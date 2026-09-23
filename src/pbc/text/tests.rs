use super::*;
use crate::circuit::Circuit;
use crate::pbc::{PauliAngle, to_pbc};

/// Test-only reader: reconstruct exactly what the exported r/m instructions say,
/// without consulting the original circuit or the Clifford lowering identities.
fn read_text(text: &str) -> PbcCircuit {
    let mut lines = text.lines();
    let n = lines
        .next()
        .unwrap()
        .strip_prefix("qubits ")
        .unwrap()
        .parse()
        .unwrap();
    let c = lines
        .next()
        .unwrap()
        .strip_prefix("registers ")
        .unwrap()
        .parse()
        .unwrap();
    let mut p = PbcCircuit::new(n, c);
    for line in lines {
        let mut words = line.split_whitespace();
        let op = words.next().unwrap();
        assert!(matches!(op, "r" | "m"));
        let k = if op == "r" {
            words.next().unwrap().parse().unwrap()
        } else {
            0
        };
        let sign = match words.next().unwrap() {
            "1" => Phase::One,
            "-1" => Phase::MinusOne,
            _ => panic!("bad sign"),
        };
        let mut axis = p.identity().as_ref().scaled(sign);
        let mut previous = None;
        for word in words.by_ref() {
            if word == "->" {
                break;
            }
            let q = word[1..].parse::<u32>().unwrap();
            assert!(previous.is_none_or(|previous| previous < q));
            previous = Some(q);
            let factor = match &word[..1] {
                "X" => Pauli::X,
                "Y" => Pauli::Y,
                "Z" => Pauli::Z,
                _ => panic!("bad Pauli"),
            };
            let single = p.single(q, factor).unwrap();
            axis = p.product(axis, single.as_ref()).unwrap();
        }
        let axis = p.hermitian_axis(axis, 100_000).unwrap();
        if op == "r" {
            p.rotate(axis, PauliAngle::new(k)).unwrap();
        } else {
            let c = words
                .next()
                .unwrap()
                .strip_prefix('c')
                .unwrap()
                .parse()
                .unwrap();
            p.measure(axis, Some(c)).unwrap();
        }
        assert!(words.next().is_none());
    }
    p
}

fn clifford_alphabet() -> Vec<Gate> {
    vec![
        Gate::h(0),
        Gate::h(1),
        Gate::x(0),
        Gate::z(1),
        Gate::s(0),
        Gate::sdg(1),
        Gate::cnot {
            control: 0,
            target: 1,
        },
        Gate::cnot {
            control: 1,
            target: 0,
        },
        Gate::cz {
            control: 0,
            target: 1,
        },
        Gate::cz {
            control: 1,
            target: 0,
        },
    ]
}

#[test]
fn exported_rotations_match_exact_unitaries_for_every_clifford_pair() {
    use crate::semantics::test_support::assert_equivalent;
    for a in clifford_alphabet() {
        for b in clifford_alphabet() {
            let input = Circuit {
                num_qubits: 2,
                num_cbits: 0,
                gates: vec![a.clone(), Gate::t(0), b, Gate::tdg(1)],
            };
            let text = to_pbc(&input).unwrap().to_text().unwrap();
            assert_equivalent(&input, &read_text(&text));
        }
    }
}

#[test]
fn exported_rotations_preserve_exact_partial_readout_channels() {
    use crate::semantics::channel::{ChannelLimits, circuit_channel, pbc_channel};
    for gate in clifford_alphabet() {
        for q in 0..2 {
            let input = Circuit {
                num_qubits: 2,
                num_cbits: 2,
                gates: vec![
                    Gate::h(0),
                    Gate::t(0),
                    gate.clone(),
                    Gate::measure { qubit: q, cbit: 1 },
                ],
            };
            let text = to_pbc(&input).unwrap().to_text().unwrap();
            let actual = read_text(&text);
            assert!(actual.output_cliffords().is_empty());
            let initial = [true, false];
            let limits = ChannelLimits::default();
            let expected = circuit_channel(&input, &initial, limits).unwrap();
            let actual = pbc_channel(&actual, &initial, limits).unwrap();
            assert_eq!(expected.compare(&actual), Ok(()), "{text}");
        }
    }
}

#[test]
fn sparse_factors_signed_angles_and_overwritten_registers() {
    let mut p = PbcCircuit::new(12, 3);
    let x = p.x(0).unwrap();
    let z = p.z(11).unwrap();
    let product = p.product(z.as_ref(), x.as_ref()).unwrap();
    let axis = p.hermitian_axis(product, 100).unwrap();
    p.rotate(axis.negated(), PauliAngle::new(-1)).unwrap();
    let y = p.y(2).unwrap();
    p.rotate(y, PauliAngle::new(2)).unwrap();
    p.measure(axis, Some(2)).unwrap();
    p.measure(z.negated(), Some(2)).unwrap();
    assert_eq!(
        p.to_text().unwrap(),
        "qubits 12\nregisters 3\nr -1 -1 X0 Z11\nr 2 1 Y2\nm 1 X0 Z11 -> c2\nm -1 Z11 -> c2\n"
    );
}

#[test]
fn shared_products_preserve_y_and_signs() {
    let mut p = PbcCircuit::new(1, 1);
    let x = p.x(0).unwrap();
    let z = p.z(0).unwrap();
    let xz = p.product(x.as_ref(), z.as_ref()).unwrap();
    let y = p.hermitian_axis(xz.scaled(Phase::I), 100).unwrap();
    p.rotate(y, PauliAngle::new(1)).unwrap();
    let square = p.product(xz, xz).unwrap();
    let minus_identity = p.hermitian_axis(square, 100).unwrap();
    p.measure(minus_identity, Some(0)).unwrap();
    assert_eq!(
        p.to_text().unwrap(),
        "qubits 1\nregisters 1\nr 1 1 Y0\nm -1 -> c0\n"
    );
}

#[test]
fn angles_are_integer_multiples_modulo_global_phase() {
    let mut p = PbcCircuit::new(1, 0);
    let z = p.z(0).unwrap();
    for k in 0..8 {
        p.rotate(z, PauliAngle::new(k)).unwrap();
    }
    assert_eq!(
        p.to_text().unwrap(),
        "qubits 1\nregisters 0\nr 0 1 Z0\nr 1 1 Z0\nr 2 1 Z0\nr 3 1 Z0\nr 4 1 Z0\nr -3 1 Z0\nr -2 1 Z0\nr -1 1 Z0\n"
    );
}

#[test]
fn quantum_outputs_are_default_classical_projection_is_explicit() {
    let c = Circuit {
        num_qubits: 1,
        num_cbits: 1,
        gates: vec![Gate::h(0), Gate::measure { qubit: 0, cbit: 0 }],
    };
    let p = to_pbc(&c).unwrap();
    assert_eq!(
        p.to_text().unwrap(),
        "qubits 1\nregisters 1\nm 1 X0 -> c0\nr 2 1 Z0\nr 2 1 X0\nr 2 1 Z0\n"
    );
    let classical = p
        .to_text_with(TextOptions {
            output_semantics: OutputSemantics::Classical,
            ..TextOptions::default()
        })
        .unwrap();
    assert_eq!(classical, "qubits 1\nregisters 1\nm 1 X0 -> c0\n");
    assert_eq!(p.output_cliffords(), &[Gate::h(0)]);
}

#[test]
fn all_suffix_gates_and_empty_program() {
    assert_eq!(
        PbcCircuit::new(0, 2).to_text().unwrap(),
        "qubits 0\nregisters 2\n"
    );
    let c = Circuit {
        num_qubits: 2,
        num_cbits: 0,
        gates: vec![
            Gate::h(0),
            Gate::x(1),
            Gate::z(0),
            Gate::s(1),
            Gate::sdg(0),
            Gate::cnot {
                control: 1,
                target: 0,
            },
            Gate::cz {
                control: 0,
                target: 1,
            },
        ],
    };
    assert_eq!(
        to_pbc(&c).unwrap().to_text().unwrap(),
        "qubits 2\nregisters 0\nr 2 1 Z0\nr 2 1 X0\nr 2 1 Z0\nr 4 1 X1\nr 4 1 Z0\nr 2 1 Z1\nr -2 1 Z0\nr 2 1 Z1\nr 2 1 X0\nr -2 1 X0 Z1\nr 2 1 Z0\nr 2 1 Z1\nr -2 1 Z0 Z1\n"
    );
}

#[test]
fn unsupported_operations_and_expansion_limit_fail_without_truncation() {
    let mut p = PbcCircuit::new(1, 1);
    let z = p.z(0).unwrap();
    let m = p.measure(z, Some(0)).unwrap();
    assert_eq!(
        p.to_text_with(TextOptions {
            max_expansion_cells: 0,
            ..TextOptions::default()
        }),
        Err(PbcError::ExpansionLimit)
    );
    p.conditional_pauli(z, m).unwrap();
    assert_eq!(
        p.to_text(),
        Err(PbcError::UnsupportedTextOperation { index: 1 })
    );
    let mut hidden = PbcCircuit::new(1, 0);
    let z = hidden.z(0).unwrap();
    hidden.measure(z, None).unwrap();
    assert_eq!(
        hidden.to_text(),
        Err(PbcError::UnsupportedTextOperation { index: 0 })
    );
    let mut late_rotation = PbcCircuit::new(1, 1);
    let z = late_rotation.z(0).unwrap();
    late_rotation.measure(z, Some(0)).unwrap();
    late_rotation.rotate(z, PauliAngle::new(1)).unwrap();
    assert_eq!(
        late_rotation.to_text().unwrap(),
        "qubits 1\nregisters 1\nm 1 Z0 -> c0\nr 1 1 Z0\n"
    );
}
