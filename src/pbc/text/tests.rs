use super::*;
use crate::circuit::{Circuit, Gate};
use crate::pbc::{Pauli, PauliAngle, to_pbc};

/// Test-only reader, independent of the original circuit and converter.
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
    let mut in_frame = false;
    let mut seen_rows = std::collections::HashSet::new();
    for line in lines {
        let mut words = line.split_whitespace();
        let op = words.next().unwrap();
        if op == "frame" {
            in_frame = true;
            let row = words.next().unwrap();
            let is_z = match &row[..1] {
                "X" => false,
                "Z" => true,
                _ => panic!("bad frame generator"),
            };
            let q: u32 = row[1..].parse().unwrap();
            assert!((q as usize) < n && seen_rows.insert((is_z, q)));
            let sign = match words.next().unwrap() {
                "1" => Phase::One,
                "-1" => Phase::MinusOne,
                _ => panic!("bad frame sign"),
            };
            let mut value = p.identity().as_ref().scaled(sign);
            let mut previous = None;
            for word in words {
                let factor_q = word[1..].parse::<u32>().unwrap();
                assert!(previous.is_none_or(|previous| previous < factor_q));
                previous = Some(factor_q);
                let factor = match &word[..1] {
                    "X" => Pauli::X,
                    "Y" => Pauli::Y,
                    "Z" => Pauli::Z,
                    _ => panic!("bad Pauli"),
                };
                let single = p.single(factor_q, factor).unwrap();
                value = p.product(value, single.as_ref()).unwrap();
            }
            let value = p.hermitian_axis(value, 100_000).unwrap().as_ref();
            if is_z {
                p.output_frame.z[q as usize] = value;
            } else {
                p.output_frame.x[q as usize] = value;
            }
            continue;
        }
        assert!(!in_frame, "frame rows must be trailing");
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
fn exported_programs_match_exact_unitaries_for_every_clifford_pair() {
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
fn exported_programs_preserve_exact_partial_readout_channels() {
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
            assert_eq!(actual.output_frame().len(), 2);
            let initial = [true, false];
            let limits = ChannelLimits::default();
            let expected = circuit_channel(&input, &initial, limits).unwrap();
            let actual = pbc_channel(&actual, &initial, limits).unwrap();
            assert_eq!(expected.compare(&actual), Ok(()), "{text}");
        }
    }
}

fn assert_full_export(input: &Circuit, initial: &[bool]) {
    use crate::semantics::channel::{ChannelLimits, circuit_channel, pbc_channel};
    let text = to_pbc(input).unwrap().to_text().unwrap();
    let exported = read_text(&text);
    let limits = ChannelLimits::default();
    let expected = circuit_channel(input, initial, limits).unwrap();
    let actual = pbc_channel(&exported, initial, limits).unwrap();
    assert_eq!(expected.compare(&actual), Ok(()), "{text}");
    assert_eq!(
        expected.classical_blocks(),
        actual.classical_blocks(),
        "{text}"
    );
}

#[test]
fn export_preserves_full_readout_channel_for_every_clifford_pair() {
    for a in clifford_alphabet() {
        for b in clifford_alphabet() {
            let input = Circuit {
                num_qubits: 2,
                num_cbits: 3,
                gates: vec![
                    a.clone(),
                    Gate::t(0),
                    b,
                    Gate::tdg(1),
                    Gate::measure { qubit: 0, cbit: 1 },
                    Gate::measure { qubit: 1, cbit: 0 },
                ],
            };
            assert_full_export(&input, &[false, true, true]);
        }
    }
}

#[test]
fn export_preserves_native_gates_and_overwritten_registers() {
    for native in [
        Gate::ccx {
            control1: 0,
            control2: 2,
            target: 1,
        },
        Gate::ccz {
            control1: 2,
            control2: 1,
            target: 0,
        },
    ] {
        for overwrite in [false, true] {
            let input = Circuit {
                num_qubits: 3,
                num_cbits: 3,
                gates: vec![
                    Gate::h(0),
                    Gate::s(0),
                    Gate::x(2),
                    native.clone(),
                    Gate::tdg(1),
                    Gate::h(1),
                    Gate::measure { qubit: 2, cbit: 2 },
                    Gate::measure { qubit: 0, cbit: 0 },
                    Gate::measure {
                        qubit: 1,
                        cbit: if overwrite { 0 } else { 1 },
                    },
                ],
            };
            assert_full_export(&input, &[true, false, true]);
        }
    }
}

#[test]
fn classical_projection_distinguishes_readout_axes_but_discards_output_frame() {
    use crate::semantics::channel::{ChannelLimits, circuit_channel, pbc_channel};
    let input = Circuit {
        num_qubits: 1,
        num_cbits: 1,
        gates: vec![Gate::h(0), Gate::measure { qubit: 0, cbit: 0 }],
    };
    let limits = ChannelLimits::default();
    let expected = circuit_channel(&input, &[false], limits).unwrap();
    let correct = pbc_channel(
        &read_text("qubits 1\nregisters 1\nm 1 X0 -> c0\n"),
        &[false],
        limits,
    )
    .unwrap();
    let wrong = pbc_channel(
        &read_text("qubits 1\nregisters 1\nm 1 Z0 -> c0\n"),
        &[false],
        limits,
    )
    .unwrap();
    assert!(
        expected.compare(&correct).is_err(),
        "quantum outputs differ"
    );
    assert_eq!(expected.classical_blocks(), correct.classical_blocks());
    assert_ne!(expected.classical_blocks(), wrong.classical_blocks());
    assert_full_export(&input, &[false]);
}

#[test]
fn export_preserves_partial_and_absent_readout_channels() {
    for measurements in [vec![], vec![Gate::measure { qubit: 1, cbit: 0 }]] {
        let mut input = Circuit {
            num_qubits: 2,
            num_cbits: 2,
            gates: vec![
                Gate::h(0),
                Gate::t(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        };
        input.gates.extend(measurements);
        assert_full_export(&input, &[true, false]);
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
fn full_readout_keeps_output_frame_and_quantum_outputs() {
    let c = Circuit {
        num_qubits: 1,
        num_cbits: 1,
        gates: vec![Gate::h(0), Gate::measure { qubit: 0, cbit: 0 }],
    };
    let p = to_pbc(&c).unwrap();
    assert_eq!(
        p.to_text().unwrap(),
        "qubits 1\nregisters 1\nm 1 X0 -> c0\nframe X0 1 Z0\nframe Z0 1 X0\n"
    );
    assert!(p.frame_matches_gates(&[Gate::h(0)]));
}

#[test]
fn all_clifford_gates_export_as_frame_and_empty_program() {
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
        "qubits 2\nregisters 0\nframe X0 -1 Y0 Z1\nframe X1 -1 Z0 X1\nframe Z0 -1 X0 Z1\nframe Z1 -1 Z1\n"
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

#[test]
fn text_budget_covers_operations_and_output_frame_together() {
    let mut p = PbcCircuit::new(1, 0);
    let z = p.z(0).unwrap();
    p.rotate(z, PauliAngle::new(1)).unwrap();
    assert_eq!(
        p.to_text_with(TextOptions {
            max_expansion_cells: 8,
        }),
        Err(PbcError::ExpansionLimit)
    );
    assert_eq!(
        p.to_text_with(TextOptions {
            max_expansion_cells: 9,
        }),
        Ok("qubits 1\nregisters 0\nr 1 1 Z0\n".into())
    );
}

/// Seeded differential fuzzing of the exporter itself: random gate prefixes
/// with random terminal readouts (subsets, orders, overwritten registers),
/// exported to text, re-read independently, and compared as exact channels.
#[test]
fn seeded_export_fuzz_preserves_exact_channels() {
    use rand::{Rng, SeedableRng, rngs::StdRng, seq::SliceRandom};
    for case in 0..48u64 {
        let mut rng = StdRng::seed_from_u64(0x5042_4354_0000 + case);
        let n = 1 + (case % 3) as usize;
        let cbits = 1 + (case / 3 % 3) as usize;
        let mut input = Circuit {
            num_qubits: n,
            num_cbits: cbits,
            gates: vec![],
        };
        for _ in 0..rng.gen_range(3..=12) {
            let mut q: Vec<u32> = (0..n as u32).collect();
            q.shuffle(&mut rng);
            input.gates.push(match rng.gen_range(0..11) {
                0 => Gate::h(q[0]),
                1 => Gate::x(q[0]),
                2 => Gate::z(q[0]),
                3 => Gate::s(q[0]),
                4 => Gate::sdg(q[0]),
                5 => Gate::t(q[0]),
                6 => Gate::tdg(q[0]),
                7 if n >= 2 => Gate::cnot {
                    control: q[0],
                    target: q[1],
                },
                8 if n >= 2 => Gate::cz {
                    control: q[0],
                    target: q[1],
                },
                9 if n >= 3 => Gate::ccx {
                    control1: q[0],
                    control2: q[1],
                    target: q[2],
                },
                10 if n >= 3 => Gate::ccz {
                    control1: q[0],
                    control2: q[1],
                    target: q[2],
                },
                _ => Gate::h(q[0]),
            });
        }
        // At most three readouts keeps the exact channel within test limits.
        for _ in 0..rng.gen_range(0..=3) {
            input.gates.push(Gate::measure {
                qubit: rng.gen_range(0..n as u32),
                cbit: rng.gen_range(0..cbits as u32),
            });
        }
        let initial: Vec<bool> = (0..cbits).map(|_| rng.gen_bool(0.5)).collect();
        assert_full_export(&input, &initial);
    }
}
