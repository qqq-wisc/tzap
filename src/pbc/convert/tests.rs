use super::*;
use crate::circuit::GateKind;
use crate::pbc::{ExpandedPauli, Pauli, Phase};
use crate::semantics::test_support::assert_equivalent;

mod fuzz;

#[test]
fn litinski_four_wire_rotation_commutation_example() {
    use Pauli::{I, X, Y, Z};

    // Figure wires q1..q4 map to 0..3. Linked Z--X means CX,
    // with Z the control. X rotations are H S^{+/-1} H.
    let rx = |q, inverse| {
        [
            Gate::h(q),
            if inverse { Gate::sdg(q) } else { Gate::s(q) },
            Gate::h(q),
        ]
    };
    let mut gates = vec![
        Gate::t(0),
        Gate::cnot {
            control: 2,
            target: 1,
        },
    ];
    gates.extend(rx(3, true));
    gates.push(Gate::cnot {
        control: 1,
        target: 0,
    });
    gates.extend(rx(2, false));
    gates.extend([
        Gate::t(3),
        Gate::cnot {
            control: 3,
            target: 0,
        },
        Gate::t(0),
        Gate::s(1),
        Gate::t(2),
        Gate::s(3),
    ]);
    for q in 0..4 {
        gates.extend(rx(q, q == 0));
    }
    let original = input(4, gates);
    let converted = check(&original);
    let expected = [
        ([Z, I, I, I], 1),
        ([I, I, I, Y], -1),
        ([Z, Z, Z, Y], -1),
        ([I, X, Y, I], 1),
    ];
    for (op, (factors, eighths)) in converted.operations().iter().zip(expected) {
        let PbcOp::Rotate { axis, angle } = op else {
            panic!("expected rotation")
        };
        let expanded = converted.expand(axis.as_ref(), 10_000).unwrap();
        assert_eq!(expanded.factors, factors);
        let normalized = match expanded.phase {
            Phase::One => *angle,
            Phase::MinusOne => -*angle,
            _ => panic!("non-Hermitian axis"),
        };
        assert_eq!(normalized, PauliAngle::new(eighths));
    }
    assert_eq!(converted.operations().len(), 4);

    // Independently construct the figure's right side, including its different
    // order of commuting rotations and the original Clifford skeleton.
    let mut pictured = PbcCircuit::new(4, 0);
    for index in [0, 3, 1, 2] {
        let (factors, eighths) = expected[index];
        let mut product = pictured.identity().as_ref();
        for (q, factor) in factors.into_iter().enumerate() {
            let single = pictured.single(q as u32, factor).unwrap();
            product = pictured.product(product, single.as_ref()).unwrap();
        }
        let axis = pictured.hermitian_axis(product, 10_000).unwrap();
        pictured.rotate(axis, PauliAngle::new(eighths)).unwrap();
    }
    let skeleton: Vec<_> = original
        .gates
        .iter()
        .filter(|g| !matches!(g, Gate::t(_) | Gate::tdg(_)))
        .cloned()
        .collect();
    assert!(converted.frame_matches_gates(&skeleton));
    for gate in &skeleton {
        pictured.push_output_clifford(gate.clone()).unwrap();
    }
    // Exact operator equality (up to global phase) before the common final
    // measurements establishes equality after those measurements as well.
    assert_equivalent(&original, &pictured);

    let mut measured = original;
    measured.num_cbits = 4;
    for q in 0..4 {
        measured.gates.push(Gate::measure { qubit: q, cbit: q });
    }
    let output = to_pbc(&measured).unwrap();
    assert_eq!(output.operations().len(), 8);
    assert!(output.frame_matches_gates(&skeleton));
    // Check terminal measurement signs and retained quantum outputs via exact
    // channels, one readout at a time to bound dense Choi storage.
    use crate::semantics::channel::{ChannelLimits, circuit_channel, pbc_channel};
    for q in 0..4 {
        let mut partial = input(4, measured.gates[..measured.gates.len() - 4].to_vec());
        partial.num_cbits = 1;
        partial.gates.push(Gate::measure { qubit: q, cbit: 0 });
        let actual = to_pbc(&partial).unwrap();
        let limits = ChannelLimits::default();
        let expected = circuit_channel(&partial, &[false], limits).unwrap();
        let actual = pbc_channel(&actual, &[false], limits).unwrap();
        assert_eq!(expected.compare(&actual), Ok(()));
    }
}

fn input(n: usize, gates: Vec<Gate>) -> Circuit {
    Circuit {
        num_qubits: n,
        num_cbits: 0,
        gates,
    }
}

fn check(input: &Circuit) -> PbcCircuit {
    let output = to_pbc(input).unwrap();
    assert_equivalent(input, &output);
    assert_linear_size(input, &output);
    output
}

fn assert_linear_size(input: &Circuit, output: &PbcCircuit) {
    assert!(output.pauli_nodes().len() <= 1 + 2 * input.num_qubits + 4 * input.gates.len());
    assert!(output.operations().len() <= 7 * input.gates.len());
    assert_eq!(output.output_frame().len(), input.num_qubits);
}

#[test]
fn empty_circuits_and_typed_pass() {
    for n in 0..4 {
        let c = input(n, vec![]);
        let output = check(&c);
        assert_eq!(output.pauli_nodes().len(), 1 + 2 * n);
        assert!(output.operations().is_empty());
        assert!(output.frame_matches_gates(&[]));
        let pass: &dyn Pass<Result<PbcCircuit, PbcError>> = &ToPbc;
        assert_eq!(pass.name(), "ToPbc");
        assert_equivalent(&c, &pass.run(&c).unwrap());
    }
}

#[test]
fn exhaustive_single_qubit_words_through_length_three() {
    let alphabet = [
        Gate::h(0),
        Gate::x(0),
        Gate::z(0),
        Gate::s(0),
        Gate::sdg(0),
        Gate::t(0),
        Gate::tdg(0),
    ];
    for length in 1..=3 {
        for mut code in 0..7usize.pow(length) {
            let mut gates = vec![];
            for _ in 0..length {
                gates.push(alphabet[code % 7].clone());
                code /= 7;
            }
            check(&input(1, gates));
        }
    }
}

#[test]
fn every_clifford_pair_probed_by_x_and_z_rotations() {
    let alphabet = [
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
    ];
    for a in &alphabet {
        for b in &alphabet {
            // Interleaved T probes retain historical axes; final H/T probes X images.
            check(&input(
                2,
                vec![
                    a.clone(),
                    Gate::t(0),
                    b.clone(),
                    Gate::tdg(1),
                    Gate::h(0),
                    Gate::t(0),
                    Gate::h(1),
                    Gate::tdg(1),
                ],
            ));
        }
    }
}

#[test]
fn native_gates_all_operand_permutations_under_nontrivial_frame() {
    for a in 0..3 {
        for b in 0..3 {
            if a != b {
                let t = 3 - a - b;
                for gate in [
                    Gate::ccx {
                        control1: a,
                        control2: b,
                        target: t,
                    },
                    Gate::ccz {
                        control1: a,
                        control2: b,
                        target: t,
                    },
                ] {
                    let prefix = vec![
                        Gate::h(0),
                        Gate::s(0),
                        Gate::h(1),
                        Gate::sdg(1),
                        Gate::cnot {
                            control: 2,
                            target: 0,
                        },
                        Gate::cz {
                            control: 0,
                            target: 1,
                        },
                    ];
                    let mut gates = prefix;
                    gates.extend([gate.clone(), Gate::t(t), Gate::h(b), gate, Gate::tdg(a)]);
                    check(&input(3, gates));
                }
            }
        }
    }
}

#[test]
fn native_lowering_agrees_with_existing_gate_decomposition() {
    use crate::decompose::DecomposeToffoli;
    let c = input(
        3,
        vec![
            Gate::h(1),
            Gate::s(0),
            Gate::ccx {
                control1: 2,
                control2: 0,
                target: 1,
            },
            Gate::t(1),
            Gate::ccz {
                control1: 1,
                control2: 2,
                target: 0,
            },
            Gate::tdg(0),
        ],
    );
    let native = check(&c);
    let decomposed = DecomposeToffoli.run(&c);
    assert_equivalent(&decomposed, &native);
    assert_equivalent(&c, &to_pbc(&decomposed).unwrap());
}

#[test]
fn pure_clifford_sequence_is_preserved_by_frame() {
    let c = input(
        3,
        vec![
            Gate::h(0),
            Gate::s(2),
            Gate::sdg(1),
            Gate::x(0),
            Gate::z(2),
            Gate::cnot {
                control: 2,
                target: 0,
            },
            Gate::cz {
                control: 1,
                target: 2,
            },
        ],
    );
    let output = check(&c);
    assert!(output.operations().is_empty());
    assert!(output.frame_matches_gates(&c.gates));
}

#[test]
fn h_then_t_keeps_h_as_an_output_frame() {
    let output = check(&input(1, vec![Gate::h(0), Gate::t(0)]));
    assert!(output.frame_matches_gates(&[Gate::h(0)]));
    assert_eq!(output.operations().len(), 1);
    let PbcOp::Rotate { axis, angle } = output.operations()[0] else {
        panic!("expected an X-axis T rotation");
    };
    assert_eq!(angle, super::PauliAngle::new(1));
    assert_eq!(
        output.expand(axis.as_ref(), 100).unwrap(),
        ExpandedPauli {
            phase: Phase::One,
            factors: vec![Pauli::X],
        }
    );
}

#[test]
fn four_qubit_nonadjacent_gates_and_idle_wire() {
    check(&input(
        4,
        vec![
            Gate::h(3),
            Gate::sdg(0),
            Gate::cnot {
                control: 3,
                target: 0,
            },
            Gate::t(0),
            Gate::ccx {
                control1: 0,
                control2: 3,
                target: 2,
            },
            Gate::h(2),
            Gate::cz {
                control: 3,
                target: 2,
            },
            Gate::tdg(3),
        ],
    ));
}

#[test]
fn repeated_conversion_is_deterministic_and_preserves_input_metadata() {
    let mut c = input(
        2,
        vec![
            Gate::h(0),
            Gate::s(1),
            Gate::cnot {
                control: 0,
                target: 1,
            },
            Gate::t(1),
        ],
    );
    // Classical register size is metadata, not an allocation in the converter.
    c.num_cbits = usize::MAX;
    let before = c.clone();
    let a = ToPbc.run(&c).unwrap();
    ToPbc.run(&input(3, vec![Gate::h(2)])).unwrap();
    let b = ToPbc.run(&c).unwrap();
    assert_eq!(a.to_ascii().unwrap(), b.to_ascii().unwrap());
    assert_eq!(a.num_cbits(), usize::MAX);
    assert_eq!(c.gates, before.gates);
    assert_eq!(c.num_qubits, before.num_qubits);
    assert_eq!(c.num_cbits, before.num_cbits);
    assert_equivalent(&c, &a);
    assert_equivalent(&c, &b);
}

#[test]
fn emitted_axes_are_snapshots_not_mutable_frame_entries() {
    let c = input(
        1,
        vec![
            Gate::t(0),
            Gate::h(0),
            Gate::t(0),
            Gate::s(0),
            Gate::h(0),
            Gate::tdg(0),
        ],
    );
    let output = check(&c);
    let axes: Vec<_> = output
        .operations()
        .iter()
        .map(|op| output.expand(op.axis().as_ref(), 100).unwrap())
        .collect();
    assert_eq!(
        axes,
        vec![
            ExpandedPauli {
                phase: Phase::One,
                factors: vec![Pauli::Z]
            },
            ExpandedPauli {
                phase: Phase::One,
                factors: vec![Pauli::X]
            },
            ExpandedPauli {
                phase: Phase::One,
                factors: vec![Pauli::Y]
            },
        ]
    );
}

#[test]
fn terminal_measurements_preserve_frame_and_immutable_outcomes() {
    let c = Circuit {
        num_qubits: 2,
        num_cbits: 1,
        gates: vec![
            Gate::h(0),
            Gate::cnot {
                control: 0,
                target: 1,
            },
            Gate::measure { qubit: 1, cbit: 0 },
            Gate::measure { qubit: 1, cbit: 0 },
        ],
    };
    let output = to_pbc(&c).unwrap();
    assert_eq!(output.measurement_count(), 2);
    assert_eq!(output.operations().len(), 2);
    assert!(output.frame_matches_gates(&c.gates[..2]));
    for (index, op) in output.operations().iter().enumerate() {
        let PbcOp::Measure {
            outcome,
            target: Some(0),
            ..
        } = op
        else {
            panic!("measurement");
        };
        assert_eq!(outcome.index(), index);
        assert_eq!(
            output.expand(op.axis().as_ref(), 100).unwrap(),
            ExpandedPauli {
                phase: Phase::One,
                factors: vec![Pauli::X, Pauli::Z]
            }
        );
    }
    assert_linear_size(&c, &output);
}

#[test]
fn reject_resets_and_any_gate_after_measurement() {
    let mut c = input(1, vec![Gate::h(0), Gate::reset(0)]);
    assert_eq!(
        to_pbc(&c).unwrap_err(),
        PbcError::InvalidInput {
            index: 1,
            cause: Box::new(PbcError::UnsupportedGate(GateKind::Reset))
        }
    );
    c.num_cbits = 1;
    for gates in [
        vec![Gate::reset(0)],
        vec![Gate::measure { qubit: 0, cbit: 0 }, Gate::reset(0)],
    ] {
        c.gates = gates;
        assert_eq!(
            to_pbc(&c).unwrap_err(),
            PbcError::InvalidInput {
                index: c.gates.len() - 1,
                cause: Box::new(PbcError::UnsupportedGate(GateKind::Reset)),
            }
        );
    }
}

/// Gates after a measurement, on the measured or any other wire, convert with
/// exactly the channel of the gate circuit. The measurement itself leaves the
/// frame unchanged, so later rotations use the same images as without it.
#[test]
fn mid_circuit_measurements_match_exact_channels() {
    use crate::semantics::channel::{ChannelLimits, circuit_channel, pbc_channel};
    let limits = ChannelLimits::default();
    let prefixes = [
        vec![],
        vec![Gate::h(0)],
        vec![
            Gate::h(0),
            Gate::cnot {
                control: 0,
                target: 1,
            },
        ],
        vec![
            Gate::s(1),
            Gate::h(1),
            Gate::cz {
                control: 1,
                target: 0,
            },
        ],
    ];
    let suffixes = [
        Gate::t(0),
        Gate::tdg(1),
        Gate::h(0),
        Gate::cnot {
            control: 1,
            target: 0,
        },
    ];
    for prefix in &prefixes {
        for q in 0..2 {
            for suffix in &suffixes {
                let mut c = input(2, prefix.clone());
                c.num_cbits = 2;
                c.gates.push(Gate::measure { qubit: q, cbit: 0 });
                // H then T on the measured wire gives a rotation that does not
                // commute with the measurement, so its position matters.
                c.gates.extend([
                    suffix.clone(),
                    Gate::h(q),
                    Gate::t(q),
                    Gate::t(1),
                    Gate::h(1),
                ]);
                c.gates.push(Gate::measure { qubit: 1, cbit: 1 });
                let output = to_pbc(&c).unwrap();
                assert_eq!(output.measurement_count(), 2);
                assert_linear_size(&c, &output);
                for initial in [[false, true], [true, false]] {
                    let expected = circuit_channel(&c, &initial, limits).unwrap();
                    let actual = pbc_channel(&output, &initial, limits).unwrap();
                    assert_eq!(expected.compare(&actual), Ok(()), "{c}");
                }
            }
        }
    }
}

#[test]
fn invalid_inputs_report_instruction_index() {
    let cases = [
        (Gate::h(2), PbcError::QubitOutOfRange(2)),
        (Gate::rz(0.0, 0), PbcError::UnsupportedGate(GateKind::Rz)),
        (
            Gate::measure { qubit: 0, cbit: 0 },
            PbcError::ClassicalBitOutOfRange(0),
        ),
        (Gate::reset(2), PbcError::UnsupportedGate(GateKind::Reset)),
        (
            Gate::cnot {
                control: 0,
                target: 0,
            },
            PbcError::RepeatedOperand,
        ),
        (
            Gate::cz {
                control: 1,
                target: 1,
            },
            PbcError::RepeatedOperand,
        ),
        (
            Gate::ccx {
                control1: 0,
                control2: 0,
                target: 1,
            },
            PbcError::RepeatedOperand,
        ),
        (
            Gate::ccz {
                control1: 0,
                control2: 1,
                target: 1,
            },
            PbcError::RepeatedOperand,
        ),
    ];
    for (gate, cause) in cases {
        assert_eq!(
            to_pbc(&input(2, vec![Gate::h(0), gate])).unwrap_err(),
            PbcError::InvalidInput {
                index: 1,
                cause: Box::new(cause)
            }
        );
    }
    assert_eq!(
        to_pbc(&input(usize::MAX, vec![])).unwrap_err(),
        PbcError::TooManyQubits
    );
}

#[test]
fn growing_parity_uses_linear_nodes_not_quadratic_strings() {
    for n in [8, 128, 4096, 32_768] {
        let mut c = input(n, vec![]);
        for q in 1..n as u32 {
            c.gates.extend([
                Gate::cnot {
                    control: q,
                    target: 0,
                },
                Gate::t(0),
            ]);
        }
        let output = to_pbc(&c).unwrap();
        assert_eq!(output.pauli_nodes().len(), 1 + 2 * n + 2 * (n - 1));
        assert_eq!(output.operations().len(), n - 1);
        assert_eq!(output.output_frame().len(), n);
        assert_linear_size(&c, &output);
        // Only the tiny instance is expanded; its supports grow 2,3,...,n.
        if n == 8 {
            for (i, op) in output.operations().iter().enumerate() {
                let axis = output.expand(op.axis().as_ref(), 1000).unwrap();
                assert_eq!(
                    axis.factors.iter().filter(|p| **p != Pauli::I).count(),
                    i + 2
                );
            }
        }
    }
}

#[test]
fn long_native_stream_has_fixed_cost_per_gate() {
    let c = input(
        3,
        vec![
            Gate::ccx {
                control1: 0,
                control2: 1,
                target: 2
            };
            10_000
        ],
    );
    let output = to_pbc(&c).unwrap();
    assert_eq!(output.pauli_nodes().len(), 7 + 4 * c.gates.len());
    assert_eq!(output.operations().len(), 7 * c.gates.len());
    assert!(output.frame_matches_gates(&[]));
    assert_linear_size(&c, &output);
}

#[test]
fn errors_display_their_instruction_and_expose_their_cause() {
    use std::error::Error;
    let err = to_pbc(&input(1, vec![Gate::h(0), Gate::rz(0.5, 0)])).unwrap_err();
    assert_eq!(
        err.to_string(),
        "input gate 1: unsupported PBC input gate: Rz"
    );
    assert_eq!(
        err.source().unwrap().to_string(),
        PbcError::UnsupportedGate(GateKind::Rz).to_string()
    );
    assert!(PbcError::ExpansionLimit.source().is_none());
}
