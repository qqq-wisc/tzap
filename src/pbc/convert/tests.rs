use super::*;
use crate::pbc::{ExpandedPauli, Pauli};
use crate::semantics::test_support::assert_equivalent;

mod fuzz;

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
    assert!(output.output_cliffords().len() <= input.gates.len());
}

#[test]
fn empty_circuits_and_typed_pass() {
    for n in 0..4 {
        let c = input(n, vec![]);
        let output = check(&c);
        assert_eq!(output.pauli_nodes().len(), 1 + 2 * n);
        assert!(output.operations().is_empty());
        assert!(output.output_cliffords().is_empty());
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
fn pure_clifford_suffix_is_preserved_in_order() {
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
    assert_eq!(output.output_cliffords(), c.gates);
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
fn measurements_and_resets_keep_frame_and_immutable_outcomes() {
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
            Gate::reset(1),
            Gate::x(0),
            Gate::measure { qubit: 1, cbit: 0 },
            Gate::t(1),
        ],
    };
    let output = to_pbc(&c).unwrap();
    assert_eq!(output.num_cbits(), 1);
    assert_eq!(output.measurement_count(), 3);
    assert_eq!(output.operations().len(), 5);
    let ops = output.operations();
    let PbcOp::Measure {
        outcome: first,
        target: Some(0),
        ..
    } = ops[0]
    else {
        panic!("measurement")
    };
    let PbcOp::Measure {
        outcome: hidden,
        target: None,
        ..
    } = ops[1]
    else {
        panic!("reset measurement")
    };
    let PbcOp::ConditionalRotate { if_one, angle, .. } = ops[2] else {
        panic!("reset correction")
    };
    let PbcOp::Measure {
        outcome: last,
        target: Some(0),
        ..
    } = ops[3]
    else {
        panic!("overwrite")
    };
    assert_eq!((first.index(), hidden.index(), last.index()), (0, 1, 2));
    assert_eq!(if_one, hidden);
    assert_eq!(angle, PauliAngle::new(4));
    for i in [0, 1, 3, 4] {
        assert_eq!(
            output.expand(ops[i].axis().as_ref(), 100).unwrap(),
            ExpandedPauli {
                phase: Phase::One,
                factors: vec![Pauli::X, Pauli::Z]
            }
        );
    }
    assert_eq!(
        output.expand(ops[2].axis().as_ref(), 100).unwrap().factors,
        vec![Pauli::I, Pauli::X]
    );
    assert_linear_size(&c, &output);
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
        (Gate::reset(2), PbcError::QubitOutOfRange(2)),
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
        assert_eq!(output.output_cliffords().len(), n - 1);
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
    assert!(output.output_cliffords().is_empty());
    assert_linear_size(&c, &output);
}
