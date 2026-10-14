use super::*;
use crate::circuit::GateKind;
use crate::pbc::{ExpandedPauli, Pauli, Phase};
use crate::semantics::test_support::assert_equivalent;

mod fuzz;
mod max_weight;

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
    let output = to_pbc(&measured, None).unwrap();
    assert_eq!(output.operations().len(), 8);
    assert!(output.frame_matches_gates(&skeleton));
    // Check terminal measurement signs and retained quantum outputs via exact
    // channels, one readout at a time to bound dense Choi storage.
    use crate::semantics::channel::{ChannelLimits, circuit_channel, pbc_channel};
    for q in 0..4 {
        let mut partial = input(4, measured.gates[..measured.gates.len() - 4].to_vec());
        partial.num_cbits = 1;
        partial.gates.push(Gate::measure { qubit: q, cbit: 0 });
        let actual = to_pbc(&partial, None).unwrap();
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
    let output = to_pbc(input, None).unwrap();
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
        let pass: &dyn Pass<Result<PbcCircuit, PbcError>> = &ToPbc::default();
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
    assert_equivalent(&c, &to_pbc(&decomposed, None).unwrap());
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
    let a = ToPbc::default().run(&c).unwrap();
    ToPbc::default().run(&input(3, vec![Gate::h(2)])).unwrap();
    let b = ToPbc::default().run(&c).unwrap();
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
    let output = to_pbc(&c, None).unwrap();
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
        to_pbc(&c, None).unwrap_err(),
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
            to_pbc(&c, None).unwrap_err(),
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
                let output = to_pbc(&c, None).unwrap();
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
        (
            Gate::rz_f64(0.1, 0).unwrap(),
            PbcError::UnsupportedGate(GateKind::Rz),
        ),
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
            to_pbc(&input(2, vec![Gate::h(0), gate]), None).unwrap_err(),
            PbcError::InvalidInput {
                index: 1,
                cause: Box::new(cause)
            }
        );
    }
    assert_eq!(
        to_pbc(&input(usize::MAX, vec![]), None).unwrap_err(),
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
        let output = to_pbc(&c, None).unwrap();
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
    let output = to_pbc(&c, None).unwrap();
    assert_eq!(output.pauli_nodes().len(), 7 + 4 * c.gates.len());
    assert_eq!(output.operations().len(), 7 * c.gates.len());
    assert!(output.frame_matches_gates(&[]));
    assert_linear_size(&c, &output);
}

#[test]
fn errors_display_their_instruction_and_expose_their_cause() {
    use std::error::Error;
    let err = to_pbc(
        &input(1, vec![Gate::h(0), Gate::rz_f64(0.5, 0).unwrap()]),
        None,
    )
    .unwrap_err();
    assert_eq!(
        err.to_string(),
        "input gate 1: unsupported PBC input gate: Rz; run DecomposeRotations before PBC conversion"
    );
    assert_eq!(
        err.source().unwrap().to_string(),
        PbcError::UnsupportedGate(GateKind::Rz).to_string()
    );
    assert!(PbcError::ExpansionLimit.source().is_none());
}

/// Unbounded and the weight bounds that flush the frame most often.
const ALL_BOUNDS: [Option<NonZeroUsize>; 3] = [None, NonZeroUsize::new(1), NonZeroUsize::new(2)];

/// The PBC unitary of `to_pbc(input)` under each weight bound, compared with
/// the native floating-point interpreter (the exact reference interpreter only
/// models the PBC input gates) up to global phase.
fn assert_pbc_matches_native(input: &Circuit, bounds: &[Option<NonZeroUsize>]) {
    for &max_weight in bounds {
        let pbc = to_pbc(input, max_weight)
            .unwrap_or_else(|e| panic!("{e} (max_weight {max_weight:?}):\n{input}"));
        let actual =
            crate::semantics::pbc_unitary(&pbc, crate::semantics::Limits::default()).unwrap();
        crate::semantics::test_support::check_numeric(&actual, input);
    }
}

/// Every statement `to_pbc` takes without a decomposition pass: each gate
/// kind on every operand order of three qubits, every quarter-turn angle of a
/// single-qubit rotation, and every Clifford angle of a controlled rotation
/// over two of its periods (CR*(2*pi) is a controlled -I, not the identity).
fn directly_convertible_statements() -> Vec<String> {
    let mut statements = Vec::new();
    for q in 0..3 {
        for name in ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "sx"] {
            statements.push(format!("{name} q[{q}]"));
        }
        for k in -8..=8 {
            for name in ["rz", "p", "rx", "ry"] {
                statements.push(format!("{name}({k}*pi/4) q[{q}]"));
            }
        }
    }
    for (a, b) in [(0, 1), (1, 0), (0, 2), (2, 0), (1, 2), (2, 1)] {
        for name in ["cx", "cz", "cy", "swap"] {
            statements.push(format!("{name} q[{a}],q[{b}]"));
        }
        for k in -4..=4 {
            for name in ["cp", "crx", "cry", "crz"] {
                statements.push(format!("{name}({k}*pi) q[{a}],q[{b}]"));
            }
        }
    }
    for (a, b, c) in [
        (0, 1, 2),
        (0, 2, 1),
        (1, 0, 2),
        (1, 2, 0),
        (2, 0, 1),
        (2, 1, 0),
    ] {
        for name in ["ccx", "ccz"] {
            statements.push(format!("{name} q[{a}],q[{b}],q[{c}]"));
        }
    }
    statements
}

#[test]
fn every_directly_convertible_gate_has_the_native_unitary() {
    // Superposed, phased, entangled inputs make every matrix entry count.
    let prefix = "qreg q[3]; h q[0]; h q[1]; t q[1]; h q[2]; s q[2]; cx q[0],q[2];";
    let statements = directly_convertible_statements();
    let mut bounded = std::collections::HashSet::new();
    for statement in &statements {
        let input = Circuit::from_qasm(&format!("{prefix} {statement};")).unwrap();
        assert!(
            !needs_decomposition(input.gates.last().unwrap()),
            "{statement}"
        );
        // Weight bounds act after lowering, so one statement per gate kind
        // exercises them; the exact interpreter is too slow for all of them.
        let name = statement.split(['(', ' ']).next().unwrap();
        let bounds = if bounded.insert(name) {
            &ALL_BOUNDS[..]
        } else {
            &ALL_BOUNDS[..1]
        };
        assert_pbc_matches_native(&input, bounds);
    }
    assert_eq!(statements.len(), 483);
}

#[test]
fn every_gate_converts_after_its_suggested_decomposition() {
    use crate::decompose::{DecomposeRotations, DecomposeToffoli};
    let prefix = "qreg q[3]; h q[0]; h q[1]; t q[1]; h q[2]; s q[2]; cx q[0],q[2];";
    for statement in [
        "ch q[0],q[1]",
        "ch q[2],q[0]",
        "cswap q[0],q[1],q[2]",
        "cswap q[2],q[0],q[1]",
        // Exact non-Clifford angles decompose without approximation.
        "crz(pi/2) q[0],q[1]",
        "crz(-3*pi/2) q[1],q[0]",
        "crx(pi/2) q[1],q[2]",
        "cry(-pi/2) q[2],q[0]",
        "cp(pi/2) q[0],q[2]",
        "cp(3*pi/2) q[2],q[1]",
    ] {
        let input = Circuit::from_qasm(&format!("{prefix} {statement};")).unwrap();
        assert!(
            needs_decomposition(input.gates.last().unwrap()),
            "{statement}"
        );
        assert!(to_pbc(&input, None).is_err(), "{statement}");
        let decomposed = if statement.starts_with("cswap") {
            DecomposeToffoli.run(&input)
        } else {
            DecomposeRotations::default().try_run(&input).unwrap()
        };
        let pbc = to_pbc(&decomposed, None).unwrap();
        let actual =
            crate::semantics::pbc_unitary(&pbc, crate::semantics::Limits::default()).unwrap();
        // Against the original native gate, not just the decomposition.
        crate::semantics::test_support::check_numeric(&actual, &input);
    }
}

#[test]
fn random_mixtures_of_convertible_gates_have_the_native_unitary() {
    // PBC_FUZZ_CASES and PBC_FUZZ_SEED extend or replay this, as for `fuzz`.
    let (count, seed) = fuzz::fuzz_range(30, 0x9e37_79b9_7f4a_7c15);
    let statements = directly_convertible_statements();
    let mut state = seed.max(1);
    let mut next = |bound: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % bound as u64) as usize
    };
    for _ in 0..count {
        let length = 4 + next(12);
        let body: Vec<_> = (0..length)
            .map(|_| statements[next(statements.len())].as_str())
            .collect();
        let input = Circuit::from_qasm(&format!("qreg q[3]; {};", body.join("; "))).unwrap();
        assert_pbc_matches_native(&input, &ALL_BOUNDS);
    }
}

#[test]
fn non_clifford_gates_name_the_decomposition_they_need() {
    for (gate, pass) in [
        ("ch q[0],q[1]", "DecomposeRotations"),
        ("rz(0.3) q[0]", "DecomposeRotations"),
        ("crz(pi/2) q[0],q[1]", "DecomposeRotations"),
        ("cp(pi/2) q[0],q[1]", "DecomposeRotations"),
        ("crx(pi/4) q[0],q[1]", "DecomposeRotations"),
        ("cswap q[0],q[1],q[2]", "DecomposeToffoli"),
    ] {
        let input = Circuit::from_qasm(&format!("qreg q[3]; h q[0]; {gate};")).unwrap();
        assert!(needs_decomposition(&input.gates[1]), "{gate}");
        let error = to_pbc(&input, None).unwrap_err();
        assert!(
            matches!(&error, PbcError::InvalidInput { index: 1, .. }),
            "{gate}: {error}"
        );
        assert!(error.to_string().contains(pass), "{gate}: {error}");
    }
    for gate in [
        Gate::y(0),
        Gate::sx(0),
        Gate::swap(0, 1),
        Gate::cy {
            control: 0,
            target: 1,
        },
        Gate::reset(0),
    ] {
        assert!(!needs_decomposition(&gate), "{gate}");
    }
}

#[test]
fn lowered_gates_still_validate_operands_even_when_they_lower_to_nothing() {
    let zero = crate::angle::Angle::pi_fraction(0, 1).unwrap();
    let pi = crate::angle::Angle::pi_fraction(1, 1).unwrap();
    for gate in [
        Gate::y(3),
        Gate::sx(3),
        Gate::swap(0, 3),
        Gate::cy {
            control: 3,
            target: 0,
        },
        Gate::rz(zero.clone(), 3),
        Gate::p(zero.clone(), 3),
        Gate::rx(zero.clone(), 3),
        Gate::ry(zero.clone(), 3),
        Gate::crz {
            theta: zero.clone(),
            control: 0,
            target: 3,
        },
        Gate::cp {
            lambda: pi,
            control: 3,
            target: 0,
        },
    ] {
        let mut input = Circuit::new(3);
        input.apply(Gate::h(0));
        input.apply(gate.clone());
        assert_eq!(
            to_pbc(&input, None).unwrap_err(),
            PbcError::InvalidInput {
                index: 1,
                cause: Box::new(PbcError::QubitOutOfRange(3)),
            },
            "{gate}"
        );
    }
}
