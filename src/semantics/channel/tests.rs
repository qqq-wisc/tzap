use super::*;
use crate::pbc::{PauliAngle, to_pbc};

mod fuzz;

fn circuit(n: usize, cbits: usize, gates: Vec<Gate>) -> Circuit {
    Circuit {
        num_qubits: n,
        num_cbits: cbits,
        gates,
    }
}
fn measure(q: u32, c: u32) -> Gate {
    Gate::measure { qubit: q, cbit: c }
}
fn cx(c: u32, t: u32) -> Gate {
    Gate::cnot {
        control: c,
        target: t,
    }
}
fn gate_channel(c: &Circuit, initial: &[bool]) -> Channel {
    circuit_channel(c, initial, ChannelLimits::default()).unwrap()
}
fn pbc(c: &PbcCircuit, initial: &[bool]) -> Channel {
    pbc_channel(c, initial, ChannelLimits::default()).unwrap()
}
fn check(c: &Circuit, initial: &[bool]) -> Channel {
    let expected = gate_channel(c, initial);
    let converted = to_pbc(c).unwrap();
    let actual = pbc(&converted, initial);
    assert_eq!(
        expected.compare(&actual),
        Ok(()),
        "input: {c}\ninitial: {initial:?}"
    );
    assert_trace_preserving(&actual);
    actual
}

/// Partial trace over output of the (unnormalized) Choi matrix is identity.
fn assert_trace_preserving(channel: &Channel) {
    let d = channel.quantum_dim;
    for a in 0..d {
        for b in 0..d {
            let mut sum = Scalar::zero();
            for block in channel.blocks.values() {
                for output in 0..d {
                    sum = sum.add(block.get(output * d + a, output * d + b));
                }
            }
            assert_eq!(sum, Scalar::integer(i64::from(a == b)));
        }
    }
}

/// Output density block for a computational-basis input, extracted directly
/// from Choi entries. Used for analytic probability/post-state checks only;
/// conversion equivalence itself compares the entire channel.
fn basis_output(channel: &Channel, classical: &[bool], input: usize) -> Matrix {
    let d = channel.quantum_dim;
    let mut rho = Matrix::zero(d);
    if let Some(block) = channel.blocks.get(classical) {
        for r in 0..d {
            for c in 0..d {
                rho.set(r, c, block.get(r * d + input, c * d + input).clone());
            }
        }
    }
    rho
}

#[test]
fn analytic_identity_and_measurement_choi_matrices() {
    let identity = check(&circuit(1, 0, vec![]), &[]);
    let mut expected = Matrix::zero(4);
    for r in [0, 3] {
        for c in [0, 3] {
            expected.set(r, c, Scalar::integer(1));
        }
    }
    assert_eq!(identity.blocks[&vec![]], expected);

    let measured = check(&circuit(1, 1, vec![measure(0, 0)]), &[false]);
    for (bit, index) in [(false, 0), (true, 3)] {
        let mut expected = Matrix::zero(4);
        expected.set(index, index, Scalar::integer(1));
        assert_eq!(measured.blocks[&vec![bit]], expected);
    }
}

#[test]
fn identity_is_not_hidden_measurement_even_though_basis_inputs_match() {
    let identity = gate_channel(&circuit(1, 0, vec![]), &[]);
    let mut hidden = PbcCircuit::new(1, 0);
    let z = hidden.z(0).unwrap();
    hidden.measure(z, None).unwrap();
    let dephased = pbc(&hidden, &[]);
    for input in 0..2 {
        assert_eq!(
            basis_output(&identity, &[], input),
            basis_output(&dephased, &[], input)
        );
    }
    assert!(identity.compare(&dephased).is_err()); // Lost off-diagonal coherence.
}

#[test]
fn bell_terminal_measurement_has_correlated_half_probability_outputs() {
    let output = check(
        &circuit(
            2,
            2,
            vec![Gate::h(0), cx(0, 1), measure(0, 0), measure(1, 1)],
        ),
        &[false, false],
    );
    for (bits, index) in [(vec![false, false], 0), (vec![true, true], 3)] {
        let mut expected = Matrix::zero(4);
        expected.set(index, index, Scalar::integer(1).half());
        assert_eq!(basis_output(&output, &bits, 0), expected);
    }
    for bits in [vec![false, true], vec![true, false]] {
        assert_eq!(basis_output(&output, &bits, 0), Matrix::zero(4));
    }
}

#[test]
fn ghz_terminal_measurements_have_only_two_outcomes_for_zero_input() {
    let output = check(
        &circuit(
            3,
            3,
            vec![
                Gate::h(0),
                cx(0, 1),
                cx(1, 2),
                measure(0, 0),
                measure(1, 1),
                measure(2, 2),
            ],
        ),
        &[false, false, false],
    );
    for bits in 0..8 {
        let classical: Vec<_> = (0..3).map(|q| bits & (1 << (2 - q)) != 0).collect();
        let mut expected = Matrix::zero(8);
        if bits == 0 || bits == 7 {
            expected.set(bits, bits, Scalar::integer(1).half());
        }
        assert_eq!(basis_output(&output, &classical, 0), expected);
    }
}

#[test]
fn four_qubit_partial_terminal_measurement_retains_unmeasured_outputs() {
    check(
        &circuit(
            4,
            1,
            vec![
                Gate::h(0),
                cx(0, 3),
                Gate::t(3),
                Gate::h(2),
                Gate::ccx {
                    control1: 3,
                    control2: 2,
                    target: 1,
                },
                measure(3, 0),
            ],
        ),
        &[true],
    );
}

#[test]
fn terminal_measurement_after_all_single_qubit_gate_pairs() {
    let gates = [
        Gate::h(0),
        Gate::x(0),
        Gate::z(0),
        Gate::s(0),
        Gate::sdg(0),
        Gate::t(0),
        Gate::tdg(0),
    ];
    for a in &gates {
        for b in &gates {
            for initial in [false, true] {
                check(
                    &circuit(1, 1, vec![a.clone(), b.clone(), measure(0, 0)]),
                    &[initial],
                );
            }
        }
    }
}

#[test]
fn terminal_measurement_subsets_orders_and_classical_destinations() {
    let prefix = vec![Gate::h(0), Gate::s(0), cx(0, 1), Gate::tdg(1), Gate::h(1)];
    for wires in [vec![0], vec![1], vec![0, 1], vec![1, 0]] {
        for swap in [false, true] {
            let mut gates = prefix.clone();
            for q in &wires {
                gates.push(measure(*q, if swap { 1 - q } else { *q }));
            }
            // c2 is untouched and must stay true.
            check(&circuit(2, 3, gates), &[false, true, true]);
        }
    }
}

#[test]
fn terminal_native_ccx_ccz_measurements_cover_all_operand_permutations() {
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
                    check(
                        &circuit(
                            3,
                            3,
                            vec![
                                Gate::h(0),
                                Gate::s(0),
                                Gate::h(1),
                                cx(1, 2),
                                gate,
                                Gate::tdg(t),
                                Gate::h(t),
                                measure(2, 0),
                                measure(0, 2),
                                measure(1, 1),
                            ],
                        ),
                        &[true, false, true],
                    );
                }
            }
        }
    }
}

#[test]
fn repeated_terminal_measurement_sums_different_histories() {
    let once = check(&circuit(1, 1, vec![measure(0, 0)]), &[true]);
    let twice = check(&circuit(1, 1, vec![measure(0, 0), measure(0, 0)]), &[true]);
    assert_eq!(once.compare(&twice), Ok(()));
    let overwritten = check(&circuit(2, 1, vec![measure(0, 0), measure(1, 0)]), &[false]);
    let mut hidden = PbcCircuit::new(2, 1);
    let z0 = hidden.z(0).unwrap();
    let z1 = hidden.z(1).unwrap();
    hidden.measure(z0, None).unwrap();
    hidden.measure(z1, Some(0)).unwrap();
    assert_eq!(overwritten.compare(&pbc(&hidden, &[false])), Ok(()));
}

#[test]
fn untouched_classical_bits_preserve_all_initial_stores() {
    let c = circuit(1, 3, vec![Gate::h(0), measure(0, 1)]);
    for bits in 0..8 {
        let initial: Vec<_> = (0..3).map(|bit| bits & (1 << bit) != 0).collect();
        let output = check(&c, &initial);
        for key in output.blocks.keys() {
            assert_eq!(key[0], initial[0]);
            assert_eq!(key[2], initial[2]);
        }
    }
}

#[test]
fn identity_measurements_and_global_phases_are_invisible() {
    for negative in [false, true] {
        let mut c = PbcCircuit::new(1, 0);
        let axis = if negative {
            c.identity().negated()
        } else {
            c.identity()
        };
        c.rotate(c.identity().negated(), PauliAngle::new(1))
            .unwrap();
        c.measure(axis, None).unwrap();
        let actual = pbc(&c, &[]);
        let identity = gate_channel(&circuit(1, 0, vec![]), &[]);
        assert_eq!(actual.compare(&identity), Ok(()));
        assert_eq!(actual.blocks.len(), 1);
    }
}

#[test]
fn detect_wrong_measurement_sign_destination_and_missing_suffix() {
    let input = circuit(1, 1, vec![Gate::x(0), measure(0, 0)]);
    let expected = check(&input, &[false]);
    let mut wrong = PbcCircuit::new(1, 1);
    let z = wrong.z(0).unwrap();
    wrong.measure(z, Some(0)).unwrap();
    wrong.push_output_clifford(Gate::x(0)).unwrap();
    assert!(expected.compare(&pbc(&wrong, &[false])).is_err());

    let expected = check(&circuit(1, 2, vec![measure(0, 0)]), &[false, false]);
    let wrong = check(&circuit(1, 2, vec![measure(0, 1)]), &[false, false]);
    assert!(matches!(
        expected.compare(&wrong),
        Err(Difference::Entry { .. })
    ));

    let expected = check(&circuit(1, 1, vec![Gate::h(0), measure(0, 0)]), &[false]);
    let mut wrong = PbcCircuit::new(1, 1);
    let x = wrong.x(0).unwrap();
    wrong.measure(x, Some(0)).unwrap();
    assert!(expected.compare(&pbc(&wrong, &[false])).is_err());
}

#[test]
fn phase_before_z_measurement_is_unobservable() {
    let plain = check(&circuit(1, 1, vec![measure(0, 0)]), &[false]);
    let phase = check(&circuit(1, 1, vec![Gate::t(0), measure(0, 0)]), &[false]);
    assert_eq!(plain.compare(&phase), Ok(()));
}

#[test]
fn reject_bad_initial_stores_operands_and_resource_limits() {
    let input = circuit(1, 1, vec![measure(0, 0)]);
    let converted = to_pbc(&input).unwrap();
    let limits = ChannelLimits::default();
    assert_eq!(
        circuit_channel(&input, &[], limits),
        Err(Error::InvalidInitialStore)
    );
    assert_eq!(
        pbc_channel(&converted, &[], limits),
        Err(Error::InvalidInitialStore)
    );
    for limits in [
        ChannelLimits {
            max_branches: 1,
            ..limits
        },
        ChannelLimits {
            max_classical_bits: 0,
            ..limits
        },
        ChannelLimits {
            matrices: Limits {
                max_qubits: 0,
                ..limits.matrices
            },
            ..limits
        },
        ChannelLimits {
            matrices: Limits {
                max_matrix_cells: 0,
                ..limits.matrices
            },
            ..limits
        },
        ChannelLimits {
            matrices: Limits {
                max_multiply_terms: 0,
                ..limits.matrices
            },
            ..limits
        },
        ChannelLimits {
            matrices: Limits {
                max_operations: 0,
                ..limits.matrices
            },
            ..limits
        },
    ] {
        assert_eq!(
            circuit_channel(&input, &[false], limits),
            Err(Error::LimitExceeded)
        );
        assert_eq!(
            pbc_channel(&converted, &[false], limits),
            Err(Error::LimitExceeded)
        );
    }
    for gate in [measure(1, 0), measure(0, 1), cx(0, 0)] {
        assert_eq!(
            circuit_channel(&circuit(1, 1, vec![gate]), &[false], limits),
            Err(Error::InvalidOperand { index: 0 })
        );
    }
    assert_eq!(
        circuit_channel(&circuit(1, 0, vec![Gate::rz(0.1, 0)]), &[], limits),
        Err(Error::UnsupportedOperation { index: 0 })
    );
    let huge = circuit(usize::MAX, 0, vec![]);
    assert_eq!(
        circuit_channel(&huge, &[], limits),
        Err(Error::LimitExceeded)
    );
    let too_many = circuit(1, 1, vec![measure(0, 0); usize::BITS as usize]);
    assert_eq!(
        circuit_channel(&too_many, &[false], limits),
        Err(Error::LimitExceeded)
    );
    assert_eq!(
        pbc_channel(&to_pbc(&too_many).unwrap(), &[false], limits),
        Err(Error::LimitExceeded)
    );
}

#[test]
fn zero_qubit_empty_channels_preserve_classical_state() {
    let output = check(&circuit(0, 2, vec![]), &[true, false]);
    assert_eq!(output.blocks[&vec![true, false]], Matrix::identity(1));
    let other = gate_channel(&circuit(1, 2, vec![]), &[true, false]);
    assert_eq!(output.compare(&other), Err(Difference::Dimensions));
}

#[test]
fn resets_and_post_measurement_gates_are_rejected() {
    let limits = ChannelLimits::default();
    for gates in [vec![Gate::reset(0)], vec![Gate::h(0), Gate::reset(0)]] {
        let index = gates.len() - 1;
        assert_eq!(
            circuit_channel(&circuit(1, 0, gates), &[], limits),
            Err(Error::UnsupportedOperation { index })
        );
    }
    for gate in [Gate::h(0), Gate::x(0), Gate::t(0), Gate::sdg(0)] {
        assert_eq!(
            circuit_channel(&circuit(1, 1, vec![measure(0, 0), gate]), &[false], limits),
            Err(Error::GateAfterMeasurement { index: 1 })
        );
    }
    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    c.measure(x, None).unwrap();
    c.rotate(x, PauliAngle::new(1)).unwrap();
    assert!(pbc_channel(&c, &[], limits).is_ok());

    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    let m = c.measure(x, None).unwrap();
    c.conditional_pauli(x, m).unwrap();
    assert_eq!(
        pbc_channel(&c, &[], limits),
        Err(Error::UnsupportedOperation { index: 1 })
    );
}

#[test]
fn partial_measurement_suffix_preserves_unmeasured_qubits() {
    let input = circuit(2, 1, vec![Gate::h(1), measure(0, 0)]);
    let converted = to_pbc(&input).unwrap();
    assert_eq!(converted.output_cliffords(), &[Gate::h(1)]);
    let expected = check(&input, &[false]);
    let mut wrong = PbcCircuit::new(2, 1);
    let z = wrong.z(0).unwrap();
    wrong.measure(z, Some(0)).unwrap();
    assert!(expected.compare(&pbc(&wrong, &[false])).is_err());
}

#[test]
fn terminal_pbc_measurements_keep_unnormalized_quarter_probabilities() {
    let mut c = PbcCircuit::new(1, 2);
    let x = c.x(0).unwrap();
    let z = c.z(0).unwrap();
    c.measure(x, Some(0)).unwrap();
    c.measure(z, Some(1)).unwrap();
    let output = pbc(&c, &[false, false]);
    for a in [false, true] {
        for b in [false, true] {
            let mut expected = Matrix::zero(2);
            expected.set(
                usize::from(b),
                usize::from(b),
                Scalar::integer(1).half().half(),
            );
            assert_eq!(basis_output(&output, &[a, b], 0), expected);
        }
    }
    assert_trace_preserving(&output);
}

#[test]
fn unitary_channels_agree_with_unitary_matrices() {
    let c = circuit(2, 0, vec![Gate::h(0), cx(0, 1), Gate::t(1), Gate::h(1)]);
    let u = circuit_unitary(&c, Limits::default()).unwrap();
    let expected = finish(
        4,
        0,
        vec![Branch {
            kraus: u,
            classical: vec![],
        }],
    );
    assert_eq!(check(&c, &[]).compare(&expected), Ok(()));
}
