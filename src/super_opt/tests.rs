use super::matrix::{IDENTITY_TOLERANCE, UnitaryMatrix, unitary_fingerprint};
use super::matrix_cache::COMPACT_KEY_MAX_GATES;
use super::murm::{LibraryGate, Murm, library_circuit_matrix, library_gates, shared_murm};
use super::*;
use crate::circuit::{GateKind, GateSet};
use crate::unitary::C as OracleScalar;

struct TestRng(u64);

impl TestRng {
    fn next(&mut self, upper: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 as usize % upper
    }

    /// [`TestRng::next`] as a qubit index.
    fn qubit(&mut self, upper: usize) -> Qubit {
        self.next(upper) as Qubit
    }
}

fn removed_indices(result: &SuperOptResult) -> Vec<Vec<usize>> {
    let mut removed: Vec<_> = result
        .rewrites
        .iter()
        .filter(|rewrite| rewrite.replacement.is_empty())
        .map(|rewrite| rewrite.gate_indices.clone())
        .collect();
    removed.sort();
    removed
}

fn union_qubits(left: &[Qubit], right: &[Qubit]) -> Vec<Qubit> {
    let mut union = Vec::with_capacity(left.len() + right.len());
    union.extend_from_slice(left);
    union.extend_from_slice(right);
    union.sort_unstable();
    union.dedup();
    union
}

fn naive_matrix(circuit: &Circuit, gate_indices: &[usize], support: &[Qubit]) -> UnitaryMatrix {
    let mut matrix = UnitaryMatrix::identity(support.len()).unwrap();
    for &gate_index in gate_indices {
        matrix
            .apply_gate_left(&circuit.gates[gate_index], support)
            .unwrap();
    }
    matrix
}

fn assert_matrix_equal(actual: &UnitaryMatrix, expected: &UnitaryMatrix) {
    assert_eq!(actual, expected);
}

fn assert_integer_entry(matrix: &UnitaryMatrix, row: usize, column: usize, value: i8) {
    assert_eq!(
        matrix.entry_coefficients(row, column),
        ([value, 0, 0, 0], 0)
    );
}

#[derive(Clone)]
struct ChannelBranch {
    amplitudes: Vec<OracleScalar>,
    cbits: Vec<bool>,
}

fn apply_unitary_to_state(
    gate: &Gate,
    num_qubits: usize,
    state: &[OracleScalar],
) -> Vec<OracleScalar> {
    let mut circuit = Circuit::new(num_qubits);
    circuit.apply(gate.clone());
    let matrix = crate::unitary::circuit_unitary(&circuit);
    (0..state.len())
        .map(|row| {
            (0..state.len()).fold(OracleScalar::ZERO, |sum, column| {
                sum + matrix[row][column] * state[column]
            })
        })
        .collect()
}

/// Independent small-circuit reference channel, represented as one density matrix
/// per final classical-bit valuation. Measurement and reset split branches;
/// unitary gates evolve each branch independently.
fn reference_channel(
    circuit: &Circuit,
) -> std::collections::BTreeMap<Vec<bool>, Vec<OracleScalar>> {
    let dim = 1 << circuit.num_qubits;
    let mut initial = vec![OracleScalar::ZERO; dim];
    initial[0] = OracleScalar::ONE;
    let mut branches = vec![ChannelBranch {
        amplitudes: initial,
        cbits: vec![false; circuit.num_cbits],
    }];

    for gate in &circuit.gates {
        let mut next = Vec::new();
        for branch in branches {
            match gate {
                Gate::measure { qubit, cbit } => {
                    let bit = 1 << (circuit.num_qubits - 1 - *qubit as usize);
                    for outcome in [false, true] {
                        let mut measured = branch.clone();
                        measured.cbits[*cbit as usize] = outcome;
                        for (basis, amplitude) in measured.amplitudes.iter_mut().enumerate() {
                            if (basis & bit != 0) != outcome {
                                *amplitude = OracleScalar::ZERO;
                            }
                        }
                        if measured
                            .amplitudes
                            .iter()
                            .any(|value| value.norm_sq() > 1e-24)
                        {
                            next.push(measured);
                        }
                    }
                }
                Gate::reset(qubit) => {
                    let bit = 1 << (circuit.num_qubits - 1 - *qubit as usize);
                    for outcome in [false, true] {
                        let mut reset = vec![OracleScalar::ZERO; dim];
                        for (basis, &amplitude) in branch.amplitudes.iter().enumerate() {
                            if (basis & bit != 0) == outcome {
                                reset[if outcome { basis ^ bit } else { basis }] = amplitude;
                            }
                        }
                        if reset.iter().any(|value| value.norm_sq() > 1e-24) {
                            next.push(ChannelBranch {
                                amplitudes: reset,
                                cbits: branch.cbits.clone(),
                            });
                        }
                    }
                }
                _ => next.push(ChannelBranch {
                    amplitudes: apply_unitary_to_state(
                        gate,
                        circuit.num_qubits,
                        &branch.amplitudes,
                    ),
                    cbits: branch.cbits,
                }),
            }
        }
        branches = next;
    }

    let mut channel = std::collections::BTreeMap::new();
    for branch in branches {
        let density = channel
            .entry(branch.cbits)
            .or_insert_with(|| vec![OracleScalar::ZERO; dim * dim]);
        for row in 0..dim {
            for column in 0..dim {
                let conjugate = branch.amplitudes[column].conj();
                density[row * dim + column] =
                    density[row * dim + column] + branch.amplitudes[row] * conjugate;
            }
        }
    }
    channel
}

fn assert_channels_equivalent(actual: &Circuit, expected: &Circuit) {
    assert_eq!(actual.num_qubits, expected.num_qubits);
    assert_eq!(actual.num_cbits, expected.num_cbits);
    let actual = reference_channel(actual);
    let expected = reference_channel(expected);
    assert_eq!(
        actual.keys().collect::<Vec<_>>(),
        expected.keys().collect::<Vec<_>>()
    );
    for (cbits, expected_density) in expected {
        let actual_density = &actual[&cbits];
        for (index, (&left, &right)) in actual_density.iter().zip(&expected_density).enumerate() {
            let delta = left - right;
            assert!(
                delta.norm_sq() <= 1e-20,
                "channel differs for classical state {cbits:?} at density entry {index}: {left:?} != {right:?}"
            );
        }
    }
}

fn naive_windows(
    circuit: &Circuit,
    max_qubits: usize,
    max_gates: usize,
) -> Vec<(Vec<usize>, Vec<Qubit>)> {
    let mut result = Vec::new();
    for anchor in 0..circuit.gates.len() {
        let mut previous_indices = Vec::new();
        for end in anchor..circuit.gates.len() {
            let mut indices = vec![anchor];
            let mut qubits: Vec<Qubit> = unique_qubits(&circuit.gates[anchor]).to_vec();
            loop {
                let mut changed = false;
                for gate_index in anchor..=end {
                    if indices.binary_search(&gate_index).is_ok() {
                        continue;
                    }
                    let gate_qubits = unique_qubits(&circuit.gates[gate_index]);
                    if gate_qubits
                        .iter()
                        .any(|qubit| qubits.binary_search(qubit).is_ok())
                    {
                        let position = indices.binary_search(&gate_index).unwrap_err();
                        indices.insert(position, gate_index);
                        qubits = union_qubits(&qubits, &gate_qubits);
                        changed = true;
                    }
                }
                if !changed {
                    break;
                }
            }

            if indices.len() > max_gates || qubits.len() > max_qubits {
                break;
            }
            if indices != previous_indices {
                previous_indices = indices.clone();
                if !indices
                    .iter()
                    .any(|&index| matches!(circuit.gates[index], Gate::rz(..)))
                {
                    result.push((indices.clone(), qubits));
                }
            }
            if indices.len() == max_gates {
                break;
            }
        }
    }
    result.sort_by(|left, right| left.0.cmp(&right.0));
    result
}

#[test]
fn disjoint_prefix_can_be_connected_by_later_gate() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(1));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });

    let result = SuperOpt::analyzer(2, 3).run(&circuit).unwrap();
    let window = result
        .subcircuits
        .iter()
        .find(|window| window.gate_indices == [0, 1, 2])
        .unwrap();
    assert_eq!(window.qubits, vec![0, 1]);
    let expected = naive_matrix(&circuit, &[0, 1, 2], &[0, 1]);
    assert_matrix_equal(&window.matrix, &expected);
}

#[test]
fn emits_one_window_per_anchor_not_all_combinations() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::s(0));
    circuit.apply(Gate::t(0));

    let result = SuperOpt::analyzer(1, 2).run(&circuit).unwrap();
    let indices: Vec<_> = result
        .subcircuits
        .iter()
        .map(|window| window.gate_indices.clone())
        .collect();
    assert_eq!(
        indices,
        vec![
            vec![0],
            vec![0, 1],
            vec![1],
            vec![1, 2],
            vec![2],
            vec![2, 3],
            vec![3],
        ]
    );
}

#[test]
fn disconnected_completed_window_is_not_emitted() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(1));

    let result = SuperOpt::analyzer(2, 2).run(&circuit).unwrap();
    assert!(
        !result
            .subcircuits
            .iter()
            .any(|window| window.gate_indices == [0, 1])
    );
}

#[test]
fn unrelated_gates_are_skipped_until_anchor_reconnects() {
    let mut circuit = Circuit::new(3);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(2));
    circuit.apply(Gate::z(2));
    circuit.apply(Gate::h(0));

    let result = SuperOpt::analyzer(1, 2).run(&circuit).unwrap();
    assert!(
        result
            .subcircuits
            .iter()
            .any(|window| window.gate_indices == [0, 3])
    );
}

#[test]
fn bridge_pulls_in_entire_intervening_component() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(1));
    circuit.apply(Gate::x(1));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });

    let three_gate = SuperOpt::analyzer(2, 3).run(&circuit).unwrap();
    assert!(
        !three_gate
            .subcircuits
            .iter()
            .any(|window| window.gate_indices == [0, 1, 3])
    );
    assert!(
        three_gate
            .subcircuits
            .iter()
            .any(|window| window.gate_indices == [1, 2, 3])
    );

    let four_gate = SuperOpt::analyzer(2, 4).run(&circuit).unwrap();
    assert!(
        four_gate
            .subcircuits
            .iter()
            .any(|window| window.gate_indices == [0, 1, 2, 3])
    );
}

#[test]
fn over_width_partial_window_is_dropped() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    circuit.apply(Gate::t(0));

    let result = SuperOpt::analyzer(1, 2).run(&circuit).unwrap();
    assert!(
        result
            .subcircuits
            .iter()
            .all(|window| window.gate_indices.len() == 1)
    );
}

#[test]
fn canonical_cache_reuses_shifted_windows() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::h(1));
    circuit.apply(Gate::x(1));

    let result = SuperOpt::analyzer(1, 2).run(&circuit).unwrap();
    assert_eq!(result.subcircuits.len(), 6);
    assert_eq!(result.cache_hits, 3);
    assert_eq!(result.cache_misses, 3);
    assert!(Arc::ptr_eq(
        &result.subcircuits[1].matrix,
        &result.subcircuits[4].matrix
    ));
}

#[test]
fn compact_cache_boundary_falls_back_without_changing_matrices() {
    let pattern = [
        Gate::h(0),
        Gate::t(0),
        Gate::s(0),
        Gate::x(0),
        Gate::tdg(0),
        Gate::z(0),
        Gate::sdg(0),
    ];
    let sequence: Vec<Gate> = pattern
        .iter()
        .cycle()
        .take(COMPACT_KEY_MAX_GATES + 1)
        .cloned()
        .collect();
    let mut circuit = Circuit::new(2);
    for gate in &sequence {
        circuit.apply(gate.clone());
    }
    for gate in &sequence {
        circuit.apply(gate.map_qubits(|_| 1));
    }

    let first_indices: Vec<_> = (0..sequence.len()).collect();
    let second_indices: Vec<_> = (sequence.len()..2 * sequence.len()).collect();
    assert!(
        compact_normalized_key(&circuit, &first_indices[..COMPACT_KEY_MAX_GATES], &[0]).is_some()
    );
    assert!(compact_normalized_key(&circuit, &first_indices, &[0]).is_none());

    let result = SuperOpt::analyzer(1, sequence.len()).run(&circuit).unwrap();
    let first = result
        .subcircuits
        .iter()
        .find(|window| window.gate_indices == first_indices)
        .unwrap();
    let second = result
        .subcircuits
        .iter()
        .find(|window| window.gate_indices == second_indices)
        .unwrap();
    assert!(Arc::ptr_eq(&first.matrix, &second.matrix));
    assert_matrix_equal(&first.matrix, &naive_matrix(&circuit, &first_indices, &[0]));
    assert_matrix_equal(
        &second.matrix,
        &naive_matrix(&circuit, &second_indices, &[1]),
    );
}

#[test]
fn compact_key_rejects_supports_wider_than_operand_encoding() {
    let mut circuit = Circuit::new(5);
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    assert!(compact_normalized_key(&circuit, &[0], &[0, 1, 2, 3]).is_some());
    assert!(compact_normalized_key(&circuit, &[0], &[0, 1, 2, 3, 4]).is_none());
}

fn murm(max_qubits: usize, max_gates: usize) -> Arc<Murm> {
    Arc::new(
        Murm::build(MurmConfig {
            max_qubits,
            max_gates,
            max_entries_per_qubit: 20_000,
            basis: BASE_GATE_SET,
        })
        .unwrap(),
    )
}

fn single_qubit_matrix(gates: &[Gate]) -> UnitaryMatrix {
    let mut matrix = UnitaryMatrix::identity(1).unwrap();
    for gate in gates {
        matrix.apply_gate_left(gate, &[0]).unwrap();
    }
    matrix
}

#[test]
fn qubit_zero_is_the_most_significant_basis_bit() {
    let mut matrix = UnitaryMatrix::identity(2).unwrap();
    matrix.apply_gate_left(&Gate::x(0), &[0, 1]).unwrap();
    // |00> -> |10>: column 0 maps to row 2.
    assert_integer_entry(&matrix, 2, 0, 1);
    assert_integer_entry(&matrix, 0, 0, 0);
}

#[test]
fn cnot_matrix_matches_truth_table() {
    let mut matrix = UnitaryMatrix::identity(2).unwrap();
    matrix
        .apply_gate_left(
            &Gate::cnot {
                control: 0,
                target: 1,
            },
            &[0, 1],
        )
        .unwrap();
    for (input, output) in [(0, 0), (1, 1), (2, 3), (3, 2)] {
        assert_integer_entry(&matrix, output, input, 1);
    }
}

#[test]
fn directional_gate_matrices_cover_every_target_position() {
    for (control, target) in [(0, 1), (1, 0)] {
        let mut matrix = UnitaryMatrix::identity(2).unwrap();
        matrix
            .apply_gate_left(&Gate::cnot { control, target }, &[0, 1])
            .unwrap();
        let control_bit = 1 << (1 - control);
        let target_bit = 1 << (1 - target);
        for input in 0..4 {
            let output = if input & control_bit == 0 {
                input
            } else {
                input ^ target_bit
            };
            assert_integer_entry(&matrix, output, input, 1);
        }
    }

    for target in 0..3 {
        let controls: Vec<_> = (0..3).filter(|&qubit| qubit != target).collect();
        let mut matrix = UnitaryMatrix::identity(3).unwrap();
        matrix
            .apply_gate_left(
                &Gate::ccx {
                    control1: controls[0],
                    control2: controls[1],
                    target,
                },
                &[0, 1, 2],
            )
            .unwrap();
        let control_mask = (1 << (2 - controls[0])) | (1 << (2 - controls[1]));
        let target_bit = 1 << (2 - target);
        for input in 0..8 {
            let output = if input & control_mask == control_mask {
                input ^ target_bit
            } else {
                input
            };
            assert_integer_entry(&matrix, output, input, 1);
        }
    }
}

#[test]
fn cz_matrix_is_symmetric_in_its_qubits() {
    let mut forward = UnitaryMatrix::identity(2).unwrap();
    forward
        .apply_gate_left(
            &Gate::cz {
                control: 0,
                target: 1,
            },
            &[0, 1],
        )
        .unwrap();
    let mut reversed = UnitaryMatrix::identity(2).unwrap();
    reversed
        .apply_gate_left(
            &Gate::cz {
                control: 1,
                target: 0,
            },
            &[0, 1],
        )
        .unwrap();
    assert_matrix_equal(&forward, &reversed);
    assert_integer_entry(&forward, 3, 3, -1);
    assert_integer_entry(&forward, 0, 0, 1);
}

#[test]
fn ccx_matrix_swaps_the_last_two_basis_states() {
    let mut matrix = UnitaryMatrix::identity(3).unwrap();
    matrix
        .apply_gate_left(
            &Gate::ccx {
                control1: 0,
                control2: 1,
                target: 2,
            },
            &[0, 1, 2],
        )
        .unwrap();
    assert_integer_entry(&matrix, 7, 6, 1);
    assert_integer_entry(&matrix, 6, 7, 1);
    for basis in 0..6 {
        assert_integer_entry(&matrix, basis, basis, 1);
    }
}

#[test]
fn ccz_matrix_negates_only_the_all_ones_state() {
    let mut matrix = UnitaryMatrix::identity(3).unwrap();
    matrix
        .apply_gate_left(
            &Gate::ccz {
                control1: 2,
                control2: 0,
                target: 1,
            },
            &[0, 1, 2],
        )
        .unwrap();

    for basis in 0..7 {
        assert_integer_entry(&matrix, basis, basis, 1);
    }
    assert_integer_entry(&matrix, 7, 7, -1);
}

#[test]
fn compact_keys_distinguish_ccx_from_ccz() {
    let mut circuit = Circuit::new(3);
    circuit.apply(Gate::ccx {
        control1: 0,
        control2: 1,
        target: 2,
    });
    circuit.apply(Gate::ccz {
        control1: 0,
        control2: 1,
        target: 2,
    });

    let ccx_key = compact_normalized_key(&circuit, &[0], &[0, 1, 2]).unwrap();
    let ccz_key = compact_normalized_key(&circuit, &[1], &[0, 1, 2]).unwrap();

    assert_ne!(ccx_key, ccz_key);
}

#[test]
fn compact_keys_canonicalize_symmetric_native_operands() {
    let key = |gate| {
        let mut circuit = Circuit::new(3);
        circuit.apply(gate);
        compact_normalized_key(&circuit, &[0], &[0, 1, 2]).unwrap()
    };

    assert_eq!(
        key(Gate::cz {
            control: 0,
            target: 2,
        }),
        key(Gate::cz {
            control: 2,
            target: 0,
        })
    );
    assert_eq!(
        key(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        }),
        key(Gate::ccx {
            control1: 1,
            control2: 0,
            target: 2,
        })
    );
    let canonical = key(Gate::ccz {
        control1: 0,
        control2: 1,
        target: 2,
    });
    for [a, b, c] in [[0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
        assert_eq!(
            canonical,
            key(Gate::ccz {
                control1: a,
                control2: b,
                target: c,
            })
        );
    }
}

#[test]
fn oversized_matrix_request_errors_instead_of_allocating() {
    for num_qubits in [40, 64, 200] {
        assert_eq!(
            UnitaryMatrix::identity(num_qubits).unwrap_err(),
            SuperOptError::MatrixTooLarge { num_qubits }
        );
    }
}

#[test]
fn global_phase_equivalence_accepts_phase_and_rejects_difference() {
    let x = single_qubit_matrix(&[Gate::x(0)]);
    let minus_x = single_qubit_matrix(&[Gate::z(0), Gate::x(0), Gate::z(0)]);
    let s = single_qubit_matrix(&[Gate::s(0)]);
    let sdg = single_qubit_matrix(&[Gate::sdg(0)]);
    assert!(x.equivalent_up_to_global_phase(&minus_x));
    assert!(!s.equivalent_up_to_global_phase(&sdg));
}

#[test]
fn fingerprint_ignores_global_phase() {
    let x = single_qubit_matrix(&[Gate::x(0)]);
    // Z X Z = -X: same unitary up to a global phase of -1.
    let minus_x = single_qubit_matrix(&[Gate::z(0), Gate::x(0), Gate::z(0)]);
    assert_eq!(unitary_fingerprint(&x), unitary_fingerprint(&minus_x));

    // T X T X = omega I, exercising a non-real eighth-root global phase.
    let identity = single_qubit_matrix(&[]);
    let omega_identity = single_qubit_matrix(&[Gate::x(0), Gate::t(0), Gate::x(0), Gate::t(0)]);
    assert!(identity.equivalent_up_to_global_phase(&omega_identity));
    assert_eq!(
        unitary_fingerprint(&identity),
        unitary_fingerprint(&omega_identity)
    );
}

#[test]
fn fingerprint_distinguishes_s_from_sdg() {
    let s = single_qubit_matrix(&[Gate::s(0)]);
    let sdg = single_qubit_matrix(&[Gate::sdg(0)]);
    assert_ne!(unitary_fingerprint(&s), unitary_fingerprint(&sdg));
}

#[test]
fn denominator_normalization_buckets_hh_with_identity() {
    let identity = single_qubit_matrix(&[]);
    let hh = single_qubit_matrix(&[Gate::h(0), Gate::h(0)]);
    assert!(identity.equivalent_up_to_global_phase(&hh));
    assert_eq!(unitary_fingerprint(&identity), unitary_fingerprint(&hh));
}

#[test]
fn library_gate_inverse_pairs() {
    assert!(LibraryGate::S(0).is_inverse_of(LibraryGate::Sdg(0)));
    assert!(LibraryGate::Tdg(1).is_inverse_of(LibraryGate::T(1)));
    assert!(LibraryGate::X(0).is_inverse_of(LibraryGate::X(0)));
    assert!(LibraryGate::Cnot(0, 1).is_inverse_of(LibraryGate::Cnot(0, 1)));
    assert!(!LibraryGate::S(0).is_inverse_of(LibraryGate::S(0)));
    assert!(!LibraryGate::T(0).is_inverse_of(LibraryGate::Tdg(1)));
    assert!(!LibraryGate::X(0).is_inverse_of(LibraryGate::X(1)));
    assert!(!LibraryGate::Cnot(0, 1).is_inverse_of(LibraryGate::Cnot(1, 0)));
}

#[test]
fn library_gate_disjointness() {
    assert!(LibraryGate::X(0).is_disjoint(LibraryGate::H(1)));
    assert!(!LibraryGate::Cnot(0, 1).is_disjoint(LibraryGate::Cnot(1, 2)));
}

#[test]
fn library_gate_counts_per_width() {
    // 7n singles + n(n-1) CNOT in the base basis.
    assert_eq!(library_gates(1, BASE_GATE_SET).len(), 7);
    assert_eq!(library_gates(2, BASE_GATE_SET).len(), 16);
    assert_eq!(library_gates(3, BASE_GATE_SET).len(), 27);
    assert_eq!(library_gates(4, BASE_GATE_SET).len(), 40);
}

#[test]
fn native_library_generators_are_canonical_and_complete() {
    let gates = library_gates(3, SUPPORTED_GATE_SET);
    assert_eq!(gates.len(), 34);
    assert_eq!(
        gates
            .iter()
            .filter(|gate| matches!(gate, LibraryGate::Cz(..)))
            .count(),
        3
    );
    assert_eq!(
        gates
            .iter()
            .filter(|gate| matches!(gate, LibraryGate::Ccx(..)))
            .count(),
        3
    );
    assert_eq!(
        gates
            .iter()
            .filter(|gate| matches!(gate, LibraryGate::Ccz(..)))
            .count(),
        1
    );
    assert!(gates.iter().all(|gate| match gate {
        LibraryGate::Cz(a, b) => a < b,
        LibraryGate::Ccx(a, b, target) => a < b && target != a && target != b,
        LibraryGate::Ccz(a, b, c) => a < b && b < c,
        _ => true,
    }));
}

#[test]
fn native_library_gates_round_trip_the_four_byte_encoding() {
    for gate in [
        LibraryGate::Cz(1, 4),
        LibraryGate::Ccx(0, 2, 1),
        LibraryGate::Ccz(0, 1, 2),
    ] {
        assert_eq!(LibraryGate::from_bytes(gate.to_bytes()), Some(gate));
    }
}

#[test]
fn native_murms_synthesize_native_representatives() {
    for (kind, gate) in [
        (
            GateKind::Cz,
            Gate::cz {
                control: 0,
                target: 1,
            },
        ),
        (
            GateKind::Ccx,
            Gate::ccx {
                control1: 0,
                control2: 1,
                target: 2,
            },
        ),
        (
            GateKind::Ccz,
            Gate::ccz {
                control1: 0,
                control2: 1,
                target: 2,
            },
        ),
    ] {
        let basis = GateSet::singleton(kind);
        let murm = Murm::build(MurmConfig::new(3, 1, 100).with_basis(basis)).unwrap();
        let mut matrix = UnitaryMatrix::identity(3).unwrap();
        matrix.apply_gate_left(&gate, &[0, 1, 2]).unwrap();
        assert_eq!(murm.synthesize(&matrix), Some(vec![gate]));
    }
}

#[test]
fn explicit_cz_murm_rewrites_h_cx_h_to_native_cz() {
    let config = MurmConfig::new(2, 1, 64).with_basis(GateSet::singleton(GateKind::Cz));
    let pass = SuperOpt::new(2, 3, config).unwrap();
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(1));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    circuit.apply(Gate::h(1));

    let result = pass.run(&circuit).unwrap().circuit;
    assert_eq!(
        result.gates,
        vec![Gate::cz {
            control: 0,
            target: 1,
        }]
    );
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &result,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn every_native_mask_builds_and_round_trips_its_own_murm() {
    let dir = tempfile::tempdir().unwrap();
    for mask in 0u8..8 {
        let optional = GateSet::from_kinds(
            OPTIONAL_GATE_KINDS
                .into_iter()
                .enumerate()
                .filter_map(|(bit, kind)| (mask & (1 << bit) != 0).then_some(kind)),
        );
        let config = MurmConfig::new(3, 1, 256).with_basis(BASE_GATE_SET.union(optional));
        let cold = Murm::build(config).unwrap();
        let path = dir.path().join(format!("murm-{mask}.bin"));
        cold.write_to_disk(&path, config).unwrap();
        let warm = Murm::read_from_disk(&path, config).unwrap();
        for width in 1..=3 {
            assert_eq!(cold.entry_count(width), warm.entry_count(width));
            assert_eq!(cold.completed_depth(width), warm.completed_depth(width));
        }
    }
}

#[test]
fn actual_murms_cover_every_single_gate_exclusion_and_mixed_basis() {
    let candidates = [
        GateKind::H,
        GateKind::X,
        GateKind::Z,
        GateKind::S,
        GateKind::Sdg,
        GateKind::T,
        GateKind::Tdg,
        GateKind::Cx,
        GateKind::Cz,
        GateKind::Ccx,
        GateKind::Ccz,
    ];
    for excluded in candidates {
        let basis = GateSet::from_kinds(
            candidates
                .into_iter()
                .filter(|candidate| *candidate != excluded),
        );
        Murm::build(MurmConfig::new(3, 1, 256).with_basis(basis)).unwrap();
    }
    for basis in [
        GateSet::from_kinds([GateKind::H, GateKind::Cz, GateKind::Ccx]),
        GateSet::from_kinds([GateKind::T, GateKind::Tdg, GateKind::Ccz]),
        GateSet::from_kinds([GateKind::Cx, GateKind::Cz]),
    ] {
        Murm::build(MurmConfig::new(3, 2, 256).with_basis(basis)).unwrap();
    }
}

#[test]
fn library_never_enumerates_toffoli_or_cz() {
    for num_qubits in 1..=4 {
        for gate in library_gates(num_qubits, BASE_GATE_SET) {
            assert!(
                !matches!(gate.to_gate(), Gate::ccx { .. } | Gate::cz { .. }),
                "library must not contain Toffoli or CZ: {gate:?}"
            );
        }
    }
}

#[test]
fn murm_does_not_synthesize_a_toffoli_representative() {
    // A Toffoli's unitary has no Clifford+T representative within the small gate
    // bound and Toffoli itself is not in the library, so it must not be found.
    let murm = murm(3, 5);
    let mut toffoli = UnitaryMatrix::identity(3).unwrap();
    toffoli
        .apply_gate_left(
            &Gate::ccx {
                control1: 0,
                control2: 1,
                target: 2,
            },
            &[0, 1, 2],
        )
        .unwrap();
    assert!(murm.synthesize(&toffoli).is_none());
}

#[test]
fn one_qubit_depth_one_murm_has_eight_distinct_entries() {
    let murm = murm(1, 1);
    // Identity plus X, H, S, Sdg, Z, T, Tdg — all distinct up to phase.
    assert_eq!(murm.entry_count(1), 8);
    assert_eq!(murm.completed_depth(1), 1);
    assert!(!murm.is_saturated(1));
}

#[test]
fn synthesize_returns_empty_circuit_for_identity() {
    let murm = murm(1, 1);
    let identity = UnitaryMatrix::identity(1).unwrap();
    assert!(murm.synthesize(&identity).unwrap().is_empty());
}

#[test]
fn synthesize_returns_none_for_unknown_width_or_depth() {
    let murm = murm(1, 1);
    let three_qubit_identity = UnitaryMatrix::identity(3).unwrap();
    assert!(murm.synthesize(&three_qubit_identity).is_none());
    // H then T is not reachable within one gate.
    let deep = single_qubit_matrix(&[Gate::h(0), Gate::t(0)]);
    assert!(murm.synthesize(&deep).is_none());
}

#[test]
fn synthesis_rejects_a_forced_fingerprint_collision() {
    let mut murm = Murm::build(MurmConfig::new(1, 1, 100)).unwrap();
    let query = single_qubit_matrix(&[Gate::s(0)]);
    let wrong_candidate = single_qubit_matrix(&[Gate::x(0)]);
    murm.inject_fingerprint_alias(&query, &wrong_candidate);

    assert!(
        murm.synthesize(&query).is_none(),
        "the exact matrix guard must reject a fingerprint hit for the wrong circuit"
    );
}

#[test]
fn hzh_synthesizes_to_x() {
    let murm = murm(1, 1);
    let hzh = single_qubit_matrix(&[Gate::h(0), Gate::z(0), Gate::h(0)]);
    let replacement = murm.synthesize(&hzh).unwrap();
    assert_eq!(replacement.len(), 1);
    assert!(matches!(replacement[0], Gate::x(0)));
}

#[test]
fn cz_unitary_synthesizes_without_emitting_cz() {
    // The CZ unitary must resolve to an H/CX-basis representative (H·CX·H),
    // never to a literal cz gate — cz is excluded from the library.
    let murm = murm(2, 3);
    let mut matrix = UnitaryMatrix::identity(2).unwrap();
    matrix
        .apply_gate_left(
            &Gate::cz {
                control: 0,
                target: 1,
            },
            &[0, 1],
        )
        .unwrap();
    let replacement = murm.synthesize(&matrix).unwrap();
    assert!(
        !replacement.iter().any(|g| matches!(g, Gate::cz { .. })),
        "synthesis must not emit cz: {replacement:?}"
    );
    assert_eq!(replacement.len(), 3, "expected H·CX·H: {replacement:?}");
}

#[test]
fn random_short_library_circuits_never_synthesize_longer() {
    let murm = murm(2, 4);
    assert!(!murm.is_saturated(2));
    let gates = library_gates(2, BASE_GATE_SET);
    let mut rng = TestRng(0x7ab1_e000_c0ff_ee00);
    for _ in 0..200 {
        let length = 1 + rng.next(4);
        let circuit: Vec<_> = (0..length).map(|_| gates[rng.next(gates.len())]).collect();
        let matrix = library_circuit_matrix(2, &circuit).unwrap().unwrap();
        let replacement = murm
            .synthesize(&matrix)
            .expect("depth-4 two-qubit MURM is complete");
        assert!(replacement.len() <= length);
    }
}

#[test]
fn empty_circuit_produces_empty_result() {
    let result = SuperOpt::analyzer(2, 4).run(&Circuit::new(2)).unwrap();
    assert!(result.subcircuits.is_empty());
    assert!(result.rewrites.is_empty());
    assert!(result.circuit.gates.is_empty());
    assert_eq!(result.cache_hits + result.cache_misses, 0);
}

#[test]
fn gate_wider_than_max_qubits_is_left_alone() {
    let mut circuit = Circuit::new(3);
    circuit.apply(Gate::ccx {
        control1: 0,
        control2: 1,
        target: 2,
    });
    let result = SuperOpt::analyzer(2, 4).run(&circuit).unwrap();
    assert!(result.subcircuits.is_empty());
    assert_eq!(result.circuit.gates.len(), 1);
}

#[test]
fn toffoli_metadata_tracks_surviving_and_removed_gates() {
    let ccx = Gate::ccx {
        control1: 0,
        control2: 1,
        target: 2,
    };

    let mut surviving = Circuit::new(3);
    surviving.apply(ccx.clone());
    let surviving = SuperOpt::analyzer(3, 2).run(&surviving).unwrap().circuit;
    assert!(surviving.has_toffoli());

    let mut cancelling = Circuit::new(3);
    cancelling.apply(ccx.clone());
    cancelling.apply(ccx);
    let removed = SuperOpt::analyzer(3, 2)
        .with_murm(murm(3, 0))
        .run(&cancelling)
        .unwrap()
        .circuit;
    assert!(removed.gates.is_empty());
    assert!(!removed.has_toffoli());
}

#[test]
fn rz_windows_are_rejected_without_matrix_lookup() {
    let murm = murm(1, 0);
    for angles in [
        vec![0.37],
        vec![0.37, -0.37],
        vec![std::f64::consts::FRAC_PI_4; 2],
        vec![f64::NAN],
    ] {
        let mut circuit = Circuit::new(1);
        for angle in angles {
            circuit.apply(Gate::rz(angle, 0));
        }
        let result = SuperOpt::analyzer(1, circuit.gates.len())
            .with_murm(Arc::clone(&murm))
            .without_subcircuits()
            .run(&circuit)
            .unwrap();
        assert!(result.rewrites.is_empty());
        assert_eq!(result.cache_hits + result.cache_misses, 0);
        assert_eq!(result.circuit.gates.len(), circuit.gates.len());
    }
}

#[test]
fn coefficient_overflow_leaves_the_window_unchanged() {
    // This one-qubit Clifford+T word reaches a numerator coefficient outside
    // -127..=127. It is deliberately longer than the CLI presets so the test
    // exercises the public configurable-window overflow path.
    let word = "HTHTHTHTHTHTHTHTHTHTHTHTHTHTHTHTHTHTHTTTHTHTHTHTHTHTHTHTH";
    let mut circuit = Circuit::new(1);
    for gate in word.bytes() {
        circuit.apply(match gate {
            b'H' => Gate::h(0),
            b'T' => Gate::t(0),
            _ => unreachable!(),
        });
    }

    let result = SuperOpt::analyzer(1, circuit.gates.len())
        .run(&circuit)
        .unwrap();
    assert_eq!(result.circuit.num_qubits, circuit.num_qubits);
    assert_eq!(result.circuit.gates, circuit.gates);
    assert!(result.rewrites.is_empty());
    assert!(
        result
            .subcircuits
            .iter()
            .all(|window| window.gate_indices.len() < word.len()),
        "the overflowing full window should be omitted"
    );
}

#[test]
fn equal_length_synthesis_is_not_applied() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::h(0));
    let pass = SuperOpt::analyzer(1, 2).with_murm(murm(1, 2));
    let result = pass.run(&circuit).unwrap();
    assert!(result.rewrites.is_empty());
    assert_eq!(result.circuit.gates.len(), 2);
}

#[test]
fn identity_removal_wins_over_synthesis() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(0));
    let pass = SuperOpt::analyzer(1, 2).with_murm(murm(1, 2));
    let result = pass.run(&circuit).unwrap();
    assert!(result.circuit.gates.is_empty());
    assert_eq!(removed_indices(&result), vec![vec![0, 1]]);
    assert!(result.rewrites.iter().any(|r| r.replacement.is_empty()));
}

#[test]
fn consecutive_synth_rewrites_do_not_double_claim() {
    let mut circuit = Circuit::new(1);
    for _ in 0..4 {
        circuit.apply(Gate::s(0));
    }
    let pass = SuperOpt::analyzer(1, 2).with_murm(murm(1, 1));
    let result = pass.run(&circuit).unwrap();
    assert_eq!(result.rewrites.len(), 2);
    assert_eq!(result.circuit.gates.len(), 2);
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &result.circuit,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn a_second_superopt_pass_can_expose_a_new_rewrite() {
    let mut circuit = Circuit::new(1);
    for _ in 0..4 {
        circuit.apply(Gate::s(0));
    }
    let pass = SuperOpt::analyzer(1, 2).with_murm(murm(1, 1));
    let first = pass.run(&circuit).unwrap().circuit;
    assert_eq!(first.gates.len(), 2);
    assert!(first.gates.iter().all(|gate| matches!(gate, Gate::z(0))));

    let second = pass.run(&first).unwrap().circuit;
    assert!(second.gates.is_empty());
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &second,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn optimizes_unitary_regions_around_measure_and_reset() {
    let mut circuit = Circuit::with_cbits(1, 1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::measure { qubit: 0, cbit: 0 });
    circuit.apply(Gate::reset(0));
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(0));

    let result = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 2))
        .run(&circuit)
        .unwrap();
    assert_eq!(result.circuit.gates.len(), 2);
    assert!(matches!(result.circuit.gates[0], Gate::measure { .. }));
    assert!(matches!(result.circuit.gates[1], Gate::reset(0)));
    assert_eq!(removed_indices(&result), vec![vec![0, 1], vec![4, 5]]);
    assert!(
        result.subcircuits.iter().all(|window| {
            !window.gate_indices.contains(&2) && !window.gate_indices.contains(&3)
        })
    );
}

#[test]
fn rejects_out_of_range_qubit() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    let error = SuperOpt::analyzer(2, 1).run(&circuit).unwrap_err();
    assert_eq!(
        error,
        SuperOptError::InvalidQubit {
            gate_index: 0,
            qubit: 1,
            num_qubits: 1,
        }
    );
}

#[test]
fn pass_trait_optimizes_unitary_regions_in_mixed_circuit() {
    let mut circuit = Circuit::with_cbits(2, 1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::measure { qubit: 1, cbit: 0 });
    circuit.apply(Gate::h(0));
    let pass = SuperOpt::analyzer(1, 2).with_murm(murm(1, 2));
    let optimized = Pass::run(&pass, &circuit);

    assert_eq!(optimized.gates.len(), 1);
    assert!(matches!(optimized.gates[0], Gate::measure { .. }));
    assert_eq!(optimized.num_qubits, circuit.num_qubits);
    assert_eq!(optimized.num_cbits, circuit.num_cbits);
    assert!(optimized.has_measurement());
}

#[test]
fn measurement_blocks_a_later_bridge_into_its_qubit() {
    let mut circuit = Circuit::with_cbits(2, 1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::measure { qubit: 1, cbit: 0 });
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });

    let result = SuperOpt::analyzer(2, 3).run(&circuit).unwrap();
    assert!(
        result.subcircuits.iter().all(|window| {
            !(window.gate_indices.contains(&0) && window.gate_indices.contains(&2))
        })
    );
    assert!(matches!(result.circuit.gates[1], Gate::measure { .. }));
}

#[test]
fn measure_and_reset_each_block_same_qubit_windows() {
    let barriers = [Gate::measure { qubit: 0, cbit: 0 }, Gate::reset(0)];
    let murm = murm(1, 2);

    for barrier in barriers {
        let mut circuit = Circuit::with_cbits(1, 1);
        circuit.apply(Gate::h(0));
        circuit.apply(barrier);
        circuit.apply(Gate::h(0));

        let result = SuperOpt::analyzer(1, 3)
            .with_murm(Arc::clone(&murm))
            .run(&circuit)
            .unwrap();
        assert!(result.rewrites.is_empty());
        assert_eq!(result.circuit.gates.len(), 3);
        assert!(
            result
                .subcircuits
                .iter()
                .all(|window| !window.gate_indices.contains(&1))
        );
        assert!(
            result
                .subcircuits
                .iter()
                .all(|window| window.gate_indices != [0, 2])
        );
    }
}

#[test]
fn disjoint_measure_and_reset_do_not_block_unitary_window() {
    let mut circuit = Circuit::with_cbits(3, 1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::measure { qubit: 1, cbit: 0 });
    circuit.apply(Gate::reset(2));
    circuit.apply(Gate::h(0));

    let result = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 2))
        .run(&circuit)
        .unwrap();
    assert_eq!(removed_indices(&result), vec![vec![0, 3]]);
    assert_eq!(result.circuit.gates.len(), 2);
    assert!(matches!(result.circuit.gates[0], Gate::measure { .. }));
    assert!(matches!(result.circuit.gates[1], Gate::reset(2)));
    assert!(
        result.subcircuits.iter().all(|window| {
            !window.gate_indices.contains(&1) && !window.gate_indices.contains(&2)
        })
    );
}

#[test]
fn reset_blocks_a_later_bridge_into_its_qubit() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::reset(1));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });

    let result = SuperOpt::analyzer(2, 3).run(&circuit).unwrap();
    assert!(
        result.subcircuits.iter().all(|window| {
            !(window.gate_indices.contains(&0) && window.gate_indices.contains(&2))
        })
    );
    assert!(matches!(result.circuit.gates[1], Gate::reset(1)));
}

#[test]
fn retroactive_multihop_closure_stops_at_measurement_and_reset() {
    for barrier in [Gate::measure { qubit: 3, cbit: 0 }, Gate::reset(3)] {
        let mut circuit = Circuit::with_cbits(4, 1);
        circuit.apply(Gate::h(0));
        circuit.apply(barrier);
        circuit.apply(Gate::h(2));
        circuit.apply(Gate::cnot {
            control: 2,
            target: 3,
        });
        circuit.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        circuit.apply(Gate::x(0));

        let result = SuperOpt::analyzer(4, 6).run(&circuit).unwrap();
        assert!(result.subcircuits.iter().all(|window| {
            !(window.gate_indices.contains(&0) && window.gate_indices.contains(&4))
        }));
        assert!(
            result
                .subcircuits
                .iter()
                .any(|window| window.gate_indices == [4, 5])
        );
    }
}

#[test]
fn mixed_circuit_rewrites_preserve_the_full_quantum_classical_channel() {
    let mut circuit = Circuit::with_cbits(2, 1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    circuit.apply(Gate::h(1));
    circuit.apply(Gate::h(1));
    circuit.apply(Gate::measure { qubit: 0, cbit: 0 });
    circuit.apply(Gate::reset(0));
    circuit.apply(Gate::x(1));
    circuit.apply(Gate::x(1));

    let optimized = SuperOpt::analyzer(2, 2)
        .with_murm(murm(2, 2))
        .run(&circuit)
        .unwrap()
        .circuit;
    assert_eq!(optimized.gates.len(), 4);
    assert_channels_equivalent(&optimized, &circuit);

    // A noncontiguous identity may span a measurement on a disjoint qubit.
    let mut disjoint = Circuit::with_cbits(2, 1);
    disjoint.apply(Gate::h(0));
    disjoint.apply(Gate::measure { qubit: 1, cbit: 0 });
    disjoint.apply(Gate::h(0));
    let optimized = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 2))
        .run(&disjoint)
        .unwrap()
        .circuit;
    assert_eq!(optimized.gates.len(), 1);
    assert_channels_equivalent(&optimized, &disjoint);
}

#[test]
fn circuit_containing_only_measure_and_reset_is_unchanged() {
    let mut circuit = Circuit::with_cbits(2, 2);
    circuit.apply(Gate::measure { qubit: 0, cbit: 1 });
    circuit.apply(Gate::reset(1));
    circuit.apply(Gate::measure { qubit: 1, cbit: 0 });

    let result = SuperOpt::analyzer(2, 8).run(&circuit).unwrap();
    assert_eq!(result.circuit.to_qasm(), circuit.to_qasm());
    assert_eq!(result.circuit.num_cbits, 2);
    assert!(result.circuit.has_measurement());
    assert!(result.subcircuits.is_empty());
    assert!(result.rewrites.is_empty());
    assert_eq!(result.cache_hits + result.cache_misses, 0);
}

#[test]
fn invalid_measure_and_reset_qubits_are_still_rejected() {
    for gate in [Gate::measure { qubit: 2, cbit: 0 }, Gate::reset(2)] {
        let mut circuit = Circuit::with_cbits(2, 1);
        circuit.apply(gate);
        let error = SuperOpt::analyzer(2, 2).run(&circuit).unwrap_err();
        assert_eq!(
            error,
            SuperOptError::InvalidQubit {
                gate_index: 0,
                qubit: 2,
                num_qubits: 2,
            }
        );
    }
}

#[test]
fn mixed_circuit_without_subcircuits_still_rewrites_unitary_windows() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::reset(1));
    circuit.apply(Gate::h(0));

    let result = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 2))
        .without_subcircuits()
        .run(&circuit)
        .unwrap();
    assert!(result.subcircuits.is_empty());
    assert_eq!(removed_indices(&result), vec![vec![0, 2]]);
    assert_eq!(result.circuit.gates.len(), 1);
    assert!(matches!(result.circuit.gates[0], Gate::reset(1)));
}

#[test]
fn error_messages_name_the_offending_values() {
    let messages = [
        SuperOptError::InvalidQubit {
            gate_index: 3,
            qubit: 9,
            num_qubits: 4,
        }
        .to_string(),
        SuperOptError::MatrixTooLarge { num_qubits: 40 }.to_string(),
    ];
    assert!(messages[0].contains('9') && messages[0].contains('4'));
    assert!(messages[1].contains("40"));
}

#[test]
fn constructor_rejects_invalid_murm_config() {
    let error = SuperOpt::new(4, 8, MurmConfig::new(0, 8, 1_000)).unwrap_err();
    assert!(matches!(error, SuperOptError::InvalidMurmConfig { .. }));

    let error = SuperOpt::new(4, 8, MurmConfig::new(6, 8, 1_000)).unwrap_err();
    assert!(matches!(error, SuperOptError::InvalidMurmConfig { .. }));

    let error = SuperOpt::new(4, 8, MurmConfig::new(4, 8, 0)).unwrap_err();
    assert!(matches!(error, SuperOptError::InvalidMurmConfig { .. }));
}

#[test]
fn cache_stats_count_every_emission() {
    let mut rng = TestRng(0xc047_7000_0000_0001);
    let mut circuit = Circuit::new(4);
    for _ in 0..40 {
        let q = rng.qubit(4);
        let q2 = (q + 1 + rng.qubit(3)) % 4;
        circuit.apply(match rng.next(3) {
            0 => Gate::h(q),
            1 => Gate::t(q),
            _ => Gate::cnot {
                control: q,
                target: q2,
            },
        });
    }
    let result = SuperOpt::analyzer(3, 5).run(&circuit).unwrap();
    assert_eq!(
        result.cache_hits + result.cache_misses,
        result.subcircuits.len()
    );
}

#[test]
fn result_lists_are_sorted() {
    let mut circuit = Circuit::new(2);
    for _ in 0..3 {
        circuit.apply(Gate::h(0));
        circuit.apply(Gate::x(1));
        circuit.apply(Gate::h(0));
        circuit.apply(Gate::x(1));
    }
    let result = SuperOpt::analyzer(1, 2).run(&circuit).unwrap();
    assert!(
        result
            .subcircuits
            .windows(2)
            .all(|pair| pair[0].gate_indices <= pair[1].gate_indices)
    );
    assert!(
        removed_indices(&result)
            .windows(2)
            .all(|pair| pair[0] <= pair[1])
    );
}

#[test]
fn without_subcircuits_matches_collected_run() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::s(1));
    circuit.apply(Gate::s(1));

    let murm = murm(2, 2);
    let collected = SuperOpt::analyzer(2, 4)
        .with_murm(Arc::clone(&murm))
        .run(&circuit)
        .unwrap();
    let skipped = SuperOpt::analyzer(2, 4)
        .with_murm(murm)
        .without_subcircuits()
        .run(&circuit)
        .unwrap();

    assert!(skipped.subcircuits.is_empty());
    assert!(!collected.subcircuits.is_empty());
    // The optimization result must be identical regardless of collection.
    assert_eq!(
        format!("{:?}", collected.circuit.gates),
        format!("{:?}", skipped.circuit.gates)
    );
    assert_eq!(removed_indices(&collected), removed_indices(&skipped));
    assert_eq!(collected.rewrites.len(), skipped.rewrites.len());
    // Cache statistics may differ: without subcircuit collection the pass skips
    // the provably-unshortenable single-gate windows, so it performs strictly
    // fewer lookups here (h, cnot, h, s, s each anchor a skipped single-gate
    // window) while reaching the same rewrites.
    let collected_lookups = collected.cache_hits + collected.cache_misses;
    let skipped_lookups = skipped.cache_hits + skipped.cache_misses;
    assert!(skipped_lookups < collected_lookups);
}

#[test]
fn window_bound_counts_gates_not_index_span() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(1));
    circuit.apply(Gate::x(1));
    circuit.apply(Gate::h(0));
    let result = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 0))
        .run(&circuit)
        .unwrap();
    assert!(removed_indices(&result).contains(&vec![0, 3]));
    assert!(removed_indices(&result).contains(&vec![1, 2]));
    assert!(result.circuit.gates.is_empty());
}

#[test]
fn library_enumeration_excludes_rotation_gates() {
    for num_qubits in 1..=4 {
        for gate in library_gates(num_qubits, BASE_GATE_SET) {
            assert!(
                !matches!(gate.to_gate(), Gate::rz(..)),
                "enumeration must stay discrete: {gate:?}"
            );
        }
    }
}

#[test]
fn murm_keeps_a_smallest_circuit() {
    let murm = murm(1, 2);
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::s(0));
    circuit.apply(Gate::s(0));
    let matrix = naive_matrix(&circuit, &[0, 1], &[0]);

    let replacement = murm.synthesize(&matrix).unwrap();
    assert_eq!(replacement.len(), 1);
    assert!(matches!(replacement[0], Gate::z(0)));
    assert!(!murm.is_saturated(1));
}

#[test]
fn replaces_subcircuit_with_shorter_synthesized_circuit() {
    let mut circuit = Circuit::new(3);
    circuit.apply(Gate::h(2));
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::x(2));
    circuit.apply(Gate::h(2));

    let pass = SuperOpt::analyzer(1, 3).with_murm(murm(1, 1));
    let result = pass.run(&circuit).unwrap();

    assert_eq!(result.rewrites.len(), 1);
    assert_eq!(result.rewrites[0].gate_indices, vec![0, 2, 3]);
    assert_eq!(result.rewrites[0].replacement.len(), 1);
    assert!(matches!(result.rewrites[0].replacement[0], Gate::z(2)));
    assert_eq!(result.circuit.gates.len(), 2);
    assert!(matches!(result.circuit.gates[0], Gate::z(2)));
    assert!(matches!(result.circuit.gates[1], Gate::x(0)));
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &result.circuit,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn h_cnot_h_window_is_not_rewritten_to_cz() {
    // H·CX·H is the CZ unitary, but cz is excluded from the library, so the
    // shortest representative is the window itself — no rewrite may fire.
    let mut circuit = Circuit::new(4);
    circuit.apply(Gate::h(3));
    circuit.apply(Gate::cnot {
        control: 1,
        target: 3,
    });
    circuit.apply(Gate::h(3));

    let pass = SuperOpt::analyzer(2, 3).with_murm(murm(2, 3));
    let result = pass.run(&circuit).unwrap();

    assert!(
        result.rewrites.is_empty(),
        "no rewrite may introduce cz: {:?}",
        result.rewrites
    );
    assert_eq!(result.circuit.gates.len(), 3);
}

#[test]
fn synthesized_cnot_preserves_direction_on_sparse_physical_qubits() {
    for (control, target) in [(1, 3), (3, 1)] {
        let mut circuit = Circuit::new(5);
        circuit.apply(Gate::h(target));
        circuit.apply(Gate::cz { control, target });
        circuit.apply(Gate::h(target));

        let result = SuperOpt::analyzer(2, 3)
            .with_murm(murm(2, 1))
            .run(&circuit)
            .unwrap();
        assert_eq!(result.circuit.gates.len(), 1);
        assert!(matches!(
            result.circuit.gates[0],
            Gate::cnot {
                control: actual_control,
                target: actual_target,
            } if actual_control == control && actual_target == target
        ));
        assert!(crate::unitary::circuits_equiv(
            &circuit,
            &result.circuit,
            IDENTITY_TOLERANCE,
        ));
    }
}

#[test]
fn overlapping_synthesized_rewrites_are_not_both_applied() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::s(0));
    circuit.apply(Gate::s(0));
    circuit.apply(Gate::s(0));

    let pass = SuperOpt::analyzer(1, 2).with_murm(murm(1, 1));
    let result = pass.run(&circuit).unwrap();

    assert_eq!(result.rewrites.len(), 1);
    assert_eq!(result.rewrites[0].gate_indices, vec![0, 1]);
    assert!(matches!(result.rewrites[0].replacement[0], Gate::z(0)));
    assert_eq!(result.circuit.gates.len(), 2);
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &result.circuit,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn murm_reports_entry_cap() {
    let murm = Murm::build(MurmConfig {
        max_qubits: 4,
        max_gates: 6,
        max_entries_per_qubit: 100,
        basis: BASE_GATE_SET,
    })
    .unwrap();
    assert_eq!(murm.entry_count(4), 100);
    assert!(murm.is_saturated(4));
}

#[test]
fn murm_handles_identity_only_and_layer_boundary_caps() {
    let identity_only = Murm::build(MurmConfig::new(2, 3, 1)).unwrap();
    for width in 1..=2 {
        assert_eq!(identity_only.entry_count(width), 1);
        assert!(identity_only.is_saturated(width));
        assert_eq!(identity_only.completed_depth(width), 0);
    }

    // One-qubit depth one contains identity plus all seven library gates.
    let complete_first_layer = Murm::build(MurmConfig::new(1, 2, 8)).unwrap();
    assert_eq!(complete_first_layer.entry_count(1), 8);
    assert!(complete_first_layer.is_saturated(1));
    assert_eq!(complete_first_layer.completed_depth(1), 1);
}

#[test]
fn disk_round_trip_reproduces_the_built_table() {
    let config = MurmConfig::new(3, 6, 5_000);
    let built = Murm::build(config).unwrap();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("murm.bin");
    built.write_to_disk(&path, config).unwrap();
    let loaded = Murm::read_from_disk(&path, config).unwrap();

    for width in 1..=3 {
        assert_eq!(built.entry_count(width), loaded.entry_count(width));
        assert_eq!(built.is_saturated(width), loaded.is_saturated(width));
        assert_eq!(built.completed_depth(width), loaded.completed_depth(width));
    }

    // Every matrix the built MURM can synthesize, the loaded MURM
    // synthesizes identically (same replacement circuit, not just "a" match).
    let mut checked = 0;
    for num_qubits in 1..=3 {
        let support: Vec<Qubit> = (0..num_qubits as Qubit).collect();
        for gate in library_gates(num_qubits, BASE_GATE_SET) {
            let mut matrix = UnitaryMatrix::identity(num_qubits).unwrap();
            matrix.apply_gate_left(&gate.to_gate(), &support).unwrap();
            assert_eq!(built.synthesize(&matrix), loaded.synthesize(&matrix));
            checked += 1;
        }
    }
    assert!(checked > 0);
}

#[test]
fn disk_read_rejects_a_mismatched_config() {
    let config = MurmConfig::new(2, 4, 1_000);
    let built = Murm::build(config).unwrap();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("murm.bin");
    built.write_to_disk(&path, config).unwrap();

    let different_entries = MurmConfig::new(2, 4, 2_000);
    assert!(Murm::read_from_disk(&path, different_entries).is_err());

    let different_qubits = MurmConfig::new(1, 4, 1_000);
    assert!(Murm::read_from_disk(&path, different_qubits).is_err());

    let different_basis = MurmConfig::new(2, 4, 1_000)
        .with_basis(BASE_GATE_SET.union(GateSet::singleton(GateKind::Cz)));
    assert!(Murm::read_from_disk(&path, different_basis).is_err());
}

#[test]
fn disk_read_rejects_a_different_crate_version() {
    let config = MurmConfig::new(1, 2, 100);
    let built = Murm::build(config).unwrap();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("murm.bin");
    built.write_to_disk(&path, config).unwrap();

    // Corrupt the length-prefixed crate-version string right after the
    // magic (4 bytes) and format-version (4 bytes) fields, simulating a
    // cache file written by a different tzap release.
    let mut bytes = std::fs::read(&path).unwrap();
    let version_len = bytes[8] as usize;
    assert!(version_len > 0, "crate version string must be non-empty");
    let last_byte = &mut bytes[8 + version_len];
    *last_byte = last_byte.wrapping_add(1);
    std::fs::write(&path, &bytes).unwrap();

    assert!(Murm::read_from_disk(&path, config).is_err());
}

#[test]
fn disk_read_rejects_a_missing_or_corrupt_file() {
    let dir = tempfile::tempdir().unwrap();
    let missing = dir.path().join("does-not-exist.bin");
    let config = MurmConfig::new(1, 2, 100);
    assert!(Murm::read_from_disk(&missing, config).is_err());

    let corrupt = dir.path().join("corrupt.bin");
    std::fs::write(&corrupt, b"not a MURM cache").unwrap();
    assert!(Murm::read_from_disk(&corrupt, config).is_err());
}

fn murm_cache_body_offset(bytes: &[u8]) -> usize {
    // magic + format + version-length byte + version + qubits + gates +
    // entry cap + basis
    27 + usize::from(bytes[8])
}

#[test]
fn disk_read_rejects_structurally_corrupt_bodies_before_using_their_lengths() {
    let config = MurmConfig::new(2, 2, 100);
    let built = Murm::build(config).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("murm.bin");
    built.write_to_disk(&path, config).unwrap();
    let original = std::fs::read(&path).unwrap();
    let body = murm_cache_body_offset(&original);

    let reject = |bytes: &[u8], expected: &str| {
        std::fs::write(&path, bytes).unwrap();
        let error = Murm::read_from_disk(&path, config).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(
            error.to_string().contains(expected),
            "expected {expected:?}, got {error}"
        );
    };

    let mut huge_width_count = original.clone();
    huge_width_count[body..body + 4].copy_from_slice(&u32::MAX.to_le_bytes());
    reject(&huge_width_count, "width count");

    // The width-zero table length immediately follows the width count. Its
    // configured bound is zero, so this must fail before any allocation.
    let mut huge_table = original.clone();
    huge_table[body + 4..body + 12].copy_from_slice(&u64::MAX.to_le_bytes());
    reject(&huge_table, "configured or file-size bound");

    // width 0: 8-byte zero length. width 1: 8-byte length, then a 16-byte
    // identity root and its first child. A child may only point backward.
    let width_one_child = body + 4 + 8 + 8 + 16;
    let mut cyclic_parent = original.clone();
    cyclic_parent[width_one_child + 8..width_one_child + 12].copy_from_slice(&1u32.to_le_bytes());
    reject(&cyclic_parent, "parent must precede");

    let mut invalid_gate = original.clone();
    invalid_gate[width_one_child + 12..width_one_child + 16].copy_from_slice(&[7, 0, 0, 0]);
    reject(&invalid_gate, "outside its width or synthesis basis");
}

#[test]
fn disk_read_rejects_checksum_mismatches_and_trailing_bytes() {
    let config = MurmConfig::new(1, 2, 100);
    let built = Murm::build(config).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("murm.bin");
    built.write_to_disk(&path, config).unwrap();
    let original = std::fs::read(&path).unwrap();
    let body = murm_cache_body_offset(&original);

    // Change a non-root node's fingerprint while leaving the structure valid;
    // the body checksum must catch corruption that shape checks cannot.
    let width_one_child = body + 4 + 8 + 8 + 16;
    let mut changed_fingerprint = original.clone();
    changed_fingerprint[width_one_child] ^= 0x80;
    std::fs::write(&path, changed_fingerprint).unwrap();
    let error = Murm::read_from_disk(&path, config).unwrap_err();
    assert!(error.to_string().contains("checksum mismatch"), "{error}");

    let mut trailing = original;
    trailing.push(0);
    std::fs::write(&path, trailing).unwrap();
    let error = Murm::read_from_disk(&path, config).unwrap_err();
    assert!(error.to_string().contains("trailing bytes"), "{error}");
}

#[test]
fn parallel_murm_build_is_deterministic_across_thread_counts() {
    let config = MurmConfig::new(2, 3, 20_000);
    let build = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| Murm::build(config).unwrap())
    };
    let sequential = build(1);
    let parallel = build(4);
    for width in 1..=2 {
        assert_eq!(sequential.entry_count(width), parallel.entry_count(width));
        assert_eq!(sequential.is_saturated(width), parallel.is_saturated(width));
        assert_eq!(
            sequential.completed_depth(width),
            parallel.completed_depth(width)
        );
    }

    let gates = library_gates(2, BASE_GATE_SET);
    for &first in &gates {
        for &second in &gates {
            for &third in &gates {
                let matrix = library_circuit_matrix(2, &[first, second, third])
                    .unwrap()
                    .unwrap();
                let left = sequential.synthesize(&matrix).unwrap();
                let right = parallel.synthesize(&matrix).unwrap();
                assert_eq!(format!("{left:?}"), format!("{right:?}"));
            }
        }
    }
}

#[test]
fn shared_murm_cache_returns_one_arc_under_concurrency() {
    let config = MurmConfig::new(1, 2, 97);
    let barrier = Arc::new(std::sync::Barrier::new(8));
    let handles: Vec<_> = (0..8)
        .map(|_| {
            let barrier = Arc::clone(&barrier);
            std::thread::spawn(move || {
                barrier.wait();
                shared_murm(config).unwrap()
            })
        })
        .collect();
    let murms: Vec<_> = handles
        .into_iter()
        .map(|handle| handle.join().unwrap())
        .collect();
    assert!(murms[1..].iter().all(|murm| Arc::ptr_eq(&murms[0], murm)));
}

#[test]
fn randomized_synthesized_rewrites_preserve_unitary() {
    let murm = murm(3, 2);
    let mut rng = TestRng(0x51a7_4e51_5eed_cafe);
    for _ in 0..25 {
        let mut circuit = Circuit::new(3);
        for _ in 0..30 {
            let q = rng.qubit(3);
            let other = (q + 1 + rng.qubit(2)) % 3;
            let gate = match rng.next(8) {
                0 => Gate::h(q),
                1 => Gate::x(q),
                2 => Gate::s(q),
                3 => Gate::t(q),
                4 => Gate::tdg(q),
                5 | 6 => Gate::cnot {
                    control: q,
                    target: other,
                },
                _ => Gate::cz {
                    control: q,
                    target: other,
                },
            };
            circuit.apply(gate);
        }

        let pass = SuperOpt::analyzer(3, 4).with_murm(Arc::clone(&murm));
        let optimized = pass.run(&circuit).unwrap().circuit;
        assert!(crate::unitary::circuits_equiv(
            &circuit,
            &optimized,
            IDENTITY_TOLERANCE,
        ));
        assert!(optimized.gates.len() <= circuit.gates.len());
    }
}

/// Independently verify every rewrite the pass selected on `circuit`:
/// disjoint claims, commutation of skipped gates past the window span,
/// support-confined replacements, and matrix equality up to global phase.
/// Local soundness of each rewrite implies whole-circuit equivalence.
fn audit_rewrites(circuit: &Circuit, result: &SuperOptResult) {
    let mut claimed = vec![false; circuit.gates.len()];
    for rewrite in &result.rewrites {
        let mut support = Vec::new();
        for &gate_index in &rewrite.gate_indices {
            assert!(!claimed[gate_index], "gate {gate_index} claimed twice");
            claimed[gate_index] = true;
            support = union_qubits(&support, &unique_qubits(&circuit.gates[gate_index]));
        }

        let first = rewrite.gate_indices[0];
        let last = *rewrite.gate_indices.last().unwrap();
        for skipped in first..=last {
            if rewrite.gate_indices.binary_search(&skipped).is_ok() {
                continue;
            }
            for qubit in unique_qubits(&circuit.gates[skipped]) {
                assert!(
                    support.binary_search(&qubit).is_err(),
                    "skipped gate {skipped} touches window qubit {qubit}"
                );
            }
        }

        assert!(
            rewrite.replacement.len() < rewrite.gate_indices.len(),
            "rewrite must strictly shrink the circuit"
        );
        let mut replacement_matrix = UnitaryMatrix::identity(support.len()).unwrap();
        for gate in &rewrite.replacement {
            assert!(
                !matches!(gate, Gate::ccx { .. }),
                "rewrite introduced a Toffoli: {:?}",
                rewrite.gate_indices
            );
            for qubit in unique_qubits(gate) {
                assert!(
                    support.binary_search(&qubit).is_ok(),
                    "replacement leaves the window support"
                );
            }
            replacement_matrix.apply_gate_left(gate, &support).unwrap();
        }
        let original = naive_matrix(circuit, &rewrite.gate_indices, &support);
        assert!(
            original.equivalent_up_to_global_phase(&replacement_matrix),
            "replacement matrix differs for gates {:?}",
            rewrite.gate_indices
        );
    }
}

#[test]
fn randomized_production_config_rewrites_are_sound() {
    let murm = Arc::new(
        Murm::build(MurmConfig {
            max_qubits: 4,
            max_gates: 3,
            max_entries_per_qubit: 2_000,
            basis: BASE_GATE_SET,
        })
        .unwrap(),
    );
    let mut rng = TestRng(0xfab1_e5ca_1e50_44d5);
    for _ in 0..10 {
        let mut circuit = Circuit::new(5);
        for gate_index in 0..60 {
            let q = rng.qubit(5);
            let q2 = (q + 1 + rng.qubit(4)) % 5;
            let mut q3 = rng.qubit(5);
            while q3 == q || q3 == q2 {
                q3 = rng.qubit(5);
            }
            let gate = match rng.next(12) {
                0 => Gate::x(q),
                1 => Gate::h(q),
                2 => Gate::s(q),
                3 => Gate::sdg(q),
                4 => Gate::z(q),
                5 => Gate::t(q),
                6 => Gate::tdg(q),
                7 => Gate::rz((gate_index + 1) as f64 / 13.0, q),
                8 | 9 => Gate::cnot {
                    control: q,
                    target: q2,
                },
                10 => Gate::cz {
                    control: q,
                    target: q2,
                },
                _ => Gate::ccx {
                    control1: q,
                    control2: q2,
                    target: q3,
                },
            };
            circuit.apply(gate);
        }

        let pass = SuperOpt::analyzer(4, 8).with_murm(Arc::clone(&murm));
        let result = pass.run(&circuit).unwrap();
        audit_rewrites(&circuit, &result);
        assert!(result.circuit.gates.len() <= circuit.gates.len());
        assert!(crate::unitary::circuits_equiv(
            &circuit,
            &result.circuit,
            IDENTITY_TOLERANCE,
        ));
    }
}

#[test]
fn removes_noncontiguous_identity_subcircuit_from_circuit() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(1));
    circuit.apply(Gate::h(0));

    let result = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 0))
        .run(&circuit)
        .unwrap();
    assert_eq!(removed_indices(&result), vec![vec![0, 2]]);
    assert_eq!(result.circuit.gates.len(), 1);
    assert!(matches!(result.circuit.gates[0], Gate::x(1)));
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &result.circuit,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn checks_identity_windows_shorter_than_gate_limit() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(0));

    let result = SuperOpt::analyzer(1, 8)
        .with_murm(murm(1, 0))
        .run(&circuit)
        .unwrap();
    assert_eq!(removed_indices(&result), vec![vec![0, 1]]);
    assert!(result.circuit.gates.is_empty());
    assert!(
        result
            .subcircuits
            .iter()
            .any(|window| window.gate_indices == [0])
    );
    assert!(
        result
            .subcircuits
            .iter()
            .any(|window| window.gate_indices == [0, 1])
    );
}

#[test]
fn removes_identity_up_to_global_phase() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::z(0));
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::z(0));

    let result = SuperOpt::analyzer(1, 4)
        .with_murm(murm(1, 0))
        .run(&circuit)
        .unwrap();
    assert_eq!(removed_indices(&result), vec![vec![0, 1, 2, 3]]);
    assert!(result.circuit.gates.is_empty());
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &result.circuit,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn overlapping_identity_windows_are_not_both_removed() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::x(0));

    let result = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 0))
        .run(&circuit)
        .unwrap();
    assert_eq!(removed_indices(&result), vec![vec![0, 1]]);
    assert_eq!(result.circuit.gates.len(), 1);
    assert!(matches!(result.circuit.gates[0], Gate::x(0)));
    assert!(crate::unitary::circuits_equiv(
        &circuit,
        &result.circuit,
        IDENTITY_TOLERANCE,
    ));
}

#[test]
fn nonidentity_window_is_preserved() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(0));

    let result = SuperOpt::analyzer(1, 2).run(&circuit).unwrap();
    assert!(removed_indices(&result).is_empty());
    assert_eq!(result.circuit.gates.len(), 2);
}

#[test]
fn implements_optimization_pass_interface() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(0));

    let pass = SuperOpt::analyzer(1, 2).with_murm(murm(1, 0));
    let optimized = Pass::run(&pass, &circuit);
    assert!(optimized.gates.is_empty());
}

#[test]
fn every_supported_unitary_gate_matches_naive_matrix() {
    let mut circuit = Circuit::new(3);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::x(0));
    circuit.apply(Gate::s(0));
    circuit.apply(Gate::sdg(0));
    circuit.apply(Gate::z(0));
    circuit.apply(Gate::t(0));
    circuit.apply(Gate::tdg(0));
    circuit.apply(Gate::rz(0.29, 0));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    circuit.apply(Gate::cz {
        control: 1,
        target: 0,
    });
    circuit.apply(Gate::ccx {
        control1: 0,
        control2: 1,
        target: 2,
    });

    let result = SuperOpt::analyzer(3, 3).run(&circuit).unwrap();
    for window in &result.subcircuits {
        let expected = naive_matrix(&circuit, &window.gate_indices, &window.qubits);
        assert_matrix_equal(&window.matrix, &expected);
    }
}

#[test]
fn randomized_results_match_naive_anchored_scan() {
    let mut rng = TestRng(0x5eed_1234_9876_abcd);
    for _ in 0..30 {
        let mut circuit = Circuit::new(4);
        for _ in 0..30 {
            let q = rng.qubit(4);
            let q2 = (q + 1 + rng.qubit(3)) % 4;
            let mut q3 = rng.qubit(4);
            while q3 == q || q3 == q2 {
                q3 = rng.qubit(4);
            }
            let gate = match rng.next(11) {
                0 => Gate::x(q),
                1 => Gate::h(q),
                2 => Gate::s(q),
                3 => Gate::sdg(q),
                4 => Gate::z(q),
                5 => Gate::t(q),
                6 => Gate::tdg(q),
                7 => Gate::rz(rng.next(100) as f64 / 17.0, q),
                8 => Gate::cnot {
                    control: q,
                    target: q2,
                },
                9 => Gate::cz {
                    control: q,
                    target: q2,
                },
                _ => Gate::ccx {
                    control1: q,
                    control2: q2,
                    target: q3,
                },
            };
            circuit.apply(gate);
        }

        for max_qubits in 1..=3 {
            for window_gates in 1..=5 {
                let result = SuperOpt::analyzer(max_qubits, window_gates)
                    .run(&circuit)
                    .unwrap();
                let expected = naive_windows(&circuit, max_qubits, window_gates);
                let actual: Vec<_> = result
                    .subcircuits
                    .iter()
                    .map(|window| (window.gate_indices.clone(), window.qubits.clone()))
                    .collect();
                assert_eq!(actual, expected);
                for window in &result.subcircuits {
                    let expected = naive_matrix(&circuit, &window.gate_indices, &window.qubits);
                    assert_matrix_equal(&window.matrix, &expected);
                }
            }
        }
    }
}

#[test]
fn rejects_zero_gate_window() {
    let error = SuperOpt::analyzer(2, 0).run(&Circuit::new(2)).unwrap_err();
    assert_eq!(error, SuperOptError::ZeroWindowGates);
}

fn load_optimizer_profile_murm() -> Arc<Murm> {
    use std::time::Instant;

    let start = Instant::now();
    let murm = shared_murm(MurmConfig::default()).unwrap();
    println!(
        "initialized optimizer MURM in {:.3} s: entries {:?}, complete depths {:?}",
        start.elapsed().as_secs_f64(),
        (1..=4)
            .map(|num_qubits| murm.entry_count(num_qubits))
            .collect::<Vec<_>>(),
        (1..=4)
            .map(|num_qubits| murm.completed_depth(num_qubits))
            .collect::<Vec<_>>(),
    );
    murm
}

#[test]
#[ignore = "manual release-mode equivalence check on small benchmarks"]
fn verify_small_benchmarks_preserve_unitary() {
    let murm = load_optimizer_profile_murm();
    for name in ["tof_3", "barenco_tof_3", "mod5_4", "hwb6"] {
        let path = format!(
            "{}/benchmarks/feynman/{name}.qasm",
            env!("CARGO_MANIFEST_DIR")
        );
        let circuit = crate::qasm::parse(&std::fs::read_to_string(&path).unwrap()).unwrap();
        let pass = SuperOpt::analyzer(4, 8)
            .with_murm(Arc::clone(&murm))
            .without_subcircuits();
        let result = pass.run(&circuit).unwrap();
        assert!(
            crate::unitary::circuits_equiv(&circuit, &result.circuit, IDENTITY_TOLERANCE),
            "{name} rewrite changed the unitary"
        );

        let default_result = Pass::run(
            &crate::phase_fold_rand::PhaseFoldRand,
            &Pass::run(&crate::cancel::CancelGates, &circuit),
        );
        let pipelined = pass.run(&default_result).unwrap();
        assert!(
            crate::unitary::circuits_equiv(&circuit, &pipelined.circuit, IDENTITY_TOLERANCE),
            "{name} default+subcircuit pipeline changed the unitary"
        );
        println!(
            "{name}: {} gates -> {} (standalone) / {} (default+subcircuit), unitary preserved",
            circuit.gates.len(),
            result.circuit.gates.len(),
            pipelined.circuit.gates.len(),
        );
    }
}

#[test]
#[ignore = "manual release-mode randomized fuzz with guaranteed rewrites"]
fn fuzz_subcircuit_rewrites_change_circuit_and_preserve_unitary() {
    let murm = load_optimizer_profile_murm();
    let mut rng = TestRng(0xdeed_beef_5eed_0001);
    let num_cases = 10_000;
    for round in 0..num_cases {
        let mut circuit = Circuit::new(6);

        // Guarantee that every fuzz case exercises the rewrite path. The
        // gadget varies across single-, two-, and three-qubit gates.
        match rng.next(8) {
            0 => {
                let q = rng.qubit(6);
                circuit.apply(Gate::h(q));
                circuit.apply(Gate::h(q));
            }
            1 => {
                let q = rng.qubit(6);
                circuit.apply(Gate::x(q));
                circuit.apply(Gate::x(q));
            }
            2 => {
                let q = rng.qubit(6);
                circuit.apply(Gate::z(q));
                circuit.apply(Gate::z(q));
            }
            3 => {
                let q = rng.qubit(6);
                circuit.apply(Gate::s(q));
                circuit.apply(Gate::sdg(q));
            }
            4 => {
                let q = rng.qubit(6);
                circuit.apply(Gate::t(q));
                circuit.apply(Gate::tdg(q));
            }
            5 => {
                circuit.apply(Gate::cnot {
                    control: 0,
                    target: 1,
                });
                circuit.apply(Gate::cnot {
                    control: 0,
                    target: 1,
                });
            }
            6 => {
                circuit.apply(Gate::cz {
                    control: 0,
                    target: 1,
                });
                circuit.apply(Gate::cz {
                    control: 1,
                    target: 0,
                });
            }
            _ => {
                circuit.apply(Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                });
                circuit.apply(Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                });
            }
        }

        for gate_index in 0..100 {
            let q = rng.qubit(6);
            let q2 = (q + 1 + rng.qubit(5)) % 6;
            let mut q3 = rng.qubit(6);
            while q3 == q || q3 == q2 {
                q3 = rng.qubit(6);
            }
            let gate = match rng.next(12) {
                0 => Gate::x(q),
                1 => Gate::h(q),
                2 => Gate::s(q),
                3 => Gate::sdg(q),
                4 => Gate::z(q),
                5 => Gate::t(q),
                6 => Gate::tdg(q),
                7 => Gate::rz((gate_index + 1) as f64 / 13.0, q),
                8 | 9 => Gate::cnot {
                    control: q,
                    target: q2,
                },
                10 => Gate::cz {
                    control: q,
                    target: q2,
                },
                _ => Gate::ccx {
                    control1: q,
                    control2: q2,
                    target: q3,
                },
            };
            circuit.apply(gate);
        }

        let pass = SuperOpt::analyzer(4, 8).with_murm(Arc::clone(&murm));
        let result = pass.run(&circuit).unwrap();
        assert!(!result.rewrites.is_empty(), "round {round} made no rewrite");
        assert!(
            result.circuit.gates.len() < circuit.gates.len(),
            "round {round} did not change the circuit"
        );
        audit_rewrites(&circuit, &result);
        assert!(
            crate::unitary::circuits_equiv(&circuit, &result.circuit, IDENTITY_TOLERANCE),
            "round {round} changed the unitary"
        );
    }
    println!(
        "{num_cases} random 6-qubit circuits: all changed, rewrites audited, unitaries preserved"
    );
}

#[test]
#[ignore = "manual release-mode rewrite audit over the full benchmark corpora"]
fn audit_all_benchmark_rewrites() {
    use crate::cancel::CancelGates;
    use crate::phase_fold_rand::PhaseFoldRand;

    let mut paths = Vec::new();
    for corpus in ["feynman", "cobble-t"] {
        let directory = format!("{}/benchmarks/{corpus}", env!("CARGO_MANIFEST_DIR"));
        paths.extend(
            std::fs::read_dir(directory)
                .unwrap()
                .map(|entry| entry.unwrap().path())
                .filter(|path| {
                    path.extension()
                        .is_some_and(|extension| extension == "qasm")
                }),
        );
    }
    paths.sort();

    let murm = load_optimizer_profile_murm();
    let pass = SuperOpt::analyzer(4, 8).with_murm(Arc::clone(&murm));

    let mut total_rewrites = 0;
    for path in paths {
        let name = path.file_stem().unwrap().to_string_lossy();
        let circuit = crate::qasm::parse(&std::fs::read_to_string(&path).unwrap()).unwrap();

        let standalone = pass.run(&circuit).unwrap();
        audit_rewrites(&circuit, &standalone);

        let default_result = Pass::run(&PhaseFoldRand, &Pass::run(&CancelGates, &circuit));
        let pipelined = pass.run(&default_result).unwrap();
        audit_rewrites(&default_result, &pipelined);

        total_rewrites += standalone.rewrites.len() + pipelined.rewrites.len();
        println!(
            "{name}: {} standalone + {} pipelined rewrites audited",
            standalone.rewrites.len(),
            pipelined.rewrites.len(),
        );
    }
    println!("TOTAL: {total_rewrites} rewrites audited");
}

#[test]
fn warm_store_reproduces_cold_run() {
    let mut circuit = Circuit::new(2);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    circuit.apply(Gate::t(1));

    let pass = SuperOpt::analyzer(2, 3).with_murm(murm(2, 3));
    let cold = pass.run(&circuit).unwrap();
    let warm = pass.run(&circuit).unwrap();

    assert!(cold.cache_misses > 0);
    // The warm run performs the same lookups, now all hits.
    assert_eq!(warm.cache_misses, 0);
    assert_eq!(warm.cache_hits, cold.cache_hits + cold.cache_misses);
    assert_eq!(removed_indices(&cold), removed_indices(&warm));
    assert_eq!(
        format!("{:?}", cold.circuit.gates),
        format!("{:?}", warm.circuit.gates)
    );
    assert_eq!(cold.subcircuits.len(), warm.subcircuits.len());
}

#[test]
fn incremental_skips_unchanged_circuit_entirely() {
    let mut circuit = Circuit::new(1);
    circuit.apply(Gate::h(0));
    circuit.apply(Gate::rz(0.3, 0));
    circuit.apply(Gate::h(0));

    let pass = SuperOpt::analyzer(1, 3)
        .with_murm(murm(1, 3))
        .without_subcircuits()
        .incremental();
    let first = pass.run(&circuit).unwrap();
    assert!(first.rewrites.is_empty());

    // An unchanged input has an empty frontier: no windows, no lookups.
    let second = pass.run(&circuit).unwrap();
    assert!(second.rewrites.is_empty());
    assert_eq!(second.cache_hits + second.cache_misses, 0);
}

#[test]
fn incremental_finds_rewrites_exposed_by_deletion() {
    // A deletion leaves no dirty gate of its own; the flanking survivors it
    // makes adjacent must be re-anchored or the new HH cancellation is lost.
    let mut before = Circuit::new(1);
    before.apply(Gate::h(0));
    before.apply(Gate::rz(0.3, 0));
    before.apply(Gate::h(0));
    let mut after = Circuit::new(1);
    after.apply(Gate::h(0));
    after.apply(Gate::h(0));

    let pass = SuperOpt::analyzer(1, 2)
        .with_murm(murm(1, 2))
        .without_subcircuits()
        .incremental();
    assert!(pass.run(&before).unwrap().rewrites.is_empty());
    let result = pass.run(&after).unwrap();
    assert!(result.circuit.gates.is_empty());
}

#[test]
fn incremental_matches_full_sweeps_on_random_circuits() {
    use crate::cancel::CancelGates;

    let mut rng = TestRng(0x1acf_1e90_b5e5_5ed1);
    let murm = murm(2, 4);
    let incremental = SuperOpt::analyzer(2, 4)
        .with_murm(Arc::clone(&murm))
        .without_subcircuits()
        .incremental();
    let full = SuperOpt::analyzer(2, 4)
        .with_murm(Arc::clone(&murm))
        .without_subcircuits();

    for _ in 0..20 {
        let mut circuit = Circuit::new(4);
        for _ in 0..50 {
            let q = rng.qubit(4);
            let q2 = (q + 1 + rng.qubit(3)) % 4;
            let gate = match rng.next(10) {
                0 => Gate::x(q),
                1 => Gate::h(q),
                2 => Gate::s(q),
                3 => Gate::sdg(q),
                4 => Gate::z(q),
                5 => Gate::t(q),
                6 => Gate::tdg(q),
                7 => Gate::rz(rng.next(100) as f64 / 17.0, q),
                8 => Gate::cnot {
                    control: q,
                    target: q2,
                },
                _ => Gate::cz {
                    control: q,
                    target: q2,
                },
            };
            circuit.apply(gate);
        }

        // Emulate the sequential fixpoint driver: the incremental instance
        // sees each successive circuit version, with CancelGates mutating it
        // between SuperOpt runs, and must match a full sweep every round.
        let mut current = circuit;
        for _ in 0..8 {
            let expected = full.run(&current).unwrap().circuit;
            let actual = incremental.run(&current).unwrap().circuit;
            assert_eq!(
                format!("{:?}", actual.gates),
                format!("{:?}", expected.gates)
            );
            let next = Pass::run(&CancelGates, &actual);
            if format!("{:?}", next.gates) == format!("{:?}", current.gates) {
                break;
            }
            current = next;
        }
    }
}
