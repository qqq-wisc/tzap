use super::*;
use crate::unitary::{C, circuit_unitary, circuits_equiv};
use std::collections::BTreeMap;

fn next(rng: &mut u64) -> u64 {
    *rng ^= *rng << 13;
    *rng ^= *rng >> 7;
    *rng ^= *rng << 17;
    *rng
}

const KINDS: [RotationKind; 4] = [
    RotationKind::X,
    RotationKind::Y,
    RotationKind::Z,
    RotationKind::Phase,
];

fn circuit(n: usize, gates: Vec<Gate>) -> Circuit {
    Circuit::with_cbits(n, n).replacing_gates(gates)
}

fn verify(c: &Circuit, expected: &[Gate]) {
    let out = phase_fold_pauli(c);
    assert_eq!(out.gates, expected, "input: {:?}", c.gates);
    assert!(circuits_equiv(c, &out, 1e-12), "{c}\n{out}");
}

#[test]
fn all_native_axes_merge_and_cancel_without_changing_spelling() {
    for kind in KINDS {
        verify(
            &circuit(
                1,
                vec![kind.numeric_gate(0.125, 0), kind.numeric_gate(0.25, 0)],
            ),
            &[kind.numeric_gate(0.375, 0)],
        );
        verify(
            &circuit(
                1,
                vec![kind.numeric_gate(0.125, 0), kind.numeric_gate(-0.125, 0)],
            ),
            &[],
        );
        for sign in [1.0, -1.0] {
            verify(
                &circuit(
                    1,
                    vec![
                        kind.numeric_gate(sign * 0.25, 0),
                        kind.numeric_gate(sign * 0.25, 0),
                    ],
                ),
                &[kind.numeric_gate(sign * 0.5, 0)],
            );
        }
    }
    for (first, last) in [
        (RotationKind::Z, RotationKind::Phase),
        (RotationKind::Phase, RotationKind::Z),
    ] {
        verify(
            &circuit(
                1,
                vec![first.numeric_gate(0.125, 0), last.numeric_gate(0.25, 0)],
            ),
            &[last.numeric_gate(0.375, 0)],
        );
    }
}

#[test]
fn cross_axis_merges_follow_h_s_and_sx_with_correct_signs() {
    for (first, clifford, last, sign) in [
        (RotationKind::X, Gate::h(0), RotationKind::Z, 1.0),
        (RotationKind::Z, Gate::h(0), RotationKind::X, 1.0),
        (RotationKind::Y, Gate::h(0), RotationKind::Y, -1.0),
        (RotationKind::X, Gate::s(0), RotationKind::Y, 1.0),
        (RotationKind::Y, Gate::s(0), RotationKind::X, -1.0),
        (RotationKind::X, Gate::sdg(0), RotationKind::Y, -1.0),
        (RotationKind::Y, Gate::sdg(0), RotationKind::X, 1.0),
        (RotationKind::Y, Gate::sx(0), RotationKind::Z, 1.0),
        (RotationKind::Z, Gate::sx(0), RotationKind::Y, -1.0),
    ] {
        verify(
            &circuit(
                1,
                vec![
                    first.numeric_gate(0.125, 0),
                    clifford.clone(),
                    last.numeric_gate(0.25, 0),
                ],
            ),
            &[clifford, last.numeric_gate(0.25 + sign * 0.125, 0)],
        );
    }
}

type MixedChain = (
    RotationKind,
    Gate,
    RotationKind,
    Gate,
    RotationKind,
    [f64; 3],
);

/// Each chain visits X, Y and Z once. Coefficients express the three angles
/// in the final gate's physical axis, derived from H/S/SX conjugation rules.
fn mixed_chains() -> [MixedChain; 6] {
    use RotationKind::{X, Y, Z};
    [
        (X, Gate::s(0), Y, Gate::sx(0), Z, [1.0, 1.0, 1.0]),
        (X, Gate::h(0), Z, Gate::sx(0), Y, [-1.0, -1.0, 1.0]),
        (Y, Gate::s(0), X, Gate::h(0), Z, [-1.0, 1.0, 1.0]),
        (Y, Gate::sx(0), Z, Gate::h(0), X, [1.0, 1.0, 1.0]),
        (Z, Gate::h(0), X, Gate::s(0), Y, [1.0, 1.0, 1.0]),
        (Z, Gate::sx(0), Y, Gate::s(0), X, [1.0, -1.0, 1.0]),
    ]
}

#[test]
fn every_rx_ry_rz_order_folds_to_the_last_axis_with_signed_angles() {
    for (first, clifford1, middle, clifford2, last, coefficients) in mixed_chains() {
        for signs in 0..8 {
            let angles: [f64; 3] = std::array::from_fn(|i| {
                let magnitude = [0.125, 0.25, 0.5][i];
                if signs & (1 << i) == 0 {
                    magnitude
                } else {
                    -magnitude
                }
            });
            let expected_angle: f64 = angles.iter().zip(coefficients).map(|(a, c)| a * c).sum();
            let c = circuit(
                1,
                vec![
                    first.numeric_gate(angles[0], 0),
                    clifford1.clone(),
                    middle.numeric_gate(angles[1], 0),
                    clifford2.clone(),
                    last.numeric_gate(angles[2], 0),
                ],
            );
            let expected = vec![
                clifford1.clone(),
                clifford2.clone(),
                last.numeric_gate(expected_angle, 0),
            ];
            verify(&c, &expected);
            // Both history implementations must reach the same conclusion,
            // even when every sketch fingerprint collides.
            for sliced in [false, true] {
                assert_eq!(
                    fold(&c, &[(0, 0)], MAX_ATTEMPT_STEPS, sliced).gates,
                    expected
                );
            }
            let once = phase_fold_pauli(&c);
            assert_eq!(phase_fold_pauli(&once).gates, expected);
        }
    }
}

#[test]
fn all_rx_ry_rz_orders_block_folds_across_noncommuting_rotations() {
    for (first, _, middle, _, last, _) in mixed_chains() {
        let c = circuit(
            1,
            vec![
                first.numeric_gate(0.125, 0),
                middle.numeric_gate(0.25, 0),
                last.numeric_gate(0.5, 0),
                first.numeric_gate(0.125, 0),
            ],
        );
        verify(&c, &c.gates);
    }
}

#[test]
fn rx_ry_rz_on_distinct_wires_fold_with_overlapping_conjugated_supports() {
    for (first, _, middle, _, last, _) in mixed_chains() {
        let cx = Gate::cnot {
            control: 0,
            target: 1,
        };
        let cy = Gate::cy {
            control: 1,
            target: 2,
        };
        let c = circuit(
            3,
            vec![
                cx.clone(),
                cy.clone(),
                first.numeric_gate(0.125, 0),
                middle.numeric_gate(0.25, 1),
                last.numeric_gate(0.5, 2),
                first.numeric_gate(0.25, 0),
                middle.numeric_gate(0.5, 1),
                last.numeric_gate(1.0, 2),
            ],
        );
        verify(
            &c,
            &[
                cx,
                cy,
                first.numeric_gate(0.375, 0),
                middle.numeric_gate(0.75, 1),
                last.numeric_gate(1.5, 2),
            ],
        );
    }
}

#[test]
fn pauli_sign_reversals_and_swapped_wires_are_retained() {
    for kind in KINDS {
        for (pauli, sign) in [
            (Gate::x(0), if kind == RotationKind::X { 1.0 } else { -1.0 }),
            (Gate::y(0), if kind == RotationKind::Y { 1.0 } else { -1.0 }),
            (
                Gate::z(0),
                if matches!(kind, RotationKind::Z | RotationKind::Phase) {
                    1.0
                } else {
                    -1.0
                },
            ),
        ] {
            verify(
                &circuit(
                    1,
                    vec![
                        kind.numeric_gate(0.125, 0),
                        pauli.clone(),
                        kind.numeric_gate(0.25, 0),
                    ],
                ),
                &[pauli, kind.numeric_gate(0.25 + sign * 0.125, 0)],
            );
        }
        verify(
            &circuit(
                2,
                vec![
                    kind.numeric_gate(0.125, 0),
                    Gate::swap(0, 1),
                    kind.numeric_gate(0.25, 1),
                ],
            ),
            &[Gate::swap(0, 1), kind.numeric_gate(0.375, 1)],
        );
        verify(
            &circuit(
                2,
                vec![
                    kind.numeric_gate(0.125, 1),
                    Gate::swap(0, 1),
                    kind.numeric_gate(0.25, 0),
                ],
            ),
            &[Gate::swap(0, 1), kind.numeric_gate(0.375, 0)],
        );
    }
}

#[test]
fn entangling_clifford_axes_fold_only_when_equal() {
    for clifford in [
        Gate::cnot {
            control: 0,
            target: 1,
        },
        Gate::cz {
            control: 0,
            target: 1,
        },
        Gate::cy {
            control: 0,
            target: 1,
        },
    ] {
        for kind in KINDS {
            for q in 0..2 {
                let c = circuit(
                    2,
                    vec![
                        kind.numeric_gate(0.125, q),
                        clifford.clone(),
                        kind.numeric_gate(0.25, q),
                    ],
                );
                let commutes = match (&clifford, kind, q) {
                    (Gate::cnot { .. }, RotationKind::X, 1) => true,
                    (Gate::cy { .. }, RotationKind::Y, 1) => true,
                    (Gate::cz { .. }, RotationKind::Z | RotationKind::Phase, _) => true,
                    (_, RotationKind::Z | RotationKind::Phase, 0) => true,
                    _ => false,
                };
                let expected = if commutes {
                    vec![clifford.clone(), kind.numeric_gate(0.375, q)]
                } else {
                    c.gates.clone()
                };
                verify(&c, &expected);
            }
        }
        for kind in KINDS {
            verify(
                &circuit(
                    2,
                    vec![
                        kind.numeric_gate(0.125, 0),
                        clifford.clone(),
                        Gate::h(1),
                        Gate::h(1),
                        clifford.clone(),
                        kind.numeric_gate(0.25, 0),
                    ],
                ),
                &[
                    clifford.clone(),
                    Gate::h(1),
                    Gate::h(1),
                    clifford.clone(),
                    kind.numeric_gate(0.375, 0),
                ],
            );
        }
    }
}

#[test]
fn anticommuting_rotations_and_opaque_gates_block_merges() {
    for first in KINDS {
        for middle in KINDS {
            let commute = first == middle
                || matches!(
                    (first, middle),
                    (RotationKind::Z, RotationKind::Phase) | (RotationKind::Phase, RotationKind::Z)
                );
            let c = circuit(
                1,
                vec![
                    first.numeric_gate(0.125, 0),
                    middle.numeric_gate(0.25, 0),
                    first.numeric_gate(0.125, 0),
                ],
            );
            let expected = if commute {
                vec![first.numeric_gate(0.5, 0)]
            } else {
                c.gates.clone()
            };
            verify(&c, &expected);
        }
        for barrier in [
            Gate::cp {
                lambda: crate::angle::Angle::from_f64(0.125).unwrap(),
                control: 0,
                target: 1,
            },
            Gate::crx {
                theta: crate::angle::Angle::from_f64(0.125).unwrap(),
                control: 0,
                target: 1,
            },
            Gate::cry {
                theta: crate::angle::Angle::from_f64(0.125).unwrap(),
                control: 0,
                target: 1,
            },
            Gate::crz {
                theta: crate::angle::Angle::from_f64(0.125).unwrap(),
                control: 0,
                target: 1,
            },
            Gate::ch {
                control: 0,
                target: 1,
            },
            Gate::cswap {
                control: 0,
                first: 1,
                second: 2,
            },
        ] {
            let c = circuit(
                3,
                vec![
                    Gate::s(0),
                    first.numeric_gate(0.125, 0),
                    barrier.clone(),
                    first.numeric_gate(0.25, 0),
                    first.numeric_gate(0.125, 0),
                ],
            );
            verify(
                &c,
                &[
                    Gate::s(0),
                    first.numeric_gate(0.125, 0),
                    barrier,
                    first.numeric_gate(0.375, 0),
                ],
            );
        }
    }
}

#[test]
fn forced_collisions_and_word_boundaries_do_not_change_axis_decisions() {
    for n in [1, 2, 32, 33, 64, 65, 129] {
        let q = (n - 1) as u32;
        for kind in KINDS {
            let c = circuit(
                n,
                vec![
                    kind.numeric_gate(0.125, q),
                    Gate::h(q),
                    Gate::h(q),
                    kind.numeric_gate(0.25, q),
                ],
            );
            for labels in [
                vec![(0, 0); n],
                (0..n)
                    .map(|i| (2 * i as u128 + 1, 2 * i as u128 + 2))
                    .collect(),
            ] {
                for sliced in [false, true]
                    .into_iter()
                    .filter(|&s| !s || n <= MAX_SLICED_QUBITS)
                {
                    let out = fold(&c, &labels, MAX_ATTEMPT_STEPS, sliced);
                    assert_eq!(
                        out.gates,
                        vec![Gate::h(q), Gate::h(q), kind.numeric_gate(0.375, q)]
                    );
                }
            }
        }
        let c = circuit(
            n,
            vec![
                Gate::rx_f64(0.125, q).unwrap(),
                Gate::ry_f64(0.25, q).unwrap(),
                Gate::rx_f64(0.125, q).unwrap(),
            ],
        );
        assert_eq!(
            fold(&c, &vec![(0, 0); n], MAX_ATTEMPT_STEPS, false).gates,
            c.gates
        );
    }
}

#[test]
fn huge_tiny_and_near_clifford_angles_are_never_reduced_or_snapped() {
    for kind in KINDS {
        for a in [
            1e16,
            -1e16,
            1e30,
            -1e30,
            1e308,
            -1e308,
            f64::from_bits(1),
            f64::MIN_POSITIVE,
            5e-7,
            PI / 4.0,
            PI / 4.0 + 5e-7,
        ] {
            let c = circuit(
                1,
                vec![
                    kind.numeric_gate(a, 0),
                    Gate::h(0),
                    Gate::h(0),
                    kind.numeric_gate(-a, 0),
                ],
            );
            verify(&c, &[Gate::h(0), Gate::h(0)]);
            let c = circuit(1, vec![kind.numeric_gate(a, 0)]);
            assert_eq!(phase_fold_pauli(&c).gates, c.gates);
        }
        let c = circuit(
            1,
            vec![kind.numeric_gate(1e308, 0), kind.numeric_gate(1e308, 0)],
        );
        assert_eq!(phase_fold_pauli(&c).gates, c.gates);
        let c = circuit(
            1,
            vec![
                kind.numeric_gate(1e16, 0),
                kind.numeric_gate(1.0, 0),
                kind.numeric_gate(-1e16, 0),
            ],
        );
        // A lossy sum is refused; the intervening small contribution survives.
        assert_eq!(phase_fold_pauli(&c).gates, c.gates);
        let c = circuit(
            1,
            vec![kind.numeric_gate(1e16, 0), kind.numeric_gate(1e16, 0)],
        );
        verify(&c, &[kind.numeric_gate(2e16, 0)]);
    }
}

// Images of every matrix unit, separated by the final classical store. This
// independently checks channels, including measurements that overwrite bits.
fn channel(c: &Circuit, initial_store: u64) -> Vec<BTreeMap<u64, Vec<Vec<C>>>> {
    let d = 1 << c.num_qubits;
    let mut images = Vec::new();
    for i in 0..d {
        for j in 0..d {
            let mut rho = vec![vec![C::ZERO; d]; d];
            rho[i][j] = C::ONE;
            let mut branches = BTreeMap::from([(initial_store, rho)]);
            for gate in &c.gates {
                let mut next: BTreeMap<u64, Vec<Vec<C>>> = BTreeMap::new();
                for (store, rho) in branches {
                    match *gate {
                        Gate::measure { qubit, cbit } => {
                            let mask = 1 << (c.num_qubits - 1 - qubit as usize);
                            for outcome in 0..2 {
                                let store = (store & !(1 << cbit)) | (outcome << cbit);
                                let out = next
                                    .entry(store)
                                    .or_insert_with(|| vec![vec![C::ZERO; d]; d]);
                                for r in 0..d {
                                    for col in 0..d {
                                        if usize::from(r & mask != 0) == outcome as usize
                                            && usize::from(col & mask != 0) == outcome as usize
                                        {
                                            out[r][col] = out[r][col] + rho[r][col];
                                        }
                                    }
                                }
                            }
                        }
                        Gate::reset(q) => {
                            let mask = 1 << (c.num_qubits - 1 - q as usize);
                            let mut out = vec![vec![C::ZERO; d]; d];
                            for r in (0..d).filter(|r| r & mask == 0) {
                                for col in (0..d).filter(|col| col & mask == 0) {
                                    out[r][col] = rho[r][col] + rho[r | mask][col | mask];
                                }
                            }
                            next.insert(store, out);
                        }
                        _ => {
                            let u = circuit_unitary(&c.replacing_gates(vec![gate.clone()]));
                            let mut out = vec![vec![C::ZERO; d]; d];
                            for r in 0..d {
                                for col in 0..d {
                                    for k in 0..d {
                                        for l in 0..d {
                                            out[r][col] = out[r][col]
                                                + u[r][k] * rho[k][l] * u[col][l].conj();
                                        }
                                    }
                                }
                            }
                            next.insert(store, out);
                        }
                    }
                }
                branches = next;
            }
            images.push(branches);
        }
    }
    images
}

#[test]
fn native_rotations_preserve_full_channels_across_measurement_and_reset() {
    for kind in KINDS {
        for q in 0..2 {
            for barrier in [Gate::measure { qubit: q, cbit: 0 }, Gate::reset(q)] {
                let c = circuit(
                    2,
                    vec![
                        kind.numeric_gate(0.125, 0),
                        barrier.clone(),
                        kind.numeric_gate(0.25, 0),
                        Gate::measure { qubit: 1, cbit: 0 },
                    ],
                );
                let out = phase_fold_pauli(&c);
                let should_fold = q == 1
                    || (matches!(barrier, Gate::measure { .. })
                        && matches!(kind, RotationKind::Z | RotationKind::Phase));
                assert_eq!(out.gates.len(), c.gates.len() - usize::from(should_fold));
                for store in 0..4 {
                    let expected = channel(&c, store);
                    let actual = channel(&out, store);
                    for (e, a) in expected.iter().zip(&actual) {
                        assert_eq!(e.keys().collect::<Vec<_>>(), a.keys().collect::<Vec<_>>());
                        for (store, e) in e {
                            for (er, ar) in e.iter().zip(&a[store]) {
                                for (&x, &y) in er.iter().zip(ar) {
                                    assert!((x - y).norm_sq() < 1e-24, "{kind:?} {barrier:?}");
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn random_native_circuits_preserve_unitaries_with_colliding_sketches() {
    let mut rng = 0x4455_ee11_7788_9234u64;
    for case in 0..1_000 {
        let n = 3;
        let mut gates = Vec::new();
        for _ in 0..25 {
            let q = (next(&mut rng) % n as u64) as u32;
            let other = (q + 1 + (next(&mut rng) % 2) as u32) % 3;
            let theta = (next(&mut rng) % 33) as f64 / 8.0 - 2.0;
            gates.push(match next(&mut rng) % 24 {
                0 => Gate::h(q),
                1 => Gate::s(q),
                2 => Gate::sdg(q),
                3 => Gate::x(q),
                4 => Gate::y(q),
                5 => Gate::z(q),
                6 => Gate::sx(q),
                7 => Gate::swap(q, other),
                8 => Gate::cy {
                    control: q,
                    target: other,
                },
                9 => Gate::cnot {
                    control: q,
                    target: other,
                },
                10 => Gate::cz {
                    control: q,
                    target: other,
                },
                11 => Gate::t(q),
                12 => Gate::tdg(q),
                13 => Gate::rx_f64(theta, q).unwrap(),
                14 => Gate::ry_f64(theta, q).unwrap(),
                15 => Gate::rz_f64(theta, q).unwrap(),
                16 => Gate::p_f64(theta, q).unwrap(),
                17 => Gate::cp {
                    lambda: crate::angle::Angle::from_f64(theta).unwrap(),
                    control: q,
                    target: other,
                },
                18 => Gate::crx {
                    theta: crate::angle::Angle::from_f64(theta).unwrap(),
                    control: q,
                    target: other,
                },
                19 => Gate::cry {
                    theta: crate::angle::Angle::from_f64(theta).unwrap(),
                    control: q,
                    target: other,
                },
                20 => Gate::crz {
                    theta: crate::angle::Angle::from_f64(theta).unwrap(),
                    control: q,
                    target: other,
                },
                21 => Gate::ch {
                    control: q,
                    target: other,
                },
                22 => Gate::cswap {
                    control: q,
                    first: other,
                    second: 3 - q - other,
                },
                _ => Gate::ccx {
                    control1: q,
                    control2: other,
                    target: 3 - q - other,
                },
            });
        }
        let kind = KINDS[case % 4];
        gates.extend([
            kind.numeric_gate(0.125, 0),
            Gate::h(0),
            Gate::h(0),
            kind.numeric_gate(0.25, 0),
        ]);
        let c = circuit(n, gates);
        let out = phase_fold_pauli(&c);
        let collided = fold(&c, &vec![(0, 0); n], MAX_ATTEMPT_STEPS, true);
        let scalar = fold(&c, &vec![(0, 0); n], MAX_ATTEMPT_STEPS, false);
        assert_eq!(out.gates, collided.gates, "case {case}");
        assert_eq!(out.gates, scalar.gates, "case {case}");
        assert!(out.gates.len() < c.gates.len());
        assert!(circuits_equiv(&c, &out, 1e-10), "case {case}\n{c}\n{out}");
    }
}

impl RotationKind {
    fn numeric_gate(self, theta: f64, q: u32) -> Gate {
        self.gate(crate::angle::Angle::from_f64(theta).unwrap(), q)
    }
}

#[test]
fn checked_native_folding_refuses_rounded_sums_and_retains_symbolic_parts() {
    for kind in KINDS {
        let c = circuit(
            1,
            vec![kind.numeric_gate(0.1, 0), kind.numeric_gate(0.2, 0)],
        );
        assert_eq!(phase_fold_pauli(&c).gates, c.gates);
        let c = circuit(
            1,
            vec![
                kind.gate(crate::angle::Angle::pi_fraction(1, 8).unwrap(), 0),
                kind.numeric_gate(0.125, 0),
            ],
        );
        let output = phase_fold_pauli(&c);
        assert!(circuits_equiv(&c, &output, 1e-12));
        for gate in &output.gates {
            if let Some(angle) = gate.angle() {
                assert!(angle.components().is_some());
            }
        }
    }
}

#[test]
fn exact_native_clifford_rotations_update_each_physical_frame() {
    for kind in KINDS {
        for n in -8..=8 {
            let c = circuit(
                1,
                vec![
                    kind.numeric_gate(0.125, 0),
                    kind.gate(crate::angle::Angle::pi_fraction(n, 2).unwrap(), 0),
                    Gate::h(0),
                    Gate::h(0),
                    kind.numeric_gate(0.25, 0),
                ],
            );
            let output = phase_fold_pauli(&c);
            assert!(circuits_equiv(&c, &output, 1e-12), "{c}\n{output}");
            for sliced in [true, false] {
                assert_eq!(
                    fold(&c, &[(0, 0)], MAX_ATTEMPT_STEPS, sliced).gates,
                    output.gates
                );
            }
        }
    }
}
