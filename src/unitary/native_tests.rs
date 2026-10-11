use super::*;

fn circuit(n: usize, gates: &[Gate]) -> Circuit {
    Circuit::with_cbits(n, 0).replacing_gates(gates.to_vec())
}

fn assert_literal(actual: &[Vec<C>], expected: &[Vec<C>]) {
    assert_eq!(actual.len(), expected.len());
    for (row, values) in actual.iter().enumerate() {
        for (col, value) in values.iter().enumerate() {
            let error = (*value - expected[row][col]).norm_sq().sqrt();
            assert!(
                error < 2e-12,
                "entry ({row}, {col}): {value:?} != {:?}",
                expected[row][col]
            );
        }
    }
}

fn assert_unitary(matrix: &[Vec<C>]) {
    for a in 0..matrix.len() {
        for b in 0..matrix.len() {
            let dot = matrix
                .iter()
                .fold(C::ZERO, |sum, row| sum + row[a].conj() * row[b]);
            let expected = if a == b { C::ONE } else { C::ZERO };
            assert!((dot - expected).norm_sq() < 1e-24);
        }
    }
}

#[test]
fn native_y_and_sx_have_literal_standard_matrices() {
    let i = C::new(0.0, 1.0);
    assert_literal(
        &circuit_unitary(&circuit(1, &[Gate::y(0)])),
        &[vec![C::ZERO, i * -1.0], vec![i, C::ZERO]],
    );
    assert_literal(
        &circuit_unitary(&circuit(1, &[Gate::sx(0)])),
        &[
            vec![C::new(0.5, 0.5), C::new(0.5, -0.5)],
            vec![C::new(0.5, -0.5), C::new(0.5, 0.5)],
        ],
    );
}

#[test]
fn native_p_preserves_its_phase_convention() {
    for theta in [-4.0 * PI, -PI, -0.37, 0.0, PI / 4.0, PI, 2.0 * PI, 0.73] {
        let matrix = circuit_unitary(&circuit(1, &[Gate::p_f64(theta, 0).unwrap()]));
        assert_literal(
            &matrix,
            &[
                vec![C::ONE, C::ZERO],
                vec![C::ZERO, C::new(theta.cos(), theta.sin())],
            ],
        );
        assert_unitary(&matrix);
        let rz = circuit_unitary(&circuit(1, &[Gate::rz_f64(theta, 0).unwrap()]));
        let phase = C::polar(1.0, theta / 2.0);
        let expected: Vec<Vec<_>> = rz
            .iter()
            .map(|row| row.iter().map(|&entry| phase * entry).collect())
            .collect();
        assert_literal(&matrix, &expected);
    }
    let p = circuit_unitary(&circuit(1, &[Gate::p_f64(PI, 0).unwrap()]));
    let rz = circuit_unitary(&circuit(1, &[Gate::rz_f64(PI, 0).unwrap()]));
    assert!((p[0][0] - rz[0][0]).norm_sq() > 1.0);
}

#[test]
fn native_rotations_use_half_angles_and_correct_signs() {
    for theta in [
        -4.0 * PI,
        -2.0 * PI,
        -PI,
        -0.37,
        0.0,
        PI / 2.0,
        PI,
        2.0 * PI,
        7.31,
    ] {
        let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
        let rx = circuit_unitary(&circuit(1, &[Gate::rx_f64(theta, 0).unwrap()]));
        let ry = circuit_unitary(&circuit(1, &[Gate::ry_f64(theta, 0).unwrap()]));
        assert_literal(
            &rx,
            &[
                vec![C::new(c, 0.0), C::new(0.0, -s)],
                vec![C::new(0.0, -s), C::new(c, 0.0)],
            ],
        );
        assert_literal(
            &ry,
            &[
                vec![C::new(c, 0.0), C::new(-s, 0.0)],
                vec![C::new(s, 0.0), C::new(c, 0.0)],
            ],
        );
        assert_unitary(&rx);
        assert_unitary(&ry);
    }
    // These are literal -I at 2*pi, not I: a projective check would miss this.
    for gate in [
        Gate::rx_f64(2.0 * PI, 0).unwrap(),
        Gate::ry_f64(2.0 * PI, 0).unwrap(),
    ] {
        let matrix = circuit_unitary(&circuit(1, &[gate]));
        assert!((matrix[0][0] - C::new(-1.0, 0.0)).norm_sq() < 1e-24);
    }
}

#[test]
fn native_single_gates_match_decompositions_on_every_wire() {
    for n in 1..=4 {
        for q in 0..n as Qubit {
            let pairs = [
                (Gate::y(q), vec![Gate::sdg(q), Gate::x(q), Gate::s(q)]),
                (Gate::sx(q), vec![Gate::h(q), Gate::s(q), Gate::h(q)]),
                (
                    Gate::rx_f64(-0.73, q).unwrap(),
                    vec![Gate::h(q), Gate::rz_f64(-0.73, q).unwrap(), Gate::h(q)],
                ),
                (
                    Gate::ry_f64(1.39, q).unwrap(),
                    vec![
                        Gate::sdg(q),
                        Gate::h(q),
                        Gate::rz_f64(1.39, q).unwrap(),
                        Gate::h(q),
                        Gate::s(q),
                    ],
                ),
            ];
            for (native, old) in pairs {
                assert_literal(
                    &circuit_unitary(&circuit(n, &[native])),
                    &circuit_unitary(&circuit(n, &old)),
                );
            }
        }
    }
}

#[test]
fn native_inverses_and_sx_squared_are_literal_identities() {
    let identity = circuit_unitary(&Circuit::new(1));
    for theta in [-13.7, -PI, 0.0, PI / 4.0, 0.37, 2.0 * PI] {
        for pair in [
            [
                Gate::p_f64(theta, 0).unwrap(),
                Gate::p_f64(-theta, 0).unwrap(),
            ],
            [
                Gate::rx_f64(theta, 0).unwrap(),
                Gate::rx_f64(-theta, 0).unwrap(),
            ],
            [
                Gate::ry_f64(theta, 0).unwrap(),
                Gate::ry_f64(-theta, 0).unwrap(),
            ],
        ] {
            assert_literal(&circuit_unitary(&circuit(1, &pair)), &identity);
        }
    }
    assert_literal(
        &circuit_unitary(&circuit(1, &[Gate::y(0), Gate::y(0)])),
        &identity,
    );
    assert_literal(
        &circuit_unitary(&circuit(1, &[Gate::sx(0), Gate::sx(0)])),
        &circuit_unitary(&circuit(1, &[Gate::x(0)])),
    );
}

#[test]
fn native_swap_truth_table_on_every_pair_and_spectator() {
    for n in 2..=4 {
        for a in 0..n as Qubit {
            for b in 0..n as Qubit {
                if a == b {
                    continue;
                }
                let matrix = circuit_unitary(&circuit(n, &[Gate::swap(a, b)]));
                for col in 0..matrix.len() {
                    let read = |q| (col >> (n - 1 - q as usize)) & 1;
                    for (row, values) in matrix.iter().enumerate() {
                        let correct = (0..n as Qubit).all(|q| {
                            let source = if q == a {
                                b
                            } else if q == b {
                                a
                            } else {
                                q
                            };
                            (row >> (n - 1 - q as usize)) & 1 == read(source)
                        });
                        let expected = if correct { C::ONE } else { C::ZERO };
                        assert!((values[col] - expected).norm_sq() < 1e-24);
                    }
                }
                let cx = |control, target| Gate::cnot { control, target };
                assert_literal(
                    &matrix,
                    &circuit_unitary(&circuit(n, &[cx(a, b), cx(b, a), cx(a, b)])),
                );
                assert_literal(
                    &circuit_unitary(&circuit(n, &[Gate::swap(a, b), Gate::swap(b, a)])),
                    &circuit_unitary(&Circuit::new(n)),
                );
                assert_unitary(&matrix);
            }
        }
    }
}

#[test]
fn native_mixed_circuits_preserve_literal_matrix_under_expansion() {
    // Noncommuting contexts, entanglement, and nonadjacent endpoints exercise
    // application to a general matrix rather than only the identity input.
    for k in -8..=8 {
        let theta = k as f64 * 0.371;
        let prefix = [
            Gate::h(0),
            Gate::cnot {
                control: 0,
                target: 2,
            },
            Gate::t(1),
        ];
        let mut native = circuit(3, &prefix);
        native.gates.extend([
            Gate::p_f64(theta, 1).unwrap(),
            Gate::rx_f64(theta, 2).unwrap(),
            Gate::y(0),
            Gate::swap(2, 0),
            Gate::ry_f64(-theta, 1).unwrap(),
            Gate::sx(2),
        ]);
        let mut expanded = circuit(3, &prefix);
        expanded.gates.extend([
            Gate::p_f64(theta, 1).unwrap(),
            Gate::h(2),
            Gate::rz_f64(theta, 2).unwrap(),
            Gate::h(2),
            Gate::sdg(0),
            Gate::x(0),
            Gate::s(0),
            Gate::cnot {
                control: 2,
                target: 0,
            },
            Gate::cnot {
                control: 0,
                target: 2,
            },
            Gate::cnot {
                control: 2,
                target: 0,
            },
            Gate::sdg(1),
            Gate::h(1),
            Gate::rz_f64(-theta, 1).unwrap(),
            Gate::h(1),
            Gate::s(1),
            Gate::h(2),
            Gate::s(2),
            Gate::h(2),
        ]);
        assert_literal(&circuit_unitary(&native), &circuit_unitary(&expanded));
        assert_unitary(&circuit_unitary(&native));
    }
}

#[test]
fn native_gates_are_preserved_safely_by_existing_passes() {
    use crate::pass::Pass;
    let passes: [&dyn Pass; 3] = [
        &crate::cancel::CancelGates,
        &crate::phase_fold_rand::PhaseFoldRand,
        &crate::cnot_min::CnotMin::default(),
    ];
    let cx = Gate::cnot {
        control: 0,
        target: 1,
    };
    for gate in [
        Gate::p_f64(0.37, 0).unwrap(),
        Gate::y(0),
        Gate::sx(0),
        Gate::rx_f64(-0.73, 0).unwrap(),
        Gate::ry_f64(1.39, 0).unwrap(),
        Gate::swap(0, 1),
        Gate::cy {
            control: 0,
            target: 1,
        },
        Gate::ch {
            control: 1,
            target: 0,
        },
        Gate::cp {
            lambda: crate::angle::Angle::from_f64(0.37).unwrap(),
            control: 0,
            target: 1,
        },
        Gate::crx {
            theta: crate::angle::Angle::from_f64(-0.73).unwrap(),
            control: 1,
            target: 0,
        },
        Gate::cry {
            theta: crate::angle::Angle::from_f64(1.39).unwrap(),
            control: 0,
            target: 1,
        },
        Gate::crz {
            theta: crate::angle::Angle::from_f64(0.37).unwrap(),
            control: 1,
            target: 0,
        },
        Gate::cswap {
            control: 0,
            first: 1,
            second: 2,
        },
    ] {
        for (prefix, suffix) in [
            (vec![Gate::t(0)], vec![Gate::tdg(0)]),
            (vec![cx.clone()], vec![cx.clone()]),
            (vec![Gate::s(0), Gate::h(1)], vec![Gate::sdg(0), Gate::t(1)]),
        ] {
            let mut input = circuit(3, &prefix);
            input.gates.push(gate.clone());
            input.gates.extend(suffix);
            for pass in passes {
                let actual = pass.run(&input);
                assert!(
                    actual.gates.contains(&gate),
                    "{} lost {gate:?}",
                    pass.name()
                );
                assert!(
                    circuits_equiv(&input, &actual, 1e-10),
                    "{} changed {gate:?}",
                    pass.name()
                );
            }
        }
        // PBC transfer rules are deferred; reject rather than silently drop.
        assert_eq!(
            crate::pbc::to_pbc(&circuit(3, &[gate.clone()]), None).unwrap_err(),
            crate::pbc::PbcError::InvalidInput {
                index: 0,
                cause: Box::new(crate::pbc::PbcError::UnsupportedGate(gate.kind()))
            }
        );
    }
}

#[test]
fn native_swap_counts_as_one_two_qubit_gate_and_occupies_both_wires() {
    let input = circuit(
        3,
        &[
            Gate::y(0),
            Gate::sx(2),
            Gate::swap(0, 2),
            Gate::rx_f64(0.37, 2).unwrap(),
        ],
    );
    assert_eq!(crate::pass::count_2q(&input), 1);
    assert_eq!(crate::pass::depth(&input), 3);
    let metrics = crate::optimize::Metrics::of(&input);
    assert_eq!(metrics.two_qubit, 1);
    assert_eq!(metrics.depth, 3);
}

#[test]
fn native_serialization_uses_qasm2_library_operations() {
    let input = circuit(
        2,
        &[
            Gate::p_f64(0.37, 0).unwrap(),
            Gate::y(1),
            Gate::sx(0),
            Gate::rx_f64(-0.73, 1).unwrap(),
            Gate::ry_f64(1.39, 0).unwrap(),
            Gate::swap(0, 1),
        ],
    );
    let output = input.to_qasm();
    for operation in [
        "u1(0.37) q[0];",
        "y q[1];",
        "h q[0];\ns q[0];\nh q[0];",
        "rx(-0.73) q[1];",
        "ry(1.39) q[0];",
        "cx q[0],q[1];\ncx q[1],q[0];\ncx q[0],q[1];",
    ] {
        assert!(
            output.contains(operation),
            "missing {operation} in {output}"
        );
    }
    let sx = circuit(1, &[Gate::sx(0)]);
    let parsed = Circuit::from_qasm(&sx.to_qasm()).unwrap();
    assert_literal(&circuit_unitary(&sx), &circuit_unitary(&parsed));
}

fn controlled_gates(theta: f64, control: Qubit, target: Qubit) -> [Gate; 6] {
    [
        Gate::cy { control, target },
        Gate::ch { control, target },
        Gate::cp {
            lambda: crate::angle::Angle::from_f64(theta).unwrap(),
            control,
            target,
        },
        Gate::crx {
            theta: crate::angle::Angle::from_f64(theta).unwrap(),
            control,
            target,
        },
        Gate::cry {
            theta: crate::angle::Angle::from_f64(theta).unwrap(),
            control,
            target,
        },
        Gate::crz {
            theta: crate::angle::Angle::from_f64(theta).unwrap(),
            control,
            target,
        },
    ]
}

#[test]
fn controlled_native_gates_have_literal_blocks_on_every_wire() {
    for theta in [-4.0 * PI, -0.73, 0.0, PI / 2.0, 2.0 * PI, 7.31] {
        let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
        let h = 1.0 / 2.0_f64.sqrt();
        let blocks = [
            [[C::ZERO, C::new(0.0, -1.0)], [C::new(0.0, 1.0), C::ZERO]],
            [
                [C::new(h, 0.0), C::new(h, 0.0)],
                [C::new(h, 0.0), C::new(-h, 0.0)],
            ],
            [[C::ONE, C::ZERO], [C::ZERO, C::polar(1.0, theta)]],
            [
                [C::new(c, 0.0), C::new(0.0, -s)],
                [C::new(0.0, -s), C::new(c, 0.0)],
            ],
            [
                [C::new(c, 0.0), C::new(-s, 0.0)],
                [C::new(s, 0.0), C::new(c, 0.0)],
            ],
            [
                [C::polar(1.0, -theta / 2.0), C::ZERO],
                [C::ZERO, C::polar(1.0, theta / 2.0)],
            ],
        ];
        for n in 2..=4 {
            for control in 0..n as Qubit {
                for target in 0..n as Qubit {
                    if control == target {
                        continue;
                    }
                    for (gate, block) in controlled_gates(theta, control, target)
                        .into_iter()
                        .zip(blocks)
                    {
                        let actual = circuit_unitary(&circuit(n, &[gate]));
                        let dim = 1 << n;
                        let cb = 1 << (n - 1 - control as usize);
                        let tb = 1 << (n - 1 - target as usize);
                        let mut expected = vec![vec![C::ZERO; dim]; dim];
                        for col in 0..dim {
                            if col & cb == 0 {
                                expected[col][col] = C::ONE;
                            } else {
                                let input_bit = usize::from(col & tb != 0);
                                expected[col & !tb][col] = block[0][input_bit];
                                expected[col | tb][col] = block[1][input_bit];
                            }
                        }
                        assert_literal(&actual, &expected);
                        assert_unitary(&actual);
                    }
                }
            }
        }
    }
}

#[test]
fn controlled_rotation_expansions_preserve_literal_branch_phases() {
    for theta in [-4.0 * PI, -2.0 * PI, -0.73, 0.0, PI / 4.0, 2.0 * PI, 7.31] {
        for (control, target) in [(0, 1), (1, 0), (0, 2), (2, 1)] {
            for gate in controlled_gates(theta, control, target).into_iter().skip(2) {
                let terms = crate::decompose::expand_controlled_rotation(&gate)
                    .unwrap()
                    .unwrap();
                assert_literal(
                    &circuit_unitary(&circuit(3, &[gate])),
                    &circuit_unitary(&circuit(3, &terms)),
                );
            }
        }
    }
    // A rotation by 2*pi contributes -I only on the active control branch.
    for gate in controlled_gates(2.0 * PI, 0, 1).into_iter().skip(3) {
        assert_literal(
            &circuit_unitary(&circuit(2, &[gate])),
            &circuit_unitary(&circuit(2, &[Gate::z(0)])),
        );
    }
    assert_literal(
        &circuit_unitary(&circuit(
            2,
            &[Gate::cp {
                lambda: crate::angle::Angle::from_f64(2.0 * PI).unwrap(),
                control: 0,
                target: 1,
            }],
        )),
        &circuit_unitary(&Circuit::new(2)),
    );
}

#[test]
fn cswap_truth_table_and_decomposition_cover_all_operand_orders() {
    for n in 3..=4 {
        for control in 0..n as Qubit {
            for first in 0..n as Qubit {
                for second in 0..n as Qubit {
                    if control == first || control == second || first == second {
                        continue;
                    }
                    let gate = Gate::cswap {
                        control,
                        first,
                        second,
                    };
                    let dim = 1 << n;
                    let mut expected = vec![vec![C::ZERO; dim]; dim];
                    for col in 0..dim {
                        let row = if col & (1 << (n - 1 - control as usize)) != 0
                            && ((col >> (n - 1 - first as usize)) & 1)
                                != ((col >> (n - 1 - second as usize)) & 1)
                        {
                            col ^ (1 << (n - 1 - first as usize)) ^ (1 << (n - 1 - second as usize))
                        } else {
                            col
                        };
                        expected[row][col] = C::ONE;
                    }
                    assert_literal(&circuit_unitary(&circuit(n, &[gate.clone()])), &expected);
                    let terms = [
                        Gate::cnot {
                            control: first,
                            target: second,
                        },
                        Gate::ccx {
                            control1: control,
                            control2: second,
                            target: first,
                        },
                        Gate::cnot {
                            control: first,
                            target: second,
                        },
                    ];
                    assert_literal(
                        &circuit_unitary(&circuit(n, &[gate])),
                        &circuit_unitary(&circuit(n, &terms)),
                    );
                }
            }
        }
    }
}
