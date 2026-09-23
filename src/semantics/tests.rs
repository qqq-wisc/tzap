use super::test_support::approximate;
use super::*;
use crate::pbc::PauliAngle;

mod pairs;

fn gates(n: usize, gates: Vec<Gate>) -> Circuit {
    Circuit {
        num_qubits: n,
        num_cbits: 0,
        gates,
    }
}

fn eval(c: &PbcCircuit) -> Matrix {
    pbc_unitary(c, Limits::default()).unwrap()
}
fn check(c: &PbcCircuit, expected: Vec<Gate>) {
    let circuit = gates(c.num_qubits(), expected);
    test_support::assert_equivalent(&circuit, c);
    // Every hand-written pair also exercises the automatic converter.
    let converted = crate::pbc::to_pbc(&circuit).unwrap();
    test_support::assert_equivalent(&circuit, &converted);
    assert!(eval(c).equivalent_up_to_global_phase(&eval(&converted)));
}

/// Small declarative fixtures: one letter per qubit, angle in units of pi/8.
fn rotations(n: usize, layers: &[(&str, i64)]) -> PbcCircuit {
    let mut c = PbcCircuit::new(n, 0);
    for &(word, k) in layers {
        assert_eq!(word.len(), n);
        let mut product = c.identity().as_ref();
        for (q, p) in word.chars().enumerate() {
            let p = match p {
                'I' => Pauli::I,
                'X' => Pauli::X,
                'Y' => Pauli::Y,
                'Z' => Pauli::Z,
                _ => panic!("invalid fixture Pauli"),
            };
            let single = c.single(q as u32, p).unwrap();
            product = c.product(product, single.as_ref()).unwrap();
        }
        let axis = c.hermitian_axis(product, 10_000).unwrap();
        c.rotate(axis, PauliAngle::new(k)).unwrap();
    }
    c
}

macro_rules! reorder_case {
    ($name:ident, $n:expr, $a:expr, $b:expr, $equivalent:expr) => {
        #[test]
        fn $name() {
            let a = eval(&rotations($n, $a));
            let b = eval(&rotations($n, $b));
            assert_eq!(a.equivalent_up_to_global_phase(&b), $equivalent);
            if $equivalent {
                assert_eq!(a, b);
            }
        }
    };
}

reorder_case!(
    commute_disjoint_qubits,
    2,
    &[("XI", 1), ("IZ", 2)],
    &[("IZ", 2), ("XI", 1)],
    true
);
reorder_case!(
    commute_same_axis,
    1,
    &[("Y", 1), ("Y", 3)],
    &[("Y", 3), ("Y", 1)],
    true
);
reorder_case!(
    commute_overlapping_z_strings,
    3,
    &[("ZZI", 1), ("IZZ", -1)],
    &[("IZZ", -1), ("ZZI", 1)],
    true
);
reorder_case!(
    commute_two_local_anticommutations,
    2,
    &[("XX", 1), ("ZZ", 2)],
    &[("ZZ", 2), ("XX", 1)],
    true
);
reorder_case!(
    commute_with_y_factors,
    2,
    &[("XY", -1), ("YX", 3)],
    &[("YX", 3), ("XY", -1)],
    true
);
reorder_case!(
    commute_whole_layers,
    3,
    &[("ZII", 1), ("IZZ", 2), ("ZZI", -1), ("IIZ", 3)],
    &[("ZZI", -1), ("IIZ", 3), ("ZII", 1), ("IZZ", 2)],
    true
);
reorder_case!(
    cannot_commute_x_and_z_rotations,
    1,
    &[("X", 1), ("Z", 1)],
    &[("Z", 1), ("X", 1)],
    false
);
reorder_case!(
    cannot_commute_one_overlap,
    3,
    &[("XXI", 1), ("IZZ", 1)],
    &[("IZZ", 1), ("XXI", 1)],
    false
);
reorder_case!(
    cannot_commute_y_and_z_rotations,
    1,
    &[("Y", 2), ("Z", 1)],
    &[("Z", 1), ("Y", 2)],
    false
);

macro_rules! gate_pbc_case {
    ($name:ident, $n:expr, $gates:expr, $rotations:expr) => {
        #[test]
        fn $name() {
            check(&rotations($n, $rotations), $gates);
        }
    };
}

gate_pbc_case!(t_is_z_rotation, 1, vec![Gate::t(0)], &[("Z", 1)]);
gate_pbc_case!(
    tdg_is_negative_z_rotation,
    1,
    vec![Gate::tdg(0)],
    &[("Z", -1)]
);
gate_pbc_case!(s_is_z_quarter_rotation, 1, vec![Gate::s(0)], &[("Z", 2)]);
gate_pbc_case!(x_is_x_half_rotation, 1, vec![Gate::x(0)], &[("X", 4)]);
gate_pbc_case!(
    h_is_three_pauli_rotations,
    1,
    vec![Gate::h(0)],
    &[("Z", 2), ("X", 2), ("Z", 2)]
);
gate_pbc_case!(
    h_t_h_is_x_rotation,
    1,
    vec![Gate::h(0), Gate::t(0), Gate::h(0)],
    &[("X", 1)]
);
gate_pbc_case!(
    sdg_h_t_h_s_is_y_rotation,
    1,
    vec![Gate::sdg(0), Gate::h(0), Gate::t(0), Gate::h(0), Gate::s(0)],
    &[("Y", 1)]
);
gate_pbc_case!(
    cz_is_three_joint_rotations,
    2,
    vec![Gate::cz {
        control: 0,
        target: 1
    }],
    &[("ZI", 2), ("IZ", 2), ("ZZ", -2)]
);
gate_pbc_case!(
    cx_is_three_joint_rotations,
    2,
    vec![Gate::cnot {
        control: 0,
        target: 1
    }],
    &[("ZI", 2), ("IX", 2), ("ZX", -2)]
);
gate_pbc_case!(
    reverse_cx_is_three_joint_rotations,
    2,
    vec![Gate::cnot {
        control: 1,
        target: 0
    }],
    &[("IZ", 2), ("XI", 2), ("XZ", -2)]
);
gate_pbc_case!(
    cx_t_cx_is_parity_rotation,
    2,
    vec![
        Gate::cnot {
            control: 0,
            target: 1
        },
        Gate::t(1),
        Gate::cnot {
            control: 0,
            target: 1
        }
    ],
    &[("ZZ", 1)]
);
gate_pbc_case!(
    ccz_is_seven_joint_rotations,
    3,
    vec![Gate::ccz {
        control1: 0,
        control2: 1,
        target: 2
    }],
    &[
        ("ZII", 1),
        ("IZI", 1),
        ("IIZ", 1),
        ("ZZI", -1),
        ("ZIZ", -1),
        ("IZZ", -1),
        ("ZZZ", 1)
    ]
);
gate_pbc_case!(
    ccx_is_seven_joint_rotations,
    3,
    vec![Gate::ccx {
        control1: 0,
        control2: 1,
        target: 2
    }],
    &[
        ("ZII", 1),
        ("IZI", 1),
        ("IIX", 1),
        ("ZZI", -1),
        ("ZIX", -1),
        ("IZX", -1),
        ("ZZX", 1)
    ]
);

#[test]
fn h_then_t_keeps_h_in_the_output_frame() {
    let mut c = rotations(1, &[("X", 1)]);
    c.push_output_clifford(Gate::h(0)).unwrap();
    check(&c, vec![Gate::h(0), Gate::t(0)]);
}

#[test]
fn cx_then_target_t_keeps_cx_in_the_output_frame() {
    let mut c = rotations(2, &[("ZZ", 1)]);
    let cx = Gate::cnot {
        control: 0,
        target: 1,
    };
    c.push_output_clifford(cx.clone()).unwrap();
    check(&c, vec![cx, Gate::t(1)]);
}

#[test]
fn commuting_layers_remain_equivalent_with_a_common_suffix() {
    let mut a = rotations(2, &[("XX", 1), ("ZZ", -1)]);
    let mut b = rotations(2, &[("ZZ", -1), ("XX", 1)]);
    for gate in [
        Gate::h(0),
        Gate::cnot {
            control: 0,
            target: 1,
        },
    ] {
        a.push_output_clifford(gate.clone()).unwrap();
        b.push_output_clifford(gate).unwrap();
    }
    assert_eq!(eval(&a), eval(&b));
}

#[test]
fn moving_a_rotation_across_a_clifford_requires_changing_its_axis() {
    let mut correct = rotations(1, &[("X", 1)]);
    correct.push_output_clifford(Gate::h(0)).unwrap();
    let mut wrong = rotations(1, &[("Z", 1)]);
    wrong.push_output_clifford(Gate::h(0)).unwrap();
    assert!(!eval(&correct).equivalent_up_to_global_phase(&eval(&wrong)));
}

#[test]
fn equal_axes_merge_and_inverse_rotations_cancel() {
    assert_eq!(
        eval(&rotations(2, &[("XY", 1), ("XY", 2)])),
        eval(&rotations(2, &[("XY", 3)]))
    );
    assert_eq!(
        eval(&rotations(2, &[("YZ", 1), ("YZ", -1)])),
        Matrix::identity(4)
    );
}

#[test]
fn exact_scalar_field_relations() {
    let one = Scalar::integer(1);
    assert_eq!(Scalar::i().mul(&Scalar::i()), one.neg());
    assert_eq!(
        Scalar::inv_sqrt_two().mul(&Scalar::inv_sqrt_two()),
        one.half()
    );
    for k in 0..8 {
        let w = Scalar::omega(k);
        assert_eq!(w.mul(&w.conj()), one);
        for j in 0..8 {
            assert_eq!(w.mul(&Scalar::omega(j)), Scalar::omega(k + j));
        }
    }
}

#[test]
fn empty_circuits_and_zero_qubit_identity_rotations() {
    for n in 0..4 {
        assert_eq!(eval(&PbcCircuit::new(n, 0)), Matrix::identity(1 << n));
    }
    let mut c = PbcCircuit::new(0, 0);
    c.rotate(c.identity(), PauliAngle::new(1)).unwrap();
    assert_eq!(eval(&c), Matrix::identity(1));
    c.rotate(c.identity().negated(), PauliAngle::new(1))
        .unwrap();
    assert_eq!(eval(&c), Matrix::identity(1).scale(&Scalar::omega(1)));
    assert!(eval(&c).equivalent_up_to_global_phase(&Matrix::identity(1)));
}

#[test]
fn every_angle_axis_and_sign_matches_gate_decompositions() {
    for p in [Pauli::X, Pauli::Y, Pauli::Z] {
        for negative in [false, true] {
            for k in 0..8 {
                let mut c = PbcCircuit::new(1, 0);
                let mut axis = c.single(0, p).unwrap();
                if negative {
                    axis = axis.negated();
                }
                c.rotate(axis, PauliAngle::new(k)).unwrap();
                let mut expected = match p {
                    Pauli::X => vec![Gate::h(0)],
                    Pauli::Y => vec![Gate::sdg(0), Gate::h(0)],
                    _ => vec![],
                };
                expected.extend((0..k).map(|_| if negative { Gate::tdg(0) } else { Gate::t(0) }));
                match p {
                    Pauli::X => expected.push(Gate::h(0)),
                    Pauli::Y => expected.extend([Gate::h(0), Gate::s(0)]),
                    _ => (),
                }
                check(&c, expected);
            }
        }
    }
}

#[test]
fn ordered_dag_products_and_reference_phases() {
    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    let z = c.z(0).unwrap();
    let xz = c.product(x.as_ref(), z.as_ref()).unwrap();
    // iXZ = Y: imaginary intermediate values must survive interpretation.
    let y = c.hermitian_axis(xz.scaled(Phase::I), 100).unwrap();
    c.rotate(y, PauliAngle::new(1)).unwrap();
    check(
        &c,
        vec![Gate::sdg(0), Gate::h(0), Gate::t(0), Gate::h(0), Gate::s(0)],
    );
    let inverse = c.hermitian_axis(xz.scaled(Phase::MinusI), 100).unwrap();
    c.rotate(inverse, PauliAngle::new(1)).unwrap();
    check(&c, vec![]);
}

#[test]
fn all_two_qubit_strings_angles_and_signs_match_parity_circuits() {
    for p0 in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
        for p1 in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
            let factors = [p0, p1];
            let support: Vec<_> = (0..2).filter(|&q| factors[q] != Pauli::I).collect();
            if support.is_empty() {
                continue;
            }
            for negative in [false, true] {
                for k in 0..8 {
                    let mut c = PbcCircuit::new(2, 0);
                    let a = c.single(0, p0).unwrap();
                    let b = c.single(1, p1).unwrap();
                    let product = c.product(a.as_ref(), b.as_ref()).unwrap();
                    let mut axis = c.hermitian_axis(product, 100).unwrap();
                    if negative {
                        axis = axis.negated();
                    }
                    c.rotate(axis, PauliAngle::new(k)).unwrap();
                    // Independent gate construction: change each basis to Z,
                    // compute parity on the last support qubit, apply T^k,
                    // then uncompute parity and the basis changes.
                    let mut prefix = vec![];
                    let mut undo = vec![];
                    for (q, p) in factors.iter().enumerate() {
                        let q = q as u32;
                        if *p == Pauli::Y {
                            prefix.push(Gate::sdg(q));
                            undo.push(Gate::s(q));
                        }
                        if matches!(p, Pauli::X | Pauli::Y) {
                            prefix.push(Gate::h(q));
                            undo.push(Gate::h(q));
                        }
                    }
                    let target = *support.last().unwrap() as u32;
                    let cx = Gate::cnot {
                        control: 0,
                        target: 1,
                    };
                    if support.len() == 2 {
                        prefix.push(cx.clone());
                    }
                    prefix.extend((0..k).map(|_| {
                        if negative {
                            Gate::tdg(target)
                        } else {
                            Gate::t(target)
                        }
                    }));
                    if support.len() == 2 {
                        prefix.push(cx);
                    }
                    prefix.extend(undo.into_iter().rev());
                    check(&c, prefix);
                }
            }
        }
    }
}

#[test]
fn joint_axes_and_suffix_order() {
    let mut c = PbcCircuit::new(3, 0);
    let z0 = c.z(0).unwrap();
    let z2 = c.z(2).unwrap();
    let product = c.product(z0.as_ref(), z2.as_ref()).unwrap();
    let axis = c.hermitian_axis(product, 100).unwrap();
    // Builder call order does not move the suffix before rotations.
    c.push_output_clifford(Gate::h(0)).unwrap();
    c.rotate(axis, PauliAngle::new(1)).unwrap();
    let cx = Gate::cnot {
        control: 0,
        target: 2,
    };
    check(&c, vec![cx.clone(), Gate::t(2), cx, Gate::h(0)]);
    let wrong =
        circuit_unitary(&gates(3, vec![Gate::h(0), Gate::t(2)]), Limits::default()).unwrap();
    assert!(!eval(&c).equivalent_up_to_global_phase(&wrong));
}

#[test]
fn negative_joint_product_and_shared_nodes() {
    let mut c = PbcCircuit::new(2, 0);
    let x0 = c.x(0).unwrap();
    let x1 = c.x(1).unwrap();
    let z0 = c.z(0).unwrap();
    let z1 = c.z(1).unwrap();
    let xx = c.product(x0.as_ref(), x1.as_ref()).unwrap();
    let zz = c.product(z0.as_ref(), z1.as_ref()).unwrap();
    let minus_yy = c.product(xx, zz).unwrap();
    let axis = c.hermitian_axis(minus_yy, 100).unwrap();
    let y0 = c.y(0).unwrap();
    let y1 = c.y(1).unwrap();
    let yy = c.product(y0.as_ref(), y1.as_ref()).unwrap();
    let inverse = c.hermitian_axis(yy, 100).unwrap();
    c.rotate(axis, PauliAngle::new(1)).unwrap();
    c.rotate(inverse, PauliAngle::new(1)).unwrap();
    check(&c, vec![]);
    let squared = c.product(minus_yy, minus_yy).unwrap();
    let identity = c.hermitian_axis(squared, 100).unwrap();
    c.rotate(identity, PauliAngle::new(3)).unwrap();
    check(&c, vec![]);
}

#[test]
fn noncommuting_rotations_keep_execution_order() {
    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    let z = c.z(0).unwrap();
    c.rotate(x, PauliAngle::new(1)).unwrap();
    c.rotate(z, PauliAngle::new(2)).unwrap();
    check(&c, vec![Gate::h(0), Gate::t(0), Gate::h(0), Gate::s(0)]);
    let reversed = circuit_unitary(
        &gates(1, vec![Gate::s(0), Gate::h(0), Gate::t(0), Gate::h(0)]),
        Limits::default(),
    )
    .unwrap();
    assert!(!eval(&c).equivalent_up_to_global_phase(&reversed));
}

#[test]
fn gate_oracle_matches_existing_independent_numeric_interpreter() {
    let mut cases = vec![];
    for q in 0..3 {
        cases.extend([
            Gate::x(q),
            Gate::z(q),
            Gate::h(q),
            Gate::s(q),
            Gate::sdg(q),
            Gate::t(q),
            Gate::tdg(q),
        ]);
    }
    for a in 0..3 {
        for b in 0..3 {
            if a != b {
                cases.extend([
                    Gate::cnot {
                        control: a,
                        target: b,
                    },
                    Gate::cz {
                        control: a,
                        target: b,
                    },
                ]);
                let t = 3 - a - b;
                cases.extend([
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
                ]);
            }
        }
    }
    for gate in cases {
        let c = gates(3, vec![Gate::h(1), gate, Gate::t(0), Gate::s(2)]);
        let exact = circuit_unitary(&c, Limits::default()).unwrap();
        let numeric = crate::unitary::circuit_unitary(&c);
        for (r, row) in numeric.iter().enumerate() {
            for (col, entry) in row.iter().enumerate() {
                let (a, b) = approximate(exact.get(r, col));
                let (x, y) = entry.components();
                assert!((a - x).abs() < 1e-12 && (b - y).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn all_suffix_gates_have_unitary_semantics() {
    let suffix = vec![
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
    ];
    let mut c = PbcCircuit::new(2, 0);
    for g in &suffix {
        c.push_output_clifford(g.clone()).unwrap();
    }
    check(&c, suffix);
}

#[test]
fn global_phase_comparison_rejects_rescaling_and_distinct_gates() {
    let i = Matrix::identity(2);
    assert!(i.equivalent_up_to_global_phase(&i.scale(&Scalar::omega(3))));
    assert!(!i.equivalent_up_to_global_phase(&i.scale(&Scalar::integer(2))));
    assert!(!i.equivalent_up_to_global_phase(&pauli(1, 0, Pauli::Z)));
    assert!(!i.equivalent_up_to_global_phase(&Matrix::identity(1)));
}

#[test]
fn nonunitary_and_unsupported_input_is_rejected() {
    let mut c = PbcCircuit::new(1, 0);
    let z = c.z(0).unwrap();
    c.rotate(z, PauliAngle::new(1)).unwrap();
    let m = c.measure(z, None).unwrap();
    c.conditional_rotate(z, PauliAngle::new(2), m).unwrap();
    assert_eq!(
        pbc_unitary(&c, Limits::default()),
        Err(Error::UnsupportedOperation { index: 1 })
    );
    for gate in [
        Gate::measure { qubit: 0, cbit: 0 },
        Gate::reset(0),
        Gate::rz(0.1, 0),
    ] {
        assert_eq!(
            circuit_unitary(&gates(1, vec![gate]), Limits::default()),
            Err(Error::UnsupportedOperation { index: 0 })
        );
    }
    for gate in [
        Gate::h(1),
        Gate::cnot {
            control: 0,
            target: 0,
        },
    ] {
        assert_eq!(
            circuit_unitary(&gates(1, vec![gate]), Limits::default()),
            Err(Error::InvalidOperand { index: 0 })
        );
    }
}

#[test]
fn limits_fail_before_dense_allocation() {
    assert_eq!(
        pbc_unitary(&PbcCircuit::new(usize::MAX, 0), Limits::default()),
        Err(Error::LimitExceeded)
    );
    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    c.rotate(x, PauliAngle::new(1)).unwrap();
    for limits in [
        Limits {
            max_matrix_cells: 0,
            ..Limits::default()
        },
        Limits {
            max_multiply_terms: 0,
            ..Limits::default()
        },
        Limits {
            max_operations: 0,
            ..Limits::default()
        },
    ] {
        assert_eq!(pbc_unitary(&c, limits), Err(Error::LimitExceeded));
    }
}
