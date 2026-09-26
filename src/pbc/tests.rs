use super::*;
use crate::circuit::Gate;

#[test]
fn angles_are_exact_modulo_global_phase() {
    for k in -64..64 {
        let angle = PauliAngle::new(k);
        assert_eq!(angle.eighths(), k.rem_euclid(8) as u8);
        assert_eq!(angle + -angle, PauliAngle::default());
        assert_eq!(
            i64::from(angle.signed_eighths()).rem_euclid(8),
            k.rem_euclid(8)
        );
    }
    assert_eq!(PauliAngle::new(i64::MIN).eighths(), 0);
    assert_eq!(PauliAngle::new(-1).to_string(), "-pi/8");
}

#[test]
fn ordered_products_keep_imaginary_phases() {
    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    let z = c.z(0).unwrap();
    let xz = c.product(x.as_ref(), z.as_ref()).unwrap();
    assert_eq!(
        c.expand(xz, 100).unwrap(),
        ExpandedPauli {
            phase: Phase::MinusI,
            factors: vec![Pauli::Y]
        }
    );
    assert_eq!(c.hermitian_axis(xz, 100), Err(PbcError::NonHermitianAxis));
    let y = c.hermitian_axis(xz.scaled(Phase::I), 100).unwrap();
    assert_eq!(c.expand(y.as_ref(), 100).unwrap().phase, Phase::One);
    let xx = c.product(x.as_ref(), x.as_ref()).unwrap();
    assert_eq!(c.expand(xx, 100).unwrap().factors, vec![Pauli::I]);
    assert_eq!(
        c.expand(y.negated().as_ref(), 100).unwrap().phase,
        Phase::MinusOne
    );
}

#[test]
fn pauli_multiplication_is_associative() {
    for a in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
        for b in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
            for c in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
                let (p, ab) = a.times(b);
                let (q, abc) = ab.times(c);
                let (r, bc) = b.times(c);
                let (s, abc2) = a.times(bc);
                assert_eq!((p * q, abc), (r * s, abc2));
            }
        }
    }
}

#[test]
fn commuting_multi_qubit_products_can_have_negative_signs() {
    let mut c = PbcCircuit::new(2, 0);
    let x0 = c.x(0).unwrap();
    let x1 = c.x(1).unwrap();
    let z0 = c.z(0).unwrap();
    let z1 = c.z(1).unwrap();
    let xx = c.product(x0.as_ref(), x1.as_ref()).unwrap();
    let zz = c.product(z0.as_ref(), z1.as_ref()).unwrap();
    let product = c.product(xx, zz).unwrap();
    let axis = c.hermitian_axis(product, 100).unwrap();
    assert_eq!(
        c.expand(axis.as_ref(), 100).unwrap(),
        ExpandedPauli {
            phase: Phase::MinusOne,
            factors: vec![Pauli::Y, Pauli::Y],
        }
    );
    c.rotate(axis, PauliAngle::new(-1)).unwrap();
    assert!(c.to_ascii().unwrap().contains("R(-pi/8,-)"));
}

#[test]
fn handles_cannot_cross_circuits() {
    let mut a = PbcCircuit::new(1, 0);
    let mut b = PbcCircuit::new(1, 0);
    let x = a.x(0).unwrap();
    let y = b.y(0).unwrap();
    let m = a.measure(x, None).unwrap();
    assert_eq!(
        b.product(x.as_ref(), y.as_ref()),
        Err(PbcError::ForeignPauli)
    );
    assert_eq!(b.rotate(x, PauliAngle::new(1)), Err(PbcError::ForeignPauli));
    assert_eq!(b.measure(x, None), Err(PbcError::ForeignPauli));
    assert_eq!(b.expand(x.as_ref(), 100), Err(PbcError::ForeignPauli));
    assert_eq!(b.conditional_pauli(y, m), Err(PbcError::UnknownMeasurement));
}

#[test]
fn measurements_have_immutable_ids_even_when_classical_bits_are_reused() {
    let mut c = PbcCircuit::new(1, 1);
    let z = c.z(0).unwrap();
    assert_eq!(
        c.measure(z, Some(1)),
        Err(PbcError::ClassicalBitOutOfRange(1))
    );
    let a = c.measure(z, Some(0)).unwrap();
    let b = c.measure(z, Some(0)).unwrap();
    let hidden = c.measure(z, None).unwrap();
    assert_eq!((a.index(), b.index(), hidden.index()), (0, 1, 2));
    c.conditional_pauli(z, a).unwrap();
    assert_eq!(c.measurement_count(), 3);
    assert_eq!(c.operations().len(), 4);
}

#[test]
fn conditional_rotations_keep_axis_angle_and_original_outcome() {
    let mut c = PbcCircuit::new(2, 1);
    let x = c.x(0).unwrap();
    let z = c.z(1).unwrap();
    let product = c.product(x.as_ref(), z.as_ref()).unwrap();
    let axis = c.hermitian_axis(product, 100).unwrap().negated();
    let m = c.measure(z, Some(0)).unwrap();
    c.measure(z, Some(0)).unwrap();
    for k in 0..8 {
        let angle = PauliAngle::new(k);
        c.conditional_rotate(axis, angle, m).unwrap();
        let op = c.operations().last().unwrap();
        assert_eq!(
            *op,
            PbcOp::ConditionalRotate {
                axis,
                angle,
                if_one: m
            }
        );
        assert_eq!(op.axis(), axis);
        assert!(
            c.to_ascii()
                .unwrap()
                .contains(&format!("R({angle},-) if m0=1"))
        );
    }
    c.conditional_pauli(axis, m).unwrap();
    assert_eq!(
        *c.operations().last().unwrap(),
        PbcOp::ConditionalRotate {
            axis,
            angle: PauliAngle::new(4),
            if_one: m,
        }
    );
}

#[test]
fn conditional_rotations_reject_invalid_handles_without_appending() {
    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    let m = c.measure(x, None).unwrap();
    let mut other = PbcCircuit::new(1, 0);
    let y = other.y(0).unwrap();
    let foreign_m = other.measure(y, None).unwrap();
    let angle = PauliAngle::new(2);
    assert_eq!(
        c.conditional_rotate(y, angle, m),
        Err(PbcError::ForeignPauli)
    );
    assert_eq!(
        c.conditional_rotate(x, angle, foreign_m),
        Err(PbcError::UnknownMeasurement)
    );
    let future_m = MeasId {
        owner: m.owner,
        index: c.measurement_count(),
    };
    assert_eq!(
        c.conditional_rotate(x, angle, future_m),
        Err(PbcError::UnknownMeasurement)
    );
    assert_eq!(c.operations().len(), 1);
}

#[test]
fn output_frame_gates_and_operands_are_validated() {
    let mut c = PbcCircuit::new(2, 0);
    assert_eq!(c.x(2), Err(PbcError::QubitOutOfRange(2)));
    assert_eq!(
        c.push_output_clifford(Gate::t(0)),
        Err(PbcError::NonCliffordFrameGate)
    );
    assert_eq!(
        c.push_output_clifford(Gate::h(2)),
        Err(PbcError::QubitOutOfRange(2))
    );
    assert_eq!(
        c.push_output_clifford(Gate::cz {
            control: 0,
            target: 0
        }),
        Err(PbcError::RepeatedOperand)
    );
    for gate in [
        Gate::h(0),
        Gate::x(0),
        Gate::z(1),
        Gate::s(0),
        Gate::sdg(1),
        Gate::cnot {
            control: 0,
            target: 1,
        },
        Gate::cz {
            control: 1,
            target: 0,
        },
    ] {
        c.push_output_clifford(gate).unwrap();
    }
    assert_eq!(c.output_frame().len(), 2);
    assert!(c.to_ascii().unwrap().contains("X0->"));
}

#[test]
fn shared_deep_dag_expands_and_drops_without_recursion() {
    let mut c = PbcCircuit::new(1, 0);
    let x = c.x(0).unwrap();
    let mut p = x.as_ref();
    for _ in 0..20_000 {
        p = c.product(p, p).unwrap();
    }
    assert_eq!(c.pauli_nodes().len(), 20_004);
    assert_eq!(c.expand(p, 20_004), Err(PbcError::ExpansionLimit));
    assert_eq!(c.expand(p, 20_005).unwrap().factors, vec![Pauli::I]);
    assert_eq!(c.expand(x.as_ref(), 4).unwrap().factors, vec![Pauli::X]);
}

#[test]
fn ascii_shows_joint_axes_outcomes_and_frame() {
    let mut c = PbcCircuit::new(3, 1);
    let x = c.x(0).unwrap();
    let z = c.z(2).unwrap();
    let p = c.product(x.as_ref(), z.as_ref()).unwrap();
    let axis = c.hermitian_axis(p, 100).unwrap();
    c.rotate(axis, PauliAngle::new(1)).unwrap();
    let m = c.measure(z.negated(), Some(0)).unwrap();
    c.conditional_pauli(x, m).unwrap();
    c.push_output_clifford(Gate::h(1)).unwrap();
    let drawing = c.to_ascii().unwrap();
    assert!(drawing.is_ascii());
    assert!(drawing.contains("R(pi/8,+)"));
    assert!(drawing.contains("M(-)->m0/c0"));
    assert!(drawing.contains("R(pi/2,+) if m0=1"));
    assert!(drawing.contains("X1->"));
    let lines: Vec<_> = drawing.lines().collect();
    let x_column = lines[1].find('X').unwrap();
    assert_eq!(lines[2].as_bytes()[x_column], b'|');
    assert_eq!(lines[3].as_bytes()[x_column], b'|');
    assert_eq!(lines[5].as_bytes()[x_column], b'Z');
}

#[test]
fn ascii_limits_and_identity_are_explicit() {
    let mut c = PbcCircuit::new(0, 0);
    assert_eq!(c.to_ascii().unwrap(), "(empty PBC circuit; 0 qubits)\n");
    c.measure(c.identity().negated(), None).unwrap();
    assert!(c.to_ascii().unwrap().contains("M(-I)->m0"));
    assert_eq!(
        c.to_ascii_with(AsciiOptions {
            max_operations: 0,
            ..AsciiOptions::default()
        }),
        Err(PbcError::DrawingLimit)
    );
    assert_eq!(
        c.to_ascii_with(AsciiOptions {
            max_expansion_cells: 0,
            ..AsciiOptions::default()
        }),
        Err(PbcError::ExpansionLimit)
    );
    assert_eq!(
        PbcCircuit::new(33, 0).to_ascii(),
        Err(PbcError::DrawingLimit)
    );
    assert!(std::panic::catch_unwind(|| PbcCircuit::new(usize::MAX, 0)).is_err());
}

#[test]
fn ascii_cx_with_t_on_each_wire_aligns_all_connectors() {
    let input = crate::circuit::Circuit {
        num_qubits: 2,
        num_cbits: 0,
        gates: vec![
            Gate::cnot {
                control: 0,
                target: 1,
            },
            Gate::t(0),
            Gate::t(1),
        ],
    };
    let drawing = to_pbc(&input).unwrap().to_ascii().unwrap();
    let rows: Vec<_> = drawing.lines().collect();
    let joint_axis = rows[1].rfind('Z').unwrap();
    let frame_x = rows[1].find('X').unwrap();
    let frame_z = rows[1].rfind('Z').unwrap();
    for column in [joint_axis, frame_x, frame_z] {
        assert_eq!(rows[2].as_bytes()[column], b'|', "{drawing}");
    }
    assert_eq!(rows[3].as_bytes()[joint_axis], b'Z');
    assert_eq!(rows[3].as_bytes()[frame_x], b'X');
    assert_eq!(rows[3].as_bytes()[frame_z], b'Z');
    assert_eq!(rows[1].len(), rows[3].len());
    assert_eq!(
        drawing,
        concat!(
            "     R(pi/8,+)  R(pi/8,+)  X0->+  Z1->+\n",
            "q0: ----[Z]--------[Z]------[X]----[Z]--\n",
            "                    |        |      |\n",
            "q1: ---------------[Z]------[X]----[Z]--\n",
        )
    );
}

#[test]
fn ascii_alignment_survives_heading_widths_and_two_digit_qubits() {
    let mut c = PbcCircuit::new(12, 1);
    let x = c.x(0).unwrap();
    let z = c.z(11).unwrap();
    let product = c.product(x.as_ref(), z.as_ref()).unwrap();
    let joint = c.hermitian_axis(product, 100).unwrap();
    for angle in 0..8 {
        c.rotate(joint, PauliAngle::new(angle)).unwrap();
    }
    let outcome = c.measure(joint.negated(), Some(0)).unwrap();
    c.conditional_rotate(joint, PauliAngle::new(-1), outcome)
        .unwrap();
    c.push_output_clifford(Gate::cnot {
        control: 11,
        target: 0,
    })
    .unwrap();
    c.push_output_clifford(Gate::cz {
        control: 0,
        target: 11,
    })
    .unwrap();
    let drawing = c.to_ascii().unwrap();
    let rows: Vec<_> = drawing.lines().collect();
    for (column, mark) in rows[1].bytes().enumerate() {
        if !matches!(mark, b'X' | b'|' | b'@') {
            continue;
        }
        for row in &rows[2..23] {
            assert_eq!(row.as_bytes()[column], b'|', "{drawing}");
        }
        assert!(
            matches!(rows[23].as_bytes()[column], b'Z' | b'|' | b'@'),
            "{drawing}"
        );
    }
}

#[test]
fn rotation_weights_count_factors_and_skip_measurements() {
    let mut c = PbcCircuit::new(3, 1);
    let x = c.x(0).unwrap();
    let z = c.z(2).unwrap();
    let product = c.product(x.as_ref(), z.as_ref()).unwrap();
    let axis = c.hermitian_axis(product, 100).unwrap();
    c.rotate(axis, PauliAngle::new(1)).unwrap();
    c.measure(z, Some(0)).unwrap();
    c.rotate(c.identity(), PauliAngle::new(1)).unwrap();
    let m = c.measure(x, None).unwrap();
    c.conditional_pauli(z, m).unwrap();
    assert_eq!(c.rotation_weights(1000).unwrap(), vec![2, 0, 1]);
}
