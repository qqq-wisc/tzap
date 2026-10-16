use tzap::{
    angle::Angle,
    circuit::{Circuit, Gate},
    qasm,
};

fn gates(angle: &Angle) -> Vec<Gate> {
    vec![
        Gate::p(angle.clone(), 1),
        Gate::rx(angle.clone(), 1),
        Gate::ry(angle.clone(), 1),
        Gate::rz(angle.clone(), 1),
        Gate::cp {
            lambda: angle.clone(),
            control: 0,
            target: 1,
        },
        Gate::crx {
            theta: angle.clone(),
            control: 0,
            target: 1,
        },
        Gate::cry {
            theta: angle.clone(),
            control: 0,
            target: 1,
        },
        Gate::crz {
            theta: angle.clone(),
            control: 0,
            target: 1,
        },
    ]
}

#[test]
fn every_angle_gate_is_numeric_with_accurate_conversion_counts() {
    for source in ["pi/7", "(pi/7)+(1e17)", "-(pi/7)"] {
        let angle = qasm::parse(&format!("qreg q[1]; rz({source}) q[0];"))
            .unwrap()
            .gates()[0]
            .angle()
            .unwrap()
            .clone();
        for gate in gates(&angle) {
            let mut input = Circuit::with_cbits(2, 3);
            input = Circuit::from_gates(
                input.num_qubits,
                input.num_cbits,
                vec![gate.clone(), Gate::measure { qubit: 1, cbit: 2 }],
            )
            .unwrap();
            let output = qasm::serialize_numeric_lossy(&input).unwrap();
            assert!(!output.qasm.contains("pi"), "{gate:?}");
            assert_eq!(output.uncertified_conversions, 1, "{gate:?}");
            let parsed = qasm::parse_numeric_sdk(&output.qasm).unwrap();
            assert_eq!((parsed.num_qubits, parsed.num_cbits), (2, 3));
            assert_eq!(parsed.gates().last(), input.gates().last());
            for g in parsed.gates() {
                if let Some(a) = g.angle() {
                    assert_eq!(a.components().unwrap().0.numerator(), 0);
                }
            }
            assert!(qasm::serialize_sdk(&input).unwrap_err().contains("gate 0"));
        }
    }
}

#[test]
fn named_expansions_are_not_uncertified_conversions() {
    for gate in gates(&Angle::pi_fraction(2, 1).unwrap()) {
        let input = Circuit::from_gates(2, 0, vec![gate]).unwrap();
        let export = qasm::serialize_numeric_lossy(&input).unwrap();
        assert_eq!(export.uncertified_conversions, 0);
        assert!(!export.qasm.contains("pi"));
        let exact = qasm::serialize_sdk(&input).unwrap();
        assert!(!exact.contains("pi"));
    }
}

#[test]
fn floats_remain_literal_and_preserved_expressions_are_rejected() {
    for value in [0.0, f64::from_bits(1), -f64::from_bits(1), 1e17, f64::MAX] {
        let input = Circuit::from_gates(2, 0, gates(&Angle::from_f64(value).unwrap())).unwrap();
        let output = qasm::serialize_numeric_lossy(&input).unwrap();
        assert_eq!(output.uncertified_conversions, 0);
        let parsed = qasm::parse_numeric_sdk(&output.qasm).unwrap();
        for gate in parsed.gates() {
            if let Some(a) = gate.angle() {
                assert_eq!(a.components().unwrap().1.get().to_bits(), value.to_bits());
            }
        }
        assert_eq!(
            qasm::parse_numeric_sdk(&qasm::serialize_sdk(&input).unwrap())
                .unwrap()
                .gates(),
            qasm::parse(&input.to_qasm()).unwrap().gates()
        );
    }
    let preserved = qasm::parse("qreg q[1]; rz(pi/9223372036854775808) q[0];")
        .unwrap()
        .gates()[0]
        .angle()
        .unwrap()
        .clone();
    for gate in gates(&preserved) {
        let input = Circuit::from_gates(2, 0, vec![Gate::h(0), gate]).unwrap();
        let error = qasm::serialize_sdk(&input).unwrap_err();
        assert!(error.contains("gate 1"));
        assert_eq!(
            qasm::serialize_numeric_lossy(&input)
                .unwrap()
                .uncertified_conversions,
            1
        );
    }
}

#[test]
fn sdk_export_splits_large_mixed_angles_without_rounding() {
    for numeric in [1e17, -1e17, f64::MAX] {
        let angle = Angle::pi_fraction(1, 4)
            .unwrap()
            .checked_add_signed(1, &Angle::from_f64(numeric).unwrap())
            .unwrap();
        for gate in gates(&angle).into_iter().take(4) {
            let input = Circuit::from_gates(2, 0, vec![gate]).unwrap();
            let output = qasm::serialize_sdk(&input).unwrap();
            assert!(output.contains("t q[1]"));
            assert!(!output.contains("pi"));
            let parsed = qasm::parse_numeric_sdk(&output).unwrap();
            let angles: Vec<_> = parsed.gates().iter().filter_map(Gate::angle).collect();
            assert_eq!(angles.len(), 1);
            assert_eq!(
                angles[0].components().unwrap().1.get().to_bits(),
                numeric.to_bits()
            );
            assert_eq!(
                qasm::serialize_numeric_lossy(&input)
                    .unwrap()
                    .uncertified_conversions,
                0
            );
        }
    }
}

#[test]
fn sdk_export_splits_nonrepresentable_dyadic_sums_exactly() {
    use num_rational::BigRational;
    let mut ctx = tzap::angle::AngleContext::default();
    let left = Angle::from_f64(0.1).unwrap();
    let right = Angle::from_f64(0.2).unwrap();
    let sum = ctx.add(&left, &right).unwrap();
    assert!(sum.to_f64_exact().is_none());
    let expected = left.as_rational().unwrap().to_big_rational()
        + right.as_rational().unwrap().to_big_rational();
    for gate in gates(&sum) {
        let input = Circuit::from_gates(2, 0, vec![gate]).unwrap();
        let text = qasm::serialize_sdk(&input).unwrap();
        let output = qasm::parse_numeric_sdk(&text).unwrap();
        let total: BigRational = output
            .gates()
            .iter()
            .filter_map(Gate::angle)
            .map(|a| a.as_rational().unwrap().to_big_rational())
            .sum();
        assert_eq!(total, expected);
        assert_eq!(
            qasm::serialize_numeric_lossy(&input)
                .unwrap()
                .uncertified_conversions,
            0
        );
    }
    let source = qasm::parse("qreg q[1]; rz(0.1+0.2) q[0];").unwrap();
    assert!(qasm::serialize_sdk(&source).is_err());
    let tiny = qasm::parse("qreg q[1]; rz(1/2^1075) q[0];").unwrap();
    assert!(qasm::serialize_sdk(&tiny).is_err());
}
