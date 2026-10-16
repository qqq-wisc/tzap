use crate::{
    angle::Angle,
    circuit::{Circuit, Gate},
    decompose::DecomposeRotations,
    optimize::{Options, PassName, optimize},
};

#[test]
fn controlled_angles_preserve_literal_turns_in_equality_and_lowering() {
    for name in ["crx", "cry", "crz"] {
        let two = Circuit::from_qasm(&format!("qreg q[2]; {name}(2*pi) q[0],q[1];")).unwrap();
        let zero = Circuit::from_qasm(&format!("qreg q[2]; {name}(0) q[0],q[1];")).unwrap();
        assert_ne!(two.gates, zero.gates);
        let output = DecomposeRotations::default().try_run(&two).unwrap();
        // A 2π controlled rotation is Z on its control, not identity.
        let expected = Circuit::from_qasm("qreg q[2]; z q[0];").unwrap();
        assert!(crate::unitary::circuits_equiv(&output, &expected, 1e-12));
        let four = Circuit::from_qasm(&format!("qreg q[2]; {name}(4*pi) q[0],q[1];")).unwrap();
        let output = DecomposeRotations::default().try_run(&four).unwrap();
        assert!(crate::unitary::circuits_equiv(
            &output,
            &Circuit::new(2),
            1e-12
        ));
        let reparsed = Circuit::from_qasm(&two.to_qasm()).unwrap();
        assert!(crate::unitary::circuits_equiv(&two, &reparsed, 1e-12));
    }
}

#[test]
fn every_parameterized_gate_preserves_pi_fractions_and_opaque_expressions() {
    for (name, operands) in [
        ("p", "q[0]"),
        ("rx", "q[0]"),
        ("ry", "q[0]"),
        ("rz", "q[0]"),
        ("cp", "q[0],q[1]"),
        ("crx", "q[0],q[1]"),
        ("cry", "q[0],q[1]"),
        ("crz", "q[0],q[1]"),
    ] {
        let input = Circuit::from_qasm(&format!("qreg q[2]; {name}(pi/7) {operands};")).unwrap();
        let angle = input.gates[0].angle().unwrap();
        let (pi, radians) = angle.components().unwrap();
        assert_eq!(
            (pi.numerator(), pi.denominator(), radians.get()),
            (1, 7, 0.0)
        );
        let serialized = input.to_qasm();
        assert!(serialized.contains("pi"));
        let parsed = Circuit::from_qasm(&serialized).unwrap();
        assert!(crate::unitary::circuits_equiv(&input, &parsed, 1e-12));

        let input = Circuit::from_qasm(&format!(
            "qreg q[2]; {name}(pi/9223372036854775808) {operands};"
        ))
        .unwrap();
        let options = Options {
            passes: Some(vec![PassName::PhaseFoldPauli]),
            ..Options::default()
        };
        let (output, report) = optimize(&input, &options).unwrap();
        assert_eq!(input.gates, output.gates);
        assert_eq!(report.numerical.preserved_expressions, 0);
        assert!(crate::unitary::circuits_equiv(
            &output,
            &Circuit::from_qasm(&output.to_qasm()).unwrap(),
            1e-12
        ));
        assert!(DecomposeRotations::default().try_run(&input).is_ok());
    }
}

#[test]
fn controlled_synthesis_certifies_unrepresentable_half_angles() {
    let angle = Angle::from_f64(f64::from_bits(1)).unwrap();
    let mut input = Circuit::new(2);
    input.apply(Gate::crx {
        theta: angle,
        control: 0,
        target: 1,
    });
    assert!(
        DecomposeRotations::default()
            .try_run(&input)
            .unwrap()
            .gates
            .is_empty()
    );
    // Export uses a CRZ basis change and never halves or drops the subnormal.
    let parsed = Circuit::from_qasm(&input.to_qasm()).unwrap();
    assert_eq!(
        parsed
            .gates
            .iter()
            .find_map(Gate::angle)
            .unwrap()
            .to_f64_lossy()
            .unwrap(),
        f64::from_bits(1)
    );
}

#[test]
fn controlled_synthesis_reports_epsilon_underflow_before_gridsynth() {
    let mut input = Circuit::new(2);
    input.apply(Gate::crz {
        theta: Angle::from_f64(0.37).unwrap(),
        control: 0,
        target: 1,
    });
    let error = DecomposeRotations {
        epsilon: f64::from_bits(1),
    }
    .try_run(&input)
    .unwrap_err();
    assert!(error.to_string().contains("epsilon underflows"));
}

/// A gate that repeats a qubit operand has no unitary meaning, so no circuit
/// can hold one, however it is built.
#[test]
fn circuits_reject_gates_with_repeated_operands() {
    use crate::circuit::RepeatedOperand;
    let degenerate = [
        Gate::cnot {
            control: 0,
            target: 0,
        },
        Gate::cz {
            control: 1,
            target: 1,
        },
        Gate::swap(2, 2),
        Gate::ccx {
            control1: 0,
            control2: 0,
            target: 1,
        },
        Gate::ccz {
            control1: 0,
            control2: 1,
            target: 1,
        },
        Gate::cswap {
            control: 2,
            first: 0,
            second: 2,
        },
    ];
    for gate in degenerate {
        let mut circuit = Circuit::new(3);
        assert!(matches!(
            circuit.try_apply(gate.clone()),
            Err(RepeatedOperand(g)) if g == gate
        ));
        assert!(circuit.gates().is_empty());
        assert!(Circuit::from_gates(3, 0, vec![Gate::h(0), gate.clone()]).is_err());
        let panic = std::panic::catch_unwind(move || Circuit::new(3).apply(gate));
        assert!(panic.is_err());
    }
    let valid = vec![
        Gate::cnot {
            control: 0,
            target: 1,
        },
        Gate::ccz {
            control1: 0,
            control2: 1,
            target: 2,
        },
    ];
    let circuit = Circuit::from_gates(3, 0, valid.clone()).unwrap();
    assert_eq!(circuit.gates(), valid);
}
