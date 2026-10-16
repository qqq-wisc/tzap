//! Concrete numerical failures from the investigation, with independent exact
//! rational phase semantics. These tests do not compare tiny errors to a tolerance.
#[path = "support/monomial.rs"]
mod monomial;
use monomial::{binary, equivalent, semantics};
use num_rational::BigRational;
use tzap::{
    angle::Angle,
    cancel::CancelGates,
    circuit::{Circuit, Gate},
    cnot_min::CnotMin,
    decompose::DecomposeRotations,
    optimize::{Level, Options, PassName, optimize},
    pass::Pass,
    phase_fold_pauli::PhaseFoldPauli,
    phase_fold_rand::PhaseFoldRand,
    qasm,
};

fn check_every_exact_pipeline(input: &Circuit) {
    let expected = semantics(input);
    let cnot = CnotMin::default();
    for pass in [
        &CancelGates as &dyn Pass,
        &PhaseFoldRand,
        &PhaseFoldPauli,
        &cnot,
    ] {
        assert!(equivalent(&expected, &pass.run(input)), "{}", pass.name());
    }
    for level in [Level::O1, Level::O2, Level::O3] {
        for parallel in [false, true] {
            let (output, _) = optimize(
                input,
                &Options {
                    level,
                    parallel,
                    ..Options::default()
                },
            )
            .unwrap();
            assert!(
                equivalent(&expected, &output),
                "{level:?}, parallel={parallel}"
            );
            let reparsed = qasm::parse(&output.to_qasm()).unwrap();
            assert!(equivalent(&expected, &reparsed), "round trip: {level:?}");
        }
    }
}

#[test]
fn all_short_diagonal_failures_across_passes_levels_and_parallel_execution() {
    for body in [
        "rz(1e16) q[0]; rz(1) q[0]; rz(-1e16) q[0];",
        "rz(1e17) q[0]; rz(pi) q[0]; rz(-1e17) q[0];",
        "t q[0]; rz(1e17) q[0]; rz(-1e17) q[0];",
        "rz(5e-10) q[0];",
        "rz(5e-7) q[0];",
        "rz(1e16) q[0];",
        "rz(28296951008113760) q[0];",
        "rz(1e308) q[0]; rz(1e308) q[0];",
        "rz(0.1) q[0]; rz(0.2) q[0];",
        "rz(0.125) q[0]; rz(0.25) q[0];",
        "rz((1e17+pi)-1e17) q[0];",
        "rz(10000000000000000*pi) q[0]; rz(pi) q[0]; rz(-10000000000000000*pi) q[0];",
        "rz(2147483647*pi/4294967295) q[0]; rz(2147483648*pi/4294967297) q[0];",
        "rz(0.125) q[0]; rz(0.25) q[0]; rz(1e17) q[0]; rz(-1e17) q[0];",
        "rz(1e16) q[0]; x q[0]; rz(1) q[0]; x q[0]; rz(-1e16) q[0];",
    ] {
        let input = qasm::parse(&format!("qreg q[1]; {body}")).unwrap();
        check_every_exact_pipeline(&input);
    }
    // Numeric pi has no symbolic provenance, even one rounding step away from S.
    let numeric = std::f64::consts::FRAC_PI_2 + 5e-10;
    let mut input = Circuit::new(2);
    input = Circuit::from_gates(
        input.num_qubits,
        input.num_cbits,
        vec![
            Gate::cnot {
                control: 0,
                target: 1
            };
            2
        ],
    )
    .unwrap();
    input.apply(Gate::rz_f64(numeric, 0).unwrap());
    check_every_exact_pipeline(&input);
}

#[test]
fn large_numeric_sums_are_preserved_for_each_single_qubit_axis() {
    for kind in ["rx", "ry", "rz", "p"] {
        for values in [
            vec![1e16, 1.0, -1e16],
            vec![0.1, 0.2],
            vec![1e308, 1e308],
            vec![5e-10],
            vec![0.125, 0.25],
        ] {
            let mut input = Circuit::new(1);
            input = Circuit::from_gates(
                input.num_qubits,
                input.num_cbits,
                values
                    .iter()
                    .map(|v| {
                        let a = Angle::from_f64(*v).unwrap();
                        match kind {
                            "rx" => Gate::rx(a, 0),
                            "ry" => Gate::ry(a, 0),
                            "rz" => Gate::rz(a, 0),
                            "p" => Gate::p(a, 0),
                            _ => unreachable!(),
                        }
                    })
                    .collect(),
            )
            .unwrap();
            let expected: BigRational = values.iter().map(|v| binary(*v)).sum();
            let sum = |output: &Circuit| -> BigRational {
                output
                    .gates()
                    .iter()
                    .map(|g| {
                        let angle = match (kind, g) {
                            ("rx", Gate::rx(a, 0)) | ("ry", Gate::ry(a, 0)) => a,
                            ("rz" | "p", Gate::rz(a, 0) | Gate::p(a, 0)) => a,
                            _ => panic!("unexpected gate for {kind}: {g:?}"),
                        };
                        angle
                            .as_rational()
                            .expect("numeric provenance")
                            .to_big_rational()
                    })
                    .sum()
            };
            for level in [Level::O1, Level::O2, Level::O3] {
                let (output, _) = optimize(
                    &input,
                    &Options {
                        level,
                        ..Options::default()
                    },
                )
                .unwrap();
                assert_eq!(sum(&output), expected, "{kind} {values:?} {level:?}");
            }
        }
    }
}

#[test]
fn decimal_chain_preserves_exact_sum_of_ten_thousand_accepted_values() {
    let mut input = Circuit::new(1);
    input = Circuit::from_gates(
        input.num_qubits,
        input.num_cbits,
        vec![Gate::rz_f64(0.0001, 0).unwrap(); 10_000],
    )
    .unwrap();
    check_every_exact_pipeline(&input);
}

fn repeated_signed_rotations(count: usize) -> Circuit {
    let mut input = Circuit::new(1);
    for _ in 0..count {
        input.extend([
            Gate::rz_f64(9e-7, 0).unwrap(),
            Gate::y(0),
            Gate::rz_f64(-9e-7, 0).unwrap(),
            Gate::y(0),
        ]);
    }
    input
}
fn repeated_phase_rotations(count: usize) -> Circuit {
    let mut input = Circuit::new(1);
    for _ in 0..count {
        input.extend([
            Gate::p_f64(3e-6, 0).unwrap(),
            Gate::p_f64(3e-6, 0).unwrap(),
            Gate::x(0),
            Gate::p_f64(-3e-6, 0).unwrap(),
            Gate::p_f64(-3e-6, 0).unwrap(),
            Gate::x(0),
        ]);
    }
    input
}
#[test]
fn repeated_tiny_rotations_survive_y_axis_sign_changes() {
    check_every_exact_pipeline(&repeated_signed_rotations(1_000));
    check_every_exact_pipeline(&repeated_phase_rotations(1_000));
}
#[test]
#[ignore = "400,000- and 600,000-gate regressions; run in release"]
fn full_hundred_thousand_repetition_regression() {
    check_every_exact_pipeline(&repeated_signed_rotations(100_000));
    check_every_exact_pipeline(&repeated_phase_rotations(100_000));
}

#[test]
fn reversed_two_qubit_operand_order_preserves_the_phase_wire() {
    // Inlined form of the BQSKit g(a) reproducer; custom gate declarations
    // remain outside tzap's parser subset.
    let input = qasm::parse("qreg q[2]; rz(0) q[1]; rz(pi) q[0];").unwrap();
    check_every_exact_pipeline(&input);
}

#[test]
fn cnot_chunk_rolls_back_after_partial_cz_phase_updates() {
    let input = qasm::parse(&format!(
        "qreg q[2]; cx q[0],q[1]; cx q[0],q[1]; rz(pi/{}) q[1]; cz q[0],q[1];",
        i64::MAX
    ))
    .unwrap();
    let output = CnotMin::default().run(&input);
    assert!(output.gates().len() <= input.gates().len());
    assert!(equivalent(&semantics(&input), &output));
}

#[test]
fn controlled_phases_use_literal_angles_before_decomposition() {
    for body in [
        "crz(2*pi) q[0],q[1];",
        "crz(4*pi) q[0],q[1];",
        "cp(2*pi) q[0],q[1];",
    ] {
        let input = qasm::parse(&format!("qreg q[2]; {body}")).unwrap();
        let output = DecomposeRotations { epsilon: 1e-20 }
            .try_run(&input)
            .unwrap();
        assert!(equivalent(&semantics(&input), &output), "{body}");
    }
    let mut tiny = Circuit::new(2);
    tiny.apply(Gate::crz {
        theta: Angle::from_f64(f64::from_bits(1)).unwrap(),
        control: 0,
        target: 1,
    });
    let output = DecomposeRotations { epsilon: 1e-20 }
        .try_run(&tiny)
        .unwrap();
    assert!(
        output.gates().is_empty(),
        "tiny gate is removed only by a certified approximate pass"
    );
    check_every_exact_pipeline(&tiny);
    // A later unsupported expression must discard an earlier successful synthesis.
    let unsupported = qasm::parse("qreg q[2]; crz(exp(2000)) q[0],q[1];").unwrap();
    tiny = Circuit::from_gates(
        tiny.num_qubits,
        tiny.num_cbits,
        vec![
            Gate::rz_f64(0.3, 0).unwrap(),
            unsupported.gates()[0].clone(),
        ],
    )
    .unwrap();
    let pass = DecomposeRotations { epsilon: 1e-20 };
    assert!(pass.try_run(&tiny).unwrap_err().0.contains("gate 1"));
    assert_eq!(pass.run(&tiny).gates(), tiny.gates());
}

#[test]
fn large_angle_synthesis_is_certified_at_tight_epsilon() {
    let input = qasm::parse("qreg q[1]; rz(1e16) q[0];").unwrap();
    let (output, report) = optimize(
        &input,
        &Options {
            passes: Some(vec![PassName::DecomposeRotations]),
            rz_epsilon: 1e-20,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(
        output
            .gates()
            .iter()
            .all(|g| matches!(g, Gate::h(_) | Gate::t(_) | Gate::s(_) | Gate::x(_)))
    );
    assert_eq!(report.numerical.certified_syntheses, 1);
    assert_eq!(report.numerical.uncertified_syntheses, 0);
    assert_eq!(report.numerical.failed_certifications, 0);
    assert!(report.numerical.certified_error_bound.unwrap() >= 1e-20);
}

#[test]
fn parser_rounding_is_separate_from_optimizer_rounding() {
    let input = qasm::parse("qreg q[1]; rz((1e17+1)-1e17) q[0];").unwrap();
    assert_eq!(input.gates()[0].angle().unwrap(), &Angle::ONE);
    let preserved_pi = qasm::parse("qreg q[1]; rz((1e17+pi)-1e17) q[0];").unwrap();
    assert_eq!(
        preserved_pi.gates()[0]
            .angle()
            .unwrap()
            .components()
            .unwrap()
            .0
            .numerator(),
        1
    );
}

#[test]
fn exactly_lowered_controlled_gates_do_not_need_a_numeric_budget() {
    let input = qasm::parse("qreg q[2]; crz(2*pi) q[0],q[1];").unwrap();
    let output = DecomposeRotations {
        epsilon: f64::from_bits(1),
    }
    .try_run(&input)
    .unwrap();
    assert!(equivalent(&semantics(&input), &output));
}

#[test]
fn invalid_operands_and_precision_limits_are_transactional_errors() {
    let valid = qasm::parse("qreg q[1]; rz(0.3) q[0];").unwrap();
    let mut invalid = valid.clone();
    invalid.apply(Gate::rz_f64(0.2, 1).unwrap());
    let pass = DecomposeRotations { epsilon: 1e-20 };
    assert!(pass.try_run(&invalid).unwrap_err().0.contains("gate 1"));
    assert_eq!(pass.run(&invalid).gates(), invalid.gates());
    assert!(
        DecomposeRotations { epsilon: 1e-100 }
            .try_run(&valid)
            .unwrap_err()
            .0
            .contains("precision limit")
    );
    // An aliased controlled rotation cannot reach the pass at all.
    let aliased = Gate::crz {
        theta: Angle::pi_fraction(2, 1).unwrap(),
        control: 0,
        target: 0,
    };
    assert!(Circuit::new(1).try_apply(aliased).is_err());
}

#[test]
fn every_native_rotation_has_a_certified_tight_epsilon_path() {
    for gate in ["rx", "ry", "rz", "p", "crx", "cry", "crz", "cp"] {
        let operands = if gate.starts_with('c') {
            "q[1],q[0]"
        } else {
            "q[0]"
        };
        let input = qasm::parse(&format!("qreg q[2]; {gate}(pi/7+0.3) {operands};")).unwrap();
        for epsilon in [1e-10, 1e-20, 1e-40] {
            let (output, report) = optimize(
                &input,
                &Options {
                    passes: Some(vec![PassName::DecomposeRotations]),
                    rz_epsilon: epsilon,
                    parallel: true,
                    ..Options::default()
                },
            )
            .unwrap();
            assert!(!output.gates().is_empty());
            assert_eq!(report.numerical.certified_syntheses, 1, "{gate}");
            assert_eq!(report.numerical.uncertified_syntheses, 0);
            assert_eq!(report.numerical.failed_certifications, 0);
            assert!(report.numerical.certified_error_bound.unwrap() >= epsilon);
        }
    }
}

#[test]
fn synthesis_bounds_accumulate_outward_and_remain_per_run() {
    let mut input = Circuit::new(1);
    input = Circuit::from_gates(
        input.num_qubits,
        input.num_cbits,
        vec![Gate::rz_f64(0.3, 0).unwrap(); 100],
    )
    .unwrap();
    let epsilon = 1e-20;
    for parallel in [false, true] {
        let options = Options {
            passes: Some(vec![PassName::DecomposeRotations]),
            rz_epsilon: epsilon,
            parallel,
            ..Options::default()
        };
        let (_, report) = optimize(&input, &options).unwrap();
        assert_eq!(report.numerical.certified_syntheses, 100);
        assert!(
            binary(report.numerical.certified_error_bound.unwrap())
                >= binary(epsilon) * BigRational::from_integer(100.into())
        );
        let named = qasm::parse("qreg q[1]; t q[0];").unwrap();
        let (_, exact_report) = optimize(&named, &options).unwrap();
        assert_eq!(exact_report.numerical.certified_syntheses, 0);
        assert_eq!(exact_report.numerical.certified_error_bound, Some(0.0));
    }
    let (_, overflow) = optimize(
        &input,
        &Options {
            passes: Some(vec![PassName::DecomposeRotations]),
            rz_epsilon: f64::MAX,
            ..Options::default()
        },
    )
    .unwrap();
    assert_eq!(overflow.numerical.certified_error_bound, None);
}

#[test]
fn sdk_export_preserves_exact_rational_phases_in_all_diagonalized_axes() {
    for value in [1e17, -1e17, f64::MAX, f64::from_bits(1)] {
        for axis in ["p", "rz", "rx", "ry"] {
            let mut input = Circuit::new(1);
            input.apply(Gate::t(0));
            if matches!(axis, "rx" | "ry") {
                input.apply(Gate::h(0));
            }
            if axis == "ry" {
                input.apply(Gate::s(0));
            }
            let angle = Angle::from_f64(value).unwrap();
            input.apply(match axis {
                "p" => Gate::p(angle, 0),
                "rx" => Gate::rx(angle, 0),
                "ry" => Gate::ry(angle, 0),
                _ => Gate::rz(angle, 0),
            });
            if axis == "ry" {
                input.apply(Gate::sdg(0));
            }
            if matches!(axis, "rx" | "ry") {
                input.apply(Gate::h(0));
            }
            let expected = semantics(&input);
            for level in [Level::O1, Level::O2, Level::O3] {
                for parallel in [false, true] {
                    let (optimized, _) = optimize(
                        &input,
                        &Options {
                            level,
                            parallel,
                            ..Options::default()
                        },
                    )
                    .unwrap();
                    let exported =
                        qasm::parse_numeric_sdk(&qasm::serialize_sdk(&optimized).unwrap()).unwrap();
                    assert!(
                        equivalent(&expected, &exported),
                        "{axis}, {value}, {level:?}, parallel={parallel}"
                    );
                }
            }
        }
    }
}

#[test]
fn certified_tiny_controlled_identity_is_counted_and_is_never_an_exact_fold() {
    for name in ["cp", "crx", "cry", "crz"] {
        for source in ["5e-324", "-5e-324", "1.5e-323", "-1.5e-323"] {
            let input = qasm::parse(&format!("qreg q[2]; {name}({source}) q[1],q[0];")).unwrap();
            for level in [Level::O1, Level::O2, Level::O3] {
                let (output, report) = optimize(
                    &input,
                    &Options {
                        level,
                        ..Options::default()
                    },
                )
                .unwrap();
                assert_eq!(output.gates(), input.gates());
                assert_eq!(report.numerical.certified_syntheses, 0);
            }
            let (output, report) = optimize(
                &input,
                &Options {
                    passes: Some(vec![PassName::DecomposeRotations]),
                    rz_epsilon: 1e-20,
                    ..Options::default()
                },
            )
            .unwrap();
            assert!(output.gates().is_empty());
            assert_eq!(report.numerical.certified_syntheses, 1);
            assert!(report.numerical.certified_error_bound.unwrap() >= 1e-20);
        }
    }
    // A named angle approximated by identity still consumes a synthesis bound.
    let input = qasm::parse("qreg q[2]; crz(2*pi) q[0],q[1];").unwrap();
    let (output, report) = optimize(
        &input,
        &Options {
            passes: Some(vec![PassName::DecomposeRotations]),
            rz_epsilon: 3.0,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(output.gates().is_empty());
    assert_eq!(report.numerical.certified_syntheses, 1);
    assert_eq!(report.numerical.uncertified_syntheses, 0);
}
