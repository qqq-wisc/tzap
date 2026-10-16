use super::*;
use num_rational::BigRational;
use num_traits::Zero;

fn rational_endpoint(n: &BigInt, ctx: &Context) -> BigRational {
    BigRational::new(n.clone(), ctx.unit.clone())
}
#[test]
fn interval_multiplication_encloses_independent_rational_endpoints() {
    let ctx = Context::new(3);
    for a in -10..10 {
        for b in -10..10 {
            let left = Interval {
                lo: a.into(),
                hi: (a + 3).into(),
            };
            let right = Interval {
                lo: b.into(),
                hi: (b + 2).into(),
            };
            let result = ctx.mul(&left, &right);
            for x in [&left.lo, &left.hi] {
                for y in [&right.lo, &right.hi] {
                    let exact = rational_endpoint(x, &ctx) * rational_endpoint(y, &ctx);
                    assert!(rational_endpoint(&result.lo, &ctx) <= exact);
                    assert!(rational_endpoint(&result.hi, &ctx) >= exact);
                }
            }
        }
    }
}
#[test]
fn mathematical_pi_and_trigonometric_extrema_are_enclosed() {
    let ctx = Context::new(256);
    let pi = ctx.pi();
    let digits = BigInt::parse_bytes(b"314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651", 10).unwrap();
    let decimal = BigRational::new(digits.clone(), BigInt::from(10).pow(110));
    let upper = BigRational::new(digits + 1, BigInt::from(10).pow(110));
    assert!(rational_endpoint(&pi.lo, &ctx) < decimal);
    assert!(rational_endpoint(&pi.hi, &ctx) > upper);
    let (sin, cos) = ctx.sin_cos(&pi.div_int(2)).unwrap();
    assert!(sin.lo <= ctx.unit && sin.hi >= ctx.unit);
    assert!(cos.lo <= BigInt::zero() && cos.hi >= BigInt::zero());
}
#[test]
fn verifier_rejects_wrong_words_and_unknown_backend_instructions() {
    let target = prepare(&Angle::from_f64(1.0).unwrap(), 256).unwrap();
    let budget = target.ctx.float(1e-20);
    assert!(!interval::verifies(&target.matrix, &[], &target.pi, &budget, &target.ctx).unwrap());
    assert!(!interval::verifies(&target.matrix, &['T'], &target.pi, &budget, &target.ctx).unwrap());
    assert!(interval::verifies(&target.matrix, &['?'], &target.pi, &budget, &target.ctx).is_err());
}
#[test]
fn verifier_accounts_for_global_phase_and_matrix_order() {
    let ctx = Context::new(256);
    let pi = ctx.pi();
    let budget = ctx.float(1e-20);
    let target = prepare(&Angle::pi_fraction(1, 4).unwrap(), 256).unwrap();
    assert!(interval::verifies(&target.matrix, &['T'], &pi, &budget, &ctx).unwrap());
    assert!(interval::verifies(&target.matrix, &['W', 'T', 'W'], &pi, &budget, &ctx).unwrap());
    let (hx, _) = interval::word_matrix(&['H', 'X'], &ctx).unwrap();
    let (xh, _) = interval::word_matrix(&['X', 'H'], &ctx).unwrap();
    assert!(hx[0][1].re.hi < 0.into());
    assert!(xh[0][1].re.lo > 0.into());
}

#[test]
fn accepted_huge_angles_match_independent_high_precision_fixtures() {
    // mpmath 1.3.0, 500 decimal digits, exact integer inputs; fixtures are
    // supplementary checks of the interval kernel, not production certificates.
    for (angle, bits, cosine, sine) in [
        (
            1e16,
            512,
            "0.432337716297637990517139024958024763417702302473549230584430093556291207557036725373332489328520567698474100922540335348",
            "0.901711760523584989956907790069803117732676668297127938633487125420365399298212213272987561401453996348484605755449745232",
        ),
        (
            f64::MAX,
            1536,
            "0.0024809850301908852035536887886555882383323977871115371060025290036519701102420897941058631037120157652075327343270151243",
            "0.999996922351904025464889675937959275072227326352092025600555207856153938612148485442266231194258125945572427389046252709",
        ),
    ] {
        let target = prepare(&Angle::from_f64(angle).unwrap(), bits).unwrap();
        for (value, spelling) in [
            (&target.matrix[0][0].re, cosine),
            (&target.matrix[1][1].im, sine),
        ] {
            let decimals = spelling.len() - 2;
            let numerator = BigInt::parse_bytes(spelling.replace('.', "").as_bytes(), 10).unwrap();
            let denominator = BigInt::from(10).pow(decimals as u32);
            let low = BigRational::new(&numerator - 1, denominator.clone());
            let high = BigRational::new(numerator + 1, denominator);
            let abs = if value.lo < BigInt::zero() {
                value.neg()
            } else {
                value.clone()
            };
            assert!(rational_endpoint(&abs.lo, &target.ctx) <= high);
            assert!(rational_endpoint(&abs.hi, &target.ctx) >= low);
            assert!(&abs.hi - &abs.lo < &target.ctx.unit / BigInt::from(10).pow(80));
        }
    }
}

#[test]
fn full_controlled_verifier_detects_lost_branch_phase() {
    let angle = Angle::pi_fraction(2, 1).unwrap();
    let crz = Gate::crz {
        theta: angle.clone(),
        control: 0,
        target: 1,
    };
    let cp = Gate::cp {
        lambda: angle,
        control: 0,
        target: 1,
    };
    assert!(verify_controlled(&crz, &[], 1e-20).is_err());
    assert!(verify_controlled(&cp, &[], 1e-20).is_ok());
}

#[test]
fn complex_distance_bounds_include_zero_crossings_and_negative_components() {
    let other = Complex {
        re: Interval::point(1.into()),
        im: Interval::point((-1).into()),
    };
    for re in -4i64..=4 {
        for im in -4i64..=4 {
            let rectangle = Complex {
                re: Interval {
                    lo: re.into(),
                    hi: (re + 3).into(),
                },
                im: Interval {
                    lo: im.into(),
                    hi: (im + 3).into(),
                },
            };
            // The extrema of this integer-endpoint rectangle occur at an
            // endpoint or zero; exhaust its integer grid independently.
            let values: Vec<_> = (re..=re + 3)
                .flat_map(|r| (im..=im + 3).map(move |i| (r - 1).pow(2) + (i + 1).pow(2)))
                .collect();
            let (lower, upper) = rectangle.difference_square_bounds(&other);
            assert_eq!(lower, BigInt::from(*values.iter().min().unwrap()));
            assert_eq!(upper, BigInt::from(*values.iter().max().unwrap()));
        }
    }
}

#[test]
fn controlled_symbolic_synthesis_retries_target_preparation() {
    use crate::{decompose::DecomposeRotations, qasm, unitary::circuits_equiv};
    // Independent mpmath fixtures at 500 decimal digits, reduced modulo 4*pi.
    // Never evaluate these huge intermediate arguments through binary64.
    for (source, reduced) in [
        ("exp(256)", "3.792252811227286521400825"),
        ("sin(exp(256))", "-0.6057118146249979501423651"),
        ("cos(exp(256))", "-0.7956841066803409008225726"),
    ] {
        let angle = crate::angle::AngleContext::default().parse(source).unwrap();
        assert!(prepare_with_period(&angle, 209, 4).is_err());
        for family in ["cp", "crx", "cry", "crz"] {
            for operands in ["q[0],q[1]", "q[1],q[0]"] {
                let input =
                    qasm::parse(&format!("qreg q[2]; {family}({source}) {operands};")).unwrap();
                let reference =
                    qasm::parse(&format!("qreg q[2]; {family}({reduced}) {operands};")).unwrap();
                let output = DecomposeRotations { epsilon: 1e-5 }
                    .try_run(&input)
                    .unwrap();
                assert!(output.gates.iter().all(|g| g.angle().is_none()));
                assert!(
                    circuits_equiv(&reference, &output, 1e-5),
                    "{family}({source}) {operands}"
                );
                assert!(verify_controlled(&input.gates[0], &[Gate::z(0)], 1e-5).is_err());
            }
        }
    }
}

#[test]
fn controlled_verifier_retries_inconclusive_distance_and_rejects_wrong_words() {
    use crate::{decompose::DecomposeRotations, qasm, unitary::circuits_equiv};
    let input = qasm::parse("qreg q[2]; crz(exp(128)) q[1],q[0];").unwrap();
    let angle = input.gates[0].angle().unwrap();
    let output = DecomposeRotations { epsilon: 1e-5 }
        .try_run(&input)
        .unwrap();
    let exact = controlled_word(&output.gates, 1, 0).unwrap();
    // Preparation succeeds at the initial precision, but the distance bound
    // remains too wide. This exercises retries after candidate evaluation.
    assert!(matches!(
        verify_controlled_at_precision(&input.gates[0], angle, &exact, 1e-5, 209).unwrap(),
        ControlledVerification::Inconclusive
    ));
    verify_controlled(&input.gates[0], &output.gates, 1e-5).unwrap();
    let reference = qasm::parse("qreg q[2]; crz(3.943723543226126377318808) q[1],q[0];").unwrap();
    assert!(circuits_equiv(&reference, &output, 1e-5));

    let wrong = qasm::parse("qreg q[2]; crz(pi) q[1],q[0];").unwrap();
    let angle = wrong.gates[0].angle().unwrap();
    assert!(matches!(
        verify_controlled_at_precision(
            &wrong.gates[0],
            angle,
            &exact_word::ExactWord::identity(),
            1e-5,
            209
        )
        .unwrap(),
        ControlledVerification::Rejected
    ));
}

#[test]
fn controlled_symbolic_synthesis_still_fails_at_resource_limits() {
    use crate::{decompose::DecomposeRotations, pass::Pass, qasm};
    let input = qasm::parse("qreg q[2]; crz(sin(exp(1025))) q[1],q[0];").unwrap();
    let pass = DecomposeRotations { epsilon: 1e-5 };
    assert!(
        pass.try_run(&input)
            .unwrap_err()
            .to_string()
            .contains("symbolic target could not be enclosed")
    );
    assert_eq!(pass.run(&input).gates, input.gates);
}

#[test]
fn numerical_limits_fail_without_accepting_unverified_output() {
    assert!(
        synthesize(&Angle::from_f64(1.0).unwrap(), 1e-100)
            .unwrap_err()
            .to_string()
            .contains("precision limit")
    );
}

#[test]
fn long_exact_words_enclose_the_independent_rational_scalar_oracle() {
    use crate::{
        circuit::Circuit,
        semantics::{self, Limits},
    };
    let ctx = Context::new(256);
    let mut circuit = Circuit::new(1);
    let mut word = Vec::new();
    for _ in 0..60 {
        word.extend(['H', 'T', 'S', 'X']);
        circuit
            .gates
            .extend([Gate::h(0), Gate::t(0), Gate::s(0), Gate::x(0)]);
    }
    let expected = semantics::circuit_unitary(&circuit, Limits::default()).unwrap();
    let (actual, _) = interval::word_matrix(&word, &ctx).unwrap();
    // Independently evaluate Q(sqrt(2), i) components using an enclosing sqrt2.
    let root = ctx.inv_sqrt_two().scale(&2.into());
    let enclose_rational = |r: &BigRational| Interval {
        lo: (r.numer() * &ctx.unit).div_floor(r.denom()),
        hi: (r.numer() * &ctx.unit).div_ceil(r.denom()),
    };
    for (i, row) in actual.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            let scalar = &expected.get(i, j).0;
            let real =
                enclose_rational(&scalar[0]).add(&ctx.mul(&enclose_rational(&scalar[1]), &root));
            let imag =
                enclose_rational(&scalar[2]).add(&ctx.mul(&enclose_rational(&scalar[3]), &root));
            for (a, b) in [(&entry.re, real), (&entry.im, imag)] {
                assert!(a.lo <= b.hi && a.hi >= b.lo);
                assert!(&a.hi - &a.lo < 100.into());
            }
        }
    }
}

#[test]
fn tiny_hadamard_rotation_is_preserved_by_each_exact_pass() {
    use crate::{
        cancel::CancelGates, circuit::Circuit, cnot_min::CnotMin, pass::Pass,
        phase_fold_pauli::PhaseFoldPauli, phase_fold_rand::PhaseFoldRand, unitary::circuits_equiv,
    };
    let mut circuit = Circuit::new(1);
    circuit.gates = vec![Gate::h(0), Gate::rz_f64(5e-10, 0).unwrap(), Gate::h(0)];
    let identity = Circuit::new(1);
    assert!(!circuits_equiv(&circuit, &identity, 1e-14));
    for pass in [
        &CancelGates as &dyn Pass,
        &PhaseFoldRand,
        &PhaseFoldPauli,
        &CnotMin::default(),
    ] {
        assert!(
            circuits_equiv(&circuit, &pass.run(&circuit), 1e-14),
            "{}",
            pass.name()
        );
    }
}

#[test]
fn exact_scaled_targets_enclose_unrepresentable_halves() {
    let tiny = Angle::from_f64(f64::from_bits(1)).unwrap();
    assert_eq!(
        tiny.checked_half()
            .unwrap()
            .as_rational()
            .unwrap()
            .to_big_rational(),
        BigRational::new(1.into(), BigInt::from(1) << 1075usize)
    );
    let bits = 1200;
    let target = prepare_scaled_with_period(&tiny, bits, 2, Scale::Half).unwrap();
    let expected = BigRational::new(1.into(), BigInt::from(1) << 1075usize);
    assert!(rational_endpoint(&target.reduced.lo, &target.ctx) <= expected);
    assert!(rational_endpoint(&target.reduced.hi, &target.ctx) >= expected);
    let negative = prepare_scaled_with_period(&tiny, bits, 2, Scale::NegativeHalf).unwrap();
    assert!(rational_endpoint(&negative.reduced.lo, &negative.ctx) <= -expected.clone());
    assert!(rational_endpoint(&negative.reduced.hi, &negative.ctx) >= -expected);
    let fraction = Angle::pi_fraction(1, i64::MAX).unwrap();
    assert_eq!(
        fraction
            .checked_half()
            .unwrap()
            .pi_coefficient()
            .unwrap()
            .to_big_rational(),
        BigRational::new(1.into(), BigInt::from(i64::MAX) * 2)
    );
    let target = prepare_scaled_with_period(&fraction, 512, 2, Scale::Half).unwrap();
    let pi = target.ctx.pi();
    let denominator = BigInt::from(i64::MAX) * 2;
    let lower = rational_endpoint(&pi.lo, &target.ctx) / &denominator;
    let upper = rational_endpoint(&pi.hi, &target.ctx) / &denominator;
    assert!(rational_endpoint(&target.reduced.lo, &target.ctx) <= lower);
    assert!(rational_endpoint(&target.reduced.hi, &target.ctx) >= upper);
}

#[test]
fn scaled_synthesis_cache_keeps_literal_sign_and_scale() {
    let angle = Angle::from_f64(0.7).unwrap();
    // The decompose entry point normally holds the process-global backend lock.
    let _guard = crate::decompose::GRIDSYNTH_BATCH_LOCK
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    for scale in [
        Scale::Whole,
        Scale::Half,
        Scale::NegativeHalf,
        Scale::Whole,
        Scale::NegativeHalf,
    ] {
        let word = synthesize_scaled(&angle, scale, 1e-10).unwrap();
        let target = prepare_scaled_with_period(&angle, 512, 2, scale).unwrap();
        assert!(
            interval::verifies(
                &target.matrix,
                &word,
                &target.pi,
                &target.ctx.float(1e-10),
                &target.ctx
            )
            .unwrap()
        );
        let wrong = prepare_scaled_with_period(
            &angle,
            512,
            2,
            if scale == Scale::Whole {
                Scale::Half
            } else {
                Scale::Whole
            },
        )
        .unwrap();
        assert!(
            !interval::verifies(
                &wrong.matrix,
                &word,
                &wrong.pi,
                &wrong.ctx.float(1e-10),
                &wrong.ctx
            )
            .unwrap()
        );
    }
}

#[test]
fn tiny_rationals_with_long_denominators_synthesize_to_the_identity() {
    let angle: Angle = "1e-600".parse().unwrap();
    assert!(angle.coefficient_bits() > MAX_BITS - 128);
    assert_eq!(synthesize(&angle, 1e-10).unwrap(), Vec::<char>::new());
    let circuit = crate::qasm::parse("qreg q[1]; rz(1e-600) q[0];").unwrap();
    let output = crate::decompose::DecomposeRotations::default()
        .try_run(&circuit)
        .unwrap();
    assert!(output.gates.is_empty(), "{:?}", output.gates);
}
