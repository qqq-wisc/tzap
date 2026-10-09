use num_bigint::BigInt;
use num_rational::BigRational;
use tzap::{
    angle::{Angle, AngleError, FiniteF64, PiFraction, exact_float_add},
    cancel::CancelGates,
    circuit::{Circuit, Gate},
    cnot_min::CnotMin,
    optimize::{Options, PassName, optimize},
    pass::Pass,
    phase_fold_pauli::PhaseFoldPauli,
    phase_fold_rand::PhaseFoldRand,
    qasm,
};
fn binary(value: f64) -> BigRational {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32;
    let mut significand = bits & ((1 << 52) - 1);
    if exponent != 0 {
        significand |= 1 << 52;
    }
    let power = if exponent == 0 {
        -1074
    } else {
        exponent - 1023 - 52
    };
    let mut n = BigInt::from(significand);
    if bits >> 63 != 0 {
        n = -n;
    }
    if power >= 0 {
        BigRational::from_integer(n << power as usize)
    } else {
        BigRational::new(n, BigInt::from(1) << (-power) as usize)
    }
}
fn numeric(value: f64) -> Angle {
    Angle::from_f64(value).unwrap()
}
fn circuit(angles: &[f64]) -> Circuit {
    let mut c = Circuit::new(1);
    for a in angles {
        c.apply(Gate::rz_f64(*a, 0).unwrap());
    }
    c
}
fn residual_sum(c: &Circuit) -> BigRational {
    c.gates
        .iter()
        .map(|g| match g {
            Gate::rz(a, _) => {
                let (p, r) = a.components().unwrap();
                assert_eq!(p.numerator(), 0);
                binary(r.get())
            }
            _ => panic!("numeric rotation must not become a named pi rotation: {g:?}"),
        })
        .sum()
}
#[test]
fn finite_inputs_and_fraction_boundaries() {
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(Angle::from_f64(x).is_err());
        assert!(Gate::rz_f64(x, 0).is_err());
    }
    assert_eq!(numeric(-0.0), numeric(0.0));
    assert_eq!(
        PiFraction::new(6, -8).unwrap(),
        PiFraction::new(-3, 4).unwrap()
    );
    assert_eq!(
        PiFraction::new(i64::MIN, i64::MIN).unwrap(),
        PiFraction::new(1, 1).unwrap()
    );
    assert!(PiFraction::new(1, i64::MIN).is_err());
    assert!(PiFraction::new(i64::MIN, 1).unwrap().checked_neg().is_err());
    assert_eq!(
        PiFraction::new(i64::MIN, i64::MAX)
            .unwrap()
            .normalized_rotation()
            .denominator(),
        i64::MAX
    );
    assert_eq!(
        Angle::pi_fraction(i64::MIN, 1).unwrap().quarter_turns(),
        Some(0)
    );
    assert_eq!(numeric(std::f64::consts::PI / 4.0).quarter_turns(), None);
    assert_eq!(Angle::pi_fraction(1, 4).unwrap().quarter_turns(), Some(1));
    assert_eq!(Angle::pi_fraction(1, 7).unwrap().quarter_turns(), None);
    for (n, d) in [
        (i64::MIN, 1),
        (i64::MIN, i64::MAX),
        (i64::MAX, i64::MAX - 1),
    ] {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(Angle::pi_fraction(n, d).unwrap(), 0));
        assert_eq!(qasm::parse(&c.to_qasm()).unwrap().gates, c.gates);
    }
}
#[test]
fn every_accepted_float_sum_matches_independent_binary_rationals() {
    let fixed = [
        0.0,
        -0.0,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::MAX,
        -f64::MAX,
        1e16,
        1.0,
        -1e16,
        1e308,
        0.1,
        0.2,
        0.3,
        0.7,
    ];
    let check = |a: f64, b: f64| {
        if let Ok(sum) = exact_float_add(FiniteF64::new(a).unwrap(), FiniteF64::new(b).unwrap()) {
            assert_eq!(binary(sum.get()), binary(a) + binary(b), "{a:?} + {b:?}");
        }
    };
    for a in fixed {
        for b in fixed {
            check(a, b);
        }
    }
    let mut seed = 0x6847487564u64;
    for _ in 0..20000 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        let a = f64::from_bits(seed);
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        let b = f64::from_bits(seed);
        if a.is_finite() && b.is_finite() {
            check(a, b);
        }
    }
    assert_eq!(
        exact_float_add(FiniteF64::new(1e16).unwrap(), FiniteF64::new(1.0).unwrap()),
        Err(AngleError::RoundedAddition)
    );
}
#[test]
fn fraction_arithmetic_matches_big_rational() {
    let mut seed = 723894723u64;
    for _ in 0..3000 {
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as i64
        };
        let a = PiFraction::new(next(), next() | 1);
        let b = PiFraction::new(next(), next() | 1);
        if let (Ok(a), Ok(b)) = (a, b) {
            let big = |p: PiFraction| {
                BigRational::new(BigInt::from(p.numerator()), BigInt::from(p.denominator()))
            };
            for (result, want) in [
                (a.checked_add(b), big(a) + big(b)),
                (a.checked_mul(b), big(a) * big(b)),
            ] {
                if let Ok(got) = result {
                    assert_eq!(big(got), want);
                }
            }
        }
    }
}
#[test]
fn parsing_pi_context_mixed_limits_fallback_and_multiline() {
    for (source, n, d) in [
        ("pi/4", 1, 4),
        ("0.1*pi", 1, 10),
        ("(1+2)*pi/8", 3, 8),
        ("(2*pi)/3", 2, 3),
        ("(pi/4)+(pi/4)", 1, 2),
        ("-pi/2", -1, 2),
    ] {
        let c = qasm::parse(&format!("OPENQASM 2.0; qreg q[1]; rz(\n{source}\n) q[0];")).unwrap();
        let Gate::rz(a, _) = &c.gates[0] else {
            panic!()
        };
        assert_eq!(a.components().unwrap().0, PiFraction::new(n, d).unwrap());
        let roundtrip = qasm::parse(&c.to_qasm()).unwrap();
        assert_eq!(roundtrip.gates, c.gates);
    }
    let c = qasm::parse("qreg q[1]; rz(pi/4+0.3) q[0];").unwrap();
    let Gate::rz(a, _) = &c.gates[0] else {
        panic!()
    };
    assert_eq!(a.components().unwrap().0, PiFraction::new(1, 4).unwrap());
    assert_eq!(a.components().unwrap().1.get(), 0.3);
    let output = c.to_qasm();
    assert!(output.contains("pi"));
    assert_eq!(qasm::parse(&output).unwrap().gates.len(), 2);
    let preserved = qasm::parse("qreg q[1]; rz(pi/9223372036854775808) q[0];").unwrap();
    assert!(matches!(&preserved.gates[0],Gate::rz(a,_) if a.is_preserved()));
    assert!(preserved.to_qasm().contains("pi/9223372036854775808"));
    let c = qasm::parse("qreg q[1]; rz(pi*pi/4) q[0];").unwrap();
    assert!(matches!(&c.gates[0],Gate::rz(a,_) if a.used_numeric_fallback()));
    assert!(qasm::parse("qreg q[1]; rz((1/0)-(1/0)) q[0];").is_err());
    assert!(qasm::parse("qreg q[1]; rz((1/0)*pi) q[0];").is_err());
    let c = qasm::parse("qreg q[1]; rz(0.1+0.2) q[0];").unwrap();
    assert!(matches!(&c.gates[0],Gate::rz(a,_) if a.components().unwrap().1.get()==0.1+0.2));
}
#[test]
fn every_numeric_folder_preserves_known_audit_regressions() {
    let passes: [&dyn Pass; 4] = [
        &CancelGates,
        &PhaseFoldRand,
        &PhaseFoldPauli,
        &CnotMin::default(),
    ];
    for values in [
        &[5e-7][..],
        &[1e16, 1.0, -1e16],
        &[1e308, 1e308],
        &[1e16],
        &[0.1, 0.2],
        &[f64::from_bits(1), f64::from_bits(1)],
    ] {
        let input = circuit(values);
        let want = residual_sum(&input);
        for pass in passes {
            let out = pass.run(&input);
            assert_eq!(residual_sum(&out), want, "{}: {values:?}", pass.name());
        }
    }
}
#[test]
fn exact_pi_and_long_named_sequences() {
    let input = qasm::parse("qreg q[1]; rz(pi/7) q[0]; rz(2*pi/7) q[0];").unwrap();
    for pass in [
        &PhaseFoldRand as &dyn Pass,
        &PhaseFoldPauli,
        &CnotMin::default(),
    ] {
        let out = pass.run(&input);
        assert_eq!(out.gates.len(), 1);
        assert!(
            matches!(&out.gates[0],Gate::rz(a,_) if a.components().unwrap().0==PiFraction::new(3,7).unwrap())
        );
    }
    let mut input = Circuit::new(1);
    input.gates = vec![Gate::t(0); 10000];
    for pass in [
        &PhaseFoldRand as &dyn Pass,
        &PhaseFoldPauli,
        &CnotMin::default(),
    ] {
        let out = pass.run(&input);
        let quarters: u32 = out
            .gates
            .iter()
            .map(|g| match g {
                Gate::t(_) => 1,
                Gate::tdg(_) => 7,
                Gate::s(_) => 2,
                Gate::sdg(_) => 6,
                Gate::z(_) => 4,
                _ => panic!("named phases must remain discrete"),
            })
            .sum();
        assert_eq!(quarters % 8, 0);
    }
    let c = qasm::parse("qreg q[1]; rz(pi/4+5e-7) q[0];").unwrap();
    assert!(
        PhaseFoldRand
            .run(&c)
            .gates
            .iter()
            .any(|g| matches!(g,Gate::rz(a,_) if a.components().unwrap().1.get()==5e-7))
    );
}
#[test]
fn reports_separate_rounding_and_random_matching_in_serial_and_parallel() {
    let input = circuit(&[1e16, 1.0, -1e16]);
    for parallel in [false, true] {
        let options = Options {
            passes: Some(vec![PassName::PhaseFoldRand]),
            parallel,
            ..Options::default()
        };
        let (out, report) = optimize(&input, &options).unwrap();
        assert_eq!(residual_sum(&out), residual_sum(&input));
        assert!(report.numerical.randomized_matching);
        if !parallel {
            assert!(report.numerical.skipped_rounded_folds > 0);
        }
    }
}
#[test]
fn pbc_accepts_only_certified_quarters() {
    let exact = qasm::parse("qreg q[1]; rz(pi/4) q[0];").unwrap();
    assert!(tzap::pbc::to_pbc(&exact, None).is_ok());
    assert!(tzap::pbc::to_pbc(&circuit(&[std::f64::consts::PI / 4.0]), None).is_err());
    let fraction = qasm::parse("qreg q[1]; rz(pi/7) q[0];").unwrap();
    assert!(tzap::pbc::to_pbc(&fraction, None).is_err());
    // Lowering a zero rotation must not hide an invalid input qubit.
    let mut invalid = Circuit::new(1);
    invalid.apply(Gate::rz(Angle::pi_fraction(0, 1).unwrap(), 1));
    assert!(tzap::pbc::to_pbc(&invalid, None).is_err());
}
#[test]
fn synthesis_errors_and_exact_lowering() {
    let preserved = qasm::parse("qreg q[1]; rz(pi/9223372036854775808) q[0];").unwrap();
    let options = Options {
        passes: Some(vec![PassName::DecomposeRz]),
        ..Options::default()
    };
    assert!(optimize(&preserved, &options).is_err());
    let exact = qasm::parse("qreg q[1]; rz(pi/4) q[0];").unwrap();
    let (out, report) = optimize(&exact, &options).unwrap();
    assert_eq!(out.gates, vec![Gate::t(0)]);
    assert_eq!(report.numerical.uncertified_syntheses, 0);
    for epsilon in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            tzap::decompose::DecomposeRz { epsilon }
                .try_run(&exact)
                .is_err()
        );
    }
    // Loose explicit requests must not underflow gridsynth's precision setup.
    assert!(
        tzap::decompose::DecomposeRz { epsilon: 10.0 }
            .try_run(&circuit(&[0.3]))
            .is_ok()
    );
}

#[test]
fn parser_complexity_is_bounded_without_changing_normal_expressions() {
    let nested = format!("{}pi{}", "(".repeat(200), ")".repeat(200));
    let chain = vec!["pi"; 300].join("+");
    for expression in [nested, chain] {
        let error = qasm::parse(&format!("qreg q[1]; rz({expression}) q[0];")).unwrap_err();
        assert!(error.to_string().contains("complexity limit"));
    }
    assert!(qasm::parse("qreg q[1]; rz((((2*pi)/3))) q[0];").is_ok());
}

#[test]
fn signed_failed_merges_and_coefficient_limits_are_transactional() {
    let sum = |c: &Circuit| {
        let mut flip = false;
        let mut radians = BigRational::from_integer(BigInt::from(0));
        for gate in &c.gates {
            match gate {
                Gate::x(_) => flip = !flip,
                Gate::rz(angle, _) => {
                    let (p, r) = angle.components().unwrap();
                    assert_eq!(p.numerator(), 0);
                    radians += if flip {
                        -binary(r.get())
                    } else {
                        binary(r.get())
                    };
                }
                _ => panic!("unexpected numeric-axis operation"),
            }
        }
        (flip, radians)
    };
    let mut input = circuit(&[1e16]);
    input.apply(Gate::x(0));
    input.apply(Gate::rz_f64(1.0, 0).unwrap());
    input.apply(Gate::x(0));
    input.apply(Gate::rz_f64(-1e16, 0).unwrap());
    for pass in [
        &CancelGates as &dyn Pass,
        &PhaseFoldRand,
        &PhaseFoldPauli,
        &CnotMin::default(),
    ] {
        assert_eq!(sum(&pass.run(&input)), sum(&input), "{}", pass.name());
    }
    let input = Circuit {
        num_qubits: 1,
        num_cbits: 0,
        gates: vec![
            Gate::rz(Angle::pi_fraction(1, i64::MAX).unwrap(), 0),
            Gate::rz(Angle::pi_fraction(1, i64::MAX - 1).unwrap(), 0),
        ],
    };
    for pass in [
        &PhaseFoldRand as &dyn Pass,
        &PhaseFoldPauli,
        &CnotMin::default(),
    ] {
        assert_eq!(pass.run(&input).gates, input.gates, "{}", pass.name());
    }
}

#[test]
fn canonical_equality_hashing_layout_and_numeric_export() {
    use std::hash::{Hash, Hasher};
    let exact = Angle::pi_fraction(1, 4).unwrap();
    let equivalent = Angle::pi_fraction(9, 4).unwrap();
    assert_eq!(exact, equivalent);
    let hash = |a: &Angle| {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        a.hash(&mut h);
        h.finish()
    };
    assert_eq!(hash(&exact), hash(&equivalent));
    println!(
        "Angle: {} bytes; Gate: {} bytes",
        std::mem::size_of::<Angle>(),
        std::mem::size_of::<Gate>()
    );
    let input = qasm::parse("qreg q[1]; rz(pi/7) q[0]; rz(pi/4) q[0]; rz(0.3) q[0];").unwrap();
    let export = qasm::serialize_numeric_lossy(&input).unwrap();
    assert_eq!(export.uncertified_conversions, 1);
    assert!(!export.qasm.contains("pi"));
    assert!(export.qasm.contains("t q[0]"));
    assert!(
        qasm::serialize_numeric_lossy(
            &qasm::parse("qreg q[1]; rz(pi/9223372036854775808) q[0];").unwrap()
        )
        .is_err()
    );
}
