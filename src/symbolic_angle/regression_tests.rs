use super::*;
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed};

fn next(seed: &mut u64) -> u64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed
}
fn exact(a: &Angle) -> BigRational {
    a.as_rational().unwrap().to_big_rational()
}
fn binary(value: f64) -> BigRational {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32;
    let mantissa = bits & ((1 << 52) - 1);
    let mut n = BigInt::from(if exponent == 0 {
        mantissa
    } else {
        mantissa + (1 << 52)
    });
    if bits >> 63 != 0 {
        n = -n;
    }
    let power = if exponent == 0 {
        -1074
    } else {
        exponent - 1075
    };
    if power < 0 {
        BigRational::new(n, BigInt::one() << -power as usize)
    } else {
        BigRational::from_integer(n << power as usize)
    }
}
fn hash(a: &impl Hash) -> u64 {
    let mut h = FxHasher::default();
    a.hash(&mut h);
    h.finish()
}

#[test]
fn binary64_pair_arithmetic_matches_independent_exact_oracle() {
    let mut ctx = AngleContext::default();
    let mut seed = 0xed12_3467_983a_7851;
    for i in 0..2000 {
        let a_bits = next(&mut seed);
        let raw = next(&mut seed);
        let b_bits = if i % 2 == 0 {
            // Concentrate half the cases around a's exponent to exercise exact sums.
            (raw & !(2047 << 52)) | (a_bits & (2047 << 52))
        } else {
            raw
        };
        let (a, b) = (f64::from_bits(a_bits), f64::from_bits(b_bits));
        if !a.is_finite() || !b.is_finite() {
            continue;
        }
        let aa = Angle::from_f64(a).unwrap();
        let bb = Angle::from_f64(b).unwrap();
        let expected = binary(a) + binary(b);
        let result = ctx.add(&aa, &bb).unwrap();
        assert_eq!(exact(&result), expected, "{a_bits:016x} + {b_bits:016x}");
        assert_eq!(ctx.add(&bb, &aa).unwrap(), result);
        assert_eq!(exact(&ctx.sub(&aa, &bb).unwrap()), binary(a) - binary(b));
        assert_eq!(exact(&ctx.neg(&aa).unwrap()), -binary(a));
        assert_eq!(exact(&ctx.half(&aa).unwrap()), binary(a) / BigInt::from(2));
        let canonical = ctx
            .big_rational(expected.numer().clone(), expected.denom().clone())
            .unwrap();
        assert_eq!(hash(&result), hash(&canonical));
        assert_eq!(result.cmp(&canonical), Ordering::Equal);
        assert_eq!(
            result.coefficient_bits(),
            result.as_rational().unwrap().bits()
        );
    }
}

#[test]
fn binary64_boundary_pairs_cover_rounding_cancellation_and_overflow() {
    let mut ctx = AngleContext::default();
    let values = [
        0.0,
        -0.0,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::from_bits((1 << 52) - 1),
        f64::MIN_POSITIVE,
        f64::from_bits((1 << 52) + 1),
        0.1,
        0.5,
        1.0,
        f64::from_bits(1.0_f64.to_bits() + 1),
        1e16,
        2.0_f64.powi(53),
        f64::MAX,
        -f64::MAX,
    ];
    for a in values {
        for b in values {
            let aa = Angle::from_f64(a).unwrap();
            let bb = Angle::from_f64(b).unwrap();
            let result = ctx.add(&aa, &bb).unwrap();
            let expected = binary(a) + binary(b);
            assert_eq!(exact(&result), expected, "{a:?} + {b:?}");
            let rounded = a + b;
            if rounded.is_finite() && binary(rounded) == expected {
                assert_eq!(result, Angle::from_f64(rounded).unwrap());
            }
        }
    }
    let tiny = Angle::from_f64(f64::from_bits(1)).unwrap();
    let half = ctx.half(&tiny).unwrap();
    assert!(!half.is_inline());
    assert_eq!(ctx.add(&half, &half).unwrap(), tiny);
}

#[test]
fn randomized_large_rationals_stay_reduced_through_arithmetic() {
    let mut ctx = AngleContext::default();
    let mut seed = 0x7863_3189_2210_4205;
    for i in 0_usize..250 {
        let a = BigInt::from(next(&mut seed)) << (65 + i % 300);
        let a = if i % 2 == 0 { a } else { -a };
        let d = BigInt::from(next(&mut seed) | 1);
        let b = BigInt::from(next(&mut seed)) << (70 + i % 200);
        let e = BigInt::from(next(&mut seed) | 1);
        let aa = ctx.big_rational(a.clone(), d.clone()).unwrap();
        let bb = ctx.big_rational(b.clone(), e.clone()).unwrap();
        let a = BigRational::new(a, d);
        let b = BigRational::new(b, e);
        for (result, expected) in [
            (ctx.add(&aa, &bb).unwrap(), &a + &b),
            (ctx.mul(&aa, &bb).unwrap(), &a * &b),
            (ctx.div(&aa, &bb).unwrap(), &a / &b),
            (ctx.add(&aa, &Angle::ONE).unwrap(), &a + BigRational::one()),
            (ctx.add(&Angle::ONE, &aa).unwrap(), &a + BigRational::one()),
        ] {
            let actual = exact(&result);
            assert_eq!(actual, expected);
            assert!(actual.denom().is_positive());
            assert_eq!(actual.numer().gcd(actual.denom()), BigInt::one());
            let rebuilt = ctx
                .big_rational(actual.numer().clone(), actual.denom().clone())
                .unwrap();
            assert_eq!(hash(&result), hash(&rebuilt));
        }
    }
}

#[test]
fn decimal_machine_boundaries_and_scientific_notation_are_exact() {
    let mut ctx = AngleContext::default();
    for literal in [
        "0",
        "-0",
        "+.0",
        "000000000.00000001",
        "1e-19",
        "1e19",
        "-9223372036854775808",
        "9223372036854775808",
        "18446744073709551615",
        "18446744073709551616",
        "-18446744073709551615",
        "10e-19",
        "0.00000000000000000001",
        "0e9223372036854775807",
        "-0e-9223372036854775808",
    ] {
        let direct = ctx.decimal(literal).unwrap();
        assert_eq!(ctx.parse(literal).unwrap(), direct, "{literal}");
        assert_eq!(direct, ctx.parse(&direct.to_string()).unwrap(), "{literal}");
    }
    let mut seed = 0xd486_8732_e687_ab24;
    for _ in 0..500 {
        let n = next(&mut seed) as i64;
        let exponent = (next(&mut seed) % 41) as i32 - 20;
        let literal = format!("{n}e{exponent}");
        let power = BigInt::from(10).pow(exponent.unsigned_abs());
        let expected = if exponent < 0 {
            BigRational::new(n.into(), power)
        } else {
            BigRational::from_integer(BigInt::from(n) * power)
        };
        assert_eq!(
            exact(&ctx.decimal(&literal).unwrap()),
            expected,
            "{literal}"
        );
    }
}

#[test]
fn trig_period_parity_and_signs_hold_on_rational_pi_grid() {
    let mut ctx = AngleContext::default();
    for d in 1..=12 {
        for n in -24_i64..=24 {
            let a = ctx.pi_fraction(n, d).unwrap();
            let neg = ctx.pi_fraction(-n, d).unwrap();
            let shifted = ctx.pi_fraction(n + 2 * d, d).unwrap();
            let sine = ctx.sin(&a).unwrap();
            let cosine = ctx.cos(&a).unwrap();
            assert_eq!(sine, ctx.sin(&shifted).unwrap());
            assert_eq!(cosine, ctx.cos(&shifted).unwrap());
            let neg_sine = ctx.neg(&sine).unwrap();
            assert_eq!(neg_sine, ctx.sin(&neg).unwrap());
            assert_eq!(cosine, ctx.cos(&neg).unwrap());
            let ss = ctx.pow(&sine, 2).unwrap();
            let cc = ctx.pow(&cosine, 2).unwrap();
            assert_eq!(ctx.add(&ss, &cc).unwrap(), Angle::ONE);
            let r = n.rem_euclid(2 * d);
            assert_eq!(
                sine.prove_zero(),
                if r == 0 || r == d {
                    Proof::Proven
                } else {
                    Proof::Disproven
                }
            );
            assert_eq!(
                sine.prove_positive(),
                if r > 0 && r < d {
                    Proof::Proven
                } else {
                    Proof::Disproven
                }
            );
            let tangent = ctx.tan(&a);
            if (2 * n).rem_euclid(2 * d) == d {
                assert_eq!(tangent, Err(Error::InvalidDomain(Function::Tan)));
            } else {
                let shifted = ctx.pi_fraction(n + d, d).unwrap();
                assert_eq!(tangent.unwrap(), ctx.tan(&shifted).unwrap());
            }
        }
    }
}

#[test]
fn random_symbolic_linear_combinations_match_independent_coefficients() {
    let mut ctx = AngleContext::default();
    let atoms: Vec<_> = ["sin(pi/7)", "cos(1/3)", "exp(1/5)", "ln(2)", "sqrt(2)"]
        .map(|s| ctx.parse(s).unwrap())
        .into();
    let mut seed = 0x1818_6677_eeff_4321;
    for _ in 0..200 {
        let mut coefficients = [0_i64; 5];
        let mut value = Angle::ZERO;
        for _ in 0..25 {
            let index = (next(&mut seed) % 5) as usize;
            let coefficient = (next(&mut seed) % 21) as i64 - 10;
            coefficients[index] += coefficient;
            let term = ctx
                .scale(&atoms[index], &Rational::integer(coefficient))
                .unwrap();
            value = ctx.add(&value, &term).unwrap();
        }
        let mut expected = Angle::ZERO;
        for index in (0..5).rev() {
            let term = ctx
                .scale(&atoms[index], &Rational::integer(coefficients[index]))
                .unwrap();
            expected = ctx.add(&expected, &term).unwrap();
        }
        assert_eq!(value, expected);
        assert_eq!(hash(&value), hash(&expected));
        let restored: Angle = value.to_expression_string(10_000).unwrap().parse().unwrap();
        assert_eq!(restored, value);
    }
}

#[test]
fn product_power_grouping_has_a_consistent_normal_form() {
    let mut ctx = AngleContext::default();
    for (a, b) in [
        ("((sin(1)*sin(1))*sin(1))*sin(1)", "sin(1)^4"),
        ("sin(1)^2*sin(1)^3", "sin(1)^5"),
        ("sin(1)^-2*sin(1)^3", "sin(1)"),
        ("3*sin(1)^2*2*sin(1)^3", "6*sin(1)^5"),
        ("sin(1)^-2*sin(1)^2", "1"),
    ] {
        let a = ctx.parse(a).unwrap();
        let b = ctx.parse(b).unwrap();
        assert_eq!(a, b);
    }
}

#[test]
fn unproved_signs_preserve_absolute_square_roots_and_partial_domains() {
    let mut ctx = AngleContext::default();
    let unknown = ctx.parse("sin(exp(10))").unwrap();
    let square = ctx.pow(&unknown, 2).unwrap();
    assert_eq!(square.prove_nonnegative(), Proof::Proven);
    assert_eq!(square.prove_positive(), Proof::Unknown);
    let root = ctx.sqrt(&square).unwrap();
    assert_ne!(root, unknown);
    assert_eq!(ctx.pow(&root, 2).unwrap(), square);
    assert_eq!(ctx.div(&Angle::ONE, &root), Err(Error::DomainUnproved));
    assert_eq!(ctx.ln(&root), Err(Error::DomainUnproved));
    assert_eq!(ctx.pow(&unknown, 0), Err(Error::DomainUnproved));
    let exponential = ctx.exp(&unknown).unwrap();
    assert_eq!(ctx.ln(&exponential).unwrap(), unknown);
    assert_eq!(ctx.div(&exponential, &exponential).unwrap(), Angle::ONE);
}

#[test]
fn interner_capacity_hits_cleanup_and_clear_preserve_live_values() {
    let mut ctx = AngleContext::with_limits(Limits {
        max_nodes: 2,
        ..Limits::default()
    });
    let sine = ctx.sin(&Angle::ONE).unwrap();
    let cosine = ctx.cos(&Angle::ONE).unwrap();
    assert!(ctx.sin(&Angle::ONE).unwrap().shares_storage(&sine));
    assert_eq!(ctx.exp(&Angle::ONE), Err(Error::Limit(Resource::Nodes)));
    drop(cosine);
    let exponential = ctx.exp(&Angle::ONE).unwrap();
    assert_eq!(ctx.stats().retained_entries, 2);
    let saved = sine.to_string();
    let misses = ctx.stats().misses;
    ctx.clear_interner();
    assert_eq!(sine.to_string(), saved);
    let rebuilt = ctx.sin(&Angle::ONE).unwrap();
    assert_eq!(rebuilt, sine);
    assert!(!rebuilt.shares_storage(&sine));
    assert!(ctx.stats().misses > misses);
    assert_eq!(exponential.prove_positive(), Proof::Proven);
}

#[test]
fn optimized_operations_still_enforce_output_bits_and_work_limits() {
    let mut ctx = AngleContext::with_limits(Limits {
        max_bits: 3,
        ..Limits::default()
    });
    let four = ctx.from_f64(4.0).unwrap();
    assert_eq!(ctx.add(&four, &four), Err(Error::Limit(Resource::Bits)));
    let mut no_work = AngleContext::with_limits(Limits {
        max_work: 0,
        ..Limits::default()
    });
    assert_eq!(
        no_work.add(&Angle::ZERO, &Angle::ONE),
        Err(Error::Limit(Resource::Work))
    );
    assert_eq!(no_work.neg(&Angle::ONE), Err(Error::Limit(Resource::Work)));
    let mut bits = AngleContext::with_limits(Limits {
        max_bits: 1075,
        ..Limits::default()
    });
    let tiny = Angle::from_f64(f64::from_bits(1)).unwrap();
    assert_eq!(bits.half(&tiny), Err(Error::Limit(Resource::Bits)));
}

#[test]
fn parser_handles_unicode_malformed_text_and_depth_exhaustion_without_panics() {
    let mut ctx = AngleContext::default();
    for source in [
        "π",
        "sin(π)",
        "1+é",
        "\0",
        "sin(1)💥",
        "1e-9223372036854775808",
        "1e9223372036854775808",
        "1e1e1",
        "--",
        "sin()",
        "sqrt(1,2)",
    ] {
        assert!(ctx.parse(source).is_err(), "{source}");
    }
    let deeply_nested = format!("{}1{}", "(".repeat(1000), ")".repeat(1000));
    assert_eq!(
        ctx.parse(&deeply_nested),
        Err(Error::Limit(Resource::Depth))
    );
    let mut seed = 0x184f_583e_5549_98c1;
    let alphabet = ['0', '1', '+', '-', '^', '(', ')', 'e', '.', 'é', 'π', '\0'];
    for _ in 0..1000 {
        let source: String = (0..next(&mut seed) % 32)
            .map(|_| alphabet[(next(&mut seed) % alphabet.len() as u64) as usize])
            .collect();
        if let Err(Error::Parse { offset, .. }) = ctx.parse(&source) {
            assert!(source.is_char_boundary(offset));
        }
    }
}

#[test]
fn power_of_two_scaling_is_exact_across_normal_and_subnormal_boundaries() {
    let mut ctx = AngleContext::default();
    let mut seed = 0x5c74_f032_e934_2851;
    for _ in 0..150 {
        let value = f64::from_bits(next(&mut seed));
        if !value.is_finite() {
            continue;
        }
        let angle = Angle::from_f64(value).unwrap();
        for n in [i64::MIN, -8, -2, -1, 1, 2, 8] {
            for d in [1, 2, 8, 1_i64 << 62] {
                let coefficient = Rational::new(n, d).unwrap();
                let result = ctx.scale(&angle, &coefficient).unwrap();
                let expected = binary(value) * BigRational::new(n.into(), d.into());
                assert_eq!(exact(&result), expected, "{value} * {coefficient}");
                assert_eq!(
                    result.coefficient_bits(),
                    result.as_rational().unwrap().bits()
                );
            }
        }
    }
}

#[test]
fn normalized_large_coefficient_keys_agree_after_reduction_and_demotion() {
    let mut ctx = AngleContext::default();
    for bits in [32, 64, 65, 128, 512, 2048] {
        for sign in [-1_i64, 1] {
            let n: BigInt = (BigInt::one() << bits as usize) * sign + 1;
            let a = ctx.big_rational(n.clone(), BigInt::from(7)).unwrap();
            let b = ctx.big_rational(&n * 9, BigInt::from(63)).unwrap();
            assert_eq!(a, b);
            assert_eq!(hash(&a), hash(&b));
            assert_eq!(a.cmp(&b), Ordering::Equal);
            let pa = ctx.big_pi_fraction(n.clone(), BigInt::from(7)).unwrap();
            let pb = ctx.big_pi_fraction(n * 9, BigInt::from(63)).unwrap();
            assert_eq!(pa, pb);
            assert_eq!(hash(&pa), hash(&pb));
        }
    }
    let huge = ctx.parse("1e100").unwrap();
    let negative = ctx.neg(&huge).unwrap();
    assert_eq!(ctx.add(&huge, &negative).unwrap(), Angle::ZERO);
    assert_eq!(ctx.div(&huge, &huge).unwrap(), Angle::ONE);
}

#[test]
fn explicit_rotation_views_match_integer_remainders_for_large_coefficients() {
    let mut ctx = AngleContext::default();
    for exponent in [0, 10, 100, 500] {
        for n in [-17_i64, -8, -1, 0, 1, 8, 17] {
            let numerator = BigInt::from(10).pow(exponent) + n;
            let angle = ctx
                .big_pi_fraction(numerator.clone(), BigInt::from(4))
                .unwrap();
            for (period, quarters) in [(RotationPeriod::TwoPi, 8_u8), (RotationPeriod::FourPi, 16)]
            {
                let expected = numerator.mod_floor(&BigInt::from(quarters));
                assert_eq!(
                    angle.quarter_turns(period),
                    Some(expected.to_string().parse().unwrap())
                );
                let reduced = ctx.rotation_view(&angle, period).unwrap();
                assert_eq!(reduced.quarter_turns(period), angle.quarter_turns(period));
                let half_period = BigInt::from(quarters / 2);
                let centered =
                    (&numerator + &half_period).mod_floor(&BigInt::from(quarters)) - half_period;
                assert_eq!(
                    reduced,
                    ctx.big_pi_fraction(centered, BigInt::from(4)).unwrap()
                );
                assert_eq!(ctx.rotation_view(&reduced, period).unwrap(), reduced);
            }
        }
    }
    let expression = ctx.parse("sin(pi/7)+pi/3").unwrap();
    assert_eq!(
        ctx.rotation_view(&expression, RotationPeriod::TwoPi)
            .unwrap(),
        expression
    );
}

#[test]
fn cached_nodes_and_fast_cancellation_respect_input_limits() {
    let mut ctx = AngleContext::default();
    let expression = ctx.parse("sin(exp(100))").unwrap();
    let mut shallow = AngleContext::with_limits(Limits {
        max_depth: 2,
        ..Limits::default()
    });
    assert_eq!(
        shallow.add(&expression, &Angle::ZERO),
        Err(Error::Limit(Resource::Depth))
    );
    assert_eq!(
        shallow.prove_equal(&expression, &expression),
        Err(Error::Limit(Resource::Depth))
    );
    let mut few = AngleContext::with_limits(Limits {
        max_nodes: 2,
        ..Limits::default()
    });
    assert_eq!(few.sin(&expression), Err(Error::Limit(Resource::Nodes)));
    assert_eq!(few.neg(&expression), Err(Error::Limit(Resource::Nodes)));
    let mut tiny = AngleContext::with_limits(Limits {
        max_bits: 4,
        ..Limits::default()
    });
    assert_eq!(tiny.decimal("1e10"), Err(Error::Limit(Resource::Bits)));
}

#[test]
fn expression_output_obeys_exact_byte_boundary_and_retains_value_on_error() {
    let mut ctx = AngleContext::default();
    for source in ["sin(pi/7)+1/10", "exp(1/3)^4", "1e100/7", "1e-400"] {
        let value = ctx.parse(source).unwrap();
        let text = value.to_string();
        assert_eq!(value.to_expression_string(text.len()).unwrap(), text);
        assert_eq!(
            value.to_expression_string(text.len() - 1),
            Err(Error::Limit(Resource::Output))
        );
        assert_eq!(value.to_string(), text);
        assert_eq!(ctx.parse(&text).unwrap(), value);
    }
}

#[test]
fn merging_powers_preserves_valid_expressions_when_exponents_exceed_i64() {
    let mut ctx = AngleContext::default();
    let x = ctx.sin(&Angle::ONE).unwrap();
    let max = ctx.pow(&x, i64::MAX).unwrap();
    let min = ctx.pow(&x, i64::MIN).unwrap();
    let over = ctx.mul(&max, &x).unwrap();
    assert_eq!(ctx.mul(&over, &min).unwrap(), Angle::ONE);
    let under = ctx.mul(&min, &min).unwrap();
    let restored = ctx.mul(&under, &over).unwrap();
    assert_eq!(restored, min);
    let text = over.to_expression_string(1000).unwrap();
    assert_eq!(ctx.parse(&text).unwrap(), over);
}

fn same(ctx: &mut AngleContext, a: &str, b: &str) {
    let x = ctx.parse(a).unwrap();
    let y = ctx.parse(b).unwrap();
    assert_eq!(x, y, "{a} => {x} but {b} => {y}");
    assert_eq!(ctx.sub(&x, &y).unwrap(), Angle::ZERO, "{a} - {b}");
}

#[test]
fn pi_coefficients_are_extracted_from_products_and_powers() {
    let mut ctx = AngleContext::default();
    let positive = ctx.parse("2*pi*sin(1)").unwrap();
    let negated = ctx.neg(&positive).unwrap();
    assert_eq!(negated, ctx.parse("-2*pi*sin(1)").unwrap());
    same(&mut ctx, "2*pi*sin(1)", "2*sin(1)*pi");
    same(&mut ctx, "2*pi*sin(1)", "sin(1)*(pi*2)");
    same(&mut ctx, "(2*pi)*(3*pi)", "6*pi^2");
    same(&mut ctx, "(-pi)^2", "pi^2");
    same(&mut ctx, "(pi/2)^3", "pi^3/8");
    same(&mut ctx, "(pi/2)^-1*pi", "2");
}

#[test]
fn simplified_powers_are_renormalized_inside_products() {
    let mut ctx = AngleContext::default();
    same(
        &mut ctx,
        "(sqrt(sin(1))*sin(1))*sqrt(sin(1))",
        "(sqrt(sin(1))*sqrt(sin(1)))*sin(1)",
    );
    same(&mut ctx, "(sqrt(sin(1))*sin(1))*sqrt(sin(1))", "sin(1)^2");
    same(
        &mut ctx,
        "(sqrt(2*sin(1))*cos(1))*sqrt(2*sin(1))",
        "2*sin(1)*cos(1)",
    );
    same(
        &mut ctx,
        "(sqrt(sin(1)*cos(1))*exp(1))*sqrt(sin(1)*cos(1))",
        "sin(1)*cos(1)*exp(1)",
    );
    same(&mut ctx, "(sqrt(2*pi)*cos(1))*sqrt(2*pi)", "2*pi*cos(1)");
    let x = ctx.parse("(sqrt(sin(1))*sin(1))*sqrt(sin(1))").unwrap();
    assert_eq!(ctx.parse(&x.to_string()).unwrap(), x);
}

#[test]
fn even_powers_of_square_roots_simplify() {
    let mut ctx = AngleContext::default();
    same(&mut ctx, "(sqrt(2)^3)*sqrt(2)", "4");
    same(&mut ctx, "sqrt(sin(1))^4", "sin(1)^2");
    same(&mut ctx, "sqrt(sin(1))^-2", "sin(1)^-1");
    same(&mut ctx, "sqrt(sin(1))^-2*sin(1)", "1");
}

#[test]
fn pythagorean_identity_applies_in_longer_sums() {
    let mut ctx = AngleContext::default();
    same(&mut ctx, "1 + sin(1)^2 + cos(1)^2", "2");
    same(&mut ctx, "sin(1)^2 + cos(1)^2 + 1", "2");
    same(&mut ctx, "exp(1) + 3*sin(1)^2 + 3*cos(1)^2", "exp(1) + 3");
    same(
        &mut ctx,
        "sin(1)^2 + cos(2)^2 + sin(2)^2 + cos(1)^2 + pi",
        "2 + pi",
    );
    let a = ctx.parse("-1").unwrap();
    let b = ctx.parse("sin(1)^2").unwrap();
    let c = ctx.parse("cos(1)^2").unwrap();
    let mut acc = PhaseAccumulator::new(RotationPeriod::TwoPi, Limits::default());
    for x in [&a, &b, &c] {
        acc.add(&mut ctx, x).unwrap();
    }
    assert_eq!(acc.finish(&mut ctx).unwrap(), Angle::ZERO);
}

#[test]
fn displayed_expressions_reparse_at_the_depth_limit() {
    let mut ctx = AngleContext::default();
    let mut source = "1".to_string();
    for _ in 0..16 {
        source = format!("sin(1+{source})");
    }
    let mut a = ctx.parse(&source).unwrap();
    assert_eq!(ctx.parse(&a.to_string()).unwrap(), a);
    // Grow to the scalar depth limit, mixing every displayed node kind.
    loop {
        let next = ctx
            .parse(&format!("-1/3*pi + 2/3*sin(-5/7 + (2*pi/3)*exp({a})^-3)"))
            .or_else(|_| ctx.parse(&format!("sin(1/2+{a})")));
        match next {
            Ok(next) => a = next,
            Err(_) => break,
        }
    }
    assert!(a.expression_depth() >= ctx.limits().max_depth - 1);
    assert_eq!(ctx.parse(&a.to_string()).unwrap(), a);
}

#[test]
fn sum_factors_have_one_representation_however_scaled() {
    let mut ctx = AngleContext::default();
    let a = ctx.parse("(sin(1)+cos(1))*(2*sin(2))").unwrap();
    let printed = a.to_string();
    let reparsed = ctx.parse(&printed).unwrap();
    assert_eq!(reparsed, a, "{printed} reparsed as {reparsed}");
    assert_eq!(ctx.prove_equal(&a, &reparsed).unwrap(), Proof::Proven);
    same(&mut ctx, "(sin(1)+cos(1))*(2*sin(2))", "2*(sin(1)+cos(1))*sin(2)");
    same(&mut ctx, "(sin(1)+cos(1))*(2*sin(2))", "(2*sin(1)+2*cos(1))*sin(2)");
    same(&mut ctx, "(-sin(1)-cos(1))*sin(2)", "-((sin(1)+cos(1))*sin(2))");
    same(&mut ctx, "(2*sin(1)+2*cos(1))^2", "4*(sin(1)+cos(1))^2");
    same(&mut ctx, "(2*sin(1)+2*cos(1))*(sin(1)+cos(1))", "2*(sin(1)+cos(1))^2");
    same(&mut ctx, "(sin(1)/3+cos(1)/3)^-1*(sin(1)+cos(1))", "3");
    for source in [
        "(sin(1)+cos(1))*(2*sin(2))",
        "(2*sin(1)+2*cos(1))^3*exp(1)",
        "-(pi/3+sin(1))*(5/7*cos(2)+exp(1))",
    ] {
        let x = ctx.parse(source).unwrap();
        assert_eq!(ctx.parse(&x.to_string()).unwrap(), x, "{source}");
    }
}

#[test]
fn decimal_size_limits_apply_to_the_exact_coefficient() {
    let mut ctx = AngleContext::default();
    let big = ctx.parse("(2^30000)^2").unwrap();
    assert_eq!(big.coefficient_bits(), 60_001);
    // 18,062 digits: four bits per digit would overestimate it past the limit.
    assert_eq!(ctx.parse(&big.to_string()).unwrap(), big);
    // The limit applies to the reduced value, not the digits spelling it.
    let one = format!("1.{}", "0".repeat(20_000));
    assert_eq!(ctx.parse(&one).unwrap(), Angle::ONE);
    assert_eq!(ctx.parse(&format!("{one}e0")).unwrap(), Angle::ONE);
    // Provably oversized magnitudes are still rejected before allocating.
    for source in ["1e20000", "1e-20000", "1e999999999", "5e-999999999"] {
        assert_eq!(ctx.parse(source), Err(Error::Limit(Resource::Bits)), "{source}");
    }
    let over = ctx.parse("2^65536");
    assert_eq!(over, Err(Error::Limit(Resource::Bits)));
}
