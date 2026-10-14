use super::*;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use std::{
    collections::{BTreeSet, HashSet},
    hash::{Hash, Hasher},
};

fn integer(n: i64) -> Angle {
    Angle::rational(n, 1).unwrap()
}
fn rational(a: &Angle) -> BigRational {
    a.as_rational().unwrap().to_big_rational()
}
fn independent_binary(value: f64) -> BigRational {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32;
    let mantissa = bits & ((1 << 52) - 1);
    let mut n = BigInt::from(mantissa | if exponent == 0 { 0 } else { 1 << 52 });
    if bits >> 63 != 0 {
        n = -n;
    }
    let power = if exponent == 0 {
        -1074
    } else {
        exponent - 1075
    };
    if power >= 0 {
        BigRational::from_integer(n << power as usize)
    } else {
        BigRational::new(n, BigInt::one() << -power as usize)
    }
}

#[test]
fn scalar_values_are_nonperiodic_and_rotation_views_are_explicit() {
    let mut ctx = AngleContext::default();
    let two = ctx.pi_fraction(2, 1).unwrap();
    let four = ctx.pi_fraction(4, 1).unwrap();
    assert_ne!(Angle::ZERO, two);
    assert_ne!(two, four);
    assert_eq!(two.prove_zero(), Proof::Disproven);
    assert_eq!(two.prove_identity(RotationPeriod::TwoPi), Proof::Proven);
    assert_eq!(two.prove_identity(RotationPeriod::FourPi), Proof::Disproven);
    assert_eq!(four.prove_identity(RotationPeriod::FourPi), Proof::Proven);
    assert_eq!(
        ctx.rotation_view(&two, RotationPeriod::TwoPi).unwrap(),
        Angle::ZERO
    );
    assert_eq!(ctx.add(&two, &two).unwrap(), four);
    assert_eq!(ctx.half(&two).unwrap(), Angle::PI);
    assert_eq!(ctx.neg(&two).unwrap(), ctx.pi_fraction(-2, 1).unwrap());
    let set: HashSet<_> = [Angle::ZERO, two, four].into_iter().collect();
    assert_eq!(set.len(), 3);
}

#[test]
fn exact_decimal_input_precedence_and_large_coefficients() {
    let mut ctx = AngleContext::default();
    for (source, n, d) in [
        ("0.1+0.2", 3, 10),
        ("(1.0e16+1)-1.0e16", 1, 1),
        ("-2^2", -4, 1),
        ("(-2)^2", 4, 1),
        ("2^3^2", 512, 1),
        ("2^-3", 1, 8),
        (".125", 1, 8),
        ("1.", 1, 1),
        ("+1e+2", 100, 1),
        ("6/-8", -3, 4),
        ("-9223372036854775808/-9223372036854775808", 1, 1),
    ] {
        assert_eq!(
            ctx.parse(source).unwrap(),
            ctx.rational(n, d).unwrap(),
            "{source}"
        );
    }
    let big = ctx.parse("1e400/3").unwrap();
    assert_eq!(
        rational(&big),
        BigRational::new(BigInt::from(10).pow(400), 3.into())
    );
    let small = ctx.parse("1e-400").unwrap();
    assert_eq!(
        rational(&small),
        BigRational::new(1.into(), BigInt::from(10).pow(400))
    );
    assert!(!big.is_inline());
    assert_eq!(
        ctx.parse("(pi/9223372036854775808)*9223372036854775808")
            .unwrap(),
        Angle::PI
    );
}

#[test]
fn rational_promotion_and_demotion_match_reference_arithmetic() {
    let mut ctx = AngleContext::default();
    let max = ctx.rational(i64::MAX, 1).unwrap();
    let twice = ctx.add(&max, &max).unwrap();
    assert_eq!(
        rational(&twice),
        BigRational::from_integer(BigInt::from(i64::MAX) * 2)
    );
    assert!(twice.as_rational().unwrap().as_small().is_none());
    assert_eq!(ctx.sub(&twice, &max).unwrap(), max);
    assert_eq!(ctx.div(&twice, &max).unwrap(), integer(2));
    let n = ctx.rational(i64::MIN, 1).unwrap();
    assert_eq!(
        rational(&ctx.neg(&n).unwrap()),
        BigRational::from_integer(-BigInt::from(i64::MIN))
    );
    let huge = ctx.parse("1e200").unwrap();
    assert_eq!(ctx.sub(&huge, &huge).unwrap(), Angle::ZERO);
    assert_eq!(ctx.div(&huge, &huge).unwrap(), Angle::ONE);
    assert!(ctx.sub(&huge, &huge).unwrap().is_inline());
}

#[test]
fn arithmetic_fuzz_matches_independent_big_rationals() {
    let mut ctx = AngleContext::default();
    let mut seed = 0x38bc_529f_7102_7631u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as i64
    };
    for _ in 0..3000 {
        let (a, b, c, d) = (next(), next() | 1, next(), next() | 1);
        let aa = ctx.rational(a, b).unwrap();
        let bb = ctx.rational(c, d).unwrap();
        let a = BigRational::new(a.into(), b.into());
        let b = BigRational::new(c.into(), d.into());
        assert_eq!(rational(&ctx.add(&aa, &bb).unwrap()), &a + &b);
        assert_eq!(rational(&ctx.sub(&aa, &bb).unwrap()), &a - &b);
        assert_eq!(rational(&ctx.mul(&aa, &bb).unwrap()), &a * &b);
        if !b.is_zero() {
            assert_eq!(rational(&ctx.div(&aa, &bb).unwrap()), &a / &b);
        }
    }
}

#[test]
fn binary64_inputs_and_nonrepresentable_sums_remain_exact() {
    let mut ctx = AngleContext::default();
    for value in [
        0.0,
        -0.0,
        0.1,
        0.2,
        1e16,
        f64::MAX,
        -f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
    ] {
        let a = ctx.from_f64(value).unwrap();
        assert!(a.is_inline());
        assert_eq!(rational(&a), independent_binary(value));
        let half = ctx.half(&a).unwrap();
        assert_eq!(rational(&half), independent_binary(value) / BigInt::from(2));
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(Angle::from_f64(value), Err(Error::NonFinite));
    }
    let a = ctx.from_f64(1e16).unwrap();
    let b = ctx.from_f64(1.0).unwrap();
    let sum = ctx.add(&a, &b).unwrap();
    assert_eq!(
        rational(&sum),
        independent_binary(1e16) + BigRational::one()
    );
    assert_eq!(ctx.sub(&sum, &a).unwrap(), Angle::ONE);
    let max = ctx.from_f64(f64::MAX).unwrap();
    assert_eq!(
        rational(&ctx.add(&max, &max).unwrap()),
        independent_binary(f64::MAX) * BigInt::from(2)
    );
    assert_ne!(ctx.parse("0.1").unwrap(), ctx.from_f64(0.1).unwrap());
}

#[test]
fn binary64_random_decoding_matches_reference() {
    let mut seed = 0x68bf_f495_55f0_6a3bu64;
    for _ in 0..5000 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        let value = f64::from_bits(seed);
        if value.is_finite() {
            let a = Angle::from_f64(value).unwrap();
            assert_eq!(rational(&a), independent_binary(value));
            assert_eq!(a.coefficient_bits(), a.as_rational().unwrap().bits());
        }
    }
}

#[test]
fn equal_numeric_representations_have_equal_hashes_and_order() {
    let mut ctx = AngleContext::default();
    for value in [0.0, 0.5, 0.125, 42.0, 0.1, f64::from_bits(1), f64::MAX] {
        let dyadic = ctx.from_f64(value).unwrap();
        let r = independent_binary(value);
        let exact = ctx
            .big_rational(r.numer().clone(), r.denom().clone())
            .unwrap();
        assert_eq!(dyadic, exact);
        assert_eq!(dyadic.cmp(&exact), Ordering::Equal);
        let hash = |a: &Angle| {
            let mut h = FxHasher::default();
            a.hash(&mut h);
            h.finish()
        };
        assert_eq!(hash(&dyadic), hash(&exact));
        let set: HashSet<_> = [dyadic.clone(), exact.clone()].into_iter().collect();
        assert_eq!(set.len(), 1);
        let set: BTreeSet<_> = [dyadic, exact].into_iter().collect();
        assert_eq!(set.len(), 1);
    }
}

#[test]
fn symbolic_linear_terms_cancel_without_evaluation() {
    let mut ctx = AngleContext::default();
    let a = ctx.parse("sin(pi/7)+1/10").unwrap();
    let b = ctx.parse("-sin(pi/7)+1/5").unwrap();
    assert_eq!(ctx.add(&a, &b).unwrap(), ctx.rational(3, 10).unwrap());
    let x = ctx.parse("sin(pi/7)+cos(1/3)").unwrap();
    let half = ctx.half(&x).unwrap();
    assert_eq!(ctx.add(&half, &half).unwrap(), x);
    assert_eq!(ctx.sub(&x, &x).unwrap(), Angle::ZERO);
    let neg = ctx.neg(&x).unwrap();
    assert_eq!(ctx.add(&x, &neg).unwrap(), Angle::ZERO);
    let product = ctx
        .parse("2*sin(pi/7)*cos(1/3)+3*sin(pi/7)*cos(1/3)")
        .unwrap();
    assert_eq!(product, ctx.parse("5*sin(pi/7)*cos(1/3)").unwrap());
    assert_eq!(
        ctx.parse("(2*sin(pi/7))*(3*cos(1/3))").unwrap(),
        ctx.parse("6*sin(pi/7)*cos(1/3)").unwrap()
    );
}

#[test]
fn special_functions_and_domain_validated_identities() {
    let mut ctx = AngleContext::default();
    for (source, expected) in [
        ("sin(pi)", "0"),
        ("cos(pi)", "-1"),
        ("tan(pi/4)", "1"),
        ("sin(pi/4)^2", "1/2"),
        ("sqrt(4/9)", "2/3"),
        ("sqrt(2)^2", "2"),
        ("sqrt((-pi)^2)", "pi"),
        ("ln(exp(sin(1)))", "sin(1)"),
        ("exp(ln(2))", "2"),
        ("sin(1/3)^2+cos(1/3)^2", "1"),
        ("sin(-pi/7)", "-sin(pi/7)"),
        ("cos(-pi/7)", "cos(pi/7)"),
        ("sin(15*pi/7)", "sin(pi/7)"),
        ("tan(8*pi/7)", "tan(pi/7)"),
    ] {
        assert_eq!(
            ctx.parse(source).unwrap(),
            ctx.parse(expected).unwrap(),
            "{source}"
        );
    }
    for source in [
        "sin(pi/7)",
        "cos(1/3)",
        "tan(1/3)",
        "exp(1/5)",
        "ln(2)",
        "sqrt(2)",
        "ln(sin(1))",
    ] {
        assert!(
            ctx.parse(source).unwrap().expression_nodes() > 1,
            "{source}"
        );
    }
}

#[test]
fn invalid_or_unknown_domains_are_never_cancelled_away() {
    let mut ctx = AngleContext::default();
    for source in [
        "1/0",
        "0*(1/0)",
        "(1/0)-(1/0)",
        "sqrt(-1)-sqrt(-1)",
        "ln(0)",
        "ln(-2)",
        "tan(pi/2)",
        "0^0",
        "0^-1",
    ] {
        assert!(ctx.parse(source).is_err(), "{source}");
    }
    assert_eq!(
        ctx.parse("sqrt(-1)"),
        Err(Error::InvalidDomain(Function::Sqrt))
    );
    assert_eq!(
        ctx.parse("tan(pi/2)"),
        Err(Error::InvalidDomain(Function::Tan))
    );
    assert_eq!(ctx.parse("1/sin(exp(10))"), Err(Error::DomainUnproved));
    assert_eq!(ctx.parse("0/sin(exp(10))"), Err(Error::DomainUnproved));
    assert_eq!(ctx.parse("sqrt(sin(exp(10)))"), Err(Error::DomainUnproved));
    let unknown = ctx.parse("sin(exp(10))").unwrap();
    assert_eq!(unknown.prove_zero(), Proof::Unknown);
    assert_eq!(unknown.prove_positive(), Proof::Unknown);
}

#[test]
fn interning_shares_nodes_and_values_survive_context_destruction() {
    let mut ctx = AngleContext::default();
    let a = ctx.parse("sin(pi/7)+cos(1/3)").unwrap();
    let b = ctx.parse("cos(1/3)+sin(pi/7)").unwrap();
    assert!(a.shares_storage(&b));
    assert!(ctx.stats().hits > 0);
    let before = ctx.stats().retained_entries;
    drop(b);
    assert_eq!(ctx.prune(), 0);
    let result = a.to_expression_string(1000).unwrap();
    drop(ctx);
    assert_eq!(a.to_expression_string(1000).unwrap(), result);
    let mut another = AngleContext::default();
    let c = another.parse("cos(1/3)+sin(pi/7)").unwrap();
    assert_eq!(a, c);
    assert!(!a.shares_storage(&c));
    assert_eq!(another.prove_equal(&a, &c).unwrap(), Proof::Proven);
    drop(a);
    drop(c);
    assert!(another.prune() > 0);
    assert!(before > 0);
}

#[test]
fn collision_buckets_and_equal_fingerprints_require_full_keys() {
    let mut ctx = AngleContext::default();
    let a = ctx.sin(&integer(1)).unwrap();
    let b = ctx.cos(&integer(1)).unwrap();
    let (Value::Shared(na), Value::Shared(nb)) = (&a.0, &b.0) else {
        panic!()
    };
    let forged = Angle(Value::Shared(Arc::new(Node {
        kind: nb.kind.clone(),
        hash: na.hash,
        meta: nb.meta,
    })));
    assert_ne!(a, forged);
    assert_ne!(a.cmp(&forged), Ordering::Equal);
    ctx.clear_interner();
    ctx.nodes.insert(nb.hash, vec![Arc::downgrade(na)]);
    ctx.stats.retained_entries = 1;
    let rebuilt = ctx.cos(&integer(1)).unwrap();
    assert_eq!(rebuilt, b);
    assert!(!rebuilt.shares_storage(&a));
}

#[test]
fn exact_expression_printing_roundtrips_across_contexts() {
    let mut ctx = AngleContext::default();
    for source in [
        "pi/7+3/10",
        "sin(pi/7)+cos(1/3)",
        "(sqrt(2)+exp(1/5))/11",
        "2*sin(pi/7)*cos(1/3)",
        "1/pi",
        "sin(1)^-2",
        "sqrt(sin(exp(10))^2)",
        "1e-400",
        "1e400/3",
    ] {
        let a = ctx.parse(source).unwrap();
        let printed = a.to_expression_string(10000).unwrap();
        let b: Angle = printed.parse().unwrap();
        assert_eq!(a, b, "{source}: {printed}");
        assert_eq!(a.to_string(), b.to_string());
    }
    let a = ctx.parse("sin(pi/7)+cos(1/3)").unwrap();
    assert_eq!(
        a.to_expression_string(3),
        Err(Error::Limit(Resource::Output))
    );
}

#[test]
fn non_integer_powers_use_the_real_principal_value() {
    let mut ctx = AngleContext::default();
    for (power, value) in [
        ("2^(1/2)", "sqrt(2)"),
        ("2^0.5", "sqrt(2)"),
        ("4^(1/2)", "2"),
        ("2^(1/3)", "exp(ln(2)/3)"),
        ("2^-0.5", "exp(-ln(2)/2)"),
        ("pi^pi", "exp(pi*ln(pi))"),
        ("0^(1/2)", "0"),
        ("-2^0.5", "-sqrt(2)"),
        ("2^3", "8"),
    ] {
        let power_angle = ctx.parse(power).unwrap();
        let value_angle = ctx.parse(value).unwrap();
        assert_eq!(
            ctx.prove_equal(&power_angle, &value_angle).unwrap(),
            Proof::Proven,
            "{power}"
        );
    }
    for source in ["(-2)^0.5", "(-1)^(1/3)", "0^(-1/2)"] {
        assert_eq!(ctx.parse(source), Err(Error::NonRealPower), "{source}");
    }
}

#[test]
fn malformed_input_and_all_budget_limits_fail_explicitly() {
    let mut ctx = AngleContext::default();
    for source in [
        "", ".", "1e", "1e+", "1..2", "pi2", "sin 1", "foo(1)", "1+", "(1", "sin(1,2)",
    ] {
        assert!(ctx.parse(source).is_err(), "{source}");
    }
    let defaults = Limits::default();
    for (limits, source, resource) in [
        (
            Limits {
                max_bits: 16,
                ..defaults
            },
            "1e100",
            Resource::Bits,
        ),
        (
            Limits {
                max_depth: 4,
                ..defaults
            },
            "sin(sin(sin(sin(1))))",
            Resource::Depth,
        ),
        (
            Limits {
                max_terms: 1,
                ..defaults
            },
            "sin(1)+cos(1)",
            Resource::Terms,
        ),
        (
            Limits {
                max_work: 2,
                ..defaults
            },
            "1+2",
            Resource::Work,
        ),
        (
            Limits {
                max_nodes: 1,
                ..defaults
            },
            "sin(1)",
            Resource::Nodes,
        ),
        (
            Limits {
                max_source_bytes: 3,
                ..defaults
            },
            "12345",
            Resource::Source,
        ),
    ] {
        assert_eq!(
            AngleContext::with_limits(limits).parse(source),
            Err(Error::Limit(resource)),
            "{source}"
        );
    }
    assert_eq!(ctx.parse("1e1000000000"), Err(Error::Limit(Resource::Bits)));
    let mut tiny = AngleContext::with_limits(Limits {
        max_bits: 8,
        ..defaults
    });
    let huge = ctx.parse("1e100").unwrap();
    assert_eq!(
        tiny.add(&huge, &Angle::ZERO),
        Err(Error::Limit(Resource::Bits))
    );
    assert_eq!(
        rational(&huge),
        BigRational::from_integer(BigInt::from(10).pow(100))
    );
}

#[test]
fn independent_workers_and_reversed_construction_have_same_output() {
    let a = std::thread::spawn(|| AngleContext::default().parse("sin(pi/7)+cos(1/3)").unwrap());
    let b = std::thread::spawn(|| AngleContext::default().parse("cos(1/3)+sin(pi/7)").unwrap());
    let a = a.join().unwrap();
    let b = b.join().unwrap();
    assert_eq!(a, b);
    assert_eq!(a.to_string(), b.to_string());
    let mut ctx = AngleContext::default();
    assert_eq!(ctx.sub(&a, &b).unwrap(), Angle::ZERO);
}

#[test]
fn inline_scalar_size_matches_the_existing_angle_size() {
    assert_eq!(
        std::mem::size_of::<Angle>(),
        std::mem::size_of::<crate::angle::Angle>()
    );
    assert!(Angle::pi_fraction(1, 7).unwrap().is_inline());
    assert!(Angle::rational(3, 10).unwrap().is_inline());
}

#[test]
fn weak_interner_churn_prunes_the_lookup_bucket() {
    let mut ctx = AngleContext::default();
    let arg = ctx.rational(1, 3).unwrap();
    for _ in 0..10_000 {
        let value = ctx.sin(&arg).unwrap();
        assert!(!value.is_inline());
    }
    assert_eq!(ctx.stats().retained_entries, 1);
    assert_eq!(ctx.stats().pruned, 9_999);
    assert_eq!(ctx.prune(), 1);
    assert_eq!(ctx.stats().retained_entries, 0);
}

#[test]
fn exponent_metadata_obeys_coefficient_limits() {
    let mut ctx = AngleContext::with_limits(Limits {
        max_bits: 8,
        ..Limits::default()
    });
    let arg = ctx.sin(&Angle::ONE).unwrap();
    assert_eq!(ctx.pow(&arg, 256), Err(Error::Limit(Resource::Bits)));
    let power = ctx.pow(&arg, 127).unwrap();
    assert_eq!(power.coefficient_bits(), 7);
}

#[test]
fn product_permutations_and_expression_roundtrips_agree() {
    let mut ctx = AngleContext::default();
    let atoms = ["sin(pi/7)", "cos(1/3)", "exp(1)"];
    for a in atoms {
        for b in atoms {
            for c in atoms {
                let first = ctx.parse(&format!("2*{a}*{b}+3*{b}*{a}+{c}")).unwrap();
                let second = ctx.parse(&format!("{c}+5*{b}*{a}")).unwrap();
                assert_eq!(first, second);
                let printed = first.to_expression_string(10_000).unwrap();
                assert_eq!(first, ctx.parse(&printed).unwrap(), "{printed}");
            }
        }
    }
}

#[test]
fn inline_addition_fast_paths_match_reference_at_integer_boundaries() {
    let mut ctx = AngleContext::default();
    for a in [i32::MIN as i64, -7, -1, 0, 1, 7, i32::MAX as i64] {
        for b in [i32::MIN as i64, -7, -1, 0, 1, 7, i32::MAX as i64] {
            for d in [1, 2, 4, 7, 1 << 31] {
                for e in [1, 2, 4, 7, 1 << 31] {
                    let expected =
                        BigRational::new(a.into(), d.into()) + BigRational::new(b.into(), e.into());
                    let x = ctx.rational(a, d).unwrap();
                    let y = ctx.rational(b, e).unwrap();
                    assert_eq!(rational(&ctx.add(&x, &y).unwrap()), expected);
                    let x = ctx.pi_fraction(a, d).unwrap();
                    let y = ctx.pi_fraction(b, e).unwrap();
                    let sum = ctx.add(&x, &y).unwrap();
                    assert_eq!(
                        sum.pi_coefficient()
                            .unwrap_or(Rational::ZERO)
                            .to_big_rational(),
                        expected
                    );
                }
            }
        }
    }
}
