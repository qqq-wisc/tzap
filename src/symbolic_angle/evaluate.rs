use super::{Angle, Error, Function, view::ExpressionRef};
use num_traits::ToPrimitive;
pub(crate) fn lossy(a: &Angle) -> Result<f64, Error> {
    let number = |r: &super::Rational| r.to_big_rational().to_f64().ok_or(Error::NonFinite);
    let value = match a.expression() {
        ExpressionRef::Rational(r) => number(&r)?,
        ExpressionRef::Pi(r) => number(&r)? * std::f64::consts::PI,
        ExpressionRef::Sum(terms) => {
            let mut sum = 0.0;
            for (c, a) in terms {
                sum += number(c)? * lossy(a)?;
            }
            sum
        }
        ExpressionRef::Product {
            coefficient,
            factors,
        } => {
            let mut p = number(coefficient)?;
            for a in factors {
                p *= lossy(a)?;
            }
            p
        }
        ExpressionRef::Power(a, n) => lossy(a)?.powf(n as f64),
        ExpressionRef::Function(fun, a) => {
            let x = lossy(a)?;
            match fun {
                Function::Sin => x.sin(),
                Function::Cos => x.cos(),
                Function::Tan => x.tan(),
                Function::Exp => x.exp(),
                Function::Ln => x.ln(),
                Function::Sqrt => x.sqrt(),
            }
        }
    };
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::NonFinite)
    }
}

// Enclosing arithmetic is shared with the independent synthesis verifier.
use crate::synthesis::interval::{Context, Interval};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Zero;
use rustc_hash::FxHashMap;

struct Evaluator<'a> {
    ctx: &'a Context,
    pi: Interval,
    cache: FxHashMap<Angle, Interval>,
    work: usize,
}
pub(crate) fn enclose(a: &Angle, ctx: &Context) -> Result<Interval, &'static str> {
    if ctx.bits > 2048 || a.expression_nodes() > 4096 || a.expression_depth() > 64 {
        return Err("symbolic evaluation resource limit");
    }
    Evaluator {
        ctx,
        pi: ctx.pi(),
        cache: FxHashMap::default(),
        work: 4_000_000,
    }
    .value(a)
}
impl Evaluator<'_> {
    fn charge(&mut self, n: usize) -> Result<(), &'static str> {
        self.work = self
            .work
            .checked_sub(n)
            .ok_or("symbolic evaluation work limit")?;
        Ok(())
    }
    fn bounded(&self, x: Interval) -> Result<Interval, &'static str> {
        if x.lo.bits().max(x.hi.bits()) > self.ctx.bits as u64 + 8192 {
            Err("symbolic evaluation magnitude limit")
        } else {
            Ok(x)
        }
    }
    fn rational(&self, r: &super::Rational) -> Interval {
        let r = r.to_big_rational();
        let n = r.numer() * &self.ctx.unit;
        Interval {
            lo: n.div_floor(r.denom()),
            hi: n.div_ceil(r.denom()),
        }
    }
    fn mul(&mut self, a: &Interval, b: &Interval) -> Result<Interval, &'static str> {
        let limbs = a.abs_upper().bits().max(b.abs_upper().bits()).div_ceil(64) as usize;
        self.charge(limbs.saturating_mul(limbs).max(1))?;
        self.bounded(self.ctx.mul(a, b))
    }
    fn div(&mut self, a: &Interval, b: &Interval) -> Result<Interval, &'static str> {
        if b.lo <= BigInt::zero() && b.hi >= BigInt::zero() {
            return Err("denominator enclosure contains zero");
        }
        self.charge(1)?;
        let mut lo = None;
        let mut hi = None;
        for n in [&a.lo, &a.hi] {
            for d in [&b.lo, &b.hi] {
                let n = n * &self.ctx.unit;
                let lower = n.div_floor(d);
                let upper = n.div_ceil(d);
                lo = Some(lo.map_or(lower.clone(), |v: BigInt| v.min(lower)));
                hi = Some(hi.map_or(upper.clone(), |v: BigInt| v.max(upper)));
            }
        }
        self.bounded(Interval {
            lo: lo.unwrap(),
            hi: hi.unwrap(),
        })
    }
    fn power(&mut self, x: &Interval, n: i64) -> Result<Interval, &'static str> {
        let mut result = self.ctx.integer(1);
        let mut base = x.clone();
        let mut e = n.unsigned_abs();
        while e != 0 {
            if e & 1 != 0 {
                result = self.mul(&result, &base)?;
            }
            e >>= 1;
            if e != 0 {
                base = self.mul(&base, &base)?;
            }
        }
        if n < 0 {
            self.div(&self.ctx.integer(1), &result)
        } else {
            Ok(result)
        }
    }
    fn trig(&mut self, x: &Interval) -> Result<(Interval, Interval), &'static str> {
        let tau = self.pi.scale(&2.into());
        let turns = (x.middle() + tau.middle() / 2u32).div_floor(&tau.middle());
        let reduced = x.sub(&tau.scale(&turns));
        self.charge((self.ctx.bits + 32).saturating_mul(self.ctx.bits.div_ceil(64).pow(2)))?;
        self.ctx.sin_cos(&reduced)
    }
    fn exp(&mut self, x: &Interval) -> Result<Interval, &'static str> {
        if x.abs_upper() > &self.ctx.unit * 1024 {
            return Err("exponential magnitude limit");
        }
        let mut reduced = x.clone();
        let mut squarings = 0;
        while reduced.abs_upper() * 2 > self.ctx.unit {
            reduced = reduced.div_int(2);
            squarings += 1;
        }
        let mut sum = self.ctx.integer(1);
        let mut term = sum.clone();
        for n in 1..self.ctx.bits + 128 {
            term = self.mul(&term, &reduced)?.div_int(n as u64);
            if n > 8 && term.abs_upper() <= BigInt::from(2) {
                // On |r|<=1/2, the tail beginning at this omitted term is
                // bounded by twice its absolute value (a geometric majorant).
                sum = sum.widen(&(term.abs_upper() * 2));
                for _ in 0..squarings {
                    sum = self.mul(&sum, &sum)?;
                }
                return self.bounded(sum);
            }
            sum = sum.add(&term);
        }
        Err("exponential precision limit")
    }
    fn ln_kernel(&mut self, x: &Interval) -> Result<Interval, &'static str> {
        let one = self.ctx.integer(1);
        let z = self.div(&x.sub(&one), &x.add(&one))?;
        if z.abs_upper() * 2 > self.ctx.unit {
            return Err("logarithm reduction unresolved");
        }
        let square = self.mul(&z, &z)?;
        let mut power = z.clone();
        let mut sum = self.ctx.integer(0);
        for n in 0..self.ctx.bits + 32 {
            let term = power.div_int((2 * n + 1) as u64);
            if n > 8 && term.abs_upper() <= BigInt::from(2) {
                // 2 * tail <= 2*|z|^(2n+1)/(2n+1)/(1-z^2) <= 3*term.
                return Ok(sum.scale(&2.into()).widen(&(term.abs_upper() * 3)));
            }
            sum = sum.add(&term);
            power = self.mul(&power, &square)?;
        }
        Err("logarithm precision limit")
    }
    fn ln_point(&mut self, n: &BigInt) -> Result<Interval, &'static str> {
        let k = n.bits() as i64 - self.ctx.bits as i64 - 1;
        if k.unsigned_abs() > 8192 {
            return Err("logarithm magnitude limit");
        }
        let x = if k >= 0 {
            let d = BigInt::from(1) << k as usize;
            Interval {
                lo: n.div_floor(&d),
                hi: n.div_ceil(&d),
            }
        } else {
            Interval::point(n << (-k) as usize)
        };
        let small = self.ln_kernel(&x)?;
        let two = self.ctx.integer(2);
        let log_two = self.ln_kernel(&two)?;
        Ok(small.add(&log_two.scale(&k.into())))
    }
    fn value(&mut self, a: &Angle) -> Result<Interval, &'static str> {
        self.charge(1)?;
        if let Some(x) = self.cache.get(a) {
            return Ok(x.clone());
        }
        if self.cache.len() >= 4096 {
            return Err("symbolic evaluation node limit");
        }
        let x = match a.expression() {
            ExpressionRef::Rational(r) => self.rational(&r),
            ExpressionRef::Pi(r) => {
                let r = self.rational(&r);
                self.mul(&r, &self.pi.clone())?
            }
            ExpressionRef::Sum(terms) => {
                let mut sum = self.ctx.integer(0);
                for (c, a) in terms {
                    let value = self.value(a)?;
                    let c = self.rational(c);
                    sum = sum.add(&self.mul(&c, &value)?);
                }
                sum
            }
            ExpressionRef::Product {
                coefficient,
                factors,
            } => {
                let mut product = self.rational(coefficient);
                for a in factors {
                    let value = self.value(a)?;
                    product = self.mul(&product, &value)?;
                }
                product
            }
            ExpressionRef::Power(a, n) => {
                let a = self.value(a)?;
                self.power(&a, n)?
            }
            ExpressionRef::Function(fun, a) => {
                let x = self.value(a)?;
                match fun {
                    Function::Sin => self.trig(&x)?.0,
                    Function::Cos => self.trig(&x)?.1,
                    Function::Tan => {
                        let (s, c) = self.trig(&x)?;
                        self.div(&s, &c)?
                    }
                    Function::Exp => self.exp(&x)?,
                    Function::Ln => {
                        if x.lo <= BigInt::zero() {
                            return Err("logarithm domain unresolved");
                        }
                        let lo = self.ln_point(&x.lo)?;
                        let hi = self.ln_point(&x.hi)?;
                        Interval {
                            lo: lo.lo,
                            hi: hi.hi,
                        }
                    }
                    Function::Sqrt => {
                        if x.lo < BigInt::zero() {
                            return Err("square root domain unresolved");
                        }
                        let lo = (&x.lo * &self.ctx.unit).sqrt();
                        let h = &x.hi * &self.ctx.unit;
                        let root = h.sqrt();
                        let hi = if &root * &root == h { root } else { root + 1 };
                        Interval { lo, hi }
                    }
                }
            }
        };
        let x = self.bounded(x)?;
        self.cache.insert(a.clone(), x.clone());
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::BigRational as R;
    use num_traits::One;
    fn q(n: i64, d: i64) -> R {
        R::new(n.into(), d.into())
    }
    fn endpoint(n: &BigInt, ctx: &Context) -> R {
        R::new(n.clone(), ctx.unit.clone())
    }
    fn check(source: &str, lower: R, upper: R) {
        check_at_precision(source, lower, upper, 192);
    }
    fn check_at_precision(source: &str, lower: R, upper: R, bits: usize) {
        let mut ctx = super::super::AngleContext::default();
        let a = ctx.parse(source).unwrap();
        let precision = Context::new(bits);
        let result = enclose(&a, &precision).unwrap();
        assert!(
            endpoint(&result.lo, &precision) <= lower,
            "lower endpoint excludes {source}"
        );
        assert!(
            endpoint(&result.hi, &precision) >= upper,
            "upper endpoint excludes {source}"
        );
        assert!(
            endpoint(&(&result.hi - &result.lo), &precision) < q(1, 1_000_000_000),
            "enclosure too loose: {source}"
        );
    }
    // Independent exact rational series with analytic remainder bounds.
    fn trig(x: R, sine: bool) -> (R, R) {
        let mut term = if sine { x.clone() } else { R::one() };
        let mut sum = term.clone();
        for n in 0..100 {
            let k = if sine { 2 * n + 2 } else { 2 * n + 1 };
            term = -term * &x * &x / q(k * (k + 1), 1);
            sum += &term;
        }
        let k = if sine { 202 } else { 201 };
        let next = -term * &x * &x / q(k * (k + 1), 1);
        let other = &sum + &next;
        (sum.clone().min(other.clone()), sum.max(other))
    }
    fn exponential(x: R) -> (R, R) {
        let mut term = R::one();
        let mut sum = term.clone();
        for n in 1..=140 {
            term = term * &x / q(n, 1);
            sum += &term;
        }
        let tail = term * &x / q(141, 1) * q(2, 1);
        (sum.clone(), sum + tail)
    }
    fn logarithm(x: R) -> (R, R) {
        let z = (&x - q(1, 1)) / (&x + q(1, 1));
        let square = &z * &z;
        let mut power = z;
        let mut sum = R::from_integer(0.into());
        for n in 0..140 {
            sum += &power / q(2 * n + 1, 1);
            power *= &square;
        }
        let lower = sum * q(2, 1);
        let upper = &lower + power / q(281, 1) * q(3, 1);
        (lower, upper)
    }
    #[test]
    fn every_function_encloses_an_independent_rational_reference() {
        let (slo, shi) = trig(q(1, 1), true);
        let (clo, chi) = trig(q(1, 1), false);
        check("sin(1)", slo.clone(), shi.clone());
        check("cos(1)", clo.clone(), chi.clone());
        check("tan(1)", &slo / &chi, &shi / &clo);
        let (lo, hi) = exponential(q(1, 1));
        check("exp(1)", lo.clone(), hi.clone());
        check("exp(-1)", q(1, 1) / hi, q(1, 1) / lo);
        for (source, x) in [("exp(1/2)", q(1, 2)), ("exp(3/2)", q(3, 2))] {
            let (lo, hi) = exponential(x);
            check(source, lo, hi);
        }
        for (source, x) in [("ln(2)", q(2, 1)), ("ln(3/2)", q(3, 2))] {
            let (lo, hi) = logarithm(x);
            check(source, lo, hi);
        }
        let (lo, hi) = logarithm(q(2, 1));
        check("ln(1/2)", -hi, -lo);
    }

    #[test]
    #[ignore = "extended independent rational-series oracle for certified functions"]
    fn certified_function_stress() {
        let cases = std::env::var("SYMBOLIC_ENCLOSURE_CASES")
            .map(|s| s.parse::<usize>().expect("decimal case count"))
            .unwrap_or(64);
        let mut seed = 0x656e_636c_6f73_7572u64;
        for case in 0..cases {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let n = (seed >> 32) as i64 % 256 + 1;
            let x = q(n, 257);
            let bits = [96, 192, 384][case % 3];
            let (slo, shi) = trig(x.clone(), true);
            let (clo, chi) = trig(x.clone(), false);
            check_at_precision(&format!("sin({n}/257)"), slo.clone(), shi.clone(), bits);
            check_at_precision(&format!("sin(-{n}/257)"), -shi.clone(), -slo.clone(), bits);
            check_at_precision(&format!("cos({n}/257)"), clo.clone(), chi.clone(), bits);
            check_at_precision(&format!("tan({n}/257)"), &slo / &chi, &shi / &clo, bits);
            let (lo, hi) = exponential(x.clone());
            check_at_precision(&format!("exp({n}/257)"), lo.clone(), hi.clone(), bits);
            check_at_precision(&format!("exp(-{n}/257)"), q(1, 1) / hi, q(1, 1) / lo, bits);
            let (lo, hi) = logarithm(q(1, 1) + &x);
            check_at_precision(
                &format!("ln({}/257)", n + 257),
                lo.clone(),
                hi.clone(),
                bits,
            );
            check_at_precision(&format!("ln(257/{})", n + 257), -hi, -lo, bits);
            let a = super::super::AngleContext::default()
                .parse(&format!("sqrt({n}/257)"))
                .unwrap();
            let ctx = Context::new(bits);
            let result = enclose(&a, &ctx).unwrap();
            let lo = endpoint(&result.lo, &ctx);
            let hi = endpoint(&result.hi, &ctx);
            assert!(&lo * &lo <= x && &hi * &hi >= x, "sqrt case={case}");
        }
        eprintln!(
            "certified function stress: {} rational-reference enclosures passed at 96/192/384 bits",
            cases * 9
        );
    }
    #[test]
    fn rational_sums_products_divisions_and_powers_enclose_literal_values() {
        for (source, want) in [
            ("0.1+0.2", q(3, 10)),
            ("(1/7)^3", q(1, 343)),
            ("(5/7)^-2", q(49, 25)),
            ("1/(2^1075)", R::new(1.into(), BigInt::from(1) << 1075usize)),
        ] {
            let mut core = super::super::AngleContext::default();
            let a = core.parse(source).unwrap();
            let ctx = Context::new(1200);
            let result = enclose(&a, &ctx).unwrap();
            assert!(endpoint(&result.lo, &ctx) <= want && endpoint(&result.hi, &ctx) >= want);
        }
    }
    #[test]
    fn sqrt_bounds_are_checked_by_independent_squared_inequalities() {
        let mut core = super::super::AngleContext::default();
        for (source, input) in [
            ("sqrt(2)", q(2, 1)),
            ("sqrt(1/10)", q(1, 10)),
            (
                "sqrt(1/2^1000)",
                R::new(1.into(), BigInt::from(1) << 1000usize),
            ),
        ] {
            let a = core.parse(source).unwrap();
            let ctx = Context::new(1200);
            let result = enclose(&a, &ctx).unwrap();
            let lo = endpoint(&result.lo, &ctx);
            let hi = endpoint(&result.hi, &ctx);
            assert!(&lo * &lo <= input && &hi * &hi >= input);
        }
    }
    #[test]
    fn shared_nodes_and_nested_functions_have_enclosing_composition() {
        let mut core = super::super::AngleContext::default();
        let a = core.parse("sin(1)+exp(1)+sqrt(2)+ln(2)").unwrap();
        let ctx = Context::new(192);
        let bounds = enclose(&a, &ctx).unwrap();
        let approx = lossy(&a).unwrap();
        use num_traits::ToPrimitive;
        assert!((endpoint(&bounds.middle(), &ctx).to_f64().unwrap() - approx).abs() < 1e-14);
        let a = core.parse("sin(exp(1))+cos(sqrt(2))").unwrap();
        let bounds = enclose(&a, &ctx).unwrap();
        assert!(
            (endpoint(&bounds.middle(), &ctx).to_f64().unwrap() - lossy(&a).unwrap()).abs() < 1e-14
        );
    }
    #[test]
    fn unresolved_domains_and_huge_arguments_fail_without_a_numeric_fallback() {
        let mut core = super::super::AngleContext::default();
        let ctx = Context::new(192);
        for source in [
            "1/sin(1/2^1000)",
            "ln(1/2^1000)",
            "exp(2000)",
            "sin(2^1000)",
        ] {
            let a = core.parse(source).unwrap();
            assert!(enclose(&a, &ctx).is_err(), "{source}");
        }
        let a = core.parse("sin(1000000000000)").unwrap();
        let result = enclose(&a, &ctx).unwrap();
        use num_traits::ToPrimitive;
        assert!(
            (endpoint(&result.middle(), &ctx).to_f64().unwrap() - 1_000_000_000_000_f64.sin())
                .abs()
                < 1e-14
        );
    }
}
