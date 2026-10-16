//! Exact scalar angles and explicit projective gate helpers.
use crate::circuit::{Circuit, Gate, Qubit};
pub use crate::symbolic_angle::{
    Angle, AngleContext, Error as AngleError, Limits, Proof, RotationPeriod,
};
use crate::symbolic_angle::{Rational, view::ExpressionRef};
use num_traits::ToPrimitive;
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PiFraction {
    numerator: i64,
    denominator: i64,
}
fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
impl PiFraction {
    pub const ZERO: Self = Self {
        numerator: 0,
        denominator: 1,
    };
    pub fn new(numerator: i64, denominator: i64) -> Result<Self, AngleError> {
        Self::wide(numerator as i128, denominator as i128)
    }
    fn wide(mut n: i128, mut d: i128) -> Result<Self, AngleError> {
        if d == 0 {
            return Err(AngleError::DivisionByZero);
        }
        if n == 0 {
            return Ok(Self::ZERO);
        }
        let g = gcd(n.unsigned_abs(), d.unsigned_abs()) as i128;
        n /= g;
        d /= g;
        if d < 0 {
            n = n.checked_neg().ok_or(AngleError::RepresentationLimit)?;
            d = d.checked_neg().ok_or(AngleError::RepresentationLimit)?;
        }
        Ok(Self {
            numerator: n.try_into().map_err(|_| AngleError::RepresentationLimit)?,
            denominator: d.try_into().map_err(|_| AngleError::RepresentationLimit)?,
        })
    }
    pub fn numerator(self) -> i64 {
        self.numerator
    }
    pub fn denominator(self) -> i64 {
        self.denominator
    }
    pub fn checked_add(self, rhs: Self) -> Result<Self, AngleError> {
        let g = gcd(self.denominator as u128, rhs.denominator as u128) as i128;
        Self::wide(
            (self.numerator as i128)
                .checked_mul(rhs.denominator as i128 / g)
                .and_then(|a| {
                    (rhs.numerator as i128)
                        .checked_mul(self.denominator as i128 / g)
                        .and_then(|b| a.checked_add(b))
                })
                .ok_or(AngleError::RepresentationLimit)?,
            (self.denominator as i128)
                .checked_mul(rhs.denominator as i128 / g)
                .ok_or(AngleError::RepresentationLimit)?,
        )
    }
    pub fn checked_neg(self) -> Result<Self, AngleError> {
        Self::wide(-(self.numerator as i128), self.denominator as i128)
    }
    pub fn checked_mul(self, rhs: Self) -> Result<Self, AngleError> {
        let g1 = gcd(
            self.numerator.unsigned_abs() as u128,
            rhs.denominator as u128,
        ) as i128;
        let g2 = gcd(
            rhs.numerator.unsigned_abs() as u128,
            self.denominator as u128,
        ) as i128;
        Self::wide(
            (self.numerator as i128 / g1)
                .checked_mul(rhs.numerator as i128 / g2)
                .ok_or(AngleError::RepresentationLimit)?,
            (self.denominator as i128 / g2)
                .checked_mul(rhs.denominator as i128 / g1)
                .ok_or(AngleError::RepresentationLimit)?,
        )
    }
    pub fn checked_div(self, rhs: Self) -> Result<Self, AngleError> {
        Self::wide(
            (self.numerator as i128)
                .checked_mul(rhs.denominator as i128)
                .ok_or(AngleError::RepresentationLimit)?,
            (self.denominator as i128)
                .checked_mul(rhs.numerator as i128)
                .ok_or(AngleError::RepresentationLimit)?,
        )
    }
    /// Rotation equivalence modulo 2*pi; never use during source arithmetic.
    pub fn normalized_rotation(self) -> Self {
        let d = self.denominator as i128;
        let n = (self.numerator as i128 + d).rem_euclid(2 * d) - d;
        Self::wide(n, d).expect("normalization stays within the coefficient bounds")
    }
    pub fn quarter_turns(self) -> Option<u8> {
        let n = 4 * self.numerator as i128;
        let d = self.denominator as i128;
        (n % d == 0).then(|| (n / d).rem_euclid(8) as u8)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FiniteF64(f64);
impl FiniteF64 {
    pub fn new(value: f64) -> Result<Self, AngleError> {
        if !value.is_finite() {
            Err(AngleError::NonFinite)
        } else {
            Ok(Self(if value.to_bits() & 0x7fff_ffff_ffff_ffff == 0 {
                0.0
            } else {
                value
            }))
        }
    }
    pub fn get(self) -> f64 {
        self.0
    }
}
impl PartialEq for FiniteF64 {
    fn eq(&self, rhs: &Self) -> bool {
        self.0.to_bits() == rhs.0.to_bits()
    }
}
impl Eq for FiniteF64 {}
impl std::hash::Hash for FiniteF64 {
    fn hash<H: std::hash::Hasher>(&self, h: &mut H) {
        self.0.to_bits().hash(h)
    }
}

/// Knuth TwoSum under IEEE round-to-nearest with gradual underflow. A nonzero
/// correction is refused rather than discarded. Rust builds must not use fast-math.
pub fn exact_float_add(a: FiniteF64, b: FiniteF64) -> Result<FiniteF64, AngleError> {
    let s = FiniteF64::new(a.0 + b.0)?.0;
    let bb = FiniteF64::new(s - a.0)?.0;
    let aa = FiniteF64::new(s - bb)?.0;
    let br = FiniteF64::new(b.0 - bb)?.0;
    let ar = FiniteF64::new(a.0 - aa)?.0;
    let correction = FiniteF64::new(ar + br)?.0;
    if correction != 0.0 || !float_sum_matches_bits(a.0, b.0, s) {
        Err(AngleError::RoundedAddition)
    } else {
        FiniteF64::new(s)
    }
}

/// Independent, allocation-free representability check. TwoSum's fast test
/// must never certify a rounded result under altered rounding/underflow modes.
/// Exponent gaps beyond 74 cannot yield a representable sum of two nonzero
/// binary64 values: the smaller significand is too far below the larger one.
fn float_sum_matches_bits(a: f64, b: f64, sum: f64) -> bool {
    fn decode(value: f64) -> (bool, u128, i32) {
        let bits = value.to_bits();
        let exponent = ((bits >> 52) & 2047) as i32;
        let fraction = bits & ((1 << 52) - 1);
        (
            bits >> 63 != 0,
            (fraction | if exponent == 0 { 0 } else { 1 << 52 }) as u128,
            if exponent == 0 {
                -1074
            } else {
                exponent - 1075
            },
        )
    }
    fn normalized(negative: bool, mantissa: u128, exponent: i32) -> (bool, u128, i32) {
        if mantissa == 0 {
            return (false, 0, 0);
        }
        let shift = mantissa.trailing_zeros();
        (negative, mantissa >> shift, exponent + shift as i32)
    }
    let (an, am, ae) = decode(a);
    let (bn, bm, be) = decode(b);
    let (sn, sm, se) = decode(sum);
    let actual = normalized(sn, sm, se);
    if am == 0 {
        return actual == normalized(bn, bm, be);
    }
    if bm == 0 {
        return actual == normalized(an, am, ae);
    }
    let exponent = ae.min(be);
    let (ashift, bshift) = ((ae - exponent) as u32, (be - exponent) as u32);
    if ashift.max(bshift) > 74 {
        return false;
    }
    let (am, bm) = (am << ashift, bm << bshift);
    let (negative, magnitude) = if an == bn {
        (an, am + bm)
    } else if am >= bm {
        (an, am - bm)
    } else {
        (bn, bm - am)
    };
    actual == normalized(negative, magnitude, exponent)
}

impl Angle {
    pub(crate) fn literal_eq(&self, rhs: &Self) -> bool {
        self == rhs
    }
    pub(crate) fn checked_half(&self) -> Result<Self, AngleError> {
        AngleContext::default().half(self)
    }
    pub(crate) fn literal_neg(&self) -> Result<Self, AngleError> {
        AngleContext::default().neg(self)
    }
    /// Literal scalar negation; callers select a gate period separately.
    pub fn checked_neg(&self) -> Result<Self, AngleError> {
        self.literal_neg()
    }
    pub(crate) fn turns(k: u8) -> Self {
        // Named P/T phases are intentionally viewed modulo 2*pi.
        Self::phase_quarters(k)
    }
    /// Compatibility classification: available only for a bounded pi coefficient
    /// and an exactly representable binary64 residual. Never used for arithmetic.
    pub fn components(&self) -> Option<(PiFraction, FiniteF64)> {
        if let Some(r) = self.to_f64_exact() {
            return Some((PiFraction::ZERO, FiniteF64(r)));
        }
        if let Some(p) = self.pi_coefficient() {
            let (n, d) = p.as_small()?;
            return Some((PiFraction::new(n, d.try_into().ok()?).ok()?, FiniteF64(0.0)));
        }
        if let ExpressionRef::Sum(terms) = self.expression() {
            let mut p = PiFraction::ZERO;
            let mut r = 0.0;
            for (coefficient, atom) in terms {
                if atom == &Self::PI {
                    let (n, d) = coefficient.as_small()?;
                    p = PiFraction::new(n, d.try_into().ok()?).ok()?;
                } else if atom == &Self::ONE {
                    let big = coefficient.to_big_rational();
                    r = big.to_f64()?;
                    if Rational::from_bits(r.to_bits()) != *coefficient {
                        return None;
                    }
                } else {
                    return None;
                }
            }
            return Some((p, FiniteF64::new(r).ok()?));
        }
        None
    }
    pub fn is_preserved(&self) -> bool {
        false
    }
    pub fn used_numeric_fallback(&self) -> bool {
        false
    }
    pub fn is_zero(&self) -> bool {
        self.prove_zero() == Proof::Proven
    }
    pub fn checked_add_signed(&self, sign: i8, rhs: &Self) -> Result<Self, AngleError> {
        self.add_signed_with(&mut AngleContext::default(), sign, rhs)
    }
    pub(crate) fn add_signed_with(
        &self,
        ctx: &mut AngleContext,
        sign: i8,
        rhs: &Self,
    ) -> Result<Self, AngleError> {
        assert!(sign == 1 || sign == -1);
        let result = if sign == 1 {
            ctx.add(self, rhs)
        } else {
            ctx.sub(rhs, self)
        };
        if let Err(error) = result {
            crate::angle_stats::failure(error);
        }
        result
    }
    #[inline]
    pub(crate) fn merge_phases_with(
        &self,
        ctx: &mut AngleContext,
        sign: i8,
        rhs: &Self,
    ) -> Result<Self, AngleError> {
        let result = ctx.merge_phases(self, sign, rhs);
        if let Err(error) = result {
            crate::angle_stats::failure(error);
        }
        result
    }
    pub(crate) fn projective_quarters(&self) -> Option<u8> {
        self.quarter_turns(RotationPeriod::TwoPi)
    }
    pub(crate) fn projective_identity(&self) -> bool {
        self.prove_identity(RotationPeriod::TwoPi) == Proof::Proven
    }
    /// Explicit, uncertified evaluation for diagnostics and lossy export only.
    pub fn to_f64_lossy(&self) -> Result<f64, AngleError> {
        crate::symbolic_angle::evaluate::lossy(self)
    }
    pub(crate) fn emitted_gate_count(&self) -> usize {
        match self.projective_quarters() {
            Some(0) => 0,
            Some(3 | 5) => 2,
            None => self
                .affine_parts()
                .map_or(1, |(p, _)| p.emitted_gate_count() + 1),
            _ => 1,
        }
    }
    pub(crate) fn try_emit(&self, q: Qubit, push: &mut impl FnMut(Gate) -> bool) -> bool {
        match self.projective_quarters() {
            Some(0) => true,
            Some(1) => push(Gate::t(q)),
            Some(2) => push(Gate::s(q)),
            Some(3) => push(Gate::s(q)) && push(Gate::t(q)),
            Some(4) => push(Gate::z(q)),
            Some(5) => push(Gate::z(q)) && push(Gate::t(q)),
            Some(6) => push(Gate::sdg(q)),
            Some(7) => push(Gate::tdg(q)),
            _ => match self.affine_parts() {
                Some((p, r)) => p.try_emit(q, push) && push(Gate::rz(r, q)),
                None => push(Gate::rz(self.clone(), q)),
            },
        }
    }
    pub(crate) fn rotation_gates(&self, q: Qubit) -> Vec<Gate> {
        let mut gates = Vec::with_capacity(self.emitted_gate_count());
        self.try_emit(q, &mut |g| {
            gates.push(g);
            true
        });
        gates
    }
    pub(crate) fn emit(&self, out: &mut Circuit, q: Qubit) {
        self.try_emit(q, &mut |g| {
            out.apply(g);
            true
        });
    }
    pub(crate) fn of_gate(g: &Gate) -> Option<Self> {
        match g {
            Gate::t(_) => Some(Self::turns(1)),
            Gate::tdg(_) => Some(Self::turns(7)),
            Gate::s(_) => Some(Self::turns(2)),
            Gate::sdg(_) => Some(Self::turns(6)),
            Gate::z(_) => Some(Self::turns(4)),
            Gate::rz(a, _) | Gate::p(a, _) | Gate::rx(a, _) | Gate::ry(a, _) => Some(a.clone()),
            _ => None,
        }
    }
}
pub(crate) fn real_token(value: f64) -> String {
    let s = value.to_string();
    let (mantissa, exp) = s
        .split_once(['e', 'E'])
        .map_or((s.as_str(), None), |(m, e)| (m, Some(e)));
    let mut out = mantissa.to_string();
    if !out.contains('.') {
        out.push_str(".0");
    }
    if let Some(e) = exp {
        out.push('e');
        out.push_str(e);
    }
    out
}
/// Lower only certified quarter-turns; preserve every other angle verbatim.
///
/// The result equals the input only up to global phase: `rz(kπ/4)` becomes
/// named phase gates (`rz(π)` becomes `z`, a factor of `-i` apart). Apply it to
/// a whole circuit, never to a subcircuit that is later controlled.
pub fn lower_exact_rotations(circuit: &Circuit) -> Circuit {
    let mut out = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
    for g in &circuit.gates {
        match g {
            Gate::rz(a, q) if a.projective_quarters().is_some() => a.emit(&mut out, *q),
            _ => out.apply(g.clone()),
        }
    }
    out
}

/// Borrow the common Clifford+T/numeric path; allocate only when lowering is needed.
pub(crate) fn lowered_if_needed(circuit: &Circuit) -> std::borrow::Cow<'_, Circuit> {
    if circuit
        .gates
        .iter()
        .any(|g| matches!(g, Gate::rz(a, _) if a.projective_quarters().is_some()))
    {
        std::borrow::Cow::Owned(lower_exact_rotations(circuit))
    } else {
        std::borrow::Cow::Borrowed(circuit)
    }
}

#[cfg(test)]
mod arithmetic_guard_tests {
    use super::float_sum_matches_bits;

    #[test]
    fn integer_guard_rejects_flushed_and_rounded_candidate_sums() {
        // Exercise altered-environment results directly, without mutating the
        // process floating-point environment used by concurrent tests.
        let tiny = f64::from_bits(1);
        assert!(!float_sum_matches_bits(tiny, tiny, 0.0));
        assert!(!float_sum_matches_bits(
            tiny,
            f64::from_bits((1 << 63) | 2),
            0.0
        ));
        assert!(float_sum_matches_bits(tiny, tiny, f64::from_bits(2)));
        assert!(float_sum_matches_bits(tiny, -tiny, -0.0));
        for rounded in [1e16, f64::from_bits(1e16_f64.to_bits() + 1)] {
            assert!(!float_sum_matches_bits(1e16, 1.0, rounded));
        }
        assert!(!float_sum_matches_bits(0.1, 0.2, 0.3));
        assert!(float_sum_matches_bits(0.125, 0.25, 0.375));
    }
}
