//! Checked angle values. Floats denote their exact binary64 value; explicit
//! pi coefficients denote rational multiples of mathematical pi.
use crate::circuit::{Circuit, Gate, Qubit};
use std::{fmt, sync::Arc};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AngleError {
    NonFinite,
    DivisionByZero,
    RepresentationLimit,
    RoundedAddition,
    PreservedExpression,
}
impl fmt::Display for AngleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::NonFinite => "angle must be finite",
            Self::DivisionByZero => "angle expression must be finite: division by zero",
            Self::RepresentationLimit => "angle coefficient exceeds the fixed-width representation",
            Self::RoundedAddition => "angle addition would round",
            Self::PreservedExpression => "preserved angle expression cannot be folded",
        })
    }
}
impl std::error::Error for AngleError {}

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
    pub(crate) fn decimal(source: &str) -> Result<Self, AngleError> {
        let (mantissa, exponent) = source.split_once(['e', 'E']).unwrap_or((source, "0"));
        let exponent: i32 = exponent
            .parse()
            .map_err(|_| AngleError::RepresentationLimit)?;
        let (whole, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        let mut digits = format!("{whole}{frac}");
        let mut scale = (frac.len() as i64) - (exponent as i64);
        while digits.ends_with('0') && digits.len() > 1 {
            digits.pop();
            scale -= 1;
        }
        let n: i128 = digits
            .parse()
            .map_err(|_| AngleError::RepresentationLimit)?;
        if n == 0 {
            return Ok(Self::ZERO);
        }
        if scale.unsigned_abs() > 38 {
            return Err(AngleError::RepresentationLimit);
        }
        let power = 10i128
            .checked_pow(scale.unsigned_abs() as u32)
            .ok_or(AngleError::RepresentationLimit)?;
        if scale >= 0 {
            Self::wide(n, power)
        } else {
            Self::wide(
                n.checked_mul(power)
                    .ok_or(AngleError::RepresentationLimit)?,
                1,
            )
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FiniteF64(f64);
impl FiniteF64 {
    pub fn new(value: f64) -> Result<Self, AngleError> {
        if !value.is_finite() {
            Err(AngleError::NonFinite)
        } else {
            Ok(Self(if value == 0.0 { 0.0 } else { value }))
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
    if correction != 0.0 {
        Err(AngleError::RoundedAddition)
    } else {
        FiniteF64::new(s)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Symbolic {
    Affine { pi: PiFraction, radians: FiniteF64 },
    Preserved { source: String, line: usize },
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Value {
    Float(FiniteF64),
    NumericFallback(FiniteF64),
    Turns(u8),
    Symbolic(Arc<Symbolic>),
}
/// Immutable, validated angle. No public unchecked payload construction.
#[derive(Clone, Debug)]
pub struct Angle(Value);
impl PartialEq for Angle {
    fn eq(&self, rhs: &Self) -> bool {
        match (self.components(), rhs.components()) {
            (Some((a, ar)), Some((b, br))) => {
                a.normalized_rotation() == b.normalized_rotation() && ar == br
            }
            (None, None) => match (&self.0, &rhs.0) {
                (Value::Symbolic(a), Value::Symbolic(b)) => match (&**a, &**b) {
                    (
                        Symbolic::Preserved { source: a, .. },
                        Symbolic::Preserved { source: b, .. },
                    ) => a == b,
                    _ => false,
                },
                _ => false,
            },
            _ => false,
        }
    }
}
impl Eq for Angle {}
impl std::hash::Hash for Angle {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        if let Some((pi, r)) = self.components() {
            0u8.hash(state);
            pi.normalized_rotation().hash(state);
            r.hash(state);
        } else if let Value::Symbolic(s) = &self.0 {
            if let Symbolic::Preserved { source, .. } = &**s {
                1u8.hash(state);
                source.hash(state);
            }
        }
    }
}
impl Angle {
    /// Compare literal values without the projective 2π rotation equivalence.
    pub(crate) fn literal_eq(&self, rhs: &Self) -> bool {
        match (self.components(), rhs.components()) {
            (Some(a), Some(b)) => a == b,
            (None, None) => self == rhs,
            _ => false,
        }
    }
    /// Exact half-angle for controlled lowering; never reduce modulo 2π.
    pub(crate) fn checked_half(&self) -> Result<Self, AngleError> {
        let (pi, r) = self.components().ok_or(AngleError::PreservedExpression)?;
        let half = r.get() / 2.0;
        if half * 2.0 != r.get() {
            return Err(AngleError::RoundedAddition);
        }
        Ok(Self::affine(
            pi.checked_div(PiFraction::new(2, 1)?)?,
            FiniteF64::new(half)?,
        ))
    }
    pub(crate) fn literal_neg(&self) -> Result<Self, AngleError> {
        let (pi, r) = self.components().ok_or(AngleError::PreservedExpression)?;
        Ok(Self::affine(pi.checked_neg()?, FiniteF64::new(-r.get())?))
    }
    pub fn from_f64(value: f64) -> Result<Self, AngleError> {
        Ok(Self(Value::Float(FiniteF64::new(value)?)))
    }
    pub fn pi_fraction(n: i64, d: i64) -> Result<Self, AngleError> {
        Ok(Self::affine(PiFraction::new(n, d)?, FiniteF64(0.0)))
    }
    pub(crate) fn turns(k: u8) -> Self {
        Self(Value::Turns(k & 7))
    }
    pub(crate) fn numeric_fallback(value: f64) -> Result<Self, AngleError> {
        Ok(Self(Value::NumericFallback(FiniteF64::new(value)?)))
    }
    pub(crate) fn preserved(source: String, line: usize) -> Self {
        Self(Value::Symbolic(Arc::new(Symbolic::Preserved {
            source,
            line,
        })))
    }
    pub(crate) fn affine(pi: PiFraction, radians: FiniteF64) -> Self {
        if pi.numerator == 0 {
            Self(Value::Float(radians))
        } else {
            Self(Value::Symbolic(Arc::new(Symbolic::Affine { pi, radians })))
        }
    }
    pub fn components(&self) -> Option<(PiFraction, FiniteF64)> {
        match &self.0 {
            Value::Float(r) | Value::NumericFallback(r) => Some((PiFraction::ZERO, *r)),
            Value::Turns(k) => Some((PiFraction::new(*k as i64, 4).unwrap(), FiniteF64(0.0))),
            Value::Symbolic(s) => match **s {
                Symbolic::Affine { pi, radians } => Some((pi, radians)),
                _ => None,
            },
        }
    }
    pub fn is_preserved(&self) -> bool {
        self.components().is_none()
    }
    pub fn used_numeric_fallback(&self) -> bool {
        matches!(self.0, Value::NumericFallback(_))
    }
    pub fn quarter_turns(&self) -> Option<u8> {
        match &self.0 {
            Value::Turns(k) => Some(*k),
            Value::Float(r) | Value::NumericFallback(r) => {
                if r.get() == 0.0 {
                    Some(0)
                } else {
                    None
                }
            }
            Value::Symbolic(s) => match &**s {
                Symbolic::Affine { pi, radians } if radians.get() == 0.0 => pi.quarter_turns(),
                _ => None,
            },
        }
    }
    pub fn is_zero(&self) -> bool {
        self.quarter_turns() == Some(0)
    }
    pub fn checked_neg(&self) -> Result<Self, AngleError> {
        if let Some(k) = self.quarter_turns() {
            return Ok(Self::turns(8u8.wrapping_sub(k)));
        }
        let (pi, r) = self.components().ok_or(AngleError::PreservedExpression)?;
        Ok(Self::affine(
            pi.normalized_rotation().checked_neg()?,
            FiniteF64::new(-r.0)?,
        )
        .normalized())
    }
    pub fn checked_add_signed(&self, sign: i8, rhs: &Self) -> Result<Self, AngleError> {
        let result = self.add_signed(sign, rhs);
        if let Err(error) = result {
            crate::angle_stats::failure(error);
        }
        result
    }
    fn add_signed(&self, sign: i8, rhs: &Self) -> Result<Self, AngleError> {
        assert!(sign == 1 || sign == -1);
        if let (Some(a), Some(b)) = (self.quarter_turns(), rhs.quarter_turns()) {
            return Ok(Self::turns(
                (sign as i16 * a as i16 + b as i16).rem_euclid(8) as u8,
            ));
        }
        let (a, ar) = self.components().ok_or(AngleError::PreservedExpression)?;
        let (b, br) = rhs.components().ok_or(AngleError::PreservedExpression)?;
        let a = a.normalized_rotation();
        let a = if sign < 0 { a.checked_neg()? } else { a };
        let pi = a
            .checked_add(b.normalized_rotation())?
            .normalized_rotation();
        let radians = exact_float_add(FiniteF64(sign as f64 * ar.0), br)?;
        Ok(Self::affine(pi, radians))
    }
    pub fn normalized(&self) -> Self {
        if matches!(
            self.0,
            Value::Turns(_) | Value::Float(_) | Value::NumericFallback(_)
        ) {
            return self.clone();
        }
        match self.components() {
            Some((p, r)) => {
                let p = p.normalized_rotation();
                if r.0 == 0.0 {
                    if let Some(k) = p.quarter_turns() {
                        return Self::turns(k);
                    }
                }
                Self::affine(p, r)
            }
            None => self.clone(),
        }
    }
    /// Explicit, uncertified numerical conversion for synthesis/test oracles.
    /// This is never used by folding or exact classification.
    pub fn to_f64_lossy(&self) -> Result<f64, AngleError> {
        let (pi, r) = self.components().ok_or(AngleError::PreservedExpression)?;
        FiniteF64::new((pi.numerator as f64 / pi.denominator as f64) * std::f64::consts::PI + r.0)
            .map(FiniteF64::get)
    }
    pub(crate) fn emitted_gate_count(&self) -> usize {
        match self.quarter_turns() {
            Some(0) => 0,
            Some(3 | 5) => 2,
            Some(_) => 1,
            None => match self.components() {
                Some((pi, r)) if pi.numerator() != 0 && r.get() != 0.0 => {
                    Self::affine(pi, FiniteF64(0.0)).emitted_gate_count() + 1
                }
                _ => 1,
            },
        }
    }
    pub(crate) fn try_emit(&self, q: Qubit, push: &mut impl FnMut(Gate) -> bool) -> bool {
        match self.quarter_turns() {
            Some(0) => true,
            Some(1) => push(Gate::t(q)),
            Some(2) => push(Gate::s(q)),
            Some(3) => push(Gate::s(q)) && push(Gate::t(q)),
            Some(4) => push(Gate::z(q)),
            Some(5) => push(Gate::z(q)) && push(Gate::t(q)),
            Some(6) => push(Gate::sdg(q)),
            Some(7) => push(Gate::tdg(q)),
            _ => match self.components() {
                Some((pi, r)) if pi.numerator != 0 && r.0 != 0.0 => {
                    Self::affine(pi.normalized_rotation(), FiniteF64(0.0)).try_emit(q, push)
                        && push(Gate::rz(Self(Value::Float(r)), q))
                }
                _ => push(Gate::rz(self.normalized(), q)),
            },
        }
    }
    pub(crate) fn rotation_gates(&self, q: Qubit) -> Vec<Gate> {
        let mut gates = Vec::with_capacity(self.emitted_gate_count());
        self.try_emit(q, &mut |gate| {
            gates.push(gate);
            true
        });
        gates
    }
    pub(crate) fn emit(&self, output: &mut Circuit, q: Qubit) {
        self.try_emit(q, &mut |gate| {
            output.apply(gate);
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
impl fmt::Display for Angle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Value::Symbolic(s) = &self.0 {
            if let Symbolic::Preserved { source, .. } = &**s {
                return f.write_str(source);
            }
        }
        let (pi, r) = self.components().unwrap();
        if pi.numerator == 0 {
            return write!(f, "{}", real_token(r.0));
        }
        write!(f, "({}*pi/{})", pi.numerator, pi.denominator)?;
        if r.0 != 0.0 {
            write!(f, "+({})", real_token(r.0))?;
        }
        Ok(())
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
pub fn lower_exact_rotations(circuit: &Circuit) -> Circuit {
    let mut out = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
    for g in &circuit.gates {
        match g {
            Gate::rz(a, q) if a.quarter_turns().is_some() => a.emit(&mut out, *q),
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
        .any(|g| matches!(g, Gate::rz(a, _) if a.quarter_turns().is_some()))
    {
        std::borrow::Cow::Owned(lower_exact_rotations(circuit))
    } else {
        std::borrow::Cow::Borrowed(circuit)
    }
}
