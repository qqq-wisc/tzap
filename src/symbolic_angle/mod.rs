//! Exact, nonperiodic scalar angles used throughout circuits and passes.
//!
//! [`crate::angle::Angle`] re-exports this type. Use an [`AngleContext`] to construct and intern closed
//! real expressions. Arithmetic never numerically evaluates a function.
//!
//! ```
//! use tzap::symbolic_angle::{AngleContext, Proof};
//! let mut ctx = AngleContext::default();
//! let a = ctx.parse("sin(pi/7) + 1/10")?;
//! let b = ctx.parse("-sin(pi/7) + 1/5")?;
//! let sum = ctx.add(&a, &b)?;
//! assert_eq!(sum, ctx.rational(3, 10)?);
//! assert_eq!(sum.prove_zero(), Proof::Disproven);
//! # Ok::<(), tzap::symbolic_angle::Error>(())
//! ```

mod accumulate;
mod context;
mod dyadic;
pub(crate) mod evaluate;
pub mod integration;
mod parse;
mod rational;
pub(crate) mod view;

pub use accumulate::PhaseAccumulator;
pub use context::{AngleContext, InternerStats};
pub use rational::Rational;

use rustc_hash::FxHasher;
use std::{
    cmp::Ordering,
    fmt,
    hash::{Hash, Hasher},
    sync::Arc,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Function {
    Sin,
    Cos,
    Tan,
    Exp,
    Ln,
    Sqrt,
}
impl Function {
    pub const fn name(self) -> &'static str {
        match self {
            Self::Sin => "sin",
            Self::Cos => "cos",
            Self::Tan => "tan",
            Self::Exp => "exp",
            Self::Ln => "ln",
            Self::Sqrt => "sqrt",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Resource {
    Bits,
    Nodes,
    Depth,
    Terms,
    Work,
    Source,
    Output,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    RepresentationLimit,
    RoundedAddition,
    PreservedExpression,
    NonFinite,
    DivisionByZero,
    ZeroPowerZero,
    NonRealPower,
    InvalidDomain(Function),
    DomainUnproved,
    Limit(Resource),
    Parse {
        offset: usize,
        message: &'static str,
    },
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RepresentationLimit => f.write_str("fixed-width compatibility coefficient limit"),
            Self::RoundedAddition => f.write_str("binary64 compatibility arithmetic would round"),
            Self::PreservedExpression => f.write_str("opaque source expression"),
            Self::NonFinite => f.write_str("angle input must be finite"),
            Self::DivisionByZero => {
                f.write_str("angle expression must be finite: division by zero")
            }
            Self::ZeroPowerZero => f.write_str("0^0 is outside the supported real domain"),
            Self::NonRealPower => {
                f.write_str("non-integer power of a non-positive base is outside the real domain")
            }
            Self::InvalidDomain(fun) => {
                write!(f, "argument is outside {}'s real domain", fun.name())
            }
            Self::DomainUnproved => {
                f.write_str("real domain could not be proved by the scalar core")
            }
            Self::Limit(resource) => write!(f, "symbolic angle resource limit: {resource:?}"),
            Self::Parse { offset, message } => {
                write!(f, "angle expression at byte {offset}: {message}")
            }
        }
    }
}
impl std::error::Error for Error {}

/// Limits apply to each operation, plus the context's retained interner entries.
/// Reaching a limit never returns a partial or approximate scalar.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    pub max_bits: usize,
    pub max_nodes: usize,
    pub max_depth: usize,
    pub max_terms: usize,
    pub max_work: usize,
    pub max_source_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_bits: 65_536,
            max_nodes: 4096,
            max_depth: 64,
            max_terms: 256,
            max_work: 4_000_000,
            max_source_bytes: 65_536,
        }
    }
}

pub(super) struct Budget {
    limits: Limits,
    remaining: usize,
}
impl Default for Budget {
    fn default() -> Self {
        Self::new(Limits::default())
    }
}
impl Budget {
    fn new(limits: Limits) -> Self {
        Self {
            remaining: limits.max_work,
            limits,
        }
    }
    fn charge(&mut self, work: usize) -> Result<(), Error> {
        self.remaining = self
            .remaining
            .checked_sub(work)
            .ok_or(Error::Limit(Resource::Work))?;
        Ok(())
    }
    fn bits(&self, bits: usize) -> Result<(), Error> {
        if bits > self.limits.max_bits {
            Err(Error::Limit(Resource::Bits))
        } else {
            Ok(())
        }
    }
    fn arithmetic(&mut self, bits: usize) -> Result<(), Error> {
        let limbs = bits.div_ceil(64).max(1);
        self.charge(
            limbs
                .checked_mul(limbs)
                .ok_or(Error::Limit(Resource::Work))?,
        )
    }
    fn value(&mut self, value: &Angle) -> Result<(), Error> {
        let meta = value.metadata();
        self.bits(meta.bits)?;
        if meta.depth > self.limits.max_depth {
            return Err(Error::Limit(Resource::Depth));
        }
        if meta.nodes > self.limits.max_nodes {
            return Err(Error::Limit(Resource::Nodes));
        }
        self.charge(1)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Proof {
    Proven,
    Disproven,
    Unknown,
}

/// An explicit gate view. Constructing and adding angles never uses a period.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RotationPeriod {
    TwoPi,
    FourPi,
}
impl RotationPeriod {
    fn quarters(self) -> u8 {
        match self {
            Self::TwoPi => 8,
            Self::FourPi => 16,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Sign {
    Negative,
    Zero,
    Positive,
    NonNegative,
    NonPositive,
    Unknown,
}
impl Sign {
    fn from_integer(sign: i8) -> Self {
        match sign {
            -1 => Self::Negative,
            0 => Self::Zero,
            _ => Self::Positive,
        }
    }
    fn neg(self) -> Self {
        match self {
            Self::Negative => Self::Positive,
            Self::Positive => Self::Negative,
            Self::NonNegative => Self::NonPositive,
            Self::NonPositive => Self::NonNegative,
            x => x,
        }
    }
    fn zero(self) -> Proof {
        match self {
            Self::Zero => Proof::Proven,
            Self::Negative | Self::Positive => Proof::Disproven,
            _ => Proof::Unknown,
        }
    }
}

/// Immutable exact literal angle. Small rationals and pi coefficients, and all
/// binary64 inputs, have inline storage. Expressions share an owned DAG.
#[derive(Clone, Debug)]
pub struct Angle(Value);
#[derive(Clone, Debug)]
enum Value {
    Rational { n: i32, d: u32 },
    Pi { n: i32, d: u32 },
    Dyadic(u64),
    Shared(Arc<Node>),
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) struct Term {
    coefficient: Rational,
    atom: Angle,
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum NodeKind {
    Rational(Rational),
    Pi(Rational),
    Sum(Box<[Term]>),
    Product {
        coefficient: Rational,
        factors: Box<[Angle]>,
    },
    Power(Angle, i64),
    Function(Function, Angle),
}
#[derive(Clone, Copy, Debug)]
pub(super) struct Metadata {
    nodes: usize,
    depth: usize,
    bits: usize,
    sign: Sign,
}
#[derive(Debug)]
pub(super) struct Node {
    kind: NodeKind,
    hash: u64,
    meta: Metadata,
}

impl Angle {
    pub const ZERO: Self = Self(Value::Rational { n: 0, d: 1 });
    pub const ONE: Self = Self(Value::Rational { n: 1, d: 1 });
    pub const PI: Self = Self(Value::Pi { n: 1, d: 1 });

    #[inline]
    pub(crate) fn phase_quarters(k: u8) -> Self {
        match k & 7 {
            0 => Self::ZERO,
            1 => Self(Value::Pi { n: 1, d: 4 }),
            2 => Self(Value::Pi { n: 1, d: 2 }),
            3 => Self(Value::Pi { n: 3, d: 4 }),
            4 => Self::PI,
            5 => Self(Value::Pi { n: 5, d: 4 }),
            6 => Self(Value::Pi { n: 3, d: 2 }),
            7 => Self(Value::Pi { n: 7, d: 4 }),
            _ => unreachable!(),
        }
    }
    /// Exact meaning of a supplied binary64 value; no pi inference or rounding.
    /// This constructor allocates nothing, including for subnormal values.
    #[inline]
    pub fn from_f64(value: f64) -> Result<Self, Error> {
        if !value.is_finite() {
            return Err(Error::NonFinite);
        }
        Ok(if value == 0.0 {
            Self::ZERO
        } else {
            Self(Value::Dyadic(value.to_bits()))
        })
    }
    pub fn rational(n: i64, d: i64) -> Result<Self, Error> {
        AngleContext::default().rational(n, d)
    }
    pub fn pi_fraction(n: i64, d: i64) -> Result<Self, Error> {
        AngleContext::default().pi_fraction(n, d)
    }
    /// Converts exact numeric storage to a rational, allocating for large dyadics.
    pub fn as_rational(&self) -> Option<Rational> {
        match &self.0 {
            Value::Rational { n, d } => Some(Rational::canonical(*n as i64, *d as u64)),
            Value::Dyadic(bits) => Some(Rational::from_bits(*bits)),
            Value::Shared(node) => match &node.kind {
                NodeKind::Rational(r) => Some(r.clone()),
                _ => None,
            },
            _ => None,
        }
    }
    /// A numeric SDK parameter only when it equals this literal scalar exactly.
    pub fn to_f64_exact(&self) -> Option<f64> {
        use num_traits::ToPrimitive;
        if let Value::Dyadic(bits) = self.0 {
            return Some(f64::from_bits(bits));
        }
        let r = self.as_rational()?;
        let value = r.to_big_rational().to_f64()?;
        if value.is_finite() && Rational::from_bits(value.to_bits()) == r {
            Some(value)
        } else {
            None
        }
    }
    /// Exact coefficients `(pi, rational)` when the scalar is affine in pi.
    /// This classification never evaluates functions or narrows coefficients.
    pub fn affine_coefficients(&self) -> Option<(Rational, Rational)> {
        if let Some(r) = self.as_rational() {
            return Some((Rational::ZERO, r));
        }
        if let Some(p) = self.pi_coefficient() {
            return Some((p, Rational::ZERO));
        }
        let NodeKind::Sum(terms) = self.kind()? else {
            return None;
        };
        let mut p = Rational::ZERO;
        let mut r = Rational::ZERO;
        for term in terms {
            if term.atom == Self::PI {
                p = term.coefficient.clone();
            } else if term.atom == Self::ONE {
                r = term.coefficient.clone();
            } else {
                return None;
            }
        }
        Some((p, r))
    }
    pub(crate) fn affine_parts(&self) -> Option<(Self, Self)> {
        let NodeKind::Sum(terms) = self.kind()? else {
            return None;
        };
        if terms.len() != 2 {
            return None;
        }
        let mut p = None;
        let mut r = None;
        for term in terms {
            if term.atom == Self::PI {
                p = Some(&term.coefficient);
            } else if term.atom == Self::ONE {
                r = Some(&term.coefficient);
            } else {
                return None;
            }
        }
        let mut ctx = AngleContext::default();
        let mut budget = Budget::default();
        Some((
            ctx.number(p?.clone(), true, &mut budget).ok()?,
            ctx.number(r?.clone(), false, &mut budget).ok()?,
        ))
    }
    pub fn pi_coefficient(&self) -> Option<Rational> {
        match &self.0 {
            Value::Pi { n, d } => Some(Rational::canonical(*n as i64, *d as u64)),
            Value::Shared(node) => match &node.kind {
                NodeKind::Pi(r) => Some(r.clone()),
                _ => None,
            },
            _ => None,
        }
    }
    #[inline]
    pub fn prove_zero(&self) -> Proof {
        self.sign().zero()
    }
    pub fn prove_positive(&self) -> Proof {
        match self.sign() {
            Sign::Positive => Proof::Proven,
            Sign::Negative | Sign::Zero | Sign::NonPositive => Proof::Disproven,
            _ => Proof::Unknown,
        }
    }
    pub fn prove_nonnegative(&self) -> Proof {
        match self.sign() {
            Sign::Positive | Sign::Zero | Sign::NonNegative => Proof::Proven,
            Sign::Negative => Proof::Disproven,
            _ => Proof::Unknown,
        }
    }
    /// Exact quarter-turn classification in an explicit rotation view.
    #[inline]
    pub fn quarter_turns(&self, period: RotationPeriod) -> Option<u8> {
        match &self.0 {
            Value::Rational { n: 0, .. } => Some(0),
            Value::Pi { n, d } => {
                let quarters = match *d {
                    1 => *n as i64 * 4,
                    2 => *n as i64 * 2,
                    4 => *n as i64,
                    _ => {
                        let n = *n as i64 * 4;
                        let d = *d as i64;
                        if n % d != 0 {
                            return None;
                        }
                        n / d
                    }
                };
                Some((quarters as u64 & (period.quarters() as u64 - 1)) as u8)
            }
            Value::Shared(node) => match &node.kind {
                NodeKind::Pi(r) => r.turns(period.quarters()),
                _ => None,
            },
            _ => None,
        }
    }
    pub fn prove_identity(&self, period: RotationPeriod) -> Proof {
        if self.prove_zero() == Proof::Proven {
            return Proof::Proven;
        }
        if let Some(k) = self.quarter_turns(period) {
            return if k == 0 {
                Proof::Proven
            } else {
                Proof::Disproven
            };
        }
        if self.category() <= 1 {
            Proof::Disproven
        } else {
            Proof::Unknown
        }
    }
    #[inline]
    pub fn is_inline(&self) -> bool {
        !matches!(self.0, Value::Shared(_))
    }
    pub fn shares_storage(&self, other: &Self) -> bool {
        match (&self.0, &other.0) {
            (Value::Shared(a), Value::Shared(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }
    pub fn expression_nodes(&self) -> usize {
        self.metadata().nodes
    }
    pub fn expression_depth(&self) -> usize {
        self.metadata().depth
    }
    pub fn coefficient_bits(&self) -> usize {
        self.metadata().bits
    }
    /// Bounds textual expansion of shared nodes. Standard Display is unbounded.
    pub fn to_expression_string(&self, max_bytes: usize) -> Result<String, Error> {
        use fmt::Write;
        struct Writer {
            out: String,
            limit: usize,
        }
        impl fmt::Write for Writer {
            fn write_str(&mut self, s: &str) -> fmt::Result {
                if s.len() > self.limit.saturating_sub(self.out.len()) {
                    return Err(fmt::Error);
                }
                self.out.push_str(s);
                Ok(())
            }
        }
        let mut writer = Writer {
            out: String::new(),
            limit: max_bytes,
        };
        write!(&mut writer, "{self}").map_err(|_| Error::Limit(Resource::Output))?;
        Ok(writer.out)
    }
    fn kind(&self) -> Option<&NodeKind> {
        match &self.0 {
            Value::Shared(n) => Some(&n.kind),
            _ => None,
        }
    }
    #[inline]
    fn sign(&self) -> Sign {
        match &self.0 {
            Value::Shared(node) => node.meta.sign,
            Value::Rational { n, .. } | Value::Pi { n, .. } => Sign::from_integer(n.signum() as i8),
            Value::Dyadic(bits) => {
                if bits >> 63 == 0 {
                    Sign::Positive
                } else {
                    Sign::Negative
                }
            }
        }
    }
    fn metadata(&self) -> Metadata {
        match &self.0 {
            Value::Shared(node) => node.meta,
            Value::Rational { n, d } | Value::Pi { n, d } => Metadata {
                nodes: 1,
                depth: 1,
                bits: (32 - n.unsigned_abs().leading_zeros()).max(32 - d.leading_zeros()) as usize,
                sign: Sign::from_integer(n.signum() as i8),
            },
            Value::Dyadic(bits) => {
                let exponent = ((bits >> 52) & 2047) as i32;
                let significand =
                    (bits & ((1 << 52) - 1)) | if exponent == 0 { 0 } else { 1 << 52 };
                let trailing = significand.trailing_zeros() as i32;
                let power = if exponent == 0 {
                    -1074
                } else {
                    exponent - 1075
                } + trailing;
                let significant = 64 - (significand >> trailing).leading_zeros() as i32;
                Metadata {
                    nodes: 1,
                    depth: 1,
                    bits: if power >= 0 {
                        (significant + power) as usize
                    } else {
                        significant.max(1 - power) as usize
                    },
                    sign: if bits >> 63 == 0 {
                        Sign::Positive
                    } else {
                        Sign::Negative
                    },
                }
            }
        }
    }
    #[inline]
    fn fingerprint(&self) -> u64 {
        match &self.0 {
            Value::Shared(node) => node.hash,
            Value::Rational { n, d } => small_number_hash(0, *n, *d),
            Value::Pi { n, d } => small_number_hash(1, *n, *d),
            Value::Dyadic(bits) => number_hash(0, &Rational::from_bits(*bits)),
        }
    }
    fn category(&self) -> u8 {
        match &self.0 {
            Value::Rational { .. } | Value::Dyadic(_) => 0,
            Value::Pi { .. } => 1,
            Value::Shared(node) => match node.kind {
                NodeKind::Rational(_) => 0,
                NodeKind::Pi(_) => 1,
                NodeKind::Sum(_) => 2,
                NodeKind::Product { .. } => 3,
                NodeKind::Power(..) => 5,
                NodeKind::Function(..) => 6,
            },
        }
    }
    fn body_cmp(&self, rhs: &Self) -> Ordering {
        let order = self.category().cmp(&rhs.category());
        if order != Ordering::Equal {
            return order;
        }
        match self.category() {
            0 => self.as_rational().unwrap().cmp(&rhs.as_rational().unwrap()),
            1 => self
                .pi_coefficient()
                .unwrap()
                .cmp(&rhs.pi_coefficient().unwrap()),
            _ => self.kind().unwrap().cmp(rhs.kind().unwrap()),
        }
    }
}

impl PartialEq for Angle {
    fn eq(&self, rhs: &Self) -> bool {
        match (&self.0, &rhs.0) {
            (Value::Rational { n: a, d: b }, Value::Rational { n: c, d })
            | (Value::Pi { n: a, d: b }, Value::Pi { n: c, d }) => a == c && b == d,
            (Value::Dyadic(a), Value::Dyadic(b)) => a == b,
            (Value::Shared(a), Value::Shared(b)) => {
                Arc::ptr_eq(a, b) || (a.hash == b.hash && a.kind == b.kind)
            }
            _ => self.category() == rhs.category() && self.body_cmp(rhs) == Ordering::Equal,
        }
    }
}
impl Eq for Angle {}
impl Hash for Angle {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.fingerprint().hash(state);
    }
}
impl Ord for Angle {
    fn cmp(&self, rhs: &Self) -> Ordering {
        if self.shares_storage(rhs) {
            return Ordering::Equal;
        }
        self.fingerprint()
            .cmp(&rhs.fingerprint())
            .then_with(|| self.body_cmp(rhs))
    }
}
impl PartialOrd for Angle {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}
impl std::str::FromStr for Angle {
    type Err = Error;
    fn from_str(source: &str) -> Result<Self, Error> {
        AngleContext::default().parse(source)
    }
}
#[inline]
fn small_number_hash(tag: u8, n: i32, d: u32) -> u64 {
    let mut h = FxHasher::default();
    tag.hash(&mut h);
    Rational::hash_small(n as i64, d as u64, &mut h);
    h.finish()
}
fn number_hash(tag: u8, r: &Rational) -> u64 {
    let mut h = FxHasher::default();
    tag.hash(&mut h);
    r.hash(&mut h);
    h.finish()
}
fn node_hash(kind: &NodeKind) -> u64 {
    match kind {
        NodeKind::Rational(r) => number_hash(0, r),
        NodeKind::Pi(r) => number_hash(1, r),
        _ => {
            let mut h = FxHasher::default();
            2u8.hash(&mut h);
            kind.hash(&mut h);
            h.finish()
        }
    }
}
impl fmt::Display for Angle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use view::ExpressionRef;
        match self.expression() {
            ExpressionRef::Rational(r) => write!(f, "{r}"),
            ExpressionRef::Pi(r) => write!(f, "({r}*pi)"),
            ExpressionRef::Sum(terms) => {
                f.write_str("(")?;
                for (i, (coefficient, atom)) in terms.enumerate() {
                    if i != 0 {
                        f.write_str("+")?;
                    }
                    write!(f, "({coefficient}*{atom})")?;
                }
                f.write_str(")")
            }
            ExpressionRef::Product {
                coefficient,
                factors,
            } => {
                write!(f, "({coefficient}")?;
                for factor in factors {
                    write!(f, "*{factor}")?;
                }
                f.write_str(")")
            }
            ExpressionRef::Power(a, e) => write!(f, "({a}^{e})"),
            ExpressionRef::Function(fun, a) => write!(f, "{}({a})", fun.name()),
        }
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod regression_tests;

#[cfg(test)]
mod pass_tests;

#[cfg(test)]
mod stress_tests;
