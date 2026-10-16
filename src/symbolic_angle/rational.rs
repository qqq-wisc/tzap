use super::{Budget, Error, Resource};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::{
    cmp::Ordering,
    fmt,
    hash::{Hash, Hasher},
    sync::Arc,
};

/// A reduced rational. Machine-sized values allocate nothing; larger values
/// share immutable arbitrary-precision storage. Construction is checked.
#[derive(Clone, Debug)]
pub struct Rational(Representation);

#[derive(Clone, Debug)]
enum Representation {
    Small { numerator: i64, denominator: u64 },
    Big(Arc<BigRational>),
}

impl Rational {
    pub const ZERO: Self = Self(Representation::Small {
        numerator: 0,
        denominator: 1,
    });
    pub const ONE: Self = Self(Representation::Small {
        numerator: 1,
        denominator: 1,
    });

    pub fn new(numerator: i64, denominator: i64) -> Result<Self, Error> {
        Self::wide(
            numerator as i128,
            denominator as i128,
            &mut Budget::default(),
        )
    }

    pub const fn integer(numerator: i64) -> Self {
        Self(Representation::Small {
            numerator,
            denominator: 1,
        })
    }
    pub(super) const fn canonical(numerator: i64, denominator: u64) -> Self {
        Self(Representation::Small {
            numerator,
            denominator,
        })
    }

    #[inline]
    pub(super) fn hash_small<H: Hasher>(n: i64, d: u64, state: &mut H) {
        0_u8.hash(state);
        n.hash(state);
        d.hash(state);
    }
    pub fn as_small(&self) -> Option<(i64, u64)> {
        match &self.0 {
            Representation::Small {
                numerator,
                denominator,
            } => Some((*numerator, *denominator)),
            Representation::Big(_) => None,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.as_small() == Some((0, 1))
    }
    pub fn is_one(&self) -> bool {
        self.as_small() == Some((1, 1))
    }
    pub fn is_integer(&self) -> bool {
        match &self.0 {
            Representation::Small { denominator, .. } => *denominator == 1,
            Representation::Big(value) => value.is_integer(),
        }
    }
    pub fn signum(&self) -> i8 {
        match &self.0 {
            Representation::Small { numerator, .. } => numerator.signum() as i8,
            Representation::Big(value) => {
                if value.is_positive() {
                    1
                } else {
                    -1
                }
            }
        }
    }
    /// Largest numerator/denominator bit length, excluding the sign.
    pub fn bits(&self) -> usize {
        match &self.0 {
            Representation::Small {
                numerator,
                denominator,
            } => (64 - numerator.unsigned_abs().leading_zeros())
                .max(64 - denominator.leading_zeros()) as usize,
            Representation::Big(value) => value.numer().bits().max(value.denom().bits()) as usize,
        }
    }
    /// Explicit conversion to an owned arbitrary-precision rational.
    pub fn to_big_rational(&self) -> BigRational {
        match &self.0 {
            Representation::Small {
                numerator,
                denominator,
            } => BigRational::new_raw((*numerator).into(), (*denominator).into()),
            Representation::Big(value) => (**value).clone(),
        }
    }

    pub(super) fn from_big(n: BigInt, d: BigInt, budget: &mut Budget) -> Result<Self, Error> {
        if d.is_zero() {
            return Err(Error::DivisionByZero);
        }
        // Charge for the reduction; the size limit applies to its result.
        budget.arithmetic(n.bits().max(d.bits()) as usize)?;
        Self::from_reduced(BigRational::new(n, d), budget)
    }
    fn from_reduced(value: BigRational, budget: &mut Budget) -> Result<Self, Error> {
        budget.bits(value.numer().bits().max(value.denom().bits()) as usize)?;
        if let (Some(n), Some(d)) = (value.numer().to_i64(), value.denom().to_u64()) {
            Ok(Self(Representation::Small {
                numerator: n,
                denominator: d,
            }))
        } else {
            Ok(Self(Representation::Big(Arc::new(value))))
        }
    }
    pub(super) fn wide(mut n: i128, mut d: i128, budget: &mut Budget) -> Result<Self, Error> {
        if d == 0 {
            return Err(Error::DivisionByZero);
        }
        if n == 0 {
            return Ok(Self::ZERO);
        }
        if d > 0 && n.unsigned_abs() <= u64::MAX as u128 && d <= u64::MAX as i128 {
            let g = (n.unsigned_abs() as u64).gcd(&(d as u64));
            n /= g as i128;
            d /= g as i128;
            if let Ok(n) = i64::try_from(n) {
                let result = Self::canonical(n, d as u64);
                budget.bits(result.bits())?;
                return Ok(result);
            }
        }
        let g = n.unsigned_abs().gcd(&d.unsigned_abs());
        // g fits i128 for all callers' checked, positive-denominator products.
        if g > i128::MAX as u128 || (d < 0 && (n == i128::MIN || d == i128::MIN)) {
            return Self::from_big(n.into(), d.into(), budget);
        }
        n /= g as i128;
        d /= g as i128;
        if d < 0 {
            n = -n;
            d = -d;
        }
        if let (Ok(n), Ok(d)) = (i64::try_from(n), u64::try_from(d)) {
            let value = Self(Representation::Small {
                numerator: n,
                denominator: d,
            });
            budget.bits(value.bits())?;
            Ok(value)
        } else {
            Self::from_big(n.into(), d.into(), budget)
        }
    }

    pub(super) fn add(&self, rhs: &Self, budget: &mut Budget) -> Result<Self, Error> {
        budget.charge(1)?;
        if self.is_zero() {
            return Ok(rhs.clone());
        }
        if rhs.is_zero() {
            return Ok(self.clone());
        }
        if let (Some((a, b)), Some((c, d))) = (self.as_small(), rhs.as_small()) {
            if b == d {
                return Self::wide(a as i128 + c as i128, b as i128, budget);
            }
            let g = b.gcd(&d) as i128;
            let n = (a as i128).checked_mul(d as i128 / g).and_then(|x| {
                (c as i128)
                    .checked_mul(b as i128 / g)
                    .and_then(|y| x.checked_add(y))
            });
            let denominator = (b as i128).checked_mul(d as i128 / g);
            if let (Some(n), Some(d)) = (n, denominator) {
                return Self::wide(n, d, budget);
            }
        }
        budget.arithmetic(self.bits().max(rhs.bits()))?;
        for (big, small) in [(self, rhs), (rhs, self)] {
            if let Representation::Big(value) = &big.0
                && let Some((n, 1)) = small.as_small()
            {
                // gcd(numerator + n*denominator, denominator) remains one.
                return Self::from_reduced(
                    BigRational::new_raw(value.numer() + value.denom() * n, value.denom().clone()),
                    budget,
                );
            }
        }
        Self::from_reduced(self.to_big_rational() + rhs.to_big_rational(), budget)
    }
    pub(super) fn neg(&self, budget: &mut Budget) -> Result<Self, Error> {
        budget.charge(1)?;
        if let Some((n, d)) = self.as_small()
            && let Some(n) = n.checked_neg()
        {
            return Ok(Self(Representation::Small {
                numerator: n,
                denominator: d,
            }));
        }
        budget.arithmetic(self.bits())?;
        Self::from_reduced(-self.to_big_rational(), budget)
    }
    pub(super) fn mul(&self, rhs: &Self, budget: &mut Budget) -> Result<Self, Error> {
        budget.charge(1)?;
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::ZERO);
        }
        if self.is_one() {
            return Ok(rhs.clone());
        }
        if rhs.is_one() {
            return Ok(self.clone());
        }
        if let (Some((a, b)), Some((c, d))) = (self.as_small(), rhs.as_small()) {
            let g1 = a.unsigned_abs().gcd(&d) as i128;
            let g2 = c.unsigned_abs().gcd(&b) as i128;
            let n = (a as i128 / g1).checked_mul(c as i128 / g2);
            let d = (b as i128 / g2).checked_mul(d as i128 / g1);
            if let (Some(n), Some(d)) = (n, d) {
                return Self::wide(n, d, budget);
            }
        }
        budget.arithmetic(self.bits().max(rhs.bits()))?;
        Self::from_reduced(self.to_big_rational() * rhs.to_big_rational(), budget)
    }
    pub(super) fn div(&self, rhs: &Self, budget: &mut Budget) -> Result<Self, Error> {
        if rhs.is_zero() {
            return Err(Error::DivisionByZero);
        }
        budget.charge(1)?;
        if self.is_zero() {
            return Ok(Self::ZERO);
        }
        if self == rhs {
            return Ok(Self::ONE);
        }
        if let (Some((a, b)), Some((c, d))) = (self.as_small(), rhs.as_small()) {
            let g1 = a.unsigned_abs().gcd(&c.unsigned_abs()) as i128;
            let g2 = d.gcd(&b) as i128;
            let n = (a as i128 / g1).checked_mul(d as i128 / g2);
            let den = (b as i128 / g2).checked_mul(c as i128 / g1);
            if let (Some(n), Some(d)) = (n, den) {
                return Self::wide(n, d, budget);
            }
        }
        budget.arithmetic(self.bits().max(rhs.bits()))?;
        Self::from_reduced(self.to_big_rational() / rhs.to_big_rational(), budget)
    }
    pub(super) fn pow(&self, exponent: i64, budget: &mut Budget) -> Result<Self, Error> {
        if self.is_zero() {
            return match exponent.cmp(&0) {
                Ordering::Greater => Ok(Self::ZERO),
                Ordering::Equal => Err(Error::ZeroPowerZero),
                Ordering::Less => Err(Error::DivisionByZero),
            };
        }
        if exponent == 0 || self.is_one() {
            return Ok(Self::ONE);
        }
        if self.as_small() == Some((-1, 1)) {
            return Ok(Self::integer(if exponent & 1 == 0 { 1 } else { -1 }));
        }
        let e =
            usize::try_from(exponent.unsigned_abs()).map_err(|_| Error::Limit(Resource::Work))?;
        budget.bits(
            self.bits()
                .checked_mul(e)
                .ok_or(Error::Limit(Resource::Bits))?,
        )?;
        let mut result = Self::ONE;
        let mut base = self.clone();
        let mut e = exponent.unsigned_abs();
        while e != 0 {
            if e & 1 != 0 {
                result = result.mul(&base, budget)?;
            }
            e >>= 1;
            if e != 0 {
                base = base.mul(&base, budget)?;
            }
        }
        if exponent < 0 {
            Self::ONE.div(&result, budget)
        } else {
            Ok(result)
        }
    }
    pub(super) fn perfect_sqrt(&self, budget: &mut Budget) -> Result<Option<Self>, Error> {
        if self.signum() < 0 {
            return Ok(None);
        }
        if let Some((n, d)) = self.as_small() {
            let a = (n as u64).isqrt();
            let b = d.isqrt();
            return if a as u128 * a as u128 == n as u128 && b as u128 * b as u128 == d as u128 {
                Ok(Some(Self::wide(a as i128, b as i128, budget)?))
            } else {
                Ok(None)
            };
        }
        budget.arithmetic(self.bits())?;
        let value = self.to_big_rational();
        let n = value.numer().sqrt();
        let d = value.denom().sqrt();
        if &n * &n == *value.numer() && &d * &d == *value.denom() {
            Ok(Some(Self::from_big(n, d, budget)?))
        } else {
            Ok(None)
        }
    }
    pub(crate) fn from_bits(bits: u64) -> Self {
        let exponent = ((bits >> 52) & 2047) as i32;
        let mut n = bits & ((1 << 52) - 1);
        if exponent != 0 {
            n |= 1 << 52;
        }
        if n == 0 {
            return Self::ZERO;
        }
        let mut power = if exponent == 0 {
            -1074
        } else {
            exponent - 1075
        };
        let trailing = n.trailing_zeros();
        n >>= trailing;
        power += trailing as i32;
        let negative = bits >> 63 != 0;
        if power < 0 && -power < 64 {
            return Self(Representation::Small {
                numerator: if negative { -(n as i64) } else { n as i64 },
                denominator: 1 << -power,
            });
        }
        if power >= 0 && (64 - n.leading_zeros()) as i32 + power <= 63 {
            let n = (n << power) as i64;
            return Self(Representation::Small {
                numerator: if negative { -n } else { n },
                denominator: 1,
            });
        }
        let mut n = BigInt::from(n);
        if negative {
            n = -n;
        }
        let value = if power >= 0 {
            BigRational::from_integer(n << power as usize)
        } else {
            BigRational::new(n, BigInt::one() << -power as usize)
        };
        // Binary64 input is bounded by 1075 bits; normalize to the machine path.
        if let (Some(n), Some(d)) = (value.numer().to_i64(), value.denom().to_u64()) {
            Self(Representation::Small {
                numerator: n,
                denominator: d,
            })
        } else {
            Self(Representation::Big(Arc::new(value)))
        }
    }
    pub(super) fn turns(&self, period: u8) -> Option<u8> {
        if let Some((n, d)) = self.as_small() {
            let n = n as i128 * 4;
            let d = d as i128;
            return (n % d == 0).then(|| (n / d).rem_euclid(period as i128) as u8);
        }
        let value = self.to_big_rational();
        let n = value.numer() * 4;
        if &n % value.denom() != BigInt::zero() {
            return None;
        }
        let k: BigInt = n / value.denom();
        k.mod_floor(&BigInt::from(period)).to_u8()
    }
    pub(super) fn reduced_pi(&self, budget: &mut Budget) -> Result<Self, Error> {
        if let Some((n, d)) = self.as_small() {
            return Self::wide(
                (n as i128 + d as i128).rem_euclid(2 * d as i128) - d as i128,
                d as i128,
                budget,
            );
        }
        budget.arithmetic(self.bits())?;
        let value = self.to_big_rational();
        let d = value.denom();
        Self::from_big(
            (value.numer() + d).mod_floor(&(d * 2)) - d,
            d.clone(),
            budget,
        )
    }
}

// Every constructor reduces and demotes its result, so exact keys can compare
// and hash the stored numerator/denominator without Ratio's repeated divisions.
impl PartialEq for Rational {
    fn eq(&self, rhs: &Self) -> bool {
        match (&self.0, &rhs.0) {
            (
                Representation::Small {
                    numerator: a,
                    denominator: b,
                },
                Representation::Small {
                    numerator: c,
                    denominator: d,
                },
            ) => a == c && b == d,
            (Representation::Big(a), Representation::Big(b)) => {
                Arc::ptr_eq(a, b) || (a.numer() == b.numer() && a.denom() == b.denom())
            }
            _ => false,
        }
    }
}
impl Eq for Rational {}
impl Hash for Rational {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match &self.0 {
            Representation::Small {
                numerator,
                denominator,
            } => {
                Self::hash_small(*numerator, *denominator, state);
            }
            Representation::Big(value) => {
                1_u8.hash(state);
                value.numer().hash(state);
                value.denom().hash(state);
            }
        }
    }
}

impl Ord for Rational {
    fn cmp(&self, rhs: &Self) -> Ordering {
        match (self.as_small(), rhs.as_small()) {
            (Some((a, b)), Some((c, d))) => (a as i128 * d as i128).cmp(&(c as i128 * b as i128)),
            _ => self.to_big_rational().cmp(&rhs.to_big_rational()),
        }
    }
}
impl PartialOrd for Rational {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}
impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Representation::Small {
                numerator,
                denominator: 1,
            } => write!(f, "{numerator}"),
            Representation::Small {
                numerator,
                denominator,
            } => write!(f, "({numerator}/{denominator})"),
            Representation::Big(value) if value.is_integer() => write!(f, "{}", value.numer()),
            Representation::Big(value) => write!(f, "({}/{})", value.numer(), value.denom()),
        }
    }
}
