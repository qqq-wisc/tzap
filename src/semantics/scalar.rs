use num_bigint::BigInt;
use num_rational::BigRational;

/// Exact coefficients of 1, sqrt(2), i, i*sqrt(2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Scalar(pub [BigRational; 4]);

impl Scalar {
    pub fn integer(n: i64) -> Self {
        let mut result = Self::zero();
        result.0[0] = BigRational::from_integer(BigInt::from(n));
        result
    }
    pub fn zero() -> Self {
        Self(std::array::from_fn(|_| {
            BigRational::from_integer(BigInt::from(0))
        }))
    }
    pub fn i() -> Self {
        let mut result = Self::zero();
        result.0[2] = BigRational::from_integer(BigInt::from(1));
        result
    }
    pub fn inv_sqrt_two() -> Self {
        let mut result = Self::zero();
        result.0[1] = BigRational::new(BigInt::from(1), BigInt::from(2));
        result
    }
    pub fn omega(k: u8) -> Self {
        let omega = Self::integer(1).add(&Self::i()).mul(&Self::inv_sqrt_two());
        (0..k % 8).fold(Self::integer(1), |a, _| a.mul(&omega))
    }
    pub fn add(&self, rhs: &Self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] + &rhs.0[i]))
    }
    pub fn neg(&self) -> Self {
        Self(std::array::from_fn(|i| -&self.0[i]))
    }
    pub fn half(&self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] / BigInt::from(2)))
    }
    pub fn conj(&self) -> Self {
        let mut result = self.clone();
        result.0[2] = -&result.0[2];
        result.0[3] = -&result.0[3];
        result
    }
    pub fn mul(&self, rhs: &Self) -> Self {
        let mut result = Self::zero();
        for a in 0..4 {
            for b in 0..4 {
                // Low bit denotes sqrt(2); high bit denotes i.
                let factor =
                    (if a & b & 1 != 0 { 2 } else { 1 }) * (if a & b & 2 != 0 { -1 } else { 1 });
                result.0[a ^ b] += &self.0[a] * &rhs.0[b] * BigInt::from(factor);
            }
        }
        result
    }
}
