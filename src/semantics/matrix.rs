use super::scalar::Scalar;

/// Row-major matrix; q0 is the most significant computational-basis bit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Matrix {
    pub dim: usize,
    pub entries: Vec<Scalar>,
}

impl Matrix {
    pub fn zero(dim: usize) -> Self {
        Self {
            dim,
            entries: vec![Scalar::zero(); dim * dim],
        }
    }
    pub fn identity(dim: usize) -> Self {
        let mut result = Self::zero(dim);
        for i in 0..dim {
            result.entries[i * dim + i] = Scalar::integer(1);
        }
        result
    }
    pub fn get(&self, row: usize, col: usize) -> &Scalar {
        &self.entries[row * self.dim + col]
    }
    pub fn set(&mut self, row: usize, col: usize, value: Scalar) {
        self.entries[row * self.dim + col] = value;
    }
    pub fn scale(&self, factor: &Scalar) -> Self {
        Self {
            dim: self.dim,
            entries: self.entries.iter().map(|x| x.mul(factor)).collect(),
        }
    }
    pub fn add(&self, rhs: &Self) -> Self {
        assert_eq!(self.dim, rhs.dim);
        Self {
            dim: self.dim,
            entries: self
                .entries
                .iter()
                .zip(&rhs.entries)
                .map(|(a, b)| a.add(b))
                .collect(),
        }
    }
    pub fn mul(&self, rhs: &Self) -> Self {
        assert_eq!(self.dim, rhs.dim);
        let mut result = Self::zero(self.dim);
        let zero = Scalar::zero();
        for r in 0..self.dim {
            for k in 0..self.dim {
                if self.get(r, k) == &zero {
                    continue;
                }
                for c in 0..self.dim {
                    if rhs.get(k, c) == &zero {
                        continue;
                    }
                    let index = r * self.dim + c;
                    result.entries[index] =
                        result.entries[index].add(&self.get(r, k).mul(rhs.get(k, c)));
                }
            }
        }
        result
    }
    pub fn adjoint(&self) -> Self {
        let mut result = Self::zero(self.dim);
        for r in 0..self.dim {
            for c in 0..self.dim {
                result.set(r, c, self.get(c, r).conj());
            }
        }
        result
    }
    /// Exact comparison without division: test a common unit-modulus phase
    /// using cross-products against a nonzero pivot.
    pub fn equivalent_up_to_global_phase(&self, rhs: &Self) -> bool {
        if self.dim != rhs.dim {
            return false;
        }
        let zero = Scalar::zero();
        let Some(i) = self.entries.iter().position(|x| x != &zero) else {
            return rhs.entries.iter().all(|x| x == &zero);
        };
        let a = &self.entries[i];
        let b = &rhs.entries[i];
        a.mul(&a.conj()) == b.mul(&b.conj())
            && self
                .entries
                .iter()
                .zip(&rhs.entries)
                .all(|(x, y)| x.mul(b) == y.mul(a))
    }
}
