//! A Clifford accumulated from Pauli rotations by pi/4 and pi/2, as the
//! images of every qubit's X, Y, and Z, for moving Clifford rotations into
//! the output frame.

use super::words::{anticommutes, product};

/// The Clifford F, as `F† P_q F` for P in X, Y, Z and every qubit q: packed
/// canonical words with phase exponents (`i^p W`). Starts as the identity.
pub(super) struct Clifford {
    l: usize,
    /// Image of X_q at index 3q, Y_q at 3q + 1, Z_q at 3q + 2; `2 l` words each.
    images: Vec<u64>,
    phases: Vec<u8>,
    /// Support signature of each X and Z image (as in `Axes::support`), to
    /// skip images that cannot anticommute with an absorbed axis.
    supports: Vec<u64>,
    /// Qubits whose images may differ from the identity's.
    touched: Vec<u64>,
    signature: fn(&[u64], usize) -> u64,
    scratch: Vec<u64>,
}

const LETTERS: usize = 3;

impl Clifford {
    pub fn identity(num_qubits: usize, l: usize, signature: fn(&[u64], usize) -> u64) -> Self {
        let stride = 2 * l;
        let mut images = vec![0u64; LETTERS * num_qubits * stride];
        let mut supports = vec![0u64; LETTERS * num_qubits];
        for q in 0..num_qubits {
            let (word, bit) = (q / 64, 1u64 << (q % 64));
            for (letter, (x, z)) in [(true, false), (true, true), (false, true)]
                .into_iter()
                .enumerate()
            {
                let start = (LETTERS * q + letter) * stride;
                if x {
                    images[start + word] = bit;
                }
                if z {
                    images[start + l + word] = bit;
                }
                supports[LETTERS * q + letter] = signature(&images[start..start + stride], l);
            }
        }
        Self {
            l,
            images,
            phases: vec![0; LETTERS * num_qubits],
            supports,
            touched: vec![0; l],
            signature,
            scratch: vec![0; stride],
        }
    }

    /// Compose a Pauli rotation `exp(-i k pi/8 B)`, k in {2, -2, 4}, executed
    /// before F: `F <- F R`. Each image I becomes `R† I R`, which for an image
    /// anticommuting with B is `i B I` (k = 2), `-i B I` (k = -2), or `-I`
    /// (k = 4); commuting images are unchanged. Y images follow the same rule,
    /// since conjugation is linear.
    pub fn absorb(&mut self, b: &[u64], b_support: u64, k: i8) {
        debug_assert!(matches!(k, 2 | -2 | 4));
        let (l, stride) = (self.l, 2 * self.l);
        for index in 0..self.phases.len() {
            if self.supports[index] & b_support == 0 {
                continue;
            }
            let start = index * stride;
            let image = &mut self.images[start..start + stride];
            if !anticommutes(image, b, l) {
                continue;
            }
            if k == 4 {
                self.phases[index] = (self.phases[index] + 2) % 4;
            } else {
                let e = product(b, image, &mut self.scratch, l);
                let turn = if k == 2 { 1 } else { 3 };
                self.phases[index] = (self.phases[index] + e + turn) % 4;
                image.copy_from_slice(&self.scratch);
                self.supports[index] = (self.signature)(image, l);
            }
            let q = index / LETTERS;
            self.touched[q / 64] |= 1 << (q % 64);
        }
    }

    /// Whether `F† W F = W` is certain because W acts on no touched qubit.
    pub fn fixes(&self, w: &[u64]) -> bool {
        let (x, z) = w.split_at(self.l);
        (0..self.l).all(|j| (x[j] | z[j]) & self.touched[j] == 0)
    }

    /// `F† W F` for a canonical word W, written to `out`; returns its phase
    /// exponent. W is the product of its single-qubit factors, whose images
    /// commute with each other, and untouched qubits map to themselves.
    pub fn conjugate(&mut self, w: &[u64], out: &mut [u64]) -> u8 {
        let (l, stride) = (self.l, 2 * self.l);
        let (x, z) = w.split_at(l);
        out.copy_from_slice(w);
        for j in 0..l {
            out[j] &= !self.touched[j];
            out[l + j] &= !self.touched[j];
        }
        let mut phase = 0u32;
        for j in 0..l {
            let mut moved = (x[j] | z[j]) & self.touched[j];
            while moved != 0 {
                let bit = moved.trailing_zeros() as usize;
                moved &= moved - 1;
                let letter = match (x[j] >> bit & 1, z[j] >> bit & 1) {
                    (1, 0) => 0,
                    (1, 1) => 1,
                    _ => 2,
                };
                let index = LETTERS * (64 * j + bit) + letter;
                let image = &self.images[index * stride..(index + 1) * stride];
                let e = product(out, image, &mut self.scratch, l);
                out.copy_from_slice(&self.scratch);
                phase += u32::from(self.phases[index]) + u32::from(e);
            }
        }
        (phase % 4) as u8
    }
}
