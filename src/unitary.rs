use std::f64::consts::{FRAC_1_SQRT_2, PI};

use crate::circuit::{Circuit, Gate};

#[derive(Clone, Copy, Debug)]
pub(crate) struct C {
    re: f64,
    im: f64,
}

impl C {
    const ZERO: C = C { re: 0.0, im: 0.0 };
    const ONE: C = C { re: 1.0, im: 0.0 };

    fn new(re: f64, im: f64) -> C {
        C { re, im }
    }

    fn polar(r: f64, theta: f64) -> C {
        C { re: r * theta.cos(), im: r * theta.sin() }
    }

    fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    fn conj(self) -> C {
        C { re: self.re, im: -self.im }
    }
}

impl std::ops::Add for C {
    type Output = C;
    fn add(self, rhs: C) -> C {
        C { re: self.re + rhs.re, im: self.im + rhs.im }
    }
}

impl std::ops::Sub for C {
    type Output = C;
    fn sub(self, rhs: C) -> C {
        C { re: self.re - rhs.re, im: self.im - rhs.im }
    }
}

impl std::ops::Mul for C {
    type Output = C;
    fn mul(self, rhs: C) -> C {
        C {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl std::ops::Mul<f64> for C {
    type Output = C;
    fn mul(self, rhs: f64) -> C {
        C { re: self.re * rhs, im: self.im * rhs }
    }
}

/// Dense 2^n × 2^n complex matrix stored in row-major order.
struct Mat {
    dim: usize,
    data: Vec<C>,
}

impl Mat {
    fn identity(num_qubits: usize) -> Mat {
        let dim = 1 << num_qubits;
        let mut data = vec![C::ZERO; dim * dim];
        for i in 0..dim {
            data[i * dim + i] = C::ONE;
        }
        Mat { dim, data }
    }

    fn get(&self, row: usize, col: usize) -> C {
        self.data[row * self.dim + col]
    }

    fn set(&mut self, row: usize, col: usize, v: C) {
        self.data[row * self.dim + col] = v;
    }

    /// Apply a single-qubit gate to every column of the matrix.
    fn apply_single(&mut self, g: [[C; 2]; 2], q: usize, num_qubits: usize) {
        let bit = 1 << (num_qubits - 1 - q);
        for col in 0..self.dim {
            for row in 0..self.dim {
                if row & bit != 0 {
                    continue;
                }
                let r0 = row;
                let r1 = row | bit;
                let a = self.get(r0, col);
                let b = self.get(r1, col);
                self.set(r0, col, g[0][0] * a + g[0][1] * b);
                self.set(r1, col, g[1][0] * a + g[1][1] * b);
            }
        }
    }

    /// Apply a controlled-NOT: flip target bit when control bit is set.
    fn apply_cnot(&mut self, control: usize, target: usize, num_qubits: usize) {
        let cb = 1 << (num_qubits - 1 - control);
        let tb = 1 << (num_qubits - 1 - target);
        for col in 0..self.dim {
            for row in 0..self.dim {
                if row & cb != 0 && row & tb == 0 {
                    let r1 = row | tb;
                    let a = self.get(row, col);
                    let b = self.get(r1, col);
                    self.set(row, col, b);
                    self.set(r1, col, a);
                }
            }
        }
    }

    /// Apply a Toffoli: flip target bit when both control bits are set.
    fn apply_ccx(&mut self, c1: usize, c2: usize, target: usize, num_qubits: usize) {
        let c1b = 1 << (num_qubits - 1 - c1);
        let c2b = 1 << (num_qubits - 1 - c2);
        let tb = 1 << (num_qubits - 1 - target);
        for col in 0..self.dim {
            for row in 0..self.dim {
                if row & c1b != 0 && row & c2b != 0 && row & tb == 0 {
                    let r1 = row | tb;
                    let a = self.get(row, col);
                    let b = self.get(r1, col);
                    self.set(row, col, b);
                    self.set(r1, col, a);
                }
            }
        }
    }
}

fn gate_matrix_x() -> [[C; 2]; 2] {
    [[C::ZERO, C::ONE], [C::ONE, C::ZERO]]
}

fn gate_matrix_h() -> [[C; 2]; 2] {
    let v = FRAC_1_SQRT_2;
    let p = C::new(v, 0.0);
    let m = C::new(-v, 0.0);
    [[p, p], [p, m]]
}

fn gate_matrix_s() -> [[C; 2]; 2] {
    [[C::ONE, C::ZERO], [C::ZERO, C::new(0.0, 1.0)]]
}

fn gate_matrix_sdg() -> [[C; 2]; 2] {
    [[C::ONE, C::ZERO], [C::ZERO, C::new(0.0, -1.0)]]
}

fn gate_matrix_z() -> [[C; 2]; 2] {
    [[C::ONE, C::ZERO], [C::ZERO, C::new(-1.0, 0.0)]]
}

fn gate_matrix_t() -> [[C; 2]; 2] {
    [[C::ONE, C::ZERO], [C::ZERO, C::polar(1.0, PI / 4.0)]]
}

fn gate_matrix_tdg() -> [[C; 2]; 2] {
    [[C::ONE, C::ZERO], [C::ZERO, C::polar(1.0, -PI / 4.0)]]
}

fn gate_matrix_rz(theta: f64) -> [[C; 2]; 2] {
    [[C::polar(1.0, -theta / 2.0), C::ZERO], [C::ZERO, C::polar(1.0, theta / 2.0)]]
}

pub(crate) fn circuit_unitary(circuit: &Circuit) -> Vec<Vec<C>> {
    let n = circuit.num_qubits;
    let mut mat = Mat::identity(n);
    for gate in &circuit.gates {
        match gate {
            Gate::x(q) => mat.apply_single(gate_matrix_x(), *q, n),
            Gate::h(q) => mat.apply_single(gate_matrix_h(), *q, n),
            Gate::s(q) => mat.apply_single(gate_matrix_s(), *q, n),
            Gate::sdg(q) => mat.apply_single(gate_matrix_sdg(), *q, n),
            Gate::z(q) => mat.apply_single(gate_matrix_z(), *q, n),
            Gate::t(q) => mat.apply_single(gate_matrix_t(), *q, n),
            Gate::tdg(q) => mat.apply_single(gate_matrix_tdg(), *q, n),
            Gate::rz(theta, q) => mat.apply_single(gate_matrix_rz(*theta), *q, n),
            Gate::cnot { control, target } => mat.apply_cnot(*control, *target, n),
            Gate::ccx { control1, control2, target } => {
                mat.apply_ccx(*control1, *control2, *target, n)
            }
        }
    }
    // convert to Vec<Vec<C>> for external use
    (0..mat.dim)
        .map(|r| (0..mat.dim).map(|c| mat.get(r, c)).collect())
        .collect()
}


/// Check if two unitaries are equal up to a global phase.
/// Finds the phase from the first nonzero entry, then checks all entries.
pub fn circuits_equiv(a: &Circuit, b: &Circuit, tol: f64) -> bool {
    assert_eq!(a.num_qubits, b.num_qubits);
    let ua = circuit_unitary(a);
    let ub = circuit_unitary(b);
    let dim = ua.len();

    // find global phase: first entry where both are nonzero
    // phase = conj(ua[r][c]) * ub[r][c] / |ua[r][c]|^2
    let mut phase = None;
    for r in 0..dim {
        for c in 0..dim {
            let nsq = ua[r][c].norm_sq();
            if nsq > tol * tol {
                // phase such that ub = phase * ua
                let p = ub[r][c] * ua[r][c].conj() * (1.0 / nsq);
                phase = Some(p);
                break;
            }
        }
        if phase.is_some() {
            break;
        }
    }

    let phase = match phase {
        Some(p) => p,
        None => return true, // both zero matrices
    };

    // check: ub[r][c] ≈ phase * ua[r][c] for all r,c
    for r in 0..dim {
        for c in 0..dim {
            let expected = phase * ua[r][c];
            let diff = ub[r][c] - expected;
            if diff.norm_sq() > tol * tol {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::Gate;

    #[test]
    fn identity_is_identity() {
        let c = Circuit::new(1);
        let u = circuit_unitary(&c);
        assert!((u[0][0].re - 1.0).abs() < 1e-10);
        assert!((u[1][1].re - 1.0).abs() < 1e-10);
        assert!(u[0][1].norm_sq() < 1e-20);
        assert!(u[1][0].norm_sq() < 1e-20);
    }

    #[test]
    fn x_squared_is_identity() {
        let mut a = Circuit::new(1);
        a.apply(Gate::x(0));
        a.apply(Gate::x(0));
        let b = Circuit::new(1);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn h_squared_is_identity() {
        let mut a = Circuit::new(1);
        a.apply(Gate::h(0));
        a.apply(Gate::h(0));
        let b = Circuit::new(1);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn t_squared_is_s() {
        let mut a = Circuit::new(1);
        a.apply(Gate::t(0));
        a.apply(Gate::t(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::s(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn t_tdg_is_identity() {
        let mut a = Circuit::new(1);
        a.apply(Gate::t(0));
        a.apply(Gate::tdg(0));
        let b = Circuit::new(1);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn s_is_rz_pi_over_2() {
        let mut a = Circuit::new(1);
        a.apply(Gate::s(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::rz(PI / 2.0, 0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn cnot_swap_decomposition() {
        // SWAP = CNOT(0,1) CNOT(1,0) CNOT(0,1)
        let mut a = Circuit::new(2);
        a.apply(Gate::cnot { control: 0, target: 1 });
        a.apply(Gate::cnot { control: 1, target: 0 });
        a.apply(Gate::cnot { control: 0, target: 1 });

        // build SWAP directly: |00⟩→|00⟩, |01⟩→|10⟩, |10⟩→|01⟩, |11⟩→|11⟩
        let u = circuit_unitary(&a);
        assert!((u[0][0].re - 1.0).abs() < 1e-10);
        assert!((u[1][2].re - 1.0).abs() < 1e-10);
        assert!((u[2][1].re - 1.0).abs() < 1e-10);
        assert!((u[3][3].re - 1.0).abs() < 1e-10);
    }

    #[test]
    fn hxh_is_z() {
        // HXH = Z
        let mut a = Circuit::new(1);
        a.apply(Gate::h(0));
        a.apply(Gate::x(0));
        a.apply(Gate::h(0));
        // Z = S^2 = Rz(pi) up to global phase
        let mut b = Circuit::new(1);
        b.apply(Gate::s(0));
        b.apply(Gate::s(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn not_equiv_h_vs_x() {
        let mut a = Circuit::new(1);
        a.apply(Gate::h(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::x(0));
        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn ccx_equiv_decomposition() {
        // Standard Toffoli decomposition into 1- and 2-qubit gates
        let mut a = Circuit::new(3);
        a.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });

        let mut b = Circuit::new(3);
        b.apply(Gate::h(2));
        b.apply(Gate::cnot { control: 1, target: 2 });
        b.apply(Gate::tdg(2));
        b.apply(Gate::cnot { control: 0, target: 2 });
        b.apply(Gate::t(2));
        b.apply(Gate::cnot { control: 1, target: 2 });
        b.apply(Gate::tdg(2));
        b.apply(Gate::cnot { control: 0, target: 2 });
        b.apply(Gate::t(1));
        b.apply(Gate::t(2));
        b.apply(Gate::h(2));
        b.apply(Gate::cnot { control: 0, target: 1 });
        b.apply(Gate::t(0));
        b.apply(Gate::tdg(1));
        b.apply(Gate::cnot { control: 0, target: 1 });

        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn ghz_3qubit_two_ways() {
        // Way 1: H(0), CNOT(0,1), CNOT(1,2)
        let mut a = Circuit::new(3);
        a.apply(Gate::h(0));
        a.apply(Gate::cnot { control: 0, target: 1 });
        a.apply(Gate::cnot { control: 1, target: 2 });

        // Way 2: H(0), CNOT(0,1), CNOT(0,2) — fan-out from q0
        let mut b = Circuit::new(3);
        b.apply(Gate::h(0));
        b.apply(Gate::cnot { control: 0, target: 1 });
        b.apply(Gate::cnot { control: 0, target: 2 });

        // These are NOT equivalent unitaries (same final state on |000⟩ but different maps)
        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn cnot_commute_disjoint_3qubit() {
        // CNOT(0,1) then CNOT(2,1) vs reversed order — disjoint controls, same target
        // These do NOT commute (both write to q1)
        let mut a = Circuit::new(3);
        a.apply(Gate::cnot { control: 0, target: 1 });
        a.apply(Gate::cnot { control: 2, target: 1 });

        let mut b = Circuit::new(3);
        b.apply(Gate::cnot { control: 2, target: 1 });
        b.apply(Gate::cnot { control: 0, target: 1 });

        // Actually XOR is commutative+associative, so these ARE equivalent
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn cnot_no_commute_overlapping() {
        // CNOT(0,1) then CNOT(1,2) vs reversed — these don't commute
        let mut a = Circuit::new(3);
        a.apply(Gate::cnot { control: 0, target: 1 });
        a.apply(Gate::cnot { control: 1, target: 2 });

        let mut b = Circuit::new(3);
        b.apply(Gate::cnot { control: 1, target: 2 });
        b.apply(Gate::cnot { control: 0, target: 1 });

        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn swap_3qubit_cycle() {
        // SWAP(0,1) then SWAP(1,2) = cyclic permutation 0→1→2→0
        // vs SWAP(1,2) then SWAP(0,1) = cyclic permutation 0→2→1→0
        fn swap(c: &mut Circuit, a: usize, b: usize) {
            c.apply(Gate::cnot { control: a, target: b });
            c.apply(Gate::cnot { control: b, target: a });
            c.apply(Gate::cnot { control: a, target: b });
        }

        let mut a = Circuit::new(3);
        swap(&mut a, 0, 1);
        swap(&mut a, 1, 2);

        let mut b = Circuit::new(3);
        swap(&mut b, 1, 2);
        swap(&mut b, 0, 1);

        // different cyclic directions — not equivalent
        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn double_ccx_is_identity() {
        // Toffoli is self-inverse
        let mut a = Circuit::new(3);
        a.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        a.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        let b = Circuit::new(3);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn t_phase_folding_3qubit() {
        // T on q0, CNOT(0,1), T on q1 is NOT the same as S on q0, CNOT(0,1)
        let mut a = Circuit::new(3);
        a.apply(Gate::t(0));
        a.apply(Gate::cnot { control: 0, target: 1 });
        a.apply(Gate::t(1));

        let mut b = Circuit::new(3);
        b.apply(Gate::s(0));
        b.apply(Gate::cnot { control: 0, target: 1 });

        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn four_qubit_parallel_cnots() {
        // CNOT(0,1) and CNOT(2,3) in either order — fully disjoint, must commute
        let mut a = Circuit::new(4);
        a.apply(Gate::cnot { control: 0, target: 1 });
        a.apply(Gate::cnot { control: 2, target: 3 });

        let mut b = Circuit::new(4);
        b.apply(Gate::cnot { control: 2, target: 3 });
        b.apply(Gate::cnot { control: 0, target: 1 });

        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn four_qubit_h_layer_order() {
        // H on all 4 qubits — order doesn't matter
        let mut a = Circuit::new(4);
        for q in 0..4 { a.apply(Gate::h(q)); }

        let mut b = Circuit::new(4);
        for q in (0..4).rev() { b.apply(Gate::h(q)); }

        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn four_qubit_cnot_ladder_not_reversible() {
        // CNOT ladder 0→1→2→3 vs 3→2→1→0 — different circuits
        let mut a = Circuit::new(4);
        for q in 0..3 {
            a.apply(Gate::cnot { control: q, target: q + 1 });
        }

        let mut b = Circuit::new(4);
        for q in (0..3).rev() {
            b.apply(Gate::cnot { control: q, target: q + 1 });
        }

        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn z_is_s_squared() {
        let mut a = Circuit::new(1);
        a.apply(Gate::z(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::s(0));
        b.apply(Gate::s(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn z_is_rz_pi() {
        let mut a = Circuit::new(1);
        a.apply(Gate::z(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::rz(PI, 0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn sdg_is_rz_neg_half_pi() {
        let mut a = Circuit::new(1);
        a.apply(Gate::sdg(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::rz(-PI / 2.0, 0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn s_sdg_is_identity() {
        let mut a = Circuit::new(1);
        a.apply(Gate::s(0));
        a.apply(Gate::sdg(0));
        let b = Circuit::new(1);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn sdg_s_is_identity() {
        let mut a = Circuit::new(1);
        a.apply(Gate::sdg(0));
        a.apply(Gate::s(0));
        let b = Circuit::new(1);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn z_z_is_identity() {
        let mut a = Circuit::new(1);
        a.apply(Gate::z(0));
        a.apply(Gate::z(0));
        let b = Circuit::new(1);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn z_is_self_inverse() {
        let mut a = Circuit::new(1);
        a.apply(Gate::z(0));
        a.apply(Gate::z(0));
        let b = Circuit::new(1);
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn sdg_sdg_is_z() {
        let mut a = Circuit::new(1);
        a.apply(Gate::sdg(0));
        a.apply(Gate::sdg(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::z(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn hzh_is_x() {
        // HZH = X
        let mut a = Circuit::new(1);
        a.apply(Gate::h(0));
        a.apply(Gate::z(0));
        a.apply(Gate::h(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::x(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn hsdgh_is_not_s() {
        // HSdgH ≠ S (different from HZH=X pattern)
        let mut a = Circuit::new(1);
        a.apply(Gate::h(0));
        a.apply(Gate::sdg(0));
        a.apply(Gate::h(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::s(0));
        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn z_commutes_with_cnot_control() {
        // Z on control commutes with CNOT
        let mut a = Circuit::new(2);
        a.apply(Gate::z(0));
        a.apply(Gate::cnot { control: 0, target: 1 });
        let mut b = Circuit::new(2);
        b.apply(Gate::cnot { control: 0, target: 1 });
        b.apply(Gate::z(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn z_does_not_commute_with_cnot_target() {
        // Z on target does NOT commute with CNOT
        let mut a = Circuit::new(2);
        a.apply(Gate::z(1));
        a.apply(Gate::cnot { control: 0, target: 1 });
        let mut b = Circuit::new(2);
        b.apply(Gate::cnot { control: 0, target: 1 });
        b.apply(Gate::z(1));
        assert!(!circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn sdg_three_times_is_s() {
        let mut a = Circuit::new(1);
        a.apply(Gate::sdg(0));
        a.apply(Gate::sdg(0));
        a.apply(Gate::sdg(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::s(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn s_three_times_is_sdg() {
        let mut a = Circuit::new(1);
        a.apply(Gate::s(0));
        a.apply(Gate::s(0));
        a.apply(Gate::s(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::sdg(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn six_t_is_sdg() {
        // 6T = 6*(π/4) = 3π/2 = Sdg
        let mut a = Circuit::new(1);
        for _ in 0..6 { a.apply(Gate::t(0)); }
        let mut b = Circuit::new(1);
        b.apply(Gate::sdg(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn z_t_is_five_quarter_pi() {
        // Z + T = Rz(π) + Rz(π/4) = Rz(5π/4)
        let mut a = Circuit::new(1);
        a.apply(Gate::z(0));
        a.apply(Gate::t(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::rz(5.0 * PI / 4.0, 0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn sdg_t_is_tdg() {
        // Sdg + T = -π/2 + π/4 = -π/4 = Tdg
        let mut a = Circuit::new(1);
        a.apply(Gate::sdg(0));
        a.apply(Gate::t(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::tdg(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }

    #[test]
    fn s_tdg_is_t() {
        // S + Tdg = π/2 - π/4 = π/4 = T
        let mut a = Circuit::new(1);
        a.apply(Gate::s(0));
        a.apply(Gate::tdg(0));
        let mut b = Circuit::new(1);
        b.apply(Gate::t(0));
        assert!(circuits_equiv(&a, &b, 1e-10));
    }
}
