use std::collections::BTreeMap;
use std::f64::consts::PI;

use indicatif::ProgressBar;

use crate::circuit::{Circuit, Gate};
use crate::pass::Pass;

/// A parity expression: XOR of a set of variables, optionally complemented.
/// Two expressions are equivalent iff they have the same variable set and
/// the same `negated` flag.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct ParityExpr {
    /// Sorted variable indices whose XOR forms this expression.
    vars: Vec<usize>,
    /// Whether the expression is bitwise-complemented (from X gates).
    negated: bool,
}

impl ParityExpr {
    fn new(var: usize) -> Self {
        ParityExpr { vars: vec![var], negated: false }
    }

    fn fresh(var: usize) -> Self {
        Self::new(var)
    }

    /// XOR two parity expressions (symmetric difference of variable sets).
    fn xor(&self, other: &ParityExpr) -> ParityExpr {
        let mut result = Vec::new();
        let (a, b) = (&self.vars, &other.vars);
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            if a[i] < b[j] {
                result.push(a[i]);
                i += 1;
            } else if a[i] > b[j] {
                result.push(b[j]);
                j += 1;
            } else {
                // Same variable in both — cancels out in XOR
                i += 1;
                j += 1;
            }
        }
        result.extend_from_slice(&a[i..]);
        result.extend_from_slice(&b[j..]);
        ParityExpr {
            vars: result,
            negated: self.negated ^ other.negated,
        }
    }

    /// Complement (NOT) — flips the negated flag.
    fn complement(&self) -> ParityExpr {
        ParityExpr {
            vars: self.vars.clone(),
            negated: !self.negated,
        }
    }

}

struct Group {
    angle: f64,
    last_idx: usize,
    last_qubit: usize,
    indices: Vec<usize>,
}

pub struct PhaseFoldGlobalExpr;

impl Pass for PhaseFoldGlobalExpr {
    fn name(&self) -> &str { "Phase folding (expr)" }
    fn run_with_progress(&self, circuit: &Circuit, pb: &ProgressBar) -> Circuit {
        phase_fold_global_expr(circuit, pb)
    }
}

pub fn phase_fold_global_expr(circuit: &Circuit, pb: &ProgressBar) -> Circuit {
    let n = circuit.num_qubits;
    let mut next_var = 0usize;
    let mut fresh = || { let v = next_var; next_var += 1; ParityExpr::fresh(v) };
    let mut qubits: Vec<ParityExpr> = (0..n).map(|_| fresh()).collect();

    let mut groups: BTreeMap<ParityExpr, Group> = BTreeMap::new();
    let mut skip = vec![false; circuit.gates.len()];
    let mut emit_at: Vec<Option<(usize, f64)>> = vec![None; circuit.gates.len()];

    for (idx, gate) in circuit.gates.iter().enumerate() {
        if idx & 0xFFF == 0 { pb.inc(0x1000); }
        match gate {
            Gate::t(q) => record_phase(&qubits, *q, PI / 4.0, idx, &mut groups),
            Gate::tdg(q) => record_phase(&qubits, *q, -PI / 4.0, idx, &mut groups),
            Gate::s(q) => record_phase(&qubits, *q, PI / 2.0, idx, &mut groups),
            Gate::sdg(q) => record_phase(&qubits, *q, -PI / 2.0, idx, &mut groups),
            Gate::z(q) => record_phase(&qubits, *q, PI, idx, &mut groups),
            Gate::rz(theta, q) => record_phase(&qubits, *q, *theta, idx, &mut groups),
            Gate::h(q) => { qubits[*q] = fresh(); }
            Gate::x(q) => { qubits[*q] = qubits[*q].complement(); }
            Gate::cnot { control, target } => {
                let ctrl = qubits[*control].clone();
                qubits[*target] = qubits[*target].xor(&ctrl);
            }
            Gate::ccx { control1: _, control2: _, target } => {
                // AND-based parity — opaque, just refresh the target.
                qubits[*target] = fresh();
            }
        }
    }

    // Finalize skip/emit from groups.
    for g in groups.values() {
        for &i in &g.indices {
            skip[i] = true;
        }
        if !angle_is_zero(g.angle) {
            skip[g.last_idx] = false;
            emit_at[g.last_idx] = Some((g.last_qubit, g.angle));
        }
    }

    // Reconstruct circuit.
    let mut output = Circuit::new(n);
    for (idx, gate) in circuit.gates.iter().enumerate() {
        if skip[idx] {
            continue;
        }
        if let Some((qubit, angle)) = emit_at[idx] {
            emit_rotation(&mut output, qubit, angle);
        } else {
            output.apply(gate.clone());
        }
    }

    output
}

fn record_phase(
    qubits: &[ParityExpr],
    q: usize,
    angle: f64,
    idx: usize,
    groups: &mut BTreeMap<ParityExpr, Group>,
) {
    let parity = &qubits[q];

    if let Some(g) = groups.get_mut(parity) {
        g.angle += angle;
        g.last_idx = idx;
        g.last_qubit = q;
        g.indices.push(idx);
        return;
    }

    groups.insert(parity.clone(), Group {
        angle,
        last_idx: idx,
        last_qubit: q,
        indices: vec![idx],
    });
}

fn angle_is_zero(angle: f64) -> bool {
    let n = angle.rem_euclid(2.0 * PI);
    n < 1e-6 || (2.0 * PI - n) < 1e-6
}

fn emit_rotation(output: &mut Circuit, q: usize, angle: f64) {
    let n = angle.rem_euclid(2.0 * PI);
    if angle_is_zero(n) {
        return;
    }
    let quarter = PI / 4.0;
    let k = (n / quarter).round();
    // Use 1e-6 tolerance to handle floating-point drift from summing many π/4 multiples.
    if (n - k * quarter).abs() < 1e-6 {
        let k = k as u32 % 8;
        match k {
            0 => {}
            1 => output.apply(Gate::t(q)),
            2 => output.apply(Gate::s(q)),
            3 => { output.apply(Gate::s(q)); output.apply(Gate::t(q)); }
            4 => output.apply(Gate::z(q)),
            5 => { output.apply(Gate::z(q)); output.apply(Gate::t(q)); }
            6 => output.apply(Gate::sdg(q)),
            7 => output.apply(Gate::tdg(q)),
            _ => unreachable!(),
        }
    } else {
        output.apply(Gate::rz(n, q));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unitary::circuits_equiv;

    const TOL: f64 = 1e-10;

    fn count_phase_gates(c: &Circuit) -> usize {
        c.gates.iter().filter(|g| matches!(g,
            Gate::t(_) | Gate::tdg(_) | Gate::s(_) | Gate::sdg(_) | Gate::z(_) | Gate::rz(..)
        )).count()
    }

    fn count_t_gates(c: &Circuit) -> usize {
        c.gates.iter().filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count()
    }

    fn run(c: &Circuit) -> Circuit {
        phase_fold_global_expr(c, &ProgressBar::hidden())
    }

    // --- Expr equivalence tests ---

    #[test]
    fn expr_fresh_not_equal() {
        let a = ParityExpr::fresh(0);
        let b = ParityExpr::fresh(1);
        assert_ne!(a, b);
    }

    #[test]
    fn expr_same_var_equal() {
        let a = ParityExpr::fresh(0);
        let b = ParityExpr::fresh(0);
        assert_eq!(a, b);
    }

    #[test]
    fn expr_xor_self_is_empty() {
        let a = ParityExpr::fresh(0);
        let r = a.xor(&a);
        assert!(r.vars.is_empty());
        assert!(!r.negated);
    }

    #[test]
    fn expr_xor_commutative() {
        let a = ParityExpr::fresh(0);
        let b = ParityExpr::fresh(1);
        let ab = a.xor(&b);
        let ba = b.xor(&a);
        assert_eq!(ab, ba);
    }

    #[test]
    fn expr_complement_not_equal() {
        let a = ParityExpr::fresh(0);
        let b = a.complement();
        assert_ne!(a, b);
    }

    #[test]
    fn expr_double_complement_equal() {
        let a = ParityExpr::fresh(0);
        let b = a.complement().complement();
        assert_eq!(a, b);
    }

    #[test]
    fn expr_xor_associative() {
        let a = ParityExpr::fresh(0);
        let b = ParityExpr::fresh(1);
        let c = ParityExpr::fresh(2);
        let ab_c = a.xor(&b).xor(&c);
        let a_bc = a.xor(&b.xor(&c));
        assert_eq!(ab_c, a_bc);
    }

    // --- Circuit-level tests ---

    #[test]
    fn two_t_merge_to_s() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::t(0));
        let out = run(&c);
        assert_eq!(count_phase_gates(&out), 1);
        assert_eq!(count_t_gates(&out), 0);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn t_tdg_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        let out = run(&c);
        assert_eq!(count_phase_gates(&out), 0);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn s_s_merge_to_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::s(0));
        let out = run(&c);
        assert_eq!(count_phase_gates(&out), 1);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn four_t_merge_to_z() {
        let mut c = Circuit::new(1);
        for _ in 0..4 { c.apply(Gate::t(0)); }
        let out = run(&c);
        assert_eq!(count_phase_gates(&out), 1);
        assert_eq!(count_t_gates(&out), 0);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn eight_t_cancel() {
        let mut c = Circuit::new(1);
        for _ in 0..8 { c.apply(Gate::t(0)); }
        let out = run(&c);
        assert_eq!(count_phase_gates(&out), 0);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn h_blocks_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = run(&c);
        assert_eq!(count_t_gates(&out), 2);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn cnot_enables_cross_qubit_merge() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        // After CNOT, q1 parity = q0 XOR q1. Not the same as q0.
        let out = run(&c);
        assert_eq!(count_t_gates(&out), 2);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn cnot_same_parity_merge() {
        // q0 = v0, q1 = v1
        // CNOT 0->1: q1 = v0 ^ v1
        // CNOT 0->1 again: q1 = v0 ^ v1 ^ v0 = v1
        // So T on q0, CNOT, CNOT, T on q1 should NOT merge (v0 != v1).
        // But: q0 = v0 and q1 after double CNOT = v1, so no merge.
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        let out = run(&c);
        // v0 and v1 are different, no merge
        assert_eq!(count_t_gates(&out), 2);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn x_flips_parity_prevents_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let out = run(&c);
        // X flips negation, so parities differ
        assert_eq!(count_t_gates(&out), 2);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn double_x_allows_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let out = run(&c);
        assert_eq!(count_t_gates(&out), 0); // merged to S
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn rz_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(0.5, 0));
        c.apply(Gate::rz(0.3, 0));
        let out = run(&c);
        assert_eq!(count_phase_gates(&out), 1);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn z_z_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::z(0));
        let out = run(&c);
        assert_eq!(count_phase_gates(&out), 0);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn ccx_refreshes_target() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(2));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::t(2));
        let out = run(&c);
        // CCX refreshes target parity, so no merge
        assert_eq!(count_t_gates(&out), 2);
    }

    #[test]
    fn phase_fold_preserves_semantics_3q() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::s(2));
        c.apply(Gate::t(0));
        let out = run(&c);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn agrees_with_hash_version_simple() {
        use crate::phase_fold_global::phase_fold_global;
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        c.apply(Gate::tdg(1));
        c.apply(Gate::s(0));
        let hash_out = phase_fold_global(&c, &ProgressBar::hidden());
        let expr_out = run(&c);
        assert_eq!(hash_out.gates.len(), expr_out.gates.len());
        assert!(circuits_equiv(&hash_out, &expr_out, TOL));
    }

    #[test]
    fn agrees_with_hash_version_3q() {
        use crate::phase_fold_global::phase_fold_global;
        let mut c = Circuit::new(3);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(0));
        c.apply(Gate::t(1));
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::s(2));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        c.apply(Gate::h(1));
        c.apply(Gate::t(1));
        let hash_out = phase_fold_global(&c, &ProgressBar::hidden());
        let expr_out = run(&c);
        // Hash version folds across X (complementary-parity lookup); expr version
        // does not. Both must be semantically equivalent to the input; the hash
        // output may be strictly smaller.
        assert!(hash_out.gates.len() <= expr_out.gates.len());
        assert!(circuits_equiv(&c, &hash_out, TOL));
        assert!(circuits_equiv(&c, &expr_out, TOL));
    }
}
