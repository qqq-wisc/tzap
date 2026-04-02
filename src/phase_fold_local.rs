//! Local phase folding pass: for each z-rotation gate, scan forward for another
//! z-rotation on the same qubit. If the subcircuit between them (on ≤5 qubits)
//! is equivalent to having the two rotations merged, replace them with a single
//! rotation. This catches cases where Rz commutes through intervening gates
//! (e.g. Rz on a CNOT control).

use std::f64::consts::PI;

use indicatif::ProgressBar;

use crate::circuit::{Circuit, Gate, Qubit, qubits_of, remap_subcircuit};
use crate::pass::Pass;
use crate::rule::{emit_z_rotation, z_angle};
use crate::unitary::circuits_equiv;

const MAX_QUBITS: usize = 5;

pub struct PhaseFoldLocal;

impl Pass for PhaseFoldLocal {
    fn name(&self) -> &str { "rotation-merging" }

    fn run_with_progress(&self, circuit: &Circuit, _pb: &ProgressBar) -> Circuit {
        let mut gates = circuit.gates.clone();
        loop {
            let prev_len = gates.len();
            if let Some((i, j, merged_angle)) = find_merge(&gates) {
                let q = z_angle(&gates[i]).unwrap().0;
                let mut new_gates = Vec::with_capacity(gates.len());
                for (k, g) in gates.iter().enumerate() {
                    if k == i {
                        // Drop the first rotation
                        continue;
                    }
                    if k == j {
                        // Replace the second rotation with the merged one (or drop if zero)
                        if let Some(mg) = nonzero_rotation(merged_angle, q) {
                            new_gates.push(mg);
                        }
                    } else {
                        new_gates.push(g.clone());
                    }
                }
                gates = new_gates;
                // Only continue if we actually reduced gate count
                if gates.len() >= prev_len {
                    break;
                }
            } else {
                break;
            }
        }
        let mut out = Circuit::new(circuit.num_qubits);
        for g in gates {
            match &g {
                Gate::rz(theta, q) => {
                    for named in emit_z_rotation(*theta, *q) {
                        out.apply(named);
                    }
                }
                _ => out.apply(g),
            }
        }
        out
    }
}

/// Emit a single z-rotation gate, or None if the angle is 0 mod 2π.
/// Uses raw Rz internally so merging passes don't re-merge decomposed pairs.
fn nonzero_rotation(angle: f64, q: Qubit) -> Option<Gate> {
    let norm = angle.rem_euclid(2.0 * PI);
    if norm.abs() < 1e-10 || (norm - 2.0 * PI).abs() < 1e-10 {
        None
    } else {
        Some(Gate::rz(angle, q))
    }
}

/// Scan for the first mergeable pair of z-rotations on the same qubit.
/// For each z-rotation at index i, extend a window forward tracking the qubit set.
/// When another z-rotation on the same qubit is found at index j, verify the merge
/// by comparing the unitary of gates[i..=j] against gates[i+1..j] + Rz(θ₁+θ₂).
fn find_merge(gates: &[Gate]) -> Option<(usize, usize, f64)> {
    for i in 0..gates.len() {
        let Some((q, angle_i)) = z_angle(&gates[i]) else { continue };
        let mut qubits = vec![q];

        for j in (i + 1)..gates.len() {
            for gq in qubits_of(&gates[j]) {
                if !qubits.contains(&gq) {
                    qubits.push(gq);
                }
            }
            if qubits.len() > MAX_QUBITS {
                break;
            }

            let Some((q_j, angle_j)) = z_angle(&gates[j]) else { continue };
            if q_j != q {
                continue;
            }

            let merged_angle = angle_i + angle_j;
            let original = remap_subcircuit(&gates[i..=j], &qubits);

            let mut repl_gates: Vec<Gate> = gates[i + 1..j].to_vec();
            if let Some(mg) = nonzero_rotation(merged_angle, q) {
                repl_gates.push(mg);
            }
            let replacement = remap_subcircuit(&repl_gates, &qubits);

            if circuits_equiv(&original, &replacement, 1e-10) {
                return Some((i, j, merged_angle));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_circuit(num_qubits: usize, gates: Vec<Gate>) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for g in gates { c.apply(g); }
        c
    }

    // --- adjacent merges ---

    #[test]
    fn t_t_merges_to_s() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::s(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn t_tdg_cancels() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::tdg(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_s_merges_to_z() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::s(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::z(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_sdg_cancels() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::sdg(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn z_z_cancels() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::z(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- merge through CNOT (z on control commutes with CNOT) ---

    #[test]
    fn t_cnot_t_merges_on_control() {
        let c = make_circuit(2, vec![
            Gate::t(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::t(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_cnot_sdg_cancels_on_control() {
        let c = make_circuit(2, vec![
            Gate::s(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::sdg(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn t_cnot_t_blocked_on_target() {
        let c = make_circuit(2, vec![
            Gate::t(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::t(1),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    // --- merge through chain of CNOTs ---

    #[test]
    fn t_two_cnots_t_merges_on_control() {
        let c = make_circuit(3, vec![
            Gate::t(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 2 },
            Gate::t(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- merge on different qubits ---

    #[test]
    fn merge_on_high_qubit() {
        let c = make_circuit(5, vec![Gate::s(4), Gate::s(4)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::z(4)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- merge past unrelated gates ---

    #[test]
    fn merge_past_gate_on_other_qubit() {
        let c = make_circuit(2, vec![Gate::t(0), Gate::h(1), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- cascading merges ---

    #[test]
    fn cascading_three_t() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::t(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2); // S + T (3π/4)
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn eight_t_cancels() {
        let c = make_circuit(1, vec![
            Gate::t(0), Gate::t(0), Gate::t(0), Gate::t(0),
            Gate::t(0), Gate::t(0), Gate::t(0), Gate::t(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- non-merge preserved ---

    #[test]
    fn h_not_merged() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::t(0), Gate::h(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    #[test]
    fn empty_circuit() {
        let c = Circuit::new(2);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 0);
    }

    #[test]
    fn single_gate_preserved() {
        let c = make_circuit(1, vec![Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
    }

    // --- qubit limit ---

    #[test]
    fn blocked_by_qubit_limit() {
        let c = make_circuit(6, vec![
            Gate::t(0),
            Gate::cnot { control: 1, target: 2 },
            Gate::cnot { control: 3, target: 4 },
            Gate::h(5),
            Gate::t(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 5);
    }

    #[test]
    fn within_qubit_limit() {
        let c = make_circuit(4, vec![
            Gate::t(0),
            Gate::cnot { control: 1, target: 2 },
            Gate::h(3),
            Gate::t(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn merge_through_cnot_sandwich() {
        let c = make_circuit(3, vec![
            Gate::t(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 2 },
            Gate::cnot { control: 0, target: 1 },
            Gate::tdg(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn surrounding_gates_preserved() {
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::t(1), Gate::t(1), Gate::h(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- all z-rotation gate combinations ---

    #[test]
    fn tdg_tdg_merges_to_sdg() {
        let c = make_circuit(1, vec![Gate::tdg(0), Gate::tdg(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::sdg(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn sdg_sdg_merges_to_z() {
        let c = make_circuit(1, vec![Gate::sdg(0), Gate::sdg(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::z(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn z_t_merges() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn z_tdg_merges() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::tdg(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2); // S + T (3π/4)
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_t_merges() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2); // S + T (3π/4)
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn sdg_t_merges_to_tdg() {
        let c = make_circuit(1, vec![Gate::sdg(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::tdg(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_tdg_merges_to_t() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::tdg(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rz_rz_merges() {
        let c = make_circuit(1, vec![Gate::rz(0.3, 0), Gate::rz(0.7, 0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rz_cancel() {
        let c = make_circuit(1, vec![Gate::rz(1.5, 0), Gate::rz(-1.5, 0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- blocking scenarios ---

    #[test]
    fn merge_through_multiple_cnots_on_control() {
        let c = make_circuit(4, vec![
            Gate::s(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 2 },
            Gate::cnot { control: 0, target: 3 },
            Gate::sdg(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn blocked_by_h_on_same_qubit() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::h(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    #[test]
    fn blocked_by_x_on_same_qubit() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::x(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    #[test]
    fn blocked_by_cnot_on_target() {
        let c = make_circuit(2, vec![
            Gate::s(0),
            Gate::cnot { control: 1, target: 0 },
            Gate::sdg(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    // --- multiple independent merges ---

    #[test]
    fn two_independent_merges() {
        let c = make_circuit(2, vec![
            Gate::t(0), Gate::t(0), Gate::t(1), Gate::t(1),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn interleaved_merges() {
        let c = make_circuit(2, vec![
            Gate::t(0), Gate::t(1), Gate::t(0), Gate::t(1),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn first_blocked_second_merges() {
        // First T can't merge past H, but the two T's after H merge.
        let c = make_circuit(1, vec![
            Gate::t(0), Gate::h(0), Gate::t(0), Gate::t(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- real circuit patterns ---

    #[test]
    fn toffoli_decomposition_has_mergeable_rotations() {
        let c = make_circuit(3, vec![
            Gate::h(2),
            Gate::cnot { control: 1, target: 2 },
            Gate::tdg(2),
            Gate::cnot { control: 0, target: 2 },
            Gate::t(2),
            Gate::cnot { control: 1, target: 2 },
            Gate::tdg(2),
            Gate::cnot { control: 0, target: 2 },
            Gate::t(1), Gate::t(2), Gate::h(2),
            Gate::cnot { control: 0, target: 1 },
            Gate::t(0), Gate::tdg(1),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = PhaseFoldLocal.run(&c);
        let orig_t = c.gates.iter().filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count();
        let new_t = r.gates.iter().filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count();
        assert!(new_t <= orig_t);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn mod5_4_phase_fold_local_reduces() {
        use crate::decompose::DecomposeToffoli;
        use crate::cancel::CancelPairs;
        use crate::phase_fold_global::PhaseFoldGlobal;
        use crate::pass::Pass;

        let qasm = std::fs::read_to_string("qasm/mod5_4.qasm").unwrap();
        let circuit = Circuit::from_qasm(&qasm).unwrap();

        let after_pf = PhaseFoldGlobal.run(&CancelPairs.run(&DecomposeToffoli.run(&circuit)));
        let t_before = after_pf.gates.iter()
            .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count();

        let after_rm = PhaseFoldLocal.run(&after_pf);
        let t_after = after_rm.gates.iter()
            .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count();

        assert!(circuits_equiv(&after_pf, &after_rm, 1e-10));
        assert!(t_after < t_before, "{t_before} -> {t_after}");
        assert_eq!(t_after, 8);
    }

    #[test]
    fn mod5_4_full_pipeline_equiv() {
        use crate::decompose::DecomposeToffoli;
        use crate::cancel::CancelPairs;
        use crate::phase_fold_global::PhaseFoldGlobal;
        use crate::pass::Pass;

        let qasm = std::fs::read_to_string("qasm/mod5_4.qasm").unwrap();
        let original = Circuit::from_qasm(&qasm).unwrap();

        let r = PhaseFoldLocal.run(
            &PhaseFoldGlobal.run(
                &CancelPairs.run(
                    &DecomposeToffoli.run(&original))));

        assert!(circuits_equiv(&original, &r, 1e-10));
    }

    // --- edge cases ---

    #[test]
    fn skip_non_phase_gate_between() {
        // CNOT q1,q0 has q0 as target — blocks the merge
        let c = make_circuit(2, vec![
            Gate::t(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 1, target: 0 },
            Gate::t(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 4);
    }

    #[test]
    fn merge_skips_intermediate_z_rotation_on_different_qubit() {
        let c = make_circuit(2, vec![Gate::t(0), Gate::s(1), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn four_qubit_merge() {
        let c = make_circuit(4, vec![
            Gate::s(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 2 },
            Gate::cnot { control: 0, target: 3 },
            Gate::s(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn five_qubit_merge() {
        let c = make_circuit(5, vec![
            Gate::t(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 2 },
            Gate::cnot { control: 0, target: 3 },
            Gate::cnot { control: 0, target: 4 },
            Gate::tdg(0),
        ]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn merge_through_z_on_same_qubit() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::t(0)]);
        let r = PhaseFoldLocal.run(&c);
        assert_eq!(r.gates.len(), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &r, 1e-10));
    }
}
