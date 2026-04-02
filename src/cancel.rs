use indicatif::ProgressBar;

use crate::circuit::{Circuit, Gate, Qubit};
use crate::pass::Pass;

/// Cancel adjacent self-inverse gate pairs (HH, XX, CNOT-CNOT) in O(n),
/// allowing commutation past gates on non-overlapping qubits.
/// Handles cascading: cancelling a pair may expose new adjacent pairs.
///
/// Uses per-qubit stacks to find the blocking gate in O(1) instead of
/// scanning backward through the entire result list.
/// Gates are tracked by index into the original slice — only surviving gates
/// are cloned at the end.
fn cancel_pairs(gates: &[Gate], num_qubits: usize, pb: &ProgressBar) -> Vec<Gate> {
    let len = gates.len();
    let mut skip = vec![false; len];
    // Per-qubit stack of gate indices — tracks the most recent gate on each qubit.
    let mut qubit_stacks: Vec<Vec<usize>> = vec![Vec::new(); num_qubits];

    for (i, gate) in gates.iter().enumerate() {
        pb.inc(1);

        if is_self_inverse(gate) {
            let (n, qs) = qubits_of(gate);
            // The blocker is the latest gate touching any of this gate's qubits.
            let mut blocker: Option<usize> = None;
            for j in 0..n {
                if let Some(&last) = qubit_stacks[qs[j]].last() {
                    blocker = Some(match blocker {
                        Some(b) => b.max(last),
                        None => last,
                    });
                }
            }
            if let Some(block_idx) = blocker {
                if gates_equal(&gates[block_idx], gate) {
                    // Cancel both gates; pop the blocker from all relevant qubit stacks.
                    skip[block_idx] = true;
                    skip[i] = true;
                    for j in 0..n {
                        debug_assert_eq!(*qubit_stacks[qs[j]].last().unwrap(), block_idx);
                        qubit_stacks[qs[j]].pop();
                    }
                    continue;
                }
            }
        }

        let (n, qs) = qubits_of(gate);
        for j in 0..n {
            qubit_stacks[qs[j]].push(i);
        }
    }

    gates.iter().enumerate()
        .filter(|(i, _)| !skip[*i])
        .map(|(_, g)| g.clone())
        .collect()
}

fn is_self_inverse(gate: &Gate) -> bool {
    matches!(gate, Gate::h(_) | Gate::x(_) | Gate::z(_)
        | Gate::cnot { .. } | Gate::ccx { .. })
}

fn gates_equal(a: &Gate, b: &Gate) -> bool {
    match (a, b) {
        (Gate::h(a), Gate::h(b))
        | (Gate::x(a), Gate::x(b))
        | (Gate::z(a), Gate::z(b)) => a == b,
        (Gate::cnot { control: ac, target: at },
         Gate::cnot { control: bc, target: bt }) => ac == bc && at == bt,
        (Gate::ccx { control1: a1, control2: a2, target: at },
         Gate::ccx { control1: b1, control2: b2, target: bt }) => a1 == b1 && a2 == b2 && at == bt,
        _ => false,
    }
}

fn qubits_of(gate: &Gate) -> (usize, [Qubit; 3]) {
    match gate {
        Gate::x(q) | Gate::h(q) | Gate::s(q) | Gate::sdg(q) | Gate::z(q)
        | Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => (1, [*q, 0, 0]),
        Gate::cnot { control, target } => (2, [*control, *target, 0]),
        Gate::ccx { control1, control2, target } => (3, [*control1, *control2, *target]),
    }
}

pub struct CancelPairs;

impl Pass for CancelPairs {
    fn name(&self) -> &str { "Pair cancellation" }
    fn run_with_progress(&self, circuit: &Circuit, pb: &ProgressBar) -> Circuit {
        let gates = cancel_pairs(&circuit.gates, circuit.num_qubits, pb);
        let has_toffoli = gates.iter().any(|g| matches!(g, Gate::ccx { .. }));
        Circuit { num_qubits: circuit.num_qubits, gates, has_toffoli }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rule::{Rule, RulePass};
    use crate::unitary::circuits_equiv;

    fn make_circuit(num_qubits: usize, gates: Vec<Gate>) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for g in gates {
            c.apply(g);
        }
        c
    }

    /// Rewrite rules only (no cancellation rules — those are handled by cancel_pairs).
    fn rewrite_rules() -> Vec<Rule> {
        vec![
            // #1: H q; S q; H q → Sdg q; H q; Sdg q
            Rule::new(
                make_circuit(1, vec![Gate::h(0), Gate::s(0), Gate::h(0)]),
                make_circuit(1, vec![Gate::sdg(0), Gate::h(0), Gate::sdg(0)]),
            ),
            // #2: H q; Sdg q; H q → S q; H q; S q
            Rule::new(
                make_circuit(1, vec![Gate::h(0), Gate::sdg(0), Gate::h(0)]),
                make_circuit(1, vec![Gate::s(0), Gate::h(0), Gate::s(0)]),
            ),
            // #3: H q1; H q2; CNOT q1 q2; H q1; H q2 → CNOT q2 q1
            Rule::new(
                make_circuit(2, vec![
                    Gate::h(0), Gate::h(1),
                    Gate::cnot { control: 0, target: 1 },
                    Gate::h(0), Gate::h(1),
                ]),
                make_circuit(2, vec![
                    Gate::cnot { control: 1, target: 0 },
                ]),
            ),
            // #4: H q2; S q2; CNOT q1 q2; Sdg q2; H q2 → Sdg q2; CNOT q1 q2; S q2
            Rule::new(
                make_circuit(2, vec![
                    Gate::h(1), Gate::s(1),
                    Gate::cnot { control: 0, target: 1 },
                    Gate::sdg(1), Gate::h(1),
                ]),
                make_circuit(2, vec![
                    Gate::sdg(1),
                    Gate::cnot { control: 0, target: 1 },
                    Gate::s(1),
                ]),
            ),
            // #5: H q2; Sdg q2; CNOT q1 q2; S q2; H q2 → S q2; CNOT q1 q2; Sdg q2
            Rule::new(
                make_circuit(2, vec![
                    Gate::h(1), Gate::sdg(1),
                    Gate::cnot { control: 0, target: 1 },
                    Gate::s(1), Gate::h(1),
                ]),
                make_circuit(2, vec![
                    Gate::s(1),
                    Gate::cnot { control: 0, target: 1 },
                    Gate::sdg(1),
                ]),
            ),
        ]
    }

    /// Run cancel_pairs + rewrite rules (for testing rules that CancelPairs no longer runs).
    fn run_with_rules(circuit: &Circuit) -> Circuit {
        let rules = rewrite_rules();
        let pass = RulePass { rules: &rules };
        let pb = ProgressBar::hidden();
        let mut gates = cancel_pairs(&circuit.gates, circuit.num_qubits, &pb);
        loop {
            let mut c = Circuit::new(circuit.num_qubits);
            for g in &gates { c.apply(g.clone()); }
            let rewritten = pass.run(&c);
            let cancelled = cancel_pairs(&rewritten.gates, circuit.num_qubits, &pb);
            if cancelled.len() >= gates.len() {
                gates = cancelled;
                break;
            }
            gates = cancelled;
        }
        let mut out = Circuit::new(circuit.num_qubits);
        for g in gates { out.apply(g); }
        out
    }

    #[test]
    fn hh_cancel() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::h(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel() {
        let c = make_circuit(1, vec![Gate::x(0), Gate::x(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel() {
        let c = make_circuit(2, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_s_h_reduction() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::s(0), Gate::h(0)]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::sdg(0)));
        assert!(matches!(&r.gates[1], Gate::h(0)));
        assert!(matches!(&r.gates[2], Gate::sdg(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_sdg_h_reduction() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::sdg(0), Gate::h(0)]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::s(0)));
        assert!(matches!(&r.gates[1], Gate::h(0)));
        assert!(matches!(&r.gates[2], Gate::s(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_flip() {
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::h(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::cnot { control: 1, target: 0 }));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_s_cnot_sdg_h_reduction() {
        let c = make_circuit(2, vec![
            Gate::h(1), Gate::s(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::sdg(1), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::sdg(1)));
        assert!(matches!(&r.gates[1], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&r.gates[2], Gate::s(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_sdg_cnot_s_h_reduction() {
        let c = make_circuit(2, vec![
            Gate::h(1), Gate::sdg(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::s(1), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::s(1)));
        assert!(matches!(&r.gates[1], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&r.gates[2], Gate::sdg(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn reduces_h_count() {
        // HH cancel followed by CNOT flip pattern
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::h(0), // HH cancel
            Gate::h(0), Gate::h(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        // HH cancel removes 2, CNOT flip removes 4 H's + replaces CNOT
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- HH cancel: varied settings ---

    #[test]
    fn hh_cancel_different_qubit() {
        let c = make_circuit(4, vec![Gate::h(3), Gate::h(3)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_skips_unrelated_gate() {
        // H q0; T q1; H q0 — T on different qubit doesn't block
        let c = make_circuit(2, vec![Gate::h(0), Gate::t(1), Gate::h(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_blocked_by_same_qubit() {
        // H q0; T q0; H q0 — T on same qubit blocks HH cancel
        let c = make_circuit(1, vec![Gate::h(0), Gate::t(0), Gate::h(0)]);
        let r = CancelPairs.run(&c);
        // Should not cancel as HH; instead rule #1/#2 may or may not apply
        // T is not S or Sdg so no hadamard reduction either — unchanged
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_multiple_pairs() {
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_parallel_qubits() {
        // H on q0 and q1 independently
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::h(1), Gate::h(0), Gate::h(1),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- XX cancel: varied settings ---

    #[test]
    fn xx_cancel_different_qubit() {
        let c = make_circuit(3, vec![Gate::x(2), Gate::x(2)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel_skips_unrelated_gate() {
        let c = make_circuit(2, vec![Gate::x(0), Gate::h(1), Gate::x(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::h(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel_blocked_by_same_qubit() {
        let c = make_circuit(1, vec![Gate::x(0), Gate::z(0), Gate::x(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel_multiple_pairs() {
        let c = make_circuit(1, vec![
            Gate::x(0), Gate::x(0), Gate::x(0), Gate::x(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- CNOT cancel: varied settings ---

    #[test]
    fn cnot_cancel_different_qubits() {
        let c = make_circuit(5, vec![
            Gate::cnot { control: 3, target: 4 },
            Gate::cnot { control: 3, target: 4 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_skips_unrelated_gate() {
        // CNOT q0,q1; T q2; CNOT q0,q1 — T on q2 doesn't interfere
        let c = make_circuit(3, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::t(2),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(2)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_blocked_by_gate_on_control() {
        let c = make_circuit(2, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::t(0),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_blocked_by_gate_on_target() {
        let c = make_circuit(2, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::t(1),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_no_match_different_direction() {
        // CNOT q0,q1 then CNOT q1,q0 — different direction, should NOT cancel
        let c = make_circuit(2, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 1, target: 0 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_multiple_pairs() {
        let c = make_circuit(2, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Rule #1: H S H → Sdg H Sdg ---

    #[test]
    fn h_s_h_different_qubit() {
        let c = make_circuit(3, vec![Gate::h(2), Gate::s(2), Gate::h(2)]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::sdg(2)));
        assert!(matches!(&r.gates[1], Gate::h(2)));
        assert!(matches!(&r.gates[2], Gate::sdg(2)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_s_h_skips_unrelated_gate() {
        // H q0; T q1; S q0; H q0 — T on q1 between H and S doesn't block
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::t(1), Gate::s(0), Gate::h(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 4); // T q1 + Sdg H Sdg
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_s_h_blocked_by_same_qubit() {
        // H q0; X q0; S q0; H q0 — X on same qubit blocks
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::x(0), Gate::s(0), Gate::h(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_s_h_with_surrounding_gates() {
        // T q0; H q0; S q0; H q0; T q0
        let c = make_circuit(1, vec![
            Gate::t(0), Gate::h(0), Gate::s(0), Gate::h(0), Gate::t(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 5);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::sdg(0)));
        assert!(matches!(&r.gates[2], Gate::h(0)));
        assert!(matches!(&r.gates[3], Gate::sdg(0)));
        assert!(matches!(&r.gates[4], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Rule #2: H Sdg H → S H S ---

    #[test]
    fn h_sdg_h_different_qubit() {
        let c = make_circuit(4, vec![Gate::h(3), Gate::sdg(3), Gate::h(3)]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::s(3)));
        assert!(matches!(&r.gates[1], Gate::h(3)));
        assert!(matches!(&r.gates[2], Gate::s(3)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_sdg_h_skips_unrelated_gate() {
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::x(1), Gate::sdg(0), Gate::h(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 4); // X q1 + S H S
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_sdg_h_blocked_by_same_qubit() {
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::x(0), Gate::sdg(0), Gate::h(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_sdg_h_with_surrounding_gates() {
        let c = make_circuit(1, vec![
            Gate::z(0), Gate::h(0), Gate::sdg(0), Gate::h(0), Gate::z(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 5);
        assert!(matches!(&r.gates[0], Gate::z(0)));
        assert!(matches!(&r.gates[1], Gate::s(0)));
        assert!(matches!(&r.gates[2], Gate::h(0)));
        assert!(matches!(&r.gates[3], Gate::s(0)));
        assert!(matches!(&r.gates[4], Gate::z(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Rule #3: H H CNOT H H → CNOT flipped ---

    #[test]
    fn cnot_flip_different_qubits() {
        let c = make_circuit(5, vec![
            Gate::h(2), Gate::h(4),
            Gate::cnot { control: 2, target: 4 },
            Gate::h(2), Gate::h(4),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::cnot { control: 4, target: 2 }));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_flip_skips_unrelated_gate() {
        // Gate on q2 between the H's and CNOT on q0,q1
        let c = make_circuit(3, vec![
            Gate::h(0), Gate::h(1),
            Gate::t(2),
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 2); // T q2 + CNOT flipped
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_flip_blocked_by_gate_on_control() {
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::h(1),
            Gate::t(0), // blocks — touches q0 which is mapped
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        // Rule #3 can't match; some sub-patterns may apply but not the full flip
        assert!(r.gates.len() > 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_flip_multiple_separated() {
        // Two independent CNOT flips on different qubit pairs
        let c = make_circuit(4, vec![
            Gate::h(0), Gate::h(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0), Gate::h(1),
            Gate::h(2), Gate::h(3),
            Gate::cnot { control: 2, target: 3 },
            Gate::h(2), Gate::h(3),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::cnot { control: 1, target: 0 }));
        assert!(matches!(&r.gates[1], Gate::cnot { control: 3, target: 2 }));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Rule #4: H S CNOT Sdg H → Sdg CNOT S ---

    #[test]
    fn rule4_different_qubits() {
        let c = make_circuit(4, vec![
            Gate::h(3), Gate::s(3),
            Gate::cnot { control: 1, target: 3 },
            Gate::sdg(3), Gate::h(3),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::sdg(3)));
        assert!(matches!(&r.gates[1], Gate::cnot { control: 1, target: 3 }));
        assert!(matches!(&r.gates[2], Gate::s(3)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rule4_skips_unrelated_gate() {
        // Gate on q2 shouldn't block rule matching on q0 (control) and q1 (target)
        let c = make_circuit(3, vec![
            Gate::h(1), Gate::s(1),
            Gate::t(2),
            Gate::cnot { control: 0, target: 1 },
            Gate::sdg(1), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 4); // T q2 + Sdg CNOT S
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rule4_blocked_by_gate_on_target() {
        let c = make_circuit(2, vec![
            Gate::h(1), Gate::s(1),
            Gate::x(1), // blocks — touches target qubit
            Gate::cnot { control: 0, target: 1 },
            Gate::sdg(1), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert!(r.gates.len() >= 5); // no reduction of the full pattern
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rule4_with_surrounding_gates() {
        let c = make_circuit(2, vec![
            Gate::t(0),
            Gate::h(1), Gate::s(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::sdg(1), Gate::h(1),
            Gate::t(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 5); // T + Sdg CNOT S + T
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::sdg(1)));
        assert!(matches!(&r.gates[2], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&r.gates[3], Gate::s(1)));
        assert!(matches!(&r.gates[4], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Rule #5: H Sdg CNOT S H → S CNOT Sdg ---

    #[test]
    fn rule5_different_qubits() {
        let c = make_circuit(4, vec![
            Gate::h(3), Gate::sdg(3),
            Gate::cnot { control: 2, target: 3 },
            Gate::s(3), Gate::h(3),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::s(3)));
        assert!(matches!(&r.gates[1], Gate::cnot { control: 2, target: 3 }));
        assert!(matches!(&r.gates[2], Gate::sdg(3)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rule5_skips_unrelated_gate() {
        let c = make_circuit(3, vec![
            Gate::h(1), Gate::sdg(1),
            Gate::t(2),
            Gate::cnot { control: 0, target: 1 },
            Gate::s(1), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 4); // T q2 + S CNOT Sdg
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rule5_blocked_by_gate_on_target() {
        let c = make_circuit(2, vec![
            Gate::h(1), Gate::sdg(1),
            Gate::x(1), // blocks
            Gate::cnot { control: 0, target: 1 },
            Gate::s(1), Gate::h(1),
        ]);
        let r = run_with_rules(&c);
        assert!(r.gates.len() >= 5);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rule5_with_surrounding_gates() {
        let c = make_circuit(2, vec![
            Gate::z(0),
            Gate::h(1), Gate::sdg(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::s(1), Gate::h(1),
            Gate::z(0),
        ]);
        let r = run_with_rules(&c);
        assert_eq!(r.gates.len(), 5);
        assert!(matches!(&r.gates[0], Gate::z(0)));
        assert!(matches!(&r.gates[1], Gate::s(1)));
        assert!(matches!(&r.gates[2], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&r.gates[3], Gate::sdg(1)));
        assert!(matches!(&r.gates[4], Gate::z(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- No-match / preservation tests ---

    #[test]
    fn no_match_preserves_circuit() {
        let c = make_circuit(2, vec![
            Gate::t(0), Gate::s(1), Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn single_h_preserved() {
        let c = make_circuit(1, vec![Gate::h(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::h(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn empty_circuit() {
        let c = Circuit::new(2);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
    }

    // --- ZZ cancel ---

    #[test]
    fn zz_cancel() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::z(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn zz_cancel_skips_unrelated_gate() {
        let c = make_circuit(2, vec![Gate::z(0), Gate::t(1), Gate::z(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn zz_cancel_blocked_by_same_qubit() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::h(0), Gate::z(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- CCX (Toffoli) cancel ---

    #[test]
    fn ccx_cancel() {
        let c = make_circuit(3, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_skips_unrelated_gate() {
        let c = make_circuit(4, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::t(3),
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(3)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_blocked_by_control1() {
        let c = make_circuit(3, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::h(0),
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_blocked_by_control2() {
        let c = make_circuit(3, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::h(1),
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_blocked_by_target() {
        let c = make_circuit(3, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::h(2),
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_no_cancel_different_controls() {
        let c = make_circuit(4, vec![
            Gate::ccx { control1: 0, control2: 1, target: 3 },
            Gate::ccx { control1: 0, control2: 2, target: 3 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_no_cancel_different_target() {
        let c = make_circuit(4, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::ccx { control1: 0, control2: 1, target: 3 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_multiple_pairs() {
        let c = make_circuit(3, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Cascading cancellation ---

    #[test]
    fn cascade_nested_h() {
        // H H H H — inner pair cancels, exposing outer pair
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_six_h() {
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_odd_count_leaves_one() {
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_nested_cnot() {
        // CNOT CNOT CNOT CNOT — fully cancels in one pass
        let c = make_circuit(2, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_mixed_self_inverse() {
        // H X H X — H cancels with H, X cancels with X? No — X blocks H.
        // H(q0) X(q0) H(q0) X(q0): X blocks H cancel, X blocks X cancel
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::x(0), Gate::h(0), Gate::x(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 4); // nothing cancels
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_different_qubits_interleaved() {
        // H(q0) H(q1) H(q1) H(q0) — inner H(q1) pair cancels, then outer H(q0) cancels
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::h(1), Gate::h(1), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_deep_nesting() {
        // H(0) X(0) Z(0) Z(0) X(0) H(0) — Z pair cancels, exposes X pair, exposes H pair
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::x(0), Gate::z(0), Gate::z(0), Gate::x(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_deep_nesting_with_residual() {
        // T(0) H(0) X(0) X(0) H(0) T(0) — X cancels, H cancels, T is not self-inverse
        let c = make_circuit(1, vec![
            Gate::t(0), Gate::h(0), Gate::x(0), Gate::x(0), Gate::h(0), Gate::t(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_cnot_nested_in_h() {
        // H(0) CNOT(0,1) CNOT(0,1) H(0) — CNOT pair cancels, H pair cancels
        let c = make_circuit(2, vec![
            Gate::h(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Commutation through unrelated qubits ---

    #[test]
    fn commute_h_past_many_unrelated() {
        // H(0); T(1); S(2); Tdg(3); H(0) — all middle gates on different qubits
        let c = make_circuit(4, vec![
            Gate::h(0), Gate::t(1), Gate::s(2), Gate::tdg(3), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn commute_cnot_past_unrelated_qubits() {
        // CNOT(0,1); H(2); T(3); CNOT(0,1)
        let c = make_circuit(4, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::h(2), Gate::t(3),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn commute_ccx_past_unrelated_qubits() {
        // CCX(0,1,2); H(3); T(4); CCX(0,1,2)
        let c = make_circuit(5, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::h(3), Gate::t(4),
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_blocked_by_gate_on_either_qubit() {
        // Gate on control blocks
        let c1 = make_circuit(3, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r1 = CancelPairs.run(&c1);
        assert_eq!(r1.gates.len(), 3);
        assert!(circuits_equiv(&c1, &r1, 1e-10));

        // Gate on target blocks
        let c2 = make_circuit(3, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::h(1),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r2 = CancelPairs.run(&c2);
        assert_eq!(r2.gates.len(), 3);
        assert!(circuits_equiv(&c2, &r2, 1e-10));

        // Gate on unrelated qubit doesn't block
        let c3 = make_circuit(3, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::h(2),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r3 = CancelPairs.run(&c3);
        assert_eq!(r3.gates.len(), 1);
        assert!(circuits_equiv(&c3, &r3, 1e-10));
    }

    // --- Non-self-inverse gates don't cancel ---

    #[test]
    fn t_t_no_cancel() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::t(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn tdg_tdg_no_cancel() {
        let c = make_circuit(1, vec![Gate::tdg(0), Gate::tdg(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_s_no_cancel() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::s(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn sdg_sdg_no_cancel() {
        let c = make_circuit(1, vec![Gate::sdg(0), Gate::sdg(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn t_tdg_no_cancel() {
        // T and Tdg are inverses of each other but cancel_pairs only handles self-inverse
        let c = make_circuit(1, vec![Gate::t(0), Gate::tdg(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_sdg_no_cancel() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::sdg(0)]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rz_rz_no_cancel() {
        let c = make_circuit(1, vec![
            Gate::rz(std::f64::consts::PI / 4.0, 0),
            Gate::rz(std::f64::consts::PI / 4.0, 0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Multi-qubit independence ---

    #[test]
    fn independent_cancellations_on_many_qubits() {
        // Each qubit has its own H-H pair, all should cancel independently
        let c = make_circuit(5, vec![
            Gate::h(0), Gate::h(1), Gate::h(2), Gate::h(3), Gate::h(4),
            Gate::h(4), Gate::h(3), Gate::h(2), Gate::h(1), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn partial_cancellation_mixed_qubits() {
        // q0: H H (cancels), q1: H T H (blocked), q2: X X (cancels)
        let c = make_circuit(3, vec![
            Gate::h(0), Gate::h(1), Gate::x(2),
            Gate::t(1),
            Gate::h(0), Gate::h(1), Gate::x(2),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3); // H(1), T(1), H(1) remain
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn interleaved_cancel_different_gate_types() {
        // H(0) X(1) H(0) X(1) — both pairs cancel through each other
        let c = make_circuit(2, vec![
            Gate::h(0), Gate::x(1), Gate::h(0), Gate::x(1),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- CNOT blocking edge cases ---

    #[test]
    fn cnot_blocks_h_on_shared_qubit() {
        // H(0); CNOT(0,1); H(0) — CNOT touches q0 so it blocks
        let c = make_circuit(2, vec![
            Gate::h(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_not_blocked_by_cnot_on_other_qubits() {
        // H(0); CNOT(1,2); H(0) — CNOT doesn't touch q0
        let c = make_circuit(3, vec![
            Gate::h(0),
            Gate::cnot { control: 1, target: 2 },
            Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_blocks_h_on_any_of_three_qubits() {
        for blocked_q in 0..3 {
            let c = make_circuit(4, vec![
                Gate::h(blocked_q),
                Gate::ccx { control1: 0, control2: 1, target: 2 },
                Gate::h(blocked_q),
            ]);
            let r = CancelPairs.run(&c);
            assert_eq!(r.gates.len(), 3, "H({}) should be blocked by CCX", blocked_q);
            assert!(circuits_equiv(&c, &r, 1e-10));
        }
        // q3 is unrelated — should cancel
        let c = make_circuit(4, vec![
            Gate::h(3),
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::h(3),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Order preservation ---

    #[test]
    fn surviving_gates_preserve_order() {
        let c = make_circuit(3, vec![
            Gate::t(0),
            Gate::h(1), Gate::h(1), // cancels
            Gate::s(0),
            Gate::x(2), Gate::x(2), // cancels
            Gate::tdg(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::s(0)));
        assert!(matches!(&r.gates[2], Gate::tdg(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn order_after_cascade() {
        // T(0) H(0) H(0) S(0) — H pair cancels, T and S remain in order
        let c = make_circuit(1, vec![
            Gate::t(0), Gate::h(0), Gate::h(0), Gate::s(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::s(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Stress / larger circuits ---

    #[test]
    fn many_pairs_single_qubit() {
        // 100 H-H pairs on q0
        let gates: Vec<Gate> = (0..200).map(|_| Gate::h(0)).collect();
        let c = make_circuit(1, gates);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn many_pairs_odd_leaves_one() {
        let gates: Vec<Gate> = (0..201).map(|_| Gate::h(0)).collect();
        let c = make_circuit(1, gates);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn alternating_cancel_no_cancel() {
        // H(0) H(0) T(0) H(0) H(0) T(0) — two H pairs cancel around T's
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::h(0), Gate::t(0),
            Gate::h(0), Gate::h(0), Gate::t(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn deeply_nested_cascade_8_layers() {
        // H X Z H H Z X H — 4 nested pairs, all cancel
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::x(0), Gate::z(0), Gate::h(0),
            Gate::h(0), Gate::z(0), Gate::x(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- has_toffoli flag ---

    #[test]
    fn has_toffoli_set_when_ccx_survives() {
        let c = make_circuit(4, vec![
            Gate::h(3), Gate::h(3), // cancels
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(r.has_toffoli);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn has_toffoli_cleared_when_ccx_cancelled() {
        let c = make_circuit(3, vec![
            Gate::ccx { control1: 0, control2: 1, target: 2 },
            Gate::ccx { control1: 0, control2: 1, target: 2 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(!r.has_toffoli);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Mixed multi-qubit gate interactions ---

    #[test]
    fn cnot_and_ccx_block_each_other_on_shared_qubit() {
        // CNOT(0,1); CCX(1,2,3); CNOT(0,1) — CCX touches q1 which blocks CNOT
        let c = make_circuit(4, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::ccx { control1: 1, control2: 2, target: 3 },
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_and_ccx_no_shared_qubit() {
        // CNOT(0,1); CCX(2,3,4); CNOT(0,1) — no shared qubits, CNOT cancels
        let c = make_circuit(5, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::ccx { control1: 2, control2: 3, target: 4 },
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn multiple_cnot_pairs_different_qubit_pairs() {
        // CNOT(0,1) CNOT(2,3) CNOT(2,3) CNOT(0,1)
        let c = make_circuit(4, vec![
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 2, target: 3 },
            Gate::cnot { control: 2, target: 3 },
            Gate::cnot { control: 0, target: 1 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Edge: single gate of each self-inverse type (no cancel) ---

    #[test]
    fn single_gate_each_type_preserved() {
        let c = make_circuit(5, vec![
            Gate::h(0),
            Gate::x(1),
            Gate::z(2),
            Gate::cnot { control: 3, target: 4 },
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Blocker is latest gate, not first ---

    #[test]
    fn blocker_is_latest_not_earliest() {
        // H(0); T(0); X(0); H(0) — blocker for second H is X (latest on q0), not T
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::t(0), Gate::x(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        // X is the blocker, not H — X != H so no cancel
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn blocker_matches_when_latest_is_same() {
        // H(0); T(0); H(0); H(0) — second H is blocked by first H? No, T is between.
        // Actually: stack is [H@0, T@1, H@2]. Third H@3 checks blocker = H@2 (latest on q0).
        // H@2 == H@3, so they cancel. Leaves [H@0, T@1].
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::t(0), Gate::h(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::h(0)));
        assert!(matches!(&r.gates[1], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Cascade after partial cancel ---

    #[test]
    fn cascade_partial_then_full() {
        // T(1) H(0) X(1) H(0) X(1) T(1)
        // H(0) pair cancels through X(1) (different qubit).
        // Then X(1) X(1) are adjacent and cancel.
        // Leaves T(1) T(1).
        let c = make_circuit(2, vec![
            Gate::t(1), Gate::h(0), Gate::x(1), Gate::h(0), Gate::x(1), Gate::t(1),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::t(1)));
        assert!(matches!(&r.gates[1], Gate::t(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Large cascade depth ---

    #[test]
    fn cascade_depth_10() {
        // 10 layers of nesting: H X Z H X Z X H Z X H Z X H ... all cancel
        // Simpler: alternating gate types nested symmetrically
        // Use q0: H(X(Z(Z(X(H()))))) = H X Z Z X H
        let c = make_circuit(1, vec![
            Gate::h(0), Gate::x(0), Gate::z(0),
            Gate::z(0), Gate::x(0), Gate::h(0),
        ]);
        let r = CancelPairs.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));

        // 5 layers deep
        let c2 = make_circuit(1, vec![
            Gate::h(0), Gate::x(0), Gate::z(0), Gate::x(0), Gate::h(0),
            Gate::h(0), Gate::x(0), Gate::z(0), Gate::x(0), Gate::h(0),
        ]);
        let r2 = CancelPairs.run(&c2);
        assert_eq!(r2.gates.len(), 0);
        assert!(circuits_equiv(&c2, &r2, 1e-10));
    }
}
