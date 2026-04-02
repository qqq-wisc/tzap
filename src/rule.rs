use crate::circuit::{Gate, Qubit};

/// If the gate is a z-axis rotation, return (qubit, angle).
pub fn z_angle(gate: &Gate) -> Option<(Qubit, f64)> {
    use std::f64::consts::PI;
    match gate {
        Gate::t(q) => Some((*q, PI / 4.0)),
        Gate::tdg(q) => Some((*q, -PI / 4.0)),
        Gate::s(q) => Some((*q, PI / 2.0)),
        Gate::sdg(q) => Some((*q, -PI / 2.0)),
        Gate::z(q) => Some((*q, PI)),
        Gate::rz(theta, q) => Some((*q, *theta)),
        _ => None,
    }
}

/// Emit a z-rotation as named gates when possible.
pub fn emit_z_rotation(angle: f64, q: Qubit) -> Vec<Gate> {
    use std::f64::consts::PI;
    let quarter = PI / 4.0;
    let n = angle.rem_euclid(2.0 * PI);
    let k = (n / quarter).round();
    if (n - k * quarter).abs() < 1e-10 {
        match k as u32 % 8 {
            0 => return vec![],
            1 => return vec![Gate::t(q)],
            2 => return vec![Gate::s(q)],
            3 => return vec![Gate::s(q), Gate::t(q)],
            4 => return vec![Gate::z(q)],
            5 => return vec![Gate::z(q), Gate::t(q)],
            6 => return vec![Gate::sdg(q)],
            7 => return vec![Gate::tdg(q)],
            _ => unreachable!(),
        }
    }
    vec![Gate::rz(angle, q)]
}

#[cfg(test)]
use crate::circuit::Circuit;
#[cfg(test)]
use crate::pass::Pass;

#[cfg(test)]
const MAX_RULE_QUBITS: usize = 3;
#[cfg(test)]
const MAX_PATTERN_GATES: usize = 8;

#[cfg(test)]
pub struct Rule {
    pub from: Circuit,
    pub to: Circuit,
}

/// Sentinel angle for wildcard Rz matching in rules.
/// Use `Gate::rz(ANY_ANGLE, q)` in `from`/`to` patterns to match any Rz angle
/// and carry the captured value through to the replacement.
#[cfg(test)]
pub const ANY_ANGLE: f64 = f64::NAN;

#[cfg(test)]
impl Rule {
    pub fn new(from: Circuit, to: Circuit) -> Self {
        assert_eq!(from.num_qubits, to.num_qubits, "rule circuits must have equal qubit counts");
        assert!(from.num_qubits <= MAX_RULE_QUBITS);
        assert!(from.gates.len() <= MAX_PATTERN_GATES);
        Rule { from, to }
    }
}

#[cfg(test)]
pub struct RulePass<'a> {
    pub rules: &'a [Rule],
}

#[cfg(test)]
impl Pass for RulePass<'_> {
    fn name(&self) -> &str { "rule-pass" }
    fn run_with_progress(&self, circuit: &Circuit, _pb: &indicatif::ProgressBar) -> Circuit {
        let mut gates = circuit.gates.clone();
        loop {
            let mut changed = false;
            let mut result = Vec::with_capacity(gates.len());
            let mut i = 0;
            while i < gates.len() {
                let mut matched = false;
                for rule in self.rules {
                    let from = &rule.from.gates;
                    if from.is_empty() { continue; }
                    if let Some(m) = try_match(&gates, from, i, rule.from.num_qubits) {
                        let last_matched = m.matched[m.match_count - 1];
                        for j in i..=last_matched {
                            if !m.matched[..m.match_count].contains(&j) {
                                result.push(gates[j].clone());
                            }
                        }
                        for to_gate in &rule.to.gates {
                            result.extend(remap_gate(to_gate, &m.qubit_map, &m.angles));
                        }
                        i = last_matched + 1;
                        matched = true;
                        changed = true;
                        break;
                    }
                }
                if !matched {
                    result.push(gates[i].clone());
                    i += 1;
                }
            }
            if !changed { break; }
            gates = result;
        }
        let mut out = Circuit::new(circuit.num_qubits);
        for g in gates {
            out.apply(g);
        }
        out
    }
}

/// Result of a successful match: qubit mapping, matched gate indices,
/// and captured Rz angles (keyed by rule qubit for NaN wildcards).
#[cfg(test)]
struct Match {
    qubit_map: [Qubit; MAX_RULE_QUBITS],
    matched: [usize; MAX_PATTERN_GATES],
    match_count: usize,
    angles: [Option<f64>; MAX_RULE_QUBITS],
}

/// Try to match `from` gates against `gates[start..]`, building a qubit mapping.
///
/// Between matched gates, any intervening gate that touches a mapped qubit
/// blocks the match (we can't commute past it).
///
/// Rz gates in `from` with NaN angle act as wildcards — they match any Rz
/// and capture the angle, keyed by rule qubit.
#[cfg(test)]
fn try_match(
    gates: &[Gate],
    from: &[Gate],
    start: usize,
    num_rule_qubits: usize,
) -> Option<Match> {
    let mut map = [None::<Qubit>; MAX_RULE_QUBITS];
    let mut matched = [0usize; MAX_PATTERN_GATES];
    let mut match_count = 0;
    let mut angles = [None::<f64>; MAX_RULE_QUBITS];
    let mut pos = start;

    for from_gate in from {
        loop {
            if pos >= gates.len() {
                return None;
            }
            if match_gate(from_gate, &gates[pos], &mut map[..num_rule_qubits], &mut angles) {
                matched[match_count] = pos;
                match_count += 1;
                pos += 1;
                break;
            }
            // This gate didn't match — check if it touches any mapped qubit
            if gate_touches_mapped(&gates[pos], &map[..num_rule_qubits]) {
                return None;
            }
            pos += 1;
        }
    }

    // Check all rule qubits were bound
    let mut qubit_map = [0usize; MAX_RULE_QUBITS];
    for i in 0..num_rule_qubits {
        qubit_map[i] = map[i]?;
    }
    Some(Match { qubit_map, matched, match_count, angles })
}

/// Try to match a single from_gate against a circuit gate, extending the qubit mapping.
/// Returns true if the gate types match and qubit bindings are consistent.
/// For Rz with NaN angle (wildcard), captures the matched angle keyed by rule qubit.
#[cfg(test)]
fn match_gate(
    from_gate: &Gate,
    circuit_gate: &Gate,
    map: &mut [Option<Qubit>],
    angles: &mut [Option<f64>],
) -> bool {
    match (from_gate, circuit_gate) {
        (Gate::x(a), Gate::x(b))
        | (Gate::h(a), Gate::h(b))
        | (Gate::s(a), Gate::s(b))
        | (Gate::sdg(a), Gate::sdg(b))
        | (Gate::z(a), Gate::z(b))
        | (Gate::t(a), Gate::t(b))
        | (Gate::tdg(a), Gate::tdg(b)) => try_bind(map, &[(*a, *b)]),

        (Gate::rz(theta_a, a), _) if theta_a.is_nan() => {
            // Wildcard: match any z-axis rotation, capture the angle
            let Some((b, angle)) = z_angle(circuit_gate) else { return false; };
            if try_bind(map, &[(*a, b)]) {
                angles[*a] = Some(angle);
                true
            } else {
                false
            }
        }

        (Gate::rz(theta_a, a), Gate::rz(theta_b, b)) => {
            (theta_a - theta_b).abs() < 1e-10 && try_bind(map, &[(*a, *b)])
        }

        (
            Gate::cnot { control: ac, target: at },
            Gate::cnot { control: bc, target: bt },
        ) => try_bind(map, &[(*ac, *bc), (*at, *bt)]),

        (
            Gate::ccx { control1: a1, control2: a2, target: at },
            Gate::ccx { control1: b1, control2: b2, target: bt },
        ) => try_bind(map, &[(*a1, *b1), (*a2, *b2), (*at, *bt)]),

        _ => false,
    }
}

/// Try to bind rule_qubit -> circuit_qubit for each pair. Injective: different rule
/// qubits must map to different circuit qubits. Rolls back on failure.
#[cfg(test)]
fn try_bind(
    map: &mut [Option<Qubit>],
    pairs: &[(Qubit, Qubit)],
) -> bool {
    let mut saved = [None::<Qubit>; MAX_RULE_QUBITS];
    saved[..map.len()].copy_from_slice(map);

    for &(rule_q, circ_q) in pairs {
        match map[rule_q] {
            Some(existing) => {
                if existing != circ_q {
                    map.copy_from_slice(&saved[..map.len()]);
                    return false;
                }
            }
            None => {
                // Check injectivity: circ_q must not already be bound to a different rule_q
                if map.iter().enumerate().any(|(i, &bound)| bound == Some(circ_q) && i != rule_q) {
                    map.copy_from_slice(&saved[..map.len()]);
                    return false;
                }
                map[rule_q] = Some(circ_q);
            }
        }
    }
    true
}

#[cfg(test)]
fn gate_touches_mapped(gate: &Gate, map: &[Option<Qubit>]) -> bool {
    let (count, qs) = qubits_of(gate);
    qs[..count].iter().any(|cq| map.iter().any(|bound| *bound == Some(*cq)))
}

#[cfg(test)]
fn qubits_of(gate: &Gate) -> (usize, [Qubit; 3]) {
    match gate {
        Gate::x(q) | Gate::h(q) | Gate::s(q) | Gate::sdg(q) | Gate::z(q)
        | Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => (1, [*q, 0, 0]),
        Gate::cnot { control, target } => (2, [*control, *target, 0]),
        Gate::ccx { control1, control2, target } => (3, [*control1, *control2, *target]),
    }
}

#[cfg(test)]
fn remap_gate(gate: &Gate, map: &[Qubit; MAX_RULE_QUBITS], angles: &[Option<f64>; MAX_RULE_QUBITS]) -> Vec<Gate> {
    match gate {
        Gate::x(q) => vec![Gate::x(map[*q])],
        Gate::h(q) => vec![Gate::h(map[*q])],
        Gate::s(q) => vec![Gate::s(map[*q])],
        Gate::sdg(q) => vec![Gate::sdg(map[*q])],
        Gate::z(q) => vec![Gate::z(map[*q])],
        Gate::t(q) => vec![Gate::t(map[*q])],
        Gate::tdg(q) => vec![Gate::tdg(map[*q])],
        Gate::rz(theta, q) => {
            if theta.is_nan() {
                let angle = angles[*q].expect("wildcard Rz in `to` has no captured angle");
                emit_z_rotation(angle, map[*q])
            } else {
                vec![Gate::rz(*theta, map[*q])]
            }
        }
        Gate::cnot { control, target } => vec![Gate::cnot { control: map[*control], target: map[*target] }],
        Gate::ccx { control1, control2, target } => vec![Gate::ccx {
            control1: map[*control1], control2: map[*control2], target: map[*target],
        }],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pass::Pass;
    use crate::unitary::circuits_equiv;

    fn make_circuit(num_qubits: usize, gates: Vec<Gate>) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for g in gates { c.apply(g); }
        c
    }

    #[test]
    fn hh_cancel() {
        let rule = Rule::new(
            make_circuit(1, vec![Gate::h(0), Gate::h(0)]),
            make_circuit(1, vec![]),
        );
        let pass = RulePass { rules: &[rule] };

        let circuit = make_circuit(1, vec![Gate::h(0), Gate::h(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 0);
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn hh_cancel_on_different_qubit() {
        let rule = Rule::new(
            make_circuit(1, vec![Gate::h(0), Gate::h(0)]),
            make_circuit(1, vec![]),
        );
        let pass = RulePass { rules: &[rule] };

        let circuit = make_circuit(3, vec![Gate::h(2), Gate::h(2)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 0);
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn hh_cancel_skips_unrelated_gates() {
        let rule = Rule::new(
            make_circuit(1, vec![Gate::h(0), Gate::h(0)]),
            make_circuit(1, vec![]),
        );
        let pass = RulePass { rules: &[rule] };

        // H q0; T q1; H q0 — T on q1 doesn't interfere
        let circuit = make_circuit(2, vec![Gate::h(0), Gate::t(1), Gate::h(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 1); // just T q1
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn hh_blocked_by_intervening_gate() {
        let rule = Rule::new(
            make_circuit(1, vec![Gate::h(0), Gate::h(0)]),
            make_circuit(1, vec![]),
        );
        let pass = RulePass { rules: &[rule] };

        // H q0; T q0; H q0 — T on same qubit blocks
        let circuit = make_circuit(1, vec![Gate::h(0), Gate::t(0), Gate::h(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 3); // unchanged
    }

    #[test]
    fn cnot_flip_rule() {
        // H q0; H q1; CNOT q0,q1; H q0; H q1 → CNOT q1,q0
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::h(0), Gate::h(1),
                Gate::cnot { control: 0, target: 1 },
                Gate::h(0), Gate::h(1),
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 1, target: 0 },
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        let circuit = make_circuit(2, vec![
            Gate::h(0), Gate::h(1),
            Gate::cnot { control: 0, target: 1 },
            Gate::h(0), Gate::h(1),
        ]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 1);
        assert!(matches!(&result.gates[0], Gate::cnot { control: 1, target: 0 }));
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn cnot_flip_on_other_qubits() {
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::h(0), Gate::h(1),
                Gate::cnot { control: 0, target: 1 },
                Gate::h(0), Gate::h(1),
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 1, target: 0 },
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // Same pattern but on qubits 3 and 5
        let circuit = make_circuit(6, vec![
            Gate::h(3), Gate::h(5),
            Gate::cnot { control: 3, target: 5 },
            Gate::h(3), Gate::h(5),
        ]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 1);
        assert!(matches!(&result.gates[0], Gate::cnot { control: 5, target: 3 }));
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn multiple_applications() {
        let rule = Rule::new(
            make_circuit(1, vec![Gate::h(0), Gate::h(0)]),
            make_circuit(1, vec![]),
        );
        let pass = RulePass { rules: &[rule] };

        let circuit = make_circuit(1, vec![Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 0);
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn no_match_preserves_circuit() {
        let rule = Rule::new(
            make_circuit(1, vec![Gate::h(0), Gate::h(0)]),
            make_circuit(1, vec![]),
        );
        let pass = RulePass { rules: &[rule] };

        let circuit = make_circuit(1, vec![Gate::t(0), Gate::h(0), Gate::s(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 3);
    }

    #[test]
    fn replacement_with_more_gates() {
        // Replace S with T T
        let rule = Rule::new(
            make_circuit(1, vec![Gate::s(0)]),
            make_circuit(1, vec![Gate::t(0), Gate::t(0)]),
        );
        let pass = RulePass { rules: &[rule] };

        let circuit = make_circuit(1, vec![Gate::s(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 2);
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn injective_mapping_enforced() {
        // Rule has 2 qubits — they must map to distinct circuit qubits
        let rule = Rule::new(
            make_circuit(2, vec![Gate::h(0), Gate::h(1)]),
            make_circuit(2, vec![]),
        );
        let pass = RulePass { rules: &[rule] };

        // H q0; H q0 — both H's on same qubit, can't map rule q0≠q1 to same qubit
        let circuit = make_circuit(1, vec![Gate::h(0), Gate::h(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 2); // no match
    }

    #[test]
    fn rz_commutes_past_cnot_control() {
        // Rz q0; CNOT q0,q1 → CNOT q0,q1; Rz q0
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::rz(ANY_ANGLE, 0),
                Gate::cnot { control: 0, target: 1 },
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 0, target: 1 },
                Gate::rz(ANY_ANGLE, 0),
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // Single application
        let circuit = make_circuit(2, vec![
            Gate::rz(0.5, 0),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 2);
        assert!(matches!(&result.gates[0], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&result.gates[1], Gate::rz(theta, 0) if (*theta - 0.5).abs() < 1e-10));
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn rz_propagates_through_cnot_chain() {
        // Rz should commute past each CNOT where it's on the control
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::rz(ANY_ANGLE, 0),
                Gate::cnot { control: 0, target: 1 },
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 0, target: 1 },
                Gate::rz(ANY_ANGLE, 0),
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // Rz q0; CNOT q0,q1; CNOT q0,q2; CNOT q0,q3
        // After one pass: CNOT q0,q1; CNOT q0,q2; CNOT q0,q3; Rz q0
        let circuit = make_circuit(4, vec![
            Gate::rz(1.23, 0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 2 },
            Gate::cnot { control: 0, target: 3 },
        ]);
        let result = pass.run(&circuit);
        assert!(circuits_equiv(&circuit, &result, 1e-10));
        assert_eq!(result.gates.len(), 4);
        // After one walk: Rz commutes past first CNOT, then past second, then past third
        // Result: CNOT q0,q1; CNOT q0,q2; CNOT q0,q3; Rz q0
        assert!(matches!(&result.gates[0], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&result.gates[1], Gate::cnot { control: 0, target: 2 }));
        assert!(matches!(&result.gates[2], Gate::cnot { control: 0, target: 3 }));
        assert!(matches!(&result.gates[3], Gate::rz(theta, 0) if (*theta - 1.23).abs() < 1e-10));
    }

    #[test]
    fn rz_stops_at_target_cnot() {
        // Rz should NOT commute past a CNOT where it's on the target
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::rz(ANY_ANGLE, 0),
                Gate::cnot { control: 0, target: 1 },
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 0, target: 1 },
                Gate::rz(ANY_ANGLE, 0),
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // Rz q0; CNOT q1,q0 — Rz is on the target, rule shouldn't match
        let circuit = make_circuit(2, vec![
            Gate::rz(0.5, 0),
            Gate::cnot { control: 1, target: 0 },
        ]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 2);
        assert!(matches!(&result.gates[0], Gate::rz(..)));
        assert!(matches!(&result.gates[1], Gate::cnot { .. }));
    }

    #[test]
    fn rz_wildcard_preserves_angle() {
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::rz(ANY_ANGLE, 0),
                Gate::cnot { control: 0, target: 1 },
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 0, target: 1 },
                Gate::rz(ANY_ANGLE, 0),
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // Two different angles on different qubits
        let circuit = make_circuit(3, vec![
            Gate::rz(0.7, 0),
            Gate::cnot { control: 0, target: 1 },
            Gate::rz(1.4, 0),
            Gate::cnot { control: 0, target: 2 },
        ]);
        let result = pass.run(&circuit);
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn t_commutes_past_cnot_control() {
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::rz(ANY_ANGLE, 0),
                Gate::cnot { control: 0, target: 1 },
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 0, target: 1 },
                Gate::rz(ANY_ANGLE, 0),
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // T q0; CNOT q0,q1 → CNOT q0,q1; T q0
        let circuit = make_circuit(2, vec![
            Gate::t(0),
            Gate::cnot { control: 0, target: 1 },
        ]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 2);
        assert!(matches!(&result.gates[0], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&result.gates[1], Gate::t(0)));
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn tdg_propagates_through_cnot_chain() {
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::rz(ANY_ANGLE, 0),
                Gate::cnot { control: 0, target: 1 },
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 0, target: 1 },
                Gate::rz(ANY_ANGLE, 0),
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // Tdg q0; CNOT q0,q1; CNOT q0,q2
        let circuit = make_circuit(3, vec![
            Gate::tdg(0),
            Gate::cnot { control: 0, target: 1 },
            Gate::cnot { control: 0, target: 2 },
        ]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 3);
        assert!(matches!(&result.gates[0], Gate::cnot { .. }));
        assert!(matches!(&result.gates[1], Gate::cnot { .. }));
        assert!(matches!(&result.gates[2], Gate::tdg(0)));
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }

    #[test]
    fn s_commutes_z_commutes() {
        let rule = Rule::new(
            make_circuit(2, vec![
                Gate::rz(ANY_ANGLE, 0),
                Gate::cnot { control: 0, target: 1 },
            ]),
            make_circuit(2, vec![
                Gate::cnot { control: 0, target: 1 },
                Gate::rz(ANY_ANGLE, 0),
            ]),
        );
        let pass = RulePass { rules: &[rule] };

        // S q0; CNOT q0,q1 → CNOT q0,q1; S q0
        let c1 = make_circuit(2, vec![Gate::s(0), Gate::cnot { control: 0, target: 1 }]);
        let r1 = pass.run(&c1);
        assert!(matches!(&r1.gates[0], Gate::cnot { .. }));
        assert!(matches!(&r1.gates[1], Gate::s(0)));
        assert!(circuits_equiv(&c1, &r1, 1e-10));

        // Z q0; CNOT q0,q1 → CNOT q0,q1; Z q0
        let c2 = make_circuit(2, vec![Gate::z(0), Gate::cnot { control: 0, target: 1 }]);
        let r2 = pass.run(&c2);
        assert!(matches!(&r2.gates[0], Gate::cnot { .. }));
        assert!(matches!(&r2.gates[1], Gate::z(0)));
        assert!(circuits_equiv(&c2, &r2, 1e-10));

        // Sdg q0; CNOT q0,q1 → CNOT q0,q1; Sdg q0
        let c3 = make_circuit(2, vec![Gate::sdg(0), Gate::cnot { control: 0, target: 1 }]);
        let r3 = pass.run(&c3);
        assert!(matches!(&r3.gates[0], Gate::cnot { .. }));
        assert!(matches!(&r3.gates[1], Gate::sdg(0)));
        assert!(circuits_equiv(&c3, &r3, 1e-10));
    }

    #[test]
    fn multiple_rules() {
        let rules = vec![
            Rule::new(
                make_circuit(1, vec![Gate::h(0), Gate::h(0)]),
                make_circuit(1, vec![]),
            ),
            Rule::new(
                make_circuit(1, vec![Gate::x(0), Gate::x(0)]),
                make_circuit(1, vec![]),
            ),
        ];
        let pass = RulePass { rules: &rules };

        let circuit = make_circuit(1, vec![Gate::h(0), Gate::h(0), Gate::x(0), Gate::x(0)]);
        let result = pass.run(&circuit);
        assert_eq!(result.gates.len(), 0);
        assert!(circuits_equiv(&circuit, &result, 1e-10));
    }
}
