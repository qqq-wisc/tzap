use crate::angle::Angle;
use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hasher};

use rustc_hash::FxHashMap;

use crate::circuit::{Circuit, Gate, Qubit};
use crate::pass::Pass;

/// A random evaluation of a Boolean polynomial in GF(2^128), together with
/// an upper bound on its total degree. Addition is XOR; multiplication is
/// carry-less multiplication modulo x^128 + x^7 + x^2 + x + 1.
#[derive(Clone, Copy)]
struct Fingerprint {
    value: u128,
    degree: u64,
}

const FIELD_ONE: u128 = 1;
/// Keeps the Schwartz–Zippel false-match bound at or below 2^-96 for each
/// comparison in the 128-bit field.
const MAX_TRACKED_DEGREE: u64 = 1 << 32;

fn fresh_fingerprint() -> Fingerprint {
    let hi = RandomState::new().build_hasher().finish() as u128;
    let lo = RandomState::new().build_hasher().finish() as u128;
    Fingerprint {
        value: (hi << 64) | lo,
        degree: 1,
    }
}

fn canonical_fingerprint(value: u128) -> (u128, bool) {
    let complement = value ^ FIELD_ONE;
    if value <= complement {
        (value, false)
    } else {
        (complement, true)
    }
}

fn field_mul(mut left: u128, mut right: u128) -> u128 {
    let mut product = 0;
    while right != 0 {
        if right & 1 != 0 {
            product ^= left;
        }
        right >>= 1;
        let carry = left >> 127;
        left <<= 1;
        if carry != 0 {
            left ^= 0x87;
        }
    }
    product
}

struct LivePhase {
    angle: Angle,
    qubit: Qubit,
    current_idx: usize,
    current_sign: bool,
    input_gates: usize,
}

/// Merges T/Rz gates that act on the same Boolean polynomial using random
/// 128-bit field evaluations rather than exact symbolic expressions. CNOT
/// remains linear; CCX multiplies its control fingerprints and introduces
/// nonlinear terms. The tracked degree bounds a false match by
/// `degree / 2^128`; tracking restarts conservatively before that bound becomes
/// too weak.
pub struct PhaseFoldRand;

impl Pass for PhaseFoldRand {
    fn name(&self) -> &str {
        "Phase folding"
    }
    fn run(&self, circuit: &Circuit) -> Circuit {
        phase_fold_rand(circuit)
    }
}

/// Free-function form of [`PhaseFoldRand`], for use outside a [`Pass`] pipeline.
pub fn phase_fold_rand(circuit: &Circuit) -> Circuit {
    let n = circuit.num_qubits;
    let mut qubits: Vec<Fingerprint> = (0..n).map(|_| fresh_fingerprint()).collect();

    let mut live: Vec<LivePhase> = Vec::new();
    let mut fingerprint_to_group: FxHashMap<u128, usize> = FxHashMap::default();
    let mut skip = vec![false; circuit.gates.len()];
    let mut emit_at: Vec<Option<(Qubit, Angle)>> = vec![None; circuit.gates.len()];

    for (idx, gate) in circuit.gates.iter().enumerate() {
        match gate {
            Gate::t(q) => record_int(
                &qubits,
                *q,
                1,
                idx,
                &mut live,
                &mut fingerprint_to_group,
                &mut skip,
            ),
            Gate::tdg(q) => record_int(
                &qubits,
                *q,
                7,
                idx,
                &mut live,
                &mut fingerprint_to_group,
                &mut skip,
            ),
            Gate::s(q) => record_int(
                &qubits,
                *q,
                2,
                idx,
                &mut live,
                &mut fingerprint_to_group,
                &mut skip,
            ),
            Gate::sdg(q) => record_int(
                &qubits,
                *q,
                6,
                idx,
                &mut live,
                &mut fingerprint_to_group,
                &mut skip,
            ),
            Gate::z(q) => record_int(
                &qubits,
                *q,
                4,
                idx,
                &mut live,
                &mut fingerprint_to_group,
                &mut skip,
            ),
            Gate::rz(theta, q) => record_angle(
                &qubits,
                *q,
                theta.clone(),
                idx,
                &mut live,
                &mut fingerprint_to_group,
                &mut skip,
            ),
            Gate::h(q) => {
                qubits[*q as usize] = fresh_fingerprint();
            }
            Gate::x(q) => {
                qubits[*q as usize].value ^= FIELD_ONE;
            }
            Gate::cnot { control, target } => {
                let control = qubits[*control as usize];
                let target = &mut qubits[*target as usize];
                target.value ^= control.value;
                target.degree = target.degree.max(control.degree);
            }
            Gate::cz { .. } | Gate::ccz { .. } => {
                // CZ and CCZ are diagonal: they change phase but leave the tracked
                // computational-basis parities unchanged. The gates themselves remain
                // in the reconstructed circuit.
            }
            Gate::ccx {
                control1,
                control2,
                target,
            } => {
                let left = qubits[*control1 as usize];
                let right = qubits[*control2 as usize];
                let target_state = qubits[*target as usize];
                qubits[*target as usize] = match left.degree.checked_add(right.degree) {
                    Some(product_degree) if product_degree <= MAX_TRACKED_DEGREE => Fingerprint {
                        value: target_state.value ^ field_mul(left.value, right.value),
                        degree: target_state.degree.max(product_degree),
                    },
                    // At extreme degree the Schwartz–Zippel collision bound is
                    // no longer useful. Treat the target's current Boolean
                    // function as a new opaque variable: subsequent equalities
                    // remain sound, while relationships to its old value are
                    // deliberately forgotten.
                    _ => fresh_fingerprint(),
                };
            }
            Gate::measure { .. } => {
                // Computational-basis measurement preserves the measured basis
                // value, so its classical transfer is x' = x. Diagonal rotations
                // commute with the measurement projectors and may still fold
                // across the measurement.
            }
            Gate::reset(qubit) => {
                // Reset overwrites the computational-basis value with |0⟩.
                qubits[*qubit as usize] = Fingerprint {
                    value: 0,
                    degree: 0,
                };
            }
        }
    }

    for lp in &live {
        let angle = if lp.current_sign {
            lp.angle
                .checked_neg()
                .expect("only mergeable angles are tracked")
        } else {
            lp.angle.clone()
        };
        if lp.input_gates == 1 && angle.emitted_gate_count() > 1 {
            continue;
        }
        if angle.is_zero() {
            skip[lp.current_idx] = true;
        } else {
            emit_at[lp.current_idx] = Some((lp.qubit, angle));
        }
    }

    // Reconstruct circuit, moving gates to avoid per-gate clones.
    let mut output = Circuit::with_cbits(n, circuit.num_cbits);
    for (idx, gate) in circuit.gates.clone().into_iter().enumerate() {
        if skip[idx] {
            continue;
        }
        if let Some((qubit, angle)) = &emit_at[idx] {
            angle.emit(&mut output, *qubit);
        } else {
            output.apply(gate);
        }
    }

    output
}

fn record_int(
    qubits: &[Fingerprint],
    q: Qubit,
    k: u8,
    idx: usize,
    live: &mut Vec<LivePhase>,
    groups: &mut FxHashMap<u128, usize>,
    skip: &mut [bool],
) {
    record_angle(qubits, q, Angle::turns(k), idx, live, groups, skip);
}
fn record_angle(
    qubits: &[Fingerprint],
    q: Qubit,
    angle: Angle,
    idx: usize,
    live: &mut Vec<LivePhase>,
    groups: &mut FxHashMap<u128, usize>,
    skip: &mut [bool],
) {
    if angle.is_preserved() {
        return;
    }
    let value = qubits[q as usize].value;
    if value == 0 || value == FIELD_ONE {
        skip[idx] = true;
        return;
    }
    let (key, complement) = canonical_fingerprint(value);
    let signed = if complement {
        angle.checked_neg().expect("validated rotation")
    } else {
        angle
    };
    if let Some(&gi) = groups.get(&key) {
        if let Ok(merged) = live[gi].angle.checked_add_signed(1, &signed) {
            if merged.emitted_gate_count() <= live[gi].input_gates + 1 {
                skip[live[gi].current_idx] = true;
                live[gi].angle = merged;
                live[gi].current_idx = idx;
                live[gi].qubit = q;
                live[gi].current_sign = complement;
                live[gi].input_gates += 1;
                return;
            }
        }
        // Leave the earlier group live for finalization and start a new one.
    }
    let gi = live.len();
    live.push(LivePhase {
        angle: signed,
        qubit: q,
        current_idx: idx,
        current_sign: complement,
        input_gates: 1,
    });
    groups.insert(key, gi);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decompose::DecomposeToffoli;
    use crate::pass::Pass;
    use crate::unitary::circuits_equiv;

    #[test]
    fn forced_sketch_collision_exposes_separate_probabilistic_matching_contract() {
        // Equal sketches are not a proof of equal Boolean functions. This
        // controlled collision documents why numerical exactness alone is
        // not an unconditional correctness guarantee for this pass.
        let qubits = [Fingerprint {
            value: 2,
            degree: 1,
        }; 2];
        let mut live = Vec::new();
        let mut groups = FxHashMap::default();
        let mut skip = [false; 2];
        record_int(&qubits, 0, 1, 0, &mut live, &mut groups, &mut skip);
        record_int(&qubits, 1, 1, 1, &mut live, &mut groups, &mut skip);
        assert!(skip[0]);
        let mut input = Circuit::new(2);
        input.gates = vec![Gate::t(0), Gate::t(1)];
        let mut output = Circuit::new(2);
        live[0].angle.emit(&mut output, live[0].qubit);
        assert!(!circuits_equiv(&input, &output, 1e-12));
    }

    fn count_gates(c: &Circuit) -> usize {
        c.gates.len()
    }

    fn count_phase_gates(c: &Circuit) -> usize {
        c.gates
            .iter()
            .filter(|g| {
                matches!(
                    g,
                    Gate::t(_)
                        | Gate::tdg(_)
                        | Gate::s(_)
                        | Gate::sdg(_)
                        | Gate::z(_)
                        | Gate::rz(..)
                )
            })
            .count()
    }

    fn count_t_gates(c: &Circuit) -> usize {
        c.gates
            .iter()
            .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_)))
            .count()
    }

    #[test]
    #[ignore] // long-running: 24-qubit pipeline benchmark
    fn large_circuit_toffoli_decompose_and_phase_fold() {
        let mut c = Circuit::new(24);
        c.apply(Gate::cnot {
            control: 3,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 8,
            target: 7,
        });
        c.apply(Gate::cnot {
            control: 14,
            target: 13,
        });
        c.apply(Gate::cnot {
            control: 21,
            target: 20,
        });
        c.apply(Gate::cnot {
            control: 3,
            target: 4,
        });
        c.apply(Gate::cnot {
            control: 8,
            target: 9,
        });
        c.apply(Gate::cnot {
            control: 14,
            target: 15,
        });
        c.apply(Gate::cnot {
            control: 21,
            target: 22,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 3,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx {
            control1: 5,
            control2: 6,
            target: 8,
        });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx {
            control1: 11,
            control2: 12,
            target: 14,
        });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx {
            control1: 18,
            control2: 19,
            target: 21,
        });
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(4));
        c.apply(Gate::h(9));
        c.apply(Gate::h(15));
        c.apply(Gate::h(22));
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx {
            control1: 2,
            control2: 3,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(9));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 9,
        });
        c.apply(Gate::h(9));
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 15,
        });
        c.apply(Gate::h(15));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::h(22));
        c.apply(Gate::ccx {
            control1: 20,
            control2: 21,
            target: 22,
        });
        c.apply(Gate::h(22));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx {
            control1: 20,
            control2: 21,
            target: 23,
        });
        c.apply(Gate::h(23));
        c.apply(Gate::cnot {
            control: 6,
            target: 5,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 11,
        });
        c.apply(Gate::cnot {
            control: 19,
            target: 18,
        });
        c.apply(Gate::cnot {
            control: 5,
            target: 8,
        });
        c.apply(Gate::cnot {
            control: 11,
            target: 14,
        });
        c.apply(Gate::cnot {
            control: 18,
            target: 21,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx {
            control1: 20,
            control2: 21,
            target: 23,
        });
        c.apply(Gate::h(23));
        c.apply(Gate::h(4));
        c.apply(Gate::h(10));
        c.apply(Gate::h(15));
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::h(17));
        c.apply(Gate::h(17));
        c.apply(Gate::ccx {
            control1: 16,
            control2: 23,
            target: 17,
        });
        c.apply(Gate::h(17));
        c.apply(Gate::h(22));
        c.apply(Gate::ccx {
            control1: 15,
            control2: 23,
            target: 22,
        });
        c.apply(Gate::h(22));
        c.apply(Gate::h(9));
        c.apply(Gate::ccx {
            control1: 4,
            control2: 10,
            target: 9,
        });
        c.apply(Gate::h(9));
        c.apply(Gate::h(17));
        c.apply(Gate::h(9));
        c.apply(Gate::h(15));
        c.apply(Gate::h(22));
        c.apply(Gate::ccx {
            control1: 9,
            control2: 17,
            target: 22,
        });
        c.apply(Gate::h(22));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx {
            control1: 9,
            control2: 16,
            target: 15,
        });
        c.apply(Gate::h(15));
        c.apply(Gate::h(15));
        c.apply(Gate::h(22));
        c.apply(Gate::h(17));
        c.apply(Gate::h(17));
        c.apply(Gate::ccx {
            control1: 16,
            control2: 23,
            target: 17,
        });
        c.apply(Gate::h(17));
        c.apply(Gate::h(17));
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx {
            control1: 20,
            control2: 21,
            target: 23,
        });
        c.apply(Gate::h(23));
        c.apply(Gate::cnot {
            control: 5,
            target: 8,
        });
        c.apply(Gate::cnot {
            control: 11,
            target: 14,
        });
        c.apply(Gate::cnot {
            control: 18,
            target: 21,
        });
        c.apply(Gate::cnot {
            control: 6,
            target: 5,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 11,
        });
        c.apply(Gate::cnot {
            control: 19,
            target: 18,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx {
            control1: 20,
            control2: 21,
            target: 23,
        });
        c.apply(Gate::h(23));
        c.apply(Gate::h(23));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 3,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx {
            control1: 5,
            control2: 6,
            target: 8,
        });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx {
            control1: 11,
            control2: 12,
            target: 14,
        });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx {
            control1: 18,
            control2: 19,
            target: 21,
        });
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::cnot {
            control: 3,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 8,
            target: 7,
        });
        c.apply(Gate::cnot {
            control: 14,
            target: 13,
        });
        c.apply(Gate::cnot {
            control: 21,
            target: 20,
        });
        c.apply(Gate::cnot {
            control: 6,
            target: 5,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 11,
        });
        c.apply(Gate::cnot {
            control: 19,
            target: 18,
        });
        c.apply(Gate::cnot {
            control: 6,
            target: 8,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 14,
        });
        c.apply(Gate::cnot {
            control: 19,
            target: 21,
        });
        c.apply(Gate::cnot {
            control: 4,
            target: 6,
        });
        c.apply(Gate::cnot {
            control: 9,
            target: 12,
        });
        c.apply(Gate::cnot {
            control: 15,
            target: 19,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 3,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx {
            control1: 5,
            control2: 6,
            target: 8,
        });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx {
            control1: 11,
            control2: 12,
            target: 14,
        });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx {
            control1: 18,
            control2: 19,
            target: 21,
        });
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::cnot {
            control: 3,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 8,
            target: 7,
        });
        c.apply(Gate::cnot {
            control: 14,
            target: 13,
        });
        c.apply(Gate::cnot {
            control: 21,
            target: 20,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 3,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx {
            control1: 5,
            control2: 6,
            target: 8,
        });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx {
            control1: 11,
            control2: 12,
            target: 14,
        });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx {
            control1: 18,
            control2: 19,
            target: 21,
        });
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::cnot {
            control: 6,
            target: 5,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 11,
        });
        c.apply(Gate::cnot {
            control: 19,
            target: 18,
        });
        c.apply(Gate::cnot {
            control: 4,
            target: 6,
        });
        c.apply(Gate::cnot {
            control: 9,
            target: 12,
        });
        c.apply(Gate::cnot {
            control: 15,
            target: 19,
        });
        c.apply(Gate::cnot {
            control: 6,
            target: 8,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 14,
        });
        c.apply(Gate::cnot {
            control: 19,
            target: 21,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::cnot {
            control: 6,
            target: 5,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 11,
        });
        c.apply(Gate::cnot {
            control: 19,
            target: 18,
        });
        c.apply(Gate::x(0));
        c.apply(Gate::x(2));
        c.apply(Gate::x(5));
        c.apply(Gate::x(7));
        c.apply(Gate::x(11));
        c.apply(Gate::x(13));
        c.apply(Gate::cnot {
            control: 3,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 8,
            target: 7,
        });
        c.apply(Gate::cnot {
            control: 14,
            target: 13,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 3,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx {
            control1: 5,
            control2: 6,
            target: 8,
        });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx {
            control1: 11,
            control2: 12,
            target: 14,
        });
        c.apply(Gate::h(14));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::cnot {
            control: 6,
            target: 5,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 11,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::cnot {
            control: 5,
            target: 8,
        });
        c.apply(Gate::cnot {
            control: 11,
            target: 14,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::h(15));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx {
            control1: 9,
            control2: 16,
            target: 15,
        });
        c.apply(Gate::h(15));
        c.apply(Gate::h(9));
        c.apply(Gate::h(9));
        c.apply(Gate::ccx {
            control1: 4,
            control2: 10,
            target: 9,
        });
        c.apply(Gate::h(9));
        c.apply(Gate::h(4));
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::cnot {
            control: 5,
            target: 8,
        });
        c.apply(Gate::cnot {
            control: 11,
            target: 14,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 10,
        });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 16,
        });
        c.apply(Gate::h(16));
        c.apply(Gate::cnot {
            control: 6,
            target: 5,
        });
        c.apply(Gate::cnot {
            control: 12,
            target: 11,
        });
        c.apply(Gate::h(9));
        c.apply(Gate::ccx {
            control1: 7,
            control2: 8,
            target: 9,
        });
        c.apply(Gate::h(9));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx {
            control1: 13,
            control2: 14,
            target: 15,
        });
        c.apply(Gate::h(15));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx {
            control1: 2,
            control2: 3,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::h(9));
        c.apply(Gate::h(10));
        c.apply(Gate::h(15));
        c.apply(Gate::h(16));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 3,
        });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx {
            control1: 5,
            control2: 6,
            target: 8,
        });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx {
            control1: 11,
            control2: 12,
            target: 14,
        });
        c.apply(Gate::h(14));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::cnot {
            control: 3,
            target: 4,
        });
        c.apply(Gate::cnot {
            control: 8,
            target: 9,
        });
        c.apply(Gate::cnot {
            control: 14,
            target: 15,
        });
        c.apply(Gate::cnot {
            control: 3,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 8,
            target: 7,
        });
        c.apply(Gate::cnot {
            control: 14,
            target: 13,
        });
        c.apply(Gate::x(0));
        c.apply(Gate::x(2));
        c.apply(Gate::x(5));
        c.apply(Gate::x(7));
        c.apply(Gate::x(11));
        c.apply(Gate::x(13));

        let original_gates = count_gates(&c);
        let original_t = count_t_gates(&c);

        let dec = DecomposeToffoli.run(&c);
        let dec_t = count_t_gates(&dec);

        let folded = phase_fold_rand(&dec);
        let folded_t = count_t_gates(&folded);

        let folded2 = phase_fold_rand(&folded);
        let folded2_t = count_t_gates(&folded2);

        println!("=== large circuit pipeline ===");
        println!(
            "original:            {} gates, {} T/Tdg",
            original_gates, original_t
        );
        println!(
            "after toffoli decomp: {} gates, {} T/Tdg",
            count_gates(&dec),
            dec_t
        );
        println!(
            "after phase fold 1:  {} gates, {} T/Tdg",
            count_gates(&folded),
            folded_t
        );
        println!(
            "after phase fold 2:  {} gates, {} T/Tdg",
            count_gates(&folded2),
            folded2_t
        );

        assert!(
            folded_t <= dec_t,
            "phase folding should not increase T count"
        );
        assert!(
            folded2_t <= folded_t,
            "second phase fold should not increase T count"
        );
    }

    #[test]
    #[ignore] // long-running: 5-qubit multi-pass pipeline
    fn small_circuit_pipeline() {
        let mut c = Circuit::new(5);
        c.apply(Gate::x(4));
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 3,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx {
            control1: 2,
            control2: 3,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::cnot {
            control: 3,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx {
            control1: 1,
            control2: 2,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::cnot {
            control: 2,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 4,
        });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::cnot {
            control: 1,
            target: 4,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 4,
        });

        let original_gates = count_gates(&c);
        let original_t = count_t_gates(&c);

        let dec = DecomposeToffoli.run(&c);
        let dec_t = count_t_gates(&dec);

        let folded = phase_fold_rand(&dec);
        let folded_t = count_t_gates(&folded);

        let folded2 = phase_fold_rand(&folded);
        let folded2_t = count_t_gates(&folded2);

        println!("=== small circuit pipeline ===");
        println!(
            "original:            {} gates, {} T/Tdg",
            original_gates, original_t
        );
        println!(
            "after toffoli decomp: {} gates, {} T/Tdg",
            count_gates(&dec),
            dec_t
        );
        println!(
            "after phase fold 1:  {} gates, {} T/Tdg",
            count_gates(&folded),
            folded_t
        );
        println!(
            "after phase fold 2:  {} gates, {} T/Tdg",
            count_gates(&folded2),
            folded2_t
        );

        assert!(folded_t <= dec_t);
    }

    #[test]
    fn two_toffoli_shared_control() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx {
            control1: 1,
            control2: 2,
            target: 0,
        });
        c.apply(Gate::cnot {
            control: 2,
            target: 0,
        });
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });

        let dec = DecomposeToffoli.run(&c);
        assert_eq!(count_t_gates(&dec), 14);

        let opt = phase_fold_rand(&dec);
        assert!(circuits_equiv(&c, &opt, 1e-10));
        assert_eq!(count_t_gates(&opt), 12);
    }

    #[test]
    fn nonlinear_ccx_fingerprint_recovers_after_inverse_pair() {
        let mut circuit = Circuit::new(3);
        circuit.apply(Gate::t(2));
        circuit.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        circuit.apply(Gate::ccx {
            control1: 1,
            control2: 0,
            target: 2,
        });
        circuit.apply(Gate::tdg(2));

        let folded = phase_fold_rand(&circuit);
        assert_eq!(
            folded.gates,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccx {
                    control1: 1,
                    control2: 0,
                    target: 2,
                },
            ]
        );
        assert!(circuits_equiv(&circuit, &folded, 1e-10));
    }

    #[test]
    fn gf128_multiplication_obeys_field_identities() {
        let samples = [0, 1, 0x1234_5678_9abc_def0, u128::MAX, 1 << 127];
        for &a in &samples {
            assert_eq!(field_mul(a, 0), 0);
            assert_eq!(field_mul(a, 1), a);
            for &b in &samples {
                assert_eq!(field_mul(a, b), field_mul(b, a));
                for &c in &samples {
                    assert_eq!(field_mul(a, b ^ c), field_mul(a, b) ^ field_mul(a, c));
                }
            }
        }
        assert_eq!(field_mul(1 << 127, 2), 0x87);
    }

    #[test]
    fn nonlinear_shared_control_dependencies_remain_equivalent() {
        let mut circuit = Circuit::new(3);
        circuit.apply(Gate::t(2));
        circuit.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        circuit.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        circuit.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        circuit.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        circuit.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        circuit.apply(Gate::tdg(2));

        let folded = phase_fold_rand(&circuit);
        assert!(circuits_equiv(&circuit, &folded, 1e-10));
    }

    #[test]
    fn extreme_polynomial_degree_conservatively_stops_phase_folding() {
        // Repeatedly update the lowest-degree wire from the other two. The
        // degree bounds then grow like Fibonacci numbers, crossing 2^32 in
        // fewer than 64 CCX gates without needing an enormous circuit.
        let mut circuit = Circuit::new(3);
        let mut degrees = [1u64; 3];
        let mut growth_gates = 0;
        let (control1, control2, target, overflow_degree) = loop {
            let target = (0..3).min_by_key(|&qubit| degrees[qubit]).unwrap();
            let controls: Vec<_> = (0..3).filter(|&qubit| qubit != target).collect();
            let product_degree = degrees[controls[0]] + degrees[controls[1]];
            if product_degree > MAX_TRACKED_DEGREE {
                break (
                    controls[0] as Qubit,
                    controls[1] as Qubit,
                    target as Qubit,
                    product_degree,
                );
            }
            circuit.apply(Gate::ccx {
                control1: controls[0] as Qubit,
                control2: controls[1] as Qubit,
                target: target as Qubit,
            });
            degrees[target] = degrees[target].max(product_degree);
            growth_gates += 1;
        };

        assert!(growth_gates < 64);
        assert!(overflow_degree > MAX_TRACKED_DEGREE);
        assert!(degrees.iter().all(|&degree| degree <= MAX_TRACKED_DEGREE));

        // The identical pair is logically the identity, so exact tracking
        // could cancel T/Tdg across it. Both CCXs exceed the safe degree bound,
        // however, and each refreshes the target fingerprint. Keeping both
        // phase gates is the required conservative result; removing them here
        // would mean the cutoff was not being used.
        circuit.apply(Gate::t(target));
        circuit.apply(Gate::ccx {
            control1,
            control2,
            target,
        });
        circuit.apply(Gate::ccx {
            control1,
            control2,
            target,
        });
        circuit.apply(Gate::tdg(target));

        let folded = phase_fold_rand(&circuit);
        assert_eq!(count_phase_gates(&folded), 2);
        assert_eq!(folded.gates, circuit.gates);
        assert!(circuits_equiv(&circuit, &folded, 1e-10));
    }

    #[test]
    #[ignore] // long-running: combinatorial CCX removal search
    fn mod5_4_remove_ccx_combos() {
        let ccx_indices: Vec<usize> = vec![1, 2, 4, 6];
        let names = [
            "ccx q0,q3,q4",
            "ccx q2,q3,q4",
            "ccx q1,q2,q4",
            "ccx q0,q1,q4",
        ];

        let build = |skip: &[usize]| -> Circuit {
            let mut c = Circuit::new(5);
            let gates: Vec<Gate> = vec![
                Gate::x(4),
                Gate::ccx {
                    control1: 0,
                    control2: 3,
                    target: 4,
                },
                Gate::ccx {
                    control1: 2,
                    control2: 3,
                    target: 4,
                },
                Gate::cnot {
                    control: 3,
                    target: 4,
                },
                Gate::ccx {
                    control1: 1,
                    control2: 2,
                    target: 4,
                },
                Gate::cnot {
                    control: 2,
                    target: 4,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 4,
                },
                Gate::cnot {
                    control: 1,
                    target: 4,
                },
                Gate::cnot {
                    control: 0,
                    target: 4,
                },
            ];
            for (j, g) in gates.into_iter().enumerate() {
                if !skip.contains(&j) {
                    c.apply(g);
                }
            }
            c
        };

        println!("\n{:<50} {:>5} {:>5}", "removed", "dec_T", "ours");
        println!("{}", "-".repeat(65));

        for keep in 0..4 {
            let skip: Vec<usize> = (0..4)
                .filter(|&i| i != keep)
                .map(|i| ccx_indices[i])
                .collect();
            let c = build(&skip);
            let dec = DecomposeToffoli.run(&c);
            let dec_t = count_t_gates(&dec);
            let opt = phase_fold_rand(&dec);
            let opt_t = count_t_gates(&opt);
            assert!(circuits_equiv(&dec, &opt, 1e-10));
            let removed: Vec<&str> = (0..4).filter(|&i| i != keep).map(|i| names[i]).collect();
            println!("{:<50} {:>5} {:>5}", removed.join(" + "), dec_t, opt_t);
        }
    }

    #[test]
    fn toffoli_decompose_then_phase_fold() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        let dec = DecomposeToffoli.run(&c);
        let opt = phase_fold_rand(&dec);
        let dec_phases = count_phase_gates(&dec);
        let opt_phases = count_phase_gates(&opt);
        assert!(opt_phases <= dec_phases);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn two_t_merge_to_s() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_and_tdg_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn four_t_merge_to_s_s() {
        let mut c = Circuit::new(1);
        for _ in 0..4 {
            c.apply(Gate::t(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn eight_t_cancel() {
        let mut c = Circuit::new(1);
        for _ in 0..8 {
            c.apply(Gate::t(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn same_parity_across_cnot() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn different_parity_no_merge() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn h_prevents_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn merge_across_x() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        // T·X·T = e^{iπ/4}·X, so both Ts fold into global phase.
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_x_tdg_folds_across_x() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        // T; X; Tdg ≡ Sdg·X up to global phase — one phase gate survives.
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_folds_across_x() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz_f64(0.3, 0).unwrap());
        c.apply(Gate::x(0));
        c.apply(Gate::rz_f64(0.7, 0).unwrap());
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_cancels_across_x() {
        // RZ(θ); X; RZ(θ) = RZ(θ)·RZ(−θ)·X = X up to global phase.
        let mut c = Circuit::new(1);
        c.apply(Gate::rz_f64(0.42, 0).unwrap());
        c.apply(Gate::x(0));
        c.apply(Gate::rz_f64(0.42, 0).unwrap());
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn triple_t_with_two_xs_folds_to_single_t() {
        // T·X·T·X·T = T up to global phase (T·X·T folds to e^{iπ/4}·X,
        // leaving X·X·T ≡ T up to global phase).
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn x_then_t_then_x_then_t_cancels_to_identity() {
        // X·T·X·T = X·(X·T) up to sign = identity up to global phase.
        let mut c = Circuit::new(1);
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        // Both Ts fold to global phase; Xs stay because this pass doesn't cancel them.
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn mixed_int_and_float_across_x() {
        // T on p, then rz(non-π/4) on ¬p — verifies int+float mixing on one group.
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::rz_f64(0.5, 0).unwrap());
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_x_z_x_cancels_to_identity() {
        // Z·X·Z·X = −I up to global phase — all phase gates fold away.
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::x(0));
        c.apply(Gate::z(0));
        c.apply(Gate::x(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cnot_target_x_sandwich_folds() {
        // T(0); CNOT(0,1); X(1); CNOT(0,1); T(0) — the CNOT·X(target)·CNOT
        // sandwich simplifies to X(1), so the two Ts fold across it as if on
        // the same parity.
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::x(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 0);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cnot_propagates_negation_from_control() {
        // X on control before CNOT negates the target's parity (!a ^ b = !(a^b)).
        // T(1); CNOT(0,1); X(0); CNOT(0,1); T(1) — second T on q1 sees
        // complement of first T's parity.
        let mut c = Circuit::new(2);
        c.apply(Gate::t(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::x(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));
        let opt = phase_fold_rand(&c);
        // Both Ts live on parities that are complements of each other,
        // so they fold (to a single phase gate on q1 — T + (−T) = 0, then X flips
        // the relationship, so the concrete count depends on the emitted rotation).
        assert!(count_t_gates(&opt) <= 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn complement_then_direct_hit() {
        // Exercises the path where the group is first reached via the
        // complement branch, then later via the direct branch.
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0)); // group[p0] int=1, sign=F
        c.apply(Gate::x(0)); // qubits[0] = !p0
        c.apply(Gate::t(0)); // complement hit → int=0, sign=T
        c.apply(Gate::x(0)); // qubits[0] = p0
        c.apply(Gate::tdg(0)); // direct hit on p0 → int=(0-1)&7=7, sign=F
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn h_still_blocks_complement_fold() {
        // H refreshes to a random hash; neither it nor its complement
        // should collide with any prior group. T on each side must survive.
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_qubit_folding() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn toffoli_decomposition_fold() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(2));
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::tdg(2));
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::t(2));
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::tdg(2));
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::t(1));
        c.apply(Gate::t(2));
        c.apply(Gate::h(2));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });

        let original_phases = count_phase_gates(&c);
        let opt = phase_fold_rand(&c);
        let optimized_phases = count_phase_gates(&opt);

        println!("toffoli: {original_phases} phase gates -> {optimized_phases} phase gates");
        assert!(optimized_phases <= original_phases);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn preserves_non_phase_structure() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::x(1));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 3);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz_f64(0.25, 0).unwrap());
        c.apply(Gate::rz_f64(0.75, 0).unwrap());
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_merge_via_cnot() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::t(0));

        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_cancel() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::tdg(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::t(0));

        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_three_way() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::rz_f64(0.5, 2).unwrap());
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::rz_f64(0.5, 2).unwrap());

        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_rz_different_qubits_same_parity() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::rz_f64(0.25, 1).unwrap());
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::rz_f64(0.75, 0).unwrap());

        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    fn print_before_after(name: &str, original: &Circuit, optimized: &Circuit) {
        println!("=== {name} ===");
        println!(
            "BEFORE ({} gates, {} phase):",
            count_gates(original),
            count_phase_gates(original)
        );
        print!("{original}");
        println!(
            "AFTER ({} gates, {} phase):",
            count_gates(optimized),
            count_phase_gates(optimized)
        );
        print!("{optimized}");
        println!();
    }

    #[test]
    fn circuit_from_diagram() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::t(2));
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::tdg(0));
        c.apply(Gate::cnot {
            control: 2,
            target: 0,
        });

        let opt = phase_fold_rand(&c);
        print_before_after("diagram circuit (T on y, T† on x)", &c, &opt);

        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cx_t_cx_cx_tdg_cx() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::tdg(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });

        let opt = phase_fold_rand(&c);
        print_before_after("cx T cx cx Tdg cx", &c, &opt);

        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_swap_h_swap_t() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::h(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));

        let opt = phase_fold_rand(&c);
        print_before_after("T, SWAP, H, SWAP, T", &c, &opt);

        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(opt.gates.last().unwrap(), Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn demo_cross_qubit_all() {
        let mut c1 = Circuit::new(2);
        c1.apply(Gate::t(0));
        c1.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c1.apply(Gate::t(1));
        c1.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c1.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c1.apply(Gate::t(0));
        let opt1 = phase_fold_rand(&c1);
        print_before_after("cross-qubit merge (T on q0 + T on q1 -> S)", &c1, &opt1);

        let mut c2 = Circuit::new(2);
        c2.apply(Gate::t(0));
        c2.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c2.apply(Gate::tdg(1));
        c2.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c2.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c2.apply(Gate::t(0));
        let opt2 = phase_fold_rand(&c2);
        print_before_after("cross-qubit cancel (T on q0 + Tdg on q1 -> 0)", &c2, &opt2);

        let mut c3 = Circuit::new(3);
        c3.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c3.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c3.apply(Gate::rz_f64(0.5, 2).unwrap());
        c3.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c3.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c3.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c3.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c3.apply(Gate::rz_f64(0.5, 2).unwrap());
        let opt3 = phase_fold_rand(&c3);
        print_before_after("3-qubit merge (Rz on q2 twice, same parity)", &c3, &opt3);

        let mut c4 = Circuit::new(2);
        c4.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c4.apply(Gate::rz_f64(0.3, 1).unwrap());
        c4.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c4.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c4.apply(Gate::rz_f64(0.7, 0).unwrap());
        let opt4 = phase_fold_rand(&c4);
        print_before_after("Rz(0.3) on q1 + Rz(0.7) on q0 -> Rz(1.0)", &c4, &opt4);
    }

    // --- z and sdg gate tests ---

    #[test]
    fn z_is_phase_gate() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_is_phase_gate() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_z_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_sdg_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_s_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::s(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_s_is_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::s(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_sdg_is_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_t_t_t_is_z() {
        let mut c = Circuit::new(1);
        for _ in 0..4 {
            c.apply(Gate::t(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_tdg_is_z_plus_t() {
        let mut c = Circuit::new(1);
        for _ in 0..3 {
            c.apply(Gate::tdg(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_t_is_z_plus_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_tdg_is_s_plus_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2); // S + T (3π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_t_is_tdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::tdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_t_folds() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2); // S + T (3π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_tdg_is_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::t(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_tdg_is_z_plus_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_s_is_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::s(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_sdg_is_s() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn six_t_is_sdg() {
        let mut c = Circuit::new(1);
        for _ in 0..6 {
            c.apply(Gate::t(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn seven_t_is_tdg() {
        let mut c = Circuit::new(1);
        for _ in 0..7 {
            c.apply(Gate::t(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::tdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_across_cnot_folds() {
        let mut c = Circuit::new(2);
        c.apply(Gate::z(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_across_cnot_folds() {
        let mut c = Circuit::new(2);
        c.apply(Gate::sdg(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::s(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_h_prevents_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::h(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_h_prevents_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::h(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_x_z_folds_across_x() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::x(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        // Z·X·Z = −X — both Zs fold into global phase.
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_z_merge() {
        let mut c = Circuit::new(2);
        c.apply(Gate::z(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        let opt = phase_fold_rand(&c);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_sdg_s_cancel() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::sdg(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::s(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_in_toffoli_decomp_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::z(2));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::z(2));
        let dec = DecomposeToffoli.run(&c);
        let pf1 = phase_fold_rand(&dec);
        let result = phase_fold_rand(&pf1);
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    #[test]
    fn sdg_in_toffoli_decomp_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::sdg(0));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::s(1));
        let dec = DecomposeToffoli.run(&c);
        let pf1 = phase_fold_rand(&dec);
        let result = phase_fold_rand(&pf1);
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    #[test]
    fn multi_qubit_z_sdg_fold() {
        let mut c = Circuit::new(2);
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(1));
        c.apply(Gate::s(0));
        c.apply(Gate::t(1));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_s_is_sdg() {
        let mut c = Circuit::new(1);
        for _ in 0..3 {
            c.apply(Gate::s(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_sdg_is_s() {
        let mut c = Circuit::new(1);
        for _ in 0..3 {
            c.apply(Gate::sdg(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn four_s_cancel() {
        let mut c = Circuit::new(1);
        for _ in 0..4 {
            c.apply(Gate::s(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn four_sdg_cancel() {
        let mut c = Circuit::new(1);
        for _ in 0..4 {
            c.apply(Gate::sdg(0));
        }
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_t_tdg_is_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn all_phase_types_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        c.apply(Gate::s(0));
        c.apply(Gate::sdg(0));
        c.apply(Gate::z(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn mixed_z_sdg_cnot_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::z(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::sdg(1));
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::t(2));
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::s(0));
        c.apply(Gate::tdg(2));
        let opt = phase_fold_rand(&c);
        assert!(circuits_equiv(&c, &opt, 1e-10));
        assert!(count_phase_gates(&opt) <= count_phase_gates(&c));
    }

    #[test]
    fn s_z_is_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_z_is_s() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_pi_folds_to_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(crate::angle_expr::parse("pi", 1).unwrap(), 0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_neg_half_pi_folds_to_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(
            crate::angle_expr::parse("-pi / 2.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_three_half_pi_folds_to_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(
            crate::angle_expr::parse("3.0 * pi / 2.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_preserves_non_phase_structure() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::z(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::sdg(1));
        c.apply(Gate::x(1));
        let opt = phase_fold_rand(&c);
        assert!(circuits_equiv(&c, &opt, 1e-10));
        let h_count = opt.gates.iter().filter(|g| matches!(g, Gate::h(_))).count();
        let x_count = opt.gates.iter().filter(|g| matches!(g, Gate::x(_))).count();
        let cx_count = opt
            .gates
            .iter()
            .filter(|g| matches!(g, Gate::cnot { .. }))
            .count();
        assert_eq!(h_count, 1);
        assert_eq!(x_count, 1);
        assert_eq!(cx_count, 1);
    }

    #[test]
    fn ccx_z_sdg_full_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(1));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::z(2));
        c.apply(Gate::s(1));
        let dec = DecomposeToffoli.run(&c);
        let pf1 = phase_fold_rand(&dec);
        let result = phase_fold_rand(&pf1);
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    // --- exact pi coefficients and finite residual accumulation tests ---
    // Named phases and explicitly symbolic Rz share exact quarter turns.
    // Numeric residuals remain separate when emitting the combined angle.

    #[test]
    fn t_plus_rz_quarter_pi_is_s() {
        // T (π/4) + rz(π/4) → 2·π/4 = S
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_plus_rz_half_pi_is_z() {
        // S (π/2) + rz(π/2) → π = Z
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 2.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_pi_plus_tdg_is_z_plus_t() {
        // rz(π) (4) + Tdg (7) → 11 & 7 = 3 → S + T
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(crate::angle_expr::parse("pi", 1).unwrap(), 0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2); // S + T
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_tdg_rz_quarter_pi_is_t() {
        // T + Tdg cancels, rz(π/4) → T
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::t(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_plus_rz_irrational_folds() {
        // Preserve the exact T and numeric residual as separate operations.
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz_f64(0.3, 0).unwrap());
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(opt.gates.iter().any(|g| matches!(g, Gate::t(0))));
        assert!(
            opt.gates
                .iter()
                .any(|g| matches!(g, Gate::rz(a, 0) if a.components().unwrap().1.get() == 0.3))
        );
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_plus_negative_exact_quarter_cancels() {
        // Explicit pi provenance proves this cancellation without a tolerance.
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("-pi / 4.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert!(circuits_equiv(&c, &opt, 1e-10));
        assert!(opt.gates.is_empty());
    }

    #[test]
    fn mixed_int_float_across_cnot() {
        // Route parity X (q0's initial) onto q1 via two cnots, so a rz(0.3)
        // on q1 merges with an earlier T on q0 (same fingerprint group → one
        // LivePhase with an exact quarter turn and a numeric residual).
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0)); // parity X, int=1
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        }); // q0 = X^Y, q1 = Y
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        }); // q1 = Y^(X^Y) = X
        c.apply(Gate::rz_f64(0.3, 1).unwrap()); // parity X → merge
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 2);
        match opt.gates.iter().find(|g| matches!(g, Gate::rz(..))) {
            Some(Gate::rz(theta, _)) => {
                let expected = 0.3;
                assert!(
                    (theta.to_f64_lossy().unwrap() - expected).abs() < 1e-9,
                    "theta={theta}, expected={expected}"
                );
            }
            other => panic!("expected combined rz gate, got {other:?}"),
        }
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_plus_two_rz_quarter_pi_cancel_to_z() {
        // S (2) + rz(π/4) (1) + rz(π/4) (1) → 4 = Z exactly.
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_rz_irrational_rz_opposite_irrational_leaves_t() {
        // T + rz(0.3) + rz(-0.3) → T (the residual cancels exactly).
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz_f64(0.3, 0).unwrap());
        c.apply(Gate::rz_f64(-0.3, 0).unwrap());
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_t_rz_pi_combine_to_sdg() {
        // S (2) + T (1) + rz(π) (4) → 7 = Tdg
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::t(0));
        c.apply(Gate::rz(crate::angle_expr::parse("pi", 1).unwrap(), 0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::tdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_t_rz_quarter_pi_cancel() {
        // Sdg (6) + T (1) + rz(π/4) (1) → 8 & 7 = 0 → no phase gate
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::t(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    // --- measurement / reset transfer tests ---

    #[test]
    fn measure_preserves_parity_for_t_t_merge() {
        // Measurement has transfer x' = x, so the two T gates fold into S.
        let mut c = Circuit::with_cbits(1, 1);
        c.apply(Gate::t(0));
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 0);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(
            opt.gates.as_slice(),
            [Gate::measure { qubit: 0, cbit: 0 }, Gate::s(0)]
        ));
        assert!(opt.has_measurement());
    }

    #[test]
    fn reset_sets_zero_and_eliminates_following_phase() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::reset(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 1);
        assert_eq!(opt.gates.len(), 2);
        assert!(matches!(&opt.gates[1], Gate::reset(0)));
    }

    #[test]
    fn measure_preserves_parity_for_t_tdg_cancel() {
        let mut c = Circuit::with_cbits(1, 1);
        c.apply(Gate::t(0));
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::tdg(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(matches!(
            opt.gates.as_slice(),
            [Gate::measure { qubit: 0, cbit: 0 }]
        ));
    }

    #[test]
    fn measure_on_other_qubit_does_not_block_merge() {
        // T q0; measure q1 -> c0; T q0 — same wire as the rotations is undisturbed.
        let mut c = Circuit::with_cbits(2, 1);
        c.apply(Gate::t(0));
        c.apply(Gate::measure { qubit: 1, cbit: 0 });
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        // T+T merge to S on q0.
        assert_eq!(count_t_gates(&opt), 0);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(opt.gates.iter().any(|g| matches!(g, Gate::s(0))));
        // Measurement is preserved.
        assert!(
            opt.gates
                .iter()
                .any(|g| matches!(g, Gate::measure { qubit: 1, cbit: 0 }))
        );
    }

    #[test]
    fn reset_on_other_qubit_does_not_block_merge() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::reset(1));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(opt.gates.iter().any(|g| matches!(g, Gate::s(0))));
        assert!(opt.gates.iter().any(|g| matches!(g, Gate::reset(1))));
    }

    #[test]
    fn measure_preserved_with_no_rotations() {
        let mut c = Circuit::with_cbits(2, 2);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::measure { qubit: 1, cbit: 1 });
        let opt = phase_fold_rand(&c);
        assert_eq!(opt.gates.len(), 4);
        assert_eq!(opt.num_cbits, 2);
        assert!(opt.has_measurement());
    }

    #[test]
    fn measure_keeps_one_phase_group_across_both_sides() {
        let mut c = Circuit::with_cbits(1, 1);
        c.apply(Gate::t(0));
        c.apply(Gate::t(0));
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::t(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(
            opt.gates.as_slice(),
            [Gate::measure { qubit: 0, cbit: 0 }, Gate::z(0)]
        ));
    }

    #[test]
    fn reset_eliminates_integer_and_float_phases_on_known_zero() {
        let mut c = Circuit::new(1);
        c.apply(Gate::reset(0));
        c.apply(Gate::t(0));
        c.apply(Gate::s(0));
        c.apply(Gate::z(0));
        c.apply(Gate::rz_f64(0.123, 0).unwrap());
        let opt = phase_fold_rand(&c);
        assert!(matches!(opt.gates.as_slice(), [Gate::reset(0)]));
    }

    #[test]
    fn reset_then_x_eliminates_phase_on_known_one() {
        let mut c = Circuit::new(1);
        c.apply(Gate::reset(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert!(matches!(opt.gates.as_slice(), [Gate::reset(0), Gate::x(0)]));
    }

    #[test]
    fn reset_zero_propagates_through_cnot() {
        let mut c = Circuit::new(2);
        c.apply(Gate::reset(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(0));
        c.apply(Gate::t(1));
        let opt = phase_fold_rand(&c);
        assert!(matches!(
            opt.gates.as_slice(),
            [
                Gate::reset(1),
                Gate::cnot {
                    control: 0,
                    target: 1
                },
                Gate::s(1)
            ]
        ));
    }

    #[test]
    fn hadamard_after_reset_makes_phase_nonconstant() {
        let mut c = Circuit::new(1);
        c.apply(Gate::reset(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 1);
        assert!(matches!(
            opt.gates.as_slice(),
            [Gate::reset(0), Gate::h(0), Gate::t(0)]
        ));
    }

    #[test]
    fn cz_preserved_without_rotations() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        let opt = phase_fold_rand(&c);
        assert!(matches!(
            opt.gates.as_slice(),
            [Gate::cz {
                control: 0,
                target: 1
            }]
        ));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_pair_folds_through_either_cz_operand() {
        for q in 0..2 {
            let mut c = Circuit::new(2);
            c.apply(Gate::t(q));
            c.apply(Gate::cz {
                control: 0,
                target: 1,
            });
            c.apply(Gate::t(q));
            let opt = phase_fold_rand(&c);
            assert_eq!(opt.gates.len(), 2, "operand {q}");
            assert!(
                matches!(
                    opt.gates[0],
                    Gate::cz {
                        control: 0,
                        target: 1
                    }
                ),
                "operand {q}"
            );
            assert!(matches!(opt.gates[1], Gate::s(p) if p == q), "operand {q}");
            assert!(circuits_equiv(&c, &opt, 1e-10), "operand {q}");
        }
    }

    #[test]
    fn t_pair_folds_through_every_ccz_operand() {
        for q in 0..3 {
            let mut c = Circuit::new(3);
            c.apply(Gate::t(q));
            c.apply(Gate::ccz {
                control1: 0,
                control2: 1,
                target: 2,
            });
            c.apply(Gate::t(q));

            let opt = phase_fold_rand(&c);

            assert!(matches!(opt.gates[0], Gate::ccz { .. }), "operand {q}");
            assert!(matches!(opt.gates[1], Gate::s(p) if p == q), "operand {q}");
            assert!(circuits_equiv(&c, &opt, 1e-10), "operand {q}");
        }
    }

    #[test]
    fn opposite_phases_cancel_through_cz_on_both_operands() {
        for q in 0..2 {
            let mut c = Circuit::new(2);
            c.apply(Gate::t(q));
            c.apply(Gate::cz {
                control: 0,
                target: 1,
            });
            c.apply(Gate::tdg(q));
            let opt = phase_fold_rand(&c);
            assert!(matches!(
                opt.gates.as_slice(),
                [Gate::cz {
                    control: 0,
                    target: 1
                }]
            ));
            assert!(circuits_equiv(&c, &opt, 1e-10));
        }
    }

    #[test]
    fn arbitrary_rz_folds_through_cz() {
        let mut c = Circuit::new(2);
        c.apply(Gate::rz_f64(0.37, 1).unwrap());
        c.apply(Gate::cz {
            control: 1,
            target: 0,
        });
        c.apply(Gate::rz_f64(-0.12, 1).unwrap());
        let opt = phase_fold_rand(&c);
        assert_eq!(opt.gates.len(), 2);
        assert!(matches!(
            opt.gates[0],
            Gate::cz {
                control: 1,
                target: 0
            }
        ));
        match &opt.gates[1] {
            Gate::rz(theta, 1) => assert!((theta.to_f64_lossy().unwrap() - 0.25).abs() < 1e-12),
            ref other => panic!("expected merged Rz, got {other:?}"),
        }
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cz_does_not_hide_real_hadamard_boundary() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 2);
        assert_eq!(opt.gates.len(), 4);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn phases_on_both_wires_fold_independently_across_cz_chain() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(1));
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cz {
            control: 1,
            target: 0,
        });
        c.apply(Gate::t(0));
        c.apply(Gate::t(1));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 0);
        assert!(opt.gates.iter().any(|g| matches!(g, Gate::s(0))));
        assert!(
            !opt.gates
                .iter()
                .any(|g| matches!(g, Gate::t(1) | Gate::tdg(1)))
        );
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cz_and_cnot_control_preserve_phase_parity_together() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(0));
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        c.apply(Gate::cz {
            control: 2,
            target: 1,
        });
        c.apply(Gate::t(0));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 0);
        assert!(matches!(opt.gates.last(), Some(Gate::s(0))));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cz_does_not_mask_cnot_target_parity_change() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(1));
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cz {
            control: 1,
            target: 0,
        });
        c.apply(Gate::t(1));
        let opt = phase_fold_rand(&c);
        assert_eq!(count_t_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn phase_fold_preserves_cz_count_and_operand_order() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(2));
        c.apply(Gate::cz {
            control: 2,
            target: 0,
        });
        c.apply(Gate::cz {
            control: 1,
            target: 2,
        });
        c.apply(Gate::t(2));
        let opt = phase_fold_rand(&c);
        let czs: Vec<_> = opt
            .gates
            .iter()
            .filter(|g| matches!(g, Gate::cz { .. }))
            .collect();
        assert_eq!(czs.len(), 2);
        assert!(matches!(
            czs[0],
            Gate::cz {
                control: 2,
                target: 0
            }
        ));
        assert!(matches!(
            czs[1],
            Gate::cz {
                control: 1,
                target: 2
            }
        ));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }
}
