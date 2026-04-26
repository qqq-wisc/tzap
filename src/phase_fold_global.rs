use std::f64::consts::PI;
use std::hash::{BuildHasher, Hasher};
use std::collections::hash_map::RandomState;

use indicatif::ProgressBar;
use rustc_hash::FxHashMap;

use crate::circuit::{Circuit, Gate};
use crate::pass::Pass;

/// Random tag representing a qubit's parity (XOR of variable set).
/// NOT is bitwise complement (!h), XOR is bitwise xor (a ^ b).
type ParityHash = u128;

fn fresh_parity() -> ParityHash {
    let hi = RandomState::new().build_hasher().finish() as u128;
    let lo = RandomState::new().build_hasher().finish() as u128;
    (hi << 64) | lo
}

/// Accumulated phase for a parity group.
/// `int_part` counts π/4 steps mod 8 (Clifford+T gates land here exclusively).
/// `float_part` holds leftover rotation for rz gates whose angle is not a π/4 multiple.
/// Pure Clifford+T circuits only ever touch `int_part`.
struct LivePhase {
    int_part: u8,
    qubit: usize,
    current_idx: usize,
    float_part: f64,
    /// True when `current_idx` sits on the complement of the group's canonical parity.
    /// On emission, the accumulated rotation is negated to compensate.
    current_sign: bool,
}

pub struct PhaseFoldGlobal;

impl Pass for PhaseFoldGlobal {
    fn name(&self) -> &str { "Phase folding" }
    fn run_with_progress(&self, circuit: &Circuit, pb: &ProgressBar) -> Circuit {
        phase_fold_global(circuit, pb)
    }
}

pub fn phase_fold_global(circuit: &Circuit, pb: &ProgressBar) -> Circuit {
    let n = circuit.num_qubits;
    let fresh = || fresh_parity();
    let mut qubits: Vec<ParityHash> = (0..n).map(|_| fresh()).collect();

    let mut live: Vec<LivePhase> = Vec::new();
    let mut parity_to_group: FxHashMap<ParityHash, usize> = FxHashMap::default();
    let mut skip = vec![false; circuit.gates.len()];
    let mut emit_at: Vec<Option<(usize, u8, f64)>> = vec![None; circuit.gates.len()];

    for (idx, gate) in circuit.gates.iter().enumerate() {
        if idx & 0xFFF == 0 { pb.inc(0x1000); }
        match gate {
            Gate::t(q)   => record_int(&qubits, *q, 1, idx, &mut live, &mut parity_to_group, &mut skip),
            Gate::tdg(q) => record_int(&qubits, *q, 7, idx, &mut live, &mut parity_to_group, &mut skip),
            Gate::s(q)   => record_int(&qubits, *q, 2, idx, &mut live, &mut parity_to_group, &mut skip),
            Gate::sdg(q) => record_int(&qubits, *q, 6, idx, &mut live, &mut parity_to_group, &mut skip),
            Gate::z(q)   => record_int(&qubits, *q, 4, idx, &mut live, &mut parity_to_group, &mut skip),
            Gate::rz(theta, q) => match classify_quarter_pi(*theta) {
                Some(k) => record_int(&qubits, *q, k, idx, &mut live, &mut parity_to_group, &mut skip),
                None    => record_float(&qubits, *q, *theta, idx, &mut live, &mut parity_to_group, &mut skip),
            },
            Gate::h(q) => { qubits[*q] = fresh(); }
            Gate::x(q) => { qubits[*q] = !qubits[*q]; }
            Gate::cnot { control, target } => {
                qubits[*target] ^= qubits[*control];
            }
            Gate::ccx { control1: _, control2: _, target } => {
                // AND-based parity — opaque, just refresh the target.
                qubits[*target] = fresh();
            }
        }
    }

    // Finalize: emit surviving groups or skip if angle cancelled to zero.
    for lp in &live {
        let (int_part, float_part) = if lp.current_sign {
            // Last seen location is on ¬canonical — negate to compensate.
            (8u8.wrapping_sub(lp.int_part) & 7, -lp.float_part)
        } else {
            (lp.int_part, lp.float_part)
        };

        if float_part == 0.0 {
            // Pure integer group — no float math.
            if int_part == 0 {
                skip[lp.current_idx] = true;
            } else {
                emit_at[lp.current_idx] = Some((lp.qubit, int_part, 0.0));
            }
        } else {
            let total = (int_part as f64) * (PI / 4.0) + float_part;
            if angle_is_zero(total) {
                skip[lp.current_idx] = true;
            } else {
                emit_at[lp.current_idx] = Some((lp.qubit, 0, total));
            }
        }
    }

    // Reconstruct circuit, moving gates to avoid per-gate clones.
    let mut output = Circuit::new(n);
    for (idx, gate) in circuit.gates.clone().into_iter().enumerate() {
        if skip[idx] {
            continue;
        }
        if let Some((qubit, int_part, float_part)) = emit_at[idx] {
            emit_rotation(&mut output, qubit, int_part, float_part);
        } else {
            output.apply(gate);
        }
    }

    output
}

fn record_int(
    qubits: &[ParityHash],
    q: usize,
    k: u8,
    idx: usize,
    live: &mut Vec<LivePhase>,
    parity_to_group: &mut FxHashMap<ParityHash, usize>,
    skip: &mut [bool],
) {
    let parity = qubits[q];

    // Direct match.
    if let Some(&gi) = parity_to_group.get(&parity) {
        skip[live[gi].current_idx] = true;
        live[gi].int_part = (live[gi].int_part + k) & 7;
        live[gi].current_idx = idx;
        live[gi].qubit = q;
        live[gi].current_sign = false;
        return;
    }

    // Complement match: RZ(θ) on ¬p ≡ RZ(−θ) on p up to global phase.
    if let Some(&gi) = parity_to_group.get(&!parity) {
        skip[live[gi].current_idx] = true;
        live[gi].int_part = live[gi].int_part.wrapping_sub(k) & 7;
        live[gi].current_idx = idx;
        live[gi].qubit = q;
        live[gi].current_sign = true;
        return;
    }

    let gi = live.len();
    live.push(LivePhase {
        int_part: k,
        float_part: 0.0,
        current_idx: idx,
        qubit: q,
        current_sign: false,
    });
    parity_to_group.insert(parity, gi);
}

fn record_float(
    qubits: &[ParityHash],
    q: usize,
    theta: f64,
    idx: usize,
    live: &mut Vec<LivePhase>,
    parity_to_group: &mut FxHashMap<ParityHash, usize>,
    skip: &mut [bool],
) {
    let parity = qubits[q];

    if let Some(&gi) = parity_to_group.get(&parity) {
        skip[live[gi].current_idx] = true;
        live[gi].float_part += theta;
        live[gi].current_idx = idx;
        live[gi].qubit = q;
        live[gi].current_sign = false;
        return;
    }

    if let Some(&gi) = parity_to_group.get(&!parity) {
        skip[live[gi].current_idx] = true;
        live[gi].float_part -= theta;
        live[gi].current_idx = idx;
        live[gi].qubit = q;
        live[gi].current_sign = true;
        return;
    }

    let gi = live.len();
    live.push(LivePhase {
        int_part: 0,
        float_part: theta,
        current_idx: idx,
        qubit: q,
        current_sign: false,
    });
    parity_to_group.insert(parity, gi);
}

/// Returns Some(k) if theta ≈ k · π/4 (mod 2π) within 1e-9, else None.
fn classify_quarter_pi(theta: f64) -> Option<u8> {
    let n = theta.rem_euclid(2.0 * PI);
    let q = PI / 4.0;
    let k = (n / q).round();
    if (n - k * q).abs() < 1e-9 {
        Some((k as u32 & 7) as u8)
    } else {
        None
    }
}

fn angle_is_zero(angle: f64) -> bool {
    let n = angle.rem_euclid(2.0 * PI);
    n < 1e-6 || (2.0 * PI - n) < 1e-6
}

fn emit_rotation(output: &mut Circuit, q: usize, int_part: u8, float_part: f64) {
    if float_part == 0.0 {
        emit_int(output, q, int_part);
    } else {
        let theta = (int_part as f64) * (PI / 4.0) + float_part;
        emit_float(output, q, theta);
    }
}

fn emit_int(output: &mut Circuit, q: usize, k: u8) {
    match k & 7 {
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
}

fn emit_float(output: &mut Circuit, q: usize, angle: f64) {
    let n = angle.rem_euclid(2.0 * PI);
    if angle_is_zero(n) {
        return;
    }
    let quarter = PI / 4.0;
    let k = (n / quarter).round();
    // Use 1e-6 tolerance to handle floating-point drift from summing many π/4 multiples.
    if (n - k * quarter).abs() < 1e-6 {
        emit_int(output, q, (k as u32 & 7) as u8);
    } else {
        output.apply(Gate::rz(n, q));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decompose::DecomposeToffoli;
    use crate::pass::Pass;
    use crate::unitary::circuits_equiv;

    fn count_gates(c: &Circuit) -> usize {
        c.gates.len()
    }

    fn count_phase_gates(c: &Circuit) -> usize {
        c.gates.iter().filter(|g| matches!(g,
            Gate::t(_) | Gate::tdg(_) | Gate::s(_) | Gate::sdg(_) | Gate::z(_) | Gate::rz(..)
        )).count()
    }

    fn count_t_gates(c: &Circuit) -> usize {
        c.gates.iter().filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count()
    }

    #[test]
    #[ignore] // long-running: 24-qubit pipeline benchmark
    fn large_circuit_toffoli_decompose_and_phase_fold() {
        let mut c = Circuit::new(24);
        c.apply(Gate::cnot { control: 3, target: 2 });
        c.apply(Gate::cnot { control: 8, target: 7 });
        c.apply(Gate::cnot { control: 14, target: 13 });
        c.apply(Gate::cnot { control: 21, target: 20 });
        c.apply(Gate::cnot { control: 3, target: 4 });
        c.apply(Gate::cnot { control: 8, target: 9 });
        c.apply(Gate::cnot { control: 14, target: 15 });
        c.apply(Gate::cnot { control: 21, target: 22 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 3 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx { control1: 5, control2: 6, target: 8 });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx { control1: 11, control2: 12, target: 14 });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx { control1: 18, control2: 19, target: 21 });
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
        c.apply(Gate::ccx { control1: 2, control2: 3, target: 4 });
        c.apply(Gate::h(4));
        c.apply(Gate::h(9));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 9 });
        c.apply(Gate::h(9));
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 15 });
        c.apply(Gate::h(15));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::h(22));
        c.apply(Gate::ccx { control1: 20, control2: 21, target: 22 });
        c.apply(Gate::h(22));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx { control1: 20, control2: 21, target: 23 });
        c.apply(Gate::h(23));
        c.apply(Gate::cnot { control: 6, target: 5 });
        c.apply(Gate::cnot { control: 12, target: 11 });
        c.apply(Gate::cnot { control: 19, target: 18 });
        c.apply(Gate::cnot { control: 5, target: 8 });
        c.apply(Gate::cnot { control: 11, target: 14 });
        c.apply(Gate::cnot { control: 18, target: 21 });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx { control1: 20, control2: 21, target: 23 });
        c.apply(Gate::h(23));
        c.apply(Gate::h(4));
        c.apply(Gate::h(10));
        c.apply(Gate::h(15));
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::h(17));
        c.apply(Gate::h(17));
        c.apply(Gate::ccx { control1: 16, control2: 23, target: 17 });
        c.apply(Gate::h(17));
        c.apply(Gate::h(22));
        c.apply(Gate::ccx { control1: 15, control2: 23, target: 22 });
        c.apply(Gate::h(22));
        c.apply(Gate::h(9));
        c.apply(Gate::ccx { control1: 4, control2: 10, target: 9 });
        c.apply(Gate::h(9));
        c.apply(Gate::h(17));
        c.apply(Gate::h(9));
        c.apply(Gate::h(15));
        c.apply(Gate::h(22));
        c.apply(Gate::ccx { control1: 9, control2: 17, target: 22 });
        c.apply(Gate::h(22));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx { control1: 9, control2: 16, target: 15 });
        c.apply(Gate::h(15));
        c.apply(Gate::h(15));
        c.apply(Gate::h(22));
        c.apply(Gate::h(17));
        c.apply(Gate::h(17));
        c.apply(Gate::ccx { control1: 16, control2: 23, target: 17 });
        c.apply(Gate::h(17));
        c.apply(Gate::h(17));
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx { control1: 20, control2: 21, target: 23 });
        c.apply(Gate::h(23));
        c.apply(Gate::cnot { control: 5, target: 8 });
        c.apply(Gate::cnot { control: 11, target: 14 });
        c.apply(Gate::cnot { control: 18, target: 21 });
        c.apply(Gate::cnot { control: 6, target: 5 });
        c.apply(Gate::cnot { control: 12, target: 11 });
        c.apply(Gate::cnot { control: 19, target: 18 });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::h(23));
        c.apply(Gate::ccx { control1: 20, control2: 21, target: 23 });
        c.apply(Gate::h(23));
        c.apply(Gate::h(23));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 3 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx { control1: 5, control2: 6, target: 8 });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx { control1: 11, control2: 12, target: 14 });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx { control1: 18, control2: 19, target: 21 });
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::cnot { control: 3, target: 2 });
        c.apply(Gate::cnot { control: 8, target: 7 });
        c.apply(Gate::cnot { control: 14, target: 13 });
        c.apply(Gate::cnot { control: 21, target: 20 });
        c.apply(Gate::cnot { control: 6, target: 5 });
        c.apply(Gate::cnot { control: 12, target: 11 });
        c.apply(Gate::cnot { control: 19, target: 18 });
        c.apply(Gate::cnot { control: 6, target: 8 });
        c.apply(Gate::cnot { control: 12, target: 14 });
        c.apply(Gate::cnot { control: 19, target: 21 });
        c.apply(Gate::cnot { control: 4, target: 6 });
        c.apply(Gate::cnot { control: 9, target: 12 });
        c.apply(Gate::cnot { control: 15, target: 19 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 3 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx { control1: 5, control2: 6, target: 8 });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx { control1: 11, control2: 12, target: 14 });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx { control1: 18, control2: 19, target: 21 });
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::cnot { control: 3, target: 2 });
        c.apply(Gate::cnot { control: 8, target: 7 });
        c.apply(Gate::cnot { control: 14, target: 13 });
        c.apply(Gate::cnot { control: 21, target: 20 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 3 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx { control1: 5, control2: 6, target: 8 });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx { control1: 11, control2: 12, target: 14 });
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::ccx { control1: 18, control2: 19, target: 21 });
        c.apply(Gate::h(21));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(21));
        c.apply(Gate::cnot { control: 6, target: 5 });
        c.apply(Gate::cnot { control: 12, target: 11 });
        c.apply(Gate::cnot { control: 19, target: 18 });
        c.apply(Gate::cnot { control: 4, target: 6 });
        c.apply(Gate::cnot { control: 9, target: 12 });
        c.apply(Gate::cnot { control: 15, target: 19 });
        c.apply(Gate::cnot { control: 6, target: 8 });
        c.apply(Gate::cnot { control: 12, target: 14 });
        c.apply(Gate::cnot { control: 19, target: 21 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::cnot { control: 6, target: 5 });
        c.apply(Gate::cnot { control: 12, target: 11 });
        c.apply(Gate::cnot { control: 19, target: 18 });
        c.apply(Gate::x(0));
        c.apply(Gate::x(2));
        c.apply(Gate::x(5));
        c.apply(Gate::x(7));
        c.apply(Gate::x(11));
        c.apply(Gate::x(13));
        c.apply(Gate::cnot { control: 3, target: 2 });
        c.apply(Gate::cnot { control: 8, target: 7 });
        c.apply(Gate::cnot { control: 14, target: 13 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::h(3));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 3 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx { control1: 5, control2: 6, target: 8 });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx { control1: 11, control2: 12, target: 14 });
        c.apply(Gate::h(14));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::cnot { control: 6, target: 5 });
        c.apply(Gate::cnot { control: 12, target: 11 });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::cnot { control: 5, target: 8 });
        c.apply(Gate::cnot { control: 11, target: 14 });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::h(15));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx { control1: 9, control2: 16, target: 15 });
        c.apply(Gate::h(15));
        c.apply(Gate::h(9));
        c.apply(Gate::h(9));
        c.apply(Gate::ccx { control1: 4, control2: 10, target: 9 });
        c.apply(Gate::h(9));
        c.apply(Gate::h(4));
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::cnot { control: 5, target: 8 });
        c.apply(Gate::cnot { control: 11, target: 14 });
        c.apply(Gate::h(10));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 10 });
        c.apply(Gate::h(10));
        c.apply(Gate::h(16));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 16 });
        c.apply(Gate::h(16));
        c.apply(Gate::cnot { control: 6, target: 5 });
        c.apply(Gate::cnot { control: 12, target: 11 });
        c.apply(Gate::h(9));
        c.apply(Gate::ccx { control1: 7, control2: 8, target: 9 });
        c.apply(Gate::h(9));
        c.apply(Gate::h(15));
        c.apply(Gate::ccx { control1: 13, control2: 14, target: 15 });
        c.apply(Gate::h(15));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx { control1: 2, control2: 3, target: 4 });
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
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 3 });
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::ccx { control1: 5, control2: 6, target: 8 });
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::ccx { control1: 11, control2: 12, target: 14 });
        c.apply(Gate::h(14));
        c.apply(Gate::h(3));
        c.apply(Gate::h(8));
        c.apply(Gate::h(14));
        c.apply(Gate::cnot { control: 3, target: 4 });
        c.apply(Gate::cnot { control: 8, target: 9 });
        c.apply(Gate::cnot { control: 14, target: 15 });
        c.apply(Gate::cnot { control: 3, target: 2 });
        c.apply(Gate::cnot { control: 8, target: 7 });
        c.apply(Gate::cnot { control: 14, target: 13 });
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

        let folded = phase_fold_global(&dec, &ProgressBar::hidden());
        let folded_t = count_t_gates(&folded);

        let folded2 = phase_fold_global(&folded, &ProgressBar::hidden());
        let folded2_t = count_t_gates(&folded2);

        println!("=== large circuit pipeline ===");
        println!("original:            {} gates, {} T/Tdg", original_gates, original_t);
        println!("after toffoli decomp: {} gates, {} T/Tdg", count_gates(&dec), dec_t);
        println!("after phase fold 1:  {} gates, {} T/Tdg", count_gates(&folded), folded_t);
        println!("after phase fold 2:  {} gates, {} T/Tdg", count_gates(&folded2), folded2_t);

        assert!(folded_t <= dec_t, "phase folding should not increase T count");
        assert!(folded2_t <= folded_t, "second phase fold should not increase T count");
    }

    #[test]
    #[ignore] // long-running: 5-qubit multi-pass pipeline
    fn small_circuit_pipeline() {
        let mut c = Circuit::new(5);
        c.apply(Gate::x(4));
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx { control1: 0, control2: 3, target: 4 });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx { control1: 2, control2: 3, target: 4 });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::cnot { control: 3, target: 4 });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx { control1: 1, control2: 2, target: 4 });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::cnot { control: 2, target: 4 });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 4 });
        c.apply(Gate::h(4));
        c.apply(Gate::h(4));
        c.apply(Gate::cnot { control: 1, target: 4 });
        c.apply(Gate::cnot { control: 0, target: 4 });

        let original_gates = count_gates(&c);
        let original_t = count_t_gates(&c);

        let dec = DecomposeToffoli.run(&c);
        let dec_t = count_t_gates(&dec);

        let folded = phase_fold_global(&dec, &ProgressBar::hidden());
        let folded_t = count_t_gates(&folded);

        let folded2 = phase_fold_global(&folded, &ProgressBar::hidden());
        let folded2_t = count_t_gates(&folded2);

        println!("=== small circuit pipeline ===");
        println!("original:            {} gates, {} T/Tdg", original_gates, original_t);
        println!("after toffoli decomp: {} gates, {} T/Tdg", count_gates(&dec), dec_t);
        println!("after phase fold 1:  {} gates, {} T/Tdg", count_gates(&folded), folded_t);
        println!("after phase fold 2:  {} gates, {} T/Tdg", count_gates(&folded2), folded2_t);

        assert!(folded_t <= dec_t);
    }

    #[test]
    fn two_toffoli_shared_control() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx { control1: 1, control2: 2, target: 0 });
        c.apply(Gate::cnot { control: 2, target: 0 });
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });

        let dec = DecomposeToffoli.run(&c);
        assert_eq!(count_t_gates(&dec), 14);

        let opt = phase_fold_global(&dec, &ProgressBar::hidden());
        assert!(circuits_equiv(&c, &opt, 1e-10));
        assert_eq!(count_t_gates(&opt), 12);
    }

    #[test]
    #[ignore] // long-running: combinatorial CCX removal search
    fn mod5_4_remove_ccx_combos() {
        let ccx_indices: Vec<usize> = vec![1, 2, 4, 6];
        let names = ["ccx q0,q3,q4", "ccx q2,q3,q4", "ccx q1,q2,q4", "ccx q0,q1,q4"];

        let build = |skip: &[usize]| -> Circuit {
            let mut c = Circuit::new(5);
            let gates: Vec<Gate> = vec![
                Gate::x(4),
                Gate::ccx { control1: 0, control2: 3, target: 4 },
                Gate::ccx { control1: 2, control2: 3, target: 4 },
                Gate::cnot { control: 3, target: 4 },
                Gate::ccx { control1: 1, control2: 2, target: 4 },
                Gate::cnot { control: 2, target: 4 },
                Gate::ccx { control1: 0, control2: 1, target: 4 },
                Gate::cnot { control: 1, target: 4 },
                Gate::cnot { control: 0, target: 4 },
            ];
            for (j, g) in gates.into_iter().enumerate() {
                if !skip.contains(&j) { c.apply(g); }
            }
            c
        };

        println!("\n{:<50} {:>5} {:>5}", "removed", "dec_T", "ours");
        println!("{}", "-".repeat(65));

        for keep in 0..4 {
            let skip: Vec<usize> = (0..4).filter(|&i| i != keep).map(|i| ccx_indices[i]).collect();
            let c = build(&skip);
            let dec = DecomposeToffoli.run(&c);
            let dec_t = count_t_gates(&dec);
            let opt = phase_fold_global(&dec, &ProgressBar::hidden());
            let opt_t = count_t_gates(&opt);
            assert!(circuits_equiv(&dec, &opt, 1e-10));
            let removed: Vec<&str> = (0..4).filter(|&i| i != keep).map(|i| names[i]).collect();
            println!("{:<50} {:>5} {:>5}", removed.join(" + "), dec_t, opt_t);
        }
    }

    #[test]
    fn toffoli_decompose_then_phase_fold() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        let dec = DecomposeToffoli.run(&c);
        let opt = phase_fold_global(&dec, &ProgressBar::hidden());
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
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_and_tdg_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn four_t_merge_to_s_s() {
        let mut c = Circuit::new(1);
        for _ in 0..4 { c.apply(Gate::t(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn eight_t_cancel() {
        let mut c = Circuit::new(1);
        for _ in 0..8 { c.apply(Gate::t(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn same_parity_across_cnot() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn different_parity_no_merge() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn h_prevents_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn merge_across_x() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        // T; X; Tdg ≡ Sdg·X up to global phase — one phase gate survives.
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_folds_across_x() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(0.3, 0));
        c.apply(Gate::x(0));
        c.apply(Gate::rz(0.7, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_cancels_across_x() {
        // RZ(θ); X; RZ(θ) = RZ(θ)·RZ(−θ)·X = X up to global phase.
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(0.42, 0));
        c.apply(Gate::x(0));
        c.apply(Gate::rz(0.42, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        c.apply(Gate::rz(0.5, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
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
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::x(1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::x(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        c.apply(Gate::t(0));       // group[p0] int=1, sign=F
        c.apply(Gate::x(0));       // qubits[0] = !p0
        c.apply(Gate::t(0));       // complement hit → int=0, sign=T
        c.apply(Gate::x(0));       // qubits[0] = p0
        c.apply(Gate::tdg(0));     // direct hit on p0 → int=(0-1)&7=7, sign=F
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_t_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_qubit_folding() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn toffoli_decomposition_fold() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(2));
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::tdg(2));
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::t(2));
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::tdg(2));
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::t(1));
        c.apply(Gate::t(2));
        c.apply(Gate::h(2));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(1));
        c.apply(Gate::cnot { control: 0, target: 1 });

        let original_phases = count_phase_gates(&c);
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        let optimized_phases = count_phase_gates(&opt);

        println!("toffoli: {original_phases} phase gates -> {optimized_phases} phase gates");
        assert!(optimized_phases <= original_phases);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn preserves_non_phase_structure() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::x(1));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 3);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(0.3, 0));
        c.apply(Gate::rz(0.7, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_merge_via_cnot() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::t(0));

        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_cancel() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::tdg(1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::t(0));

        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_three_way() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::rz(0.5, 2));
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::rz(0.5, 2));

        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_rz_different_qubits_same_parity() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::rz(0.3, 1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::rz(0.7, 0));

        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    fn print_before_after(name: &str, original: &Circuit, optimized: &Circuit) {
        println!("=== {name} ===");
        println!("BEFORE ({} gates, {} phase):", count_gates(original), count_phase_gates(original));
        print!("{original}");
        println!("AFTER ({} gates, {} phase):", count_gates(optimized), count_phase_gates(optimized));
        print!("{optimized}");
        println!();
    }

    #[test]
    fn circuit_from_diagram() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::t(2));
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::tdg(0));
        c.apply(Gate::cnot { control: 2, target: 0 });

        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        print_before_after("diagram circuit (T on y, T† on x)", &c, &opt);

        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cx_t_cx_cx_tdg_cx() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::tdg(0));
        c.apply(Gate::cnot { control: 0, target: 1 });

        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        print_before_after("cx T cx cx Tdg cx", &c, &opt);

        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_swap_h_swap_t() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::h(1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(1));

        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        print_before_after("T, SWAP, H, SWAP, T", &c, &opt);

        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(opt.gates.last().unwrap(), Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn demo_cross_qubit_all() {
        let mut c1 = Circuit::new(2);
        c1.apply(Gate::t(0));
        c1.apply(Gate::cnot { control: 0, target: 1 });
        c1.apply(Gate::t(1));
        c1.apply(Gate::cnot { control: 0, target: 1 });
        c1.apply(Gate::cnot { control: 1, target: 0 });
        c1.apply(Gate::t(0));
        let opt1 = phase_fold_global(&c1, &ProgressBar::hidden());
        print_before_after("cross-qubit merge (T on q0 + T on q1 -> S)", &c1, &opt1);

        let mut c2 = Circuit::new(2);
        c2.apply(Gate::t(0));
        c2.apply(Gate::cnot { control: 0, target: 1 });
        c2.apply(Gate::tdg(1));
        c2.apply(Gate::cnot { control: 0, target: 1 });
        c2.apply(Gate::cnot { control: 1, target: 0 });
        c2.apply(Gate::t(0));
        let opt2 = phase_fold_global(&c2, &ProgressBar::hidden());
        print_before_after("cross-qubit cancel (T on q0 + Tdg on q1 -> 0)", &c2, &opt2);

        let mut c3 = Circuit::new(3);
        c3.apply(Gate::cnot { control: 0, target: 2 });
        c3.apply(Gate::cnot { control: 1, target: 2 });
        c3.apply(Gate::rz(0.5, 2));
        c3.apply(Gate::cnot { control: 1, target: 2 });
        c3.apply(Gate::cnot { control: 0, target: 2 });
        c3.apply(Gate::cnot { control: 1, target: 2 });
        c3.apply(Gate::cnot { control: 0, target: 2 });
        c3.apply(Gate::rz(0.5, 2));
        let opt3 = phase_fold_global(&c3, &ProgressBar::hidden());
        print_before_after("3-qubit merge (Rz on q2 twice, same parity)", &c3, &opt3);

        let mut c4 = Circuit::new(2);
        c4.apply(Gate::cnot { control: 0, target: 1 });
        c4.apply(Gate::rz(0.3, 1));
        c4.apply(Gate::cnot { control: 0, target: 1 });
        c4.apply(Gate::cnot { control: 1, target: 0 });
        c4.apply(Gate::rz(0.7, 0));
        let opt4 = phase_fold_global(&c4, &ProgressBar::hidden());
        print_before_after("Rz(0.3) on q1 + Rz(0.7) on q0 -> Rz(1.0)", &c4, &opt4);
    }

    // --- z and sdg gate tests ---

    #[test]
    fn z_is_phase_gate() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_is_phase_gate() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_z_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_sdg_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_s_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::s(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_s_is_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::s(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_sdg_is_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_t_t_t_is_z() {
        let mut c = Circuit::new(1);
        for _ in 0..4 { c.apply(Gate::t(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_tdg_is_z_plus_t() {
        let mut c = Circuit::new(1);
        for _ in 0..3 { c.apply(Gate::tdg(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_t_is_z_plus_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_tdg_is_s_plus_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2); // S + T (3π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_t_is_tdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::tdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_t_folds() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::t(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2); // S + T (3π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_tdg_is_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::t(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_tdg_is_z_plus_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2); // Z + T (5π/4)
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_s_is_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::s(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_sdg_is_s() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn six_t_is_sdg() {
        let mut c = Circuit::new(1);
        for _ in 0..6 { c.apply(Gate::t(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn seven_t_is_tdg() {
        let mut c = Circuit::new(1);
        for _ in 0..7 { c.apply(Gate::t(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::tdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_across_cnot_folds() {
        let mut c = Circuit::new(2);
        c.apply(Gate::z(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::z(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_across_cnot_folds() {
        let mut c = Circuit::new(2);
        c.apply(Gate::sdg(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::s(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_h_prevents_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::h(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_h_prevents_merge() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::h(0));
        c.apply(Gate::sdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_x_z_folds_across_x() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::x(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        // Z·X·Z = −X — both Zs fold into global phase.
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_z_merge() {
        let mut c = Circuit::new(2);
        c.apply(Gate::z(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn cross_qubit_sdg_s_cancel() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::sdg(1));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::s(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_in_toffoli_decomp_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::z(2));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::z(2));
        let dec = DecomposeToffoli.run(&c);
        let pf1 = phase_fold_global(&dec, &ProgressBar::hidden());
        let result = phase_fold_global(&pf1, &ProgressBar::hidden());
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    #[test]
    fn sdg_in_toffoli_decomp_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::sdg(0));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::s(1));
        let dec = DecomposeToffoli.run(&c);
        let pf1 = phase_fold_global(&dec, &ProgressBar::hidden());
        let result = phase_fold_global(&pf1, &ProgressBar::hidden());
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    #[test]
    fn multi_qubit_z_sdg_fold() {
        let mut c = Circuit::new(2);
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(1));
        c.apply(Gate::s(0));
        c.apply(Gate::t(1));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_s_is_sdg() {
        let mut c = Circuit::new(1);
        for _ in 0..3 { c.apply(Gate::s(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn three_sdg_is_s() {
        let mut c = Circuit::new(1);
        for _ in 0..3 { c.apply(Gate::sdg(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn four_s_cancel() {
        let mut c = Circuit::new(1);
        for _ in 0..4 { c.apply(Gate::s(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn four_sdg_cancel() {
        let mut c = Circuit::new(1);
        for _ in 0..4 { c.apply(Gate::sdg(0)); }
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_t_tdg_is_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn mixed_z_sdg_cnot_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::z(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::sdg(1));
        c.apply(Gate::cnot { control: 1, target: 2 });
        c.apply(Gate::t(2));
        c.apply(Gate::cnot { control: 0, target: 2 });
        c.apply(Gate::s(0));
        c.apply(Gate::tdg(2));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert!(circuits_equiv(&c, &opt, 1e-10));
        assert!(count_phase_gates(&opt) <= count_phase_gates(&c));
    }

    #[test]
    fn s_z_is_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn sdg_z_is_s() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        c.apply(Gate::z(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_pi_folds_to_z() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_neg_half_pi_folds_to_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(-PI / 2.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_three_half_pi_folds_to_sdg() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(3.0 * PI / 2.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::sdg(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn z_preserves_non_phase_structure() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::z(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::sdg(1));
        c.apply(Gate::x(1));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert!(circuits_equiv(&c, &opt, 1e-10));
        let h_count = opt.gates.iter().filter(|g| matches!(g, Gate::h(_))).count();
        let x_count = opt.gates.iter().filter(|g| matches!(g, Gate::x(_))).count();
        let cx_count = opt.gates.iter().filter(|g| matches!(g, Gate::cnot { .. })).count();
        assert_eq!(h_count, 1);
        assert_eq!(x_count, 1);
        assert_eq!(cx_count, 1);
    }

    #[test]
    fn ccx_z_sdg_full_pipeline() {
        let mut c = Circuit::new(3);
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(1));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::z(2));
        c.apply(Gate::s(1));
        let dec = DecomposeToffoli.run(&c);
        let pf1 = phase_fold_global(&dec, &ProgressBar::hidden());
        let result = phase_fold_global(&pf1, &ProgressBar::hidden());
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    // --- mixed int + float (π/4-multiple) accumulation tests ---
    // These exercise the interaction between int_part (T/Tdg/S/Sdg/Z) and
    // rz(·) gates, both when the rz angle is a π/4 multiple (int_part path)
    // and when it is arbitrary (float_part path).

    #[test]
    fn t_plus_rz_quarter_pi_is_s() {
        // T (π/4) + rz(π/4) → 2·π/4 = S
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz(PI / 4.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::s(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_plus_rz_half_pi_is_z() {
        // S (π/2) + rz(π/2) → π = Z
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::rz(PI / 2.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn rz_pi_plus_tdg_is_z_plus_t() {
        // rz(π) (4) + Tdg (7) → 11 & 7 = 3 → S + T
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI, 0));
        c.apply(Gate::tdg(0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 2); // S + T
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_tdg_rz_quarter_pi_is_t() {
        // T + Tdg cancels, rz(π/4) → T
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(0));
        c.apply(Gate::rz(PI / 4.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::t(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_plus_rz_irrational_folds() {
        // T (π/4) + rz(0.3) — float_part nonzero, should emit an rz with
        // the combined angle π/4 + 0.3 (not a π/4 multiple).
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz(0.3, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        match &opt.gates[0] {
            Gate::rz(theta, _) => {
                let expected = (PI / 4.0 + 0.3).rem_euclid(2.0 * PI);
                assert!((theta - expected).abs() < 1e-9);
            }
            g => panic!("expected rz, got {g:?}"),
        }
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_rz_irrational_cancels_to_zero_after_second_rz() {
        // T + rz(-π/4) — float_part(-π/4) doesn't match π/4-multiple at
        // classify time (0.7853... may differ from negated). Check that
        // the integer and float parts combine correctly to cancel.
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz(-PI / 4.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        // Either cancels entirely (both routed to int_part) or emits a
        // near-identity rz; equivalence check is authoritative.
        assert!(circuits_equiv(&c, &opt, 1e-10));
        assert!(count_phase_gates(&opt) <= 1);
    }

    #[test]
    fn mixed_int_float_across_cnot() {
        // Route parity X (q0's initial) onto q1 via two cnots, so a rz(0.3)
        // on q1 merges with an earlier T on q0 (same parity group → one
        // LivePhase with int_part=1, float_part=0.3).
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));                                // parity X, int=1
        c.apply(Gate::cnot { control: 1, target: 0 });      // q0 = X^Y, q1 = Y
        c.apply(Gate::cnot { control: 0, target: 1 });      // q1 = Y^(X^Y) = X
        c.apply(Gate::rz(0.3, 1));                          // parity X → merge
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        match opt.gates.iter().find(|g| matches!(g, Gate::rz(..))) {
            Some(Gate::rz(theta, _)) => {
                let expected = (PI / 4.0 + 0.3).rem_euclid(2.0 * PI);
                assert!((theta - expected).abs() < 1e-9, "theta={theta}, expected={expected}");
            }
            other => panic!("expected combined rz gate, got {other:?}"),
        }
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_plus_two_rz_quarter_pi_cancel_to_z() {
        // S (2) + rz(π/4) (1) + rz(π/4) (1) → 4 = Z (all in int_part)
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::rz(PI / 4.0, 0));
        c.apply(Gate::rz(PI / 4.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(matches!(&opt.gates[0], Gate::z(_)));
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn t_rz_irrational_rz_opposite_irrational_leaves_t() {
        // T + rz(0.3) + rz(-0.3) → T (float_part cancels, int_part stays)
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::rz(0.3, 0));
        c.apply(Gate::rz(-0.3, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 1);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }

    #[test]
    fn s_t_rz_pi_combine_to_sdg() {
        // S (2) + T (1) + rz(π) (4) → 7 = Tdg
        let mut c = Circuit::new(1);
        c.apply(Gate::s(0));
        c.apply(Gate::t(0));
        c.apply(Gate::rz(PI, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
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
        c.apply(Gate::rz(PI / 4.0, 0));
        let opt = phase_fold_global(&c, &ProgressBar::hidden());
        assert_eq!(count_phase_gates(&opt), 0);
        assert!(circuits_equiv(&c, &opt, 1e-10));
    }
}
