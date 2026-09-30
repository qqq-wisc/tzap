//! Circuit-level Pauli folding with exact commutation checks.
//!
//! A randomized linear image of each Clifford-frame row finds possible
//! repeated rotation axes. Axis equality and intervening commutation are
//! checked exactly before any rewrite is made.

use std::collections::hash_map::RandomState;
use std::f64::consts::PI;
use std::hash::{BuildHasher, Hasher};

use rustc_hash::{FxHashMap, FxHashSet};

use crate::circuit::{Circuit, Gate, qubits_of};
use crate::pass::Pass;
use crate::pbc::{packed_anticommutes, support_signature};

/// Maximum retained exact axis storage, in 64-bit words (256 MiB). Reaching
/// it stops further folding; folds already found are kept.
const MAX_AXIS_WORDS: usize = 1 << 25;
/// Work allowed for one fold attempt (candidate lookup plus the blocker
/// check). An attempt that runs out is treated as blocked, and the scan goes
/// on, so one expensive candidate cannot end the pass.
const MAX_ATTEMPT_STEPS: usize = 1 << 16;
/// Folds T and Rz rotations across Clifford gates, including Hadamards.
/// Can run independently or after [`crate::phase_fold_rand::PhaseFoldRand`].
pub struct PauliFoldRand;

impl Pass for PauliFoldRand {
    fn name(&self) -> &str {
        "Pauli folding"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        pauli_fold_rand(circuit)
    }
}

/// Circuit rewrite up to global phase. Random fingerprints only select
/// candidates; exact axis and commutation checks authorize every rewrite.
pub fn pauli_fold_rand(circuit: &Circuit) -> Circuit {
    let labels: Vec<_> = (0..circuit.num_qubits)
        .map(|_| (fresh_label(), fresh_label()))
        .collect();
    pauli_fold_with_labels(circuit, &labels)
}

fn fresh_label() -> u128 {
    let hi = RandomState::new().build_hasher().finish() as u128;
    let lo = RandomState::new().build_hasher().finish() as u128;
    (hi << 64) | lo
}

/// Linear projection of the unsigned symplectic vector of each frame row.
struct SketchFrame {
    x: Vec<u128>,
    z: Vec<u128>,
}

impl SketchFrame {
    fn new(labels: &[(u128, u128)]) -> Self {
        Self {
            x: labels.iter().map(|&(x, _)| x).collect(),
            z: labels.iter().map(|&(_, z)| z).collect(),
        }
    }

    fn apply(&mut self, gate: &Gate) {
        match *gate {
            Gate::h(q) => {
                std::mem::swap(&mut self.x[q as usize], &mut self.z[q as usize]);
            }
            Gate::s(q) | Gate::sdg(q) => self.x[q as usize] ^= self.z[q as usize],
            Gate::cnot { control, target } => {
                let (c, t) = (control as usize, target as usize);
                self.x[c] ^= self.x[t];
                self.z[t] ^= self.z[c];
            }
            Gate::cz { control, target } => {
                let (c, t) = (control as usize, target as usize);
                self.x[c] ^= self.z[t];
                self.x[t] ^= self.z[c];
            }
            Gate::x(_) | Gate::z(_) => {} // only the exact row sign changes
            _ => unreachable!("supported Clifford gate"),
        }
    }
}

/// `i^phase X^x Z^z`; the first `l` words are x, the next `l` are z.
struct Row {
    words: Vec<u64>,
    phase: u8,
}

impl Row {
    fn single(nwords: usize, q: usize, x: bool) -> Self {
        let mut words = vec![0; 2 * nwords];
        let offset = if x { 0 } else { nwords };
        words[offset + q / 64] = 1 << (q % 64);
        Self { words, phase: 0 }
    }

    /// Ordered Pauli product `self * rhs`.
    fn multiply(&mut self, rhs: &Self, l: usize) {
        let mut parity = 0u32;
        for w in 0..l {
            parity ^= (self.words[l + w] & rhs.words[w]).count_ones() & 1;
        }
        self.phase = (self.phase + rhs.phase + 2 * parity as u8) & 3;
        for (a, b) in self.words.iter_mut().zip(&rhs.words) {
            *a ^= *b;
        }
    }

    fn sign(&self, l: usize) -> i8 {
        let y = (0..l)
            .map(|w| (self.words[w] & self.words[l + w]).count_ones())
            .sum::<u32>() as u8
            & 3;
        match (self.phase + 4 - y) & 3 {
            0 => 1,
            2 => -1,
            _ => unreachable!("Clifford conjugation preserves Hermiticity"),
        }
    }
}

fn two_mut<T>(rows: &mut [T], a: usize, b: usize) -> (&mut T, &mut T) {
    debug_assert_ne!(a, b);
    if a < b {
        let (left, right) = rows.split_at_mut(b);
        (&mut left[a], &mut right[0])
    } else {
        let (left, right) = rows.split_at_mut(a);
        (&mut right[0], &mut left[b])
    }
}

struct ExactFrame {
    x: Vec<Row>,
    z: Vec<Row>,
    l: usize,
}

impl ExactFrame {
    fn new(n: usize) -> Self {
        let l = n.div_ceil(64).max(1);
        Self {
            x: (0..n).map(|q| Row::single(l, q, true)).collect(),
            z: (0..n).map(|q| Row::single(l, q, false)).collect(),
            l,
        }
    }

    fn apply(&mut self, gate: &Gate) {
        let l = self.l;
        match *gate {
            Gate::h(q) => {
                let q = q as usize;
                std::mem::swap(&mut self.x[q], &mut self.z[q]);
            }
            Gate::x(q) => self.z[q as usize].phase = (self.z[q as usize].phase + 2) & 3,
            Gate::z(q) => self.x[q as usize].phase = (self.x[q as usize].phase + 2) & 3,
            Gate::s(q) | Gate::sdg(q) => {
                let q = q as usize;
                self.x[q].multiply(&self.z[q], l);
                self.x[q].phase =
                    (self.x[q].phase + if matches!(gate, Gate::s(_)) { 3 } else { 1 }) & 3;
            }
            Gate::cnot { control, target } => {
                let (c, t) = (control as usize, target as usize);
                let (xc, xt) = two_mut(&mut self.x, c, t);
                xc.multiply(xt, l);
                let (zt, zc) = two_mut(&mut self.z, t, c);
                zt.multiply(zc, l); // these two rows commute
            }
            Gate::cz { control, target } => {
                let (c, t) = (control as usize, target as usize);
                self.x[c].multiply(&self.z[t], l);
                self.x[t].multiply(&self.z[c], l); // these two rows commute
            }
            _ => unreachable!("supported Clifford gate"),
        }
    }
}

/// A live or folded rotation. Rotations that cannot be rewritten (the seven
/// Pauli rotations a CCX or CCZ implies) are recorded as blockers with no
/// gate: they are never merge candidates, only obstacles.
#[derive(Clone, Copy)]
struct Event {
    gate: Option<usize>,
    offset: usize,
    sign: i8,
    angle: Angle,
    prev_hash: Option<usize>,
    live: bool,
}

/// Follow an axis fingerprint's chain past rotations already folded away.
/// Compressing dead links prevents long runs of successful folds from making
/// each subsequent lookup revisit all earlier dead events.
fn live_hash(events: &mut [Event], start: Option<usize>) -> Option<usize> {
    let mut cursor = start;
    while let Some(id) = cursor {
        if events[id].live {
            break;
        }
        cursor = events[id].prev_hash;
    }
    let result = cursor;
    let mut cursor = start;
    while let Some(id) = cursor {
        if events[id].live {
            break;
        }
        cursor = events[id].prev_hash;
        events[id].prev_hash = result;
    }
    result
}

/// The rotations and blockers recorded so far, with an index by support.
struct History {
    l: usize,
    /// Work allowed per fold attempt.
    attempt_steps: usize,
    events: Vec<Event>,
    /// Exact unsigned axes, `2 l` words per event.
    axes: Vec<u64>,
    /// Per support-signature bit, the events whose axis touches it, in order.
    /// Only these can anticommute with an axis whose signature has that bit.
    touching: [Vec<u32>; 64],
    /// Newest event per axis fingerprint (merge candidates only).
    last: FxHashMap<u128, usize>,
}

impl History {
    fn new(l: usize, attempt_steps: usize) -> Self {
        Self {
            l,
            attempt_steps,
            events: Vec::new(),
            axes: Vec::new(),
            touching: std::array::from_fn(|_| Vec::new()),
            last: FxHashMap::default(),
        }
    }

    fn axis(&self, id: usize) -> &[u64] {
        let offset = self.events[id].offset;
        &self.axes[offset..offset + 2 * self.l]
    }

    /// Record a rotation; `hash` makes it a merge candidate. `None` when the
    /// axis storage cap is reached.
    fn push(&mut self, event: Event, axis: &[u64], hash: Option<u128>) -> Option<()> {
        if self.axes.len().saturating_add(axis.len()) > MAX_AXIS_WORDS {
            return None;
        }
        let id = self.events.len();
        let mut event = Event {
            offset: self.axes.len(),
            prev_hash: None,
            ..event
        };
        if let Some(h) = hash {
            event.prev_hash = self.last.insert(h, id);
        }
        self.axes.extend_from_slice(axis);
        for bit in bits(support_signature(axis, self.l)) {
            self.touching[bit].push(id as u32);
        }
        self.events.push(event);
        Some(())
    }

    /// Record a fixed obstacle: a later fold across it must commute with
    /// `axis`. `None` when the axis storage cap is reached.
    fn block(&mut self, axis: &[u64]) -> Option<()> {
        let blocker = Event {
            gate: None,
            offset: 0,
            sign: 1,
            angle: Angle {
                eighths: 1,
                residual: 0.0,
            },
            prev_hash: None,
            live: true,
        };
        self.push(blocker, axis, None)
    }

    /// The newest live candidate with exactly this axis, within the attempt
    /// budget (`None` if it runs out).
    fn candidate(&mut self, hash: u128, axis: &[u64], steps: &mut usize) -> Option<Option<usize>> {
        let mut cursor = live_hash(&mut self.events, self.last.get(&hash).copied());
        while let Some(id) = cursor {
            *steps += 1;
            if *steps > self.attempt_steps {
                return None;
            }
            if self.axis(id) == axis {
                return Some(Some(id));
            }
            let prev = self.events[id].prev_hash;
            cursor = live_hash(&mut self.events, prev);
        }
        Some(None)
    }

    /// Whether `axis` commutes with every live rotation recorded after event
    /// `after`. Only rotations sharing a support-signature bit can fail to
    /// commute, so only those are checked. `None` if the attempt budget runs
    /// out first.
    fn commutes_after(&self, after: usize, axis: &[u64], steps: &mut usize) -> Option<bool> {
        for bit in bits(support_signature(axis, self.l)) {
            let list = &self.touching[bit];
            let start = list.partition_point(|&id| id as usize <= after);
            for &id in &list[start..] {
                *steps += 1;
                if *steps > self.attempt_steps {
                    return None;
                }
                let id = id as usize;
                if self.events[id].live && packed_anticommutes(self.axis(id), axis, self.l) {
                    return Some(false);
                }
            }
        }
        Some(true)
    }
}

fn bits(mut w: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        (w != 0).then(|| {
            let bit = w.trailing_zeros() as usize;
            w &= w - 1;
            bit
        })
    })
}

/// Keep exact multiples of pi/4 separate from arbitrary floating-point angles.
/// This preserves exact Clifford+T folding even after many merges.
#[derive(Clone, Copy)]
struct Angle {
    eighths: u8,
    residual: f64,
}

/// Only turn an Rz into an exact multiple of pi/4 when its stored f64 value is
/// itself that value. A tolerance here could erase a small but real rotation.
fn exact_eighths(theta: f64) -> Option<u8> {
    let q = PI / 4.0;
    let k = (theta / q).round();
    if k.abs() <= 1024.0 && theta == k * q {
        Some((k as i32).rem_euclid(8) as u8)
    } else {
        None
    }
}

impl Angle {
    fn from_gate(gate: &Gate) -> Self {
        match *gate {
            Gate::t(_) => Self {
                eighths: 1,
                residual: 0.0,
            },
            Gate::tdg(_) => Self {
                eighths: 7,
                residual: 0.0,
            },
            Gate::rz(theta, _) => match exact_eighths(theta) {
                Some(eighths) => Self {
                    eighths,
                    residual: 0.0,
                },
                None => Self {
                    eighths: 0,
                    residual: theta,
                },
            },
            _ => unreachable!("rotation gate"),
        }
    }

    /// Express two rotations about signed copies of one axis at the later gate.
    fn merge(prior: Self, prior_sign: i8, current: Self, current_sign: i8) -> Option<Self> {
        let relative_sign = i32::from(prior_sign * current_sign);
        let eighths = (relative_sign * i32::from(prior.eighths) + i32::from(current.eighths))
            .rem_euclid(8) as u8;
        let residual = f64::from(prior_sign * current_sign) * prior.residual + current.residual;
        if !residual.is_finite() {
            return None;
        }
        if residual == 0.0 {
            return Some(Self {
                eighths,
                residual: 0.0,
            });
        }
        let total = f64::from(eighths) * (PI / 4.0) + residual;
        if !total.is_finite() {
            return None;
        }
        Some(match exact_eighths(residual) {
            Some(extra) => Self {
                eighths: (eighths + extra) & 7,
                residual: 0.0,
            },
            None => Self { eighths, residual },
        })
    }

    fn clifford_gate(self, q: u32) -> Option<Gate> {
        if self.residual != 0.0 {
            return None;
        }
        match self.eighths {
            2 => Some(Gate::s(q)),
            4 => Some(Gate::z(q)),
            6 => Some(Gate::sdg(q)),
            _ => None,
        }
    }

    fn is_clifford(self) -> bool {
        self.residual == 0.0 && self.eighths & 1 == 0
    }

    fn emit(self, output: &mut Circuit, q: u32) {
        if self.residual != 0.0 {
            let theta = f64::from(self.eighths) * (PI / 4.0) + self.residual;
            output.apply(Gate::rz(theta.rem_euclid(2.0 * PI), q));
            return;
        }
        match self.eighths {
            0 => {}
            1 => output.apply(Gate::t(q)),
            2 => output.apply(Gate::s(q)),
            3 => {
                output.apply(Gate::s(q));
                output.apply(Gate::t(q));
            }
            4 => output.apply(Gate::z(q)),
            5 => {
                output.apply(Gate::z(q));
                output.apply(Gate::t(q));
            }
            6 => output.apply(Gate::sdg(q)),
            7 => output.apply(Gate::tdg(q)),
            _ => unreachable!(),
        }
    }
}

/// Every gate kind is handled; the input must only be well formed (qubits
/// in range and distinct, finite Rz angles). Anything else is left as is.
fn well_formed(circuit: &Circuit) -> bool {
    let n = circuit.num_qubits;
    circuit.gates.iter().all(|g| {
        if matches!(g, Gate::rz(theta, _) if !theta.is_finite()) {
            return false;
        }
        let mut qubits = qubits_of(g);
        qubits.sort_unstable();
        qubits.iter().all(|&q| (q as usize) < n) && qubits.windows(2).all(|w| w[0] != w[1])
    })
}

/// The unsigned axes of the seven Pauli rotations a CCX (`target_x`) or CCZ
/// implies, from the current frame rows of its operands: products of the Z
/// images of the controls and the X (CCX) or Z (CCZ) image of the target, as
/// in `to_pbc`. They are fixed obstacles to folding, not candidates.
fn toffoli_axes(frame: &ExactFrame, a: usize, b: usize, t: usize, target_x: bool) -> Vec<Vec<u64>> {
    let rows = [
        &frame.z[a].words,
        &frame.z[b].words,
        if target_x {
            &frame.x[t].words
        } else {
            &frame.z[t].words
        },
    ];
    (1..8u8)
        .map(|mask| {
            let mut axis = vec![0u64; rows[0].len()];
            for (i, row) in rows.iter().enumerate() {
                if mask >> i & 1 != 0 {
                    axis.iter_mut().zip(row.iter()).for_each(|(a, r)| *a ^= r);
                }
            }
            axis
        })
        .collect()
}

fn pauli_fold_with_labels(circuit: &Circuit, labels: &[(u128, u128)]) -> Circuit {
    pauli_fold_with(circuit, labels, MAX_ATTEMPT_STEPS)
}

/// [`pauli_fold_with_labels`] with a given per-attempt budget.
fn pauli_fold_with(circuit: &Circuit, labels: &[(u128, u128)], attempt_steps: usize) -> Circuit {
    if labels.len() != circuit.num_qubits || !well_formed(circuit) {
        return circuit.clone();
    }
    let num_rotations = circuit
        .gates
        .iter()
        .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_) | Gate::rz(..)))
        .count();
    if num_rotations < 2 {
        return circuit.clone();
    }

    // No exact tableau or stored axes on the common no-candidate path.
    // Measurements, resets, CCX and CCZ change no frame row; rotations may
    // fold across the first two (see the scan below), so they do not reset
    // the candidates.
    let mut sketch = SketchFrame::new(labels);
    let mut seen = FxHashSet::default();
    let mut possible = false;
    for gate in &circuit.gates {
        match *gate {
            Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => {
                let angle = Angle::from_gate(gate);
                if matches!(gate, Gate::rz(..)) && angle.is_clifford() {
                    if let Some(clifford) = angle.clifford_gate(q) {
                        sketch.apply(&clifford);
                    }
                    continue;
                }
                if !seen.insert(sketch.z[q as usize]) {
                    possible = true;
                    break;
                }
            }
            Gate::measure { .. } | Gate::reset(_) | Gate::ccx { .. } | Gate::ccz { .. } => {}
            _ => sketch.apply(gate),
        }
    }
    if !possible {
        return circuit.clone();
    }

    // The exact tableau has two rows per qubit. Bound its quadratic storage
    // before allocating it, even for a circuit with very few rotations.
    let l = circuit.num_qubits.div_ceil(64).max(1);
    if circuit
        .num_qubits
        .checked_mul(4 * l)
        .is_none_or(|words| words > MAX_AXIS_WORDS)
    {
        return circuit.clone();
    }

    let mut sketch = SketchFrame::new(labels);
    let mut exact = ExactFrame::new(circuit.num_qubits);
    let mut history = History::new(l, attempt_steps);
    // 0: keep; 1: delete; >=2: index into replacements.
    let mut edits = vec![0u32; circuit.gates.len()];
    let mut replacements: Vec<Angle> = Vec::new();
    let mut rewrites = 0;

    'scan: for (idx, gate) in circuit.gates.iter().enumerate() {
        match *gate {
            Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => {
                let angle = Angle::from_gate(gate);
                if matches!(gate, Gate::rz(..)) && angle.is_clifford() {
                    if let Some(clifford) = angle.clifford_gate(q) {
                        sketch.apply(&clifford);
                        exact.apply(&clifford);
                    }
                    continue;
                }
                let h = sketch.z[q as usize];
                let row = &exact.z[q as usize];
                let sign = row.sign(l);
                let axis = row.words.clone();
                // One attempt: find the newest candidate with this exact axis
                // and check that nothing live since it anticommutes. Running
                // out of budget counts as blocked.
                let mut steps = 0;
                let reachable = history
                    .candidate(h, &axis, &mut steps)
                    .flatten()
                    .filter(|&id| history.commutes_after(id, &axis, &mut steps) == Some(true));
                let mut event = Event {
                    gate: Some(idx),
                    offset: 0,
                    sign,
                    angle,
                    prev_hash: None,
                    live: true,
                };
                if let Some(id) = reachable {
                    // With P_i = s_i W and P_j = s_j W, commuting P_i through
                    // the live interval gives an angle s_i*theta_i +
                    // s_j*theta_j on W. At j, the physical Z gate needs the
                    // angle multiplied by s_j.
                    let prior = history.events[id];
                    if let Some(merged) = Angle::merge(prior.angle, prior.sign, angle, sign) {
                        if replacements.len() >= (u32::MAX - 2) as usize {
                            break 'scan;
                        }
                        let source = prior.gate.expect("candidates are rewritable rotations");
                        edits[source] = 1;
                        history.events[id].live = false;
                        rewrites += 1;
                        edits[idx] = replacements.len() as u32 + 2;
                        replacements.push(merged);
                        if merged.is_clifford() {
                            // S, S† or Z now sits here and changes every later
                            // axis; the identity leaves nothing.
                            if let Some(emitted) = merged.clifford_gate(q) {
                                sketch.apply(&emitted);
                                exact.apply(&emitted);
                            }
                            continue;
                        }
                        // A non-Clifford merged rotation stays live and may
                        // itself merge with a later rotation on this axis.
                        event.angle = merged;
                    }
                }
                if history.push(event, &axis, Some(h)).is_none() {
                    break 'scan;
                }
            }
            Gate::ccx {
                control1,
                control2,
                target,
            }
            | Gate::ccz {
                control1,
                control2,
                target,
            } => {
                // Non-Clifford, so the frame is unchanged; its rotations stay
                // in the circuit and block folds they anticommute with.
                let target_x = matches!(gate, Gate::ccx { .. });
                let axes = toffoli_axes(
                    &exact,
                    control1 as usize,
                    control2 as usize,
                    target as usize,
                    target_x,
                );
                for axis in axes {
                    if history.block(&axis).is_none() {
                        break 'scan;
                    }
                }
            }
            // A rotation about P commutes with a Z-basis measurement of q,
            // whose projectors are (1 ± Z_q)/2, exactly when P commutes with
            // Z_q. It commutes with a reset of q, whose Kraus maps are
            // |0⟩⟨b|, when P acts trivially on q: P commutes with Z_q and X_q.
            // So each is a blocker on the input-frame images of those
            // Paulis, and only folds that cross it with an anticommuting axis
            // are refused. The frame stays a valid reference for later axes:
            // folds only move a rotation forward, and both axes are compared
            // after the same Clifford prefix.
            Gate::measure { qubit, .. } => {
                if history.block(&exact.z[qubit as usize].words).is_none() {
                    break 'scan;
                }
            }
            Gate::reset(q) => {
                let q = q as usize;
                if history.block(&exact.z[q].words).is_none()
                    || history.block(&exact.x[q].words).is_none()
                {
                    break 'scan;
                }
            }
            _ => {
                sketch.apply(gate);
                exact.apply(gate);
            }
        }
    }

    if rewrites == 0 {
        return circuit.clone();
    }
    let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
    output.gates.reserve(circuit.gates.len() - rewrites);
    for (gate, edit) in circuit.gates.iter().cloned().zip(edits) {
        match edit {
            0 => output.apply(gate),
            1 => {}
            edit => {
                let q = match gate {
                    Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => q,
                    _ => unreachable!(),
                };
                replacements[(edit - 2) as usize].emit(&mut output, q);
            }
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cancel::CancelGates;
    use crate::pass::{count_rz, count_t};
    use crate::pbc::{Pauli, Phase, to_pbc};
    use crate::phase_fold_rand::phase_fold_rand;
    use crate::qasm;
    use crate::unitary::{C, circuit_unitary, circuits_equiv};

    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn up_to(&mut self, n: usize) -> usize {
            self.next() as usize % n
        }
    }

    #[test]
    fn also_folds_without_h_when_used_standalone() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(out.gates, vec![Gate::s(0)]);
    }

    #[test]
    fn folds_across_h_with_commuting_rotation() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(1));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_t(&out), 1);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn anticommuting_rotation_blocks_fold() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_t(&out), 3);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn opposite_sign_axes_cancel() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::h(1));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_t(&out), 0);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn merges_rz_and_t_across_commuting_rotations() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::rz(0.37, 1));
        c.apply(Gate::h(0));
        c.apply(Gate::rz(0.21, 0));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_t(&out), 0);
        assert_eq!(count_rz(&out), 2);
        assert_eq!(out.gates.len(), c.gates.len() - 1);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn opposite_sign_rz_and_t_cancel() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI / 4.0, 0));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(out.gates, vec![Gate::x(0)]);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn rz_can_fold_to_clifford_and_change_the_frame() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI / 4.0, 0));
        c.apply(Gate::rz(PI / 4.0, 0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_rz(&out), 0);
        assert_eq!(out.gates.first(), Some(&Gate::s(0)));
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn anticommuting_rz_blocks_t_fold() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::rz(0.3, 0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(out.gates, c.gates);
    }

    #[test]
    fn clifford_rz_updates_frame_and_zero_rz_does_not_block() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::rz(0.0, 0));
        c.apply(Gate::h(0));
        c.apply(Gate::rz(PI / 2.0, 0));
        c.apply(Gate::t(0));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_t(&out), 0);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn near_quarter_rz_is_not_rounded_away() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI / 4.0 + 1e-10, 0));
        c.apply(Gate::rz(PI / 4.0, 0));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_rz(&out), 1);
        assert_eq!(count_t(&out), 0);
        assert!(circuits_equiv(&c, &out, 1e-12));
    }

    #[test]
    fn rz_axis_across_packed_word_boundary() {
        let mut c = Circuit::new(65);
        let cx = Gate::cnot {
            control: 0,
            target: 64,
        };
        c.apply(Gate::rz(0.3, 64));
        c.apply(cx.clone());
        c.apply(Gate::rz(0.2, 64));
        c.apply(cx);
        c.apply(Gate::rz(0.4, 64));
        let out = pauli_fold_rand(&c);
        assert_eq!(count_rz(&out), 2);
        assert_eq!(out.gates.len(), 4);
        assert!(
            matches!(out.gates.last(), Some(Gate::rz(theta, 64)) if (*theta - 0.7).abs() < 1e-12)
        );
    }

    #[test]
    fn forced_fingerprint_collisions_cannot_authorize_rewrites() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = pauli_fold_with_labels(&c, &[(0, 0)]);
        assert_eq!(count_t(&out), 2);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn mod5_4_improves_on_phase_fold_and_preserves_unitary() {
        let input = qasm::parse(include_str!("../benchmarks/feynman/mod5_4.qasm")).unwrap();
        let cancelled = CancelGates.run(&input);
        let baseline = phase_fold_rand(&cancelled);
        let improved = pauli_fold_rand(&baseline);
        assert!(count_t(&improved) < count_t(&baseline));
        assert!(circuits_equiv(&baseline, &improved, 1e-9));
    }

    #[test]
    fn exact_frame_agrees_with_pbc_frame() {
        let mut rng = Rng(0x7e57_1234_9abc_def0);
        let mut c = Circuit::new(4);
        let mut frame = ExactFrame::new(4);
        let labels: Vec<_> = (0..4)
            .map(|_| (rng.next() as u128, rng.next() as u128))
            .collect();
        let mut sketch = SketchFrame::new(&labels);
        for _ in 0..100 {
            let q = rng.up_to(4) as u32;
            let other = (q + 1 + rng.up_to(3) as u32) % 4;
            let gate = match rng.up_to(8) {
                0 => Gate::h(q),
                1 => Gate::x(q),
                2 => Gate::z(q),
                3 => Gate::s(q),
                4 => Gate::sdg(q),
                5 => Gate::cnot {
                    control: q,
                    target: other,
                },
                _ => Gate::cz {
                    control: q,
                    target: other,
                },
            };
            c.apply(gate.clone());
            frame.apply(&gate);
            sketch.apply(&gate);
            let pbc = to_pbc(&c, None).unwrap();
            for r in 0..4u32 {
                for (row, hash, reference) in [
                    (
                        &frame.x[r as usize],
                        sketch.x[r as usize],
                        pbc.output_frame().x(r),
                    ),
                    (
                        &frame.z[r as usize],
                        sketch.z[r as usize],
                        pbc.output_frame().z(r),
                    ),
                ] {
                    let expanded = pbc.expand(reference, 100_000).unwrap();
                    assert_eq!(
                        row.sign(frame.l),
                        if expanded.phase == Phase::One { 1 } else { -1 }
                    );
                    let mut projected = 0;
                    for (bit, factor) in expanded.factors.iter().enumerate() {
                        let x = (row.words[bit / 64] >> (bit % 64)) & 1 != 0;
                        let z = (row.words[frame.l + bit / 64] >> (bit % 64)) & 1 != 0;
                        if x {
                            projected ^= labels[bit].0;
                        }
                        if z {
                            projected ^= labels[bit].1;
                        }
                        assert_eq!(
                            (x, z),
                            match factor {
                                Pauli::I => (false, false),
                                Pauli::X => (true, false),
                                Pauli::Y => (true, true),
                                Pauli::Z => (false, true),
                            }
                        );
                    }
                    assert_eq!(hash, projected);
                }
            }
        }
    }

    #[test]
    fn random_circuits_preserve_unitary_and_baseline_counts() {
        let mut rng = Rng(0x9d61_a45e_1754_3cba);
        for case in 0..2_000 {
            let n = if case % 4 == 0 { 4 } else { 3 };
            let mut c = Circuit::new(n);
            for _ in 0..40 {
                let q = rng.up_to(n) as u32;
                let other = (q + 1 + rng.up_to(n - 1) as u32) % n as u32;
                let gate = match rng.up_to(12) {
                    0 => Gate::h(q),
                    1 => Gate::x(q),
                    2 => Gate::z(q),
                    3 => Gate::s(q),
                    4 => Gate::sdg(q),
                    5..=7 => Gate::t(q),
                    8 => Gate::tdg(q),
                    9 | 10 => Gate::cnot {
                        control: q,
                        target: other,
                    },
                    _ => Gate::cz {
                        control: q,
                        target: other,
                    },
                };
                c.apply(gate);
            }
            let baseline = phase_fold_rand(&c);
            let out = pauli_fold_rand(&baseline);
            assert!(count_t(&out) <= count_t(&baseline));
            assert!(out.gates.len() <= baseline.gates.len());
            assert!(circuits_equiv(&baseline, &out, 1e-9));
            if case < 100 {
                let collided = pauli_fold_with_labels(&baseline, &vec![(0, 0); n]);
                assert!(circuits_equiv(&baseline, &collided, 1e-9));
            }
        }
    }

    #[test]
    fn random_rz_circuits_preserve_unitary_and_rotation_count() {
        let mut rng = Rng(0x76a5_371e_c118_abd9);
        for case in 0..800 {
            let n = if case % 3 == 0 { 4 } else { 3 };
            let mut c = Circuit::new(n);
            for _ in 0..35 {
                let q = rng.up_to(n) as u32;
                let other = (q + 1 + rng.up_to(n - 1) as u32) % n as u32;
                let gate = match rng.up_to(16) {
                    0 => Gate::h(q),
                    1 => Gate::x(q),
                    2 => Gate::z(q),
                    3 => Gate::s(q),
                    4 => Gate::sdg(q),
                    5..=6 => Gate::t(q),
                    7 => Gate::tdg(q),
                    8..=10 => Gate::rz((rng.up_to(17) as f64 - 8.0) / 13.0, q),
                    11 => Gate::rz((rng.up_to(9) as f64 - 4.0) * (PI / 4.0), q),
                    12 | 13 => Gate::cnot {
                        control: q,
                        target: other,
                    },
                    _ => Gate::cz {
                        control: q,
                        target: other,
                    },
                };
                c.apply(gate);
            }
            let out = pauli_fold_rand(&c);
            assert!(count_t(&out) + count_rz(&out) <= count_t(&c) + count_rz(&c));
            assert!(out.gates.len() <= c.gates.len());
            assert!(circuits_equiv(&c, &out, 1e-9), "case {case}");
            if case < 100 {
                let labels = vec![(0, 0); n];
                let collided = pauli_fold_with_labels(&c, &labels);
                assert!(circuits_equiv(&c, &collided, 1e-9), "collision case {case}");
            }
        }
    }

    fn ccx(a: u32, b: u32, t: u32) -> Gate {
        Gate::ccx {
            control1: a,
            control2: b,
            target: t,
        }
    }

    fn ccz(a: u32, b: u32, t: u32) -> Gate {
        Gate::ccz {
            control1: a,
            control2: b,
            target: t,
        }
    }

    /// A CCX's rotations include X on its target, which anticommutes with Z:
    /// T gates on the target cannot fold across it.
    #[test]
    fn ccx_blocks_folds_on_its_target() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(2));
        c.apply(ccx(0, 1, 2));
        c.apply(Gate::t(2));
        let out = pauli_fold_rand(&c);
        assert_eq!(out.gates, c.gates);
    }

    /// Its control axes are Z only, so T gates on a control fold across it,
    /// and a CCZ is diagonal, so T gates fold across it on every operand.
    #[test]
    fn folds_across_toffolis_where_the_axes_commute() {
        for (gate, q) in [(ccx(0, 1, 2), 0), (ccx(0, 1, 2), 1), (ccz(0, 1, 2), 2)] {
            let mut c = Circuit::new(3);
            c.apply(Gate::h(2));
            c.apply(Gate::t(q));
            c.apply(gate.clone());
            c.apply(Gate::t(q));
            let out = pauli_fold_rand(&c);
            assert_eq!(count_t(&out), 0, "{gate:?} on q{q}");
            assert!(out.gates.contains(&gate));
            assert!(circuits_equiv(&c, &out, 1e-10));
        }
    }

    /// The Toffoli's blocking axes use the frame at its position: after an H
    /// on the target, a CCX's target axis is Z, so T gates there fold across.
    #[test]
    fn toffoli_axes_follow_the_frame() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(2));
        c.apply(Gate::t(2));
        c.apply(Gate::h(2));
        c.apply(ccx(0, 1, 2));
        c.apply(Gate::h(2));
        c.apply(Gate::t(2));
        let out = pauli_fold_rand(&c);
        // Both T gates act on X2 in the input frame; the CCX's rotations
        // commute with X2 exactly when its target axis is X2 -- here it is.
        assert!(circuits_equiv(&c, &out, 1e-10));
        assert_eq!(count_t(&out), 0);
    }

    fn circuit_of(n: usize, gates: &[Gate]) -> Circuit {
        let mut c = Circuit::with_cbits(n, n);
        for gate in gates {
            c.apply(gate.clone());
        }
        c
    }

    /// A rotation folds across a measurement or reset it commutes with, on
    /// another qubit or, for a measurement, about Z on the measured qubit.
    #[test]
    fn folds_cross_commuting_measurements_and_resets() {
        let cnot = |control, target| Gate::cnot { control, target };
        for (name, gates) in [
            (
                "measure other",
                vec![
                    Gate::t(0),
                    Gate::measure { qubit: 1, cbit: 1 },
                    Gate::tdg(0),
                ],
            ),
            (
                "measure same",
                vec![
                    Gate::t(0),
                    Gate::measure { qubit: 0, cbit: 0 },
                    Gate::tdg(0),
                ],
            ),
            (
                "reset other",
                vec![Gate::t(0), Gate::reset(1), Gate::tdg(0)],
            ),
            (
                "measure Z0 under Z0Z1",
                vec![
                    cnot(0, 1),
                    Gate::t(1),
                    cnot(0, 1),
                    Gate::measure { qubit: 0, cbit: 0 },
                    cnot(0, 1),
                    Gate::tdg(1),
                    cnot(0, 1),
                ],
            ),
        ] {
            let c = circuit_of(2, &gates);
            let out = pauli_fold_rand(&c);
            assert_eq!(count_t(&out), 0, "{name}: {:?}", out.gates);
        }
    }

    /// A fold whose axis anticommutes with a measured Z, or touches a reset
    /// qubit, is refused; the rotations on either side still fold.
    #[test]
    fn measurements_and_resets_block_anticommuting_folds() {
        let cnot = |control, target| Gate::cnot { control, target };
        for (name, gates) in [
            (
                "X0 across measure q0",
                vec![
                    Gate::h(0),
                    Gate::t(0),
                    Gate::h(0),
                    Gate::measure { qubit: 0, cbit: 0 },
                    Gate::h(0),
                    Gate::tdg(0),
                    Gate::h(0),
                ],
            ),
            (
                "Z0 across reset q0",
                vec![Gate::t(0), Gate::reset(0), Gate::tdg(0)],
            ),
            (
                "Z0Z1 across reset q0",
                vec![
                    cnot(0, 1),
                    Gate::t(1),
                    cnot(0, 1),
                    Gate::reset(0),
                    cnot(0, 1),
                    Gate::tdg(1),
                    cnot(0, 1),
                ],
            ),
        ] {
            let c = circuit_of(2, &gates);
            let out = pauli_fold_rand(&c);
            assert_eq!(count_t(&out), 2, "{name}: {:?}", out.gates);
        }
        // T T T | reset q0 | T T: S T before, S after.
        let c = circuit_of(
            1,
            &[
                Gate::t(0),
                Gate::t(0),
                Gate::t(0),
                Gate::reset(0),
                Gate::t(0),
                Gate::t(0),
            ],
        );
        let out = pauli_fold_rand(&c);
        assert_eq!(count_t(&out), 1, "{:?}", out.gates);
        let at = out.gates.iter().position(|g| *g == Gate::reset(0)).unwrap();
        assert_eq!(
            count_t(&Circuit {
                gates: out.gates[at..].to_vec(),
                ..c.clone()
            }),
            0
        );
    }

    /// An attempt whose interval exceeds its budget is treated as blocked;
    /// later folds are still found.
    #[test]
    fn an_exhausted_attempt_does_not_stop_the_pass() {
        const BUDGET: usize = 64;
        let n = 20;
        let mut c = Circuit::new(n);
        c.apply(Gate::t(0));
        // More distinct rotations Z0 Za Zb, all commuting with Z0 and
        // touching qubit 0's support bit, than one attempt may scan.
        let pairs = (1..n as u32).flat_map(|a| ((a + 1)..n as u32).map(move |b| (a, b)));
        for (a, b) in pairs.take(2 * BUDGET) {
            for gate in [
                Gate::cnot {
                    control: a,
                    target: 0,
                },
                Gate::cnot {
                    control: b,
                    target: 0,
                },
                Gate::t(0),
                Gate::cnot {
                    control: b,
                    target: 0,
                },
                Gate::cnot {
                    control: a,
                    target: 0,
                },
            ] {
                c.apply(gate);
            }
        }
        c.apply(Gate::t(0));
        c.apply(Gate::t(1));
        c.apply(Gate::t(1));
        let labels: Vec<_> = (0..n)
            .map(|q| (q as u128 * 2 + 1, q as u128 * 2 + 2))
            .collect();
        let out = pauli_fold_with(&c, &labels, BUDGET);
        // The Z0 pair is out of reach within one attempt; the Z1 pair folds.
        assert_eq!(count_t(&out), count_t(&c) - 2);
        assert_eq!(out.gates.last(), Some(&Gate::s(1)));
        // With the default budget the Z0 pair folds too.
        let full = pauli_fold_with_labels(&c, &labels);
        assert_eq!(count_t(&full), count_t(&c) - 4);
    }

    fn random_clifford(rng: &mut Rng, n: usize, len: usize) -> Vec<Gate> {
        (0..len)
            .map(|_| {
                let q = rng.up_to(n) as u32;
                let r = (q + 1 + rng.up_to(n - 1) as u32) % n as u32;
                match rng.up_to(5) {
                    0 | 1 => Gate::h(q),
                    2 => Gate::s(q),
                    3 => Gate::cnot {
                        control: q,
                        target: r,
                    },
                    _ => Gate::cz {
                        control: q,
                        target: r,
                    },
                }
            })
            .collect()
    }

    fn inverse(gates: &[Gate]) -> Vec<Gate> {
        gates
            .iter()
            .rev()
            .map(|g| match *g {
                Gate::s(q) => Gate::sdg(q),
                Gate::sdg(q) => Gate::s(q),
                ref g => g.clone(),
            })
            .collect()
    }

    /// `prefix · T_q · D · barrier · D⁻¹ · T_q^± · suffix`: the two rotations
    /// share an axis, and the barrier decides whether they may fold.
    fn sandwich(rng: &mut Rng, n: usize, barrier: Gate) -> Circuit {
        let q = rng.up_to(n) as u32;
        let (before, middle, after) = (rng.up_to(4), 1 + rng.up_to(5), rng.up_to(4));
        let prefix = random_clifford(rng, n, before);
        let d = random_clifford(rng, n, middle);
        let suffix = random_clifford(rng, n, after);
        let second = if rng.up_to(2) == 0 {
            Gate::t(q)
        } else {
            Gate::tdg(q)
        };
        let mut c = Circuit::with_cbits(n, n);
        let pieces = [
            prefix,
            vec![Gate::t(q)],
            d.clone(),
            vec![barrier],
            inverse(&d),
            vec![second],
            suffix,
        ];
        for gate in pieces.into_iter().flatten() {
            c.apply(gate);
        }
        c
    }

    fn random_gate(rng: &mut Rng, n: usize, measured: bool) -> Gate {
        let q = rng.up_to(n) as u32;
        let others: Vec<u32> = (0..n as u32).filter(|&r| r != q).collect();
        let r = others[rng.up_to(others.len())];
        let s = others.iter().copied().find(|&s| s != r).unwrap_or(r);
        match rng.up_to(if measured { 15 } else { 14 }) {
            0 => Gate::h(q),
            1 => Gate::x(q),
            2 => Gate::z(q),
            3 => Gate::s(q),
            4 => Gate::sdg(q),
            5..=7 => Gate::t(q),
            8 => Gate::tdg(q),
            9 | 10 => Gate::cnot {
                control: q,
                target: r,
            },
            11 => Gate::cz {
                control: q,
                target: r,
            },
            12 if n >= 3 => ccx(q, r, s),
            13 if n >= 3 => ccz(q, r, s),
            14 => Gate::measure { qubit: q, cbit: q },
            _ => Gate::t(q),
        }
    }

    /// Random circuits with native CCX and CCZ keep their unitary, including
    /// with forced fingerprint collisions.
    #[test]
    fn random_toffoli_circuits_preserve_unitary() {
        let mut rng = Rng(0x7f0f_f011_ccc0_0001);
        let mut folded = 0;
        for case in 0..1_000 {
            let n = 3 + case % 2;
            let mut c = Circuit::new(n);
            for _ in 0..40 {
                c.apply(random_gate(&mut rng, n, false));
            }
            let out = pauli_fold_rand(&c);
            assert!(count_t(&out) <= count_t(&c));
            assert!(circuits_equiv(&c, &out, 1e-9), "case {case}\n{c}");
            folded += usize::from(out.gates != c.gates);
            if case < 100 {
                let collided = pauli_fold_with_labels(&c, &vec![(0, 0); n]);
                assert!(circuits_equiv(&c, &collided, 1e-9), "collision case {case}");
            }
        }
        assert!(folded > 500, "only {folded} cases changed");
    }

    /// Random circuits with mid-circuit measurements keep their exact
    /// quantum-classical channel.
    #[test]
    fn random_measured_circuits_preserve_the_channel() {
        use crate::semantics::channel::{ChannelLimits, circuit_channel};
        let mut rng = Rng(0x3ea5_0000_1234_5678);
        let mut folded = 0;
        for case in 0..450 {
            let n = 2 + case % 2;
            let c = if case < 150 {
                let mut c = Circuit::with_cbits(n, n);
                for _ in 0..16 {
                    c.apply(random_gate(&mut rng, n, true));
                }
                c
            } else {
                let m = rng.up_to(n) as u32;
                sandwich(&mut rng, n, Gate::measure { qubit: m, cbit: m })
            };
            let out = pauli_fold_rand(&c);
            folded += usize::from(out.gates != c.gates);
            let store = vec![false; n];
            let limits = ChannelLimits::default();
            let expected = circuit_channel(&c, &store, limits).unwrap();
            let actual = circuit_channel(&out, &store, limits).unwrap();
            assert_eq!(expected.compare(&actual), Ok(()), "case {case}\n{c}");
        }
        assert!(folded > 150, "only {folded} cases changed");
    }

    /// The channel of a circuit whose only non-unitary gates are resets, as
    /// its images of the matrix units |i⟩⟨j|, from dense unitary segments and
    /// the reset map ρ ↦ Σ_b |0⟩⟨b| ρ |b⟩⟨0|.
    fn reset_channel(c: &Circuit) -> Vec<Vec<Vec<C>>> {
        let n = c.num_qubits;
        let d = 1usize << n;
        let mut steps: Vec<Result<Vec<Vec<C>>, usize>> = Vec::new();
        let mut segment = Circuit::new(n);
        for gate in &c.gates {
            if let Gate::reset(q) = *gate {
                steps.push(Ok(circuit_unitary(&segment)));
                steps.push(Err(q as usize));
                segment = Circuit::new(n);
            } else {
                segment.apply(gate.clone());
            }
        }
        steps.push(Ok(circuit_unitary(&segment)));
        let mut images = Vec::with_capacity(d * d);
        for i in 0..d {
            for j in 0..d {
                let mut rho = vec![vec![C::ZERO; d]; d];
                rho[i][j] = C::ONE;
                for step in &steps {
                    rho = match step {
                        Ok(u) => {
                            let mut left = vec![vec![C::ZERO; d]; d];
                            for r in 0..d {
                                for k in 0..d {
                                    for col in 0..d {
                                        left[r][col] = left[r][col] + u[r][k] * rho[k][col];
                                    }
                                }
                            }
                            let mut out = vec![vec![C::ZERO; d]; d];
                            for r in 0..d {
                                for col in 0..d {
                                    for k in 0..d {
                                        out[r][col] = out[r][col] + left[r][k] * u[col][k].conj();
                                    }
                                }
                            }
                            out
                        }
                        Err(q) => {
                            let bit = 1 << (n - 1 - q);
                            let mut out = vec![vec![C::ZERO; d]; d];
                            for r in (0..d).filter(|r| r & bit == 0) {
                                for col in (0..d).filter(|col| col & bit == 0) {
                                    out[r][col] = rho[r][col] + rho[r | bit][col | bit];
                                }
                            }
                            out
                        }
                    };
                }
                images.push(rho);
            }
        }
        images
    }

    /// Random circuits with resets, including Toffolis and forced
    /// fingerprint collisions, keep their exact channel.
    #[test]
    fn random_reset_circuits_preserve_the_channel() {
        let mut rng = Rng(0x7e5e_7000_9abc_def1);
        let mut folded = 0;
        for case in 0..400 {
            let n = 2 + case % 2;
            let c = if case < 200 {
                let mut c = Circuit::new(n);
                for _ in 0..20 {
                    if rng.up_to(8) == 0 {
                        c.apply(Gate::reset(rng.up_to(n) as u32));
                    } else {
                        c.apply(random_gate(&mut rng, n, false));
                    }
                }
                c
            } else {
                let q = rng.up_to(n) as u32;
                sandwich(&mut rng, n, Gate::reset(q))
            };
            let expected = reset_channel(&c);
            let mut outputs = vec![pauli_fold_rand(&c)];
            if case < 50 {
                outputs.push(pauli_fold_with_labels(&c, &vec![(0, 0); n]));
            }
            folded += usize::from(outputs[0].gates != c.gates);
            for out in outputs {
                let actual = reset_channel(&out);
                for (e, a) in expected.iter().zip(&actual) {
                    for (er, ar) in e.iter().zip(a) {
                        for (&x, &y) in er.iter().zip(ar) {
                            assert!((x - y).norm_sq() < 1e-18, "case {case}\n{c}\n{out}");
                        }
                    }
                }
            }
        }
        assert!(folded > 100, "only {folded} cases changed");
    }
}
