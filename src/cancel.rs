//! The `CancelGates` pass: removes adjacent self-inverse gate pairs.

use crate::circuit::{
    Circuit, Gate, Qubit, canonical_ccx, canonical_pair, canonical_triple, qubit_operands,
};
use crate::pass::Pass;
use crate::phase_fold_rand::classify_quarter_pi;

/// Working buffers shared by every sweep of one [`CancelGates::run`].
///
/// The three sweep kinds each want a per-qubit index of gate positions, a
/// circuit-sized deletion mark, and (for Hadamard reduction) a circuit-sized
/// replacement table. Allocating those per sweep was the pass's dominant cost
/// on large inputs: a fixpoint over cobble's `ols-ridge-d12.qasm`
/// runs 76 self-inverse, 55 Hadamard, and 66 commuting-pair sweeps over 33M
/// gate positions in total, and `Gate` is 32 bytes — so the replacement table
/// alone was ~290 MB of freshly allocated, freshly zeroed, freshly
/// page-faulted memory. Sizes only ever shrink within a run (every sweep
/// deletes or rewrites in place, none grows the circuit), so one allocation
/// at the input's size serves all of them.
struct Scratch {
    /// Per-qubit ordered list of positions of the gates touching that qubit.
    tracks: Vec<Vec<usize>>,
    /// Gate-list length `tracks` was built for, or `None` when they are stale.
    /// The Hadamard and commuting-pair sweeps index gates the same way, and a
    /// sweep that reaches its fixpoint returns the list it was given
    /// unchanged — which is the common case — so the next sweep can very often
    /// reuse the index instead of rebuilding it over millions of gates.
    tracks_for: Option<usize>,
    /// Per-qubit stack of live gate positions, for [`cancel_pairs`].
    stacks: Vec<Vec<usize>>,
    /// Deletion marks, held as sweep stamps rather than booleans so that
    /// starting a sweep is a counter bump instead of clearing a
    /// circuit-sized bitmap. `delete[i] == stamp` means "deleted this sweep".
    delete: Vec<u32>,
    stamp: u32,
    /// Pending gate replacements, indexed by gate position. Every sweep that
    /// writes one also `take`s it back out again while emitting its output
    /// (replaced gates are never also deleted), so this returns to all-`None`
    /// on its own and likewise never needs re-zeroing.
    replace: Vec<Option<Gate>>,
    /// The previous sweep's gate buffer, reused as the next sweep's output
    /// instead of allocating and freeing one gate list per sweep.
    spare: Vec<Gate>,
}

impl Scratch {
    fn new(num_qubits: usize, gate_count: usize) -> Self {
        Self {
            tracks: vec![Vec::new(); num_qubits],
            tracks_for: None,
            stacks: vec![Vec::new(); num_qubits],
            delete: vec![0; gate_count],
            stamp: 0,
            replace: vec![None; gate_count],
            spare: Vec::new(),
        }
    }

    /// Begin a sweep, returning the stamp that marks a gate deleted by it.
    /// Starts at 1 so it never collides with the initial all-zero `delete`.
    fn begin(&mut self) -> u32 {
        self.stamp += 1;
        self.stamp
    }

    /// Make the per-qubit position index describe `gates`, rebuilding it only
    /// when a previous sweep invalidated it. Rebuilding keeps the capacity each
    /// track reached on earlier sweeps.
    fn ensure_tracks(&mut self, gates: &[Gate]) {
        if self.tracks_for == Some(gates.len()) {
            return;
        }
        for track in &mut self.tracks {
            track.clear();
        }
        for (i, gate) in gates.iter().enumerate() {
            let (n, qs) = qubit_operands(gate);
            for &q in &qs[..n] {
                self.tracks[q as usize].push(i);
            }
        }
        self.tracks_for = Some(gates.len());
    }

    /// Mark the position index stale, because a sweep is about to hand back a
    /// gate list the index no longer describes.
    fn invalidate_tracks(&mut self) {
        self.tracks_for = None;
    }

    /// An empty gate buffer carrying a previous sweep's capacity.
    fn out_buffer(&mut self) -> Vec<Gate> {
        let mut out = std::mem::take(&mut self.spare);
        out.clear();
        out
    }
}

/// Cancel adjacent self-inverse gate pairs (HH, XX, CNOT-CNOT) in O(n),
/// allowing commutation past gates on non-overlapping qubits.
/// Handles cascading: cancelling a pair may expose new adjacent pairs.
///
/// Uses per-qubit stacks to find the blocking gate in O(1) instead of
/// scanning backward through the entire result list.
/// Gates are tracked by index into the original slice — only surviving gates
/// are cloned at the end.
fn cancel_pairs(gates: &[Gate], scratch: &mut Scratch) -> Vec<Gate> {
    let stamp = scratch.begin();
    let Scratch { stacks, delete, .. } = &mut *scratch;
    for stack in stacks.iter_mut() {
        stack.clear();
    }

    for (i, gate) in gates.iter().enumerate() {
        if is_self_inverse(gate) {
            let (n, qs) = qubit_operands(gate);
            // The blocker is the latest gate touching any of this gate's qubits.
            let mut blocker: Option<usize> = None;
            for j in 0..n {
                if let Some(&last) = stacks[qs[j] as usize].last() {
                    blocker = Some(match blocker {
                        Some(b) => b.max(last),
                        None => last,
                    });
                }
            }
            if let Some(block_idx) = blocker
                && gates_equal(&gates[block_idx], gate)
            {
                // Cancel both gates; pop the blocker from all relevant qubit stacks.
                delete[block_idx] = stamp;
                delete[i] = stamp;
                for j in 0..n {
                    debug_assert_eq!(*stacks[qs[j] as usize].last().unwrap(), block_idx);
                    stacks[qs[j] as usize].pop();
                }
                continue;
            }
        }

        let (n, qs) = qubit_operands(gate);
        for j in 0..n {
            stacks[qs[j] as usize].push(i);
        }
    }

    scratch.invalidate_tracks();
    let mut out = scratch.out_buffer();
    let delete = &scratch.delete;
    out.extend(
        gates
            .iter()
            .enumerate()
            .filter(|(i, _)| delete[*i] != stamp)
            .map(|(_, g)| g.clone()),
    );
    out
}

fn is_self_inverse(gate: &Gate) -> bool {
    matches!(
        gate,
        Gate::h(_)
            | Gate::x(_)
            | Gate::z(_)
            | Gate::cnot { .. }
            | Gate::cz { .. }
            | Gate::ccx { .. }
            | Gate::ccz { .. }
    )
}

fn gates_equal(a: &Gate, b: &Gate) -> bool {
    match (a, b) {
        (Gate::h(a), Gate::h(b)) | (Gate::x(a), Gate::x(b)) | (Gate::z(a), Gate::z(b)) => a == b,
        (
            Gate::cnot {
                control: ac,
                target: at,
            },
            Gate::cnot {
                control: bc,
                target: bt,
            },
        ) => ac == bc && at == bt,
        (
            Gate::cz {
                control: ac,
                target: at,
            },
            Gate::cz {
                control: bc,
                target: bt,
            },
        ) => canonical_pair(*ac, *at) == canonical_pair(*bc, *bt),
        (
            Gate::ccx {
                control1: a1,
                control2: a2,
                target: at,
            },
            Gate::ccx {
                control1: b1,
                control2: b2,
                target: bt,
            },
        ) => canonical_ccx(*a1, *a2, *at) == canonical_ccx(*b1, *b2, *bt),
        (
            Gate::ccz {
                control1: a1,
                control2: a2,
                target: a3,
            },
            Gate::ccz {
                control1: b1,
                control2: b2,
                target: b3,
            },
        ) => canonical_triple(*a1, *a2, *a3) == canonical_triple(*b1, *b2, *b3),
        _ => false,
    }
}

fn is_h(g: &Gate, q: Qubit) -> bool {
    matches!(g, Gate::h(p) if *p == q)
}
fn is_x(g: &Gate, q: Qubit) -> bool {
    matches!(g, Gate::x(p) if *p == q)
}

/// For a single-qubit diagonal gate on `q`, classifies its rotation:
///   `Some(Some(k))` — diagonal, rotates by k·π/4
///   `Some(None)`    — diagonal `rz` whose angle is not a π/4 multiple
///   `None`          — not a single-qubit diagonal gate on `q`
fn diagonal_k(g: &Gate, q: Qubit) -> Option<Option<u32>> {
    let (k, p) = match g {
        Gate::t(p) => (1, p),
        Gate::tdg(p) => (7, p),
        Gate::s(p) => (2, p),
        Gate::sdg(p) => (6, p),
        Gate::z(p) => (4, p),
        Gate::rz(theta, p) => {
            return if *p == q {
                Some(classify_quarter_pi(*theta).map(|k| k as u32))
            } else {
                None
            };
        }
        _ => return None,
    };
    if *p == q { Some(Some(k)) } else { None }
}

/// Local Clifford identities that strictly lower the Hadamard count:
///
///   H·H          = I
///   H·X·H        = Z
///   H·D·H        = X            when the diagonal run D ≡ Z   (mod 2π)
///   H·D·H        = Sdg·H·Sdg    when D ≡ S
///   H·D·H        = S·H·S        when D ≡ Sdg
///   H·D·H        = I            when D ≡ I
///
/// `D` is a maximal run of diagonal gates with no other gate touching the
/// qubit in between, so each identity is local to that wire and sound
/// regardless of the rest of the circuit. Fewer Hadamards means longer
/// Hadamard-free sections for the downstream phase-folding pass to merge
/// rotations across. Every rule removes at least one `h`, so the fixpoint
/// loop terminates. The S-rules differ by a global phase, which the rest of
/// the pipeline (and `circuits_equiv`) ignores.
///
/// Returns the rewritten gate list and whether any rewrite fired.
/// Takes the gate list by value: it used to be cloned on entry purely to have
/// something owned to iterate on, which cost a full copy of the circuit per
/// call even on the (common) call that reaches the fixpoint immediately.
fn reduce_hadamards(mut gates: Vec<Gate>, scratch: &mut Scratch) -> (Vec<Gate>, bool) {
    let mut changed = false;
    while let Some(next) = reduce_hadamards_pass(&gates, scratch) {
        scratch.spare = std::mem::replace(&mut gates, next);
        changed = true;
    }
    (gates, changed)
}

/// One sweep of [`reduce_hadamards`]. Returns the rewritten gate list when at
/// least one rewrite fired, or `None` at the fixpoint.
fn reduce_hadamards_pass(gates: &[Gate], scratch: &mut Scratch) -> Option<Vec<Gate>> {
    let stamp = scratch.begin();
    scratch.ensure_tracks(gates);
    // Per-index edit: a gate is dropped, replaced, or kept. Rewrites on
    // different qubits never share an index, and within a qubit the scan
    // skips past consumed indices, so the edits never conflict.
    let Scratch {
        tracks,
        delete,
        replace,
        ..
    } = &mut *scratch;
    let mut changed = false;

    for (q, track) in tracks.iter().enumerate() {
        let q = q as Qubit;
        let mut p = 0;
        while p < track.len() {
            let io = track[p];
            if !is_h(&gates[io], q) {
                p += 1;
                continue;
            }
            // H·X·H = Z
            if p + 2 < track.len() && is_x(&gates[track[p + 1]], q) && is_h(&gates[track[p + 2]], q)
            {
                delete[io] = stamp;
                replace[track[p + 1]] = Some(Gate::z(q));
                delete[track[p + 2]] = stamp;
                changed = true;
                p += 3;
                continue;
            }
            // H · (maximal diagonal run) · H
            let mut k = 0u32;
            let mut dirty = false;
            let mut j = p + 1;
            while j < track.len() {
                match diagonal_k(&gates[track[j]], q) {
                    Some(Some(v)) => {
                        k = (k + v) & 7;
                        j += 1;
                    }
                    Some(None) => {
                        dirty = true;
                        j += 1;
                    }
                    None => break,
                }
            }
            if j < track.len() && is_h(&gates[track[j]], q) && !dirty && k.is_multiple_of(2) {
                let ic = track[j];
                let run = &track[p + 1..j];
                match k {
                    0 => {
                        delete[io] = stamp;
                        delete[ic] = stamp;
                        for &r in run {
                            delete[r] = stamp;
                        }
                    }
                    4 => {
                        delete[io] = stamp;
                        delete[ic] = stamp;
                        replace[run[0]] = Some(Gate::x(q));
                        for &r in &run[1..] {
                            delete[r] = stamp;
                        }
                    }
                    2 | 6 => {
                        let outer = if k == 2 { Gate::sdg(q) } else { Gate::s(q) };
                        replace[io] = Some(outer.clone());
                        replace[ic] = Some(outer);
                        replace[run[0]] = Some(Gate::h(q));
                        for &r in &run[1..] {
                            delete[r] = stamp;
                        }
                    }
                    _ => unreachable!(),
                }
                changed = true;
                p = j + 1;
            } else if j < track.len() && is_h(&gates[track[j]], q) {
                // Closing H exists but the run carries a T (or a non-π/4 rz):
                // not a Clifford, so leave it — but that H may still open the
                // next triple.
                p = j;
            } else {
                p = j + 1;
            }
        }
    }

    if !changed {
        return None;
    }
    scratch.invalidate_tracks();
    let mut out = scratch.out_buffer();
    let Scratch {
        delete, replace, ..
    } = &mut *scratch;
    out.reserve(gates.len());
    for (i, g) in gates.iter().enumerate() {
        if delete[i] == stamp {
            continue;
        }
        out.push(replace[i].take().unwrap_or_else(|| g.clone()));
    }
    Some(out)
}

#[derive(Clone, Copy)]
enum CommutingPair {
    Cnot(Qubit, Qubit),
    Cz(Qubit, Qubit),
    Ccx(Qubit, Qubit, Qubit),
    Ccz(Qubit, Qubit, Qubit),
}

impl CommutingPair {
    fn from_gate(gate: &Gate) -> Option<Self> {
        match gate {
            Gate::cnot { control, target } => Some(Self::Cnot(*control, *target)),
            Gate::cz { control, target } => {
                let (a, b) = canonical_pair(*control, *target);
                Some(Self::Cz(a, b))
            }
            Gate::ccx {
                control1,
                control2,
                target,
            } => {
                let (a, b, target) = canonical_ccx(*control1, *control2, *target);
                Some(Self::Ccx(a, b, target))
            }
            Gate::ccz {
                control1,
                control2,
                target,
            } => {
                let (a, b, c) = canonical_triple(*control1, *control2, *target);
                Some(Self::Ccz(a, b, c))
            }
            _ => None,
        }
    }

    fn operands(self) -> (usize, [Qubit; 3]) {
        match self {
            Self::Cnot(a, b) | Self::Cz(a, b) => (2, [a, b, 0]),
            Self::Ccx(a, b, c) | Self::Ccz(a, b, c) => (3, [a, b, c]),
        }
    }

    fn matches(self, gate: &Gate) -> bool {
        match self {
            Self::Cnot(c, t) => {
                matches!(gate, Gate::cnot { control, target } if *control == c && *target == t)
            }
            Self::Cz(a, b) => matches!(gate, Gate::cz { control, target }
                if canonical_pair(*control, *target) == (a, b)),
            Self::Ccx(a, b, t) => matches!(gate, Gate::ccx { control1, control2, target }
                if canonical_ccx(*control1, *control2, *target) == (a, b, t)),
            Self::Ccz(a, b, c) => match gate {
                Gate::ccz {
                    control1,
                    control2,
                    target,
                } => canonical_triple(*control1, *control2, *target) == (a, b, c),
                _ => false,
            },
        }
    }

    fn commutes_with(self, gate: &Gate) -> bool {
        match self {
            Self::Cnot(c, t) => commutes_past_cnot(gate, c, t),
            Self::Cz(a, b) => commutes_past_cz(gate, a, b),
            Self::Ccx(a, b, t) => commutes_past_ccx(gate, a, b, t),
            Self::Ccz(a, b, c) => commutes_past_ccz(gate, a, b, c),
        }
    }
}

/// Cancels matching CNOT and CZ pairs across gates that commute with them.
/// Directional CNOT matching, symmetric CZ matching, and their distinct
/// commutation rules share one track build and lookahead walk.
/// Takes the gate list by value, for the reason given on [`reduce_hadamards`].
fn cancel_commuting_pairs(mut gates: Vec<Gate>, scratch: &mut Scratch) -> (Vec<Gate>, bool) {
    let mut changed = false;
    while let Some(next) = cancel_commuting_pairs_pass(&gates, scratch) {
        scratch.spare = std::mem::replace(&mut gates, next);
        changed = true;
    }
    (gates, changed)
}

/// One sweep of two-qubit lookahead cancellation. Per-qubit tracks skip gates
/// that cannot interact with either operand.
fn cancel_commuting_pairs_pass(gates: &[Gate], scratch: &mut Scratch) -> Option<Vec<Gate>> {
    let stamp = scratch.begin();
    scratch.ensure_tracks(gates);
    let Scratch { tracks, delete, .. } = &mut *scratch;
    let mut fired = false;

    for i in 0..gates.len() {
        if delete[i] == stamp {
            continue;
        }
        let pair = match CommutingPair::from_gate(&gates[i]) {
            Some(pair) => pair,
            _ => continue,
        };
        let (arity, operands) = pair.operands();
        let mut cursors = [0usize; 3];
        for operand in 0..arity {
            let track = &tracks[operands[operand] as usize];
            cursors[operand] = track
                .binary_search(&i)
                .expect("gate missing from operand track")
                + 1;
        }
        let mut cancel_at: Option<usize> = None;

        loop {
            let mut next = [None; 3];
            for operand in 0..arity {
                let track = &tracks[operands[operand] as usize];
                while cursors[operand] < track.len() && delete[track[cursors[operand]]] == stamp {
                    cursors[operand] += 1;
                }
                next[operand] = track.get(cursors[operand]).copied();
            }
            let Some(j) = next[..arity].iter().flatten().copied().min() else {
                break;
            };

            if pair.matches(&gates[j]) {
                cancel_at = Some(j);
                break;
            }
            if !pair.commutes_with(&gates[j]) {
                break;
            }
            for operand in 0..arity {
                if next[operand] == Some(j) {
                    cursors[operand] += 1;
                }
            }
        }

        if let Some(j) = cancel_at {
            delete[i] = stamp;
            delete[j] = stamp;
            fired = true;
        }
    }

    if !fired {
        return None;
    }
    scratch.invalidate_tracks();
    let mut out = scratch.out_buffer();
    let delete = &scratch.delete;
    out.reserve(gates.len());
    out.extend(
        gates
            .iter()
            .enumerate()
            .filter(|(i, _)| delete[*i] != stamp)
            .map(|(_, g)| g.clone()),
    );
    Some(out)
}

/// True if `g` commutes past `CNOT(c, t)` so that the CNOT can hop over
/// it without changing the unitary. Covers:
///   - X on the target wire (X(t)·CNOT = CNOT·X(t)).
///   - Diagonal on the control wire (D(c)·CNOT = CNOT·D(c)).
///   - Anything on qubits disjoint from {c, t}.
///   - Another CNOT(c2, t2) provided c2 ≠ t and t2 ≠ c (covers same
///     control / same target / fully-disjoint cases).
///
/// Hadamard, X on control, diagonal on target, CCX touching c or t,
/// measurement, and reset all block.
fn commutes_past_cnot(g: &Gate, c: Qubit, t: Qubit) -> bool {
    match g {
        Gate::x(q) => *q != c,
        Gate::h(q) => *q != c && *q != t,
        Gate::s(q) | Gate::sdg(q) | Gate::z(q) | Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => {
            *q != t
        }
        Gate::cnot {
            control: c2,
            target: t2,
        } => *c2 != t && *t2 != c,
        Gate::cz { control, target } => *control != t && *target != t,
        Gate::ccx {
            control1,
            control2,
            target,
        } => ![*control1, *control2, *target]
            .iter()
            .any(|&q| q == c || q == t),
        Gate::ccz {
            control1,
            control2,
            target,
        } => *target != t && *control1 != t && *control2 != t,
        Gate::measure { qubit, .. } => *qubit != c && *qubit != t,
        Gate::reset(q) => *q != c && *q != t,
    }
}

/// True when `g` commutes with CZ(a, b).
///
/// Diagonal gates and other CZs always commute. Classical controlled gates
/// commute when they do not modify either CZ operand: for CNOT/CCX, only the
/// target is modified. X, H, measurement, and reset block when they act on an
/// operand.
fn commutes_past_cz(g: &Gate, a: Qubit, b: Qubit) -> bool {
    match g {
        Gate::x(q) | Gate::h(q) => *q != a && *q != b,
        Gate::s(_)
        | Gate::sdg(_)
        | Gate::z(_)
        | Gate::t(_)
        | Gate::tdg(_)
        | Gate::rz(..)
        | Gate::cz { .. }
        | Gate::ccz { .. } => true,
        Gate::cnot { target, .. } | Gate::ccx { target, .. } => *target != a && *target != b,
        Gate::measure { qubit, .. } => *qubit != a && *qubit != b,
        Gate::reset(q) => *q != a && *q != b,
    }
}

/// True when `g` commutes with CCX(a, b, t). These rules follow the
/// classical dependency graph: neither operation may modify the other's
/// controls. Diagonal gates commute only on CCX controls, and X commutes on
/// the CCX target.
fn commutes_past_ccx(g: &Gate, a: Qubit, b: Qubit, t: Qubit) -> bool {
    let touched = |q| q == a || q == b || q == t;
    match g {
        Gate::x(q) => *q == t || !touched(*q),
        Gate::h(q) => !touched(*q),
        Gate::s(q) | Gate::sdg(q) | Gate::z(q) | Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => {
            *q != t
        }
        Gate::cnot { control, target } => *target != a && *target != b && *control != t,
        Gate::cz { control, target } => *control != t && *target != t,
        Gate::ccx {
            control1,
            control2,
            target,
        } => *target != a && *target != b && *control1 != t && *control2 != t,
        Gate::ccz {
            control1,
            control2,
            target,
        } => ![*control1, *control2, *target].contains(&t),
        Gate::measure { qubit, .. } => !touched(*qubit),
        Gate::reset(q) => !touched(*q),
    }
}

/// True when `g` commutes with the diagonal CCZ on `{a,b,c}`.
fn commutes_past_ccz(g: &Gate, a: Qubit, b: Qubit, c: Qubit) -> bool {
    let touched = |q| q == a || q == b || q == c;
    match g {
        Gate::x(q) | Gate::h(q) => !touched(*q),
        Gate::s(_)
        | Gate::sdg(_)
        | Gate::z(_)
        | Gate::t(_)
        | Gate::tdg(_)
        | Gate::rz(..)
        | Gate::cz { .. }
        | Gate::ccz { .. } => true,
        Gate::cnot { target, .. } | Gate::ccx { target, .. } => !touched(*target),
        Gate::measure { qubit, .. } => !touched(*qubit),
        Gate::reset(q) => !touched(*q),
    }
}

/// Removes adjacent self-inverse gate pairs (HH, XX, CNOT-CNOT, etc.),
/// commuting gates past non-overlapping operands to expose more pairs.
pub struct CancelGates;

impl Pass for CancelGates {
    fn name(&self) -> &str {
        "Gate cancellation"
    }
    fn run(&self, circuit: &Circuit) -> Circuit {
        let n = circuit.num_qubits;
        // Cancel self-inverse pairs, shrink Hadamard barriers, and cancel
        // CNOT/CZ pairs across commuting gates — alternated to a combined
        // fixpoint. Each step can expose work for the others: dropping
        // gates between two H's or two CNOTs exposes new reducible runs,
        // and a rewrite that emits an X or Z exposes a new cancellable
        // pair.
        // One set of working buffers for every sweep below; see `Scratch`.
        let mut scratch = Scratch::new(n, circuit.gates.len());
        let mut gates = cancel_pairs(&circuit.gates, &mut scratch);
        loop {
            let (reduced, reduce_changed) = reduce_hadamards(gates, &mut scratch);
            let (pair_reduced, pair_changed) = cancel_commuting_pairs(reduced, &mut scratch);
            let before = pair_reduced.len();
            gates = cancel_pairs(&pair_reduced, &mut scratch);
            scratch.spare = pair_reduced;
            let cancel_changed = gates.len() != before;
            if !reduce_changed && !pair_changed && !cancel_changed {
                break;
            }
        }
        let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
        for gate in gates {
            output.apply(gate);
        }
        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unitary::circuits_equiv;

    fn make_circuit(num_qubits: usize, gates: Vec<Gate>) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for g in gates {
            c.apply(g);
        }
        c
    }

    #[test]
    fn hh_cancel() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel() {
        let c = make_circuit(1, vec![Gate::x(0), Gate::x(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel() {
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- HH cancel: varied settings ---

    #[test]
    fn hh_cancel_different_qubit() {
        let c = make_circuit(4, vec![Gate::h(3), Gate::h(3)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_skips_unrelated_gate() {
        // H q0; T q1; H q0 — T on different qubit doesn't block
        let c = make_circuit(2, vec![Gate::h(0), Gate::t(1), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_blocked_by_same_qubit() {
        // H q0; T q0; H q0 — T on same qubit blocks HH cancel
        let c = make_circuit(1, vec![Gate::h(0), Gate::t(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        // Should not cancel as HH; instead rule #1/#2 may or may not apply
        // T is not S or Sdg so no hadamard reduction either — unchanged
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_multiple_pairs() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hh_cancel_parallel_qubits() {
        // H on q0 and q1 independently
        let c = make_circuit(2, vec![Gate::h(0), Gate::h(1), Gate::h(0), Gate::h(1)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- XX cancel: varied settings ---

    #[test]
    fn xx_cancel_different_qubit() {
        let c = make_circuit(3, vec![Gate::x(2), Gate::x(2)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel_skips_unrelated_gate() {
        let c = make_circuit(2, vec![Gate::x(0), Gate::h(1), Gate::x(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::h(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel_blocked_by_same_qubit() {
        let c = make_circuit(1, vec![Gate::x(0), Gate::z(0), Gate::x(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn xx_cancel_multiple_pairs() {
        let c = make_circuit(1, vec![Gate::x(0), Gate::x(0), Gate::x(0), Gate::x(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- CNOT cancel: varied settings ---

    #[test]
    fn cnot_cancel_different_qubits() {
        let c = make_circuit(
            5,
            vec![
                Gate::cnot {
                    control: 3,
                    target: 4,
                },
                Gate::cnot {
                    control: 3,
                    target: 4,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_skips_unrelated_gate() {
        // CNOT q0,q1; T q2; CNOT q0,q1 — T on q2 doesn't interfere
        let c = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::t(2),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(2)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancels_through_diagonal_on_control() {
        // T on the control commutes with CNOT(0,1), so the two CNOTs cancel
        // and only the T survives. Caught by the lookahead pass, not by
        // adjacent pair cancellation.
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::t(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(r.gates[0], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- cancel_cnot (lookahead CNOT cancellation) tests ---

    #[test]
    fn cnots_cancel_through_diagonal_run_on_control() {
        // CNOT · T · S · Tdg · CNOT — all diagonals on q0 (control) commute
        // past the CNOT, so the two CNOTs annihilate and T·S·Tdg survives.
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::t(0),
                Gate::s(0),
                Gate::tdg(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        // Result has no CNOTs.
        assert!(!r.gates.iter().any(|g| matches!(g, Gate::cnot { .. })));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnots_cancel_through_x_on_target() {
        // X on target commutes with CNOT.
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::x(1),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(r.gates[0], Gate::x(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnots_cancel_through_sharing_control_cnot() {
        // CNOT(0,1) · CNOT(0,2) · CNOT(0,1): the middle CNOT shares the
        // control with the outer ones, so they all commute past it. The
        // two outer CNOTs cancel.
        let c = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(
            r.gates[0],
            Gate::cnot {
                control: 0,
                target: 2
            }
        ));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnots_cancel_through_sharing_target_cnot() {
        // CNOT(0,2) · CNOT(1,2) · CNOT(0,2): the middle CNOT shares the
        // target, so the outer ones cancel through it.
        let c = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                Gate::cnot {
                    control: 1,
                    target: 2,
                },
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(
            r.gates[0],
            Gate::cnot {
                control: 1,
                target: 2
            }
        ));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnots_blocked_by_cnot_with_swapped_endpoints() {
        // CNOT(0,1) · CNOT(1,0) · CNOT(0,1) — the middle CNOT has q1 as
        // control and q0 as target, which overlaps with the propagating
        // CNOT's qubits in a way that doesn't commute.
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 1,
                    target: 0,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnots_blocked_by_diagonal_on_target() {
        // T on the target doesn't commute with CNOT.
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::t(1),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnots_blocked_by_x_on_control() {
        // X on control doesn't commute with CNOT (it'd add an X on target).
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::x(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnots_cancel_through_long_mixed_run() {
        // A mix of commuting things in between: Rz on control, X on target,
        // disjoint-qubit gate, sharing-control CNOT. Outer two should cancel.
        let c = make_circuit(
            4,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::t(0),
                Gate::x(1),
                Gate::h(3),
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                Gate::rz(0.31, 0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert!(!r.gates.iter().any(|g| matches!(
            g,
            Gate::cnot {
                control: 0,
                target: 1
            }
        )));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_lookahead_idempotent() {
        // A second run after CancelGates changes nothing.
        let c = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::t(0),
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        let r2 = CancelGates.run(&r);
        assert_eq!(r.gates.len(), r2.gates.len());
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_blocked_by_hadamard_on_control() {
        // H is not diagonal — does not commute past CNOT on either wire.
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_blocked_by_gate_on_target() {
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::t(1),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_no_match_different_direction() {
        // CNOT q0,q1 then CNOT q1,q0 — different direction, should NOT cancel
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 1,
                    target: 0,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_cancel_multiple_pairs() {
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- No-match / preservation tests ---

    #[test]
    fn no_match_preserves_circuit() {
        let c = make_circuit(
            2,
            vec![
                Gate::t(0),
                Gate::s(1),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn single_h_preserved() {
        let c = make_circuit(1, vec![Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::h(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn empty_circuit() {
        let c = Circuit::new(2);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
    }

    // --- ZZ cancel ---

    #[test]
    fn zz_cancel() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::z(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn zz_cancel_skips_unrelated_gate() {
        let c = make_circuit(2, vec![Gate::z(0), Gate::t(1), Gate::z(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn zz_cancel_blocked_by_same_qubit() {
        let c = make_circuit(1, vec![Gate::z(0), Gate::h(0), Gate::z(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- CCX (Toffoli) cancel ---

    #[test]
    fn ccx_cancel() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccz_cancel_is_symmetric_in_all_operands() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccz {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccz {
                    control1: 2,
                    control2: 0,
                    target: 1,
                },
            ],
        );

        let r = CancelGates.run(&c);

        assert!(r.gates.is_empty());
        assert!(!r.has_toffoli());
        assert!(!r.has_ccz());
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccz_cancellation_is_blocked_on_each_operand() {
        for blocker in 0..3 {
            let c = make_circuit(
                3,
                vec![
                    Gate::ccz {
                        control1: 0,
                        control2: 1,
                        target: 2,
                    },
                    Gate::h(blocker),
                    Gate::ccz {
                        control1: 2,
                        control2: 0,
                        target: 1,
                    },
                ],
            );

            let r = CancelGates.run(&c);

            assert_eq!(r.gates.len(), 3, "blocker q{blocker}");
            assert!(r.has_ccz(), "blocker q{blocker}");
            assert!(circuits_equiv(&c, &r, 1e-10), "blocker q{blocker}");
        }
    }

    #[test]
    fn ccz_cancels_across_disjoint_gate() {
        let c = make_circuit(
            4,
            vec![
                Gate::ccz {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::t(3),
                Gate::ccz {
                    control1: 2,
                    control2: 1,
                    target: 0,
                },
            ],
        );

        let r = CancelGates.run(&c);

        assert!(matches!(r.gates.as_slice(), [Gate::t(3)]));
        assert!(!r.has_ccz());
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_skips_unrelated_gate() {
        let c = make_circuit(
            4,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::t(3),
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::t(3)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_blocked_by_control1() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::h(0),
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_blocked_by_control2() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::h(1),
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_blocked_by_target() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::h(2),
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_no_cancel_different_controls() {
        let c = make_circuit(
            4,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 3,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 2,
                    target: 3,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_no_cancel_different_target() {
        let c = make_circuit(
            4,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 3,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_cancel_multiple_pairs() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Cascading cancellation ---

    #[test]
    fn cascade_nested_h() {
        // H H H H — inner pair cancels, exposing outer pair
        let c = make_circuit(1, vec![Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_six_h() {
        let c = make_circuit(
            1,
            vec![
                Gate::h(0),
                Gate::h(0),
                Gate::h(0),
                Gate::h(0),
                Gate::h(0),
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_odd_count_leaves_one() {
        let c = make_circuit(
            1,
            vec![Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0), Gate::h(0)],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_nested_cnot() {
        // CNOT CNOT CNOT CNOT — fully cancels in one pass
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_mixed_self_inverse() {
        // H X H X — no adjacent self-inverse pair cancels, but the Hadamard
        // pass rewrites H·X·H = Z, leaving Z·X.
        let c = make_circuit(1, vec![Gate::h(0), Gate::x(0), Gate::h(0), Gate::x(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(r.gates[0], Gate::z(0)));
        assert!(matches!(r.gates[1], Gate::x(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_different_qubits_interleaved() {
        // H(q0) H(q1) H(q1) H(q0) — inner H(q1) pair cancels, then outer H(q0) cancels
        let c = make_circuit(2, vec![Gate::h(0), Gate::h(1), Gate::h(1), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_deep_nesting() {
        // H(0) X(0) Z(0) Z(0) X(0) H(0) — Z pair cancels, exposes X pair, exposes H pair
        let c = make_circuit(
            1,
            vec![
                Gate::h(0),
                Gate::x(0),
                Gate::z(0),
                Gate::z(0),
                Gate::x(0),
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_deep_nesting_with_residual() {
        // T(0) H(0) X(0) X(0) H(0) T(0) — X cancels, H cancels, T is not self-inverse
        let c = make_circuit(
            1,
            vec![
                Gate::t(0),
                Gate::h(0),
                Gate::x(0),
                Gate::x(0),
                Gate::h(0),
                Gate::t(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cascade_cnot_nested_in_h() {
        // H(0) CNOT(0,1) CNOT(0,1) H(0) — CNOT pair cancels, H pair cancels
        let c = make_circuit(
            2,
            vec![
                Gate::h(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Commutation through unrelated qubits ---

    #[test]
    fn commute_h_past_many_unrelated() {
        // H(0); T(1); S(2); Tdg(3); H(0) — all middle gates on different qubits
        let c = make_circuit(
            4,
            vec![Gate::h(0), Gate::t(1), Gate::s(2), Gate::tdg(3), Gate::h(0)],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn commute_cnot_past_unrelated_qubits() {
        // CNOT(0,1); H(2); T(3); CNOT(0,1)
        let c = make_circuit(
            4,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(2),
                Gate::t(3),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn commute_ccx_past_unrelated_qubits() {
        // CCX(0,1,2); H(3); T(4); CCX(0,1,2)
        let c = make_circuit(
            5,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::h(3),
                Gate::t(4),
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_blocked_by_gate_on_either_qubit() {
        // Gate on control blocks
        let c1 = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r1 = CancelGates.run(&c1);
        assert_eq!(r1.gates.len(), 3);
        assert!(circuits_equiv(&c1, &r1, 1e-10));

        // Gate on target blocks
        let c2 = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(1),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r2 = CancelGates.run(&c2);
        assert_eq!(r2.gates.len(), 3);
        assert!(circuits_equiv(&c2, &r2, 1e-10));

        // Gate on unrelated qubit doesn't block
        let c3 = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(2),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r3 = CancelGates.run(&c3);
        assert_eq!(r3.gates.len(), 1);
        assert!(circuits_equiv(&c3, &r3, 1e-10));
    }

    // --- Non-self-inverse gates don't cancel ---

    #[test]
    fn t_t_no_cancel() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::t(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn tdg_tdg_no_cancel() {
        let c = make_circuit(1, vec![Gate::tdg(0), Gate::tdg(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_s_no_cancel() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::s(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn sdg_sdg_no_cancel() {
        let c = make_circuit(1, vec![Gate::sdg(0), Gate::sdg(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn t_tdg_no_cancel() {
        // T and Tdg are inverses of each other but cancel_pairs only handles self-inverse
        let c = make_circuit(1, vec![Gate::t(0), Gate::tdg(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn s_sdg_no_cancel() {
        let c = make_circuit(1, vec![Gate::s(0), Gate::sdg(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn rz_rz_no_cancel() {
        let c = make_circuit(
            1,
            vec![
                Gate::rz(std::f64::consts::PI / 4.0, 0),
                Gate::rz(std::f64::consts::PI / 4.0, 0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Multi-qubit independence ---

    #[test]
    fn independent_cancellations_on_many_qubits() {
        // Each qubit has its own H-H pair, all should cancel independently
        let c = make_circuit(
            5,
            vec![
                Gate::h(0),
                Gate::h(1),
                Gate::h(2),
                Gate::h(3),
                Gate::h(4),
                Gate::h(4),
                Gate::h(3),
                Gate::h(2),
                Gate::h(1),
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn partial_cancellation_mixed_qubits() {
        // q0: H H (cancels), q1: H T H (blocked), q2: X X (cancels)
        let c = make_circuit(
            3,
            vec![
                Gate::h(0),
                Gate::h(1),
                Gate::x(2),
                Gate::t(1),
                Gate::h(0),
                Gate::h(1),
                Gate::x(2),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3); // H(1), T(1), H(1) remain
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn interleaved_cancel_different_gate_types() {
        // H(0) X(1) H(0) X(1) — both pairs cancel through each other
        let c = make_circuit(2, vec![Gate::h(0), Gate::x(1), Gate::h(0), Gate::x(1)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- CNOT blocking edge cases ---

    #[test]
    fn cnot_blocks_h_on_shared_qubit() {
        // H(0); CNOT(0,1); H(0) — CNOT touches q0 so it blocks
        let c = make_circuit(
            2,
            vec![
                Gate::h(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn h_not_blocked_by_cnot_on_other_qubits() {
        // H(0); CNOT(1,2); H(0) — CNOT doesn't touch q0
        let c = make_circuit(
            3,
            vec![
                Gate::h(0),
                Gate::cnot {
                    control: 1,
                    target: 2,
                },
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn ccx_blocks_h_on_any_of_three_qubits() {
        for blocked_q in 0..3 {
            let c = make_circuit(
                4,
                vec![
                    Gate::h(blocked_q),
                    Gate::ccx {
                        control1: 0,
                        control2: 1,
                        target: 2,
                    },
                    Gate::h(blocked_q),
                ],
            );
            let r = CancelGates.run(&c);
            assert_eq!(
                r.gates.len(),
                3,
                "H({}) should be blocked by CCX",
                blocked_q
            );
            assert!(circuits_equiv(&c, &r, 1e-10));
        }
        // q3 is unrelated — should cancel
        let c = make_circuit(
            4,
            vec![
                Gate::h(3),
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::h(3),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Order preservation ---

    #[test]
    fn surviving_gates_preserve_order() {
        let c = make_circuit(
            3,
            vec![
                Gate::t(0),
                Gate::h(1),
                Gate::h(1), // cancels
                Gate::s(0),
                Gate::x(2),
                Gate::x(2), // cancels
                Gate::tdg(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::s(0)));
        assert!(matches!(&r.gates[2], Gate::tdg(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn order_after_cascade() {
        // T(0) H(0) H(0) S(0) — H pair cancels, T and S remain in order
        let c = make_circuit(1, vec![Gate::t(0), Gate::h(0), Gate::h(0), Gate::s(0)]);
        let r = CancelGates.run(&c);
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
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn many_pairs_odd_leaves_one() {
        let gates: Vec<Gate> = (0..201).map(|_| Gate::h(0)).collect();
        let c = make_circuit(1, gates);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn alternating_cancel_no_cancel() {
        // H(0) H(0) T(0) H(0) H(0) T(0) — two H pairs cancel around T's
        let c = make_circuit(
            1,
            vec![
                Gate::h(0),
                Gate::h(0),
                Gate::t(0),
                Gate::h(0),
                Gate::h(0),
                Gate::t(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::t(0)));
        assert!(matches!(&r.gates[1], Gate::t(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn deeply_nested_cascade_8_layers() {
        // H X Z H H Z X H — 4 nested pairs, all cancel
        let c = make_circuit(
            1,
            vec![
                Gate::h(0),
                Gate::x(0),
                Gate::z(0),
                Gate::h(0),
                Gate::h(0),
                Gate::z(0),
                Gate::x(0),
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- has_toffoli flag ---

    #[test]
    fn has_toffoli_set_when_ccx_survives() {
        let c = make_circuit(
            4,
            vec![
                Gate::h(3),
                Gate::h(3), // cancels
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(r.has_toffoli());
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn has_toffoli_cleared_when_ccx_cancelled() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(!r.has_toffoli());
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Mixed multi-qubit gate interactions ---

    #[test]
    fn cnot_and_ccx_block_each_other_on_shared_qubit() {
        // CNOT(0,1); CCX(1,2,3); CNOT(0,1) — CCX touches q1 which blocks CNOT
        let c = make_circuit(
            4,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::ccx {
                    control1: 1,
                    control2: 2,
                    target: 3,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_and_ccx_no_shared_qubit() {
        // CNOT(0,1); CCX(2,3,4); CNOT(0,1) — no shared qubits, CNOT cancels
        let c = make_circuit(
            5,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::ccx {
                    control1: 2,
                    control2: 3,
                    target: 4,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn multiple_cnot_pairs_different_qubit_pairs() {
        // CNOT(0,1) CNOT(2,3) CNOT(2,3) CNOT(0,1)
        let c = make_circuit(
            4,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 2,
                    target: 3,
                },
                Gate::cnot {
                    control: 2,
                    target: 3,
                },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Edge: single gate of each self-inverse type (no cancel) ---

    #[test]
    fn single_gate_each_type_preserved() {
        let c = make_circuit(
            5,
            vec![
                Gate::h(0),
                Gate::x(1),
                Gate::z(2),
                Gate::cnot {
                    control: 3,
                    target: 4,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Blocker is latest gate, not first ---

    #[test]
    fn blocker_is_latest_not_earliest() {
        // H(0); T(0); X(0); H(0) — blocker for second H is X (latest on q0), not T
        let c = make_circuit(1, vec![Gate::h(0), Gate::t(0), Gate::x(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        // X is the blocker, not H — X != H so no cancel
        assert_eq!(r.gates.len(), 4);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn blocker_matches_when_latest_is_same() {
        // H(0); T(0); H(0); H(0) — second H is blocked by first H? No, T is between.
        // Actually: stack is [H@0, T@1, H@2]. Third H@3 checks blocker = H@2 (latest on q0).
        // H@2 == H@3, so they cancel. Leaves [H@0, T@1].
        let c = make_circuit(1, vec![Gate::h(0), Gate::t(0), Gate::h(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
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
        let c = make_circuit(
            2,
            vec![
                Gate::t(1),
                Gate::h(0),
                Gate::x(1),
                Gate::h(0),
                Gate::x(1),
                Gate::t(1),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
        assert!(matches!(&r.gates[0], Gate::t(1)));
        assert!(matches!(&r.gates[1], Gate::t(1)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- Large cascade depth ---

    // --- measurement / reset barrier tests ---

    fn make_circuit_with_cbits(num_qubits: usize, num_cbits: usize, gates: Vec<Gate>) -> Circuit {
        let mut c = Circuit::with_cbits(num_qubits, num_cbits);
        for g in gates {
            c.apply(g);
        }
        c
    }

    #[test]
    fn hh_blocked_by_measure_on_same_qubit() {
        let c = make_circuit_with_cbits(
            1,
            1,
            vec![Gate::h(0), Gate::measure { qubit: 0, cbit: 0 }, Gate::h(0)],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(r.has_measurement());
        assert_eq!(r.num_cbits, 1);
    }

    #[test]
    fn hh_blocked_by_reset_on_same_qubit() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::reset(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(r.has_measurement());
    }

    #[test]
    fn hh_allowed_past_measure_on_other_qubit() {
        let c = make_circuit_with_cbits(
            2,
            1,
            vec![Gate::h(0), Gate::measure { qubit: 1, cbit: 0 }, Gate::h(0)],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::measure { qubit: 1, cbit: 0 }));
    }

    #[test]
    fn hh_allowed_past_reset_on_other_qubit() {
        let c = make_circuit(2, vec![Gate::h(0), Gate::reset(1), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::reset(1)));
    }

    #[test]
    fn xx_blocked_by_measure() {
        let c = make_circuit_with_cbits(
            1,
            1,
            vec![Gate::x(0), Gate::measure { qubit: 0, cbit: 0 }, Gate::x(0)],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    #[test]
    fn cnot_blocked_by_measure_on_control() {
        let c = make_circuit_with_cbits(
            2,
            1,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::measure { qubit: 0, cbit: 0 },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    #[test]
    fn cnot_blocked_by_measure_on_target() {
        let c = make_circuit_with_cbits(
            2,
            1,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::measure { qubit: 1, cbit: 0 },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    #[test]
    fn cnot_allowed_past_measure_on_disjoint_qubit() {
        let c = make_circuit_with_cbits(
            3,
            1,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::measure { qubit: 2, cbit: 0 },
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(&r.gates[0], Gate::measure { qubit: 2, cbit: 0 }));
    }

    #[test]
    fn cnot_blocked_by_reset_on_target() {
        let c = make_circuit(
            2,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::reset(1),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
    }

    #[test]
    fn ccx_blocked_by_measure_on_any_qubit() {
        for blocked_q in 0..3 {
            let c = make_circuit_with_cbits(
                3,
                1,
                vec![
                    Gate::ccx {
                        control1: 0,
                        control2: 1,
                        target: 2,
                    },
                    Gate::measure {
                        qubit: blocked_q,
                        cbit: 0,
                    },
                    Gate::ccx {
                        control1: 0,
                        control2: 1,
                        target: 2,
                    },
                ],
            );
            let r = CancelGates.run(&c);
            assert_eq!(
                r.gates.len(),
                3,
                "measure q{} should block CCX pair",
                blocked_q
            );
        }
    }

    #[test]
    fn measure_reset_pair_does_not_cancel() {
        // measure and reset are NOT self-inverse, so two adjacent ones must stay.
        let c = make_circuit_with_cbits(
            1,
            1,
            vec![
                Gate::measure { qubit: 0, cbit: 0 },
                Gate::measure { qubit: 0, cbit: 0 },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);

        let c = make_circuit(1, vec![Gate::reset(0), Gate::reset(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 2);
    }

    #[test]
    fn num_cbits_preserved() {
        let c = make_circuit_with_cbits(1, 3, vec![Gate::h(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert_eq!(r.num_cbits, 3);
    }

    #[test]
    fn cascade_depth_10() {
        // 10 layers of nesting: H X Z H X Z X H Z X H Z X H ... all cancel
        // Simpler: alternating gate types nested symmetrically
        // Use q0: H(X(Z(Z(X(H()))))) = H X Z Z X H
        let c = make_circuit(
            1,
            vec![
                Gate::h(0),
                Gate::x(0),
                Gate::z(0),
                Gate::z(0),
                Gate::x(0),
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));

        // 5 layers deep
        let c2 = make_circuit(
            1,
            vec![
                Gate::h(0),
                Gate::x(0),
                Gate::z(0),
                Gate::x(0),
                Gate::h(0),
                Gate::h(0),
                Gate::x(0),
                Gate::z(0),
                Gate::x(0),
                Gate::h(0),
            ],
        );
        let r2 = CancelGates.run(&c2);
        assert_eq!(r2.gates.len(), 0);
        assert!(circuits_equiv(&c2, &r2, 1e-10));
    }

    // --- Hadamard-reduction tests ---
    // CancelGates also collapses H barriers via local Clifford identities, so
    // the downstream phase-folding pass sees longer Hadamard-free sections.

    fn count_h(c: &Circuit) -> usize {
        c.gates.iter().filter(|g| matches!(g, Gate::h(_))).count()
    }

    #[test]
    fn hxh_becomes_z() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::x(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(r.gates[0], Gate::z(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hzh_becomes_x() {
        let c = make_circuit(1, vec![Gate::h(0), Gate::z(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 1);
        assert!(matches!(r.gates[0], Gate::x(0)));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hsh_drops_one_hadamard() {
        // H·S·H = Sdg·H·Sdg — still three gates, but one fewer Hadamard.
        let c = make_circuit(1, vec![Gate::h(0), Gate::s(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(count_h(&r), 1);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hsh_cascade_eliminates_all_hadamards() {
        // H·S·H·S·H reduces to a single Sdg, with no Hadamards left.
        let c = make_circuit(
            1,
            vec![Gate::h(0), Gate::s(0), Gate::h(0), Gate::s(0), Gate::h(0)],
        );
        let r = CancelGates.run(&c);
        assert_eq!(count_h(&r), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hadamard_run_summing_to_identity_collapses() {
        // H·Sdg·S·H — the run Sdg·S is the identity, so both H's vanish.
        let c = make_circuit(1, vec![Gate::h(0), Gate::sdg(0), Gate::s(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hxh_reduced_despite_other_qubit_gate_between() {
        // A gate on another wire interleaved between the H's must not block
        // the (wire-local) H·X·H = Z rewrite.
        let c = make_circuit(
            3,
            vec![
                Gate::h(0),
                Gate::x(0),
                Gate::cnot {
                    control: 1,
                    target: 2,
                },
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(count_h(&r), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_inside_hadamards_blocks_reduction() {
        // A CNOT touching the wire between the two H's makes the run
        // non-Clifford — nothing may be rewritten.
        let c = make_circuit(
            2,
            vec![
                Gate::h(0),
                Gate::cnot {
                    control: 0,
                    target: 1,
                },
                Gate::h(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn hadamard_run_with_t_is_not_reducible() {
        // H·T·H carries a genuine T — not a Clifford run, both H's survive.
        let c = make_circuit(1, vec![Gate::h(0), Gate::t(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(count_h(&r), 2);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn lone_hadamard_is_kept() {
        let c = make_circuit(1, vec![Gate::t(0), Gate::h(0), Gate::t(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn measure_between_hadamards_blocks_reduction() {
        let mut c = Circuit::with_cbits(1, 1);
        c.apply(Gate::h(0));
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::h(0));
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(matches!(r.gates[1], Gate::measure { .. }));
    }

    #[test]
    fn reduction_exposes_pair_cancellation() {
        // H·X·H = Z, and the emitted Z then cancels with the leading Z — only
        // the combined cancel+reduce fixpoint catches this.
        let c = make_circuit(1, vec![Gate::z(0), Gate::h(0), Gate::x(0), Gate::h(0)]);
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn pair_cancellation_exposes_reduction() {
        // The inner X·X cancels, leaving H·H which then also cancels — the
        // pair pass feeds the Hadamard pass.
        let c = make_circuit(
            1,
            vec![
                Gate::t(0),
                Gate::h(0),
                Gate::x(0),
                Gate::x(0),
                Gate::h(0),
                Gate::t(0),
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(count_h(&r), 0);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    // --- CZ cancellation ---

    #[test]
    fn adjacent_cz_pair_cancels() {
        let c = make_circuit(
            2,
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::cz {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert!(r.gates.is_empty());
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn reversed_operand_cz_pair_cancels() {
        let c = make_circuit(
            3,
            vec![
                Gate::cz {
                    control: 0,
                    target: 2,
                },
                Gate::cz {
                    control: 2,
                    target: 0,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert!(r.gates.is_empty());
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cz_pair_cancels_through_diagonals_on_both_operands() {
        let c = make_circuit(
            2,
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::t(0),
                Gate::sdg(1),
                Gate::rz(0.31, 0),
                Gate::z(1),
                Gate::cz {
                    control: 1,
                    target: 0,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 4);
        assert!(!r.gates.iter().any(|g| matches!(g, Gate::cz { .. })));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cz_pair_lookahead_commutation_table() {
        let cases = [
            (
                "other CZ",
                3,
                Gate::cz {
                    control: 1,
                    target: 2,
                },
                true,
                1,
            ),
            (
                "CNOT with external target",
                3,
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                true,
                0,
            ),
            (
                "CNOT targeting operand",
                3,
                Gate::cnot {
                    control: 2,
                    target: 1,
                },
                false,
                2,
            ),
            (
                "CCX with external target",
                4,
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 3,
                },
                true,
                0,
            ),
            (
                "CCX targeting operand",
                3,
                Gate::ccx {
                    control1: 1,
                    control2: 2,
                    target: 0,
                },
                false,
                2,
            ),
            ("X on operand", 2, Gate::x(0), false, 2),
            ("H on operand", 2, Gate::h(1), false, 2),
            ("X on other qubit", 3, Gate::x(2), true, 0),
        ];

        for (name, num_qubits, middle, should_cancel, expected_cz_count) in cases {
            let c = make_circuit(
                num_qubits,
                vec![
                    Gate::cz {
                        control: 0,
                        target: 1,
                    },
                    middle,
                    Gate::cz {
                        control: 1,
                        target: 0,
                    },
                ],
            );
            let r = CancelGates.run(&c);
            let cz_count = r
                .gates
                .iter()
                .filter(|g| matches!(g, Gate::cz { .. }))
                .count();
            assert_eq!(cz_count, expected_cz_count, "{name}");
            assert_eq!(r.gates.len(), if should_cancel { 1 } else { 3 }, "{name}");
            assert!(circuits_equiv(&c, &r, 1e-10), "{name}");
        }
    }

    #[test]
    fn measurement_on_operand_blocks_cz_cancellation() {
        let mut c = Circuit::with_cbits(2, 1);
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::cz {
            control: 1,
            target: 0,
        });
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(r.has_measurement());
    }

    #[test]
    fn reset_on_operand_blocks_cz_cancellation() {
        let c = make_circuit(
            2,
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::reset(1),
                Gate::cz {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(r.has_measurement());
    }

    #[test]
    fn cnot_pair_cancels_through_commuting_cz() {
        let c = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert!(matches!(
            r.gates.as_slice(),
            [Gate::cz {
                control: 0,
                target: 1
            }]
        ));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn cnot_pair_blocked_by_cz_on_target() {
        let c = make_circuit(
            3,
            vec![
                Gate::cnot {
                    control: 2,
                    target: 1,
                },
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::cnot {
                    control: 2,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(r.gates.len(), 3);
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn alternating_cz_pairs_cancel_to_fixpoint() {
        let c = make_circuit(
            3,
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::cz {
                    control: 1,
                    target: 2,
                },
                Gate::t(0),
                Gate::cz {
                    control: 1,
                    target: 0,
                },
                Gate::cz {
                    control: 2,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert!(matches!(r.gates.as_slice(), [Gate::t(0)]));
        assert!(circuits_equiv(&c, &r, 1e-10));
    }

    #[test]
    fn measurement_on_other_qubit_does_not_block_cz_cancellation() {
        let mut c = Circuit::with_cbits(3, 1);
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::measure { qubit: 2, cbit: 0 });
        c.apply(Gate::cz {
            control: 1,
            target: 0,
        });
        let r = CancelGates.run(&c);
        assert!(matches!(
            r.gates.as_slice(),
            [Gate::measure { qubit: 2, cbit: 0 }]
        ));
        assert!(r.has_measurement());
        assert_eq!(r.num_cbits, 1);
    }

    #[test]
    fn reset_on_other_qubit_does_not_block_cz_cancellation() {
        let c = make_circuit(
            3,
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::reset(2),
                Gate::cz {
                    control: 0,
                    target: 1,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert!(matches!(r.gates.as_slice(), [Gate::reset(2)]));
        assert!(r.has_measurement());
    }

    #[test]
    fn cz_cancellation_through_long_mixed_commuting_run() {
        let c = make_circuit(
            5,
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::t(0),
                Gate::rz(0.19, 1),
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 3,
                },
                Gate::cz {
                    control: 2,
                    target: 4,
                },
                Gate::x(4),
                Gate::cz {
                    control: 1,
                    target: 0,
                },
            ],
        );
        let r = CancelGates.run(&c);
        assert_eq!(
            r.gates
                .iter()
                .filter(|g| matches!(
                    g,
                    Gate::cz {
                        control: 0,
                        target: 1
                    } | Gate::cz {
                        control: 1,
                        target: 0
                    }
                ))
                .count(),
            0
        );
        assert_eq!(
            r.gates
                .iter()
                .filter(|g| matches!(
                    g,
                    Gate::cz {
                        control: 2,
                        target: 4
                    }
                ))
                .count(),
            1
        );
        assert!(circuits_equiv(&c, &r, 1e-9));
    }

    #[test]
    fn native_cz_cancellation_is_structurally_idempotent() {
        let c = make_circuit(
            4,
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                },
                Gate::t(0),
                Gate::cz {
                    control: 2,
                    target: 3,
                },
                Gate::cz {
                    control: 1,
                    target: 0,
                },
                Gate::cz {
                    control: 3,
                    target: 2,
                },
            ],
        );
        let once = CancelGates.run(&c);
        let twice = CancelGates.run(&once);
        assert_eq!(once.to_qasm(), twice.to_qasm());
        assert!(circuits_equiv(&c, &twice, 1e-10));
    }

    #[test]
    fn ccx_swapped_controls_cancel_across_commuting_control_phase() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::t(0),
                Gate::ccx {
                    control1: 1,
                    control2: 0,
                    target: 2,
                },
            ],
        );
        let result = CancelGates.run(&c);
        assert_eq!(result.gates, vec![Gate::t(0)]);
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    #[test]
    fn ccx_lookahead_is_blocked_by_a_target_phase() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::t(2),
                Gate::ccx {
                    control1: 1,
                    control2: 0,
                    target: 2,
                },
            ],
        );
        assert_eq!(CancelGates.run(&c).gates.len(), 3);
    }

    #[test]
    fn ccz_lookahead_cancels_across_diagonal_gates() {
        let c = make_circuit(
            3,
            vec![
                Gate::ccz {
                    control1: 0,
                    control2: 1,
                    target: 2,
                },
                Gate::rz(0.37, 1),
                Gate::ccz {
                    control1: 2,
                    control2: 0,
                    target: 1,
                },
            ],
        );
        let result = CancelGates.run(&c);
        assert_eq!(result.gates, vec![Gate::rz(0.37, 1)]);
        assert!(circuits_equiv(&c, &result, 1e-10));
    }

    #[test]
    fn every_positive_ccx_and_ccz_commutation_rule_is_unitary_sound() {
        let mut candidates = Vec::new();
        for q in 0..4 {
            candidates.extend([
                Gate::x(q),
                Gate::h(q),
                Gate::s(q),
                Gate::sdg(q),
                Gate::z(q),
                Gate::t(q),
                Gate::tdg(q),
                Gate::rz(0.37, q),
            ]);
        }
        for a in 0..4 {
            for b in 0..4 {
                if a == b {
                    continue;
                }
                candidates.push(Gate::cnot {
                    control: a,
                    target: b,
                });
                candidates.push(Gate::cz {
                    control: a,
                    target: b,
                });
                for c in 0..4 {
                    if c == a || c == b {
                        continue;
                    }
                    candidates.push(Gate::ccx {
                        control1: a,
                        control2: b,
                        target: c,
                    });
                    candidates.push(Gate::ccz {
                        control1: a,
                        control2: b,
                        target: c,
                    });
                }
            }
        }

        let ccx = Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        };
        let ccz = Gate::ccz {
            control1: 0,
            control2: 1,
            target: 2,
        };
        for gate in candidates {
            for (native, commutes) in [
                (&ccx, commutes_past_ccx(&gate, 0, 1, 2)),
                (&ccz, commutes_past_ccz(&gate, 0, 1, 2)),
            ] {
                if !commutes {
                    continue;
                }
                let left = make_circuit(4, vec![native.clone(), gate.clone()]);
                let right = make_circuit(4, vec![gate.clone(), native.clone()]);
                assert!(
                    circuits_equiv(&left, &right, 1e-10),
                    "unsound commutation: {native:?} across {gate:?}"
                );
            }
        }
    }
}
