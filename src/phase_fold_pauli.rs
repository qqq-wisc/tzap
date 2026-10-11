//! Circuit-level Pauli folding with exact commutation checks.
//!
//! The pass walks the circuit once and keeps the Clifford frame: for the
//! Clifford prefix `U` seen so far, the row of each generator `P` (`X_q` or
//! `Z_q`) is `U† P U`. A rotation about `P` after `U` equals a rotation about
//! `U† P U` before it. RX uses the X row, RZ/P/T/T† the Z row, and RY the
//! signed Hermitian product `i X Z`. Two rotations about the same unsigned
//! axis merge, at the later one's position, when every live rotation recorded
//! between them commutes with that axis. CCX and CCZ gates, measurements and
//! resets are recorded as *blockers*: they never merge, but a fold across one
//! must commute with it.
//!
//! The scan is made fast by these tricks. `docs/tableau-abstraction.tex`,
//! section "Making the scan fast", has the details and measurements.
//!
//! - **Fingerprints.** Each X and Z generator gets a random 128-bit label. A
//!   row's fingerprint is the XOR of the labels of its factors, so a second,
//!   *sketch* frame keeps every row's fingerprint at O(1) per gate. A map from
//!   fingerprint to the newest rotation, chained to the older ones, proposes
//!   merge candidates. Only exact axes and an exact commutation check
//!   authorize a rewrite, so the randomness affects speed, not output.
//! - **Pre-check.** A sketch-only scan looks for a repeated fingerprint. If
//!   there is none, nothing can fold, and the exact frame is never built.
//! - **Support lists.** Two axes can anticommute only if they share a qubit.
//!   Each event is listed under each bit of its 64-bit support signature, so
//!   the commutation check visits only events that share support with the
//!   axis.
//! - **Compaction.** Most rotations fold away. A support list drops its dead
//!   entries once a fixed fraction of it is dead, or later scans would walk
//!   past them.
//! - **Bit slicing.** Up to [`MAX_SLICED_QUBITS`] qubits, each support list
//!   is stored as per-qubit X and Z masks over blocks of 64 entries, and a
//!   whole block is tested against an axis with a few word XORs.
//! - **Budget.** Each fold attempt may do a fixed amount of work. Running out
//!   counts as blocked, which bounds the pass to linear time.

use std::collections::hash_map::RandomState;
#[cfg(test)]
use std::f64::consts::PI;
use std::hash::{BuildHasher, Hasher};

use rustc_hash::{FxHashMap, FxHashSet};

use crate::circuit::{Circuit, Gate, qubit_operands};
use crate::pass::Pass;
use crate::pbc::{packed_anticommutes, support_signature};

/// Maximum retained exact axis storage, in 64-bit words (256 MiB). Reaching
/// it stops further folding; folds already found are kept. Each event stores
/// at least two words, so this also keeps event ids below 2^24, small enough
/// for the `u32` ids of the support lists.
const MAX_AXIS_WORDS: usize = 1 << 25;

/// Work allowed for one fold attempt (candidate lookup plus the commutation
/// check). An attempt that runs out is treated as blocked, and the scan goes
/// on, so one expensive candidate cannot end the pass.
const MAX_ATTEMPT_STEPS: usize = 1 << 20;

/// Widest circuit that uses bit-sliced support lists. A rotation with support
/// `w` costs `O(w²)` to record there (it sets `w` bits in each of `w` lists),
/// which dense wide circuits such as QFT do not recover in faster scans.
const MAX_SLICED_QUBITS: usize = 32;

/// Folds T/T†, P and RZ/RX/RY across Clifford gates, including Hadamards.
/// Numeric angles never snap to Clifford angles or undergo floating modulo
/// reduction. Checked angle arithmetic retains exact pi fractions and declines
/// numeric merges that would round, overflow, or exceed representation limits.
/// Can run independently or after [`crate::phase_fold_rand::PhaseFoldRand`].
pub struct PhaseFoldPauli;

impl Pass for PhaseFoldPauli {
    fn name(&self) -> &str {
        "Pauli folding"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        phase_fold_pauli(circuit)
    }
}

/// Fold rotations using checked angle arithmetic and exact Pauli commutation.
/// Random fingerprints only select candidates; exact checks authorize rewrites.
pub fn phase_fold_pauli(circuit: &Circuit) -> Circuit {
    let lowered = crate::angle::lowered_if_needed(circuit);
    let circuit = &lowered;
    let n = circuit.num_qubits;
    if !frame_fits(n) || !well_formed(circuit) {
        return (**circuit).clone();
    }
    let labels: Vec<_> = (0..circuit.num_qubits)
        .map(|_| (random_label(), random_label()))
        .collect();
    let sliced = circuit.num_qubits <= MAX_SLICED_QUBITS;
    fold(circuit, &labels, MAX_ATTEMPT_STEPS, sliced)
}

fn frame_fits(n: usize) -> bool {
    let l = n.div_ceil(64).max(1);
    n.checked_mul(4 * l)
        .is_some_and(|words| words <= MAX_AXIS_WORDS)
}

fn random_label() -> u128 {
    let hi = RandomState::new().build_hasher().finish() as u128;
    let lo = RandomState::new().build_hasher().finish() as u128;
    (hi << 64) | lo
}

/// The pass, with fingerprint labels `(X_q, Z_q)` for each qubit `q`, a
/// per-attempt budget, and bit-sliced (`sliced`) or scalar support lists. The
/// two kinds of list reach the same verdicts; only the budget accounting
/// differs.
fn fold(circuit: &Circuit, labels: &[(u128, u128)], attempt_steps: usize, sliced: bool) -> Circuit {
    let n = circuit.num_qubits;
    let l = n.div_ceil(64).max(1);
    // The exact frame has 2n rows of 2l words. Bound it before allocating it,
    // even for a circuit with very few rotations.
    if labels.len() != n || !frame_fits(n) || !well_formed(circuit) || !may_fold(circuit, labels) {
        return circuit.clone();
    }
    let mut folder = Folder::new(labels, n, l, attempt_steps, sliced, circuit.gates.len());
    for (idx, gate) in circuit.gates.iter().enumerate() {
        if folder.step(idx, gate).is_none() {
            break; // axis storage is full
        }
    }
    folder.rewrite(circuit)
}

/// Every gate kind is handled; the input must only be well formed (qubits
/// in range and distinct, finite angles). Anything else is left as is.
fn well_formed(circuit: &Circuit) -> bool {
    let n = circuit.num_qubits;
    circuit.gates.iter().all(|g| {
        let (len, mut qubits) = qubit_operands(g);
        let qubits = &mut qubits[..len];
        qubits.sort_unstable();
        qubits.iter().all(|&q| (q as usize) < n) && qubits.windows(2).all(|w| w[0] != w[1])
    })
}

/// The pre-check: whether two non-Clifford rotations of the unfolded circuit
/// have the same axis fingerprint. If none do, no two have the same axis.
/// The first fold needs such a pair, so then nothing folds.
fn may_fold(circuit: &Circuit, labels: &[(u128, u128)]) -> bool {
    let mut rotations = circuit
        .gates
        .iter()
        .filter(|g| RotationKind::of(g).is_some());
    if rotations.nth(1).is_none() {
        return false;
    }
    let mut sketch = SketchFrame::new(labels);
    let mut seen = FxHashSet::default();
    for gate in &circuit.gates {
        match *gate {
            Gate::t(q)
            | Gate::tdg(q)
            | Gate::rz(_, q)
            | Gate::rx(_, q)
            | Gate::ry(_, q)
            | Gate::p(_, q) => {
                let kind = RotationKind::of(gate).unwrap();
                let angle = Angle::of(gate);
                if angle.is_clifford() {
                    if let Some(clifford) = angle.clifford_gate(q) {
                        if matches!(kind, RotationKind::X | RotationKind::Y) {
                            if kind == RotationKind::Y {
                                sketch.apply(&Gate::sdg(q));
                            }
                            sketch.apply(&Gate::h(q));
                        }
                        sketch.apply(&clifford);
                        if matches!(kind, RotationKind::X | RotationKind::Y) {
                            sketch.apply(&Gate::h(q));
                            if kind == RotationKind::Y {
                                sketch.apply(&Gate::s(q));
                            }
                        }
                    }
                } else if !seen.insert(sketch.axis(kind, q)) {
                    return true;
                }
            }
            // No frame change, as in the main scan.
            Gate::measure { .. } | Gate::reset(_) | Gate::ccx { .. } | Gate::ccz { .. } => {}
            Gate::cp { .. }
            | Gate::crx { .. }
            | Gate::cry { .. }
            | Gate::crz { .. }
            | Gate::ch { .. }
            | Gate::cswap { .. } => seen.clear(),
            _ => sketch.apply(gate),
        }
    }
    false
}

/// The state of the main scan.
struct Folder {
    /// Fingerprints of the frame rows.
    sketch: SketchFrame,
    /// The frame rows themselves.
    exact: ExactFrame,
    history: History,
    /// Rotations merged into a later one, by gate index.
    deleted: BitSet,
    /// The merged rotations that replace gates, in gate order.
    replacements: Vec<Replacement>,
    /// The current axis, a buffer reused across gates.
    axis: Vec<u64>,
}

/// A merged rotation emitted in place of the rotation at gate `gate`.
struct Replacement {
    gate: usize,
    qubit: u32,
    angle: Angle,
    kind: RotationKind,
}

impl Folder {
    fn new(
        labels: &[(u128, u128)],
        n: usize,
        l: usize,
        attempt_steps: usize,
        sliced: bool,
        num_gates: usize,
    ) -> Self {
        Self {
            sketch: SketchFrame::new(labels),
            exact: ExactFrame::new(n, l),
            history: History::new(n, l, attempt_steps, sliced),
            deleted: BitSet::new(num_gates),
            replacements: Vec::new(),
            axis: Vec::with_capacity(2 * l),
        }
    }

    /// Scan gate `idx`. `None` once axis storage is full.
    fn step(&mut self, idx: usize, gate: &Gate) -> Option<()> {
        match *gate {
            Gate::t(q)
            | Gate::tdg(q)
            | Gate::rz(_, q)
            | Gate::rx(_, q)
            | Gate::ry(_, q)
            | Gate::p(_, q) => {
                self.rotation(idx, q, RotationKind::of(gate).unwrap(), Angle::of(gate))
            }
            Gate::ccx {
                control1,
                control2,
                target,
            } => self.toffoli([control1, control2], target, true),
            Gate::ccz {
                control1,
                control2,
                target,
            } => self.toffoli([control1, control2], target, false),
            // A rotation about P commutes with a Z-basis measurement of q,
            // whose projectors are (1 ± Z_q)/2, exactly when P commutes with
            // Z_q. It commutes with a reset of q, whose Kraus maps are
            // |0⟩⟨b|, when P acts trivially on q: P commutes with Z_q and X_q.
            // So each is a blocker on the frame rows of those Paulis, and
            // only folds that cross it with an anticommuting axis are refused.
            // The frame stays a valid reference for later axes: folds only
            // move a rotation forward, and both axes are compared after the
            // same Clifford prefix.
            Gate::measure { qubit, .. } => self
                .history
                .record_blocker(&self.exact.z[qubit as usize].words),
            Gate::reset(q) => {
                let q = q as usize;
                self.history.record_blocker(&self.exact.z[q].words)?;
                self.history.record_blocker(&self.exact.x[q].words)
            }
            Gate::cp { .. }
            | Gate::crx { .. }
            | Gate::cry { .. }
            | Gate::crz { .. }
            | Gate::ch { .. }
            | Gate::cswap { .. } => {
                // Start a new segment. The retained frame is merely a common
                // change of Pauli basis for axes within that segment; no
                // candidate can cross the unsupported operation.
                self.history.newest.clear();
                Some(())
            }
            _ => {
                self.apply_clifford(gate);
                Some(())
            }
        }
    }

    /// A single-qubit rotation: merge it into the newest earlier
    /// rotation it can reach, or record it for later ones to merge into.
    fn rotation(&mut self, idx: usize, q: u32, kind: RotationKind, angle: Angle) -> Option<()> {
        if angle.is_clifford() {
            // An Rz by a multiple of pi/2 is part of the frame.
            self.apply_clifford_angle(angle, q, kind);
            return Some(());
        }
        let hash = self.sketch.axis(kind, q);
        let sign = self.exact.axis(kind, q, &mut self.axis);
        let mut angle = angle;
        let mut input_gates = 1;
        if let Some(id) = self.history.partner(hash, &self.axis) {
            // With P_i = s_i W and P_j = s_j W, commuting P_i through the
            // interval gives an angle s_i θ_i + s_j θ_j on W. At j, the
            // physical Z rotation needs that angle multiplied by s_j.
            let prior = self.history.rotation(id);
            if let Some(merged) =
                Angle::merge(prior.angle.clone(), prior.sign * sign, angle.clone())
                && merged.0.emitted_gate_count() <= prior.input_gates + 1
            {
                input_gates = prior.input_gates + 1;
                self.history.kill(id);
                self.deleted.insert(prior.gate);
                self.replacements.push(Replacement {
                    gate: idx,
                    qubit: q,
                    angle: merged.clone(),
                    kind,
                });
                if merged.is_clifford() {
                    // S, S† or Z now sits here and changes every later
                    // axis; the identity leaves nothing.
                    self.apply_clifford_angle(merged.clone(), q, kind);
                    return Some(());
                }
                // A non-Clifford merged rotation stays live and may itself
                // merge with a later rotation on this axis.
                angle = merged;
            }
        }
        let rotation = Rotation {
            gate: idx,
            sign,
            angle,
            input_gates,
        };
        self.history.record_rotation(rotation, hash, &self.axis)
    }

    /// A CCX (`target_x`) or CCZ is not Clifford, so the frame is unchanged.
    /// As in `to_pbc`, it is seven Pauli rotations, about the products of the
    /// nonempty subsets of its controls' Z rows and its target's X (CCX) or Z
    /// (CCZ) row. They stay in the circuit and block the folds they
    /// anticommute with.
    fn toffoli(&mut self, controls: [u32; 2], target: u32, target_x: bool) -> Option<()> {
        let frame = &self.exact;
        let [a, b] = controls.map(|c| &frame.z[c as usize].words);
        let target_rows = if target_x { &frame.x } else { &frame.z };
        let t = &target_rows[target as usize].words;
        // Bits 0, 1 and 2 of `subset` select a, b and t.
        for subset in 1..8 {
            self.axis.clear();
            self.axis.resize(t.len(), 0);
            for (i, row) in [a, b, t].into_iter().enumerate() {
                if subset >> i & 1 != 0 {
                    self.axis.iter_mut().zip(row).for_each(|(w, r)| *w ^= r);
                }
            }
            self.history.record_blocker(&self.axis)?;
        }
        Some(())
    }

    fn apply_clifford(&mut self, gate: &Gate) {
        self.sketch.apply(gate);
        self.exact.apply(gate);
    }

    /// Apply the S, S† or Z a Clifford angle amounts to (nothing for zero).
    fn apply_clifford_angle(&mut self, angle: Angle, q: u32, kind: RotationKind) {
        if matches!(kind, RotationKind::X | RotationKind::Y) {
            if kind == RotationKind::Y {
                self.apply_clifford(&Gate::sdg(q));
            }
            self.apply_clifford(&Gate::h(q));
        }
        if let Some(gate) = angle.clifford_gate(q) {
            self.apply_clifford(&gate);
        }
        if matches!(kind, RotationKind::X | RotationKind::Y) {
            self.apply_clifford(&Gate::h(q));
            if kind == RotationKind::Y {
                self.apply_clifford(&Gate::s(q));
            }
        }
    }

    /// The circuit with the merges applied: merged-away rotations dropped and
    /// the others replaced by their merged angles. An Rz by an exact multiple
    /// of pi/4 is written as the T, S and Z gates it amounts to (none for zero),
    /// as `PhaseFoldRand` does, whether or not it folded.
    fn rewrite(self, circuit: &Circuit) -> Circuit {
        let exact_rz =
            |g: &Gate| matches!(g, Gate::rz(theta, _) if theta.quarter_turns().is_some());
        if self.replacements.is_empty() && !circuit.gates.iter().any(exact_rz) {
            return circuit.clone();
        }
        let mut output = circuit.empty_like();
        output
            .gates
            .reserve(circuit.gates.len() - self.replacements.len());
        let mut replacements = self.replacements.iter().peekable();
        for (idx, gate) in circuit.gates.iter().enumerate() {
            // A replacement may itself have merged into a later rotation.
            let replacement = replacements.next_if(|r| r.gate == idx);
            if self.deleted.contains(idx) {
                continue;
            }
            match replacement {
                Some(r) => r.angle.emit(&mut output, r.qubit, r.kind),
                None => match gate {
                    Gate::rz(_, q) => Angle::of(gate).emit(&mut output, *q, RotationKind::Z),
                    _ => output.apply(gate.clone()),
                },
            }
        }
        output
    }
}

/// The events recorded so far, rotations and blockers, with ids in order of
/// recording, and the indexes a fold attempt searches.
struct History {
    /// Words per bit plane of an axis.
    l: usize,
    /// Work allowed per fold attempt.
    attempt_steps: usize,
    events: Vec<Event>,
    /// Exact unsigned axes, `2 l` words per event, so event `id`'s axis
    /// starts at `2 l id`.
    axes: Vec<u64>,
    /// Support signature of each event.
    sigs: Vec<u64>,
    /// The live events: blockers always, rotations until they merge away.
    /// Kept apart from `events` so scans read one bit per event.
    live: BitSet,
    /// Newest rotation per axis fingerprint. Older ones are chained through
    /// [`Event::older`].
    newest: FxHashMap<u128, usize>,
    lists: SupportLists,
}

struct Event {
    /// `None` for a blocker, which is never a merge candidate.
    rotation: Option<Rotation>,
    /// The next older rotation with the same fingerprint.
    older: Option<usize>,
}

#[derive(Clone)]
struct Rotation {
    /// Index of its gate in the input circuit.
    gate: usize,
    /// Its axis is `sign` times the stored unsigned axis.
    sign: i8,
    angle: Angle,
    input_gates: usize,
}

impl History {
    fn new(n: usize, l: usize, attempt_steps: usize, sliced: bool) -> Self {
        let lists = if sliced {
            assert!(n <= MAX_SLICED_QUBITS);
            SupportLists::Sliced((0..n).map(|_| SlicedList::default()).collect())
        } else {
            SupportLists::Scalar(Box::new(std::array::from_fn(|_| ScalarList::default())))
        };
        Self {
            l,
            attempt_steps,
            events: Vec::new(),
            axes: Vec::new(),
            sigs: Vec::new(),
            live: BitSet::new(0),
            newest: FxHashMap::default(),
            lists,
        }
    }

    fn axis(&self, id: usize) -> &[u64] {
        let words = 2 * self.l;
        &self.axes[words * id..words * (id + 1)]
    }

    fn rotation(&self, id: usize) -> Rotation {
        self.events[id]
            .rotation
            .clone()
            .expect("candidates are rotations")
    }

    /// Record a rotation as a merge candidate. `None` once axis storage is
    /// full.
    fn record_rotation(&mut self, rotation: Rotation, hash: u128, axis: &[u64]) -> Option<()> {
        let id = self.push(axis, Some(rotation))?;
        self.events[id].older = self.newest.insert(hash, id);
        Some(())
    }

    /// Record a fixed obstacle: a later fold across it must commute with
    /// `axis`. `None` once axis storage is full.
    fn record_blocker(&mut self, axis: &[u64]) -> Option<()> {
        self.push(axis, None).map(drop)
    }

    fn push(&mut self, axis: &[u64], rotation: Option<Rotation>) -> Option<usize> {
        debug_assert_eq!(axis.len(), 2 * self.l);
        if self.axes.len() + axis.len() > MAX_AXIS_WORDS {
            return None;
        }
        let id = self.events.len();
        let signature = support_signature(axis, self.l);
        self.lists.push(id as u32, axis, signature);
        self.axes.extend_from_slice(axis);
        self.sigs.push(signature);
        self.live.insert(id);
        self.events.push(Event {
            rotation,
            older: None,
        });
        Some(id)
    }

    /// Mark a rotation that merged away dead.
    fn kill(&mut self, id: usize) {
        self.live.remove(id);
        self.lists
            .kill(id as u32, self.sigs[id], &self.live, &self.axes);
    }

    /// The rotation `axis` can merge into, if one fold attempt finds it
    /// within the budget: the newest live rotation with exactly this axis,
    /// provided every live event since commutes with `axis`. (Any older one
    /// would be blocked too, since the newest also commutes with `axis`.)
    fn partner(&mut self, hash: u128, axis: &[u64]) -> Option<usize> {
        let mut budget = Budget(self.attempt_steps);
        let id = self.candidate(hash, axis, &mut budget)?;
        self.commutes_after(id, axis, &mut budget).then_some(id)
    }

    /// The newest live rotation with exactly `axis`, on the chain of
    /// fingerprint `hash`. One step per axis compared.
    fn candidate(&mut self, hash: u128, axis: &[u64], budget: &mut Budget) -> Option<usize> {
        let mut cursor = self.newest_live(self.newest.get(&hash).copied());
        while let Some(id) = cursor {
            if !budget.spend(1) {
                return None;
            }
            if self.axis(id) == axis {
                return Some(id);
            }
            cursor = self.newest_live(self.events[id].older);
        }
        None
    }

    /// The first live rotation on the fingerprint chain from `start`. The
    /// dead links passed on the way are pointed at it, so a long run of
    /// merges does not make every later lookup walk all of them again.
    fn newest_live(&mut self, start: Option<usize>) -> Option<usize> {
        let mut found = start;
        while let Some(id) = found.filter(|&id| !self.live.contains(id)) {
            found = self.events[id].older;
        }
        let mut cursor = start;
        while cursor != found {
            let id = cursor.expect("`found` is on the chain");
            cursor = std::mem::replace(&mut self.events[id].older, found);
        }
        found
    }

    /// Whether `axis` commutes with every live event recorded after event
    /// `after`. Running out of budget counts as not commuting.
    ///
    /// Only events that share a support-signature bit with `axis` can fail to
    /// commute, so only the lists of those bits are scanned. Each list is in
    /// id order, so the events after `after` are a suffix, found by binary
    /// search, and scanned oldest first.
    fn commutes_after(&self, after: usize, axis: &[u64], budget: &mut Budget) -> bool {
        match &self.lists {
            SupportLists::Scalar(lists) => self.commutes_scalar(lists, after, axis, budget),
            SupportLists::Sliced(lists) => commutes_sliced(lists, after, axis, budget),
        }
    }

    /// [`History::commutes_after`] on scalar lists, one step per entry
    /// visited. An event that shares several bits with `axis` is in several
    /// of the lists, but is tested only in the one of the lowest shared bit.
    fn commutes_scalar(
        &self,
        lists: &[ScalarList; 64],
        after: usize,
        axis: &[u64],
        budget: &mut Budget,
    ) -> bool {
        let signature = support_signature(axis, self.l);
        for bit in bits(signature) {
            let ids = &lists[bit].ids;
            let start = ids.partition_point(|&id| id as usize <= after);
            for &id in &ids[start..] {
                if !budget.spend(1) {
                    return false;
                }
                let id = id as usize;
                if self.live.contains(id)
                    && (self.sigs[id] & signature).trailing_zeros() as usize == bit
                    && packed_anticommutes(self.axis(id), axis, self.l)
                {
                    return false;
                }
            }
        }
        true
    }
}

/// [`History::commutes_after`] on bit-sliced lists, one step per live entry
/// of each block scanned. An event in several of the lists is tested in each.
///
/// An axis here is one X word and one Z word, and entry `e` anticommutes with
/// `axis = (ax, az)` when `(ax ∧ z_e) ⊕ (az ∧ x_e)` has odd weight. Grouping
/// that parity by qubit, it is, for all 64 entries of a block at once, the
/// XOR of the block's Z masks of the qubits where `axis` has X and its X masks
/// of the qubits where `axis` has Z.
fn commutes_sliced(lists: &[SlicedList], after: usize, axis: &[u64], budget: &mut Budget) -> bool {
    let (ax, az) = (axis[0], axis[1]);
    let n = lists.len();
    // The mask words to XOR in each block (see [`SlicedList`] for the layout).
    let mut words = [0u8; 2 * MAX_SLICED_QUBITS];
    let mut count = 0;
    for p in bits(ax) {
        words[count] = (2 * p + 1) as u8;
        count += 1;
    }
    for p in bits(az) {
        words[count] = (2 * p) as u8;
        count += 1;
    }
    let words = &words[..count];
    for q in bits(ax | az) {
        let list = &lists[q];
        let start = list.ids.partition_point(|&id| id as usize <= after);
        // The entries before `start` in its block are not after `after`.
        let mut first = !0u64 << (start % 64);
        for block in start / 64..list.live.len() {
            if !budget.spend(list.live[block].count_ones() as usize) {
                return false;
            }
            let masks = list.block(block, n);
            let mut anticommuting = 0;
            for &w in words {
                anticommuting ^= masks[w as usize];
            }
            if anticommuting & list.live[block] & first != 0 {
                return false;
            }
            first = !0;
        }
    }
    true
}

/// For each support-signature bit, the events whose axis has it, in id order.
enum SupportLists {
    /// Plain lists of ids, for circuits wider than [`MAX_SLICED_QUBITS`].
    Scalar(Box<[ScalarList; 64]>),
    /// Bit-sliced lists, one per qubit. A circuit this narrow has one word
    /// per plane, and the signature bit of qubit `q` is `q`.
    Sliced(Vec<SlicedList>),
}

impl SupportLists {
    fn push(&mut self, id: u32, axis: &[u64], signature: u64) {
        match self {
            Self::Scalar(lists) => {
                for bit in bits(signature) {
                    lists[bit].ids.push(id);
                }
            }
            Self::Sliced(lists) => {
                let n = lists.len();
                for q in bits(signature) {
                    lists[q].push(id, axis[0], axis[1], n);
                }
            }
        }
    }

    /// Mark event `id`, already removed from `live`, dead in its lists.
    fn kill(&mut self, id: u32, signature: u64, live: &BitSet, axes: &[u64]) {
        match self {
            Self::Scalar(lists) => {
                for bit in bits(signature) {
                    lists[bit].note_dead(live);
                }
            }
            Self::Sliced(lists) => {
                let n = lists.len();
                for q in bits(signature) {
                    lists[q].kill(id, live, axes, n);
                }
            }
        }
    }
}

/// Event ids, increasing. Dead ids stay until the list is compacted.
#[derive(Default)]
struct ScalarList {
    ids: Vec<u32>,
    /// Dead ids still in the list.
    dead: usize,
}

impl ScalarList {
    /// Count one more dead id, and drop the dead ids once a quarter of the
    /// list is dead. At least a quarter of the list died since the previous
    /// compaction, so each costs O(1) per death.
    fn note_dead(&mut self, live: &BitSet) {
        self.dead += 1;
        if 4 * self.dead > self.ids.len() {
            self.ids.retain(|&id| live.contains(id as usize));
            self.dead = 0;
        }
    }
}

/// One qubit's support list, bit-sliced in blocks of 64 entries. Entry `i` of
/// a block owns bit `i` of each of the block's `2 n` masks: mask `2 p` marks
/// the entries with X or Y on qubit `p`, and mask `2 p + 1` those with Z or Y.
#[derive(Default)]
struct SlicedList {
    /// Event ids, increasing.
    ids: Vec<u32>,
    /// `2 n` words per block.
    masks: Vec<u64>,
    /// The live entries, one word per block.
    live: Vec<u64>,
    /// Dead entries still in the list.
    dead: usize,
}

impl SlicedList {
    fn block(&self, block: usize, n: usize) -> &[u64] {
        &self.masks[2 * n * block..2 * n * (block + 1)]
    }

    /// Append an event with axis `(x, z)`.
    fn push(&mut self, id: u32, x: u64, z: u64, n: usize) {
        let pos = self.ids.len();
        if pos.is_multiple_of(64) {
            self.masks.resize(self.masks.len() + 2 * n, 0);
            self.live.push(0);
        }
        let (block, bit) = (pos / 64, 1u64 << (pos % 64));
        let masks = &mut self.masks[2 * n * block..2 * n * (block + 1)];
        for p in bits(x) {
            masks[2 * p] |= bit;
        }
        for p in bits(z) {
            masks[2 * p + 1] |= bit;
        }
        self.live[block] |= bit;
        self.ids.push(id);
    }

    /// Mark event `id`, already removed from `live`, dead, and rebuild the
    /// list from its live entries once half of it is dead.
    fn kill(&mut self, id: u32, live: &BitSet, axes: &[u64], n: usize) {
        // Rotations usually merge away soon after they are recorded, so
        // gallop back from the end.
        let ids = &self.ids;
        let (mut lo, mut span) = (ids.len() - 1, 1);
        while lo > 0 && ids[lo] > id {
            lo = lo.saturating_sub(span);
            span *= 2;
        }
        let pos = lo + ids[lo..].partition_point(|&e| e < id);
        debug_assert_eq!(ids[pos], id, "a live event is in its lists");
        self.live[pos / 64] &= !(1 << (pos % 64));
        self.dead += 1;
        if 2 * self.dead > self.ids.len() {
            self.rebuild(live, axes, n);
        }
    }

    fn rebuild(&mut self, live: &BitSet, axes: &[u64], n: usize) {
        let ids = std::mem::take(&mut self.ids);
        self.masks.clear();
        self.live.clear();
        self.dead = 0;
        for id in ids {
            let i = id as usize;
            if live.contains(i) {
                self.push(id, axes[2 * i], axes[2 * i + 1], n);
            }
        }
    }
}

/// The work left in one fold attempt.
struct Budget(usize);

impl Budget {
    /// Spend `steps`, or return `false` if fewer are left.
    fn spend(&mut self, steps: usize) -> bool {
        match self.0.checked_sub(steps) {
            Some(left) => {
                self.0 = left;
                true
            }
            None => false,
        }
    }
}

/// A growable set of small integers, one bit each.
struct BitSet(Vec<u64>);

impl BitSet {
    /// An empty set with room for `0..len`.
    fn new(len: usize) -> Self {
        Self(vec![0; len.div_ceil(64)])
    }

    fn insert(&mut self, i: usize) {
        let word = i / 64;
        if word >= self.0.len() {
            self.0.resize(word + 1, 0);
        }
        self.0[word] |= 1 << (i % 64);
    }

    fn remove(&mut self, i: usize) {
        self.0[i / 64] &= !(1 << (i % 64));
    }

    fn contains(&self, i: usize) -> bool {
        self.0[i / 64] >> (i % 64) & 1 != 0
    }
}

/// The positions of the set bits of `w`, lowest first.
fn bits(mut w: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        (w != 0).then(|| {
            let bit = w.trailing_zeros() as usize;
            w &= w - 1;
            bit
        })
    })
}

/// The fingerprints of the frame rows: each row's XOR of the labels of its
/// factors (`Y_q` counts as `X_q` and `Z_q`). Products of rows XOR their
/// fingerprints, and signs do not change them, so this follows
/// [`ExactFrame::apply`] at O(1) per gate.
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

    fn axis(&self, kind: RotationKind, q: u32) -> u128 {
        let q = q as usize;
        match kind {
            RotationKind::X => self.x[q],
            RotationKind::Y => self.x[q] ^ self.z[q],
            RotationKind::Z | RotationKind::Phase => self.z[q],
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
            Gate::swap(a, b) => {
                self.x.swap(a as usize, b as usize);
                self.z.swap(a as usize, b as usize);
            }
            Gate::sx(q) => {
                self.apply(&Gate::h(q));
                self.apply(&Gate::s(q));
                self.apply(&Gate::h(q));
            }
            Gate::cy { control, target } => {
                self.apply(&Gate::sdg(target));
                self.apply(&Gate::cnot { control, target });
                self.apply(&Gate::s(target));
            }
            Gate::x(_) | Gate::y(_) | Gate::z(_) => {} // only the exact row sign changes
            _ => unreachable!("supported Clifford gate"),
        }
    }
}

/// The frame rows `U† X_q U` and `U† Z_q U` for the Clifford prefix `U`.
struct ExactFrame {
    x: Vec<Row>,
    z: Vec<Row>,
}

impl ExactFrame {
    /// The identity frame on `n` qubits, with `l` words per bit plane.
    fn new(n: usize, l: usize) -> Self {
        Self {
            x: (0..n).map(|q| Row::generator(l, q, true)).collect(),
            z: (0..n).map(|q| Row::generator(l, q, false)).collect(),
        }
    }

    /// Populate a reusable axis buffer and return its Hermitian sign.
    fn axis(&self, kind: RotationKind, q: u32, axis: &mut Vec<u64>) -> i8 {
        let q = q as usize;
        let row = match kind {
            RotationKind::X | RotationKind::Y => &self.x[q],
            RotationKind::Z | RotationKind::Phase => &self.z[q],
        };
        axis.clone_from(&row.words);
        if kind != RotationKind::Y {
            return row.sign();
        }
        // Y = i X Z. Multiplication must retain its phase before converting
        // to the canonical Hermitian word: XOR alone loses the axis sign.
        let z = &self.z[q];
        let l = axis.len() / 2;
        let crossings = row.words[l..]
            .iter()
            .zip(&z.words[..l])
            .fold(0, |acc, (a, b)| acc ^ (a & b));
        let phase = (1 + row.phase + z.phase + 2 * (crossings.count_ones() & 1) as u8) & 3;
        for (a, b) in axis.iter_mut().zip(&z.words) {
            *a ^= *b;
        }
        Row::word_sign(axis, phase)
    }

    fn apply(&mut self, gate: &Gate) {
        match *gate {
            Gate::h(q) => {
                let q = q as usize;
                std::mem::swap(&mut self.x[q], &mut self.z[q]);
            }
            Gate::x(q) => self.z[q as usize].negate(),
            Gate::z(q) => self.x[q as usize].negate(),
            Gate::y(q) => {
                self.x[q as usize].negate();
                self.z[q as usize].negate();
            }
            Gate::swap(a, b) => {
                self.x.swap(a as usize, b as usize);
                self.z.swap(a as usize, b as usize);
            }
            Gate::sx(q) => {
                self.apply(&Gate::h(q));
                self.apply(&Gate::s(q));
                self.apply(&Gate::h(q));
            }
            Gate::cy { control, target } => {
                self.apply(&Gate::sdg(target));
                self.apply(&Gate::cnot { control, target });
                self.apply(&Gate::s(target));
            }
            Gate::s(q) | Gate::sdg(q) => {
                // S† X S = -Y = i³ X Z and S X S† = Y = i X Z.
                let q = q as usize;
                self.x[q].multiply(&self.z[q]);
                let phase = if matches!(gate, Gate::s(_)) { 3 } else { 1 };
                self.x[q].phase = (self.x[q].phase + phase) & 3;
            }
            Gate::cnot { control, target } => {
                let (c, t) = (control as usize, target as usize);
                let [xc, xt] = self.x.get_disjoint_mut([c, t]).expect("distinct operands");
                xc.multiply(xt);
                let [zt, zc] = self.z.get_disjoint_mut([t, c]).expect("distinct operands");
                zt.multiply(zc); // these two rows commute
            }
            Gate::cz { control, target } => {
                let (c, t) = (control as usize, target as usize);
                self.x[c].multiply(&self.z[t]);
                self.x[t].multiply(&self.z[c]); // these two rows commute
            }
            _ => unreachable!("supported Clifford gate"),
        }
    }
}

/// The Pauli string `i^phase X^x Z^z`; `words` holds the `x` plane and then
/// the `z` plane, `l` words each.
struct Row {
    words: Vec<u64>,
    phase: u8,
}

impl Row {
    /// `X_q` if `x`, else `Z_q`.
    fn generator(l: usize, q: usize, x: bool) -> Self {
        let mut words = vec![0; 2 * l];
        let plane = if x { 0 } else { l };
        words[plane + q / 64] = 1 << (q % 64);
        Self { words, phase: 0 }
    }

    fn negate(&mut self) {
        self.phase = (self.phase + 2) & 3;
    }

    /// `self ← self · rhs`.
    fn multiply(&mut self, rhs: &Self) {
        let l = self.words.len() / 2;
        // Moving rhs's X factors left past self's Z factors flips the sign
        // once per qubit where both are set.
        let crossings = self.words[l..]
            .iter()
            .zip(&rhs.words[..l])
            .fold(0, |acc, (z, x)| acc ^ (z & x));
        self.phase = (self.phase + rhs.phase + 2 * (crossings.count_ones() & 1) as u8) & 3;
        for (a, b) in self.words.iter_mut().zip(&rhs.words) {
            *a ^= *b;
        }
    }

    /// The sign `s` with `self = s W` for the Hermitian Pauli string `W`
    /// with these planes (whose `Y = i X Z` factors carry the phase
    /// `i^|x ∧ z|`).
    fn sign(&self) -> i8 {
        Self::word_sign(&self.words, self.phase)
    }

    fn word_sign(words: &[u64], phase: u8) -> i8 {
        let (x, z) = words.split_at(words.len() / 2);
        let ys: u32 = x.iter().zip(z).map(|(x, z)| (x & z).count_ones()).sum();
        match (phase + 4 - (ys & 3) as u8) & 3 {
            0 => 1,
            2 => -1,
            _ => unreachable!("Clifford conjugation preserves Hermiticity"),
        }
    }
}

/// Physical rotation at the later gate's position. P shares the Z axis but
/// retains its native spelling; its scalar phase is irrelevant to this pass.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RotationKind {
    X,
    Y,
    Z,
    Phase,
}

impl RotationKind {
    fn of(gate: &Gate) -> Option<Self> {
        match gate {
            Gate::rx(..) => Some(Self::X),
            Gate::ry(..) => Some(Self::Y),
            Gate::rz(..) | Gate::t(_) | Gate::tdg(_) => Some(Self::Z),
            Gate::p(..) => Some(Self::Phase),
            _ => None,
        }
    }

    fn gate(self, theta: crate::angle::Angle, q: u32) -> Gate {
        match self {
            Self::X => Gate::rx(theta, q),
            Self::Y => Gate::ry(theta, q),
            Self::Z => Gate::rz(theta, q),
            Self::Phase => Gate::p(theta, q),
        }
    }
}

#[derive(Clone)]
struct Angle(crate::angle::Angle);
impl Angle {
    fn of(gate: &Gate) -> Self {
        Self(crate::angle::Angle::of_gate(gate).expect("rotation gate"))
    }
    fn merge(prior: Self, sign: i8, current: Self) -> Option<Self> {
        prior.0.checked_add_signed(sign, &current.0).ok().map(Self)
    }
    fn is_clifford(&self) -> bool {
        self.0.quarter_turns().is_some_and(|k| k & 1 == 0)
    }
    fn clifford_gate(&self, q: u32) -> Option<Gate> {
        match self.0.quarter_turns() {
            Some(2) => Some(Gate::s(q)),
            Some(4) => Some(Gate::z(q)),
            Some(6) => Some(Gate::sdg(q)),
            _ => None,
        }
    }
    fn emit(&self, output: &mut Circuit, q: u32, kind: RotationKind) {
        if kind == RotationKind::Z {
            self.0.emit(output, q);
        } else if !self.0.is_zero() {
            output.apply(kind.gate(self.0.normalized(), q));
        }
    }
}

#[cfg(test)]
mod rotation_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cancel::CancelGates;
    use crate::pass::{count_rz, count_t};
    use crate::pbc::{Pauli, Phase, to_pbc};
    use crate::phase_fold_rand::phase_fold_rand;
    use crate::qasm;
    use crate::unitary::{C, circuit_unitary, circuits_equiv};

    /// The pass as [`phase_fold_pauli`] runs it, with the given labels.
    fn fold_with_labels(c: &Circuit, labels: &[(u128, u128)]) -> Circuit {
        fold(
            c,
            labels,
            MAX_ATTEMPT_STEPS,
            c.num_qubits <= MAX_SLICED_QUBITS,
        )
    }

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
        let out = phase_fold_pauli(&c);
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
        let out = phase_fold_pauli(&c);
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
        let out = phase_fold_pauli(&c);
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
        let out = phase_fold_pauli(&c);
        assert_eq!(count_t(&out), 0);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn numeric_and_named_angles_merge_across_commuting_rotations() {
        let mut c = Circuit::new(2);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::rz_f64(0.37, 1).unwrap());
        c.apply(Gate::h(0));
        c.apply(Gate::rz_f64(0.21, 0).unwrap());
        let out = phase_fold_pauli(&c);
        assert_eq!(count_t(&out), 1);
        assert_eq!(count_rz(&out), 2);
        assert_eq!(out.gates.len(), c.gates.len());
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn numeric_pi_can_cancel_named_t_in_a_mixed_fold() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::x(0));
        c.apply(Gate::t(0));
        let out = phase_fold_pauli(&c);
        assert_eq!(out.gates, vec![Gate::x(0)]);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn named_t_can_fold_to_clifford_and_change_the_frame() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = phase_fold_pauli(&c);
        assert_eq!(count_rz(&out), 0);
        assert_eq!(out.gates.first(), Some(&Gate::s(0)));
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn anticommuting_rz_blocks_t_fold() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::rz_f64(0.3, 0).unwrap());
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = phase_fold_pauli(&c);
        assert_eq!(out.gates, c.gates);
    }

    #[test]
    fn clifford_rz_updates_frame_and_zero_rz_does_not_block() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::rz_f64(0.0, 0).unwrap());
        c.apply(Gate::h(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 2.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::t(0));
        let out = phase_fold_pauli(&c);
        assert_eq!(count_t(&out), 0);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn near_quarter_rz_is_not_rounded_away() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz_f64(PI / 4.0 + 1e-10, 0).unwrap());
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        let out = phase_fold_pauli(&c);
        assert_eq!(count_rz(&out), 1);
        assert_eq!(count_t(&out), 1);
        assert!(circuits_equiv(&c, &out, 1e-12));
    }

    #[test]
    fn rz_axis_across_packed_word_boundary() {
        let mut c = Circuit::new(65);
        let cx = Gate::cnot {
            control: 0,
            target: 64,
        };
        c.apply(Gate::rz_f64(0.25, 64).unwrap());
        c.apply(cx.clone());
        c.apply(Gate::rz_f64(0.2, 64).unwrap());
        c.apply(cx);
        c.apply(Gate::rz_f64(0.5, 64).unwrap());
        let out = phase_fold_pauli(&c);
        assert_eq!(count_rz(&out), 2);
        assert_eq!(out.gates.len(), 4);
        assert!(
            matches!(out.gates.last(), Some(Gate::rz(theta, 64)) if (theta.to_f64_lossy().unwrap() - 0.75).abs() < 1e-12)
        );
    }

    #[test]
    fn forced_fingerprint_collisions_cannot_authorize_rewrites() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::h(0));
        c.apply(Gate::t(0));
        let out = fold_with_labels(&c, &[(0, 0)]);
        assert_eq!(count_t(&out), 2);
        assert!(circuits_equiv(&c, &out, 1e-10));
    }

    #[test]
    fn mod5_4_improves_on_phase_fold_and_preserves_unitary() {
        let input = qasm::parse(include_str!("../benchmarks/feynman/mod5_4.qasm")).unwrap();
        let cancelled = CancelGates.run(&input);
        let baseline = phase_fold_rand(&cancelled);
        let improved = phase_fold_pauli(&baseline);
        assert!(count_t(&improved) < count_t(&baseline));
        assert!(circuits_equiv(&baseline, &improved, 1e-9));
    }

    #[test]
    fn exact_frame_agrees_with_pbc_frame() {
        let mut rng = Rng(0x7e57_1234_9abc_def0);
        let mut c = Circuit::new(4);
        let mut frame = ExactFrame::new(4, 1);
        let l = 1;
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
                        row.sign(),
                        if expanded.phase == Phase::One { 1 } else { -1 }
                    );
                    let mut projected = 0;
                    for (bit, factor) in expanded.factors.iter().enumerate() {
                        let x = (row.words[bit / 64] >> (bit % 64)) & 1 != 0;
                        let z = (row.words[l + bit / 64] >> (bit % 64)) & 1 != 0;
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
            let out = phase_fold_pauli(&baseline);
            assert!(count_t(&out) <= count_t(&baseline));
            assert!(out.gates.len() <= baseline.gates.len());
            assert!(circuits_equiv(&baseline, &out, 1e-9));
            if case < 100 {
                let collided = fold_with_labels(&baseline, &vec![(0, 0); n]);
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
                    8..=10 => Gate::rz_f64((rng.up_to(17) as f64 - 8.0) / 13.0, q).unwrap(),
                    11 => Gate::rz_f64((rng.up_to(9) as f64 - 4.0) * (PI / 4.0), q).unwrap(),
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
            let out = phase_fold_pauli(&c);
            assert!(count_t(&out) + count_rz(&out) <= count_t(&c) + count_rz(&c));
            assert!(out.gates.len() <= c.gates.len());
            assert!(circuits_equiv(&c, &out, 1e-9), "case {case}");
            if case < 100 {
                let labels = vec![(0, 0); n];
                let collided = fold_with_labels(&c, &labels);
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
        let out = phase_fold_pauli(&c);
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
            let out = phase_fold_pauli(&c);
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
        let out = phase_fold_pauli(&c);
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
            let out = phase_fold_pauli(&c);
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
            let out = phase_fold_pauli(&c);
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
        let out = phase_fold_pauli(&c);
        assert_eq!(count_t(&out), 1, "{:?}", out.gates);
        let at = out.gates.iter().position(|g| *g == Gate::reset(0)).unwrap();
        assert_eq!(count_t(&c.replacing_gates(out.gates[at..].to_vec())), 0);
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
        let out = fold(&c, &labels, BUDGET, true);
        // The Z0 pair is out of reach within one attempt; the Z1 pair folds.
        assert_eq!(count_t(&out), count_t(&c) - 2);
        assert_eq!(out.gates.last(), Some(&Gate::s(1)));
        // With the default budget the Z0 pair folds too.
        let full = fold_with_labels(&c, &labels);
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

    /// The bit-sliced and scalar support lists reach the same verdicts, so
    /// with a budget neither can exhaust they rewrite identically.
    #[test]
    fn sliced_and_scalar_lists_agree() {
        let mut rng = Rng(0x51ce_d0ff_5ca1_a400);
        let labels: Vec<_> = (0..8)
            .map(|q| (q as u128 * 7 + 1, q as u128 * 13 + 5))
            .collect();
        for case in 0..600 {
            let n = 2 + case % 7;
            let mut c = Circuit::with_cbits(n, n);
            for _ in 0..(20 + rng.up_to(120)) {
                c.apply(random_gate(&mut rng, n, case % 3 == 0));
            }
            let sliced = fold(&c, &labels[..n], usize::MAX, true);
            let scalar = fold(&c, &labels[..n], usize::MAX, false);
            assert_eq!(sliced.gates, scalar.gates, "case {case}");
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
            let out = phase_fold_pauli(&c);
            assert!(count_t(&out) <= count_t(&c));
            assert!(circuits_equiv(&c, &out, 1e-9), "case {case}\n{c}");
            folded += usize::from(out.gates != c.gates);
            if case < 100 {
                let collided = fold_with_labels(&c, &vec![(0, 0); n]);
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
            let out = phase_fold_pauli(&c);
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
            let mut outputs = vec![phase_fold_pauli(&c)];
            if case < 50 {
                outputs.push(fold_with_labels(&c, &vec![(0, 0); n]));
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
