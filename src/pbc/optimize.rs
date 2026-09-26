//! T-count optimization of PBC rotations: commuting-rotation merging, moving
//! Clifford rotations into the output frame, and generalized
//! multiproduct-commutation (MCR) group swaps.
//!
//! See `docs/pbc-phase-folding.tex` for the design and `docs/pbc-optimize.md`
//! for how this implementation differs and what it achieves. Every axis
//! (rotations, measurements, conditional rotations, and the output frame's
//! images) is packed once into canonical bit planes and interned. Then rounds
//! of two steps run until a round removes no T:
//!
//! 1. **Stream.** One pass over all operations. A rotation merges into the
//!    latest earlier rotation with the same axis when every rotation between
//!    them commutes with it, at most `lookback` rotations back; an
//!    interned-axis table finds that candidate in O(1), and only rotations
//!    sharing a support-signature bit are checked as blockers. With
//!    `clifford_to_frame`, Clifford rotations (given, or produced by a merge)
//!    move into a Clifford F carried to the end, and every later axis is
//!    conjugated by F; F is finally composed into the output frame.
//!    Measurements and conditional rotations end merging but not F.
//! 2. **MCR swaps.** Partition each segment greedily into internally
//!    commuting runs of at most `window` rotations. For neighboring runs
//!    A|B|C (and up to `candidates` variants with a shorter A or C) that share
//!    an axis between A and C, swap A,B (or B,C) when the generators commute,
//!    `[H_A, H_B] = 0`, and keep the swap only if the circuit then shrinks.
//!
//! `Strategy::Litinski` replaces both steps with Litinski's greedy layering,
//! as a baseline (see `litinski.rs`). Output rotations reuse existing arena
//! handles where their axis is unchanged.

mod clifford;
mod litinski;
mod words;

use super::pauli::{PauliArena, PauliNode, PauliRef};
use super::{Pauli, PauliAngle, PauliAxis, PbcCircuit, PbcError, PbcOp, Phase};
use clifford::Clifford;
use words::{Axes, Rot, anticommutes, normalize, pack, product};

/// Limits for [`PbcCircuit::optimize_rotations`].
#[derive(Clone, Copy, Debug)]
pub struct OptimizeOptions {
    /// Rotations a merge candidate may be separated by.
    pub lookback: usize,
    /// Maximum length of each group in an MCR swap. 0 disables swaps.
    pub window: usize,
    /// Per run triple, how many shorter group variants to try besides the
    /// full runs: suffixes of A and prefixes of C.
    pub candidates: usize,
    /// Maximum rounds of streaming merge and MCR swaps; also bounds the
    /// merge sweeps of `Strategy::Litinski`.
    pub rounds: usize,
    /// Move Clifford (even-angle) rotations into the output frame, which can
    /// unblock further merges.
    pub clifford_to_frame: bool,
    /// Accept a certified MCR swap whenever it lets rotations merge, even if
    /// the T count does not drop immediately; otherwise only when it does.
    pub eager_swaps: bool,
    /// Keep merged Clifford rotations as merge targets until an
    /// anticommuting rotation, a barrier, or the end forces them into the
    /// frame, instead of moving them into the frame as soon as they appear.
    pub lazy_cliffords: bool,
    /// Which algorithm to run.
    pub strategy: Strategy,
    /// Also report the T depth before and after (costs one layering pass
    /// each).
    pub measure_depth: bool,
    /// Bound on packed axis storage, in u64 words.
    pub max_packed_words: usize,
}

impl Default for OptimizeOptions {
    fn default() -> Self {
        Self {
            lookback: 1 << 20,
            window: 8,
            candidates: 16,
            rounds: 8,
            clifford_to_frame: true,
            eager_swaps: true,
            lazy_cliffords: false,
            strategy: Strategy::Merge,
            measure_depth: false,
            max_packed_words: 64 << 20,
        }
    }
}

/// The optimization algorithm.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Strategy {
    /// Streaming merge with Clifford moves, plus MCR swaps (this module).
    #[default]
    Merge,
    /// Litinski's greedy layering ("A Game of Surface Codes", Sec. 1),
    /// as a baseline: partition rotations into layers of mutually commuting
    /// rotations, repeatedly move rotations to the previous layer when they
    /// commute with all of it, and combine equal rotations that meet in a
    /// layer, commuting the resulting Clifford to the end.
    Litinski,
}

/// What [`PbcCircuit::optimize_rotations`] did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OptimizeStats {
    pub t_before: usize,
    pub t_after: usize,
    pub rotations_before: usize,
    pub rotations_after: usize,
    /// Pairs of rotations merged into one.
    pub merges: usize,
    /// Accepted MCR group swaps.
    pub swaps: usize,
    /// Clifford rotations moved into the output frame.
    pub cliffords_to_frame: usize,
    /// T depth (layers of mutually commuting odd-angle rotations, placed as
    /// early as possible), when `measure_depth` is set.
    pub t_depth_before: usize,
    pub t_depth_after: usize,
    /// Total Pauli weight of the rotation axes: the number of non-identity
    /// single-qubit factors, summed over all rotations.
    pub weight_before: usize,
    pub weight_after: usize,
}

/// An operation during optimization: a rotation, or a barrier (measurement
/// or conditional rotation, by index into the original operations) with its
/// current axis, which evaluates to `sign * W`.
#[derive(Clone, Copy, Debug)]
enum Item {
    Rot(Rot),
    Barrier { op: usize, axis: u32, sign: i8 },
}

impl PbcCircuit {
    /// Odd-angle rotations, conditional or not: each costs one T (pi/8)
    /// rotation. Identity-axis rotations are counted too; the optimizer
    /// removes them.
    pub fn t_count(&self) -> usize {
        self.operations
            .iter()
            .filter(|op| match op {
                PbcOp::Rotate { angle, .. } | PbcOp::ConditionalRotate { angle, .. } => {
                    angle.eighths() % 2 == 1
                }
                PbcOp::Measure { .. } => false,
            })
            .count()
    }

    /// Reduce the T count by merging and MCR-swapping rotations, preserving
    /// the circuit's quantum-classical channel. On error the circuit is
    /// unchanged.
    pub fn optimize_rotations(
        &mut self,
        options: OptimizeOptions,
    ) -> Result<OptimizeStats, PbcError> {
        let rotations = |c: &Self| {
            c.operations
                .iter()
                .filter(|op| matches!(op, PbcOp::Rotate { .. }))
                .count()
        };
        let mut stats = OptimizeStats {
            t_before: self.t_count(),
            rotations_before: rotations(self),
            ..OptimizeStats::default()
        };
        // Roots: every operation's axis, then the output frame's images.
        let n = self.num_qubits;
        let mut roots: Vec<PauliRef> = self.operations.iter().map(|op| op.axis().0).collect();
        for q in 0..n as u32 {
            roots.extend([self.output_frame.x(q), self.output_frame.z(q)]);
        }
        let (axes, records) = pack(&self.arena, n, &roots, options.max_packed_words)?;
        let mut items: Vec<Item> = self
            .operations
            .iter()
            .zip(&records)
            .enumerate()
            .filter_map(|(op, (operation, &(axis, sign)))| match operation {
                PbcOp::Rotate { angle, .. } => {
                    let k = normalize(i32::from(sign) * i32::from(angle.eighths()));
                    (k != 0 && !axes.is_identity(axis)).then(|| {
                        Item::Rot(Rot {
                            axis,
                            k,
                            support: axes.support(axis),
                        })
                    })
                }
                _ => Some(Item::Barrier { op, axis, sign }),
            })
            .collect();
        let mut frame: Vec<(u32, i8)> = records[self.operations.len()..].to_vec();

        let mut optimizer = Optimizer {
            axes,
            options,
            last: Vec::new(),
            touching: std::array::from_fn(|_| Vec::new()),
            scratch: Vec::new(),
            word: Vec::new(),
            stats: &mut stats,
        };
        if options.measure_depth {
            optimizer.stats.t_depth_before = optimizer.t_depth(&items);
        }
        // Counted over the original rotations, including any the packing
        // step dropped (identity axes weigh nothing anyway).
        optimizer.stats.weight_before = optimizer.weight(&items);
        let mut frame_changed = false;
        if options.strategy == Strategy::Litinski {
            frame_changed = optimizer.litinski(&mut items, &mut frame)?;
        }
        // Rounds of streaming merge (moving Cliffords into the frame as they
        // appear) and MCR swaps, until a round removes no T or the rounds run
        // out. Swaps can leave Clifford rotations, so a final stream moves
        // those into the frame too.
        let mut t = usize::MAX;
        let rounds = if options.strategy == Strategy::Merge {
            options.rounds.max(1)
        } else {
            0
        };
        for _ in 0..rounds {
            frame_changed |= optimizer.stream(&mut items, &mut frame, true)? > 0;
            let swaps = optimizer.swap_segments(&mut items);
            let now = t_count_items(&items);
            if now >= t && swaps == 0 {
                break;
            }
            t = now;
        }
        if options.strategy == Strategy::Merge
            && options.clifford_to_frame
            && items
                .iter()
                .any(|item| matches!(item, Item::Rot(r) if r.k % 2 == 0))
        {
            frame_changed |= optimizer.stream(&mut items, &mut frame, true)? > 0;
        }
        if options.measure_depth {
            optimizer.stats.t_depth_after = optimizer.t_depth(&items);
        }
        optimizer.stats.weight_after = optimizer.weight(&items);

        // Install: build handles for synthesized axes, sharing leaves.
        let mut axes = optimizer.axes;
        let mut leaves = vec![None; 3 * n];
        let mut handle = |arena: &mut PauliArena, id: u32| -> PauliRef {
            if let Some(h) = axes.handles[id as usize] {
                return h;
            }
            let h = build_handle(arena, &mut leaves, axes.get(id), axes.l);
            axes.handles[id as usize] = Some(h);
            h
        };
        let signed = |h: PauliRef, sign: i8| {
            if sign < 0 {
                h.scaled(Phase::MinusOne)
            } else {
                h
            }
        };
        let mut output = Vec::with_capacity(items.len());
        for item in items {
            output.push(match item {
                Item::Rot(rot) => PbcOp::Rotate {
                    axis: PauliAxis(handle(&mut self.arena, rot.axis)),
                    angle: PauliAngle::new(i64::from(rot.k)),
                },
                Item::Barrier { op, axis, sign } => {
                    let h = handle(&mut self.arena, axis);
                    match self.operations[op] {
                        PbcOp::Measure {
                            outcome, target, ..
                        } => PbcOp::Measure {
                            axis: PauliAxis(signed(h, sign)),
                            outcome,
                            target,
                        },
                        PbcOp::ConditionalRotate { angle, if_one, .. } => {
                            PbcOp::ConditionalRotate {
                                axis: PauliAxis(h),
                                angle: PauliAngle::new(
                                    i64::from(sign) * i64::from(angle.eighths()),
                                ),
                                if_one,
                            }
                        }
                        PbcOp::Rotate { .. } => unreachable!("rotations are not barriers"),
                    }
                }
            });
        }
        if frame_changed {
            for q in 0..n {
                let (x, sx) = frame[2 * q];
                let (z, sz) = frame[2 * q + 1];
                self.output_frame.x[q] = signed(handle(&mut self.arena, x), sx);
                self.output_frame.z[q] = signed(handle(&mut self.arena, z), sz);
            }
        }
        self.operations = output;
        stats.t_after = self.t_count();
        stats.rotations_after = rotations(self);
        Ok(stats)
    }
}

/// An arena expression for the canonical word `w`: the product of its
/// single-qubit factors, which is exactly `+W`. Leaves are shared through
/// `leaves` (index `3q + letter`).
fn build_handle(
    arena: &mut PauliArena,
    leaves: &mut [Option<PauliRef>],
    w: &[u64],
    l: usize,
) -> PauliRef {
    let (x, z) = w.split_at(l);
    let mut acc: Option<PauliRef> = None;
    for j in 0..l {
        let mut support = x[j] | z[j];
        while support != 0 {
            let bit = support.trailing_zeros() as usize;
            support &= support - 1;
            let q = 64 * j + bit;
            let (letter, pauli) = match (x[j] >> bit & 1, z[j] >> bit & 1) {
                (1, 0) => (0, Pauli::X),
                (1, 1) => (1, Pauli::Y),
                _ => (2, Pauli::Z),
            };
            let leaf = *leaves[3 * q + letter].get_or_insert_with(|| {
                arena.push(PauliNode::Single {
                    qubit: q as u32,
                    pauli,
                })
            });
            acc = Some(match acc {
                None => leaf,
                Some(acc) => arena.push(PauliNode::Product(acc, leaf)),
            });
        }
    }
    acc.unwrap_or_else(|| arena.reference(0))
}

const NONE: u32 = u32::MAX;

/// The indices of the set bits of `word`.
fn bits(mut word: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        (word != 0).then(|| {
            let bit = word.trailing_zeros() as usize;
            word &= word - 1;
            bit
        })
    })
}

fn t_count_items(items: &[Item]) -> usize {
    items
        .iter()
        .filter(|item| matches!(item, Item::Rot(r) if r.k % 2 != 0))
        .count()
}

fn t_count(rots: &[Rot]) -> usize {
    rots.iter().filter(|r| r.k % 2 != 0).count()
}

struct Optimizer<'a> {
    axes: Axes,
    options: OptimizeOptions,
    /// Per axis: position of its latest live rotation in the segment being
    /// merged, or NONE. Reset after every sweep.
    last: Vec<u32>,
    /// Per support-signature bit: positions, in increasing order, of the
    /// rotations in the segment being merged whose signature has that bit.
    /// Only these can block a merge of an axis with the bit. Reset per sweep.
    touching: [Vec<u32>; 64],
    /// Scratch for certificate products: packed words, then coefficients.
    scratch: Vec<u64>,
    /// Scratch for one product word.
    word: Vec<u64>,
    stats: &'a mut OptimizeStats,
}

impl Optimizer<'_> {
    /// One streaming pass over all operations, in order. With
    /// `clifford_to_frame`, keep the accumulated Clifford F (initially the
    /// identity) that is being moved to the end: conjugate each axis by F,
    /// and absorb every unconditional Clifford rotation into F, whether given
    /// or produced by a merge. A merged rotation at position j commutes with
    /// everything after it (that is what allowed the merge), so absorbing it
    /// from there is sound. Each rotation merges into the latest earlier
    /// rotation with the same axis that it can reach. Barriers end merging
    /// but not F. Finally the output frame's images are conjugated by F.
    /// Returns the number of Clifford rotations absorbed.
    fn stream(
        &mut self,
        items: &mut Vec<Item>,
        frame: &mut [(u32, i8)],
        merge: bool,
    ) -> Result<usize, PbcError> {
        let num_qubits = frame.len() / 2;
        let l = self.axes.l;
        let to_frame = self.options.clifford_to_frame;
        let mut f = Clifford::identity(if to_frame { num_qubits } else { 0 }, l, Axes::signature);
        let mut word = vec![0u64; 2 * l];
        // F† (sign W) F as (id, sign).
        let mut conjugate = |axes: &mut Axes, f: &mut Clifford, id: u32, sign: i8| {
            if !to_frame || f.fixes(axes.get(id)) {
                return (id, sign);
            }
            let phase = f.conjugate(axes.get(id), &mut word);
            debug_assert!(phase.is_multiple_of(2), "conjugation preserves Hermiticity");
            let new = axes.intern(&word, None);
            (new, if phase == 0 { sign } else { -sign })
        };
        let (mut merges, mut absorbed) = (0, 0);
        let mut out: Vec<Item> = Vec::with_capacity(items.len());
        let mut segment_start = 0;
        // With `lazy_cliffords`, Clifford rotations stay in `out` as merge
        // targets until something forces them into F, so a later rotation can
        // still merge into them (T + T + T stays one rotation instead of S in
        // F plus a T). Positions; entries go stale once the angle turns odd or
        // zero. Without it, Cliffords go into F at once and this stays empty.
        let lazy = self.options.lazy_cliffords;
        let mut pending: Vec<u32> = Vec::new();
        for &item in items.iter() {
            match item {
                Item::Rot(rot) => {
                    let (mut axis, mut sign) = conjugate(&mut self.axes, &mut f, rot.axis, 1);
                    self.last.resize(self.axes.len(), NONE);
                    // A pending Clifford that anticommutes with this rotation
                    // cannot stay in front of it: move it into F first.
                    if to_frame && !pending.is_empty() {
                        let blocking = Some((axis, self.axes.support(axis)));
                        let n = self.absorb_pending(&mut out, &mut pending, &mut f, blocking);
                        if n > 0 {
                            absorbed += n;
                            (axis, sign) = conjugate(&mut self.axes, &mut f, rot.axis, 1);
                            self.last.resize(self.axes.len(), NONE);
                        }
                    }
                    let k = normalize(i32::from(sign) * i32::from(rot.k));
                    if to_frame && !lazy && k % 2 == 0 {
                        f.absorb(self.axes.get(axis), self.axes.support(axis), k);
                        absorbed += 1;
                        continue;
                    }
                    let rot = Rot {
                        axis,
                        k,
                        support: self.axes.support(axis),
                    };
                    let j = if merge {
                        self.last[axis as usize]
                    } else {
                        NONE
                    };
                    if j != NONE && self.can_reach(&out, j as usize, rot) {
                        merges += 1;
                        let Item::Rot(target) = &mut out[j as usize] else {
                            unreachable!("merge targets are rotations")
                        };
                        target.k = normalize(i32::from(target.k) + i32::from(rot.k));
                        if to_frame && target.k != 0 && target.k % 2 == 0 {
                            if lazy {
                                pending.push(j);
                            } else {
                                f.absorb(self.axes.get(axis), rot.support, target.k);
                                absorbed += 1;
                                target.k = 0;
                            }
                        }
                        if target.k == 0 {
                            self.last[axis as usize] = NONE;
                        }
                        continue;
                    }
                    let position = out.len() as u32;
                    self.last[axis as usize] = position;
                    for bit in bits(rot.support) {
                        self.touching[bit].push(position);
                    }
                    if to_frame && rot.k % 2 == 0 {
                        pending.push(position);
                    }
                    out.push(Item::Rot(rot));
                }
                Item::Barrier { op, axis, sign } => {
                    absorbed += self.absorb_pending(&mut out, &mut pending, &mut f, None);
                    let (axis, sign) = conjugate(&mut self.axes, &mut f, axis, sign);
                    self.end_segment(&out[segment_start..]);
                    out.push(Item::Barrier { op, axis, sign });
                    segment_start = out.len();
                }
            }
        }
        absorbed += self.absorb_pending(&mut out, &mut pending, &mut f, None);
        self.end_segment(&out[segment_start..]);
        if absorbed > 0 {
            for image in frame.iter_mut() {
                *image = conjugate(&mut self.axes, &mut f, image.0, image.1);
            }
        }
        if self.axes.words_len() > self.options.max_packed_words {
            return Err(PbcError::ExpansionLimit);
        }
        out.retain(|item| !matches!(item, Item::Rot(r) if r.k == 0));
        *items = out;
        self.stats.merges += merges;
        self.stats.cliffords_to_frame += absorbed;
        Ok(absorbed)
    }

    /// Move pending Clifford rotations from `out` into F: all of them, or
    /// only those anticommuting with `blocking` (an axis and its support).
    /// Sound because each pending Clifford commutes with everything after it
    /// in `out`: any rotation pushed after it that anticommuted would have
    /// absorbed it first. Returns the number absorbed.
    fn absorb_pending(
        &mut self,
        out: &mut [Item],
        pending: &mut Vec<u32>,
        f: &mut Clifford,
        blocking: Option<(u32, u64)>,
    ) -> usize {
        let mut absorbed = 0;
        let mut keep = 0;
        for i in 0..pending.len() {
            let p = pending[i];
            let Item::Rot(r) = &mut out[p as usize] else {
                unreachable!("pending entries are rotations")
            };
            if r.k % 2 != 0 || r.k == 0 {
                continue; // stale: merged into an odd angle, or cancelled
            }
            if let Some((axis, support)) = blocking
                && (r.support & support == 0 || !self.axes.anticommute(r.axis, axis))
            {
                pending[keep] = p;
                keep += 1;
                continue;
            }
            f.absorb(self.axes.get(r.axis), r.support, r.k);
            self.last[r.axis as usize] = NONE;
            r.k = 0;
            absorbed += 1;
        }
        pending.truncate(keep);
        absorbed
    }

    /// Total Pauli weight of the rotations' axes.
    fn weight(&self, items: &[Item]) -> usize {
        items
            .iter()
            .map(|item| match item {
                Item::Rot(rot) => self.axes.weight(rot.axis),
                Item::Barrier { .. } => 0,
            })
            .sum()
    }

    /// Reset the merge state for a finished segment.
    fn end_segment(&mut self, segment: &[Item]) {
        for item in segment {
            if let Item::Rot(rot) = item {
                self.last[rot.axis as usize] = NONE;
            }
        }
        for list in &mut self.touching {
            list.clear();
        }
    }

    /// MCR swaps within each run of rotations between barriers.
    fn swap_segments(&mut self, items: &mut Vec<Item>) -> usize {
        if self.options.window == 0 {
            return 0;
        }
        let mut swaps = 0;
        let mut output = Vec::with_capacity(items.len());
        let mut segment = Vec::new();
        for &item in items.iter() {
            match item {
                Item::Rot(rot) => segment.push(rot),
                Item::Barrier { .. } => {
                    swaps += self.swaps(&mut segment);
                    output.extend(segment.drain(..).map(Item::Rot));
                    output.push(item);
                }
            }
        }
        swaps += self.swaps(&mut segment);
        output.extend(segment.drain(..).map(Item::Rot));
        *items = output;
        swaps
    }

    /// Whether `rot`'s axis commutes with every live rotation after position
    /// `j`, within the lookback limit. A blocker must share a qubit, and so a
    /// signature bit, with the axis, so only those bits' position lists are
    /// visited, newest first, where blockers are likeliest.
    fn can_reach(&self, out: &[Item], j: usize, rot: Rot) -> bool {
        out.len() - j - 1 <= self.options.lookback
            && bits(rot.support).all(|bit| {
                self.touching[bit]
                    .iter()
                    .rev()
                    .take_while(|&&p| p as usize > j)
                    .all(|&p| match out[p as usize] {
                        Item::Rot(r) => r.k == 0 || !self.axes.anticommute(r.axis, rot.axis),
                        Item::Barrier { .. } => unreachable!("barriers end segments"),
                    })
            })
    }

    /// End of the greedy internally commuting run starting at `start`.
    fn run_end(&self, rots: &[Rot], start: usize) -> usize {
        let mut end = start + 1;
        while end < rots.len()
            && end - start < self.options.window
            && rots[start..end]
                .iter()
                .all(|r| !self.axes.anticommute(r.axis, rots[end].axis))
        {
            end += 1;
        }
        end
    }

    /// One pass of MCR swaps over the segment; returns the number accepted.
    /// Accepted rewrites cover disjoint, increasing ranges and are applied in
    /// one rebuild; the next round sees their results.
    fn swaps(&mut self, rots: &mut Vec<Rot>) -> usize {
        let mut edits = Vec::new();
        let mut start = 0;
        while start < rots.len() {
            let a = self.run_end(rots, start);
            if a == rots.len() {
                break;
            }
            let b = self.run_end(rots, a);
            if b == rots.len() {
                break;
            }
            let c = self.run_end(rots, b);
            match self.try_triple(rots, start, a, b, c) {
                Some((range, trial)) => {
                    start = range.end;
                    edits.push((range, trial));
                }
                None => start = a,
            }
        }
        if edits.is_empty() {
            return 0;
        }
        let accepted = edits.len();
        let mut out = Vec::with_capacity(rots.len());
        let mut next = 0;
        for (range, trial) in edits {
            out.extend_from_slice(&rots[next..range.start]);
            out.extend(trial);
            next = range.end;
        }
        out.extend_from_slice(&rots[next..]);
        *rots = out;
        self.stats.swaps += accepted;
        accepted
    }

    /// Try swaps for runs A = [start, a), B = [a, b), C = [b, c): first the
    /// full runs, then up to `candidates` variants with A shortened to a
    /// suffix or C to a prefix. Returns the replaced range and its rewrite for
    /// the first trial that lowers T.
    fn try_triple(
        &mut self,
        rots: &[Rot],
        start: usize,
        a: usize,
        b: usize,
        c: usize,
    ) -> Option<(std::ops::Range<usize>, Vec<Rot>)> {
        let variants = (start..a)
            .flat_map(|s| (b + 1..=c).rev().map(move |e| (s, e)))
            .take(1 + self.options.candidates);
        // Certificates depend only on A (for A,B) or only on C (for B,C), so
        // each is computed at most once across the variants.
        let mut ab: Vec<Option<bool>> = vec![None; a - start];
        let mut bc: Vec<Option<bool>> = vec![None; c - b];
        for (s, e) in variants {
            let (ga, gb, gc) = (&rots[s..a], &rots[a..b], &rots[b..e]);
            // A swap helps only by bringing A and C together. Both groups
            // commute internally, so after either swap every axis they share
            // merges. That removes T only if both angles are odd; eagerly, any
            // merge is worth it (a later round or the frame may profit).
            let eager = self.options.eager_swaps;
            let useful = |r: &&Rot| eager || r.k % 2 != 0;
            if !ga
                .iter()
                .filter(useful)
                .any(|x| gc.iter().filter(useful).any(|y| x.axis == y.axis))
            {
                continue;
            }
            let baseline = &rots[s..e];
            if *ab[s - start].get_or_insert_with(|| self.generators_commute(ga, gb))
                && let Some(trial) = self.improves([gb, ga, gc], baseline)
            {
                return Some((s..e, trial));
            }
            if *bc[e - b - 1].get_or_insert_with(|| self.generators_commute(gb, gc))
                && let Some(trial) = self.improves([ga, gc, gb], baseline)
            {
                return Some((s..e, trial));
            }
        }
        None
    }

    /// The groups concatenated and locally merged, if that beats `baseline`:
    /// fewer T, or, with eager swaps, no more T and fewer rotations. Either
    /// way each accepted swap strictly shrinks the circuit.
    fn improves(&self, groups: [&[Rot]; 3], baseline: &[Rot]) -> Option<Vec<Rot>> {
        let mut trial: Vec<Rot> = groups.concat();
        for i in 0..trial.len() {
            for j in (0..i).rev() {
                if trial[j].k == 0 {
                    continue;
                }
                if trial[j].axis == trial[i].axis {
                    trial[j].k = normalize(i32::from(trial[j].k) + i32::from(trial[i].k));
                    trial[i].k = 0;
                    break;
                }
                if self.axes.anticommute(trial[j].axis, trial[i].axis) {
                    break;
                }
            }
        }
        trial.retain(|r| r.k != 0);
        let (t, t_baseline) = (t_count(&trial), t_count(baseline));
        let better = t < t_baseline
            || (self.options.eager_swaps && t <= t_baseline && trial.len() < baseline.len());
        better.then_some(trial)
    }

    /// The MCR certificate for internally commuting groups A and B:
    /// `[H_A, H_B] = 0` for `H = sum k P`. Each anticommuting pair with
    /// `P Q = i^e S` contributes `+-k_P k_Q` to the coefficient of S (+ for
    /// e = 1, - for e = 3); the generators commute iff every total is zero.
    fn generators_commute(&mut self, a: &[Rot], b: &[Rot]) -> bool {
        let l = self.axes.l;
        let stride = 2 * l;
        // Entries: packed word followed by its coefficient, stride + 1 u64s.
        self.scratch.clear();
        let mut word = std::mem::take(&mut self.word);
        word.resize(stride, 0);
        for p in a {
            for q in b {
                let (pw, qw) = (self.axes.get(p.axis), self.axes.get(q.axis));
                if !anticommutes(pw, qw, l) {
                    continue;
                }
                let e = product(pw, qw, &mut word, l);
                let sign = if e == 1 { 1 } else { -1 };
                let delta = sign * i64::from(p.k) * i64::from(q.k);
                match self
                    .scratch
                    .chunks_exact_mut(stride + 1)
                    .find(|entry| entry[..stride] == word[..])
                {
                    Some(entry) => entry[stride] = (entry[stride] as i64 + delta) as u64,
                    None => {
                        self.scratch.extend_from_slice(&word);
                        self.scratch.push(delta as u64);
                    }
                }
            }
        }
        self.word = word;
        self.scratch
            .chunks_exact(stride + 1)
            .all(|entry| entry[stride] == 0)
    }
}

#[cfg(test)]
mod tests;
