//! CNOT minimization by phase-polynomial resynthesis.
//!
//! A port of Feynman's `-cnotmin` pass (`minCNOT`, the algorithm of Amy,
//! Azimzadeh and Mosca, *Q. Sci. Tech.* 2017) into tzap's pass pipeline.
//!
//! # What it does
//!
//! Between the gates it cannot interpret, a circuit is *CNOT-dihedral*: every
//! gate is either a CNOT/X, which permutes computational-basis states, or a
//! diagonal rotation. Such a block is fully described by two things, and
//! nothing else:
//!
//! - the **linear map** it applies to the qubits (each qubit ends up holding
//!   some XOR of the values the qubits held on entry), and
//! - the **phase polynomial** it applies: a set of (parity, angle) pairs, each
//!   meaning "rotate by `angle` on the subspace where this XOR of inputs is 1".
//!
//! Two blocks with the same pair are the same operator, whatever gates they
//! are written with. So a block can be thrown away and *re-synthesized* from
//! its (linear map, phase polynomial) pair — and if the synthesis is cleverer
//! than what the block happened to be written as, the rewrite is a win. That
//! is the whole pass: recover the pair, re-synthesize, keep whichever of the
//! two is smaller.
//!
//! # Structure
//!
//! 1. A forward scan ([`Chunk::feed`]) tracks each qubit's parity as a
//!    bitvector over the values the qubits held at the start of the current
//!    block, folding rotations into a phase map keyed by parity.
//! 2. A gate outside the fragment (H, ccx, ccz, measure, reset) ends the
//!    block. The block is synthesized, the gate is emitted, and a fresh block
//!    starts with the basis renormalized so that every qubit is once again its
//!    own variable — which is what keeps parities `n`-dimensional forever
//!    rather than growing a variable per Hadamard.
//! 3. Synthesis ([`synthesize`]) is Gray-code phase synthesis
//!    ([`gray_synth`]) for the rotations, then Gaussian elimination
//!    ([`linear_synth`]) to fix up whatever linear map that left behind into
//!    the one the block actually needs.
//!
//! # Monotonicity
//!
//! Re-synthesis is not always a win: on a block that was already written
//! optimally, Gaussian elimination can easily emit more CNOTs than it
//! replaces. Every block is therefore synthesized *and* compared against the
//! gates it came from, and the original is kept unless the replacement is
//! strictly better with no more gates overall (see [`Chunk::flush`]). The pass
//! can consequently never make a circuit longer, which is what lets it sit in
//! a fixpoint loop without oscillating.

use std::f64::consts::PI;

use rustc_hash::FxHashMap;

use crate::circuit::{Circuit, Gate, Qubit};
use crate::pass::Pass;

/// A parity: bit `i` is set when local qubit `i`'s starting value is XOR-ed
/// into it. One `u128` bounds a block at [`MAX_CHUNK_QUBITS`] qubits, which is
/// what makes parities `Copy` and hashable without allocating.
type Parity = u128;

/// Largest number of distinct qubits one block may span. A block that would
/// outgrow this is flushed and a fresh one started, so the cap costs
/// optimization across the split but never correctness.
///
/// This is also the pass's main speed control. Synthesis is superlinear in
/// the block's qubit count (Gaussian elimination is cubic in it, and
/// [`gray_synth`]'s basis rewrite is charged once per emitted CNOT), so
/// halving this cap cuts the worst-case cost of a block by roughly eight.
/// 128 is the ceiling [`Parity`] can represent at all.
///
/// Wide blocks are where the pass spends its time and not where it wins.
/// Measured over all 8.67 million blocks built across `benchmarks/` at -O3,
/// blocks of 49 qubits and wider were 75% of total synthesis time and won
/// zero times, and a sweep of this cap from 32 to 128 left output quality
/// unchanged while runtime grew monotonically (1.17x the no-CnotMin runtime
/// at 32, against 1.47x at 128). [`Budget`] is what keeps that tail
/// affordable at the full width: it abandons a synthesis as soon as the
/// result can no longer be kept, so a wide block costs its own size rather
/// than the size of the circuit Gauss-Jordan would have produced for it.
pub const MAX_CHUNK_QUBITS: usize = 128;

/// Largest number of distinct rotation parities one block may carry, for the
/// same reason as [`MAX_CHUNK_QUBITS`]: [`gray_synth`] splits the parity set
/// recursively and rewrites the surviving parities on every CNOT it emits, so
/// its cost grows faster than linearly in this.
///
/// 512 was measured, not guessed. Raising it to 4096 leaves the output
/// byte-identical on all 135 benchmarks in `evaluation/comparison/bench` --
/// no real block here carries anything like that many distinct parities --
/// while costing 3.5x on the adversarial shape the cap exists for (a single
/// huge block whose rotations sit on thousands of unrelated parities, where
/// synthesis always loses to the original anyway; see the `term_cap_sweep`
/// timing test).
pub const MAX_CHUNK_TERMS: usize = 512;

/// Rotations closer to a multiple of 2π than this are dropped as identity,
/// and angles within it of a multiple of π/4 are emitted as the corresponding
/// Clifford+T gate rather than as an `rz`. Matches the tolerance the phase
/// folding passes use, so a rotation that survives one survives the other.
const ANGLE_TOL: f64 = 1e-9;

/// Replaces CNOT-dihedral blocks with re-synthesized equivalents that use
/// fewer two-qubit gates. See the module documentation.
pub struct CnotMin {
    /// See [`MAX_CHUNK_QUBITS`].
    pub max_qubits: usize,
    /// See [`MAX_CHUNK_TERMS`].
    pub max_terms: usize,
}

impl Default for CnotMin {
    fn default() -> Self {
        CnotMin {
            max_qubits: MAX_CHUNK_QUBITS,
            max_terms: MAX_CHUNK_TERMS,
        }
    }
}

impl Pass for CnotMin {
    fn name(&self) -> &str {
        "CNOT minimization"
    }
    fn run(&self, circuit: &Circuit) -> Circuit {
        cnot_min_with(circuit, self.max_qubits, self.max_terms)
    }
}

/// Free-function form of [`CnotMin`] at its default bounds.
pub fn cnot_min(circuit: &Circuit) -> Circuit {
    cnot_min_with(circuit, MAX_CHUNK_QUBITS, MAX_CHUNK_TERMS)
}

/// [`cnot_min`] with explicit block bounds. `max_qubits` is clamped to what a
/// [`Parity`] can address, and both bounds to at least one, since a block that
/// can hold nothing at all would never make progress.
pub fn cnot_min_with(circuit: &Circuit, max_qubits: usize, max_terms: usize) -> Circuit {
    let max_qubits = max_qubits.clamp(1, MAX_CHUNK_QUBITS);
    let max_terms = max_terms.max(1);
    let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
    let mut chunk = Chunk::new(circuit.num_qubits, max_qubits, max_terms);

    for gate in &circuit.gates {
        if let Some(passthrough) = chunk.feed(gate, &mut output) {
            // Outside the CNOT-dihedral fragment: the block ends here, and the
            // gate itself is copied over untouched.
            output.apply(passthrough.clone());
        }
    }
    chunk.flush(&mut output);
    output
}

/// [`cnot_min`] with the synthesis budget disabled: every block is
/// synthesized to completion and only then compared. Slow on wide blocks --
/// which is the entire reason the budget exists -- but it must produce
/// byte-identical output, which is what
/// `budget_is_a_pure_early_exit` checks.
#[cfg(test)]
fn cnot_min_unbounded(circuit: &Circuit) -> Circuit {
    let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
    let mut chunk = Chunk::new(circuit.num_qubits, MAX_CHUNK_QUBITS, MAX_CHUNK_TERMS);
    chunk.bounded = false;
    for gate in &circuit.gates {
        if let Some(passthrough) = chunk.feed(gate, &mut output) {
            output.apply(passthrough.clone());
        }
    }
    chunk.flush(&mut output);
    output
}

/// One CNOT-dihedral block under construction.
///
/// Local indices number only the qubits the block has actually touched, so a
/// block confined to three wires of a 200-qubit circuit is synthesized as a
/// three-qubit problem. `parity[i]` is a XOR of the values the block's qubits
/// held on entry, and `consts[i]` records an odd number of X gates on top.
struct Chunk {
    max_qubits: usize,
    max_terms: usize,
    /// Local index -> circuit qubit.
    qubits: Vec<Qubit>,
    /// Circuit qubit -> local index, or `usize::MAX` when untouched.
    local: Vec<usize>,
    parity: Vec<Parity>,
    consts: Vec<bool>,
    /// Rotation angle per parity. A parity is never zero (the qubit parities
    /// stay a basis, see [`Chunk::renormalize`]), so no entry here is a pure
    /// global phase.
    terms: FxHashMap<Parity, f64>,
    /// The gates this block was built from, kept so [`Chunk::flush`] can fall
    /// back to them when synthesis does no better.
    original: Vec<Gate>,
    /// Whether [`Chunk::flush`] gives synthesis a [`Budget`] bounded by the
    /// original block. Always `true` in production; the tests turn it off to
    /// check that the bound is a pure early exit -- that abandoning a losing
    /// synthesis leaves exactly the output that finishing and rejecting it
    /// would have.
    bounded: bool,
    /// Buffers handed to synthesis and taken back afterwards, so that a pass
    /// over a circuit with millions of blocks allocates once rather than once
    /// per block. Purely a cache: their contents never carry across a flush.
    scratch: Scratch,
}

/// Reusable allocations for [`Chunk::flush`]. Every field is cleared before
/// use, so nothing here is state -- it exists only to keep synthesis off the
/// allocator.
#[derive(Default)]
struct Scratch {
    /// [`Budget`]'s gate buffer.
    gates: Vec<Gate>,
    /// Retired parity-set buffers, for [`gray_synth`] to hand back out.
    pool: Pool,
    /// [`gray_synth`]'s explicit recursion stack.
    stack: Vec<Pt>,
    /// The parity each local qubit holds as [`gray_synth`] proceeds.
    state: Vec<Parity>,
    /// [`invert`]'s working copy of the matrix being reduced, and the inverse
    /// it accumulates.
    elim: Vec<Parity>,
    inverse: Vec<Parity>,
    /// [`linear_synth`]'s residual map and the row operations reducing it.
    m: Vec<Parity>,
    ops: Vec<(usize, usize)>,
}

impl Chunk {
    fn new(num_qubits: usize, max_qubits: usize, max_terms: usize) -> Self {
        Chunk {
            max_qubits,
            max_terms,
            qubits: Vec::new(),
            local: vec![usize::MAX; num_qubits],
            parity: Vec::new(),
            consts: Vec::new(),
            terms: FxHashMap::default(),
            original: Vec::new(),
            bounded: true,
            scratch: Scratch::default(),
        }
    }

    /// The local index for `q`, which [`Chunk::admits`] must already have
    /// confirmed there is room for.
    fn index_of(&mut self, q: Qubit) -> usize {
        let existing = self.local[q as usize];
        if existing != usize::MAX {
            return existing;
        }
        let i = self.qubits.len();
        self.qubits.push(q);
        self.local[q as usize] = i;
        // A qubit joins holding its own starting value, by definition of the
        // basis this block is expressed in.
        self.parity.push(1u128 << i);
        self.consts.push(false);
        i
    }

    /// Whether `operands` fit under the qubit cap without evicting anything.
    /// Checked before any mutation, so a gate whose second operand does not
    /// fit cannot leave the first half-admitted.
    fn admits(&self, operands: &[Qubit]) -> bool {
        let mut fresh = 0;
        for (i, &q) in operands.iter().enumerate() {
            if self.local[q as usize] == usize::MAX && !operands[..i].contains(&q) {
                fresh += 1;
            }
        }
        self.qubits.len() + fresh <= self.max_qubits
    }

    /// Absorb `gate`. Returns `Some(gate)` when the caller should emit it
    /// verbatim instead: either it is outside the CNOT-dihedral fragment (the
    /// block having been flushed first), or it is wider than the qubit cap
    /// allows any block to be.
    fn feed<'g>(&mut self, gate: &'g Gate, output: &mut Circuit) -> Option<&'g Gate> {
        let (rotation, operands): (Option<f64>, &[Qubit]) = match gate {
            Gate::t(q) => (Some(PI / 4.0), std::slice::from_ref(q)),
            Gate::tdg(q) => (Some(-PI / 4.0), std::slice::from_ref(q)),
            Gate::s(q) => (Some(PI / 2.0), std::slice::from_ref(q)),
            Gate::sdg(q) => (Some(-PI / 2.0), std::slice::from_ref(q)),
            Gate::z(q) => (Some(PI), std::slice::from_ref(q)),
            Gate::rz(theta, q) => (Some(*theta), std::slice::from_ref(q)),
            Gate::x(q) => (None, std::slice::from_ref(q)),
            // A degenerate two-qubit gate (control == target) is not
            // well-formed, but the corpora contain them, and interpreting one
            // would zero a qubit's parity -- destroying the basis this pass's
            // whole representation rests on, and with it the invertibility
            // `linear_synth` needs. Treat it the way an uninterpretable gate
            // is treated: end the block and copy it through untouched, so
            // whatever the rest of tzap takes it to mean is preserved exactly.
            Gate::cnot { control, target } | Gate::cz { control, target } if control == target => {
                self.flush(output);
                return Some(gate);
            }
            Gate::cnot { control, target } | Gate::cz { control, target } => {
                (None, &[*control, *target])
            }
            // Everything else leaves the fragment. `ccz` is diagonal but its
            // phase is cubic in the inputs, not a parity, so it gets the same
            // treatment as the genuinely non-diagonal gates.
            Gate::h(_)
            | Gate::ccx { .. }
            | Gate::ccz { .. }
            | Gate::measure { .. }
            | Gate::reset(_) => {
                self.flush(output);
                return Some(gate);
            }
        };

        if !self.admits(operands) {
            // An empty block that still cannot take the gate never will:
            // flushing again would free nothing and loop forever. Let it
            // through untouched, which is sound precisely because there is no
            // block in progress for it to be reordered against.
            if self.qubits.is_empty() {
                return Some(gate);
            }
            self.flush(output);
            return self.feed(gate, output);
        }

        // The parity a rotation would land on has to be known before the term
        // cap can be checked, and the check has to happen before the block is
        // mutated -- hence resolving the index first and only then deciding.
        if let Some(angle) = rotation {
            let i = self.index_of(operands[0]);
            let (p, k) = (self.parity[i], self.consts[i]);
            if !self.terms.contains_key(&p) && self.terms.len() >= self.max_terms {
                self.flush(output);
                return self.feed(gate, output);
            }
            self.add_term(p, k, angle);
            self.original.push(gate.clone());
            return None;
        }

        match gate {
            Gate::x(q) => {
                let i = self.index_of(*q);
                self.consts[i] = !self.consts[i];
            }
            Gate::cnot { control, target } => {
                let (c, t) = (self.index_of(*control), self.index_of(*target));
                self.parity[t] ^= self.parity[c];
                self.consts[t] ^= self.consts[c];
            }
            // CZ is diagonal, so it contributes phase without disturbing any
            // parity: pi on the AND of the two wires, which as a polynomial
            // over parities is pi/2 on each and -pi/2 on their XOR.
            Gate::cz { control, target } => {
                let (c, t) = (self.index_of(*control), self.index_of(*target));
                let (pc, pt) = (self.parity[c], self.parity[t]);
                let (kc, kt) = (self.consts[c], self.consts[t]);
                self.add_term(pc, kc, PI / 2.0);
                self.add_term(pt, kt, PI / 2.0);
                self.add_term(pc ^ pt, kc ^ kt, -PI / 2.0);
            }
            _ => unreachable!("non-rotation gate kinds are exhausted above"),
        }
        self.original.push(gate.clone());
        None
    }

    /// Fold a rotation of `angle` on parity `p` (complemented when `k`) into
    /// the phase map.
    ///
    /// A rotation sitting on the complement of `p` is the same as one on `p`
    /// with the opposite angle, plus a global phase — which is dropped, since
    /// everything downstream compares circuits up to global phase.
    fn add_term(&mut self, p: Parity, k: bool, angle: f64) {
        let signed = if k { -angle } else { angle };
        *self.terms.entry(p).or_insert(0.0) += signed;
    }

    /// Emit the block into `output` and start a fresh one.
    ///
    /// The synthesized form replaces the original only when it uses fewer
    /// two-qubit gates and no more gates in total, so neither count can ever
    /// rise. Either way the qubits end in the same state, so the next block
    /// starts from the same place regardless of which was chosen.
    fn flush(&mut self, output: &mut Circuit) {
        if self.original.is_empty() {
            self.reset();
            return;
        }
        let (old_2q, old_all) = (count_2q(&self.original), self.original.len());
        // The budget bounds synthesis by the original's own size. It never
        // decides anything: a replacement that grows past the original in
        // either count has already lost the comparison below, so abandoning
        // it there yields the same output as finishing and rejecting it --
        // just without doing the work. (Equality is left in bounds, since
        // matching the two-qubit count still wins on a lower total.)
        let (max_two_q, max_gates) = match self.bounded {
            true => (old_2q, old_all),
            false => (usize::MAX, usize::MAX),
        };
        let mut budget = Budget::with_buffer(
            std::mem::take(&mut self.scratch.gates),
            max_two_q,
            max_gates,
        );
        // `push` stops one gate past `max_gates`, so a bounded synthesis never
        // needs more room than this -- and the buffer is reused, so this
        // reserve is a no-op after the first few blocks.
        budget.gates.reserve(old_all.saturating_add(1));
        let complete = synthesize(
            &self.parity,
            &self.consts,
            &self.terms,
            &self.qubits,
            &mut budget,
            &mut self.scratch,
        );
        let (new_2q, new_all) = (budget.two_q, budget.gates.len());
        let keep = complete && new_all <= old_all && (new_2q, new_all) < (old_2q, old_all);
        let mut buffer = budget.gates;
        for gate in if keep { &buffer } else { &self.original } {
            output.apply(gate.clone());
        }
        buffer.clear();
        self.scratch.gates = buffer;
        self.reset();
    }

    /// Renormalize: forget the block, and let the state the qubits are now in
    /// be the basis the next block is expressed over.
    ///
    /// This is what bounds the representation. Feynman grows the variable set
    /// at every Hadamard and carries parities in the widening space; because
    /// the qubit parities are always a full-rank basis here (X gates do not
    /// touch the linear part, and CNOT is invertible), the state can instead
    /// be renamed back to one variable per qubit at every block boundary, and
    /// a parity never needs more bits than the circuit has qubits.
    fn reset(&mut self) {
        for &q in &self.qubits {
            self.local[q as usize] = usize::MAX;
        }
        self.qubits.clear();
        self.parity.clear();
        self.consts.clear();
        self.terms.clear();
        self.original.clear();
    }
}

/// A gate buffer for a synthesized block that refuses to grow past the size
/// of the block it would replace.
///
/// [`Chunk::flush`] keeps the synthesized form only when it uses no more
/// gates overall and lexicographically fewer `(two-qubit, total)`. Both
/// counts only ever grow as synthesis proceeds, so once either has passed the
/// original's, no continuation can win, and every gate emitted after that
/// point is thrown away. Stopping there is *exactly* as correct as running to
/// completion and losing the comparison -- `budget_is_a_pure_early_exit`
/// pins that down -- but it bounds the work by the size of the block rather
/// than by the size of what synthesis would have produced.
///
/// That gap is the pass's whole cost on wide blocks: [`linear_synth`] reduces
/// the block's `n`x`n` map by Gauss-Jordan and emits a CNOT per row
/// operation, which for a dense map is order `n^2` however few CNOTs the
/// block was actually written with. Measured on `gf2^128_mult`, blocks of
/// 65-128 qubits were 97% of the pass's runtime and won zero times, the worst
/// of them expanding 632 gates into 12,286 before being discarded.
struct Budget {
    gates: Vec<Gate>,
    two_q: usize,
    /// Ceilings from the original block. Exceeding either decides the
    /// comparison against the replacement, so synthesis stops.
    max_two_q: usize,
    max_gates: usize,
}

impl Budget {
    #[cfg(test)]
    fn new(max_two_q: usize, max_gates: usize) -> Self {
        Budget::with_buffer(Vec::new(), max_two_q, max_gates)
    }

    /// [`Budget::new`], reusing `gates` as the buffer. The buffer is emptied
    /// first, so this differs from `new` only in the allocation it avoids.
    fn with_buffer(mut gates: Vec<Gate>, max_two_q: usize, max_gates: usize) -> Self {
        gates.clear();
        Budget {
            gates,
            two_q: 0,
            max_two_q,
            max_gates,
        }
    }

    /// Append `gate`, returning `false` once the budget is spent -- at which
    /// point the caller must abandon synthesis, since nothing it emits from
    /// here on can be kept.
    #[must_use]
    fn push(&mut self, gate: Gate) -> bool {
        if matches!(gate, Gate::cnot { .. } | Gate::cz { .. }) {
            self.two_q += 1;
            if self.two_q > self.max_two_q {
                return false;
            }
        }
        self.gates.push(gate);
        self.gates.len() <= self.max_gates
    }

    /// Whether `extra` more two-qubit gates would still fit, so a caller that
    /// knows in advance how many it owes can give up before emitting them.
    fn admits_2q(&self, extra: usize) -> bool {
        self.two_q + extra <= self.max_two_q && self.gates.len() + extra <= self.max_gates
    }
}

fn count_2q(gates: &[Gate]) -> usize {
    gates
        .iter()
        .filter(|g| matches!(g, Gate::cnot { .. } | Gate::cz { .. }))
        .count()
}

/// Synthesize a block: rotations by Gray-code synthesis, then a linear fix-up
/// onto the map the block actually applies, then the X gates its affine part
/// calls for.
///
/// `parity`/`consts` describe the block's target state in local indices,
/// `terms` its phase polynomial, and `qubits` maps local indices back to
/// circuit qubits.
fn synthesize(
    parity: &[Parity],
    consts: &[bool],
    terms: &FxHashMap<Parity, f64>,
    qubits: &[Qubit],
    budget: &mut Budget,
    scratch: &mut Scratch,
) -> bool {
    let n = qubits.len();
    let mut phases = scratch.pool.pop().unwrap_or_default();
    phases.extend(
        terms
            .iter()
            .filter(|&(_, &a)| !angle_is_zero(a))
            .map(|(&p, &a)| (p, a)),
    );
    // Iteration order of a hash map is not deterministic across runs; the
    // synthesized circuit must be.
    phases.sort_unstable_by_key(|&(p, _)| p);

    // `state[i]` is the parity local qubit `i` holds as synthesis proceeds,
    // starting from the basis the block is expressed in. Held outside
    // `scratch` for the call so that synthesis can borrow the rest of it.
    let mut state = std::mem::take(&mut scratch.state);
    state.clear();
    state.extend((0..n).map(|i| (1 as Parity) << i));
    let synthesized = gray_synth(n, phases, &mut state, budget, qubits, scratch)
        && linear_synth(&state, parity, qubits, budget, scratch);
    scratch.state = state;
    if !synthesized {
        return false;
    }
    for (i, &k) in consts.iter().enumerate() {
        if k && !budget.push(Gate::x(qubits[i])) {
            return false;
        }
    }
    true
}

/// One node of the Gray-code recursion: a set of parities still to be
/// produced, the columns not yet split on, the qubit they are being
/// accumulated onto, and a CNOT owed before the node can proceed.
struct Pt {
    /// The columns not yet split on, as a bit set -- the same shape a
    /// [`Parity`] uses, and for the same reason: a block spans at most
    /// [`MAX_CHUNK_QUBITS`] columns, so the set fits in a register and a node
    /// costs no allocation to narrow.
    candidates: Parity,
    target: Option<usize>,
    pending: Option<usize>,
    vectors: ParitySet,
}

/// A set of parities with the angle each carries, as [`gray_synth`] splits it.
type ParitySet = Vec<(Parity, f64)>;

/// Retired [`Pt::vectors`] buffers. The recursion splits one parity set into
/// two and drops the parent, so without a pool every node would allocate;
/// with one, the whole recursion runs on the buffers its first few nodes
/// established.
type Pool = Vec<ParitySet>;

/// Hand `buffer` back for a later node to fill.
fn recycle(pool: &mut Pool, mut buffer: ParitySet) {
    if buffer.capacity() > 0 {
        buffer.clear();
        pool.push(buffer);
    }
}

/// Gray-code phase synthesis: emit CNOTs that walk the qubits through every
/// parity that carries a rotation, applying each rotation as its parity comes
/// up.
///
/// The recursion splits the parity set on the column that separates it most
/// evenly — parities agreeing on that column share the CNOT that establishes
/// it, so the split that shares the most work is the one to take. A node with
/// no columns left holds parities that agree everywhere, i.e. one parity, and
/// its rotation is applied to the qubit the node accumulated onto.
///
/// `state` is updated to the linear map the emitted CNOTs leave behind, which
/// is generally *not* the one the block needs; [`linear_synth`] corrects it.
///
/// Returns `false` when `budget` ran out, meaning synthesis has already lost
/// to the block it would replace and the caller should abandon it.
fn gray_synth(
    n: usize,
    phases: ParitySet,
    state: &mut [Parity],
    budget: &mut Budget,
    qubits: &[Qubit],
    scratch: &mut Scratch,
) -> bool {
    let Scratch { pool, stack, .. } = scratch;
    // A synthesis abandoned mid-way leaves nodes behind; reclaim their
    // buffers rather than dropping them.
    for node in stack.drain(..) {
        recycle(pool, node.vectors);
    }
    if phases.is_empty() {
        recycle(pool, phases);
        return true;
    }
    let all_columns = match n >= Parity::BITS as usize {
        true => Parity::MAX,
        false => ((1 as Parity) << n) - 1,
    };
    stack.push(Pt {
        candidates: all_columns,
        target: None,
        pending: None,
        vectors: phases,
    });

    while let Some(mut node) = stack.pop() {
        if node.vectors.is_empty() {
            recycle(pool, node.vectors);
            continue;
        }
        if let (Some(t), Some(p)) = (node.target, node.pending) {
            if !budget.push(Gate::cnot {
                control: qubits[p],
                target: qubits[t],
            }) {
                return false;
            }
            state[t] ^= state[p];
            // Qubit `t` now holds the old `t` XOR `p`. Any parity still to be
            // produced that used the old `t` can keep using the new one as
            // long as the `p` it dragged in is cancelled back out.
            for other in stack.iter_mut() {
                for (bv, _) in other.vectors.iter_mut() {
                    if *bv >> t & 1 == 1 {
                        *bv ^= 1u128 << p;
                    }
                }
            }
            node.pending = None;
            stack.push(node);
            continue;
        }
        if node.candidates == 0 {
            // Out of columns: the node's parities are all equal, so there is
            // at most one, and it is sitting on `target` right now.
            if let (Some(t), 1) = (node.target, node.vectors.len())
                && !emit_rotation(budget, qubits[t], node.vectors[0].1)
            {
                return false;
            }
            recycle(pool, node.vectors);
            continue;
        }
        let target = node.target;
        let col = best_column(node.candidates, &node.vectors);
        let rest = node.candidates & !((1 as Parity) << col);
        let (zeros, ones) = split_on(pool, &node.vectors, col);
        recycle(pool, node.vectors);
        let zero_node = Pt {
            candidates: rest,
            target,
            pending: None,
            vectors: zeros,
        };
        let one_node = match target {
            // Already accumulating onto `t`: reaching the parities that carry
            // this column costs a CNOT from it.
            Some(_) => Pt {
                candidates: rest,
                target,
                pending: Some(col),
                vectors: ones,
            },
            // Nothing accumulated yet, so this column's own qubit already
            // holds it and becomes the target for free.
            None => Pt {
                candidates: rest,
                target: Some(col),
                pending: None,
                vectors: ones,
            },
        };
        // Popped LIFO, so the zero side runs first.
        stack.push(one_node);
        stack.push(zero_node);
    }
    true
}

/// Pick the column in `candidates` whose 0/1 split of `vectors` leaves the
/// largest side: parities agreeing on that column share the CNOT that
/// establishes it, so the most lopsided split shares the most work.
///
/// `candidates` is non-empty, so there is always a column to return. Ties go
/// to the highest column, which is what the ascending scan leaves behind.
fn best_column(candidates: Parity, vectors: &[(Parity, f64)]) -> usize {
    let mut best = 0;
    let mut best_score = usize::MIN;
    let mut rest = candidates;
    while rest != 0 {
        let c = rest.trailing_zeros() as usize;
        rest &= rest - 1;
        let ones = vectors.iter().filter(|(bv, _)| bv >> c & 1 == 1).count();
        let score = ones.max(vectors.len() - ones);
        if score >= best_score {
            best_score = score;
            best = c;
        }
    }
    best
}

/// Split `vectors` on `col` into the parities that clear it and the ones that
/// carry it, in that order, reusing buffers from `pool`.
fn split_on(pool: &mut Pool, vectors: &[(Parity, f64)], col: usize) -> (ParitySet, ParitySet) {
    let mut zeros = pool.pop().unwrap_or_default();
    let mut ones = pool.pop().unwrap_or_default();
    for &entry in vectors {
        match entry.0 >> col & 1 == 1 {
            true => ones.push(entry),
            false => zeros.push(entry),
        }
    }
    (zeros, ones)
}

/// Emit CNOTs taking the qubits from parities `from` to parities `to`.
///
/// Both are invertible, so `m = to * from^-1` is the map still to be applied.
/// Reducing `m` to the identity by row operations expresses it as a product of
/// elementary row additions, each of which *is* a CNOT; replaying them in
/// reverse applies `m` itself.
/// Returns `false` when `budget` ran out, on the same terms as
/// [`gray_synth`].
fn linear_synth(
    from: &[Parity],
    to: &[Parity],
    qubits: &[Qubit],
    budget: &mut Budget,
    scratch: &mut Scratch,
) -> bool {
    let n = from.len();
    if from == to {
        return true;
    }
    let Scratch {
        elim,
        inverse,
        m,
        ops,
        ..
    } = scratch;
    if !invert_into(from, n, elim, inverse) {
        // Unreachable: the qubit parities stay a basis for as long as a block
        // lives. Giving up keeps the caller's comparison honest -- synthesis
        // simply loses to the original.
        debug_assert!(false, "block linear map was singular");
        return false;
    }
    m.clear();
    m.extend(to.iter().map(|&row| row_times_matrix(row, inverse)));

    ops.clear();
    for col in 0..n {
        // Every op below becomes exactly one CNOT, so the budget can be
        // checked against the running op count without emitting anything.
        if !budget.admits_2q(ops.len()) {
            return false;
        }
        if m[col] >> col & 1 == 0 {
            // Columns left of `col` are already cleared to the identity, so a
            // donor row must come from below to keep them that way; one always
            // exists because `m` is invertible.
            let Some(r) = (col + 1..n).find(|&r| m[r] >> col & 1 == 1) else {
                debug_assert!(false, "block linear map was singular");
                return false;
            };
            m[col] ^= m[r];
            ops.push((r, col));
        }
        for r in 0..n {
            if r != col && m[r] >> col & 1 == 1 {
                m[r] ^= m[col];
                ops.push((col, r));
            }
        }
    }
    for &(control, target) in ops.iter().rev() {
        if !budget.push(Gate::cnot {
            control: qubits[control],
            target: qubits[target],
        }) {
            return false;
        }
    }
    true
}

/// XOR together the rows of `m` selected by the set bits of `row`.
fn row_times_matrix(row: Parity, m: &[Parity]) -> Parity {
    let mut out = 0;
    let mut bits = row;
    while bits != 0 {
        let i = bits.trailing_zeros() as usize;
        bits &= bits - 1;
        out ^= m[i];
    }
    out
}

/// Invert an `n`x`n` bit matrix by Gauss-Jordan into `inv`, returning `false`
/// if it is singular. `a` is scratch space for the reduction.
fn invert_into(rows: &[Parity], n: usize, a: &mut Vec<Parity>, inv: &mut Vec<Parity>) -> bool {
    a.clear();
    a.extend_from_slice(&rows[..n]);
    inv.clear();
    inv.extend((0..n).map(|i| (1 as Parity) << i));
    for col in 0..n {
        let Some(pivot) = (col..n).find(|&r| a[r] >> col & 1 == 1) else {
            return false;
        };
        a.swap(col, pivot);
        inv.swap(col, pivot);
        for r in 0..n {
            if r != col && a[r] >> col & 1 == 1 {
                a[r] ^= a[col];
                inv[r] ^= inv[col];
            }
        }
    }
    true
}

/// [`invert_into`] with its scratch space owned, for the tests that check the
/// inversion itself rather than the pass around it.
#[cfg(test)]
fn invert(rows: &[Parity], n: usize) -> Option<Vec<Parity>> {
    let (mut a, mut inv) = (Vec::new(), Vec::new());
    invert_into(rows, n, &mut a, &mut inv).then_some(inv)
}

fn angle_is_zero(angle: f64) -> bool {
    let n = angle.rem_euclid(2.0 * PI);
    n < ANGLE_TOL || (2.0 * PI - n) < ANGLE_TOL
}

/// Append `angle` on `qubit`, as Clifford+T gates when it is a multiple of
/// pi/4 and an `rz` otherwise.
/// Returns `false` when `budget` ran out, on the same terms as
/// [`gray_synth`].
#[must_use]
fn emit_rotation(budget: &mut Budget, qubit: Qubit, angle: f64) -> bool {
    let n = angle.rem_euclid(2.0 * PI);
    if angle_is_zero(n) {
        return true;
    }
    let quarter = PI / 4.0;
    let k = (n / quarter).round();
    if (n - k * quarter).abs() >= ANGLE_TOL {
        return budget.push(Gate::rz(n, qubit));
    }
    match k as u32 % 8 {
        0 => true,
        1 => budget.push(Gate::t(qubit)),
        2 => budget.push(Gate::s(qubit)),
        3 => budget.push(Gate::s(qubit)) && budget.push(Gate::t(qubit)),
        4 => budget.push(Gate::z(qubit)),
        5 => budget.push(Gate::z(qubit)) && budget.push(Gate::t(qubit)),
        6 => budget.push(Gate::sdg(qubit)),
        7 => budget.push(Gate::tdg(qubit)),
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unitary::circuits_equiv;

    const TOL: f64 = 1e-9;

    pub(super) struct TestRng(pub(super) u64);

    impl TestRng {
        pub(super) fn next(&mut self, upper: usize) -> usize {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0 as usize % upper
        }

        /// [`TestRng::next`] as a qubit index.
        pub(super) fn qubit(&mut self, upper: usize) -> Qubit {
            self.next(upper) as Qubit
        }
    }

    fn two_qubit(c: &Circuit) -> usize {
        count_2q(&c.gates)
    }

    /// A budget that never runs out, for the tests that exercise synthesis
    /// itself rather than the give-up path.
    fn unbounded() -> Budget {
        Budget::new(usize::MAX, usize::MAX)
    }

    pub(super) fn count_t(gates: &[Gate]) -> usize {
        gates
            .iter()
            .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_)))
            .count()
    }

    /// A random circuit mixing everything the pass interprets with the gates
    /// that end a block, on up to `max_qubits` wires and `max_len` gates.
    pub(super) fn random_mixed_circuit(
        rng: &mut TestRng,
        max_qubits: usize,
        max_len: usize,
    ) -> Circuit {
        let n = 1 + rng.next(max_qubits);
        let mut c = Circuit::new(n);
        for _ in 0..rng.next(max_len) {
            let q = rng.qubit(n);
            match rng.next(9) {
                0 => c.apply(Gate::t(q)),
                1 => c.apply(Gate::tdg(q)),
                2 => c.apply(Gate::s(q)),
                3 => c.apply(Gate::sdg(q)),
                4 => c.apply(Gate::z(q)),
                5 => c.apply(Gate::x(q)),
                6 => c.apply(Gate::h(q)),
                7 if n > 1 => {
                    let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                    c.apply(Gate::cz {
                        control: q,
                        target: t,
                    });
                }
                _ if n > 1 => {
                    let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                    c.apply(Gate::cnot {
                        control: q,
                        target: t,
                    });
                }
                _ => c.apply(Gate::z(q)),
            }
        }
        c
    }

    /// Apply a CNOT list to a parity state, the way the hardware would.
    fn replay(gates: &[Gate], qubits: &[Qubit], n: usize) -> Vec<Parity> {
        let mut state: Vec<Parity> = (0..n).map(|i| 1u128 << i).collect();
        for gate in gates {
            if let Gate::cnot { control, target } = gate {
                let c = qubits.iter().position(|q| q == control).unwrap();
                let t = qubits.iter().position(|q| q == target).unwrap();
                state[t] ^= state[c];
            }
        }
        state
    }

    // --- bit-matrix primitives ---

    #[test]
    fn invert_identity_is_identity() {
        let id: Vec<Parity> = (0..5).map(|i| 1u128 << i).collect();
        assert_eq!(invert(&id, 5).unwrap(), id);
    }

    #[test]
    fn invert_composed_with_original_is_identity() {
        let mut rng = TestRng(0x9e37_79b9_7f4a_7c15);
        for _ in 0..200 {
            let n = 1 + rng.next(8);
            // A random invertible matrix, built as a product of elementary
            // row additions so it is invertible by construction.
            let mut m: Vec<Parity> = (0..n).map(|i| 1u128 << i).collect();
            for _ in 0..(3 * n) {
                let a = rng.next(n);
                let b = rng.next(n);
                if a != b {
                    m[a] ^= m[b];
                }
            }
            let inv = invert(&m, n).unwrap();
            for (i, &row) in m.iter().enumerate() {
                assert_eq!(row_times_matrix(row, &inv), 1u128 << i, "n={n}");
            }
        }
    }

    #[test]
    fn invert_rejects_singular() {
        // Two equal rows cannot be a basis.
        let m: Vec<Parity> = vec![0b01, 0b01, 0b100];
        assert!(invert(&m, 3).is_none());
    }

    // --- linear synthesis ---

    #[test]
    fn linear_synth_realizes_random_maps() {
        let mut rng = TestRng(0x243f_6a88_85a3_08d3);
        for _ in 0..300 {
            let n = 1 + rng.next(8);
            let qubits: Vec<Qubit> = (0..n as Qubit).collect();
            let from: Vec<Parity> = (0..n).map(|i| 1u128 << i).collect();
            let mut to = from.clone();
            for _ in 0..(3 * n) {
                let a = rng.next(n);
                let b = rng.next(n);
                if a != b {
                    to[a] ^= to[b];
                }
            }
            let mut b = unbounded();
            assert!(linear_synth(
                &from,
                &to,
                &qubits,
                &mut b,
                &mut Scratch::default()
            ));
            assert_eq!(replay(&b.gates, &qubits, n), to, "n={n}");
        }
    }

    #[test]
    fn linear_synth_emits_nothing_for_identity() {
        let from: Vec<Parity> = (0..4).map(|i| 1u128 << i).collect();
        let mut b = unbounded();
        assert!(linear_synth(
            &from,
            &from.clone(),
            &[0, 1, 2, 3],
            &mut b,
            &mut Scratch::default(),
        ));
        assert!(b.gates.is_empty());
    }

    #[test]
    fn linear_synth_realizes_maps_from_a_nonidentity_start() {
        let mut rng = TestRng(0xb7e1_5162_8aed_2a6b);
        for _ in 0..200 {
            let n = 2 + rng.next(6);
            let qubits: Vec<Qubit> = (0..n as Qubit).collect();
            let mut from: Vec<Parity> = (0..n).map(|i| 1u128 << i).collect();
            for _ in 0..(2 * n) {
                let (a, b) = (rng.next(n), rng.next(n));
                if a != b {
                    from[a] ^= from[b];
                }
            }
            let mut to = from.clone();
            for _ in 0..(2 * n) {
                let (a, b) = (rng.next(n), rng.next(n));
                if a != b {
                    to[a] ^= to[b];
                }
            }
            let mut b = unbounded();
            assert!(linear_synth(
                &from,
                &to,
                &qubits,
                &mut b,
                &mut Scratch::default()
            ));
            // Replaying starts from the identity, so compose onto `from`.
            let applied = replay(&b.gates, &qubits, n);
            let composed: Vec<Parity> = applied
                .iter()
                .map(|&row| row_times_matrix(row, &from))
                .collect();
            assert_eq!(composed, to, "n={n}");
        }
    }

    // --- Gray-code synthesis ---

    #[test]
    fn gray_synth_places_every_rotation_on_its_own_parity() {
        // The invariant the whole pass rests on: when a rotation is emitted,
        // the qubit it lands on is holding exactly that rotation's parity.
        let mut rng = TestRng(0xf3bc_c908_b2fb_1366);
        for _ in 0..200 {
            let n = 1 + rng.next(6);
            let count = 1 + rng.next(8);
            let mut wanted: Vec<Parity> = Vec::new();
            for _ in 0..count {
                let p = (rng.next(1 << n) as u128) & ((1u128 << n) - 1);
                if p != 0 && !wanted.contains(&p) {
                    wanted.push(p);
                }
            }
            if wanted.is_empty() {
                continue;
            }
            let phases: Vec<(Parity, f64)> = wanted.iter().map(|&p| (p, PI / 4.0)).collect();
            let qubits: Vec<Qubit> = (0..n as Qubit).collect();
            let mut state: Vec<Parity> = (0..n).map(|i| 1u128 << i).collect();
            let mut b = unbounded();
            assert!(gray_synth(
                n,
                phases,
                &mut state,
                &mut b,
                &qubits,
                &mut Scratch::default(),
            ));
            let gates = b.gates;

            // Walk the emitted gates, checking each rotation sits on its parity.
            let mut live: Vec<Parity> = (0..n).map(|i| 1u128 << i).collect();
            let mut seen = Vec::new();
            for gate in &gates {
                match gate {
                    Gate::cnot { control, target } => {
                        live[*target as usize] ^= live[*control as usize]
                    }
                    Gate::t(q) => seen.push(live[*q as usize]),
                    other => panic!("unexpected gate {other:?}"),
                }
            }
            assert_eq!(live, state, "tracked state disagrees with gray_synth's");
            seen.sort_unstable();
            wanted.sort_unstable();
            assert_eq!(seen, wanted, "n={n}");
        }
    }

    #[test]
    fn gray_synth_on_no_phases_emits_nothing() {
        let mut state: Vec<Parity> = (0..4).map(|i| 1u128 << i).collect();
        let mut b = unbounded();
        assert!(gray_synth(
            4,
            Vec::new(),
            &mut state,
            &mut b,
            &[0, 1, 2, 3],
            &mut Scratch::default(),
        ));
        assert!(b.gates.is_empty());
        assert_eq!(state, (0..4).map(|i| 1u128 << i).collect::<Vec<_>>());
    }

    // --- end-to-end equivalence ---

    #[test]
    fn preserves_a_cnot_ladder() {
        let mut c = Circuit::new(4);
        for i in 0..3 {
            c.apply(Gate::cnot {
                control: i,
                target: i + 1,
            });
        }
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn cancels_a_repeated_cnot_pair() {
        // Two identical CNOTs are the identity; resynthesis should see the
        // linear map is trivial and emit nothing at all.
        let mut c = Circuit::new(2);
        for _ in 0..2 {
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
        }
        let out = cnot_min(&c);
        assert_eq!(two_qubit(&out), 0);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn shrinks_a_redundant_linear_map() {
        // Six CNOTs that compose to a map reachable in fewer.
        let mut c = Circuit::new(3);
        let pairs = [(0, 1), (1, 2), (0, 1), (1, 2), (0, 1), (1, 2)];
        for (control, target) in pairs {
            c.apply(Gate::cnot { control, target });
        }
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
        assert!(two_qubit(&out) < 6, "got {}", two_qubit(&out));
    }

    #[test]
    fn merges_rotations_on_one_parity() {
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
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
        assert!(out.gates.len() <= c.gates.len());
    }

    #[test]
    fn hadamard_splits_blocks_but_preserves_semantics() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));
        c.apply(Gate::h(1));
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::t(2));
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
        assert!(out.gates.iter().any(|g| matches!(g, Gate::h(1))));
    }

    #[test]
    fn x_gates_are_preserved_through_resynthesis() {
        let mut c = Circuit::new(3);
        c.apply(Gate::x(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(1));
        c.apply(Gate::x(1));
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::tdg(2));
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn cz_is_absorbed_and_preserved() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::cnot {
            control: 1,
            target: 2,
        });
        c.apply(Gate::cz {
            control: 0,
            target: 2,
        });
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn ccz_and_ccx_act_as_block_boundaries() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::ccz {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::t(0));
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
        assert_eq!(
            out.gates
                .iter()
                .filter(|g| matches!(g, Gate::ccz { .. } | Gate::ccx { .. }))
                .count(),
            2
        );
    }

    #[test]
    fn rz_rotations_survive() {
        let mut c = Circuit::new(2);
        c.apply(Gate::rz(0.37, 0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::rz(-1.2, 1));
        let out = cnot_min(&c);
        assert!(circuits_equiv(&c, &out, TOL));
    }

    #[test]
    fn measurement_and_reset_are_preserved_and_bound_blocks() {
        let mut c = Circuit::with_cbits(2, 1);
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::reset(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        let out = cnot_min(&c);
        assert!(out.has_measurement());
        assert_eq!(out.num_cbits, 1);
        assert!(matches!(out.gates[1], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(matches!(out.gates[2], Gate::reset(1)));
    }

    // --- randomized equivalence, the real safety net ---

    /// Random CNOT-dihedral circuits: every gate is one this pass interprets,
    /// so the whole circuit is a single block and gets fully resynthesized.
    #[test]
    fn random_cnot_dihedral_circuits_are_preserved() {
        let mut rng = TestRng(0x0123_4567_89ab_cdef);
        for case in 0..400 {
            let n = 1 + rng.next(4);
            let len = rng.next(24);
            let mut c = Circuit::new(n);
            for _ in 0..len {
                let q = rng.qubit(n);
                match rng.next(8) {
                    0 => c.apply(Gate::t(q)),
                    1 => c.apply(Gate::tdg(q)),
                    2 => c.apply(Gate::s(q)),
                    3 => c.apply(Gate::z(q)),
                    4 => c.apply(Gate::x(q)),
                    5 if n > 1 => {
                        let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cz {
                            control: q,
                            target: t,
                        });
                    }
                    _ if n > 1 => {
                        let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cnot {
                            control: q,
                            target: t,
                        });
                    }
                    _ => c.apply(Gate::sdg(q)),
                }
            }
            let out = cnot_min(&c);
            assert!(circuits_equiv(&c, &out, TOL), "case {case}");
        }
    }

    /// The same, with block-splitting gates mixed in.
    #[test]
    fn random_mixed_circuits_are_preserved() {
        let mut rng = TestRng(0xdead_beef_cafe_1234);
        for case in 0..400 {
            let n = 1 + rng.next(4);
            let len = rng.next(24);
            let mut c = Circuit::new(n);
            for _ in 0..len {
                let q = rng.qubit(n);
                match rng.next(10) {
                    0 => c.apply(Gate::t(q)),
                    1 => c.apply(Gate::tdg(q)),
                    2 => c.apply(Gate::h(q)),
                    3 => c.apply(Gate::s(q)),
                    4 => c.apply(Gate::x(q)),
                    5 => c.apply(Gate::sdg(q)),
                    6 => c.apply(Gate::rz(0.1 + 0.3 * (rng.next(7) as f64), q)),
                    7 if n > 2 => c.apply(Gate::ccz {
                        control1: 0,
                        control2: 1,
                        target: 2,
                    }),
                    _ if n > 1 => {
                        let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cnot {
                            control: q,
                            target: t,
                        });
                    }
                    _ => c.apply(Gate::z(q)),
                }
            }
            let out = cnot_min(&c);
            assert!(circuits_equiv(&c, &out, TOL), "case {case}");
        }
    }

    /// Blocks forced to split mid-stream by a tiny qubit cap must still
    /// compose back to the same operator.
    #[test]
    fn tight_qubit_cap_preserves_semantics() {
        let mut rng = TestRng(0x5555_aaaa_3333_cccc);
        for cap in 1..=4 {
            for case in 0..120 {
                let n = 4;
                let mut c = Circuit::new(n);
                for _ in 0..rng.next(20) {
                    let q = rng.qubit(n);
                    match rng.next(5) {
                        0 => c.apply(Gate::t(q)),
                        1 => c.apply(Gate::x(q)),
                        2 => c.apply(Gate::h(q)),
                        _ => {
                            let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                            c.apply(Gate::cnot {
                                control: q,
                                target: t,
                            });
                        }
                    }
                }
                let out = cnot_min_with(&c, cap, MAX_CHUNK_TERMS);
                assert!(circuits_equiv(&c, &out, TOL), "cap {cap} case {case}");
            }
        }
    }

    /// Likewise for a term cap that forces a flush mid-rotation-sequence.
    #[test]
    fn tight_term_cap_preserves_semantics() {
        let mut rng = TestRng(0x1111_2222_3333_4444);
        for cap in 1..=3 {
            for case in 0..120 {
                let n = 3;
                let mut c = Circuit::new(n);
                for _ in 0..rng.next(20) {
                    let q = rng.qubit(n);
                    match rng.next(4) {
                        0 => c.apply(Gate::t(q)),
                        1 => c.apply(Gate::tdg(q)),
                        2 => c.apply(Gate::s(q)),
                        _ => {
                            let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                            c.apply(Gate::cnot {
                                control: q,
                                target: t,
                            });
                        }
                    }
                }
                let out = cnot_min_with(&c, MAX_CHUNK_QUBITS, cap);
                assert!(circuits_equiv(&c, &out, TOL), "cap {cap} case {case}");
            }
        }
    }

    // --- monotonicity and determinism ---

    #[test]
    fn never_increases_gate_or_two_qubit_count() {
        let mut rng = TestRng(0x7777_8888_9999_aaaa);
        for case in 0..400 {
            let n = 1 + rng.next(5);
            let mut c = Circuit::new(n);
            for _ in 0..rng.next(30) {
                let q = rng.qubit(n);
                match rng.next(7) {
                    0 => c.apply(Gate::t(q)),
                    1 => c.apply(Gate::h(q)),
                    2 => c.apply(Gate::x(q)),
                    3 => c.apply(Gate::s(q)),
                    _ if n > 1 => {
                        let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cnot {
                            control: q,
                            target: t,
                        });
                    }
                    _ => c.apply(Gate::tdg(q)),
                }
            }
            let out = cnot_min(&c);
            assert!(
                out.gates.len() <= c.gates.len(),
                "case {case}: {} -> {}",
                c.gates.len(),
                out.gates.len()
            );
            assert!(
                two_qubit(&out) <= two_qubit(&c),
                "case {case}: {} -> {} two-qubit",
                two_qubit(&c),
                two_qubit(&out)
            );
        }
    }

    /// The budget must be a pure *early exit*: abandoning a synthesis that
    /// can no longer win has to leave exactly the circuit that running it to
    /// completion and rejecting it would have left. If this ever diverges,
    /// the bound is changing results rather than just saving work.
    /// Every multiple of pi/4 must come out as the Clifford+T gates that
    /// realize it -- never as an `rz`, and never as more gates than the
    /// angle needs. Rotations are the pass's whole output apart from CNOTs,
    /// so a wrong entry here would be a silent miscompile, and an `rz` where
    /// a `t` belongs would leak an Rz into a Clifford+T circuit.
    #[test]
    fn quarter_turn_angles_emit_clifford_t_gates() {
        // Fewest gates realizing k*pi/4: 0 is identity; +-pi/4, +-pi/2 and pi
        // are single gates; 3pi/4 and 5pi/4 are not any single gate in the
        // basis, so they cost two.
        let expected_len = [0, 1, 1, 2, 1, 2, 1, 1];
        for k in 0..8u32 {
            let angle = f64::from(k) * PI / 4.0;
            let mut budget = unbounded();
            assert!(emit_rotation(&mut budget, 0, angle), "k={k}");
            assert_eq!(budget.gates.len(), expected_len[k as usize], "k={k}");
            assert!(
                !budget.gates.iter().any(|g| matches!(g, Gate::rz(..))),
                "k={k} fell back to rz"
            );

            let mut got = Circuit::new(1);
            for gate in &budget.gates {
                got.apply(gate.clone());
            }
            let mut want = Circuit::new(1);
            if k != 0 {
                want.apply(Gate::rz(angle, 0));
            }
            assert!(circuits_equiv(&want, &got, TOL), "k={k} is the wrong gate");
        }
    }

    /// The same angles reached the way the pass actually reaches them --
    /// summed into a block's phase map -- rather than handed to
    /// `emit_rotation` directly, so accumulated floating-point error is
    /// covered too.
    #[test]
    fn accumulated_quarter_turns_stay_clifford_t() {
        for k in 1..8usize {
            // k T gates on one wire, wrapped in CNOTs so the block is worth
            // re-synthesizing and the rotation is actually re-emitted.
            let mut c = Circuit::new(2);
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
            for _ in 0..k {
                c.apply(Gate::t(1));
            }
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
            let out = cnot_min(&c);
            assert!(circuits_equiv(&c, &out, TOL), "k={k}");
            assert!(
                !out.gates.iter().any(|g| matches!(g, Gate::rz(..))),
                "k={k}: Clifford+T input produced an rz"
            );
        }
    }

    /// An angle that is genuinely not a multiple of pi/4 has to stay an `rz`.
    #[test]
    fn non_quarter_turn_angles_stay_rz() {
        let mut budget = unbounded();
        assert!(emit_rotation(&mut budget, 0, 0.37));
        assert!(matches!(budget.gates.as_slice(), [Gate::rz(..)]));
    }

    #[test]
    fn budget_is_a_pure_early_exit() {
        let mut rng = TestRng(0xc0ff_ee00_1234_5678);
        for case in 0..600 {
            let n = 1 + rng.next(6);
            let mut c = Circuit::new(n);
            for _ in 0..rng.next(40) {
                let q = rng.qubit(n);
                match rng.next(9) {
                    0 => c.apply(Gate::t(q)),
                    1 => c.apply(Gate::tdg(q)),
                    2 => c.apply(Gate::s(q)),
                    3 => c.apply(Gate::h(q)),
                    4 => c.apply(Gate::x(q)),
                    5 => c.apply(Gate::rz(0.1 + 0.3 * (rng.next(7) as f64), q)),
                    6 if n > 1 => {
                        let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cz {
                            control: q,
                            target: t,
                        });
                    }
                    _ if n > 1 => {
                        let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cnot {
                            control: q,
                            target: t,
                        });
                    }
                    _ => c.apply(Gate::z(q)),
                }
            }
            assert_eq!(
                cnot_min(&c).gates,
                cnot_min_unbounded(&c).gates,
                "case {case}: budget changed the result"
            );
        }
    }

    /// A wide, shallow block is the shape the budget exists for: Gaussian
    /// elimination would emit far more CNOTs than the block was written with,
    /// so synthesis must be abandoned and the original kept.
    #[test]
    fn a_wide_shallow_block_is_abandoned_not_expanded() {
        let n = 60;
        let mut rng = TestRng(0xabcd_0987_6543_210f);
        let mut c = Circuit::new(n);
        for _ in 0..(2 * n) {
            let a = rng.qubit(n);
            let b = (a + 1 + rng.qubit(n - 1)) % n as Qubit;
            c.apply(Gate::cnot {
                control: a,
                target: b,
            });
            c.apply(Gate::t(rng.qubit(n)));
        }
        let out = cnot_min(&c);
        assert_eq!(out.gates, cnot_min_unbounded(&c).gates);
        assert!(
            two_qubit(&out) <= two_qubit(&c),
            "{} -> {}",
            two_qubit(&c),
            two_qubit(&out)
        );
    }

    /// T-count must not rise either. [`Chunk::flush`] only compares
    /// `(two-qubit, total)`, so nothing in the accept test looks at T
    /// directly -- it holds because the phase map merges every rotation on a
    /// parity into one, and [`emit_rotation`] spends at most one T on it.
    /// Worth pinning precisely because it is a consequence rather than a
    /// check: T is the expensive resource these circuits are measured in.
    #[test]
    fn never_increases_t_count() {
        let mut rng = TestRng(0x1234_5678_9abc_def0);
        for case in 0..2000 {
            let c = random_mixed_circuit(&mut rng, 5, 30);
            let out = cnot_min(&c);
            assert!(
                count_t(&out.gates) <= count_t(&c.gates),
                "case {case}: T {} -> {}",
                count_t(&c.gates),
                count_t(&out.gates)
            );
        }
    }

    /// The pass is not idempotent -- a second application can find more,
    /// since the first one reshapes the blocks the second sees -- but it must
    /// *settle*, and never grow on the way. That is what lets it sit in a
    /// fixpoint loop: an oscillating pass would spin forever, and a growing
    /// one would undo the loop's progress.
    ///
    /// Measured over 200k random circuits, 0.16% were not idempotent after
    /// one extra application and every one settled within three.
    #[test]
    fn repeated_application_reaches_a_fixed_point() {
        let mut rng = TestRng(0xfeed_beef_1234_5678);
        for case in 0..2000 {
            let c = random_mixed_circuit(&mut rng, 6, 40);
            let mut cur = cnot_min(&c);
            let mut settled = false;
            for _ in 0..8 {
                let next = cnot_min(&cur);
                assert!(
                    next.gates.len() <= cur.gates.len(),
                    "case {case}: grew on re-application"
                );
                assert!(count_2q(&next.gates) <= count_2q(&cur.gates), "case {case}");
                if next.gates == cur.gates {
                    settled = true;
                    break;
                }
                cur = next;
            }
            assert!(settled, "case {case}: no fixed point within 8 applications");
            assert!(circuits_equiv(&c, &cur, TOL), "case {case}");
        }
    }

    /// A circuit wider than [`MAX_CHUNK_QUBITS`] exercises the cap's
    /// flush-and-retry path and the escape hatch for a gate that no empty
    /// block could ever admit. Too wide to compare unitaries, so this checks
    /// what can be checked: it terminates, stays monotone, and touches no
    /// qubit the input did not.
    #[test]
    fn circuits_wider_than_the_qubit_cap_are_handled() {
        let n = MAX_CHUNK_QUBITS + 72;
        let mut rng = TestRng(0x0fed_cba9_8765_4321);
        for case in 0..20 {
            let mut c = Circuit::new(n);
            for _ in 0..4000 {
                let a = rng.qubit(n);
                match rng.next(6) {
                    0 => c.apply(Gate::t(a)),
                    1 => c.apply(Gate::h(a)),
                    2 => c.apply(Gate::x(a)),
                    _ => {
                        let b = (a + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cnot {
                            control: a,
                            target: b,
                        });
                    }
                }
            }
            let out = cnot_min(&c);
            assert!(out.gates.len() <= c.gates.len(), "case {case}");
            assert!(count_2q(&out.gates) <= count_2q(&c.gates), "case {case}");
            assert_eq!(out.num_qubits, n);
        }
    }

    #[test]
    fn is_deterministic() {
        let mut rng = TestRng(0xfeed_face_0bad_f00d);
        for _ in 0..50 {
            let n = 1 + rng.next(4);
            let mut c = Circuit::new(n);
            for _ in 0..rng.next(20) {
                let q = rng.qubit(n);
                match rng.next(4) {
                    0 => c.apply(Gate::t(q)),
                    1 => c.apply(Gate::h(q)),
                    _ if n > 1 => {
                        let t = (q + 1 + rng.qubit(n - 1)) % n as Qubit;
                        c.apply(Gate::cnot {
                            control: q,
                            target: t,
                        });
                    }
                    _ => c.apply(Gate::s(q)),
                }
            }
            let a = cnot_min(&c);
            let b = cnot_min(&c);
            assert_eq!(a.gates, b.gates);
        }
    }

    #[test]
    fn empty_and_single_gate_circuits_round_trip() {
        let empty = Circuit::new(3);
        assert_eq!(cnot_min(&empty).gates.len(), 0);

        let mut one = Circuit::new(2);
        one.apply(Gate::h(1));
        assert_eq!(cnot_min(&one).gates, one.gates);

        let mut t_only = Circuit::new(1);
        t_only.apply(Gate::t(0));
        let out = cnot_min(&t_only);
        assert!(circuits_equiv(&t_only, &out, TOL));
    }

    #[test]
    fn rotations_summing_to_identity_disappear() {
        let mut c = Circuit::new(1);
        for _ in 0..8 {
            c.apply(Gate::t(0));
        }
        let out = cnot_min(&c);
        assert_eq!(out.gates.len(), 0);
    }

    #[test]
    fn degenerate_two_qubit_gates_pass_through_untouched() {
        // `cx q,q` is ill-formed but appears in the Feynman corpus
        // (cycle_17_3). Interpreting it would zero a parity and leave the
        // block's linear map singular, so it must be copied through verbatim
        // and act as a block boundary.
        for degenerate in [
            Gate::cnot {
                control: 1,
                target: 1,
            },
            Gate::cz {
                control: 1,
                target: 1,
            },
        ] {
            let mut c = Circuit::new(3);
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
            c.apply(Gate::t(1));
            c.apply(degenerate.clone());
            c.apply(Gate::cnot {
                control: 1,
                target: 2,
            });
            c.apply(Gate::t(2));
            let out = cnot_min(&c);
            assert_eq!(
                out.gates.iter().filter(|g| **g == degenerate).count(),
                1,
                "{degenerate:?} was not preserved"
            );
        }
    }

    #[test]
    fn a_circuit_of_only_degenerate_gates_is_unchanged() {
        let mut c = Circuit::new(2);
        for _ in 0..3 {
            c.apply(Gate::cnot {
                control: 0,
                target: 0,
            });
        }
        assert_eq!(cnot_min(&c).gates, c.gates);
    }

    #[test]
    fn qubits_outside_a_block_are_untouched() {
        // A 40-qubit circuit whose activity is confined to three wires must
        // not acquire gates anywhere else.
        let mut c = Circuit::new(40);
        c.apply(Gate::cnot {
            control: 3,
            target: 7,
        });
        c.apply(Gate::t(7));
        c.apply(Gate::cnot {
            control: 7,
            target: 11,
        });
        let out = cnot_min(&c);
        for gate in &out.gates {
            let (n, qs) = crate::circuit::qubit_operands(gate);
            for &q in &qs[..n] {
                assert!(matches!(q, 3 | 7 | 11), "stray gate on qubit {q}");
            }
        }
    }
}

/// Long-running randomized checks. Not part of the default run -- they take
/// minutes, and the fast tests cover the same ground at smaller scale -- but
/// they are the pass's real safety net when it changes.
///
///     cargo test --release --lib cnot_min::long_running -- --ignored --nocapture
#[cfg(test)]
mod long_running {
    use super::tests::{TestRng, count_t, random_mixed_circuit};
    use super::*;
    use crate::unitary::circuits_equiv;

    const TOL: f64 = 1e-9;

    /// Randomized equivalence over far more circuits, and wider ones, than
    /// the fast tests can afford. Unitary comparison is exponential in the
    /// qubit count, which is what bounds the width here.
    #[test]
    #[ignore] // long-running: 200k random circuits with unitary equivalence checks
    fn randomized_equivalence_at_scale() {
        let mut rng = TestRng(0x5eed_1234_abcd_9876);
        let cases = 200_000;
        for case in 0..cases {
            let c = random_mixed_circuit(&mut rng, 7, 60);
            let out = cnot_min(&c);
            assert!(circuits_equiv(&c, &out, TOL), "case {case}: {:?}", c.gates);
            assert!(out.gates.len() <= c.gates.len(), "case {case}");
            assert!(count_2q(&out.gates) <= count_2q(&c.gates), "case {case}");
            assert!(count_t(&out.gates) <= count_t(&c.gates), "case {case}");
            if case % 20_000 == 0 {
                println!("  ..{case}/{cases}");
            }
        }
        println!("{cases} random circuits: equivalent, and monotone in gates, 2q and T");
    }

    /// The same, with the block caps driven to values that force splits
    /// mid-block -- the paths a normally-sized circuit never reaches.
    #[test]
    #[ignore] // long-running: 80k circuits across 20 block-cap combinations
    fn cap_split_equivalence_at_scale() {
        let mut rng = TestRng(0x9876_5432_10fe_dcba);
        let mut checked = 0;
        for qubit_cap in [1, 2, 3, 5, 8] {
            for term_cap in [1, 2, 4, 16] {
                for case in 0..4_000 {
                    let c = random_mixed_circuit(&mut rng, 6, 40);
                    let out = cnot_min_with(&c, qubit_cap, term_cap);
                    assert!(
                        circuits_equiv(&c, &out, TOL),
                        "qubit_cap {qubit_cap} term_cap {term_cap} case {case}: {:?}",
                        c.gates
                    );
                    assert!(out.gates.len() <= c.gates.len());
                    checked += 1;
                }
            }
        }
        println!("{checked} circuits across 20 cap combinations: equivalent and monotone");
    }

    /// Every synthesis the budget abandons must leave exactly what running it
    /// to completion and rejecting it would have left, at scale.
    #[test]
    #[ignore] // long-running: 100k circuits, budgeted against unbudgeted synthesis
    fn budget_early_exit_at_scale() {
        let mut rng = TestRng(0xc0de_f00d_5555_aaaa);
        let cases = 100_000;
        for case in 0..cases {
            let c = random_mixed_circuit(&mut rng, 7, 60);
            assert_eq!(
                cnot_min(&c).gates,
                cnot_min_unbounded(&c).gates,
                "case {case}: budget changed the result on {:?}",
                c.gates
            );
        }
        println!("{cases} random circuits: budget is a pure early exit");
    }
}

#[cfg(test)]
mod cap_bench {
    use super::*;
    use std::time::Instant;

    /// Not an assertion so much as a record: how the term cap trades runtime
    /// against quality on an adversarial block (one huge CNOT-dihedral region
    /// whose rotations sit on thousands of distinct parities).
    #[test]
    #[ignore = "timing sweep, run manually with --ignored --nocapture"]
    fn term_cap_sweep() {
        let n = 100;
        let mut seed = 7u64;
        let mut rng = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        let mut c = Circuit::new(n);
        for i in 0..200_000 {
            if i % 3 == 0 {
                c.apply(Gate::t((rng() % n) as Qubit));
            } else {
                let a = (rng() % n) as Qubit;
                let b = (a + 1 + (rng() % (n - 1)) as Qubit) % n as Qubit;
                c.apply(Gate::cnot {
                    control: a,
                    target: b,
                });
            }
        }
        for cap in [128, 512, 1024, 2048, 4096] {
            let start = Instant::now();
            let out = cnot_min_with(&c, MAX_CHUNK_QUBITS, cap);
            println!(
                "terms<={cap:<5} {:>7.3}s  gates {} -> {}  2q {} -> {}",
                start.elapsed().as_secs_f64(),
                c.gates.len(),
                out.gates.len(),
                count_2q(&c.gates),
                count_2q(&out.gates),
            );
        }
    }
}
