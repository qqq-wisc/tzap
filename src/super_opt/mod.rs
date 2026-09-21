//! Peephole superoptimization over anchored circuit windows.
//!
//! The pass makes one forward scan over the circuit, carves out small
//! *windows* — causally connected groups of gates on a few qubits — computes
//! each window's unitary matrix, and asks a precomputed minimal unitary
//! representative map (MURM) whether the same unitary is reachable with fewer gates. Where
//! it is, the window is replaced. Because replacement is by matrix
//! equivalence rather than by syntactic rule, any identity the MURM can
//! express is discovered without ever being written down.
//!
//! # Windows
//!
//! Every gate *anchors* a window: the connected component of that gate among
//! all gates observed since it, connectivity meaning "shares a qubit". As the
//! scan advances, each new gate extends every window that touches one of its
//! qubits. Two invariants make this cheap and exhaustive:
//!
//! - A window always holds the *full* connected component of its anchor over
//!   the scanned region: gates on a window qubit between the anchor and the
//!   scan head are all members, so per-qubit gate history only needs to be
//!   consulted when a gate bridges a **new** qubit into the window
//!   (`expand_component_closure`). Such a bridge can retroactively pull in
//!   an entire previously unrelated component.
//! - Unrelated gates in between are simply absent, so a window is a
//!   *subsequence* of the circuit, not a contiguous slice. Rewriting it in
//!   place is still sound because every skipped gate commutes past the
//!   window: it shares no qubit with any member between anchor and head.
//!
//! A window dies when it exceeds `max_qubits` or `window_gates`, when Rz,
//! measurement, or reset touches one of its qubits (Rz is outside SuperOpt's
//! exact Clifford+T matrix domain; measurement/reset are non-unitary), or when
//! one of its gates is claimed by a selected rewrite. Each window is analyzed
//! after every extension, so every intermediate size is considered, not just
//! the final one.
//!
//! # Rewrites
//!
//! Windows are looked up as they grow; the first strictly smaller equivalent
//! claims its gates greedily in scan order (`RewriteSet`), later windows
//! overlapping a claimed gate are discarded, and all selected rewrites are
//! applied in one reconstruction pass at the end.
//!
//! # Caching
//!
//! Three layers keep repeated work off the hot path:
//!
//! 1. The MURM is built once per configuration and shared process-wide
//!    (`murm::shared_murm`), and persisted to disk
//!    (under the `murm/` cache subdirectory, one file per distinct config) so
//!    later processes load it instead of rebuilding it.
//! 2. The matrix store (`MatrixStore`) interns each canonical window shape
//!    with its matrix and synthesis outcome — including the negative one —
//!    and is carried across runs of a pass instance, so later fixpoint
//!    sweeps replay recurring shapes as hash hits. It also caches shape
//!    *transitions*, so extending a window by one gate usually resolves from
//!    the pair (previous shape, appended gate) without a key being built.
//! 3. Incremental mode ([`SuperOpt::incremental`]) diffs each input against
//!    the instance's previous input and only anchors windows near the
//!    changes; everything else was analyzed before and selected nothing
//!    (`incremental.rs`).

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use smallvec::{SmallVec, smallvec};

use crate::circuit::{Circuit, Gate, GateKind, GateSet, Qubit};
use crate::pass::Pass;

/// Inline storage for a window's qubit support (bounded by `max_qubits`, at most
/// four in practice) and its gate indices (bounded by `window_gates`), so the
/// per-gate window bookkeeping stays off the heap.
type QubitVec = SmallVec<[Qubit; 4]>;
type IndexVec = SmallVec<[usize; 8]>;

mod config;
mod error;
mod incremental;
mod matrix;
mod matrix_cache;
mod murm;
mod synthesis_arena;

pub use config::MurmConfig;
pub use error::SuperOptError;

/// Ordinary SuperOpt synthesis library.
pub const BASE_GATE_SET: GateSet = GateSet::from_bits_const(
    (1 << GateKind::H as u8)
        | (1 << GateKind::X as u8)
        | (1 << GateKind::Z as u8)
        | (1 << GateKind::S as u8)
        | (1 << GateKind::Sdg as u8)
        | (1 << GateKind::T as u8)
        | (1 << GateKind::Tdg as u8)
        | (1 << GateKind::Cx as u8),
);

/// Native controlled gates SuperOpt may add to its base library.
pub const OPTIONAL_GATE_KINDS: [GateKind; 3] = [GateKind::Cz, GateKind::Ccx, GateKind::Ccz];
pub const OPTIONAL_GATE_SET: GateSet = GateSet::from_bits_const(
    (1 << GateKind::Cz as u8) | (1 << GateKind::Ccx as u8) | (1 << GateKind::Ccz as u8),
);
pub const SUPPORTED_GATE_SET: GateSet =
    GateSet::from_bits_const(BASE_GATE_SET.bits() | OPTIONAL_GATE_SET.bits());

use matrix::UnitaryMatrix;
#[cfg(test)]
use matrix_cache::compact_normalized_key;
use matrix_cache::{CachedMatrix, GateCode, MatrixStore, gate_code};
use murm::{Murm, shared_murm};

/// Whether a MURM matching `config` is already cached on disk —
/// a hint for callers wanting to report whether the next `SuperOpt::new`
/// call will be a fast cache load or a fresh, slow build. Purely
/// informational: `SuperOpt::new` re-validates independently, so this is
/// never load-bearing for correctness.
pub fn murm_is_cached(config: MurmConfig) -> bool {
    murm::disk_cache_exists(config)
}

/// Size in bytes of the on-disk MURM cache for `config`, if one
/// exists — a reporting aid alongside `murm_is_cached`, never load-bearing
/// for correctness.
pub fn murm_cache_size_bytes(config: MurmConfig) -> Option<u64> {
    murm::disk_cache_size_bytes(config)
}

/// One MURM cached on disk: where it lives and how big it is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CacheEntry {
    pub path: std::path::PathBuf,
    pub bytes: u64,
}

/// Point tzap's on-disk caches at `dir` for the rest of the process,
/// overriding `$TZAP_CACHE_DIR` and the XDG lookup. Returns `Err` with the
/// root already in force if one was installed earlier: this must be called
/// once, before any MURM is built or loaded, or a single run could read
/// half its MURMs from one directory and half from another.
pub fn set_cache_dir(dir: std::path::PathBuf) -> Result<(), std::path::PathBuf> {
    murm::set_cache_root(dir)
}

/// The cache root currently in force — `--cache-dir`, `$TZAP_CACHE_DIR`,
/// `$XDG_CACHE_HOME/tzap`, `$HOME/.cache/tzap`, `%LOCALAPPDATA%\tzap`, or
/// `%USERPROFILE%\.cache\tzap`, in that order. `None` when none of them
/// resolve, in which case tzap simply doesn't cache to disk.
pub fn cache_dir() -> Option<std::path::PathBuf> {
    murm::cache_root()
}

/// The MURMs currently cached on disk, oldest name first. An
/// unreadable or absent cache directory reads as empty rather than as an
/// error: a cache nobody can list is, for every purpose tzap has, a cache
/// with nothing in it.
pub fn cache_entries() -> Vec<CacheEntry> {
    let Some(dir) = murm::cache_dir() else {
        return Vec::new();
    };
    let Ok(read_dir) = std::fs::read_dir(&dir) else {
        return Vec::new();
    };
    let mut entries: Vec<CacheEntry> = read_dir
        .flatten()
        .filter(|entry| entry.path().extension().is_some_and(|ext| ext == "bin"))
        .filter_map(|entry| {
            let bytes = entry.metadata().ok().filter(|m| m.is_file())?.len();
            Some(CacheEntry {
                path: entry.path(),
                bytes,
            })
        })
        .collect();
    entries.sort_by(|a, b| a.path.cmp(&b.path));
    entries
}

/// Delete every cached MURM, returning the entries removed.
/// Only the MURM files themselves are touched — never the directory, and
/// never anything else inside it — so a cache root shared with other tools
/// (`$XDG_CACHE_HOME/tzap`, say) can't lose data that isn't tzap's to
/// delete. The first failure stops the walk and is returned, with the
/// entries removed up to that point discarded from the error path.
pub fn clear_cache() -> std::io::Result<Vec<CacheEntry>> {
    let mut removed = Vec::new();
    for entry in cache_entries() {
        std::fs::remove_file(&entry.path)?;
        removed.push(entry);
    }
    Ok(removed)
}

/// Matrix and location information for one completed anchored window.
#[derive(Clone, Debug)]
pub struct SuperOptWindow {
    /// Gate positions in chronological order; unrelated intervening gates are absent.
    pub gate_indices: Vec<usize>,
    /// Sorted physical qubits corresponding to the matrix's local qubit order.
    pub qubits: Vec<Qubit>,
    /// Shared when another canonical gate sequence has the same matrix construction.
    pub matrix: Arc<UnitaryMatrix>,
}

/// One selected semantics-preserving peephole rewrite, in input coordinates.
#[derive(Clone, Debug)]
pub struct SuperOptRewrite {
    /// Indices, into the input circuit, of the gates this rewrite replaces.
    pub gate_indices: Vec<usize>,
    /// Replacement gates on the original circuit's physical qubits.
    pub replacement: Vec<Gate>,
}

/// Results and matrix-cache statistics from [`SuperOpt::run`].
#[derive(Clone, Debug)]
pub struct SuperOptResult {
    /// Input circuit with a non-overlapping set of strictly smaller rewrites applied.
    pub circuit: Circuit,
    /// Completed windows and their matrices; empty when built with
    /// [`SuperOpt::without_subcircuits`].
    pub subcircuits: Vec<SuperOptWindow>,
    /// All selected rewrites, including identity removals with empty replacements.
    pub rewrites: Vec<SuperOptRewrite>,
    /// Reused matrices for completed closed components.
    pub cache_hits: usize,
    /// Completed components with previously unseen canonical gate sequences.
    pub cache_misses: usize,
}

/// Configuration for the connected anchored-window analysis.
#[derive(Clone, Debug)]
pub struct SuperOpt {
    /// Maximum number of distinct qubits in a tracked window.
    pub max_qubits: usize,
    /// Maximum number of connected gates in a reported window.
    pub window_gates: usize,
    collect_subcircuits: bool,
    incremental: bool,
    murm: Option<Arc<Murm>>,
    /// Matrix cache carried across runs of this pass instance (and its
    /// clones), so repeated fixpoint sweeps skip re-deriving recurring window
    /// shapes. Reuse returns exactly what a cold run would recompute; see
    /// `MatrixStore`. `Rc<RefCell<_>>`, not `Arc<Mutex<_>>`: a `SuperOpt`
    /// instance (and any clones sharing this cache) is meant to stay on one
    /// thread — parallel callers construct one fresh instance per worker
    /// rather than sharing or cloning one across threads — so this is never
    /// actually touched from more than one thread at a time; `Pass` doesn't
    /// require `Sync` for exactly this reason.
    store: Rc<RefCell<MatrixStore>>,
    /// The input the previous run saw, diffed against the next input to
    /// bound where new windows can anchor when `incremental` is set.
    prev_input: Rc<RefCell<Option<Circuit>>>,
}

#[derive(Debug)]
struct ActiveWindow {
    gate_indices: IndexVec,
    qubits: QubitVec,
    // Index of the store entry interned for this window's current canonical
    // shape, which `MatrixStore::resolve` uses as the source state of a shape
    // transition when the next gate extends the window without changing its
    // support. `None` before the window has been resolved once, and after any
    // extension whose result the store declined to intern. Replaced the
    // window's own copy of its canonical key: with transitions, the key is
    // only ever needed on the miss path, where it is rebuilt into one shared
    // scratch buffer.
    state: Option<u32>,
    // Tracked incrementally by `expand_component_closure` (the only place an
    // Rz gate can enter `gate_indices`, via a bridged-in qubit's history), so
    // `analyze_window` can reject an Rz window in O(1) instead of rescanning
    // the whole (monotonically growing) gate list on every touch.
    contains_rz: bool,
    // How many of `gate_indices` were already found unclaimed, and the
    // `RewriteSet` version that finding is valid at — see
    // `RewriteSet::window_claims_any`.
    checked_len: u32,
    checked_version: u32,
    // Anchor order: the number of windows created before this one. Slot
    // indices in `slots` get recycled (see `WindowArena`), so a slot index is
    // no longer a proxy for age — but rewrite selection is greedy in anchor
    // order, so the per-gate touched-window list must still be processed
    // oldest-first. This is the key it sorts by.
    //
    // It doubles as the slot's liveness token: `WindowArena::DEAD` marks a
    // retired slot, so a per-qubit registration naming this slot is stale
    // exactly when the seq it recorded no longer matches (see
    // `Registration`).
    seq: u32,
}

/// One window's entry in a qubit's list of live windows: its anchor order in
/// the high 32 bits, its arena slot in the low 32.
///
/// Packing the two serves three purposes at once. Ordering by anchor — what
/// greedy rewrite selection consumes — is a plain `u64` comparison over a
/// contiguous array, with no indirection into the arena and no key closure.
/// The recorded seq validates the registration against the slot's *current*
/// tenant, which is what lets a window be retired without touching the
/// per-qubit lists at all. And a registration is 8 bytes rather than 16.
type Registration = u64;

fn register(seq: u32, slot: u32) -> Registration {
    (u64::from(seq) << 32) | u64::from(slot)
}

fn registered_slot(entry: Registration) -> usize {
    (entry & u64::from(u32::MAX)) as usize
}

fn registered_seq(entry: Registration) -> u32 {
    (entry >> 32) as u32
}

/// Slot storage for the live windows of one scan.
///
/// Nearly every gate anchors a window, but a window dies within
/// `window_gates` gates of its anchor, so only a few dozen are ever live at
/// once (measured on `benchmarks/cobble-t/ols-ridge.qasm` at `-O3`: 35–68 live,
/// against 536k anchored). Indexing slots by a monotonically increasing
/// counter therefore grew `slots` to one slot per gate — tens of megabytes,
/// reallocated and page-faulted on every fixpoint sweep — to hold a working
/// set of a few kilobytes. Retiring a window returns its slot to `free`, so
/// the arena instead stays at high-water-mark size.
///
/// Slot indices are consequently *not* anchor order; `ActiveWindow::seq` is.
///
/// Windows live in their slots and are mutated in place. They used to be moved
/// out and back on every touch — 136 bytes copied twice, tens of millions of
/// times per scan — purely so that `retire` could take `&mut self`; retiring
/// now needs nothing but the slot index.
#[derive(Debug, Default)]
struct WindowArena {
    slots: Vec<ActiveWindow>,
    free: Vec<u32>,
    /// Stamped with `gate_index + 1` when a slot is added to the per-gate
    /// touched list, to dedup it without sorting. 0 means "not yet touched".
    touch_stamp: Vec<usize>,
    next_seq: u32,
}

impl WindowArena {
    /// The `seq` of a retired slot. Never a real anchor order, which is
    /// bounded by the circuit's gate count, so a per-qubit registration for a
    /// window that has since died can never validate against its old slot.
    const DEAD: u32 = u32::MAX;

    fn window(&mut self, slot: usize) -> &mut ActiveWindow {
        &mut self.slots[slot]
    }

    /// Start a window anchored at `anchor` over `qubits`, returning its slot.
    ///
    /// A recycled slot is overwritten field by field rather than assigned a
    /// whole new window, so the previous tenant's `gate_indices` allocation is
    /// reused — `IndexVec` spills to the heap past eight gates, and
    /// `window_gates` is 25 at `-O3`, so replacing the window wholesale freed
    /// and reallocated that buffer once per anchored gate.
    fn anchor(&mut self, anchor: usize, qubits: &[Qubit]) -> u32 {
        let seq = self.next_seq;
        debug_assert_ne!(seq, Self::DEAD, "anchor order exhausted");
        self.next_seq += 1;
        match self.free.pop() {
            Some(slot) => {
                let window = &mut self.slots[slot as usize];
                window.gate_indices.clear();
                window.gate_indices.push(anchor);
                window.qubits.clear();
                window.qubits.extend_from_slice(qubits);
                window.state = None;
                window.contains_rz = false;
                window.checked_len = 0;
                window.checked_version = 0;
                window.seq = seq;
                slot
            }
            None => {
                self.slots.push(ActiveWindow {
                    gate_indices: smallvec![anchor],
                    qubits: QubitVec::from_slice(qubits),
                    state: None,
                    contains_rz: false,
                    checked_len: 0,
                    checked_version: 0,
                    seq,
                });
                self.touch_stamp.push(0);
                let slot = self.slots.len() - 1;
                debug_assert!(slot < u32::MAX as usize, "window slots exhausted");
                slot as u32
            }
        }
    }

    /// Retire the window in `slot` and return the slot for reuse.
    ///
    /// The per-qubit index is left alone: marking the slot dead invalidates
    /// every registration naming it, and the next walk of each list drops what
    /// it finds (see the collection loop in [`SuperOpt::run`]). Eagerly
    /// unregistering meant a linear search of one list per window qubit on
    /// every window death, and a window dies for nearly every gate.
    fn retire(&mut self, slot: usize) {
        self.slots[slot].seq = Self::DEAD;
        // Clear the stamp so a future tenant of this slot can never inherit a
        // touch recorded against the previous one.
        self.touch_stamp[slot] = 0;
        self.free.push(slot as u32);
    }
}

impl SuperOpt {
    /// An optimizing pass: windows are checked against a MURM
    /// (built on first use per `murm_config`, then shared process-wide) and
    /// rewritten when a smaller equivalent exists.
    pub fn new(
        max_qubits: usize,
        window_gates: usize,
        murm_config: MurmConfig,
    ) -> Result<Self, SuperOptError> {
        Ok(Self::analyzer(max_qubits, window_gates).with_murm(shared_murm(murm_config)?))
    }

    /// An analysis-only instance: windows and their matrices are reported in
    /// `result.subcircuits`, but with no MURM nothing is rewritten.
    pub fn analyzer(max_qubits: usize, window_gates: usize) -> Self {
        Self {
            max_qubits,
            window_gates,
            collect_subcircuits: true,
            incremental: false,
            murm: None,
            store: Rc::default(),
            prev_input: Rc::default(),
        }
    }

    fn with_murm(mut self, murm: Arc<Murm>) -> Self {
        self.murm = Some(murm);
        self
    }

    /// Anchor windows only near gates that changed since the previous run's
    /// input, skipping regions whose windows were already analyzed and
    /// selected nothing. The output is identical to a full sweep as long as
    /// successive runs see successive versions of the same circuit (a
    /// sequential fixpoint driver); do not combine with parallel chunking,
    /// which feeds one instance interleaved unrelated circuits, or with
    /// subcircuit collection, which would come back incomplete.
    pub fn incremental(mut self) -> Self {
        self.incremental = true;
        self
    }

    /// Skip accumulating per-window [`SuperOptWindow`] diagnostics.
    ///
    /// Rewrites still apply; only `result.subcircuits` stays empty. Large
    /// circuits emit millions of windows, so optimization-only callers should
    /// opt out of retaining them.
    pub const fn without_subcircuits(mut self) -> Self {
        self.collect_subcircuits = false;
        self
    }

    /// Run one forward scan while maintaining one closed unitary component
    /// per anchor (see the module documentation for the algorithm).
    ///
    /// For every gate, in order: (1) windows touching a measurement or reset
    /// die; (2) every live window touching the gate is extended, re-closed,
    /// and analyzed; (3) the gate anchors a fresh window of its own.
    pub fn run(&self, circuit: &Circuit) -> Result<SuperOptResult, SuperOptError> {
        if self.window_gates == 0 {
            return Err(SuperOptError::ZeroWindowGates);
        }
        let gates_per_qubit = validate_circuit(circuit)?;

        let frontier = self.take_anchor_frontier(circuit);

        // Scan state: `arena` owns the live windows, recycling slots as they
        // die; `windows_by_qubit` inverts it so a gate finds the windows it
        // touches without scanning them all; `gates_by_qubit` is the full
        // per-qubit gate history that window closure consults.
        let mut arena = WindowArena::default();
        let mut windows_by_qubit: Vec<Vec<Registration>> = vec![Vec::new(); circuit.num_qubits];
        let mut gates_by_qubit: Vec<Vec<usize>> = gates_per_qubit
            .iter()
            .map(|&count| Vec::with_capacity(count))
            .collect();
        // Take the persistent store for the duration of this run; an early
        // error drops it, which only costs the next run a cold start.
        let mut store = MatrixStore::take_from(&self.store);
        let mut subcircuits = Vec::new();
        let mut rewrites = RewriteSet::new(circuit.gates.len());
        let mut touched_windows = Vec::new();

        for (gate_index, gate) in circuit.gates.iter().enumerate() {
            let gate_qubits = unique_qubits(gate);

            // Collect the windows this gate touches, deduped by slot via
            // `touch_stamp` (a window registered on two of the gate's qubits
            // appears in both lists) and ordered by anchor age, which is what
            // greedy rewrite selection consumes.
            touched_windows.clear();
            let stamp = gate_index + 1;
            for &qubit in &gate_qubits {
                // Walking a qubit's list also compacts it: a retired window is
                // not unregistered when it dies, only marked dead, so the
                // registrations to drop are exactly the ones whose recorded
                // anchor order no longer matches their slot's tenant.
                let live = &mut windows_by_qubit[qubit as usize];
                let mut kept = 0;
                for index in 0..live.len() {
                    let entry = live[index];
                    let slot = registered_slot(entry);
                    if arena.slots[slot].seq != registered_seq(entry) {
                        continue;
                    }
                    live[kept] = entry;
                    kept += 1;
                    // A window registered on two of the gate's qubits appears
                    // in both lists; stamping the slot dedups it.
                    if arena.touch_stamp[slot] != stamp {
                        arena.touch_stamp[slot] = stamp;
                        touched_windows.push(entry);
                    }
                }
                live.truncate(kept);
                gates_by_qubit[qubit as usize].push(gate_index);
            }
            // Registrations carry anchor order in their high bits, so this
            // sorts into the order greedy rewrite selection consumes without a
            // key closure — `sort_unstable_by_key` would re-evaluate its key
            // on every comparison — and without reading the arena at all.
            touched_windows.sort_unstable();

            // Rz is outside the exact Clifford+T matrix domain; measurement
            // and reset are non-unitary. All three terminate windows touching
            // their qubit. Keep the gate in per-qubit history so a disjoint
            // window cannot later bridge across this barrier.
            if matches!(gate, Gate::rz(..) | Gate::measure { .. } | Gate::reset(_)) {
                for &entry in &touched_windows {
                    arena.retire(registered_slot(entry));
                }
                continue;
            }

            // The gate's kind and operands, resolved once rather than once per
            // window it touches; only the support-local operand positions
            // differ per window. The barrier check above has already excluded
            // the one gate kind `gate_code` cannot encode.
            let gate_parts = gate_code(gate);

            for &entry in &touched_windows {
                let window_id = registered_slot(entry);
                if rewrites.is_claimed(gate_index)
                    || rewrites.window_claims_any(arena.window(window_id))
                {
                    arena.retire(window_id);
                    continue;
                }
                let (within_bounds, added_qubits) = expand_component_closure(
                    circuit,
                    arena.window(window_id),
                    gate_index,
                    &gate_qubits,
                    &gates_by_qubit,
                    self.max_qubits,
                    self.window_gates,
                );
                if !within_bounds {
                    arena.retire(window_id);
                    continue;
                }

                // An extension that added no qubit left every earlier member's
                // support-local encoding intact, so the store can reach this
                // window's new shape as a transition from its previous one. A
                // bridged-in qubit renumbers the support, so it cannot.
                let appended = added_qubits.is_empty().then_some(()).and(gate_parts);

                let at_gate_limit = arena.window(window_id).gate_indices.len() == self.window_gates;
                let selected = self.analyze_window(
                    circuit,
                    arena.window(window_id),
                    appended,
                    &mut store,
                    &mut rewrites,
                    &mut subcircuits,
                )?;

                if at_gate_limit || selected {
                    arena.retire(window_id);
                } else {
                    for &qubit in &added_qubits {
                        windows_by_qubit[qubit as usize].push(entry);
                    }
                }
            }

            // The current gate anchors a new one-gate closed component.
            if frontier.as_ref().is_none_or(|f| f[gate_index])
                && gate_qubits.len() <= self.max_qubits
                && !rewrites.is_claimed(gate_index)
            {
                // The window is built in its slot: anchoring builds no
                // canonical key at all now (the first touch resolves the
                // window's shape, and every touch after that steps the
                // transition cache from it), so there is nothing to compute
                // before a slot is available.
                let window_id = arena.anchor(gate_index, &gate_qubits);
                // A single non-identity Clifford+T gate cannot be rewritten to
                // the empty circuit. Analyze it only when diagnostics were
                // requested; Rz windows are rejected inside `analyze_window`.
                if self.collect_subcircuits {
                    self.analyze_window(
                        circuit,
                        arena.window(window_id as usize),
                        None,
                        &mut store,
                        &mut rewrites,
                        &mut subcircuits,
                    )?;
                }

                if self.window_gates > 1 && !rewrites.is_claimed(gate_index) {
                    let entry = register(arena.slots[window_id as usize].seq, window_id);
                    for &qubit in &gate_qubits {
                        windows_by_qubit[qubit as usize].push(entry);
                    }
                } else {
                    arena.retire(window_id as usize);
                }
            }
        }

        subcircuits.sort_by(|left, right| left.gate_indices.cmp(&right.gate_indices));
        let (optimized, rewrites) = rewrites.apply(circuit);
        let (cache_hits, cache_misses) = (store.hits, store.misses);
        store.store_back(&self.store);
        Ok(SuperOptResult {
            circuit: optimized,
            subcircuits,
            rewrites,
            cache_hits,
            cache_misses,
        })
    }

    /// In incremental mode, the bitmap of gates allowed to anchor a window
    /// (`None` anchors everywhere); also remembers this input for the next
    /// run's diff. See [`incremental::anchor_frontier`].
    fn take_anchor_frontier(&self, circuit: &Circuit) -> Option<Vec<bool>> {
        if !self.incremental {
            return None;
        }
        let mut prev = self.prev_input.borrow_mut();
        let frontier = incremental::anchor_frontier(circuit, prev.as_ref(), self.window_gates);
        *prev = Some(circuit.clone());
        frontier
    }

    /// Resolve one window emission: intern its matrix and synthesis outcome
    /// in the store, offer the outcome to the rewrite selection, and record
    /// the window when diagnostics are collected. Returns whether a rewrite
    /// was selected.
    /// `appended` is the code parts of the single gate this extension added,
    /// present when the extension left the window's support unchanged, which
    /// lets the store step its shape transition cache instead of rebuilding a
    /// canonical key; `None` forces a rebuild (a bridged-in qubit renumbered
    /// the support, or the window has no resolved shape yet).
    fn analyze_window(
        &self,
        circuit: &Circuit,
        window: &mut ActiveWindow,
        appended: Option<GateCode>,
        store: &mut MatrixStore,
        rewrites: &mut RewriteSet,
        subcircuits: &mut Vec<SuperOptWindow>,
    ) -> Result<bool, SuperOptError> {
        // Exact SuperOpt matrices represent Clifford+T only. Never construct
        // a matrix or accept a rewrite for a window containing Rz; phase
        // folding/decomposition is responsible for those rotations.
        // `contains_rz` is tracked incrementally by `expand_component_closure`,
        // so this is O(1) rather than a rescan of the (monotonically growing)
        // gate list on every touch.
        if window.contains_rz {
            return Ok(false);
        }
        let code = appended.and_then(|code| code.encode(&window.qubits));
        let resolved = store.resolve(
            window.state,
            code,
            circuit,
            &window.gate_indices,
            &window.qubits,
            self.murm.as_deref(),
        )?;
        window.state = resolved;
        let Some(entry_index) = resolved else {
            // The exact i8 numerator bound was exceeded. Leaving this window
            // untouched is conservative and preserves rewrite soundness.
            return Ok(false);
        };
        // Unspill the SmallVecs once; this runs per emitted window.
        let gate_indices: &[usize] = &window.gate_indices;
        let qubits: &[Qubit] = &window.qubits;
        let cached = store.entry(entry_index);
        let selected = rewrites.consider(cached, gate_indices, qubits);
        if self.collect_subcircuits {
            subcircuits.push(SuperOptWindow {
                gate_indices: gate_indices.to_vec(),
                qubits: qubits.to_vec(),
                matrix: Arc::clone(&cached.matrix),
            });
        }
        Ok(selected)
    }
}

impl Pass for SuperOpt {
    fn name(&self) -> &str {
        "SuperOpt"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        match SuperOpt::run(self, circuit) {
            Ok(result) => result.circuit,
            Err(error) => panic!("SuperOpt failed: {error}"),
        }
    }
}

/// Reject gates whose qubit operands fall outside the circuit, so the scan can
/// index per-qubit state unchecked — and, from the same walk, count how many
/// gates touch each qubit.
///
/// The counts size `gates_by_qubit` exactly. Growing those lists on demand
/// reallocated and recopied them all the way up to one entry per gate operand
/// on every sweep, which is the pass's largest remaining memory-traffic cost;
/// this walk already visits every gate and every operand, so the counts are
/// free.
fn validate_circuit(circuit: &Circuit) -> Result<Vec<usize>, SuperOptError> {
    let mut gates_per_qubit = vec![0usize; circuit.num_qubits];
    for (gate_index, gate) in circuit.gates.iter().enumerate() {
        for qubit in unique_qubits(gate) {
            if qubit as usize >= circuit.num_qubits {
                return Err(SuperOptError::InvalidQubit {
                    gate_index,
                    qubit,
                    num_qubits: circuit.num_qubits,
                });
            }
            gates_per_qubit[qubit as usize] += 1;
        }
    }
    Ok(gates_per_qubit)
}

/// Extend `window` with `current_gate` and re-close it, scanning only the
/// qubits this step newly introduces.
///
/// The window already holds the connected closure of its anchor over earlier
/// gates, and every gate touching a window qubit expands the window when it is
/// processed, so all history on qubits already in the window is present. Only a
/// gate that bridges in a *new* qubit requires scanning, so the BFS queue is
/// seeded with new qubits alone rather than the whole support.
///
/// Returns whether the window stayed within `max_gates` and `max_qubits`, and
/// the qubits inserted into `window.qubits` by this call (on both paths), none
/// of which have been registered yet.
fn expand_component_closure(
    circuit: &Circuit,
    window: &mut ActiveWindow,
    current_gate: usize,
    current_qubits: &[Qubit],
    gates_by_qubit: &[Vec<usize>],
    max_qubits: usize,
    max_gates: usize,
) -> (bool, QubitVec) {
    /// Admit `qubit` into the sorted support if absent, recording it in
    /// `added` and queueing it for a history scan. False when the support
    /// grows past `max_qubits`.
    fn admit(
        window: &mut ActiveWindow,
        added: &mut QubitVec,
        pending: &mut QubitVec,
        qubit: Qubit,
        max_qubits: usize,
    ) -> bool {
        if let Err(position) = window.qubits.binary_search(&qubit) {
            window.qubits.insert(position, qubit);
            added.push(qubit);
            pending.push(qubit);
            return window.qubits.len() <= max_qubits;
        }
        true
    }

    let mut added = QubitVec::new();
    let mut pending = QubitVec::new();

    let anchor = window.gate_indices[0];
    window.gate_indices.push(current_gate);
    if window.gate_indices.len() > max_gates {
        return (false, added);
    }

    for &qubit in current_qubits {
        if !admit(window, &mut added, &mut pending, qubit, max_qubits) {
            return (false, added);
        }
    }

    while let Some(qubit) = pending.pop() {
        let history = &gates_by_qubit[qubit as usize];
        let start = history.partition_point(|&gate_index| gate_index < anchor);
        for &gate_index in &history[start..] {
            if gate_index > current_gate {
                break;
            }
            if matches!(
                circuit.gates[gate_index],
                Gate::measure { .. } | Gate::reset(_)
            ) {
                return (false, added);
            }
            let Err(position) = window.gate_indices.binary_search(&gate_index) else {
                continue;
            };
            window.gate_indices.insert(position, gate_index);
            // A bridged-in gate lands anywhere in the list, so the "first
            // `checked_len` members are unclaimed" record no longer names a
            // prefix of the members it was taken over.
            window.checked_len = 0;
            if matches!(circuit.gates[gate_index], Gate::rz(..)) {
                window.contains_rz = true;
            }
            if window.gate_indices.len() > max_gates {
                return (false, added);
            }

            for gate_qubit in unique_qubits(&circuit.gates[gate_index]) {
                if !admit(window, &mut added, &mut pending, gate_qubit, max_qubits) {
                    return (false, added);
                }
            }
        }
    }

    (true, added)
}

/// Greedy, non-overlapping rewrite selection over the input's gate positions.
///
/// The first window offered with a strictly smaller replacement claims all
/// its gates; anything overlapping a claimed gate is refused afterwards, so
/// selected rewrites never conflict and [`RewriteSet::apply`] can splice them
/// all in a single reconstruction pass.
struct RewriteSet {
    claimed: Vec<bool>,
    /// Index into `selected` for the rewrite anchored at each gate position.
    anchored: Vec<Option<usize>>,
    selected: Vec<SuperOptRewrite>,
}

impl RewriteSet {
    fn new(gate_count: usize) -> Self {
        Self {
            claimed: vec![false; gate_count],
            anchored: vec![None; gate_count],
            selected: Vec::new(),
        }
    }

    fn is_claimed(&self, gate_index: usize) -> bool {
        self.claimed[gate_index]
    }

    fn claims_any(&self, gate_indices: &[usize]) -> bool {
        gate_indices.iter().any(|&index| self.claimed[index])
    }

    /// Whether any of `window`'s gates is already claimed, re-examining only
    /// the members added since the last time this window was cleared.
    ///
    /// Claimed bits are only ever set by [`RewriteSet::consider`], which also
    /// pushes to `selected` — so `selected.len()` is a version stamp for
    /// `claimed`, and a window cleared at the current version can trust every
    /// member it was cleared for. Selections are rare next to window touches,
    /// so this turns a per-touch scan of the whole (growing) gate list, at
    /// scattered indices into a circuit-sized bitmap, into a look at the one
    /// gate the touch just added. `expand_component_closure` resets
    /// `checked_len` when it inserts a bridged-in gate *before* the end, since
    /// then a prefix length no longer names the members already checked.
    fn window_claims_any(&self, window: &mut ActiveWindow) -> bool {
        let version = self.selected.len() as u32;
        let start = if window.checked_version == version {
            window.checked_len as usize
        } else {
            0
        };
        if window.gate_indices[start..]
            .iter()
            .any(|&index| self.claimed[index])
        {
            return true;
        }
        window.checked_version = version;
        window.checked_len = window.gate_indices.len() as u32;
        false
    }

    /// Claim this window's rewrite if it has a strictly smaller replacement
    /// and overlaps nothing already claimed. Returns whether it was selected.
    ///
    /// The replacement is stored on the window's *physical* qubits: the
    /// synthesized circuit uses support-local qubits `0..n`, and `qubits[q]`
    /// is exactly where local qubit `q` lives in the input circuit.
    fn consider(
        &mut self,
        cached: &CachedMatrix,
        gate_indices: &[usize],
        qubits: &[Qubit],
    ) -> bool {
        // Cheap rejections first. Most windows have no smaller equivalent at
        // all, and `claims_any` is a scan of the window's gate list indexing a
        // circuit-sized bitmap at scattered positions — running it before
        // knowing there is anything to select paid that scan on every emitted
        // window. It has no side effects, so the order is free to change.
        let Some(local) = cached.synthesized_replacement.as_ref() else {
            return false;
        };
        if local.len() >= gate_indices.len() {
            return false;
        }
        if self.claims_any(gate_indices) {
            return false;
        }

        let replacement: Vec<_> = local
            .iter()
            .map(|gate| gate.map_qubits(|q| qubits[q as usize]))
            .collect();
        for &index in gate_indices {
            self.claimed[index] = true;
        }
        self.anchored[gate_indices[0]] = Some(self.selected.len());
        self.selected.push(SuperOptRewrite {
            gate_indices: gate_indices.to_vec(),
            replacement,
        });
        true
    }

    /// Rebuild the circuit with every selected rewrite spliced in: each
    /// replacement is emitted at its window's anchor position, claimed gates
    /// are dropped, and everything else is copied through unchanged. Sound
    /// because a window's members all commute past the unclaimed gates
    /// between them (see the module documentation).
    fn apply(self, circuit: &Circuit) -> (Circuit, Vec<SuperOptRewrite>) {
        let mut optimized = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
        for (index, gate) in circuit.gates.iter().enumerate() {
            if let Some(rewrite) = self.anchored[index] {
                for gate in &self.selected[rewrite].replacement {
                    optimized.apply(gate.clone());
                }
            }
            if !self.claimed[index] {
                optimized.apply(gate.clone());
            }
        }
        (optimized, self.selected)
    }
}

/// A gate's qubit operands, sorted and deduplicated — the window code relies
/// on supports being sorted sets.
/// Sorts and dedups per arity rather than running the generic
/// `sort_unstable` + `dedup` over every arm: this is called once per gate of
/// every scan, and the one- and two-operand arms — everything left after
/// Toffoli decomposition — need no comparison sort at all, where the generic
/// path reached `slice::sort`'s insertion-sort small-case each time.
fn unique_qubits(gate: &Gate) -> QubitVec {
    match gate {
        Gate::x(q)
        | Gate::h(q)
        | Gate::s(q)
        | Gate::sdg(q)
        | Gate::z(q)
        | Gate::t(q)
        | Gate::tdg(q)
        | Gate::rz(_, q)
        | Gate::reset(q) => smallvec![*q],
        Gate::measure { qubit, .. } => smallvec![*qubit],
        Gate::cnot { control, target } | Gate::cz { control, target } => {
            match (*control).cmp(target) {
                std::cmp::Ordering::Less => smallvec![*control, *target],
                std::cmp::Ordering::Greater => smallvec![*target, *control],
                std::cmp::Ordering::Equal => smallvec![*control],
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
            let mut qubits: QubitVec = smallvec![*control1, *control2, *target];
            qubits.sort_unstable();
            qubits.dedup();
            qubits
        }
    }
}

#[cfg(test)]
mod tests;
