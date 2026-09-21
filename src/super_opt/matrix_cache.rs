//! Canonical window keys and the interned matrix store.
//!
//! Two windows with the same gate kinds in the same order, acting on the
//! same *support-local* qubit positions, have the same unitary — regardless
//! of which physical qubits or circuit positions they sit on. Keying on that
//! canonical shape lets one matrix construction and one MURM
//! probe serve every recurrence of a shape, across the whole circuit and
//! across runs.
//!
//! Keys come in two forms: a fixed-size, non-allocating [`CompactKey`] fast
//! path covering the common case (at most [`COMPACT_KEY_MAX_GATES`] gates on
//! at most four qubits), and a general, heap-allocated [`NormalizedGate`]
//! sequence for wider or longer Clifford+T windows.

use std::cell::RefCell;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::circuit::{Circuit, Gate, Qubit, canonical_ccx, canonical_pair, canonical_triple};

use super::SuperOptError;
use super::matrix::UnitaryMatrix;
use super::murm::Murm;

/// Everything the pass knows about one canonical window shape: its unitary,
/// and the smallest support-local replacement the MURM offered
/// (`None` caches the negative outcome, so misses are never retried).
#[derive(Clone, Debug)]
pub(super) struct CachedMatrix {
    pub(super) matrix: Arc<UnitaryMatrix>,
    pub(super) synthesized_replacement: Option<Vec<Gate>>,
}

/// Interned canonical-window matrices. A lookup returns an entry *index* and
/// builds its key, when it needs one at all, into a reused scratch buffer, so
/// the per-emission hot path never allocates on a cache hit.
///
/// Both key forms use support-local qubit indices and each entry is resolved
/// against the pass instance's fixed MURM, so entries are valid for
/// any window with the same shape — including windows in a different circuit.
/// That is what makes carrying the store across runs of one pass instance
/// sound (see [`MatrixStore::take_from`]).
#[derive(Debug, Default)]
pub(super) struct MatrixStore {
    // FxHash: this is probed once per emitted window (millions of times on
    // large circuits), and the keys are short gate sequences where SipHash's
    // per-lookup overhead dominates.
    cache: FxHashMap<Box<[NormalizedGate]>, usize>,
    compact_cache: FxHashMap<CompactKey, u32>,
    /// Successor cache over canonical shapes. A window that grows by appending
    /// one gate *without* changing its support keeps every earlier member's
    /// support-local encoding, so the grown window's canonical shape — and
    /// hence its interned entry — is a function of the pair (shape before,
    /// appended gate code) alone. Keying on `state << 16 | code` therefore
    /// answers the common extension from a two-integer probe, instead of
    /// rebuilding, rehashing, and `memcmp`ing a canonical key that is up to 64
    /// gates long. Entries are interned per distinct key, so an entry index
    /// names a shape uniquely, which is what makes the pair sufficient.
    transitions: FxHashMap<u64, u32>,
    entries: Vec<CachedMatrix>,
    scratch: Vec<NormalizedGate>,
    /// Reused buffer for the compact key built on a transition miss.
    scratch_key: CompactKey,
    pub(super) hits: usize,
    pub(super) misses: usize,
}

impl MatrixStore {
    /// Take the persistent store out of `slot` for one run, leaving the slot
    /// empty. The hit/miss counters describe a single run and are reset here.
    pub(super) fn take_from(slot: &RefCell<MatrixStore>) -> MatrixStore {
        let mut store = std::mem::take(&mut *slot.borrow_mut());
        store.hits = 0;
        store.misses = 0;
        store
    }

    /// Return the store to `slot` so the pass's next run starts warm. If two
    /// clones of the same pass instance (sharing this `Rc`) both ran and are
    /// now returning their store, the one with more interned entries wins
    /// and the other is simply dropped — clones only ever run one at a time
    /// on the single thread they're confined to, so this is an ordering
    /// choice, not a data race.
    pub(super) fn store_back(self, slot: &RefCell<MatrixStore>) {
        let mut guard = slot.borrow_mut();
        if guard.entries.len() < self.entries.len() {
            *guard = self;
        }
    }

    /// The interned entry for a window's canonical shape.
    ///
    /// `state` is the entry index for the window's shape *before* this
    /// extension and `code` the appended gate's encoding, both present only
    /// when the extension left the window's support unchanged; given the pair,
    /// the answer comes from `transitions` without a key being built at all.
    /// Otherwise the canonical key is rebuilt and interned as usual.
    ///
    /// Exactly one of `hits`/`misses` is counted per call, with the same
    /// meaning as before the transition cache existed: a transition hit can
    /// only have been recorded once the successor shape was interned, so it is
    /// a shape the store had already seen.
    pub(super) fn resolve(
        &mut self,
        state: Option<u32>,
        code: Option<u16>,
        circuit: &Circuit,
        gate_indices: &[usize],
        qubits: &[Qubit],
        murm: Option<&Murm>,
    ) -> Result<Option<u32>, SuperOptError> {
        let step = match (state, code) {
            (Some(state), Some(code)) => Some((u64::from(state) << 16) | u64::from(code)),
            _ => None,
        };
        if let Some(step) = step
            && let Some(&next) = self.transitions.get(&step)
        {
            self.hits += 1;
            return Ok(Some(next));
        }

        let resolved = self.intern(circuit, gate_indices, qubits, murm)?;
        if let (Some(step), Some(next)) = (step, resolved) {
            self.transitions.insert(step, next);
        }
        Ok(resolved)
    }

    /// The cached matrix and synthesis outcome for an index from [`Self::resolve`].
    pub(super) fn entry(&self, index: u32) -> &CachedMatrix {
        &self.entries[index as usize]
    }

    /// Build the window's canonical key, returning the index of its interned
    /// entry — or `None` when the exact matrix exceeded the representable
    /// numerator range, which leaves the window conservatively un-rewritten
    /// and un-interned.
    fn intern(
        &mut self,
        circuit: &Circuit,
        gate_indices: &[usize],
        qubits: &[Qubit],
        murm: Option<&Murm>,
    ) -> Result<Option<u32>, SuperOptError> {
        let compact =
            compact_normalized_key_into(circuit, gate_indices, qubits, &mut self.scratch_key);
        if compact {
            if let Some(&entry_index) = self.compact_cache.get(&self.scratch_key) {
                self.hits += 1;
                return Ok(Some(entry_index));
            }
        } else {
            normalized_gate_key(circuit, gate_indices, qubits, &mut self.scratch);
            if let Some(&entry_index) = self.cache.get(self.scratch.as_slice()) {
                self.hits += 1;
                return Ok(Some(entry_index as u32));
            }
        }

        self.misses += 1;
        let mut matrix = UnitaryMatrix::identity(qubits.len())?;
        for &gate_index in gate_indices {
            if matrix
                .apply_gate_left(&circuit.gates[gate_index], qubits)
                .is_err()
            {
                return Ok(None);
            }
        }
        let synthesized_replacement = murm.and_then(|murm| murm.synthesize(&matrix));
        let entry_index = self.entries.len();
        self.entries.push(CachedMatrix {
            synthesized_replacement,
            matrix: Arc::new(matrix),
        });
        if compact {
            let key = self.scratch_key;
            self.compact_cache.insert(key, entry_index as u32);
        } else {
            self.cache
                .insert(self.scratch.as_slice().into(), entry_index);
        }
        Ok(Some(entry_index as u32))
    }
}

/// One gate of a general canonical key on support-local qubit positions. The
/// fallback when the packed `u128` form below cannot represent the window.
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
enum NormalizedGate {
    X(usize),
    H(usize),
    S(usize),
    Sdg(usize),
    Z(usize),
    T(usize),
    Tdg(usize),
    Cnot(usize, usize),
    Cz(usize, usize),
    Ccx(usize, usize, usize),
    Ccz(usize, usize, usize),
}

/// Longest window the fixed-size compact key can represent without falling
/// back to the heap-allocated general key. Comfortably covers `-Osuper`'s
/// `window_gates=40` (see `main.rs`), with headroom for the hidden
/// `--superopt-window-gates` experimentation flag.
pub(super) const COMPACT_KEY_MAX_GATES: usize = 64;

/// Exact non-allocating normalized key for the common Clifford+T window: one
/// `u16` gate code (see [`GateCode`]) per slot, with only the first `len`
/// slots meaningful. Hashing and equality look at just that prefix (see the
/// trait impls below), so cost tracks the window's actual length rather than
/// the array's fixed capacity.
#[derive(Clone, Copy, Debug)]
pub(super) struct CompactKey {
    len: u8,
    gates: [u16; COMPACT_KEY_MAX_GATES],
}

impl Default for CompactKey {
    fn default() -> Self {
        Self::EMPTY
    }
}

impl CompactKey {
    const EMPTY: Self = Self {
        len: 0,
        gates: [0; COMPACT_KEY_MAX_GATES],
    };

    fn as_slice(&self) -> &[u16] {
        &self.gates[..usize::from(self.len)]
    }

    /// Append one gate's code in place. Returns `false` (leaving `self`
    /// unchanged) if the window would exceed `COMPACT_KEY_MAX_GATES` gates or
    /// the gate/support isn't compactly encodable; the caller then falls back
    /// to the general key.
    fn try_push(&mut self, gate: &Gate, support: &[Qubit]) -> bool {
        let len = usize::from(self.len);
        if len >= COMPACT_KEY_MAX_GATES {
            return false;
        }
        let Some(code) = gate_code(gate).and_then(|code| code.encode(support)) else {
            return false;
        };
        self.gates[len] = code;
        self.len = (len + 1) as u8;
        true
    }
}

impl PartialEq for CompactKey {
    fn eq(&self, other: &Self) -> bool {
        self.as_slice() == other.as_slice()
    }
}

impl Eq for CompactKey {}

impl Hash for CompactKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_slice().hash(state);
    }
}

/// Rebuild `out` as the canonical key of `gate_indices` over `support`,
/// reporting whether the window is compactly encodable at all. Writes in
/// place into a reused buffer, so a rebuild neither allocates nor zeroes the
/// key's fixed-size tail — only the `len` prefix is ever meaningful.
fn compact_normalized_key_into(
    circuit: &Circuit,
    gate_indices: &[usize],
    support: &[Qubit],
    out: &mut CompactKey,
) -> bool {
    out.len = 0;
    for &gate_index in gate_indices {
        if !out.try_push(&circuit.gates[gate_index], support) {
            return false;
        }
    }
    true
}

/// Owned-return form of [`compact_normalized_key_into`], for tests that assert
/// on encodability and key distinctness rather than running the scan.
#[cfg(test)]
pub(super) fn compact_normalized_key(
    circuit: &Circuit,
    gate_indices: &[usize],
    support: &[Qubit],
) -> Option<CompactKey> {
    let mut key = CompactKey::EMPTY;
    compact_normalized_key_into(circuit, gate_indices, support, &mut key).then_some(key)
}

/// A gate's compact encoding, split into the part that is the same for every
/// window — its kind tag and its operands, in operand order — and the part
/// that is not, each operand's position within the window's support.
///
/// The scan needs a gate's code once per live window the gate touches, a few
/// dozen times per gate, and the kind dispatch answers identically every time.
/// Resolving the kind once per gate leaves only the support-local positions on
/// the per-window path.
#[derive(Clone, Copy, Debug)]
pub(super) struct GateCode {
    tag: u16,
    operands: [Qubit; 3],
    arity: u8,
}

/// Split `gate` into its window-independent code parts. `None` for `rz`, which
/// a compact key cannot express — an Rz gate can reach a window through a
/// bridged-in qubit's history, so this is a real case, not an assertion.
pub(super) fn gate_code(gate: &Gate) -> Option<GateCode> {
    let (tag, operands, arity) = match gate {
        Gate::x(q) => (0, [*q, 0, 0], 1),
        Gate::h(q) => (1, [*q, 0, 0], 1),
        Gate::s(q) => (2, [*q, 0, 0], 1),
        Gate::sdg(q) => (3, [*q, 0, 0], 1),
        Gate::z(q) => (4, [*q, 0, 0], 1),
        Gate::t(q) => (5, [*q, 0, 0], 1),
        Gate::tdg(q) => (6, [*q, 0, 0], 1),
        Gate::cnot { control, target } => (7, [*control, *target, 0], 2),
        Gate::cz { control, target } => {
            let (a, b) = canonical_pair(*control, *target);
            (8, [a, b, 0], 2)
        }
        Gate::ccx {
            control1,
            control2,
            target,
        } => {
            let (a, b, target) = canonical_ccx(*control1, *control2, *target);
            (9, [a, b, target], 3)
        }
        Gate::ccz {
            control1,
            control2,
            target,
        } => {
            let (a, b, c) = canonical_triple(*control1, *control2, *target);
            (10, [a, b, c], 3)
        }
        Gate::rz(..) => return None,
        Gate::measure { .. } | Gate::reset(_) => {
            unreachable!("measurement and reset are window barriers")
        }
    };
    Some(GateCode {
        tag,
        operands,
        arity,
    })
}

impl GateCode {
    /// This gate's code over `support`. `None` when the support is wider than
    /// the two-bit operand positions can address — wider analyzer-only windows
    /// must use the general normalized key rather than silently alias local
    /// qubits.
    pub(super) fn encode(self, support: &[Qubit]) -> Option<u16> {
        if support.len() > 4 {
            return None;
        }
        let mut code = self.tag;
        for index in 0..usize::from(self.arity) {
            // `support` is sorted, so an operand's local position is just how
            // many support qubits come before it: a count over at most four
            // entries, with no branching and no panicking path, rather than a
            // `binary_search(..).expect(..)`.
            let operand = self.operands[index];
            let mut position = 0u16;
            for &qubit in support {
                position += u16::from(qubit < operand);
            }
            code |= position << (4 + 2 * index);
        }
        Some(code)
    }
}

/// Write the window's general canonical key into `key` (reused scratch, so
/// the hot path allocates nothing on a hit).
fn normalized_gate_key(
    circuit: &Circuit,
    gate_indices: &[usize],
    support: &[Qubit],
    key: &mut Vec<NormalizedGate>,
) {
    let local = |q| {
        support
            .binary_search(&q)
            .expect("window qubit is in support")
    };
    key.clear();
    key.extend(
        gate_indices
            .iter()
            .map(|&gate_index| match &circuit.gates[gate_index] {
                Gate::x(q) => NormalizedGate::X(local(*q)),
                Gate::h(q) => NormalizedGate::H(local(*q)),
                Gate::s(q) => NormalizedGate::S(local(*q)),
                Gate::sdg(q) => NormalizedGate::Sdg(local(*q)),
                Gate::z(q) => NormalizedGate::Z(local(*q)),
                Gate::t(q) => NormalizedGate::T(local(*q)),
                Gate::tdg(q) => NormalizedGate::Tdg(local(*q)),
                Gate::rz(..) => unreachable!("Rz gates are SuperOpt window barriers"),
                Gate::cnot { control, target } => {
                    NormalizedGate::Cnot(local(*control), local(*target))
                }
                Gate::cz { control, target } => {
                    let (a, b) = canonical_pair(local(*control), local(*target));
                    NormalizedGate::Cz(a, b)
                }
                Gate::ccx {
                    control1,
                    control2,
                    target,
                } => {
                    let (a, b, target) =
                        canonical_ccx(local(*control1), local(*control2), local(*target));
                    NormalizedGate::Ccx(a, b, target)
                }
                Gate::ccz {
                    control1,
                    control2,
                    target,
                } => {
                    let (a, b, c) =
                        canonical_triple(local(*control1), local(*control2), local(*target));
                    NormalizedGate::Ccz(a, b, c)
                }
                Gate::measure { .. } | Gate::reset(_) => {
                    unreachable!("measurement and reset are SuperOpt window barriers")
                }
            }),
    );
}
