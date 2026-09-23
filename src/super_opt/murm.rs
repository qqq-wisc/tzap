//! A bounded minimal unitary representative map (MURM): breadth-first
//! enumeration of configured library-gate circuits keyed by unitary
//! fingerprint, plus the process-wide cache that shares built MURMs across
//! passes.
//!
//! For each width the enumeration grows circuits one gate at a time, layer
//! by layer, recording each unitary the first time it appears. Because
//! layers are visited in gate-count order, the first circuit to reach a
//! unitary is a smallest one — so a MURM hit *is* the synthesis answer, no
//! search needed at lookup time. Two prunes keep the frontier tractable
//! without losing any unitary: a child never follows its parent's inverse
//! (the product would revisit the grandparent's unitary), and among
//! qubit-disjoint neighbors only the canonically ordered interleaving is
//! expanded (the swapped one has the same product).

use std::collections::HashMap;
use std::io::{self, Read, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use rayon::prelude::*;

use crate::circuit::{
    Gate, GateKind, GateSet, Qubit, canonical_ccx, canonical_pair, canonical_triple,
};

use super::matrix::{UnitaryFingerprint, UnitaryMatrix, unitary_fingerprint};
use super::synthesis_arena::WidthTable;
use super::{MurmConfig, SUPPORTED_GATE_SET, SuperOptError};

/// On-disk MURM cache format identifier and version. Bump
/// `CACHE_FORMAT_VERSION` whenever the byte layout below changes; a mismatch
/// (or a missing/corrupt file) simply falls back to rebuilding, never to a
/// misread MURM.
///
/// The crate version is checked too (see `CACHE_CRATE_VERSION`), separately
/// from the byte layout: MURM *construction* (pruning rules, the library
/// gate set, etc.) can change between releases without touching how a MURM
/// is serialized, and such a change must still invalidate old caches even
/// though `CACHE_FORMAT_VERSION` didn't move. Tying invalidation to the crate
/// version means that never has to be caught by hand — every release gets a
/// fresh cache namespace for free.
const CACHE_MAGIC: &[u8; 4] = b"MURM";
// Version 2 introduced exact cyclotomic fingerprints, version 3 the compact
// i8 coefficient bound, version 4 stores 64-bit rather than 128-bit keys,
// version 5 stores the configurable gate basis, and version 6 validates the
// complete body with a checksum and strict structural bounds.
const CACHE_FORMAT_VERSION: u32 = 6;
const CACHE_CRATE_VERSION: &str = env!("CARGO_PKG_VERSION");
static TEMP_FILE_SEQUENCE: AtomicU64 = AtomicU64::new(0);
const SERIALIZED_NODE_BYTES: u64 = 16;
const CHECKSUM_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
const CHECKSUM_PRIME: u64 = 0x0000_0100_0000_01b3;

fn update_checksum(mut checksum: u64, bytes: &[u8]) -> u64 {
    for &byte in bytes {
        checksum ^= u64::from(byte);
        checksum = checksum.wrapping_mul(CHECKSUM_PRIME);
    }
    checksum
}

struct ChecksummedWriter<W> {
    inner: W,
    checksum: u64,
}

impl<W> ChecksummedWriter<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            checksum: CHECKSUM_OFFSET_BASIS,
        }
    }

    fn checksum(&self) -> u64 {
        self.checksum
    }
}

impl<W: Write> Write for ChecksummedWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let written = self.inner.write(bytes)?;
        self.checksum = update_checksum(self.checksum, &bytes[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

struct ChecksummedReader<R> {
    inner: R,
    checksum: u64,
}

impl<R> ChecksummedReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            checksum: CHECKSUM_OFFSET_BASIS,
        }
    }

    fn checksum(&self) -> u64 {
        self.checksum
    }
}

impl<R: Read> Read for ChecksummedReader<R> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(bytes)?;
        self.checksum = update_checksum(self.checksum, &bytes[..read]);
        Ok(read)
    }
}

/// Reads and validates a cache file's header — magic, format version, and
/// config fields — against `config`. Shared by `read_from_disk` (which
/// continues on to read the MURM body) and `disk_cache_exists` (which only
/// needs to know the header matches).
fn read_cache_header(input: &mut impl Read, config: MurmConfig) -> io::Result<()> {
    let invalid = |msg: &str| io::Error::new(io::ErrorKind::InvalidData, msg.to_owned());

    let mut magic = [0u8; 4];
    input.read_exact(&mut magic)?;
    if magic != *CACHE_MAGIC {
        return Err(invalid("not a SuperOpt MURM cache file"));
    }
    let mut version_buf = [0u8; 4];
    input.read_exact(&mut version_buf)?;
    if u32::from_le_bytes(version_buf) != CACHE_FORMAT_VERSION {
        return Err(invalid("cache format version mismatch"));
    }

    let mut crate_version_len = [0u8; 1];
    input.read_exact(&mut crate_version_len)?;
    let mut crate_version_buf = vec![0u8; crate_version_len[0] as usize];
    input.read_exact(&mut crate_version_buf)?;
    if crate_version_buf != CACHE_CRATE_VERSION.as_bytes() {
        return Err(invalid("cache crate version mismatch"));
    }

    let mut qubits_buf = [0u8; 4];
    input.read_exact(&mut qubits_buf)?;
    let mut gates_buf = [0u8; 4];
    input.read_exact(&mut gates_buf)?;
    let mut entries_buf = [0u8; 8];
    input.read_exact(&mut entries_buf)?;
    let mut basis_buf = [0u8; 2];
    input.read_exact(&mut basis_buf)?;
    let stored_entries = u64::from_le_bytes(entries_buf);
    if usize::try_from(u32::from_le_bytes(qubits_buf)).ok() != Some(config.max_qubits)
        || usize::try_from(u32::from_le_bytes(gates_buf)).ok() != Some(config.max_gates)
        || u64::try_from(config.max_entries_per_qubit).ok() != Some(stored_entries)
        || u16::from_le_bytes(basis_buf) != config.basis.bits()
    {
        return Err(invalid("cache config mismatch"));
    }
    Ok(())
}

/// A gate the MURM enumerates over support-local qubits.
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum LibraryGate {
    H(u8),
    X(u8),
    Z(u8),
    S(u8),
    Sdg(u8),
    T(u8),
    Tdg(u8),
    Cnot(u8, u8),
    Cz(u8, u8),
    Ccx(u8, u8, u8),
    Ccz(u8, u8, u8),
}

impl LibraryGate {
    fn cz(a: u8, b: u8) -> Self {
        let (a, b) = canonical_pair(a, b);
        Self::Cz(a, b)
    }

    fn ccx(control1: u8, control2: u8, target: u8) -> Self {
        let (control1, control2, target) = canonical_ccx(control1, control2, target);
        Self::Ccx(control1, control2, target)
    }

    fn ccz(a: u8, b: u8, c: u8) -> Self {
        let (a, b, c) = canonical_triple(a, b, c);
        Self::Ccz(a, b, c)
    }

    pub(super) fn to_gate(self) -> Gate {
        match self {
            Self::X(q) => Gate::x(q.into()),
            Self::H(q) => Gate::h(q.into()),
            Self::S(q) => Gate::s(q.into()),
            Self::Sdg(q) => Gate::sdg(q.into()),
            Self::Z(q) => Gate::z(q.into()),
            Self::T(q) => Gate::t(q.into()),
            Self::Tdg(q) => Gate::tdg(q.into()),
            Self::Cnot(control, target) => Gate::cnot {
                control: control.into(),
                target: target.into(),
            },
            Self::Cz(control, target) => Gate::cz {
                control: control.into(),
                target: target.into(),
            },
            Self::Ccx(control1, control2, target) => Gate::ccx {
                control1: control1.into(),
                control2: control2.into(),
                target: target.into(),
            },
            Self::Ccz(control1, control2, target) => Gate::ccz {
                control1: control1.into(),
                control2: control2.into(),
                target: target.into(),
            },
        }
    }

    pub(super) fn qubits(self) -> [Option<u8>; 3] {
        match self {
            Self::X(q)
            | Self::H(q)
            | Self::S(q)
            | Self::Sdg(q)
            | Self::Z(q)
            | Self::T(q)
            | Self::Tdg(q) => [Some(q), None, None],
            Self::Cnot(left, right) | Self::Cz(left, right) => [Some(left), Some(right), None],
            Self::Ccx(a, b, c) | Self::Ccz(a, b, c) => [Some(a), Some(b), Some(c)],
        }
    }

    /// Whether a serialized gate is a canonical member of `basis` over the
    /// given support width. Cache files are inputs, so do not assume their
    /// operands are in range or distinct merely because the tag is known.
    pub(super) fn is_valid_for(self, num_qubits: usize, basis: GateSet) -> bool {
        if !basis.contains(self.to_gate().kind()) {
            return false;
        }
        let in_range = |q: u8| usize::from(q) < num_qubits;
        match self {
            Self::H(q)
            | Self::X(q)
            | Self::Z(q)
            | Self::S(q)
            | Self::Sdg(q)
            | Self::T(q)
            | Self::Tdg(q) => in_range(q),
            Self::Cnot(control, target) => {
                control != target && in_range(control) && in_range(target)
            }
            Self::Cz(a, b) => a < b && in_range(a) && in_range(b),
            Self::Ccx(control1, control2, target) => {
                control1 < control2
                    && control1 != target
                    && control2 != target
                    && in_range(control1)
                    && in_range(control2)
                    && in_range(target)
            }
            Self::Ccz(a, b, c) => a < b && b < c && in_range(a) && in_range(b) && in_range(c),
        }
    }

    pub(super) fn is_disjoint(self, other: Self) -> bool {
        let left = self.qubits();
        let right = other.qubits();
        left.into_iter()
            .flatten()
            .all(|qubit| !right.contains(&Some(qubit)))
    }

    pub(super) fn is_inverse_of(self, other: Self) -> bool {
        match (self, other) {
            (Self::S(q), Self::Sdg(r))
            | (Self::Sdg(q), Self::S(r))
            | (Self::T(q), Self::Tdg(r))
            | (Self::Tdg(q), Self::T(r)) => q == r,
            _ => {
                self == other
                    && matches!(
                        self,
                        Self::X(_)
                            | Self::H(_)
                            | Self::Z(_)
                            | Self::Cnot(..)
                            | Self::Cz(..)
                            | Self::Ccx(..)
                            | Self::Ccz(..)
                    )
            }
        }
    }

    /// Fixed 4-byte encoding: a tag plus up to three qubit operands.
    pub(super) fn to_bytes(self) -> [u8; 4] {
        match self {
            Self::H(q) => [0, q, 0, 0],
            Self::X(q) => [1, q, 0, 0],
            Self::Z(q) => [2, q, 0, 0],
            Self::S(q) => [3, q, 0, 0],
            Self::Sdg(q) => [4, q, 0, 0],
            Self::T(q) => [5, q, 0, 0],
            Self::Tdg(q) => [6, q, 0, 0],
            Self::Cnot(control, target) => [7, control, target, 0],
            Self::Cz(a, b) => [8, a, b, 0],
            Self::Ccx(a, b, c) => [9, a, b, c],
            Self::Ccz(a, b, c) => [10, a, b, c],
        }
    }

    pub(super) fn from_bytes(bytes: [u8; 4]) -> Option<Self> {
        let [tag, a, b, c] = bytes;
        Some(match tag {
            0 => Self::H(a),
            1 => Self::X(a),
            2 => Self::Z(a),
            3 => Self::S(a),
            4 => Self::Sdg(a),
            5 => Self::T(a),
            6 => Self::Tdg(a),
            7 => Self::Cnot(a, b),
            8 => Self::cz(a, b),
            9 => Self::ccx(a, b, c),
            10 => Self::ccz(a, b, c),
            _ => return None,
        })
    }
}

/// Breadth-first map from a unitary fingerprint to the smallest circuit found.
#[derive(Clone, Debug)]
pub(super) struct Murm {
    // Only `entries` serves lookups; the rest is bookkeeping read by tests.
    entries: Vec<WidthTable>,
    #[cfg_attr(not(test), allow(dead_code))]
    saturated: Vec<bool>,
    #[cfg_attr(not(test), allow(dead_code))]
    completed_depth: Vec<usize>,
}

impl Murm {
    pub(super) fn build(config: MurmConfig) -> Result<Self, SuperOptError> {
        if !(1..=5).contains(&config.max_qubits) {
            return Err(SuperOptError::InvalidMurmConfig {
                reason: format!("max_qubits must be in 1..=5, got {}", config.max_qubits),
            });
        }
        if config.max_entries_per_qubit == 0 {
            return Err(SuperOptError::InvalidMurmConfig {
                reason: "max_entries_per_qubit must be greater than zero".to_owned(),
            });
        }
        if config.basis.is_empty() || !config.basis.is_subset(SUPPORTED_GATE_SET) {
            return Err(SuperOptError::InvalidMurmConfig {
                reason: format!("basis must be a non-empty subset of {}", SUPPORTED_GATE_SET),
            });
        }

        let mut entries = vec![WidthTable::default(); config.max_qubits + 1];
        let mut saturated = vec![false; config.max_qubits + 1];
        let mut completed_depth = vec![0; config.max_qubits + 1];
        for num_qubits in 1..=config.max_qubits {
            let identity = UnitaryMatrix::identity(num_qubits)?;
            entries[num_qubits] = WidthTable::with_identity(unitary_fingerprint(&identity));
            let gates = library_gates(num_qubits, config.basis);
            // A valid basis can have no generators at a smaller width (for
            // example, the explicit basis `ccx` at widths one and two). The
            // identity-only MURM for that width is complete at depth zero.
            if gates.is_empty() {
                continue;
            }
            let support: Vec<Qubit> = (0..num_qubits as Qubit).collect();
            let mut frontier = vec![(0, identity)];

            // Parents per parallel batch: enough candidates to spread across
            // threads while a batch's survivor list stays small in memory.
            let batch_parents = (65_536 / gates.len()).max(1);

            'depths: for depth in 1..=config.max_gates {
                // Accepted children this layer as (frontier position, node, gate).
                let mut accepted = Vec::new();
                for (batch_index, batch) in frontier.chunks(batch_parents).enumerate() {
                    // Matrix products and fingerprints dominate the build, so
                    // candidates are generated in parallel against a read-only
                    // view of the table. Survivors are then inserted serially
                    // in enumeration order, which keeps the table (and the
                    // exact saturation point) identical to a sequential build;
                    // candidates already present at batch start would have been
                    // skipped by the sequential scan too, so pre-filtering them
                    // in the parallel phase changes nothing.
                    let table = &entries[num_qubits];
                    let batch_survivors: Vec<Vec<(LibraryGate, UnitaryFingerprint)>> = batch
                        .par_iter()
                        .map(|(parent, base)| {
                            let last = table.nodes[*parent].gate;
                            let mut scratch = base.clone();
                            let mut survivors = Vec::new();
                            for &gate in &gates {
                                if let Some(last) = last
                                    && (last.is_inverse_of(gate)
                                        || (last.is_disjoint(gate) && gate < last))
                                {
                                    continue;
                                }
                                scratch.copy_from(base);
                                if scratch.apply_gate_left(&gate.to_gate(), &support).is_err() {
                                    continue;
                                }
                                let fingerprint = unitary_fingerprint(&scratch);
                                if !table.contains_key(&fingerprint) {
                                    survivors.push((gate, fingerprint));
                                }
                            }
                            survivors
                        })
                        .collect();

                    let table = &mut entries[num_qubits];
                    for (offset, survivors) in batch_survivors.into_iter().enumerate() {
                        let position = batch_index * batch_parents + offset;
                        let parent = frontier[position].0;
                        for (gate, fingerprint) in survivors {
                            if table.contains_key(&fingerprint) {
                                continue;
                            }
                            if table.len() >= config.max_entries_per_qubit {
                                saturated[num_qubits] = true;
                                break 'depths;
                            }
                            let node = table.insert_child(fingerprint, parent, gate);
                            accepted.push((position, node, gate));
                        }
                    }
                }
                completed_depth[num_qubits] = depth;
                if accepted.is_empty() {
                    break;
                }
                // Re-deriving each child from its parent repeats one gate
                // application per accepted node, in exchange for never holding
                // matrices for the (mostly duplicate) rejected candidates.
                let next_frontier = accepted
                    .into_par_iter()
                    .map(|(position, node, gate)| {
                        let mut matrix = frontier[position].1.clone();
                        matrix
                            .apply_gate_left(&gate.to_gate(), &support)
                            .expect("an accepted table child remains representable");
                        (node, matrix)
                    })
                    .collect();
                frontier = next_frontier;
            }
        }

        Ok(Self {
            entries,
            saturated,
            completed_depth,
        })
    }

    /// Write this MURM to `path` for reuse by a later process, tagged with
    /// `config` so a mismatched config on read is rejected rather than
    /// silently misinterpreted. Written to a sibling temp file and renamed
    /// into place, so a reader never observes a partially written cache file
    /// (concurrent writers each rename their own complete file; last one
    /// wins, which is fine since every writer for the same `config` builds
    /// byte-identical content).
    pub(super) fn write_to_disk(
        &self,
        path: &std::path::Path,
        config: MurmConfig,
    ) -> io::Result<()> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let sequence = TEMP_FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let file_name = path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("murm");
        let tmp_path = path.with_file_name(format!(
            ".{file_name}.tmp.{}.{}",
            std::process::id(),
            sequence
        ));
        let write_result = (|| {
            let file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&tmp_path)?;
            let mut out = io::BufWriter::new(file);
            out.write_all(CACHE_MAGIC)?;
            out.write_all(&CACHE_FORMAT_VERSION.to_le_bytes())?;
            out.write_all(&[CACHE_CRATE_VERSION.len() as u8])?;
            out.write_all(CACHE_CRATE_VERSION.as_bytes())?;
            out.write_all(&(config.max_qubits as u32).to_le_bytes())?;
            out.write_all(&(config.max_gates as u32).to_le_bytes())?;
            out.write_all(&(config.max_entries_per_qubit as u64).to_le_bytes())?;
            out.write_all(&config.basis.bits().to_le_bytes())?;
            let checksum = {
                let mut body = ChecksummedWriter::new(&mut out);
                body.write_all(&(self.entries.len() as u32).to_le_bytes())?;
                for width_table in &self.entries {
                    width_table.write_to(&mut body)?;
                }
                for &saturated in &self.saturated {
                    body.write_all(&[u8::from(saturated)])?;
                }
                for &depth in &self.completed_depth {
                    body.write_all(&(depth as u32).to_le_bytes())?;
                }
                body.checksum()
            };
            out.write_all(&checksum.to_le_bytes())?;
            out.flush()?;
            out.get_ref().sync_all()
        })();
        if let Err(error) = write_result {
            let _ = std::fs::remove_file(&tmp_path);
            return Err(error);
        }
        match std::fs::rename(&tmp_path, path) {
            Ok(()) => Ok(()),
            // On platforms where rename cannot replace a file, another
            // process winning the same deterministic cache race is success
            // only when the published file is complete and valid.
            Err(_) if Self::read_from_disk(path, config).is_ok() => {
                let _ = std::fs::remove_file(&tmp_path);
                Ok(())
            }
            Err(first_error) => {
                let retry = if path.exists() {
                    std::fs::remove_file(path).and_then(|()| std::fs::rename(&tmp_path, path))
                } else {
                    Err(first_error)
                };
                if retry.is_err() {
                    let _ = std::fs::remove_file(&tmp_path);
                }
                retry
            }
        }
    }

    /// Read a MURM previously written by `write_to_disk`, rejecting it
    /// (with an `io::Error`) unless its header matches `config` exactly and
    /// the format version is one this build understands. Any error here —
    /// missing file, truncated write, config mismatch, version bump — should
    /// be treated by the caller as "no usable cache", not as a hard failure.
    pub(super) fn read_from_disk(path: &std::path::Path, config: MurmConfig) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let file_len = file.metadata()?.len();
        let mut input = io::BufReader::new(file);
        read_cache_header(&mut input, config)?;
        let invalid = |msg: &str| io::Error::new(io::ErrorKind::InvalidData, msg.to_owned());

        let (entries, saturated, completed_depth, computed_checksum) = {
            let mut body = ChecksummedReader::new(&mut input);
            let mut width_count_buf = [0u8; 4];
            body.read_exact(&mut width_count_buf)?;
            let width_count = usize::try_from(u32::from_le_bytes(width_count_buf))
                .map_err(|_| invalid("MURM width count does not fit usize"))?;
            if width_count != config.max_qubits + 1 {
                return Err(invalid("MURM width count does not match max_qubits"));
            }

            // A node occupies exactly 16 bytes on disk. Bounding the sum by
            // both the configured entry cap and the actual file size prevents
            // corrupt length fields from driving disproportionate allocation.
            let mut remaining_nodes =
                usize::try_from(file_len / SERIALIZED_NODE_BYTES).unwrap_or(usize::MAX);
            let mut entries = Vec::with_capacity(width_count);
            for num_qubits in 0..width_count {
                let identity = if num_qubits == 0 {
                    None
                } else {
                    let matrix = UnitaryMatrix::identity(num_qubits)
                        .map_err(|error| invalid(&format!("invalid MURM width: {error}")))?;
                    Some(unitary_fingerprint(&matrix))
                };
                let max_nodes = if num_qubits == 0 {
                    0
                } else {
                    config.max_entries_per_qubit.min(remaining_nodes)
                };
                let table = WidthTable::read_from(
                    &mut body,
                    num_qubits,
                    max_nodes,
                    config.basis,
                    identity,
                )?;
                remaining_nodes = remaining_nodes.saturating_sub(table.nodes.len());
                entries.push(table);
            }

            let mut saturated = Vec::with_capacity(width_count);
            for num_qubits in 0..width_count {
                let mut byte = [0u8; 1];
                body.read_exact(&mut byte)?;
                if byte[0] > 1 || (num_qubits == 0 && byte[0] != 0) {
                    return Err(invalid("invalid MURM saturation marker"));
                }
                saturated.push(byte[0] != 0);
            }
            let mut completed_depth = Vec::with_capacity(width_count);
            for num_qubits in 0..width_count {
                let mut depth_buf = [0u8; 4];
                body.read_exact(&mut depth_buf)?;
                let depth = usize::try_from(u32::from_le_bytes(depth_buf))
                    .map_err(|_| invalid("MURM depth does not fit usize"))?;
                if depth > config.max_gates || (num_qubits == 0 && depth != 0) {
                    return Err(invalid("invalid MURM completed depth"));
                }
                completed_depth.push(depth);
            }
            (entries, saturated, completed_depth, body.checksum())
        };

        let mut checksum_buf = [0u8; 8];
        input.read_exact(&mut checksum_buf)?;
        if u64::from_le_bytes(checksum_buf) != computed_checksum {
            return Err(invalid("MURM body checksum mismatch"));
        }
        let mut trailing = [0u8; 1];
        if input.read(&mut trailing)? != 0 {
            return Err(invalid("trailing bytes after MURM cache body"));
        }

        Ok(Self {
            entries,
            saturated,
            completed_depth,
        })
    }

    #[cfg(test)]
    pub(super) fn entry_count(&self, num_qubits: usize) -> usize {
        self.entries.get(num_qubits).map_or(0, WidthTable::len)
    }

    #[cfg(test)]
    pub(super) fn is_saturated(&self, num_qubits: usize) -> bool {
        self.saturated.get(num_qubits).copied().unwrap_or(false)
    }

    /// Largest gate count whose entire breadth-first layer was enumerated.
    #[cfg(test)]
    pub(super) fn completed_depth(&self, num_qubits: usize) -> usize {
        self.completed_depth.get(num_qubits).copied().unwrap_or(0)
    }

    /// Test seam for exercising the release-mode fingerprint collision guard.
    #[cfg(test)]
    pub(super) fn inject_fingerprint_alias(
        &mut self,
        query: &UnitaryMatrix,
        wrong_candidate: &UnitaryMatrix,
    ) {
        assert_eq!(query.num_qubits(), wrong_candidate.num_qubits());
        let width = query.num_qubits();
        let candidate = self.entries[width]
            .node_for(&unitary_fingerprint(wrong_candidate))
            .expect("wrong candidate is present in the test table");
        self.entries[width].insert_fingerprint_alias(unitary_fingerprint(query), candidate);
    }

    /// A smallest known library circuit implementing `matrix` up to global
    /// phase, on local qubits `0..matrix.num_qubits()`.
    pub(super) fn synthesize(&self, matrix: &UnitaryMatrix) -> Option<Vec<Gate>> {
        let table = self.entries.get(matrix.num_qubits())?;
        let node = table.node_for(&unitary_fingerprint(matrix))?;
        let circuit = table.circuit(node);
        // A fingerprint is still a finite hash of the exact matrix. This
        // exact comparison is the release-mode collision guard that makes
        // accepting a rewrite sound; it is not a redundant post-rewrite audit.
        let candidate = library_circuit_matrix(matrix.num_qubits(), &circuit)
            .ok()
            .flatten()?;
        matrix
            .equivalent_up_to_global_phase(&candidate)
            .then(|| circuit.into_iter().map(LibraryGate::to_gate).collect())
    }
}

type SharedMurm = Result<Arc<Murm>, SuperOptError>;
type MurmCache = HashMap<MurmConfig, SharedMurm>;

pub(super) fn shared_murm(config: MurmConfig) -> SharedMurm {
    static MURMS: OnceLock<Mutex<MurmCache>> = OnceLock::new();

    let murms = MURMS.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(murm) = lock_cache(murms).get(&config) {
        return murm.clone();
    }

    // Built with the cache lock *released*, which is load-bearing rather than
    // merely tidy: `Murm::build` is rayon-parallel, so holding
    // a process-wide mutex across it deadlocks. A rayon worker running an
    // unrelated stolen job can call back in here — the map-reduce path builds
    // a `SuperOpt` per chunk, and `SuperOpt::new` is public, so any caller can
    // do this from inside their own parallel iterator — and block on the lock,
    // while the builder blocks waiting for workers to drain its own parallel
    // iterators. Once every worker is parked on the lock, nothing can make
    // progress. Two threads racing the same cold config may now each build a
    // MURM, but only the first insertion is kept, so all callers still share
    // one `Arc` per config.
    let murm = build_or_load_from_disk(config).map(Arc::new);
    lock_cache(murms).entry(config).or_insert(murm).clone()
}

fn lock_cache(murms: &Mutex<MurmCache>) -> std::sync::MutexGuard<'_, MurmCache> {
    murms
        .lock()
        .expect("SuperOpt MURM cache mutex was poisoned")
}

/// Process-wide override for the cache root, set once at startup by
/// `super_opt::set_cache_dir` (the CLI's `--cache-dir`). Takes precedence
/// over every environment-derived location below.
static CACHE_ROOT_OVERRIDE: OnceLock<std::path::PathBuf> = OnceLock::new();

/// Basename of the MURM subdirectory inside whichever cache root wins.
const MURM_SUBDIR: &str = "murm";

/// Install a process-wide cache root, overriding `$TZAP_CACHE_DIR` and the
/// XDG lookup. Returns `Err` with the already-installed root if one was set
/// before — callers set this once, before any MURM is built or loaded, so a
/// second call would silently mean half a run read from one directory and
/// half from another.
pub(super) fn set_cache_root(dir: std::path::PathBuf) -> Result<(), std::path::PathBuf> {
    CACHE_ROOT_OVERRIDE.set(dir).map_err(|_| {
        CACHE_ROOT_OVERRIDE
            .get()
            .cloned()
            .unwrap_or_else(std::path::PathBuf::new)
    })
}

/// tzap's cache root: `--cache-dir`, then `$TZAP_CACHE_DIR`, then the XDG
/// Base Directory locations (`$XDG_CACHE_HOME/tzap`, then the spec's default
/// of `$HOME/.cache/tzap`), then the Windows ones (`%LOCALAPPDATA%\tzap`,
/// then `%USERPROFILE%\.cache\tzap`). `None` when none of them resolve, in
/// which case callers just skip disk caching entirely — it is always a pure
/// speed optimization, never required for correctness.
///
/// The Windows names are read on every platform rather than behind
/// `cfg(windows)`, so the resolution order is one list that any platform's
/// tests can exercise; nothing sets them on a Unix machine. They come last
/// because `$HOME` is what a Unix-shell environment on Windows (MSYS,
/// Git Bash) provides, and a user who has been caching there should keep
/// reading the same MURMs from either shell. Without them a native Windows
/// process — which gets no `$HOME` — cached nothing at all and rebuilt its
/// MURM on every run.
///
/// An empty `$XDG_CACHE_HOME` counts as unset, per the spec, and a relative
/// one is ignored the same way: the spec requires absolute paths, and
/// honoring a relative one would scatter caches wherever tzap happened to be
/// invoked from.
pub(super) fn cache_root() -> Option<std::path::PathBuf> {
    if let Some(dir) = CACHE_ROOT_OVERRIDE.get() {
        return Some(dir.clone());
    }
    if let Some(dir) = non_empty_env("TZAP_CACHE_DIR") {
        return Some(std::path::PathBuf::from(dir));
    }
    if let Some(dir) = non_empty_env("XDG_CACHE_HOME") {
        let dir = std::path::PathBuf::from(dir);
        if dir.is_absolute() {
            return Some(dir.join("tzap"));
        }
    }
    if let Some(home) = non_empty_env("HOME") {
        return Some(std::path::Path::new(&home).join(".cache").join("tzap"));
    }
    if let Some(local) = non_empty_env("LOCALAPPDATA") {
        return Some(std::path::Path::new(&local).join("tzap"));
    }
    let profile = non_empty_env("USERPROFILE")?;
    Some(std::path::Path::new(&profile).join(".cache").join("tzap"))
}

fn non_empty_env(name: &str) -> Option<std::ffi::OsString> {
    std::env::var_os(name).filter(|value| !value.is_empty())
}

/// Directory holding on-disk MURM caches.
pub(super) fn cache_dir() -> Option<std::path::PathBuf> {
    Some(cache_root()?.join(MURM_SUBDIR))
}

/// Filename for `config`'s MURM. One file per distinct `config`, since
/// different bounds produce different MURMs; the format version and crate
/// version are in the name too so a bump of either can't collide with (and
/// doesn't need to explicitly invalidate) old files — they're simply never
/// looked up again and can be cleaned up with `tzap --clear-cache`.
fn cache_file_name(config: MurmConfig) -> String {
    format!(
        "q{}_g{}_e{}_b{:04x}.v{CACHE_FORMAT_VERSION}.tzap{CACHE_CRATE_VERSION}.bin",
        config.max_qubits,
        config.max_gates,
        config.max_entries_per_qubit,
        config.basis.bits(),
    )
}

/// Where `config`'s MURM is written, and the first place it's looked for.
fn cache_file_path(config: MurmConfig) -> Option<std::path::PathBuf> {
    Some(cache_dir()?.join(cache_file_name(config)))
}

/// Whether a matching on-disk cache candidate for `config` exists, checked by
/// reading just the header rather than loading and validating the full MURM
/// body twice. Purely informational — for callers wanting to report whether a
/// `SuperOpt::new` call is likely to do a fast cache load or a slow fresh
/// build. `build_or_load_from_disk` is the sole source of truth and validates
/// the full body independently, so corruption or a race can only make this
/// hint optimistic; it can never cause a bad load.
pub(super) fn disk_cache_exists(config: MurmConfig) -> bool {
    existing_cache_file(config).is_some()
}

/// A readable cache candidate for `config`, matched by its header without the
/// more expensive full-body validation performed by `read_from_disk`.
fn existing_cache_file(config: MurmConfig) -> Option<std::path::PathBuf> {
    let path = cache_file_path(config)?;
    std::fs::File::open(&path)
        .is_ok_and(|file| read_cache_header(&mut io::BufReader::new(file), config).is_ok())
        .then_some(path)
}

/// Size in bytes of the on-disk cache file for `config`, if a valid one
/// exists — a reporting aid alongside `disk_cache_exists`, never
/// load-bearing for correctness.
pub(super) fn disk_cache_size_bytes(config: MurmConfig) -> Option<u64> {
    let path = existing_cache_file(config)?;
    std::fs::metadata(&path).ok().map(|metadata| metadata.len())
}

/// Load a matching MURM from disk if one exists, else build it and (on a
/// best-effort basis) write it back for the next process to reuse. A disk
/// read/write failure of any kind — missing file, corrupt content, a
/// read-only cache directory — is swallowed here: it can only make this run
/// as slow as a cold run would have been anyway, never wrong.
fn build_or_load_from_disk(config: MurmConfig) -> Result<Murm, SuperOptError> {
    if let Some(path) = cache_file_path(config)
        && let Ok(murm) = Murm::read_from_disk(&path, config)
    {
        return Ok(murm);
    }

    let murm = Murm::build(config)?;
    if let Some(path) = cache_file_path(config) {
        let _ = murm.write_to_disk(&path, config);
    }
    Ok(murm)
}

/// Build a library circuit's exact matrix. The outer `Result` reports an
/// unusably large dense matrix; `Ok(None)` means its coefficients exceeded the
/// bounded i8 representation and the candidate must be skipped.
pub(super) fn library_circuit_matrix(
    num_qubits: usize,
    circuit: &[LibraryGate],
) -> Result<Option<UnitaryMatrix>, SuperOptError> {
    let support: Vec<Qubit> = (0..num_qubits as Qubit).collect();
    let mut matrix = UnitaryMatrix::identity(num_qubits)?;
    for &gate in circuit {
        if matrix.apply_gate_left(&gate.to_gate(), &support).is_err() {
            return Ok(None);
        }
    }
    Ok(Some(matrix))
}

pub(super) fn library_gates(num_qubits: usize, basis: GateSet) -> Vec<LibraryGate> {
    let mut gates = Vec::new();
    for kind in GateKind::ALL {
        if !basis.contains(kind) {
            continue;
        }
        match kind {
            GateKind::H
            | GateKind::X
            | GateKind::Z
            | GateKind::S
            | GateKind::Sdg
            | GateKind::T
            | GateKind::Tdg => {
                for q in 0..num_qubits as u8 {
                    gates.push(match kind {
                        GateKind::H => LibraryGate::H(q),
                        GateKind::X => LibraryGate::X(q),
                        GateKind::Z => LibraryGate::Z(q),
                        GateKind::S => LibraryGate::S(q),
                        GateKind::Sdg => LibraryGate::Sdg(q),
                        GateKind::T => LibraryGate::T(q),
                        GateKind::Tdg => LibraryGate::Tdg(q),
                        _ => unreachable!(),
                    });
                }
            }
            GateKind::Cx => {
                for control in 0..num_qubits as u8 {
                    for target in 0..num_qubits as u8 {
                        if control != target {
                            gates.push(LibraryGate::Cnot(control, target));
                        }
                    }
                }
            }
            GateKind::Cz => {
                for a in 0..num_qubits as u8 {
                    for b in a + 1..num_qubits as u8 {
                        gates.push(LibraryGate::cz(a, b));
                    }
                }
            }
            GateKind::Ccx => {
                for a in 0..num_qubits as u8 {
                    for b in a + 1..num_qubits as u8 {
                        for target in 0..num_qubits as u8 {
                            if target != a && target != b {
                                gates.push(LibraryGate::ccx(a, b, target));
                            }
                        }
                    }
                }
            }
            GateKind::Ccz => {
                for a in 0..num_qubits as u8 {
                    for b in a + 1..num_qubits as u8 {
                        for c in b + 1..num_qubits as u8 {
                            gates.push(LibraryGate::ccz(a, b, c));
                        }
                    }
                }
            }
            _ => unreachable!("unsupported MURM gate kind"),
        }
    }
    gates
}
