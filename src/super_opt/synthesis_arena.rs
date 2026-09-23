//! Prefix-sharing storage for the MURM's circuits.
//!
//! Breadth-first enumeration only ever extends an existing circuit by one
//! gate, so the hundreds of thousands of stored circuits form a tree: each
//! node records just its final gate and its parent. A full circuit is
//! recovered by walking to the root and reversing — done only on a MURM
//! hit, never during enumeration.

use std::io::{self, Read, Write};

use rustc_hash::FxHashMap;

use crate::circuit::GateSet;

use super::matrix::UnitaryFingerprint;
use super::murm::LibraryGate;

/// Sentinel gate tag marking the root node (no gate, no parent) in the
/// on-disk format; distinct from every real `LibraryGate::to_bytes` tag.
const NO_GATE_TAG: u8 = 255;

/// One stored circuit: its last gate plus the node holding the rest. The
/// root (the empty circuit, i.e. the identity) has neither.
#[derive(Clone, Copy, Debug)]
pub(super) struct CircuitNode {
    pub(super) parent: Option<usize>,
    pub(super) gate: Option<LibraryGate>,
}

/// One MURM width stored as a prefix-sharing circuit arena.
#[derive(Clone, Debug, Default)]
pub(super) struct WidthTable {
    fingerprints: FxHashMap<UnitaryFingerprint, usize>,
    pub(super) nodes: Vec<CircuitNode>,
}

impl WidthTable {
    pub(super) fn with_identity(fingerprint: UnitaryFingerprint) -> Self {
        let mut table = Self::default();
        table.nodes.push(CircuitNode {
            parent: None,
            gate: None,
        });
        table.fingerprints.insert(fingerprint, 0);
        table
    }

    pub(super) fn len(&self) -> usize {
        self.fingerprints.len()
    }

    pub(super) fn contains_key(&self, fingerprint: &UnitaryFingerprint) -> bool {
        self.fingerprints.contains_key(fingerprint)
    }

    pub(super) fn node_for(&self, fingerprint: &UnitaryFingerprint) -> Option<usize> {
        self.fingerprints.get(fingerprint).copied()
    }

    #[cfg(test)]
    pub(super) fn insert_fingerprint_alias(
        &mut self,
        fingerprint: UnitaryFingerprint,
        node: usize,
    ) {
        self.fingerprints.insert(fingerprint, node);
    }

    pub(super) fn insert_child(
        &mut self,
        fingerprint: UnitaryFingerprint,
        parent: usize,
        gate: LibraryGate,
    ) -> usize {
        let node = self.nodes.len();
        self.nodes.push(CircuitNode {
            parent: Some(parent),
            gate: Some(gate),
        });
        self.fingerprints.insert(fingerprint, node);
        node
    }

    pub(super) fn circuit(&self, mut node: usize) -> Vec<LibraryGate> {
        let mut circuit = Vec::new();
        loop {
            let entry = self.nodes[node];
            if let Some(gate) = entry.gate {
                circuit.push(gate);
            }
            let Some(parent) = entry.parent else {
                break;
            };
            node = parent;
        }
        circuit.reverse();
        circuit
    }

    /// Persist this width's table verbatim: each node's fingerprint (8
    /// bytes), parent index (`u32`, `u32::MAX` for none), and gate (4-byte
    /// encoding). Node order is preserved so parent indices stay valid on
    /// read; the fingerprint map is rebuilt from the same pairs on load,
    /// so no BFS re-derivation is needed.
    pub(super) fn write_to(&self, out: &mut impl Write) -> io::Result<()> {
        out.write_all(&(self.nodes.len() as u64).to_le_bytes())?;
        let mut fingerprint_by_node = vec![None; self.nodes.len()];
        for (&fingerprint, &node) in &self.fingerprints {
            fingerprint_by_node[node] = Some(fingerprint);
        }
        for (node, fingerprint) in self.nodes.iter().zip(&fingerprint_by_node) {
            let fingerprint = fingerprint
                .expect("every stored node was inserted alongside exactly one fingerprint");
            out.write_all(&fingerprint.to_bits().to_le_bytes())?;
            let parent = node.parent.map_or(u32::MAX, |p| p as u32);
            out.write_all(&parent.to_le_bytes())?;
            let gate_bytes = node
                .gate
                .map_or([NO_GATE_TAG, 0, 0, 0], LibraryGate::to_bytes);
            out.write_all(&gate_bytes)?;
        }
        Ok(())
    }

    pub(super) fn read_from(
        input: &mut impl Read,
        num_qubits: usize,
        max_nodes: usize,
        basis: GateSet,
        identity_fingerprint: Option<UnitaryFingerprint>,
    ) -> io::Result<Self> {
        let invalid = |msg: &str| io::Error::new(io::ErrorKind::InvalidData, msg.to_owned());
        let mut len_buf = [0u8; 8];
        input.read_exact(&mut len_buf)?;
        let len = usize::try_from(u64::from_le_bytes(len_buf))
            .map_err(|_| invalid("MURM width table length does not fit usize"))?;
        if len > max_nodes {
            return Err(invalid(
                "MURM width table exceeds its configured or file-size bound",
            ));
        }
        if (num_qubits == 0 && len != 0) || (num_qubits > 0 && len == 0) {
            return Err(invalid("MURM width table has an invalid root count"));
        }

        let mut nodes = Vec::with_capacity(len);
        let mut fingerprints = FxHashMap::with_capacity_and_hasher(len, Default::default());
        for index in 0..len {
            let mut fingerprint_buf = [0u8; 8];
            input.read_exact(&mut fingerprint_buf)?;
            let fingerprint = UnitaryFingerprint::from_bits(u64::from_le_bytes(fingerprint_buf));

            let mut parent_buf = [0u8; 4];
            input.read_exact(&mut parent_buf)?;
            let parent_raw = u32::from_le_bytes(parent_buf);
            let parent = (parent_raw != u32::MAX).then_some(parent_raw as usize);

            let mut gate_buf = [0u8; 4];
            input.read_exact(&mut gate_buf)?;
            let gate = if gate_buf == [NO_GATE_TAG, 0, 0, 0] {
                None
            } else if gate_buf[0] == NO_GATE_TAG {
                return Err(invalid("invalid MURM root-gate encoding"));
            } else {
                let gate = LibraryGate::from_bytes(gate_buf).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "unknown library gate tag")
                })?;
                if gate.to_bytes() != gate_buf {
                    return Err(invalid("non-canonical MURM library gate encoding"));
                }
                Some(gate)
            };

            if index == 0 {
                if parent.is_some() || gate.is_some() || Some(fingerprint) != identity_fingerprint {
                    return Err(invalid("MURM width table has an invalid identity root"));
                }
            } else {
                let Some(parent) = parent else {
                    return Err(invalid("non-root MURM node has no parent"));
                };
                let Some(gate) = gate else {
                    return Err(invalid("non-root MURM node has no gate"));
                };
                if parent >= index {
                    return Err(invalid("MURM parent must precede its child"));
                }
                if !gate.is_valid_for(num_qubits, basis) {
                    return Err(invalid("MURM gate is outside its width or synthesis basis"));
                }
            }

            if fingerprints.insert(fingerprint, index).is_some() {
                return Err(invalid("duplicate MURM fingerprint"));
            }
            nodes.push(CircuitNode { parent, gate });
        }
        Ok(Self {
            fingerprints,
            nodes,
        })
    }
}
