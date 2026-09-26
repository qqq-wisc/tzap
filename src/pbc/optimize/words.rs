//! Packed canonical Pauli words and the evaluator that produces them.
//!
//! A canonical word on n qubits is stored as two bit planes of `l` u64 words
//! each, x then z: I = (0,0), X = (1,0), Z = (0,1), Y = (1,1). Qubit q is bit
//! q % 64 of word q / 64. A word carries no phase; a value `i^p W` keeps its
//! phase exponent `p` (mod 4) beside the planes.

use std::hash::{Hash, Hasher};

use rustc_hash::{FxHashMap, FxHasher};

use super::super::pauli::{PauliArena, PauliNode, PauliRef};
use super::super::{Pauli, PbcError};

/// Whether two canonical words anticommute: the parity of positions where
/// both are non-identity and differ.
pub(super) fn anticommutes(p: &[u64], q: &[u64], l: usize) -> bool {
    let (px, pz) = p.split_at(l);
    let (qx, qz) = q.split_at(l);
    let mut acc = 0;
    for j in 0..l {
        acc ^= (px[j] & qz[j]) ^ (pz[j] & qx[j]);
    }
    acc.count_ones() & 1 != 0
}

/// Write the canonical word of `p * q` into `out` and return its phase
/// exponent `e`, so that `p * q = i^e * out`.
///
/// With `W = i^nu(W) X^x Z^z` and `nu(W) = popcount(x & z)`,
/// `e = nu(p) + nu(q) - nu(out) + 2 popcount(z_p & x_q) (mod 4)`.
pub(super) fn product(p: &[u64], q: &[u64], out: &mut [u64], l: usize) -> u8 {
    let (px, pz) = p.split_at(l);
    let (qx, qz) = q.split_at(l);
    let (ox, oz) = out.split_at_mut(l);
    let (mut nu_p, mut nu_q, mut nu_out, mut cross) = (0u32, 0u32, 0u32, 0u64);
    for j in 0..l {
        ox[j] = px[j] ^ qx[j];
        oz[j] = pz[j] ^ qz[j];
        nu_p += (px[j] & pz[j]).count_ones();
        nu_q += (qx[j] & qz[j]).count_ones();
        nu_out += (ox[j] & oz[j]).count_ones();
        cross ^= pz[j] & qx[j];
    }
    let e = nu_p + nu_q + 4 - (nu_out % 4) + 2 * (cross.count_ones() & 1);
    (e % 4) as u8
}

fn hash_words(words: &[u64]) -> u64 {
    let mut hasher = FxHasher::default();
    words.hash(&mut hasher);
    hasher.finish()
}

/// Interned canonical axes. Equal words always get the same id. Each packed
/// id keeps an arena handle that evaluates to `+W`; ids synthesized during
/// optimization have none until one is built for output.
pub(super) struct Axes {
    pub l: usize,
    words: Vec<u64>,
    pub handles: Vec<Option<PauliRef>>,
    supports: Vec<u64>,
    /// Hash -> most recently interned id with that hash; `chain` links older
    /// ids with the same hash, so collisions are resolved by exact comparison.
    index: FxHashMap<u64, u32>,
    chain: Vec<u32>,
}

const NONE: u32 = u32::MAX;

impl Axes {
    fn new(l: usize) -> Self {
        Self {
            l,
            words: Vec::new(),
            handles: Vec::new(),
            supports: Vec::new(),
            index: FxHashMap::default(),
            chain: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.handles.len()
    }

    pub fn get(&self, id: u32) -> &[u64] {
        let start = id as usize * 2 * self.l;
        &self.words[start..start + 2 * self.l]
    }

    /// Packed storage in use, in u64 words.
    pub fn words_len(&self) -> usize {
        self.words.len()
    }

    /// Number of non-identity single-qubit factors.
    pub fn weight(&self, id: u32) -> usize {
        let (x, z) = self.get(id).split_at(self.l);
        x.iter()
            .zip(z)
            .map(|(x, z)| (x | z).count_ones() as usize)
            .sum()
    }

    pub fn is_identity(&self, id: u32) -> bool {
        self.get(id).iter().all(|&w| w == 0)
    }

    pub fn anticommute(&self, a: u32, b: u32) -> bool {
        anticommutes(self.get(a), self.get(b), self.l)
    }

    /// Signature of the axis's support: bit q / l for every qubit q it acts
    /// on, so each bit covers l consecutive qubits (exactly one when n <= 64).
    /// Contiguous ranges suit circuits whose gates act on nearby qubits.
    pub fn support(&self, id: u32) -> u64 {
        self.supports[id as usize]
    }

    pub fn signature(words: &[u64], l: usize) -> u64 {
        let (x, z) = words.split_at(l);
        let mut signature = 0;
        for (j, (x, z)) in x.iter().zip(z).enumerate() {
            let mut w = x | z;
            while w != 0 {
                let q = 64 * j + w.trailing_zeros() as usize;
                signature |= 1 << (q / l);
                w &= w - 1;
            }
        }
        signature
    }

    pub fn intern(&mut self, words: &[u64], handle: Option<PauliRef>) -> u32 {
        let hash = hash_words(words);
        let mut id = self.index.get(&hash).copied().unwrap_or(NONE);
        while id != NONE {
            if self.get(id) == words {
                return id;
            }
            id = self.chain[id as usize];
        }
        let id = self.handles.len() as u32;
        self.words.extend_from_slice(words);
        self.handles.push(handle);
        self.supports.push(Self::signature(words, self.l));
        self.chain.push(self.index.insert(hash, id).unwrap_or(NONE));
        id
    }
}

/// A rotation record: an interned canonical axis and a signed angle `k`,
/// normalized to -3..=4, for `exp(-i k pi/8 W)`. `support` is the axis's
/// signature (`Axes::support`): disjoint signatures prove that two axes
/// commute without touching their packed words.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Rot {
    pub axis: u32,
    pub k: i8,
    pub support: u64,
}

pub(super) fn normalize(k: i32) -> i8 {
    let k = k.rem_euclid(8) as i8;
    if k > 4 { k - 8 } else { k }
}

/// Pack each root expression. Returns the interned axes and, per root, its
/// axis id and sign: the root evaluates to `sign * W`.
///
/// Nodes are evaluated once, in arena order (children are always older), and
/// each value is freed after its last use. `max_words` bounds the live and
/// interned packed storage, in u64 words.
pub(super) fn pack(
    arena: &PauliArena,
    num_qubits: usize,
    root_refs: &[PauliRef],
    max_words: usize,
) -> Result<(Axes, Vec<(u32, i8)>), PbcError> {
    let l = num_qubits.div_ceil(64).max(1);
    let stride = 2 * l;
    let mut axes = Axes::new(l);
    let mut records = vec![(0, 1); root_refs.len()];

    // Roots sorted by node, so each is packed as soon as its node is evaluated.
    let mut roots: Vec<(usize, usize)> = root_refs
        .iter()
        .enumerate()
        .map(|(i, r)| (r.node, i))
        .collect();
    roots.sort_unstable();

    // Use counts over the reachable subgraph: one per root reference and one
    // per product edge out of a reachable node.
    let mut uses = vec![0u32; arena.nodes.len()];
    let mut stack: Vec<usize> = Vec::new();
    for &(node, _) in &roots {
        if uses[node] == 0 {
            stack.push(node);
        }
        uses[node] += 1;
        while let Some(n) = stack.pop() {
            if let PauliNode::Product(a, b) = arena.nodes[n] {
                for child in [a.node, b.node] {
                    if uses[child] == 0 {
                        stack.push(child);
                    }
                    uses[child] += 1;
                }
            }
        }
    }

    // Slot storage for live node values: planes plus a phase exponent.
    let mut slots: Vec<u64> = Vec::new();
    let mut phases: Vec<u8> = Vec::new();
    let mut free: Vec<u32> = Vec::new();
    let mut node_slot = vec![NONE; arena.nodes.len()];
    let mut live = 0usize;
    let mut scratch = vec![0u64; stride];

    let check = |live: usize, axes: &Axes| {
        let words = (live + axes.len())
            .checked_mul(stride)
            .ok_or(PbcError::ExpansionLimit)?;
        if words > max_words {
            return Err(PbcError::ExpansionLimit);
        }
        Ok(())
    };

    let mut next_root = 0;
    for node in 0..arena.nodes.len() {
        if uses[node] == 0 {
            continue;
        }
        // Evaluate into scratch, then move into a slot.
        let phase = match arena.nodes[node] {
            PauliNode::Identity => {
                scratch.fill(0);
                0
            }
            PauliNode::Single { qubit, pauli } => {
                scratch.fill(0);
                let (word, bit) = (qubit as usize / 64, 1u64 << (qubit % 64));
                if matches!(pauli, Pauli::X | Pauli::Y) {
                    scratch[word] |= bit;
                }
                if matches!(pauli, Pauli::Z | Pauli::Y) {
                    scratch[l + word] |= bit;
                }
                0
            }
            PauliNode::Product(a, b) => {
                let (sa, sb) = (node_slot[a.node] as usize, node_slot[b.node] as usize);
                let e = product(
                    &slots[sa * stride..(sa + 1) * stride],
                    &slots[sb * stride..(sb + 1) * stride],
                    &mut scratch,
                    l,
                );
                let p = e as usize
                    + phases[sa] as usize
                    + a.phase as usize
                    + phases[sb] as usize
                    + b.phase as usize;
                for child in [a.node, b.node] {
                    uses[child] -= 1;
                    if uses[child] == 0 {
                        free.push(node_slot[child]);
                        node_slot[child] = NONE;
                        live -= 1;
                    }
                }
                (p % 4) as u8
            }
        };
        let slot = match free.pop() {
            Some(slot) => slot,
            None => {
                check(live + 1, &axes)?;
                slots.resize(slots.len() + stride, 0);
                phases.push(0);
                phases.len() as u32 - 1
            }
        };
        let s = slot as usize;
        slots[s * stride..(s + 1) * stride].copy_from_slice(&scratch);
        phases[s] = phase;
        node_slot[node] = slot;
        live += 1;

        // Pack the operations rooted here.
        while next_root < roots.len() && roots[next_root].0 == node {
            let root = roots[next_root].1;
            next_root += 1;
            let reference = root_refs[root];
            let total = (phase as usize + reference.phase as usize) % 4;
            let (sign, handle) = match total {
                0 => (1, reference),
                2 => (-1, reference.scaled(super::super::Phase::MinusOne)),
                _ => return Err(PbcError::NonHermitianAxis),
            };
            check(live, &axes)?;
            let id = axes.intern(&slots[s * stride..(s + 1) * stride], Some(handle));
            records[root] = (id, sign);
            uses[node] -= 1;
            if uses[node] == 0 {
                free.push(slot);
                node_slot[node] = NONE;
                live -= 1;
            }
        }
    }
    Ok((axes, records))
}

#[cfg(test)]
mod tests;
