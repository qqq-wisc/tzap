//! Demand-driven sparse evaluation. No allocation is indexed by total qubits
//! or arena length. Dead child maps are moved into parents, not copied.
use super::*;
use std::collections::{HashMap, hash_map::Entry};

#[derive(Clone, Debug)]
pub(crate) struct SparsePauli {
    pub phase: Phase,
    pub factors: HashMap<Qubit, Pauli>,
}

impl SparsePauli {
    /// Fixed-width radix ordering keeps large sparse output linear in support.
    pub fn sorted_factors(&self) -> Vec<(Qubit, Pauli)> {
        let mut factors: Vec<_> = self.factors.iter().map(|(&q, &p)| (q, p)).collect();
        if factors.len() <= 32 {
            factors.sort_unstable_by_key(|&(q, _)| q);
            return factors;
        }
        let mut scratch = vec![(0, Pauli::I); factors.len()];
        for shift in [0, 8, 16, 24] {
            let mut offsets = [0usize; 256];
            for &(q, _) in &factors {
                offsets[((q >> shift) & 255) as usize] += 1;
            }
            let mut total = 0;
            for count in &mut offsets {
                let next = total + *count;
                *count = total;
                total = next;
            }
            for &(q, p) in &factors {
                let offset = &mut offsets[((q >> shift) & 255) as usize];
                scratch[*offset] = (q, p);
                *offset += 1;
            }
            std::mem::swap(&mut factors, &mut scratch);
        }
        factors
    }
}

/// Deterministic work counts, also used to test scaling without timing noise.
#[derive(Default, Debug)]
pub(crate) struct ExpansionStats {
    pub nodes: usize,
    pub copied_factors: usize,
    pub multiplied_factors: usize,
    pub output_factors: usize,
    pub work: usize,
}

impl ExpansionStats {
    fn charge(&mut self, amount: usize, limit: usize) -> Result<(), PbcError> {
        self.work = self
            .work
            .checked_add(amount)
            .ok_or(PbcError::ExpansionLimit)?;
        if self.work > limit {
            return Err(PbcError::ExpansionLimit);
        }
        Ok(())
    }
}

impl PauliArena {
    /// Evaluate reachable nodes once and visit roots in the supplied order.
    /// Expected work is O(reachable nodes + roots + copied/merged/output
    /// factors). Shared large intermediates can still make this superlinear in
    /// arena size; max_work bounds that work before each allocation or merge.
    pub fn materialize(
        &self,
        roots: &[PauliRef],
        max_work: usize,
        mut visit: impl FnMut(PauliRef, &SparsePauli) -> Result<(), PbcError>,
    ) -> Result<ExpansionStats, PbcError> {
        let mut stats = ExpansionStats::default();
        let mut uses = HashMap::<usize, usize>::new();
        let mut pending = Vec::new();
        // Count graph edges once, plus each external root use. This enables
        // last-use ownership transfer even with diamonds and repeated roots.
        for &root in roots {
            self.check(root)?;
            stats.charge(1, max_work)?;
            *uses.entry(root.node).or_default() += 1;
            pending.push(root.node);
        }
        let mut discovered = std::collections::HashSet::new();
        while let Some(node) = pending.pop() {
            if !discovered.insert(node) {
                continue;
            }
            stats.charge(1, max_work)?;
            stats.nodes += 1;
            if let PauliNode::Product(a, b) = self.nodes[node] {
                for child in [a, b] {
                    *uses.entry(child.node).or_default() += 1;
                    pending.push(child.node);
                }
            }
        }
        drop(discovered);
        let mut values = HashMap::<usize, SparsePauli>::new();
        for &root in roots {
            pending.push(root.node);
            while let Some(&node) = pending.last() {
                if values.contains_key(&node) {
                    pending.pop();
                    continue;
                }
                let value = match self.nodes[node] {
                    PauliNode::Identity => SparsePauli {
                        phase: Phase::One,
                        factors: HashMap::new(),
                    },
                    PauliNode::Single { qubit, pauli } => SparsePauli {
                        phase: Phase::One,
                        factors: if pauli == Pauli::I {
                            HashMap::new()
                        } else {
                            HashMap::from([(qubit, pauli)])
                        },
                    },
                    PauliNode::Product(a, b) => {
                        if !values.contains_key(&a.node) {
                            pending.push(a.node);
                            continue;
                        }
                        if !values.contains_key(&b.node) {
                            pending.push(b.node);
                            continue;
                        }
                        let mut left = take(a.node, &mut uses, &mut values, &mut stats, max_work)?;
                        let mut right = take(b.node, &mut uses, &mut values, &mut stats, max_work)?;
                        let mut phase = left.phase.times(a.phase).times(right.phase).times(b.phase);
                        let swapped = left.factors.len() < right.factors.len();
                        if swapped {
                            std::mem::swap(&mut left, &mut right);
                        }
                        stats.charge(right.factors.len(), max_work)?;
                        stats.multiplied_factors += right.factors.len();
                        for (q, other) in right.factors {
                            match left.factors.entry(q) {
                                Entry::Vacant(entry) => {
                                    entry.insert(other);
                                }
                                Entry::Occupied(mut entry) => {
                                    let base = *entry.get();
                                    let (local, factor) = if swapped {
                                        other.times(base)
                                    } else {
                                        base.times(other)
                                    };
                                    phase = phase.times(local);
                                    if factor == Pauli::I {
                                        entry.remove();
                                    } else {
                                        *entry.get_mut() = factor;
                                    }
                                }
                            }
                        }
                        // HashMap iteration costs capacity, not length. Keep
                        // cancellation from leaving oversized sparse maps.
                        if left.factors.len() < left.factors.capacity() / 4 {
                            left.factors.shrink_to_fit();
                        }
                        left.phase = phase;
                        left
                    }
                };
                values.insert(node, value);
                pending.pop();
            }
            let value = &values[&root.node];
            stats.charge(value.factors.len(), max_work)?;
            stats.output_factors += value.factors.len();
            visit(root, value)?;
            let remaining = uses.get_mut(&root.node).unwrap();
            *remaining -= 1;
            if *remaining == 0 {
                values.remove(&root.node);
            }
        }
        Ok(stats)
    }
}

#[cfg(test)]
mod tests;

fn take(
    node: usize,
    uses: &mut HashMap<usize, usize>,
    values: &mut HashMap<usize, SparsePauli>,
    stats: &mut ExpansionStats,
    limit: usize,
) -> Result<SparsePauli, PbcError> {
    let remaining = uses.get_mut(&node).unwrap();
    *remaining -= 1;
    if *remaining == 0 {
        Ok(values.remove(&node).unwrap())
    } else {
        let value = &values[&node];
        stats.charge(value.factors.len(), limit)?;
        stats.copied_factors += value.factors.len();
        Ok(value.clone())
    }
}
