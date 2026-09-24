//! Demand-driven sparse evaluation. No allocation is indexed by total qubits
//! or arena length, and a child's value is moved into its last user, not copied.
use super::*;
use std::collections::btree_map::{BTreeMap, Entry};
use std::collections::{HashMap, HashSet};

/// Non-identity factors of a Pauli product, in increasing qubit order.
pub(crate) type Factors = BTreeMap<Qubit, Pauli>;

#[derive(Clone, Debug, Default)]
struct SparsePauli {
    phase: Phase,
    factors: Factors,
}

impl SparsePauli {
    fn single(qubit: Qubit, pauli: Pauli) -> Self {
        let mut value = Self::default();
        if pauli != Pauli::I {
            value.factors.insert(qubit, pauli);
        }
        value
    }

    /// The ordered product `self * rhs`, merging the smaller support into the
    /// larger one.
    fn mul(mut self, mut rhs: Self) -> Self {
        let swapped = self.factors.len() < rhs.factors.len();
        if swapped {
            std::mem::swap(&mut self, &mut rhs);
        }
        self.phase = self.phase * rhs.phase;
        for (q, other) in rhs.factors {
            match self.factors.entry(q) {
                Entry::Vacant(entry) => {
                    entry.insert(other);
                }
                Entry::Occupied(mut entry) => {
                    let base = *entry.get();
                    let (phase, factor) = if swapped {
                        other.times(base)
                    } else {
                        base.times(other)
                    };
                    self.phase = self.phase * phase;
                    if factor == Pauli::I {
                        entry.remove();
                    } else {
                        *entry.get_mut() = factor;
                    }
                }
            }
        }
        self
    }
}

/// Deterministic work counts, also used to test scaling without timing noise.
#[derive(Debug)]
pub(crate) struct ExpansionStats {
    limit: usize,
    pub nodes: usize,
    pub copied_factors: usize,
    pub multiplied_factors: usize,
    pub output_factors: usize,
    pub work: usize,
}

impl ExpansionStats {
    fn charge(&mut self, amount: usize) -> Result<(), PbcError> {
        self.work = self
            .work
            .checked_add(amount)
            .filter(|&work| work <= self.limit)
            .ok_or(PbcError::ExpansionLimit)?;
        Ok(())
    }
}

/// Evaluated nodes, each kept only until its last remaining use.
struct Values {
    uses: HashMap<usize, usize>,
    values: HashMap<usize, SparsePauli>,
}

impl Values {
    /// Consume one use of `reference`, moving its value out on the last use.
    fn take(
        &mut self,
        reference: PauliRef,
        stats: &mut ExpansionStats,
    ) -> Result<SparsePauli, PbcError> {
        let mut value = if self.release(reference.node) {
            self.values.remove(&reference.node).unwrap()
        } else {
            let value = &self.values[&reference.node];
            stats.charge(value.factors.len())?;
            stats.copied_factors += value.factors.len();
            value.clone()
        };
        value.phase = value.phase * reference.phase;
        Ok(value)
    }

    /// Consume one use of `node`, returning whether it was the last.
    fn release(&mut self, node: usize) -> bool {
        let remaining = self.uses.get_mut(&node).unwrap();
        *remaining -= 1;
        *remaining == 0
    }
}

impl PauliArena {
    /// Evaluate reachable nodes once and visit roots in the supplied order,
    /// passing each root's index, overall phase, and factors.
    ///
    /// Work is O(reachable nodes + roots + copied/merged/output factors), up
    /// to logarithmic map costs. Shared large intermediates can still make
    /// this superlinear in arena size; `max_work` bounds that work before
    /// each copy, merge, or output.
    pub fn materialize(
        &self,
        roots: &[PauliRef],
        max_work: usize,
        mut visit: impl FnMut(usize, Phase, &Factors) -> Result<(), PbcError>,
    ) -> Result<ExpansionStats, PbcError> {
        let mut stats = ExpansionStats {
            limit: max_work,
            nodes: 0,
            copied_factors: 0,
            multiplied_factors: 0,
            output_factors: 0,
            work: 0,
        };
        // Count graph edges once, plus each external root use. This enables
        // last-use ownership transfer even with diamonds and repeated roots.
        let mut values = Values {
            uses: HashMap::new(),
            values: HashMap::new(),
        };
        let mut pending = Vec::with_capacity(roots.len());
        for &root in roots {
            self.check(root)?;
            stats.charge(1)?;
            *values.uses.entry(root.node).or_default() += 1;
            pending.push(root.node);
        }
        let mut discovered = HashSet::new();
        while let Some(node) = pending.pop() {
            if !discovered.insert(node) {
                continue;
            }
            stats.charge(1)?;
            stats.nodes += 1;
            if let PauliNode::Product(a, b) = self.nodes[node] {
                for child in [a, b] {
                    *values.uses.entry(child.node).or_default() += 1;
                    pending.push(child.node);
                }
            }
        }
        drop(discovered);

        // Iterative post-order evaluation, so deep DAGs cannot overflow the stack.
        for (index, &root) in roots.iter().enumerate() {
            pending.push(root.node);
            while let Some(&node) = pending.last() {
                if values.values.contains_key(&node) {
                    pending.pop();
                    continue;
                }
                let value = match self.nodes[node] {
                    PauliNode::Identity => SparsePauli::default(),
                    PauliNode::Single { qubit, pauli } => SparsePauli::single(qubit, pauli),
                    PauliNode::Product(a, b) => {
                        if let Some(child) = [a, b]
                            .into_iter()
                            .find(|child| !values.values.contains_key(&child.node))
                        {
                            pending.push(child.node);
                            continue;
                        }
                        let left = values.take(a, &mut stats)?;
                        let right = values.take(b, &mut stats)?;
                        let merged = left.factors.len().min(right.factors.len());
                        stats.charge(merged)?;
                        stats.multiplied_factors += merged;
                        left.mul(right)
                    }
                };
                values.values.insert(node, value);
                pending.pop();
            }
            let value = &values.values[&root.node];
            stats.charge(value.factors.len())?;
            stats.output_factors += value.factors.len();
            visit(index, value.phase * root.phase, &value.factors)?;
            if values.release(root.node) {
                values.values.remove(&root.node);
            }
        }
        Ok(stats)
    }
}

#[cfg(test)]
mod tests;
