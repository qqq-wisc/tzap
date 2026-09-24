//! Exact quantum-classical channels for a fixed, explicit initial classical
//! store and arbitrary quantum input (including inputs entangled with a reference).
//!
//! Branches remain unnormalized: rho -> K rho K†. Measurement bit 0 is the +1
//! eigenvalue. All output qubits and final user classical bits are observable;
//! internal outcome IDs are not. Hidden histories are summed incoherently into
//! Choi blocks, never compared as individual Kraus operators or sampled.
//! Gate inputs require terminal measurements. PBC may include rotations after
//! measurements to restore quantum outputs. Resets and feed-forward are rejected.

use super::*;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug)]
pub(crate) struct ChannelLimits {
    pub matrices: Limits,
    pub max_branches: usize,
    pub max_classical_bits: usize,
}

impl Default for ChannelLimits {
    fn default() -> Self {
        Self {
            matrices: Limits::default(),
            max_branches: 64,
            max_classical_bits: 64,
        }
    }
}

/// Row-major vec(K): index = output_basis * quantum_dim + input_basis.
/// Each block is sum_h vec(K_h) vec(K_h)† for one final classical store.
/// Zero blocks are omitted canonically. Classical bits are in c0,c1,... order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Channel {
    pub quantum_dim: usize,
    pub num_cbits: usize,
    pub blocks: BTreeMap<Vec<bool>, Matrix>,
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum Difference {
    Dimensions,
    Entry {
        classical: Vec<bool>,
        row: usize,
        col: usize,
        values: Box<(Scalar, Scalar)>,
    },
}

impl Channel {
    /// Trace out quantum outputs from each Choi block, retaining the classical
    /// channel for arbitrary input states (not just computational-basis inputs).
    /// With our row-major vec(K), these matrices are transposed POVM effects.
    pub fn classical_blocks(&self) -> BTreeMap<Vec<bool>, Matrix> {
        self.blocks
            .iter()
            .map(|(store, block)| {
                let dim = self.quantum_dim;
                let mut effect = Matrix::zero(dim);
                for row in 0..dim {
                    for col in 0..dim {
                        let mut value = Scalar::zero();
                        for output in 0..dim {
                            value = value.add(block.get(output * dim + row, output * dim + col));
                        }
                        effect.set(row, col, value);
                    }
                }
                (store.clone(), effect)
            })
            .collect()
    }

    /// Compare exact entries, treating a missing block as zero. Returns the
    /// first differing classical store and matrix position, not a huge dump.
    pub fn compare(&self, other: &Self) -> Result<(), Difference> {
        if self.quantum_dim != other.quantum_dim || self.num_cbits != other.num_cbits {
            return Err(Difference::Dimensions);
        }
        let keys: BTreeSet<_> = self.blocks.keys().chain(other.blocks.keys()).collect();
        let zero = Scalar::zero();
        let width = self.quantum_dim * self.quantum_dim;
        for key in keys {
            let a = self.blocks.get(key);
            let b = other.blocks.get(key);
            for row in 0..width {
                for col in 0..width {
                    let x = a.map_or(&zero, |m| m.get(row, col));
                    let y = b.map_or(&zero, |m| m.get(row, col));
                    if x != y {
                        return Err(Difference::Entry {
                            classical: key.clone(),
                            row,
                            col,
                            values: Box::new((x.clone(), y.clone())),
                        });
                    }
                }
            }
        }
        Ok(())
    }
}

struct Branch {
    kraus: Matrix,
    classical: Vec<bool>,
}

struct Footprint {
    operations: usize,
    nodes: usize,
    splits: usize,
    destinations: usize,
}

/// Conservative preflight: account for branch splitting, retained DAG matrices,
/// dense Choi blocks, and scalar-product work before allocating any matrices.
fn preflight(
    n: usize,
    cbits: usize,
    initial: &[bool],
    footprint: Footprint,
    limits: ChannelLimits,
) -> Result<usize, Error> {
    let Footprint {
        operations,
        nodes,
        splits,
        destinations,
    } = footprint;
    if cbits != initial.len() {
        return Err(Error::InvalidInitialStore);
    }
    if cbits > limits.max_classical_bits {
        return Err(Error::LimitExceeded);
    }
    let fail = Error::LimitExceeded;
    let branches = 1usize
        .checked_shl(u32::try_from(splits).map_err(|_| fail)?)
        .ok_or(fail)?;
    if branches > limits.max_branches {
        return Err(fail);
    }
    let dim = dimensions(n, SCRATCH_MATRICES, 0, operations, limits.matrices)?;
    let square = dim.checked_mul(dim).ok_or(fail)?;
    let fourth = square.checked_mul(square).ok_or(fail)?;
    // At most 2^distinct-written-bits observable stores. Hidden measurements
    // do not create additional visible blocks.
    let visible = 1usize
        .checked_shl(u32::try_from(destinations.min(splits)).map_err(|_| fail)?)
        .ok_or(fail)?;
    let retained = branches
        .checked_mul(3)
        .and_then(|x| x.checked_add(nodes))
        .and_then(|x| x.checked_add(SCRATCH_MATRICES))
        .and_then(|x| x.checked_mul(square))
        .ok_or(fail)?;
    let cells = visible
        .checked_mul(fourth)
        .and_then(|x| x.checked_add(retained))
        .ok_or(fail)?;
    let multiplications = operations
        .checked_mul(branches)
        .and_then(|x| x.checked_mul(2))
        .and_then(|x| x.checked_add(nodes))
        .ok_or(fail)?;
    let terms = multiplications
        .checked_mul(square)
        .and_then(|x| x.checked_mul(dim))
        .and_then(|x| branches.checked_mul(fourth).and_then(|y| x.checked_add(y)))
        .ok_or(fail)?;
    if cells > limits.matrices.max_matrix_cells || terms > limits.matrices.max_multiply_terms {
        return Err(fail);
    }
    Ok(dim)
}

fn initial_branch(dim: usize, classical: &[bool]) -> Vec<Branch> {
    vec![Branch {
        kraus: Matrix::identity(dim),
        classical: classical.to_vec(),
    }]
}

fn apply(branches: &mut [Branch], unitary: &Matrix) {
    for branch in branches {
        branch.kraus = unitary.mul(&branch.kraus);
    }
}

fn split(branches: Vec<Branch>, operators: [Matrix; 2], target: Option<u32>) -> Vec<Branch> {
    let mut result = Vec::with_capacity(branches.len() * 2);
    for branch in branches {
        for (outcome, operator) in operators.iter().enumerate() {
            let mut next = Branch {
                kraus: operator.mul(&branch.kraus),
                classical: branch.classical.clone(),
            };
            let bit = outcome == 1;
            if let Some(cbit) = target {
                next.classical[cbit as usize] = bit;
            }
            result.push(next);
        }
    }
    result
}

fn finish(dim: usize, cbits: usize, branches: Vec<Branch>) -> Channel {
    let mut blocks = BTreeMap::<Vec<bool>, Matrix>::new();
    let zero = Scalar::zero();
    for branch in branches {
        // Impossible outcomes carry zero probability and must not introduce
        // spurious observable blocks (e.g. the -1 outcome of measuring +I).
        let entries = &branch.kraus.entries;
        if entries.iter().all(|x| x == &zero) {
            continue;
        }
        let block = blocks
            .entry(branch.classical)
            .or_insert_with(|| Matrix::zero(dim * dim));
        for (row, a) in entries.iter().enumerate() {
            if a == &zero {
                continue;
            }
            for (col, b) in entries.iter().enumerate() {
                if b == &zero {
                    continue;
                }
                block.set(row, col, block.get(row, col).add(&a.mul(&b.conj())));
            }
        }
    }
    Channel {
        quantum_dim: dim,
        num_cbits: cbits,
        blocks,
    }
}

/// Gate semantics from native matrices and direct terminal-measurement Kraus
/// operators. `initial` contains c0,c1,...; arbitrary Rz remains unsupported.
pub(crate) fn circuit_channel(
    circuit: &Circuit,
    initial: &[bool],
    limits: ChannelLimits,
) -> Result<Channel, Error> {
    if circuit.gates.len() > limits.matrices.max_operations {
        return Err(Error::LimitExceeded);
    }
    let mut splits = 0;
    let mut destinations = BTreeSet::new();
    for (index, gate) in circuit.gates.iter().enumerate() {
        let (n, qs) = qubit_operands(gate);
        for (i, &q) in qs[..n].iter().enumerate() {
            if q as usize >= circuit.num_qubits || qs[..i].contains(&q) {
                return Err(Error::InvalidOperand { index });
            }
        }
        match gate {
            Gate::measure { cbit, .. } => {
                if *cbit as usize >= circuit.num_cbits {
                    return Err(Error::InvalidOperand { index });
                }
                splits += 1;
                destinations.insert(*cbit);
            }
            Gate::reset(_) | Gate::rz(..) => return Err(Error::UnsupportedOperation { index }),
            _ if splits > 0 => return Err(Error::GateAfterMeasurement { index }),
            _ => (),
        }
    }
    let dim = preflight(
        circuit.num_qubits,
        circuit.num_cbits,
        initial,
        Footprint {
            operations: circuit.gates.len(),
            nodes: 0,
            splits,
            destinations: destinations.len(),
        },
        limits,
    )?;
    let mut branches = initial_branch(dim, initial);
    for (index, gate) in circuit.gates.iter().enumerate() {
        match *gate {
            Gate::measure { qubit, cbit } => {
                // Direct computational-basis projectors, independent of the
                // PBC interpreter's Pauli-projector formula.
                let mut operators = [Matrix::zero(dim), Matrix::zero(dim)];
                let mask = 1 << (circuit.num_qubits - 1 - qubit as usize);
                for col in 0..dim {
                    let bit = usize::from(col & mask != 0);
                    operators[bit].set(col, col, Scalar::integer(1));
                }
                branches = split(branches, operators, Some(cbit));
            }
            _ => apply(
                &mut branches,
                &gate_matrix(circuit.num_qubits, gate, index)?,
            ),
        }
    }
    Ok(finish(dim, circuit.num_cbits, branches))
}

/// PBC semantics for rotations and measurements, followed by any internal suffix.
/// Preserves all quantum outputs; only internal measurement histories are hidden.
pub(crate) fn pbc_channel(
    circuit: &PbcCircuit,
    initial: &[bool],
    limits: ChannelLimits,
) -> Result<Channel, Error> {
    let operations = circuit
        .operations()
        .len()
        .checked_add(circuit.output_cliffords().len())
        .ok_or(Error::LimitExceeded)?;
    if operations > limits.matrices.max_operations {
        return Err(Error::LimitExceeded);
    }
    for (index, op) in circuit.operations().iter().enumerate() {
        if matches!(op, PbcOp::ConditionalRotate { .. }) {
            return Err(Error::UnsupportedOperation { index });
        }
    }
    let destinations: BTreeSet<_> = circuit
        .operations()
        .iter()
        .filter_map(|op| {
            if let PbcOp::Measure { target, .. } = op {
                *target
            } else {
                None
            }
        })
        .collect();
    let nodes = referenced_nodes(circuit);
    let dim = preflight(
        circuit.num_qubits(),
        circuit.num_cbits(),
        initial,
        Footprint {
            operations,
            nodes,
            splits: circuit.measurement_count(),
            destinations: destinations.len(),
        },
        limits,
    )?;
    let values = pauli_matrices(circuit, nodes, dim);
    let mut branches = initial_branch(dim, initial);
    for op in circuit.operations() {
        let axis = reference(&values, op.axis().as_ref());
        match *op {
            PbcOp::Rotate { angle, .. } => apply(&mut branches, &rotation(&axis, angle.eighths())),
            PbcOp::Measure { target, .. } => {
                let identity = Matrix::identity(dim);
                let half = Scalar::integer(1).half();
                let positive = identity.add(&axis).scale(&half);
                let negative = identity.add(&axis.scale(&Scalar::integer(-1))).scale(&half);
                branches = split(branches, [positive, negative], target);
            }
            PbcOp::ConditionalRotate { .. } => unreachable!("validated above"),
        }
    }
    for (i, gate) in circuit.output_cliffords().iter().enumerate() {
        apply(
            &mut branches,
            &gate_matrix(circuit.num_qubits(), gate, circuit.operations().len() + i)?,
        );
    }
    Ok(finish(dim, circuit.num_cbits(), branches))
}

#[cfg(test)]
mod tests;
