//! Exact, bounded reference semantics for small quantum circuits.
//!
//! Test-only: this is an independent oracle, not part of optimization. PBC
//! rotations use V(P,k) = (I+P)/2 + omega^k (I-P)/2, omega = exp(i*pi/4).
//! Thus the returned matrix is a unitary representative up to global phase,
//! not necessarily the literal exp(-i*k*pi*P/8). All arithmetic is exact in
//! Q(sqrt(2), i). Qubit 0 is the most significant basis bit.
//! [`channel`] extends this to exact quantum-classical maps with terminal
//! measurements, for an explicit initial classical store and arbitrary quantum
//! input. Gate inputs reject resets and mid-circuit measurements; PBC outputs
//! may have rotations after measurement to restore their quantum output frame.

pub(crate) mod channel;
mod matrix;
mod scalar;
pub(crate) mod test_support;
#[cfg(test)]
mod tests;

use crate::circuit::{Circuit, Gate, qubit_operands};
use crate::pbc::{Pauli, PauliNode, PauliRef, PbcCircuit, PbcOp, Phase};
pub(crate) use matrix::Matrix;
use scalar::Scalar;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Limits {
    pub max_qubits: usize,
    /// Conservative bound on simultaneously retained dense scalar entries.
    pub max_matrix_cells: usize,
    /// Bound on dense multiplication terms, including intermediate products.
    pub max_multiply_terms: usize,
    pub max_operations: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_qubits: 5,
            max_matrix_cells: 262_144,
            max_multiply_terms: 10_000_000,
            max_operations: 1024,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Error {
    LimitExceeded,
    UnsupportedOperation { index: usize },
    InvalidOperand { index: usize },
    InvalidInitialStore,
    GateAfterMeasurement { index: usize },
}

// Conservatively covers the accumulator, identity, gate, and intermediate
// matrices in either interpreter, in addition to retained Pauli DAG values.
const SCRATCH_MATRICES: usize = 8;

fn dimensions(
    n: usize,
    matrices: usize,
    multiplications: usize,
    operations: usize,
    limits: Limits,
) -> Result<usize, Error> {
    let fail = Error::LimitExceeded;
    if n > limits.max_qubits || operations > limits.max_operations {
        return Err(fail);
    }
    let dim = 1usize
        .checked_shl(u32::try_from(n).map_err(|_| fail)?)
        .ok_or(fail)?;
    let square = dim.checked_mul(dim).ok_or(fail)?;
    let cells = square.checked_mul(matrices).ok_or(fail)?;
    let terms = square
        .checked_mul(dim)
        .and_then(|x| x.checked_mul(multiplications))
        .ok_or(fail)?;
    if cells > limits.max_matrix_cells || terms > limits.max_multiply_terms {
        return Err(fail);
    }
    Ok(dim)
}

fn phase(p: Phase) -> Scalar {
    match p {
        Phase::One => Scalar::integer(1),
        Phase::MinusOne => Scalar::integer(-1),
        Phase::I => Scalar::i(),
        Phase::MinusI => Scalar::i().neg(),
    }
}

fn single(n: usize, q: u32, entries: [[Scalar; 2]; 2]) -> Matrix {
    let dim = 1 << n;
    let mask = 1 << (n - 1 - q as usize);
    let mut result = Matrix::zero(dim);
    for col in 0..dim {
        let bit = usize::from(col & mask != 0);
        result.set(col & !mask, col, entries[0][bit].clone());
        result.set(col | mask, col, entries[1][bit].clone());
    }
    result
}

fn pauli(n: usize, q: u32, p: Pauli) -> Matrix {
    let zero = Scalar::zero();
    let one = Scalar::integer(1);
    match p {
        Pauli::I => Matrix::identity(1 << n),
        Pauli::X => single(n, q, [[zero.clone(), one.clone()], [one, zero]]),
        Pauli::Y => single(
            n,
            q,
            [[zero.clone(), Scalar::i().neg()], [Scalar::i(), zero]],
        ),
        Pauli::Z => single(n, q, [[one, zero.clone()], [zero, Scalar::integer(-1)]]),
    }
}

/// Direct truth-table semantics for controlled X/Z, with validated operands.
/// Control bits must all be 1; X flips the target and Z negates its |1> state.
fn controlled(n: usize, controls: &[u32], target: u32, flip: bool) -> Matrix {
    let controls = controls
        .iter()
        .fold(0, |mask, &q| mask | (1 << (n - 1 - q as usize)));
    let target = 1 << (n - 1 - target as usize);
    let mut matrix = Matrix::zero(1 << n);
    for col in 0..matrix.dim {
        let active = col & controls == controls;
        let row = if flip && active { col ^ target } else { col };
        let sign = if !flip && active && col & target != 0 {
            -1
        } else {
            1
        };
        matrix.set(row, col, Scalar::integer(sign));
    }
    matrix
}

/// Reject out-of-range or repeated qubit operands of instruction `index`.
fn check_operands(n: usize, gate: &Gate, index: usize) -> Result<(), Error> {
    let (count, qs) = qubit_operands(gate);
    for (i, q) in qs[..count].iter().enumerate() {
        if *q as usize >= n || qs[..i].contains(q) {
            return Err(Error::InvalidOperand { index });
        }
    }
    Ok(())
}

fn gate_matrix(n: usize, gate: &Gate, index: usize) -> Result<Matrix, Error> {
    check_operands(n, gate, index)?;
    let zero = Scalar::zero();
    let one = Scalar::integer(1);
    Ok(match *gate {
        Gate::x(q) => pauli(n, q, Pauli::X),
        Gate::z(q) => pauli(n, q, Pauli::Z),
        Gate::h(q) => {
            let a = Scalar::inv_sqrt_two();
            single(n, q, [[a.clone(), a.clone()], [a.clone(), a.neg()]])
        }
        Gate::s(q) | Gate::sdg(q) | Gate::t(q) | Gate::tdg(q) => {
            let k = match gate {
                Gate::s(_) => 2,
                Gate::sdg(_) => 6,
                Gate::t(_) => 1,
                _ => 7,
            };
            single(n, q, [[one, zero.clone()], [zero, Scalar::omega(k)]])
        }
        Gate::cnot { control, target } | Gate::cz { control, target } => {
            controlled(n, &[control], target, matches!(gate, Gate::cnot { .. }))
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
        } => controlled(
            n,
            &[control1, control2],
            target,
            matches!(gate, Gate::ccx { .. }),
        ),
        Gate::rz(..) | Gate::measure { .. } | Gate::reset(_) => {
            return Err(Error::UnsupportedOperation { index });
        }
    })
}

pub(crate) fn circuit_unitary(circuit: &Circuit, limits: Limits) -> Result<Matrix, Error> {
    let dim = dimensions(
        circuit.num_qubits,
        SCRATCH_MATRICES,
        circuit.gates.len(),
        circuit.gates.len(),
        limits,
    )?;
    let mut result = Matrix::identity(dim);
    for (i, gate) in circuit.gates.iter().enumerate() {
        result = gate_matrix(circuit.num_qubits, gate, i)?.mul(&result);
    }
    Ok(result)
}

fn reference(values: &[Matrix], r: PauliRef) -> Matrix {
    values[r.node_index()].scale(&phase(r.phase()))
}

fn referenced_nodes(circuit: &PbcCircuit) -> usize {
    circuit
        .operations()
        .iter()
        .map(|op| op.axis().as_ref().node_index() + 1)
        .max()
        .unwrap_or(0)
}

/// Independent dense interpretation, not the production Pauli materializer.
fn pauli_matrices(circuit: &PbcCircuit, nodes: usize, dim: usize) -> Vec<Matrix> {
    let mut values = Vec::with_capacity(nodes);
    for node in &circuit.pauli_nodes()[..nodes] {
        values.push(match *node {
            PauliNode::Identity => Matrix::identity(dim),
            PauliNode::Single { qubit, pauli: p } => pauli(circuit.num_qubits(), qubit, p),
            PauliNode::Product(a, b) => reference(&values, a).mul(&reference(&values, b)),
        });
    }
    values
}

fn rotation(axis: &Matrix, eighths: u8) -> Matrix {
    let omega = Scalar::omega(eighths);
    let a = Scalar::integer(1).add(&omega).half();
    let b = Scalar::integer(1).add(&omega.neg()).half();
    Matrix::identity(axis.dim).scale(&a).add(&axis.scale(&b))
}

pub(crate) fn pbc_unitary(circuit: &PbcCircuit, limits: Limits) -> Result<Matrix, Error> {
    // Reject dynamic operations rather than silently ignoring a measurement or
    // guessing the classical branch of a conditional rotation.
    for (i, op) in circuit.operations().iter().enumerate() {
        if !matches!(op, PbcOp::Rotate { .. }) {
            return Err(Error::UnsupportedOperation { index: i });
        }
    }
    let nodes = referenced_nodes(circuit);
    let operations = circuit
        .operations()
        .len()
        .checked_add(circuit.output_cliffords().len())
        .ok_or(Error::LimitExceeded)?;
    let dim = dimensions(
        circuit.num_qubits(),
        nodes
            .checked_add(SCRATCH_MATRICES)
            .ok_or(Error::LimitExceeded)?,
        nodes.checked_add(operations).ok_or(Error::LimitExceeded)?,
        operations,
        limits,
    )?;
    let values = pauli_matrices(circuit, nodes, dim);
    let mut result = Matrix::identity(dim);
    for op in circuit.operations() {
        let PbcOp::Rotate { axis, angle } = op else {
            unreachable!("validated above")
        };
        let rotation = rotation(&reference(&values, axis.as_ref()), angle.eighths());
        result = rotation.mul(&result);
    }
    for (i, gate) in circuit.output_cliffords().iter().enumerate() {
        result =
            gate_matrix(circuit.num_qubits(), gate, circuit.operations().len() + i)?.mul(&result);
    }
    Ok(result)
}
