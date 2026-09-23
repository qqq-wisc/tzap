//! Logical Pauli-based circuits with shared, exact Pauli expressions.
//!
//! Circuit construction does not simulate outcomes. The Clifford suffix is
//! executed after the PBC operations and preserves the quantum output frame.
//! Rotation angles use `exp(-i * k*pi/8 * P)`, unlike QASM's half-angle convention.
//!
//! ```
//! use tzap::pbc::{PauliAngle, PbcCircuit};
//! let mut pbc = PbcCircuit::new(2, 1);
//! let x = pbc.x(0)?;
//! let z = pbc.z(1)?;
//! let product = pbc.product(x.as_ref(), z.as_ref())?;
//! let axis = pbc.hermitian_axis(product, 1024)?;
//! pbc.rotate(axis, PauliAngle::new(1))?;
//! let outcome = pbc.measure(z, Some(0))?;
//! pbc.conditional_pauli(x, outcome)?;
//! println!("{}", pbc.to_ascii()?);
//! # Ok::<(), tzap::pbc::PbcError>(())
//! ```

mod ascii;
mod convert;
mod pauli;

pub use ascii::AsciiOptions;
pub use convert::{ToPbc, to_pbc};
pub use pauli::{ExpandedPauli, Pauli, PauliAxis, PauliNode, PauliRef, Phase};

use crate::circuit::{CBit, Gate, GateKind, Qubit, qubit_operands};
use pauli::PauliArena;
use std::fmt;

/// Exact multiples of pi/8 modulo pi, ignoring global phase.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PauliAngle(u8);

impl PauliAngle {
    pub fn new(eighths: i64) -> Self {
        Self(eighths.rem_euclid(8) as u8)
    }
    pub fn eighths(self) -> u8 {
        self.0
    }
    pub fn inverse(self) -> Self {
        Self::new(-i64::from(self.0))
    }
    pub fn plus(self, other: Self) -> Self {
        Self((self.0 + other.0) % 8)
    }
}

impl fmt::Display for PauliAngle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(
            [
                "0", "pi/8", "pi/4", "3pi/8", "pi/2", "-3pi/8", "-pi/4", "-pi/8",
            ][self.0 as usize],
        )
    }
}

/// Immutable measurement identity, independent of the mutable user classical bits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MeasId {
    owner: u64,
    index: usize,
}

impl MeasId {
    pub fn index(self) -> usize {
        self.index
    }
}
impl fmt::Display for MeasId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "m{}", self.index)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PbcOp {
    Rotate {
        axis: PauliAxis,
        angle: PauliAngle,
    },
    /// +1 eigenvalue yields bit 0; -1 yields bit 1. None keeps the outcome
    /// internal, as required by reset lowering.
    Measure {
        axis: PauliAxis,
        outcome: MeasId,
        target: Option<CBit>,
    },
    /// Apply exp(-i * angle * axis) when the recorded outcome is 1;
    /// otherwise act as identity. The condition is classical, not coherent.
    ConditionalRotate {
        axis: PauliAxis,
        angle: PauliAngle,
        if_one: MeasId,
    },
}

impl PbcOp {
    pub fn axis(&self) -> PauliAxis {
        match *self {
            Self::Rotate { axis, .. }
            | Self::Measure { axis, .. }
            | Self::ConditionalRotate { axis, .. } => axis,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PbcError {
    /// Invalid source instruction; index is zero-based.
    InvalidInput {
        index: usize,
        cause: Box<PbcError>,
    },
    UnsupportedGate(GateKind),
    TooManyQubits,
    QubitOutOfRange(Qubit),
    ClassicalBitOutOfRange(CBit),
    ForeignPauli,
    UnknownMeasurement,
    NonHermitianAxis,
    NonCliffordSuffix,
    RepeatedOperand,
    ExpansionLimit,
    DrawingLimit,
}

impl fmt::Display for PbcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { index, cause } => write!(f, "input gate {index}: {cause}"),
            Self::UnsupportedGate(gate) => write!(f, "unsupported PBC input gate: {gate:?}"),
            Self::TooManyQubits => f.write_str("PBC qubit count exceeds the supported u32 range"),
            Self::QubitOutOfRange(q) => write!(f, "PBC qubit {q} is out of range"),
            Self::ClassicalBitOutOfRange(c) => write!(f, "PBC classical bit {c} is out of range"),
            Self::ForeignPauli => f.write_str("Pauli reference belongs to another circuit"),
            Self::UnknownMeasurement => f.write_str("measurement is not defined in this circuit"),
            Self::NonHermitianAxis => {
                f.write_str("rotation and measurement axes must be Hermitian")
            }
            Self::NonCliffordSuffix => f.write_str("output suffix only accepts Clifford gates"),
            Self::RepeatedOperand => f.write_str("gate operands must be distinct"),
            Self::ExpansionLimit => f.write_str("Pauli expansion exceeds the cell budget"),
            Self::DrawingLimit => f.write_str("circuit exceeds the ASCII drawing limits"),
        }
    }
}
impl std::error::Error for PbcError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidInput { cause, .. } => Some(cause.as_ref()),
            _ => None,
        }
    }
}

/// An owned PBC program. Handles are scoped to this circuit. Read-only slices
/// expose the IR; checked methods maintain operand and outcome validity.
///
/// No Clone implementation is provided: duplicating an arena requires remapping
/// its scoped handles, otherwise separately edited clones could alias node IDs.
#[derive(Debug)]
pub struct PbcCircuit {
    num_qubits: usize,
    num_cbits: usize,
    arena: PauliArena,
    operations: Vec<PbcOp>,
    output_cliffords: Vec<Gate>,
    measurements: usize,
}

impl PbcCircuit {
    pub fn new(num_qubits: usize, num_cbits: usize) -> Self {
        Self {
            num_qubits,
            num_cbits,
            arena: PauliArena::new(),
            operations: Vec::new(),
            output_cliffords: Vec::new(),
            measurements: 0,
        }
    }
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }
    pub fn num_cbits(&self) -> usize {
        self.num_cbits
    }
    pub fn operations(&self) -> &[PbcOp] {
        &self.operations
    }
    pub fn pauli_nodes(&self) -> &[PauliNode] {
        &self.arena.nodes
    }
    pub fn output_cliffords(&self) -> &[Gate] {
        &self.output_cliffords
    }
    pub fn measurement_count(&self) -> usize {
        self.measurements
    }
    pub fn identity(&self) -> PauliAxis {
        PauliAxis(self.arena.reference(0))
    }

    fn check_qubit(&self, q: Qubit) -> Result<(), PbcError> {
        if q as usize >= self.num_qubits {
            return Err(PbcError::QubitOutOfRange(q));
        }
        Ok(())
    }

    pub fn single(&mut self, qubit: Qubit, pauli: Pauli) -> Result<PauliAxis, PbcError> {
        self.check_qubit(qubit)?;
        if pauli == Pauli::I {
            return Ok(self.identity());
        }
        Ok(PauliAxis(
            self.arena.push(PauliNode::Single { qubit, pauli }),
        ))
    }
    pub fn x(&mut self, q: Qubit) -> Result<PauliAxis, PbcError> {
        self.single(q, Pauli::X)
    }
    pub fn y(&mut self, q: Qubit) -> Result<PauliAxis, PbcError> {
        self.single(q, Pauli::Y)
    }
    pub fn z(&mut self, q: Qubit) -> Result<PauliAxis, PbcError> {
        self.single(q, Pauli::Z)
    }

    /// Append one ordered product without expanding either operand. Amortized O(1).
    pub fn product(&mut self, left: PauliRef, right: PauliRef) -> Result<PauliRef, PbcError> {
        self.arena.check(left)?;
        self.arena.check(right)?;
        Ok(self.arena.push(PauliNode::Product(left, right)))
    }

    /// Explicitly inspect an expression. Cost and storage are O(D * Q) for the
    /// arena prefix of D nodes; max_cells bounds D * max(Q, 1).
    pub fn expand(&self, reference: PauliRef, max_cells: usize) -> Result<ExpandedPauli, PbcError> {
        self.arena.check(reference)?;
        let mut values = self
            .arena
            .expand(self.num_qubits, reference.node, max_cells)?;
        let mut result = values.pop().expect("every reference has a node");
        result.phase = result.phase.times(reference.phase);
        Ok(result)
    }

    /// Check an arbitrary expression before using it as a physical axis.
    /// This explicit materialization is intended for manual construction. A
    /// converter constructs axes from its Clifford-frame invariants without
    /// expanding them.
    pub fn hermitian_axis(
        &self,
        reference: PauliRef,
        max_cells: usize,
    ) -> Result<PauliAxis, PbcError> {
        let phase = self.expand(reference, max_cells)?.phase;
        if !matches!(phase, Phase::One | Phase::MinusOne) {
            return Err(PbcError::NonHermitianAxis);
        }
        Ok(PauliAxis(reference))
    }

    pub fn rotate(&mut self, axis: PauliAxis, angle: PauliAngle) -> Result<(), PbcError> {
        self.arena.check(axis.0)?;
        self.operations.push(PbcOp::Rotate { axis, angle });
        Ok(())
    }

    pub fn measure(&mut self, axis: PauliAxis, target: Option<CBit>) -> Result<MeasId, PbcError> {
        self.arena.check(axis.0)?;
        if let Some(c) = target
            && c as usize >= self.num_cbits
        {
            return Err(PbcError::ClassicalBitOutOfRange(c));
        }
        let outcome = MeasId {
            owner: self.arena.owner,
            index: self.measurements,
        };
        self.operations.push(PbcOp::Measure {
            axis,
            outcome,
            target,
        });
        self.measurements += 1;
        Ok(outcome)
    }

    /// Append a classically conditioned Pauli-string rotation. Angles are
    /// multiples of pi/8, modulo an unobservable phase within the outcome branch.
    pub fn conditional_rotate(
        &mut self,
        axis: PauliAxis,
        angle: PauliAngle,
        if_one: MeasId,
    ) -> Result<(), PbcError> {
        self.arena.check(axis.0)?;
        if if_one.owner != self.arena.owner || if_one.index >= self.measurements {
            return Err(PbcError::UnknownMeasurement);
        }
        self.operations.push(PbcOp::ConditionalRotate {
            axis,
            angle,
            if_one,
        });
        Ok(())
    }

    /// Apply a Pauli string conditionally via a pi/2 rotation. This gives -iP,
    /// equivalent to P up to an unobservable phase in the classical branch.
    pub fn conditional_pauli(&mut self, axis: PauliAxis, if_one: MeasId) -> Result<(), PbcError> {
        self.conditional_rotate(axis, PauliAngle::new(4), if_one)
    }

    /// Append a gate to the suffix executed after all PBC operations, regardless
    /// of the order in which the builder's methods are called.
    pub fn push_output_clifford(&mut self, gate: Gate) -> Result<(), PbcError> {
        if !matches!(
            gate,
            Gate::h(_)
                | Gate::x(_)
                | Gate::z(_)
                | Gate::s(_)
                | Gate::sdg(_)
                | Gate::cnot { .. }
                | Gate::cz { .. }
        ) {
            return Err(PbcError::NonCliffordSuffix);
        }
        let (n, qs) = qubit_operands(&gate);
        for (i, &q) in qs[..n].iter().enumerate() {
            self.check_qubit(q)?;
            if qs[..i].contains(&q) {
                return Err(PbcError::RepeatedOperand);
            }
        }
        self.output_cliffords.push(gate);
        Ok(())
    }
}

#[cfg(test)]
mod tests;
