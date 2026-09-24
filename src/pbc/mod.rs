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
mod text;

pub use ascii::AsciiOptions;
pub use convert::{ToPbc, to_pbc};
pub use pauli::{ExpandedPauli, Pauli, PauliAxis, PauliNode, PauliRef, Phase};
pub use text::TextOptions;

use crate::circuit::{CBit, Gate, GateKind, Qubit, qubit_operands};
use pauli::{Factors, PauliArena};
use std::fmt;
use std::ops::{Add, Neg};

/// Exact multiples of pi/8 modulo pi, ignoring global phase.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PauliAngle(u8);

impl PauliAngle {
    pub fn new(eighths: i64) -> Self {
        Self(eighths.rem_euclid(8) as u8)
    }
    /// The angle as `k*pi/8` with `0 <= k < 8`.
    pub fn eighths(self) -> u8 {
        self.0
    }
    /// The angle as `k*pi/8` with `-3 <= k <= 4`.
    pub fn signed_eighths(self) -> i8 {
        let k = self.0 as i8;
        if k > 4 { k - 8 } else { k }
    }
}

impl Add for PauliAngle {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self((self.0 + rhs.0) % 8)
    }
}

impl Neg for PauliAngle {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-i64::from(self.0))
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
    /// internal, without writing a user classical bit.
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
    GateAfterMeasurement,
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
    UnsupportedTextOperation {
        index: usize,
    },
}

impl fmt::Display for PbcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { index, cause } => write!(f, "input gate {index}: {cause}"),
            Self::UnsupportedGate(gate) => write!(f, "unsupported PBC input gate: {gate:?}"),
            Self::GateAfterMeasurement => {
                f.write_str("measurements must form the final input block")
            }
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
            Self::UnsupportedTextOperation { index } => write!(
                f,
                "PBC operation {index} cannot be exported: expected rotations or measurements into classical registers"
            ),
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

    pub fn single(&mut self, qubit: Qubit, pauli: Pauli) -> Result<PauliAxis, PbcError> {
        if qubit as usize >= self.num_qubits {
            return Err(PbcError::QubitOutOfRange(qubit));
        }
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

    /// Inspect one expression using sparse reachable dependencies, then allocate
    /// the requested dense Q-qubit result. The budget charges Q output cells
    /// plus sparse evaluation work; unrelated arena nodes are never expanded.
    pub fn expand(&self, reference: PauliRef, max_cells: usize) -> Result<ExpandedPauli, PbcError> {
        self.arena.check(reference)?;
        let budget = max_cells
            .checked_sub(self.num_qubits)
            .ok_or(PbcError::ExpansionLimit)?;
        let mut result = ExpandedPauli {
            phase: Phase::One,
            factors: Vec::new(),
        };
        self.arena
            .materialize(&[reference], budget, |_, phase, factors| {
                result.phase = phase;
                result.factors = vec![Pauli::I; self.num_qubits];
                for (&q, &p) in factors {
                    result.factors[q as usize] = p;
                }
                Ok(())
            })?;
        Ok(result)
    }

    /// Check an arbitrary expression before using it as a physical axis.
    /// This explicit materialization is intended for manual construction;
    /// `to_pbc` constructs axes that are Hermitian by construction.
    pub fn hermitian_axis(
        &self,
        reference: PauliRef,
        max_cells: usize,
    ) -> Result<PauliAxis, PbcError> {
        self.arena
            .materialize(&[reference], max_cells, |_, phase, _| match phase {
                Phase::One | Phase::MinusOne => Ok(()),
                Phase::I | Phase::MinusI => Err(PbcError::NonHermitianAxis),
            })?;
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
        if !is_suffix_gate(&gate) {
            return Err(PbcError::NonCliffordSuffix);
        }
        check_operands(&gate, self.num_qubits)?;
        self.output_cliffords.push(gate);
        Ok(())
    }

    /// Materialize every operation's axis within `max_work`, visiting the
    /// operations in order with the axis's overall phase and sparse factors.
    fn visit_axes(
        &self,
        max_work: usize,
        mut visit: impl FnMut(&PbcOp, Phase, &Factors) -> Result<(), PbcError>,
    ) -> Result<(), PbcError> {
        let roots: Vec<_> = self.operations.iter().map(|op| op.axis().0).collect();
        self.arena
            .materialize(&roots, max_work, |index, phase, factors| {
                visit(&self.operations[index], phase, factors)
            })?;
        Ok(())
    }
}

/// The Clifford gates allowed in the output suffix.
fn is_suffix_gate(gate: &Gate) -> bool {
    use GateKind::*;
    matches!(gate.kind(), H | X | Z | S | Sdg | Cx | Cz)
}

/// Check that a gate's qubit operands are in range and distinct.
fn check_operands(gate: &Gate, num_qubits: usize) -> Result<(), PbcError> {
    let (n, qs) = qubit_operands(gate);
    for (i, &q) in qs[..n].iter().enumerate() {
        if q as usize >= num_qubits {
            return Err(PbcError::QubitOutOfRange(q));
        }
        if qs[..i].contains(&q) {
            return Err(PbcError::RepeatedOperand);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
