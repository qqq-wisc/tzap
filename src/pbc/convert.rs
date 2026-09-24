use crate::circuit::{Circuit, Gate};
use crate::pass::Pass;

use super::{
    PauliAngle, PauliAxis, PauliNode, PauliRef, PbcCircuit, PbcError, PbcOp, check_operands,
};

/// Terminal gate-circuit to logical-PBC conversion. Ordinary gate optimization
/// pipelines retain the default `Pass<Circuit>` type.
#[derive(Clone, Copy, Debug, Default)]
pub struct ToPbc;

impl Pass<Result<PbcCircuit, PbcError>> for ToPbc {
    fn name(&self) -> &str {
        "ToPbc"
    }

    fn run(&self, circuit: &Circuit) -> Result<PbcCircuit, PbcError> {
        to_pbc(circuit)
    }
}

/// Convert Clifford+T, CZ, CCX, CCZ, and terminal measurements to logical PBC.
/// Resets and gates after measurements are rejected. Arbitrary Rz must be
/// decomposed by the caller first. The output frame is retained to preserve
/// quantum outputs, including circuits with no or partial measurements.
///
/// Time and space are O(G + Q): each gate creates at most four shared Pauli
/// nodes and seven PBC operations, while updating only O(1) frame references.
/// No axis is expanded, no algebraic equality is checked, and the frame is
/// never rescanned. The initial arena contains one identity and two leaves per
/// qubit. Classical-bit counts
/// are metadata, not a bit-count-sized allocation.
///
/// ```
/// use tzap::circuit::{Circuit, Gate};
/// use tzap::pbc::to_pbc;
/// let mut input = Circuit::new(1);
/// input.apply(Gate::h(0));
/// input.apply(Gate::t(0));
/// let pbc = to_pbc(&input)?; // X rotation, followed by an H frame
/// assert_eq!(pbc.operations().len(), 1);
/// assert_eq!(pbc.expand(pbc.output_frame().x(0), 100)?.factors, vec![tzap::pbc::Pauli::Z]);
/// # Ok::<(), tzap::pbc::PbcError>(())
/// ```
pub fn to_pbc(circuit: &Circuit) -> Result<PbcCircuit, PbcError> {
    u32::try_from(circuit.num_qubits).map_err(|_| PbcError::TooManyQubits)?;
    // Validate everything up front, so conversion itself cannot fail or panic.
    let mut measuring = false;
    for (index, gate) in circuit.gates.iter().enumerate() {
        let invalid = |cause| PbcError::InvalidInput {
            index,
            cause: Box::new(cause),
        };
        validate(circuit, gate).map_err(invalid)?;
        let is_measure = matches!(gate, Gate::measure { .. });
        if measuring && !is_measure {
            return Err(invalid(PbcError::GateAfterMeasurement));
        }
        measuring |= is_measure;
    }
    let output = PbcCircuit::new(circuit.num_qubits, circuit.num_cbits);
    let mut converter = Converter { output };
    for gate in &circuit.gates {
        converter.gate(gate);
    }
    Ok(converter.output)
}

fn validate(circuit: &Circuit, gate: &Gate) -> Result<(), PbcError> {
    if matches!(gate, Gate::rz(..) | Gate::reset(_)) {
        return Err(PbcError::UnsupportedGate(gate.kind()));
    }
    check_operands(gate, circuit.num_qubits)?;
    match *gate {
        Gate::measure { cbit, .. } if cbit as usize >= circuit.num_cbits => {
            Err(PbcError::ClassicalBitOutOfRange(cbit))
        }
        _ => Ok(()),
    }
}

struct Converter {
    output: PbcCircuit,
}

impl Converter {
    fn product(&mut self, a: PauliRef, b: PauliRef) -> PauliRef {
        self.output.arena.push(PauliNode::Product(a, b))
    }

    /// Emit an unchecked rotation. Every axis is a conjugated Pauli or a
    /// product of commuting images of distinct input wires, so it is Hermitian.
    fn rotate(&mut self, axis: PauliRef, eighths: i64) {
        self.output.operations.push(PbcOp::Rotate {
            axis: PauliAxis(axis),
            angle: PauliAngle::new(eighths),
        });
    }

    /// A doubly controlled Pauli, as rotations about the products of the
    /// images of Z on both controls and of `d` on the target.
    fn controlled_controlled(&mut self, a: PauliRef, b: PauliRef, d: PauliRef) {
        let ab = self.product(a, b);
        let ad = self.product(a, d);
        let bd = self.product(b, d);
        let abd = self.product(ab, d);
        for (axis, k) in [
            (a, 1),
            (b, 1),
            (d, 1),
            (ab, -1),
            (ad, -1),
            (bd, -1),
            (abd, 1),
        ] {
            self.rotate(axis, k);
        }
    }

    fn gate(&mut self, gate: &Gate) {
        let x = |q: u32| self.output.output_frame.x(q);
        let z = |q: u32| self.output.output_frame.z(q);
        match *gate {
            Gate::t(q) => self.rotate(z(q), 1),
            Gate::tdg(q) => self.rotate(z(q), -1),
            Gate::ccx {
                control1,
                control2,
                target,
            } => self.controlled_controlled(z(control1), z(control2), x(target)),
            Gate::ccz {
                control1,
                control2,
                target,
            } => self.controlled_controlled(z(control1), z(control2), z(target)),
            Gate::measure { qubit, cbit } => {
                self.output
                    .measure(PauliAxis(z(qubit)), Some(cbit))
                    .expect("measurement operands are validated");
            }
            Gate::rz(..) | Gate::reset(_) => unreachable!("input validated before conversion"),
            _ => {
                self.output
                    .output_frame
                    .append(&mut self.output.arena, gate);
            }
        }
    }
}

#[cfg(test)]
mod tests;
