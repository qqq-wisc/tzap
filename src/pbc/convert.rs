use crate::circuit::{Circuit, Gate};
use crate::pass::Pass;

use super::{
    Pauli, PauliAngle, PauliAxis, PauliNode, PauliRef, PbcCircuit, PbcError, PbcOp, Phase,
    check_operands,
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
/// decomposed by the caller first. The suffix is always retained to preserve
/// quantum outputs, including circuits with no or partial measurements.
///
/// Time and space are O(G + Q): each gate creates at most four shared Pauli
/// nodes, seven PBC operations, and one suffix gate. No axis is expanded, no
/// algebraic equality is checked, and the suffix is never rescanned. The initial
/// arena contains one identity and two leaves per qubit. Classical-bit counts
/// are metadata, not a bit-count-sized allocation.
///
/// ```
/// use tzap::circuit::{Circuit, Gate};
/// use tzap::pbc::to_pbc;
/// let mut input = Circuit::new(1);
/// input.apply(Gate::h(0));
/// input.apply(Gate::t(0));
/// let pbc = to_pbc(&input)?; // X rotation, followed by H
/// assert_eq!(pbc.operations().len(), 1);
/// assert_eq!(pbc.output_cliffords(), &[Gate::h(0)]);
/// # Ok::<(), tzap::pbc::PbcError>(())
/// ```
pub fn to_pbc(circuit: &Circuit) -> Result<PbcCircuit, PbcError> {
    let qubits = u32::try_from(circuit.num_qubits).map_err(|_| PbcError::TooManyQubits)?;
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
    let mut output = PbcCircuit::new(circuit.num_qubits, circuit.num_cbits);
    let frame = (0..qubits)
        .map(|q| Frame {
            x: output.arena.push(PauliNode::Single {
                qubit: q,
                pauli: Pauli::X,
            }),
            z: output.arena.push(PauliNode::Single {
                qubit: q,
                pauli: Pauli::Z,
            }),
        })
        .collect();
    let mut converter = Converter { output, frame };
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

/// Images C† Xq C and C† Zq C for the postponed Clifford suffix C.
/// At every prefix, input = C * emitted-PBC, up to branch-global phase.
#[derive(Clone, Copy)]
struct Frame {
    x: PauliRef,
    z: PauliRef,
}

struct Converter {
    output: PbcCircuit,
    frame: Vec<Frame>,
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
        let frame = |q: u32| self.frame[q as usize];
        match *gate {
            Gate::t(q) => self.rotate(frame(q).z, 1),
            Gate::tdg(q) => self.rotate(frame(q).z, -1),
            Gate::ccx {
                control1,
                control2,
                target,
            } => self.controlled_controlled(frame(control1).z, frame(control2).z, frame(target).x),
            Gate::ccz {
                control1,
                control2,
                target,
            } => self.controlled_controlled(frame(control1).z, frame(control2).z, frame(target).z),
            Gate::measure { qubit, cbit } => {
                self.output
                    .measure(PauliAxis(frame(qubit).z), Some(cbit))
                    .expect("measurement operands are validated");
            }
            Gate::rz(..) | Gate::reset(_) => unreachable!("input validated before conversion"),
            _ => {
                self.push_clifford(gate);
                // Append once; the suffix is never simplified or rescanned.
                self.output.output_cliffords.push(gate.clone());
            }
        }
    }

    /// Append a Clifford gate G to the suffix. The suffix operator becomes
    /// G·C (G runs after C), so the image of each Pauli P becomes the image of
    /// G† P G, a product of existing images.
    fn push_clifford(&mut self, gate: &Gate) {
        match *gate {
            Gate::h(q) => {
                let f = &mut self.frame[q as usize];
                std::mem::swap(&mut f.x, &mut f.z);
            }
            Gate::x(q) => {
                let f = &mut self.frame[q as usize];
                f.z = f.z.scaled(Phase::MinusOne);
            }
            Gate::z(q) => {
                let f = &mut self.frame[q as usize];
                f.x = f.x.scaled(Phase::MinusOne);
            }
            Gate::s(q) | Gate::sdg(q) => {
                // G† X G is -Y for S and +Y for Sdg; Y = iXZ.
                let phase = if matches!(gate, Gate::s(_)) {
                    Phase::MinusI
                } else {
                    Phase::I
                };
                let f = self.frame[q as usize];
                self.frame[q as usize].x = self.product(f.x, f.z).scaled(phase);
            }
            Gate::cnot { control, target } => {
                let c = self.frame[control as usize];
                let t = self.frame[target as usize];
                self.frame[control as usize].x = self.product(c.x, t.x);
                self.frame[target as usize].z = self.product(c.z, t.z);
            }
            Gate::cz { control, target } => {
                let c = self.frame[control as usize];
                let t = self.frame[target as usize];
                self.frame[control as usize].x = self.product(c.x, t.z);
                self.frame[target as usize].x = self.product(c.z, t.x);
            }
            _ => unreachable!("not a suffix Clifford: {gate:?}"),
        }
    }
}

#[cfg(test)]
mod tests;
