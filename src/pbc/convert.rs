use crate::circuit::{Circuit, Gate, GateKind, qubit_operands};
use crate::pass::Pass;

use super::{PauliAngle, PauliAxis, PauliRef, PbcCircuit, PbcError, PbcOp, Phase};

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
    // Validate before building output so invalid public Circuit data cannot
    // cause indexing panics or silently lose unsupported operations.
    let mut measuring = false;
    for (index, gate) in circuit.gates.iter().enumerate() {
        validate(circuit, gate)
            .and_then(|()| {
                if matches!(gate, Gate::measure { .. }) {
                    measuring = true;
                } else if measuring {
                    return Err(PbcError::GateAfterMeasurement);
                }
                Ok(())
            })
            .map_err(|cause| PbcError::InvalidInput {
                index,
                cause: Box::new(cause),
            })?;
    }
    let mut output = PbcCircuit::new(circuit.num_qubits, circuit.num_cbits);
    let mut frame = Vec::with_capacity(circuit.num_qubits);
    for q in 0..qubits {
        frame.push(Frame {
            x: output.x(q)?.as_ref(),
            z: output.z(q)?.as_ref(),
        });
    }
    let mut converter = Converter { output, frame };
    for gate in &circuit.gates {
        converter.gate(gate)?;
    }
    Ok(converter.output)
}

fn validate(circuit: &Circuit, gate: &Gate) -> Result<(), PbcError> {
    if matches!(gate, Gate::rz(..)) {
        return Err(PbcError::UnsupportedGate(GateKind::Rz));
    }
    if matches!(gate, Gate::reset(_)) {
        return Err(PbcError::UnsupportedGate(GateKind::Reset));
    }
    let (n, qs) = qubit_operands(gate);
    for (i, &q) in qs[..n].iter().enumerate() {
        if q as usize >= circuit.num_qubits {
            return Err(PbcError::QubitOutOfRange(q));
        }
        if qs[..i].contains(&q) {
            return Err(PbcError::RepeatedOperand);
        }
    }
    if let Gate::measure { cbit, .. } = gate
        && *cbit as usize >= circuit.num_cbits
    {
        return Err(PbcError::ClassicalBitOutOfRange(*cbit));
    }
    Ok(())
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
        // All handles originate in this converter's arena and are immutable.
        self.output.arena.push(super::PauliNode::Product(a, b))
    }

    fn rotate(&mut self, axis: PauliRef, eighths: i64) {
        // Trusted Hermitian construction: axes are conjugated Paulis or
        // products of mutually commuting images of distinct input wires.
        self.output.operations.push(PbcOp::Rotate {
            axis: PauliAxis(axis),
            angle: PauliAngle::new(eighths),
        });
    }

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

    fn gate(&mut self, gate: &Gate) -> Result<(), PbcError> {
        match *gate {
            Gate::h(q) => {
                let f = &mut self.frame[q as usize];
                std::mem::swap(&mut f.x, &mut f.z);
            }
            Gate::x(q) => {
                self.frame[q as usize].z = self.frame[q as usize].z.scaled(Phase::MinusOne)
            }
            Gate::z(q) => {
                self.frame[q as usize].x = self.frame[q as usize].x.scaled(Phase::MinusOne)
            }
            Gate::s(q) | Gate::sdg(q) => {
                let f = self.frame[q as usize];
                // G† X G is -Y for S and +Y for Sdg; Y = iXZ.
                let phase = if matches!(gate, Gate::s(_)) {
                    Phase::MinusI
                } else {
                    Phase::I
                };
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
            Gate::t(q) | Gate::tdg(q) => {
                self.rotate(
                    self.frame[q as usize].z,
                    if matches!(gate, Gate::t(_)) { 1 } else { -1 },
                );
                return Ok(());
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
            } => {
                let t = self.frame[target as usize];
                let d = if matches!(gate, Gate::ccx { .. }) {
                    t.x
                } else {
                    t.z
                };
                self.controlled_controlled(
                    self.frame[control1 as usize].z,
                    self.frame[control2 as usize].z,
                    d,
                );
                return Ok(());
            }
            Gate::measure { qubit, cbit } => {
                self.output
                    .measure(PauliAxis(self.frame[qubit as usize].z), Some(cbit))?;
                return Ok(());
            }
            Gate::rz(..) | Gate::reset(_) => unreachable!("input validated before conversion"),
        }
        // Only Clifford arms reach here. Append once; do not simplify or scan.
        self.output.output_cliffords.push(gate.clone());
        Ok(())
    }
}

#[cfg(test)]
mod tests;
