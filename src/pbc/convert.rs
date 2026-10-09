use std::num::NonZeroUsize;

use crate::circuit::{Circuit, Gate, Qubit};
use crate::pass::Pass;

use super::{
    PauliAngle, PauliAxis, PauliNode, PauliRef, PbcCircuit, PbcError, PbcOp, check_operands,
};

/// Terminal gate-circuit to logical-PBC conversion. Ordinary gate optimization
/// pipelines retain the default `Pass<Circuit>` type.
#[derive(Clone, Copy, Debug, Default)]
pub struct ToPbc {
    /// Weight bound passed to [`to_pbc`].
    pub max_weight: Option<NonZeroUsize>,
}

impl Pass<Result<PbcCircuit, PbcError>> for ToPbc {
    fn name(&self) -> &str {
        "ToPbc"
    }

    fn run(&self, circuit: &Circuit) -> Result<PbcCircuit, PbcError> {
        to_pbc(circuit, self.max_weight)
    }
}

/// Convert Clifford+T, CZ, CCX, CCZ, and measurements to logical PBC.
/// Measurements may appear anywhere: each measures the current image of Z
/// and leaves the frame unchanged. Resets are rejected. Arbitrary Rz must be
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
/// `max_weight` bounds the weight (number of qubits acted on) of every T-type
/// rotation and measurement. The Clifford gates folded into the frame widen
/// later axes; when an axis would exceed the bound, the Clifford gates since
/// the last such point are emitted in place as pi/4 and pi/2 Pauli rotations
/// (`r 2`/`r 4`, weight at most 2) and the frame restarts from the identity.
/// Those Clifford rotations thus have weight at most `max(max_weight, 2)`:
/// an entangling gate needs weight 2. With a bound below 3, CCX and CCZ are decomposed into Clifford+T first,
/// since their own axes have weight 3. The bound also tracks each frame
/// image's support as packed bits, costing O(Q/64) per gate and O(Q) per
/// flush. `None` places every Clifford in the output frame.
///
/// ```
/// use tzap::circuit::{Circuit, Gate};
/// use tzap::pbc::to_pbc;
/// let mut input = Circuit::new(1);
/// input.apply(Gate::h(0));
/// input.apply(Gate::t(0));
/// let pbc = to_pbc(&input, None)?; // X rotation, followed by an H frame
/// assert_eq!(pbc.operations().len(), 1);
/// assert_eq!(pbc.expand(pbc.output_frame().x(0), 100)?.factors, vec![tzap::pbc::Pauli::Z]);
/// # Ok::<(), tzap::pbc::PbcError>(())
/// ```
pub fn to_pbc(circuit: &Circuit, max_weight: Option<NonZeroUsize>) -> Result<PbcCircuit, PbcError> {
    // Validate the caller's indices before lowering can remove a zero Rz or
    // expand one instruction. Certified angles are validated as one-qubit gates.
    u32::try_from(circuit.num_qubits).map_err(|_| PbcError::TooManyQubits)?;
    for (index, gate) in circuit.gates.iter().enumerate() {
        let validation_gate = match gate {
            Gate::rz(a, q) if a.quarter_turns().is_some() => Gate::z(*q),
            _ => gate.clone(),
        };
        validate(circuit, &validation_gate).map_err(|cause| PbcError::InvalidInput {
            index,
            cause: Box::new(cause),
        })?;
    }
    let lowered = crate::angle::lowered_if_needed(circuit);
    let circuit = &lowered;
    if max_weight.is_some_and(|w| w.get() < 3)
        && circuit
            .gates
            .iter()
            .any(|g| matches!(g, Gate::ccx { .. } | Gate::ccz { .. }))
    {
        // Validate first, so errors refer to the caller's instruction indices.
        validate_all(circuit)?;
        let decomposed = crate::decompose::DecomposeToffoli.run(circuit);
        return convert(&decomposed, max_weight);
    }
    convert(circuit, max_weight)
}

fn validate_all(circuit: &Circuit) -> Result<(), PbcError> {
    u32::try_from(circuit.num_qubits).map_err(|_| PbcError::TooManyQubits)?;
    // Validate everything up front, so conversion itself cannot fail or panic.
    for (index, gate) in circuit.gates.iter().enumerate() {
        validate(circuit, gate).map_err(|cause| PbcError::InvalidInput {
            index,
            cause: Box::new(cause),
        })?;
    }
    Ok(())
}

fn convert(circuit: &Circuit, max_weight: Option<NonZeroUsize>) -> Result<PbcCircuit, PbcError> {
    validate_all(circuit)?;
    let n = circuit.num_qubits;
    let output = PbcCircuit::new(n, circuit.num_cbits);
    let cap = max_weight.map(|max_weight| {
        let identity: Vec<(PauliRef, PauliRef)> = (0..n)
            .map(|q| (output.output_frame.x[q], output.output_frame.z[q]))
            .collect();
        Cap {
            max_weight: max_weight.get(),
            supports: Supports::identity(n),
            pending: Vec::new(),
            identity,
        }
    });
    let mut converter = Converter { output, cap };
    for gate in &circuit.gates {
        converter.gate(gate);
    }
    Ok(converter.output)
}

/// Supports of the frame's images, as packed x/z planes (phases are not
/// needed for weights). The image of X_q is at `2q`, of Z_q at `2q + 1`.
struct Supports {
    l: usize,
    planes: Vec<u64>,
}

impl Supports {
    fn identity(n: usize) -> Self {
        let l = n.div_ceil(64).max(1);
        let mut planes = vec![0u64; 2 * n * 2 * l];
        for q in 0..n {
            let (word, bit) = (q / 64, 1u64 << (q % 64));
            planes[2 * q * 2 * l + word] = bit; // X_q: x plane
            planes[(2 * q + 1) * 2 * l + l + word] = bit; // Z_q: z plane
        }
        Self { l, planes }
    }

    fn image(&self, index: usize) -> &[u64] {
        &self.planes[index * 2 * self.l..(index + 1) * 2 * self.l]
    }

    /// `target ^= source`, image-wise: the support of a product of images.
    fn xor_into(&mut self, target: usize, source: usize) {
        let stride = 2 * self.l;
        for j in 0..stride {
            let value = self.planes[source * stride + j];
            self.planes[target * stride + j] ^= value;
        }
    }

    /// Mirror `CliffordFrame::append` on supports.
    fn apply(&mut self, gate: &Gate) {
        let (x, z) = (|q: Qubit| 2 * q as usize, |q: Qubit| 2 * q as usize + 1);
        match *gate {
            Gate::h(q) => {
                let stride = 2 * self.l;
                let (a, b) = (x(q) * stride, z(q) * stride);
                for j in 0..stride {
                    self.planes.swap(a + j, b + j);
                }
            }
            Gate::x(_) | Gate::z(_) => {} // signs only
            Gate::s(q) | Gate::sdg(q) => self.xor_into(x(q), z(q)),
            Gate::cnot { control, target } => {
                self.xor_into(x(control), x(target));
                self.xor_into(z(target), z(control));
            }
            Gate::cz { control, target } => {
                self.xor_into(x(control), z(target));
                self.xor_into(x(target), z(control));
            }
            _ => unreachable!("only Clifford gates update the frame"),
        }
    }

    /// Weight of the product of the given images.
    fn weight(&self, images: &[usize]) -> usize {
        let l = self.l;
        (0..l)
            .map(|j| {
                let (mut x, mut z) = (0u64, 0u64);
                for &i in images {
                    let w = self.image(i);
                    x ^= w[j];
                    z ^= w[l + j];
                }
                (x | z).count_ones() as usize
            })
            .sum()
    }
}

/// State for a weight bound: the frame images' supports, the Clifford gates
/// folded into the frame since the last flush, and the identity frame.
struct Cap {
    max_weight: usize,
    supports: Supports,
    pending: Vec<Gate>,
    /// Arena handles of X_q and Z_q, for emitting flushed gates and
    /// restarting the frame.
    identity: Vec<(PauliRef, PauliRef)>,
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
    cap: Option<Cap>,
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

    /// With a weight bound, flush the frame first if any of the products of
    /// images listed would exceed it. Images are `2q` (X_q) and `2q + 1`
    /// (Z_q).
    fn bound(&mut self, products: &[&[usize]]) {
        let Some(cap) = &self.cap else { return };
        if products
            .iter()
            .any(|images| cap.supports.weight(images) > cap.max_weight)
        {
            self.flush();
        }
    }

    /// Emit the Clifford gates folded into the frame since the last flush as
    /// Pauli rotations about plain Paulis, in order, and restart the frame
    /// from the identity. Before: input = C * emitted, with C those gates;
    /// after: input = I * (emitted ; C).
    fn flush(&mut self) {
        let cap = self.cap.as_mut().expect("flush needs a weight bound");
        let gates = std::mem::take(&mut cap.pending);
        let identity = cap.identity.clone();
        cap.supports = Supports::identity(self.output.num_qubits);
        let (x, z) = (
            |q: Qubit| identity[q as usize].0,
            |q: Qubit| identity[q as usize].1,
        );
        for gate in gates {
            // Each gate as exp(-i k pi/8 P) rotations, up to global phase.
            match gate {
                Gate::h(q) => {
                    self.rotate(z(q), 2);
                    self.rotate(x(q), 2);
                    self.rotate(z(q), 2);
                }
                Gate::x(q) => self.rotate(x(q), 4),
                Gate::z(q) => self.rotate(z(q), 4),
                Gate::s(q) => self.rotate(z(q), 2),
                Gate::sdg(q) => self.rotate(z(q), -2),
                Gate::cnot { control, target } => {
                    let joint = self.product(z(control), x(target));
                    self.rotate(z(control), 2);
                    self.rotate(x(target), 2);
                    self.rotate(joint, -2);
                }
                Gate::cz { control, target } => {
                    let joint = self.product(z(control), z(target));
                    self.rotate(z(control), 2);
                    self.rotate(z(target), 2);
                    self.rotate(joint, -2);
                }
                _ => unreachable!("only Clifford gates are pending"),
            }
        }
        for (q, &(xq, zq)) in identity.iter().enumerate() {
            self.output.output_frame.x[q] = xq;
            self.output.output_frame.z[q] = zq;
        }
    }

    fn gate(&mut self, gate: &Gate) {
        // Keep T-type rotations and measurements within the weight bound.
        let (xi, zi) = (|q: Qubit| 2 * q as usize, |q: Qubit| 2 * q as usize + 1);
        match *gate {
            Gate::t(q) | Gate::tdg(q) | Gate::measure { qubit: q, .. } => self.bound(&[&[zi(q)]]),
            Gate::ccx {
                control1: a,
                control2: b,
                target,
            }
            | Gate::ccz {
                control1: a,
                control2: b,
                target,
            } => {
                let d = if matches!(gate, Gate::ccx { .. }) {
                    xi(target)
                } else {
                    zi(target)
                };
                let (a, b) = (zi(a), zi(b));
                self.bound(&[&[a], &[b], &[d], &[a, b], &[a, d], &[b, d], &[a, b, d]]);
            }
            _ => {}
        }
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
                if let Some(cap) = &mut self.cap {
                    cap.supports.apply(gate);
                    cap.pending.push(gate.clone());
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;
