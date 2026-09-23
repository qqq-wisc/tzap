//! Manually constructed logical PBC; see `to_pbc` for gate-circuit conversion.
use tzap::circuit::Gate;
use tzap::pbc::{PauliAngle, PbcCircuit};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut circuit = PbcCircuit::new(3, 1);
    let x = circuit.x(0)?;
    let z = circuit.z(2)?;
    let xz = circuit.product(x.as_ref(), z.as_ref())?;
    let axis = circuit.hermitian_axis(xz, 1024)?;
    circuit.rotate(axis, PauliAngle::new(1))?;
    let m = circuit.measure(z, Some(0))?;
    circuit.conditional_pauli(x, m)?;
    circuit.conditional_rotate(axis, PauliAngle::new(2), m)?;
    circuit.push_output_clifford(Gate::h(1))?;
    circuit.push_output_clifford(Gate::cnot {
        control: 0,
        target: 2,
    })?;
    print!("{}", circuit.to_ascii()?);
    Ok(())
}
