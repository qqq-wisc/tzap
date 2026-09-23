//! Convert a gate circuit and render its PBC rotations and Clifford suffix.
use tzap::circuit::{Circuit, Gate};
use tzap::pbc::to_pbc;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = Circuit::new(2);
    input.apply(Gate::h(0));
    input.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    input.apply(Gate::t(1));
    let output = to_pbc(&input)?;
    print!("{}", output.to_ascii()?);
    Ok(())
}
