//! Generate the HWB-like circuit family used by the POPL 2027 RQ3 scaling
//! experiment.
//!
//! Hadamards mint fresh parity variables, while CNOTs accumulate them on a
//! qubit that never receives a Hadamard. After `s` stages, the accumulator's
//! parity contains `s + 1` variables and the circuit contains `5s` gates.

use std::env;
use std::fs;
use std::process;

use tzap::circuit::{Circuit, Gate, Qubit};

fn generate(num_qubits: usize, stages: usize) -> Circuit {
    assert!(
        num_qubits >= 2,
        "need at least 2 qubits (1 accumulator + 1 worker)"
    );

    let mut circuit = Circuit::new(num_qubits);
    let accumulator: Qubit = 0;
    for stage in 0..stages {
        let worker = Qubit::try_from(1 + stage % (num_qubits - 1))
            .expect("qubit index exceeds the supported range");
        circuit.apply(Gate::h(worker));
        circuit.apply(Gate::cnot {
            control: worker,
            target: accumulator,
        });
        circuit.apply(Gate::t(accumulator));
        circuit.apply(Gate::cnot {
            control: accumulator,
            target: worker,
        });
        circuit.apply(Gate::tdg(worker));
    }
    circuit
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let (Some(num_qubits), Some(stages)) = (
        args.get(1).and_then(|value| value.parse::<usize>().ok()),
        args.get(2).and_then(|value| value.parse::<usize>().ok()),
    ) else {
        eprintln!("usage: gen_hwb_like <num_qubits> <stages> [output.qasm]");
        process::exit(1);
    };

    let qasm = generate(num_qubits, stages).to_qasm();
    match args.get(3) {
        Some(path) => fs::write(path, qasm).unwrap_or_else(|error| {
            eprintln!("error writing {path}: {error}");
            process::exit(1);
        }),
        None => print!("{qasm}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_stage_emits_the_five_gate_rq3_pattern() {
        let circuit = generate(4, 7);
        assert_eq!(circuit.num_qubits, 4);
        assert_eq!(circuit.gates().len(), 35);
        assert!(matches!(circuit.gates()[0], Gate::h(1)));
        assert!(matches!(circuit.gates()[5], Gate::h(2)));
        assert!(matches!(circuit.gates()[10], Gate::h(3)));
        assert!(matches!(circuit.gates()[15], Gate::h(1)));
    }
}
