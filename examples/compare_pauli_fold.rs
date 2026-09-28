//! Compare the existing phase fold with circuit-level Pauli folding.
//! Run with `cargo run --release --example compare_pauli_fold -- [--direct] DIR...`.
//! `--direct` runs each pass independently on the same parsed input.

use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use tzap::cancel::CancelGates;
use tzap::pass::{Pass, count_rz, count_t};
use tzap::pauli_fold_rand::pauli_fold_rand;
use tzap::phase_fold_rand::phase_fold_rand;
use tzap::qasm;

fn collect(path: &Path, files: &mut Vec<PathBuf>) -> io::Result<()> {
    if path.is_dir() {
        for entry in fs::read_dir(path)? {
            collect(&entry?.path(), files)?;
        }
    } else if path.extension().is_some_and(|e| e == "qasm") {
        files.push(path.to_owned());
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut paths: Vec<_> = env::args_os().skip(1).collect();
    let direct = paths.first().is_some_and(|arg| arg == "--direct");
    if direct {
        paths.remove(0);
    }
    let mut files = Vec::new();
    for path in paths {
        collect(Path::new(&path), &mut files)?;
    }
    files.sort();
    if files.is_empty() {
        return Err("provide at least one directory or QASM file".into());
    }
    println!(
        "file,qubits,input_gates,prepared_gates,phase_gates,phase_t,pauli_gates,pauli_t,phase_ms,pauli_ms,eligible,mode,phase_rz,pauli_rz,phase_rotations,pauli_rotations"
    );
    for path in files {
        let source = fs::read_to_string(&path)?;
        let circuit =
            qasm::parse(&source).map_err(|error| format!("{}: {error}", path.display()))?;
        drop(source);
        let cancelled = (!direct).then(|| CancelGates.run(&circuit));
        let phase_input = cancelled.as_ref().unwrap_or(&circuit);
        let started = Instant::now();
        let baseline = phase_fold_rand(phase_input);
        let phase_ms = started.elapsed().as_secs_f64() * 1e3;
        let pauli_input = if direct { &circuit } else { &baseline };
        // Every gate kind is handled (CCX/CCZ block, measurements and resets
        // are barriers), so every circuit is eligible.
        let eligible = true;
        let started = Instant::now();
        let pauli = eligible.then(|| pauli_fold_rand(pauli_input));
        let pauli_ms = started.elapsed().as_secs_f64() * 1e3;
        let result = pauli.as_ref().unwrap_or(pauli_input);
        assert!(count_t(result) + count_rz(result) <= count_t(pauli_input) + count_rz(pauli_input));
        assert!(result.gates.len() <= pauli_input.gates.len());
        println!(
            "{},{},{},{},{},{},{},{},{:.3},{:.3},{},{},{},{},{},{}",
            path.display(),
            circuit.num_qubits,
            circuit.gates.len(),
            phase_input.gates.len(),
            baseline.gates.len(),
            count_t(&baseline),
            result.gates.len(),
            count_t(result),
            phase_ms,
            pauli_ms,
            eligible,
            if direct { "direct" } else { "additive" },
            count_rz(&baseline),
            count_rz(result),
            count_t(&baseline) + count_rz(&baseline),
            count_t(result) + count_rz(result),
        );
        io::stdout().flush()?;
    }
    Ok(())
}
