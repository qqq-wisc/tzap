//! Experimental driver for phase folding over quantum WHILE programs.
//!
//!     tzap-prog <pass> <program.qasm> [--print]
//!
//! `<pass>` is `none`, `pauli` (PhaseFoldPauli over programs; `pauli-disj`
//! with up to 100 disjuncts, `pauli-zero` with eigenstate facts from resets and
//! measurements, `pauli-full` with both, `pauli-full<k>` with both and up to k
//! disjuncts), or
//! `statefold1`, `statefold2`, `statefold0` (StateFold over programs; 0 is
//! unbounded degree). Prints the static T-count before and after, and the
//! optimized program with `--print`.

use std::time::Instant;

use tzap::program::{self, Program};
use tzap::program::pauli_fold::Options;

fn run(pass: &str, prog: &Program) -> Result<Program, String> {
    Ok(match pass {
        "none" => prog.clone(),
        "pauli" => program::pauli_fold::fold(prog),
        "pauli-disj" => program::pauli_fold::fold_with(prog, Options { disjuncts: 100, zero_facts: false }),
        "pauli-zero" => program::pauli_fold::fold_with(prog, Options { disjuncts: 1, zero_facts: true }),
        "pauli-full" => program::pauli_fold::fold_with(prog, Options { disjuncts: 100, zero_facts: true }),
        "statefold1" => program::state_fold::fold(prog, Some(1)),
        "statefold2" => program::state_fold::fold(prog, Some(2)),
        "statefold0" => program::state_fold::fold(prog, None),
        // `pauli-full<k>`: eigenstate facts and up to k disjuncts.
        other => match other.strip_prefix("pauli-full").and_then(|k| k.parse().ok()) {
            Some(disjuncts) => program::pauli_fold::fold_with(prog, Options { disjuncts, zero_facts: true }),
            None => return Err(format!("unknown pass {other}")),
        },
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: tzap-prog <none|pauli|pauli-disj|pauli-zero|pauli-full[<k>]|statefold1|statefold2|statefold0> <program.qasm> [--print]");
        std::process::exit(2);
    }
    let src = std::fs::read_to_string(&args[2]).unwrap_or_else(|e| {
        eprintln!("{}: {e}", args[2]);
        std::process::exit(1)
    });
    let prog = program::parse(&src).unwrap_or_else(|e| {
        eprintln!("{}: {e}", args[2]);
        std::process::exit(1)
    });
    let start = Instant::now();
    let out = run(&args[1], &prog).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(2)
    });
    let secs = start.elapsed().as_secs_f64();
    println!("{{\"qubits\": {}, \"t_in\": {}, \"t_out\": {}, \"seconds\": {secs}}}", prog.num_qubits, prog.t_count(), out.t_count());
    if args.iter().any(|a| a == "--print") {
        print!("{}", out.to_qasm3());
    }
}
