#[path = "support/mod.rs"]
mod support;
use support::Tzap;

fn qasm(n: usize, body: &str) -> String {
    format!("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{n}];\ncreg c[{n}];\n{body}\n")
}

#[test]
fn documented_measurement_examples_match_cli_output() {
    let doc = include_str!("../docs/pbc.md");
    for (n, body, expected) in [
        (
            1,
            "h q[0];\nmeasure q[0] -> c[0];",
            "qubits 1\nregisters 1\nm 1 X0 -> c0\n",
        ),
        (
            1,
            "x q[0];\nmeasure q[0] -> c[0];",
            "qubits 1\nregisters 1\nm -1 Z0 -> c0\n",
        ),
        (
            1,
            "h q[0];\nt q[0];\nh q[0];\nmeasure q[0] -> c[0];",
            "qubits 1\nregisters 1\nr 1 1 X0\nm 1 Z0 -> c0\n",
        ),
        (
            2,
            "h q[0];\ncx q[0],q[1];\nmeasure q[0] -> c[0];\nmeasure q[1] -> c[1];",
            "qubits 2\nregisters 2\nm 1 X0 -> c0\nm 1 X0 Z1 -> c1\n",
        ),
    ] {
        assert!(doc.contains(&format!("```text\n{expected}```")));
        let run = Tzap::new(&[
            "-",
            "-o",
            "-",
            "--to-pbc",
            "--passes",
            "CancelGates",
            "--quiet",
        ])
        .stdin(&qasm(n, body))
        .run()
        .ok("documented example");
        assert_eq!(run.stdout, expected);
    }
}

#[test]
fn stdout_has_only_pbc_and_full_readout_omits_suffix() {
    let run = Tzap::new(&["-", "-o", "-", "--to-pbc", "--passes", "CancelGates"])
        .stdin(&qasm(1, "h q[0];\nmeasure q[0] -> c[0];"))
        .run()
        .ok("full readout");
    assert_eq!(run.stdout, "qubits 1\nregisters 1\nm 1 X0 -> c0\n");
    assert!(run.stderr.contains("classical outputs only"));
}

#[test]
fn partial_repeated_and_absent_readout_lower_remaining_cliffords() {
    for body in [
        "h q[0];",
        "h q[0];\nmeasure q[0] -> c[0];",
        "h q[0];\nmeasure q[0] -> c[0];\nmeasure q[0] -> c[1];",
    ] {
        let run = Tzap::new(&[
            "-",
            "-o",
            "-",
            "--to-pbc",
            "--passes",
            "CancelGates",
            "--quiet",
        ])
        .stdin(&qasm(2, body))
        .run()
        .ok("partial readout");
        assert!(run.stdout.ends_with("r 2 1 Z0\nr 2 1 X0\nr 2 1 Z0\n"));
        assert!(!run.stdout.contains("suffix"));
        assert!(run.stderr.is_empty());
    }
}

#[test]
fn conversion_runs_after_optimization_and_custom_decomposition() {
    let run = Tzap::new(&[
        "-",
        "-o",
        "-",
        "--to-pbc",
        "--passes",
        "CancelGates,DecomposeCz",
    ])
    .stdin(&qasm(2, "h q[0];\nh q[0];\ncz q[0],q[1];"))
    .run()
    .ok("final transformation");
    assert_eq!(
        run.stdout,
        "qubits 2\nregisters 2\nr 2 1 Z1\nr 2 1 X1\nr 2 1 Z1\nr 2 1 Z0\nr 2 1 X1\nr -2 1 Z0 X1\nr 2 1 Z1\nr 2 1 X1\nr 2 1 Z1\n"
    );
}

#[test]
fn default_pipeline_decomposes_rz_before_conversion() {
    let run = Tzap::new(&["-", "-o", "-", "--to-pbc", "--decompose-rz", "-O1"])
        .stdin(&qasm(1, "h q[0];\nrz(pi/4) q[0];"))
        .run()
        .ok("Rz decomposition");
    assert!(run.stdout.starts_with("qubits 1\nregisters 1\n"));
    assert!(run.stdout.contains("r "));
    assert!(!run.stdout.contains("rz"));
}

#[test]
fn native_ccx_ccz_and_cz_export_without_decomposition_flags() {
    let run = Tzap::new(&["-", "-o", "-", "--to-pbc", "--passes", "CancelGates"])
        .stdin(&qasm(
            3,
            "ccx q[0],q[1],q[2];\nccz q[0],q[1],q[2];\ncz q[0],q[1];",
        ))
        .run()
        .ok("native multi-qubit gates");
    assert_eq!(
        run.stdout
            .lines()
            .filter(|line| line.starts_with("r "))
            .count(),
        17
    );
    assert!(run.stdout.ends_with("r 2 1 Z0\nr 2 1 Z1\nr -2 1 Z0 Z1\n"));
    assert!(run.stdout.contains("r 1 1 Z0 Z1 X2\n"));
    assert!(run.stdout.contains("r 1 1 Z0 Z1 Z2\n"));
}

#[test]
fn unsupported_inputs_fail_even_without_output_destination() {
    for (body, error) in [
        ("reset q[0];", "Reset"),
        ("rz(pi/5) q[0];", "--decompose-rz"),
        ("measure q[0] -> c[0];\nh q[0];", "final input block"),
    ] {
        let run = Tzap::new(&["-", "--to-pbc", "--passes", "CancelGates"])
            .stdin(&qasm(1, body))
            .run()
            .failed("invalid PBC input");
        assert!(run.stderr.contains(error), "{}", run.stderr);
    }
}

#[test]
fn file_output_json_and_errors_preserve_stream_contract() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("output.pbc");
    let run = Tzap::new(&[
        "-",
        "-o",
        path.to_str().unwrap(),
        "--to-pbc",
        "--passes",
        "CancelGates",
        "--json",
    ])
    .stdin(&qasm(1, "x q[0];\nmeasure q[0] -> c[0];"))
    .run()
    .ok("PBC file and JSON");
    assert!(run.stdout.trim_start().starts_with('{'));
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "qubits 1\nregisters 1\nm -1 Z0 -> c0\n"
    );
    Tzap::new(&[
        "-",
        "-o",
        path.to_str().unwrap(),
        "--to-pbc",
        "--passes",
        "CancelGates",
    ])
    .stdin(&qasm(1, "reset q[0];"))
    .run()
    .failed("no overwrite on conversion failure");
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "qubits 1\nregisters 1\nm -1 Z0 -> c0\n"
    );
    Tzap::new(&["-", "-o", "-", "--to-pbc", "--json"])
        .run()
        .failed("conflicting stdout writers");
}
