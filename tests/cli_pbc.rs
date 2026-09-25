#[path = "support/mod.rs"]
mod support;
use support::Tzap;

fn qasm(n: usize, body: &str) -> String {
    qasm_with_cbits(n, n, body)
}

fn qasm_with_cbits(n: usize, cbits: usize, body: &str) -> String {
    format!("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{n}];\ncreg c[{cbits}];\n{body}\n")
}

/// Convert without optimization beyond gate cancellation, so outputs are the
/// direct transformation the docs describe.
fn convert(source: &str) -> String {
    Tzap::new(&[
        "-",
        "-o",
        "-",
        "--to-pbc",
        "--passes",
        "CancelGates",
        "--quiet",
    ])
    .stdin(source)
    .run()
    .ok("PBC conversion")
    .stdout
}

#[test]
fn documented_measurement_examples_match_cli_output() {
    let doc = include_str!("../docs/pbc.md");
    let readme = include_str!("../README.md");
    for (n, body, expected) in [
        (
            1,
            "h q[0];\nmeasure q[0] -> c[0];",
            "qubits 1\nregisters 1\nm 1 X0 -> c0\nf X0 1 Z0\nf Z0 1 X0\n",
        ),
        (
            1,
            "x q[0];\nmeasure q[0] -> c[0];",
            "qubits 1\nregisters 1\nm -1 Z0 -> c0\nf Z0 -1 Z0\n",
        ),
        (
            1,
            "h q[0];\nt q[0];\nh q[0];\nmeasure q[0] -> c[0];",
            "qubits 1\nregisters 1\nr 1 1 X0\nm 1 Z0 -> c0\n",
        ),
        (
            2,
            "h q[0];\ncx q[0],q[1];\nmeasure q[0] -> c[0];\nmeasure q[1] -> c[1];",
            "qubits 2\nregisters 2\nm 1 X0 -> c0\nm 1 X0 Z1 -> c1\nf X0 1 Z0 X1\nf Z0 1 X0\nf Z1 1 X0 Z1\n",
        ),
    ] {
        assert!(doc.contains(&format!("```text\n{expected}```")));
        assert_eq!(convert(&qasm(n, body)), expected);
    }
    for (body, expected, text) in [
        (
            "h q[0];\ncx q[0],q[1];\nmeasure q[1] -> c[0];",
            "qubits 2\nregisters 1\nm 1 X0 Z1 -> c0\nf X0 1 Z0 X1\nf Z0 1 X0\nf Z1 1 X0 Z1\n",
            doc,
        ),
        (
            "h q[0];\ncx q[0],q[1];\nt q[1];\nmeasure q[1] -> c[0];",
            "qubits 2\nregisters 1\nr 1 1 X0 Z1\nm 1 X0 Z1 -> c0\nf X0 1 Z0 X1\nf Z0 1 X0\nf Z1 1 X0 Z1\n",
            readme,
        ),
    ] {
        assert!(text.contains(&format!("```text\n{expected}```")));
        assert_eq!(convert(&qasm_with_cbits(2, 1, body)), expected);
    }
}

/// Documented: registers from several QASM declarations are numbered
/// consecutively, in declaration order.
#[test]
fn multiple_registers_are_numbered_in_declaration_order() {
    let source = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg a[1];\nqreg b[2];\n\
                  creg x[1];\ncreg y[2];\nh b[1];\nmeasure b[1] -> y[1];\nmeasure a[0] -> x[0];\n";
    assert_eq!(
        convert(source),
        "qubits 3\nregisters 3\nm 1 X2 -> c2\nm 1 Z0 -> c0\nf X2 1 Z2\nf Z2 1 X2\n"
    );
}

#[test]
fn stdout_has_only_pbc_and_full_readout_retains_frame() {
    let run = Tzap::new(&["-", "-o", "-", "--to-pbc", "--passes", "CancelGates"])
        .stdin(&qasm(1, "h q[0];\nmeasure q[0] -> c[0];"))
        .run()
        .ok("full readout");
    assert_eq!(
        run.stdout,
        "qubits 1\nregisters 1\nm 1 X0 -> c0\nf X0 1 Z0\nf Z0 1 X0\n"
    );
    assert!(
        run.stderr
            .contains("retaining quantum and classical outputs")
    );
}

#[test]
fn partial_repeated_and_absent_readout_keep_frames() {
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
        assert!(run.stdout.ends_with("f X0 1 Z0\nf Z0 1 X0\n"));
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
        "qubits 2\nregisters 2\nf X0 1 X0 Z1\nf X1 1 Z0 X1\n"
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
        14
    );
    assert!(run.stdout.ends_with("f X0 1 X0 Z1\nf X1 1 Z0 X1\n"));
    assert!(run.stdout.contains("r 1 1 Z0 Z1 X2\n"));
    assert!(run.stdout.contains("r 1 1 Z0 Z1 Z2\n"));
}

#[test]
fn unsupported_inputs_fail_even_without_output_destination() {
    for (body, error) in [
        ("reset q[0];", "Reset"),
        ("rz(pi/5) q[0];", "--decompose-rz"),
    ] {
        let run = Tzap::new(&["-", "--to-pbc", "--passes", "CancelGates"])
            .stdin(&qasm(2, body))
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
        "qubits 1\nregisters 1\nm -1 Z0 -> c0\nf Z0 -1 Z0\n"
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
        "qubits 1\nregisters 1\nm -1 Z0 -> c0\nf Z0 -1 Z0\n"
    );
    Tzap::new(&["-", "-o", "-", "--to-pbc", "--json"])
        .run()
        .failed("conflicting stdout writers");
}

/// Every mid-circuit example in docs/pbc.md, through the CLI.
#[test]
fn documented_mid_circuit_examples_match_cli_output() {
    let doc = include_str!("../docs/pbc.md");
    for (n, body, expected) in [
        (
            1,
            "h q[0];\nt q[0];\nmeasure q[0] -> c[0];\nh q[0];\nt q[0];\nmeasure q[0] -> c[1];",
            "qubits 1\nregisters 2\nr 1 1 X0\nm 1 X0 -> c0\nr 1 1 Z0\nm 1 Z0 -> c1\n",
        ),
        (
            2,
            "h q[0];\ncx q[0],q[1];\nmeasure q[0] -> c[0];\nt q[1];\nh q[1];\nmeasure q[1] -> c[1];",
            "qubits 2\nregisters 2\nm 1 X0 -> c0\nr 1 1 X0 Z1\nm 1 X1 -> c1\n\
             f X0 1 Z0 X1\nf X1 1 X0 Z1\nf Z0 1 X0\nf Z1 1 X1\n",
        ),
        (
            3,
            "h q[0];\nt q[0];\ncx q[0],q[2];\ncx q[1],q[2];\nmeasure q[2] -> c[0];\n\
             h q[0];\nt q[0];\nt q[1];\nmeasure q[0] -> c[1];",
            "qubits 3\nregisters 2\nr 1 1 X0\nm 1 X0 Z1 Z2 -> c0\nr 1 1 Z0 X2\nr 1 1 Z1\n\
             m 1 Z0 X2 -> c1\nf X1 1 X1 X2\nf Z0 1 Z0 X2\nf Z2 1 X0 Z1 Z2\n",
        ),
        (
            3,
            "h q[0];\nh q[1];\nccx q[0],q[1],q[2];\nmeasure q[2] -> c[0];\ntdg q[0];\n\
             cx q[0],q[1];\nmeasure q[1] -> c[1];",
            "qubits 3\nregisters 2\nr 1 1 X0\nr 1 1 X1\nr 1 1 X2\nr -1 1 X0 X1\n\
             r -1 1 X0 X2\nr -1 1 X1 X2\nr 1 1 X0 X1 X2\nm 1 Z2 -> c0\nr -1 1 X0\n\
             m 1 X0 X1 -> c1\nf X0 1 Z0 Z1\nf X1 1 Z1\nf Z0 1 X0\nf Z1 1 X0 X1\n",
        ),
    ] {
        assert!(
            doc.contains(&format!("```text\n{expected}```")),
            "{expected}"
        );
        assert_eq!(convert(&qasm_with_cbits(n, 2, body)), expected);
    }
}

/// The default pipeline optimizes around mid-circuit measurements, and the
/// converted output keeps a measurement followed by further rotations.
#[test]
fn default_pipeline_converts_mid_circuit_measurements() {
    let run = Tzap::new(&["-", "-o", "-", "--to-pbc", "--quiet"])
        .stdin(&qasm(
            2,
            "h q[0];\nt q[0];\nmeasure q[0] -> c[0];\nh q[0];\nt q[0];\ncx q[0],q[1];\nmeasure q[1] -> c[1];",
        ))
        .run()
        .ok("mid-circuit measurement");
    let lines: Vec<_> = run.stdout.lines().collect();
    let first = lines.iter().position(|l| l.starts_with("m ")).unwrap();
    assert!(
        lines[first + 1..].iter().any(|l| l.starts_with("r ")),
        "{}",
        run.stdout
    );
}
