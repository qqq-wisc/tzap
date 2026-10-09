#[path = "support/mod.rs"]
mod support;
use support::Tzap;

fn qasm(n: usize, body: &str) -> String {
    qasm_with_cbits(n, n, body)
}

fn qasm_with_cbits(n: usize, cbits: usize, body: &str) -> String {
    format!("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{n}];\ncreg c[{cbits}];\n{body}\n")
}

fn pbc_guide() -> String {
    include_str!("../docs/pbc.md").replace("\r\n", "\n")
}

/// Convert without optimization beyond gate cancellation, so outputs are the
/// direct transformation the docs describe.
fn convert(source: &str) -> String {
    Tzap::new(&[
        "-",
        "-o",
        "-",
        "--to-pbc",
        "--pbc-no-opt",
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
    let doc = pbc_guide();
    // (qubits, body, expected output, whether docs/pbc.md shows it)
    for (n, body, expected, in_doc) in [
        (
            1,
            "h q[0];\nmeasure q[0] -> c[0];",
            "pbc 0.1\nqubits 1\nregisters 1\nm 1 X0 -> c0\nf X0 1 Z0\nf Z0 1 X0\n",
            true,
        ),
        (
            1,
            "x q[0];\nmeasure q[0] -> c[0];",
            "pbc 0.1\nqubits 1\nregisters 1\nm -1 Z0 -> c0\nf Z0 -1 Z0\n",
            false,
        ),
        (
            1,
            "h q[0];\nt q[0];\nh q[0];\nmeasure q[0] -> c[0];",
            "pbc 0.1\nqubits 1\nregisters 1\nr 1 1 X0\nm 1 Z0 -> c0\n",
            false,
        ),
        (
            2,
            "h q[0];\ncx q[0],q[1];\nmeasure q[0] -> c[0];\nmeasure q[1] -> c[1];",
            "pbc 0.1\nqubits 2\nregisters 2\nm 1 X0 -> c0\nm 1 X0 Z1 -> c1\nf X0 1 Z0 X1\nf Z0 1 X0\nf Z1 1 X0 Z1\n",
            true,
        ),
    ] {
        assert!(!in_doc || doc.contains(&format!("```text\n{expected}```")));
        assert_eq!(convert(&qasm(n, body)), expected);
    }
    for (body, expected) in [
        (
            "h q[0];\ncx q[0],q[1];\nmeasure q[1] -> c[0];",
            "pbc 0.1\nqubits 2\nregisters 1\nm 1 X0 Z1 -> c0\nf X0 1 Z0 X1\nf Z0 1 X0\nf Z1 1 X0 Z1\n",
        ),
        (
            "h q[0];\ncx q[0],q[1];\nt q[1];\nmeasure q[1] -> c[0];",
            "pbc 0.1\nqubits 2\nregisters 1\nr 1 1 X0 Z1\nm 1 X0 Z1 -> c0\nf X0 1 Z0 X1\nf Z0 1 X0\nf Z1 1 X0 Z1\n",
        ),
    ] {
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
        "pbc 0.1\nqubits 3\nregisters 3\nm 1 X2 -> c2\nm 1 Z0 -> c0\nf X2 1 Z2\nf Z2 1 X2\n"
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
        "pbc 0.1\nqubits 1\nregisters 1\nm 1 X0 -> c0\nf X0 1 Z0\nf Z0 1 X0\n"
    );
    assert!(
        run.stderr
            .contains("s\n\t├─ 0 π/8 rotations · 1 measurement\n"),
        "{}",
        run.stderr
    );
    assert!(
        run.stderr
            .contains("\t└─ measurement weight min/median/max: 1/1/1\n"),
        "{}",
        run.stderr
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
        "pbc 0.1\nqubits 2\nregisters 2\nf X0 1 X0 Z1\nf X1 1 Z0 X1\n"
    );
}

#[test]
fn default_pipeline_decomposes_rz_before_conversion() {
    let run = Tzap::new(&["-", "-o", "-", "--to-pbc", "--decompose-rz", "-O1"])
        .stdin(&qasm(1, "h q[0];\nrz(pi/4) q[0];"))
        .run()
        .ok("Rz decomposition");
    assert!(run.stdout.starts_with("pbc 0.1\nqubits 1\nregisters 1\n"));
    assert!(run.stdout.contains("r "));
    assert!(!run.stdout.contains("rz"));
}

#[test]
fn native_ccx_ccz_and_cz_export_without_decomposition_flags() {
    let run = Tzap::new(&[
        "-",
        "-o",
        "-",
        "--to-pbc",
        "--pbc-no-opt",
        "--passes",
        "CancelGates",
    ])
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
        ("reset q[0];", "gate 1 (reset q0): PBC has no resets"),
        (
            "rz(pi/5) q[0];",
            "gate 1 (rz(0.6283) q0): Rz needs --decompose-rz",
        ),
        (
            "h q[1];\ncx q[0],q[0];",
            "line 6: cx expects distinct qubit operands, qubit 0 is repeated",
        ),
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
        "pbc 0.1\nqubits 1\nregisters 1\nm -1 Z0 -> c0\nf Z0 -1 Z0\n"
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
        "pbc 0.1\nqubits 1\nregisters 1\nm -1 Z0 -> c0\nf Z0 -1 Z0\n"
    );
    Tzap::new(&["-", "-o", "-", "--to-pbc", "--json"])
        .run()
        .failed("conflicting stdout writers");
}

/// Mid-circuit examples through the CLI; the one docs/pbc.md shows must match.
#[test]
fn documented_mid_circuit_examples_match_cli_output() {
    let doc = pbc_guide();
    for (n, body, expected) in [
        (
            1,
            "h q[0];\nt q[0];\nmeasure q[0] -> c[0];\nh q[0];\nt q[0];\nmeasure q[0] -> c[1];",
            "pbc 0.1\nqubits 1\nregisters 2\nr 1 1 X0\nm 1 X0 -> c0\nr 1 1 Z0\nm 1 Z0 -> c1\n",
        ),
        (
            2,
            "h q[0];\ncx q[0],q[1];\nmeasure q[0] -> c[0];\nt q[1];\nh q[1];\nmeasure q[1] -> c[1];",
            "pbc 0.1\nqubits 2\nregisters 2\nm 1 X0 -> c0\nr 1 1 X0 Z1\nm 1 X1 -> c1\n\
             f X0 1 Z0 X1\nf X1 1 X0 Z1\nf Z0 1 X0\nf Z1 1 X1\n",
        ),
        (
            3,
            "h q[0];\nt q[0];\ncx q[0],q[2];\ncx q[1],q[2];\nmeasure q[2] -> c[0];\n\
             h q[0];\nt q[0];\nt q[1];\nmeasure q[0] -> c[1];",
            "pbc 0.1\nqubits 3\nregisters 2\nr 1 1 X0\nm 1 X0 Z1 Z2 -> c0\nr 1 1 Z0 X2\nr 1 1 Z1\n\
             m 1 Z0 X2 -> c1\nf X1 1 X1 X2\nf Z0 1 Z0 X2\nf Z2 1 X0 Z1 Z2\n",
        ),
        (
            3,
            "h q[0];\nh q[1];\nccx q[0],q[1],q[2];\nmeasure q[2] -> c[0];\ntdg q[0];\n\
             cx q[0],q[1];\nmeasure q[1] -> c[1];",
            "pbc 0.1\nqubits 3\nregisters 2\nr 1 1 X0\nr 1 1 X1\nr 1 1 X2\nr -1 1 X0 X1\n\
             r -1 1 X0 X2\nr -1 1 X1 X2\nr 1 1 X0 X1 X2\nm 1 Z2 -> c0\nr -1 1 X0\n\
             m 1 X0 X1 -> c1\nf X0 1 Z0 Z1\nf X1 1 Z1\nf Z0 1 X0\nf Z1 1 X0 X1\n",
        ),
    ] {
        let in_doc = n == 2;
        assert!(
            !in_doc || doc.contains(&format!("```text\n{expected}```")),
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

/// --to-pbc optimizes the PBC rotations by default and reports the T counts;
/// --pbc-no-opt skips that, and requires PBC output.
#[test]
fn to_pbc_optimizes_rotations_unless_pbc_no_opt() {
    // T on q0, then CX (its Z image is unchanged), then T on q0 again: the two
    // Z0 rotations merge into one S rotation, which moves into the frame.
    let source = qasm(2, "t q[0];\ncx q[0],q[1];\nt q[0];");
    let run = Tzap::new(&["-", "-o", "-", "--to-pbc"])
        .stdin(&source)
        .run()
        .ok("PBC optimization");
    assert!(!run.stdout.contains("\nr "), "{}", run.stdout);
    assert!(run.stdout.contains("f X0 -1 Y0 X1\n"), "{}", run.stdout);
    // Both T rotations merge and the resulting S moves into the frame.
    assert!(
        run.stderr
            .contains("s\n\t└─ 2 → 0 π/8 rotations (↓100.0%)\n"),
        "{}",
        run.stderr
    );
    assert!(!run.stderr.contains('┌'), "{}", run.stderr);
    let plain = Tzap::new(&["-", "-o", "-", "--to-pbc", "--pbc-no-opt"])
        .stdin(&source)
        .run()
        .ok("plain conversion");
    assert!(
        plain.stderr.contains("s\n\t├─ 2 π/8 rotations\n"),
        "{}",
        plain.stderr
    );
    assert!(
        plain
            .stderr
            .contains("\t└─ π/8 weight min/median/max: 1/1/1\n"),
        "{}",
        plain.stderr
    );
    assert!(!plain.stderr.contains("Optimized PBC"));
    Tzap::new(&["-", "--pbc-no-opt"])
        .stdin(&source)
        .run()
        .failed("--pbc-no-opt without --to-pbc");
}

#[test]
fn invalid_pbc_option_with_long_stdin_reports_the_option_error() {
    let source = qasm(1, &"t q[0];\n".repeat(100_000));
    let run = Tzap::new(&["-", "--pbc-no-opt"])
        .stdin(&source)
        .run()
        .failed("--pbc-no-opt without --to-pbc");
    assert!(run.stderr.contains("--pbc-no-opt"), "{}", run.stderr);
}

/// --to-pbc alone turns gate-level optimization off and says so; -O* or
/// --passes turn it back on. Requested decompositions still run.
#[test]
fn to_pbc_turns_gate_optimization_off_unless_requested() {
    // A cancelling H pair and a T/Tdg pair: gate passes remove all four.
    let source = qasm(1, "h q[0];\nh q[0];\nt q[0];\ntdg q[0];");
    let notice = "Gate-level optimization is off with --to-pbc";
    let off = Tzap::new(&["-", "-o", "-", "--to-pbc"])
        .stdin(&source)
        .run()
        .ok("--to-pbc");
    assert!(off.stderr.contains(notice), "{}", off.stderr);
    // No gate reduction box: the circuit's figures, then the PBC's.
    assert!(
        off.stderr
            .contains("\t├─ 1 qubit · 4 gates\n\t├─ 0 2q gates · 2 T/Tdg · 4 depth\n"),
        "{}",
        off.stderr
    );
    assert!(!off.stderr.contains("Final result"), "{}", off.stderr);
    assert!(
        off.stderr.contains("s\n\t├─ 2 π/8 rotations"),
        "{}",
        off.stderr
    );
    for flags in [&["-O3"][..], &["--passes", "CancelGates,PhaseFoldRand"]] {
        let mut args = vec!["-", "-o", "-", "--to-pbc"];
        args.extend_from_slice(flags);
        let on = Tzap::new(&args).stdin(&source).run().ok("explicit");
        assert!(!on.stderr.contains(notice), "{flags:?}: {}", on.stderr);
        assert!(on.stderr.contains("4 → 0"), "{flags:?}: {}", on.stderr);
    }
    // Without --to-pbc, the default O3 still runs.
    let qasm_out = Tzap::new(&["-", "-o", "-"]).stdin(&source).run().ok("O3");
    assert!(!qasm_out.stderr.contains(notice));
    assert!(qasm_out.stderr.contains("4 → 0"), "{}", qasm_out.stderr);
    // Rz still needs --decompose-rz, which still runs.
    let rz = Tzap::new(&["-", "-o", "-", "--to-pbc", "--decompose-rz"])
        .stdin(&qasm(1, "h q[0];\nrz(pi/4) q[0];"))
        .run()
        .ok("Rz decomposition");
    assert!(rz.stdout.starts_with("pbc 0.1\nqubits 1\nregisters 1\n"));
    assert!(rz.stdout.contains("r "), "{}", rz.stdout);
}

/// --visualize-pbc writes an SVG drawing, alone or together with --to-pbc.
#[test]
fn visualize_pbc_writes_an_svg() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("circuit.svg");
    let run = Tzap::new(&[
        "-",
        "--visualize-pbc",
        path.to_str().unwrap(),
        "--passes",
        "CancelGates",
    ])
    .stdin(&qasm(
        2,
        "h q[0];\ncx q[0],q[1];\nt q[1];\nmeasure q[1] -> c[0];",
    ))
    .run()
    .ok("drawing");
    assert!(run.stderr.contains("wrote PBC drawing"), "{}", run.stderr);
    let svg = std::fs::read_to_string(&path).unwrap();
    assert!(svg.starts_with("<svg") && svg.contains("</svg>"));
    // One green T rotation box and one blue measurement box.
    assert!(svg.contains("#E3FFA1") && svg.contains("#70B3F5"));
    Tzap::new(&["-", "--visualize-pbc"])
        .run()
        .failed("missing path");
}

/// --pbc-max-weight emits the Cliffords that would widen an axis as
/// rotations; it requires N >= 1 and --to-pbc or --visualize-pbc.
#[test]
fn pbc_max_weight_flushes_cliffords_as_rotations() {
    let doc = pbc_guide();
    let source = qasm(2, "h q[0];\ncx q[0],q[1];\nt q[1];\nmeasure q[1] -> c[0];");
    let expected = "pbc 0.1\nqubits 2\nregisters 2\nr 2 1 Z0\nr 2 1 X0\nr 4 1 Z0\nr 2 1 X1\n\
                    r -2 1 Z0 X1\nr 1 1 Z1\nm 1 Z1 -> c0\n";
    assert!(
        doc.contains(expected),
        "docs/pbc.md is missing:\n{expected}"
    );
    let run = Tzap::new(&[
        "-",
        "-o",
        "-",
        "--to-pbc",
        "--pbc-max-weight",
        "1",
        "--passes",
        "CancelGates",
    ])
    .stdin(&source)
    .run()
    .ok("bounded conversion");
    assert_eq!(run.stdout, expected);
    assert!(
        run.stderr
            .contains("s\n\t├─ 1 π/8 rotation · 6 Clifford rotations · 1 measurement\n"),
        "{}",
        run.stderr
    );
    assert!(
        run.stderr.contains(
            "s\n\t├─ 1 → 1 π/8 rotation (↓0.0%)\n\t├─ 6 → 5 Clifford rotations (↓16.7%)\n\t\
             ├─ π/8 weight min/median/max: 1/1/1\n\t└─ Clifford weight min/median/max: 1/1/2\n"
        ),
        "{}",
        run.stderr
    );
    // Bound 2 admits the unbounded conversion unchanged.
    let wide = Tzap::new(&[
        "-",
        "-o",
        "-",
        "--to-pbc",
        "--pbc-max-weight",
        "2",
        "--passes",
        "CancelGates",
    ])
    .stdin(&source)
    .run()
    .ok("bound 2");
    assert_eq!(wide.stdout, convert(&source));
    for bad in [
        &["-", "--to-pbc", "--pbc-max-weight", "0"][..],
        &["-", "--to-pbc", "--pbc-max-weight"],
        &["-", "--pbc-max-weight", "2"],
    ] {
        Tzap::new(bad)
            .stdin(&source)
            .run()
            .failed("invalid --pbc-max-weight");
    }
}

/// ToPbc and PbcOpt are passes: listed after the gate passes, they match
/// --to-pbc (with and without --pbc-no-opt), and may run with no gate passes.
#[test]
fn to_pbc_and_pbc_opt_are_passes() {
    let source = qasm(2, "t q[0];\ncx q[0],q[1];\nt q[0];\nh q[1];\nh q[1];");
    let run = |args: &[&str]| {
        let mut all = vec!["-", "-o", "-"];
        all.extend_from_slice(args);
        Tzap::new(&all).stdin(&source).run().ok("PBC passes")
    };
    let flags = run(&["--passes", "CancelGates", "--to-pbc"]);
    let passes = run(&["--passes", "CancelGates,ToPbc,PbcOpt"]);
    assert_eq!(passes.stdout, flags.stdout);
    assert!(
        passes
            .stderr
            .contains("s\n\t└─ 2 → 0 π/8 rotations (↓100.0%)"),
        "{}",
        passes.stderr
    );
    assert_eq!(
        run(&["--passes", "CancelGates", "--to-pbc", "--pbc-no-opt"]).stdout,
        run(&["--passes", "CancelGates,ToPbc"]).stdout
    );
    // No gate passes: the uncancelled H pair folds into the frame, where it
    // cancels, leaving the two unmerged T rotations.
    let bare = run(&["--passes", "ToPbc"]);
    assert_eq!(
        bare.stdout,
        "pbc 0.1\nqubits 2\nregisters 2\nr 1 1 Z0\nr 1 1 Z0\nf X0 1 X0 X1\nf Z1 1 Z0 Z1\n"
    );
    let optimized = run(&["--passes", "ToPbc,PbcOpt"]);
    assert!(!optimized.stdout.contains("\nr "), "{}", optimized.stdout);
    for bad in [
        "PbcOpt",
        "ToPbc,CancelGates",
        "ToPbc,ToPbc",
        "ToPbc,PbcOpt,PbcOpt",
    ] {
        Tzap::new(&["-", "--passes", bad])
            .stdin(&source)
            .run()
            .failed(bad);
    }
    for flag in ["--to-pbc", "--pbc-no-opt"] {
        Tzap::new(&["-", "--passes", "ToPbc", flag])
            .stdin(&source)
            .run()
            .failed(flag);
    }
}

/// Invalid gate operands are rejected by the parser before any PBC or
/// gate-level work; --fixpoint needs gate passes to repeat.
#[test]
fn pbc_input_is_checked_before_optimization() {
    let run = Tzap::new(&["-", "--to-pbc", "-O3"])
        .stdin(&qasm(2, "h q[0];\nt q[0];\ncx q[1],q[1];"))
        .run()
        .failed("repeated operand");
    assert!(
        run.stderr
            .contains("line 7: cx expects distinct qubit operands, qubit 1 is repeated"),
        "{}",
        run.stderr
    );
    assert!(!run.stderr.contains("Gate-level"), "{}", run.stderr);
    // Rz is fine when a pass decomposes it.
    Tzap::new(&["-", "--to-pbc", "--decompose-rz"])
        .stdin(&qasm(1, "rz(pi/4) q[0];"))
        .run()
        .ok("decomposed Rz");
    let fixpoint = Tzap::new(&["-", "--to-pbc", "--fixpoint"])
        .stdin(&qasm(1, "t q[0];"))
        .run()
        .failed("--fixpoint without gate passes");
    assert!(fixpoint.stderr.contains("use -O3"), "{}", fixpoint.stderr);
    Tzap::new(&["-", "--to-pbc", "--passes", "CancelGates", "--fixpoint"])
        .stdin(&qasm(1, "t q[0];"))
        .run()
        .ok("--fixpoint with --passes");
}

/// The report: one line when the optimizer changes nothing, notices that
/// match the flags, and the main output written before the drawing.
#[test]
fn pbc_report_lines_match_the_run() {
    let source = qasm(1, "h q[0];\nt q[0];\nmeasure q[0] -> c[0];");
    let run = Tzap::new(&["-", "--to-pbc", "--pbc-no-opt", "--superopt-gates", "base"])
        .stdin(&source)
        .run()
        .ok("no optimizer");
    assert!(
        run.stderr
            .contains("Gate-level optimization is off with --to-pbc; pass -O3"),
        "{}",
        run.stderr
    );
    assert!(
        run.stderr
            .contains("--superopt-* has no effect without gate-level optimization"),
        "{}",
        run.stderr
    );
    let optimized = Tzap::new(&["-", "--to-pbc"])
        .stdin(&source)
        .run()
        .ok("nothing to optimize");
    assert!(
        optimized.stderr.contains("(the PBC is optimized instead)"),
        "{}",
        optimized.stderr
    );
    assert!(
        optimized.stderr.contains("s\n\t└─ no reduction"),
        "{}",
        optimized.stderr
    );
    assert!(!optimized.stderr.contains('┌'), "{}", optimized.stderr);

    let dir = tempfile::tempdir().unwrap();
    let pbc = dir.path().join("out.pbc");
    let bad_svg = dir.path().join("missing").join("x.svg");
    Tzap::new(&[
        "-",
        "--to-pbc",
        "-o",
        pbc.to_str().unwrap(),
        "--visualize-pbc",
        bad_svg.to_str().unwrap(),
    ])
    .stdin(&source)
    .run()
    .failed("unwritable drawing");
    assert!(
        std::fs::read_to_string(&pbc)
            .unwrap()
            .starts_with("pbc 0.1\nqubits 1\n")
    );
}

/// --json describes the PBC: counts, weights, and the optimizer's work.
#[test]
fn json_reports_the_pbc() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.pbc");
    let source = qasm(2, "t q[0];\ncx q[0],q[1];\nt q[0];\nmeasure q[1] -> c[0];");
    let json = |args: &[&str]| {
        let mut all = vec!["-", "-o", path.to_str().unwrap(), "--json"];
        all.extend_from_slice(args);
        Tzap::new(&all).stdin(&source).run().ok("--json").stdout
    };
    let off = json(&["--to-pbc"]);
    for expected in [
        "\"level\": null",
        "\"gate_set\": null",
        "\"pi8_rotations\": 2",
        "\"pi8_rotations\": 0",
        "\"measurements\": 1",
        "\"merges\": 1",
        "\"cliffords_to_frame\": 1",
        "\"max_weight\": null",
    ] {
        assert!(off.contains(expected), "missing {expected}:\n{off}");
    }
    let passes = json(&["--passes", "CancelGates,ToPbc"]);
    assert!(
        passes.contains("\"CancelGates\",\n      \"ToPbc\""),
        "{passes}"
    );
    assert!(passes.contains("\"optimization\": null"), "{passes}");
    let qasm_only = json(&[]);
    assert!(qasm_only.contains("\"pbc\": null"), "{qasm_only}");
}

/// The report's counts and weights agree however they are obtained: from
/// the optimizer, from the exported text, or by materializing for a drawing.
#[test]
fn report_weights_agree_across_sources() {
    let source = qasm(
        3,
        "h q[0];\ncx q[0],q[1];\nt q[1];\ncx q[1],q[2];\nt q[2];\nh q[2];\nt q[2];\n\
         measure q[2] -> c[0];\nt q[0];",
    );
    let converted = |stderr: &str| -> String {
        let start = stderr.find("Converted to PBC").expect(stderr);
        let block = &stderr[start..];
        let end = block.find("\n  ").unwrap_or(block.len());
        // Drop the timing on the heading line.
        block[..end].split_once('\n').unwrap().1.to_string()
    };
    let dir = tempfile::tempdir().unwrap();
    let svg = dir.path().join("c.svg");
    let run = |args: &[&str]| {
        let mut all = vec!["-"];
        all.extend_from_slice(args);
        Tzap::new(&all).stdin(&source).run().ok("report").stderr
    };
    let optimizer = converted(&run(&["--to-pbc"]));
    let text = converted(&run(&["--to-pbc", "--pbc-no-opt"]));
    let drawing = converted(&run(&[
        "--visualize-pbc",
        svg.to_str().unwrap(),
        "--pbc-no-opt",
        "--passes",
        "CancelGates",
    ]));
    assert_eq!(optimizer, text);
    assert_eq!(text, drawing);
    assert!(text.contains("measurement weight min/median/max"), "{text}");
}

/// --pbc-expansion-budget bounds writing, and needs PBC output.
#[test]
fn expansion_budget_flag() {
    let source = qasm(2, "t q[0];\nh q[0];\nt q[0];\ncx q[0],q[1];\nt q[1];");
    let over = Tzap::new(&["-", "--to-pbc", "-o", "-", "--pbc-expansion-budget", "1"])
        .stdin(&source)
        .run()
        .failed("tiny budget");
    assert!(
        over.stderr.contains("raise it with --pbc-expansion-budget"),
        "{}",
        over.stderr
    );
    for bad in [
        &["-", "--pbc-expansion-budget", "5"][..],
        &["-", "--to-pbc", "--pbc-expansion-budget", "lots"],
    ] {
        Tzap::new(bad).stdin(&source).run().failed("bad flag use");
    }
}
