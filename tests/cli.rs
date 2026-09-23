use std::fs;
use std::process::Command;

#[path = "support/mod.rs"]
mod support;

use support::Json;

const TEST_QASM: &str = "tests/fixtures/test.qasm";
const TWO_CCX_QASM: &str = "tests/fixtures/two_ccx.qasm";
const MOD5_4_QASM: &str = "benchmarks/feynman/mod5_4.qasm";

fn tzap() -> Command {
    Command::new(env!("CARGO_BIN_EXE_tzap"))
}

fn tzap_run(args: &[&str]) -> std::process::Output {
    tzap()
        .args(args)
        .args([
            "--superopt-qubits",
            "3",
            "--superopt-window-gates",
            "4",
            "--superopt-murm-entries",
            "500",
        ])
        .output()
        .expect("failed to run tzap")
}

/// The `--json` report for a run. How these tests check what a run *did*
/// internally — how many rounds it took, whether it converged — now that the
/// progress box those facts used to be read out of is drawn only on a live
/// terminal. See `tests/cli_json.rs` for the report's own tests.
fn tzap_report(args: &[&str]) -> Json {
    let mut full = vec!["--json", "-q"];
    full.extend_from_slice(args);
    let out = tzap_run(&full);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    Json::parse(&String::from_utf8_lossy(&out.stdout))
}

#[test]
fn no_args_prints_usage() {
    let out = tzap_run(&[]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("tzap"),
        "expected usage message, got: {stderr}"
    );
}

#[test]
fn missing_file_errors() {
    let out = tzap_run(&["nonexistent.qasm"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Error reading"),
        "expected error message, got: {stderr}"
    );
}

#[test]
fn parsed_circuit_gate_set_uses_canonical_order() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("all-gates.qasm");
    fs::write(
        &input,
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\ncreg c[1];\n\
reset q[2];\nmeasure q[1] -> c[0];\nccz q[0],q[1],q[2];\nccx q[0],q[1],q[2];\n\
cz q[0],q[1];\ncx q[0],q[1];\nrz(pi/7) q[0];\ntdg q[0];\nt q[0];\nsdg q[0];\n\
s q[0];\nz q[0];\nx q[0];\nh q[0];\n",
    )
    .unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "--passes", "CancelGates"]);
    assert_success(&out, "canonical circuit gate set");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains(
            "Circuit gates: {h, x, z, s, sdg, t, tdg, rz, cx, cz, ccx, ccz, measure, reset}"
        ),
        "got:\n{stderr}"
    );
}

#[test]
fn murm_reports_auto_base_and_explicit_synthesis_bases() {
    for (argument, expected) in [
        (
            "auto",
            "Synthesis basis: {h, x, z, s, sdg, t, tdg, cx, ccx}",
        ),
        ("base", "Synthesis basis: {h, x, z, s, sdg, t, tdg, cx}"),
        ("h,t,cx,cz", "Synthesis basis: {h, t, cx, cz}"),
    ] {
        let out = tzap_run(&[TWO_CCX_QASM, "-O2", "--superopt-gates", argument]);
        assert_success(&out, argument);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains("Loaded MURM in "), "got:\n{stderr}");
        assert!(stderr.contains(expected), "{argument}: got:\n{stderr}");
    }
}

/// A circuit big enough to trip the automatic parallel threshold, written to
/// `path`. One qubit and one gate per line keeps it cheap to generate and to
/// parse; what matters here is the gate count.
fn write_big_circuit(path: &std::path::Path, gates: usize) {
    let mut qasm = String::with_capacity(gates * 8 + 64);
    qasm.push_str("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\n");
    for _ in 0..gates {
        qasm.push_str("h q[0];\n");
    }
    fs::write(path, qasm).unwrap();
}

/// A big circuit goes parallel without being asked: the size alone decides,
/// the run says so, and the report agrees.
#[test]
fn a_large_circuit_is_optimized_in_parallel_by_default() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("big.qasm");
    write_big_circuit(&input, 1_000_000);
    let input = input.to_str().unwrap();

    let out = tzap_run(&[input, "--passes", "CancelGates"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Optimizing chunks in parallel"),
        "expected the automatic-parallel note:\n{stderr}"
    );
    assert!(
        tzap_report(&[input, "--passes", "CancelGates"])
            .at("options/parallel")
            .as_bool()
    );
}

/// `--no-parallel` overrides that, and stays quiet about it.
#[test]
fn no_parallel_keeps_a_large_circuit_sequential() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("big.qasm");
    write_big_circuit(&input, 1_000_000);
    let input = input.to_str().unwrap();

    let out = tzap_run(&[input, "--no-parallel", "--passes", "CancelGates"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("parallel"),
        "--no-parallel must not announce a parallel run:\n{stderr}"
    );
    assert!(
        !tzap_report(&[input, "--no-parallel", "--passes", "CancelGates"])
            .at("options/parallel")
            .as_bool()
    );
}

/// A circuit below the threshold is sequential, and nothing is said about it.
#[test]
fn a_small_circuit_stays_sequential_and_silent() {
    let out = tzap_run(&[MOD5_4_QASM]);
    assert!(out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!stderr.contains("parallel"), "got:\n{stderr}");
    assert!(!tzap_report(&[MOD5_4_QASM]).at("options/parallel").as_bool());
}

#[test]
fn optimizes_test_qasm_to_file() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&[TEST_QASM, "-o", out_path.to_str().unwrap()]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let content = fs::read_to_string(&out_path).unwrap();
    assert!(content.starts_with("OPENQASM 2.0;"));
    assert!(content.contains("qreg q[3];"));

    // A test harness captures stderr through a pipe, where tzap emits no
    // escapes at all — so the banner arrives as plain text (see
    // `cli_streams.rs` for that rule and `src/main.rs` for the styled form).
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains(&format!("tzap v{}", env!("CARGO_PKG_VERSION"))));
    assert!(!stderr.contains('\x1b'), "got: {stderr}");
    assert!(stderr.contains("test.qasm"));
    assert!(stderr.contains("gates"));
    assert!(stderr.contains("T/Tdg"));
}

#[test]
fn output_is_valid_qasm() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&[TEST_QASM, "-o", out_path.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&out_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "OPENQASM 2.0;");
    assert_eq!(lines[1], "include \"qelib1.inc\";");
    assert!(lines[2].starts_with("qreg q["));

    // Every gate line should end with semicolon
    for line in &lines[3..] {
        if !line.is_empty() {
            assert!(line.ends_with(';'), "gate line missing semicolon: {line}");
        }
    }
}

#[test]
fn writes_to_output_file() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&[TEST_QASM, out_path.to_str().unwrap()]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // stdout should be empty when writing to file
    assert!(
        out.stdout.is_empty(),
        "stdout should be empty when output file given"
    );

    let content = fs::read_to_string(&out_path).unwrap();
    assert!(content.starts_with("OPENQASM 2.0;"));

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("wrote"));
}

#[test]
fn roundtrip_preserves_qasm_structure() {
    // Run tzap, then feed the output back through tzap — should be stable
    let dir = tempfile::tempdir().unwrap();
    let pass1 = dir.path().join("pass1.qasm");
    let pass2 = dir.path().join("pass2.qasm");

    let out1 = tzap_run(&[TEST_QASM, pass1.to_str().unwrap()]);
    assert!(out1.status.success());

    let out2 = tzap_run(&[pass1.to_str().unwrap(), pass2.to_str().unwrap()]);
    assert!(out2.status.success());

    let content1 = fs::read_to_string(&pass1).unwrap();
    let content2 = fs::read_to_string(&pass2).unwrap();
    assert_eq!(content1, content2, "second pass should be idempotent");
}

#[test]
fn toffoli_decomposition_increases_gate_count() {
    // two_ccx.qasm has 3 gates (2 CCX + 1 CX) — output should have more after decomposition
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&[
        TWO_CCX_QASM,
        "-o",
        out_path.to_str().unwrap(),
        "--decompose-ccx",
    ]);
    assert!(out.status.success());

    let content = fs::read_to_string(&out_path).unwrap();
    let gates = gate_lines_from(&content);
    assert!(
        gates.len() > 3,
        "decomposed circuit should have more than 3 gates, got {}",
        gates.len()
    );
    // No CCX gates should remain
    assert!(
        !gates.iter().any(|g| g.starts_with("ccx ")),
        "no ccx gates should remain after decomposition"
    );
}

#[test]
fn explicit_superopt_basis_cannot_reintroduce_requested_decompositions() {
    let dir = tempfile::tempdir().unwrap();
    for (name, gate, decompose, forbidden) in [
        ("cz", "cz q[0],q[1];", "--decompose-cz", "cz "),
        ("ccx", "ccx q[0],q[1],q[2];", "--decompose-ccx", "ccx "),
        ("ccz", "ccz q[0],q[1],q[2];", "--decompose-ccx", "ccz "),
    ] {
        let input = dir.path().join(format!("{name}.qasm"));
        let output = dir.path().join(format!("{name}-out.qasm"));
        fs::write(
            &input,
            format!("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\n{gate}\n"),
        )
        .unwrap();
        let run = tzap_run(&[
            input.to_str().unwrap(),
            "-O2",
            decompose,
            "--superopt-gates",
            name,
            "-o",
            output.to_str().unwrap(),
        ]);
        assert!(
            run.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let gates = gate_lines_from(&fs::read_to_string(output).unwrap());
        assert!(
            !gates.iter().any(|line| line.starts_with(forbidden)),
            "{name} was reintroduced: {gates:?}"
        );
    }
}

/// The decomposition line must show both sides of the gate count. The final
/// result banner measures its reduction against the *decomposed* circuit, so
/// printing only the post-decomposition figure left readers no way to see
/// where it came from — a 3-gate input reported as "31 → 21 gates" with no
/// mention of the 3.
#[test]
fn decomposition_line_shows_both_gate_counts() {
    let out = tzap_run(&[TWO_CCX_QASM, "--decompose-ccx"]);
    assert!(out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert!(stderr.contains("Toffoli decomposition"), "got: {stderr}");
    assert!(
        stderr.contains("3 → "),
        "expected the decomposition line to show the input gate count (3) \
         alongside the decomposed one:\n{stderr}"
    );
}

/// The result banner's heading carries the headline gate reduction, so the
/// number most readers want is visible without reading the rows.
#[test]
fn final_result_heading_leads_with_the_gate_reduction() {
    let out = tzap_run(&[MOD5_4_QASM]);
    assert!(out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert!(
        stderr.contains("Final result · ") && stderr.contains("% fewer gates · "),
        "expected the heading to lead with the gate reduction, ahead of the \
         elapsed time:\n{stderr}"
    );
}

/// A sub-megabyte input must not report its size as "0.0 MB", which reads as
/// an empty file.
#[test]
fn small_input_size_is_not_reported_as_zero_megabytes() {
    let out = tzap_run(&[TWO_CCX_QASM]);
    assert!(out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert!(
        !stderr.contains("0.0 MB"),
        "a tiny input should report bytes, not 0.0 MB:\n{stderr}"
    );
    assert!(stderr.contains(" B) in "), "got: {stderr}");
}

#[test]
fn ccz_is_native_by_default() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[3];
ccz q[0],q[1],q[2];
",
    );

    assert_eq!(gates, ["ccz q[0],q[1],q[2];"]);
}

#[test]
fn mod5_4_reduces_t_count() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&[MOD5_4_QASM, "-o", out_path.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&out_path).unwrap();
    let gates = gate_lines_from(&content);

    let t_count = gates
        .iter()
        .filter(|g| g.starts_with("t ") || g.starts_with("tdg "))
        .count();

    assert_eq!(
        t_count, 16,
        "mod5_4 should optimize to 16 T/Tdg, got {t_count}"
    );
    assert!(gates.len() < 79, "mod5_4 gate count did not decrease");
}

#[test]
fn mod5_4_idempotent() {
    let dir = tempfile::tempdir().unwrap();
    let pass1 = dir.path().join("pass1.qasm");
    let pass2 = dir.path().join("pass2.qasm");

    let out1 = tzap_run(&[MOD5_4_QASM, pass1.to_str().unwrap()]);
    assert!(out1.status.success());

    let out2 = tzap_run(&[pass1.to_str().unwrap(), pass2.to_str().unwrap()]);
    assert!(out2.status.success());

    let c1 = fs::read_to_string(&pass1).unwrap();
    let c2 = fs::read_to_string(&pass2).unwrap();
    assert_eq!(c1, c2, "mod5_4 output should be stable on second pass");
}

#[test]
fn mod5_4_no_rz_in_output() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");
    let out = tzap_run(&[MOD5_4_QASM, "-o", out_path.to_str().unwrap()]);
    assert!(out.status.success());
    let content = fs::read_to_string(&out_path).unwrap();
    for line in content.lines() {
        assert!(
            !line.starts_with("rz("),
            "mod5_4 output should not contain raw rz gates, found: {line}"
        );
    }
}

#[test]
fn inline_qasm_optimization() {
    // Write a known circuit to a temp file and verify T+T folds to S
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("tt.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(
        &input,
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
",
    )
    .unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "-o", output.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&output).unwrap();
    let gate_lines: Vec<&str> = content
        .lines()
        .filter(|l| {
            !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg")
        })
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(
        gate_lines.len(),
        1,
        "T+T should fold to S, got: {gate_lines:?}"
    );
    assert!(
        gate_lines[0].starts_with("s "),
        "T+T should fold to S, got: {gate_lines:?}"
    );
}

#[test]
fn t_t_folds_to_s() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("tt.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(
        &input,
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
",
    )
    .unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "-o", output.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&output).unwrap();
    let gate_lines: Vec<&str> = content
        .lines()
        .filter(|l| {
            !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg")
        })
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(
        gate_lines.len(),
        1,
        "T+T should fold to single gate, got: {gate_lines:?}"
    );
    assert_eq!(gate_lines[0], "s q[0];", "T+T should fold to S");
}

#[test]
fn t_tdg_cancels() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("t_tdg.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(
        &input,
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
tdg q[0];
",
    )
    .unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "-o", output.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&output).unwrap();
    let gate_lines: Vec<&str> = content
        .lines()
        .filter(|l| {
            !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg")
        })
        .filter(|l| !l.is_empty())
        .collect();
    assert!(
        gate_lines.is_empty(),
        "T+Tdg should cancel, got: {gate_lines:?}"
    );
}

#[test]
fn phase_fold_through_cnot_control() {
    // T q0; CNOT q0,q1; T q0 should merge to CNOT + S
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
t q[0];
cx q[0],q[1];
t q[0];
",
    );
    assert_eq!(
        gates.len(),
        2,
        "should merge T through CNOT control: {gates:?}"
    );
    assert!(
        gates.iter().any(|g| g.starts_with("cx ")),
        "CNOT should remain"
    );
    assert!(gates.iter().any(|g| g == "s q[0];"), "T+T should become S");
}

#[test]
fn phase_fold_cancel_through_cnot_control() {
    // S q0; CNOT q0,q1; Sdg q0 — should cancel rotations
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
s q[0];
cx q[0],q[1];
sdg q[0];
",
    );
    assert_eq!(
        gates.len(),
        1,
        "S and Sdg should cancel through CNOT: {gates:?}"
    );
    assert!(gates[0].starts_with("cx "), "only CNOT should remain");
}

#[test]
fn phase_fold_blocked_on_cnot_target() {
    // T q[1]; CNOT q[0],q[1]; T q[1] — q1 is target, can't merge
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
t q[1];
cx q[0],q[1];
t q[1];
",
    );
    // Phase fold may or may not touch this, but the two T's should not merge
    assert!(
        gates.len() >= 3,
        "should not merge T on CNOT target: {gates:?}"
    );
}

fn rz_qasm(theta_expr: &str) -> String {
    format!("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz({theta_expr}) q[0];\n")
}

#[test]
fn decompose_rz_with_epsilon_produces_cliffordt() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "1e-3",
        "--superopt-gates",
        "base",
    ]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let content = fs::read_to_string(&output).unwrap();
    let gates = gate_lines_from(&content);
    assert!(
        !gates.iter().any(|g| g.starts_with("rz(")),
        "no rz gates should remain, got: {gates:?}"
    );
    assert!(
        !gates.is_empty(),
        "output should have gates after decomposition"
    );
}

#[test]
fn epsilon_accepts_scientific_notation_variants() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    for eps in &["1e-3", "1E-3", "1.5e-2", "0.001"] {
        let out = tzap_run(&[input.to_str().unwrap(), "--decompose-rz", "--epsilon", eps]);
        assert!(
            out.status.success(),
            "--epsilon {eps} should be accepted, stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
}

#[test]
fn invalid_epsilon_exits_with_error() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "not-a-number",
    ]);
    assert!(!out.status.success(), "invalid epsilon should fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--epsilon"),
        "error should mention --epsilon, got: {stderr}"
    );
}

#[test]
fn epsilon_without_decompose_rz_is_ignored() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--epsilon",
        "1e-3",
    ]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let content = fs::read_to_string(&output).unwrap();
    let gates = gate_lines_from(&content);
    assert!(
        gates.iter().any(|g| g.starts_with("rz(")),
        "rz gate should be preserved without --decompose-rz, got: {gates:?}"
    );
}

#[test]
fn coarser_epsilon_produces_fewer_t_gates_than_finer() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    let out_fine = dir.path().join("fine.qasm");
    let out_coarse = dir.path().join("coarse.qasm");

    let r1 = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        out_fine.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "1e-4",
    ]);
    let r2 = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        out_coarse.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "1e-2",
    ]);
    assert!(r1.status.success());
    assert!(r2.status.success());

    let t_fine = gate_lines_from(&fs::read_to_string(&out_fine).unwrap())
        .iter()
        .filter(|g| g.starts_with("t ") || g.starts_with("tdg "))
        .count();
    let t_coarse = gate_lines_from(&fs::read_to_string(&out_coarse).unwrap())
        .iter()
        .filter(|g| g.starts_with("t ") || g.starts_with("tdg "))
        .count();

    assert!(
        t_coarse <= t_fine,
        "coarser epsilon should need <= T gates: coarse={t_coarse}, fine={t_fine}"
    );
}

// --- --passes pipeline override ---

const HHTT_QASM: &str =
    "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0];\nh q[0];\nt q[0];\nt q[0];\n";

#[test]
fn passes_runs_only_the_listed_pass() {
    // CancelGates alone removes the H·H pair but does not fold T·T (that is the
    // phase-folding pass), so the two T gates survive.
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, HHTT_QASM).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--passes",
        "CancelGates",
    ]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let gates = gate_lines_from(&fs::read_to_string(&output).unwrap());
    assert_eq!(gates, vec!["t q[0];", "t q[0];"]);
}

#[test]
fn passes_run_in_the_given_order() {
    // CancelGates removes H·H leaving T·T, then PhaseFoldRand folds T·T into S.
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, HHTT_QASM).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--passes",
        "CancelGates,PhaseFoldRand",
    ]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let gates = gate_lines_from(&fs::read_to_string(&output).unwrap());
    assert_eq!(gates, vec!["s q[0];"]);
}

#[test]
fn passes_allows_decompose_rz_with_epsilon() {
    // --epsilon is permitted with --passes and configures the DecomposeRz pass.
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--passes",
        "DecomposeRz,CancelGates,PhaseFoldRand",
        "--epsilon",
        "1e-3",
    ]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let gates = gate_lines_from(&fs::read_to_string(&output).unwrap());
    assert!(
        !gates.iter().any(|g| g.starts_with("rz(")),
        "no rz gates should remain, got: {gates:?}"
    );
    assert!(!gates.is_empty());
}

#[test]
fn passes_conflicts_with_decompose_rz() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    fs::write(&input, HHTT_QASM).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "--passes",
        "CancelGates",
        "--decompose-rz",
    ]);
    assert!(
        !out.status.success(),
        "should reject --passes with --decompose-rz"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("cannot be combined"), "got: {stderr}");
}

#[test]
fn passes_conflicts_with_decompose_cz() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    fs::write(&input, HHTT_QASM).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "--passes",
        "CancelGates",
        "--decompose-cz",
    ]);
    assert!(
        !out.status.success(),
        "should reject --passes with --decompose-cz"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("cannot be combined"), "got: {stderr}");
}

#[test]
fn passes_unknown_name_errors_with_valid_list() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    fs::write(&input, HHTT_QASM).unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "--passes", "Foo,CancelGates"]);
    assert!(!out.status.success(), "unknown pass should fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("Unknown pass 'Foo'"), "got: {stderr}");
    assert!(
        stderr.contains("CancelGates"),
        "error should list valid passes, got: {stderr}"
    );
}

#[test]
fn superopt_pass_builds_murm_behind_cli() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, HHTT_QASM).unwrap();

    let out = tzap_run(&[
        "--passes",
        "PhaseFoldRand,",
        "SuperOpt",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        gate_lines_from(&fs::read_to_string(output).unwrap()),
        ["s q[0];"]
    );
}

// --- optimization levels ---

#[test]
fn o3_is_the_default_pipeline() {
    let dir = tempfile::tempdir().unwrap();
    let default_output = dir.path().join("default.qasm");
    let o3_output = dir.path().join("o3.qasm");
    let default = tzap_run(&[TEST_QASM, "-o", default_output.to_str().unwrap()]);
    let o3 = tzap_run(&[TEST_QASM, "-o", o3_output.to_str().unwrap(), "-O3"]);
    assert!(default.status.success());
    assert!(o3.status.success());
    assert_eq!(
        fs::read_to_string(default_output).unwrap(),
        fs::read_to_string(o3_output).unwrap()
    );

    let stderr = String::from_utf8_lossy(&o3.stderr);
    assert!(
        stderr.contains("Gates") && stderr.contains("T/Tdg"),
        "expected the result box's metric rows:\n{stderr}"
    );
    assert!(
        stderr.contains("Loaded MURM"),
        "O3 uses SuperOpt, so it should build/load the MURM:\n{stderr}"
    );
    assert!(
        stderr.contains("Converged after"),
        "O3 runs to a true fixpoint:\n{stderr}"
    );
    assert!(
        tzap_report(&[TEST_QASM, "-O3"]).get("fixpoints").arr()[0]
            .get("converged")
            .as_bool(),
        "and reports that it did"
    );
}

#[test]
fn o2_uses_superopt_pass_capped_at_two_rounds() {
    let out = tzap_run(&[TEST_QASM, "-O2"]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("Loaded MURM"), "got: {stderr}");

    // O2 may converge sooner, but it must never exceed its two-round cap.
    let report = tzap_report(&[MOD5_4_QASM, "-O2"]);
    assert!(
        report.get("fixpoints").arr()[0].get("rounds").as_usize() <= 2,
        "O2 exceeded its two-round cap"
    );
}

#[test]
fn o3_reports_native_and_post_rz_fixpoints() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-O3",
        "--decompose-rz",
        "--epsilon",
        "1e-3",
    ]);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Gates") && stderr.contains("T/Tdg") && stderr.contains('%'),
        "the result box should carry the gate/T reduction rows:\n{stderr}"
    );
    assert!(
        !stderr.contains("Gate cancellation") && !stderr.contains("Phase folding"),
        "fixpoint progress should not print optimization-pass logs:\n{stderr}"
    );
    assert!(stderr.contains("Rz → Clifford+T decomposition"));
    assert!(stderr.contains("Converged after"), "got: {stderr}");

    let report = tzap_report(&[
        input.to_str().unwrap(),
        "-O3",
        "--decompose-rz",
        "--epsilon",
        "1e-3",
        "--superopt-gates",
        "base",
    ]);
    let fixpoints = report.get("fixpoints").arr();
    assert_eq!(fixpoints.len(), 2, "native and post-decomposition stages");
    assert_eq!(fixpoints[0].get("stage").as_usize(), 0);
    assert_eq!(fixpoints[1].get("stage").as_usize(), 2);
    assert!(fixpoints[1].get("rounds").as_usize() >= 1);

    let gates = gate_lines_from(&fs::read_to_string(output).unwrap());
    assert!(
        !gates.iter().any(|g| g.starts_with("rz(")),
        "no rz gates should remain, got: {gates:?}"
    );
}

#[test]
fn optimization_levels_conflict_with_passes_and_fixpoint() {
    for level in ["-O1", "-O2", "-O3"] {
        for conflicting in [
            ["--passes", "CancelGates"].as_slice(),
            ["--fixpoint"].as_slice(),
        ] {
            let mut args = vec![TEST_QASM, level];
            args.extend_from_slice(conflicting);
            let out = tzap_run(&args);
            assert!(
                !out.status.success(),
                "should reject {level} with {conflicting:?}"
            );
            let stderr = String::from_utf8_lossy(&out.stderr);
            assert!(
                stderr.contains("cannot be combined with --passes or --fixpoint"),
                "got: {stderr}"
            );
        }
    }
}

#[test]
fn optimization_levels_are_mutually_exclusive() {
    let out = tzap_run(&[TEST_QASM, "-O1", "-O2"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("-O1, -O2, -O3, and -Osuper cannot be combined"),
        "got: {stderr}"
    );
}

// --- pass and mode combinations ---

/// Every pass selectable via `--passes`, in default pipeline order.
const ALL_PASS_NAMES: [&str; 7] = [
    "DecomposeToffoli",
    "DecomposeCz",
    "DecomposeRz",
    "CancelGates",
    "SuperOpt",
    "PhaseFoldRand",
    "CnotMin",
];

/// Exercises every decomposition pass: a Toffoli, a CZ, and an Rz rotation,
/// plus cancellable H·H and foldable T·T pairs.
const MIXED_QASM: &str = "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[3];
ccx q[0],q[1],q[2];
cz q[0],q[1];
rz(pi/5) q[2];
h q[0];
h q[0];
t q[1];
t q[1];
";

/// Read an output file, assert it parses as QASM, and return its gate lines.
fn read_valid_qasm(path: &std::path::Path) -> Vec<String> {
    let content = fs::read_to_string(path).unwrap();
    assert!(
        tzap::qasm::parse(&content).is_ok(),
        "output must be valid QASM:\n{content}"
    );
    gate_lines_from(&content)
}

fn assert_success(out: &std::process::Output, context: &str) {
    assert!(
        out.status.success(),
        "{context} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn every_pass_runs_standalone() {
    let dir = tempfile::tempdir().unwrap();
    let mixed = dir.path().join("mixed.qasm");
    fs::write(&mixed, MIXED_QASM).unwrap();

    for pass in ALL_PASS_NAMES {
        // Decomposition passes get gates to decompose; the optimization
        // passes run on the Clifford+T fixture.
        let input = if pass.starts_with("Decompose") {
            mixed.to_str().unwrap().to_string()
        } else {
            TEST_QASM.to_string()
        };
        let output = dir.path().join(format!("{pass}.qasm"));
        let mut args = vec![
            input.as_str(),
            "-o",
            output.to_str().unwrap(),
            "--passes",
            pass,
        ];
        if pass == "DecomposeRz" {
            args.extend_from_slice(&["--epsilon", "1e-3"]);
        }

        let out = tzap_run(&args);
        assert_success(&out, &format!("--passes {pass}"));

        let gates = read_valid_qasm(&output);
        let survives = |prefix: &str| gates.iter().any(|g| g.starts_with(prefix));
        match pass {
            "DecomposeToffoli" => assert!(!survives("ccx "), "{pass} left ccx: {gates:?}"),
            "DecomposeCz" => assert!(!survives("cz "), "{pass} left cz: {gates:?}"),
            "DecomposeRz" => assert!(!survives("rz("), "{pass} left rz: {gates:?}"),
            _ => {}
        }
    }
}

#[test]
fn all_passes_combined_pipeline() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("mixed.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, MIXED_QASM).unwrap();

    let list = ALL_PASS_NAMES.join(",");
    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--passes",
        &list,
        "--superopt-gates",
        "base",
        "--epsilon",
        "1e-3",
    ]);
    assert_success(&out, "--passes with all passes");

    // A base-only MURM cannot reintroduce native controlled gates after the
    // explicit decomposition passes.
    let gates = read_valid_qasm(&output);
    for prefix in ["ccx ", "ccz ", "cz ", "rz("] {
        assert!(
            !gates.iter().any(|g| g.starts_with(prefix)),
            "combined pipeline left {prefix}: {gates:?}"
        );
    }
}

#[test]
fn passes_with_fixpoint_reaches_fixpoint() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, HHTT_QASM).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--passes",
        "CancelGates,PhaseFoldRand",
        "--fixpoint",
    ]);
    assert_success(&out, "--passes with --fixpoint");

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("Converged after"), "got: {stderr}");
    assert_eq!(read_valid_qasm(&output), vec!["s q[0];"]);
}

#[test]
fn passes_with_parallel_produces_valid_output() {
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.qasm");

    let out = tzap_run(&[
        MOD5_4_QASM,
        "-o",
        output.to_str().unwrap(),
        "--passes",
        "CancelGates,PhaseFoldRand",
        "--parallel",
    ]);
    assert_success(&out, "--passes with --parallel");
    assert!(!read_valid_qasm(&output).is_empty());
    assert!(
        tzap_report(&[
            MOD5_4_QASM,
            "--passes",
            "CancelGates,PhaseFoldRand",
            "--parallel",
        ])
        .at("options/parallel")
        .as_bool()
    );
}

#[test]
fn passes_with_parallel_and_fixpoint() {
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.qasm");

    let out = tzap_run(&[
        MOD5_4_QASM,
        "-o",
        output.to_str().unwrap(),
        "--passes",
        "CancelGates,PhaseFoldRand",
        "--fixpoint",
        "--parallel",
    ]);
    assert_success(&out, "--passes with --fixpoint --parallel");
    assert!(!read_valid_qasm(&output).is_empty());
}

#[test]
fn optimization_levels_run_with_parallel() {
    let dir = tempfile::tempdir().unwrap();
    for level in ["-O1", "-O2", "-O3"] {
        let output = dir.path().join(format!("{level}.qasm"));
        let out = tzap_run(&[
            MOD5_4_QASM,
            "-o",
            output.to_str().unwrap(),
            level,
            "--parallel",
        ]);
        assert_success(&out, &format!("{level} --parallel"));

        let stderr = String::from_utf8_lossy(&out.stderr);
        if level != "-O1" {
            assert!(stderr.contains("Loaded MURM"), "got: {stderr}");
        }
        assert!(!read_valid_qasm(&output).is_empty());
    }
}

#[test]
fn default_pipeline_with_parallel() {
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.qasm");

    let out = tzap_run(&[TEST_QASM, "-o", output.to_str().unwrap(), "--parallel"]);
    assert_success(&out, "default --parallel");
    read_valid_qasm(&output);
}

#[test]
fn default_pipeline_with_fixpoint() {
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.qasm");

    let out = tzap_run(&[TEST_QASM, "-o", output.to_str().unwrap(), "--fixpoint"]);
    assert_success(&out, "default --fixpoint");

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("Converged after"), "got: {stderr}");
    read_valid_qasm(&output);
}

#[test]
fn decompose_rz_with_parallel() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, rz_qasm("pi/5")).unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "1e-3",
        "--parallel",
    ]);
    assert_success(&out, "--decompose-rz --parallel");

    let gates = read_valid_qasm(&output);
    assert!(
        !gates.iter().any(|g| g.starts_with("rz(")),
        "no rz gates should remain, got: {gates:?}"
    );
}

fn gate_lines_from(stdout: &str) -> Vec<String> {
    stdout
        .lines()
        .filter(|l| {
            !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg")
        })
        .filter(|l| !l.is_empty())
        .map(|l| l.to_string())
        .collect()
}

fn run_qasm(qasm: &str) -> (Vec<String>, String) {
    run_qasm_with_args(qasm, &[])
}

fn run_qasm_with_args(qasm: &str, extra_args: &[&str]) -> (Vec<String>, String) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.qasm");
    let output = dir.path().join("output.qasm");
    fs::write(&input, qasm).unwrap();
    let mut args = vec![input.to_str().unwrap(), "-o", output.to_str().unwrap()];
    args.extend_from_slice(extra_args);
    let out = tzap_run(&args);
    assert!(
        out.status.success(),
        "tzap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let content = fs::read_to_string(&output).unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).to_string();
    (gate_lines_from(&content), stderr)
}

#[test]
fn s_sdg_cancels() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
sdg q[0];
",
    );
    assert!(gates.is_empty(), "S+Sdg should cancel, got: {gates:?}");
}

#[test]
fn sdg_s_cancels() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
sdg q[0];
s q[0];
",
    );
    assert!(gates.is_empty(), "Sdg+S should cancel, got: {gates:?}");
}

#[test]
fn z_z_cancels() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
z q[0];
z q[0];
",
    );
    assert!(gates.is_empty(), "Z+Z should cancel, got: {gates:?}");
}

#[test]
fn s_s_becomes_z() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
s q[0];
",
    );
    assert_eq!(
        gates.len(),
        1,
        "S+S should fold to one gate, got: {gates:?}"
    );
    assert_eq!(gates[0], "z q[0];");
}

#[test]
fn sdg_sdg_becomes_z() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
sdg q[0];
sdg q[0];
",
    );
    assert_eq!(
        gates.len(),
        1,
        "Sdg+Sdg should fold to one gate, got: {gates:?}"
    );
    assert_eq!(gates[0], "z q[0];");
}

#[test]
fn z_s_becomes_sdg() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
z q[0];
s q[0];
",
    );
    assert_eq!(gates.len(), 1, "Z+S should fold to Sdg, got: {gates:?}");
    assert_eq!(gates[0], "sdg q[0];");
}

#[test]
fn z_sdg_becomes_s() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
z q[0];
sdg q[0];
",
    );
    assert_eq!(gates.len(), 1, "Z+Sdg should fold to S, got: {gates:?}");
    assert_eq!(gates[0], "s q[0];");
}

#[test]
fn sdg_t_becomes_tdg() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
sdg q[0];
t q[0];
",
    );
    assert_eq!(gates.len(), 1, "Sdg+T should fold to Tdg, got: {gates:?}");
    assert_eq!(gates[0], "tdg q[0];");
}

#[test]
fn s_tdg_becomes_t() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
tdg q[0];
",
    );
    assert_eq!(gates.len(), 1, "S+Tdg should fold to T, got: {gates:?}");
    assert_eq!(gates[0], "t q[0];");
}

#[test]
fn six_t_becomes_sdg() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
",
    );
    assert_eq!(gates.len(), 1, "6T should fold to Sdg, got: {gates:?}");
    assert_eq!(gates[0], "sdg q[0];");
}

#[test]
fn four_t_becomes_z() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
t q[0];
t q[0];
",
    );
    assert_eq!(gates.len(), 1, "4T should fold to Z, got: {gates:?}");
    assert_eq!(gates[0], "z q[0];");
}

#[test]
fn four_s_cancels() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
s q[0];
s q[0];
s q[0];
",
    );
    assert!(
        gates.is_empty(),
        "4S should cancel to identity, got: {gates:?}"
    );
}

#[test]
fn z_roundtrip_qasm() {
    // z gate should survive a round trip
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
z q[0];
cx q[0],q[1];
",
    );
    assert!(
        gates.iter().any(|g| g == "z q[0];"),
        "z gate should be in output, got: {gates:?}"
    );
    assert!(
        gates.iter().any(|g| g.starts_with("cx ")),
        "cx should be in output, got: {gates:?}"
    );
}

#[test]
fn sdg_roundtrip_qasm() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
sdg q[0];
cx q[0],q[1];
",
    );
    assert!(
        gates.iter().any(|g| g == "sdg q[0];"),
        "sdg gate should be in output, got: {gates:?}"
    );
}

#[test]
fn no_rz_for_known_angles() {
    // After optimization, known angles should never produce raw rz
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
",
    );
    // 16T = 2*2π = identity
    assert!(gates.is_empty(), "16T should cancel, got: {gates:?}");
}

#[test]
fn mixed_z_sdg_pipeline() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
z q[0];
sdg q[1];
cx q[0],q[1];
s q[0];
t q[1];
",
    );
    // z + s = 3π/2 = sdg on q0
    // sdg + t = -π/4 = tdg on q1
    // Should not contain any rz gates
    for g in &gates {
        assert!(
            !g.starts_with("rz("),
            "should not have raw rz, got: {gates:?}"
        );
    }
}

#[test]
fn output_idempotent_with_z_sdg() {
    // First pass output with z/sdg should be stable on second pass
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    let pass1 = dir.path().join("p1.qasm");
    let pass2 = dir.path().join("p2.qasm");
    fs::write(
        &input,
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[3];
z q[0];
sdg q[1];
s q[2];
cx q[0],q[1];
cx q[1],q[2];
t q[0];
tdg q[2];
",
    )
    .unwrap();

    let out1 = tzap_run(&[input.to_str().unwrap(), pass1.to_str().unwrap()]);
    assert!(out1.status.success());
    let out2 = tzap_run(&[pass1.to_str().unwrap(), pass2.to_str().unwrap()]);
    assert!(out2.status.success());

    let c1 = fs::read_to_string(&pass1).unwrap();
    let c2 = fs::read_to_string(&pass2).unwrap();
    assert_eq!(c1, c2, "output should be idempotent with z/sdg gates");
}

#[test]
fn native_cz_is_preserved_by_default() {
    let (gates, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
cz q[0],q[1];
",
    );
    assert_eq!(gates, vec!["cz q[0],q[1];"]);
}

const CZ_CHAIN_QASM: &str = "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[3];
cz q[0],q[1];
cz q[1],q[2];
";

#[test]
fn decompose_cz_flag_decomposes_the_default_pipeline_input() {
    let (gates, _) = run_qasm_with_args(CZ_CHAIN_QASM, &["--decompose-cz"]);

    assert!(!gates.iter().any(|gate| gate.starts_with("cz ")));
    assert_eq!(
        gates.iter().filter(|gate| gate.starts_with("cx ")).count(),
        2
    );
    // Both CZs decompose sandwiching q[1] in H...H; since the default
    // pipeline (-O3) further cancels gates, the two separate H-CX-H
    // sandwiches collapse to a single shared H...H around both CXs.
    assert_eq!(
        gates.iter().filter(|gate| gate.starts_with("h ")).count(),
        2
    );
}

#[test]
fn decompose_cz_flag_is_a_noop_when_no_cz_is_present() {
    let (without_flag, _) = run_qasm(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
h q[0];
cx q[0],q[1];
t q[1];
",
    );
    let (with_flag, _) = run_qasm_with_args(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
h q[0];
cx q[0],q[1];
t q[1];
",
        &["--decompose-cz"],
    );

    assert_eq!(with_flag, without_flag);
}

#[test]
fn phase_fold_through_native_cz_on_second_operand() {
    let (gates, _) = run_qasm_with_args(
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
t q[1];
cz q[0],q[1];
t q[1];
",
        &["-O1"],
    );
    assert_eq!(
        gates.len(),
        2,
        "expected native CZ plus folded S: {gates:?}"
    );
    assert_eq!(gates[0], "cz q[0],q[1];");
    assert_eq!(gates[1], "s q[1];");
}

#[test]
fn explicit_decompose_cz_pass_emits_h_cx_h() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("cz.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(
        &input,
        "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
cz q[1],q[0];
",
    )
    .unwrap();

    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--passes",
        "DecomposeCz",
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let content = fs::read_to_string(output).unwrap();
    assert!(!content.lines().any(|line| line.starts_with("cz ")));
    let gates: Vec<_> = content
        .lines()
        .filter(|line| line.starts_with("h ") || line.starts_with("cx "))
        .collect();
    assert_eq!(gates, vec!["h q[0];", "cx q[1],q[0];", "h q[0];"]);
}
