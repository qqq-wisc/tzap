use std::fs;
use std::process::Command;

fn tzap() -> Command {
    Command::new(env!("CARGO_BIN_EXE_tzap"))
}

fn tzap_run(args: &[&str]) -> std::process::Output {
    tzap().args(args).output().expect("failed to run tzap")
}

#[test]
fn no_args_prints_usage() {
    let out = tzap_run(&[]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("tzap"), "expected usage message, got: {stderr}");
}

#[test]
fn missing_file_errors() {
    let out = tzap_run(&["nonexistent.qasm"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("Error reading"), "expected error message, got: {stderr}");
}

#[test]
fn optimizes_test_qasm_to_file() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&["qasm/test.qasm", "-o", out_path.to_str().unwrap()]);
    assert!(out.status.success(), "tzap failed: {}", String::from_utf8_lossy(&out.stderr));

    let content = fs::read_to_string(&out_path).unwrap();
    assert!(content.starts_with("OPENQASM 2.0;"));
    assert!(content.contains("qreg q[3];"));

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("test.qasm"));
    assert!(stderr.contains("gates"));
    assert!(stderr.contains("T/Tdg"));
}

#[test]
fn output_is_valid_qasm() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&["qasm/test.qasm", "-o", out_path.to_str().unwrap()]);
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

    let out = tzap_run(&["qasm/test.qasm", out_path.to_str().unwrap()]);
    assert!(out.status.success(), "tzap failed: {}", String::from_utf8_lossy(&out.stderr));

    // stdout should be empty when writing to file
    assert!(out.stdout.is_empty(), "stdout should be empty when output file given");

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

    let out1 = tzap_run(&["qasm/test.qasm", pass1.to_str().unwrap()]);
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

    let out = tzap_run(&["qasm/two_ccx.qasm", "-o", out_path.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&out_path).unwrap();
    let gates = gate_lines_from(&content);
    assert!(gates.len() > 3, "decomposed circuit should have more than 3 gates, got {}", gates.len());
    // No CCX gates should remain
    assert!(!gates.iter().any(|g| g.starts_with("ccx ")),
        "no ccx gates should remain after decomposition");
}

#[test]
fn mod5_4_reduces_t_count() {
    let dir = tempfile::tempdir().unwrap();
    let out_path = dir.path().join("out.qasm");

    let out = tzap_run(&["qasm/mod5_4.qasm", "-o", out_path.to_str().unwrap(), "--cancel"]);
    assert!(out.status.success());

    let content = fs::read_to_string(&out_path).unwrap();
    let gates = gate_lines_from(&content);

    let t_count = gates.iter()
        .filter(|g| g.starts_with("t ") || g.starts_with("tdg "))
        .count();

    assert_eq!(t_count, 16, "mod5_4 should optimize to 16 T/Tdg, got {t_count}");
    assert_eq!(gates.len(), 57, "mod5_4 should optimize to 57 gates, got {}", gates.len());
}

#[test]
fn mod5_4_idempotent() {
    let dir = tempfile::tempdir().unwrap();
    let pass1 = dir.path().join("pass1.qasm");
    let pass2 = dir.path().join("pass2.qasm");

    let out1 = tzap_run(&["qasm/mod5_4.qasm", pass1.to_str().unwrap()]);
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
    let out = tzap_run(&["qasm/mod5_4.qasm", "-o", out_path.to_str().unwrap()]);
    assert!(out.status.success());
    let content = fs::read_to_string(&out_path).unwrap();
    for line in content.lines() {
        assert!(!line.starts_with("rz("),
            "mod5_4 output should not contain raw rz gates, found: {line}");
    }
}

#[test]
fn inline_qasm_optimization() {
    // Write a known circuit to a temp file and verify T+T folds to S
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("tt.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
").unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "-o", output.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&output).unwrap();
    let gate_lines: Vec<&str> = content.lines()
        .filter(|l| !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg"))
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(gate_lines.len(), 1, "T+T should fold to S, got: {gate_lines:?}");
    assert!(gate_lines[0].starts_with("s "), "T+T should fold to S, got: {gate_lines:?}");
}

#[test]
fn t_t_folds_to_s() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("tt.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
").unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "-o", output.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&output).unwrap();
    let gate_lines: Vec<&str> = content.lines()
        .filter(|l| !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg"))
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(gate_lines.len(), 1, "T+T should fold to single gate, got: {gate_lines:?}");
    assert_eq!(gate_lines[0], "s q[0];", "T+T should fold to S");
}

#[test]
fn t_tdg_cancels() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("t_tdg.qasm");
    let output = dir.path().join("out.qasm");
    fs::write(&input, "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
tdg q[0];
").unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "-o", output.to_str().unwrap()]);
    assert!(out.status.success());

    let content = fs::read_to_string(&output).unwrap();
    let gate_lines: Vec<&str> = content.lines()
        .filter(|l| !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg"))
        .filter(|l| !l.is_empty())
        .collect();
    assert!(gate_lines.is_empty(), "T+Tdg should cancel, got: {gate_lines:?}");
}

#[test]
fn phase_fold_local_through_cnot() {
    // T q0; CNOT q0,q1; T q0 should merge to CNOT + S
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
t q[0];
cx q[0],q[1];
t q[0];
");
    assert_eq!(gates.len(), 2, "should merge T through CNOT control: {gates:?}");
    assert!(gates.iter().any(|g| g.starts_with("cx ")), "CNOT should remain");
    assert!(gates.iter().any(|g| g == "s q[0];"), "T+T should become S");
}

#[test]
fn phase_fold_local_cancel_through_cnot() {
    // S q0; CNOT q0,q1; Sdg q0 — should cancel rotations
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
s q[0];
cx q[0],q[1];
sdg q[0];
");
    assert_eq!(gates.len(), 1, "S and Sdg should cancel through CNOT: {gates:?}");
    assert!(gates[0].starts_with("cx "), "only CNOT should remain");
}

#[test]
fn phase_fold_local_blocked_on_target() {
    // T q[1]; CNOT q[0],q[1]; T q[1] — q1 is target, can't merge
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
t q[1];
cx q[0],q[1];
t q[1];
");
    // Phase fold may or may not touch this, but the two T's should not merge
    assert!(gates.len() >= 3, "should not merge T on CNOT target: {gates:?}");
}

fn gate_lines_from(stdout: &str) -> Vec<String> {
    stdout.lines()
        .filter(|l| !l.starts_with("OPENQASM") && !l.starts_with("include") && !l.starts_with("qreg"))
        .filter(|l| !l.is_empty())
        .map(|l| l.to_string())
        .collect()
}

fn run_qasm(qasm: &str) -> (Vec<String>, String) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.qasm");
    let output = dir.path().join("output.qasm");
    fs::write(&input, qasm).unwrap();
    let out = tzap_run(&[input.to_str().unwrap(), "-o", output.to_str().unwrap()]);
    assert!(out.status.success(), "tzap failed: {}", String::from_utf8_lossy(&out.stderr));
    let content = fs::read_to_string(&output).unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).to_string();
    (gate_lines_from(&content), stderr)
}

#[test]
fn s_sdg_cancels() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
sdg q[0];
");
    assert!(gates.is_empty(), "S+Sdg should cancel, got: {gates:?}");
}

#[test]
fn sdg_s_cancels() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
sdg q[0];
s q[0];
");
    assert!(gates.is_empty(), "Sdg+S should cancel, got: {gates:?}");
}

#[test]
fn z_z_cancels() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
z q[0];
z q[0];
");
    assert!(gates.is_empty(), "Z+Z should cancel, got: {gates:?}");
}

#[test]
fn s_s_becomes_z() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
s q[0];
");
    assert_eq!(gates.len(), 1, "S+S should fold to one gate, got: {gates:?}");
    assert_eq!(gates[0], "z q[0];");
}

#[test]
fn sdg_sdg_becomes_z() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
sdg q[0];
sdg q[0];
");
    assert_eq!(gates.len(), 1, "Sdg+Sdg should fold to one gate, got: {gates:?}");
    assert_eq!(gates[0], "z q[0];");
}

#[test]
fn z_s_becomes_sdg() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
z q[0];
s q[0];
");
    assert_eq!(gates.len(), 1, "Z+S should fold to Sdg, got: {gates:?}");
    assert_eq!(gates[0], "sdg q[0];");
}

#[test]
fn z_sdg_becomes_s() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
z q[0];
sdg q[0];
");
    assert_eq!(gates.len(), 1, "Z+Sdg should fold to S, got: {gates:?}");
    assert_eq!(gates[0], "s q[0];");
}

#[test]
fn sdg_t_becomes_tdg() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
sdg q[0];
t q[0];
");
    assert_eq!(gates.len(), 1, "Sdg+T should fold to Tdg, got: {gates:?}");
    assert_eq!(gates[0], "tdg q[0];");
}

#[test]
fn s_tdg_becomes_t() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
tdg q[0];
");
    assert_eq!(gates.len(), 1, "S+Tdg should fold to T, got: {gates:?}");
    assert_eq!(gates[0], "t q[0];");
}

#[test]
fn six_t_becomes_sdg() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
t q[0];
");
    assert_eq!(gates.len(), 1, "6T should fold to Sdg, got: {gates:?}");
    assert_eq!(gates[0], "sdg q[0];");
}

#[test]
fn four_t_becomes_z() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
t q[0];
t q[0];
t q[0];
t q[0];
");
    assert_eq!(gates.len(), 1, "4T should fold to Z, got: {gates:?}");
    assert_eq!(gates[0], "z q[0];");
}

#[test]
fn four_s_cancels() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[1];
s q[0];
s q[0];
s q[0];
s q[0];
");
    assert!(gates.is_empty(), "4S should cancel to identity, got: {gates:?}");
}

#[test]
fn z_roundtrip_qasm() {
    // z gate should survive a round trip
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
z q[0];
cx q[0],q[1];
");
    assert!(gates.iter().any(|g| g == "z q[0];"), "z gate should be in output, got: {gates:?}");
    assert!(gates.iter().any(|g| g.starts_with("cx ")), "cx should be in output, got: {gates:?}");
}

#[test]
fn sdg_roundtrip_qasm() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
sdg q[0];
cx q[0],q[1];
");
    assert!(gates.iter().any(|g| g == "sdg q[0];"), "sdg gate should be in output, got: {gates:?}");
}

#[test]
fn no_rz_for_known_angles() {
    // After optimization, known angles should never produce raw rz
    let (gates, _) = run_qasm("\
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
");
    // 16T = 2*2π = identity
    assert!(gates.is_empty(), "16T should cancel, got: {gates:?}");
}

#[test]
fn mixed_z_sdg_pipeline() {
    let (gates, _) = run_qasm("\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
z q[0];
sdg q[1];
cx q[0],q[1];
s q[0];
t q[1];
");
    // z + s = 3π/2 = sdg on q0
    // sdg + t = -π/4 = tdg on q1
    // Should not contain any rz gates
    for g in &gates {
        assert!(!g.starts_with("rz("), "should not have raw rz, got: {gates:?}");
    }
}

#[test]
fn output_idempotent_with_z_sdg() {
    // First pass output with z/sdg should be stable on second pass
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.qasm");
    let pass1 = dir.path().join("p1.qasm");
    let pass2 = dir.path().join("p2.qasm");
    fs::write(&input, "\
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
").unwrap();

    let out1 = tzap_run(&[input.to_str().unwrap(), pass1.to_str().unwrap()]);
    assert!(out1.status.success());
    let out2 = tzap_run(&[pass1.to_str().unwrap(), pass2.to_str().unwrap()]);
    assert!(out2.status.success());

    let c1 = fs::read_to_string(&pass1).unwrap();
    let c2 = fs::read_to_string(&pass2).unwrap();
    assert_eq!(c1, c2, "output should be idempotent with z/sdg gates");
}




