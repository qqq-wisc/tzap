//! Every erroneous input tzap can be handed — malformed QASM, bad flags,
//! missing files, conflicting options — must be caught with a clear,
//! actionable error message on stderr and a non-zero exit. None of it may
//! ever surface as a raw Rust panic/backtrace: that's an unclear "error
//! message" by definition, no matter how descriptive the panic text is.

use std::fs;
use std::path::Path;
use std::process::{Command, Output};

fn tzap() -> Command {
    Command::new(env!("CARGO_BIN_EXE_tzap"))
}

fn tzap_run(args: &[&str]) -> Output {
    tzap().args(args).output().expect("failed to run tzap")
}

fn write_qasm(dir: &Path, name: &str, contents: &str) -> std::path::PathBuf {
    let path = dir.join(name);
    fs::write(&path, contents).unwrap();
    path
}

/// Assert `out` is how tzap must fail: non-zero exit, nothing on stdout,
/// and a stderr message that unambiguously reads as an error (never a raw
/// panic/backtrace). Returns stderr for further content checks.
fn assert_clear_error(out: &Output, context: &str) -> String {
    assert!(
        !out.status.success(),
        "{context}: expected tzap to fail, but it exited successfully"
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.trim().is_empty(),
        "{context}: expected no stdout on error, got:\n{stdout}"
    );
    let stderr = String::from_utf8_lossy(&out.stderr).to_string();
    let lower = stderr.to_lowercase();
    assert!(
        !lower.contains("panicked"),
        "{context}: tzap panicked instead of returning a clean error:\n{stderr}"
    );
    assert!(
        !lower.contains("rust_backtrace"),
        "{context}: tzap surfaced a raw backtrace hint instead of a clean error:\n{stderr}"
    );
    assert!(
        lower.contains("error"),
        "{context}: message doesn't clearly read as an error:\n{stderr}"
    );
    stderr
}

// ---------------------------------------------------------------------
// Missing / unreadable files
// ---------------------------------------------------------------------

#[test]
fn missing_input_file() {
    let out = tzap_run(&["does-not-exist.qasm"]);
    let stderr = assert_clear_error(&out, "missing input file");
    assert!(stderr.contains("Error reading"), "got: {stderr}");
    assert!(stderr.contains("does-not-exist.qasm"), "got: {stderr}");
}

#[test]
fn directory_as_input_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = tzap_run(&[dir.path().to_str().unwrap()]);
    let stderr = assert_clear_error(&out, "directory as input");
    assert!(stderr.contains("Error reading"), "got: {stderr}");
}

#[test]
fn directory_as_output_path() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out_dir = dir.path().join("subdir");
    fs::create_dir(&out_dir).unwrap();

    let out = tzap_run(&[input.to_str().unwrap(), "-o", out_dir.to_str().unwrap()]);
    let stderr = assert_clear_error(&out, "directory as output path");
    assert!(stderr.contains("Error writing"), "got: {stderr}");
}

#[test]
fn output_path_in_nonexistent_directory() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let bad_out = dir.path().join("no-such-subdir").join("out.qasm");

    let out = tzap_run(&[input.to_str().unwrap(), "-o", bad_out.to_str().unwrap()]);
    let stderr = assert_clear_error(&out, "output path in missing directory");
    assert!(stderr.contains("Error writing"), "got: {stderr}");
}

const TRIVIAL_QASM: &str = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0];\n";

// ---------------------------------------------------------------------
// Malformed QASM: every error path in src/qasm.rs
// ---------------------------------------------------------------------

fn qasm_error(dir: &Path, name: &str, contents: &str) -> String {
    let path = write_qasm(dir, name, contents);
    let out = tzap_run(&[path.to_str().unwrap()]);
    assert_clear_error(&out, name)
}

#[test]
fn qreg_after_gate_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "qreg_after_gate.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0];\nqreg r[1];\n",
    );
    assert!(
        stderr.contains("qreg declaration after gate"),
        "got: {stderr}"
    );
}

#[test]
fn creg_after_gate_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "creg_after_gate.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0];\ncreg c[1];\n",
    );
    assert!(
        stderr.contains("creg declaration after gate"),
        "got: {stderr}"
    );
}

#[test]
fn bad_qreg_size_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "bad_qreg.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[abc];\n",
    );
    assert!(stderr.contains("bad qreg size"), "got: {stderr}");
}

#[test]
fn bad_creg_size_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "bad_creg.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[abc];\n",
    );
    assert!(stderr.contains("bad creg size"), "got: {stderr}");
}

#[test]
fn unsupported_gate_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "unsupported.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ny q[0];\n",
    );
    assert!(stderr.contains("unsupported: y"), "got: {stderr}");
}

#[test]
fn garbage_input_is_rejected_not_panicked() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "garbage.qasm",
        "this is not qasm at all\n{}[]<>\n",
    );
    assert!(stderr.contains("Error parsing"), "got: {stderr}");
}

#[test]
fn empty_angle_expression_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "empty_angle.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz() q[0];\n",
    );
    assert!(stderr.contains("empty angle expression"), "got: {stderr}");
}

#[test]
fn unexpected_token_in_angle_expression_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "bad_angle.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi pi) q[0];\n",
    );
    assert!(stderr.contains("unexpected token"), "got: {stderr}");
}

#[test]
fn unexpected_character_in_angle_expression_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "bad_char.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi & 2) q[0];\n",
    );
    assert!(stderr.contains("unexpected character"), "got: {stderr}");
}

#[test]
fn dangling_operator_in_angle_expression_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "dangling_op.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi+) q[0];\n",
    );
    assert!(
        stderr.contains("unexpected end of angle expression"),
        "got: {stderr}"
    );
}

#[test]
fn unclosed_parenthesis_in_angle_expression_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "unclosed_paren_expr.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz((pi+1) q[0];\n",
    );
    // The outer rz(...) itself never finds its matching close paren either,
    // since the extra '(' shifts what find_matching_paren treats as depth 0.
    assert!(stderr.contains("rz missing closing"), "got: {stderr}");
}

#[test]
fn rz_missing_closing_paren_is_rejected_not_dropped() {
    // Regression guard: this used to silently drop the gate (no error, no
    // gate emitted) instead of reporting the malformed line.
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "rz_unclosed.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi/5 q[0];\n",
    );
    assert!(stderr.contains("rz missing closing ')'"), "got: {stderr}");
}

#[test]
fn measure_missing_arrow_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "measure_no_arrow.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[1];\nmeasure q[0] c[0];\n",
    );
    assert!(stderr.contains("measure missing '->'"), "got: {stderr}");
}

#[test]
fn measure_operand_size_mismatch_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "measure_mismatch.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncreg c[1];\nmeasure q -> c;\n",
    );
    assert!(
        stderr.contains("measure operand size mismatch"),
        "got: {stderr}"
    );
}

#[test]
fn unknown_register_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "unknown_reg.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh r[0];\n",
    );
    assert!(stderr.contains("unknown register 'r'"), "got: {stderr}");
}

#[test]
fn unknown_classical_register_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "unknown_creg.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[1];\nmeasure q[0] -> d[0];\n",
    );
    assert!(
        stderr.contains("unknown classical register 'd'"),
        "got: {stderr}"
    );
}

#[test]
fn bad_qubit_index_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "bad_qubit_index.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nh q[x];\n",
    );
    assert!(stderr.contains("bad qubit index"), "got: {stderr}");
}

#[test]
fn bad_cbit_index_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "bad_cbit_index.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[1];\nmeasure q[0] -> c[x];\n",
    );
    assert!(stderr.contains("bad cbit index"), "got: {stderr}");
}

#[test]
fn qubit_index_out_of_range_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "qubit_oob.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nh q[5];\n",
    );
    assert!(
        stderr.contains("index 5 out of range for register 'q' (size 2)"),
        "got: {stderr}"
    );
}

#[test]
fn cbit_index_out_of_range_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "cbit_oob.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[1];\nmeasure q[0] -> c[5];\n",
    );
    assert!(
        stderr.contains("index 5 out of range for classical register 'c' (size 1)"),
        "got: {stderr}"
    );
}

/// Regression guard: `cx`/`ccx` used to index straight into the resolved
/// qubit list (`qubits[0]`, `qubits[1]`, ...) with no arity check, so the
/// wrong operand count panicked with a raw index-out-of-bounds instead of
/// a parse error.
#[test]
fn cx_wrong_arity_is_rejected_not_panicked() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "cx_arity.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncx q[0];\n",
    );
    assert!(
        stderr.contains("cx expects 2 qubit operands, got 1"),
        "got: {stderr}"
    );
}

#[test]
fn cx_too_many_operands_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "cx_arity3.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\ncx q[0],q[1],q[2];\n",
    );
    assert!(
        stderr.contains("cx expects 2 qubit operands, got 3"),
        "got: {stderr}"
    );
}

#[test]
fn ccx_wrong_arity_is_rejected_not_panicked() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "ccx_arity.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nccx q[0],q[1];\n",
    );
    assert!(
        stderr.contains("ccx expects 3 qubit operands, got 2"),
        "got: {stderr}"
    );
}

#[test]
fn cz_wrong_arity_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "cz_arity.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\ncz q[0],q[1],q[2];\n",
    );
    assert!(
        stderr.contains("cz expects 2 qubit operands, got 3"),
        "got: {stderr}"
    );
}

#[test]
fn ccz_wrong_arity_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "ccz_arity.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nccz q[0],q[1];\n",
    );
    assert!(
        stderr.contains("ccz expects 3 qubit operands, got 2"),
        "got: {stderr}"
    );
}

/// Regression guard: single-qubit gates (`h`, `x`, `s`, `sdg`, `z`, `t`,
/// `tdg`, `rz`) used to index `resolve_qubits(..)[0]` directly, panicking on
/// an empty or multi-qubit operand instead of returning a parse error.
#[test]
fn single_qubit_gates_reject_wrong_arity_not_panic() {
    let dir = tempfile::tempdir().unwrap();
    for gate in ["h", "x", "s", "sdg", "z", "t", "tdg"] {
        let stderr = qasm_error(
            dir.path(),
            &format!("{gate}_two_operands.qasm"),
            &format!("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\n{gate} q[0],q[1];\n"),
        );
        assert!(
            stderr.contains(&format!("{gate} expects 1 qubit operand, got 2")),
            "gate {gate} got: {stderr}"
        );
    }
}

#[test]
fn rz_wrong_arity_is_rejected_not_panicked() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "rz_arity.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nrz(pi/5) q[0],q[1];\n",
    );
    assert!(
        stderr.contains("rz expects 1 qubit operand, got 2"),
        "got: {stderr}"
    );
}

#[test]
fn negative_index_qubit_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let stderr = qasm_error(
        dir.path(),
        "negative_index.qasm",
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\nh q[-1];\n",
    );
    assert!(stderr.contains("bad qubit index"), "got: {stderr}");
}

// ---------------------------------------------------------------------
// CLI argument errors
// ---------------------------------------------------------------------

#[test]
fn no_args_is_rejected() {
    let out = tzap_run(&[]);
    let stderr = assert_clear_error(&out, "no args");
    assert!(
        stderr.contains("missing required <input.qasm> argument"),
        "got: {stderr}"
    );
    assert!(
        stderr.contains("--help"),
        "should hint at --help, got: {stderr}"
    );
}

#[test]
fn unknown_flag_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "--not-a-real-flag"]);
    let stderr = assert_clear_error(&out, "unknown flag");
    assert!(
        stderr.contains("unknown flag '--not-a-real-flag'"),
        "got: {stderr}"
    );
    assert!(
        stderr.contains("--help"),
        "should hint at --help, got: {stderr}"
    );
}

#[test]
fn dash_o_with_no_value_is_rejected() {
    // Regression guard: this used to silently fall back to "no output file"
    // instead of erroring.
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "-o"]);
    let stderr = assert_clear_error(&out, "-o with no value");
    assert!(
        stderr.contains("-o requires an output file path"),
        "got: {stderr}"
    );
}

#[test]
fn extra_positional_argument_is_rejected() {
    // Regression guard: a third positional argument used to be silently
    // dropped instead of erroring.
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "out.qasm", "unexpected-extra.qasm"]);
    let stderr = assert_clear_error(&out, "extra positional argument");
    assert!(
        stderr.contains("unexpected extra argument 'unexpected-extra.qasm'"),
        "got: {stderr}"
    );
}

#[test]
fn epsilon_with_no_value_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "--epsilon"]);
    let stderr = assert_clear_error(&out, "--epsilon with no value");
    assert!(
        stderr.contains("--epsilon requires a number"),
        "got: {stderr}"
    );
}

#[test]
fn epsilon_non_number_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "--epsilon", "banana"]);
    let stderr = assert_clear_error(&out, "--epsilon non-number");
    assert!(
        stderr.contains("--epsilon requires a number"),
        "got: {stderr}"
    );
}

/// Regression guard: `--epsilon 0` used to reach rsgridsynth and panic
/// there with a raw dashu-int backtrace ("logarithm is not defined for 0").
#[test]
fn epsilon_zero_is_rejected_not_panicked() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(
        dir.path(),
        "rz.qasm",
        "\
OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi/5) q[0];\n",
    );
    let out = tzap_run(&[input.to_str().unwrap(), "--decompose-rz", "--epsilon", "0"]);
    let stderr = assert_clear_error(&out, "--epsilon 0");
    assert!(
        stderr.contains("--epsilon must be a positive"),
        "got: {stderr}"
    );
}

/// Regression guard: negative epsilon used to be silently accepted and fed
/// straight to gridsynth instead of being rejected as meaningless.
#[test]
fn epsilon_negative_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(
        dir.path(),
        "rz.qasm",
        "\
OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi/5) q[0];\n",
    );
    let out = tzap_run(&[
        input.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "-1e-3",
    ]);
    let stderr = assert_clear_error(&out, "--epsilon negative");
    assert!(
        stderr.contains("--epsilon must be a positive"),
        "got: {stderr}"
    );
}

#[test]
fn epsilon_nan_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(
        dir.path(),
        "rz.qasm",
        "\
OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi/5) q[0];\n",
    );
    let out = tzap_run(&[
        input.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "NaN",
    ]);
    let stderr = assert_clear_error(&out, "--epsilon NaN");
    assert!(
        stderr.contains("--epsilon must be a positive"),
        "got: {stderr}"
    );
}

#[test]
fn epsilon_infinity_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(
        dir.path(),
        "rz.qasm",
        "\
OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi/5) q[0];\n",
    );
    let out = tzap_run(&[
        input.to_str().unwrap(),
        "--decompose-rz",
        "--epsilon",
        "inf",
    ]);
    let stderr = assert_clear_error(&out, "--epsilon inf");
    assert!(
        stderr.contains("--epsilon must be a positive"),
        "got: {stderr}"
    );
}

#[test]
fn passes_with_no_value_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "--passes"]);
    let stderr = assert_clear_error(&out, "--passes with no value");
    assert!(
        stderr.contains("--passes requires a comma-separated list"),
        "got: {stderr}"
    );
}

#[test]
fn passes_empty_list_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "--passes", ",,,"]);
    let stderr = assert_clear_error(&out, "--passes empty list");
    assert!(
        stderr.contains("--passes requires at least one pass name"),
        "got: {stderr}"
    );
}

#[test]
fn passes_unknown_name_lists_valid_names() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "--passes", "NotARealPass"]);
    let stderr = assert_clear_error(&out, "--passes unknown name");
    assert!(
        stderr.contains("Unknown pass 'NotARealPass'"),
        "got: {stderr}"
    );
    assert!(stderr.contains("Available passes:"), "got: {stderr}");
    assert!(stderr.contains("CancelGates"), "got: {stderr}");
}

#[test]
fn optimization_level_conflicts_with_passes() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "-O2", "--passes", "CancelGates"]);
    let stderr = assert_clear_error(&out, "-O2 with --passes");
    assert!(
        stderr.contains("cannot be combined with --passes or --fixpoint"),
        "got: {stderr}"
    );
}

#[test]
fn optimization_levels_mutually_exclusive_all_pairs() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let levels = ["-O1", "-O2", "-O3", "-Osuper"];
    for a in levels {
        for b in levels {
            if a == b {
                continue;
            }
            let out = tzap_run(&[input.to_str().unwrap(), a, b]);
            let stderr = assert_clear_error(&out, &format!("{a} + {b}"));
            assert!(
                stderr.contains("cannot be combined"),
                "{a}+{b} got: {stderr}"
            );
        }
    }
}

#[test]
fn passes_conflicts_with_decomposition_flags() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    for flag in ["--decompose-rz", "--decompose-cz", "--decompose-ccx"] {
        let out = tzap_run(&[input.to_str().unwrap(), "--passes", "CancelGates", flag]);
        let stderr = assert_clear_error(&out, &format!("--passes + {flag}"));
        assert!(stderr.contains("cannot be combined"), "got: {stderr}");
    }
}

#[test]
fn superopt_gates_rejects_missing_empty_and_unsupported_values() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let input = input.to_str().unwrap();
    for (args, expected) in [
        (vec![input, "--superopt-gates"], "--superopt-gates requires"),
        (
            vec![input, "--superopt-gates="],
            "--superopt-gates requires",
        ),
        (
            vec![input, "--superopt-gates", "rz"],
            "cannot be emitted by SuperOpt",
        ),
        (
            vec![input, "--superopt-gates", "h,nope"],
            "unsupported SuperOpt gate",
        ),
    ] {
        let out = tzap_run(&args);
        let stderr = assert_clear_error(&out, &format!("{args:?}"));
        assert!(stderr.contains(expected), "got: {stderr}");
    }
}

// --- hidden --superopt-* flags ---

#[test]
fn superopt_qubits_non_integer_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-O2",
        "--superopt-qubits",
        "banana",
    ]);
    let stderr = assert_clear_error(&out, "--superopt-qubits non-integer");
    assert!(
        stderr.contains("--superopt-qubits requires a positive integer"),
        "got: {stderr}"
    );
}

#[test]
fn superopt_qubits_with_no_value_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "-O2", "--superopt-qubits"]);
    let stderr = assert_clear_error(&out, "--superopt-qubits with no value");
    assert!(
        stderr.contains("--superopt-qubits requires a positive integer"),
        "got: {stderr}"
    );
}

#[test]
fn superopt_qubits_negative_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "-O2", "--superopt-qubits", "-1"]);
    let stderr = assert_clear_error(&out, "--superopt-qubits negative");
    assert!(
        stderr.contains("--superopt-qubits requires a positive integer"),
        "got: {stderr}"
    );
}

/// Regression guard: `--superopt-window-gates 0` used to pass CLI parsing
/// (0 is a valid `usize`) and only panic later, mid-run, inside the
/// `Pass::run` trait impl ("SuperOpt failed: window_gates must be greater
/// than zero").
#[test]
fn superopt_window_gates_zero_is_rejected_not_panicked() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-O2",
        "--superopt-window-gates",
        "0",
    ]);
    let stderr = assert_clear_error(&out, "--superopt-window-gates 0");
    assert!(
        stderr.contains("--superopt-window-gates requires a positive integer, got 0"),
        "got: {stderr}"
    );
}

#[test]
fn superopt_qubits_zero_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[input.to_str().unwrap(), "-O2", "--superopt-qubits", "0"]);
    let stderr = assert_clear_error(&out, "--superopt-qubits 0");
    assert!(
        stderr.contains("--superopt-qubits requires a positive integer, got 0"),
        "got: {stderr}"
    );
}

#[test]
fn superopt_murm_entries_zero_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let out = tzap_run(&[
        input.to_str().unwrap(),
        "-O2",
        "--superopt-murm-entries",
        "0",
    ]);
    let stderr = assert_clear_error(&out, "--superopt-murm-entries 0");
    assert!(
        stderr.contains("--superopt-murm-entries requires a positive integer, got 0"),
        "got: {stderr}"
    );
}

// ---------------------------------------------------------------------
// Table-driven sweep: a large batch of erroneous invocations, checked in
// bulk against the same "reads as a clear error, never panics" bar. This
// is deliberately redundant with the targeted tests above — the point is
// breadth, so a future change that reintroduces a raw panic on any one of
// these has nowhere to hide.
// ---------------------------------------------------------------------

#[test]
fn bulk_sweep_of_malformed_qasm_never_panics() {
    let dir = tempfile::tempdir().unwrap();
    let cases: &[(&str, &str)] = &[
        ("gate_before_qreg", "OPENQASM 2.0;\nh q[0];\n"),
        (
            "reset_unknown_register",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nreset r[0];\n",
        ),
        (
            "reset_bad_index",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nreset q[9];\n",
        ),
        (
            "ccz_wrong_arity_one",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nccz q[0];\n",
        ),
        (
            "ccx_wrong_arity_four",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[4];\nccx q[0],q[1],q[2],q[3];\n",
        ),
        (
            "cz_self_and_unknown_mix",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncz q[0],r[1];\n",
        ),
        (
            "rz_bad_index",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi/4) q[9];\n",
        ),
        (
            "rz_double_open_paren",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz((pi/4) q[0];\n",
        ),
        (
            "angle_unary_minus_only",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(-) q[0];\n",
        ),
        (
            "angle_trailing_operator",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nrz(pi*) q[0];\n",
        ),
        (
            "measure_bad_qubit_register",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[1];\nmeasure r[0] -> c[0];\n",
        ),
        (
            "measure_bad_cbit_register",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\ncreg c[1];\nmeasure q[0] -> d[0];\n",
        ),
        (
            "h_no_operand",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh;\n",
        ),
        (
            "cx_no_operands",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncx;\n",
        ),
        (
            "unknown_word_as_statement",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nfrobnicate q[0];\n",
        ),
    ];

    for (name, qasm) in cases {
        let path = write_qasm(dir.path(), &format!("{name}.qasm"), qasm);
        let out = tzap_run(&[path.to_str().unwrap(), "-O1"]);
        assert_clear_error(&out, name);
    }
}

#[test]
fn bulk_sweep_of_bad_cli_invocations_never_panics() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_qasm(dir.path(), "in.qasm", TRIVIAL_QASM);
    let input = input.to_str().unwrap();

    let cases: Vec<Vec<&str>> = vec![
        vec![],
        vec!["--help-me"],
        vec![input, "-x"],
        vec![input, "--"],
        vec![input, "-O5"],
        vec![input, "--passes"],
        vec![input, "--passes", ""],
        vec![input, "--passes", "cancelgates"], // wrong case
        vec![input, "--epsilon"],
        vec![input, "--epsilon", ""],
        vec![input, "--epsilon", "1e"],
        vec![input, "-o"],
        vec!["-O1", "-O2"],
        vec![input, input, input],
        vec!["--superopt-qubits", "1"],
    ];

    for args in cases {
        let out = tzap_run(&args);
        assert_clear_error(&out, &format!("{args:?}"));
    }
}

/// Degenerate inputs that are *not* errors — an empty program, a header
/// with no gates, a statement with no trailing `;` on its own line, an
/// empty statement from `;;` — must still succeed cleanly rather than
/// being rejected or (worse) crashing.
#[test]
fn degenerate_but_valid_qasm_succeeds() {
    let dir = tempfile::tempdir().unwrap();
    let cases: &[(&str, &str)] = &[
        ("empty_program", ""),
        ("only_header", "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n"),
        (
            "unterminated_line",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0]",
        ),
        (
            "double_semicolon",
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];;\nh q[0];\n",
        ),
    ];
    for (name, qasm) in cases {
        let path = write_qasm(dir.path(), &format!("{name}.qasm"), qasm);
        let out = tzap_run(&[path.to_str().unwrap()]);
        assert!(
            out.status.success(),
            "{name} should succeed, stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
}
