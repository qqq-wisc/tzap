use std::{fs, process::Command};
use tzap::{circuit::Gate, qasm};

#[test]
fn cli_optimizes_across_custom_call_boundaries_and_emits_flat_qasm() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.qasm");
    let output = dir.path().join("output.qasm");
    fs::write(
        &input,
        "OPENQASM 2.0; gate pulse a { t a; } gate inverse a { tdg a; }
        qreg q[2]; pulse q; inverse q;",
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_tzap"))
        .arg(&input)
        .args(["-O1", "-q", "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let flat = fs::read_to_string(output).unwrap();
    assert!(!flat.contains("gate "));
    assert!(qasm::parse(&flat).unwrap().gates().is_empty());
}

#[test]
fn malformed_definition_does_not_overwrite_requested_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.qasm");
    let output = dir.path().join("output.qasm");
    fs::write(
        &input,
        "gate g(theta) a { rz(1/theta) a; } qreg q[1]; g(0) q[0];",
    )
    .unwrap();
    fs::write(&output, "existing output").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_tzap"))
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("division by zero"));
    assert_eq!(fs::read_to_string(output).unwrap(), "existing output");
}

#[test]
fn configurable_limits_apply_to_expansion_but_not_flat_circuits() {
    let limits = qasm::ParseLimits {
        expanded_gates: 1,
        ..Default::default()
    };
    assert!(qasm::parse_with_limits("gate g a { x a; x a; } qreg q[1]; g q[0];", limits).is_err());
    let flat = qasm::parse_with_limits("qreg q[1]; x q[0]; x q[0];", limits).unwrap();
    assert_eq!(flat.gates(), vec![Gate::x(0), Gate::x(0)]);
}
