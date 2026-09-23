import subprocess
import sys
from types import SimpleNamespace

import pytest
from qiskit import qasm2
from tzap import _cli

QASM = """\
OPENQASM 2.0;
include "qelib1.inc";
qreg q[1];
h q[0];
h q[0];
"""


def run_module(*args):
    return subprocess.run(
        [sys.executable, "-m", "tzap", *map(str, args)],
        check=False,
        capture_output=True,
        text=True,
    )


def fake_result(qasm=QASM):
    metrics = SimpleNamespace(gates=0, t=0)
    report = SimpleNamespace(baseline=metrics, output=metrics)
    return SimpleNamespace(qasm=qasm, report=report)


def test_module_version():
    completed = run_module("--version")

    assert completed.returncode == 0
    assert completed.stdout.startswith("tzap ")
    assert completed.stderr == ""


def test_module_help():
    completed = run_module("--help")

    assert completed.returncode == 0
    assert "fast Clifford+T/Rz circuit optimizer" in completed.stdout
    assert "--decompose-rz" in completed.stdout
    assert "--decompose-ccx" in completed.stdout
    assert "--superopt-gates" in completed.stdout


def test_module_optimizes_to_output_flag(tmp_path):
    source = tmp_path / "input.qasm"
    output = tmp_path / "output.qasm"
    source.write_text(QASM, encoding="utf-8")

    completed = run_module(source, "-O1", "-o", output)

    assert completed.returncode == 0
    assert "2 -> 0 gates" in completed.stderr
    assert "wrote" in completed.stderr
    assert output.exists()
    assert len(qasm2.load(output).data) == 0


def test_module_supports_positional_output(tmp_path):
    source = tmp_path / "input.qasm"
    output = tmp_path / "output.qasm"
    source.write_text(QASM, encoding="utf-8")

    completed = run_module(source, output, "-O1")

    assert completed.returncode == 0
    assert output.exists()


def test_module_custom_pass_pipeline(tmp_path):
    source = tmp_path / "input.qasm"
    output = tmp_path / "output.qasm"
    source.write_text(QASM, encoding="utf-8")

    completed = run_module(
        source,
        "-o",
        output,
        "--passes",
        "CancelGates,PhaseFoldRand",
    )

    assert completed.returncode == 0
    assert len(qasm2.load(output).data) == 0


def test_module_missing_input_is_an_argparse_error():
    completed = run_module()

    assert completed.returncode == 2
    assert "required" in completed.stderr


def test_module_missing_file_is_reported_without_traceback(tmp_path):
    completed = run_module(tmp_path / "missing.qasm", "-O1")

    assert completed.returncode == 1
    assert "Error:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_module_bad_qasm_is_reported_without_traceback(tmp_path):
    source = tmp_path / "bad.qasm"
    source.write_text(
        "OPENQASM 2.0; qreg q[1]; y q[0];",
        encoding="utf-8",
    )

    completed = run_module(source, "-O1")

    assert completed.returncode == 1
    assert "unsupported" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_module_rejects_invalid_epsilon_cleanly(tmp_path):
    source = tmp_path / "input.qasm"
    source.write_text(QASM, encoding="utf-8")

    completed = run_module(source, "-O1", "--epsilon", "0")

    assert completed.returncode == 1
    assert "positive, finite" in completed.stderr


def test_module_reports_output_write_failures_without_traceback(tmp_path):
    source = tmp_path / "input.qasm"
    source.write_text(QASM, encoding="utf-8")

    completed = run_module(source, "-O1", "-o", tmp_path)

    assert completed.returncode == 1
    assert "Error:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_main_forwards_every_optimizer_flag(monkeypatch, tmp_path):
    source = tmp_path / "input.qasm"
    source.write_text(QASM, encoding="utf-8")
    captured = {}

    def optimize_qasm(qasm, **options):
        captured["qasm"] = qasm
        captured["options"] = options
        return fake_result()

    monkeypatch.setattr(_cli, "optimize_qasm", optimize_qasm)

    exit_code = _cli.main(
        [
            str(source),
            "-Osuper",
            "--fixpoint",
            "--decompose-rz",
            "--decompose-cz",
            "--decompose-ccx",
            "--superopt-gates",
            "base",
            "--epsilon",
            "0.001",
            "--parallel",
        ]
    )

    assert exit_code == 0
    assert captured == {
        "qasm": QASM,
        "options": {
            "level": "Osuper",
            "passes": None,
            "fixpoint": True,
            "decompose_rz": True,
            "decompose_cz": True,
            "decompose_ccx": True,
            "superopt_gates": "base",
            "rz_epsilon": 0.001,
            "parallel": True,
        },
    }


@pytest.mark.parametrize(
    ("flag", "expected_level"),
    [("-O1", "O1"), ("-O2", "O2"), ("-O3", "O3"), ("-Osuper", "Osuper")],
)
def test_main_forwards_every_optimization_level(
    monkeypatch, tmp_path, flag, expected_level
):
    source = tmp_path / "input.qasm"
    source.write_text(QASM, encoding="utf-8")
    captured = {}

    def optimize_qasm(_qasm, **options):
        captured.update(options)
        return fake_result()

    monkeypatch.setattr(_cli, "optimize_qasm", optimize_qasm)

    assert _cli.main([str(source), flag]) == 0
    assert captured["level"] == expected_level


def test_main_normalizes_comma_separated_passes(monkeypatch, tmp_path):
    source = tmp_path / "input.qasm"
    source.write_text(QASM, encoding="utf-8")
    captured = {}

    def optimize_qasm(_qasm, **options):
        captured.update(options)
        return fake_result()

    monkeypatch.setattr(_cli, "optimize_qasm", optimize_qasm)

    assert (
        _cli.main(
            [
                str(source),
                "--passes",
                " CancelGates, , PhaseFoldRand ",
            ]
        )
        == 0
    )
    assert captured["passes"] == ["CancelGates", "PhaseFoldRand"]
