import dataclasses
import math
import re

import pytest
import tzap


def make_qasm(body, *, qubits=1, cbits=0):
    declarations = [f"qreg q[{qubits}];"]
    if cbits:
        declarations.append(f"creg c[{cbits}];")
    return "\n".join(
        [
            "OPENQASM 2.0;",
            'include "qelib1.inc";',
            *declarations,
            body.strip(),
            "",
        ]
    )


QASM = make_qasm(
    """
    h q[0];
    h q[0];
    t q[0];
    t q[0];
    """
)


def test_optimize_qasm_returns_qasm_and_metrics():
    result = tzap.optimize_qasm(QASM, level="O1")

    assert "h q[0]" not in result.qasm
    assert "s q[0]" in result.qasm
    assert result.report.input.gates == 4
    assert result.report.output.gates == 1
    assert result.report.input.t == 2
    assert result.report.output.t == 0


def test_public_version_is_semver_and_matches_native_module():
    from tzap import _native

    assert re.fullmatch(r"\d+\.\d+\.\d+(?:[-+].+)?", tzap.__version__)
    assert tzap.__version__ == _native.__version__


def test_public_all_exports_resolve():
    expected = {
        "Metrics",
        "OptimizationError",
        "OptimizationReport",
        "OptimizationResult",
        "QasmError",
        "TzapError",
        "__version__",
        "optimize_qasm",
    }

    assert set(tzap.__all__) == expected
    assert all(hasattr(tzap, name) for name in tzap.__all__)


def test_public_exception_hierarchy_and_modules():
    assert issubclass(tzap.QasmError, tzap.TzapError)
    assert issubclass(tzap.OptimizationError, tzap.TzapError)
    assert issubclass(tzap.TzapError, Exception)
    assert tzap.QasmError.__module__ == "tzap._native"
    assert tzap.OptimizationError.__module__ == "tzap._native"


def test_result_and_report_are_typed_frozen_dataclasses():
    result = tzap.optimize_qasm(QASM, level="O1")

    assert isinstance(result, tzap.OptimizationResult)
    assert isinstance(result.report, tzap.OptimizationReport)
    assert isinstance(result.report.output, tzap.Metrics)
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.report.output.gates = 100
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.qasm = "changed"


def test_all_metrics_are_populated():
    circuit = make_qasm(
        """
        h q[0];
        cx q[0],q[1];
        rz(0.3) q[1];
        """,
        qubits=2,
    )

    result = tzap.optimize_qasm(circuit, passes=["CancelGates"])

    assert result.report.input == tzap.Metrics(gates=3, two_qubit=1, depth=3, t=0, rz=1)
    assert result.report.baseline == result.report.input
    assert result.report.output == result.report.input


@pytest.mark.parametrize("level", ["O1", "o1", "1"])
def test_o1_level_spellings(level):
    result = tzap.optimize_qasm(QASM, level=level)

    assert result.report.output.gates == 1


@pytest.mark.parametrize("level", ["O2", "o2", "2", "O3", "o3", "3"])
def test_o2_and_o3_level_spellings(level):
    result = tzap.optimize_qasm(make_qasm("x q[0]; x q[0];"), level=level)

    assert result.report.output.gates == 0


def test_super_level_spelling_with_small_bounds():
    result = tzap.optimize_qasm(
        make_qasm("x q[0]; x q[0];"),
        level="super",
        superopt_qubits=1,
        superopt_window_gates=3,
        superopt_murm_entries=64,
    )

    assert result.report.output.gates == 0


def test_passes_accept_a_one_shot_generator():
    passes = (name for name in ["CancelGates", "PhaseFoldRand"])

    result = tzap.optimize_qasm(QASM, passes=passes)

    assert result.report.output.gates == 1


@pytest.mark.parametrize(
    "pass_name",
    [
        "DecomposeToffoli",
        "DecomposeCz",
        "DecomposeRz",
        "CancelGates",
        "SuperOpt",
        "PhaseFoldRand",
        "CnotMin",
    ],
)
def test_every_public_pass_name_is_accepted(pass_name):
    kwargs = {}
    if pass_name == "SuperOpt":
        kwargs = {
            "superopt_qubits": 1,
            "superopt_window_gates": 3,
            "superopt_murm_entries": 64,
        }

    result = tzap.optimize_qasm(make_qasm("x q[0];"), passes=[pass_name], **kwargs)

    assert result.report.input.gates == 1


def test_explicit_pipeline_reports_input_as_baseline():
    circuit = make_qasm("ccx q[0],q[1],q[2];", qubits=3)

    result = tzap.optimize_qasm(circuit, passes=["CancelGates"])

    assert result.report.input == result.report.baseline
    assert "ccx q[0],q[1],q[2]" in result.qasm


def test_default_pipeline_preserves_native_toffoli_and_baseline():
    circuit = make_qasm("ccx q[0],q[1],q[2];", qubits=3)

    result = tzap.optimize_qasm(circuit, level="O1")

    assert result.report.input.gates == 1
    assert result.report.baseline == result.report.input
    assert "ccx q[0],q[1],q[2]" in result.qasm


def test_decompose_ccx_removes_ccx_and_keeps_original_baseline():
    circuit = make_qasm("ccx q[0],q[1],q[2];", qubits=3)

    result = tzap.optimize_qasm(circuit, level="O1", decompose_ccx=True)

    assert result.report.input.gates == 1
    assert result.report.baseline == result.report.input
    assert "ccx " not in result.qasm


def test_decompose_cz_keeps_original_baseline_and_removes_cz():
    circuit = make_qasm("cz q[0],q[1];", qubits=2)

    result = tzap.optimize_qasm(circuit, level="O1", decompose_cz=True)

    assert result.report.input.gates == 1
    assert result.report.baseline == result.report.input
    assert "cz " not in result.qasm
    assert "cx q[0],q[1]" in result.qasm


def test_decompose_rz_removes_all_rz_gates():
    circuit = make_qasm("rz(0.321) q[0];")

    result = tzap.optimize_qasm(
        circuit,
        level="O1",
        decompose_rz=True,
        rz_epsilon=1e-3,
    )

    assert result.report.input.rz == 1
    assert result.report.output.rz == 0
    assert "rz(" not in result.qasm


def test_custom_pipeline_fixpoint_is_accepted():
    circuit = make_qasm("h q[0]; x q[0]; h q[0]; z q[0];")

    result = tzap.optimize_qasm(
        circuit,
        passes=["CancelGates"],
        fixpoint=True,
    )

    assert result.report.output.gates < result.report.input.gates


def test_parallel_mode_preserves_classical_register_and_measurement():
    circuit = make_qasm(
        """
        measure q[0] -> c[1];
        reset q[1];
        """,
        qubits=2,
        cbits=2,
    )

    result = tzap.optimize_qasm(circuit, level="O1", parallel=True)

    assert "creg c[2];" in result.qasm
    assert "measure q[0] -> c[1];" in result.qasm
    assert "reset q[1];" in result.qasm
    assert result.report.output.gates == 2


def test_multiple_quantum_and_classical_registers_are_flattened_correctly():
    circuit = """\
OPENQASM 2.0;
include "qelib1.inc";
qreg left[1];
qreg right[2];
creg first[1];
creg second[2];
cx left[0],right[1];
measure right[0] -> second[1];
"""

    result = tzap.optimize_qasm(circuit, passes=["CancelGates"])

    assert "qreg q[3];" in result.qasm
    assert "creg c[3];" in result.qasm
    assert "cx q[0],q[2];" in result.qasm
    assert "measure q[1] -> c[2];" in result.qasm


def test_empty_circuit_round_trip():
    result = tzap.optimize_qasm(make_qasm(""), level="O1")

    assert result.report.input == tzap.Metrics(0, 0, 0, 0, 0)
    assert result.report.output == result.report.input
    assert result.qasm.endswith("qreg q[1];\n")


def test_comments_and_angle_expressions_are_accepted():
    circuit = make_qasm(
        """
        // a line comment
        rz(pi / 2) q[0]; /* a block comment */
        """
    )

    result = tzap.optimize_qasm(circuit, passes=["CancelGates"])

    assert result.report.input.rz == 1
    assert "rz(" in result.qasm


def test_invalid_qasm_raises_specific_exception():
    with pytest.raises(tzap.QasmError, match="unsupported"):
        tzap.optimize_qasm("OPENQASM 2.0; qreg q[1]; y q[0];", level="O1")


@pytest.mark.parametrize(
    ("source", "message"),
    [
        ("OPENQASM 2.0; qreg q[1]; y q[0];", "unsupported"),
        ("OPENQASM 2.0; qreg q[1]; h q[2];", "out of range"),
        ("OPENQASM 2.0; qreg q[1]; cx q[0];", "expects 2"),
        ("OPENQASM 2.0; qreg q[1]; rz() q[0];", "angle"),
        ("OPENQASM 2.0; qreg q[1]; measure q[0];", "measure"),
    ],
)
def test_qasm_parse_failures_never_escape_as_generic_errors(source, message):
    with pytest.raises(tzap.QasmError, match=message):
        tzap.optimize_qasm(source, level="O1")


@pytest.mark.parametrize("level", ["nope", "O4", "", "O 1"])
def test_invalid_level(level):
    with pytest.raises(ValueError, match="optimization level"):
        tzap.optimize_qasm(QASM, level=level)


def test_empty_pass_list_is_rejected():
    with pytest.raises(ValueError, match="at least one"):
        tzap.optimize_qasm(QASM, passes=[])


def test_unknown_pass_lists_available_passes():
    with pytest.raises(ValueError, match="unknown pass.*CancelGates"):
        tzap.optimize_qasm(QASM, passes=["NotAPass"])


@pytest.mark.parametrize("option", ["decompose_rz", "decompose_cz", "decompose_ccx"])
def test_explicit_passes_reject_conflicting_options(option):
    with pytest.raises(ValueError, match="passes cannot be combined"):
        tzap.optimize_qasm(QASM, passes=["CancelGates"], **{option: True})


@pytest.mark.parametrize("epsilon", [0.0, -1.0, math.inf, -math.inf, math.nan])
def test_invalid_rz_epsilon_is_rejected(epsilon):
    with pytest.raises(ValueError, match="positive, finite"):
        tzap.optimize_qasm(QASM, level="O1", rz_epsilon=epsilon)


@pytest.mark.parametrize(
    "option",
    [
        "superopt_qubits",
        "superopt_window_gates",
        "superopt_murm_entries",
    ],
)
def test_zero_superopt_bound_is_rejected(option):
    with pytest.raises(ValueError, match=option):
        tzap.optimize_qasm(QASM, level="O1", **{option: 0})


@pytest.mark.parametrize("basis", ["", "rz", "h,nope"])
def test_invalid_superopt_gate_basis_is_rejected(basis):
    with pytest.raises(ValueError, match="superopt|SuperOpt|basis"):
        tzap.optimize_qasm(QASM, level="O1", superopt_gates=basis)


@pytest.mark.parametrize("basis", ["auto", "base", "h,h,t,cx,cz"])
def test_superopt_gate_modes_are_accepted(basis):
    result = tzap.optimize_qasm(QASM, level="O1", superopt_gates=basis)
    assert result.report.output.gates == 1


def test_non_string_qasm_is_rejected_by_binding_type_check():
    with pytest.raises(TypeError):
        tzap.optimize_qasm(b"not a Python string")


def test_non_string_pass_item_is_rejected_by_binding_type_check():
    with pytest.raises(TypeError):
        tzap.optimize_qasm(QASM, passes=["CancelGates", 123])
