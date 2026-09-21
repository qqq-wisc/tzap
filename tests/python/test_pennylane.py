import math
from pathlib import Path

import numpy as np
import pennylane as qml
import pytest
from pennylane import numpy as pnp
from tzap import optimize_qasm
from tzap import pennylane as tzap_pennylane
from tzap.pennylane import (
    PennyLaneError,
    _concrete_angle,
    _output_operations,
    _tape_to_qasm,
    optimize,
    tzap_transform,
)


def transform_tape(tape, **options):
    batch, postprocessing = optimize(tape, **options)
    assert len(batch) == 1
    return batch[0], postprocessing


def operation_names(tape):
    return [operation.name for operation in tape.operations]


def test_transform_optimizes_a_quantum_script():
    tape = qml.tape.QuantumScript(
        [
            qml.Hadamard("wire"),
            qml.Hadamard("wire"),
            qml.T(7),
            qml.T(7),
        ],
        [qml.probs(wires=["wire", 7])],
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert operation_names(transformed) == ["S"]
    assert list(transformed.operations[0].wires) == [7]


def test_transform_alias_and_exports():
    from tzap import pennylane

    assert tzap_transform is optimize
    assert set(pennylane.__all__) == {
        "PennyLaneError",
        "optimize",
        "tzap_transform",
    }


def test_null_postprocessing_returns_the_single_execution_result():
    tape = qml.tape.QuantumScript([qml.PauliX(0)], [qml.probs(wires=0)])

    _, postprocessing = transform_tape(tape, level="O1")

    result = np.array([0.0, 1.0])
    assert postprocessing([result]) is result


def test_terminal_measurements_and_observables_are_preserved():
    measurements = [
        qml.expval(qml.PauliZ("alice")),
        qml.var(qml.PauliX("bob")),
        qml.probs(wires=["alice", "bob"]),
    ]
    tape = qml.tape.QuantumScript(
        [qml.Hadamard("alice"), qml.Hadamard("alice")],
        measurements,
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert transformed.measurements == tape.measurements
    assert operation_names(transformed) == []


def test_shots_are_preserved():
    tape = qml.tape.QuantumScript(
        [qml.PauliX(0)],
        [qml.sample(wires=0)],
        shots=[(10, 2), 25],
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert transformed.shots == tape.shots


def test_arbitrary_hashable_wire_labels_round_trip():
    labels = ("left", -3, ("ancilla", 1))
    tape = qml.tape.QuantumScript(
        [
            qml.CNOT(wires=[labels[0], labels[1]]),
            qml.CZ(wires=[labels[1], labels[2]]),
        ],
        [qml.probs(wires=labels)],
    )

    transformed, _ = transform_tape(tape, passes=["CancelGates"])

    assert [tuple(operation.wires) for operation in transformed.operations] == [
        (labels[0], labels[1]),
        (labels[1], labels[2]),
    ]


def test_measurement_only_wires_are_included_in_native_mapping():
    tape = qml.tape.QuantumScript([], [qml.probs(wires=["unused", 9])])

    qasm, wires, _ = _tape_to_qasm(tape)
    transformed, _ = transform_tape(tape, level="O1")

    assert wires == ("unused", 9)
    assert "qreg q[2];" in qasm
    assert transformed.measurements == tape.measurements


def test_wireless_empty_tape_round_trip():
    tape = qml.tape.QuantumScript()

    transformed, _ = transform_tape(tape, level="O1")

    assert transformed.operations == []
    assert transformed.measurements == []


def test_every_supported_operation_crosses_the_bridge():
    tape = qml.tape.QuantumScript(
        [
            qml.PauliX(0),
            qml.Hadamard(1),
            qml.S(2),
            qml.adjoint(qml.S(3)),
            qml.PauliZ(4),
            qml.T(5),
            qml.adjoint(qml.T(6)),
            qml.RZ(0.321, wires=7),
            qml.CNOT(wires=[8, 9]),
            qml.CZ(wires=[10, 11]),
            qml.Toffoli(wires=[0, 1, 2]),
            qml.CCZ(wires=[3, 4, 5]),
        ],
        [qml.probs(wires=range(12))],
    )

    transformed, _ = transform_tape(tape, passes=["CancelGates"])

    assert sorted(operation_names(transformed)) == sorted(operation_names(tape))


def test_native_qasm_bridge_uses_expected_gate_names():
    tape = qml.tape.QuantumScript(
        [
            qml.PauliX("a"),
            qml.adjoint(qml.S("a")),
            qml.adjoint(qml.T("a")),
            qml.RZ(-0.25, wires="a"),
            qml.CNOT(wires=["a", "b"]),
            qml.CCZ(wires=["a", "b", "c"]),
        ]
    )

    qasm, wires, phases = _tape_to_qasm(tape)

    assert wires == ("a", "b", "c")
    assert phases == ()
    assert "x q[0];" in qasm
    assert "sdg q[0];" in qasm
    assert "tdg q[0];" in qasm
    assert "rz(-0.25) q[0];" in qasm
    assert "cx q[0],q[1];" in qasm
    assert "ccz q[0],q[1],q[2];" in qasm


def test_output_bridge_reconstructs_every_gate():
    qasm = """\
OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
x q[0];
h q[0];
s q[0];
sdg q[0];
z q[0];
t q[0];
tdg q[0];
rz(0.125) q[0];
cx q[0],q[1];
cz q[0],q[1];
ccx q[0],q[1],q[2];
ccz q[0],q[1],q[2];
"""

    operations = _output_operations(qasm, ("a", "b", "c"))

    assert [operation.name for operation in operations] == [
        "PauliX",
        "Hadamard",
        "S",
        "Adjoint(S)",
        "PauliZ",
        "T",
        "Adjoint(T)",
        "RZ",
        "CNOT",
        "CZ",
        "Toffoli",
        "CCZ",
    ]
    assert list(operations[-1].wires) == ["a", "b", "c"]


def test_existing_global_phase_is_preserved():
    phase = qml.GlobalPhase(0.432)
    tape = qml.tape.QuantumScript(
        [qml.Hadamard(0), phase, qml.Hadamard(0)],
        [qml.probs(wires=0)],
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert len(transformed.operations) == 1
    assert transformed.operations[0] is phase


def test_trainable_global_phase_is_preserved_without_converting_to_float():
    phase = pnp.array(0.2, requires_grad=True)
    tape = qml.tape.QuantumScript(
        [qml.GlobalPhase(phase), qml.PauliX(0)],
        [qml.probs(wires=0)],
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert transformed.operations[0].data[0] is phase


def test_bound_non_trainable_rz_round_trips():
    angle = np.float64(-0.123456789)
    tape = qml.tape.QuantumScript([qml.RZ(angle, wires="q")])

    transformed, _ = transform_tape(tape, passes=["CancelGates"])

    assert transformed.operations[0].data[0] == pytest.approx(float(angle))
    assert list(transformed.operations[0].wires) == ["q"]


def test_trainable_rz_is_rejected_to_protect_autodiff():
    angle = pnp.array(0.2, requires_grad=True)
    tape = qml.tape.QuantumScript([qml.RZ(angle, wires=0)])

    with pytest.raises(PennyLaneError, match="non-trainable.*autodiff"):
        optimize(tape, level="O1")


@pytest.mark.parametrize("angle", [math.inf, -math.inf, math.nan])
def test_non_finite_rz_is_rejected(angle):
    tape = qml.tape.QuantumScript([qml.RZ(angle, wires=0)])

    with pytest.raises(PennyLaneError, match="finite RZ"):
        optimize(tape, level="O1")


def test_non_scalar_rz_is_rejected():
    tape = qml.tape.QuantumScript([qml.RZ(np.array([0.1, 0.2]), wires=0)])

    with pytest.raises(PennyLaneError, match="real scalars"):
        optimize(tape, level="O1")


def test_rz_without_an_angle_is_rejected():
    class MissingAngle:
        data = ()

    with pytest.raises(PennyLaneError, match="exactly one angle"):
        _concrete_angle(MissingAngle())


def test_supported_operation_with_wrong_arity_is_rejected(monkeypatch):
    tape = qml.tape.QuantumScript([qml.PauliX(0)])
    monkeypatch.setattr(tzap_pennylane, "_operation_name", lambda _operation: "cx")

    with pytest.raises(PennyLaneError, match="1 wires, expected 2"):
        _tape_to_qasm(tape)


@pytest.mark.parametrize(
    "operation",
    [
        qml.RX(0.2, wires=0),
        qml.RY(0.2, wires=0),
        qml.SWAP(wires=[0, 1]),
        qml.Identity(wires=0),
        qml.adjoint(qml.Hadamard(0)),
        qml.BasisState(np.array([1]), wires=0),
    ],
)
def test_unsupported_operations_have_actionable_errors(operation):
    tape = qml.tape.QuantumScript([operation])

    with pytest.raises(PennyLaneError, match="does not support PennyLane operation"):
        optimize(tape, level="O1")


def test_unconditioned_mid_circuit_measurement_is_preserved():
    operation = qml.measurements.MidMeasureMP(wires=0)
    tape = qml.tape.QuantumScript(
        [qml.T(0), qml.T(0), operation],
        [qml.probs(wires=0)],
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert operation_names(transformed) == ["S", "MidMeasureMP"]
    assert transformed.operations[-1] is operation


def test_output_bridge_reconstructs_measurement_without_original_tape_metadata():
    operations = _output_operations(
        "measure q[0] -> c[0];",
        ("wire",),
    )

    assert len(operations) == 1
    assert isinstance(operations[0], qml.measurements.MidMeasureMP)
    assert list(operations[0].wires) == ["wire"]


def test_output_bridge_reconstructs_standalone_reset():
    operations = _output_operations("reset q[0];", ("wire",))

    assert len(operations) == 1
    assert isinstance(operations[0], qml.measurements.MidMeasureMP)
    assert operations[0].reset
    assert list(operations[0].wires) == ["wire"]


def test_output_bridge_rejects_unknown_native_operation():
    with pytest.raises(PennyLaneError, match="unsupported.*'y'"):
        _output_operations("y q[0];", ("wire",))


def test_qnode_unconditioned_mid_circuit_measurement_executes():
    # Analytic execution defers the measurement onto an auxiliary wire.
    device = qml.device("default.qubit", wires=3)

    @optimize(level="O1")
    @qml.qnode(device)
    def circuit():
        qml.Hadamard(0)
        qml.measure(0)
        return qml.probs(wires=[0, 1])

    assert np.allclose(circuit(), np.array([0.5, 0.0, 0.5, 0.0]))


def test_reset_is_preserved_and_executes():
    reset = qml.measurements.MidMeasureMP(wires=0, reset=True)
    tape = qml.tape.QuantumScript(
        [qml.PauliX(0), reset],
        [qml.probs(wires=0)],
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert transformed.operations[-1] is reset
    assert transformed.operations[-1].reset

    # Analytic execution defers the reset's measurement onto an auxiliary wire.
    device = qml.device("default.qubit", wires=2)

    @optimize(level="O1")
    @qml.qnode(device)
    def circuit():
        qml.PauliX(0)
        qml.measure(0, reset=True)
        return qml.probs(wires=0)

    assert np.allclose(circuit(), np.array([1.0, 0.0]))


def test_postselection_is_rejected_as_dynamic():
    tape = qml.tape.QuantumScript(
        [qml.measurements.MidMeasureMP(wires=0, postselect=1)],
        [qml.probs(wires=0)],
    )

    with pytest.raises(PennyLaneError, match="dynamic circuits.*postselection"):
        optimize(tape, level="O1")


def test_qnode_classical_feed_forward_is_rejected():
    device = qml.device("default.qubit", wires=2)

    @optimize(level="O1")
    @qml.qnode(device)
    def circuit():
        measurement = qml.measure(0)
        qml.cond(measurement, qml.PauliX)(1)
        return qml.probs(wires=[0, 1])

    with pytest.raises(PennyLaneError, match="dynamic circuits"):
        circuit()


def test_decompose_cz_option_reaches_native_optimizer():
    tape = qml.tape.QuantumScript([qml.CZ(wires=["a", "b"])])

    transformed, _ = transform_tape(tape, level="O1", decompose_cz=True)

    assert operation_names(transformed) == ["Hadamard", "CNOT", "Hadamard"]


def test_default_pipeline_preserves_toffoli_and_ccz():
    tape = qml.tape.QuantumScript(
        [
            qml.Toffoli(wires=[0, 1, 2]),
            qml.CCZ(wires=[0, 1, 2]),
        ]
    )

    transformed, _ = transform_tape(tape, level="O1")

    assert operation_names(transformed) == ["Toffoli", "CCZ"]


def test_decompose_ccx_option_removes_toffoli_and_ccz():
    tape = qml.tape.QuantumScript(
        [
            qml.Toffoli(wires=[0, 1, 2]),
            qml.CCZ(wires=[0, 1, 2]),
        ]
    )

    transformed, _ = transform_tape(
        tape,
        level="O1",
        decompose_ccx=True,
        superopt_gates="base",
    )

    assert "Toffoli" not in operation_names(transformed)
    assert "CCZ" not in operation_names(transformed)
    assert len(transformed.operations) > 2


def test_decompose_rz_option_removes_rz():
    tape = qml.tape.QuantumScript([qml.RZ(0.321, wires=0)])

    transformed, _ = transform_tape(
        tape,
        level="O1",
        decompose_rz=True,
        rz_epsilon=1e-3,
    )

    assert "RZ" not in operation_names(transformed)
    assert transformed.operations


def test_custom_pipeline_accepts_pass_generator():
    tape = qml.tape.QuantumScript([qml.T(0), qml.T(0)])
    passes = (name for name in ["CancelGates", "PhaseFoldRand"])

    transformed, _ = transform_tape(tape, passes=passes)

    assert operation_names(transformed) == ["S"]


@pytest.mark.parametrize("level", ["O1", "O2", "O3"])
def test_standard_optimization_levels(level):
    tape = qml.tape.QuantumScript([qml.PauliX(0), qml.PauliX(0)])

    transformed, _ = transform_tape(tape, level=level)

    assert transformed.operations == []


def test_invalid_optimizer_options_propagate():
    tape = qml.tape.QuantumScript([qml.PauliX(0)])

    with pytest.raises(ValueError, match="optimization level"):
        optimize(tape, level="invalid")


def test_parallel_transform_path():
    tape = qml.tape.QuantumScript(
        [qml.PauliX(0), qml.PauliX(0)],
        [qml.probs(wires=0)],
    )

    transformed, _ = transform_tape(tape, level="O1", parallel=True)

    # Map-reduce chunking is permitted to optimize independently; this checks
    # that the transformed tape remains executable and measurements survive.
    assert transformed.measurements == tape.measurements
    assert all(
        isinstance(operation, qml.PauliX) for operation in transformed.operations
    )


def test_qnode_decorator_executes_optimized_circuit():
    device = qml.device("default.qubit", wires=2)

    @optimize(level="O1")
    @qml.qnode(device)
    def circuit():
        qml.Hadamard(0)
        qml.Hadamard(0)
        qml.T(1)
        qml.T(1)
        return qml.probs(wires=[0, 1])

    assert np.allclose(circuit(), np.array([1.0, 0.0, 0.0, 0.0]))
    assert "S" in qml.draw(circuit)()
    assert "H" not in qml.draw(circuit)()


def test_qnode_transform_preserves_expectation_value():
    device = qml.device("default.qubit", wires=2)

    @qml.qnode(device)
    def original():
        qml.Hadamard(0)
        qml.CNOT(wires=[0, 1])
        qml.T(1)
        qml.T(1)
        return qml.expval(qml.PauliZ(0) @ qml.PauliZ(1))

    transformed = optimize(original, level="O1")

    assert transformed() == pytest.approx(original())


def test_qnode_with_trainable_rz_fails_clearly_at_execution():
    device = qml.device("default.qubit", wires=1)

    @optimize(level="O1")
    @qml.qnode(device, interface="autograd")
    def circuit(angle):
        qml.RZ(angle, wires=0)
        return qml.expval(qml.PauliX(0))

    with pytest.raises(PennyLaneError, match="autodiff"):
        circuit(pnp.array(0.2, requires_grad=True))


def test_transformed_and_original_tapes_have_same_unitary():
    tape = qml.tape.QuantumScript(
        [
            qml.Hadamard(0),
            qml.CNOT(wires=[0, 1]),
            qml.T(1),
            qml.T(1),
            qml.CZ(wires=[0, 1]),
        ]
    )

    transformed, _ = transform_tape(tape, level="O3")

    original_matrix = qml.matrix(tape, wire_order=[0, 1])
    transformed_matrix = qml.matrix(transformed, wire_order=[0, 1])
    phase = (
        original_matrix.flat[np.argmax(np.abs(original_matrix))]
        / (transformed_matrix.flat[np.argmax(np.abs(original_matrix))])
    )
    assert np.allclose(original_matrix, phase * transformed_matrix)


def test_laplacian_filter_matches_native_pipeline():
    benchmark = (
        Path(__file__).parents[2] / "benchmarks" / "cobble-t" / "laplacian-filter.qasm"
    )
    source_qasm = benchmark.read_text(encoding="utf-8")
    tape = qml.tape.QuantumScript(_output_operations(source_qasm, tuple(range(11))))

    transformed, _ = transform_tape(tape, level="O3")
    native = optimize_qasm(source_qasm, level="O3")
    native_operations = _output_operations(native.qasm, tuple(range(11)))

    assert _operation_signatures(transformed.operations) == _operation_signatures(
        native_operations
    )
    assert len(transformed.operations) == 26_402


def _operation_signatures(operations):
    return [
        (
            operation.name,
            tuple(operation.wires),
            tuple(float(parameter) for parameter in operation.data),
        )
        for operation in operations
    ]
