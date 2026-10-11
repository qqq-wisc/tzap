"""Randomized end-to-end equivalence checks for the Python integrations."""

from __future__ import annotations

import math
import random
from typing import Any

import pennylane as qml
import pytest
from qiskit import QuantumCircuit
from tzap.pennylane import optimize as optimize_pennylane
from tzap.qiskit import optimize as optimize_qiskit

try:
    from mqt import qcec
    from mqt.qcec.pyqcec import EquivalenceCriterion
except ImportError:
    pytest.skip(
        "MQT QCEC is optional; install tzap[qcec] to run equivalence checks",
        allow_module_level=True,
    )


Gate = tuple[str, tuple[int, ...], float | None]

ACCEPTED = {
    EquivalenceCriterion.equivalent,
    EquivalenceCriterion.equivalent_up_to_global_phase,
}
SEEDS = range(8)
OPTIMIZER_CONFIGS = (
    pytest.param({"level": "O1"}, id="O1"),
    pytest.param({"level": "O2"}, id="O2"),
    pytest.param({"level": "O3"}, id="O3"),
    pytest.param({"level": "Osuper"}, id="Osuper"),
    pytest.param({"level": "O3", "fixpoint": True}, id="O3-fixpoint"),
    pytest.param({"level": "O3", "parallel": True}, id="O3-parallel"),
    pytest.param({"level": "O1", "decompose_cz": True}, id="O1-decompose-cz"),
    pytest.param({"level": "O1", "decompose_ccx": True}, id="O1-decompose-ccx"),
    pytest.param(
        {
            "level": "O2",
            "superopt_gates": "base",
            "superopt_qubits": 3,
            "superopt_window_gates": 4,
            "superopt_murm_entries": 500,
        },
        id="O2-superopt-base",
    ),
    pytest.param(
        {
            "level": "O2",
            "superopt_gates": "h,t,tdg,cx,cz,ccx,ccz",
            "superopt_qubits": 3,
            "superopt_window_gates": 4,
            "superopt_murm_entries": 500,
        },
        id="O2-superopt-explicit-native",
    ),
)


@pytest.mark.qcec
@pytest.mark.parametrize(
    "adapter",
    [
        pytest.param("qiskit", id="qiskit"),
        pytest.param("pennylane", id="pennylane"),
    ],
)
@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("options", OPTIMIZER_CONFIGS)
def test_randomized_python_adapter_equivalence(adapter, seed, options):
    """Exercise every pipeline through both third-party circuit adapters.

    Eight deterministic programs, ten optimizer configurations, and two
    adapters produce 160 independent QCEC proofs per Python version.
    """

    num_qubits, program = _random_program(seed)
    original = _qiskit_circuit(program, num_qubits)
    optimized = _optimize_with(adapter, program, num_qubits, options)

    result = qcec.verify(original, optimized)

    assert result.equivalence in ACCEPTED, (
        f"{adapter} seed={seed} options={options} was {result.equivalence}"
    )


@pytest.mark.qcec
def test_qcec_negative_control_rejects_a_changed_unitary():
    """Ensure the external checker is capable of failing this test suite."""

    original = QuantumCircuit(3)
    original.h(0)
    original.cx(0, 1)
    original.t(2)
    changed = original.copy()
    changed.x(2)

    assert qcec.verify(original, changed).equivalence not in ACCEPTED


def _optimize_with(
    adapter: str,
    program: list[Gate],
    num_qubits: int,
    options: dict[str, Any],
) -> QuantumCircuit:
    if adapter == "qiskit":
        return optimize_qiskit(_qiskit_circuit(program, num_qubits), **options)

    tape = qml.tape.QuantumScript(_pennylane_operations(program))
    batch, _ = optimize_pennylane(tape, **options)
    assert len(batch) == 1
    return _pennylane_to_qiskit(batch[0], num_qubits)


def _random_program(seed: int) -> tuple[int, list[Gate]]:
    rng = random.Random(seed)
    num_qubits = 3 + seed % 3
    program: list[Gate] = [
        # Guaranteed rewrite opportunities ensure the adapters exercise both
        # gate removal and gate reconstruction rather than merely round-trip.
        ("h", (0,), None),
        ("h", (0,), None),
        ("x", (1,), None),
        ("x", (1,), None),
        ("t", (2,), None),
        ("t", (2,), None),
        ("cx", (0, 1), None),
        ("cx", (0, 1), None),
        ("cz", (1, 2), None),
        ("cz", (1, 2), None),
    ]

    single_qubit = ("x", "h", "s", "sdg", "z", "t", "tdg")
    for _ in range(48 + seed * 2):
        kind = rng.randrange(10)
        if kind < 5:
            program.append(
                (rng.choice(single_qubit), (rng.randrange(num_qubits),), None)
            )
        elif kind == 5:
            # Irrational-looking and Clifford-aligned rotations cover both
            # numeric serialization and phase-folding behavior.
            angle = rng.choice(
                (
                    -math.pi,
                    -math.pi / 4,
                    math.pi / 8,
                    math.pi / 2,
                    0.123456789,
                    -0.987654321,
                )
            )
            program.append(("rz", (rng.randrange(num_qubits),), angle))
        elif kind < 9:
            left, right = rng.sample(range(num_qubits), 2)
            program.append((rng.choice(("cx", "cz")), (left, right), None))
        else:
            controls_and_target = tuple(rng.sample(range(num_qubits), 3))
            program.append((rng.choice(("ccx", "ccz")), controls_and_target, None))

    return num_qubits, program


def _qiskit_circuit(program: list[Gate], num_qubits: int) -> QuantumCircuit:
    circuit = QuantumCircuit(num_qubits)
    for name, wires, angle in program:
        if angle is not None:
            getattr(circuit, name)(angle, *wires)
        else:
            getattr(circuit, name)(*wires)
    return circuit


def _pennylane_operations(program: list[Gate]):
    constructors = {
        "x": qml.PauliX,
        "h": qml.Hadamard,
        "s": qml.S,
        "z": qml.PauliZ,
        "t": qml.T,
        "rz": qml.RZ,
        "rx": qml.RX,
        "ry": qml.RY,
        "p": qml.PhaseShift,
        "y": qml.PauliY,
        "sx": qml.SX,
        "swap": qml.SWAP,
        "cy": qml.CY,
        "ch": qml.CH,
        "cp": qml.ControlledPhaseShift,
        "crx": qml.CRX,
        "cry": qml.CRY,
        "crz": qml.CRZ,
        "cswap": qml.CSWAP,
        "cx": qml.CNOT,
        "cz": qml.CZ,
        "ccx": qml.Toffoli,
        "ccz": qml.CCZ,
    }
    operations = []
    with qml.QueuingManager.stop_recording():
        for name, wires, angle in program:
            operation_wires: int | tuple[int, ...]
            operation_wires = wires[0] if len(wires) == 1 else wires
            if name == "sdg":
                operations.append(qml.adjoint(qml.S(wires=operation_wires)))
            elif name == "tdg":
                operations.append(qml.adjoint(qml.T(wires=operation_wires)))
            elif angle is not None:
                operations.append(constructors[name](angle, wires=operation_wires))
            else:
                operations.append(constructors[name](wires=operation_wires))
    return operations


def _pennylane_to_qiskit(tape, num_qubits: int) -> QuantumCircuit:
    names = {
        "PauliX": "x",
        "Hadamard": "h",
        "S": "s",
        "Adjoint(S)": "sdg",
        "PauliZ": "z",
        "T": "t",
        "Adjoint(T)": "tdg",
        "RZ": "rz",
        "RX": "rx",
        "RY": "ry",
        "PhaseShift": "p",
        "PauliY": "y",
        "SX": "sx",
        "SWAP": "swap",
        "CY": "cy",
        "CH": "ch",
        "ControlledPhaseShift": "cp",
        "CRX": "crx",
        "CRY": "cry",
        "CRZ": "crz",
        "CSWAP": "cswap",
        "CNOT": "cx",
        "CZ": "cz",
        "Toffoli": "ccx",
        "CCZ": "ccz",
    }
    circuit = QuantumCircuit(num_qubits)
    for operation in tape.operations:
        name = names[operation.name]
        wires = tuple(int(wire) for wire in operation.wires)
        if operation.data:
            getattr(circuit, name)(float(operation.data[0]), *wires)
        else:
            getattr(circuit, name)(*wires)
    return circuit


@pytest.mark.qcec
@pytest.mark.parametrize("adapter", ["qiskit", "pennylane"])
@pytest.mark.parametrize("theta", [0.375, 2 * math.pi])
@pytest.mark.parametrize("level", ["O1", "O2"])
@pytest.mark.parametrize(
    "gate",
    [
        "p",
        "y",
        "sx",
        "rx",
        "ry",
        "swap",
        "cy",
        "ch",
        "cp",
        "crx",
        "cry",
        "crz",
        "cswap",
    ],
)
def test_native_gate_adapter_equivalence(adapter, theta, level, gate):
    arity = (
        3
        if gate == "cswap"
        else 2
        if gate in {"swap", "cy", "ch", "cp", "crx", "cry", "crz"}
        else 1
    )
    angle = theta if gate in {"p", "rx", "ry", "cp", "crx", "cry", "crz"} else None
    program = [
        ("h", (0,), None),
        ("h", (1,), None),
        ("h", (2,), None),
        ("cx", (0, 2), None),
        (gate, tuple(range(arity)), angle),
        ("t", (1,), None),
        ("h", (1,), None),
        ("cx", (0, 2), None),
    ]
    original = _qiskit_circuit(program, 3)
    optimized = _optimize_with(
        adapter,
        program,
        3,
        {
            "level": level,
            "superopt_qubits": 3,
            "superopt_window_gates": 4,
            "superopt_murm_entries": 500,
        },
    )
    assert qcec.verify(original, optimized).equivalence in ACCEPTED
