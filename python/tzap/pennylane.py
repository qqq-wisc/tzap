"""PennyLane transform integration for tzap."""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence

import pennylane as qml
from pennylane.tape import QuantumScript

from ._core import optimize_qasm


class PennyLaneError(ValueError):
    """A PennyLane circuit cannot be represented safely by tzap."""


_SIMPLE_OPERATIONS = (
    (qml.PauliX, "x"),
    (qml.Hadamard, "h"),
    (qml.S, "s"),
    (qml.PauliZ, "z"),
    (qml.T, "t"),
    (qml.RZ, "rz"),
    (qml.CNOT, "cx"),
    (qml.CZ, "cz"),
    (qml.Toffoli, "ccx"),
    (qml.CCZ, "ccz"),
)

_ARITIES = {
    "x": 1,
    "h": 1,
    "s": 1,
    "sdg": 1,
    "z": 1,
    "t": 1,
    "tdg": 1,
    "rz": 1,
    "cx": 2,
    "cz": 2,
    "ccx": 3,
    "ccz": 3,
}


def _operation_name(operation) -> str:
    if isinstance(operation, qml.ops.op_math.Conditional):
        raise PennyLaneError(
            "tzap does not support dynamic circuits with classical feed-forward"
        )

    if isinstance(operation, qml.ops.op_math.Adjoint):
        if isinstance(operation.base, qml.S):
            return "sdg"
        if isinstance(operation.base, qml.T):
            return "tdg"
        raise PennyLaneError(
            f"tzap does not support PennyLane operation {operation.name!r}; only adjoints of "
            "S and T are supported"
        )

    for operation_type, name in _SIMPLE_OPERATIONS:
        if isinstance(operation, operation_type):
            return name
    raise PennyLaneError(
        f"tzap does not support PennyLane operation {operation.name!r}; decompose to "
        "[PauliX, Hadamard, S, Adjoint(S), PauliZ, T, Adjoint(T), RZ, "
        "CNOT, CZ, Toffoli, CCZ] before applying tzap"
    )


def _concrete_angle(operation) -> float:
    try:
        angle = operation.data[0]
    except (AttributeError, IndexError) as error:
        raise PennyLaneError("RZ must have exactly one angle") from error

    if qml.math.requires_grad(angle) or qml.math.is_abstract(angle):
        raise PennyLaneError(
            "tzap requires RZ angles to be concrete and non-trainable; "
            "optimizing a trainable or traced angle would break autodiff"
        )
    try:
        value = float(angle)
    except (TypeError, ValueError) as error:
        raise PennyLaneError("tzap requires RZ angles to be real scalars") from error
    if not math.isfinite(value):
        raise PennyLaneError("tzap requires finite RZ angles")
    return value


def _tape_to_qasm(
    tape: QuantumScript,
) -> tuple[str, Sequence[object], Sequence[object]]:
    wires = tuple(tape.wires)
    # QuantumScript derives wire order by first appearance. For conventional
    # integer or string labels, a canonical ordering avoids changing tzap's
    # greedy window choices merely because a high-numbered wire happened to
    # receive the first operation. Mixed/custom labels may not be orderable,
    # in which case PennyLane's own order remains the deterministic fallback.
    try:
        wires = tuple(sorted(wires))
    except TypeError:
        pass
    wire_indices: dict[object, int] = {wire: index for index, wire in enumerate(wires)}
    global_phases = []
    lines = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f"qreg q[{len(wires)}];",
    ]
    mid_measurements = [
        operation
        for operation in tape.operations
        if isinstance(operation, qml.measurements.MidMeasureMP)
    ]
    measurement_indices = {
        id(operation): index for index, operation in enumerate(mid_measurements)
    }
    if mid_measurements:
        lines.append(f"creg c[{len(mid_measurements)}];")

    for operation in tape.operations:
        if isinstance(operation, qml.GlobalPhase):
            global_phases.append(operation)
            continue

        if isinstance(operation, qml.measurements.MidMeasureMP):
            if operation.postselect is not None:
                raise PennyLaneError(
                    "tzap does not support dynamic circuits with postselection"
                )
            wire = wire_indices[operation.wires[0]]
            cbit = measurement_indices[id(operation)]
            lines.append(f"measure q[{wire}] -> c[{cbit}];")
            if operation.reset:
                lines.append(f"reset q[{wire}];")
            continue

        name = _operation_name(operation)
        operation_wires = tuple(operation.wires)
        if len(operation_wires) != _ARITIES[name]:
            raise PennyLaneError(
                f"operation {operation.name!r} has {len(operation_wires)} wires, expected {_ARITIES[name]}"
            )
        operands = ",".join(f"q[{wire_indices[wire]}]" for wire in operation_wires)
        if name == "rz":
            lines.append(f"rz({_concrete_angle(operation)!r}) {operands};")
        else:
            lines.append(f"{name} {operands};")

    return "\n".join(lines) + "\n", wires, tuple(global_phases)


def _index(operand: str) -> int:
    return int(operand[operand.index("[") + 1 : operand.index("]")])


def _indices(operands: str) -> list[int]:
    return [_index(operand.strip()) for operand in operands.split(",")]


def _output_operations(
    qasm: str,
    wires: Sequence[object],
    mid_measurements: Sequence[object] = (),
):
    operations = []
    resets_to_skip: dict[object, int] = {}
    with qml.QueuingManager.stop_recording():
        for raw_line in qasm.splitlines():
            line = raw_line.strip()
            if not line or line.startswith(("OPENQASM", "include ", "qreg ", "creg ")):
                continue

            statement = line.removesuffix(";")
            if statement.startswith("measure "):
                source, target = statement[len("measure ") :].split(" -> ")
                wire = wires[_index(source)]
                cbit = _index(target)
                if cbit < len(mid_measurements):
                    measurement = mid_measurements[cbit]
                else:
                    measurement = qml.measurements.MidMeasureMP(wires=wire)
                operations.append(measurement)
                if measurement.reset:
                    resets_to_skip[wire] = resets_to_skip.get(wire, 0) + 1
                continue

            name, operands = statement.split(" ", 1)
            wire_operands = [wires[index] for index in _indices(operands)]
            if name.startswith("rz("):
                operations.append(
                    qml.RZ(float(name[3 : name.rfind(")")]), wires=wire_operands[0])
                )
            elif name == "x":
                operations.append(qml.PauliX(wires=wire_operands[0]))
            elif name == "h":
                operations.append(qml.Hadamard(wires=wire_operands[0]))
            elif name == "s":
                operations.append(qml.S(wires=wire_operands[0]))
            elif name == "sdg":
                operations.append(qml.adjoint(qml.S(wires=wire_operands[0])))
            elif name == "z":
                operations.append(qml.PauliZ(wires=wire_operands[0]))
            elif name == "t":
                operations.append(qml.T(wires=wire_operands[0]))
            elif name == "tdg":
                operations.append(qml.adjoint(qml.T(wires=wire_operands[0])))
            elif name == "cx":
                operations.append(qml.CNOT(wires=wire_operands))
            elif name == "cz":
                operations.append(qml.CZ(wires=wire_operands))
            elif name == "ccx":
                operations.append(qml.Toffoli(wires=wire_operands))
            elif name == "ccz":
                operations.append(qml.CCZ(wires=wire_operands))
            elif name == "reset":
                wire = wire_operands[0]
                if resets_to_skip.get(wire, 0):
                    resets_to_skip[wire] -= 1
                else:
                    operations.append(
                        qml.measurements.MidMeasureMP(wires=wire, reset=True)
                    )
            else:
                raise PennyLaneError(
                    "tzap returned an operation unsupported by the PennyLane "
                    f"adapter: {name!r}"
                )
    return operations


@qml.transform
def _optimize_transform(
    tape: QuantumScript,
    *,
    level: str = "O3",
    passes: Iterable[str] | None = None,
    fixpoint: bool = False,
    decompose_rz: bool = False,
    decompose_cz: bool = False,
    decompose_ccx: bool = False,
    rz_epsilon: float = 1e-10,
    parallel: bool = False,
    superopt_qubits: int | None = None,
    superopt_window_gates: int | None = None,
    superopt_murm_entries: int | None = None,
    superopt_gates: str = "auto",
):
    """Optimize a PennyLane tape, quantum function, or QNode with tzap.

    Terminal measurements and observables are retained by copying the input
    ``QuantumScript`` with only its operation list replaced. Existing
    ``GlobalPhase`` operations are preserved. RZ parameters must be concrete,
    finite, and non-trainable.
    """

    qasm, wires, global_phases = _tape_to_qasm(tape)
    pass_list = None if passes is None else tuple(passes)
    result = optimize_qasm(
        qasm,
        level=level,
        passes=pass_list,
        fixpoint=fixpoint,
        decompose_rz=decompose_rz,
        decompose_cz=decompose_cz,
        decompose_ccx=decompose_ccx,
        rz_epsilon=rz_epsilon,
        parallel=parallel,
        superopt_qubits=superopt_qubits,
        superopt_window_gates=superopt_window_gates,
        superopt_murm_entries=superopt_murm_entries,
        superopt_gates=superopt_gates,
    )
    optimized_operations = [
        *global_phases,
        *_output_operations(
            result.qasm,
            wires,
            tuple(
                operation
                for operation in tape.operations
                if isinstance(operation, qml.measurements.MidMeasureMP)
            ),
        ),
    ]
    transformed_tape = tape.copy(operations=optimized_operations)

    def null_postprocessing(results):
        return results[0]

    return (transformed_tape,), null_postprocessing


def optimize(
    tape=None,
    *,
    level: str = "O3",
    passes: Iterable[str] | None = None,
    fixpoint: bool = False,
    decompose_rz: bool = False,
    decompose_cz: bool = False,
    decompose_ccx: bool = False,
    rz_epsilon: float = 1e-10,
    parallel: bool = False,
    superopt_qubits: int | None = None,
    superopt_window_gates: int | None = None,
    superopt_murm_entries: int | None = None,
    superopt_gates: str = "auto",
):
    """Optimize a PennyLane tape, quantum function, or QNode with tzap.

    Called without a circuit-like first argument, this returns a decorator.
    This wrapper keeps parameterized-decorator syntax consistent across
    PennyLane versions whose transform dispatchers differ on that syntax.
    """

    options = {
        "level": level,
        "passes": None if passes is None else tuple(passes),
        "fixpoint": fixpoint,
        "decompose_rz": decompose_rz,
        "decompose_cz": decompose_cz,
        "decompose_ccx": decompose_ccx,
        "rz_epsilon": rz_epsilon,
        "parallel": parallel,
        "superopt_qubits": superopt_qubits,
        "superopt_window_gates": superopt_window_gates,
        "superopt_murm_entries": superopt_murm_entries,
        "superopt_gates": superopt_gates,
    }

    def apply_transform(target):
        return _optimize_transform(target, **options)

    if tape is None:
        return apply_transform
    return apply_transform(tape)


tzap_transform = optimize

__all__ = ["PennyLaneError", "optimize", "tzap_transform"]
