"""Qiskit transpiler integration for tzap."""

from __future__ import annotations

import math
from collections.abc import Iterable
from typing import Any

from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.transpiler import PassManager, TransformationPass
from qiskit.transpiler.exceptions import TranspilerError

from ._core import optimize_qasm

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
    "measure": 1,
    "reset": 1,
}


def _dag_to_qasm(dag: DAGCircuit) -> str:
    """Serialize exactly tzap's gate subset, with flat stable bit indices."""

    qubit_indices = {bit: index for index, bit in enumerate(dag.qubits)}
    cbit_indices = {bit: index for index, bit in enumerate(dag.clbits)}
    lines = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f"qreg q[{len(dag.qubits)}];",
    ]
    if dag.clbits:
        lines.append(f"creg c[{len(dag.clbits)}];")

    # Qiskit's default topological tie-breaker sorts independent operations by
    # their bit operands. That is semantically valid, but it changes the
    # linearized gate stream passed to tzap and can alter its greedy window
    # choices compared with optimizing the original QASM through the CLI.
    # Prefer the DAG's operation-node order while still asking Qiskit for a
    # valid topological traversal. The fallback rank applies to input/output
    # nodes, which the topological sorter also presents to the key function.
    insertion_rank = {node: index for index, node in enumerate(dag.op_nodes())}
    nodes = dag.topological_op_nodes(
        key=lambda node: f"{insertion_rank.get(node, -1):020d}"
    )
    for node in nodes:
        operation = node.op
        name = operation.name
        if name not in _ARITIES:
            raise TranspilerError(
                f"tzap does not support Qiskit operation {name!r}; transpile to the "
                "basis [x, h, s, sdg, z, t, tdg, rz, cx, cz, ccx, ccz, "
                "measure, reset] before running TzapPass"
            )
        if getattr(operation, "condition", None) is not None:
            raise TranspilerError(
                "tzap does not support classically conditioned operations"
            )
        if len(node.qargs) != _ARITIES[name]:
            raise TranspilerError(
                f"operation {name!r} has {len(node.qargs)} qubits, expected {_ARITIES[name]}"
            )

        qubits = [qubit_indices[bit] for bit in node.qargs]
        operands = ",".join(f"q[{index}]" for index in qubits)
        if name == "rz":
            try:
                angle = float(operation.params[0])
            except (IndexError, TypeError, ValueError) as error:
                raise TranspilerError(
                    "tzap requires rz angles to be bound real numbers"
                ) from error
            if not math.isfinite(angle):
                raise TranspilerError("tzap requires finite rz angles")
            lines.append(f"rz({angle!r}) {operands};")
        elif name == "measure":
            if len(node.cargs) != 1:
                raise TranspilerError("measure must target exactly one classical bit")
            lines.append(f"measure {operands} -> c[{cbit_indices[node.cargs[0]]}];")
        else:
            if node.cargs:
                raise TranspilerError(
                    f"operation {name!r} unexpectedly has classical operands"
                )
            lines.append(f"{name} {operands};")

    return "\n".join(lines) + "\n"


def _rebuild_on_original_bits(
    original: QuantumCircuit, optimized_qasm: str
) -> QuantumCircuit:
    # tzap's serializer emits only primitive one-line statements. Parsing this
    # known output ourselves avoids Qiskit's QASM exporter/importer inventing
    # custom gate definitions (notably for CCZ).
    flat = QuantumCircuit(original.num_qubits, original.num_clbits)
    for raw_line in optimized_qasm.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(("OPENQASM", "include ", "qreg ", "creg ")):
            continue
        statement = line.removesuffix(";")
        if statement.startswith("measure "):
            source, target = statement[len("measure ") :].split(" -> ")
            flat.measure(_index(source), _classical_index(target))
            continue

        name, operands = statement.split(" ", 1)
        if name.startswith("rz("):
            close = name.rfind(")")
            angle = float(name[3:close])
            flat.rz(angle, _indices(operands)[0])
        elif name == "reset":
            flat.reset(_indices(operands)[0])
        else:
            getattr(flat, name)(*_indices(operands))

    rebuilt = original.copy_empty_like()
    for instruction in flat.data:
        qargs = [rebuilt.qubits[flat.find_bit(bit).index] for bit in instruction.qubits]
        cargs = [rebuilt.clbits[flat.find_bit(bit).index] for bit in instruction.clbits]
        rebuilt.append(instruction.operation, qargs, cargs)
    return rebuilt


def _index(operand: str) -> int:
    return int(operand[operand.index("[") + 1 : operand.index("]")])


def _classical_index(operand: str) -> int:
    return _index(operand)


def _indices(operands: str) -> list[int]:
    return [_index(operand.strip()) for operand in operands.split(",")]


class TzapPass(TransformationPass):
    """Run tzap as a Qiskit transformation pass.

    Input operations must already be in tzap's supported basis. Circuit name,
    metadata, global phase, quantum/classical registers, and bit identities
    are retained. Like tzap itself, circuit equivalence is considered up to
    global phase.
    """

    def __init__(
        self,
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
    ) -> None:
        super().__init__()
        self._options: dict[str, Any] = {
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

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        source_qasm = _dag_to_qasm(dag)
        original = dag_to_circuit(dag, copy_operations=False)
        result = optimize_qasm(source_qasm, **self._options)
        rebuilt = _rebuild_on_original_bits(original, result.qasm)
        return circuit_to_dag(rebuilt, copy_operations=False)


def optimize(circuit: QuantumCircuit, **options: Any) -> QuantumCircuit:
    """Optimize a Qiskit circuit with a one-pass Qiskit ``PassManager``."""

    return PassManager([TzapPass(**options)]).run(circuit)


__all__ = ["TzapPass", "optimize"]
