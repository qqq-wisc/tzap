import math
from pathlib import Path

import pytest
import tzap.qiskit as tzap_qiskit
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qasm2
from qiskit.circuit import Instruction, Parameter
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.transpiler.exceptions import TranspilerError
from tzap import optimize_qasm
from tzap.qiskit import TzapPass, _dag_to_qasm, optimize


def test_pass_has_qiskit_transformation_pass_contract():
    pass_ = TzapPass(level="O1")

    assert isinstance(pass_, TransformationPass)
    assert pass_.is_transformation_pass
    assert not pass_.is_analysis_pass
    assert tzap_qiskit.__all__ == ["TzapPass", "optimize"]


def test_transformation_pass_optimizes_and_preserves_structure():
    qreg = QuantumRegister(2, "logical")
    creg = ClassicalRegister(1, "readout")
    circuit = QuantumCircuit(
        qreg,
        creg,
        name="kept-name",
        global_phase=math.pi / 7,
        metadata={"purpose": "test"},
    )
    circuit.h(qreg[0])
    circuit.h(qreg[0])
    circuit.t(qreg[1])
    circuit.t(qreg[1])
    circuit.measure(qreg[1], creg[0])

    optimized = PassManager([TzapPass(level="O1")]).run(circuit)

    assert optimized.count_ops() == {"s": 1, "measure": 1}
    assert optimized.name == "kept-name"
    assert optimized.metadata == {"purpose": "test"}
    assert float(optimized.global_phase) == pytest.approx(math.pi / 7)
    assert [register.name for register in optimized.qregs] == ["logical"]
    assert [register.name for register in optimized.cregs] == ["readout"]


def test_pass_and_convenience_function():
    circuit = QuantumCircuit(1)
    circuit.x(0)
    circuit.x(0)

    assert isinstance(TzapPass(level="O1"), TransformationPass)
    assert len(optimize(circuit, level="O1").data) == 0


def test_convenience_function_does_not_mutate_input():
    circuit = QuantumCircuit(1)
    circuit.h(0)
    circuit.h(0)

    optimized = optimize(circuit, level="O1")

    assert len(circuit.data) == 2
    assert len(optimized.data) == 0


def test_run_accepts_and_returns_a_dag():
    circuit = QuantumCircuit(1)
    circuit.x(0)
    circuit.x(0)

    output = TzapPass(level="O1").run(circuit_to_dag(circuit))

    assert isinstance(output, DAGCircuit)
    assert len(output.op_nodes()) == 0


def test_run_serializes_before_no_copy_conversions(monkeypatch):
    circuit = QuantumCircuit(1)
    circuit.x(0)
    circuit.x(0)
    dag = circuit_to_dag(circuit)
    events = []

    def record_dag_to_qasm(input_dag):
        events.append(("dag_to_qasm", None))
        return _dag_to_qasm(input_dag)

    def record_dag_to_circuit(input_dag, *, copy_operations=True):
        events.append(("dag_to_circuit", copy_operations))
        return dag_to_circuit(input_dag, copy_operations=copy_operations)

    def record_circuit_to_dag(input_circuit, *, copy_operations=True):
        events.append(("circuit_to_dag", copy_operations))
        return circuit_to_dag(input_circuit, copy_operations=copy_operations)

    monkeypatch.setattr(tzap_qiskit, "_dag_to_qasm", record_dag_to_qasm)
    monkeypatch.setattr(tzap_qiskit, "dag_to_circuit", record_dag_to_circuit)
    monkeypatch.setattr(tzap_qiskit, "circuit_to_dag", record_circuit_to_dag)

    output = TzapPass(level="O1").run(dag)

    assert len(output.op_nodes()) == 0
    assert events == [
        ("dag_to_qasm", None),
        ("dag_to_circuit", False),
        ("circuit_to_dag", False),
    ]


def test_empty_qiskit_circuit_round_trip():
    circuit = QuantumCircuit(3, name="empty", metadata={"empty": True})

    optimized = optimize(circuit, level="O1")

    assert optimized.num_qubits == 3
    assert len(optimized.data) == 0
    assert optimized.name == "empty"
    assert optimized.metadata == {"empty": True}


def test_every_supported_qiskit_operation_round_trips():
    circuit = QuantumCircuit(14, 1)
    circuit.x(0)
    circuit.h(1)
    circuit.s(2)
    circuit.sdg(3)
    circuit.z(4)
    circuit.t(5)
    circuit.tdg(6)
    circuit.rz(0.321, 7)
    circuit.cx(8, 9)
    circuit.cz(10, 11)
    circuit.ccx(0, 1, 2)
    circuit.ccz(3, 4, 5)
    circuit.reset(12)
    circuit.measure(13, 0)

    optimized = optimize(circuit, passes=["CancelGates"])

    assert dict(optimized.count_ops()) == dict(circuit.count_ops())
    assert _canonical_operations(optimized) == _canonical_operations(circuit)


def test_bound_numeric_rz_round_trips():
    circuit = QuantumCircuit(1)
    circuit.rz(-0.123456789, 0)

    optimized = optimize(circuit, passes=["CancelGates"])

    assert float(optimized.data[0].operation.params[0]) == pytest.approx(-0.123456789)


def test_decompose_cz_option_reaches_native_optimizer():
    circuit = QuantumCircuit(2)
    circuit.cz(0, 1)

    optimized = optimize(circuit, level="O1", decompose_cz=True)

    assert "cz" not in optimized.count_ops()
    assert optimized.count_ops() == {"h": 2, "cx": 1}


def test_default_qiskit_pipeline_preserves_native_toffoli_gates():
    circuit = QuantumCircuit(3)
    circuit.ccx(0, 1, 2)
    circuit.ccz(0, 1, 2)

    optimized = optimize(circuit, level="O1")

    assert optimized.count_ops() == {"ccx": 1, "ccz": 1}


def test_qiskit_decompose_ccx_option_removes_toffoli_gates():
    circuit = QuantumCircuit(3)
    circuit.ccx(0, 1, 2)
    circuit.ccz(0, 1, 2)

    optimized = optimize(
        circuit,
        level="O1",
        decompose_ccx=True,
        superopt_gates="base",
    )

    assert "ccx" not in optimized.count_ops()
    assert "ccz" not in optimized.count_ops()
    assert len(optimized.data) > 2


def test_qiskit_decompose_rz_option_removes_rz():
    circuit = QuantumCircuit(1)
    circuit.rz(0.321, 0)

    optimized = optimize(
        circuit,
        level="O1",
        decompose_rz=True,
        rz_epsilon=1e-3,
    )

    assert "rz" not in optimized.count_ops()
    assert len(optimized.data) > 0


def test_passes_accept_generator_at_construction():
    circuit = QuantumCircuit(1)
    circuit.t(0)
    circuit.t(0)
    passes = (name for name in ["CancelGates", "PhaseFoldRand"])

    optimized = optimize(circuit, passes=passes)

    assert optimized.count_ops() == {"s": 1}


def test_multiple_qiskit_registers_and_measurement_mapping_are_preserved():
    left = QuantumRegister(1, "left")
    right = QuantumRegister(2, "right")
    first = ClassicalRegister(1, "first")
    second = ClassicalRegister(2, "second")
    circuit = QuantumCircuit(left, right, first, second)
    circuit.cx(left[0], right[1])
    circuit.measure(right[0], second[1])

    optimized = optimize(circuit, passes=["CancelGates"])

    assert [register.name for register in optimized.qregs] == ["left", "right"]
    assert [register.name for register in optimized.cregs] == ["first", "second"]
    assert _canonical_operations(optimized) == _canonical_operations(circuit)


def test_measurement_and_reset_order_on_same_wire_is_preserved():
    circuit = QuantumCircuit(1, 2)
    circuit.measure(0, 1)
    circuit.reset(0)
    circuit.measure(0, 0)

    optimized = optimize(circuit, passes=["CancelGates"])

    assert [instruction.operation.name for instruction in optimized.data] == [
        "measure",
        "reset",
        "measure",
    ]
    assert optimized.find_bit(optimized.data[0].clbits[0]).index == 1
    assert optimized.find_bit(optimized.data[2].clbits[0]).index == 0


def test_unsupported_gate_has_actionable_error():
    circuit = QuantumCircuit(1)
    circuit.y(0)

    with pytest.raises(TranspilerError, match="does not support.*'y'"):
        PassManager([TzapPass(level="O1")]).run(circuit)


@pytest.mark.parametrize("operation", ["barrier", "swap", "rx", "delay"])
def test_other_unsupported_operations_name_the_operation(operation):
    circuit = QuantumCircuit(2)
    if operation == "barrier":
        circuit.barrier()
    elif operation == "swap":
        circuit.swap(0, 1)
    elif operation == "rx":
        circuit.rx(0.2, 0)
    else:
        circuit.delay(10, 0)

    with pytest.raises(TranspilerError, match=f"'{operation}'"):
        optimize(circuit, level="O1")


def test_unbound_rz_is_rejected():
    circuit = QuantumCircuit(1)
    circuit.rz(Parameter("theta"), 0)

    with pytest.raises(TranspilerError, match="bound real numbers"):
        PassManager([TzapPass(level="O1")]).run(circuit)


@pytest.mark.parametrize("angle", [math.inf, -math.inf, math.nan])
def test_non_finite_rz_is_rejected(angle):
    circuit = QuantumCircuit(1)
    circuit.rz(angle, 0)

    with pytest.raises(TranspilerError, match="finite rz"):
        optimize(circuit, level="O1")


def test_control_flow_is_rejected_as_unsupported():
    circuit = QuantumCircuit(1, 1)
    with circuit.if_test((circuit.clbits[0], True)):
        circuit.x(0)

    with pytest.raises(TranspilerError, match="'if_else'"):
        optimize(circuit, level="O1")


def test_classically_conditioned_supported_gate_is_rejected():
    circuit = QuantumCircuit(1, 1)
    operation = Instruction("x", 1, 0, [])
    operation.condition = (circuit.cregs[0], 1)
    circuit.append(operation, [0])

    with pytest.raises(TranspilerError, match="classically conditioned"):
        _dag_to_qasm(circuit_to_dag(circuit))


def test_supported_gate_with_wrong_qubit_arity_is_rejected():
    circuit = QuantumCircuit(2)
    circuit.append(Instruction("x", 2, 0, []), [0, 1])

    with pytest.raises(TranspilerError, match="2 qubits, expected 1"):
        _dag_to_qasm(circuit_to_dag(circuit))


def test_measurement_without_classical_target_is_rejected():
    circuit = QuantumCircuit(1)
    circuit.append(Instruction("measure", 1, 0, []), [0])

    with pytest.raises(TranspilerError, match="exactly one classical bit"):
        _dag_to_qasm(circuit_to_dag(circuit))


def test_supported_quantum_gate_with_classical_operand_is_rejected():
    circuit = QuantumCircuit(1, 1)
    circuit.append(Instruction("x", 1, 1, []), [0], [0])

    with pytest.raises(TranspilerError, match="classical operands"):
        _dag_to_qasm(circuit_to_dag(circuit))


def test_invalid_native_options_propagate_through_direct_pass_run():
    circuit = QuantumCircuit(1)

    with pytest.raises(ValueError, match="optimization level"):
        TzapPass(level="invalid").run(circuit_to_dag(circuit))


def test_native_ccz_reset_and_measure_round_trip():
    circuit = QuantumCircuit(3, 1)
    circuit.ccz(0, 1, 2)
    circuit.reset(0)
    circuit.measure(2, 0)

    optimized = PassManager([TzapPass(passes=["CancelGates"])]).run(circuit)

    assert optimized.count_ops() == {"ccz": 1, "reset": 1, "measure": 1}


def test_tzap_input_prefers_original_order_for_independent_operations():
    circuit = QuantumCircuit(2)
    circuit.s(1)
    circuit.t(0)

    qasm = _dag_to_qasm(circuit_to_dag(circuit))

    assert qasm.index("s q[1]") < qasm.index("t q[0]")


def test_private_qasm_bridge_emits_flat_registers_and_all_operand_kinds():
    circuit = QuantumCircuit(3, 1)
    circuit.rz(0.25, 0)
    circuit.cx(0, 1)
    circuit.ccz(0, 1, 2)
    circuit.measure(2, 0)

    bridge = _dag_to_qasm(circuit_to_dag(circuit))

    assert "qreg q[3];" in bridge
    assert "creg c[1];" in bridge
    assert "rz(0.25) q[0];" in bridge
    assert "cx q[0],q[1];" in bridge
    assert "ccz q[0],q[1],q[2];" in bridge
    assert "measure q[2] -> c[0];" in bridge


def test_parallel_qiskit_run_preserves_measurements():
    circuit = QuantumCircuit(2, 2)
    circuit.measure(0, 1)
    circuit.reset(1)

    optimized = optimize(circuit, level="O1", parallel=True)

    assert optimized.count_ops() == {"measure": 1, "reset": 1}


@pytest.mark.parametrize("fixture_name", ["test.qasm", "two_ccx.qasm"])
def test_qiskit_matches_native_pipeline_on_repository_fixtures(fixture_name):
    fixture = Path(__file__).parents[1] / "fixtures" / fixture_name
    source_qasm = fixture.read_text(encoding="utf-8")

    native = optimize_qasm(source_qasm, level="O3")
    qiskit_output = optimize(qasm2.loads(source_qasm), level="O3")
    native_output = qasm2.loads(native.qasm)

    assert _canonical_operations(qiskit_output) == _canonical_operations(native_output)


def test_laplacian_filter_matches_native_cli_pipeline():
    benchmark = (
        Path(__file__).parents[2] / "benchmarks" / "cobble-t" / "laplacian-filter.qasm"
    )
    source_qasm = benchmark.read_text(encoding="utf-8")

    native = optimize_qasm(source_qasm, level="O3")
    qiskit_output = PassManager([TzapPass(level="O3")]).run(qasm2.loads(source_qasm))
    native_output = qasm2.loads(native.qasm)

    assert len(qiskit_output.data) == native.report.output.gates
    assert qiskit_output.depth() == native.report.output.depth
    assert _canonical_operations(qiskit_output) == _canonical_operations(native_output)


def _canonical_operations(circuit):
    circuit = dag_to_circuit(circuit_to_dag(circuit))
    return [
        (
            item.operation.name,
            tuple(float(parameter) for parameter in item.operation.params),
            tuple(circuit.find_bit(bit).index for bit in item.qubits),
            tuple(circuit.find_bit(bit).index for bit in item.clbits),
        )
        for item in circuit.data
    ]
