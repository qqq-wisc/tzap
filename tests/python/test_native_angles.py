import math

import numpy as np
import pennylane as qml
import pytest
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
from tzap._angles import angle_radians
from tzap.pennylane import optimize as optimize_tape
from tzap.qiskit import optimize


@pytest.mark.parametrize("axis", ["rx", "ry", "p"])
def test_qiskit_rebuilds_symbolic_native_angles(axis):
    circuit = QuantumCircuit(1)
    circuit.t(0)
    if axis in ("rx", "ry"):
        circuit.h(0)
    if axis == "ry":
        circuit.s(0)
    getattr(circuit, axis)(0.125, 0)
    result = optimize(circuit, passes=["PhaseFoldPauli"])
    assert Operator(circuit).equiv(Operator(result))
    assert len(result.data) < len(circuit.data)


@pytest.mark.parametrize("axis", ["rx", "ry", "p"])
def test_pennylane_rebuilds_symbolic_native_angles(axis):
    gates = [qml.T(0)]
    if axis in ("rx", "ry"):
        gates.append(qml.Hadamard(0))
    if axis == "ry":
        gates.append(qml.S(0))
    kind = {"rx": qml.RX, "ry": qml.RY, "p": qml.PhaseShift}[axis]
    gates.append(kind(0.125, wires=0))
    tape = qml.tape.QuantumScript(gates)
    (result,), _ = optimize_tape(tape, passes=["PhaseFoldPauli"])
    before, after = qml.matrix(tape), qml.matrix(result)
    overlap = np.vdot(before, after)
    assert np.allclose(after, before * overlap / abs(overlap), atol=1e-12)
    assert len(result.operations) < len(tape.operations)


def test_emitted_angle_conversion_accepts_only_constant_arithmetic():
    assert angle_radians("(1*pi/4)+(0.125)") == math.pi / 4 + 0.125
    assert angle_radians("-(pi/7)") == -math.pi / 7
    for expression in ["abs(1)", "__import__('os')", "pi.real", "2**3", "'1'"]:
        with pytest.raises(ValueError):
            angle_radians(expression)
