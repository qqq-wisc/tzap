# PennyLane API

tzap integrates with PennyLane as a transform for quantum functions, QNodes,
and `QuantumScript` tapes. It replaces only the circuit's operations and
retains the execution-facing tape data.

## Installation

Install the Python library:

```bash
pip install tzap
```

This installs PennyLane, Qiskit, the native tzap Python extension, and both
framework adapters. Python 3.10 or later is required. A Rust compiler is not
required when a tzap wheel is available for your platform.

## QNode decorator

Place `optimize` above `qml.qnode` so tzap receives the constructed QNode:

```python
import pennylane as qml
from tzap.pennylane import optimize

device = qml.device("default.qubit", wires=2)

@optimize(level="O3")
@qml.qnode(device)
def circuit():
    qml.Hadamard(0)
    qml.Hadamard(0)
    qml.T(1)
    qml.T(1)
    return qml.probs(wires=[0, 1])

print(qml.draw(circuit)())
print(circuit())
```

The transform can also be applied functionally:

```python
optimized_qnode = optimize(qnode, level="O1")
```

## Quantum functions and tapes

The same decorator works on a quantum function:

```python
@optimize(level="O3")
def ansatz():
    qml.Hadamard(0)
    qml.Hadamard(0)
    qml.CNOT(wires=[0, 1])
```

On a `QuantumScript`, PennyLane's transform protocol returns a batch and a
postprocessing function:

```python
tape = qml.tape.QuantumScript(
    [
        qml.Hadamard(0),
        qml.Hadamard(0),
        qml.T(1),
        qml.T(1),
    ],
    [qml.probs(wires=[0, 1])],
)

optimized_tapes, postprocess = optimize(tape, level="O3")
optimized_tape = optimized_tapes[0]
```

`tzap_transform` is an alias for `optimize`.

## Supported operations

The input circuit must use the gates listed under
[Circuit support](../README.md#circuit-support). Decompose other PennyLane
operations into that basis before applying tzap.

`RZ` angles must be concrete, finite, real scalars. Trainable or traced angles
are rejected because discrete circuit rewriting cannot safely preserve their
autodiff behavior.

For example, put trainable rotations outside the tzap-transformed portion of a
hybrid circuit, or apply tzap only after binding all rotation values.

## Measurements and circuit data

The transform preserves:

- arbitrary hashable wire labels;
- shots and shot vectors;
- terminal measurements and observables;
- existing `GlobalPhase` operations;
- unconditioned mid-circuit measurements and reset.

Postselection and classical feed-forward are rejected as unsupported dynamic
circuits. The optimized circuit is equivalent up to global phase.

## API and options

```python
optimize(tape=None, **options)
tzap_transform(tape=None, **options)
```

Common keyword options are:

| Option | Default | Meaning |
|---|---:|---|
| `level` | `"O3"` | `"O1"`, `"O2"`, `"O3"`, or `"Osuper"` |
| `passes` | `None` | Explicit ordered pass names; replaces `level`'s pipeline |
| `fixpoint` | `False` | Run `O1` or an explicit pass pipeline until gate count stops decreasing |
| `decompose_rz` | `False` | Approximate `RZ` operations with Clifford+T |
| `decompose_cz` | `False` | Lower `CZ` to `Hadamard` + `CNOT` + `Hadamard` |
| `decompose_ccx` | `False` | Lower both `Toffoli` and `CCZ` to Clifford+T |
| `rz_epsilon` | `1e-10` | Approximation tolerance used by `decompose_rz` |
| `parallel` | `False` | Enable the native parallel optimizer |
| `superopt_gates` | `"auto"` | MURM basis: `auto`, `base`, or an exact comma-separated gate list |

Explicit pass names are the same as the CLI's `--passes` values:
`DecomposeToffoli`, `DecomposeCz`, `DecomposeRz`, `CancelGates`, `SuperOpt`,
`PhaseFoldRand`, and `CnotMin`.

For example:

```python
optimized_qnode = optimize(
    qnode,
    passes=["CancelGates", "PhaseFoldRand"],
    fixpoint=True,
)
```

`Toffoli`, `CCZ`, `CZ`, and `RZ` remain native by default. Requested
decompositions run between native and post-decomposition optimization stages.

An explicit `passes` list cannot be combined with `decompose_rz`,
`decompose_cz`, or `decompose_ccx`; put the corresponding pass directly in the list
instead. `O2` always runs two rounds, while `O3` and `Osuper` already run to a
fixpoint, so `fixpoint` only changes `O1` or an explicit pipeline.

## Error handling

Unsupported PennyLane circuits raise `PennyLaneError`, a `ValueError`
subclass. Native parsing or optimization failures use tzap's `QasmError` or
`OptimizationError`:

```python
from tzap import OptimizationError, QasmError
from tzap.pennylane import PennyLaneError, optimize

try:
    optimized_qnode = optimize(qnode)
except PennyLaneError as error:
    print(f"unsupported PennyLane circuit: {error}")
except (QasmError, OptimizationError) as error:
    print(f"tzap failed: {error}")
```

See the repository [README](../README.md) for installation alternatives and
CLI usage.
