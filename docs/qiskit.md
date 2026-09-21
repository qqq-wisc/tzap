# Qiskit API

tzap integrates with Qiskit as a transformation pass. It sends a circuit
through the same native optimizer and optimization levels as the tzap CLI,
then returns a Qiskit circuit.

## Installation

Install the Python library:

```bash
pip install tzap
```

This installs Qiskit, PennyLane, the native tzap Python extension, and both
framework adapters. Python 3.10 or later is required. A Rust compiler is not
required when a wheel is available for your platform.

## Quick start

The convenience function runs a one-pass `PassManager` and returns a new
`QuantumCircuit`. It does not mutate the input circuit.

```python
from qiskit import QuantumCircuit
from tzap.qiskit import optimize

circuit = QuantumCircuit(2)
circuit.h(0)
circuit.h(0)
circuit.t(1)
circuit.t(1)

optimized = optimize(circuit, level="O3")

print(optimized.count_ops())
```

Use `TzapPass` when composing tzap with other transpiler passes:

```python
from qiskit.transpiler import PassManager
from tzap.qiskit import TzapPass

manager = PassManager([
    TzapPass(level="O3"),
])
optimized = manager.run(circuit)
```

## Preparing a circuit

The input circuit must use the gates listed under
[Circuit support](../README.md#circuit-support). Use Qiskit's transpiler to
translate other operations into that basis before running tzap.

`rz` angles must be bound, finite real numbers. Symbolic parameters,
classically conditioned gates, control-flow operations, and operations
outside the supported basis raise Qiskit's `TranspilerError`.

## API and options

```python
optimize(circuit, **options) -> QuantumCircuit
TzapPass(**options)
```

Common keyword options are:

| Option | Default | Meaning |
|---|---:|---|
| `level` | `"O3"` | `"O1"`, `"O2"`, `"O3"`, or `"Osuper"` |
| `passes` | `None` | Explicit ordered pass names; replaces `level`'s pipeline |
| `fixpoint` | `False` | Run `O1` or an explicit pass pipeline until gate count stops decreasing |
| `decompose_rz` | `False` | Approximate `rz` operations with Clifford+T |
| `decompose_cz` | `False` | Lower `cz` to `h` + `cx` + `h` |
| `decompose_ccx` | `False` | Lower both `ccx` and `ccz` to Clifford+T |
| `rz_epsilon` | `1e-10` | Approximation tolerance used by `decompose_rz` |
| `parallel` | `False` | Enable the native parallel optimizer |
| `superopt_gates` | `"auto"` | MURM basis: `auto`, `base`, or an exact comma-separated gate list |

Explicit pass names are the same as the CLI's `--passes` values:
`DecomposeToffoli`, `DecomposeCz`, `DecomposeRz`, `CancelGates`, `SuperOpt`,
`PhaseFoldRand`, and `CnotMin`.

For example:

```python
optimized = optimize(
    circuit,
    passes=["CancelGates", "PhaseFoldRand"],
    fixpoint=True,
)
```

`ccx`, `ccz`, `cz`, and `rz` remain native by default. Requested
decompositions run between native and post-decomposition optimization stages.

An explicit `passes` list cannot be combined with `decompose_rz`,
`decompose_cz`, or `decompose_ccx`; put the corresponding pass directly in the list
instead. `O2` always runs two rounds, while `O3` and `Osuper` already run to a
fixpoint, so `fixpoint` only changes `O1` or an explicit pipeline.

## What the pass preserves

The adapter retains:

- quantum and classical registers and their bit identities;
- circuit name and metadata;
- the circuit's existing global phase;
- measurements, resets, and their classical-bit mappings.

Like tzap itself, circuit equivalence is considered up to global phase.
Independent operations may have a different valid topological order after a
Qiskit pass-manager round trip.

## Error handling

Adapter validation errors use `qiskit.transpiler.exceptions.TranspilerError`.
Native parsing or optimization failures use tzap's `QasmError` or
`OptimizationError`:

```python
from qiskit.transpiler.exceptions import TranspilerError
from tzap import OptimizationError, QasmError

try:
    optimized = optimize(circuit)
except TranspilerError as error:
    print(f"unsupported Qiskit circuit: {error}")
except (QasmError, OptimizationError) as error:
    print(f"tzap failed: {error}")
```

See the repository [README](../README.md) for installation alternatives and
CLI usage.
