# tzap API

## Circuits

A `Circuit` holds a list of gates over a fixed number of qubits.

```rust
use tzap::circuit::{Circuit, Gate};

let mut circuit = Circuit::new(2);
circuit.apply(Gate::h(0));
circuit.apply(Gate::cnot { control: 0, target: 1 });
circuit.apply(Gate::t(0));
```

### Supported gates

| Gate | Constructor |
|------|------------|
| X | `Gate::x(qubit)` |
| H | `Gate::h(qubit)` |
| S | `Gate::s(qubit)` |
| Sdg | `Gate::sdg(qubit)` |
| Z | `Gate::z(qubit)` |
| T | `Gate::t(qubit)` |
| Tdg | `Gate::tdg(qubit)` |
| Rz | `Gate::rz(angle, qubit)` |
| CNOT | `Gate::cnot { control, target }` |
| Toffoli | `Gate::ccx { control1, control2, target }` |

### QASM I/O

Parse from and convert to OpenQASM 2.0:

```rust
let circuit = Circuit::from_qasm("
    OPENQASM 2.0;
    include \"qelib1.inc\";
    qreg q[2];
    h q[0];
    cx q[0],q[1];
");

let qasm_string = circuit.to_qasm();
```

## Passes

Every pass implements the `Pass` trait and returns a new circuit:

```rust
use tzap::pass::Pass;
```

### Available passes

| Pass | Import | Description |
|------|--------|-------------|
| `DecomposeToffoli` | `tzap::decompose` | Breaks Toffoli gates into CNOT+T/Tdg |
| `DecomposeRz` | `tzap::decompose_rz` | Decomposes Rz gates into Clifford+T via gridsynth |
| `CancelPairs` | `tzap::cancel` | Removes adjacent self-inverse gate pairs (HH, XX, etc.) |
| `PhaseFoldGlobal` | `tzap::phase_fold_global` | Merges T/Rz gates across the circuit via symbolic parity tracking |
| `PhaseFoldGlobalExpr` | `tzap::phase_fold_global_expr` | Like `PhaseFoldGlobal` but uses exact symbolic expressions instead of hashes |

### Running passes

Run a single pass:

```rust
use tzap::decompose::DecomposeToffoli;

let optimized = DecomposeToffoli.run(&circuit);
```

Run a pipeline:

```rust
use tzap::decompose::DecomposeToffoli;
use tzap::cancel::CancelPairs;
use tzap::phase_fold_global::PhaseFoldGlobal;
use tzap::pass::run_passes;

let passes: Vec<&dyn Pass> = vec![
    &DecomposeToffoli,
    &CancelPairs,
    &PhaseFoldGlobal,
];

let result = run_passes(&circuit, &passes);
println!("{} gates, {} T", result.circuit.gates.len(), result.circuit.gates.iter()
    .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count());
```

### DecomposeRz epsilon

Control the approximation precision with the `epsilon` field (default `1e-10`):

```rust
use tzap::decompose_rz::DecomposeRz;

let pass = DecomposeRz { epsilon: 1e-6 };
let cliffordt = pass.run(&circuit);
```
