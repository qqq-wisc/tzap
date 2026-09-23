# PBC exchange format

Export an optimized circuit with:

```bash
tzap input.qasm --to-pbc -o output.pbc
```

Conversion runs after optimization and requested decompositions. Gate-circuit
inputs must have no resets and measurements only at the end. For Rz gates, also
use `--decompose-rz`.

## Syntax

The first two lines declare the number of qubits and classical registers:

```text
qubits 3
registers 2
```

Qubit IDs and register IDs start at zero. The remaining lines are instructions,
executed from top to bottom:

```text
r <k> <sign> <Pauli factors>
m <sign> <Pauli factors> -> c<register ID>
```

- `r` applies `exp(-i * k*pi/8 * P)`, where `k` is an integer and `P` is the
  signed Pauli product. The exporter uses `k` from -3 through 4; angles differing
  by 8 are equivalent up to global phase.
- `m` measures `P`: eigenvalue +1 writes bit 0; eigenvalue -1 writes bit 1.
- The sign is `1` or `-1`. Signs `i` and `-i` are not allowed: rotation and
  measurement axes must be Hermitian.
- Factors are `X0`, `Y1`, `Z2`, etc., listed in increasing qubit-ID order, with
  each qubit appearing at most once. Omitted qubits carry identity; an empty
  factor list means identity.
- All factors on a line form **one joint operation**, not separate operations.
- Measurements write directly to registers such as `c0`. A later write to the
  same register overwrites its value. There are no separate measurement IDs.

For example, `r 1 -1 X0 Z2` rotates by pi/8 about `-X0 Z2`, and
`m -1 Z1 -> c0` measures `-Z1`, reversing the outcome labels of a Z measurement.

## Output semantics

When every input qubit is measured at the end, CLI export preserves the classical
result distribution and discards quantum outputs. The remaining Clifford frame
is absorbed into measurement axes and signs; no final Clifford operations are
needed.

Otherwise, the current CLI preserves **all** quantum outputs, including measured
qubits. Remaining Cliffords become ordinary `r` instructions, which may follow
measurements to restore the quantum output state. There is no `suffix` section
and no separate syntax for H, CX, or other Clifford gates. Discarding only the
measured qubits while retaining unmeasured quantum outputs is not currently an
automatic simplification.

## Examples ending in measurement

These examples give the direct transformation; optimization may simplify a
circuit further. Every input qubit is measured, so only classical outputs matter.

### H followed by measurement

Input: `H q0; measure q0 -> c0`.

```text
qubits 1
registers 1
m 1 X0 -> c0
```

H changes the measurement axis from Z to X.

### X followed by measurement

Input: `X q0; measure q0 -> c0`.

```text
qubits 1
registers 1
m -1 Z0 -> c0
```

X reverses the Z-measurement result.

### A T rotation between two H gates

Input: `H q0; T q0; H q0; measure q0 -> c0`.

```text
qubits 1
registers 1
r 1 1 X0
m 1 Z0 -> c0
```

The T gate becomes an X-axis rotation, followed by Z measurement.

### Bell preparation and readout

Input: `H q0; CX q0 -> q1; measure q0 -> c0; measure q1 -> c1`.

```text
qubits 2
registers 2
m 1 X0 -> c0
m 1 X0 Z1 -> c1
```

For input `|00>`, the classical outputs are `00` and `11`, each with probability
one half. The second instruction measures the joint product `X0 Z1`.

## Limits

This is an expanded exchange format, not the converter's shared internal
representation. Export can exceed linear size and currently has a budget of
16 million expansion cells (arena nodes times qubit count, at least one cell per
node). Exceeding the budget returns an error rather than truncating the circuit.
The Rust API exposes this budget through `TextOptions`.
