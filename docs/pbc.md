# PBC exchange format

Export an optimized circuit with:

```bash
tzap input.qasm --to-pbc -o output.pbc
```

Conversion runs after optimization and requested decompositions. Gate-circuit
inputs must have no resets, and measurements must form a final block: once any
qubit is measured, only measurements may follow, on any qubit. For Rz gates,
also use `--decompose-rz`. CCX and CCZ convert natively, into seven rotations
each, and need no decomposition flag.

## Syntax

The first two lines declare the number of qubits and of classical registers,
each of which holds one bit:

```text
qubits 3
registers 2
```

Qubit IDs and register IDs start at zero. Multiple QASM `qreg`s or `creg`s are
numbered consecutively in declaration order, so with `creg a[1]; creg b[2];`,
`b[1]` becomes `c2`. The remaining lines contain operations in execution
order, followed by frame records:

```text
r <k> <sign> <Pauli factors>
m <sign> <Pauli factors> -> c<register ID>
frame <Xq or Zq> <sign> <Pauli factors>
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
- `frame` records specify the complete output Clifford C by giving its images
  `C† Xq C` and `C† Zq C`. They are not individual gates. The Clifford C acts
  after all `r` and `m` operations and is unique up to global phase.
  An omitted row means the identity image (`Xq` or `Zq` with sign `1`).
  Export writes changed X rows in increasing q order, then changed Z rows in
  increasing q order. The rows must collectively describe a valid Clifford
  action; their signs are `1` or `-1`.

For example, `r 1 -1 X0 Z2` rotates by pi/8 about `-X0 Z2` (equivalently, by
-pi/8 about `X0 Z2`), and
`m -1 Z1 -> c0` measures `-Z1`, reversing the outcome labels of a Z measurement.

## Output semantics

PBC conversion preserves the **full quantum–classical channel** of the circuit
being converted: the joint classical-result distribution, the conditional quantum
output states on every wire (including measured wires), and their correlations.
This holds for arbitrary input states, including inputs entangled with an external
system. Global phase is unobservable. Any approximation introduced earlier by
Rz decomposition is separate from this exact conversion.

Measurements are nondestructive projective measurements. Conversion conjugates
their axes by the accumulated Clifford frame, then retains the **entire** frame
as terminal `frame` records. Applying the represented Clifford restores the
quantum output states; discarding it generally preserves only classical
probabilities, not the full channel. No qubits or frame information are
automatically discarded, even when every input qubit is measured. Circuits
without measurements retain their frame too.

The terminal-measurement restriction applies to the gate-circuit input. The
output Clifford represented by the frame acts after measurements. Resets,
hidden measurement outcomes, and conditional operations are not supported by
this exchange format.

## Examples with terminal input measurements

These examples give the direct transformation; optimization may simplify a
circuit further. The exported frame preserves quantum outputs as well as the
classical results.

### H followed by measurement

Input: `H q0; measure q0 -> c0`.

```text
qubits 1
registers 1
m 1 X0 -> c0
frame X0 1 Z0
frame Z0 1 X0
```

H changes the measurement axis from Z to X; the frame represents that final H.

### X followed by measurement

Input: `X q0; measure q0 -> c0`.

```text
qubits 1
registers 1
m -1 Z0 -> c0
frame Z0 -1 Z0
```

X reverses the Z-measurement result; the frame represents that final X.

### A T rotation between two H gates

Input: `H q0; T q0; H q0; measure q0 -> c0`.

```text
qubits 1
registers 1
r 1 1 X0
m 1 Z0 -> c0
```

The T gate becomes an X-axis rotation, followed by Z measurement. The two
H gates cancel in the frame, so no `frame` records are needed.

### Bell preparation and readout

Input: `H q0; CX q0 -> q1; measure q0 -> c0; measure q1 -> c1`.

```text
qubits 2
registers 2
m 1 X0 -> c0
m 1 X0 Z1 -> c1
frame X0 1 Z0 X1
frame Z0 1 X0
frame Z1 1 X0 Z1
```

For input `|00>`, the classical outputs are `00` and `11`, each with probability
one half. The second instruction measures the joint product `X0 Z1`.
The frame restores the corresponding quantum outputs `|00>` and `|11>`.

### Partial readout

Input: `H q0; CX q0 -> q1; measure q1 -> c0` (one classical register).

```text
qubits 2
registers 1
m 1 X0 Z1 -> c0
frame X0 1 Z0 X1
frame Z0 1 X0
frame Z1 1 X0 Z1
```

The frame represents both Cliffords, including the entangling CX. The state on unmeasured q0,
the state on measured q1, and their correlations with c0 are preserved.
