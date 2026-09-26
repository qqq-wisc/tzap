# Worked example: why two pipelines give the same T count but different Pauli weight

This note walks through one 8-gate circuit, step by step. It shows how
tzap's gate-level optimization and the PBC rotation optimizer reach the same
T count but produce rotation axes of different total weight. Every
intermediate state below was produced by the tools, not derived by hand.

## The two pipelines

- **Pipeline A: tzap -O1, then conversion.** tzap optimizes the gate
  circuit at level O1 (phase folding and gate cancellation), and the result
  is converted to PBC:

  ```bash
  tzap input.qasm -O1 -o optimized.qasm
  tzap optimized.qasm --to-pbc --passes CancelGates -o out.pbc
  ```

- **Pipeline B: conversion, then "the pass".** The unoptimized circuit is
  converted to PBC, and then the PBC rotation optimizer runs on the
  converted circuit:

  ```bash
  tzap input.qasm --to-pbc --pbc-opt --passes CancelGates -o out.pbc
  ```

  (`--passes CancelGates` only stops tzap's default gate optimization from
  running first; on this circuit CancelGates changes nothing.)

**"The pass"** means the PBC rotation optimizer,
`PbcCircuit::optimize_rotations` in `src/pbc/optimize.rs`, run by
`--pbc-opt`. It works on the PBC form, not on gates. It walks the list of
Pauli rotations, merges rotations about the same axis when everything
between them commutes with that axis, and moves any Clifford rotation that
a merge produces into the output frame. It is described in
`docs/pbc-optimize.md`. Here it runs with its default settings
(`clifford_to_frame` on, `lazy_cliffords` off).

## Notation

A PBC program is a list of Pauli rotations followed by a Clifford frame:

- `r k s P` is the rotation `exp(-i · k·π/8 · s·P)` about the Pauli product
  P with sign s. Odd k costs one T gate (k = 1 is a T, k = -1 a T†). Even k
  is a Clifford: k = 2 is an S about P.
- `f Xq s P` and `f Zq s P` describe the output Clifford C, applied after
  all rotations, by its images `C† Xq C` and `C† Zq C`. Omitted lines map a
  Pauli to itself.
- The **weight** of an axis is its number of non-identity factors: `Z0 Z1`
  has weight 2, `Y0` has weight 1. The total weight sums over all rotations.

Conversion keeps the Clifford gates in a running frame instead of emitting
them. A T gate on qubit q becomes a rotation about the frame's current
image of Zq. That axis is Zq "pulled back" through every Clifford gate seen
so far. So an axis depends on where the Clifford gates are, and that is what
this example is about.

Two Pauli identities used below: `Z X = i Y` on one qubit, and `Z Z = I`.

## The circuit

```
cx q[0],q[1];
t q[1];
t q[1];
tdg q[1];
cx q[0],q[1];
h q[0];
cx q[1],q[0];
t q[0];
```

The three gates `t; t; tdg` on q1 sit between two CX gates on the same
qubits, so they act on the parity of q0 and q1, i.e. on the Pauli Z0Z1.
Together they equal a single T on that parity: 1 + 1 − 1 = 1.

## Step 1: converting the unoptimized circuit

Conversion processes the gates in order:

| gate | what happens | frame afterwards (non-trivial images) |
|---|---|---|
| `cx q[0],q[1]` | Clifford: update the frame | X0 → X0X1, Z1 → Z0Z1 |
| `t q[1]` | rotation about the image of Z1: `r 1 1 Z0 Z1` | unchanged |
| `t q[1]` | `r 1 1 Z0 Z1` | unchanged |
| `tdg q[1]` | `r -1 1 Z0 Z1` | unchanged |
| `cx q[0],q[1]` | undoes the first CX | identity |
| `h q[0]` | swaps X0 and Z0 | X0 → Z0, Z0 → X0 |
| `cx q[1],q[0]` | CX with control q1, target q0 | X0 → Z0, X1 → Z0X1, Z0 → X0Z1 |
| `t q[0]` | rotation about the image of Z0: `r 1 1 X0 Z1` | unchanged |

Result (no optimization at all):

```
r 1 1 Z0 Z1
r 1 1 Z0 Z1
r -1 1 Z0 Z1
r 1 1 X0 Z1
f X0 1 Z0
f X1 1 Z0 X1
f Z0 1 X0 Z1
```

T count 4, total weight 2 + 2 + 2 + 2 = 8. This is the starting point for
pipeline B.

## Step 2, pipeline A: tzap -O1 first

tzap's phase folding sees that the three T-type gates on q1 act on the same
parity and combines them: T · T · T† = T. It outputs

```
cx q[0],q[1]; t q[1]; cx q[0],q[1]; h q[0]; cx q[1],q[0]; t q[0];
```

No Clifford gate is created: the combined angle is odd, so a single `t`
expresses it. Converting this circuit follows the same frame updates as
step 1, with one T in place of three:

```
r 1 1 Z0 Z1
r 1 1 X0 Z1
f X0 1 Z0
f X1 1 Z0 X1
f Z0 1 X0 Z1
```

**T count 2, total weight 2 + 2 = 4.** Running the pass on this output
changes nothing: there is nothing left to merge (0 merges).

## Step 3, pipeline B: the pass on the converted circuit

The pass starts from the four rotations of step 1 and walks them in order.

**What F is.** Clifford rotations cost no T, so the pass does not keep them
in the rotation list. It moves them to the end of the program, where they
join the output Clifford (the `f` lines). F is the Clifford the pass is
currently moving toward the end. It starts as the identity.

**The one rule.** Moving a Clifford G later past a rotation about Q is
allowed, but it changes that rotation. Running G and then a rotation about
Q is the same as running a rotation about `G† Q G` and then G. For
G = S about P (`r 2` about P):

- if Q commutes with P, `G† Q G = Q`: the rotation is unchanged;
- if Q anticommutes with P, `G† Q G = i P Q`: the rotation gets a new
  axis.

Either way the rotation keeps its angle, so a T stays exactly one T.

### The walk, step by step

At every moment the program reads, left to right:

> processed rotations ; F ; rotations not yet processed ; C

where C is the output Clifford from conversion. F sits between the part
already processed and the rest. Processing the next rotation moves F one
position to the right, past that rotation. Every row describes the same
program, i.e. the same unitary.

| step | next input rotation | what the pass does | processed rotations | F | not yet processed |
|---|---|---|---|---|---|
| start | | | (none) | identity | T(Z0Z1), T(Z0Z1), T†(Z0Z1), T(X0Z1) |
| 1 | T(Z0Z1) | F is the identity, so the axis is unchanged. Nothing to merge with: keep it. | T(Z0Z1) | identity | T(Z0Z1), T†(Z0Z1), T(X0Z1) |
| 2 | T(Z0Z1) | Axis unchanged (F is still the identity). Same axis as the processed T(Z0Z1), with nothing in between: merge, 1 + 1 = 2. Angle 2 is S(Z0Z1), a Clifford: take it out of the list and make it F. | (none) | **S(Z0Z1)** | T†(Z0Z1), T(X0Z1) |
| 3 | T†(Z0Z1) | Move F past it. Z0Z1 commutes with Z0Z1, so the axis is unchanged. There is no earlier Z0Z1 rotation left to merge with (it became F): keep it. | T†(Z0Z1) | S(Z0Z1) | T(X0Z1) |
| 4 | T(X0Z1) | Move F past it. X0Z1 anticommutes with Z0Z1 (X0 vs Z0 anticommute, Z1 vs Z1 commute), so the axis becomes i·(Z0Z1)(X0Z1) = i·(Z0X0)(Z1Z1) = i·(iY0) = **−Y0**. The shared Z1 cancels: weight 2 → 1. Keep it. | T†(Z0Z1), T(−Y0) | S(Z0Z1) | (none) |
| end | | Nothing left: merge F into C. The new output Clifford is "S(Z0Z1), then C". | T†(Z0Z1), T(−Y0) | (in C) | (none) |

Changes to the rotations, compared with the input:

| input rotation | output rotation | why |
|---|---|---|
| T(Z0Z1) | removed | merged with the next T into S(Z0Z1), which moved into the output Clifford |
| T(Z0Z1) | removed | merged, as above |
| T†(Z0Z1) | T†(Z0Z1), unchanged | commutes with S(Z0Z1) |
| T(X0Z1) | T(−Y0), written `r -1 1 Y0` | anticommutes with S(Z0Z1): axis multiplied by Z0Z1 |

### Merging F into the output Clifford

The `f` lines give C's images `C† P C` for each single-qubit Pauli P. The
new output Clifford is S(Z0Z1) followed by C, so each image gets the same
rule applied (moved past S(Z0Z1)):

| image | before (C) | commutes with Z0Z1? | after (S(Z0Z1), then C) |
|---|---|---|---|
| X0 | Z0 | yes | Z0 |
| X1 | Z0 X1 | no (X1 vs Z1) | i·(Z0Z1)(Z0X1) = i·(Z1X1) = i·(iY1) = −Y1 |
| Z0 | X0 Z1 | no (X0 vs Z0) | i·(Z0Z1)(X0Z1) = −Y0 |
| Z1 | Z1 (implicit) | yes | Z1 (implicit) |

### Result

```
r -1 1 Z0 Z1      T†(Z0Z1)
r -1 1 Y0         T(−Y0)
f X0 1 Z0
f X1 -1 Y1
f Z0 -1 Y0
```

**T count 2, total weight 2 + 1 = 3.**

## Comparing the results

| | rotations | T | total weight |
|---|---|---|---|
| unoptimized | Z0Z1, Z0Z1, −Z0Z1, X0Z1 | 4 | 8 |
| A: tzap -O1, then conversion | Z0Z1, X0Z1 | 2 | 4 |
| B: conversion, then the pass | −Z0Z1, −Y0 | 2 | **3** |

Both outputs implement the same unitary, up to global phase. The test
suite checks the pass's rewrites exactly.

**Why B is narrower.** Both pipelines reduce the three T-type gates on Z0Z1
to one odd rotation. They do it differently:

- **tzap** adds the angles into one gate, 1 + 1 − 1 = 1, and leaves no
  Clifford behind.
- **The pass** merges the first two into a Clifford (1 + 1 = 2), moves it
  into the frame, and keeps the third as a separate T†. The Clifford in the
  frame then conjugates the last rotation.

Here the conjugation *narrows* the last axis, because X0Z1 and Z0Z1 overlap
on Z1 and that factor cancels.

## Verification

Every claim above is machine-checked by
`src/pbc/optimize/tests/weight_example.rs`, using the exact-arithmetic
reference semantics in `src/semantics/`. That code computes unitaries with
exact numbers of the form a + b√2 + i(c + d√2), with a, b, c, d rational,
and compares them up to a global phase, with no floating point. The tests
check that each of the following implements the same unitary as the input
circuit:

- **Step 1:** the input's direct conversion.
- **Pipeline A:** tzap -O1's output circuit, and its conversion.
- **Pipeline B:** the pass's output.
- **The walk table:** every row, rebuilt as the program *processed
  rotations ; F ; not yet processed ; C*, including the final merge of F
  into the output Clifford.
- **Each individual move of S(Z0Z1)** past an axis, as claimed in the
  tables: running S and then a rotation about Q equals running a rotation
  about Q' and then S.

They also pin the exact PBC text, T counts, and weights quoted here. As a
check that the checks can fail, deliberately wrong claims (the step-4 axis
as +Y0 instead of −Y0, or T†(Z0Z1) changing when moved past S) make the
tests fail.

## The general lesson

Moving a Clifford into the frame multiplies every later anticommuting axis
Q by the Clifford's axis P, up to phase. The weight of P·Q depends on how
the two overlap:

- **Disjoint supports widen.** Qubits appear in P·Q that were not in Q.
  This is what usually happens on large circuits, and why tzap's phase
  folding raises weight: its merges create S gates that conversion puts into
  the frame (see "Why tzap's optimizations raise weight" in
  `docs/pbc-optimize.md`).
- **Overlapping supports can narrow.** Shared factors cancel, as with Z1
  here.

So neither pipeline is always lighter. The reverse effect, where tzap's
output is wider, shows up in this 7-gate variant:

```
cx q[0],q[1]; t q[1]; t q[1]; t q[1]; cx q[0],q[1]; h q[0]; t q[0];
```

Here the three gates add up to T³, which tzap must write as `s; t`. The S
goes into the frame and widens the last axis from X0 to Y0Z1: T 2, weight
4. The pass with default settings does the same: its first merge also
leaves an S in the frame, so it also reaches weight 4. With
`lazy_cliffords: true`, the pass keeps a merged Clifford as a merge target
until something forces it into the frame. The third T then merges into it,
T³ stays one `r 3` rotation, no Clifford reaches the frame, and the last
axis stays X0: T 2, weight 3.

(In the 8-gate circuit above, lazy mode instead gives weight 4. Without the
S in the frame, nothing narrows the last axis. Which choice is lighter
depends on the overlaps.)
