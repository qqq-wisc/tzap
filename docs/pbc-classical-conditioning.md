# Classical conditioning in tzap's PBC pipeline

This note lists every place where classical conditioning (an operation that
depends on a measurement outcome) arises around PBC conversion. For each case
it records what tzap does today and how the case could be supported.

Throughout, `C` is the postponed output Clifford. At every point in the
conversion, the gate prefix equals `C · (emitted PBC operations)`, and the
converter tracks the images `C† Xq C` and `C† Zq C`. Rotations use
`R(P, k) = exp(-i·kπ/8·P)`: T is `k = 1`, S is `k = 2`, and a Pauli is `k = 4`.

## Summary

| # | Source | Classical conditioning | Today |
|---|--------|------------------------|-------|
| 1 | Mid-circuit measurement | None (measurement only) | **Supported** |
| 2 | Reset | Conditional Pauli on a hidden outcome | Rejected by `to_pbc` |
| 3 | QASM `if (c==n) gate;` | Conditional gate on a register value | Not parsed |
| 4 | Measurement-based uncomputation | Conditional Clifford (e.g. CZ) | Not produced |
| 5 | Teleportation and syndrome decoding | Conditional Pauli on a parity of outcomes | Not representable |
| 6 | Commuting conditional Paulis to the end | Conditional rotation signs, flipped outcomes | Not done |
| 7 | Magic-state injection of π/8 rotations | Conditional π/4 correction | Out of scope (logical level) |

What the PBC IR already has (`src/pbc/mod.rs`):

- **Hidden measurements.** `measure(axis, None)` records an outcome with an
  immutable `MeasId` but writes no classical register.
- **Conditional rotations.** `PbcOp::ConditionalRotate { axis, angle, if_one }`
  and `conditional_pauli(axis, m)`, which is `k = 4`. They apply only when the
  outcome `if_one` is 1. The condition is a single outcome, not a register
  value or a parity.
- **ASCII rendering.** Both are drawn, for example `R(pi/2,+) if m0=1`.

What rejects them:

- **`to_pbc`.** Rejects `reset`. (Gates after measurements are accepted;
  see case 1.)
- **Text format.** `to_text` returns `UnsupportedTextOperation` for hidden
  measurements and conditional rotations.
- **Exact channel oracle.** `pbc_channel` supports hidden measurements, since
  their histories are summed incoherently, but rejects `ConditionalRotate`.
  `circuit_channel` rejects resets.
- **QASM parser and gate IR.** They have no `if` and no conditional gate.

The optimizer treats `measure` and `reset` as barriers on their own qubits
only (`src/cancel.rs`). Gates on other qubits may move across a measurement,
which is correct because they commute with it. The extended fuzzer
`fuzz_mid_circuit_measurements` checks this: it optimizes random mid-circuit
circuits at every level, sequentially and in parallel, then converts them,
and compares exact channels.

**Key invariant for every future case: conditioning never changes the output
frame.** A conditional operation is emitted as conditional rotations about the
current images, and `C` stays unconditional. Conversion stays linear, and the
terminal `f` records keep their meaning.

---

## 1. Mid-circuit measurement

**Conditioning:** none. Measurements in PBC are nondestructive, so later gates
simply continue to act on the post-measurement state.

**Today: supported.** `to_pbc` emits `measure(C† Zq C, Some(cbit))`, leaves
the frame unchanged, and continues converting. The text format expresses it
directly, because `m` lines may appear anywhere among the `r` lines, and
`circuit_channel` evaluates gates after measurements on each branch. See the
mid-circuit examples in `docs/pbc.md`.

**Future optimization impact:** a measurement of `P` blocks merging or
commuting a later rotation about `P'` past it when `P'` anticommutes with
`P`. Rotations that commute with `P` can pass.

## 2. Reset

**Conditioning:** reset = measure Z, then apply X if the outcome is 1. The
outcome is internal, not a user-visible bit.

**Today:** the parser accepts `reset`, the optimizer treats it as a per-qubit
barrier, and `to_pbc` rejects it with `UnsupportedGate(Reset)`.

**Future (conversion):**

```
m    = measure(C† Zq C, None)         // hidden outcome
conditional_pauli(C† Xq C, m)         // R(C† Xq C, 4) if m = 1
```

The two images anticommute, so this prepares the +1 eigenstate of `C† Zq C`.
After the suffix `C`, qubit q is in |0>. The frame is unchanged. The cost is
two operations, with no expansion.

**Future (format):** there are two equivalent choices, which is Abtin's point.

- **Expose the conditional.** Add hidden outcome IDs and an `if` clause
  (see [Text format proposal](#text-format-proposal)).
- **Reset as a primitive.** Add a line `p <P> ; <Q>`: "prepare the +1
  eigenstate of P, using the anticommuting Q as the flip". It hides the
  conditional exactly as a logical π/8 rotation hides its injection correction
  (case 7).

Either way, the channel oracle needs `ConditionalRotate` support (see
[Semantics and testing](#semantics-and-testing)).

**Policy question:** QASM exporters such as Qiskit often emit leading resets.
Under tzap's semantics (arbitrary input state), a leading reset is not a
no-op. Dropping it is valid only under a "|0…0> input" assumption, and that
should be an opt-in flag, not the default.

## 3. Classically conditioned gates (QASM `if`)

**Conditioning:** OpenQASM 2 `if (creg == n) gate;` compares the *whole
register* to an integer at that point in the program.

**Today:** not parsed, and the gate IR has no conditional gate.

**Future (conversion):** a gate `G` conditioned on `c` satisfies
`G^c · C · P = C · (C† G C)^c · P`. So emit the rotations that implement `G`,
about the current images, each with the same condition. The frame is
unchanged.

- **Pauli X/Z:** one conditional `R(image, 4)`.
- **S/Sdg:** conditional `R(C† Zq C, ±2)`.
- **T/Tdg:** conditional `R(C† Zq C, ±1)`.
- **H:** `R(Z,2) R(X,2) R(Z,2)` on the images, all conditional.
- **CZ:** `R(Za,2) R(Zb,2) R(ZaZb,-2)`.
- **CX:** `R(Zc,2) R(Xt,2) R(ZcXt,-2)`.
- **CCX/CCZ:** the seven native rotations, all conditional.

**Needed IR change:** `if_one: MeasId` is not enough. A QASM condition is a
conjunction over the bits of a register. Each bit is either the latest outcome
written to it, or its initial value if it was never written. Options:

- Generalize the condition to a conjunction of literals over `MeasId`s and
  initial bits.
- Assume the OpenQASM 2 initial value of 0. Then unwritten bits are constants
  and fold away statically. This is consistent with QASM, but narrower than
  the oracle's arbitrary initial classical store.

**Optimizer impact:** every pass would need to treat conditioned gates, at
minimum as barriers on their qubits. Alternatively, run `to_pbc` on the
unconditioned segments only. Adding `if` to `Circuit` touches every pass, so
it is the most invasive item in this note.

## 4. Measurement-based uncomputation

**Conditioning:** Gidney's AND/Toffoli uncomputation measures the target in
the X basis and applies CZ to the controls if the outcome is 1. This saves
all of the uncomputation's T gates.

**Today:** tzap neither produces nor accepts it.

**Future:** once cases 2 and 3 exist, it is a hidden X-basis measurement
followed by a conditional CZ (case 3's three conditional π/4 rotations). It
is also an *optimization opportunity*: a pass could recognize compute/uncompute
Toffoli pairs and rewrite the uncompute this way, trading T gates for
measurements and conditional Cliffords.

## 5. Parity-conditioned Paulis (teleportation, syndromes)

**Conditioning:** corrections such as `X^{m2} Z^{m1}`, or Pauli corrections
from a decoder, depend on *XORs* of several outcomes.

**Today:** not representable, because a condition names exactly one `MeasId`.

**Future:** use a parity (a set of `MeasId`s) as the condition type for
conditional Paulis. Conditional Paulis with parity conditions form a group:
two conditional Paulis on the same axis combine by XOR-ing their conditions,
so they compose and cancel without ever expanding into branches. This is also
the natural representation for case 6.

## 6. Commuting conditional Paulis to the end (Pauli frame tracking)

This is the origin of the "conditional rotation layers" in the question.

**Conditioning:** moving a conditional Pauli `Q^m` later in the program:

- **Past a rotation `R(P, k)` that anticommutes with `Q`:** the rotation
  becomes `R(P, (-1)^m k)`, so later rotation signs depend on outcomes. For
  `k = ±1`, the sign flip is a conditional `R(P, ∓2)`, a conditional Clifford.
- **Past a measurement of `P` that anticommutes with `Q`:** the outcome flips.
  Its label becomes `outcome XOR m`. This change is purely classical and free.
- **Past a commuting operation:** nothing changes.

**Today:** not done. Nothing produces conditional Paulis yet.

**Future:** the default should be to **leave conditional Paulis in place**.
This is exact, linear, and matches hardware that tracks the Pauli frame in
software. Commuting them to the end is worthwhile only for a specific goal,
such as reducing non-Clifford depth. When it is done:

- **Measurements:** fold the flips into outcome parities (case 5).
- **Rotations:** materialize the flips as conditional π/4 rotations, or keep
  a per-rotation sign-parity annotation.

Either choice should be a separate, optional pass, measured apart from the
linear converter, as `docs/pbc-todo.md` already requires for optimizations.

## 7. Magic-state injection of π/8 rotations

**Conditioning:** on fault-tolerant hardware, `R(P, ±1)` is applied by
injection:

1. Measure `P ⊗ Z_anc` jointly with a magic-state ancilla.
2. Measure `X_anc`.
3. Apply a π/4 correction `R(P, ±2)` conditioned on the outcomes.

In Litinski's auto-corrected variant, the correction becomes a conditional
*choice of measurement basis* for a second ancilla (X or Y). Either way, every
π/8 rotation carries classical conditioning with it.

**Today:** out of scope. tzap emits *logical* PBC, in which π/8 rotations are
primitive.

**Future:** this applies only if tzap gains a lowering to injection-level PBC
or to lattice surgery. The conditions would be outcome parities (case 5), and
the correction π/4 rotations would be conditional rotations, or,
equivalently, Clifford updates to a conditional Pauli frame (case 6).

**Consistency argument (Abtin):** treating logical π/8 rotations as primitives
already hides a conditional operation. The format can hide reset's
conditional Pauli in the same way (case 2, "reset as a primitive"). So there
is no reason to support one and not the other.

---

## Text format proposal

This is a minimal extension that covers cases 1–5 and keeps today's files
valid.

```text
m <sign> <factors> [-> c<k>]          # every m gets implicit ID m0, m1, ... in order
r <k> <sign> <factors> [if <cond>]
<cond> := <lit> | <cond> ^ <cond>     # parity of outcomes
<lit>  := m<id> | !m<id>
```

- **Hidden measurements:** an `m` line without `-> c<k>`.
- **Reset:** a hidden `m` followed by `r 4 <sign> <Q> if m<id>`.
- **Register conditions (case 3):** resolved at export into outcome literals
  and constants. A conjunction of several bits needs either one conditional
  rotation per satisfying branch, or a richer `<cond>` that allows `&`.
- **`frame` records:** they stay terminal and unconditional.

## Semantics and testing

**Oracle:** `pbc_channel` must evaluate conditional rotations. Each branch
carries its outcome history, as a vector indexed by `MeasId`, and applies the
rotation when the condition holds. Hidden outcomes are then summed as they are
today. `circuit_channel` needs reset, gates after measurements, and
conditioned gates, all driven by the explicit initial classical store it
already takes. The branch budget stays `2^measurements`, so exact tests are
limited to about 6 measurements.

**Tests to add with each case:**

- **Reset:** reset agrees with measure + conditional X, including on inputs
  entangled with a reference system.
- **Mid-circuit measurement:** for each Clifford probe, a mid-circuit
  measurement followed by more gates matches the gate circuit's channel.
- **Conditioned gates:** `if` on an unwritten register depends on the initial
  store, and on a written one depends on the latest write.
- **Parity conditions:** conditional Paulis compose by XOR.
- **Fuzzing:** extend the seeded export fuzz (`src/pbc/text/tests.rs`) to
  random resets and conditionals.

## Suggested order

1. ~~**Mid-circuit measurements (case 1).**~~ Done.
2. **Reset (case 2) plus `ConditionalRotate` in the oracle and text format.**
   This settles the "expose vs. primitive" question.
3. **Parity conditions (case 5).** This is the prerequisite for cases 6 and 7.
4. **QASM `if` (case 3).** This is the largest change, because it touches
   the parser, the gate IR, and every optimization pass.
5. **Measurement-based uncomputation (case 4)** and **Pauli frame tracking
   (case 6)**, as optional optimization passes.
6. **Injection-level lowering (case 7)**, only if tzap targets fault-tolerant
   PBC below the logical level.
