# tzap-lean

`tzap-lean` is a Lean 4 port of [tzap](https://github.com/qqq-wisc/tzap), the Rust quantum-circuit optimizer. It formalizes the circuit representation and channel semantics, then implements optimizer passes together with proofs that they preserve those semantics.

Every executable optimizer pass is a `Pass`. Its `run` function takes an indexed `Circuit.Checked n m`, which fixes the register sizes and carries the operand-distinctness invariant needed by the proofs. The verified executable core composes the selected passes and fixpoint loop; `runConfigured_correct` proves its result equivalent to its input. Checked serialization reparses generated OpenQASM before it is written, and `runConfigured_checkedOutput_correct` connects accepted input, optimization, and successful output emission.

## Installation

Install [Lean through `elan`](https://lean-lang.org/install/), then clone and build the project:

```sh
git clone https://github.com/qqq-wisc/tzap.git
cd tzap/lean
lake exe cache get
lake build
```

## Running

Optimize an OpenQASM 2.0 circuit and write the result to a file:

```sh
lake exe tzap-lean input.qasm optimized.qasm
```

## The obligation

```lean
structure Pass where
  name : String
  run : ∀ {n m}, Circuit.Checked n m → Circuit.Checked n m
  correct : ∀ {n m} (c : Circuit.Checked n m),
    (run c).Equivalent c
```

`Pass.comp` and `Pass.runAll` compose that obligation, so any pipeline of passes is correct by construction (`Pass.correct_runAll`).

Raw circuits remain parser-facing so malformed QASM can receive diagnostics. `Circuit.Checked` is the optimizer boundary: its indices prevent register-size changes and its `Wf` field prevents repeated multi-qubit operands. `WellFormed` and `FlagsOk` are deliberately not part of semantic pass correctness; `serializeChecked` remains the fail-closed output boundary.

`Circuit.Wf` requires distinct operands in every multi-qubit gate. This is semantically necessary: `cnot q q` is idempotent rather than self-inverse, so cancelling two such gates would be unsound.

## What is trusted

| Component | Role in the trust boundary | Consequence if wrong |
|---|---|---|
| Formal circuit semantics and specification | Define what “correct” means: equality of the modeled quantum channels. | The development could prove preservation of an unintended model. |
| OpenQASM output boundary | `serializeChecked` reparses emitted text and refuses to write it unless it reconstructs the optimized circuit. | A failed check stops output rather than emitting a changed circuit. |
| Lean compiler and runtime | Execute the proved definitions and IO shell. | Outside the logic’s proof guarantee, as for any compiled verification artifact. |
