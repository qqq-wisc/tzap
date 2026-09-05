# tzap-lean

`tzap-lean` is a Lean 4 port of [tzap](https://github.com/qqq-wisc/tzap), the Rust quantum-circuit optimizer. It formalizes the circuit representation and channel semantics, then implements optimizer passes together with proofs that they preserve those semantics.

Every executable optimizer pass is a `Pass`, so it carries its structural and semantic proof obligations. The verified executable core composes the selected passes and fixpoint loop; `runConfigured_correct` proves its result equivalent to its input. Checked serialization reparses generated OpenQASM before it is written, and `runConfigured_checkedOutput_correct` connects accepted input, optimization, and successful output emission.

```sh
lake exe cache get
lake build
```

## The obligation

```lean
structure Pass where
  name : String
  run : Circuit → Circuit
  numQubits_run : ∀ c, (run c).numQubits = c.numQubits
  numCbits_run : ∀ c, (run c).numCbits = c.numCbits
  wf_run : ∀ c, c.Wf → (run c).Wf
  wellFormed_run : ∀ c, c.Wf → c.WellFormed → (run c).WellFormed
  flagsOk_run : ∀ c, c.FlagsOk → (run c).FlagsOk
  correct : ∀ c, c.Wf → Equivalent c.numQubits c.numCbits (run c).gates c.gates
```

`Pass.comp` and `Pass.runAll` compose that obligation, so any pipeline of passes is correct by construction (`Pass.correct_runAll`).

The final two structural fields ensure that output remains serializable: `WellFormed` means every operand is in range, and `FlagsOk` means cached circuit flags match the gates. `Qasm.parse_valid` establishes these conditions, together with `Wf`, for accepted OpenQASM input.

`Circuit.Wf` requires distinct operands in every multi-qubit gate. This is semantically necessary: `cnot q q` is idempotent rather than self-inverse, so cancelling two such gates would be unsound.

## What is trusted

| Component | Role in the trust boundary | Consequence if wrong |
|---|---|---|
| Formal circuit semantics and specification | Define what “correct” means: equality of the modeled quantum channels. | The development could prove preservation of an unintended model. |
| OpenQASM output boundary | `serializeChecked` reparses emitted text and refuses to write it unless it reconstructs the optimized circuit. | A failed check stops output rather than emitting a changed circuit. |
| Lean compiler and runtime | Execute the proved definitions and IO shell. | Outside the logic’s proof guarantee, as for any compiled verification artifact. |
