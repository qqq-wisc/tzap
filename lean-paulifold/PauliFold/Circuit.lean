import TzapLean.Semantics

/-!
# A Clifford+T circuit model and its concrete domain

Circuits are lists of `CTGate`s on `n` qubits, with qubits as `Fin n` and two-qubit gates on
distinct wires. Each gate maps to a `TzapLean.Gate`, so the concrete semantics is the one of
`TzapLean/Semantics.lean`.

The *concrete domain* of rotation folding is the operator being transported: a rotation
about `P` moved forward through a gate `U` becomes a rotation about `U P U†` (the transport
lemma of `docs/pauli-axis-abstract-interpretation.md`, Section 1). Concrete values are
therefore `2ⁿ × 2ⁿ` matrices (`TzapLean.Density n`), and a gate's concrete transformer is
conjugation by its unitary.
-/

namespace PauliFold

open TzapLean

/-- The Clifford+T gates: the Cliffords H, X, Z, S, S†, CX, CZ and the rotations T, T†.
Two-qubit gates carry a proof that their operands differ. -/
inductive CTGate (n : Nat) where
  | h (q : Fin n)
  | x (q : Fin n)
  | z (q : Fin n)
  | s (q : Fin n)
  | sdg (q : Fin n)
  | t (q : Fin n)
  | tdg (q : Fin n)
  | cx (c t : Fin n) (hne : c ≠ t)
  | cz (c t : Fin n) (hne : c ≠ t)

namespace CTGate

variable {n : Nat}

/-- The corresponding gate of `TzapLean`. -/
def toGate : CTGate n → Gate
  | .h q => .h q.val
  | .x q => .x q.val
  | .z q => .z q.val
  | .s q => .s q.val
  | .sdg q => .sdg q.val
  | .t q => .t q.val
  | .tdg q => .tdg q.val
  | .cx ctl tgt _ => .cnot ctl.val tgt.val
  | .cz ctl tgt _ => .cz ctl.val tgt.val

end CTGate

/-- Concrete values: operators on `n` qubits. -/
abbrev Op (n : Nat) := Density n

/-- The concrete transformer of a gate: `A ↦ U A U†`. -/
noncomputable def CTGate.concrete {n : Nat} (g : CTGate n) (A : Op n) : Op n :=
  conj (gateUnitary n g.toGate) A

/-- The concrete transformer of a gate list, first gate first. -/
noncomputable def concreteRun {n : Nat} (gs : List (CTGate n)) (A : Op n) : Op n :=
  gs.foldl (fun A g => g.concrete A) A

end PauliFold
