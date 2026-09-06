import TzapLean.Equivalence

/-!
# Passes

The Rust `Pass` trait is a name and a `Circuit → Circuit` function:

```rust
pub trait Pass {
    fn name(&self) -> &str;
    fn run(&self, circuit: &Circuit) -> Circuit;
}
```

with "returns an equivalent one" a comment. Here the transformation acts directly on
`Circuit n m`, whose indices fix the register sizes and whose value carries the `Wf`
precondition. An unverified transformation is not a `Pass`, so composition does not need
separate preservation obligations.
-/

namespace TzapLean

noncomputable section

/-- An optimization pass: a circuit transformation carrying its own correctness proof. -/
structure Pass where
  /-- The pass's name, as in the Rust trait. -/
  name : String
  /-- The transformation on a validated, size-indexed circuit. -/
  run : ∀ {n m}, Circuit n m → Circuit n m
  /-- **The correctness obligation**: the output denotes the same channel as the input. -/
  correct : ∀ {n m} (c : Circuit n m),
    (run c).Equivalent c

namespace Pass

/-- Running one pass after another. -/
def comp (p q : Pass) : Pass where
  name := q.name ++ " ∘ " ++ p.name
  run := fun c => q.run (p.run c)
  correct c := Equivalent.trans (q.correct (p.run c)) (p.correct c)

/-- Run a list of passes in order, as the Rust `run_passes` does. -/
def runAll : List Pass → Circuit n m → Circuit n m
  | [], c => c
  | p :: ps, c => runAll ps (p.run c)

/-- **Composed correctness**: any pipeline of passes preserves the semantics. -/
theorem correct_runAll (ps : List Pass) (c : Circuit n m) :
    (runAll ps c).Equivalent c := by
  induction ps generalizing c with
  | nil => exact Equivalent.refl _ _ _
  | cons p ps ih =>
      exact Equivalent.trans (ih (p.run c)) (p.correct c)

end Pass

/-! ## Circuit statistics (the Rust `pass.rs` helpers) -/

/-- Number of `t`/`tdg` gates. -/
def countT (c : RawCircuit) : Nat :=
  c.gates.countP fun g => match g with | .t _ | .tdg _ => true | _ => false

/-- Number of `rz` gates. -/
def countRz (c : RawCircuit) : Nat :=
  c.gates.countP fun g => match g with | .rz _ _ => true | _ => false

/-- Number of two-qubit `cnot`/`cz` gates. -/
def count2q (c : RawCircuit) : Nat :=
  c.gates.countP fun g => match g with | .cnot .. | .cz .. => true | _ => false

/-- Per-wire next-free layer after scheduling `gs`, greedily as early as possible. -/
def depthAux : (Qubit → Nat) → List Gate → (Qubit → Nat)
  | next, [] => next
  | next, g :: gs =>
      let layer := (g.qubitsOf.map next).foldl max 0 + 1
      depthAux (fun q => if g.qubitsOf.contains q then layer else next q) gs

/-- Circuit depth: gates on disjoint qubits share a layer. -/
def depth (c : RawCircuit) : Nat :=
  ((List.range c.numQubits).map (depthAux (fun _ => 0) c.gates)).foldl max 0

end
end TzapLean
