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

with "returns an equivalent one" a comment. Here it is a field: a `Pass` *is* a function
together with proofs that it preserves the circuit's shape, its well-formedness, and — the
point of the exercise — its semantics. An unverified transformation is not a `Pass`, so
`runPasses` composes correctness for free.

The obligation is conditional on `Circuit.Wf` (distinct operands for multi-qubit gates,
see `Gate.Wf`): that is exactly the class of circuits the QASM front end produces, and
`cnot q q` really would break gate cancellation. Every pass must also *preserve* `Wf`,
which is what makes the conditional obligation compose.
-/

namespace TzapLean

noncomputable section

/-- An optimization pass: a circuit transformation carrying its own correctness proof. -/
structure Pass where
  /-- The pass's name, as in the Rust trait. -/
  name : String
  /-- The transformation. -/
  run : Circuit → Circuit
  /-- Passes never change the number of qubits. -/
  numQubits_run : ∀ c, (run c).numQubits = c.numQubits
  /-- Passes never change the number of classical bits. -/
  numCbits_run : ∀ c, (run c).numCbits = c.numCbits
  /-- Passes preserve well-formedness, so the obligation below composes. -/
  wf_run : ∀ c, c.Wf → (run c).Wf
  /-- Passes keep every operand in range: they never invent a wire the register does not
  have, so the output is a circuit the QASM back end is allowed to print. Conditional on `Wf`
  for the same reason `correct` is — `SuperOpt` reads the operands of the gates it absorbs
  into a window, and a repeated operand would not name the wire it appears to. -/
  wellFormed_run : ∀ c, c.Wf → c.WellFormed → (run c).WellFormed
  /-- Passes leave the cached `has*` flags describing the gates that came out, not the ones
  that went in. Every pass here rebuilds them outright, so the hypothesis is unused — it is
  there only so that the identity is a pass. -/
  flagsOk_run : ∀ c, c.FlagsOk → (run c).FlagsOk
  /-- **The correctness obligation**: the output denotes the same channel as the input. -/
  correct : ∀ c, c.Wf → Equivalent c.numQubits c.numCbits (run c).gates c.gates

namespace Pass

/-- Running one pass after another. -/
def comp (p q : Pass) : Pass where
  name := q.name ++ " ∘ " ++ p.name
  run := q.run ∘ p.run
  numQubits_run c := by simp [q.numQubits_run, p.numQubits_run]
  numCbits_run c := by simp [q.numCbits_run, p.numCbits_run]
  wf_run c hc := q.wf_run _ (p.wf_run c hc)
  wellFormed_run c hwf hc := q.wellFormed_run _ (p.wf_run c hwf) (p.wellFormed_run c hwf hc)
  flagsOk_run c hc := q.flagsOk_run _ (p.flagsOk_run c hc)
  correct c hc := by
    have h₁ : Equivalent c.numQubits c.numCbits (p.run c).gates c.gates := p.correct c hc
    have h₂ : Equivalent (p.run c).numQubits (p.run c).numCbits
        (q.run (p.run c)).gates (p.run c).gates := q.correct _ (p.wf_run c hc)
    rw [p.numQubits_run, p.numCbits_run] at h₂
    exact Equivalent.trans (by simpa using h₂) h₁

/-- Run a list of passes in order, as the Rust `run_passes` does. -/
def runAll : List Pass → Circuit → Circuit
  | [], c => c
  | p :: ps, c => runAll ps (p.run c)

theorem numQubits_runAll (ps : List Pass) (c : Circuit) :
    (runAll ps c).numQubits = c.numQubits := by
  induction ps generalizing c with
  | nil => rfl
  | cons p ps ih => simp [runAll, ih, p.numQubits_run]

theorem numCbits_runAll (ps : List Pass) (c : Circuit) :
    (runAll ps c).numCbits = c.numCbits := by
  induction ps generalizing c with
  | nil => rfl
  | cons p ps ih => simp [runAll, ih, p.numCbits_run]

theorem wf_runAll (ps : List Pass) (c : Circuit) (hc : c.Wf) : (runAll ps c).Wf := by
  induction ps generalizing c with
  | nil => exact hc
  | cons p ps ih => exact ih _ (p.wf_run c hc)

theorem wellFormed_runAll (ps : List Pass) (c : Circuit) (hwf : c.Wf) (hc : c.WellFormed) :
    (runAll ps c).WellFormed := by
  induction ps generalizing c with
  | nil => exact hc
  | cons p ps ih => exact ih _ (p.wf_run c hwf) (p.wellFormed_run c hwf hc)

theorem flagsOk_runAll (ps : List Pass) (c : Circuit) (hc : c.FlagsOk) :
    (runAll ps c).FlagsOk := by
  induction ps generalizing c with
  | nil => exact hc
  | cons p ps ih => exact ih _ (p.flagsOk_run c hc)

/-- **A pipeline of passes returns a circuit one can print and re-parse.** Together with
`Qasm.parse_valid`, which establishes the same three properties for whatever the front end
accepts, this closes the loop: the CLI can only ever write out a valid circuit. -/
theorem structural_runAll (ps : List Pass) (c : Circuit) (hwf : c.Wf) (hs : c.Structural) :
    (runAll ps c).Wf ∧ (runAll ps c).Structural :=
  ⟨wf_runAll ps c hwf, wellFormed_runAll ps c hwf hs.1, flagsOk_runAll ps c hs.2⟩

/-- **Composed correctness**: any pipeline of passes preserves the semantics. -/
theorem correct_runAll (ps : List Pass) (c : Circuit) (hc : c.Wf) :
    Equivalent c.numQubits c.numCbits (runAll ps c).gates c.gates := by
  induction ps generalizing c with
  | nil => exact Equivalent.refl _ _ _
  | cons p ps ih =>
      have h₁ : Equivalent c.numQubits c.numCbits (p.run c).gates c.gates := p.correct c hc
      have h₂ := ih (p.run c) (p.wf_run c hc)
      rw [p.numQubits_run, p.numCbits_run] at h₂
      exact Equivalent.trans h₂ h₁

end Pass

/-! ## Circuit statistics (the Rust `pass.rs` helpers) -/

/-- Number of `t`/`tdg` gates. -/
def countT (c : Circuit) : Nat :=
  c.gates.countP fun g => match g with | .t _ | .tdg _ => true | _ => false

/-- Number of `rz` gates. -/
def countRz (c : Circuit) : Nat :=
  c.gates.countP fun g => match g with | .rz _ _ => true | _ => false

/-- Number of two-qubit `cnot`/`cz` gates. -/
def count2q (c : Circuit) : Nat :=
  c.gates.countP fun g => match g with | .cnot .. | .cz .. => true | _ => false

/-- Per-wire next-free layer after scheduling `gs`, greedily as early as possible. -/
def depthAux : (Qubit → Nat) → List Gate → (Qubit → Nat)
  | next, [] => next
  | next, g :: gs =>
      let layer := (g.qubitsOf.map next).foldl max 0 + 1
      depthAux (fun q => if g.qubitsOf.contains q then layer else next q) gs

/-- Circuit depth: gates on disjoint qubits share a layer. -/
def depth (c : Circuit) : Nat :=
  ((List.range c.numQubits).map (depthAux (fun _ => 0) c.gates)).foldl max 0

end
end TzapLean
