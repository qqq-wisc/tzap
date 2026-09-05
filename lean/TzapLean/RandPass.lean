import Mathlib.Probability.ProbabilityMassFunction.Constructions
import TzapLean.Pass

/-!
# Randomized Passes

`Pass` demands unconditional correctness, which a randomized optimizer cannot sign: fixing a
seed fixes the transformation, and no fixed seed can be right on *every* circuit. A pass
whose parities are `k`-bit random tags is already defeated by a circuit on `k+1` wires,
which has more distinct parities than there are tags.

So the primitive here is `RandPass`: a transformation together with a distribution on its
seed and a *bound on the probability that its output is wrong*.

```
correct : ∀ c, c.Wf → Pr_{s ← dist c} [ ⟦run c s⟧ ≠ ⟦c⟧ ] ≤ error c
```

A deterministic pass is the `error = 0` case, with a one-point seed: `Pass.toRand` turns any
`Pass` into a `RandPass` whose failure event is literally empty. Nothing about the existing
`CancelGates` and `CnotMin` proofs changes; they are reused verbatim.

Composition is where the design pays off. `RandPass.comp` draws the second pass's seed
*after* seeing the first pass's output — `PMF.bind`, so the seed space of the composite is a
sigma type and no independence argument is needed. The error adds:

```
error (comp p q) c = p.error c + ⨆ s, q.error (p.run c s)
```

and a pipeline of deterministic passes collapses to `0`, recovering `Pass.correct_runAll`.
Note that the bound needs no independence *between* the two failure events — it is a union
bound, which holds regardless.
-/

namespace TzapLean

open scoped ENNReal

noncomputable section

/-- An optimization pass that may consult randomness, carrying a bound on how often its
output can fail to denote the same channel as its input. -/
structure RandPass where
  /-- The pass's name. -/
  name : String
  /-- The randomness the pass consumes, as a function of the circuit it is given. -/
  Seed : Circuit → Type
  /-- The distribution the seed is drawn from (`PMF.uniformOfFintype` in practice). -/
  dist : (c : Circuit) → PMF (Seed c)
  /-- The transformation, for a given seed. -/
  run : (c : Circuit) → Seed c → Circuit
  /-- The failure probability this pass is allowed. -/
  error : Circuit → ℝ≥0∞
  /-- Passes never change the number of qubits. -/
  numQubits_run : ∀ c s, (run c s).numQubits = c.numQubits
  /-- Passes never change the number of classical bits. -/
  numCbits_run : ∀ c s, (run c s).numCbits = c.numCbits
  /-- Passes preserve well-formedness, whatever the seed. -/
  wf_run : ∀ c s, c.Wf → (run c s).Wf
  /-- Passes keep every operand in range, whatever the seed. -/
  wellFormed_run : ∀ c s, c.Wf → c.WellFormed → (run c s).WellFormed
  /-- Passes leave the cached `has*` flags describing the gates that came out. -/
  flagsOk_run : ∀ c s, c.FlagsOk → (run c s).FlagsOk
  /-- **The correctness obligation**: the output denotes the same channel as the input,
  except on a set of seeds of probability at most `error c`. -/
  correct : ∀ c, c.Wf →
    (dist c).toOuterMeasure
        {s | ¬ Equivalent c.numQubits c.numCbits (run c s).gates c.gates} ≤ error c

namespace RandPass

/-- The failure event of a pass on a circuit. -/
def failure (p : RandPass) (c : Circuit) : Set (p.Seed c) :=
  {s | ¬ Equivalent c.numQubits c.numCbits (p.run c s).gates c.gates}

/-! ## Deterministic passes are the `error = 0` case -/

/-- Any `Pass` is a `RandPass` that ignores its seed and never fails. -/
def _root_.TzapLean.Pass.toRand (p : Pass) : RandPass where
  name := p.name
  Seed := fun _ => Unit
  dist := fun _ => PMF.pure ()
  run := fun c _ => p.run c
  error := fun _ => 0
  numQubits_run c _ := p.numQubits_run c
  numCbits_run c _ := p.numCbits_run c
  wf_run c _ hc := p.wf_run c hc
  wellFormed_run c _ hwf hc := p.wellFormed_run c hwf hc
  flagsOk_run c _ hc := p.flagsOk_run c hc
  correct c hc := by
    have hempty : {s : Unit | ¬ Equivalent c.numQubits c.numCbits (p.run c).gates c.gates} = ∅ := by
      ext s
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact p.correct c hc
    rw [hempty]
    simp

@[simp] theorem toRand_error (p : Pass) (c : Circuit) : (Pass.toRand p).error c = 0 := rfl

@[simp] theorem toRand_run (p : Pass) (c : Circuit) (s : Unit) :
    (Pass.toRand p).run c s = p.run c := rfl

/-- The identity pass. -/
def id : RandPass where
  name := "id"
  Seed := fun _ => Unit
  dist := fun _ => PMF.pure ()
  run := fun c _ => c
  error := fun _ => 0
  numQubits_run _ _ := rfl
  numCbits_run _ _ := rfl
  wf_run _ _ hc := hc
  wellFormed_run _ _ _ hc := hc
  flagsOk_run _ _ hc := hc
  correct c _ := by
    have hempty : {s : Unit | ¬ Equivalent c.numQubits c.numCbits c.gates c.gates} = ∅ := by
      ext s
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact Equivalent.refl _ _ _
    rw [hempty]
    simp

/-- **Zero error collapses to deterministic correctness.** A `RandPass` with `error c = 0` is
right on *every* seed its distribution can produce — the `Pass` notion, recovered from the
randomized one rather than sitting beside it. -/
theorem correct_of_error_eq_zero (p : RandPass) (c : Circuit) (hc : c.Wf)
    (h : p.error c = 0) {s : p.Seed c} (hs : s ∈ (p.dist c).support) :
    Equivalent c.numQubits c.numCbits (p.run c s).gates c.gates := by
  have hzero : (p.dist c).toOuterMeasure (p.failure c) = 0 :=
    le_antisymm (le_of_le_of_eq (p.correct c hc) h) (by simp)
  have hdisj := (PMF.toOuterMeasure_apply_eq_zero_iff _ _).mp hzero
  by_contra hne
  exact (Set.disjoint_left.mp hdisj hs) hne

/-! ## Composition -/

/-- Run `p`, then run `q` on its output *if* `cond` says to, drawing `q`'s seed after seeing
that output.

The condition sees both circuits, which is what the optimizer's round loop needs: it repeats
the pipeline exactly while the gate count keeps falling. `comp` is the unconditional case,
and the two share this proof rather than having one each.

The seed is drawn for `q` either way. That costs nothing — a `RandPass` is a specification,
never run — and it keeps the seed space a plain sigma type, so the measure argument below
needs no transport. -/
def compWhen (p q : RandPass) (cond : Circuit → Circuit → Bool) : RandPass where
  name := q.name ++ " ∘? " ++ p.name
  Seed := fun c => Σ s : p.Seed c, q.Seed (p.run c s)
  dist := fun c => (p.dist c).bind fun s =>
    (q.dist (p.run c s)).map (fun s₂ => (⟨s, s₂⟩ : Σ s : p.Seed c, q.Seed (p.run c s)))
  run := fun c s =>
    if cond c (p.run c s.1) then q.run (p.run c s.1) s.2 else p.run c s.1
  error := fun c => p.error c + ⨆ s : p.Seed c, q.error (p.run c s)
  numQubits_run c s := by
    split
    · rw [q.numQubits_run, p.numQubits_run]
    · rw [p.numQubits_run]
  numCbits_run c s := by
    split
    · rw [q.numCbits_run, p.numCbits_run]
    · rw [p.numCbits_run]
  wf_run c s hc := by
    split
    · exact q.wf_run _ _ (p.wf_run c s.1 hc)
    · exact p.wf_run c s.1 hc
  wellFormed_run c s hwf hc := by
    split
    · exact q.wellFormed_run _ _ (p.wf_run c s.1 hwf) (p.wellFormed_run c s.1 hwf hc)
    · exact p.wellFormed_run c s.1 hwf hc
  flagsOk_run c s hc := by
    split
    · exact q.flagsOk_run _ _ (p.flagsOk_run c s.1 hc)
    · exact p.flagsOk_run c s.1 hc
  correct c hc := by
    set E : ℝ≥0∞ := ⨆ s : p.Seed c, q.error (p.run c s) with hE
    set F : Set (Σ s : p.Seed c, q.Seed (p.run c s)) :=
      {s | ¬ Equivalent c.numQubits c.numCbits
        (if cond c (p.run c s.1) then q.run (p.run c s.1) s.2
          else p.run c s.1).gates c.gates} with hF
    -- the inner measure, for a fixed first-stage seed
    have hinner : ∀ s : p.Seed c,
        (q.dist (p.run c s)).toOuterMeasure
            ((fun s₂ => (⟨s, s₂⟩ : Σ s : p.Seed c, q.Seed (p.run c s))) ⁻¹' F) ≤
          Set.indicator (p.failure c) (fun _ => (1 : ℝ≥0∞)) s + E := by
      intro s
      by_cases hgood : Equivalent c.numQubits c.numCbits (p.run c s).gates c.gates
      · -- `p` succeeded here, so any final failure is a failure of `q`
        have hsub : (fun s₂ => (⟨s, s₂⟩ : Σ s : p.Seed c, q.Seed (p.run c s))) ⁻¹' F ⊆
            q.failure (p.run c s) := by
          intro s₂ hs₂
          simp only [hF, Set.mem_preimage, Set.mem_ofPred_eq] at hs₂
          intro hq
          refine hs₂ ?_
          -- when the condition fails the composite *is* `p`'s output, which is already good
          split
          · have hq' : Equivalent c.numQubits c.numCbits
                (q.run (p.run c s) s₂).gates (p.run c s).gates := by
              rw [p.numQubits_run, p.numCbits_run] at hq
              exact hq
            exact Equivalent.trans hq' hgood
          · exact hgood
        calc (q.dist (p.run c s)).toOuterMeasure
              ((fun s₂ => (⟨s, s₂⟩ : Σ s : p.Seed c, q.Seed (p.run c s))) ⁻¹' F)
            ≤ (q.dist (p.run c s)).toOuterMeasure (q.failure (p.run c s)) :=
              (q.dist (p.run c s)).toOuterMeasure_mono (by
                intro x hx; exact hsub hx.1)
          _ ≤ q.error (p.run c s) := q.correct _ (p.wf_run c s hc)
          _ ≤ E := le_iSup (fun s => q.error (p.run c s)) s
          _ ≤ Set.indicator (p.failure c) (fun _ => (1 : ℝ≥0∞)) s + E := le_add_self
      · -- `p` already failed here; bound the inner probability by one
        have hone : Set.indicator (p.failure c) (fun _ => (1 : ℝ≥0∞)) s = 1 := by
          rw [Set.indicator_of_mem]
          exact hgood
        rw [hone]
        refine le_add_right ?_
        calc (q.dist (p.run c s)).toOuterMeasure _
            ≤ (q.dist (p.run c s)).toOuterMeasure Set.univ :=
              (q.dist (p.run c s)).toOuterMeasure_mono (by intro x _; exact Set.mem_univ x)
          _ = 1 := by rw [PMF.toOuterMeasure_apply]; simpa using (q.dist (p.run c s)).tsum_coe
    calc ((p.dist c).bind fun s =>
            (q.dist (p.run c s)).map
              (fun s₂ => (⟨s, s₂⟩ : Σ s : p.Seed c, q.Seed (p.run c s)))).toOuterMeasure F
        = ∑' s, (p.dist c) s *
            ((q.dist (p.run c s)).map
              (fun s₂ => (⟨s, s₂⟩ : Σ s : p.Seed c, q.Seed (p.run c s)))).toOuterMeasure F := by
          rw [PMF.toOuterMeasure_bind_apply]
      _ = ∑' s, (p.dist c) s * (q.dist (p.run c s)).toOuterMeasure
            ((fun s₂ => (⟨s, s₂⟩ : Σ s : p.Seed c, q.Seed (p.run c s))) ⁻¹' F) := by
          refine tsum_congr fun s => ?_
          rw [PMF.toOuterMeasure_map_apply]
      _ ≤ ∑' s, (p.dist c) s *
            (Set.indicator (p.failure c) (fun _ => (1 : ℝ≥0∞)) s + E) := by
          refine ENNReal.tsum_le_tsum fun s => ?_
          exact mul_le_mul_right (hinner s) _
      _ = (∑' s, (p.dist c) s * Set.indicator (p.failure c) (fun _ => (1 : ℝ≥0∞)) s) +
            (∑' s, (p.dist c) s * E) := by
          rw [← ENNReal.tsum_add]
          exact tsum_congr fun s => by ring
      _ = (p.dist c).toOuterMeasure (p.failure c) + E := by
          have h1 : (∑' s, (p.dist c) s * Set.indicator (p.failure c) (fun _ => (1 : ℝ≥0∞)) s)
              = (p.dist c).toOuterMeasure (p.failure c) := by
            rw [PMF.toOuterMeasure_apply]
            refine tsum_congr fun s => ?_
            by_cases hs : s ∈ p.failure c <;> simp [hs]
          have h2 : (∑' _s : p.Seed c, (p.dist c) _s * E) = E := by
            rw [ENNReal.tsum_mul_right, PMF.tsum_coe, one_mul]
          rw [h1, h2]
      _ ≤ p.error c + E := add_le_add (p.correct c hc) le_rfl

/-- Run `p`, then `q` on its output, drawing `q`'s seed after seeing that output. -/
def comp (p q : RandPass) : RandPass := p.compWhen q (fun _ _ => true)

@[simp] theorem comp_run (p q : RandPass) (c : Circuit) (s : (p.comp q).Seed c) :
    (p.comp q).run c s = q.run (p.run c s.1) s.2 := rfl

@[simp] theorem comp_error (p q : RandPass) (c : Circuit) :
    (p.comp q).error c = p.error c + ⨆ s : p.Seed c, q.error (p.run c s) := rfl

@[simp] theorem compWhen_error (p q : RandPass) (cond : Circuit → Circuit → Bool) (c : Circuit) :
    (p.compWhen q cond).error c = p.error c + ⨆ s : p.Seed c, q.error (p.run c s) := rfl

/-- What one step of a conditional composition computes — the round loop's rule, as a
rewrite. -/
theorem compWhen_run (p q : RandPass) (cond : Circuit → Circuit → Bool) (c : Circuit)
    (s : (p.compWhen q cond).Seed c) :
    (p.compWhen q cond).run c s =
      (if cond c (p.run c s.1) then q.run (p.run c s.1) s.2 else p.run c s.1) := rfl

/-- A pipeline, run left to right: the head runs first, on the original circuit. -/
def pipeline : List RandPass → RandPass
  | [] => RandPass.id
  | p :: ps => p.comp (pipeline ps)

@[simp] theorem pipeline_nil : pipeline [] = RandPass.id := rfl

@[simp] theorem pipeline_cons (p : RandPass) (ps : List RandPass) :
    pipeline (p :: ps) = p.comp (pipeline ps) := rfl

/-- The gate count fell: the optimizer's rule for going round again. -/
def Shrank (c c' : Circuit) : Bool := decide (c'.gates.length < c.gates.length)

/-- **Repeat `p` while it keeps shrinking the circuit**, at most `fuel` times — the driver's
round loop, as a pass rather than as an `IO` loop.

`fuel` bounds the recursion, and `c.gates.length + 1` is always enough: a round that does not
remove a gate ends the loop, so at most `c.gates.length` rounds can continue. -/
def fixpointShrink (p : RandPass) : Nat → RandPass
  | 0 => RandPass.id
  | n + 1 => p.compWhen (p.fixpointShrink n) Shrank

@[simp] theorem fixpointShrink_zero (p : RandPass) : p.fixpointShrink 0 = RandPass.id := rfl

@[simp] theorem fixpointShrink_succ (p : RandPass) (n : Nat) :
    p.fixpointShrink (n + 1) = p.compWhen (p.fixpointShrink n) Shrank := rfl

/-- One turn of the loop, as a rewrite: run the round, and go again exactly when it shrank. -/
theorem fixpointShrink_run (p : RandPass) (n : Nat) (c : Circuit)
    (s : (p.fixpointShrink (n + 1)).Seed c) :
    (p.fixpointShrink (n + 1)).run c s =
      (if (p.run c s.1).gates.length < c.gates.length then
        (p.fixpointShrink n).run (p.run c s.1) s.2 else p.run c s.1) := by
  have hiff : (Shrank c (p.run c s.1) = true) ↔
      ((p.run c s.1).gates.length < c.gates.length) := by simp [Shrank]
  exact if_congr hiff rfl rfl

/-- **The union bound over rounds.** `fuel` rounds are wrong with probability at most `fuel`
times one round's, whatever the rounds do to each other's inputs — the composition draws each
round's seed after seeing the previous round's output, so nothing here assumes independence
of the *events*, only that the draws are fresh. -/
theorem fixpointShrink_error_le (p : RandPass) (B : ℝ≥0∞) (hB : ∀ c, p.error c ≤ B) :
    ∀ (n : Nat) (c : Circuit), (p.fixpointShrink n).error c ≤ n * B := by
  intro n
  induction n with
  | zero => intro c; simp [RandPass.id]
  | succ n ih =>
      intro c
      have hsup : (⨆ s : p.Seed c, (p.fixpointShrink n).error (p.run c s)) ≤ (n : ℝ≥0∞) * B :=
        iSup_le fun s => ih _
      calc (p.fixpointShrink (n + 1)).error c
          = p.error c + ⨆ s : p.Seed c, (p.fixpointShrink n).error (p.run c s) := rfl
        _ ≤ B + (n : ℝ≥0∞) * B := add_le_add (hB c) hsup
        _ = ((n : ℝ≥0∞) + 1) * B := by ring
        _ = ((n + 1 : Nat) : ℝ≥0∞) * B := by push_cast; ring

/-- The same union bound along a pipeline. -/
theorem pipeline_error_le (B : ℝ≥0∞) :
    ∀ (ps : List RandPass), (∀ p ∈ ps, ∀ c, p.error c ≤ B) →
      ∀ c, (pipeline ps).error c ≤ ps.length * B := by
  intro ps
  induction ps with
  | nil => intro _ c; simp [pipeline, RandPass.id]
  | cons p ps ih =>
      intro hB c
      have hsup : (⨆ s : p.Seed c, (pipeline ps).error (p.run c s)) ≤ (ps.length : ℝ≥0∞) * B :=
        iSup_le fun s => ih (fun q hq => hB q (by simp [hq])) _
      calc (pipeline (p :: ps)).error c
          = p.error c + ⨆ s : p.Seed c, (pipeline ps).error (p.run c s) := rfl
        _ ≤ B + (ps.length : ℝ≥0∞) * B := add_le_add (hB p (by simp) c) hsup
        _ = ((ps.length : ℝ≥0∞) + 1) * B := by ring
        _ = ((p :: ps).length : ℝ≥0∞) * B := by
              rw [List.length_cons]; push_cast; ring

/-- **A pipeline of deterministic passes has error zero** — the `Pass` world, recovered. -/
theorem pipeline_error_eq_zero (ps : List Pass) (c : Circuit) :
    (pipeline (ps.map Pass.toRand)).error c = 0 := by
  induction ps generalizing c with
  | nil => rfl
  | cons p ps ih =>
      simp only [List.map_cons, pipeline_cons, comp_error, Pass.toRand]
      simp [ih]

end RandPass

end
end TzapLean
