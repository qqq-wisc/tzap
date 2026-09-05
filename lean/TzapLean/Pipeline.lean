import TzapLean.RandPass
import TzapLean.CnotMinProof
import TzapLean.PhaseFoldRand
import TzapLean.SuperOptProof

/-!
# The Passes, in the Randomized World

`RandPass` is the common currency. Three of tzap's four passes are deterministic, so they
enter at `error = 0` with a one-point seed; phase folding is the one that consumes randomness.
This file is only that embedding, plus what it costs (nothing).

The pipeline itself — which passes, in which order, repeated how often — is *not* here. It is
`passOf` / `tzapRound` / `tzapRun` in `TzapLean/Optimize.lean`, sitting next to `verifiedStep`
and `runToFixpoint`, the driver code that runs it. Keeping the two adjacent is deliberate: a
pipeline modelled in one file and executed in another is a pipeline that drifts.
-/

namespace TzapLean

open scoped ENNReal

noncomputable section

/-- Embed a concrete deterministic transformation in the legacy randomized theory.  This is
kept separate from `Pass`: the executable `Pass` interface is indexed by `Circuit.Checked`. -/
def deterministicRand (name : String) (f : Circuit → Circuit)
    (hn : ∀ c, (f c).numQubits = c.numQubits)
    (hm : ∀ c, (f c).numCbits = c.numCbits)
    (hwf : ∀ c, c.Wf → (f c).Wf)
    (hform : ∀ c, c.Wf → c.WellFormed → (f c).WellFormed)
    (hflags : ∀ c, c.FlagsOk → (f c).FlagsOk)
    (hcorrect : ∀ c, c.Wf → Equivalent c.numQubits c.numCbits (f c).gates c.gates) : RandPass where
  name := name
  Seed := fun _ => Unit
  dist := fun _ => PMF.pure ()
  run := fun c _ => f c
  error := fun _ => 0
  numQubits_run c _ := hn c
  numCbits_run c _ := hm c
  wf_run c _ := hwf c
  wellFormed_run c _ := hform c
  flagsOk_run c _ := hflags c
  correct c hc := by
    have : {s : Unit | ¬ Equivalent c.numQubits c.numCbits (f c).gates c.gates} = ∅ := by
      ext s; simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact hcorrect c hc
    rw [this]; simp

/-- `CancelGates` as a zero-error randomized pass. -/
def CancelGatesR : RandPass := deterministicRand "Gate cancellation" cancelGatesCircuit
  (by intro; rfl) (by intro; rfl) (fun _ => cancelGates_wf) (fun _ _ => cancelGates_inRange)
  (fun c _ => Circuit.flagsOk_withGates _ _) (fun c => cancelGates_correct c.gates)

/-- `CnotMin` as a zero-error randomized pass. -/
def CnotMinR : RandPass := deterministicRand "CNOT minimization" cnotMinCircuit
  (by intro; rfl) (by intro; rfl)
  (fun c => cnotMinGates_wf _ _ c.gates) (fun c _ => cnotMinGates_inRange _ _ c.gates)
  (fun c _ => Circuit.flagsOk_withGates _ _) (fun c => cnotMinGates_correct _ _ c.gates)

/-- `SuperOpt` as a zero-error randomized pass: it verifies each rewrite by exact matrix
comparison, so despite the search inside it there is nothing probabilistic about it. -/
def SuperOptR (cfg : SuperOptConfig) (tbl : SynthTable) : RandPass :=
  deterministicRand "Superoptimization" (superOpt cfg tbl)
    (by intro; rfl) (by intro; rfl)
    (fun c => superOptGates_wf cfg tbl c.gates) (fun c _ => superOptGates_inRange cfg tbl c.gates)
    (fun c _ => Circuit.flagsOk_withGates _ _) (fun c _ => superOptGates_correct cfg tbl c.gates)

@[simp] theorem CancelGatesR_error (c : Circuit) : CancelGatesR.error c = 0 := rfl
@[simp] theorem CnotMinR_error (c : Circuit) : CnotMinR.error c = 0 := rfl

@[simp] theorem SuperOptR_error (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit) :
    (SuperOptR cfg tbl).error c = 0 := rfl

@[simp] theorem CancelGatesR_run (c : Circuit) (s : CancelGatesR.Seed c) :
    CancelGatesR.run c s = cancelGatesCircuit c := rfl

@[simp] theorem CnotMinR_run (c : Circuit) (s : CnotMinR.Seed c) :
    CnotMinR.run c s = cnotMinCircuit c := rfl

@[simp] theorem SuperOptR_run (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit)
    (s : (SuperOptR cfg tbl).Seed c) : (SuperOptR cfg tbl).run c s = superOpt cfg tbl c := rfl

/-- Phase folding's bound, in closed form, for reference from the pipeline: `t` compared
parities collide with probability at most `C(t,2)·2⁻ᵏ`, so doubling the tag width squares the
odds against it. -/
theorem phaseFold_error (k : Nat) (c : Circuit) :
    (PhaseFoldRand k).error c =
      ((relevantForms c).length.choose 2 : ℝ≥0∞) * ((2 : ℝ≥0∞)⁻¹) ^ k :=
  PhaseFoldRand_error k c

end
end TzapLean
