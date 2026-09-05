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
`passOf` / `tzapRound` / `tzapRun` in `TzapLean/Optimize.lean`, sitting next to `stepOf` and
`runToFixpoint`, the driver code that runs it. Keeping the two adjacent is deliberate: a
pipeline modelled in one file and executed in another is a pipeline that drifts.
-/

namespace TzapLean

open scoped ENNReal

noncomputable section

/-- `CancelGates` as a zero-error randomized pass. -/
def CancelGatesR : RandPass := Pass.toRand CancelGates

/-- `CnotMin` as a zero-error randomized pass. -/
def CnotMinR : RandPass := Pass.toRand CnotMin

/-- `SuperOpt` as a zero-error randomized pass: it verifies each rewrite by exact matrix
comparison, so despite the search inside it there is nothing probabilistic about it. -/
def SuperOptR (cfg : SuperOptConfig) (tbl : SynthTable) : RandPass :=
  Pass.toRand (SuperOpt cfg tbl)

@[simp] theorem CancelGatesR_error (c : Circuit) : CancelGatesR.error c = 0 := rfl
@[simp] theorem CnotMinR_error (c : Circuit) : CnotMinR.error c = 0 := rfl

@[simp] theorem SuperOptR_error (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit) :
    (SuperOptR cfg tbl).error c = 0 := rfl

@[simp] theorem CancelGatesR_run (c : Circuit) (s : CancelGatesR.Seed c) :
    CancelGatesR.run c s = CancelGates.run c := rfl

@[simp] theorem CnotMinR_run (c : Circuit) (s : CnotMinR.Seed c) :
    CnotMinR.run c s = CnotMin.run c := rfl

@[simp] theorem SuperOptR_run (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit)
    (s : (SuperOptR cfg tbl).Seed c) : (SuperOptR cfg tbl).run c s = superOpt cfg tbl c := rfl

/-- **Nothing is given up by moving to the randomized setting.** A pipeline built only from
`Pass`es carries error `0`, and `correct_of_error_eq_zero` turns that back into the
unconditional statement `Pass.correct_runAll` already made. -/
theorem pipeline_toRand_correct (ps : List Pass) (c : Circuit) (hc : c.Wf)
    {s : (RandPass.pipeline (ps.map Pass.toRand)).Seed c}
    (hs : s ∈ ((RandPass.pipeline (ps.map Pass.toRand)).dist c).support) :
    Equivalent c.numQubits c.numCbits
      ((RandPass.pipeline (ps.map Pass.toRand)).run c s).gates c.gates :=
  RandPass.correct_of_error_eq_zero _ c hc (RandPass.pipeline_error_eq_zero ps c) hs

/-- Phase folding's bound, in closed form, for reference from the pipeline: `t` compared
parities collide with probability at most `C(t,2)·2⁻ᵏ`, so doubling the tag width squares the
odds against it. -/
theorem phaseFold_error (k : Nat) (c : Circuit) :
    (PhaseFoldRand k).error c =
      ((relevantForms c).length.choose 2 : ℝ≥0∞) * ((2 : ℝ≥0∞)⁻¹) ^ k :=
  PhaseFoldRand_error k c

end
end TzapLean
