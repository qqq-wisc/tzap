import TzapLean.RandPass
import TzapLean.CnotMinProof
import TzapLean.SuperOptProof
import TzapLean.Decompose

/-!
# The Passes, in the Randomized World

`RandPass` is the common currency. The other pipeline passes are deterministic, so they
enter at `error = 0` with a one-point seed; phase folding is the one that consumes randomness.
This file is only that embedding, plus what it costs (nothing).

The pipeline itself — which passes, in which order, repeated how often — is *not* here. It is
`passOf` / `tzapRound` / `tzapRun` in `TzapLean/Optimize.lean`, sitting next to
`executableStep`, the driver code that runs it. Keeping the two adjacent is deliberate: a
pipeline modelled in one file and executed in another is a pipeline that drifts.
-/

namespace TzapLean

open scoped ENNReal

noncomputable section

/-- Embed a deterministic verified pass in the randomized theory with a one-point seed and
zero error. Both interfaces act on the same checked circuit type, so register-size and `Wf`
obligations are not duplicated; only the output-boundary range and cache invariants remain. -/
def deterministicRand (p : Pass)
    (hform : ∀ {n m} (c : Circuit n m), c.raw.WellFormed → (p.run c).raw.WellFormed)
    (hflags : ∀ {n m} (c : Circuit n m), c.raw.FlagsOk → (p.run c).raw.FlagsOk) : RandPass where
  name := p.name
  Seed := fun _ => Unit
  dist := fun _ => PMF.pure ()
  run := fun c _ => p.run c
  error := fun _ => 0
  wellFormed_run c _ := hform c
  flagsOk_run c _ := hflags c
  correct c := by
    have : {s : Unit | ¬ (p.run c).Equivalent c} = ∅ := by
      ext s; simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_not]
      exact p.correct c
    rw [this]; simp

/-- `CancelGates` as a zero-error randomized pass. -/
def CancelGatesR : RandPass := deterministicRand CancelGates
  (fun _ => cancelGates_inRange)
  (fun _ _ => RawCircuit.flagsOk_withGates _ _)

/-- `CnotMin` as a zero-error randomized pass. -/
def CnotMinR : RandPass := deterministicRand CnotMin
  (fun c => cnotMinGates_inRange _ _ c.raw.gates)
  (fun _ _ => RawCircuit.flagsOk_withGates _ _)

/-- `SuperOpt` as a zero-error randomized pass: it verifies each rewrite by exact matrix
comparison, so despite the search inside it there is nothing probabilistic about it. -/
def SuperOptR (cfg : SuperOptConfig) (murm : Murm) : RandPass :=
  deterministicRand (SuperOpt cfg murm)
    (fun c => superOptGates_inRange cfg murm c.raw.gates)
    (fun _ _ => RawCircuit.flagsOk_withGates _ _)

/-- Exact CCX/CCZ lowering as a zero-error randomized pass. -/
def DecomposeToffoliR : RandPass := deterministicRand DecomposeToffoli
  decomposeToffoliChecked_wellFormed decomposeToffoliChecked_flagsOk

/-- Exact CZ lowering as a zero-error randomized pass. -/
def DecomposeCzR : RandPass := deterministicRand DecomposeCz
  decomposeCzChecked_wellFormed decomposeCzChecked_flagsOk

@[simp] theorem CancelGatesR_error (c : Circuit n m) : CancelGatesR.error c = 0 := rfl
@[simp] theorem CnotMinR_error (c : Circuit n m) : CnotMinR.error c = 0 := rfl
@[simp] theorem DecomposeToffoliR_error (c : Circuit n m) : DecomposeToffoliR.error c = 0 := rfl
@[simp] theorem DecomposeCzR_error (c : Circuit n m) : DecomposeCzR.error c = 0 := rfl

@[simp] theorem SuperOptR_error (cfg : SuperOptConfig) (murm : Murm) (c : Circuit n m) :
    (SuperOptR cfg murm).error c = 0 := rfl

@[simp] theorem CancelGatesR_run (c : Circuit n m) (s : CancelGatesR.Seed c) :
    CancelGatesR.run c s = CancelGates.run c := rfl

@[simp] theorem CnotMinR_run (c : Circuit n m) (s : CnotMinR.Seed c) :
    CnotMinR.run c s = CnotMin.run c := rfl

@[simp] theorem SuperOptR_run (cfg : SuperOptConfig) (murm : Murm) (c : Circuit n m)
    (s : (SuperOptR cfg murm).Seed c) :
    (SuperOptR cfg murm).run c s = (SuperOpt cfg murm).run c := rfl

end
end TzapLean
