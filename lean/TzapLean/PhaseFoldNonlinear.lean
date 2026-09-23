import TzapLean.PhaseFoldProof
import TzapLean.GF128

/-!
# Nonlinear executable state for `PhaseFoldRand`

This module is the exact executable transfer system used by the Rust pass.  A wire carries a
packed `GF(2^128)` evaluation and a total-degree upper bound.  CCX evaluates
`target + control₁ * control₂`; if the product degree would exceed `2^32`, the target is
replaced by a fresh opaque degree-one value.  Reset writes the field zero.

The affine state in `PhaseFold.lean` remains available to its existing proof. The packed-field
refinement and nonlinear collision theorem are connected in `GF128Bridge.lean`.
-/

namespace TzapLean

/-- Per-wire nonlinear fingerprints and the next unused random draw. -/
structure NState where
  fingerprints : List Fingerprint
  fresh : Nat

namespace NState

/-- Wire `q`'s fingerprint, or the constant zero outside the represented register. -/
def fpOf (st : NState) (q : Qubit) : Fingerprint :=
  st.fingerprints.getD q Fingerprint.zero

/-- The packed field value used as the phase-group key. -/
def tagOf (st : NState) (q : Qubit) : Tag := (st.fpOf q).value

/-- Each input wire begins as an independent degree-one variable. -/
def initial (draws : Nat → Tag) (n : Nat) : NState where
  fingerprints := (List.range n).map (fun i => Fingerprint.fresh (draws i))
  fresh := n

/-- Rust-compatible nonlinear transfer functions. -/
def step (draws : Nat → Tag) (st : NState) (g : Gate) : NState :=
  match g with
  | .x q =>
      { st with fingerprints := st.fingerprints.set q ((st.fpOf q).add Fingerprint.one) }
  | .cnot c t =>
      { st with fingerprints := st.fingerprints.set t ((st.fpOf t).add (st.fpOf c)) }
  | .h q =>
      { fingerprints := st.fingerprints.set q (Fingerprint.fresh (draws st.fresh))
        fresh := st.fresh + 1 }
  | .ccx c₁ c₂ t =>
      match (st.fpOf c₁).mul? (st.fpOf c₂) with
      | some product =>
          { st with fingerprints := st.fingerprints.set t ((st.fpOf t).add product) }
      | none =>
          { fingerprints := st.fingerprints.set t (Fingerprint.fresh (draws st.fresh))
            fresh := st.fresh + 1 }
  | .reset q =>
      { st with fingerprints := st.fingerprints.set q Fingerprint.zero }
  | _ => st

def steps (draws : Nat → Tag) (st : NState) : List Gate → NState
  | [] => st
  | g :: gs => steps draws (st.step draws g) gs

end NState

/-- Equality/complement matching in the extension field.  The Boolean constant one is the
field element `1`, not the all-ones packed word used by the old affine bit-vector model. -/
def matchFingerprint (pending later : Tag) : Option Bool :=
  if later == pending then some false
  else if later == GF128.add pending 1 then some true
  else none

def mergeIntoNonlinear (draws : Nat → Tag) (st : NState) (tag : Tag) (θ : ℚ) :
    List Gate → Option (List Gate)
  | [] => none
  | g :: gs =>
      if g.isUnitary || g.isMeasurement then
        match rotAngle g with
        | some (φ, q') =>
            match matchFingerprint tag (st.tagOf q') with
            | some sign => some (Gate.rz (φ + signedAngle sign θ) q' :: gs)
            | none => (mergeIntoNonlinear draws (st.step draws g) tag θ gs).map (g :: ·)
        | none => (mergeIntoNonlinear draws (st.step draws g) tag θ gs).map (g :: ·)
      else none

theorem mergeIntoNonlinear_length (draws : Nat → Tag) (tag : Tag) (θ : ℚ) :
    ∀ (gs gs' : List Gate) (st : NState),
      mergeIntoNonlinear draws st tag θ gs = some gs' → gs'.length = gs.length := by
  intro gs
  induction gs with
  | nil => intro gs' st h; simp [mergeIntoNonlinear] at h
  | cons g gs ih =>
      intro gs' st h
      simp only [mergeIntoNonlinear] at h
      split at h
      · split at h
        · split at h
          · simp only [Option.some.injEq] at h
            subst h
            simp
          · rcases Option.map_eq_some_iff.1 h with ⟨t, ht, rfl⟩
            simp [ih t _ ht]
        · rcases Option.map_eq_some_iff.1 h with ⟨t, ht, rfl⟩
          simp [ih t _ ht]
      · exact absurd h (by simp)

def canonicalFingerprint (value : Tag) : Tag :=
  min value (GF128.add value 1)

/-- Precompute which rotations have a later member of the same nonlinear fingerprint group. -/
def nonlinearMergeTargets (draws : Nat → Tag) (n : Nat) (gs : List Gate) : Array Bool :=
  Id.run do
    let arr := gs.toArray
    let mut canons : Array (Option Tag) := Array.replicate arr.size none
    let mut st := NState.initial draws n
    for h : i in [0 : arr.size] do
      let g := arr[i]
      match rotAngle g with
      | some (_, q) => canons := canons.set! i (some (canonicalFingerprint (st.tagOf q)))
      | none => pure ()
      st := st.step draws g
    let mut result : Array Bool := Array.replicate arr.size false
    let mut seen : Std.HashSet Tag := ∅
    for i in [0 : arr.size] do
      let j := arr.size - 1 - i
      let g := arr[j]!
      if !(g.isUnitary || g.isMeasurement) then
        seen := ∅
      else
        match canons[j]! with
        | some key =>
            result := result.set! j (seen.contains key)
            seen := seen.insert key
        | none => pure ()
    return result

def foldFromNonlinear (draws : Nat → Tag) (targets : Array Bool) (st : NState) :
    Nat → List Gate → List Gate
  | _, [] => []
  | at_, g :: gs =>
      match rotAngle g with
      | some (θ, q) =>
          if st.tagOf q == 0 || st.tagOf q == 1 then
            foldFromNonlinear draws targets (st.step draws g) (at_ + 1) gs
          else
            if targets[at_]?.getD true then
              match hm : mergeIntoNonlinear draws st (st.tagOf q) θ gs with
              | some gs' => foldFromNonlinear draws targets st (at_ + 1) gs'
              | none => g :: foldFromNonlinear draws targets (st.step draws g) (at_ + 1) gs
            else g :: foldFromNonlinear draws targets (st.step draws g) (at_ + 1) gs
      | none => g :: foldFromNonlinear draws targets (st.step draws g) (at_ + 1) gs
  termination_by _ gs => gs.length
  decreasing_by
    all_goals
      first
        | (rw [mergeIntoNonlinear_length draws (st.tagOf q) θ gs gs' st hm]; simp)
        | simp

/-- Executable nonlinear phase folding, matching Rust's state transitions. -/
def phaseFoldGatesNonlinear (draws : Nat → Tag) (n : Nat) (gs : List Gate) : List Gate :=
  emitAll (foldFromNonlinear draws (nonlinearMergeTargets draws n gs)
    (NState.initial draws n) 0 gs)

def phaseFoldNonlinear (draws : Nat → Tag) (c : RawCircuit) : RawCircuit :=
  c.withGates (phaseFoldGatesNonlinear draws c.numQubits c.gates)

/-! ## Structural preservation

These lemmas deliberately do not use any algebraic property of the fingerprints.  The fold
only removes rotations or replaces a rotation by another rotation on the same wire, so gate
well-formedness and operand ranges are independent of the collision proof. -/

theorem nonlinear_emitRotation_wf (q : Qubit) (a : ℚ) : ∀ g ∈ emitRotation q a, g.Wf := by
  by_cases h0 : BlockState.angleMod a = 0
  · rw [emitRotation_eq_nil h0]; simp
  · cases hcl : classifyQuarterPi (BlockState.angleMod a) with
    | some j =>
        rw [emitRotation_eq_diagRun h0 hcl]
        match j with
        | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 => intro g hg; fin_cases hg <;> trivial
        | (i + 8) => intro g hg; simp [diagRun] at hg
    | none =>
        rw [emitRotation_eq_rz h0 hcl]
        intro g hg
        rw [List.mem_singleton.1 hg]
        trivial

theorem nonlinear_emitAll_wf {gs : List Gate} (h : ∀ g ∈ gs, g.Wf) :
    ∀ g ∈ emitAll gs, g.Wf := by
  induction gs with
  | nil => intro g hg; simp [emitAll] at hg
  | cons x xs ih =>
      intro g hg
      rw [emitAll] at hg
      rcases List.mem_append.1 hg with hg | hg
      · cases hrot : rotAngle x with
        | some p =>
            obtain ⟨a, q⟩ := p
            rw [hrot] at hg
            exact nonlinear_emitRotation_wf q a g hg
        | none =>
            rw [hrot] at hg
            rw [List.mem_singleton.1 hg]
            exact h x (by simp)
      · exact ih (fun y hy => h y (by simp [hy])) g hg

theorem mergeIntoNonlinear_wf (draws : Nat → Tag) (tag : Tag) (θ : ℚ) :
    ∀ (gs gs' : List Gate) (st : NState),
      (∀ g ∈ gs, g.Wf) → mergeIntoNonlinear draws st tag θ gs = some gs' →
      ∀ g ∈ gs', g.Wf := by
  intro gs
  induction gs with
  | nil => intro gs' st _ hm; simp [mergeIntoNonlinear] at hm
  | cons x xs ih =>
      intro gs' st hwf hm
      have hx : x.Wf := hwf x (by simp)
      simp only [mergeIntoNonlinear] at hm
      split at hm
      · split at hm
        · split at hm
          · simp only [Option.some.injEq] at hm
            subst gs'
            intro g hg
            rcases List.mem_cons.1 hg with rfl | hg
            · trivial
            · exact hwf g (by simp [hg])
          · rcases Option.map_eq_some_iff.1 hm with ⟨ys, hys, rfl⟩
            intro g hg
            rcases List.mem_cons.1 hg with rfl | hg
            · exact hx
            · exact ih ys _ (fun y hy => hwf y (by simp [hy])) hys g hg
        · rcases Option.map_eq_some_iff.1 hm with ⟨ys, hys, rfl⟩
          intro g hg
          rcases List.mem_cons.1 hg with rfl | hg
          · exact hx
          · exact ih ys _ (fun y hy => hwf y (by simp [hy])) hys g hg
      · simp at hm

theorem foldFromNonlinear_wf (draws : Nat → Tag) (targets : Array Bool) :
    ∀ (N : Nat) (gs : List Gate), gs.length ≤ N → ∀ (at_ : Nat) (st : NState),
      (∀ g ∈ gs, g.Wf) → ∀ g ∈ foldFromNonlinear draws targets st at_ gs, g.Wf := by
  intro N
  induction N with
  | zero =>
      intro gs hlen at_ st _ g hg
      rw [List.eq_nil_of_length_eq_zero (Nat.le_zero.1 hlen)] at hg
      simp [foldFromNonlinear] at hg
  | succ N ih =>
      intro gs hlen at_ st hwf
      cases gs with
      | nil => intro g hg; simp [foldFromNonlinear] at hg
      | cons x xs =>
          have hxs : xs.length ≤ N := by simp only [List.length_cons] at hlen; omega
          have keep : ∀ g ∈ x :: foldFromNonlinear draws targets (st.step draws x)
                (at_ + 1) xs, g.Wf := by
            intro g hg
            rcases List.mem_cons.1 hg with rfl | hg
            · exact hwf _ (by simp)
            · exact ih xs hxs _ _ (fun y hy => hwf y (by simp [hy])) g hg
          simp only [foldFromNonlinear]
          split
          · split
            · exact ih xs hxs _ _ (fun y hy => hwf y (by simp [hy]))
            · split
              · split
                · rename_i θ q hrot hconstant htarget gs' hm
                  have hlen' : gs'.length ≤ N := by
                    rw [mergeIntoNonlinear_length draws (st.tagOf q) θ xs gs' st hm]
                    exact hxs
                  exact ih gs' hlen' (at_ + 1) st
                    (mergeIntoNonlinear_wf draws (st.tagOf q) θ xs gs' st
                      (fun g hg => hwf g (by simp [hg])) hm)
                · exact keep
              · exact keep
          · exact keep

theorem phaseFoldGatesNonlinear_wf (draws : Nat → Tag) {n : Nat} {gs : List Gate}
    (h : ∀ g ∈ gs, g.Wf) : ∀ g ∈ phaseFoldGatesNonlinear draws n gs, g.Wf :=
  nonlinear_emitAll_wf (foldFromNonlinear_wf draws _ gs.length gs le_rfl 0 _ h)

@[simp] theorem phaseFoldNonlinear_numQubits (draws : Nat → Tag) (c : RawCircuit) :
    (phaseFoldNonlinear draws c).numQubits = c.numQubits := rfl

@[simp] theorem phaseFoldNonlinear_numCbits (draws : Nat → Tag) (c : RawCircuit) :
    (phaseFoldNonlinear draws c).numCbits = c.numCbits := rfl

end TzapLean
