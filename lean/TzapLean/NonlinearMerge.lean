import TzapLean.NonlinearRefinement
import TzapLean.NonlinearCompared
import TzapLean.Merge

/-!
# Nonlinear phase-hop soundness

The local rewrite at the heart of phase folding works for formal Boolean polynomials,
including CCX products and fresh variables introduced by the degree cutoff. This theorem
does not assume that the polynomials are affine.
-/

namespace TzapLean

noncomputable section

/-- A pair of rotations can be merged across a unitary segment whenever the nonlinear
analyzer proves that the later wire is the earlier wire or its complement. -/
theorem nonlinear_merge_equivalent {n m : Nat} (segment : List Gate)
    (hunitary : ∀ gate ∈ segment, gate.isUnitary = true)
    (st : NonlinearAState) (hbounded : st.VariablesBounded)
    (hlen : st.wires.length = n) (hgeneric : st.Generic n)
    (θ φ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((st.steps segment).wireOf q').polynomial =
      if sign then (st.wireOf q).polynomial.flip else (st.wireOf q).polynomial) :
    Equivalent n m
      (Gate.rz θ q :: segment ++ [Gate.rz φ q'])
      (segment ++ [Gate.rz (φ + signedAngle sign θ) q']) := by
  have hpath : ∀ b k b' : Basis n,
      unitary n [] k b ≠ 0 → unitary n segment b' k ≠ 0 →
        k.get q = (b'.get q' != sign) := by
    intro _ k b' _ hseg
    exact NonlinearAState.path_of_polynomial_eq hunitary hbounded hlen hgeneric hpoly
      k b' hseg
  simpa using merge_equivalent (n := n) (m := m) [] segment θ φ q q' sign
    (by simp) hunitary hpath

/-- The executable packed comparison licenses the local nonlinear rewrite whenever the
three formal polynomials it may compare are collision-free. -/
theorem nonlinear_merge_equivalent_of_match {n m : Nat} (segment : List Gate)
    (hunitary : ∀ gate ∈ segment, gate.isUnitary = true)
    (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (st : NonlinearAState) (runtime : NState)
    (hsim : NonlinearSim draws evaluation st runtime)
    (hbounded : st.VariablesBounded) (hlen : st.wires.length = n)
    (hgeneric : st.Generic n)
    (θ φ : ℚ) (q q' : Qubit) (sign : Bool)
    (hfaithful : PolynomialFaithful evaluation
      [(st.wireOf q).polynomial, (st.wireOf q).polynomial.flip,
        ((st.steps segment).wireOf q').polynomial])
    (hmatch : matchFingerprint (runtime.tagOf q)
      ((runtime.steps draws segment).tagOf q') = some sign) :
    Equivalent n m
      (Gate.rz θ q :: segment ++ [Gate.rz φ q'])
      (segment ++ [Gate.rz (φ + signedAngle sign θ) q']) := by
  have hsimMiddle := nonlinearSim_steps hsim segment
  have hmatch' : matchFingerprint
      (evaluation.eval (st.wireOf q).polynomial)
      (evaluation.eval ((st.steps segment).wireOf q').polynomial) = some sign := by
    rw [← nonlinearSim_tagOf hsim q, ← nonlinearSim_tagOf hsimMiddle q']
    exact hmatch
  have hpoly : ((st.steps segment).wireOf q').polynomial =
      if sign then (st.wireOf q).polynomial.flip else (st.wireOf q).polynomial :=
    matchFingerprint_exact hfaithful (by cases sign <;> simp)
      (by simp) hmatch'
  exact nonlinear_merge_equivalent segment hunitary st hbounded hlen hgeneric
    θ φ q q' sign hpoly

/-- The usual circuit-wide comparison set supplies the local faithfulness premise. -/
theorem nonlinear_merge_equivalent_of_relevantFaithful {n m : Nat}
    (segment : List Gate)
    (hunitary : ∀ gate ∈ segment, gate.isUnitary = true)
    (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (st : NonlinearAState) (runtime : NState)
    (hsim : NonlinearSim draws evaluation st runtime)
    (hbounded : st.VariablesBounded) (hlen : st.wires.length = n)
    (hgeneric : st.Generic n)
    (θ φ : ℚ) (q q' : Qubit) (sign : Bool)
    (hfaithful : PolynomialFaithful evaluation
      (nonlinearRelevant n st (segment ++ [Gate.rz φ q'])))
    (hmatch : matchFingerprint (runtime.tagOf q)
      ((runtime.steps draws segment).tagOf q') = some sign) :
    Equivalent n m
      (Gate.rz θ q :: segment ++ [Gate.rz φ q'])
      (segment ++ [Gate.rz (φ + signedAngle sign θ) q']) := by
  have hpending : (st.wireOf q).polynomial ∈
      nonlinearVisited n st (segment ++ [Gate.rz φ q']) :=
    nonlinearWire_mem_visited hlen q _
  have hlater : ((st.steps segment).wireOf q').polynomial ∈
      nonlinearVisited n st (segment ++ [Gate.rz φ q']) :=
    nonlinearVisited_append_mem segment [Gate.rz φ q']
      (nonlinearWire_mem_visited (by rw [NonlinearAState.length_steps, hlen]) q' _)
  have hsimMiddle := nonlinearSim_steps hsim segment
  have hmatch' : matchFingerprint
      (evaluation.eval (st.wireOf q).polynomial)
      (evaluation.eval ((st.steps segment).wireOf q').polynomial) = some sign := by
    rw [← nonlinearSim_tagOf hsim q, ← nonlinearSim_tagOf hsimMiddle q']
    exact hmatch
  have hcomparison : (if sign then (st.wireOf q).polynomial.flip
      else (st.wireOf q).polynomial) ∈
      nonlinearRelevant n st (segment ++ [Gate.rz φ q']) := by
    cases sign with
    | false => exact List.mem_append.2 (Or.inl hpending)
    | true => exact List.mem_append.2 (Or.inr
        (List.mem_map.2 ⟨_, hpending, rfl⟩))
  have hpoly := matchFingerprint_exact hfaithful hcomparison
    (List.mem_append.2 (Or.inl hlater)) hmatch'
  exact nonlinear_merge_equivalent segment hunitary st hbounded hlen hgeneric
    θ φ q q' sign hpoly

/-- Every successful executable search identifies one later rotation and the intervening
segment, with the exact packed comparison used for its sign. -/
theorem mergeIntoNonlinear_spec (draws : Nat → Tag) (tag : Tag) (θ : ℚ) :
    ∀ (gates result : List Gate) (runtime : NState),
      mergeIntoNonlinear draws runtime tag θ gates = some result →
      ∃ (middle rest : List Gate) (later : Gate) (φ : ℚ) (q' : Qubit) (sign : Bool),
        gates = middle ++ later :: rest ∧
        result = middle ++ Gate.rz (φ + signedAngle sign θ) q' :: rest ∧
        (∀ gate ∈ middle, gate.isUnitary = true ∨ gate.isMeasurement = true) ∧
        rotAngle later = some (φ, q') ∧
        matchFingerprint tag ((runtime.steps draws middle).tagOf q') = some sign := by
  intro gates
  induction gates with
  | nil =>
      intro result runtime h
      simp [mergeIntoNonlinear] at h
  | cons gate gates ih =>
      intro result runtime h
      simp only [mergeIntoNonlinear] at h
      split at h
      · rename_i hallowed
        split at h
        · rename_i φ q' hrot
          split at h
          · rename_i sign hmatch
            simp only [Option.some.injEq] at h
            subst result
            exact ⟨[], gates, gate, φ, q', sign, rfl, rfl, by simp, hrot, hmatch⟩
          · rcases Option.map_eq_some_iff.1 h with ⟨tail, htail, rfl⟩
            obtain ⟨middle, rest, later, φ', q'', sign, h1, h2, h3, h4, h5⟩ :=
              ih tail _ htail
            exact ⟨gate :: middle, rest, later, φ', q'', sign,
              by rw [h1]; rfl, by rw [h2]; rfl,
              by intro x hx
                 rcases List.mem_cons.1 hx with rfl | hx
                 · simp at hallowed ⊢
                 · exact h3 x hx,
              h4, h5⟩
        · rcases Option.map_eq_some_iff.1 h with ⟨tail, htail, rfl⟩
          obtain ⟨middle, rest, later, φ', q'', sign, h1, h2, h3, h4, h5⟩ :=
            ih tail _ htail
          exact ⟨gate :: middle, rest, later, φ', q'', sign,
            by rw [h1]; rfl, by rw [h2]; rfl,
            by intro x hx
               rcases List.mem_cons.1 hx with rfl | hx
               · simp at hallowed ⊢
               · exact h3 x hx,
            h4, h5⟩
      · exact absurd h (by simp)

end

end TzapLean
