import TzapLean.NonlinearSoundness

/-!
# Nonlinear fingerprint refinement

This isolates the algebra needed to relate executable packed words to formal polynomial
lifts. `GF128Bridge.lean` constructs the concrete `PackedEvaluation` and uses this
refinement in the executable probability theorem.
-/

namespace TzapLean

noncomputable section

/-- The ring-homomorphism laws required of packed polynomial evaluation at a draw stream. -/
structure PackedEvaluation (draws : Nat → Tag) where
  eval : BoolPolynomial → Tag
  zero : eval 0 = 0
  one : eval 1 = 1
  var_eq : ∀ i, eval (BoolPolynomial.var i) = GF128.normalize (draws i)
  add : ∀ p q, eval (p + q) = GF128.add (eval p) (eval q)
  mul : ∀ p q, eval (p * q) = GF128.mul (eval p) (eval q)

theorem PackedEvaluation.flip {draws : Nat → Tag} (evaluation : PackedEvaluation draws)
    (polynomial : BoolPolynomial) :
    evaluation.eval polynomial.flip = GF128.add (evaluation.eval polynomial) 1 := by
  rw [BoolPolynomial.flip, evaluation.add, evaluation.one]

/-- The runtime and symbolic states carry the same degree and evaluated value on every wire. -/
def NonlinearSim (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (symbolic : NonlinearAState) (runtime : NState) : Prop :=
  runtime.fingerprints.length = symbolic.wires.length ∧
  runtime.fresh = symbolic.fresh ∧
  ∀ q, (runtime.fpOf q).value = evaluation.eval (symbolic.wireOf q).polynomial ∧
    (runtime.fpOf q).degree = (symbolic.wireOf q).degree

theorem nonlinearSim_initial (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (n : Nat) :
    NonlinearSim draws evaluation (NonlinearAState.initial n) (NState.initial draws n) := by
  refine ⟨by simp [NState.initial, NonlinearAState.initial], rfl, fun q => ?_⟩
  simp only [NState.fpOf, NState.initial, NonlinearAState.wireOf,
    NonlinearAState.initial, getD_map_range]
  by_cases hq : q < n
  · simp [hq, Fingerprint.fresh, TrackedPolynomial.fresh, evaluation.var_eq]
  · simp [hq, Fingerprint.zero, TrackedPolynomial.zero, evaluation.zero]

/-- Updating one wire preserves refinement when the new packed value and degree agree. -/
theorem nonlinearSim_set {draws : Nat → Tag} {evaluation : PackedEvaluation draws}
    {symbolic : NonlinearAState} {runtime : NState}
    (hsim : NonlinearSim draws evaluation symbolic runtime)
    (q : Qubit) (polynomial : TrackedPolynomial) (fingerprint : Fingerprint)
    (hvalue : fingerprint.value = evaluation.eval polynomial.polynomial)
    (hdegree : fingerprint.degree = polynomial.degree)
    {symbolicFresh runtimeFresh : Nat} (hfresh : runtimeFresh = symbolicFresh) :
    NonlinearSim draws evaluation
      ⟨symbolic.wires.set q polynomial, symbolicFresh⟩
      ⟨runtime.fingerprints.set q fingerprint, runtimeFresh⟩ := by
  obtain ⟨hlen, -, hwire⟩ := hsim
  refine ⟨by simp [hlen], hfresh, fun r => ?_⟩
  simp only [NState.fpOf, NonlinearAState.wireOf]
  by_cases hr : r = q
  · subst r
    by_cases hq : q < runtime.fingerprints.length
    · rw [List.getD_eq_getElem?_getD, List.getElem?_set_self hq,
        List.getD_eq_getElem?_getD, List.getElem?_set_self (hlen ▸ hq)]
      exact ⟨hvalue, hdegree⟩
    · rw [List.getD_eq_getElem?_getD,
        List.getElem?_eq_none (by simpa using Nat.le_of_not_lt hq),
        List.getD_eq_getElem?_getD,
        List.getElem?_eq_none (by simpa [hlen] using Nat.le_of_not_lt hq)]
      simpa [NState.fpOf, NonlinearAState.wireOf, List.getD_eq_getElem?_getD,
        List.getElem?_eq_none (by simpa using Nat.le_of_not_lt hq),
        List.getElem?_eq_none (by simpa [hlen] using Nat.le_of_not_lt hq)] using hwire q
  · simp only [List.getD_eq_getElem?_getD, List.getElem?_set,
      if_neg (by simpa using Ne.symm hr)]
    exact hwire r

/-- The nonlinear executable transfer follows formal polynomial evaluation, including
the CCX product branch and its fresh-variable cutoff branch. -/
theorem nonlinearSim_step {draws : Nat → Tag} {evaluation : PackedEvaluation draws}
    {symbolic : NonlinearAState} {runtime : NState}
    (hsim : NonlinearSim draws evaluation symbolic runtime) (gate : Gate) :
    NonlinearSim draws evaluation (symbolic.step gate) (runtime.step draws gate) := by
  have hfresh : runtime.fresh = symbolic.fresh := hsim.2.1
  have hwire := hsim.2.2
  cases gate with
  | x q =>
      apply nonlinearSim_set hsim q _ _
      · simp [Fingerprint.add, Fingerprint.one, TrackedPolynomial.add,
          TrackedPolynomial.one, evaluation.add,
          evaluation.one, (hwire q).1]
      · simp [Fingerprint.add, Fingerprint.one, TrackedPolynomial.add,
          TrackedPolynomial.one, (hwire q).2]
      · exact hfresh
  | cnot c t =>
      apply nonlinearSim_set hsim t _ _
      · simp [Fingerprint.add, TrackedPolynomial.add, evaluation.add,
          (hwire t).1, (hwire c).1]
      · simp [Fingerprint.add, TrackedPolynomial.add, (hwire t).2, (hwire c).2]
      · exact hfresh
  | h q =>
      apply nonlinearSim_set hsim q _ _
      · simp [Fingerprint.fresh, TrackedPolynomial.fresh, evaluation.var_eq, hfresh]
      · rfl
      · omega
  | ccx c₁ c₂ t =>
      have hcutoff := NonlinearAState.ccx_cutoff_agrees symbolic runtime c₁ c₂
        (hwire c₁).2.symm (hwire c₂).2.symm
      simp only [NonlinearAState.step, NState.step]
      cases hs : (symbolic.wireOf c₁).mul? (symbolic.wireOf c₂) with
      | some product =>
          cases hr : (runtime.fpOf c₁).mul? (runtime.fpOf c₂) with
          | none => simp [hs, hr] at hcutoff
          | some fingerprint =>
              have hprod : fingerprint.value = evaluation.eval product.polynomial ∧
                  fingerprint.degree = product.degree := by
                simp only [TrackedPolynomial.mul?, Fingerprint.mul?] at hs hr
                split at hs
                · split at hr
                  · cases hs; cases hr
                    exact ⟨by simpa [evaluation.mul] using
                      congrArg₂ GF128.mul (hwire c₁).1 (hwire c₂).1, by
                        simp [Fingerprint.productDegree, (hwire c₁).2, (hwire c₂).2]⟩
                  · simp at hr
                · simp at hs
              apply nonlinearSim_set hsim t _ _
              · simp [Fingerprint.add, TrackedPolynomial.add, evaluation.add,
                  (hwire t).1, hprod.1]
              · simp [Fingerprint.add, TrackedPolynomial.add, (hwire t).2, hprod.2]
              · exact hfresh
      | none =>
          cases hr : (runtime.fpOf c₁).mul? (runtime.fpOf c₂) with
          | some fingerprint => simp [hs, hr] at hcutoff
          | none =>
              apply nonlinearSim_set hsim t _ _
              · simp [Fingerprint.fresh, TrackedPolynomial.fresh, evaluation.var_eq, hfresh]
              · rfl
              · omega
  | reset q =>
      apply nonlinearSim_set hsim q _ _
      · simpa [Fingerprint.zero, TrackedPolynomial.zero] using evaluation.zero.symm
      · rfl
      · exact hfresh
  | _ => exact hsim

theorem nonlinearSim_steps {draws : Nat → Tag} {evaluation : PackedEvaluation draws}
    {symbolic : NonlinearAState} {runtime : NState}
    (hsim : NonlinearSim draws evaluation symbolic runtime) (gates : List Gate) :
    NonlinearSim draws evaluation (symbolic.steps gates) (runtime.steps draws gates) := by
  induction gates generalizing symbolic runtime with
  | nil => exact hsim
  | cons gate gates ih => exact ih (nonlinearSim_step hsim gate)

theorem nonlinearSim_tagOf {draws : Nat → Tag} {evaluation : PackedEvaluation draws}
    {symbolic : NonlinearAState} {runtime : NState}
    (hsim : NonlinearSim draws evaluation symbolic runtime) (q : Qubit) :
    runtime.tagOf q = evaluation.eval (symbolic.wireOf q).polynomial :=
  (hsim.2.2 q).1

/-- When the packed comparison reports a match, its two keys are equal, or differ by
the evaluation of the Boolean complement. -/
theorem matchFingerprint_sound {draws : Nat → Tag} (evaluation : PackedEvaluation draws)
    {pending later : BoolPolynomial} {sign : Bool}
    (hmatch : matchFingerprint (evaluation.eval pending) (evaluation.eval later) =
      some sign) :
    evaluation.eval later = evaluation.eval
      (if sign then pending.flip else pending) := by
  simp only [matchFingerprint] at hmatch
  split at hmatch
  · rename_i heq
    have hs : sign = false := by simpa using hmatch.symm
    subst sign
    simpa using of_decide_eq_true heq
  · split at hmatch
    · rename_i heq
      have hs : sign = true := by simpa using hmatch.symm
      subst sign
      simpa [evaluation.flip] using of_decide_eq_true heq
    · simp at hmatch

/-- Injectivity is required only on the finite set of wire polynomials and complements
that this circuit actually compares. -/
def PolynomialFaithful {draws : Nat → Tag} (evaluation : PackedEvaluation draws)
    (polynomials : List BoolPolynomial) : Prop :=
  ∀ p ∈ polynomials, ∀ q ∈ polynomials,
    evaluation.eval p = evaluation.eval q → p = q

theorem matchFingerprint_exact {draws : Nat → Tag}
    {evaluation : PackedEvaluation draws} {polynomials : List BoolPolynomial}
    (hfaithful : PolynomialFaithful evaluation polynomials)
    {pending later : BoolPolynomial} {sign : Bool}
    (hpending : (if sign then pending.flip else pending) ∈ polynomials)
    (hlater : later ∈ polynomials)
    (hmatch : matchFingerprint (evaluation.eval pending) (evaluation.eval later) =
      some sign) :
    later = if sign then pending.flip else pending :=
  hfaithful later hlater _ hpending (matchFingerprint_sound evaluation hmatch)

end

end TzapLean
