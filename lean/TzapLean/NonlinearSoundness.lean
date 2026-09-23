import TzapLean.NonlinearAnalysis
import TzapLean.Analysis

/-!
# Semantic soundness of nonlinear CCX tracking

The executable fingerprint is an evaluation of the formal state developed in
`NonlinearAnalysis`.  This module proves the semantic half of that statement: on Boolean
valuations, every below-cutoff CCX update is exactly Toffoli, while the cutoff fallback and
Hadamard introduce a fresh variable that can explain the observed output branch.

The packed-field refinement and its collision probability are separate obligations; no field
property is assumed here.
-/

namespace TzapLean

noncomputable section

namespace NonlinearAState

/-- Every wire polynomial mentions only variables already allocated by the state. -/
def VariablesBounded (st : NonlinearAState) : Prop :=
  ∀ q : Qubit, BoolPolynomial.Bounded st.fresh (st.wireOf q).polynomial

theorem fresh_le_step (st : NonlinearAState) (g : Gate) : st.fresh ≤ (st.step g).fresh := by
  cases g <;> simp [step]
  split <;> simp

theorem variablesBounded_set {st : NonlinearAState} {q : Qubit}
    {value : TrackedPolynomial} {fresh : Nat}
    (hst : ∀ r, BoolPolynomial.Bounded fresh (st.wireOf r).polynomial)
    (hvalue : BoolPolynomial.Bounded fresh value.polynomial) :
    (⟨st.wires.set q value, fresh⟩ : NonlinearAState).VariablesBounded := by
  intro r
  simp only [wireOf, List.getD_eq_getElem?_getD]
  by_cases hr : r = q
  · subst r
    by_cases hq : q < st.wires.length
    · rw [List.getElem?_set_self hq]
      exact hvalue
    · rw [List.getElem?_eq_none]
      · simpa [TrackedPolynomial.zero] using BoolPolynomial.bounded_zero fresh
      · simpa using Nat.le_of_not_lt hq
  · rw [List.getElem?_set, if_neg (by simpa using Ne.symm hr)]
    exact hst r

theorem variablesBounded_initial (n : Nat) : (initial n).VariablesBounded := by
  intro q
  simp only [wireOf, initial, getD_map_range]
  split
  · exact BoolPolynomial.bounded_var (by assumption)
  · simpa [TrackedPolynomial.zero] using BoolPolynomial.bounded_zero n

theorem variablesBounded_step {st : NonlinearAState} (hst : st.VariablesBounded) (g : Gate) :
    (st.step g).VariablesBounded := by
  cases g with
  | x q =>
      exact variablesBounded_set hst
        (BoolPolynomial.bounded_add (hst q) (by
          simpa [TrackedPolynomial.one] using BoolPolynomial.bounded_one st.fresh))
  | cnot c t =>
      exact variablesBounded_set hst (BoolPolynomial.bounded_add (hst t) (hst c))
  | h q =>
      have hmono : ∀ r, BoolPolynomial.Bounded (st.fresh + 1) (st.wireOf r).polynomial :=
        fun r => BoolPolynomial.bounded_mono (by omega) (hst r)
      exact variablesBounded_set hmono (BoolPolynomial.bounded_var (by omega))
  | ccx c₁ c₂ t =>
      simp only [step]
      cases hm : (st.wireOf c₁).mul? (st.wireOf c₂) with
      | some product =>
          have hp : product.polynomial =
              (st.wireOf c₁).polynomial * (st.wireOf c₂).polynomial := by
            simp only [TrackedPolynomial.mul?] at hm
            split at hm
            · cases hm; rfl
            · simp at hm
          exact variablesBounded_set hst
            (BoolPolynomial.bounded_add (hst t)
              (hp ▸ BoolPolynomial.bounded_mul (hst c₁) (hst c₂)))
      | none =>
          have hmono : ∀ r, BoolPolynomial.Bounded (st.fresh + 1)
              (st.wireOf r).polynomial :=
            fun r => BoolPolynomial.bounded_mono (by omega) (hst r)
          exact variablesBounded_set hmono (BoolPolynomial.bounded_var (by omega))
  | reset q =>
      exact variablesBounded_set hst (by
        simpa [TrackedPolynomial.zero] using BoolPolynomial.bounded_zero st.fresh)
  | _ => exact hst

/-- A valuation explains a basis state when each wire equals its formal polynomial on that
valuation. -/
def Consistent (n : Nat) (st : NonlinearAState) (valuation : Nat → Bool) (basis : Basis n) :
    Prop := ∀ q : Qubit, basis.get q = BoolPolynomial.evalB valuation (st.wireOf q).polynomial

theorem consistent_set {n : Nat} {st : NonlinearAState} {fresh : Nat} {q : Qubit}
    {value : TrackedPolynomial} {valuation : Nat → Bool} {basis : Basis n}
    (hlen : st.wires.length = n)
    (hq : q < n → basis.get q = BoolPolynomial.evalB valuation value.polynomial)
    (hother : ∀ r, r ≠ q →
      basis.get r = BoolPolynomial.evalB valuation (st.wireOf r).polynomial) :
    Consistent n ⟨st.wires.set q value, fresh⟩ valuation basis := by
  intro r
  simp only [wireOf, List.getD_eq_getElem?_getD]
  by_cases hr : r = q
  · subst r
    by_cases hlt : q < st.wires.length
    · rw [List.getElem?_set_self hlt]
      exact hq (hlen ▸ hlt)
    · rw [List.getElem?_eq_none]
      · change basis.get q = false
        exact basis_get_of_ge basis (hlen ▸ hlt)
      · simpa using Nat.le_of_not_lt hlt
  · rw [List.getElem?_set, if_neg (by simpa using Ne.symm hr)]
    exact hother r hr

theorem consistent_initial {n : Nat} (basis : Basis n) :
    Consistent n (initial n) (fun i => basis.get i) basis := by
  intro q
  simp only [initial, wireOf, getD_map_range]
  split
  · simp [TrackedPolynomial.fresh]
  · change basis.get q = false
    exact basis_get_of_ge basis (by assumption)

/-- Reset preserves the explanatory valuation while setting the analyzed wire to zero.
This is the support invariant needed after the analyzer ceases to be generic. -/
theorem consistent_reset {n : Nat} {st : NonlinearAState}
    (hlen : st.wires.length = n) {valuation : Nat → Bool} {before : Basis n}
    (hbefore : Consistent n st valuation before) (q : Qubit) :
    Consistent n (st.step (Gate.reset q)) valuation (before.set q false) := by
  change Consistent n
    ⟨st.wires.set q TrackedPolynomial.zero, st.fresh⟩ valuation (before.set q false)
  refine consistent_set hlen ?_ ?_
  · intro hq
    simp [Basis.get_set_same _ _ _ hq, TrackedPolynomial.zero,
      BoolPolynomial.evalB, BoolPolynomial.evalF₂, unbit]
  · intro r hr
    rw [Basis.get_set_ne _ _ _ _ hr]
    exact hbefore r

/-- Measurement does not alter the quantum wire polynomials. Each supported basis branch
therefore retains the same explanatory valuation. -/
theorem consistent_measure {n : Nat} {st : NonlinearAState}
    {valuation : Nat → Bool} {basis : Basis n}
    (hbefore : Consistent n st valuation basis) (q : Qubit) (c : CBit) :
    Consistent n (st.step (Gate.measure q c)) valuation basis := by
  simpa [NonlinearAState.step] using hbefore

/-- One nonlinear analysis step explains every nonzero unitary basis transition.  In the
ordinary CCX branch this uses `target + control₁ * control₂` directly; only an overflowing
CCX target is replaced by a fresh opaque variable. -/
theorem step_sound {n : Nat} {st : NonlinearAState} (hbounded : st.VariablesBounded)
    (hlen : st.wires.length = n) {g : Gate} (hunitary : g.isUnitary = true)
    {valuation : Nat → Bool} {before after : Basis n}
    (hbefore : Consistent n st valuation before)
    (hne : gateUnitary n g after before ≠ 0) :
    ∃ valuation', (∀ i, i < st.fresh → valuation' i = valuation i) ∧
      Consistent n (st.step g) valuation' after := by
  have diagonalCase : ∀ (U : Density n), IsDiagonal U → gateUnitary n g = U → st.step g = st →
      ∃ valuation', (∀ i, i < st.fresh → valuation' i = valuation i) ∧
        Consistent n (st.step g) valuation' after := by
    intro U hU hgU hstep
    refine ⟨valuation, fun _ _ => rfl, ?_⟩
    have hsame : after = before := AState.diag_entry hU (by rw [← hgU]; exact hne)
    rw [hstep, hsame]
    exact hbefore
  cases g with
  | x q =>
      refine ⟨valuation, fun _ _ => rfl, ?_⟩
      have hafter : after = before.set q (!before.get q) :=
        AState.perm_entry (by rw [← embed1_x2_eq_perm]; exact hne)
      refine consistent_set hlen ?_ ?_
      · intro hqn
        rw [hafter, Basis.get_set_same _ _ _ hqn, TrackedPolynomial.evalB_add,
          TrackedPolynomial.evalB_one, ← hbefore q]
        simp
      · intro r hr
        rw [hafter, Basis.get_set_ne _ _ _ _ hr, hbefore r]
  | cnot c t =>
      refine ⟨valuation, fun _ _ => rfl, ?_⟩
      have hafter : after = before.set t (before.get t != before.get c) :=
        AState.perm_entry
          (n := n) (b' := after) (b := before)
          (σ := fun b : Basis n => b.set t (b.get t != b.get c)) hne
      refine consistent_set hlen ?_ ?_
      · intro htn
        rw [hafter, Basis.get_set_same _ _ _ htn, TrackedPolynomial.evalB_add,
          ← hbefore t, ← hbefore c]
      · intro r hr
        rw [hafter, Basis.get_set_ne _ _ _ _ hr, hbefore r]
  | ccx c₁ c₂ t =>
      have hafter : after = before.set t
          (before.get t != (before.get c₁ && before.get c₂)) :=
        AState.perm_entry
          (n := n) (b' := after) (b := before)
          (σ := fun b : Basis n => b.set t (b.get t != (b.get c₁ && b.get c₂))) hne
      cases hm : (st.wireOf c₁).mul? (st.wireOf c₂) with
      | some product =>
          refine ⟨valuation, fun _ _ => rfl, ?_⟩
          simp only [step, hm]
          refine consistent_set hlen ?_ ?_
          · intro htn
            have hp : product.polynomial =
                (st.wireOf c₁).polynomial * (st.wireOf c₂).polynomial := by
              simp only [TrackedPolynomial.mul?] at hm
              split at hm
              · cases hm; rfl
              · simp at hm
            rw [hafter, Basis.get_set_same _ _ _ htn,
              TrackedPolynomial.evalB_add, hp, BoolPolynomial.evalB_mul, ← hbefore t,
              ← hbefore c₁, ← hbefore c₂]
          · intro r hr
            rw [hafter, Basis.get_set_ne _ _ _ _ hr, hbefore r]
      | none =>
          let valuation' := Function.update valuation st.fresh (after.get t)
          refine ⟨valuation', fun i hi => Function.update_of_ne (by omega) _ _, ?_⟩
          simp only [step, hm]
          refine consistent_set hlen ?_ ?_
          · intro _
            simp [valuation', TrackedPolynomial.fresh]
          · intro r hr
            rw [hafter, Basis.get_set_ne _ _ _ _ hr, hbefore r]
            exact (BoolPolynomial.evalB_congr (hbounded r)
              (fun i hi => Function.update_of_ne (by omega) _ _)).symm
  | h q =>
      let valuation' := Function.update valuation st.fresh (after.get q)
      refine ⟨valuation', fun i hi => Function.update_of_ne (by omega) _ _, ?_⟩
      refine consistent_set hlen ?_ ?_
      · intro _
        simp [valuation', TrackedPolynomial.fresh]
      · intro r hr
        rw [AState.embed1_entry hne hr, hbefore r]
        exact (BoolPolynomial.evalB_congr (hbounded r)
          (fun i hi => Function.update_of_ne (by omega) _ _)).symm
  | s q => exact diagonalCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | sdg q => exact diagonalCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | z q => exact diagonalCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | t q => exact diagonalCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | tdg q => exact diagonalCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | rz angle q => exact diagonalCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | cz c t => exact diagonalCase _ (isDiagonal_phaseMatrix _) rfl rfl
  | ccz c₁ c₂ t => exact diagonalCase _ (isDiagonal_phaseMatrix _) rfl rfl
  | measure q c => simp [Gate.isUnitary, Gate.isMeasurement] at hunitary
  | reset q => simp [Gate.isUnitary, Gate.isMeasurement] at hunitary

theorem variablesBounded_steps {st : NonlinearAState} (hbounded : st.VariablesBounded) :
    ∀ gates, (st.steps gates).VariablesBounded := by
  intro gates
  induction gates generalizing st with
  | nil => exact hbounded
  | cons gate gates ih => exact ih (variablesBounded_step hbounded gate)

theorem length_steps (st : NonlinearAState) (gates : List Gate) :
    (st.steps gates).wires.length = st.wires.length := by
  induction gates generalizing st with
  | nil => rfl
  | cons gate gates ih => rw [steps, ih, length_step]

/-- The nonlinear polynomial analysis explains every nonzero transition through a unitary
fragment. This is the circuit-level CCX simulation theorem needed by phase folding. -/
theorem analyze_sound {n : Nat} {gates : List Gate}
    (hunitary : ∀ gate ∈ gates, gate.isUnitary = true)
    {st : NonlinearAState} (hbounded : st.VariablesBounded)
    (hlen : st.wires.length = n) {valuation : Nat → Bool} {before after : Basis n}
    (hbefore : Consistent n st valuation before)
    (hne : unitary n gates after before ≠ 0) :
    ∃ valuation', (∀ i, i < st.fresh → valuation' i = valuation i) ∧
      Consistent n (st.steps gates) valuation' after := by
  induction gates generalizing st valuation before with
  | nil =>
      have hsame : after = before := AState.one_entry (by simpa using hne)
      refine ⟨valuation, fun _ _ => rfl, ?_⟩
      rw [steps, hsame]
      exact hbefore
  | cons gate gates ih =>
      rw [unitary_cons, Matrix.mul_apply] at hne
      obtain ⟨middle, -, hmiddle⟩ := Finset.exists_ne_zero_of_sum_ne_zero hne
      have hsuffix : unitary n gates after middle ≠ 0 := left_ne_zero_of_mul hmiddle
      have hgate : gateUnitary n gate middle before ≠ 0 := right_ne_zero_of_mul hmiddle
      obtain ⟨valuation₁, hvaluation₁, hconsistent₁⟩ :=
        step_sound hbounded hlen (hunitary gate (by simp)) hbefore hgate
      obtain ⟨valuation', hvaluation', hconsistent'⟩ :=
        ih (fun later hlater => hunitary later (by simp [hlater]))
          (variablesBounded_step hbounded gate) (by rw [length_step, hlen])
          hconsistent₁ hsuffix
      refine ⟨valuation', fun i hi => ?_, ?_⟩
      rw [hvaluation' i (lt_of_lt_of_le hi (fresh_le_step st gate)), hvaluation₁ i hi]
      change Consistent n ((st.step gate).steps gates) valuation' after
      exact hconsistent'

/-! ## Branch paths through measurement and reset -/

/-- A possible basis-state branch of one gate. Measurement keeps the measured basis value;
reset discards it and writes zero. Other gates use their nonzero matrix entries. -/
def GatePossible (n : Nat) (gate : Gate) (before after : Basis n) : Prop :=
  match gate with
  | .measure _ _ => after = before
  | .reset q => after = before.set q false
  | _ => gateUnitary n gate after before ≠ 0

theorem measurementKraus_possible {n : Nat} (q : Qubit) (c : CBit) (b : Bool)
    {before after : Basis n} (hne : proj n q b after before ≠ 0) :
    GatePossible n (Gate.measure q c) before after := by
  change after = before
  exact AState.diag_entry (isDiagonal_phaseMatrix _) hne

theorem resetKraus_possible {n : Nat} (q : Qubit) (b : Bool)
    {before after : Basis n} (hne : resetKraus n q b after before ≠ 0) :
    GatePossible n (Gate.reset q) before after := by
  change after = before.set q false
  simp only [resetKraus] at hne
  split at hne
  · exact (by assumption : before.get q = b ∧ after = before.set q false).2
  · simp at hne

/-- A basis-state branch through a gate list, abstracting away the classical memory value. -/
inductive NonlinearPath (n : Nat) : List Gate → Basis n → Basis n → Prop where
  | nil (basis : Basis n) : NonlinearPath n [] basis basis
  | cons {gate gates before middle after} :
      GatePossible n gate before middle →
      NonlinearPath n gates middle after →
      NonlinearPath n (gate :: gates) before after

theorem step_sound_possible {n : Nat} {st : NonlinearAState}
    (hbounded : st.VariablesBounded) (hlen : st.wires.length = n)
    {gate : Gate} {valuation : Nat → Bool} {before after : Basis n}
    (hbefore : Consistent n st valuation before)
    (hpossible : GatePossible n gate before after) :
    ∃ valuation', (∀ i, i < st.fresh → valuation' i = valuation i) ∧
      Consistent n (st.step gate) valuation' after := by
  cases gate with
  | measure q c =>
      change after = before at hpossible
      subst after
      exact ⟨valuation, fun _ _ => rfl, consistent_measure hbefore q c⟩
  | reset q =>
      change after = before.set q false at hpossible
      subst after
      exact ⟨valuation, fun _ _ => rfl, consistent_reset hlen hbefore q⟩
  | h q => exact step_sound hbounded hlen rfl hbefore hpossible
  | x q => exact step_sound hbounded hlen rfl hbefore hpossible
  | z q => exact step_sound hbounded hlen rfl hbefore hpossible
  | s q => exact step_sound hbounded hlen rfl hbefore hpossible
  | sdg q => exact step_sound hbounded hlen rfl hbefore hpossible
  | t q => exact step_sound hbounded hlen rfl hbefore hpossible
  | tdg q => exact step_sound hbounded hlen rfl hbefore hpossible
  | rz angle q => exact step_sound hbounded hlen rfl hbefore hpossible
  | cnot c t => exact step_sound hbounded hlen rfl hbefore hpossible
  | cz c t => exact step_sound hbounded hlen rfl hbefore hpossible
  | ccx c₁ c₂ t => exact step_sound hbounded hlen rfl hbefore hpossible
  | ccz c₁ c₂ t => exact step_sound hbounded hlen rfl hbefore hpossible

/-- The formal nonlinear analyzer explains every possible branch of the entire circuit,
including CCX, measurement, and reset. Unlike `Generic`, this invariant survives reset. -/
theorem analyze_path_sound {n : Nat} {gates : List Gate}
    {st : NonlinearAState} (hbounded : st.VariablesBounded)
    (hlen : st.wires.length = n) {valuation : Nat → Bool}
    {before after : Basis n}
    (hbefore : Consistent n st valuation before)
    (hpath : NonlinearPath n gates before after) :
    ∃ valuation', (∀ i, i < st.fresh → valuation' i = valuation i) ∧
      Consistent n (st.steps gates) valuation' after := by
  induction hpath generalizing st valuation with
  | nil _ =>
      exact ⟨valuation, fun _ _ => rfl, hbefore⟩
  | @cons gate gates before middle after hgate htail ih =>
      obtain ⟨valuation₁, hagree₁, hmiddle⟩ :=
        step_sound_possible hbounded hlen hbefore hgate
      obtain ⟨valuation₂, hagree₂, hafter⟩ :=
        ih (variablesBounded_step hbounded gate) (by rw [length_step, hlen]) hmiddle
      refine ⟨valuation₂, ?_, ?_⟩
      · intro i hi
        rw [hagree₂ i (lt_of_lt_of_le hi (fresh_le_step st gate)), hagree₁ i hi]
      · exact hafter

/-- Formal equality (or complement) transports a wire bit along any supported branch,
including branches containing measurement or reset. The starting branch must be consistent
with the symbolic state; no global genericity assumption is needed. -/
theorem path_of_polynomial_eq_supported {n : Nat} {gates : List Gate}
    {st : NonlinearAState} (hbounded : st.VariablesBounded)
    (hlen : st.wires.length = n)
    {valuation : Nat → Bool} {before after : Basis n}
    (hbefore : Consistent n st valuation before)
    {q q' : Qubit} {sign : Bool}
    (hpoly : ((st.steps gates).wireOf q').polynomial =
      if sign then (st.wireOf q).polynomial.flip else (st.wireOf q).polynomial)
    (hpath : NonlinearPath n gates before after) :
    before.get q = (after.get q' != sign) := by
  obtain ⟨valuation', hagree, hafter⟩ :=
    analyze_path_sound hbounded hlen hbefore hpath
  have hbeforeq := hbefore q
  have hafterq := hafter q'
  rw [hpoly] at hafterq
  have hsame := BoolPolynomial.evalB_congr (hbounded q) hagree
  cases sign with
  | false =>
      simp only [Bool.false_eq_true, ↓reduceIte] at hafterq
      rw [hbeforeq, ← hsame, ← hafterq]
      simp
  | true =>
      simp only [↓reduceIte, BoolPolynomial.evalB_flip] at hafterq
      rw [hbeforeq, ← hsame, hafterq]
      simp

/-- Every branch from a circuit input has a Boolean valuation explaining its final
symbolic state, even after reset has destroyed genericity. -/
theorem path_from_initial_consistent {n : Nat} {gates : List Gate}
    {before after : Basis n} (hpath : NonlinearPath n gates before after) :
    ∃ valuation : Nat → Bool,
      Consistent n ((initial n).steps gates) valuation after := by
  obtain ⟨valuation, -, hafter⟩ := analyze_path_sound
    (variablesBounded_initial n) (by simp [initial])
    (consistent_initial before) hpath
  exact ⟨valuation, hafter⟩

/-- The branch-sensitive phase-hop premise for an arbitrary circuit prefix. In particular,
the prefix may contain reset, and the middle segment may contain CCX or measurement. -/
theorem path_of_polynomial_eq_after_prefix {n : Nat}
    {pre middle : List Gate} {input before after : Basis n}
    (hprefix : NonlinearPath n pre input before)
    (hmiddle : NonlinearPath n middle before after)
    {q q' : Qubit} {sign : Bool}
    (hpoly : ((((initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((initial n).steps pre).wireOf q).polynomial.flip
      else (((initial n).steps pre).wireOf q).polynomial) :
    before.get q = (after.get q' != sign) := by
  obtain ⟨valuation, hbefore⟩ := path_from_initial_consistent hprefix
  exact path_of_polynomial_eq_supported
    (variablesBounded_steps (variablesBounded_initial n) pre)
    (by rw [length_steps]; simp [initial]) hbefore hpoly hmiddle

/-- Every basis state can be described by some valuation of a generic polynomial state. -/
def Generic (n : Nat) (st : NonlinearAState) : Prop :=
  ∀ basis : Basis n, ∃ valuation, Consistent n st valuation basis

theorem generic_initial (n : Nat) : Generic n (initial n) :=
  fun basis => ⟨fun i => basis.get i, consistent_initial basis⟩

private theorem generic_fresh {n : Nat} {st : NonlinearAState}
    (hbounded : st.VariablesBounded) (hlen : st.wires.length = n)
    (hgeneric : Generic n st) (q : Qubit) :
    Generic n ⟨st.wires.set q (TrackedPolynomial.fresh st.fresh), st.fresh + 1⟩ := by
  intro after
  obtain ⟨valuation, hbefore⟩ := hgeneric after
  let valuation' := Function.update valuation st.fresh (after.get q)
  refine ⟨valuation', consistent_set hlen ?_ ?_⟩
  · intro _
    simp [valuation', TrackedPolynomial.fresh]
  · intro r _
    rw [hbefore r]
    exact (BoolPolynomial.evalB_congr (hbounded r)
      (fun i hi => Function.update_of_ne (by omega) _ _)).symm

/-- Unitary nonlinear transfers preserve coverage of all basis states. Measurement leaves
the wire polynomials unchanged. Reset is excluded because it fixes one wire to zero. -/
theorem generic_step {n : Nat} {st : NonlinearAState}
    (hbounded : st.VariablesBounded) (hlen : st.wires.length = n)
    (hgeneric : Generic n st) {gate : Gate} (hwf : gate.Wf)
    (hreset : ∀ q, gate ≠ Gate.reset q) : Generic n (st.step gate) := by
  cases gate with
  | x q =>
      intro after
      obtain ⟨valuation, hbefore⟩ := hgeneric (after.set q (!after.get q))
      refine ⟨valuation, consistent_set hlen ?_ ?_⟩
      · intro hq
        rw [TrackedPolynomial.evalB_add, TrackedPolynomial.evalB_one,
          ← hbefore q, Basis.get_set_same _ _ _ hq]
        simp
      · intro r hr
        rw [← hbefore r, Basis.get_set_ne _ _ _ _ hr]
  | cnot c t =>
      intro after
      obtain ⟨valuation, hbefore⟩ :=
        hgeneric (after.set t (after.get t != after.get c))
      have hct : c ≠ t := by simpa [Gate.Wf] using hwf
      refine ⟨valuation, consistent_set hlen ?_ ?_⟩
      · intro ht
        rw [TrackedPolynomial.evalB_add, ← hbefore t, ← hbefore c,
          Basis.get_set_same _ _ _ ht, Basis.get_set_ne _ _ _ _ hct]
        cases after.get t <;> cases after.get c <;> rfl
      · intro r hr
        rw [← hbefore r, Basis.get_set_ne _ _ _ _ hr]
  | h q => exact generic_fresh hbounded hlen hgeneric q
  | ccx c₁ c₂ t =>
      cases hm : (st.wireOf c₁).mul? (st.wireOf c₂) with
      | none =>
          simpa [step, hm] using generic_fresh hbounded hlen hgeneric t
      | some product =>
          intro after
          obtain ⟨valuation, hbefore⟩ := hgeneric
            (after.set t (after.get t != (after.get c₁ && after.get c₂)))
          have hct₁ : c₁ ≠ t := (show c₁ ≠ c₂ ∧ c₁ ≠ t ∧ c₂ ≠ t from hwf).2.1
          have hct₂ : c₂ ≠ t := (show c₁ ≠ c₂ ∧ c₁ ≠ t ∧ c₂ ≠ t from hwf).2.2
          simp only [step, hm]
          refine ⟨valuation, consistent_set hlen ?_ ?_⟩
          · intro ht
            have hp : product.polynomial =
                (st.wireOf c₁).polynomial * (st.wireOf c₂).polynomial := by
              simp only [TrackedPolynomial.mul?] at hm
              split at hm
              · cases hm; rfl
              · simp at hm
            rw [TrackedPolynomial.evalB_add, hp, BoolPolynomial.evalB_mul,
              ← hbefore t, ← hbefore c₁, ← hbefore c₂,
              Basis.get_set_same _ _ _ ht,
              Basis.get_set_ne _ _ _ _ hct₁,
              Basis.get_set_ne _ _ _ _ hct₂]
            cases after.get t <;> cases after.get c₁ <;> cases after.get c₂ <;> rfl
          · intro r hr
            rw [← hbefore r, Basis.get_set_ne _ _ _ _ hr]
  | reset q => exact False.elim (hreset q rfl)
  | _ => exact hgeneric

theorem generic_steps {n : Nat} {st : NonlinearAState}
    (hbounded : st.VariablesBounded) (hlen : st.wires.length = n)
    (hgeneric : Generic n st) {gates : List Gate}
    (hwf : ∀ gate ∈ gates, gate.Wf)
    (hreset : ∀ gate ∈ gates, ∀ q, gate ≠ Gate.reset q) :
    Generic n (st.steps gates) := by
  induction gates generalizing st with
  | nil => exact hgeneric
  | cons gate gates ih =>
      exact ih (variablesBounded_step hbounded gate)
        (by rw [length_step, hlen])
        (generic_step hbounded hlen hgeneric (hwf gate (by simp))
          (hreset gate (by simp)))
        (fun later hlater => hwf later (by simp [hlater]))
        (fun later hlater => hreset later (by simp [hlater]))

/-- Equality of the formal wire polynomials transports the wire value along every
nonzero path through a unitary segment, including segments containing CCX. -/
theorem path_of_polynomial_eq {n : Nat} {gates : List Gate}
    (hunitary : ∀ gate ∈ gates, gate.isUnitary = true)
    {st : NonlinearAState} (hbounded : st.VariablesBounded)
    (hlen : st.wires.length = n) (hgeneric : Generic n st)
    {q q' : Qubit} {sign : Bool}
    (hpoly : ((st.steps gates).wireOf q').polynomial =
      if sign then (st.wireOf q).polynomial.flip else (st.wireOf q).polynomial) :
    ∀ before after : Basis n, unitary n gates after before ≠ 0 →
      before.get q = (after.get q' != sign) := by
  intro before after hpath
  obtain ⟨valuation, hbefore⟩ := hgeneric before
  obtain ⟨valuation', hagree, hafter⟩ :=
    analyze_sound hunitary hbounded hlen hbefore hpath
  have hbeforeq := hbefore q
  have hafterq := hafter q'
  rw [hpoly] at hafterq
  have hsame := BoolPolynomial.evalB_congr (hbounded q) hagree
  cases sign with
  | false =>
      simp only [Bool.false_eq_true, ↓reduceIte] at hafterq
      rw [hbeforeq, ← hsame, ← hafterq]
      simp
  | true =>
      simp only [↓reduceIte, BoolPolynomial.evalB_flip] at hafterq
      rw [hbeforeq, ← hsame, hafterq]
      simp

end NonlinearAState

end

end TzapLean
