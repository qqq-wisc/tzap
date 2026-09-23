import TzapLean.NonlinearRefinement

/-!
# Polynomials compared by nonlinear phase folding

The finite comparison set contains both constants, every wire polynomial at every circuit
position, and the complements used by `matchFingerprint`. The degree and variable bounds
below are the inputs needed for a whole-run Schwartz–Zippel union bound.
-/

namespace TzapLean

noncomputable section

def nonlinearPolynomialsOf (n : Nat) (st : NonlinearAState) : List BoolPolynomial :=
  0 :: 1 :: (List.range n).map (fun q => (st.wireOf q).polynomial)

def nonlinearVisited (n : Nat) (st : NonlinearAState) : List Gate → List BoolPolynomial
  | [] => nonlinearPolynomialsOf n st
  | gate :: gates =>
      nonlinearPolynomialsOf n st ++ nonlinearVisited n (st.step gate) gates

def nonlinearRelevant (n : Nat) (st : NonlinearAState) (gates : List Gate) :
    List BoolPolynomial :=
  nonlinearVisited n st gates ++ (nonlinearVisited n st gates).map BoolPolynomial.flip

theorem nonlinearWire_mem_polynomialsOf {n : Nat} {st : NonlinearAState}
    (hlen : st.wires.length = n) (q : Qubit) :
    (st.wireOf q).polynomial ∈ nonlinearPolynomialsOf n st := by
  by_cases hq : q < n
  · exact List.mem_cons_of_mem _ (List.mem_cons_of_mem _
      (List.mem_map.2 ⟨q, List.mem_range.2 hq, rfl⟩))
  · have hzero : (st.wireOf q).polynomial = 0 := by
      have hge : st.wires.length ≤ q := by rw [hlen]; exact Nat.le_of_not_lt hq
      simp [NonlinearAState.wireOf, List.getD_eq_getElem?_getD,
        List.getElem?_eq_none hge, TrackedPolynomial.zero]
    rw [hzero]
    exact List.mem_cons_self

theorem nonlinearWire_mem_visited {n : Nat} {st : NonlinearAState}
    (hlen : st.wires.length = n) (q : Qubit) (gates : List Gate) :
    (st.wireOf q).polynomial ∈ nonlinearVisited n st gates := by
  cases gates with
  | nil => exact nonlinearWire_mem_polynomialsOf hlen q
  | cons gate gates =>
      exact List.mem_append.2 (Or.inl (nonlinearWire_mem_polynomialsOf hlen q))

theorem nonlinearVisited_append_mem {n : Nat} {st : NonlinearAState}
    {polynomial : BoolPolynomial} (pre rest : List Gate)
    (h : polynomial ∈ nonlinearVisited n (st.steps pre) rest) :
    polynomial ∈ nonlinearVisited n st (pre ++ rest) := by
  induction pre generalizing st with
  | nil => exact h
  | cons gate pre ih => exact List.mem_append.2 (Or.inr (ih h))

theorem nonlinear_step_of_rotAngle {st : NonlinearAState} {gate : Gate}
    {θ : ℚ} {q : Qubit} (hrot : rotAngle gate = some (θ, q)) :
    st.step gate = st := by
  cases gate <;> simp [rotAngle] at hrot ⊢ <;> rfl

theorem nonlinearVisited_append_sub (n : Nat) (pre before after : List Gate)
    (st : NonlinearAState)
    (hsub : ∀ p ∈ nonlinearVisited n (st.steps pre) before,
      p ∈ nonlinearVisited n (st.steps pre) after) :
    ∀ p ∈ nonlinearVisited n st (pre ++ before),
      p ∈ nonlinearVisited n st (pre ++ after) := by
  induction pre generalizing st with
  | nil => exact hsub
  | cons gate pre ih =>
      intro p hp
      rcases List.mem_append.1 hp with hp | hp
      · exact List.mem_append.2 (Or.inl hp)
      · exact List.mem_append.2 (Or.inr (ih (st.step gate) hsub p hp))

theorem nonlinearVisited_cons_rot {n : Nat} {st : NonlinearAState}
    {gate gate' : Gate} {θ θ' : ℚ} {q q' : Qubit}
    (hrot : rotAngle gate = some (θ, q))
    (hrot' : rotAngle gate' = some (θ', q')) (gates : List Gate) :
    nonlinearVisited n st (gate :: gates) = nonlinearVisited n st (gate' :: gates) := by
  rw [nonlinearVisited, nonlinearVisited, nonlinear_step_of_rotAngle hrot,
    nonlinear_step_of_rotAngle hrot']

/-- Replacing a later rotation does not introduce a new comparison polynomial, even
when the intervening segment contains nonlinear CCX transfers. -/
theorem nonlinearVisited_merged_subset {n : Nat} (st : NonlinearAState)
    (first later : Gate) (middle rest : List Gate)
    (θ φ : ℚ) (q q' : Qubit) (sign : Bool)
    (hfirst : rotAngle first = some (θ, q))
    (hlater : rotAngle later = some (φ, q')) :
    ∀ p ∈ nonlinearVisited n st
      (middle ++ Gate.rz (φ + signedAngle sign θ) q' :: rest),
      p ∈ nonlinearVisited n st (first :: middle ++ later :: rest) := by
  intro p hp
  refine List.mem_append.2 (Or.inr ?_)
  rw [nonlinear_step_of_rotAngle hfirst]
  refine nonlinearVisited_append_sub n middle
    (Gate.rz (φ + signedAngle sign θ) q' :: rest)
    (later :: rest) st ?_ p hp
  intro r hr
  rw [nonlinearVisited_cons_rot (gates := rest) (hrot := rfl) (hrot' := hlater)] at hr
  exact hr

theorem nonlinearFaithful_of_sub {draws : Nat → Tag}
    {evaluation : PackedEvaluation draws} {n : Nat}
    {st st' : NonlinearAState} {gates smaller : List Gate}
    (hfaithful : PolynomialFaithful evaluation (nonlinearRelevant n st gates))
    (hsub : ∀ p ∈ nonlinearVisited n st' smaller,
      p ∈ nonlinearVisited n st gates) :
    PolynomialFaithful evaluation (nonlinearRelevant n st' smaller) := by
  intro p hp q hq heq
  apply hfaithful p ?_ q ?_ heq
  · rcases List.mem_append.1 hp with hp | hp
    · exact List.mem_append.2 (Or.inl (hsub p hp))
    · rcases List.mem_map.1 hp with ⟨r, hr, rfl⟩
      exact List.mem_append.2 (Or.inr (List.mem_map.2 ⟨r, hsub r hr, rfl⟩))
  · rcases List.mem_append.1 hq with hq | hq
    · exact List.mem_append.2 (Or.inl (hsub q hq))
    · rcases List.mem_map.1 hq with ⟨r, hr, rfl⟩
      exact List.mem_append.2 (Or.inr (List.mem_map.2 ⟨r, hsub r hr, rfl⟩))

theorem nonlinearPolynomialsOf_variablesBounded {n m : Nat} {st : NonlinearAState}
    (hbounded : st.VariablesBounded) (hfresh : st.fresh ≤ m) :
    ∀ polynomial ∈ nonlinearPolynomialsOf n st,
      BoolPolynomial.Bounded m polynomial := by
  intro polynomial hmem
  rcases List.mem_cons.1 hmem with rfl | hmem
  · exact BoolPolynomial.bounded_zero m
  rcases List.mem_cons.1 hmem with rfl | hmem
  · exact BoolPolynomial.bounded_one m
  rcases List.mem_map.1 hmem with ⟨q, -, rfl⟩
  exact BoolPolynomial.bounded_mono hfresh (hbounded q)

theorem nonlinearPolynomialsOf_degreeBounded {n : Nat} {st : NonlinearAState}
    (hdegree : st.DegreeBounded) (hcutoff : st.WithinCutoff) :
    ∀ polynomial ∈ nonlinearPolynomialsOf n st,
      polynomial.totalDegree ≤ Fingerprint.maxDegree := by
  intro polynomial hmem
  rcases List.mem_cons.1 hmem with rfl | hmem
  · simp [Fingerprint.maxDegree]
  rcases List.mem_cons.1 hmem with rfl | hmem
  · simp [Fingerprint.maxDegree]
  rcases List.mem_map.1 hmem with ⟨q, -, rfl⟩
  exact TrackedPolynomial.totalDegree_le_cutoff (hdegree q) (hcutoff q)

theorem nonlinear_step_fresh_le (st : NonlinearAState) (gate : Gate) :
    (st.step gate).fresh ≤ st.fresh + if Gate.allocates gate then 1 else 0 := by
  cases gate <;> simp [NonlinearAState.step, Gate.allocates]
  all_goals split <;> simp

theorem nonlinearVisited_variablesBounded {n m : Nat} :
    ∀ (gates : List Gate) (st : NonlinearAState), st.VariablesBounded →
      st.fresh + gates.countP Gate.allocates ≤ m →
      ∀ polynomial ∈ nonlinearVisited n st gates,
        BoolPolynomial.Bounded m polynomial := by
  intro gates
  induction gates with
  | nil =>
      intro st hbounded hfresh polynomial hmem
      exact nonlinearPolynomialsOf_variablesBounded hbounded (by simpa using hfresh)
        polynomial hmem
  | cons gate gates ih =>
      intro st hbounded hfresh polynomial hmem
      rw [nonlinearVisited] at hmem
      rcases List.mem_append.1 hmem with hmem | hmem
      · exact nonlinearPolynomialsOf_variablesBounded hbounded (by
          simp only [List.countP_cons] at hfresh
          omega) polynomial hmem
      · exact ih (st.step gate) (NonlinearAState.variablesBounded_step hbounded gate)
          (by
            have hstep := nonlinear_step_fresh_le st gate
            simp only [List.countP_cons] at hfresh
            omega) polynomial hmem

theorem nonlinearVisited_degreeBounded {n : Nat} :
    ∀ (gates : List Gate) (st : NonlinearAState),
      st.DegreeBounded → st.WithinCutoff →
      ∀ polynomial ∈ nonlinearVisited n st gates,
        polynomial.totalDegree ≤ Fingerprint.maxDegree := by
  intro gates
  induction gates with
  | nil =>
      intro st hdegree hcutoff polynomial hmem
      exact nonlinearPolynomialsOf_degreeBounded hdegree hcutoff polynomial hmem
  | cons gate gates ih =>
      intro st hdegree hcutoff polynomial hmem
      rw [nonlinearVisited] at hmem
      rcases List.mem_append.1 hmem with hmem | hmem
      · exact nonlinearPolynomialsOf_degreeBounded hdegree hcutoff polynomial hmem
      · exact ih (st.step gate) (NonlinearAState.degreeBounded_step hdegree gate)
          (NonlinearAState.withinCutoff_step hcutoff gate) polynomial hmem

theorem nonlinearRelevant_variablesBounded {n m : Nat} {st : NonlinearAState}
    {gates : List Gate} (hbounded : st.VariablesBounded)
    (hfresh : st.fresh + gates.countP Gate.allocates ≤ m) :
    ∀ polynomial ∈ nonlinearRelevant n st gates,
      BoolPolynomial.Bounded m polynomial := by
  intro polynomial hmem
  rcases List.mem_append.1 hmem with hmem | hmem
  · exact nonlinearVisited_variablesBounded gates st hbounded hfresh polynomial hmem
  · rcases List.mem_map.1 hmem with ⟨prior, hprior, rfl⟩
    exact BoolPolynomial.bounded_add
      (nonlinearVisited_variablesBounded gates st hbounded hfresh prior hprior)
      (BoolPolynomial.bounded_one m)

theorem nonlinearRelevant_degreeBounded {n : Nat} {st : NonlinearAState}
    {gates : List Gate} (hdegree : st.DegreeBounded) (hcutoff : st.WithinCutoff) :
    ∀ polynomial ∈ nonlinearRelevant n st gates,
      polynomial.totalDegree ≤ Fingerprint.maxDegree := by
  intro polynomial hmem
  rcases List.mem_append.1 hmem with hmem | hmem
  · exact nonlinearVisited_degreeBounded gates st hdegree hcutoff polynomial hmem
  · rcases List.mem_map.1 hmem with ⟨prior, hprior, rfl⟩
    rw [BoolPolynomial.flip]
    exact (BoolPolynomial.totalDegree_add_le prior 1).trans (max_le
      (nonlinearVisited_degreeBounded gates st hdegree hcutoff prior hprior)
      (by simp [Fingerprint.maxDegree]))

/-- A seed-independent finite superset of every polynomial comparison in one circuit run. -/
def nonlinearRelevantCircuit (circuit : RawCircuit) : List BoolPolynomial :=
  nonlinearRelevant circuit.numQubits
    (NonlinearAState.initial circuit.numQubits) circuit.gates

theorem nonlinearRelevantCircuit_variablesBounded (circuit : RawCircuit) :
    ∀ polynomial ∈ nonlinearRelevantCircuit circuit,
      BoolPolynomial.Bounded (varBound circuit) polynomial := by
  exact nonlinearRelevant_variablesBounded
    (NonlinearAState.variablesBounded_initial circuit.numQubits)
    (by simp [varBound, NonlinearAState.initial])

theorem nonlinearRelevantCircuit_degreeBounded (circuit : RawCircuit) :
    ∀ polynomial ∈ nonlinearRelevantCircuit circuit,
      polynomial.totalDegree ≤ Fingerprint.maxDegree := by
  exact nonlinearRelevant_degreeBounded
    (NonlinearAState.degreeBounded_initial circuit.numQubits)
    (NonlinearAState.withinCutoff_initial circuit.numQubits)

end

end TzapLean
