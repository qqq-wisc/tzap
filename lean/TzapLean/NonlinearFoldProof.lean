import TzapLean.NonlinearChannel
import TzapLean.NonlinearMerge
import TzapLean.PhaseFoldRand

/-!
# Nonlinear executable fold correctness

The channel-level phase-hop and constant-phase rules are applied here to comparisons
made by the executable packed-state scan. The packed-field interpretation and the
whole-run induction remain separate obligations.
-/

namespace TzapLean

noncomputable section

theorem nonlinear_steps_append (st : NonlinearAState) (pre post : List Gate) :
    st.steps (pre ++ post) = (st.steps pre).steps post := by
  induction pre generalizing st with
  | nil => rfl
  | cons gate pre ih => exact ih (st.step gate)

theorem nonlinear_runtime_step_of_rotAngle {st : NState} {gate : Gate}
    {θ : ℚ} {q : Qubit} (draws : Nat → Tag)
    (hrot : rotAngle gate = some (θ, q)) : st.step draws gate = st := by
  cases gate <;> simp [rotAngle] at hrot ⊢ <;> rfl

theorem nonlinear_constant_of_faithful_tag {n : Nat} (gates : List Gate)
    (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (st : NonlinearAState) (runtime : NState) (q : Qubit) (bit : Bool)
    (hsim : NonlinearSim draws evaluation st runtime)
    (hlen : st.wires.length = n)
    (hfaithful : PolynomialFaithful evaluation (nonlinearRelevant n st gates))
    (htag : runtime.tagOf q = if bit then 1 else 0) :
    (st.wireOf q).polynomial = BoolPolynomial.const bit := by
  have hp : (st.wireOf q).polynomial ∈ nonlinearRelevant n st gates :=
    List.mem_append.2 (Or.inl (nonlinearWire_mem_visited hlen q gates))
  have hc : BoolPolynomial.const bit ∈ nonlinearRelevant n st gates := by
    apply List.mem_append.2
    apply Or.inl
    cases gates <;> cases bit <;>
      simp [nonlinearVisited, nonlinearPolynomialsOf, BoolPolynomial.const, bit]
  apply hfaithful _ hp _ hc
  cases bit <;> simp [BoolPolynomial.const, bit, evaluation.zero,
    evaluation.one, ← nonlinearSim_tagOf hsim q] at htag ⊢
  all_goals exact htag

theorem nonlinear_drop_constant_channel_of_tag {n m : Nat}
    (pre rest : List Gate) (first : Gate) (θ : ℚ) (q : Qubit) (bit : Bool)
    (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (st : NonlinearAState) (runtime : NState)
    (hprefix : st = (NonlinearAState.initial n).steps pre)
    (hsim : NonlinearSim draws evaluation st runtime)
    (hlen : st.wires.length = n)
    (hfirst : rotAngle first = some (θ, q))
    (hfaithful : PolynomialFaithful evaluation
      (nonlinearRelevant n st (first :: rest)))
    (htag : runtime.tagOf q = if bit then 1 else 0) :
    Equivalent n m (pre ++ first :: rest) (pre ++ rest) := by
  have hpoly := nonlinear_constant_of_faithful_tag (first :: rest) draws evaluation
    st runtime q bit hsim hlen hfaithful htag
  have hpolyInitial : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit := by simpa [hprefix] using hpoly
  have hrotEq : Equivalent n m (pre ++ first :: rest)
      (pre ++ Gate.rz θ q :: rest) := by
    simpa using Equivalent.window (n := n) (m := m) pre rest
      (equivalent_rot_rz hfirst)
  have hconstant := nonlinear_constant_phase_channel (n := n) (m := m)
    pre θ q bit hpolyInitial
  have hrest := Equivalent.append_right (n := n) (m := m) rest hconstant
  exact hrotEq.trans (by simpa [List.append_assoc] using hrest)

/-- A successful executable scan licenses its local merge even when the prefix or
intervening segment contains measurement and reset branches. -/
theorem nonlinear_merge_channel_of_relevantFaithful {n m : Nat}
    (pre gates result : List Gate) (first : Gate) (θ : ℚ) (q : Qubit)
    (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (st : NonlinearAState) (runtime : NState)
    (hprefix : st = (NonlinearAState.initial n).steps pre)
    (hsim : NonlinearSim draws evaluation st runtime)
    (hlen : st.wires.length = n)
    (hfirst : rotAngle first = some (θ, q))
    (hfaithful : PolynomialFaithful evaluation
      (nonlinearRelevant n st (first :: gates)))
    (hmerge : mergeIntoNonlinear draws runtime (runtime.tagOf q) θ gates = some result) :
    Equivalent n m (pre ++ first :: gates) (pre ++ result) := by
  obtain ⟨middle, rest, later, φ, q', sign, hgates, hresult, -, hlater, hmatch⟩ :=
    mergeIntoNonlinear_spec draws (runtime.tagOf q) θ gates result runtime hmerge
  have hpending : (st.wireOf q).polynomial ∈
      nonlinearVisited n st (first :: gates) :=
    nonlinearWire_mem_visited hlen q _
  have hlaterPoly : ((st.steps middle).wireOf q').polynomial ∈
      nonlinearVisited n st (first :: gates) := by
    refine List.mem_append.2 (Or.inr ?_)
    rw [nonlinear_step_of_rotAngle hfirst, hgates]
    exact nonlinearVisited_append_mem middle (later :: rest)
      (nonlinearWire_mem_visited (by rw [NonlinearAState.length_steps, hlen]) q' _)
  have hcomparison : (if sign then (st.wireOf q).polynomial.flip
      else (st.wireOf q).polynomial) ∈ nonlinearRelevant n st (first :: gates) := by
    cases sign with
    | false => exact List.mem_append.2 (Or.inl hpending)
    | true => exact List.mem_append.2 (Or.inr
        (List.mem_map.2 ⟨_, hpending, rfl⟩))
  have hsimMiddle := nonlinearSim_steps hsim middle
  have hmatch' : matchFingerprint
      (evaluation.eval (st.wireOf q).polynomial)
      (evaluation.eval ((st.steps middle).wireOf q').polynomial) = some sign := by
    rw [← nonlinearSim_tagOf hsim q, ← nonlinearSim_tagOf hsimMiddle q']
    exact hmatch
  have hpoly : ((st.steps middle).wireOf q').polynomial =
      if sign then (st.wireOf q).polynomial.flip else (st.wireOf q).polynomial :=
    matchFingerprint_exact hfaithful hcomparison
      (List.mem_append.2 (Or.inl hlaterPoly)) hmatch'
  have hpolyInitial :
      ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
        if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
        else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial := by
    simpa [hprefix] using hpoly
  have hfirstEq : Equivalent n m (pre ++ first :: gates)
      (pre ++ Gate.rz θ q :: gates) := by
    simpa using Equivalent.window (n := n) (m := m) pre gates
      (equivalent_rot_rz hfirst)
  have hlaterEq : Equivalent n m
      (pre ++ Gate.rz θ q :: middle ++ later :: rest)
      (pre ++ Gate.rz θ q :: middle ++ Gate.rz φ q' :: rest) := by
    simpa using Equivalent.window (n := n) (m := m)
      (pre ++ Gate.rz θ q :: middle) rest (equivalent_rot_rz hlater)
  have hhop : Equivalent n m
      (pre ++ Gate.rz θ q :: middle ++ Gate.rz φ q' :: rest)
      (pre ++ middle ++ Gate.rz (φ + signedAngle sign θ) q' :: rest) := by
    have hlocal := nonlinear_phase_hop_channel (n := n) (m := m)
      pre middle θ q q' sign hpolyInitial
    have htail := Equivalent.append_right (n := n) (m := m)
      (Gate.rz φ q' :: rest) hlocal
    have hcombine : Equivalent n m
        (pre ++ middle ++ Gate.rz (signedAngle sign θ) q' :: Gate.rz φ q' :: rest)
        (pre ++ middle ++ Gate.rz (φ + signedAngle sign θ) q' :: rest) := by
      have hrot := merge_equivalent (n := n) (m := m) [] []
        (signedAngle sign θ) φ q' q' false (by simp) (by simp)
        (by
          intro b k b' hA hB
          have hk : k = b := AState.one_entry hA
          have hb : b' = k := AState.one_entry hB
          subst b
          subst b'
          simp)
      have hwindow := Equivalent.window (n := n) (m := m)
        (pre ++ middle) rest hrot
      simpa [signedAngle] using hwindow
    have htail' : Equivalent n m
        (pre ++ Gate.rz θ q :: middle ++ Gate.rz φ q' :: rest)
        (pre ++ middle ++ Gate.rz (signedAngle sign θ) q' :: Gate.rz φ q' :: rest) := by
      simpa [List.append_assoc] using htail
    exact htail'.trans hcombine
  rw [hgates, hresult]
  rw [hgates] at hfirstEq
  have hfirstEq' : Equivalent n m
      (pre ++ first :: middle ++ later :: rest)
      (pre ++ Gate.rz θ q :: middle ++ later :: rest) := by
    simpa [List.append_assoc] using hfirstEq
  simpa [List.append_assoc] using hfirstEq'.trans (hlaterEq.trans hhop)

/-- Under collision-free packed comparisons, the executable nonlinear fold preserves
the full classical-quantum channel, including CCX and measurement. -/
theorem foldFromNonlinear_correct {n m : Nat}
    (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (targets : Array Bool) :
    ∀ (N : Nat) (gates : List Gate), gates.length ≤ N →
      ∀ (pre : List Gate) (at_ : Nat) (st : NonlinearAState) (runtime : NState),
        st = (NonlinearAState.initial n).steps pre →
        NonlinearSim draws evaluation st runtime →
        st.wires.length = n →
        PolynomialFaithful evaluation (nonlinearRelevant n st gates) →
        Equivalent n m
          (pre ++ foldFromNonlinear draws targets runtime at_ gates)
          (pre ++ gates) := by
  intro N
  induction N with
  | zero =>
      intro gates hlen pre at_ st runtime _ _ _ _
      have hnil : gates = [] := List.eq_nil_of_length_eq_zero (Nat.le_zero.1 hlen)
      subst gates
      simp [foldFromNonlinear, Equivalent.refl]
  | succ N ih =>
      intro gates hlen pre at_ st runtime hprefix hsim hstlen hfaithful
      cases gates with
      | nil => simp [foldFromNonlinear, Equivalent.refl]
      | cons gate tail =>
          have htailLen : tail.length ≤ N := by
            simp only [List.length_cons] at hlen
            omega
          have hprefixStep : st.step gate =
              (NonlinearAState.initial n).steps (pre ++ [gate]) := by
            rw [hprefix, nonlinear_steps_append]
            rfl
          have hfaithTail : PolynomialFaithful evaluation
              (nonlinearRelevant n (st.step gate) tail) :=
            nonlinearFaithful_of_sub hfaithful (by
              intro p hp
              exact List.mem_append.2 (Or.inr hp))
          have hkeep : Equivalent n m
              (pre ++ gate :: foldFromNonlinear draws targets
                (runtime.step draws gate) (at_ + 1) tail)
              (pre ++ gate :: tail) := by
            have hIH := ih tail htailLen (pre ++ [gate]) (at_ + 1)
              (st.step gate) (runtime.step draws gate) hprefixStep
              (nonlinearSim_step hsim gate)
              (by rw [NonlinearAState.length_step, hstlen]) hfaithTail
            simpa [List.append_assoc] using hIH
          simp only [foldFromNonlinear]
          split
          · split
            · rename_i θ q hrot hconstant
              have hIH : Equivalent n m
                  (pre ++ foldFromNonlinear draws targets
                    (runtime.step draws gate) (at_ + 1) tail)
                  (pre ++ tail) := by
                have hprefixSame : st.step gate =
                    (NonlinearAState.initial n).steps pre := by
                  rw [nonlinear_step_of_rotAngle hrot, hprefix]
                exact ih tail htailLen pre (at_ + 1) (st.step gate)
                  (runtime.step draws gate) hprefixSame
                  (nonlinearSim_step hsim gate)
                  (by rw [NonlinearAState.length_step, hstlen]) hfaithTail
              have htag : runtime.tagOf q = 0 ∨ runtime.tagOf q = 1 := by
                simpa using hconstant
              have hdrop : Equivalent n m (pre ++ gate :: tail) (pre ++ tail) := by
                rcases htag with hzero | hone
                · exact nonlinear_drop_constant_channel_of_tag pre tail gate θ q false
                    draws evaluation st runtime hprefix hsim hstlen hrot hfaithful
                    (by simpa using hzero)
                · exact nonlinear_drop_constant_channel_of_tag pre tail gate θ q true
                    draws evaluation st runtime hprefix hsim hstlen hrot hfaithful
                    (by simpa using hone)
              exact hIH.trans hdrop.symm
            · split
              · split
                · rename_i θ q hrot _ _ gs' hmerge
                  have hmergedLen : gs'.length ≤ N := by
                    rw [mergeIntoNonlinear_length draws (runtime.tagOf q) θ tail
                      gs' runtime hmerge]
                    exact htailLen
                  obtain ⟨middle, rest, later, φ, q', sign, htailEq,
                    hmergedEq, -, hlater, -⟩ :=
                    mergeIntoNonlinear_spec draws (runtime.tagOf q) θ tail
                      gs' runtime hmerge
                  have hfaithMerged : PolynomialFaithful evaluation
                      (nonlinearRelevant n st gs') :=
                    nonlinearFaithful_of_sub hfaithful (by
                      rw [hmergedEq, htailEq]
                      exact nonlinearVisited_merged_subset (n := n) st gate later middle rest
                        θ φ q q' sign hrot hlater)
                  have hIH := ih gs' hmergedLen pre (at_ + 1) st runtime
                    hprefix hsim hstlen hfaithMerged
                  have hlocal := nonlinear_merge_channel_of_relevantFaithful
                    (n := n) (m := m) pre tail gs' gate θ q draws evaluation
                    st runtime hprefix hsim hstlen hrot hfaithful hmerge
                  exact hIH.trans hlocal.symm
                · exact hkeep
              · exact hkeep
          · exact hkeep

/-- The exact executable nonlinear gate transformation is channel-equivalent to its
input whenever its packed polynomial evaluations are collision-free on the visited set. -/
theorem phaseFoldGatesNonlinear_correct {n m : Nat}
    (draws : Nat → Tag) (evaluation : PackedEvaluation draws)
    (gates : List Gate)
    (hfaithful : PolynomialFaithful evaluation
      (nonlinearRelevant n (NonlinearAState.initial n) gates)) :
    Equivalent n m (phaseFoldGatesNonlinear draws n gates) gates := by
  apply Equivalent.trans (emitAll_correct _)
  exact foldFromNonlinear_correct draws evaluation
    (nonlinearMergeTargets draws n gates)
    gates.length gates le_rfl [] 0 (NonlinearAState.initial n)
    (NState.initial draws n) rfl (nonlinearSim_initial draws evaluation n)
    (by simp [NonlinearAState.initial]) hfaithful

/-- The sampled executable circuit is correct whenever its concrete packed evaluation
exists and is collision-free on this circuit's comparison set. -/
theorem phaseFoldNonlinearWithSample_correct {n m k : Nat}
    (circuit : Circuit n m) (sample : Sample (varBound circuit.raw) k)
    (evaluation : PackedEvaluation (wordsOf k (liftSample sample)))
    (hfaithful : PolynomialFaithful evaluation
      (nonlinearRelevantCircuit circuit.raw)) :
    Circuit.Equivalent (phaseFoldNonlinearWithSample k circuit sample) circuit := by
  rcases circuit with ⟨raw, hn, hm, hwf⟩
  subst n
  subst m
  simpa [Circuit.Equivalent, phaseFoldNonlinearWithSample,
    phaseFoldNonlinear, nonlinearRelevantCircuit] using
    (phaseFoldGatesNonlinear_correct (n := raw.numQubits)
      (m := raw.numCbits) (wordsOf k (liftSample sample)) evaluation
      raw.gates hfaithful)

end

end TzapLean
