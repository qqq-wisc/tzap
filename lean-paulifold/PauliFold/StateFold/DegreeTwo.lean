import PauliFold.StateFold.Strengthen

/-!
# A fold that StateFold finds and PauliFold misses

The 21-gate, 11-T circuit `t q0; G; tdg q0; G; t q0` on three qubits, where

  `G = h q0; tdg q0; cx q2,q0; t q0; cx q1,q0; tdg q0; cx q2,q0; t q0; h q0`

is a relative-phase Toffoli: target `q0`, controls `q1` and `q2`. Its four Ts put the phases
`-[y] + [y ⊕ x₂] - [y ⊕ x₁ ⊕ x₂] + [y ⊕ x₁]` on the path variable `y` of its first H. These
add up to `π · y x₁ x₂`, plus a term free of `y`. The middle T† acts on the
Toffoli-computed value `x₀ ⊕ x₁x₂`, and the second `G` uncomputes it.

**StateFold at degree 2 folds the outer Ts** (`midE_proves`). It takes two HH steps:

1. sum out `y₃` and substitute `y₄ := x₀ ⊕ x₁x₂` (degree 2);
2. sum out `y₅` and substitute `y₆ := x₀`.

After them, the last T's predicate `y₆` has become `x₀`, the first T's.

**PauliFold** reaches `⊤` at the first inner T† (`midE_pauli`): after the H, the axis is
`X₀`. So on this segment degree-2 StateFold is strictly more precise (`midE_ssubset`), for
every `d ≥ 2`, and so is Strengthen (`midE_ssubsetS`).

The tools behave the same way:

| tool | T-count, 11 in |
|---|---|
| Feynman `-statefold 2` or `-statefold 0` | 9 (the outer Ts merge into an S) |
| Feynman `-statefold 1` | 11 |
| tzap `--passes PhaseFoldPauli,CancelGates --fixpoint` | 11 |
| tzap `-O3` | 11 |

**Minimality.** The circuit was found by a random search over Toffoli circuits, followed by
local search: deleting one or two gates, and replacing gates. No smaller circuit in that
neighborhood keeps the separation. It is also minimal for its construction:

- A block's Ts must produce `π · y x₁ x₂`. Writing `[x₁ ∧ x₂]` over parities needs all four
  characters `1, x₁, x₂, x₁ ⊕ x₂`, so each block needs at least four Ts.
- One block leaves the wire nonlinear, so a second block is needed to bring it back.
- The middle rotation must be a T: without it the pass pairs the two blocks' Ts and clears
  them, and it sees through any Clifford.

That gives 2 · 4 + 2 + 1 = 11 Ts.
-/

namespace PauliFold.SF

open TzapLean Matrix

noncomputable section

/-- `e^{iπa} = e^{iπb}` when `a - b` is an even integer. -/
theorem ep_eq_of_sub_int {a b : ℚ} (k : ℤ) (h : a = b + 2 * k) : ep a = ep b := by
  rw [h, ep_add]
  have : ep (2 * k) = 1 := by
    unfold ep
    have e : ((Real.pi * ((2 * k : ℚ) : ℝ) : ℝ) : ℂ) * Complex.I = k * (2 * Real.pi * Complex.I) := by
      push_cast; ring
    rw [e, Complex.exp_int_mul_two_pi_mul_I]
  rw [this, mul_one]

macro "ep_num" : tactic => `(tactic| first
  | rfl
  | (refine ep_eq_of_sub_int (0) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (1) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (-1) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (2) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (-2) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (3) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (-3) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (4) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (-4) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (5) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (-5) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (6) ?_; norm_num; done)
  | (refine ep_eq_of_sub_int (-6) ?_; norm_num; done))

/-- The relative-phase Toffoli block. -/
def blockG : List (CTGate 3) :=
  [.h 0, .tdg 0, .cx 2 0 (by decide), .t 0, .cx 1 0 (by decide), .tdg 0, .cx 2 0 (by decide),
   .t 0, .h 0]

/-- The segment between the two outer Ts. -/
def midE : List (CTGate 3) := blockG ++ [.tdg 0] ++ blockG

theorem midE_temps (y : ℕ) : y ∈ (run midE).temps ↔ 3 ≤ y ∧ y ≤ 6 := by
  simp [run, midE, blockG, runFrom, step, rot, init]
  omega

theorem midE_ket (ν : Val) : ev ν ((run midE).ket 0) = ν 6 ∧
    ev ν ((run midE).ket 1) = ν 1 ∧ ev ν ((run midE).ket 2) = ν 2 := by
  simp [run, midE, blockG, runFrom, step, rot, init]

theorem midE_pauli : PauliAbs.run midE (.one (zString {0} false)) = .top := by
  simp [PauliAbs.run, midE, blockG, PauliAbs.step, SPauli.conjClifford, SPauli.apply1,
    SPauli.apply2, zString, Letter.conjH, Letter.conjCX, Letter.anticommutesZ]

/-! ## The two HH steps -/

abbrev A0 : State 3 := run midE

/-- `x₀ ⊕ x₁x₂`, degree 2: the Toffoli's output. -/
abbrev R1 : Poly := BoolPolynomial.var 0 + BoolPolynomial.var 1 * BoolPolynomial.var 2
/-- `x₀`. -/
abbrev R2 : Poly := BoolPolynomial.var 0

abbrev Φ1 : Val → ℚ := fun ν => A0.phase (Function.update ν 3 false)
abbrev A1 : State 3 := hhResult A0 3 4 R1 Φ1
abbrev Φ2 : Val → ℚ := fun ν => A1.phase (Function.update ν 5 false)
abbrev A2 : State 3 := hhResult A1 5 6 R2 Φ2

set_option maxHeartbeats 4000000 in
theorem split1 (ν : Val) :
    ep (A0.phase ν) = ep (Φ1 ν + b2q (ν 3 && (ν 4 != ev ν R1))) := by
  simp [Φ1, A0, run, midE, blockG, runFrom, step, rot, init, ev_add, ev_mul, ev_var, Function.update_apply]
  generalize ν 0 = a0; generalize ν 1 = a1; generalize ν 2 = a2
  generalize ν 3 = a3; generalize ν 4 = a4; generalize ν 5 = a5; generalize ν 6 = a6
  cases a0 <;> cases a1 <;> cases a2 <;> cases a3 <;> cases a4 <;> cases a5 <;> cases a6
  all_goals simp only [b2q, Bool.bne_true, Bool.bne_false, Bool.not_true,
      Bool.not_false, Bool.and_true, Bool.and_false, Bool.true_and, Bool.false_and, if_true,
      if_false, Bool.false_eq_true, Bool.true_eq_false]
  all_goals ep_num

set_option maxHeartbeats 4000000 in
theorem split2 (ν : Val) :
    ep (A1.phase ν) = ep (Φ2 ν + b2q (ν 5 && (ν 6 != ev ν R2))) := by
  simp [Φ2, A1, Φ1, A0, hhResult, run, midE, blockG, runFrom, step, rot, init, ev_add, ev_mul, ev_var, Function.update_apply]
  generalize ν 0 = a0; generalize ν 1 = a1; generalize ν 2 = a2
  generalize ν 5 = a5; generalize ν 6 = a6
  cases a0 <;> cases a1 <;> cases a2 <;> cases a5 <;> cases a6
  all_goals simp only [b2q, Bool.bne_true, Bool.bne_false, Bool.not_true,
      Bool.not_false, Bool.and_true, Bool.and_false, Bool.true_and, Bool.false_and, if_true,
      if_false, Bool.false_eq_true, Bool.true_eq_false]
  all_goals ep_num

theorem supp_phase_update {n : ℕ} {A : State n} (hA : A.WF) (y : ℕ) :
    Supp (fun ν => A.phase (Function.update ν y false)) (A.dom \ {y}) := fun ν ν' h =>
  hA.phase _ _ fun v hv => by
    by_cases hvy : v = y
    · subst hvy; simp
    · rw [Function.update_of_ne hvy, Function.update_of_ne hvy]; exact h v ⟨hv, hvy⟩

theorem mem_dom_input {A : State 3} {v : ℕ} (hv : v < 3) (S : Set ℕ) (hS : v ∉ S) :
    v ∈ A.dom \ S := ⟨Or.inl hv, hS⟩

theorem A1_temps (y : ℕ) : y ∈ A1.temps ↔ 5 ≤ y ∧ y ≤ 6 := by
  show y ∈ (A0.temps.erase 3).erase 4 ↔ _
  simp only [Finset.mem_erase, midE_temps]; omega

theorem A1_ket (ν : Val) : ev ν (A1.ket 0) = ν 6 ∧ ev ν (A1.ket 1) = ν 1 ∧
    ev ν (A1.ket 2) = ν 2 := by
  obtain ⟨h0, h1, h2⟩ := midE_ket (Function.update ν 4 (ev ν R1))
  refine ⟨?_, ?_, ?_⟩
  · show ev ν (subst 4 R1 (A0.ket 0)) = _; rw [ev_subst, h0]; simp [Function.update_apply]
  · show ev ν (subst 4 R1 (A0.ket 1)) = _; rw [ev_subst, h1]; simp [Function.update_apply]
  · show ev ν (subst 4 R1 (A0.ket 2)) = _; rw [ev_subst, h2]; simp [Function.update_apply]

theorem hh1 : HH 2 A0 A1 3 4 R1 Φ1 where
  hy := (midE_temps 3).2 (by norm_num)
  hz := (midE_temps 4).2 (by norm_num)
  hyz := by decide
  ket q ν ν' h := by
    have e6 : ν 6 = ν' 6 := h 6 ⟨Or.inr ((midE_temps 6).2 (by norm_num)), by simp⟩
    have e1 : ν 1 = ν' 1 := h 1 (mem_dom_input (by norm_num) _ (by simp))
    have e2 : ν 2 = ν' 2 := h 2 (mem_dom_input (by norm_num) _ (by simp))
    fin_cases q
    · show ev ν (A0.ket 0) = ev ν' (A0.ket 0); rw [(midE_ket ν).1, (midE_ket ν').1, e6]
    · show ev ν (A0.ket 1) = ev ν' (A0.ket 1); rw [(midE_ket ν).2.1, (midE_ket ν').2.1, e1]
    · show ev ν (A0.ket 2) = ev ν' (A0.ket 2); rw [(midE_ket ν).2.2, (midE_ket ν').2.2, e2]
  deg := (MvPolynomial.totalDegree_add _ _).trans (max_le
    (by simp [BoolPolynomial.var, MvPolynomial.totalDegree_X])
    ((MvPolynomial.totalDegree_mul _ _).trans (by simp [BoolPolynomial.var, MvPolynomial.totalDegree_X])))
  suppR ν ν' h := by
    show ev ν R1 = ev ν' R1
    simp only [R1, ev_add, ev_mul, ev_var]
    rw [h 0 (mem_dom_input (by norm_num) _ (by simp)), h 1 (mem_dom_input (by norm_num) _ (by simp)),
      h 2 (mem_dom_input (by norm_num) _ (by simp))]
  supp₀ := supp_phase_update (wf_run midE) 3
  split := split1
  result := rfl

theorem wf_A1 : A1.WF := wf_hh hh1 (wf_run midE)

theorem hh2 : HH 2 A1 A2 5 6 R2 Φ2 where
  hy := (A1_temps 5).2 (by norm_num)
  hz := (A1_temps 6).2 (by norm_num)
  hyz := by decide
  ket q ν ν' h := by
    have e6 : ν 6 = ν' 6 := h 6 ⟨Or.inr ((A1_temps 6).2 (by norm_num)), by simp⟩
    have e1 : ν 1 = ν' 1 := h 1 (mem_dom_input (by norm_num) _ (by simp))
    have e2 : ν 2 = ν' 2 := h 2 (mem_dom_input (by norm_num) _ (by simp))
    fin_cases q
    · show ev ν (A1.ket 0) = ev ν' (A1.ket 0); rw [(A1_ket ν).1, (A1_ket ν').1, e6]
    · show ev ν (A1.ket 1) = ev ν' (A1.ket 1); rw [(A1_ket ν).2.1, (A1_ket ν').2.1, e1]
    · show ev ν (A1.ket 2) = ev ν' (A1.ket 2); rw [(A1_ket ν).2.2, (A1_ket ν').2.2, e2]
  deg := by simp [R2, BoolPolynomial.var, MvPolynomial.totalDegree_X]
  suppR ν ν' h := by
    simp only [R2, ev_var]; exact h 0 (mem_dom_input (by norm_num) _ (by simp))
  supp₀ := supp_phase_update wf_A1 5
  split := split2
  result := rfl

/-- Degree-2 StateFold proves that the axis `Z₀` of the first T returns to `Z₀` at the last. -/
theorem midE_proves : Proves 2 (run midE) (BoolPolynomial.var 0) 0 false := by
  refine ⟨_, _, _, .hh hh1 ?_ ?_ (.hh hh2 ?_ ?_ (.refl _ _ _)), fun ν => ?_⟩
  · intro ν ν' h; simp only [ev_var]; exact h 0 (mem_dom_input (by norm_num) _ (by simp))
  · intro ν ν' h
    show ev ν (A0.ket 0) = ev ν' (A0.ket 0)
    rw [(midE_ket ν).1, (midE_ket ν').1]
    exact h 6 ⟨Or.inr ((midE_temps 6).2 (by norm_num)), by simp⟩
  · intro ν ν' h; simp only [ev_subst, ev_var, Function.update_apply]
    exact h 0 (mem_dom_input (by norm_num) _ (by simp))
  · intro ν ν' h
    show ev ν (A1.ket 0) = ev ν' (A1.ket 0)
    rw [(A1_ket ν).1, (A1_ket ν').1]
    exact h 6 ⟨Or.inr ((A1_temps 6).2 (by norm_num)), by simp⟩
  · show ev ν (subst 6 R2 (A1.ket 0)) =
      (false != ev ν (subst 6 R2 (subst 4 R1 (BoolPolynomial.var 0))))
    rw [ev_subst, (A1_ket _).1]
    simp [ev_subst, Function.update_apply]

/-- **A segment where StateFold beats PauliFold at degree 2.** -/
theorem midE_ssubset (d : ℕ) (hd : 2 ≤ d) : sfCandidate d midE 0 ⊂ pauliCandidate midE 0 := by
  have hP : pauliCandidate midE 0 = Set.univ := by
    rw [pauliCandidate, midE_pauli]; rfl
  rw [hP]
  refine Set.ssubset_univ_iff.2 fun h => ?_
  have h0 : (0 : Op 3) ∈ sfCandidate 2 midE 0 :=
    sfCandidate_anti hd midE 0 (h ▸ Set.mem_univ _)
  exact zString_toMatrix_ne_zero _ _ (h0 0 false midE_proves).symm

/-- The same holds with Strengthen. -/
theorem midE_ssubsetS (d : ℕ) (hd : 2 ≤ d) : sfCandidateS d midE 0 ⊂ pauliCandidate midE 0 :=
  lt_of_le_of_lt (sfCandidateS_subset d midE 0) (midE_ssubset d hd)

end

end PauliFold.SF
