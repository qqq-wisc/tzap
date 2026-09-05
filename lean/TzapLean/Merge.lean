import TzapLean.Analysis
import TzapLean.Cancel

/-!
# The Merge Lemma

Phase folding's one rewrite: a rotation may be deleted from where it stands and its angle
added to a later rotation, provided the two sit on the same parity. This file proves that
rewrite correct.

The argument is entrywise. Write the prefix as a matrix `A`, the segment between the two
rotations as `B`. A rotation is a diagonal matrix whose entry at basis state `b` depends
only on one wire of `b`, so

```
    (B · D_θ(q) · A) b' b = ∑ₖ B b' k · phase_θ(k q) · A k b
    (D_θ(q') · B · A) b' b = phase_θ(b' q') · ∑ₖ B b' k · A k b
```

agree term by term as soon as `k q = b' q'` on every path with nonzero amplitude — which is
exactly what `AState.analyze_sound` delivers when the analysis assigns the two sites the same
parity. A parity matching the *complement* works the same way with the angle negated.
-/

namespace TzapLean

/-! ## Rotations as phase matrices -/

/-- A wire outside the register reads as `false`. -/
theorem basis_get_of_ge {n : Nat} (b : Basis n) {q : Qubit} (hq : ¬ q < n) : b.get q = false := by
  simp [Basis.get, hq]


/-- The phase an `rz θ` on wire `q` applies to a basis state. -/
noncomputable def rzPhase (n : Nat) (θ : ℚ) (q : Qubit) : Basis n → ℂ :=
  fun b => if b.get q then ep (θ / 2) else ep (-θ / 2)

theorem gateUnitary_rz_eq (n : Nat) (θ : ℚ) (q : Qubit) (hq : q < n) :
    gateUnitary n (Gate.rz θ q) = phaseMatrix (rzPhase n θ q) :=
  embed1_diag2_eq_phaseMatrix _ _ q hq

theorem rzPhase_add (n : Nat) (θ φ : ℚ) (q : Qubit) (b : Basis n) :
    rzPhase n θ q b * rzPhase n φ q b = rzPhase n (θ + φ) q b := by
  simp only [rzPhase]
  by_cases h : b.get q
  · rw [if_pos h, if_pos h, if_pos h, ← ep_add]
    ring_nf
  · rw [if_neg h, if_neg h, if_neg h, ← ep_add]
    ring_nf

/-- Two rotations on one wire add. -/
theorem gateUnitary_rz_mul (n : Nat) (θ φ : ℚ) (q : Qubit) (hq : q < n) :
    gateUnitary n (Gate.rz θ q) * gateUnitary n (Gate.rz φ q)
      = gateUnitary n (Gate.rz (θ + φ) q) := by
  rw [gateUnitary_rz_eq n θ q hq, gateUnitary_rz_eq n φ q hq, gateUnitary_rz_eq n _ q hq,
    phaseMatrix_mul]
  congr 1
  funext b
  exact rzPhase_add n θ φ q b

/-- A phase matrix scales rows. -/
theorem phaseMatrix_mul_apply {n : Nat} (f : Basis n → ℂ) (A : Density n) (b' b : Basis n) :
    (phaseMatrix f * A) b' b = f b' * A b' b := by
  rw [Matrix.mul_apply, Finset.sum_eq_single b']
  · simp [phaseMatrix]
  · intro k _ hk
    simp [phaseMatrix, Ne.symm hk]
  · simp

/-! ## The hop -/

/-- The angle a rotation carries after hopping onto a wire holding the complementary value. -/
def signedAngle (sign : Bool) (θ : ℚ) : ℚ := if sign then -θ else θ

@[simp] theorem signedAngle_false (θ : ℚ) : signedAngle false θ = θ := rfl
@[simp] theorem signedAngle_true (θ : ℚ) : signedAngle true θ = -θ := rfl

/-- A rotation is its phase matrix, up to a unit factor — on a wire outside the register it
is the identity, which is that phase up to the global factor the register cannot see. -/
theorem unitary_rot_smul (n : Nat) (θ : ℚ) (q : Qubit) :
    ∃ c : ℂ, c * star c = 1 ∧
      gateUnitary n (Gate.rz θ q) = c • phaseMatrix (rzPhase n θ q) := by
  by_cases hq : q < n
  · exact ⟨1, by simp, by rw [gateUnitary_rz_eq n θ q hq, one_smul]⟩
  · refine ⟨ep (θ / 2), ep_mul_star _, ?_⟩
    have hconst : phaseMatrix (rzPhase n θ q) = ep (-θ / 2) • (1 : Density n) := by
      funext out inp
      by_cases h : out = inp <;>
        simp [phaseMatrix, h, rzPhase, Basis.get, hq, Matrix.one_apply]
    rw [hconst, smul_smul, gateUnitary, embed1_eq_one _ hq, ← ep_add,
      show θ / 2 + -θ / 2 = 0 by ring, ep_zero, one_smul]

/-- **The hop.** A rotation on wire `q` before `B` is a rotation on wire `q'` after `B`, as
soon as every path with nonzero amplitude carries the value of `q` to the value of `q'`
(`sign = true`: to its complement, which negates the angle). -/
theorem phase_hop {n : Nat} (A B : Density n) (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (h : ∀ b k b' : Basis n, A k b ≠ 0 → B b' k ≠ 0 → k.get q = (b'.get q' != sign)) :
    B * (phaseMatrix (rzPhase n θ q) * A)
      = phaseMatrix (rzPhase n (signedAngle sign θ) q') * (B * A) := by
  funext b' b
  rw [Matrix.mul_apply, phaseMatrix_mul_apply, Matrix.mul_apply, Finset.mul_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [phaseMatrix_mul_apply]
  by_cases hB : B b' k = 0
  · rw [hB]; ring
  by_cases hA : A k b = 0
  · rw [hA]; ring
  have hkq : k.get q = (b'.get q' != sign) := h b k b' hA hB
  have hphase : rzPhase n θ q k = rzPhase n (signedAngle sign θ) q' b' := by
    simp only [rzPhase, hkq]
    cases sign <;> cases hb : b'.get q' <;> simp [signedAngle, hb]
  rw [hphase]; ring

/-! ## From analysis to paths -/

theorem length_initial_par (n : Nat) : (AState.initial n).par.length = n := by
  simp [AState.initial]

/-- The parity condition the analysis checks implies the path condition the hop needs: from a
generic state, whatever value the first wire holds, the second holds it (or its complement)
after the segment. -/
theorem path_of_generic {n : Nat} {M : List Gate} (hM : ∀ g ∈ M, g.isUnitary = true)
    {st : AState} (hst : st.Bounded) (hlen : st.par.length = n) (hgen : AState.Generic n st)
    {q q' : Qubit} {sign : Bool}
    (hpar : (st.steps M).parOf q' =
      if sign then (st.parOf q).flip else st.parOf q) :
    ∀ b k b' : Basis n, unitary n [] k b ≠ 0 → unitary n M b' k ≠ 0 →
      k.get q = (b'.get q' != sign) := by
  intro b k b' _hA hB
  obtain ⟨v, hv⟩ := hgen k
  obtain ⟨v', hv', hc'⟩ := AState.analyze_sound hM hst hlen hv hB
  have hkq : k.get q = (st.parOf q).evalB v := hv q
  have hagree : (st.parOf q).evalB v' = (st.parOf q).evalB v := Form.evalB_congr (hst q) hv'
  have hb' := hc' q'
  rw [hpar] at hb'
  cases sign with
  | false =>
      simp only [if_neg (by simp : ¬ (false = true))] at hb'
      rw [hkq, ← hagree, ← hb']
      simp
  | true =>
      simp only [if_pos rfl, Form.evalB_flip] at hb'
      rw [hkq, ← hagree, hb']
      simp

/-! ## The merge -/

/-- **Merging two rotations.** With `P` the prefix, `M` the segment between the two rotation
sites, both measurement-free, and the analysis assigning the second site the first site's
parity (or its complement, when `sign`), the first rotation may be deleted and its angle
added to the second. -/
theorem merge_equivalent {n m : Nat} (P M : List Gate) (θ φ : ℚ) (q q' : Qubit) (sign : Bool)
    (hP : ∀ g ∈ P, g.isUnitary = true) (hM : ∀ g ∈ M, g.isUnitary = true)
    (hpath : ∀ b k b' : Basis n, unitary n P k b ≠ 0 → unitary n M b' k ≠ 0 →
      k.get q = (b'.get q' != sign)) :
    Equivalent n m (P ++ Gate.rz θ q :: (M ++ [Gate.rz φ q']))
      (P ++ (M ++ [Gate.rz (φ + signedAngle sign θ) q'])) := by
  have hrz : ∀ (a : ℚ) (r : Qubit), (Gate.rz a r).isUnitary = true := by
    intro a r; rfl
  obtain ⟨c₁, hc₁, e₁⟩ := unitary_rot_smul n θ q
  obtain ⟨c₂, hc₂, e₂⟩ := unitary_rot_smul n φ q'
  obtain ⟨c₃, hc₃, e₃⟩ := unitary_rot_smul n (φ + signedAngle sign θ) q'
  have hunit : (c₁ * c₂ * star c₃) * star (c₁ * c₂ * star c₃) = 1 := by
    have : star (c₁ * c₂ * star c₃) = star c₁ * star c₂ * c₃ := by
      simp [star_mul, mul_comm]
    rw [this]
    calc c₁ * c₂ * star c₃ * (star c₁ * star c₂ * c₃)
        = (c₁ * star c₁) * (c₂ * star c₂) * (c₃ * star c₃) := by ring
      _ = 1 := by rw [hc₁, hc₂, hc₃]; ring
  refine equivalent_of_unitary_smul ?_ ?_ (c₁ * c₂ * star c₃) hunit ?_
  · intro g hg
    rcases List.mem_append.1 hg with h | h
    · exact hP g h
    · rcases List.mem_cons.1 h with h | h
      · subst h; exact hrz _ _
      · rcases List.mem_append.1 h with h | h
        · exact hM g h
        · rw [List.mem_singleton.1 h]; exact hrz _ _
  · intro g hg
    rcases List.mem_append.1 hg with h | h
    · exact hP g h
    · rcases List.mem_append.1 h with h | h
      · exact hM g h
      · rw [List.mem_singleton.1 h]; exact hrz _ _
  · have hadd : (fun b : Basis n => rzPhase n φ q' b * rzPhase n (signedAngle sign θ) q' b)
        = rzPhase n (φ + signedAngle sign θ) q' := by
      funext b
      exact rzPhase_add n φ (signedAngle sign θ) q' b
    simp only [unitary_append, unitary_cons, unitary_nil, one_mul, Matrix.mul_assoc, e₁, e₂, e₃,
      Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    rw [phase_hop (unitary n P) (unitary n M) θ q q' sign hpath, ← Matrix.mul_assoc,
      phaseMatrix_mul, hadd]
    congr 1
    calc c₁ * c₂ = c₁ * c₂ * (star c₃ * c₃) := by rw [mul_comm (star c₃) c₃, hc₃]; ring
      _ = c₁ * c₂ * star c₃ * c₃ := by ring

end TzapLean
