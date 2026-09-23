import TzapLean.GF128Field
import TzapLean.NonlinearFoldProof
import TzapLean.RandPass
import Mathlib.Data.Nat.Bitwise

/-!
# Packed words and the abstract fingerprint field

The coefficient at `X^i` is bit `i` of the packed word. This is the representation
boundary between the executable `Nat` arithmetic and `AdjoinRoot modulus`.
-/

namespace TzapLean.GF128Proof

open Polynomial

set_option maxRecDepth 10000
set_option maxHeartbeats 1000000

/-- The low 128 bits of a packed word, read as coefficients over `F₂`. -/
noncomputable def packedPolynomial (x : Nat) : F₂[X] :=
  ∑ i : Fin 128, if x.testBit i then X ^ (i : Nat) else 0

theorem coeff_packedPolynomial (x n : Nat) :
    (packedPolynomial x).coeff n =
      if n < 128 then if x.testBit n then 1 else 0 else 0 := by
  classical
  change (lcoeff F₂ n) (∑ i : Fin 128, if x.testBit i then X ^ (i : Nat) else 0) = _
  rw [map_sum]
  simp only [lcoeff_apply]
  have hterm (i : Fin 128) :
      (if x.testBit i then X ^ (i : Nat) else (0 : F₂[X])).coeff n =
        if x.testBit i then if (i : Nat) = n then 1 else 0 else 0 := by
    by_cases hb : x.testBit i
    · by_cases he : (i : Nat) = n
      · subst n
        simp [hb]
      · have hne : n ≠ (i : Nat) := Ne.symm he
        simp [hb, he, hne]
    · simp [hb]
  simp only [hterm]
  by_cases hn : n < 128
  · simp only [hn]
    have hfin (i : Fin 128) : ((i : Nat) = n) ↔ i = ⟨n, hn⟩ := by
      simp [Fin.ext_iff]
    simp_rw [hfin]
    have hswap (i : Fin 128) :
        (if x.testBit i then if i = ⟨n, hn⟩ then (1 : F₂) else 0 else 0) =
          if i = ⟨n, hn⟩ then if x.testBit i then 1 else 0 else 0 := by
      split_ifs <;> rfl
    simp_rw [hswap]
    simp
  · simp only [hn]
    apply Finset.sum_eq_zero
    intro i _
    have hi : (i : Nat) ≠ n := by omega
    simp [hi]

theorem packedPolynomial_normalize (x : Nat) :
    packedPolynomial (GF128.normalize x) = packedPolynomial x := by
  apply Polynomial.ext
  intro n
  simp only [coeff_packedPolynomial]
  by_cases hn : n < 128
  · simp only [hn, ↓reduceIte, GF128.normalize, GF128.mask,
      Nat.testBit_land, Nat.testBit_two_pow_sub_one]
    simp [GF128.width, hn]
  · simp [hn]

theorem packedPolynomial_xor (x y : Nat) :
    packedPolynomial (x ^^^ y) = packedPolynomial x + packedPolynomial y := by
  apply Polynomial.ext
  intro n
  simp only [coeff_packedPolynomial, coeff_add]
  by_cases hn : n < 128
  · simp only [hn, ↓reduceIte, Nat.testBit_xor]
    cases x.testBit n <;> cases y.testBit n <;> decide
  · simp [hn]

@[simp] theorem packedPolynomial_zero : packedPolynomial 0 = 0 := by
  apply Polynomial.ext
  intro n
  simp [coeff_packedPolynomial]

@[simp] theorem packedPolynomial_one : packedPolynomial 1 = 1 := by
  apply Polynomial.ext
  intro n
  by_cases hn : n = 0
  · subst n
    simp [coeff_packedPolynomial]
  · have hb : (1 : Nat).testBit n = false := by
      cases h : (1 : Nat).testBit n with
      | false => rfl
      | true => exact (hn (Nat.testBit_one_eq_true_iff_self_eq_zero.mp h)).elim
    simp [coeff_packedPolynomial, hb, coeff_one, hn]

/-- Interpret the packed coefficients in the irreducible quotient field. -/
noncomputable def packedToField (x : Nat) : Field128 :=
  AdjoinRoot.mk modulus (packedPolynomial x)

@[simp] theorem packedToField_normalize (x : Nat) :
    packedToField (GF128.normalize x) = packedToField x := by
  simp [packedToField, packedPolynomial_normalize]

@[simp] theorem packedToField_zero : packedToField 0 = 0 := by
  simp [packedToField]

@[simp] theorem packedToField_one : packedToField 1 = 1 := by
  simp [packedToField]

theorem packedToField_xor (x y : Nat) :
    packedToField (x ^^^ y) = packedToField x + packedToField y := by
  simp [packedToField, packedPolynomial_xor]

theorem packedToField_add (x y : Nat) :
    packedToField (GF128.add x y) = packedToField x + packedToField y := by
  simp [packedToField, GF128.add, packedPolynomial_xor,
    packedPolynomial_normalize]

theorem packedPolynomial_degree_lt (x : Nat) :
    (packedPolynomial x).degree < 128 := by
  apply (degree_lt_iff_coeff_zero _ _).2
  intro n hn
  simp [coeff_packedPolynomial, Nat.not_lt.mpr hn]

theorem packedPolynomial_injective_on_words {x y : Nat}
    (hx : x < 2 ^ 128) (hy : y < 2 ^ 128)
    (heq : packedPolynomial x = packedPolynomial y) : x = y := by
  apply Nat.eq_of_testBit_eq
  intro i
  by_cases hi : i < 128
  · have hc := congrArg (fun p : F₂[X] => p.coeff i) heq
    simp only [coeff_packedPolynomial, hi, ↓reduceIte] at hc
    cases hxi : x.testBit i <;> cases hyi : y.testBit i <;> simp_all
  · have hxi : x.testBit i = false :=
      Nat.testBit_eq_false_of_lt
        (lt_of_lt_of_le hx (Nat.pow_le_pow_right (by decide) (by omega)))
    have hyi : y.testBit i = false :=
      Nat.testBit_eq_false_of_lt
        (lt_of_lt_of_le hy (Nat.pow_le_pow_right (by decide) (by omega)))
    simp [hxi, hyi]

theorem packedToField_injective_on_words {x y : Nat}
    (hx : x < 2 ^ 128) (hy : y < 2 ^ 128)
    (heq : packedToField x = packedToField y) : x = y := by
  have hpoly : packedPolynomial x = packedPolynomial y := by
    have hm := congrArg (AdjoinRoot.modByMonicHom modulus_monic) heq
    simp only [packedToField, AdjoinRoot.modByMonicHom_mk] at hm
    have hdegree (z : Nat) : (packedPolynomial z).degree < modulus.degree := by
      rw [degree_eq_natDegree modulus_monic.ne_zero, modulus_natDegree]
      exact packedPolynomial_degree_lt z
    rw [(modByMonic_eq_self_iff modulus_monic).mpr (hdegree x),
      (modByMonic_eq_self_iff modulus_monic).mpr (hdegree y)] at hm
    exact hm
  exact packedPolynomial_injective_on_words hx hy hpoly

noncomputable def packedToFieldFin : Fin (2 ^ 128) → Field128 :=
  fun x => packedToField x

theorem packedToFieldFin_bijective : Function.Bijective packedToFieldFin := by
  apply (Fintype.bijective_iff_injective_and_card packedToFieldFin).2
  constructor
  · intro x y h
    apply Fin.ext
    exact packedToField_injective_on_words x.isLt y.isLt h
  · simp [field128_card]

/-- A left shift has one possible overflow coefficient, at degree 128. -/
theorem packedPolynomial_shiftLeft_one (x : Nat) :
    X * packedPolynomial x =
      packedPolynomial (x <<< 1) +
        (if x.testBit 127 then X ^ 128 else (0 : F₂[X])) := by
  apply Polynomial.ext
  intro n
  cases n with
  | zero =>
      by_cases hb : x.testBit 127 <;>
        simp [hb, coeff_packedPolynomial, coeff_X_pow]
  | succ m =>
      rw [coeff_X_mul, coeff_add]
      by_cases hm : m < 127
      · have hn : m + 1 < 128 := by omega
        have hm' : m < 128 := by omega
        by_cases hb : x.testBit 127 <;>
          simp [hb, coeff_packedPolynomial, hm', hn, Nat.testBit_shiftLeft,
            coeff_X_pow]; omega
      · by_cases he : m = 127
        · subst m
          by_cases hb : x.testBit 127 <;>
            simp [hb, coeff_packedPolynomial, coeff_X_pow]
        · have hm' : ¬m < 128 := by omega
          have hn : ¬m + 1 < 128 := by omega
          by_cases hb : x.testBit 127 <;>
            simp [hb, coeff_packedPolynomial, hm', hn, he, coeff_X_pow]

private theorem packedPolynomial_reduction :
    packedPolynomial GF128.reduction = X ^ 7 + X ^ 2 + X + 1 := by
  apply Polynomial.ext
  intro n
  rw [coeff_packedPolynomial]
  by_cases hn : n < 8
  · interval_cases n <;>
      norm_num [GF128.reduction, coeff_add, coeff_X_pow, coeff_X, coeff_one,
        Nat.testBit] <;> decide
  · have hb : (135 : Nat).testBit n = false :=
      Nat.testBit_eq_false_of_lt
        (lt_of_lt_of_le (by decide : 135 < 2 ^ 8)
          (Nat.pow_le_pow_right (by decide) (by omega)))
    have h0 : n ≠ 0 := by omega
    have h1 : n ≠ 1 := by omega
    have h1' : 1 ≠ n := Ne.symm h1
    have h2 : n ≠ 2 := by omega
    have h7 : n ≠ 7 := by omega
    simp [GF128.reduction, hb, coeff_add, coeff_X_pow, coeff_X, coeff_one,
      h0, h1', h2, h7]

private theorem modulus_root_relation :
    AdjoinRoot.mk modulus (X ^ 128) =
      AdjoinRoot.mk modulus (X ^ 7 + X ^ 2 + X + 1) := by
  have hmod : modulus = X ^ 128 + (X ^ 7 + X ^ 2 + X + 1) := by
    unfold modulus
    ring
  have hzero : AdjoinRoot.mk modulus modulus = 0 := AdjoinRoot.mk_self
  have hzero' : AdjoinRoot.mk modulus (X ^ 128 + (X ^ 7 + X ^ 2 + X + 1)) = 0 := by
    rw [← hmod]
    exact hzero
  rw [map_add] at hzero'
  apply sub_eq_zero.mp
  simpa [CharTwo.sub_eq_add] using hzero'

theorem packedToField_xtime (x : Nat) :
    packedToField (GF128.xtime x) = AdjoinRoot.root modulus * packedToField x := by
  let z := GF128.normalize x
  have hshift := congrArg (AdjoinRoot.mk modulus) (packedPolynomial_shiftLeft_one z)
  have hreduce : packedToField GF128.reduction =
      AdjoinRoot.mk modulus (X ^ 128) := by
    rw [packedToField, packedPolynomial_reduction, modulus_root_relation]
  have hxtime : packedToField (GF128.xtime x) =
      packedToField (z <<< 1) +
        (if z.testBit 127 then packedToField GF128.reduction else 0) := by
    change packedToField (GF128.normalize (z <<< 1) ^^^
      (if z.testBit 127 then GF128.reduction else 0)) = _
    rw [packedToField_xor, packedToField_normalize]
    by_cases hb : z.testBit 127 <;> simp [hb]
  rw [hxtime, hreduce]
  simp only [map_mul, map_add, AdjoinRoot.mk_X] at hshift
  have hcase : (AdjoinRoot.mk modulus)
      (if z.testBit 127 then X ^ 128 else (0 : F₂[X])) =
        if z.testBit 127 then (AdjoinRoot.mk modulus) (X ^ 128) else 0 := by
    split_ifs <;> simp
  rw [hcase] at hshift
  change AdjoinRoot.root modulus * packedToField z =
    packedToField (z <<< 1) +
      (if z.testBit 127 then (AdjoinRoot.mk modulus) (X ^ 128) else 0) at hshift
  rw [← hshift]
  change AdjoinRoot.root modulus * packedToField (GF128.normalize x) = _
  rw [packedToField_normalize]

/-- Peeling the low bit of a canonical word gives its polynomial recursion. -/
private theorem packedPolynomial_split (b : Nat) (hb : b < 2 ^ 128) :
    packedPolynomial b =
      (if b.testBit 0 then (1 : F₂[X]) else 0) + X * packedPolynomial (b >>> 1) := by
  apply Polynomial.ext
  intro n
  cases n with
  | zero =>
      by_cases hbit : b.testBit 0 <;>
        simp [hbit, coeff_packedPolynomial, coeff_one]
  | succ m =>
      rw [coeff_add, coeff_X_mul]
      have hshift : (b >>> 1).testBit m = b.testBit (m + 1) := by
        simp [Nat.testBit_shiftRight, Nat.add_comm]
      by_cases hm : m < 127
      · have hn : m + 1 < 128 := by omega
        by_cases hbit : b.testBit 0 <;>
          simp [hbit, coeff_packedPolynomial, hn, hshift, coeff_one] <;> omega
      · by_cases he : m = 127
        · subst m
          have htop : b.testBit 128 = false := Nat.testBit_eq_false_of_lt hb
          by_cases hbit : b.testBit 0 <;>
            simp [hbit, coeff_packedPolynomial, hshift, htop, coeff_one]
        · have hm' : ¬m < 128 := by omega
          have hn : ¬m + 1 < 128 := by omega
          by_cases hbit : b.testBit 0 <;>
            simp [hbit, coeff_packedPolynomial, hn, hm', coeff_one]

theorem packedToField_split (b : Nat) (hb : b < 2 ^ 128) :
    packedToField b =
      (if b.testBit 0 then 1 else 0) +
        AdjoinRoot.root modulus * packedToField (b >>> 1) := by
  have h := congrArg (AdjoinRoot.mk modulus) (packedPolynomial_split b hb)
  simp only [map_add, map_mul, AdjoinRoot.mk_X] at h
  by_cases hbit : b.testBit 0 <;>
    simpa [hbit, packedToField] using h

/-- The executable multiply loop preserves the corresponding field expression. -/
theorem packedToField_mulAux :
    ∀ fuel a b acc, fuel ≤ 128 → b < 2 ^ fuel →
      packedToField (GF128.mulAux fuel a b acc) =
        packedToField acc + packedToField a * packedToField b := by
  intro fuel
  induction fuel with
  | zero =>
      intro a b acc _ hb
      have hz : b = 0 := by simpa using hb
      subst b
      simp [GF128.mulAux]
  | succ fuel ih =>
      intro a b acc hf hb
      have hb128 : b < 2 ^ 128 :=
        lt_of_lt_of_le hb (Nat.pow_le_pow_right (by decide) hf)
      have htail : b >>> 1 < 2 ^ fuel := by
        have hpow : 2 ^ (fuel + 1) = 2 ^ fuel * 2 := by simp [pow_succ]
        have hdiv : b / 2 < 2 ^ fuel := by omega
        simpa [Nat.shiftRight_eq_div_pow] using hdiv
      have hsplit := packedToField_split b hb128
      simp only [GF128.mulAux]
      by_cases hbit : b.testBit 0
      · simp only [hbit, ↓reduceIte]
        rw [ih _ _ _ (by omega) htail, packedToField_xor, packedToField_xtime]
        simp [hbit] at hsplit
        rw [hsplit]
        ring
      · simp only [hbit]
        rw [ih _ _ _ (by omega) htail, packedToField_xtime]
        simp only [hbit] at hsplit
        rw [hsplit]
        simp at hsplit ⊢
        ring

theorem packedToField_mul (x y : Nat) :
    packedToField (GF128.mul x y) = packedToField x * packedToField y := by
  have hy : GF128.normalize y < 2 ^ GF128.width := by
    simp only [GF128.normalize, GF128.mask, Nat.and_two_pow_sub_one_eq_mod]
    exact Nat.mod_lt _ (by positivity)
  rw [GF128.mul, packedToField_mulAux GF128.width _ _ _ (by decide) hy]
  simp [packedToField_normalize]

theorem normalize_lt (x : Nat) : GF128.normalize x < 2 ^ 128 := by
  simp only [GF128.normalize, GF128.mask, GF128.width,
    Nat.and_two_pow_sub_one_eq_mod]
  exact Nat.mod_lt _ (by positivity)

theorem add_lt (x y : Nat) : GF128.add x y < 2 ^ 128 := by
  exact Nat.xor_lt_two_pow (normalize_lt x) (normalize_lt y)

private theorem mulAux_lt (fuel a b acc : Nat) :
    GF128.mulAux fuel a b acc < 2 ^ 128 := by
  induction fuel generalizing a b acc with
  | zero => exact normalize_lt acc
  | succ fuel ih => exact ih _ _ _

theorem mul_lt (x y : Nat) : GF128.mul x y < 2 ^ 128 :=
  mulAux_lt _ _ _ _

/-- The canonical packed words are exactly the elements of the quotient field. -/
noncomputable def packedEquiv : Fin (2 ^ 128) ≃ Field128 :=
  Equiv.ofBijective packedToFieldFin packedToFieldFin_bijective

noncomputable def fieldToPacked (v : Field128) : Nat :=
  (packedEquiv.symm v).val

theorem fieldToPacked_lt (v : Field128) : fieldToPacked v < 2 ^ 128 :=
  (packedEquiv.symm v).isLt

@[simp] theorem packedToField_fieldToPacked (v : Field128) :
    packedToField (fieldToPacked v) = v := by
  change packedEquiv (packedEquiv.symm v) = v
  exact packedEquiv.apply_symm_apply v

@[simp] theorem fieldToPacked_packedToField (x : Nat) :
    fieldToPacked (packedToField x) = GF128.normalize x := by
  apply packedToField_injective_on_words
    (fieldToPacked_lt _) (normalize_lt x)
  simp

theorem fieldToPacked_add (u v : Field128) :
    fieldToPacked (u + v) = GF128.add (fieldToPacked u) (fieldToPacked v) := by
  apply packedToField_injective_on_words (fieldToPacked_lt _) (add_lt _ _)
  rw [packedToField_fieldToPacked, packedToField_add]
  simp

theorem fieldToPacked_mul (u v : Field128) :
    fieldToPacked (u * v) = GF128.mul (fieldToPacked u) (fieldToPacked v) := by
  apply packedToField_injective_on_words (fieldToPacked_lt _) (mul_lt _ _)
  rw [packedToField_fieldToPacked, packedToField_mul]
  simp

@[simp] theorem fieldToPacked_zero : fieldToPacked 0 = 0 := by
  simpa using fieldToPacked_packedToField 0

@[simp] theorem fieldToPacked_one : fieldToPacked 1 = 1 := by
  have hnorm : GF128.normalize 1 = 1 := by native_decide
  simpa [hnorm] using fieldToPacked_packedToField 1

/-- The concrete evaluation required by the nonlinear runtime refinement. -/
noncomputable def packedEvaluation (draws : Nat → Tag) : PackedEvaluation draws where
  eval p := fieldToPacked (BoolPolynomial.evalAt (fun i => packedToField (draws i)) p)
  zero := by simp [BoolPolynomial.evalAt]
  one := by simp [BoolPolynomial.evalAt]
  var_eq i := by simp
  add p q := by simp [fieldToPacked_add]
  mul p q := by simp [fieldToPacked_mul]

/-- The executable state follows the formal polynomial state for concrete packed arithmetic. -/
theorem packed_nonlinearSim_steps (draws : Nat → Tag) (n : Nat) (gates : List Gate) :
    NonlinearSim draws (packedEvaluation draws)
      ((NonlinearAState.initial n).steps gates)
      ((NState.initial draws n).steps draws gates) :=
  nonlinearSim_steps (nonlinearSim_initial draws (packedEvaluation draws) n) gates

/-- Pad a finite sequence of canonical 128-bit words with zero draws. -/
def packedDraws (m : Nat) (words : Fin m → Fin (2 ^ 128)) : Nat → Tag :=
  fun i => if h : i < m then (words ⟨i, h⟩).val else 0

/-- The abstract field valuation induced by those very same packed words. -/
noncomputable def packedValuation (m : Nat) (words : Fin m → Fin (2 ^ 128)) :
    Fin m → Field128 := fun i => packedEquiv (words i)

theorem packedDraws_field (m : Nat) (words : Fin m → Fin (2 ^ 128)) (i : Nat) :
    packedToField (packedDraws m words i) =
      if h : i < m then packedValuation m words ⟨i, h⟩ else 0 := by
  by_cases hi : i < m
  · simp [packedDraws, packedValuation, hi, packedEquiv, packedToFieldFin]
  · simp [packedDraws, hi]

/-- Concrete polynomial evaluation agrees with field evaluation for every polynomial. -/
theorem packedEvaluation_agrees (m : Nat) (words : Fin m → Fin (2 ^ 128))
    (p : BoolPolynomial) :
    (packedEvaluation (packedDraws m words)).eval p =
      fieldToPacked (BoolPolynomial.evalAt
        (fun i => if h : i < m then packedValuation m words ⟨i, h⟩ else 0) p) := by
  simp only [packedEvaluation, packedDraws_field]

/-- A collision-free field valuation makes the executable comparisons exact. -/
theorem packedEvaluation_faithful (m : Nat) (words : Fin m → Fin (2 ^ 128))
    (polynomials : List BoolPolynomial)
    (hgood : packedValuation m words ∉
      nonlinearBadValuations (K := Field128) m polynomials) :
    PolynomialFaithful (packedEvaluation (packedDraws m words)) polynomials := by
  apply nonlinearGoodValuation_packedFaithful m polynomials
    (packedValuation m words) (packedEvaluation (packedDraws m words)) fieldToPacked
  · intro x y h
    exact (packedToField_fieldToPacked x).symm.trans
      ((congrArg packedToField h).trans (packedToField_fieldToPacked y))
  · intro p _
    exact packedEvaluation_agrees m words p
  · exact hgood

/-- Coordinatewise, canonical packed words and field elements have the same sample space. -/
noncomputable def packedValuationEquiv (m : Nat) :
    (Fin m → Fin (2 ^ 128)) ≃ (Fin m → Field128) :=
  Equiv.piCongrRight (fun _ => packedEquiv)

theorem packedValuationEquiv_apply (m : Nat) (words : Fin m → Fin (2 ^ 128)) :
    packedValuationEquiv m words = packedValuation m words := rfl

/-- The bad-seed event, stated directly on the executable's canonical word samples. -/
def packedBadValuations (m : Nat) (polynomials : List BoolPolynomial) :
    Set (Fin m → Fin (2 ^ 128)) :=
  {words | packedValuation m words ∈
    nonlinearBadValuations (K := Field128) m polynomials}

/-- A bijective encoding preserves the uniform probability of the bad-seed event. -/
theorem packedBadProbability_eq (m : Nat) (polynomials : List BoolPolynomial) :
    (PMF.uniformOfFintype (Fin m → Fin (2 ^ 128))).toOuterMeasure
      (packedBadValuations m polynomials) =
    (PMF.uniformOfFintype (Fin m → Field128)).toOuterMeasure
      (↑(nonlinearBadValuations (K := Field128) m polynomials) :
        Set (Fin m → Field128)) := by
  classical
  rw [PMF.toOuterMeasure_uniformOfFintype_apply,
    PMF.toOuterMeasure_uniformOfFintype_apply]
  have hdenom : Fintype.card (Fin m → Fin (2 ^ 128)) =
      Fintype.card (Fin m → Field128) :=
    Fintype.card_congr (packedValuationEquiv m)
  have hnumer : Fintype.card
      {words : Fin m → Fin (2 ^ 128) // words ∈ packedBadValuations m polynomials} =
      Fintype.card {valuation : Fin m → Field128 //
        valuation ∈ nonlinearBadValuations (K := Field128) m polynomials} := by
    apply Fintype.card_congr
    apply (packedValuationEquiv m).subtypeEquiv
    intro words
    simp only [packedBadValuations, Set.mem_ofPred_eq, packedValuationEquiv_apply]
  rw [hnumer, hdenom]
  simp

/-- The concrete packed-word distribution inherits the sharp field collision bound. -/
theorem packedCircuit_badProbability_bound (circuit : RawCircuit) :
    (PMF.uniformOfFintype (Fin (varBound circuit) → Fin (2 ^ 128))).toOuterMeasure
      (packedBadValuations (varBound circuit)
        (nonlinearRelevantCircuit circuit)) ≤
      ((nonlinearComparisonPairs (nonlinearRelevantCircuit circuit)).card : ENNReal) *
        ((Fingerprint.maxDegree : ENNReal) / ((2 ^ 128 : Nat) : ENNReal)) := by
  rw [packedBadProbability_eq]
  exact field128_badProbability_bound circuit

/-- Outside that event, the actual nonlinear folding transformation is channel-equivalent
to its input, including CCX products and the degree-cutoff fallback. -/
theorem packedPhaseFold_correct {n m : Nat} (circuit : Circuit n m)
    (words : Fin (varBound circuit.raw) → Fin (2 ^ 128))
    (hgood : words ∉ packedBadValuations (varBound circuit.raw)
      (nonlinearRelevantCircuit circuit.raw)) :
    Circuit.Equivalent
      ⟨phaseFoldNonlinear (packedDraws (varBound circuit.raw) words) circuit.raw,
        (phaseFoldNonlinear_numQubits _ circuit.raw).trans circuit.numQubits_eq,
        (phaseFoldNonlinear_numCbits _ circuit.raw).trans circuit.numCbits_eq,
        phaseFoldGatesNonlinear_wf _ circuit.wf⟩ circuit := by
  have hfaithful := packedEvaluation_faithful (varBound circuit.raw) words
    (nonlinearRelevantCircuit circuit.raw) hgood
  rcases circuit with ⟨raw, hn, hm, hwf⟩
  subst n
  subst m
  simpa [Circuit.Equivalent, phaseFoldNonlinear, nonlinearRelevantCircuit] using
    (phaseFoldGatesNonlinear_correct (n := raw.numQubits)
      (m := raw.numCbits) (packedDraws (varBound raw) words)
      (packedEvaluation (packedDraws (varBound raw) words)) raw.gates hfaithful)

/-- Packing a 128-bit sample produces a canonical word. -/
theorem bitsToWord_lt (bits : BitString 128) : bitsToWord bits < 2 ^ 128 := by
  apply Nat.lt_of_testBit 128
  · simp [bitsToWord, testBit_bitsToWordAux]
  · decide
  · intro j hj
    have hj' : ¬j < 128 := by omega
    have hleft : (bitsToWord bits).testBit j = false := by
      simp [bitsToWord, testBit_bitsToWordAux, hj']
    rw [hleft]
    exact (Nat.testBit_two_pow_of_ne (by omega : 128 ≠ j)).symm

/-- The bit-vector sample representation and the packed-word representation are equivalent. -/
noncomputable def bitStringWordEquiv : BitString 128 ≃ Fin (2 ^ 128) where
  toFun bits := ⟨bitsToWord bits, bitsToWord_lt bits⟩
  invFun word := wordToBits word.val
  left_inv bits := wordToBits_bitsToWord bits
  right_inv word := by
    apply Fin.ext
    apply Nat.eq_of_testBit_eq
    intro i
    by_cases hi : i < 128
    · simp only [bitsToWord, testBit_bitsToWordAux, hi, decide_true,
        Bool.true_and, wordToBits]
      exact unbit_bit _
    · have hzero : word.val.testBit i = false := by
        exact Nat.testBit_eq_false_of_lt (lt_of_lt_of_le word.isLt (pow_le_pow_right₀ (by omega) (by omega)))
      simp [bitsToWord, testBit_bitsToWordAux, hi, hzero]

/-- The actual 128-bit sample type used by `phaseFoldNonlinearWithSample` is coordinatewise
equivalent to the canonical packed-word sample type. -/
noncomputable def packedSampleEquiv (m : Nat) :
    Sample m 128 ≃ (Fin m → Fin (2 ^ 128)) :=
  Equiv.piCongrRight (fun _ => bitStringWordEquiv)

private theorem bitsToWord_zero : bitsToWord (0 : BitString 128) = 0 := by
  have h := bitStringWordEquiv.right_inv (0 : Fin (2 ^ 128))
  simpa [bitStringWordEquiv, wordToBits_zero] using congrArg Fin.val h

/-- Packing the executable sample and padding it yields exactly the draw stream modeled
by `packedDraws`. -/
theorem wordsOf_liftSample_eq_packedDraws (m : Nat) (sample : Sample m 128) :
    wordsOf 128 (liftSample sample) =
      packedDraws m (packedSampleEquiv m sample) := by
  funext i
  by_cases hi : i < m
  · simp [wordsOf, liftSample, packedDraws, packedSampleEquiv, bitStringWordEquiv, hi]
  · simp [wordsOf, liftSample, packedDraws, hi, bitsToWord_zero]

/-- The channel-correctness theorem applies to the executable's exact sample representation. -/
theorem phaseFoldNonlinearWithSample128_correct {n m : Nat} (circuit : Circuit n m)
    (sample : Sample (varBound circuit.raw) 128)
    (hgood : packedSampleEquiv (varBound circuit.raw) sample ∉
      packedBadValuations (varBound circuit.raw)
        (nonlinearRelevantCircuit circuit.raw)) :
    Circuit.Equivalent (phaseFoldNonlinearWithSample 128 circuit sample) circuit := by
  have h := packedPhaseFold_correct circuit
    (packedSampleEquiv (varBound circuit.raw) sample) hgood
  simpa [phaseFoldNonlinearWithSample, wordsOf_liftSample_eq_packedDraws] using h

/-- Bad seeds in the actual executable sample format. -/
def nonlinearBadSamples (circuit : RawCircuit) : Set (Sample (varBound circuit) 128) :=
  {sample | packedSampleEquiv (varBound circuit) sample ∈
    packedBadValuations (varBound circuit)
      (nonlinearRelevantCircuit circuit)}

/-- A uniform 128-bit bit-string sample has the same bad-event probability as uniform
canonical packed words. -/
theorem nonlinearBadSamples_probability_eq (circuit : RawCircuit) :
    (PMF.uniformOfFintype (Sample (varBound circuit) 128)).toOuterMeasure
      (nonlinearBadSamples circuit) =
    (PMF.uniformOfFintype
      (Fin (varBound circuit) → Fin (2 ^ 128))).toOuterMeasure
      (packedBadValuations (varBound circuit)
        (nonlinearRelevantCircuit circuit)) := by
  classical
  rw [PMF.toOuterMeasure_uniformOfFintype_apply,
    PMF.toOuterMeasure_uniformOfFintype_apply]
  have hdenom : Fintype.card (Sample (varBound circuit) 128) =
      Fintype.card (Fin (varBound circuit) → Fin (2 ^ 128)) :=
    Fintype.card_congr (packedSampleEquiv (varBound circuit))
  have hnumer : Fintype.card
      {sample : Sample (varBound circuit) 128 //
        sample ∈ nonlinearBadSamples circuit} =
      Fintype.card {words : Fin (varBound circuit) → Fin (2 ^ 128) //
        words ∈ packedBadValuations (varBound circuit)
          (nonlinearRelevantCircuit circuit)} := by
    apply Fintype.card_congr
    apply (packedSampleEquiv (varBound circuit)).subtypeEquiv
    intro sample
    rfl
  rw [hnumer, hdenom]

/-- The sharp nonlinear failure bound for the exact finite sample format of the executable. -/
theorem phaseFoldNonlinearWithSample128_failure_bound {n m : Nat}
    (circuit : Circuit n m) :
    (PMF.uniformOfFintype (Sample (varBound circuit.raw) 128)).toOuterMeasure
      {sample | ¬(phaseFoldNonlinearWithSample 128 circuit sample).Equivalent circuit} ≤
      ((nonlinearComparisonPairs
        (nonlinearRelevantCircuit circuit.raw)).card : ENNReal) *
        ((Fingerprint.maxDegree : ENNReal) / ((2 ^ 128 : Nat) : ENNReal)) := by
  calc
    _ ≤ (PMF.uniformOfFintype (Sample (varBound circuit.raw) 128)).toOuterMeasure
        (nonlinearBadSamples circuit.raw) := by
          apply (PMF.uniformOfFintype _).toOuterMeasure_mono
          intro sample hfail
          by_contra hgood
          exact hfail.1 (phaseFoldNonlinearWithSample128_correct circuit sample hgood)
    _ = (PMF.uniformOfFintype
          (Fin (varBound circuit.raw) → Fin (2 ^ 128))).toOuterMeasure
          (packedBadValuations (varBound circuit.raw)
            (nonlinearRelevantCircuit circuit.raw)) :=
          nonlinearBadSamples_probability_eq circuit.raw
    _ ≤ _ := packedCircuit_badProbability_bound circuit.raw

/-- A successful nonlinear merge only replaces a later rotation on its existing wire. -/
theorem mergeIntoNonlinear_inRange {n m : Nat} (draws : Nat → Tag)
    (st : NState) (tag : Tag) (angle : ℚ) (gates result : List Gate)
    (hin : ∀ gate ∈ gates, gate.InRange n m)
    (hmerge : mergeIntoNonlinear draws st tag angle gates = some result) :
    ∀ gate ∈ result, gate.InRange n m := by
  obtain ⟨middle, rest, later, φ, q, sign, hgates, hresult, -, hrot, -⟩ :=
    mergeIntoNonlinear_spec draws tag angle gates result st hmerge
  intro gate hgate
  rw [hresult] at hgate
  rcases List.mem_append.1 hgate with hmiddle | htail
  · exact hin gate (by rw [hgates]; simp [hmiddle])
  · rcases List.mem_cons.1 htail with rfl | hrest
    · exact Gate.InRange.onWire (hin later (by rw [hgates]; simp))
        (rotAngle_mem hrot) rfl rfl
    · exact hin gate (by rw [hgates]; simp [hrest])

/-- The nonlinear fold preserves every operand bound, independently of the samples. -/
theorem foldFromNonlinear_inRange {n m : Nat} (draws : Nat → Tag)
    (targets : Array Bool) :
    ∀ (fuel : Nat) (gates : List Gate), gates.length ≤ fuel →
      ∀ (at_ : Nat) (st : NState),
        (∀ gate ∈ gates, gate.InRange n m) →
          ∀ gate ∈ foldFromNonlinear draws targets st at_ gates,
            gate.InRange n m := by
  intro fuel
  induction fuel with
  | zero =>
      intro gates hlen at_ st hin gate hgate
      rw [List.eq_nil_of_length_eq_zero (Nat.le_zero.1 hlen)] at hgate
      simp [foldFromNonlinear] at hgate
  | succ fuel ih =>
      intro gates hlen at_ st hin
      cases gates with
      | nil => intro gate hgate; simp [foldFromNonlinear] at hgate
      | cons head tail =>
          have htail : tail.length ≤ fuel := by
            simp only [List.length_cons] at hlen
            omega
          have hkeep : ∀ gate ∈ head ::
              foldFromNonlinear draws targets (st.step draws head) (at_ + 1) tail,
              gate.InRange n m := by
            intro gate hgate
            rcases List.mem_cons.1 hgate with rfl | hgate
            · exact hin _ (by simp)
            · exact ih tail htail _ _ (fun g hg => hin g (by simp [hg])) gate hgate
          simp only [foldFromNonlinear]
          split
          · split
            · exact ih tail htail _ _ (fun g hg => hin g (by simp [hg]))
            · split
              · split
                · rename_i angle q hrot hconstant htarget result hmerge
                  have hresult : result.length ≤ fuel := by
                    rw [mergeIntoNonlinear_length draws (st.tagOf q) angle tail result st
                      hmerge]
                    exact htail
                  exact ih result hresult (at_ + 1) st
                    (mergeIntoNonlinear_inRange draws st (st.tagOf q) angle tail result
                      (fun g hg => hin g (by simp [hg])) hmerge)
                · exact hkeep
              · exact hkeep
          · exact hkeep

theorem phaseFoldGatesNonlinear_inRange {n m : Nat} (draws : Nat → Tag)
    (qubits : Nat) {gates : List Gate}
    (hin : ∀ gate ∈ gates, gate.InRange n m) :
    ∀ gate ∈ phaseFoldGatesNonlinear draws qubits gates, gate.InRange n m :=
  emitAll_inRange (foldFromNonlinear_inRange draws _ gates.length gates le_rfl 0 _ hin)

/-- The verified nonlinear pass now runs exactly the sampled transformation used by the
executable. Uniformity of the OS-backed draws remains an environmental assumption. -/
noncomputable def PhaseFoldRandNonlinear128 : RandPass where
  name := "Phase folding"
  Seed := fun circuit => Sample (varBound circuit.raw) 128
  dist := fun _ => PMF.uniformOfFintype _
  run := phaseFoldNonlinearWithSample 128
  error := fun circuit =>
    ((nonlinearComparisonPairs
      (nonlinearRelevantCircuit circuit.raw)).card : ENNReal) *
      ((Fingerprint.maxDegree : ENNReal) / ((2 ^ 128 : Nat) : ENNReal))
  wellFormed_run _circuit sample hin :=
    phaseFoldGatesNonlinear_inRange (wordsOf 128 (liftSample sample)) _ hin
  flagsOk_run _circuit _sample _ := RawCircuit.flagsOk_withGates _ _
  correct circuit := phaseFoldNonlinearWithSample128_failure_bound circuit

theorem PhaseFoldRandNonlinear128_run (circuit : Circuit n m)
    (sample : (PhaseFoldRandNonlinear128).Seed circuit) :
    (PhaseFoldRandNonlinear128).run circuit sample =
      phaseFoldNonlinearWithSample 128 circuit sample := rfl

end TzapLean.GF128Proof

namespace TzapLean

/-- The verified 128-bit nonlinear phase-folding pass used by the optimizer. -/
noncomputable abbrev PhaseFoldRand : RandPass := GF128Proof.PhaseFoldRandNonlinear128

@[simp] theorem PhaseFoldRand_run (circuit : Circuit n m)
    (sample : PhaseFoldRand.Seed circuit) :
    PhaseFoldRand.run circuit sample =
      phaseFoldNonlinearWithSample 128 circuit sample := rfl

@[simp] theorem PhaseFoldRand_error (circuit : Circuit n m) :
    PhaseFoldRand.error circuit =
      ((nonlinearComparisonPairs
        (nonlinearRelevantCircuit circuit.raw)).card : ENNReal) *
        ((Fingerprint.maxDegree : ENNReal) / ((2 ^ 128 : Nat) : ENNReal)) := rfl

end TzapLean
