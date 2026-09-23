import Mathlib.Algebra.MvPolynomial.SchwartzZippel
import TzapLean.NonlinearCompared
import TzapLean.Hash

/-!
# Finite-field collision bound for nonlinear fingerprints

Schwartz–Zippel bounds each formal polynomial collision. A finite union bound then controls
all comparisons in one circuit, and a good valuation is faithful on that comparison set.
The concrete packed representation is proved separately.
-/

namespace TzapLean

noncomputable section

open Finset Fintype

set_option maxHeartbeats 0

private theorem totalDegree_rename_injective {σ τ R : Type*} [CommRing R]
    [DecidableEq τ] (f : σ → τ) (hf : Function.Injective f)
    (polynomial : MvPolynomial σ R) :
    (MvPolynomial.rename f polynomial).totalDegree = polynomial.totalDegree := by
  classical
  simp only [MvPolynomial.totalDegree, MvPolynomial.support_rename_of_injective hf,
    Finset.sup_image]
  apply Finset.sup_congr rfl
  intro monomial _
  simp [Finsupp.sum_mapDomain_index]

private theorem totalDegree_map_injective {σ R S : Type*} [CommRing R] [CommRing S]
    (f : R →+* S) (hf : Function.Injective f) (polynomial : MvPolynomial σ R) :
    (MvPolynomial.map f polynomial).totalDegree = polynomial.totalDegree := by
  classical
  simp only [MvPolynomial.totalDegree, MvPolynomial.support_map_of_injective _ hf]

/-- A polynomial whose variables are below `m` can be viewed as a polynomial in `Fin m`
variables without changing its total degree. -/
theorem boundedPolynomial_fin {m : Nat} {polynomial : BoolPolynomial}
    (hbounded : BoolPolynomial.Bounded m polynomial) :
    ∃ finitePolynomial : MvPolynomial (Fin m) F₂,
      MvPolynomial.rename (fun i : Fin m => (i : Nat)) finitePolynomial = polynomial ∧
        finitePolynomial.totalDegree = polynomial.totalDegree := by
  classical
  obtain ⟨finitePolynomial, hrename⟩ :=
    MvPolynomial.exists_rename_eq_of_vars_subset_range polynomial
      (fun i : Fin m => (i : Nat)) Fin.val_injective (by
        intro i hi
        exact ⟨⟨i, hbounded i hi⟩, rfl⟩)
  refine ⟨finitePolynomial, hrename, ?_⟩
  rw [← hrename, totalDegree_rename_injective _ Fin.val_injective]

theorem evalAt_rename_fin {K : Type*} [CommRing K] [CharP K 2]
    {m : Nat} (finitePolynomial : MvPolynomial (Fin m) F₂)
    (valuation : Fin m → K) :
    BoolPolynomial.evalAt
        (fun i => if h : i < m then valuation ⟨i, h⟩ else 0)
        (MvPolynomial.rename (fun i : Fin m => (i : Nat)) finitePolynomial) =
      MvPolynomial.eval valuation
        (MvPolynomial.map (ZMod.castHom dvd_rfl K) finitePolynomial) := by
  unfold BoolPolynomial.evalAt
  rw [MvPolynomial.eval₂Hom_rename]
  have hvaluation :
      (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) ∘
        (fun i : Fin m => (i : Nat)) =
        valuation := by
    funext i
    simp [i.isLt]
  change MvPolynomial.eval₂Hom (ZMod.castHom dvd_rfl K)
    ((fun i => if h : i < m then valuation ⟨i, h⟩ else 0) ∘
      fun i : Fin m => (i : Nat)) finitePolynomial = _
  rw [hvaluation]
  exact MvPolynomial.eval₂_eq_eval_map (ZMod.castHom dvd_rfl K) valuation finitePolynomial

/-- For two different formal polynomials, the fraction of uniform field evaluations
at which they agree is bounded by their largest total degree divided by field size. -/
theorem finiteField_polynomial_collision {K : Type*} [Field K] [Fintype K]
    [DecidableEq K]
    (m : Nat) (p q : MvPolynomial (Fin m) K) (hne : p ≠ q) :
    (#{valuation ∈ (piFinset fun _ : Fin m => (Finset.univ : Finset K)) |
        MvPolynomial.eval valuation p = MvPolynomial.eval valuation q} : ℚ≥0) /
        (Fintype.card K : ℚ≥0) ^ m ≤
      (max p.totalDegree q.totalDegree : ℚ≥0) / (Fintype.card K : ℚ≥0) := by
  classical
  have hnonzero : p - q ≠ 0 := sub_ne_zero.mpr hne
  have hsz := MvPolynomial.schwartz_zippel_totalDegree hnonzero
    (Finset.univ : Finset K)
  have hdegree : (p - q).totalDegree ≤ max p.totalDegree q.totalDegree :=
    MvPolynomial.totalDegree_sub p q
  have hbase :
      (#{valuation ∈ (piFinset fun _ : Fin m => (Finset.univ : Finset K)) |
          MvPolynomial.eval valuation p = MvPolynomial.eval valuation q} : ℚ≥0) /
          (Fintype.card K : ℚ≥0) ^ m ≤
        ((p - q).totalDegree : ℚ≥0) / (Fintype.card K : ℚ≥0) := by
    simpa [MvPolynomial.eval_sub, sub_eq_zero] using hsz
  exact hbase.trans (by gcongr; exact_mod_cast hdegree)

/-- The same bound for the actual, naturally indexed Boolean polynomials of the analyzer.
The `Fin m` presentation makes the sample space finite without changing evaluations or degree. -/
theorem boundedPolynomial_collision {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] (m : Nat) (p q : BoolPolynomial)
    (hp : BoolPolynomial.Bounded m p) (hq : BoolPolynomial.Bounded m q)
    (hne : p ≠ q) :
    (#{valuation ∈ (piFinset fun _ : Fin m => (Finset.univ : Finset K)) |
        BoolPolynomial.evalAt
          (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) p =
        BoolPolynomial.evalAt
          (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) q} : ℚ≥0) /
        (Fintype.card K : ℚ≥0) ^ m ≤
      (max p.totalDegree q.totalDegree : ℚ≥0) / (Fintype.card K : ℚ≥0) := by
  classical
  obtain ⟨pf, hpf, hpd⟩ := boundedPolynomial_fin hp
  obtain ⟨qf, hqf, hqd⟩ := boundedPolynomial_fin hq
  let f : F₂ →+* K := ZMod.castHom dvd_rfl K
  have hmapne : MvPolynomial.map f pf ≠ MvPolynomial.map f qf := by
    intro heq
    have hpq : pf = qf := MvPolynomial.map_injective f
      (ZMod.castHom_injective (R := K)) heq
    exact hne (hpf ▸ hqf ▸ congrArg (MvPolynomial.rename
      (fun i : Fin m => (i : Nat))) hpq)
  have hbound := finiteField_polynomial_collision m
    (MvPolynomial.map f pf) (MvPolynomial.map f qf) hmapne
  have hdegrees :
      (MvPolynomial.map f pf).totalDegree = p.totalDegree ∧
        (MvPolynomial.map f qf).totalDegree = q.totalDegree := by
    constructor
    · exact (totalDegree_map_injective f (ZMod.castHom_injective (R := K)) pf).trans hpd
    · exact (totalDegree_map_injective f (ZMod.castHom_injective (R := K)) qf).trans hqd
  simpa only [← hpf, ← hqf, evalAt_rename_fin, hdegrees.1, hdegrees.2] using hbound

/-- Pairs of distinct polynomials whose evaluations must not collide. -/
def nonlinearComparisonPairs (polynomials : List BoolPolynomial) :
    Finset (BoolPolynomial × BoolPolynomial) :=
  ((polynomials.toFinset.product polynomials.toFinset).filter fun pair => pair.1 ≠ pair.2)

/-- The finite set of sampled valuations on which one comparison is a false positive. -/
def nonlinearCollisionValuations {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] (m : Nat) (pair : BoolPolynomial × BoolPolynomial) :
    Finset (Fin m → K) :=
  (piFinset fun _ : Fin m => (Finset.univ : Finset K)).filter fun valuation =>
    BoolPolynomial.evalAt
      (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) pair.1 =
    BoolPolynomial.evalAt
      (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) pair.2

/-- A seed is bad if any two different polynomials in the circuit's comparison set collide. -/
def nonlinearBadValuations {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] (m : Nat) (polynomials : List BoolPolynomial) :
    Finset (Fin m → K) :=
  (nonlinearComparisonPairs polynomials).biUnion
    (nonlinearCollisionValuations (K := K) m)

theorem nonlinearBadValuations_bound {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] (m degreeBound : Nat)
    (polynomials : List BoolPolynomial)
    (hvars : ∀ p ∈ polynomials, BoolPolynomial.Bounded m p)
    (hdegree : ∀ p ∈ polynomials, p.totalDegree ≤ degreeBound) :
    ((nonlinearBadValuations (K := K) m polynomials).card : ℚ≥0) /
        (Fintype.card K : ℚ≥0) ^ m ≤
      (nonlinearComparisonPairs polynomials).card *
        ((degreeBound : ℚ≥0) / (Fintype.card K : ℚ≥0)) := by
  classical
  let pairs := nonlinearComparisonPairs polynomials
  let collisions := nonlinearCollisionValuations (K := K) m
  have hcount :
      ((pairs.biUnion collisions).card : ℚ≥0) ≤
        ∑ pair ∈ pairs, (collisions pair).card := by
    exact_mod_cast (Finset.card_biUnion_le (s := pairs) (t := collisions))
  calc
    ((nonlinearBadValuations (K := K) m polynomials).card : ℚ≥0) /
        (Fintype.card K : ℚ≥0) ^ m =
      ((pairs.biUnion collisions).card : ℚ≥0) /
        (Fintype.card K : ℚ≥0) ^ m := rfl
    _ ≤ (∑ pair ∈ pairs, ((collisions pair).card : ℚ≥0)) /
        (Fintype.card K : ℚ≥0) ^ m := by
      exact div_le_div_of_nonneg_right (by simpa only [Nat.cast_sum] using hcount)
        (by positivity)
    _ = ∑ pair ∈ pairs, (((collisions pair).card : ℚ≥0) /
        (Fintype.card K : ℚ≥0) ^ m) := by rw [Finset.sum_div]
    _ ≤ ∑ _pair ∈ pairs, ((degreeBound : ℚ≥0) / (Fintype.card K : ℚ≥0)) := by
      apply Finset.sum_le_sum
      intro pair hpair
      have hpair' : (pair.1 ∈ polynomials ∧ pair.2 ∈ polynomials) ∧ pair.1 ≠ pair.2 := by
        simpa [pairs, nonlinearComparisonPairs] using hpair
      obtain ⟨⟨hp, hq⟩, hne⟩ := hpair'
      have hcollision := boundedPolynomial_collision (K := K) m pair.1 pair.2
        (hvars _ hp) (hvars _ hq) hne
      have hmax : max pair.1.totalDegree pair.2.totalDegree ≤ degreeBound :=
        max_le (hdegree _ hp) (hdegree _ hq)
      have hcollision' :
          (((collisions pair).card : ℚ≥0) / (Fintype.card K : ℚ≥0) ^ m) ≤
            (max pair.1.totalDegree pair.2.totalDegree : ℚ≥0) /
              (Fintype.card K : ℚ≥0) := by
        simpa only [collisions, nonlinearCollisionValuations] using hcollision
      exact hcollision'.trans (by gcongr; exact_mod_cast hmax)
    _ = (nonlinearComparisonPairs polynomials).card *
        ((degreeBound : ℚ≥0) / (Fintype.card K : ℚ≥0)) := by
      simp [pairs]

/-- Whole-circuit collision bound, including CCX products and fresh cutoff variables. -/
theorem nonlinearCircuit_badValuations_bound {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] (circuit : RawCircuit) :
    ((nonlinearBadValuations (K := K) (varBound circuit)
      (nonlinearRelevantCircuit circuit)).card : ℚ≥0) /
        (Fintype.card K : ℚ≥0) ^ (varBound circuit) ≤
      (nonlinearComparisonPairs (nonlinearRelevantCircuit circuit)).card *
        ((Fingerprint.maxDegree : ℚ≥0) / (Fintype.card K : ℚ≥0)) := by
  exact nonlinearBadValuations_bound (varBound circuit) Fingerprint.maxDegree
    (nonlinearRelevantCircuit circuit)
    (nonlinearRelevantCircuit_variablesBounded circuit)
    (nonlinearRelevantCircuit_degreeBounded circuit)

theorem nonlinearCircuit_badProbability_bound {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] (circuit : RawCircuit) :
    (PMF.uniformOfFintype (Fin (varBound circuit) → K)).toOuterMeasure
      (↑(nonlinearBadValuations (K := K) (varBound circuit)
        (nonlinearRelevantCircuit circuit)) : Set (Fin (varBound circuit) → K)) ≤
      ((nonlinearComparisonPairs (nonlinearRelevantCircuit circuit)).card : ENNReal) *
        ((Fingerprint.maxDegree : ENNReal) / (Fintype.card K : ENNReal)) := by
  classical
  rw [PMF.toOuterMeasure_uniformOfFintype_apply, Fintype.card_subtype]
  have hcard : Fintype.card (Fin (varBound circuit) → K) =
      (Fintype.card K) ^ (varBound circuit) := by simp
  rw [hcard]
  have hset :
      #{x : Fin (varBound circuit) → K |
        x ∈ (↑(nonlinearBadValuations (K := K) (varBound circuit)
          (nonlinearRelevantCircuit circuit)) : Set (Fin (varBound circuit) → K))} =
        (nonlinearBadValuations (K := K) (varBound circuit)
          (nonlinearRelevantCircuit circuit)).card := by simp
  rw [hset]
  have hrat := nonlinearCircuit_badValuations_bound (K := K) circuit
  have hreal :
      (((nonlinearBadValuations (K := K) (varBound circuit)
        (nonlinearRelevantCircuit circuit)).card : NNReal) /
        (Fintype.card K : NNReal) ^ (varBound circuit)) ≤
      ((nonlinearComparisonPairs (nonlinearRelevantCircuit circuit)).card : NNReal) *
        ((Fingerprint.maxDegree : NNReal) / (Fintype.card K : NNReal)) := by
    simpa only [NNRat.cast_div, NNRat.cast_mul, NNRat.cast_natCast,
      NNRat.cast_pow, Nat.cast_pow] using
      (NNRat.cast_mono (K := NNReal) hrat)
  have hK : (Fintype.card K : NNReal) ≠ 0 := by positivity
  have hcast := ENNReal.coe_le_coe.mpr hreal
  simpa only [ENNReal.coe_div (by positivity :
      (Fintype.card K : NNReal) ^ varBound circuit ≠ 0),
    ENNReal.coe_div hK, ENNReal.coe_mul, ENNReal.coe_pow,
    ENNReal.coe_natCast, Nat.cast_pow] using hcast

/-- Outside the bad set, a field evaluation is injective on every polynomial the circuit
may compare. -/
theorem nonlinearGoodValuation_faithful {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] (m : Nat) (polynomials : List BoolPolynomial)
    (valuation : Fin m → K)
    (hgood : valuation ∉ nonlinearBadValuations (K := K) m polynomials) :
    ∀ p ∈ polynomials, ∀ q ∈ polynomials,
      BoolPolynomial.evalAt
        (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) p =
      BoolPolynomial.evalAt
        (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) q → p = q := by
  intro p hp q hq heq
  by_contra hne
  apply hgood
  simp only [nonlinearBadValuations, Finset.mem_biUnion]
  refine ⟨(p, q), ?_, ?_⟩
  · simp [nonlinearComparisonPairs, hp, hq, hne]
  · simp [nonlinearCollisionValuations, heq]

/-- A correct, injective packed encoding transfers the field-level good-seed statement to
the precise faithfulness hypothesis used by `matchFingerprint_exact`. -/
theorem nonlinearGoodValuation_packedFaithful {K : Type*} [Field K] [Fintype K]
    [DecidableEq K] [CharP K 2] {draws : Nat → Tag}
    (m : Nat) (polynomials : List BoolPolynomial) (valuation : Fin m → K)
    (evaluation : PackedEvaluation draws) (encode : K → Tag)
    (hencode : Function.Injective encode)
    (hagree : ∀ p ∈ polynomials,
      evaluation.eval p = encode (BoolPolynomial.evalAt
        (fun i => if h : i < m then valuation ⟨i, h⟩ else 0) p))
    (hgood : valuation ∉ nonlinearBadValuations (K := K) m polynomials) :
    PolynomialFaithful evaluation polynomials := by
  intro p hp q hq heq
  apply nonlinearGoodValuation_faithful m polynomials valuation hgood p hp q hq
  apply hencode
  rw [← hagree p hp, ← hagree q hq]
  exact heq

end

end TzapLean
