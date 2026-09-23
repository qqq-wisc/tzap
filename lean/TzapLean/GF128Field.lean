import TzapLean.GF128Certificate
import TzapLean.NonlinearProbability
import Mathlib.FieldTheory.Finiteness

/-!
# The abstract 128-bit fingerprint field

The Rabin certificate makes the quotient by `X^128 + X^7 + X^2 + X + 1` a field.
Its cardinality is exactly `2^128`. `GF128Bridge.lean` relates this quotient to the
packed `Nat` operations in `GF128.lean`.
-/

namespace TzapLean.GF128Proof

noncomputable instance : Fact (Irreducible modulus) := ⟨modulus_irreducible⟩

abbrev Field128 := AdjoinRoot modulus

noncomputable instance : CharP Field128 2 :=
  CharP.of_ringHom_of_ne_zero (AdjoinRoot.of modulus) 2 (by decide)

noncomputable instance : Module.Finite F₂ Field128 := modulus_monic.finite_adjoinRoot
noncomputable instance : Finite Field128 := Module.finite_of_finite F₂

noncomputable instance : Fintype Field128 := Fintype.ofFinite Field128
noncomputable instance : DecidableEq Field128 := Classical.decEq Field128

theorem field128_card : Fintype.card Field128 = 2 ^ 128 := by
  rw [Module.card_eq_pow_finrank (K := F₂) (V := Field128)]
  have hfin : Module.finrank F₂ Field128 = modulus.natDegree :=
    finrank_quotient_span_eq_natDegree
  rw [hfin, modulus_natDegree]
  norm_num [F₂, Field128]

/-- The sharp field-level union bound for one nonlinear phase-folding run. The per-pair
factor is `2^32 / 2^128 = 2^-96`. -/
theorem field128_badProbability_bound (circuit : RawCircuit) :
    (PMF.uniformOfFintype (Fin (varBound circuit) → Field128)).toOuterMeasure
      (↑(nonlinearBadValuations (K := Field128) (varBound circuit)
        (nonlinearRelevantCircuit circuit)) : Set (Fin (varBound circuit) → Field128)) ≤
      ((nonlinearComparisonPairs (nonlinearRelevantCircuit circuit)).card : ENNReal) *
        ((Fingerprint.maxDegree : ENNReal) / ((2 ^ 128 : Nat) : ENNReal)) := by
  simpa only [field128_card] using
    nonlinearCircuit_badProbability_bound (K := Field128) circuit

end TzapLean.GF128Proof
