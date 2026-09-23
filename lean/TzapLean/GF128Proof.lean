import Mathlib.FieldTheory.Finite.Extension
import Mathlib.RingTheory.Polynomial.UniqueFactorization
import Mathlib.RingTheory.Coprime.Basic
import Mathlib.Tactic.ComputeDegree

/-!
# Algebraic certificate boundary for the packed fingerprint field

The packed implementation reduces by `x^128 + x^7 + x^2 + x + 1`.  The theorem below
specializes Rabin's irreducibility criterion to degree 128: because the only prime divisor of
128 is 2, it is enough to certify the Frobenius-128 split and one gcd at Frobenius-64.

The remaining certificate is deliberately computational: repeated squaring modulo `modulus`
must prove the two hypotheses.  Once supplied, `irreducible_of_rabin_128` yields the quotient
field used by the Schwartz–Zippel proof.
-/

set_option maxRecDepth 10000

namespace TzapLean.GF128Proof

open Polynomial

abbrev F₂ := ZMod 2

noncomputable def modulus : F₂[X] := X ^ 128 + X ^ 7 + X ^ 2 + X + 1

theorem modulus_monic : modulus.Monic := by
  unfold modulus
  monicity!

theorem modulus_natDegree : modulus.natDegree = 128 := by
  unfold modulus
  compute_degree!

theorem irreducible_of_rabin_128 {k : Type*} [Field k] [Finite k]
    {f : k[X]} (hmonic : f.Monic) (hdeg : f.natDegree = 128)
    (hsplit : f ∣ X ^ (Nat.card k) ^ 128 - X)
    (hcop : IsCoprime f (X ^ (Nat.card k) ^ 64 - X)) : Irreducible f := by
  rw [hmonic.irreducible_iff_lt_natDegree_lt]
  · intro q hqmonic hqdeg hqf
    have hqbounds : 0 < q.natDegree ∧ q.natDegree ≤ f.natDegree / 2 := by
      simpa only [Finset.mem_Ioc] using hqdeg
    obtain ⟨p, hp, hpq⟩ := q.exists_irreducible_of_natDegree_pos hqbounds.1
    have hpdeg_le_q : p.natDegree ≤ q.natDegree :=
      natDegree_le_of_dvd hpq hqmonic.ne_zero
    have hpdeg_le : p.natDegree ≤ 64 := by
      rw [hdeg] at hqbounds
      omega
    have hpf : p ∣ f := dvd_trans hpq hqf
    have hp_split : p ∣ X ^ (Nat.card k) ^ 128 - X := dvd_trans hpf hsplit
    have hpdeg_dvd128 : p.natDegree ∣ 128 :=
      hp.natDegree_dvd_of_dvd_X_pow_card_pow_sub_X hp_split
    have hpdeg_dvd64 : p.natDegree ∣ 64 := by
      have hd : p.natDegree ∣ 2 ^ 7 := by
        norm_num
        exact hpdeg_dvd128
      obtain ⟨e, he7, heq⟩ := (Nat.dvd_prime_pow Nat.prime_two).mp hd
      have he6 : e ≤ 6 := by
        by_contra h
        have he : e = 7 := by omega
        subst e
        norm_num at heq
        omega
      have hd64 : p.natDegree ∣ 2 ^ 6 :=
        (Nat.dvd_prime_pow Nat.prime_two).mpr ⟨e, he6, heq⟩
      norm_num at hd64 ⊢
      exact hd64
    have hp64 : p ∣ X ^ (Nat.card k) ^ 64 - X :=
      (hp.natDegree_dvd_iff_dvd_X_pow_card_pow_sub_X).1 hpdeg_dvd64
    exact hp.not_isUnit (hcop.isUnit_of_dvd' hpf hp64)
  · intro hf
    rw [hf] at hdeg
    simp at hdeg

end TzapLean.GF128Proof
