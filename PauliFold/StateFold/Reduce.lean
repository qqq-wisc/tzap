import PauliFold.StateFold.Exact
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Algebra.MvPolynomial.Monad

/-!
# StateFold's reductions

After the gates, StateFold simplifies the path sum by eliminating path variables that appear
in no wire (Feynman's `applyReductions`):

* **HH.** If `Φ = Φ₀ + π · y · (z ⊕ R)` with `Φ₀` independent of `y` and `R` independent of
  `y` and `z`, then summing over `y` gives `2 · [z = R]`: eliminate `y`, and substitute
  `z := R` everywhere. `statefold d` requires `deg R ≤ d`.
* **ω.** If `Φ = Φ₀ + π · y · (1/2 + R)`, then summing over `y` gives
  `√2 · e^{iπ/4} · e^{-iπR/2}`: eliminate `y` and add `1/4 - R/2` to the phase.

Both are stated as relations between states and proved to preserve the denotation. The degree
bound appears only as a side condition, so every `d` is sound; `d` affects only which folds are
found.
-/

namespace PauliFold.SF

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Substitution -/

/-- Substitute `R` for the variable `z`. -/
def subst (z : ℕ) (R p : Poly) : Poly :=
  MvPolynomial.bind₁ (fun v => if v = z then R else MvPolynomial.X v) p

theorem ev_subst (ν : Val) (z : ℕ) (R p : Poly) :
    ev ν (subst z R p) = ev (Function.update ν z (ev ν R)) p := by
  unfold ev subst BoolPolynomial.evalB BoolPolynomial.evalF₂
  congr 1
  rw [show MvPolynomial.eval (fun i => bit (ν i))
      (MvPolynomial.bind₁ (fun v => if v = z then R else MvPolynomial.X v) p) =
    MvPolynomial.eval (fun v => MvPolynomial.eval (fun i => bit (ν i))
      (if v = z then R else MvPolynomial.X v)) p from MvPolynomial.eval₂Hom_bind₁ _ _ _ p]
  have key : (fun v => MvPolynomial.eval (fun i => bit (ν i))
      (if v = z then R else MvPolynomial.X v)) =
      fun i => bit (Function.update ν z (unbit (MvPolynomial.eval (fun i => bit (ν i)) R)) i) := by
    funext v
    by_cases hv : v = z
    · subst hv
      simp
    · simp [hv, Function.update_of_ne hv]
  rw [key]

/-! ## The rules -/

/-- The result of an HH reduction. -/
def hhResult (A : State n) (y z : ℕ) (R : Poly) (Φ₀ : Val → ℚ) : State n where
  ket := fun q => subst z R (A.ket q)
  phase := fun ν => Φ₀ (Function.update ν z (ev ν R))
  temps := (A.temps.erase y).erase z
  next := A.next
  sites := A.sites.map fun s => { s with pred := subst z R s.pred }

/-- An HH reduction at degree `d`: eliminate `y`, substituting `z := R`. -/
structure HH (d : ℕ) (A B : State n) (y z : ℕ) (R : Poly) (Φ₀ : Val → ℚ) : Prop where
  hy : y ∈ A.temps
  hz : z ∈ A.temps
  hyz : y ≠ z
  ket : ∀ q, Supp (fun ν => ev ν (A.ket q)) (A.dom \ {y})
  deg : R.totalDegree ≤ d
  suppR : Supp (fun ν => ev ν R) (A.dom \ {y, z})
  supp₀ : Supp Φ₀ (A.dom \ {y})
  split : ∀ ν, ep (A.phase ν) = ep (Φ₀ ν + b2q (ν y && (ν z != ev ν R)))
  result : B = hhResult A y z R Φ₀

/-- The result of an ω reduction. -/
def omegaResult (A : State n) (y : ℕ) (R : Poly) (Φ₀ : Val → ℚ) : State n :=
  { A with
    phase := fun ν => Φ₀ ν + 1 / 4 - 1 / 2 * b2q (ev ν R)
    temps := A.temps.erase y }

/-- An ω reduction: eliminate `y`. -/
structure Omega (A B : State n) (y : ℕ) (R : Poly) (Φ₀ : Val → ℚ) : Prop where
  hy : y ∈ A.temps
  ket : ∀ q, Supp (fun ν => ev ν (A.ket q)) (A.dom \ {y})
  suppR : Supp (fun ν => ev ν R) (A.dom \ {y})
  supp₀ : Supp Φ₀ (A.dom \ {y})
  split : ∀ ν, ep (A.phase ν) = ep (Φ₀ ν + b2q (ν y) * (1 / 2 + b2q (ev ν R)))
  result : B = omegaResult A y R Φ₀

/-! ## Arithmetic -/

theorem invSqrt2_sq : invSqrt2 * invSqrt2 = 1 / 2 := by
  have hs := sqrt2_mul_self
  have h0 : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := by
    intro h; rw [h, zero_mul] at hs; norm_num at hs
  unfold invSqrt2
  rw [← mul_inv, hs]
  norm_num

theorem ep_quarter : ep (1 / 4) = invSqrt2 * (1 + Complex.I) := by
  have : ((Real.pi * ((1 / 4 : ℚ) : ℝ) : ℝ) : ℂ) * Complex.I =
      (Real.pi / 4 : ℝ) * Complex.I := by push_cast; ring
  rw [ep, this, Complex.exp_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin,
    Real.cos_pi_div_four, Real.sin_pi_div_four]
  have hs := sqrt2_mul_self
  have h0 : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := by
    intro h; rw [h, zero_mul] at hs; norm_num at hs
  unfold invSqrt2
  push_cast
  field_simp
  linear_combination (1 + Complex.I) * hs

theorem ep_neg_quarter : ep (-1 / 4) = invSqrt2 * (1 - Complex.I) := by
  have h : ep (-1 / 4) = star (ep (1 / 4)) := by rw [star_ep]; norm_num
  rw [h, ep_quarter]
  simp [invSqrt2, Complex.conj_ofReal]
  ring

@[simp] theorem b2q_true : b2q true = 1 := rfl
@[simp] theorem b2q_false : b2q false = 0 := rfl

/-- The ω identity: `(1 + e^{iπ(1/2 + r)}) / √2 = e^{iπ(1/4 - r/2)}`. -/
theorem omega_sum (r : Bool) :
    invSqrt2 * (1 + ep (1 / 2 + b2q r)) = ep (1 / 4 - 1 / 2 * b2q r) := by
  cases r
  · simp only [b2q_false, add_zero, mul_zero, sub_zero]
    rw [ep_half, ep_quarter]
  · simp only [b2q_true, mul_one]
    rw [show (1 / 2 : ℚ) + 1 = 1 + 1 / 2 by norm_num, ep_add, ep_one, ep_half,
      show (1 / 4 : ℚ) - 1 / 2 = -1 / 4 by norm_num, ep_neg_quarter]
    ring

/-- The HH identity: `1 + (-1)^c = 2 [c = 0]`. -/
theorem hh_sum (c : Bool) : 1 + ep (b2q c) = if c then 0 else 2 := by
  cases c <;> simp [b2q, ep_one] <;> norm_num

/-! ## Splitting a path sum over one variable -/

/-- The terms of a path sum depend on the path variables only through the valuation. -/
theorem sum_split (U : Finset ℕ) (y : ℕ) (hy : y ∉ U) (hyn : n ≤ y) (i : Basis n)
    (g : Val → ℂ) :
    ∑ S ∈ (insert y U).powerset, g (val i S) =
      ∑ S ∈ U.powerset, (g (val i S) + g (Function.update (val i S) y true)) := by
  rw [Finset.sum_powerset_insert hy, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun S _ => ?_
  rw [val_insert i S y hyn]

/-- The summand of a path sum, as a function of the valuation. -/
def State.term (A : State n) (o : Basis n) (ν : Val) : ℂ :=
  if ∀ q, ev ν (A.ket q) = o q then ep (A.phase ν) else 0

theorem State.denote_eq (A : State n) (o i : Basis n) :
    A.denote o i = invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset, A.term o (val i S) :=
  rfl

theorem mem_dom_diff {A : State n} {v y : ℕ} (hv : v ∈ A.dom) (hvy : v ≠ y) :
    v ∈ A.dom \ {y} := ⟨hv, hvy⟩

/-! ## Soundness of HH -/

theorem denote_hh {d : ℕ} {A B : State n} {y z : ℕ} {R : Poly} {Φ₀ : Val → ℚ}
    (h : HH d A B y z R Φ₀) (hA : A.WF) : B.denote = A.denote := by
  obtain ⟨hy, hz, hyz, hket, -, hR, h₀, hsplit, rfl⟩ := h
  funext o i
  set U := (A.temps.erase y).erase z with hU
  have hzU : z ∉ U := Finset.notMem_erase z _
  have hzy : z ∈ A.temps.erase y := Finset.mem_erase.mpr ⟨Ne.symm hyz, hz⟩
  have hyU : y ∉ insert z U := by
    rw [Finset.mem_insert, not_or]
    exact ⟨hyz, fun h => by simp [hU] at h⟩
  have hT : A.temps = insert y (insert z U) := by
    rw [hU, Finset.insert_erase hzy, Finset.insert_erase hy]
  have hyn := hA.temps_ge y hy
  have hzn := hA.temps_ge z hz
  rw [State.denote_eq, State.denote_eq, hT, Finset.card_insert_of_notMem hyU,
    Finset.card_insert_of_notMem hzU, sum_split _ y hyU hyn,
    sum_split U z hzU hzn i (fun ν => A.term o ν + A.term o (Function.update ν y true))]
  show invSqrt2 ^ U.card * _ = _
  rw [pow_succ, pow_succ, mul_assoc, mul_assoc]
  congr 1
  rw [Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun S hS => ?_
  have hSU : S ⊆ U := Finset.mem_powerset.mp hS
  have hyS : y ∉ S := fun h => hyU (Finset.mem_insert_of_mem (hSU h))
  have hzS : z ∉ S := fun h => hzU (hSU h)
  set ν := val i S
  have hνy : ν y = false := val_not_mem i S y hyn hyS
  have hνz : ν z = false := val_not_mem i S z hzn hzS
  -- `y` is not in the domain the ket, `R`, and `Φ₀` depend on.
  have hyD : y ∉ A.dom \ {y} := fun h => h.2 rfl
  have hyzD : ∀ w, w ∈ ({y, z} : Set ℕ) → w ∉ A.dom \ {y, z} := fun w hw h => h.2 hw
  set r := ev ν R
  have hRy : ∀ μ b, ev (Function.update μ y b) R = ev μ R := fun μ b =>
    hR.update (hyzD y (by simp)) μ b
  have hRz : ∀ μ b, ev (Function.update μ z b) R = ev μ R := fun μ b =>
    hR.update (hyzD z (by simp)) μ b
  -- The term at `(y, z) = (a, b)`.
  have hterm : ∀ b : Bool, ∀ a : Bool,
      A.term o (Function.update (Function.update ν z b) y a) =
        (if ∀ q, ev (Function.update ν z b) (A.ket q) = o q then
          ep (Φ₀ (Function.update ν z b)) * ep (b2q (a && (b != r))) else 0) := by
    intro b a
    unfold State.term
    have hk : ∀ q, ev (Function.update (Function.update ν z b) y a) (A.ket q) =
        ev (Function.update ν z b) (A.ket q) := fun q => (hket q).update hyD _ _
    simp only [hk]
    split_ifs
    · rw [hsplit, h₀.update hyD, ep_add]
      have hyv : Function.update (Function.update ν z b) y a y = a := Function.update_self _ _ _
      have hzv : Function.update (Function.update ν z b) y a z = b := by
        rw [Function.update_of_ne (Ne.symm hyz), Function.update_self]
      have hRv : ev (Function.update (Function.update ν z b) y a) R = r := by
        rw [hRy, hRz]
      rw [hyv, hzv, hRv]
    · rfl
  have e1 : Function.update ν z false = ν := by
    rw [← hνz]; exact Function.update_eq_self _ _
  have e2 : Function.update ν y false = ν := by
    rw [← hνy]; exact Function.update_eq_self _ _
  have e3 : Function.update (Function.update ν z true) y false = Function.update ν z true := by
    have : Function.update ν z true y = false := by rw [Function.update_of_ne hyz, hνy]
    rw [← this]; exact Function.update_eq_self _ _
  have t00 := hterm false false
  rw [e1, e2] at t00
  have t10 := hterm false true
  rw [e1] at t10
  have t01 := hterm true false
  rw [e3] at t01
  have t11 := hterm true true
  rw [t00, t10, t01, t11]
  -- The `B` term: substitution `z := R`.
  have hB : (hhResult A y z R Φ₀).term o ν =
      if ∀ q, ev (Function.update ν z r) (A.ket q) = o q then
        ep (Φ₀ (Function.update ν z r)) else 0 := by
    unfold State.term hhResult
    simp only [ev_subst]
    rfl
  rw [hB]
  have hνz0 : Function.update ν z false = ν := e1
  rw [← mul_assoc, invSqrt2_sq]
  cases hr : r
  · simp only [hνz0, Bool.and_false, Bool.and_true, Bool.bne_false, Bool.bne_true,
      Bool.not_false, b2q_false, b2q_true, ep_zero, ep_one]
    split_ifs <;> ring
  · simp only [hνz0, Bool.and_false, Bool.and_true, Bool.bne_false, Bool.bne_true,
      Bool.not_true, Bool.not_false, b2q_false, b2q_true, ep_zero, ep_one]
    split_ifs <;> ring

/-! ## Soundness of ω -/

theorem denote_omega {A B : State n} {y : ℕ} {R : Poly} {Φ₀ : Val → ℚ}
    (h : Omega A B y R Φ₀) (hA : A.WF) : B.denote = A.denote := by
  obtain ⟨hy, hket, hR, h₀, hsplit, rfl⟩ := h
  funext o i
  set U := A.temps.erase y with hU
  have hyU : y ∉ U := Finset.notMem_erase y _
  have hT : A.temps = insert y U := (Finset.insert_erase hy).symm
  have hyn := hA.temps_ge y hy
  rw [State.denote_eq, State.denote_eq, hT, Finset.card_insert_of_notMem hyU,
    sum_split _ y hyU hyn]
  show invSqrt2 ^ U.card * _ = _
  rw [pow_succ, mul_assoc]
  congr 1
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun S hS => ?_
  have hyS : y ∉ S := fun h => hyU (Finset.mem_powerset.mp hS h)
  set ν := val i S
  have hνy : ν y = false := val_not_mem i S y hyn hyS
  have hyD : y ∉ A.dom \ {y} := fun h => h.2 rfl
  have hk : ∀ q, ev (Function.update ν y true) (A.ket q) = ev ν (A.ket q) :=
    fun q => (hket q).update hyD _ _
  have hR1 : ev (Function.update ν y true) R = ev ν R := hR.update hyD _ _
  have hΦ1 : Φ₀ (Function.update ν y true) = Φ₀ ν := h₀.update hyD _ _
  unfold State.term omegaResult
  simp only [hk]
  split_ifs
  · rw [hsplit, hsplit, hΦ1, hR1, hνy, Function.update_self, b2q_false, b2q_true, zero_mul,
      add_zero, one_mul, show Φ₀ ν + 1 / 4 - 1 / 2 * b2q (ev ν R) =
        Φ₀ ν + (1 / 4 - 1 / 2 * b2q (ev ν R)) by ring, ep_add, ep_add, ← omega_sum]
    ring
  · simp

end

end PauliFold.SF
