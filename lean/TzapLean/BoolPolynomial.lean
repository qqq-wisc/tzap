import Mathlib.Algebra.MvPolynomial.Variables
import TzapLean.Hash

/-!
# Boolean functions through formal polynomial lifts

The nonlinear phase analysis does not evaluate in the Boolean quotient.  A wire carries an
ordinary polynomial in `F₂[x₀,x₁,…]`; restricting its variables to `0` and `1` gives the
Boolean value of the wire.  This matters because random evaluation happens in a larger field,
where the Boolean relation `x² = x` is deliberately not imposed.
-/

namespace TzapLean

noncomputable section

/-- Formal polynomial lifts of Boolean wire values. -/
abbrev BoolPolynomial := MvPolynomial Nat F₂

namespace BoolPolynomial

def const (b : Bool) : BoolPolynomial := MvPolynomial.C (bit b)

def var (i : Nat) : BoolPolynomial := MvPolynomial.X i

/-- Evaluation of a formal lift at a Boolean valuation, still represented in `F₂`. -/
def evalF₂ (valuation : Nat → Bool) (p : BoolPolynomial) : F₂ :=
  MvPolynomial.eval (fun i => bit (valuation i)) p

/-- Boolean interpretation of a formal lift. -/
def evalB (valuation : Nat → Bool) (p : BoolPolynomial) : Bool :=
  unbit (evalF₂ valuation p)

/-- Evaluate the formal lift in any characteristic-two commutative ring.  In particular, the
nonlinear collision proof will instantiate `K` with the 128-bit fingerprint field. -/
def evalAt {K : Type*} [CommRing K] [CharP K 2]
    (valuation : Nat → K) (p : BoolPolynomial) : K :=
  MvPolynomial.eval₂Hom (ZMod.castHom dvd_rfl K) valuation p

@[simp] theorem evalAt_const {K : Type*} [CommRing K] [CharP K 2]
    (v : Nat → K) (b : Bool) :
    evalAt v (const b) = if b then 1 else 0 := by
  cases b <;> simp [evalAt, const, bit]

@[simp] theorem evalAt_var {K : Type*} [CommRing K] [CharP K 2]
    (v : Nat → K) (i : Nat) : evalAt v (var i) = v i := by
  simp [evalAt, var]

@[simp] theorem evalAt_add {K : Type*} [CommRing K] [CharP K 2]
    (v : Nat → K) (p q : BoolPolynomial) :
    evalAt v (p + q) = evalAt v p + evalAt v q := by
  simp [evalAt]

@[simp] theorem evalAt_mul {K : Type*} [CommRing K] [CharP K 2]
    (v : Nat → K) (p q : BoolPolynomial) :
    evalAt v (p * q) = evalAt v p * evalAt v q := by
  simp [evalAt]

@[simp] theorem evalF₂_const (v : Nat → Bool) (b : Bool) : evalF₂ v (const b) = bit b := by
  simp [evalF₂, const]

@[simp] theorem evalF₂_var (v : Nat → Bool) (i : Nat) : evalF₂ v (var i) = bit (v i) := by
  simp [evalF₂, var]

@[simp] theorem evalF₂_add (v : Nat → Bool) (p q : BoolPolynomial) :
    evalF₂ v (p + q) = evalF₂ v p + evalF₂ v q := by
  simp [evalF₂]

@[simp] theorem evalF₂_mul (v : Nat → Bool) (p q : BoolPolynomial) :
    evalF₂ v (p * q) = evalF₂ v p * evalF₂ v q := by
  simp [evalF₂]

@[simp] theorem evalB_const (v : Nat → Bool) (b : Bool) : evalB v (const b) = b := by
  simp [evalB]

@[simp] theorem evalB_var (v : Nat → Bool) (i : Nat) : evalB v (var i) = v i := by
  simp [evalB]

@[simp] theorem evalB_add (v : Nat → Bool) (p q : BoolPolynomial) :
    evalB v (p + q) = (evalB v p != evalB v q) := by
  apply Bool.eq_iff_iff.mpr
  simp only [evalB, evalF₂_add]
  generalize evalF₂ v p = a
  generalize evalF₂ v q = b
  revert a b
  decide

@[simp] theorem evalB_mul (v : Nat → Bool) (p q : BoolPolynomial) :
    evalB v (p * q) = (evalB v p && evalB v q) := by
  apply Bool.eq_iff_iff.mpr
  simp only [evalB, evalF₂_mul]
  generalize evalF₂ v p = a
  generalize evalF₂ v q = b
  revert a b
  decide

/-- Complement a Boolean polynomial without quotienting by Boolean identities. -/
def flip (p : BoolPolynomial) : BoolPolynomial := p + 1

@[simp] theorem evalB_flip (v : Nat → Bool) (p : BoolPolynomial) :
    evalB v p.flip = !evalB v p := by
  rw [flip, evalB_add]
  have hone : unbit (1 : F₂) = true := by decide
  simp [evalB, evalF₂, hone]

/-- A polynomial mentions only variables below `m`. -/
def Bounded (m : Nat) (p : BoolPolynomial) : Prop :=
  ∀ i ∈ p.vars, i < m

theorem bounded_const (m : Nat) (b : Bool) : Bounded m (const b) := by
  intro i hi
  simp [const] at hi

theorem bounded_zero (m : Nat) : Bounded m (0 : BoolPolynomial) := by
  intro i hi
  simp at hi

theorem bounded_one (m : Nat) : Bounded m (1 : BoolPolynomial) := by
  intro i hi
  simp at hi

theorem bounded_var {m i : Nat} (h : i < m) : Bounded m (var i) := by
  intro j hj
  have : j = i := by simpa [var, MvPolynomial.vars_X] using hj
  simpa [this] using h

theorem bounded_add {m : Nat} {p q : BoolPolynomial} (hp : Bounded m p) (hq : Bounded m q) :
    Bounded m (p + q) := by
  intro i hi
  have hi' : i ∈ p.vars ∪ q.vars := MvPolynomial.vars_add_subset p q hi
  rcases Finset.mem_union.mp hi' with hi | hi
  · exact hp i hi
  · exact hq i hi

theorem bounded_mul {m : Nat} {p q : BoolPolynomial} (hp : Bounded m p) (hq : Bounded m q) :
    Bounded m (p * q) := by
  intro i hi
  have hi' : i ∈ p.vars ∪ q.vars := MvPolynomial.vars_mul p q hi
  rcases Finset.mem_union.mp hi' with hi | hi
  · exact hp i hi
  · exact hq i hi

theorem bounded_mono {m m' : Nat} {p : BoolPolynomial} (h : m ≤ m') (hp : Bounded m p) :
    Bounded m' p := fun i hi => lt_of_lt_of_le (hp i hi) h

/-- A bounded polynomial has the same Boolean value under valuations that agree on all
variables it can mention. -/
theorem evalB_congr {m : Nat} {p : BoolPolynomial} (hp : Bounded m p)
    {v v' : Nat → Bool} (h : ∀ i, i < m → v' i = v i) :
    evalB v' p = evalB v p := by
  apply congrArg unbit
  simp only [evalF₂]
  apply MvPolynomial.eval₂_congr
  intro i exponent hi hcoeff
  apply congrArg bit
  exact h i (hp i ((MvPolynomial.mem_vars_iff_mem_support i).2
    ⟨exponent, MvPolynomial.mem_support_iff.mpr hcoeff, hi⟩))

theorem totalDegree_const (b : Bool) : (const b).totalDegree = 0 := by
  simp [const]

theorem totalDegree_var (i : Nat) : (var i).totalDegree = 1 := by
  simp [var]

theorem totalDegree_add_le (p q : BoolPolynomial) :
    (p + q).totalDegree ≤ max p.totalDegree q.totalDegree :=
  MvPolynomial.totalDegree_add p q

theorem totalDegree_mul_le (p q : BoolPolynomial) :
    (p * q).totalDegree ≤ p.totalDegree + q.totalDegree :=
  MvPolynomial.totalDegree_mul p q

/-- Formal equality is sufficient for equality on every Boolean valuation.  The converse is
intentionally not claimed: `x²` and `x` are different lifts of the same Boolean function. -/
theorem evalB_eq_of_eq {p q : BoolPolynomial} (h : p = q) (v : Nat → Bool) :
    evalB v p = evalB v q := by rw [h]

end BoolPolynomial

end

end TzapLean
