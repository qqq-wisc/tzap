import TzapLean.BoolPolynomial
import TzapLean.GF128
import TzapLean.Circuit
import TzapLean.PhaseFoldNonlinear

/-!
# Symbolic nonlinear analysis for `PhaseFoldRand`

This is the formal counterpart of `NState`.  Each wire carries an ordinary polynomial over
`F₂` plus a checked upper bound on its total degree.  It is intentionally not a polynomial in
the Boolean quotient: `x²` and `x` denote the same Boolean function but remain different formal
polynomials for Schwartz–Zippel evaluation.

When a CCX product would exceed `2^32`, the target is replaced by a fresh degree-one variable.
That is a conservative abstraction of the target's whole current Boolean function and exactly
matches the executable fallback.
-/

namespace TzapLean

noncomputable section

structure TrackedPolynomial where
  polynomial : BoolPolynomial
  degree : Nat

namespace TrackedPolynomial

def zero : TrackedPolynomial := ⟨0, 0⟩
def one : TrackedPolynomial := ⟨1, 0⟩
def fresh (i : Nat) : TrackedPolynomial := ⟨BoolPolynomial.var i, 1⟩

def add (x y : TrackedPolynomial) : TrackedPolynomial :=
  ⟨x.polynomial + y.polynomial, max x.degree y.degree⟩

def mul? (x y : TrackedPolynomial) : Option TrackedPolynomial :=
  let degree := x.degree + y.degree
  if degree ≤ Fingerprint.maxDegree then
    some ⟨x.polynomial * y.polynomial, degree⟩
  else none

def ccx (target left right : TrackedPolynomial) (freshVariable : Nat) : TrackedPolynomial :=
  match left.mul? right with
  | some product => target.add product
  | none => fresh freshVariable

/-- The stored number is a sound upper bound on formal total degree. -/
def DegreeBounded (x : TrackedPolynomial) : Prop := x.polynomial.totalDegree ≤ x.degree

/-- Every retained formal degree stays within the Schwartz–Zippel budget. -/
def WithinCutoff (x : TrackedPolynomial) : Prop := x.degree ≤ Fingerprint.maxDegree

theorem degreeBounded_zero : zero.DegreeBounded := by simp [DegreeBounded, zero]
theorem degreeBounded_one : one.DegreeBounded := by simp [DegreeBounded, one]
theorem degreeBounded_fresh (i : Nat) : (fresh i).DegreeBounded := by
  simp [DegreeBounded, fresh, BoolPolynomial.totalDegree_var]

theorem withinCutoff_zero : zero.WithinCutoff := by norm_num [WithinCutoff, zero, Fingerprint.maxDegree]
theorem withinCutoff_one : one.WithinCutoff := by norm_num [WithinCutoff, one, Fingerprint.maxDegree]
theorem withinCutoff_fresh (i : Nat) : (fresh i).WithinCutoff := by
  norm_num [WithinCutoff, fresh, Fingerprint.maxDegree]

/-- Every retained formal polynomial has total degree at most the global cutoff. -/
theorem totalDegree_le_cutoff {x : TrackedPolynomial}
    (hdegree : x.DegreeBounded) (hcutoff : x.WithinCutoff) :
    x.polynomial.totalDegree ≤ Fingerprint.maxDegree :=
  hdegree.trans hcutoff

theorem withinCutoff_add {x y : TrackedPolynomial}
    (hx : x.WithinCutoff) (hy : y.WithinCutoff) : (x.add y).WithinCutoff := by
  exact max_le hx hy

theorem withinCutoff_mul {x y product : TrackedPolynomial}
    (h : x.mul? y = some product) : product.WithinCutoff := by
  simp only [mul?, WithinCutoff] at h ⊢
  split at h
  · cases h
    assumption
  · simp at h

theorem degreeBounded_add {x y : TrackedPolynomial}
    (hx : x.DegreeBounded) (hy : y.DegreeBounded) : (x.add y).DegreeBounded := by
  apply le_trans (BoolPolynomial.totalDegree_add_le x.polynomial y.polynomial)
  exact max_le_max hx hy

theorem degreeBounded_mul {x y product : TrackedPolynomial}
    (hx : x.DegreeBounded) (hy : y.DegreeBounded) (h : x.mul? y = some product) :
    product.DegreeBounded := by
  simp only [mul?] at h
  split at h
  · cases h
    apply le_trans (BoolPolynomial.totalDegree_mul_le x.polynomial y.polynomial)
    exact Nat.add_le_add hx hy
  · simp at h

theorem degreeBounded_ccx {target left right : TrackedPolynomial} (freshVariable : Nat)
    (ht : target.DegreeBounded) (hl : left.DegreeBounded) (hr : right.DegreeBounded) :
    (ccx target left right freshVariable).DegreeBounded := by
  simp only [ccx]
  cases h : left.mul? right with
  | some product => exact degreeBounded_add ht (degreeBounded_mul hl hr h)
  | none => exact degreeBounded_fresh _

theorem withinCutoff_ccx {target left right : TrackedPolynomial} (freshVariable : Nat)
    (ht : target.WithinCutoff) : (ccx target left right freshVariable).WithinCutoff := by
  simp only [ccx]
  cases h : left.mul? right with
  | some product => exact withinCutoff_add ht (withinCutoff_mul h)
  | none => exact withinCutoff_fresh _

@[simp] theorem evalB_zero (v : Nat → Bool) : BoolPolynomial.evalB v zero.polynomial = false := by
  simp [zero, BoolPolynomial.evalB, BoolPolynomial.evalF₂, unbit]

@[simp] theorem evalB_one (v : Nat → Bool) : BoolPolynomial.evalB v one.polynomial = true := by
  simp [one, BoolPolynomial.evalB, BoolPolynomial.evalF₂, unbit]

@[simp] theorem evalB_fresh (v : Nat → Bool) (i : Nat) :
    BoolPolynomial.evalB v (fresh i).polynomial = v i := by
  simp [fresh]

@[simp] theorem evalB_add (v : Nat → Bool) (x y : TrackedPolynomial) :
    BoolPolynomial.evalB v (x.add y).polynomial =
      (BoolPolynomial.evalB v x.polynomial != BoolPolynomial.evalB v y.polynomial) := by
  simp [add]

/-- Below the cutoff, symbolic CCX is exactly the Boolean Toffoli update. -/
theorem evalB_ccx_of_mul {target left right product : TrackedPolynomial}
    (freshVariable : Nat) (h : left.mul? right = some product) (v : Nat → Bool) :
    BoolPolynomial.evalB v (ccx target left right freshVariable).polynomial =
      (BoolPolynomial.evalB v target.polynomial !=
        (BoolPolynomial.evalB v left.polynomial && BoolPolynomial.evalB v right.polynomial)) := by
  rw [ccx, h]
  have hp : product.polynomial = left.polynomial * right.polynomial := by
    simp only [mul?] at h
    split at h
    · cases h; rfl
    · simp at h
  simp [hp]

end TrackedPolynomial

/-- One tracked polynomial per wire and the next opaque-variable index. -/
structure NonlinearAState where
  wires : List TrackedPolynomial
  fresh : Nat

namespace NonlinearAState

def wireOf (st : NonlinearAState) (q : Qubit) : TrackedPolynomial :=
  st.wires.getD q TrackedPolynomial.zero

def initial (n : Nat) : NonlinearAState where
  wires := (List.range n).map TrackedPolynomial.fresh
  fresh := n

/-- Symbolic transfer functions in lock-step with `NState.step`. -/
def step (st : NonlinearAState) (g : Gate) : NonlinearAState :=
  match g with
  | .x q => { st with wires := st.wires.set q ((st.wireOf q).add TrackedPolynomial.one) }
  | .cnot c t => { st with wires := st.wires.set t ((st.wireOf t).add (st.wireOf c)) }
  | .h q =>
      { wires := st.wires.set q (TrackedPolynomial.fresh st.fresh), fresh := st.fresh + 1 }
  | .ccx c₁ c₂ t =>
      match (st.wireOf c₁).mul? (st.wireOf c₂) with
      | some product => { st with wires := st.wires.set t ((st.wireOf t).add product) }
      | none =>
          { wires := st.wires.set t (TrackedPolynomial.fresh st.fresh), fresh := st.fresh + 1 }
  | .reset q => { st with wires := st.wires.set q TrackedPolynomial.zero }
  | _ => st

def steps (st : NonlinearAState) : List Gate → NonlinearAState
  | [] => st
  | g :: gs => steps (st.step g) gs

def DegreeBounded (st : NonlinearAState) : Prop :=
  ∀ q : Qubit, (st.wireOf q).DegreeBounded

def WithinCutoff (st : NonlinearAState) : Prop :=
  ∀ q : Qubit, (st.wireOf q).WithinCutoff

theorem degreeBounded_set {st : NonlinearAState} {q : Qubit} {value : TrackedPolynomial}
    (hst : st.DegreeBounded) (hvalue : value.DegreeBounded) (fresh : Nat) :
    (⟨st.wires.set q value, fresh⟩ : NonlinearAState).DegreeBounded := by
  intro r
  simp only [wireOf, List.getD_eq_getElem?_getD]
  by_cases hr : r = q
  · subst r
    by_cases hq : q < st.wires.length
    · rw [List.getElem?_set_self hq]
      exact hvalue
    · rw [List.getElem?_eq_none]
      · exact TrackedPolynomial.degreeBounded_zero
      · simpa using Nat.le_of_not_lt hq
  · rw [List.getElem?_set, if_neg (by simpa using Ne.symm hr)]
    exact hst r

theorem withinCutoff_set {st : NonlinearAState} {q : Qubit} {value : TrackedPolynomial}
    (hst : st.WithinCutoff) (hvalue : value.WithinCutoff) (fresh : Nat) :
    (⟨st.wires.set q value, fresh⟩ : NonlinearAState).WithinCutoff := by
  intro r
  simp only [wireOf, List.getD_eq_getElem?_getD]
  by_cases hr : r = q
  · subst r
    by_cases hq : q < st.wires.length
    · rw [List.getElem?_set_self hq]
      exact hvalue
    · rw [List.getElem?_eq_none]
      · exact TrackedPolynomial.withinCutoff_zero
      · simpa using Nat.le_of_not_lt hq
  · rw [List.getElem?_set, if_neg (by simpa using Ne.symm hr)]
    exact hst r

theorem degreeBounded_initial (n : Nat) : (initial n).DegreeBounded := by
  intro q
  simp only [wireOf, initial, getD_map_range]
  split
  · exact TrackedPolynomial.degreeBounded_fresh _
  · exact TrackedPolynomial.degreeBounded_zero

theorem withinCutoff_initial (n : Nat) : (initial n).WithinCutoff := by
  intro q
  simp only [wireOf, initial, getD_map_range]
  split
  · exact TrackedPolynomial.withinCutoff_fresh _
  · exact TrackedPolynomial.withinCutoff_zero

theorem degreeBounded_step {st : NonlinearAState} (hst : st.DegreeBounded) (g : Gate) :
    (st.step g).DegreeBounded := by
  cases g with
  | x q =>
      exact degreeBounded_set hst
        (TrackedPolynomial.degreeBounded_add (hst q) TrackedPolynomial.degreeBounded_one) _
  | cnot c t =>
      exact degreeBounded_set hst
        (TrackedPolynomial.degreeBounded_add (hst t) (hst c)) _
  | h q =>
      exact degreeBounded_set hst (TrackedPolynomial.degreeBounded_fresh _) _
  | ccx c₁ c₂ t =>
      simp only [step]
      cases h : (st.wireOf c₁).mul? (st.wireOf c₂) with
      | some product =>
          exact degreeBounded_set hst (TrackedPolynomial.degreeBounded_add (hst t)
            (TrackedPolynomial.degreeBounded_mul (hst c₁) (hst c₂) h)) _
      | none => exact degreeBounded_set hst (TrackedPolynomial.degreeBounded_fresh _) _
  | reset q => exact degreeBounded_set hst TrackedPolynomial.degreeBounded_zero _
  | _ => exact hst

theorem withinCutoff_step {st : NonlinearAState} (hst : st.WithinCutoff) (g : Gate) :
    (st.step g).WithinCutoff := by
  cases g with
  | x q =>
      exact withinCutoff_set hst
        (TrackedPolynomial.withinCutoff_add (hst q) TrackedPolynomial.withinCutoff_one) _
  | cnot c t =>
      exact withinCutoff_set hst
        (TrackedPolynomial.withinCutoff_add (hst t) (hst c)) _
  | h q => exact withinCutoff_set hst (TrackedPolynomial.withinCutoff_fresh _) _
  | ccx c₁ c₂ t =>
      simp only [step]
      cases h : (st.wireOf c₁).mul? (st.wireOf c₂) with
      | some product =>
          exact withinCutoff_set hst (TrackedPolynomial.withinCutoff_add (hst t)
            (TrackedPolynomial.withinCutoff_mul h)) _
      | none => exact withinCutoff_set hst (TrackedPolynomial.withinCutoff_fresh _) _
  | reset q => exact withinCutoff_set hst TrackedPolynomial.withinCutoff_zero _
  | _ => exact hst

theorem degreeBounded_steps {st : NonlinearAState} (hst : st.DegreeBounded) :
    ∀ gates, (st.steps gates).DegreeBounded := by
  intro gates
  induction gates generalizing st with
  | nil => exact hst
  | cons g gates ih => exact ih (degreeBounded_step hst g)

theorem withinCutoff_steps {st : NonlinearAState} (hst : st.WithinCutoff) :
    ∀ gates, (st.steps gates).WithinCutoff := by
  intro gates
  induction gates generalizing st with
  | nil => exact hst
  | cons g gates ih => exact ih (withinCutoff_step hst g)

theorem length_step (st : NonlinearAState) (g : Gate) :
    (st.step g).wires.length = st.wires.length := by
  cases g <;> simp [step]
  split <;> simp

/-- Both sides make the cutoff decision from the same degree metadata. -/
theorem ccx_cutoff_agrees (st : NonlinearAState) (runtime : NState)
    (c₁ c₂ : Qubit)
    (h₁ : (st.wireOf c₁).degree = (runtime.fpOf c₁).degree)
    (h₂ : (st.wireOf c₂).degree = (runtime.fpOf c₂).degree) :
    ((st.wireOf c₁).mul? (st.wireOf c₂)).isSome =
      ((runtime.fpOf c₁).mul? (runtime.fpOf c₂)).isSome := by
  simp only [TrackedPolynomial.mul?, Fingerprint.mul?, Fingerprint.productDegree]
  rw [h₁, h₂]
  by_cases h : (runtime.fpOf c₁).degree + (runtime.fpOf c₂).degree ≤ Fingerprint.maxDegree
  · simp [h]
  · simp [h]

end NonlinearAState

end

end TzapLean
