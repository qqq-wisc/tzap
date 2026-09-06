import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Finsupp.Basic
import Mathlib.Algebra.CharP.Two
import Mathlib.Algebra.BigOperators.Finsupp.Basic
import Mathlib.Probability.Distributions.Uniform
import Mathlib.GroupTheory.Index
import Mathlib.Tactic.LinearCombination

/-!
# Affine Parities and Their Hashes

`PhaseFoldRand` tracks each wire's parity — the XOR of computational-basis values it
currently holds — not symbolically but as a random `k`-bit tag. This file is the algebra
and the probability behind that: what a parity *is*, what hashing one means, and how often
two different parities hash the same.

* `Form` — a Boolean affine function in normal form: a constant in `𝔽₂` plus a finitely
  supported coefficient vector. Equality of `Form`s is equality of functions, and it is
  decidable, which is what makes an exact parity comparison executable.
* `Sample m k` — one uniform `k`-bit string per variable, the finite probability space.
* `output p sample` — the hash of `p`: the 𝔽₂-sum of the draws of its variables, plus its
  constant. This is exactly the Rust tag arithmetic (`^` for XOR, `!` for negation).
* `affine_collision_bound` — **distinct parities collide with probability at most `2⁻ᵏ`**,
  by writing collision as membership in one fiber of a surjective 𝔽₂-linear map.
* `collides_probability_le` — the union bound over the parities a run compares.

This module combines the affine-form, collision, and finite-probability developments directly
on `Form`; the analysis does not need a separate syntactic parity-expression type.
-/

namespace TzapLean

noncomputable section

open scoped ENNReal

/-! ## `𝔽₂` -/

/-- The two-element field, the coefficient field for affine parities. -/
abbrev F₂ := ZMod 2

/-- Embed a Boolean into `𝔽₂`, turning xor into addition. -/
def bit (b : Bool) : F₂ := if b then 1 else 0

theorem bit_xor (a b : Bool) : bit (a != b) = bit a + bit b := by
  cases a <;> cases b <;> simp [bit]
  exact (CharTwo.add_self_eq_zero (1 : F₂)).symm

/-- The inverse of `bit`: `𝔽₂` back to `Bool`. Arithmetic in `ZMod 2` goes through a ring
instance that cannot be specialized, which costs ~600× what the corresponding `Bool`
operation costs — measured. The algorithm therefore stores `Bool`s and uses these to talk to
the theory. -/
def unbit (x : F₂) : Bool := decide (x = 1)

@[simp] theorem bit_unbit (x : F₂) : bit (unbit x) = x := by revert x; decide

@[simp] theorem unbit_bit (b : Bool) : unbit (bit b) = b := by cases b <;> simp [unbit, bit]

/-! ## Affine forms -/

/-- A Boolean affine function in normal form: `c ⊕ ⨁_{i ∈ support} vᵢ`. Since coefficients
live in `𝔽₂` the representation is canonical — two forms are equal as structures iff they
denote the same function. -/
structure Form where
  /-- The constant term. -/
  constant : F₂
  /-- One coefficient per variable; `1` exactly when the variable occurs. -/
  coefficients : Nat →₀ F₂
deriving DecidableEq

namespace Form

@[ext] theorem ext {p q : Form} (hc : p.constant = q.constant)
    (hl : p.coefficients = q.coefficients) : p = q := by
  cases p; cases q; simp_all

instance : Zero Form := ⟨⟨0, 0⟩⟩
instance : Add Form := ⟨fun p q => ⟨p.constant + q.constant, p.coefficients + q.coefficients⟩⟩
instance : Sub Form := ⟨fun p q => ⟨p.constant - q.constant, p.coefficients - q.coefficients⟩⟩

@[simp] theorem add_constant (p q : Form) : (p + q).constant = p.constant + q.constant := rfl
@[simp] theorem add_coefficients (p q : Form) :
    (p + q).coefficients = p.coefficients + q.coefficients := rfl
@[simp] theorem sub_constant (p q : Form) :
    (p - q).constant = p.constant - q.constant := rfl
@[simp] theorem sub_coefficients (p q : Form) :
    (p - q).coefficients = p.coefficients - q.coefficients := rfl

/-- The constant parity. -/
def const (b : Bool) : Form := ⟨bit b, 0⟩

/-- The parity consisting of the single variable `i`. -/
def var (i : Nat) : Form := ⟨0, Finsupp.single i 1⟩

/-- Complementing a parity (the `X` gate, and Rust's `!h`). -/
def flip (p : Form) : Form := p + const true

/-- Evaluate a form at a valuation of the variables. -/
def eval (valuation : Nat → F₂) (p : Form) : F₂ :=
  p.constant + p.coefficients.sum (fun i coefficient => coefficient * valuation i)

@[simp] theorem eval_add (v) (p q : Form) : eval v (p + q) = eval v p + eval v q := by
  simp [eval, add_assoc, add_left_comm, add_comm, add_mul, Finsupp.sum_add_index]

@[simp] theorem eval_const (v) (b : Bool) : eval v (const b) = bit b := by
  simp [eval, const]

@[simp] theorem eval_var (v) (i : Nat) : eval v (var i) = v i := by
  simp [eval, var]

@[simp] theorem eval_flip (v) (p : Form) : eval v p.flip = eval v p + 1 := by
  simp [flip, eval_add, bit]

/-- Evaluate at a Boolean valuation, back in `Bool`. -/
def evalB (valuation : Nat → Bool) (p : Form) : Bool :=
  eval (fun i => bit (valuation i)) p = 1

theorem bit_eq_one_iff (b : Bool) : bit b = 1 ↔ b = true := by
  cases b <;> simp [bit]

@[simp] theorem evalB_const (v : Nat → Bool) (b : Bool) : evalB v (const b) = b := by
  simp [evalB, bit_eq_one_iff]

@[simp] theorem evalB_var (v : Nat → Bool) (i : Nat) : evalB v (var i) = v i := by
  simp [evalB, bit_eq_one_iff]

@[simp] theorem evalB_add (v : Nat → Bool) (p q : Form) :
    evalB v (p + q) = (evalB v p != evalB v q) := by
  simp only [evalB, eval_add]
  generalize eval (fun i => bit (v i)) p = a
  generalize eval (fun i => bit (v i)) q = b
  revert a b
  decide

@[simp] theorem evalB_flip (v : Nat → Bool) (p : Form) : evalB v p.flip = !evalB v p := by
  simp [flip]

/-- A form is `Bounded m` when every variable it mentions is below `m`. -/
def Bounded (m : Nat) (p : Form) : Prop := ∀ i ∈ p.coefficients.support, i < m

theorem bounded_const (m : Nat) (b : Bool) : Bounded m (const b) := by
  intro i hi; simp [const] at hi

theorem bounded_var {m i : Nat} (h : i < m) : Bounded m (var i) := by
  intro j hj
  simp only [var, Finsupp.mem_support_iff, Finsupp.single_apply] at hj
  split at hj
  · omega
  · exact absurd rfl hj

theorem bounded_add {m : Nat} {p q : Form} (hp : Bounded m p) (hq : Bounded m q) :
    Bounded m (p + q) := by
  intro i hi
  simp only [add_coefficients, Finsupp.mem_support_iff, Finsupp.add_apply] at hi
  by_cases hpi : p.coefficients i = 0
  · have hqi : q.coefficients i ≠ 0 := by
      intro hq0; simp [hpi, hq0] at hi
    exact hq i (Finsupp.mem_support_iff.mpr hqi)
  · exact hp i (Finsupp.mem_support_iff.mpr hpi)

theorem bounded_flip {m : Nat} {p : Form} (hp : Bounded m p) : Bounded m p.flip :=
  bounded_add hp (bounded_const m true)

theorem bounded_mono {m m' : Nat} {p : Form} (h : m ≤ m') (hp : Bounded m p) :
    Bounded m' p := fun i hi => lt_of_lt_of_le (hp i hi) h

/-- Only the variables a form mentions matter: valuations agreeing below its bound agree on it. -/
theorem evalB_congr {m : Nat} {p : Form} (hp : Bounded m p) {v v' : Nat → Bool}
    (h : ∀ i, i < m → v' i = v i) : evalB v' p = evalB v p := by
  have heval : eval (fun i => bit (v' i)) p = eval (fun i => bit (v i)) p := by
    simp only [eval]
    congr 1
    exact Finsupp.sum_congr fun i hi => by rw [h i (hp i hi)]
  simp [evalB, heval]

theorem sub_bounded {m : Nat} {p q : Form} (hp : Bounded m p) (hq : Bounded m q) :
    Bounded m (p - q) := by
  intro i hi
  simp only [sub_coefficients, Finsupp.mem_support_iff, Finsupp.sub_apply] at hi ⊢
  by_cases hpi : p.coefficients i = 0
  · have hqi : q.coefficients i ≠ 0 := by
      intro hq0; simp [hpi, hq0] at hi
    exact hq i (Finsupp.mem_support_iff.mpr hqi)
  · exact hp i (Finsupp.mem_support_iff.mpr hpi)

end Form

/-! ## Hashing -/

/-- A `k`-bit string; XOR is pointwise addition. -/
abbrev BitString (k : Nat) := Fin k → F₂

/-- One `k`-bit draw per variable. -/
abbrev Draws (k : Nat) := Nat → BitString k

/-- The finite sample space: one draw per variable below `m`. -/
abbrev Sample (m k : Nat) := Fin m → BitString k

/-- Pad a finite sample out to a full stream. Bounded forms never look past `m`. -/
def liftSample {m k : Nat} (sample : Sample m k) : Draws k :=
  fun i => if h : i < m then sample ⟨i, h⟩ else 0

/-- The hash of a parity under a draw stream: the 𝔽₂-sum of its variables' draws plus its
constant — Rust's tag arithmetic. -/
def hash {k : Nat} (draws : Draws k) (p : Form) : BitString k :=
  fun j => Form.eval (fun i => draws i j) p

/-- The hash under a finite sample. -/
def output {m k : Nat} (p : Form) (sample : Sample m k) : BitString k :=
  hash (liftSample sample) p

@[simp] theorem hash_add {k} (draws : Draws k) (p q : Form) :
    hash draws (p + q) = hash draws p + hash draws q := by
  funext j; simp [hash]

@[simp] theorem hash_var {k} (draws : Draws k) (i : Nat) : hash draws (Form.var i) = draws i := by
  funext j; simp [hash]

@[simp] theorem hash_const {k} (draws : Draws k) (b : Bool) :
    hash draws (Form.const b) = fun _ => bit b := by
  funext j; simp [hash]

/-! ## The collision bound

Ported from `Collision.lean` and `FiniteProbability.lean`. -/

/-- Every fiber of a surjective homomorphism between finite additive groups has probability
`|B|⁻¹` under the uniform PMF on the domain. -/
theorem uniform_fiber_of_surjective {A : Type u} {B : Type v} [AddGroup A] [Fintype A]
    [AddGroup B] [Fintype B] [DecidableEq B] (f : A →+ B) (hf : Function.Surjective f) (y : B) :
    (PMF.uniformOfFintype A).toOuterMeasure {x | f x = y} = (Fintype.card B : ℝ≥0∞)⁻¹ := by
  classical
  let fiberCard := (Finset.univ.filter fun x : A => f x = y).card
  have hfib (b : B) : (Finset.univ.filter fun x : A => f x = b).card = fiberCard :=
    AddMonoidHom.card_fiber_eq_of_mem_range f (hf b) (hf y)
  have hcard : Fintype.card A = Fintype.card B * fiberCard := by
    calc
      Fintype.card A = ∑ b ∈ (Finset.univ : Finset B),
          (Finset.univ.filter fun x : A => f x = b).card := by
            rw [← Finset.card_univ]
            exact Finset.card_eq_sum_card_fiberwise (fun _ _ => Finset.mem_univ _)
      _ = ∑ _b ∈ (Finset.univ : Finset B), fiberCard :=
            Finset.sum_congr rfl fun b _ => hfib b
      _ = Fintype.card B * fiberCard := by simp
  have hfiber_pos : 0 < fiberCard := by
    rcases hf y with ⟨x, hx⟩
    exact Finset.card_pos.mpr ⟨x, by simp [hx]⟩
  rw [PMF.toOuterMeasure_uniformOfFintype_apply, Fintype.card_subtype]
  change (fiberCard : ℝ≥0∞) / Fintype.card A = (Fintype.card B : ℝ≥0∞)⁻¹
  rw [hcard, Nat.cast_mul, ENNReal.div_eq_inv_mul]
  have hBtop : (Fintype.card B : ℝ≥0∞) ≠ ⊤ := ENNReal.natCast_ne_top _
  have hFtop : (fiberCard : ℝ≥0∞) ≠ ⊤ := ENNReal.natCast_ne_top _
  have hFzero : (fiberCard : ℝ≥0∞) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hfiber_pos
  rw [ENNReal.mul_inv (Or.inr hFtop) (Or.inl hBtop), mul_assoc,
    ENNReal.inv_mul_cancel hFzero hFtop, mul_one]

theorem sum_liftSample_eq {m k : Nat} (p : Form) (hp : Form.Bounded m p)
    (sample : Sample m k) (j : Fin k) :
    p.coefficients.sum (fun i coefficient => coefficient * liftSample sample i j) =
      ∑ i : Fin m, p.coefficients i.val * sample i j := by
  rw [Finsupp.sum_of_support_subset p.coefficients
    (s := Finset.range m) (fun i hi => Finset.mem_range.mpr (hp i hi))]
  · calc
      (∑ i ∈ Finset.range m, p.coefficients i * liftSample sample i j) =
          ∑ i : Fin m, p.coefficients i.val * liftSample sample i.val j :=
            (Fin.sum_univ_eq_sum_range
              (fun i => p.coefficients i * liftSample sample i j) m).symm
      _ = ∑ i : Fin m, p.coefficients i.val * sample i j := by
            refine Finset.sum_congr rfl fun i _ => ?_
            simp [liftSample, i.isLt]
  · intro i _; simp

/-- The linear part of hashing, as an additive homomorphism. -/
def linearMap (coefficients : Nat →₀ F₂) (m k : Nat) : Sample m k →+ BitString k where
  toFun sample := fun j => ∑ i : Fin m, coefficients i.val * sample i j
  map_zero' := by funext j; simp
  map_add' left right := by funext j; simp [mul_add, Finset.sum_add_distrib]

theorem output_eq_constant_add_linear {m k : Nat} (p : Form) (hp : Form.Bounded m p)
    (sample : Sample m k) :
    output p sample = fun j => p.constant + linearMap p.coefficients m k sample j := by
  funext j
  change p.constant +
      p.coefficients.sum (fun i coefficient => coefficient * liftSample sample i j) =
    p.constant + ∑ i : Fin m, p.coefficients i.val * sample i j
  rw [sum_liftSample_eq p hp sample j]

theorem linearMap_surjective_of_coeff_ne_zero {m k : Nat} (coefficients : Nat →₀ F₂)
    {pivot : Fin m} (hpivot : coefficients pivot.val ≠ 0) :
    Function.Surjective (linearMap coefficients m k) := by
  intro target
  let sample : Sample m k := fun i j =>
    if i = pivot then (coefficients pivot.val)⁻¹ * target j else 0
  refine ⟨sample, ?_⟩
  funext j
  change (∑ i : Fin m, coefficients i.val * sample i j) = target j
  rw [Finset.sum_eq_single pivot]
  · simp only [sample, if_pos]
    have hone : coefficients pivot.val = (1 : F₂) := by
      apply (ZMod.val_eq_one (by norm_num) _).mp
      have hvne : (coefficients pivot.val).val ≠ 0 := by
        intro hv
        exact hpivot ((ZMod.val_eq_zero _).mp hv)
      have hvlt := ZMod.val_lt (coefficients pivot.val)
      omega
    simp [hone]
  · intro i _ hne; simp [sample, hne]
  · simp

theorem linearMap_sub (p q : Form) (m k : Nat) (sample : Sample m k) :
    linearMap (p - q).coefficients m k sample =
      linearMap p.coefficients m k sample - linearMap q.coefficients m k sample := by
  funext j
  simp [linearMap, Finset.sum_sub_distrib, sub_mul]

theorem output_eq_iff_linear_eq {m k : Nat} {p q : Form}
    (hp : Form.Bounded m p) (hq : Form.Bounded m q) (sample : Sample m k) :
    output p sample = output q sample ↔
      linearMap (p - q).coefficients m k sample = fun _ => -(p - q).constant := by
  rw [output_eq_constant_add_linear p hp, output_eq_constant_add_linear q hq]
  constructor
  · intro h
    funext j
    have hj := congrFun h j
    change p.constant + linearMap p.coefficients m k sample j =
      q.constant + linearMap q.coefficients m k sample j at hj
    rw [linearMap_sub]
    change linearMap p.coefficients m k sample j -
      linearMap q.coefficients m k sample j = -(p.constant - q.constant)
    linear_combination hj
  · intro h
    funext j
    have hj := congrFun h j
    rw [linearMap_sub] at hj
    change linearMap p.coefficients m k sample j -
      linearMap q.coefficients m k sample j = -(p.constant - q.constant) at hj
    change p.constant + linearMap p.coefficients m k sample j =
      q.constant + linearMap q.coefficients m k sample j
    linear_combination hj

theorem card_bitString (k : Nat) : Fintype.card (BitString k) = 2 ^ k := by
  simp [BitString, ZMod.card]

theorem inv_card_bitString (k : Nat) :
    (Fintype.card (BitString k) : ℝ≥0∞)⁻¹ = ((2 : ℝ≥0∞)⁻¹) ^ k := by
  rw [card_bitString]
  push_cast
  exact ENNReal.inv_pow

/-- **Two distinct parities collide with probability at most `2⁻ᵏ`.** -/
theorem affine_collision_bound {m k : Nat} (p q : Form)
    (hp : Form.Bounded m p) (hq : Form.Bounded m q) (hne : p ≠ q) :
    (PMF.uniformOfFintype (Sample m k)).toOuterMeasure
        {sample | output p sample = output q sample} ≤ ((2 : ℝ≥0∞)⁻¹) ^ k := by
  by_cases hk : k = 0
  · subst k
    calc
      (PMF.uniformOfFintype (Sample m 0)).toOuterMeasure
          {sample | output p sample = output q sample} ≤
          (PMF.uniformOfFintype (Sample m 0)).toOuterMeasure Set.univ :=
        (PMF.uniformOfFintype (Sample m 0)).toOuterMeasure.mono (Set.subset_univ _)
      _ = 1 := (PMF.toOuterMeasure_apply_eq_one_iff _ _).2 (Set.subset_univ _)
      _ = ((2 : ℝ≥0∞)⁻¹) ^ 0 := by simp
  let d := p - q
  have hdBounded : Form.Bounded m d := Form.sub_bounded hp hq
  by_cases hlinear : d.coefficients = 0
  · have hconstant : d.constant ≠ 0 := by
      intro hc
      refine hne (Form.ext ?_ ?_)
      · exact sub_eq_zero.mp (by simpa [d] using hc)
      · exact sub_eq_zero.mp (by simpa [d] using hlinear)
    have hno (sample : Sample m k) : output p sample ≠ output q sample := by
      intro heq
      have hlin := (output_eq_iff_linear_eq hp hq sample).mp heq
      have j : Fin k := ⟨0, Nat.pos_of_ne_zero hk⟩
      have hj := congrFun hlin j
      have hz : linearMap d.coefficients m k sample j = 0 := by
        simp [hlinear, linearMap]
      rw [hz] at hj
      exact hconstant (neg_eq_zero.mp hj.symm)
    have hevent : {sample : Sample m k | output p sample = output q sample} = ∅ := by
      ext sample; simp [hno sample]
    rw [hevent]; simp
  · have hsupp : d.coefficients.support.Nonempty := Finsupp.support_nonempty_iff.mpr hlinear
    obtain ⟨i, hiSupport⟩ := hsupp
    have hi : d.coefficients i ≠ 0 := Finsupp.mem_support_iff.mp hiSupport
    let pivot : Fin m := ⟨i, hdBounded i hiSupport⟩
    have hpivot : d.coefficients pivot.val ≠ 0 := hi
    have hsurj : Function.Surjective (linearMap d.coefficients m k) :=
      linearMap_surjective_of_coeff_ne_zero d.coefficients hpivot
    have hevent : {sample : Sample m k | output p sample = output q sample} =
        {sample | linearMap d.coefficients m k sample = fun _ => -d.constant} := by
      ext sample
      simpa [d] using output_eq_iff_linear_eq hp hq sample
    rw [hevent, uniform_fiber_of_surjective _ hsurj, inv_card_bitString]

/-! ## The union bound over one run's comparisons -/

/-- The bad event: two of the compared parities are distinct but hash the same. -/
def Collides {m k : Nat} (ps : List Form) (sample : Sample m k) : Prop :=
  ∃ p ∈ ps, ∃ q ∈ ps, p ≠ q ∧ output p sample = output q sample

/-- **The union bound**: a run comparing `t` parities is unfaithful with probability at most
`C(t,2)·2⁻ᵏ`. -/
theorem collides_probability_le {m k : Nat} (ps : List Form)
    (hps : ∀ p ∈ ps, Form.Bounded m p) :
    (PMF.uniformOfFintype (Sample m k)).toOuterMeasure {sample | Collides ps sample} ≤
      (ps.length.choose 2 : ℝ≥0∞) * ((2 : ℝ≥0∞)⁻¹) ^ k := by
  classical
  let μ := (PMF.uniformOfFintype (Sample m k)).toOuterMeasure
  let get : Nat → Form := fun i => ps.getD i (Form.const false)
  let E : Nat → Nat → Set (Sample m k) := fun i j =>
    {sample | get i ≠ get j ∧ output (get i) sample = output (get j) sample}
  have hgetD : ∀ i (h : i < ps.length), get i = ps[i] := fun i h => by
    simp [get, List.getD_eq_getElem?_getD, List.getElem?_eq_getElem h]
  have hsub : {sample : Sample m k | Collides ps sample} ⊆
      ⋃ i ∈ Finset.range ps.length, ⋃ j ∈ Finset.range i, E i j := by
    rintro sample ⟨p, hp, p', hp', hne, hout⟩
    rcases List.mem_iff_getElem.mp hp with ⟨i, hi, hpi⟩
    rcases List.mem_iff_getElem.mp hp' with ⟨j, hj, hpj⟩
    have hgi : get i = p := (hgetD i hi).trans hpi
    have hgj : get j = p' := (hgetD j hj).trans hpj
    have hij : i ≠ j := by
      rintro rfl
      exact hne (hpi.symm.trans hpj)
    rcases Nat.lt_or_ge i j with hlt | hge
    · refine Set.mem_biUnion (Finset.mem_range.mpr hj)
        (Set.mem_biUnion (Finset.mem_range.mpr hlt) ?_)
      exact ⟨by rw [hgi, hgj]; exact hne.symm, by rw [hgi, hgj]; exact hout.symm⟩
    · refine Set.mem_biUnion (Finset.mem_range.mpr hi)
        (Set.mem_biUnion (Finset.mem_range.mpr (lt_of_le_of_ne hge hij.symm)) ?_)
      exact ⟨by rw [hgi, hgj]; exact hne, by rw [hgi, hgj]; exact hout⟩
  have hpair : ∀ i < ps.length, ∀ j < i, μ (E i j) ≤ ((2 : ℝ≥0∞)⁻¹) ^ k := by
    intro i hi j hj
    have hj' : j < ps.length := hj.trans hi
    have hmem : ∀ {a : Nat}, a < ps.length → get a ∈ ps := fun {a} ha =>
      (hgetD a ha) ▸ List.getElem_mem ha
    by_cases heq : get i = get j
    · have hempty : E i j = ∅ := by ext sample; simp [E, heq]
      rw [hempty]; simp
    · calc
        μ (E i j) ≤ μ {sample : Sample m k |
            output (get i) sample = output (get j) sample} := μ.mono fun sample hs => hs.2
        _ ≤ ((2 : ℝ≥0∞)⁻¹) ^ k :=
          affine_collision_bound _ _ (hps _ (hmem hi)) (hps _ (hmem hj')) heq
  calc
    (PMF.uniformOfFintype (Sample m k)).toOuterMeasure {sample | Collides ps sample} ≤
        μ (⋃ i ∈ Finset.range ps.length, ⋃ j ∈ Finset.range i, E i j) := μ.mono hsub
    _ ≤ ∑ i ∈ Finset.range ps.length, μ (⋃ j ∈ Finset.range i, E i j) :=
      MeasureTheory.measure_biUnion_finset_le _ _
    _ ≤ ∑ i ∈ Finset.range ps.length, ∑ j ∈ Finset.range i, μ (E i j) :=
      Finset.sum_le_sum fun i _ => MeasureTheory.measure_biUnion_finset_le _ _
    _ ≤ ∑ i ∈ Finset.range ps.length, ∑ _j ∈ Finset.range i, ((2 : ℝ≥0∞)⁻¹) ^ k :=
      Finset.sum_le_sum fun i hi => Finset.sum_le_sum fun j hj =>
        hpair i (Finset.mem_range.mp hi) j (Finset.mem_range.mp hj)
    _ = ((∑ i ∈ Finset.range ps.length, i : ℕ) : ℝ≥0∞) * ((2 : ℝ≥0∞)⁻¹) ^ k := by
      simp only [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
      rw [← Finset.sum_mul]
      push_cast
      rfl
    _ = (ps.length.choose 2 : ℝ≥0∞) * ((2 : ℝ≥0∞)⁻¹) ^ k := by
      rw [Finset.sum_range_id, ← Nat.choose_two_right]

end
end TzapLean
