import TzapLean.Support

/-!
# Gate Matrix Algebra

Products of gate matrices, in the form the rewrites of an optimization pass need:

* `embed1_mul` — two one-wire gates on the same wire multiply as their `2 × 2` matrices;
* `permMatrix_mul` / `phaseMatrix_mul` — permutations compose, diagonals multiply
  pointwise (so all diagonal gates commute);
* `gateUnitary_mul_self` — the seven self-inverse gates square to the identity, which is
  what licenses `CancelGates`' first sweep.

`Gate.Wf` records the side condition the last of these needs: a multi-qubit gate's operands
must be distinct. `cnot q q` is not a gate QASM can express and the Rust implementation
never builds one; semantically it would be the map `b ↦ b[q := 0]`, which is idempotent
rather than self-inverse, so cancelling a pair of them would be unsound.
-/

namespace TzapLean

open Matrix

noncomputable section

/-- The gates `CancelGates` treats as self-inverse: `h`, `x`, `z`, `cnot`, `cz`, `ccx`,
`ccz` (the Rust `is_self_inverse`). -/
def Gate.isSelfInverse : Gate → Bool
  | .h _ | .x _ | .z _ | .cnot .. | .cz .. | .ccx .. | .ccz .. => true
  | _ => false

/-! ## One-wire gates -/

/-- The `2 × 2` matrix product, in the `M output input` convention. -/
def mul2 (M M' : Bool → Bool → ℂ) : Bool → Bool → ℂ :=
  fun out inp => ∑ k : Bool, M out k * M' k inp

/-- The `2 × 2` identity. -/
def id2 : Bool → Bool → ℂ := fun out inp => if out = inp then 1 else 0

@[simp] theorem embed1_id2 (n : Nat) (q : Qubit) : embed1 n id2 q = 1 := by
  by_cases h : q < n
  · funext out inp
    rw [embed1_apply_of_lt _ h]
    simp only [Matrix.one_apply, id2]
    by_cases hall : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
    · rw [if_pos hall]
      by_cases hq : out ⟨q, h⟩ = inp ⟨q, h⟩
      · have heq : out = inp := by
          funext r
          by_cases hr : (r : Nat) = q
          · have : r = ⟨q, h⟩ := Fin.ext hr
            rw [this]; exact hq
          · exact hall r hr
        simp [heq]
      · have : out ≠ inp := fun hc => hq (by rw [hc])
        simp [hq, this]
    · rw [if_neg hall]
      have : out ≠ inp := fun hc => hall fun r _ => by rw [hc]
      simp [this]
  · simp [embed1_eq_one _ h]

/-- Two gates on the same wire multiply as their `2 × 2` matrices. -/
theorem embed1_mul (n : Nat) (M M' : Bool → Bool → ℂ) (q : Qubit) :
    embed1 n M q * embed1 n M' q = embed1 n (mul2 M M') q := by
  by_cases h : q < n
  · funext out inp
    rw [Matrix.mul_apply]
    simp only [embed1_apply_of_lt _ h]
    set kf : Basis n := out.set q false with hkf
    set kt : Basis n := out.set q true with hkt
    have hkf_off : ∀ r : Fin n, (r : Nat) ≠ q → kf r = out r := by
      intro r hr; simp [hkf, Basis.set, hr]
    have hkt_off : ∀ r : Fin n, (r : Nat) ≠ q → kt r = out r := by
      intro r hr; simp [hkt, Basis.set, hr]
    have hkf_q : kf ⟨q, h⟩ = false := by simp [hkf, Basis.set]
    have hkt_q : kt ⟨q, h⟩ = true := by simp [hkt, Basis.set]
    have hne : kf ≠ kt := by
      intro hc
      have hb : kf ⟨q, h⟩ = kt ⟨q, h⟩ := by rw [hc]
      rw [hkf_q, hkt_q] at hb
      exact Bool.false_ne_true hb
    by_cases hoff : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
    · set vf : ℂ := M (out ⟨q, h⟩) false * M' false (inp ⟨q, h⟩) with hvf
      set vt : ℂ := M (out ⟨q, h⟩) true * M' true (inp ⟨q, h⟩) with hvt
      have key : ∀ k : Basis n,
          ((if ∀ r : Fin n, (r : Nat) ≠ q → out r = k r then
              M (out ⟨q, h⟩) (k ⟨q, h⟩) else 0) *
           (if ∀ r : Fin n, (r : Nat) ≠ q → k r = inp r then
              M' (k ⟨q, h⟩) (inp ⟨q, h⟩) else 0))
          = (if k = kf then vf else 0) + (if k = kt then vt else 0) := by
        intro k
        by_cases hk : ∀ r : Fin n, (r : Nat) ≠ q → out r = k r
        · have hk2 : ∀ r : Fin n, (r : Nat) ≠ q → k r = inp r := by
            intro r hr; rw [← hk r hr]; exact hoff r hr
          have hsplit : k = kf ∨ k = kt := by
            cases hkq : k ⟨q, h⟩
            · left
              funext r
              by_cases hr : (r : Nat) = q
              · have : r = ⟨q, h⟩ := Fin.ext hr
                rw [this, hkq, hkf_q]
              · rw [hkf_off r hr, ← hk r hr]
            · right
              funext r
              by_cases hr : (r : Nat) = q
              · have : r = ⟨q, h⟩ := Fin.ext hr
                rw [this, hkq, hkt_q]
              · rw [hkt_off r hr, ← hk r hr]
          rcases hsplit with hs | hs
          · subst hs
            rw [if_pos hk, if_pos hk2, if_pos rfl, if_neg hne, hkf_q]
            simp [hvf]
          · subst hs
            rw [if_pos hk, if_pos hk2, if_neg (Ne.symm hne), if_pos rfl, hkt_q]
            simp [hvt]
        · have hkf' : k ≠ kf := by
            intro hc; exact hk (fun r hr => by rw [hc, hkf_off r hr])
          have hkt' : k ≠ kt := by
            intro hc; exact hk (fun r hr => by rw [hc, hkt_off r hr])
          rw [if_neg hk, if_neg hkf', if_neg hkt']
          simp
      rw [Finset.sum_congr rfl fun k _ => key k, Finset.sum_add_distrib,
        Finset.sum_ite_eq' Finset.univ, Finset.sum_ite_eq' Finset.univ]
      simp only [Finset.mem_univ, if_pos]
      rw [if_pos hoff]
      simp [mul2, hvf, hvt, add_comm]
    · rw [if_neg hoff]
      refine Finset.sum_eq_zero fun k _ => ?_
      by_cases hk : ∀ r : Fin n, (r : Nat) ≠ q → out r = k r
      · have : ¬ ∀ r : Fin n, (r : Nat) ≠ q → k r = inp r := by
          intro hc
          exact hoff fun r hr => (hk r hr).trans (hc r hr)
        rw [if_neg this]
        ring
      · rw [if_neg hk]
        ring
  · simp [embed1_eq_one _ h]

/-- Scaling a one-wire matrix scales the embedded operator (in range). -/
theorem embed1_smul {n : Nat} (c : ℂ) (M : Bool → Bool → ℂ) (q : Qubit) (hq : q < n) :
    embed1 n (fun o i => c * M o i) q = c • embed1 n M q := by
  funext out inp
  simp only [embed1_apply_of_lt _ hq, Matrix.smul_apply, smul_eq_mul]
  by_cases h : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r <;> simp [h]

/-! ## Permutations and diagonals -/

theorem permMatrix_mul {n : Nat} (σ τ : Basis n → Basis n) :
    permMatrix σ * permMatrix τ = permMatrix (σ ∘ τ) := by
  funext out inp
  rw [Matrix.mul_apply]
  simp only [permMatrix, Function.comp_apply]
  rw [Finset.sum_eq_single (τ inp)]
  · simp
  · intro k _ hk; simp [hk]
  · simp

@[simp] theorem permMatrix_id {n : Nat} : permMatrix (fun b : Basis n => b) = 1 := by
  funext out inp
  simp [permMatrix, Matrix.one_apply]

theorem phaseMatrix_mul {n : Nat} (f g : Basis n → ℂ) :
    phaseMatrix f * phaseMatrix g = phaseMatrix (fun b => f b * g b) := by
  funext out inp
  rw [Matrix.mul_apply]
  simp only [phaseMatrix]
  rw [Finset.sum_eq_single inp]
  · by_cases h : out = inp <;> simp [h]
  · intro k _ hk; simp [hk]
  · simp

/-! ## Squares of the self-inverse gates -/

/-- Writing a wire twice keeps only the second value. -/
@[simp] theorem Basis.set_set {n : Nat} (b : Basis n) (q : Qubit) (v v' : Bool) :
    (b.set q v).set q v' = b.set q v' := by
  funext r
  by_cases hr : (r : Nat) = q <;> simp [Basis.set, hr]

/-- Writing a wire back the value it already holds changes nothing. -/
theorem Basis.set_get_self {n : Nat} (b : Basis n) (q : Qubit) : b.set q (b.get q) = b := by
  funext r
  by_cases hr : (r : Nat) = q
  · have hq : q < n := hr ▸ r.isLt
    have : r = ⟨q, hq⟩ := Fin.ext hr
    subst this
    simp [Basis.set, Basis.get, hq]
  · simp [Basis.set, hr]

@[simp] theorem mul2_diag2 (a b c d : ℂ) : mul2 (diag2 a b) (diag2 c d) = diag2 (a * c) (b * d) := by
  funext out inp
  cases out <;> cases inp <;> simp [mul2, diag2]

@[simp] theorem diag2_one_one : diag2 1 1 = id2 := by
  funext out inp
  cases out <;> cases inp <;> simp [diag2, id2]

@[simp] theorem mul2_x2 : mul2 x2 x2 = id2 := by
  funext out inp
  cases out <;> cases inp <;> simp [mul2, x2, id2]

theorem ofReal_sqrt_two_mul_self :
    ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := by
  rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by norm_num : (0:ℝ) ≤ 2)]
  norm_num

theorem ofReal_sqrt_two_ne_zero : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := by
  simp only [ne_eq, Complex.ofReal_eq_zero]
  positivity

@[simp] theorem mul2_h2 : mul2 h2 h2 = id2 := by
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) ^ 2 = 2 := by
    rw [pow_two]; exact ofReal_sqrt_two_mul_self
  have hne := ofReal_sqrt_two_ne_zero
  funext out inp
  cases out <;> cases inp <;>
    · simp only [mul2, h2, id2, Fintype.sum_bool]
      norm_num
      field_simp
      rw [hsq]
      ring

/-- `e^{2πi} = 1`. -/
@[simp] theorem ep_two : ep 2 = 1 := by
  simp only [ep]
  rw [show ((Real.pi * ((2 : ℚ) : ℝ) : ℝ) : ℂ) * Complex.I = 2 * (Real.pi : ℂ) * Complex.I by
    push_cast; ring]
  exact Complex.exp_two_pi_mul_I

/-- `e^{iπθ}` has period 2 in `θ`. -/
@[simp] theorem ep_add_two (θ : ℚ) : ep (θ + 2) = ep θ := by
  rw [ep_add, ep_two, mul_one]

/-- The seven gates the pass treats as self-inverse really are involutions.
The `Wf` hypothesis rules out `cnot q q` and `ccx` with a repeated operand. -/
theorem gateUnitary_mul_self (n : Nat) (g : Gate) (hwf : g.Wf) (hs : g.isSelfInverse = true) :
    gateUnitary n g * gateUnitary n g = 1 := by
  cases g with
  | h q => simp [gateUnitary, embed1_mul]
  | x q => simp [gateUnitary, embed1_mul]
  | z q =>
      have : mul2 (diag2 1 (ep 1)) (diag2 1 (ep 1)) = id2 := by
        rw [mul2_diag2, ← ep_add]
        norm_num
      simp [gateUnitary, embed1_mul, this]
  | cnot c tgt =>
      have hct : c ≠ tgt := hwf
      rw [gateUnitary, permMatrix_mul]
      have : (fun b : Basis n => b.set tgt (b.get tgt != b.get c)) ∘
          (fun b : Basis n => b.set tgt (b.get tgt != b.get c)) = fun b => b := by
        funext b
        simp only [Function.comp_apply]
        by_cases htn : tgt < n
        · rw [Basis.get_set_same _ _ _ htn, Basis.get_set_ne _ _ _ _ hct]
          have : ((b.get tgt != b.get c) != b.get c) = b.get tgt := by
            cases b.get tgt <;> cases b.get c <;> rfl
          rw [this, Basis.set_set, Basis.set_get_self]
        · simp [Basis.set_out_of_range _ _ _ htn]
      rw [this, permMatrix_id]
  | ccx c₁ c₂ tgt =>
      obtain ⟨h12, h1t, h2t⟩ := hwf
      rw [gateUnitary, permMatrix_mul]
      have : (fun b : Basis n => b.set tgt (b.get tgt != (b.get c₁ && b.get c₂))) ∘
          (fun b : Basis n => b.set tgt (b.get tgt != (b.get c₁ && b.get c₂))) = fun b => b := by
        funext b
        simp only [Function.comp_apply]
        by_cases htn : tgt < n
        · rw [Basis.get_set_same _ _ _ htn, Basis.get_set_ne _ _ _ _ h1t,
            Basis.get_set_ne _ _ _ _ h2t]
          have : ((b.get tgt != (b.get c₁ && b.get c₂)) != (b.get c₁ && b.get c₂)) =
              b.get tgt := by
            cases b.get tgt <;> cases b.get c₁ <;> cases b.get c₂ <;> rfl
          rw [this, Basis.set_set, Basis.set_get_self]
        · simp [Basis.set_out_of_range _ _ _ htn]
      rw [this, permMatrix_id]
  | cz c tgt =>
      rw [gateUnitary, phaseMatrix_mul]
      have : (fun b : Basis n => (if b.get c && b.get tgt then (-1 : ℂ) else 1) *
          (if b.get c && b.get tgt then (-1 : ℂ) else 1)) = fun _ => 1 := by
        funext b
        by_cases hb : b.get c && b.get tgt <;> simp [hb]
      rw [this, phaseMatrix_one]
  | ccz c₁ c₂ tgt =>
      rw [gateUnitary, phaseMatrix_mul]
      have : (fun b : Basis n =>
          (if b.get c₁ && b.get c₂ && b.get tgt then (-1 : ℂ) else 1) *
          (if b.get c₁ && b.get c₂ && b.get tgt then (-1 : ℂ) else 1)) = fun _ => 1 := by
        funext b
        by_cases hb : b.get c₁ && b.get c₂ && b.get tgt <;> simp [hb]
      rw [this, phaseMatrix_one]
  | s q => simp [Gate.isSelfInverse] at hs
  | sdg q => simp [Gate.isSelfInverse] at hs
  | t q => simp [Gate.isSelfInverse] at hs
  | tdg q => simp [Gate.isSelfInverse] at hs
  | rz θ q => simp [Gate.isSelfInverse] at hs
  | measure q c => simp [Gate.isSelfInverse] at hs
  | reset q => simp [Gate.isSelfInverse] at hs

/-! ## Operand symmetries -/

theorem gateUnitary_cz_comm (n : Nat) (c tgt : Qubit) :
    gateUnitary n (.cz c tgt) = gateUnitary n (.cz tgt c) := by
  simp only [gateUnitary]
  congr 1
  funext b
  rw [Bool.and_comm]

theorem gateUnitary_ccz_swap₁₂ (n : Nat) (a b c : Qubit) :
    gateUnitary n (.ccz a b c) = gateUnitary n (.ccz b a c) := by
  simp only [gateUnitary]
  congr 1
  funext s
  rw [Bool.and_comm (s.get a)]

theorem gateUnitary_ccz_swap₂₃ (n : Nat) (a b c : Qubit) :
    gateUnitary n (.ccz a b c) = gateUnitary n (.ccz a c b) := by
  simp only [gateUnitary]
  congr 1
  funext s
  rw [Bool.and_assoc, Bool.and_comm (s.get b), ← Bool.and_assoc]

/-- Writes to different wires commute. -/
theorem Basis.set_comm {n : Nat} (b : Basis n) {q r : Qubit} (h : q ≠ r) (v w : Bool) :
    (b.set q v).set r w = (b.set r w).set q v := by
  funext i
  by_cases hi : (i : Nat) = q <;> by_cases hi' : (i : Nat) = r <;>
    simp_all [Basis.set]

/-- A one-wire embedding at a wire that does not exist is the identity map on states. -/
theorem Basis.set_eq_self_of_ge {n : Nat} (b : Basis n) (q : Qubit) (v : Bool) (h : ¬ q < n) :
    b.set q v = b := Basis.set_out_of_range b q v h

/-! ## Commutation shapes

Three ways two gate matrices can commute without their wires being disjoint. Together with
`SupportedOn.mul_comm` these discharge every entry of the `CancelGates` commutation tables.
-/

/-- The matrix is diagonal in the computational basis. -/
def IsDiagonal {n : Nat} (U : Density n) : Prop := ∀ out inp, out ≠ inp → U out inp = 0

theorem isDiagonal_phaseMatrix {n : Nat} (f : Basis n → ℂ) : IsDiagonal (phaseMatrix f) := by
  intro out inp h; simp [phaseMatrix, h]

theorem isDiagonal_one {n : Nat} : IsDiagonal (1 : Density n) := by
  intro out inp h; simp [h]

theorem isDiagonal_embed1_diag2 {n : Nat} (a b : ℂ) (q : Qubit) :
    IsDiagonal (embed1 n (diag2 a b) q) := by
  by_cases hq : q < n
  · intro out inp hne
    rw [embed1_apply_of_lt _ hq]
    by_cases hall : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
    · rw [if_pos hall]
      have : out ⟨q, hq⟩ ≠ inp ⟨q, hq⟩ := by
        intro hc
        refine hne (funext fun r => ?_)
        by_cases hr : (r : Nat) = q
        · have : r = ⟨q, hq⟩ := Fin.ext hr
          rw [this]; exact hc
        · exact hall r hr
      simp [diag2, this]
    · simp [hall]
  · simpa [embed1_eq_one _ hq] using isDiagonal_one (n := n)

/-- **Diagonal matrices commute.** -/
theorem IsDiagonal.mul_comm {n : Nat} {U V : Density n} (hU : IsDiagonal U) (hV : IsDiagonal V) :
    U * V = V * U := by
  have key : ∀ A B : Density n, IsDiagonal A → IsDiagonal B → ∀ out inp : Basis n,
      (A * B) out inp = if out = inp then A out out * B out out else 0 := by
    intro A B hA hB out inp
    simp only [Matrix.mul_apply]
    by_cases h : out = inp
    · subst h
      rw [if_pos rfl, Finset.sum_eq_single out]
      · intro k _ hk
        rw [hA out k (fun hc => hk hc.symm), zero_mul]
      · simp
    · rw [if_neg h]
      refine Finset.sum_eq_zero fun k _ => ?_
      by_cases hk : out = k
      · subst hk
        rw [hB out inp h, mul_zero]
      · rw [hA out k hk, zero_mul]
  funext out inp
  rw [key U V hU hV, key V U hV hU]
  by_cases h : out = inp <;> simp [h, _root_.mul_comm]

/-- **A permutation commutes with a diagonal it does not disturb.** -/
theorem comm_perm_phase {n : Nat} (σ : Basis n → Basis n) (f : Basis n → ℂ)
    (h : ∀ b, f (σ b) = f b) :
    permMatrix σ * phaseMatrix f = phaseMatrix f * permMatrix σ := by
  funext out inp
  have hL : (permMatrix σ * phaseMatrix f) out inp = (if out = σ inp then 1 else 0) * f inp := by
    simp only [Matrix.mul_apply, permMatrix, phaseMatrix]
    rw [Finset.sum_eq_single inp]
    · simp
    · intro k _ hk; simp [hk]
    · simp
  have hR : (phaseMatrix f * permMatrix σ) out inp = f out * (if out = σ inp then 1 else 0) := by
    simp only [Matrix.mul_apply, permMatrix, phaseMatrix]
    rw [Finset.sum_eq_single out]
    · simp
    · intro k _ hk; simp [Ne.symm hk]
    · simp
  rw [hL, hR]
  by_cases hc : out = σ inp
  · subst hc; simp [h inp]
  · simp [hc]

/-- `X` as a basis permutation. -/
theorem embed1_x2_eq_perm (n : Nat) (q : Qubit) :
    embed1 n x2 q = permMatrix (fun b : Basis n => b.set q (!b.get q)) := by
  by_cases hq : q < n
  · funext out inp
    simp only [embed1_apply_of_lt _ hq, permMatrix, x2]
    by_cases hall : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
    · rw [if_pos hall]
      by_cases hb : out ⟨q, hq⟩ = !inp ⟨q, hq⟩
      · have : out = inp.set q (!inp.get q) := by
          funext r
          by_cases hr : (r : Nat) = q
          · have hrq : r = ⟨q, hq⟩ := Fin.ext hr
            subst hrq
            simp [Basis.set, Basis.get, hq, hb]
          · simp [Basis.set, hr, hall r hr]
        rw [if_pos hb, if_pos this]
      · have : out ≠ inp.set q (!inp.get q) := by
          intro hc
          refine hb ?_
          have := congrFun hc ⟨q, hq⟩
          simpa [Basis.set, Basis.get, hq] using this
        rw [if_neg hb, if_neg this]
    · rw [if_neg hall]
      have : out ≠ inp.set q (!inp.get q) := by
        intro hc
        refine hall fun r hr => ?_
        rw [hc]
        simp [Basis.set, hr]
      simp [this]
  · rw [embed1_eq_one _ hq]
    funext out inp
    simp [permMatrix, Basis.set_out_of_range _ _ _ hq, Matrix.one_apply]

/-! ## Concrete phase values -/

theorem ep_eq_cos_sin (θ : ℚ) :
    ep θ = (Real.cos (Real.pi * (θ : ℝ)) : ℂ) + (Real.sin (Real.pi * (θ : ℝ)) : ℂ) * Complex.I := by
  simp only [ep]
  rw [Complex.exp_mul_I]
  simp

@[simp] theorem ep_one : ep 1 = -1 := by
  rw [ep_eq_cos_sin]
  norm_num

@[simp] theorem ep_half : ep (1/2) = Complex.I := by
  rw [ep_eq_cos_sin, show Real.pi * ((1/2 : ℚ) : ℝ) = Real.pi / 2 by push_cast; ring]
  simp

@[simp] theorem ep_neg_half : ep (-1/2) = -Complex.I := by
  rw [ep_eq_cos_sin]
  rw [show Real.pi * ((-1/2 : ℚ) : ℝ) = -(Real.pi / 2) by push_cast; ring]
  simp

theorem ep_quarter : ep (1/4) = ((Real.sqrt 2 / 2 : ℝ) : ℂ) * (1 + Complex.I) := by
  rw [ep_eq_cos_sin]
  rw [show Real.pi * ((1/4 : ℚ) : ℝ) = Real.pi / 4 by push_cast; ring]
  rw [Real.cos_pi_div_four, Real.sin_pi_div_four]
  push_cast
  ring

theorem ep_neg_quarter : ep (-1/4) = ((Real.sqrt 2 / 2 : ℝ) : ℂ) * (1 - Complex.I) := by
  rw [ep_eq_cos_sin]
  rw [show Real.pi * ((-1/4 : ℚ) : ℝ) = -(Real.pi / 4) by push_cast; ring]
  rw [Real.cos_neg, Real.sin_neg, Real.cos_pi_div_four, Real.sin_pi_div_four]
  push_cast
  ring

/-! ## The Hadamard-reduction identities

The local Clifford identities the `reduce_hadamards` sweep of `src/cancel.rs` applies. Each
is a `2 × 2` fact; `embed1_mul` lifts it to the full register. The `S` rules hold only up to
the global phase `e^{±iπ/4}`, which is invisible to the density-matrix semantics — see
`equivalent_of_unitary_smul`. -/

/-- `H·X·H = Z`. -/
theorem mul2_hxh : mul2 (mul2 h2 x2) h2 = diag2 1 (ep 1) := by
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) ^ 2 = 2 := by
    rw [pow_two]; exact ofReal_sqrt_two_mul_self
  have hne := ofReal_sqrt_two_ne_zero
  funext out inp
  cases out <;> cases inp <;>
    · simp only [mul2, h2, x2, diag2, Fintype.sum_bool, ep_one]
      norm_num
      try field_simp
      try rw [hsq]
      try ring

/-- `H·Z·H = X`. -/
theorem mul2_hzh : mul2 (mul2 h2 (diag2 1 (ep 1))) h2 = x2 := by
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) ^ 2 = 2 := by
    rw [pow_two]; exact ofReal_sqrt_two_mul_self
  have hne := ofReal_sqrt_two_ne_zero
  funext out inp
  cases out <;> cases inp <;>
    · simp only [mul2, h2, x2, diag2, Fintype.sum_bool, ep_one]
      norm_num
      try field_simp
      try rw [hsq]
      try ring

/-- `H·S·H = e^{iπ/4} · Sdg·H·Sdg`. -/
theorem mul2_hsh : mul2 (mul2 h2 (diag2 1 (ep (1/2)))) h2 =
    fun o i => ep (1/4) * mul2 (mul2 (diag2 1 (ep (-1/2))) h2) (diag2 1 (ep (-1/2))) o i := by
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) ^ 2 = 2 := by
    rw [pow_two]; exact ofReal_sqrt_two_mul_self
  have hne := ofReal_sqrt_two_ne_zero
  have hI2 : Complex.I ^ 2 = -1 := Complex.I_sq
  have hI3 : Complex.I ^ 3 = -Complex.I := by rw [pow_succ, Complex.I_sq]; ring
  funext out inp
  cases out <;> cases inp <;>
    · simp only [mul2, h2, diag2, Fintype.sum_bool, ep_half, ep_neg_half, ep_quarter]
      norm_num
      try field_simp
      try ring_nf
      try rw [hsq]
      try simp only [hI3, hI2]
      try ring

/-- `H·Sdg·H = e^{-iπ/4} · S·H·S`. -/
theorem mul2_hsdgh : mul2 (mul2 h2 (diag2 1 (ep (-1/2)))) h2 =
    fun o i => ep (-1/4) * mul2 (mul2 (diag2 1 (ep (1/2))) h2) (diag2 1 (ep (1/2))) o i := by
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) ^ 2 = 2 := by
    rw [pow_two]; exact ofReal_sqrt_two_mul_self
  have hne := ofReal_sqrt_two_ne_zero
  have hI2 : Complex.I ^ 2 = -1 := Complex.I_sq
  have hI3 : Complex.I ^ 3 = -Complex.I := by rw [pow_succ, Complex.I_sq]; ring
  funext out inp
  cases out <;> cases inp <;>
    · simp only [mul2, h2, diag2, Fintype.sum_bool, ep_half, ep_neg_half, ep_neg_quarter]
      norm_num
      try field_simp
      try ring_nf
      try rw [hsq]
      try simp only [hI3, hI2]
      try ring

end
end TzapLean
