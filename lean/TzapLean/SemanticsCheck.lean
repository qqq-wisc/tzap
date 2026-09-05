import TzapLean.Cancel

/-!
# Independent Validation of the Semantics

The semantics of `TzapLean/Semantics.lean` is the foundation everything else stands on, so
this file checks it from angles the pass proofs never exercise:

1. **Unitarity.** Every gate matrix satisfies `UᴴU = 1`. A single mistyped entry — a wrong
   sign, a missing `1/√2`, a swapped control and target — breaks this.
2. **Trace preservation, unconditionally.** With unitarity in hand, the `totalTrace` results
   of `Semantics.lean` lose their hypothesis: every gate, unitary or not, is a channel.
3. **Concrete states.** Explicit amplitudes for explicit circuits: `X|0⟩ = |1⟩`, the Bell
   state, `CNOT` flipping the target and not the control, measurement statistics of `|+⟩`,
   `reset` returning `|1⟩` to `|0⟩`.
4. **The `src/unitary.rs` suite.** Rust's own semantics tests, ported — including the
   *negative* ones, which are what rule out a degenerate semantics that equates everything.
-/

namespace TzapLean

open Matrix

noncomputable section

/-! ## Every gate matrix is unitary -/

theorem embed1_conjTranspose (n : Nat) (M : Bool → Bool → ℂ) (q : Qubit) :
    (embed1 n M q)ᴴ = embed1 n (fun o i => star (M i o)) q := by
  by_cases hq : q < n
  · funext out inp
    simp only [Matrix.conjTranspose_apply, embed1_apply_of_lt _ hq]
    by_cases h : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
    · rw [if_pos h, if_pos fun r hr => (h r hr).symm]
    · rw [if_neg h, if_neg fun hc => h fun r hr => (hc r hr).symm]
      simp
  · simp [embed1_eq_one _ hq]

/-- A permutation matrix of an involution is self-adjoint. -/
theorem permMatrix_conjTranspose_of_involutive {n : Nat} {σ : Basis n → Basis n}
    (h : ∀ b, σ (σ b) = b) : (permMatrix σ)ᴴ = permMatrix σ := by
  funext out inp
  simp only [Matrix.conjTranspose_apply, permMatrix]
  by_cases hc : out = σ inp
  · rw [if_pos hc, if_pos (by rw [hc, h inp])]
    simp
  · rw [if_neg hc, if_neg ?_]
    · simp
    · intro hcon
      exact hc (by rw [hcon, h out])

/-- A real diagonal is self-adjoint. -/
theorem phaseMatrix_conjTranspose_of_real {n : Nat} {f : Basis n → ℂ}
    (h : ∀ b, star (f b) = f b) : (phaseMatrix f)ᴴ = phaseMatrix f := by
  funext out inp
  simp only [Matrix.conjTranspose_apply, phaseMatrix]
  by_cases hc : out = inp
  · subst hc; simp [h]
  · rw [if_neg hc, if_neg (Ne.symm hc)]
    simp

/-- **Every gate's matrix is unitary.** For the two- and three-qubit permutation gates this
needs `Gate.Wf`: `cnot q q` really is not unitary. -/
theorem gateUnitary_unitary (n : Nat) (g : Gate) (hwf : g.Wf) :
    (gateUnitary n g)ᴴ * gateUnitary n g = 1 := by
  have embed : ∀ (M : Bool → Bool → ℂ) (q : Qubit),
      mul2 (fun o i => star (M i o)) M = id2 →
      (embed1 n M q)ᴴ * embed1 n M q = 1 := by
    intro M q hM
    rw [embed1_conjTranspose, embed1_mul, hM, embed1_id2]
  have starmul : ∀ a : ℚ, star (ep a) * ep a = 1 := by
    intro a; rw [mul_comm]; exact ep_mul_star a
  have diagstar2 : ∀ x y : ℚ,
      mul2 (fun o i => star (diag2 (ep x) (ep y) i o)) (diag2 (ep x) (ep y)) = id2 := by
    intro x y
    funext o i
    cases o <;> cases i <;>
      · simp only [mul2, diag2, id2, Fintype.sum_bool]
        norm_num
        try exact starmul x
        try exact starmul y
  have diagstar : ∀ a : ℚ, mul2 (fun o i => star (diag2 1 (ep a) i o)) (diag2 1 (ep a)) = id2 := by
    intro a
    have := diagstar2 0 a
    simpa using this
  cases g with
  | x q =>
      refine embed x2 q ?_
      funext o i
      cases o <;> cases i <;> simp [mul2, x2, id2]
  | h q =>
      refine embed h2 q ?_
      have hsq : ((Real.sqrt 2 : ℝ) : ℂ) ^ 2 = 2 := by
        rw [pow_two]; exact ofReal_sqrt_two_mul_self
      have hne := ofReal_sqrt_two_ne_zero
      funext o i
      cases o <;> cases i <;>
        · simp only [mul2, h2, id2, Fintype.sum_bool, RCLike.star_def, map_div₀,
            Complex.conj_ofReal]
          norm_num
          try field_simp
          try exact hsq.symm
          try (rw [pow_two, hsq])
          try rw [hsq]
          try ring
  | s q => exact embed _ q (diagstar (1/2))
  | sdg q => exact embed _ q (diagstar (-1/2))
  | z q => exact embed _ q (diagstar 1)
  | t q => exact embed _ q (diagstar (1/4))
  | tdg q => exact embed _ q (diagstar (-1/4))
  | rz a q => exact embed _ q (diagstar2 (-a/2) (a/2))
  | cnot c t =>
      have hinv : ∀ b : Basis n,
          (fun b : Basis n => b.set t (b.get t != b.get c))
            ((fun b : Basis n => b.set t (b.get t != b.get c)) b) = b := by
        have := gateUnitary_mul_self n (Gate.cnot c t) hwf rfl
        intro b
        have hcomp : (fun b : Basis n => b.set t (b.get t != b.get c)) ∘
            (fun b : Basis n => b.set t (b.get t != b.get c)) = fun b => b := by
          have hct : c ≠ t := hwf
          funext b
          simp only [Function.comp_apply]
          by_cases htn : t < n
          · rw [Basis.get_set_same _ _ _ htn, Basis.get_set_ne _ _ _ _ hct]
            have : ((b.get t != b.get c) != b.get c) = b.get t := by
              cases b.get t <;> cases b.get c <;> rfl
            rw [this, Basis.set_set, Basis.set_get_self]
          · simp [Basis.set_out_of_range _ _ _ htn]
        exact congrFun hcomp b
      rw [gateUnitary, permMatrix_conjTranspose_of_involutive hinv]
      exact gateUnitary_mul_self n (Gate.cnot c t) hwf rfl
  | ccx c₁ c₂ t =>
      obtain ⟨h12, h1t, h2t⟩ := hwf
      have hinv : ∀ b : Basis n,
          (fun b : Basis n => b.set t (b.get t != (b.get c₁ && b.get c₂)))
            ((fun b : Basis n => b.set t (b.get t != (b.get c₁ && b.get c₂))) b) = b := by
        intro b
        simp only
        by_cases htn : t < n
        · rw [Basis.get_set_same _ _ _ htn, Basis.get_set_ne _ _ _ _ h1t,
            Basis.get_set_ne _ _ _ _ h2t]
          have : ((b.get t != (b.get c₁ && b.get c₂)) != (b.get c₁ && b.get c₂)) = b.get t := by
            cases b.get t <;> cases b.get c₁ <;> cases b.get c₂ <;> rfl
          rw [this, Basis.set_set, Basis.set_get_self]
        · simp [Basis.set_out_of_range _ _ _ htn]
      rw [gateUnitary, permMatrix_conjTranspose_of_involutive hinv]
      exact gateUnitary_mul_self n (Gate.ccx c₁ c₂ t) ⟨h12, h1t, h2t⟩ rfl
  | cz c t =>
      have hreal : ∀ b : Basis n,
          star (if b.get c && b.get t then (-1 : ℂ) else 1) =
            if b.get c && b.get t then (-1 : ℂ) else 1 := by
        intro b; by_cases h : b.get c && b.get t <;> simp [h]
      rw [gateUnitary, phaseMatrix_conjTranspose_of_real hreal]
      exact gateUnitary_mul_self n (Gate.cz c t) hwf rfl
  | ccz c₁ c₂ t =>
      have hreal : ∀ b : Basis n,
          star (if b.get c₁ && b.get c₂ && b.get t then (-1 : ℂ) else 1) =
            if b.get c₁ && b.get c₂ && b.get t then (-1 : ℂ) else 1 := by
        intro b; by_cases h : b.get c₁ && b.get c₂ && b.get t <;> simp [h]
      rw [gateUnitary, phaseMatrix_conjTranspose_of_real hreal]
      exact gateUnitary_mul_self n (Gate.ccz c₁ c₂ t) hwf rfl
  | measure q c => simp [gateUnitary]
  | reset q => simp [gateUnitary]

/-! ## Trace preservation, with no hypothesis left -/

/-- Every gate — unitary, `measure` or `reset` — preserves total probability. -/
theorem totalTrace_step {n m : Nat} (g : Gate) (hwf : g.Wf) (ρ : CQState n m) :
    (step g ρ).totalTrace = ρ.totalTrace := by
  cases g with
  | measure q c => exact totalTrace_step_measure q c ρ
  | reset q => exact totalTrace_step_reset q ρ
  | _ => exact totalTrace_step_of_isUnitary rfl ρ (gateUnitary_unitary n _ hwf)

/-- And therefore so does every circuit. -/
theorem totalTrace_denote {n m : Nat} : ∀ (gs : List Gate), (∀ g ∈ gs, g.Wf) →
    ∀ ρ : CQState n m, (denote gs ρ).totalTrace = ρ.totalTrace := by
  intro gs
  induction gs with
  | nil => intro _ ρ; rfl
  | cons g gs ih =>
      intro hwf ρ
      rw [denote_cons, ih (fun g' hg' => hwf g' (by simp [hg'])),
        totalTrace_step g (hwf g (by simp))]

/-! ## Evaluating concrete circuits

Three lemmas turn a matrix product into arithmetic: a one-wire gate on the left sums over
that wire's two values (this is exactly the Rust `Mat::apply_single`), a permutation gate
reads one entry, and a diagonal scales one. -/

theorem embed1_mul_apply {n : Nat} (M : Bool → Bool → ℂ) (q : Qubit) (hq : q < n)
    (A : Density n) (out inp : Basis n) :
    (embed1 n M q * A) out inp = ∑ v : Bool, M (out ⟨q, hq⟩) v * A (out.set q v) inp := by
  simp only [Matrix.mul_apply, embed1_apply_of_lt _ hq]
  set kf : Basis n := out.set q false with hkf
  set kt : Basis n := out.set q true with hkt
  have hkf_q : kf ⟨q, hq⟩ = false := by simp [hkf, Basis.set]
  have hkt_q : kt ⟨q, hq⟩ = true := by simp [hkt, Basis.set]
  have hne : kf ≠ kt := by
    intro hc
    have hb : kf ⟨q, hq⟩ = kt ⟨q, hq⟩ := by rw [hc]
    rw [hkf_q, hkt_q] at hb
    exact Bool.false_ne_true hb
  have key : ∀ k : Basis n,
      (if ∀ r : Fin n, (r : Nat) ≠ q → out r = k r then M (out ⟨q, hq⟩) (k ⟨q, hq⟩) else 0) *
        A k inp
      = (if k = kf then M (out ⟨q, hq⟩) false * A kf inp else 0) +
        (if k = kt then M (out ⟨q, hq⟩) true * A kt inp else 0) := by
    intro k
    by_cases hk : ∀ r : Fin n, (r : Nat) ≠ q → out r = k r
    · have hsplit : k = kf ∨ k = kt := by
        cases hkq : k ⟨q, hq⟩
        · left
          funext r
          by_cases hr : (r : Nat) = q
          · have : r = ⟨q, hq⟩ := Fin.ext hr
            rw [this, hkq, hkf_q]
          · rw [hkf, Basis.set, if_neg hr, ← hk r hr]
        · right
          funext r
          by_cases hr : (r : Nat) = q
          · have : r = ⟨q, hq⟩ := Fin.ext hr
            rw [this, hkq, hkt_q]
          · rw [hkt, Basis.set, if_neg hr, ← hk r hr]
      rcases hsplit with rfl | rfl
      · rw [if_pos hk, if_pos rfl, if_neg hne, hkf_q]
        simp
      · rw [if_pos hk, if_neg (Ne.symm hne), if_pos rfl, hkt_q]
        simp
    · have hkf' : k ≠ kf := by
        intro hc
        exact hk fun r hr => by rw [hc, hkf, Basis.set, if_neg hr]
      have hkt' : k ≠ kt := by
        intro hc
        exact hk fun r hr => by rw [hc, hkt, Basis.set, if_neg hr]
      rw [if_neg hk, if_neg hkf', if_neg hkt']
      simp
  rw [Finset.sum_congr rfl fun k _ => key k, Finset.sum_add_distrib,
    Finset.sum_ite_eq' Finset.univ, Finset.sum_ite_eq' Finset.univ, Fintype.sum_bool]
  simp only [Finset.mem_univ, if_true, hkf, hkt]
  ring

theorem permMatrix_mul_apply {n : Nat} {σ : Basis n → Basis n} (hinv : ∀ b, σ (σ b) = b)
    (A : Density n) (out inp : Basis n) : (permMatrix σ * A) out inp = A (σ out) inp := by
  simp only [Matrix.mul_apply, permMatrix]
  rw [Finset.sum_eq_single (σ out)]
  · rw [hinv out, if_pos rfl, one_mul]
  · intro k _ hk
    rw [if_neg, zero_mul]
    intro hc
    exact hk (by rw [hc, hinv])
  · simp

theorem phaseMatrix_mul_apply {n : Nat} (f : Basis n → ℂ) (A : Density n) (out inp : Basis n) :
    (phaseMatrix f * A) out inp = f out * A out inp := by
  simp only [Matrix.mul_apply, phaseMatrix]
  rw [Finset.sum_eq_single out]
  · simp
  · intro k _ hk; simp [Ne.symm hk]
  · simp

/-- `|b⟩⟨b'|`. -/
def ketBra {n : Nat} (b b' : Basis n) : Density n :=
  fun out inp => if out = b ∧ inp = b' then 1 else 0

/-- Conjugating `|b⟩⟨b'|` reads off one column of `U` and one of its conjugate. -/
theorem conj_ketBra {n : Nat} (U : Density n) (b b' out inp : Basis n) :
    conj U (ketBra b b') out inp = U out b * star (U inp b') := by
  simp only [conj, Matrix.mul_apply, Matrix.conjTranspose_apply, ketBra]
  rw [Finset.sum_eq_single b']
  · rw [Finset.sum_eq_single b]
    · simp
    · intro k _ hk; simp [hk]
    · simp
  · intro l _ hl
    rw [Finset.sum_eq_zero, zero_mul]
    intro k _
    simp [hl]
  · simp

/-! ## Concrete circuits on one and two wires -/

theorem basis_one_ext {b b' : Basis 1} (h : b 0 = b' 0) : b = b' := by
  funext r
  have : r = 0 := Fin.ext (by omega)
  rw [this]; exact h

theorem basis_two_ext {b b' : Basis 2} (h0 : b 0 = b' 0) (h1 : b 1 = b' 1) : b = b' := by
  funext r
  match r with
  | 0 => exact h0
  | 1 => exact h1

theorem basis_one_eq_true_iff (b : Basis 1) : b = ![true] ↔ b 0 = true := by
  constructor
  · intro h; rw [h]; rfl
  · intro h; exact basis_one_ext (by rw [h]; rfl)

theorem basis_one_eq_false_iff (b : Basis 1) : b = ![false] ↔ b 0 = false := by
  constructor
  · intro h; rw [h]; rfl
  · intro h; exact basis_one_ext (by rw [h]; rfl)

theorem embed1_apply_one (M : Bool → Bool → ℂ) (out inp : Basis 1) :
    embed1 1 M 0 out inp = M (out 0) (inp 0) := by
  simp only [embed1_apply_of_lt _ (by norm_num : (0:Nat) < 1)]
  rw [if_pos]
  · rfl
  · intro r hr
    exact absurd (Fin.val_eq_zero r) hr

theorem embed1_apply_two_zero (M : Bool → Bool → ℂ) (out inp : Basis 2) :
    embed1 2 M 0 out inp = if out 1 = inp 1 then M (out 0) (inp 0) else 0 := by
  simp only [embed1_apply_of_lt _ (by norm_num : (0:Nat) < 2)]
  by_cases h : out 1 = inp 1
  · rw [if_pos, if_pos h]
    · rfl
    · intro r hr
      match r with
      | 0 => exact absurd rfl hr
      | 1 => exact h
  · rw [if_neg h, if_neg]
    intro hall
    exact h (hall 1 (by norm_num))

theorem embed1_apply_two_one (M : Bool → Bool → ℂ) (out inp : Basis 2) :
    embed1 2 M 1 out inp = if out 0 = inp 0 then M (out 1) (inp 1) else 0 := by
  simp only [embed1_apply_of_lt _ (by norm_num : (1:Nat) < 2)]
  by_cases h : out 0 = inp 0
  · rw [if_pos, if_pos h]
    · rfl
    · intro r hr
      match r with
      | 0 => exact h
      | 1 => exact absurd rfl hr
  · rw [if_neg h, if_neg]
    intro hall
    exact h (hall 0 (by norm_num))

/-! ### `X` flips, `T` phases, and the list order is execution order -/

/-- `X|0⟩ = |1⟩`. -/
theorem x_maps_zero_to_one :
    conj (gateUnitary 1 (.x 0)) (ketBra ![false] ![false]) = ketBra ![true] ![true] := by
  funext out inp
  rw [conj_ketBra]
  have hcol : ∀ b : Basis 1, gateUnitary 1 (.x 0) b ![false] = if b 0 = true then 1 else 0 := by
    intro b
    rw [gateUnitary, embed1_apply_one]
    simp [x2]
  rw [hcol, hcol]
  simp only [ketBra, basis_one_eq_true_iff]
  by_cases h1 : out 0 = true <;> by_cases h2 : inp 0 = true <;> simp [h1, h2]

/-- `T` multiplies the `|1⟩` branch by `e^{iπ/4}` — and not by its conjugate. -/
theorem t_phases_one_branch :
    conj (gateUnitary 1 (.t 0)) (ketBra ![true] ![false]) ![true] ![false] = ep (1/4) := by
  rw [conj_ketBra, gateUnitary, embed1_apply_one, embed1_apply_one]
  simp [diag2]

/-- The list is executed left to right: `[x, h]` differs from `[h, x]`. -/
theorem list_order_is_execution_order :
    unitary 1 [Gate.x 0, Gate.h 0] ![true] ![false] ≠
      unitary 1 [Gate.h 0, Gate.x 0] ![true] ![false] := by
  have hsqrt : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := ofReal_sqrt_two_ne_zero
  have hxh : unitary 1 [Gate.x 0, Gate.h 0] ![true] ![false] = -1 / (Real.sqrt 2 : ℂ) := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary, embed1_mul]
    rw [embed1_apply_one]
    simp [mul2, h2, x2]
  have hhx : unitary 1 [Gate.h 0, Gate.x 0] ![true] ![false] = 1 / (Real.sqrt 2 : ℂ) := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary, embed1_mul]
    rw [embed1_apply_one]
    simp [mul2, h2, x2]
  rw [hxh, hhx]
  intro hc
  field_simp at hc
  exact absurd hc (by norm_num)

/-! ### Two wires: `CNOT`'s direction, `CZ`'s phase, and the Bell state -/

/-- The two components of the `CNOT` permutation. -/
theorem cnot_perm_zero (b : Basis 2) : (b.set 1 (b.get 1 != b.get 0)) 0 = b 0 := by
  simp [Basis.set]

theorem cnot_perm_one (b : Basis 2) : (b.set 1 (b.get 1 != b.get 0)) 1 = (b 1 != b 0) := by
  simp [Basis.set, Basis.get]

theorem cnot_perm_apply (b : Basis 2) :
    b.set 1 (b.get 1 != b.get 0) = ![b 0, (b 1 != b 0)] :=
  basis_two_ext (by rw [cnot_perm_zero]; rfl) (by rw [cnot_perm_one]; rfl)

/-- `CNOT 0 1` flips wire 1 when wire 0 is set — and not the other way round. -/
theorem cnot_flips_target :
    conj (gateUnitary 2 (.cnot 0 1)) (ketBra ![true, false] ![true, false]) =
      ketBra ![true, true] ![true, true] := by
  have hσ : Basis.set (![true, false] : Basis 2) 1
      (Basis.get (![true, false] : Basis 2) 1 != Basis.get (![true, false] : Basis 2) 0) =
      ![true, true] := by
    rw [cnot_perm_apply]; rfl
  funext out inp
  rw [conj_ketBra]
  have hcol : ∀ b : Basis 2,
      gateUnitary 2 (.cnot 0 1) b ![true, false] = if b = ![true, true] then 1 else 0 := by
    intro b
    simp only [gateUnitary, permMatrix, hσ]
  rw [hcol, hcol]
  simp only [ketBra]
  by_cases h1 : out = ![true, true] <;> by_cases h2 : inp = ![true, true] <;> simp [h1, h2]

/-- With the control clear, `CNOT` does nothing. -/
theorem cnot_leaves_clear_control :
    conj (gateUnitary 2 (.cnot 0 1)) (ketBra ![false, true] ![false, true]) =
      ketBra ![false, true] ![false, true] := by
  have hσ : Basis.set (![false, true] : Basis 2) 1
      (Basis.get (![false, true] : Basis 2) 1 != Basis.get (![false, true] : Basis 2) 0) =
      ![false, true] := by
    rw [cnot_perm_apply]; rfl
  funext out inp
  rw [conj_ketBra]
  have hcol : ∀ b : Basis 2,
      gateUnitary 2 (.cnot 0 1) b ![false, true] = if b = ![false, true] then 1 else 0 := by
    intro b
    simp only [gateUnitary, permMatrix, hσ]
  rw [hcol, hcol]
  simp only [ketBra]
  by_cases h1 : out = ![false, true] <;> by_cases h2 : inp = ![false, true] <;> simp [h1, h2]

/-- `CZ` is `-1` exactly on `|11⟩`. -/
theorem cz_phase (out : Basis 2) :
    gateUnitary 2 (.cz 0 1) out out = if out 0 && out 1 then -1 else 1 := by
  simp [gateUnitary, phaseMatrix, Basis.get]

/-- **The Bell state.** `H` on wire 0 then `CNOT 0 1`, applied to `|00⟩`, gives the density
matrix with `1/2` on the four entries linking `|00⟩` and `|11⟩` and `0` everywhere else.
This one check pins `H`'s normalization, `CNOT`'s direction, and the fact that a gate list
is executed left to right. -/
theorem bell_state (out inp : Basis 2) :
    conj (unitary 2 [Gate.h 0, Gate.cnot 0 1]) (ketBra ![false, false] ![false, false]) out inp =
      if out 0 = out 1 ∧ inp 0 = inp 1 then 1/2 else 0 := by
  have hsqrt : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := ofReal_sqrt_two_ne_zero
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := ofReal_sqrt_two_mul_self
  have hinv : ∀ b : Basis 2,
      (fun s : Basis 2 => s.set 1 (s.get 1 != s.get 0))
        ((fun s : Basis 2 => s.set 1 (s.get 1 != s.get 0)) b) = b := by
    intro b
    simp only
    refine basis_two_ext ?_ ?_
    · rw [cnot_perm_zero, cnot_perm_zero]
    · rw [cnot_perm_one, cnot_perm_one, cnot_perm_zero]
      cases b 0 <;> cases b 1 <;> rfl
  have hcol : ∀ b : Basis 2, unitary 2 [Gate.h 0, Gate.cnot 0 1] b ![false, false] =
      if b 0 = b 1 then 1 / (Real.sqrt 2 : ℂ) else 0 := by
    intro b
    simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary]
    rw [permMatrix_mul_apply hinv, embed1_apply_two_zero]
    simp only [cnot_perm_zero, cnot_perm_one, Matrix.cons_val_zero, Matrix.cons_val_one,
      ]
    by_cases h : b 0 = b 1
    · rw [if_pos (by rw [← h]; simp), if_pos h]
      cases hb : b 0 <;> simp [h2]
    · rw [if_neg (by simp; intro hc; exact h hc.symm), if_neg h]
  rw [conj_ketBra, hcol, hcol]
  have hstar : star ((1 : ℂ) / ((Real.sqrt 2 : ℝ) : ℂ)) = 1 / ((Real.sqrt 2 : ℝ) : ℂ) := by
    simp [Complex.conj_ofReal]
  by_cases h1 : out 0 = out 1 <;> by_cases h2 : inp 0 = inp 1
  · rw [if_pos h1, if_pos h2, if_pos (And.intro h1 h2), hstar, div_mul_div_comm, one_mul, hsq]
  · rw [if_pos h1, if_neg h2,
      if_neg (show ¬(out 0 = out 1 ∧ inp 0 = inp 1) from fun hc => h2 hc.2)]
    simp
  · rw [if_neg h1, if_pos h2,
      if_neg (show ¬(out 0 = out 1 ∧ inp 0 = inp 1) from fun hc => h1 hc.1)]
    simp
  · rw [if_neg h1, if_neg h2,
      if_neg (show ¬(out 0 = out 1 ∧ inp 0 = inp 1) from fun hc => h1 hc.1)]
    simp

/-! ### Measurement and reset

The clauses of `step` that are not conjugation by a unitary: the statistics of measuring
`|+⟩`, the coherence it destroys, and `reset` returning `|1⟩` to `|0⟩`. -/

theorem sum_basis_one (f : Basis 1 → ℂ) : ∑ b, f b = f ![true] + f ![false] := by
  rw [← Equiv.sum_comp
    ({ toFun := fun v => ![v], invFun := fun b => b 0, left_inv := fun _ => rfl,
       right_inv := fun b => basis_one_ext rfl } : Bool ≃ Basis 1) f, Fintype.sum_bool]
  rfl

theorem conj_phaseMatrix {n : Nat} (f : Basis n → ℂ) (ρ : Density n) (out inp : Basis n) :
    conj (phaseMatrix f) ρ out inp = f out * ρ out inp * star (f inp) := by
  simp only [conj]
  rw [Matrix.mul_apply, Finset.sum_eq_single inp]
  · rw [phaseMatrix_mul_apply]
    simp [Matrix.conjTranspose_apply, phaseMatrix]
  · intro l _ hl
    rw [phaseMatrix_mul_apply]
    simp [Matrix.conjTranspose_apply, phaseMatrix, Ne.symm hl]
  · simp

/-- `|+⟩⟨+|`: every entry is `1/2`. -/
theorem plus_entries (out inp : Basis 1) :
    conj (gateUnitary 1 (.h 0)) (ketBra ![false] ![false]) out inp = 1/2 := by
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := ofReal_sqrt_two_mul_self
  have hcol : ∀ b : Basis 1,
      gateUnitary 1 (.h 0) b ![false] = 1 / ((Real.sqrt 2 : ℝ) : ℂ) := by
    intro b
    rw [gateUnitary, embed1_apply_one]
    simp [h2]
  rw [conj_ketBra, hcol, hcol]
  have hstar : star ((1 : ℂ) / ((Real.sqrt 2 : ℝ) : ℂ)) = 1 / ((Real.sqrt 2 : ℝ) : ℂ) := by
    simp [Complex.conj_ofReal]
  rw [hstar, div_mul_div_comm, one_mul, hsq]

/-- The measured state, with its classical bit cleared. -/
private def plusState : CQState 1 1 :=
  CQState.ofDensity (conj (gateUnitary 1 (.h 0)) (ketBra ![false] ![false]))

private theorem plus_block (w : Memory 1) :
    plusState (w.write 0 false) + plusState (w.write 0 true) =
      conj (gateUnitary 1 (.h 0)) (ketBra ![false] ![false]) := by
  have hmem : w.write 0 false = (fun _ => false) := by
    funext i
    have : i = 0 := Fin.ext (by omega)
    simp [Memory.write, this]
  have hmem' : w.write 0 true ≠ (fun _ => false) := by
    intro hc
    have := congrFun hc 0
    simp [Memory.write] at this
  simp only [plusState, CQState.ofDensity]
  rw [if_pos hmem, if_neg hmem', add_zero]

/-- **Measuring `|+⟩` gives each outcome probability `1/2`.** -/
theorem measure_plus_outcome (v : Bool) :
    (step (n := 1) (m := 1) (.measure 0 0) plusState).outcome ![v] = 1/2 := by
  have hlt : (0 : Nat) < 1 := by norm_num
  have hread : Memory.read (![v] : Memory 1) 0 = v := by simp [Memory.read]
  rw [step_measure_of_lt 0 0 hlt]
  simp only [CQState.outcome]
  rw [plus_block, hread]
  simp only [Matrix.trace, Matrix.diag, proj]
  rw [sum_basis_one, conj_phaseMatrix, conj_phaseMatrix, plus_entries, plus_entries]
  have h1 : Basis.get (![true] : Basis 1) 0 = true := by simp [Basis.get]
  have h0 : Basis.get (![false] : Basis 1) 0 = false := by simp [Basis.get]
  rw [h1, h0]
  cases v <;> norm_num

/-- **Measurement destroys the coherence of `|+⟩`**: the off-diagonal entry is gone. -/
theorem measure_kills_coherence (v : Bool) :
    (step (n := 1) (m := 1) (.measure 0 0) plusState) ![v] ![true] ![false] = 0 := by
  have hlt : (0 : Nat) < 1 := by norm_num
  have hread : Memory.read (![v] : Memory 1) 0 = v := by simp [Memory.read]
  rw [step_measure_of_lt 0 0 hlt]
  simp only
  rw [plus_block, hread, proj, conj_phaseMatrix]
  have h1 : Basis.get (![true] : Basis 1) 0 = true := by simp [Basis.get]
  have h0 : Basis.get (![false] : Basis 1) 0 = false := by simp [Basis.get]
  rw [h1, h0]
  cases v <;> norm_num

/-- **`reset` returns `|1⟩` to `|0⟩`.** -/
theorem reset_one_to_zero :
    step (n := 1) (m := 0) (.reset 0) (fun _ => ketBra ![true] ![true])
      = fun _ => ketBra ![false] ![false] := by
  funext w out inp
  have hset : Basis.set (![true] : Basis 1) 0 false = ![false] :=
    basis_one_ext (by simp [Basis.set])
  have hK : ∀ (b : Bool) (x : Basis 1),
      resetKraus 1 0 b x ![true] = if b = true ∧ x = ![false] then 1 else 0 := by
    intro b x
    simp only [resetKraus, Basis.get, dif_pos (by norm_num : (0:Nat) < 1), hset]
    cases b <;> simp
  have e : ∀ b : Bool,
      (resetKraus 1 0 b * (ketBra ![true] ![true] : Density 1) * (resetKraus 1 0 b)ᴴ) out inp
        = resetKraus 1 0 b out ![true] * star (resetKraus 1 0 b inp ![true]) :=
    fun b => conj_ketBra (resetKraus 1 0 b) ![true] ![true] out inp
  rw [step_reset]
  simp only
  rw [Matrix.sum_apply, Fintype.sum_bool, e true, e false, hK, hK, hK, hK]
  simp only [ketBra]
  by_cases h1 : out = ![false] <;> by_cases h2 : inp = ![false] <;> simp [h1, h2]

/-! ## The `src/unitary.rs` suite

Rust's own semantics tests. Each `circuits_equiv(&a, &b, 1e-10)` becomes `Equivalent`, and
each `!circuits_equiv(..)` becomes a refutation — those are what rule out a semantics that
equates too much.

The one-wire rotation tests all reduce to the same modular arithmetic: a diagonal gate is
`diagRun k` for its `k·π/4`, and runs add mod 8. -/

/-- Two diagonal gates on one wire merge. -/
theorem two_diag {n m : Nat} {q : Qubit} {g₁ g₂ : Gate} {k₁ k₂ : Nat}
    (h₁ : diagonalK g₁ q = some (some k₁)) (h₂ : diagonalK g₂ q = some (some k₂)) :
    Equivalent n m [g₁, g₂] (diagRun ((k₁ + k₂) % 8) q) := by
  have e₁ : Equivalent n m ([g₁] ++ [g₂]) (diagRun k₁ q ++ [g₂]) :=
    Equivalent.append_right _ (Equivalent.diagonal_gate h₁)
  have e₂ : Equivalent n m (diagRun k₁ q ++ [g₂]) (diagRun k₁ q ++ diagRun k₂ q) :=
    Equivalent.append_left _ (Equivalent.diagonal_gate h₂)
  exact Equivalent.trans (by simpa using e₁)
    (Equivalent.trans e₂ (Equivalent.diagRun_append k₁ k₂ q (diagonalK_lt h₁) (diagonalK_lt h₂)))

/-- A repeated diagonal gate. -/
theorem diag_replicate {n m : Nat} {q : Qubit} {g : Gate} {k : Nat}
    (h : diagonalK g q = some (some k)) :
    ∀ j : Nat, Equivalent n m (List.replicate j g) (diagRun ((j * k) % 8) q) := by
  intro j
  induction j with
  | zero => simpa [diagRun] using Equivalent.refl n m []
  | succ j ih =>
      have e₁ : Equivalent n m (g :: List.replicate j g) ([g] ++ diagRun ((j * k) % 8) q) := by
        simpa using Equivalent.append_left [g] ih
      have e₂ : Equivalent n m ([g] ++ diagRun ((j * k) % 8) q)
          (diagRun k q ++ diagRun ((j * k) % 8) q) :=
        Equivalent.append_right _ (Equivalent.diagonal_gate h)
      have e₃ : Equivalent n m (diagRun k q ++ diagRun ((j * k) % 8) q)
          (diagRun ((k + (j * k) % 8) % 8) q) :=
        Equivalent.diagRun_append _ _ q (diagonalK_lt h) (Nat.mod_lt _ (by norm_num))
      have harith : (k + (j * k) % 8) % 8 = ((j + 1) * k) % 8 := by
        conv_rhs => rw [Nat.add_mul, one_mul, Nat.add_comm]
        omega
      rw [List.replicate_succ]
      rw [← harith]
      exact Equivalent.trans e₁ (Equivalent.trans e₂ e₃)


/-! ### One-wire rotation identities -/

private theorem dk_t (q : Qubit) : diagonalK (Gate.t q) q = some (some 1) := by simp [diagonalK]
private theorem dk_tdg (q : Qubit) : diagonalK (Gate.tdg q) q = some (some 7) := by
  simp [diagonalK]
private theorem dk_s (q : Qubit) : diagonalK (Gate.s q) q = some (some 2) := by simp [diagonalK]
private theorem dk_sdg (q : Qubit) : diagonalK (Gate.sdg q) q = some (some 6) := by
  simp [diagonalK]
private theorem dk_z (q : Qubit) : diagonalK (Gate.z q) q = some (some 4) := by simp [diagonalK]
private theorem dk_rz (q : Qubit) {θ : ℚ} {k : Nat} (h : classifyQuarterPi θ = some k) :
    diagonalK (Gate.rz θ q) q = some (some k) := by simp [diagonalK, h]

/-- `x_squared_is_identity` -/
theorem x_squared_is_identity {n m : Nat} (q : Qubit) : Equivalent n m [Gate.x q, Gate.x q] [] :=
  Equivalent.selfInverse_pair _ trivial rfl

/-- `h_squared_is_identity` -/
theorem h_squared_is_identity {n m : Nat} (q : Qubit) : Equivalent n m [Gate.h q, Gate.h q] [] :=
  Equivalent.selfInverse_pair _ trivial rfl

/-- `z_is_self_inverse`, `z_z_is_identity` -/
theorem z_z_is_identity {n m : Nat} (q : Qubit) : Equivalent n m [Gate.z q, Gate.z q] [] :=
  Equivalent.selfInverse_pair _ trivial rfl

/-- `t_squared_is_s` -/
theorem t_squared_is_s {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.t q, Gate.t q] [Gate.s q] := by
  have := two_diag (n := n) (m := m) (dk_t q) (dk_t q)
  norm_num [diagRun] at this
  exact this

/-- `t_tdg_is_identity` -/
theorem t_tdg_is_identity {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.t q, Gate.tdg q] [] := by
  have := two_diag (n := n) (m := m) (dk_t q) (dk_tdg q)
  norm_num [diagRun] at this
  exact this

/-- `s_sdg_is_identity` -/
theorem s_sdg_is_identity {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.s q, Gate.sdg q] [] := by
  have := two_diag (n := n) (m := m) (dk_s q) (dk_sdg q)
  norm_num [diagRun] at this
  exact this

/-- `sdg_s_is_identity` -/
theorem sdg_s_is_identity {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.sdg q, Gate.s q] [] := by
  have := two_diag (n := n) (m := m) (dk_sdg q) (dk_s q)
  norm_num [diagRun] at this
  exact this

/-- `z_is_s_squared` -/
theorem z_is_s_squared {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.z q] [Gate.s q, Gate.s q] := by
  have := two_diag (n := n) (m := m) (dk_s q) (dk_s q)
  norm_num [diagRun] at this
  exact this.symm

/-- `sdg_sdg_is_z` -/
theorem sdg_sdg_is_z {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.sdg q, Gate.sdg q] [Gate.z q] := by
  have := two_diag (n := n) (m := m) (dk_sdg q) (dk_sdg q)
  norm_num [diagRun] at this
  exact this

/-- `sdg_t_is_tdg` -/
theorem sdg_t_is_tdg {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.sdg q, Gate.t q] [Gate.tdg q] := by
  have := two_diag (n := n) (m := m) (dk_sdg q) (dk_t q)
  norm_num [diagRun] at this
  exact this

/-- `s_tdg_is_t` -/
theorem s_tdg_is_t {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.s q, Gate.tdg q] [Gate.t q] := by
  have := two_diag (n := n) (m := m) (dk_s q) (dk_tdg q)
  norm_num [diagRun] at this
  exact this

/-- `s_three_times_is_sdg` -/
theorem s_three_times_is_sdg {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.s q, Gate.s q, Gate.s q] [Gate.sdg q] := by
  have := diag_replicate (n := n) (m := m) (dk_s q) 3
  norm_num [diagRun, List.replicate] at this
  exact this

/-- `sdg_three_times_is_s` -/
theorem sdg_three_times_is_s {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.sdg q, Gate.sdg q, Gate.sdg q] [Gate.s q] := by
  have := diag_replicate (n := n) (m := m) (dk_sdg q) 3
  norm_num [diagRun, List.replicate] at this
  exact this

/-- `six_t_is_sdg` -/
theorem six_t_is_sdg {n m : Nat} (q : Qubit) :
    Equivalent n m (List.replicate 6 (Gate.t q)) [Gate.sdg q] := by
  have := diag_replicate (n := n) (m := m) (dk_t q) 6
  norm_num [diagRun] at this
  exact this

/-- `s_is_rz_pi_over_2` -/
theorem s_is_rz_pi_over_2 {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.s q] [Gate.rz (1/2) q] := by
  have hcl : classifyQuarterPi (1/2) = some 2 := by
    norm_num [classifyQuarterPi]
    rfl
  have := Equivalent.diagonal_gate (n := n) (m := m) (dk_rz q hcl)
  norm_num [diagRun] at this
  exact this.symm

/-- `z_is_rz_pi` -/
theorem z_is_rz_pi {n m : Nat} (q : Qubit) : Equivalent n m [Gate.z q] [Gate.rz 1 q] := by
  have hcl : classifyQuarterPi 1 = some 4 := by
    norm_num [classifyQuarterPi]
    rfl
  have := Equivalent.diagonal_gate (n := n) (m := m) (dk_rz q hcl)
  norm_num [diagRun] at this
  exact this.symm

/-- `sdg_is_rz_neg_half_pi` -/
theorem sdg_is_rz_neg_half_pi {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.sdg q] [Gate.rz (-(1/2)) q] := by
  have hcl : classifyQuarterPi (-(1/2)) = some 6 := by
    norm_num [classifyQuarterPi]
    rfl
  have := Equivalent.diagonal_gate (n := n) (m := m) (dk_rz q hcl)
  norm_num [diagRun] at this
  exact this.symm

/-- `z_t_is_five_quarter_pi` -/
theorem z_t_is_five_quarter_pi {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.z q, Gate.t q] [Gate.rz (5/4) q] := by
  have hcl : classifyQuarterPi (5/4) = some 5 := by
    norm_num [classifyQuarterPi]
    rfl
  have h₁ := two_diag (n := n) (m := m) (dk_z q) (dk_t q)
  have h₂ := Equivalent.diagonal_gate (n := n) (m := m) (dk_rz q hcl)
  norm_num [diagRun] at h₁ h₂
  exact Equivalent.trans h₁ h₂.symm

/-- `hxh_is_z` and `hzh_is_x` -/
theorem hxh_is_z {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.h q, Gate.x q, Gate.h q] [Gate.z q] := Equivalent.hxh q

theorem hzh_is_x {n m : Nat} (q : Qubit) :
    Equivalent n m [Gate.h q, Gate.z q, Gate.h q] [Gate.x q] := Equivalent.hzh q

/-! ### Multi-wire identities -/

/-- `double_ccx_is_identity` -/
theorem double_ccx_is_identity {n m : Nat} :
    Equivalent n m [Gate.ccx 0 1 2, Gate.ccx 0 1 2] [] :=
  Equivalent.selfInverse_pair _ (by norm_num [Gate.Wf]) rfl

/-- `cz_squared_is_identity` — with the operands swapped, as the Rust test writes it. -/
theorem cz_squared_is_identity {n m : Nat} : Equivalent n m [Gate.cz 0 1, Gate.cz 1 0] [] :=
  Equivalent.gatesEqual_cancel (by decide) (by norm_num [Gate.Wf]) rfl

/-- `ccz_squared_is_identity` -/
theorem ccz_squared_is_identity {n m : Nat} :
    Equivalent n m [Gate.ccz 0 1 2, Gate.ccz 2 0 1] [] :=
  Equivalent.gatesEqual_cancel (by decide) (by norm_num [Gate.Wf]) rfl

/-- `cz_is_symmetric` -/
theorem cz_is_symmetric {n m : Nat} (a b : Qubit) :
    Equivalent n m [Gate.cz a b] [Gate.cz b a] :=
  equivalent_of_unitary_eq (by intro g hg; fin_cases hg; rfl) (by intro g hg; fin_cases hg; rfl)
    (by simp [unitary_cons, gateUnitary_cz_comm n a b])

/-- `ccz_is_symmetric_in_all_operands` — the five nontrivial reorderings. -/
theorem ccz_is_symmetric (a b c a' b' c' : Qubit)
    (h : gatesEqual (Gate.ccz a b c) (Gate.ccz a' b' c') = true) {n m : Nat} :
    Equivalent n m [Gate.ccz a b c] [Gate.ccz a' b' c'] :=
  equivalent_of_unitary_eq (by intro g hg; fin_cases hg; rfl) (by intro g hg; fin_cases hg; rfl)
    (by simp [unitary_cons, gateUnitary_of_gatesEqual n h])

example {n m : Nat} : Equivalent n m [Gate.ccz 0 1 2] [Gate.ccz 0 2 1] :=
  ccz_is_symmetric 0 1 2 0 2 1 (by decide)
example {n m : Nat} : Equivalent n m [Gate.ccz 0 1 2] [Gate.ccz 1 0 2] :=
  ccz_is_symmetric 0 1 2 1 0 2 (by decide)
example {n m : Nat} : Equivalent n m [Gate.ccz 0 1 2] [Gate.ccz 1 2 0] :=
  ccz_is_symmetric 0 1 2 1 2 0 (by decide)
example {n m : Nat} : Equivalent n m [Gate.ccz 0 1 2] [Gate.ccz 2 0 1] :=
  ccz_is_symmetric 0 1 2 2 0 1 (by decide)
example {n m : Nat} : Equivalent n m [Gate.ccz 0 1 2] [Gate.ccz 2 1 0] :=
  ccz_is_symmetric 0 1 2 2 1 0 (by decide)

/-- `four_qubit_parallel_cnots` -/
theorem four_qubit_parallel_cnots {n m : Nat} :
    Equivalent n m [Gate.cnot 0 1, Gate.cnot 2 3] [Gate.cnot 2 3, Gate.cnot 0 1] :=
  Equivalent.swap_of_disjoint rfl (disjoint_support_of (by decide))

/-- `cnot_commute_disjoint_3qubit`: two `CNOT`s onto the same target commute. -/
theorem cnot_commute_same_target {n m : Nat} :
    Equivalent n m [Gate.cnot 0 1, Gate.cnot 2 1] [Gate.cnot 2 1, Gate.cnot 0 1] :=
  Equivalent.pairCommutes_swap rfl (by decide)

/-- `z_commutes_with_cnot_control` -/
theorem z_commutes_with_cnot_control {n m : Nat} :
    Equivalent n m [Gate.z 0, Gate.cnot 0 1] [Gate.cnot 0 1, Gate.z 0] :=
  (Equivalent.pairCommutes_swap (p := Gate.cnot 0 1) rfl (by decide)).symm

/-- `cz_commutes_with_cnot_when_target_is_outside_operands` -/
theorem cz_commutes_with_cnot_outside {n m : Nat} :
    Equivalent n m [Gate.cz 0 1, Gate.cnot 0 2] [Gate.cnot 0 2, Gate.cz 0 1] :=
  Equivalent.pairCommutes_swap rfl (by decide)

/-- A list of pairwise-disjoint unitary gates may be reversed. -/
theorem Equivalent.reverse_of_pairwise {n m : Nat} : ∀ {gs : List Gate},
    (∀ g ∈ gs, g.isUnitary = true) →
    gs.Pairwise (fun g g' => Wires.Disjoint g.support g'.support) →
    Equivalent n m gs gs.reverse := by
  intro gs
  induction gs with
  | nil => intro _ _; exact Equivalent.refl _ _ _
  | cons g gs ih =>
      intro hu hp
      have hu' : ∀ g' ∈ gs, g'.isUnitary = true := fun g' hg' => hu g' (by simp [hg'])
      have hhead : ∀ g' ∈ gs, Wires.Disjoint g.support g'.support :=
        (List.pairwise_cons.mp hp).1
      have htail := (List.pairwise_cons.mp hp).2
      have step₁ : Equivalent n m (g :: gs) (g :: gs.reverse) :=
        by simpa using Equivalent.append_left [g] (ih hu' htail)
      have step₂ : Equivalent n m (g :: gs.reverse) (gs.reverse ++ [g]) :=
        Equivalent.move_back (hu g (by simp)) gs.reverse
          (fun g' hg' => hhead g' (List.mem_reverse.mp hg'))
      simpa using Equivalent.trans step₁ step₂

/-- `four_qubit_h_layer_order` -/
theorem four_qubit_h_layer_order {n m : Nat} :
    Equivalent n m [Gate.h 0, Gate.h 1, Gate.h 2, Gate.h 3]
      [Gate.h 3, Gate.h 2, Gate.h 1, Gate.h 0] := by
  have hd : ∀ a b : Qubit, a ≠ b → Wires.Disjoint (Gate.h a).support (Gate.h b).support := by
    intro a b hab
    exact disjoint_support_of (by
      intro q hq q' hq'
      simp only [Gate.qubitsOf, List.mem_singleton] at hq hq'
      subst hq; subst hq'; exact hab)
  have hpw : List.Pairwise (fun g g' => Wires.Disjoint g.support g'.support)
      [Gate.h 0, Gate.h 1, Gate.h 2, Gate.h 3] := by
    refine List.Pairwise.cons (fun g hg => ?_) (List.Pairwise.cons (fun g hg => ?_)
      (List.Pairwise.cons (fun g hg => ?_) (List.Pairwise.cons (fun g hg => ?_)
        List.Pairwise.nil)))
    · simp only [List.mem_cons, List.not_mem_nil, or_false] at hg
      rcases hg with rfl | rfl | rfl <;> exact hd _ _ (by decide)
    · simp only [List.mem_cons, List.not_mem_nil, or_false] at hg
      rcases hg with rfl | rfl <;> exact hd _ _ (by decide)
    · simp only [List.mem_cons, List.not_mem_nil, or_false] at hg
      rcases hg with rfl
      exact hd _ _ (by decide)
    · simp at hg
  exact Equivalent.reverse_of_pairwise (by intro g hg; fin_cases hg <;> rfl) hpw

/-! ### The negative tests

`!circuits_equiv(..)` in Rust. These are what rule out a semantics that equates too much —
a `denote` that returned a constant would satisfy every positive test above. -/

/-- To refute an equivalence it is enough to exhibit one entry of one density matrix where
the two circuits disagree. -/
theorem not_equivalent_of_entry {n : Nat} {gs hs : List Gate}
    (hg : ∀ g ∈ gs, g.isUnitary = true) (hh : ∀ g ∈ hs, g.isUnitary = true)
    (b b' out inp : Basis n)
    (hne : unitary n gs out b * star (unitary n gs inp b') ≠
      unitary n hs out b * star (unitary n hs inp b')) : ¬ Equivalent n 0 gs hs := by
  intro heq
  have h1 := heq (fun _ => ketBra b b')
  rw [denote_eq_conj_unitary gs hg, denote_eq_conj_unitary hs hh] at h1
  have h2 := congrFun h1 (fun i => i.elim0)
  have h3 := congrFun (congrFun h2 out) inp
  rw [conj_ketBra, conj_ketBra] at h3
  exact hne h3

/-- Two permutation circuits sending a basis state to different places are not equivalent. -/
theorem not_equivalent_of_perm {n : Nat} {σ τ : Basis n → Basis n} {gs hs : List Gate}
    (hg : ∀ g ∈ gs, g.isUnitary = true) (hh : ∀ g ∈ hs, g.isUnitary = true)
    (hgs : unitary n gs = permMatrix σ) (hhs : unitary n hs = permMatrix τ)
    (b : Basis n) (hne : σ b ≠ τ b) : ¬ Equivalent n 0 gs hs := by
  refine not_equivalent_of_entry hg hh b b (σ b) (σ b) ?_
  rw [hgs, hhs]
  simp only [permMatrix, if_neg hne]
  norm_num

/-- `not_equiv_h_vs_x` -/
theorem not_equiv_h_vs_x : ¬ Equivalent 1 0 [Gate.h 0] [Gate.x 0] := by
  have hsqrt : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := ofReal_sqrt_two_ne_zero
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := ofReal_sqrt_two_mul_self
  refine not_equivalent_of_entry (by intro g hg; fin_cases hg; rfl)
    (by intro g hg; fin_cases hg; rfl) ![false] ![false] ![false] ![false] ?_
  have hh : unitary 1 [Gate.h 0] ![false] ![false] = 1 / ((Real.sqrt 2 : ℝ) : ℂ) := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary]
    rw [embed1_apply_one]
    simp [h2]
  have hx : unitary 1 [Gate.x 0] ![false] ![false] = 0 := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary]
    rw [embed1_apply_one]
    simp [x2]
  rw [hh, hx]
  have hstar : star ((1 : ℂ) / ((Real.sqrt 2 : ℝ) : ℂ)) = 1 / ((Real.sqrt 2 : ℝ) : ℂ) := by
    simp [Complex.conj_ofReal]
  rw [hstar, div_mul_div_comm, one_mul, hsq]
  norm_num

/-- `hsdgh_is_not_s` -/
theorem hsdgh_is_not_s : ¬ Equivalent 1 0 [Gate.h 0, Gate.sdg 0, Gate.h 0] [Gate.s 0] := by
  have hsqrt : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := ofReal_sqrt_two_ne_zero
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := ofReal_sqrt_two_mul_self
  refine not_equivalent_of_entry (by intro g hg; fin_cases hg <;> rfl)
    (by intro g hg; fin_cases hg; rfl) ![false] ![false] ![false] ![false] ?_
  have ha : unitary 1 [Gate.h 0, Gate.sdg 0, Gate.h 0] ![false] ![false] =
      (1 + -Complex.I) / (((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ)) := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary, embed1_mul]
    rw [embed1_apply_one]
    simp [mul2, h2, diag2, ep_neg_half]
    ring
  have hb : unitary 1 [Gate.s 0] ![false] ![false] = 1 := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary]
    rw [embed1_apply_one]
    simp [diag2]
  rw [ha, hb, hsq]
  have hstar : star ((1 + -Complex.I) / 2) = (1 + Complex.I) / 2 := by
    simp
  rw [hstar]
  intro hc
  have h2' : (1 + -Complex.I) * (1 + Complex.I) = 2 := by
    have := Complex.I_sq
    ring_nf
    rw [Complex.I_sq]
    ring
  field_simp at hc
  rw [h2'] at hc
  norm_num at hc

/-- `cnot_no_commute_overlapping`: `CNOT(0,1)` and `CNOT(1,2)` do not commute. -/
theorem cnot_no_commute_overlapping :
    ¬ Equivalent 3 0 [Gate.cnot 0 1, Gate.cnot 1 2] [Gate.cnot 1 2, Gate.cnot 0 1] := by
  refine not_equivalent_of_perm
    (σ := (fun b : Basis 3 => b.set 2 (b.get 2 != b.get 1)) ∘
      (fun b : Basis 3 => b.set 1 (b.get 1 != b.get 0)))
    (τ := (fun b : Basis 3 => b.set 1 (b.get 1 != b.get 0)) ∘
      (fun b : Basis 3 => b.set 2 (b.get 2 != b.get 1)))
    (by intro g hg; fin_cases hg <;> rfl) (by intro g hg; fin_cases hg <;> rfl)
    (by simp [unitary_cons, gateUnitary, permMatrix_mul])
    (by simp [unitary_cons, gateUnitary, permMatrix_mul])
    ![true, false, false] ?_
  intro hc
  have h2 := congrFun hc 2
  simp [Function.comp_apply, Basis.set, Basis.get] at h2

/-- `four_qubit_cnot_ladder_not_reversible` -/
theorem four_qubit_cnot_ladder_not_reversible :
    ¬ Equivalent 4 0 [Gate.cnot 0 1, Gate.cnot 1 2, Gate.cnot 2 3]
      [Gate.cnot 2 3, Gate.cnot 1 2, Gate.cnot 0 1] := by
  refine not_equivalent_of_perm
    (σ := ((fun b : Basis 4 => b.set 3 (b.get 3 != b.get 2)) ∘
      (fun b : Basis 4 => b.set 2 (b.get 2 != b.get 1))) ∘
      (fun b : Basis 4 => b.set 1 (b.get 1 != b.get 0)))
    (τ := ((fun b : Basis 4 => b.set 1 (b.get 1 != b.get 0)) ∘
      (fun b : Basis 4 => b.set 2 (b.get 2 != b.get 1))) ∘
      (fun b : Basis 4 => b.set 3 (b.get 3 != b.get 2)))
    (by intro g hg; fin_cases hg <;> rfl) (by intro g hg; fin_cases hg <;> rfl)
    (by simp [unitary_cons, gateUnitary, permMatrix_mul])
    (by simp [unitary_cons, gateUnitary, permMatrix_mul])
    ![true, false, false, false] ?_
  intro hc
  have h2 := congrFun hc 2
  simp [Function.comp_apply, Basis.set, Basis.get] at h2

theorem gateUnitary_cnot_perm (n : Nat) (c t : Qubit) :
    gateUnitary n (.cnot c t) = permMatrix (fun b : Basis n => b.set t (b.get t != b.get c)) :=
  rfl

theorem gateUnitary_cz_phase (n : Nat) (c t : Qubit) :
    gateUnitary n (.cz c t) =
      phaseMatrix (fun b : Basis n => if b.get c && b.get t then (-1 : ℂ) else 1) := rfl

private theorem cnot_involutive_two : ∀ b : Basis 2,
    (fun s : Basis 2 => s.set 1 (s.get 1 != s.get 0))
      ((fun s : Basis 2 => s.set 1 (s.get 1 != s.get 0)) b) = b := by
  intro b
  simp only
  refine basis_two_ext ?_ ?_
  · rw [cnot_perm_zero, cnot_perm_zero]
  · rw [cnot_perm_one, cnot_perm_one, cnot_perm_zero]
    cases b 0 <;> cases b 1 <;> rfl

/-- `z_does_not_commute_with_cnot_target`: a `Z` on the *target* does not commute with the
`CNOT`, while `z_commutes_with_cnot_control` above shows one on the control does. -/
theorem z_does_not_commute_with_cnot_target :
    ¬ Equivalent 2 0 [Gate.z 1, Gate.cnot 0 1] [Gate.cnot 0 1, Gate.z 1] := by
  have hz : gateUnitary 2 (.z 1) =
      phaseMatrix (fun s : Basis 2 => if s.get 1 then ep 1 else 1) := by
    rw [gateUnitary, embed1_diag2_eq_phaseMatrix _ _ _ (by norm_num : (1:Nat) < 2)]
  refine not_equivalent_of_entry (by intro g hg; fin_cases hg <;> rfl)
    (by intro g hg; fin_cases hg <;> rfl)
    ![true, false] ![false, false] ![true, true] ![false, false] ?_
  have hσ₁ : Basis.set (![true, true] : Basis 2) 1
      (Basis.get (![true, true] : Basis 2) 1 != Basis.get (![true, true] : Basis 2) 0)
      = ![true, false] := by rw [cnot_perm_apply]; rfl
  have hσ₂ : Basis.set (![false, false] : Basis 2) 1
      (Basis.get (![false, false] : Basis 2) 1 != Basis.get (![false, false] : Basis 2) 0)
      = ![false, false] := by rw [cnot_perm_apply]; rfl
  have hσ₃ : Basis.set (![true, false] : Basis 2) 1
      (Basis.get (![true, false] : Basis 2) 1 != Basis.get (![true, false] : Basis 2) 0)
      = ![true, true] := by rw [cnot_perm_apply]; rfl
  have ha : ∀ out b : Basis 2, unitary 2 [Gate.z 1, Gate.cnot 0 1] out b =
      if Basis.set out 1 (Basis.get out 1 != Basis.get out 0) = b then
        (if Basis.get b 1 then ep 1 else 1) else 0 := by
    intro out b
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [hz, gateUnitary_cnot_perm, permMatrix_mul_apply cnot_involutive_two]
    simp only [phaseMatrix]
  have hb : ∀ out b : Basis 2, unitary 2 [Gate.cnot 0 1, Gate.z 1] out b =
      (if Basis.get out 1 then ep 1 else 1) *
        (if out = Basis.set b 1 (Basis.get b 1 != Basis.get b 0) then 1 else 0) := by
    intro out b
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [hz, gateUnitary_cnot_perm, phaseMatrix_mul_apply]
    simp only [permMatrix]
  rw [ha, ha, hb, hb, hσ₁, hσ₂, hσ₃]
  have hg1 : Basis.get (![true, false] : Basis 2) 1 = false := by simp [Basis.get]
  have hg2 : Basis.get (![false, false] : Basis 2) 1 = false := by simp [Basis.get]
  have hg3 : Basis.get (![true, true] : Basis 2) 1 = true := by simp [Basis.get]
  rw [hg1, hg2, hg3, ep_one]
  norm_num

/-- `cz_does_not_commute_with_x_on_operand` -/
theorem cz_does_not_commute_with_x_on_operand :
    ¬ Equivalent 2 0 [Gate.x 0, Gate.cz 0 1] [Gate.cz 0 1, Gate.x 0] := by
  have hx : gateUnitary 2 (.x 0) = permMatrix (fun b : Basis 2 => b.set 0 (!b.get 0)) := by
    rw [gateUnitary, embed1_x2_eq_perm]
  have hxinv : ∀ b : Basis 2,
      (fun s : Basis 2 => s.set 0 (!s.get 0)) ((fun s : Basis 2 => s.set 0 (!s.get 0)) b) = b := by
    intro b
    simp only
    refine basis_two_ext ?_ ?_ <;> simp [Basis.set, Basis.get]
  refine not_equivalent_of_entry (by intro g hg; fin_cases hg <;> rfl)
    (by intro g hg; fin_cases hg <;> rfl)
    ![false, true] ![false, false] ![true, true] ![true, false] ?_
  have hs1 : Basis.set (![true, true] : Basis 2) 0 (!Basis.get (![true, true] : Basis 2) 0)
      = ![false, true] := basis_two_ext (by simp [Basis.set, Basis.get]) (by simp [Basis.set])
  have hs2 : Basis.set (![true, false] : Basis 2) 0 (!Basis.get (![true, false] : Basis 2) 0)
      = ![false, false] := basis_two_ext (by simp [Basis.set, Basis.get]) (by simp [Basis.set])
  have ha : ∀ out b : Basis 2, unitary 2 [Gate.x 0, Gate.cz 0 1] out b =
      (if Basis.get out 0 && Basis.get out 1 then (-1 : ℂ) else 1) *
        (if out = Basis.set b 0 (!Basis.get b 0) then 1 else 0) := by
    intro out b
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [hx, gateUnitary_cz_phase, phaseMatrix_mul_apply]
    simp only [permMatrix]
  have hb : ∀ out b : Basis 2, unitary 2 [Gate.cz 0 1, Gate.x 0] out b =
      if Basis.set out 0 (!Basis.get out 0) = b then
        (if Basis.get b 0 && Basis.get b 1 then (-1 : ℂ) else 1) else 0 := by
    intro out b
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [hx, gateUnitary_cz_phase, permMatrix_mul_apply hxinv]
    simp only [phaseMatrix]
  have hs3 : Basis.set (![false, false] : Basis 2) 0 true = ![true, false] :=
    basis_two_ext (by simp [Basis.set]) (by simp [Basis.set])
  have hs4 : Basis.set (![false, true] : Basis 2) 0 true = ![true, true] :=
    basis_two_ext (by simp [Basis.set]) (by simp [Basis.set])
  rw [ha, ha, hb, hb, hs1, hs2]
  norm_num [Basis.get, hs3, hs4]

/-! ### `CZ` by conjugating `CNOT`, and the remaining commutation tests -/

/-- `identity_is_identity` -/
theorem identity_is_identity (n : Nat) : unitary n [] = 1 := rfl

/-- `cz_equals_h_cnot_h`: `H` on the target turns `CNOT` into `CZ`. This is the check that
exercises composition *across* wires with a non-permutation gate. -/
theorem cz_equals_h_cnot_h :
    unitary 2 [Gate.h 1, Gate.cnot 0 1, Gate.h 1] = gateUnitary 2 (.cz 0 1) := by
  have hsqrt : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := ofReal_sqrt_two_ne_zero
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := ofReal_sqrt_two_mul_self
  funext out inp
  have hentry : unitary 2 [Gate.h 1, Gate.cnot 0 1, Gate.h 1] out inp =
      ∑ v : Bool, h2 (out 1) v *
        (if out 0 = inp 0 then h2 (v != out 0) (inp 1) else 0) := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [show gateUnitary 2 (Gate.h 1) = embed1 2 h2 1 from rfl, Matrix.mul_assoc,
      embed1_mul_apply h2 1 (by norm_num) _ out inp]
    refine Finset.sum_congr rfl fun v _ => ?_
    congr 1
    rw [gateUnitary_cnot_perm, permMatrix_mul_apply cnot_involutive_two, embed1_apply_two_one]
    have h0 : (Basis.set (out.set 1 v) 1
        (Basis.get (out.set 1 v) 1 != Basis.get (out.set 1 v) 0)) 0 = out 0 := by
      rw [cnot_perm_zero]
      simp [Basis.set]
    have h1 : (Basis.set (out.set 1 v) 1
        (Basis.get (out.set 1 v) 1 != Basis.get (out.set 1 v) 0)) 1 = (v != out 0) := by
      rw [cnot_perm_one]
      simp [Basis.set]
    rw [h0, h1]
  have g0 : ∀ b : Basis 2, Basis.get b 0 = b 0 := fun b => by simp [Basis.get]
  have g1 : ∀ b : Basis 2, Basis.get b 1 = b 1 := fun b => by simp [Basis.get]
  rw [hentry, gateUnitary_cz_phase]
  simp only [phaseMatrix, Fintype.sum_bool, g0, g1]
  by_cases h0 : out 0 = inp 0
  · rw [if_pos h0]
    by_cases hout : out = inp
    · subst hout
      rw [if_pos rfl]
      cases hb0 : out 0 <;> cases hb1 : out 1 <;>
        · try rw [hb0, hb1]
          simp only [h2]
          norm_num
          try field_simp
          try exact hsq.symm
          try (rw [pow_two, hsq])
          try rw [hsq]
          try ring
    · rw [if_neg hout]
      have h1 : out 1 ≠ inp 1 := by
        intro hc
        exact hout (basis_two_ext h0 hc)
      rw [← h0]
      cases hb0 : out 0 <;> cases hb1 : out 1 <;> cases hi1 : inp 1 <;>
        · first
            | (exfalso; rw [hb1, hi1] at h1; exact h1 rfl)
            | (try rw [hb0, hb1, hi1]
               simp only [h2]
               norm_num
               try field_simp
               try exact hsq.symm
               try (rw [pow_two, hsq])
               try rw [hsq]
               try ring)
  · simp only [if_neg h0, if_neg (show ¬ out = inp from fun hc => h0 (by rw [hc])),
      mul_zero, add_zero]

/-- `cz_commutes_with_diagonals_on_both_operands` -/
theorem cz_commutes_with_diagonals {n m : Nat} :
    Equivalent n m [Gate.t 0, Gate.sdg 1, Gate.cz 0 1] [Gate.cz 1 0, Gate.t 0, Gate.sdg 1] := by
  have swap₁ : Equivalent n m [Gate.sdg 1, Gate.cz 0 1] [Gate.cz 0 1, Gate.sdg 1] :=
    (Equivalent.pairCommutes_swap (p := Gate.cz 0 1) rfl (by decide)).symm
  have swap₂ : Equivalent n m [Gate.t 0, Gate.cz 0 1] [Gate.cz 0 1, Gate.t 0] :=
    (Equivalent.pairCommutes_swap (p := Gate.cz 0 1) rfl (by decide)).symm
  have step₁ : Equivalent n m [Gate.t 0, Gate.sdg 1, Gate.cz 0 1]
      [Gate.t 0, Gate.cz 0 1, Gate.sdg 1] := by
    simpa using Equivalent.append_left [Gate.t 0] swap₁
  have step₂ : Equivalent n m [Gate.t 0, Gate.cz 0 1, Gate.sdg 1]
      [Gate.cz 0 1, Gate.t 0, Gate.sdg 1] := by
    simpa using Equivalent.append_right [Gate.sdg 1] swap₂
  have step₃ : Equivalent n m [Gate.cz 0 1, Gate.t 0, Gate.sdg 1]
      [Gate.cz 1 0, Gate.t 0, Gate.sdg 1] := by
    simpa using Equivalent.append_right [Gate.t 0, Gate.sdg 1] (cz_is_symmetric (n := n) 0 1)
  exact Equivalent.trans step₁ (Equivalent.trans step₂ step₃)

/-- `cnot_swap_decomposition`: three `CNOT`s make a `SWAP`. -/
theorem cnot_swap_decomposition :
    unitary 2 [Gate.cnot 0 1, Gate.cnot 1 0, Gate.cnot 0 1] =
      permMatrix (fun b : Basis 2 => ![b 1, b 0]) := by
  simp only [unitary_cons, unitary_nil, Matrix.one_mul, gateUnitary_cnot_perm, permMatrix_mul]
  congr 1
  funext b
  refine basis_two_ext ?_ ?_ <;>
    simp [Function.comp_apply, Basis.set, Basis.get]
  all_goals (cases b 0 <;> cases b 1 <;> rfl)

/-- `cz_does_not_commute_with_cnot_targeting_operand` -/
theorem cz_does_not_commute_with_cnot_targeting_operand :
    ¬ Equivalent 3 0 [Gate.cz 0 1, Gate.cnot 2 1] [Gate.cnot 2 1, Gate.cz 0 1] := by
  have hinv : ∀ b : Basis 3,
      (fun s : Basis 3 => s.set 1 (s.get 1 != s.get 2))
        ((fun s : Basis 3 => s.set 1 (s.get 1 != s.get 2)) b) = b := by
    intro b
    funext r
    by_cases h1 : (r : Nat) = 1
    · have hr : r = ⟨1, by norm_num⟩ := Fin.ext h1
      subst hr
      simp [Basis.set, Basis.get]
    · simp [Basis.set, h1]
  refine not_equivalent_of_entry (by intro g hg; fin_cases hg <;> rfl)
    (by intro g hg; fin_cases hg <;> rfl)
    ![true, false, true] ![false, false, false] ![true, true, true] ![false, false, false] ?_
  have ha : ∀ out b : Basis 3, unitary 3 [Gate.cz 0 1, Gate.cnot 2 1] out b =
      if Basis.set out 1 (Basis.get out 1 != Basis.get out 2) = b then
        (if Basis.get b 0 && Basis.get b 1 then (-1 : ℂ) else 1) else 0 := by
    intro out b
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [gateUnitary_cz_phase, gateUnitary_cnot_perm, permMatrix_mul_apply hinv]
    simp only [phaseMatrix]
  have hb : ∀ out b : Basis 3, unitary 3 [Gate.cnot 2 1, Gate.cz 0 1] out b =
      (if Basis.get out 0 && Basis.get out 1 then (-1 : ℂ) else 1) *
        (if out = Basis.set b 1 (Basis.get b 1 != Basis.get b 2) then 1 else 0) := by
    intro out b
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [gateUnitary_cz_phase, gateUnitary_cnot_perm, phaseMatrix_mul_apply]
    simp only [permMatrix]
  rw [ha, ha, hb, hb]
  have e1 : Basis.set (![true, true, true] : Basis 3) 1
      (Basis.get (![true, true, true] : Basis 3) 1 !=
        Basis.get (![true, true, true] : Basis 3) 2) = ![true, false, true] := by
    funext r
    match r with
    | 0 => simp [Basis.set]
    | 1 => simp [Basis.set, Basis.get]
    | 2 => simp [Basis.set]
  have e2 : Basis.set (![false, false, false] : Basis 3) 1
      (Basis.get (![false, false, false] : Basis 3) 1 !=
        Basis.get (![false, false, false] : Basis 3) 2) = ![false, false, false] := by
    funext r
    match r with
    | 0 => simp [Basis.set]
    | 1 => simp [Basis.set, Basis.get]
    | 2 => simp [Basis.set]
  have e3 : Basis.set (![true, false, true] : Basis 3) 1
      (Basis.get (![true, false, true] : Basis 3) 1 !=
        Basis.get (![true, false, true] : Basis 3) 2) = ![true, true, true] := by
    funext r
    match r with
    | 0 => simp [Basis.set]
    | 1 => simp [Basis.set, Basis.get]
    | 2 => simp [Basis.set]
  rw [e1, e2, e3]
  norm_num [Basis.get]

/-- Entries of a one-wire gate, at index pairs that agree off the wire. -/
theorem embed1_apply_agree {n : Nat} (M : Bool → Bool → ℂ) {q : Qubit} (hq : q < n)
    (out inp : Basis n) (h : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r) :
    embed1 n M q out inp = M (out ⟨q, hq⟩) (inp ⟨q, hq⟩) := by
  simp only [embed1_apply_of_lt _ hq, if_pos h]

theorem embed1_apply_disagree {n : Nat} (M : Bool → Bool → ℂ) {q : Qubit} (hq : q < n)
    (out inp : Basis n) (h : ¬ ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r) :
    embed1 n M q out inp = 0 := by
  simp only [embed1_apply_of_lt _ hq, if_neg h]

theorem gateUnitary_ccx_perm (n : Nat) (c₁ c₂ t : Qubit) :
    gateUnitary n (.ccx c₁ c₂ t) =
      permMatrix (fun b : Basis n => b.set t (b.get t != (b.get c₁ && b.get c₂))) := rfl

theorem gateUnitary_ccz_phase (n : Nat) (c₁ c₂ t : Qubit) :
    gateUnitary n (.ccz c₁ c₂ t) =
      phaseMatrix (fun b : Basis n =>
        if b.get c₁ && b.get c₂ && b.get t then (-1 : ℂ) else 1) := rfl

/-- `ccz_equals_h_ccx_h`: `H` on the target turns `CCX` into `CCZ` — the three-wire version
of `cz_equals_h_cnot_h`. -/
theorem ccz_equals_h_ccx_h :
    unitary 3 [Gate.h 2, Gate.ccx 0 1 2, Gate.h 2] = gateUnitary 3 (.ccz 0 1 2) := by
  have hsqrt : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := ofReal_sqrt_two_ne_zero
  have hsq : ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := ofReal_sqrt_two_mul_self
  have h2lt : (2 : Nat) < 3 := by norm_num
  have g0 : ∀ b : Basis 3, Basis.get b 0 = b 0 := fun b => by simp [Basis.get]
  have g1 : ∀ b : Basis 3, Basis.get b 1 = b 1 := fun b => by simp [Basis.get]
  have g2 : ∀ b : Basis 3, Basis.get b 2 = b 2 := fun b => by simp [Basis.get]
  have hinv : ∀ b : Basis 3,
      (fun s : Basis 3 => s.set 2 (s.get 2 != (s.get 0 && s.get 1)))
        ((fun s : Basis 3 => s.set 2 (s.get 2 != (s.get 0 && s.get 1))) b) = b := by
    intro b
    funext r
    by_cases hr : (r : Nat) = 2
    · have hrq : r = ⟨2, h2lt⟩ := Fin.ext hr
      subst hrq
      simp [Basis.set, Basis.get]
    · simp [Basis.set, hr]
  funext out inp
  have hentry : unitary 3 [Gate.h 2, Gate.ccx 0 1 2, Gate.h 2] out inp =
      ∑ v : Bool, h2 (out 2) v *
        (if (∀ r : Fin 3, (r : Nat) ≠ 2 → out r = inp r) then
          h2 (v != (out 0 && out 1)) (inp 2) else 0) := by
    simp only [unitary_cons, unitary_nil, Matrix.one_mul]
    rw [show gateUnitary 3 (Gate.h 2) = embed1 3 h2 2 from rfl, Matrix.mul_assoc,
      embed1_mul_apply h2 2 h2lt _ out inp]
    refine Finset.sum_congr rfl fun v _ => ?_
    congr 1
    rw [gateUnitary_ccx_perm, permMatrix_mul_apply hinv]
    set k : Basis 3 := Basis.set (out.set 2 v) 2
      (Basis.get (out.set 2 v) 2 != (Basis.get (out.set 2 v) 0 && Basis.get (out.set 2 v) 1))
      with hk
    have hk0 : k 0 = out 0 := by rw [hk]; simp [Basis.set]
    have hk1 : k 1 = out 1 := by rw [hk]; simp [Basis.set]
    have hk2 : k ⟨2, h2lt⟩ = (v != (out 0 && out 1)) := by
      rw [hk]; simp [Basis.set, Basis.get]
    by_cases hag : ∀ r : Fin 3, (r : Nat) ≠ 2 → out r = inp r
    · have hkag : ∀ r : Fin 3, (r : Nat) ≠ 2 → k r = inp r := by
        intro r hr
        match r with
        | 0 => rw [hk0]; exact hag 0 hr
        | 1 => rw [hk1]; exact hag 1 hr
        | 2 => exact absurd rfl hr
      rw [if_pos hag, embed1_apply_agree h2 h2lt k inp hkag, hk2]
      rfl
    · have hkdis : ¬ ∀ r : Fin 3, (r : Nat) ≠ 2 → k r = inp r := by
        intro hcon
        refine hag fun r hr => ?_
        match r with
        | 0 => rw [← hk0]; exact hcon 0 hr
        | 1 => rw [← hk1]; exact hcon 1 hr
        | 2 => exact absurd rfl hr
      rw [if_neg hag, embed1_apply_disagree h2 h2lt k inp hkdis]
  rw [hentry, gateUnitary_ccz_phase]
  simp only [phaseMatrix, Fintype.sum_bool, g0, g1, g2]
  by_cases hag : ∀ r : Fin 3, (r : Nat) ≠ 2 → out r = inp r
  · simp only [if_pos hag]
    by_cases hout : out = inp
    · subst hout
      rw [if_pos rfl]
      cases hb0 : out 0 <;> cases hb1 : out 1 <;> cases hb2 : out 2 <;>
        · simp only [h2]
          norm_num
          try field_simp
          try exact hsq.symm
          try (rw [pow_two, hsq])
          try rw [hsq]
          try ring
    · rw [if_neg hout]
      have h2ne : out 2 ≠ inp 2 := by
        intro hc
        refine hout (funext fun r => ?_)
        by_cases hr : (r : Nat) = 2
        · have : r = ⟨2, h2lt⟩ := Fin.ext hr
          rw [this]; exact hc
        · exact hag r hr
      cases hb0 : out 0 <;> cases hb1 : out 1 <;> cases hb2 : out 2 <;> cases hi2 : inp 2 <;>
        · first
            | (exfalso; rw [hb2, hi2] at h2ne; exact h2ne rfl)
            | (simp only [h2]
               norm_num
               try field_simp
               try exact hsq.symm
               try (rw [pow_two, hsq])
               try rw [hsq]
               try ring)
  · simp only [if_neg hag, mul_zero, add_zero,
      if_neg (show ¬ out = inp from fun hc => hag fun r _ => by rw [hc])]

/-! ## What is not ported from `src/unitary.rs`

Seven of the fifty tests have no theorem here. Three of them *cannot*:

* `equiv_panics_on_measurement_left`, `equiv_panics_on_measurement_right`,
  `equiv_panics_on_reset` assert that Rust's `circuits_equiv` **panics** when handed a
  circuit containing `measure` or `reset` — it compares unitaries and those gates have none.
  The semantics here is a channel semantics, so those circuits have a meaning:
  `measure_plus_outcome`, `measure_kills_coherence` and `reset_one_to_zero` above compute it.

The other four are simply not done, and are the honest gaps in this file:

* `ccx_equiv_decomposition` — the fifteen-gate Toffoli decomposition, an `8 × 8` identity
  with `T` gates on three wires.
* `ghz_3qubit_two_ways`, `swap_3qubit_cycle`, `t_phase_folding_3qubit` — three-wire
  refutations. The machinery for them is here (`not_equivalent_of_entry`,
  `not_equivalent_of_perm`, `embed1_mul_apply`); each is a page of arithmetic that was not
  written. Their *kind* is covered: `cnot_no_commute_overlapping` and
  `four_qubit_cnot_ladder_not_reversible` are the same shape at three and four wires.
-/

end
end TzapLean
