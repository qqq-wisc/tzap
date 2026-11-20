import PauliFold.Paper.Circuits

/-!
# Paper: three claims of Sections 2 and 4 about non-Clifford behaviour

* Section 2.2: the `T` gate takes `X_q` to `(X_q + Y_q)/√2`, a sum of Pauli strings that is
  not a Pauli string (`T_conj_X`, `overview_T_not_pauli`, and `soundness_example_X2_T` for the circuit of
  Example 4.3).
* Example 4.7: `γ(α(T_q))` contains operators that are not multiples of `T_q`
  (`example_4_7_not_only_T`).
* Section 2.2: `T₁; H₁; T₁; H₁` and `H₁; T₁; H₁; T₁` have the same abstraction up to the order of
  the axes, but are not equivalent up to a global phase (`overview_not_equiv`).

The paper's `T` gate is `diag(1, e^{iπ/4})`; here it is `Circuit.T q = R_Z^q(π/4)`, which differs
from it by the global phase `e^{iπ/8}` (`gateUnitary_t`). Neither conjugation nor equivalence up
to a global phase sees the difference. Qubits are 0-indexed.
-/

namespace PauliFold.Paper

open TzapLean Matrix Circuit

noncomputable section

variable {n : ℕ}

/-! ## Products of `X_q`, `Y_q`, `Z_q` -/

theorem Xq_sq (q : Fin n) : (Xq q).toMatrix * (Xq q).toMatrix = 1 := by
  rw [← PP.toMatrix_mul, PP.mul_single_eq, ← PP.toMatrix_one]
  congr 1
  apply PP.ext'
  · rfl
  · intro r; simp [mulLetter, PP.one]

theorem XZ_mat (q : Fin n) :
    (Xq q).toMatrix * (Zq q).toMatrix = (-Complex.I) • (Yq q).toMatrix := by
  rw [← PP.toMatrix_mul, XZ_eq_Y, PP.toMatrix_scale, ipow_three]

theorem ZX_mat (q : Fin n) :
    (Zq q).toMatrix * (Xq q).toMatrix = Complex.I • (Yq q).toMatrix := by
  rw [← PP.toMatrix_mul, ZX_eq_Y, PP.toMatrix_scale, ipow_one]

theorem XZ_anticomm (q : Fin n) :
    (Xq q).toMatrix * (Zq q).toMatrix = -((Zq q).toMatrix * (Xq q).toMatrix) := by
  rw [XZ_mat, ZX_mat, neg_smul]

/-- `X_q R_Z^q(θ) = R_Z^q(-θ) X_q`. -/
theorem Xq_rzMat (q : Fin n) (θ : ℝ) :
    (Xq q).toMatrix * rzMat q θ = rzMat q (-θ) * (Xq q).toMatrix := by
  rw [rzMat_eq, rzMat_eq, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_smul, Matrix.mul_smul,
    Matrix.smul_mul, Matrix.smul_mul, Matrix.mul_one, Matrix.one_mul, XZ_anticomm, neg_div,
    Real.cos_neg, Real.sin_neg, smul_neg, sub_neg_eq_add, sub_eq_add_neg, ← neg_smul]
  push_cast
  congr 2
  ring

/-- Conjugating `X_q` by `R_Z^q(θ)` gives `R_Z^q(2θ) X_q = cos θ X_q + sin θ Y_q`
(`overview_rz_conj_X` on any number of qubits). -/
theorem rzMat_conj_X (q : Fin n) (θ : ℝ) :
    rzMat q θ * (Xq q).toMatrix * (rzMat q θ)ᴴ = rzMat q (2 * θ) * (Xq q).toMatrix := by
  rw [rzMat_conjTranspose, Matrix.mul_assoc, Xq_rzMat, neg_neg, ← Matrix.mul_assoc, rzMat_mul,
    two_mul]

theorem rzMat_mul_X (q : Fin n) (θ : ℝ) :
    rzMat q (2 * θ) * (Xq q).toMatrix =
      (Real.cos θ : ℂ) • (Xq q).toMatrix + (Real.sin θ : ℂ) • (Yq q).toMatrix := by
  rw [rzMat_eq, Matrix.sub_mul, Matrix.smul_mul, Matrix.smul_mul, Matrix.one_mul, ZX_mat,
    smul_smul, mul_div_cancel_left₀ θ two_ne_zero, sub_eq_add_neg, ← neg_smul]
  congr 2
  ring_nf; rw [Complex.I_sq]; ring

/-- **Paper, Section 2.2**: the `T` gate on qubit `q` takes `X_q` to
`(X_q + Y_q)/√2`, a sum of Pauli strings. -/
theorem T_conj_X (q : Fin n) :
    sem (T q) * (Xq q).toMatrix * (sem (T q))ᴴ =
      (1 / Real.sqrt 2 : ℂ) • ((Xq q).toMatrix + (Yq q).toMatrix) := by
  rw [sem_rz, rzMat_conj_X, show 2 * (Real.pi / 4) = Real.pi / 2 by ring,
    show Real.pi / 2 = 2 * (Real.pi / 4) by ring, rzMat_mul_X, Real.cos_pi_div_four,
    Real.sin_pi_div_four, smul_add]
  have h2 : (Real.sqrt 2 : ℂ) ≠ 0 := by exact_mod_cast (Real.sqrt_pos.2 two_pos).ne'
  have : ((Real.sqrt 2 / 2 : ℝ) : ℂ) = 1 / (Real.sqrt 2 : ℂ) := by
    push_cast
    field_simp
    rw [← Complex.ofReal_pow, Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)]
    norm_num
  rw [this]

theorem one_ne_zero_density : (1 : Density n) ≠ 0 := by
  intro h
  have := congrFun (congrFun h (fun _ => false)) (fun _ => false)
  simp at this

theorem Zq_ne_zero (q : Fin n) : (Zq q).toMatrix ≠ 0 := by
  intro h
  have := congrFun (congrFun h (fun _ => false)) (fun _ => false)
  simp [Zq_apply] at this

/-- `R_Z^q(θ)` is a Pauli string only if `cos(θ/2) = 0` or `sin(θ/2) = 0`. -/
theorem rzMat_not_mem_Paulis (q : Fin n) (θ : ℝ) (hc : Real.cos (θ / 2) ≠ 0)
    (hs : Real.sin (θ / 2) ≠ 0) : rzMat q θ ∉ Paulis n := by
  rintro ⟨P, hP⟩
  have hXR := Xq_rzMat q θ
  have cancel : ∀ A B : Density n, A * (Xq q).toMatrix = B * (Xq q).toMatrix → A = B := by
    intro A B hAB
    have := congrArg (· * (Xq q).toMatrix) hAB
    simpa only [Matrix.mul_assoc, Xq_sq, Matrix.mul_one] using this
  rcases PP.comm_or_anticomm P (Xq q) with h | h
  · -- `R X = X R`, so `R(θ) = R(-θ)`, so `2 i sin(θ/2) Z = 0`.
    rw [hP] at h
    have e := cancel _ _ (h.trans hXR)
    rw [rzMat_eq, rzMat_eq, neg_div, Real.cos_neg, Real.sin_neg, Complex.ofReal_neg] at e
    have e2 : (2 * Complex.I * Real.sin (θ / 2)) • (Zq q).toMatrix = 0 := by
      have := sub_eq_zero.2 e
      rw [← neg_eq_zero, ← this]; module
    rcases smul_eq_zero.1 e2 with h0 | h0
    · rcases mul_eq_zero.1 h0 with h1 | h1
      · simp [Complex.I_ne_zero] at h1
      · exact hs (Complex.ofReal_eq_zero.1 h1)
    · exact Zq_ne_zero q h0
  · -- `R X = -X R`, so `R(θ) = -R(-θ)`, so `2 cos(θ/2) I = 0`.
    rw [hP, hXR, ← Matrix.neg_mul] at h
    have e := cancel _ _ h
    rw [rzMat_eq, rzMat_eq, neg_div, Real.cos_neg, Real.sin_neg, Complex.ofReal_neg] at e
    have e2 : (2 * (Real.cos (θ / 2) : ℂ)) • (1 : Density n) = 0 := by
      have := sub_eq_zero.2 e
      rw [← this]; module
    rcases smul_eq_zero.1 e2 with h0 | h0
    · rcases mul_eq_zero.1 h0 with h1 | h1
      · norm_num at h1
      · exact hc (Complex.ofReal_eq_zero.1 h1)
    · exact one_ne_zero_density h0

/-- **Paper, Section 2.2**: the `T` gate takes `X_q` to a matrix that is not a
Pauli string, `T X_q T† ∉ 𝒫_n`. -/
theorem overview_T_not_pauli (q : Fin n) :
    sem (T q) * (Xq q).toMatrix * (sem (T q))ᴴ ∉ Paulis n := by
  rw [sem_rz, rzMat_conj_X]
  rintro ⟨P, hP⟩
  apply rzMat_not_mem_Paulis q (2 * (Real.pi / 4))
  · rw [show 2 * (Real.pi / 4) / 2 = Real.pi / 4 by ring, Real.cos_pi_div_four]; positivity
  · rw [show 2 * (Real.pi / 4) / 2 = Real.pi / 4 by ring, Real.sin_pi_div_four]; positivity
  refine ⟨P.mul (Xq q), ?_⟩
  rw [PP.toMatrix_mul, hP, Matrix.mul_assoc, Xq_sq, Matrix.mul_one]

/-- **Paper, Section 4.4**, for the relation `σ` of Example 4.3 (`H₁; CX_{1,2}; T₂`, Lean qubits
`0, 1`): `X₂ ∉ 𝓕(σ)`, and the `T₂` gate takes `X₂` to `(X₂ + Y₂)/√2`, which is not a Pauli
string. -/
theorem soundness_example_X2_T :
    Xq 1 ∉ Examples.exSigma.F ∧
    sem (T (1 : Fin 2)) * (Xq 1).toMatrix * (sem (T 1))ᴴ =
      (1 / Real.sqrt 2 : ℂ) • ((Xq (1 : Fin 2)).toMatrix + (Yq 1).toMatrix) ∧
    sem (T (1 : Fin 2)) * (Xq 1).toMatrix * (sem (T 1))ᴴ ∉ Paulis 2 :=
  ⟨Examples.soundness_example_X2.2.2, T_conj_X 1, overview_T_not_pauli 1⟩

/-! ## Example 4.7 -/

/-- **Paper, Example 4.7**: `T_q` is in `γ(α(T_q))`, "but not the only one": the identity, the
Z-rotation `R_Z^q(0)`, is in `γ(α(T_q))` and is not a scalar multiple of `T_q`. (Every other
Z-rotation of `q` is in `γ(α(T_q))` too: `Examples.example_4_7_rz`.) -/
theorem example_4_7_not_only_T (q : Fin n) :
    sem (T q) ∈ (alpha (T q)).gamma ∧ rzMat q 0 ∈ (alpha (T q)).gamma ∧
      ¬ ∃ c : ℂ, rzMat q 0 = c • sem (T q) := by
  refine ⟨soundness _, Examples.example_4_7_rz q 0, ?_⟩
  rintro ⟨c, hc⟩
  rw [rzMat_zero, sem_rz, rzMat_eq] at hc
  set o : Basis n := fun _ => false
  set o' : Basis n := Function.update o q true
  have e₀ := congrFun (congrFun hc o) o
  have e₁ := congrFun (congrFun hc o') o'
  simp only [Matrix.smul_apply, Matrix.sub_apply, Matrix.one_apply_eq, Zq_apply, if_pos,
    smul_eq_mul, o, o', Function.update_self] at e₀ e₁
  simp only [Bool.false_eq_true, if_false, mul_one] at e₀ e₁
  have hs : Real.sin (Real.pi / 4 / 2) ≠ 0 :=
    (Real.sin_pos_of_pos_of_lt_pi (by positivity) (by linarith [Real.pi_pos])).ne'
  have : c * (2 * Complex.I * Real.sin (Real.pi / 4 / 2)) = 0 := by
    linear_combination e₀ - e₁
  rcases mul_eq_zero.1 this with h | h
  · rw [h] at e₀; simp at e₀
  · rcases mul_eq_zero.1 h with h1 | h1
    · simp [Complex.I_ne_zero] at h1
    · exact hs (Complex.ofReal_eq_zero.1 h1)

/-! ## Section 2.2: same abstraction, different unitaries -/

/-- `H Z H = X` and `H H = 1` for the Hadamard gate. -/
theorem gH_conj_Z (q : Fin n) :
    (gH q).mat * (Zq q).toMatrix * (gH q).mat = (Xq q).toMatrix := by
  have h := (gH q).zs_spec q
  have hHt : (gH q).matᴴ = (gH q).mat := (gateUnitary_inv (CTGate.h q)).symm
  rw [hHt] at h
  rw [h]; simp [gH, ofCT]

theorem HTH_eq (q : Fin n) (a b : ℂ) :
    (gH q).mat * (a • (1 : Density n) + b • (Zq q).toMatrix) * (gH q).mat =
      a • (1 : Density n) + b • (Xq q).toMatrix := by
  rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.mul_smul, Matrix.smul_mul,
    Matrix.smul_mul, Matrix.mul_one, gH_mul_self, gH_conj_Z]

theorem Xq_diag (q : Fin n) (o : Basis n) : (Xq q).toMatrix o o = 0 := by
  rw [Xq_apply, if_neg]
  intro h
  have := congrFun h q
  simp [flipQ] at this

theorem Zq_mul_diag (q : Fin n) (M : Density n) (o : Basis n) :
    ((Zq q).toMatrix * M) o o = (if o q then -1 else 1) * M o o := by
  rw [Matrix.mul_apply, Finset.sum_eq_single o]
  · rw [Zq_apply, if_pos rfl]
  · intro k _ hk; rw [Zq_apply, if_neg (Ne.symm hk), zero_mul]
  · simp

theorem mul_Zq_diag (q : Fin n) (M : Density n) (o : Basis n) :
    (M * (Zq q).toMatrix) o o = M o o * (if o q then -1 else 1) := by
  rw [Matrix.mul_apply, Finset.sum_eq_single o]
  · rw [Zq_apply, if_pos rfl]
  · intro k _ hk; rw [Zq_apply, if_neg hk, mul_zero]
  · simp

/-- With `M = a I + b Z_q` and `N = a I + b X_q`, `N M ≡ M N` (up to a global phase) forces
`a(a + b) = 0` or `b = 0`. -/
theorem combo_not_equiv (q : Fin n) {a b : ℂ} (hab : a * (a + b) ≠ 0) (hb : b ≠ 0) :
    ¬ EquivP ((a • (1 : Density n) + b • (Xq q).toMatrix) * (a • 1 + b • (Zq q).toMatrix))
      ((a • (1 : Density n) + b • (Zq q).toMatrix) * (a • 1 + b • (Xq q).toMatrix)) := by
  rintro ⟨φ, h⟩
  set e := Complex.exp (φ * Complex.I)
  have expand : ∀ A B : Density n, (a • (1 : Density n) + b • A) * (a • 1 + b • B) =
      (a * a) • (1 : Density n) + (a * b) • B + (a * b) • A + (b * b) • (A * B) := by
    intro A B
    simp only [Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul, Matrix.mul_smul, Matrix.one_mul,
      Matrix.mul_one]
    module
  rw [expand, expand] at h
  -- The `(o, o)` entry, `o = |0⋯0⟩`: `a² + ab = e (a² + ab)`.
  set o : Basis n := fun _ => false
  have hd := congrFun (congrFun h o) o
  simp only [Matrix.add_apply, Matrix.smul_apply, Matrix.one_apply_eq, smul_eq_mul, Xq_diag,
    Zq_apply, if_pos, Zq_mul_diag, mul_Zq_diag, o, Bool.false_eq_true, if_false] at hd
  have he : e = 1 := by
    have : (a * (a + b)) * (e - 1) = 0 := by linear_combination -hd
    rcases mul_eq_zero.1 this with h0 | h0
    · exact absurd h0 hab
    · exact sub_eq_zero.1 h0
  rw [he, one_smul] at h
  -- Then `b² X Z = b² Z X = -b² X Z`, so `2 b² X Z = 0`, and `X Z Z X = 1`.
  have h2 : (2 * (b * b)) • ((Xq q).toMatrix * (Zq q).toMatrix) = 0 := by
    have := sub_eq_zero.2 h
    rw [XZ_anticomm] at this ⊢
    rw [← this]; module
  have h3 : (2 * (b * b)) • (1 : Density n) = 0 := by
    have := congrArg (· * ((Zq q).toMatrix * (Xq q).toMatrix)) h2
    simp only [Matrix.smul_mul, Matrix.zero_mul] at this
    rwa [Matrix.mul_assoc, ← Matrix.mul_assoc (Zq q).toMatrix, Zq_sq, Matrix.one_mul, Xq_sq]
      at this
  rcases smul_eq_zero.1 h3 with h0 | h0
  · simp [hb] at h0
  · exact one_ne_zero_density h0

/-- **Paper, Section 2.2**: `T₁; H₁; T₁; H₁` and `H₁; T₁; H₁; T₁` (one qubit) have the same
abstraction, the tableau of the identity with the axes `Z₁, X₁` (in a different order, as the
axes are a list), but they are not equivalent up to a global phase. -/
theorem overview_not_equiv :
    alpha (T 0 ⨾ gate (gH (0 : Fin 1)) ⨾ T 0 ⨾ gate (gH 0)) = .rel ⟨Cons.id 1, [Zq 0, Xq 0]⟩ ∧
    alpha (gate (gH (0 : Fin 1)) ⨾ T 0 ⨾ gate (gH 0) ⨾ T 0) = .rel ⟨Cons.id 1, [Xq 0, Zq 0]⟩ ∧
    ¬ (T 0 ⨾ gate (gH (0 : Fin 1)) ⨾ T 0 ⨾ gate (gH 0)).Equiv
      (gate (gH (0 : Fin 1)) ⨾ T 0 ⨾ gate (gH 0) ⨾ T 0) := by
  refine ⟨Examples.overview_same_abstraction.1, Examples.overview_same_abstraction.2, ?_⟩
  set a : ℂ := (Real.cos (Real.pi / 4 / 2) : ℂ)
  set b : ℂ := -(Complex.I * Real.sin (Real.pi / 4 / 2))
  have hT : rzMat (0 : Fin 1) (Real.pi / 4) = a • (1 : Density 1) + b • (Zq 0).toMatrix := by
    rw [rzMat_eq, sub_eq_add_neg, ← neg_smul]
  have hc : 0 < Real.cos (Real.pi / 4 / 2) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith [Real.pi_pos]⟩
  have hs : 0 < Real.sin (Real.pi / 4 / 2) :=
    Real.sin_pos_of_pos_of_lt_pi (by positivity) (by linarith [Real.pi_pos])
  have hab : a * (a + b) ≠ 0 := by
    apply mul_ne_zero
    · exact Complex.ofReal_ne_zero.2 hc.ne'
    · intro h
      have := congrArg Complex.re h
      simp only [a, b, Complex.add_re, Complex.neg_re, Complex.mul_re, Complex.I_re,
        Complex.I_im, Complex.ofReal_re, Complex.ofReal_im, Complex.zero_re] at this
      linarith
  have hb : b ≠ 0 := by
    intro h
    have := congrArg Complex.im h
    simp only [b, Complex.neg_im, Complex.mul_im, Complex.I_re, Complex.I_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.zero_im] at this
    linarith
  intro hE
  apply combo_not_equiv (0 : Fin 1) hab hb
  unfold Circuit.Equiv at hE
  simp only [sem_seq, sem_gate, sem_rz, hT] at hE
  rw [← HTH_eq (0 : Fin 1) a b]
  simpa only [Matrix.mul_assoc] using hE

end

end PauliFold.Paper
