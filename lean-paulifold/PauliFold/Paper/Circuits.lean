import PauliFold.Paper.Relation

/-!
# Paper, Section 4: abstracting circuits, soundness, and the examples

* `Circuit.alpha`: the abstraction of Definitions 4.8, 4.10 and 4.13.
* `Circuit.alpha_spec`: `α(C)` is never `⊤`, all its constraints are known, and the Clifford
  part `Circuit.cliff C` of the circuit is a model.
* Lemma 4.14 (Figure 2): `lemma_4_14_*`.
* Theorem 4.15 (Soundness): `Circuit.soundness`.
* "The tableau loses nothing": for a circuit of Clifford gates only, `γ(α(C))` is exactly the
  multiples of `⟦C⟧` (`Circuit.gamma_alpha_of_isClifford`), with the instances `H` and `H; H`
  of Sections 2.1 and 7.1.
* Examples 4.3, 4.7, 4.9 and 4.12, the soundness example of Section 4.4, and the Section 2 claims
  they repeat.

Qubits are 0-indexed: the paper's qubits 1 and 2 are the Lean qubits `0` and `1`.
-/

namespace PauliFold.Paper

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## The abstraction (Definitions 4.8, 4.10, 4.13) -/

/-- Paper, Definition 4.8: the constraints `Γ_G` of a Clifford gate,
`Gᴴ X_q G ↦ X_q` and `Gᴴ Z_q G ↦ Z_q` (the gate's tables). -/
def CliffordGate.cons (G : CliffordGate n) : Cons n := Cons.ofTables G.xs G.zs

/-- The constraints of the identity gate are `Γ_I`. -/
theorem CliffordGate.one_cons : (CliffordGate.one n).cons = Cons.id n := rfl

/-- `Γ_G` is exactly Definition 4.8: each constraint is `Gᴴ X_q G`, respectively `Gᴴ Z_q G`. -/
theorem CliffordGate.cons_spec (G : CliffordGate n) (q : Fin n) :
    (∃ P, G.cons.x q = some P ∧ P.toMatrix = G.matᴴ * (Xq q).toMatrix * G.mat) ∧
    (∃ P, G.cons.z q = some P ∧ P.toMatrix = G.matᴴ * (Zq q).toMatrix * G.mat) :=
  ⟨⟨_, rfl, (G.xs_spec q).symm⟩, ⟨_, rfl, (G.zs_spec q).symm⟩⟩

namespace Circuit

/-- Paper, Definitions 4.8, 4.10 and 4.13 (Abstraction): `α(G) = (Γ_G, ∅)`,
`α(R_Z^q(θ)) = (Γ_I, {Z_q})`, `α(C₁; C₂) = α(C₂) ∘ α(C₁)`. -/
def alpha : Circuit n → Dom n
  | gate G => .rel ⟨G.cons, []⟩
  | rz q _ => .rel ⟨Cons.id n, [Zq q]⟩
  | seq C₁ C₂ => alpha C₂ ⊚ alpha C₁

@[simp] theorem alpha_gate (G : CliffordGate n) : alpha (gate G) = .rel ⟨G.cons, []⟩ := rfl
@[simp] theorem alpha_rz (q : Fin n) (θ : ℝ) : alpha (rz q θ) = .rel ⟨Cons.id n, [Zq q]⟩ := rfl
@[simp] theorem alpha_seq (C₁ C₂ : Circuit n) : alpha (C₁ ⨾ C₂) = alpha C₂ ⊚ alpha C₁ := rfl

/-- The Clifford part of a circuit: its gates, with the rotations removed. -/
def cliff : Circuit n → Density n
  | gate G => G.mat
  | rz _ _ => 1
  | seq C₁ C₂ => cliff C₂ * cliff C₁

/-- **`α(C)` is an abstract relation of Section 4**: it is not `⊤`, all its constraints are
known, and the Clifford part of `C` is a model (Section 4.3). -/
theorem alpha_spec : ∀ C : Circuit n,
    ∃ σ, alpha C = .rel σ ∧ σ.cons.AllKnown ∧ Models (cliff C) σ.cons
  | gate G => ⟨_, rfl, Cons.ofTables_allKnown _ _, G.models⟩
  | rz q _ => ⟨_, rfl, Cons.ofTables_allKnown _ _, models_id⟩
  | seq C₁ C₂ => by
    obtain ⟨σ₁, h₁, k₁, m₁⟩ := alpha_spec C₁
    obtain ⟨σ₂, h₂, k₂, m₂⟩ := alpha_spec C₂
    refine ⟨⟨σ₂.cons.comp σ₁.cons, σ₁.axes ++ σ₂.axes.filterMap σ₁.cons.wp⟩, ?_,
      Cons.AllKnown.comp k₁ k₂, m₁.comp m₂⟩
    rw [alpha_seq, h₁, h₂, Dom.rel_comp_rel, ARel.comp_of_allKnown k₁]

/-- `α(C)` as an abstract relation (it is never `⊤`). -/
def alphaRel (C : Circuit n) : ARel n :=
  match alpha C with
  | .rel σ => σ
  | .top => ⟨Cons.id n, []⟩

theorem alpha_eq_rel (C : Circuit n) : alpha C = .rel (alphaRel C) := by
  obtain ⟨σ, h, -⟩ := alpha_spec C
  rw [alphaRel, h]

theorem alphaRel_allKnown (C : Circuit n) : (alphaRel C).cons.AllKnown := by
  obtain ⟨σ, h, hk, -⟩ := alpha_spec C
  rw [alphaRel, h]; exact hk

theorem models_alphaRel (C : Circuit n) : Models (cliff C) (alphaRel C).cons := by
  obtain ⟨σ, h, -, hm⟩ := alpha_spec C
  rw [alphaRel, h]; exact hm

theorem alpha_sat (C : Circuit n) : (alpha C).Sat := by
  rw [alpha_eq_rel]; exact ⟨_, models_alphaRel C⟩

@[simp] theorem alphaRel_gate (G : CliffordGate n) : alphaRel (gate G) = ⟨G.cons, []⟩ := rfl
@[simp] theorem alphaRel_rz (q : Fin n) (θ : ℝ) :
    alphaRel (rz q θ) = ⟨Cons.id n, [Zq q]⟩ := rfl

/-- `α(C₁; C₂)` in terms of the parts. -/
theorem alphaRel_seq (C₁ C₂ : Circuit n) :
    alphaRel (C₁ ⨾ C₂) = ⟨(alphaRel C₂).cons.comp (alphaRel C₁).cons,
      (alphaRel C₁).axes ++ (alphaRel C₂).axes.filterMap (alphaRel C₁).cons.wp⟩ := by
  have h := alpha_eq_rel (C₁ ⨾ C₂)
  rw [alpha_seq, alpha_eq_rel C₁, alpha_eq_rel C₂, Dom.rel_comp_rel,
    ARel.comp_of_allKnown (alphaRel_allKnown C₁)] at h
  exact (Dom.rel.inj h).symm

/-- Lemma 4.5 for `α(C)`: every weakest precondition is known, and is the conjugate by the
Clifford part. -/
theorem wp_alpha (C : Circuit n) (P : PP n) :
    ∃ B, (alphaRel C).cons.wp P = some B ∧ B.toMatrix = (cliff C)ᴴ * P.toMatrix * cliff C :=
  Cons.wp_conj (alphaRel_allKnown C) (models_alphaRel C) P

/-- The weakest precondition through `α(C₁; C₂)` is the one through `α(C₂)`, then `α(C₁)`. -/
theorem wp_alphaRel_seq (C₁ C₂ : Circuit n) (P : PP n) :
    (alphaRel (C₁ ⨾ C₂)).cons.wp P =
      ((alphaRel C₂).cons.wp P).bind (alphaRel C₁).cons.wp := by
  rw [alphaRel_seq]
  exact Cons.wp_comp_eq (models_alphaRel C₁) (alphaRel_allKnown C₁) P

end Circuit

/-! ## Lemma 4.14: the abstractions of the standard gates (Figure 2) -/

open Circuit

/-- Paper, Lemma 4.14 for `H_q`: `Z_q ↦ X_q`, `X_q ↦ Z_q`, the other constraints those of the
identity, no axes; each constraint is `Hᴴ X_r H`, respectively `Hᴴ Z_r H`. The matrix is tzap's
`H`. -/
theorem lemma_4_14_H (q : Fin n) :
    (gH q).mat = gateUnitary n (.h q.val) ∧
    alpha (gate (gH q)) = .rel ⟨⟨fun r => some (if r = q then Zq q else Xq r),
      fun r => some (if r = q then Xq q else Zq r)⟩, []⟩ ∧
    ∀ r, (gH q).matᴴ * (Xq r).toMatrix * (gH q).mat = (if r = q then Zq q else Xq r).toMatrix ∧
      (gH q).matᴴ * (Zq r).toMatrix * (gH q).mat = (if r = q then Xq q else Zq r).toMatrix :=
  ⟨rfl, rfl, fun r => ⟨(gH q).xs_spec r, (gH q).zs_spec r⟩⟩

/-- Paper, Lemma 4.14 for `S_q`: `-Y_q ↦ X_q`. -/
theorem lemma_4_14_S (q : Fin n) :
    (gS q).mat = gateUnitary n (.s q.val) ∧
    alpha (gate (gS q)) = .rel ⟨⟨fun r => some (if r = q then (Yq q).neg else Xq r),
      fun r => some (Zq r)⟩, []⟩ ∧
    ∀ r, (gS q).matᴴ * (Xq r).toMatrix * (gS q).mat =
        (if r = q then (Yq q).neg else Xq r).toMatrix ∧
      (gS q).matᴴ * (Zq r).toMatrix * (gS q).mat = (Zq r).toMatrix :=
  ⟨rfl, rfl, fun r => ⟨(gS q).xs_spec r, (gS q).zs_spec r⟩⟩

/-- Paper, Lemma 4.14 for `X_q`: `-Z_q ↦ Z_q`. -/
theorem lemma_4_14_X (q : Fin n) :
    (gX q).mat = gateUnitary n (.x q.val) ∧
    alpha (gate (gX q)) = .rel ⟨⟨fun r => some (Xq r),
      fun r => some (if r = q then (Zq q).neg else Zq r)⟩, []⟩ ∧
    ∀ r, (gX q).matᴴ * (Xq r).toMatrix * (gX q).mat = (Xq r).toMatrix ∧
      (gX q).matᴴ * (Zq r).toMatrix * (gX q).mat =
        (if r = q then (Zq q).neg else Zq r).toMatrix :=
  ⟨rfl, rfl, fun r => ⟨(gX q).xs_spec r, (gX q).zs_spec r⟩⟩

/-- Paper, Lemma 4.14 for `Z_q`: `-X_q ↦ X_q`. -/
theorem lemma_4_14_Z (q : Fin n) :
    (gZ q).mat = gateUnitary n (.z q.val) ∧
    alpha (gate (gZ q)) = .rel ⟨⟨fun r => some (if r = q then (Xq q).neg else Xq r),
      fun r => some (Zq r)⟩, []⟩ ∧
    ∀ r, (gZ q).matᴴ * (Xq r).toMatrix * (gZ q).mat =
        (if r = q then (Xq q).neg else Xq r).toMatrix ∧
      (gZ q).matᴴ * (Zq r).toMatrix * (gZ q).mat = (Zq r).toMatrix :=
  ⟨rfl, rfl, fun r => ⟨(gZ q).xs_spec r, (gZ q).zs_spec r⟩⟩

/-- Paper, Lemma 4.14 for `CX_{c,t}`: `X_c X_t ↦ X_c`, `Z_c Z_t ↦ Z_t`. -/
theorem lemma_4_14_CX (c t : Fin n) (hne : c ≠ t) :
    (gCX c t hne).mat = gateUnitary n (.cnot c.val t.val) ∧
    alpha (gate (gCX c t hne)) = .rel ⟨⟨fun r => some (if r = c then (Xq c).mul (Xq t) else Xq r),
      fun r => some (if r = t then (Zq c).mul (Zq t) else Zq r)⟩, []⟩ ∧
    ∀ r, (gCX c t hne).matᴴ * (Xq r).toMatrix * (gCX c t hne).mat =
        (if r = c then (Xq c).mul (Xq t) else Xq r).toMatrix ∧
      (gCX c t hne).matᴴ * (Zq r).toMatrix * (gCX c t hne).mat =
        (if r = t then (Zq c).mul (Zq t) else Zq r).toMatrix :=
  ⟨rfl, rfl, fun r => ⟨(gCX c t hne).xs_spec r, (gCX c t hne).zs_spec r⟩⟩

/-- Lemma 4.14 for `S†_q` (not in Figure 2): `Y_q ↦ X_q`. -/
theorem lemma_4_14_Sdg (q : Fin n) :
    (gSdg q).mat = gateUnitary n (.sdg q.val) ∧
    alpha (gate (gSdg q)) = .rel ⟨⟨fun r => some (if r = q then Yq q else Xq r),
      fun r => some (Zq r)⟩, []⟩ ∧
    ∀ r, (gSdg q).matᴴ * (Xq r).toMatrix * (gSdg q).mat = (if r = q then Yq q else Xq r).toMatrix ∧
      (gSdg q).matᴴ * (Zq r).toMatrix * (gSdg q).mat = (Zq r).toMatrix :=
  ⟨rfl, rfl, fun r => ⟨(gSdg q).xs_spec r, (gSdg q).zs_spec r⟩⟩

/-- Lemma 4.14 for `CZ_{c,t}` (not in Figure 2): `X_c Z_t ↦ X_c`, `Z_c X_t ↦ X_t`. -/
theorem lemma_4_14_CZ (c t : Fin n) (hne : c ≠ t) :
    (gCZ c t hne).mat = gateUnitary n (.cz c.val t.val) ∧
    alpha (gate (gCZ c t hne)) = .rel ⟨⟨fun r => some
        (if r = c then (Xq c).mul (Zq t) else if r = t then (Zq c).mul (Xq t) else Xq r),
      fun r => some (Zq r)⟩, []⟩ ∧
    ∀ r, (gCZ c t hne).matᴴ * (Xq r).toMatrix * (gCZ c t hne).mat =
        (if r = c then (Xq c).mul (Zq t) else if r = t then (Zq c).mul (Xq t) else Xq r).toMatrix ∧
      (gCZ c t hne).matᴴ * (Zq r).toMatrix * (gCZ c t hne).mat = (Zq r).toMatrix :=
  ⟨rfl, rfl, fun r => ⟨(gCZ c t hne).xs_spec r, (gCZ c t hne).zs_spec r⟩⟩

/-- Paper, Lemma 4.14 for `R_Z^q(θ)`: no constraint changes, one axis `Z_q`. -/
theorem lemma_4_14_rz (q : Fin n) (θ : ℝ) : alpha (rz q θ) = .rel ⟨Cons.id n, [Zq q]⟩ := rfl

/-! ## Soundness (Theorem 4.15) -/

namespace Circuit

/-- A Clifford gate takes `wp(P)` to `P` exactly: `G ∈ γ(α(G))`. -/
theorem gate_mem_gamma (G : CliffordGate n) : G.mat ∈ (alpha (gate G)).gamma := by
  intro P B hB _
  rw [Cons.wp_toMatrix G.models hB]
  simp only [← Matrix.mul_assoc]
  rw [G.mul_conj, Matrix.one_mul]

/-- Every Z-rotation of qubit `q`, whatever its angle, is in `γ(α(R_Z^q(θ')))`: it commutes
with every Pauli string that commutes with `Z_q`. -/
theorem rz_mem_gamma (q : Fin n) (θ θ' : ℝ) : rzMat q θ ∈ (alpha (rz q θ')).gamma := by
  intro P B hB hc
  rw [Cons.wp_id] at hB
  cases hB
  exact rzMat_comm (hc _ (List.mem_singleton_self _)) θ

/-- **Paper, Theorem 4.15 (Soundness)**: `⟦C⟧ ∈ γ(α(C))` for every circuit `C`. -/
theorem soundness : ∀ C : Circuit n, sem C ∈ (alpha C).gamma
  | gate G => gate_mem_gamma G
  | rz q θ => rz_mem_gamma q θ θ
  | seq C₁ C₂ => gamma_comp (alpha_sat C₁) (soundness C₁) (soundness C₂)

/-- Soundness, as a statement about tracked strings: if `P ∈ 𝓕(α(C))`, then `⟦C⟧` takes
`wp(P)` to `P`. -/
theorem takes_of_tracked (C : Circuit n) {P B : PP n} (hB : (alphaRel C).cons.wp P = some B)
    (hc : ∀ A ∈ (alphaRel C).axes, PP.Comm B A) : sem C * B.toMatrix = P.toMatrix * sem C := by
  have h := soundness C
  rw [alpha_eq_rel] at h
  exact h P B hB hc

/-! ## The tableau loses nothing -/

theorem cliff_eq_sem : ∀ {C : Circuit n}, C.IsClifford → cliff C = sem C
  | gate _, _ => rfl
  | rz _ _, h => h.elim
  | seq _ _, h => by rw [cliff, sem, cliff_eq_sem h.1, cliff_eq_sem h.2]

theorem axes_nil_of_isClifford : ∀ {C : Circuit n}, C.IsClifford → (alphaRel C).axes = []
  | gate _, _ => rfl
  | rz _ _, h => h.elim
  | seq C₁ C₂, h => by
    rw [alphaRel_seq, axes_nil_of_isClifford h.1, axes_nil_of_isClifford h.2]; rfl

/-- **The tableau loses nothing.** For a circuit made only of Clifford gates, `γ(α(C))` is
exactly the set of multiples of `⟦C⟧`. -/
theorem gamma_alpha_of_isClifford {C : Circuit n} (hC : C.IsClifford) :
    (alpha C).gamma = {K | ∃ c : ℂ, K = c • sem C} := by
  ext K
  constructor
  · intro hK
    rw [alpha_eq_rel] at hK
    have hm := models_alphaRel C
    rw [cliff_eq_sem hC] at hm
    set V := sem C
    have hax := axes_nil_of_isClifford hC
    have key : ∀ P : PP n, (K * Vᴴ) * P.toMatrix = P.toMatrix * (K * Vᴴ) := by
      intro P
      obtain ⟨B, hB, hBV⟩ := Cons.wp_conj (alphaRel_allKnown C) hm P
      have e := hK P B hB (by rw [hax]; simp)
      rw [hBV] at e
      calc K * Vᴴ * P.toMatrix = K * (Vᴴ * P.toMatrix * V) * Vᴴ := by
            simp only [Matrix.mul_assoc]; rw [hm.mul_conj, Matrix.mul_one]
        _ = P.toMatrix * (K * Vᴴ) := by rw [e, Matrix.mul_assoc]
    obtain ⟨c, hc⟩ := eq_smul_one_of_comm_all _ key
    refine ⟨c, ?_⟩
    calc K = K * Vᴴ * V := by rw [Matrix.mul_assoc, hm.conj_mul, Matrix.mul_one]
      _ = c • V := by rw [hc, Matrix.smul_mul, Matrix.one_mul]
  · rintro ⟨c, rfl⟩
    exact Dom.smul_mem_gamma c (soundness C)

theorem gH_mul_self (q : Fin n) : (gH q).mat * (gH q).mat = 1 := by
  have h : (gH q).matᴴ = (gH q).mat := (gateUnitary_inv (CTGate.h q)).symm
  calc (gH q).mat * (gH q).mat = (gH q).mat * (gH q).matᴴ := by rw [h]
    _ = 1 := (gH q).mul_conj

/-- Sections 2.1 and 7.1: `γ(α(H))` is exactly the multiples of `H` (here on any number of
qubits). -/
theorem gamma_alpha_H (q : Fin n) :
    (alpha (gate (gH q))).gamma = {K | ∃ c : ℂ, K = c • (gH q).mat} :=
  gamma_alpha_of_isClifford trivial

/-- Section 7.1: `γ(α(H; H))` is exactly the multiples of the identity. -/
theorem gamma_alpha_HH (q : Fin n) :
    (alpha (gate (gH q) ⨾ gate (gH q))).gamma = {K | ∃ c : ℂ, K = c • (1 : Density n)} := by
  rw [gamma_alpha_of_isClifford (C := gate (gH q) ⨾ gate (gH q)) ⟨trivial, trivial⟩, sem_seq, sem_gate, gH_mul_self]

/-- Section 7.1: on one qubit, the tableau of `H; H` is the tableau of the identity. -/
theorem alpha_HH : alpha (gate (gH (0 : Fin 1)) ⨾ gate (gH 0)) = .rel ⟨Cons.id 1, []⟩ := by
  decide

/-- Section 2.1: on one qubit, `γ(α(H))` consists of the `K` with `K Z = X K` and `K X = Z K`,
which are exactly the multiples of `H`. -/
theorem gamma_alpha_H_one (K : Density 1) :
    K ∈ (alpha (gate (gH (0 : Fin 1)))).gamma ↔
      K * (Zq 0).toMatrix = (Xq 0).toMatrix * K ∧ K * (Xq 0).toMatrix = (Zq 0).toMatrix * K := by
  constructor
  · intro hK
    exact ⟨hK (Xq 0) (Zq 0) (by decide) (by simp), hK (Zq 0) (Xq 0) (by decide) (by simp)⟩
  · rintro ⟨h1, h2⟩
    rw [gamma_alpha_H]
    -- `K H` commutes with `X` and `Z`.
    set H := (gH (0 : Fin 1)).mat
    have hX := ((lemma_4_14_H (0 : Fin 1)).2.2 0).1
    have hZ := ((lemma_4_14_H (0 : Fin 1)).2.2 0).2
    simp only [↓reduceIte] at hX hZ
    have hHH : H * H = 1 := gH_mul_self 0
    have hHt : Hᴴ = H := (gateUnitary_inv (CTGate.h (0 : Fin 1))).symm
    rw [hHt] at hX hZ
    -- `H X = Z H` and `H Z = X H`.
    have hXH : H * (Xq 0).toMatrix = (Zq 0).toMatrix * H := by
      rw [← hX, Matrix.mul_assoc _ H H, hHH, Matrix.mul_one]
    have hZH : H * (Zq 0).toMatrix = (Xq 0).toMatrix * H := by
      rw [← hZ, Matrix.mul_assoc _ H H, hHH, Matrix.mul_one]
    obtain ⟨c, hc⟩ := eq_smul_one_of_comm (K * H)
      (fun q => by
        rw [Fin.fin_one_eq_zero q, Matrix.mul_assoc, hXH, ← Matrix.mul_assoc, h1,
          Matrix.mul_assoc])
      (fun q => by
        rw [Fin.fin_one_eq_zero q, Matrix.mul_assoc, hZH, ← Matrix.mul_assoc, h2,
          Matrix.mul_assoc])
    exact ⟨c, by rw [← Matrix.mul_one K, ← hHH, ← Matrix.mul_assoc, hc, Matrix.smul_mul,
      Matrix.one_mul]⟩

end Circuit

/-! ## Examples (qubits `1, 2` of the paper are `0, 1`) -/

namespace Examples

open Circuit

/-- `X₁ Z₂` (Lean qubits `0, 1`). -/
abbrev X1Z2 : PP 2 := (Xq 0).mul (Zq 1)
/-- `Z₁ X₂`. -/
abbrev Z1X2 : PP 2 := (Zq 0).mul (Xq 1)

/-- The circuit `H₁; CX_{1,2}; T₂` of Section 2.2 and Example 4.3. -/
def exCirc : Circuit 2 := gate (gH 0) ⨾ gate (gCX 0 1 (by decide)) ⨾ T 1

/-- The tableau of `H₁; CX_{1,2}` (Section 2.1, Example 4.12):
`Z₁X₂ ↦ X₁, X₁ ↦ Z₁, X₂ ↦ X₂, X₁Z₂ ↦ Z₂`. -/
def exCons : Cons 2 := ⟨![some Z1X2, some (Xq 1)], ![some (Xq 0), some X1Z2]⟩

/-- The abstract relation of Example 4.3: `exCons` with the axis `X₁Z₂`. -/
def exSigma : ARel 2 := ⟨exCons, [X1Z2]⟩

/-- Paper, Example 4.3: `α(H₁; CX_{1,2}; T₂)` is `σ`. -/
theorem example_4_3 : alpha exCirc = .rel exSigma := by decide

/-- Paper, Example 4.3: the Clifford part `V = CX_{1,2} H₁` is a model of the constraints. -/
theorem example_4_3_model :
    Models ((gCX (0 : Fin 2) 1 (by decide)).mat * (gH 0).mat) exSigma.cons := by
  have h := models_alphaRel exCirc
  have e : alphaRel exCirc = exSigma := by
    have := alpha_eq_rel exCirc
    rw [example_4_3] at this
    exact (Dom.rel.inj this).symm
  rw [e] at h
  simpa [exCirc, cliff] using h

/-- Paper, Example 4.3: `Vᴴ X₁ V = H₁ X₁X₂ H₁ = Z₁X₂`. -/
theorem example_4_3_conj :
    ((gCX (0 : Fin 2) 1 (by decide)).mat * (gH 0).mat)ᴴ * (Xq 0).toMatrix *
        ((gCX (0 : Fin 2) 1 (by decide)).mat * (gH 0).mat) = Z1X2.toMatrix ∧
      (gH (0 : Fin 2)).matᴴ * ((Xq 0).mul (Xq 1)).toMatrix * (gH 0).mat = Z1X2.toMatrix := by
  refine ⟨(example_4_3_model.x 0 Z1X2 rfl).symm, ?_⟩
  have h : (gH (0 : Fin 2)).cons.wp ((Xq 0).mul (Xq 1)) = some Z1X2 := by decide
  exact (Cons.wp_toMatrix (gH 0).models h).symm

/-- Paper, Section 4.4 (for Example 4.3): `wp(X₁) = Z₁X₂` commutes with the axis `X₁Z₂`, so `X₁ ∈ 𝓕(σ)`. -/
theorem soundness_example_X1 : exSigma.cons.wp (Xq 0) = some Z1X2 ∧ PP.Comm Z1X2 X1Z2 ∧
    Xq 0 ∈ exSigma.F := by decide

/-- Paper, Section 4.4 (for Example 4.3): every `K ∈ γ(σ)` satisfies `K Z₁X₂ = X₁ K`. -/
theorem soundness_example_gamma (K : Density 2) (hK : K ∈ exSigma.gamma) :
    K * Z1X2.toMatrix = (Xq 0).toMatrix * K :=
  hK (Xq 0) Z1X2 (by decide) (by decide)

/-- Paper, Section 4.4 (for Example 4.3): `wp(X₂) = X₂` anticommutes with the axis, so `X₂ ∉ 𝓕(σ)`. -/
theorem soundness_example_X2 : exSigma.cons.wp (Xq 1) = some (Xq 1) ∧
    (Xq 1).toMatrix * X1Z2.toMatrix = -(X1Z2.toMatrix * (Xq 1).toMatrix) ∧
    Xq 1 ∉ exSigma.F := by
  refine ⟨by decide, ?_, by decide⟩
  rw [← PP.toMatrix_mul, ← PP.toMatrix_mul,
    show (Xq (1 : Fin 2)).mul X1Z2 = (X1Z2.mul (Xq 1)).neg by decide, PP.toMatrix_neg]

/-- Section 4.4 (and 2.2): the circuit of Example 4.3 takes `Z₁X₂` to `X₁`. -/
theorem exCirc_takes : sem exCirc * Z1X2.toMatrix = (Xq 0).toMatrix * sem exCirc := by
  have h := soundness exCirc
  rw [example_4_3] at h
  exact soundness_example_gamma _ h

/-- Paper, Example 4.9: `α(H₁)` has `Z₁ ↦ X₁`, `X₁ ↦ Z₁`, and no axes. -/
theorem example_4_9_H :
    alpha (gate (gH (0 : Fin 2))) = .rel ⟨⟨![some (Zq 0), some (Xq 1)],
      ![some (Xq 0), some (Zq 1)]⟩, []⟩ := by decide

/-- Paper, Example 4.9: `α(CX_{1,2})` has `X₁X₂ ↦ X₁`, `Z₁Z₂ ↦ Z₂`, and no axes. -/
theorem example_4_9_CX :
    alpha (gate (gCX (0 : Fin 2) 1 (by decide))) = .rel ⟨⟨![some ((Xq 0).mul (Xq 1)), some (Xq 1)],
      ![some (Zq 0), some ((Zq 0).mul (Zq 1))]⟩, []⟩ := by decide

/-- Paper, Example 4.7: `α(T_q) = (Γ_I, {Z_q})`. -/
theorem example_4_7 (q : Fin n) : alpha (T q) = .rel ⟨Cons.id n, [Zq q]⟩ := rfl

/-- Paper, Example 4.7: the strings tracked by `α(T_q)` are those that commute with `Z_q`. -/
theorem example_4_7_F (q : Fin n) : (alpha (T q)).F = {P | PP.Comm P (Zq q)} := by
  ext P
  simp [ARel.F, ARel.Tracked]

/-- Paper, Example 4.7: `γ(α(T_q))` consists of the operators that commute with every Pauli
string commuting with `Z_q`. -/
theorem example_4_7_gamma (q : Fin n) :
    (alpha (T q)).gamma = {K | ∀ P : PP n, PP.Comm P (Zq q) → K * P.toMatrix = P.toMatrix * K} := by
  ext K
  simp only [alpha_rz, Dom.gamma_rel, ARel.gamma, Cons.wp_id, Option.some.injEq,
    List.mem_singleton, forall_eq, Set.mem_ofPred_eq]
  constructor
  · intro h P hP; exact h P P rfl hP
  · rintro h P B rfl hP; exact h P hP

/-- Paper, Example 4.7, on one qubit: `γ(α(T₁))` consists of the operators `K` with
`K Z₁ = Z₁ K`. -/
theorem example_4_7_one_qubit (K : Density 1) :
    K ∈ (alpha (T (0 : Fin 1))).gamma ↔ K * (Zq 0).toMatrix = (Zq 0).toMatrix * K := by
  rw [example_4_7_gamma]
  refine ⟨fun h => h (Zq 0) (by decide), fun h P hP => ?_⟩
  obtain ⟨k, L⟩ := P
  have hPk : (⟨k, L⟩ : PP 1) = (PP.single 0 (L 0)).scale k := by
    simp only [PP.single, PP.scale, zero_add, PP.mk.injEq, true_and]
    funext r
    rw [Subsingleton.elim r 0]; simp
  rw [hPk] at hP ⊢
  rw [PP.toMatrix_scale, Matrix.mul_smul, Matrix.smul_mul]
  congr 1
  cases hl : L 0 with
  | I => rw [Tab.single_I, PP.toMatrix_one, Matrix.mul_one, Matrix.one_mul]
  | Z => exact h
  | X =>
    have key : ∀ k : ZMod 4, ¬ ((PP.single (0 : Fin 1) .X).scale k).Comm (Zq 0) := by decide
    exact absurd (hl ▸ hP) (key k)
  | Y =>
    have key : ∀ k : ZMod 4, ¬ ((PP.single (0 : Fin 1) .Y).scale k).Comm (Zq 0) := by decide
    exact absurd (hl ▸ hP) (key k)

/-- Paper, Example 4.7: `T_q` and every other Z-rotation of qubit `q` are in `γ(α(T_q))`. -/
theorem example_4_7_rz (q : Fin n) (θ : ℝ) : rzMat q θ ∈ (alpha (T q)).gamma :=
  rz_mem_gamma q θ _

/-- Paper, Example 4.12: `α(CX_{1,2}) ∘ α(H₁)` is the tableau of `H₁; CX_{1,2}`, in which
the constraint `X₁X₂ ↦ X₁` becomes `Z₁X₂ ↦ X₁`, since `wp_{Γ_{H₁}}(X₁X₂) = wp(X₁) wp(X₂) = Z₁X₂`. -/
theorem example_4_12 :
    alpha (gate (gCX (0 : Fin 2) 1 (by decide))) ⊚ alpha (gate (gH 0)) = .rel ⟨exCons, []⟩ ∧
    (gH (0 : Fin 2)).cons.wp ((Xq 0).mul (Xq 1)) = some Z1X2 ∧
    (gH (0 : Fin 2)).cons.wp (Xq 0) = some (Zq 0) ∧ (gH (0 : Fin 2)).cons.wp (Xq 1) = some (Xq 1) :=
  by decide

/-- Paper, Example 4.12: composing with `α(T₂)` keeps `Γ` and adds the axis
`wp_Γ(Z₂) = X₁Z₂`, giving the relation of Example 4.3. -/
theorem example_4_12_T : alpha (T 1) ⊚ .rel ⟨exCons, []⟩ = .rel exSigma ∧
    exCons.wp (Zq 1) = some X1Z2 := by decide

/-! ### Section 2 claims -/

/-- Section 2.1: `X₁` and `Z₁` anticommute, and `X₁X₂` and `Z₁Z₂` commute. -/
theorem overview_comm : ¬ PP.Comm (Xq (0 : Fin 2)) (Zq 0) ∧
    PP.Comm ((Xq (0 : Fin 2)).mul (Xq 1)) ((Zq 0).mul (Zq 1)) := by decide

/-- Section 2.1: for the tableau of `H`, `wp(Y₁) = i Z₁X₁ = -Y₁`. -/
theorem overview_wp_Y : (gH (0 : Fin 1)).cons.wp (Yq 0) = some (Yq 0).neg := by decide

/-- Section 2.2: adding `T₁` after cut 2 adds the axis `wp(Z₁) = X₁`, giving
`𝒜 = {X₁Z₂, X₁}`. -/
theorem overview_T1 : alpha (exCirc ⨾ T 0) = .rel ⟨exCons, [X1Z2, Xq 0]⟩ := by decide

/-- Section 2.2: `T₁; H₁; T₁; H₁` and `H₁; T₁; H₁; T₁` have the same abstraction, the tableau
of the identity with the axes `{Z₁, X₁}` (as lists, in different orders). -/
theorem overview_same_abstraction :
    alpha (T 0 ⨾ gate (gH (0 : Fin 1)) ⨾ T 0 ⨾ gate (gH 0)) = .rel ⟨Cons.id 1, [Zq 0, Xq 0]⟩ ∧
    alpha (gate (gH (0 : Fin 1)) ⨾ T 0 ⨾ gate (gH 0) ⨾ T 0) = .rel ⟨Cons.id 1, [Xq 0, Zq 0]⟩ := by
  decide

end Examples

end

end PauliFold.Paper
