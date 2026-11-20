import PauliFold.TableauContainment
import PauliFold.FoldRule

/-!
# Paper, Section 3: Pauli strings, Clifford gates, Z-rotations, circuits

The paper-faithful development of "Pauli Abstractions for Quantum Circuits" starts here.

* Pauli matrices are `Letter.mat` (Definition 3.1), Pauli strings `i^k P₁ ⊗ ⋯ ⊗ Pₙ` are `PP n`
  (Definition 3.2) with `PP.toMatrix`. Qubits are `Fin n`, so the paper's qubit `q` is the
  Lean qubit `q - 1`.
* A Clifford gate (Definition 3.3) is a `CliffordGate`: a unitary matrix together with its
  tableau, the Pauli strings `Gᴴ X_q G` and `Gᴴ Z_q G`, as witnesses.
* The standard gates `H, S, S†, X, Z, CX, CZ` (Section 3.2) use tzap's `gateUnitary` matrices;
  their tables are the ones of Figure 2.
* Z-rotations `R_Z(θ) = e^{-iθZ/2}` with a real angle (Definition 3.4); how they conjugate `X` (Section 2.2).
* Circuits (Definition 3.5), their semantics (Definition 3.6) and equivalence up to a global
  phase (Definition 3.7).
-/

namespace PauliFold.PP

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-- `-P`. -/
def neg (P : PP n) : PP n := P.scale 2

theorem toMatrix_neg (P : PP n) : P.neg.toMatrix = -P.toMatrix := by
  rw [neg, PP.toMatrix_scale, ipow_two, neg_one_smul]

theorem scale_scale (P : PP n) (a b : ZMod 4) : (P.scale a).scale b = P.scale (a + b) := by
  simp [PP.scale, add_assoc]

@[simp] theorem scale_zero (P : PP n) : P.scale 0 = P := by
  cases P; simp [PP.scale]

@[simp] theorem letter_scale (P : PP n) (k : ZMod 4) : (P.scale k).letter = P.letter := rfl

theorem one_mul (P : PP n) : PP.one.mul P = P := by
  apply PP.ext'
  · simp [PP.mul, PP.one, PP.mulLetter_I_left]
  · intro r; simp [PP.mul, PP.one, PP.mulLetter_I_left]

theorem mul_one (P : PP n) : P.mul PP.one = P := by
  apply PP.ext'
  · simp [PP.mul, PP.one, PP.mulLetter_I_right]
  · intro r; simp [PP.mul, PP.one, PP.mulLetter_I_right]

/-- The sum over the qubits of the phases of a product with a single-qubit string. -/
theorem sum_mulLetter_single_right (P : PP n) (r : Fin n) (L : Letter) :
    ∑ s, (mulLetter (P.letter s) ((PP.single r L).letter s)).1 = (mulLetter (P.letter r) L).1 := by
  rw [Finset.sum_eq_single r]
  · simp [PP.single_letter]
  · intro s _ hs; simp [PP.single_letter, hs, PP.mulLetter_I_right]
  · simp

theorem sum_mulLetter_single_left (P : PP n) (r : Fin n) (L : Letter) :
    ∑ s, (mulLetter ((PP.single r L).letter s) (P.letter s)).1 = (mulLetter L (P.letter r)).1 := by
  rw [Finset.sum_eq_single r]
  · simp [PP.single_letter]
  · intro s _ hs; simp [PP.single_letter, hs, PP.mulLetter_I_left]
  · simp

/-- `P` commutes with a single-qubit string `L_r` exactly when the two letters at `r` commute. -/
theorem comm_single (P : PP n) (r : Fin n) (L : Letter) :
    PP.Comm P (PP.single r L) ↔ (mulLetter (P.letter r) L).1 = (mulLetter L (P.letter r)).1 := by
  unfold PP.Comm
  simp only [PP.mul, sum_mulLetter_single_right, sum_mulLetter_single_left]
  simp only [PP.single, zero_add, add_zero]
  constructor
  · intro h; linear_combination h
  · intro h; rw [h]

/-- **The matrix of a Pauli string determines the string.** -/
theorem toMatrix_injective : Function.Injective (PP.toMatrix (n := n)) := by
  intro P Q h
  have hc : ∀ A : PP n, PP.Comm P A ↔ PP.Comm Q A := fun A => by
    rw [PP.comm_iff, PP.comm_iff, h]
  have hL : P.letter = Q.letter := by
    funext r
    have hx := hc (PP.single r .X)
    have hz := hc (PP.single r .Z)
    rw [PP.comm_single, PP.comm_single] at hx hz
    revert hx hz
    cases P.letter r <;> cases Q.letter r <;> decide
  have hm := h
  rw [PP.toMatrix_letters P, PP.toMatrix_letters Q, hL] at hm
  have hne := PP.toMatrix_base_ne_zero Q.letter
  have hp : ipow P.phase = ipow Q.phase := by
    by_contra hne'
    apply hne
    have : (ipow P.phase - ipow Q.phase) • (⟨0, Q.letter⟩ : PP n).toMatrix = 0 := by
      rw [sub_smul, hm, sub_self]
    exact (smul_eq_zero.1 this).resolve_left (sub_ne_zero.2 hne')
  exact PP.ext' (ipow_inj hp) (fun r => congrFun hL r)

theorem toMatrix_inj {P Q : PP n} : P.toMatrix = Q.toMatrix ↔ P = Q :=
  PP.toMatrix_injective.eq_iff

theorem zmod4_two_mul (x : ZMod 4) (h : 2 * x = 0) : x = 0 ∨ x = 2 := by
  revert x; decide

/-- Section 3.1: **two Pauli strings either commute or anticommute.** -/
theorem comm_or_anticomm (P Q : PP n) :
    P.toMatrix * Q.toMatrix = Q.toMatrix * P.toMatrix ∨
      P.toMatrix * Q.toMatrix = -(Q.toMatrix * P.toMatrix) := by
  set D : ZMod 4 := ∑ r, ((mulLetter (P.letter r) (Q.letter r)).1 -
    (mulLetter (Q.letter r) (P.letter r)).1)
  have hD : 2 * D = 0 := by
    simp only [D, Finset.mul_sum]
    refine Finset.sum_eq_zero fun r _ => ?_
    cases P.letter r <;> cases Q.letter r <;> decide
  have hPQ : P.mul Q = (Q.mul P).scale D := by
    apply PP.ext'
    · simp only [PP.mul, PP.scale, D, Finset.sum_sub_distrib]; ring
    · intro r; simp [PP.mul, PP.scale, PP.mulLetter_snd_comm]
  rw [← PP.toMatrix_mul, ← PP.toMatrix_mul, hPQ, PP.toMatrix_scale]
  rcases zmod4_two_mul D hD with h | h <;> rw [h]
  · left; simp
  · right; simp

end

end PauliFold.PP

namespace PauliFold.Paper

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Pauli matrices (Definition 3.1) -/

/-- Paper, Definition 3.1 (Pauli matrices), entry by entry. `Letter.mat L o i` is the entry in
row `o`, column `i`, with `false = |0⟩` first. -/
theorem pauli_matrices (o i : Bool) :
    Letter.I.mat o i = (if o = i then 1 else 0) ∧
    Letter.X.mat o i = (if o = i then 0 else 1) ∧
    Letter.Y.mat o i = (if o = i then 0 else if o then Complex.I else -Complex.I) ∧
    Letter.Z.mat o i = (if o = i then (if o then -1 else 1) else 0) := by
  cases o <;> cases i <;> simp [Letter.mat]

/-- Paper, Definition 3.1: each Pauli matrix is its own inverse, `P² = I`. -/
theorem letter_mat_sq (L : Letter) (o i : Bool) :
    ∑ c : Bool, L.mat o c * L.mat c i = Letter.I.mat o i := by
  rw [mulLetter_mat]
  cases L <;> simp [mulLetter]

/-! ## Pauli strings (Definition 3.2) -/

/-- Paper, Definition 3.2: the set `𝒫_n` of Pauli strings, as matrices. -/
def Paulis (n : ℕ) : Set (Density n) := Set.range (PP.toMatrix (n := n))

/-- `X_q`: the letter `X` at qubit `q`, the identity elsewhere. -/
abbrev Xq (q : Fin n) : PP n := PP.single q .X
/-- `Y_q`. -/
abbrev Yq (q : Fin n) : PP n := PP.single q .Y
/-- `Z_q`. -/
abbrev Zq (q : Fin n) : PP n := PP.single q .Z

/-- The product of Pauli strings is a Pauli string (Section 3.1). -/
theorem mul_mem_Paulis {M N : Density n} (hM : M ∈ Paulis n) (hN : N ∈ Paulis n) :
    M * N ∈ Paulis n := by
  obtain ⟨P, rfl⟩ := hM
  obtain ⟨Q, rfl⟩ := hN
  exact ⟨P.mul Q, PP.toMatrix_mul P Q⟩

/-- An operator `K` *takes* `P` to `Q` when `K P = Q K` (Section 3.1). -/
def Takes (K : Density n) (P Q : PP n) : Prop := K * P.toMatrix = Q.toMatrix * K

/-! ## Unitaries -/

/-- The unitaries `𝒰_n`. -/
def Unitaries (n : ℕ) : Set (Density n) := {U | Uᴴ * U = 1}

theorem Unitaries.mul_conj {U : Density n} (h : U ∈ Unitaries n) : U * Uᴴ = 1 :=
  mul_eq_one_comm.mp h

theorem Unitaries.mul {U V : Density n} (hU : U ∈ Unitaries n) (hV : V ∈ Unitaries n) :
    U * V ∈ Unitaries n := by
  show (U * V)ᴴ * (U * V) = 1
  rw [conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc Uᴴ, hU, Matrix.one_mul, hV]

theorem Unitaries.one : (1 : Density n) ∈ Unitaries n := by
  show (1 : Density n)ᴴ * 1 = 1; simp

/-! ## Clifford unitaries and Clifford gates (Definition 3.3) -/

/-- Paper, Definition 3.3 (Clifford unitaries): the unitaries `G` with `Gᴴ X_q G` and
`Gᴴ Z_q G` Pauli strings for every qubit `q`. -/
def Cliffords (n : ℕ) : Set (Density n) :=
  {G | G ∈ Unitaries n ∧ ∀ q : Fin n,
    Gᴴ * (Xq q).toMatrix * G ∈ Paulis n ∧ Gᴴ * (Zq q).toMatrix * G ∈ Paulis n}

/-- A Clifford gate: a matrix `G ∈ 𝒞_n` together with its tableau, the Pauli strings
`xs q = Gᴴ X_q G` and `zs q = Gᴴ Z_q G`, as witnesses of Definition 3.3. -/
structure CliffordGate (n : ℕ) where
  /-- The matrix of the gate. -/
  mat : Density n
  unitary : mat ∈ Unitaries n
  /-- `xs q = Gᴴ X_q G`. -/
  xs : Fin n → PP n
  /-- `zs q = Gᴴ Z_q G`. -/
  zs : Fin n → PP n
  xs_spec : ∀ q, matᴴ * (Xq q).toMatrix * mat = (xs q).toMatrix
  zs_spec : ∀ q, matᴴ * (Zq q).toMatrix * mat = (zs q).toMatrix

namespace CliffordGate

/-- The matrix of a Clifford gate is a Clifford unitary. -/
theorem mat_mem_cliffords (G : CliffordGate n) : G.mat ∈ Cliffords n :=
  ⟨G.unitary, fun q => ⟨⟨G.xs q, (G.xs_spec q).symm⟩, ⟨G.zs q, (G.zs_spec q).symm⟩⟩⟩

/-- Every Clifford unitary is the matrix of a Clifford gate (its tables chosen). -/
def ofMem (G : Density n) (hG : G ∈ Cliffords n) : CliffordGate n where
  mat := G
  unitary := hG.1
  xs q := ((hG.2 q).1).choose
  zs q := ((hG.2 q).2).choose
  xs_spec q := ((hG.2 q).1).choose_spec.symm
  zs_spec q := ((hG.2 q).2).choose_spec.symm

@[simp] theorem ofMem_mat (G : Density n) (hG : G ∈ Cliffords n) : (ofMem G hG).mat = G := rfl

theorem mul_conj (G : CliffordGate n) : G.mat * G.matᴴ = 1 := Unitaries.mul_conj G.unitary

/-- The identity gate. -/
def one (n : ℕ) : CliffordGate n where
  mat := 1
  unitary := Unitaries.one
  xs q := Xq q
  zs q := Zq q
  xs_spec q := by simp
  zs_spec q := by simp

end CliffordGate

/-! ## The standard Clifford gates (Section 3.2, Figure 2) -/

theorem XZ_eq_Y (q : Fin n) : (Xq q).mul (Zq q) = (Yq q).scale 3 := by
  rw [PP.mul_single_eq]
  apply PP.ext'
  · rfl
  · intro r; by_cases h : r = q <;> simp [h, PP.scale, PP.single_letter, mulLetter]

/-- A standard gate built from a Clifford `CTGate`: its matrix is tzap's `gateUnitary`, and its
tables are given explicitly and checked against the pull-back lemmas of `Tableau.lean`. -/
def ofCT (g : CTGate n) (hg : g.IsClifford) (xs zs : Fin n → PP n)
    (hx : ∀ r, pullG g (Xq r) = xs r) (hz : ∀ r, pullG g (Zq r) = zs r) : CliffordGate n where
  mat := gateUnitary n g.toGate
  unitary := gate_adj_mul g
  xs := xs
  zs := zs
  xs_spec r := by rw [← hx, pullG_toMatrix g hg]
  zs_spec r := by rw [← hz, pullG_toMatrix g hg]

/-- `H_q`. Figure 2: `Z_q ↦ X_q`, `X_q ↦ Z_q`. -/
def gH (q : Fin n) : CliffordGate n :=
  ofCT (.h q) trivial (fun r => if r = q then Zq q else Xq r)
    (fun r => if r = q then Xq q else Zq r) (fun r => pull_h_X r q) (fun r => pull_h_Z r q)

/-- `S_q`. Figure 2: `-Y_q ↦ X_q`. -/
def gS (q : Fin n) : CliffordGate n :=
  ofCT (.s q) trivial (fun r => if r = q then (Yq q).neg else Xq r) (fun r => Zq r)
    (fun r => by
      rw [pull_s_X]; split_ifs
      · rw [XZ_eq_Y, PP.scale_scale, PP.neg]; rfl
      · rfl)
    (fun r => pull_s_Z r q)

/-- `S†_q`: `Y_q ↦ X_q`. -/
def gSdg (q : Fin n) : CliffordGate n :=
  ofCT (.sdg q) trivial (fun r => if r = q then Yq q else Xq r) (fun r => Zq r)
    (fun r => by
      rw [pull_sdg_X]; split_ifs
      · rw [XZ_eq_Y, PP.scale_scale]; exact PP.scale_zero _
      · rfl)
    (fun r => pull_sdg_Z r q)

/-- `X_q`. Figure 2: `-Z_q ↦ Z_q`. -/
def gX (q : Fin n) : CliffordGate n :=
  ofCT (.x q) trivial (fun r => Xq r) (fun r => if r = q then (Zq q).neg else Zq r)
    (fun r => pull_x_X r q) (fun r => pull_x_Z r q)

/-- `Z_q`. Figure 2: `-X_q ↦ X_q`. -/
def gZ (q : Fin n) : CliffordGate n :=
  ofCT (.z q) trivial (fun r => if r = q then (Xq q).neg else Xq r) (fun r => Zq r)
    (fun r => pull_z_X r q) (fun r => pull_z_Z r q)

/-- `CX_{c,t}`. Figure 2: `X_c X_t ↦ X_c`, `Z_c Z_t ↦ Z_t`. -/
def gCX (c t : Fin n) (hne : c ≠ t) : CliffordGate n :=
  ofCT (.cx c t hne) trivial (fun r => if r = c then (Xq c).mul (Xq t) else Xq r)
    (fun r => if r = t then (Zq c).mul (Zq t) else Zq r)
    (fun r => pull_cx_X r c t hne) (fun r => pull_cx_Z r c t hne)

/-- `CZ_{c,t}`: `X_c Z_t ↦ X_c`, `Z_c X_t ↦ X_t`. -/
def gCZ (c t : Fin n) (hne : c ≠ t) : CliffordGate n :=
  ofCT (.cz c t hne) trivial
    (fun r => if r = c then (Xq c).mul (Zq t) else if r = t then (Zq c).mul (Xq t) else Xq r)
    (fun r => Zq r) (fun r => pull_cz_X r c t hne) (fun r => pull_cz_Z r c t hne)

/-! ## Z-rotations (Definition 3.4) -/

/-- Paper, Definition 3.4: `R_Z(θ) = e^{-iθZ/2}` on qubit `q`, for a real angle `θ`:
`diag(e^{-iθ/2}, e^{iθ/2})` on `q`, the identity elsewhere. -/
def rzMat (q : Fin n) (θ : ℝ) : Density n :=
  embed1 n (diag2 (Complex.exp (-(θ / 2 : ℝ) * Complex.I)) (Complex.exp ((θ / 2 : ℝ) * Complex.I)))
    q.val

/-- A diagonal one-qubit gate is `((a + b)/2) I + ((a - b)/2) Z_q`. -/
theorem embed1_diag2_ab (a b : ℂ) (q : Fin n) :
    embed1 n (diag2 a b) q.val = ((a + b) / 2) • (1 : Density n) + ((a - b) / 2) • (Zq q).toMatrix := by
  rw [PP.toMatrix_single]
  funext o i
  rw [Matrix.add_apply, Matrix.smul_apply, Matrix.smul_apply, embed1_apply_of_lt _ q.isLt,
    embed1_apply_of_lt _ q.isLt, Matrix.one_apply, Fin.eta]
  by_cases hc : ∀ r : Fin n, (r : ℕ) ≠ q.val → o r = i r
  · rw [if_pos hc, if_pos hc]
    have hoi : o = i ↔ o q = i q := ⟨fun h => h ▸ rfl, fun h => funext fun r => by
      by_cases hr : r = q
      · subst hr; exact h
      · exact hc r (Fin.val_ne_of_ne hr)⟩
    simp only [hoi, smul_eq_mul]
    cases o q <;> cases i q <;> simp [diag2, Letter.mat] <;> ring
  · rw [if_neg hc, if_neg hc, if_neg (fun h => hc fun r _ => by rw [h])]
    simp

/-- Paper, Definition 3.4: `R_Z(θ) = cos(θ/2) I - i sin(θ/2) Z`. -/
theorem rzMat_eq (q : Fin n) (θ : ℝ) :
    rzMat q θ = (Real.cos (θ / 2) : ℂ) • (1 : Density n) -
      (Complex.I * Real.sin (θ / 2)) • (Zq q).toMatrix := by
  have hneg : Complex.exp (-((θ / 2 : ℝ) : ℂ) * Complex.I) =
      Complex.cos ((θ / 2 : ℝ) : ℂ) - Complex.sin ((θ / 2 : ℝ) : ℂ) * Complex.I := by
    rw [Complex.exp_mul_I, Complex.cos_neg, Complex.sin_neg]; ring
  have hpos := Complex.exp_mul_I ((θ / 2 : ℝ) : ℂ)
  have ha : (Complex.exp (-((θ / 2 : ℝ) : ℂ) * Complex.I) +
      Complex.exp (((θ / 2 : ℝ) : ℂ) * Complex.I)) / 2 = (Real.cos (θ / 2) : ℂ) := by
    rw [hneg, hpos, Complex.ofReal_cos]; ring
  have hb : (Complex.exp (-((θ / 2 : ℝ) : ℂ) * Complex.I) -
      Complex.exp (((θ / 2 : ℝ) : ℂ) * Complex.I)) / 2 = -(Complex.I * Real.sin (θ / 2)) := by
    rw [hneg, hpos, Complex.ofReal_sin]; ring
  rw [rzMat, embed1_diag2_ab, ha, hb, neg_smul, ← sub_eq_add_neg]

/-- A Z-rotation is a combination `a I + b Z_q`. -/
theorem rzMat_combo (q : Fin n) (θ : ℝ) :
    ∃ a b : ℂ, rzMat q θ = a • (1 : Density n) + b • (Zq q).toMatrix :=
  ⟨_, _, by rw [rzMat_eq, sub_eq_add_neg, ← neg_smul]⟩

theorem Zq_sq (q : Fin n) : (Zq q).toMatrix * (Zq q).toMatrix = 1 := single_Z_sq q

theorem Zq_conjTranspose (q : Fin n) : (Zq q).toMatrixᴴ = (Zq q).toMatrix := by
  rw [PP.toMatrix_single, embed1_conjTranspose]
  congr 1
  funext o i
  cases o <;> cases i <;> simp [Letter.mat]

/-- Rotations about `Z_q` add their angles. -/
theorem rzMat_mul (q : Fin n) (θ φ : ℝ) : rzMat q θ * rzMat q φ = rzMat q (θ + φ) := by
  rw [rzMat, rzMat, rzMat, embed1_mul]
  congr 1
  funext o i
  cases o <;> cases i <;> simp [mul2, diag2, ← Complex.exp_add] <;> ring_nf

@[simp] theorem rzMat_zero (q : Fin n) : rzMat q 0 = 1 := by
  rw [rzMat_eq]; simp

theorem rzMat_conjTranspose (q : Fin n) (θ : ℝ) : (rzMat q θ)ᴴ = rzMat q (-θ) := by
  rw [rzMat_eq, rzMat_eq, conjTranspose_sub, conjTranspose_smul, conjTranspose_smul,
    conjTranspose_one, Zq_conjTranspose]
  have h1 : star ((Real.cos (θ / 2) : ℝ) : ℂ) = ((Real.cos (-θ / 2) : ℝ) : ℂ) := by
    rw [Complex.star_def, Complex.conj_ofReal, neg_div, Real.cos_neg]
  have h2 : star (Complex.I * ((Real.sin (θ / 2) : ℝ) : ℂ)) =
      Complex.I * ((Real.sin (-θ / 2) : ℝ) : ℂ) := by
    rw [Complex.star_def, map_mul, Complex.conj_ofReal, Complex.conj_I, neg_div, Real.sin_neg]
    push_cast; ring
  rw [h1, h2]

theorem rzMat_unitary (q : Fin n) (θ : ℝ) : rzMat q θ ∈ Unitaries n := by
  show (rzMat q θ)ᴴ * rzMat q θ = 1
  rw [rzMat_conjTranspose, rzMat_mul, neg_add_cancel, rzMat_zero]

/-- tzap's `T` gate `diag(1, e^{iπ/4})` is `R_Z(π/4)` up to the global phase `e^{iπ/8}`. -/
theorem gateUnitary_t (q : Fin n) :
    gateUnitary n (CTGate.t q).toGate = Complex.exp ((Real.pi / 8 : ℝ) * Complex.I) • rzMat q (Real.pi / 4) := by
  show embed1 n (diag2 1 (ep (1 / 4))) q.val = _
  rw [rzMat, ← embed1_smul _ _ _ q.isLt]
  congr 1
  funext o i
  cases o <;> cases i <;> simp [diag2, ep, ← Complex.exp_add] <;> (ring_nf; try simp)

/-- tzap's `T†` gate is `R_Z(-π/4)` up to the global phase `e^{-iπ/8}`. -/
theorem gateUnitary_tdg (q : Fin n) :
    gateUnitary n (CTGate.tdg q).toGate =
      Complex.exp ((-(Real.pi / 8) : ℝ) * Complex.I) • rzMat q (-(Real.pi / 4)) := by
  show embed1 n (diag2 1 (ep (-1 / 4))) q.val = _
  rw [rzMat, ← embed1_smul _ _ _ q.isLt]
  congr 1
  funext o i
  cases o <;> cases i <;> simp [diag2, ep, ← Complex.exp_add] <;> (ring_nf; try simp)

/-! ## Z-rotations conjugate `X` out of the Pauli strings (Section 2.2) -/

theorem ZX_eq_Y (q : Fin n) : (Zq q).mul (Xq q) = (Yq q).scale 1 := by
  rw [PP.mul_single_eq]
  apply PP.ext'
  · rfl
  · intro r; by_cases h : r = q <;> simp [h, PP.scale, PP.single_letter, mulLetter]

/-- A combination `a I + b Z` conjugates `X`, on one qubit. -/
theorem rot_conj_X (a b a' b' : ℂ) :
    (a • (1 : Density 1) + b • (Zq 0).toMatrix) * (Xq (0 : Fin 1)).toMatrix * (a' • 1 + b' • (Zq 0).toMatrix) =
      (a * a' - b * b') • (Xq (0 : Fin 1)).toMatrix + (Complex.I * (b * a' - a * b')) • (Yq 0).toMatrix := by
  have hXZ : (Xq (0 : Fin 1)).toMatrix * (Zq 0).toMatrix = (-Complex.I) • (Yq 0).toMatrix := by
    rw [← PP.toMatrix_mul, show (Xq (0 : Fin 1)).mul (Zq 0) = (Yq 0).scale 3 by decide,
      PP.toMatrix_scale, ipow_three]
  have hZX : (Zq (0 : Fin 1)).toMatrix * (Xq 0).toMatrix = Complex.I • (Yq 0).toMatrix := by
    rw [← PP.toMatrix_mul, show (Zq (0 : Fin 1)).mul (Xq 0) = (Yq 0).scale 1 by decide,
      PP.toMatrix_scale, ipow_one]
  have hYZ : (Yq (0 : Fin 1)).toMatrix * (Zq 0).toMatrix = Complex.I • (Xq 0).toMatrix := by
    rw [← PP.toMatrix_mul, show (Yq (0 : Fin 1)).mul (Zq 0) = (Xq 0).scale 1 by decide,
      PP.toMatrix_scale, ipow_one]
  simp only [add_mul, mul_add, smul_mul_assoc, mul_smul_comm, Matrix.one_mul, Matrix.mul_one,
    hXZ, hZX, hYZ, smul_smul]
  match_scalars <;> (ring_nf; try simp only [Complex.I_sq]; try ring)

/-- Paper, Section 2.2: on one qubit, `R_Z(θ) X R_Z(θ)† = cos θ X + sin θ Y`. -/
theorem overview_rz_conj_X (θ : ℝ) :
    rzMat (0 : Fin 1) θ * (Xq 0).toMatrix * (rzMat 0 θ)ᴴ =
      (Real.cos θ : ℂ) • (Xq (0 : Fin 1)).toMatrix + (Real.sin θ : ℂ) • (Yq 0).toMatrix := by
  rw [rzMat_conjTranspose, rzMat_eq, rzMat_eq, sub_eq_add_neg, sub_eq_add_neg, ← neg_smul,
    ← neg_smul, rot_conj_X]
  have hc : Real.cos θ = Real.cos (θ / 2) ^ 2 - Real.sin (θ / 2) ^ 2 := by
    conv_lhs => rw [show θ = 2 * (θ / 2) by ring]
    rw [Real.cos_two_mul]; nlinarith [Real.sin_sq_add_cos_sq (θ / 2)]
  have hs : Real.sin θ = 2 * Real.sin (θ / 2) * Real.cos (θ / 2) := by
    conv_lhs => rw [show θ = 2 * (θ / 2) by ring]
    rw [Real.sin_two_mul]
  rw [hc, hs, neg_div, Real.cos_neg, Real.sin_neg]
  congr 1 <;> congr 1 <;> (push_cast; ring_nf; try simp [Complex.I_sq]; try ring)

/-- Paper, Section 2.2: for `θ = π/4` (the `T` gate up to a global phase), the conjugate of `X`
is `(X + Y)/√2`. -/
theorem overview_T_conj_X :
    rzMat (0 : Fin 1) (Real.pi / 4) * (Xq 0).toMatrix * (rzMat 0 (Real.pi / 4))ᴴ =
      (1 / Real.sqrt 2 : ℂ) • ((Xq (0 : Fin 1)).toMatrix + (Yq 0).toMatrix) := by
  rw [overview_rz_conj_X, Real.cos_pi_div_four, Real.sin_pi_div_four, smul_add]
  have h2 : (Real.sqrt 2 : ℂ) ≠ 0 := by exact_mod_cast (Real.sqrt_pos.2 two_pos).ne'
  have : ((Real.sqrt 2 / 2 : ℝ) : ℂ) = 1 / (Real.sqrt 2 : ℂ) := by
    push_cast
    field_simp
    rw [← Complex.ofReal_pow, Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)]
    norm_num
  rw [this]

/-- Paper, Section 2.2: `cos θ X + sin θ Y` is a Pauli string only when `θ` is a multiple
of `π/2`. -/
theorem overview_rz_conj_X_pauli (θ : ℝ)
    (h : (Real.cos θ : ℂ) • (Xq (0 : Fin 1)).toMatrix + (Real.sin θ : ℂ) • (Yq 0).toMatrix ∈
      Paulis 1) : ∃ m : ℤ, θ = m * (Real.pi / 2) := by
  obtain ⟨P, hP⟩ := h
  have e := congrFun (congrFun hP (fun _ => true)) (fun _ => false)
  simp [PP.toMatrix, PP.single, Letter.mat] at e
  have hexp := Complex.exp_mul_I (θ : ℂ)
  have hk4 : ipow P.phase ^ 4 = 1 := by
    unfold ipow; rw [← pow_mul, mul_comm, pow_mul, Complex.I_pow_four, one_pow]
  have h4 : Complex.exp (θ * Complex.I) ^ 4 = 1 := by
    rw [hexp, ← e]
    cases hL : P.letter 0 <;> simp [mul_pow, hk4, Complex.I_pow_four]
    all_goals
      exfalso
      have := Complex.exp_ne_zero (θ * Complex.I)
      rw [hexp, ← e, hL] at this
      simp at this
  rw [← Complex.exp_nat_mul, Complex.exp_eq_one_iff] at h4
  obtain ⟨m, hm⟩ := h4
  refine ⟨m, ?_⟩
  have := congrArg Complex.im hm
  simp at this
  nlinarith [Real.pi_pos]

/-! ## Circuits (Definitions 3.5–3.7) -/

/-- Paper, Definition 3.5 (Circuits): `C ::= G | R_Z^q(θ) | C ; C` with `G ∈ 𝒞_n`. -/
inductive Circuit (n : ℕ) where
  | gate (G : CliffordGate n)
  | rz (q : Fin n) (θ : ℝ)
  /-- `C₁ ; C₂`: run `C₁`, then `C₂`. -/
  | seq (C₁ C₂ : Circuit n)

@[inherit_doc] scoped infixl:65 " ⨾ " => Circuit.seq

namespace Circuit

/-- Paper, Definition 3.6 (Semantics): `⟦G⟧ = G`, `⟦R_Z^q(θ)⟧ = R_Z^q(θ)`,
`⟦C₁; C₂⟧ = ⟦C₂⟧ ⟦C₁⟧`. -/
def sem : Circuit n → Density n
  | gate G => G.mat
  | rz q θ => rzMat q θ
  | seq C₁ C₂ => sem C₂ * sem C₁

@[simp] theorem sem_gate (G : CliffordGate n) : sem (gate G) = G.mat := rfl
@[simp] theorem sem_rz (q : Fin n) (θ : ℝ) : sem (rz q θ) = rzMat q θ := rfl
@[simp] theorem sem_seq (C₁ C₂ : Circuit n) : sem (C₁ ⨾ C₂) = sem C₂ * sem C₁ := rfl

/-- The semantics of a circuit is unitary. -/
theorem sem_unitary : ∀ C : Circuit n, sem C ∈ Unitaries n
  | gate G => G.unitary
  | rz q θ => rzMat_unitary q θ
  | seq C₁ C₂ => Unitaries.mul (sem_unitary C₂) (sem_unitary C₁)

/-- `T_q`, the rotation `R_Z^q(π/4)` (the paper's `T`, up to a global phase). -/
abbrev T (q : Fin n) : Circuit n := rz q (Real.pi / 4)

/-- `T†_q = R_Z^q(-π/4)`. -/
abbrev Tdg (q : Fin n) : Circuit n := rz q (-(Real.pi / 4))

/-- A circuit made only of Clifford gates. -/
def IsClifford : Circuit n → Prop
  | gate _ => True
  | rz _ _ => False
  | seq C₁ C₂ => IsClifford C₁ ∧ IsClifford C₂

end Circuit

/-- Paper, Definition 3.7: `U ≡ U'` when `U = e^{iφ} U'` for some real `φ`. -/
def EquivP (U U' : Density n) : Prop := ∃ φ : ℝ, U = Complex.exp (φ * Complex.I) • U'

theorem EquivP.refl (U : Density n) : EquivP U U := ⟨0, by simp⟩

theorem EquivP.symm {U U' : Density n} (h : EquivP U U') : EquivP U' U := by
  obtain ⟨φ, rfl⟩ := h
  refine ⟨-φ, ?_⟩
  rw [smul_smul, ← Complex.exp_add]
  simp

theorem EquivP.trans {U U' U'' : Density n} (h : EquivP U U') (h' : EquivP U' U'') :
    EquivP U U'' := by
  obtain ⟨φ, rfl⟩ := h
  obtain ⟨ψ, rfl⟩ := h'
  refine ⟨φ + ψ, ?_⟩
  rw [smul_smul, ← Complex.exp_add]
  push_cast; ring_nf

theorem EquivP.of_eq {U U' : Density n} (h : U = U') : EquivP U U' := h ▸ EquivP.refl U

/-- Paper, Definition 3.7: circuits are equivalent when their semantics are. -/
def Circuit.Equiv (C₁ C₂ : Circuit n) : Prop := EquivP C₁.sem C₂.sem

end

end PauliFold.Paper
