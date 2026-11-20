import PauliFold.Paper.Basic
import PauliFold.Join

/-!
# Paper, Section 4 (and 6.3): abstract relations, weakest preconditions, concretization,
sequential composition

One type serves Sections 4 and 6:

* `Cons n`: constraints `Γ`, one left-hand side `Γ(W)` for each generator `W = X_q, Z_q`, each known
  (`some P`) or unknown (`none`, the paper's `⋆`, Definition 6.4). The relations of Section 4
  (Definition 4.2) are those with every constraint known (`Cons.AllKnown`).
* `ARel n`: a pair `(Γ, 𝒜)` with a list of axes (the paper's finite set).
* `Dom n`: an abstract relation or `⊤` (Definition 6.4).

A *model* (Definitions 4.1 and 6.4) is a unitary `V` with `Γ(X_q) = Vᴴ X_q V` and
`Γ(Z_q) = Vᴴ Z_q V` for every *known* constraint. The paper requires every abstract relation to
have a model; here this is the predicate `Cons.Sat` (`Dom.Sat`), and the abstraction is shown to
produce only relations with a model.

The weakest precondition (Definition 4.4, Section 6.2) is the syntactic computation
`wp(i^k P₁ ⋯ Pₙ) = i^k w₁(P₁) ⋯ wₙ(Pₙ)`, unknown as soon as a needed constraint is; it is
`PTab.back`. Lemma 4.5: for every model `V`, a known `wp(P)` is `Vᴴ P V`.
-/

namespace PauliFold.Paper

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Constraints and models (Definitions 4.1 and 6.4) -/

/-- Constraints `Γ`: the left-hand sides of the constraints `Γ(X_q) ↦ X_q` and `Γ(Z_q) ↦ Z_q` for every qubit `q`; `none` is the
unknown `⋆` of Definition 6.4. -/
@[ext] structure Cons (n : ℕ) where
  /-- `Γ(X_q)`, the left-hand side of the constraint on `X_q`. -/
  x : Fin n → Option (PP n)
  /-- `Γ(Z_q)`, the left-hand side of the constraint on `Z_q`. -/
  z : Fin n → Option (PP n)

instance : DecidableEq (Cons n) := fun a b =>
  decidable_of_iff (a.x = b.x ∧ a.z = b.z) (by cases a; cases b; simp)

namespace Cons

/-- The constraints as a partial tableau, to reuse `PTab.back`. -/
def toPTab (Γ : Cons n) : PTab n := ⟨Γ.x, Γ.z, []⟩

/-- Every constraint is known (the abstract relations of Definition 4.2). -/
def AllKnown (Γ : Cons n) : Prop := ∀ q, (Γ.x q).isSome ∧ (Γ.z q).isSome

instance (Γ : Cons n) : Decidable Γ.AllKnown := by unfold AllKnown; infer_instance

/-- The constraints with the given left-hand sides, all known. -/
def ofTables (xs zs : Fin n → PP n) : Cons n := ⟨fun q => some (xs q), fun q => some (zs q)⟩

/-- `Γ_I`, the constraints of the identity: `X_q ↦ X_q`, `Z_q ↦ Z_q`. -/
def id (n : ℕ) : Cons n := ofTables (fun q => Xq q) (fun q => Zq q)

theorem ofTables_allKnown (xs zs : Fin n → PP n) : (ofTables xs zs).AllKnown :=
  fun _ => ⟨rfl, rfl⟩

end Cons

/-- Paper, Definitions 4.1 and 6.4: `V` is a *model* of `Γ`, `V ⊨ Γ`, when `V` is unitary and
`Γ(X_q) = Vᴴ X_q V`, `Γ(Z_q) = Vᴴ Z_q V` for every known constraint. -/
structure Models (V : Density n) (Γ : Cons n) : Prop where
  unitary : V ∈ Unitaries n
  x : ∀ q P, Γ.x q = some P → P.toMatrix = Vᴴ * (Xq q).toMatrix * V
  z : ∀ q P, Γ.z q = some P → P.toMatrix = Vᴴ * (Zq q).toMatrix * V

/-- `Γ` has a model. -/
def Cons.Sat (Γ : Cons n) : Prop := ∃ V, Models V Γ

theorem Models.mul_conj {V : Density n} {Γ : Cons n} (h : Models V Γ) : V * Vᴴ = 1 :=
  Unitaries.mul_conj h.unitary

theorem Models.conj_mul {V : Density n} {Γ : Cons n} (h : Models V Γ) : Vᴴ * V = 1 := h.unitary

/-- The identity is a model of `Γ_I`. -/
theorem models_id : Models (1 : Density n) (Cons.id n) :=
  ⟨Unitaries.one, fun q P h => by cases h; simp, fun q P h => by cases h; simp⟩

/-- A Clifford gate is a model of its tables. -/
theorem CliffordGate.models (G : CliffordGate n) : Models G.mat (Cons.ofTables G.xs G.zs) :=
  ⟨G.unitary, fun q P h => by cases h; exact (G.xs_spec q).symm,
    fun q P h => by cases h; exact (G.zs_spec q).symm⟩

/-! ## Weakest preconditions (Definition 4.4, Lemma 4.5) -/

namespace Cons

/-- `w_q(L)` of Definition 4.4: `w_q(I) = I`, `w_q(X) = Γ(X_q)`, `w_q(Z) = Γ(Z_q)`,
`w_q(Y) = i Γ(X_q) Γ(Z_q)`; unknown if a needed constraint is. -/
def w (Γ : Cons n) (q : Fin n) (L : Letter) : Option (PP n) := Γ.toPTab.rowOf q L

/-- Paper, Definition 4.4 (Weakest precondition):
`wp_Γ(i^k P₁ ⊗ ⋯ ⊗ Pₙ) = i^k w₁(P₁) ⋯ wₙ(Pₙ)`, unknown (`none`) as soon as a factor is
(Section 6.2). This is `PTab.back`. -/
def wp (Γ : Cons n) (P : PP n) : Option (PP n) := Γ.toPTab.back P

/-- Definition 4.4, unfolded. -/
theorem wp_def (Γ : Cons n) (P : PP n) :
    Γ.wp P = ((List.finRange n).foldr (fun r acc => mulO (Γ.w r (P.letter r)) acc)
      (some PP.one)).map (·.scale P.phase) := rfl

@[simp] theorem w_I (Γ : Cons n) (q : Fin n) : Γ.w q .I = some PP.one := rfl
@[simp] theorem w_X (Γ : Cons n) (q : Fin n) : Γ.w q .X = Γ.x q := rfl
@[simp] theorem w_Z (Γ : Cons n) (q : Fin n) : Γ.w q .Z = Γ.z q := rfl
@[simp] theorem w_Y (Γ : Cons n) (q : Fin n) :
    Γ.w q .Y = (mulO (Γ.x q) (Γ.z q)).map (·.scale 1) := rfl

theorem mulO_some_one (a : Option (PP n)) : mulO (some PP.one) a = a := by
  cases a <;> simp [mulO, PP.one_mul]

theorem mulO_one_some (a : Option (PP n)) : mulO a (some PP.one) = a := by
  cases a <;> simp [mulO, PP.mul_one]

theorem mulO_isSome {a b : Option (PP n)} : (mulO a b).isSome ↔ a.isSome ∧ b.isSome := by
  cases a <;> cases b <;> simp [mulO]

theorem backL_one (Γ : Cons n) :
    ∀ l : List (Fin n), Γ.toPTab.backL l PP.one = some PP.one
  | [] => rfl
  | r :: l => by
    show mulO (some PP.one) (Γ.toPTab.backL l PP.one) = some PP.one
    rw [backL_one Γ l]; exact mulO_some_one _

/-- Lemma 4.5: `wp(I) = I`. -/
@[simp] theorem wp_one (Γ : Cons n) : Γ.wp PP.one = some PP.one := by
  simp only [wp, PTab.back, backL_one]
  simp [PP.one]

/-- `wp(i^k P) = i^k wp(P)`. -/
theorem wp_scale (Γ : Cons n) (P : PP n) (k : ZMod 4) :
    Γ.wp (P.scale k) = (Γ.wp P).map (·.scale k) := by
  show Option.map (fun x : PP n => x.scale (P.phase + k)) (Γ.toPTab.backL _ P) =
    Option.map (fun x : PP n => x.scale k)
      (Option.map (fun x : PP n => x.scale P.phase) (Γ.toPTab.backL _ P))
  rw [Option.map_map]
  congr 1
  funext x
  exact (PP.scale_scale x _ _).symm

/-- Lemma 4.5: `wp(i^k I) = i^k I`. -/
theorem wp_phase (Γ : Cons n) (k : ZMod 4) : Γ.wp (PP.one.scale k) = some (PP.one.scale k) := by
  rw [wp_scale, wp_one]; rfl

theorem backL_single (Γ : Cons n) (q : Fin n) (L : Letter) :
    ∀ l : List (Fin n), l.Nodup →
      Γ.toPTab.backL l (PP.single q L) = if q ∈ l then Γ.w q L else some PP.one
  | [], _ => rfl
  | r :: l, hnd => by
    have hr : r ∉ l := (List.nodup_cons.1 hnd).1
    have ih := backL_single Γ q L l (List.nodup_cons.1 hnd).2
    simp only [PTab.backL, List.foldr_cons] at ih ⊢
    rw [ih, PP.single_letter]
    by_cases hrq : r = q
    · subst hrq
      simp only [if_true, List.mem_cons, true_or, if_neg hr]
      exact mulO_one_some _
    · rw [if_neg hrq]
      have : (q ∈ r :: l) ↔ q ∈ l := by simp [Ne.symm hrq]
      rw [if_congr this rfl rfl]
      exact mulO_some_one _

/-- `wp(L_q) = w_q(L)` for a single-qubit string. -/
theorem wp_single (Γ : Cons n) (q : Fin n) (L : Letter) : Γ.wp (PP.single q L) = Γ.w q L := by
  simp only [wp, PTab.back]
  rw [backL_single Γ q L _ (List.nodup_finRange n), if_pos (List.mem_finRange q)]
  show Option.map (fun x : PP n => x.scale 0) _ = _
  simp

/-- Lemma 4.5: `wp(X_q) = Γ(X_q)`. -/
@[simp] theorem wp_Xq (Γ : Cons n) (q : Fin n) : Γ.wp (Xq q) = Γ.x q := wp_single Γ q .X

/-- Lemma 4.5: `wp(Z_q) = Γ(Z_q)`. -/
@[simp] theorem wp_Zq (Γ : Cons n) (q : Fin n) : Γ.wp (Zq q) = Γ.z q := wp_single Γ q .Z

/-- Section 4.2: `wp(Y_q) = i Γ(X_q) Γ(Z_q)`. -/
theorem wp_Yq (Γ : Cons n) (q : Fin n) :
    Γ.wp (Yq q) = (mulO (Γ.x q) (Γ.z q)).map (·.scale 1) := wp_single Γ q .Y

/-! ### Which weakest preconditions are known -/

theorem w_isSome_iff (Γ : Cons n) (q : Fin n) (L : Letter) :
    (Γ.w q L).isSome ↔ match L with
      | .I => True | .X => (Γ.x q).isSome | .Z => (Γ.z q).isSome
      | .Y => (Γ.x q).isSome ∧ (Γ.z q).isSome := by
  cases L <;> simp [mulO_isSome]

theorem backL_isSome_iff (Γ : Cons n) (P : PP n) :
    ∀ l : List (Fin n), (Γ.toPTab.backL l P).isSome ↔ ∀ r ∈ l, (Γ.w r (P.letter r)).isSome
  | [] => by simp [PTab.backL]
  | r :: l => by
    simp only [PTab.backL, List.foldr_cons] at *
    rw [mulO_isSome, ← PTab.backL, backL_isSome_iff Γ P l]
    simp [w]

/-- `wp(P)` is known exactly when every needed constraint is. -/
theorem wp_isSome_iff (Γ : Cons n) (P : PP n) :
    (Γ.wp P).isSome ↔ ∀ r, (Γ.w r (P.letter r)).isSome := by
  simp only [wp, PTab.back, Option.isSome_map]
  rw [backL_isSome_iff]
  simp

theorem w_mul_isSome (Γ : Cons n) (r : Fin n) {a b : Letter} (ha : (Γ.w r a).isSome)
    (hb : (Γ.w r b).isSome) : (Γ.w r (mulLetter a b).2).isSome := by
  rw [w_isSome_iff] at ha hb ⊢
  cases a <;> cases b <;> simp_all [mulLetter]

/-- A product of strings with known weakest preconditions has a known weakest precondition. -/
theorem wp_mul_isSome (Γ : Cons n) {P Q : PP n} (hP : (Γ.wp P).isSome) (hQ : (Γ.wp Q).isSome) :
    (Γ.wp (P.mul Q)).isSome := by
  rw [wp_isSome_iff] at *
  exact fun r => w_mul_isSome Γ r (hP r) (hQ r)

/-- With every constraint known, every weakest precondition is known. -/
theorem wp_isSome_of_allKnown {Γ : Cons n} (h : Γ.AllKnown) (P : PP n) : (Γ.wp P).isSome := by
  rw [wp_isSome_iff]
  intro r
  rw [w_isSome_iff]
  cases P.letter r <;> simp [(h r).1, (h r).2]

/-! ### Weakest preconditions are conjugations (Lemma 4.5) -/

theorem w_toMatrix {Γ : Cons n} {V : Density n} (hV : Models V Γ) (r : Fin n) (L : Letter)
    {A : PP n} (h : Γ.w r L = some A) : A.toMatrix = Vᴴ * (PP.single r L).toMatrix * V := by
  cases L with
  | I => cases h; rw [Tab.single_I, PP.toMatrix_one, Matrix.mul_one, hV.conj_mul]
  | X => exact hV.x r A h
  | Z => exact hV.z r A h
  | Y =>
    simp only [w_Y, Option.map_eq_some_iff] at h
    obtain ⟨B, hB, rfl⟩ := h
    obtain ⟨x, y, hx, hy, rfl⟩ := PTab.mulO_eq_some hB
    rw [Tab.single_Y, PP.toMatrix_scale, PP.toMatrix_scale, PP.toMatrix_mul, PP.toMatrix_mul,
      hV.x r x hx, hV.z r y hy, Matrix.mul_smul, Matrix.smul_mul]
    congr 1
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc V, hV.mul_conj, Matrix.one_mul]

theorem backL_toMatrix {Γ : Cons n} {V : Density n} (hV : Models V Γ) (P : PP n) :
    ∀ l : List (Fin n), l.Nodup → ∀ {B : PP n}, Γ.toPTab.backL l P = some B →
      B.toMatrix = Vᴴ * (P.restrict l).toMatrix * V
  | [], _, B, h => by
    cases h; rw [PP.restrict_nil, PP.toMatrix_one, Matrix.mul_one, hV.conj_mul]
  | r :: l, hnd, B, h => by
    have hr : r ∉ l := (List.nodup_cons.1 hnd).1
    simp only [PTab.backL, List.foldr_cons] at h
    obtain ⟨x, y, hx, hy, rfl⟩ := PTab.mulO_eq_some h
    have ih := backL_toMatrix hV P l (List.nodup_cons.1 hnd).2 hy
    rw [PP.toMatrix_mul, ih, w_toMatrix hV r _ hx, PP.restrict_cons P r l hr, PP.toMatrix_mul]
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc V, hV.mul_conj, Matrix.one_mul]

/-- Paper, Lemma 4.5 (Weakest preconditions are conjugations): if `V ⊨ Γ` and `wp_Γ(P)` is
known, it is `Vᴴ P V`. -/
theorem wp_toMatrix {Γ : Cons n} {V : Density n} (hV : Models V Γ) {P B : PP n}
    (h : Γ.wp P = some B) : B.toMatrix = Vᴴ * P.toMatrix * V := by
  simp only [wp, PTab.back, Option.map_eq_some_iff] at h
  obtain ⟨B', hB', rfl⟩ := h
  rw [PP.toMatrix_scale, backL_toMatrix hV P _ (List.nodup_finRange n) hB']
  conv_rhs => rw [← PP.restrict_all P, PP.toMatrix_scale]
  rw [Matrix.mul_smul, Matrix.smul_mul]

/-- Paper, Lemma 4.5, for the relations of Section 4 (every constraint known): `wp_Γ(P)` is
the Pauli string `Vᴴ P V`, for every model `V`. -/
theorem wp_conj {Γ : Cons n} (hk : Γ.AllKnown) {V : Density n} (hV : Models V Γ) (P : PP n) :
    ∃ B, Γ.wp P = some B ∧ B.toMatrix = Vᴴ * P.toMatrix * V := by
  obtain ⟨B, hB⟩ := Option.isSome_iff_exists.1 (wp_isSome_of_allKnown hk P)
  exact ⟨B, hB, wp_toMatrix hV hB⟩

/-- Paper, Lemma 4.5: `wp(PQ) = wp(P) wp(Q)`, given a model. -/
theorem wp_mul {Γ : Cons n} {V : Density n} (hV : Models V Γ) {P Q A B : PP n}
    (hP : Γ.wp P = some A) (hQ : Γ.wp Q = some B) : Γ.wp (P.mul Q) = some (A.mul B) := by
  obtain ⟨C, hC⟩ := Option.isSome_iff_exists.1
    (wp_mul_isSome Γ (Option.isSome_iff_exists.2 ⟨A, hP⟩) (Option.isSome_iff_exists.2 ⟨B, hQ⟩))
  rw [hC]
  congr 1
  apply PP.toMatrix_injective
  rw [wp_toMatrix hV hC, PP.toMatrix_mul, PP.toMatrix_mul, wp_toMatrix hV hP, wp_toMatrix hV hQ]
  simp only [Matrix.mul_assoc]
  rw [← Matrix.mul_assoc V, hV.mul_conj, Matrix.one_mul]

/-- For the identity constraints, `wp(P) = P`. -/
@[simp] theorem wp_id (P : PP n) : (Cons.id n).wp P = some P := by
  obtain ⟨B, hB, hm⟩ := wp_conj (Γ := Cons.id n) (ofTables_allKnown _ _) models_id P
  rw [hB]
  congr 1
  apply PP.toMatrix_injective
  rw [hm]; simp

end Cons

/-- The Pauli strings `Gᴴ P G`: for a Clifford gate, `Gᴴ P G` is a Pauli string for every Pauli
string `P` (Section 3.2, after Definition 3.3). -/
theorem CliffordGate.conj_mem_Paulis (G : CliffordGate n) (P : PP n) :
    G.matᴴ * P.toMatrix * G.mat ∈ Paulis n := by
  obtain ⟨B, -, hB⟩ := Cons.wp_conj (Cons.ofTables_allKnown _ _) G.models P
  exact ⟨B, hB⟩

/-- Every Clifford unitary maps every Pauli string to a Pauli string (Section 3.2). -/
theorem conj_mem_Paulis_of_mem_cliffords {G : Density n} (hG : G ∈ Cliffords n) (P : PP n) :
    Gᴴ * P.toMatrix * G ∈ Paulis n :=
  (CliffordGate.ofMem G hG).conj_mem_Paulis P

/-! ## Commutation, as matrices -/

/-- Conjugation by a unitary preserves and reflects commutation. -/
theorem comm_conj_iff {V M A : Density n} (hV : V ∈ Unitaries n) :
    (Vᴴ * M * V) * (Vᴴ * A * V) = (Vᴴ * A * V) * (Vᴴ * M * V) ↔ M * A = A * M := by
  have hVV : V * Vᴴ = 1 := Unitaries.mul_conj hV
  have key : ∀ X Y : Density n, (Vᴴ * X * V) * (Vᴴ * Y * V) = Vᴴ * (X * Y) * V := by
    intro X Y
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc V, hVV, Matrix.one_mul]
  have cancel : ∀ X : Density n, V * (Vᴴ * X * V) * Vᴴ = X := by
    intro X
    simp only [← Matrix.mul_assoc, hVV, Matrix.one_mul]
    rw [Matrix.mul_assoc, hVV, Matrix.mul_one]
  rw [key, key]
  constructor
  · intro h; rw [← cancel (M * A), h, cancel]
  · intro h; rw [h]

/-- A combination `a I + b A` commutes with every Pauli string that commutes with `A`. -/
theorem combo_comm {A P : PP n} (h : PP.Comm P A) (a b : ℂ) :
    (a • (1 : Density n) + b • A.toMatrix) * P.toMatrix =
      P.toMatrix * (a • (1 : Density n) + b • A.toMatrix) := by
  rw [Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul, Matrix.smul_mul, Matrix.mul_smul,
    Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, h.toMatrix]

/-- A Z-rotation of qubit `q` commutes with every Pauli string that commutes with `Z_q`. -/
theorem rzMat_comm {q : Fin n} {P : PP n} (h : PP.Comm P (Zq q)) (θ : ℝ) :
    rzMat q θ * P.toMatrix = P.toMatrix * rzMat q θ := by
  obtain ⟨a, b, hab⟩ := rzMat_combo q θ
  rw [hab]; exact combo_comm h a b

/-! ## Matrices commuting with every Pauli string are scalars -/

/-- Flip qubit `q` of a basis state. -/
def flipQ (q : Fin n) (b : Basis n) : Basis n := Function.update b q (!b q)

theorem Zq_apply (q : Fin n) (o i : Basis n) :
    (Zq q).toMatrix o i = if o = i then (if o q then -1 else 1) else 0 := by
  simp only [PP.toMatrix, PP.single, ipow_zero, one_mul]
  by_cases h : o = i
  · subst h
    rw [if_pos rfl, Finset.prod_eq_single q]
    · cases hq : o q <;> simp [Letter.mat]
    · intro r _ hr; simp [Function.update_of_ne hr, Letter.mat]
    · simp
  · rw [if_neg h]
    obtain ⟨r, hr⟩ : ∃ r, o r ≠ i r := by
      by_contra hc; push Not at hc; exact h (funext hc)
    refine Finset.prod_eq_zero (Finset.mem_univ r) ?_
    by_cases hrq : r = q
    · subst hrq; simp [Letter.mat, hr]
    · simp [Function.update_of_ne hrq, Letter.mat, hr]

theorem Xq_apply (q : Fin n) (o i : Basis n) :
    (Xq q).toMatrix o i = if o = flipQ q i then 1 else 0 := by
  simp only [PP.toMatrix, PP.single, ipow_zero, one_mul]
  split_ifs with h
  · subst h
    refine Finset.prod_eq_one fun r _ => ?_
    by_cases hrq : r = q
    · subst hrq; simp [flipQ, Letter.mat]
    · simp [Function.update_of_ne hrq, flipQ, Letter.mat]
  · obtain ⟨r, hr⟩ : ∃ r, o r ≠ flipQ q i r := by
      by_contra hc; push Not at hc; exact h (funext hc)
    refine Finset.prod_eq_zero (Finset.mem_univ r) ?_
    by_cases hrq : r = q
    · subst hrq; simp only [flipQ, Function.update_self] at hr; simp [Letter.mat, hr]
    · simp only [flipQ, Function.update_of_ne hrq] at hr
      simp [Function.update_of_ne hrq, Letter.mat, hr]

theorem flipQ_flipQ (q : Fin n) (b : Basis n) : flipQ q (flipQ q b) = b := by
  funext r; by_cases h : r = q
  · subst h; simp [flipQ]
  · simp [flipQ, Function.update_of_ne h]

/-- **A matrix that commutes with every `X_q` and `Z_q` is a multiple of the identity.**
Commuting with the `Z_q` makes it diagonal; commuting with the `X_q` makes the diagonal
constant. -/
theorem eq_smul_one_of_comm (M : Density n)
    (hX : ∀ q, M * (Xq q).toMatrix = (Xq q).toMatrix * M)
    (hZ : ∀ q, M * (Zq q).toMatrix = (Zq q).toMatrix * M) :
    ∃ c : ℂ, M = c • (1 : Density n) := by
  -- Off-diagonal entries vanish.
  have hoff : ∀ o i : Basis n, o ≠ i → M o i = 0 := by
    intro o i hne
    obtain ⟨q, hq⟩ : ∃ q, o q ≠ i q := by
      by_contra hc; push Not at hc; exact hne (funext hc)
    have e := congrFun (congrFun (hZ q) o) i
    simp only [Matrix.mul_apply, Zq_apply] at e
    rw [Finset.sum_eq_single i (fun k _ hk => by rw [if_neg hk, mul_zero]) (by simp),
      Finset.sum_eq_single o (fun k _ hk => by rw [if_neg (Ne.symm hk), zero_mul]) (by simp),
      if_pos rfl, if_pos rfl] at e
    have : M o i * ((if i q then -1 else 1) - (if o q then -1 else 1)) = 0 := by
      rw [mul_sub, e]; ring
    rcases mul_eq_zero.1 this with h | h
    · exact h
    · exfalso; revert h hq; cases o q <;> cases i q <;> norm_num
  -- The diagonal is invariant under flips.
  have hflip : ∀ (q : Fin n) (o : Basis n), M (flipQ q o) (flipQ q o) = M o o := by
    intro q o
    have e := congrFun (congrFun (hX q) o) (flipQ q o)
    simp only [Matrix.mul_apply, Xq_apply] at e
    rw [Finset.sum_eq_single (flipQ q (flipQ q o))
        (fun k _ hk => by rw [if_neg hk, mul_zero]) (by simp),
      Finset.sum_eq_single (flipQ q o) (fun k _ hk => by
        rw [if_neg (fun h => hk (by rw [h, flipQ_flipQ])), zero_mul]) (by simp),
      if_pos rfl, if_pos (flipQ_flipQ q o).symm, flipQ_flipQ] at e
    simpa using e.symm
  -- Hence constant.
  have hdiag : ∀ S : Finset (Fin n),
      M (fun r => decide (r ∈ S)) (fun r => decide (r ∈ S)) = M (fun _ => false) (fun _ => false) := by
    intro S
    induction S using Finset.induction_on with
    | empty => simp
    | insert a S ha ih =>
      have : (fun r => decide (r ∈ insert a S)) = flipQ a (fun r => decide (r ∈ S)) := by
        funext r
        by_cases hr : r = a
        · subst hr; simp [flipQ, ha]
        · simp [flipQ, hr]
      rw [this, hflip, ih]
  refine ⟨M (fun _ => false) (fun _ => false), ?_⟩
  funext o i
  rw [Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
  by_cases h : o = i
  · subst h
    have ho : o = fun r => decide (r ∈ Finset.univ.filter fun r => o r = true) := by
      funext r; simp
    rw [if_pos rfl, mul_one, ho, hdiag]
  · rw [if_neg h, mul_zero, hoff o i h]

/-- A matrix commuting with every Pauli string is a multiple of the identity. -/
theorem eq_smul_one_of_comm_all (M : Density n)
    (h : ∀ P : PP n, M * P.toMatrix = P.toMatrix * M) : ∃ c : ℂ, M = c • (1 : Density n) :=
  eq_smul_one_of_comm M (fun _ => h _) (fun _ => h _)

/-- With every constraint known, all models are equal up to a scalar (the remark after
Definition 4.4 in earlier drafts of the paper). -/
theorem Models.unique {Γ : Cons n} (hk : Γ.AllKnown) {V W : Density n} (hV : Models V Γ)
    (hW : Models W Γ) : ∃ c : ℂ, W = c • V := by
  have key : ∀ P : PP n, (W * Vᴴ) * P.toMatrix = P.toMatrix * (W * Vᴴ) := by
    intro P
    obtain ⟨B, hB, hBV⟩ := Cons.wp_conj hk hV P
    have hBW := Cons.wp_toMatrix hW hB
    -- `Vᴴ P V = Wᴴ P W`, so `Vᴴ P = Wᴴ P W Vᴴ` and `W Vᴴ P = P W Vᴴ`.
    have e : W * (Vᴴ * P.toMatrix * V) * Vᴴ = W * (Wᴴ * P.toMatrix * W) * Vᴴ := by
      rw [← hBV, ← hBW]
    simp only [Matrix.mul_assoc] at e
    rw [hV.mul_conj, Matrix.mul_one, ← Matrix.mul_assoc W Wᴴ, hW.mul_conj, Matrix.one_mul] at e
    simp only [Matrix.mul_assoc]; exact e
  obtain ⟨c, hc⟩ := eq_smul_one_of_comm_all _ key
  refine ⟨c, ?_⟩
  calc W = W * Vᴴ * V := by rw [Matrix.mul_assoc, hV.conj_mul, Matrix.mul_one]
    _ = c • V := by rw [hc, Matrix.smul_mul, Matrix.one_mul]

/-! ## Abstract relations (Definitions 4.2 and 6.4) -/

/-- Paper, Definition 4.2 (Abstract relations), with the unknowns of Definition 6.4: a pair
`(Γ, 𝒜)` of constraints and axes. The axes are a list (the paper's finite set). Having a model
is the separate predicate `Cons.Sat`. -/
structure ARel (n : ℕ) where
  /-- The constraints `Γ`. -/
  cons : Cons n
  /-- The axes `𝒜`. -/
  axes : List (PP n)

instance : DecidableEq (ARel n) := fun a b =>
  decidable_of_iff (a.cons = b.cons ∧ a.axes = b.axes) (by cases a; cases b; simp)

/-- Paper, Definition 6.4: an abstract relation is `⊤` or a pair `(Γ, 𝒜)`. -/
inductive Dom (n : ℕ) where
  | top
  | rel (σ : ARel n)

instance : DecidableEq (Dom n) := fun a b =>
  match a, b with
  | .top, .top => isTrue rfl
  | .top, .rel _ => isFalse (by simp)
  | .rel _, .top => isFalse (by simp)
  | .rel σ, .rel τ => decidable_of_iff (σ = τ) (by simp)

instance : Top (Dom n) := ⟨.top⟩

@[simp] theorem Dom.top_def : (⊤ : Dom n) = .top := rfl

/-! ## Concretization (Definitions 4.6 and 6.5) -/

namespace ARel

/-- `P` is *tracked*: `wp(P)` is known and commutes with every axis. -/
def Tracked (σ : ARel n) (P : PP n) : Prop :=
  ∃ B, σ.cons.wp P = some B ∧ ∀ A ∈ σ.axes, PP.Comm B A

instance (σ : ARel n) (P : PP n) : Decidable (σ.Tracked P) :=
  match h : σ.cons.wp P with
  | none => isFalse (by rintro ⟨B, hB, -⟩; rw [h] at hB; cases hB)
  | some B => decidable_of_iff (∀ A ∈ σ.axes, PP.Comm B A)
      ⟨fun hc => ⟨B, h, hc⟩, fun ⟨B', hB', hc⟩ => by rw [h] at hB'; cases hB'; exact hc⟩

/-- Paper, Definitions 4.6 and 6.5: the tracked Pauli strings
`𝓕(σ) = {P | wp(P) ≠ ⋆ ∧ [wp(P), A] = 0 for all A ∈ 𝒜}`. -/
def F (σ : ARel n) : Set (PP n) := {P | σ.Tracked P}

instance (σ : ARel n) : DecidablePred (· ∈ σ.F) := fun P =>
  inferInstanceAs (Decidable (σ.Tracked P))

/-- `𝓕(σ)` with the commutator of matrices, as in the paper. -/
theorem mem_F_iff (σ : ARel n) (P : PP n) :
    P ∈ σ.F ↔ ∃ B, σ.cons.wp P = some B ∧
      ∀ A ∈ σ.axes, B.toMatrix * A.toMatrix - A.toMatrix * B.toMatrix = 0 := by
  simp only [F, Set.mem_ofPred_eq, Tracked, PP.comm_iff, sub_eq_zero]

/-- Paper, Definitions 4.6 and 6.5 (Concretization):
`γ(σ) = {K ∈ 𝒪_n | K wp(P) = P K for all P ∈ 𝓕(σ)}`, over **all** matrices. -/
def gamma (σ : ARel n) : Set (Density n) :=
  {K | ∀ P B, σ.cons.wp P = some B → (∀ A ∈ σ.axes, PP.Comm B A) →
    K * B.toMatrix = P.toMatrix * K}

/-- The concretization, stated with `𝓕(σ)`. -/
theorem mem_gamma_iff (σ : ARel n) (K : Density n) :
    K ∈ σ.gamma ↔ ∀ P ∈ σ.F, ∀ B, σ.cons.wp P = some B → K * B.toMatrix = P.toMatrix * K := by
  constructor
  · rintro hK P ⟨B, hB, hc⟩ B' hB'
    rw [hB] at hB'; cases hB'
    exact hK P B hB hc
  · intro hK P B hB hc
    exact hK P ⟨B, hB, hc⟩ B hB

theorem zero_mem_gamma (σ : ARel n) : (0 : Density n) ∈ σ.gamma := by
  intro P B _ _; simp

theorem add_mem_gamma {σ : ARel n} {K L : Density n} (hK : K ∈ σ.gamma) (hL : L ∈ σ.gamma) :
    K + L ∈ σ.gamma := by
  intro P B hB hc
  rw [Matrix.add_mul, Matrix.mul_add, hK P B hB hc, hL P B hB hc]

theorem smul_mem_gamma {σ : ARel n} (c : ℂ) {K : Density n} (hK : K ∈ σ.gamma) :
    c • K ∈ σ.gamma := by
  intro P B hB hc
  rw [Matrix.smul_mul, Matrix.mul_smul, hK P B hB hc]

end ARel

namespace Dom

/-- Paper, Definition 6.5: `γ(⊤) = 𝒪_n`, and `γ(σ)` as in Definition 4.6 otherwise. -/
def gamma : Dom n → Set (Density n)
  | top => Set.univ
  | rel σ => σ.gamma

/-- The tracked strings; `⊤` tracks nothing (a convention; the paper defines `𝓕` only on
pairs). -/
def F : Dom n → Set (PP n)
  | top => ∅
  | rel σ => σ.F

/-- The constraints of an abstract relation have a model (`⊤` trivially). -/
def Sat : Dom n → Prop
  | top => True
  | rel σ => σ.cons.Sat

@[simp] theorem gamma_top : (top : Dom n).gamma = Set.univ := rfl
@[simp] theorem gamma_rel (σ : ARel n) : (rel σ).gamma = σ.gamma := rfl
@[simp] theorem F_rel (σ : ARel n) : (rel σ).F = σ.F := rfl

theorem smul_mem_gamma {σ : Dom n} (c : ℂ) {K : Density n} (hK : K ∈ σ.gamma) :
    c • K ∈ σ.gamma := by
  cases σ with
  | top => trivial
  | rel σ => exact ARel.smul_mem_gamma c hK

theorem add_mem_gamma {σ : Dom n} {K L : Density n} (hK : K ∈ σ.gamma) (hL : L ∈ σ.gamma) :
    K + L ∈ σ.gamma := by
  cases σ with
  | top => trivial
  | rel σ => exact ARel.add_mem_gamma hK hL

end Dom

/-! ## Sequential composition (Definition 4.11, Section 6.2) -/

namespace Cons

/-- Paper, Definition 4.11: `(Γ₂ ∘ Γ₁)(W) = wp_{Γ₁}(Γ₂(W))`, unknown if `Γ₂(W)` is
(Section 6.2: `wp(⋆) = ⋆`). -/
def comp (Γ₂ Γ₁ : Cons n) : Cons n :=
  ⟨fun q => (Γ₂.x q).bind Γ₁.wp, fun q => (Γ₂.z q).bind Γ₁.wp⟩

/-- Paper, Section 4.3: if `V₁ ⊨ Γ₁` and `V₂ ⊨ Γ₂` then `V₂ V₁ ⊨ Γ₂ ∘ Γ₁`. -/
theorem _root_.PauliFold.Paper.Models.comp {Γ₁ Γ₂ : Cons n} {V₁ V₂ : Density n} (h₁ : Models V₁ Γ₁)
    (h₂ : Models V₂ Γ₂) : Models (V₂ * V₁) (Γ₂.comp Γ₁) := by
  have key : ∀ (P : PP n) (o : Option (PP n)) (L : Density n),
      (∀ Q, o = some Q → Q.toMatrix = V₂ᴴ * L * V₂) → o.bind Γ₁.wp = some P →
        P.toMatrix = (V₂ * V₁)ᴴ * L * (V₂ * V₁) := by
    intro P o L hL h
    obtain ⟨Q, hQ, hP⟩ := Option.bind_eq_some_iff.1 h
    rw [wp_toMatrix h₁ hP, hL Q hQ, conjTranspose_mul]
    simp only [Matrix.mul_assoc]
  exact ⟨Unitaries.mul h₂.unitary h₁.unitary,
    fun q P h => key P _ _ (h₂.x q) h, fun q P h => key P _ _ (h₂.z q) h⟩

theorem Sat.comp {Γ₁ Γ₂ : Cons n} (h₁ : Γ₁.Sat) (h₂ : Γ₂.Sat) : (Γ₂.comp Γ₁).Sat := by
  obtain ⟨V₁, h₁⟩ := h₁
  obtain ⟨V₂, h₂⟩ := h₂
  exact ⟨_, h₁.comp h₂⟩

theorem AllKnown.comp {Γ₁ Γ₂ : Cons n} (h₁ : Γ₁.AllKnown) (h₂ : Γ₂.AllKnown) :
    (Γ₂.comp Γ₁).AllKnown := by
  intro q
  obtain ⟨P, hP⟩ := Option.isSome_iff_exists.1 (h₂ q).1
  obtain ⟨Q, hQ⟩ := Option.isSome_iff_exists.1 (h₂ q).2
  show ((Γ₂.x q).bind Γ₁.wp).isSome ∧ ((Γ₂.z q).bind Γ₁.wp).isSome
  rw [hP, hQ]
  exact ⟨wp_isSome_of_allKnown h₁ P, wp_isSome_of_allKnown h₁ Q⟩

theorem bind_wp_isSome {Γ₁ : Cons n} (hk : Γ₁.AllKnown) {o : Option (PP n)} (h : o.isSome) :
    (o.bind Γ₁.wp).isSome := by
  obtain ⟨P, rfl⟩ := Option.isSome_iff_exists.1 h
  exact wp_isSome_of_allKnown hk P

/-- The constraint for a letter in a composition is the weakest precondition of the second
constraint's. -/
theorem w_comp {Γ₁ Γ₂ : Cons n} {V₁ : Density n} (h₁ : Models V₁ Γ₁) (r : Fin n) (L : Letter)
    {B : PP n} (h : (Γ₂.comp Γ₁).w r L = some B) :
    ∃ R, Γ₂.w r L = some R ∧ Γ₁.wp R = some B := by
  cases L with
  | I => cases h; exact ⟨PP.one, rfl, wp_one Γ₁⟩
  | X => exact Option.bind_eq_some_iff.1 h
  | Z => exact Option.bind_eq_some_iff.1 h
  | Y =>
    simp only [w_Y, Option.map_eq_some_iff] at h
    obtain ⟨C, hC, rfl⟩ := h
    obtain ⟨a, b, ha, hb, rfl⟩ := PTab.mulO_eq_some hC
    obtain ⟨Ra, hRa, hwa⟩ := Option.bind_eq_some_iff.1 ha
    obtain ⟨Rb, hRb, hwb⟩ := Option.bind_eq_some_iff.1 hb
    refine ⟨(Ra.mul Rb).scale 1, by simp [w_Y, hRa, hRb, mulO], ?_⟩
    rw [wp_scale, wp_mul h₁ hwa hwb]; rfl

theorem backL_comp {Γ₁ Γ₂ : Cons n} {V₁ : Density n} (h₁ : Models V₁ Γ₁) (P : PP n) :
    ∀ (l : List (Fin n)) {B : PP n}, (Γ₂.comp Γ₁).toPTab.backL l P = some B →
      ∃ M, Γ₂.toPTab.backL l P = some M ∧ Γ₁.wp M = some B
  | [], B, h => by cases h; exact ⟨PP.one, rfl, wp_one Γ₁⟩
  | r :: l, B, h => by
    simp only [PTab.backL, List.foldr_cons] at h ⊢
    obtain ⟨x, y, hx, hy, rfl⟩ := PTab.mulO_eq_some h
    obtain ⟨M, hM, hwM⟩ := backL_comp h₁ P l hy
    obtain ⟨R, hR, hwR⟩ := w_comp h₁ r _ hx
    refine ⟨R.mul M, ?_, wp_mul h₁ hwR hwM⟩
    rw [← PTab.backL, hM]
    show mulO (Γ₂.w r (P.letter r)) (some M) = _
    rw [hR]; rfl

/-- **The weakest precondition of a composition, in terms of the parts.** If `V₁ ⊨ Γ₁` and
`wp_{Γ₂ ∘ Γ₁}(P) = B`, then `wp_{Γ₂}(P) = M` for some `M` with `wp_{Γ₁}(M) = B`. -/
theorem wp_comp {Γ₁ Γ₂ : Cons n} {V₁ : Density n} (h₁ : Models V₁ Γ₁) {P B : PP n}
    (h : (Γ₂.comp Γ₁).wp P = some B) : ∃ M, Γ₂.wp P = some M ∧ Γ₁.wp M = some B := by
  simp only [wp, PTab.back, Option.map_eq_some_iff] at h
  obtain ⟨B', hB', rfl⟩ := h
  obtain ⟨M, hM, hwM⟩ := backL_comp h₁ P _ hB'
  refine ⟨M.scale P.phase, by simp [wp, PTab.back, hM], ?_⟩
  rw [wp_scale, hwM]; rfl

/-- With `Γ₁` total, `wp_{Γ₂ ∘ Γ₁} = wp_{Γ₁} ∘ wp_{Γ₂}`. -/
theorem wp_comp_eq {Γ₁ Γ₂ : Cons n} {V₁ : Density n} (h₁ : Models V₁ Γ₁) (hk : Γ₁.AllKnown)
    (P : PP n) : (Γ₂.comp Γ₁).wp P = (Γ₂.wp P).bind Γ₁.wp := by
  cases h : (Γ₂.comp Γ₁).wp P with
  | some B =>
    obtain ⟨M, hM, hwM⟩ := wp_comp h₁ h
    rw [hM]; exact hwM.symm
  | none =>
    cases hM : Γ₂.wp P with
    | none => rfl
    | some M =>
      exfalso
      have hs : ((Γ₂.comp Γ₁).wp P).isSome := by
        rw [wp_isSome_iff]
        have hM' := (wp_isSome_iff Γ₂ P).1 (by rw [hM]; rfl)
        intro r
        have := hM' r
        rw [w_isSome_iff] at this ⊢
        cases hL : P.letter r <;> rw [hL] at this <;> simp only at this ⊢
        all_goals first
          | trivial
          | exact bind_wp_isSome hk this
          | exact ⟨bind_wp_isSome hk this.1, bind_wp_isSome hk this.2⟩
      rw [h] at hs; cases hs

end Cons

namespace ARel

/-- Paper, Definition 4.11 with Section 6.2:
`σ₂ ∘ σ₁ = (Γ₂ ∘ Γ₁, 𝒜₁ ∪ wp_{Γ₁}(𝒜₂))`, and `⊤` if some `wp_{Γ₁}(A)`, `A ∈ 𝒜₂`, is unknown. -/
def comp (σ₂ σ₁ : ARel n) : Dom n :=
  if ∀ A ∈ σ₂.axes, (σ₁.cons.wp A).isSome then
    .rel ⟨σ₂.cons.comp σ₁.cons, σ₁.axes ++ σ₂.axes.filterMap σ₁.cons.wp⟩
  else .top

theorem comp_of_allKnown {σ₁ : ARel n} (h : σ₁.cons.AllKnown) (σ₂ : ARel n) :
    σ₂.comp σ₁ = .rel ⟨σ₂.cons.comp σ₁.cons, σ₁.axes ++ σ₂.axes.filterMap σ₁.cons.wp⟩ :=
  if_pos fun A _ => Cons.wp_isSome_of_allKnown h A

end ARel

/-- Paper, Definition 4.11 and Section 6.2: sequential composition `σ₂ ∘ σ₁` of abstract
relations, with `⊤ ∘ σ = σ ∘ ⊤ = ⊤`. -/
def Dom.comp : Dom n → Dom n → Dom n
  | .rel σ₂, .rel σ₁ => σ₂.comp σ₁
  | _, _ => .top

@[inherit_doc] scoped infixr:80 " ⊚ " => Dom.comp

@[simp] theorem Dom.top_comp (σ : Dom n) : .top ⊚ σ = .top := rfl
@[simp] theorem Dom.comp_top (σ : Dom n) : σ ⊚ .top = .top := by cases σ <;> rfl
@[simp] theorem Dom.rel_comp_rel (σ₂ σ₁ : ARel n) : .rel σ₂ ⊚ .rel σ₁ = σ₂.comp σ₁ := rfl

/-- Paper, Section 4.3: the composition of abstract relations with models has a model, so it
is an abstract relation. -/
theorem Dom.Sat.comp {σ₁ σ₂ : Dom n} (h₁ : σ₁.Sat) (h₂ : σ₂.Sat) : (σ₂ ⊚ σ₁).Sat := by
  cases σ₂ with
  | top => trivial
  | rel σ₂ =>
    cases σ₁ with
    | top => trivial
    | rel σ₁ =>
      show (σ₂.comp σ₁).Sat
      unfold ARel.comp
      split_ifs
      · exact Cons.Sat.comp h₁ h₂
      · trivial

/-- **Compositionality.** If `K₁ ∈ γ(σ₁)` and `K₂ ∈ γ(σ₂)`, then `K₂ K₁ ∈ γ(σ₂ ∘ σ₁)`. Only
`σ₁` needs a model; `K₁` and `K₂` are arbitrary matrices. -/
theorem gamma_comp {σ₁ σ₂ : Dom n} (h₁ : σ₁.Sat) {K₁ K₂ : Density n}
    (hK₁ : K₁ ∈ σ₁.gamma) (hK₂ : K₂ ∈ σ₂.gamma) : K₂ * K₁ ∈ (σ₂ ⊚ σ₁).gamma := by
  cases σ₂ with
  | top => trivial
  | rel σ₂ =>
    cases σ₁ with
    | top => trivial
    | rel σ₁ =>
      obtain ⟨V₁, hV₁⟩ := h₁
      show K₂ * K₁ ∈ (σ₂.comp σ₁).gamma
      unfold ARel.comp
      split_ifs with hax
      · intro P B hB hc
        obtain ⟨M, hM, hwM⟩ := Cons.wp_comp hV₁ hB
        have hcM : ∀ A ∈ σ₂.axes, PP.Comm M A := by
          intro A hA
          obtain ⟨A', hA'⟩ := Option.isSome_iff_exists.1 (hax A hA)
          have hBA : PP.Comm B A' :=
            hc A' (List.mem_append_right _ (List.mem_filterMap.2 ⟨A, hA, hA'⟩))
          rw [PP.comm_iff] at hBA ⊢
          rw [Cons.wp_toMatrix hV₁ hwM, Cons.wp_toMatrix hV₁ hA'] at hBA
          exact (comm_conj_iff hV₁.unitary).1 hBA
        have e₂ := hK₂ P M hM hcM
        have e₁ := hK₁ M B hwM (fun A hA => hc A (List.mem_append_left _ hA))
        rw [Matrix.mul_assoc, e₁, ← Matrix.mul_assoc, e₂, Matrix.mul_assoc]
      · trivial

end

end PauliFold.Paper
