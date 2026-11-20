import PauliFold.Paper.Circuits

/-!
# Paper, Section 6: handling control flow

* Projectors `Π^q_b = (I + (-1)^b Z_q)/2` (`proj`), Definition 6.1 (`alphaMeas`, `alphaReset`),
  Lemma 6.2 (`lemma_6_2_*`), Example 6.3 (`Examples.example_6_3*`).
* Definitions 6.4 and 6.5 (`⊤`, composition with `⊤`) are the core's `Dom`, `Dom.gamma` and
  `Dom.comp` (`Paper/Relation.lean`).
* Section 4.2 (partial tableaux): the example without a constraint on `Z₂` (`Examples.unknowns_example`)
  and the claim that Lemma 4.5 still holds where the weakest precondition is defined (`lemma_4_5_unknowns`, which is
  the core's `Cons.wp_toMatrix`).
* Definition 6.6 (join, `Dom.join`, notation `σ₁ ⊔ σ₂`), Lemma 6.7 (`lemma_6_7`), the join has a
  model (`Dom.Sat.join_left`, `Dom.Sat.join_right`), Example 6.8 (`Examples.example_6_8*`).
* Definition 6.9 (programs, `Program`), Definition 6.10 (semantics, `Program.psem`), and
  `𝒪⟦C⟧ = {⟦C⟧}` for circuits (`Program.psem_ofCircuit`).
* Definition 6.11 (abstraction of programs, `Program.alpha`, with `loopIter`), Lemma 6.12
  (`lemma_6_12`), Example 6.13 (`Examples.example_6_13*`), Theorem 6.14 (`theorem_6_14`).
* The disjunctive refinement of Section 8.4 (RQ2, "Configurations"): branches and loops keep up to
  `N` abstract relations before joining them (`Disj.run`, soundness `Disj.program_sound`).

Qubits are 0-indexed: the paper's qubits 1 and 2 are the Lean qubits `0` and `1`.
-/

namespace PauliFold.Paper

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Measurement and reset (Section 6.1) -/

/-- Paper, Section 6.1: the projector `Π^q_b = ½ (I + (-1)^b Z_q)` onto outcome `b` of measuring
qubit `q`. -/
def proj (q : Fin n) (b : Fin 2) : Density n :=
  (1 / 2 : ℂ) • ((1 : Density n) + ((-1 : ℂ) ^ (b : ℕ)) • (Zq q).toMatrix)

/-- `Π^q_b` is a combination `a I + c Z_q`. -/
theorem proj_combo (q : Fin n) (b : Fin 2) :
    proj q b = (1 / 2 : ℂ) • (1 : Density n) + ((1 / 2 : ℂ) * (-1 : ℂ) ^ (b : ℕ)) • (Zq q).toMatrix := by
  rw [proj, smul_add, smul_smul]

/-- `Π^q_b` is a projector: `Π^q_b Π^q_b = Π^q_b`. -/
theorem proj_mul_self (q : Fin n) (b : Fin 2) : proj q b * proj q b = proj q b := by
  have hs : (-1 : ℂ) ^ ((b : ℕ) * 2) = 1 := by rw [mul_comm, pow_mul]; simp
  simp only [proj, Matrix.smul_mul, Matrix.mul_smul, Matrix.add_mul, Matrix.mul_add,
    Matrix.one_mul, Matrix.mul_one, smul_smul, Zq_sq]
  rw [smul_add, smul_add, smul_add, smul_add, smul_add]
  match_scalars <;> (ring_nf; try rw [hs]; try ring)

/-- The two outcomes are complete: `Π^q_0 + Π^q_1 = I`. -/
theorem proj_add (q : Fin n) : proj q 0 + proj q 1 = 1 := by
  simp only [proj, Fin.val_zero, Fin.val_one, pow_zero, pow_one, one_smul, neg_smul, ← smul_add]
  rw [add_add_add_comm, add_neg_cancel, add_zero, ← two_smul ℂ, smul_smul]
  norm_num

/-- Paper, Definition 6.1: `α(measure q) = (Γ_I, {Z_q})`. -/
def alphaMeas (q : Fin n) : Dom n := .rel ⟨Cons.id n, [Zq q]⟩

/-- Paper, Definition 6.1: `α(reset q) = (Γ_I, {X_q, Z_q})`. -/
def alphaReset (q : Fin n) : Dom n := .rel ⟨Cons.id n, [Xq q, Zq q]⟩

/-- The concretization of `(Γ_I, 𝒜)`: the operators commuting with every Pauli string that
commutes with all the axes. -/
theorem mem_gamma_id_iff (axes : List (PP n)) (K : Density n) :
    K ∈ (Dom.rel ⟨Cons.id n, axes⟩).gamma ↔
      ∀ P : PP n, (∀ A ∈ axes, PP.Comm P A) → K * P.toMatrix = P.toMatrix * K := by
  simp only [Dom.gamma_rel, ARel.gamma, Cons.wp_id, Option.some.injEq, Set.mem_ofPred_eq]
  constructor
  · intro h P hP; exact h P P rfl hP
  · rintro h P B rfl hP; exact h P hP

/-- **Paper, Lemma 6.2** (measurement): `Π^q_b ∈ γ(α(measure q))` for `b ∈ {0, 1}`. -/
theorem lemma_6_2_measure (q : Fin n) (b : Fin 2) : proj q b ∈ (alphaMeas q).gamma := by
  rw [alphaMeas, mem_gamma_id_iff]
  intro P hP
  rw [proj_combo]
  exact (combo_comm (hP _ (List.mem_singleton_self _)) _ _).symm ▸ rfl

/-- **Paper, Lemma 6.2** (reset, outcome 0): `Π^q_0 ∈ γ(α(reset q))`. -/
theorem lemma_6_2_reset_zero (q : Fin n) : proj q 0 ∈ (alphaReset q).gamma := by
  rw [alphaReset, mem_gamma_id_iff]
  intro P hP
  rw [proj_combo]
  exact combo_comm (hP _ (by simp)) _ _

/-- **Paper, Lemma 6.2** (reset, outcome 1): `X_q Π^q_1 ∈ γ(α(reset q))`. -/
theorem lemma_6_2_reset_one (q : Fin n) : (Xq q).toMatrix * proj q 1 ∈ (alphaReset q).gamma := by
  rw [alphaReset, mem_gamma_id_iff]
  intro P hP
  have hX := (hP (Xq q) (by simp)).toMatrix
  have hZ : proj q 1 * P.toMatrix = P.toMatrix * proj q 1 := by
    rw [proj_combo]; exact combo_comm (hP _ (by simp)) _ _
  rw [Matrix.mul_assoc, hZ, ← Matrix.mul_assoc, ← hX, Matrix.mul_assoc]

/-- Paper, Lemma 6.2, all three statements. -/
theorem lemma_6_2 (q : Fin n) :
    (∀ b : Fin 2, proj q b ∈ (alphaMeas q).gamma) ∧ proj q 0 ∈ (alphaReset q).gamma ∧
      (Xq q).toMatrix * proj q 1 ∈ (alphaReset q).gamma :=
  ⟨lemma_6_2_measure q, lemma_6_2_reset_zero q, lemma_6_2_reset_one q⟩

/-! ## Partial tableaux (Sections 4.1 and 4.2) -/

/-- Paper, Definition 4.4: the weakest precondition is undefined as soon as one of the constraints it
uses is dropped (and only then). -/
theorem wp_eq_none_iff (Γ : Cons n) (P : PP n) :
    Γ.wp P = none ↔ ∃ r, Γ.w r (P.letter r) = none := by
  rw [← Option.not_isSome_iff_eq_none, Cons.wp_isSome_iff]
  push Not
  simp only [ne_eq, Option.not_isSome_iff_eq_none]

/-- Paper, Definition 4.11: composition keeps a dropped constraint dropped. -/
@[simp] theorem wp_unknown (Γ : Cons n) : (none : Option (PP n)).bind Γ.wp = none := rfl

/-- **Paper, Lemma 4.5, where the weakest precondition is defined.** If `V ⊨ Γ` and `wp_Γ(P)` is defined, then `wp_Γ(P) = Vᴴ P V`. (This is the core's `Cons.wp_toMatrix`.) -/
theorem lemma_4_5_unknowns {Γ : Cons n} {V : Density n} (hV : Models V Γ) {P B : PP n}
    (h : Γ.wp P = some B) : B.toMatrix = Vᴴ * P.toMatrix * V :=
  Cons.wp_toMatrix hV h

/-- Paper, Definition 6.5: `σ₂ ∘ σ₁ = ⊤` when `wp_{Γ₁}(A)` is undefined for some `A ∈ 𝒜₂`, and `⊤ ∘ σ = σ ∘ ⊤ = ⊤`
(the core's `ARel.comp` and `Dom.comp`). -/
theorem comp_eq_top {σ₁ σ₂ : ARel n} (h : ∃ A ∈ σ₂.axes, σ₁.cons.wp A = none) :
    (.rel σ₂ : Dom n) ⊚ .rel σ₁ = .top ∧ ∀ σ : Dom n, .top ⊚ σ = .top ∧ σ ⊚ .top = .top := by
  refine ⟨?_, fun σ => ⟨Dom.top_comp σ, Dom.comp_top σ⟩⟩
  rw [Dom.rel_comp_rel, ARel.comp, if_neg]
  obtain ⟨A, hA, hn⟩ := h
  exact fun hall => by have := hall A hA; rw [hn] at this; cases this

/-! ## Weaker constraints and the order on abstract relations -/

namespace Cons

/-- `Γ'` is weaker than `Γ`: every known constraint of `Γ'` is the same constraint of `Γ`. -/
def Le (Γ Γ' : Cons n) : Prop :=
  (∀ q, Γ'.x q = none ∨ Γ'.x q = Γ.x q) ∧ (∀ q, Γ'.z q = none ∨ Γ'.z q = Γ.z q)

instance (Γ Γ' : Cons n) : Decidable (Γ.Le Γ') := by unfold Le; infer_instance

theorem Le.refl (Γ : Cons n) : Γ.Le Γ := ⟨fun _ => Or.inr rfl, fun _ => Or.inr rfl⟩

theorem Le.trans {Γ₁ Γ₂ Γ₃ : Cons n} (h₁ : Γ₁.Le Γ₂) (h₂ : Γ₂.Le Γ₃) : Γ₁.Le Γ₃ := by
  refine ⟨fun q => ?_, fun q => ?_⟩
  · rcases h₂.1 q with h | h
    · exact Or.inl h
    · rcases h₁.1 q with h' | h'
      · exact Or.inl (h.trans h')
      · exact Or.inr (h.trans h')
  · rcases h₂.2 q with h | h
    · exact Or.inl h
    · rcases h₁.2 q with h' | h'
      · exact Or.inl (h.trans h')
      · exact Or.inr (h.trans h')

theorem Le.x {Γ Γ' : Cons n} (h : Γ.Le Γ') {q : Fin n} {P : PP n} (hP : Γ'.x q = some P) :
    Γ.x q = some P := by
  rcases h.1 q with h' | h'
  · rw [h'] at hP; cases hP
  · rw [← h', hP]

theorem Le.z {Γ Γ' : Cons n} (h : Γ.Le Γ') {q : Fin n} {P : PP n} (hP : Γ'.z q = some P) :
    Γ.z q = some P := by
  rcases h.2 q with h' | h'
  · rw [h'] at hP; cases hP
  · rw [← h', hP]

theorem w_of_le {Γ Γ' : Cons n} (h : Γ.Le Γ') (r : Fin n) (L : Letter) {A : PP n}
    (hA : Γ'.w r L = some A) : Γ.w r L = some A := by
  cases L with
  | I => exact hA
  | X => exact h.x hA
  | Z => exact h.z hA
  | Y =>
    simp only [w_Y, Option.map_eq_some_iff] at hA ⊢
    obtain ⟨B, hB, rfl⟩ := hA
    obtain ⟨x, y, hx, hy, rfl⟩ := PTab.mulO_eq_some hB
    exact ⟨x.mul y, by rw [h.x hx, h.z hy]; rfl, rfl⟩

theorem backL_of_le {Γ Γ' : Cons n} (h : Γ.Le Γ') (P : PP n) :
    ∀ (l : List (Fin n)) {B : PP n}, Γ'.toPTab.backL l P = some B → Γ.toPTab.backL l P = some B
  | [], B, hB => hB
  | r :: l, B, hB => by
    simp only [PTab.backL, List.foldr_cons] at hB ⊢
    obtain ⟨x, y, hx, hy, rfl⟩ := PTab.mulO_eq_some hB
    have ih := backL_of_le h P l (B := y) hy
    rw [← PTab.backL, ih]
    show mulO (Γ.w r (P.letter r)) (some y) = _
    rw [w_of_le h r _ hx]; rfl

/-- A known weakest precondition through weaker constraints is the same through stronger ones. -/
theorem wp_of_le {Γ Γ' : Cons n} (h : Γ.Le Γ') {P B : PP n} (hB : Γ'.wp P = some B) :
    Γ.wp P = some B := by
  simp only [wp, PTab.back, Option.map_eq_some_iff] at hB ⊢
  obtain ⟨B', hB', rfl⟩ := hB
  exact ⟨B', backL_of_le h P _ hB', rfl⟩

/-- A model of constraints is a model of weaker constraints. -/
theorem _root_.PauliFold.Paper.Models.of_le {Γ Γ' : Cons n} {V : Density n} (hV : Models V Γ)
    (h : Γ.Le Γ') : Models V Γ' :=
  ⟨hV.unitary, fun q P hP => hV.x q P (h.x hP), fun q P hP => hV.z q P (h.z hP)⟩

end Cons

/-- `σ ⊑ τ`: `τ`'s constraints are weaker and its axes include `σ`'s. -/
def ARel.Le (σ τ : ARel n) : Prop := σ.cons.Le τ.cons ∧ ∀ A ∈ σ.axes, A ∈ τ.axes

instance (σ τ : ARel n) : Decidable (σ.Le τ) := by unfold ARel.Le; infer_instance

/-- `σ ⊑ τ` on abstract relations with `⊤`; everything is below `⊤`. -/
def Dom.Le : Dom n → Dom n → Prop
  | _, .top => True
  | .top, .rel _ => False
  | .rel σ, .rel τ => σ.Le τ

instance : (σ τ : Dom n) → Decidable (σ.Le τ)
  | _, .top => isTrue trivial
  | .top, .rel _ => isFalse id
  | .rel σ, .rel τ => inferInstanceAs (Decidable (σ.Le τ))

namespace Dom

theorem Le.refl : ∀ σ : Dom n, σ.Le σ
  | .top => trivial
  | .rel _ => ⟨Cons.Le.refl _, fun _ h => h⟩

theorem Le.trans : ∀ {σ₁ σ₂ σ₃ : Dom n}, σ₁.Le σ₂ → σ₂.Le σ₃ → σ₁.Le σ₃
  | _, _, .top, _, _ => trivial
  | .top, .rel _, .rel _, h, _ => h.elim
  | _, .top, .rel _, _, h => h.elim
  | .rel _, .rel _, .rel _, h₁, h₂ => ⟨h₁.1.trans h₂.1, fun A hA => h₂.2 A (h₁.2 A hA)⟩

/-- A weaker abstract relation describes more operators. -/
theorem gamma_mono : ∀ {σ τ : Dom n}, σ.Le τ → σ.gamma ⊆ τ.gamma
  | _, .top, _ => fun _ _ => trivial
  | .top, .rel _, h => h.elim
  | .rel σ, .rel τ, h => by
    intro K hK P B hB hc
    exact hK P B (Cons.wp_of_le h.1 hB) (fun A hA => hc A (h.2 A hA))

/-- A weaker abstract relation still has a model. -/
theorem Sat.mono : ∀ {σ τ : Dom n}, σ.Le τ → σ.Sat → τ.Sat
  | _, .top, _, _ => trivial
  | .top, .rel _, h, _ => h.elim
  | .rel _, .rel _, h, ⟨V, hV⟩ => ⟨V, hV.of_le h.1⟩

end Dom

/-! ## Join (Definition 6.6, Lemma 6.7) -/

/-- Paper, Definition 6.6, for one generator `W`: `Γ₁(W)` if `Γ₁(W) = Γ₂(W)`, and dropped otherwise. -/
def joinV (a b : Option (PP n)) : Option (PP n) := if a = b then a else none

/-- Paper, Definition 6.6: the join of constraints, generator by generator. -/
def Cons.join (Γ₁ Γ₂ : Cons n) : Cons n :=
  ⟨fun q => joinV (Γ₁.x q) (Γ₂.x q), fun q => joinV (Γ₁.z q) (Γ₂.z q)⟩

/-- The union `𝒜₁ ∪ 𝒜₂` of two lists of axes: `𝒜₁`, then the axes of `𝒜₂` not in `𝒜₁`. -/
def axUnion (A₁ A₂ : List (PP n)) : List (PP n) := A₁ ++ A₂.filter (fun A => A ∉ A₁)

theorem mem_axUnion {A : PP n} {A₁ A₂ : List (PP n)} : A ∈ axUnion A₁ A₂ ↔ A ∈ A₁ ∨ A ∈ A₂ := by
  simp only [axUnion, List.mem_append, List.mem_filter, decide_eq_true_eq]
  tauto

theorem axUnion_of_subset {A₁ A₂ : List (PP n)} (h : ∀ A ∈ A₂, A ∈ A₁) : axUnion A₁ A₂ = A₁ := by
  rw [axUnion, List.filter_eq_nil_iff.2 (fun A hA => by simpa using h A hA), List.append_nil]

/-- Paper, Definition 6.6: `σ₁ ⊔ σ₂ = (Γ, 𝒜₁ ∪ 𝒜₂)` with `Γ` the join of the constraints. -/
def ARel.join (σ₁ σ₂ : ARel n) : ARel n := ⟨σ₁.cons.join σ₂.cons, axUnion σ₁.axes σ₂.axes⟩

/-- **Paper, Definition 6.6 (Join)**: `⊤ ⊔ σ = σ ⊔ ⊤ = ⊤`, and `ARel.join` on pairs. -/
def Dom.join : Dom n → Dom n → Dom n
  | .rel σ₁, .rel σ₂ => .rel (σ₁.join σ₂)
  | _, _ => .top

/-- The notation `σ₁ ⊔ σ₂` for the join. -/
instance : Max (Dom n) := ⟨Dom.join⟩

@[simp] theorem Dom.top_join (σ : Dom n) : (.top : Dom n) ⊔ σ = .top := rfl
@[simp] theorem Dom.join_top (σ : Dom n) : σ ⊔ (.top : Dom n) = .top := by cases σ <;> rfl
@[simp] theorem Dom.rel_join_rel (σ₁ σ₂ : ARel n) :
    (.rel σ₁ : Dom n) ⊔ .rel σ₂ = .rel (σ₁.join σ₂) := rfl

theorem Cons.le_join_left (Γ₁ Γ₂ : Cons n) : Γ₁.Le (Γ₁.join Γ₂) := by
  refine ⟨fun q => ?_, fun q => ?_⟩ <;> simp only [join, joinV] <;> split_ifs <;> simp

theorem Cons.le_join_right (Γ₁ Γ₂ : Cons n) : Γ₂.Le (Γ₁.join Γ₂) := by
  refine ⟨fun q => ?_, fun q => ?_⟩ <;> simp only [join, joinV] <;> split_ifs with h <;> simp [h]

theorem Dom.le_join_left : ∀ σ₁ σ₂ : Dom n, σ₁.Le (σ₁ ⊔ σ₂)
  | .top, _ => trivial
  | .rel _, .top => trivial
  | .rel _, .rel _ => ⟨Cons.le_join_left _ _, fun _ hA => mem_axUnion.2 (Or.inl hA)⟩

theorem Dom.le_join_right : ∀ σ₁ σ₂ : Dom n, σ₂.Le (σ₁ ⊔ σ₂)
  | .top, _ => by rw [Dom.top_join]; trivial
  | .rel _, .top => trivial
  | .rel _, .rel _ => ⟨Cons.le_join_right _ _, fun _ hA => mem_axUnion.2 (Or.inr hA)⟩

/-- **Paper, Lemma 6.7**: `γ(σ₁) ∪ γ(σ₂) ⊆ γ(σ₁ ⊔ σ₂)`. -/
theorem lemma_6_7 (σ₁ σ₂ : Dom n) : σ₁.gamma ∪ σ₂.gamma ⊆ (σ₁ ⊔ σ₂).gamma := by
  rintro K (hK | hK)
  · exact Dom.gamma_mono (Dom.le_join_left σ₁ σ₂) hK
  · exact Dom.gamma_mono (Dom.le_join_right σ₁ σ₂) hK

/-- Paper, Section 6.3: **the join has a model** (every model of `σ₁` is one), so it is an
abstract relation. -/
theorem Dom.Sat.join_left {σ₁ : Dom n} (h : σ₁.Sat) (σ₂ : Dom n) : (σ₁ ⊔ σ₂).Sat :=
  Dom.Sat.mono (Dom.le_join_left σ₁ σ₂) h

/-- The join has a model (every model of `σ₂` is one). -/
theorem Dom.Sat.join_right (σ₁ : Dom n) {σ₂ : Dom n} (h : σ₂.Sat) : (σ₁ ⊔ σ₂).Sat :=
  Dom.Sat.mono (Dom.le_join_right σ₁ σ₂) h

/-- Every model of `Γ₁` is a model of `Γ₁ ⊔ Γ₂`. -/
theorem Models.join_left {Γ₁ : Cons n} {V : Density n} (h : Models V Γ₁) (Γ₂ : Cons n) :
    Models V (Γ₁.join Γ₂) := h.of_le (Cons.le_join_left _ _)

/-! ## Programs and their semantics (Definitions 6.9 and 6.10) -/

/-- Paper, Definition 6.9 (Programs):
`C ::= G | R_Z^q(θ) | measure q | reset q | C ; C | if (*) {C} else {C} | while (*) {C}`. -/
inductive Program (n : ℕ) where
  | gate (G : CliffordGate n)
  | rz (q : Fin n) (θ : ℝ)
  | meas (q : Fin n)
  | reset (q : Fin n)
  /-- `C₁ ; C₂`. -/
  | seq (C₁ C₂ : Program n)
  /-- `if (*) {C₁} else {C₂}`. -/
  | choice (C₁ C₂ : Program n)
  /-- `while (*) {C}`. -/
  | loop (C : Program n)

namespace Program

/-- The identity program `I` (the identity gate). -/
def skip (n : ℕ) : Program n := gate (CliffordGate.one n)

/-- `if (*) {C}`: a conditional without an `else` branch has the identity as its `else`
branch. -/
def ifThen (C : Program n) : Program n := choice C (skip n)

/-- `𝒪⟦C⟧^k` of Definition 6.10: `S⁰ = {I}`, `S^{k+1} = {K' K | K ∈ S^k, K' ∈ S}`. -/
def powSem (S : Set (Density n)) : ℕ → Set (Density n)
  | 0 => {1}
  | k + 1 => {M | ∃ K ∈ powSem S k, ∃ K' ∈ S, M = K' * K}

/-- **Paper, Definition 6.10 (Semantics of programs)**: the set `𝒪⟦C⟧ ⊆ 𝒪_n` of the operators
of the runs of `C`. -/
def psem : Program n → Set (Density n)
  | gate G => {G.mat}
  | rz q θ => {rzMat q θ}
  | meas q => {proj q 0, proj q 1}
  | reset q => {proj q 0, (Xq q).toMatrix * proj q 1}
  | seq C₁ C₂ => {M | ∃ K₁ ∈ psem C₁, ∃ K₂ ∈ psem C₂, M = K₂ * K₁}
  | choice C₁ C₂ => psem C₁ ∪ psem C₂
  | loop C => ⋃ k, powSem (psem C) k

/-- A circuit as a program. -/
def ofCircuit : Circuit n → Program n
  | .gate G => gate G
  | .rz q θ => rz q θ
  | .seq C₁ C₂ => seq (ofCircuit C₁) (ofCircuit C₂)

/-- Paper, Section 6.4: for a circuit `C`, `𝒪⟦C⟧ = {⟦C⟧}`. -/
theorem psem_ofCircuit : ∀ C : Circuit n, psem (ofCircuit C) = {C.sem}
  | .gate _ => rfl
  | .rz _ _ => rfl
  | .seq C₁ C₂ => by
    ext M
    simp only [ofCircuit, psem, psem_ofCircuit C₁, psem_ofCircuit C₂, Set.mem_singleton_iff,
      Set.mem_ofPred_eq, exists_eq_left, Circuit.sem_seq]

end Program

/-! ## The abstraction of programs (Definition 6.11) -/

/-- `α(I) = (Γ_I, ∅)`, the abstraction of the identity. -/
def alphaI (n : ℕ) : Dom n := .rel ⟨Cons.id n, []⟩

/-- One loop iteration of Definition 6.11: `σ ↦ σ ⊔ (τ ∘ σ)`, with `τ = α(C)`. -/
def loopStep (τ σ : Dom n) : Dom n := σ ⊔ (τ ⊚ σ)

/-- Paper, Definition 6.11: the iterates `σ₀ = α(I)`, `σ_{i+1} = σ_i ⊔ (τ ∘ σ_i)` of a loop
whose body has the abstraction `τ`. -/
def loopIter (τ : Dom n) : ℕ → Dom n
  | 0 => alphaI n
  | i + 1 => loopIter τ i ⊔ (τ ⊚ loopIter τ i)

namespace Program

/-- **Paper, Definition 6.11 (Abstraction of programs)**: Definitions 4.8, 4.10, 4.13 and 6.1,
with `α(if (*) {C₁} else {C₂}) = α(C₁) ⊔ α(C₂)` and `α(while (*) {C}) = σ_{2n+1}`. -/
def alpha : Program n → Dom n
  | gate G => .rel ⟨G.cons, []⟩
  | rz q _ => .rel ⟨Cons.id n, [Zq q]⟩
  | meas q => alphaMeas q
  | reset q => alphaReset q
  | seq C₁ C₂ => alpha C₂ ⊚ alpha C₁
  | choice C₁ C₂ => alpha C₁ ⊔ alpha C₂
  | loop C => loopIter (alpha C) (2 * n + 1)

@[simp] theorem alpha_seq (C₁ C₂ : Program n) : alpha (seq C₁ C₂) = alpha C₂ ⊚ alpha C₁ := rfl
@[simp] theorem alpha_choice (C₁ C₂ : Program n) :
    alpha (choice C₁ C₂) = alpha C₁ ⊔ alpha C₂ := rfl
@[simp] theorem alpha_loop (C : Program n) : alpha (loop C) = loopIter (alpha C) (2 * n + 1) := rfl

theorem alpha_skip : alpha (skip n) = alphaI n := rfl

/-- On circuits, the abstraction of programs is the abstraction of circuits. -/
theorem alpha_ofCircuit : ∀ C : Circuit n, alpha (ofCircuit C) = C.alpha
  | .gate _ => rfl
  | .rz _ _ => rfl
  | .seq C₁ C₂ => by rw [ofCircuit, alpha_seq, alpha_ofCircuit C₁, alpha_ofCircuit C₂]; rfl

end Program

/-! ## Lemma 6.12: the loop iterates stabilize -/

namespace Cons

/-- The number of known constraints. -/
def known (Γ : Cons n) : ℕ := ∑ q, ((Γ.x q).isSome.toNat + (Γ.z q).isSome.toNat)

theorem known_le (Γ : Cons n) : Γ.known ≤ 2 * n := by
  calc Γ.known ≤ ∑ _q : Fin n, 2 := Finset.sum_le_sum fun q _ => by
        have h1 := Bool.toNat_le (Γ.x q).isSome
        have h2 := Bool.toNat_le (Γ.z q).isSome
        omega
    _ = 2 * n := by simp [mul_comm]

theorem joinV_toNat_le (a b : Option (PP n)) : (joinV a b).isSome.toNat ≤ a.isSome.toNat := by
  unfold joinV; split_ifs <;> simp

theorem joinV_toNat_lt {a b : Option (PP n)} (h : joinV a b ≠ a) :
    (joinV a b).isSome.toNat < a.isSome.toNat := by
  unfold joinV at *
  split_ifs at * with hab
  · exact absurd rfl h
  · cases a
    · exact absurd rfl h
    · simp

theorem known_join_le (Γ Γ' : Cons n) : (Γ.join Γ').known ≤ Γ.known :=
  Finset.sum_le_sum fun _ _ => Nat.add_le_add (joinV_toNat_le _ _) (joinV_toNat_le _ _)

/-- The join can only drop constraints: if it changes the constraints, it decreases their
number. -/
theorem known_join_lt {Γ Γ' : Cons n} (h : Γ.join Γ' ≠ Γ) : (Γ.join Γ').known < Γ.known := by
  have : ∃ q, joinV (Γ.x q) (Γ'.x q) ≠ Γ.x q ∨ joinV (Γ.z q) (Γ'.z q) ≠ Γ.z q := by
    by_contra hc
    push Not at hc
    exact h (Cons.ext (funext fun q => (hc q).1) (funext fun q => (hc q).2))
  obtain ⟨q, hq⟩ := this
  apply Finset.sum_lt_sum (fun q _ => Nat.add_le_add (joinV_toNat_le _ _) (joinV_toNat_le _ _))
  refine ⟨q, Finset.mem_univ _, ?_⟩
  rcases hq with hq | hq
  · exact Nat.add_lt_add_of_lt_of_le (joinV_toNat_lt hq) (joinV_toNat_le _ _)
  · exact Nat.add_lt_add_of_le_of_lt (joinV_toNat_le _ _) (joinV_toNat_lt hq)

end Cons

@[simp] theorem loopStep_top (τ : Dom n) : loopStep τ .top = .top := rfl

theorem loopIter_succ (τ : Dom n) (i : ℕ) : loopIter τ (i + 1) = loopStep τ (loopIter τ i) := rfl

theorem loopIter_eq_iterate (τ : Dom n) : ∀ i, loopIter τ i = (loopStep τ)^[i] (alphaI n)
  | 0 => rfl
  | i + 1 => by rw [loopIter_succ, loopIter_eq_iterate τ i, Function.iterate_succ_apply']

/-- One iteration from a relation: the constraints are joined with something, and if they did not
change, the next iteration changes nothing (the axes have caught up). -/
theorem loopStep_rel {τ : Dom n} {s s' : ARel n} (h : loopStep τ (.rel s) = .rel s') :
    (∃ Γ, s'.cons = s.cons.join Γ) ∧ (s'.cons = s.cons → loopStep τ (.rel s') = .rel s') := by
  cases τ with
  | top => cases h
  | rel t =>
    simp only [loopStep, Dom.rel_comp_rel, ARel.comp] at h
    split_ifs at h with hc
    · rw [Dom.rel_join_rel] at h
      cases h
      refine ⟨⟨_, rfl⟩, fun hcons => ?_⟩
      simp only [loopStep, Dom.rel_comp_rel, ARel.comp]
      simp only [ARel.join] at hcons ⊢
      simp only [hcons]
      rw [if_pos hc, Dom.rel_join_rel]
      simp only [ARel.join, hcons]
      congr 2
      apply axUnion_of_subset
      intro A hA
      rcases List.mem_append.1 hA with hA | hA
      · exact hA
      · exact mem_axUnion.2 (Or.inr (List.mem_append_right _ hA))
    · rw [Dom.join_top] at h; cases h

/-- After `m + 1` iterations from a relation with at most `m` known constraints, the iteration
has reached a fixed point. -/
theorem loopStep_iterate_fixed (τ : Dom n) :
    ∀ (m : ℕ) (σ : Dom n), (∀ s, σ = .rel s → s.cons.known ≤ m) →
      loopStep τ ((loopStep τ)^[m + 1] σ) = (loopStep τ)^[m + 1] σ := by
  intro m
  induction m with
  | zero =>
    intro σ hσ
    cases σ with
    | top => rw [Function.iterate_fixed (loopStep_top τ)]; rfl
    | rel s =>
      rw [Function.iterate_succ_apply, Function.iterate_zero_apply]
      cases h1 : loopStep τ (.rel s) with
      | top => rfl
      | rel s' =>
        obtain ⟨⟨Γ, hΓ⟩, hfix⟩ := loopStep_rel h1
        by_cases hc : s'.cons = s.cons
        · exact hfix hc
        · have := Cons.known_join_lt (Γ' := Γ) (by rw [← hΓ]; exact hc)
          have := hσ s rfl
          omega
  | succ m ih =>
    intro σ hσ
    cases σ with
    | top => rw [Function.iterate_fixed (loopStep_top τ)]; rfl
    | rel s =>
      rw [Function.iterate_succ_apply]
      cases h1 : loopStep τ (.rel s) with
      | top => rw [Function.iterate_fixed (loopStep_top τ)]; rfl
      | rel s' =>
        obtain ⟨⟨Γ, hΓ⟩, hfix⟩ := loopStep_rel h1
        by_cases hc : s'.cons = s.cons
        · rw [Function.iterate_fixed (hfix hc)]; exact hfix hc
        · have hlt := Cons.known_join_lt (Γ' := Γ) (by rw [← hΓ]; exact hc)
          have := hσ s rfl
          exact ih (.rel s') (fun s'' h => by cases h; rw [hΓ]; omega)

/-- **Paper, Lemma 6.12**, for the iterates of any loop body abstraction `τ`:
`σ_i = σ_{2n+1}` for every `i ≥ 2n + 1`. -/
theorem loopIter_stable (τ : Dom n) {i : ℕ} (hi : 2 * n + 1 ≤ i) :
    loopIter τ i = loopIter τ (2 * n + 1) := by
  obtain ⟨j, rfl⟩ := Nat.exists_eq_add_of_le' hi
  rw [loopIter_eq_iterate, loopIter_eq_iterate, Function.iterate_add_apply]
  exact Function.iterate_fixed (loopStep_iterate_fixed τ (2 * n) (alphaI n)
    (fun s h => by cases h; exact Cons.known_le _)) j

/-- **Paper, Lemma 6.12**: in Definition 6.11, `σ_i = σ_{2n+1}` for every `i ≥ 2n + 1`. -/
theorem lemma_6_12 (C : Program n) {i : ℕ} (hi : 2 * n + 1 ≤ i) :
    loopIter C.alpha i = (Program.loop C).alpha :=
  loopIter_stable C.alpha hi

/-- Once two consecutive iterates agree, all later ones agree too. -/
theorem loopIter_eq_of_succ {τ : Dom n} {k : ℕ} (h : loopIter τ (k + 1) = loopIter τ k) :
    ∀ i, k ≤ i → loopIter τ i = loopIter τ k := by
  intro i hi
  induction i, hi using Nat.le_induction with
  | base => rfl
  | succ i _ ih => rw [loopIter_succ, ih, ← loopIter_succ, h]

/-! ## Soundness for programs (Theorem 6.14) -/

theorem loopIter_sat (τ : Dom n) : ∀ i, (loopIter τ i).Sat
  | 0 => ⟨1, models_id⟩
  | i + 1 => Dom.Sat.join_left (loopIter_sat τ i) _

theorem loopIter_le (τ : Dom n) {i j : ℕ} (h : i ≤ j) : (loopIter τ i).Le (loopIter τ j) := by
  induction j, h using Nat.le_induction with
  | base => exact Dom.Le.refl _
  | succ j _ ih => exact ih.trans (Dom.le_join_left _ _)

theorem one_mem_alphaI : (1 : Density n) ∈ (alphaI n).gamma := by
  rw [alphaI, mem_gamma_id_iff]; intro P _; simp

namespace Program

/-- Paper, Section 6.4: every `α(C)` has a model, so it is an abstract relation. -/
theorem alpha_sat : ∀ C : Program n, (alpha C).Sat
  | gate G => ⟨_, G.models⟩
  | rz _ _ => ⟨1, models_id⟩
  | meas _ => ⟨1, models_id⟩
  | reset _ => ⟨1, models_id⟩
  | seq C₁ C₂ => Dom.Sat.comp (alpha_sat C₁) (alpha_sat C₂)
  | choice C₁ _ => Dom.Sat.join_left (alpha_sat C₁) _
  | loop _ => loopIter_sat _ _

/-- The runs of `k` iterations of a loop body are in `γ(σ_k)`. -/
theorem powSem_subset {S : Set (Density n)} {τ : Dom n} (hS : S ⊆ τ.gamma) :
    ∀ k, powSem S k ⊆ (loopIter τ k).gamma
  | 0 => by rintro K rfl; exact one_mem_alphaI
  | k + 1 => by
    rintro M ⟨K, hK, K', hK', rfl⟩
    exact Dom.gamma_mono (Dom.le_join_right _ _)
      (gamma_comp (loopIter_sat τ k) (powSem_subset hS k hK) (hS hK'))

/-- **Paper, Theorem 6.14 (Soundness for programs)**: `𝒪⟦C⟧ ⊆ γ(α(C))` for every program `C`. -/
theorem soundness : ∀ C : Program n, psem C ⊆ (alpha C).gamma
  | gate G => by rintro K rfl; exact Circuit.gate_mem_gamma G
  | rz q θ => by rintro K rfl; exact Circuit.rz_mem_gamma q θ θ
  | meas q => by
    rintro K (rfl | rfl)
    · exact lemma_6_2_measure q 0
    · exact lemma_6_2_measure q 1
  | reset q => by
    rintro K (rfl | rfl)
    · exact lemma_6_2_reset_zero q
    · exact lemma_6_2_reset_one q
  | seq C₁ C₂ => by
    rintro M ⟨K₁, h₁, K₂, h₂, rfl⟩
    exact gamma_comp (alpha_sat C₁) (soundness C₁ h₁) (soundness C₂ h₂)
  | choice C₁ C₂ => fun K hK => lemma_6_7 _ _ (Or.imp (fun h => soundness C₁ h)
      (fun h => soundness C₂ h) hK)
  | loop C => by
    intro K hK
    obtain ⟨k, hk⟩ := Set.mem_iUnion.1 hK
    have h := powSem_subset (soundness C) k hk
    rcases le_total k (2 * n + 1) with hle | hle
    · exact Dom.gamma_mono (loopIter_le _ hle) h
    · rw [alpha_loop, ← loopIter_stable _ hle]; exact h

end Program

/-- **Paper, Theorem 6.14 (Soundness for programs)**: `𝒪⟦C⟧ ⊆ γ(α(C))`. -/
theorem theorem_6_14 (C : Program n) : C.psem ⊆ C.alpha.gamma := Program.soundness C

/-! ## The disjunctive refinement (Section 8.4, RQ2 "Configurations")

The implementation keeps up to `N` separate abstract relations at branches and loops before
joining them (`N = 100` in the evaluation). Here an abstract state is a list of abstract relations
(each from the program's input), its concretization is the union, and the analysis runs forward:

* an atomic operation `C` composes `α(C)` with every relation of the list;
* a conditional runs both branches from the same list, concatenates the results, and joins them
  into a single relation with `⊔` if there are more than `N` (`collapse`);
* a loop iterates `S ↦ collapse N (S ++ new)`, where `new` are the relations of one more
  iteration of the body not already below a relation of `S`, until one more iteration adds nothing
  (a post-fixpoint, checked decidably); after `fuel` rounds without one it gives up with `[⊤]`.
-/

namespace Disj

/-- The bound `N` used in the evaluation (Section 8). -/
def implBound : ℕ := 100

/-- The concretization of a list of abstract relations: the union of their concretizations. -/
def gamma (L : List (Dom n)) : Set (Density n) := {K | ∃ σ ∈ L, K ∈ σ.gamma}

/-- Every relation of the list has a model. -/
def AllSat (L : List (Dom n)) : Prop := ∀ σ ∈ L, σ.Sat

/-- `L ⊑ L'`: every relation of `L` is below one of `L'`. -/
def LeL (L L' : List (Dom n)) : Prop := ∀ σ ∈ L, ∃ τ ∈ L', σ.Le τ

instance (L L' : List (Dom n)) : Decidable (LeL L L') := by unfold LeL; infer_instance

theorem LeL.refl (L : List (Dom n)) : LeL L L := fun σ h => ⟨σ, h, Dom.Le.refl σ⟩

theorem LeL.trans {L₁ L₂ L₃ : List (Dom n)} (h₁ : LeL L₁ L₂) (h₂ : LeL L₂ L₃) : LeL L₁ L₃ := by
  intro σ hσ
  obtain ⟨τ, hτ, hle⟩ := h₁ σ hσ
  obtain ⟨ρ, hρ, hle'⟩ := h₂ τ hτ
  exact ⟨ρ, hρ, hle.trans hle'⟩

theorem LeL.of_subset {L L' : List (Dom n)} (h : ∀ σ ∈ L, σ ∈ L') : LeL L L' :=
  fun σ hσ => ⟨σ, h σ hσ, Dom.Le.refl σ⟩

theorem gamma_mono {L L' : List (Dom n)} (h : LeL L L') : gamma L ⊆ gamma L' := by
  rintro K ⟨σ, hσ, hK⟩
  obtain ⟨τ, hτ, hle⟩ := h σ hσ
  exact ⟨τ, hτ, Dom.gamma_mono hle hK⟩

/-- The join of all the relations of a nonempty list (`⊤` for the empty list, never used). -/
def joinAll : List (Dom n) → Dom n
  | [] => .top
  | σ :: L => L.foldl (· ⊔ ·) σ

theorem le_foldl_join : ∀ (L : List (Dom n)) (acc : Dom n),
    acc.Le (L.foldl (· ⊔ ·) acc) ∧ ∀ σ ∈ L, σ.Le (L.foldl (· ⊔ ·) acc)
  | [], acc => ⟨Dom.Le.refl _, by simp⟩
  | τ :: L, acc => by
    obtain ⟨h1, h2⟩ := le_foldl_join L (acc ⊔ τ)
    refine ⟨(Dom.le_join_left _ _).trans h1, ?_⟩
    intro σ hσ
    rcases List.mem_cons.1 hσ with rfl | hσ
    · exact (Dom.le_join_right _ _).trans h1
    · exact h2 σ hσ

theorem le_joinAll : ∀ (L : List (Dom n)), ∀ σ ∈ L, σ.Le (joinAll L)
  | [], _, h => by cases h
  | τ :: L, σ, hσ => by
    obtain ⟨h1, h2⟩ := le_foldl_join L τ
    rcases List.mem_cons.1 hσ with rfl | hσ
    · exact h1
    · exact h2 σ hσ

/-- Keep up to `N` relations; join them into one when there are more. -/
def collapse (N : ℕ) (L : List (Dom n)) : List (Dom n) :=
  if L.length ≤ N then L else [joinAll L]

theorem le_collapse (N : ℕ) (L : List (Dom n)) : LeL L (collapse N L) := by
  unfold collapse
  split_ifs
  · exact LeL.refl L
  · exact fun σ hσ => ⟨_, List.mem_singleton_self _, le_joinAll L σ hσ⟩

theorem allSat_collapse (N : ℕ) {L : List (Dom n)} (h : AllSat L) : AllSat (collapse N L) := by
  unfold collapse
  split_ifs
  · exact h
  · intro σ hσ
    rw [List.mem_singleton] at hσ
    subst hσ
    cases L with
    | nil => trivial
    | cons τ L => exact Dom.Sat.mono (le_joinAll _ τ (List.mem_cons_self ..)) (h τ (List.mem_cons_self ..))

theorem allSat_append {L L' : List (Dom n)} (h : AllSat L) (h' : AllSat L') : AllSat (L ++ L') :=
  fun σ hσ => (List.mem_append.1 hσ).elim (h σ) (h' σ)

/-- The relations of `T` that are not below a relation of `S`. -/
def fresh (S T : List (Dom n)) : List (Dom n) :=
  T.filter (fun σ => decide (¬ ∃ τ ∈ S, σ.Le τ))

theorem mem_of_mem_fresh {S T : List (Dom n)} {σ : Dom n} (h : σ ∈ fresh S T) : σ ∈ T :=
  (List.mem_filter.1 h).1

/-- Loop iteration to a post-fixpoint: stop when one more iteration of the body `f` adds nothing
new; otherwise add the new relations (collapsing above `N`) and try again; give up with `[⊤]`
when the fuel runs out. -/
def loopFix (N : ℕ) (f : List (Dom n) → List (Dom n)) : ℕ → List (Dom n) → List (Dom n)
  | 0, _ => [.top]
  | k + 1, S =>
    if LeL (f S) S then S
    else loopFix N f k (collapse N (S ++ fresh S (f S)))

theorem loopFix_spec (N : ℕ) (f : List (Dom n) → List (Dom n))
    (hf : ∀ S, AllSat S → AllSat (f S)) :
    ∀ (k : ℕ) (S : List (Dom n)), AllSat S →
      AllSat (loopFix N f k S) ∧ LeL S (loopFix N f k S) ∧
        LeL (f (loopFix N f k S)) (loopFix N f k S)
  | 0, S, _ =>
    ⟨fun σ hσ => by rw [List.mem_singleton.1 hσ]; trivial,
      fun _ _ => ⟨.top, List.mem_singleton_self _, trivial⟩,
      fun _ _ => ⟨.top, List.mem_singleton_self _, trivial⟩⟩
  | k + 1, S, hS => by
    simp only [loopFix]
    split_ifs with hpost
    · exact ⟨hS, LeL.refl S, hpost⟩
    · have hS' : AllSat (collapse N (S ++ fresh S (f S))) :=
        allSat_collapse N (allSat_append hS (fun σ hσ => hf S hS σ (mem_of_mem_fresh hσ)))
      obtain ⟨h1, h2, h3⟩ := loopFix_spec N f hf k _ hS'
      refine ⟨h1, LeL.trans ?_ h2, h3⟩
      exact (LeL.of_subset fun σ hσ => List.mem_append_left _ hσ).trans (le_collapse N _)

/-- The fuel of the loop iteration (enough rounds of up to `N` steps for each of the `2n`
constraints to be dropped; soundness does not depend on it). -/
def fuel (N n : ℕ) : ℕ := (2 * n + 3) * (N + 2)

/-- **The disjunctive analysis** with bound `N`: the abstract state after running `C` from the
list of relations `S`. -/
def run (N : ℕ) : Program n → List (Dom n) → List (Dom n)
  | .gate G, S => S.map ((Program.gate G).alpha ⊚ ·)
  | .rz q θ, S => S.map ((Program.rz q θ).alpha ⊚ ·)
  | .meas q, S => S.map ((Program.meas q).alpha ⊚ ·)
  | .reset q, S => S.map ((Program.reset q).alpha ⊚ ·)
  | .seq C₁ C₂, S => run N C₂ (run N C₁ S)
  | .choice C₁ C₂, S => collapse N (run N C₁ S ++ run N C₂ S)
  | .loop C, S => loopFix N (run N C) (fuel N n) S

/-- An atomic step `S ↦ map (α(C) ∘ ·) S` is sound. -/
theorem map_sound (C : Program n) {S : List (Dom n)} (hS : AllSat S) :
    AllSat (S.map (C.alpha ⊚ ·)) ∧
      ∀ K ∈ C.psem, ∀ M ∈ gamma S, K * M ∈ gamma (S.map (C.alpha ⊚ ·)) := by
  refine ⟨fun σ hσ => ?_, ?_⟩
  · obtain ⟨τ, hτ, rfl⟩ := List.mem_map.1 hσ
    exact Dom.Sat.comp (hS τ hτ) (Program.alpha_sat C)
  · rintro K hK M ⟨σ, hσ, hM⟩
    exact ⟨_, List.mem_map_of_mem hσ, gamma_comp (hS σ hσ) hM (Program.soundness C hK)⟩

/-- **Soundness of the disjunctive analysis**, for programs with measurement and reset: if every
relation of `S` has a model, so does every relation of `run N C S`, and every run `K` of `C`
takes every operator `M` described by `S` to an operator `K M` described by `run N C S`. -/
theorem run_sound (N : ℕ) (C : Program n) :
    ∀ S : List (Dom n), AllSat S →
      AllSat (run N C S) ∧ ∀ K ∈ C.psem, ∀ M ∈ gamma S, K * M ∈ gamma (run N C S) := by
  induction C with
  | gate G => exact fun S hS => map_sound (.gate G) hS
  | rz q θ => exact fun S hS => map_sound (.rz q θ) hS
  | meas q => exact fun S hS => map_sound (.meas q) hS
  | reset q => exact fun S hS => map_sound (.reset q) hS
  | seq C₁ C₂ ih₁ ih₂ =>
    intro S hS
    obtain ⟨h1, h1'⟩ := ih₁ S hS
    obtain ⟨h2, h2'⟩ := ih₂ _ h1
    refine ⟨h2, ?_⟩
    rintro K ⟨K₁, hK₁, K₂, hK₂, rfl⟩ M hM
    rw [Matrix.mul_assoc]
    exact h2' K₂ hK₂ _ (h1' K₁ hK₁ M hM)
  | choice C₁ C₂ ih₁ ih₂ =>
    intro S hS
    obtain ⟨h1, h1'⟩ := ih₁ S hS
    obtain ⟨h2, h2'⟩ := ih₂ S hS
    refine ⟨allSat_collapse N (allSat_append h1 h2), ?_⟩
    rintro K (hK | hK) M hM
    · exact gamma_mono ((LeL.of_subset fun σ hσ => List.mem_append_left _ hσ).trans
        (le_collapse N _)) (h1' K hK M hM)
    · exact gamma_mono ((LeL.of_subset fun σ hσ => List.mem_append_right _ hσ).trans
        (le_collapse N _)) (h2' K hK M hM)
  | loop C ih =>
    intro S hS
    obtain ⟨hR, hSR, hpost⟩ := loopFix_spec N (run N C) (fun S hS => (ih S hS).1) (fuel N n) S hS
    refine ⟨hR, ?_⟩
    set R := loopFix N (run N C) (fuel N n) S
    have key : ∀ k, ∀ K ∈ Program.powSem C.psem k, ∀ M ∈ gamma S, K * M ∈ gamma R := by
      intro k
      induction k with
      | zero =>
        rintro K rfl M hM
        rw [Matrix.one_mul]; exact gamma_mono hSR hM
      | succ k ihk =>
        rintro _ ⟨K, hK, K', hK', rfl⟩ M hM
        rw [Matrix.mul_assoc]
        exact gamma_mono hpost ((ih R hR).2 K' hK' _ (ihk K hK M hM))
    intro K hK M hM
    obtain ⟨k, hk⟩ := Set.mem_iUnion.1 hK
    exact key k K hk M hM

/-- **Soundness of the disjunctive refinement**: every run of a program is described by one of
the relations the disjunctive analysis computes from `[α(I)]`, for every bound `N`. -/
theorem program_sound (N : ℕ) (C : Program n) : C.psem ⊆ gamma (run N C [alphaI n]) := by
  intro K hK
  have h := (run_sound N C [alphaI n] (fun σ hσ => by
    rw [List.mem_singleton.1 hσ]; exact ⟨1, models_id⟩)).2 K hK 1
    ⟨_, List.mem_singleton_self _, one_mem_alphaI⟩
  rwa [Matrix.mul_one] at h

/-- With the evaluation's bound `N = 100`. -/
theorem program_sound_impl (C : Program n) : C.psem ⊆ gamma (run implBound C [alphaI n]) :=
  program_sound implBound C

/-- Every relation the disjunctive analysis computes has a model. -/
theorem run_allSat (N : ℕ) (C : Program n) : AllSat (run N C [alphaI n]) :=
  (run_sound N C [alphaI n] (fun σ hσ => by
    rw [List.mem_singleton.1 hσ]; exact ⟨1, models_id⟩)).1

end Disj

/-! ## Examples (qubits `1, 2` of the paper are `0, 1`) -/

/-- Conjugation by a model preserves commutation of weakest preconditions:
`[wp(P), wp(A)] = 0 ↔ [P, A] = 0`. -/
theorem comm_wp_iff {Γ : Cons n} {V : Density n} (hV : Models V Γ) {P A B A' : PP n}
    (hP : Γ.wp P = some B) (hA : Γ.wp A = some A') : PP.Comm B A' ↔ PP.Comm P A := by
  rw [PP.comm_iff, PP.comm_iff, Cons.wp_toMatrix hV hP, Cons.wp_toMatrix hV hA]
  exact comm_conj_iff hV.unitary

/-- A Pauli string commutes with `X_q` and `Z_q` exactly when it does not touch `q`. -/
theorem comm_X_Z_iff (P : PP n) (q : Fin n) :
    (PP.Comm P (Xq q) ∧ PP.Comm P (Zq q)) ↔ P.letter q = .I := by
  rw [PP.comm_single, PP.comm_single]
  cases P.letter q <;> decide

/-- An operator commuting with `Z_q` commutes with every Z-rotation of `q`. -/
theorem comm_rzMat_of_comm {K : Density n} {q : Fin n}
    (h : K * (Zq q).toMatrix = (Zq q).toMatrix * K) (θ : ℝ) : K * rzMat q θ = rzMat q θ * K := by
  obtain ⟨a, b, hab⟩ := rzMat_combo q θ
  rw [hab, Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.mul_smul, Matrix.smul_mul,
    Matrix.smul_mul, Matrix.mul_one, Matrix.one_mul, h]

namespace Examples

open Program

/-- The tableau of `H₁` on two qubits: `Z₁ ↦ X₁`, `X₁ ↦ Z₁`, `X₂ ↦ X₂`, `Z₂ ↦ Z₂`. -/
theorem consH_eq : (gH (0 : Fin 2)).cons =
    ⟨![some (Zq 0), some (Xq 1)], ![some (Xq 0), some (Zq 1)]⟩ := by decide

/-- The program `H₁; measure 1` of Section 2.4. -/
def hMeas : Program 2 := seq (gate (gH 0)) (meas 0)

/-- The program `H₁; reset 1`. -/
def hReset : Program 2 := seq (gate (gH 0)) (reset 0)

/-- **Paper, Example 6.3**: `α(H₁; measure 1)` keeps the tableau of `H₁` and adds the axis
`wp_{Γ_{H₁}}(Z₁) = X₁`. -/
theorem example_6_3 : hMeas.alpha = .rel ⟨(gH 0).cons, [Xq 0]⟩ ∧
    (gH (0 : Fin 2)).cons.wp (Zq 0) = some (Xq 0) := by decide

/-- Paper, Example 6.3: `Z₁` is tracked (its weakest precondition `X₁` commutes with the axis),
`X₁` is not (its weakest precondition `Z₁` anticommutes with the axis). -/
theorem example_6_3_F : Zq 0 ∈ (⟨(gH 0).cons, [Xq 0]⟩ : ARel 2).F ∧
    (gH (0 : Fin 2)).cons.wp (Xq 0) = some (Zq 0) ∧ ¬ PP.Comm (Zq (0 : Fin 2)) (Xq 0) ∧
    Xq 0 ∉ (⟨(gH 0).cons, [Xq 0]⟩ : ARel 2).F := by decide

/-- Paper, Example 6.3: on both outcomes, the program takes `X₁` to `Z₁`. -/
theorem example_6_3_takes : ∀ K ∈ hMeas.psem, K * (Xq 0).toMatrix = (Zq 0).toMatrix * K := by
  intro K hK
  have h := theorem_6_14 hMeas hK
  rw [example_6_3.1] at h
  exact h (Zq 0) (Xq 0) (by decide) (by decide)

/-- Paper, Example 6.3 (reset): the axes are `{Z₁, X₁}`. -/
theorem example_6_3_reset : hReset.alpha = .rel ⟨(gH 0).cons, [Zq 0, Xq 0]⟩ := by decide

/-- Paper, Example 6.3 (reset): the abstraction tracks exactly the Pauli strings that do not
touch `q₁`. -/
theorem example_6_3_reset_F :
    (⟨(gH 0).cons, [Zq 0, Xq 0]⟩ : ARel 2).F = {P | P.letter 0 = .I} := by
  ext P
  obtain ⟨B, hB, -⟩ := Cons.wp_conj (Cons.ofTables_allKnown _ _) (gH (0 : Fin 2)).models P
  have hZ : (gH (0 : Fin 2)).cons.wp (Zq 0) = some (Xq 0) := by decide
  have hX : (gH (0 : Fin 2)).cons.wp (Xq 0) = some (Zq 0) := by decide
  have hm := (gH (0 : Fin 2)).models
  simp only [ARel.F, ARel.Tracked, Set.mem_ofPred_eq, List.mem_cons, List.not_mem_nil,
    or_false, forall_eq_or_imp, forall_eq]
  constructor
  · rintro ⟨B', hB', h1, h2⟩
    exact (comm_X_Z_iff P 0).1 ⟨(comm_wp_iff hm hB' hX).1 h1, (comm_wp_iff hm hB' hZ).1 h2⟩
  · intro h
    obtain ⟨h1, h2⟩ := (comm_X_Z_iff P 0).2 h
    exact ⟨B, hB, (comm_wp_iff hm hB hX).2 h1, (comm_wp_iff hm hB hZ).2 h2⟩

/-- Paper, Section 4.2: no constraint on `Z₂`, and every other constraint that of the identity. -/
def consUnk : Cons 2 := ⟨![some (Xq 0), some (Xq 1)], ![some (Zq 0), none]⟩

/-- **Paper, Section 4.2**: `wp(X₂) = X₂`, but `wp(Z₁Z₂)` is undefined. -/
theorem unknowns_example :
    consUnk.wp (Xq 1) = some (Xq 1) ∧ consUnk.wp ((Zq 0).mul (Zq 1)) = none := by decide

/-- The program `H₁; measure 1; if (*) {X₂}` of Section 2.4 and Example 6.8. -/
def hMeasIf : Program 2 := seq hMeas (ifThen (gate (gX 1)))

/-- The relation after the conditional (Example 6.8 and Section 2.4):
`Z₁ ↦ X₁, X₁ ↦ Z₁, X₂ ↦ X₂` (no constraint on `Z₂`), `𝒜 = {X₁}`. -/
def sigmaJoin : ARel 2 := ⟨⟨![some (Zq 0), some (Xq 1)], ![some (Xq 0), none]⟩, [Xq 0]⟩

/-- **Paper, Example 6.8**: the `if` branch changes `Z₂ ↦ Z₂` to `-Z₂ ↦ Z₂`, the `else` branch
changes nothing, and their join is `sigmaJoin`; it is the abstraction of the program. -/
theorem example_6_8 :
    (gate (gX (1 : Fin 2))).alpha ⊚ hMeas.alpha =
      .rel ⟨⟨![some (Zq 0), some (Xq 1)], ![some (Xq 0), some (Zq 1).neg]⟩, [Xq 0]⟩ ∧
    (skip 2).alpha ⊚ hMeas.alpha = hMeas.alpha ∧
    ((gate (gX (1 : Fin 2))).alpha ⊚ hMeas.alpha) ⊔ ((skip 2).alpha ⊚ hMeas.alpha) =
      .rel sigmaJoin ∧
    hMeasIf.alpha = .rel sigmaJoin := by decide

/-- Paper, Example 6.8: `X₂` is still tracked, `Z₂` is not. -/
theorem example_6_8_F : Xq 1 ∈ sigmaJoin.F ∧ Zq 1 ∉ sigmaJoin.F := by decide

/-- Paper, Example 6.8: on either branch, the program takes `X₂` to `X₂`. -/
theorem example_6_8_takes : ∀ K ∈ hMeasIf.psem, K * (Xq 1).toMatrix = (Xq 1).toMatrix * K := by
  intro K hK
  have h := theorem_6_14 hMeasIf hK
  rw [example_6_8.2.2.2] at h
  exact h (Xq 1) (Xq 1) (by decide) (by decide)

/-- `CX_{1,2}`. -/
abbrev cx12 : CliffordGate 2 := gCX 0 1 (by decide)

/-- The loop `while (*) {CX_{1,2}}` of Loop-simple. -/
def loopCX : Program 2 := loop (gate cx12)

/-- Loop-simple of Amy and Lunderville, `T₁; while (*) {CX_{1,2}}; T₁†`. -/
def loopSimple : Program 2 := seq (seq (rz 0 (Real.pi / 4)) loopCX) (rz 0 (-(Real.pi / 4)))

/-- `σ₁` of Example 6.13: no constraints on `X₁` and `Z₂`, the other constraints those of the
identity, no axes. -/
def sigmaLoop : Dom 2 := .rel ⟨⟨![none, some (Xq 1)], ![some (Zq 0), none]⟩, []⟩

/-- **Paper, Example 6.13**: `σ₀ = α(I)`; `α(CX_{1,2})` disagrees with the identity on the constraint on `X₁` and
the constraint on `Z₂`, so `σ₁` drops them; and `σ₂ = σ₁`. -/
theorem example_6_13_iterates :
    loopIter (gate cx12).alpha 0 = alphaI 2 ∧
    (gate cx12).alpha = .rel ⟨⟨![some ((Xq 0).mul (Xq 1)), some (Xq 1)],
      ![some (Zq 0), some ((Zq 0).mul (Zq 1))]⟩, []⟩ ∧
    loopIter (gate cx12).alpha 1 = sigmaLoop ∧
    loopIter (gate cx12).alpha 2 = loopIter (gate cx12).alpha 1 := by decide

/-- Paper, Example 6.13: the iteration has stabilized, `σ_i = σ₁` for every `i ≥ 1`, so
`α(while (*) {CX_{1,2}}) = σ₁`. -/
theorem example_6_13_stable : (∀ i, 1 ≤ i → loopIter (gate cx12).alpha i = sigmaLoop) ∧
    loopCX.alpha = sigmaLoop := by
  have h : ∀ i, 1 ≤ i → loopIter (gate cx12).alpha i = sigmaLoop := fun i hi => by
    rw [loopIter_eq_of_succ example_6_13_iterates.2.2.2 i hi, example_6_13_iterates.2.2.1]
  exact ⟨h, h _ (by omega)⟩

/-- Paper, Example 6.13: throughout, `Z₁ ↦ Z₁` and the loop adds no axes. -/
theorem example_6_13_z1 : ∀ i, ∃ s, loopIter (gate cx12).alpha i = .rel s ∧
    s.cons.z 0 = some (Zq 0) ∧ s.axes = [] := by
  intro i
  rcases Nat.eq_zero_or_pos i with rfl | hi
  · exact ⟨_, rfl, rfl, rfl⟩
  · exact ⟨_, example_6_13_stable.1 i hi, rfl, rfl⟩

/-- Paper, Example 6.13: however many times the loop runs, it takes `Z₁` to `Z₁`. -/
theorem example_6_13_takes :
    ∀ K ∈ loopCX.psem, K * (Zq 0).toMatrix = (Zq 0).toMatrix * K := by
  intro K hK
  have h := theorem_6_14 loopCX hK
  rw [example_6_13_stable.2] at h
  exact h (Zq 0) (Zq 0) (by decide) (by simp)

/-- Paper, Example 6.13: the `T₁†` after the loop cancels the `T₁` before it: for every run `K`
of the loop, `T₁† K T₁ = K`. -/
theorem example_6_13_cancel :
    ∀ K ∈ loopCX.psem, rzMat 0 (-(Real.pi / 4)) * K * rzMat 0 (Real.pi / 4) = K := by
  intro K hK
  rw [Matrix.mul_assoc, comm_rzMat_of_comm (example_6_13_takes K hK), ← Matrix.mul_assoc,
    rzMat_mul, neg_add_cancel, rzMat_zero, Matrix.one_mul]

/-- Paper, Example 6.13: so the runs of Loop-simple are exactly the runs of the loop alone. -/
theorem example_6_13_runs : loopSimple.psem = loopCX.psem := by
  ext M
  constructor
  · rintro ⟨_, ⟨_, rfl, K, hK, rfl⟩, _, rfl, rfl⟩
    rw [← Matrix.mul_assoc, example_6_13_cancel K hK]; exact hK
  · intro hM
    exact ⟨M * rzMat 0 (Real.pi / 4), ⟨_, rfl, M, hM, rfl⟩, _, rfl,
      by rw [← Matrix.mul_assoc, example_6_13_cancel M hM]⟩

/-! ### The disjunctive refinement on Loop-simple -/

/-- With the bound `N = 100`, the disjunctive analysis of `while (*) {CX_{1,2}}` keeps the two
relations `α(I)` and `α(CX_{1,2})` apart instead of joining them. -/
theorem disj_loopCX :
    Disj.run Disj.implBound loopCX [alphaI 2] = [alphaI 2, (gate cx12).alpha] := by decide

/-- The single relation `σ₁` does not track `X₁` (it has no constraint on `X₁`). -/
theorem sigmaLoop_X1 : ∀ s, sigmaLoop = .rel s → Xq 0 ∉ s.F := by
  rintro s ⟨⟩; decide

/-- The disjunctive analysis does: every run of the loop takes `X₁` to `X₁` or `X₁X₂` to
`X₁`. -/
theorem disj_loopCX_takes : ∀ K ∈ loopCX.psem,
    K * (Xq 0).toMatrix = (Xq 0).toMatrix * K ∨
      K * ((Xq 0).mul (Xq 1)).toMatrix = (Xq 0).toMatrix * K := by
  intro K hK
  obtain ⟨σ, hσ, hK⟩ := Disj.program_sound_impl loopCX hK
  rw [disj_loopCX] at hσ
  rcases List.mem_pair.1 hσ with rfl | rfl
  · exact Or.inl (hK (Xq 0) (Xq 0) (by decide) (by simp))
  · exact Or.inr (hK (Xq 0) ((Xq 0).mul (Xq 1)) (by decide) (by simp))

end Examples

end

end PauliFold.Paper
