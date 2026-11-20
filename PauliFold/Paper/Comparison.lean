import PauliFold.Paper.Circuits
import PauliFold.Alpha
import Mathlib.Algebra.Field.ZMod
import Mathlib.LinearAlgebra.Dual.Lemmas

/-!
# Paper, Section 7: comparison to other abstractions

* **Clifford+T circuits as paper circuits.** `CTCirc n` is the syntax of circuits over the
  gates of `CTGate` (the paper's `H, S, X, Z, CX, T, T†`, plus `S†` and `CZ`), with sequential
  composition. `CTCirc.toPaper` embeds it into the paper's `Circuit` (Cliffords ↦ the standard
  Clifford gates of `Circuits.lean`, `T ↦ R_Z(π/4)`, `T† ↦ R_Z(-π/4)`), and `CTCirc.toList` into
  the gate lists of the earlier development. The semantics agree up to a global phase with
  `circ` (`CTCirc.sem_toPaper`), and the paper's abstraction has exactly the rows and the
  (reversed) axes of the earlier `Tab.run` (`CTCirc.alpha_rows_eq_tab`,
  `CTCirc.alpha_axes_eq_tab`).
* **The linear abstraction (Section 7.1).** Parities `Par(V)` (running text), linear
  abstract states and their formula `Φ_ρ` (Definition 7.1), the transformers (running text),
  support (Definition 7.2), and the concretization `γ_L` (Definition 7.3). `γ_L(ρ)` is exactly the
  earlier `affGamma`, the support-based concretization by the relations without Hadamard
  variables (`gammaL_eq_affGamma`); the proof is linear-algebra duality over `𝔽₂`.
* **Theorem 7.4** (Pauli refines linear: `pauli_refines_linear`,
  `pauli_refines_linear_strict`), the linear soundness claim (`linear_soundness`), and the
  claims after Theorem 7.4: a parity relation holds on the support of `K` exactly when `K`
  takes `(-1)^c Z_T` to `Z_S` (`parity_rel_iff_takes`), the entailed relations are the derived
  ones (`LinState.entails_iff_rel`), and the Pauli abstraction tracks `Z_S` with the same
  weakest precondition (`CTCirc.tracks_Zset_of_rel`).

Numbering follows the current `main.aux`: Definitions 7.1–7.3 and Theorem 7.4.
* **Section 7.2**: the two separating circuits against Amy and Lunderville's polynomial
  abstraction `alphaGammaC d` (`polynomial_more_precise`, `pauli_more_precise`) and the
  incomparability theorem `pauli_polynomial_incomparable`.

Qubits are 0-indexed: the paper's qubits 1 and 2 are the Lean qubits `0` and `1`.
-/

namespace PauliFold.Paper

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

open Circuit

/-! ## Clifford+T circuits as paper circuits -/

/-- The paper circuit of a single `CTGate`: Cliffords are the standard Clifford gates of
Section 3.2, `T_q = R_Z^q(π/4)` and `T†_q = R_Z^q(-π/4)`. -/
def embedG : CTGate n → Circuit n
  | .h q => .gate (gH q)
  | .x q => .gate (gX q)
  | .z q => .gate (gZ q)
  | .s q => .gate (gS q)
  | .sdg q => .gate (gSdg q)
  | .t q => Circuit.T q
  | .tdg q => Circuit.Tdg q
  | .cx c t hne => .gate (gCX c t hne)
  | .cz c t hne => .gate (gCZ c t hne)

/-- Circuits over the gates `H, S, S†, X, Z, CX, CZ, T, T†`, with sequential composition. -/
inductive CTCirc (n : ℕ) where
  | gate (g : CTGate n)
  /-- `C₁ ; C₂`. -/
  | seq (C₁ C₂ : CTCirc n)

namespace CTCirc

/-- The gate list of the earlier development (first gate first). -/
def toList : CTCirc n → List (CTGate n)
  | gate g => [g]
  | seq C₁ C₂ => toList C₁ ++ toList C₂

/-- The paper circuit. -/
def toPaper : CTCirc n → Circuit n
  | gate g => embedG g
  | seq C₁ C₂ => toPaper C₁ ⨾ toPaper C₂

@[simp] theorem toList_gate (g : CTGate n) : (gate g).toList = [g] := rfl
@[simp] theorem toList_seq (C₁ C₂ : CTCirc n) : (seq C₁ C₂).toList = C₁.toList ++ C₂.toList := rfl
@[simp] theorem toPaper_gate (g : CTGate n) : (gate g).toPaper = embedG g := rfl
@[simp] theorem toPaper_seq (C₁ C₂ : CTCirc n) : (seq C₁ C₂).toPaper = C₁.toPaper ⨾ C₂.toPaper :=
  rfl

end CTCirc

/-- A Clifford `CTGate` embeds as a standard Clifford gate whose matrix is tzap's. -/
theorem embedG_clifford {g : CTGate n} (hg : g.IsClifford) :
    ∃ G : CliffordGate n, embedG g = .gate G ∧ G.mat = gateUnitary n g.toGate := by
  cases g with
  | t q => exact hg.elim
  | tdg q => exact hg.elim
  | h q => exact ⟨_, rfl, rfl⟩
  | x q => exact ⟨_, rfl, rfl⟩
  | z q => exact ⟨_, rfl, rfl⟩
  | s q => exact ⟨_, rfl, rfl⟩
  | sdg q => exact ⟨_, rfl, rfl⟩
  | cx c t hne => exact ⟨_, rfl, rfl⟩
  | cz c t hne => exact ⟨_, rfl, rfl⟩

/-- A `T` or `T†` embeds as a Z-rotation. -/
theorem embedG_rot {g : CTGate n} (hg : ¬ g.IsClifford) :
    ∃ q θ, embedG g = .rz q θ ∧ (g = .t q ∨ g = .tdg q) := by
  cases g with
  | t q => exact ⟨q, _, rfl, Or.inl rfl⟩
  | tdg q => exact ⟨q, _, rfl, Or.inr rfl⟩
  | _ => exact (hg trivial).elim

/-- The embedding of a gate agrees with tzap's matrix up to a global phase. -/
theorem embedG_sem (g : CTGate n) : EquivP (gateUnitary n g.toGate) (embedG g).sem := by
  cases g with
  | t q => exact ⟨_, gateUnitary_t q⟩
  | tdg q => exact ⟨_, gateUnitary_tdg q⟩
  | _ => exact EquivP.refl _

theorem EquivP.mul {A A' B B' : Density n} (hA : EquivP A A') (hB : EquivP B B') :
    EquivP (A * B) (A' * B') := by
  obtain ⟨φ, rfl⟩ := hA
  obtain ⟨ψ, rfl⟩ := hB
  refine ⟨φ + ψ, ?_⟩
  rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, ← Complex.exp_add]
  push_cast; ring_nf

/-- **The embedding preserves the semantics up to a global phase**: `⟦C⟧ ≡ circ C`. -/
theorem CTCirc.sem_toPaper : ∀ C : CTCirc n, EquivP (circ C.toList) C.toPaper.sem
  | .gate g => by
    simpa [circ] using embedG_sem g
  | .seq C₁ C₂ => by
    simp only [toList_seq, toPaper_seq, Circuit.sem_seq, circ, List.map_append, unitary_append]
    exact EquivP.mul (sem_toPaper C₂) (sem_toPaper C₁)

/-! ## The Pauli abstraction of an embedded circuit tracks the forward Pauli runs -/

theorem pauliRun_top (gs : List (CTGate n)) : PauliAbs.run gs .top = .top := by
  induction gs with
  | nil => rfl
  | cons g gs ih => exact ih

theorem pauliRun_append (l₁ l₂ : List (CTGate n)) (p : PauliAbs n) :
    PauliAbs.run (l₁ ++ l₂) p = PauliAbs.run l₂ (PauliAbs.run l₁ p) := by
  simp [PauliAbs.run, List.foldl_append]

/-- **Every forward run of the per-string Pauli transformer is a tracked string of the paper's
abstraction.** If the transformer of `Pauli.lean` takes `B` through `C` to `P` (never `⊤`),
then `wp(P) = B` in `α(C)`, and `B` commutes with every axis of `α(C)`. -/
theorem CTCirc.tracked_of_pauliRun : ∀ (C : CTCirc n) (B P : SPauli n),
    PauliAbs.run C.toList (.one B) = .one P →
    (alphaRel C.toPaper).cons.wp (PP.ofS P) = some (PP.ofS B) ∧
      ∀ A ∈ (alphaRel C.toPaper).axes, PP.Comm (PP.ofS B) A
  | .gate g, B, P, h => by
    by_cases hg : g.IsClifford
    · obtain ⟨G, hG, hGm⟩ := embedG_clifford hg
      have hstep : PauliAbs.step g (.one B) = .one (B.conjClifford g) := by
        cases g <;> simp_all [PauliAbs.step, CTGate.IsClifford]
      have hP : P = B.conjClifford g := by
        have : PauliAbs.run [g] (.one B) = PauliAbs.step g (.one B) := rfl
        rw [toList_gate, this, hstep] at h
        cases h; rfl
      subst hP
      have hs := step_sound g (.one B) B.toMatrix (Set.mem_singleton _)
      rw [hstep] at hs
      simp only [PauliAbs.gamma, Set.mem_singleton_iff, CTGate.concrete, conj] at hs
      rw [toPaper_gate, hG, alphaRel_gate]
      refine ⟨?_, fun A hA => absurd hA List.not_mem_nil⟩
      obtain ⟨B', hB', hm⟩ := Cons.wp_conj (Cons.ofTables_allKnown _ _) G.models
        (PP.ofS (B.conjClifford g))
      refine hB'.trans (congrArg some (PP.toMatrix_injective ?_))
      rw [hm, PP.toMatrix_ofS, PP.toMatrix_ofS, ← hs, hGm]
      have hGG := gate_adj_mul g
      simp only [← Matrix.mul_assoc]
      rw [hGG, Matrix.one_mul, Matrix.mul_assoc, hGG, Matrix.mul_one]
    · obtain ⟨q, θ, hE, hq⟩ := embedG_rot hg
      have hstep : PauliAbs.step g (.one B) =
          if (B.letter q).anticommutesZ then .top else .one B := by
        rcases hq with rfl | rfl <;> rfl
      have : PauliAbs.run [g] (.one B) = PauliAbs.step g (.one B) := rfl
      rw [toList_gate, this, hstep] at h
      split_ifs at h with hanti
      cases h
      rw [toPaper_gate, hE, alphaRel_rz, Cons.wp_id]
      refine ⟨rfl, fun A hA => ?_⟩
      rw [List.mem_singleton] at hA
      subst hA
      rw [PP.comm_iff, PP.toMatrix_ofS]
      exact comm_Z_of_letter q B (letter_of_not_anti _ hanti)
  | .seq C₁ C₂, B, P, h => by
    rw [toList_seq, pauliRun_append] at h
    cases h₁ : PauliAbs.run C₁.toList (.one B) with
    | top => rw [h₁, pauliRun_top] at h; cases h
    | one B₁ =>
      rw [h₁] at h
      obtain ⟨hw₂, hc₂⟩ := tracked_of_pauliRun C₂ B₁ P h
      obtain ⟨hw₁, hc₁⟩ := tracked_of_pauliRun C₁ B B₁ h₁
      rw [toPaper_seq, Circuit.wp_alphaRel_seq, hw₂, Option.bind_some, hw₁, Circuit.alphaRel_seq]
      refine ⟨rfl, fun A hA => ?_⟩
      rcases List.mem_append.1 hA with hA | hA
      · exact hc₁ A hA
      · obtain ⟨A₂, hA₂, hwA⟩ := List.mem_filterMap.1 hA
        have hm := Circuit.models_alphaRel C₁.toPaper
        have e1 := Cons.wp_toMatrix hm hw₁
        have e2 := Cons.wp_toMatrix hm hwA
        rw [PP.comm_iff, e1, e2, comm_conj_iff hm.unitary, ← PP.comm_iff]
        exact hc₂ A₂ hA₂

/-! ## The paper's abstraction is the tableau of `Tableau.lean` -/

/-- The Clifford part of a gate: its matrix, or `1` for a rotation. -/
def cliffG (g : CTGate n) : Density n :=
  match g with
  | .t _ | .tdg _ => 1
  | g => gateUnitary n g.toGate

/-- The Clifford part of a gate list. -/
def cliffL : List (CTGate n) → Density n
  | [] => 1
  | g :: gs => cliffL gs * cliffG g

theorem cliffL_append (l₁ l₂ : List (CTGate n)) : cliffL (l₁ ++ l₂) = cliffL l₂ * cliffL l₁ := by
  induction l₁ with
  | nil => simp [cliffL]
  | cons g l ih => simp [cliffL, ih, Matrix.mul_assoc]

theorem cliffG_unitary (g : CTGate n) : cliffG g * (cliffG g)ᴴ = 1 := by
  cases g with
  | t _ | tdg _ => simp [cliffG]
  | _ => exact gate_mul_adj _

theorem cliffL_unitary : ∀ l : List (CTGate n), cliffL l * (cliffL l)ᴴ = 1
  | [] => by simp [cliffL]
  | g :: l => by
    rw [cliffL, conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc (cliffG g),
      cliffG_unitary, Matrix.one_mul, cliffL_unitary l]

/-- The Clifford part of the embedded circuit is the Clifford part of the gate list. -/
theorem CTCirc.cliff_toPaper : ∀ C : CTCirc n, cliff C.toPaper = cliffL C.toList
  | .gate g => by
    cases g <;> simp [embedG, cliff, cliffL, cliffG] <;> rfl
  | .seq C₁ C₂ => by
    rw [toPaper_seq, cliff, toList_seq, cliffL_append, cliff_toPaper C₁, cliff_toPaper C₂]

theorem Tab.init_rowsOK : (Tab.init n).RowsOK 1 :=
  fun _ => ⟨by simp [Tab.init], by simp [Tab.init]⟩

/-- The rows of `Tab.run` represent the Clifford part of the gate list. -/
theorem Tab.rowsOK_run : ∀ (l : List (CTGate n)) (T : Tab n) (C₀ : Density n),
    C₀ * C₀ᴴ = 1 → T.RowsOK C₀ → (T.run l).RowsOK (cliffL l * C₀)
  | [], T, C₀, _, h => by simpa [cliffL, Tab.run] using h
  | g :: l, T, C₀, hCC, hrows => by
    show ((T.step g).run l).RowsOK _
    by_cases hg : g.IsClifford
    · have hG : cliffG g = gateUnitary n g.toGate := by cases g <;> first | rfl | exact hg.elim
      have hu : (gateUnitary n g.toGate * C₀) * (gateUnitary n g.toGate * C₀)ᴴ = 1 := by
        rw [conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc C₀, hCC, Matrix.one_mul,
          gate_mul_adj]
      have := Tab.rowsOK_run l _ _ hu (Tab.rowsOK_step hCC hrows hg)
      rwa [cliffL, hG, Matrix.mul_assoc]
    · obtain ⟨q, rfl | rfl⟩ : ∃ q, g = .t q ∨ g = .tdg q := by
        cases g <;> simp_all [CTGate.IsClifford]
      · simpa [cliffL, cliffG] using Tab.rowsOK_run l (T.step (.t q)) C₀ hCC hrows
      · simpa [cliffL, cliffG] using Tab.rowsOK_run l (T.step (.tdg q)) C₀ hCC hrows

/-- **The rows of `α(C)` are the rows of `Tab.run`**: the paper's constraints on `X_q` and `Z_q` for
the embedded circuit are exactly the rows the earlier tableau abstraction computes. -/
theorem CTCirc.alpha_rows_eq_tab (C : CTCirc n) (q : Fin n) :
    (alphaRel C.toPaper).cons.x q = some (((Tab.init n).run C.toList).rx q) ∧
      (alphaRel C.toPaper).cons.z q = some (((Tab.init n).run C.toList).rz q) := by
  have hm := models_alphaRel C.toPaper
  rw [cliff_toPaper] at hm
  have hr := Tab.rowsOK_run C.toList (Tab.init n) 1 (by simp) Tab.init_rowsOK q
  rw [Matrix.mul_one] at hr
  have hk := alphaRel_allKnown C.toPaper q
  obtain ⟨P, hP⟩ := Option.isSome_iff_exists.1 hk.1
  obtain ⟨Q, hQ⟩ := Option.isSome_iff_exists.1 hk.2
  rw [hP, hQ]
  exact ⟨congrArg some (PP.toMatrix_injective (by rw [hm.x q P hP, hr.1])),
    congrArg some (PP.toMatrix_injective (by rw [hm.z q Q hQ, hr.2]))⟩

theorem forall₂_filterMap {f : PP n → Option (PP n)} {V : Density n}
    (hf : ∀ A, ∃ B, f A = some B ∧ B.toMatrix = Vᴴ * A.toMatrix * V) :
    ∀ l : List (PP n),
      List.Forall₂ (fun A B => B.toMatrix = Vᴴ * A.toMatrix * V) l (l.filterMap f)
  | [] => .nil
  | A :: l => by
    obtain ⟨B, hB, hm⟩ := hf A
    rw [List.filterMap_cons_some hB]
    exact .cons hm (forall₂_filterMap hf l)

theorem forall₂_compose {α β γ : Type*} {W : α → β → Prop} {R₂ : α → γ → Prop}
    {R : β → γ → Prop} (h : ∀ a b c, W a b → R₂ a c → R b c) :
    ∀ {l : List α} {l' : List β} {L : List γ}, List.Forall₂ W l l' → List.Forall₂ R₂ l L →
      List.Forall₂ R l' L
  | _, _, _, .nil, .nil => .nil
  | _, _, _, .cons hw tw, .cons hr tr => .cons (h _ _ _ hw hr) (forall₂_compose h tw tr)

theorem Tab.run_append (l₁ l₂ : List (CTGate n)) (T : Tab n) :
    T.run (l₁ ++ l₂) = (T.run l₁).run l₂ := by
  simp [Tab.run, List.foldl_append]

theorem CTCirc.tab_axes_run : ∀ (C : CTCirc n) (T : Tab n) (C₀ : Density n),
    C₀ * C₀ᴴ = 1 → T.RowsOK C₀ →
    ∃ L, (T.run C.toList).axes = L ++ T.axes ∧
      List.Forall₂ (fun A B => B.toMatrix = C₀ᴴ * A.toMatrix * C₀) (alphaRel C.toPaper).axes.reverse L
  | .gate g, T, C₀, _, hrows => by
    by_cases hg : g.IsClifford
    · obtain ⟨G, hG, -⟩ := embedG_clifford hg
      refine ⟨[], ?_, ?_⟩
      · exact Tab.axes_step_clifford hg
      · rw [toPaper_gate, hG, alphaRel_gate]; exact .nil
    · obtain ⟨q, θ, hE, hq⟩ := embedG_rot hg
      refine ⟨[T.rz q], ?_, ?_⟩
      · show (T.step g).axes = _
        rcases hq with rfl | rfl <;> rfl
      · rw [toPaper_gate, hE, alphaRel_rz]
        exact .cons (hrows q).2 .nil
  | .seq C₁ C₂, T, C₀, hCC, hrows => by
    obtain ⟨L₁, h₁, f₁⟩ := tab_axes_run C₁ T C₀ hCC hrows
    set V₁ := cliffL C₁.toList
    have hu : (V₁ * C₀) * (V₁ * C₀)ᴴ = 1 := by
      rw [conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc C₀, hCC, Matrix.one_mul,
        cliffL_unitary]
    obtain ⟨L₂, h₂, f₂⟩ := tab_axes_run C₂ (T.run C₁.toList) (V₁ * C₀) hu
      (Tab.rowsOK_run _ T C₀ hCC hrows)
    refine ⟨L₂ ++ L₁, ?_, ?_⟩
    · rw [toList_seq, Tab.run_append, h₂, h₁, List.append_assoc]
    · rw [toPaper_seq, alphaRel_seq, List.reverse_append]
      refine List.rel_append ?_ f₁
      have hm := models_alphaRel C₁.toPaper
      rw [cliff_toPaper] at hm
      have hw := forall₂_filterMap (f := (alphaRel C₁.toPaper).cons.wp) (V := V₁)
        (fun A => Cons.wp_conj (alphaRel_allKnown _) hm A) (alphaRel C₂.toPaper).axes
      rw [← List.forall₂_reverse_iff] at hw
      refine forall₂_compose (fun A A' B hW hR => ?_) hw f₂
      rw [hR, hW, conjTranspose_mul]
      simp only [Matrix.mul_assoc]

/-- **The axes of `α(C)` are the axes of `Tab.run`, in reverse order** (the tableau conses each
new axis onto its list; the paper appends it). -/
theorem CTCirc.alpha_axes_eq_tab (C : CTCirc n) :
    ((Tab.init n).run C.toList).axes = (alphaRel C.toPaper).axes.reverse := by
  obtain ⟨L, hL, hf⟩ := C.tab_axes_run (Tab.init n) 1 (by simp) Tab.init_rowsOK
  rw [hL, show (Tab.init n).axes = [] from rfl, List.append_nil]
  have : List.Forall₂ (· = ·) (alphaRel C.toPaper).axes.reverse L :=
    hf.imp fun A B h => (PP.toMatrix_injective (by rw [h]; simp)).symm
  rw [List.forall₂_eq_eq_eq] at this
  exact this.symm

/-! ## The linear abstraction (Section 7.1) -/

/-- Paper, Section 7.1 (Parities): a parity `c ⊕ ⊕_{v ∈ S} v` over the variables
`V = {v₁, …, vₙ} ⊎ V_H`. Variables are natural numbers: `q < n` is the input variable of qubit
`q`, and `n, n + 1, …` are the Hadamard variables. A parity is the pair of `S` (as a finitely
supported `𝔽₂`-valued function) and `c`; this is the type `Affine` of `Phase.lean`. -/
abbrev Parity := Affine

/-- The value of a parity under an assignment `ν` of the variables. -/
def Parity.eval (ν : ℕ → ZMod 2) (p : Parity) : ZMod 2 :=
  p.2 + Finsupp.linearCombination (ZMod 2) ν p.1

theorem Parity.eval_add (ν : ℕ → ZMod 2) (p p' : Parity) :
    Parity.eval ν (p + p') = Parity.eval ν p + Parity.eval ν p' := by
  simp only [Parity.eval, Prod.fst_add, Prod.snd_add, map_add]; ring

theorem Parity.eval_sum {ι : Type*} (ν : ℕ → ZMod 2) (s : Finset ι) (f : ι → Parity) :
    Parity.eval ν (∑ i ∈ s, f i) = ∑ i ∈ s, Parity.eval ν (f i) := by
  simp only [Parity.eval, Prod.fst_sum, Prod.snd_sum, map_sum, Finset.sum_add_distrib]

theorem lc_add (ν μ : ℕ → ZMod 2) (ℓ : Lin) :
    Finsupp.linearCombination (ZMod 2) (ν + μ) ℓ =
      Finsupp.linearCombination (ZMod 2) ν ℓ + Finsupp.linearCombination (ZMod 2) μ ℓ := by
  simp only [Finsupp.linearCombination_apply, Pi.add_apply, smul_add, Finsupp.sum_add]

theorem lc_smul (c : ZMod 2) (μ : ℕ → ZMod 2) (ℓ : Lin) :
    Finsupp.linearCombination (ZMod 2) (c • μ) ℓ = c • Finsupp.linearCombination (ZMod 2) μ ℓ := by
  simp only [Finsupp.linearCombination_apply, Pi.smul_apply, Finsupp.smul_sum]
  refine Finset.sum_congr rfl fun v _ => ?_
  simp only [smul_eq_mul]; ring

theorem lc_indicator (ℓ : Lin) (v : ℕ) :
    Finsupp.linearCombination (ZMod 2) (fun x => if x = v then 1 else 0) ℓ = ℓ v := by
  rw [Finsupp.linearCombination_apply, Finsupp.sum]
  simp only [smul_eq_mul, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finsupp.mem_support_iff]
  split_ifs with h
  · rfl
  · push Not at h; exact h.symm

/-- The value `0` or `1` of a bit. -/
def bit (b : Bool) : ZMod 2 := if b then 1 else 0

/-- Paper, Definition 7.1 (Linear abstract states, `def:linear-state`; the paper's numbering in the
current `main.aux`): the parity `ρ(q)` that each qubit
currently carries. The state also carries the counter of the next fresh Hadamard variable;
it is `AffAbs` of `Affine.lean`. -/
abbrev LinState (n : ℕ) := AffAbs n

/-- Paper, Definition 7.1: the formula `Φ_ρ = ⋀_q (y_q = ρ(q))`, as a predicate of an
assignment `ν` of `V` and of the values `y` of the output variables. -/
def LinState.Phi (ρ : LinState n) (ν : ℕ → ZMod 2) (y : Fin n → ZMod 2) : Prop :=
  ∀ q, y q = Parity.eval ν (ρ.wires q)

/-- Section 7.1: the initial state `ρ₀(q) = v_q`, with the first fresh variable `n`. -/
theorem LinState.init_wires (q : Fin n) :
    (AffAbs.init n).wires q = (Finsupp.single q.val 1, 0) ∧ (AffAbs.init n).fresh = n :=
  ⟨rfl, rfl⟩

/-- Section 7.1, the transformer of `X_q`: `ρ'(q) = ρ(q) ⊕ 1`. -/
theorem LinState.step_x (ρ : LinState n) (q : Fin n) :
    (AffAbs.step (.x q) ρ).wires = Function.update ρ.wires q (ρ.wires q + (0, 1)) := by
  funext r
  by_cases h : r = q
  · subst h; simp [AffAbs.step]
  · simp [AffAbs.step, h]

/-- Section 7.1, the transformer of `CX_{c,t}`: `ρ'(t) = ρ(t) ⊕ ρ(c)`. -/
theorem LinState.step_cx (ρ : LinState n) (c t : Fin n) (hne : c ≠ t) :
    (AffAbs.step (.cx c t hne) ρ).wires = Function.update ρ.wires t (ρ.wires t + ρ.wires c) := by
  funext r
  by_cases h : r = t
  · subst h; simp [AffAbs.step]
  · simp [AffAbs.step, h]

/-- Section 7.1, the transformer of `H_q`: `ρ'(q) = u` for the fresh variable `u`, which is a
Hadamard variable (`n ≤ u`) when the state is reachable (`LinState.fresh_ge`). -/
theorem LinState.step_h (ρ : LinState n) (q : Fin n) :
    (AffAbs.step (.h q) ρ).wires = Function.update ρ.wires q (Finsupp.single ρ.fresh 1, 0) ∧
      (AffAbs.step (.h q) ρ).fresh = ρ.fresh + 1 := by
  refine ⟨?_, rfl⟩
  funext r
  by_cases h : r = q
  · subst h; simp [AffAbs.step]
  · simp [AffAbs.step, h]

/-- Section 7.1, the transformers of the diagonal gates `S_q, Z_q, T_q, T†_q` (and `S†_q`,
`CZ`): `ρ' = ρ`. -/
theorem LinState.step_diag (ρ : LinState n) (g : CTGate n)
    (hg : (∃ q, g = .s q ∨ g = .z q ∨ g = .t q ∨ g = .tdg q ∨ g = .sdg q) ∨
      ∃ c t hne, g = .cz c t hne) :
    AffAbs.step g ρ = ρ := by
  rcases hg with ⟨q, rfl | rfl | rfl | rfl | rfl⟩ | ⟨c, t, hne, rfl⟩ <;> rfl

/-- The fresh variables are Hadamard variables: on every reachable state, `n ≤ fresh`, and no
wire mentions a variable `≥ fresh`. -/
theorem LinState.fresh_ge (gs : List (CTGate n)) :
    n ≤ ((AffAbs.init n).run gs).fresh ∧
      ∀ q v, ((AffAbs.init n).run gs).fresh ≤ v → (((AffAbs.init n).run gs).wires q).1 v = 0 := by
  suffices h : ∀ (gs : List (CTGate n)) (ρ : LinState n), n ≤ ρ.fresh →
      (∀ q v, ρ.fresh ≤ v → (ρ.wires q).1 v = 0) →
      n ≤ (ρ.run gs).fresh ∧ ∀ q v, (ρ.run gs).fresh ≤ v → ((ρ.run gs).wires q).1 v = 0 from
    h gs _ le_rfl fun q v hv => by
      show (Finsupp.single q.val (1 : ZMod 2)) v = 0
      rw [Finsupp.single_apply, if_neg]; have := q.isLt; change n ≤ v at hv; omega
  intro gs
  induction gs with
  | nil => intro ρ h1 h2; exact ⟨h1, h2⟩
  | cons g gs ih =>
    intro ρ h1 h2
    refine ih (AffAbs.step g ρ) ?_ ?_
    · cases g <;> simp [AffAbs.step] <;> omega
    · intro q v hv
      cases g with
      | x r =>
        have := h2 q v hv
        by_cases h : q = r
        · subst h; simp_all [AffAbs.step]
        · simp_all [AffAbs.step]
      | cx c t hne =>
        have := h2 q v hv
        have := h2 c v hv
        by_cases h : q = t
        · subst h; simp_all [AffAbs.step]
        · simp_all [AffAbs.step]
      | h r =>
        change ρ.fresh + 1 ≤ v at hv
        have := h2 q v (by omega)
        by_cases h : q = r
        · subst h; simp [AffAbs.step, Finsupp.single_apply]; omega
        · simp_all [AffAbs.step]
      | _ => exact h2 q v hv

/-- Section 7.1: `ρ_C`, the linear abstract state of a circuit, from `ρ₀` gate by gate. -/
def CTCirc.rho (C : CTCirc n) : LinState n := (AffAbs.init n).run C.toList

theorem AffAbs.run_append (l₁ l₂ : List (CTGate n)) (a : AffAbs n) :
    a.run (l₁ ++ l₂) = (a.run l₁).run l₂ := by
  simp [AffAbs.run, List.foldl_append]

/-- Paper, Definition 7.2 (Support): `supp(K) = {(s, t) | ⟨t|K|s⟩ ≠ 0}`. -/
def supp (K : Density n) : Set (Basis n × Basis n) := {st | K st.2 st.1 ≠ 0}

/-- The assignment of `V` with the inputs `v_q := s_q` and the Hadamard variables given by
`u`. -/
def assign (s : Basis n) (u : ℕ → ZMod 2) : ℕ → ZMod 2 :=
  fun v => if h : v < n then bit (s ⟨v, h⟩) else u v

/-- `Φ_ρ[v_q := s_q, y_q := t_q]_q` is satisfiable: some value of the Hadamard variables
satisfies it. -/
def LinState.SatAt (ρ : LinState n) (s t : Basis n) : Prop :=
  ∃ u : ℕ → ZMod 2, ρ.Phi (assign s u) (fun q => bit (t q))

/-- Paper, Definition 7.3 (Concretization of linear abstract states):
`γ_L(ρ) = {K | Φ_ρ[v_q := s_q, y_q := t_q]_q is satisfiable for every (s, t) ∈ supp(K)}`. -/
def gammaL (ρ : LinState n) : Set (Density n) := {K | ∀ st ∈ supp K, ρ.SatAt st.1 st.2}

/-! ### `γ_L` is the concretization by the relations without Hadamard variables -/

theorem parity_eq_sum_bit (S : Finset (Fin n)) (x : Basis n) :
    parity S x = ∑ i ∈ S, bit (x i) := rfl

theorem assign_lt (s : Basis n) (u : ℕ → ZMod 2) (i : Fin n) : assign s u i.val = bit (s i) := by
  simp [assign, i.isLt]

theorem eval_inputs (s : Basis n) (u : ℕ → ZMod 2) (T : Finset (Fin n)) (c : ZMod 2) :
    Parity.eval (assign s u) (AffAbs.inputs T, c) = c + parity T s := by
  simp only [Parity.eval, AffAbs.inputs, map_sum, Finsupp.linearCombination_single, one_smul,
    assign_lt, parity_eq_sum_bit]

theorem zmod2_sum_mul (s : Finset (Fin n)) (x a : Fin n → ZMod 2) :
    ∑ i ∈ s, x i * a i = ∑ i ∈ s.filter (fun i => a i = 1), x i := by
  rw [Finset.sum_filter]
  refine Finset.sum_congr rfl fun i _ => ?_
  rcases zmod2_cases (a i) with h | h <;> simp [h]

/-- **`Φ_ρ`, with the Hadamard variables projected away, is the conjunction of the relations
without Hadamard variables it entails.** The pair `(s, t)` satisfies `Φ_ρ` for some value of
the Hadamard variables exactly when it satisfies every relation
`⊕_{q ∈ S} y_q = c ⊕ ⊕_{q ∈ T} v_q` with `∑_{q ∈ S} ρ(q) = c ⊕ ⊕_{q ∈ T} v_q` (the relations
`AffAbs.Rel` of `Affine.lean`). This is linear-algebra duality over `𝔽₂`. -/
theorem LinState.satAt_iff_rel (ρ : LinState n) (s t : Basis n) :
    ρ.SatAt s t ↔ ∀ r ∈ ρ.Rel, parity r.1 t = parity r.2.1 s + r.2.2 := by
  constructor
  · rintro ⟨u, hu⟩ ⟨S, T, c⟩ hr
    have hr' : ∑ i ∈ S, ρ.wires i = (AffAbs.inputs T, c) := hr
    simp only
    rw [parity_eq_sum_bit, Finset.sum_congr rfl fun i _ => hu i, ← Parity.eval_sum, hr',
      eval_inputs, add_comm]
  · intro hrel
    -- The Hadamard part of the wires, as a linear map of the assignment.
    let E : (ℕ → ZMod 2) →ₗ[ZMod 2] (Fin n → ZMod 2) :=
      { toFun := fun μ q => Finsupp.linearCombination (ZMod 2) μ (ρ.wires q).1
        map_add' := fun μ μ' => by funext q; exact lc_add μ μ' _
        map_smul' := fun c μ => by funext q; exact lc_smul c μ _ }
    -- The assignments of the Hadamard variables alone.
    let H : Submodule (ZMod 2) (ℕ → ZMod 2) :=
      { carrier := {μ | ∀ v < n, μ v = 0}
        add_mem' := fun ha hb v hv => by simp [ha v hv, hb v hv]
        zero_mem' := fun v _ => rfl
        smul_mem' := fun c μ hμ v hv => by simp [hμ v hv] }
    let β := assign s 0
    let w : Fin n → ZMod 2 := fun q => bit (t q) - Parity.eval β (ρ.wires q)
    have hw : w ∈ H.map E := by
      refine (Subspace.forall_mem_dualAnnihilator_apply_eq_zero_iff (H.map E) w).1 fun φ hφ => ?_
      rw [Submodule.mem_dualAnnihilator] at hφ
      have hφx : ∀ x : Fin n → ZMod 2, φ x =
          ∑ i ∈ Finset.univ.filter (fun i => φ (fun j => if i = j then 1 else 0) = 1), x i := by
        intro x
        rw [LinearMap.pi_apply_eq_sum_univ φ x]
        simp only [smul_eq_mul]
        exact zmod2_sum_mul _ x _
      set S := Finset.univ.filter (fun i => φ (fun j => if i = j then (1 : ZMod 2) else 0) = 1)
      set ℓ : Lin := ∑ i ∈ S, (ρ.wires i).1
      -- `ℓ` mentions no Hadamard variable.
      have hℓH : ∀ v, n ≤ v → ℓ v = 0 := by
        intro v hv
        have hμ : (fun x => if x = v then (1 : ZMod 2) else 0) ∈ H := fun x hx => by
          simp only [ite_eq_right_iff]; omega
        have := hφ _ (Submodule.mem_map_of_mem hμ)
        rw [hφx] at this
        simp only [E, LinearMap.coe_mk, AddHom.coe_mk, lc_indicator] at this
        simp only [ℓ, Finsupp.finsetSum_apply]
        exact this
      set T := Finset.univ.filter (fun i : Fin n => ℓ i.val = 1)
      have hℓ : ℓ = AffAbs.inputs T := by
        ext v
        by_cases hv : v < n
        · have := sum_singles_apply T ⟨v, hv⟩
          simp only [AffAbs.inputs]
          rw [this]
          simp only [T, Finset.mem_filter, Finset.mem_univ, true_and]
          rcases zmod2_cases (ℓ v) with h | h <;> simp [h]
        · rw [hℓH v (by omega)]
          simp only [AffAbs.inputs, Finsupp.finsetSum_apply, Finsupp.single_apply]
          exact (Finset.sum_eq_zero fun i _ => if_neg (by have := i.isLt; omega)).symm
      set c := ∑ i ∈ S, (ρ.wires i).2
      have hsum : ∑ i ∈ S, ρ.wires i = (AffAbs.inputs T, c) := by
        rw [← hℓ]; exact Prod.ext (Prod.fst_sum ..) (Prod.snd_sum ..)
      have hr := hrel (S, T, c) hsum
      simp only at hr
      rw [hφx]
      simp only [w, Finset.sum_sub_distrib, ← parity_eq_sum_bit, ← Parity.eval_sum, hsum,
        hr]
      rw [eval_inputs s 0]
      ring
    obtain ⟨μ, hμ, hEμ⟩ := Submodule.mem_map.1 hw
    refine ⟨μ, fun q => ?_⟩
    have hass : assign s μ = β + μ := by
      funext v
      by_cases hv : v < n
      · simp [β, assign, hv, hμ v hv]
      · simp [β, assign, hv]
    have hq := congrFun hEμ q
    simp only [E, LinearMap.coe_mk, AddHom.coe_mk, w] at hq
    rw [hass, Parity.eval, lc_add, ← add_assoc, ← Parity.eval, hq]
    ring

/-- Paper, Section 7.1, after Theorem 7.4: "the relations entailed by `Φ_ρ` once the Hadamard
variables are projected away have the form `⊕_{q ∈ S} y_q = c ⊕ ⊕_{q ∈ T} v_q`". A relation of
this form is entailed exactly when it is one of the derived relations `AffAbs.Rel`, i.e. when
`∑_{q ∈ S} ρ(q) = c ⊕ ⊕_{q ∈ T} v_q` syntactically. -/
theorem LinState.entails_iff_rel (ρ : LinState n) (S T : Finset (Fin n)) (c : ZMod 2) :
    (∀ s t, ρ.SatAt s t → parity S t = parity T s + c) ↔ (S, T, c) ∈ ρ.Rel := by
  have hbit : ∀ x : ZMod 2, bit (decide (x = 1)) = x := fun x => by
    rcases zmod2_cases x with rfl | rfl <;> simp [bit]
  constructor
  · intro h
    have key : ∀ ν : ℕ → ZMod 2,
        Parity.eval ν (∑ i ∈ S, ρ.wires i) = Parity.eval ν (AffAbs.inputs T, c) := by
      intro ν
      let s : Basis n := fun i => decide (ν i.val = 1)
      let t : Basis n := fun q => decide (Parity.eval ν (ρ.wires q) = 1)
      have hass : assign s ν = ν := by
        funext v
        by_cases hv : v < n
        · simp [assign, hv, s, hbit]
        · simp [assign, hv]
      have hsat : ρ.SatAt s t := ⟨ν, fun q => by rw [hass]; simp only [t, hbit]⟩
      have e := h s t hsat
      rw [parity_eq_sum_bit] at e
      simp only [t, hbit] at e
      rw [Parity.eval_sum, e]
      conv_rhs => rw [← hass]
      rw [eval_inputs, add_comm]
    have h0 := key 0
    simp [Parity.eval] at h0
    show ∑ i ∈ S, ρ.wires i = (AffAbs.inputs T, c)
    refine Prod.ext ?_ h0
    ext v
    have h1 := key (fun x => if x = v then 1 else 0)
    simp only [Parity.eval, lc_indicator] at h1
    rw [h0] at h1
    exact add_left_cancel h1
  · intro hr s t hst
    exact (ρ.satAt_iff_rel s t).1 hst _ hr

/-- **`γ_L` is `affGamma`.** The paper's concretization of a linear abstract state
(Definition 7.3) is the concretization `affGamma` of `Affine.lean` by the derived relations. -/
theorem gammaL_eq_affGamma (gs : List (CTGate n)) :
    gammaL ((AffAbs.init n).run gs) = affGamma gs := by
  ext K
  simp only [gammaL, supp, affGamma, Set.mem_ofPred_eq, LinState.satAt_iff_rel, Prod.forall]

theorem CTCirc.gammaL_rho (C : CTCirc n) : gammaL C.rho = affGamma C.toList :=
  gammaL_eq_affGamma _

/-! ### Parity relations and `Z` strings -/

/-- `Z_T`, the product of the strings `Z_q` with `q ∈ T`. -/
def Zset (T : Finset (Fin n)) : PP n := PP.ofS (zString T false)

/-- `(-1)^c Z_T`. -/
def sZset (c : ZMod 2) (T : Finset (Fin n)) : PP n := PP.ofS (zString T (decide (c = 1)))

theorem Zset_singleton (q : Fin n) : Zset {q} = Zq q := by
  apply PP.ext'
  · rfl
  · intro r; by_cases h : r = q <;> simp [Zset, PP.ofS, zString, PP.single_letter, h]

theorem sgn_decide (c : ZMod 2) : SPauli.sgn (decide (c = 1)) = zsgn c := by
  rcases zmod2_cases c with rfl | rfl <;> simp [SPauli.sgn, zsgn]

theorem sZset_toMatrix (c : ZMod 2) (T : Finset (Fin n)) :
    (sZset c T).toMatrix = zsgn c • (Zset T).toMatrix := by
  funext o i
  simp only [sZset, Zset, PP.toMatrix_ofS, Matrix.smul_apply, zString_apply, smul_eq_mul]
  split_ifs
  · rcases zmod2_cases c with rfl | rfl <;> simp [SPauli.sgn, zsgn]
  · simp

theorem zsgn_inj {a b : ZMod 2} (h : zsgn a = zsgn b) : a = b := by
  rcases zmod2_cases a with rfl | rfl <;> rcases zmod2_cases b with rfl | rfl <;>
    simp_all [zsgn] <;> norm_num at h

/-- The two sides of `K (±Z_T) = (±Z_S) K`, entry by entry. -/
theorem takes_zString_iff (K : Density n) (S T : Finset (Fin n)) (b b' : Bool) :
    K * (zString T b).toMatrix = (zString S b').toMatrix * K ↔
      ∀ s t, K t s ≠ 0 → parity S t + bit b' = parity T s + bit b := by
  have hl : ∀ t s, (K * (zString T b).toMatrix) t s = K t s * (SPauli.sgn b * zsgn (parity T s)) := by
    intro t s
    rw [Matrix.mul_apply, Finset.sum_eq_single s (fun k _ hk => by rw [zString_apply, if_neg hk,
      mul_zero]) (by simp), zString_apply, if_pos rfl]
  have hr : ∀ t s, ((zString S b').toMatrix * K) t s = SPauli.sgn b' * zsgn (parity S t) * K t s := by
    intro t s
    rw [Matrix.mul_apply, Finset.sum_eq_single t (fun k _ hk => by rw [zString_apply,
      if_neg (Ne.symm hk), zero_mul]) (by simp), zString_apply, if_pos rfl]
  have hsb : ∀ b : Bool, SPauli.sgn b = zsgn (bit b) := fun b => by
    cases b <;> simp [SPauli.sgn, zsgn, bit]
  constructor
  · intro h s t hne
    have e := congrFun (congrFun h t) s
    rw [hl, hr, mul_comm (K t s)] at e
    have e' := mul_right_cancel₀ hne e
    rw [hsb, hsb, ← zsgn_add, ← zsgn_add] at e'
    have := zsgn_inj e'
    rw [add_comm (bit b'), add_comm (bit b)] at this
    exact this.symm
  · intro h
    ext t s
    rw [hl, hr]
    by_cases hne : K t s = 0
    · simp [hne]
    · rw [hsb, hsb, ← zsgn_add, ← zsgn_add, add_comm (bit b), add_comm (bit b'), h s t hne]
      ring

/-- **Paper, Section 7.1, after Theorem 7.4**: for every operator `K` (not only unitaries),
the relation `⊕_{q ∈ S} y_q = c ⊕ ⊕_{q ∈ T} v_q` holds on the support of `K` exactly when `K`
takes `(-1)^c Z_T` to `Z_S`. -/
theorem parity_rel_iff_takes (K : Density n) (S T : Finset (Fin n)) (c : ZMod 2) :
    (∀ st ∈ supp K, parity S st.2 = c + parity T st.1) ↔ Takes K (sZset c T) (Zset S) := by
  rw [Takes, sZset, Zset, PP.toMatrix_ofS, PP.toMatrix_ofS, takes_zString_iff]
  have hb : bit (decide (c = 1)) = c := by
    rcases zmod2_cases c with rfl | rfl <;> simp [bit]
  simp only [hb, show bit false = 0 from rfl, add_zero, supp, Set.mem_ofPred_eq, Prod.forall]
  exact ⟨fun h s t hne => by rw [h s t hne, add_comm], fun h s t hne => by rw [h s t hne, add_comm]⟩

/-! ## Theorem 7.4: the Pauli abstraction refines the linear abstraction -/

/-- The proof of Theorem 7.4: whenever the linear abstraction derives a relation
`⊕_{q ∈ S} y_q = c ⊕ ⊕_{q ∈ T} v_q`, the Pauli abstraction tracks `(-1)^c Z_S` with weakest
precondition `Z_T`. -/
theorem CTCirc.tracked_of_rel (C : CTCirc n) {S T : Finset (Fin n)} {c : ZMod 2}
    (hr : (S, T, c) ∈ C.rho.Rel) :
    (alphaRel C.toPaper).cons.wp (sZset c S) = some (Zset T) ∧
      ∀ A ∈ (alphaRel C.toPaper).axes, PP.Comm (Zset T) A :=
  C.tracked_of_pauliRun _ _ (rel_pauliFact C.toList hr)

theorem sZset_eq_scale (c : ZMod 2) (T : Finset (Fin n)) :
    sZset c T = (Zset T).scale (if c = 1 then 2 else 0) := by
  apply PP.ext'
  · rcases zmod2_cases c with rfl | rfl <;> simp [sZset, Zset, PP.ofS, PP.scale, zString]
  · intro r; rfl

theorem Zset_eq_scale (c : ZMod 2) (T : Finset (Fin n)) :
    Zset T = (sZset c T).scale (if c = 1 then 2 else 0) := by
  rw [sZset_eq_scale, PP.scale_scale]
  rcases zmod2_cases c with rfl | rfl
  · simp
  · exact (PP.scale_zero _).symm

theorem PP.comm_scale_left (P A : PP n) (k : ZMod 4) : PP.Comm (P.scale k) A ↔ PP.Comm P A := by
  rw [PP.comm_iff, PP.comm_iff, PP.toMatrix_scale, Matrix.smul_mul, Matrix.mul_smul]
  exact ⟨fun h => smul_right_injective _ (ipow_ne_zero k) h, fun h => by rw [h]⟩

/-- Section 7.1, after Theorem 7.4: "whenever the linear abstraction derives such a relation,
the Pauli abstraction tracks `Z_S` and derives the same weakest precondition",
`wp(Z_S) = (-1)^c Z_T`. -/
theorem CTCirc.tracks_Zset_of_rel (C : CTCirc n) {S T : Finset (Fin n)} {c : ZMod 2}
    (hr : (S, T, c) ∈ C.rho.Rel) :
    (alphaRel C.toPaper).cons.wp (Zset S) = some (sZset c T) ∧ Zset S ∈ (alpha C.toPaper).F := by
  obtain ⟨hw, hc⟩ := C.tracked_of_rel hr
  have hw' : (alphaRel C.toPaper).cons.wp (Zset S) = some (sZset c T) := by
    rw [Zset_eq_scale c S, Cons.wp_scale, hw, Option.map_some, sZset_eq_scale]
  refine ⟨hw', ?_⟩
  rw [alpha_eq_rel, Dom.F_rel]
  exact ⟨_, hw', fun A hA => by rw [sZset_eq_scale, PP.comm_scale_left]; exact hc A hA⟩

/-- **Paper, Theorem 7.4 (Pauli refines linear)**, against `affGamma`: for every circuit `C`
over `H, S, X, Z, CX, T, T†` (and `S†`, `CZ`), `γ(α(C)) ⊆ affGamma C`, with `γ` the paper's
concretization over all operators. -/
theorem pauli_refines_linear_aff (C : CTCirc n) : (alpha C.toPaper).gamma ⊆ affGamma C.toList := by
  intro K hK x x' hne ⟨S, T, c⟩ hr
  obtain ⟨hw, hc⟩ := C.tracked_of_rel hr
  rw [alpha_eq_rel] at hK
  have h := hK _ _ hw hc
  rw [sZset, Zset, PP.toMatrix_ofS, PP.toMatrix_ofS, takes_zString_iff] at h
  have e := h x x' hne
  have hb : bit (decide (c = 1)) = c := by
    rcases zmod2_cases c with rfl | rfl <;> simp [bit]
  simp only [hb, show bit false = 0 from rfl, add_zero] at e
  calc parity S x' = parity S x' + c + c := by rw [add_assoc, ← two_mul, show (2 : ZMod 2) = 0 from rfl, zero_mul, add_zero]
    _ = parity T x + c := by rw [← e]

/-- **Paper, Theorem 7.4 (Pauli refines linear)**: for every circuit `C` over the gates
`H, S, X, Z, CX, T, T†` (and also `S†`, `CZ`), `γ(α(C)) ⊆ γ_L(ρ_C)`. -/
theorem pauli_refines_linear (C : CTCirc n) : (alpha C.toPaper).gamma ⊆ gammaL C.rho := by
  rw [CTCirc.gammaL_rho]; exact pauli_refines_linear_aff C

/-- Paper, Section 7.1: the linear abstraction is sound, `⟦C⟧ ∈ γ_L(ρ_C)`. -/
theorem linear_soundness (C : CTCirc n) : C.toPaper.sem ∈ gammaL C.rho :=
  pauli_refines_linear C (soundness _)

/-- The linear abstraction is sound for tzap's semantics of the gate list as well. -/
theorem linear_soundness_circ (C : CTCirc n) : circ C.toList ∈ gammaL C.rho := by
  rw [CTCirc.gammaL_rho]; exact aff_sound _

/-! ### The examples of Section 7.1 -/

/-- Section 7.1: after `CX_{1,2}`, the formula is `(y₁ = v₁) ∧ (y₂ = v₁ ⊕ v₂)`. -/
theorem linear_CX :
    (CTCirc.gate (.cx (0 : Fin 2) 1 (by decide))).rho.wires 0 = (Finsupp.single 0 1, 0) ∧
      (CTCirc.gate (.cx (0 : Fin 2) 1 (by decide))).rho.wires 1 =
        (Finsupp.single 1 1 + Finsupp.single 0 1, 0) := by
  refine ⟨?_, ?_⟩ <;> simp [CTCirc.rho, AffAbs.run, AffAbs.step, AffAbs.init]

/-- The circuit `H₁; CX_{1,2}`. -/
def HCX : CTCirc 2 := .seq (.gate (.h 0)) (.gate (.cx 0 1 (by decide)))

/-- Section 7.1: for `H₁; CX_{1,2}`, the formula at cut 2 is `(y₁ = u₁) ∧ (y₂ = u₁ ⊕ v₂)` (the
fresh variable `u₁` is the variable `2`); it entails `y₁ ⊕ y₂ = v₂`, and it says nothing about
`y₁` alone. -/
theorem linear_HCX :
    HCX.rho.wires 0 = (Finsupp.single 2 1, 0) ∧
      HCX.rho.wires 1 = (Finsupp.single 1 1 + Finsupp.single 2 1, 0) ∧
      ({0, 1}, {1}, 0) ∈ HCX.rho.Rel ∧ ∀ T c, ({0}, T, c) ∉ HCX.rho.Rel := by
  have h0 : HCX.rho.wires 0 = (Finsupp.single 2 1, 0) := by
    simp [CTCirc.rho, HCX, AffAbs.run, AffAbs.step, AffAbs.init]
  have h1 : HCX.rho.wires 1 = (Finsupp.single 1 1 + Finsupp.single 2 1, 0) := by
    simp [CTCirc.rho, HCX, AffAbs.run, AffAbs.step, AffAbs.init]
  refine ⟨h0, h1, ?_, ?_⟩
  · show ∑ i ∈ ({0, 1} : Finset (Fin 2)), HCX.rho.wires i = (AffAbs.inputs {1}, 0)
    rw [Finset.sum_pair (by decide), h0, h1]
    refine Prod.ext ?_ rfl
    ext v
    simp only [Prod.fst_add, AffAbs.inputs, Finset.sum_singleton, Finsupp.coe_add, Pi.add_apply,
      Finsupp.single_apply]
    split_ifs <;> first | rfl | omega
  · intro T c hr
    have hr' : ∑ i ∈ ({0} : Finset (Fin 2)), HCX.rho.wires i = (AffAbs.inputs T, c) := hr
    rw [Finset.sum_singleton, h0] at hr'
    have := congrArg (fun p : Affine => p.1 2) hr'
    simp only [Finsupp.single_eq_same, AffAbs.inputs, Finsupp.finsetSum_apply,
      Finsupp.single_apply] at this
    rw [Finset.sum_eq_zero fun i _ => if_neg (by have := i.isLt; omega)] at this
    exact zmod2_one_ne_zero this

/-- The circuit `H₁; H₁` (Lean qubit `0`). -/
def HH : CTCirc 1 := .seq (.gate (.h 0)) (.gate (.h 0))

theorem HH_toPaper : HH.toPaper = gate (gH 0) ⨾ gate (gH 0) := rfl

/-- Section 7.1: on `H₁; H₁` the linear abstraction has `ρ(1) = u₂`, a Hadamard variable,
so its formula relates the output to nothing and `γ_L` allows every operator. -/
theorem gammaL_HH : HH.rho.wires 0 = (Finsupp.single 2 1, 0) ∧ gammaL HH.rho = Set.univ := by
  have hw : HH.rho.wires 0 = (Finsupp.single 2 1, 0) := by
    simp [CTCirc.rho, HH, AffAbs.run, AffAbs.step, AffAbs.init]
  refine ⟨hw, Set.eq_univ_of_forall fun K st _ => ⟨fun _ => bit (st.2 0), fun q => ?_⟩⟩
  rw [Fin.fin_one_eq_zero q, hw]
  simp [Parity.eval, assign]

/-- **Paper, Theorem 7.4, strictness**: on `C = H₁; H₁`, the inclusion is strict; the matrix
`H` is in `γ_L(ρ_C)` but not in `γ(α(C))`, which holds only the multiples of the identity. -/
theorem pauli_refines_linear_strict :
    (alpha HH.toPaper).gamma ⊂ gammaL HH.rho ∧
      (gH (0 : Fin 1)).mat ∈ gammaL HH.rho ∧ (gH (0 : Fin 1)).mat ∉ (alpha HH.toPaper).gamma := by
  have hnot : (gH (0 : Fin 1)).mat ∉ (alpha HH.toPaper).gamma := by
    rw [HH_toPaper, gamma_alpha_HH]
    rintro ⟨c, hc⟩
    have hx := (gH (0 : Fin 1)).xs_spec 0
    rw [hc, conjTranspose_smul, conjTranspose_one, Matrix.smul_mul, Matrix.one_mul,
      Matrix.mul_smul, Matrix.mul_one, smul_smul] at hx
    have e := congrFun (congrFun hx (fun _ => false)) (fun _ => false)
    simp [gH, ofCT, PP.toMatrix, PP.single, Letter.mat] at e
  have hin : (gH (0 : Fin 1)).mat ∈ gammaL HH.rho := by rw [gammaL_HH.2]; trivial
  exact ⟨⟨pauli_refines_linear HH, fun h => hnot (h hin)⟩, hin, hnot⟩

/-! ## Section 7.2: polynomial abstractions -/

theorem PP.toMatrix_ne_zero (P : PP n) : P.toMatrix ≠ 0 := by
  rw [PP.toMatrix_letters]
  exact smul_ne_zero (ipow_ne_zero _) (PP.toMatrix_base_ne_zero _)

/-- If `Q` anticommutes with `P`, then `Q` does not take `P` to `P`. -/
theorem not_takes_of_anticomm {P Q : PP n}
    (h : Q.toMatrix * P.toMatrix = -(P.toMatrix * Q.toMatrix)) : ¬ Takes Q.toMatrix P P := by
  intro ht
  rw [Takes] at ht
  have h2 : (2 : ℂ) • (P.toMatrix * Q.toMatrix) = 0 := by
    rw [two_smul]
    nth_rewrite 1 [← ht]
    rw [h, neg_add_cancel]
  rw [← PP.toMatrix_mul] at h2
  exact PP.toMatrix_ne_zero _ ((smul_eq_zero.1 h2).resolve_left two_ne_zero)

theorem Xq_mat_unitary (q : Fin n) : (Xq q).toMatrix ∈ Unitaries n := by
  show (Xq q).toMatrixᴴ * (Xq q).toMatrix = 1
  have hX : (Xq q).toMatrix = gateUnitary n (CTGate.x q).toGate := by
    rw [PP.toMatrix_single]
    show embed1 n Letter.X.mat q.val = embed1 n x2 q.val
    congr 1
  rw [hX]
  exact gate_adj_mul (CTGate.x q)

theorem Xq_Zq_anticomm (q : Fin n) :
    (Xq q).toMatrix * (Zq q).toMatrix = -((Zq q).toMatrix * (Xq q).toMatrix) := by
  rw [← PP.toMatrix_mul, ← PP.toMatrix_mul, XZ_eq_Y, ZX_eq_Y, PP.toMatrix_scale,
    PP.toMatrix_scale, ipow_three, ipow_one, neg_smul]

/-- `X_q` takes `Z_q` to `-Z_q`. -/
theorem Xq_takes_Zq_neg (q : Fin n) : Takes (Xq q).toMatrix (Zq q) (Zq q).neg := by
  rw [Takes, PP.toMatrix_neg, Matrix.neg_mul]; exact Xq_Zq_anticomm q

/-! ### The Pauli abstraction is more precise: `CX₁₂; H₁; H₂; CX₂₁; H₂` -/

/-- The Clifford circuit `CX_{1,2}; H₁; H₂; CX_{2,1}; H₂` of Section 7.2 (`SF.cexMid`). -/
def cexMidC : CTCirc 2 :=
  .seq (.seq (.seq (.seq (.gate (.cx 0 1 (by decide))) (.gate (.h 0))) (.gate (.h 1)))
    (.gate (.cx 1 0 (by decide)))) (.gate (.h 1))

theorem cexMidC_toList : cexMidC.toList = SF.cexMid := rfl

theorem cexMidC_toPaper : cexMidC.toPaper =
    gate (gCX 0 1 (by decide)) ⨾ gate (gH 0) ⨾ gate (gH 1) ⨾ gate (gCX 1 0 (by decide)) ⨾
      gate (gH 1) := rfl

/-- Section 7.2: the Pauli abstraction of the Clifford circuit (its tableau) has no axes and
`wp(Z₂) = Z₂`, so it tracks `Z₂`, and the circuit takes `Z₂` to `Z₂`. -/
theorem cexMid_pauli :
    (alphaRel cexMidC.toPaper).axes = [] ∧
      (alphaRel cexMidC.toPaper).cons.wp (Zq 1) = some (Zq 1) ∧
      Zq 1 ∈ (alpha cexMidC.toPaper).F ∧
      Takes cexMidC.toPaper.sem (Zq 1) (Zq 1) := by
  have hax : (alphaRel cexMidC.toPaper).axes = [] := by decide
  have hwp : (alphaRel cexMidC.toPaper).cons.wp (Zq 1) = some (Zq 1) := by decide
  refine ⟨hax, hwp, ?_, takes_of_tracked _ hwp (by simp [hax])⟩
  rw [alpha_eq_rel, Dom.F_rel]
  exact ⟨_, hwp, by simp [hax]⟩

/-- Section 7.2: the exact relation of the circuit contains `y₂ = v₂`. -/
theorem cexMid_exact_rel : ∀ st ∈ supp cexMidC.toPaper.sem, st.2 1 = st.1 1 := by
  have h := (parity_rel_iff_takes cexMidC.toPaper.sem {1} {1} 0).2 (by
    rw [show sZset (0 : ZMod 2) {(1 : Fin 2)} = Zq 1 by
      rw [← Zset_singleton]; rfl, Zset_singleton]
    exact cexMid_pauli.2.2.2)
  intro st hst
  have e := h st hst
  simp only [parity, Finset.sum_singleton, zero_add] at e
  revert e; cases st.2 1 <;> cases st.1 1 <;> simp [zmod2_one_ne_zero.symm]

/-- Section 7.2: the path sum of the circuit has three path variables (`2, 3, 4`, after the
inputs `0, 1`), its output bits are `y₁ = u₁ ⊕ u₂` and `y₂ = u₃`, and rewriting cannot
simplify it at any degree. -/
theorem cexMid_pathSum (d : ℕ) :
    (∀ y, y ∈ (SF.run SF.cexMid).temps ↔ y = 2 ∨ y = 3 ∨ y = 4) ∧
      (∀ ν, SF.ev ν ((SF.run SF.cexMid).ket 0) = (ν 2 != ν 3)) ∧
      (∀ ν, SF.ev ν ((SF.run SF.cexMid).ket 1) = ν 4) ∧
      ∀ C, SF.Reduces d (SF.run SF.cexMid) C → C = SF.run SF.cexMid :=
  ⟨SF.mid_temps, SF.mid_ket0, SF.mid_ket1, fun _ h => SF.cexMid_reduces h⟩

/-- **Paper, Section 7.2 ("The Pauli abstraction is more precise")**: on
`CX_{1,2}; H₁; H₂; CX_{2,1}; H₂`, at every degree `d`, the polynomial abstraction relates
every input to every output, so its concretization contains every unitary, in particular
`X₂`; `X₂` takes `Z₂` to `-Z₂`, so it is not in `γ(α(C))`. -/
theorem pauli_more_precise (d : ℕ) :
    (∀ U ∈ Unitaries 2, U ∈ alphaGammaC d cexMidC.toList) ∧
      Takes (Xq (1 : Fin 2)).toMatrix (Zq 1) (Zq 1).neg ∧
      (Xq (1 : Fin 2)).toMatrix ∈ Unitaries 2 ∧
      (Xq (1 : Fin 2)).toMatrix ∈ alphaGammaC d cexMidC.toList ∧
      (Xq (1 : Fin 2)).toMatrix ∉ (alpha cexMidC.toPaper).gamma := by
  have hall : ∀ U ∈ Unitaries 2, U ∈ alphaGammaC d cexMidC.toList := fun U hU =>
    alpha_cexMid d U (Unitaries.mul_conj hU)
  refine ⟨hall, Xq_takes_Zq_neg 1, Xq_mat_unitary 1, hall _ (Xq_mat_unitary 1), fun hX => ?_⟩
  obtain ⟨hax, hwp, -⟩ := cexMid_pauli
  rw [alpha_eq_rel] at hX
  exact not_takes_of_anticomm (Xq_Zq_anticomm 1) (hX _ _ hwp (by simp [hax]))

/-! ### The polynomial abstraction is more precise: `H; T; T†; H` -/

/-- The circuit `H₁; T₁; T†₁; H₁` of Section 7.2 (`SF.midA`). -/
def midAC : CTCirc 1 :=
  .seq (.seq (.seq (.gate (.h 0)) (.gate (.t 0))) (.gate (.tdg 0))) (.gate (.h 0))

theorem midAC_toList : midAC.toList = SF.midA := rfl

theorem midAC_toPaper : midAC.toPaper = gate (gH 0) ⨾ T 0 ⨾ Tdg 0 ⨾ gate (gH 0) := rfl

/-- Section 7.2: the Pauli abstraction of `H; T; T†; H` records the axis `X₁` for both
rotations, and its constraints are those of the identity. Since `Z₁` does not commute with
`X₁`, it does not track `Z₁`. -/
theorem midA_alpha :
    alpha midAC.toPaper = .rel ⟨Cons.id 1, [Xq 0, Xq 0]⟩ ∧ Zq 0 ∉ (alpha midAC.toPaper).F := by
  have h : alpha midAC.toPaper = .rel ⟨Cons.id 1, [Xq 0, Xq 0]⟩ := by decide
  refine ⟨h, ?_⟩
  rw [h, Dom.F_rel]
  decide

/-- Section 7.2: `γ(α(H; T; T†; H))` contains `X₁` (which takes `Z₁` to `-Z₁`). -/
theorem midA_X_mem : (Xq (0 : Fin 1)).toMatrix ∈ (alpha midAC.toPaper).gamma := by
  rw [midA_alpha.1, Dom.gamma_rel]
  intro P B hB hc
  rw [Cons.wp_id] at hB
  cases hB
  exact (hc _ (by simp)).toMatrix.symm

/-- Section 7.2: rewriting lets the polynomial abstraction derive `y₁ = v₁`, at every degree
`d ≥ 1`: every operator in its concretization takes `Z₁` to `Z₁`, and its support satisfies
`y₁ = v₁`. -/
theorem midA_polynomial (d : ℕ) (hd : 1 ≤ d) (U : Density 1) (hU : U ∈ alphaGammaC d midAC.toList) :
    Takes U (Zq 0) (Zq 0) ∧ ∀ st ∈ supp U, st.2 0 = st.1 0 := by
  have hmem : (zString {(0 : Fin 1)} false, zString {0} false) ∈ sfFacts d SF.midA :=
    ⟨0, 0, false, rfl, rfl, SF.midA_proves d hd⟩
  have hf := alpha_refines_sf d SF.midA hU _ hmem
  simp only [Holds, conj] at hf
  have hUU : Uᴴ * U = 1 := mul_eq_one_comm.mp hU.1
  have hc : U * (zString {(0 : Fin 1)} false).toMatrix = (zString {0} false).toMatrix * U := by
    calc U * (zString {(0 : Fin 1)} false).toMatrix
        = U * (zString {(0 : Fin 1)} false).toMatrix * (Uᴴ * U) := by rw [hUU, Matrix.mul_one]
      _ = _ := by rw [← Matrix.mul_assoc, hf]
  have hZ : (zString {(0 : Fin 1)} false).toMatrix = (Zq 0).toMatrix := by
    rw [← PP.toMatrix_ofS, ← Zset_singleton]; rfl
  refine ⟨by rw [Takes, ← hZ]; exact hc, fun st hst => ?_⟩
  have e := (takes_zString_iff U {0} {0} false false).1 hc st.1 st.2 hst
  simp only [parity, Finset.sum_singleton, bit, Bool.false_eq_true, if_false, add_zero] at e
  revert e; cases st.2 0 <;> cases st.1 0 <;> simp [zmod2_one_ne_zero.symm]

/-- **Paper, Section 7.2 ("The polynomial abstraction is more precise")**: on `H; T; T†; H`,
at every degree `d ≥ 1`, the unitary `X₁` is in `γ(α(C))` but not in the polynomial
abstraction's concretization. -/
theorem polynomial_more_precise (d : ℕ) (hd : 1 ≤ d) :
    (Xq (0 : Fin 1)).toMatrix ∈ Unitaries 1 ∧
      (Xq (0 : Fin 1)).toMatrix ∈ (alpha midAC.toPaper).gamma ∧
      Takes (Xq (0 : Fin 1)).toMatrix (Zq 0) (Zq 0).neg ∧
      (Xq (0 : Fin 1)).toMatrix ∉ alphaGammaC d midAC.toList :=
  ⟨Xq_mat_unitary 0, midA_X_mem, Xq_takes_Zq_neg 0, fun h =>
    not_takes_of_anticomm (Xq_Zq_anticomm 0) (midA_polynomial d hd _ h).1⟩

/-- **Paper, Section 7.2: "Neither abstraction refines the other, at any degree `d ≥ 1`."**
For every `d ≥ 1` there is a circuit with a unitary in the polynomial abstraction's
concretization `alphaGammaC d` but not in the Pauli abstraction's `γ(α(C))`
(`CX_{1,2}; H₁; H₂; CX_{2,1}; H₂`, with `X₂`), and one with a unitary in `γ(α(C))` but not in
`alphaGammaC d` (`H; T; T†; H`, with `X₁`).

The witnesses are unitaries: `alphaGammaC` contains only unitaries while `γ` ranges over all
operators, so a non-unitary operator such as `0 ∈ γ(α(C))` would separate them trivially. -/
theorem pauli_polynomial_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (C : CTCirc n) (U : Density n), U ∈ Unitaries n ∧
        U ∈ alphaGammaC d C.toList ∧ U ∉ (alpha C.toPaper).gamma) ∧
      (∃ (n : ℕ) (C : CTCirc n) (U : Density n), U ∈ Unitaries n ∧
        U ∈ (alpha C.toPaper).gamma ∧ U ∉ alphaGammaC d C.toList) := by
  obtain ⟨-, -, hu, hin, hout⟩ := pauli_more_precise d
  obtain ⟨hu', hin', -, hout'⟩ := polynomial_more_precise d hd
  exact ⟨⟨2, cexMidC, _, hu, hin, hout⟩, ⟨1, midAC, _, hu', hin', hout'⟩⟩

/-- The same, as non-inclusions of the concretizations restricted to unitaries. -/
theorem pauli_polynomial_incomparable' (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (C : CTCirc n), ¬ alphaGammaC d C.toList ⊆ (alpha C.toPaper).gamma) ∧
      (∃ (n : ℕ) (C : CTCirc n),
        ¬ (alpha C.toPaper).gamma ∩ Unitaries n ⊆ alphaGammaC d C.toList) := by
  obtain ⟨⟨n, C, U, -, hin, hout⟩, ⟨n', C', U', hu', hin', hout'⟩⟩ :=
    pauli_polynomial_incomparable d hd
  exact ⟨⟨n, C, fun h => hout (h hin)⟩, ⟨n', C', fun h => hout' (h ⟨hin', hu'⟩)⟩⟩

end

end PauliFold.Paper
