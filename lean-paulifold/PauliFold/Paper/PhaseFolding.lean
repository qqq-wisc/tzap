import PauliFold.Paper.Circuits

/-!
# Paper, Section 5: phase folding

* Theorem 5.1 (Phase folding): `phase_folding` (and `phase_folding_eq`, as an equality), with
  the matrix identity behind it, `fold_sem`.
* Algorithm 1: `Gate`, `ofList` (the circuit `I; g₁; ⋯; g_m`), the rotations of a gate list
  (`numRots`, `splitRot`), the axis `A_j` of the `j`th rotation (`axis`), the fold condition of
  line 3 (`Folds`), the choice of `j` (`Choice`, `Choice.Valid`), one step of the loop (`step`),
  and the algorithm (`phaseFold`).
* "So `𝒜 = {A_1, …, A_r}`": `axes_ofList`.
* "The condition on line 3 is the hypothesis of Theorem 5.1 for the part of `D` after its `j`th
  rotation": `folds_iff`.
* Theorem 5.2 (Correctness): `phaseFold_correct`.
* A valid choice function (the latest `j` that folds): `latest`, `latest_valid`.
* The worked example of Figure 1 / Section 2.3 (namespace `Figure1`): `cuts`, `figure1`,
  `figure1_folds`, `figure1_phaseFold`, `figure1_latest`.

Qubits are 0-indexed: the paper's qubits 1 and 2 are the Lean qubits `0` and `1`. Rotations of
`D` are also counted from `0`: the paper's `j`th rotation is the Lean rotation `j - 1`.
-/

namespace PauliFold.Paper

open TzapLean Matrix Circuit

noncomputable section

variable {n : ℕ}

/-! ## Signs -/

/-- `s P` for a sign `s ∈ {1, -1}`: `P` when `s = 1`, and `-P` otherwise. -/
def signed (s : ℤ) (P : PP n) : PP n := if s = 1 then P else P.neg

theorem signed_toMatrix {s : ℤ} (hs : s = 1 ∨ s = -1) (P : PP n) :
    (signed s P).toMatrix = (s : ℂ) • P.toMatrix := by
  rcases hs with rfl | rfl
  · simp [signed]
  · simp [signed, PP.toMatrix_neg]

/-! ## Theorem 5.1 -/

/-- If `U (s Z_q) = Z_{q'} U`, then a rotation about `Z_{q'}` after `U` is a rotation about
`Z_q` before `U`, by the angle times `s`: `R_{q'}(θ') U = U R_q(sθ')`. -/
theorem rz_slide {U : Density n} {q q' : Fin n} {s : ℤ} (hs : s = 1 ∨ s = -1)
    (h : U * (signed s (Zq q)).toMatrix = (Zq q').toMatrix * U) (θ' : ℝ) :
    rzMat q' θ' * U = U * rzMat q (s * θ') := by
  rw [signed_toMatrix hs] at h
  rw [rzMat_eq, rzMat_eq, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, ← h, Matrix.mul_smul,
    smul_smul]
  rcases hs with rfl | rfl
  · simp
  · congr 1
    · congr 1; push_cast; rw [neg_one_mul, neg_div, Complex.cos_neg]
    · congr 1; push_cast; rw [neg_one_mul, neg_div, Complex.sin_neg]; ring

/-- The matrix identity behind Theorem 5.1: under its hypotheses,
`R_{q'}(θ') ⟦R_q(θ); C⟧ = ⟦R_q(θ + sθ'); C⟧`. -/
theorem fold_sem (C : Circuit n) (q q' : Fin n) (s : ℤ) (hs : s = 1 ∨ s = -1) (θ : ℝ)
    (σ : ARel n) (hσ : alpha (rz q θ ⨾ C) = .rel σ)
    (hF : Zq q' ∈ σ.F) (hwp : σ.cons.wp (Zq q') = some (signed s (Zq q))) (θ' : ℝ) :
    rzMat q' θ' * sem (rz q θ ⨾ C) = sem (rz q (θ + s * θ') ⨾ C) := by
  have hsound := soundness (rz q θ ⨾ C)
  rw [hσ] at hsound
  obtain ⟨B, hB, hc⟩ := hF
  rw [hwp] at hB
  cases hB
  have htakes := hsound (Zq q') (signed s (Zq q)) hwp hc
  rw [rz_slide hs htakes, sem_seq, sem_seq, sem_rz, sem_rz, Matrix.mul_assoc, rzMat_mul]

/-- **Paper, Theorem 5.1 (Phase folding).** Let `C` be a circuit, `q, q'` qubits, `s ∈ {1, -1}`,
`θ ∈ ℝ`, and `σ = (Γ, 𝒜) = α(R_q(θ); C)`. If `Z_{q'} ∈ 𝓕(σ)` and `wp_Γ(Z_{q'}) = s Z_q`, then
for all circuits `C₁, C₂` and every `θ' ∈ ℝ`,
`C₁; R_q(θ); C; R_{q'}(θ'); C₂ ≡ C₁; R_q(θ + sθ'); C; C₂`.
(`α` returns an element of `Dom`; `hσ` names the abstract relation it returns, which is never
`⊤`. The two sides are in fact equal, `phase_folding_eq`.) -/
theorem phase_folding (C : Circuit n) (q q' : Fin n) (s : ℤ) (hs : s = 1 ∨ s = -1) (θ : ℝ)
    (σ : ARel n) (hσ : alpha (rz q θ ⨾ C) = .rel σ)
    (hF : Zq q' ∈ σ.F) (hwp : σ.cons.wp (Zq q') = some (signed s (Zq q)))
    (C₁ C₂ : Circuit n) (θ' : ℝ) :
    (C₁ ⨾ rz q θ ⨾ C ⨾ rz q' θ' ⨾ C₂).Equiv (C₁ ⨾ rz q (θ + s * θ') ⨾ C ⨾ C₂) := by
  refine EquivP.of_eq ?_
  have h := fold_sem C q q' s hs θ σ hσ hF hwp θ'
  simp only [sem_seq, sem_rz] at h ⊢
  calc _ = C₂.sem * ((rzMat q' θ' * (C.sem * rzMat q θ)) * C₁.sem) := by
        simp only [Matrix.mul_assoc]
    _ = _ := by rw [h]; simp only [Matrix.mul_assoc]

/-- Theorem 5.1, as an equality of matrices. -/
theorem phase_folding_eq (C : Circuit n) (q q' : Fin n) (s : ℤ) (hs : s = 1 ∨ s = -1) (θ : ℝ)
    (σ : ARel n) (hσ : alpha (rz q θ ⨾ C) = .rel σ)
    (hF : Zq q' ∈ σ.F) (hwp : σ.cons.wp (Zq q') = some (signed s (Zq q)))
    (C₁ C₂ : Circuit n) (θ' : ℝ) :
    sem (C₁ ⨾ rz q θ ⨾ C ⨾ rz q' θ' ⨾ C₂) = sem (C₁ ⨾ rz q (θ + s * θ') ⨾ C ⨾ C₂) := by
  have h := fold_sem C q q' s hs θ σ hσ hF hwp θ'
  simp only [sem_seq, sem_rz] at h ⊢
  calc _ = C₂.sem * ((rzMat q' θ' * (C.sem * rzMat q θ)) * C₁.sem) := by
        simp only [Matrix.mul_assoc]
    _ = _ := by rw [h]; simp only [Matrix.mul_assoc]

/-! ## Gate lists -/

/-- A gate of Algorithm 1: a Clifford gate or a Z-rotation. -/
inductive Gate (n : ℕ) where
  | cliff (G : CliffordGate n)
  | rz (q : Fin n) (θ : ℝ)

/-- A gate as a circuit. -/
def Gate.toCircuit : Gate n → Circuit n
  | .cliff G => .gate G
  | .rz q θ => .rz q θ

/-- `C₀; g₁; ⋯; g_m`. -/
def ofListFrom (C₀ : Circuit n) (l : List (Gate n)) : Circuit n :=
  l.foldl (fun D g => D ⨾ g.toCircuit) C₀

/-- The circuit `I; g₁; ⋯; g_m` of a list of gates (Algorithm 1 starts from `D = I`). -/
def ofList (l : List (Gate n)) : Circuit n := ofListFrom (gate (CliffordGate.one n)) l

/-- The gates `g₁, …, g_m` of a circuit `C = g₁; ⋯; g_m`. -/
def Circuit.gates : Circuit n → List (Gate n)
  | gate G => [.cliff G]
  | rz q θ => [.rz q θ]
  | seq C₁ C₂ => C₁.gates ++ C₂.gates

@[simp] theorem ofListFrom_nil (C₀ : Circuit n) : ofListFrom C₀ [] = C₀ := rfl
@[simp] theorem ofListFrom_cons (C₀ : Circuit n) (g : Gate n) (l : List (Gate n)) :
    ofListFrom C₀ (g :: l) = ofListFrom (C₀ ⨾ g.toCircuit) l := rfl

theorem ofListFrom_append (C₀ : Circuit n) (l₁ l₂ : List (Gate n)) :
    ofListFrom C₀ (l₁ ++ l₂) = ofListFrom (ofListFrom C₀ l₁) l₂ := List.foldl_append

theorem ofList_concat (l : List (Gate n)) (g : Gate n) :
    ofList (l ++ [g]) = ofList l ⨾ g.toCircuit := by
  rw [ofList, ofListFrom_append]; rfl

theorem sem_ofListFrom (C₀ : Circuit n) (l : List (Gate n)) :
    sem (ofListFrom C₀ l) = sem (ofList l) * sem C₀ := by
  induction l generalizing C₀ with
  | nil => simp [ofList, CliffordGate.one]
  | cons g l ih =>
    rw [ofListFrom_cons, ofList, ofListFrom_cons, ih, ih (gate _ ⨾ _)]
    simp [CliffordGate.one, Matrix.mul_assoc]

theorem sem_ofList_append (l₁ l₂ : List (Gate n)) :
    sem (ofList (l₁ ++ l₂)) = sem (ofList l₂) * sem (ofList l₁) := by
  rw [ofList, ofListFrom_append, sem_ofListFrom]; rfl

theorem sem_ofList_cons (g : Gate n) (l : List (Gate n)) :
    sem (ofList (g :: l)) = sem (ofList l) * sem g.toCircuit := by
  rw [ofList, ofListFrom_cons, sem_ofListFrom]; simp [CliffordGate.one]

theorem cliff_ofListFrom (C₀ : Circuit n) (l : List (Gate n)) :
    cliff (ofListFrom C₀ l) = cliff (ofList l) * cliff C₀ := by
  induction l generalizing C₀ with
  | nil => simp [ofList, cliff, CliffordGate.one]
  | cons g l ih =>
    rw [ofListFrom_cons, ofList, ofListFrom_cons, ih, ih (gate _ ⨾ _)]
    simp [cliff, CliffordGate.one, Matrix.mul_assoc]

theorem cliff_ofList_append (l₁ l₂ : List (Gate n)) :
    cliff (ofList (l₁ ++ l₂)) = cliff (ofList l₂) * cliff (ofList l₁) := by
  rw [ofList, ofListFrom_append, cliff_ofListFrom]; rfl

theorem cliff_ofList_cons (g : Gate n) (l : List (Gate n)) :
    cliff (ofList (g :: l)) = cliff (ofList l) * cliff g.toCircuit := by
  rw [ofList, ofListFrom_cons, cliff_ofListFrom]; simp [cliff, CliffordGate.one]

/-- `I; g₁; ⋯; g_m` is equal to `C = g₁; ⋯; g_m`. -/
theorem sem_ofList_gates : ∀ C : Circuit n, sem (ofList C.gates) = sem C
  | gate G => by simp [Circuit.gates, ofList, Gate.toCircuit, CliffordGate.one]
  | rz q θ => by simp [Circuit.gates, ofList, Gate.toCircuit, CliffordGate.one]
  | seq C₁ C₂ => by rw [Circuit.gates, sem_ofList_append, sem_ofList_gates, sem_ofList_gates]; rfl

/-- The abstraction of `C₀; g₁; ⋯; g_m` depends on `C₀` only through `α(C₀)`. -/
theorem alpha_ofListFrom_congr {C C' : Circuit n} (h : alpha C = alpha C') (l : List (Gate n)) :
    alpha (ofListFrom C l) = alpha (ofListFrom C' l) := by
  induction l generalizing C C' with
  | nil => exact h
  | cons g l ih => exact ih (by rw [alpha_seq, alpha_seq, h])

/-! ## The rotations of a gate list -/

/-- The number `r` of rotations. -/
def numRots : List (Gate n) → ℕ
  | [] => 0
  | .cliff _ :: l => numRots l
  | .rz _ _ :: l => numRots l + 1

/-- The `j`th rotation `R_{q_j}(θ_j)` of a gate list (counted from `0`), with the gates before
it and the gates after it: `splitRot D j = some (pre, q_j, θ_j, post)` when
`D = pre ++ [R_{q_j}(θ_j)] ++ post` and `pre` has `j` rotations. -/
def splitRot : List (Gate n) → ℕ → Option (List (Gate n) × Fin n × ℝ × List (Gate n))
  | [], _ => none
  | .cliff G :: l, j => (splitRot l j).map fun x => (.cliff G :: x.1, x.2)
  | .rz q θ :: l, 0 => some ([], q, θ, l)
  | .rz q θ :: l, j + 1 => (splitRot l j).map fun x => (.rz q θ :: x.1, x.2)

theorem splitRot_spec : ∀ {l : List (Gate n)} {j : ℕ} {p : List (Gate n)} {q : Fin n} {θ : ℝ}
    {r : List (Gate n)}, splitRot l j = some (p, q, θ, r) →
      l = p ++ .rz q θ :: r ∧ numRots p = j
  | [], _, _, _, _, _, h => by cases h
  | .cliff G :: l, j, p, q, θ, r, h => by
    simp only [splitRot, Option.map_eq_some_iff] at h
    obtain ⟨⟨p', q', θ', r'⟩, h', e⟩ := h
    simp only [Prod.mk.injEq] at e
    obtain ⟨rfl, rfl, rfl, rfl⟩ := e
    obtain ⟨rfl, rfl⟩ := splitRot_spec h'
    exact ⟨rfl, rfl⟩
  | .rz q' θ' :: l, 0, p, q, θ, r, h => by
    simp only [splitRot, Option.some.injEq, Prod.mk.injEq] at h
    obtain ⟨rfl, rfl, rfl, rfl⟩ := h
    exact ⟨rfl, rfl⟩
  | .rz q' θ' :: l, j + 1, p, q, θ, r, h => by
    simp only [splitRot, Option.map_eq_some_iff] at h
    obtain ⟨⟨p', q'', θ'', r'⟩, h', e⟩ := h
    simp only [Prod.mk.injEq] at e
    obtain ⟨rfl, rfl, rfl, rfl⟩ := e
    obtain ⟨rfl, rfl⟩ := splitRot_spec h'
    exact ⟨rfl, rfl⟩

theorem splitRot_isSome : ∀ (l : List (Gate n)) (j : ℕ), (splitRot l j).isSome ↔ j < numRots l
  | [], j => by simp [splitRot, numRots]
  | .cliff G :: l, j => by simp [splitRot, numRots, splitRot_isSome l j]
  | .rz q θ :: l, 0 => by simp [splitRot, numRots]
  | .rz q θ :: l, j + 1 => by simp [splitRot, numRots, splitRot_isSome l j]

theorem splitRot_append_left : ∀ (l₁ l₂ : List (Gate n)) (j : ℕ), j < numRots l₁ →
    splitRot (l₁ ++ l₂) j = (splitRot l₁ j).map fun x => (x.1, x.2.1, x.2.2.1, x.2.2.2 ++ l₂)
  | [], _, j, h => by simp [numRots] at h
  | .cliff G :: l, l₂, j, h => by
    simp only [List.cons_append, splitRot, Option.map_map]
    rw [splitRot_append_left l l₂ j h, Option.map_map]; rfl
  | .rz q θ :: l, l₂, 0, _ => rfl
  | .rz q θ :: l, l₂, j + 1, h => by
    simp only [List.cons_append, splitRot, Option.map_map]
    rw [splitRot_append_left l l₂ j (by simpa [numRots] using h), Option.map_map]; rfl

theorem splitRot_append_right : ∀ (l₁ l₂ : List (Gate n)) (k : ℕ),
    splitRot (l₁ ++ l₂) (numRots l₁ + k) = (splitRot l₂ k).map fun x => (l₁ ++ x.1, x.2)
  | [], l₂, k => by simp [numRots]
  | .cliff G :: l, l₂, k => by
    simp only [List.cons_append, splitRot, numRots]
    rw [splitRot_append_right l l₂ k, Option.map_map]; rfl
  | .rz q θ :: l, l₂, k => by
    simp only [List.cons_append, numRots, show numRots l + 1 + k = (numRots l + k) + 1 by omega,
      splitRot]
    rw [splitRot_append_right l l₂ k, Option.map_map]; rfl

/-- Add `φ` to the angle of the `j`th rotation: `θ_j ← θ_j + φ`. -/
def addAngle (D : List (Gate n)) (j : ℕ) (φ : ℝ) : List (Gate n) :=
  match splitRot D j with
  | some (p, q, θ, r) => p ++ .rz q (θ + φ) :: r
  | none => D

/-! ## Axes -/

/-- The axis `A_j = wp_{Γ_j}(Z_{q_j})` of the `j`th rotation `R_{q_j}(θ_j)` of `D`, where
`Γ_j` are the constraints of the abstraction of the part of `D` before the `j`th rotation
(`none` if `D` has no `j`th rotation). -/
def axis (D : List (Gate n)) (j : ℕ) : Option (PP n) :=
  (splitRot D j).bind fun x => (alphaRel (ofList x.1)).cons.wp (Zq x.2.1)

theorem axis_of_splitRot {D p r : List (Gate n)} {j : ℕ} {q : Fin n} {θ : ℝ}
    (h : splitRot D j = some (p, q, θ, r)) :
    ∃ A, axis D j = some A ∧
      A.toMatrix = (cliff (ofList p))ᴴ * (Zq q).toMatrix * cliff (ofList p) := by
  obtain ⟨A, hA, hm⟩ := wp_alpha (ofList p) (Zq q)
  exact ⟨A, by rw [axis, h]; exact hA, hm⟩

theorem axis_isSome (D : List (Gate n)) (j : ℕ) : (axis D j).isSome ↔ j < numRots D := by
  rw [← splitRot_isSome]
  cases h : splitRot D j with
  | none => simp [axis, h]
  | some x =>
    obtain ⟨p, q, θ, r⟩ := x
    obtain ⟨A, hA, -⟩ := axis_of_splitRot h
    simp [hA]

theorem axis_append_left (l₁ l₂ : List (Gate n)) (j : ℕ) (h : j < numRots l₁) :
    axis (l₁ ++ l₂) j = axis l₁ j := by
  rw [axis, axis, splitRot_append_left l₁ l₂ j h, Option.bind_map]; rfl

/-- `α(D)` for `D` with one more gate. -/
theorem alphaRel_ofList_concat (l : List (Gate n)) (g : Gate n) :
    alphaRel (ofList (l ++ [g])) = ⟨(alphaRel g.toCircuit).cons.comp (alphaRel (ofList l)).cons,
      (alphaRel (ofList l)).axes ++
        (alphaRel g.toCircuit).axes.filterMap (alphaRel (ofList l)).cons.wp⟩ := by
  rw [ofList_concat, alphaRel_seq]

/-- **Paper, Section 5: "So `𝒜 = {A_1, …, A_r}`."** The axes of `α(D)` (a list) are exactly the
axes `A_0, …, A_{r-1}` of the rotations of `D`, in order. -/
theorem axes_ofList (D : List (Gate n)) :
    (alphaRel (ofList D)).axes.length = numRots D ∧
      ∀ k, (alphaRel (ofList D)).axes[k]? = axis D k := by
  induction D using List.reverseRecOn with
  | nil => exact ⟨rfl, fun k => by simp [ofList, axis, splitRot]⟩
  | append_singleton l g ih =>
    obtain ⟨hlen, hget⟩ := ih
    have hnum : ∀ k, k < numRots l → axis (l ++ [g]) k = axis l k :=
      fun k hk => axis_append_left l [g] k hk
    have hright : ∀ k, axis (l ++ [g]) (numRots l + k) =
        (splitRot [g] k).bind fun x => (alphaRel (ofList (l ++ x.1))).cons.wp (Zq x.2.1) := by
      intro k; rw [axis, splitRot_append_right, Option.bind_map]; rfl
    cases g with
    | cliff G =>
      have hax : (alphaRel (ofList (l ++ [Gate.cliff G]))).axes = (alphaRel (ofList l)).axes := by
        rw [alphaRel_ofList_concat]; simp [Gate.toCircuit]
      have hn : numRots (l ++ [Gate.cliff G]) = numRots l := by
        rw [show numRots (l ++ [Gate.cliff G]) = numRots l from ?_]
        clear hlen hget hnum hright hax
        induction l with
        | nil => rfl
        | cons g l ih => cases g <;> simp [numRots, ih]
      refine ⟨by rw [hax, hlen, hn], fun k => ?_⟩
      rw [hax, hget]
      by_cases hk : k < numRots l
      · rw [hnum k hk]
      · obtain ⟨k', rfl⟩ := Nat.exists_eq_add_of_le (not_lt.1 hk)
        rw [hright]
        have : axis l (numRots l + k') = none := by
          rw [← Option.not_isSome_iff_eq_none, axis_isSome]; omega
        rw [this]; cases k' <;> rfl
    | rz q θ =>
      obtain ⟨A, hA, -⟩ := wp_alpha (ofList l) (Zq q)
      have hax : (alphaRel (ofList (l ++ [Gate.rz q θ]))).axes = (alphaRel (ofList l)).axes ++ [A] := by
        rw [alphaRel_ofList_concat]; simp [Gate.toCircuit, hA]
      have hn : numRots (l ++ [Gate.rz q θ]) = numRots l + 1 := by
        clear hlen hget hnum hright hax hA
        induction l with
        | nil => rfl
        | cons g l ih => cases g <;> simp [numRots, ih]
      refine ⟨by rw [hax, List.length_append, hlen, hn]; rfl, fun k => ?_⟩
      rw [hax]
      by_cases hk : k < numRots l
      · rw [hnum k hk, List.getElem?_append_left (by omega), hget]
      · obtain ⟨k', rfl⟩ := Nat.exists_eq_add_of_le (not_lt.1 hk)
        rw [hright, List.getElem?_append_right (by omega), hlen, Nat.add_sub_cancel_left]
        cases k' with
        | zero =>
          simp only [splitRot, Option.bind_some, List.append_nil, List.getElem?_cons_zero]
          exact hA.symm
        | succ k' => simp [splitRot]

theorem mem_axes_ofList (D : List (Gate n)) (A : PP n) :
    A ∈ (alphaRel (ofList D)).axes ↔ ∃ k, axis D k = some A := by
  rw [List.mem_iff_getElem?]
  simp only [(axes_ofList D).2]

/-! ## Weakest preconditions through `α(D)`, as conjugations -/

/-- A weakest precondition through `α(C)` is determined by its matrix. -/
theorem wp_alpha_eq (C : Circuit n) {P B : PP n}
    (h : B.toMatrix = (cliff C)ᴴ * P.toMatrix * cliff C) : (alphaRel C).cons.wp P = some B := by
  obtain ⟨B', hB', hm⟩ := wp_alpha C P
  rw [hB', PP.toMatrix_injective (hm.trans h.symm)]

theorem wp_alpha_congr {C C' : Circuit n} (h : cliff C = cliff C') (P : PP n) :
    (alphaRel C).cons.wp P = (alphaRel C').cons.wp P := by
  obtain ⟨B, hB, hm⟩ := wp_alpha C P
  rw [hB, wp_alpha_eq C' (by rw [hm, h])]

theorem wp_alpha_inj (C : Circuit n) {P Q A : PP n} (hP : (alphaRel C).cons.wp P = some A)
    (hQ : (alphaRel C).cons.wp Q = some A) : P = Q := by
  have hV := (models_alphaRel C).mul_conj
  have e := (Cons.wp_toMatrix (models_alphaRel C) hP).symm.trans
    (Cons.wp_toMatrix (models_alphaRel C) hQ)
  apply PP.toMatrix_injective
  have e' := congrArg (fun M => cliff C * M * (cliff C)ᴴ) e
  simp only [← Matrix.mul_assoc, hV, Matrix.one_mul] at e'
  simpa only [Matrix.mul_assoc, hV, Matrix.mul_one] using e'

/-- Weakest preconditions through `α(C)` preserve and reflect commutation. -/
theorem comm_wp_iff_alpha (C : Circuit n) {P Q A B : PP n} (hA : (alphaRel C).cons.wp P = some A)
    (hB : (alphaRel C).cons.wp Q = some B) : PP.Comm A B ↔ PP.Comm P Q := by
  rw [PP.comm_iff, PP.comm_iff, Cons.wp_toMatrix (models_alphaRel C) hA,
    Cons.wp_toMatrix (models_alphaRel C) hB]
  exact comm_conj_iff (models_alphaRel C).unitary

theorem wp_alpha_signed (C : Circuit n) {s : ℤ} (hs : s = 1 ∨ s = -1) {P A : PP n}
    (h : (alphaRel C).cons.wp P = some A) :
    (alphaRel C).cons.wp (signed s P) = some (signed s A) := by
  apply wp_alpha_eq
  rw [signed_toMatrix hs, signed_toMatrix hs, Cons.wp_toMatrix (models_alphaRel C) h,
    Matrix.mul_smul, Matrix.smul_mul]

theorem comm_signed_iff {s : ℤ} (P Q : PP n) : PP.Comm (signed s P) Q ↔ PP.Comm P Q := by
  unfold signed; split_ifs
  · rfl
  · exact PP.comm_scale P Q 2

/-- The weakest precondition through `I; g₁; ⋯; g_m` is computed backwards through the two
halves of the list. -/
theorem wp_ofList_append (l₁ l₂ : List (Gate n)) (P : PP n) :
    (alphaRel (ofList (l₁ ++ l₂))).cons.wp P =
      ((alphaRel (ofList l₂)).cons.wp P).bind (alphaRel (ofList l₁)).cons.wp := by
  obtain ⟨B₂, h₂, m₂⟩ := wp_alpha (ofList l₂) P
  obtain ⟨B₁, h₁, m₁⟩ := wp_alpha (ofList l₁) B₂
  rw [h₂, Option.bind_some, h₁]
  apply wp_alpha_eq
  rw [m₁, m₂, cliff_ofList_append, conjTranspose_mul]
  simp only [Matrix.mul_assoc]

theorem wp_ofList_rz_cons (q : Fin n) (θ : ℝ) (l : List (Gate n)) (P : PP n) :
    (alphaRel (ofList (.rz q θ :: l))).cons.wp P = (alphaRel (ofList l)).cons.wp P :=
  wp_alpha_congr (by rw [cliff_ofList_cons]; simp [Gate.toCircuit, cliff]) P

theorem numRots_append : ∀ l₁ l₂ : List (Gate n), numRots (l₁ ++ l₂) = numRots l₁ + numRots l₂
  | [], _ => by simp [numRots]
  | .cliff _ :: l, l₂ => by simp [numRots, numRots_append l l₂]
  | .rz _ _ :: l, l₂ => by simp [numRots, numRots_append l l₂]; omega

/-- The axes of the rotations after the `j`th, seen from the input of `D`, are the axes of the
part of `D` after the `j`th rotation, translated through `wp_{Γ_j}`: `A_{j+1+k} = wp_{Γ_j}(B_k)`. -/
theorem axis_shift {D p r : List (Gate n)} {j : ℕ} {qj : Fin n} {θj : ℝ}
    (hsplit : splitRot D j = some (p, qj, θj, r)) (k : ℕ) :
    axis D (j + 1 + k) = (axis r k).bind (alphaRel (ofList p)).cons.wp := by
  obtain ⟨rfl, rfl⟩ := splitRot_spec hsplit
  rw [axis, Nat.add_assoc, splitRot_append_right, Nat.add_comm 1 k]
  simp only [splitRot, Option.map_map, Option.bind_map, axis]
  cases splitRot r k with
  | none => rfl
  | some x =>
    simp only [Option.bind_some, Function.comp]
    rw [show p ++ Gate.rz qj θj :: x.1 = p ++ (Gate.rz qj θj :: x.1) from rfl, wp_ofList_append,
      wp_ofList_rz_cons]

/-! ## Algorithm 1 -/

/-- **Paper, Algorithm 1, line 3**: the rotation `R_q(θ)` *folds* into the `j`th rotation of `D`
(counted from `0`) with sign `s ∈ {1, -1}`, where `σ = (Γ, 𝒜)` is the abstraction of `D`:
`wp_Γ(Z_q) = s A_j` and `[A_j, A_l] = 0` for all later rotations `l > j` of `D`. The
commutator is `PP.Comm` (equivalently, of the matrices: `PP.comm_iff`). It never holds for
`σ = ⊤`, which the algorithm never produces. -/
def Folds (σ : Dom n) (D : List (Gate n)) (q : Fin n) (j : ℕ) (s : ℤ) : Prop :=
  match σ with
  | .top => False
  | .rel σ => (s = 1 ∨ s = -1) ∧ ∃ Aj, axis D j = some Aj ∧
      σ.cons.wp (Zq q) = some (signed s Aj) ∧
      ∀ l < numRots D, j < l → ∀ Al, axis D l = some Al → PP.Comm Aj Al

instance (σ : Dom n) (D : List (Gate n)) (q : Fin n) (j : ℕ) (s : ℤ) :
    Decidable (Folds σ D q j s) :=
  match σ with
  | .top => isFalse id
  | .rel σ =>
    match h : axis D j with
    | none => isFalse fun ⟨_, A, hA, _⟩ => by rw [h] at hA; cases hA
    | some Aj => decidable_of_iff ((s = 1 ∨ s = -1) ∧ σ.cons.wp (Zq q) = some (signed s Aj) ∧
        ∀ l < numRots D, j < l → ∀ Al ∈ axis D l, PP.Comm Aj Al)
      ⟨fun ⟨hs, hw, hc⟩ => ⟨hs, Aj, h, hw, hc⟩,
        fun ⟨hs, A, hA, hw, hc⟩ => by
          rw [h, Option.some.injEq] at hA; subst hA; exact ⟨hs, hw, hc⟩⟩

/-- A choice of `j` (and `s`) on line 3, given `σ`, `D` and the rotation `R_q(θ)`; `none` means
"no fold". -/
abbrev Choice (n : ℕ) := Dom n → List (Gate n) → Fin n → ℝ → Option (ℕ × ℤ)

/-- A choice is *valid* when every `j` and `s` it returns satisfy the condition of line 3. -/
def Choice.Valid (ch : Choice n) : Prop :=
  ∀ σ D q θ j s, ch σ D q θ = some (j, s) → Folds σ D q j s

/-- One iteration of the loop of Algorithm 1 on the state `(D, σ)`: if `g = R_q(θ)` folds into
the `j`th rotation of `D` with sign `s` (as chosen by `ch`), add `sθ` to `θ_j`; otherwise
`D ← D; g` and `σ ← α(g) ∘ σ`. -/
def step (ch : Choice n) (st : List (Gate n) × Dom n) : Gate n → List (Gate n) × Dom n
  | .rz q θ =>
    match ch st.2 st.1 q θ with
    | some (j, s) => (addAngle st.1 j (s * θ), st.2)
    | none => (st.1 ++ [.rz q θ], alpha (.rz q θ) ⊚ st.2)
  | .cliff G => (st.1 ++ [.cliff G], alpha (.gate G) ⊚ st.2)

/-- The loop of Algorithm 1, from `D = I` (the empty list) and `σ = α(I)`. -/
def run (ch : Choice n) (gs : List (Gate n)) : List (Gate n) × Dom n :=
  gs.foldl (step ch) ([], alpha (gate (CliffordGate.one n)))

/-- **Paper, Algorithm 1 (Phase folding)**, for the choice `ch` of `j` on line 3: run the loop
on the gates `g₁, …, g_m` of `C` and return `D`. -/
def phaseFold (ch : Choice n) (C : Circuit n) : Circuit n := ofList (run ch C.gates).1

/-! ## The condition on line 3 is the hypothesis of Theorem 5.1 -/

/-- **Paper, Section 5: "The condition on line 3 is the hypothesis of Theorem 5.1 for the part
of `D` after its `j`th rotation, with both Pauli strings expressed at the input of `D` through
`wp_Γ`."** Let `D = pre; R_{q_j}(θ_j); post`. Then `R_q(θ)` folds into the `j`th rotation of `D`
with sign `s` exactly when `Z_q ∈ 𝓕(σ')` and `wp_{Γ'}(Z_q) = s Z_{q_j}` for
`σ' = (Γ', 𝒜') = α(R_{q_j}(θ_j); post)` (here `post` is the circuit `I; post`). -/
theorem folds_iff {D p r : List (Gate n)} {j : ℕ} {qj : Fin n} {θj : ℝ}
    (hsplit : splitRot D j = some (p, qj, θj, r)) (q : Fin n) (s : ℤ) (hs : s = 1 ∨ s = -1) :
    Folds (alpha (ofList D)) D q j s ↔
      Zq q ∈ (alphaRel (rz qj θj ⨾ ofList r)).F ∧
        (alphaRel (rz qj θj ⨾ ofList r)).cons.wp (Zq q) = some (signed s (Zq qj)) := by
  obtain ⟨hD, hj⟩ := splitRot_spec hsplit
  have hax' : (alphaRel (rz qj θj ⨾ ofList r)).axes = Zq qj :: (alphaRel (ofList r)).axes := by
    rw [alphaRel_seq]; simp
  have hwp' : ∀ P, (alphaRel (rz qj θj ⨾ ofList r)).cons.wp P = (alphaRel (ofList r)).cons.wp P :=
    wp_alpha_congr (by simp [cliff])
  have hwpD : ∀ P, (alphaRel (ofList D)).cons.wp P =
      ((alphaRel (ofList r)).cons.wp P).bind (alphaRel (ofList p)).cons.wp := by
    intro P; rw [hD, wp_ofList_append, wp_ofList_rz_cons]
  obtain ⟨Aj, hAj, -⟩ := axis_of_splitRot hsplit
  have hAj' : (alphaRel (ofList p)).cons.wp (Zq qj) = some Aj := by
    rw [axis, hsplit] at hAj; exact hAj
  obtain ⟨B', hB', -⟩ := wp_alpha (ofList r) (Zq q)
  have key1 : (alphaRel (ofList D)).cons.wp (Zq q) = some (signed s Aj) ↔
      B' = signed s (Zq qj) := by
    rw [hwpD, hB', Option.bind_some]
    exact ⟨fun h => wp_alpha_inj _ h (wp_alpha_signed _ hs hAj'),
      fun h => h ▸ wp_alpha_signed _ hs hAj'⟩
  have hnum : numRots D = j + 1 + numRots r := by
    rw [hD, numRots_append, hj]; simp [numRots]; omega
  rw [alpha_eq_rel, hwp', hB']
  simp only [Folds, Option.some.injEq]
  constructor
  · rintro ⟨-, A, hA, hw, hc⟩
    rw [hAj] at hA; cases hA
    have hB'eq := key1.1 hw
    refine ⟨⟨B', by rw [hwp', hB'], ?_⟩, hB'eq⟩
    intro A hA
    rw [hax'] at hA
    rcases List.mem_cons.1 hA with rfl | hA
    · rw [hB'eq, comm_signed_iff]; rfl
    · obtain ⟨k, hk⟩ := (mem_axes_ofList r A).1 hA
      have hklt : k < numRots r := (axis_isSome r k).1 (by rw [hk]; rfl)
      obtain ⟨Al, hAl, -⟩ := wp_alpha (ofList p) A
      have hc' := hc (j + 1 + k) (by omega) (by omega) Al
        (by rw [axis_shift hsplit, hk, Option.bind_some, hAl])
      rw [hB'eq, comm_signed_iff, ← comm_wp_iff_alpha (ofList p) hAj' hAl]
      exact hc'
  · rintro ⟨⟨B, hB, hc⟩, hw⟩
    rw [hwp', hB'] at hB; cases hB
    refine ⟨hs, Aj, hAj, key1.2 hw, ?_⟩
    intro l hl hjl Al hAl
    obtain ⟨k, rfl⟩ : ∃ k, l = j + 1 + k := ⟨l - (j + 1), by omega⟩
    rw [axis_shift hsplit] at hAl
    obtain ⟨Bk, hBk, hAlBk⟩ := Option.bind_eq_some_iff.1 hAl
    have hcB := hc Bk (by rw [hax']; exact List.mem_cons_of_mem _ ((mem_axes_ofList r Bk).2 ⟨k, hBk⟩))
    rw [hw, comm_signed_iff] at hcB
    exact (comm_wp_iff_alpha (ofList p) hAj' hAlBk).2 hcB

/-! ## Theorem 5.2 -/

/-- A fold changes the angle of the `j`th rotation only, so it keeps the abstraction, and by
Theorem 5.1 (through `folds_iff`) it multiplies the semantics by the folded rotation. -/
theorem fold_step {D : List (Gate n)} {q : Fin n} {j : ℕ} {s : ℤ}
    (h : Folds (alpha (ofList D)) D q j s) (θ : ℝ) :
    alpha (ofList (addAngle D j (s * θ))) = alpha (ofList D) ∧
      sem (ofList (addAngle D j (s * θ))) = rzMat q θ * sem (ofList D) := by
  have h' := h
  rw [alpha_eq_rel] at h'
  obtain ⟨hs, Aj, hAj, -⟩ := h'
  obtain ⟨⟨p, qj, θj, r⟩, hsplit, -⟩ := Option.bind_eq_some_iff.1 hAj
  obtain ⟨hD, -⟩ := splitRot_spec hsplit
  obtain ⟨hF, hwp⟩ := (folds_iff hsplit q s hs).1 h
  have hadd : addAngle D j (s * θ) = p ++ .rz qj (θj + s * θ) :: r := by
    rw [addAngle, hsplit]
  rw [hadd, hD]
  refine ⟨?_, ?_⟩
  · simp only [ofList, ofListFrom_append, ofListFrom_cons]
    apply alpha_ofListFrom_congr; rfl
  · have e := fold_sem (ofList r) qj q s hs θj _ (alpha_eq_rel _) hF hwp θ
    simp only [sem_seq, sem_rz] at e
    rw [sem_ofList_append, sem_ofList_append, sem_ofList_cons, sem_ofList_cons]
    simp only [Gate.toCircuit, sem_rz]
    rw [← e]; simp only [Matrix.mul_assoc]

/-- The loop keeps `σ = α(D)` and multiplies `⟦D⟧` by `⟦g⟧`. -/
theorem step_spec {ch : Choice n} (hch : ch.Valid) {D : List (Gate n)} {σ : Dom n}
    (hσ : σ = alpha (ofList D)) (g : Gate n) :
    (step ch (D, σ) g).2 = alpha (ofList (step ch (D, σ) g).1) ∧
      sem (ofList (step ch (D, σ) g).1) = sem g.toCircuit * sem (ofList D) := by
  have hcat : ∀ g : Gate n, alpha g.toCircuit ⊚ σ = alpha (ofList (D ++ [g])) ∧
      sem (ofList (D ++ [g])) = sem g.toCircuit * sem (ofList D) := fun g =>
    ⟨by rw [ofList_concat, alpha_seq, hσ], by rw [ofList_concat, sem_seq]⟩
  cases g with
  | cliff G => exact hcat (.cliff G)
  | rz q θ =>
    simp only [step]
    cases hc : ch σ D q θ with
    | none => exact hcat (.rz q θ)
    | some js =>
      obtain ⟨j, s⟩ := js
      have hf := hch σ D q θ j s hc
      rw [hσ] at hf
      obtain ⟨h1, h2⟩ := fold_step hf θ
      exact ⟨by rw [h1, hσ], by rw [h2]; rfl⟩

theorem run_spec {ch : Choice n} (hch : ch.Valid) (gs : List (Gate n)) :
    ∀ (D : List (Gate n)) (σ : Dom n), σ = alpha (ofList D) →
      (gs.foldl (step ch) (D, σ)).2 = alpha (ofList (gs.foldl (step ch) (D, σ)).1) ∧
        sem (ofList (gs.foldl (step ch) (D, σ)).1) = sem (ofList gs) * sem (ofList D) := by
  induction gs with
  | nil => intro D σ hσ; simp [hσ, ofList, CliffordGate.one]
  | cons g gs ih =>
    intro D σ hσ
    obtain ⟨h1, h2⟩ := step_spec hch hσ g
    rw [List.foldl_cons]
    obtain ⟨h3, h4⟩ := ih _ _ h1
    refine ⟨h3, ?_⟩
    rw [h4, h2, sem_ofList_cons, Matrix.mul_assoc]

/-- The invariant of Algorithm 1: `σ = α(D)`. -/
theorem run_alpha {ch : Choice n} (hch : ch.Valid) (gs : List (Gate n)) :
    (run ch gs).2 = alpha (ofList (run ch gs).1) :=
  (run_spec hch gs [] (alpha (gate (CliffordGate.one n))) rfl).1

/-- Algorithm 1 preserves the semantics exactly. -/
theorem phaseFold_sem {ch : Choice n} (hch : ch.Valid) (C : Circuit n) :
    sem (phaseFold ch C) = sem C := by
  rw [phaseFold, run, (run_spec hch C.gates [] (alpha (gate (CliffordGate.one n))) rfl).2,
    sem_ofList_gates]
  simp [ofList, CliffordGate.one]

/-- **Paper, Theorem 5.2 (Correctness).** For every circuit `C` and every valid choice of `j`
on line 3, Algorithm 1 returns a circuit `D ≡ C`. -/
theorem phaseFold_correct {ch : Choice n} (hch : ch.Valid) (C : Circuit n) :
    (phaseFold ch C).Equiv C :=
  EquivP.of_eq (phaseFold_sem hch C)

/-! ## A valid choice: the latest rotation that folds -/

/-- The candidates `(j, s)`, latest rotation first. -/
def candidates (D : List (Gate n)) : List (ℕ × ℤ) :=
  (List.range (numRots D)).reverse.flatMap fun j => [(j, 1), (j, -1)]

/-- A choice function: the latest rotation of `D` into which `R_q(θ)` folds, with `s = 1` tried
first. -/
def latest : Choice n := fun σ D q _ =>
  (candidates D).find? fun js => decide (Folds σ D q js.1 js.2)

theorem latest_valid : (latest : Choice n).Valid := by
  intro σ D q θ j s h
  simpa using List.find?_some h

/-- `latest` folds whenever some fold is possible. -/
theorem latest_ne_none {σ : Dom n} {D : List (Gate n)} {q : Fin n} {θ : ℝ} {j : ℕ} {s : ℤ}
    (h : Folds σ D q j s) : latest σ D q θ ≠ none := by
  intro hn
  rw [latest, List.find?_eq_none] at hn
  have hmem : (j, s) ∈ candidates D := by
    cases σ with
    | top => exact h.elim
    | rel σ =>
      obtain ⟨hs, A, hA, -⟩ := h
      have hj : j < numRots D := (axis_isSome D j).1 (by rw [hA]; rfl)
      simp only [candidates, List.mem_flatMap, List.mem_reverse, List.mem_range]
      exact ⟨j, hj, by rcases hs with rfl | rfl <;> simp⟩
  exact hn _ hmem (by simpa using h)

theorem not_folds_nil (σ : Dom n) (q : Fin n) (j : ℕ) (s : ℤ) : ¬ Folds σ [] q j s := by
  cases σ with
  | top => exact id
  | rel σ => rintro ⟨-, A, hA, -⟩; cases hA

/-! ## The worked example: Figure 1 and Section 2.3 -/

/-- `S_q` is `R_Z^q(π/2)` up to the global phase `e^{iπ/4}`. -/
theorem gS_mat (q : Fin n) :
    (gS q).mat = Complex.exp ((Real.pi / 4 : ℝ) * Complex.I) • rzMat q (Real.pi / 2) := by
  show embed1 n (diag2 1 (ep (1 / 2))) q.val = _
  rw [rzMat, ← embed1_smul _ _ _ q.isLt]
  congr 1
  funext o i
  cases o <;> cases i <;> simp [diag2, ep, ← Complex.exp_add] <;> (ring_nf; try simp)

namespace Figure1

/-- `CX_{1,2}` (Lean qubits `0, 1`). -/
abbrev cx : CliffordGate 2 := gCX 0 1 (by decide)

/-- `CX_{2,1}` (Lean qubits `1, 0`). -/
abbrev cx' : CliffordGate 2 := gCX 1 0 (by decide)

/-- The Clifford gates between the two `T` gates of Figure 1: `CX_{1,2}; H₁; H₂; CX_{2,1}; H₂`. -/
def mid : Circuit 2 := gate cx ⨾ gate (gH 0) ⨾ gate (gH 1) ⨾ gate cx' ⨾ gate (gH 1)

/-- Figure 1(a), "Before": `T₂; CX_{1,2}; H₁; H₂; CX_{2,1}; H₂; T₂` (Lean qubits `0, 1`). -/
def before : Circuit 2 :=
  T 1 ⨾ gate cx ⨾ gate (gH 0) ⨾ gate (gH 1) ⨾ gate cx' ⨾ gate (gH 1) ⨾ T 1

/-- Figure 1(b), "After": `S₂; CX_{1,2}; H₁; H₂; CX_{2,1}; H₂`. -/
def after : Circuit 2 :=
  gate (gS 1) ⨾ gate cx ⨾ gate (gH 0) ⨾ gate (gH 1) ⨾ gate cx' ⨾ gate (gH 1)

/-- Section 2.3: going backwards, the gates between the two `T` gates turn `Z₂` into `X₂`, then
`X₁X₂`, then `Z₁Z₂`, and finally back into `Z₂`. -/
theorem chain :
    (alphaRel (gate (gH (1 : Fin 2)))).cons.wp (Zq 1) = some (Xq 1) ∧
    (alphaRel (gate cx' ⨾ gate (gH 1))).cons.wp (Zq 1) = some ((Xq 0).mul (Xq 1)) ∧
    (alphaRel (gate (gH 0) ⨾ gate (gH 1) ⨾ gate cx' ⨾ gate (gH 1))).cons.wp (Zq 1) =
      some ((Zq 0).mul (Zq 1)) ∧
    (alphaRel mid).cons.wp (Zq 1) = some (Zq 1) := by
  decide

/-- Section 2.3: at cut 1 the abstraction has the axis `Z₂`; at cut 2 (after the Clifford
gates), `wp(Z₂) = Z₂`, no other axis has been added, and `Z₂` is tracked. These are the
hypotheses of Theorem 5.1 with `s = 1`. -/
theorem cuts :
    alpha (T (1 : Fin 2)) = .rel ⟨Cons.id 2, [Zq 1]⟩ ∧
    (alphaRel (T 1 ⨾ mid)).cons.wp (Zq 1) = some (Zq 1) ∧
    (alphaRel (T 1 ⨾ mid)).axes = [Zq 1] ∧
    Zq 1 ∈ (alphaRel (T 1 ⨾ mid)).F := by
  decide

/-- Figure 1, via Theorem 5.1: the second `T` merges into the first. -/
theorem before_eq : sem before = sem (rz 1 (Real.pi / 4 + Real.pi / 4) ⨾ mid) := by
  have h := fold_sem mid 1 1 1 (Or.inl rfl) (Real.pi / 4) _ (alpha_eq_rel _)
    cuts.2.2.2 (by rw [cuts.2.1]; rfl) (Real.pi / 4)
  simp only [Int.cast_one, one_mul] at h
  rw [← h]
  simp only [before, mid, sem_seq, sem_gate, sem_rz, Matrix.mul_assoc]

theorem after_equiv : (rz 1 (Real.pi / 4 + Real.pi / 4) ⨾ mid).Equiv after := by
  refine ⟨-(Real.pi / 4), ?_⟩
  simp only [after, mid, sem_seq, sem_gate, sem_rz, gS_mat, Matrix.mul_smul, smul_smul,
    ← Complex.exp_add, Matrix.mul_assoc]
  rw [show Real.pi / 4 + Real.pi / 4 = Real.pi / 2 by ring]
  push_cast; simp

/-- **Paper, Figure 1 / Section 2.3**: `T₂; CX_{1,2}; H₁; H₂; CX_{2,1}; H₂; T₂ ≡
S₂; CX_{1,2}; H₁; H₂; CX_{2,1}; H₂` (up to a global phase). -/
theorem figure1 : before.Equiv after :=
  (EquivP.of_eq before_eq).trans after_equiv

/-- `D` when Algorithm 1 reaches the second `T` gate: the first six gates of Figure 1(a). -/
def D₆ : List (Gate 2) :=
  [.rz 1 (Real.pi / 4), .cliff cx, .cliff (gH 0), .cliff (gH 1), .cliff cx', .cliff (gH 1)]

/-- The output of Algorithm 1 on Figure 1(a): the first rotation now has angle `π/4 + π/4`. -/
def folded : List (Gate 2) :=
  [.rz 1 (Real.pi / 4 + Real.pi / 4), .cliff cx, .cliff (gH 0), .cliff (gH 1), .cliff cx',
    .cliff (gH 1)]

theorem before_gates : before.gates = D₆ ++ [.rz 1 (Real.pi / 4)] := rfl

/-- **Paper, Section 5: "On Figure 1(a), the second `T` gate folds into the first with
`s = 1`."** When Algorithm 1 reaches the second `T₂`, it folds into the first rotation of `D` (Lean
rotation `0`) with `s = 1`, and into no other rotation with any other sign. -/
theorem figure1_folds :
    Folds (alpha (ofList D₆)) D₆ 1 0 1 ∧
      ∀ j s, Folds (alpha (ofList D₆)) D₆ 1 j s → j = 0 ∧ s = 1 := by
  refine ⟨by decide, fun j s h => ?_⟩
  have h' := h
  rw [alpha_eq_rel] at h'
  obtain ⟨hs, A, hA, -⟩ := h'
  have hj : j < numRots D₆ := (axis_isSome D₆ j).1 (by rw [hA]; rfl)
  have hj0 : j = 0 := by simp [numRots, D₆] at hj; omega
  subst hj0
  rcases hs with rfl | rfl
  · exact ⟨rfl, rfl⟩
  · exact absurd h (by decide)

/-- Algorithm 1 reaches the second `T` gate with `D = D₆` and `σ = α(D₆)`, whatever the
(valid) choice function: the first `T` cannot fold, since `D` has no rotation yet. -/
theorem run_D₆ {ch : Choice 2} (hch : ch.Valid) : run ch D₆ = (D₆, alpha (ofList D₆)) := by
  have h1 : ch (alpha (gate (CliffordGate.one 2))) [] 1 (Real.pi / 4) = none := by
    cases h : ch (alpha (gate (CliffordGate.one 2))) [] 1 (Real.pi / 4) with
    | none => rfl
    | some js => exact absurd (hch _ _ _ _ _ _ h) (not_folds_nil _ _ _ _)
  simp only [run, D₆, List.foldl, step, h1]
  rfl

/-- **Paper, Section 5 and Figure 1**: Algorithm 1, with any valid choice that folds when it
can (e.g. `latest`), turns Figure 1(a) into `R_Z^2(π/4 + π/4); CX_{1,2}; H₁; H₂; CX_{2,1}; H₂`,
that is, Figure 1(b) up to a global phase. -/
theorem figure1_phaseFold {ch : Choice 2} (hch : ch.Valid)
    (hfold : ch (alpha (ofList D₆)) D₆ 1 (Real.pi / 4) ≠ none) :
    phaseFold ch before = ofList folded ∧ (phaseFold ch before).Equiv after := by
  obtain ⟨⟨j, s⟩, h7⟩ := Option.ne_none_iff_exists'.1 hfold
  obtain ⟨rfl, rfl⟩ := figure1_folds.2 j s (hch _ _ _ _ _ _ h7)
  have hrun : run ch before.gates = (folded, alpha (ofList D₆)) := by
    rw [before_gates, run, List.foldl_append, ← run, run_D₆ hch, List.foldl_cons, List.foldl_nil,
      step]
    simp only [h7]
    simp [addAngle, splitRot, D₆, folded]
  have hout : phaseFold ch before = ofList folded := by rw [phaseFold, hrun]
  exact ⟨hout, (phaseFold_correct hch before).trans figure1⟩

/-- `latest` folds the second `T` gate of Figure 1(a) into the first. -/
theorem figure1_latest :
    phaseFold latest before = ofList folded ∧ (phaseFold latest before).Equiv after :=
  figure1_phaseFold latest_valid (latest_ne_none figure1_folds.1)

end Figure1

end

end PauliFold.Paper
