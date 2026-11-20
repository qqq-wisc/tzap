import PauliFold.Join

/-!
# A bounded disjunctive domain

The per-row join of `Join.lean` forgets every row on which two branches disagree. After
`while ★ { swap a,b }` it forgets all four rows, although zero and one iteration agree on
`Z_a Z_b`. A *disjunctive* domain keeps the branches apart: an abstract value is a list of
partial tableaux, read as "the program is described by one of them".

- **Concretization.** The union of the disjuncts' concretizations. A fact holds for the
  program when it holds in every disjunct.
- **Normalization.** Disjuncts with the same rows are merged, keeping the rows and taking
  the union of the axes. This lets loops with rotations in the body stabilize.
- **Join.** Concatenate and normalize. If more than `bound` disjuncts remain, collapse them
  with the per-row join into one.
- **Transformers.** Applied disjunct by disjunct. If one disjunct reaches `⊤`, the whole
  value is `⊤`.
- **Loops.** As in `Join.lean`: iterate `S ← S ⊔ body♯(S)` until `body♯(S) ⊑ S`.

Soundness (`absRunD_sound`, `programD_sound`) reuses the per-disjunct lemmas of
`Join.lean`. Examples checked by `decide`:

- `while ★ { swap a,b }` maps `Z_a Z_b` to itself (`loopSwap_ZZ`), which the per-row join
  cannot show;
- `while ★ { swap a,b; swap b,c }` does not preserve `Z_b` (the value of `b` cycles).
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-- The maximum number of disjuncts before collapsing. -/
def bound : ℕ := 8

namespace PTab

theorem leP_refl (X : PTab n) : LeP X X := ⟨fun _ => Or.inr rfl, fun _ => Or.inr rfl, fun _ h => h⟩

theorem leP_trans {X Y Z : PTab n} (h₁ : LeP X Y) (h₂ : LeP Y Z) : LeP X Z := by
  obtain ⟨hx1, hz1, ha1⟩ := h₁
  obtain ⟨hx2, hz2, ha2⟩ := h₂
  refine ⟨fun r => ?_, fun r => ?_, fun A hA => ha2 A (ha1 A hA)⟩
  · rcases hx2 r with h | h
    · exact Or.inl h
    · rcases hx1 r with h' | h'
      · exact Or.inl (h.trans h')
      · exact Or.inr (h.trans h')
  · rcases hz2 r with h | h
    · exact Or.inl h
    · rcases hz1 r with h' | h'
      · exact Or.inl (h.trans h')
      · exact Or.inr (h.trans h')

theorem leP_joinP_left (X Y : PTab n) : LeP X (joinP X Y) := le_join_left (some X) (some Y)
theorem leP_joinP_right (X Y : PTab n) : LeP Y (joinP X Y) := le_join_right (some X) (some Y)

end PTab

/-! ## Normalizing a list of disjuncts -/

/-- Add a disjunct, merging it with one that has the same rows. -/
def insertD (σ : PTab n) : List (PTab n) → List (PTab n)
  | [] => [σ]
  | τ :: L =>
    if σ.rx = τ.rx ∧ σ.rz = τ.rz then ⟨τ.rx, τ.rz, τ.axes ++ σ.axes⟩ :: L else τ :: insertD σ L

theorem insertD_new (σ : PTab n) :
    ∀ L : List (PTab n), ∃ τ ∈ insertD σ L, LeP σ τ
  | [] => ⟨σ, List.mem_singleton_self _, PTab.leP_refl σ⟩
  | τ :: L => by
    unfold insertD
    split_ifs with h
    · exact ⟨_, List.mem_cons_self .., fun r => Or.inr (by rw [h.1]), fun r => Or.inr (by rw [h.2]),
        fun A hA => List.mem_append_right _ hA⟩
    · obtain ⟨τ', hτ', hle⟩ := insertD_new σ L
      exact ⟨τ', List.mem_cons_of_mem _ hτ', hle⟩

theorem insertD_old (σ : PTab n) :
    ∀ (L : List (PTab n)) (τ : PTab n), τ ∈ L → ∃ τ' ∈ insertD σ L, LeP τ τ'
  | [], τ, h => absurd h List.not_mem_nil
  | ρ :: L, τ, h => by
    unfold insertD
    split_ifs with hs
    · rcases List.mem_cons.1 h with rfl | h
      · exact ⟨_, List.mem_cons_self .., fun r => Or.inr rfl, fun r => Or.inr rfl,
          fun A hA => List.mem_append_left _ hA⟩
      · exact ⟨τ, List.mem_cons_of_mem _ h, PTab.leP_refl τ⟩
    · rcases List.mem_cons.1 h with rfl | h
      · exact ⟨τ, List.mem_cons_self .., PTab.leP_refl τ⟩
      · obtain ⟨τ', hτ', hle⟩ := insertD_old σ L τ h
        exact ⟨τ', List.mem_cons_of_mem _ hτ', hle⟩

/-- Merge all disjuncts with equal rows. -/
def mergeAll (L : List (PTab n)) : List (PTab n) := L.foldr insertD []

theorem mergeAll_le : ∀ (L : List (PTab n)) (σ : PTab n), σ ∈ L → ∃ τ ∈ mergeAll L, LeP σ τ
  | [], σ, h => absurd h List.not_mem_nil
  | ρ :: L, σ, h => by
    show ∃ τ ∈ insertD ρ (mergeAll L), LeP σ τ
    rcases List.mem_cons.1 h with rfl | h
    · exact insertD_new σ (mergeAll L)
    · obtain ⟨τ, hτ, hle⟩ := mergeAll_le L σ h
      obtain ⟨τ', hτ', hle'⟩ := insertD_old ρ (mergeAll L) τ hτ
      exact ⟨τ', hτ', PTab.leP_trans hle hle'⟩

theorem foldl_joinP_le (σ : PTab n) :
    ∀ (L : List (PTab n)) (τ : PTab n), (τ = σ ∨ τ ∈ L) → LeP τ (L.foldl joinP σ)
  | [], τ, h => by
    rcases h with rfl | h
    · exact PTab.leP_refl _
    · exact absurd h List.not_mem_nil
  | ρ :: L, τ, h => by
    show LeP τ (L.foldl joinP (joinP σ ρ))
    rcases h with rfl | h
    · exact PTab.leP_trans (PTab.leP_joinP_left τ ρ) (foldl_joinP_le _ L _ (Or.inl rfl))
    · rcases List.mem_cons.1 h with rfl | h
      · exact PTab.leP_trans (PTab.leP_joinP_right σ τ) (foldl_joinP_le _ L _ (Or.inl rfl))
      · exact foldl_joinP_le _ L τ (Or.inr h)

/-- Merge equal rows; if more than `bound` disjuncts remain, collapse them with the per-row
join. -/
def normalize (L : List (PTab n)) : List (PTab n) :=
  match mergeAll L with
  | [] => []
  | σ :: M => if (σ :: M).length ≤ bound then σ :: M else [M.foldl joinP σ]

theorem normalize_le (L : List (PTab n)) (σ : PTab n) (h : σ ∈ L) :
    ∃ τ ∈ normalize L, LeP σ τ := by
  obtain ⟨τ, hτ, hle⟩ := mergeAll_le L σ h
  unfold normalize
  split
  · rename_i hm; rw [hm] at hτ; exact absurd hτ List.not_mem_nil
  · rename_i ρ M hm
    rw [hm] at hτ
    split_ifs
    · exact ⟨τ, hτ, hle⟩
    · refine ⟨_, List.mem_singleton_self _, PTab.leP_trans hle ?_⟩
      rcases List.mem_cons.1 hτ with rfl | h
      · exact foldl_joinP_le _ M _ (Or.inl rfl)
      · exact foldl_joinP_le _ M _ (Or.inr h)

/-! ## The domain -/

/-- A list of disjuncts, or `⊤` (`none`). -/
abbrev DAbs (n : ℕ) := Option (List (PTab n))

/-- Sound for `U`: `⊤`, or some disjunct is sound for `U`. -/
def SoundD (U : Density n) : DAbs n → Prop
  | none => True
  | some L => ∃ σ ∈ L, Sound U (some σ)

/-- Concretization: the union of the disjuncts'. -/
def gammaD : DAbs n → Set (Density n)
  | none => Set.univ
  | some L => {U | ∃ σ ∈ L, U ∈ gammaP σ}

theorem gammaD_of_sound {U : Density n} {a : DAbs n} (h : SoundD U a) : U ∈ gammaD a := by
  cases a with
  | none => trivial
  | some L => obtain ⟨σ, hσ, hs⟩ := h; exact ⟨σ, hσ, gammaA_of_sound hs⟩

/-- Every disjunct is below some disjunct of the other side. -/
def LeD : DAbs n → DAbs n → Prop
  | _, none => True
  | none, some _ => False
  | some L, some M => ∀ σ ∈ L, ∃ τ ∈ M, LeP σ τ

instance (a b : DAbs n) : Decidable (LeD a b) := by
  cases a <;> cases b <;> unfold LeD <;> infer_instance

theorem soundD_of_le {U : Density n} {a b : DAbs n} (h : SoundD U a) (hle : LeD a b) :
    SoundD U b := by
  cases b with
  | none => trivial
  | some M =>
    cases a with
    | none => exact hle.elim
    | some L =>
      obtain ⟨σ, hσ, hs⟩ := h
      obtain ⟨τ, hτ, hle'⟩ := hle σ hσ
      exact ⟨τ, hτ, sound_of_le hs hle'⟩

/-- **The join**: union of the disjuncts, normalized. -/
def joinD : DAbs n → DAbs n → DAbs n
  | some L, some M => some (normalize (L ++ M))
  | _, _ => none

theorem le_joinD_left (a b : DAbs n) : LeD a (joinD a b) := by
  cases a with
  | none => cases b <;> trivial
  | some L =>
    cases b with
    | none => trivial
    | some M => exact fun σ hσ => normalize_le _ σ (List.mem_append_left _ hσ)

theorem le_joinD_right (a b : DAbs n) : LeD b (joinD a b) := by
  cases a with
  | none => cases b <;> trivial
  | some L =>
    cases b with
    | none => trivial
    | some M => exact fun σ hσ => normalize_le _ σ (List.mem_append_right _ hσ)

/-- Apply a gate to every disjunct; `⊤` if any disjunct reaches `⊤`. -/
def stepList (g : CTGate n) : List (PTab n) → Option (List (PTab n))
  | [] => some []
  | σ :: L => (PTab.stepP g σ).bind fun σ' => (stepList g L).map (σ' :: ·)

theorem stepList_mem (g : CTGate n) :
    ∀ {L L' : List (PTab n)}, stepList g L = some L' → ∀ σ ∈ L,
      ∃ σ', PTab.stepP g σ = some σ' ∧ σ' ∈ L'
  | [], _, _, σ, h => absurd h List.not_mem_nil
  | ρ :: L, L', hs, σ, hσ => by
    simp only [stepList, Option.bind_eq_some_iff, Option.map_eq_some_iff] at hs
    obtain ⟨ρ', hρ', M, hM, rfl⟩ := hs
    rcases List.mem_cons.1 hσ with rfl | hσ
    · exact ⟨ρ', hρ', List.mem_cons_self ..⟩
    · obtain ⟨σ', h1, h2⟩ := stepList_mem g hM σ hσ
      exact ⟨σ', h1, List.mem_cons_of_mem _ h2⟩

def stepD (g : CTGate n) (a : DAbs n) : DAbs n := a.bind (stepList g)

theorem stepD_sound {U : Density n} {a : DAbs n} (h : SoundD U a) (g : CTGate n) :
    SoundD (gateUnitary n g.toGate * U) (stepD g a) := by
  cases a with
  | none => trivial
  | some L =>
    obtain ⟨σ, hσ, hs⟩ := h
    show SoundD _ (stepList g L)
    cases hL : stepList g L with
    | none => trivial
    | some L' =>
      obtain ⟨σ', h1, h2⟩ := stepList_mem g hL σ hσ
      have := stepA_sound hs g
      simp only [stepA, Option.bind_some, h1] at this
      exact ⟨σ', h2, this⟩

/-! ## Programs -/

def loopFixD (f : DAbs n → DAbs n) : ℕ → DAbs n → DAbs n
  | 0, _ => none
  | k + 1, a =>
    if LeD (f (joinD a (f a))) (joinD a (f a)) then joinD a (f a) else loopFixD f k (joinD a (f a))

theorem loopFixD_spec (f : DAbs n → DAbs n) :
    ∀ (k : ℕ) (a : DAbs n), (∀ U, SoundD U a → SoundD U (loopFixD f k a)) ∧
      LeD (f (loopFixD f k a)) (loopFixD f k a)
  | 0, a => ⟨fun _ _ => trivial, by cases f none <;> trivial⟩
  | k + 1, a => by
    simp only [loopFixD]
    split_ifs with hle
    · exact ⟨fun U h => soundD_of_le h (le_joinD_left _ _), hle⟩
    · obtain ⟨h1, h2⟩ := loopFixD_spec f k (joinD a (f a))
      exact ⟨fun U h => h1 U (soundD_of_le h (le_joinD_left _ _)), h2⟩

/-- The abstract semantics over the disjunctive domain. -/
def absRunD : Prog n → DAbs n → DAbs n
  | .gate g, a => stepD g a
  | .seq p q, a => absRunD q (absRunD p a)
  | .choice p q, a => joinD (absRunD p a) (absRunD q a)
  | .loop p, a => loopFixD (absRunD p) (2 * n + 3 + bound) a

/-- **Soundness of the disjunctive abstract semantics.** -/
theorem absRunD_sound {p : Prog n} {U : Density n} (hU : Exec p U) :
    ∀ {a : DAbs n} {U₀ : Density n}, SoundD U₀ a → SoundD (U * U₀) (absRunD p a) := by
  induction hU with
  | gate g => intro a U₀ h; exact stepD_sound h g
  | seq _ _ ihp ihq =>
    intro a U₀ h
    rw [Matrix.mul_assoc]
    exact ihq (ihp h)
  | choiceL _ ih => intro a U₀ h; exact soundD_of_le (ih h) (le_joinD_left _ _)
  | choiceR _ ih => intro a U₀ h; exact soundD_of_le (ih h) (le_joinD_right _ _)
  | @loopZero p =>
    intro a U₀ h
    rw [Matrix.one_mul]
    exact (loopFixD_spec (absRunD p) _ a).1 _ h
  | @loopStep p U V _ _ ihloop ihbody =>
    intro a U₀ h
    have hI : SoundD (U * U₀) (loopFixD (absRunD p) (2 * n + 3 + bound) a) := ihloop h
    rw [Matrix.mul_assoc]
    exact soundD_of_le (ihbody hI) (loopFixD_spec (absRunD p) _ a).2

theorem soundD_init : SoundD (1 : Density n) (some [PTab.init n]) :=
  ⟨PTab.init n, List.mem_singleton_self _, sound_init⟩

/-- **Programs are soundly abstracted** by the disjunctive domain. -/
theorem programD_sound {p : Prog n} {U : Density n} (hU : Exec p U) :
    U ∈ gammaD (absRunD p (some [PTab.init n])) := by
  have := absRunD_sound hU soundD_init
  rw [Matrix.mul_one] at this
  exact gammaD_of_sound this

/-- A fact that holds in every disjunct of the result holds for every execution. -/
theorem factD_of_run {p : Prog n} {U : Density n} (hU : Exec p U) {L : List (PTab n)}
    (hL : absRunD p (some [PTab.init n]) = some L) {P B : PP n}
    (hf : ∀ σ ∈ L, σ.back P = some B ∧ ∀ A ∈ σ.axes, PP.Comm B A) :
    conj U B.toMatrix = P.toMatrix := by
  have h := programD_sound hU
  rw [hL] at h
  obtain ⟨σ, hσ, hg⟩ := h
  exact hg P B (hf σ hσ).1 (hf σ hσ).2

/-- The check of `factD_of_run`, as a Boolean. -/
def factD (P : PP n) (L : List (PTab n)) : Bool :=
  L.all fun σ => σ.back P = some P && decide (∀ A ∈ σ.axes, PP.Comm P A)

/-! ## Examples -/

/-- The swap, as three CNOTs. -/
def swapP (a b : Fin n) (h : a ≠ b) : Prog n :=
  .seq (.gate (.cx a b h)) (.seq (.gate (.cx b a h.symm)) (.gate (.cx a b h)))

/-- `while ★ { swap a,b }`, the loop of `loop-swap`. -/
def loopSwapSeg : Prog 2 := .loop (swapP 0 1 (by decide))

theorem loopSwapSeg_some : (absRunD loopSwapSeg (some [PTab.init 2])).isSome := by decide

/-- **The swap loop preserves `Z_a Z_b`.** The disjuncts after the loop are the rows of the
identity and of the swap, and both pull `Z_a Z_b` back to itself. So the `T` and `T†` of
`loop-swap`, which act on this parity, merge. -/
theorem loopSwap_ZZ {U : Density 2} (hU : Exec loopSwapSeg U) :
    conj U ((PP.single (0 : Fin 2) .Z).mul (PP.single 1 .Z)).toMatrix =
      ((PP.single (0 : Fin 2) .Z).mul (PP.single 1 .Z)).toMatrix := by
  refine factD_of_run hU (Option.some_get loopSwapSeg_some).symm fun σ hσ => ?_
  revert σ
  decide

/-- The two disjuncts: zero and one iteration. -/
example : ((absRunD loopSwapSeg (some [PTab.init 2])).get loopSwapSeg_some).length = 2 := by
  decide

/-- `while ★ { swap a,b; swap b,c }` cycles the value of `b`, so `Z_b` is not preserved in
every disjunct. -/
def loopCycleSeg : Prog 3 := .loop (.seq (swapP 0 1 (by decide)) (swapP 1 2 (by decide)))

example : ((absRunD loopCycleSeg (some [PTab.init 3])).map fun L =>
    L.all fun σ => σ.back (PP.single 1 .Z) = some (PP.single 1 .Z)) = some false := by
  decide

end

end PauliFold
