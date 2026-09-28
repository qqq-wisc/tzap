import PauliFold.StateFold.Fold

/-!
# A fold StateFold misses at every degree

The circuit

  `t q1; cx q0,q1; h q0; h q1; cx q1,q0; h q1; tdg q1`

has its two rotations on the same Pauli axis: the Pauli domain carries `+Z₁` from the T to
the T† (`cex_pauli`), so the pair folds to the identity. StateFold's analysis ends with the
wires `y₂ ⊕ y₃` and `y₄`, which between them mention every path variable. Neither the HH
nor the ω rule may sum out a variable a wire mentions, so the only trace is the empty one,
and the two predicates `x₁` and `y₄` are neither equal nor complementary
(`cex_not_foldable`), for every degree `d`.
-/

namespace PauliFold.SF

open TzapLean

noncomputable section

variable {n : ℕ}

/-- If every path variable appears in some wire, no reduction applies: every trace is empty. -/
theorem trace_eq_of_exposed {d : ℕ} {A C : State n} {f g f' g' : Poly}
    (hexp : ∀ y ∈ A.temps, ∃ q, ¬ Supp (fun ν => ev ν (A.ket q)) (A.dom \ {y}))
    (htr : Trace d A f g C f' g') : f' = f ∧ g' = g := by
  cases htr with
  | refl => exact ⟨rfl, rfl⟩
  | hh h _ _ _ =>
    obtain ⟨q, hq⟩ := hexp _ h.hy
    exact absurd (h.ket q) hq
  | omega h _ _ _ =>
    obtain ⟨q, hq⟩ := hexp _ h.hy
    exact absurd (h.ket q) hq

/-- A wire whose value flips with `y` (from the all-false valuation) mentions `y`. -/
theorem not_supp_of_flip {p : Poly} {D : Set ℕ} (y : ℕ)
    (h : ev (fun _ => false) p ≠ ev (Function.update (fun _ => false) y true) p) :
    ¬ Supp (fun ν => ev ν p) (D \ {y}) := fun hs =>
  h (hs _ _ fun v hv => by
    rw [Function.update_of_ne (by simpa using hv.2)])

/-- The segment between the two rotations. -/
def cexMid : List (CTGate 2) :=
  [.cx 0 1 (by decide), .h 0, .h 1, .cx 1 0 (by decide), .h 1]

/-- The full circuit, in the shape of `fold_sound`. -/
def cex : List (CTGate 2) := [] ++ .t 1 :: cexMid ++ .tdg 1 :: []

/-- The Pauli domain carries the T's axis `+Z₁` unchanged to the T†. -/
theorem cex_pauli :
    PauliAbs.run cexMid (.one (zString {1} false)) = .one (zString {1} false) := by
  simp only [PauliAbs.run, cexMid, List.foldl, PauliAbs.step]
  congr 1
  simp only [zString, SPauli.conjClifford]
  simp only [SPauli.apply1, SPauli.apply2]
  congr 1
  funext i; fin_cases i <;> decide

theorem cex_temps (y : ℕ) : y ∈ (run cex).temps ↔ y = 2 ∨ y = 3 ∨ y = 4 := by
  simp [run, cex, cexMid, runFrom, step, rot, init]
  omega

theorem cex_ket0 (ν : Val) : ev ν ((run cex).ket 0) = (ν 2 != ν 3) := by
  simp [run, cex, cexMid, runFrom, step, rot, init]

theorem cex_ket1 (ν : Val) : ev ν ((run cex).ket 1) = ν 4 := by
  simp [run, cex, cexMid, runFrom, step, rot, init]

theorem cex_pred₂ (ν : Val) : ev ν (predAfter₂ [] cexMid 1) = ν 4 := by
  simp [predAfter₂, cexMid, runFrom, step, init]

/-- StateFold cannot fold the pair, at any degree. -/
theorem cex_not_foldable (d : ℕ) (s : Bool) {C : State 2} {f' g' : Poly}
    (htr : Trace d (run cex) (predAfter ([] : List (CTGate 2)) 1) (predAfter₂ [] cexMid 1) C f' g') :
    ¬ ∀ ν, ev ν g' = (s != ev ν f') := by
  have hexp : ∀ y ∈ (run cex).temps, ∃ q,
      ¬ Supp (fun ν => ev ν ((run cex).ket q)) ((run cex).dom \ {y}) := by
    intro y hy
    rcases (cex_temps y).1 hy with rfl | rfl | rfl
    · exact ⟨0, not_supp_of_flip _ (by simp [cex_ket0])⟩
    · exact ⟨0, not_supp_of_flip _ (by simp [cex_ket0])⟩
    · exact ⟨1, not_supp_of_flip _ (by simp [cex_ket1])⟩
  obtain ⟨rfl, rfl⟩ := trace_eq_of_exposed hexp htr
  intro heq
  have h₁ := heq fun _ => false
  have h₂ := heq (Function.update (fun _ => false) 4 true)
  simp only [cex_pred₂, predAfter, runFrom, init, ev_var] at h₁ h₂
  cases s <;> simp at h₁ h₂

end

end PauliFold.SF
