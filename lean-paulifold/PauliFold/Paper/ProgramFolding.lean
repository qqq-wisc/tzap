import PauliFold.Paper.PhaseFolding
import PauliFold.Paper.Programs

/-!
# Paper, Section 5.4: phase folding for programs

* Corollary 5.13 (Phase folding for programs): `program_folding`. By soundness, every run of
  `C` maps `±Z_q` to `Z_{q'}` (`maps_of_mem_gamma`), so the two programs have the same runs
  (`program_folding_of_runs`).
* The same merge with the disjunctive refinement of Section 7 (RQ2): `program_folding_disj`.
* "Since `⟦·⟧` is compositional, the merge applies anywhere inside a larger program":
  `psem_seq_congr`, `psem_choice_congr`, `psem_loop_congr`.
* Example 5.11 (Loop-simple) is an instance: `loopSimple_folds`.

Unlike Theorem 4.1, the two sides have the same runs, with no global phase.
Qubits are 0-indexed: the paper's qubit 1 is the Lean qubit `0`.
-/

namespace PauliFold.Paper

open TzapLean Matrix Program

noncomputable section

variable {n : ℕ}

/-- The core of Corollary 5.13: if every run of `C` maps `s Z_q` to `Z_{q'}`, then
`⟦R_q(θ); C; R_{q'}(θ')⟧ = ⟦R_q(θ + sθ'); C⟧`. -/
theorem program_folding_of_runs (C : Program n) (q q' : Fin n) (s : ℤ) (hs : s = 1 ∨ s = -1)
    (hrun : ∀ K ∈ C.psem, K * (signed s (Zq q)).toMatrix = (Zq q').toMatrix * K) (θ θ' : ℝ) :
    (seq (seq (rz q θ) C) (rz q' θ')).psem = (seq (rz q (θ + s * θ')) C).psem := by
  have key : ∀ K ∈ C.psem, rzMat q' θ' * (K * rzMat q θ) = K * rzMat q (θ + s * θ') := by
    intro K hK
    rw [← Matrix.mul_assoc, rz_slide hs (hrun K hK) θ', Matrix.mul_assoc, rzMat_mul, add_comm]
  ext M
  constructor
  · rintro ⟨_, ⟨_, rfl, K, hK, rfl⟩, _, rfl, rfl⟩
    exact ⟨_, rfl, K, hK, key K hK⟩
  · rintro ⟨_, rfl, K, hK, rfl⟩
    exact ⟨_, ⟨_, rfl, K, hK, rfl⟩, _, rfl, (key K hK).symm⟩

/-- The fold condition in one abstract relation: `wp_Γ(Z_{q'}) = s Z_q`, and `Z_q` commutes with
every axis. Every operator of `γ(σ)` then maps `s Z_q` to `Z_{q'}`. -/
theorem maps_of_mem_gamma {σ : ARel n} {q q' : Fin n} {s : ℤ}
    (hwp : σ.cons.wp (Zq q') = some (signed s (Zq q))) (hc : ∀ A ∈ σ.axes, PP.Comm (Zq q) A)
    {K : Density n} (hK : K ∈ σ.gamma) :
    K * (signed s (Zq q)).toMatrix = (Zq q').toMatrix * K :=
  hK (Zq q') _ hwp (fun A hA => (comm_signed_iff _ _).2 (hc A hA))

/-- **Paper, Corollary 5.13 (Phase folding for programs).** Let `C` be a program with
`α(C) = (Γ, 𝒜)`. If `wp_Γ(Z_{q'}) = s Z_q` with `s = ±1`, and `Z_q` commutes with every
`A ∈ 𝒜`, then `⟦R_q(θ); C; R_{q'}(θ')⟧ = ⟦R_q(θ + sθ'); C⟧`. -/
theorem program_folding (C : Program n) (q q' : Fin n) (s : ℤ) (hs : s = 1 ∨ s = -1)
    (σ : ARel n) (hσ : C.alpha = .rel σ) (hwp : σ.cons.wp (Zq q') = some (signed s (Zq q)))
    (hc : ∀ A ∈ σ.axes, PP.Comm (Zq q) A) (θ θ' : ℝ) :
    (seq (seq (rz q θ) C) (rz q' θ')).psem = (seq (rz q (θ + s * θ')) C).psem := by
  refine program_folding_of_runs C q q' s hs (fun K hK => ?_) θ θ'
  have h := Program.soundness C hK
  rw [hσ] at h
  exact maps_of_mem_gamma hwp hc h

/-- Corollary 5.13 with the disjunctive refinement of Section 7 (RQ2, "Configurations"): if the
fold condition holds, with the same sign, in every relation that the disjunctive analysis keeps
for `C` (with any bound `N`), the merge is correct. -/
theorem program_folding_disj (N : ℕ) (C : Program n) (q q' : Fin n) (s : ℤ)
    (hs : s = 1 ∨ s = -1)
    (hL : ∀ τ ∈ Disj.run N C [alphaI n], ∃ σ : ARel n, τ = .rel σ ∧
      σ.cons.wp (Zq q') = some (signed s (Zq q)) ∧ ∀ A ∈ σ.axes, PP.Comm (Zq q) A)
    (θ θ' : ℝ) :
    (seq (seq (rz q θ) C) (rz q' θ')).psem = (seq (rz q (θ + s * θ')) C).psem := by
  refine program_folding_of_runs C q q' s hs (fun K hK => ?_) θ θ'
  obtain ⟨τ, hτ, hK⟩ := Disj.program_sound N C hK
  obtain ⟨σ, rfl, hwp, hc⟩ := hL τ hτ
  exact maps_of_mem_gamma hwp hc hK

/-! ## The semantics of programs is compositional -/

theorem psem_seq_congr {C₁ C₁' C₂ C₂' : Program n} (h₁ : C₁.psem = C₁'.psem)
    (h₂ : C₂.psem = C₂'.psem) : (seq C₁ C₂).psem = (seq C₁' C₂').psem := by
  simp only [psem, h₁, h₂]

theorem psem_choice_congr {C₁ C₁' C₂ C₂' : Program n} (h₁ : C₁.psem = C₁'.psem)
    (h₂ : C₂.psem = C₂'.psem) : (choice C₁ C₂).psem = (choice C₁' C₂').psem := by
  simp only [psem, h₁, h₂]

theorem psem_loop_congr {C C' : Program n} (h : C.psem = C'.psem) :
    (loop C).psem = (loop C').psem := by
  simp only [psem, h]

/-! ## Example 5.11 as an instance -/

/-- Paper, Example 5.11: Corollary 5.13 merges the `T₁†` after the loop into the `T₁` before
it, leaving a rotation by `π/4 - π/4 = 0`. -/
theorem loopSimple_folds :
    Examples.loopSimple.psem =
      (seq (rz 0 (Real.pi / 4 + ((1 : ℤ) : ℝ) * (-(Real.pi / 4)))) Examples.loopCX).psem :=
  program_folding Examples.loopCX 0 0 1 (Or.inl rfl)
    ⟨⟨![none, some (Xq 1)], ![some (Zq 0), none]⟩, []⟩ Examples.example_6_13_stable.2
    (by decide) (by simp) _ _

end

end PauliFold.Paper
