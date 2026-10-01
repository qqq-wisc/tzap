import PauliFold.Domains

/-!
# The tableau domain refines the affine (phase-folding) domain

The containment theorems of `Domains.lean` are stated for `pauliGamma`, the concretization
of the per-string Pauli transformer. This file transfers them to the tableau abstraction
of `Tableau.lean`, the one the Rust `PhaseFoldPauli` pass computes.

The bridge is `tab_gamma_subset_pauli`: every fact of the Pauli domain is also a fact of
the tableau. The proof runs the two analyses side by side (`joint_run`). For a fixed input
Pauli `B`, while the forward Pauli run is at some `P`:

- the tableau rows represent the Clifford part `C`, and `C B C† = P`;
- `B` commutes with every axis the tableau has recorded.

A Clifford gate conjugates both sides. A T on `q` keeps the Pauli run alive only if `P`
has `I` or `Z` at `q`. Then `P` commutes with `Z_q`, so `B = C† P C` commutes with the new
axis `C† Z_q C`.

Consequences:

- `tableau_refines_phase`: the tableau refines phase folding on every circuit.
- `tableau_refines_phase_strict`: strictly, on `h; h`.
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Commutation, combinatorially and as matrices -/

theorem ipow_ne_zero (k : ZMod 4) : ipow k ≠ 0 := pow_ne_zero _ Complex.I_ne_zero

theorem ipow_eq_one {k : ZMod 4} (h : ipow k = 1) : k = 0 := by
  have hk := ZMod.val_lt k
  unfold ipow at h
  interval_cases hv : k.val
  · exact (ZMod.val_eq_zero k).1 hv
  · have := congrArg Complex.im h; simp at this
  · simp [Complex.I_sq] at h; norm_num at h
  · rw [pow_succ, Complex.I_sq] at h
    have := congrArg Complex.im h
    simp at this

theorem ipow_inj {a b : ZMod 4} (h : ipow a = ipow b) : a = b := by
  have hb := ipow_ne_zero b
  have : ipow (a - b) = 1 := by
    have e := ipow_add (a - b) b
    rw [sub_add_cancel, h] at e
    exact (mul_eq_right₀ hb).1 e.symm
  exact sub_eq_zero.1 (ipow_eq_one this)

namespace PP

theorem toMatrix_letters (P : PP n) : P.toMatrix = ipow P.phase • (⟨0, P.letter⟩ : PP n).toMatrix := by
  funext o i
  simp [toMatrix]

/-- A phase-0 string has a nonzero entry. -/
theorem toMatrix_base_ne_zero (L : Fin n → Letter) : (⟨0, L⟩ : PP n).toMatrix ≠ 0 := by
  intro h
  have h' := congrFun (congrFun h fun r => decide (L r = .X ∨ L r = .Y)) fun _ => false
  simp only [toMatrix, ipow_zero, one_mul, Matrix.zero_apply] at h'
  rw [Finset.prod_eq_zero_iff] at h'
  obtain ⟨r, -, hr⟩ := h'
  cases hL : L r <;> simp [hL, Letter.mat] at hr

/-- `Comm` is commutation of the matrices. -/
theorem comm_iff (P A : PP n) : Comm P A ↔ P.toMatrix * A.toMatrix = A.toMatrix * P.toMatrix := by
  refine ⟨Comm.toMatrix, fun h => ?_⟩
  rw [← toMatrix_mul, ← toMatrix_mul, toMatrix_letters (P.mul A), toMatrix_letters (A.mul P)] at h
  have hL : (P.mul A).letter = (A.mul P).letter := by
    funext r; exact mulLetter_snd_comm _ _
  rw [hL] at h
  have h0 := toMatrix_base_ne_zero (A.mul P).letter
  have : (ipow (P.mul A).phase - ipow (A.mul P).phase) •
      (⟨0, (A.mul P).letter⟩ : PP n).toMatrix = 0 := by
    rw [sub_smul, h, sub_self]
  rcases smul_eq_zero.1 this with h1 | h1
  · exact ipow_inj (sub_eq_zero.1 h1)
  · exact absurd h1 h0

end PP

/-! ## Running the Pauli domain and the tableau together -/

theorem gateZ_eq (q : Fin n) : embed1 n (diag2 1 (ep 1)) q.val = (PP.single q .Z).toMatrix := by
  rw [PP.toMatrix_single]
  congr 1
  funext o i
  cases o <;> cases i <;> simp [diag2, Letter.mat, PauliFold.ep_one]

/-- A Pauli with `I` or `Z` at `q` commutes with `Z_q`. -/
theorem comm_Z_of_letter (q : Fin n) (P : SPauli n) (hL : P.letter q = .I ∨ P.letter q = .Z) :
    P.toMatrix * (PP.single q .Z).toMatrix = (PP.single q .Z).toMatrix * P.toMatrix := by
  have h := conj_diag_fix 1 q P hL
  rw [gateZ_eq] at h
  set Z := (PP.single q .Z).toMatrix
  have hZ : Zᴴ * Z = 1 := by
    have := gate_adj_mul (n := n) (.z q)
    rwa [show gateUnitary n (CTGate.z q).toGate = Z from gateZ_eq q] at this
  unfold conj at h
  calc P.toMatrix * Z = (Z * P.toMatrix * Zᴴ) * Z := by rw [h]
    _ = Z * P.toMatrix := by rw [Matrix.mul_assoc, hZ, Matrix.mul_one]

theorem letter_of_not_anti (L : Letter) (h : ¬ L.anticommutesZ = true) : L = .I ∨ L = .Z := by
  cases L <;> simp_all [Letter.anticommutesZ]

/-- The joint invariant: the rows represent `C`, and while the Pauli run from `B` is at some
`P`, `C B C† = P` and `B` commutes with every axis. -/
def Joint (B : SPauli n) (T : Tab n) (C : Density n) (p : PauliAbs n) : Prop :=
  C * Cᴴ = 1 ∧ T.RowsOK C ∧
    ∀ P, p = .one P → C * B.toMatrix * Cᴴ = P.toMatrix ∧
      ∀ A ∈ T.axes, B.toMatrix * A.toMatrix = A.toMatrix * B.toMatrix

theorem joint_step {B : SPauli n} {T : Tab n} {C : Density n} {p : PauliAbs n}
    (h : Joint B T C p) (g : CTGate n) :
    ∃ C', Joint B (T.step g) C' (PauliAbs.step g p) := by
  obtain ⟨hCC, hrows, hP⟩ := h
  have hC : Cᴴ * C = 1 := mul_eq_one_comm.mp hCC
  by_cases hT : ∃ q, g = .t q ∨ g = .tdg q
  · obtain ⟨q, hq⟩ := hT
    have hstep : T.step g = { T with axes := T.rz q :: T.axes } := by
      rcases hq with rfl | rfl <;> rfl
    refine ⟨C, hCC, by rw [hstep]; exact hrows, fun P' hP' => ?_⟩
    cases p with
    | top => rcases hq with rfl | rfl <;> simp [PauliAbs.step] at hP'
    | one P =>
      have hkeep : PauliAbs.step g (.one P) = .one P' →
          P' = P ∧ (P.letter q = .I ∨ P.letter q = .Z) := by
        intro he
        rcases hq with rfl | rfl <;>
        · simp only [PauliAbs.step] at he
          split_ifs at he with hanti <;> cases he
          exact ⟨rfl, letter_of_not_anti _ hanti⟩
      obtain ⟨rfl, hL⟩ := hkeep hP'
      obtain ⟨hconj, haxes⟩ := hP P' rfl
      refine ⟨hconj, fun A hA => ?_⟩
      rw [hstep] at hA
      rcases List.mem_cons.1 hA with rfl | hA
      · -- `B = C† P C` commutes with `z[q] = C† Z_q C`.
        have hB : B.toMatrix = Cᴴ * P'.toMatrix * C := by
          rw [← hconj]; simp only [Matrix.mul_assoc]; rw [hC, Matrix.mul_one, ← Matrix.mul_assoc,
            hC, Matrix.one_mul]
        rw [(hrows q).2, hB]
        have hPZ := comm_Z_of_letter q P' hL
        simp only [Matrix.mul_assoc]
        rw [← Matrix.mul_assoc C Cᴴ, hCC, Matrix.one_mul, ← Matrix.mul_assoc C Cᴴ, hCC,
          Matrix.one_mul, ← Matrix.mul_assoc P'.toMatrix, hPZ, Matrix.mul_assoc]
      · exact haxes A hA
  · push Not at hT
    have hcl : g.IsClifford := by
      cases g with
      | t q => exact ((hT q).1 rfl).elim
      | tdg q => exact ((hT q).2 rfl).elim
      | _ => trivial
    set G := gateUnitary n g.toGate
    have hGG : G * Gᴴ = 1 := gate_mul_adj g
    refine ⟨G * C, ?_, Tab.rowsOK_step hCC hrows hcl, fun P' hP' => ?_⟩
    · rw [Matrix.conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc C, hCC,
        Matrix.one_mul, hGG]
    · cases p with
      | top => cases g <;> simp_all [PauliAbs.step, CTGate.IsClifford]
      | one P =>
        have hstepP : PauliAbs.step g (.one P) = .one (P.conjClifford g) := by
          cases g <;> simp_all [PauliAbs.step, CTGate.IsClifford]
        rw [hstepP] at hP'
        cases hP'
        obtain ⟨hconj, haxes⟩ := hP P rfl
        refine ⟨?_, ?_⟩
        · have hs := step_sound g (.one P) P.toMatrix (Set.mem_singleton _)
          rw [hstepP] at hs
          simp only [PauliAbs.gamma, Set.mem_singleton_iff, CTGate.concrete, conj] at hs
          rw [← hs, ← hconj, Matrix.conjTranspose_mul]
          simp only [Matrix.mul_assoc]; rfl
        · intro A hA
          rw [Tab.axes_step_clifford hcl] at hA
          exact haxes A hA

theorem joint_run (B : SPauli n) (gs : List (CTGate n)) :
    ∀ {T : Tab n} {C : Density n} {p : PauliAbs n}, Joint B T C p →
      ∃ C', Joint B (T.run gs) C' (PauliAbs.run gs p) := by
  induction gs with
  | nil => intro T C p h; exact ⟨C, h⟩
  | cons g gs ih =>
    intro T C p h
    obtain ⟨C₁, h₁⟩ := joint_step h g
    exact ih h₁

theorem joint_init (B : SPauli n) : Joint B (Tab.init n) 1 (.one B) :=
  ⟨by simp, fun r => ⟨by simp [Tab.init], by simp [Tab.init]⟩,
    fun P hP => by cases hP; exact ⟨by simp, fun A hA => absurd hA List.not_mem_nil⟩⟩

/-! ## Containment -/

/-- **Every Pauli fact is a tableau fact**, so the tableau concretization is contained in the
Pauli one. -/
theorem tab_gamma_subset_pauli (gs : List (CTGate n)) : Tab.gamma gs ⊆ pauliGamma gs := by
  intro U hU f hf
  obtain ⟨C, hCC, hrows, hP⟩ := joint_run f.1 gs (joint_init f.1)
  have hC : Cᴴ * C = 1 := mul_eq_one_comm.mp hCC
  obtain ⟨hconj, haxes⟩ := hP f.2 hf
  set T := (Tab.init n).run gs
  -- The tableau's pull-back of `f.2` is `f.1`.
  have hback : (T.back (PP.ofS f.2)).toMatrix = f.1.toMatrix := by
    rw [Tab.back_ok hC hCC hrows, PP.toMatrix_ofS, ← hconj]
    simp only [Matrix.mul_assoc]
    rw [hC, Matrix.mul_one, ← Matrix.mul_assoc, hC, Matrix.one_mul]
  have hfact := hU (PP.ofS f.2) fun A hA => by
    rw [PP.comm_iff, hback]
    exact haxes A hA
  rw [hback, PP.toMatrix_ofS] at hfact
  exact hfact

/-- **The tableau refines phase folding** on every circuit. -/
theorem tableau_refines_phase (gs : List (CTGate n)) : Tab.gamma gs ⊆ phaseGamma gs :=
  (tab_gamma_subset_pauli gs).trans (pauli_refines_phase gs)

/-- **Strictly:** on `h; h` the tableau proves `Z ↦ Z` and phase folding proves nothing. -/
theorem tableau_refines_phase_strict : Tab.gamma hh ⊂ phaseGamma hh :=
  ⟨tableau_refines_phase hh, fun h =>
    pauli_refines_phase_strict.2 (h.trans (tab_gamma_subset_pauli hh))⟩

end

end PauliFold
