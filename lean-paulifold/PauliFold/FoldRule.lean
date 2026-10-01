import PauliFold.Tableau
import PauliFold.StateFold.Incomparable

/-!
# The fold theorem

The rewrite `PhaseFoldPauli` performs, justified by one fact of the tableau abstraction.

Let `m` be a segment of a circuit and `σ` its abstract state from the identity tableau.
Suppose that
- (a) the pull-back of `Z_{q'}` through `m` is `± Z_q`, and
- (b) `Z_q` commutes with every axis of `σ`.

Then a `Z_q` rotation just before `m` slides through `m` and becomes a `Z_{q'}` rotation
just after it (`fold_pos`, `fold_neg`). There it merges with a rotation on `q'`
(`fold_merge_pos`, `fold_merge_neg`). For gates: `T_q; m; T_{q'}` equals `m; S_{q'}` when the
sign is `+` (`fold_t_t`), and equals `m` up to a global phase when the sign is `−`
(`fold_t_t_neg`).

Rotations are phase-exact here: `rot θ P = diag(1, e^{iπθ})` about `P`, so `T = rot (1/4)`
and `S = rot (1/2)` exactly. The paper's `R_P(θ) = e^{-iθP/2}` differs by a global phase.
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Rotations about a Pauli string -/

/-- The rotation `diag(1, e^{iπθ})` about `P`: `((1 + e^{iπθ})/2)·I + ((1 − e^{iπθ})/2)·P`. -/
def rot (θ : ℚ) (P : PP n) : Density n :=
  ((1 + ep θ) / 2) • (1 : Density n) + ((1 - ep θ) / 2) • P.toMatrix

theorem rot_zero (P : PP n) : rot 0 P = 1 := by
  simp [rot]

/-- Rotations about a fixed involution add their angles. -/
theorem rot_mul {P : PP n} (hP : P.toMatrix * P.toMatrix = 1) (θ φ : ℚ) :
    rot θ P * rot φ P = rot (θ + φ) P := by
  simp only [rot, add_mul, mul_add, smul_mul_assoc, mul_smul_comm, one_mul, mul_one, hP,
    ep_add]
  module

theorem single_Z_sq (q : Fin n) :
    (PP.single q .Z).toMatrix * (PP.single q .Z).toMatrix = 1 := by
  rw [← PP.toMatrix_mul, PP.mul_single_eq, ← PP.toMatrix_one]
  congr 1
  apply PP.ext'
  · rfl
  · intro r; simp [mulLetter, PP.one]

theorem gate_t (q : Fin n) : gateUnitary n (CTGate.t q).toGate = rot (1 / 4) (PP.single q .Z) := by
  rw [rot, PP.toMatrix_single]; exact embed1_diag2 (ep (1 / 4)) q

theorem gate_tdg (q : Fin n) :
    gateUnitary n (CTGate.tdg q).toGate = rot (-1 / 4) (PP.single q .Z) := by
  rw [rot, PP.toMatrix_single]; exact embed1_diag2 (ep (-1 / 4)) q

theorem gate_s (q : Fin n) : gateUnitary n (CTGate.s q).toGate = rot (1 / 2) (PP.single q .Z) := by
  rw [rot, PP.toMatrix_single]; exact embed1_diag2 (ep (1 / 2)) q

/-! ## Sliding a rotation through a segment -/

theorem PP.comm_scale (P A : PP n) (k : ZMod 4) : PP.Comm (P.scale k) A ↔ PP.Comm P A := by
  simp only [PP.Comm, PP.mul, PP.scale]
  constructor <;> intro h <;> linear_combination h

theorem ipow_zero' : ipow 0 = 1 := by simp [ipow]

/-- **The fact.** Under (a) and (b), the segment maps `Z_q` to `i^{-k} Z_{q'}`. -/
theorem segment_Z (m : List (CTGate n)) (q q' : Fin n) (k : ZMod 4)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = (PP.single q .Z).scale k)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A) :
    circ m * (PP.single q .Z).toMatrix = ipow (-k) • (PP.single q' .Z).toMatrix * circ m := by
  have h := Tab.tableau_sound m (PP.single q' .Z)
    (by rw [ha]; exact fun A hA => (PP.comm_scale _ _ _).2 (hb A hA))
  rw [ha, PP.toMatrix_scale, conj] at h
  have hU : (circ m)ᴴ * circ m = 1 := mul_eq_one_comm.mp (SF.unitary_mul_conjTranspose m)
  have h2 : ipow k • (circ m * (PP.single q .Z).toMatrix) =
      (PP.single q' .Z).toMatrix * circ m := by
    have : circ m * (ipow k • (PP.single q .Z).toMatrix) * (circ m)ᴴ * circ m =
        (PP.single q' .Z).toMatrix * circ m := by rw [h]
    rwa [Matrix.mul_assoc, hU, Matrix.mul_one, Matrix.mul_smul] at this
  calc circ m * (PP.single q .Z).toMatrix
      = (ipow (-k) * ipow k) • (circ m * (PP.single q .Z).toMatrix) := by
        rw [← ipow_add, neg_add_cancel, ipow_zero', one_smul]
    _ = ipow (-k) • (PP.single q' .Z).toMatrix * circ m := by
        rw [← smul_smul, h2, Matrix.smul_mul]

/-- A combination `a·I + b·Z_q` slides through the segment. -/
theorem slide (m : List (CTGate n)) (q q' : Fin n) (k : ZMod 4)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = (PP.single q .Z).scale k)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A) (a b : ℂ) :
    circ m * (a • 1 + b • (PP.single q .Z).toMatrix) =
      (a • 1 + (b * ipow (-k)) • (PP.single q' .Z).toMatrix) * circ m := by
  rw [Matrix.mul_add, Matrix.mul_smul, Matrix.mul_smul, Matrix.mul_one,
    segment_Z m q q' k ha hb, Matrix.add_mul]
  simp only [Matrix.smul_mul, Matrix.one_mul, smul_smul]

/-- **Fold theorem, sign `+`.** If `β(Z_{q'}) = Z_q` and `Z_q` commutes with every axis, a
rotation about `Z_q` before the segment equals the same rotation about `Z_{q'}` after it. -/
theorem fold_pos (m : List (CTGate n)) (q q' : Fin n)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = PP.single q .Z)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A) (θ : ℚ) :
    circ m * rot θ (PP.single q .Z) = rot θ (PP.single q' .Z) * circ m := by
  have ha' : ((Tab.init n).run m).back (PP.single q' .Z) = (PP.single q .Z).scale 0 := by
    rw [ha]; rfl
  rw [rot, slide m q q' 0 ha' hb, neg_zero, ipow_zero', mul_one, rot]

theorem ipow_neg_two : ipow (-2) = -1 := by
  show Complex.I ^ (2 : ℕ) = -1
  exact Complex.I_sq

/-- **Fold theorem, sign `−`.** If `β(Z_{q'}) = −Z_q`, the rotation comes out with the
opposite angle, up to the global phase `e^{iπθ}`. -/
theorem fold_neg (m : List (CTGate n)) (q q' : Fin n)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = (PP.single q .Z).scale 2)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A) (θ : ℚ) :
    circ m * rot θ (PP.single q .Z) = (ep θ • rot (-θ) (PP.single q' .Z)) * circ m := by
  rw [rot, slide m q q' 2 ha hb, ipow_neg_two]
  have hab : ep θ * ep (-θ) = 1 := by rw [← ep_add, add_neg_cancel, ep_zero]
  congr 1
  rw [rot, smul_add, smul_smul, smul_smul]
  congr 2
  · linear_combination (-(1 / 2 : ℂ)) * hab
  · linear_combination (1 / 2 : ℂ) * hab

/-! ## Merging -/

/-- **Merge, sign `+`** (Theorem 9.1 of the notes): `ℓ; R_q(θ); m; R_{q'}(θ'); k` equals
`ℓ; m; R_{q'}(θ + θ'); k`. -/
theorem fold_merge_pos (m : List (CTGate n)) (q q' : Fin n)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = PP.single q .Z)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A)
    (L K : Density n) (θ θ' : ℚ) :
    K * rot θ' (PP.single q' .Z) * circ m * rot θ (PP.single q .Z) * L =
      K * rot (θ + θ') (PP.single q' .Z) * circ m * L := by
  calc K * rot θ' (PP.single q' .Z) * circ m * rot θ (PP.single q .Z) * L
      = K * rot θ' (PP.single q' .Z) * (circ m * rot θ (PP.single q .Z)) * L := by
        simp only [Matrix.mul_assoc]
    _ = K * (rot θ' (PP.single q' .Z) * rot θ (PP.single q' .Z)) * circ m * L := by
        rw [fold_pos m q q' ha hb]; simp only [Matrix.mul_assoc]
    _ = K * rot (θ + θ') (PP.single q' .Z) * circ m * L := by
        rw [rot_mul (single_Z_sq q'), add_comm]

/-- **Merge, sign `−`**: the angles subtract, up to a global phase. -/
theorem fold_merge_neg (m : List (CTGate n)) (q q' : Fin n)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = (PP.single q .Z).scale 2)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A)
    (L K : Density n) (θ θ' : ℚ) :
    K * rot θ' (PP.single q' .Z) * circ m * rot θ (PP.single q .Z) * L =
      ep θ • (K * rot (-θ + θ') (PP.single q' .Z) * circ m * L) := by
  calc K * rot θ' (PP.single q' .Z) * circ m * rot θ (PP.single q .Z) * L
      = K * rot θ' (PP.single q' .Z) * (circ m * rot θ (PP.single q .Z)) * L := by
        simp only [Matrix.mul_assoc]
    _ = ep θ • (K * (rot θ' (PP.single q' .Z) * rot (-θ) (PP.single q' .Z)) * circ m * L) := by
        rw [fold_neg m q q' ha hb]
        simp only [Matrix.mul_assoc, Matrix.smul_mul, Matrix.mul_smul]
    _ = ep θ • (K * rot (-θ + θ') (PP.single q' .Z) * circ m * L) := by
        rw [rot_mul (single_Z_sq q'), add_comm θ']

/-! ## Gate-level corollaries -/

theorem circ_split (ℓ m k : List (CTGate n)) (g g' : CTGate n) :
    circ (ℓ ++ g :: m ++ g' :: k) =
      circ k * gateUnitary n g'.toGate * circ m * gateUnitary n g.toGate * circ ℓ := by
  simp only [circ, List.map_append, List.map_cons, unitary_append, unitary_cons,
    Matrix.mul_assoc]

/-- `ℓ; T_q; m; T_{q'}; k` equals `ℓ; m; S_{q'}; k` when the segment's fact has sign `+`. -/
theorem fold_t_t (ℓ m k : List (CTGate n)) (q q' : Fin n)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = PP.single q .Z)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A) :
    circ (ℓ ++ .t q :: m ++ .t q' :: k) = circ (ℓ ++ m ++ .s q' :: k) := by
  rw [circ_split, gate_t, gate_t, fold_merge_pos m q q' ha hb]
  simp only [circ, List.map_append, List.map_cons, unitary_append, unitary_cons,
    Matrix.mul_assoc, gate_s]
  norm_num

/-- `ℓ; T_q; m; T_{q'}; k` equals `ℓ; m; k` up to the global phase `e^{iπ/4}` when the
segment's fact has sign `−`: the two `T`s cancel. -/
theorem fold_t_t_neg (ℓ m k : List (CTGate n)) (q q' : Fin n)
    (ha : ((Tab.init n).run m).back (PP.single q' .Z) = (PP.single q .Z).scale 2)
    (hb : ∀ A ∈ ((Tab.init n).run m).axes, PP.Comm (PP.single q .Z) A) :
    circ (ℓ ++ .t q :: m ++ .t q' :: k) = ep (1 / 4) • circ (ℓ ++ m ++ k) := by
  rw [circ_split, gate_t, gate_t, fold_merge_neg m q q' ha hb, neg_add_cancel, rot_zero]
  simp only [circ, List.map_append, unitary_append, Matrix.mul_one, Matrix.mul_assoc]

/-- Example: the segment `cx 0 1; t 1; cx 0 1` of the worked example (gates 2–4) satisfies
the hypotheses with `q = q' = 0`, so gates 1 and 5 merge into an `S`. -/
example : circ ([] ++ .t (0 : Fin 2) :: [.cx 0 1 (by decide), .t 1, .cx 0 1 (by decide)] ++
      .t 0 :: []) = circ ([] ++ [.cx 0 1 (by decide), .t 1, .cx 0 1 (by decide)] ++ [.s 0]) :=
  fold_t_t [] _ [] 0 0 (by decide) (by decide)

end

end PauliFold
