import PauliFold.Soundness
import TzapLean.SemanticsCheck
import Mathlib.LinearAlgebra.Matrix.SemiringInverse

/-!
# A tableau abstraction of Clifford+T circuits

The abstraction the Rust `PhaseFoldPauli` pass computes, for whole circuits.

**The state.** A Clifford tableau and a list of axes:

- The rows record, for each qubit `q`, the Paulis `x[q] = C† X_q C` and `z[q] = C† Z_q C`,
  where `C` is the Clifford part of the circuit so far.
- The axes are input-frame Paulis, one per T or T† gate.

**The transformers.**

- A Clifford gate updates the rows and leaves the axes alone.
- A T or T† on `q` leaves the rows alone and records the current row `z[q]` as an axis.

**Why the T step is sound.** A T is a linear combination `α·I + β·Z_q`. Pulled back to the
input frame it is `α·I + β·z[q]`, so it commutes with every Pauli that commutes with
`z[q]`.

**The invariant.** The circuit's unitary is `C · V`, where the rows are correct for the
Clifford `C` and `V` commutes with every Pauli that commutes with all the axes.

**Soundness** (`tableau_sound`). If the pull-back `C† P C` of a Pauli `P` commutes with every
axis, the circuit maps it back to `P`: `U (C† P C) U† = P`.

Pauli strings here carry a phase `i^k`, so that products of strings are strings again.
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

/-! ## Phases `i^k` -/

/-- `i^k` for `k : ZMod 4`. -/
def ipow (k : ZMod 4) : ℂ := Complex.I ^ k.val

theorem I_pow_mod (m : ℕ) : Complex.I ^ (m % 4) = Complex.I ^ m := by
  conv_rhs => rw [← Nat.div_add_mod m 4, pow_add, pow_mul, Complex.I_pow_four, one_pow, one_mul]

theorem ipow_add (a b : ZMod 4) : ipow (a + b) = ipow a * ipow b := by
  unfold ipow
  rw [ZMod.val_add, I_pow_mod, pow_add]

@[simp] theorem ipow_zero : ipow 0 = 1 := by simp [ipow]
@[simp] theorem ipow_one : ipow 1 = Complex.I := by
  unfold ipow; rw [show (1 : ZMod 4).val = 1 from rfl, pow_one]
@[simp] theorem ipow_two : ipow 2 = -1 := by
  unfold ipow; rw [show (2 : ZMod 4).val = 2 from rfl, Complex.I_sq]
@[simp] theorem ipow_three : ipow 3 = -Complex.I := by
  unfold ipow; rw [show (3 : ZMod 4).val = 3 from rfl, pow_succ, Complex.I_sq]; ring

theorem prod_ipow {ι : Type*} (s : Finset ι) (f : ι → ZMod 4) :
    ∏ r ∈ s, ipow (f r) = ipow (∑ r ∈ s, f r) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih => rw [Finset.prod_insert ha, Finset.sum_insert ha, ih, ipow_add]

/-! ## Letter products -/

/-- `a · b = i^k · c`. -/
def mulLetter : Letter → Letter → ZMod 4 × Letter
  | .I, b => (0, b)
  | a, .I => (0, a)
  | .X, .X | .Y, .Y | .Z, .Z => (0, .I)
  | .X, .Y => (1, .Z)
  | .Y, .X => (3, .Z)
  | .Y, .Z => (1, .X)
  | .Z, .Y => (3, .X)
  | .Z, .X => (1, .Y)
  | .X, .Z => (3, .Y)

theorem mulLetter_mat (a b : Letter) (o i : Bool) :
    ∑ c : Bool, a.mat o c * b.mat c i = ipow (mulLetter a b).1 * (mulLetter a b).2.mat o i := by
  cases a <;> cases b <;> cases o <;> cases i <;>
    simp [Letter.mat, mulLetter, Fintype.sum_bool] <;> ring_nf <;> simp [Complex.I_sq]

/-! ## Pauli strings with a phase -/

/-- `i^phase · ⊗_r letter r`. -/
structure PP (n : ℕ) where
  phase : ZMod 4
  letter : Fin n → Letter

namespace PP

variable {n : ℕ}

instance : DecidableEq (PP n) := fun a b =>
  decidable_of_iff (a.phase = b.phase ∧ a.letter = b.letter) (by cases a; cases b; simp)

def toMatrix (P : PP n) : Density n :=
  fun o i => ipow P.phase * ∏ r, (P.letter r).mat (o r) (i r)

/-- The product of two strings. -/
def mul (P Q : PP n) : PP n :=
  ⟨P.phase + Q.phase + ∑ r, (mulLetter (P.letter r) (Q.letter r)).1,
    fun r => (mulLetter (P.letter r) (Q.letter r)).2⟩

/-- **The matrix of a product is the product of the matrices.** -/
theorem toMatrix_mul (P Q : PP n) : (P.mul Q).toMatrix = P.toMatrix * Q.toMatrix := by
  funext o i
  rw [Matrix.mul_apply]
  simp only [toMatrix, mul]
  have hsum : ∑ k : Basis n, ipow P.phase * (∏ r, (P.letter r).mat (o r) (k r)) *
        (ipow Q.phase * ∏ r, (Q.letter r).mat (k r) (i r)) =
      ipow P.phase * ipow Q.phase *
        ∏ r, ∑ c : Bool, (P.letter r).mat (o r) c * (Q.letter r).mat c (i r) := by
    rw [Fintype.prod_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Finset.prod_mul_distrib]
    ring
  rw [hsum]
  simp only [mulLetter_mat, Finset.prod_mul_distrib, prod_ipow, ipow_add]
  ring

/-- Multiply by `i^k`. -/
def scale (k : ZMod 4) (P : PP n) : PP n := ⟨P.phase + k, P.letter⟩

theorem toMatrix_scale (k : ZMod 4) (P : PP n) : (P.scale k).toMatrix = ipow k • P.toMatrix := by
  funext o i
  simp only [toMatrix, scale, ipow_add, Matrix.smul_apply, smul_eq_mul]
  ring

/-- The identity string. -/
def one : PP n := ⟨0, fun _ => .I⟩

theorem toMatrix_one : (one : PP n).toMatrix = 1 := by
  funext o i
  simp only [toMatrix, one, ipow_zero, one_mul, Matrix.one_apply]
  by_cases h : o = i
  · subst h; simp [Letter.mat]
  · rw [if_neg h]
    obtain ⟨r, hr⟩ : ∃ r, o r ≠ i r := by
      by_contra hc; push Not at hc; exact h (funext hc)
    exact Finset.prod_eq_zero (Finset.mem_univ r) (by simp [Letter.mat, hr])

/-- The letter `L` on qubit `q`, the identity elsewhere. -/
def single (q : Fin n) (L : Letter) : PP n := ⟨0, Function.update (fun _ => .I) q L⟩

theorem toMatrix_single (q : Fin n) (L : Letter) :
    (single q L).toMatrix = embed1 n L.mat q.val := by
  funext o i
  rw [embed1_apply_of_lt _ q.isLt]
  simp only [toMatrix, single, ipow_zero, one_mul]
  rw [prod_split _ q, Function.update_self, Fin.eta]
  have hrest : ∏ r ∈ Finset.univ.erase q,
      (Function.update (fun _ => Letter.I) q L r).mat (o r) (i r) =
      if ∀ r : Fin n, (r : ℕ) ≠ q.val → o r = i r then 1 else 0 := by
    split_ifs with h
    · refine Finset.prod_eq_one fun r hr => ?_
      have hr' : r ≠ q := Finset.ne_of_mem_erase hr
      simp [hr', Letter.mat, h r (Fin.val_ne_of_ne hr')]
    · push Not at h
      obtain ⟨r, hr, hne⟩ := h
      have hr' : r ≠ q := fun h => hr (congrArg Fin.val h)
      exact Finset.prod_eq_zero (Finset.mem_erase.2 ⟨hr', Finset.mem_univ r⟩)
        (by simp [hr', Letter.mat, hne])
  rw [hrest]
  split_ifs <;> simp

/-! ## Decomposition into single-qubit factors -/

/-- The part of `P`'s letters on the qubits in `l`, with phase 0. -/
def restrict (P : PP n) (l : List (Fin n)) : PP n :=
  ⟨0, fun r => if r ∈ l then P.letter r else .I⟩

theorem mulLetter_I_left (b : Letter) : mulLetter .I b = (0, b) := by cases b <;> rfl
theorem mulLetter_I_right (a : Letter) : mulLetter a .I = (0, a) := by cases a <;> rfl

theorem restrict_nil (P : PP n) : P.restrict [] = one := rfl

theorem restrict_cons (P : PP n) (r : Fin n) (l : List (Fin n)) (hr : r ∉ l) :
    P.restrict (r :: l) = (single r (P.letter r)).mul (P.restrict l) := by
  have hmul : ∀ s : Fin n, mulLetter ((single r (P.letter r)).letter s) ((P.restrict l).letter s) =
      (0, (P.restrict (r :: l)).letter s) := fun s => by
    by_cases hs : s = r
    · subst hs
      simp [single, restrict, hr, mulLetter_I_right]
    · simp [single, restrict, hs, Function.update_of_ne hs, mulLetter_I_left]
  simp only [mul, hmul, Finset.sum_const_zero, add_zero]
  rfl

theorem restrict_all (P : PP n) : (P.restrict (List.finRange n)).scale P.phase = P := by
  cases P
  simp [restrict, scale]

/-! ## Commutation -/

/-- `P` and `A` commute: their two products have the same phase (their letters always
agree). -/
def Comm (P A : PP n) : Prop := (P.mul A).phase = (A.mul P).phase

instance (P A : PP n) : Decidable (Comm P A) := inferInstanceAs (Decidable (_ = _))

theorem mulLetter_snd_comm (a b : Letter) : (mulLetter a b).2 = (mulLetter b a).2 := by
  cases a <;> cases b <;> rfl

theorem Comm.toMatrix {P A : PP n} (h : Comm P A) :
    P.toMatrix * A.toMatrix = A.toMatrix * P.toMatrix := by
  rw [← toMatrix_mul, ← toMatrix_mul]
  have : P.mul A = A.mul P := by
    cases hPA : P.mul A with
    | mk k L =>
      cases hAP : A.mul P with
      | mk k' L' =>
        have hk : k = k' := by
          have := h; unfold Comm at this; rw [hPA, hAP] at this; exact this
        have hL : L = L' := by
          funext r
          have e1 := congrArg (fun X => X.letter r) hPA
          have e2 := congrArg (fun X => X.letter r) hAP
          simp only [mul] at e1 e2
          rw [← e1, ← e2, mulLetter_snd_comm]
        rw [hk, hL]
  rw [this]

end PP

/-! ## Adjoints of the gates -/

/-- The inverse of a gate. -/
def CTGate.inv {n : ℕ} : CTGate n → CTGate n
  | .s q => .sdg q
  | .sdg q => .s q
  | .t q => .tdg q
  | .tdg q => .t q
  | g => g

theorem star_diag2_ep (θ : ℚ) :
    (fun o i => star (diag2 1 (ep θ) i o)) = diag2 1 (ep (-θ)) := by
  funext o i
  cases o <;> cases i <;> simp [diag2, star_ep]

theorem gateUnitary_inv {n : ℕ} (g : CTGate n) :
    gateUnitary n g.inv.toGate = (gateUnitary n g.toGate)ᴴ := by
  cases g with
  | h q =>
    show embed1 n h2 q.val = (embed1 n h2 q.val)ᴴ
    rw [embed1_conjTranspose]
    congr 1
    funext o i
    cases o <;> cases i <;> simp [h2, Complex.conj_ofReal]
  | x q =>
    show embed1 n x2 q.val = (embed1 n x2 q.val)ᴴ
    rw [embed1_conjTranspose]
    congr 1
    funext o i
    cases o <;> cases i <;> simp [x2]
  | z q =>
    show embed1 n (diag2 1 (ep 1)) q.val = (embed1 n (diag2 1 (ep 1)) q.val)ᴴ
    rw [embed1_conjTranspose, star_diag2_ep]
    have h : ep (-1) = -1 := by rw [← star_ep, ep_one]; simp
    rw [h, ep_one]
  | s q =>
    show embed1 n (diag2 1 (ep (-1 / 2))) q.val = (embed1 n (diag2 1 (ep (1 / 2))) q.val)ᴴ
    rw [embed1_conjTranspose, star_diag2_ep]
    norm_num
  | sdg q =>
    show embed1 n (diag2 1 (ep (1 / 2))) q.val = (embed1 n (diag2 1 (ep (-1 / 2))) q.val)ᴴ
    rw [embed1_conjTranspose, star_diag2_ep]
    norm_num
  | t q =>
    show embed1 n (diag2 1 (ep (-1 / 4))) q.val = (embed1 n (diag2 1 (ep (1 / 4))) q.val)ᴴ
    rw [embed1_conjTranspose, star_diag2_ep]
    norm_num
  | tdg q =>
    show embed1 n (diag2 1 (ep (1 / 4))) q.val = (embed1 n (diag2 1 (ep (-1 / 4))) q.val)ᴴ
    rw [embed1_conjTranspose, star_diag2_ep]
    norm_num
  | cx c t hne =>
    show gateUnitary n (CTGate.cx c t hne).toGate = (gateUnitary n (CTGate.cx c t hne).toGate)ᴴ
    rw [gateUnitary_cx, permMatrix_conjTranspose_of_involutive (cxPerm_involutive c t hne)]
  | cz c t hne =>
    show gateUnitary n (CTGate.cz c t hne).toGate = (gateUnitary n (CTGate.cz c t hne).toGate)ᴴ
    rw [gateUnitary_cz, phaseMatrix_conjTranspose_of_real]
    intro b
    simp only [czPhase]
    split_ifs <;> simp

theorem ctgate_wf {n : ℕ} : ∀ g : CTGate n, g.toGate.Wf
  | .cx _ _ hne => show _ ≠ _ from Fin.val_ne_of_ne hne
  | .cz _ _ hne => show _ ≠ _ from Fin.val_ne_of_ne hne
  | .h _ | .x _ | .z _ | .s _ | .sdg _ | .t _ | .tdg _ => trivial

theorem gate_mul_adj {n : ℕ} (g : CTGate n) :
    gateUnitary n g.toGate * (gateUnitary n g.toGate)ᴴ = 1 :=
  mul_eq_one_comm.mp (gateUnitary_unitary n _ (ctgate_wf g))

theorem gate_adj_mul {n : ℕ} (g : CTGate n) :
    (gateUnitary n g.toGate)ᴴ * gateUnitary n g.toGate = 1 :=
  gateUnitary_unitary n _ (ctgate_wf g)


/-! ## Pulling a Pauli back through one Clifford gate -/

namespace PP

variable {n : ℕ}

/-- A signed Pauli as a phased one. -/
def ofS (S : SPauli n) : PP n := ⟨if S.neg then 2 else 0, S.letter⟩

theorem toMatrix_ofS (S : SPauli n) : (ofS S).toMatrix = S.toMatrix := by
  funext o i
  simp only [toMatrix, ofS, SPauli.toMatrix]
  cases S.neg <;> simp [SPauli.sgn]

theorem toMatrix_eq_smul (P : PP n) :
    P.toMatrix = ipow P.phase • (⟨false, P.letter⟩ : SPauli n).toMatrix := by
  funext o i
  simp [toMatrix, SPauli.toMatrix, SPauli.sgn]

end PP

/-- The Clifford gates. -/
def CTGate.IsClifford {n : ℕ} : CTGate n → Prop
  | .t _ | .tdg _ => False
  | _ => True

/-- `g† P g`, computed with the forward table of `g⁻¹`. -/
def pullG {n : ℕ} (g : CTGate n) (P : PP n) : PP n :=
  (PP.ofS (SPauli.conjClifford g.inv ⟨false, P.letter⟩)).scale P.phase

theorem pullG_toMatrix {n : ℕ} (g : CTGate n) (hg : g.IsClifford) (P : PP n) :
    (pullG g P).toMatrix =
      (gateUnitary n g.toGate)ᴴ * P.toMatrix * gateUnitary n g.toGate := by
  set S : SPauli n := ⟨false, P.letter⟩
  have hstep : PauliAbs.step g.inv (.one S) = .one (SPauli.conjClifford g.inv S) := by
    cases g <;> simp_all [CTGate.IsClifford, CTGate.inv, PauliAbs.step]
  have hs := step_sound g.inv (.one S) S.toMatrix (Set.mem_singleton _)
  rw [hstep] at hs
  simp only [PauliAbs.gamma, Set.mem_singleton_iff, CTGate.concrete, conj,
    gateUnitary_inv, Matrix.conjTranspose_conjTranspose] at hs
  rw [pullG, PP.toMatrix_scale, PP.toMatrix_ofS, ← hs, PP.toMatrix_eq_smul P, Matrix.mul_smul,
    Matrix.smul_mul]

/-! ## The tableau -/

/-- Rows `x[q] = C† X_q C` and `z[q] = C† Z_q C`, and the input-frame axes of the T gates. -/
structure Tab (n : ℕ) where
  rx : Fin n → PP n
  rz : Fin n → PP n
  axes : List (PP n)

namespace Tab

variable {n : ℕ}

/-- The empty circuit: each row is its own generator, and there are no axes. -/
def init (n : ℕ) : Tab n := ⟨fun q => PP.single q .X, fun q => PP.single q .Z, []⟩

/-- The pull-back of the single-qubit letter `L` on `r`. `Y = i · X Z`. -/
def rowOf (T : Tab n) (r : Fin n) : Letter → PP n
  | .I => PP.one
  | .X => T.rx r
  | .Z => T.rz r
  | .Y => ((T.rx r).mul (T.rz r)).scale 1

/-- The pull-back of `P`'s letters on the qubits in `l`, as a product of rows. -/
def backL (T : Tab n) (l : List (Fin n)) (P : PP n) : PP n :=
  l.foldr (fun r acc => (T.rowOf r (P.letter r)).mul acc) PP.one

/-- The pull-back `C† P C` of a Pauli, from the rows. -/
def back (T : Tab n) (P : PP n) : PP n := (T.backL (List.finRange n) P).scale P.phase

/-- The transformer of a gate. A Clifford gate updates the rows by the explicit rules below;
`step_rows` proves that each new row is the old pull-back `β(G† X_q G)`, respectively
`β(G† Z_q G)`. A T or T† records the current `z[q]` as an axis. -/
def step (T : Tab n) : CTGate n → Tab n
  | .h q => { T with rx := Function.update T.rx q (T.rz q), rz := Function.update T.rz q (T.rx q) }
  | .x q => { T with rz := Function.update T.rz q ((T.rz q).scale 2) }
  | .z q => { T with rx := Function.update T.rx q ((T.rx q).scale 2) }
  | .s q => { T with rx := Function.update T.rx q (((T.rx q).mul (T.rz q)).scale 3) }
  | .sdg q => { T with rx := Function.update T.rx q (((T.rx q).mul (T.rz q)).scale 1) }
  | .cx c t _ => { T with
      rx := Function.update T.rx c ((T.rx c).mul (T.rx t))
      rz := Function.update T.rz t ((T.rz c).mul (T.rz t)) }
  | .cz c t _ => { T with
      rx := Function.update (Function.update T.rx c ((T.rx c).mul (T.rz t))) t
        ((T.rz c).mul (T.rx t)) }
  | .t q | .tdg q => { T with axes := T.rz q :: T.axes }

/-- The transformer of a circuit, first gate first. -/
def run (gs : List (CTGate n)) (T : Tab n) : Tab n := gs.foldl step T

/-- The rows are correct for the Clifford `C`. -/
def RowsOK (T : Tab n) (C : Density n) : Prop :=
  ∀ r, (T.rx r).toMatrix = Cᴴ * (PP.single r .X).toMatrix * C ∧
    (T.rz r).toMatrix = Cᴴ * (PP.single r .Z).toMatrix * C

theorem single_I (r : Fin n) : PP.single r .I = (PP.one : PP n) := by
  simp only [PP.single, PP.one]
  congr 1
  funext s
  by_cases hs : s = r <;> simp [hs]

theorem single_Y (r : Fin n) :
    PP.single r .Y = ((PP.single r .X).mul (PP.single r .Z)).scale 1 := by
  have hl : ∀ s, mulLetter ((PP.single r .X).letter s) ((PP.single r .Z).letter s) =
      (if s = r then 3 else 0, (PP.single r .Y).letter s) := fun s => by
    by_cases hs : s = r
    · subst hs; simp [PP.single, mulLetter]
    · simp [PP.single, hs, Function.update_of_ne hs, mulLetter]
  simp only [PP.scale, PP.mul, hl, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  congr 1
  simp [PP.single]
  decide

theorem rowOf_ok {T : Tab n} {C : Density n} (hC : Cᴴ * C = 1) (hCC : C * Cᴴ = 1)
    (h : T.RowsOK C) (r : Fin n) (L : Letter) :
    (T.rowOf r L).toMatrix = Cᴴ * (PP.single r L).toMatrix * C := by
  cases L with
  | I => rw [single_I, rowOf, PP.toMatrix_one, Matrix.mul_one, hC]
  | X => exact (h r).1
  | Z => exact (h r).2
  | Y =>
    rw [rowOf, single_Y, PP.toMatrix_scale, PP.toMatrix_scale, PP.toMatrix_mul, PP.toMatrix_mul,
      (h r).1, (h r).2, Matrix.mul_smul, Matrix.smul_mul]
    congr 1
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc C, hCC, Matrix.one_mul]

theorem backL_ok {T : Tab n} {C : Density n} (hC : Cᴴ * C = 1) (hCC : C * Cᴴ = 1)
    (h : T.RowsOK C) (P : PP n) :
    ∀ l : List (Fin n), l.Nodup →
      (T.backL l P).toMatrix = Cᴴ * (P.restrict l).toMatrix * C
  | [], _ => by rw [backL, List.foldr_nil, PP.restrict_nil, PP.toMatrix_one, Matrix.mul_one, hC]
  | r :: l, hnd => by
    have hr : r ∉ l := (List.nodup_cons.1 hnd).1
    have ih := backL_ok hC hCC h P l (List.nodup_cons.1 hnd).2
    rw [backL, List.foldr_cons, ← backL, PP.toMatrix_mul, ih, rowOf_ok hC hCC h,
      PP.restrict_cons P r l hr, PP.toMatrix_mul]
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc C, hCC, Matrix.one_mul]

/-- **The pull-back is correct**: `back P = C† P C`. -/
theorem back_ok {T : Tab n} {C : Density n} (hC : Cᴴ * C = 1) (hCC : C * Cᴴ = 1)
    (h : T.RowsOK C) (P : PP n) : (T.back P).toMatrix = Cᴴ * P.toMatrix * C := by
  rw [back, PP.toMatrix_scale, backL_ok hC hCC h P _ (List.nodup_finRange n)]
  conv_rhs => rw [← PP.restrict_all P, PP.toMatrix_scale]
  rw [Matrix.mul_smul, Matrix.smul_mul]

end Tab

/-! ## The Clifford row updates are pull-backs -/

variable {n : ℕ}

theorem PP.ext' {P Q : PP n} (h1 : P.phase = Q.phase) (h2 : ∀ r, P.letter r = Q.letter r) :
    P = Q := by
  cases P; cases Q; simp only at h1 h2; subst h1; congr; funext r; exact h2 r

theorem PP.single_letter (q r : Fin n) (L : Letter) :
    (PP.single q L).letter r = if r = q then L else .I := by
  simp [PP.single, Function.update_apply]

/-- The product of single-qubit strings on distinct qubits. -/
theorem PP.mul_single_ne {a b : Fin n} (hab : a ≠ b) (L M : Letter) :
    (PP.single a L).mul (PP.single b M) =
      ⟨0, fun r => if r = a then L else if r = b then M else .I⟩ := by
  apply PP.ext'
  · simp only [PP.mul, PP.single_letter]
    simp only [PP.single, zero_add]
    refine Finset.sum_eq_zero fun r _ => ?_
    by_cases ha : r = a
    · subst ha; simp [hab, mulLetter_I_right]
    · by_cases hb : r = b
      · subst hb; simp [ha, mulLetter_I_left]
      · simp [ha, hb, mulLetter]
  · intro r
    simp only [PP.mul, PP.single_letter]
    by_cases ha : r = a
    · subst ha; simp [hab, mulLetter_I_right]
    · by_cases hb : r = b
      · subst hb; simp [ha, mulLetter_I_left]
      · simp [ha, hb, mulLetter]

/-- The product of two single-qubit strings on the same qubit. -/
theorem PP.mul_single_eq (a : Fin n) (L M : Letter) :
    (PP.single a L).mul (PP.single a M) =
      ⟨(mulLetter L M).1, fun r => if r = a then (mulLetter L M).2 else .I⟩ := by
  apply PP.ext'
  · simp only [PP.mul, PP.single_letter]
    simp only [PP.single, zero_add]
    rw [Finset.sum_eq_single a (fun r _ hr => by simp [hr, mulLetter]) (by simp)]
    simp
  · intro r
    simp only [PP.mul, PP.single_letter]
    by_cases ha : r = a <;> simp [ha, mulLetter]

/-- The matrix computations shared by all the row updates. -/
theorem row_mul {A B : PP n} {C A₀ B₀ : Density n} (hCC : C * Cᴴ = 1)
    (hA : A.toMatrix = Cᴴ * A₀ * C) (hB : B.toMatrix = Cᴴ * B₀ * C) :
    (A.mul B).toMatrix = Cᴴ * (A₀ * B₀) * C := by
  rw [PP.toMatrix_mul, hA, hB]
  simp only [Matrix.mul_assoc]
  rw [← Matrix.mul_assoc C, hCC, Matrix.one_mul]

theorem row_scale {A : PP n} {C A₀ : Density n} (k : ZMod 4) (hA : A.toMatrix = Cᴴ * A₀ * C) :
    (A.scale k).toMatrix = Cᴴ * (ipow k • A₀) * C := by
  rw [PP.toMatrix_scale, hA, Matrix.mul_smul, Matrix.smul_mul]

/-! The pulled-back generators `G† X_r G` and `G† Z_r G`, as products of single-qubit strings:
the middle column of the row-update table. -/

section pull

variable (r : Fin n)

local macro "pull_tac" : tactic => `(tactic|
  (apply PP.ext' <;>
    (try intro s) <;>
    simp [pullG, PP.ofS, PP.scale, PP.single, SPauli.conjClifford, CTGate.inv, SPauli.apply1,
      SPauli.apply2, Letter.conjH, Letter.conjS, Letter.conjSdg, Letter.conjX, Letter.conjZ,
      Letter.conjCX, Letter.conjCZ, Function.update_apply, PP.mul_single_ne, PP.mul_single_eq,
      mulLetter, *] <;>
    (try split_ifs) <;> (try simp_all) <;> (try decide)))

theorem pull_h_X (q : Fin n) :
    pullG (.h q) (PP.single r .X) = if r = q then PP.single q .Z else PP.single r .X := by
  split_ifs with hr
  · subst hr; pull_tac
  · pull_tac

theorem pull_h_Z (q : Fin n) :
    pullG (.h q) (PP.single r .Z) = if r = q then PP.single q .X else PP.single r .Z := by
  split_ifs with hr
  · subst hr; pull_tac
  · pull_tac

theorem pull_x_X (q : Fin n) : pullG (.x q) (PP.single r .X) = PP.single r .X := by
  by_cases hr : r = q
  · subst hr; pull_tac
  · pull_tac

theorem pull_x_Z (q : Fin n) :
    pullG (.x q) (PP.single r .Z) = if r = q then (PP.single q .Z).scale 2 else PP.single r .Z := by
  split_ifs with hr
  · subst hr; pull_tac
  · pull_tac

theorem pull_z_X (q : Fin n) :
    pullG (.z q) (PP.single r .X) = if r = q then (PP.single q .X).scale 2 else PP.single r .X := by
  split_ifs with hr
  · subst hr; pull_tac
  · pull_tac

theorem pull_z_Z (q : Fin n) : pullG (.z q) (PP.single r .Z) = PP.single r .Z := by
  by_cases hr : r = q
  · subst hr; pull_tac
  · pull_tac

theorem pull_s_X (q : Fin n) :
    pullG (.s q) (PP.single r .X) =
      if r = q then ((PP.single q .X).mul (PP.single q .Z)).scale 3 else PP.single r .X := by
  split_ifs with hr
  · subst hr; rw [PP.mul_single_eq]; pull_tac
  · pull_tac

theorem pull_s_Z (q : Fin n) : pullG (.s q) (PP.single r .Z) = PP.single r .Z := by
  by_cases hr : r = q
  · subst hr; pull_tac
  · pull_tac

theorem pull_sdg_X (q : Fin n) :
    pullG (.sdg q) (PP.single r .X) =
      if r = q then ((PP.single q .X).mul (PP.single q .Z)).scale 1 else PP.single r .X := by
  split_ifs with hr
  · subst hr; rw [PP.mul_single_eq]; pull_tac
  · pull_tac

theorem pull_sdg_Z (q : Fin n) : pullG (.sdg q) (PP.single r .Z) = PP.single r .Z := by
  by_cases hr : r = q
  · subst hr; pull_tac
  · pull_tac

theorem pull_cx_X (c t : Fin n) (hne : c ≠ t) :
    pullG (.cx c t hne) (PP.single r .X) =
      if r = c then (PP.single c .X).mul (PP.single t .X) else PP.single r .X := by
  split_ifs with hr
  · subst hr; rw [PP.mul_single_ne hne]; pull_tac
  · by_cases ht : r = t
    · subst ht; pull_tac
    · pull_tac

theorem pull_cx_Z (c t : Fin n) (hne : c ≠ t) :
    pullG (.cx c t hne) (PP.single r .Z) =
      if r = t then (PP.single c .Z).mul (PP.single t .Z) else PP.single r .Z := by
  split_ifs with hr
  · subst hr; rw [PP.mul_single_ne hne]; pull_tac
  · by_cases hc : r = c
    · subst hc; pull_tac
    · pull_tac

theorem pull_cz_X (c t : Fin n) (hne : c ≠ t) :
    pullG (.cz c t hne) (PP.single r .X) =
      if r = c then (PP.single c .X).mul (PP.single t .Z)
      else if r = t then (PP.single c .Z).mul (PP.single t .X) else PP.single r .X := by
  split_ifs with hc ht
  · subst hc; rw [PP.mul_single_ne hne]; pull_tac
  · subst ht; rw [PP.mul_single_ne hne]; pull_tac
  · pull_tac

theorem pull_cz_Z (c t : Fin n) (hne : c ≠ t) :
    pullG (.cz c t hne) (PP.single r .Z) = PP.single r .Z := by
  by_cases hc : r = c
  · subst hc; pull_tac
  · by_cases ht : r = t
    · subst ht; pull_tac
    · pull_tac

end pull

namespace Tab

theorem axes_step_clifford {T : Tab n} {g : CTGate n} (hcl : g.IsClifford) :
    (T.step g).axes = T.axes := by
  cases g <;> first | rfl | exact hcl.elim

/-- **The row updates are pull-backs.** For a Clifford gate `g` with matrix `G`, each new row
is the old pull-back of the pulled-back generator: `x'[r] = C† (G† X_r G) C` and
`z'[r] = C† (G† Z_r G) C`. By `back_ok`, the right-hand sides are `β(G† X_r G)` and
`β(G† Z_r G)`. -/
theorem step_rows {T : Tab n} {C : Density n} (hCC : C * Cᴴ = 1) (hrows : T.RowsOK C)
    {g : CTGate n} (hcl : g.IsClifford) (r : Fin n) :
    ((T.step g).rx r).toMatrix = Cᴴ * (pullG g (PP.single r .X)).toMatrix * C ∧
    ((T.step g).rz r).toMatrix = Cᴴ * (pullG g (PP.single r .Z)).toMatrix * C := by
  have hX := fun s => (hrows s).1
  have hZ := fun s => (hrows s).2
  cases g with
  | t q => exact hcl.elim
  | tdg q => exact hcl.elim
  | h q =>
    rw [pull_h_X, pull_h_Z]
    by_cases hr : r = q
    · subst hr; simp [step, hX, hZ]
    · simp [step, hr, hX, hZ]
  | x q =>
    rw [pull_x_X, pull_x_Z]
    by_cases hr : r = q
    · subst hr
      simp only [step, Function.update_self, if_true]
      exact ⟨hX r, (row_scale 2 (hZ r)).trans (by rw [PP.toMatrix_scale])⟩
    · simp [step, hr, hX, hZ]
  | z q =>
    rw [pull_z_X, pull_z_Z]
    by_cases hr : r = q
    · subst hr
      simp only [step, Function.update_self, if_true]
      exact ⟨(row_scale 2 (hX r)).trans (by rw [PP.toMatrix_scale]), hZ r⟩
    · simp [step, hr, hX, hZ]
  | s q =>
    rw [pull_s_X, pull_s_Z]
    by_cases hr : r = q
    · subst hr
      simp only [step, Function.update_self, if_true]
      refine ⟨?_, hZ r⟩
      rw [row_scale 3 (row_mul hCC (hX r) (hZ r)), PP.toMatrix_scale, PP.toMatrix_mul]
    · simp [step, hr, hX, hZ]
  | sdg q =>
    rw [pull_sdg_X, pull_sdg_Z]
    by_cases hr : r = q
    · subst hr
      simp only [step, Function.update_self, if_true]
      refine ⟨?_, hZ r⟩
      rw [row_scale 1 (row_mul hCC (hX r) (hZ r)), PP.toMatrix_scale, PP.toMatrix_mul]
    · simp [step, hr, hX, hZ]
  | cx c t hne =>
    rw [pull_cx_X, pull_cx_Z]
    refine ⟨?_, ?_⟩
    · by_cases hr : r = c
      · subst hr
        simp only [step, Function.update_self, if_true]
        rw [row_mul hCC (hX r) (hX t), PP.toMatrix_mul]
      · simp [step, hr, hX]
    · by_cases hr : r = t
      · subst hr
        simp only [step, Function.update_self, if_true]
        rw [row_mul hCC (hZ c) (hZ r), PP.toMatrix_mul]
      · simp [step, hr, hZ]
  | cz c t hne =>
    rw [pull_cz_X, pull_cz_Z]
    refine ⟨?_, by simp [step, hZ]⟩
    by_cases hc : r = c
    · subst hc
      simp only [step, Function.update_of_ne hne, Function.update_self, if_true]
      rw [row_mul hCC (hX r) (hZ t), PP.toMatrix_mul]
    · by_cases ht : r = t
      · subst ht
        simp only [step, Function.update_self, if_neg hc, if_true]
        rw [row_mul hCC (hZ c) (hX r), PP.toMatrix_mul]
      · simp [step, hc, ht, hX]

/-- A Clifford gate keeps the rows correct, for the Clifford `G C`. -/
theorem rowsOK_step {T : Tab n} {C : Density n} (hCC : C * Cᴴ = 1) (hrows : T.RowsOK C)
    {g : CTGate n} (hcl : g.IsClifford) :
    (T.step g).RowsOK (gateUnitary n g.toGate * C) := fun r => by
  obtain ⟨h1, h2⟩ := step_rows hCC hrows hcl r
  refine ⟨?_, ?_⟩
  · rw [h1, pullG_toMatrix g hcl, Matrix.conjTranspose_mul]; simp only [Matrix.mul_assoc]
  · rw [h2, pullG_toMatrix g hcl, Matrix.conjTranspose_mul]; simp only [Matrix.mul_assoc]

/-- The same, stated with the pull-back map: `x'[r] = β(G† X_r G)`. -/
theorem step_rows_back {T : Tab n} {C : Density n} (hCC : C * Cᴴ = 1) (hrows : T.RowsOK C)
    {g : CTGate n} (hcl : g.IsClifford) (r : Fin n) :
    ((T.step g).rx r).toMatrix = (T.back (pullG g (PP.single r .X))).toMatrix ∧
    ((T.step g).rz r).toMatrix = (T.back (pullG g (PP.single r .Z))).toMatrix := by
  have hC : Cᴴ * C = 1 := mul_eq_one_comm.mp hCC
  rw [back_ok hC hCC hrows, back_ok hC hCC hrows]
  exact step_rows hCC hrows hcl r

end Tab

/-! ## T as a combination of I and Z -/

theorem embed1_diag2 {n : ℕ} (w : ℂ) (q : Fin n) :
    embed1 n (diag2 1 w) q.val =
      ((1 + w) / 2) • (1 : Density n) + ((1 - w) / 2) • embed1 n Letter.Z.mat q.val := by
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

/-- A T or T† is `α·I + β·Z_q`. -/
theorem rotation_decomp {n : ℕ} {g : CTGate n} {q : Fin n} (hg : g = .t q ∨ g = .tdg q) :
    ∃ a b : ℂ, gateUnitary n g.toGate = a • 1 + b • (PP.single q .Z).toMatrix := by
  rw [PP.toMatrix_single]
  rcases hg with rfl | rfl
  · exact ⟨_, _, embed1_diag2 (ep (1 / 4)) q⟩
  · exact ⟨_, _, embed1_diag2 (ep (-1 / 4)) q⟩

/-! ## The invariant and soundness -/

/-- The unitary of a circuit. -/
abbrev circ {n : ℕ} (gs : List (CTGate n)) : Density n := unitary n (gs.map CTGate.toGate)

namespace Tab

variable {n : ℕ}

/-- `U = C · V`: the rows are correct for the Clifford `C`, and `V` commutes with every
Pauli that commutes with all the axes. -/
def Inv (U : Density n) (T : Tab n) : Prop :=
  ∃ C V : Density n, U = C * V ∧ C * Cᴴ = 1 ∧ V * Vᴴ = 1 ∧ T.RowsOK C ∧
    ∀ B : PP n, (∀ A ∈ T.axes, PP.Comm B A) → V * B.toMatrix = B.toMatrix * V

theorem inv_init : Inv 1 (init n) :=
  ⟨1, 1, by simp, by simp, by simp,
    fun r => ⟨by simp [init], by simp [init]⟩, fun _ _ => by simp⟩

theorem inv_step {U : Density n} {T : Tab n} (h : Inv U T) (g : CTGate n) :
    Inv (gateUnitary n g.toGate * U) (T.step g) := by
  obtain ⟨C, V, rfl, hCC, hVV, hrows, hcomm⟩ := h
  have hC : Cᴴ * C = 1 := mul_eq_one_comm.mp hCC
  set G := gateUnitary n g.toGate
  by_cases hT : ∃ q, g = .t q ∨ g = .tdg q
  · -- A rotation: record its axis, keep the Clifford.
    obtain ⟨q, hq⟩ := hT
    obtain ⟨a, b, hG⟩ := rotation_decomp hq
    have hstep : T.step g = { T with axes := T.rz q :: T.axes } := by
      rcases hq with rfl | rfl <;> rfl
    rw [hstep]
    refine ⟨C, Cᴴ * G * C * V, ?_, hCC, ?_, hrows, fun B hB => ?_⟩
    · simp only [← Matrix.mul_assoc]; rw [hCC, Matrix.one_mul]
    · have hGG : G * Gᴴ = 1 := gate_mul_adj g
      simp only [Matrix.conjTranspose_mul, Matrix.conjTranspose_conjTranspose, Matrix.mul_assoc]
      rw [← Matrix.mul_assoc V, hVV, Matrix.one_mul, ← Matrix.mul_assoc C, hCC, Matrix.one_mul,
        ← Matrix.mul_assoc G, hGG, Matrix.one_mul, hC]
    · -- `C† G C = a·I + b·z[q]` commutes with `B`.
      have hA : PP.Comm B (T.rz q) := hB _ (List.mem_cons_self ..)
      have hold : ∀ A ∈ T.axes, PP.Comm B A := fun A hA => hB A (List.mem_cons_of_mem _ hA)
      have hmid : Cᴴ * G * C = a • 1 + b • (T.rz q).toMatrix := by
        rw [(hrows q).2, show G = a • 1 + b • (PP.single q .Z).toMatrix from hG, Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul,
          Matrix.mul_one, hC, Matrix.mul_smul, Matrix.smul_mul]
      have hBA := hA.toMatrix
      rw [hmid, Matrix.mul_assoc _ V, hcomm B hold, ← Matrix.mul_assoc, ← Matrix.mul_assoc]
      congr 1
      rw [Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul, Matrix.smul_mul, Matrix.mul_smul,
        Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, hBA]
  · -- A Clifford: update the rows.
    push Not at hT
    have hcl : g.IsClifford := by
      cases g with
      | t q => exact ((hT q).1 rfl).elim
      | tdg q => exact ((hT q).2 rfl).elim
      | _ => trivial
    have hGG : G * Gᴴ = 1 := gate_mul_adj g
    refine ⟨G * C, V, by rw [Matrix.mul_assoc], ?_, hVV, rowsOK_step hCC hrows hcl, ?_⟩
    · rw [Matrix.conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc C, hCC,
        Matrix.one_mul, hGG]
    · rw [axes_step_clifford hcl]; exact hcomm

theorem inv_run (gs : List (CTGate n)) :
    ∀ {U : Density n} {T : Tab n}, Inv U T → Inv (circ gs * U) (T.run gs) := by
  induction gs with
  | nil => intro U T h; simpa [run, circ] using h
  | cons g gs ih =>
    intro U T h
    have := ih (inv_step h g)
    simp only [run, List.foldl_cons, circ, List.map_cons, unitary_cons] at this ⊢
    rwa [Matrix.mul_assoc]

/-- **Soundness of the tableau abstraction.** After the circuit `gs`, let `B = C† P C` be the
pull-back of a Pauli `P` computed from the rows. If `B` commutes with every recorded axis,
the circuit maps `B` exactly to `P`: `U B U† = P`. -/
theorem tableau_sound (gs : List (CTGate n)) (P : PP n)
    (h : ∀ A ∈ ((init n).run gs).axes, PP.Comm (((init n).run gs).back P) A) :
    conj (circ gs) (((init n).run gs).back P).toMatrix = P.toMatrix := by
  have hI := inv_run gs (inv_init (n := n))
  rw [Matrix.mul_one] at hI
  obtain ⟨C, V, hU, hCC, hVV, hrows, hcomm⟩ := hI
  have hC : Cᴴ * C = 1 := mul_eq_one_comm.mp hCC
  have hB := hcomm _ h
  rw [hU, conj, Matrix.conjTranspose_mul]
  simp only [Matrix.mul_assoc]
  rw [← Matrix.mul_assoc V, hB, Matrix.mul_assoc, ← Matrix.mul_assoc V Vᴴ, hVV, Matrix.one_mul,
    back_ok hC hCC hrows]
  simp only [Matrix.mul_assoc]
  rw [hCC, Matrix.mul_one, ← Matrix.mul_assoc, hCC, Matrix.one_mul]


/-- The tableau concretization of a circuit: the unitaries that map `C† P C` to `P` for every
Pauli `P` whose pull-back commutes with all the recorded axes. -/
def gamma (gs : List (CTGate n)) : Set (Density n) :=
  {U | ∀ P : PP n, (∀ A ∈ ((init n).run gs).axes, PP.Comm (((init n).run gs).back P) A) →
    conj U (((init n).run gs).back P).toMatrix = P.toMatrix}

/-- The circuit's unitary is in its tableau concretization. -/
theorem circ_mem_gamma (gs : List (CTGate n)) : circ gs ∈ gamma gs :=
  fun P h => tableau_sound gs P h

end Tab

/-- Example: on `h; t; h`, the T's axis is `X` in the input frame, and `X` commutes with it,
so the tableau proves `U X U† = X`. -/
example : conj (circ [CTGate.h (0 : Fin 1), .t 0, .h 0]) (PP.single (0 : Fin 1) .X).toMatrix =
    (PP.single (0 : Fin 1) .X).toMatrix := by
  have h := Tab.tableau_sound [CTGate.h (0 : Fin 1), .t 0, .h 0] (PP.single 0 .X) (by decide)
  have hb : ((Tab.init 1).run [CTGate.h (0 : Fin 1), .t 0, .h 0]).back (PP.single 0 .X) =
      PP.single 0 .X := by decide
  rwa [hb] at h

end

end PauliFold
