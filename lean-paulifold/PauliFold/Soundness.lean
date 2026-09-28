import PauliFold.Containment

/-!
# Soundness against the concrete semantics

The concrete transformer of a gate conjugates an operator by the gate's unitary
(`TzapLean.gateUnitary`). This file proves that the Pauli-fold transformer is sound for it:
if an operator lies in the concretization before a gate, its conjugate lies in the
concretization after. In particular every conjugation table in `Pauli.lean` is a proved
identity of matrices, not a numerical fact.

The proof has three layers.

1. *Conjugation lemmas.* Conjugating a Pauli string by `embed1 M q` changes only its factor
   at `q`, into `∑ M L M†`; by a permutation matrix it evaluates the string at permuted basis
   states; by a phase matrix it multiplies entries by phases.
2. *Letter identities.* Finite checks on letters: 2 × 2 identities for the single-qubit
   gates and two-letter identities for CX and CZ.
3. *Assembly.* Soundness of `PauliAbs.step` and `PauliAbs.run`, and, through `containment`,
   of the phase-folding domain.
-/

namespace PauliFold

open TzapLean Matrix
open SPauli (sgn)

variable {n : Nat}

/-! ## Basis states and qubit indices -/

theorem basis_get_val (b : Basis n) (q : Fin n) : b.get q.val = b q := by
  simp [Basis.get]

theorem basis_set_val (b : Basis n) (q : Fin n) (v : Bool) :
    b.set q.val v = Function.update b q v := by
  funext r
  by_cases h : r = q
  · subst h; simp [Basis.set]
  · simp [Basis.set, Function.update, h, Fin.val_inj]

/-! ## Layer 1: conjugation lemmas -/

theorem embed1_val_apply (M : Bool → Bool → ℂ) (q : Fin n) (o a : Basis n) :
    embed1 n M q.val o a =
      if ∀ r : Fin n, r ≠ q → o r = a r then M (o q) (a q) else 0 := by
  rw [embed1_apply_of_lt M q.isLt]
  simp only [Fin.eta, ne_eq, ← Fin.val_inj]

/-- A basis state agreeing with `o` off `q` is `o` updated at `q`. -/
theorem sum_agree (o : Basis n) (q : Fin n) (g : Basis n → Bool → ℂ) :
    (∑ a : Basis n, if ∀ r : Fin n, r ≠ q → o r = a r then g a (a q) else 0) =
      ∑ α : Bool, g (Function.update o q α) α := by
  have key : ∀ a : Basis n,
      (if ∀ r : Fin n, r ≠ q → o r = a r then g a (a q) else 0) =
        ∑ α : Bool, if a = Function.update o q α then g a α else 0 := by
    intro a
    by_cases h : ∀ r : Fin n, r ≠ q → o r = a r
    · have ha : a = Function.update o q (a q) := by
        funext r
        by_cases hr : r = q
        · subst hr; simp
        · simp [Function.update, hr, h r hr]
      rw [if_pos h, Fintype.sum_bool]
      cases hq : a q
      · have h1 : a = Function.update o q false := by rw [hq] at ha; exact ha
        have h2 : a ≠ Function.update o q true := fun h' => by
          have := congrFun h' q; simp [hq] at this
        rw [if_neg h2, if_pos h1, zero_add]
      · have h1 : a ≠ Function.update o q false := fun h' => by
          have := congrFun h' q; simp [hq] at this
        have h2 : a = Function.update o q true := by rw [hq] at ha; exact ha
        rw [if_pos h2, if_neg h1, add_zero]
    · rw [if_neg h]
      symm
      apply Finset.sum_eq_zero
      intro α _
      rw [if_neg]
      rintro rfl
      exact h fun r hr => by simp [Function.update, hr]
  rw [Finset.sum_congr rfl fun a _ => key a, Finset.sum_comm]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Finset.sum_ite_eq']
  simp

theorem embed1_mul_apply (M : Bool → Bool → ℂ) (q : Fin n) (B : Op n) (o i : Basis n) :
    (embed1 n M q.val * B) o i =
      ∑ α : Bool, M (o q) α * B (Function.update o q α) i := by
  rw [Matrix.mul_apply]
  simp only [embed1_val_apply, ite_mul, zero_mul]
  exact sum_agree o q (fun a α => M (o q) α * B a i)

theorem mul_embed1_conjTranspose_apply (M : Bool → Bool → ℂ) (q : Fin n) (B : Op n)
    (o i : Basis n) :
    (B * (embed1 n M q.val)ᴴ) o i =
      ∑ β : Bool, B o (Function.update i q β) * star (M (i q) β) := by
  rw [Matrix.mul_apply]
  simp only [Matrix.conjTranspose_apply, embed1_val_apply]
  have := sum_agree i q (fun b β => B o b * star (M (i q) β))
  rw [← this]
  refine Finset.sum_congr rfl fun b _ => ?_
  split_ifs <;> simp

/-- Split a product over the qubits at `q`. -/
theorem prod_split (F : Fin n → ℂ) (q : Fin n) :
    ∏ r, F r = F q * ∏ r ∈ Finset.univ.erase q, F r :=
  (Finset.mul_prod_erase _ _ (Finset.mem_univ q)).symm

/-- The entries of a Pauli string at states updated at `q`. -/
theorem toMatrix_update (P : SPauli n) (q : Fin n) (o i : Basis n) (α β : Bool) :
    P.toMatrix (Function.update o q α) (Function.update i q β) =
      sgn P.neg * ((P.letter q).mat α β *
        ∏ r ∈ Finset.univ.erase q, (P.letter r).mat (o r) (i r)) := by
  unfold SPauli.toMatrix
  rw [prod_split _ q]
  congr 2
  · simp
  · refine Finset.prod_congr rfl fun r hr => ?_
    have hr' : r ≠ q := Finset.ne_of_mem_erase hr
    simp [Function.update, hr']

/-- **Single-qubit conjugation.** If `M L M† = ±L'` for the letter `L` at `q`, conjugating
the string by `embed1 M q` replaces `L` by `L'` and multiplies the sign. -/
theorem conj_embed1_toMatrix (M : Bool → Bool → ℂ) (q : Fin n) (P : SPauli n) (s : Bool)
    (L' : Letter)
    (hM : ∀ x y : Bool,
      ∑ β : Bool, ∑ α : Bool, M x α * (P.letter q).mat α β * star (M y β) =
        sgn s * L'.mat x y) :
    conj (embed1 n M q.val) P.toMatrix =
      (⟨xor P.neg s, Function.update P.letter q L'⟩ : SPauli n).toMatrix := by
  funext o i
  unfold conj
  rw [mul_embed1_conjTranspose_apply]
  simp only [embed1_mul_apply, toMatrix_update, Finset.sum_mul]
  have hrest : ∀ x y : ℂ, ∀ α β : Bool,
      M (o q) α * (sgn P.neg * ((P.letter q).mat α β * x)) * y =
        sgn P.neg * x * (M (o q) α * (P.letter q).mat α β * y) := by
    intros; ring
  simp only [hrest, ← Finset.mul_sum]
  rw [hM (o q) (i q)]
  unfold SPauli.toMatrix
  rw [prod_split _ q]
  simp only [Function.update_self]
  have hsame : ∏ r ∈ Finset.univ.erase q, (Function.update P.letter q L' r).mat (o r) (i r) =
      ∏ r ∈ Finset.univ.erase q, (P.letter r).mat (o r) (i r) :=
    Finset.prod_congr rfl fun r hr => by
      simp [Function.update, Finset.ne_of_mem_erase hr]
  rw [hsame]
  cases P.neg <;> cases s <;> simp [sgn] <;> ring

theorem permMatrix_mul_apply (σ : Basis n → Basis n) (hσ : Function.Involutive σ)
    (B : Op n) (o i : Basis n) : (permMatrix σ * B) o i = B (σ o) i := by
  rw [Matrix.mul_apply]
  have : ∀ a, (o = σ a) ↔ (a = σ o) := fun a =>
    ⟨fun h => by rw [h, hσ], fun h => by rw [h, hσ]⟩
  simp only [permMatrix, ite_mul, one_mul, zero_mul, this]
  simp

theorem mul_permMatrix_conjTranspose_apply (σ : Basis n → Basis n)
    (hσ : Function.Involutive σ) (B : Op n) (o i : Basis n) :
    (B * (permMatrix σ)ᴴ) o i = B o (σ i) := by
  rw [Matrix.mul_apply]
  have : ∀ b, (i = σ b) ↔ (b = σ i) := fun b =>
    ⟨fun h => by rw [h, hσ], fun h => by rw [h, hσ]⟩
  simp only [Matrix.conjTranspose_apply, permMatrix, this]
  simp

theorem phaseMatrix_mul_apply (f : Basis n → ℂ) (B : Op n) (o i : Basis n) :
    (phaseMatrix f * B) o i = f o * B o i := by
  rw [Matrix.mul_apply]
  simp [phaseMatrix, ite_mul]

theorem mul_phaseMatrix_conjTranspose_apply (f : Basis n → ℂ) (B : Op n) (o i : Basis n) :
    (B * (phaseMatrix f)ᴴ) o i = B o i * star (f i) := by
  rw [Matrix.mul_apply, Finset.sum_eq_single i]
  · simp [phaseMatrix, Matrix.conjTranspose_apply]
  · intro j _ hj
    simp [phaseMatrix, Matrix.conjTranspose_apply, Ne.symm hj]
  · simp

/-- Split a product over the qubits at two distinct qubits. -/
theorem prod_split2 (F : Fin n → ℂ) (c t : Fin n) (hne : c ≠ t) :
    ∏ r, F r = F c * (F t * ∏ r ∈ (Finset.univ.erase c).erase t, F r) := by
  rw [prod_split F c, Finset.mul_prod_erase _ _ (Finset.mem_erase.mpr ⟨Ne.symm hne,
    Finset.mem_univ t⟩)]

/-! ## Layer 2: letter identities

Each identity is a finite check over letters and bits. -/

theorem sqrt2_mul_self : ((Real.sqrt 2 : ℝ) : ℂ) * ((Real.sqrt 2 : ℝ) : ℂ) = 2 := by
  rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by norm_num)]
  norm_num

theorem ep_one : ep 1 = -1 := by
  simp [ep, Complex.exp_pi_mul_I]

theorem ep_half : ep (1 / 2) = Complex.I := by
  have : ((Real.pi * ((1 / 2 : ℚ) : ℝ) : ℝ) : ℂ) * Complex.I = (Real.pi / 2 : ℝ) * Complex.I := by
    push_cast; ring
  rw [ep, this, Complex.exp_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin,
    Real.cos_pi_div_two, Real.sin_pi_div_two]
  simp

theorem ep_neg_half : ep (-1 / 2) = -Complex.I := by
  have h : ep (-1 / 2) = star (ep (1 / 2)) := by rw [star_ep]; norm_num
  rw [h, ep_half]
  simp

/-- The 2 × 2 identity `M L M† = ±L'` behind a single-qubit table entry. -/
def Table1 (M : Bool → Bool → ℂ) (L : Letter) (r : Bool × Letter) : Prop :=
  ∀ x y : Bool, ∑ β : Bool, ∑ α : Bool, M x α * L.mat α β * star (M y β) =
    sgn r.1 * r.2.mat x y

theorem table_x (L : Letter) : Table1 x2 L (Letter.conjX L) := by
  intro x y
  cases L <;> cases x <;> cases y <;> simp [x2, Letter.mat, Letter.conjX, sgn]

theorem table_diag (a b : ℂ) (ha : a * star a = 1) (hb : b * star b = 1) (L : Letter)
    (hL : L = .I ∨ L = .Z) : Table1 (diag2 a b) L (false, L) := by
  intro x y
  have ha' : a * (starRingEnd ℂ) a = 1 := ha
  have hb' : b * (starRingEnd ℂ) b = 1 := hb
  rcases hL with rfl | rfl <;> cases x <;> cases y <;>
    simp [diag2, Letter.mat, sgn, ha', hb']

theorem ep_mul_star' (θ : ℚ) : ep θ * star (ep θ) = 1 := ep_mul_star θ

theorem table_z (L : Letter) : Table1 (diag2 1 (ep 1)) L (Letter.conjZ L) := by
  rw [ep_one]
  intro x y
  cases L <;> cases x <;> cases y <;>
    simp [diag2, Letter.mat, Letter.conjZ, sgn]

theorem table_s (L : Letter) : Table1 (diag2 1 (ep (1 / 2))) L (Letter.conjS L) := by
  rw [ep_half]
  intro x y
  cases L <;> cases x <;> cases y <;>
    simp [diag2, Letter.mat, Letter.conjS, sgn]

theorem table_sdg (L : Letter) : Table1 (diag2 1 (ep (-1 / 2))) L (Letter.conjSdg L) := by
  rw [ep_neg_half]
  intro x y
  cases L <;> cases x <;> cases y <;>
    simp [diag2, Letter.mat, Letter.conjSdg, sgn]

theorem table_h (L : Letter) : Table1 h2 L (Letter.conjH L) := by
  intro x y
  have hs := sqrt2_mul_self
  have h0 : ((Real.sqrt 2 : ℝ) : ℂ) ≠ 0 := by
    intro h; rw [h, zero_mul] at hs; norm_num at hs
  have hs2 : ((Real.sqrt 2 : ℝ) : ℂ) ^ 2 = 2 := by rw [sq]; exact hs
  cases L <;> cases x <;> cases y <;>
    simp [h2, Letter.mat, Letter.conjH, sgn, Complex.conj_ofReal] <;>
    field_simp <;> ring_nf <;> (try simp only [hs2])

/-- The two-letter identity behind a CX table entry, with the control first. -/
theorem table_cx (lc lt : Letter) (oc ic ot it : Bool) :
    lc.mat oc ic * lt.mat (ot != oc) (it != ic) =
      sgn (Letter.conjCX lc lt).1 *
        ((Letter.conjCX lc lt).2.1.mat oc ic * (Letter.conjCX lc lt).2.2.mat ot it) := by
  cases lc <;> cases lt <;> cases oc <;> cases ic <;> cases ot <;> cases it <;>
    simp [Letter.mat, Letter.conjCX, sgn]

/-- The two-letter identity behind a CZ table entry, with the control first. -/
theorem table_cz (lc lt : Letter) (oc ic ot it : Bool) :
    (if oc && ot then -1 else 1) * (lc.mat oc ic * lt.mat ot it) *
        star (if ic && it then (-1 : ℂ) else 1) =
      sgn (Letter.conjCZ lc lt).1 *
        ((Letter.conjCZ lc lt).2.1.mat oc ic * (Letter.conjCZ lc lt).2.2.mat ot it) := by
  cases lc <;> cases lt <;> cases oc <;> cases ic <;> cases ot <;> cases it <;>
    simp [Letter.mat, Letter.conjCZ, sgn]

/-! ## Layer 3: the transformers are sound -/

theorem sgn_xor (a b : Bool) : sgn (xor a b) = sgn a * sgn b := by
  cases a <;> cases b <;> simp [sgn]

/-- A single-qubit table entry, proved as a 2 × 2 identity, gives the exact conjugate. -/
theorem conj_apply1 (M : Bool → Bool → ℂ) (f : Letter → Bool × Letter) (q : Fin n)
    (P : SPauli n) (hM : Table1 M (P.letter q) (f (P.letter q))) :
    conj (embed1 n M q.val) P.toMatrix = (SPauli.apply1 f q P).toMatrix :=
  conj_embed1_toMatrix M q P _ _ hM

/-- CX's basis permutation. -/
def cxPerm (c t : Fin n) (b : Basis n) : Basis n := Function.update b t (b t != b c)

theorem cxPerm_involutive (c t : Fin n) (hne : c ≠ t) : Function.Involutive (cxPerm c t) := by
  intro b
  funext r
  by_cases hr : r = t
  · subst hr; simp [cxPerm, Function.update, hne]
  · simp [cxPerm, Function.update, hr]

theorem gateUnitary_cx (c t : Fin n) (hne : c ≠ t) :
    gateUnitary n (CTGate.cx c t hne).toGate = permMatrix (cxPerm c t) := by
  show permMatrix _ = _
  congr 1
  funext b
  simp [cxPerm, basis_set_val, basis_get_val]

theorem conj_cx (c t : Fin n) (hne : c ≠ t) (P : SPauli n) :
    conj (gateUnitary n (CTGate.cx c t hne).toGate) P.toMatrix =
      (SPauli.apply2 Letter.conjCX c t P).toMatrix := by
  rw [gateUnitary_cx]
  funext o i
  unfold conj
  rw [mul_permMatrix_conjTranspose_apply _ (cxPerm_involutive c t hne),
    permMatrix_mul_apply _ (cxPerm_involutive c t hne)]
  unfold SPauli.toMatrix SPauli.apply2
  rw [prod_split2 _ c t hne, prod_split2 _ c t hne]
  have hrest : ∏ r ∈ (Finset.univ.erase c).erase t,
        (P.letter r).mat (cxPerm c t o r) (cxPerm c t i r) =
      ∏ r ∈ (Finset.univ.erase c).erase t,
        (Function.update (Function.update P.letter c (Letter.conjCX (P.letter c)
          (P.letter t)).2.1) t (Letter.conjCX (P.letter c) (P.letter t)).2.2 r).mat (o r) (i r) := by
    refine Finset.prod_congr rfl fun r hr => ?_
    have hrt : r ≠ t := Finset.ne_of_mem_erase hr
    have hrc : r ≠ c := Finset.ne_of_mem_erase (Finset.mem_of_mem_erase hr)
    simp [cxPerm, Function.update, hrt, hrc]
  rw [hrest]
  simp only [cxPerm, Function.update_self, Function.update_of_ne hne, sgn_xor]
  have key := table_cx (P.letter c) (P.letter t) (o c) (i c) (o t) (i t)
  linear_combination (sgn P.neg * ∏ r ∈ (Finset.univ.erase c).erase t,
    (Function.update (Function.update P.letter c (Letter.conjCX (P.letter c)
      (P.letter t)).2.1) t (Letter.conjCX (P.letter c) (P.letter t)).2.2 r).mat (o r) (i r)) * key

theorem star_pm (p : Prop) [Decidable p] :
    star (if p then (-1 : ℂ) else 1) = if p then -1 else 1 := by
  split_ifs <;> simp

/-- CZ's phase. -/
def czPhase (c t : Fin n) (b : Basis n) : ℂ := if b c && b t then -1 else 1

theorem gateUnitary_cz (c t : Fin n) (hne : c ≠ t) :
    gateUnitary n (CTGate.cz c t hne).toGate = phaseMatrix (czPhase c t) := by
  show phaseMatrix _ = _
  congr 1
  funext b
  simp [czPhase, basis_get_val]

theorem conj_cz (c t : Fin n) (hne : c ≠ t) (P : SPauli n) :
    conj (gateUnitary n (CTGate.cz c t hne).toGate) P.toMatrix =
      (SPauli.apply2 Letter.conjCZ c t P).toMatrix := by
  rw [gateUnitary_cz]
  funext o i
  unfold conj
  rw [mul_phaseMatrix_conjTranspose_apply, phaseMatrix_mul_apply]
  unfold SPauli.toMatrix SPauli.apply2
  rw [prod_split2 _ c t hne, prod_split2 _ c t hne]
  have hrest : ∏ r ∈ (Finset.univ.erase c).erase t, (P.letter r).mat (o r) (i r) =
      ∏ r ∈ (Finset.univ.erase c).erase t,
        (Function.update (Function.update P.letter c (Letter.conjCZ (P.letter c)
          (P.letter t)).2.1) t (Letter.conjCZ (P.letter c) (P.letter t)).2.2 r).mat (o r) (i r) := by
    refine Finset.prod_congr rfl fun r hr => ?_
    have hrt : r ≠ t := Finset.ne_of_mem_erase hr
    have hrc : r ≠ c := Finset.ne_of_mem_erase (Finset.mem_of_mem_erase hr)
    simp [Function.update, hrt, hrc]
  rw [← hrest]
  simp only [czPhase, Function.update_self, Function.update_of_ne hne, sgn_xor,
    star_pm]
  have key := table_cz (P.letter c) (P.letter t) (o c) (i c) (o t) (i t)
  simp only [star_pm] at key
  linear_combination (sgn P.neg * ∏ r ∈ (Finset.univ.erase c).erase t,
    (P.letter r).mat (o r) (i r)) * key

/-- A T or T† fixes a string whose letter at the qubit is `I` or `Z`. -/
theorem conj_diag_fix (θ : ℚ) (q : Fin n) (P : SPauli n)
    (hL : P.letter q = .I ∨ P.letter q = .Z) :
    conj (embed1 n (diag2 1 (ep θ)) q.val) P.toMatrix = P.toMatrix := by
  have h := conj_embed1_toMatrix (diag2 1 (ep θ)) q P false (P.letter q)
    (table_diag 1 (ep θ) (by simp) (ep_mul_star θ) (P.letter q) hL)
  rw [h]
  congr 1
  cases P
  simp

/-- **Soundness of the Pauli transformer.** -/
theorem step_sound (g : CTGate n) (p : PauliAbs n) (A : Op n) (hA : A ∈ p.gamma) :
    g.concrete A ∈ (PauliAbs.step g p).gamma := by
  cases p with
  | top => cases g <;> trivial
  | one P =>
    simp only [PauliAbs.gamma, Set.mem_singleton_iff] at hA
    subst hA
    unfold CTGate.concrete
    cases g with
    | h q => exact conj_apply1 h2 Letter.conjH q P (table_h _)
    | x q => exact conj_apply1 x2 Letter.conjX q P (table_x _)
    | z q => exact conj_apply1 _ Letter.conjZ q P (table_z _)
    | s q => exact conj_apply1 _ Letter.conjS q P (table_s _)
    | sdg q => exact conj_apply1 _ Letter.conjSdg q P (table_sdg _)
    | cx c t hne => exact conj_cx c t hne P
    | cz c t hne => exact conj_cz c t hne P
    | t q | tdg q =>
      simp only [PauliAbs.step]
      split_ifs with hanti
      · trivial
      · have hL : P.letter q = .I ∨ P.letter q = .Z := by
          cases h : P.letter q <;> simp_all [Letter.anticommutesZ]
        exact conj_diag_fix _ q P hL

/-- **Soundness over a segment.** If an operator lies in the Pauli concretization before a
Clifford+T segment, its transport through the segment lies in the concretization after. -/
theorem run_sound (gs : List (CTGate n)) (p : PauliAbs n) (A : Op n) (hA : A ∈ p.gamma) :
    concreteRun gs A ∈ (PauliAbs.run gs p).gamma := by
  induction gs generalizing p A with
  | nil => exact hA
  | cons g gs ih => exact ih _ _ (step_sound g p A hA)

/-- **Soundness of both domains for a candidate rotation.** Transporting the axis `Z_q` of a
rotation on qubit `q` through any Clifford+T segment gives an operator in the Pauli
concretization, and therefore, by `containment`, in the phase-folding concretization. -/
theorem candidate_sound (q : Fin n) (gs : List (CTGate n)) :
    concreteRun gs (zString {q} false).toMatrix ∈
        (PauliAbs.run gs (.one (zString {q} false))).gamma ∧
      concreteRun gs (zString {q} false).toMatrix ∈
        (PhaseAbs.run gs (PhaseAbs.init n q)).gamma := by
  have h := run_sound gs (.one (zString {q} false)) _
    (show _ ∈ (PauliAbs.one (zString {q} false)).gamma from Set.mem_singleton _)
  exact ⟨h, containment q gs h⟩

end PauliFold
