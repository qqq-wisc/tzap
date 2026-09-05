import TzapLean.CnotMin
import TzapLean.Cancel

/-!
# Soundness of the `CnotMin` Analysis

The forward scan of `TzapLean/CnotMin.lean` claims to recover, from a CNOT-dihedral gate
list, the two things that determine it: the affine map each wire ends up holding, and the
phase polynomial. This file proves that claim against the density-matrix semantics, and
draws the conclusion the pass rests on:

> two gate lists whose analyses agree denote the same channel.

Everything is phrased with a *monomial* matrix — a basis permutation composed with a
diagonal phase, `permMatrix σ * phaseMatrix f`. Every interpretable gate is one, monomials
are closed under multiplication, and a block's analysis names exactly the monomial it
denotes (up to a global phase, which the channel semantics discards).
-/

namespace TzapLean

open Matrix

noncomputable section

/-! ## Evaluating parities -/

/-- The value parity `p` takes on the starting values of the block's first `k` wires. -/
def parityValAux {n : Nat} (b : Basis n) (qs : List Qubit) (p : Parity) : Nat → Bool
  | 0 => false
  | i + 1 => xor (if p.testBit i then b.get (qs.getD i 0) else false) (parityValAux b qs p i)

/-- The value parity `p` takes on the block's starting values, read off `b`. -/
def parityVal {n : Nat} (qs : List Qubit) (p : Parity) (b : Basis n) : Bool :=
  parityValAux b qs p qs.length

theorem parityValAux_xor {n : Nat} (b : Basis n) (qs : List Qubit) (p q : Parity) (k : Nat) :
    parityValAux b qs (p ^^^ q) k = xor (parityValAux b qs p k) (parityValAux b qs q k) := by
  induction k with
  | zero => simp [parityValAux]
  | succ k ih =>
      simp only [parityValAux, Nat.testBit_xor, ih]
      cases hp : p.testBit k <;> cases hq : q.testBit k <;>
        cases hv : b.get (qs.getD k 0) <;>
        cases h1 : parityValAux b qs p k <;> cases h2 : parityValAux b qs q k <;> rfl

theorem parityVal_xor {n : Nat} (qs : List Qubit) (p q : Parity) (b : Basis n) :
    parityVal qs (p ^^^ q) b = xor (parityVal qs p b) (parityVal qs q b) :=
  parityValAux_xor b qs p q qs.length

theorem parityValAux_two_pow {n : Nat} (b : Basis n) (qs : List Qubit) (i k : Nat) :
    parityValAux b qs (2 ^ i) k = if i < k then b.get (qs.getD i 0) else false := by
  induction k with
  | zero => simp [parityValAux]
  | succ k ih =>
      simp only [parityValAux, ih, Nat.testBit_two_pow]
      by_cases h : i = k
      · subst h
        simp
      · have h' : ¬ (i < k + 1) ↔ ¬ (i < k) := by omega
        by_cases hlt : i < k
        · simp [h, hlt, Nat.lt_succ_of_lt hlt]
        · have : ¬ i < k + 1 := by omega
          simp [h, hlt, this]

theorem parityVal_two_pow {n : Nat} (qs : List Qubit) (b : Basis n) {i : Nat}
    (hi : i < qs.length) : parityVal qs (2 ^ i) b = b.get (qs.getD i 0) := by
  simp [parityVal, parityValAux_two_pow, hi]

/-! ## The monomial a block denotes -/

/-- The affine map the block applies: each of its wires ends up holding the parity it
tracks, complemented when an odd number of `x` gates landed on it. -/
def blockPerm {n : Nat} (qs : List Qubit) (st : BlockState) (b : Basis n) : Basis n :=
  fun r =>
    match localIdx qs (r : Nat) with
    | some i => xor (parityVal qs (st.parity.getD i 0) b) (st.consts.getD i false)
    | none => b r

/-- The phase polynomial's value on `b`. -/
def blockPhase {n : Nat} (qs : List Qubit) (st : BlockState) (b : Basis n) : ℂ :=
  st.terms.foldr (fun t acc => (if parityVal qs t.1 b then ep t.2 else 1) * acc) 1

/-- The monomial matrix a block denotes. -/
def blockMatrix (n : Nat) (qs : List Qubit) (st : BlockState) : Density n :=
  permMatrix (blockPerm qs st) * phaseMatrix (blockPhase qs st)

/-! ## Monomials multiply -/

theorem phaseMatrix_mul_permMatrix {n : Nat} (f : Basis n → ℂ) (σ : Basis n → Basis n) :
    phaseMatrix f * permMatrix σ = permMatrix σ * phaseMatrix (fun b => f (σ b)) := by
  funext out inp
  have hL : (phaseMatrix f * permMatrix σ) out inp = f out * (if out = σ inp then 1 else 0) := by
    simp only [Matrix.mul_apply, phaseMatrix, permMatrix]
    rw [Finset.sum_eq_single out]
    · simp
    · intro k _ hk; simp [Ne.symm hk]
    · simp
  have hR : (permMatrix σ * phaseMatrix (fun b => f (σ b))) out inp
      = (if out = σ inp then 1 else 0) * f (σ inp) := by
    simp only [Matrix.mul_apply, phaseMatrix, permMatrix]
    rw [Finset.sum_eq_single inp]
    · simp
    · intro k _ hk; simp [hk]
    · simp
  rw [hL, hR]
  by_cases h : out = σ inp
  · subst h; simp
  · simp [h]

/-- Left-multiplying a block's monomial by a permutation. -/
theorem permMatrix_mul_blockMatrix {n : Nat} (σ : Basis n → Basis n) (qs : List Qubit)
    (st : BlockState) :
    permMatrix σ * blockMatrix n qs st =
      permMatrix (σ ∘ blockPerm qs st) * phaseMatrix (blockPhase qs st) := by
  rw [blockMatrix, ← Matrix.mul_assoc, permMatrix_mul]

/-- Left-multiplying a block's monomial by a diagonal. -/
theorem phaseMatrix_mul_blockMatrix {n : Nat} (f : Basis n → ℂ) (qs : List Qubit)
    (st : BlockState) :
    phaseMatrix f * blockMatrix n qs st =
      permMatrix (blockPerm qs st) *
        phaseMatrix (fun b => f (blockPerm qs st b) * blockPhase qs st b) := by
  rw [blockMatrix, ← Matrix.mul_assoc, phaseMatrix_mul_permMatrix, Matrix.mul_assoc,
    phaseMatrix_mul]

/-! ## Reading the block's wires -/

/-- The block's bookkeeping has one entry per wire. -/
structure Shaped (qs : List Qubit) (st : BlockState) : Prop where
  parity : st.parity.length = qs.length
  consts : st.consts.length = qs.length

theorem localIdx_lt {qs : List Qubit} : ∀ {q : Qubit} {i : Nat},
    localIdx qs q = some i → i < qs.length := by
  induction qs with
  | nil => intro q i h; simp [localIdx] at h
  | cons a as ih =>
      intro q i h
      simp only [localIdx, List.findIdx?_cons] at h
      by_cases hq : (a == q) = true
      · rw [if_pos hq] at h
        cases h; simp
      · rw [if_neg (by simpa using hq)] at h
        rcases hmap : as.findIdx? (· == q) with _ | j
        · rw [hmap] at h; simp at h
        · rw [hmap] at h
          simp only [Option.map_some, Option.some.injEq] at h
          subst h
          have := ih (q := q) (i := j) hmap
          simpa using this

theorem localIdx_getD {qs : List Qubit} : ∀ {q : Qubit} {i : Nat},
    localIdx qs q = some i → qs.getD i 0 = q := by
  induction qs with
  | nil => intro q i h; simp [localIdx] at h
  | cons a as ih =>
      intro q i h
      simp only [localIdx, List.findIdx?_cons] at h
      by_cases hq : (a == q) = true
      · rw [if_pos hq] at h
        cases h
        simpa using (by simpa using hq : a = q)
      · rw [if_neg (by simpa using hq)] at h
        rcases hmap : as.findIdx? (· == q) with _ | j
        · rw [hmap] at h; simp at h
        · rw [hmap] at h
          simp only [Option.map_some, Option.some.injEq] at h
          subst h
          simpa using ih (q := q) (i := j) hmap

theorem blockPerm_get {n : Nat} (qs : List Qubit) (st : BlockState) (b : Basis n)
    {q : Qubit} {i : Nat} (hq : q < n) (hi : localIdx qs q = some i) :
    (blockPerm qs st b).get q =
      xor (parityVal qs (st.parity.getD i 0) b) (st.consts.getD i false) := by
  simp only [Basis.get, dif_pos hq, blockPerm]
  simp [hi]

theorem blockPerm_get_of_none {n : Nat} (qs : List Qubit) (st : BlockState) (b : Basis n)
    {q : Qubit} (hq : q < n) (hi : localIdx qs q = none) :
    (blockPerm qs st b).get q = b.get q := by
  simp only [Basis.get, dif_pos hq, blockPerm]
  simp [hi]

/-! ## Folding a rotation into the phase polynomial -/

/-- The phase polynomial's value, as a function of the term list alone. -/
def phaseOf {n : Nat} (qs : List Qubit) (b : Basis n) (ts : List (Parity × ℚ)) : ℂ :=
  ts.foldr (fun t acc => (if parityVal qs t.1 b then ep t.2 else 1) * acc) 1

theorem blockPhase_eq {n : Nat} (qs : List Qubit) (st : BlockState) (b : Basis n) :
    blockPhase qs st b = phaseOf qs b st.terms := rfl

theorem phaseOf_insertTerm {n : Nat} (qs : List Qubit) (b : Basis n) (ts : List (Parity × ℚ))
    (p : Parity) (s : ℚ) :
    phaseOf qs b (BlockState.insertTerm ts p s) =
      phaseOf qs b ts * (if parityVal qs p b then ep s else 1) := by
  induction ts with
  | nil => simp [BlockState.insertTerm, phaseOf]
  | cons t ts ih =>
      obtain ⟨q, a⟩ := t
      simp only [BlockState.insertTerm]
      by_cases hq : (q == p) = true
      · have hqp : q = p := by simpa using hq
        subst hqp
        rw [if_pos hq]
        simp only [phaseOf, List.foldr_cons]
        by_cases hpv : parityVal qs q b
        · rw [if_pos hpv, if_pos hpv, if_pos hpv, ep_add]
          ring
        · simp only [if_neg hpv]
          ring
      · rw [if_neg (by simpa using hq)]
        simp only [phaseOf, List.foldr_cons] at *
        rw [ih]
        ring

theorem blockPhase_addTerm {n : Nat} (qs : List Qubit) (st : BlockState) (p : Parity)
    (k : Bool) (θ : ℚ) (b : Basis n) :
    blockPhase qs (st.addTerm p k θ) b =
      blockPhase qs st b * (if parityVal qs p b then ep (if k then -θ else θ) else 1) := by
  simp only [BlockState.addTerm, blockPhase_eq]
  exact phaseOf_insertTerm qs b st.terms p _

/-- The rotation a gate applies, expressed on the *tracked* parity: complementing the
parity flips the angle and contributes a global phase. -/
theorem rotation_factor {n : Nat} (qs : List Qubit) (st : BlockState) (p : Parity) (k : Bool)
    (θ : ℚ) (b : Basis n) :
    (if xor (parityVal qs p b) k then ep θ else 1) * blockPhase qs st b =
      (if k then ep θ else 1) * blockPhase qs (st.addTerm p k θ) b := by
  rw [blockPhase_addTerm]
  cases k
  · simp only [Bool.xor_false]
    by_cases h : parityVal qs p b <;> simp [h]
    all_goals ring
  · simp only [Bool.xor_true, if_true]
    by_cases h : parityVal qs p b
    · simp only [h, Bool.not_true, if_true]
      rw [show ep θ * (blockPhase qs st b * ep (-θ)) = blockPhase qs st b * (ep θ * ep (-θ)) by ring,
        ← ep_add]
      simp
    · simp only [h, Bool.not_false, if_true]
      simp

/-- The `CZ` phase, written on the three parities the analysis records. -/
theorem cz_phase_identity (u v : Bool) :
    (if u && v then (-1 : ℂ) else 1) =
      (if u then ep (1/2) else 1) * (if v then ep (1/2) else 1) *
        (if xor u v then ep (-1/2) else 1) := by
  have h1 : ep (1/2) = Complex.I := ep_half
  have h2 : ep (-1/2) = -Complex.I := ep_neg_half
  cases u <;> cases v <;>
    simp only [h1, h2, Bool.and_true, Bool.and_false, Bool.xor_true, Bool.xor_false,
      if_true, if_false, Bool.false_eq_true] <;>
    norm_num [Complex.I_mul_I]

/-! ## Each interpretable gate is a monomial -/

theorem phaseMatrix_const_smul {n : Nat} (c : ℂ) (f : Basis n → ℂ) :
    phaseMatrix (fun b => c * f b) = c • phaseMatrix f := by
  funext out inp
  by_cases h : out = inp <;> simp [phaseMatrix, h]

/-- A diagonal gate is the phase `e^{iπθ}` on its wire, up to a global factor. -/
theorem rotAngle_gateUnitary {n : Nat} {g : Gate} {θ : ℚ} {q : Qubit}
    (h : rotAngle g = some (θ, q)) (hq : q < n) :
    ∃ c : ℂ, c * star c = 1 ∧
      gateUnitary n g = c • phaseMatrix (fun s : Basis n => if s.get q then ep θ else 1) := by
  have hbase : ∀ a : ℚ, embed1 n (diag2 1 (ep a)) q =
      phaseMatrix (fun s : Basis n => if s.get q then ep a else 1) := by
    intro a
    rw [embed1_diag2_eq_phaseMatrix _ _ _ hq]
  cases g with
  | t p =>
      simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, rfl⟩ := h
      exact ⟨1, by simp, by simpa [gateUnitary] using hbase (1/4)⟩
  | tdg p =>
      simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, rfl⟩ := h
      exact ⟨1, by simp, by simpa [gateUnitary] using hbase (-1/4)⟩
  | s p =>
      simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, rfl⟩ := h
      exact ⟨1, by simp, by simpa [gateUnitary] using hbase (1/2)⟩
  | sdg p =>
      simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, rfl⟩ := h
      exact ⟨1, by simp, by simpa [gateUnitary] using hbase (-1/2)⟩
  | z p =>
      simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, rfl⟩ := h
      exact ⟨1, by simp, by simpa [gateUnitary] using hbase 1⟩
  | rz a p =>
      simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, rfl⟩ := h
      refine ⟨ep (-a/2), ep_mul_star _, ?_⟩
      have hdiag : diag2 (ep (-a/2)) (ep (a/2)) =
          fun o i => ep (-a/2) * diag2 1 (ep a) o i := by
        funext o i
        have hmul : ep (-a/2) * ep a = ep (a/2) := by
          rw [← ep_add]; congr 1; ring
        cases o <;> cases i <;> simp [diag2, hmul]
      rw [gateUnitary, hdiag, embed1_smul _ _ _ hq, hbase]
  | _ => simp [rotAngle] at h

theorem localIdx_mem {qs : List Qubit} {q : Qubit} {i : Nat} (h : localIdx qs q = some i) :
    q ∈ qs := by
  have hlt := localIdx_lt h
  have hget := localIdx_getD h
  rw [← hget]
  simp [List.getD, List.getElem?_eq_getElem hlt]

/-! ## The analysis is sound -/

theorem feedGate_shaped {qs : List Qubit} {st st' : BlockState} (hsh : Shaped qs st)
    {g : Gate} (h : feedGate qs st g = some st') : Shaped qs st' := by
  unfold feedGate at h
  rcases hrot : rotAngle g with _ | ⟨θ, q⟩
  · rw [hrot] at h
    cases g with
    | x q =>
        simp only at h
        rcases hidx : localIdx qs q with _ | i
        · simp only [hidx] at h; simp at h
        · simp only [hidx] at h
          simp only [Option.some.injEq] at h
          subst h
          exact ⟨hsh.parity, by simpa using hsh.consts⟩
    | cnot c t =>
        simp only at h
        by_cases hct : (c == t) = true
        · simp [hct] at h
        · rw [if_neg (by simpa using hct)] at h
          rcases hc : localIdx qs c with _ | ci <;> rcases ht : localIdx qs t with _ | ti <;>
            simp only [hc, ht] at h
          · simp at h
          · simp at h
          · simp at h
          · simp only [Option.some.injEq] at h
            subst h
            exact ⟨by simpa using hsh.parity, by simpa using hsh.consts⟩
    | cz c t =>
        simp only at h
        by_cases hct : (c == t) = true
        · simp [hct] at h
        · rw [if_neg (by simpa using hct)] at h
          rcases hc : localIdx qs c with _ | ci <;> rcases ht : localIdx qs t with _ | ti <;>
            simp only [hc, ht] at h
          · simp at h
          · simp at h
          · simp at h
          · simp only [Option.some.injEq] at h
            subst h
            exact ⟨hsh.parity, hsh.consts⟩
    | _ => simp at h
  · rw [hrot] at h
    simp only at h
    rcases hidx : localIdx qs q with _ | i
    · simp only [hidx] at h; simp at h
    · simp only [hidx] at h
      rcases hp : st.parity[i]? with _ | p <;> rcases hk : st.consts[i]? with _ | k <;>
        simp only [hp, hk] at h
      · simp at h
      · simp at h
      · simp at h
      · simp only [Option.some.injEq] at h
        subst h
        exact ⟨hsh.parity, hsh.consts⟩

theorem unit_mul {c d : ℂ} (hc : c * star c = 1) (hd : d * star d = 1) :
    (c * d) * star (c * d) = 1 := by
  rw [star_mul']
  calc c * d * (star c * star d) = (c * star c) * (d * star d) := by ring
    _ = 1 := by rw [hc, hd, one_mul]

theorem unit_ite (k : Bool) (θ : ℚ) :
    ((if k then ep θ else 1) : ℂ) * star (if k then ep θ else 1) = 1 := by
  cases k
  · simp
  · simpa using ep_mul_star θ

/-! ### List surgery -/

theorem getD_set_self {α : Type*} (l : List α) {i : Nat} (v d : α) (h : i < l.length) :
    (l.set i v).getD i d = v := by
  rw [List.getD_eq_getElem?_getD, List.getElem?_set_self]
  · rfl
  · exact h

theorem getD_set_ne {α : Type*} (l : List α) {i j : Nat} (h : j ≠ i) (v d : α) :
    (l.set i v).getD j d = l.getD j d := by
  simp [List.getD_eq_getElem?_getD, Ne.symm h]

theorem getElem!_eq_getD {α : Type*} [Inhabited α] (l : List α) (i : Nat) :
    l[i]! = l.getD i default := by
  rw [List.getElem!_eq_getElem?_getD, List.getD_eq_getElem?_getD]

theorem getElem!_bool (l : List Bool) (i : Nat) : l[i]! = l.getD i false := getElem!_eq_getD l i

theorem getElem!_nat (l : List Nat) (i : Nat) : l[i]! = l.getD i 0 := getElem!_eq_getD l i

/-- Distinct wires have distinct local indices. -/
theorem localIdx_inj {qs : List Qubit} {q r : Qubit} {i j : Nat}
    (hq : localIdx qs q = some i) (hr : localIdx qs r = some j) (hne : r ≠ q) : j ≠ i := by
  intro hji
  subst hji
  exact hne ((localIdx_getD hr).symm.trans (localIdx_getD hq))

/-- The `x`, `cnot` and `cz` cases all leave the phase polynomial's *shape* and the
parities of untouched wires alone; this is the pointwise reading of that. -/
theorem blockPerm_other {n : Nat} (qs : List Qubit) (st st' : BlockState) (b : Basis n)
    (r : Fin n) {q : Qubit} {i : Nat} (hidx : localIdx qs q = some i) (hr : (r : Nat) ≠ q)
    (hpar : ∀ j, j ≠ i → st'.parity.getD j 0 = st.parity.getD j 0)
    (hcon : ∀ j, j ≠ i → st'.consts.getD j false = st.consts.getD j false) :
    blockPerm qs st' b r = blockPerm qs st b r := by
  simp only [blockPerm]
  rcases hj : localIdx qs (r : Nat) with _ | j
  · rfl
  · have hji : j ≠ i := localIdx_inj hidx hj hr
    dsimp only
    rw [hpar j hji, hcon j hji]

/-- The `CZ` phase, folded into the three terms the analysis records. -/
theorem cz_factor {n : Nat} (qs : List Qubit) (st : BlockState) (pc pt : Parity)
    (kc kt : Bool) (b : Basis n) :
    (if (xor (parityVal qs pc b) kc && xor (parityVal qs pt b) kt) then (-1 : ℂ) else 1)
        * blockPhase qs st b
      = ((if kc then ep (1/2) else 1) * (if kt then ep (1/2) else 1) *
          (if xor kc kt then ep (-1/2) else 1))
        * blockPhase qs (((st.addTerm pc kc (1/2)).addTerm pt kt (1/2)).addTerm
            (pc ^^^ pt) (xor kc kt) (-1/2)) b := by
  set u := xor (parityVal qs pc b) kc with hu
  set v := xor (parityVal qs pt b) kt with hv
  set st₁ := st.addTerm pc kc (1/2) with hst₁
  set st₂ := st₁.addTerm pt kt (1/2) with hst₂
  have hxor : xor u v = xor (parityVal qs (pc ^^^ pt) b) (xor kc kt) := by
    rw [parityVal_xor, hu, hv]
    cases parityVal qs pc b <;> cases parityVal qs pt b <;> cases kc <;> cases kt <;> rfl
  have h1 := rotation_factor qs st pc kc (1/2) b
  have h2 := rotation_factor qs st₁ pt kt (1/2) b
  have h3 := rotation_factor qs st₂ (pc ^^^ pt) (xor kc kt) (-1/2) b
  rw [← hu] at h1
  rw [← hv] at h2
  rw [← hxor] at h3
  rw [cz_phase_identity u v]
  calc (if u then ep (1/2) else 1) * (if v then ep (1/2) else 1) *
          (if xor u v then ep (-1/2) else 1) * blockPhase qs st b
      = (if v then ep (1/2) else 1) * ((if xor u v then ep (-1/2) else 1) *
          ((if u then ep (1/2) else 1) * blockPhase qs st b)) := by ring
    _ = (if v then ep (1/2) else 1) * ((if xor u v then ep (-1/2) else 1) *
          ((if kc then ep (1/2) else 1) * blockPhase qs st₁ b)) := by rw [h1]
    _ = (if kc then ep (1/2) else 1) * ((if xor u v then ep (-1/2) else 1) *
          ((if v then ep (1/2) else 1) * blockPhase qs st₁ b)) := by ring
    _ = (if kc then ep (1/2) else 1) * ((if xor u v then ep (-1/2) else 1) *
          ((if kt then ep (1/2) else 1) * blockPhase qs st₂ b)) := by rw [h2]
    _ = (if kc then ep (1/2) else 1) * ((if kt then ep (1/2) else 1) *
          ((if xor u v then ep (-1/2) else 1) * blockPhase qs st₂ b)) := by ring
    _ = (if kc then ep (1/2) else 1) * ((if kt then ep (1/2) else 1) *
          ((if xor kc kt then ep (-1/2) else 1) *
            blockPhase qs (st₂.addTerm (pc ^^^ pt) (xor kc kt) (-1/2)) b)) := by rw [h3]
    _ = _ := by ring

/-- **The analysis is sound**: feeding a gate multiplies the block's monomial by that gate's
matrix, up to a global phase. -/
theorem feedGate_sound {n : Nat} {qs : List Qubit} (hrange : ∀ q ∈ qs, q < n)
    {st st' : BlockState} (hsh : Shaped qs st) {g : Gate} (h : feedGate qs st g = some st') :
    ∃ c : ℂ, c * star c = 1 ∧
      gateUnitary n g * blockMatrix n qs st = c • blockMatrix n qs st' := by
  unfold feedGate at h
  rcases hrot : rotAngle g with _ | ⟨θ, q⟩
  · rw [hrot] at h
    cases g with
    | x q =>
        simp only at h
        rcases hidx : localIdx qs q with _ | i
        · simp only [hidx] at h; simp at h
        · simp only [hidx] at h
          simp only [Option.some.injEq] at h
          have hq : q < n := hrange q (localIdx_mem hidx)
          have hilt : i < st.consts.length := by
            rw [hsh.consts]; exact localIdx_lt hidx
          subst h
          refine ⟨1, by simp, ?_⟩
          have hperm : (fun b : Basis n => b.set q (!b.get q)) ∘ blockPerm qs st
              = blockPerm qs { st with consts := st.consts.set i (!st.consts[i]!) } := by
            funext b r
            by_cases hr : (r : Nat) = q
            · have hrq : r = ⟨q, hq⟩ := Fin.ext hr
              have hget : (blockPerm qs st b).get q = (blockPerm qs st b) r := by
                rw [hrq]; simp [Basis.get, hq]
              simp only [Function.comp_apply, Basis.set, if_pos hr, hget]
              simp only [blockPerm, hr, hidx]
              rw [getD_set_self _ _ _ hilt, getElem!_bool]
              cases parityVal qs (st.parity.getD i 0) b <;>
                cases st.consts.getD i false <;> rfl
            · simp only [Function.comp_apply, Basis.set, if_neg hr]
              exact (blockPerm_other qs st
                { st with consts := st.consts.set i (!st.consts[i]!) } b r hidx hr
                (fun j _ => rfl) (fun j hj => getD_set_ne _ hj _ _)).symm
          rw [gateUnitary, embed1_x2_eq_perm, permMatrix_mul_blockMatrix, one_smul, blockMatrix,
            hperm]
          rfl
    | cnot c t =>
        simp only at h
        by_cases hct : (c == t) = true
        · simp [hct] at h
        · rw [if_neg (by simpa using hct)] at h
          rcases hc : localIdx qs c with _ | ci <;> rcases ht : localIdx qs t with _ | ti <;>
            simp only [hc, ht] at h
          · simp at h
          · simp at h
          · simp at h
          · simp only [Option.some.injEq] at h
            have hqc : c < n := hrange c (localIdx_mem hc)
            have hqt : t < n := hrange t (localIdx_mem ht)
            have hct' : c ≠ t := by simpa using hct
            have hplt : ti < st.parity.length := by rw [hsh.parity]; exact localIdx_lt ht
            have hclt : ti < st.consts.length := by rw [hsh.consts]; exact localIdx_lt ht
            subst h
            refine ⟨1, by simp, ?_⟩
            have hperm : (fun b : Basis n => b.set t (xor (b.get t) (b.get c))) ∘ blockPerm qs st
                = blockPerm qs { st with
                    parity := st.parity.set ti (st.parity[ti]! ^^^ st.parity[ci]!)
                    consts := st.consts.set ti (xor (st.consts[ti]!) (st.consts[ci]!)) } := by
              funext b r
              by_cases hr : (r : Nat) = t
              · have hrq : r = ⟨t, hqt⟩ := Fin.ext hr
                have hgt : (blockPerm qs st b).get t = (blockPerm qs st b) r := by
                  rw [hrq]; simp [Basis.get, hqt]
                simp only [Function.comp_apply, Basis.set, if_pos hr, hgt]
                rw [blockPerm_get qs st b hqc hc]
                simp only [blockPerm, hr, ht]
                rw [getD_set_self _ _ _ hplt, getD_set_self _ _ _ hclt, getElem!_nat,
                  getElem!_nat, getElem!_bool, getElem!_bool]
                rw [parityVal_xor]
                cases parityVal qs (st.parity.getD ti 0) b <;>
                  cases parityVal qs (st.parity.getD ci 0) b <;>
                  cases st.consts.getD ti false <;>
                  cases st.consts.getD ci false <;> rfl
              · simp only [Function.comp_apply, Basis.set, if_neg hr]
                exact (blockPerm_other qs st
                  { st with
                    parity := st.parity.set ti (st.parity[ti]! ^^^ st.parity[ci]!)
                    consts := st.consts.set ti (xor (st.consts[ti]!) (st.consts[ci]!)) } b r ht hr
                  (fun j hj => getD_set_ne _ hj _ _) (fun j hj => getD_set_ne _ hj _ _)).symm
            rw [gateUnitary, permMatrix_mul_blockMatrix, one_smul, blockMatrix]
            rw [show (fun b : Basis n => b.set t (b.get t != b.get c)) =
                (fun b : Basis n => b.set t (xor (b.get t) (b.get c))) from rfl, hperm]
            rfl
    | cz c t =>
        simp only at h
        by_cases hct : (c == t) = true
        · simp [hct] at h
        · rw [if_neg (by simpa using hct)] at h
          rcases hc : localIdx qs c with _ | ci <;> rcases ht : localIdx qs t with _ | ti <;>
            simp only [hc, ht] at h
          · simp at h
          · simp at h
          · simp at h
          · simp only [Option.some.injEq] at h
            have hqc : c < n := hrange c (localIdx_mem hc)
            have hqt : t < n := hrange t (localIdx_mem ht)
            have hpc : st.parity.getD ci 0 = st.parity[ci]! := (getElem!_nat _ _).symm
            have hpt : st.parity.getD ti 0 = st.parity[ti]! := (getElem!_nat _ _).symm
            have hkc : st.consts.getD ci false = st.consts[ci]! := (getElem!_bool _ _).symm
            have hkt : st.consts.getD ti false = st.consts[ti]! := (getElem!_bool _ _).symm
            subst h
            refine ⟨(if st.consts[ci]! then ep (1/2) else 1) *
                (if st.consts[ti]! then ep (1/2) else 1) *
                (if xor (st.consts[ci]!) (st.consts[ti]!) then ep (-1/2) else 1),
              unit_mul (unit_mul (unit_ite _ _) (unit_ite _ _)) (unit_ite _ _), ?_⟩
            rw [gateUnitary, phaseMatrix_mul_blockMatrix]
            have hfun : (fun b : Basis n =>
                  (if ((blockPerm qs st b).get c && (blockPerm qs st b).get t) then (-1 : ℂ)
                    else 1) * blockPhase qs st b)
                = fun b => ((if st.consts[ci]! then ep (1/2) else 1) *
                    (if st.consts[ti]! then ep (1/2) else 1) *
                    (if xor (st.consts[ci]!) (st.consts[ti]!) then ep (-1/2) else 1)) *
                  blockPhase qs (((st.addTerm st.parity[ci]! st.consts[ci]! (1/2)).addTerm
                    st.parity[ti]! st.consts[ti]! (1/2)).addTerm
                      (st.parity[ci]! ^^^ st.parity[ti]!)
                      (xor (st.consts[ci]!) (st.consts[ti]!)) (-1/2)) b := by
              funext b
              rw [blockPerm_get qs st b hqc hc, blockPerm_get qs st b hqt ht, hpc, hpt, hkc, hkt]
              exact cz_factor qs st _ _ _ _ b
            rw [hfun, phaseMatrix_const_smul, Matrix.mul_smul, blockMatrix]
            rfl
    | _ => simp at h
  · rw [hrot] at h
    simp only at h
    rcases hidx : localIdx qs q with _ | i
    · simp only [hidx] at h; simp at h
    · simp only [hidx] at h
      rcases hp : st.parity[i]? with _ | p <;> rcases hk : st.consts[i]? with _ | k <;>
        simp only [hp, hk] at h
      · simp at h
      · simp at h
      · simp at h
      · simp only [Option.some.injEq] at h
        have hq : q < n := hrange q (localIdx_mem hidx)
        have hpD : st.parity.getD i 0 = p := by
          rw [List.getD_eq_getElem?_getD, hp]; rfl
        have hkD : st.consts.getD i false = k := by
          rw [List.getD_eq_getElem?_getD, hk]; rfl
        obtain ⟨c₀, hc₀, hgu⟩ := rotAngle_gateUnitary hrot hq
        subst h
        refine ⟨c₀ * (if k then ep θ else 1), unit_mul hc₀ (unit_ite k θ), ?_⟩
        rw [hgu, Matrix.smul_mul, phaseMatrix_mul_blockMatrix]
        have hfun : (fun b : Basis n =>
              (if (blockPerm qs st b).get q then ep θ else 1) * blockPhase qs st b)
            = fun b => (if k then ep θ else 1) * blockPhase qs (st.addTerm p k θ) b := by
          funext b
          rw [blockPerm_get qs st b hq hidx, hpD, hkD]
          exact rotation_factor qs st p k θ b
        rw [hfun, phaseMatrix_const_smul, Matrix.mul_smul, smul_smul, blockMatrix]
        rfl

/-! ## From gates to blocks -/

theorem feedGate_isUnitary {qs : List Qubit} {st st' : BlockState} {g : Gate}
    (h : feedGate qs st g = some st') : g.isUnitary = true := by
  cases g <;> first | rfl | (simp [feedGate, rotAngle] at h)

theorem blockMatrix_initial {n : Nat} (qs : List Qubit) (_hrange : ∀ q ∈ qs, q < n) :
    blockMatrix n qs (BlockState.initial qs.length) = 1 := by
  have hperm : blockPerm qs (BlockState.initial qs.length) = fun b : Basis n => b := by
    funext b r
    simp only [blockPerm]
    rcases hj : localIdx qs (r : Nat) with _ | j
    · rfl
    · have hlt : j < qs.length := localIdx_lt hj
      have hpar : (BlockState.initial qs.length).parity.getD j 0 = 2 ^ j := by
        simp only [BlockState.initial]
        rw [List.getD_eq_getElem?_getD, List.getElem?_map, List.getElem?_range hlt]
        rfl
      have hcon : (BlockState.initial qs.length).consts.getD j false = false := by
        simp only [BlockState.initial, List.getD_eq_getElem?_getD]
        rcases h : (List.replicate qs.length false)[j]? with _ | v
        · rfl
        · have hv : v = false := by
            have := List.getElem?_replicate (a := false) (n := qs.length) (i := j)
            rw [this] at h
            split at h <;> simp_all
          simp [hv]
      dsimp only
      rw [hpar, hcon, parityVal_two_pow qs b hlt, localIdx_getD hj]
      have : (r : Nat) < n := r.isLt
      simp [Basis.get, this]
  have hphase : blockPhase qs (BlockState.initial qs.length) = fun _ : Basis n => (1 : ℂ) := by
    funext b; rfl
  rw [blockMatrix, hperm, hphase, permMatrix_id, phaseMatrix_one, Matrix.one_mul]

theorem foldl_bind_none (qs : List Qubit) :
    ∀ (gs : List Gate), gs.foldl (fun s g => s.bind fun s => feedGate qs s g) none = none := by
  intro gs
  induction gs with
  | nil => rfl
  | cons g gs ih => simpa using ih

/-- **Soundness of the analysis over a whole block.** -/
theorem analyzeFold_sound {n : Nat} {qs : List Qubit} (hrange : ∀ q ∈ qs, q < n) :
    ∀ (gs : List Gate) (st st' : BlockState), Shaped qs st →
      gs.foldl (fun s g => s.bind fun s => feedGate qs s g) (some st) = some st' →
      Shaped qs st' ∧ (∀ g ∈ gs, g.isUnitary = true) ∧
        ∃ c : ℂ, c * star c = 1 ∧
          unitary n gs * blockMatrix n qs st = c • blockMatrix n qs st' := by
  intro gs
  induction gs with
  | nil =>
      intro st st' hsh h
      simp only [List.foldl_nil, Option.some.injEq] at h
      subst h
      exact ⟨hsh, by simp, 1, by simp, by simp⟩
  | cons g gs ih =>
      intro st st' hsh h
      simp only [List.foldl_cons, Option.bind_some] at h
      rcases hfeed : feedGate qs st g with _ | st₁
      · rw [hfeed, foldl_bind_none] at h
        exact absurd h (by simp)
      · rw [hfeed] at h
        obtain ⟨hsh', hu', c', hc', heq'⟩ := ih st₁ st' (feedGate_shaped hsh hfeed) h
        obtain ⟨c, hc, heq⟩ := feedGate_sound hrange hsh hfeed
        refine ⟨hsh', ?_, c * c', unit_mul hc hc', ?_⟩
        · intro g' hg'
          rcases List.mem_cons.mp hg' with rfl | hg'
          · exact feedGate_isUnitary hfeed
          · exact hu' g' hg'
        · rw [unitary_cons, Matrix.mul_assoc, heq, Matrix.mul_smul, heq', smul_smul]

/-- The analysis names the block's monomial, up to a global phase. -/
theorem analyzeGates_sound {n : Nat} {qs : List Qubit} (hrange : ∀ q ∈ qs, q < n)
    {gs : List Gate} {st : BlockState} (h : analyzeGates qs gs = some st) :
    (∀ g ∈ gs, g.isUnitary = true) ∧
      ∃ c : ℂ, c * star c = 1 ∧ unitary n gs = c • blockMatrix n qs st := by
  obtain ⟨-, hu, c, hc, heq⟩ :=
    analyzeFold_sound hrange gs (BlockState.initial qs.length) st
      ⟨by simp [BlockState.initial], by simp [BlockState.initial]⟩ h
  refine ⟨hu, c, hc, ?_⟩
  rwa [blockMatrix_initial qs hrange, Matrix.mul_one] at heq

/-! ## Normalization does not move the block -/

theorem ep_of_angleIsZero {θ : ℚ} (h : BlockState.angleIsZero θ = true) : ep θ = 1 := by
  simp only [BlockState.angleIsZero, beq_iff_eq, sub_eq_zero] at h
  rw [h]
  simpa using ep_two_mul_int ⌊θ / 2⌋

theorem ep_angleMod (θ : ℚ) : ep (BlockState.angleMod θ) = ep θ := by
  simp only [BlockState.angleMod]
  have : θ - 2 * (⌊θ / 2⌋ : ℚ) = θ + 2 * (-⌊θ / 2⌋ : ℤ) := by push_cast; ring
  rw [this, ep_add]
  rw [ep_two_mul_int (-⌊θ / 2⌋), mul_one]

theorem phaseOf_cons {n : Nat} (qs : List Qubit) (b : Basis n) (t : Parity × ℚ)
    (ts : List (Parity × ℚ)) :
    phaseOf qs b (t :: ts) = (if parityVal qs t.1 b then ep t.2 else 1) * phaseOf qs b ts := rfl

theorem phaseOf_eq_prod {n : Nat} (qs : List Qubit) (b : Basis n) (ts : List (Parity × ℚ)) :
    phaseOf qs b ts = (ts.map fun t => if parityVal qs t.1 b then ep t.2 else 1).prod := by
  induction ts with
  | nil => rfl
  | cons t ts ih => rw [phaseOf_cons, ih, List.map_cons, List.prod_cons]

theorem phaseOf_perm {n : Nat} (qs : List Qubit) (b : Basis n) {ts us : List (Parity × ℚ)}
    (h : ts.Perm us) : phaseOf qs b ts = phaseOf qs b us := by
  rw [phaseOf_eq_prod, phaseOf_eq_prod]
  exact (h.map _).prod_eq

theorem phaseOf_filter {n : Nat} (qs : List Qubit) (b : Basis n) (ts : List (Parity × ℚ)) :
    phaseOf qs b (ts.filter fun t => !BlockState.angleIsZero t.2) = phaseOf qs b ts := by
  induction ts with
  | nil => rfl
  | cons t ts ih =>
      rw [List.filter_cons]
      by_cases h : BlockState.angleIsZero t.2 = true
      · rw [if_neg (by simp [h]), ih, phaseOf_cons, ep_of_angleIsZero h]
        simp
      · rw [if_pos (by simp [h]), phaseOf_cons, phaseOf_cons, ih]

theorem phaseOf_map_mod {n : Nat} (qs : List Qubit) (b : Basis n) (ts : List (Parity × ℚ)) :
    phaseOf qs b (ts.map fun t => (t.1, BlockState.angleMod t.2)) = phaseOf qs b ts := by
  induction ts with
  | nil => rfl
  | cons t ts ih =>
      rw [List.map_cons, phaseOf_cons, phaseOf_cons, ih, ep_angleMod]

/-- Normalization changes the term list but not the operator it denotes. -/
theorem blockMatrix_normalize {n : Nat} (qs : List Qubit) (st : BlockState) :
    blockMatrix n qs st.normalize = blockMatrix n qs st := by
  have hperm : blockPerm qs st.normalize = blockPerm (n := n) qs st := rfl
  have hphase : blockPhase (n := n) qs st.normalize = blockPhase qs st := by
    funext b
    simp only [blockPhase_eq, BlockState.normalize]
    rw [phaseOf_perm qs b (List.mergeSort_perm _ _), phaseOf_map_mod, phaseOf_filter]
  rw [blockMatrix, blockMatrix, hperm, hphase]

/-! ## The conclusion the pass rests on -/

/-- **Blocks with the same analysis are the same channel.** -/
theorem equivalent_of_analyze {n m : Nat} {qs : List Qubit} (hrange : ∀ q ∈ qs, q < n)
    {gs hs : List Gate} {d₁ d₂ : BlockState}
    (h₁ : analyzeGates qs gs = some d₁) (h₂ : analyzeGates qs hs = some d₂)
    (heq : d₁.normalize = d₂.normalize) : Equivalent n m gs hs := by
  obtain ⟨hu₁, c₁, hc₁, he₁⟩ := analyzeGates_sound hrange h₁
  obtain ⟨hu₂, c₂, hc₂, he₂⟩ := analyzeGates_sound hrange h₂
  have hB : blockMatrix n qs d₁ = blockMatrix n qs d₂ := by
    rw [← blockMatrix_normalize qs d₁, ← blockMatrix_normalize qs d₂, heq]
  have hc₂' : c₂ ≠ 0 := by
    intro hzero
    rw [hzero] at hc₂
    simp at hc₂
  refine equivalent_of_unitary_smul hu₁ hu₂ (c₁ * star c₂) (unit_mul hc₁ ?_) ?_
  · rw [star_star, mul_comm]
    exact hc₂
  · rw [he₁, he₂, hB, smul_smul]
    congr 1
    calc c₁ = c₁ * (c₂ * star c₂) := by rw [hc₂, mul_one]
      _ = c₁ * star c₂ * c₂ := by ring

/-! ## Well-formedness from the certificate -/

theorem feedGate_wf {qs : List Qubit} {st st' : BlockState} {g : Gate}
    (h : feedGate qs st g = some st') : g.Wf := by
  cases g with
  | cnot c t =>
      simp only [feedGate, rotAngle] at h
      by_cases hct : (c == t) = true
      · simp [hct] at h
      · simpa [Gate.Wf] using (by simpa using hct : ¬ c = t)
  | cz c t =>
      simp only [feedGate, rotAngle] at h
      by_cases hct : (c == t) = true
      · simp [hct] at h
      · simpa [Gate.Wf] using (by simpa using hct : ¬ c = t)
  | _ => trivial

theorem analyzeFold_wf (qs : List Qubit) :
    ∀ (gs : List Gate) (st st' : BlockState),
      gs.foldl (fun s g => s.bind fun s => feedGate qs s g) (some st) = some st' →
      ∀ g ∈ gs, g.Wf := by
  intro gs
  induction gs with
  | nil => intro st st' _ g hg; simp at hg
  | cons g gs ih =>
      intro st st' h g' hg'
      simp only [List.foldl_cons, Option.bind_some] at h
      rcases hfeed : feedGate qs st g with _ | st₁
      · rw [hfeed, foldl_bind_none] at h; exact absurd h (by simp)
      · rw [hfeed] at h
        rcases List.mem_cons.mp hg' with rfl | hg'
        · exact feedGate_wf hfeed
        · exact ih st₁ st' h g' hg'

theorem analyzeGates_wf {qs : List Qubit} {gs : List Gate} {d : BlockState}
    (h : analyzeGates qs gs = some d) : ∀ g ∈ gs, g.Wf :=
  analyzeFold_wf qs gs _ d h

/-! ## Operand ranges from the certificate

The `Wf` argument above and this one both read off the same fact: `feedGate` returns `none`
for anything it cannot place on the block's wires, so a gate list the analysis *accepts* has
every operand in `ch.qubits` — and `Chunk.Ok` says those are wires of the register. That
covers the synthesized replacement, which is the only gate list `CnotMin` invents. -/

theorem rotAngle_shape {g : Gate} {θ : ℚ} {q : Qubit} (h : rotAngle g = some (θ, q)) :
    g.qubitsOf = [q] ∧ g.cbitsOf = [] := by
  cases g <;> simp_all [rotAngle, Gate.qubitsOf, Gate.cbitsOf]

/-- A gate the block accepted lives on the block's wires and writes no classical bit. -/
theorem feedGate_mem {qs : List Qubit} {st st' : BlockState} {g : Gate}
    (h : feedGate qs st g = some st') : (∀ q ∈ g.qubitsOf, q ∈ qs) ∧ g.cbitsOf = [] := by
  unfold feedGate at h
  rcases hrot : rotAngle g with _ | ⟨θ, q⟩
  · rw [hrot] at h
    cases g with
    | x q =>
        simp only at h
        rcases hidx : localIdx qs q with _ | i
        · simp only [hidx] at h; simp at h
        · refine ⟨fun r hr => ?_, rfl⟩
          simp only [Gate.qubitsOf, List.mem_singleton] at hr
          exact hr ▸ localIdx_mem hidx
    | cnot c t =>
        simp only at h
        by_cases hct : (c == t) = true
        · simp [hct] at h
        · rw [if_neg (by simpa using hct)] at h
          rcases hc : localIdx qs c with _ | ci <;> rcases ht : localIdx qs t with _ | ti <;>
            simp only [hc, ht] at h
          · simp at h
          · simp at h
          · simp at h
          · refine ⟨fun r hr => ?_, rfl⟩
            simp only [Gate.qubitsOf, List.mem_cons,
              List.not_mem_nil, or_false] at hr
            rcases hr with rfl | rfl
            · exact localIdx_mem hc
            · exact localIdx_mem ht
    | cz c t =>
        simp only at h
        by_cases hct : (c == t) = true
        · simp [hct] at h
        · rw [if_neg (by simpa using hct)] at h
          rcases hc : localIdx qs c with _ | ci <;> rcases ht : localIdx qs t with _ | ti <;>
            simp only [hc, ht] at h
          · simp at h
          · simp at h
          · simp at h
          · refine ⟨fun r hr => ?_, rfl⟩
            simp only [Gate.qubitsOf, List.mem_cons,
              List.not_mem_nil, or_false] at hr
            rcases hr with rfl | rfl
            · exact localIdx_mem hc
            · exact localIdx_mem ht
    | _ => simp at h
  · rw [hrot] at h
    simp only at h
    rcases hidx : localIdx qs q with _ | i
    · simp only [hidx] at h; simp at h
    · obtain ⟨hq, hc⟩ := rotAngle_shape hrot
      refine ⟨fun r hr => ?_, hc⟩
      rw [hq, List.mem_singleton] at hr
      exact hr ▸ localIdx_mem hidx

theorem feedGate_inRange {n m : Nat} {qs : List Qubit} (hqs : ∀ q ∈ qs, q < n)
    {st st' : BlockState} {g : Gate} (h : feedGate qs st g = some st') : g.InRange n m := by
  obtain ⟨hq, hc⟩ := feedGate_mem h
  exact ⟨fun r hr => hqs r (hq r hr), fun b hb => absurd hb (by rw [hc]; simp)⟩

theorem analyzeFold_inRange {n m : Nat} (qs : List Qubit) (hqs : ∀ q ∈ qs, q < n) :
    ∀ (gs : List Gate) (st st' : BlockState),
      gs.foldl (fun s g => s.bind fun s => feedGate qs s g) (some st) = some st' →
      ∀ g ∈ gs, g.InRange n m := by
  intro gs
  induction gs with
  | nil => intro st st' _ g hg; simp at hg
  | cons g gs ih =>
      intro st st' h g' hg'
      simp only [List.foldl_cons, Option.bind_some] at h
      rcases hfeed : feedGate qs st g with _ | st₁
      · rw [hfeed, foldl_bind_none] at h; exact absurd h (by simp)
      · rw [hfeed] at h
        rcases List.mem_cons.mp hg' with rfl | hg'
        · exact feedGate_inRange hqs hfeed
        · exact ih st₁ st' h g' hg'

theorem analyzeGates_inRange {n m : Nat} {qs : List Qubit} (hqs : ∀ q ∈ qs, q < n)
    {gs : List Gate} {d : BlockState} (h : analyzeGates qs gs = some d) :
    ∀ g ∈ gs, g.InRange n m :=
  analyzeFold_inRange qs hqs gs _ d h

/-! ## The chunk invariant -/

/-- A chunk only ever holds wires the circuit has. -/
def Chunk.Ok (ch : Chunk) : Prop := ∀ q ∈ ch.qubits, q < ch.numQubits

theorem Chunk.extend_numQubits (ch : Chunk) (q : Qubit) :
    (ch.extend q).numQubits = ch.numQubits := by
  unfold Chunk.extend; split <;> rfl

theorem Chunk.extend_ok {ch : Chunk} {q : Qubit} (hok : ch.Ok) (hq : q < ch.numQubits) :
    (ch.extend q).Ok := by
  unfold Chunk.Ok Chunk.extend
  split
  · exact hok
  · intro r hr
    rcases List.mem_append.mp hr with hr | hr
    · exact hok r hr
    · simp only [List.mem_singleton] at hr; subst hr; exact hq

theorem Chunk.foldl_extend_numQubits : ∀ (ops : List Qubit) (ch : Chunk),
    (ops.foldl Chunk.extend ch).numQubits = ch.numQubits := by
  intro ops
  induction ops with
  | nil => intro ch; rfl
  | cons q ops ih => intro ch; rw [List.foldl_cons, ih, Chunk.extend_numQubits]

theorem Chunk.foldl_extend_original : ∀ (ops : List Qubit) (ch : Chunk),
    (ops.foldl Chunk.extend ch).original = ch.original := by
  intro ops
  induction ops with
  | nil => intro ch; rfl
  | cons q ops ih =>
      intro ch
      rw [List.foldl_cons, ih]
      unfold Chunk.extend; split <;> rfl

theorem Chunk.foldl_extend_ok : ∀ (ops : List Qubit) (ch : Chunk), ch.Ok →
    (∀ q ∈ ops, q < ch.numQubits) → (ops.foldl Chunk.extend ch).Ok := by
  intro ops
  induction ops with
  | nil => intro ch hok _; exact hok
  | cons q ops ih =>
      intro ch hok hlt
      rw [List.foldl_cons]
      refine ih _ (Chunk.extend_ok hok (hlt q (by simp))) ?_
      intro r hr
      rw [Chunk.extend_numQubits]
      exact hlt r (by simp [hr])

theorem Chunk.feedInto_spec {ch ch' : Chunk} {g : Gate} {out : List Gate}
    (h : ch.feedInto g = some (out, ch')) :
    out = [] ∧ ch' = { ch with state := ch'.state, original := ch.original ++ [g] } := by
  unfold Chunk.feedInto at h
  by_cases hcap : ch.capOk g = true
  · rw [if_neg (by simp [hcap])] at h
    rcases hfeed : feedGateFast ch.qubits ch.state g with _ | st
    · rw [hfeed] at h; simp at h
    · rw [hfeed] at h
      simp only [Option.some.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, rfl⟩ := h
      exact ⟨rfl, rfl⟩
  · rw [if_pos (by simpa using hcap)] at h; simp at h

theorem Chunk.feedTry_basic {ch ch' : Chunk} {g : Gate} {out : List Gate}
    (h : ch.feedTry g = some (out, ch')) :
    out = [] ∧ ch'.numQubits = ch.numQubits ∧ ch'.original = ch.original ++ [g] := by
  unfold Chunk.feedTry at h
  rcases hint : Chunk.interpretable g with _ | ops
  · rw [hint] at h; simp at h
  · rw [hint] at h
    dsimp only at h
    by_cases hrange : (ops.all fun q => decide (q < ch.numQubits)) = true
    · rw [if_neg (by simp [hrange])] at h
      by_cases hadm : ch.admits ops = true
      · rw [if_neg (by simp [hadm])] at h
        obtain ⟨rfl, hch'⟩ := Chunk.feedInto_spec h
        refine ⟨rfl, ?_, ?_⟩
        · rw [hch']; exact Chunk.foldl_extend_numQubits ops ch
        · rw [hch']
          simp only
          rw [Chunk.foldl_extend_original]
      · rw [if_pos (by simpa using hadm)] at h; simp at h
    · rw [if_pos (by simpa using hrange)] at h; simp at h

theorem Chunk.feedTry_ok {ch ch' : Chunk} {g : Gate} {out : List Gate}
    (hok : ch.Ok) (h : ch.feedTry g = some (out, ch')) : ch'.Ok := by
  unfold Chunk.feedTry at h
  rcases hint : Chunk.interpretable g with _ | ops
  · rw [hint] at h; simp at h
  · rw [hint] at h
    dsimp only at h
    by_cases hrange : (ops.all fun q => decide (q < ch.numQubits)) = true
    · rw [if_neg (by simp [hrange])] at h
      by_cases hadm : ch.admits ops = true
      · rw [if_neg (by simp [hadm])] at h
        have hops : ∀ q ∈ ops, q < ch.numQubits := by
          intro q hq
          have := List.all_eq_true.mp hrange q hq
          simpa using this
        have hnq := Chunk.foldl_extend_numQubits ops ch
        have hok₁ : (ops.foldl Chunk.extend ch).Ok := Chunk.foldl_extend_ok ops ch hok hops
        obtain ⟨rfl, hch'⟩ := Chunk.feedInto_spec h
        unfold Chunk.Ok
        rw [hch']
        intro r hr
        simpa [hnq] using hok₁ r hr
      · rw [if_pos (by simpa using hadm)] at h; simp at h
    · rw [if_pos (by simpa using hrange)] at h; simp at h

theorem Chunk.feedTry_spec {ch ch' : Chunk} {g : Gate} {out : List Gate}
    (hok : ch.Ok) (h : ch.feedTry g = some (out, ch')) :
    out = [] ∧ ch'.numQubits = ch.numQubits ∧ ch'.Ok ∧ ch'.original = ch.original ++ [g] :=
  let ⟨h1, h2, h3⟩ := Chunk.feedTry_basic h
  ⟨h1, h2, Chunk.feedTry_ok hok h, h3⟩

/-! ## Flushing a block -/

/-- **A flushed block denotes the block it replaces.** Either the original gates come back
unchanged, or the synthesized ones do — and then only because re-analysing them reproduced
the original's linear map and phase polynomial. -/
theorem flush_correct {n m : Nat} (ch : Chunk) (hn : ch.numQubits = n) (hok : ch.Ok) :
    Equivalent n m ch.flush ch.original := by
  have hrange : ∀ q ∈ ch.qubits, q < n := by
    intro q hq; rw [← hn]; exact hok q hq
  unfold Chunk.flush
  split
  · rename_i horig
    rw [List.isEmpty_iff.mp horig]
  · split
    · exact Equivalent.refl _ _ _
    · rename_i synth hsynth
      split
      · rename_i d₁ d₂ h₁ h₂
        split
        · rename_i hcond
          have hnorm : d₂.normalize = d₁.normalize := by
            have h' := (Bool.and_eq_true _ _).mp hcond
            exact (beq_iff_eq.mp h'.2).symm
          exact equivalent_of_analyze (m := m) hrange h₂ h₁ hnorm
        · exact Equivalent.refl _ _ _
      all_goals exact Equivalent.refl _ _ _

theorem flush_wf (ch : Chunk) (horig : ∀ g ∈ ch.original, g.Wf) : ∀ g ∈ ch.flush, g.Wf := by
  unfold Chunk.flush
  split
  · simp
  · split
    · exact horig
    · rename_i synth hsynth
      split
      · rename_i d₁ d₂ h₁ h₂
        split
        · exact analyzeGates_wf h₂
        · exact horig
      all_goals exact horig

theorem flush_inRange {n m : Nat} (ch : Chunk) (hn : ch.numQubits = n) (hok : ch.Ok)
    (horig : ∀ g ∈ ch.original, g.InRange n m) : ∀ g ∈ ch.flush, g.InRange n m := by
  have hqs : ∀ q ∈ ch.qubits, q < n := fun q hq => hn ▸ hok q hq
  unfold Chunk.flush
  split
  · simp
  · split
    · exact horig
    · rename_i synth hsynth
      split
      · rename_i d₁ d₂ h₁ h₂
        split
        · exact analyzeGates_inRange hqs h₂
        · exact horig
      all_goals exact horig

/-! ## Running the chunker -/

theorem Chunk.reset_ok (ch : Chunk) : ch.reset.Ok := by
  intro q hq; simp [Chunk.reset, Chunk.empty] at hq

theorem Chunk.reset_numQubits (ch : Chunk) : ch.reset.numQubits = ch.numQubits := rfl

theorem Chunk.reset_original (ch : Chunk) : ch.reset.original = [] := rfl

/-- **Feeding one gate preserves the circuit.** What has been emitted, followed by what the
block still holds, denotes what was consumed. -/
theorem feed_correct {n m : Nat} (ch : Chunk) (hn : ch.numQubits = n) (hok : ch.Ok) (g : Gate) :
    Equivalent n m ((ch.feed g).1 ++ (ch.feed g).2.original) (ch.original ++ [g]) := by
  unfold Chunk.feed
  split
  · rename_i res h
    obtain ⟨out, ch'⟩ := res
    obtain ⟨rfl, hnq, hok', horig⟩ := Chunk.feedTry_spec hok h
    simp only [horig, List.nil_append]
    exact Equivalent.refl _ _ _
  · rename_i h
    split
    · rename_i hempty
      have horig : ch.original = [] := by
        have h' := (Bool.and_eq_true _ _).mp hempty
        exact List.isEmpty_iff.mp h'.2
      simp only [horig, List.nil_append]
      simp only [List.append_nil]
      exact Equivalent.refl _ _ _
    · split
      · rename_i out ch' hfresh
        obtain ⟨rfl, hnq, hok', horig⟩ := Chunk.feedTry_spec (Chunk.reset_ok ch) hfresh
        rw [horig, Chunk.reset_original]
        simpa using Equivalent.append_right (n := n) (m := m) [g] (flush_correct ch hn hok)
      · rw [Chunk.reset_original]
        simpa using Equivalent.append_right (n := n) (m := m) [g] (flush_correct ch hn hok)

theorem feed_ok (ch : Chunk) (hok : ch.Ok) (g : Gate) :
    (ch.feed g).2.numQubits = ch.numQubits ∧ (ch.feed g).2.Ok := by
  unfold Chunk.feed
  split
  · rename_i res h
    obtain ⟨out, ch'⟩ := res
    obtain ⟨-, hnq, hok', -⟩ := Chunk.feedTry_spec hok h
    exact ⟨hnq, hok'⟩
  · split
    · exact ⟨rfl, hok⟩
    · split
      · rename_i out ch' hfresh
        obtain ⟨-, hnq, hok', -⟩ := Chunk.feedTry_spec (Chunk.reset_ok ch) hfresh
        exact ⟨by rw [hnq]; rfl, hok'⟩
      · exact ⟨rfl, Chunk.reset_ok ch⟩

theorem feed_wf (ch : Chunk) (g : Gate) (horig : ∀ g' ∈ ch.original, g'.Wf) (hg : g.Wf) :
    (∀ g' ∈ (ch.feed g).1, g'.Wf) ∧ (∀ g' ∈ (ch.feed g).2.original, g'.Wf) := by
  unfold Chunk.feed
  split
  · rename_i res h
    obtain ⟨out, ch'⟩ := res
    obtain ⟨rfl, -, horig'⟩ := Chunk.feedTry_basic h
    refine ⟨by simp, ?_⟩
    rw [horig']
    intro g' hg'
    rcases List.mem_append.mp hg' with h' | h'
    · exact horig g' h'
    · simp only [List.mem_singleton] at h'
      exact h' ▸ hg
  · split
    · exact ⟨by simpa using hg, horig⟩
    · split
      · rename_i out ch' hfresh
        obtain ⟨rfl, -, horig'⟩ := Chunk.feedTry_basic hfresh
        refine ⟨by simpa using flush_wf ch horig, ?_⟩
        rw [horig', Chunk.reset_original]
        intro g' hg'
        simp only [List.nil_append, List.mem_singleton] at hg'
        exact hg' ▸ hg
      · refine ⟨?_, by rw [Chunk.reset_original]; simp⟩
        intro g' hg'
        rcases List.mem_append.mp hg' with h' | h'
        · exact flush_wf ch horig g' h'
        · simp only [List.mem_singleton] at h'; exact h' ▸ hg

theorem feed_inRange {n m : Nat} (ch : Chunk) (hn : ch.numQubits = n) (hok : ch.Ok) (g : Gate)
    (horig : ∀ g' ∈ ch.original, g'.InRange n m) (hg : g.InRange n m) :
    (∀ g' ∈ (ch.feed g).1, g'.InRange n m) ∧
      (∀ g' ∈ (ch.feed g).2.original, g'.InRange n m) := by
  unfold Chunk.feed
  split
  · rename_i res h
    obtain ⟨out, ch'⟩ := res
    obtain ⟨rfl, -, horig'⟩ := Chunk.feedTry_basic h
    refine ⟨by simp, ?_⟩
    rw [horig']
    intro g' hg'
    rcases List.mem_append.mp hg' with h' | h'
    · exact horig g' h'
    · simp only [List.mem_singleton] at h'
      exact h' ▸ hg
  · split
    · exact ⟨by simpa using hg, horig⟩
    · split
      · rename_i out ch' hfresh
        obtain ⟨rfl, -, horig'⟩ := Chunk.feedTry_basic hfresh
        refine ⟨by simpa using flush_inRange ch hn hok horig, ?_⟩
        rw [horig', Chunk.reset_original]
        intro g' hg'
        simp only [List.nil_append, List.mem_singleton] at hg'
        exact hg' ▸ hg
      · refine ⟨?_, by rw [Chunk.reset_original]; simp⟩
        intro g' hg'
        rcases List.mem_append.mp hg' with h' | h'
        · exact flush_inRange ch hn hok horig g' h'
        · simp only [List.mem_singleton] at h'; exact h' ▸ hg

/-! ## The whole pass -/

theorem runGates_spec {n m : Nat} : ∀ (gs : List Gate) (ch : Chunk),
    ch.numQubits = n → ch.Ok → (∀ g ∈ ch.original, g.Wf) → (∀ g ∈ gs, g.Wf) →
    Equivalent n m ((runGates ch gs).1 ++ (runGates ch gs).2.original) (ch.original ++ gs) ∧
      (runGates ch gs).2.numQubits = ch.numQubits ∧ (runGates ch gs).2.Ok ∧
      (∀ g ∈ (runGates ch gs).1, g.Wf) ∧ (∀ g ∈ (runGates ch gs).2.original, g.Wf) := by
  intro gs
  induction gs with
  | nil =>
      intro ch hn hok hwf _
      refine ⟨?_, rfl, hok, by simp [runGates], by simpa [runGates] using hwf⟩
      simp only [runGates, List.nil_append, List.append_nil]
      exact Equivalent.refl _ _ _
  | cons g gs ih =>
      intro ch hn hok hwf hgs
      have hg : g.Wf := hgs g (by simp)
      have hgs' : ∀ g' ∈ gs, g'.Wf := fun g' hg' => hgs g' (by simp [hg'])
      obtain ⟨hnq, hok'⟩ := feed_ok ch hok g
      obtain ⟨hwf₁, hwf₂⟩ := feed_wf ch g hwf hg
      have hfeed := feed_correct (m := m) ch hn hok g
      obtain ⟨hih, hnq', hok'', hwf₃, hwf₄⟩ :=
        ih (ch.feed g).2 (by rw [hnq, hn]) hok' hwf₂ hgs'
      rcases hfe : ch.feed g with ⟨e, ch'⟩
      rcases hru : runGates ch' gs with ⟨rest, ch''⟩
      rw [hfe] at hnq hok' hwf₁ hwf₂ hfeed
      rw [hfe, hru] at hih hnq' hok'' hwf₃ hwf₄
      simp only at hnq hok' hwf₁ hwf₂ hfeed hih hnq' hok'' hwf₃ hwf₄
      have hrun : runGates ch (g :: gs) = (e ++ rest, ch'') := by
        rw [runGates, hfe]
        dsimp only
        rw [hru]
      rw [hrun]
      refine ⟨?_, by rw [hnq'] at *; simp [hnq], hok'', ?_, hwf₄⟩
      · have step₁ : Equivalent n m (e ++ (rest ++ ch''.original)) (e ++ (ch'.original ++ gs)) :=
          Equivalent.append_left e (by simpa using hih)
        have step₂ : Equivalent n m ((e ++ ch'.original) ++ gs) ((ch.original ++ [g]) ++ gs) :=
          Equivalent.append_right gs hfeed
        refine Equivalent.trans (by simpa [List.append_assoc] using step₁) ?_
        simpa [List.append_assoc] using step₂
      · intro g' hg'
        rcases List.mem_append.mp hg' with h' | h'
        · exact hwf₁ g' h'
        · exact hwf₃ g' h'

theorem runGates_inRange {n m : Nat} : ∀ (gs : List Gate) (ch : Chunk),
    ch.numQubits = n → ch.Ok → (∀ g ∈ ch.original, g.InRange n m) →
    (∀ g ∈ gs, g.InRange n m) →
    (∀ g ∈ (runGates ch gs).1, g.InRange n m) ∧
      (∀ g ∈ (runGates ch gs).2.original, g.InRange n m) := by
  intro gs
  induction gs with
  | nil => intro ch _ _ hin _; exact ⟨by simp [runGates], by simpa [runGates] using hin⟩
  | cons g gs ih =>
      intro ch hn hok hin hgs
      have hg : g.InRange n m := hgs g (by simp)
      have hgs' : ∀ g' ∈ gs, g'.InRange n m := fun g' hg' => hgs g' (by simp [hg'])
      obtain ⟨hnq, hok'⟩ := feed_ok ch hok g
      obtain ⟨hin₁, hin₂⟩ := feed_inRange (m := m) ch hn hok g hin hg
      obtain ⟨hin₃, hin₄⟩ := ih (ch.feed g).2 (by rw [hnq, hn]) hok' hin₂ hgs'
      rcases hfe : ch.feed g with ⟨e, ch'⟩
      rcases hru : runGates ch' gs with ⟨rest, ch''⟩
      rw [hfe] at hin₁ hin₂
      rw [hfe, hru] at hin₃ hin₄
      simp only at hin₁ hin₂ hin₃ hin₄
      have hrun : runGates ch (g :: gs) = (e ++ rest, ch'') := by
        rw [runGates, hfe]
        dsimp only
        rw [hru]
      rw [hrun]
      refine ⟨?_, hin₄⟩
      intro g' hg'
      rcases List.mem_append.mp hg' with h' | h'
      · exact hin₁ g' h'
      · exact hin₃ g' h'

/-- The chunker's bookkeeping alone: the wire count and `Chunk.Ok` survive a run, with no
assumption about the gates. `runGates_spec` proves this too, but only alongside the
equivalence, which needs hypotheses the range argument has no reason to carry. -/
theorem runGates_ok : ∀ (gs : List Gate) (ch : Chunk), ch.Ok →
    (runGates ch gs).2.numQubits = ch.numQubits ∧ (runGates ch gs).2.Ok := by
  intro gs
  induction gs with
  | nil => intro ch hok; exact ⟨rfl, hok⟩
  | cons g gs ih =>
      intro ch hok
      obtain ⟨hnq, hok'⟩ := feed_ok ch hok g
      obtain ⟨hnq', hok''⟩ := ih (ch.feed g).2 hok'
      rcases hfe : ch.feed g with ⟨e, ch'⟩
      rcases hru : runGates ch' gs with ⟨rest, ch''⟩
      rw [hfe] at hnq hok'
      rw [hfe, hru] at hnq' hok''
      simp only at hnq hok' hnq' hok''
      have hrun : runGates ch (g :: gs) = (e ++ rest, ch'') := by
        rw [runGates, hfe]
        dsimp only
        rw [hru]
      rw [hrun]
      exact ⟨hnq'.trans hnq, hok''⟩

/-- **`CnotMin` keeps every operand in range.** -/
theorem cnotMinGates_inRange {n m : Nat} (maxQ maxT : Nat) (gs : List Gate)
    (hin : ∀ g ∈ gs, g.InRange n m) :
    ∀ g ∈ cnotMinGates n maxQ maxT gs, g.InRange n m := by
  have hok : (Chunk.empty n maxQ maxT).Ok := by intro q hq; simp [Chunk.empty] at hq
  obtain ⟨hnq, hok'⟩ := runGates_ok gs (Chunk.empty n maxQ maxT) hok
  obtain ⟨hin₁, hin₂⟩ :=
    runGates_inRange (m := m) gs (Chunk.empty n maxQ maxT) rfl hok (by simp [Chunk.empty]) hin
  intro g hg
  simp only [cnotMinGates] at hg
  rcases List.mem_append.mp hg with h' | h'
  · exact hin₁ g h'
  · exact flush_inRange _ (by rw [hnq]) hok' hin₂ g h'

/-- **`CnotMin` preserves the circuit.** -/
theorem cnotMinGates_correct {n m : Nat} (maxQ maxT : Nat) (gs : List Gate)
    (hwf : ∀ g ∈ gs, g.Wf) : Equivalent n m (cnotMinGates n maxQ maxT gs) gs := by
  have hok : (Chunk.empty n maxQ maxT).Ok := by intro q hq; simp [Chunk.empty] at hq
  obtain ⟨hcorr, hnq, hok', hwf₁, hwf₂⟩ :=
    runGates_spec (n := n) (m := m) gs (Chunk.empty n maxQ maxT) rfl hok (by simp [Chunk.empty]) hwf
  set res := runGates (Chunk.empty n maxQ maxT) gs with hres
  have hflush : Equivalent n m res.2.flush res.2.original :=
    flush_correct res.2 (by rw [hnq]; rfl) hok'
  have : Equivalent n m (res.1 ++ res.2.flush) (res.1 ++ res.2.original) :=
    Equivalent.append_left res.1 hflush
  refine Equivalent.trans this ?_
  simpa [Chunk.empty] using hcorr

theorem cnotMinGates_wf {n : Nat} (maxQ maxT : Nat) (gs : List Gate) (hwf : ∀ g ∈ gs, g.Wf) :
    ∀ g ∈ cnotMinGates n maxQ maxT gs, g.Wf := by
  have hok : (Chunk.empty n maxQ maxT).Ok := by intro q hq; simp [Chunk.empty] at hq
  obtain ⟨-, -, -, hwf₁, hwf₂⟩ :=
    runGates_spec (n := n) (m := 0) gs (Chunk.empty n maxQ maxT) rfl hok (by simp [Chunk.empty]) hwf
  intro g hg
  simp only [cnotMinGates] at hg
  rcases List.mem_append.mp hg with h' | h'
  · exact hwf₁ g h'
  · exact flush_wf _ hwf₂ g h'

/-- Run CNOT minimization on a raw circuit. -/
def cnotMinCircuit (c : Circuit) : Circuit :=
  c.withGates (cnotMinGates c.numQubits maxChunkQubits maxChunkTerms c.gates)

/-- **The `CnotMin` pass**, with its correctness proof as a field. -/
def CnotMin : Pass where
  name := "CNOT minimization"
  run := fun c => ⟨cnotMinCircuit c.raw, c.numQubits_eq, c.numCbits_eq,
    cnotMinGates_wf _ _ c.raw.gates c.wf⟩
  correct := by
    intro n m c
    rcases c with ⟨c, rfl, rfl, hc⟩
    exact cnotMinGates_correct _ _ c.gates hc

end
end TzapLean
