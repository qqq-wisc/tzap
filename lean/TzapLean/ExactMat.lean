import TzapLean.Cyclotomic
import TzapLean.Merge

/-!
# Exact Clifford+T Matrices

A window's meaning, computed without floating point: a `2ⁿ × 2ⁿ` matrix of cyclotomic
integers over one shared `√2^den`. Gates act by the row operations of
`src/super_opt/matrix.rs` — Hadamard adds and subtracts row pairs and bumps `den`, phase
gates scale rows by `ω^p`, controlled-`X` permutes rows — and each is proved to be exactly
left multiplication by that gate's matrix.

The index type is `Basis n`, not a packed row number. Rust has to pick a bit order (qubit
zero is the most significant bit) and translate; here `out.get q` *is* the bit, so the gate
lemmas are direct.

`matrixOf` builds a gate list's matrix, returning `none` on the gates SuperOpt treats as
window barriers (`rz`, `measure`, `reset`). `matrixOf_sound` is what the pass rests on:

```lean
matrixOf n gs = some M  →  M.interp = unitary n gs
```
-/

namespace TzapLean

open Cyc

/-- `√2`, as a complex number. -/
noncomputable def sq2 : ℂ := ((Real.sqrt 2 : ℝ) : ℂ)

theorem sq2_ne_zero : sq2 ≠ 0 := ofReal_sqrt_two_ne_zero

@[simp] theorem sq2_mul_self : sq2 * sq2 = 2 := ofReal_sqrt_two_mul_self

/-- A `2ⁿ × 2ⁿ` exact matrix: cyclotomic-integer entries over a shared `√2^den`. -/
structure ExactMat (n : Nat) where
  /-- The shared denominator exponent: every entry is divided by `√2^den`. -/
  den : Nat
  /-- The numerators. -/
  get : Basis n → Basis n → Cyc

namespace ExactMat

/-- The complex matrix an exact matrix denotes. -/
noncomputable def interp {n : Nat} (M : ExactMat n) : Density n :=
  fun out inp => (M.get out inp).interp / sq2 ^ M.den

/-- The identity, denominator zero. -/
def id (n : Nat) : ExactMat n where
  den := 0
  get := fun out inp => if out = inp then 1 else 0

@[simp] theorem interp_id (n : Nat) : (id n).interp = 1 := by
  funext out inp
  by_cases h : out = inp <;> simp [interp, id, h, Matrix.one_apply]

/-! ## Gate actions -/

/-- Scale the rows where `q` is set by `ω^p` — the phase gates. -/
def applyPhase {n : Nat} (q : Qubit) (p : Nat) (M : ExactMat n) : ExactMat n where
  den := M.den
  get := fun out inp => if out.get q then (M.get out inp).timesOmega p else M.get out inp

/-- Scale the rows where every wire in `qs` is set by `ω^p` — `z`, `cz`, `ccz`. -/
def applyPhaseMask {n : Nat} (qs : List Qubit) (p : Nat) (M : ExactMat n) : ExactMat n where
  den := M.den
  get := fun out inp =>
    if qs.all (fun q => out.get q) then (M.get out inp).timesOmega p else M.get out inp

/-- Permute rows by flipping `tgt` where every control is set — `x`, `cnot`, `ccx`. -/
def applyX {n : Nat} (ctrls : List Qubit) (tgt : Qubit) (M : ExactMat n) : ExactMat n where
  den := M.den
  get := fun out inp =>
    M.get (out.set tgt (out.get tgt != ctrls.all (fun q => out.get q))) inp

/-- Add and subtract row pairs across `q`, and bump the denominator — the Hadamard. -/
def applyH {n : Nat} (q : Qubit) (M : ExactMat n) : ExactMat n :=
  if q < n then
    { den := M.den + 1
      get := fun out inp =>
        if out.get q then M.get (out.set q false) inp - M.get out inp
        else M.get out inp + M.get (out.set q true) inp }
  else M

/-! ## Each action is a left multiplication -/

/-- A permutation of basis states acts on a matrix by permuting rows. -/
theorem permMatrix_mul_apply {n : Nat} {σ : Basis n → Basis n} (hinv : ∀ b, σ (σ b) = b)
    (M : Density n) (out inp : Basis n) : (permMatrix σ * M) out inp = M (σ out) inp := by
  rw [Matrix.mul_apply, Finset.sum_eq_single (σ out)]
  · rw [permMatrix, if_pos (hinv out).symm, one_mul]
  · intro k _ hk
    rw [permMatrix, if_neg, zero_mul]
    intro hc
    exact hk (by rw [hc, hinv])
  · simp

theorem interp_applyPhase {n : Nat} (q : Qubit) (p : Nat) (M : ExactMat n) :
    (applyPhase q p M).interp = phaseMatrix (fun b : Basis n => if b.get q then ω ^ p else 1) *
      M.interp := by
  funext out inp
  rw [phaseMatrix_mul_apply]
  by_cases h : out.get q = true
  · simp only [interp, applyPhase, if_pos h, Cyc.interp_timesOmega]
    rw [mul_div_assoc]
  · simp only [interp, applyPhase, if_neg h, one_mul]

theorem interp_applyPhaseMask {n : Nat} (qs : List Qubit) (p : Nat) (M : ExactMat n) :
    (applyPhaseMask qs p M).interp =
      phaseMatrix (fun b : Basis n => if qs.all (fun q => b.get q) then ω ^ p else 1) *
        M.interp := by
  funext out inp
  rw [phaseMatrix_mul_apply]
  by_cases h : qs.all (fun q => out.get q) = true
  · simp only [interp, applyPhaseMask, if_pos h, Cyc.interp_timesOmega]
    rw [mul_div_assoc]
  · simp only [interp, applyPhaseMask, if_neg h, one_mul]

theorem interp_applyX {n : Nat} (ctrls : List Qubit) (tgt : Qubit) (M : ExactMat n)
    {σ : Basis n → Basis n}
    (hσ : σ = fun b => b.set tgt (b.get tgt != ctrls.all (fun q => b.get q)))
    (hinv : ∀ b, σ (σ b) = b) :
    (applyX ctrls tgt M).interp = permMatrix σ * M.interp := by
  funext out inp
  rw [permMatrix_mul_apply hinv, hσ]
  rfl

/-! ### The Hadamard -/

theorem h2_apply (o i : Bool) : h2 o i = (if o && i then -1 else 1) / sq2 := rfl

theorem embed1_apply_set {n : Nat} (M : Bool → Bool → ℂ) {q : Qubit} (hq : q < n)
    (out : Basis n) (v : Bool) : embed1 n M q out (out.set q v) = M (out.get q) v := by
  simp only [embed1_apply_of_lt _ hq]
  rw [if_pos (fun r hr => by simp [Basis.set, hr])]
  congr 1
  · simp [Basis.get, hq]
  · simp [Basis.set, hq]

theorem embed1_apply_eq_zero {n : Nat} (M : Bool → Bool → ℂ) {q : Qubit} (hq : q < n)
    (out k : Basis n) (h : k ≠ out.set q false) (h' : k ≠ out.set q true) :
    embed1 n M q out k = 0 := by
  simp only [embed1_apply_of_lt _ hq]
  refine if_neg (fun hall => ?_)
  have hk : k = out.set q (k ⟨q, hq⟩) := by
    funext r
    by_cases hr : (r : Nat) = q
    · have hrq : r = ⟨q, hq⟩ := Fin.ext hr
      rw [hrq]; simp [Basis.set, hq]
    · rw [Basis.set]
      simp only [hr, if_neg]
      exact (hall r hr).symm
  cases hkq : k ⟨q, hq⟩
  · exact h (by rw [hk, hkq])
  · exact h' (by rw [hk, hkq])

theorem embed1_row_sum {n : Nat} (M : Bool → Bool → ℂ) {q : Qubit} (hq : q < n)
    (out : Basis n) (f : Basis n → ℂ) :
    (∑ k : Basis n, embed1 n M q out k * f k)
      = M (out.get q) false * f (out.set q false) + M (out.get q) true * f (out.set q true) := by
  have hne : out.set q false ≠ out.set q true := by
    intro hc
    have := congrFun hc ⟨q, hq⟩
    simp [Basis.set, hq] at this
  rw [Finset.sum_eq_add_of_mem (out.set q false) (out.set q true) (Finset.mem_univ _)
    (Finset.mem_univ _) hne ?_]
  · rw [embed1_apply_set M hq out false, embed1_apply_set M hq out true]
  · intro k _ hk
    rw [embed1_apply_eq_zero M hq out k hk.1 hk.2, zero_mul]

theorem interp_applyH {n : Nat} (q : Qubit) (M : ExactMat n) :
    (applyH q M).interp = embed1 n h2 q * M.interp := by
  by_cases hq : q < n
  · have hne : sq2 ≠ 0 := sq2_ne_zero
    funext out inp
    have hsum : (embed1 n h2 q * M.interp) out inp =
        h2 (out.get q) false * M.interp (out.set q false) inp +
          h2 (out.get q) true * M.interp (out.set q true) inp := by
      rw [Matrix.mul_apply]
      exact embed1_row_sum h2 hq out (fun k => M.interp k inp)
    rw [hsum]
    by_cases hb : out.get q = true
    · have hself : out.set q true = out := by
        rw [← hb]; exact Basis.set_get_self out q
      simp only [applyH, if_pos hq, interp, hb, if_true, hself, h2_apply, Bool.true_and,
        Bool.and_true, if_false, if_true, Cyc.interp_sub, pow_succ, reduceIte]
      field_simp
      norm_num [Cyc.interp_add, Cyc.interp_sub]
      try ring
    · have hbf : out.get q = false := by simpa using hb
      have hself : out.set q false = out := by
        rw [← hbf]; exact Basis.set_get_self out q
      simp only [applyH, if_pos hq, interp, hbf, if_false, hself, h2_apply, Bool.false_and,
        if_false, Cyc.interp_add, pow_succ, reduceIte]
      field_simp
      norm_num [Cyc.interp_add, Cyc.interp_sub]
      try ring
  · rw [applyH, if_neg hq, embed1_eq_one _ hq, Matrix.one_mul]

/-! ## Building a gate list's matrix -/

theorem omega_pow (p : Nat) : ω ^ p = ep ((p : ℚ) / 4) := by
  induction p with
  | zero => simp
  | succ p ih =>
      rw [pow_succ, ih, ω, ← ep_add]
      congr 1
      push_cast
      ring

theorem ep_add_two (θ : ℚ) : ep (θ + 2) = ep θ := by
  rw [ep_add, show (2 : ℚ) = 1 + 1 by norm_num, ep_add, ep_one]
  ring

/-- Controls that avoid the target make the row permutation an involution. -/
theorem all_get_set_ne {n : Nat} (ctrls : List Qubit) {tgt : Qubit} (v : Bool) (b : Basis n)
    (h : ∀ q ∈ ctrls, q ≠ tgt) :
    ctrls.all (fun q => (b.set tgt v).get q) = ctrls.all (fun q => b.get q) := by
  induction ctrls with
  | nil => rfl
  | cons c cs ih =>
      simp only [List.all_cons]
      rw [Basis.get_set_ne _ _ _ _ (h c (by simp)), ih (fun q hq => h q (by simp [hq]))]

theorem flipPerm_involutive {n : Nat} {ctrls : List Qubit} {tgt : Qubit}
    (h : ∀ q ∈ ctrls, q ≠ tgt) (b : Basis n) :
    (fun b : Basis n => b.set tgt (b.get tgt != ctrls.all (fun q => b.get q)))
      ((fun b : Basis n => b.set tgt (b.get tgt != ctrls.all (fun q => b.get q))) b) = b := by
  simp only
  rw [all_get_set_ne ctrls _ b h]
  by_cases htn : tgt < n
  · rw [Basis.get_set_same _ _ _ htn, Basis.set_set]
    have : ((b.get tgt != ctrls.all (fun q => b.get q)) != ctrls.all (fun q => b.get q))
        = b.get tgt := by
      cases b.get tgt <;> cases ctrls.all (fun q => b.get q) <;> rfl
    rw [this, Basis.set_get_self]
  · rw [Basis.set_eq_self_of_ge _ _ _ htn, Basis.set_eq_self_of_ge _ _ _ htn]

/-- One gate's action, or `none` for the gates SuperOpt treats as window barriers. -/
def applyGate {n : Nat} (g : Gate) (M : ExactMat n) : Option (ExactMat n) :=
  match g with
  | .x q => some (applyX [] q M)
  | .h q => some (applyH q M)
  | .t q => some (applyPhase q 1 M)
  | .tdg q => some (applyPhase q 7 M)
  | .s q => some (applyPhase q 2 M)
  | .sdg q => some (applyPhase q 6 M)
  | .z q => some (applyPhase q 4 M)
  | .cnot c t => some (applyX [c] t M)
  | .ccx c₁ c₂ t => some (applyX [c₁, c₂] t M)
  | .cz c t => some (applyPhaseMask [c, t] 4 M)
  | .ccz c₁ c₂ t => some (applyPhaseMask [c₁, c₂, t] 4 M)
  | _ => none

/-- A phase gate's matrix, as a phase matrix on its wire — in range or not. -/
theorem gateUnitary_phase {n : Nat} (q : Qubit) (a : ℚ) (p : Nat) (hp : ω ^ p = ep a) :
    embed1 n (diag2 1 (ep a)) q =
      phaseMatrix (fun b : Basis n => if b.get q then ω ^ p else 1) := by
  by_cases hq : q < n
  · rw [embed1_diag2_eq_phaseMatrix _ _ _ hq, hp]
  · rw [embed1_eq_one _ hq, ← phaseMatrix_one]
    congr 1
    funext b
    rw [basis_get_of_ge b hq, if_neg (by simp)]

theorem interp_applyGate {n : Nat} {g : Gate} {M M' : ExactMat n} (hwf : g.Wf)
    (h : applyGate g M = some M') : M'.interp = gateUnitary n g * M.interp := by
  cases g with
  | x q =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyX [] q M rfl (flipPerm_involutive (by simp)), gateUnitary,
        embed1_x2_eq_perm]
      congr 1
      funext b
      simp
  | h q =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      exact interp_applyH q M
  | t q =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyPhase, gateUnitary, gateUnitary_phase q (1/4) 1 (by rw [omega_pow]; norm_num)]
  | tdg q =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyPhase, gateUnitary,
        gateUnitary_phase q (-1/4) 7 (by rw [omega_pow, show ((7 : Nat) : ℚ)/4 = -1/4 + 2 by
          norm_num, ep_add_two])]
  | s q =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyPhase, gateUnitary, gateUnitary_phase q (1/2) 2 (by rw [omega_pow]; norm_num)]
  | sdg q =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyPhase, gateUnitary,
        gateUnitary_phase q (-1/2) 6 (by rw [omega_pow, show ((6 : Nat) : ℚ)/4 = -1/2 + 2 by
          norm_num, ep_add_two])]
  | z q =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyPhase, gateUnitary, gateUnitary_phase q 1 4 (by rw [omega_pow]; norm_num)]
  | cnot c t =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      have hct : c ≠ t := by simpa [Gate.Wf] using hwf
      rw [interp_applyX [c] t M rfl (flipPerm_involutive (by simpa using hct)), gateUnitary]
      congr 1
      funext b
      simp
  | ccx c₁ c₂ t =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      have hwf' : c₁ ≠ c₂ ∧ c₁ ≠ t ∧ c₂ ≠ t := by simpa [Gate.Wf] using hwf
      rw [interp_applyX [c₁, c₂] t M rfl (flipPerm_involutive (by
        intro q hq
        rcases List.mem_cons.1 hq with rfl | hq
        · exact hwf'.2.1
        · rw [List.mem_singleton.1 hq]; exact hwf'.2.2)), gateUnitary]
      congr 1
      funext b
      simp
  | cz c t =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyPhaseMask, gateUnitary]
      congr 1
      funext b
      simp only [List.all_cons, List.all_nil, Bool.and_true, ω_pow_four]
  | ccz c₁ c₂ t =>
      rw [applyGate, Option.some.injEq] at h
      subst h
      rw [interp_applyPhaseMask, gateUnitary]
      congr 1
      funext b
      simp only [List.all_cons, List.all_nil, Bool.and_true, ω_pow_four, Bool.and_assoc]
  | rz θ q => simp [applyGate] at h
  | measure q c => simp [applyGate] at h
  | reset q => simp [applyGate] at h

/-- The matrix of a gate list, applied on the left of `M`. -/
def matrixOfFrom {n : Nat} (M : ExactMat n) : List Gate → Option (ExactMat n)
  | [] => some M
  | g :: gs => (applyGate g M).bind (fun M' => matrixOfFrom M' gs)

/-- The matrix of a gate list. `none` when it contains a window barrier. -/
def matrixOf (n : Nat) (gs : List Gate) : Option (ExactMat n) := matrixOfFrom (id n) gs

theorem interp_matrixOfFrom {n : Nat} : ∀ (gs : List Gate) (M M' : ExactMat n),
    (∀ g ∈ gs, g.Wf) → matrixOfFrom M gs = some M' → M'.interp = unitary n gs * M.interp := by
  intro gs
  induction gs with
  | nil =>
      intro M M' _ h
      rw [matrixOfFrom, Option.some.injEq] at h
      subst h
      rw [unitary_nil, Matrix.one_mul]
  | cons g gs ih =>
      intro M M' hwf h
      rw [matrixOfFrom] at h
      cases hg : applyGate g M with
      | none => rw [hg] at h; exact absurd h (by simp)
      | some M₁ =>
          rw [hg, Option.bind_some] at h
          rw [ih M₁ M' (fun x hx => hwf x (by simp [hx])) h,
            interp_applyGate (hwf g (by simp)) hg, unitary_cons, Matrix.mul_assoc]

/-- **The exact matrix is the circuit's matrix.** -/
theorem matrixOf_sound {n : Nat} {gs : List Gate} {M : ExactMat n} (hwf : ∀ g ∈ gs, g.Wf)
    (h : matrixOf n gs = some M) : M.interp = unitary n gs := by
  rw [interp_matrixOfFrom gs (id n) M hwf h, interp_id, Matrix.mul_one]

/-! ## Normalizing the denominator -/

/-- Every `Bool` vector of a given length. -/
def boolVecs : Nat → List (List Bool)
  | 0 => [[]]
  | k + 1 => (boolVecs k).flatMap fun l => [false :: l, true :: l]

theorem mem_boolVecs : ∀ l : List Bool, l ∈ boolVecs l.length
  | [] => by simp [boolVecs]
  | b :: l => by
      refine List.mem_flatMap.2 ⟨l, mem_boolVecs l, ?_⟩
      cases b <;> simp

/-- Every basis state, in a fixed order. Deciding a property of a matrix by folding this
list is far cheaper at runtime than deciding a quantifier over the function type, and
`mem_basisList` makes the two interchangeable in proofs. -/
def basisList (n : Nat) : List (Basis n) :=
  (boolVecs n).map fun l => fun j => l.getD (j : Nat) false

theorem mem_basisList {n : Nat} (b : Basis n) : b ∈ basisList n := by
  have hlen : ((List.finRange n).map b).length = n := by simp
  refine List.mem_map.2 ⟨(List.finRange n).map b, ?_, ?_⟩
  · have := mem_boolVecs ((List.finRange n).map b)
    rwa [hlen] at this
  · funext j
    rw [List.getD_eq_getElem?_getD, List.getElem?_map,
      List.getElem?_eq_getElem (by simp)]
    simp

/-- The matrix's entries in row-major order. -/
def entries {n : Nat} (M : ExactMat n) : List Cyc :=
  (basisList n).flatMap fun r => (basisList n).map fun c => M.get r c

theorem mem_entries {n : Nat} (M : ExactMat n) (out inp : Basis n) :
    M.get out inp ∈ M.entries :=
  List.mem_flatMap.2 ⟨out, mem_basisList out,
    List.mem_map.2 ⟨inp, mem_basisList inp, rfl⟩⟩

/-- Whether every entry is divisible by `√2`. -/
def allDivisible {n : Nat} (M : ExactMat n) : Bool :=
  M.entries.all fun x => x.divisibleBySqrt2

theorem allDivisible_spec {n : Nat} {M : ExactMat n} (h : M.allDivisible = true)
    (out inp : Basis n) : (M.get out inp).divisibleBySqrt2 = true :=
  (List.all_eq_true.1 h) _ (mem_entries M out inp)

/-- Strip common factors of `√2`, as Rust does after every Hadamard. -/
def normalizeAux {n : Nat} : Nat → ExactMat n → ExactMat n
  | 0, M => M
  | fuel + 1, M =>
      match hd : M.den with
      | 0 => M
      | d + 1 =>
          if allDivisible M then
            normalizeAux fuel ⟨d, fun out inp => (M.get out inp).divSqrt2⟩
          else M

/-- The matrix with its denominator reduced as far as it exactly can be. -/
def normalize {n : Nat} (M : ExactMat n) : ExactMat n := normalizeAux M.den M

theorem interp_normalizeAux {n : Nat} : ∀ (fuel : Nat) (M : ExactMat n),
    (normalizeAux fuel M).interp = M.interp := by
  intro fuel
  induction fuel with
  | zero => intro M; rfl
  | succ fuel ih =>
      intro M
      rw [normalizeAux]
      split
      · rfl
      · rename_i d hd
        split
        · rename_i hdiv
          rw [ih]
          funext out inp
          have hdiv' : (M.get out inp).divisibleBySqrt2 = true := allDivisible_spec hdiv out inp
          have hspec : sq2 * (M.get out inp).divSqrt2.interp = (M.get out inp).interp :=
            Cyc.divSqrt2_spec hdiv'
          simp only [interp, hd]
          rw [← hspec, pow_succ]
          rw [mul_comm sq2 ((M.get out inp).divSqrt2.interp), mul_div_mul_right _ _ sq2_ne_zero]
        · rfl

/-- **Normalizing does not change the matrix it denotes.** -/
@[simp] theorem interp_normalize {n : Nat} (M : ExactMat n) : (normalize M).interp = M.interp :=
  interp_normalizeAux M.den M

/-! ## Comparing two matrices up to global phase -/

/-- The power of `ω` carrying `M` to `N`, if one does — Rust canonicalizes both phases and
compares; searching the eight powers directly is the same test. -/
def phaseMatch {n : Nat} (M N : ExactMat n) : Option Nat :=
  if M.den = N.den then
    (List.range 8).find? fun p =>
      decide (∀ out inp : Basis n, (M.get out inp).timesOmega p = N.get out inp)
  else none

/-- **The phase test is sound**: a match means the matrices differ by that global phase. -/
theorem phaseMatch_sound {n : Nat} {M N : ExactMat n} {p : Nat} (h : phaseMatch M N = some p) :
    N.interp = ω ^ p • M.interp := by
  rw [phaseMatch] at h
  split at h
  · rename_i hden
    have hp := List.find?_some h
    have hall : ∀ out inp : Basis n, (M.get out inp).timesOmega p = N.get out inp :=
      of_decide_eq_true hp
    funext out inp
    simp only [interp, ← hall out inp, Cyc.interp_timesOmega, hden, Matrix.smul_apply,
      smul_eq_mul, mul_div_assoc]
  · exact absurd h (by simp)

/-- `ω^p` is a unit, so a phase match is invisible to the semantics. -/
theorem omega_pow_unit (p : Nat) : (ω ^ p) * star (ω ^ p) = 1 := by
  rw [omega_pow]
  exact ep_mul_star _

end ExactMat

end TzapLean
