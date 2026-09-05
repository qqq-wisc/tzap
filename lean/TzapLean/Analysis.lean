import TzapLean.Hash
import TzapLean.Cancel

/-!
# The Parity Analysis and Its Soundness

`PhaseFoldRand` merges two rotations when the wires they act on carry the same *parity* —
the same XOR of computational-basis values. This file defines the analysis that computes
those parities and proves the fact the whole pass rests on:

> **`analyze_sound`.** If a gate list has a nonzero amplitude from `b` to `b'`, then the
> output values `b'` are exactly what the analysis's parities predict, at *some* valuation
> of the variables extending the one describing `b`.

The variables are the initial wire values plus one fresh variable per `h` and per `ccx`
target — the two places where a wire's new value is not an affine function of the old ones.
Allocating a fresh variable there is what lets the analysis stay affine while remaining
sound: the valuation simply records whichever branch the amplitude took.

The transfer functions mirror `src/phase_fold_rand.rs`: `x` complements, `cnot` XORs, `h`
and `ccx` allocate, `measure` and every diagonal gate leave parities alone. `reset` allocates
too, where Rust records the constant `0`: the weaker fact keeps every analysis state
*generic* (any basis state is explainable), which is what the induction in `PhaseFoldProof`
needs, and phase folding never merges across a `reset` here anyway.
-/

namespace TzapLean

noncomputable section

/-- The analysis state: one parity per wire, and the next unused variable. -/
structure AState where
  /-- Parity of each wire, indexed by wire number. -/
  par : List Form
  /-- Next free variable index. -/
  fresh : Nat

namespace AState

/-- The parity of wire `q` (the constant `0` for a wire the state does not cover). -/
def parOf (st : AState) (q : Qubit) : Form := st.par.getD q (Form.const false)

/-- The state a circuit starts in: wire `i` holds variable `i`. -/
def initial (n : Nat) : AState where
  par := (List.range n).map Form.var
  fresh := n

/-- The gates that allocate a fresh variable: `h`, `ccx` and `reset` are the ones whose
output parity the affine analysis cannot express in terms of the input's. Counting them
rather than every gate is what keeps the drawn sample — one tag per variable — to the size
the analysis actually needs. -/
def _root_.TzapLean.Gate.allocates : Gate → Bool
  | .h _ | .ccx .. | .reset _ => true
  | _ => false

/-- One step of the analysis — the Rust transfer functions. -/
def step (st : AState) (g : Gate) : AState :=
  match g with
  | .x q => { st with par := st.par.set q (st.parOf q).flip }
  | .cnot c t => { st with par := st.par.set t ((st.parOf t) + (st.parOf c)) }
  | .h q => { par := st.par.set q (Form.var st.fresh), fresh := st.fresh + 1 }
  | .ccx _ _ t => { par := st.par.set t (Form.var st.fresh), fresh := st.fresh + 1 }
  | .reset q => { par := st.par.set q (Form.var st.fresh), fresh := st.fresh + 1 }
  | _ => st

/-- The analysis of a gate list, from a starting state. -/
def steps (st : AState) : List Gate → AState
  | [] => st
  | g :: gs => steps (st.step g) gs

@[simp] theorem steps_nil (st : AState) : st.steps [] = st := rfl
@[simp] theorem steps_cons (st : AState) (g : Gate) (gs : List Gate) :
    st.steps (g :: gs) = (st.step g).steps gs := rfl

theorem steps_append (st : AState) (gs hs : List Gate) :
    st.steps (gs ++ hs) = (st.steps gs).steps hs := by
  induction gs generalizing st with
  | nil => rfl
  | cons g gs ih => simp [ih]

/-! ## Invariants -/

/-- Every parity mentions only variables the state has allocated. -/
def Bounded (st : AState) : Prop := ∀ q : Qubit, Form.Bounded st.fresh (st.parOf q)

theorem length_step (st : AState) (g : Gate) : (st.step g).par.length = st.par.length := by
  cases g <;> simp [step]

theorem length_steps (st : AState) (gs : List Gate) :
    (st.steps gs).par.length = st.par.length := by
  induction gs generalizing st with
  | nil => rfl
  | cons g gs ih => rw [steps_cons, ih, length_step]

theorem fresh_le_step (st : AState) (g : Gate) : st.fresh ≤ (st.step g).fresh := by
  cases g <;> simp [step]

theorem fresh_le_steps (st : AState) (gs : List Gate) : st.fresh ≤ (st.steps gs).fresh := by
  induction gs generalizing st with
  | nil => exact le_rfl
  | cons g gs ih => exact le_trans (fresh_le_step st g) (ih (st.step g))

theorem parOf_mk (l : List Form) (fr : Nat) (q : Qubit) :
    (⟨l, fr⟩ : AState).parOf q = l.getD q (Form.const false) := rfl

theorem getD_set_self (l : List Form) (q : Qubit) (f : Form) (h : q < l.length) :
    (l.set q f).getD q (Form.const false) = f := by
  rw [List.getD_eq_getElem?_getD, List.getElem?_set_self]
  · rfl
  · exact h

theorem getD_set_ne (l : List Form) {q r : Qubit} (f : Form) (h : r ≠ q) :
    (l.set q f).getD r (Form.const false) = l.getD r (Form.const false) := by
  simp [List.getD_eq_getElem?_getD, Ne.symm h]

theorem getD_set_out (l : List Form) (q : Qubit) (f : Form) (h : ¬ q < l.length) :
    (l.set q f).getD q (Form.const false) = Form.const false := by
  rw [List.getD_eq_getElem?_getD, List.getElem?_eq_none]
  · rfl
  · simpa using Nat.le_of_not_lt h

/-- Writing one wire's parity keeps the state bounded, given the new parity is. -/
theorem bounded_set {st : AState} {q : Qubit} {f : Form} {fr : Nat}
    (hst : ∀ r, Form.Bounded fr (st.parOf r)) (hf : Form.Bounded fr f) :
    (⟨st.par.set q f, fr⟩ : AState).Bounded := by
  intro r
  rw [parOf_mk]
  by_cases hr : r = q
  · subst hr
    by_cases hlen : r < st.par.length
    · rw [getD_set_self _ _ _ hlen]; exact hf
    · rw [getD_set_out _ _ _ hlen]; exact Form.bounded_const _ false
  · rw [getD_set_ne _ _ hr]; exact hst r

theorem bounded_step {st : AState} (hst : st.Bounded) (g : Gate) : (st.step g).Bounded := by
  have hmono : ∀ q, Form.Bounded (st.step g).fresh (st.parOf q) := fun q =>
    Form.bounded_mono (fresh_le_step st g) (hst q)
  cases g with
  | x p => exact bounded_set hmono (Form.bounded_flip (hmono p))
  | cnot c t => exact bounded_set hmono (Form.bounded_add (hmono t) (hmono c))
  | h p => exact bounded_set hmono (Form.bounded_var (by simp))
  | ccx c₁ c₂ t => exact bounded_set hmono (Form.bounded_var (by simp))
  | reset p => exact bounded_set hmono (Form.bounded_var (by simp))
  | _ => exact hmono

theorem bounded_steps {st : AState} (hst : st.Bounded) (gs : List Gate) :
    (st.steps gs).Bounded := by
  induction gs generalizing st with
  | nil => exact hst
  | cons g gs ih => exact ih (bounded_step hst g)

theorem bounded_initial (n : Nat) : (initial n).Bounded := by
  intro q
  simp only [initial, parOf]
  by_cases hq : q < n
  · have : ((List.range n).map Form.var).getD q (Form.const false) = Form.var q := by
      rw [List.getD_eq_getElem?_getD, List.getElem?_map, List.getElem?_range hq]
      rfl
    rw [this]
    exact Form.bounded_var hq
  · have : ((List.range n).map Form.var).getD q (Form.const false) = Form.const false := by
      rw [List.getD_eq_getElem?_getD, List.getElem?_map,
        List.getElem?_eq_none (by simpa using Nat.le_of_not_lt hq)]
      rfl
    rw [this]
    exact Form.bounded_const _ false

/-! ## Consistency: what the parities mean -/

/-- A valuation `v` *explains* the basis state `b` at analysis state `st` when every wire's
value is the value of its parity. Out-of-range wires hold `false` on both sides, so the
condition is stated for all wires. -/
def Consistent (n : Nat) (st : AState) (v : Nat → Bool) (b : Basis n) : Prop :=
  ∀ q : Qubit, b.get q = (st.parOf q).evalB v

theorem initial_parOf (n : Nat) (q : Qubit) :
    (initial n).parOf q = if q < n then Form.var q else Form.const false := by
  simp only [initial, parOf, List.getD_eq_getElem?_getD, List.getElem?_map]
  by_cases hq : q < n
  · rw [List.getElem?_range hq, if_pos hq]; rfl
  · rw [List.getElem?_eq_none (by simpa using Nat.le_of_not_lt hq), if_neg hq]; rfl

theorem consistent_initial {n : Nat} (b : Basis n) :
    Consistent n (initial n) (fun i => b.get i) b := by
  intro q
  rw [initial_parOf]
  by_cases hq : q < n
  · rw [if_pos hq, Form.evalB_var]
  · rw [if_neg hq, Form.evalB_const, Basis.get, dif_neg hq]

/-- Consistency after writing one wire's parity. -/
theorem consistent_set {n : Nat} {st : AState} {fr : Nat} {q : Qubit} {f : Form}
    {v : Nat → Bool} {b' : Basis n} (hlen : st.par.length = n)
    (hq : q < n → b'.get q = f.evalB v)
    (hother : ∀ r, r ≠ q → b'.get r = (st.parOf r).evalB v) :
    Consistent n ⟨st.par.set q f, fr⟩ v b' := by
  intro r
  rw [parOf_mk]
  by_cases hr : r = q
  · subst hr
    by_cases hlt : r < st.par.length
    · rw [getD_set_self _ _ _ hlt]
      exact hq (hlen ▸ hlt)
    · rw [getD_set_out _ _ _ hlt, Form.evalB_const, Basis.get, dif_neg (hlen ▸ hlt)]
  · rw [getD_set_ne _ _ hr]
    exact hother r hr

/-! ## Entrywise facts about gate matrices -/

theorem perm_entry {n : Nat} {σ : Basis n → Basis n} {b' b : Basis n}
    (h : permMatrix σ b' b ≠ 0) : b' = σ b := by
  by_contra hc
  exact h (by simp [permMatrix, hc])

theorem diag_entry {n : Nat} {U : Density n} (hU : IsDiagonal U) {b' b : Basis n}
    (h : U b' b ≠ 0) : b' = b := by
  by_contra hc
  exact h (hU _ _ hc)

theorem one_entry {n : Nat} {b' b : Basis n} (h : (1 : Density n) b' b ≠ 0) : b' = b :=
  diag_entry isDiagonal_one h

/-- A one-wire gate leaves every other wire's value alone on any nonzero amplitude. -/
theorem embed1_entry {n : Nat} {M : Bool → Bool → ℂ} {q : Qubit} {b' b : Basis n}
    (h : embed1 n M q b' b ≠ 0) {r : Qubit} (hr : r ≠ q) : b'.get r = b.get r := by
  by_cases hq : q < n
  · simp only [embed1_apply_of_lt _ hq] at h
    by_cases hall : ∀ s : Fin n, (s : Nat) ≠ q → b' s = b s
    · by_cases hrn : r < n
      · have := hall ⟨r, hrn⟩ (by simpa using hr)
        simpa [Basis.get, hrn] using this
      · simp [Basis.get, hrn]
    · rw [if_neg hall] at h
      exact absurd rfl h
  · rw [embed1_eq_one _ hq] at h
    rw [one_entry h]

/-! ## Soundness of one step -/

/-- **One step is sound.** Any basis transition a unitary gate can make is explained by a
valuation extending the current one; only `h` and `ccx` need a new variable, and they get a
fresh one. -/
theorem step_sound {n : Nat} {st : AState} (hst : st.Bounded) (hlen : st.par.length = n)
    {g : Gate} (hu : g.isUnitary = true) {v : Nat → Bool} {b b' : Basis n}
    (hb : Consistent n st v b) (hne : gateUnitary n g b' b ≠ 0) :
    ∃ v', (∀ i, i < st.fresh → v' i = v i) ∧ Consistent n (st.step g) v' b' := by
  -- The diagonal gates: the state does not move and neither does the analysis.
  have diagCase : ∀ (U : Density n), IsDiagonal U → gateUnitary n g = U → st.step g = st →
      ∃ v', (∀ i, i < st.fresh → v' i = v i) ∧ Consistent n (st.step g) v' b' := by
    intro U hU hgU hstep
    refine ⟨v, fun i _ => rfl, ?_⟩
    have hbb : b' = b := diag_entry hU (by rw [← hgU]; exact hne)
    rw [hstep, hbb]
    exact hb
  cases g with
  | x q =>
      refine ⟨v, fun i _ => rfl, ?_⟩
      have hb' : b' = b.set q (!b.get q) :=
        perm_entry (σ := fun b : Basis n => b.set q (!b.get q))
          (by rw [← embed1_x2_eq_perm]; exact hne)
      show Consistent n ⟨st.par.set q (st.parOf q).flip, st.fresh⟩ v b'
      refine consistent_set hlen ?_ ?_
      · intro hqn
        rw [hb', Basis.get_set_same _ _ _ hqn, Form.evalB_flip, ← hb q]
      · intro r hr
        rw [hb', Basis.get_set_ne _ _ _ _ hr, hb r]
  | cnot c t =>
      refine ⟨v, fun i _ => rfl, ?_⟩
      have hb' : b' = b.set t (b.get t != b.get c) :=
        perm_entry (σ := fun b : Basis n => b.set t (b.get t != b.get c)) hne
      show Consistent n ⟨st.par.set t ((st.parOf t) + (st.parOf c)), st.fresh⟩ v b'
      refine consistent_set hlen ?_ ?_
      · intro htn
        rw [hb', Basis.get_set_same _ _ _ htn, Form.evalB_add, ← hb t, ← hb c]
      · intro r hr
        rw [hb', Basis.get_set_ne _ _ _ _ hr, hb r]
  | ccx c₁ c₂ t =>
      have hb' : b' = b.set t (b.get t != (b.get c₁ && b.get c₂)) :=
        perm_entry (σ := fun b : Basis n => b.set t (b.get t != (b.get c₁ && b.get c₂))) hne
      refine ⟨Function.update v st.fresh (b'.get t), fun i hi => Function.update_of_ne
        (by omega) _ _, ?_⟩
      show Consistent n ⟨st.par.set t (Form.var st.fresh), st.fresh + 1⟩ _ b'
      refine consistent_set hlen ?_ ?_
      · intro _
        rw [Form.evalB_var, Function.update_self]
      · intro r hr
        rw [hb', Basis.get_set_ne _ _ _ _ hr, hb r]
        exact (Form.evalB_congr (hst r) (fun i hi => Function.update_of_ne (by omega) _ _)).symm
  | h q =>
      refine ⟨Function.update v st.fresh (b'.get q), fun i hi => Function.update_of_ne
        (by omega) _ _, ?_⟩
      show Consistent n ⟨st.par.set q (Form.var st.fresh), st.fresh + 1⟩ _ b'
      refine consistent_set hlen ?_ ?_
      · intro _
        rw [Form.evalB_var, Function.update_self]
      · intro r hr
        rw [embed1_entry hne hr, hb r]
        exact (Form.evalB_congr (hst r) (fun i hi => Function.update_of_ne (by omega) _ _)).symm
  | s q => exact diagCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | sdg q => exact diagCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | z q => exact diagCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | t q => exact diagCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | tdg q => exact diagCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | rz θ q => exact diagCase _ (isDiagonal_embed1_diag2 _ _ q) rfl rfl
  | cz c t => exact diagCase _ (isDiagonal_phaseMatrix _) rfl rfl
  | ccz c₁ c₂ t => exact diagCase _ (isDiagonal_phaseMatrix _) rfl rfl
  | measure q c => simp [Gate.isUnitary, Gate.isMeasurement] at hu
  | reset q => simp [Gate.isUnitary, Gate.isMeasurement] at hu

/-- **The analysis is sound.** Every basis transition a measurement-free fragment can make is
explained by a valuation of the analysis' variables that extends the one describing the input.
This is the single semantic fact `PhaseFoldRand` rests on. -/
theorem analyze_sound {n : Nat} {gs : List Gate} (hu : ∀ g ∈ gs, g.isUnitary = true)
    {st : AState} (hst : st.Bounded) (hlen : st.par.length = n)
    {v : Nat → Bool} {b b' : Basis n}
    (hb : Consistent n st v b) (hne : unitary n gs b' b ≠ 0) :
    ∃ v', (∀ i, i < st.fresh → v' i = v i) ∧ Consistent n (st.steps gs) v' b' := by
  induction gs generalizing st v b with
  | nil =>
      have hbb : b' = b := one_entry (by simpa using hne)
      exact ⟨v, fun i _ => rfl, by rw [steps_nil, hbb]; exact hb⟩
  | cons g gs ih =>
      rw [unitary_cons, Matrix.mul_apply] at hne
      obtain ⟨k, -, hk⟩ := Finset.exists_ne_zero_of_sum_ne_zero hne
      have hk₁ : unitary n gs b' k ≠ 0 := left_ne_zero_of_mul hk
      have hk₂ : gateUnitary n g k b ≠ 0 := right_ne_zero_of_mul hk
      obtain ⟨v₁, hv₁, hc₁⟩ :=
        step_sound hst hlen (hu g (by simp)) hb hk₂
      obtain ⟨v', hv', hc'⟩ :=
        ih (fun g' hg' => hu g' (by simp [hg']))
          (bounded_step hst g) (by rw [length_step, hlen]) hc₁ hk₁
      refine ⟨v', fun i hi => ?_, by simpa using hc'⟩
      rw [hv' i (lt_of_lt_of_le hi (fresh_le_step st g)), hv₁ i hi]

/-! ## Genericity -/

/-- A state is *generic* when it explains every basis state: no wire's value is pinned down.
Every state reachable from `initial` is generic, which is what lets the fold's induction drop
the prefix it has already processed. -/
def Generic (n : Nat) (st : AState) : Prop := ∀ b : Basis n, ∃ v, Consistent n st v b

theorem generic_initial (n : Nat) : Generic n (initial n) :=
  fun b => ⟨fun i => b.get i, consistent_initial b⟩

/-- Allocating a fresh variable on one wire keeps a state generic. -/
theorem generic_fresh {n : Nat} {st : AState} (hst : st.Bounded) (hlen : st.par.length = n)
    (hgen : Generic n st) (q : Qubit) :
    Generic n ⟨st.par.set q (Form.var st.fresh), st.fresh + 1⟩ := by
  intro b'
  obtain ⟨v, hv⟩ := hgen b'
  refine ⟨Function.update v st.fresh (b'.get q), consistent_set hlen ?_ ?_⟩
  · intro _
    rw [Form.evalB_var, Function.update_self]
  · intro r _
    rw [hv r]
    exact (Form.evalB_congr (hst r) (fun i hi => Function.update_of_ne (by omega) _ _)).symm

theorem generic_step {n : Nat} {st : AState} (hst : st.Bounded) (hlen : st.par.length = n)
    (hgen : Generic n st) {g : Gate} (hwf : g.Wf) : Generic n (st.step g) := by
  cases g with
  | x q =>
      intro b'
      obtain ⟨v, hv⟩ := hgen (b'.set q (!b'.get q))
      refine ⟨v, ?_⟩
      show Consistent n ⟨st.par.set q (st.parOf q).flip, st.fresh⟩ v b'
      refine consistent_set hlen ?_ ?_
      · intro hqn
        rw [Form.evalB_flip, ← hv q, Basis.get_set_same _ _ _ hqn, Bool.not_not]
      · intro r hr
        rw [← hv r, Basis.get_set_ne _ _ _ _ hr]
  | cnot c t =>
      intro b'
      obtain ⟨v, hv⟩ := hgen (b'.set t (b'.get t != b'.get c))
      have hct : c ≠ t := by
        simpa [Gate.Wf] using hwf
      refine ⟨v, ?_⟩
      show Consistent n ⟨st.par.set t ((st.parOf t) + (st.parOf c)), st.fresh⟩ v b'
      refine consistent_set hlen ?_ ?_
      · intro htn
        rw [Form.evalB_add, ← hv t, ← hv c, Basis.get_set_same _ _ _ htn,
          Basis.get_set_ne _ _ _ _ hct]
        cases b'.get t <;> cases b'.get c <;> rfl
      · intro r hr
        rw [← hv r, Basis.get_set_ne _ _ _ _ hr]
  | h q => exact generic_fresh hst hlen hgen q
  | ccx c₁ c₂ t => exact generic_fresh hst hlen hgen t
  | reset q => exact generic_fresh hst hlen hgen q
  | measure q c => exact hgen
  | s q => exact hgen
  | sdg q => exact hgen
  | z q => exact hgen
  | t q => exact hgen
  | tdg q => exact hgen
  | rz θ q => exact hgen
  | cz c t => exact hgen
  | ccz c₁ c₂ t => exact hgen

end AState

end
end TzapLean
