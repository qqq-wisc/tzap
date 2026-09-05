import TzapLean.PhaseFold
import TzapLean.CnotMinProof

/-!
# `PhaseFoldRand`: correctness

The pass compares tags; the proof compares parities. This file ties the two together and
then proves that folding preserves the circuit's meaning **whenever the tags are faithful** —
that is, whenever no two distinct parities the run compares happen to hash alike.

* `Sim` — the tag state mirrors the symbolic state.
* `Faithful` — on the forms this run compares, equal hashes mean equal forms.
* `foldFrom_correct` — the fold is meaning-preserving under `Faithful`.

Everything random is confined to `Faithful`; `TzapLean.PhaseFoldRand` bounds its failure
probability with the collision bound from `TzapLean.Hash`.
-/

namespace TzapLean

open Form

/-! ## Tags simulate parities -/

@[simp] theorem hash_const_false {k : Nat} (draws : Draws k) :
    hash draws (Form.const false) = 0 := by
  funext j
  simp [hash, Form.eval, Form.const, bit]

/-- The tag state and the symbolic state agree, wire by wire. -/
def Sim {k : Nat} (draws : Draws k) (st : AState) (ts : TState k) : Prop :=
  ts.tags.length = st.par.length ∧ ts.fresh = st.fresh ∧
    ∀ q : Qubit, wordToBits (k := k) (ts.tagOf q) = hash draws (st.parOf q)

theorem getD_map_range {α : Type*} (f : Nat → α) (n q : Nat) (d : α) :
    ((List.range n).map f).getD q d = if q < n then f q else d := by
  by_cases hq : q < n
  · have hget : ((List.range n).map f)[q]? = some (f q) := by
      rw [List.getElem?_map, List.getElem?_range hq]; rfl
    rw [List.getD_eq_getElem?_getD, hget, if_pos hq]; rfl
  · have hget : ((List.range n).map f)[q]? = none := by
      refine List.getElem?_eq_none ?_
      rw [List.length_map, List.length_range]
      omega
    rw [List.getD_eq_getElem?_getD, hget, if_neg hq]; rfl

theorem hash_flip {k : Nat} (draws : Draws k) (p : Form) :
    hash draws p.flip = hash draws p + ones k := by
  funext j
  simp only [hash, Form.flip, Form.const, bit, ones, Form.eval, Pi.add_apply,
    Form.add_constant, Form.add_coefficients, if_pos]
  ring_nf

theorem sim_initial {k : Nat} {draws : Draws k} {wdraws : Nat → Tag}
    (hw : ∀ i, wordToBits (k := k) (wdraws i) = draws i) (n : Nat) :
    Sim draws (AState.initial n) (TState.initial (k := k) wdraws n) := by
  refine ⟨by simp [TState.initial, AState.initial], rfl, fun q => ?_⟩
  rw [AState.initial_parOf, TState.tagOf, TState.initial, getD_map_range]
  by_cases hq : q < n
  · rw [if_pos hq, if_pos hq, hw, hash_var]
  · rw [if_neg hq, if_neg hq, wordToBits_zero, hash_const_false]

/-- Writing one wire keeps the simulation, given the tag written is the parity's hash. -/
theorem sim_set {k : Nat} {draws : Draws k} {st : AState} {ts : TState k}
    (hsim : Sim draws st ts) (q : Qubit) (f : Form) (t : Tag)
    (ht : wordToBits (k := k) t = hash draws f) (fr fr' : Nat) (hfr : fr = fr') :
    Sim draws ⟨st.par.set q f, fr'⟩ ⟨ts.tags.set q t, fr⟩ := by
  obtain ⟨hlen, -, htag⟩ := hsim
  refine ⟨by simp [hlen], hfr, fun r => ?_⟩
  simp only [TState.tagOf, AState.parOf]
  by_cases hr : r = q
  · subst hr
    by_cases hlt : r < ts.tags.length
    · rw [List.getD_eq_getElem?_getD, List.getElem?_set_self hlt,
        AState.getD_set_self _ _ _ (hlen ▸ hlt)]
      exact ht
    · rw [List.getD_eq_getElem?_getD, List.getElem?_eq_none (by simpa using Nat.le_of_not_lt hlt),
        AState.getD_set_out _ _ _ (hlen ▸ hlt), hash_const_false]
      exact wordToBits_zero
  · rw [List.getD_eq_getElem?_getD, List.getElem?_set, if_neg (by simpa using Ne.symm hr),
      ← List.getD_eq_getElem?_getD, AState.getD_set_ne _ _ hr]
    exact htag r

theorem sim_step {k : Nat} {draws : Draws k} {wdraws : Nat → Tag} {st : AState} {ts : TState k}
    (hw : ∀ i, wordToBits (k := k) (wdraws i) = draws i)
    (hsim : Sim draws st ts) (g : Gate) : Sim draws (st.step g) (ts.step wdraws g) := by
  have hfresh : ts.fresh = st.fresh := hsim.2.1
  have htag := hsim.2.2
  cases g with
  | x q =>
      refine sim_set hsim q _ _ ?_ _ _ hfresh
      rw [wordToBits_xor, htag q, wordToBits_onesTag, hash_flip]
  | cnot c t =>
      refine sim_set hsim t _ _ ?_ _ _ hfresh
      rw [wordToBits_xor, htag t, htag c, hash_add]
  | h q => exact sim_set hsim q _ _ (by rw [hw, hash_var, hfresh]) _ _ (by rw [hfresh])
  | ccx c₁ c₂ t => exact sim_set hsim t _ _ (by rw [hw, hash_var, hfresh]) _ _ (by rw [hfresh])
  | reset q => exact sim_set hsim q _ _ (by rw [hw, hash_var, hfresh]) _ _ (by rw [hfresh])
  | _ => exact hsim

theorem sim_steps {k : Nat} {draws : Draws k} {wdraws : Nat → Tag} {st : AState} {ts : TState k}
    (hw : ∀ i, wordToBits (k := k) (wdraws i) = draws i)
    (hsim : Sim draws st ts) (gs : List Gate) :
    Sim draws (st.steps gs) (ts.steps wdraws gs) := by
  induction gs generalizing st ts with
  | nil => exact hsim
  | cons g gs ih => exact ih (sim_step hw hsim g)

/-! ## The forms a run compares -/

/-- Every wire's parity at one point in the run, plus the constant parity that stands for
wires outside the register. -/
noncomputable def formsOf (n : Nat) (st : AState) : List Form :=
  Form.const false :: (List.range n).map st.parOf

/-- Every wire's parity at every point in the run. -/
noncomputable def visited (n : Nat) (st : AState) : List Gate → List Form
  | [] => formsOf n st
  | g :: gs => formsOf n st ++ visited n (st.step g) gs

/-- The forms the fold may compare: the parities it sees and their complements. -/
noncomputable def relevant (n : Nat) (st : AState) (gs : List Gate) : List Form :=
  visited n st gs ++ (visited n st gs).map Form.flip

/-- The randomness assumption: on the forms this run compares, equal tags mean equal
parities. Its failure probability is bounded in `TzapLean.PhaseFoldRand`. -/
def Faithful {k : Nat} (draws : Draws k) (ps : List Form) : Prop :=
  ∀ p ∈ ps, ∀ q ∈ ps, hash draws p = hash draws q → p = q

theorem Faithful.mono {k : Nat} {draws : Draws k} {ps qs : List Form}
    (h : Faithful draws ps) (hsub : ∀ p ∈ qs, p ∈ ps) : Faithful draws qs :=
  fun p hp q hq hpq => h p (hsub p hp) q (hsub q hq) hpq

theorem parOf_of_ge {st : AState} {q : Qubit} (h : ¬ q < st.par.length) :
    st.parOf q = Form.const false := by
  rw [AState.parOf, List.getD_eq_getElem?_getD,
    List.getElem?_eq_none (Nat.le_of_not_lt h)]
  rfl

theorem mem_formsOf {n : Nat} {st : AState} (hlen : st.par.length = n) (q : Qubit) :
    st.parOf q ∈ formsOf n st := by
  by_cases hq : q < n
  · exact List.mem_cons_of_mem _ (List.mem_map.2 ⟨q, List.mem_range.2 hq, rfl⟩)
  · rw [parOf_of_ge (by rw [hlen]; exact hq)]
    exact List.mem_cons_self

theorem formsOf_sub_visited {n : Nat} {st : AState} {p : Form} (gs : List Gate)
    (h : p ∈ formsOf n st) : p ∈ visited n st gs := by
  cases gs with
  | nil => exact h
  | cons g gs => exact List.mem_append.2 (Or.inl h)

theorem mem_visited_append {n : Nat} {st : AState} {p : Form} (as bs : List Gate)
    (h : p ∈ visited n (st.steps as) bs) : p ∈ visited n st (as ++ bs) := by
  induction as generalizing st with
  | nil => exact h
  | cons a as ih => exact List.mem_append.2 (Or.inr (ih h))

theorem visited_append_sub {n : Nat} (as bs cs : List Gate) (st : AState)
    (h : ∀ p ∈ visited n (st.steps as) bs, p ∈ visited n (st.steps as) cs) :
    ∀ p ∈ visited n st (as ++ bs), p ∈ visited n st (as ++ cs) := by
  induction as generalizing st with
  | nil => exact h
  | cons a as ih =>
      intro p hp
      rcases List.mem_append.1 hp with hp | hp
      · exact List.mem_append.2 (Or.inl hp)
      · exact List.mem_append.2 (Or.inr (ih (st.step a) h p hp))

/-- A rotation gate does not move the analysis. -/
theorem step_of_rotAngle {st : AState} {g : Gate} {θ : ℚ} {q : Qubit}
    (h : rotAngle g = some (θ, q)) : st.step g = st := by
  cases g <;> simp [rotAngle] at h ⊢ <;> rfl

/-! ## Rotations, re-emitted -/

theorem phaseMatrix_trivial {n : Nat} {q : Qubit} (hq : ¬ q < n) (z : ℂ) :
    phaseMatrix (fun s : Basis n => if s.get q then z else 1) = 1 := by
  rw [← phaseMatrix_one]
  congr 1
  funext s
  rw [basis_get_of_ge s hq, if_neg (by simp)]

/-- Any of the gates phase folding tracks is a phase on its wire, up to a unit factor —
in range or not. -/
theorem phase_of_rotAngle {n : Nat} {g : Gate} {θ : ℚ} {q : Qubit}
    (h : rotAngle g = some (θ, q)) :
    ∃ c : ℂ, c * star c = 1 ∧
      gateUnitary n g = c • phaseMatrix (fun s : Basis n => if s.get q then ep θ else 1) := by
  by_cases hq : q < n
  · exact rotAngle_gateUnitary h hq
  · refine ⟨1, by simp, ?_⟩
    rw [phaseMatrix_trivial hq, one_smul]
    cases g with
    | t p => simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
             obtain ⟨rfl, rfl⟩ := h
             simp [gateUnitary, embed1, hq]
    | tdg p => simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
               obtain ⟨rfl, rfl⟩ := h
               simp [gateUnitary, embed1, hq]
    | s p => simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
             obtain ⟨rfl, rfl⟩ := h
             simp [gateUnitary, embed1, hq]
    | sdg p => simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
               obtain ⟨rfl, rfl⟩ := h
               simp [gateUnitary, embed1, hq]
    | z p => simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
             obtain ⟨rfl, rfl⟩ := h
             simp [gateUnitary, embed1, hq]
    | rz b p => simp only [rotAngle, Option.some.injEq, Prod.mk.injEq] at h
                obtain ⟨rfl, rfl⟩ := h
                simp [gateUnitary, embed1, hq]
    | _ => simp [rotAngle] at h

theorem emitRotation_eq_nil {q : Qubit} {a : ℚ} (h0 : BlockState.angleMod a = 0) :
    emitRotation q a = [] := by
  simp [emitRotation, h0]

theorem emitRotation_eq_diagRun {q : Qubit} {a : ℚ} {j : Nat}
    (h0 : BlockState.angleMod a ≠ 0) (hcl : classifyQuarterPi (BlockState.angleMod a) = some j) :
    emitRotation q a = diagRun j q := by
  have hden : (4 * BlockState.angleMod a).den = 1 := by
    by_contra hc
    simp [classifyQuarterPi, hc] at hcl
  have hj : ((((4 * BlockState.angleMod a).num % 8 + 8) % 8).toNat) = j := by
    simpa [classifyQuarterPi, hden] using hcl
  have h8 : j < 8 := classifyQuarterPi_lt hcl
  simp only [emitRotation, beq_iff_eq, if_neg h0, hden, hj]
  interval_cases j <;> rfl

theorem emitRotation_eq_rz {q : Qubit} {a : ℚ}
    (h0 : BlockState.angleMod a ≠ 0) (hcl : classifyQuarterPi (BlockState.angleMod a) = none) :
    emitRotation q a = [Gate.rz (BlockState.angleMod a) q] := by
  have hden : (4 * BlockState.angleMod a).den ≠ 1 := by
    by_contra hc
    simp [classifyQuarterPi, hc] at hcl
  simp [emitRotation, h0, hden]

/-- Re-emitting a rotation reproduces its phase, up to a unit factor. -/
theorem unitary_emitRotation {n : Nat} (q : Qubit) (a : ℚ) :
    ∃ c : ℂ, c * star c = 1 ∧
      unitary n (emitRotation q a) =
        c • phaseMatrix (fun s : Basis n => if s.get q then ep a else 1) := by
  have hmod : ep (BlockState.angleMod a) = ep a := ep_angleMod a
  by_cases h0 : BlockState.angleMod a = 0
  · have hep : ep a = 1 := by rw [← hmod, h0, ep_zero]
    refine ⟨1, by simp, ?_⟩
    rw [emitRotation_eq_nil h0, unitary_nil, one_smul, hep, ← phaseMatrix_one]
    congr 1
    funext s
    by_cases hs : s.get q <;> simp [hs]
  · cases hcl : classifyQuarterPi (BlockState.angleMod a) with
    | some j =>
        have h8 : j < 8 := classifyQuarterPi_lt hcl
        refine ⟨1, by simp, ?_⟩
        rw [emitRotation_eq_diagRun h0 hcl, unitary_diagRun n q j h8, one_smul]
        by_cases hq : q < n
        · rw [embed1_diag2_eq_phaseMatrix _ _ _ hq, ← ep_of_classifyQuarterPi hcl, hmod]
        · rw [phaseMatrix_trivial hq]
          simp [embed1, hq]
    | none =>
        rw [emitRotation_eq_rz h0 hcl]
        obtain ⟨c, hc, hgu⟩ :=
          phase_of_rotAngle (n := n) (g := Gate.rz (BlockState.angleMod a) q)
            (θ := BlockState.angleMod a) (q := q) rfl
        exact ⟨c, hc, by rw [unitary_cons, unitary_nil, Matrix.one_mul, hgu, hmod]⟩

theorem rotAngle_isUnitary {g : Gate} {a : ℚ} {q : Qubit} (h : rotAngle g = some (a, q)) :
    g.isUnitary = true := by
  cases g <;> simp [rotAngle] at h ⊢ <;> rfl

theorem diagRun_isUnitary (j : Nat) (q : Qubit) : ∀ g ∈ diagRun j q, g.isUnitary = true := by
  match j with
  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 => intro g hg; fin_cases hg <;> rfl
  | (i + 8) => intro g hg; simp [diagRun] at hg

theorem emitRotation_isUnitary (q : Qubit) (a : ℚ) :
    ∀ g ∈ emitRotation q a, g.isUnitary = true := by
  by_cases h0 : BlockState.angleMod a = 0
  · rw [emitRotation_eq_nil h0]; simp
  · cases hcl : classifyQuarterPi (BlockState.angleMod a) with
    | some j => rw [emitRotation_eq_diagRun h0 hcl]; exact diagRun_isUnitary j q
    | none => rw [emitRotation_eq_rz h0 hcl]; intro g hg; rw [List.mem_singleton.1 hg]; rfl

/-- The pass may replace any tracked rotation with its re-emission. -/
theorem equivalent_emitRotation {n m : Nat} {g : Gate} {a : ℚ} {q : Qubit}
    (h : rotAngle g = some (a, q)) : Equivalent n m (emitRotation q a) [g] := by
  obtain ⟨c, hc, hgu⟩ := phase_of_rotAngle (n := n) h
  obtain ⟨d, hd, hem⟩ := unitary_emitRotation (n := n) q a
  have hcs : star c * star (star c) = 1 := by rw [star_star, mul_comm]; exact hc
  refine equivalent_of_unitary_smul (emitRotation_isUnitary q a) ?_ (d * star c)
    (unit_mul hd hcs) ?_
  · intro x hx
    rw [List.mem_singleton.1 hx]
    exact rotAngle_isUnitary h
  · rw [unitary_cons, unitary_nil, Matrix.one_mul, hem, hgu, smul_smul]
    congr 1
    rw [mul_assoc, show star c * c = 1 from by rw [mul_comm]; exact hc, mul_one]

/-- Re-emitting every rotation preserves meaning. -/
theorem emitAll_correct {n m : Nat} (gs : List Gate) : Equivalent n m (emitAll gs) gs := by
  induction gs with
  | nil => exact fun _ => rfl
  | cons g gs ih =>
      have hhead : Equivalent n m (match rotAngle g with
          | some (a, q) => emitRotation q a
          | none => [g]) [g] := by
        cases hrot : rotAngle g with
        | some p =>
            obtain ⟨a, q⟩ := p
            simpa [hrot] using equivalent_emitRotation (n := n) (m := m) hrot
        | none => simp; exact fun _ => rfl
      show Equivalent n m ((match rotAngle g with
          | some (a, q) => emitRotation q a
          | none => [g]) ++ emitAll gs) (g :: gs)
      exact Equivalent.trans (Equivalent.append_right _ hhead)
        (by simpa using Equivalent.append_left [g] ih)

/-! ## What a successful merge means -/

theorem mergeInto_spec {k : Nat} (wdraws : Nat → Tag) (tag : Tag) (θ : ℚ) :
    ∀ (gs gs' : List Gate) (ts : TState k), mergeInto wdraws ts tag θ gs = some gs' →
      ∃ (M rest : List Gate) (g' : Gate) (φ : ℚ) (q' : Qubit) (sign : Bool),
        gs = M ++ g' :: rest ∧
        gs' = M ++ Gate.rz (φ + signedAngle sign θ) q' :: rest ∧
        (∀ g ∈ M, g.isUnitary = true) ∧
        rotAngle g' = some (φ, q') ∧
        matchTag k tag ((ts.steps wdraws M).tagOf q') = some sign := by
  intro gs
  induction gs with
  | nil => intro gs' ts h; simp [mergeInto] at h
  | cons g gs ih =>
      intro gs' ts h
      simp only [mergeInto] at h
      split at h
      · rename_i hu
        split at h
        · rename_i φ q' hrot
          split at h
          · rename_i sign hmatch
            simp only [Option.some.injEq] at h
            subst h
            exact ⟨[], gs, g, φ, q', sign, rfl, rfl, by simp, hrot, hmatch⟩
          · rcases Option.map_eq_some_iff.1 h with ⟨t, ht, rfl⟩
            obtain ⟨M, rest, g', φ', q'', sign, h1, h2, h3, h4, h5⟩ := ih t _ ht
            exact ⟨g :: M, rest, g', φ', q'', sign, by rw [h1]; rfl, by rw [h2]; rfl,
              by intro x hx; rcases List.mem_cons.1 hx with rfl | hx
                 · exact hu
                 · exact h3 x hx,
              h4, h5⟩
        · rcases Option.map_eq_some_iff.1 h with ⟨t, ht, rfl⟩
          obtain ⟨M, rest, g', φ', q'', sign, h1, h2, h3, h4, h5⟩ := ih t _ ht
          exact ⟨g :: M, rest, g', φ', q'', sign, by rw [h1]; rfl, by rw [h2]; rfl,
            by intro x hx; rcases List.mem_cons.1 hx with rfl | hx
               · exact hu
               · exact h3 x hx,
            h4, h5⟩
      · exact absurd h (by simp)

/-! ## Unfolding the fold -/

@[simp] theorem foldFrom_nil {k : Nat} (wdraws : Nat → Tag) (targets : Array Bool)
    (ts : TState k) (at_ : Nat) : foldFrom wdraws targets ts at_ [] = [] := by
  rw [foldFrom]

theorem foldFrom_cons_merge {k : Nat} {wdraws : Nat → Tag} {targets : Array Bool}
    {ts : TState k} {at_ : Nat} {g : Gate} {θ : ℚ} {q : Qubit} {gs gs' : List Gate}
    (hrot : rotAngle g = some (θ, q)) (hsel : targets[at_]?.getD true = true)
    (hm : mergeInto wdraws ts (ts.tagOf q) θ gs = some gs') :
    foldFrom wdraws targets ts at_ (g :: gs) = foldFrom wdraws targets ts (at_ + 1) gs' := by
  rw [foldFrom]
  simp only [hrot, hsel, if_true]
  split
  · rename_i x heq
    rw [hm] at heq
    cases heq
    rfl
  · rename_i heq
    rw [hm] at heq
    exact absurd heq (by simp)

theorem foldFrom_cons_none {k : Nat} {wdraws : Nat → Tag} {targets : Array Bool}
    {ts : TState k} {at_ : Nat} {g : Gate} {gs : List Gate} (hrot : rotAngle g = none) :
    foldFrom wdraws targets ts at_ (g :: gs)
      = g :: foldFrom wdraws targets (ts.step wdraws g) (at_ + 1) gs := by
  rw [foldFrom]
  simp only [hrot]

/-- The two ways the fold declines to merge: the filter said there is nothing to merge with,
or the scan looked and found nothing. -/
theorem foldFrom_cons_keep {k : Nat} {wdraws : Nat → Tag} {targets : Array Bool}
    {ts : TState k} {at_ : Nat} {g : Gate} {θ : ℚ} {q : Qubit} {gs : List Gate}
    (hrot : rotAngle g = some (θ, q))
    (hm : targets[at_]?.getD true = false ∨ mergeInto wdraws ts (ts.tagOf q) θ gs = none) :
    foldFrom wdraws targets ts at_ (g :: gs)
      = g :: foldFrom wdraws targets (ts.step wdraws g) (at_ + 1) gs := by
  rw [foldFrom]
  simp only [hrot]
  by_cases hsel : targets[at_]?.getD true = true
  · rw [if_pos hsel]
    have hmm : mergeInto wdraws ts (ts.tagOf q) θ gs = none := by
      rcases hm with h | h
      · rw [h] at hsel; exact absurd hsel (by simp)
      · exact h
    split
    · rename_i x heq
      rw [hmm] at heq
      exact absurd heq (by simp)
    · rfl
  · rw [if_neg hsel]

/-! ## Rotations are `rz`s -/

/-- Every gate phase folding tracks is an `rz` of its angle, up to global phase. -/
theorem equivalent_rot_rz {n m : Nat} {g : Gate} {θ : ℚ} {q : Qubit}
    (h : rotAngle g = some (θ, q)) : Equivalent n m [g] [Gate.rz θ q] := by
  obtain ⟨c, hc, hgu⟩ := phase_of_rotAngle (n := n) h
  obtain ⟨d, hd, hrz⟩ := phase_of_rotAngle (n := n) (g := Gate.rz θ q) (θ := θ) (q := q) rfl
  have hds : star d * star (star d) = 1 := by rw [star_star, mul_comm]; exact hd
  refine equivalent_of_unitary_smul ?_ ?_ (c * star d) (unit_mul hc hds) ?_
  · intro x hx
    rw [List.mem_singleton.1 hx]
    exact rotAngle_isUnitary h
  · intro x hx
    rw [List.mem_singleton.1 hx]
    rfl
  · rw [unitary_cons, unitary_nil, Matrix.one_mul, unitary_cons, unitary_nil, Matrix.one_mul,
      hgu, hrz, smul_smul]
    congr 1
    rw [mul_assoc, show star d * d = 1 from by rw [mul_comm]; exact hd, mul_one]

/-! ## Faithfulness bookkeeping -/

theorem mem_visited_of_parOf {n : Nat} {st : AState} (hlen : st.par.length = n)
    (q : Qubit) (gs : List Gate) : st.parOf q ∈ visited n st gs :=
  formsOf_sub_visited gs (mem_formsOf hlen q)

theorem visited_cons_rot {n : Nat} {st : AState} {g g' : Gate} {θ θ' : ℚ} {q q' : Qubit}
    (h : rotAngle g = some (θ, q)) (h' : rotAngle g' = some (θ', q')) (gs : List Gate) :
    visited n st (g :: gs) = visited n st (g' :: gs) := by
  rw [visited, visited, step_of_rotAngle h, step_of_rotAngle h']

theorem faithful_of_sub {k n : Nat} {draws : Draws k} {st st' : AState} {gs hs : List Gate}
    (hf : Faithful draws (relevant n st gs))
    (hsub : ∀ p ∈ visited n st' hs, p ∈ visited n st gs) :
    Faithful draws (relevant n st' hs) := by
  refine hf.mono ?_
  intro p hp
  rcases List.mem_append.1 hp with hp | hp
  · exact List.mem_append.2 (Or.inl (hsub p hp))
  · rcases List.mem_map.1 hp with ⟨r, hr, rfl⟩
    exact List.mem_append.2 (Or.inr (List.mem_map.2 ⟨r, hsub r hr, rfl⟩))

/-! ## The fold is correct -/

/-- **The fold preserves meaning.** Under `Faithful` — no two distinct parities compared by
this run hash alike — folding a gate list yields an equivalent one. -/
theorem foldFrom_correct {n m k : Nat} {draws : Draws k} {wdraws : Nat → Tag}
    (hw : ∀ i, wordToBits (k := k) (wdraws i) = draws i) :
    ∀ (N : Nat) (gs : List Gate), gs.length ≤ N → ∀ (at_ : Nat) (st : AState) (ts : TState k),
      Sim draws st ts → st.Bounded → st.par.length = n → AState.Generic n st →
      (∀ g ∈ gs, g.Wf) → Faithful draws (relevant n st gs) →
      Equivalent n m (foldFrom wdraws targets ts at_ gs) gs := by
  intro N
  induction N with
  | zero =>
      intro gs hgs at_ st ts _ _ _ _ _ _
      have : gs = [] := List.eq_nil_of_length_eq_zero (Nat.le_zero.1 hgs)
      subst this
      simp only [foldFrom_nil]
      exact fun _ => rfl
  | succ N ih =>
      intro gs hlen' at_ st ts hsim hbd hlen hgen hwf hfaith
      cases gs with
      | nil => simp only [foldFrom_nil]; exact fun _ => rfl
      | cons g gs =>
          have hlenN : gs.length ≤ N := by
            simp only [List.length_cons] at hlen'
            omega
          -- The two "keep the gate" branches share a proof.
          have keep : foldFrom wdraws targets ts at_ (g :: gs)
                = g :: foldFrom wdraws targets (ts.step wdraws g) (at_ + 1) gs →
              Equivalent n m (foldFrom wdraws targets ts at_ (g :: gs)) (g :: gs) := by
            intro heq
            rw [heq]
            have hsub : ∀ p ∈ visited n (st.step g) gs, p ∈ visited n st (g :: gs) := by
              intro p hp
              exact List.mem_append.2 (Or.inr hp)
            have := ih gs hlenN (at_ + 1) (st.step g) (ts.step wdraws g) (sim_step hw hsim g)
              (AState.bounded_step hbd g) (by rw [AState.length_step, hlen])
              (AState.generic_step hbd hlen hgen (hwf g (by simp)))
              (fun x hx => hwf x (by simp [hx])) (faithful_of_sub hfaith hsub)
            simpa using Equivalent.append_left [g] this
          cases hrot : rotAngle g with
          | none => exact keep (foldFrom_cons_none hrot)
          | some p =>
              obtain ⟨θ, q⟩ := p
              -- The filter can decline outright; otherwise the scan decides.
              by_cases hsel : targets[at_]?.getD true = true
              case neg => exact keep (foldFrom_cons_keep hrot (Or.inl (by simpa using hsel)))
              case pos =>
              cases hm : mergeInto wdraws ts (ts.tagOf q) θ gs with
              | none => exact keep (foldFrom_cons_keep hrot (Or.inr hm))
              | some gs' =>
                  obtain ⟨M, rest, g', φ, q', sign, hgseq, hgs'eq, hMu, hrot', hmatch⟩ :=
                    mergeInto_spec wdraws (ts.tagOf q) θ gs gs' ts hm
                  -- The recursive call sees the same state and a no-longer list.
                  have hlen'' : gs'.length ≤ N := by
                    have := mergeInto_length wdraws (ts.tagOf q) θ gs gs' ts hm
                    omega
                  have hwf' : ∀ x ∈ gs', x.Wf := by
                    intro x hx
                    rw [hgs'eq] at hx
                    rcases List.mem_append.1 hx with hx | hx
                    · exact hwf x (by rw [hgseq]; simp [hx])
                    · rcases List.mem_cons.1 hx with rfl | hx
                      · trivial
                      · exact hwf x (by rw [hgseq]; simp [hx])
                  have hvis : ∀ p ∈ visited n st gs', p ∈ visited n st (g :: gs) := by
                    intro p hp
                    refine List.mem_append.2 (Or.inr ?_)
                    rw [step_of_rotAngle hrot]
                    rw [hgs'eq] at hp
                    rw [hgseq]
                    refine visited_append_sub M _ _ st ?_ p hp
                    intro r hr
                    rwa [visited_cons_rot (g := Gate.rz (φ + signedAngle sign θ) q')
                      (g' := g') rfl hrot'] at hr
                  have hIH : Equivalent n m (foldFrom wdraws targets ts (at_ + 1) gs') gs' :=
                    ih gs' hlen'' (at_ + 1) st ts hsim hbd hlen hgen hwf'
                      (faithful_of_sub hfaith hvis)
                  -- The parity condition, from the tag match and faithfulness.
                  have hsimM : Sim draws (st.steps M) (ts.steps wdraws M) := sim_steps hw hsim M
                  have hmemq : st.parOf q ∈ visited n st (g :: gs) :=
                    mem_visited_of_parOf hlen q (g :: gs)
                  have hmemq' : (st.steps M).parOf q' ∈ visited n st (g :: gs) := by
                    refine List.mem_append.2 (Or.inr ?_)
                    rw [step_of_rotAngle hrot, hgseq]
                    exact mem_visited_append M (g' :: rest)
                      (mem_visited_of_parOf (by rw [AState.length_steps, hlen]) q' (g' :: rest))
                  have hpar : (st.steps M).parOf q' =
                      if sign then (st.parOf q).flip else st.parOf q := by
                    simp only [matchTag] at hmatch
                    split at hmatch
                    · rename_i heq
                      have hsign : sign = false := by simpa using hmatch.symm
                      subst hsign
                      rw [if_neg (by simp)]
                      refine hfaith _ (List.mem_append.2 (Or.inl hmemq')) _
                        (List.mem_append.2 (Or.inl hmemq)) ?_
                      rw [← hsimM.2.2 q', ← hsim.2.2 q]
                      exact wordToBits_congr (by simpa using heq)
                    · split at hmatch
                      · rename_i heq
                        have hsign : sign = true := by simpa using hmatch.symm
                        subst hsign
                        rw [if_pos rfl]
                        refine hfaith _ (List.mem_append.2 (Or.inl hmemq')) _
                          (List.mem_append.2 (Or.inr (List.mem_map.2 ⟨_, hmemq, rfl⟩))) ?_
                        rw [hash_flip, ← hsimM.2.2 q', ← hsim.2.2 q, ← wordToBits_onesTag k,
                          ← wordToBits_xor]
                        exact wordToBits_congr (by simpa using heq)
                      · exact absurd hmatch (by simp)
                  have hpath := path_of_generic hMu hbd hlen hgen hpar
                  -- The merge itself.
                  have hmerge : Equivalent n m gs' (g :: gs) := by
                    have hstep₁ : Equivalent n m (g :: (M ++ g' :: rest))
                        (Gate.rz θ q :: (M ++ g' :: rest)) := by
                      simpa using
                        Equivalent.window (n := n) (m := m) [] (M ++ g' :: rest)
                          (equivalent_rot_rz hrot)
                    have hstep₂ : Equivalent n m (Gate.rz θ q :: (M ++ g' :: rest))
                        (Gate.rz θ q :: (M ++ Gate.rz φ q' :: rest)) := by
                      have := Equivalent.window (n := n) (m := m) (Gate.rz θ q :: M) rest
                        (equivalent_rot_rz (g := g') hrot')
                      simpa using this
                    have hstep₃ : Equivalent n m (Gate.rz θ q :: (M ++ Gate.rz φ q' :: rest))
                        (M ++ Gate.rz (φ + signedAngle sign θ) q' :: rest) := by
                      have hmg := merge_equivalent (n := n) (m := m) [] M θ φ q q' sign
                        (by simp) hMu hpath
                      have := Equivalent.append_right (n := n) (m := m) rest hmg
                      simpa using this
                    rw [hgs'eq, hgseq]
                    exact (hstep₁.trans (hstep₂.trans hstep₃)).symm
                  rw [foldFrom_cons_merge hrot hsel hm]
                  exact hIH.trans hmerge

/-- **Phase folding preserves meaning**, under `Faithful`. -/
theorem phaseFoldGates_correct {n m k : Nat} {draws : Draws k} {wdraws : Nat → Tag}
    (hw : ∀ i, wordToBits (k := k) (wdraws i) = draws i) (gs : List Gate)
    (hwf : ∀ g ∈ gs, g.Wf)
    (hf : Faithful draws (relevant n (AState.initial n) gs)) :
    Equivalent n m (phaseFoldGates k wdraws n gs) gs := by
  refine Equivalent.trans (emitAll_correct _) ?_
  exact foldFrom_correct hw gs.length gs le_rfl 0 (AState.initial n)
    (TState.initial (k := k) wdraws n)
    (sim_initial hw n) (AState.bounded_initial n) (length_initial_par n)
    (AState.generic_initial n) hwf hf

end TzapLean
