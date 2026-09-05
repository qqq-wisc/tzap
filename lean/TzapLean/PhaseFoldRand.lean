import TzapLean.PhaseFoldProof
import TzapLean.RandPass

/-!
# `PhaseFoldRand` as a `RandPass`

Everything is in place: `phaseFoldGates_correct` says the pass is right whenever its tags are
faithful, and `collides_probability_le` says unfaithful tags are unlikely. Putting the two
together gives the obligation `RandPass` demands,

```
Pr_{s ← uniform} [ ⟦phaseFold s c⟧ ≠ ⟦c⟧ ]  ≤  C(L, 2) · 2⁻ᵏ
```

where `L` is the number of parities (and complements) this circuit makes the pass compare.
The seed is one ideal uniform `k`-bit tag per variable — `Sample (varBound c) k`.
`phaseFoldIO_run` records the pure correspondence for callers that choose this probabilistic
variant. The CLI uses the collision-free `PhaseFoldExact` defined below instead.
-/

namespace TzapLean

open scoped ENNReal

open Form

/-! ## Well-formedness is preserved -/

theorem emitRotation_wf (q : Qubit) (a : ℚ) : ∀ g ∈ emitRotation q a, g.Wf := by
  by_cases h0 : BlockState.angleMod a = 0
  · rw [emitRotation_eq_nil h0]; simp
  · cases hcl : classifyQuarterPi (BlockState.angleMod a) with
    | some j =>
        rw [emitRotation_eq_diagRun h0 hcl]
        match j with
        | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 => intro g hg; fin_cases hg <;> trivial
        | (i + 8) => intro g hg; simp [diagRun] at hg
    | none =>
        rw [emitRotation_eq_rz h0 hcl]
        intro g hg
        rw [List.mem_singleton.1 hg]
        trivial

theorem emitAll_wf {gs : List Gate} (h : ∀ g ∈ gs, g.Wf) : ∀ g ∈ emitAll gs, g.Wf := by
  induction gs with
  | nil => intro g hg; simp [emitAll] at hg
  | cons x gs ih =>
      intro g hg
      rw [emitAll] at hg
      rcases List.mem_append.1 hg with hg | hg
      · cases hrot : rotAngle x with
        | some p =>
            obtain ⟨a, q⟩ := p
            rw [hrot] at hg
            exact emitRotation_wf q a g hg
        | none =>
            rw [hrot] at hg
            rw [List.mem_singleton.1 hg]
            exact h x (by simp)
      · exact ih (fun y hy => h y (by simp [hy])) g hg

theorem foldFrom_wf {k : Nat} (wdraws : Nat → Tag) (targets : Array Bool) :
    ∀ (N : Nat) (gs : List Gate), gs.length ≤ N → ∀ (at_ : Nat) (ts : TState k),
      (∀ g ∈ gs, g.Wf) → ∀ g ∈ foldFrom wdraws targets ts at_ gs, g.Wf := by
  intro N
  induction N with
  | zero =>
      intro gs hgs at_ ts _ g hg
      rw [List.eq_nil_of_length_eq_zero (Nat.le_zero.1 hgs)] at hg
      simp at hg
  | succ N ih =>
      intro gs hlen at_ ts hwf
      cases gs with
      | nil => intro g hg; simp at hg
      | cons x gs =>
          have hlenN : gs.length ≤ N := by
            simp only [List.length_cons] at hlen
            omega
          have keep : foldFrom wdraws targets ts at_ (x :: gs)
                = x :: foldFrom wdraws targets (ts.step wdraws x) (at_ + 1) gs →
              ∀ g ∈ foldFrom wdraws targets ts at_ (x :: gs), g.Wf := by
            intro heq g hg
            rw [heq] at hg
            rcases List.mem_cons.1 hg with rfl | hg
            · exact hwf g (by simp)
            · exact ih gs hlenN _ _ (fun y hy => hwf y (by simp [hy])) g hg
          cases hrot : rotAngle x with
          | none => exact keep (foldFrom_cons_none hrot)
          | some p =>
              obtain ⟨θ, q⟩ := p
              by_cases hsel : targets[at_]?.getD true = true
              case neg => exact keep (foldFrom_cons_keep hrot (Or.inl (by simpa using hsel)))
              case pos =>
              cases hm : mergeInto wdraws ts (ts.tagOf q) θ gs with
              | none => exact keep (foldFrom_cons_keep hrot (Or.inr hm))
              | some gs' =>
                  obtain ⟨M, rest, g', φ, q', sign, hgseq, hgs'eq, -, -, -⟩ :=
                    mergeInto_spec wdraws (ts.tagOf q) θ gs gs' ts hm
                  have hlen'' : gs'.length ≤ N := by
                    have := mergeInto_length wdraws (ts.tagOf q) θ gs gs' ts hm
                    omega
                  have hwf' : ∀ y ∈ gs', y.Wf := by
                    intro y hy
                    rw [hgs'eq] at hy
                    rcases List.mem_append.1 hy with hy | hy
                    · exact hwf y (by rw [hgseq]; simp [hy])
                    · rcases List.mem_cons.1 hy with rfl | hy
                      · trivial
                      · exact hwf y (by rw [hgseq]; simp [hy])
                  intro g hg
                  rw [foldFrom_cons_merge hrot hsel hm] at hg
                  exact ih gs' hlen'' (at_ + 1) ts hwf' g hg

theorem phaseFoldGates_wf {k n : Nat} (wdraws : Nat → Tag) {gs : List Gate}
    (h : ∀ g ∈ gs, g.Wf) : ∀ g ∈ phaseFoldGates k wdraws n gs, g.Wf :=
  emitAll_wf (foldFrom_wf (k := k) wdraws _ gs.length gs le_rfl 0 _ h)

/-! ## Operand ranges are preserved

`Wf` above is about *distinctness*; this is about *range*, and together they are what
`Pass.wf_run` and `Pass.wellFormed_run` ask for. The two arguments have the same shape
because the pass only ever invents one kind of gate: a diagonal rotation on a wire the gate
it replaced already used. -/

theorem diagRun_shape (j : Nat) (q : Qubit) :
    ∀ g ∈ diagRun j q, g.qubitsOf = [q] ∧ g.cbitsOf = [] := by
  match j with
  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 =>
      intro g hg; fin_cases hg <;> exact ⟨rfl, rfl⟩
  | (i + 8) => intro g hg; simp [diagRun] at hg

theorem emitRotation_shape (q : Qubit) (a : ℚ) :
    ∀ g ∈ emitRotation q a, g.qubitsOf = [q] ∧ g.cbitsOf = [] := by
  by_cases h0 : BlockState.angleMod a = 0
  · rw [emitRotation_eq_nil h0]; simp
  · cases hcl : classifyQuarterPi (BlockState.angleMod a) with
    | some j => rw [emitRotation_eq_diagRun h0 hcl]; exact diagRun_shape j q
    | none =>
        rw [emitRotation_eq_rz h0 hcl]
        intro g hg
        rw [List.mem_singleton.1 hg]
        exact ⟨rfl, rfl⟩

/-- A gate the folder can re-emit lives on one of its own wires. -/
theorem rotAngle_mem {g : Gate} {a : ℚ} {q : Qubit} (h : rotAngle g = some (a, q)) :
    q ∈ g.qubitsOf := by
  cases g <;> simp_all [rotAngle, Gate.qubitsOf]

theorem emitRotation_inRange {n m : Nat} {g : Gate} {q : Qubit} {a : ℚ}
    (hg : g.InRange n m) (hq : q ∈ g.qubitsOf) : ∀ g' ∈ emitRotation q a, g'.InRange n m := by
  intro g' hg'
  obtain ⟨h₁, h₂⟩ := emitRotation_shape q a g' hg'
  exact Gate.InRange.onWire hg hq h₁ h₂

theorem emitAll_inRange {n m : Nat} {gs : List Gate} (h : ∀ g ∈ gs, g.InRange n m) :
    ∀ g ∈ emitAll gs, g.InRange n m := by
  induction gs with
  | nil => intro g hg; simp [emitAll] at hg
  | cons x gs ih =>
      intro g hg
      rw [emitAll] at hg
      rcases List.mem_append.1 hg with hg | hg
      · cases hrot : rotAngle x with
        | some p =>
            obtain ⟨a, q⟩ := p
            rw [hrot] at hg
            exact emitRotation_inRange (h x (by simp)) (rotAngle_mem hrot) g hg
        | none =>
            rw [hrot] at hg
            rw [List.mem_singleton.1 hg]
            exact h x (by simp)
      · exact ih (fun y hy => h y (by simp [hy])) g hg

theorem foldFrom_inRange {k n m : Nat} (wdraws : Nat → Tag) (targets : Array Bool) :
    ∀ (N : Nat) (gs : List Gate), gs.length ≤ N → ∀ (at_ : Nat) (ts : TState k),
      (∀ g ∈ gs, g.InRange n m) →
        ∀ g ∈ foldFrom wdraws targets ts at_ gs, g.InRange n m := by
  intro N
  induction N with
  | zero =>
      intro gs hgs at_ ts _ g hg
      rw [List.eq_nil_of_length_eq_zero (Nat.le_zero.1 hgs)] at hg
      simp at hg
  | succ N ih =>
      intro gs hlen at_ ts hin
      cases gs with
      | nil => intro g hg; simp at hg
      | cons x gs =>
          have hlenN : gs.length ≤ N := by
            simp only [List.length_cons] at hlen
            omega
          have keep : foldFrom wdraws targets ts at_ (x :: gs)
                = x :: foldFrom wdraws targets (ts.step wdraws x) (at_ + 1) gs →
              ∀ g ∈ foldFrom wdraws targets ts at_ (x :: gs), g.InRange n m := by
            intro heq g hg
            rw [heq] at hg
            rcases List.mem_cons.1 hg with rfl | hg
            · exact hin g (by simp)
            · exact ih gs hlenN _ _ (fun y hy => hin y (by simp [hy])) g hg
          cases hrot : rotAngle x with
          | none => exact keep (foldFrom_cons_none hrot)
          | some p =>
              obtain ⟨θ, q⟩ := p
              by_cases hsel : targets[at_]?.getD true = true
              case neg => exact keep (foldFrom_cons_keep hrot (Or.inl (by simpa using hsel)))
              case pos =>
              cases hm : mergeInto wdraws ts (ts.tagOf q) θ gs with
              | none => exact keep (foldFrom_cons_keep hrot (Or.inr hm))
              | some gs' =>
                  obtain ⟨M, rest, g', φ, q', sign, hgseq, hgs'eq, -, hrot', -⟩ :=
                    mergeInto_spec wdraws (ts.tagOf q) θ gs gs' ts hm
                  have hlen'' : gs'.length ≤ N := by
                    have := mergeInto_length wdraws (ts.tagOf q) θ gs gs' ts hm
                    omega
                  have hin' : ∀ y ∈ gs', y.InRange n m := by
                    intro y hy
                    rw [hgs'eq] at hy
                    rcases List.mem_append.1 hy with hy | hy
                    · exact hin y (by rw [hgseq]; simp [hy])
                    · rcases List.mem_cons.1 hy with rfl | hy
                      · -- the merged rotation sits on the wire of the gate it replaced
                        exact Gate.InRange.onWire
                          (hin g' (by rw [hgseq]; simp)) (rotAngle_mem hrot') rfl rfl
                      · exact hin y (by rw [hgseq]; simp [hy])
                  intro g hg
                  rw [foldFrom_cons_merge hrot hsel hm] at hg
                  exact ih gs' hlen'' (at_ + 1) ts hin' g hg

/-- **Phase folding keeps every operand in range.** -/
theorem phaseFoldGates_inRange {k n' n m : Nat} (wdraws : Nat → Tag) {gs : List Gate}
    (h : ∀ g ∈ gs, g.InRange n m) :
    ∀ g ∈ phaseFoldGates k wdraws n' gs, g.InRange n m :=
  emitAll_inRange (foldFrom_inRange (k := k) wdraws _ gs.length gs le_rfl 0 _ h)

/-! ## The compared parities are bounded -/

theorem fresh_steps_le (st : AState) (gs : List Gate) :
    (st.steps gs).fresh ≤ st.fresh + gs.countP Gate.allocates := by
  induction gs generalizing st with
  | nil => simp
  | cons g gs ih =>
      have := ih (st.step g)
      rw [AState.steps_cons, List.countP_cons]
      by_cases hall : Gate.allocates g = true
      · rw [if_pos hall]
        have hstep : (st.step g).fresh ≤ st.fresh + 1 := by cases g <;> simp [AState.step]
        omega
      · rw [if_neg hall]
        have hstep : (st.step g).fresh = st.fresh := by
          cases g <;> simp_all [AState.step, Gate.allocates]
        omega

theorem bounded_formsOf {n : Nat} {st : AState} (hst : st.Bounded) {m : Nat}
    (h : st.fresh ≤ m) : ∀ p ∈ formsOf n st, Form.Bounded m p := by
  intro p hp
  rcases List.mem_cons.1 hp with rfl | hp
  · exact Form.bounded_const _ false
  · rcases List.mem_map.1 hp with ⟨q, -, rfl⟩
    exact Form.bounded_mono h (hst q)

theorem bounded_visited {n : Nat} : ∀ (gs : List Gate) (st : AState), st.Bounded →
    ∀ {m : Nat}, st.fresh + gs.countP Gate.allocates ≤ m →
      ∀ p ∈ visited n st gs, Form.Bounded m p := by
  intro gs
  induction gs with
  | nil =>
      intro st hst m h p hp
      exact bounded_formsOf hst (by simpa using h) p hp
  | cons g gs ih =>
      intro st hst m h p hp
      rw [List.countP_cons] at h
      by_cases hall : Gate.allocates g = true
      · rw [if_pos hall] at h
        have hstep : (st.step g).fresh ≤ st.fresh + 1 := by cases g <;> simp [AState.step]
        rcases List.mem_append.1 hp with hp | hp
        · exact bounded_formsOf hst (by omega) p hp
        · exact ih (st.step g) (AState.bounded_step hst g) (by omega) p hp
      · rw [if_neg hall] at h
        have hstep : (st.step g).fresh = st.fresh := by
          cases g <;> simp_all [AState.step, Gate.allocates]
        rcases List.mem_append.1 hp with hp | hp
        · exact bounded_formsOf hst (by omega) p hp
        · exact ih (st.step g) (AState.bounded_step hst g) (by omega) p hp

/-- The forms one run of the pass can compare. -/
noncomputable def relevantForms (c : Circuit) : List Form :=
  relevant c.numQubits (AState.initial c.numQubits) c.gates

theorem bounded_relevantForms (c : Circuit) :
    ∀ p ∈ relevantForms c, Form.Bounded (varBound c) p := by
  intro p hp
  have hbase : ∀ r ∈ visited c.numQubits (AState.initial c.numQubits) c.gates,
      Form.Bounded (varBound c) r := by
    refine bounded_visited c.gates (AState.initial c.numQubits) (AState.bounded_initial _) ?_
    simp [varBound, AState.initial]
  rcases List.mem_append.1 hp with hp | hp
  · exact hbase p hp
  · rcases List.mem_map.1 hp with ⟨r, hr, rfl⟩
    exact Form.bounded_flip (hbase r hr)

/-! ## Faithful, unless the tags collide -/

theorem faithful_of_not_collides {m k : Nat} {ps : List Form} {sample : Sample m k}
    (h : ¬ Collides ps sample) : Faithful (liftSample sample) ps := by
  intro p hp q hq hpq
  by_contra hne
  exact h ⟨p, hp, q, hq, hne, hpq⟩

/-! ## A collision-free executable seed

The probabilistic development above remains useful for analysing fixed-width hashing, but the
executable uses this seed.  It gives variable `i` the `i`-th basis vector in `m + 1` bits; the
last bit separates the affine constant.  Consequently hashing is injective on forms bounded by
`m`, with no assumption about a runtime random-number generator. -/

/-- Collision-free draws for affine forms over variables below `m`. -/
def exactDraws (m : Nat) : Draws (m + 1) :=
  fun i j => if j.val = i then 1 else 0

theorem hash_exactDraws_last (m : Nat) (p : Form) (hp : Form.Bounded m p) :
    hash (exactDraws m) p ⟨m, Nat.lt_succ_self m⟩ = p.constant := by
  unfold hash Form.eval
  have hz : p.coefficients.sum
      (fun i coefficient => coefficient * exactDraws m i ⟨m, Nat.lt_succ_self m⟩) = 0 := by
    calc
      _ = p.coefficients.sum (fun _ (_ : F₂) => (0 : F₂)) := Finsupp.sum_congr fun i hi => by
        have hil := hp i hi
        simp [exactDraws]
        omega
      _ = 0 := by simp
  rw [hz, add_zero]

theorem hash_exactDraws_coord {m : Nat} (p : Form) {i : Nat} (hi : i < m) :
    hash (exactDraws m) p ⟨i, Nat.lt_succ_of_lt hi⟩ = p.constant + p.coefficients i := by
  unfold hash Form.eval
  congr 1
  rw [Finsupp.sum_eq_single i]
  · simp [exactDraws]
  · intro b _ hbi
    simp [exactDraws, Ne.symm hbi]
  · simp [exactDraws]

theorem exactDraws_injective {m : Nat} {p q : Form} (hp : Form.Bounded m p)
    (hq : Form.Bounded m q) (h : hash (exactDraws m) p = hash (exactDraws m) q) : p = q := by
  apply Form.ext
  · simpa [hash_exactDraws_last m p hp, hash_exactDraws_last m q hq] using
      congrFun h ⟨m, Nat.lt_succ_self m⟩
  · apply Finsupp.ext
    intro i
    by_cases hi : i < m
    · have hc : p.constant = q.constant := by
        simpa [hash_exactDraws_last m p hp, hash_exactDraws_last m q hq] using
          congrFun h ⟨m, Nat.lt_succ_self m⟩
      have hcoord := congrFun h ⟨i, Nat.lt_succ_of_lt hi⟩
      rw [hash_exactDraws_coord p hi, hash_exactDraws_coord q hi, hc] at hcoord
      exact add_left_cancel hcoord
    · have hp0 : p.coefficients i = 0 := by
        by_contra hn
        exact hi (hp i (Finsupp.mem_support_iff.mpr hn))
      have hq0 : q.coefficients i = 0 := by
        by_contra hn
        exact hi (hq i (Finsupp.mem_support_iff.mpr hn))
      rw [hp0, hq0]

theorem exactDraws_faithful (c : Circuit) :
    Faithful (exactDraws (varBound c)) (relevantForms c) := by
  intro p hp q hq heq
  exact exactDraws_injective (bounded_relevantForms c p hp) (bounded_relevantForms c q hq) heq

@[simp] theorem wordToBits_pow (m i : Nat) :
    wordToBits (k := m + 1) (2 ^ i) = exactDraws m i := by
  funext j
  simp [wordToBits, exactDraws, bit, Nat.testBit_two_pow, eq_comm]

/-- The executable collision-free phase folder. -/
def phaseFoldExact (c : Circuit) : Circuit :=
  phaseFold (varBound c + 1) (fun i => 2 ^ i) c

/-- Phase folding with a collision-free seed.  Unlike `PhaseFoldRand`, this is a deterministic
`Pass`: every output is correct. -/
def PhaseFoldExact : Pass where
  name := "Phase folding"
  run := phaseFoldExact
  numQubits_run _ := rfl
  numCbits_run _ := rfl
  wf_run c hc := phaseFoldGates_wf _ hc
  wellFormed_run c _ hc := phaseFoldGates_inRange _ hc
  flagsOk_run c _ := Circuit.flagsOk_withGates _ _
  correct c hc := phaseFoldGates_correct (wordToBits_pow (varBound c))
    c.gates hc (exactDraws_faithful c)

/-! ## The pass -/

noncomputable section

/-- **Phase folding, as a randomized pass.** The seed is one uniform `k`-bit tag per
variable; the failure probability is the chance that two of the parities this circuit makes
the pass compare hash alike. -/
def PhaseFoldRand (k : Nat) : RandPass where
  name := "Phase folding"
  Seed := fun c => Sample (varBound c) k
  dist := fun _ => PMF.uniformOfFintype _
  run := fun c s => phaseFold k (wordsOf k (liftSample s)) c
  error := fun c => ((relevantForms c).length.choose 2 : ℝ≥0∞) * ((2 : ℝ≥0∞)⁻¹) ^ k
  numQubits_run _ _ := rfl
  numCbits_run _ _ := rfl
  wf_run c s hc := phaseFoldGates_wf (wordsOf k (liftSample s)) hc
  wellFormed_run c s _ hc := phaseFoldGates_inRange (wordsOf k (liftSample s)) hc
  flagsOk_run c _ _ := Circuit.flagsOk_withGates _ _
  correct c hc := by
    refine le_trans ((PMF.uniformOfFintype (Sample (varBound c) k)).toOuterMeasure_mono ?_)
      (collides_probability_le (relevantForms c) (bounded_relevantForms c))
    intro s hs
    by_contra hcol
    exact hs.1 (phaseFoldGates_correct (wordToBits_wordsOf k (liftSample s)) c.gates hc
      (faithful_of_not_collides hcol))

@[simp] theorem PhaseFoldRand_run (k : Nat) (c : Circuit) (s : (PhaseFoldRand k).Seed c) :
    (PhaseFoldRand k).run c s = phaseFold k (wordsOf k (liftSample s)) c := rfl

/-- **What `phaseFoldIO` computes is what the bound is about.** This optional runner draws an element
of `Sample (varBound c) k` — the space `correct` above takes a measure over — and applies the
pass at it; this says the two are the same function, by definition and not by resemblance.
The one thing left unproved is that the draw is uniform, which is a fact about `IO.rand` and
not about any Lean term. -/
theorem phaseFoldIO_run (k : Nat) (c : Circuit) (s : Sample (varBound c) k) :
    phaseFold k (wordsOf k (padSample s)) c = (PhaseFoldRand k).run c s := rfl

/-- The failure bound in closed form: with `t` compared parities the pass is wrong with
probability at most `C(t,2)·2⁻ᵏ`, so doubling the tag width squares the odds against it. -/
theorem PhaseFoldRand_error (k : Nat) (c : Circuit) :
    (PhaseFoldRand k).error c =
      ((relevantForms c).length.choose 2 : ℝ≥0∞) * ((2 : ℝ≥0∞)⁻¹) ^ k := rfl

end

end TzapLean
