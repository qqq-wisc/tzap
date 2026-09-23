import TzapLean.PhaseFoldProof
import TzapLean.PhaseFoldNonlinear
import TzapLean.ExecutableRandPass

/-!
# Affine phase-folding reference and nonlinear sampled runner

The affine `phaseFoldGates_correct` says folding is right whenever its tags are faithful,
and `collides_probability_le` says unfaithful tags are unlikely. Together they give the
reference bound

```
Pr_{s ← uniform} [ ⟦phaseFold s c⟧ ≠ ⟦c⟧ ]  ≤  C(L, 2) · 2⁻ᵏ
```

where `L` is the number of parities (and complements) this circuit makes the pass compare.
The seed is one ideal uniform `k`-bit tag per affine variable — `Sample (varBound c) k`.
The CLI uses the nonlinear OS-backed runner with 128-bit samples.

The affine proof covers the conservative CCX fallback. `GF128Bridge.lean` proves the
nonlinear pass's 128-bit sampled probability bound and supplies the sole verified
`PhaseFoldRand` wired into the optimizer.
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
the optimizer's checked representation and output boundary require. The two arguments have the same shape
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
noncomputable def relevantForms (c : RawCircuit) : List Form :=
  relevant c.numQubits (AState.initial c.numQubits) c.gates

theorem bounded_relevantForms (c : RawCircuit) :
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

/-! ## Executable randomized phase folding -/

/-- The concrete nonlinear transformation used by the executable.  Its random words are
normalized to 128 bits by `Fingerprint.fresh`, and the CLI instantiates this with `k = 128`. -/
def phaseFoldNonlinearWithSample (k : Nat) (c : Circuit n m)
    (s : Sample (varBound c.raw) k) : Circuit n m :=
  ⟨phaseFoldNonlinear (wordsOf k (liftSample s)) c.raw,
    (phaseFoldNonlinear_numQubits _ c.raw).trans c.numQubits_eq,
    (phaseFoldNonlinear_numCbits _ c.raw).trans c.numCbits_eq,
    phaseFoldGatesNonlinear_wf _ c.wf⟩

/-- The runtime nonlinear phase-folding pass. Every invocation obtains a fresh 128-bit sample
from `IO.getRandomBytes` and calls the same pure transformation as the verified pass. -/
def PhaseFoldRandExec : ExecutableRandPass where
  name := "Phase folding"
  run := fun c => do
    let s ← randomSample (varBound c.raw) 128
    return phaseFoldNonlinearWithSample 128 c s

end TzapLean
