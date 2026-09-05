import TzapLean.GateAlgebra

/-!
# Circuit Equivalence

`Equivalent n m gs hs` says two gate lists denote the same channel on `n` qubits and `m`
classical bits: `⟦gs⟧ ρ = ⟦hs⟧ ρ` for every classical-quantum state. This is the obligation
every `Pass` has to discharge.

Two features of the density-matrix semantics matter here:

* **Global phase is free.** `conj (c • U) = conj U` whenever `‖c‖ = 1`, since the phase and
  its conjugate cancel. Rewrites that are only correct up to global phase — the `H·S·H`
  identities of the Hadamard sweep, for instance — are therefore *exactly* equivalent as
  channels, with no "up to phase" caveat anywhere in the development.
* **`measure` and `reset` are ordinary members of the list.** Commuting a unitary gate past
  them is proved once and for all in `step_comm_of_disjoint`, so a pass may move a gate
  across a measurement on another wire, exactly as the Rust `cancel_pairs` does.

The main results are `step_comm_of_disjoint` (a unitary gate commutes with any gate on
disjoint wires), `Equivalent.selfInverse_pair` (a self-inverse gate cancels its twin), and
`Equivalent.move_past` (a unitary gate may be moved to the front of a run of gates it is
disjoint from) — between them, everything `CancelGates`' first sweep does.
-/

namespace TzapLean

open Matrix

noncomputable section

/-- Two gate lists denote the same channel on `n` qubits and `m` classical bits. -/
def Equivalent (n m : Nat) (gs hs : List Gate) : Prop :=
  ∀ ρ : CQState n m, denote (n := n) (m := m) gs ρ = denote hs ρ

namespace Equivalent

@[refl] theorem refl (n m : Nat) (gs : List Gate) : Equivalent n m gs gs := fun _ => rfl

theorem symm {n m : Nat} {gs hs : List Gate} (h : Equivalent n m gs hs) :
    Equivalent n m hs gs := fun ρ => (h ρ).symm

theorem trans {n m : Nat} {gs hs ks : List Gate}
    (h₁ : Equivalent n m gs hs) (h₂ : Equivalent n m hs ks) : Equivalent n m gs ks :=
  fun ρ => (h₁ ρ).trans (h₂ ρ)

/-- Equivalent tails give equivalent lists under a common prefix. -/
theorem append_left {n m : Nat} (pre : List Gate) {gs hs : List Gate}
    (h : Equivalent n m gs hs) : Equivalent n m (pre ++ gs) (pre ++ hs) := by
  intro ρ
  rw [denote_append, denote_append]
  exact h _

/-- Equivalent heads give equivalent lists under a common suffix. -/
theorem append_right {n m : Nat} (suf : List Gate) {gs hs : List Gate}
    (h : Equivalent n m gs hs) : Equivalent n m (gs ++ suf) (hs ++ suf) := by
  intro ρ
  rw [denote_append, denote_append, h ρ]

/-- Rewriting a contiguous window of a circuit. -/
theorem window {n m : Nat} (pre suf : List Gate) {gs hs : List Gate}
    (h : Equivalent n m gs hs) : Equivalent n m (pre ++ gs ++ suf) (pre ++ hs ++ suf) := by
  simpa [List.append_assoc] using append_left pre (append_right suf h)

end Equivalent

namespace Circuit.Checked

/-- Two checked circuits with the same register indices denote the same channel. Cached
metadata is intentionally excluded because it does not contribute to circuit semantics. -/
abbrev Equivalent (c d : Circuit.Checked n m) : Prop :=
  TzapLean.Equivalent n m c.raw.gates d.raw.gates

end Circuit.Checked

/-! ## Global phase is invisible -/

theorem conj_smul {n : Nat} (c : ℂ) (U ρ : Density n) (hc : c * star c = 1) :
    conj (c • U) ρ = conj U ρ := by
  simp only [conj, Matrix.conjTranspose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  rw [mul_comm] at hc
  rw [hc, one_smul]

theorem conj_add {n : Nat} (U ρ σ : Density n) : conj U (ρ + σ) = conj U ρ + conj U σ := by
  simp [conj, Matrix.mul_add, Matrix.add_mul]

theorem conj_sum {n : Nat} {ι : Type*} (s : Finset ι) (U : Density n) (f : ι → Density n) :
    conj U (∑ i ∈ s, f i) = ∑ i ∈ s, conj U (f i) := by
  simp [conj, Finset.mul_sum, Finset.sum_mul]

/-! ## Commutation -/

/-- Matrices with disjoint supports commute after conjugate transposition too. -/
theorem conjTranspose_comm_of_mul_comm {n : Nat} {U V : Density n} (h : U * V = V * U) :
    Uᴴ * Vᴴ = Vᴴ * Uᴴ := by
  rw [← Matrix.conjTranspose_mul, ← Matrix.conjTranspose_mul, h]

theorem conj_comm {n : Nat} {U V : Density n} (h : U * V = V * U) (ρ : Density n) :
    conj U (conj V ρ) = conj V (conj U ρ) := by
  rw [← conj_mul, ← conj_mul, h]

/-! ### Equation lemmas for the non-unitary steps -/

theorem step_measure_of_lt {n m : Nat} (q : Qubit) (c : CBit) (hc : c < m) (ρ : CQState n m) :
    step (.measure q c) ρ =
      fun w => conj (proj n q (w.read c)) (ρ (w.write c false) + ρ (w.write c true)) := by
  simp only [step, if_pos hc]

theorem step_measure_of_ge {n m : Nat} (q : Qubit) (c : CBit) (hc : ¬ c < m) (ρ : CQState n m) :
    step (.measure q c) ρ = fun w => ∑ b : Bool, conj (proj n q b) (ρ w) := by
  simp only [step, if_neg hc]

theorem step_reset {n m : Nat} (q : Qubit) (ρ : CQState n m) :
    step (.reset q) ρ =
      fun w => ∑ b : Bool, resetKraus n q b * ρ w * (resetKraus n q b)ᴴ := rfl

/-- Two unitary gates whose matrices commute commute as channel steps. -/
theorem step_comm_unitary {n m : Nat} {g₁ g₂ : Gate}
    (h₁ : g₁.isUnitary = true) (h₂ : g₂.isUnitary = true)
    (hcomm : gateUnitary n g₁ * gateUnitary n g₂ = gateUnitary n g₂ * gateUnitary n g₁)
    (ρ : CQState n m) : step g₂ (step g₁ ρ) = step g₁ (step g₂ ρ) := by
  have e1 : ∀ σ : CQState n m, step g₁ σ = fun w => conj (gateUnitary n g₁) (σ w) :=
    fun σ => step_of_isUnitary h₁ σ
  have e2 : ∀ σ : CQState n m, step g₂ σ = fun w => conj (gateUnitary n g₂) (σ w) :=
    fun σ => step_of_isUnitary h₂ σ
  simp only [e1, e2]
  funext w
  exact conj_comm hcomm.symm _

/-- **A unitary gate commutes with any gate on disjoint wires** — including `measure` and
`reset`. This is the single fact behind every "commute past unrelated gates" step in the
Rust pass. -/
theorem step_comm_of_disjoint {n m : Nat} {g₁ g₂ : Gate} (h₁ : g₁.isUnitary = true)
    (hdisj : Wires.Disjoint g₁.support g₂.support) (ρ : CQState n m) :
    step g₂ (step g₁ ρ) = step g₁ (step g₂ ρ) := by
  have hU : SupportedOn (gateUnitary n g₁) g₁.support := gateUnitary_supportedOn n g₁
  have e1 : ∀ σ : CQState n m, step g₁ σ = fun w => conj (gateUnitary n g₁) (σ w) :=
    fun σ => step_of_isUnitary h₁ σ
  cases g₂ with
  | measure q c =>
      have hproj : ∀ b : Bool,
          gateUnitary n g₁ * proj n q b = proj n q b * gateUnitary n g₁ := by
        intro b
        refine hU.mul_comm (proj_supportedOn n q b) fun r hr => ?_
        simpa [Gate.support] using hdisj r hr
      by_cases hc : c < m
      · rw [step_measure_of_lt _ _ hc, step_measure_of_lt _ _ hc, e1, e1]
        funext w
        simp only
        rw [← conj_add, ← conj_mul, ← conj_mul, hproj]
      · rw [step_measure_of_ge _ _ hc, step_measure_of_ge _ _ hc, e1, e1]
        funext w
        simp only
        rw [conj_sum]
        exact Finset.sum_congr rfl fun b _ => by rw [← conj_mul, ← conj_mul, hproj]
  | reset q =>
      have hk : ∀ b : Bool,
          gateUnitary n g₁ * resetKraus n q b = resetKraus n q b * gateUnitary n g₁ := by
        intro b
        refine hU.mul_comm (resetKraus_supportedOn n q b) fun r hr => ?_
        simpa [Gate.support] using hdisj r hr
      rw [step_reset, step_reset, e1, e1]
      funext w
      simp only
      rw [conj_sum]
      refine Finset.sum_congr rfl fun b _ => ?_
      have hk' : (resetKraus n q b)ᴴ * (gateUnitary n g₁)ᴴ =
          (gateUnitary n g₁)ᴴ * (resetKraus n q b)ᴴ := by
        rw [← Matrix.conjTranspose_mul, ← Matrix.conjTranspose_mul, hk b]
      simp only [conj]
      calc resetKraus n q b * (gateUnitary n g₁ * ρ w * (gateUnitary n g₁)ᴴ) *
              (resetKraus n q b)ᴴ
          = (resetKraus n q b * gateUnitary n g₁) * ρ w *
              ((gateUnitary n g₁)ᴴ * (resetKraus n q b)ᴴ) := by
            simp [Matrix.mul_assoc]
        _ = (gateUnitary n g₁ * resetKraus n q b) * ρ w *
              ((resetKraus n q b)ᴴ * (gateUnitary n g₁)ᴴ) := by rw [← hk b, ← hk']
        _ = gateUnitary n g₁ * (resetKraus n q b * ρ w * (resetKraus n q b)ᴴ) *
              (gateUnitary n g₁)ᴴ := by simp [Matrix.mul_assoc]
  | _ =>
      exact step_comm_unitary h₁ rfl
        (hU.mul_comm (gateUnitary_supportedOn n _) hdisj) ρ

/-! ## Equivalences a pass can use -/

/-- **Matrix identities lift to circuit equivalences, up to global phase.** For
measurement-free fragments this reduces a rewrite's correctness to a `2ⁿ × 2ⁿ` identity that
may differ by any unit scalar. -/
theorem equivalent_of_unitary_smul {n m : Nat} {gs hs : List Gate}
    (hgs : ∀ g ∈ gs, g.isUnitary = true) (hhs : ∀ g ∈ hs, g.isUnitary = true)
    (c : ℂ) (hc : c * star c = 1) (h : unitary n gs = c • unitary n hs) :
    Equivalent n m gs hs := by
  intro ρ
  rw [denote_eq_conj_unitary gs hgs, denote_eq_conj_unitary hs hhs, h]
  funext w
  exact conj_smul c _ _ hc

/-- The same with the two matrices literally equal. -/
theorem equivalent_of_unitary_eq {n m : Nat} {gs hs : List Gate}
    (hgs : ∀ g ∈ gs, g.isUnitary = true) (hhs : ∀ g ∈ hs, g.isUnitary = true)
    (h : unitary n gs = unitary n hs) : Equivalent n m gs hs :=
  equivalent_of_unitary_smul hgs hhs 1 (by simp) (by simpa using h)

/-- A self-inverse gate cancels its own repetition — the rewrite `CancelGates`' first sweep
is built on. -/
theorem Equivalent.selfInverse_pair {n m : Nat} (g : Gate) (hwf : g.Wf)
    (hs : g.isSelfInverse = true) : Equivalent n m [g, g] [] := by
  have hu : g.isUnitary = true := by
    cases g <;> simp_all [Gate.isSelfInverse, Gate.isUnitary, Gate.isMeasurement]
  refine equivalent_of_unitary_eq (fun g' hg' => by simp at hg'; subst hg'; exact hu)
    (by simp) ?_
  simp only [unitary_cons, unitary_nil, Matrix.one_mul]
  exact gateUnitary_mul_self n g hwf hs

/-- Two gates on disjoint wires, one of them unitary, may be swapped. -/
theorem Equivalent.swap_of_disjoint {n m : Nat} {g₁ g₂ : Gate} (h₁ : g₁.isUnitary = true)
    (hdisj : Wires.Disjoint g₁.support g₂.support) : Equivalent n m [g₁, g₂] [g₂, g₁] := by
  intro ρ
  simp only [denote_cons, denote_nil]
  exact step_comm_of_disjoint h₁ hdisj ρ

/-- **A unitary gate may be moved to the front of a run of gates it does not touch.** This
is the "commute past unrelated gates" step of `cancel_pairs`, in the form the pass needs:
the gate reaches back over an arbitrary window of disjoint gates. -/
theorem Equivalent.move_past {n m : Nat} {g : Gate} (hu : g.isUnitary = true)
    (gs : List Gate) (hdisj : ∀ g' ∈ gs, Wires.Disjoint g.support g'.support) :
    Equivalent n m (gs ++ [g]) (g :: gs) := by
  induction gs with
  | nil => exact Equivalent.refl _ _ _
  | cons g' gs ih =>
      have hg' : Wires.Disjoint g.support g'.support := hdisj g' (by simp)
      have htail : ∀ g'' ∈ gs, Wires.Disjoint g.support g''.support :=
        fun g'' hg'' => hdisj g'' (by simp [hg''])
      intro ρ
      calc denote ((g' :: gs) ++ [g]) ρ
          = denote (gs ++ [g]) (step g' ρ) := rfl
        _ = denote (g :: gs) (step g' ρ) := ih htail (step g' ρ)
        _ = denote gs (step g (step g' ρ)) := rfl
        _ = denote gs (step g' (step g ρ)) := by rw [step_comm_of_disjoint hu hg' ρ]
        _ = denote (g :: g' :: gs) ρ := rfl

/-- Dually: a unitary gate may be moved to the back of a run of gates it does not touch. -/
theorem Equivalent.move_back {n m : Nat} {g : Gate} (hu : g.isUnitary = true)
    (gs : List Gate) (hdisj : ∀ g' ∈ gs, Wires.Disjoint g.support g'.support) :
    Equivalent n m (g :: gs) (gs ++ [g]) :=
  (Equivalent.move_past hu gs hdisj).symm

/-- **Disjoint runs commute wholesale.** A block of unitary gates may be moved across a
block of gates it shares no wire with. -/
theorem Equivalent.append_comm {n m : Nat} (A B : List Gate) (hu : ∀ g ∈ A, g.isUnitary = true)
    (hdisj : ∀ g ∈ A, ∀ g' ∈ B, Wires.Disjoint g.support g'.support) :
    Equivalent n m (A ++ B) (B ++ A) := by
  induction A with
  | nil => simp [Equivalent.refl]
  | cons a A ih =>
      have hu' : ∀ g ∈ A, g.isUnitary = true := fun g hg => hu g (by simp [hg])
      have hd' : ∀ g ∈ A, ∀ g' ∈ B, Wires.Disjoint g.support g'.support :=
        fun g hg => hdisj g (by simp [hg])
      have step₁ : Equivalent n m (a :: (A ++ B)) (a :: (B ++ A)) :=
        Equivalent.append_left [a] (ih hu' hd')
      have step₂ : Equivalent n m (a :: (B ++ A)) ((B ++ [a]) ++ A) := by
        have h := Equivalent.append_right (n := n) (m := m) A
          (Equivalent.move_back (n := n) (m := m) (hu a (by simp)) B (hdisj a (by simp)))
        simpa using h
      have : Equivalent n m (a :: (A ++ B)) (B ++ (a :: A)) := by
        refine Equivalent.trans step₁ (Equivalent.trans step₂ ?_)
        simp [Equivalent.refl]
      simpa using this

end
end TzapLean
