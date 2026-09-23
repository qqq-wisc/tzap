import TzapLean.NonlinearBranches

/-!
# Channel lemmas for nonlinear phase folding

The path argument is branchwise, while circuit semantics is a linear channel on
classical-quantum states. These lemmas supply the linearity needed to sum Kraus branches
and show that the nonlinear phase hop survives conjugation on each branch.
-/

namespace TzapLean

noncomputable section

/-- Every gate acts additively on classical-quantum states, including measurement with
classical-memory overwrite and reset. -/
theorem nonlinear_step_add {n m : Nat} (gate : Gate) (ρ σ : CQState n m) :
    step gate (ρ + σ) = step gate ρ + step gate σ := by
  funext w
  cases gate with
  | measure q c =>
      by_cases hc : c < m
      · simp [step, hc, conj_add, add_left_comm, add_comm]
      · simp [step, hc, conj_add, Finset.sum_add_distrib]
  | reset q =>
      simp only [step, Matrix.mul_add, Matrix.add_mul, Finset.sum_add_distrib,
        Pi.add_apply]
  | h q => simp [step, conj_add]
  | x q => simp [step, conj_add]
  | z q => simp [step, conj_add]
  | s q => simp [step, conj_add]
  | sdg q => simp [step, conj_add]
  | t q => simp [step, conj_add]
  | tdg q => simp [step, conj_add]
  | rz angle q => simp [step, conj_add]
  | cnot c t => simp [step, conj_add]
  | cz c t => simp [step, conj_add]
  | ccx c₁ c₂ t => simp [step, conj_add]
  | ccz c₁ c₂ t => simp [step, conj_add]

/-- The whole circuit denotes an additive channel. -/
theorem nonlinear_denote_add {n m : Nat} (gates : List Gate) (ρ σ : CQState n m) :
    denote gates (ρ + σ) = denote gates ρ + denote gates σ := by
  induction gates generalizing ρ σ with
  | nil => rfl
  | cons gate gates ih =>
      rw [denote_cons, nonlinear_step_add, ih]
      rfl

theorem nonlinear_step_zero {n m : Nat} (gate : Gate) :
    step gate (0 : CQState n m) = 0 := by
  funext w
  cases gate with
  | measure q c =>
      by_cases hc : c < m <;> simp [step, hc, conj]
  | reset q => simp [step]
  | h q => simp [step, conj]
  | x q => simp [step, conj]
  | z q => simp [step, conj]
  | s q => simp [step, conj]
  | sdg q => simp [step, conj]
  | t q => simp [step, conj]
  | tdg q => simp [step, conj]
  | rz angle q => simp [step, conj]
  | cnot c t => simp [step, conj]
  | cz c t => simp [step, conj]
  | ccx c₁ c₂ t => simp [step, conj]
  | ccz c₁ c₂ t => simp [step, conj]

theorem nonlinear_denote_zero {n m : Nat} (gates : List Gate) :
    denote gates (0 : CQState n m) = 0 := by
  induction gates with
  | nil => rfl
  | cons gate gates ih =>
      rw [denote_cons, nonlinear_step_zero, ih]

/-- Circuit semantics distributes over finite sums of Kraus contributions. -/
theorem nonlinear_denote_sum {n m : Nat} {ι : Type*} (s : Finset ι)
    (gates : List Gate) (f : ι → CQState n m) :
    denote gates (∑ i ∈ s, f i) = ∑ i ∈ s, denote gates (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [nonlinear_denote_zero]
  | @insert i s hi ih =>
      simp [hi, nonlinear_denote_add, ih]

theorem nonlinear_step_sum {n m : Nat} {ι : Type*} (s : Finset ι)
    (gate : Gate) (f : ι → CQState n m) :
    step gate (∑ i ∈ s, f i) = ∑ i ∈ s, step gate (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [nonlinear_step_zero]
  | @insert i s hi ih => simp [hi, nonlinear_step_add, ih]

theorem nonlinear_step_list_sum {n m : Nat} (gate : Gate)
    (xs : List (CQState n m)) :
    step gate xs.sum = (xs.map (step gate)).sum := by
  induction xs with
  | nil => simp [nonlinear_step_zero]
  | cons x xs ih => simp [nonlinear_step_add, ih]

theorem nonlinear_denote_list_sum {n m : Nat} (gates : List Gate)
    (xs : List (CQState n m)) :
    denote gates xs.sum = (xs.map (denote gates)).sum := by
  induction xs with
  | nil => simp [nonlinear_denote_zero]
  | cons x xs ih => simp [nonlinear_denote_add, ih]

/-- A density matrix supported at exactly one classical-memory value. -/
def nonlinearSingleMemory {n m : Nat} (w₀ : Memory m) (D : Density n) : CQState n m :=
  fun w => if w = w₀ then D else 0

theorem nonlinear_sum_singleMemory {n m : Nat} (ρ : CQState n m) :
    (∑ w : Memory m, nonlinearSingleMemory w (ρ w)) = ρ := by
  funext w
  simp [nonlinearSingleMemory]

/-- One selected Kraus outcome, with the classical-memory update of measurement.
For a unitary gate only the `false` choice is active. -/
def nonlinearBranchStep {n m : Nat} (gate : Gate) (choice : Bool)
    (ρ : CQState n m) : CQState n m :=
  match gate with
  | .measure q c =>
      if c < m then
        fun w => if w.read c = choice then
          conj (proj n q choice) (ρ (w.write c false) + ρ (w.write c true)) else 0
      else fun w => conj (proj n q choice) (ρ w)
  | .reset q => fun w => conj (resetKraus n q choice) (ρ w)
  | gate => if choice then 0 else fun w => conj (gateUnitary n gate) (ρ w)

def nonlinearBranchMemoryOutput {m : Nat} (gate : Gate) (choice : Bool)
    (w : Memory m) : Memory m :=
  match gate with
  | .measure _ c => w.write c choice
  | _ => w

/-- The unused `true` branch of a unitary gate contributes the zero matrix. -/
def nonlinearEffectiveBranchOp (n : Nat) (gate : Gate) (choice : Bool) : Density n :=
  if gate.isUnitary && choice then 0 else nonlinearBranchOp n gate choice

theorem nonlinear_memory_write_preimage {m : Nat} (w w₀ : Memory m)
    (c : CBit) (choice : Bool) (hc : c < m) (hread : w.read c = choice) :
    w.write c (w₀.read c) = w₀ ↔ w = w₀.write c choice := by
  constructor
  · intro h
    calc w = w.write c (w.read c) := (Memory.write_read_self w c hc).symm
      _ = w.write c choice := by rw [hread]
      _ = (w.write c (w₀.read c)).write c choice := by rw [Memory.write_write]
      _ = w₀.write c choice := by rw [h]
  · intro h
    rw [h, Memory.write_write, Memory.write_read_self w₀ c hc]

theorem nonlinear_memory_write_other_ne {m : Nat} (w w₀ : Memory m)
    (c : CBit) (hc : c < m) : w.write c (!w₀.read c) ≠ w₀ := by
  intro h
  have hr := congrArg (fun v : Memory m => v.read c) h
  rw [Memory.read_write_same _ _ _ hc] at hr
  cases hb : w₀.read c <;> simp [hb] at hr

/-- A selected gate branch sends a single classical-memory block to a single block.
The matrix on that block is the branch's effective Kraus operator. -/
theorem nonlinearBranchStep_singleMemory {n m : Nat} (gate : Gate) (choice : Bool)
    (w₀ : Memory m) (D : Density n) :
    nonlinearBranchStep gate choice (nonlinearSingleMemory w₀ D) =
      nonlinearSingleMemory (nonlinearBranchMemoryOutput gate choice w₀)
        (conj (nonlinearEffectiveBranchOp n gate choice) D) := by
  funext w
  cases gate with
  | measure q c =>
      by_cases hc : c < m
      · have hgoal :
            (if w.read c = choice then
              conj (proj n q choice)
                ((if w.write c false = w₀ then D else 0) +
                  if w.write c true = w₀ then D else 0)
              else 0) =
            if w = w₀.write c choice then conj (proj n q choice) D else 0 := by
          by_cases hr : w.read c = choice
          · have hpre := nonlinear_memory_write_preimage w w₀ c choice hc hr
            have hother := nonlinear_memory_write_other_ne w w₀ c hc
            cases hb : w₀.read c with
            | false =>
                simp only [hb, Bool.not_false] at hpre hother
                simp [hr, hpre, hother, conj]
            | true =>
                simp only [hb, Bool.not_true] at hpre hother
                simp [hr, hpre, hother, conj]
          · have hne : w ≠ w₀.write c choice := by
              intro heq
              apply hr
              rw [heq, Memory.read_write_same _ _ _ hc]
            simp [hr, hne]
        simpa [nonlinearBranchStep, nonlinearSingleMemory,
          nonlinearBranchMemoryOutput, nonlinearEffectiveBranchOp,
          nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, hc] using hgoal
      · have hout : w₀.write c choice = w₀ := by
          funext i
          have hi : (i : Nat) ≠ c := by
            intro heq
            exact hc (heq ▸ i.isLt)
          simp [Memory.write, hi]
        simp [nonlinearBranchStep, nonlinearSingleMemory,
          nonlinearBranchMemoryOutput, nonlinearEffectiveBranchOp,
          nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, hc, hout, conj]
  | reset q =>
      simp [nonlinearBranchStep, nonlinearSingleMemory,
        nonlinearBranchMemoryOutput, nonlinearEffectiveBranchOp,
        nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | _ =>
      cases choice <;> simp [nonlinearBranchStep, nonlinearSingleMemory,
        nonlinearBranchMemoryOutput, nonlinearEffectiveBranchOp,
        nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]

theorem nonlinearBranchStep_add {n m : Nat} (gate : Gate) (choice : Bool)
    (ρ σ : CQState n m) :
    nonlinearBranchStep gate choice (ρ + σ) =
      nonlinearBranchStep gate choice ρ + nonlinearBranchStep gate choice σ := by
  funext w
  cases gate with
  | measure q c =>
      by_cases hc : c < m
      · by_cases hw : w.read c = choice <;>
          simp [nonlinearBranchStep, hc, hw, conj_add, add_left_comm, add_comm]
      · simp [nonlinearBranchStep, hc, conj_add]
  | reset q => simp [nonlinearBranchStep, conj_add]
  | _ => cases choice <;> simp [nonlinearBranchStep, conj_add]

theorem nonlinearBranchStep_zero {n m : Nat} (gate : Gate) (choice : Bool) :
    nonlinearBranchStep gate choice (0 : CQState n m) = 0 := by
  funext w
  cases gate with
  | measure q c =>
      by_cases hc : c < m <;> simp [nonlinearBranchStep, hc, conj]
  | reset q => simp [nonlinearBranchStep, conj]
  | _ => cases choice <;> simp [nonlinearBranchStep, conj]

theorem nonlinearBranchStep_sum {n m : Nat} {ι : Type*} (s : Finset ι)
    (gate : Gate) (choice : Bool) (f : ι → CQState n m) :
    nonlinearBranchStep gate choice (∑ i ∈ s, f i) =
      ∑ i ∈ s, nonlinearBranchStep gate choice (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [nonlinearBranchStep_zero]
  | @insert i s hi ih => simp [hi, nonlinearBranchStep_add, ih]

/-- The two branch outcomes sum to the exact single-gate channel. -/
theorem nonlinear_step_eq_branches {n m : Nat} (gate : Gate) (ρ : CQState n m) :
    step gate ρ = nonlinearBranchStep gate false ρ +
      nonlinearBranchStep gate true ρ := by
  funext w
  cases gate with
  | measure q c =>
      by_cases hc : c < m
      · cases hw : w.read c <;> simp [step, nonlinearBranchStep, hc, hw]
      · simp [step, nonlinearBranchStep, hc, add_comm]
  | reset q =>
      simp [step, nonlinearBranchStep, conj, add_comm]
  | h q => simp [step, nonlinearBranchStep]
  | x q => simp [step, nonlinearBranchStep]
  | z q => simp [step, nonlinearBranchStep]
  | s q => simp [step, nonlinearBranchStep]
  | sdg q => simp [step, nonlinearBranchStep]
  | t q => simp [step, nonlinearBranchStep]
  | tdg q => simp [step, nonlinearBranchStep]
  | rz angle q => simp [step, nonlinearBranchStep]
  | cnot c t => simp [step, nonlinearBranchStep]
  | cz c t => simp [step, nonlinearBranchStep]
  | ccx c₁ c₂ t => simp [step, nonlinearBranchStep]
  | ccz c₁ c₂ t => simp [step, nonlinearBranchStep]

/-- Enumerate all selected-Kraus branches of a circuit, retaining zero branches so
the recursion has a fixed binary shape. -/
def nonlinearBranchFamily {n m : Nat} : List Gate → CQState n m → List (CQState n m)
  | [], ρ => [ρ]
  | gate :: gates, ρ =>
      nonlinearBranchFamily gates (nonlinearBranchStep gate false ρ) ++
        nonlinearBranchFamily gates (nonlinearBranchStep gate true ρ)

/-- Summing all branch cq-states recovers exact circuit semantics. -/
theorem nonlinear_denote_eq_branchFamily_sum {n m : Nat} (gates : List Gate)
    (ρ : CQState n m) : denote gates ρ = (nonlinearBranchFamily gates ρ).sum := by
  induction gates generalizing ρ with
  | nil => simp [nonlinearBranchFamily]
  | cons gate gates ih =>
      rw [denote_cons, nonlinear_step_eq_branches, nonlinear_denote_add,
        ih, ih]
      simp [nonlinearBranchFamily]

/-- Run one fixed sequence of Kraus choices. -/
def nonlinearBranchDenote {n m : Nat} :
    List Gate → (Nat → Bool) → CQState n m → CQState n m
  | [], _, ρ => ρ
  | gate :: gates, choices, ρ =>
      nonlinearBranchDenote gates (fun i => choices (i + 1))
        (nonlinearBranchStep gate (choices 0) ρ)

def nonlinearEffectiveBranchMatrix (n : Nat) : List Gate → (Nat → Bool) → Density n
  | [], _ => 1
  | gate :: gates, choices =>
      nonlinearEffectiveBranchMatrix n gates (fun i => choices (i + 1)) *
        nonlinearEffectiveBranchOp n gate (choices 0)

theorem nonlinearEffectiveBranchMatrix_eq_zero_or_branch (n : Nat)
    (gates : List Gate) (choices : Nat → Bool) :
    nonlinearEffectiveBranchMatrix n gates choices = 0 ∨
      nonlinearEffectiveBranchMatrix n gates choices =
        nonlinearBranchMatrix n gates choices := by
  induction gates generalizing choices with
  | nil => exact Or.inr rfl
  | cons gate gates ih =>
      rcases ih (fun i => choices (i + 1)) with hzero | heq
      · exact Or.inl (by simp [nonlinearEffectiveBranchMatrix, hzero])
      · by_cases hbad : (gate.isUnitary && choices 0) = true
        · have hop : nonlinearEffectiveBranchOp n gate (choices 0) = 0 := by
            simp only [nonlinearEffectiveBranchOp, hbad, ite_true]
          exact Or.inl (by rw [nonlinearEffectiveBranchMatrix, hop, mul_zero])
        · have hop : nonlinearEffectiveBranchOp n gate (choices 0) =
              nonlinearBranchOp n gate (choices 0) := by
            simp only [nonlinearEffectiveBranchOp, if_neg hbad]
          exact Or.inr (by rw [nonlinearEffectiveBranchMatrix, nonlinearBranchMatrix,
            hop, heq])

def nonlinearBranchMemory {m : Nat} :
    List Gate → (Nat → Bool) → Memory m → Memory m
  | [], _, w => w
  | gate :: gates, choices, w =>
      nonlinearBranchMemory gates (fun i => choices (i + 1))
        (nonlinearBranchMemoryOutput gate (choices 0) w)

/-- A fixed branch of any circuit, including measurements that write classical bits,
maps one input memory block to one output block by its composed Kraus matrix. -/
theorem nonlinearBranchDenote_singleMemory {n m : Nat} (gates : List Gate)
    (choices : Nat → Bool) (w₀ : Memory m) (D : Density n) :
    nonlinearBranchDenote gates choices (nonlinearSingleMemory w₀ D) =
      nonlinearSingleMemory (nonlinearBranchMemory gates choices w₀)
        (conj (nonlinearEffectiveBranchMatrix n gates choices) D) := by
  induction gates generalizing choices w₀ D with
  | nil => simp [nonlinearBranchDenote, nonlinearBranchMemory,
      nonlinearEffectiveBranchMatrix, conj_one]
  | cons gate gates ih =>
      rw [nonlinearBranchDenote,
        nonlinearBranchStep_singleMemory gate (choices 0) w₀ D,
        ih]
      simp only [nonlinearBranchMemory, nonlinearEffectiveBranchMatrix]
      rw [conj_mul]

theorem nonlinearBranchDenote_add {n m : Nat} (gates : List Gate)
    (choices : Nat → Bool) (ρ σ : CQState n m) :
    nonlinearBranchDenote gates choices (ρ + σ) =
      nonlinearBranchDenote gates choices ρ + nonlinearBranchDenote gates choices σ := by
  induction gates generalizing choices ρ σ with
  | nil => rfl
  | cons gate gates ih =>
      rw [nonlinearBranchDenote, nonlinearBranchStep_add, ih]
      rfl

theorem nonlinearBranchDenote_zero {n m : Nat} (gates : List Gate)
    (choices : Nat → Bool) :
    nonlinearBranchDenote gates choices (0 : CQState n m) = 0 := by
  induction gates generalizing choices with
  | nil => rfl
  | cons gate gates ih =>
      rw [nonlinearBranchDenote, nonlinearBranchStep_zero, ih]

theorem nonlinearBranchDenote_sum {n m : Nat} {ι : Type*} (s : Finset ι)
    (gates : List Gate) (choices : Nat → Bool) (f : ι → CQState n m) :
    nonlinearBranchDenote gates choices (∑ i ∈ s, f i) =
      ∑ i ∈ s, nonlinearBranchDenote gates choices (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [nonlinearBranchDenote_zero]
  | @insert i s hi ih => simp [hi, nonlinearBranchDenote_add, ih]

theorem nonlinearBranchDenote_list_sum {n m : Nat} (gates : List Gate)
    (choices : Nat → Bool) (xs : List (CQState n m)) :
    nonlinearBranchDenote gates choices xs.sum =
      (xs.map (nonlinearBranchDenote gates choices)).sum := by
  induction xs with
  | nil => simp [nonlinearBranchDenote_zero]
  | cons x xs ih => simp [nonlinearBranchDenote_add, ih]

/-- Exact selected-branch semantics for arbitrary classical memory, including writes. -/
theorem nonlinearBranchDenote_memory_sum {n m : Nat} (gates : List Gate)
    (choices : Nat → Bool) (ρ : CQState n m) :
    nonlinearBranchDenote gates choices ρ =
      ∑ w₀ : Memory m,
        nonlinearSingleMemory (nonlinearBranchMemory gates choices w₀)
          (conj (nonlinearEffectiveBranchMatrix n gates choices) (ρ w₀)) := by
  conv_lhs => rw [← nonlinear_sum_singleMemory ρ]
  rw [nonlinearBranchDenote_sum]
  simp_rw [nonlinearBranchDenote_singleMemory]

def nonlinearBranchChoices : List Gate → List (Nat → Bool)
  | [] => [fun _ => false]
  | _ :: gates =>
      (nonlinearBranchChoices gates).map (fun tail =>
        fun | 0 => false | i + 1 => tail i) ++
      (nonlinearBranchChoices gates).map (fun tail =>
        fun | 0 => true | i + 1 => tail i)

/-- The exact channel is the sum of all fixed-choice branch channels, including
classical-memory writes. -/
theorem nonlinear_denote_eq_branchChoices_sum {n m : Nat} (gates : List Gate)
    (ρ : CQState n m) :
    denote gates ρ =
      ((nonlinearBranchChoices gates).map
        (fun choices => nonlinearBranchDenote gates choices ρ)).sum := by
  induction gates generalizing ρ with
  | nil => simp [nonlinearBranchChoices, nonlinearBranchDenote]
  | cons gate gates ih =>
      rw [denote_cons, nonlinear_step_eq_branches, nonlinear_denote_add,
        ih, ih]
      simp only [nonlinearBranchChoices, List.map_append, List.sum_append, List.map_map]
      rfl

/-- Without classical writes, a selected branch acts independently on every memory
block by its branch matrix. This includes reset and measurements whose outcomes are
discarded. -/
theorem nonlinearBranchDenote_matrix {n m : Nat} (gates : List Gate)
    (choices : Nat → Bool) (ρ : CQState n m)
    (hdiscard : ∀ gate ∈ gates, ∀ q c, gate = Gate.measure q c → m ≤ c)
    (hchoices : ∀ i (hi : i < gates.length),
      (gates[i]).isUnitary = true → choices i = false) :
    nonlinearBranchDenote gates choices ρ =
      fun w => conj (nonlinearBranchMatrix n gates choices) (ρ w) := by
  induction gates generalizing choices ρ with
  | nil =>
      funext w
      simp [nonlinearBranchDenote, nonlinearBranchMatrix, conj_one]
  | cons gate gates ih =>
      have htail : ∀ g ∈ gates, ∀ q c, g = Gate.measure q c → m ≤ c :=
        fun g hg q c heq => hdiscard g (by simp [hg]) q c heq
      have hhead : ∀ q c, gate = Gate.measure q c → m ≤ c :=
        hdiscard gate (by simp)
      have htailChoices : ∀ i (hi : i < gates.length),
          (gates[i]).isUnitary = true → choices (i + 1) = false := by
        intro i hi hu
        exact hchoices (i + 1) (by simp; omega) (by simpa using hu)
      have hheadChoice : gate.isUnitary = true → choices 0 = false := by
        intro hu
        exact hchoices 0 (by simp) (by simpa using hu)
      rw [nonlinearBranchDenote, ih _ _ htail htailChoices]
      funext w
      have hstep : nonlinearBranchStep gate (choices 0) ρ w =
          conj (nonlinearBranchOp n gate (choices 0)) (ρ w) := by
        cases gate with
        | measure q c =>
            simp [nonlinearBranchStep, nonlinearBranchOp,
              Nat.not_lt.mpr (hhead q c rfl)]
        | reset q => rfl
        | h q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | x q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | z q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | s q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | sdg q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | t q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | tdg q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | rz angle q => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | cnot c t => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | cz c t => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | ccx c₁ c₂ t => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
        | ccz c₁ c₂ t => simp [nonlinearBranchStep, nonlinearBranchOp, hheadChoice rfl]
      rw [hstep, nonlinearBranchMatrix, conj_mul]

theorem nonlinearBranchStep_matrix {n m : Nat} (gate : Gate) (choice : Bool)
    (ρ : CQState n m) (hdiscard : ∀ q c, gate = Gate.measure q c → m ≤ c) :
    nonlinearBranchStep gate choice ρ =
      fun w => conj (nonlinearEffectiveBranchOp n gate choice) (ρ w) := by
  funext w
  cases gate with
  | measure q c =>
      simp [nonlinearBranchStep, nonlinearEffectiveBranchOp, nonlinearBranchOp,
        Gate.isUnitary, Gate.isMeasurement,
        Nat.not_lt.mpr (hdiscard q c rfl)]
  | reset q =>
      simp [nonlinearBranchStep, nonlinearEffectiveBranchOp, nonlinearBranchOp,
        Gate.isUnitary, Gate.isMeasurement]
  | h q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | x q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | z q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | s q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | sdg q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | t q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | tdg q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | rz angle q => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | cnot c t => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | cz c t => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | ccx c₁ c₂ t => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]
  | ccz c₁ c₂ t => cases choice <;> simp [nonlinearBranchStep, nonlinearEffectiveBranchOp,
      nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement, conj]

def nonlinearBranchMatrices (n : Nat) : List Gate → List (Density n)
  | [] => [1]
  | gate :: gates =>
      (nonlinearBranchMatrices n gates).map
        (· * nonlinearEffectiveBranchOp n gate false) ++
      (nonlinearBranchMatrices n gates).map
        (· * nonlinearEffectiveBranchOp n gate true)

/-- For circuits without classical writes, every branch in the channel sum is
conjugation by its composed Kraus matrix. -/
theorem nonlinearBranchFamily_matrix {n m : Nat} (gates : List Gate)
    (ρ : CQState n m)
    (hdiscard : ∀ gate ∈ gates, ∀ q c, gate = Gate.measure q c → m ≤ c) :
    nonlinearBranchFamily gates ρ =
      (nonlinearBranchMatrices n gates).map (fun M w => conj M (ρ w)) := by
  induction gates generalizing ρ with
  | nil => simp [nonlinearBranchFamily, nonlinearBranchMatrices, conj_one]
  | cons gate gates ih =>
      have htail : ∀ g ∈ gates, ∀ q c, g = Gate.measure q c → m ≤ c :=
        fun g hg q c heq => hdiscard g (by simp [hg]) q c heq
      have hhead : ∀ q c, gate = Gate.measure q c → m ≤ c :=
        hdiscard gate (by simp)
      rw [nonlinearBranchFamily, ih _ htail, ih _ htail,
        nonlinearBranchStep_matrix gate false ρ hhead,
        nonlinearBranchStep_matrix gate true ρ hhead]
      simp only [nonlinearBranchMatrices, List.map_append, List.map_map]
      congr 1 <;> apply List.map_congr_left <;> intro M _ <;>
        funext w <;> simp [conj_mul]

/-- Exact Kraus expansion for circuits with no classical-memory writes. In particular,
it applies to arbitrary circuits when `m = 0`, including measurement and reset. -/
theorem nonlinear_denote_eq_branchMatrices_sum {n m : Nat} (gates : List Gate)
    (ρ : CQState n m)
    (hdiscard : ∀ gate ∈ gates, ∀ q c, gate = Gate.measure q c → m ≤ c) :
    denote gates ρ =
      ((nonlinearBranchMatrices n gates).map (fun M w => conj M (ρ w))).sum := by
  rw [nonlinear_denote_eq_branchFamily_sum, nonlinearBranchFamily_matrix gates ρ hdiscard]

/-- Every nonzero matrix in the channel expansion is one of the path matrices used by
the nonlinear polynomial phase-hop theorem. -/
theorem nonlinearBranchMatrices_path_matrix (n : Nat) (gates : List Gate) :
    ∀ M ∈ nonlinearBranchMatrices n gates,
      M = 0 ∨ ∃ choices : Nat → Bool, M = nonlinearBranchMatrix n gates choices := by
  induction gates with
  | nil =>
      intro M hM
      have : M = 1 := by simpa [nonlinearBranchMatrices] using hM
      exact Or.inr ⟨fun _ => false, by simp [this, nonlinearBranchMatrix]⟩
  | cons gate gates ih =>
      intro M hM
      simp only [nonlinearBranchMatrices, List.mem_append, List.mem_map] at hM
      rcases hM with ⟨tail, htail, rfl⟩ | ⟨tail, htail, rfl⟩
      · rcases ih tail htail with hzero | ⟨tailChoices, htailMatrix⟩
        · exact Or.inl (by rw [hzero, zero_mul])
        · refine Or.inr ⟨(fun | 0 => false | i + 1 => tailChoices i), ?_⟩
          simp [htailMatrix, nonlinearEffectiveBranchOp, nonlinearBranchMatrix,
            Gate.isUnitary]
      · rcases ih tail htail with hzero | ⟨tailChoices, htailMatrix⟩
        · exact Or.inl (by rw [hzero, zero_mul])
        · by_cases hu : gate.isUnitary = true
          · exact Or.inl (by simp only [nonlinearEffectiveBranchOp, Bool.and_true,
              hu, ite_true, mul_zero])
          · refine Or.inr ⟨(fun | 0 => true | i + 1 => tailChoices i), ?_⟩
            have hu' : gate.isUnitary = false := Bool.eq_false_iff.mpr hu
            have heff : nonlinearEffectiveBranchOp n gate true =
                nonlinearBranchOp n gate true := by
              simp only [nonlinearEffectiveBranchOp, Bool.and_true, hu', Bool.false_eq_true,
                ite_false]
            rw [heff, htailMatrix]
            rfl

theorem nonlinear_denote_matrix_sum_apply {n m : Nat} (gates : List Gate)
    (ρ : CQState n m) (w : Memory m)
    (hdiscard : ∀ gate ∈ gates, ∀ q c, gate = Gate.measure q c → m ≤ c) :
    denote gates ρ w =
      ((nonlinearBranchMatrices n gates).map (fun M => conj M (ρ w))).sum := by
  rw [nonlinear_denote_eq_branchMatrices_sum gates ρ hdiscard]
  induction nonlinearBranchMatrices n gates with
  | nil => rfl
  | cons M Ms ih => simp [ih]

/-- The nonlinear polynomial relation gives the phase-hop identity simultaneously
for every Kraus matrix contributing to the prefix and middle channels. -/
theorem nonlinearBranchMatrices_phase_hop {n : Nat} (pre middle : List Gate)
    (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial) :
    ∀ A ∈ nonlinearBranchMatrices n pre,
      ∀ B ∈ nonlinearBranchMatrices n middle,
        B * (phaseMatrix (rzPhase n θ q) * A) =
          phaseMatrix (rzPhase n (signedAngle sign θ) q') * (B * A) := by
  intro A hA B hB
  rcases nonlinearBranchMatrices_path_matrix n pre A hA with hzeroA | ⟨choicesA, hAeq⟩
  · simp [hzeroA]
  rcases nonlinearBranchMatrices_path_matrix n middle B hB with hzeroB | ⟨choicesB, hBeq⟩
  · simp [hzeroB]
  rw [hAeq, hBeq]
  exact nonlinearBranch_phase_hop pre middle choicesA choicesB θ q q' sign hpoly

theorem nonlinear_conj_list_sum {n : Nat} (U : Density n) (xs : List (Density n)) :
    conj U xs.sum = (xs.map (conj U)).sum := by
  induction xs with
  | nil => simp [conj]
  | cons x xs ih => simp [conj_add, ih]

theorem nonlinear_step_rz_phase {n m : Nat} (θ : ℚ) (q : Qubit)
    (ρ : CQState n m) :
    step (Gate.rz θ q) ρ = fun w => conj (phaseMatrix (rzPhase n θ q)) (ρ w) := by
  obtain ⟨c, hc, hphase⟩ := unitary_rot_smul n θ q
  funext w
  simp only [step, hphase]
  exact conj_smul c _ _ hc

theorem nonlinear_step_rz_singleMemory {n m : Nat} (θ : ℚ) (q : Qubit)
    (w₀ : Memory m) (D : Density n) :
    step (Gate.rz θ q) (nonlinearSingleMemory w₀ D) =
      nonlinearSingleMemory w₀ (conj (phaseMatrix (rzPhase n θ q)) D) := by
  rw [nonlinear_step_rz_phase]
  funext w
  by_cases hw : w = w₀ <;> simp [nonlinearSingleMemory, hw, conj]

theorem nonlinearEffectiveBranch_phase_hop {n : Nat} (pre middle : List Gate)
    (preChoices middleChoices : Nat → Bool) (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial) :
    nonlinearEffectiveBranchMatrix n middle middleChoices *
        (phaseMatrix (rzPhase n θ q) * nonlinearEffectiveBranchMatrix n pre preChoices) =
      phaseMatrix (rzPhase n (signedAngle sign θ) q') *
        (nonlinearEffectiveBranchMatrix n middle middleChoices *
          nonlinearEffectiveBranchMatrix n pre preChoices) := by
  rcases nonlinearEffectiveBranchMatrix_eq_zero_or_branch n pre preChoices with hzero | heq
  · simp [hzero]
  rcases nonlinearEffectiveBranchMatrix_eq_zero_or_branch n middle middleChoices with hzero | heq'
  · simp [hzero]
  rw [heq, heq']
  exact nonlinearBranch_phase_hop pre middle preChoices middleChoices θ q q' sign hpoly

/-- The hop preserves each selected branch's output memory as well as its density matrix. -/
theorem nonlinearBranchDenote_phase_hop_singleMemory {n m : Nat}
    (pre middle : List Gate) (preChoices middleChoices : Nat → Bool)
    (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial)
    (w₀ : Memory m) (D : Density n) :
    nonlinearBranchDenote middle middleChoices
        (step (Gate.rz θ q)
          (nonlinearBranchDenote pre preChoices (nonlinearSingleMemory w₀ D))) =
      step (Gate.rz (signedAngle sign θ) q')
        (nonlinearBranchDenote middle middleChoices
          (nonlinearBranchDenote pre preChoices (nonlinearSingleMemory w₀ D))) := by
  rw [nonlinearBranchDenote_singleMemory pre preChoices w₀ D,
    nonlinear_step_rz_singleMemory,
    nonlinearBranchDenote_singleMemory middle middleChoices,
    nonlinearBranchDenote_singleMemory middle middleChoices,
    nonlinear_step_rz_singleMemory]
  congr 1
  rw [← conj_mul, ← conj_mul, ← conj_mul, ← conj_mul,
    Matrix.mul_assoc (nonlinearEffectiveBranchMatrix n middle middleChoices)
      (phaseMatrix (rzPhase n θ q))
      (nonlinearEffectiveBranchMatrix n pre preChoices),
    Matrix.mul_assoc (phaseMatrix (rzPhase n (signedAngle sign θ) q'))
      (nonlinearEffectiveBranchMatrix n middle middleChoices)
      (nonlinearEffectiveBranchMatrix n pre preChoices),
    nonlinearEffectiveBranch_phase_hop pre middle preChoices middleChoices θ q q' sign hpoly]

theorem nonlinearBranchDenote_phase_hop {n m : Nat}
    (pre middle : List Gate) (preChoices middleChoices : Nat → Bool)
    (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial)
    (ρ : CQState n m) :
    nonlinearBranchDenote middle middleChoices
        (step (Gate.rz θ q) (nonlinearBranchDenote pre preChoices ρ)) =
      step (Gate.rz (signedAngle sign θ) q')
        (nonlinearBranchDenote middle middleChoices
          (nonlinearBranchDenote pre preChoices ρ)) := by
  conv_lhs => rw [← nonlinear_sum_singleMemory ρ]
  conv_rhs => rw [← nonlinear_sum_singleMemory ρ]
  simp only [nonlinearBranchDenote_sum, nonlinear_step_sum]
  apply Finset.sum_congr rfl
  intro w₀ _
  exact nonlinearBranchDenote_phase_hop_singleMemory pre middle preChoices
    middleChoices θ q q' sign hpoly w₀ (ρ w₀)

/-- The nonlinear phase hop preserves the full classical-quantum channel, including
measurements that overwrite classical bits and reset. -/
theorem nonlinear_phase_hop_channel {n m : Nat}
    (pre middle : List Gate) (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial) :
    Equivalent n m (pre ++ Gate.rz θ q :: middle)
      (pre ++ middle ++ [Gate.rz (signedAngle sign θ) q']) := by
  intro ρ
  simp only [denote_append, denote_cons, denote_nil]
  simp only [nonlinear_denote_eq_branchChoices_sum pre ρ]
  conv_lhs =>
    rw [nonlinear_step_list_sum, nonlinear_denote_list_sum]
  conv_rhs =>
    rw [nonlinear_denote_list_sum, nonlinear_step_list_sum]
  simp only [List.map_map]
  apply congrArg List.sum
  apply List.map_congr_left
  intro preChoices _
  change denote middle (step (Gate.rz θ q)
      (nonlinearBranchDenote pre preChoices ρ)) =
    step (Gate.rz (signedAngle sign θ) q')
      (denote middle (nonlinearBranchDenote pre preChoices ρ))
  rw [nonlinear_denote_eq_branchChoices_sum middle,
    nonlinear_denote_eq_branchChoices_sum middle, nonlinear_step_list_sum]
  apply congrArg List.sum
  rw [List.map_map]
  apply List.map_congr_left
  intro middleChoices _
  exact nonlinearBranchDenote_phase_hop pre middle preChoices middleChoices
    θ q q' sign hpoly ρ

theorem nonlinearEffectiveBranch_constant_phase {n : Nat} (pre : List Gate)
    (choices : Nat → Bool) (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) (ρ : Density n) :
    conj (phaseMatrix (rzPhase n θ q) *
      nonlinearEffectiveBranchMatrix n pre choices) ρ =
      conj (nonlinearEffectiveBranchMatrix n pre choices) ρ := by
  rcases nonlinearEffectiveBranchMatrix_eq_zero_or_branch n pre choices with hzero | heq
  · simp [hzero, conj]
  · rw [heq, nonlinearBranch_constant_phase pre choices θ q bit hpoly]
    split <;> exact conj_smul _ _ _ (ep_mul_star _)

theorem nonlinearBranchDenote_constant_phase_singleMemory {n m : Nat}
    (pre : List Gate) (choices : Nat → Bool) (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) (w₀ : Memory m) (D : Density n) :
    step (Gate.rz θ q)
        (nonlinearBranchDenote pre choices (nonlinearSingleMemory w₀ D)) =
      nonlinearBranchDenote pre choices (nonlinearSingleMemory w₀ D) := by
  rw [nonlinearBranchDenote_singleMemory,
    nonlinear_step_rz_singleMemory]
  congr 1
  rw [← conj_mul]
  exact nonlinearEffectiveBranch_constant_phase pre choices θ q bit hpoly D

theorem nonlinearBranchDenote_constant_phase {n m : Nat}
    (pre : List Gate) (choices : Nat → Bool) (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) (ρ : CQState n m) :
    step (Gate.rz θ q) (nonlinearBranchDenote pre choices ρ) =
      nonlinearBranchDenote pre choices ρ := by
  conv_lhs => rw [← nonlinear_sum_singleMemory ρ]
  conv_rhs => rw [← nonlinear_sum_singleMemory ρ]
  simp only [nonlinearBranchDenote_sum, nonlinear_step_sum]
  apply Finset.sum_congr rfl
  intro w₀ _
  exact nonlinearBranchDenote_constant_phase_singleMemory pre choices θ q bit
    hpoly w₀ (ρ w₀)

/-- Constant-phase elimination preserves the full channel, even after measurements
that write classical memory. -/
theorem nonlinear_constant_phase_channel {n m : Nat}
    (pre : List Gate) (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) :
    Equivalent n m (pre ++ [Gate.rz θ q]) pre := by
  intro ρ
  rw [denote_append, denote_cons, denote_nil,
    nonlinear_denote_eq_branchChoices_sum pre ρ,
    nonlinear_step_list_sum]
  apply congrArg List.sum
  rw [List.map_map]
  apply List.map_congr_left
  intro choices _
  exact nonlinearBranchDenote_constant_phase pre choices θ q bit hpoly ρ

theorem nonlinear_conj_phase_hop_sum {n : Nat} (As Bs : List (Density n))
    (D D' ρ : Density n)
    (hhop : ∀ A ∈ As, ∀ B ∈ Bs, B * (D * A) = D' * (B * A)) :
    (Bs.map (fun B => conj B (conj D ((As.map (fun A => conj A ρ)).sum)))).sum =
      conj D' ((Bs.map (fun B => conj B ((As.map (fun A => conj A ρ)).sum))).sum) := by
  rw [nonlinear_conj_list_sum D' _]
  apply congrArg List.sum
  rw [List.map_map]
  apply List.map_congr_left
  intro B hB
  change conj B (conj D ((As.map (fun A => conj A ρ)).sum)) =
    conj D' (conj B ((As.map (fun A => conj A ρ)).sum))
  conv_lhs =>
    rw [nonlinear_conj_list_sum D _]
    rw [nonlinear_conj_list_sum B _]
  conv_rhs =>
    rw [nonlinear_conj_list_sum B _]
    rw [nonlinear_conj_list_sum D' _]
  simp only [List.map_map]
  apply congrArg List.sum
  apply List.map_congr_left
  intro A hA
  change conj B (conj D (conj A ρ)) = conj D' (conj B (conj A ρ))
  rw [← conj_mul, ← conj_mul, ← conj_mul, ← conj_mul,
    Matrix.mul_assoc B D A, Matrix.mul_assoc D' B A, hhop A hA B hB]

theorem nonlinear_denote_composed_matrix_sum {n m : Nat} (pre middle : List Gate)
    (ρ : CQState n m) (w : Memory m)
    (hpre : ∀ gate ∈ pre, ∀ q c, gate = Gate.measure q c → m ≤ c)
    (hmiddle : ∀ gate ∈ middle, ∀ q c, gate = Gate.measure q c → m ≤ c) :
    denote (pre ++ middle) ρ w =
      ((nonlinearBranchMatrices n middle).map (fun B => conj B
        (((nonlinearBranchMatrices n pre).map (fun A => conj A (ρ w))).sum))).sum := by
  rw [denote_append, nonlinear_denote_matrix_sum_apply middle _ w hmiddle,
    nonlinear_denote_matrix_sum_apply pre ρ w hpre]

theorem nonlinear_denote_phase_middle_matrix_sum {n m : Nat} (pre middle : List Gate)
    (ρ : CQState n m) (w : Memory m) (θ : ℚ) (q : Qubit)
    (hpre : ∀ gate ∈ pre, ∀ q c, gate = Gate.measure q c → m ≤ c)
    (hmiddle : ∀ gate ∈ middle, ∀ q c, gate = Gate.measure q c → m ≤ c) :
    denote (pre ++ Gate.rz θ q :: middle) ρ w =
      ((nonlinearBranchMatrices n middle).map (fun B => conj B
        (conj (phaseMatrix (rzPhase n θ q))
          (((nonlinearBranchMatrices n pre).map (fun A => conj A (ρ w))).sum)))).sum := by
  rw [denote_append, denote_cons, nonlinear_denote_matrix_sum_apply middle _ w hmiddle,
    nonlinear_step_rz_phase θ q]
  dsimp only
  rw [nonlinear_denote_matrix_sum_apply pre ρ w hpre]

theorem nonlinear_denote_middle_phase_matrix_sum {n m : Nat} (pre middle : List Gate)
    (ρ : CQState n m) (w : Memory m) (θ : ℚ) (q : Qubit)
    (hpre : ∀ gate ∈ pre, ∀ q c, gate = Gate.measure q c → m ≤ c)
    (hmiddle : ∀ gate ∈ middle, ∀ q c, gate = Gate.measure q c → m ≤ c) :
    denote (pre ++ middle ++ [Gate.rz θ q]) ρ w =
      conj (phaseMatrix (rzPhase n θ q))
        (((nonlinearBranchMatrices n middle).map (fun B => conj B
          (((nonlinearBranchMatrices n pre).map (fun A => conj A (ρ w))).sum))).sum) := by
  rw [denote_append, denote_cons, denote_nil, nonlinear_step_rz_phase θ q]
  dsimp only
  rw [nonlinear_denote_composed_matrix_sum pre middle ρ w hpre hmiddle]

/-- Nonlinear phase folding is channel-correct across CCX, reset, and measurements
whose outcomes are discarded. This is the full phase-hop rule for circuits without
classical-memory writes; in particular it applies with `m = 0`. -/
theorem nonlinear_phase_hop_channel_no_classical_write {n m : Nat}
    (pre middle : List Gate) (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpre : ∀ gate ∈ pre, ∀ r c, gate = Gate.measure r c → m ≤ c)
    (hmiddle : ∀ gate ∈ middle, ∀ r c, gate = Gate.measure r c → m ≤ c)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial) :
    Equivalent n m (pre ++ Gate.rz θ q :: middle)
      (pre ++ middle ++ [Gate.rz (signedAngle sign θ) q']) := by
  intro ρ
  funext w
  rw [nonlinear_denote_phase_middle_matrix_sum pre middle ρ w θ q hpre hmiddle,
    nonlinear_denote_middle_phase_matrix_sum pre middle ρ w
      (signedAngle sign θ) q' hpre hmiddle]
  exact nonlinear_conj_phase_hop_sum
    (nonlinearBranchMatrices n pre) (nonlinearBranchMatrices n middle)
    (phaseMatrix (rzPhase n θ q))
    (phaseMatrix (rzPhase n (signedAngle sign θ) q')) (ρ w)
    (nonlinearBranchMatrices_phase_hop pre middle θ q q' sign hpoly)

/-- The constant-polynomial rule is valid on every contributing Kraus branch,
including branches created by reset. -/
theorem nonlinearBranchMatrices_constant_phase {n : Nat} (pre : List Gate)
    (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) :
    ∀ A ∈ nonlinearBranchMatrices n pre, ∀ ρ : Density n,
      conj (phaseMatrix (rzPhase n θ q) * A) ρ = conj A ρ := by
  intro A hA ρ
  rcases nonlinearBranchMatrices_path_matrix n pre A hA with hzero | ⟨choices, hAeq⟩
  · simp [hzero, conj]
  · rw [hAeq]
    rw [nonlinearBranch_constant_phase pre choices θ q bit hpoly]
    split <;> exact conj_smul _ _ _ (ep_mul_star _)

/-- A rotation on a formally constant wire has no effect on the channel, for circuits
with no classical-memory writes. -/
theorem nonlinear_constant_phase_channel_no_classical_write {n m : Nat}
    (pre : List Gate) (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpre : ∀ gate ∈ pre, ∀ r c, gate = Gate.measure r c → m ≤ c)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) :
    Equivalent n m (pre ++ [Gate.rz θ q]) pre := by
  intro ρ
  funext w
  rw [denote_append, denote_cons, denote_nil, nonlinear_step_rz_phase θ q]
  dsimp only
  rw [nonlinear_denote_matrix_sum_apply pre ρ w hpre,
    nonlinear_conj_list_sum]
  apply congrArg List.sum
  rw [List.map_map]
  apply List.map_congr_left
  intro A hA
  change conj (phaseMatrix (rzPhase n θ q)) (conj A (ρ w)) = conj A (ρ w)
  rw [← conj_mul]
  exact nonlinearBranchMatrices_constant_phase pre θ q bit hpoly A hA (ρ w)

/-- Regression: a phase hops through a cancelling pair of nonlinear CCX gates. -/
example (θ : ℚ) : Equivalent 3 0
    ([Gate.rz θ 2, Gate.ccx 0 1 2, Gate.ccx 0 1 2])
    ([Gate.ccx 0 1 2, Gate.ccx 0 1 2, Gate.rz θ 2]) := by
  apply nonlinear_phase_hop_channel_no_classical_write []
    [Gate.ccx 0 1 2, Gate.ccx 0 1 2] θ 2 2 false
  · simp
  · simp
  · simp [NonlinearAState.steps, NonlinearAState.step,
      NonlinearAState.wireOf, NonlinearAState.initial,
      TrackedPolynomial.add, TrackedPolynomial.mul?, TrackedPolynomial.fresh,
      Fingerprint.maxDegree]
    have htwo : (2 : BoolPolynomial) = 0 := by
      rw [← one_add_one_eq_two, ← MvPolynomial.C_1, ← MvPolynomial.C_add]
      simp [CharTwo.add_self_eq_zero]
    rw [add_assoc, ← two_mul, htwo, zero_mul, add_zero]

/-- Regression: reset makes the wire constant, so its following phase is global. -/
example (θ : ℚ) : Equivalent 1 0
    [Gate.reset 0, Gate.rz θ 0] [Gate.reset 0] := by
  apply nonlinear_constant_phase_channel_no_classical_write [Gate.reset 0] θ 0 false
  · simp
  · simp [NonlinearAState.steps, NonlinearAState.step, NonlinearAState.wireOf,
      NonlinearAState.initial, TrackedPolynomial.zero, BoolPolynomial.const, bit]

/-- Regression: moving a phase past measurement preserves the recorded outcome. -/
example (θ : ℚ) : Equivalent 1 1
    [Gate.rz θ 0, Gate.measure 0 0]
    [Gate.measure 0 0, Gate.rz θ 0] := by
  apply nonlinear_phase_hop_channel [] [Gate.measure 0 0] θ 0 0 false
  simp [NonlinearAState.steps, NonlinearAState.step]

/-- Regression: constant-phase removal remains valid after a measurement write. -/
example (θ : ℚ) : Equivalent 1 1
    [Gate.measure 0 0, Gate.reset 0, Gate.rz θ 0]
    [Gate.measure 0 0, Gate.reset 0] := by
  apply nonlinear_constant_phase_channel
    [Gate.measure 0 0, Gate.reset 0] θ 0 false
  simp [NonlinearAState.steps, NonlinearAState.step, NonlinearAState.wireOf,
    NonlinearAState.initial, TrackedPolynomial.zero, BoolPolynomial.const, bit]

/-- A phase hop established on one Kraus branch is preserved when that branch acts on
an arbitrary density matrix. -/
theorem nonlinearBranch_phase_hop_conj {n : Nat} (pre middle : List Gate)
    (preChoices middleChoices : Nat → Bool) (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial)
    (ρ : Density n) :
    conj (nonlinearBranchMatrix n middle middleChoices *
      (phaseMatrix (rzPhase n θ q) * nonlinearBranchMatrix n pre preChoices)) ρ =
    conj (phaseMatrix (rzPhase n (signedAngle sign θ) q') *
      (nonlinearBranchMatrix n middle middleChoices *
        nonlinearBranchMatrix n pre preChoices)) ρ := by
  rw [nonlinearBranch_phase_hop pre middle preChoices middleChoices θ q q' sign hpoly]

/-- A phase on a formally constant wire disappears from a Kraus branch after
conjugation, because its common scalar has unit norm. -/
theorem nonlinearBranch_constant_phase_conj {n : Nat} (pre : List Gate)
    (choices : Nat → Bool) (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) (ρ : Density n) :
    conj (phaseMatrix (rzPhase n θ q) * nonlinearBranchMatrix n pre choices) ρ =
      conj (nonlinearBranchMatrix n pre choices) ρ := by
  rw [nonlinearBranch_constant_phase pre choices θ q bit hpoly]
  split <;> exact conj_smul _ _ _ (ep_mul_star _)

end

end TzapLean
