import TzapLean.NonlinearSoundness
import TzapLean.Merge

/-!
# Kraus-branch paths for nonlinear phase folding

One selected Kraus operator per gate gives a matrix branch of the channel. A nonzero
entry of that matrix determines a `NonlinearPath`, including measurement and reset.
The polynomial phase-hop relation therefore applies to each branch matrix.
-/

namespace TzapLean

noncomputable section

def nonlinearBranchOp (n : Nat) (gate : Gate) (choice : Bool) : Density n :=
  match gate with
  | .measure q _ => proj n q choice
  | .reset q => resetKraus n q choice
  | _ => gateUnitary n gate

def nonlinearBranchMatrix (n : Nat) : List Gate → (Nat → Bool) → Density n
  | [], _ => 1
  | gate :: gates, choices =>
      nonlinearBranchMatrix n gates (fun i => choices (i + 1)) *
        nonlinearBranchOp n gate (choices 0)

theorem nonlinearBranchOp_unitary {n : Nat} {gate : Gate} (choice : Bool)
    (hunitary : gate.isUnitary = true) :
    nonlinearBranchOp n gate choice = gateUnitary n gate := by
  cases gate <;> simp [nonlinearBranchOp, Gate.isUnitary, Gate.isMeasurement] at hunitary ⊢

theorem nonlinearBranchMatrix_unitary {n : Nat} (gates : List Gate)
    (choices : Nat → Bool) (hunitary : ∀ gate ∈ gates, gate.isUnitary = true) :
    nonlinearBranchMatrix n gates choices = unitary n gates := by
  induction gates generalizing choices with
  | nil => rfl
  | cons gate gates ih =>
      rw [nonlinearBranchMatrix, unitary_cons,
        nonlinearBranchOp_unitary (choices 0) (hunitary gate (by simp))]
      rw [ih _ (fun later hlater => hunitary later (by simp [hlater]))]

theorem nonlinearBranchOp_possible {n : Nat} (gate : Gate) (choice : Bool)
    {before after : Basis n}
    (hne : nonlinearBranchOp n gate choice after before ≠ 0) :
    NonlinearAState.GatePossible n gate before after := by
  cases gate with
  | measure q c => exact NonlinearAState.measurementKraus_possible q c choice hne
  | reset q => exact NonlinearAState.resetKraus_possible q choice hne
  | h q => exact hne
  | x q => exact hne
  | z q => exact hne
  | s q => exact hne
  | sdg q => exact hne
  | t q => exact hne
  | tdg q => exact hne
  | rz angle q => exact hne
  | cnot c t => exact hne
  | cz c t => exact hne
  | ccx c₁ c₂ t => exact hne
  | ccz c₁ c₂ t => exact hne

theorem nonlinearBranchMatrix_path {n : Nat} (gates : List Gate)
    (choices : Nat → Bool) {before after : Basis n}
    (hne : nonlinearBranchMatrix n gates choices after before ≠ 0) :
    NonlinearAState.NonlinearPath n gates before after := by
  induction gates generalizing choices before with
  | nil =>
      have heq : after = before := AState.one_entry (by simpa [nonlinearBranchMatrix] using hne)
      subst after
      exact .nil before
  | cons gate gates ih =>
      simp only [nonlinearBranchMatrix, Matrix.mul_apply] at hne
      obtain ⟨middle, -, hmiddle⟩ := Finset.exists_ne_zero_of_sum_ne_zero hne
      have htail :
          nonlinearBranchMatrix n gates (fun i => choices (i + 1)) after middle ≠ 0 :=
        left_ne_zero_of_mul hmiddle
      have hhead : nonlinearBranchOp n gate (choices 0) middle before ≠ 0 :=
        right_ne_zero_of_mul hmiddle
      exact .cons (nonlinearBranchOp_possible gate (choices 0) hhead) (ih _ htail)

/-- The phase-hop matrix identity holds on each selected Kraus branch of a circuit,
including prefixes with reset and middle segments with measurement. -/
theorem nonlinearBranch_phase_hop {n : Nat} (pre middle : List Gate)
    (preChoices middleChoices : Nat → Bool) (θ : ℚ) (q q' : Qubit) (sign : Bool)
    (hpoly : ((((NonlinearAState.initial n).steps pre).steps middle).wireOf q').polynomial =
      if sign then (((NonlinearAState.initial n).steps pre).wireOf q).polynomial.flip
      else (((NonlinearAState.initial n).steps pre).wireOf q).polynomial) :
    nonlinearBranchMatrix n middle middleChoices *
        (phaseMatrix (rzPhase n θ q) * nonlinearBranchMatrix n pre preChoices) =
      phaseMatrix (rzPhase n (signedAngle sign θ) q') *
        (nonlinearBranchMatrix n middle middleChoices *
          nonlinearBranchMatrix n pre preChoices) := by
  apply phase_hop
  intro input before after hpre hmiddle
  exact NonlinearAState.path_of_polynomial_eq_after_prefix
    (nonlinearBranchMatrix_path pre preChoices hpre)
    (nonlinearBranchMatrix_path middle middleChoices hmiddle) hpoly

/-- A formally constant wire makes an `rz` a common scalar on every nonzero entry of
each selected branch matrix. This is the branch-level justification for constant-phase
elimination, including constants created by reset. -/
theorem nonlinearBranch_constant_phase {n : Nat} (pre : List Gate)
    (choices : Nat → Bool) (θ : ℚ) (q : Qubit) (bit : Bool)
    (hpoly : (((NonlinearAState.initial n).steps pre).wireOf q).polynomial =
      BoolPolynomial.const bit) :
    phaseMatrix (rzPhase n θ q) * nonlinearBranchMatrix n pre choices =
      (if bit then ep (θ / 2) else ep (-θ / 2)) •
        nonlinearBranchMatrix n pre choices := by
  ext out inp
  rw [phaseMatrix_mul_apply, Matrix.smul_apply]
  by_cases hentry : nonlinearBranchMatrix n pre choices out inp = 0
  · simp [hentry]
  obtain ⟨valuation, hconsistent⟩ := NonlinearAState.path_from_initial_consistent
    (nonlinearBranchMatrix_path pre choices hentry)
  have hbit : out.get q = bit := by
    have hq := hconsistent q
    rw [hpoly] at hq
    simpa using hq
  simp [rzPhase, hbit]

end

end TzapLean
