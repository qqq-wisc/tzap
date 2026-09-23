import TzapLean.Pass
import TzapLean.SuperOptProof

/-!
# Exact native controlled-gate decompositions

The executable rewrites are kept separate from their `Pass` packaging so the local matrix
identities and the list-level composition proof remain visible.  Rz decomposition is
deliberately absent from tzap-lean.
-/

namespace TzapLean

/-- The 15-gate Clifford+T Toffoli decomposition used by the Rust implementation. -/
def ccxDecomposition (control₁ control₂ target : Qubit) : List Gate :=
  [ .h target,
    .cnot control₂ target,
    .tdg target,
    .cnot control₁ target,
    .t target,
    .cnot control₂ target,
    .tdg target,
    .cnot control₁ target,
    .t control₂,
    .t target,
    .h target,
    .cnot control₁ control₂,
    .t control₁,
    .tdg control₂,
    .cnot control₁ control₂ ]

/-- Lower one CCX or CCZ; every other gate is preserved. -/
def decomposeToffoliGate : Gate → List Gate
  | .ccx a b t => ccxDecomposition a b t
  | .ccz a b t => [.h t] ++ ccxDecomposition a b t ++ [.h t]
  | gate => [gate]

def decomposeToffoliGates (gates : List Gate) : List Gate :=
  gates.flatMap decomposeToffoliGate

/-- Lower one CZ to `H(target); CX(control,target); H(target)`. -/
def decomposeCzGate : Gate → List Gate
  | .cz control target => [.h target, .cnot control target, .h target]
  | gate => [gate]

def decomposeCzGates (gates : List Gate) : List Gate :=
  gates.flatMap decomposeCzGate

def decomposeToffoli (circuit : RawCircuit) : RawCircuit :=
  circuit.withGates (decomposeToffoliGates circuit.gates)

def decomposeCz (circuit : RawCircuit) : RawCircuit :=
  circuit.withGates (decomposeCzGates circuit.gates)

/-! ## Canonical local matrix identities -/

private def ccxTarget : ExactMat 3 :=
  ExactMat.applyX [0, 1] 2 (ExactMat.id 3)

private def cczTarget : ExactMat 3 :=
  ExactMat.applyPhaseMask [0, 1, 2] 4 (ExactMat.id 3)

private def czTarget : ExactMat 2 :=
  ExactMat.applyPhaseMask [0, 1] 4 (ExactMat.id 2)

private theorem ccxCandidateAccepted :
    accepts ccxTarget (ccxDecomposition 0 1 2) = true := by native_decide

private theorem cczCandidateAccepted :
    accepts cczTarget ([.h 2] ++ ccxDecomposition 0 1 2 ++ [.h 2]) = true := by native_decide

private theorem czCandidateAccepted :
    accepts czTarget [.h 1, .cnot 0 1, .h 1] = true := by native_decide

/-- The exact checker certifies the canonical Toffoli lowering up to an irrelevant global
phase. -/
theorem ccxDecomposition_local :
    ∃ phase : Nat, unitary 3 (ccxDecomposition 0 1 2) =
      ω ^ phase • gateUnitary 3 (.ccx 0 1 2) := by
  obtain ⟨_, _, matrix, phase, hmatrix, hphase⟩ := accepts_spec ccxCandidateAccepted
  refine ⟨phase, ?_⟩
  have hcand := ExactMat.matrixOf_sound
    (n := 3) (gs := ccxDecomposition 0 1 2) (by native_decide) hmatrix
  have hsource := ExactMat.matrixOf_sound
    (n := 3) (gs := [.ccx 0 1 2]) (by native_decide)
    (show ExactMat.matrixOf 3 [.ccx 0 1 2] = some ccxTarget from rfl)
  have hmatch := ExactMat.phaseMatch_sound hphase
  rw [ExactMat.interp_normalize, hcand, hsource] at hmatch
  simpa using hmatch

theorem cczDecomposition_local :
    ∃ phase : Nat, unitary 3 ([.h 2] ++ ccxDecomposition 0 1 2 ++ [.h 2]) =
      ω ^ phase • gateUnitary 3 (.ccz 0 1 2) := by
  obtain ⟨_, _, matrix, phase, hmatrix, hphase⟩ := accepts_spec cczCandidateAccepted
  refine ⟨phase, ?_⟩
  have hcand := ExactMat.matrixOf_sound
    (n := 3) (gs := [.h 2] ++ ccxDecomposition 0 1 2 ++ [.h 2]) (by native_decide) hmatrix
  have hsource := ExactMat.matrixOf_sound
    (n := 3) (gs := [.ccz 0 1 2]) (by native_decide)
    (show ExactMat.matrixOf 3 [.ccz 0 1 2] = some cczTarget from rfl)
  have hmatch := ExactMat.phaseMatch_sound hphase
  rw [ExactMat.interp_normalize, hcand, hsource] at hmatch
  simpa using hmatch

theorem czDecomposition_local :
    ∃ phase : Nat, unitary 2 [.h 1, .cnot 0 1, .h 1] =
      ω ^ phase • gateUnitary 2 (.cz 0 1) := by
  obtain ⟨_, _, matrix, phase, hmatrix, hphase⟩ := accepts_spec czCandidateAccepted
  refine ⟨phase, ?_⟩
  have hcand := ExactMat.matrixOf_sound
    (n := 2) (gs := [.h 1, .cnot 0 1, .h 1]) (by native_decide) hmatrix
  have hsource := ExactMat.matrixOf_sound
    (n := 2) (gs := [.cz 0 1]) (by native_decide)
    (show ExactMat.matrixOf 2 [.cz 0 1] = some czTarget from rfl)
  have hmatch := ExactMat.phaseMatch_sound hphase
  rw [ExactMat.interp_normalize, hcand, hsource] at hmatch
  simpa using hmatch

/-! ## Lifting canonical identities to physical wires -/

theorem localizeGates_map_globalize {support : List Qubit} (hnd : support.Nodup)
    {gates : List Gate} (hbound : ∀ g ∈ gates, ∀ q ∈ g.qubitsOf, q < support.length) :
    localizeGates support (gates.map (globalizeGate support)) = gates := by
  rw [localizeGates, List.map_map]
  refine Eq.trans (List.map_congr_left ?_) (List.map_id gates)
  intro gate hgate
  exact mapQubits_comp fun q hq =>
    localIdxD_eq (localIdx_getD_self hnd (hbound gate hgate q hq))

theorem globalized_qubits_mem {support : List Qubit} {gates : List Gate}
    (hbound : ∀ g ∈ gates, ∀ q ∈ g.qubitsOf, q < support.length) :
    ∀ g ∈ gates.map (globalizeGate support), ∀ q ∈ g.qubitsOf, q ∈ support := by
  intro gate hgate q hq
  rcases List.mem_map.1 hgate with ⟨localGate, hlocal, rfl⟩
  rw [globalizeGate, qubitsOf_mapQubits] at hq
  rcases List.mem_map.1 hq with ⟨i, hi, rfl⟩
  exact getD_mem (hbound localGate hlocal i hi)

/-- A local Clifford+T identity lifts to any distinct in-range physical support. -/
theorem equivalent_map_globalize_of_local_smul {n m : Nat} {support : List Qubit}
    (hnd : support.Nodup) (hrange : ∀ q ∈ support, q < n)
    {source replacement : List Gate}
    (hsourceBound : ∀ g ∈ source, ∀ q ∈ g.qubitsOf, q < support.length)
    (hreplacementBound : ∀ g ∈ replacement, ∀ q ∈ g.qubitsOf, q < support.length)
    (hsourceUnitary : ∀ g ∈ source, g.isUnitary = true)
    (hreplacementUnitary : ∀ g ∈ replacement, g.isUnitary = true)
    (phase : ℂ) (hphase : phase * star phase = 1)
    (hlocal : unitary support.length replacement = phase • unitary support.length source) :
    Equivalent n m (replacement.map (globalizeGate support))
      (source.map (globalizeGate support)) := by
  apply equivalent_of_local_smul hnd hrange
  · exact globalized_qubits_mem hsourceBound
  · exact globalized_qubits_mem hreplacementBound
  · intro g hg
    rcases List.mem_map.1 hg with ⟨localGate, hlocalGate, rfl⟩
    rw [globalizeGate, isUnitary_mapQubits]
    exact hsourceUnitary localGate hlocalGate
  · intro g hg
    rcases List.mem_map.1 hg with ⟨localGate, hlocalGate, rfl⟩
    rw [globalizeGate, isUnitary_mapQubits]
    exact hreplacementUnitary localGate hlocalGate
  · exact hphase
  · rw [localizeGates_map_globalize hnd hreplacementBound,
        localizeGates_map_globalize hnd hsourceBound]
    exact hlocal

@[simp] theorem ccxDecomposition_globalize (a b t : Qubit) :
    (ccxDecomposition 0 1 2).map (globalizeGate [a, b, t]) = ccxDecomposition a b t := by
  rfl

@[simp] theorem cczDecomposition_globalize (a b t : Qubit) :
    (([Gate.h 2] ++ ccxDecomposition 0 1 2 ++ [Gate.h 2]).map
      (globalizeGate [a, b, t])) =
      [Gate.h t] ++ ccxDecomposition a b t ++ [Gate.h t] := by
  rfl

@[simp] theorem czDecomposition_globalize (a b : Qubit) :
    ([Gate.h 1, Gate.cnot 0 1, Gate.h 1].map (globalizeGate [a, b])) =
      [Gate.h b, Gate.cnot a b, Gate.h b] := by
  rfl

@[simp] theorem ccxGate_globalize (a b t : Qubit) :
    globalizeGate [a, b, t] (.ccx 0 1 2) = .ccx a b t := by
  rfl

@[simp] theorem cczGate_globalize (a b t : Qubit) :
    globalizeGate [a, b, t] (.ccz 0 1 2) = .ccz a b t := by
  rfl

@[simp] theorem czGate_globalize (a b : Qubit) :
    globalizeGate [a, b] (.cz 0 1) = .cz a b := by
  rfl

@[simp] theorem hTwo_globalize (a b t : Qubit) :
    globalizeGate [a, b, t] (.h 2) = .h t := by
  rfl

@[simp] theorem hOne_globalize (a b : Qubit) :
    globalizeGate [a, b] (.h 1) = .h b := by
  rfl

@[simp] theorem cnotZeroOne_globalize (a b : Qubit) :
    globalizeGate [a, b] (.cnot 0 1) = .cnot a b := by
  rfl

private theorem ccxSourceBound :
    ∀ g ∈ [Gate.ccx 0 1 2], ∀ q ∈ g.qubitsOf, q < 3 := by native_decide

private theorem cczSourceBound :
    ∀ g ∈ [Gate.ccz 0 1 2], ∀ q ∈ g.qubitsOf, q < 3 := by native_decide

private theorem ccxReplacementBound :
    ∀ g ∈ ccxDecomposition 0 1 2, ∀ q ∈ g.qubitsOf, q < 3 := by native_decide

private theorem cczReplacementBound :
    ∀ g ∈ [Gate.h 2] ++ ccxDecomposition 0 1 2 ++ [Gate.h 2],
      ∀ q ∈ g.qubitsOf, q < 3 := by native_decide

private theorem czSourceBound :
    ∀ g ∈ [Gate.cz 0 1], ∀ q ∈ g.qubitsOf, q < 2 := by native_decide

private theorem czReplacementBound :
    ∀ g ∈ [Gate.h 1, Gate.cnot 0 1, Gate.h 1], ∀ q ∈ g.qubitsOf, q < 2 := by
  native_decide

theorem ccxDecomposition_equivalent {n m : Nat} {a b t : Qubit}
    (hwf : (Gate.ccx a b t).Wf) (hrange : (Gate.ccx a b t).InRange n m) :
    Equivalent n m (ccxDecomposition a b t) [.ccx a b t] := by
  obtain ⟨hab, hat, hbt⟩ := hwf
  have ha : a < n := hrange.qubits a (by simp [Gate.qubitsOf])
  have hb : b < n := hrange.qubits b (by simp [Gate.qubitsOf])
  have ht : t < n := hrange.qubits t (by simp [Gate.qubitsOf])
  obtain ⟨phase, hlocal⟩ := ccxDecomposition_local
  have hlift := equivalent_map_globalize_of_local_smul (n := n) (m := m)
    (support := [a, b, t]) (source := [.ccx 0 1 2])
    (replacement := ccxDecomposition 0 1 2)
    (by simp [hab, hat, hbt])
    (by intro q hq; simp only [List.mem_cons, List.not_mem_nil, or_false] at hq;
        rcases hq with rfl | rfl | rfl <;> assumption)
    (by simpa using ccxSourceBound) (by simpa using ccxReplacementBound)
    (by native_decide) (by native_decide) (ω ^ phase) (ExactMat.omega_pow_unit phase)
    (by simpa [unitary] using hlocal)
  simpa using hlift

theorem cczDecomposition_equivalent {n m : Nat} {a b t : Qubit}
    (hwf : (Gate.ccz a b t).Wf) (hrange : (Gate.ccz a b t).InRange n m) :
    Equivalent n m ([.h t] ++ ccxDecomposition a b t ++ [.h t]) [.ccz a b t] := by
  obtain ⟨hab, hat, hbt⟩ := hwf
  have ha : a < n := hrange.qubits a (by simp [Gate.qubitsOf])
  have hb : b < n := hrange.qubits b (by simp [Gate.qubitsOf])
  have ht : t < n := hrange.qubits t (by simp [Gate.qubitsOf])
  obtain ⟨phase, hlocal⟩ := cczDecomposition_local
  have hlift := equivalent_map_globalize_of_local_smul (n := n) (m := m)
    (support := [a, b, t]) (source := [.ccz 0 1 2])
    (replacement := [.h 2] ++ ccxDecomposition 0 1 2 ++ [.h 2])
    (by simp [hab, hat, hbt])
    (by intro q hq; simp only [List.mem_cons, List.not_mem_nil, or_false] at hq;
        rcases hq with rfl | rfl | rfl <;> assumption)
    (by simpa using cczSourceBound) (by simpa using cczReplacementBound)
    (by native_decide) (by native_decide) (ω ^ phase) (ExactMat.omega_pow_unit phase)
    (by simpa [unitary] using hlocal)
  simpa using hlift

theorem czDecomposition_equivalent {n m : Nat} {a b : Qubit}
    (hwf : (Gate.cz a b).Wf) (hrange : (Gate.cz a b).InRange n m) :
    Equivalent n m [.h b, .cnot a b, .h b] [.cz a b] := by
  obtain ⟨phase, hlocal⟩ := czDecomposition_local
  have ha : a < n := hrange.qubits a (by simp [Gate.qubitsOf])
  have hb : b < n := hrange.qubits b (by simp [Gate.qubitsOf])
  have hlift := equivalent_map_globalize_of_local_smul (n := n) (m := m)
    (support := [a, b]) (source := [.cz 0 1])
    (replacement := [.h 1, .cnot 0 1, .h 1])
    (by simpa [Gate.Wf] using hwf)
    (by intro q hq; simp only [List.mem_cons, List.not_mem_nil, or_false] at hq;
        rcases hq with rfl | rfl <;> assumption)
    (by simpa using czSourceBound) (by simpa using czReplacementBound)
    (by native_decide) (by native_decide) (ω ^ phase) (ExactMat.omega_pow_unit phase)
    (by simpa [unitary] using hlocal)
  simpa using hlift

/-- The Toffoli lowering introduces only valid gates when its source gate is valid. -/
theorem decomposeToffoliGate_wf {gate : Gate} (h : gate.Wf) :
    ∀ out ∈ decomposeToffoliGate gate, out.Wf := by
  cases gate with
  | ccx a b t | ccz a b t =>
      simp only [Gate.Wf] at h
      obtain ⟨hab, hat, hbt⟩ := h
      simp [decomposeToffoliGate, ccxDecomposition, Gate.Wf, hab, hat, hbt]
  | _ => simpa [decomposeToffoliGate] using h

theorem decomposeCzGate_wf {gate : Gate} (h : gate.Wf) :
    ∀ out ∈ decomposeCzGate gate, out.Wf := by
  cases gate with
  | cz a b => simpa [decomposeCzGate, Gate.Wf] using h
  | _ => simpa [decomposeCzGate] using h

private theorem ccxDecomposition_inRange {n m : Nat} {a b t : Qubit}
    (h : (Gate.ccx a b t).InRange n m) :
    ∀ out ∈ ccxDecomposition a b t, out.InRange n m := by
  intro out hout
  simp only [ccxDecomposition, List.mem_cons, List.not_mem_nil, or_false] at hout
  rcases hout with (rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl | rfl | rfl | rfl | rfl | rfl)
  all_goals
    constructor
    · intro q hq
      apply h.qubits q
      simp_all [Gate.qubitsOf] <;> tauto
    · simp [Gate.cbitsOf]

theorem decomposeToffoliGate_inRange {n m : Nat} {gate : Gate} (h : gate.InRange n m) :
    ∀ out ∈ decomposeToffoliGate gate, out.InRange n m := by
  cases gate with
  | ccx a b t => simpa [decomposeToffoliGate] using ccxDecomposition_inRange h
  | ccz a b t =>
      intro out hout
      change out ∈ ([Gate.h t] ++ ccxDecomposition a b t ++ [Gate.h t]) at hout
      rcases List.mem_append.mp hout with hleft | hright
      · rcases List.mem_append.mp hleft with hpref | hccx
        · have : out = Gate.h t := by simpa using hpref
          subst out
          constructor
          · intro q hq
            apply h.qubits q
            simp_all [Gate.qubitsOf]
          · simp [Gate.cbitsOf]
        · have hccxRange : (Gate.ccx a b t).InRange n m := ⟨h.qubits, h.cbits⟩
          exact ccxDecomposition_inRange hccxRange out hccx
      · have : out = Gate.h t := by simpa using hright
        subst out
        constructor
        · intro q hq
          apply h.qubits q
          simp_all [Gate.qubitsOf]
        · simp [Gate.cbitsOf]
  | _ => simpa [decomposeToffoliGate] using h

theorem decomposeCzGate_inRange {n m : Nat} {gate : Gate} (h : gate.InRange n m) :
    ∀ out ∈ decomposeCzGate gate, out.InRange n m := by
  cases gate with
  | cz a b =>
      intro out hout
      simp only [decomposeCzGate, List.mem_cons, List.not_mem_nil, or_false] at hout
      rcases hout with (rfl | rfl | rfl)
      all_goals
        constructor
        · intro q hq
          apply h.qubits q
          simp_all [Gate.qubitsOf]
        · simp [Gate.cbitsOf]
  | _ => simpa [decomposeCzGate] using h

theorem decomposeToffoliGates_wf {gates : List Gate} (h : ∀ g ∈ gates, g.Wf) :
    ∀ g ∈ decomposeToffoliGates gates, g.Wf := by
  intro g hg
  rw [decomposeToffoliGates, List.mem_flatMap] at hg
  obtain ⟨source, hs, hg⟩ := hg
  exact decomposeToffoliGate_wf (h source hs) g hg

theorem decomposeCzGates_wf {gates : List Gate} (h : ∀ g ∈ gates, g.Wf) :
    ∀ g ∈ decomposeCzGates gates, g.Wf := by
  intro g hg
  rw [decomposeCzGates, List.mem_flatMap] at hg
  obtain ⟨source, hs, hg⟩ := hg
  exact decomposeCzGate_wf (h source hs) g hg

theorem decomposeToffoliGates_inRange {n m : Nat} {gates : List Gate}
    (h : ∀ g ∈ gates, g.InRange n m) :
    ∀ g ∈ decomposeToffoliGates gates, g.InRange n m := by
  intro g hg
  rw [decomposeToffoliGates, List.mem_flatMap] at hg
  obtain ⟨source, hs, hg⟩ := hg
  exact decomposeToffoliGate_inRange (h source hs) g hg

theorem decomposeCzGates_inRange {n m : Nat} {gates : List Gate}
    (h : ∀ g ∈ gates, g.InRange n m) :
    ∀ g ∈ decomposeCzGates gates, g.InRange n m := by
  intro g hg
  rw [decomposeCzGates, List.mem_flatMap] at hg
  obtain ⟨source, hs, hg⟩ := hg
  exact decomposeCzGate_inRange (h source hs) g hg

theorem decomposeToffoliGate_equivalent {n m : Nat} {gate : Gate} (hwf : gate.Wf)
    (hrange : gate.InRange n m) : Equivalent n m (decomposeToffoliGate gate) [gate] := by
  cases gate with
  | ccx a b t => exact ccxDecomposition_equivalent hwf hrange
  | ccz a b t => exact cczDecomposition_equivalent hwf hrange
  | _ => exact Equivalent.refl _ _ _

theorem decomposeCzGate_equivalent {n m : Nat} {gate : Gate} (hwf : gate.Wf)
    (hrange : gate.InRange n m) : Equivalent n m (decomposeCzGate gate) [gate] := by
  cases gate with
  | cz a b => exact czDecomposition_equivalent hwf hrange
  | _ => exact Equivalent.refl _ _ _

theorem decomposeToffoliGates_correct {n m : Nat} : ∀ {gates : List Gate},
    (∀ g ∈ gates, g.Wf) → (∀ g ∈ gates, g.InRange n m) →
      Equivalent n m (decomposeToffoliGates gates) gates := by
  intro gates hwf hrange
  induction gates with
  | nil => exact Equivalent.refl _ _ _
  | cons gate gates ih =>
      rw [decomposeToffoliGates, List.flatMap_cons]
      have hgate := decomposeToffoliGate_equivalent
        (hwf gate (by simp)) (hrange gate (by simp))
      have htail := ih (fun g hg => hwf g (by simp [hg]))
        (fun g hg => hrange g (by simp [hg]))
      have htail' : Equivalent n m (List.flatMap decomposeToffoliGate gates) gates := by
        simpa [decomposeToffoliGates] using htail
      exact Equivalent.trans (Equivalent.append_right _ hgate)
        (by simpa using Equivalent.append_left [gate] htail')

theorem decomposeCzGates_correct {n m : Nat} : ∀ {gates : List Gate},
    (∀ g ∈ gates, g.Wf) → (∀ g ∈ gates, g.InRange n m) →
      Equivalent n m (decomposeCzGates gates) gates := by
  intro gates hwf hrange
  induction gates with
  | nil => exact Equivalent.refl _ _ _
  | cons gate gates ih =>
      rw [decomposeCzGates, List.flatMap_cons]
      have hgate := decomposeCzGate_equivalent (hwf gate (by simp)) (hrange gate (by simp))
      have htail := ih (fun g hg => hwf g (by simp [hg]))
        (fun g hg => hrange g (by simp [hg]))
      have htail' : Equivalent n m (List.flatMap decomposeCzGate gates) gates := by
        simpa [decomposeCzGates] using htail
      exact Equivalent.trans (Equivalent.append_right _ hgate)
        (by simpa using Equivalent.append_left [gate] htail')

/-- Checked lowering is intentionally a no-op for an internally constructed circuit whose
operands are out of range. Parsed QASM always takes the transforming branch. -/
def decomposeToffoliChecked (c : Circuit n m) : Circuit n m :=
  if h : c.raw.WellFormed then
    ⟨decomposeToffoli c.raw,
      by simpa [decomposeToffoli] using c.numQubits_eq,
      by simpa [decomposeToffoli] using c.numCbits_eq,
      decomposeToffoliGates_wf c.wf⟩
  else c

def decomposeCzChecked (c : Circuit n m) : Circuit n m :=
  if h : c.raw.WellFormed then
    ⟨decomposeCz c.raw,
      by simpa [decomposeCz] using c.numQubits_eq,
      by simpa [decomposeCz] using c.numCbits_eq,
      decomposeCzGates_wf c.wf⟩
  else c

theorem decomposeToffoliChecked_wellFormed (c : Circuit n m) (h : c.raw.WellFormed) :
    (decomposeToffoliChecked c).raw.WellFormed := by
  unfold decomposeToffoliChecked
  simp only [h, ↓reduceDIte]
  simpa [decomposeToffoli, RawCircuit.WellFormed] using
    (decomposeToffoliGates_inRange h)

theorem decomposeCzChecked_wellFormed (c : Circuit n m) (h : c.raw.WellFormed) :
    (decomposeCzChecked c).raw.WellFormed := by
  unfold decomposeCzChecked
  simp only [h, ↓reduceDIte]
  simpa [decomposeCz, RawCircuit.WellFormed] using (decomposeCzGates_inRange h)

theorem decomposeToffoliChecked_flagsOk (c : Circuit n m) (h : c.raw.FlagsOk) :
    (decomposeToffoliChecked c).raw.FlagsOk := by
  unfold decomposeToffoliChecked
  split
  · exact RawCircuit.flagsOk_withGates _ _
  · exact h

theorem decomposeCzChecked_flagsOk (c : Circuit n m) (h : c.raw.FlagsOk) :
    (decomposeCzChecked c).raw.FlagsOk := by
  unfold decomposeCzChecked
  split
  · exact RawCircuit.flagsOk_withGates _ _
  · exact h

def DecomposeToffoli : Pass where
  name := "Toffoli decomposition"
  run := decomposeToffoliChecked
  correct := by
    intro n m c
    unfold decomposeToffoliChecked
    split
    · rename_i h
      have hin : ∀ g ∈ c.raw.gates, g.InRange n m := by
        intro g hg
        have hg' := h g hg
        rw [c.numQubits_eq, c.numCbits_eq] at hg'
        exact hg'
      exact decomposeToffoliGates_correct c.wf hin
    · exact Equivalent.refl _ _ _

def DecomposeCz : Pass where
  name := "CZ decomposition"
  run := decomposeCzChecked
  correct := by
    intro n m c
    unfold decomposeCzChecked
    split
    · rename_i h
      have hin : ∀ g ∈ c.raw.gates, g.InRange n m := by
        intro g hg
        have hg' := h g hg
        rw [c.numQubits_eq, c.numCbits_eq] at hg'
        exact hg'
      exact decomposeCzGates_correct c.wf hin
    · exact Equivalent.refl _ _ _

end TzapLean
