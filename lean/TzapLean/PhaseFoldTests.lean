import TzapLean.PhaseFoldRand
import TzapLean.GF128
import TzapLean.Decompose

/-!
# `PhaseFoldRand`: the Rust test suite, ported

The practical behavioural tests from `src/phase_fold_rand.rs`, including the
Toffoli-decomposition pipelines, are reproduced below. The three ignored Rust cases are
replayed separately in `PhaseFoldSlowTests.lean`. Rust checks gate *counts* plus a numerical
`circuits_equiv`; these
`#guard`s pin down the exact gate list or count the executable pass produces. The
Clifford+T subset also has exact matrix-equivalence checks. Arbitrary-Rz examples use
rational angles in units of π, so Rust's per-case floating-point equivalence assertions
are not reproduced literally. `GF128Bridge.lean` supplies the nonlinear pass's general
packed-field probability bound.

Angles are rationals in units of `π`, so Rust's radian constants become the `π`-fraction with
the same classification: `0.3` (not a quarter turn) becomes `3/10`, `PI/4` becomes `1/4`.

The draws here are a fixed splitmix stream, so the results are reproducible. The optimizer
itself obtains fresh 128-bit samples from `IO.getRandomBytes`; the theorem models those as
independent uniform samples, while the platform RNG assumption remains explicit.
-/

namespace TzapLean

open Gate

/-! ### Nonlinear fingerprint arithmetic

These checks pin down the packed arithmetic, nonlinear transfer functions, and degree cutoff
used by both Rust and the Lean executable. The quotient-field refinement and sharp
probabilistic theorem are proved in `GF128Bridge.lean`. -/

#guard GF128.add 0x1234 0 == 0x1234
#guard GF128.mul 0 0xDEADBEEF == 0
#guard GF128.mul 1 0xDEADBEEF == 0xDEADBEEF
#guard GF128.mul (2 ^ 127) 2 == 0x87
#guard GF128.mul 0x123456789ABCDEF 0xFEDCBA987654321 ==
  GF128.mul 0xFEDCBA987654321 0x123456789ABCDEF
#guard GF128.mul 0x12345 (GF128.add 0xABC 0xDEF) ==
  GF128.add (GF128.mul 0x12345 0xABC) (GF128.mul 0x12345 0xDEF)

def gf128Samples : List Nat := [0, 1, 0x123456789ABCDEF0, 2 ^ 128 - 1, 2 ^ 127]

-- The complete finite sample matrix from Rust: zero/one, commutativity, and distributivity.
#guard gf128Samples.all fun a =>
  GF128.mul a 0 == 0 && GF128.mul a 1 == GF128.normalize a &&
    gf128Samples.all fun b =>
      GF128.mul a b == GF128.mul b a &&
        gf128Samples.all fun c =>
          GF128.mul a (GF128.add b c) == GF128.add (GF128.mul a b) (GF128.mul a c)

def fpA : Fingerprint := ⟨0x123456, 9⟩
def fpB : Fingerprint := ⟨0xABCDEF, 12⟩
def fpTarget : Fingerprint := ⟨0x55AA, 7⟩

-- Applying the same in-range CCX update twice restores the target field value.
#guard
  let once := Fingerprint.ccx fpTarget fpA fpB 0
  let twice := Fingerprint.ccx once fpA fpB 0
  twice.value == fpTarget.value

/-- Successive nonlinear products can grow like Fibonacci numbers. -/
def fibonacciDegrees : Nat → Nat × Nat
  | 0 => (1, 1)
  | n + 1 =>
      let (a, b) := fibonacciDegrees n
      (b, a + b)

-- The synthetic sequence stays within the budget at step 45 and crosses `2^32` at step 46.
#guard (fibonacciDegrees 45).2 ≤ Fingerprint.maxDegree
#guard (fibonacciDegrees 46).2 > Fingerprint.maxDegree
#guard
  let a : Fingerprint := { value := 3, degree := (fibonacciDegrees 45).1 }
  let b : Fingerprint := { value := 5, degree := (fibonacciDegrees 45).2 }
  let freshValue := 0xC0FFEE
  Fingerprint.ccx fpTarget a b freshValue == Fingerprint.fresh freshValue

/-- A fixed draw stream: splitmix64 bit mixing, one 63-bit tag per variable. Reproducible,
and *only* for these tests — see the module docstring. -/
def seedWords (k : Nat) (seed : Nat) : Nat → Tag := fun i =>
  let x : UInt64 := (seed.toUInt64 + i.toUInt64 + 1) * 0x9E3779B97F4A7C15
  let x := (x ^^^ (x >>> 30)) * 0xBF58476D1CE4E5B9
  let x := (x ^^^ (x >>> 27)) * 0x94D049BB133111EB
  let x := x ^^^ (x >>> 31)
  x.toNat % 2 ^ k

/-- A fixed draw stream for the tests below. -/
def testWords : Nat → Tag := seedWords 63 0

-- OS entropy is packed little-endian before the tag is exposed as bits.
#guard natOfBytes (ByteArray.mk #[0x01, 0x02, 0x03]) 0 3 == 0x030201

/-- Phase folding with those draws. -/
def pf (n : Nat) (gs : List Gate) : List Gate := phaseFoldGatesNonlinear testWords n gs

/-! ### Nonlinear CCX tracking -/

-- An identical CCX pair restores the target polynomial, so the surrounding T/Tdg cancel.
#guard pf 3 [.t 2, .ccx 0 1 2, .ccx 1 0 2, .tdg 2] ==
  [.ccx 0 1 2, .ccx 1 0 2]

-- One CCX changes the target to a genuinely nonlinear fingerprint and blocks the merge.
#guard ((pf 3 [.t 2, .ccx 0 1 2, .tdg 2]).filter fun g => (rotAngle g).isSome).length == 2

-- The transfer function uses target + left*right and the product-degree upper bound.
#guard
  let st := NState.initial testWords 3
  let after := st.step testWords (.ccx 0 1 2)
  (after.fpOf 2).value == GF128.add (st.fpOf 2).value
    (GF128.mul (st.fpOf 0).value (st.fpOf 1).value)
#guard
  let st := NState.initial testWords 3
  let after := st.step testWords (.ccx 0 1 2)
  (after.fpOf 2).degree == 2

-- Reset is the Boolean constant zero and does not consume a random variable.
#guard
  let st := (NState.initial testWords 1).step testWords (.reset 0)
  st.fpOf 0 == Fingerprint.zero && st.fresh == 1

-- Shared-control dependencies are not mistaken for an affine cancellation.
#guard pf 3
    [.t 2, .cnot 0 1, .ccx 0 1 2, .cnot 0 1, .ccx 0 1 2, .cnot 0 2, .tdg 2] ==
  [.t 2, .cnot 0 1, .ccx 0 1 2, .cnot 0 1, .ccx 0 1 2, .cnot 0 2, .tdg 2]

structure DegreeGrowth where
  degrees : Array Nat
  gates : List Gate
  overflowControls : Qubit × Qubit
  overflowTarget : Qubit
  crossed : Bool
deriving Repr

def leastDegreeWire (degrees : Array Nat) : Qubit :=
  if degrees[0]!.min degrees[1]! ≤ degrees[2]! then
    if degrees[0]! ≤ degrees[1]! then 0 else 1
  else 2

def otherWires : Qubit → Qubit × Qubit
  | 0 => (1, 2)
  | 1 => (0, 2)
  | _ => (0, 1)

/-- Build Rust's synthetic Fibonacci-degree circuit, with a hard recursion bound that is
itself checked below. -/
def growDegrees : Nat → DegreeGrowth → DegreeGrowth
  | 0, result => result
  | fuel + 1, result =>
      if result.crossed then result
      else
        let target := leastDegreeWire result.degrees
        let controls := otherWires target
        let productDegree := result.degrees[controls.1]! + result.degrees[controls.2]!
        if productDegree > Fingerprint.maxDegree then
          { result with overflowControls := controls, overflowTarget := target, crossed := true }
        else
          growDegrees fuel {
            degrees := result.degrees.set! target (max result.degrees[target]! productDegree)
            gates := result.gates ++ [.ccx controls.1 controls.2 target]
            overflowControls := controls
            overflowTarget := target
            crossed := false }

def extremeDegreeGrowth : DegreeGrowth :=
  growDegrees 64 {
    degrees := #[1, 1, 1]
    gates := []
    overflowControls := (0, 1)
    overflowTarget := 2
    crossed := false }

#guard extremeDegreeGrowth.crossed
#guard extremeDegreeGrowth.gates.length < 64
#guard extremeDegreeGrowth.degrees.all (· ≤ Fingerprint.maxDegree)

def extremeDegreeCircuit : List Gate :=
  let growth := extremeDegreeGrowth
  growth.gates ++ [.t growth.overflowTarget,
    .ccx growth.overflowControls.1 growth.overflowControls.2 growth.overflowTarget,
    .ccx growth.overflowControls.1 growth.overflowControls.2 growth.overflowTarget,
    .tdg growth.overflowTarget]

-- Each overflowing CCX refreshes the target independently, so folding stops conservatively.
#guard pf 3 extremeDegreeCircuit == extremeDegreeCircuit

/-! ### Replayable native-gate circuit sweep

Every generated circuit contains CZ, CCX, and CCZ as well as phase and classical-update gates.
The seed rotates gate and operand choices, giving a deterministic regression corpus. -/

def nativePhaseGate (seed i : Nat) : Gate :=
  let q := (seed + i) % 3
  match (seed * 7 + i) % 10 with
  | 0 => .t q
  | 1 => .tdg q
  | 2 => .x q
  | 3 => .h q
  | 4 => .cnot q ((q + 1) % 3)
  | 5 => .cz q ((q + 1) % 3)
  | 6 => .ccx q ((q + 1) % 3) ((q + 2) % 3)
  | 7 => .ccz q ((q + 1) % 3) ((q + 2) % 3)
  | 8 => .rz (3/10) q
  | _ => .s q

def nativePhaseCircuit (seed : Nat) : List Gate :=
  [.cz 0 1, .ccx 0 1 2, .ccz 0 1 2] ++
    (List.range 37).map (nativePhaseGate seed)

def nativePhaseFuzzCase (seed : Nat) : Bool :=
  let input := nativePhaseCircuit seed
  let output := pf 3 input
  input.any (fun g => g.kind == .ccx) && input.any (fun g => g.kind == .ccz) &&
    output.all Gate.Wf

#guard (List.range 256).all nativePhaseFuzzCase

/-- Exact semantic differential checks on a replayable Clifford+T subset of the generated
native-gate corpus. `ExactMat` intentionally treats arbitrary `rz` as a barrier. -/
def nativePhaseMatrixCircuit (seed : Nat) : List Gate :=
  [.cz 0 1, .ccx 0 1 2, .ccz 0 1 2] ++
    (List.range 20).map (fun i =>
      match nativePhaseGate seed i with
      | .rz _ q => .tdg q
      | gate => gate)

/-- Rust's `count_phase_gates`. -/
def countPhaseGates (gs : List Gate) : Nat := (gs.filter fun g => (rotAngle g).isSome).length

/-- Rust's `count_t_gates`. -/
def countTGates (gs : List Gate) : Nat :=
  (gs.filter fun g => match g with | .t _ | .tdg _ => true | _ => false).length

/-- Exact Clifford+T equivalence up to global phase, rather than Rust's floating-point
`circuits_equiv` tolerance. The checker is sound by `matrixOf_sound` and `phaseMatch_sound`. -/
def exactMatrixEquivalent (n : Nat) (source output : List Gate) : Bool :=
  match ExactMat.matrixOf n source, ExactMat.matrixOf n output with
  | some sourceMat, some outputMat =>
      (ExactMat.phaseMatch sourceMat.normalize outputMat.normalize).isSome
  | _, _ => false

/-- The executable exact-matrix check implies equality of the modeled channels. -/
theorem exactMatrixEquivalent_sound {n m : Nat} {source output : List Gate}
    (hsourceWf : ∀ gate ∈ source, gate.Wf)
    (houtputWf : ∀ gate ∈ output, gate.Wf)
    (hsourceUnitary : ∀ gate ∈ source, gate.isUnitary = true)
    (houtputUnitary : ∀ gate ∈ output, gate.isUnitary = true)
    (h : exactMatrixEquivalent n source output = true) :
    Equivalent n m output source := by
  unfold exactMatrixEquivalent at h
  cases hsource : ExactMat.matrixOf n source with
  | none => simp [hsource] at h
  | some sourceMat =>
      cases houtput : ExactMat.matrixOf n output with
      | none => simp [hsource, houtput] at h
      | some outputMat =>
          cases hphase : ExactMat.phaseMatch sourceMat.normalize outputMat.normalize with
          | none => simp [hsource, houtput, hphase] at h
          | some phase =>
              have hmat : unitary n output = ω ^ phase • unitary n source := by
                have hs := ExactMat.matrixOf_sound hsourceWf hsource
                have ho := ExactMat.matrixOf_sound houtputWf houtput
                have hp := ExactMat.phaseMatch_sound hphase
                simpa [hs, ho] using hp
              exact equivalent_of_unitary_smul houtputUnitary hsourceUnitary
                (ω ^ phase) (ExactMat.omega_pow_unit phase) hmat

def nativePhaseMatrixCase (seed : Nat) : Bool :=
  let input := nativePhaseMatrixCircuit seed
  exactMatrixEquivalent 3 input (pf 3 input)

#guard (List.range 64).all nativePhaseMatrixCase

-- The CCX inverse-pair regression is exact-equivalent, not merely equal on the test seed.
#guard
  let input : List Gate := [.t 2, .ccx 0 1 2, .ccx 1 0 2, .tdg 2]
  exactMatrixEquivalent 3 input (pf 3 input)

-- The shared-control regression keeps its original channel.
#guard
  let input : List Gate :=
    [.t 2, .cnot 0 1, .ccx 0 1 2, .cnot 0 1, .ccx 0 1 2, .cnot 0 2, .tdg 2]
  exactMatrixEquivalent 3 input (pf 3 input)

/-! ### Decompose-then-fold pipelines -/

def sharedControlToffolis : List Gate :=
  [.ccx 1 2 0, .cnot 2 0, .ccx 0 1 2]

def smallToffoliPipeline : List Gate :=
  [.x 4, .h 4, .h 4, .ccx 0 3 4, .h 4, .h 4, .ccx 2 3 4, .h 4, .h 4,
    .cnot 3 4, .h 4, .h 4, .ccx 1 2 4, .h 4, .h 4, .cnot 2 4, .h 4, .h 4,
    .ccx 0 1 4, .h 4, .h 4, .cnot 1 4, .cnot 0 4]

-- `toffoli_decompose_then_phase_fold`: phase folding never raises the T count.
#guard
  let decomposed := decomposeToffoliGates [.ccx 0 1 2]
  countTGates (pf 3 decomposed) ≤ countTGates decomposed
#guard exactMatrixEquivalent 3 [.ccx 0 1 2]
  (pf 3 (decomposeToffoliGates [.ccx 0 1 2]))

-- `two_toffoli_shared_control`: port the exact Rust regression, 14 T/Tdg gates become 12.
#guard countTGates (decomposeToffoliGates sharedControlToffolis) == 14
#guard countTGates (pf 3 (decomposeToffoliGates sharedControlToffolis)) == 12
#guard exactMatrixEquivalent 3 sharedControlToffolis
  (pf 3 (decomposeToffoliGates sharedControlToffolis))

-- `small_circuit_pipeline`: both folding rounds are monotone in T/Tdg count.
#guard
  let decomposed := decomposeToffoliGates smallToffoliPipeline
  let once := pf 5 decomposed
  let twice := pf 5 once
  countTGates once ≤ countTGates decomposed && countTGates twice ≤ countTGates once

-- Integer phases surrounding decomposed Toffolis remain stable under a second fold.
#guard
  let once := pf 3 (decomposeToffoliGates [.z 2, .ccx 0 1 2, .z 2])
  pf 3 once == once
#guard
  let original : List Gate := [.z 2, .ccx 0 1 2, .z 2]
  let once := pf 3 (decomposeToffoliGates original)
  exactMatrixEquivalent 3 original (pf 3 once)
#guard
  let once := pf 3 (decomposeToffoliGates [.sdg 0, .ccx 0 1 2, .s 1])
  pf 3 once == once
#guard
  let original : List Gate := [.sdg 0, .ccx 0 1 2, .s 1]
  let once := pf 3 (decomposeToffoliGates original)
  exactMatrixEquivalent 3 original (pf 3 once)
#guard
  let once := pf 3 (decomposeToffoliGates [.z 0, .sdg 1, .ccx 0 1 2, .z 2, .s 1])
  pf 3 once == once
#guard
  let original : List Gate := [.z 0, .sdg 1, .ccx 0 1 2, .z 2, .s 1]
  let once := pf 3 (decomposeToffoliGates original)
  exactMatrixEquivalent 3 original (pf 3 once)


/-! ### Basic merges (`two_t_merge_to_s` … `h_prevents_merge`) -/

-- two_t_merge_to_s
#guard pf 1 [t 0, t 0] == [Gate.s 0]

-- t_and_tdg_cancel
#guard pf 1 [t 0, tdg 0] == []

-- four_t_merge_to_s_s
#guard pf 1 [t 0, t 0, t 0, t 0] == [Gate.z 0]

-- eight_t_cancel
#guard pf 1 (List.replicate 8 (t 0)) == []

-- same_parity_across_cnot
#guard pf 2 [t 0, cnot 0 1, t 0] == [Gate.cnot 0 1, Gate.s 0]

-- different_parity_no_merge
#guard pf 2 [t 0, cnot 0 1, t 1] == [Gate.t 0, Gate.cnot 0 1, Gate.t 1]

-- h_prevents_merge
#guard pf 1 [t 0, h 0, t 0] == [Gate.t 0, Gate.h 0, Gate.t 0]


/-! ### Folding across `x`: the complemented parity (`merge_across_x` … `z_x_z_x_identity`) -/

-- merge_across_x
#guard pf 1 [t 0, x 0, t 0] == [Gate.x 0]

-- t_x_tdg_folds_across_x
#guard pf 1 [t 0, x 0, tdg 0] == [Gate.x 0, Gate.sdg 0]

-- rz_folds_across_x
#guard pf 1 [rz (3/10) 0, x 0, rz (7/10) 0] == [Gate.x 0, Gate.rz ((2 : Rat)/5) 0]

-- rz_cancels_across_x
#guard pf 1 [rz (21/50) 0, x 0, rz (21/50) 0] == [Gate.x 0]

-- triple_t_with_two_xs
#guard pf 1 [t 0, x 0, t 0, x 0, t 0] == [Gate.x 0, Gate.x 0, Gate.t 0]

-- x_t_x_t_identity
#guard pf 1 [x 0, t 0, x 0, t 0] == [Gate.x 0, Gate.x 0]

-- mixed_int_and_float_across_x
#guard pf 1 [t 0, x 0, rz (3/10) 0] == [Gate.x 0, Gate.rz ((1 : Rat)/20) 0]

-- z_x_z_x_identity
#guard pf 1 [z 0, x 0, z 0, x 0] == [Gate.x 0, Gate.x 0]


/-! ### Folding across `cnot` (`cnot_target_x_sandwich` … `t_swap_h_swap_t`) -/

-- cnot_target_x_sandwich
#guard pf 2 [t 0, cnot 0 1, x 1, cnot 0 1, t 0] == [Gate.cnot 0 1, Gate.x 1, Gate.cnot 0 1, Gate.s 0]

-- cnot_propagates_negation
#guard pf 2 [t 1, cnot 0 1, x 0, cnot 0 1, t 1] == [Gate.cnot 0 1, Gate.x 0, Gate.cnot 0 1]

-- complement_then_direct_hit
#guard pf 1 [t 0, x 0, t 0, x 0, tdg 0] == [Gate.x 0, Gate.x 0, Gate.tdg 0]

-- h_still_blocks_complement
#guard pf 1 [t 0, h 0, x 0, t 0] == [Gate.t 0, Gate.h 0, Gate.x 0, Gate.t 0]

-- three_qubit_folding
#guard pf 3 [t 0, cnot 0 1, cnot 0 2, t 0, tdg 0] == [Gate.cnot 0 1, Gate.cnot 0 2, Gate.t 0]

-- toffoli_decomposition_fold
#guard pf 3 [h 2, cnot 1 2, tdg 2, cnot 0 2, t 2, cnot 1 2, tdg 2, cnot 0 2, t 1, t 2, h 2, cnot 0 1, t 0, tdg 1, cnot 0 1] == [Gate.h 2,  Gate.cnot 1 2,  Gate.tdg 2,  Gate.cnot 0 2,  Gate.t 2,  Gate.cnot 1 2,  Gate.tdg 2,  Gate.cnot 0 2,  Gate.t 1,  Gate.t 2,  Gate.h 2,  Gate.cnot 0 1,  Gate.t 0,  Gate.tdg 1,  Gate.cnot 0 1]

-- preserves_non_phase_structure
#guard pf 2 [h 0, cnot 0 1, x 1] == [Gate.h 0, Gate.cnot 0 1, Gate.x 1]

-- rz_merge
#guard pf 1 [rz (3/10) 0, rz (7/10) 0] == [Gate.z 0]

-- cross_qubit_merge_via_cnot
#guard pf 2 [t 0, cnot 0 1, t 1, cnot 0 1, cnot 1 0, t 0] == [Gate.t 0, Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 1 0, Gate.s 0]

-- cross_qubit_cancel
#guard pf 2 [t 0, cnot 0 1, tdg 1, cnot 0 1, cnot 1 0, t 0] == [Gate.t 0, Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 1 0]

-- cross_qubit_three_way
#guard pf 3 [cnot 0 2, cnot 1 2, rz (1/2) 2, cnot 1 2, cnot 0 2, cnot 1 2, cnot 0 2, rz (1/2) 2] == [Gate.cnot 0 2,  Gate.cnot 1 2,  Gate.cnot 1 2,  Gate.cnot 0 2,  Gate.cnot 1 2,  Gate.cnot 0 2,  Gate.z 2]

-- cross_qubit_rz_same_parity
#guard pf 2 [cnot 0 1, rz (3/10) 1, cnot 0 1, cnot 1 0, rz (7/10) 0] == [Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 1 0, Gate.z 0]

-- circuit_from_diagram
#guard pf 3 [cnot 0 2, t 2, cnot 1 2, cnot 1 0, tdg 0, cnot 2 0] == [Gate.cnot 0 2,  Gate.t 2,  Gate.cnot 1 2,  Gate.cnot 1 0,  Gate.tdg 0,  Gate.cnot 2 0]

-- cx_t_cx_cx_tdg_cx
#guard pf 2 [cnot 0 1, t 1, cnot 0 1, cnot 1 0, tdg 0, cnot 0 1] == [Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 1 0, Gate.cnot 0 1]

-- t_swap_h_swap_t
#guard pf 2 [t 1, cnot 0 1, cnot 1 0, cnot 0 1, h 1, cnot 0 1, cnot 1 0, cnot 0 1, t 1] == [Gate.cnot 0 1,  Gate.cnot 1 0,  Gate.cnot 0 1,  Gate.h 1,  Gate.cnot 0 1,  Gate.cnot 1 0,  Gate.cnot 0 1,  Gate.s 1]


/-! ### The Clifford phase gates (`z_is_phase_gate` … `sdg_z_is_s`) -/

-- z_is_phase_gate
#guard pf 1 [z 0] == [Gate.z 0]

-- sdg_is_phase_gate
#guard pf 1 [sdg 0] == [Gate.sdg 0]

-- z_z_cancel
#guard pf 1 [z 0, z 0] == []

-- s_sdg_cancel
#guard pf 1 [s 0, sdg 0] == []

-- sdg_s_cancel
#guard pf 1 [sdg 0, s 0] == []

-- s_s_is_z
#guard pf 1 [s 0, s 0] == [Gate.z 0]

-- sdg_sdg_is_z
#guard pf 1 [sdg 0, sdg 0] == [Gate.z 0]

-- three_tdg_is_z_plus_t
#guard pf 1 [tdg 0, tdg 0, tdg 0] == [Gate.z 0, Gate.t 0]

-- z_t_is_z_plus_t
#guard pf 1 [z 0, t 0] == [Gate.z 0, Gate.t 0]

-- z_tdg_is_s_plus_t
#guard pf 1 [z 0, tdg 0] == [Gate.s 0, Gate.t 0]

-- sdg_t_is_tdg
#guard pf 1 [sdg 0, t 0] == [Gate.tdg 0]

-- s_t_folds
#guard pf 1 [s 0, t 0] == [Gate.s 0, Gate.t 0]

-- s_tdg_is_t
#guard pf 1 [s 0, tdg 0] == [Gate.t 0]

-- sdg_tdg_is_z_plus_t
#guard pf 1 [sdg 0, tdg 0] == [Gate.z 0, Gate.t 0]

-- z_s_is_sdg
#guard pf 1 [z 0, s 0] == [Gate.sdg 0]

-- z_sdg_is_s
#guard pf 1 [z 0, sdg 0] == [Gate.s 0]

-- six_t_is_sdg
#guard pf 1 (List.replicate 6 (t 0)) == [Gate.sdg 0]

-- seven_t_is_tdg
#guard pf 1 (List.replicate 7 (t 0)) == [Gate.tdg 0]

-- z_across_cnot_folds
#guard pf 2 [z 0, cnot 0 1, z 0] == [Gate.cnot 0 1]

-- sdg_across_cnot_folds
#guard pf 2 [sdg 0, cnot 0 1, s 0] == [Gate.cnot 0 1]

-- z_h_prevents_merge
#guard pf 1 [z 0, h 0, z 0] == [Gate.z 0, Gate.h 0, Gate.z 0]

-- sdg_h_prevents_merge
#guard pf 1 [sdg 0, h 0, sdg 0] == [Gate.sdg 0, Gate.h 0, Gate.sdg 0]

-- z_x_z_folds_across_x
#guard pf 1 [z 0, x 0, z 0] == [Gate.x 0]

-- cross_qubit_z_merge
#guard pf 2 [z 0, cnot 0 1, cnot 0 1, cnot 1 0] == [Gate.z 0, Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 1 0]

-- cross_qubit_sdg_s_cancel
#guard pf 2 [cnot 0 1, sdg 1, cnot 0 1, cnot 1 0, s 0] == [Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 1 0]

-- multi_qubit_z_sdg_fold
#guard pf 2 [z 0, sdg 1, s 0, t 1] == [Gate.sdg 0, Gate.tdg 1]

-- three_s_is_sdg
#guard pf 1 (List.replicate 3 (s 0)) == [Gate.sdg 0]

-- three_sdg_is_s
#guard pf 1 (List.replicate 3 (sdg 0)) == [Gate.s 0]

-- four_s_cancel
#guard pf 1 (List.replicate 4 (s 0)) == []

-- four_sdg_cancel
#guard pf 1 (List.replicate 4 (sdg 0)) == []

-- z_t_tdg_is_z
#guard pf 1 [z 0, t 0, tdg 0] == [Gate.z 0]

-- all_phase_types_cancel
#guard pf 1 [t 0, tdg 0, s 0, sdg 0, z 0, z 0] == []

-- mixed_z_sdg_cnot_pipeline
#guard pf 3 [z 0, cnot 0 1, sdg 1, cnot 1 2, t 2, cnot 0 2, s 0, tdg 2] == [Gate.cnot 0 1,  Gate.sdg 1,  Gate.cnot 1 2,  Gate.t 2,  Gate.cnot 0 2,  Gate.sdg 0,  Gate.tdg 2]

-- s_z_is_sdg
#guard pf 1 [s 0, z 0] == [Gate.sdg 0]

-- sdg_z_is_s
#guard pf 1 [sdg 0, z 0] == [Gate.s 0]


/-! ### `rz` angles, quarter-turn and not (`rz_pi_folds_to_z` … `sdg_t_rz_quarter_pi_cancel`) -/

-- rz_pi_folds_to_z
#guard pf 1 [rz 1 0] == [Gate.z 0]

-- rz_neg_half_pi_folds_to_sdg
#guard pf 1 [rz (-1/2) 0] == [Gate.sdg 0]

-- rz_three_half_pi_folds_to_sdg
#guard pf 1 [rz (3/2) 0] == [Gate.sdg 0]

-- z_preserves_non_phase_structure
#guard pf 2 [h 0, z 0, cnot 0 1, sdg 1, x 1] == [Gate.h 0, Gate.z 0, Gate.cnot 0 1, Gate.sdg 1, Gate.x 1]

-- t_plus_rz_quarter_pi_is_s
#guard pf 1 [t 0, rz (1/4) 0] == [Gate.s 0]

-- s_plus_rz_half_pi_is_z
#guard pf 1 [s 0, rz (1/2) 0] == [Gate.z 0]

-- rz_pi_plus_tdg
#guard pf 1 [rz 1 0, tdg 0] == [Gate.s 0, Gate.t 0]

-- t_tdg_rz_quarter_pi_is_t
#guard pf 1 [t 0, tdg 0, rz (1/4) 0] == [Gate.t 0]

-- t_plus_rz_irrational_folds
#guard pf 1 [t 0, rz (3/10) 0] == [Gate.rz ((11 : Rat)/20) 0]

-- t_rz_neg_quarter_cancels
#guard pf 1 [t 0, rz (-1/4) 0] == []

-- mixed_int_float_across_cnot
#guard pf 2 [t 0, cnot 1 0, cnot 0 1, rz (3/10) 1] == [Gate.cnot 1 0, Gate.cnot 0 1, Gate.rz ((11 : Rat)/20) 1]

-- s_plus_two_rz_quarter_pi
#guard pf 1 [s 0, rz (1/4) 0, rz (1/4) 0] == [Gate.z 0]

-- t_rz_irr_rz_opposite
#guard pf 1 [t 0, rz (3/10) 0, rz (-3/10) 0] == [Gate.t 0]

-- s_t_rz_pi_combine
#guard pf 1 [s 0, t 0, rz 1 0] == [Gate.tdg 0]

-- sdg_t_rz_quarter_pi_cancel
#guard pf 1 [sdg 0, t 0, rz (1/4) 0] == []


/-! ### Measurement and reset

The nonlinear executable matches Rust here: a measurement preserves the value it measures,
and a `reset` pins a wire to the field constants `0` or (after X) `1`. Rotations on either
known computational-basis constant are observationally irrelevant and are removed. -/

-- measure_t_t
#guard pf 1 [t 0, measure 0 0, t 0] == [Gate.measure 0 0, Gate.s 0]

-- reset_then_phase
#guard pf 1 [t 0, reset 0, t 0] == [Gate.t 0, Gate.reset 0]

-- measure_t_tdg
#guard pf 1 [t 0, measure 0 0, tdg 0] == [Gate.measure 0 0]

-- measure_other_qubit
#guard pf 2 [t 0, measure 1 0, t 0] == [Gate.measure 1 0, Gate.s 0]

-- reset_other_qubit
#guard pf 2 [t 0, reset 1, t 0] == [Gate.reset 1, Gate.s 0]

-- measure_no_rotations
#guard pf 2 [h 0, cnot 0 1, measure 0 0, measure 1 1] == [Gate.h 0, Gate.cnot 0 1, Gate.measure 0 0, Gate.measure 1 1]

-- measure_both_sides
#guard pf 1 [t 0, t 0, measure 0 0, t 0, t 0] == [Gate.measure 0 0, Gate.z 0]

-- reset_then_phases
#guard pf 1 [reset 0, t 0, s 0, z 0, rz (123/1000) 0] == [Gate.reset 0]

-- reset_then_x_then_t
#guard pf 1 [reset 0, x 0, t 0] == [Gate.reset 0, Gate.x 0]

-- reset_zero_through_cnot
#guard pf 2 [reset 1, cnot 0 1, t 0, t 1] == [Gate.reset 1, Gate.cnot 0 1, Gate.s 1]

-- hadamard_after_reset
#guard pf 1 [reset 0, h 0, t 0] == [Gate.reset 0, Gate.h 0, Gate.t 0]


/-! ### Diagonal two- and three-qubit gates are transparent (`cz_preserved` … `cz_no_hadamard_hiding`) -/

-- cz_preserved
#guard pf 2 [cz 0 1] == [Gate.cz 0 1]

-- t_pair_through_cz_q0
#guard pf 2 [t 0, cz 0 1, t 0] == [Gate.cz 0 1, Gate.s 0]

-- t_pair_through_cz_q1
#guard pf 2 [t 1, cz 0 1, t 1] == [Gate.cz 0 1, Gate.s 1]

-- t_pair_through_ccz_q0
#guard pf 3 [t 0, ccz 0 1 2, t 0] == [Gate.ccz 0 1 2, Gate.s 0]

-- t_pair_through_ccz_q1
#guard pf 3 [t 1, ccz 0 1 2, t 1] == [Gate.ccz 0 1 2, Gate.s 1]

-- t_pair_through_ccz_q2
#guard pf 3 [t 2, ccz 0 1 2, t 2] == [Gate.ccz 0 1 2, Gate.s 2]

-- opposite_through_cz_q0
#guard pf 2 [t 0, cz 0 1, tdg 0] == [Gate.cz 0 1]

-- opposite_through_cz_q1
#guard pf 2 [t 1, cz 0 1, tdg 1] == [Gate.cz 0 1]

-- arbitrary_rz_through_cz
#guard pf 2 [rz (37/100) 1, cz 1 0, rz (-17/100) 1] == [Gate.cz 1 0, Gate.rz ((1 : Rat)/5) 1]

-- cz_no_hadamard_hiding
#guard pf 2 [t 0, cz 0 1, h 0, t 0] == [Gate.t 0, Gate.cz 0 1, Gate.h 0, Gate.t 0]

-- phases_on_both_wires_fold_independently_across_cz_chain
#guard pf 2 [t 0, tdg 1, cz 0 1, cz 1 0, t 0, t 1] ==
  [Gate.cz 0 1, Gate.cz 1 0, Gate.s 0]

-- cz_and_cnot_control_preserve_phase_parity_together
#guard pf 3 [t 0, cz 0 1, cnot 0 2, cz 2 1, t 0] ==
  [Gate.cz 0 1, Gate.cnot 0 2, Gate.cz 2 1, Gate.s 0]

-- cz_does_not_mask_cnot_target_parity_change
#guard pf 2 [t 1, cz 0 1, cnot 0 1, cz 1 0, t 1] ==
  [Gate.t 1, Gate.cz 0 1, Gate.cnot 0 1, Gate.cz 1 0, Gate.t 1]

-- phase_fold_preserves_cz_count_and_operand_order
#guard pf 3 [t 2, cz 2 0, cz 1 2, t 2] ==
  [Gate.cz 2 0, Gate.cz 1 2, Gate.s 2]

/-! ### The count assertions Rust states as inequalities -/

-- `cnot_propagates_negation_from_control`: at most one `t` survives.
#guard countTGates (pf 2 [t 1, cnot 0 1, x 0, cnot 0 1, t 1]) ≤ 1

-- `toffoli_decomposition_fold`: never more phase gates than we started with.
#guard countPhaseGates (pf 3 [h 2, cnot 1 2, tdg 2, cnot 0 2, t 2, cnot 1 2, tdg 2, cnot 0 2,
    t 1, t 2, h 2, cnot 0 1, t 0, tdg 1, cnot 0 1]) ≤
  countPhaseGates [h 2, cnot 1 2, tdg 2, cnot 0 2, t 2, cnot 1 2, tdg 2, cnot 0 2,
    t 1, t 2, h 2, cnot 0 1, t 0, tdg 1, cnot 0 1]

-- `mixed_z_sdg_cnot_pipeline`: likewise on a mixed circuit.
#guard countPhaseGates (pf 3 [z 0, cnot 0 1, sdg 1, cnot 1 2, t 2, cnot 0 2, s 0, tdg 2]) ≤
  countPhaseGates [z 0, cnot 0 1, sdg 1, cnot 1 2, t 2, cnot 0 2, s 0, tdg 2]

-- `z_preserves_non_phase_structure`: one `h`, one `x`, one `cnot`, whatever the phases do.
#guard (pf 2 [h 0, z 0, cnot 0 1, sdg 1, x 1]).filter (fun g => (rotAngle g).isNone) ==
  [h 0, cnot 0 1, x 1]

/-! ### The pass is idempotent on its own output -/

#guard pf 1 (pf 1 [t 0, t 0, t 0]) == pf 1 [t 0, t 0, t 0]
#guard pf 2 (pf 2 [t 0, cnot 0 1, t 0, x 1, t 1]) == pf 2 [t 0, cnot 0 1, t 0, x 1, t 1]
#guard pf 3 (pf 3 [t 0, cnot 0 2, tdg 2, h 1, t 1, cnot 1 2, s 2]) ==
  pf 3 [t 0, cnot 0 2, tdg 2, h 1, t 1, cnot 1 2, s 2]

end TzapLean
