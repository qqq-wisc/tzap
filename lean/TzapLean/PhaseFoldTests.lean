import TzapLean.PhaseFoldRand

/-!
# `PhaseFoldRand`: the Rust test suite, ported

Every behavioural test from `src/phase_fold_rand.rs` that does not depend on a pass this
development has not ported (the Toffoli decomposition). Rust checks gate *counts* plus a
numerical `circuits_equiv`; here the equivalence is already a theorem for every input
(`phaseFoldGates_correct`), so these `#guard`s pin down the stronger fact — the exact gate
list the pass produces — and the counts follow.

Angles are rationals in units of `π`, so Rust's radian constants become the `π`-fraction with
the same classification: `0.3` (not a quarter turn) becomes `3/10`, `PI/4` becomes `1/4`.

The draws here are a fixed splitmix stream, so the results are reproducible. The optimizer
itself never uses one: `phaseFoldIO` draws a uniform `Sample`, which is the space
`PhaseFoldRand.correct` bounds a measure over, and nothing relates a splitmix stream to it.
-/

namespace TzapLean

open Gate

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

/-- Phase folding with those draws. -/
def pf (n : Nat) (gs : List Gate) : List Gate := phaseFoldGates 63 testWords n gs

/-- Rust's `count_phase_gates`. -/
def countPhaseGates (gs : List Gate) : Nat := (gs.filter fun g => (rotAngle g).isSome).length

/-- Rust's `count_t_gates`. -/
def countTGates (gs : List Gate) : Nat :=
  (gs.filter fun g => match g with | .t _ | .tdg _ => true | _ => false).length


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

Rust folds *through* both: a measurement preserves the value it measures, and a `reset`
pins a wire's parity to the constant `0`. This port stops the lookahead at either gate, so
it merges strictly less — never wrongly. Four of these agree with Rust exactly; the rest
keep gates Rust removes, and the expectations below record that difference rather than
paper over it. Lifting it needs the merge lemma restated on channels (a diagonal operator
commutes with the measurement projectors) instead of on unitaries. -/

-- measure_t_t
#guard pf 1 [t 0, measure 0 0, t 0] == [Gate.t 0, Gate.measure 0 0, Gate.t 0]

-- reset_then_phase
#guard pf 1 [t 0, reset 0, t 0] == [Gate.t 0, Gate.reset 0, Gate.t 0]

-- measure_t_tdg
#guard pf 1 [t 0, measure 0 0, tdg 0] == [Gate.t 0, Gate.measure 0 0, Gate.tdg 0]

-- measure_other_qubit
#guard pf 2 [t 0, measure 1 0, t 0] == [Gate.t 0, Gate.measure 1 0, Gate.t 0]

-- reset_other_qubit
#guard pf 2 [t 0, reset 1, t 0] == [Gate.t 0, Gate.reset 1, Gate.t 0]

-- measure_no_rotations
#guard pf 2 [h 0, cnot 0 1, measure 0 0, measure 1 1] == [Gate.h 0, Gate.cnot 0 1, Gate.measure 0 0, Gate.measure 1 1]

-- measure_both_sides
#guard pf 1 [t 0, t 0, measure 0 0, t 0, t 0] == [Gate.s 0, Gate.measure 0 0, Gate.s 0]

-- reset_then_phases
#guard pf 1 [reset 0, t 0, s 0, z 0, rz (123/1000) 0] == [Gate.reset 0, Gate.rz ((1873 : Rat)/1000) 0]

-- reset_then_x_then_t
#guard pf 1 [reset 0, x 0, t 0] == [Gate.reset 0, Gate.x 0, Gate.t 0]

-- reset_zero_through_cnot
#guard pf 2 [reset 1, cnot 0 1, t 0, t 1] == [Gate.reset 1, Gate.cnot 0 1, Gate.t 0, Gate.t 1]

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
