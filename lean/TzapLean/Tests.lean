import TzapLean.Cancel
import TzapLean.CnotMinProof

/-!
# Ported Tests

The Rust test suites of `src/circuit.rs`, `src/pass.rs`, and `src/cancel.rs`, as `#guard`
checks — they run at compile time, so `lake build` is the test run.

One systematic difference: every Rust test also asserts
`circuits_equiv(&c, &r, 1e-10)`, a numerical check that the pass preserved the unitary on
that one input. Those assertions have no counterpart here because
`TzapLean.cancelGates_correct` proves the equivalence for *every* well-formed input, which
is strictly stronger. What remains worth testing is the *structure* of the output — which
gates survive, and how many — and that is what these checks cover.

`rz` angles are in units of `π` (see `TzapLean.Gate`), so the Rust `Gate::rz(PI/4.0, q)`
is `Gate.rz (1/4) q`. Where a Rust test used a non-`π` angle such as `0.31`, the port keeps
the literal: those tests only care that the gate is a non-Clifford diagonal, which
`0.31 · π` equally is.
-/

namespace TzapLean

/-- Number of gates of a given kind, for the count-based assertions. -/
def countKind (p : Gate → Bool) (c : Circuit) : Nat := c.gates.countP p

def isHGate : Gate → Bool | .h _ => true | _ => false
def isCnotGate : Gate → Bool | .cnot .. => true | _ => false
def isCzGate : Gate → Bool | .cz .. => true | _ => false

/-- `count_h` from the Rust tests. -/
def countH (c : Circuit) : Nat := countKind isHGate c

/-- The pass under test, as a plain function. -/
def runCancel (c : Circuit) : Circuit := cancelGatesCircuit c

/-! ## `src/cancel.rs` -/

-- hh_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0])).gates.length = 0

-- xx_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.x 0, Gate.x 0])).gates.length = 0

-- cnot_cancel
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.cnot 0 1])).gates.length = 0

-- hh_cancel_different_qubit
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 3, Gate.h 3])).gates.length = 0

-- hh_cancel_skips_unrelated_gate
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.t 1, Gate.h 0])).gates.length = 1
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.t 1, Gate.h 0])).gates[0]? = some (Gate.t 1)

-- hh_cancel_blocked_by_same_qubit
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.t 0, Gate.h 0])).gates.length = 3

-- hh_cancel_multiple_pairs
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0, Gate.h 0, Gate.h 0])).gates.length = 0

-- hh_cancel_parallel_qubits
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.h 1, Gate.h 0, Gate.h 1])).gates.length = 0

-- xx_cancel_different_qubit
#guard (runCancel (Circuit.ofGates 3 0 [Gate.x 2, Gate.x 2])).gates.length = 0

-- xx_cancel_skips_unrelated_gate
#guard (runCancel (Circuit.ofGates 2 0 [Gate.x 0, Gate.h 1, Gate.x 0])).gates.length = 1
#guard (runCancel (Circuit.ofGates 2 0 [Gate.x 0, Gate.h 1, Gate.x 0])).gates[0]? = some (Gate.h 1)

-- xx_cancel_blocked_by_same_qubit
#guard (runCancel (Circuit.ofGates 1 0 [Gate.x 0, Gate.z 0, Gate.x 0])).gates.length = 3

-- xx_cancel_multiple_pairs
#guard (runCancel (Circuit.ofGates 1 0 [Gate.x 0, Gate.x 0, Gate.x 0, Gate.x 0])).gates.length = 0

-- cnot_cancel_different_qubits
#guard (runCancel (Circuit.ofGates 5 0 [Gate.cnot 3 4, Gate.cnot 3 4])).gates.length = 0

-- cnot_cancel_skips_unrelated_gate
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cnot 0 1, Gate.t 2, Gate.cnot 0 1])).gates.length = 1
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cnot 0 1, Gate.t 2, Gate.cnot 0 1])).gates[0]? = some (Gate.t 2)

-- cnot_cancels_through_diagonal_on_control
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.t 0, Gate.cnot 0 1])).gates.length = 1

-- cnots_cancel_through_x_on_target
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.x 1, Gate.cnot 0 1])).gates.length = 1

-- cnots_cancel_through_sharing_control_cnot
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cnot 0 1, Gate.cnot 0 2, Gate.cnot 0 1])).gates.length = 1

-- cnots_cancel_through_sharing_target_cnot
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cnot 0 2, Gate.cnot 1 2, Gate.cnot 0 2])).gates.length = 1

-- cnots_blocked_by_cnot_with_swapped_endpoints
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.cnot 1 0, Gate.cnot 0 1])).gates.length = 3

-- cnots_blocked_by_diagonal_on_target
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.t 1, Gate.cnot 0 1])).gates.length = 3

-- cnots_blocked_by_x_on_control
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.x 0, Gate.cnot 0 1])).gates.length = 3

-- cnot_cancel_blocked_by_hadamard_on_control
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.h 0, Gate.cnot 0 1])).gates.length = 3

-- cnot_cancel_blocked_by_gate_on_target
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.t 1, Gate.cnot 0 1])).gates.length = 3

-- cnot_cancel_no_match_different_direction
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.cnot 1 0])).gates.length = 2

-- cnot_cancel_multiple_pairs
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 0 1])).gates.length = 0

-- no_match_preserves_circuit
#guard (runCancel (Circuit.ofGates 2 0 [Gate.t 0, Gate.s 1, Gate.cnot 0 1])).gates.length = 3

-- single_h_preserved
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0])).gates.length = 1
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0])).gates[0]? = some (Gate.h 0)

-- zz_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.z 0, Gate.z 0])).gates.length = 0

-- zz_cancel_skips_unrelated_gate
#guard (runCancel (Circuit.ofGates 2 0 [Gate.z 0, Gate.t 1, Gate.z 0])).gates.length = 1
#guard (runCancel (Circuit.ofGates 2 0 [Gate.z 0, Gate.t 1, Gate.z 0])).gates[0]? = some (Gate.t 1)

-- zz_cancel_blocked_by_same_qubit
#guard (runCancel (Circuit.ofGates 1 0 [Gate.z 0, Gate.h 0, Gate.z 0])).gates.length = 3

-- ccx_cancel
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccx 0 1 2, Gate.ccx 0 1 2])).gates.length = 0

-- ccz_cancel_is_symmetric_in_all_operands
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccz 0 1 2, Gate.ccz 2 0 1])).hasToffoli = false

-- ccx_cancel_skips_unrelated_gate
#guard (runCancel (Circuit.ofGates 4 0 [Gate.ccx 0 1 2, Gate.t 3, Gate.ccx 0 1 2])).gates.length = 1
#guard (runCancel (Circuit.ofGates 4 0 [Gate.ccx 0 1 2, Gate.t 3, Gate.ccx 0 1 2])).gates[0]? = some (Gate.t 3)

-- ccx_cancel_blocked_by_control1
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccx 0 1 2, Gate.h 0, Gate.ccx 0 1 2])).gates.length = 3

-- ccx_cancel_blocked_by_control2
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccx 0 1 2, Gate.h 1, Gate.ccx 0 1 2])).gates.length = 3

-- ccx_cancel_blocked_by_target
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccx 0 1 2, Gate.h 2, Gate.ccx 0 1 2])).gates.length = 3

-- ccx_no_cancel_different_controls
#guard (runCancel (Circuit.ofGates 4 0 [Gate.ccx 0 1 3, Gate.ccx 0 2 3])).gates.length = 2

-- ccx_no_cancel_different_target
#guard (runCancel (Circuit.ofGates 4 0 [Gate.ccx 0 1 2, Gate.ccx 0 1 3])).gates.length = 2

-- ccx_cancel_multiple_pairs
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccx 0 1 2, Gate.ccx 0 1 2, Gate.ccx 0 1 2, Gate.ccx 0 1 2])).gates.length = 0

-- cascade_nested_h
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0, Gate.h 0, Gate.h 0])).gates.length = 0

-- cascade_six_h
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0, Gate.h 0, Gate.h 0, Gate.h 0, Gate.h 0])).gates.length = 0

-- cascade_odd_count_leaves_one
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0, Gate.h 0, Gate.h 0, Gate.h 0])).gates.length = 1

-- cascade_nested_cnot
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 0 1])).gates.length = 0

-- cascade_mixed_self_inverse
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.x 0, Gate.h 0, Gate.x 0])).gates.length = 2

-- cascade_different_qubits_interleaved
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.h 1, Gate.h 1, Gate.h 0])).gates.length = 0

-- cascade_deep_nesting
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.x 0, Gate.z 0, Gate.z 0, Gate.x 0, Gate.h 0])).gates.length = 0

-- cascade_deep_nesting_with_residual
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.h 0, Gate.x 0, Gate.x 0, Gate.h 0, Gate.t 0])).gates.length = 2
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.h 0, Gate.x 0, Gate.x 0, Gate.h 0, Gate.t 0])).gates[0]? = some (Gate.t 0)
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.h 0, Gate.x 0, Gate.x 0, Gate.h 0, Gate.t 0])).gates[1]? = some (Gate.t 0)

-- cascade_cnot_nested_in_h
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.cnot 0 1, Gate.cnot 0 1, Gate.h 0])).gates.length = 0

-- commute_h_past_many_unrelated
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 0, Gate.t 1, Gate.s 2, Gate.tdg 3, Gate.h 0])).gates.length = 3

-- commute_cnot_past_unrelated_qubits
#guard (runCancel (Circuit.ofGates 4 0 [Gate.cnot 0 1, Gate.h 2, Gate.t 3, Gate.cnot 0 1])).gates.length = 2

-- commute_ccx_past_unrelated_qubits
#guard (runCancel (Circuit.ofGates 5 0 [Gate.ccx 0 1 2, Gate.h 3, Gate.t 4, Gate.ccx 0 1 2])).gates.length = 2

-- t_t_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.t 0])).gates.length = 2

-- tdg_tdg_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.tdg 0, Gate.tdg 0])).gates.length = 2

-- s_s_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.s 0, Gate.s 0])).gates.length = 2

-- sdg_sdg_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.sdg 0, Gate.sdg 0])).gates.length = 2

-- t_tdg_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.tdg 0])).gates.length = 2

-- s_sdg_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.s 0, Gate.sdg 0])).gates.length = 2

-- rz_rz_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.rz (1/4) 0, Gate.rz (1/4) 0])).gates.length = 2

-- independent_cancellations_on_many_qubits
#guard (runCancel (Circuit.ofGates 5 0 [Gate.h 0, Gate.h 1, Gate.h 2, Gate.h 3, Gate.h 4, Gate.h 4, Gate.h 3, Gate.h 2, Gate.h 1, Gate.h 0])).gates.length = 0

-- partial_cancellation_mixed_qubits
#guard (runCancel (Circuit.ofGates 3 0 [Gate.h 0, Gate.h 1, Gate.x 2, Gate.t 1, Gate.h 0, Gate.h 1, Gate.x 2])).gates.length = 3

-- interleaved_cancel_different_gate_types
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.x 1, Gate.h 0, Gate.x 1])).gates.length = 0

-- cnot_blocks_h_on_shared_qubit
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.cnot 0 1, Gate.h 0])).gates.length = 3

-- h_not_blocked_by_cnot_on_other_qubits
#guard (runCancel (Circuit.ofGates 3 0 [Gate.h 0, Gate.cnot 1 2, Gate.h 0])).gates.length = 1

-- surviving_gates_preserve_order
#guard (runCancel (Circuit.ofGates 3 0 [Gate.t 0, Gate.h 1, Gate.h 1, Gate.s 0, Gate.x 2, Gate.x 2, Gate.tdg 0])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 0 [Gate.t 0, Gate.h 1, Gate.h 1, Gate.s 0, Gate.x 2, Gate.x 2, Gate.tdg 0])).gates[0]? = some (Gate.t 0)
#guard (runCancel (Circuit.ofGates 3 0 [Gate.t 0, Gate.h 1, Gate.h 1, Gate.s 0, Gate.x 2, Gate.x 2, Gate.tdg 0])).gates[1]? = some (Gate.s 0)
#guard (runCancel (Circuit.ofGates 3 0 [Gate.t 0, Gate.h 1, Gate.h 1, Gate.s 0, Gate.x 2, Gate.x 2, Gate.tdg 0])).gates[2]? = some (Gate.tdg 0)

-- order_after_cascade
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.h 0, Gate.h 0, Gate.s 0])).gates.length = 2
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.h 0, Gate.h 0, Gate.s 0])).gates[0]? = some (Gate.t 0)
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.h 0, Gate.h 0, Gate.s 0])).gates[1]? = some (Gate.s 0)

-- alternating_cancel_no_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0, Gate.t 0, Gate.h 0, Gate.h 0, Gate.t 0])).gates.length = 2
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0, Gate.t 0, Gate.h 0, Gate.h 0, Gate.t 0])).gates[0]? = some (Gate.t 0)
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0, Gate.t 0, Gate.h 0, Gate.h 0, Gate.t 0])).gates[1]? = some (Gate.t 0)

-- deeply_nested_cascade_8_layers
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.x 0, Gate.z 0, Gate.h 0, Gate.h 0, Gate.z 0, Gate.x 0, Gate.h 0])).gates.length = 0

-- has_toffoli_set_when_ccx_survives
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 3, Gate.h 3, Gate.ccx 0 1 2])).gates.length = 1
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 3, Gate.h 3, Gate.ccx 0 1 2])).hasToffoli = true

-- has_toffoli_cleared_when_ccx_cancelled
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccx 0 1 2, Gate.ccx 0 1 2])).gates.length = 0
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccx 0 1 2, Gate.ccx 0 1 2])).hasToffoli = false

-- cnot_and_ccx_block_each_other_on_shared_qubit
#guard (runCancel (Circuit.ofGates 4 0 [Gate.cnot 0 1, Gate.ccx 1 2 3, Gate.cnot 0 1])).gates.length = 3

-- cnot_and_ccx_no_shared_qubit
#guard (runCancel (Circuit.ofGates 5 0 [Gate.cnot 0 1, Gate.ccx 2 3 4, Gate.cnot 0 1])).gates.length = 1

-- multiple_cnot_pairs_different_qubit_pairs
#guard (runCancel (Circuit.ofGates 4 0 [Gate.cnot 0 1, Gate.cnot 2 3, Gate.cnot 2 3, Gate.cnot 0 1])).gates.length = 0

-- single_gate_each_type_preserved
#guard (runCancel (Circuit.ofGates 5 0 [Gate.h 0, Gate.x 1, Gate.z 2, Gate.cnot 3 4])).gates.length = 4

-- blocker_is_latest_not_earliest
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.t 0, Gate.x 0, Gate.h 0])).gates.length = 4

-- blocker_matches_when_latest_is_same
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.t 0, Gate.h 0, Gate.h 0])).gates.length = 2
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.t 0, Gate.h 0, Gate.h 0])).gates[0]? = some (Gate.h 0)
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.t 0, Gate.h 0, Gate.h 0])).gates[1]? = some (Gate.t 0)

-- cascade_partial_then_full
#guard (runCancel (Circuit.ofGates 2 0 [Gate.t 1, Gate.h 0, Gate.x 1, Gate.h 0, Gate.x 1, Gate.t 1])).gates.length = 2
#guard (runCancel (Circuit.ofGates 2 0 [Gate.t 1, Gate.h 0, Gate.x 1, Gate.h 0, Gate.x 1, Gate.t 1])).gates[0]? = some (Gate.t 1)
#guard (runCancel (Circuit.ofGates 2 0 [Gate.t 1, Gate.h 0, Gate.x 1, Gate.h 0, Gate.x 1, Gate.t 1])).gates[1]? = some (Gate.t 1)

-- hh_blocked_by_reset_on_same_qubit
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.reset 0, Gate.h 0])).gates.length = 3

-- hh_allowed_past_reset_on_other_qubit
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.reset 1, Gate.h 0])).gates.length = 1
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.reset 1, Gate.h 0])).gates[0]? = some (Gate.reset 1)

-- cnot_blocked_by_reset_on_target
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.reset 1, Gate.cnot 0 1])).gates.length = 3

-- measure_reset_pair_does_not_cancel
#guard (runCancel (Circuit.ofGates 1 0 [Gate.reset 0, Gate.reset 0])).gates.length = 2

-- hxh_becomes_z
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.x 0, Gate.h 0])).gates.length = 1

-- hzh_becomes_x
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.z 0, Gate.h 0])).gates.length = 1

-- hadamard_run_summing_to_identity_collapses
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.sdg 0, Gate.s 0, Gate.h 0])).gates.length = 0

-- cnot_inside_hadamards_blocks_reduction
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.cnot 0 1, Gate.h 0])).gates.length = 3

-- lone_hadamard_is_kept
#guard (runCancel (Circuit.ofGates 1 0 [Gate.t 0, Gate.h 0, Gate.t 0])).gates.length = 3

-- reduction_exposes_pair_cancellation
#guard (runCancel (Circuit.ofGates 1 0 [Gate.z 0, Gate.h 0, Gate.x 0, Gate.h 0])).gates.length = 0

-- cz_pair_cancels_through_diagonals_on_both_operands
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cz 0 1, Gate.t 0, Gate.sdg 1, Gate.rz (0.31) 0, Gate.z 1, Gate.cz 1 0])).gates.length = 4

-- reset_on_operand_blocks_cz_cancellation
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cz 0 1, Gate.reset 1, Gate.cz 0 1])).gates.length = 3

-- cnot_pair_blocked_by_cz_on_target
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cnot 2 1, Gate.cz 0 1, Gate.cnot 2 1])).gates.length = 3


/-! ### Tests needing richer assertions than the mechanical port -/

-- cnots_cancel_through_diagonal_run_on_control
#guard countKind isCnotGate
  (runCancel (Circuit.ofGates 2 0
    [Gate.cnot 0 1, Gate.t 0, Gate.s 0, Gate.tdg 0, Gate.cnot 0 1])) = 0

-- cnots_cancel_through_long_mixed_run
#guard countKind isCnotGate
  (runCancel (Circuit.ofGates 4 0
    [Gate.cnot 0 1, Gate.t 0, Gate.x 1, Gate.h 3, Gate.cnot 0 2, Gate.rz (31/100) 0,
     Gate.cnot 0 1])) = 1

-- cnot_lookahead_idempotent
#guard (runCancel (runCancel (Circuit.ofGates 3 0
    [Gate.cnot 0 1, Gate.t 0, Gate.cnot 0 2, Gate.cnot 0 1]))).gates.length
  = (runCancel (Circuit.ofGates 3 0
    [Gate.cnot 0 1, Gate.t 0, Gate.cnot 0 2, Gate.cnot 0 1])).gates.length

-- empty_circuit
#guard (runCancel (Circuit.new 2)).gates.length = 0

-- ccz_cancellation_is_blocked_on_each_operand
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.ccz 0 1 2, Gate.h 0, Gate.ccz 2 0 1])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccz 0 1 2, Gate.h 0, Gate.ccz 2 0 1])).hasCcz = true
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.ccz 0 1 2, Gate.h 1, Gate.ccz 2 0 1])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccz 0 1 2, Gate.h 1, Gate.ccz 2 0 1])).hasCcz = true
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.ccz 0 1 2, Gate.h 2, Gate.ccz 2 0 1])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 0 [Gate.ccz 0 1 2, Gate.h 2, Gate.ccz 2 0 1])).hasCcz = true

-- ccz_cancels_across_disjoint_gate
#guard (runCancel (Circuit.ofGates 4 0
  [Gate.ccz 0 1 2, Gate.t 3, Gate.ccz 2 1 0])).gates = [Gate.t 3]
#guard (runCancel (Circuit.ofGates 4 0 [Gate.ccz 0 1 2, Gate.t 3, Gate.ccz 2 1 0])).hasCcz = false

-- cnot_blocked_by_gate_on_either_qubit
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cnot 0 1, Gate.h 0, Gate.cnot 0 1])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cnot 0 1, Gate.h 1, Gate.cnot 0 1])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cnot 0 1, Gate.h 2, Gate.cnot 0 1])).gates.length = 1

-- ccx_blocks_h_on_any_of_three_qubits
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 0, Gate.ccx 0 1 2, Gate.h 0])).gates.length = 3
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 1, Gate.ccx 0 1 2, Gate.h 1])).gates.length = 3
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 2, Gate.ccx 0 1 2, Gate.h 2])).gates.length = 3
#guard (runCancel (Circuit.ofGates 4 0 [Gate.h 3, Gate.ccx 0 1 2, Gate.h 3])).gates.length = 1

-- many_pairs_single_qubit
#guard (runCancel (Circuit.ofGates 1 0 (List.replicate 200 (Gate.h 0)))).gates.length = 0

-- many_pairs_odd_leaves_one
#guard (runCancel (Circuit.ofGates 1 0 (List.replicate 201 (Gate.h 0)))).gates.length = 1

-- hh_blocked_by_measure_on_same_qubit
#guard (runCancel (Circuit.ofGates 1 1
  [Gate.h 0, Gate.measure 0 0, Gate.h 0])).gates.length = 3
#guard (runCancel (Circuit.ofGates 1 1
  [Gate.h 0, Gate.measure 0 0, Gate.h 0])).hasMeasurement = true
#guard (runCancel (Circuit.ofGates 1 1 [Gate.h 0, Gate.measure 0 0, Gate.h 0])).numCbits = 1

-- hh_blocked_by_reset_on_same_qubit
#guard (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.reset 0, Gate.h 0])).gates.length = 3

-- hh_allowed_past_measure_on_other_qubit
#guard (runCancel (Circuit.ofGates 2 1
  [Gate.h 0, Gate.measure 1 0, Gate.h 0])).gates = [Gate.measure 1 0]

-- hh_allowed_past_reset_on_other_qubit
#guard (runCancel (Circuit.ofGates 2 0 [Gate.h 0, Gate.reset 1, Gate.h 0])).gates = [Gate.reset 1]

-- xx_blocked_by_measure
#guard (runCancel (Circuit.ofGates 1 1
  [Gate.x 0, Gate.measure 0 0, Gate.x 0])).gates.length = 3

-- cnot_blocked_by_measure_on_control
#guard (runCancel (Circuit.ofGates 2 1
  [Gate.cnot 0 1, Gate.measure 0 0, Gate.cnot 0 1])).gates.length = 3

-- cnot_blocked_by_measure_on_target
#guard (runCancel (Circuit.ofGates 2 1
  [Gate.cnot 0 1, Gate.measure 1 0, Gate.cnot 0 1])).gates.length = 3

-- cnot_allowed_past_measure_on_disjoint_qubit
#guard (runCancel (Circuit.ofGates 3 1
  [Gate.cnot 0 1, Gate.measure 2 0, Gate.cnot 0 1])).gates = [Gate.measure 2 0]

-- cnot_blocked_by_reset_on_target
#guard (runCancel (Circuit.ofGates 2 0
  [Gate.cnot 0 1, Gate.reset 1, Gate.cnot 0 1])).gates.length = 3

-- ccx_blocked_by_measure_on_any_qubit
#guard (runCancel (Circuit.ofGates 3 1
  [Gate.ccx 0 1 2, Gate.measure 0 0, Gate.ccx 0 1 2])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 1
  [Gate.ccx 0 1 2, Gate.measure 1 0, Gate.ccx 0 1 2])).gates.length = 3
#guard (runCancel (Circuit.ofGates 3 1
  [Gate.ccx 0 1 2, Gate.measure 2 0, Gate.ccx 0 1 2])).gates.length = 3

-- measure_reset_pair_does_not_cancel
#guard (runCancel (Circuit.ofGates 1 1
  [Gate.measure 0 0, Gate.measure 0 0])).gates.length = 2
#guard (runCancel (Circuit.ofGates 1 0 [Gate.reset 0, Gate.reset 0])).gates.length = 2

-- num_cbits_preserved
#guard (runCancel (Circuit.ofGates 1 3 [Gate.h 0, Gate.h 0])).gates.length = 0
#guard (runCancel (Circuit.ofGates 1 3 [Gate.h 0, Gate.h 0])).numCbits = 3

-- cascade_depth_10
#guard (runCancel (Circuit.ofGates 1 0
  [Gate.h 0, Gate.x 0, Gate.z 0, Gate.z 0, Gate.x 0, Gate.h 0])).gates.length = 0
#guard (runCancel (Circuit.ofGates 1 0
  [Gate.h 0, Gate.x 0, Gate.z 0, Gate.x 0, Gate.h 0,
   Gate.h 0, Gate.x 0, Gate.z 0, Gate.x 0, Gate.h 0])).gates.length = 0

-- hsh_drops_one_hadamard
#guard countH (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.s 0, Gate.h 0])) = 1

-- hsh_cascade_eliminates_all_hadamards
#guard countH (runCancel (Circuit.ofGates 1 0
  [Gate.h 0, Gate.s 0, Gate.h 0, Gate.s 0, Gate.h 0])) = 0

-- hxh_reduced_despite_other_qubit_gate_between
#guard countH (runCancel (Circuit.ofGates 3 0
  [Gate.h 0, Gate.x 0, Gate.cnot 1 2, Gate.h 0])) = 0

-- hadamard_run_with_t_is_not_reducible
#guard countH (runCancel (Circuit.ofGates 1 0 [Gate.h 0, Gate.t 0, Gate.h 0])) = 2

-- measure_between_hadamards_blocks_reduction
#guard (runCancel (Circuit.ofGates 1 1
  [Gate.h 0, Gate.measure 0 0, Gate.h 0])).gates[1]? = some (Gate.measure 0 0)

-- pair_cancellation_exposes_reduction
#guard countH (runCancel (Circuit.ofGates 1 0
  [Gate.t 0, Gate.h 0, Gate.x 0, Gate.x 0, Gate.h 0, Gate.t 0])) = 0

-- adjacent_cz_pair_cancels
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cz 0 1, Gate.cz 0 1])).gates = []

-- reversed_operand_cz_pair_cancels
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cz 0 2, Gate.cz 2 0])).gates = []

-- cz_pair_lookahead_commutation_table
#guard countKind isCzGate (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.cz 1 2, Gate.cz 1 0])) = 1
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cz 0 1, Gate.cz 1 2, Gate.cz 1 0])).gates.length = 1
#guard countKind isCzGate (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.cnot 0 2, Gate.cz 1 0])) = 0
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cz 0 1, Gate.cnot 0 2, Gate.cz 1 0])).gates.length = 1
#guard countKind isCzGate (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.cnot 2 1, Gate.cz 1 0])) = 2
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cz 0 1, Gate.cnot 2 1, Gate.cz 1 0])).gates.length = 3
#guard countKind isCzGate (runCancel (Circuit.ofGates 4 0
  [Gate.cz 0 1, Gate.ccx 0 1 3, Gate.cz 1 0])) = 0
#guard (runCancel (Circuit.ofGates 4 0
  [Gate.cz 0 1, Gate.ccx 0 1 3, Gate.cz 1 0])).gates.length = 1
#guard countKind isCzGate (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.ccx 1 2 0, Gate.cz 1 0])) = 2
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.ccx 1 2 0, Gate.cz 1 0])).gates.length = 3
#guard countKind isCzGate (runCancel (Circuit.ofGates 2 0
  [Gate.cz 0 1, Gate.x 0, Gate.cz 1 0])) = 2
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cz 0 1, Gate.x 0, Gate.cz 1 0])).gates.length = 3
#guard countKind isCzGate (runCancel (Circuit.ofGates 2 0
  [Gate.cz 0 1, Gate.h 1, Gate.cz 1 0])) = 2
#guard (runCancel (Circuit.ofGates 2 0 [Gate.cz 0 1, Gate.h 1, Gate.cz 1 0])).gates.length = 3
#guard countKind isCzGate (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.x 2, Gate.cz 1 0])) = 0
#guard (runCancel (Circuit.ofGates 3 0 [Gate.cz 0 1, Gate.x 2, Gate.cz 1 0])).gates.length = 1

-- measurement_on_operand_blocks_cz_cancellation
#guard (runCancel (Circuit.ofGates 2 1
  [Gate.cz 0 1, Gate.measure 0 0, Gate.cz 1 0])).gates.length = 3
#guard (runCancel (Circuit.ofGates 2 1
  [Gate.cz 0 1, Gate.measure 0 0, Gate.cz 1 0])).hasMeasurement = true

-- reset_on_operand_blocks_cz_cancellation
#guard (runCancel (Circuit.ofGates 2 0
  [Gate.cz 0 1, Gate.reset 0, Gate.cz 1 0])).gates.length = 3

-- cnot_pair_cancels_through_commuting_cz
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cnot 0 2, Gate.cz 0 1, Gate.cnot 0 2])).gates = [Gate.cz 0 1]

-- cnot_pair_blocked_by_cz_on_target
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cnot 0 2, Gate.cz 2 1, Gate.cnot 0 2])).gates.length = 3

-- alternating_cz_pairs_cancel_to_fixpoint
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.cz 1 2, Gate.t 0, Gate.cz 1 0, Gate.cz 2 1])).gates = [Gate.t 0]

-- measurement_on_other_qubit_does_not_block_cz_cancellation
#guard (runCancel (Circuit.ofGates 3 1
  [Gate.cz 0 1, Gate.measure 2 0, Gate.cz 1 0])).gates = [Gate.measure 2 0]
#guard (runCancel (Circuit.ofGates 3 1
  [Gate.cz 0 1, Gate.measure 2 0, Gate.cz 1 0])).numCbits = 1

-- reset_on_other_qubit_does_not_block_cz_cancellation
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.reset 2, Gate.cz 0 1])).gates = [Gate.reset 2]
#guard (runCancel (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.reset 2, Gate.cz 0 1])).hasMeasurement = true

-- cz_cancellation_through_long_mixed_commuting_run
#guard countKind isCzGate (runCancel (Circuit.ofGates 5 0
  [Gate.cz 0 1, Gate.t 0, Gate.rz (19/100) 1, Gate.cnot 0 2, Gate.ccx 0 1 3, Gate.cz 2 4,
   Gate.x 4, Gate.cz 1 0])) = 1

-- native_cz_cancellation_is_structurally_idempotent
#guard (runCancel (runCancel (Circuit.ofGates 4 0
    [Gate.cz 0 1, Gate.t 0, Gate.cz 2 3, Gate.cz 1 0, Gate.cz 3 2]))).gates
  = (runCancel (Circuit.ofGates 4 0
    [Gate.cz 0 1, Gate.t 0, Gate.cz 2 3, Gate.cz 1 0, Gate.cz 3 2])).gates


/-! ## `src/circuit.rs` -/

-- bell_pair
#guard (Circuit.ofGates 2 0 [Gate.h 0, Gate.cnot 0 1]).gates.length = 2
#guard ((Circuit.ofGates 2 0 [Gate.h 0, Gate.cnot 0 1]).toString.splitOn "h q0").length = 2
#guard ((Circuit.ofGates 2 0 [Gate.h 0, Gate.cnot 0 1]).toString.splitOn "cnot q0, q1").length = 2

-- ghz_state
#guard (Circuit.ofGates 4 0
  [Gate.h 0, Gate.cnot 0 1, Gate.cnot 1 2, Gate.cnot 2 3]).gates.length = 4

-- t_gate_decomposition_of_rz
#guard (Circuit.ofGates 1 0 [Gate.t 0, Gate.s 0, Gate.rz (1/4) 0]).gates.length = 3

-- ccx_gate
#guard (Circuit.ofGates 3 0 [Gate.h 2, Gate.ccx 0 1 2, Gate.h 2]).gates.length = 3
#guard (Circuit.ofGates 3 0 [Gate.h 2, Gate.ccx 0 1 2, Gate.h 2]).hasToffoli = true

-- ccz_gate_display_and_metadata
#guard ((Circuit.ofGates 3 0 [Gate.ccz 2 0 1]).toString.splitOn "ccz q2, q0, q1").length = 2
#guard (Circuit.ofGates 3 0 [Gate.ccz 2 0 1]).hasToffoli = false
#guard (Circuit.ofGates 3 0 [Gate.ccz 2 0 1]).hasCcz = true
#guard (Gate.ccz 2 0 1).qubitsOf = [2, 0, 1]

-- ccz_remap
#guard remapGate (Gate.ccz 8 2 5) [2, 5, 8] = Gate.ccz 2 0 1

-- qft_3qubit
#guard (Circuit.ofGates 3 0
  [Gate.h 0, Gate.rz (1/2) 0, Gate.cnot 1 0, Gate.rz (1/4) 0, Gate.cnot 2 0,
   Gate.h 1, Gate.rz (1/2) 1, Gate.cnot 2 1, Gate.h 2]).gates.length = 9

-- z_gate_display
-- sdg_gate_display
#guard ((Circuit.ofGates 1 0 [Gate.z 0]).toString.splitOn "z q0").length = 2
#guard ((Circuit.ofGates 1 0 [Gate.sdg 0]).toString.splitOn "sdg q0").length = 2

-- cz_gate_display_and_metadata
#guard ((Circuit.ofGates 3 0 [Gate.cz 2 0]).toString.splitOn "cz q2, q0").length = 2
#guard (Circuit.ofGates 3 0 [Gate.cz 2 0]).hasToffoli = false
#guard (Circuit.ofGates 3 0 [Gate.cz 2 0]).hasMeasurement = false

-- cz_qubits_of_preserves_operand_order
#guard (Gate.cz 4 1).qubitsOf = [4, 1]

-- cz_remap
#guard remapGate (Gate.cz 7 3) [3, 7] = Gate.cz 1 0

-- cz_remap_subcircuit
#guard (remapSubcircuit [Gate.t 8, Gate.cz 8 2] [2, 8]).numQubits = 2
#guard (remapSubcircuit [Gate.t 8, Gate.cz 8 2] [2, 8]).gates = [Gate.t 1, Gate.cz 1 0]

-- measure_gate_display
#guard ((Circuit.ofGates 1 1 [Gate.measure 0 0]).toString.splitOn "measure q0 -> c0").length = 2
#guard (Circuit.ofGates 1 1 [Gate.measure 0 0]).hasMeasurement = true

-- reset_gate_display
#guard ((Circuit.ofGates 1 0 [Gate.reset 0]).toString.splitOn "reset q0").length = 2
#guard (Circuit.ofGates 1 0 [Gate.reset 0]).hasMeasurement = true

-- measure_qubits_of
-- reset_qubits_of
#guard (Gate.measure 3 7).qubitsOf = [3]
#guard (Gate.reset 2).qubitsOf = [2]

-- measure_remap
-- reset_remap
#guard remapGate (Gate.measure 5 2) [5] = Gate.measure 0 2
#guard remapGate (Gate.reset 5) [5] = Gate.reset 0

-- with_cbits_default_fields
#guard (Circuit.withCbits 2 3).numQubits = 2
#guard (Circuit.withCbits 2 3).numCbits = 3
#guard (Circuit.withCbits 2 3).hasMeasurement = false
#guard (Circuit.withCbits 2 3).hasToffoli = false
#guard (Circuit.withCbits 2 3).hasCcz = false
#guard (Circuit.withCbits 2 3).gates.length = 0

-- has_measurement_flag_set_by_reset_alone
#guard (Circuit.new 1).hasMeasurement = false
#guard ((Circuit.new 1).apply (Gate.reset 0)).hasMeasurement = true

/-! ## `src/pass.rs` -/

-- count_2q_counts_cnot_and_cz_only
#guard count2q (Circuit.ofGates 3 0 [Gate.h 0, Gate.cnot 0 1, Gate.cz 1 2, Gate.t 0]) = 2

-- count_rz_counts_only_rz_gates
#guard countRz (Circuit.ofGates 1 0 [Gate.rz (1/4) 0, Gate.t 0, Gate.rz (-1/2) 0]) = 2

-- depth_layers_disjoint_gates_and_serializes_shared_qubits
#guard depth (Circuit.ofGates 3 0 [Gate.h 0, Gate.h 1, Gate.cnot 0 2, Gate.t 1]) = 2

-- depth_of_empty_circuit_is_zero
#guard depth (Circuit.new 4) = 0

-- depth_serializes_gates_on_one_qubit
#guard depth (Circuit.ofGates 1 0 [Gate.h 0, Gate.t 0, Gate.z 0]) = 3

-- depth_propagates_through_multi_qubit_dependencies
#guard depth (Circuit.ofGates 3 0 [Gate.cnot 0 1, Gate.x 1, Gate.cz 1 2]) = 3

/-! ## The pass API itself -/

-- Composition of two passes is again a pass, correctness included.
#guard ((Pass.comp CancelGates CancelGates).run
  (Circuit.Checked.of (Circuit.ofGates 1 0 [Gate.h 0, Gate.h 0]) (by decide))).raw.gates.length = 0
#guard (Pass.runAll [CancelGates, CancelGates]
  (Circuit.Checked.of (Circuit.ofGates 2 0
    [Gate.cnot 0 1, Gate.h 0, Gate.h 0, Gate.cnot 0 1]) (by decide))).raw.gates.length = 0


/-! ## `src/cnot_min.rs`

`cnotMin` needs the circuit's qubit count (out-of-range operands are not interpretable), so
the Rust `cnot_min(&c)` becomes `cnotMin c.numQubits c.gates`. As above, the
`circuits_equiv` assertions have no counterpart: `TzapLean.cnotMinGates_correct` proves the
equivalence for every well-formed input. -/

/-- The pass under test, as a plain function. -/
def runCnotMin (c : Circuit) : Circuit := cnotMinCircuit c

-- invert_identity_is_identity
#guard invert ((Array.range 5).map (2 ^ ·)) 5 = some ((Array.range 5).map (2 ^ ·))

-- invert_rejects_singular
#guard invert #[1, 1, 4] 3 = none

-- linear_synth_emits_nothing_for_identity
#guard linearSynth ((Array.range 4).map (2 ^ ·)) ((Array.range 4).map (2 ^ ·)) #[0,1,2,3] = some []

-- gray_synth_on_no_phases_emits_nothing
#guard graySynth 4 [] ((Array.range 4).map (2 ^ ·)) #[0,1,2,3]
  = ([], (Array.range 4).map (2 ^ ·))

-- The native parity path covers both halves of the 128-wire chunk bound.
#guard (FastParity.basis 63).testBit 63
#guard (FastParity.basis 64).testBit 64
#guard !((FastParity.basis 63).xor (FastParity.basis 64)).testBit 62
#guard ((FastParity.basis 63).xor (FastParity.basis 64)).testBit 63
#guard ((FastParity.basis 63).xor (FastParity.basis 64)).testBit 64

-- A spent synthesis budget is a pure early rejection.
#guard (SynthBudget.push { maxCount := 0, maxTwoQ := 0 } (Gate.cnot 0 1)).isNone
#guard (SynthBudget.push { maxCount := 1, maxTwoQ := 1 } (Gate.cnot 0 1)).isSome

-- gray_synth_places_every_rotation_on_its_own_parity (one rotation per singleton parity)
#guard (graySynth 2 [(1, (1:ℚ)/4), (2, (1:ℚ)/4)] #[1, 2] #[0,1]).1.length ≤ 4

-- preserves_a_cnot_ladder
#guard (runCnotMin (Circuit.ofGates 4 0
  [Gate.cnot 0 1, Gate.cnot 1 2, Gate.cnot 2 3])).gates.length = 3

-- cancels_a_repeated_cnot_pair
#guard count2q (runCnotMin (Circuit.ofGates 2 0 [Gate.cnot 0 1, Gate.cnot 0 1])) = 0

-- shrinks_a_redundant_linear_map
#guard count2q (runCnotMin (Circuit.ofGates 3 0
  [Gate.cnot 0 1, Gate.cnot 1 2, Gate.cnot 0 1, Gate.cnot 1 2, Gate.cnot 0 1,
   Gate.cnot 1 2])) < 6

-- merges_rotations_on_one_parity
#guard (runCnotMin (Circuit.ofGates 2 0
  [Gate.cnot 0 1, Gate.t 1, Gate.cnot 0 1, Gate.cnot 0 1, Gate.t 1, Gate.cnot 0 1])).gates.length
  ≤ 6

-- hadamard_splits_blocks_but_preserves_semantics
#guard (runCnotMin (Circuit.ofGates 3 0
  [Gate.cnot 0 1, Gate.t 1, Gate.h 1, Gate.cnot 1 2, Gate.t 2])).gates.contains (Gate.h 1)

-- x_gates_are_preserved_through_resynthesis (the X's affine part survives)
#guard (runCnotMin (Circuit.ofGates 3 0
  [Gate.x 0, Gate.cnot 0 1, Gate.t 1, Gate.x 1, Gate.cnot 1 2, Gate.tdg 2])).gates.length ≤ 6

-- cz_is_absorbed_and_preserved
#guard (runCnotMin (Circuit.ofGates 3 0
  [Gate.cz 0 1, Gate.cnot 1 2, Gate.cz 0 2])).gates.length ≤ 3

-- ccz_and_ccx_act_as_block_boundaries
#guard countKind (fun g => match g with | .ccz .. | .ccx .. => true | _ => false)
  (runCnotMin (Circuit.ofGates 3 0
    [Gate.cnot 0 1, Gate.ccz 0 1 2, Gate.cnot 0 1, Gate.ccx 0 1 2, Gate.t 0])) = 2

-- rz_rotations_survive
#guard (runCnotMin (Circuit.ofGates 2 0
  [Gate.rz (37/100) 0, Gate.cnot 0 1, Gate.rz (-12/10) 1])).gates.length = 3

-- measurement_and_reset_are_preserved_and_bound_blocks
#guard (runCnotMin (Circuit.ofGates 2 1
  [Gate.cnot 0 1, Gate.measure 0 0, Gate.reset 1, Gate.cnot 0 1])).gates[1]?
  = some (Gate.measure 0 0)
#guard (runCnotMin (Circuit.ofGates 2 1
  [Gate.cnot 0 1, Gate.measure 0 0, Gate.reset 1, Gate.cnot 0 1])).gates[2]?
  = some (Gate.reset 1)
#guard (runCnotMin (Circuit.ofGates 2 1
  [Gate.cnot 0 1, Gate.measure 0 0, Gate.reset 1, Gate.cnot 0 1])).hasMeasurement = true
#guard (runCnotMin (Circuit.ofGates 2 1
  [Gate.cnot 0 1, Gate.measure 0 0, Gate.reset 1, Gate.cnot 0 1])).numCbits = 1

-- quarter_turn_angles_emit_clifford_t_gates
#guard (List.range 8).map (fun k => (emitRotation 0 (k/4 : ℚ)).length) = [0, 1, 1, 2, 1, 2, 1, 1]
#guard ((List.range 8).all fun k =>
  (emitRotation 0 (k/4 : ℚ)).all fun g => match g with | .rz .. => false | _ => true)

-- non_quarter_turn_angles_stay_rz
#guard emitRotation 0 (37/100) = [Gate.rz (37/100) 0]

-- accumulated_quarter_turns_stay_clifford_t
#guard ((List.range 7).all fun k =>
  (runCnotMin (Circuit.ofGates 2 0
    ([Gate.cnot 0 1] ++ List.replicate (k+1) (Gate.t 1) ++
      [Gate.cnot 0 1, Gate.cnot 0 1, Gate.cnot 0 1]))).gates.all
    fun g => match g with | .rz .. => false | _ => true)

-- rotations_summing_to_identity_disappear
#guard (runCnotMin (Circuit.ofGates 1 0 (List.replicate 8 (Gate.t 0)))).gates.length = 0

-- empty_and_single_gate_circuits_round_trip
#guard (runCnotMin (Circuit.new 3)).gates.length = 0
#guard (runCnotMin (Circuit.ofGates 2 0 [Gate.h 1])).gates = [Gate.h 1]

-- degenerate_two_qubit_gates_pass_through_untouched
#guard (runCnotMin (Circuit.ofGates 3 0
  [Gate.cnot 0 1, Gate.t 1, Gate.cnot 1 1, Gate.cnot 1 2, Gate.t 2])).gates.count
    (Gate.cnot 1 1) = 1
#guard (runCnotMin (Circuit.ofGates 3 0
  [Gate.cnot 0 1, Gate.t 1, Gate.cz 1 1, Gate.cnot 1 2, Gate.t 2])).gates.count
    (Gate.cz 1 1) = 1

-- a_circuit_of_only_degenerate_gates_is_unchanged
#guard (runCnotMin (Circuit.ofGates 2 0 (List.replicate 3 (Gate.cnot 0 0)))).gates
  = List.replicate 3 (Gate.cnot 0 0)

-- qubits_outside_a_block_are_untouched
#guard (runCnotMin (Circuit.ofGates 40 0
  [Gate.cnot 3 7, Gate.t 7, Gate.cnot 7 11])).gates.all fun g =>
    g.qubitsOf.all fun q => q == 3 || q == 7 || q == 11

-- is_deterministic (the pass is a function, so this is `rfl`)
#guard (runCnotMin (Circuit.ofGates 3 0 [Gate.cnot 0 1, Gate.t 1, Gate.h 0, Gate.cnot 1 2])).gates
  = (runCnotMin (Circuit.ofGates 3 0
      [Gate.cnot 0 1, Gate.t 1, Gate.h 0, Gate.cnot 1 2])).gates

/-! ### A randomized sweep, with the Rust test RNG

`never_increases_gate_or_two_qubit_count`, `never_increases_t_count` and
`repeated_application_reaches_a_fixed_point`, ported with the same xorshift generator the
Rust tests use. -/

/-- The Rust `TestRng`. -/
def rngNext (s : UInt64) (upper : Nat) : Nat × UInt64 :=
  let s := s ^^^ (s <<< 13)
  let s := s ^^^ (s >>> 7)
  let s := s ^^^ (s <<< 17)
  (s.toNat % upper, s)

/-- A random Clifford+T circuit over `n` wires. -/
def randomGates (s : UInt64) (n : Nat) : Nat → List Gate × UInt64
  | 0 => ([], s)
  | k + 1 =>
      let (q, s) := rngNext s n
      let (kind, s) := rngNext s 7
      let (g, s) :=
        if kind == 0 then (Gate.t q, s)
        else if kind == 1 then (Gate.h q, s)
        else if kind == 2 then (Gate.x q, s)
        else if kind == 3 then (Gate.s q, s)
        else if n > 1 then
          let (d, s) := rngNext s (n - 1)
          (Gate.cnot q ((q + 1 + d) % n), s)
        else (Gate.tdg q, s)
      let (rest, s) := randomGates s n k
      (g :: rest, s)

/-- Iterate the pass, as `repeated_application_reaches_a_fixed_point` does: each round must
not grow the gate or two-qubit count, and the output must settle within eight rounds. -/
def settles : Nat → Circuit → Bool
  | 0, _ => false
  | k + 1, c =>
      let next := runCnotMin c
      if next.gates.length ≤ c.gates.length && count2q next ≤ count2q c then
        if next.gates == c.gates then true else settles k next
      else false

/-- One sweep case: the pass never grows the gate count, the two-qubit count or the T-count,
and iterating it reaches a fixed point. -/
def sweepCase (s : UInt64) : Bool × UInt64 :=
  let (n, s) := rngNext s 5
  let n := n + 1
  let (len, s) := rngNext s 24
  let (gs, s) := randomGates s n len
  let c := Circuit.ofGates n 0 gs
  let out := runCnotMin c
  (out.gates.length ≤ c.gates.length && count2q out ≤ count2q c &&
     countT out ≤ countT c && settles 8 c, s)

def sweep : Nat → UInt64 → Bool
  | 0, _ => true
  | k + 1, s => let (ok, s) := sweepCase s; ok && sweep k s

-- never_increases_gate_or_two_qubit_count
-- never_increases_t_count
-- repeated_application_reaches_a_fixed_point
-- random_cnot_dihedral_circuits_are_preserved
-- random_mixed_circuits_are_preserved
-- randomized_equivalence_at_scale
-- cap_split_equivalence_at_scale
#guard sweep 300 0x7777888899993333
#guard sweep 300 0x1234567_9abcdef0


/-! ### `cnot_min` internals

Two Rust tests have no counterpart here and are deliberately not ported:
`budget_is_a_pure_early_exit` and `budget_early_exit_at_scale` pin the Rust `Budget` — a
performance device that abandons a synthesis once it can no longer win — to producing
byte-identical output to the unbounded synthesis. This port *is* the unbounded synthesis, so
there is nothing to compare it against. -/

/-- Apply a `CNOT` list to a linear map, as the Rust `linear_synth` tests do. -/
def applyCnotList (state : List Parity) : List Gate → List Parity
  | [] => state
  | g :: gs =>
      match g with
      | .cnot c t => applyCnotList (state.set t (state[t]! ^^^ state[c]!)) gs
      | _ => applyCnotList state gs

-- invert_composed_with_original_is_identity
#guard (let m : Array Parity := #[3, 2, 4]
        match invert m 3 with
        | some inv => (List.range 3).map (fun i => rowTimesMatrix m[i]! inv)
        | none => []) = (List.range 3).map (2 ^ ·)

-- linear_synth_realizes_random_maps (fixed maps, checked by replaying the CNOTs)
#guard (let n := 4
        let frm : Array Parity := (Array.range n).map (2 ^ ·)
        let tgt : Array Parity := #[1, 3, 5, 8]
        match linearSynth frm tgt (Array.range n) with
        | some gs => applyCnotList frm.toList gs == tgt.toList
        | none => false)

-- linear_synth_realizes_maps_from_a_nonidentity_start
#guard (let frm : Array Parity := #[3, 2, 4]
        let tgt : Array Parity := #[1, 3, 4]
        match linearSynth frm tgt #[0, 1, 2] with
        | some gs => applyCnotList frm.toList gs == tgt.toList
        | none => false)

/-! ### Caps and sweeps

`never_increases_gate_or_two_qubit_count`, `never_increases_t_count`,
`repeated_application_reaches_a_fixed_point`, `is_deterministic`,
`random_cnot_dihedral_circuits_are_preserved`, `random_mixed_circuits_are_preserved`,
`randomized_equivalence_at_scale` and `cap_split_equivalence_at_scale` are all the same
shape: run the pass over randomly generated circuits and check it never grows a count, that
iterating settles, and that the result is equivalent. The equivalence half is
`cnotMinGates_correct`; the counting half is the sweep below, run at the default caps and at
the tightest ones. -/

/-- The pass at explicit caps. -/
def runCnotMinAt (maxQ maxT : Nat) (c : Circuit) : Circuit :=
  c.withGates (cnotMinGates c.numQubits maxQ maxT c.gates)

def settlesAt (maxQ maxT : Nat) : Nat → Circuit → Bool
  | 0, _ => false
  | k + 1, c =>
      let next := runCnotMinAt maxQ maxT c
      if next.gates.length ≤ c.gates.length && count2q next ≤ count2q c then
        if next.gates == c.gates then true else settlesAt maxQ maxT k next
      else false

def sweepCaseAt (maxQ maxT : Nat) (nQubits : Nat) (s : UInt64) : Bool × UInt64 :=
  let (n, s) := rngNext s nQubits
  let n := n + 1
  let (len, s) := rngNext s 24
  let (gs, s) := randomGates s n len
  let c := Circuit.ofGates n 0 gs
  let out := runCnotMinAt maxQ maxT c
  (out.gates.length ≤ c.gates.length && count2q out ≤ count2q c &&
     countT out ≤ countT c && settlesAt maxQ maxT 8 c, s)

def sweepAt (maxQ maxT nQubits : Nat) : Nat → UInt64 → Bool
  | 0, _ => true
  | k + 1, s => let (ok, s) := sweepCaseAt maxQ maxT nQubits s; ok && sweepAt maxQ maxT nQubits k s

-- tight_qubit_cap_preserves_semantics (the counting half; equivalence is proved)
#guard sweepAt 2 512 5 120 0x51de51de51de51de

-- tight_term_cap_preserves_semantics
#guard sweepAt 128 1 5 120 0x7e2c7e2c7e2c7e2c

-- term_cap_sweep
#guard ([1, 2, 4, 8].all fun t => sweepAt 128 t 5 40 0x1111222233334444)

-- circuits_wider_than_the_qubit_cap_are_handled
#guard sweepAt 4 512 12 60 0x0fedcba987654321

-- a_wide_shallow_block_is_abandoned_not_expanded: a wide block never grows the circuit
#guard (let n := 12
        let gs := (List.range (n-1)).map (fun i => Gate.cnot i (i+1)) ++
          (List.range n).map (fun i => Gate.t i)
        let c := Circuit.ofGates n 0 gs
        (runCnotMinAt 4 512 c).gates.length ≤ c.gates.length)

end TzapLean
