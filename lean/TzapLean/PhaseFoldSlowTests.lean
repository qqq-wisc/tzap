import TzapLean.PhaseFoldTests

/-!
# Slow PhaseFoldRand parity cases

Run explicitly with `lake build TzapLean.PhaseFoldSlowTests`. These replay the three Rust
`#[ignore]` tests without placing their larger circuits in the normal build. These replay
the inputs and T-count assertions. The four `mod5` numerical equivalence assertions are
replaced by exact Clifford+T matrix checks, evaluated with `native_decide`. The much longer
five-qubit pipeline remains covered by the general good-sample theorem, not an exact matrix
evaluation here.
-/

namespace TzapLean

/-- The exact 330-gate input from Rust's `large_circuit_toffoli_decompose_and_phase_fold`. -/
def largeRustPipeline : List Gate :=
  [ .cnot 3 2
  , .cnot 8 7
  , .cnot 14 13
  , .cnot 21 20
  , .cnot 3 4
  , .cnot 8 9
  , .cnot 14 15
  , .cnot 21 22
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .h 3
  , .ccx 0 1 3
  , .h 3
  , .h 8
  , .ccx 5 6 8
  , .h 8
  , .h 14
  , .ccx 11 12 14
  , .h 14
  , .h 21
  , .ccx 18 19 21
  , .h 21
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .h 4
  , .h 9
  , .h 15
  , .h 22
  , .h 10
  , .h 16
  , .h 23
  , .h 4
  , .ccx 2 3 4
  , .h 4
  , .h 9
  , .ccx 7 8 9
  , .h 9
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 15
  , .ccx 13 14 15
  , .h 15
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .h 22
  , .ccx 20 21 22
  , .h 22
  , .h 23
  , .ccx 20 21 23
  , .h 23
  , .cnot 6 5
  , .cnot 12 11
  , .cnot 19 18
  , .cnot 5 8
  , .cnot 11 14
  , .cnot 18 21
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .h 23
  , .ccx 20 21 23
  , .h 23
  , .h 4
  , .h 10
  , .h 15
  , .h 16
  , .h 23
  , .h 17
  , .h 17
  , .ccx 16 23 17
  , .h 17
  , .h 22
  , .ccx 15 23 22
  , .h 22
  , .h 9
  , .ccx 4 10 9
  , .h 9
  , .h 17
  , .h 9
  , .h 15
  , .h 22
  , .ccx 9 17 22
  , .h 22
  , .h 15
  , .ccx 9 16 15
  , .h 15
  , .h 15
  , .h 22
  , .h 17
  , .h 17
  , .ccx 16 23 17
  , .h 17
  , .h 17
  , .h 10
  , .h 16
  , .h 23
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .h 23
  , .ccx 20 21 23
  , .h 23
  , .cnot 5 8
  , .cnot 11 14
  , .cnot 18 21
  , .cnot 6 5
  , .cnot 12 11
  , .cnot 19 18
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .h 23
  , .ccx 20 21 23
  , .h 23
  , .h 23
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .h 3
  , .ccx 0 1 3
  , .h 3
  , .h 8
  , .ccx 5 6 8
  , .h 8
  , .h 14
  , .ccx 11 12 14
  , .h 14
  , .h 21
  , .ccx 18 19 21
  , .h 21
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .cnot 3 2
  , .cnot 8 7
  , .cnot 14 13
  , .cnot 21 20
  , .cnot 6 5
  , .cnot 12 11
  , .cnot 19 18
  , .cnot 6 8
  , .cnot 12 14
  , .cnot 19 21
  , .cnot 4 6
  , .cnot 9 12
  , .cnot 15 19
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .h 3
  , .ccx 0 1 3
  , .h 3
  , .h 8
  , .ccx 5 6 8
  , .h 8
  , .h 14
  , .ccx 11 12 14
  , .h 14
  , .h 21
  , .ccx 18 19 21
  , .h 21
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .cnot 3 2
  , .cnot 8 7
  , .cnot 14 13
  , .cnot 21 20
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .h 3
  , .ccx 0 1 3
  , .h 3
  , .h 8
  , .ccx 5 6 8
  , .h 8
  , .h 14
  , .ccx 11 12 14
  , .h 14
  , .h 21
  , .ccx 18 19 21
  , .h 21
  , .h 3
  , .h 8
  , .h 14
  , .h 21
  , .cnot 6 5
  , .cnot 12 11
  , .cnot 19 18
  , .cnot 4 6
  , .cnot 9 12
  , .cnot 15 19
  , .cnot 6 8
  , .cnot 12 14
  , .cnot 19 21
  , .cnot 1 0
  , .cnot 6 5
  , .cnot 12 11
  , .cnot 19 18
  , .x 0
  , .x 2
  , .x 5
  , .x 7
  , .x 11
  , .x 13
  , .cnot 3 2
  , .cnot 8 7
  , .cnot 14 13
  , .h 3
  , .h 8
  , .h 14
  , .h 3
  , .ccx 0 1 3
  , .h 3
  , .h 8
  , .ccx 5 6 8
  , .h 8
  , .h 14
  , .ccx 11 12 14
  , .h 14
  , .h 3
  , .h 8
  , .h 14
  , .cnot 6 5
  , .cnot 12 11
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .cnot 5 8
  , .cnot 11 14
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .h 10
  , .h 16
  , .h 15
  , .h 15
  , .ccx 9 16 15
  , .h 15
  , .h 9
  , .h 9
  , .ccx 4 10 9
  , .h 9
  , .h 4
  , .h 10
  , .h 16
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .cnot 5 8
  , .cnot 11 14
  , .h 10
  , .ccx 7 8 10
  , .h 10
  , .h 16
  , .ccx 13 14 16
  , .h 16
  , .cnot 6 5
  , .cnot 12 11
  , .h 9
  , .ccx 7 8 9
  , .h 9
  , .h 15
  , .ccx 13 14 15
  , .h 15
  , .h 4
  , .ccx 2 3 4
  , .h 4
  , .h 4
  , .h 9
  , .h 10
  , .h 15
  , .h 16
  , .h 3
  , .h 8
  , .h 14
  , .h 3
  , .ccx 0 1 3
  , .h 3
  , .h 8
  , .ccx 5 6 8
  , .h 8
  , .h 14
  , .ccx 11 12 14
  , .h 14
  , .h 3
  , .h 8
  , .h 14
  , .cnot 3 4
  , .cnot 8 9
  , .cnot 14 15
  , .cnot 3 2
  , .cnot 8 7
  , .cnot 14 13
  , .x 0
  , .x 2
  , .x 5
  , .x 7
  , .x 11
  , .x 13 ]

#guard
  let decomposed := decomposeToffoliGates largeRustPipeline
  let once := pf 24 decomposed
  let twice := pf 24 once
  largeRustPipeline.length == 330 &&
    countTGates once ≤ countTGates decomposed &&
    countTGates twice ≤ countTGates once

/-- Rust's four-CCX `mod5_4_remove_ccx_combos` input. -/
def mod5Base : List Gate :=
  [.x 4, .ccx 0 3 4, .ccx 2 3 4, .cnot 3 4, .ccx 1 2 4,
    .cnot 2 4, .ccx 0 1 4, .cnot 1 4, .cnot 0 4]

def mod5Case (keep : Nat) : List Gate :=
  (mod5Base.zipIdx).filterMap fun (gate, i) =>
    if ([1, 2, 4, 6] : List Nat).contains i && i != ([1, 2, 4, 6] : List Nat).getD keep 1
    then none else some gate

#guard (List.range 4).all fun keep =>
  let input := mod5Case keep
  let decomposed := decomposeToffoliGates input
  let folded := pf 5 decomposed
  input.length == 6 && countTGates folded ≤ countTGates decomposed

private theorem mod5Cases_exact :
    (List.range 4).all (fun keep =>
      let decomposed := decomposeToffoliGates (mod5Case keep)
      exactMatrixEquivalent 5 decomposed (pf 5 decomposed)) = true := by
  native_decide

/- Rust's ignored 5-qubit multi-pass case. -/
#guard
  let decomposed := decomposeToffoliGates smallToffoliPipeline
  let once := pf 5 decomposed
  let twice := pf 5 once
  countTGates once ≤ countTGates decomposed &&
    countTGates twice ≤ countTGates once

end TzapLean
