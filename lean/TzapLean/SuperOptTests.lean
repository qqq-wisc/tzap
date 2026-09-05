import TzapLean.SuperOptProof

/-!
# `SuperOpt`: tests

The behavioural half of `src/super_opt/tests.rs`. Rust's remaining tests cover machinery this
port does not have — the on-disk cache, the matrix-store arena, incremental mode, the
subcircuit/rewrite reports — or assert `circuits_equiv`, which here is a theorem for every
input (`superOptGates_correct`).

The exact-matrix checks come first: they are what a floating-point implementation would get
subtly wrong, and what the cyclotomic representation exists to get right.

The pass tests share **one** table. Each `#guard` is evaluated independently, so a table
built inside each would be rebuilt each time; `passFailures` therefore builds it once and
runs every case against it. When a case fails, `#eval passFailures` names it.
-/

namespace TzapLean

open Gate

/-! ## Exact matrices -/

/-- Two gate lists with the same exact matrix, up to one of the eight Clifford+T phases. -/
def eqUpToPhase (k : Nat) (gs hs : List Gate) : Bool :=
  match ExactMat.matrixOf k gs, ExactMat.matrixOf k hs with
  | some A, some B => (ExactMat.phaseMatch A.normalize B.normalize).isSome
  | _, _ => false

/-! `x` and `−x` are the same operator; `s` and `sdg` are not. A representation that rounded,
or that canonicalized phase by dividing, could get either of these wrong. -/

-- `phase_equivalence_ignores_global_phase`: `z·x·z = -x`.
#guard eqUpToPhase 1 [x 0] [z 0, x 0, z 0]
#guard !eqUpToPhase 1 [s 0] [sdg 0]
-- `identity_matches_omega_identity`, `hh_is_identity`.
#guard eqUpToPhase 1 [] [h 0, h 0]
#guard eqUpToPhase 1 [] [x 0, x 0]
#guard eqUpToPhase 1 [] [s 0, s 0, s 0, s 0]
-- `ccx_and_ccz_have_different_keys`.
#guard !eqUpToPhase 3 [ccx 0 1 2] [ccz 0 1 2]
-- Control and target are not interchangeable.
#guard !eqUpToPhase 2 [cnot 0 1] [cnot 1 0]
-- Textbook identities, as matrix facts.
#guard eqUpToPhase 2 [h 1, cnot 0 1, h 1] [cz 0 1]
#guard eqUpToPhase 1 [t 0, t 0, t 0, t 0] [z 0]
#guard eqUpToPhase 1 [h 0, z 0, h 0] [x 0]
#guard eqUpToPhase 2 [x 0, cnot 0 1, x 0] [x 1, cnot 0 1]
#guard eqUpToPhase 2 [t 0, cnot 0 1, tdg 0] [cnot 0 1]
#guard eqUpToPhase 2 [cnot 0 1, cnot 1 0, cnot 0 1] [cnot 1 0, cnot 0 1, cnot 1 0]

/-! ## The builder's fast paths agree with the definitions they replace -/

/-- Sample coefficient tuples, including negatives and zeros. -/
def sampleCycs : List Cyc :=
  [⟨0,0,0,0⟩, ⟨1,0,0,0⟩, ⟨0,1,0,0⟩, ⟨0,0,1,0⟩, ⟨0,0,0,1⟩,
   ⟨1,2,3,4⟩, ⟨-1,2,-3,4⟩, ⟨127,-127,3,-9⟩, ⟨-5,-6,-7,-8⟩]

-- The one-step `ω^p` the table builder uses is the recursive one, for every `p` mod 8.
#guard sampleCycs.all fun x => (List.range 24).all fun p => x.timesOmegaFast p == x.timesOmega p

/-! ## Canonical window-shape cache -/

-- Physical wire names disappear from the key, but operand order does not.
#guard compactShapeKey [3, 7] #[h 3, cnot 3 7, t 7] [0, 1, 2] ==
  compactShapeKey [11, 2] #[h 11, cnot 11 2, t 2] [0, 1, 2]
#guard compactShapeKey [3, 7] #[cnot 3 7] [0] !=
  compactShapeKey [3, 7] #[cnot 7 3] [0]

/-! ## Batched rewrite certificates -/

def sampleTagged : List Tagged :=
  [(h 0, some 0), (x 2, none), (t 0, some 0),
   (h 1, some 1), (z 2, none), (h 1, some 1)]

def sampleRepl : Nat → List Gate
  | 0 => [s 0]
  | 1 => []
  | _ => []

#guard groupedClaim sampleTagged 0 = claimedBy 0 sampleTagged
#guard groupedClaim sampleTagged 1 = claimedBy 1 sampleTagged
#guard applyAllLinear sampleRepl [] sampleTagged = applyAll sampleRepl sampleTagged
#guard sepAllB (fun w => if w = 0 then [0] else [1]) sampleTagged =
  sepB (fun w => if w = 0 then [0] else [1]) [] sampleTagged

def blockedTagged : List Tagged := [(h 0, some 0), (x 0, none), (h 0, some 0)]
#guard !sepAllB (fun _ => [0]) blockedTagged
#guard sepAllB (fun _ => [0]) blockedTagged = sepB (fun _ => [0]) [] blockedTagged

/-! ## The library gate set

`library_gates(k)` has `7k` one-wire gates and `k(k−1)` `CNOT`s — Rust asserts exactly these
counts. `ccx` and `cz` are absent by design, so a rewrite never introduces them. -/

#guard (libGates 1).length == 7
#guard (libGates 2).length == 16
#guard (libGates 3).length == 27
#guard (libGates 4).length == 40

/-! ## The table

Rust asserts a depth-1 width-1 table holds 8 unitaries and reports depth 1. -/

/-- A depth-1 width-1 table, for the counts Rust checks. -/
def tinyTable : WidthTable := buildWidth 1 { maxQubits := 1, maxGates := 1 }

#guard tinyTable.size == 8
#guard tinyTable.depth == 1
#guard !tinyTable.saturated

-- The entry cap stops the build and says so.
#guard (buildWidth 2 { maxQubits := 2, maxGates := 4, maxEntriesPerQubit := 20 }).saturated

/-! ## The pass -/

/-- Two wires, replacements of up to three gates. -/
def tcfg : SuperOptTableConfig := { maxQubits := 2, maxGates := 3 }

/-- Windows of up to six gates on up to two wires. -/
def cfg : SuperOptConfig := { maxQubits := 2, maxWindow := 6 }

/-- Each case is a name, a wire count, an input, and the expected output. -/
def passCases : List (String × Nat × List Gate × List Gate) :=
  [ -- `empty_circuit_is_unchanged`, `single_gate_is_unchanged`
    ("empty", 1, [], []),
    ("single gate", 1, [h 0], [h 0]),
    -- self-inverse pairs, found by lookup rather than by rule
    ("hh", 1, [h 0, h 0], []),
    ("xx", 1, [x 0, x 0], []),
    ("t tdg", 1, [t 0, tdg 0], []),
    ("ssss", 2, [s 0, s 0, s 0, s 0], []),
    ("cx cx", 2, [cnot 0 1, cnot 0 1], []),
    -- rotations fold
    ("tt", 1, [t 0, t 0], [s 0]),
    ("ss", 1, [s 0, s 0], [z 0]),
    ("tttt", 1, [t 0, t 0, t 0, t 0], [z 0]),
    ("s tttt", 1, [s 0, t 0, t 0, t 0, t 0], [sdg 0]),
    -- conjugations the pass discovers
    ("h z h", 1, [h 0, z 0, h 0], [x 0]),
    ("h x h", 1, [h 0, x 0, h 0], [z 0]),
    ("(hs)^3", 1, [h 0, s 0, h 0, s 0, h 0, s 0], []),
    ("x on control", 2, [x 0, cnot 0 1, x 0], [x 1, cnot 0 1]),
    ("t through control", 2, [t 0, cnot 0 1, tdg 0], [cnot 0 1]),
    -- physical support order must not reverse a synthesized CNOT
    ("sparse cnot forward", 5, [h 3, cz 1 3, h 3], [cnot 1 3]),
    ("sparse cnot reverse", 5, [h 1, cz 3 1, h 1], [cnot 3 1]),
    -- CZ is not in the library, so `h·cx·h` has no replacement the pass may emit
    ("h cx h stays", 2, [h 1, cnot 0 1, h 1], [h 1, cnot 0 1, h 1]),
    -- a SWAP is already optimal
    ("swap stays", 2, [cnot 0 1, cnot 1 0, cnot 0 1], [cnot 0 1, cnot 1 0, cnot 0 1]),
    -- windows are subsequences: a gate on other wires is skipped and re-emitted
    ("skip disjoint", 2, [t 0, h 1, t 0], [s 0, h 1]),
    ("skip reset", 2, [h 0, reset 1, h 0], [reset 1]),
    ("skip measure", 2, [h 0, measure 1 0, h 0], [measure 1 0]),
    -- …but a barrier on the window's own wire stops it, and both sides still optimize
    ("measure splits", 1, [h 0, h 0, measure 0 0, h 0, h 0], [measure 0 0]),
    ("rz blocks", 1, [h 0, rz (1/3) 0, h 0], [h 0, rz (1/3) 0, h 0]),
    ("rz then pair", 1, [rz (1/3) 0, h 0, h 0], [rz (1/3) 0]),
    -- a three-wire gate is out of a two-wire window's reach
    ("ccx out of reach", 3, [ccx 0 1 2, ccx 0 1 2], [ccx 0 1 2, ccx 0 1 2]),
    -- the Toffoli decomposition has no two-wire win: its `T`s sit on three-wire parities
    ("toffoli decomp", 3,
      [h 2, cnot 1 2, tdg 2, cnot 0 2, t 2, cnot 1 2, tdg 2, cnot 0 2, t 1, t 2, h 2],
      [h 2, cnot 1 2, tdg 2, cnot 0 2, t 2, cnot 1 2, tdg 2, cnot 0 2, t 1, t 2, h 2]) ]

/-! ## Retroactive closure

A gate that brings a new wire into the window pulls in the earlier gates on that wire too.
Before the window was closed over its wires, a skipped gate touching the newly-bridged wire
killed the window outright, and every rewrite that has to be discovered that way was
systematically missed. Each of these needs exactly that: the `h q1` is skipped while the
window sits on wire 0, and only becomes reachable when the `cnot` bridges the two. -/
def bridgeCases : List (String × Nat × List Gate × List Gate) :=
  [ ("bridge x", 2, [x 0, h 1, cnot 0 1, x 0], [h 1, x 1, cnot 0 1]),
    ("bridge t", 2, [t 0, h 1, cnot 0 1, tdg 0], [h 1, cnot 0 1]) ]

def bridgeFailures : List String :=
  let tbl := buildTable tcfg
  bridgeCases.filterMap fun (name, n, inp, want) =>
    if superOptGates cfg tbl n inp == want then none else some name

#guard bridgeFailures.isEmpty

/-! ## Greedy completion order

The `ss` window on wire 1 completes before the older `x` window on wire 0 reaches it through
the final bridge. Rust therefore claims the two `s` gates first. This distinguishes the
simultaneous-window scan from the old earliest-anchor scan, which looked through the bridge
before giving the `s` anchor a turn. -/
def completionOrderOk : Bool :=
  let inp := [x 0, s 1, s 1, cnot 0 1]
  let tbl := buildTable tcfg
  let st := proposeRewrites cfg tbl 2 inp.toArray
  st.tags.toList == [none, some 0, some 0, none] &&
    superOptGates cfg tbl 2 inp == [x 0, z 1, cnot 0 1]

#guard completionOrderOk

/-- The second pair has the same support-local shape as the first pair on another physical
wire, so its synthesis answer comes from the scan-local cache. -/
def shapeCacheReusesAnswer : Bool :=
  let tbl := buildTable tcfg
  let inp := #[h 0, h 0, h 1, h 1]
  let st := proposeRewrites cfg tbl 2 inp
  st.shapeHits > 0 && st.shapeMisses > 0 &&
    superOptGates cfg tbl 2 inp.toList == []

#guard shapeCacheReusesAnswer

/-- Iterate the scan to a fixpoint, as the driver does.

`superOptGates` is *one forward scan*, because `SuperOpt::run` is one forward scan: a case
needing two rewrites on the same wires needs two of them, and it is the pipeline that repeats
the pass, not the pass that repeats itself. -/
def superOptFix (tbl : SynthTable) (n : Nat) : Nat → List Gate → List Gate
  | 0, gs => gs
  | fuel + 1, gs =>
      let gs' := superOptGates cfg tbl n gs
      if gs'.length < gs.length then superOptFix tbl n fuel gs' else gs'

/-- Cases whose fixpoint differs from what is expected — empty when all pass. One table serves
every case. -/
def passFailures : List String :=
  let tbl := buildTable tcfg
  passCases.filterMap fun (name, n, inp, want) =>
    if superOptFix tbl n inp.length inp == want then none else some name

#guard passFailures.isEmpty

/-- The cases above where a *single* scan does not already finish, and what it leaves for the
driver's next round. Recorded rather than hidden: this is the whole difference between one
invocation of the pass and the pipeline that repeats it. -/
def oneScanCases : List (String × Nat × List Gate × List Gate) :=
  [ ("ssss", 2, [s 0, s 0, s 0, s 0], [z 0, z 0]),
    ("tttt", 1, [t 0, t 0, t 0, t 0], [s 0, s 0]),
    ("s tttt", 1, [s 0, t 0, t 0, t 0, t 0], [z 0, s 0]),
    ("(hs)^3", 1, [h 0, s 0, h 0, s 0, h 0, s 0], [sdg 0, h 0, h 0, s 0]) ]

def oneScanFailures : List String :=
  let tbl := buildTable tcfg
  oneScanCases.filterMap fun (name, n, inp, want) =>
    if superOptGates cfg tbl n inp == want then none else some name

#guard oneScanFailures.isEmpty

/-- Every other case *is* finished by a single scan. -/
def singleScanSuffices : List String :=
  let tbl := buildTable tcfg
  passCases.filterMap fun (name, n, inp, want) =>
    if oneScanCases.any (·.1 == name) then none
    else if superOptGates cfg tbl n inp == want then none else some name

#guard singleScanSuffices.isEmpty

/-- Iterating reaches a fixpoint: one more scan over the settled circuit changes nothing. -/
def fixpointFailures : List String :=
  let tbl := buildTable tcfg
  passCases.filterMap fun (name, n, inp, _) =>
    let done := superOptFix tbl n inp.length inp
    if superOptGates cfg tbl n done == done then none else some name

#guard fixpointFailures.isEmpty

/-! ## Three-wire windows

At width 3 the table reaches `ccx` and `ccz` windows, which two-wire windows cannot. -/

/-- Cases needing a three-wire table. -/
def wide3Failures : List String :=
  let tbl := buildTable { maxQubits := 3, maxGates := 1 }
  let cfg3 : SuperOptConfig := { maxQubits := 3, maxWindow := 6 }
  [("ccx pair", [ccx 0 1 2, ccx 0 1 2]),
   ("ccz pair", [ccz 0 1 2, ccz 0 1 2]),
   ("cz pair", [Gate.cz 0 1, Gate.cz 0 1])].filterMap fun (name, inp) =>
    if superOptGates cfg3 tbl 3 inp == ([] : List Gate) then none else some name

#guard wide3Failures.isEmpty

end TzapLean
