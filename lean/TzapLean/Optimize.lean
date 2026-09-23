import TzapLean.SuperOptProof
import TzapLean.PhaseFoldRand
import TzapLean.GF128Bridge
import TzapLean.Qasm
import TzapLean.MurmCache
import TzapLean.Pipeline

/-!
# The Optimizer Driver

A port of the parts of `src/optimize.rs` this development can support: the optimization
levels, the pass names selectable by `--passes`, the fixpoint loop, and the metrics the CLI
reports.

The deterministic passes are `Pass`es and carry unconditional proofs. Phase folding is
executed with fresh 128-bit tags, like Rust, and is modelled by `RandPass`: its probability of
failure is machine-checked under ideal independent uniform sampling. `ExecutableRandPass`
keeps the OS entropy boundary explicit; no theorem here claims that the platform RNG itself is
uniform.

Rz decomposition, Clifford re-synthesis, and parallel chunking are not ported. `Level.O3` is
therefore `O2` run to a true fixpoint, without Rust's one-shot Clifford re-synthesis at the end.
-/

namespace TzapLean

/-! ## Metrics -/

/-- The counters the CLI reports, from one walk of the gate list. -/
structure Metrics where
  /-- Total gates. -/
  gates : Nat := 0
  /-- `cnot` and `cz` gates. -/
  twoQubit : Nat := 0
  /-- Circuit depth. -/
  depth : Nat := 0
  /-- `t` and `tdg` gates. -/
  t : Nat := 0
  /-- `rz` gates. -/
  rz : Nat := 0
deriving Repr, Inhabited, DecidableEq

/-- Measure a structurally valid circuit in one gate-list walk.

`TzapLean.depth` is the readable specification: it builds a chain of functions and then
queries that chain once per qubit, which costs `O(gates × qubits)`. The driver has already
established that every operand is below `numQubits`, so its reporting path can keep the same
next-free-layer state in an array and compute all five counters alongside it in `O(gates)`.
Metrics are presentation only; no optimizer decision or correctness proof depends on them. -/
def Metrics.of (c : RawCircuit) : Metrics := Id.run do
  let mut next : Array Nat := Array.replicate c.numQubits 0
  let mut gates := 0
  let mut twoQubit := 0
  let mut depth := 0
  let mut t := 0
  let mut rz := 0
  for g in c.gates do
    gates := gates + 1
    match g with
    | .cnot .. | .cz .. => twoQubit := twoQubit + 1
    | .t _ | .tdg _ => t := t + 1
    | .rz _ _ => rz := rz + 1
    | _ => pure ()
    let mut layer := 1
    for q in g.qubitsOf do
      if q < next.size then layer := max layer (next[q]! + 1)
    for q in g.qubitsOf do
      if q < next.size then next := next.set! q layer
    depth := max depth layer
  return ⟨gates, twoQubit, depth, t, rz⟩

/-! ## Levels and passes -/

/-- Which default pipeline to run. -/
inductive Level where
  /-- Randomized phase folding + gate cancellation. Fastest. -/
  | O1
  /-- Adds `SuperOpt`, capped at two rounds. -/
  | O2
  /-- Like `O2`, run to a true fixpoint. The default. -/
  | O3
deriving Repr, DecidableEq, Inhabited

/-- A pass selectable by name in `--passes`. -/
inductive PassName where
  /-- Lower CCX and CCZ exactly to Clifford+T. -/
  | DecomposeToffoli
  /-- Lower CZ exactly to H/CX/H. -/
  | DecomposeCz
  /-- Cancel adjacent self-inverse pairs and reduce Hadamards. -/
  | CancelGates
  /-- Re-synthesize CNOT-dihedral blocks. -/
  | CnotMin
  /-- Peephole superoptimization against a MURM. -/
  | SuperOpt
  /-- Randomized nonlinear phase folding with fresh 128-bit field tags. -/
  | PhaseFoldRand
deriving Repr, DecidableEq, Inhabited

namespace PassName

/-- All passes — name, variant, and what it carries — in a stable order for listing. -/
def all : List (String × PassName × String) :=
  [ ("DecomposeToffoli", .DecomposeToffoli,
     "Decompose CCX and CCZ exactly into Clifford+T"),
    ("DecomposeCz", .DecomposeCz,
     "Decompose CZ exactly into H, CX, H"),
    ("CancelGates", .CancelGates,
     "Cancel adjacent self-inverse gate pairs and reduce Hadamards"),
    ("CnotMin", .CnotMin,
     "Re-synthesize CNOT-dihedral blocks to cut CNOT count"),
    ("SuperOpt", .SuperOpt,
     "Peephole superoptimization against an exact MURM"),
    ("PhaseFoldRand", .PhaseFoldRand,
     "Merge rotations on the same phase polynomial (randomized 128-bit tags)") ]

/-- Parse a pass name. -/
def parse (s : String) : Option PassName :=
  (all.find? (·.1 == s)).map (·.2.1)

/-- Every pass name, comma-separated — for error messages. -/
def allNames : String := String.intercalate ", " (all.map (·.1))

/-- Whether a pass has an unconditional proof for every run. Phase folding instead carries
a proved failure-probability bound. -/
def verified : PassName → Bool
  | .PhaseFoldRand => false
  | _ => true

end PassName

/-! ## Options -/

/-- Which gate families SuperOpt may emit. -/
inductive SuperOptGates where
  /-- Base Clifford+T plus CZ/CCX/CCZ present in the current stage. -/
  | auto
  /-- The fixed Clifford+T base basis. -/
  | base
  /-- Exactly this nonempty subset of the supported synthesis gates. -/
  | explicit (gates : GateSet)
deriving Repr, Inhabited, DecidableEq

namespace SuperOptGates

/-- Resolve a basis from the circuit entering this particular SuperOpt occurrence. -/
def effective (mode : SuperOptGates) (stage : GateSet) : GateSet :=
  match mode with
  | .auto => baseGateSet.union (stage.inter optionalGateSet)
  | .base => baseGateSet
  | .explicit gates => gates

/-- Parse `auto`, `base`, or an exact comma-separated gate list. -/
def parse (raw : String) : Except String SuperOptGates :=
  if raw == "auto" then .ok .auto
  else if raw == "base" then .ok .base
  else
    let names := (raw.splitOn ",").map (·.trimAscii.toString) |>.filter (· ≠ "")
    if names.isEmpty then .error "--superopt-gates requires auto, base, or a nonempty gate list"
    else do
      let mut gates := GateSet.empty
      for name in names do
        let some kind := GateKind.parse name |
          throw s!"unsupported SuperOpt gate '{name}'"
        if !supportedGateSet.contains kind then
          throw s!"gate '{name}' cannot be used for SuperOpt synthesis"
        gates := gates.insert kind
      if gates.isEmpty then .error "--superopt-gates requires a nonempty gate list"
      else .ok (.explicit gates)

end SuperOptGates

/-- `SuperOpt` window and MURM bounds, `none` meaning "whatever the level implies". -/
structure SuperOptBounds where
  /-- Widest window and MURM, in wires. -/
  qubits : Option Nat := none
  /-- Longest window, in gates. -/
  windowGates : Option Nat := none
  /-- Cap on stored unitaries per MURM width. -/
  murmEntries : Option Nat := none
deriving Repr, Inhabited

/-- Everything the driver needs. -/
structure Options where
  /-- Which default pipeline, when `passes` is absent. -/
  level : Level := .O3
  /-- An explicit pipeline, overriding `level`. -/
  passes : Option (List PassName) := none
  /-- Repeat the pipeline until the gate count stops decreasing. -/
  fixpoint : Bool := false
  /-- `SuperOpt` bounds overrides. -/
  superopt : SuperOptBounds := {}
  /-- How the synthesis basis is selected at each SuperOpt stage. -/
  superoptGates : SuperOptGates := .auto
  /-- Lower CCX and CCZ between the two optimization stages. -/
  decomposeCcx : Bool := false
  /-- Lower CZ between the two optimization stages. -/
  decomposeCz : Bool := false
  /-- Compatibility value accepted from `--seed`; OS-backed phase folding ignores it. -/
  seed : Option Nat := none
  /-- Print detailed input and MURM-loading information. -/
  verbose : Bool := false
deriving Repr, Inhabited

/-- Optional native families that a requested decomposition stage must not allow SuperOpt to
reintroduce afterward. -/
def Options.decomposedGateSet (o : Options) : GateSet :=
  let set := if o.decomposeCcx then GateSet.ofKinds [.ccx, .ccz] else GateSet.empty
  if o.decomposeCz then set.insert .cz else set

/-- The window/MURM bounds a level implies: Rust's own — 3 wires, 25-gate windows, a
200,000-entry MURM, with the same `table_gates = window_gates - 1` mapping.

That MURM takes about 76 seconds to build here against Rust's parallel builder, which is
affordable only because it is built once and cached (`MurmCache`): a warm run loads its
549,456 unitaries in 0.07 s. `--superopt-qubits`, `--superopt-window-gates` and
`--superopt-table-entries` override any of the three. -/
def Level.bounds : Level → Nat × Nat × Nat
  | _ => (3, 25, 200000)

/-- Resolve the bounds for a run: level preset, then any explicit override. -/
def resolveBounds (o : Options) (stage : GateSet) (forbidden : GateSet := GateSet.empty) :
    SuperOptConfig × MurmConfig :=
  let (q, w, e) := o.level.bounds
  let q := o.superopt.qubits.getD q
  let w := o.superopt.windowGates.getD w
  let e := o.superopt.murmEntries.getD e
  let basis := (o.superoptGates.effective stage).diff forbidden
  -- A MURM entry only ever replaces a strictly larger window, so the MURM never needs to
  -- be deeper than `windowGates - 1`.
  ({ maxQubits := q, maxWindow := w }, { maxQubits := q, maxGates := w - 1,
                                         maxEntriesPerQubit := e, basis })

/-- Whether a level's pipeline includes `SuperOpt`, and so pays for a MURM. -/
def Level.usesSuperOpt : Level → Bool
  | .O1 => false
  | _ => true

/-! ## Idealized randomized model

This section is the ideal probabilistic model of the executable pipeline. The runtime obtains
its samples from the operating system rather than evaluating the `PMF`; the bound therefore
applies under the explicit assumption that those samples realize the model's independent
uniform draws. -/

/-- Tag width for phase folding, matching Rust's `u128` fingerprints. -/
def tagBits : Nat := 128

/-- **The verified object a pass name denotes.** The three deterministic passes enter at
error `0` with a one-point seed; phase folding is the one that consumes randomness. -/
noncomputable def passOf (cfg : SuperOptConfig) (murm : Murm) : PassName → RandPass
  | .DecomposeToffoli => DecomposeToffoliR
  | .DecomposeCz => DecomposeCzR
  | .CancelGates => CancelGatesR
  | .CnotMin => CnotMinR
  | .SuperOpt => SuperOptR cfg murm
  | .PhaseFoldRand => PhaseFoldRand

/-- The idealized pipeline applies the nonlinear transformation to its sampled circuit. -/
theorem passOf_phaseFoldRand_run (cfg : SuperOptConfig) (murm : Murm)
    (circuit : Circuit n m)
    (sample : (passOf cfg murm .PhaseFoldRand).Seed circuit) :
    (passOf cfg murm .PhaseFoldRand).run circuit sample =
      phaseFoldNonlinearWithSample 128 circuit sample := rfl

/-- The pipeline a level runs, when `--passes` is absent.

`CnotMin` leads the sweep: it re-synthesizes whole CNOT-dihedral blocks, reshaping the circuit
far more than the peephole rewriter does, and the passes after it work on the result. -/
def Level.pipeline : Level → List PassName
  | .O1 => [.CancelGates, .PhaseFoldRand]
  | _ => [.CnotMin, .CancelGates, .SuperOpt, .PhaseFoldRand]

/-- One round of the idealized randomized pipeline, in executable pass order. -/
noncomputable def tzapRound (cfg : SuperOptConfig) (murm : Murm) (names : List PassName) :
    RandPass := RandPass.pipeline (names.map (passOf cfg murm))

/-- The idealized randomized round, repeated while it keeps removing gates, at most `fuel`
times. -/
noncomputable def tzapRun (cfg : SuperOptConfig) (murm : Murm) (names : List PassName)
    (fuel : Nat) : RandPass := (tzapRound cfg murm names).fixpointShrink fuel

/-! ### What the run is worth

Two statements, and between them they say what the optimizer guarantees.

`tzapRun_correct` is the general one: the output denotes the same channel as the input except
on a set of seeds whose measure is at most `error`, which by `fixpointShrink_error_le` and
`pipeline_error_le` is at most the sum of the nonlinear comparison bounds across phase-fold
invocations. This needs no independence *between* rounds' failure events — a union bound never does —
only that each round's tags are drawn afresh, which is why `PhaseFoldRandExec` draws per call.

`tzapRun_exact` is the special one: drop `PhaseFoldRand` from `--passes` and the bound is
`0`, so *every* run is right, and the randomized machinery gives back exactly the
unconditional `Pass` guarantee. -/

theorem passOf_error_eq_zero {nm : PassName} (h : nm ≠ .PhaseFoldRand) (cfg : SuperOptConfig)
    (murm : Murm) (c : Circuit n m) : (passOf cfg murm nm).error c = 0 := by
  cases nm <;> simp_all [passOf, CancelGatesR, CnotMinR, SuperOptR, deterministicRand]

theorem tzapRound_error_eq_zero {names : List PassName} (h : PassName.PhaseFoldRand ∉ names)
    (cfg : SuperOptConfig) (murm : Murm) (c : Circuit n m) :
    (tzapRound cfg murm names).error c = 0 := by
  refine le_antisymm ?_ (by simp)
  have := RandPass.pipeline_error_le 0 (names.map (passOf cfg murm)) ?_ c
  · simpa [tzapRound] using this
  · intro p hp n' m' c
    obtain ⟨nm, hnm, rfl⟩ := List.mem_map.1 hp
    exact le_of_eq (passOf_error_eq_zero (by rintro rfl; exact h hnm) cfg murm c)

theorem tzapRun_error_eq_zero {names : List PassName} (h : PassName.PhaseFoldRand ∉ names)
    (cfg : SuperOptConfig) (murm : Murm) (fuel : Nat) (c : Circuit n m) :
    (tzapRun cfg murm names fuel).error c = 0 := by
  refine le_antisymm ?_ (by simp)
  have := RandPass.fixpointShrink_error_le (tzapRound cfg murm names) 0
    (fun c => le_of_eq (tzapRound_error_eq_zero h cfg murm c)) fuel c
  simpa [tzapRun] using this

/-- **The optimizer is correct.** For a well-formed circuit, the pipeline's output denotes the
same channel as its input, except on a set of seeds of measure at most `error`. -/
theorem tzapRun_correct (cfg : SuperOptConfig) (murm : Murm) (names : List PassName)
    (fuel : Nat) (c : Circuit n m) :
    ((tzapRun cfg murm names fuel).dist c).toOuterMeasure
        {s | ¬ ((tzapRun cfg murm names fuel).run c s).Equivalent c}
      ≤ (tzapRun cfg murm names fuel).error c :=
  (tzapRun cfg murm names fuel).correct c

/-- **…and exactly correct without the randomized pass.** -/
theorem tzapRun_exact {names : List PassName} (h : PassName.PhaseFoldRand ∉ names)
    (cfg : SuperOptConfig) (murm : Murm) (fuel : Nat) (c : Circuit n m)
    {s : (tzapRun cfg murm names fuel).Seed c}
    (hs : s ∈ ((tzapRun cfg murm names fuel).dist c).support) :
    ((tzapRun cfg murm names fuel).run c s).Equivalent c :=
  RandPass.correct_of_error_eq_zero _ c (tzapRun_error_eq_zero h cfg murm fuel c) hs

/-- **The run returns a circuit the back end may print**, for any seed: operands in range and
honest `has*` flags, from `RandPass`'s structural obligations. With `Qasm.parse_valid`, which
establishes the same of whatever the front end accepts, this holds from parse to emit. -/
theorem tzapRun_structural (cfg : SuperOptConfig) (murm : Murm) (names : List PassName)
    (fuel : Nat) (c : Circuit n m) (hc : c.raw.Structural)
    (s : (tzapRun cfg murm names fuel).Seed c) :
    ((tzapRun cfg murm names fuel).run c s).raw.Structural :=
  ⟨(tzapRun cfg murm names fuel).wellFormed_run c s hc.1,
   (tzapRun cfg murm names fuel).flagsOk_run c s hc.2⟩

/-! ## Randomized executable core -/

/-- The runtime interpretation of a pass name. Deterministic passes are lifted into `IO`;
phase folding obtains a fresh OS-random sample on every invocation. -/
def executableStep (cfg : SuperOptConfig) (murm : Murm) : PassName → ExecutableRandPass
  | .DecomposeToffoli => ExecutableRandPass.ofPass DecomposeToffoli
  | .DecomposeCz => ExecutableRandPass.ofPass DecomposeCz
  | .CancelGates => ExecutableRandPass.ofPass CancelGates
  | .CnotMin => ExecutableRandPass.ofPass CnotMin
  | .SuperOpt => ExecutableRandPass.ofPass (SuperOpt cfg murm)
  | .PhaseFoldRand => PhaseFoldRandExec

/-- Execution draws a sample, then applies the very same pure nonlinear pass as the
probability model. The theorem does not assert that OS bytes are uniformly distributed. -/
theorem executablePhaseFoldRand_run (cfg : SuperOptConfig) (murm : Murm)
    (circuit : Circuit n m) :
    (executableStep cfg murm .PhaseFoldRand).run circuit = do
      let sample ← randomSample (varBound circuit.raw) 128
      return (passOf cfg murm .PhaseFoldRand).run circuit sample := rfl

/-- How many fixpoint rounds a level allows: `O2` is the cheap bounded tier, the rest run out
fully. -/
def Level.maxRounds : Level → Option Nat
  | .O2 => some 2
  | _ => none

/-! ## Formatting -/

/-- Thousands separators, as Rust's `fmt_num`. -/
def fmtNum (n : Nat) : String :=
  let ds := (toString n).toList
  let grouped :=
    ds.reverse.foldl (fun (acc, i) c =>
      (if i ≠ 0 && i % 3 == 0 then c :: ',' :: acc else c :: acc, i + 1)) ([], 0) |>.1
  String.ofList grouped

/-- Seconds to three decimal places. -/
def fmtSecs (nanos : Nat) : String :=
  let ms := nanos / 1000000
  s!"{ms / 1000}.{String.ofList ((toString (ms % 1000)).toList.reverse.take 3 |>.reverse)
      |> fun x => (String.ofList (List.replicate (3 - x.length) '0')) ++ x}"

/-- A percentage reduction from `before` to `after`, one decimal place. -/
def fmtPct (before after : Nat) : String :=
  if before == 0 then "0.0"
  else
    let tenths := ((before - min before after) * 1000) / before
    s!"{tenths / 10}.{tenths % 10}"

/-! ## The run -/

/-- Force a `Nat` before reading the clock.

Lean's `let` is lazy, so a timing that brackets an unforced binding measures the cost of
allocating a thunk and nothing else. `IO.lazyPure` evaluates its thunk when the action runs,
and a `Nat` in weak head normal form is fully evaluated — so forcing a sum of counters forces
the work that produced them. (Branching on the value instead does *not* work: with both arms
equal, the compiler drops the test.) -/
def force (n : Unit → Nat) : IO Unit := do
  let _ ← IO.lazyPure n
  pure ()

/-- How many rounds a run may take.

`Level.maxRounds` when the level caps them; otherwise `gates + 1`, which is the whole loop:
a round that removes no gate ends it, so no more than `gates` rounds can continue. This is
the `fuel` `tzapRun` is indexed by. -/
def roundFuel (maxRounds : Option Nat) (c : RawCircuit) : Nat :=
  maxRounds.getD (c.gates.length + 1)

/-- The result of a run: the counts the banner compares. -/
structure Report where
  /-- The circuit as the pipeline received it. -/
  baseline : Metrics
  /-- The circuit the pipeline returned. -/
  output : Metrics
deriving Repr, Inhabited

/-- The checked randomized core used by the CLI. Conditional composition ensures that each
phase-folding invocation draws after seeing the circuit produced by the preceding passes and
rounds. -/
def runConfiguredChecked (cfg : SuperOptConfig) (murm : Murm)
    (c : Circuit n m) (o : Options) : IO (Circuit n m) := do
  let names := o.passes.getD o.level.pipeline
  let round := ExecutableRandPass.pipeline (names.map (executableStep cfg murm))
  if o.passes.isSome then
    if o.fixpoint then (round.fixpointShrink (roundFuel none c.raw)).run c
    else round.run c
  else if o.level.usesSuperOpt then
    (round.fixpointShrink (roundFuel o.level.maxRounds c.raw)).run c
  else if o.fixpoint then
    (round.fixpointShrink (roundFuel none c.raw)).run c
  else round.run c

/-- Raw randomized API boundary. Malformed internal circuits are left unchanged; parsed QASM
always takes the checked branch. -/
def runConfigured (cfg : SuperOptConfig) (murm : Murm)
    (c : RawCircuit) (o : Options) : IO RawCircuit := do
  if hc : c.Wf then
    return (← runConfiguredChecked cfg murm (Circuit.of c hc) o).raw
  else return c

/-- Run one verified deterministic pass at the raw API boundary. -/
def runPassRaw (pass : Pass) (c : RawCircuit) : RawCircuit :=
  if hc : c.Wf then (pass.run (Circuit.of c hc)).raw else c

/-- Load exactly the MURM needed by one stage.  An empty synthesis basis disables SuperOpt
for that stage instead of constructing a useless map. -/
def prepareMurm (c : RawCircuit) (o : Options) (forbidden : GateSet)
    (names : List PassName) : IO (SuperOptConfig × MurmConfig × Murm) := do
  let (cfg, murmCfg) := resolveBounds o c.gateSet forbidden
  let needsMurm := names.contains .SuperOpt && !murmCfg.basis.isEmpty
  let murm ← if needsMurm then do
      let cached ← MurmCache.isCached murmCfg
      if !cached then
        IO.eprintln "  🔧 Building MURM (one-time — cached for future use)..."
      let t0 ← IO.monoNanosNow
      let (murm, fromCache) ← MurmCache.loadOrBuild murmCfg
      let total := (List.range (murmCfg.maxQubits + 1)).foldl
        (fun acc k => acc + (murm.widths[k]?.map WidthMurm.size |>.getD 0)) 0
      force fun _ => total
      let t1 ← IO.monoNanosNow
      if o.verbose || !fromCache then
        let verb := if fromCache then "Loaded" else "Built"
        IO.eprintln s!"  {verb} MURM ({fmtNum total} unitaries) in \
                       {fmtSecs (t1 - t0)}s"
        IO.eprintln s!"    └─ Synthesis basis: {murmCfg.basis}"
        IO.eprintln ""
      pure murm
    else pure default
  return (cfg, murmCfg, murm)

/-- One default optimization stage, with its synthesis basis resolved from the circuit that
actually enters the stage. -/
def runDefaultStage (c : RawCircuit) (o : Options) (forbidden : GateSet) : IO RawCircuit := do
  let names := o.level.pipeline
  let (cfg, _, murm) ← prepareMurm c o forbidden names
  runConfigured cfg murm c { o with passes := none }

/-- Run an explicit pipeline once.  SuperOpt is handled at the pass boundary so `auto`
observes every preceding transformation, including an explicit decomposition. -/
def runExplicitSweep (c : RawCircuit) (o : Options) : List PassName → IO RawCircuit
  | [] => pure c
  | name :: names => do
      let (cfg, murmCfg, murm) ← prepareMurm c o GateSet.empty [name]
      let next ← if name == .SuperOpt && murmCfg.basis.isEmpty then pure c
        else runConfigured cfg murm c { o with passes := some [name], fixpoint := false }
      runExplicitSweep next o names

/-- Repeat an explicit pipeline while it strictly shrinks the circuit. -/
def runExplicitFixpoint : Nat → RawCircuit → Options → List PassName → IO RawCircuit
  | 0, c, _, _ => pure c
  | fuel + 1, c, o, names => do
      let next ← runExplicitSweep c o names
      if next.gates.length < c.gates.length then runExplicitFixpoint fuel next o names
      else pure next

/-- Run the optimizer. Loads the stage's MURM when the pipeline needs one, then
delegates the circuit transformation to the OS-randomized executable pipeline. -/
def optimize (c : RawCircuit) (o : Options) : IO (RawCircuit × Report) := do
  let baseline := Metrics.of c
  let result ← match o.passes with
    | some names =>
        if o.fixpoint then runExplicitFixpoint (roundFuel none c) c o names
        else runExplicitSweep c o names
    | none => do
        -- The default workflow deliberately optimizes before lowering native gates.
        let inputOptimized ← runDefaultStage c o GateSet.empty
        let mut decomposed := false
        let mut result := inputOptimized
        let ccxKinds := GateSet.ofKinds [.ccx, .ccz]
        if o.decomposeCcx && !(result.gateSet.inter ccxKinds).isEmpty then
          result := runPassRaw DecomposeToffoli result
          decomposed := true
        if o.decomposeCz && result.gateSet.contains .cz then
          result := runPassRaw DecomposeCz result
          decomposed := true
        if decomposed then runDefaultStage result o o.decomposedGateSet else pure result
  return (result, ⟨baseline, Metrics.of result⟩)

end TzapLean
