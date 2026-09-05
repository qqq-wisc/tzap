import TzapLean.SuperOptProof
import TzapLean.PhaseFoldRand
import TzapLean.Qasm
import TzapLean.TableCache
import TzapLean.Pipeline

/-!
# The Optimizer Driver

A port of the parts of `src/optimize.rs` this development can support: the optimization
levels, the pass names selectable by `--passes`, the fixpoint loop, and the metrics the CLI
reports.

All four executable passes are `Pass`es and carry unconditional proofs.  Phase folding uses a
collision-free affine encoding; the fixed-width randomized construction remains in
`PhaseFoldRand.lean` as a separately proved algorithm, but is not in the executable's trusted
path. `runPipeline`, `runToFixpoint`, and `runConfigured` are pure executable definitions with
correctness theorems, so there is no hand-matched IO loop.

Not ported: `DecomposeToffoli`, `DecomposeCz`, `DecomposeRz` (gridsynth), `CliffordResynth`,
and parallel chunking. `Level.O3` is therefore `O2` run to a true fixpoint, without Rust's
one-shot Clifford re-synthesis at the end.
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
def Metrics.of (c : Circuit) : Metrics := Id.run do
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
  /-- Exact phase folding + gate cancellation. Fastest. -/
  | O1
  /-- Adds `SuperOpt`, capped at two rounds. -/
  | O2
  /-- Like `O2`, run to a true fixpoint. The default. -/
  | O3
deriving Repr, DecidableEq, Inhabited

/-- A pass selectable by name in `--passes`. -/
inductive PassName where
  /-- Cancel adjacent self-inverse pairs and reduce Hadamards. -/
  | CancelGates
  /-- Re-synthesize CNOT-dihedral blocks. -/
  | CnotMin
  /-- Peephole superoptimization against the synthesis table. -/
  | SuperOpt
  /-- Collision-free phase folding (the historical CLI name is retained). -/
  | PhaseFoldRand
deriving Repr, DecidableEq, Inhabited

namespace PassName

/-- All passes — name, variant, and what it carries — in a stable order for listing. -/
def all : List (String × PassName × String) :=
  [ ("CancelGates", .CancelGates,
     "Cancel adjacent self-inverse gate pairs and reduce Hadamards"),
    ("CnotMin", .CnotMin,
     "Re-synthesize CNOT-dihedral blocks to cut CNOT count"),
    ("SuperOpt", .SuperOpt,
     "Peephole superoptimization against the exact synthesis table"),
    ("PhaseFoldRand", .PhaseFoldRand,
     "Merge rotations on the same parity (collision-free exact tags)") ]

/-- Parse a pass name. -/
def parse (s : String) : Option PassName :=
  (all.find? (·.1 == s)).map (·.2.1)

/-- Every pass name, comma-separated — for error messages. -/
def allNames : String := String.intercalate ", " (all.map (·.1))

/-- Every executable pass carries an unconditional correctness proof. -/
def verified (_ : PassName) : Bool := true

end PassName

/-! ## Options -/

/-- `SuperOpt` window and table bounds, `none` meaning "whatever the level implies". -/
structure SuperOptBounds where
  /-- Widest window and table, in wires. -/
  qubits : Option Nat := none
  /-- Longest window, in gates. -/
  windowGates : Option Nat := none
  /-- Cap on stored unitaries per table width. -/
  tableEntries : Option Nat := none
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
  /-- Compatibility value accepted from `--seed`; the exact executable pass ignores it. -/
  seed : Option Nat := none
  /-- Print detailed input and synthesis-table loading information. -/
  verbose : Bool := false
deriving Repr, Inhabited

/-- The window/table bounds a level implies: Rust's own — 3 wires, 25-gate windows, a
200,000-entry table, with the same `table_gates = window_gates - 1` mapping.

That table takes about 76 seconds to build here against Rust's parallel builder, which is
affordable only because it is built once and cached (`TableCache`): a warm run loads its
549,456 unitaries in 0.07 s. `--superopt-qubits`, `--superopt-window-gates` and
`--superopt-table-entries` override any of the three. -/
def Level.bounds : Level → Nat × Nat × Nat
  | _ => (3, 25, 200000)

/-- Resolve the bounds for a run: level preset, then any explicit override. -/
def resolveBounds (o : Options) : SuperOptConfig × SuperOptTableConfig :=
  let (q, w, e) := o.level.bounds
  let q := o.superopt.qubits.getD q
  let w := o.superopt.windowGates.getD w
  let e := o.superopt.tableEntries.getD e
  -- A table entry only ever replaces a strictly larger window, so the table never needs to
  -- be deeper than `windowGates - 1`.
  ({ maxQubits := q, maxWindow := w }, { maxQubits := q, maxGates := w - 1,
                                         maxEntriesPerQubit := e })

/-- Whether a level's pipeline includes `SuperOpt`, and so pays for a table. -/
def Level.usesSuperOpt : Level → Bool
  | .O1 => false
  | _ => true

/-! ## Idealized randomized model

This section retains the fixed-width probabilistic pipeline and its theorem for library callers.
It is not the executable pipeline; the executable core below selects `PhaseFoldExact`. -/

/-- Tag width for phase folding: 63 bits, so one round misleads the pass with probability at
most `C(t,2)·2⁻⁶³`. -/
def tagBits : Nat := 63

/-- **The verified object a pass name denotes.** The three deterministic passes enter at
error `0` with a one-point seed; phase folding is the one that consumes randomness. -/
noncomputable def passOf (cfg : SuperOptConfig) (tbl : SynthTable) : PassName → RandPass
  | .CancelGates => CancelGatesR
  | .CnotMin => CnotMinR
  | .SuperOpt => SuperOptR cfg tbl
  | .PhaseFoldRand => PhaseFoldRand tagBits

/-- The pipeline a level runs, when `--passes` is absent.

`CnotMin` leads the sweep: it re-synthesizes whole CNOT-dihedral blocks, reshaping the circuit
far more than the peephole rewriter does, and the passes after it work on the result. -/
def Level.pipeline : Level → List PassName
  | .O1 => [.CancelGates, .PhaseFoldRand]
  | _ => [.CnotMin, .CancelGates, .SuperOpt, .PhaseFoldRand]

/-- One round of the idealized randomized pipeline, in executable pass order. -/
noncomputable def tzapRound (cfg : SuperOptConfig) (tbl : SynthTable) (names : List PassName) :
    RandPass := RandPass.pipeline (names.map (passOf cfg tbl))

/-- The idealized randomized round, repeated while it keeps removing gates, at most `fuel`
times. -/
noncomputable def tzapRun (cfg : SuperOptConfig) (tbl : SynthTable) (names : List PassName)
    (fuel : Nat) : RandPass := (tzapRound cfg tbl names).fixpointShrink fuel

/-! ### What the run is worth

Two statements, and between them they say what the optimizer guarantees.

`tzapRun_correct` is the general one: the output denotes the same channel as the input except
on a set of seeds whose measure is at most `error`, which by `fixpointShrink_error_le` and
`pipeline_error_le` is at most (rounds × passes) times one phase fold's `C(t,2)·2⁻ᵏ`. Note
that this needs no independence *between* rounds' failure events — a union bound never does —
only that each round's tags are drawn afresh, which is why `phaseFoldIO` draws per call.

`tzapRun_exact` is the special one: drop `PhaseFoldRand` from `--passes` and the bound is
`0`, so *every* run is right, and the randomized machinery gives back exactly the
unconditional `Pass` guarantee. -/

theorem passOf_error_eq_zero {nm : PassName} (h : nm ≠ .PhaseFoldRand) (cfg : SuperOptConfig)
    (tbl : SynthTable) (c : Circuit) : (passOf cfg tbl nm).error c = 0 := by
  cases nm <;> simp_all [passOf, CancelGatesR, CnotMinR, SuperOptR, deterministicRand]

theorem tzapRound_error_eq_zero {names : List PassName} (h : PassName.PhaseFoldRand ∉ names)
    (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit) :
    (tzapRound cfg tbl names).error c = 0 := by
  refine le_antisymm ?_ (by simp)
  have := RandPass.pipeline_error_le 0 (names.map (passOf cfg tbl)) ?_ c
  · simpa [tzapRound] using this
  · intro p hp c
    obtain ⟨nm, hnm, rfl⟩ := List.mem_map.1 hp
    exact le_of_eq (passOf_error_eq_zero (by rintro rfl; exact h hnm) cfg tbl c)

theorem tzapRun_error_eq_zero {names : List PassName} (h : PassName.PhaseFoldRand ∉ names)
    (cfg : SuperOptConfig) (tbl : SynthTable) (fuel : Nat) (c : Circuit) :
    (tzapRun cfg tbl names fuel).error c = 0 := by
  refine le_antisymm ?_ (by simp)
  have := RandPass.fixpointShrink_error_le (tzapRound cfg tbl names) 0
    (fun c => le_of_eq (tzapRound_error_eq_zero h cfg tbl c)) fuel c
  simpa [tzapRun] using this

/-- **The optimizer is correct.** For a well-formed circuit, the pipeline's output denotes the
same channel as its input, except on a set of seeds of measure at most `error`. -/
theorem tzapRun_correct (cfg : SuperOptConfig) (tbl : SynthTable) (names : List PassName)
    (fuel : Nat) (c : Circuit) (hc : c.Wf) :
    ((tzapRun cfg tbl names fuel).dist c).toOuterMeasure
        {s | ¬ Equivalent c.numQubits c.numCbits
          ((tzapRun cfg tbl names fuel).run c s).gates c.gates}
      ≤ (tzapRun cfg tbl names fuel).error c :=
  (tzapRun cfg tbl names fuel).correct c hc

/-- **…and exactly correct without the randomized pass.** -/
theorem tzapRun_exact {names : List PassName} (h : PassName.PhaseFoldRand ∉ names)
    (cfg : SuperOptConfig) (tbl : SynthTable) (fuel : Nat) (c : Circuit) (hc : c.Wf)
    {s : (tzapRun cfg tbl names fuel).Seed c}
    (hs : s ∈ ((tzapRun cfg tbl names fuel).dist c).support) :
    Equivalent c.numQubits c.numCbits ((tzapRun cfg tbl names fuel).run c s).gates c.gates :=
  RandPass.correct_of_error_eq_zero _ c hc (tzapRun_error_eq_zero h cfg tbl fuel c) hs

/-- **The run returns a circuit the back end may print**, for any seed: operands in range and
honest `has*` flags, from `RandPass`'s structural obligations. With `Qasm.parse_valid`, which
establishes the same of whatever the front end accepts, this holds from parse to emit. -/
theorem tzapRun_structural (cfg : SuperOptConfig) (tbl : SynthTable) (names : List PassName)
    (fuel : Nat) (c : Circuit) (hwf : c.Wf) (hc : c.Structural)
    (s : (tzapRun cfg tbl names fuel).Seed c) :
    ((tzapRun cfg tbl names fuel).run c s).Structural :=
  ⟨(tzapRun cfg tbl names fuel).wellFormed_run c s hwf hc.1,
   (tzapRun cfg tbl names fuel).flagsOk_run c s hc.2⟩

/-! ## Verified executable core -/

/-- The verified pass behind an executable step. -/
def verifiedStep (cfg : SuperOptConfig) (tbl : SynthTable) : PassName → Pass
  | .CancelGates => CancelGates
  | .CnotMin => CnotMin
  | .SuperOpt => SuperOpt cfg tbl
  | .PhaseFoldRand => PhaseFoldExact

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

/-- Run every pass once, in order. The checked type makes this pure function both the
executable and the object of its correctness theorem. -/
def runPipeline (passes : List Pass) (c : Circuit.Checked n m) : Circuit.Checked n m :=
  Pass.runAll passes c

/-- How many rounds a run may take.

`Level.maxRounds` when the level caps them; otherwise `gates + 1`, which is the whole loop:
a round that removes no gate ends it, so no more than `gates` rounds can continue. This is
the `fuel` `tzapRun` is indexed by. -/
def roundFuel (maxRounds : Option Nat) (c : Circuit) : Nat :=
  maxRounds.getD (c.gates.length + 1)

/-- Repeat the pipeline while it keeps removing gates, at most the supplied number of rounds. -/
def runToFixpoint (passes : List Pass) : Nat → Circuit.Checked n m → Circuit.Checked n m
  | 0, c => c
  | fuel + 1, c =>
      let out := runPipeline passes c
      if out.raw.gates.length < c.raw.gates.length then runToFixpoint passes fuel out else out

theorem runToFixpoint_correct (passes : List Pass) : ∀ fuel (c : Circuit.Checked n m),
    (runToFixpoint passes fuel c).Equivalent c := by
  intro fuel
  induction fuel with
  | zero => intro c; exact Equivalent.refl _ _ _
  | succ fuel ih =>
      intro c
      have hround := Pass.correct_runAll passes c
      simp only [runToFixpoint]
      split
      · exact Equivalent.trans (ih (runPipeline passes c)) hround
      · exact hround

/-- The result of a run: the counts the banner compares. -/
structure Report where
  /-- The circuit as the pipeline received it. -/
  baseline : Metrics
  /-- The circuit the pipeline returned. -/
  output : Metrics
deriving Repr, Inhabited

/-- The pure checked optimization core used by the IO shell. -/
def runConfiguredChecked (cfg : SuperOptConfig) (tbl : SynthTable)
    (c : Circuit.Checked n m) (o : Options) : Circuit.Checked n m :=
  let names := o.passes.getD o.level.pipeline
  let passes := names.map (verifiedStep cfg tbl)
  if o.passes.isSome then
    if o.fixpoint then runToFixpoint passes (roundFuel none c.raw) c
    else runPipeline passes c
  else if o.level.usesSuperOpt then
    runToFixpoint passes (roundFuel o.level.maxRounds c.raw) c
  else if o.fixpoint then runToFixpoint passes (roundFuel none c.raw) c
  else runPipeline passes c

/-- The raw API boundary. Malformed internal circuits are left unchanged; parsed QASM always
takes the checked branch. -/
def runConfigured (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit) (o : Options) : Circuit :=
  if hc : c.Wf then (runConfiguredChecked cfg tbl (Circuit.Checked.of c hc) o).raw else c

/-- Correctness of the checked optimization core. -/
theorem runConfiguredChecked_correct (cfg : SuperOptConfig) (tbl : SynthTable)
    (c : Circuit.Checked n m) (o : Options) :
    (runConfiguredChecked cfg tbl c o).Equivalent c := by
  unfold runConfiguredChecked
  split <;> split
  · exact runToFixpoint_correct _ _ _
  · exact Pass.correct_runAll _ _
  · exact runToFixpoint_correct _ _ _
  · split
    · exact runToFixpoint_correct _ _ _
    · exact Pass.correct_runAll _ _

/-- **End-to-end correctness of the executable optimization core.** -/
theorem runConfigured_correct (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit)
    (o : Options) (hc : c.Wf) :
    Equivalent c.numQubits c.numCbits (runConfigured cfg tbl c o).gates c.gates := by
  rw [runConfigured, dif_pos hc]
  exact runConfiguredChecked_correct cfg tbl (Circuit.Checked.of c hc) o

/-- **End-to-end checked-output theorem.** If the front end accepts `input` and the checked
serializer accepts the optimized circuit, then the emitted source parses back to that exact
circuit (apart from rebuilt cache flags), and its gates denote the same channel as the input. -/
theorem runConfigured_checkedOutput_correct (cfg : SuperOptConfig) (tbl : SynthTable)
    (o : Options) {input output : String} {c : Circuit}
    (hinput : Qasm.parse input = .ok c)
    (houtput : Qasm.serializeChecked (runConfigured cfg tbl c o) = .ok output) :
    Qasm.parse output = .ok ((runConfigured cfg tbl c o).withGates
      (runConfigured cfg tbl c o).gates) ∧
    Equivalent c.numQubits c.numCbits (runConfigured cfg tbl c o).gates c.gates :=
  ⟨(Qasm.serializeChecked_sound houtput).2,
   runConfigured_correct cfg tbl c o (Qasm.parse_wf hinput)⟩

/-- Run the optimizer. Builds the synthesis table first when the pipeline needs one, then
delegates the circuit transformation exactly to `runConfigured`. -/
def optimize (c : Circuit) (o : Options) : IO (Circuit × Report) := do
  let names := o.passes.getD o.level.pipeline
  let (cfg, tcfg) := resolveBounds o
  -- Only pay for a table if some selected pass will consult it.
  let needsTable := names.contains .SuperOpt
  let tbl ← if needsTable then do
      -- Captured before the load below can create the file, so a cold run says so.
      let cached ← TableCache.isCached tcfg
      -- A cold build takes tens of seconds, so it is announced even quietly; a warm load is
      -- fast enough to stay silent unless asked for.
      if !cached then
        IO.eprintln "  🔧 Building superoptimizer table (one-time — cached for future use)..."
      let t0 ← IO.monoNanosNow
      let (tbl, fromCache) ← TableCache.loadOrBuild tcfg
      let total := (List.range (tcfg.maxQubits + 1)).foldl
        (fun acc k => acc + (tbl.widths[k]?.map WidthTable.size |>.getD 0)) 0
      force fun _ => total
      let t1 ← IO.monoNanosNow
      if o.verbose || !fromCache then
        let verb := if fromCache then "Loaded" else "Built"
        IO.eprintln s!"  {verb} superoptimizer table ({fmtNum total} unitaries) in \
                       {fmtSecs (t1 - t0)}s"
        IO.eprintln ""
      pure tbl
    else pure default
  let baseline := Metrics.of c
  let result := runConfigured cfg tbl c o
  return (result, ⟨baseline, Metrics.of result⟩)

end TzapLean
