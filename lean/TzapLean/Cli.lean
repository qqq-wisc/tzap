import TzapLean.Optimize

/-!
# Command-line Interface

A port of `src/cli.rs` and `src/main.rs`: flag parsing, `--help`, and the run banner. Flags
keep Rust's names and meanings so a command line transfers between the two.

What is absent, and says so when asked for: `--parallel` (deliberately not ported),
`--decompose-rz` / `--decompose-cz` / `--epsilon` (no decomposition passes here), and the
`DecomposeToffoli` / `DecomposeCz` / `DecomposeRz` / `CliffordResynth` pass names. Rather than
report those as unknown flags — which would read as a typo — each gets an error saying it is
not in this build.

`--seed` is accepted for command-line compatibility; exact phase folding no longer consumes it.
-/

namespace TzapLean

namespace Cli

/-- Parsed options: the paths the CLI owns, plus the optimizer configuration. -/
structure Opts where
  /-- Input QASM path. -/
  inputPath : String
  /-- Output QASM path, if requested. -/
  outputPath : Option String := none
  /-- What to run. -/
  options : Options := {}
deriving Inhabited

/-- Print `Error: {msg}` and exit 1 — the single entry point for every argument failure, so
all CLI errors share one prefix. -/
def argError {α : Type} (msg : String) : IO α := do
  IO.eprintln s!"Error: {msg}"
  IO.Process.exit 1

/-- Parse an argument that must be a positive integer. -/
def parseUsize (args : List String) (i : Nat) (flag : String) : IO Nat := do
  let some raw := args[i]? | argError s!"{flag} requires a positive integer"
  let some v := raw.toNat? | argError s!"{flag} requires a positive integer"
  if v == 0 then argError s!"{flag} requires a positive integer, got 0" else pure v

/-- Flags Rust has that this build does not, each with the reason. -/
def unsupportedFlag (flag : String) : Option String :=
  match flag with
  | "--parallel" =>
      some "--parallel is not implemented in this build"
  | "--decompose-rz" | "--decompose-cz" | "--epsilon" =>
      some s!"{flag} is not available in this build — no decomposition passes are ported \
              (angles here are exact rationals in units of π, and gridsynth is not ported)"
  | _ => none

/-- Pass names Rust has that this build does not. -/
def unsupportedPass (name : String) : Bool :=
  ["DecomposeToffoli", "DecomposeCz", "DecomposeRz", "CliffordResynth"].contains name

/-- Parse a comma-separated pass list. -/
def parsePassList (list : String) : IO (List PassName) := do
  let names := (list.splitOn ",").map (·.trimAscii.toString) |>.filter (· ≠ "")
  if names.isEmpty then
    argError "--passes requires at least one pass name (e.g. --passes CancelGates,SuperOpt)"
  names.mapM fun name =>
    match PassName.parse name with
    | some p => pure p
    | none =>
        if unsupportedPass name then
          argError s!"Pass '{name}' is not ported to this build. \
                      Available passes: {PassName.allNames}"
        else
          argError s!"Unknown pass '{name}'. Available passes: {PassName.allNames}"

/-- Whether a token continues a `--passes` list rather than starting a new flag. -/
def looksLikePassFragment (token : String) : Bool :=
  let parts := (token.splitOn ",").map (·.trimAscii.toString) |>.filter (· ≠ "")
  !parts.isEmpty && parts.all fun n => (PassName.parse n).isSome

/-- The `--help` text. -/
def printHelp : IO Unit := do
  IO.println ""
  IO.println "  \x1b[1m⚡️ tzap-lean\x1b[0m  —  verified quantum circuit optimizer"
  IO.println ""
  IO.println "  \x1b[1;33mUSAGE\x1b[0m"
  IO.println "    tzap-lean <input.qasm> [output.qasm] [options]"
  IO.println ""
  IO.println "  Every executable pass carries an unconditional machine-checked proof"
  IO.println "  that its output denotes the same channel as its input."
  IO.println ""
  IO.println "  \x1b[1;33mARGS\x1b[0m"
  IO.println "    \x1b[1m<input.qasm>\x1b[0m     Input OpenQASM 2.0 file"
  IO.println "    \x1b[1m[output.qasm]\x1b[0m    Output file (no output if omitted)"
  IO.println ""
  IO.println "  \x1b[1;33mOPTIONS\x1b[0m"
  IO.println "    \x1b[1m-o\x1b[0m <file>        Write output to <file>"
  IO.println "    \x1b[1m--passes\x1b[0m <list>  Run these passes in order, overriding the default pipeline"
  IO.println "                     (see PASSES)"
  IO.println "    \x1b[1m--fixpoint\x1b[0m       Repeat the pipeline until gate count stops decreasing"
  IO.println "    \x1b[1m--seed\x1b[0m <n>       Accepted for compatibility (exact tags need no seed)"
  IO.println "    \x1b[1m--verbose\x1b[0m        Report detailed input and table-loading information"
  IO.println "    \x1b[1m-O1\x1b[0m              Fastest: phase folding + gate cancellation only"
  IO.println "    \x1b[1m-O2\x1b[0m              Adds a superoptimization pass to O1 (2 rounds)"
  IO.println "    \x1b[1m-O3\x1b[0m              Runs -O2 to a fixpoint (default)"
  IO.println "                     (slower first run: the table is built from scratch)"
  IO.println "    \x1b[1m-h, --help\x1b[0m       Print this help message"
  IO.println "    \x1b[1m-v, --version\x1b[0m    Print the version"
  IO.println ""
  IO.println "  \x1b[1;33mPASSES\x1b[0m (names for --passes)"
  for (name, p, desc) in PassName.all do
    let mark := if p.verified then "proved  " else "bounded "
    IO.println s!"    \x1b[1m{name}\x1b[0m{String.ofList (List.replicate (16 - name.length) ' ')}\x1b[2m{mark}\x1b[0m {desc}"
  IO.println ""
  IO.println "  \x1b[2mproved = output proved equivalent to input, unconditionally\x1b[0m"
  IO.println ""

/-- Parse the command line. -/
partial def parseArgs (args : List String) : IO Opts := do
  let mut inputPath : Option String := none
  let mut outputPath : Option String := none
  let mut passes : Option (List PassName) := none
  let mut fixpoint := false
  let mut level : Option Level := none
  let mut seed : Option Nat := none
  let mut verbose := false
  let mut soQubits : Option Nat := none
  let mut soWindow : Option Nat := none
  let mut soEntries : Option Nat := none
  let mut i := 0
  while h : i < args.length do
    let arg := args[i]
    match arg with
    | "--help" | "-h" => printHelp; IO.Process.exit 0
    | "--version" | "-v" => IO.println "tzap-lean 0.1.0"; IO.Process.exit 0
    | "--fixpoint" => fixpoint := true
    | "--verbose" => verbose := true
    | "--passes" =>
        i := i + 1
        let some first := args[i]? |
          argError "--passes requires a comma-separated list of pass names \
                    (e.g. --passes CancelGates,SuperOpt)"
        let mut list := first
        -- Rust allows a space-separated continuation, so `--passes A, B` works.
        while (args[i + 1]?.map looksLikePassFragment).getD false do
          list := list ++ "," ++ args[i + 1]!
          i := i + 1
        passes := some (← parsePassList list)
    | "--seed" =>
        i := i + 1
        seed := some (← parseUsize args i "--seed")
    | "-O1" | "-O2" | "-O3" =>
        if level.isSome then
          argError "-O1, -O2 and -O3 cannot be combined — pick exactly one"
        level := some (match arg with
          | "-O1" => .O1 | "-O2" => .O2 | _ => .O3)
    | "-o" =>
        i := i + 1
        let some p := args[i]? | argError "-o requires an output file path"
        outputPath := some p
    -- Hidden: not listed in --help, for experimenting with SuperOpt's bounds.
    | "--superopt-qubits" =>
        i := i + 1
        soQubits := some (← parseUsize args i "--superopt-qubits")
    | "--superopt-window-gates" =>
        i := i + 1
        soWindow := some (← parseUsize args i "--superopt-window-gates")
    | "--superopt-table-entries" =>
        i := i + 1
        soEntries := some (← parseUsize args i "--superopt-table-entries")
    | _ =>
        if let some reason := unsupportedFlag arg then
          argError reason
        else if arg.startsWith "-" then
          argError s!"unknown flag '{arg}'. Run `tzap-lean --help` for the list of valid options"
        else if inputPath.isNone then inputPath := some arg
        else if outputPath.isNone then outputPath := some arg
        else
          argError s!"unexpected extra argument '{arg}' — tzap-lean takes at most \
                      <input.qasm> and [output.qasm]"
    i := i + 1
  let some inPath := inputPath |
    argError "missing required <input.qasm> argument\n\n  \
              Usage: tzap-lean <input.qasm> [-o output.qasm] [-O1|-O2|-O3] \
              [--passes <list>] [--fixpoint]\n  \
              Run `tzap-lean --help` for the full option list."
  if level.isSome && (passes.isSome || fixpoint) then
    argError "-O1, -O2 and -O3 cannot be combined with --passes or --fixpoint"
  return { inputPath := inPath, outputPath,
           options := { level := level.getD .O3, passes, fixpoint, seed, verbose,
                        superopt := ⟨soQubits, soWindow, soEntries⟩ } }

/-! ## The run -/

/-- One row of the result banner. -/
def resultRow (label : String) (before after : Nat) : String :=
  let pad := String.ofList (List.replicate (8 - min 8 label.length) ' ')
  s!"    {label}{pad}  {fmtNum before} → {fmtNum after}  ({fmtPct before after}% reduction)"

/-- Read and parse the input, exiting with Rust's message shape on failure. -/
def readCircuit (verbose : Bool) (path : String) : IO Circuit := do
  let t0 ← IO.monoNanosNow
  let source ← try IO.FS.readFile path
    catch e => do
      IO.eprintln s!"\nError reading {path}: {e}"
      IO.Process.exit 1
  match Qasm.parse source with
  | .error e => do
      IO.eprintln s!"\nError parsing {path}: {e}"
      IO.Process.exit 1
  | .ok c => do
      force fun _ => c.gates.length + c.numQubits
      let t1 ← IO.monoNanosNow
      if verbose then
        IO.eprintln s!"  Parsed {path} in {fmtSecs (t1 - t0)}s"
        IO.eprintln s!"\t└─ {fmtNum c.numQubits} qubits · {fmtNum c.gates.length} gates"
        IO.eprintln ""
      return c

/-- The entry point. -/
def main (argv : List String) : IO UInt32 := do
  let t0 ← IO.monoNanosNow
  let opts ← parseArgs argv
  if opts.options.verbose then IO.eprintln "\x1b[1m⚡️ tzap-lean\x1b[0m"
  let circuit ← readCircuit opts.options.verbose opts.inputPath
  let (result, report) ← optimize circuit opts.options
  let t1 ← IO.monoNanosNow
  IO.eprintln s!"  \x1b[1mResult\x1b[0m  ({fmtSecs (t1 - t0)}s)"
  IO.eprintln (resultRow "Gates" report.baseline.gates report.output.gates)
  IO.eprintln (resultRow "2q gates" report.baseline.twoQubit report.output.twoQubit)
  IO.eprintln (resultRow "T/Tdg" report.baseline.t report.output.t)
  if report.baseline.rz > 0 || report.output.rz > 0 then
    IO.eprintln (resultRow "Rz" report.baseline.rz report.output.rz)
  IO.eprintln (resultRow "Depth" report.baseline.depth report.output.depth)
  IO.eprintln ""
  match opts.outputPath with
  | none => pure ()
  | some p => do
      try
        let source ← match Qasm.serializeChecked result with
          | .ok source => pure source
          | .error e => throw (IO.userError e)
        IO.FS.writeFile p source
        IO.eprintln s!"  wrote {p}"
      catch e => do
        IO.eprintln s!"Error writing {p}: {e}"
        IO.Process.exit 1
  return 0

end Cli

end TzapLean
