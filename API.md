# tzap API

## Circuits

A `Circuit` holds a list of gates over a fixed number of qubits.

```rust
use tzap::circuit::{Circuit, Gate};

let mut circuit = Circuit::new(2);
circuit.apply(Gate::h(0));
circuit.apply(Gate::cnot { control: 0, target: 1 });
circuit.apply(Gate::t(0));
```

To use `measure` gates, allocate classical bits with
`Circuit::with_cbits(num_qubits, num_cbits)` instead of `Circuit::new`.

### Supported gates

| Gate | Constructor |
|------|------------|
| X | `Gate::x(qubit)` |
| H | `Gate::h(qubit)` |
| S | `Gate::s(qubit)` |
| Sdg | `Gate::sdg(qubit)` |
| Z | `Gate::z(qubit)` |
| T | `Gate::t(qubit)` |
| Tdg | `Gate::tdg(qubit)` |
| Rz | `Gate::rz(angle, qubit)` |
| CNOT | `Gate::cnot { control, target }` |
| CZ | `Gate::cz { control, target }` |
| Toffoli | `Gate::ccx { control1, control2, target }` |
| CCZ | `Gate::ccz { control1, control2, target }` |
| Measure | `Gate::measure { qubit, cbit }` |
| Reset | `Gate::reset(qubit)` |

Operands are `Qubit` and `CBit`, both aliases for `u32` — a `Gate` is 16
bytes, which matters on circuits of tens of millions of gates. Integer
literals still work unannotated; a `usize` index needs `as u32`.

### QASM I/O

Parse from and convert to OpenQASM 2.0. `from_qasm` returns
`Result<Circuit, String>`:

```rust
use tzap::circuit::Circuit;

let circuit = Circuit::from_qasm("
    OPENQASM 2.0;
    include \"qelib1.inc\";
    qreg q[2];
    h q[0];
    cx q[0],q[1];
").expect("invalid QASM");

let qasm_string = circuit.to_qasm();
```

The QASM parser accepts `ccz` as a native circuit gate. `DecomposeToffoli`
lowers both `ccx` and `ccz` to Clifford+T.

## Optimizing

`tzap::optimize` runs the same staged workflow as the `tzap` CLI: optimize
the input circuit, apply requested decompositions, then optimize the changed
circuit again. There is no need to assemble one pass at a time to get the
CLI's results.

```rust,ignore
use tzap::circuit::Circuit;
use tzap::optimize::{Options, optimize};

let circuit = Circuit::from_qasm(qasm)?;
let (optimized, report) = optimize(&circuit, &Options::default())?;

println!(
    "{} → {} gates, {} → {} T",
    report.baseline.gates, report.output.gates,
    report.baseline.t, report.output.t,
);
# Ok::<(), Box<dyn std::error::Error>>(())
```

`Options::default()` is the CLI's default: `-O3`, sequential, and no
CCX/CCZ, CZ, or Rz decomposition.

| Field | Default | Meaning |
|-------|---------|---------|
| `level` | `Level::O3` | `O1` (cancel + phase-fold), `O2` (adds SuperOpt, 2 rounds), `O3` (same, to a fixpoint), `Osuper` (`O3` with the bigger SuperOpt bounds) |
| `passes` | `None` | An explicit `Vec<PassName>` pipeline, replacing `level`'s |
| `fixpoint` | `false` | Repeat until the gate count stops falling. Only consulted for pipelines that aren't already fixpoint loops (`passes`, or `O1`) |
| `decompose_rz` | `false` | Decompose Rz into Clifford+T via gridsynth |
| `decompose_cz` | `false` | Decompose CZ into H+CX+H before optimizing |
| `decompose_ccx` | `false` | Decompose CCX and CCZ into Clifford+T |
| `rz_epsilon` | `1e-10` | Approximation epsilon for `decompose_rz` |
| `parallel` | `false` | Optimize gate-contiguous chunks concurrently, then concatenate |
| `superopt` | all `None` | Per-run overrides for the SuperOpt window/MURM bounds |
| `superopt_gates` | `Auto` | Base MURM basis plus native CZ/CCX/CCZ gates present in the current stage |

`Report` carries three sets of `Metrics` (`gates`, `two_qubit`, `depth`, `t`,
`rz`): `input` as handed in, `baseline` (the same original input comparison
point), and `output`. Requested decompositions are opt-in middle stages:
optimize the input circuit, decompose CCX/CCZ then CZ then Rz, and optimize
again when a decomposition changed the circuit.

The post-decomposition MURM basis excludes every gate family requested for
decomposition, including gates named in an explicit `superopt_gates` basis.
Within an explicit `passes` pipeline, `Auto` is resolved independently at each
SuperOpt occurrence from the circuit produced by the preceding passes.

### Reporting progress

`optimize` is silent. To report progress, implement `Observer` (every method
defaults to doing nothing) and call `optimize_with`:

```rust,ignore
use tzap::circuit::Circuit;
use tzap::optimize::{Metrics, Observer, Options, optimize_with};

struct Log;

impl Observer for Log {
    fn progress_update(&self, round: Option<usize>, current: &Circuit, baseline: Metrics) {
        eprintln!("round {round:?}: {} → {} gates", baseline.gates, current.gates.len());
    }
}

let (optimized, _) = optimize_with(&circuit, &Options::default(), &Log)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

Events fire from whichever thread reaches them, so an `Observer` must be
`Sync`; under `parallel`, `chunk_done` is called concurrently from rayon
workers. The chunk workers' own pipelines are always observed by `Silent`,
since their events would otherwise interleave. Set `tracks_chunks` to `true`
to receive the `chunks_start`/`chunk_done`/`chunks_end` events — they're
skipped by default, along with the whole-circuit stitch needed to compute
their metrics.

## Passes

The passes below are the building blocks `tzap::optimize` composes. Reach for
them directly to build a pipeline it doesn't offer.

Every pass implements the `Pass` trait:

```rust,ignore
use tzap::pass::Pass;

pub trait Pass {
    fn name(&self) -> &str;
    fn run(&self, circuit: &Circuit) -> Circuit;
}
```

A custom pass only needs to supply `name` and `run`.

### Available passes

| Pass | Import | Description |
|------|--------|-------------|
| `DecomposeToffoli` | `tzap::decompose` | Breaks CCX and CCZ gates into Clifford+T |
| `DecomposeCz` | `tzap::decompose` | Explicitly lowers CZ gates to H+CX+H |
| `DecomposeRz` | `tzap::decompose` | Decomposes Rz gates into Clifford+T via gridsynth |
| `CancelGates` | `tzap::cancel` | Removes adjacent self-inverse gate pairs (HH, XX, etc.) |
| `SuperOpt` | `tzap::super_opt` | Replaces small windows using its shared MURM |
| `PhaseFoldRand` | `tzap::phase_fold_rand` | Merges T/Rz gates across the circuit via randomized parity tracking |
| `CnotMin` | `tzap::cnot_min` | Re-synthesizes CNOT-dihedral blocks to cut two-qubit gates |

### Running passes

Run a single pass:

```rust,ignore
use tzap::decompose::DecomposeToffoli;

let optimized = DecomposeToffoli.run(&circuit);
```

Run a pipeline:

```rust,ignore
use tzap::decompose::DecomposeToffoli;
use tzap::cancel::CancelGates;
use tzap::phase_fold_rand::PhaseFoldRand;
use tzap::pass::{Pass, PassResult, run_passes, count_t};

let passes: Vec<&dyn Pass> = vec![
    &DecomposeToffoli,
    &CancelGates,
    &PhaseFoldRand,
];

let result: PassResult = run_passes(&circuit, &passes);
println!("{} gates, {} T", result.circuit.gates.len(), count_t(&result.circuit));
```

`run_passes` returns a `PassResult`:

```rust,ignore
pub struct PassResult {
    pub circuit: Circuit,
    pub t_after_first: usize,       // T-count after only the first pass
    pub gates_after_first: usize,   // gate count after only the first pass
}
```

The `t_after_first` / `gates_after_first` fields are useful for
attributing reductions to the leading decomposition pass when reporting
end-to-end numbers. Helpers `count_t` and `count_rz` are also exposed
from `tzap::pass`.

### SuperOpt

`SuperOpt` is a peephole pass. It scans each maximal connected subcircuit window
and replaces it with the smallest equivalent circuit from a precomputed
MURM, applying a rewrite only when it strictly reduces the
gate count. Every replacement is verified by matrix equality up to global phase
before use, so rewrites are always semantics-preserving. Matrices use exact
Clifford+T arithmetic; Rz gates act as window barriers and are left unchanged.
The pass accepts unitary circuits only.

```rust,ignore
use tzap::super_opt::{SuperOpt, MurmConfig};

let pass = SuperOpt::new(3, 10, MurmConfig::default())?;
let result = pass.run(&circuit)?;
println!("{} rewrites", result.rewrites.len());
# Ok::<(), Box<dyn std::error::Error>>(())
```

Parameters:

- `max_qubits` — maximum distinct qubits in a scanned window.
- `window_gates` — maximum gates in a scanned window.
- `MurmConfig::new(max_qubits, max_gates, max_entries_per_qubit)` — bounds
  for the MURM, independent of the window size; `default()` is
  `(3, 8, 200_000)`. A MURM entry can only ever be used when it's strictly
  smaller than the window it would replace, so `max_gates` never needs to exceed
  `window_gates - 1`.

For a materially more thorough (but slower to build) configuration — the CLI's
`-Osuper` uses exactly this — try:

```rust,ignore
use tzap::super_opt::{SuperOpt, MurmConfig};

let pass = SuperOpt::new(5, 30, MurmConfig::new(5, 29, 5_000_000))?;
# Ok::<(), tzap::super_opt::SuperOptError>(())
```

**MURM construction and caching.** Building the MURM is the
expensive part — breadth-first enumeration over the gate library, bounded by
`max_gates` and `max_entries_per_qubit`. The exact basis is part of
`MurmConfig` and the cache identity. MURMs are cached two ways:

1. **Per-process, in-memory.** Every `SuperOpt::new` call with the same
   `MurmConfig` shares one already-built MURM for the life of the
   process (`Arc`-backed, keyed by config).
2. **On disk, across processes.** The built MURM is also persisted under
   `<cache root>/murm/` (one file per distinct config), so a later
   process with the same config loads it in well under a second instead of
   rebuilding it. The reader bounds every stored length, validates roots,
   parent ordering, gate operands/basis, metadata, exact EOF, and a full-body
   checksum before accepting a MURM. A missing, stale, or corrupt cache file
   is never a hard error — it triggers a fresh build, which then gets cached
   for next time. Call `tzap::super_opt::murm_is_cached(config)` to check up
   front whether a matching cache candidate exists (useful for deciding
   whether the next `SuperOpt::new` is likely to be slow); the subsequent load
   remains authoritative and may rebuild a candidate whose body is invalid.

   The cache root follows the XDG Base Directory Specification:
   `$XDG_CACHE_HOME/tzap`, falling back to `$HOME/.cache/tzap`, with
   `$TZAP_CACHE_DIR` overriding both. A native Windows process has none of
   those, so `%LOCALAPPDATA%\tzap` and then `%USERPROFILE%\.cache\tzap` are
   tried after them. `tzap::super_opt::set_cache_dir(dir)`
   overrides it for the process (call it once, before any MURM is built);
   `cache_dir()` reports the location in force, `cache_entries()` lists what
   is cached, and `clear_cache()` deletes it. Pre-MURM synthesis tables under
   `~/.tzap/superopt-tables/` use an incompatible format and are not read as
   MURMs. `clear_cache()` only removes current MURM files; obsolete tables may
   be deleted manually.

`SuperOpt` also implements `tzap::pass::Pass`. Chain `.without_subcircuits()` when
only the optimized circuit is needed, to skip retaining per-window diagnostics.
Chain `.incremental()` when repeatedly re-running the same pass instance on
successive versions of one evolving circuit (e.g. inside a fixpoint loop) — it
anchors new windows only near what changed since the previous `run` call,
which is unsound if the instance ever sees unrelated circuits or concurrent
chunks, so don't share an incremental instance across parallel workers.

### DecomposeRz epsilon

Control the approximation precision with the `epsilon` field (default `1e-10`):

```rust,ignore
use tzap::decompose::DecomposeRz;

let pass = DecomposeRz { epsilon: 1e-6 };
let cliffordt = pass.run(&circuit);
```
