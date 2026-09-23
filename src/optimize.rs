//! The optimization driver: optimization levels, pass pipelines, the fixpoint
//! loop, `SuperOpt` construction, and map-reduce parallelism.
//!
//! This is everything the `tzap` binary used to own itself, so that every
//! frontend — the CLI, a Rust caller, a language binding — runs the exact same
//! pipeline rather than reimplementing `-O3`. The CLI keeps only argument
//! parsing, file I/O, and terminal rendering; it plugs the latter in through
//! [`Observer`].

use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::cancel::CancelGates;
use crate::circuit::{Circuit, Gate, GateKind, GateSet, qubit_operands};
use crate::cnot_min::CnotMin;
use crate::decompose::{DecomposeCz, DecomposeRz, DecomposeToffoli};
use crate::pass::Pass;
use crate::phase_fold_rand::PhaseFoldRand;
#[cfg(test)]
use crate::super_opt::OPTIONAL_GATE_KINDS;
use crate::super_opt::{
    BASE_GATE_SET, MurmConfig, OPTIONAL_GATE_SET, SUPPORTED_GATE_SET, SuperOpt, SuperOptError,
    murm_is_cached,
};

/// Map-reduce chunks per logical core. Deliberately more than one thread per
/// core (see [`num_threads`]): chunks cost varies (some hit more SuperOpt
/// rewrites than others), so more chunks than threads lets rayon's
/// work-stealing load-balance that unevenness across a right-sized pool.
const CHUNK_MULTIPLIER: usize = 2;

/// Default approximation epsilon for [`Options::rz_epsilon`].
pub const DEFAULT_RZ_EPSILON: f64 = 1e-10;

/// Default SuperOpt window/MURM bounds, overridable per-run via
/// [`SuperOptBounds`].
///
/// The window and MURM share both a qubit bound and a gate-count bound: the
/// [`SuperOpt`] pass itself allows window and MURM bounds to differ on either
/// axis (e.g. a window wider or deeper than the MURM backing it, to exercise
/// window mechanics beyond what the MURM can synthesize replacements for —
/// see `super_opt::tests`), but the driver has no everyday use case for that,
/// so it only exposes one knob per axis.
///
/// `window_gates=10` leaves real T-count on the table suite-wide; the T
/// floor is reached by `window_gates≈15` and gate-count keeps improving
/// slowly beyond that, so 25 is used as a deliberately more thorough
/// default. `qubits` and `murm_entries` showed no benefit worth their
/// added cost at this tier and were left alone.
pub const DEFAULT_SUPEROPT_QUBITS: usize = 3;
/// See [`DEFAULT_SUPEROPT_QUBITS`].
pub const DEFAULT_SUPEROPT_WINDOW_GATES: usize = 25;
/// See [`DEFAULT_SUPEROPT_QUBITS`].
pub const DEFAULT_SUPEROPT_MURM_ENTRIES: usize = 200_000;

/// SuperOpt bounds for [`Level::Osuper`]: a materially bigger window/MURM
/// than the default. Confirmed (by direct comparison against
/// `DEFAULT_SUPEROPT_*` across the full feynman+cobble benchmark suite) to be
/// a real, zero-regression improvement — concentrated in circuits with long
/// single-qubit runs — at the cost of a slower one-time MURM build (still
/// cached to disk after the first run). `window_gates` was swept 15→20→30,
/// each step a further zero-or-near-zero-regression win (feynman gains
/// saturated by 20; cobble kept improving through 30 with zero regressions).
/// Two other axes stopped helping, though: bigger qubit widths hit an
/// out-of-memory wall during MURM construction, and a bigger entries cap
/// starts *regressing* output (a bigger MURM is a strict fingerprint
/// superset, but the greedy, non-backtracking rewrite-selection rule doesn't
/// turn "more matches available" into "better output" monotonically).
/// `window_gates=40` is the best gate-count point found at this
/// qubit/entries setting; T-count is unchanged from `window_gates=30`.
pub const SUPER_SUPEROPT_QUBITS: usize = 5;
/// See [`SUPER_SUPEROPT_QUBITS`].
pub const SUPER_SUPEROPT_WINDOW_GATES: usize = 40;
/// See [`SUPER_SUPEROPT_QUBITS`].
pub const SUPER_SUPEROPT_MURM_ENTRIES: usize = 5_000_000;

/// An optimization level: which default pipeline [`optimize`] runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    /// Randomized phase folding + gate cancellation. Fastest.
    O1,
    /// Adds a SuperOpt pass to [`Level::O1`], capped at 2 rounds rather than
    /// run to a true fixpoint — see `optimize_default`'s `max_rounds`. With
    /// `decompose_rz` the cap allows one extra round, so Rz synthesis lands in
    /// the same place it does at the uncapped levels (see [`run_to_fixpoint`]).
    O2,
    /// Like [`Level::O2`], but run to a true fixpoint instead of capped at 2
    /// rounds. The default.
    O3,
    /// Like [`Level::O3`], but with the `SUPER_SUPEROPT_*` bounds.
    Osuper,
}

/// A pass selectable by name in [`Options::passes`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PassName {
    DecomposeToffoli,
    DecomposeCz,
    DecomposeRz,
    CancelGates,
    SuperOpt,
    PhaseFoldRand,
    CnotMin,
}

/// A structured top-level stage in an optimization run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StageKind {
    ExplicitPipeline,
    InputOptimization,
    DecomposeCcx,
    DecomposeCz,
    DecomposeRz,
    PostDecompositionOptimization,
}

impl StageKind {
    pub const fn name(self) -> &'static str {
        match self {
            Self::ExplicitPipeline => "Running explicit pipeline",
            Self::InputOptimization => "Optimizing input circuit",
            Self::DecomposeCcx => "Decomposing CCX/CCZ",
            Self::DecomposeCz => "Decomposing CZ",
            Self::DecomposeRz => "Decomposing Rz",
            Self::PostDecompositionOptimization => "Optimizing decomposed circuit",
        }
    }
}

impl PassName {
    /// All passes — `(name, variant, description)` — in a stable order
    /// suitable for listing to a user.
    pub const ALL: [(&'static str, PassName, &'static str); 7] = [
        (
            "DecomposeToffoli",
            PassName::DecomposeToffoli,
            "Decompose ccx (Toffoli) and ccz gates into Clifford+T",
        ),
        (
            "DecomposeCz",
            PassName::DecomposeCz,
            "Decompose cz gates into H+CX+H",
        ),
        (
            "DecomposeRz",
            PassName::DecomposeRz,
            "Decompose Rz gates into Clifford+T (gridsynth; see --epsilon)",
        ),
        (
            "CancelGates",
            PassName::CancelGates,
            "Cancel adjacent self-inverse gate pairs and reduce Hadamards",
        ),
        (
            "SuperOpt",
            PassName::SuperOpt,
            "Replace small subcircuit windows using a MURM",
        ),
        (
            "PhaseFoldRand",
            PassName::PhaseFoldRand,
            "Merge T/Rz rotations via randomized parity tracking",
        ),
        (
            "CnotMin",
            PassName::CnotMin,
            "Re-synthesize CNOT-dihedral blocks to cut two-qubit gates",
        ),
    ];

    /// Look up a pass by its exact name, as listed in [`PassName::ALL`].
    pub fn parse(s: &str) -> Option<PassName> {
        Self::ALL
            .iter()
            .find(|(n, _, _)| *n == s)
            .map(|(_, p, _)| *p)
    }

    /// Comma-separated list of every valid name (for help / error messages).
    pub fn all_names() -> String {
        Self::ALL
            .iter()
            .map(|(n, _, _)| *n)
            .collect::<Vec<_>>()
            .join(", ")
    }
}

/// Per-run overrides for the SuperOpt window/MURM bounds. `None` means "use
/// whichever preset the optimization level implies" (`DEFAULT_SUPEROPT_*`, or
/// `SUPER_SUPEROPT_*` under [`Level::Osuper`]).
#[derive(Clone, Copy, Debug, Default)]
pub struct SuperOptBounds {
    pub qubits: Option<usize>,
    pub window_gates: Option<usize>,
    pub murm_entries: Option<usize>,
}

/// Which gates SuperOpt may use in MURM representatives.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum SuperOptGates {
    /// Base Clifford+T/CX plus optional native gates present in this stage.
    #[default]
    Auto,
    /// Exactly the built-in Clifford+T/CX basis.
    Base,
    /// Exactly the requested non-empty supported basis.
    Explicit(GateSet),
}

impl SuperOptGates {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "auto" => return Ok(Self::Auto),
            "base" => return Ok(Self::Base),
            "" => {
                return Err(
                    "--superopt-gates requires 'auto', 'base', or a non-empty gate list".to_owned(),
                );
            }
            _ => {}
        }
        let mut basis = GateSet::EMPTY;
        for name in value.split(',') {
            let Some(kind) = GateKind::parse(name) else {
                return Err(format!(
                    "unsupported SuperOpt gate {name:?}; expected a subset of {}",
                    SUPPORTED_GATE_SET
                ));
            };
            if !SUPPORTED_GATE_SET.contains(kind) {
                return Err(format!(
                    "gate {name:?} cannot be emitted by SuperOpt; expected a subset of {}",
                    SUPPORTED_GATE_SET
                ));
            }
            basis.insert(kind);
        }
        if basis.is_empty() {
            return Err("the explicit SuperOpt basis cannot be empty".to_owned());
        }
        Ok(Self::Explicit(basis))
    }

    pub const fn effective(&self, stage_gates: GateSet) -> GateSet {
        match self {
            Self::Auto => BASE_GATE_SET.union(stage_gates.intersection(OPTIONAL_GATE_SET)),
            Self::Base => BASE_GATE_SET,
            Self::Explicit(basis) => *basis,
        }
    }

    pub fn as_argument(&self) -> String {
        match self {
            Self::Auto => "auto".to_owned(),
            Self::Base => "base".to_owned(),
            Self::Explicit(basis) => basis.names().collect::<Vec<_>>().join(","),
        }
    }
}

impl SuperOptBounds {
    /// These bounds with every unset field filled in from `level`'s preset —
    /// the exact `(qubits, window_gates, murm_entries)` a run at `level`
    /// will use. Shared by [`initialize_superopt`] and by callers that want
    /// to *report* the effective configuration (the CLI's `--json` and
    /// `--cache-info`) without duplicating the fallback rules.
    pub fn resolved(self, level: Level) -> (usize, usize, usize) {
        let (qubits, window_gates, murm_entries) = if level == Level::Osuper {
            (
                SUPER_SUPEROPT_QUBITS,
                SUPER_SUPEROPT_WINDOW_GATES,
                SUPER_SUPEROPT_MURM_ENTRIES,
            )
        } else {
            (
                DEFAULT_SUPEROPT_QUBITS,
                DEFAULT_SUPEROPT_WINDOW_GATES,
                DEFAULT_SUPEROPT_MURM_ENTRIES,
            )
        };
        (
            self.qubits.unwrap_or(qubits),
            self.window_gates.unwrap_or(window_gates),
            self.murm_entries.unwrap_or(murm_entries),
        )
    }
}

/// Everything [`optimize`] needs to know about how to optimize a circuit.
///
/// `Default` is the CLI's default: `-O3`, sequential, with all decomposition
/// options disabled.
#[derive(Clone, Debug)]
pub struct Options {
    /// Which default pipeline to run. Ignored when `passes` is set.
    pub level: Level,
    /// An explicit, ordered pass pipeline, replacing the `level` pipeline.
    pub passes: Option<Vec<PassName>>,
    /// Repeat the pipeline until the gate count stops decreasing. Only
    /// consulted for pipelines that aren't already a fixpoint loop (`passes`,
    /// or [`Level::O1`]).
    pub fixpoint: bool,
    /// Decompose Rz gates into Clifford+T via gridsynth.
    pub decompose_rz: bool,
    /// Decompose CZ gates into H+CX+H before optimizing.
    pub decompose_cz: bool,
    /// Decompose CCX and CCZ gates into Clifford+T.
    pub decompose_ccx: bool,
    /// Approximation epsilon for `decompose_rz`.
    pub rz_epsilon: f64,
    /// Optimize gate-contiguous chunks of the circuit in parallel, then
    /// concatenate the results (see [`optimize`]).
    pub parallel: bool,
    /// SuperOpt window/MURM bounds.
    pub superopt: SuperOptBounds,
    /// Gate basis SuperOpt may emit.
    pub superopt_gates: SuperOptGates,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            level: Level::O3,
            passes: None,
            fixpoint: false,
            decompose_rz: false,
            decompose_cz: false,
            decompose_ccx: false,
            rz_epsilon: DEFAULT_RZ_EPSILON,
            parallel: false,
            superopt: SuperOptBounds::default(),
            superopt_gates: SuperOptGates::Auto,
        }
    }
}

/// The metrics tzap reports on a circuit, all counted in one place so callers
/// and progress renderers agree on what "2q gates" or "depth" means.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Metrics {
    pub gates: usize,
    pub two_qubit: usize,
    pub depth: usize,
    pub t: usize,
    pub rz: usize,
}

impl Metrics {
    /// Every counter from one walk of the gate list rather than four.
    ///
    /// This runs after every pass of every fixpoint round, purely to drive the
    /// progress display, so on a multi-million-gate circuit the difference
    /// between one traversal and four is a measurable share of total runtime.
    /// `count_2q`/`count_t`/`count_rz`/`depth` remain as published, separately
    /// callable API; this just doesn't route through them.
    pub fn of(circuit: &Circuit) -> Metrics {
        let mut two_qubit = 0;
        let mut t = 0;
        let mut rz = 0;
        // Depth, computed exactly as `pass::depth` does: a gate lands one layer
        // past the deepest layer already occupied by any of its operands.
        let mut next_layer = vec![0usize; circuit.num_qubits];
        let mut depth = 0;
        for gate in &circuit.gates {
            match gate {
                Gate::cnot { .. } | Gate::cz { .. } => two_qubit += 1,
                Gate::t(_) | Gate::tdg(_) => t += 1,
                Gate::rz(..) => rz += 1,
                _ => {}
            }
            let (arity, operands) = qubit_operands(gate);
            let layer = operands[..arity]
                .iter()
                .map(|&qubit| next_layer[qubit as usize])
                .max()
                .unwrap_or(0)
                + 1;
            for &qubit in &operands[..arity] {
                next_layer[qubit as usize] = layer;
            }
            depth = depth.max(layer);
        }
        Metrics {
            gates: circuit.gates.len(),
            two_qubit,
            depth,
            t,
            rz,
        }
    }

    /// These metrics with `before`'s contribution swapped out for `after`'s —
    /// how a parallel run reports partial progress: the whole circuit's
    /// baseline, adjusted by the chunks finished so far (see
    /// [`run_map_reduce`]).
    ///
    /// Exact for the counting metrics, which are plain sums over the gate
    /// list. Depth is carried through untouched, because it is *not* a sum:
    /// concatenated chunks share layers, so a chunk halving its own depth may
    /// barely move the circuit's. Adjusting it the same way came out 32% low
    /// on `qft_q010_d14171` — a bar that would overstate the reduction all
    /// run and then snap back once the real number arrived. Getting it right
    /// means measuring the whole partially optimized circuit, which is
    /// exactly the O(chunks x circuit) work this exists to avoid.
    ///
    /// Saturating, so a progress bar can never panic a real run.
    fn adjusted(self, before: Metrics, after: Metrics) -> Metrics {
        let adjust =
            |base: usize, before: usize, after: usize| (base + after).saturating_sub(before);
        Metrics {
            gates: adjust(self.gates, before.gates, after.gates),
            two_qubit: adjust(self.two_qubit, before.two_qubit, after.two_qubit),
            depth: self.depth,
            t: adjust(self.t, before.t, after.t),
            rz: adjust(self.rz, before.rz, after.rz),
        }
    }

    /// Accumulate the metrics that are additive across independent chunks.
    fn add_counts(&mut self, other: Metrics) {
        self.gates += other.gates;
        self.two_qubit += other.two_qubit;
        self.t += other.t;
        self.rz += other.rz;
    }
}

/// A consistent snapshot of the completed chunks' contributions to progress.
#[derive(Default)]
struct ChunkProgress {
    done: usize,
    before: Metrics,
    after: Metrics,
}

impl ChunkProgress {
    fn record(&mut self, before: Metrics, after: Metrics, baseline: Metrics) -> (usize, Metrics) {
        self.before.add_counts(before);
        self.after.add_counts(after);
        self.done += 1;
        (self.done, baseline.adjusted(self.before, self.after))
    }
}

/// What [`optimize`] achieved.
#[derive(Clone, Copy, Debug)]
pub struct Report {
    /// The circuit as handed in.
    pub input: Metrics,
    /// Stable comparison metrics for reduction reporting. This is the
    /// original input, including when opt-in decompositions run between the
    /// two optimization stages.
    pub baseline: Metrics,
    /// The optimized circuit.
    pub output: Metrics,
}

/// Anything that can go wrong in [`optimize`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// The MURM backing [`PassName::SuperOpt`] could not be built.
    SuperOpt(SuperOptError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::SuperOpt(error) => write!(f, "failed to initialize SuperOpt: {error}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<SuperOptError> for Error {
    fn from(error: SuperOptError) -> Error {
        Error::SuperOpt(error)
    }
}

/// A sink for a run's progress events, so the driver can report what it's
/// doing without knowing anything about terminals. Every method defaults to
/// doing nothing; implement only the events you care about.
///
/// Events fire from whichever thread reached them: in a `parallel` run,
/// `chunk_done` is called from rayon workers, in completion order. An
/// `Observer` must still be `Sync` because it is shared across workers and
/// other callbacks can run from different threads. The chunk workers
/// themselves are always observed by [`Silent`] — their pipelines run
/// concurrently, so their progress events would interleave into garbage.
pub trait Observer: Sync {
    /// Whether this observer consumes the per-chunk events. When false (the
    /// default), a parallel run skips `chunks_start`/`chunk_done`/`chunks_end`
    /// *and* the per-chunk metric walks that feed them — pure overhead for a
    /// run nobody is watching.
    fn tracks_chunks(&self) -> bool {
        false
    }

    /// A top-level optimization/decomposition stage is starting.
    fn stage_start(&self, _stage: StageKind) {}

    /// A whole-circuit decomposition pass finished.
    /// Both sides are reported: a decomposition *grows* the circuit, so a
    /// renderer that only showed `result` would leave its counts looking
    /// unexplained next to the input's.
    fn pass_done(&self, _name: &str, _input: &Circuit, _result: &Circuit, _elapsed: Duration) {}

    /// The SuperOpt MURM is about to be loaded from disk
    /// (`cached`) or built from scratch.
    fn murm_load_start(&self, _cached: bool, _basis: GateSet) {}

    /// The SuperOpt MURM is ready.
    fn murm_load_done(&self, _cached: bool, _basis: GateSet, _elapsed: Duration) {}

    /// A pass pipeline is starting, against `baseline`. Paired with
    /// [`Observer::progress_end`].
    fn progress_start(&self, _baseline: Metrics) {}

    /// A pass within that pipeline finished. `round` is the (1-based)
    /// fixpoint iteration, or `None` for a pipeline that runs once.
    fn progress_update(&self, _round: Option<usize>, _current: &Circuit, _baseline: Metrics) {}

    /// The pipeline finished. Paired with [`Observer::progress_start`].
    fn progress_end(&self, _baseline: Metrics) {}

    /// The fixpoint driver stopped after `rounds` sweeps, either because a
    /// sweep stopped reducing the gate count (`reached_fixpoint`) or because
    /// it hit its round cap.
    fn fixpoint_done(&self, _rounds: usize, _reached_fixpoint: bool) {}

    /// A parallel run split the circuit into `total` chunks. Only called when
    /// [`Observer::tracks_chunks`] is true.
    fn chunks_start(&self, _total: usize, _baseline: Metrics) {}

    /// `done` of `total` chunks have been optimized; `current` is the whole
    /// circuit's metrics counting finished chunks as optimized and pending
    /// ones as-is.
    ///
    /// The counts are exact. `current.depth` is *not* tracked — it stays at
    /// `baseline.depth` for the whole run, because depth can only be measured
    /// on the assembled circuit (see [`Metrics::adjusted`]). Only called when
    /// [`Observer::tracks_chunks`] is true.
    fn chunk_done(&self, _done: usize, _total: usize, _current: Metrics, _baseline: Metrics) {}

    /// Every chunk is done. Only called when [`Observer::tracks_chunks`] is
    /// true.
    fn chunks_end(&self, _baseline: Metrics) {}
}

/// The observer that reports nothing — what [`optimize`] uses, and what every
/// map-reduce chunk worker gets.
pub struct Silent;

impl Observer for Silent {}

/// Collects the only chunk-local event that has whole-stage meaning. Parallel
/// chunks can take different numbers of fixpoint sweeps, so the stage reports
/// the maximum round count and is converged only when every reporting chunk
/// converged. Other chunk events remain intentionally silent.
struct ParallelFixpointAggregate {
    reports: AtomicUsize,
    max_rounds: AtomicUsize,
    all_converged: AtomicBool,
}

impl ParallelFixpointAggregate {
    fn new() -> Self {
        Self {
            reports: AtomicUsize::new(0),
            max_rounds: AtomicUsize::new(0),
            all_converged: AtomicBool::new(true),
        }
    }

    fn report_to(&self, observer: &dyn Observer) {
        if self.reports.load(Ordering::Relaxed) != 0 {
            observer.fixpoint_done(
                self.max_rounds.load(Ordering::Relaxed),
                self.all_converged.load(Ordering::Relaxed),
            );
        }
    }
}

impl Observer for ParallelFixpointAggregate {
    fn fixpoint_done(&self, rounds: usize, reached_fixpoint: bool) {
        self.max_rounds.fetch_max(rounds, Ordering::Relaxed);
        self.all_converged
            .fetch_and(reached_fixpoint, Ordering::Relaxed);
        self.reports.fetch_add(1, Ordering::Relaxed);
    }
}

/// Number of logical cores, for sizing the rayon thread pool. CPU-bound work
/// like ours gets no benefit from oversubscribing OS threads beyond the core
/// count — it only adds context-switch and cache-thrashing overhead.
fn num_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}

/// Number of map-reduce chunks to split the circuit into. See
/// [`CHUNK_MULTIPLIER`] for why this is more than [`num_threads`].
fn num_par_chunks() -> usize {
    num_threads() * CHUNK_MULTIPLIER
}

/// Build the global rayon pool, sized to the number of logical cores. A
/// no-op if already built (it is process-global).
fn init_global_pool() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads())
        .build_global()
        .ok();
}

/// Run one pass over the whole circuit, timed, reporting it as
/// [`Observer::pass_done`].
fn run_logged(pass: &dyn Pass, circuit: &Circuit, observer: &dyn Observer) -> Circuit {
    let start = Instant::now();
    let c = pass.run(circuit);
    observer.pass_done(pass.name(), circuit, &c, start.elapsed());
    c
}

/// Split a circuit into `num_chunks` gate-contiguous pieces for map-reduce
/// parallelism. Each chunk is later optimized completely independently (see
/// [`run_map_reduce`]) and the results concatenated back in order by
/// [`stitch`] — chunk boundaries are fixed up front and never revisited.
/// `max(1)` guards the empty case where `slice::chunks(0)` would panic.
fn chunk_circuit(circuit: &Circuit, num_chunks: usize) -> Vec<Circuit> {
    let chunk_size = circuit.gates.len().div_ceil(num_chunks).max(1);
    circuit
        .gates
        .chunks(chunk_size)
        .map(|slice| {
            let mut c = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
            for g in slice {
                c.apply(g.clone());
            }
            c
        })
        .collect()
}

/// Concatenate optimized chunks back into a single circuit, in order.
fn stitch(num_qubits: usize, num_cbits: usize, chunks: &[Circuit]) -> Circuit {
    let mut out = Circuit::with_cbits(num_qubits, num_cbits);
    for c in chunks {
        for g in &c.gates {
            out.apply(g.clone());
        }
    }
    out
}

/// Run a pass pipeline once over `circuit`, in order, reporting each pass to
/// `observer` with no round number — unlike the fixpoint driver, this only
/// ever makes one pass over `passes`.
fn run_pipeline(circuit: &Circuit, passes: &[&dyn Pass], observer: &dyn Observer) -> Circuit {
    let baseline = Metrics::of(circuit);
    let mut c = circuit.clone();
    observer.progress_start(baseline);
    observer.progress_update(None, &c, baseline);
    for p in passes {
        c = p.run(&c);
        observer.progress_update(None, &c, baseline);
    }
    observer.progress_end(baseline);
    c
}

/// Run one fixpoint sweep over `circuit`, reporting each pass under
/// `iteration`.
fn run_fixpoint_sweep(
    circuit: &Circuit,
    passes: &[&dyn Pass],
    iteration: usize,
    observer: &dyn Observer,
    baseline: Metrics,
) -> Circuit {
    let mut c = circuit.clone();
    observer.progress_update(Some(iteration), &c, baseline);
    for pass in passes {
        c = pass.run(&c);
        observer.progress_update(Some(iteration), &c, baseline);
    }
    c
}

/// Sweep `passes` until one fails to reduce the gate count, or until
/// `max_rounds` sweeps have run, whichever comes first. Rounds are numbered
/// from `first_round` so a caller running two phases (see [`run_to_fixpoint`])
/// reports one continuous sequence to the observer. Returns the circuit, the
/// last round number used, and whether it converged (as opposed to stopping on
/// the cap).
fn run_fixpoint_phase(
    circuit: &Circuit,
    passes: &[&dyn Pass],
    observer: &dyn Observer,
    baseline: Metrics,
    first_round: usize,
    max_rounds: Option<usize>,
) -> (Circuit, usize, bool) {
    let mut c = circuit.clone();
    let mut round = first_round - 1;
    let mut reduced;
    let mut swept = 0;
    loop {
        round += 1;
        swept += 1;
        let before = c.gates.len();
        c = run_fixpoint_sweep(&c, passes, round, observer, baseline);
        reduced = c.gates.len() < before;
        if !reduced || max_rounds.is_some_and(|m| swept >= m) {
            break;
        }
    }
    (c, round, !reduced)
}

/// Repeatedly run `passes` until a sweep fails to reduce the gate count, or
/// (when `max_rounds` is given) until that many sweeps have run, whichever
/// comes first.
///
/// When `rz_decompose` is given it runs exactly once, after the pre-synthesis
/// sweeps have converged — identically at every SuperOpt level (O2, O3,
/// Osuper), which all reach this through `optimize_default`.
///
/// That placement matters a lot: what gridsynth expands into is what
/// `SuperOpt`'s greedy, non-backtracking window selection has to work with
/// afterwards, and converging first can shrink that circuit by several-fold (on
/// cobble's ols-ridge, 180k gates after one sweep vs 25k at the fixpoint).
/// Synthesizing into the smaller circuit measured 2-22% fewer T across cobble,
/// and *faster* despite running more rounds, since every post-synthesis round
/// then sweeps far fewer gates.
///
/// Afterwards, if there were Rz gates to decompose, the sweeps resume on the
/// synthesized circuit. `max_rounds` is a budget shared across both phases,
/// except that the post-synthesis phase always gets at least one sweep —
/// freshly synthesized Clifford+T sequences are never left unoptimized just
/// because the pre-synthesis phase used up the cap. So O2, the one capped
/// level, runs up to `max_rounds + 1` sweeps when synthesis intervenes; paying
/// that extra sweep is what buys O2 the same Rz placement as the uncapped
/// levels rather than silently degrading to synthesize-after-one-sweep.
fn run_to_fixpoint(
    circuit: &Circuit,
    passes: &[&dyn Pass],
    rz_decompose: Option<&dyn Pass>,
    observer: &dyn Observer,
    max_rounds: Option<usize>,
) -> Circuit {
    let baseline = Metrics::of(circuit);
    observer.progress_start(baseline);

    let (mut c, mut round, mut converged) =
        run_fixpoint_phase(circuit, passes, observer, baseline, 1, max_rounds);

    if let Some(pass) = rz_decompose {
        let had_rz = c.gates.iter().any(|g| matches!(g, Gate::rz(..)));
        c = pass.run(&c);
        observer.progress_update(Some(round), &c, baseline);
        if had_rz {
            let spent = round;
            let post_cap = max_rounds.map(|m| m.saturating_sub(spent).max(1));
            (c, round, converged) =
                run_fixpoint_phase(&c, passes, observer, baseline, round + 1, post_cap);
        }
    }

    observer.progress_end(baseline);
    observer.fixpoint_done(round, converged);
    c
}

/// Build a fresh `SuperOpt` instance. Callers must construct one of these
/// per map-reduce chunk (never share or reuse one instance across chunks):
/// each instance owns its own matrix cache and incremental-diff state, so a
/// fresh instance per chunk means `.incremental()` is always sound here —
/// every instance only ever sees successive versions of the one circuit it
/// was built for. `level` selects which bounds preset an unset
/// [`SuperOptBounds`] field falls back to — `SUPER_SUPEROPT_*` under
/// [`Level::Osuper`], `DEFAULT_SUPEROPT_*` otherwise.
fn initialize_superopt(
    options: &Options,
    level: Level,
    basis: GateSet,
    observer: &dyn Observer,
) -> Result<SuperOpt, Error> {
    let (qubits, window_gates, murm_entries) = options.superopt.resolved(level);

    // A MURM entry needs strictly fewer gates than the window it replaces
    // (see `ActiveWindow::consider`'s `local.len() >= gate_indices.len()`
    // rejection), and no window ever exceeds `window_gates`. So a stored
    // representative at exactly `window_gates` depth could never be strictly
    // smaller than the largest possible window — `window_gates - 1` is the
    // deepest depth any MURM entry can ever be used at.
    let murm_gates = window_gates.saturating_sub(1);
    let murm_config = MurmConfig::new(qubits, murm_gates, murm_entries).with_basis(basis);
    // Captured before the build/load below can create the cache file (which
    // would make a second `murm_is_cached` call always say "cached").
    let cached = murm_is_cached(murm_config);
    observer.murm_load_start(cached, basis);
    let start = Instant::now();
    let pass = SuperOpt::new(qubits, window_gates, murm_config)?;
    observer.murm_load_done(cached, basis, start.elapsed());
    Ok(pass.without_subcircuits().incremental())
}

/// Run `optimize` once on the whole circuit when sequential, or independently
/// on each of `num_chunks` chunks in parallel (map), recombining the results
/// in order (reduce). Each chunk is optimized completely independently: no
/// state — not even a MURM's matrix cache — is shared between
/// chunks, so `optimize` must construct every stateful pass (`SuperOpt`)
/// fresh on every call.
fn run_map_reduce(
    circuit: &Circuit,
    parallel: bool,
    num_chunks: usize,
    observer: &dyn Observer,
    optimize: impl Fn(&Circuit, &dyn Observer) -> Result<Circuit, Error> + Sync + Send,
) -> Result<Circuit, Error> {
    if !parallel {
        return optimize(circuit, observer);
    }
    let chunks = chunk_circuit(circuit, num_chunks);
    let total = chunks.len();
    let tracking = observer.tracks_chunks();

    // Whole-circuit baseline, which each finished chunk then adjusts by its
    // own before/after difference: pending chunks (not yet optimized)
    // contribute their original metrics, completed ones their current metrics.
    //
    // Measuring the partially optimized circuit directly instead — stitching
    // every chunk back together under a lock on each completion — cost
    // O(chunks x circuit) and serialized the workers behind it: ~20% of a
    // parallel run on a 4M-gate circuit, all of it to move a progress bar.
    let baseline = Metrics::of(circuit);
    let progress = Mutex::new(ChunkProgress::default());
    let fixpoint = ParallelFixpointAggregate::new();

    if tracking {
        observer.chunks_start(total, baseline);
    }
    let optimized: Vec<Result<Circuit, Error>> = chunks
        .par_iter()
        .map(|chunk| {
            // Walked only when someone is watching: two passes over the chunk
            // that a `Silent` run has no use for.
            let before = tracking.then(|| Metrics::of(chunk));
            let result = optimize(chunk, &fixpoint)?;
            if let Some(before) = before {
                // Keep the totals, completion number, and callback in one
                // critical section. Otherwise workers can report mixed
                // before/after snapshots or deliver callbacks out of order.
                let after = Metrics::of(&result);
                let mut progress = progress.lock().unwrap();
                let (n, current) = progress.record(before, after, baseline);
                observer.chunk_done(n, total, current, baseline);
            }
            Ok(result)
        })
        .collect();
    if tracking {
        observer.chunks_end(baseline);
    }
    let optimized = optimized.into_iter().collect::<Result<Vec<_>, _>>()?;
    fixpoint.report_to(observer);
    Ok(stitch(circuit.num_qubits, circuit.num_cbits, &optimized))
}

/// Run one named pass over the current whole-circuit stage. SuperOpt resolves
/// `auto` lazily from that exact circuit, after every earlier named pass. In a
/// parallel run the whole-circuit basis is fixed before chunking, then each
/// worker constructs its own stateful pass instance.
fn run_explicit_pass(
    circuit: &Circuit,
    name: PassName,
    options: &Options,
    num_chunks: usize,
    observer: &dyn Observer,
) -> Result<Circuit, Error> {
    macro_rules! map_pass {
        ($make:expr) => {
            run_map_reduce(
                circuit,
                options.parallel,
                num_chunks,
                observer,
                |chunk, _| {
                    let pass = $make;
                    Ok(Pass::run(&pass, chunk))
                },
            )
        };
    }

    match name {
        PassName::DecomposeToffoli => map_pass!(DecomposeToffoli),
        PassName::DecomposeCz => map_pass!(DecomposeCz),
        PassName::DecomposeRz => map_pass!(DecomposeRz {
            epsilon: options.rz_epsilon,
        }),
        PassName::CancelGates => map_pass!(CancelGates),
        PassName::PhaseFoldRand => map_pass!(PhaseFoldRand),
        PassName::CnotMin => map_pass!(CnotMin::default()),
        PassName::SuperOpt => {
            let basis = options.superopt_gates.effective(circuit.gate_set());
            if basis.is_empty() {
                return Ok(circuit.clone());
            }
            if options.parallel {
                initialize_superopt(options, Level::O1, basis, observer)?;
            }
            run_map_reduce(
                circuit,
                options.parallel,
                num_chunks,
                observer,
                |chunk, chunk_observer| {
                    let pass = initialize_superopt(options, Level::O1, basis, chunk_observer)?;
                    Ok(Pass::run(&pass, chunk))
                },
            )
        }
    }
}

/// Run one sweep of an explicit pipeline, resolving stage-sensitive state at
/// each pass boundary.
fn run_explicit_sweep(
    circuit: &Circuit,
    names: &[PassName],
    options: &Options,
    num_chunks: usize,
    observer: &dyn Observer,
    round: Option<usize>,
    baseline: Metrics,
) -> Result<Circuit, Error> {
    let mut current = circuit.clone();
    if !options.parallel {
        observer.progress_update(round, &current, baseline);
    }
    for &name in names {
        current = run_explicit_pass(&current, name, options, num_chunks, observer)?;
        if !options.parallel {
            observer.progress_update(round, &current, baseline);
        }
    }
    Ok(current)
}

/// Run the explicit [`Options::passes`] pipeline. Unlike the default preset,
/// every pass boundary is semantically visible: a SuperOpt after a
/// decomposition resolves `auto` from the decomposed circuit, and a later
/// SuperOpt may therefore use a different MURM again.
fn optimize_explicit(
    circuit: &Circuit,
    options: &Options,
    num_chunks: usize,
    observer: &dyn Observer,
) -> Result<Circuit, Error> {
    let names = options
        .passes
        .as_ref()
        .expect("only called when `passes` is set");
    let baseline = Metrics::of(circuit);
    if !options.parallel {
        observer.progress_start(baseline);
    }

    let mut current = circuit.clone();
    if options.fixpoint {
        let mut round = 0;
        let converged = loop {
            round += 1;
            let before = current.gates.len();
            current = run_explicit_sweep(
                &current,
                names,
                options,
                num_chunks,
                observer,
                Some(round),
                baseline,
            )?;
            if current.gates.len() >= before {
                break true;
            }
        };
        observer.fixpoint_done(round, converged);
    } else {
        current = run_explicit_sweep(
            &current, names, options, num_chunks, observer, None, baseline,
        )?;
    }

    if !options.parallel {
        observer.progress_end(baseline);
    }
    Ok(current)
}

/// Run the default pipeline for [`Options::level`] on `circuit`, constructing
/// a fresh `SuperOpt` whenever the level uses one. Used as one map-reduce
/// worker per chunk.
fn optimize_default(
    circuit: &Circuit,
    options: &Options,
    superopt_basis: GateSet,
    observer: &dyn Observer,
) -> Result<Circuit, Error> {
    let cancel_pass = CancelGates;
    let global = PhaseFoldRand;
    let cnot_min_pass = CnotMin::default();

    if level_uses_superopt(options.level) {
        // O2 runs a fixed 2 rounds rather than to a true fixpoint — O3 and
        // Osuper are the "run it out fully" tiers; O2 is the cheap, bounded
        // one.
        let max_rounds = (options.level == Level::O2).then_some(2);
        let superopt_pass = (!superopt_basis.is_empty())
            .then(|| initialize_superopt(options, options.level, superopt_basis, observer))
            .transpose()?;
        // CnotMin leads the sweep: it re-synthesizes whole CNOT-dihedral
        // blocks, which reshapes the circuit far more than the peephole
        // rewriter does, and the passes after it then work on the result.
        let mut passes: Vec<&dyn Pass> = vec![&cnot_min_pass, &cancel_pass];
        if let Some(superopt_pass) = &superopt_pass {
            passes.push(superopt_pass);
        }
        passes.push(&global);
        Ok(run_to_fixpoint(
            circuit, &passes, None, observer, max_rounds,
        ))
    } else {
        let optimization_passes: Vec<&dyn Pass> = vec![&cancel_pass, &global];
        Ok(if options.fixpoint {
            run_to_fixpoint(circuit, &optimization_passes, None, observer, None)
        } else {
            run_pipeline(circuit, &optimization_passes, observer)
        })
    }
}

#[derive(Clone, Copy)]
enum OptimizationStage {
    Explicit,
    Preset { forbidden: GateSet },
}

/// Execute one structured optimization stage. This is the single boundary at
/// which whole-circuit gate inspection, MURM policy, parallel prewarming, and
/// observer stage identity meet.
fn run_optimization_stage(
    circuit: &Circuit,
    kind: StageKind,
    stage: OptimizationStage,
    options: &Options,
    num_chunks: usize,
    observer: &dyn Observer,
) -> Result<Circuit, Error> {
    observer.stage_start(kind);
    match stage {
        OptimizationStage::Explicit => optimize_explicit(circuit, options, num_chunks, observer),
        OptimizationStage::Preset { forbidden } => {
            let basis = options
                .superopt_gates
                .effective(circuit.gate_set())
                .difference(forbidden);
            if options.parallel && level_uses_superopt(options.level) && !basis.is_empty() {
                initialize_superopt(options, options.level, basis, observer)?;
            }
            run_map_reduce(
                circuit,
                options.parallel,
                num_chunks,
                observer,
                |chunk, chunk_observer| optimize_default(chunk, options, basis, chunk_observer),
            )
        }
    }
}

/// Whether `level`'s pipeline includes a SuperOpt pass (and so pays for a
/// MURM).
fn level_uses_superopt(level: Level) -> bool {
    matches!(level, Level::O2 | Level::O3 | Level::Osuper)
}

/// Assert output gate contracts independently of unitary equivalence. The
/// latter cannot detect a decomposed gate being synthesized back into the
/// output by a later optimizer stage.
fn check_output_invariants(input: &Circuit, result: &Circuit, decomposed: GateSet) {
    let input_gates = input.gate_set();
    let output_gates = result.gate_set();
    if output_gates.contains(GateKind::Rz) && !input_gates.contains(GateKind::Rz) {
        panic!("BUG: output contains Rz gates but input did not");
    }
    for kind in decomposed.iter() {
        assert!(
            !output_gates.contains(kind),
            "BUG: output contains {} gates after their requested decomposition",
            kind.name()
        );
    }
}

/// Optimize `circuit`, reporting nothing along the way.
///
/// ```rust,ignore
/// use tzap::circuit::Circuit;
/// use tzap::optimize::{Level, Options, optimize};
///
/// let circuit = Circuit::from_qasm(qasm)?;
/// let (optimized, report) = optimize(&circuit, &Options::default())?;
/// println!("{} → {} T", report.baseline.t, report.output.t);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn optimize(circuit: &Circuit, options: &Options) -> Result<(Circuit, Report), Error> {
    optimize_with(circuit, options, &Silent)
}

/// Optimize `circuit`, reporting progress to `observer`.
///
/// The default workflow optimizes the input circuit, applies requested
/// CCX/CCZ, CZ, and Rz decompositions in that order, then optimizes again if
/// any decomposition changed the circuit. [`Options::passes`] overrides that
/// with an explicit, user-ordered pipeline. Under [`Options::parallel`], each
/// optimization stage runs map-reduce style (see [`run_map_reduce`]).
pub fn optimize_with(
    circuit: &Circuit,
    options: &Options,
    observer: &dyn Observer,
) -> Result<(Circuit, Report), Error> {
    let input = Metrics::of(circuit);
    let num_chunks = num_par_chunks();

    // Explicit pipeline via `passes`: run exactly what the caller listed, in
    // order, with the original input as the reporting baseline.
    if let Some(names) = &options.passes {
        let uses_rz = names.iter().any(|p| matches!(p, PassName::DecomposeRz));
        if options.parallel || uses_rz {
            init_global_pool();
        }
        let result = run_optimization_stage(
            circuit,
            StageKind::ExplicitPipeline,
            OptimizationStage::Explicit,
            options,
            num_chunks,
            observer,
        )?;
        check_output_invariants(circuit, &result, GateSet::EMPTY);
        let report = Report {
            input,
            baseline: input,
            output: Metrics::of(&result),
        };
        return Ok((result, report));
    }

    if options.parallel || options.decompose_rz {
        init_global_pool();
    }

    // Native optimization always comes first. Requested decompositions form
    // an opt-in middle stage, followed by the same optimizer when any pass
    // actually changed the circuit.
    let native = run_optimization_stage(
        circuit,
        StageKind::InputOptimization,
        OptimizationStage::Preset {
            forbidden: GateSet::EMPTY,
        },
        options,
        num_chunks,
        observer,
    )?;
    let mut result = native;
    let mut forbidden = GateSet::EMPTY;
    if options.decompose_ccx {
        forbidden = forbidden.union(GateSet::from_bits_const(
            (1 << GateKind::Ccx as u8) | (1 << GateKind::Ccz as u8),
        ));
    }
    if options.decompose_cz {
        forbidden = forbidden.union(GateSet::singleton(GateKind::Cz));
    }
    if options.decompose_rz {
        forbidden = forbidden.union(GateSet::singleton(GateKind::Rz));
    }
    let mut decomposed = false;

    let decompositions: [(bool, GateSet, StageKind, &dyn Pass); 3] = [
        (
            options.decompose_ccx,
            GateSet::from_bits_const((1 << GateKind::Ccx as u8) | (1 << GateKind::Ccz as u8)),
            StageKind::DecomposeCcx,
            &DecomposeToffoli,
        ),
        (
            options.decompose_cz,
            GateSet::singleton(GateKind::Cz),
            StageKind::DecomposeCz,
            &DecomposeCz,
        ),
        (
            options.decompose_rz,
            GateSet::singleton(GateKind::Rz),
            StageKind::DecomposeRz,
            &DecomposeRz {
                epsilon: options.rz_epsilon,
            },
        ),
    ];
    for (requested, kinds, stage, pass) in decompositions {
        if requested && !result.gate_set().intersection(kinds).is_empty() {
            observer.stage_start(stage);
            result = run_logged(pass, &result, observer);
            decomposed = true;
        }
    }
    if decomposed {
        result = run_optimization_stage(
            &result,
            StageKind::PostDecompositionOptimization,
            OptimizationStage::Preset { forbidden },
            options,
            num_chunks,
            observer,
        )?;
    }

    check_output_invariants(circuit, &result, forbidden);
    let report = Report {
        input,
        baseline: input,
        output: Metrics::of(&result),
    };
    Ok((result, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qasm;

    /// Parallel (map-reduce) optimization of a measured circuit must
    /// round-trip to valid QASM. Regression guard: the stitched
    /// reconstruction must use the circuit's `num_cbits` (not `Circuit::new`,
    /// which zeroes it and drops the `creg` declaration).
    #[test]
    fn parallel_mode_preserves_creg() {
        let mut c = Circuit::with_cbits(1, 1);
        c.apply(Gate::h(0));
        c.apply(Gate::measure { qubit: 0, cbit: 0 });

        let cancel = CancelGates;
        let par = run_map_reduce(&c, true, 4, &Silent, |chunk, obs| {
            Ok(run_pipeline(chunk, &[&cancel], obs))
        })
        .expect("no SuperOpt pass, so nothing can fail");

        assert_eq!(
            par.num_cbits, 1,
            "parallel reconstruction must preserve num_cbits"
        );

        let out = qasm::serialize(&par);
        assert!(
            out.contains("creg"),
            "parallel output has measurements but no creg declaration:\n{out}"
        );
        assert!(
            qasm::parse(&out).is_ok(),
            "parallel output must round-trip to valid QASM:\n{out}"
        );
    }

    /// Optimizing an empty circuit in parallel must not panic. Regression guard:
    /// chunk_size must never be 0 (else `slice::chunks(0)` panics).
    #[test]
    fn parallel_mode_handles_empty_circuit() {
        let empty = Circuit::new(1);
        let out = run_map_reduce(&empty, true, 4, &Silent, |chunk, obs| {
            Ok(run_pipeline(chunk, &[], obs))
        })
        .expect("no SuperOpt pass, so nothing can fail");
        assert!(out.gates.is_empty());
    }

    #[test]
    fn parallel_fixpoint_reports_the_maximum_chunk_round() {
        #[derive(Default)]
        struct FixpointObserver(std::sync::Mutex<Vec<(usize, bool)>>);
        impl Observer for FixpointObserver {
            fn fixpoint_done(&self, rounds: usize, converged: bool) {
                self.0.lock().unwrap().push((rounds, converged));
            }
        }

        let mut circuit = Circuit::new(4);
        for qubit in 0..4 {
            circuit.apply(Gate::h(qubit));
        }
        let observer = FixpointObserver::default();
        run_map_reduce(&circuit, true, 4, &observer, |chunk, chunk_observer| {
            let qubit = match chunk.gates[0] {
                Gate::h(qubit) => qubit,
                _ => unreachable!(),
            };
            chunk_observer.fixpoint_done(qubit as usize + 1, qubit != 2);
            Ok(chunk.clone())
        })
        .unwrap();

        assert_eq!(
            *observer.0.lock().unwrap(),
            vec![(4, false)],
            "parallel stages report one record: max rounds and all-chunks convergence"
        );
    }

    /// The parallel progress numbers are derived from per-chunk deltas rather
    /// than measured off the partially stitched circuit, so the counting
    /// metrics must still land exactly on the real output's once every chunk
    /// has reported.
    #[test]
    fn chunk_progress_counts_match_the_finished_circuit() {
        #[derive(Default)]
        struct LastChunk(std::sync::Mutex<Option<Metrics>>);
        impl Observer for LastChunk {
            fn tracks_chunks(&self) -> bool {
                true
            }
            fn chunk_done(&self, _done: usize, _total: usize, current: Metrics, _: Metrics) {
                *self.0.lock().unwrap() = Some(current);
            }
        }

        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\n".to_string()
            + &"h q[0];\ncx q[0],q[1];\nt q[1];\ncx q[1],q[2];\ntdg q[2];\nrz(0.3) q[0];\n"
                .repeat(20);
        let c = Circuit::from_qasm(&qasm).unwrap();
        let cancel = CancelGates;
        let observer = LastChunk::default();
        let out = run_map_reduce(&c, true, 4, &observer, |chunk, obs| {
            Ok(run_pipeline(chunk, &[&cancel], obs))
        })
        .unwrap();

        let reported = observer.0.lock().unwrap().expect("every chunk reports");
        let actual = Metrics::of(&out);
        assert_eq!(reported.gates, actual.gates);
        assert_eq!(reported.two_qubit, actual.two_qubit);
        assert_eq!(reported.t, actual.t);
        assert_eq!(reported.rz, actual.rz);
    }

    /// A slow first callback must not let later workers publish their progress
    /// ahead of it, even when their optimization has already finished.
    #[test]
    fn chunk_progress_callbacks_are_ordered() {
        #[derive(Default)]
        struct OrderedChunks {
            reported: Mutex<Vec<usize>>,
        }
        impl Observer for OrderedChunks {
            fn tracks_chunks(&self) -> bool {
                true
            }
            fn chunk_done(&self, done: usize, _total: usize, _: Metrics, _: Metrics) {
                if done == 1 {
                    std::thread::sleep(Duration::from_millis(50));
                }
                self.reported.lock().unwrap().push(done);
            }
        }

        let mut circuit = Circuit::new(16);
        for qubit in 0..16 {
            circuit.apply(Gate::h(qubit));
        }
        let observer = OrderedChunks::default();
        rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| run_map_reduce(&circuit, true, 16, &observer, |chunk, _| Ok(chunk.clone())))
            .unwrap();

        assert_eq!(
            *observer.reported.lock().unwrap(),
            (1..=16).collect::<Vec<_>>()
        );
    }

    /// Depth is not a sum over chunks, so a chunk reporting its own depth
    /// collapsing must leave the circuit's reported depth alone rather than
    /// subtracting from it.
    #[test]
    fn adjusted_leaves_depth_at_the_baseline() {
        let baseline = Metrics {
            gates: 10,
            two_qubit: 4,
            depth: 6,
            t: 2,
            rz: 0,
        };
        // A chunk that locally held more depth than the whole circuit did,
        // and optimized all of it away.
        let before = Metrics {
            depth: 9,
            ..baseline
        };
        let adjusted = baseline.adjusted(before, Metrics::default());
        assert_eq!(adjusted.depth, baseline.depth);
        assert_eq!(adjusted.gates, 0, "the counts still adjust exactly");
    }

    /// Each map-reduce chunk must get its own independent `SuperOpt`
    /// instance — the whole point of this design (no shared `MatrixStore`,
    /// no `Arc`, no cross-chunk locking). Regression guard: running the same
    /// SuperOpt-using pipeline on a circuit split into several chunks must
    /// not panic or deadlock, and must produce a valid, gate-count-bounded
    /// output (each chunk's `SuperOpt` instance only ever rewrites within
    /// its own chunk).
    #[test]
    fn parallel_mode_gives_each_chunk_its_own_superopt() {
        let mut c = Circuit::new(2);
        for _ in 0..8 {
            c.apply(Gate::h(0));
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
        }

        let options = Options {
            level: Level::O2,
            parallel: true,
            ..Options::default()
        };

        let out = run_map_reduce(&c, true, 4, &Silent, |chunk, obs| {
            optimize_default(chunk, &options, BASE_GATE_SET, obs)
        })
        .expect("the default SuperOpt MURM must build");
        assert!(out.gates.len() <= c.gates.len());
    }

    /// Rz synthesis is an opt-in middle stage and its output is optimized
    /// again at every SuperOpt level.
    #[test]
    fn rz_synthesis_waits_for_the_presynthesis_fixpoint_at_every_level() {
        // An H·H / CNOT·CNOT body (so sweep 1 finds real reduction and the
        // pre-synthesis phase runs past round 1) plus a single non-π/4 Rz for
        // gridsynth to synthesize.
        let mut c = Circuit::new(3);
        for _ in 0..12 {
            c.apply(Gate::h(0));
            c.apply(Gate::h(0));
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
            c.apply(Gate::cnot {
                control: 0,
                target: 1,
            });
        }
        c.apply(Gate::rz(0.37, 2));

        for level in [Level::O2, Level::O3, Level::Osuper] {
            let options = Options {
                level,
                decompose_rz: true,
                rz_epsilon: 1e-3,
                // The placement logic is bound-independent, so use the tiny
                // CI MURM while still exercising each level's driver path.
                superopt: tiny_superopt_bounds(),
                ..Options::default()
            };
            let (out, _) = optimize(&c, &options).expect("pipeline must run");

            assert!(
                !out.gates.iter().any(|g| matches!(g, Gate::rz(..))),
                "{level:?}: decompose_rz must leave no Rz behind"
            );
        }
    }

    /// Native gates remain native unless their decomposition flag is set.
    #[test]
    fn ccx_decomposition_is_opt_in() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });

        let options = Options {
            level: Level::O1,
            ..Options::default()
        };
        let (native, report) = optimize(&c, &options).expect("O1 builds no MURM");

        assert_eq!(report.input.gates, 1);
        assert_eq!(report.baseline, report.input);
        assert!(native.gate_set().contains(GateKind::Ccx));

        let options = Options {
            level: Level::O1,
            decompose_ccx: true,
            ..Options::default()
        };
        let (decomposed, report) = optimize(&c, &options).expect("O1 builds no MURM");
        assert!(!decomposed.gate_set().contains(GateKind::Ccx));
        assert_eq!(report.baseline, report.input);
    }

    /// A circuit needing no decomposition reports `input == baseline`.
    #[test]
    fn report_baseline_equals_input_without_decomposition() {
        let mut c = Circuit::new(1);
        c.apply(Gate::h(0));
        c.apply(Gate::h(0));

        let options = Options {
            level: Level::O1,
            ..Options::default()
        };
        let (out, report) = optimize(&c, &options).expect("O1 builds no SuperOpt MURM");

        assert_eq!(report.input, report.baseline);
        assert_eq!(out.gates.len(), 0, "HH must cancel");
        assert_eq!(report.output.gates, 0);
    }

    #[test]
    fn every_explicit_superopt_basis_mask_round_trips() {
        let candidates = [
            GateKind::H,
            GateKind::X,
            GateKind::Z,
            GateKind::S,
            GateKind::Sdg,
            GateKind::T,
            GateKind::Tdg,
            GateKind::Cx,
            GateKind::Cz,
            GateKind::Ccx,
            GateKind::Ccz,
        ];
        assert!(SuperOptGates::parse("").is_err());
        for mask in 1u16..1 << candidates.len() {
            let expected = GateSet::from_kinds(
                candidates
                    .into_iter()
                    .enumerate()
                    .filter_map(|(bit, kind)| (mask & (1 << bit) != 0).then_some(kind)),
            );
            let argument = expected.names().collect::<Vec<_>>().join(",");
            assert_eq!(
                SuperOptGates::parse(&argument).unwrap(),
                SuperOptGates::Explicit(expected),
                "mask {mask:#05x}"
            );
        }
        assert_eq!(
            SuperOptGates::parse("h,h,cz").unwrap(),
            SuperOptGates::Explicit(GateSet::from_kinds([GateKind::H, GateKind::Cz]))
        );
        assert!(SuperOptGates::parse("rz").is_err());
    }

    #[test]
    fn auto_superopt_basis_covers_all_native_gate_subsets() {
        for mask in 0u8..8 {
            let mut stage = GateSet::from_kinds([GateKind::H, GateKind::Rz]);
            for (bit, kind) in OPTIONAL_GATE_KINDS.into_iter().enumerate() {
                if mask & (1 << bit) != 0 {
                    stage.insert(kind);
                }
            }
            let expected = BASE_GATE_SET.union(stage.intersection(OPTIONAL_GATE_SET));
            assert_eq!(SuperOptGates::Auto.effective(stage), expected);
            assert_eq!(SuperOptGates::Base.effective(stage), BASE_GATE_SET);
        }
    }

    #[test]
    fn all_native_and_decomposition_subsets_follow_the_staged_contract() {
        use crate::unitary::circuits_equiv;

        for native_mask in 0u8..8 {
            let mut input = Circuit::new(3);
            input.apply(Gate::h(0));
            input.apply(Gate::rz(0.37, 2));
            if native_mask & 1 != 0 {
                input.apply(Gate::cz {
                    control: 0,
                    target: 1,
                });
            }
            if native_mask & 2 != 0 {
                input.apply(Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                });
            }
            if native_mask & 4 != 0 {
                input.apply(Gate::ccz {
                    control1: 0,
                    control2: 1,
                    target: 2,
                });
            }

            for decompose_mask in 0u8..8 {
                let options = Options {
                    level: Level::O1,
                    decompose_ccx: decompose_mask & 1 != 0,
                    decompose_cz: decompose_mask & 2 != 0,
                    decompose_rz: decompose_mask & 4 != 0,
                    rz_epsilon: 1e-3,
                    ..Options::default()
                };
                let (output, _) = optimize(&input, &options).unwrap();
                if options.decompose_ccx {
                    assert!(!output.gate_set().contains(GateKind::Ccx));
                    assert!(!output.gate_set().contains(GateKind::Ccz));
                }
                if options.decompose_cz {
                    assert!(!output.gate_set().contains(GateKind::Cz));
                }
                if options.decompose_rz {
                    assert!(!output.gate_set().contains(GateKind::Rz));
                }
                assert!(Circuit::from_qasm(&output.to_qasm()).is_ok());
                assert!(
                    circuits_equiv(&input, &output, 5e-3),
                    "native mask {native_mask:#05b}, decomposition mask {decompose_mask:#05b}"
                );
            }
        }
    }

    fn native_subset_circuit(mask: u8) -> Circuit {
        let mut circuit = Circuit::new(3);
        circuit.apply(Gate::h(0));
        circuit.apply(Gate::t(1));
        circuit.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        circuit.apply(Gate::rz(0.37, 2));
        if mask & 1 != 0 {
            circuit.apply(Gate::cz {
                control: 0,
                target: 1,
            });
        }
        if mask & 2 != 0 {
            circuit.apply(Gate::ccx {
                control1: 0,
                control2: 1,
                target: 2,
            });
        }
        if mask & 4 != 0 {
            circuit.apply(Gate::ccz {
                control1: 0,
                control2: 1,
                target: 2,
            });
        }
        circuit
    }

    fn tiny_superopt_bounds() -> SuperOptBounds {
        SuperOptBounds {
            qubits: Some(3),
            window_gates: Some(4),
            murm_entries: Some(500),
        }
    }

    #[derive(Default)]
    struct BasisObserver(std::sync::Mutex<Vec<GateSet>>);

    impl Observer for BasisObserver {
        fn murm_load_done(&self, _cached: bool, basis: GateSet, _elapsed: Duration) {
            self.0.lock().unwrap().push(basis);
        }
    }

    #[derive(Default)]
    struct StageBasisState {
        current: Option<StageKind>,
        loads: Vec<(StageKind, GateSet)>,
    }

    #[derive(Default)]
    struct StageBasisObserver(std::sync::Mutex<StageBasisState>);

    impl Observer for StageBasisObserver {
        fn stage_start(&self, stage: StageKind) {
            self.0.lock().unwrap().current = Some(stage);
        }

        fn murm_load_done(&self, _cached: bool, basis: GateSet, _elapsed: Duration) {
            let mut state = self.0.lock().unwrap();
            let stage = state.current.expect("MURM load belongs to a stage");
            state.loads.push((stage, basis));
        }
    }

    #[test]
    fn auto_basis_is_recomputed_after_decomposition() {
        let mut input = Circuit::new(3);
        input.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        let observer = StageBasisObserver::default();
        let options = Options {
            level: Level::O2,
            decompose_ccx: true,
            superopt: tiny_superopt_bounds(),
            ..Options::default()
        };

        let (output, _) = optimize_with(&input, &options, &observer).unwrap();
        assert!(!output.gate_set().contains(GateKind::Ccx));
        assert_eq!(
            observer.0.lock().unwrap().loads,
            vec![
                (
                    StageKind::InputOptimization,
                    BASE_GATE_SET.union(GateSet::singleton(GateKind::Ccx)),
                ),
                (StageKind::PostDecompositionOptimization, BASE_GATE_SET),
            ]
        );
    }

    #[test]
    fn requested_decomposition_constrains_an_explicit_post_basis() {
        for kind in [GateKind::Cz, GateKind::Ccx, GateKind::Ccz] {
            let mut input = Circuit::new(3);
            match kind {
                GateKind::Cz => input.apply(Gate::cz {
                    control: 0,
                    target: 1,
                }),
                GateKind::Ccx => input.apply(Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                }),
                GateKind::Ccz => input.apply(Gate::ccz {
                    control1: 0,
                    control2: 1,
                    target: 2,
                }),
                _ => unreachable!(),
            }
            let observer = StageBasisObserver::default();
            let options = Options {
                level: Level::O2,
                decompose_ccx: matches!(kind, GateKind::Ccx | GateKind::Ccz),
                decompose_cz: kind == GateKind::Cz,
                superopt: tiny_superopt_bounds(),
                superopt_gates: SuperOptGates::Explicit(GateSet::singleton(kind)),
                ..Options::default()
            };

            let (output, _) = optimize_with(&input, &options, &observer).unwrap();
            assert!(!output.gate_set().contains(kind), "{kind:?}");
            assert_eq!(
                observer.0.lock().unwrap().loads,
                vec![(StageKind::InputOptimization, GateSet::singleton(kind))],
                "{kind:?}"
            );
        }
    }

    #[test]
    fn absent_requested_gate_is_still_forbidden_from_the_post_basis() {
        let mut input = Circuit::new(3);
        input.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        let observer = StageBasisObserver::default();
        let explicit = BASE_GATE_SET.union(GateSet::from_kinds([
            GateKind::Cz,
            GateKind::Ccx,
            GateKind::Ccz,
        ]));
        let options = Options {
            level: Level::O2,
            decompose_ccx: true,
            decompose_cz: true,
            superopt: tiny_superopt_bounds(),
            superopt_gates: SuperOptGates::Explicit(explicit),
            ..Options::default()
        };

        let (output, _) = optimize_with(&input, &options, &observer).unwrap();
        assert!(output.gate_set().intersection(OPTIONAL_GATE_SET).is_empty());
        assert_eq!(
            observer.0.lock().unwrap().loads,
            vec![
                (StageKind::InputOptimization, explicit),
                (StageKind::PostDecompositionOptimization, BASE_GATE_SET),
            ]
        );
    }

    #[test]
    fn explicit_pipeline_recomputes_auto_at_every_superopt_boundary() {
        let mut ccx = Circuit::new(3);
        ccx.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        let mut cz = Circuit::new(2);
        cz.apply(Gate::cz {
            control: 0,
            target: 1,
        });

        let cases = [
            (
                &ccx,
                vec![PassName::DecomposeToffoli, PassName::SuperOpt],
                vec![BASE_GATE_SET],
            ),
            (
                &cz,
                vec![PassName::DecomposeCz, PassName::SuperOpt],
                vec![BASE_GATE_SET],
            ),
            (
                &ccx,
                vec![
                    PassName::SuperOpt,
                    PassName::DecomposeToffoli,
                    PassName::SuperOpt,
                ],
                vec![
                    BASE_GATE_SET.union(GateSet::singleton(GateKind::Ccx)),
                    BASE_GATE_SET,
                ],
            ),
        ];

        for parallel in [false, true] {
            for (input, passes, expected) in &cases {
                let observer = BasisObserver::default();
                let options = Options {
                    passes: Some(passes.clone()),
                    parallel,
                    superopt: tiny_superopt_bounds(),
                    ..Options::default()
                };
                optimize_with(input, &options, &observer).unwrap();
                assert_eq!(
                    *observer.0.lock().unwrap(),
                    *expected,
                    "{passes:?}, {parallel}"
                );
            }
        }
    }

    /// Eight input profiles × (`auto`, `base`, and eight explicit profiles).
    #[test]
    fn superopt_selection_matrix_covers_80_end_to_end_cases() {
        use crate::unitary::circuits_equiv;

        for input_mask in 0u8..8 {
            let input = native_subset_circuit(input_mask);
            let mut modes = vec![SuperOptGates::Auto, SuperOptGates::Base];
            for explicit_mask in 0u8..8 {
                let optional =
                    GateSet::from_kinds(OPTIONAL_GATE_KINDS.into_iter().enumerate().filter_map(
                        |(bit, kind)| (explicit_mask & (1 << bit) != 0).then_some(kind),
                    ));
                modes.push(SuperOptGates::Explicit(BASE_GATE_SET.union(optional)));
            }

            for (mode_index, mode) in modes.into_iter().enumerate() {
                let expected_basis = mode.effective(input.gate_set());
                let observer = BasisObserver::default();
                let options = Options {
                    level: Level::O2,
                    superopt: tiny_superopt_bounds(),
                    superopt_gates: mode,
                    ..Options::default()
                };
                let (output, _) = optimize_with(&input, &options, &observer).unwrap();
                assert_eq!(
                    *observer.0.lock().unwrap(),
                    vec![expected_basis],
                    "input mask {input_mask:#05b}, mode {mode_index}"
                );

                // SuperOpt may preserve optional gates already in the input,
                // or emit ones explicitly enabled by the basis, but cannot
                // introduce any other optional native kind.
                let allowed_optional = input
                    .gate_set()
                    .union(expected_basis)
                    .intersection(OPTIONAL_GATE_SET);
                assert!(
                    output
                        .gate_set()
                        .intersection(OPTIONAL_GATE_SET)
                        .is_subset(allowed_optional),
                    "input mask {input_mask:#05b}, mode {mode_index}: {}",
                    output.gate_set()
                );
                assert!(Circuit::from_qasm(&output.to_qasm()).is_ok());
                assert!(
                    circuits_equiv(&input, &output, 1e-9),
                    "input mask {input_mask:#05b}, mode {mode_index}"
                );
            }
        }
    }

    /// Four levels × all eight optional-native input profiles.
    #[test]
    fn optimization_level_matrix_covers_32_end_to_end_cases() {
        use crate::unitary::circuits_equiv;

        for level in [Level::O1, Level::O2, Level::O3, Level::Osuper] {
            for native_mask in 0u8..8 {
                let input = native_subset_circuit(native_mask);
                let options = Options {
                    level,
                    superopt: tiny_superopt_bounds(),
                    ..Options::default()
                };
                let (output, _) = optimize(&input, &options).unwrap();
                assert!(Circuit::from_qasm(&output.to_qasm()).is_ok());
                assert!(
                    circuits_equiv(&input, &output, 1e-9),
                    "level {level:?}, native mask {native_mask:#05b}"
                );
            }
        }
    }

    /// Eight input profiles × auto/base/representative-explicit in parallel.
    #[test]
    fn parallel_superopt_matrix_covers_24_end_to_end_cases() {
        use crate::unitary::circuits_equiv;

        let explicit = SuperOptGates::Explicit(
            BASE_GATE_SET.union(GateSet::from_kinds([GateKind::Cz, GateKind::Ccz])),
        );
        for native_mask in 0u8..8 {
            let input = native_subset_circuit(native_mask);
            for mode in [SuperOptGates::Auto, SuperOptGates::Base, explicit.clone()] {
                let sequential = Options {
                    level: Level::O2,
                    superopt: tiny_superopt_bounds(),
                    superopt_gates: mode.clone(),
                    ..Options::default()
                };
                let parallel = Options {
                    parallel: true,
                    ..sequential.clone()
                };
                let (seq_output, _) = optimize(&input, &sequential).unwrap();
                let (par_output, _) = optimize(&input, &parallel).unwrap();
                assert!(circuits_equiv(&input, &seq_output, 1e-9));
                assert!(circuits_equiv(&input, &par_output, 1e-9));
                assert!(
                    circuits_equiv(&seq_output, &par_output, 1e-9),
                    "native mask {native_mask:#05b}, mode {mode:?}"
                );
            }
        }
    }

    /// Bounded deterministic fuzzing: every optional-native input subset and
    /// every decomposition subset gets eight reproducible random circuits.
    /// Basis and execution modes cycle across the samples, and failure
    /// messages include the seed for replay.
    #[test]
    fn bounded_native_pipeline_fuzz_covers_512_cross_feature_cases() {
        use crate::pass::Pass;
        use crate::unitary::circuits_equiv;

        fn next(seed: &mut u64) -> u64 {
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            *seed
        }

        let candidates = [
            GateKind::H,
            GateKind::X,
            GateKind::Z,
            GateKind::S,
            GateKind::Sdg,
            GateKind::T,
            GateKind::Tdg,
            GateKind::Cx,
            GateKind::Cz,
            GateKind::Ccx,
            GateKind::Ccz,
        ];

        for native_mask in 0u8..8 {
            for decompose_mask in 0u8..8 {
                for sample in 0u64..8 {
                    let replay_seed = 0x4e41_5449_5645_0000u64
                        | (sample << 12)
                        | ((native_mask as u64) << 8)
                        | decompose_mask as u64;
                    let mut seed = replay_seed;
                    let mut input = native_subset_circuit(native_mask);
                    for _ in 0..24 {
                        let q = (next(&mut seed) % 3) as u32;
                        match next(&mut seed) % 12 {
                            0 => input.apply(Gate::h(q)),
                            1 => input.apply(Gate::x(q)),
                            2 => input.apply(Gate::z(q)),
                            3 => input.apply(Gate::s(q)),
                            4 => input.apply(Gate::sdg(q)),
                            5 => input.apply(Gate::t(q)),
                            6 => input.apply(Gate::tdg(q)),
                            7 => input.apply(Gate::rz((next(&mut seed) % 13 + 1) as f64 / 17.0, q)),
                            8 => input.apply(Gate::cnot {
                                control: q,
                                target: (q + 1) % 3,
                            }),
                            9 if native_mask & 1 != 0 => input.apply(Gate::cz {
                                control: q,
                                target: (q + 1) % 3,
                            }),
                            10 if native_mask & 2 != 0 => input.apply(Gate::ccx {
                                control1: q,
                                control2: (q + 1) % 3,
                                target: (q + 2) % 3,
                            }),
                            11 if native_mask & 4 != 0 => input.apply(Gate::ccz {
                                control1: q,
                                control2: (q + 1) % 3,
                                target: (q + 2) % 3,
                            }),
                            _ => input.apply(Gate::h(q)),
                        }
                    }

                    let mode = match sample % 3 {
                        0 => SuperOptGates::Auto,
                        1 => SuperOptGates::Base,
                        _ => {
                            let mask = (next(&mut seed) as u16) & 0x07ff;
                            let nonempty = if mask == 0 { 1 } else { mask };
                            SuperOptGates::Explicit(GateSet::from_kinds(
                                candidates
                                    .into_iter()
                                    .enumerate()
                                    .filter_map(|(bit, kind)| {
                                        (nonempty & (1 << bit) != 0).then_some(kind)
                                    }),
                            ))
                        }
                    };
                    let options = Options {
                        level: Level::O2,
                        decompose_ccx: decompose_mask & 1 != 0,
                        decompose_cz: decompose_mask & 2 != 0,
                        decompose_rz: decompose_mask & 4 != 0,
                        rz_epsilon: 1e-3,
                        parallel: sample & 1 != 0,
                        superopt: SuperOptBounds {
                            qubits: Some(3),
                            window_gates: Some(4),
                            murm_entries: Some(128),
                        },
                        superopt_gates: mode,
                        ..Options::default()
                    };
                    let (output, _) = optimize(&input, &options).unwrap();
                    assert!(
                        circuits_equiv(&input, &output, 2e-2),
                        "seed {replay_seed:#018x}, native {native_mask:#05b}, decomposition {decompose_mask:#05b}"
                    );
                    assert!(Circuit::from_qasm(&output.to_qasm()).is_ok());
                    if options.decompose_ccx {
                        assert!(!output.gate_set().contains(GateKind::Ccx));
                        assert!(!output.gate_set().contains(GateKind::Ccz));
                    }
                    if options.decompose_cz {
                        assert!(!output.gate_set().contains(GateKind::Cz));
                    }
                    if options.decompose_rz {
                        assert!(!output.gate_set().contains(GateKind::Rz));
                    }

                    let once = CancelGates.run(&output);
                    let twice = CancelGates.run(&once);
                    assert_eq!(
                        once.gates, twice.gates,
                        "CancelGates not idempotent for seed {replay_seed:#018x}"
                    );
                }
            }
        }
    }

    /// The literal Cartesian product is intentionally manual/nightly: the
    /// factored tests above cover the same axes cheaply in normal CI.
    #[test]
    #[ignore = "extended 5,120-case native-gate matrix"]
    fn extended_full_native_configuration_matrix_has_5120_cases() {
        use crate::unitary::circuits_equiv;

        let mut cases = 0usize;
        for input_mask in 0u8..8 {
            let input = native_subset_circuit(input_mask);
            let mut modes = vec![SuperOptGates::Auto, SuperOptGates::Base];
            for explicit_mask in 0u8..8 {
                let optional =
                    GateSet::from_kinds(OPTIONAL_GATE_KINDS.into_iter().enumerate().filter_map(
                        |(bit, kind)| (explicit_mask & (1 << bit) != 0).then_some(kind),
                    ));
                modes.push(SuperOptGates::Explicit(BASE_GATE_SET.union(optional)));
            }
            for mode in modes {
                for decompose_mask in 0u8..8 {
                    for level in [Level::O1, Level::O2, Level::O3, Level::Osuper] {
                        for parallel in [false, true] {
                            let options = Options {
                                level,
                                decompose_ccx: decompose_mask & 1 != 0,
                                decompose_cz: decompose_mask & 2 != 0,
                                decompose_rz: decompose_mask & 4 != 0,
                                rz_epsilon: 1e-3,
                                parallel,
                                superopt: tiny_superopt_bounds(),
                                superopt_gates: mode.clone(),
                                ..Options::default()
                            };
                            let (output, _) = optimize(&input, &options).unwrap();
                            assert!(circuits_equiv(&input, &output, 5e-3));
                            if options.decompose_ccx {
                                assert!(!output.gate_set().contains(GateKind::Ccx));
                                assert!(!output.gate_set().contains(GateKind::Ccz));
                            }
                            if options.decompose_cz {
                                assert!(!output.gate_set().contains(GateKind::Cz));
                            }
                            if options.decompose_rz {
                                assert!(!output.gate_set().contains(GateKind::Rz));
                            }
                            cases += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(cases, 5_120);
    }
}
