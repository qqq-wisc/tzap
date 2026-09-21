use std::env;
use std::fs;
use std::io::{self, Read};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use tzap::circuit::{Circuit, GateSet};
use tzap::optimize::{Metrics, Observer, Report, StageKind, optimize_with};
#[cfg(test)]
use tzap::super_opt::BASE_GATE_SET;

mod cli;
mod json;
mod progress;
mod ui;

use cli::{AUTO_PARALLEL_GATES, Action, Opts, Run, STREAM_PATH, arg_error, parse_args};
use json::{FixpointRecord, MurmRecord, PassRecord, Recording, RunInfo};
use progress::{box_lines, fmt_num, fmt_size};
use ui::Ui;

/// Renders a run's progress to the terminal: the per-pass result lines, the
/// SuperOpt MURM-load status, and the live redrawn progress boxes. The whole
/// of the CLI's output during optimization, and the only thing standing
/// between `tzap::optimize` and a silent run.
///
/// Every write goes through [`Ui`], which decides whether this run may color
/// and whether it may redraw in place at all — so the same observer serves a
/// terminal, a pipe, and `--quiet`. It doubles as the recorder for `--json`:
/// the events it renders are exactly the ones that report needs, which is why
/// a run that renders nothing still reports everything.
struct Terminal {
    ui: Ui,
    /// Filled in only when `--json` will consume it — a run nobody asked a
    /// report from shouldn't pay for the bookkeeping, and `pass_done` fires
    /// under the same lock every parallel chunk contends for.
    recording: Option<Mutex<Recording>>,
}

impl Terminal {
    fn new(ui: Ui, json: bool) -> Terminal {
        Terminal {
            ui,
            recording: json.then(|| Mutex::new(Recording::default())),
        }
    }

    /// Add to the `--json` recording, if one is being kept.
    fn record(&self, f: impl FnOnce(&mut Recording)) {
        if let Some(recording) = &self.recording {
            f(&mut recording.lock().expect("--json recording mutex poisoned"));
        }
    }

    /// Whether this run will draw progress boxes at all, which needs a live
    /// terminal to redraw on.
    ///
    /// When false, the driver skips the progress events *and* the per-pass
    /// and per-chunk `Metrics` walks that feed them — pure waste for a piped
    /// or quiet run that would render none of it. Nothing in the `--json`
    /// report comes from these events, so skipping them costs it nothing.
    fn draws_progress(&self) -> bool {
        self.ui.live()
    }

    /// Take the recording this run accumulated, leaving an empty one behind
    /// (and returning an empty one when `--json` wasn't asked for). Takes
    /// `&self` rather than consuming: the observer also owns the [`Ui`] the
    /// report is written through.
    fn take_recording(&self) -> Recording {
        self.recording
            .as_ref()
            .map_or_else(Recording::default, |recording| {
                std::mem::take(&mut *recording.lock().expect("--json recording mutex poisoned"))
            })
    }
}

/// Number of bar rows in a reduction progress box — the Rz row only appears
/// for circuits that have Rz gates to report on.
fn reduction_rows(baseline: Metrics) -> usize {
    if baseline.rz > 0 { 5 } else { 4 }
}

/// Number of bar rows in the parallel chunk progress box: [`reduction_rows`]
/// with a Chunks row added on top and the Depth row dropped (a parallel run
/// can't track depth cheaply — see `Metrics::adjusted`), which happens to
/// leave the two boxes the same height.
fn chunk_rows(baseline: Metrics) -> usize {
    reduction_rows(baseline)
}

impl Observer for Terminal {
    fn tracks_chunks(&self) -> bool {
        self.draws_progress()
    }

    fn stage_start(&self, stage: StageKind) {
        self.record(|recording| recording.stages.push(stage.name().to_string()));
        self.ui.info(&format!("  {}", stage.name()));
    }

    /// Report one pass with timing and a result line, followed by a blank
    /// separator line — this and the SuperOpt MURM-load message each own a
    /// trailing blank line, so a live progress box that follows never needs to
    /// print one itself. [`read_circuit`] deliberately does *not* trail with a
    /// blank: it should stay flush with whatever comes right after it, whether
    /// that's this, the MURM message, or a box directly.
    fn pass_done(&self, name: &str, input: &Circuit, result: &Circuit, elapsed: Duration) {
        let metrics = Metrics::of(result);
        self.record(|recording| {
            let stage = recording.current_stage();
            recording.passes.push(PassRecord {
                stage,
                name: name.to_string(),
                input_gates: input.gates.len(),
                output_gates: metrics.gates,
                seconds: elapsed.as_secs_f64(),
            })
        });
        let rz_report = (metrics.rz > 0).then(|| format!(" · {} Rz", fmt_num(metrics.rz)));
        // Gate count shows both sides: a decomposition grows the circuit, and
        // the final result banner measures its reduction against *this* count,
        // not the parsed one. Printing only the post-decomposition figure left
        // readers to guess where it came from.
        self.ui.info(&format!(
            "  {}\n\t└─ {} → {} gates · {} 2q gates · {} T/Tdg{} · {} depth · {:.3}s",
            name,
            fmt_num(input.gates.len()),
            fmt_num(metrics.gates),
            fmt_num(metrics.two_qubit),
            fmt_num(metrics.t),
            rz_report.as_deref().unwrap_or(""),
            fmt_num(metrics.depth),
            elapsed.as_secs_f64()
        ));
        self.ui.blank();
    }

    /// One name for the synthesis artifact in every message.
    fn murm_load_start(&self, cached: bool, _basis: GateSet) {
        if cached {
            // Reading a large cached MURM off disk can itself take a
            // moment, so say so before it starts; overwritten in place with
            // the completed message below rather than left as its own line.
            self.ui.start_inline("  Loading MURM...");
        } else {
            self.ui
                .info("  🔧 Building MURM (one-time — cached for future use)...");
        }
    }

    fn murm_load_done(&self, cached: bool, basis: GateSet, elapsed: Duration) {
        self.record(|recording| {
            let stage = recording.current_stage();
            recording.murms.push(MurmRecord {
                stage,
                cached,
                basis,
                seconds: elapsed.as_secs_f64(),
            })
        });
        let message = format!("  Loaded MURM in {:.3}s", elapsed.as_secs_f64());
        if cached {
            self.ui.finish_inline(&message);
        } else {
            self.ui.info(&message);
        }
        self.ui.info(&format!("\t└─ Synthesis basis: {basis}"));
        self.ui.blank();
    }

    fn progress_start(&self, baseline: Metrics) {
        self.ui
            .start_progress_block(box_lines(reduction_rows(baseline)));
    }

    fn progress_update(&self, round: Option<usize>, current: &Circuit, baseline: Metrics) {
        if !self.draws_progress() {
            return;
        }
        let m = Metrics::of(current);
        match round {
            Some(round) => self.ui.update_fixpoint_progress(
                round,
                m.gates,
                m.two_qubit,
                m.depth,
                m.t,
                baseline.gates,
                baseline.two_qubit,
                baseline.depth,
                baseline.t,
                m.rz,
                baseline.rz,
            ),
            None => self.ui.update_reduction_progress(
                "% reduction so far",
                m.gates,
                m.two_qubit,
                m.depth,
                m.t,
                baseline.gates,
                baseline.two_qubit,
                baseline.depth,
                baseline.t,
                m.rz,
                baseline.rz,
            ),
        }
    }

    fn progress_end(&self, baseline: Metrics) {
        self.ui
            .end_progress_block(box_lines(reduction_rows(baseline)));
    }

    fn fixpoint_done(&self, rounds: usize, reached_fixpoint: bool) {
        self.record(|recording| {
            let stage = recording.current_stage();
            recording.fixpoints.push(FixpointRecord {
                stage,
                rounds,
                converged: reached_fixpoint,
            })
        });
        if reached_fixpoint {
            let plural = if rounds == 1 { "round" } else { "rounds" };
            self.ui
                .info(&format!("  Converged after {rounds} {plural}"));
            self.ui.blank();
        }
    }

    fn chunks_start(&self, total: usize, baseline: Metrics) {
        self.ui
            .start_progress_block(box_lines(chunk_rows(baseline)));
        self.chunk_done(0, total, baseline, baseline);
    }

    fn chunk_done(&self, done: usize, total: usize, current: Metrics, baseline: Metrics) {
        self.ui.update_chunk_progress(
            done,
            total,
            baseline.gates,
            current.gates,
            baseline.two_qubit,
            current.two_qubit,
            baseline.t,
            current.t,
            baseline.rz,
            current.rz,
        );
    }

    fn chunks_end(&self, baseline: Metrics) {
        self.ui.end_progress_block(box_lines(chunk_rows(baseline)));
    }
}

/// Write a circuit to the output destination (if any), exiting on error.
/// `-` writes it to stdout, so tzap can sit in the middle of a pipeline;
/// every message tzap prints goes to stderr, which is what keeps that
/// stdout stream clean enough to pipe into a parser.
fn write_output(ui: &Ui, run: &Run, circuit: &Circuit) {
    let Some(path) = &run.output_path else {
        return;
    };
    let output = circuit.to_qasm();
    if run.writes_stdout() {
        ui.write_stdout(&output);
        return;
    }
    fs::write(path, &output).unwrap_or_else(|e| ui.abort(&format!("Error writing {path}: {e}")));
    ui.info(&format!("  wrote {path}"));
}

/// A parsed input circuit and what the CLI knows about where it came from.
struct Parsed {
    circuit: Circuit,
    /// Size of the input, in bytes; `None` when it couldn't be determined.
    bytes: Option<u64>,
    seconds: f64,
}

/// Read and parse the input into a circuit, logging parse stats. Exits on
/// error. Overwrites its own "Parsing..." line with "Parsed..." in place
/// (see [`Ui::start_inline`]/[`Ui::finish_inline`]) once done; without a live
/// terminal only the completed line is printed.
fn read_circuit(ui: &Ui, path: &str) -> Parsed {
    let stdin = path == STREAM_PATH;
    let label = if stdin { "<stdin>" } else { path };
    // A file's size is known before the read, so it can be shown while the
    // read is still in flight; stdin's isn't known until it has all arrived.
    let known_size = (!stdin)
        .then(|| fs::metadata(path).map(|m| m.len()).ok())
        .flatten();
    let parse_start = Instant::now();
    ui.start_inline(&match known_size {
        Some(size) => format!("  Parsing {label} ({})", fmt_size(size)),
        None => format!("  Parsing {label}"),
    });

    let qasm = if stdin {
        let mut buffer = String::new();
        io::stdin()
            .read_to_string(&mut buffer)
            .unwrap_or_else(|e| ui.abort(&format!("Error reading stdin: {e}")));
        buffer
    } else {
        fs::read_to_string(path).unwrap_or_else(|e| ui.abort(&format!("Error reading {path}: {e}")))
    };
    let bytes = known_size.or(Some(qasm.len() as u64));
    let circuit = Circuit::from_qasm(&qasm)
        .unwrap_or_else(|e| ui.abort(&format!("Error parsing {label}: {e}")));
    let seconds = parse_start.elapsed().as_secs_f64();
    ui.finish_inline(&match bytes {
        Some(size) => format!("  Parsed {label} ({}) in {seconds:.3}s", fmt_size(size)),
        None => format!("  Parsed {label} in {seconds:.3}s"),
    });
    ui.info(&format!(
        "\t├─ {} qubits · {} gates",
        fmt_num(circuit.num_qubits),
        fmt_num(circuit.gates.len()),
    ));
    ui.info(&format!("\t└─ Circuit gates: {}", circuit.gate_set()));
    Parsed {
        circuit,
        bytes,
        seconds,
    }
}

/// `--cache-info`: where tzap's on-disk MURMs live and what they
/// cost. A query, so it answers on stdout — including under `--quiet`, which
/// silences commentary, not the thing that was asked for.
fn print_cache_info(ui: &Ui, json: bool) {
    let entries = tzap::super_opt::cache_entries();
    if json {
        ui.write_stdout(&json::render_cache_info(&entries));
        return;
    }
    let Some(dir) = tzap::super_opt::cache_dir() else {
        ui.write_stdout(
            "No cache directory: none of --cache-dir, $TZAP_CACHE_DIR, \
             $XDG_CACHE_HOME, $HOME, %LOCALAPPDATA%, or %USERPROFILE% is set, \
             so MURMs are rebuilt every run.\n",
        );
        return;
    };
    let total: u64 = entries.iter().map(|entry| entry.bytes).sum();
    let plural = if entries.len() == 1 { "MURM" } else { "MURMs" };
    let mut out = format!(
        "Cache directory: {}\n{} cached {plural} · {}\n",
        dir.display(),
        entries.len(),
        fmt_size(total)
    );
    for entry in &entries {
        let name = entry
            .path
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| entry.path.display().to_string());
        out.push_str(&format!("  {name}  {}\n", fmt_size(entry.bytes)));
    }
    ui.write_stdout(&out);
}

/// `--clear-cache`: delete every cached MURM. The summary is
/// commentary on an action rather than a queried result, so it goes to
/// stderr and `--quiet` silences it; `--json` puts the machine-readable list
/// on stdout as usual.
fn clear_cache(ui: &Ui, json: bool) {
    let removed = tzap::super_opt::clear_cache()
        .unwrap_or_else(|e| ui.abort(&format!("Error clearing the MURM cache: {e}")));
    if json {
        ui.write_stdout(&json::render_cache_info(&removed));
        return;
    }
    let total: u64 = removed.iter().map(|entry| entry.bytes).sum();
    let plural = if removed.len() == 1 { "MURM" } else { "MURMs" };
    ui.info(&format!(
        "  Removed {} cached {plural} · {} freed",
        removed.len(),
        fmt_size(total)
    ));
}

/// Print the result banner against the original input baseline and write the
/// output file (if requested).
fn finish(ui: &Ui, report: &Report, result: &Circuit, run: &Run, start: Instant) {
    ui.print_result(
        report.baseline.gates,
        report.output.gates,
        report.baseline.two_qubit,
        report.output.two_qubit,
        report.baseline.depth,
        report.output.depth,
        report.baseline.t,
        report.output.t,
        report.baseline.rz,
        report.output.rz,
        start.elapsed().as_secs_f64(),
    );

    write_output(ui, run, result);
}

fn main() {
    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    let opts = parse_args(&args);
    let Opts { action, ui, json } = opts;

    let mut run = match action {
        Action::CacheInfo => return print_cache_info(&ui, json),
        Action::ClearCache => return clear_cache(&ui, json),
        Action::Optimize(run) => run,
    };

    ui.info(&format!(
        "{}⚡\u{FE0F} tzap{} {}v{}{}",
        ui.sgr("\x1b[1m"),
        ui.reset(),
        ui.sgr("\x1b[2m"),
        env!("CARGO_PKG_VERSION"),
        ui.reset()
    ));
    let parsed = read_circuit(&ui, &run.input_path);
    // Said out loud only when the size decided it: a run nobody asked to
    // parallelize otherwise reaches the chunked progress box with no
    // explanation, and a slightly different gate count at the end.
    if run.resolve_parallel(parsed.circuit.gates.len()) {
        ui.info(&format!(
            "  Optimizing chunks in parallel ({}+ gates) · --no-parallel to disable",
            fmt_num(AUTO_PARALLEL_GATES)
        ));
    }
    let observer = Terminal::new(ui, json);
    let (result, report) =
        optimize_with(&parsed.circuit, &run.options, &observer).unwrap_or_else(|e| arg_error(e));
    finish(&observer.ui, &report, &result, &run, start);

    if json {
        let info = RunInfo {
            input_path: (!run.reads_stdin()).then_some(run.input_path.as_str()),
            input_bytes: parsed.bytes,
            input_qubits: parsed.circuit.num_qubits,
            input_gate_set: parsed.circuit.gate_set(),
            parse_seconds: parsed.seconds,
            output_path: run.output_path.as_deref(),
            output_gate_set: result.gate_set(),
            seconds: start.elapsed().as_secs_f64(),
        };
        let recording = observer.take_recording();
        observer
            .ui
            .write_stdout(&json::render(&info, &run.options, &report, &recording));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cli::Action;
    use tzap::optimize::{Level, SuperOptBounds};

    /// The `Run` a successful `parse_args` produced, or a panic naming what
    /// it produced instead — every test here is about optimization runs.
    fn parse_run(args: &[&str]) -> Run {
        let args: Vec<String> = std::iter::once("tzap".to_string())
            .chain(args.iter().map(|s| s.to_string()))
            .collect();
        match parse_args(&args).action {
            Action::Optimize(run) => run,
            _ => panic!("expected an optimization run for {args:?}"),
        }
    }

    /// The Rz row appears in a progress box exactly when the circuit has Rz
    /// gates to report on — the box's height must match, or the live redraw
    /// leaves a stray line behind.
    #[test]
    fn rz_row_only_counted_when_present() {
        let with_rz = Metrics {
            rz: 3,
            ..Metrics::default()
        };
        assert_eq!(reduction_rows(Metrics::default()), 4);
        assert_eq!(reduction_rows(with_rz), 5);
        assert_eq!(chunk_rows(Metrics::default()), 4);
        assert_eq!(chunk_rows(with_rz), 5);
    }

    /// An absent `-O` flag must behave exactly like `-O3`, and the hidden
    /// `--superopt-*` flags must reach the optimizer.
    #[test]
    fn parse_args_defaults_to_o3() {
        let run = parse_run(&["in.qasm"]);
        assert_eq!(run.options.level, Level::O3);
        assert!(run.options.passes.is_none());
        assert!(!run.options.parallel);

        let run = parse_run(&["in.qasm", "--superopt-qubits", "4"]);
        let SuperOptBounds { qubits, .. } = run.options.superopt;
        assert_eq!(qubits, Some(4));
    }

    /// Parallelism follows the circuit's size unless it was asked for
    /// outright: the flags win at any size, and with neither of them the
    /// threshold decides.
    #[test]
    fn parallelism_follows_the_circuit_size_by_default() {
        let mut run = parse_run(&["in.qasm"]);
        assert!(!run.resolve_parallel(AUTO_PARALLEL_GATES - 1));
        assert!(!run.options.parallel, "below the threshold: sequential");

        let mut run = parse_run(&["in.qasm"]);
        assert!(
            run.resolve_parallel(AUTO_PARALLEL_GATES),
            "the threshold is inclusive, and the caller is told why"
        );
        assert!(run.options.parallel);

        // Asked for outright: honored at any size, and never announced as a
        // decision the size made.
        let mut run = parse_run(&["in.qasm", "--parallel"]);
        assert!(!run.resolve_parallel(1));
        assert!(run.options.parallel);

        let mut run = parse_run(&["in.qasm", "--no-parallel"]);
        assert!(!run.resolve_parallel(AUTO_PARALLEL_GATES * 10));
        assert!(!run.options.parallel);
    }

    /// `--parallel` and `--no-parallel` are a pair, so the last one on the
    /// line wins rather than one of them being privileged.
    #[test]
    fn the_last_parallel_flag_wins() {
        let mut run = parse_run(&["in.qasm", "--parallel", "--no-parallel"]);
        run.resolve_parallel(AUTO_PARALLEL_GATES);
        assert!(!run.options.parallel);

        let mut run = parse_run(&["in.qasm", "--no-parallel", "--parallel"]);
        run.resolve_parallel(1);
        assert!(run.options.parallel);
    }

    /// `-` names the standard streams on either side, and is never mistaken
    /// for a flag.
    #[test]
    fn dash_selects_the_standard_streams() {
        let run = parse_run(&["-", "-o", "-"]);
        assert!(run.reads_stdin());
        assert!(run.writes_stdout());

        let run = parse_run(&["in.qasm", "-"]);
        assert!(!run.reads_stdin());
        assert!(run.writes_stdout());

        let run = parse_run(&["in.qasm"]);
        assert!(!run.reads_stdin());
        assert!(!run.writes_stdout());
    }

    /// The live rendering path — the progress boxes and their cursor motion —
    /// can only be reached with a terminal on stderr, which no test has: an
    /// integration test's streams are both pipes, and a piped run correctly
    /// renders none of this. So drive the whole `Observer` surface here
    /// against a `Ui` that claims a terminal, which at minimum pins that
    /// every event renders, that the box heights the block bracketing
    /// reserves match what gets drawn into them, and that no arithmetic in
    /// the bar/box layout can panic on a real run's numbers.
    #[test]
    fn the_live_observer_renders_every_event() {
        let circuit = Circuit::from_qasm(
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\n\
             h q[0];\ncx q[0],q[1];\nt q[1];\nrz(0.3) q[2];\ntdg q[1];\n",
        )
        .expect("fixture parses");
        let baseline = Metrics::of(&circuit);
        let observer = Terminal::new(Ui::live_for_tests(), true);
        assert!(observer.draws_progress());
        assert!(observer.tracks_chunks());

        let elapsed = Duration::from_millis(7);
        observer.pass_done("Toffoli decomposition", &circuit, &circuit, elapsed);
        for cached in [true, false] {
            observer.murm_load_start(cached, BASE_GATE_SET);
            observer.murm_load_done(cached, BASE_GATE_SET, elapsed);
        }

        // A sequential pipeline, then a fixpoint one, then a parallel run —
        // each bracketed the way the driver brackets it.
        observer.progress_start(baseline);
        observer.progress_update(None, &circuit, baseline);
        for round in 1..=3 {
            observer.progress_update(Some(round), &circuit, baseline);
        }
        observer.progress_end(baseline);
        observer.fixpoint_done(3, true);
        observer.fixpoint_done(2, false);

        observer.chunks_start(4, baseline);
        for done in 1..=4 {
            observer.chunk_done(done, 4, baseline, baseline);
        }
        observer.chunks_end(baseline);

        observer.ui.print_result(
            baseline.gates,
            baseline.gates / 2,
            baseline.two_qubit,
            baseline.two_qubit,
            baseline.depth,
            baseline.depth - 1,
            baseline.t,
            0,
            baseline.rz,
            baseline.rz,
            0.125,
        );

        // The pass and fixpoint events still feed `--json` from the live path
        // exactly as they do from the silent one.
        let recording = observer.take_recording();
        assert_eq!(recording.passes.len(), 1);
        assert_eq!(
            recording
                .fixpoints
                .iter()
                .map(|f| f.rounds)
                .collect::<Vec<_>>(),
            vec![3, 2]
        );
        assert_eq!(recording.murms.len(), 2);
    }

    /// A circuit with no gates at all still renders: every bar is a division
    /// by a zero baseline, and every box has to survive it.
    #[test]
    fn the_live_observer_renders_an_empty_circuit() {
        let circuit = Circuit::from_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\n")
            .expect("fixture parses");
        let baseline = Metrics::of(&circuit);
        assert_eq!(baseline.gates, 0);

        let observer = Terminal::new(Ui::live_for_tests(), false);
        observer.progress_start(baseline);
        observer.progress_update(Some(1), &circuit, baseline);
        observer.progress_end(baseline);
        observer.chunks_start(0, baseline);
        observer.chunks_end(baseline);
        observer.ui.print_result(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0);
    }

    /// Without a live terminal the driver is told not to bother: the chunk
    /// events, and the per-chunk metric walks behind them, are skipped
    /// outright rather than computed and discarded.
    #[test]
    fn a_piped_observer_asks_for_no_progress_work() {
        let observer = Terminal::new(Ui::plain(), false);
        assert!(!observer.draws_progress());
        assert!(!observer.tracks_chunks());
    }

    /// A recording is kept only for the runs that will consume one.
    #[test]
    fn json_recording_is_kept_only_when_asked_for() {
        let observer = Terminal::new(Ui::plain(), false);
        observer.record(|recording| {
            recording.fixpoints.push(FixpointRecord {
                stage: None,
                rounds: 1,
                converged: true,
            })
        });
        assert!(observer.take_recording().fixpoints.is_empty());

        let observer = Terminal::new(Ui::plain(), true);
        observer.record(|recording| {
            recording.fixpoints.push(FixpointRecord {
                stage: None,
                rounds: 2,
                converged: false,
            })
        });
        let recording = observer.take_recording();
        assert_eq!(recording.fixpoints.first().map(|f| f.rounds), Some(2));
        assert!(
            observer.take_recording().fixpoints.is_empty(),
            "the recording is taken, not copied"
        );
    }
}
