use std::env;
use std::fs;
use std::io::{self, Read};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use tzap::circuit::{Circuit, Gate, GateSet};
use tzap::optimize::{Metrics, Observer, Report, StageKind, optimize_with};
use tzap::pbc::{AxisKind, AxisWeights};
#[cfg(test)]
use tzap::super_opt::BASE_GATE_SET;

mod cli;
mod json;
mod progress;
mod ui;

use cli::{AUTO_PARALLEL_GATES, Action, Opts, Run, STREAM_PATH, arg_error, parse_args};
use json::{FixpointRecord, MurmRecord, PassRecord, Recording, RunInfo};
use progress::{box_lines, fmt_num, fmt_size};
use ui::{Ui, Verbosity};

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
            "  {}\n\t{} {} → {} gates · {} 2q gates · {} T/Tdg{} · {} depth · {:.3}s",
            name,
            self.ui.elbow(true),
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
        self.ui.info(&format!(
            "\t{} Synthesis basis: {basis}",
            self.ui.elbow(true)
        ));
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

/// A heading (none if empty) with its figures on `├─`/`└─` lines below.
fn tree(ui: &Ui, heading: &str, lines: &[String]) {
    if !heading.is_empty() {
        ui.info(&format!("  {heading}"));
    }
    for (i, line) in lines.iter().enumerate() {
        ui.info(&format!("\t{} {line}", ui.elbow(i + 1 == lines.len())));
    }
}

/// `1 gate`, `2 gates`: a count with thousands separators and a plural noun.
fn count(n: usize, noun: &str) -> String {
    if n == 1 {
        format!("1 {noun}")
    } else {
        format!("{} {noun}s", fmt_num(n))
    }
}

/// Operation counts and axis weights of a PBC circuit, and whether the
/// weights are known (they are all zero when materializing exceeded its
/// budget).
#[derive(Clone, Copy)]
struct PbcStats {
    axes: AxisWeights,
    weighed: bool,
}

impl PbcStats {
    /// By materializing the axes; counts only when that exceeds its budget.
    fn materialize(pbc: &tzap::pbc::PbcCircuit, budget: usize) -> Self {
        match pbc.axis_weight_summary(budget) {
            Ok(axes) => Self {
                axes,
                weighed: true,
            },
            Err(_) => Self {
                axes: AxisWeights::of(
                    pbc.operations()
                        .iter()
                        .filter_map(|op| Some((AxisKind::of(op)?, 0))),
                ),
                weighed: false,
            },
        }
    }

    /// From exported PBC text, which lists every axis already expanded:
    /// each `r` or `m` line's weight is its number of factors.
    fn from_text(text: &str) -> Self {
        let axes = AxisWeights::of(text.lines().filter_map(|line| {
            let mut fields = line.split_whitespace();
            match fields.next()? {
                "r" => {
                    let k: i64 = fields.next()?.parse().ok()?;
                    let kind = match k.rem_euclid(8) {
                        0 => return None,
                        k if k % 2 == 1 => AxisKind::Pi8,
                        _ => AxisKind::Clifford,
                    };
                    Some((kind, fields.skip(1).count()))
                }
                "m" => Some((
                    AxisKind::Measurement,
                    fields.skip(1).take_while(|f| *f != "->").count(),
                )),
                _ => None,
            }
        }));
        Self {
            axes,
            weighed: true,
        }
    }

    /// `N π/8 rotations · N Clifford rotations · N measurements`, omitting
    /// empty Clifford and measurement counts.
    fn counts(&self) -> String {
        let a = &self.axes;
        let mut parts = vec![count(a.pi8.count, "π/8 rotation")];
        if a.clifford.count > 0 {
            parts.push(count(a.clifford.count, "Clifford rotation"));
        }
        if a.measurements.count > 0 {
            parts.push(count(a.measurements.count, "measurement"));
        }
        parts.join(" · ")
    }

    /// One `π/8 weight min/median/max: 1/2/3` line per non-empty kind; none
    /// when weights are unknown.
    fn weights(&self, measurements: bool) -> Vec<String> {
        if !self.weighed {
            return Vec::new();
        }
        let a = &self.axes;
        let mut kinds = vec![("π/8", a.pi8), ("Clifford", a.clifford)];
        if measurements {
            kinds.push(("measurement", a.measurements));
        }
        kinds
            .into_iter()
            .filter(|(_, w)| w.count > 0)
            .map(|(label, w)| {
                format!(
                    "{label} weight min/median/max: {}/{}/{}",
                    w.min,
                    fmt_median(w.median()),
                    w.max
                )
            })
            .collect()
    }

    fn json(&self) -> json::PbcCounts {
        json::PbcCounts {
            axes: self.axes,
            weighed: self.weighed,
        }
    }
}

/// A median, which is a whole number or halfway between two.
fn fmt_median(median: f64) -> String {
    if median.fract() == 0.0 {
        format!("{median:.0}")
    } else {
        format!("{median:.1}")
    }
}

fn toffolis(circuit: &Circuit) -> usize {
    circuit
        .gates
        .iter()
        .filter(|g| matches!(g, Gate::ccx { .. } | Gate::ccz { .. }))
        .count()
}

/// Reject input PBC conversion cannot take before any work is done, naming
/// the offending instruction: resets, gates repeating a qubit, and Rz gates
/// that no requested pass decomposes.
fn check_pbc_input(ui: &Ui, run: &Run, circuit: &Circuit) {
    let decomposes_rz = run.options.decompose_rotations_enabled()
        || run
            .options
            .passes
            .as_ref()
            .is_some_and(|p| p.contains(&tzap::optimize::PassName::DecomposeRz));
    for (index, gate) in circuit.gates.iter().enumerate() {
        let problem = match gate {
            Gate::reset(_) => "PBC has no resets",
            Gate::rz(..)
            | Gate::p(..)
            | Gate::rx(..)
            | Gate::ry(..)
            | Gate::cp { .. }
            | Gate::crx { .. }
            | Gate::cry { .. }
            | Gate::crz { .. }
                if !decomposes_rz =>
            {
                "Rotations need --decompose-rotations (or DecomposeRotations in --passes) before PBC conversion"
            }
            _ => {
                let mut qubits = tzap::circuit::qubits_of(gate);
                qubits.sort_unstable();
                if qubits.windows(2).any(|w| w[0] == w[1]) {
                    "a gate's qubits must be distinct"
                } else {
                    continue;
                }
            }
        };
        ui.abort(&format!(
            "Error converting to PBC: gate {} ({gate}): {problem}",
            fmt_num(index + 1)
        ));
    }
}

/// What a run produces: the text for `-o`, an SVG drawing to write after
/// it, and the PBC report for `--json`.
struct Prepared {
    output: Option<String>,
    drawing: Option<(String, String)>,
    pbc: Option<json::PbcReport>,
}

/// Prepare the final representation before reporting success or writing a file.
fn prepare_output(ui: &Ui, run: &Run, circuit: &Circuit) -> Prepared {
    if !(run.to_pbc || run.visualize_pbc.is_some()) {
        let output = run.output_path.as_ref().map(|_| circuit.to_qasm());
        return Prepared {
            output,
            drawing: None,
            pbc: None,
        };
    }
    // Conversion is terminal: optimization and requested decompositions have
    // already finished. Validate even when no output destination was requested.
    let convert_start = Instant::now();
    let mut pbc = tzap::pbc::to_pbc(circuit, run.pbc_max_weight)
        .unwrap_or_else(|e| ui.abort(&format!("Error converting to PBC: {e}")));
    let convert_seconds = convert_start.elapsed().as_secs_f64();
    // The optimizer reports weights before and after from its own packed
    // axes; otherwise they come from the exported text. Only a drawing
    // without either materializes the axes just for this report.
    let mut optimized = None;
    if run.pbc_opt {
        let start = Instant::now();
        // Moving Cliffords into the frame is the only step that widens
        // axes; merges and swaps keep every axis, so a bound still holds.
        let options = tzap::pbc::OptimizeOptions {
            clifford_to_frame: run.pbc_max_weight.is_none(),
            ..Default::default()
        };
        match pbc.optimize_rotations(options) {
            Ok(stats) => optimized = Some((stats, start.elapsed().as_secs_f64())),
            // The pass leaves the circuit unchanged on failure.
            Err(e) => ui.note(&format!("  PBC optimization skipped: {e}")),
        }
    }
    let output = if run.to_pbc {
        Some(
            pbc.to_text_with(tzap::pbc::TextOptions {
                max_expansion_cells: run.pbc_expansion_budget,
            })
            .unwrap_or_else(|e| {
                ui.abort(&format!(
                    "Error exporting PBC: {e} ({} work units; raise it with \
                     --pbc-expansion-budget)",
                    fmt_num(run.pbc_expansion_budget)
                ))
            }),
        )
    } else {
        run.output_path.as_ref().map(|_| circuit.to_qasm())
    };
    let stats = |axes| PbcStats {
        axes,
        weighed: true,
    };
    let before = match (&optimized, &output) {
        (Some((s, _)), _) => stats(s.axes_before),
        (None, Some(text)) if run.to_pbc => PbcStats::from_text(text),
        _ => PbcStats::materialize(&pbc, run.pbc_expansion_budget),
    };
    let mut lines = vec![before.counts()];
    lines.extend(before.weights(true));
    tree(
        ui,
        &format!("Converted to PBC in {convert_seconds:.3}s"),
        &lines,
    );
    let optimization = optimized.map(|(s, seconds)| {
        let after = stats(s.axes_after);
        let (b, a) = (&before.axes, &after.axes);
        let t = (b.pi8.count, a.pi8.count);
        let clifford = (b.clifford.count, a.clifford.count);
        if t.0 == t.1 && clifford.0 == clifford.1 {
            tree(
                ui,
                &format!("Optimized PBC in {seconds:.3}s"),
                &["no reduction".to_string()],
            );
        } else {
            // `a → b N π/8 rotations (↓x%)`.
            let change = |(before, after): (usize, usize), noun: &str| {
                let reduction = if before > 0 {
                    (before as f64 - after as f64) / before as f64 * 100.0
                } else {
                    0.0
                };
                let noun = if after == 1 {
                    noun.to_string()
                } else {
                    format!("{noun}s")
                };
                format!(
                    "{} → {} {noun} (↓{reduction:.1}%)",
                    fmt_num(before),
                    fmt_num(after)
                )
            };
            let mut lines = vec![change(t, "π/8 rotation")];
            if clifford != (0, 0) {
                lines.push(change(clifford, "Clifford rotation"));
            }
            lines.extend(after.weights(false));
            tree(ui, &format!("Optimized PBC in {seconds:.3}s"), &lines);
        }
        json::PbcOptimization {
            output: after.json(),
            merges: s.merges,
            mcr_swaps: s.swaps,
            cliffords_to_frame: s.cliffords_to_frame,
            seconds,
        }
    });
    let drawing = run.visualize_pbc.as_ref().map(|path| {
        let options = tzap::pbc::SvgOptions {
            // Every operation, however large the circuit.
            max_operations: usize::MAX,
            max_expansion_cells: run.pbc_expansion_budget,
            ..Default::default()
        };
        let svg = pbc
            .to_svg_with(options)
            .unwrap_or_else(|e| ui.abort(&format!("Error drawing PBC: {e}")));
        (path.clone(), svg)
    });
    let report = json::PbcReport {
        max_weight: run.pbc_max_weight.map(|w| w.get()),
        converted: before.json(),
        convert_seconds,
        optimization,
    };
    Prepared {
        output,
        drawing,
        pbc: Some(report),
    }
}

/// Write to the requested file or stdout; diagnostics remain on stderr.
fn write_output(ui: &Ui, run: &Run, output: Option<String>) {
    let Some(path) = &run.output_path else {
        return;
    };
    let output = output.expect("requested output was prepared");
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
/// With `metrics`, also report the circuit's gate metrics (for runs that
/// convert to PBC with no gate-level pass, which report nothing else).
fn read_circuit(ui: &Ui, path: &str, metrics: bool) -> Parsed {
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
        "\t{} {} · {}",
        ui.elbow(false),
        count(circuit.num_qubits, "qubit"),
        count(circuit.gates.len(), "gate"),
    ));
    if metrics {
        let m = Metrics::of(&circuit);
        let mut parts = vec![
            count(m.two_qubit, "2q gate"),
            format!("{} T/Tdg", fmt_num(m.t)),
        ];
        if m.rz > 0 {
            parts.push(format!("{} Rz", fmt_num(m.rz)));
        }
        let toffolis = toffolis(&circuit);
        if toffolis > 0 {
            // Each carries seven T gates, which the T/Tdg count leaves out.
            parts.push(format!("{} CCX/CCZ", fmt_num(toffolis)));
        }
        parts.push(format!("{} depth", fmt_num(m.depth)));
        ui.info(&format!("\t{} {}", ui.elbow(false), parts.join(" · ")));
    }
    ui.info(&format!(
        "\t{} Circuit gates: {}",
        ui.elbow(true),
        circuit.gate_set()
    ));
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
/// Report the result and write the outputs: the main output first, so a
/// failing drawing path cannot lose it. Returns the PBC report for `--json`.
fn finish(
    ui: &Ui,
    report: &Report,
    result: &Circuit,
    run: &Run,
    start: Instant,
) -> Option<json::PbcReport> {
    // With PBC output, the gate pipeline was already summarized in one line.
    if !run.to_pbc {
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
    }
    let prepared = prepare_output(ui, run, result);
    write_output(ui, run, prepared.output);
    if let Some((path, svg)) = &prepared.drawing {
        fs::write(path, svg).unwrap_or_else(|e| ui.abort(&format!("Error writing {path}: {e}")));
        ui.info(&format!("  wrote PBC drawing {path}"));
    }
    if run.to_pbc {
        ui.info(&format!("  Done in {:.3}s", start.elapsed().as_secs_f64()));
    }
    prepared.pbc
}

fn main() {
    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    let opts = parse_args(&args);
    let Opts { action, ui, json } = opts;

    let mut run = match action {
        Action::CacheInfo => return print_cache_info(&ui, json),
        Action::ClearCache => return clear_cache(&ui, json),
        Action::Optimize(run) => *run,
    };

    ui.info(&format!(
        "{}⚡\u{FE0F} tzap{} {}v{}{}",
        ui.sgr("\x1b[1m"),
        ui.reset(),
        ui.sgr("\x1b[2m"),
        env!("CARGO_PKG_VERSION"),
        ui.reset()
    ));
    let bare_pbc = run.to_pbc && run.options.passes.as_ref().is_some_and(Vec::is_empty);
    let parsed = read_circuit(&ui, &run.input_path, bare_pbc);
    let preserved = parsed
        .circuit
        .gates
        .iter()
        .filter(|g| matches!(g, Gate::rz(a, _) if a.is_preserved()))
        .count();
    if preserved != 0 {
        ui.note(&format!("Warning: {preserved} angle expressions retained unchanged because exact representation or folding exceeded its limits"));
    }
    // Said out loud only when the size decided it: a run nobody asked to
    // parallelize otherwise reaches the chunked progress box with no
    // explanation, and a slightly different gate count at the end.
    if run.to_pbc || run.visualize_pbc.is_some() {
        check_pbc_input(&ui, &run, &parsed.circuit);
    }
    if run.resolve_parallel(parsed.circuit.gates.len()) && !run.gate_opts_off {
        ui.note(&format!(
            "  Optimizing chunks in parallel ({}+ gates) · --no-parallel to disable",
            fmt_num(AUTO_PARALLEL_GATES)
        ));
    }
    if run.gate_opts_off {
        let instead = if run.pbc_opt {
            " (the PBC is optimized instead)"
        } else {
            ""
        };
        ui.note(&format!(
            "  Gate-level optimization is off with --to-pbc{instead}; pass -O3 or --passes to \
             enable it"
        ));
        for flag in &run.ignored {
            ui.note(&format!(
                "  {flag} has no effect without gate-level optimization"
            ));
        }
    }
    // With PBC output, the gate pipeline is a preliminary: it runs silently
    // and is summarized in one line, keeping the report on the PBC. The
    // observer still records for --json.
    let (observer, pbc_ui) = if run.to_pbc {
        (Terminal::new(Ui::new(Verbosity::Quiet), json), Some(ui))
    } else {
        (Terminal::new(ui, json), None)
    };
    let ui = pbc_ui.as_ref().unwrap_or(&observer.ui);
    let (result, report) = if run.options.passes.as_ref().is_some_and(Vec::is_empty) {
        // `--to-pbc` alone with no decompositions: nothing for the gate
        // pipeline to do.
        let metrics = Metrics::of(&parsed.circuit);
        let report = Report {
            input: metrics,
            baseline: metrics,
            output: metrics,
            numerical: tzap::angle_stats::NumericalReport::input(&parsed.circuit, &run.options),
        };
        (parsed.circuit.clone(), report)
    } else if run.to_pbc {
        let (what, pending) = if run.gate_opts_off {
            ("Decomposition".to_string(), "Decomposing".to_string())
        } else if run.options.passes.is_some() {
            (
                "Gate-level passes".to_string(),
                "Running gate-level passes".to_string(),
            )
        } else {
            let level = format!("{:?}", run.options.level);
            (
                format!("Gate-level {level}"),
                format!("Optimizing gate-level circuit ({level})"),
            )
        };
        let gate_start = Instant::now();
        ui.start_inline(&format!("  {pending}..."));
        let (result, report) = optimize_with(&parsed.circuit, &run.options, &observer)
            .unwrap_or_else(|e| arg_error(e));
        let (b, o) = (report.baseline, report.output);
        let arrow = |before: usize, after: usize, noun: &str| {
            let noun = if after == 1 {
                noun.to_string()
            } else {
                format!("{noun}s")
            };
            format!("{} → {} {noun}", fmt_num(before), fmt_num(after))
        };
        let mut parts = vec![
            arrow(b.gates, o.gates, "gate"),
            arrow(b.two_qubit, o.two_qubit, "2q gate"),
            format!("{} → {} T/Tdg", fmt_num(b.t), fmt_num(o.t)),
        ];
        if b.rz > 0 || o.rz > 0 {
            parts.push(format!("{} → {} Rz", fmt_num(b.rz), fmt_num(o.rz)));
        }
        let (tb, to) = (toffolis(&parsed.circuit), toffolis(&result));
        if tb > 0 || to > 0 {
            parts.push(format!("{} → {} CCX/CCZ", fmt_num(tb), fmt_num(to)));
        }
        parts.push(format!("{} → {} depth", fmt_num(b.depth), fmt_num(o.depth)));
        ui.finish_inline(&format!(
            "  {what} in {:.3}s",
            gate_start.elapsed().as_secs_f64()
        ));
        tree(ui, "", &parts);
        (result, report)
    } else {
        optimize_with(&parsed.circuit, &run.options, &observer).unwrap_or_else(|e| arg_error(e))
    };
    if report.numerical.uncertified_syntheses != 0 {
        ui.note(&format!(
            "{} Rz syntheses used numeric targets; epsilon applies per rotation, and total conversion and circuit error are uncertified",
            report.numerical.uncertified_syntheses
        ));
    }
    let pbc = finish(ui, &report, &result, &run, start);

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
            gate_opts_off: run.gate_opts_off,
            output_is_pbc: run.to_pbc,
            pbc_passes: &run.pbc_passes,
            pbc: pbc.as_ref(),
        };
        let recording = observer.take_recording();
        ui.write_stdout(&json::render(&info, &run.options, &report, &recording));
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
            Action::Optimize(run) => *run,
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
