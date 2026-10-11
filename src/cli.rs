//! Command-line argument parsing, `Opts`, and `--help` text. The pass and
//! optimization-level enums, and everything they select, live in
//! `tzap::optimize` — this module only maps flags onto an
//! [`Options`](tzap::optimize::Options).

use std::process;

use tzap::optimize::{DEFAULT_RZ_EPSILON, Level, Options, PassName, SuperOptBounds, SuperOptGates};

use crate::ui::{Ui, Verbosity};

/// The path spelling that means "the standard stream" rather than a file, for
/// both the input and the output — the usual Unix convention, so tzap can sit
/// in the middle of a pipeline.
pub(crate) const STREAM_PATH: &str = "-";

/// Default `--pbc-expansion-budget`: sparse work for expanding PBC axes into
/// explicit Pauli strings, a guard against runaway memory (roughly 8 bytes a
/// unit at worst). The largest benchmark, gf2^256_mult, needs about 10M.
pub(crate) const DEFAULT_PBC_EXPANSION_BUDGET: usize = 256_000_000;

/// Passes that run on the PBC circuit, after the gate passes: `(name,
/// description)`. `ToPbc` converts; `PbcOpt` then optimizes the rotations.
pub(crate) const PBC_PASSES: [(&str, &str); 2] = [
    (
        "ToPbc",
        "Convert to PBC; outputs PBC (see --pbc-max-weight). Must follow the gate passes",
    ),
    (
        "PbcOpt",
        "Merge same-axis PBC rotations to lower the T count. Must follow ToPbc",
    ),
];

fn is_pass_name(name: &str) -> bool {
    PassName::parse(name).is_some() || PBC_PASSES.iter().any(|(n, _)| *n == name)
}

/// A parsed `--passes` list: gate passes, then optional PBC passes.
struct PassList {
    gates: Vec<PassName>,
    to_pbc: bool,
    pbc_opt: bool,
}

fn parse_pass_list(list: &str) -> PassList {
    let mut parsed = PassList {
        gates: Vec::new(),
        to_pbc: false,
        pbc_opt: false,
    };
    let names: Vec<&str> = list
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    if names.is_empty() {
        arg_error(
            "--passes requires at least one pass name \
             (e.g. --passes CancelGates,PhaseFoldRand)",
        );
    }
    for name in names {
        match name {
            "ToPbc" if parsed.to_pbc => arg_error("--passes lists ToPbc more than once"),
            "ToPbc" => parsed.to_pbc = true,
            "PbcOpt" if !parsed.to_pbc => {
                arg_error("PbcOpt optimizes PBC — list it after ToPbc in --passes")
            }
            "PbcOpt" if parsed.pbc_opt => arg_error("--passes lists PbcOpt more than once"),
            "PbcOpt" => parsed.pbc_opt = true,
            _ => {
                let pass = PassName::parse(name).unwrap_or_else(|| {
                    arg_error(format!(
                        "Unknown pass '{name}'. Available passes: {}, {}",
                        PassName::all_names(),
                        PBC_PASSES.map(|(n, _)| n).join(", ")
                    ))
                });
                if parsed.to_pbc {
                    arg_error(format!(
                        "{name} is a gate pass and cannot run after ToPbc — list gate passes first"
                    ));
                }
                parsed.gates.push(pass);
            }
        }
    }
    parsed
}

fn looks_like_pass_list_fragment(token: &str) -> bool {
    token
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .all(is_pass_name)
}

/// What this invocation is for. Optimizing a circuit is the whole point of
/// tzap; the cache actions are maintenance on the artifact an `-Osuper` run
/// leaves behind, and take no input circuit.
pub(crate) enum Action {
    Optimize(Box<Run>),
    /// `--cache-info`: report where the on-disk MURMs live and
    /// what they cost.
    CacheInfo,
    /// `--clear-cache`: delete them.
    ClearCache,
}

/// Gate count at and above which a run goes parallel on its own. Map-reduce
/// chunking costs a little quality — chunk boundaries are walls no pass can
/// optimize across, worth a few hundredths of a percent of the gate count —
/// and buys several-fold wall-clock. That trade is the wrong way round for
/// the circuits people iterate on and clearly right for the ones where a
/// serial run takes tens of seconds, so the default follows the size.
///
/// Measured against the circuit *as parsed*, which is the number tzap has
/// just printed and the one a user would predict from the file. A later
/// opt-in decomposition does not retroactively change this decision.
pub(crate) const AUTO_PARALLEL_GATES: usize = 1_000_000;

/// One optimization run: the file paths the CLI itself owns, plus the
/// optimizer configuration to hand to `tzap::optimize`.
pub(crate) struct Run {
    /// The input circuit's path, or [`STREAM_PATH`] for stdin.
    pub(crate) input_path: String,
    /// Where to write the optimized circuit: a path, [`STREAM_PATH`] for
    /// stdout, or `None` to discard it.
    pub(crate) output_path: Option<String>,
    pub(crate) to_pbc: bool,
    /// Run the PBC rotation optimizer after conversion: on by default,
    /// off with `--pbc-no-opt`; under `--passes … ToPbc`, only if `PbcOpt` is
    /// listed.
    pub(crate) pbc_opt: bool,
    /// Gate-level optimization was turned off because `--to-pbc` was given
    /// without `-O*` or `--passes`; only requested decompositions run.
    pub(crate) gate_opts_off: bool,
    /// Flags given that have no effect because gate-level optimization is off.
    pub(crate) ignored: Vec<&'static str>,
    /// The PBC passes listed in `--passes`, for the `--json` report.
    pub(crate) pbc_passes: Vec<&'static str>,
    /// Write an SVG drawing of the PBC circuit here (`--visualize-pbc`).
    pub(crate) visualize_pbc: Option<String>,
    /// Sparse-expansion work allowed for writing, drawing, or measuring PBC
    /// axes (`--pbc-expansion-budget`).
    pub(crate) pbc_expansion_budget: usize,
    /// Bound on PBC rotation and measurement weight (`--pbc-max-weight`).
    pub(crate) pbc_max_weight: Option<std::num::NonZeroUsize>,
    /// `--parallel`/`--no-parallel` as asked for, or `None` to decide from
    /// the circuit's size (see [`Run::resolve_parallel`]). Distinct from
    /// `options.parallel`, which is the answer rather than the request.
    parallel: Option<bool>,
    pub(crate) options: Options,
}

impl Run {
    pub(crate) fn reads_stdin(&self) -> bool {
        self.input_path == STREAM_PATH
    }

    pub(crate) fn writes_stdout(&self) -> bool {
        self.output_path.as_deref() == Some(STREAM_PATH)
    }

    /// Settle whether this run is parallel, now that the circuit has been
    /// read and `gates` is known, and record it in `options` so the
    /// optimizer, the progress box, and the `--json` report all describe the
    /// same run. Returns whether the size alone decided it, which is the one
    /// case worth mentioning on the way past: nobody asked for parallel, so
    /// nobody expects the chunked progress box.
    pub(crate) fn resolve_parallel(&mut self, gates: usize) -> bool {
        let by_size = self.parallel.is_none() && gates >= AUTO_PARALLEL_GATES;
        self.options.parallel = self.parallel.unwrap_or(by_size);
        by_size
    }
}

/// Parsed command-line options.
pub(crate) struct Opts {
    pub(crate) action: Action,
    /// Terminal capabilities and verbosity for this run, decided once from
    /// `--quiet` and the streams themselves.
    pub(crate) ui: Ui,
    /// `--json`: write a machine-readable report of the run to stdout.
    pub(crate) json: bool,
}

/// Print `Error: {msg}` and exit 1. The single entry point for every
/// argument-parsing failure, so all CLI errors share one unmistakable
/// prefix instead of some being phrased as errors and others not.
pub(crate) fn arg_error(msg: impl std::fmt::Display) -> ! {
    eprintln!("Error: {msg}");
    process::exit(1);
}

/// Parse the next argument as a `usize`, exiting with `flag_name` in the
/// error message on failure, if there is no next argument, or if it parses
/// to 0 — every caller of this (the hidden `--superopt-*` bounds) feeds a
/// count or width that must be at least 1, and the message already promises
/// "positive integer", so 0 must be rejected too rather than silently
/// accepted as a valid `usize`.
fn parse_usize_arg(args: &[String], i: usize, flag_name: &str) -> usize {
    let value = args
        .get(i)
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| arg_error(format!("{flag_name} requires a positive integer")));
    if value == 0 {
        arg_error(format!("{flag_name} requires a positive integer, got 0"));
    }
    value
}

/// Take the next argument as a string value for `flag_name`, exiting if
/// there isn't one or if it's empty.
fn parse_string_arg(args: &[String], i: usize, flag_name: &str, what: &str) -> String {
    match args.get(i) {
        Some(value) if !value.is_empty() => value.clone(),
        _ => arg_error(format!("{flag_name} requires {what}")),
    }
}

/// Every long flag that takes a separate value, and so may also be written
/// `--flag=value` (see [`split_flag_values`]).
const VALUE_FLAGS: [&str; 10] = [
    "--epsilon",
    "--pbc-expansion-budget",
    "--pbc-max-weight",
    "--visualize-pbc",
    "--passes",
    "--cache-dir",
    "--superopt-qubits",
    "--superopt-window-gates",
    "--superopt-murm-entries",
    "--superopt-gates",
];

/// Rewrite `--flag=value` into two arguments, so every value-taking long
/// flag accepts both spellings — `--epsilon=1e-6` is as natural to type as
/// `--epsilon 1e-6`, and a parser that took only one of the two forms would
/// reject what someone reaches for first.
///
/// Only the known [`VALUE_FLAGS`] are split, so nothing else containing an
/// `=` is disturbed: a pass list, or a file whose name happens to contain
/// one, arrives untouched.
fn split_flag_values(args: &[String]) -> Vec<String> {
    let mut out = Vec::with_capacity(args.len());
    for arg in args {
        let split = arg
            .split_once('=')
            .filter(|(flag, _)| VALUE_FLAGS.contains(flag));
        match split {
            Some((flag, value)) => {
                out.push(flag.to_string());
                out.push(value.to_string());
            }
            None => out.push(arg.clone()),
        }
    }
    out
}

pub(crate) fn parse_args(args: &[String]) -> Opts {
    let args = split_flag_values(args);
    let mut input_path: Option<String> = None;
    let mut output_path: Option<String> = None;
    let mut to_pbc = false;
    let mut pbc_no_opt = false;
    let mut visualize_pbc: Option<String> = None;
    let mut pbc_expansion_budget: Option<usize> = None;
    let mut pbc_max_weight: Option<std::num::NonZeroUsize> = None;
    let mut decompose_rotations = false;
    let mut decompose_cz = false;
    let mut decompose_ccx = false;
    let mut rz_epsilon: f64 = DEFAULT_RZ_EPSILON;
    let mut parallel: Option<bool> = None;
    let mut passes: Option<Vec<PassName>> = None;
    // (ToPbc, PbcOpt) listed in --passes.
    let mut pbc_passes = (false, false);
    let mut fixpoint = false;
    let mut optimization_level = None;
    let mut superopt_qubits: Option<usize> = None;
    let mut superopt_window_gates: Option<usize> = None;
    let mut superopt_murm_entries: Option<usize> = None;
    let mut superopt_gates = SuperOptGates::Auto;
    let mut quiet = false;
    let mut json = false;
    let mut cache_dir: Option<String> = None;
    let mut cache_info = false;
    let mut clear_cache = false;
    // Help and version are honored after the whole line is parsed rather than
    // on the spot, so a flag that changes how they print — `--quiet`, today —
    // works on either side of them.
    let mut help = false;
    let mut version = false;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--help" | "-h" => help = true,
            "--version" | "-v" | "-V" => version = true,
            "--decompose-rotations" | "--decompose-rz" => decompose_rotations = true,
            "--to-pbc" => to_pbc = true,
            "--pbc-no-opt" => pbc_no_opt = true,
            "--decompose-cz" => decompose_cz = true,
            "--decompose-ccx" => decompose_ccx = true,
            "--superopt-gates" => {
                i += 1;
                let value =
                    parse_string_arg(&args, i, "--superopt-gates", "a basis mode or gate list");
                superopt_gates =
                    SuperOptGates::parse(&value).unwrap_or_else(|error| arg_error(error));
            }
            "--epsilon" => {
                i += 1;
                let value: f64 = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_else(|| arg_error("--epsilon requires a number (e.g. 1e-10)"));
                if !(value.is_finite() && value > 0.0) {
                    arg_error(format!(
                        "--epsilon must be a positive, finite number, got {value} \
                         (e.g. 1e-10) — zero or negative values make Rz synthesis undefined"
                    ));
                }
                rz_epsilon = value;
            }
            "--passes" => {
                i += 1;
                let list = args.get(i).unwrap_or_else(|| {
                    arg_error(
                        "--passes requires a comma-separated list of pass names \
                         (e.g. --passes CancelGates,PhaseFoldRand)",
                    )
                });
                let mut list = list.clone();
                while let Some(next) = args.get(i + 1) {
                    if next.starts_with('-') || !looks_like_pass_list_fragment(next) {
                        break;
                    }
                    list.push(',');
                    list.push_str(next);
                    i += 1;
                }
                let parsed = parse_pass_list(&list);
                passes = Some(parsed.gates);
                pbc_passes = (parsed.to_pbc, parsed.pbc_opt);
            }
            // Last spelling on the line wins, as for any --x/--no-x pair, so
            // a wrapper script's default can be overridden by appending to it.
            "--parallel" => parallel = Some(true),
            "--no-parallel" => parallel = Some(false),
            "--fixpoint" => fixpoint = true,
            "-q" | "--quiet" => quiet = true,
            "--json" => json = true,
            "--cache-dir" => {
                i += 1;
                cache_dir = Some(parse_string_arg(
                    &args,
                    i,
                    "--cache-dir",
                    "a directory path",
                ));
            }
            "--cache-info" => cache_info = true,
            "--clear-cache" => clear_cache = true,
            "-O1" | "-O2" | "-O3" | "-Osuper" => {
                if optimization_level.is_some() {
                    arg_error("-O1, -O2, -O3, and -Osuper cannot be combined — pick exactly one");
                }
                optimization_level = Some(match args[i].as_str() {
                    "-O1" => Level::O1,
                    "-O2" => Level::O2,
                    "-O3" => Level::O3,
                    "-Osuper" => Level::Osuper,
                    _ => unreachable!(),
                });
            }
            "--visualize-pbc" => {
                i += 1;
                visualize_pbc = Some(args.get(i).cloned().unwrap_or_else(|| {
                    arg_error("--visualize-pbc requires an output SVG file path")
                }));
            }
            "--pbc-expansion-budget" => {
                i += 1;
                let value = args.get(i).map_or("", |s| s.as_str());
                pbc_expansion_budget = Some(value.parse().unwrap_or_else(|_| {
                    arg_error(format!(
                        "--pbc-expansion-budget requires a number of work units, got '{value}'"
                    ))
                }));
            }
            "--pbc-max-weight" => {
                i += 1;
                let value = args.get(i).map_or("", |s| s.as_str());
                pbc_max_weight = Some(value.parse().unwrap_or_else(|_| {
                    arg_error(format!(
                        "--pbc-max-weight requires an integer of at least 1, got '{value}'"
                    ))
                }));
            }
            "-o" => {
                i += 1;
                output_path = Some(
                    args.get(i)
                        .cloned()
                        .unwrap_or_else(|| arg_error("-o requires an output file path")),
                );
            }
            // Hidden: not listed in --help, for experimentation with SuperOpt's
            // window/MURM bounds without a rebuild.
            "--superopt-qubits" => {
                i += 1;
                superopt_qubits = Some(parse_usize_arg(&args, i, "--superopt-qubits"));
            }
            "--superopt-window-gates" => {
                i += 1;
                superopt_window_gates = Some(parse_usize_arg(&args, i, "--superopt-window-gates"));
            }
            "--superopt-murm-entries" => {
                i += 1;
                superopt_murm_entries = Some(parse_usize_arg(&args, i, "--superopt-murm-entries"));
            }
            // A bare "-" is a positional, not a flag: stdin as the input,
            // stdout as the output. Checked before the unknown-flag arm,
            // which every other leading dash falls into.
            _ if args[i].starts_with('-') && args[i] != STREAM_PATH => {
                arg_error(format!(
                    "unknown flag '{}'. Run `tzap --help` for the list of valid options",
                    args[i]
                ));
            }
            _ => {
                if input_path.is_none() {
                    input_path = Some(args[i].clone());
                } else if output_path.is_none() {
                    output_path = Some(args[i].clone());
                } else {
                    arg_error(format!(
                        "unexpected extra argument '{}' — tzap takes at most \
                         <input.qasm> and [output.qasm]",
                        args[i]
                    ));
                }
            }
        }
        i += 1;
    }

    let verbosity = if quiet {
        Verbosity::Quiet
    } else {
        Verbosity::Normal
    };
    let ui = Ui::new(verbosity);

    if help {
        print_help(&ui);
        process::exit(0);
    }
    if version {
        ui.write_stdout(&format!("tzap {}\n", env!("CARGO_PKG_VERSION")));
        process::exit(0);
    }

    // Installed before anything can touch a MURM, so `--cache-dir` governs
    // this whole process — including the cache actions below, which report on
    // and clear whichever directory is in force.
    if let Some(dir) = &cache_dir
        && let Err(existing) = tzap::super_opt::set_cache_dir(dir.into())
    {
        arg_error(format!(
            "the cache directory is already fixed at {} for this process",
            existing.display()
        ));
    }

    if cache_info && clear_cache {
        arg_error("--cache-info and --clear-cache cannot be combined — pick exactly one");
    }
    if cache_info || clear_cache {
        let flag = if cache_info {
            "--cache-info"
        } else {
            "--clear-cache"
        };
        if let Some(path) = input_path {
            arg_error(format!(
                "{flag} takes no input circuit, but '{path}' was given — \
                 it only reports on or clears tzap's on-disk cache"
            ));
        }
        return Opts {
            action: if cache_info {
                Action::CacheInfo
            } else {
                Action::ClearCache
            },
            ui,
            json,
        };
    }

    let Some(input_path) = input_path else {
        arg_error(
            "missing required <input.qasm> argument\n\n  \
             Usage: tzap <input.qasm> [-o output.qasm] [-O1|-O2|-O3|-Osuper] \
             [--decompose-ccx] [--decompose-cz] [--decompose-rotations] [--passes <list>] \
             [--parallel|--no-parallel] [--fixpoint]\n  \
             Pass - to read the circuit from stdin.\n  \
             Run `tzap --help` for the full option list.",
        );
    };

    if optimization_level.is_some() && (passes.is_some() || fixpoint) {
        arg_error("-O1, -O2, -O3, and -Osuper cannot be combined with --passes or --fixpoint");
    }
    if passes.is_some() && (decompose_rotations || decompose_cz || decompose_ccx) {
        arg_error(
            "--passes cannot be combined with --decompose-rotations, --decompose-cz, or --decompose-ccx \
             — list the corresponding decomposition passes instead",
        );
    }
    // --to-pbc is shorthand for ending the pipeline with ToPbc,PbcOpt (ToPbc
    // alone with --pbc-no-opt); mixing the two spellings would make the order
    // ambiguous.
    if pbc_passes.0 && to_pbc {
        arg_error("--passes already lists ToPbc — drop --to-pbc");
    }
    if pbc_passes.0 && pbc_no_opt {
        arg_error(
            "--pbc-no-opt does not apply to --passes — list PbcOpt after ToPbc to optimize, \
             or leave it out",
        );
    }
    if to_pbc && fixpoint && passes.is_none() {
        arg_error(
            "--fixpoint repeats gate-level passes, which --to-pbc turns off — use -O3 \
             (a fixpoint) or --passes <list> --fixpoint",
        );
    }
    if pbc_no_opt && !to_pbc && visualize_pbc.is_none() {
        arg_error(
            "--pbc-no-opt applies to PBC output — combine it with --to-pbc or --visualize-pbc",
        );
    }
    let pbc_opt = if pbc_passes.0 {
        pbc_passes.1
    } else {
        !pbc_no_opt
    };
    // Gate-level passes tend to widen PBC axes and cost time, and PBC has its
    // own optimizer: --to-pbc alone runs only the requested decompositions.
    let gate_opts_off = to_pbc && optimization_level.is_none() && passes.is_none() && !fixpoint;
    if gate_opts_off {
        let requested = [
            (decompose_ccx, PassName::DecomposeToffoli),
            (decompose_cz, PassName::DecomposeCz),
            (decompose_rotations, PassName::DecomposeRotations),
        ];
        passes = Some(requested.iter().filter(|r| r.0).map(|r| r.1).collect());
        (decompose_ccx, decompose_cz, decompose_rotations) = (false, false, false);
    }
    to_pbc |= pbc_passes.0;
    if pbc_expansion_budget.is_some() && !to_pbc && visualize_pbc.is_none() {
        arg_error(
            "--pbc-expansion-budget applies to PBC output — combine it with --to-pbc or \
             --visualize-pbc",
        );
    }
    if pbc_max_weight.is_some() && !to_pbc && visualize_pbc.is_none() {
        arg_error(
            "--pbc-max-weight bounds PBC output — combine it with --to-pbc or --visualize-pbc",
        );
    }
    // Two writers, one stream: whichever won, the other's output would be
    // interleaved into it and neither would parse. Better to say so than to
    // emit a QASM file with a JSON object spliced through it.
    if json && output_path.as_deref() == Some(STREAM_PATH) {
        let file = if to_pbc { "out.pbc" } else { "out.qasm" };
        arg_error(format!(
            "--json and `-o -` both write to stdout — write the circuit to a \
             file (-o {file}) and keep --json on stdout, or drop --json"
        ));
    }
    // Flags that only configure gate-level optimization, which is off.
    let mut ignored = Vec::new();
    if gate_opts_off {
        if parallel.is_some() {
            ignored.push("--parallel/--no-parallel");
        }
        if superopt_qubits.is_some()
            || superopt_window_gates.is_some()
            || superopt_murm_entries.is_some()
            || !matches!(superopt_gates, SuperOptGates::Auto)
        {
            ignored.push("--superopt-*");
        }
    }

    Opts {
        action: Action::Optimize(Box::new(Run {
            input_path,
            output_path,
            to_pbc,
            pbc_opt,
            gate_opts_off,
            ignored,
            pbc_passes: [(pbc_passes.0, "ToPbc"), (pbc_passes.1, "PbcOpt")]
                .into_iter()
                .filter_map(|(listed, name)| listed.then_some(name))
                .collect(),
            visualize_pbc,
            pbc_expansion_budget: pbc_expansion_budget.unwrap_or(DEFAULT_PBC_EXPANSION_BUDGET),
            pbc_max_weight,
            parallel,
            // An absent `-O` flag means O3 too; the distinction only ever
            // mattered for the validation above, which has already run.
            options: Options {
                level: optimization_level.unwrap_or(Level::O3),
                passes,
                fixpoint,
                decompose_rotations,
                decompose_rz: false,
                decompose_cz,
                decompose_ccx,
                rz_epsilon,
                // Only what was asked for outright; [`Run::resolve_parallel`]
                // settles the rest once the circuit's size is known.
                parallel: parallel == Some(true),
                // Hidden (undocumented in `--help`) bounds overrides; `None`
                // means "use whichever preset the optimization level implies".
                superopt: SuperOptBounds {
                    qubits: superopt_qubits,
                    window_gates: superopt_window_gates,
                    murm_entries: superopt_murm_entries,
                },
                superopt_gates,
            },
        })),
        ui,
        json,
    }
}

/// Print the help text.
///
/// Help is requested output, so it goes to stdout and takes its styling from
/// stdout — `tzap --help | less` must not paint escapes into the pager.
/// Built as one string and written once through [`Ui::write_stdout`], so a
/// reader that closes early (`tzap --help | head`) ends the pipeline rather
/// than panicking inside a `println!`.
fn print_help(ui: &Ui) {
    let mut out = String::new();
    let bold = ui.out_sgr("\x1b[1m");
    let dim = ui.out_sgr("\x1b[2m");
    let heading = ui.out_sgr("\x1b[1;33m");
    let reset = ui.out_reset();
    out.push('\n');
    out.push_str(&format!(
        "  {bold}⚡\u{FE0F} tzap{reset}  —  fast quantum circuit optimizer  {dim}v{}{reset}\n",
        env!("CARGO_PKG_VERSION")
    ));
    out.push('\n');
    out.push_str(&format!("  {heading}USAGE{reset}\n"));
    out.push_str("    tzap <input.qasm> [output] [options]\n");
    out.push('\n');
    out.push_str(&format!("  {heading}ARGS{reset}\n"));
    out.push_str(&format!(
        "    {bold}<input.qasm>{reset}     Input OpenQASM 2.0 file, or - for stdin\n"
    ));
    out.push_str(&format!(
        "    {bold}[output]{reset}         Output file, or - for stdout (no output if omitted)\n"
    ));
    out.push('\n');
    out.push_str(&format!("  {heading}OPTIONS{reset}\n"));
    out.push_str(&format!(
        "    {bold}-o{reset} <file>        Write output to <file> (- for stdout)\n"
    ));
    out.push_str(&format!(
        "    {bold}--decompose-ccx{reset}  Decompose CCX and CCZ gates into Clifford+T\n"
    ));
    out.push_str(&format!(
        "    {bold}--decompose-cz{reset}   Decompose CZ gates into H+CX+H\n"
    ));
    out.push_str(&format!(
        "    {bold}--decompose-rotations{reset}   Decompose P, Rz/Rx/Ry, CP, CRz/CRx/CRy into Clifford+T (gridsynth)\n"
    ));
    out.push_str("                     --decompose-rz is a compatibility alias.\n");
    out.push_str(&format!("    {bold}--epsilon{reset} <eps>  Approximation epsilon for --decompose-rotations (default: 1e-10)\n"));
    out.push_str(&format!(
        "    {bold}--superopt-gates{reset} <basis>  MURM basis: auto (default), base, or a\n"
    ));
    out.push_str("                     comma-separated gate list (for example h,t,tdg,cx,cz)\n");
    out.push_str(&format!(
        "    {bold}--parallel{reset}       Optimize chunks of the circuit concurrently. On by\n"
    ));
    out.push_str(&format!(
        "                     default from {} gates up\n",
        crate::progress::fmt_num(AUTO_PARALLEL_GATES)
    ));
    out.push_str(&format!(
        "    {bold}--no-parallel{reset}    Keep the run sequential at any size\n"
    ));
    out.push_str(&format!("    {bold}--passes{reset} <list>  Run these passes in order, overriding the default pipeline\n"));
    out.push_str("                     (see PASSES). Excludes --decompose-* — list the\n");
    out.push_str("                     corresponding decomposition pass names directly.\n");
    out.push_str("                     --epsilon still configures DecomposeRotations.\n");
    out.push_str(&format!(
        "    {bold}--fixpoint{reset}       Repeat the pipeline until gate count stops decreasing\n"
    ));
    out.push_str(&format!(
        "    {bold}-O1{reset}              Fastest: phase folding + gate cancellation only\n"
    ));
    out.push_str(&format!(
        "    {bold}-O2{reset}              Adds a superoptimization pass to O1 (2 rounds)\n"
    ));
    out.push_str(&format!("    {bold}-O3{reset}              Like -O2, run to a fixpoint instead of 2 rounds (default)\n"));
    out.push_str(&format!(
        "    {bold}-Osuper{reset}          Like -O3, with a larger SuperOpt window/MURM (slower\n"
    ));
    out.push_str("                     first run; the MURM is cached to disk afterward)\n");
    out.push_str(&format!(
        "    {bold}--json{reset}           Write a machine-readable report of the run to stdout\n"
    ));
    out.push_str(&format!(
        "    {bold}-q, --quiet{reset}      Print nothing but errors (output is still written)\n"
    ));
    out.push_str(&format!(
        "    {bold}--cache-dir{reset} <d>  Keep the MURM cache under <d>\n"
    ));
    out.push_str(&format!(
        "    {bold}--cache-info{reset}     Report the cache location and MURMs in it\n"
    ));
    out.push_str(&format!(
        "    {bold}--clear-cache{reset}    Delete every cached MURM\n"
    ));
    out.push_str(&format!(
        "    {bold}-h, --help{reset}       Print this help message\n"
    ));
    out.push_str(&format!(
        "    {bold}-V, --version{reset}    Print the version\n"
    ));
    out.push('\n');
    out.push_str(&format!("  {heading}PBC OUTPUT{reset}\n"));
    out.push_str(&format!(
        "    {bold}--to-pbc{reset}         Convert to Pauli-based computation (PBC), then optimize\n"
    ));
    out.push_str("                     its rotations; write it with -o <file.pbc> or -o -.\n");
    out.push_str(
        "                     Gate-level optimization is off unless -O1/-O2/-O3/-Osuper\n",
    );
    out.push_str("                     or --passes is given. Same as ending --passes with\n");
    out.push_str("                     ToPbc,PbcOpt (ToPbc alone with --pbc-no-opt).\n");
    out.push_str(&format!(
        "    {bold}--pbc-no-opt{reset}     Skip the PBC rotation optimizer\n"
    ));
    out.push_str(&format!("    {bold}--pbc-max-weight{reset} <N>\n"));
    out.push_str("                     Bound the weight of π/8 rotations and measurements\n");
    out.push_str("                     (N >= 1); Cliffords that would widen them become π/4\n");
    out.push_str("                     and π/2 rotations of weight <= 2\n");
    out.push_str(&format!("    {bold}--visualize-pbc{reset} <file.svg>\n"));
    out.push_str("                     Draw the PBC circuit as SVG, in Litinski's style\n");
    out.push_str(&format!("    {bold}--pbc-expansion-budget{reset} <N>\n"));
    out.push_str("                     Work allowed for expanding Pauli axes when writing,\n");
    out.push_str(&format!(
        "                     drawing, or measuring PBC (default {})\n",
        crate::progress::fmt_num(DEFAULT_PBC_EXPANSION_BUDGET)
    ));
    out.push('\n');
    out.push_str(&format!("  {heading}PASSES{reset} (names for --passes)\n"));
    for (name, _pass, desc) in PassName::ALL {
        out.push_str(&format!("    {bold}{name:<19}{reset}  {desc}\n"));
    }
    for (name, desc) in PBC_PASSES {
        out.push_str(&format!("    {bold}{name:<19}{reset}  {desc}\n"));
    }
    out.push('\n');
    out.push_str(&format!("  {heading}OUTPUT{reset}\n"));
    out.push_str("    Progress and messages go to stderr; the circuit and --json go to\n");
    out.push_str("    stdout, so tzap composes in a pipeline:\n");
    out.push_str(&format!(
        "      {dim}tzap in.qasm -o - | other-tool{reset}\n"
    ));
    out.push_str(&format!(
        "      {dim}tzap in.qasm -O3 --json -q | jq .reduction_percent{reset}\n"
    ));
    out.push_str("    Color and live progress bars are used only when stderr is a terminal.\n");
    out.push('\n');
    out.push_str(&format!("  {heading}ENVIRONMENT{reset}\n"));
    out.push_str(&format!(
        "    {bold}TZAP_CACHE_DIR{reset}   Cache root (overridden by --cache-dir)\n"
    ));
    out.push_str(&format!(
        "    {bold}XDG_CACHE_HOME{reset}   Cache root falls back to $XDG_CACHE_HOME/tzap,\n"
    ));
    out.push_str("                     then $HOME/.cache/tzap, then the Windows\n");
    out.push_str("                     %LOCALAPPDATA% and %USERPROFILE% locations\n");
    out.push('\n');
    ui.write_stdout(&out);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_flag_values_may_be_written_with_an_equals_sign() {
        let args: Vec<String> = ["tzap", "--passes=CancelGates", "in.qasm", "--epsilon=1e-8"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(
            split_flag_values(&args),
            vec![
                "tzap",
                "--passes",
                "CancelGates",
                "in.qasm",
                "--epsilon",
                "1e-8"
            ]
        );
    }

    /// Only the known value-taking flags are split, so an `=` anywhere else —
    /// a file name, a pass list — survives intact.
    #[test]
    fn other_arguments_containing_equals_are_left_alone() {
        let args: Vec<String> = ["tzap", "a=b.qasm", "--parallel", "--unknown=x"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(split_flag_values(&args), args);
    }
}
