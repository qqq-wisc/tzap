//! Terminal etiquette and stream discipline: what tzap is allowed to put on
//! stderr when nobody is watching it through a terminal, what belongs on
//! stdout, and what `--quiet` silences.
//!
//! Every test here runs tzap as a child process, so *both* its streams are
//! pipes rather than terminals — which is precisely the case these tests
//! exist to pin down. There is no way to ask for styling anyway: a terminal
//! gets it and nothing else does. (The live rendering path is covered by unit
//! tests in `src/main.rs`, which can hand it a `Ui` that claims a terminal.)

#[path = "support/mod.rs"]
mod support;

use std::fs;

use support::{
    Json, Tzap, assert_plain, assert_valid_qasm, gate_lines, read, tzap as raw_tzap,
    without_timings,
};

const TEST_QASM: &str = "tests/fixtures/test.qasm";
const TWO_CCX_QASM: &str = "tests/fixtures/two_ccx.qasm";
const MOD5_4_QASM: &str = "benchmarks/feynman/mod5_4.qasm";

/// A circuit with an Rz gate, so the Rz progress row and `--decompose-rz`
/// have something to act on.
const RZ_QASM: &str = "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
h q[0];
rz(pi/5) q[0];
cx q[0],q[1];
t q[1];
";

/// Every optimization level, including the one whose bounds are only
/// affordable in a test with the hidden overrides applied.
const LEVELS: [&[&str]; 5] = [
    &[
        "--superopt-qubits",
        "3",
        "--superopt-window-gates",
        "4",
        "--superopt-murm-entries",
        "500",
    ],
    &["-O1"],
    &[
        "-O2",
        "--superopt-qubits",
        "3",
        "--superopt-window-gates",
        "4",
        "--superopt-murm-entries",
        "500",
    ],
    &[
        "-O3",
        "--superopt-qubits",
        "3",
        "--superopt-window-gates",
        "4",
        "--superopt-murm-entries",
        "500",
    ],
    &[
        "-Osuper",
        "--superopt-qubits",
        "2",
        "--superopt-window-gates",
        "4",
        "--superopt-murm-entries",
        "500",
    ],
];

/// SuperOpt bounds small enough to build and cache a real MURM in
/// milliseconds, for the test that needs to do it twice (see `LEVELS` above,
/// which spells out the same thing for `-Osuper`).
const TINY_SUPEROPT: [&str; 6] = [
    "--superopt-qubits",
    "2",
    "--superopt-window-gates",
    "4",
    "--superopt-murm-entries",
    "500",
];

/// Keep tests of stream behavior independent of the production MURM size.
fn tzap(args: &[&str]) -> support::Run {
    let mut full = args.to_vec();
    if !args.contains(&"--superopt-qubits") {
        full.extend_from_slice(&TINY_SUPEROPT);
    }
    raw_tzap(&full)
}

fn level_name(level: &[&str]) -> String {
    if level.is_empty() {
        "default".to_string()
    } else {
        level[0].to_string()
    }
}

// ---------------------------------------------------------------------------
// Color and cursor motion (clig.dev: "disable color if not in a terminal")
// ---------------------------------------------------------------------------

/// The headline rule: with stderr piped, nothing tzap prints may contain an
/// escape sequence or a carriage return — across every level, sequential and
/// parallel, with and without a decomposition.
#[test]
fn a_piped_run_emits_no_escapes_at_any_level() {
    for level in LEVELS {
        for parallel in [&[][..], &["--parallel"][..]] {
            for input in [TEST_QASM, TWO_CCX_QASM, MOD5_4_QASM] {
                let mut args = vec![input];
                args.extend_from_slice(level);
                args.extend_from_slice(parallel);
                let context = format!("{input} {} {parallel:?}", level_name(level));
                let run = tzap(&args).ok(&context);
                assert_plain(&run.stderr, &format!("{context} stderr"));
                assert_plain(&run.stdout, &format!("{context} stdout"));
            }
        }
    }
}

/// The same rule for the flags that change what gets printed rather than
/// what gets optimized.
#[test]
fn a_piped_run_emits_no_escapes_under_any_verbosity() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    fs::write(&input, RZ_QASM).unwrap();
    let input = input.to_str().unwrap();

    let variants: [&[&str]; 7] = [
        &[],
        &["--quiet"],
        &["-q"],
        &["--json"],
        &["--json", "-q"],
        &["--decompose-rz", "--epsilon", "1e-2"],
        &["--parallel", "--passes", "CancelGates"],
    ];
    for variant in variants {
        let mut args = vec![input];
        args.extend_from_slice(variant);
        let run = tzap(&args).ok(&format!("{variant:?}"));
        assert_plain(&run.stderr, &format!("{variant:?} stderr"));
        assert_plain(&run.stdout, &format!("{variant:?} stdout"));
    }
}

/// Styling has no override: the flags that used to force or suppress it are
/// gone, and are rejected as the unknown flags they now are.
#[test]
fn there_is_no_way_to_ask_for_color() {
    for args in [
        vec![TEST_QASM, "--color", "always"],
        vec![TEST_QASM, "--color=never"],
        vec![TEST_QASM, "--color"],
        vec![TEST_QASM, "--no-color"],
    ] {
        let run = tzap(&args).failed(&format!("{args:?}"));
        assert!(
            run.stderr.contains("unknown flag"),
            "{args:?}: got {}",
            run.stderr
        );
    }
}

/// ...and no environment variable changes it either. A run with each of the
/// usual color conventions set is byte-for-byte the run without them: whether
/// the stream is a terminal is the whole of the decision.
#[test]
fn the_environment_does_not_change_what_is_printed() {
    let plain = tzap(&[MOD5_4_QASM, "-O1"]).ok("no environment");
    for (name, value) in [
        ("NO_COLOR", "1"),
        ("CLICOLOR_FORCE", "1"),
        ("CLICOLOR", "1"),
        ("TERM", "dumb"),
        ("TERM", "xterm-256color"),
    ] {
        let run = Tzap::new(&[MOD5_4_QASM, "-O1"])
            .env(name, value)
            .run()
            .ok(&format!("{name}={value}"));
        assert_plain(&run.stderr, &format!("{name}={value}"));
        // Timings are the one thing two identical runs may legitimately
        // disagree on, and they are printed to the millisecond.
        assert_eq!(
            without_timings(&run.stderr),
            without_timings(&plain.stderr),
            "{name}={value} changed what was printed"
        );
    }
}

/// The normalization the test above leans on, pinned so a comparison can't
/// pass by erasing everything: timings go, and nothing else does.
#[test]
fn only_timings_are_normalized_away() {
    assert_eq!(
        without_timings("Parsed a.qasm (871 B) in 0.000s\n↓27.8% · 79 → 57 · 12.345s\n"),
        "Parsed a.qasm (871 B) in N.NNNs\n↓27.8% · 79 → 57 · N.NNNs\n"
    );
}

/// Help is requested output: it goes to stdout, and takes its color from
/// stdout — which is a pipe here, so it arrives plain.
#[test]
fn help_goes_to_stdout_and_is_plain_when_piped() {
    for flag in ["--help", "-h"] {
        let run = tzap(&[flag]).ok(flag);
        assert_plain(&run.stdout, flag);
        assert!(
            run.stderr.is_empty(),
            "{flag} wrote to stderr: {}",
            run.stderr
        );
        for expected in [
            "USAGE",
            "OPTIONS",
            "PASSES",
            "OUTPUT",
            "ENVIRONMENT",
            "--json",
            "--quiet",
            "--cache-dir",
            "--cache-info",
            "--clear-cache",
            "TZAP_CACHE_DIR",
            "XDG_CACHE_HOME",
        ] {
            assert!(
                run.stdout.contains(expected),
                "{flag} should document {expected}:\n{}",
                run.stdout
            );
        }
    }
}

/// The whole line is parsed before help prints, so a flag that changes how it
/// prints works on either side of the `--help` that consumes it.
#[test]
fn help_is_printed_whichever_side_its_flags_are_on() {
    let plain = tzap(&["--help"]).ok("--help").stdout;
    for args in [vec!["--help", "-q"], vec!["-q", "--help"]] {
        let run = tzap(&args).ok(&format!("{args:?}"));
        assert_eq!(
            run.stdout, plain,
            "{args:?}: help is requested output, so --quiet doesn't suppress it"
        );
    }
}

#[test]
fn version_goes_to_stdout_under_every_spelling() {
    for flag in ["--version", "-v", "-V"] {
        let run = tzap(&[flag]).ok(flag);
        assert_eq!(
            run.stdout.trim(),
            format!("tzap {}", env!("CARGO_PKG_VERSION"))
        );
        assert!(run.stderr.is_empty(), "{flag} wrote to stderr");
    }
}

// ---------------------------------------------------------------------------
// stdout as the circuit's destination (clig.dev: primary output to stdout)
// ---------------------------------------------------------------------------

/// `-o -` sends the circuit to stdout, byte for byte what the same run would
/// have written to a file — with the commentary still on stderr, which is
/// what makes the stream pipeable.
#[test]
fn dash_writes_the_circuit_to_stdout() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("out.qasm");

    let to_file = tzap(&[MOD5_4_QASM, "-o", file.to_str().unwrap(), "-O1"]).ok("-o file");
    let to_stdout = tzap(&[MOD5_4_QASM, "-o", "-", "-O1"]).ok("-o -");

    assert_eq!(to_stdout.stdout, read(&file));
    assert_valid_qasm(&to_stdout.stdout, "-o -");
    assert!(
        to_file.stdout.is_empty(),
        "writing to a file must leave stdout empty, got:\n{}",
        to_file.stdout
    );
    // The commentary is unchanged; only the destination moved.
    assert!(to_stdout.stderr.contains("Final result"));
    assert!(
        !to_stdout.stderr.contains("OPENQASM"),
        "the circuit must not be duplicated onto stderr:\n{}",
        to_stdout.stderr
    );
    // No "wrote <path>" line for a stream that has no path.
    assert!(
        !to_stdout.stderr.contains("wrote"),
        "got: {}",
        to_stdout.stderr
    );
    assert!(to_file.stderr.contains("wrote"));
}

/// The second positional argument takes `-` too, not just `-o`.
#[test]
fn a_positional_dash_writes_to_stdout() {
    let run = tzap(&[MOD5_4_QASM, "-", "-O1"]).ok("positional -");
    assert_valid_qasm(&run.stdout, "positional -");
}

/// `-` as the input reads the circuit from stdin, and says so rather than
/// naming a file that doesn't exist.
#[test]
fn dash_reads_the_circuit_from_stdin() {
    let run = Tzap::new(&["-", "-O1", "-o", "-"])
        .stdin(&fs::read_to_string(MOD5_4_QASM).unwrap())
        .run()
        .ok("stdin to stdout");

    let piped = assert_valid_qasm(&run.stdout, "stdin to stdout");
    let from_file = tzap(&[MOD5_4_QASM, "-o", "-", "-O1"]).ok("file to stdout");
    assert_eq!(run.stdout, from_file.stdout);
    assert!(!piped.is_empty());
    assert!(
        run.stderr.contains("<stdin>"),
        "a stdin run should name its input <stdin>:\n{}",
        run.stderr
    );
    // Size is reported from what actually arrived, not from a stat that
    // couldn't have happened.
    assert!(run.stderr.contains(" KB) in ") || run.stderr.contains(" B) in "));
}

/// Reading and writing streams works at every level and in parallel mode —
/// the full pipeline shape, `cat x.qasm | tzap - -o - | ...`.
#[test]
fn streaming_works_at_every_level() {
    let qasm = fs::read_to_string(TWO_CCX_QASM).unwrap();
    for level in LEVELS {
        for parallel in [&[][..], &["--parallel"][..]] {
            let mut args = vec!["-", "-o", "-"];
            args.extend_from_slice(level);
            args.extend_from_slice(parallel);
            let context = format!("{} {parallel:?} streamed", level_name(level));
            let run = Tzap::new(&args).stdin(&qasm).run().ok(&context);
            let gates = assert_valid_qasm(&run.stdout, &context);
            assert!(!gates.is_empty(), "{context}: no gates in the output");
        }
    }

    let run = Tzap::new(&["-", "-o", "-", "-O1", "--decompose-ccx"])
        .stdin(&qasm)
        .run()
        .ok("streamed --decompose-ccx");
    let gates = assert_valid_qasm(&run.stdout, "streamed --decompose-ccx");
    assert!(
        !gates.iter().any(|gate| gate.starts_with("ccx ")),
        "the opt-in decomposition must remove Toffolis"
    );
}

/// Malformed input on stdin fails the way a malformed file does: a clear
/// error, nothing on stdout, non-zero exit.
#[test]
fn malformed_stdin_fails_cleanly() {
    let run = Tzap::new(&["-", "-o", "-"])
        .stdin("this is not qasm at all\n")
        .run()
        .failed("malformed stdin");
    assert!(
        run.stderr.contains("<stdin>"),
        "the error should name the input as <stdin>:\n{}",
        run.stderr
    );
    // The message must not begin with a stray blank line left over from an
    // in-progress "Parsing..." line that was never printed.
    assert!(
        !run.stderr.starts_with('\n'),
        "leading blank line in:\n{:?}",
        run.stderr
    );
}

/// Empty input is an empty circuit, not an error — and stdin agrees with a
/// file, which is the property that matters: `-` must be the same code path,
/// not a second one with its own behavior.
#[test]
fn empty_stdin_matches_an_empty_file() {
    let dir = tempfile::tempdir().unwrap();
    let empty = dir.path().join("empty.qasm");
    fs::write(&empty, "").unwrap();

    let from_file = tzap(&[empty.to_str().unwrap(), "-o", "-", "-O1"]).ok("empty file");
    let from_stdin = Tzap::new(&["-", "-o", "-", "-O1"])
        .stdin("")
        .run()
        .ok("empty stdin");
    assert_eq!(from_stdin.stdout, from_file.stdout);
    assert!(from_stdin.stdout.starts_with("OPENQASM 2.0;"));
    assert!(from_stdin.stderr.contains("0 qubits · 0 gates"));
}

/// With no output argument at all, stdout stays completely empty — the run is
/// a measurement, and its numbers are commentary on stderr.
#[test]
fn no_output_argument_leaves_stdout_empty() {
    for level in LEVELS {
        let mut args = vec![TEST_QASM];
        args.extend_from_slice(level);
        let run = tzap(&args).ok(&level_name(level));
        assert!(
            run.stdout.is_empty(),
            "{}: expected empty stdout, got:\n{}",
            level_name(level),
            run.stdout
        );
    }
}

// ---------------------------------------------------------------------------
// --quiet
// ---------------------------------------------------------------------------

/// `--quiet` silences the commentary entirely — and only the commentary. The
/// circuit is still written, wherever it was asked to go.
#[test]
fn quiet_silences_stderr_but_not_the_output() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("out.qasm");

    for flag in ["-q", "--quiet"] {
        let run = tzap(&[MOD5_4_QASM, "-o", file.to_str().unwrap(), "-O1", flag]).ok(flag);
        assert_eq!(run.stderr, "", "{flag} should print nothing to stderr");
        assert!(run.stdout.is_empty());
        assert_valid_qasm(&read(&file), flag);

        let run = tzap(&[MOD5_4_QASM, "-o", "-", "-O1", flag]).ok(flag);
        assert_eq!(run.stderr, "", "{flag} should print nothing to stderr");
        assert_valid_qasm(&run.stdout, flag);
    }
}

/// Quiet at every level and in parallel mode, including the ones that load a
/// MURM and converge over several rounds — each of which has its
/// own message to suppress.
#[test]
fn quiet_is_silent_at_every_level() {
    for level in LEVELS {
        for parallel in [&[][..], &["--parallel"][..]] {
            let mut args = vec![TEST_QASM, "-q"];
            args.extend_from_slice(level);
            args.extend_from_slice(parallel);
            let context = format!("{} {parallel:?} -q", level_name(level));
            let run = tzap(&args).ok(&context);
            assert_eq!(run.stderr, "", "{context} should be silent");
        }
    }
}

/// Quiet is not a way to lose an error: a failing run still explains itself
/// and still exits non-zero.
#[test]
fn quiet_still_reports_errors() {
    let run = tzap(&["definitely-missing.qasm", "-q"]).failed("-q on a missing file");
    assert!(run.stderr.contains("Error reading"), "got: {}", run.stderr);

    let dir = tempfile::tempdir().unwrap();
    let bad = dir.path().join("bad.qasm");
    fs::write(&bad, "not qasm").unwrap();
    tzap(&[bad.to_str().unwrap(), "--quiet"]).failed("-q on a bad circuit");

    tzap(&["--quiet", "--passes", "NoSuchPass"]).failed("-q with a bad pass name");
}

/// A piped run stays quiet *during* optimization and still reports its
/// result: the progress boxes are for a terminal, and nothing replaces them
/// elsewhere.
#[test]
fn a_piped_run_reports_its_result_without_per_round_noise() {
    let run = tzap(&[MOD5_4_QASM, "-O3"]).ok("-O3");
    assert!(run.stderr.contains("Final result"), "got:\n{}", run.stderr);
    assert!(
        run.stderr.contains("Converged after"),
        "got:\n{}",
        run.stderr
    );
    for absent in ["Iteration 1", "% reduction so far", "Parallel optimization"] {
        assert!(
            !run.stderr.contains(absent),
            "{absent:?} is a live-terminal frame and must not be appended to a pipe:\n{}",
            run.stderr
        );
    }

    let parallel = tzap(&[MOD5_4_QASM, "-O1", "--parallel"]).ok("-O1 --parallel");
    assert!(
        !parallel.stderr.contains("Parallel optimization"),
        "got:\n{}",
        parallel.stderr
    );
    assert!(parallel.stderr.contains("Final result"));
}

/// There is no verbosity above the default, so `--verbose` is an unknown flag
/// rather than a silently accepted no-op.
#[test]
fn there_is_no_verbose_flag() {
    for flag in ["--verbose", "--verbosity", "-vv"] {
        let run = tzap(&[TEST_QASM, flag]).failed(flag);
        assert!(
            run.stderr.contains("unknown flag"),
            "{flag}: got {}",
            run.stderr
        );
    }
}

/// The parse line collapses to a single completed line when there's nothing
/// to overwrite, rather than printing both halves.
#[test]
fn the_parse_line_is_printed_once_when_piped() {
    let run = tzap(&[TEST_QASM, "-O1"]).ok("parse line");
    assert_eq!(
        run.stderr.matches("test.qasm").count(),
        1,
        "expected exactly one line naming the input:\n{}",
        run.stderr
    );
    assert!(run.stderr.contains("Parsed "), "got:\n{}", run.stderr);
    assert!(
        !run.stderr.contains("Parsing "),
        "the in-progress half has nothing to overwrite and must be skipped:\n{}",
        run.stderr
    );
}

/// The same for the MURM-load line, which uses the same overwrite-in-place
/// mechanism — on both of its paths: a cold run announces the build it is
/// about to do, a warm one reports the load, and neither leaves the
/// in-progress half on screen.
///
/// The cache is pinned to a directory this test owns, with bounds small
/// enough to build in milliseconds. Reading the default location instead made
/// the result depend on whether this machine happened to have that MURM
/// cached already: warm elsewhere, and always cold on Windows, where none of
/// the cache locations resolved at all.
#[test]
fn the_murm_line_is_printed_once_when_piped() {
    let cache = tempfile::tempdir().unwrap();
    let mut args = vec![
        TEST_QASM,
        "-O2",
        "--cache-dir",
        cache.path().to_str().unwrap(),
    ];
    args.extend_from_slice(&TINY_SUPEROPT);

    let cold = tzap(&args).ok("cold MURM line");
    assert!(
        cold.stderr.contains("Building MURM"),
        "a cold run must say the MURM is being built:\n{}",
        cold.stderr
    );

    let warm = tzap(&args).ok("warm MURM line");
    assert_eq!(
        warm.stderr.matches("MURM").count(),
        1,
        "expected exactly one MURM-load line:\n{}",
        warm.stderr
    );
    assert!(warm.stderr.contains("Loaded MURM"), "got:\n{}", warm.stderr);

    for (path, run) in [("cold", &cold), ("warm", &warm)] {
        assert!(
            !run.stderr.contains("Loading MURM"),
            "{path}: the in-progress half has nothing to overwrite and must \
             be skipped:\n{}",
            run.stderr
        );
    }
}

// ---------------------------------------------------------------------------
// Interactions between the two new destinations
// ---------------------------------------------------------------------------

/// `--json` and `-o -` both want stdout; tzap says so instead of interleaving
/// a JSON object through a QASM file.
#[test]
fn json_and_stdout_output_cannot_be_combined() {
    for args in [
        vec![TEST_QASM, "--json", "-o", "-"],
        vec![TEST_QASM, "-o", "-", "--json"],
        vec![TEST_QASM, "-", "--json"],
    ] {
        let run = tzap(&args).failed(&format!("{args:?}"));
        assert!(
            run.stderr.contains("--json") && run.stderr.contains("stdout"),
            "the error should explain the conflict:\n{}",
            run.stderr
        );
    }
}

/// `--json` alongside a file destination is the useful combination, and both
/// arrive intact.
#[test]
fn json_and_a_file_destination_coexist() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("out.qasm");
    let run = tzap(&[MOD5_4_QASM, "-o", file.to_str().unwrap(), "-O1", "--json"]).ok("--json -o");

    let report = Json::parse(&run.stdout);
    let written = assert_valid_qasm(&read(&file), "--json -o");
    assert_eq!(
        report.at("metrics/output/gates").as_usize(),
        written.len(),
        "the report's gate count must match the circuit it wrote"
    );
    assert!(run.stderr.contains("Final result"));
}

/// A quiet JSON run is the scriptable one: stdout is exactly the report,
/// stderr exactly nothing.
#[test]
fn quiet_json_puts_the_report_alone_on_stdout() {
    let run = tzap(&[MOD5_4_QASM, "-O1", "--json", "-q"]).ok("--json -q");
    assert_eq!(run.stderr, "");
    let report = Json::parse(&run.stdout);
    assert_eq!(report.get("tzap").as_str(), env!("CARGO_PKG_VERSION"));
}

/// Output is written before the report, and neither corrupts the other when
/// both streams are captured together.
#[test]
fn the_two_streams_stay_separable() {
    let run = tzap(&[MOD5_4_QASM, "-o", "-", "-O1"]).ok("-o -");
    let combined = run.both();
    assert!(combined.contains("OPENQASM 2.0;"));
    // stdout on its own is parseable QASM with no commentary spliced in.
    for line in gate_lines(&run.stdout) {
        assert!(
            !line.contains("Final result") && !line.contains("tzap"),
            "commentary leaked into the circuit: {line}"
        );
    }
}

// ---------------------------------------------------------------------------
// Pipeline manners
// ---------------------------------------------------------------------------

/// A reader that closes early (`tzap in.qasm -o - | head -1`) ends the
/// pipeline the way every Unix tool does: quietly, with a zero exit — not
/// with a broken-pipe error, and above all not with a Rust panic out of a
/// failed `println!`.
///
/// Needs an output larger than a pipe buffer for the write to fail at all,
/// hence the largest benchmark rather than a fixture.
#[test]
fn a_closed_reader_ends_the_pipeline_quietly() {
    use std::io::Read;
    use std::process::{Command, Stdio};

    const BIG_QASM: &str = "benchmarks/qft/qft_q020_d32421.qasm";
    if !std::path::Path::new(BIG_QASM).exists() {
        return;
    }

    for args in [
        vec![BIG_QASM, "-o", "-", "--passes", "CancelGates", "-q"],
        vec![BIG_QASM, "--json", "-q", "--passes", "CancelGates"],
        vec!["--help"],
    ] {
        let mut child = Command::new(env!("CARGO_BIN_EXE_tzap"))
            .args(&args)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("failed to spawn tzap");

        // Read one small chunk, then drop the pipe while tzap is still
        // writing.
        let mut stdout = child.stdout.take().expect("piped stdout");
        let mut first = [0u8; 16];
        let _ = stdout.read(&mut first);
        drop(stdout);

        let mut stderr = String::new();
        child
            .stderr
            .take()
            .expect("piped stderr")
            .read_to_string(&mut stderr)
            .expect("failed to read stderr");
        let status = child.wait().expect("failed to wait for tzap");

        assert!(
            !stderr.to_lowercase().contains("panic"),
            "{args:?}: panicked on a closed reader:\n{stderr}"
        );
        assert!(
            !stderr.to_lowercase().contains("broken pipe"),
            "{args:?}: a closed reader is not an error to report:\n{stderr}"
        );
        assert!(
            status.success(),
            "{args:?}: expected a clean exit, got {:?}\n{stderr}",
            status.code()
        );
    }
}

/// Every value-taking long flag accepts the `--flag=value` spelling, which is
/// what most people reach for first.
#[test]
fn value_flags_accept_the_equals_spelling() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    fs::write(&input, RZ_QASM).unwrap();
    let input = input.to_str().unwrap();
    let out = dir.path().join("out.qasm");
    let out_flag = format!("-o={}", out.to_str().unwrap());

    // The long flags, in both spellings, must agree.
    let pairs: [(Vec<&str>, Vec<&str>); 4] = [
        (
            vec!["--epsilon=1e-2", "--decompose-rz"],
            vec!["--epsilon", "1e-2", "--decompose-rz"],
        ),
        (
            vec!["--passes=CancelGates,PhaseFoldRand"],
            vec!["--passes", "CancelGates,PhaseFoldRand"],
        ),
        (
            vec!["--superopt-qubits=2", "--passes=SuperOpt"],
            vec!["--superopt-qubits", "2", "--passes", "SuperOpt"],
        ),
        (
            vec!["--superopt-murm-entries=400", "--passes=SuperOpt"],
            vec!["--superopt-murm-entries", "400", "--passes", "SuperOpt"],
        ),
    ];
    for (equals, spaced) in pairs {
        let mut with_equals = vec![input, "-o", "-", "-q"];
        with_equals.extend_from_slice(&equals);
        let mut with_space = vec![input, "-o", "-", "-q"];
        with_space.extend_from_slice(&spaced);

        let a = tzap(&with_equals).ok(&format!("{equals:?}"));
        let b = tzap(&with_space).ok(&format!("{spaced:?}"));
        assert_eq!(a.stdout, b.stdout, "{equals:?} differs from {spaced:?}");
        assert!(!a.stdout.is_empty());
    }

    // The splitting is deliberately limited to the known long flags, so the
    // short `-o=file` is rejected as the unknown flag it is rather than
    // silently writing to a file named "=file".
    let run = tzap(&[input, &out_flag, "-q"]).failed("-o=file");
    assert!(run.stderr.contains("unknown flag"), "got: {}", run.stderr);
    assert!(!out.exists(), "nothing should have been written");
}

/// A file name containing an `=` is a positional argument, not a flag with a
/// value, and survives untouched.
#[test]
fn a_filename_containing_an_equals_sign_is_left_alone() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("a=b.qasm");
    fs::write(&input, RZ_QASM).unwrap();

    let run = tzap(&[input.to_str().unwrap(), "-o", "-", "-q"]).ok("a=b.qasm");
    assert_valid_qasm(&run.stdout, "a=b.qasm");
}
