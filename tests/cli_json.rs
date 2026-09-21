//! The `--json` report: that it is real JSON, that its numbers agree with
//! the circuit tzap actually wrote and with the human banner, and that every
//! documented key is present with the right type for every shape of run.
//!
//! The parser in `support` is written independently of `src/json.rs`, so
//! these tests check the output rather than the writer's idea of it.

#[path = "support/mod.rs"]
mod support;

use std::fs;

use support::{Json, Tzap, assert_valid_qasm, read, tzap};

const TEST_QASM: &str = "tests/fixtures/test.qasm";
const TWO_CCX_QASM: &str = "tests/fixtures/two_ccx.qasm";
const MOD5_4_QASM: &str = "benchmarks/feynman/mod5_4.qasm";

const RZ_QASM: &str = "\
OPENQASM 2.0;
include \"qelib1.inc\";
qreg q[2];
h q[0];
rz(pi/5) q[0];
cx q[0],q[1];
t q[1];
";

/// Small `-Osuper` bounds, so the level can be exercised without building its
/// real five-million-entry MURM.
const SMALL_SUPER: [&str; 6] = [
    "--superopt-qubits",
    "2",
    "--superopt-window-gates",
    "4",
    "--superopt-murm-entries",
    "500",
];

const TINY_MURM: [&str; 6] = [
    "--superopt-qubits",
    "3",
    "--superopt-window-gates",
    "4",
    "--superopt-murm-entries",
    "500",
];

/// Run tzap with `--json -q` and parse the report. Quiet, so any stray
/// commentary landing on stdout would break the parse loudly.
fn json_report(args: &[&str]) -> Json {
    let mut full = vec!["--json", "-q"];
    full.extend_from_slice(args);
    if !args.contains(&"--superopt-qubits") {
        full.extend_from_slice(&TINY_MURM);
    }
    let run = tzap(&full).ok(&format!("{args:?}"));
    assert_eq!(run.stderr, "", "--json -q should be silent on stderr");
    Json::parse(&run.stdout)
}

/// Assert every key the report documents is present, and typed as promised.
fn assert_schema(report: &Json, context: &str) {
    assert_eq!(
        report.keys(),
        vec![
            "tzap",
            "input",
            "output",
            "options",
            "metrics",
            "reduction_percent",
            "stages",
            "passes",
            "table",
            "murms",
            "fixpoint",
            "fixpoints",
            "cache_dir",
            "seconds",
        ],
        "{context}: unexpected top-level shape"
    );

    report.get("tzap").as_str();
    report.at("input/stdin").as_bool();
    report.at("input/qubits").as_usize();
    report.at("input/gate_set").strings();
    report.at("input/parse_seconds").as_f64();
    report.at("output/stdout").as_bool();
    report.at("output/gate_set").strings();
    report.at("seconds").as_f64();

    let options = report.get("options");
    assert_eq!(
        options.keys(),
        vec![
            "level",
            "passes",
            "fixpoint",
            "decompose_rz",
            "decompose_cz",
            "decompose_ccx",
            "rz_epsilon",
            "parallel",
            "superopt",
        ],
        "{context}: unexpected options shape"
    );
    options.get("level").as_str();
    options.get("fixpoint").as_bool();
    options.get("decompose_rz").as_bool();
    options.get("decompose_cz").as_bool();
    options.get("decompose_ccx").as_bool();
    options.get("parallel").as_bool();
    assert!(
        options.get("rz_epsilon").as_f64() > 0.0,
        "{context}: epsilon must be positive"
    );
    for key in ["qubits", "window_gates", "murm_entries"] {
        assert!(
            options.at(&format!("superopt/{key}")).as_usize() > 0,
            "{context}: superopt/{key} must be a positive integer"
        );
    }
    options.at("superopt/gates").as_str();

    for stage in ["input", "baseline", "output"] {
        let metrics = report.get("metrics").get(stage);
        assert_eq!(
            metrics.keys(),
            vec!["gates", "two_qubit", "t", "rz", "depth"],
            "{context}: unexpected metrics shape for {stage}"
        );
        for key in metrics.keys() {
            metrics.get(key).as_usize();
        }
    }
    for key in report.get("reduction_percent").keys() {
        report.get("reduction_percent").get(key).as_f64();
    }
    report.get("stages").strings();
    for pass in report.get("passes").arr() {
        assert_eq!(
            pass.keys(),
            vec!["stage", "name", "input_gates", "output_gates", "seconds"],
            "{context}: unexpected pass shape"
        );
        pass.get("stage").as_usize();
        pass.get("name").as_str();
        pass.get("input_gates").as_usize();
        pass.get("output_gates").as_usize();
        pass.get("seconds").as_f64();
    }
    for murm in report.get("murms").arr() {
        assert_eq!(murm.keys(), vec!["stage", "cached", "basis", "seconds"]);
        murm.get("stage").as_usize();
        murm.get("cached").as_bool();
        murm.get("basis").strings();
        murm.get("seconds").as_f64();
    }
    let table = report.get("table");
    if !table.is_null() {
        assert_eq!(table.keys(), vec!["cached", "seconds"]);
        table.get("cached").as_bool();
        table.get("seconds").as_f64();
    }
    for fixpoint in report.get("fixpoints").arr() {
        assert_eq!(fixpoint.keys(), vec!["stage", "rounds", "converged"]);
        fixpoint.get("stage").as_usize();
        fixpoint.get("rounds").as_usize();
        fixpoint.get("converged").as_bool();
    }
    let fixpoint = report.get("fixpoint");
    if !fixpoint.is_null() {
        assert_eq!(fixpoint.keys(), vec!["rounds", "converged"]);
        fixpoint.get("rounds").as_usize();
        fixpoint.get("converged").as_bool();
    }
}

/// Every shape of run produces a report with the same schema.
#[test]
fn the_schema_holds_across_every_kind_of_run() {
    let dir = tempfile::tempdir().unwrap();
    let rz = dir.path().join("rz.qasm");
    fs::write(&rz, RZ_QASM).unwrap();
    let rz = rz.to_str().unwrap();
    let out = dir.path().join("out.qasm");
    let out = out.to_str().unwrap();

    let variants: [Vec<&str>; 12] = [
        vec![TEST_QASM],
        vec![TEST_QASM, "-O1"],
        vec![TEST_QASM, "-O2"],
        vec![TEST_QASM, "-O3"],
        {
            let mut args = vec![TEST_QASM, "-Osuper"];
            args.extend_from_slice(&SMALL_SUPER);
            args
        },
        vec![TEST_QASM, "--parallel"],
        vec![TEST_QASM, "-O1", "--parallel"],
        vec![TEST_QASM, "--passes", "CancelGates"],
        vec![
            TEST_QASM,
            "--passes",
            "CancelGates,PhaseFoldRand",
            "--fixpoint",
        ],
        vec![rz, "--decompose-rz", "--epsilon", "1e-2"],
        vec![rz, "--decompose-cz", "-O1"],
        vec![MOD5_4_QASM, "-O1", "-o", out],
    ];
    for variant in variants {
        assert_schema(&json_report(&variant), &format!("{variant:?}"));
    }
}

#[test]
fn stages_and_superopt_gate_mode_are_reported() {
    let report = json_report(&[
        TWO_CCX_QASM,
        "-O2",
        "--decompose-ccx",
        "--superopt-gates",
        "h,h,t,cx,cz",
    ]);
    assert_eq!(
        report.at("options/superopt/gates").as_str(),
        "h,t,cx,cz",
        "duplicate explicit gates should be canonicalized"
    );
    assert_eq!(
        report.get("stages").strings(),
        vec![
            "Optimizing input circuit",
            "Decomposing CCX/CCZ",
            "Optimizing decomposed circuit",
        ]
    );
    let fixpoints = report.get("fixpoints").arr();
    assert_eq!(fixpoints.len(), 2);
    assert_eq!(fixpoints[0].get("stage").as_usize(), 0);
    assert_eq!(fixpoints[1].get("stage").as_usize(), 2);
    let murms = report.get("murms").arr();
    assert_eq!(murms.len(), 2);
    assert_eq!(murms[0].get("stage").as_usize(), 0);
    assert_eq!(murms[1].get("stage").as_usize(), 2);
    let table = report.get("table");
    assert_eq!(
        table.get("cached").as_bool(),
        murms[1].get("cached").as_bool()
    );
    assert_eq!(
        table.get("seconds").as_f64(),
        murms[1].get("seconds").as_f64()
    );
    assert_eq!(report.get("passes").arr()[0].get("stage").as_usize(), 1);
}

/// The report is exactly one JSON document, newline-terminated, with nothing
/// before or after it.
#[test]
fn the_report_is_one_complete_document() {
    let run = tzap(&[TEST_QASM, "-O1", "--json", "-q"]).ok("--json -q");
    assert!(run.stdout.starts_with('{'), "got: {}", run.stdout);
    assert!(run.stdout.ends_with("}\n"), "got: {:?}", run.stdout);
    assert_eq!(run.stdout.matches("\"tzap\":").count(), 1);
    // The parser rejects trailing content, so this is the real check.
    Json::parse(&run.stdout);
}

/// The reported counts are the circuit's, not an approximation of it: gate
/// count matches the file tzap wrote, gate by gate.
#[test]
fn the_metrics_match_the_circuit_that_was_written() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("out.qasm");
    let path = file.to_str().unwrap();

    for level in [&["-O1"][..], &["-O3"][..], &["-O1", "--parallel"][..]] {
        let mut args = vec![MOD5_4_QASM, "-o", path];
        args.extend_from_slice(level);
        let report = json_report(&args);
        let written = assert_valid_qasm(&read(&file), &format!("{level:?}"));

        assert_eq!(
            report.at("metrics/output/gates").as_usize(),
            written.len(),
            "{level:?}: gate count disagrees with the written circuit"
        );
        let t_gates = written
            .iter()
            .filter(|gate| gate.starts_with("t ") || gate.starts_with("tdg "))
            .count();
        assert_eq!(report.at("metrics/output/t").as_usize(), t_gates);
        let two_qubit = written
            .iter()
            .filter(|gate| gate.starts_with("cx ") || gate.starts_with("cz "))
            .count();
        assert_eq!(report.at("metrics/output/two_qubit").as_usize(), two_qubit);
        let rz = written
            .iter()
            .filter(|gate| gate.starts_with("rz("))
            .count();
        assert_eq!(report.at("metrics/output/rz").as_usize(), rz);
    }
}

/// `reduction_percent` is measured against `baseline` and is the same number
/// the human banner leads with when both invocations use the same MURM bounds.
#[test]
fn the_reduction_agrees_with_the_human_banner() {
    let report = json_report(&[MOD5_4_QASM]);
    let baseline = report.at("metrics/baseline/gates").as_usize() as f64;
    let output = report.at("metrics/output/gates").as_usize() as f64;
    let reported = report.at("reduction_percent/gates").as_f64();
    assert!(
        (reported - (baseline - output) / baseline * 100.0).abs() < 1e-9,
        "reduction_percent/gates doesn't match its own metrics"
    );

    let mut human_args = vec![MOD5_4_QASM];
    human_args.extend_from_slice(&TINY_MURM);
    let human = tzap(&human_args).ok("human banner");
    assert!(
        human
            .stderr
            .contains(&format!("Final result · {reported:.1}% fewer gates")),
        "the banner should lead with the same figure ({reported:.1}%):\n{}",
        human.stderr
    );
}

/// A pass that grows the circuit reports a negative reduction rather than
/// clamping — `--passes DecomposeRz` expands one rotation into a long
/// Clifford+T sequence.
#[test]
fn a_circuit_that_grew_reports_a_negative_reduction() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    fs::write(&input, RZ_QASM).unwrap();

    let report = json_report(&[
        input.to_str().unwrap(),
        "--passes",
        "DecomposeRz",
        "--epsilon",
        "1e-2",
    ]);
    assert!(
        report.at("reduction_percent/gates").as_f64() < 0.0,
        "expected a negative reduction, got {}",
        report.at("reduction_percent/gates").as_f64()
    );
    assert!(
        report.at("metrics/output/gates").as_usize()
            > report.at("metrics/baseline/gates").as_usize()
    );
    assert_eq!(report.at("metrics/output/rz").as_usize(), 0);
}

/// The baseline remains the original input now that decomposition is an
/// opt-in middle stage between two optimization stages.
#[test]
fn input_and_baseline_remain_the_original_circuit() {
    let report = json_report(&[MOD5_4_QASM, "-O1"]);
    assert_eq!(
        report.at("metrics/input/gates").as_usize(),
        report.at("metrics/baseline/gates").as_usize(),
        "mod5_4 has no Toffolis to decompose"
    );

    let report = json_report(&[TWO_CCX_QASM, "-O1", "--decompose-ccx"]);
    assert_eq!(
        report.at("metrics/baseline/gates").as_usize(),
        report.at("metrics/input/gates").as_usize()
    );
    assert_eq!(report.at("metrics/input/gates").as_usize(), 3);
}

/// Whole-circuit passes are recorded with both gate counts and a timing.
#[test]
fn the_pass_list_records_requested_decompositions() {
    let report = json_report(&[TWO_CCX_QASM, "-O1", "--decompose-ccx"]);
    let passes = report.get("passes").arr();
    assert!(!passes.is_empty(), "expected the Toffoli decomposition");
    let toffoli = &passes[0];
    assert!(
        toffoli.get("name").as_str().contains("Toffoli"),
        "got: {}",
        toffoli.get("name").as_str()
    );
    assert_eq!(toffoli.get("input_gates").as_usize(), 3);
    assert!(toffoli.get("output_gates").as_usize() > 3);
    assert!(toffoli.get("seconds").as_f64() >= 0.0);

    // Nothing to decompose, nothing recorded.
    assert!(
        json_report(&[MOD5_4_QASM, "-O1"])
            .get("passes")
            .arr()
            .is_empty()
    );
}

/// `fixpoints` reports every stage's rounds and whether it converged — including the
/// `-O2` case, which stops on its round cap rather than at a true fixpoint.
#[test]
fn fixpoint_reports_rounds_and_convergence() {
    let o3 = json_report(&[MOD5_4_QASM, "-O3"]);
    let o3_fixpoints = o3.get("fixpoints").arr();
    assert_eq!(o3_fixpoints.len(), 1);
    assert!(
        o3_fixpoints[0].get("converged").as_bool(),
        "-O3 runs to a fixpoint"
    );
    assert!(o3_fixpoints[0].get("rounds").as_usize() >= 2);
    assert_eq!(
        o3.get("fixpoint").get("rounds").as_usize(),
        o3_fixpoints[0].get("rounds").as_usize()
    );
    assert_eq!(
        o3.get("fixpoint").get("converged").as_bool(),
        o3_fixpoints[0].get("converged").as_bool()
    );

    let o2 = json_report(&[MOD5_4_QASM, "-O2"]);
    assert!(o2.get("fixpoints").arr()[0].get("rounds").as_usize() <= 2);

    // A single-shot pipeline has no fixpoint to report.
    assert!(
        json_report(&[MOD5_4_QASM, "--passes", "CancelGates"])
            .get("fixpoints")
            .arr()
            .is_empty()
    );

    // ...and --fixpoint gives it one.
    let looped = json_report(&[MOD5_4_QASM, "--passes", "CancelGates", "--fixpoint"]);
    assert!(looped.get("fixpoints").arr()[0].get("rounds").as_usize() >= 1);
    assert!(looped.get("fixpoints").arr()[0].get("converged").as_bool());
}

#[test]
fn parallel_fixpoint_reports_one_maximum_round_record() {
    let report = json_report(&[MOD5_4_QASM, "-O2", "--parallel"]);
    let fixpoints = report.get("fixpoints").arr();
    assert_eq!(fixpoints.len(), 1);
    let rounds = fixpoints[0].get("rounds").as_usize();
    assert!((1..=2).contains(&rounds));
    assert_eq!(report.get("fixpoint").get("rounds").as_usize(), rounds);
}

/// `murms` records each MURM load for levels that use SuperOpt.
#[test]
fn the_murm_records_appear_only_for_levels_that_load_one() {
    assert!(
        json_report(&[TEST_QASM, "-O1"])
            .get("murms")
            .arr()
            .is_empty(),
        "-O1 doesn't use SuperOpt"
    );

    for level in ["-O2", "-O3"] {
        let report = json_report(&[TEST_QASM, level]);
        let murms = report.get("murms").arr();
        assert_eq!(murms.len(), 1);
        murms[0].get("cached").as_bool();
        murms[0].get("basis").strings();
        assert!(murms[0].get("seconds").as_f64() >= 0.0, "{level}");
    }
}

/// The echoed options are the effective ones: defaults filled in, presets
/// resolved, overrides applied.
#[test]
fn the_options_echo_what_the_run_actually_used() {
    let run = tzap(&["--json", "-q", TEST_QASM, "--passes", "CancelGates"])
        .ok("production default options");
    let default = Json::parse(&run.stdout);
    assert_eq!(default.at("options/level").as_str(), "O3");
    assert_eq!(default.at("options/superopt/qubits").as_usize(), 3);
    assert_eq!(default.at("options/superopt/window_gates").as_usize(), 25);
    assert_eq!(
        default.at("options/superopt/murm_entries").as_usize(),
        200_000
    );
    assert_eq!(default.at("options/rz_epsilon").as_f64(), 1e-10);
    assert!(!default.at("options/parallel").as_bool());

    for level in ["O1", "O2", "O3"] {
        let report = json_report(&[TEST_QASM, &format!("-{level}")]);
        assert_eq!(report.at("options/level").as_str(), level);
    }

    // -Osuper's preset shows through even though no bound was named...
    let mut args = vec![TEST_QASM, "-Osuper"];
    args.extend_from_slice(&SMALL_SUPER);
    let overridden = json_report(&args);
    assert_eq!(overridden.at("options/level").as_str(), "Osuper");
    assert_eq!(overridden.at("options/superopt/qubits").as_usize(), 2);
    assert_eq!(overridden.at("options/superopt/window_gates").as_usize(), 4);
    assert_eq!(
        overridden.at("options/superopt/murm_entries").as_usize(),
        500
    );

    // ...and a partial override leaves the rest at the preset.
    let partial = json_report(&[TEST_QASM, "-O1", "--superopt-qubits", "2"]);
    assert_eq!(partial.at("options/superopt/qubits").as_usize(), 2);
    assert_eq!(partial.at("options/superopt/window_gates").as_usize(), 25);
}

/// A `--passes` pipeline is echoed by name, in order, in the spelling the
/// flag itself accepts.
#[test]
fn the_pass_pipeline_round_trips_through_the_report() {
    let report = json_report(&[
        MOD5_4_QASM,
        "--passes",
        "CancelGates,PhaseFoldRand,CancelGates",
    ]);
    assert_eq!(
        report.at("options/passes").strings(),
        vec!["CancelGates", "PhaseFoldRand", "CancelGates"]
    );
    assert_eq!(report.at("options/level").as_str(), "O3");

    // Feeding the reported names back in is accepted.
    let names = report.at("options/passes").strings().join(",");
    tzap(&[MOD5_4_QASM, "--passes", &names, "-q"]).ok("round-tripped pass list");
}

#[test]
fn the_decomposition_and_epsilon_flags_are_echoed() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rz.qasm");
    fs::write(&input, RZ_QASM).unwrap();
    let input = input.to_str().unwrap();

    let report = json_report(&[
        input,
        "--decompose-rz",
        "--decompose-cz",
        "--epsilon",
        "1e-4",
    ]);
    assert!(report.at("options/decompose_rz").as_bool());
    assert!(report.at("options/decompose_cz").as_bool());
    assert_eq!(report.at("options/rz_epsilon").as_f64(), 1e-4);
    assert_eq!(
        report.at("metrics/output/rz").as_usize(),
        0,
        "--decompose-rz should leave no rotations behind"
    );

    let report = json_report(&[input, "--parallel", "-O1"]);
    assert!(report.at("options/parallel").as_bool());
}

/// Where the circuit came from and went, as the report describes it.
#[test]
fn the_report_names_its_input_and_output() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("out.qasm");
    let path = file.to_str().unwrap();

    let to_file = json_report(&[MOD5_4_QASM, "-o", path, "-O1"]);
    assert_eq!(to_file.at("input/path").as_str(), MOD5_4_QASM);
    assert!(!to_file.at("input/stdin").as_bool());
    assert_eq!(to_file.at("output/path").as_str(), path);
    assert!(!to_file.at("output/stdout").as_bool());
    assert_eq!(
        to_file.at("input/bytes").as_usize(),
        fs::metadata(MOD5_4_QASM).unwrap().len() as usize
    );

    let discarded = json_report(&[MOD5_4_QASM, "-O1"]);
    assert!(discarded.at("output/path").is_null());
    assert!(!discarded.at("output/stdout").as_bool());
}

/// A stdin run reports `stdin: true` with a null path, and sizes the input
/// from what actually arrived.
#[test]
fn a_stdin_run_is_reported_as_such() {
    let qasm = fs::read_to_string(MOD5_4_QASM).unwrap();
    let run = Tzap::new(&["-", "-O1", "--json", "-q"])
        .stdin(&qasm)
        .run()
        .ok("stdin --json");
    let report = Json::parse(&run.stdout);

    assert!(report.at("input/stdin").as_bool());
    assert!(report.at("input/path").is_null());
    assert_eq!(report.at("input/bytes").as_usize(), qasm.len());
    assert_eq!(
        report.at("metrics/input/gates").as_usize(),
        report_from_file().at("metrics/input/gates").as_usize(),
        "the same circuit through either door gives the same numbers"
    );
}

fn report_from_file() -> Json {
    json_report(&[MOD5_4_QASM, "-O1"])
}

/// A path that needs escaping produces valid JSON, not a broken document —
/// quotes, backslashes and spaces all survive the round trip.
#[test]
fn paths_needing_escapes_still_produce_valid_json() {
    let dir = tempfile::tempdir().unwrap();
    for name in [
        "a \"quoted\" name.qasm",
        "with spaces.qasm",
        "back\\slash.qasm",
    ] {
        let path = dir.path().join(name);
        if fs::write(&path, RZ_QASM).is_err() {
            // Some filesystems reject these names outright; the escaping is
            // what's under test, not the filesystem's tolerance.
            continue;
        }
        let report = json_report(&[path.to_str().unwrap(), "-O1"]);
        assert_eq!(
            report.at("input/path").as_str(),
            path.to_str().unwrap(),
            "the path must survive escaping intact"
        );
    }
}

/// Timings are present and consistent: the whole run is at least as long as
/// the parse it contains.
#[test]
fn the_timings_are_plausible() {
    let report = json_report(&[MOD5_4_QASM, "-O3"]);
    let total = report.get("seconds").as_f64();
    let parse = report.at("input/parse_seconds").as_f64();
    assert!(total > 0.0 && parse >= 0.0, "{total} {parse}");
    assert!(
        total >= parse,
        "the run ({total}s) can't be shorter than its own parse ({parse}s)"
    );
    assert!(
        total >= report.get("murms").arr()[0].get("seconds").as_f64(),
        "the run can't be shorter than its MURM load"
    );
}

/// `--json` without `-q` keeps the human report on stderr and the machine one
/// on stdout — the whole point of the split.
#[test]
fn json_and_human_output_occupy_different_streams() {
    let run = tzap(&[MOD5_4_QASM, "-O1", "--json"]).ok("--json");
    Json::parse(&run.stdout);
    assert!(run.stderr.contains("Final result"), "got: {}", run.stderr);
    assert!(
        !run.stderr.contains("\"tzap\":"),
        "the report must not also go to stderr:\n{}",
        run.stderr
    );
    assert!(
        !run.stdout.contains("Final result"),
        "the banner must not go to stdout:\n{}",
        run.stdout
    );
}

/// The report is what a run says about itself, and does not depend on
/// whether anything was rendered for a human: the same numbers arrive with
/// the commentary on and off.
#[test]
fn the_report_is_the_same_with_and_without_the_human_output() {
    let with_banner = tzap(&[MOD5_4_QASM, "-O3", "--json"]).ok("--json");
    let quiet = tzap(&[MOD5_4_QASM, "-O3", "--json", "-q"]).ok("--json -q");

    let a = Json::parse(&with_banner.stdout);
    let b = Json::parse(&quiet.stdout);
    assert_schema(&a, "--json");
    assert_schema(&b, "--json -q");
    for path in [
        "metrics/output/gates",
        "metrics/output/t",
        "metrics/baseline/gates",
    ] {
        assert_eq!(a.at(path).as_usize(), b.at(path).as_usize(), "{path}");
    }
    assert_eq!(
        a.get("fixpoints").arr()[0].get("rounds").as_usize(),
        b.get("fixpoints").arr()[0].get("rounds").as_usize()
    );
    assert_eq!(
        a.get("fixpoints").arr()[0].get("converged").as_bool(),
        b.get("fixpoints").arr()[0].get("converged").as_bool()
    );
}

/// A failing run writes no report at all — an exit code and a message on
/// stderr, with stdout left clean for the caller to distinguish.
#[test]
fn a_failed_run_writes_no_report() {
    for args in [
        vec!["definitely-missing.qasm", "--json"],
        vec![TEST_QASM, "--json", "--passes", "NoSuchPass"],
        vec![TEST_QASM, "--json", "--epsilon", "0"],
    ] {
        tzap(&args).failed(&format!("{args:?}"));
    }
}

/// `cache_dir` reports the location actually in force, so a report carries
/// enough to reproduce the run.
#[test]
fn the_report_names_the_cache_directory() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    let report = json_report(&[TEST_QASM, "-O1", "--cache-dir", path]);
    assert_eq!(report.get("cache_dir").as_str(), path);

    let run = Tzap::new(&[TEST_QASM, "-O1", "--json", "-q"])
        .env("TZAP_CACHE_DIR", path)
        .run()
        .ok("TZAP_CACHE_DIR --json");
    assert_eq!(Json::parse(&run.stdout).get("cache_dir").as_str(), path);
}
