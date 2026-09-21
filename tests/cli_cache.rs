//! Where tzap keeps its on-disk MURMs, and the two actions that
//! inspect and clear them.
//!
//! Every test here pins the cache to a directory it owns: the resolution
//! order is the thing under test, so a test that leaked into the developer's
//! real cache would be both wrong and destructive. `--clear-cache` in
//! particular is never invoked without an explicit location.

#[path = "support/mod.rs"]
mod support;

use std::fs;
use std::path::{Path, PathBuf};

use support::{Json, Tzap, assert_plain, tzap};

const TEST_QASM: &str = "tests/fixtures/test.qasm";

/// Tiny SuperOpt bounds: a real MURM, built and cached in milliseconds.
const TINY: [&str; 6] = [
    "--superopt-qubits",
    "2",
    "--superopt-window-gates",
    "4",
    "--superopt-murm-entries",
    "300",
];

/// The subdirectory the MURMs themselves live in, inside the cache root.
const MURMS: &str = "murm";

/// A run with no environment-derived cache location at all, so only what the
/// test sets can be found. `HOME` is removed too — several tests are about
/// what happens when it is the only thing left.
fn isolated(args: &[&str]) -> Tzap {
    Tzap::new(args).env_remove("HOME")
}

/// Optimize something with SuperOpt in the pipeline, so a MURM is built or
/// loaded, against the cache root `env` describes.
fn warm(dir_flag: Option<&Path>, env: &[(&str, &str)]) -> support::Run {
    let mut args = vec![TEST_QASM, "-q", "--passes", "SuperOpt"];
    args.extend_from_slice(&TINY);
    let dir = dir_flag.map(|d| d.to_str().unwrap().to_string());
    if let Some(dir) = &dir {
        args.push("--cache-dir");
        args.push(dir);
    }
    let mut command = isolated(&args);
    for (name, value) in env {
        command = command.env(name, value);
    }
    command.run().ok("warming the cache")
}

/// `--cache-info`'s reported directory, as the JSON form gives it.
fn reported_dir(env: &[(&str, &str)], args: &[&str]) -> Option<PathBuf> {
    let mut full = vec!["--cache-info", "--json"];
    full.extend_from_slice(args);
    let mut command = isolated(&full);
    for (name, value) in env {
        command = command.env(name, value);
    }
    let run = command.run().ok("--cache-info --json");
    let dir = Json::parse(&run.stdout).get("cache_dir").clone();
    match dir {
        Json::Null => None,
        dir => Some(PathBuf::from(dir.as_str())),
    }
}

fn murm_files(root: &Path) -> Vec<PathBuf> {
    let dir = root.join(MURMS);
    let Ok(entries) = fs::read_dir(&dir) else {
        return Vec::new();
    };
    let mut files: Vec<PathBuf> = entries
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "bin"))
        .collect();
    files.sort();
    files
}

// ---------------------------------------------------------------------------
// Resolution order: --cache-dir > TZAP_CACHE_DIR > XDG_CACHE_HOME > HOME
// ---------------------------------------------------------------------------

#[test]
fn the_cache_root_follows_the_documented_precedence() {
    let dir = tempfile::tempdir().unwrap();
    let flag = dir.path().join("from-flag");
    let env = dir.path().join("from-env");
    let xdg = dir.path().join("from-xdg");
    let home = dir.path().join("from-home");
    for path in [&flag, &env, &xdg, &home] {
        fs::create_dir_all(path).unwrap();
    }
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();

    // HOME alone: the XDG spec's default.
    assert_eq!(
        reported_dir(&[("HOME", &s(&home))], &[]),
        Some(home.join(".cache").join("tzap"))
    );
    // XDG_CACHE_HOME beats HOME.
    assert_eq!(
        reported_dir(&[("HOME", &s(&home)), ("XDG_CACHE_HOME", &s(&xdg))], &[]),
        Some(xdg.join("tzap"))
    );
    // TZAP_CACHE_DIR beats XDG, and is used as-is (it names tzap already).
    assert_eq!(
        reported_dir(
            &[
                ("HOME", &s(&home)),
                ("XDG_CACHE_HOME", &s(&xdg)),
                ("TZAP_CACHE_DIR", &s(&env)),
            ],
            &[]
        ),
        Some(env.clone())
    );
    // --cache-dir beats everything.
    assert_eq!(
        reported_dir(
            &[
                ("HOME", &s(&home)),
                ("XDG_CACHE_HOME", &s(&xdg)),
                ("TZAP_CACHE_DIR", &s(&env)),
            ],
            &["--cache-dir", &s(&flag)]
        ),
        Some(flag.clone())
    );
}

/// The Windows locations come after the XDG ones: a native Windows process
/// has no `$HOME` at all, and used to end up with no cache — rebuilding its
/// MURM on every single run — while being told it was cached for
/// future use.
#[test]
fn the_windows_locations_are_the_last_resort() {
    let dir = tempfile::tempdir().unwrap();
    let home = dir.path().join("from-home");
    let local = dir.path().join("from-localappdata");
    let profile = dir.path().join("from-userprofile");
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();

    // USERPROFILE alone: the same `.cache` layout as $HOME.
    assert_eq!(
        reported_dir(&[("USERPROFILE", &s(&profile))], &[]),
        Some(profile.join(".cache").join("tzap"))
    );
    // LOCALAPPDATA beats it, and is joined directly — it already names a
    // per-application directory.
    assert_eq!(
        reported_dir(
            &[("USERPROFILE", &s(&profile)), ("LOCALAPPDATA", &s(&local)),],
            &[]
        ),
        Some(local.join("tzap"))
    );
    // ...and $HOME beats both, so a Unix shell on Windows (MSYS, Git Bash)
    // keeps reading the MURMs it has already cached.
    assert_eq!(
        reported_dir(
            &[
                ("USERPROFILE", &s(&profile)),
                ("LOCALAPPDATA", &s(&local)),
                ("HOME", &s(&home)),
            ],
            &[]
        ),
        Some(home.join(".cache").join("tzap"))
    );
}

/// Per the XDG spec, `$XDG_CACHE_HOME` must be absolute; a relative one is
/// ignored rather than resolved against the working directory, which would
/// scatter caches wherever tzap happened to be run from.
#[test]
fn a_relative_xdg_cache_home_is_ignored() {
    let dir = tempfile::tempdir().unwrap();
    let home = dir.path().to_str().unwrap();
    assert_eq!(
        reported_dir(&[("HOME", home), ("XDG_CACHE_HOME", "relative/cache")], &[]),
        Some(dir.path().join(".cache").join("tzap"))
    );
}

/// An empty variable is unset, not a location of "".
#[test]
fn empty_environment_variables_are_treated_as_unset() {
    let dir = tempfile::tempdir().unwrap();
    let home = dir.path().to_str().unwrap();
    assert_eq!(
        reported_dir(
            &[
                ("HOME", home),
                ("XDG_CACHE_HOME", ""),
                ("TZAP_CACHE_DIR", "")
            ],
            &[]
        ),
        Some(dir.path().join(".cache").join("tzap"))
    );
}

/// With nothing to resolve, caching is simply off — and the run still works,
/// because the cache is a speed optimization and never load-bearing.
#[test]
fn a_run_with_no_cache_location_still_works() {
    assert_eq!(reported_dir(&[], &[]), None);

    let run = isolated(&["--cache-info"])
        .run()
        .ok("--cache-info with no HOME");
    assert!(
        run.stdout.contains("No cache directory"),
        "got: {}",
        run.stdout
    );

    let mut args = vec![TEST_QASM, "--passes", "SuperOpt"];
    args.extend_from_slice(&TINY);
    let run = isolated(&args)
        .run()
        .ok("optimizing with no cache location");
    assert!(run.stderr.contains("Building MURM"));
}

/// `--cache-dir` is where the MURM actually lands, not just what gets
/// reported.
#[test]
fn the_murm_is_written_under_the_chosen_root() {
    let dir = tempfile::tempdir().unwrap();
    assert!(murm_files(dir.path()).is_empty());

    warm(Some(dir.path()), &[]);
    let files = murm_files(dir.path());
    assert_eq!(files.len(), 1, "expected exactly one MURM, got {files:?}");
    let name = files[0].file_name().unwrap().to_string_lossy().into_owned();
    assert!(
        name.contains("q2_g3_e300"),
        "the file name should carry the bounds it was built for: {name}"
    );
    assert!(
        name.ends_with(&format!("tzap{}.bin", env!("CARGO_PKG_VERSION"))),
        "the file name should carry the crate version: {name}"
    );
}

/// The same for each environment-derived location.
#[test]
fn the_murm_is_written_under_each_environment_root() {
    let dir = tempfile::tempdir().unwrap();

    let env_root = dir.path().join("env");
    warm(None, &[("TZAP_CACHE_DIR", env_root.to_str().unwrap())]);
    assert_eq!(murm_files(&env_root).len(), 1);

    let xdg = dir.path().join("xdg");
    warm(None, &[("XDG_CACHE_HOME", xdg.to_str().unwrap())]);
    assert_eq!(murm_files(&xdg.join("tzap")).len(), 1);

    let home = dir.path().join("home");
    fs::create_dir_all(&home).unwrap();
    warm(None, &[("HOME", home.to_str().unwrap())]);
    assert_eq!(murm_files(&home.join(".cache").join("tzap")).len(), 1);
}

/// The second run reads the MURM instead of rebuilding it — the entire
/// reason the cache exists.
#[test]
fn a_warm_cache_is_loaded_rather_than_rebuilt() {
    let dir = tempfile::tempdir().unwrap();
    let mut args = vec![TEST_QASM, "--passes", "SuperOpt", "--cache-dir"];
    let path = dir.path().to_str().unwrap();
    args.push(path);
    args.extend_from_slice(&TINY);

    let cold = isolated(&args).run().ok("cold run");
    assert!(
        cold.stderr.contains("Building MURM"),
        "got: {}",
        cold.stderr
    );

    let warm = isolated(&args).run().ok("warm run");
    assert!(
        !warm.stderr.contains("Building MURM"),
        "the second run should not rebuild:\n{}",
        warm.stderr
    );
    assert!(warm.stderr.contains("Loaded MURM"), "got: {}", warm.stderr);
}

#[test]
fn a_structurally_valid_cache_with_body_corruption_is_rebuilt() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    let mut args = vec![TEST_QASM, "--passes", "SuperOpt", "--cache-dir", path];
    args.extend_from_slice(&TINY);

    isolated(&args).run().ok("initial cache build");
    let cache = murm_files(dir.path()).pop().expect("one cache file");
    let original = fs::read(&cache).unwrap();
    let mut corrupt = original.clone();
    let body_byte = corrupt.len() - 16;
    corrupt[body_byte] ^= 0x80;
    fs::write(&cache, &corrupt).unwrap();

    isolated(&args).run().ok("corrupt cache fallback");
    assert_eq!(
        fs::read(&cache).unwrap(),
        original,
        "a corrupt body should be replaced by a deterministic fresh MURM"
    );
    isolated(&args).run().ok("rebuilt cache warm read");
}

#[test]
fn concurrent_cold_writers_leave_one_valid_warm_murm() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    let mut args = vec![TEST_QASM, "--passes", "SuperOpt", "--cache-dir", path];
    args.extend_from_slice(&TINY);

    let first = isolated(&args).spawn();
    let second = isolated(&args).spawn();
    let first = support::Run::from(first.wait_with_output().unwrap()).ok("first concurrent writer");
    let second =
        support::Run::from(second.wait_with_output().unwrap()).ok("second concurrent writer");
    assert!(first.status.success() && second.status.success());
    assert_eq!(murm_files(dir.path()).len(), 1);
    let leftovers: Vec<_> = fs::read_dir(dir.path().join(MURMS))
        .unwrap()
        .flatten()
        .filter(|entry| entry.file_name().to_string_lossy().contains(".tmp."))
        .collect();
    assert!(
        leftovers.is_empty(),
        "temporary files leaked after the race"
    );

    let warm = isolated(&args)
        .run()
        .ok("warm read after concurrent writers");
    assert!(!warm.stderr.contains("Building MURM"), "{}", warm.stderr);
}

/// Pre-MURM caches used a different format and directory. They are never
/// mistaken for current MURMs; a current cache is built in the XDG location.
#[test]
fn obsolete_superopt_tables_are_not_read_as_murms() {
    let dir = tempfile::tempdir().unwrap();
    let home = dir.path();
    let legacy = home.join(".tzap").join("superopt-tables");
    fs::create_dir_all(&legacy).unwrap();
    fs::write(legacy.join("old-table.bin"), b"obsolete format").unwrap();

    let mut args = vec![TEST_QASM, "--passes", "SuperOpt"];
    args.extend_from_slice(&TINY);
    let run = isolated(&args)
        .env("HOME", home.to_str().unwrap())
        .run()
        .ok("obsolete cache ignored");
    assert!(
        run.stderr.contains("Building MURM"),
        "the obsolete table must not be treated as a MURM:\n{}",
        run.stderr
    );
    assert_eq!(murm_files(&home.join(".cache").join("tzap")).len(), 1);
}

// ---------------------------------------------------------------------------
// --cache-info
// ---------------------------------------------------------------------------

#[test]
fn cache_info_reports_the_murms_it_finds() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();

    let empty = isolated(&["--cache-info", "--cache-dir", path])
        .run()
        .ok("empty --cache-info");
    assert!(empty.stdout.contains(path), "got: {}", empty.stdout);
    assert!(
        empty.stdout.contains("0 cached MURMs"),
        "got: {}",
        empty.stdout
    );
    assert_plain(&empty.stdout, "--cache-info");
    assert!(empty.stderr.is_empty(), "got: {}", empty.stderr);

    warm(Some(dir.path()), &[]);
    let warm_info = isolated(&["--cache-info", "--cache-dir", path])
        .run()
        .ok("warm --cache-info");
    assert!(
        warm_info.stdout.contains("1 cached MURM ·"),
        "the count should be singular for one MURM:\n{}",
        warm_info.stdout
    );
    assert!(
        warm_info.stdout.contains("q2_g3_e300"),
        "each MURM should be listed by name:\n{}",
        warm_info.stdout
    );
}

#[test]
fn cache_info_json_lists_every_murm_with_its_size() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    warm(Some(dir.path()), &[]);

    let run = isolated(&["--cache-info", "--json", "--cache-dir", path])
        .run()
        .ok("--cache-info --json");
    let report = Json::parse(&run.stdout);
    assert_eq!(
        report.keys(),
        vec!["tzap", "cache_dir", "tables", "murms", "total_bytes"]
    );
    assert_eq!(report.get("tzap").as_str(), env!("CARGO_PKG_VERSION"));
    assert_eq!(report.get("cache_dir").as_str(), path);

    let murms = report.get("murms").arr();
    let tables = report.get("tables").arr();
    assert_eq!(murms.len(), 1);
    assert_eq!(tables.len(), murms.len());
    assert_eq!(murms[0].keys(), vec!["path", "bytes"]);
    let listed = PathBuf::from(murms[0].get("path").as_str());
    assert!(listed.exists(), "the listed path should exist: {listed:?}");
    assert_eq!(
        murms[0].get("bytes").as_usize(),
        fs::metadata(&listed).unwrap().len() as usize
    );
    assert_eq!(
        report.get("total_bytes").as_usize(),
        murms[0].get("bytes").as_usize()
    );
    assert_eq!(
        tables[0].get("path").as_str(),
        murms[0].get("path").as_str()
    );
    assert_eq!(
        tables[0].get("bytes").as_usize(),
        murms[0].get("bytes").as_usize()
    );
}

/// `--cache-info` is a query, so it answers on stdout even under `--quiet` —
/// quiet silences commentary, not the thing that was asked for.
#[test]
fn cache_info_answers_even_when_quiet() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    let run = isolated(&["--cache-info", "--cache-dir", path, "-q"])
        .run()
        .ok("--cache-info -q");
    assert!(run.stdout.contains(path), "got: {:?}", run.stdout);
    assert!(run.stderr.is_empty());
}

/// The cache actions take no circuit, and say so rather than silently
/// ignoring one.
#[test]
fn the_cache_actions_reject_an_input_circuit() {
    for flag in ["--cache-info", "--clear-cache"] {
        let dir = tempfile::tempdir().unwrap();
        let run = isolated(&[flag, TEST_QASM, "--cache-dir", dir.path().to_str().unwrap()])
            .run()
            .failed(&format!("{flag} with an input"));
        assert!(
            run.stderr.contains(flag) && run.stderr.contains(TEST_QASM),
            "the error should name both:\n{}",
            run.stderr
        );
    }
}

#[test]
fn the_two_cache_actions_are_mutually_exclusive() {
    let dir = tempfile::tempdir().unwrap();
    let run = isolated(&[
        "--cache-info",
        "--clear-cache",
        "--cache-dir",
        dir.path().to_str().unwrap(),
    ])
    .run()
    .failed("--cache-info --clear-cache");
    assert!(
        run.stderr.contains("--cache-info") && run.stderr.contains("--clear-cache"),
        "got: {}",
        run.stderr
    );
}

#[test]
fn cache_dir_requires_a_value() {
    tzap(&["--cache-info", "--cache-dir"]).failed("--cache-dir with no value");
    tzap(&["--cache-info", "--cache-dir="]).failed("--cache-dir with an empty value");
}

/// `--cache-dir` accepts the `=` spelling like every other value-taking flag.
#[test]
fn cache_dir_accepts_the_equals_spelling() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    let run = isolated(&[&format!("--cache-dir={path}"), "--cache-info"])
        .run()
        .ok("--cache-dir=");
    assert!(run.stdout.contains(path), "got: {}", run.stdout);
}

// ---------------------------------------------------------------------------
// --clear-cache
// ---------------------------------------------------------------------------

#[test]
fn clear_cache_removes_the_murms_and_reports_what_it_freed() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    warm(Some(dir.path()), &[]);
    assert_eq!(murm_files(dir.path()).len(), 1);

    let run = isolated(&["--clear-cache", "--cache-dir", path])
        .run()
        .ok("--clear-cache");
    assert!(
        run.stderr.contains("Removed 1 cached MURM"),
        "got: {}",
        run.stderr
    );
    assert!(run.stderr.contains("freed"), "got: {}", run.stderr);
    assert!(murm_files(dir.path()).is_empty(), "MURMs should be gone");

    // Clearing an already-empty cache is a no-op, not an error.
    let again = isolated(&["--clear-cache", "--cache-dir", path])
        .run()
        .ok("--clear-cache twice");
    assert!(again.stderr.contains("Removed 0 cached MURMs"));
}

/// A cache root can be shared with other tools (`$XDG_CACHE_HOME/tzap`, or
/// whatever an operator points `--cache-dir` at), so clearing touches only
/// the MURM files — never the directory, and never anything else in it.
#[test]
fn clear_cache_leaves_everything_that_is_not_a_murm_alone() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    warm(Some(dir.path()), &[]);

    let murms_dir = dir.path().join(MURMS);
    let stray = murms_dir.join("notes.txt");
    fs::write(&stray, "not a MURM").unwrap();
    let sibling = dir.path().join("something-else.bin");
    fs::write(&sibling, "another tool's cache").unwrap();

    isolated(&["--clear-cache", "--cache-dir", path])
        .run()
        .ok("--clear-cache with strays");

    assert!(murm_files(dir.path()).is_empty());
    assert!(stray.exists(), "a non-MURM file must survive");
    assert!(
        sibling.exists(),
        "a file outside the MURMs dir must survive"
    );
    assert!(murms_dir.exists(), "the directory itself must survive");
}

/// `--clear-cache --json` reports what it removed, for a script that wants
/// to know.
#[test]
fn clear_cache_json_lists_what_it_removed() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    warm(Some(dir.path()), &[]);
    let before = fs::metadata(&murm_files(dir.path())[0]).unwrap().len() as usize;

    let run = isolated(&["--clear-cache", "--json", "--cache-dir", path])
        .run()
        .ok("--clear-cache --json");
    let report = Json::parse(&run.stdout);
    assert_eq!(report.get("murms").arr().len(), 1);
    assert_eq!(report.get("total_bytes").as_usize(), before);
    assert_eq!(report.get("cache_dir").as_str(), path);
    assert!(murm_files(dir.path()).is_empty());
}

/// Its summary is commentary on an action, so `--quiet` silences it — while
/// the action itself still happens.
#[test]
fn clear_cache_is_silent_when_quiet() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    warm(Some(dir.path()), &[]);

    let run = isolated(&["--clear-cache", "--cache-dir", path, "-q"])
        .run()
        .ok("--clear-cache -q");
    assert_eq!(run.stderr, "");
    assert!(run.stdout.is_empty());
    assert!(
        murm_files(dir.path()).is_empty(),
        "it should still have run"
    );
}

/// Clearing the cache costs speed, never correctness: the next run rebuilds
/// and produces the same circuit.
#[test]
fn clearing_the_cache_does_not_change_the_output() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    let out = dir.path().join("out.qasm");

    let mut args = vec![
        TEST_QASM,
        "-o",
        out.to_str().unwrap(),
        "-q",
        "--passes",
        "SuperOpt",
        "--cache-dir",
        path,
    ];
    args.extend_from_slice(&TINY);

    isolated(&args).run().ok("first run");
    let first = fs::read_to_string(&out).unwrap();

    isolated(&["--clear-cache", "--cache-dir", path])
        .run()
        .ok("clear");
    isolated(&args).run().ok("run after clearing");
    assert_eq!(first, fs::read_to_string(&out).unwrap());
}

/// The cache location a run reports and the one it writes to are the same
/// thing — checked through `--json`, which is where a run states the location
/// it used.
#[test]
fn a_run_reports_the_same_cache_directory_it_uses() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().to_str().unwrap();
    let mut args = vec![
        TEST_QASM,
        "--json",
        "-q",
        "--passes",
        "SuperOpt",
        "--cache-dir",
    ];
    args.push(path);
    args.extend_from_slice(&TINY);

    let run = isolated(&args).run().ok("--json --cache-dir");
    assert_eq!(Json::parse(&run.stdout).get("cache_dir").as_str(), path);
    assert_eq!(murm_files(dir.path()).len(), 1);

    // With nothing to resolve, the report says so rather than naming a
    // directory that was never used.
    let mut args = vec![TEST_QASM, "--json", "-q", "--passes", "SuperOpt"];
    args.extend_from_slice(&TINY);
    let run = isolated(&args).run().ok("--json with no cache location");
    assert!(Json::parse(&run.stdout).get("cache_dir").is_null());
}
