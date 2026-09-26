//! Measure the PBC rotation optimizer's T-count reduction and speed.
//!
//! ```text
//! cargo run --release --example pbc_opt_eval -- [--level none|O1|O2|O3]
//!     [--decompose-ccx] [--no-frame] [--no-eager] [--lazy-cliffords] [--litinski] [--depth] [--window W]
//!     [--lookback L] [--rounds R] [--candidates C] <file.qasm | directory>...
//! ```
//!
//! For each circuit: optimize with tzap at `--level` (default none; with
//! `--decompose-ccx`, CCX/CCZ are decomposed before optimizing), convert
//! to PBC, count T, run the rotation optimizer, and count again. Circuits
//! with Rz or resets are skipped. Prints one CSV row per circuit, then totals.

use std::path::{Path, PathBuf};
use std::time::Instant;

use tzap::circuit::{Circuit, Gate};
use tzap::optimize::{Level, Options, optimize};
use tzap::pbc::{OptimizeOptions, Strategy, to_pbc};

fn collect(path: &Path, out: &mut Vec<PathBuf>) {
    if path.is_dir() {
        let mut entries: Vec<_> = std::fs::read_dir(path)
            .unwrap()
            .map(|e| e.unwrap().path())
            .collect();
        entries.sort();
        for entry in entries {
            collect(&entry, out);
        }
    } else if path.extension().is_some_and(|e| e == "qasm") {
        out.push(path.to_path_buf());
    }
}

fn main() {
    let mut level = None;
    let mut decompose_ccx = false;
    let mut options = OptimizeOptions::default();
    let mut files = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = || args.next().expect("flag value");
        match arg.as_str() {
            "--level" => {
                level = match value().as_str() {
                    "none" => None,
                    "O1" => Some(Level::O1),
                    "O2" => Some(Level::O2),
                    "O3" => Some(Level::O3),
                    other => panic!("unknown level {other}"),
                }
            }
            "--decompose-ccx" => decompose_ccx = true,
            "--no-frame" => options.clifford_to_frame = false,
            "--no-eager" => options.eager_swaps = false,
            "--lazy-cliffords" => options.lazy_cliffords = true,
            "--litinski" => options.strategy = Strategy::Litinski,
            "--depth" => options.measure_depth = true,
            "--window" => options.window = value().parse().unwrap(),
            "--lookback" => options.lookback = value().parse().unwrap(),
            "--rounds" => options.rounds = value().parse().unwrap(),
            "--candidates" => options.candidates = value().parse().unwrap(),
            path => collect(Path::new(path), &mut files),
        }
    }

    println!(
        "file,qubits,gates,t_gate_circuit,pbc_t_before,pbc_t_after,reduction_pct,merges,swaps,cliffords_to_frame,tzap_s,convert_s,optimize_s,t_depth_before,t_depth_after,weight_before,weight_after,weights_before,weights_after"
    );
    let (mut before, mut after, mut convert_time, mut optimize_time, mut count) =
        (0usize, 0usize, 0f64, 0f64, 0usize);
    let mut tzap_time = 0f64;
    for file in files {
        let text = std::fs::read_to_string(&file).unwrap();
        let Ok(mut circuit) = Circuit::from_qasm(&text) else {
            eprintln!("skip (parse): {}", file.display());
            continue;
        };
        if circuit
            .gates
            .iter()
            .any(|g| matches!(g, Gate::rz(..) | Gate::reset(_)))
        {
            eprintln!("skip (rz/reset): {}", file.display());
            continue;
        }
        let start = Instant::now();
        if let Some(level) = level {
            let options = Options {
                level,
                decompose_ccx,
                ..Options::default()
            };
            circuit = optimize(&circuit, &options).unwrap().0;
        }
        let tzap = start.elapsed().as_secs_f64();
        let t_gates = circuit
            .gates
            .iter()
            .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_)))
            .count();
        let start = Instant::now();
        let mut pbc = match to_pbc(&circuit) {
            Ok(pbc) => pbc,
            Err(e) => {
                eprintln!("skip (convert: {e}): {}", file.display());
                continue;
            }
        };
        let convert = start.elapsed().as_secs_f64();
        let weights_before = histogram(&pbc.rotation_weights(WEIGHT_BUDGET).unwrap_or_default());
        let start = Instant::now();
        let stats = match pbc.optimize_rotations(options) {
            Ok(stats) => stats,
            Err(e) => {
                eprintln!("skip ({e}): {}", file.display());
                continue;
            }
        };
        let elapsed = start.elapsed().as_secs_f64();
        let weights_after = histogram(&pbc.rotation_weights(WEIGHT_BUDGET).unwrap_or_default());
        let pct = if stats.t_before == 0 {
            0.0
        } else {
            100.0 * (stats.t_before - stats.t_after) as f64 / stats.t_before as f64
        };
        println!(
            "{},{},{},{},{},{},{:.2},{},{},{},{:.4},{:.4},{:.4},{},{},{},{},{},{}",
            file.display(),
            circuit.num_qubits,
            circuit.gates.len(),
            t_gates,
            stats.t_before,
            stats.t_after,
            pct,
            stats.merges,
            stats.swaps,
            stats.cliffords_to_frame,
            tzap,
            convert,
            elapsed,
            stats.t_depth_before,
            stats.t_depth_after,
            stats.weight_before,
            stats.weight_after,
            weights_before,
            weights_after
        );
        before += stats.t_before;
        after += stats.t_after;
        convert_time += convert;
        tzap_time += tzap;
        optimize_time += elapsed;
        count += 1;
    }
    let pct = if before == 0 {
        0.0
    } else {
        100.0 * (before - after) as f64 / before as f64
    };
    eprintln!(
        "TOTAL {count} circuits: PBC T {before} -> {after} ({pct:.2}% reduction), \
         tzap {tzap_time:.3}s, convert {convert_time:.3}s, optimize {optimize_time:.3}s"
    );
}

/// Sparse work budget for measuring per-rotation weights.
const WEIGHT_BUDGET: usize = 1 << 30;

/// A weight histogram as `weight:count` pairs joined by `;`.
fn histogram(weights: &[usize]) -> String {
    let mut counts = std::collections::BTreeMap::new();
    for &w in weights {
        *counts.entry(w).or_insert(0usize) += 1;
    }
    counts
        .iter()
        .map(|(w, c)| format!("{w}:{c}"))
        .collect::<Vec<_>>()
        .join(";")
}
