use std::env;
use std::fs;
use std::io::{BufReader, BufWriter, Write};
use std::process;

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use tzap::circuit::{Circuit, Gate};
use tzap::decompose::DecomposeToffoli;
use tzap::decompose_rz::DecomposeRz;
use tzap::cancel::CancelPairs;
use tzap::pass::{Pass, count_t};
use tzap::phase_fold_global::PhaseFoldGlobal;
use tzap::phase_fold_global_expr::PhaseFoldGlobalExpr;
use tzap::qasm::{StreamingReader, serialize_gates};

fn fmt_num<N: std::fmt::Display>(n: N) -> String {
    let s = n.to_string();
    let is_negative = s.starts_with('-');
    let num_part = if is_negative { &s[1..] } else { &s[..] };
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    if is_negative { result.push('-'); }

    let rem = num_part.len() % 3;
    for (i, c) in num_part.chars().enumerate() {
        if i > 0 && i % 3 == rem {
            result.push(',');
        }
        result.push(c);
    }
    result
}

fn make_progress_bar(total: u64, name: &str) -> ProgressBar {
    let pb = ProgressBar::new(total);
    pb.set_style(ProgressStyle::with_template(
        &format!("  {name:<17} {{bar:40.green/dim}} {{percent:>3}}% ({{eta}})")
    ).unwrap().progress_chars("━╸─"));
    pb
}

fn run_passes_parallel(circuit: &Circuit, passes: &[&dyn Pass], num_chunks: usize) -> Circuit {
    let chunk_size = (circuit.gates.len() + num_chunks - 1) / num_chunks;
    let mut chunks: Vec<Circuit> = circuit.gates
        .chunks(chunk_size)
        .map(|slice| {
            let mut c = Circuit::new(circuit.num_qubits);
            for g in slice { c.apply(g.clone()); }
            c
        })
        .collect();
    for p in passes {
        chunks = chunks.par_iter()
            .map(|chunk| p.run(chunk))
            .collect();
    }
    let mut out = Circuit::new(circuit.num_qubits);
    for c in &chunks {
        for g in &c.gates { out.apply(g.clone()); }
    }
    out
}

fn run_streaming(
    input_path: &str,
    output_path: Option<&str>,
    batch_size: usize,
    cancel: bool,
    no_global: bool,
    expr: bool,
    decompose_rz: bool,
) {
    let file = fs::File::open(input_path).unwrap_or_else(|e| {
        eprintln!("Error reading {input_path}: {e}");
        process::exit(1);
    });
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);

    let mut stream = StreamingReader::new(reader).unwrap_or_else(|e| {
        eprintln!("Error parsing header: {e}");
        process::exit(1);
    });
    let num_qubits = stream.num_qubits;
    eprintln!("\t└─ {} qubits (streaming, batch {})\n", fmt_num(num_qubits), fmt_num(batch_size));

    let mut writer: Option<BufWriter<fs::File>> = output_path.map(|p| {
        let out_file = fs::File::create(p).unwrap_or_else(|e| {
            eprintln!("Error creating {p}: {e}");
            process::exit(1);
        });
        let mut w = BufWriter::new(out_file);
        writeln!(w, "OPENQASM 2.0;").unwrap();
        writeln!(w, "include \"qelib1.inc\";").unwrap();
        writeln!(w, "qreg q[{num_qubits}];").unwrap();
        w
    });

    let decompose = DecomposeToffoli;
    let cancel_pass = CancelPairs;
    let global = PhaseFoldGlobal;
    let global_expr = PhaseFoldGlobalExpr;
    let rz_decompose = DecomposeRz::default();
    let cancel_pass2 = CancelPairs;
    let global2 = PhaseFoldGlobal;

    let num_par_chunks = std::thread::available_parallelism()
        .map(|n| n.get() * 4)
        .unwrap_or(8);
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_par_chunks)
        .build_global()
        .ok();

    let start = std::time::Instant::now();
    let mut chunk_num = 0usize;
    let mut total_in_gates = 0usize;
    let mut total_out_gates = 0usize;
    let mut total_in_t = 0usize;
    let mut total_out_t = 0usize;

    loop {
        let batch = match stream.next_batch(batch_size) {
            Ok(Some(b)) => b,
            Ok(None) => break,
            Err(e) => { eprintln!("Error parsing: {e}"); process::exit(1); }
        };

        chunk_num += 1;
        let mut circuit = Circuit::new(num_qubits);
        for g in batch { circuit.apply(g); }
        eprint!("  Chunk {}: {} gates", fmt_num(chunk_num), fmt_num(circuit.gates.len()));

        let circuit = if circuit.has_toffoli { decompose.run(&circuit) } else { circuit };

        total_in_gates += circuit.gates.len();
        total_in_t += count_t(&circuit);

        let mut passes: Vec<&dyn Pass> = vec![];
        if cancel { passes.push(&cancel_pass); }
        if !no_global {
            if expr { passes.push(&global_expr); } else { passes.push(&global); }
        }
        let c = run_passes_parallel(&circuit, &passes, num_par_chunks);

        let c = if decompose_rz && c.gates.iter().any(|g| matches!(g, Gate::rz(..))) {
            let rz_passes: Vec<&dyn Pass> = vec![&rz_decompose, &cancel_pass2, &global2];
            run_passes_parallel(&c, &rz_passes, num_par_chunks)
        } else {
            c
        };

        total_out_gates += c.gates.len();
        total_out_t += count_t(&c);
        eprintln!(" → {} gates · {} T", fmt_num(c.gates.len()), fmt_num(count_t(&c)));

        if let Some(ref mut w) = writer {
            let gate_text = serialize_gates(&c.gates);
            w.write_all(gate_text.as_bytes()).unwrap_or_else(|e| {
                eprintln!("Error writing output: {e}"); process::exit(1);
            });
        }
    }

    if let Some(ref mut w) = writer {
        w.flush().unwrap();
    }
    let total = start.elapsed();

    let gate_pct = if total_in_gates > 0 {
        ((total_in_gates as f64 - total_out_gates as f64) / total_in_gates as f64) * 100.0
    } else { 0.0 };
    let t_pct = if total_in_t > 0 {
        ((total_in_t as f64 - total_out_t as f64) / total_in_t as f64) * 100.0
    } else { 0.0 };

    eprintln!("\n\x1b[1m  ⚡\u{FE0F} Result\x1b[0m (streaming, {} chunks)", fmt_num(chunk_num));
    eprintln!("\t├─ Gates  {} → {} (↓{:.1}%)", fmt_num(total_in_gates), fmt_num(total_out_gates), gate_pct);
    eprintln!("\t├─ T/Tdg  {} → {} (↓{:.1}%)", fmt_num(total_in_t), fmt_num(total_out_t), t_pct);
    eprintln!("\t└─ Time   {:.3}s", total.as_secs_f64());
    if let Some(p) = output_path {
        eprintln!("  wrote {p}");
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut input_path = None;
    let mut output_path = None;
    let mut cancel = false;
    let mut no_global = false;
    let mut expr = false;
    let mut decompose_rz = false;
    let mut to_cliffordt = false;
    let mut parallel = None;
    let mut streaming = false;
    let mut batch_size: usize = 20_000_000;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--help" | "-h" => {
                println!();
                println!("  \x1b[1m⚡\u{FE0F} tzap\x1b[0m  —  fast T-gate optimizer for Clifford+T circuits");
                println!();
                println!("  Decomposes Toffoli (ccx) gates into Clifford+T by default.");
                println!("  Pass --decompose-rz to also decompose Rz gates via gridsynth.");
                println!();
                println!("  \x1b[1;33mUSAGE\x1b[0m");
                println!("    tzap <input.qasm> [output.qasm] [options]");
                println!();
                println!("  \x1b[1;33mARGS\x1b[0m");
                println!("    \x1b[1m<input.qasm>\x1b[0m     Input OpenQASM 2.0 file");
                println!("    \x1b[1m[output.qasm]\x1b[0m    Output file (no output if omitted)");
                println!();
                println!("  \x1b[1;33mOPTIONS\x1b[0m");
                println!("    \x1b[1m-o\x1b[0m <file>        Write output to <file>");
                println!("    \x1b[1m--decompose-rz\x1b[0m   Decompose Rz gates into Clifford+T (gridsynth)");
                println!("    \x1b[1m--to-cliffordt\x1b[0m   Decompose ccx + Rz to Clifford+T, no optimization");
                println!("    \x1b[1m--cancel\x1b[0m         Enable the gate cancellation pass");
                println!("    \x1b[1m--expr\x1b[0m           Use expr-based phase folding (exact parity)");
                println!("    \x1b[1m--no-global\x1b[0m      Skip the global phase folding pass");
                println!("    \x1b[1m--parallel\x1b[0m       Force parallel mode");
                println!("    \x1b[1m--no-parallel\x1b[0m    Force sequential mode");
                println!("    \x1b[1m--streaming\x1b[0m      Stream-process in chunks (low memory, for very large files)");
                println!("    \x1b[1m--batch-size\x1b[0m <N> Gates per streaming chunk (default: 20000000)");
                println!("    \x1b[1m-h, --help\x1b[0m       Print this help message");
                println!();
                process::exit(0);
            }
            "--cancel" => cancel = true,
            "--no-global" => no_global = true,
            "--expr" => expr = true,
            "--decompose-rz" => decompose_rz = true,
            "--to-cliffordt" => to_cliffordt = true,
            "--parallel" => parallel = Some(true),
            "--no-parallel" => parallel = Some(false),
            "--streaming" => streaming = true,
            "--batch-size" => {
                i += 1;
                batch_size = args.get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_else(|| { eprintln!("--batch-size requires a number"); process::exit(1); });
            }
            "-o" => {
                i += 1;
                output_path = args.get(i).map(|s| s.as_str());
            }
            _ if args[i].starts_with('-') => {
                eprintln!("Unknown flag: {}", args[i]);
                process::exit(1);
            }
            _ => {
                if input_path.is_none() {
                    input_path = Some(args[i].as_str());
                } else if output_path.is_none() {
                    output_path = Some(args[i].as_str());
                }
            }
        }
        i += 1;
    }

    let Some(input_path) = input_path else {
        eprintln!("\x1b[1m⚡\u{FE0F} tzap\x1b[0m <input.qasm> [-o output.qasm] [--decompose-rz] [--to-cliffordt] [--cancel] [--no-global] [--expr] [--parallel] [--no-parallel]");
        process::exit(1);
    };

    eprintln!("\x1b[1m⚡\u{FE0F} tzap\x1b[0m");

    let file_size = fs::metadata(input_path).map(|m| m.len()).unwrap_or(0);
    eprintln!("  Parsing {input_path} ({:.1} MB)", file_size as f64 / (1024.0 * 1024.0));

    if streaming {
        run_streaming(input_path, output_path, batch_size, cancel, no_global, expr, decompose_rz);
        return;
    }

    let qasm = fs::read_to_string(input_path).unwrap_or_else(|e| {
        eprintln!("Error reading {input_path}: {e}");
        process::exit(1);
    });

    let parse_start = std::time::Instant::now();
    let circuit = Circuit::from_qasm(&qasm).unwrap_or_else(|e| {
        eprintln!("Error parsing {input_path}: {e}");
        process::exit(1);
    });
    let parse_time = parse_start.elapsed();
    let initial_gates = circuit.gates.len();
    let initial_t = count_t(&circuit);
    let num_qubits = circuit.num_qubits;
    eprintln!("\t└─ {} qubits · {} gates · {} T/Tdg · {:.3}s\n",
        fmt_num(num_qubits),
        fmt_num(initial_gates),
        fmt_num(initial_t),
        parse_time.as_secs_f64()
    );

    if to_cliffordt {
        let num_par_chunks = std::thread::available_parallelism()
            .map(|n| n.get() * 4)
            .unwrap_or(8);
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_par_chunks)
            .build_global()
            .ok();

        let start = std::time::Instant::now();
        let circuit = if circuit.has_toffoli {
            let pb = make_progress_bar(circuit.gates.len() as u64, "Toffoli → Clifford+T");
            let c = DecomposeToffoli.run_with_progress(&circuit, &pb);
            pb.finish_and_clear();
            eprintln!("  Toffoli → Clifford+T\n\t└─ {} gates", fmt_num(c.gates.len()));
            c
        } else {
            circuit
        };
        let circuit = if circuit.gates.iter().any(|g| matches!(g, Gate::rz(..))) {
            let rz = DecomposeRz::default();
            let pb = make_progress_bar(circuit.gates.len() as u64, rz.name());
            let c = rz.run_with_progress(&circuit, &pb);
            pb.finish_and_clear();
            eprintln!("  {}\n\t└─ {} gates · {} T", rz.name(), fmt_num(c.gates.len()), fmt_num(count_t(&c)));
            c
        } else {
            circuit
        };
        let total = start.elapsed();
        eprintln!("\n  {} gates · {:.3}s", fmt_num(circuit.gates.len()), total.as_secs_f64());

        if let Some(output_path) = output_path {
            let output = circuit.to_qasm();
            fs::write(output_path, &output).unwrap_or_else(|e| {
                eprintln!("Error writing {output_path}: {e}");
                process::exit(1);
            });
            eprintln!("  wrote {output_path}");
        }
        return;
    }

    let parallel = parallel.unwrap_or(initial_gates > 1_000_000);

    let num_par_chunks = std::thread::available_parallelism()
        .map(|n| n.get() * 4)
        .unwrap_or(8);

    if parallel || decompose_rz {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_par_chunks)
            .build_global()
            .ok();
    }

    let decompose = DecomposeToffoli;
    let cancel_pass = CancelPairs;
    let global = PhaseFoldGlobal;
    let global_expr = PhaseFoldGlobalExpr;
    let rz_decompose = DecomposeRz::default();

    // Run toffoli decomposition eagerly so we can use post-decomp counts as baseline.
    let circuit = if circuit.has_toffoli {
        let pb = make_progress_bar(circuit.gates.len() as u64, decompose.name());
        let pass_start = std::time::Instant::now();
        let c = decompose.run_with_progress(&circuit, &pb);
        let elapsed = pass_start.elapsed();
        pb.finish_and_clear();
        let t = count_t(&c);
        eprintln!("  {}\n\t└─ {} gates · {} T · {:.3}s", decompose.name(), fmt_num(c.gates.len()), fmt_num(t), elapsed.as_secs_f64());
        c
    } else {
        circuit
    };

    // Decompose Rz gates early so cancel/global passes can optimize the result.
    let circuit = if decompose_rz && circuit.gates.iter().any(|g| matches!(g, Gate::rz(..))) {
        let pb = make_progress_bar(circuit.gates.len() as u64, rz_decompose.name());
        let pass_start = std::time::Instant::now();
        let c = rz_decompose.run_with_progress(&circuit, &pb);
        let elapsed = pass_start.elapsed();
        pb.finish_and_clear();
        let t = count_t(&c);
        eprintln!("  {}\n\t└─ {} gates · {} T · {:.3}s", rz_decompose.name(), fmt_num(c.gates.len()), fmt_num(t), elapsed.as_secs_f64());
        c
    } else {
        circuit
    };
    let baseline_gates = circuit.gates.len();
    let baseline_t = count_t(&circuit);

    let mut passes: Vec<&dyn Pass> = vec![];
    if cancel { passes.push(&cancel_pass); }
    if !no_global {
        if expr {
            passes.push(&global_expr);
        } else {
            passes.push(&global);
        }
    }
    let mut total_pass_time = std::time::Duration::ZERO;
    let result = if parallel {
        eprintln!("  Parallel mode: {} chunks / {} threads\n",
            fmt_num(num_par_chunks), fmt_num(rayon::current_num_threads()));

        // Start with one chunk per slice of the circuit.
        let chunk_size = (circuit.gates.len() + num_par_chunks - 1) / num_par_chunks;
        let mut chunks: Vec<Circuit> = circuit.gates
            .chunks(chunk_size)
            .map(|slice| {
                let mut c = Circuit::new(circuit.num_qubits);
                for g in slice { c.apply(g.clone()); }
                c
            })
            .collect();

        // Run each pass across all chunks before moving to the next.
        for p in &passes {
            let total_gates: u64 = chunks.iter().map(|c| c.gates.len() as u64).sum();
            let pb = make_progress_bar(total_gates, p.name());
            let pass_start = std::time::Instant::now();
            chunks = chunks.par_iter()
                .map(|chunk| p.run_with_progress(chunk, &pb))
                .collect();
            let elapsed = pass_start.elapsed();
            total_pass_time += elapsed;
            pb.finish_and_clear();

            let total_gates: usize = chunks.iter().map(|c| c.gates.len()).sum();
            let total_t: usize = chunks.iter().map(|c| count_t(c)).sum();
            eprintln!("  {}\n\t└─ {} gates · {} T · {:.3}s", p.name(), fmt_num(total_gates), fmt_num(total_t), elapsed.as_secs_f64());
        }

        // Stitch final result.
        let mut out = Circuit::new(circuit.num_qubits);
        for c in &chunks {
            for g in &c.gates { out.apply(g.clone()); }
        }
        out
    } else {
        // Sequential: run each pass on the full circuit with logging.
        let mut c = circuit.clone();
        for p in &passes {
            let pb = make_progress_bar(c.gates.len() as u64, p.name());
            let pass_start = std::time::Instant::now();
            c = p.run_with_progress(&c, &pb);
            let elapsed = pass_start.elapsed();
            total_pass_time += elapsed;
            pb.finish_and_clear();
            let t = count_t(&c);
            eprintln!("  {}\n\t└─ {} gates · {} T · {:.3}s", p.name(), fmt_num(c.gates.len()), fmt_num(t), elapsed.as_secs_f64());
        }
        c
    };

    let final_gates = result.gates.len();
    let final_t = count_t(&result);

    let gate_pct = if baseline_gates > 0 {
        ((baseline_gates as f64 - final_gates as f64) / baseline_gates as f64) * 100.0
    } else {
        0.0
    };
    let t_pct = if baseline_t > 0 {
        ((baseline_t as f64 - final_t as f64) / baseline_t as f64) * 100.0
    } else {
        0.0
    };
    eprintln!("\n\x1b[1m  ⚡\u{FE0F} Result\x1b[0m");
    eprintln!("\t├─ Gates  {} → {} (↓{:.1}%)", fmt_num(baseline_gates), fmt_num(final_gates), gate_pct);
    eprintln!("\t├─ T/Tdg  {} → {} (↓{:.1}%)", fmt_num(baseline_t), fmt_num(final_t), t_pct);
    eprintln!("\t└─ Time   {:.3}s", total_pass_time.as_secs_f64());

    let input_has_rz = circuit.gates.iter().any(|g| matches!(g, Gate::rz(..)));
    let output_has_rz = result.gates.iter().any(|g| matches!(g, Gate::rz(..)));
    if output_has_rz && !input_has_rz {
        panic!("BUG: output contains Rz gates but input did not");
    }
    if output_has_rz && decompose_rz {
        panic!("BUG: output contains Rz gates after --decompose-rz");
    }

    if let Some(output_path) = output_path {
        let output = result.to_qasm();
        fs::write(output_path, &output).unwrap_or_else(|e| {
            eprintln!("Error writing {output_path}: {e}");
            process::exit(1);
        });
        eprintln!("  wrote {output_path}");
    }
}

