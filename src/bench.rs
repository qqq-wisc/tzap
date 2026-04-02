#[cfg(test)]
mod tests {
    use std::time::Instant;

    use indicatif::ProgressBar;

    use crate::circuit::{Circuit, Gate};
    use crate::cancel::CancelPairs;
    use crate::decompose::DecomposeToffoli;
    use crate::pass::{Pass, count_t};
    use crate::phase_fold_global::phase_fold_global;
    use crate::unitary::circuits_equiv;
    struct Rng(u64);

    impl Rng {
        fn new(seed: u64) -> Self { Rng(seed) }

        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn range(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn random_circuit(rng: &mut Rng, num_qubits: usize, num_gates: usize) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for _ in 0..num_gates {
            let q = rng.range(num_qubits);
            let kind = if num_qubits < 3 {
                // skip ccx (needs 3 qubits) and cnot (needs 2 qubits) slots
                let k = if num_qubits < 2 { rng.range(14) + 8 } else { rng.range(22) };
                if k < 8 { k } else { k + 2 } // skip 8..=9 (ccx)
            } else {
                rng.range(24)
            };
            match kind {
                0..=5 => {
                    let t = (q + 1 + rng.range(num_qubits - 1)) % num_qubits;
                    c.apply(Gate::cnot { control: q, target: t });
                }
                6..=7 => c.apply(Gate::rz(0.1 * rng.range(60) as f64, q)),
                8..=9 => {
                    let mut qs = [q, 0, 0];
                    qs[1] = (q + 1 + rng.range(num_qubits - 1)) % num_qubits;
                    loop {
                        qs[2] = rng.range(num_qubits);
                        if qs[2] != qs[0] && qs[2] != qs[1] { break; }
                    }
                    c.apply(Gate::ccx { control1: qs[0], control2: qs[1], target: qs[2] });
                }
                10..=13 => c.apply(Gate::t(q)),
                14..=16 => c.apply(Gate::tdg(q)),
                17 | 18 => c.apply(Gate::h(q)),
                19 | 20 => c.apply(Gate::s(q)),
                21 => c.apply(Gate::sdg(q)),
                22 => c.apply(Gate::x(q)),
                23 => c.apply(Gate::z(q)),
                _ => unreachable!(),
            }
        }
        c
    }

    /// Deterministic build for profiling (not random).
    fn build_circuit(num_qubits: usize, num_gates: usize) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for i in 0..num_gates {
            let q = i % num_qubits;
            match i % 24 {
                0..=5 => {
                    let target = (q + 1) % num_qubits;
                    c.apply(Gate::cnot { control: q, target });
                }
                6..=7 => c.apply(Gate::rz(0.123 * (i as f64), q)),
                8..=9 => {
                    let c1 = q;
                    let c2 = (q + 1) % num_qubits;
                    let t = (q + 2) % num_qubits;
                    c.apply(Gate::ccx { control1: c1, control2: c2, target: t });
                }
                10..=13 => c.apply(Gate::t(q)),
                14..=16 => c.apply(Gate::tdg(q)),
                17 | 18 => c.apply(Gate::h(q)),
                19 | 20 => c.apply(Gate::s(q)),
                21 => c.apply(Gate::sdg(q)),
                22 => c.apply(Gate::x(q)),
                23 => c.apply(Gate::z(q)),
                _ => unreachable!(),
            }
        }
        c
    }

    fn time<F: FnOnce() -> T, T>(label: &str, f: F) -> T {
        let start = Instant::now();
        let result = f();
        let elapsed = start.elapsed();
        println!("  {label}: {elapsed:?}");
        result
    }

    #[test]
    #[ignore] // long-running: 1M gate benchmark
    fn profile_1000_gates() {
        let num_qubits = 10;
        let num_gates = 1_000_000;

        println!("\n=== Profile: {num_gates} gates on {num_qubits} qubits ===");

        let circuit = time("build circuit", || build_circuit(num_qubits, num_gates));
        println!("  circuit: {} gates", circuit.gates.len());

        let optimized = time("phase_fold", || {
            phase_fold_global(&circuit, &ProgressBar::hidden())
        });
        println!(
            "  result: {} -> {} gates",
            circuit.gates.len(),
            optimized.gates.len()
        );

        println!("=== done ===\n");
    }

    #[test]
    #[ignore] // long-running: 10k random circuits with unitary equivalence checks
    fn fuzz_phase_fold_global() {
        let mut rng = Rng::new(0xDEAD_BEEF);
        let num_cases = 10_000;
        let mut total_t_before = 0;
        let mut total_t_after = 0;
        let mut reductions: Vec<(usize, usize, usize, usize, usize)> = Vec::new(); // (qubits, t_before, t_after, gates_after, case)

        for i in 0..num_cases {
            let num_qubits = rng.range(6) + 1; // 1..=6
            let num_gates = rng.range(991) + 10; // 10..=1000
            let circuit = random_circuit(&mut rng, num_qubits, num_gates);
            let decomposed = DecomposeToffoli.run(&circuit);
            let cancelled = CancelPairs.run(&decomposed);
            let optimized = phase_fold_global(&cancelled, &ProgressBar::hidden());

            assert!(
                circuits_equiv(&circuit, &optimized, 1e-10),
                "MISMATCH on case {i}: {num_qubits} qubits, {num_gates} gates\n{circuit}"
            );

            let t_before = count_t(&decomposed);
            let t_after = count_t(&optimized);
            total_t_before += t_before;
            total_t_after += t_after;
            if t_before != t_after {
                reductions.push((num_qubits, t_before, t_after, optimized.gates.len(), i));
            }
        }

        reductions.sort_by(|a, b| {
            let pct_a = (a.1 - a.2) as f64 / a.1 as f64;
            let pct_b = (b.1 - b.2) as f64 / b.1 as f64;
            pct_b.partial_cmp(&pct_a).unwrap()
        });

        println!("\n{:>5} {:>4}q {:>6} {:>6} {:>7}", "case", "", "T before", "T after", "reduced");
        println!("{}", "-".repeat(40));
        for (q, t_before, t_after, _gates, case) in &reductions {
            let removed = t_before - t_after;
            let pct = removed as f64 / *t_before as f64 * 100.0;
            println!("{case:>5} {q:>4}q {t_before:>6} {t_after:>6} {removed:>6} ({pct:.0}%)");
        }
        println!("{}", "-".repeat(40));
        println!(
            "{} cases with T reductions out of {num_cases} ({:.0}%)",
            reductions.len(),
            reductions.len() as f64 / num_cases as f64 * 100.0,
        );
        println!(
            "total T: {} -> {} ({:.1}% reduction)",
            total_t_before,
            total_t_after,
            (1.0 - total_t_after as f64 / total_t_before as f64) * 100.0,
        );
    }

    #[test]
    #[ignore] // long-running: 100 mutation detection tests
    fn fuzz_mutation_detected() {
        let mut rng = Rng::new(0xFEED_FACE);
        let num_cases = 100;
        let mut caught = 0;

        for _ in 0..num_cases {
            let num_qubits = rng.range(5) + 2;
            let num_gates = rng.range(91) + 10;
            let circuit = random_circuit(&mut rng, num_qubits, num_gates);
            let mut optimized = phase_fold_global(&circuit, &ProgressBar::hidden());

            // inject a bug: append X q0
            optimized.apply(Gate::x(0));

            if !circuits_equiv(&circuit, &optimized, 1e-10) {
                caught += 1;
            }
        }

        println!("\nmutation test: {caught}/{num_cases} mutations detected");
        assert_eq!(caught, num_cases, "some mutations were not detected");
    }

    #[test]
    #[ignore] // long-running: 1000 random circuit pairs
    fn fuzz_inequivalent_circuits() {
        let mut rng = Rng::new(0xCAFE_BABE);
        let num_cases = 1000;
        let mut false_equiv = 0;

        for i in 0..num_cases {
            let num_qubits = rng.range(5) + 2; // 2..=6
            let num_gates = rng.range(41) + 10; // 10..=50

            let a = random_circuit(&mut rng, num_qubits, num_gates);
            let b = random_circuit(&mut rng, num_qubits, num_gates);

            if circuits_equiv(&a, &b, 1e-10) {
                // two independent random circuits happened to be equivalent — rare but possible
                false_equiv += 1;
                continue;
            }

            // optimize both independently
            let opt_a = phase_fold_global(&a, &ProgressBar::hidden());
            let opt_b = phase_fold_global(&b, &ProgressBar::hidden());

            assert!(
                !circuits_equiv(&opt_a, &opt_b, 1e-10),
                "BUG: independently optimized circuits became equivalent on case {i}\n\
                 a ({num_qubits}q, {num_gates}g):\n{a}\nopt_a:\n{opt_a}\n\
                 b ({num_qubits}q, {num_gates}g):\n{b}\nopt_b:\n{opt_b}"
            );

            // also check each optimization is correct
            assert!(
                circuits_equiv(&a, &opt_a, 1e-10),
                "opt_a != a on case {i}"
            );
            assert!(
                circuits_equiv(&b, &opt_b, 1e-10),
                "opt_b != b on case {i}"
            );
        }

        println!(
            "\nfuzz inequiv: {num_cases} cases, {false_equiv} coincidentally equivalent, {} confirmed inequivalent",
            num_cases - false_equiv
        );
    }

}
