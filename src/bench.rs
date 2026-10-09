#[cfg(test)]
mod tests {
    use std::time::Instant;

    use crate::cancel::CancelGates;
    use crate::circuit::{Circuit, Gate, Qubit};
    use crate::decompose::{DecomposeCz, DecomposeToffoli};
    use crate::pass::{Pass, count_t};
    use crate::phase_fold_rand::phase_fold_rand;
    use crate::unitary::circuits_equiv;
    struct Rng(u64);

    impl Rng {
        fn new(seed: u64) -> Self {
            Rng(seed)
        }

        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn range(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }

        /// [`Rng::range`] as a qubit index.
        fn qubit(&mut self, n: usize) -> Qubit {
            self.range(n) as Qubit
        }
    }

    fn other_qubit(rng: &mut Rng, q: Qubit, num_qubits: usize) -> Qubit {
        debug_assert!(num_qubits >= 2);
        (q + 1 + rng.qubit(num_qubits - 1)) % num_qubits as Qubit
    }

    fn random_circuit(rng: &mut Rng, num_qubits: usize, num_gates: usize) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for _ in 0..num_gates {
            let q = rng.qubit(num_qubits);
            let kind = if num_qubits < 3 {
                // skip ccx (needs 3 qubits) and cnot (needs 2 qubits) slots
                let k = if num_qubits < 2 {
                    rng.range(14) + 8
                } else {
                    rng.range(22)
                };
                if k < 8 { k } else { k + 2 } // skip 8..=9 (ccx)
            } else {
                rng.range(24)
            };
            match kind {
                0..=4 => {
                    let t = other_qubit(rng, q, num_qubits);
                    c.apply(Gate::cnot {
                        control: q,
                        target: t,
                    });
                }
                5 => {
                    let other = other_qubit(rng, q, num_qubits);
                    c.apply(Gate::cz {
                        control: q,
                        target: other,
                    });
                }
                6..=7 => c.apply(Gate::rz_f64(0.1 * rng.range(60) as f64, q).unwrap()),
                8..=9 => {
                    let mut qs = [q, 0, 0];
                    qs[1] = other_qubit(rng, q, num_qubits);
                    loop {
                        qs[2] = rng.qubit(num_qubits);
                        if qs[2] != qs[0] && qs[2] != qs[1] {
                            break;
                        }
                    }
                    c.apply(Gate::ccx {
                        control1: qs[0],
                        control2: qs[1],
                        target: qs[2],
                    });
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

    /// Random circuit dominated by single-qubit H/X/S/Z gates — the regime
    /// that drives the cancel pass's pair cancellation and Hadamard reduction.
    /// A small tail of sdg/t/tdg/cnot also exercises the non-reducible
    /// branches (a `t` makes a run non-Clifford, a `cnot` breaks a run).
    fn random_hxsz_circuit(rng: &mut Rng, num_qubits: usize, num_gates: usize) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for _ in 0..num_gates {
            let q = rng.qubit(num_qubits);
            match rng.range(24) {
                0..=4 => c.apply(Gate::h(q)),
                5..=9 => c.apply(Gate::x(q)),
                10..=14 => c.apply(Gate::s(q)),
                15..=19 => c.apply(Gate::z(q)),
                20 => c.apply(Gate::sdg(q)),
                21 => c.apply(Gate::t(q)),
                22 => c.apply(Gate::tdg(q)),
                23 => {
                    if num_qubits >= 2 {
                        let t = other_qubit(rng, q, num_qubits);
                        if rng.range(2) == 0 {
                            c.apply(Gate::cnot {
                                control: q,
                                target: t,
                            });
                        } else {
                            c.apply(Gate::cz {
                                control: q,
                                target: t,
                            });
                        }
                    } else {
                        c.apply(Gate::z(q));
                    }
                }
                _ => unreachable!(),
            }
        }
        c
    }

    /// Random circuits with a deliberately high CZ density. This stresses
    /// symmetric CZ matching, lookahead cancellation through commuting gates,
    /// conservative blockers, phase folding through CZ, and explicit lowering.
    fn random_cz_circuit(rng: &mut Rng, num_qubits: usize, num_gates: usize) -> Circuit {
        debug_assert!(num_qubits >= 2);
        let mut c = Circuit::new(num_qubits);
        for _ in 0..num_gates {
            let q = rng.qubit(num_qubits);
            let other = other_qubit(rng, q, num_qubits);
            match rng.range(24) {
                0..=7 => {
                    // Reverse half the emitted operand orders to exercise CZ symmetry.
                    if rng.range(2) == 0 {
                        c.apply(Gate::cz {
                            control: q,
                            target: other,
                        });
                    } else {
                        c.apply(Gate::cz {
                            control: other,
                            target: q,
                        });
                    }
                }
                8..=10 => c.apply(Gate::cnot {
                    control: q,
                    target: other,
                }),
                11 | 12 => c.apply(Gate::h(q)),
                13 | 14 => c.apply(Gate::x(q)),
                15 | 16 => c.apply(Gate::t(q)),
                17 => c.apply(Gate::tdg(q)),
                18 => c.apply(Gate::s(q)),
                19 => c.apply(Gate::sdg(q)),
                20 => c.apply(Gate::z(q)),
                21..=23 => c.apply(Gate::rz_f64(0.05 * rng.range(80) as f64, q).unwrap()),
                _ => unreachable!(),
            }
        }
        c
    }

    /// Deterministic build for profiling (not random).
    fn build_circuit(num_qubits: usize, num_gates: usize) -> Circuit {
        let mut c = Circuit::new(num_qubits);
        for i in 0..num_gates {
            let q = (i % num_qubits) as Qubit;
            let width = num_qubits as Qubit;
            match i % 24 {
                0..=5 => {
                    let target = (q + 1) % width;
                    c.apply(Gate::cnot { control: q, target });
                }
                6..=7 => c.apply(Gate::rz_f64(0.123 * (i as f64), q).unwrap()),
                8..=9 => {
                    let c1 = q;
                    let c2 = (q + 1) % width;
                    let t = (q + 2) % width;
                    c.apply(Gate::ccx {
                        control1: c1,
                        control2: c2,
                        target: t,
                    });
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

        let optimized = time("phase_fold", || phase_fold_rand(&circuit));
        println!(
            "  result: {} -> {} gates",
            circuit.gates.len(),
            optimized.gates.len()
        );

        println!("=== done ===\n");
    }

    #[test]
    #[ignore] // long-running: 10k random circuits with unitary equivalence checks
    fn fuzz_phase_fold_rand() {
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
            let cancelled = CancelGates.run(&decomposed);
            let optimized = phase_fold_rand(&cancelled);

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

        println!(
            "\n{:>5} {:>4}q {:>6} {:>6} {:>7}",
            "case", "", "T before", "T after", "reduced"
        );
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
    #[ignore] // long-running: 10k single-qubit-heavy circuits through the cancel pass
    fn fuzz_cancel_hadamard() {
        // Stress the cancel pass (pair cancellation + Hadamard reduction) on
        // circuits dominated by single-qubit H/X/S/Z sequences. Few qubits so
        // gates pile up per wire, producing long reducible/cancellable runs.
        let mut rng = Rng::new(0x4861_6441);
        let num_cases = 10_000;
        let count_h = |c: &Circuit| c.gates.iter().filter(|g| matches!(g, Gate::h(_))).count();
        let mut gates_before = 0;
        let mut gates_after = 0;
        let mut h_before = 0;
        let mut h_after = 0;

        for i in 0..num_cases {
            let num_qubits = rng.range(4) + 1; // 1..=4
            let num_gates = rng.range(191) + 10; // 10..=200
            let circuit = random_hxsz_circuit(&mut rng, num_qubits, num_gates);
            let cancelled = CancelGates.run(&circuit);

            assert!(
                circuits_equiv(&circuit, &cancelled, 1e-9),
                "MISMATCH on case {i}: {num_qubits} qubits, {num_gates} gates\n{circuit}"
            );
            // The pass never grows the circuit or its Hadamard count.
            assert!(
                cancelled.gates.len() <= circuit.gates.len(),
                "case {i}: gate count grew"
            );
            assert!(
                count_h(&cancelled) <= count_h(&circuit),
                "case {i}: Hadamard count grew"
            );
            // The pass runs to a fixpoint — a second run changes nothing.
            let twice = CancelGates.run(&cancelled);
            assert_eq!(
                twice.gates.len(),
                cancelled.gates.len(),
                "case {i}: cancel pass is not idempotent"
            );

            gates_before += circuit.gates.len();
            gates_after += cancelled.gates.len();
            h_before += count_h(&circuit);
            h_after += count_h(&cancelled);
        }

        println!("\nfuzz cancel (H/X/S/Z): {num_cases} cases, all equivalent");
        println!(
            "gates: {gates_before} -> {gates_after} ({:.1}% removed)",
            (1.0 - gates_after as f64 / gates_before as f64) * 100.0,
        );
        println!(
            "H gates: {h_before} -> {h_after} ({:.1}% removed)",
            (1.0 - h_after as f64 / h_before as f64) * 100.0,
        );
    }

    #[test]
    #[ignore] // long-running: 10k CZ-dense circuits through the full native-CZ pipeline
    fn fuzz_cz_pipeline() {
        let mut rng = Rng::new(0xC2_C2_C2_C2);
        let num_cases = 10_000;
        let mut total_cz_before = 0;
        let mut total_cz_after = 0;

        for i in 0..num_cases {
            let num_qubits = rng.range(4) + 2; // 2..=5
            let num_gates = rng.range(111) + 10; // 10..=120
            let circuit = random_cz_circuit(&mut rng, num_qubits, num_gates);
            let label = format!("case {i}: {num_qubits} qubits, {num_gates} gates");
            let cancelled = check_cz_pipeline_case(&circuit, &label);

            total_cz_before += count_cz(&circuit);
            total_cz_after += count_cz(&cancelled);
        }

        println!(
            "\nfuzz CZ pipeline: {num_cases} cases, CZ {total_cz_before} -> {total_cz_after}, all equivalent"
        );
    }

    #[test]
    fn cz_pipeline_fuzz_smoke() {
        let mut rng = Rng::new(0xC2_5A_0C_E5);
        for i in 0..256 {
            let num_qubits = rng.range(3) + 2; // 2..=4
            let num_gates = rng.range(31) + 10; // 10..=40
            let circuit = random_cz_circuit(&mut rng, num_qubits, num_gates);
            check_cz_pipeline_case(&circuit, &format!("smoke case {i}"));
        }
    }

    fn count_cz(circuit: &Circuit) -> usize {
        circuit
            .gates
            .iter()
            .filter(|g| matches!(g, Gate::cz { .. }))
            .count()
    }

    fn check_cz_pipeline_case(circuit: &Circuit, label: &str) -> Circuit {
        let cancelled = CancelGates.run(circuit);
        let optimized = phase_fold_rand(&cancelled);
        let decomposed = DecomposeCz.run(&optimized);

        for (stage, result) in [
            ("cancellation", &cancelled),
            ("phase folding", &optimized),
            ("CZ decomposition", &decomposed),
        ] {
            assert!(
                circuits_equiv(circuit, result, 1e-9),
                "{stage} mismatch on {label}\n{circuit}"
            );
        }
        assert!(
            cancelled.gates.len() <= circuit.gates.len(),
            "{label}: cancellation grew the circuit"
        );
        assert!(
            count_cz(&cancelled) <= count_cz(circuit),
            "{label}: cancellation increased the CZ count"
        );
        assert!(
            !decomposed
                .gates
                .iter()
                .any(|g| matches!(g, Gate::cz { .. })),
            "{label}: DecomposeCz left a native CZ"
        );

        let twice = CancelGates.run(&cancelled);
        assert_eq!(
            twice.to_qasm(),
            cancelled.to_qasm(),
            "{label}: CZ cancellation is not idempotent"
        );
        cancelled
    }

    #[test]
    #[ignore] // long-running: unitary equivalence checks on benchmark circuits
    fn verify_benchmark_circuits() {
        // End-to-end soundness check on real (measurement-free) benchmark
        // circuits, capped at qubit counts the full-unitary check can handle.
        let names = [
            "tof_3",
            "tof_4",
            "tof_5",
            "barenco_tof_3",
            "barenco_tof_4",
            "barenco_tof_5",
            "mod5_4",
            "mod_mult_55",
            "vbe_adder_3",
            "grover_5",
        ];
        for name in names {
            let path = format!(
                "{}/benchmarks/feynman/{name}.qasm",
                env!("CARGO_MANIFEST_DIR")
            );
            let src = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
            let circuit = crate::qasm::parse(&src).unwrap_or_else(|e| panic!("parse {name}: {e}"));
            let decomposed = DecomposeToffoli.run(&circuit);
            let cancelled = CancelGates.run(&decomposed);
            let optimized = phase_fold_rand(&cancelled);
            assert!(
                circuits_equiv(&decomposed, &optimized, 1e-9),
                "{name}: optimized circuit is not equivalent to the input"
            );
            assert!(
                count_t(&optimized) <= count_t(&decomposed),
                "{name}: T count increased"
            );
            println!(
                "{name}: {} -> {} T",
                count_t(&decomposed),
                count_t(&optimized)
            );
        }
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
            let mut optimized = phase_fold_rand(&circuit);

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
            let opt_a = phase_fold_rand(&a);
            let opt_b = phase_fold_rand(&b);

            assert!(
                !circuits_equiv(&opt_a, &opt_b, 1e-10),
                "BUG: independently optimized circuits became equivalent on case {i}\n\
                 a ({num_qubits}q, {num_gates}g):\n{a}\nopt_a:\n{opt_a}\n\
                 b ({num_qubits}q, {num_gates}g):\n{b}\nopt_b:\n{opt_b}"
            );

            // also check each optimization is correct
            assert!(circuits_equiv(&a, &opt_a, 1e-10), "opt_a != a on case {i}");
            assert!(circuits_equiv(&b, &opt_b, 1e-10), "opt_b != b on case {i}");
        }

        println!(
            "\nfuzz inequiv: {num_cases} cases, {false_equiv} coincidentally equivalent, {} confirmed inequivalent",
            num_cases - false_equiv
        );
    }
}
