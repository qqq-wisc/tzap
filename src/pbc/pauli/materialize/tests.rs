use super::*;
use crate::circuit::{Circuit, Gate};
use crate::pbc::to_pbc;
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::time::{Duration, Instant};

// Independent dense evaluator retained only as a small-test oracle and an
// explicitly timed baseline for the previous prefix-wide algorithm.
fn dense(arena: &PauliArena, qubits: usize, last: usize) -> Vec<ExpandedPauli> {
    let mut values: Vec<ExpandedPauli> = Vec::new();
    for node in &arena.nodes[..=last] {
        let mut value = ExpandedPauli {
            phase: Phase::One,
            factors: vec![Pauli::I; qubits],
        };
        match *node {
            PauliNode::Identity => (),
            PauliNode::Single { qubit, pauli } => value.factors[qubit as usize] = pauli,
            PauliNode::Product(a, b) => {
                value.phase = values[a.node]
                    .phase
                    .times(a.phase)
                    .times(values[b.node].phase)
                    .times(b.phase);
                for (q, factor) in value.factors.iter_mut().enumerate() {
                    let (phase, p) = values[a.node].factors[q].times(values[b.node].factors[q]);
                    value.phase = value.phase.times(phase);
                    *factor = p;
                }
            }
        }
        values.push(value);
    }
    values
}

fn chain(n: usize) -> (PauliArena, PauliRef) {
    let mut arena = PauliArena::new();
    let mut root = arena.push(PauliNode::Single {
        qubit: 0,
        pauli: Pauli::X,
    });
    for q in 1..n {
        let leaf = arena.push(PauliNode::Single {
            qubit: q as u32,
            pauli: Pauli::X,
        });
        root = arena.push(PauliNode::Product(root, leaf));
    }
    (arena, root)
}

#[test]
fn sparse_materialization_matches_dense_random_shared_dags() {
    for seed in 0..200 {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut arena = PauliArena::new();
        let phases = [Phase::One, Phase::I, Phase::MinusOne, Phase::MinusI];
        for _ in 0..120 {
            let node = if rng.gen_bool(0.3) {
                PauliNode::Single {
                    qubit: rng.gen_range(0..5),
                    pauli: [Pauli::X, Pauli::Y, Pauli::Z][rng.gen_range(0..3)],
                }
            } else {
                let a = arena
                    .reference(rng.gen_range(0..arena.nodes.len()))
                    .scaled(phases[rng.gen_range(0..4)]);
                let b = arena
                    .reference(rng.gen_range(0..arena.nodes.len()))
                    .scaled(phases[rng.gen_range(0..4)]);
                PauliNode::Product(a, b)
            };
            arena.push(node);
        }
        let expected = dense(&arena, 5, arena.nodes.len() - 1);
        // Arbitrary order, duplicate roots, and scalar phases on root handles.
        let roots: Vec<_> = (0..40)
            .map(|_| {
                arena
                    .reference(rng.gen_range(0..arena.nodes.len()))
                    .scaled(phases[rng.gen_range(0..4)])
            })
            .collect();
        arena
            .materialize(&roots, 100_000, |reference, value| {
                let target = &expected[reference.node];
                assert_eq!(
                    value.phase.times(reference.phase),
                    target.phase.times(reference.phase),
                    "seed {seed}"
                );
                let factors: Vec<_> = target
                    .factors
                    .iter()
                    .enumerate()
                    .filter(|(_, p)| **p != Pauli::I)
                    .map(|(q, &p)| (q as u32, p))
                    .collect();
                assert_eq!(value.sorted_factors(), factors, "seed {seed}");
                Ok(())
            })
            .unwrap();
    }
}

#[test]
fn late_sparse_leaf_export_is_independent_of_arena_prefix_and_width() {
    for n in [1, 100, 10_000, 100_000] {
        let circuit = Circuit {
            num_qubits: n,
            num_cbits: 0,
            gates: vec![Gate::t((n - 1) as u32)],
        };
        let pbc = to_pbc(&circuit).unwrap();
        let root = pbc.operations()[0].axis().as_ref();
        let stats = pbc
            .arena
            .materialize(&[root], 3, |_, value| {
                assert_eq!(value.factors.len(), 1);
                Ok(())
            })
            .unwrap();
        assert_eq!((stats.nodes, stats.work), (1, 3));
        assert_eq!(
            pbc.to_text().unwrap(),
            format!("qubits {n}\nregisters 0\nr 1 1 Z{}\n", n - 1)
        );
    }
}

#[test]
fn growing_chain_moves_storage_with_exact_linear_work() {
    for n in [1, 32, 1_000, 10_000] {
        let (arena, root) = chain(n);
        let stats = arena
            .materialize(&[root], 4 * n, |_, value| {
                assert_eq!(
                    value.sorted_factors(),
                    (0..n as u32).map(|q| (q, Pauli::X)).collect::<Vec<_>>()
                );
                Ok(())
            })
            .unwrap();
        assert_eq!(stats.nodes, 2 * n - 1);
        assert_eq!(stats.copied_factors, 0);
        assert_eq!(stats.multiplied_factors, n - 1);
        assert_eq!(stats.work, 4 * n - 1);
    }
}

#[test]
fn repeated_roots_are_evaluated_once_and_budgeted_per_output() {
    let (arena, root) = chain(100);
    let roots = vec![root; 100];
    let stats = arena.materialize(&roots, 20_000, |_, _| Ok(())).unwrap();
    assert_eq!(stats.nodes, 199);
    assert_eq!(stats.copied_factors, 0);
    assert_eq!(stats.output_factors, 10_000);
    assert!(matches!(
        arena.materialize(&roots, stats.work - 1, |_, _| Ok(())),
        Err(PbcError::ExpansionLimit)
    ));
    assert!(arena.materialize(&roots, stats.work, |_, _| Ok(())).is_ok());
}

#[test]
fn growing_prefix_outputs_are_budgeted_for_their_genuinely_quadratic_size() {
    let (arena, _) = chain(100);
    let roots: Vec<_> = (1..arena.nodes.len())
        .filter(|&node| node == 1 || node % 2 == 1)
        .map(|node| arena.reference(node))
        .collect();
    let stats = arena.materialize(&roots, 30_000, |_, _| Ok(())).unwrap();
    assert_eq!(stats.nodes, 199);
    assert_eq!(stats.output_factors, 100 * 101 / 2);
    assert!(matches!(
        arena.materialize(&roots, 400, |_, _| Ok(())),
        Err(PbcError::ExpansionLimit)
    ));
}

#[test]
fn radix_order_handles_all_qubit_bytes() {
    let factors = (0..100u32)
        .map(|i| (i.wrapping_mul(0x8123_4567), Pauli::Z))
        .collect();
    let value = SparsePauli {
        phase: Phase::One,
        factors,
    };
    let mut expected: Vec<_> = value.factors.iter().map(|(&q, &p)| (q, p)).collect();
    expected.sort_unstable_by_key(|&(q, _)| q);
    assert_eq!(value.sorted_factors(), expected);
}

#[test]
fn large_support_cancellation_keeps_phase_and_releases_capacity() {
    let (mut arena, root) = chain(1_000);
    let squared = arena.push(PauliNode::Product(root.scaled(Phase::I), root));
    arena
        .materialize(&vec![squared; 100], 20_000, |_, value| {
            assert_eq!(value.phase, Phase::I);
            assert!(value.factors.is_empty());
            assert_eq!(value.factors.capacity(), 0);
            Ok(())
        })
        .unwrap();
}

fn median(mut run: impl FnMut(), repetitions: usize) -> Duration {
    run(); // Warm-up is excluded.
    let mut elapsed = Vec::new();
    for _ in 0..repetitions {
        let start = Instant::now();
        run();
        elapsed.push(start.elapsed());
    }
    elapsed.sort_unstable();
    elapsed[elapsed.len() / 2]
}

/// Reproduce with cargo test --release --lib benchmark_sparse_materialization
/// -- --ignored --nocapture --test-threads=1. Timings are descriptive; exact
/// work-count tests above enforce complexity without flaky timing thresholds.
#[test]
#[ignore = "release-mode scaling benchmark with a dense baseline"]
fn benchmark_sparse_materialization() {
    for n in [1_000, 2_000, 4_000] {
        let (arena, root) = chain(n);
        let sparse = median(
            || {
                std::hint::black_box(
                    arena
                        .materialize(&[root], usize::MAX, |_, value| {
                            std::hint::black_box(value.sorted_factors());
                            Ok(())
                        })
                        .unwrap(),
                );
            },
            5,
        );
        let legacy = median(
            || {
                std::hint::black_box(dense(&arena, n, root.node));
            },
            5,
        );
        println!(
            "chain n={n} sparse_us={} dense_us={} dense_cells={}",
            sparse.as_micros(),
            legacy.as_micros(),
            (root.node + 1) * n
        );
    }
    for n in [10_000, 20_000, 40_000, 80_000] {
        let (arena, root) = chain(n);
        let sparse = median(
            || {
                std::hint::black_box(
                    arena
                        .materialize(&[root], usize::MAX, |_, value| {
                            std::hint::black_box(value.sorted_factors());
                            Ok(())
                        })
                        .unwrap(),
                );
            },
            7,
        );
        let input = Circuit {
            num_qubits: n,
            num_cbits: 0,
            gates: vec![Gate::t((n - 1) as u32)],
        };
        let pbc = to_pbc(&input).unwrap();
        let export = median(
            || {
                std::hint::black_box(pbc.to_text().unwrap());
            },
            101,
        );
        println!(
            "scale n={n} chain_us={} sparse_leaf_export_ns={}",
            sparse.as_micros(),
            export.as_nanos()
        );
    }
    for n in [10_000, 20_000, 40_000, 80_000] {
        let gates = (0..n)
            .map(|i| {
                let q = (i / 4 % 32) as u32;
                let target = (q + 1) % 32;
                match i % 4 {
                    0 => Gate::h(q),
                    1 => Gate::cnot { control: q, target },
                    2 => Gate::t(q),
                    _ => Gate::tdg(target),
                }
            })
            .collect();
        let input = Circuit {
            num_qubits: 32,
            num_cbits: 0,
            gates,
        };
        let converted = to_pbc(&input).unwrap();
        let bytes = converted.to_text().unwrap().len();
        let roots: Vec<_> = converted
            .operations()
            .iter()
            .map(|op| op.axis().as_ref())
            .collect();
        let stats = converted
            .arena
            .materialize(&roots, 16_000_000, |_, _| Ok(()))
            .unwrap();
        let export = median(
            || {
                std::hint::black_box(converted.to_text().unwrap());
            },
            5,
        );
        let end_to_end = median(
            || {
                std::hint::black_box(to_pbc(&input).unwrap().to_text().unwrap());
            },
            5,
        );
        println!(
            "entangled gates={n} export_us={} convert_export_us={} bytes={bytes} work={}",
            export.as_micros(),
            end_to_end.as_micros(),
            stats.work
        );
    }
}
