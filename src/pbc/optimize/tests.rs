use super::*;
use crate::circuit::{Circuit, Gate};
use crate::pbc::{Pauli, Strategy, to_pbc};
use crate::semantics::{Limits, pbc_unitary};
use rand::{Rng, SeedableRng, rngs::StdRng};

mod weight_example;

/// A PBC circuit from rotations given as words ("XIZ") and angles in pi/8.
fn rotations(n: usize, layers: &[(&str, i64)]) -> PbcCircuit {
    let mut c = PbcCircuit::new(n, 0);
    for &(word, k) in layers {
        let (sign, word) = match word.strip_prefix('-') {
            Some(rest) => (Phase::MinusOne, rest),
            None => (Phase::One, word),
        };
        let mut product = c.identity().as_ref().scaled(sign);
        for (q, p) in word.chars().enumerate() {
            let p = match p {
                'I' => Pauli::I,
                'X' => Pauli::X,
                'Y' => Pauli::Y,
                'Z' => Pauli::Z,
                _ => panic!("bad Pauli"),
            };
            let single = c.single(q as u32, p).unwrap();
            product = c.product(product, single.as_ref()).unwrap();
        }
        let axis = c.hermitian_axis(product, 1 << 20).unwrap();
        c.rotate(axis, PauliAngle::new(k)).unwrap();
    }
    c
}

use crate::pbc::Phase;

fn optimize(c: &mut PbcCircuit, options: OptimizeOptions) -> OptimizeStats {
    let before = pbc_unitary(c, Limits::default()).unwrap();
    let stats = c.optimize_rotations(options).unwrap();
    let after = pbc_unitary(c, Limits::default()).unwrap();
    assert!(
        before.equivalent_up_to_global_phase(&after),
        "optimization changed the unitary"
    );
    assert_eq!(stats.t_after, c.t_count());
    assert!(stats.t_after <= stats.t_before);
    stats
}

#[test]
fn commuting_rotations_merge_and_anticommuting_ones_block() {
    // Z0 T, X1 T, Z0 T: X1 commutes with Z0, so the two T merge into S.
    let mut c = rotations(2, &[("ZI", 1), ("IX", 1), ("ZI", 1)]);
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!((stats.t_before, stats.t_after, stats.merges), (3, 1, 1));
    // X0 anticommutes with Z0: no merge.
    let mut c = rotations(1, &[("Z", 1), ("X", 1), ("Z", 1)]);
    assert_eq!(optimize(&mut c, OptimizeOptions::default()).t_after, 3);
}

#[test]
fn signs_are_absorbed_so_opposite_rotations_cancel() {
    let mut c = rotations(2, &[("ZZ", 1), ("-ZZ", 1)]);
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!((stats.t_after, stats.rotations_after), (0, 0));
    assert!(c.operations().is_empty());
}

#[test]
fn identity_axes_are_removed() {
    let mut c = rotations(1, &[("I", 1), ("Z", 3)]);
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!(
        (stats.t_before, stats.t_after, stats.rotations_after),
        (2, 1, 1)
    );
}

/// The design note's first example: A = Z0,Z1; B = X0X1,Y0Y1; C = Z0,Z1.
/// The certified A,B swap exposes two merges: T count 6 -> 2.
#[test]
fn design_note_mcr_example_reduces_six_to_two() {
    let layers = [
        ("ZI", 1),
        ("IZ", 1),
        ("XX", 1),
        ("YY", 1),
        ("ZI", 1),
        ("IZ", 1),
    ];
    let mut c = rotations(2, &layers);
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!((stats.t_before, stats.t_after, stats.swaps), (6, 2, 1));
    // Without swaps, ordinary merging cannot cross the XX/YY block.
    let mut c = rotations(2, &layers);
    let no_swaps = OptimizeOptions {
        window: 0,
        ..OptimizeOptions::default()
    };
    assert_eq!(optimize(&mut c, no_swaps).t_after, 6);
}

/// Inverting the YY rotation breaks the certificate: no swap.
#[test]
fn design_note_negative_variant_is_rejected() {
    let mut c = rotations(
        2,
        &[
            ("ZI", 1),
            ("IZ", 1),
            ("XX", 1),
            ("YY", -1),
            ("ZI", 1),
            ("IZ", 1),
        ],
    );
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!((stats.t_after, stats.swaps), (6, 0));
}

/// The four-versus-four example: every cross pair anticommutes, yet the
/// generators commute. T count 12 -> 4.
#[test]
fn design_note_four_versus_four_example_reduces_twelve_to_four() {
    let a = [("ZII", 1), ("IZI", 1), ("ZIZ", 1), ("IZZ", 1)];
    let b = [("XXI", 1), ("YYI", 1), ("XXZ", 1), ("YYZ", 1)];
    let layers: Vec<_> = a.iter().chain(&b).chain(&a).copied().collect();
    let mut c = rotations(3, &layers);
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!((stats.t_before, stats.t_after, stats.swaps), (12, 4, 1));
}

#[test]
fn lookback_bounds_merge_distance() {
    let mut layers = vec![("ZII", 1)];
    layers.extend([("IXI", 2); 5]);
    layers.push(("ZII", 1));
    let mut c = rotations(3, &layers);
    let short = OptimizeOptions {
        lookback: 0,
        clifford_to_frame: false,
        ..OptimizeOptions::default()
    };
    // The five X1 rotations merge among themselves (adjacent), but Z0 cannot
    // look back past the survivor (which, as a Clifford, would otherwise move
    // into the frame and unblock it).
    assert_eq!(optimize(&mut c, short).t_after, 2);
    let mut c = rotations(3, &layers);
    assert_eq!(optimize(&mut c, OptimizeOptions::default()).t_after, 0);
}

#[test]
fn failure_leaves_the_circuit_unchanged() {
    let mut c = rotations(1, &[("Z", 1), ("Z", 1)]);
    let tiny = OptimizeOptions {
        max_packed_words: 0,
        ..OptimizeOptions::default()
    };
    assert_eq!(c.optimize_rotations(tiny), Err(PbcError::ExpansionLimit));
    assert_eq!(c.operations().len(), 2);
}

/// Measurements are barriers: rotations never merge across them, and the
/// exact channel is preserved.
#[test]
fn measurements_are_barriers_and_channels_are_preserved() {
    use crate::semantics::channel::{ChannelLimits, pbc_channel};
    let input = Circuit {
        num_qubits: 2,
        num_cbits: 1,
        gates: vec![
            Gate::t(0),
            Gate::h(1),
            Gate::t(1),
            Gate::measure { qubit: 1, cbit: 0 },
            Gate::t(0),
            Gate::h(1),
            Gate::t(1),
            Gate::t(0),
        ],
    };
    let mut c = to_pbc(&input).unwrap();
    let limits = ChannelLimits::default();
    let before = pbc_channel(&c, &[false], limits).unwrap();
    let stats = c.optimize_rotations(OptimizeOptions::default()).unwrap();
    // T0 before the measurement cannot merge with the two after it; those two
    // merge with each other (Z0 commutes with the rotation between them).
    assert_eq!((stats.t_before, stats.t_after), (5, 3));
    assert_eq!(pbc_channel(&c, &[false], limits).unwrap(), before);
    assert_eq!(c.measurement_count(), 1);
}

fn random_word(n: usize, rng: &mut StdRng) -> String {
    let mut w: String = (0..n)
        .map(|_| ['I', 'X', 'Y', 'Z'][rng.gen_range(0..4)])
        .collect();
    if rng.gen_bool(0.3) {
        w.insert(0, '-');
    }
    w
}

/// Random rotation sequences over a small axis pool (so merges and MCR
/// patterns occur), checked exactly for every window and lookback setting.
#[test]
fn seeded_fuzz_preserves_unitaries_and_never_increases_t() {
    let (mut reductions, mut swaps) = (0, 0);
    for case in 0..24u64 {
        let mut rng = StdRng::seed_from_u64(0x4d43_5200 + case);
        let n = 1 + (case % 3) as usize;
        let pool: Vec<String> = (0..rng.gen_range(2..=5))
            .map(|_| random_word(n, &mut rng))
            .collect();
        let layers: Vec<(String, i64)> = (0..rng.gen_range(4..=24))
            .map(|_| {
                (
                    pool[rng.gen_range(0..pool.len())].clone(),
                    rng.gen_range(-3..=4),
                )
            })
            .collect();
        let layers: Vec<(&str, i64)> = layers.iter().map(|(w, k)| (w.as_str(), *k)).collect();
        for options in [
            OptimizeOptions::default(),
            OptimizeOptions {
                window: 0,
                ..OptimizeOptions::default()
            },
            OptimizeOptions {
                window: 2,
                lookback: 1,
                rounds: 1,
                ..OptimizeOptions::default()
            },
            OptimizeOptions {
                strategy: Strategy::Litinski,
                ..OptimizeOptions::default()
            },
            OptimizeOptions {
                lazy_cliffords: true,
                ..OptimizeOptions::default()
            },
        ] {
            let mut c = rotations(n, &layers);
            let stats = optimize(&mut c, options);
            reductions += stats.t_before - stats.t_after;
            swaps += stats.swaps;
        }
    }
    eprintln!("fuzz: {reductions} T removed, {swaps} swaps");
    assert!(reductions > 0);
}

/// Converted gate circuits, including CCX/CCZ, through the whole pipeline.
#[test]
fn seeded_fuzz_on_converted_circuits() {
    for case in 0..64u64 {
        let mut rng = StdRng::seed_from_u64(0x4d43_5300 + case);
        let n = 1 + (case % 3) as usize;
        let mut input = Circuit {
            num_qubits: n,
            num_cbits: 0,
            gates: vec![],
        };
        for _ in 0..rng.gen_range(4..=30) {
            let q = rng.gen_range(0..n as u32);
            let r = (q + 1) % n as u32;
            let s = (q + 2) % n as u32;
            input.gates.push(match rng.gen_range(0..9) {
                0 => Gate::h(q),
                1 => Gate::s(q),
                2 => Gate::t(q),
                3 => Gate::tdg(q),
                4 | 5 if n >= 2 => Gate::cnot {
                    control: q,
                    target: r,
                },
                6 if n >= 2 => Gate::cz {
                    control: q,
                    target: r,
                },
                7 if n >= 3 => Gate::ccx {
                    control1: q,
                    control2: r,
                    target: s,
                },
                _ => Gate::t(q),
            });
        }
        for strategy in [Strategy::Merge, Strategy::Litinski] {
            let mut c = to_pbc(&input).unwrap();
            let options = OptimizeOptions {
                strategy,
                measure_depth: true,
                ..OptimizeOptions::default()
            };
            let stats = optimize(&mut c, options);
            crate::semantics::test_support::assert_equivalent(&input, &c);
            assert!(stats.t_after <= stats.t_before);
        }
    }
}

/// Plant the design note's MCR pattern (A = Z_a, Z_b; B = X_aX_b, Y_aY_b;
/// C = A) on random qubits, with a random letter permutation, random common
/// angles, an optional spectator factor, and random noise rotations around
/// it, so the swap path is exercised and checked exactly.
#[test]
fn seeded_fuzz_on_planted_mcr_patterns() {
    let mut swaps = 0;
    for case in 0..32u64 {
        let mut rng = StdRng::seed_from_u64(0x4d43_5400 + case);
        let n = 2 + (case % 2) as usize;
        let mut letters = ['X', 'Y', 'Z'];
        for i in (1..3).rev() {
            letters.swap(i, rng.gen_range(0..=i));
        }
        let [lx, ly, lz] = letters;
        let (qa, qb) = (0, 1);
        let spectator = (n == 3 && rng.gen_bool(0.5)).then_some(2);
        let word = |pairs: &[(usize, char)]| -> String {
            let mut w = vec!['I'; n];
            for &(q, p) in pairs {
                w[q] = p;
            }
            if let Some(s) = spectator {
                w[s] = lz;
            }
            w.into_iter().collect()
        };
        let (a, b) = (rng.gen_range(-3..=4), rng.gen_range(-3..=4));
        let group_a = [(word(&[(qa, lz)]), a), (word(&[(qb, lz)]), a)];
        let group_b = [
            (word(&[(qa, lx), (qb, lx)]), b),
            (word(&[(qa, ly), (qb, ly)]), b),
        ];
        let mut layers: Vec<(String, i64)> = Vec::new();
        for _ in 0..rng.gen_range(0..3) {
            layers.push((random_word(n, &mut rng), rng.gen_range(-3..=4)));
        }
        layers.extend(group_a.iter().cloned());
        layers.extend(group_b.iter().cloned());
        layers.extend(group_a.iter().cloned());
        for _ in 0..rng.gen_range(0..3) {
            layers.push((random_word(n, &mut rng), rng.gen_range(-3..=4)));
        }
        let layers: Vec<(&str, i64)> = layers.iter().map(|(w, k)| (w.as_str(), *k)).collect();
        let mut c = rotations(n, &layers);
        optimize(&mut c, OptimizeOptions::default());
        // Without frame moves, even-angle patterns still need the swap.
        let mut c = rotations(n, &layers);
        let no_frame = OptimizeOptions {
            clifford_to_frame: false,
            ..OptimizeOptions::default()
        };
        swaps += optimize(&mut c, no_frame).swaps;
    }
    eprintln!("planted: {swaps} swaps");
    assert!(swaps >= 8, "only {swaps} swaps");
}

/// Z0 T; X0 S; Y0 T. Moving the S (a Clifford) into the frame conjugates the
/// later Y0 into -Z0, which then cancels the first T: T count 2 -> 0. The
/// frame absorbs the S, so the unitary is unchanged.
#[test]
fn cliffords_move_to_the_frame_and_unblock_merges() {
    let layers = [("Z", 1), ("X", 2), ("Y", 1)];
    let mut c = rotations(1, &layers);
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!(
        (stats.t_before, stats.t_after, stats.cliffords_to_frame),
        (2, 0, 1)
    );
    assert!(c.operations().is_empty());
    let mut c = rotations(1, &layers);
    let kept = OptimizeOptions {
        clifford_to_frame: false,
        ..OptimizeOptions::default()
    };
    assert_eq!(optimize(&mut c, kept).t_after, 2);
}

/// Every Clifford rotation leaves the operation list, and measurement axes
/// are conjugated so the exact channel is preserved.
#[test]
fn clifford_to_frame_conjugates_measurements() {
    use crate::semantics::channel::{ChannelLimits, pbc_channel};
    let mut c = rotations(2, &[("XZ", 2), ("ZI", 1), ("YY", 4)]);
    let z = c.z(0).unwrap();
    c.measure(z, None).unwrap();
    let x = c.x(1).unwrap();
    c.rotate(x, PauliAngle::new(2)).unwrap();
    let limits = ChannelLimits::default();
    let before = pbc_channel(&c, &[], limits).unwrap();
    let stats = c.optimize_rotations(OptimizeOptions::default()).unwrap();
    assert_eq!(stats.cliffords_to_frame, 3);
    assert!(c.operations().iter().all(|op| match op {
        PbcOp::Rotate { angle, .. } => angle.eighths() % 2 == 1,
        _ => true,
    }));
    assert_eq!(pbc_channel(&c, &[], limits).unwrap(), before);
}

/// A = Z0 T, Z1 T; B = XX T, YY T; C = Z0 S, Z1 S. The swap merges A with C
/// without lowering T (T + S is still odd). Eager swaps take it; otherwise
/// the swap is skipped.
#[test]
fn eager_swaps_accept_merges_that_keep_t() {
    let layers = [
        ("ZI", 1),
        ("IZ", 1),
        ("XX", 1),
        ("YY", 1),
        ("ZI", 2),
        ("IZ", 2),
    ];
    let no_frame = OptimizeOptions {
        clifford_to_frame: false,
        ..OptimizeOptions::default()
    };
    let mut c = rotations(2, &layers);
    let eager = optimize(&mut c, no_frame);
    assert_eq!(
        (eager.swaps, eager.rotations_after, eager.t_after),
        (1, 4, 4)
    );
    let mut c = rotations(2, &layers);
    let strict = optimize(
        &mut c,
        OptimizeOptions {
            eager_swaps: false,
            ..no_frame
        },
    );
    assert_eq!(
        (strict.swaps, strict.rotations_after, strict.t_after),
        (0, 6, 4)
    );
}

fn depth(n: usize, layers: &[(&str, i64)]) -> (usize, usize) {
    let mut c = rotations(n, layers);
    let stats = c
        .optimize_rotations(OptimizeOptions {
            measure_depth: true,
            ..OptimizeOptions::default()
        })
        .unwrap();
    (stats.t_depth_before, stats.t_depth_after)
}

/// T depth places each odd rotation one layer past the latest anticommuting
/// one: Z0, Z1 share a layer; Z0, X0 do not; in Z0, X0, Z1 the Z1 joins the
/// first layer. Even rotations do not count.
#[test]
fn t_depth_is_as_soon_as_possible_layering() {
    assert_eq!(depth(2, &[("ZI", 1), ("IZ", 1)]).0, 1);
    assert_eq!(depth(1, &[("Z", 1), ("X", 1)]).0, 2);
    assert_eq!(depth(2, &[("ZI", 1), ("XI", 1), ("IZ", 1)]).0, 2);
    assert_eq!(depth(2, &[("ZI", 1), ("XI", 2), ("ZI", -3)]).0, 1);
    assert_eq!(depth(1, &[("Z", 1), ("X", 1), ("Y", 1)]).0, 3);
}

/// Litinski's layering combines equal rotations that share a layer into a
/// Clifford, which goes to the frame; a rotation blocked by an anticommuting
/// layer does not move.
#[test]
fn litinski_combines_equal_rotations_meeting_in_a_layer() {
    let layers = [("ZI", 1), ("IX", 1), ("XI", 1), ("ZX", 1), ("ZI", 1)];
    let mut c = rotations(2, &layers);
    let options = OptimizeOptions {
        strategy: Strategy::Litinski,
        measure_depth: true,
        ..OptimizeOptions::default()
    };
    let stats = optimize(&mut c, options);
    // Z0 and X0 anticommute, so the last Z0 cannot pass X0: nothing combines.
    assert_eq!(stats.t_after, 5);
    let mut c = rotations(2, &[("ZI", 1), ("IX", 1), ("ZI", 1)]);
    let stats = optimize(&mut c, options);
    assert_eq!(
        (stats.t_before, stats.t_after, stats.t_depth_after),
        (3, 1, 1)
    );
}

/// Weight counts non-identity factors over all rotation axes: X0Y1, Z1, Z1
/// weigh 4. The two Z1 merge into a Clifford that moves to the frame; it
/// passes no later rotation, so only X0Y1 remains, weight 2.
#[test]
fn weight_counts_non_identity_factors() {
    let mut c = rotations(2, &[("XY", 1), ("IZ", 1), ("IZ", 1)]);
    let stats = optimize(&mut c, OptimizeOptions::default());
    assert_eq!((stats.weight_before, stats.weight_after), (4, 2));
}

/// Three T about one axis: greedily, the first two become an S in the frame
/// and the third stays a T; lazily, all three merge into one r 3 before
/// anything reaches the frame. Either way one T remains.
#[test]
fn lazy_cliffords_merge_before_reaching_the_frame() {
    let layers = [("ZZ", 1), ("ZZ", 1), ("ZZ", 1), ("XI", 1)];
    let mut c = rotations(2, &layers);
    let greedy = optimize(&mut c, OptimizeOptions::default());
    assert_eq!((greedy.t_after, greedy.cliffords_to_frame), (2, 1));
    let mut c = rotations(2, &layers);
    let lazy = OptimizeOptions {
        lazy_cliffords: true,
        ..OptimizeOptions::default()
    };
    let stats = optimize(&mut c, lazy);
    assert_eq!((stats.t_after, stats.cliffords_to_frame), (2, 0));
    assert!(c.operations().iter().any(|op| matches!(
        op,
        PbcOp::Rotate { angle, .. } if angle.eighths() == 3
    )));
}
