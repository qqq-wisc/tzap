//! Seeded differential fuzzing, like the bounded/ignored fuzzers elsewhere in
//! tzap. Failures print the complete input and an independently replayable seed.
use super::*;
use rand::{Rng, SeedableRng, rngs::StdRng, seq::SliceRandom};

fn random_gate(kind: GateKind, n: usize, rng: &mut StdRng) -> Gate {
    let mut wires: Vec<u32> = (0..n as u32).collect();
    wires.shuffle(rng);
    let q = wires[0];
    match kind {
        GateKind::H => Gate::h(q),
        GateKind::X => Gate::x(q),
        GateKind::Z => Gate::z(q),
        GateKind::S => Gate::s(q),
        GateKind::Sdg => Gate::sdg(q),
        GateKind::T => Gate::t(q),
        GateKind::Tdg => Gate::tdg(q),
        GateKind::Cx => Gate::cnot {
            control: q,
            target: wires[1],
        },
        GateKind::Cz => Gate::cz {
            control: q,
            target: wires[1],
        },
        GateKind::Ccx => Gate::ccx {
            control1: q,
            control2: wires[1],
            target: wires[2],
        },
        GateKind::Ccz => Gate::ccz {
            control1: q,
            control2: wires[1],
            target: wires[2],
        },
        _ => unreachable!("unitary fuzz palette"),
    }
}

fn fuzz_case(seed: u64) -> usize {
    // Width and optional-native subset are encoded in the seed so replay only
    // needs one value. Every permitted gate kind is forced into the input.
    let n = (seed % 4 + 1) as usize;
    let native_mask = (seed >> 2) & 7;
    let mut palette = vec![
        GateKind::H,
        GateKind::X,
        GateKind::Z,
        GateKind::S,
        GateKind::Sdg,
        GateKind::T,
        GateKind::Tdg,
    ];
    if n >= 2 {
        palette.push(GateKind::Cx);
    }
    if n >= 2 && native_mask & 1 != 0 {
        palette.push(GateKind::Cz);
    }
    if n >= 3 && native_mask & 2 != 0 {
        palette.push(GateKind::Ccx);
    }
    if n >= 3 && native_mask & 4 != 0 {
        palette.push(GateKind::Ccz);
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut c = input(
        n,
        palette
            .iter()
            .map(|&kind| random_gate(kind, n, &mut rng))
            .collect(),
    );
    c.gates.shuffle(&mut rng);
    for _ in 0..rng.gen_range(4..=20) {
        let kind = palette[rng.gen_range(0..palette.len())];
        c.gates.push(random_gate(kind, n, &mut rng));
    }
    let result = std::panic::catch_unwind(|| check(&c));
    assert!(
        result.is_ok(),
        "PBC fuzz mismatch: seed={seed:#x}, qubits={n}, native_mask={native_mask:#05b}\n\
        Replay: PBC_FUZZ_SEED={seed} PBC_FUZZ_CASES=1 cargo test --release fuzz_to_pbc -- --ignored --nocapture\n{c}"
    );
    n
}

#[test]
fn bounded_fuzz_to_pbc_covers_64_seeded_circuits() {
    let mut coverage = [[0; 8]; 4];
    for case in 0..64 {
        let seed = 0x5042_4300 + case;
        let n = fuzz_case(seed);
        coverage[n - 1][((seed >> 2) & 7) as usize] += 1;
    }
    // Each width/native-mask combination is exercised twice. Native gates are
    // enabled only when the selected width can accommodate their operands.
    assert_eq!(coverage, [[2; 8]; 4]);
}

/// Extended runs and single-seed replay:
/// PBC_FUZZ_CASES=1000 cargo test --release fuzz_to_pbc -- --ignored --nocapture
#[test]
#[ignore = "extended exact-unitary fuzzing; set PBC_FUZZ_CASES and PBC_FUZZ_SEED to replay"]
fn fuzz_to_pbc() {
    let count: u64 = std::env::var("PBC_FUZZ_CASES")
        .map(|s| s.parse().expect("decimal case count"))
        .unwrap_or(1000);
    let seed: u64 = std::env::var("PBC_FUZZ_SEED")
        .map(|s| {
            if let Some(hex) = s.strip_prefix("0x") {
                u64::from_str_radix(hex, 16).expect("hex seed")
            } else {
                s.parse().expect("decimal seed")
            }
        })
        .unwrap_or(0x5042_4400);
    for case in 0..count {
        fuzz_case(seed.wrapping_add(case));
    }
    eprintln!("PBC fuzz: {count} cases passed, starting seed {seed:#x}");
}
