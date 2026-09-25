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

/// A random unitary circuit and its RNG, for callers that extend it. Width and
/// optional-native subset are encoded in the seed so replay only needs one
/// value. Every permitted gate kind is forced into the input.
fn random_circuit(seed: u64) -> (Circuit, StdRng) {
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
    (c, rng)
}

fn fuzz_case(seed: u64) -> usize {
    let (c, _) = random_circuit(seed);
    let (n, native_mask) = (c.num_qubits, (seed >> 2) & 7);
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
    let (count, seed) = fuzz_range(1000, 0x5042_4400);
    for case in 0..count {
        fuzz_case(seed.wrapping_add(case));
    }
    eprintln!("PBC fuzz: {count} cases passed, starting seed {seed:#x}");
}

/// `PBC_FUZZ_CASES` (decimal) and `PBC_FUZZ_SEED` (decimal or 0x-hex), with
/// the given defaults.
fn fuzz_range(default_count: u64, default_seed: u64) -> (u64, u64) {
    let count = std::env::var("PBC_FUZZ_CASES")
        .map(|s| s.parse().expect("decimal case count"))
        .unwrap_or(default_count);
    let seed = std::env::var("PBC_FUZZ_SEED")
        .map(|s| match s.strip_prefix("0x") {
            Some(hex) => u64::from_str_radix(hex, 16).expect("hex seed"),
            None => s.parse().expect("decimal seed"),
        })
        .unwrap_or(default_seed);
    (count, seed)
}

/// A random circuit from [`random_circuit`] with one to three measurements at
/// random positions, so most are followed by more gates. Targets include
/// overwritten registers; the last register is never written, so its initial
/// value must survive. Four-qubit circuits get one measurement, keeping exact
/// Choi blocks within the oracle's limits.
fn random_measured_circuit(seed: u64) -> (Circuit, Vec<bool>) {
    let (mut c, mut rng) = random_circuit(seed);
    let n = c.num_qubits;
    c.num_cbits = n + 1;
    let count = if n == 4 { 1 } else { rng.gen_range(1..=3) };
    for _ in 0..count {
        let position = rng.gen_range(0..=c.gates.len());
        let measure = Gate::measure {
            qubit: rng.gen_range(0..n as u32),
            cbit: rng.gen_range(0..n as u32),
        };
        c.gates.insert(position, measure);
    }
    let initial = (0..c.num_cbits).map(|_| rng.gen_bool(0.5)).collect();
    (c, initial)
}

/// Check conversion, optionally after the optimizer at `levels`, against the
/// exact channel of the input circuit.
fn mid_circuit_case(seed: u64, levels: &[crate::optimize::Level]) {
    use crate::optimize::{Options, optimize};
    use crate::semantics::channel::{ChannelLimits, circuit_channel, pbc_channel};
    let (c, initial) = random_measured_circuit(seed);
    let limits = ChannelLimits::default();
    let replay = format!(
        "seed={seed:#x}, initial={initial:?}\n\
         Replay: PBC_FUZZ_SEED={seed} PBC_FUZZ_CASES=1 cargo test --release \
         fuzz_mid_circuit_measurements -- --ignored --nocapture\n{c}"
    );
    let expected = circuit_channel(&c, &initial, limits).unwrap();
    let converted = to_pbc(&c).unwrap_or_else(|e| panic!("conversion failed: {e}, {replay}"));
    assert_linear_size(&c, &converted);
    let actual = pbc_channel(&converted, &initial, limits).unwrap();
    assert_eq!(
        expected.compare(&actual),
        Ok(()),
        "mid-circuit mismatch: {replay}"
    );
    for &level in levels {
        for parallel in [false, true] {
            let options = Options {
                level,
                parallel,
                ..Options::default()
            };
            let (optimized, _) = optimize(&c, &options).unwrap();
            let actual = pbc_channel(&to_pbc(&optimized).unwrap(), &initial, limits).unwrap();
            assert_eq!(
                expected.compare(&actual),
                Ok(()),
                "optimize {level:?} (parallel={parallel}) then convert mismatch: {replay}\n--> {optimized}"
            );
        }
    }
}

/// Widths 1-3 only: exact four-qubit channels take seconds each in debug
/// builds, so they run in the extended fuzzer below.
#[test]
fn bounded_fuzz_mid_circuit_measurements_covers_24_seeded_circuits() {
    let seeds = (0x4d43_4d00..).filter(|seed| seed % 4 != 3);
    for seed in seeds.take(24) {
        mid_circuit_case(seed, &[]);
    }
}

/// Extended runs, also optimizing at every level (sequential and parallel)
/// before conversion. Single-seed replay:
/// PBC_FUZZ_CASES=200 cargo test --release fuzz_mid_circuit_measurements -- --ignored --nocapture
#[test]
#[ignore = "extended mid-circuit fuzzing through the optimizer; set PBC_FUZZ_CASES and PBC_FUZZ_SEED"]
fn fuzz_mid_circuit_measurements() {
    use crate::optimize::Level;
    let (count, seed) = fuzz_range(200, 0x4d43_4e00);
    let levels = [Level::O1, Level::O2, Level::O3, Level::Osuper];
    for case in 0..count {
        mid_circuit_case(seed.wrapping_add(case), &levels);
    }
    eprintln!("PBC mid-circuit fuzz: {count} cases passed, starting seed {seed:#x}");
}
