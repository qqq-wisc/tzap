use super::*;
use rand::{Rng, SeedableRng, rngs::StdRng, seq::SliceRandom};

fn random_unitary(n: usize, rng: &mut StdRng) -> Gate {
    let mut q: Vec<_> = (0..n as u32).collect();
    q.shuffle(rng);
    match rng.gen_range(0..11) {
        0 => Gate::h(q[0]),
        1 => Gate::x(q[0]),
        2 => Gate::z(q[0]),
        3 => Gate::s(q[0]),
        4 => Gate::sdg(q[0]),
        5 => Gate::t(q[0]),
        6 => Gate::tdg(q[0]),
        7 if n >= 2 => cx(q[0], q[1]),
        8 if n >= 2 => Gate::cz {
            control: q[0],
            target: q[1],
        },
        9 if n >= 3 => Gate::ccx {
            control1: q[0],
            control2: q[1],
            target: q[2],
        },
        10 if n >= 3 => Gate::ccz {
            control1: q[0],
            control2: q[1],
            target: q[2],
        },
        _ => Gate::h(q[0]),
    }
}

#[test]
fn seeded_measurement_channel_fuzz_96_cases() {
    for case in 0..96 {
        let seed = 0x4348_414e_0000 + case;
        let mut rng = StdRng::seed_from_u64(seed);
        let n = 1 + (case % 3) as usize;
        let initial: Vec<bool> = (0..n + 1).map(|_| rng.gen_bool(0.5)).collect();
        let mut gates: Vec<_> = (0..8).map(|_| random_unitary(n, &mut rng)).collect();
        // Full, partial, repeated/overwritten, mid-circuit, or no terminal
        // measurements. The extra classical bit is untouched. No resets.
        let q = rng.gen_range(0..n) as u32;
        gates.extend((0..3).map(|_| random_unitary(n, &mut rng)));
        match case % 4 {
            0 | 2 => {
                if case % 4 == 2 {
                    gates.push(measure(q, 0));
                }
                for q in (0..n as u32).rev() {
                    gates.push(measure(q, (n as u32 - 1) - q));
                }
            }
            1 => gates.push(measure(q, 0)),
            _ => {
                // Mid-circuit readouts, followed by more gates.
                gates.insert(rng.gen_range(0..gates.len()), measure(q, 0));
                let p = rng.gen_range(0..n) as u32;
                gates.insert(rng.gen_range(0..gates.len()), measure(p, n as u32 - 1));
            }
        }
        let input = circuit(n, n + 1, gates);
        let result = std::panic::catch_unwind(|| check(&input, &initial));
        assert!(
            result.is_ok(),
            "channel fuzz failed seed={seed:#x}, case={case}, initial={initial:?}\n{input}"
        );
    }
}
