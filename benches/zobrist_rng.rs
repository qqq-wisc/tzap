use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hasher};
use std::time::Instant;

use rand::Rng;

const N: usize = 10_000_000;
const ROUNDS: usize = 5;

fn bench_random_state() -> u64 {
    let mut sink = 0u64;
    for _ in 0..N {
        sink ^= RandomState::new().build_hasher().finish();
    }
    sink
}

fn bench_rand_thread_rng() -> u64 {
    let mut rng = rand::thread_rng();
    let mut sink = 0u64;
    for _ in 0..N {
        sink ^= rng.r#gen::<u64>();
    }
    sink
}

fn main() {
    println!("Generating {} random u64 values × {} rounds\n", N, ROUNDS);

    for name_fn in [
        ("RandomState::new()", bench_random_state as fn() -> u64),
        ("rand::thread_rng()", bench_rand_thread_rng),
    ] {
        let (name, f) = name_fn;
        let mut best = f64::MAX;
        for _ in 0..ROUNDS {
            let start = Instant::now();
            std::hint::black_box(f());
            let elapsed = start.elapsed().as_secs_f64();
            if elapsed < best { best = elapsed; }
        }
        let ns_per = best / N as f64 * 1e9;
        println!("  {:<25} {:>8.3}s  ({:.1} ns/call)", name, best, ns_per);
    }
}
