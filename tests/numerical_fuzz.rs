//! Numerical stress testing with exact basis-state phases, without evaluating
//! sin/cos or using the optimizer's angle arithmetic as an equivalence oracle.
//!
//! Replay any case reported by a failure:
//! NUMERICAL_FUZZ_SEED=123 NUMERICAL_FUZZ_START=17 NUMERICAL_FUZZ_CASES=1 \
//!   cargo test --release --test numerical_fuzz stress -- --ignored --nocapture
//!
//! This oracle covers monomial circuits (X, CNOT, CZ, Toffoli, CCZ and phases).
//! It exhausts computational-basis inputs and compares their exact relative
//! phases as well as the output permutation, proving equivalence on arbitrary
//! superpositions up to global phase. H, measurement and reset are outside this
//! oracle's domain; their channel tests live in semantics::channel.

use std::{collections::BTreeMap, env};

use num_bigint::BigInt;
use num_rational::BigRational;
use tzap::{
    angle::Angle,
    cancel::CancelGates,
    circuit::{Circuit, Gate},
    cnot_min::CnotMin,
    optimize::{Options, PassName, optimize},
    pass::Pass,
    phase_fold_pauli::PhaseFoldPauli,
    phase_fold_rand::PhaseFoldRand,
    qasm,
};

const SEED: u64 = 0x24b0_1695_d7a8_3ef1;

fn zero() -> BigRational {
    BigRational::from_integer(0.into())
}

/// Decode IEEE bits independently: no arithmetic on the input floats.
fn binary(value: f64) -> BigRational {
    assert!(value.is_finite());
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mantissa = bits & ((1 << 52) - 1);
    let significand = if exponent == 0 {
        mantissa
    } else {
        mantissa | (1 << 52)
    };
    let power = if exponent == 0 {
        -1074
    } else {
        exponent - 1075
    };
    let signed = BigInt::from(significand) * if bits >> 63 == 0 { 1 } else { -1 };
    if power >= 0 {
        BigRational::from_integer(signed << power as usize)
    } else {
        BigRational::new(signed, BigInt::from(1) << (-power) as usize)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Phase {
    pi: BigRational,
    radians: BigRational,
    // Preserved expressions are independent formal parameters. Comparing
    // their coefficients proves a stronger property than assuming a value.
    opaque: BTreeMap<String, BigInt>,
}

impl Phase {
    fn zero() -> Self {
        Self {
            pi: zero(),
            radians: zero(),
            opaque: BTreeMap::new(),
        }
    }

    fn add(&mut self, rhs: &Self) {
        self.pi += &rhs.pi;
        self.radians += &rhs.radians;
        for (source, coefficient) in &rhs.opaque {
            *self.opaque.entry(source.clone()).or_default() += coefficient;
        }
    }

    fn relative_to(&self, origin: &Self) -> Self {
        let mut result = self.clone();
        result.pi -= &origin.pi;
        result.radians -= &origin.radians;
        for (source, coefficient) in &origin.opaque {
            *result.opaque.entry(source.clone()).or_default() -= coefficient;
        }
        result.opaque.retain(|_, n| *n != BigInt::from(0));
        result
    }

    fn equivalent(&self, other: &Self) -> bool {
        let delta = self.relative_to(other);
        delta.radians == zero()
            && delta.opaque.is_empty()
            && delta.pi.numer() % (delta.pi.denom() * 2) == BigInt::from(0)
    }
}

fn phase(gate: &Gate) -> Option<(u32, Phase)> {
    let mut p = Phase::zero();
    let (q, quarters) = match gate {
        Gate::t(q) => (*q, 1),
        Gate::tdg(q) => (*q, -1),
        Gate::s(q) => (*q, 2),
        Gate::sdg(q) => (*q, -2),
        Gate::z(q) => (*q, 4),
        Gate::rz(angle, q) => {
            if let Some((pi, residual)) = angle.components() {
                p.pi = BigRational::new(pi.numerator().into(), pi.denominator().into());
                p.radians = binary(residual.get());
            } else {
                p.opaque.insert(angle.to_string(), 1.into());
            }
            return Some((*q, p));
        }
        _ => return None,
    };
    p.pi = BigRational::new(quarters.into(), 4.into());
    Some((q, p))
}

fn semantics(circuit: &Circuit) -> Vec<(usize, Phase)> {
    assert_eq!(circuit.num_cbits, 0);
    assert!(circuit.num_qubits <= 4);
    // Decode each phase once; the state interpreter only uses exact integers.
    let phases: Vec<_> = circuit.gates.iter().map(phase).collect();
    (0..1usize << circuit.num_qubits)
        .map(|input| {
            let mut bits = input;
            let mut accumulated = Phase::zero();
            for (gate, p) in circuit.gates.iter().zip(&phases) {
                let bit = |q: u32| {
                    assert!((q as usize) < circuit.num_qubits);
                    bits & (1 << q) != 0
                };
                if let Some((q, p)) = p {
                    if bit(*q) {
                        accumulated.add(p);
                    }
                    continue;
                }
                match gate {
                    Gate::x(q) => {
                        bit(*q);
                        bits ^= 1 << q;
                    }
                    Gate::cnot { control, target } => {
                        assert_ne!(control, target);
                        let active = bit(*control);
                        bit(*target);
                        if active {
                            bits ^= 1 << target;
                        }
                    }
                    Gate::cz { control, target } => {
                        assert_ne!(control, target);
                        if bit(*control) && bit(*target) {
                            accumulated.pi += BigRational::from_integer(1.into());
                        }
                    }
                    Gate::ccx {
                        control1,
                        control2,
                        target,
                    } => {
                        assert!(control1 != control2 && control1 != target && control2 != target);
                        let active = bit(*control1) && bit(*control2);
                        bit(*target);
                        if active {
                            bits ^= 1 << target;
                        }
                    }
                    Gate::ccz {
                        control1,
                        control2,
                        target,
                    } => {
                        assert!(control1 != control2 && control1 != target && control2 != target);
                        if bit(*control1) && bit(*control2) && bit(*target) {
                            accumulated.pi += BigRational::from_integer(1.into());
                        }
                    }
                    unsupported => panic!("oracle domain violated by {unsupported:?}"),
                }
            }
            (bits, accumulated)
        })
        .collect()
}

fn equivalent(expected: &[(usize, Phase)], actual: &Circuit) -> bool {
    let actual = semantics(actual);
    expected.len() == actual.len()
        && expected.iter().zip(&actual).all(|((a, p), (b, q))| {
            a == b
                && p.relative_to(&expected[0].1)
                    .equivalent(&q.relative_to(&actual[0].1))
        })
}

struct Random(u64);
impl Random {
    fn next(&mut self) -> u64 {
        // SplitMix64: a case's seed determines its input independently of the
        // number of random draws made in earlier cases, including seed zero.
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut n = self.0;
        n = (n ^ (n >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        n = (n ^ (n >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        n ^ (n >> 31)
    }
    fn index(&mut self, count: usize) -> usize {
        self.next() as usize % count
    }
    fn float(&mut self) -> f64 {
        const EDGES: &[f64] = &[
            0.0,
            -0.0,
            f64::from_bits(1),
            -f64::from_bits(1),
            f64::from_bits((1 << 52) - 1),
            f64::MIN_POSITIVE,
            f64::MAX,
            -f64::MAX,
            1e308,
            -1e308,
            1e16,
            -1e16,
            1.0,
            -1.0,
            0.1,
            0.2,
            0.3,
            0.7,
            5e-7,
            -5e-7,
            1.570796,
            -1.570796,
            std::f64::consts::PI / 4.0,
        ];
        match self.index(4) {
            0 | 1 => EDGES[self.index(EDGES.len())],
            2 => {
                let center = (std::f64::consts::PI / 4.0).to_bits();
                f64::from_bits(center + self.index(5) as u64 - 2)
            }
            _ => loop {
                let value = f64::from_bits(self.next());
                if value.is_finite() {
                    break value;
                }
            },
        }
    }
    fn angle(&mut self) -> Angle {
        match self.index(5) {
            0 | 1 => Angle::from_f64(self.float()).unwrap(),
            2 => {
                let n = [i64::MIN, i64::MAX, -7, -1, 0, 1, 3, 8][self.index(8)];
                let d = [1, 4, 7, 8, i64::MAX, i64::MAX - 1][self.index(6)];
                Angle::pi_fraction(n, d).unwrap()
            }
            _ => {
                let expression = if self.index(4) == 0 {
                    // Both a coefficient-limit and a rounded affine residual
                    // must remain formal, unconverted expressions.
                    ["pi/9223372036854775808", "pi/4+0.1+0.2"][self.index(2)].to_string()
                } else {
                    format!("pi/{}+({})", [4, 7, 8][self.index(3)], self.float())
                };
                let mut parsed =
                    qasm::parse(&format!("qreg q[1]; rz({expression}) q[0];")).unwrap();
                match parsed.gates.pop().unwrap() {
                    Gate::rz(a, _) => a,
                    _ => unreachable!(),
                }
            }
        }
    }
}

fn generate(seed: u64, case: usize) -> Circuit {
    let mut r = Random(seed ^ (case as u64).wrapping_mul(0xd134_2543_de82_ef95));
    let n = 1 + r.index(4);
    let mut circuit = Circuit::new(n);
    let numeric = |v| Gate::rz_f64(v, 0).unwrap();
    // Mandatory adversarial prefixes force repeated failed merges rather than
    // relying on uniformly random exponent bits to create difficult sums.
    circuit.gates = match case % 6 {
        0 => vec![numeric(1e16), numeric(1.0), numeric(-1e16)],
        1 => vec![numeric(1e308), numeric(1e308), numeric(-1e308)],
        2 => vec![
            Gate::rz(Angle::pi_fraction(1, i64::MAX).unwrap(), 0),
            Gate::rz(Angle::pi_fraction(1, i64::MAX - 1).unwrap(), 0),
        ],
        3 => vec![numeric(5e-7), numeric(std::f64::consts::PI / 4.0)],
        4 => vec![
            numeric(1e16),
            Gate::x(0),
            numeric(1.0),
            Gate::x(0),
            numeric(-1e16),
        ],
        _ => vec![
            numeric(f64::MAX),
            numeric(-f64::MAX),
            numeric(f64::from_bits(1)),
        ],
    };
    for _ in 0..12 + r.index(24) {
        let q = r.index(n) as u32;
        let mut operands: Vec<_> = (0..n as u32).collect();
        for i in (1..n).rev() {
            let j = r.index(i + 1);
            operands.swap(i, j);
        }
        match r.index(12) {
            0 => circuit.apply(Gate::x(q)),
            1 | 2 if n > 1 => circuit.apply(Gate::cnot {
                control: operands[0],
                target: operands[1],
            }),
            3 if n > 1 => circuit.apply(Gate::cz {
                control: operands[0],
                target: operands[1],
            }),
            4 if n > 2 => circuit.apply(Gate::ccx {
                control1: operands[0],
                control2: operands[1],
                target: operands[2],
            }),
            5 if n > 2 => circuit.apply(Gate::ccz {
                control1: operands[0],
                control2: operands[1],
                target: operands[2],
            }),
            6 => circuit.apply(Gate::t(q)),
            7 => circuit.apply(Gate::tdg(q)),
            8 => {
                let a = r.angle();
                circuit.gates.extend([
                    Gate::rz(a.clone(), q),
                    Gate::x(q),
                    Gate::rz(a, q),
                    Gate::x(q),
                ]);
            }
            _ => circuit.apply(Gate::rz(r.angle(), q)),
        }
    }
    circuit
}

#[derive(Default)]
struct Coverage {
    rounded: usize,
    nonfinite: usize,
    coefficient: usize,
    preserved: usize,
    reduced: usize,
}

fn run_cases(seed: u64, start: usize, cases: usize) -> Coverage {
    let mut coverage = Coverage::default();
    for case in start..start.checked_add(cases).expect("case range overflow") {
        // Annotate production panics and oracle-domain failures too, rather
        // than only failed equivalence assertions, with a replayable case.
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let input = generate(seed, case);
            let expected = semantics(&input);
            let check = |output: &Circuit, context: &str| {
                assert_eq!(output.num_qubits, input.num_qubits);
                assert!(
                    equivalent(&expected, output),
                    "seed={seed} case={case} profile={context}\ninput:\n{}\noutput:\n{}",
                    input.to_qasm(),
                    output.to_qasm()
                );
                let reparsed = qasm::parse(&output.to_qasm()).unwrap_or_else(|e| {
                    panic!(
                        "seed={seed} case={case} profile={context}: output does not parse: {e}\n{}",
                        output.to_qasm()
                    )
                });
                assert!(
                    equivalent(&expected, &reparsed),
                    "seed={seed} case={case} profile={context}: serialization changed semantics\ninput:\n{}\noutput:\n{}",
                    input.to_qasm(),
                    output.to_qasm()
                );
            };
            check(&input, "input round trip");
            let cnot = CnotMin::default();
            let passes: [&dyn Pass; 4] = [&CancelGates, &PhaseFoldRand, &PhaseFoldPauli, &cnot];
            for pass in passes {
                check(&pass.run(&input), pass.name());
            }
            for parallel in [false, true] {
                let options = Options {
                    passes: Some(vec![
                        PassName::PhaseFoldRand,
                        PassName::CnotMin,
                        PassName::PhaseFoldPauli,
                        PassName::CancelGates,
                    ]),
                    fixpoint: true,
                    parallel,
                    ..Options::default()
                };
                let (output, report) = optimize(&input, &options).unwrap();
                check(
                    &output,
                    if parallel {
                        "parallel fixpoint"
                    } else {
                        "serial fixpoint"
                    },
                );
                coverage.rounded += report.numerical.skipped_rounded_folds;
                coverage.nonfinite += report.numerical.skipped_nonfinite_folds;
                coverage.coefficient += report.numerical.skipped_coefficient_limit_folds;
                coverage.preserved += report.numerical.preserved_expressions;
                coverage.reduced += usize::from(output.gates.len() < input.gates.len());
            }
            // A different composition and repeated application exercise changed
            // axes and replacement history after a previous pass's successful fold.
            let mut output = input.clone();
            for _ in 0..2 {
                for pass in passes.into_iter().rev() {
                    output = pass.run(&output);
                }
            }
            check(&output, "reverse composition twice");
        }));
        if let Err(reason) = outcome {
            let message = reason
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| reason.downcast_ref::<&str>().copied())
                .unwrap_or("non-string panic");
            panic!("numerical fuzz failed: seed={seed} case={case}\n{message}");
        }
    }
    eprintln!(
        "numerical fuzz seed={seed} start={start} cases={cases}: rounded={} nonfinite={} coefficient={} preserved={} reduced={}",
        coverage.rounded,
        coverage.nonfinite,
        coverage.coefficient,
        coverage.preserved,
        coverage.reduced
    );
    coverage
}

#[test]
fn exact_oracle_detects_numeric_errors_and_ignores_only_global_phase() {
    let mut identity = Circuit::new(1);
    let expected = semantics(&identity);
    identity.apply(Gate::rz_f64(f64::from_bits(1), 0).unwrap());
    assert!(!equivalent(&expected, &identity)); // A floating norm would miss this.
    let mut large = Circuit::new(1);
    for a in [1e16, 1.0, -1e16] {
        large.apply(Gate::rz_f64(a, 0).unwrap());
    }
    assert!(!equivalent(&expected, &large));
    let mut quarter = Circuit::new(1);
    quarter.apply(Gate::rz_f64(std::f64::consts::PI / 4.0, 0).unwrap());
    let mut t = Circuit::new(1);
    t.apply(Gate::t(0));
    assert!(!equivalent(&semantics(&quarter), &t));
    let mut global = Circuit::new(1);
    global.gates = vec![
        Gate::x(0),
        Gate::rz_f64(0.3, 0).unwrap(),
        Gate::x(0),
        Gate::rz_f64(0.3, 0).unwrap(),
    ];
    assert!(equivalent(&expected, &global));
    let mut two_qubits = Circuit::new(2);
    let expected = semantics(&two_qubits);
    two_qubits.apply(Gate::cnot {
        control: 0,
        target: 1,
    });
    assert!(!equivalent(&expected, &two_qubits));
}

#[test]
fn bounded_numerical_fuzz_covers_exact_arithmetic_and_rollback() {
    let coverage = run_cases(SEED, 0, 32);
    assert!(
        coverage.rounded > 0
            && coverage.nonfinite > 0
            && coverage.coefficient > 0
            && coverage.preserved > 0
            && coverage.reduced > 0,
        "the bounded sample must exercise successful rewrites and every failure category"
    );
}

fn setting(name: &str, default: u64) -> u64 {
    match env::var(name) {
        Ok(value) => value
            .parse()
            .unwrap_or_else(|_| panic!("{name} must be an unsigned decimal integer")),
        Err(env::VarError::NotPresent) => default,
        Err(e) => panic!("invalid {name}: {e}"),
    }
}

#[test]
#[ignore = "numerical stress fuzzer; run in release with --ignored --nocapture"]
fn numerical_stress_fuzz() {
    let seed = setting("NUMERICAL_FUZZ_SEED", SEED);
    let start = usize::try_from(setting("NUMERICAL_FUZZ_START", 0)).unwrap();
    let cases = usize::try_from(setting("NUMERICAL_FUZZ_CASES", 2048)).unwrap();
    assert!(cases > 0, "NUMERICAL_FUZZ_CASES must be positive");
    run_cases(seed, start, cases);
}
