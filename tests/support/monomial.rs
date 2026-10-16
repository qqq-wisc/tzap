use num_bigint::BigInt;
use num_rational::BigRational;
use std::collections::BTreeMap;
use tzap::circuit::{Circuit, Gate};

fn zero() -> BigRational {
    BigRational::from_integer(0.into())
}

/// Decode IEEE bits independently: no arithmetic on the input floats.
pub fn binary(value: f64) -> BigRational {
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
pub struct Phase {
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
        Gate::rz(angle, q) | Gate::p(angle, q) => {
            if let Some((pi, residual)) = angle.affine_coefficients() {
                p.pi = pi.to_big_rational();
                p.radians = if let Some(f) = angle.to_f64_exact() {
                    binary(f)
                } else {
                    residual.to_big_rational()
                };
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

/// Independently recognize local basis changes that make RX/RY diagonal.
fn diagonalize(gates: &[Gate]) -> Vec<Gate> {
    let mut current = gates.to_vec();
    loop {
        let mut out = Vec::new();
        let mut i = 0;
        while i < current.len() {
            let rest = &current[i..];
            if let [
                Gate::h(a),
                Gate::s(b),
                Gate::ry(angle, c),
                Gate::sdg(d),
                Gate::h(e),
                ..,
            ] = rest
            {
                if a == b && b == c && c == d && d == e {
                    out.push(Gate::rz(angle.clone(), *a));
                    i += 5;
                    continue;
                }
            }
            if let [Gate::h(a), Gate::rx(angle, b), Gate::h(c), ..] = rest {
                if a == b && b == c {
                    out.push(Gate::rz(angle.clone(), *a));
                    i += 3;
                    continue;
                }
            }
            if let [a, b, ..] = rest {
                let inverse = matches!((a,b), (Gate::h(q), Gate::h(r)) | (Gate::s(q), Gate::sdg(r)) | (Gate::sdg(q), Gate::s(r)) if q == r);
                if inverse {
                    i += 2;
                    continue;
                }
            }
            out.push(current[i].clone());
            i += 1;
        }
        if out.len() == current.len() {
            return out;
        }
        current = out;
    }
}

pub fn semantics(circuit: &Circuit) -> Vec<(usize, Phase)> {
    assert_eq!(circuit.num_cbits, 0);
    assert!(circuit.num_qubits <= 4);
    // Decode each phase once; the state interpreter only uses exact integers.
    let gates = diagonalize(circuit.gates());
    let phases: Vec<_> = gates.iter().map(phase).collect();
    (0..1usize << circuit.num_qubits)
        .map(|input| {
            let mut bits = input;
            let mut accumulated = Phase::zero();
            for (gate, p) in gates.iter().zip(&phases) {
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
                    Gate::y(q) => {
                        if bit(*q) {
                            accumulated.pi += BigRational::from_integer(1.into());
                        }
                        bits ^= 1 << q;
                    }
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
                    Gate::swap(a, b) => {
                        if bit(*a) != bit(*b) {
                            bits ^= (1 << a) | (1 << b);
                        }
                    }
                    Gate::cswap {
                        control,
                        first,
                        second,
                    } => {
                        if bit(*control) && bit(*first) != bit(*second) {
                            bits ^= (1 << first) | (1 << second);
                        }
                    }
                    Gate::cp {
                        lambda,
                        control,
                        target,
                    }
                    | Gate::crz {
                        theta: lambda,
                        control,
                        target,
                    } => {
                        let (pi, residual) =
                            lambda.affine_coefficients().expect("affine oracle domain");
                        let mut angle = Phase::zero();
                        angle.pi = pi.to_big_rational();
                        angle.radians = residual.to_big_rational();
                        if bit(*control) {
                            if matches!(gate, Gate::crz { .. }) {
                                let factor = BigRational::new(
                                    if bit(*target) { 1.into() } else { (-1).into() },
                                    2.into(),
                                );
                                angle.pi *= &factor;
                                angle.radians *= factor;
                                accumulated.add(&angle);
                            } else if bit(*target) {
                                accumulated.add(&angle);
                            }
                        }
                    }
                    unsupported => panic!("oracle domain violated by {unsupported:?}"),
                }
            }
            (bits, accumulated)
        })
        .collect()
}

pub fn equivalent(expected: &[(usize, Phase)], actual: &Circuit) -> bool {
    let actual = semantics(actual);
    expected.len() == actual.len()
        && expected.iter().zip(&actual).all(|((a, p), (b, q))| {
            a == b
                && p.relative_to(&expected[0].1)
                    .equivalent(&q.relative_to(&actual[0].1))
        })
}
