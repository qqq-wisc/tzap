//! Experimental: quantum WHILE programs, for phase folding across control flow.
//!
//! A [`Program`] is a statement tree over [`Gate`]s, resets and measurements,
//! with non-deterministic branches and loops (`if ★`, `while ★`), as in Amy and
//! Lunderville, *Linear and non-linear relational analyses for quantum program
//! optimization* (POPL 2025). [`parse`] reads the OpenQASM 3 subset their
//! benchmarks use: branch and loop conditions are ignored (treated as ★),
//! `for` loops over constant ranges are unrolled, gate definitions are inlined,
//! and Toffoli-family gates are expanded into Clifford+T.

mod parse;
pub mod pauli_fold;
pub mod state_fold;

use std::fmt::Write as _;

use crate::circuit::{Gate, Qubit};

pub use parse::{parse, ParseError};

/// A program statement.
#[derive(Clone, Debug, PartialEq)]
pub enum Stmt {
    Gate(Gate),
    Reset(Qubit),
    Measure(Qubit),
    Seq(Vec<Stmt>),
    /// A non-deterministic branch: either side may run.
    If(Box<Stmt>, Box<Stmt>),
    /// A non-deterministic loop: the body runs any number of times.
    While(Box<Stmt>),
}

/// A program over `num_qubits` qubits, named for output.
#[derive(Clone, Debug)]
pub struct Program {
    pub num_qubits: usize,
    pub body: Stmt,
}

impl Stmt {
    /// Every gate, in program order, with nesting flattened.
    pub fn gates(&self) -> Vec<&Gate> {
        let mut out = Vec::new();
        self.visit_gates(&mut |g| out.push(g));
        out
    }

    fn visit_gates<'a>(&'a self, f: &mut impl FnMut(&'a Gate)) {
        match self {
            Stmt::Gate(g) => f(g),
            Stmt::Reset(_) | Stmt::Measure(_) => {}
            Stmt::Seq(xs) => xs.iter().for_each(|x| x.visit_gates(f)),
            Stmt::If(a, b) => {
                a.visit_gates(f);
                b.visit_gates(f);
            }
            Stmt::While(b) => b.visit_gates(f),
        }
    }
}

/// Whether `theta` is a multiple of pi/2 (a Clifford rotation).
fn clifford_angle(theta: f64) -> bool {
    let k = theta / std::f64::consts::FRAC_PI_2;
    (k - k.round()).abs() < 1e-9
}

impl Program {
    /// The static count of non-Clifford rotations: T, T-dagger, and Rz by a
    /// non-Clifford angle, each counted once however often it may run.
    pub fn t_count(&self) -> usize {
        self.body
            .gates()
            .iter()
            .filter(|g| match g {
                Gate::t(_) | Gate::tdg(_) => true,
                Gate::rz(theta, _) => !clifford_angle(*theta),
                _ => false,
            })
            .count()
    }

    /// The program as OpenQASM 3 over one register `q`.
    pub fn to_qasm3(&self) -> String {
        let mut out = String::from("OPENQASM 3.0;\ninclude \"stdgates.inc\";\n");
        let _ = writeln!(out, "qubit[{}] q;", self.num_qubits);
        write_stmt(&mut out, &self.body, 0);
        out
    }
}

fn write_stmt(out: &mut String, s: &Stmt, depth: usize) {
    let pad = "  ".repeat(depth);
    match s {
        Stmt::Gate(g) => {
            let _ = writeln!(out, "{pad}{};", gate_text(g));
        }
        Stmt::Reset(q) => {
            let _ = writeln!(out, "{pad}reset q[{q}];");
        }
        Stmt::Measure(q) => {
            let _ = writeln!(out, "{pad}measure q[{q}];");
        }
        Stmt::Seq(xs) => xs.iter().for_each(|x| write_stmt(out, x, depth)),
        Stmt::If(a, b) => {
            let _ = writeln!(out, "{pad}if (true) {{");
            write_stmt(out, a, depth + 1);
            let _ = writeln!(out, "{pad}}} else {{");
            write_stmt(out, b, depth + 1);
            let _ = writeln!(out, "{pad}}}");
        }
        Stmt::While(b) => {
            let _ = writeln!(out, "{pad}while (true) {{");
            write_stmt(out, b, depth + 1);
            let _ = writeln!(out, "{pad}}}");
        }
    }
}

fn gate_text(g: &Gate) -> String {
    match *g {
        Gate::h(q) => format!("h q[{q}]"),
        Gate::x(q) => format!("x q[{q}]"),
        Gate::z(q) => format!("z q[{q}]"),
        Gate::s(q) => format!("s q[{q}]"),
        Gate::sdg(q) => format!("sdg q[{q}]"),
        Gate::t(q) => format!("t q[{q}]"),
        Gate::tdg(q) => format!("tdg q[{q}]"),
        Gate::rz(theta, q) => format!("rz({theta}) q[{q}]"),
        Gate::cnot { control, target } => format!("cx q[{control}], q[{target}]"),
        Gate::cz { control, target } => format!("cz q[{control}], q[{target}]"),
        ref other => format!("{other:?}"),
    }
}

#[cfg(test)]
pub(crate) mod testing {
    //! Path semantics for checking program rewrites: every path (loops
    //! unrolled up to a bound, both sides of each branch, each measurement and
    //! reset outcome) must keep its operator, up to a phase of its own.

    use super::Stmt;
    use crate::circuit::{Circuit, Gate};
    use crate::unitary::{C, circuit_unitary};

    pub type Mat = Vec<Vec<C>>;

    #[derive(Clone)]
    pub enum Action {
        Gate(Gate),
        /// Post-select qubit on outcome.
        Assume(u32, bool),
    }

    /// The paths of `s`, with loops run at most `bound` times.
    pub fn paths(s: &Stmt, bound: usize) -> Vec<Vec<Action>> {
        match s {
            Stmt::Gate(g) => vec![vec![Action::Gate(g.clone())]],
            Stmt::Measure(q) => vec![vec![Action::Assume(*q, false)], vec![Action::Assume(*q, true)]],
            Stmt::Reset(q) => vec![
                vec![Action::Assume(*q, false)],
                vec![Action::Assume(*q, true), Action::Gate(Gate::x(*q))],
            ],
            Stmt::Seq(xs) => {
                let mut acc = vec![Vec::new()];
                for x in xs {
                    let ps = paths(x, bound);
                    acc = acc
                        .iter()
                        .flat_map(|a| ps.iter().map(move |p| [a.clone(), p.clone()].concat()))
                        .collect();
                }
                acc
            }
            Stmt::If(a, b) => [paths(a, bound), paths(b, bound)].concat(),
            Stmt::While(b) => {
                let body = paths(b, bound);
                let mut out = vec![Vec::new()];
                let mut layer = vec![Vec::new()];
                for _ in 0..bound {
                    layer = layer
                        .iter()
                        .flat_map(|a: &Vec<Action>| body.iter().map(move |p| [a.clone(), p.clone()].concat()))
                        .collect();
                    out.extend(layer.iter().cloned());
                }
                out
            }
        }
    }

    /// The number of paths, without building them.
    pub fn count(s: &Stmt, bound: usize) -> usize {
        match s {
            Stmt::Gate(_) => 1,
            Stmt::Measure(_) | Stmt::Reset(_) => 2,
            Stmt::Seq(xs) => xs.iter().map(|x| count(x, bound)).fold(1usize, |a, b| a.saturating_mul(b)),
            Stmt::If(a, b) => count(a, bound) + count(b, bound),
            Stmt::While(b) => {
                let c = count(b, bound);
                (0..=bound as u32).map(|k| c.saturating_pow(k)).fold(0usize, |a, b| a.saturating_add(b))
            }
        }
    }

    fn mul(a: &Mat, b: &Mat) -> Mat {
        let n = a.len();
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| (0..n).fold(C::ZERO, |acc, k| acc + a[i][k] * b[k][j]))
                    .collect()
            })
            .collect()
    }

    fn gate_matrix(n: usize, g: &Gate) -> Mat {
        let mut c = Circuit::with_cbits(n, 0);
        c.apply(g.clone());
        circuit_unitary(&c)
    }

    /// The operator of a path.
    pub fn operator(n: usize, path: &[Action]) -> Mat {
        let dim = 1 << n;
        let mut m: Mat = (0..dim)
            .map(|i| (0..dim).map(|j| if i == j { C::ONE } else { C::ZERO }).collect())
            .collect();
        for a in path {
            let step = match a {
                Action::Gate(g) => gate_matrix(n, g),
                Action::Assume(q, b) => {
                    // (I ± Z_q) / 2.
                    let z = gate_matrix(n, &Gate::z(*q));
                    let sign = if *b { -1.0 } else { 1.0 };
                    (0..dim)
                        .map(|i| {
                            (0..dim)
                                .map(|j| {
                                    let id = if i == j { C::ONE } else { C::ZERO };
                                    (id + z[i][j] * sign) * 0.5
                                })
                                .collect()
                        })
                        .collect()
                }
            };
            m = mul(&step, &m);
        }
        m
    }

    /// Equal up to a phase (both may be zero).
    pub fn equiv(a: &Mat, b: &Mat, tol: f64) -> bool {
        let n = a.len();
        let mut phase: Option<C> = None;
        for i in 0..n {
            for j in 0..n {
                let (x, y) = (a[i][j], b[i][j]);
                if x.norm_sq() > tol && phase.is_none() {
                    let r = x.conj() * y * (1.0 / x.norm_sq());
                    phase = Some(r);
                }
            }
        }
        let Some(p) = phase else {
            return b.iter().flatten().all(|y| y.norm_sq() < tol);
        };
        if (p.norm_sq() - 1.0).abs() > 1e-6 {
            return false;
        }
        (0..n).all(|i| (0..n).all(|j| {
            let d = a[i][j] * p + b[i][j] * -1.0;
            d.norm_sq() < tol
        }))
    }

    /// A random program over `n` qubits.
    pub fn random_stmt(rng: &mut u64, n: u32, len: usize, depth: usize) -> Stmt {
        let mut next = |m: u64| -> u64 {
            *rng ^= *rng << 13;
            *rng ^= *rng >> 7;
            *rng ^= *rng << 17;
            *rng % m
        };
        let mut out = Vec::new();
        for _ in 0..len {
            let q = next(n as u64) as u32;
            let r = (q + 1 + next(n as u64 - 1) as u32) % n;
            let k = next(if depth > 0 { 18 } else { 15 });
            out.push(match k {
                0..=2 => Stmt::Gate(Gate::h(q)),
                3..=5 => Stmt::Gate(Gate::t(q)),
                6 => Stmt::Gate(Gate::tdg(q)),
                7 => Stmt::Gate(Gate::s(q)),
                8 => Stmt::Gate(Gate::x(q)),
                9 | 10 => Stmt::Gate(Gate::cnot { control: q, target: r }),
                11 => Stmt::Gate(Gate::rz(0.3, q)),
                12 | 13 => Stmt::Measure(q),
                14 => Stmt::Reset(q),
                15 => {
                    let seed = next(1 << 30);
                    let mut s2 = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
                    let a = random_stmt(&mut s2, n, 1 + (seed % 3) as usize, depth - 1);
                    let b = random_stmt(&mut s2, n, (seed % 2) as usize, depth - 1);
                    Stmt::If(Box::new(a), Box::new(b))
                }
                _ => {
                    let seed = next(1 << 30);
                    let mut s2 = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
                    Stmt::While(Box::new(random_stmt(&mut s2, n, 1 + (seed % 3) as usize, depth - 1)))
                }
            });
        }
        Stmt::Seq(out)
    }
}
