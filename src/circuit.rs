//! Circuit representation: gates, qubits, and display.

use std::fmt;

pub type Qubit = usize;

#[allow(non_camel_case_types)]
#[derive(Clone, Debug)]
pub enum Gate {
    x(Qubit),
    h(Qubit),
    s(Qubit),
    sdg(Qubit),
    z(Qubit),
    t(Qubit),
    tdg(Qubit),
    rz(f64, Qubit),
    cnot { control: Qubit, target: Qubit },
    ccx { control1: Qubit, control2: Qubit, target: Qubit },
}

#[derive(Clone, Debug)]
pub struct Circuit {
    pub num_qubits: usize,
    pub gates: Vec<Gate>,
    pub has_toffoli: bool,
}

impl Circuit {
    pub fn new(num_qubits: usize) -> Self {
        Circuit { num_qubits, gates: Vec::new(), has_toffoli: false }
    }

    pub fn apply(&mut self, gate: Gate) {
        if matches!(gate, Gate::ccx { .. }) {
            self.has_toffoli = true;
        }
        self.gates.push(gate);
    }

    pub fn to_qasm(&self) -> String {
        crate::qasm::serialize(self)
    }

    pub fn from_qasm(qasm: &str) -> Result<Self, String> {
        crate::qasm::parse(qasm)
    }
}

impl fmt::Display for Gate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Gate::x(q) => write!(f, "x q{q}"),
            Gate::h(q) => write!(f, "h q{q}"),
            Gate::s(q) => write!(f, "s q{q}"),
            Gate::sdg(q) => write!(f, "sdg q{q}"),
            Gate::z(q) => write!(f, "z q{q}"),
            Gate::t(q) => write!(f, "t q{q}"),
            Gate::tdg(q) => write!(f, "tdg q{q}"),
            Gate::rz(theta, q) => write!(f, "rz({theta:.4}) q{q}"),
            Gate::cnot { control, target } => write!(f, "cnot q{control}, q{target}"),
            Gate::ccx { control1, control2, target } => {
                write!(f, "ccx q{control1}, q{control2}, q{target}")
            }
        }
    }
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Circuit ({} qubits, {} gates):", self.num_qubits, self.gates.len())?;
        for (i, gate) in self.gates.iter().enumerate() {
            writeln!(f, "  {i}: {gate}")?;
        }
        Ok(())
    }
}

/// Return the qubits a gate acts on.
pub fn qubits_of(gate: &Gate) -> Vec<Qubit> {
    match gate {
        Gate::x(q) | Gate::h(q) | Gate::s(q) | Gate::sdg(q) | Gate::z(q)
        | Gate::t(q) | Gate::tdg(q) | Gate::rz(_, q) => vec![*q],
        Gate::cnot { control, target } => vec![*control, *target],
        Gate::ccx { control1, control2, target } => vec![*control1, *control2, *target],
    }
}

/// Remap a gate's qubits through a lookup table: qubit i becomes its index in `qubits`.
pub fn remap_gate(gate: &Gate, qubits: &[Qubit]) -> Gate {
    let m = |q: &Qubit| qubits.iter().position(|&x| x == *q).unwrap();
    match gate {
        Gate::x(q) => Gate::x(m(q)),
        Gate::h(q) => Gate::h(m(q)),
        Gate::s(q) => Gate::s(m(q)),
        Gate::sdg(q) => Gate::sdg(m(q)),
        Gate::z(q) => Gate::z(m(q)),
        Gate::t(q) => Gate::t(m(q)),
        Gate::tdg(q) => Gate::tdg(m(q)),
        Gate::rz(theta, q) => Gate::rz(*theta, m(q)),
        Gate::cnot { control, target } => Gate::cnot { control: m(control), target: m(target) },
        Gate::ccx { control1, control2, target } => Gate::ccx {
            control1: m(control1),
            control2: m(control2),
            target: m(target),
        },
    }
}

/// Build a compact circuit with qubits remapped to 0..n.
pub fn remap_subcircuit(gates: &[Gate], qubits: &[Qubit]) -> Circuit {
    let n = qubits.len();
    let mut c = Circuit::new(n);
    for g in gates {
        c.apply(remap_gate(g, qubits));
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn bell_pair() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        assert_eq!(c.gates.len(), 2);
        let s = format!("{c}");
        assert!(s.contains("h q0"));
        assert!(s.contains("cnot q0, q1"));
        println!("{c}");
    }

    #[test]
    fn ghz_state() {
        let n = 4;
        let mut c = Circuit::new(n);
        c.apply(Gate::h(0));
        for i in 0..n - 1 {
            c.apply(Gate::cnot { control: i, target: i + 1 });
        }
        assert_eq!(c.gates.len(), 4);
        println!("{c}");
    }

    #[test]
    fn t_gate_decomposition_of_rz() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::s(0));
        c.apply(Gate::rz(PI / 4.0, 0));
        let s = format!("{c}");
        assert!(s.contains("t q0"));
        assert!(s.contains("s q0"));
        assert!(s.contains("rz(0.7854) q0"));
        println!("{c}");
    }

    #[test]
    fn ccx_gate() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(2));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::h(2));
        assert_eq!(c.gates.len(), 3);
        let s = format!("{c}");
        assert!(s.contains("ccx q0, q1, q2"));
        println!("{c}");
    }

    #[test]
    fn qft_3qubit() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(0));
        c.apply(Gate::rz(PI / 2.0, 0));
        c.apply(Gate::cnot { control: 1, target: 0 });
        c.apply(Gate::rz(PI / 4.0, 0));
        c.apply(Gate::cnot { control: 2, target: 0 });
        c.apply(Gate::h(1));
        c.apply(Gate::rz(PI / 2.0, 1));
        c.apply(Gate::cnot { control: 2, target: 1 });
        c.apply(Gate::h(2));
        assert_eq!(c.num_qubits, 3);
        assert_eq!(c.gates.len(), 9);
        println!("{c}");
    }

    #[test]
    fn z_gate_display() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        let s = format!("{c}");
        assert!(s.contains("z q0"));
    }

    #[test]
    fn sdg_gate_display() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        let s = format!("{c}");
        assert!(s.contains("sdg q0"));
    }
}
