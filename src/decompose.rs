use indicatif::ProgressBar;
use crate::circuit::{Circuit, Gate};
use crate::pass::Pass;

pub struct DecomposeToffoli;

impl Pass for DecomposeToffoli {
    fn name(&self) -> &str { "Toffoli decomposition" }
    fn run_with_progress(&self, circuit: &Circuit, pb: &ProgressBar) -> Circuit {
        let mut output = Circuit::new(circuit.num_qubits);
        for gate in &circuit.gates {
            pb.inc(1);
            match gate {
                Gate::ccx { control1, control2, target } => {
                    let c0 = *control1;
                    let c1 = *control2;
                    let t = *target;
                    output.apply(Gate::h(t));
                    output.apply(Gate::cnot { control: c1, target: t });
                    output.apply(Gate::tdg(t));
                    output.apply(Gate::cnot { control: c0, target: t });
                    output.apply(Gate::t(t));
                    output.apply(Gate::cnot { control: c1, target: t });
                    output.apply(Gate::tdg(t));
                    output.apply(Gate::cnot { control: c0, target: t });
                    output.apply(Gate::t(c1));
                    output.apply(Gate::t(t));
                    output.apply(Gate::h(t));
                    output.apply(Gate::cnot { control: c0, target: c1 });
                    output.apply(Gate::t(c0));
                    output.apply(Gate::tdg(c1));
                    output.apply(Gate::cnot { control: c0, target: c1 });
                }
                other => output.apply(other.clone()),
            }
        }
        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unitary::circuits_equiv;

    #[test]
    fn single() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 15);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::ccx { .. })));
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn preserves_non_ccx() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(2));
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 3);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn multiple() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 30);
        assert!(circuits_equiv(&c, &dec, 1e-10));
        let identity = Circuit::new(3);
        assert!(circuits_equiv(&c, &identity, 1e-10));
    }

    #[test]
    fn different_qubits() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx { control1: 2, control2: 0, target: 1 });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 15);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn mixed_circuit() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(2));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::ccx { control1: 1, control2: 2, target: 0 });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 32);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn four_qubit() {
        let mut c = Circuit::new(4);
        c.apply(Gate::h(3));
        c.apply(Gate::ccx { control1: 0, control2: 1, target: 2 });
        c.apply(Gate::cnot { control: 2, target: 3 });
        let dec = DecomposeToffoli.run(&c);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn empty() {
        let c = Circuit::new(3);
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 0);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn preserves_z_sdg() {
        let mut c = Circuit::new(3);
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(1));
        c.apply(Gate::s(2));
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 3);
        assert!(matches!(&dec.gates[0], Gate::z(0)));
        assert!(matches!(&dec.gates[1], Gate::sdg(1)));
        assert!(matches!(&dec.gates[2], Gate::s(2)));
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }
}
