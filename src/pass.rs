use crate::circuit::{Circuit, Gate, qubit_operands};

/// A circuit pass. The default output keeps optimization pipelines gate-based;
/// terminal conversions can select a different output type.
///
/// Implementors are not required to be `Sync`/`Send`: a pass may cache state
/// behind interior mutability that isn't safe to share across threads (see
/// `SuperOpt`'s `MatrixStore`, which uses `Rc<RefCell<_>>` since it's never
/// actually accessed from more than one thread — each parallel worker
/// constructs and owns its own pass instance).
pub trait Pass<Output = Circuit> {
    fn name(&self) -> &str;
    fn run(&self, circuit: &Circuit) -> Output;
}

/// The outcome of [`run_passes`].
pub struct PassResult {
    pub circuit: Circuit,
    /// T-count after only the first pass, for attributing reductions to it.
    pub t_after_first: usize,
    /// Gate count after only the first pass, for attributing reductions to it.
    pub gates_after_first: usize,
}

/// Run `passes` in order, feeding each pass's output to the next.
pub fn run_passes(circuit: &Circuit, passes: &[&dyn Pass]) -> PassResult {
    let mut c = circuit.clone();
    let mut t_after_first = 0;
    let mut gates_after_first = 0;
    for (i, p) in passes.iter().enumerate() {
        c = p.run(&c);
        if i == 0 {
            t_after_first = count_t(&c);
            gates_after_first = c.gates.len();
        }
    }
    PassResult {
        circuit: c,
        t_after_first,
        gates_after_first,
    }
}

/// Number of `t`/`tdg` gates in the circuit.
pub fn count_t(c: &Circuit) -> usize {
    c.gates
        .iter()
        .filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_)))
        .count()
}

/// Number of `rz` gates in the circuit.
pub fn count_rz(c: &Circuit) -> usize {
    c.gates.iter().filter(|g| matches!(g, Gate::rz(..))).count()
}

/// Number of two-qubit `cnot`/`cz` gates in the circuit.
pub fn count_2q(c: &Circuit) -> usize {
    c.gates
        .iter()
        .filter(|g| matches!(g, Gate::cnot { .. } | Gate::cz { .. }))
        .count()
}

/// Circuit depth, where gates on disjoint qubits may occupy the same layer.
pub fn depth(c: &Circuit) -> usize {
    let mut next_layer = vec![0; c.num_qubits];
    for gate in &c.gates {
        let (n, qs) = qubit_operands(gate);
        let layer = qs[..n]
            .iter()
            .map(|&q| next_layer[q as usize])
            .max()
            .unwrap_or(0)
            + 1;
        for &q in &qs[..n] {
            next_layer[q as usize] = layer;
        }
    }
    next_layer.into_iter().max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::{count_2q, count_rz, depth};
    use crate::circuit::{Circuit, Gate};

    #[test]
    fn count_2q_counts_cnot_and_cz_only() {
        let mut circuit = Circuit::new(3);
        circuit.apply(Gate::h(0));
        circuit.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        circuit.apply(Gate::cz {
            control: 1,
            target: 2,
        });
        circuit.apply(Gate::t(0));

        assert_eq!(count_2q(&circuit), 2);
    }

    #[test]
    fn count_rz_counts_only_rz_gates() {
        let mut circuit = Circuit::new(1);
        circuit.apply(Gate::rz(0.25, 0));
        circuit.apply(Gate::t(0));
        circuit.apply(Gate::rz(-0.5, 0));

        assert_eq!(count_rz(&circuit), 2);
    }

    #[test]
    fn depth_layers_disjoint_gates_and_serializes_shared_qubits() {
        let mut circuit = Circuit::new(3);
        circuit.apply(Gate::h(0));
        circuit.apply(Gate::h(1));
        circuit.apply(Gate::cnot {
            control: 0,
            target: 2,
        });
        circuit.apply(Gate::t(1));

        assert_eq!(depth(&circuit), 2);
    }

    #[test]
    fn depth_of_empty_circuit_is_zero() {
        assert_eq!(depth(&Circuit::new(4)), 0);
    }

    #[test]
    fn depth_serializes_gates_on_one_qubit() {
        let mut circuit = Circuit::new(1);
        circuit.apply(Gate::h(0));
        circuit.apply(Gate::t(0));
        circuit.apply(Gate::z(0));

        assert_eq!(depth(&circuit), 3);
    }

    #[test]
    fn depth_propagates_through_multi_qubit_dependencies() {
        let mut circuit = Circuit::new(3);
        circuit.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        circuit.apply(Gate::x(1));
        circuit.apply(Gate::cz {
            control: 1,
            target: 2,
        });

        assert_eq!(depth(&circuit), 3);
    }
}
