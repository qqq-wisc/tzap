use crate::circuit::{Circuit, Gate, qubit_operands};

/// Per-pipeline exact arithmetic state. Workers own independent contexts.
/// Stored circuit angles own their nodes and outlive this state.
pub struct PassState {
    pub angles: crate::angle::AngleContext,
}
impl Default for PassState {
    fn default() -> Self {
        Self {
            angles: crate::angle::AngleContext::default().with_run_budget(),
        }
    }
}
impl PassState {
    pub fn with_limits(limits: crate::angle::Limits) -> Self {
        Self {
            angles: crate::angle::AngleContext::with_limits(limits).with_run_budget(),
        }
    }
}

/// A circuit pass. The default output keeps optimization pipelines gate-based;
/// terminal conversions can select a different output type. Whole-circuit gate
/// passes retain the input qubit and classical-bit counts.
///
/// Implementors are not required to be `Sync`/`Send`: a pass may cache state
/// behind interior mutability that isn't safe to share across threads (see
/// `SuperOpt`'s `MatrixStore`, which uses `Rc<RefCell<_>>` since it's never
/// actually accessed from more than one thread — each parallel worker
/// constructs and owns its own pass instance).
pub trait Pass<Output = Circuit> {
    fn name(&self) -> &str;
    fn run(&self, circuit: &Circuit) -> Output;
    fn run_with_state(&self, circuit: &Circuit, _state: &mut PassState) -> Output {
        self.run(circuit)
    }
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
    let mut state = PassState::default();
    for (i, p) in passes.iter().enumerate() {
        c = p.run_with_state(&c, &mut state);
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

/// T/Tdg cost of named gates and exact phase/axis rotation terms.
/// Controlled rotations contribute two half-angle terms (CP contributes three).
/// Angles requiring synthesis contribute to `count_rz`, not an estimated T cost.
pub fn count_t(c: &Circuit) -> usize {
    c.gates.iter().map(|g| rotation_counts(g).0).sum()
}

/// Number of rotation terms that still require arbitrary-angle synthesis.
/// P/Rz/Rx/Ry contribute one term; controlled rotations contribute two,
/// or three for CP. Exact Clifford+T terms are counted separately by `count_t`.
pub fn count_rz(c: &Circuit) -> usize {
    c.gates.iter().map(|g| rotation_counts(g).1).sum()
}

/// Shared by the public counters and the single-traversal progress metrics.
pub(crate) fn rotation_counts(gate: &Gate) -> (usize, usize) {
    fn angle_counts(angle: &crate::angle::Angle) -> (usize, usize) {
        match angle.projective_quarters() {
            Some(k) => (usize::from(k % 2 == 1), 0),
            None => match angle.affine_parts() {
                Some((phase, _)) => {
                    let (t, rz) = angle_counts(&phase);
                    (t, rz + 1)
                }
                None => (0, 1),
            },
        }
    }
    match gate {
        Gate::t(_) | Gate::tdg(_) => (1, 0),
        Gate::rz(a, _) | Gate::p(a, _) | Gate::rx(a, _) | Gate::ry(a, _) => angle_counts(a),
        Gate::cp { lambda: a, .. }
        | Gate::crx { theta: a, .. }
        | Gate::cry { theta: a, .. }
        | Gate::crz { theta: a, .. } => {
            let terms = if matches!(gate, Gate::cp { .. }) {
                3
            } else {
                2
            };
            let (t, rz) = a.checked_half().map_or((0, 1), |half| angle_counts(&half));
            (terms * t, terms * rz)
        }
        Gate::ch { .. } => (2, 0),
        _ => (0, 0),
    }
}

/// Number of two-qubit gates in the circuit.
pub fn count_2q(c: &Circuit) -> usize {
    c.gates
        .iter()
        .filter(|g| {
            matches!(
                g,
                Gate::cnot { .. }
                    | Gate::cz { .. }
                    | Gate::swap(..)
                    | Gate::cy { .. }
                    | Gate::ch { .. }
                    | Gate::cp { .. }
                    | Gate::crx { .. }
                    | Gate::cry { .. }
                    | Gate::crz { .. }
            )
        })
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
    use super::{count_2q, count_rz, count_t, depth};
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
    fn count_rz_counts_numeric_rotation_terms() {
        let mut circuit = Circuit::new(1);
        circuit.apply(Gate::rz_f64(0.25, 0).unwrap());
        circuit.apply(Gate::t(0));
        circuit.apply(Gate::rz_f64(-0.5, 0).unwrap());

        assert_eq!(count_rz(&circuit), 2);
    }

    #[test]
    fn rotation_resource_counts_include_native_gates_and_exact_angles() {
        for (gate, t, rz) in [
            ("p(pi/4) q[0]", 1, 0),
            ("rz(pi/4) q[0]", 1, 0),
            ("rx(pi/4) q[0]", 1, 0),
            ("ry(-3*pi/4) q[0]", 1, 0),
            ("p(pi/2) q[0]", 0, 0),
            ("rz(0) q[0]", 0, 0),
            ("cp(pi/4) q[0],q[1]", 0, 3),
            ("cp(pi/2) q[0],q[1]", 3, 0),
            ("crz(pi/2) q[0],q[1]", 2, 0),
            ("crx(-pi/2) q[0],q[1]", 2, 0),
            ("cry(3*pi/2) q[0],q[1]", 2, 0),
            ("crz(2*pi) q[0],q[1]", 0, 0),
            ("crz(4*pi) q[0],q[1]", 0, 0),
            ("ch q[0],q[1]", 2, 0),
            ("rx(0.7853981633974483) q[0]", 0, 1),
            ("ry(0.3) q[0]", 0, 1),
            ("cp(0.3) q[0],q[1]", 0, 3),
            ("crz(0.3) q[0],q[1]", 0, 2),
            ("rz(pi/4 + 0.3) q[0]", 1, 1),
        ] {
            let circuit = Circuit::from_qasm(&format!("qreg q[2]; {gate};")).unwrap();
            assert_eq!((count_t(&circuit), count_rz(&circuit)), (t, rz), "{gate}");
            let metrics = crate::optimize::Metrics::of(&circuit);
            assert_eq!((metrics.t, metrics.rz), (t, rz), "{gate}");
        }
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
