use indicatif::ProgressBar;
use rayon::prelude::*;

use crate::circuit::{Circuit, Gate};
use crate::pass::Pass;
use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gridsynth::gridsynth_gates;

pub struct DecomposeRz {
    pub epsilon: f64,
}

impl Default for DecomposeRz {
    fn default() -> Self {
        Self { epsilon: 1e-10 }
    }
}

impl Pass for DecomposeRz {
    fn name(&self) -> &str {
        "Rz → Clifford+T decomposition"
    }

    fn run_with_progress(&self, circuit: &Circuit, pb: &ProgressBar) -> Circuit {
        let epsilon = self.epsilon;

        // Synthesize all Rz gates in parallel.
        let expanded: Vec<Vec<Gate>> = circuit.gates.par_iter().map(|gate| {
            let result = match gate {
                Gate::rz(theta, q) => {
                    let q = *q;
                    let chars = synthesize_rz(*theta, epsilon);
                    let mut gates = Vec::with_capacity(chars.len());
                    for g in chars {
                        match g {
                            'H' => gates.push(Gate::h(q)),
                            'T' => gates.push(Gate::t(q)),
                            'S' => gates.push(Gate::s(q)),
                            'X' => gates.push(Gate::x(q)),
                            'I' | 'W' => {} // identity / global phase, skip
                            c => eprintln!("warning: unknown gridsynth gate '{c}'"),
                        }
                    }
                    gates
                }
                other => vec![other.clone()],
            };
            pb.inc(1);
            result
        }).collect();

        let mut output = Circuit::new(circuit.num_qubits);
        for gates in expanded {
            for g in gates {
                output.apply(g);
            }
        }
        output
    }
}

fn synthesize_rz(theta: f64, epsilon: f64) -> Vec<char> {
    let mut config = config_from_theta_epsilon(theta, epsilon, 0, false, true);
    let result = gridsynth_gates(&mut config);
    result.gates.chars().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn decomposes_rz_into_clifford_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI / 5.0, 0));
        let dec = DecomposeRz { epsilon: 1e-3 }.run(&c);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::rz(..))));
        assert!(!dec.gates.is_empty());
        for g in &dec.gates {
            assert!(matches!(g, Gate::h(_) | Gate::t(_) | Gate::s(_) | Gate::x(_)));
        }
    }

    #[test]
    fn preserves_non_rz() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::t(0));
        let dec = DecomposeRz::default().run(&c);
        assert_eq!(dec.gates.len(), 3);
    }

    #[test]
    fn mixed_circuit() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::rz(PI / 3.0, 0));
        c.apply(Gate::cnot { control: 0, target: 1 });
        c.apply(Gate::rz(PI / 7.0, 1));
        let dec = DecomposeRz { epsilon: 1e-3 }.run(&c);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::rz(..))));
    }

    #[test]
    fn empty_circuit() {
        let c = Circuit::new(1);
        let dec = DecomposeRz::default().run(&c);
        assert_eq!(dec.gates.len(), 0);
    }

    #[test]
    fn default_epsilon_is_1e_10() {
        assert_eq!(DecomposeRz::default().epsilon, 1e-10);
    }

    #[test]
    fn coarser_epsilon_produces_fewer_or_equal_t_gates() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI / 5.0, 0));
        let fine = DecomposeRz { epsilon: 1e-4 }.run(&c);
        let coarse = DecomposeRz { epsilon: 1e-2 }.run(&c);
        let t_fine = fine.gates.iter().filter(|g| matches!(g, Gate::t(_))).count();
        let t_coarse = coarse.gates.iter().filter(|g| matches!(g, Gate::t(_))).count();
        assert!(t_coarse <= t_fine,
            "coarser epsilon should not require more T gates ({t_coarse} > {t_fine})");
    }
}
