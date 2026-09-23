//! Gate-decomposition passes:
//!   * `DecomposeToffoli` — rewrites every `ccx` and `ccz` into Clifford+T.
//!   * `DecomposeCz` — rewrites every `cz` into `H · CNOT · H`.
//!   * `DecomposeRz` — synthesizes each `rz(θ)` into Clifford+T via
//!     gridsynth (`rsgridsynth`).

use std::sync::Mutex;

use crate::circuit::{Circuit, Gate, Qubit};
use crate::pass::Pass;
use rsgridsynth::config::config_from_theta_epsilon;
use rsgridsynth::gridsynth::gridsynth_gates;

// --- Toffoli decomposition --------------------------------------------------

/// Rewrites every `ccx` and `ccz` gate into Clifford+T.
pub struct DecomposeToffoli;

fn emit_ccx_decomposition(output: &mut Circuit, c0: Qubit, c1: Qubit, t: Qubit) {
    output.apply(Gate::h(t));
    output.apply(Gate::cnot {
        control: c1,
        target: t,
    });
    output.apply(Gate::tdg(t));
    output.apply(Gate::cnot {
        control: c0,
        target: t,
    });
    output.apply(Gate::t(t));
    output.apply(Gate::cnot {
        control: c1,
        target: t,
    });
    output.apply(Gate::tdg(t));
    output.apply(Gate::cnot {
        control: c0,
        target: t,
    });
    output.apply(Gate::t(c1));
    output.apply(Gate::t(t));
    output.apply(Gate::h(t));
    output.apply(Gate::cnot {
        control: c0,
        target: c1,
    });
    output.apply(Gate::t(c0));
    output.apply(Gate::tdg(c1));
    output.apply(Gate::cnot {
        control: c0,
        target: c1,
    });
}

impl Pass for DecomposeToffoli {
    fn name(&self) -> &str {
        "Toffoli decomposition"
    }
    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
        for gate in &circuit.gates {
            match gate {
                Gate::ccx {
                    control1,
                    control2,
                    target,
                } => {
                    emit_ccx_decomposition(&mut output, *control1, *control2, *target);
                }
                Gate::ccz {
                    control1,
                    control2,
                    target,
                } => {
                    // CCZ = H(target) · CCX · H(target). Keep this explicit;
                    // the cancellation pass removes the adjacent Hadamards when
                    // it follows decomposition in an optimization pipeline.
                    output.apply(Gate::h(*target));
                    emit_ccx_decomposition(&mut output, *control1, *control2, *target);
                    output.apply(Gate::h(*target));
                }
                other => output.apply(other.clone()),
            }
        }
        output
    }
}

// --- CZ decomposition --------------------------------------------------------

/// Explicitly lower CZ gates for backends that only accept H and CNOT.
/// CZ remains native unless this pass is requested.
pub struct DecomposeCz;

impl Pass for DecomposeCz {
    fn name(&self) -> &str {
        "CZ decomposition"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
        for gate in &circuit.gates {
            match gate {
                Gate::cz { control, target } => {
                    output.apply(Gate::h(*target));
                    output.apply(Gate::cnot {
                        control: *control,
                        target: *target,
                    });
                    output.apply(Gate::h(*target));
                }
                other => output.apply(other.clone()),
            }
        }
        output
    }
}

// --- Rz decomposition via gridsynth -----------------------------------------

/// Synthesizes each `rz(θ)` into Clifford+T via gridsynth (`rsgridsynth`).
pub struct DecomposeRz {
    /// Approximation precision; smaller is more accurate but produces a
    /// larger decomposition. Defaults to `1e-10`.
    pub epsilon: f64,
}

// rsgridsynth stores working precision in process-global state. Serialize
// configuration and synthesis so one call cannot change another's precision.
static GRIDSYNTH_BATCH_LOCK: Mutex<()> = Mutex::new(());

impl Default for DecomposeRz {
    fn default() -> Self {
        Self { epsilon: 1e-10 }
    }
}

impl Pass for DecomposeRz {
    fn name(&self) -> &str {
        "Rz → Clifford+T decomposition"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let epsilon = self.epsilon;
        let _batch_guard = GRIDSYNTH_BATCH_LOCK
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());

        // Keep gridsynth calls serial: the dependency has process-global state.
        let expanded: Vec<Vec<Gate>> = circuit
            .gates
            .iter()
            .map(|gate| {
                match gate {
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
                }
            })
            .collect();

        let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
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

// --- Tests ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unitary::circuits_equiv;
    use std::f64::consts::PI;

    // --- Toffoli decomposition tests ---

    #[test]
    fn toffoli_single() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 15);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::ccx { .. })));
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn ccz_single() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccz {
            control1: 0,
            control2: 1,
            target: 2,
        });

        let dec = DecomposeToffoli.run(&c);

        assert_eq!(dec.gates.len(), 17);
        assert!(
            !dec.gates
                .iter()
                .any(|g| matches!(g, Gate::ccx { .. } | Gate::ccz { .. }))
        );
        assert!(!dec.has_toffoli());
        assert!(!dec.has_ccz());
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn decomposes_mixed_ccx_and_ccz() {
        let mut c = Circuit::new(4);
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::ccz {
            control1: 3,
            control2: 1,
            target: 0,
        });

        let dec = DecomposeToffoli.run(&c);

        assert_eq!(dec.gates.len(), 32);
        assert!(
            !dec.gates
                .iter()
                .any(|g| matches!(g, Gate::ccx { .. } | Gate::ccz { .. }))
        );
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn ccz_decomposition_handles_every_operand_order() {
        for [control1, control2, target] in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let mut c = Circuit::new(3);
            c.apply(Gate::ccz {
                control1,
                control2,
                target,
            });

            let dec = DecomposeToffoli.run(&c);

            assert!(circuits_equiv(&c, &dec, 1e-10));
        }
    }

    #[test]
    fn ccz_decomposition_is_idempotent() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccz {
            control1: 0,
            control2: 1,
            target: 2,
        });

        let once = DecomposeToffoli.run(&c);
        let twice = DecomposeToffoli.run(&once);

        assert_eq!(once.to_qasm(), twice.to_qasm());
    }

    #[test]
    fn ccz_decomposition_preserves_measurement_metadata() {
        let mut c = Circuit::with_cbits(3, 1);
        c.apply(Gate::reset(0));
        c.apply(Gate::ccz {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::measure { qubit: 2, cbit: 0 });

        let dec = DecomposeToffoli.run(&c);

        assert_eq!(dec.num_cbits, 1);
        assert!(dec.has_measurement());
        assert!(matches!(dec.gates.first(), Some(Gate::reset(0))));
        assert!(matches!(
            dec.gates.last(),
            Some(Gate::measure { qubit: 2, cbit: 0 })
        ));
    }

    #[test]
    fn toffoli_preserves_non_ccx() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(2));
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 3);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn toffoli_multiple() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 30);
        assert!(circuits_equiv(&c, &dec, 1e-10));
        let identity = Circuit::new(3);
        assert!(circuits_equiv(&c, &identity, 1e-10));
    }

    #[test]
    fn toffoli_different_qubits() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccx {
            control1: 2,
            control2: 0,
            target: 1,
        });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 15);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn toffoli_mixed_circuit() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(2));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::ccx {
            control1: 1,
            control2: 2,
            target: 0,
        });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 32);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn toffoli_four_qubit() {
        let mut c = Circuit::new(4);
        c.apply(Gate::h(3));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::cnot {
            control: 2,
            target: 3,
        });
        let dec = DecomposeToffoli.run(&c);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn toffoli_empty() {
        let c = Circuit::new(3);
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 0);
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn toffoli_preserves_z_sdg() {
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

    #[test]
    fn toffoli_preserves_measure() {
        let mut c = Circuit::with_cbits(2, 2);
        c.apply(Gate::h(0));
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::measure { qubit: 1, cbit: 1 });
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 3);
        assert!(matches!(&dec.gates[1], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(matches!(&dec.gates[2], Gate::measure { qubit: 1, cbit: 1 }));
    }

    #[test]
    fn toffoli_preserves_reset() {
        let mut c = Circuit::new(2);
        c.apply(Gate::reset(0));
        c.apply(Gate::h(1));
        c.apply(Gate::reset(1));
        let dec = DecomposeToffoli.run(&c);
        assert_eq!(dec.gates.len(), 3);
        assert!(matches!(&dec.gates[0], Gate::reset(0)));
        assert!(matches!(&dec.gates[2], Gate::reset(1)));
    }

    #[test]
    fn ccx_decomposes_with_surrounding_measure() {
        // Toffoli decomposes; surrounding measurements / resets pass through unchanged.
        let mut c = Circuit::with_cbits(3, 1);
        c.apply(Gate::reset(0));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::measure { qubit: 2, cbit: 0 });
        let dec = DecomposeToffoli.run(&c);
        // 1 (reset) + 15 (CCX decomposition) + 1 (measure) = 17
        assert_eq!(dec.gates.len(), 17);
        assert!(matches!(&dec.gates[0], Gate::reset(0)));
        assert!(matches!(
            dec.gates.last().unwrap(),
            Gate::measure { qubit: 2, cbit: 0 }
        ));
        assert_eq!(dec.num_cbits, 1);
    }

    // --- CZ decomposition tests ---

    #[test]
    fn cz_decomposes_to_h_cnot_h() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        let dec = DecomposeCz.run(&c);
        assert_eq!(dec.gates.len(), 3);
        assert!(matches!(dec.gates[0], Gate::h(1)));
        assert!(matches!(
            dec.gates[1],
            Gate::cnot {
                control: 0,
                target: 1
            }
        ));
        assert!(matches!(dec.gates[2], Gate::h(1)));
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn cz_decomposition_respects_operand_order() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cz {
            control: 2,
            target: 0,
        });
        let dec = DecomposeCz.run(&c);
        assert!(matches!(dec.gates[0], Gate::h(0)));
        assert!(matches!(
            dec.gates[1],
            Gate::cnot {
                control: 2,
                target: 0
            }
        ));
        assert!(matches!(dec.gates[2], Gate::h(0)));
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn cz_decomposes_multiple_and_preserves_other_gates() {
        let mut c = Circuit::new(3);
        c.apply(Gate::t(0));
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::x(2));
        c.apply(Gate::cz {
            control: 2,
            target: 1,
        });
        let dec = DecomposeCz.run(&c);
        assert_eq!(dec.gates.len(), 8);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::cz { .. })));
        assert!(matches!(dec.gates[0], Gate::t(0)));
        assert!(matches!(dec.gates[4], Gate::x(2)));
        assert!(circuits_equiv(&c, &dec, 1e-10));
    }

    #[test]
    fn cz_decomposition_preserves_measurement_metadata() {
        let mut c = Circuit::with_cbits(2, 1);
        c.apply(Gate::reset(0));
        c.apply(Gate::cz {
            control: 0,
            target: 1,
        });
        c.apply(Gate::measure { qubit: 1, cbit: 0 });
        let dec = DecomposeCz.run(&c);
        assert_eq!(dec.num_cbits, 1);
        assert!(dec.has_measurement());
        assert!(matches!(dec.gates[0], Gate::reset(0)));
        assert!(matches!(
            dec.gates.last(),
            Some(Gate::measure { qubit: 1, cbit: 0 })
        ));
    }

    #[test]
    fn other_decomposers_preserve_native_cz() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cz {
            control: 0,
            target: 2,
        });
        let toffoli = DecomposeToffoli.run(&c);
        let rz = DecomposeRz::default().run(&c);
        assert!(matches!(
            toffoli.gates.as_slice(),
            [Gate::cz {
                control: 0,
                target: 2
            }]
        ));
        assert!(matches!(
            rz.gates.as_slice(),
            [Gate::cz {
                control: 0,
                target: 2
            }]
        ));
    }

    #[test]
    fn cz_decomposition_empty_circuit() {
        let c = Circuit::new(2);
        let dec = DecomposeCz.run(&c);
        assert!(dec.gates.is_empty());
        assert_eq!(dec.num_qubits, 2);
    }

    #[test]
    fn cz_decomposition_is_idempotent() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cz {
            control: 0,
            target: 2,
        });
        c.apply(Gate::t(1));
        c.apply(Gate::cz {
            control: 2,
            target: 1,
        });
        let once = DecomposeCz.run(&c);
        let twice = DecomposeCz.run(&once);
        assert_eq!(once.to_qasm(), twice.to_qasm());
        assert!(circuits_equiv(&c, &twice, 1e-10));
    }

    // --- Rz decomposition tests ---

    #[test]
    fn rz_decomposes_into_clifford_t() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI / 5.0, 0));
        let dec = DecomposeRz { epsilon: 1e-3 }.run(&c);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::rz(..))));
        assert!(!dec.gates.is_empty());
        for g in &dec.gates {
            assert!(matches!(
                g,
                Gate::h(_) | Gate::t(_) | Gate::s(_) | Gate::x(_)
            ));
        }
    }

    #[test]
    fn rz_multiple_angles_do_not_reuse_an_invalid_gridsynth_solution() {
        // rsgridsynth 0.2.0 cached these by coefficients alone and could reuse
        // a result with the wrong denominator exponent for the second angle.
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(0.02, 0));
        c.apply(Gate::rz(0.03, 0));

        let dec = DecomposeRz { epsilon: 1e-3 }.run(&c);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::rz(..))));
        assert!(circuits_equiv(&c, &dec, 2e-3));
    }

    #[test]
    fn rz_preserves_non_rz() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::t(0));
        let dec = DecomposeRz::default().run(&c);
        assert_eq!(dec.gates.len(), 3);
    }

    #[test]
    fn rz_mixed_circuit() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::rz(PI / 3.0, 0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::rz(PI / 7.0, 1));
        let dec = DecomposeRz { epsilon: 1e-3 }.run(&c);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::rz(..))));
    }

    #[test]
    fn rz_empty_circuit() {
        let c = Circuit::new(1);
        let dec = DecomposeRz::default().run(&c);
        assert_eq!(dec.gates.len(), 0);
    }

    #[test]
    fn rz_default_epsilon_is_1e_10() {
        assert_eq!(DecomposeRz::default().epsilon, 1e-10);
    }

    #[test]
    fn rz_preserves_measure_and_reset() {
        let mut c = Circuit::with_cbits(2, 1);
        c.apply(Gate::reset(0));
        c.apply(Gate::rz(PI / 5.0, 0));
        c.apply(Gate::measure { qubit: 1, cbit: 0 });
        let dec = DecomposeRz { epsilon: 1e-3 }.run(&c);
        // No rz survives; reset and measure both still there.
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::rz(..))));
        assert!(dec.gates.iter().any(|g| matches!(g, Gate::reset(0))));
        assert!(
            dec.gates
                .iter()
                .any(|g| matches!(g, Gate::measure { qubit: 1, cbit: 0 }))
        );
        assert_eq!(dec.num_cbits, 1);
        assert!(dec.has_measurement());
    }

    #[test]
    fn rz_coarser_epsilon_produces_fewer_or_equal_t_gates() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(PI / 5.0, 0));
        let fine = DecomposeRz { epsilon: 1e-4 }.run(&c);
        let coarse = DecomposeRz { epsilon: 1e-2 }.run(&c);
        let t_fine = fine
            .gates
            .iter()
            .filter(|g| matches!(g, Gate::t(_)))
            .count();
        let t_coarse = coarse
            .gates
            .iter()
            .filter(|g| matches!(g, Gate::t(_)))
            .count();
        assert!(
            t_coarse <= t_fine,
            "coarser epsilon should not require more T gates ({t_coarse} > {t_fine})"
        );
    }
}
