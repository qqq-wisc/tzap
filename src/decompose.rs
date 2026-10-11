//! Gate-decomposition passes:
//!   * `DecomposeToffoli` — rewrites every `ccx` and `ccz` into Clifford+T.
//!   * `DecomposeCz` — rewrites every `cz` into `H · CNOT · H`.
//!   * `DecomposeRotations` — lowers phase and ordinary/controlled rotations via
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

/// Synthesizes phase and ordinary/controlled axis rotations into Clifford+T.
pub struct DecomposeRotations {
    /// Approximation precision; smaller is more accurate but produces a
    /// larger decomposition. Defaults to `1e-10`.
    pub epsilon: f64,
}

pub type DecomposeRz = DecomposeRotations;

// rsgridsynth stores working precision in process-global state. Serialize
// configuration and synthesis so one call cannot change another's precision.
static GRIDSYNTH_BATCH_LOCK: Mutex<()> = Mutex::new(());

impl Default for DecomposeRotations {
    fn default() -> Self {
        Self { epsilon: 1e-10 }
    }
}

impl Pass for DecomposeRotations {
    fn name(&self) -> &str {
        "Rotations → Clifford+T decomposition"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        self.try_run(circuit).unwrap_or_else(|_| circuit.clone())
    }
}

/// A synthesis failure; preserved expressions have no permitted numeric value.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SynthesisError(pub String);
impl std::fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for SynthesisError {}
impl DecomposeRotations {
    /// Explicit approximate synthesis. The epsilon applies to gridsynth's
    /// numeric target, not to conversion error or the whole circuit.
    pub fn try_run(&self, circuit: &Circuit) -> Result<Circuit, SynthesisError> {
        let epsilon = self.epsilon;
        if !epsilon.is_finite() || epsilon <= 0.0 {
            return Err(SynthesisError(
                "Rz synthesis epsilon must be positive and finite".into(),
            ));
        }
        let _batch_guard = GRIDSYNTH_BATCH_LOCK
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());

        // Keep gridsynth calls serial: the dependency has process-global state.
        let expanded: Result<Vec<Vec<Gate>>, SynthesisError> = circuit
            .gates
            .iter()
            .map(|gate| {
                if let Some(terms) = expand_controlled_rotation(gate)? {
                    let count = terms
                        .iter()
                        .filter(|g| matches!(g, Gate::rz(..) | Gate::p(..)))
                        .count();
                    let term_epsilon = epsilon / count as f64;
                    if term_epsilon == 0.0 {
                        return Err(SynthesisError(
                            "rotation synthesis epsilon underflows when split across terms".into(),
                        ));
                    }
                    terms
                        .iter()
                        .map(|g| synthesize_single_rotation(g, term_epsilon))
                        .collect::<Result<Vec<_>, _>>()
                        .map(|groups| groups.into_iter().flatten().collect())
                } else {
                    synthesize_single_rotation(gate, epsilon)
                }
            })
            .collect();

        let mut output = Circuit::with_cbits(circuit.num_qubits, circuit.num_cbits);
        for gates in expanded? {
            for g in gates {
                output.apply(g);
            }
        }
        Ok(output)
    }
}

/// Expand controlled rotations exactly before projective single-qubit synthesis.
/// Half-angles retain the literal 4π period of controlled RX/RY/RZ.
pub(crate) fn expand_controlled_rotation(gate: &Gate) -> Result<Option<Vec<Gate>>, SynthesisError> {
    let (theta, control, target) = match gate {
        Gate::cp {
            lambda,
            control,
            target,
        } => (lambda, *control, *target),
        Gate::crx {
            theta,
            control,
            target,
        }
        | Gate::cry {
            theta,
            control,
            target,
        }
        | Gate::crz {
            theta,
            control,
            target,
        } => (theta, *control, *target),
        _ => return Ok(None),
    };
    let half = theta
        .checked_half()
        .map_err(|e| SynthesisError(e.to_string()))?;
    let negative = half
        .literal_neg()
        .map_err(|e| SynthesisError(e.to_string()))?;
    let mut terms = Vec::with_capacity(9);
    if matches!(gate, Gate::cp { .. }) {
        terms.push(Gate::p(half.clone(), control));
    }
    if matches!(gate, Gate::cry { .. }) {
        terms.push(Gate::sdg(target));
    }
    if matches!(gate, Gate::crx { .. } | Gate::cry { .. }) {
        terms.push(Gate::h(target));
    }
    let cx = Gate::cnot { control, target };
    terms.extend([
        Gate::rz(half, target),
        cx.clone(),
        Gate::rz(negative, target),
        cx,
    ]);
    if matches!(gate, Gate::crx { .. } | Gate::cry { .. }) {
        terms.push(Gate::h(target));
    }
    if matches!(gate, Gate::cry { .. }) {
        terms.push(Gate::s(target));
    }
    Ok(Some(terms))
}

fn synthesize_single_rotation(gate: &Gate, epsilon: f64) -> Result<Vec<Gate>, SynthesisError> {
    let (theta, q) = match gate {
        Gate::rz(a, q) | Gate::p(a, q) | Gate::rx(a, q) | Gate::ry(a, q) => (a, *q),
        _ => return Ok(vec![gate.clone()]),
    };
    let mut gates = Vec::new();
    if matches!(gate, Gate::ry(..)) {
        gates.push(Gate::sdg(q));
    }
    if matches!(gate, Gate::rx(..) | Gate::ry(..)) {
        gates.push(Gate::h(q));
    }
    if theta.quarter_turns().is_some() {
        gates.extend(theta.rotation_gates(q));
    } else {
        let numeric = theta
            .to_f64_lossy()
            .map_err(|e| SynthesisError(e.to_string()))?;
        crate::angle_stats::synthesis();
        for c in synthesize_rz(numeric, epsilon) {
            match c {
                'H' => gates.push(Gate::h(q)),
                'T' => gates.push(Gate::t(q)),
                'S' => gates.push(Gate::s(q)),
                'X' => gates.push(Gate::x(q)),
                'I' | 'W' => {}
                _ => return Err(SynthesisError(format!("unknown gridsynth gate {c:?}"))),
            }
        }
    }
    if matches!(gate, Gate::rx(..) | Gate::ry(..)) {
        gates.push(Gate::h(q));
    }
    if matches!(gate, Gate::ry(..)) {
        gates.push(Gate::s(q));
    }
    Ok(gates)
}

fn synthesize_rz(theta: f64, epsilon: f64) -> Vec<char> {
    // The dependency subtracts unsigned decimal exponents when configuring
    // precision. A looser request can safely use the stricter bound of one.
    let epsilon = epsilon.min(1.0);
    let mut config =
        config_from_theta_epsilon(bounded_synthesis_angle(theta), epsilon, 0, false, true);
    let result = gridsynth_gates(&mut config);
    result.gates.chars().collect()
}

fn bounded_synthesis_angle(theta: f64) -> f64 {
    // Preserve ordinary synthesis, including literal 2pi and 4pi rotations.
    // Beyond this modest interval the dependency's fixed working precision
    // cannot reliably reduce large arguments, and its decimal conversion can
    // change the phase of a binary64 angle by an arbitrarily large amount.
    if theta.abs() <= 2.0 * std::f64::consts::TAU {
        return theta;
    }
    // These are the same eigenphases used by the native Rz reference matrix.
    // Unlike remainder with a rounded 2pi, native sin_cos performs argument
    // reduction before atan2 returns a bounded phase in [-pi, pi]. This occurs
    // only after controlled gates have expanded into unconditional rotations;
    // any discarded synthesis phase therefore remains whole-circuit global.
    let (sin, cos) = (theta / 2.0).sin_cos();
    2.0 * sin.atan2(cos)
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
        let rz = DecomposeRotations::default().run(&c);
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
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 5.0", 1).unwrap(),
            0,
        ));
        let dec = DecomposeRotations { epsilon: 1e-3 }.run(&c);
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
        c.apply(Gate::rz_f64(0.02, 0).unwrap());
        c.apply(Gate::rz_f64(0.03, 0).unwrap());

        let dec = DecomposeRotations { epsilon: 1e-3 }.run(&c);
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
        let dec = DecomposeRotations::default().run(&c);
        assert_eq!(dec.gates.len(), 3);
    }

    #[test]
    fn rz_mixed_circuit() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 3.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 7.0", 1).unwrap(),
            1,
        ));
        let dec = DecomposeRotations { epsilon: 1e-3 }.run(&c);
        assert!(!dec.gates.iter().any(|g| matches!(g, Gate::rz(..))));
    }

    #[test]
    fn rz_empty_circuit() {
        let c = Circuit::new(1);
        let dec = DecomposeRotations::default().run(&c);
        assert_eq!(dec.gates.len(), 0);
    }

    #[test]
    fn rz_default_epsilon_is_1e_10() {
        assert_eq!(DecomposeRotations::default().epsilon, 1e-10);
    }

    #[test]
    fn rz_preserves_measure_and_reset() {
        let mut c = Circuit::with_cbits(2, 1);
        c.apply(Gate::reset(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 5.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::measure { qubit: 1, cbit: 0 });
        let dec = DecomposeRotations { epsilon: 1e-3 }.run(&c);
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
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 5.0", 1).unwrap(),
            0,
        ));
        let fine = DecomposeRotations { epsilon: 1e-4 }.run(&c);
        let coarse = DecomposeRotations { epsilon: 1e-2 }.run(&c);
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

    fn spectral_error(input: &Circuit, output: &Circuit) -> f64 {
        use crate::unitary::{C, circuit_unitary};
        let (a, b) = (circuit_unitary(input), circuit_unitary(output));
        let n = a.len();
        let overlap = a
            .iter()
            .flatten()
            .zip(b.iter().flatten())
            .fold(C::ZERO, |sum, (&a, &b)| sum + a.conj() * b);
        assert!(overlap.norm_sq().is_finite() && overlap.norm_sq() > 0.0);
        let phase = overlap * (1.0 / overlap.norm_sq().sqrt());
        let diff = (0..n)
            .map(|row| {
                (0..n)
                    .map(|col| b[row][col] - phase * a[row][col])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let gram = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        (0..n).fold(C::ZERO, |sum, row| sum + diff[row][i].conj() * diff[row][j])
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut maximum: f64 = 0.0;
        for seed in 0..n {
            let mut vector = vec![C::ZERO; n];
            vector[seed] = C::ONE;
            for _ in 0..64 {
                let next = gram
                    .iter()
                    .map(|row| {
                        row.iter()
                            .zip(&vector)
                            .fold(C::ZERO, |sum, (&entry, &v)| sum + entry * v)
                    })
                    .collect::<Vec<_>>();
                let norm = next.iter().map(|entry| entry.norm_sq()).sum::<f64>().sqrt();
                assert!(norm.is_finite());
                if norm == 0.0 {
                    break;
                }
                vector = next.into_iter().map(|entry| entry * (1.0 / norm)).collect();
            }
            let error = diff
                .iter()
                .map(|row| {
                    row.iter()
                        .zip(&vector)
                        .fold(C::ZERO, |sum, (&entry, &v)| sum + entry * v)
                        .norm_sq()
                })
                .sum::<f64>()
                .sqrt();
            maximum = maximum.max(error);
        }
        maximum
    }

    #[test]
    fn bounded_synthesis_preserves_ordinary_angles_and_large_native_eigenphases() {
        let limit = 2.0 * std::f64::consts::TAU;
        for theta in [-limit, -7.31, -PI, -0.0, 0.0, PI, 7.31, limit] {
            assert_eq!(bounded_synthesis_angle(theta).to_bits(), theta.to_bits());
        }
        for theta in [
            f64::from_bits(limit.to_bits() + 1),
            1e16,
            1e20,
            1e30,
            f64::MAX,
        ] {
            for theta in [theta, -theta] {
                let reduced = bounded_synthesis_angle(theta);
                assert!(reduced.is_finite() && reduced.abs() <= std::f64::consts::TAU);
                let (old_sin, old_cos) = (theta / 2.0).sin_cos();
                let (new_sin, new_cos) = (reduced / 2.0).sin_cos();
                assert!((new_sin - old_sin).abs() < 1e-15);
                assert!((new_cos - old_cos).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn large_and_extreme_native_rotations_preserve_per_gate_synthesis_precision() {
        let epsilon = 1e-3;
        let template = Circuit::new(2);
        let mut checked = 0;
        for magnitude in [
            1e16,
            1e20,
            1e30,
            f64::MAX,
            f64::MIN_POSITIVE,
            f64::from_bits(1),
            2.0 * PI,
            4.0 * PI,
        ] {
            for theta in [magnitude, -magnitude] {
                for (control, target) in [(0, 1), (1, 0)] {
                    for gate in [
                        Gate::p_f64(theta, target).unwrap(),
                        Gate::rz_f64(theta, target).unwrap(),
                        Gate::rx_f64(theta, target).unwrap(),
                        Gate::ry_f64(theta, target).unwrap(),
                        Gate::cp {
                            lambda: crate::angle::Angle::from_f64(theta).unwrap(),
                            control,
                            target,
                        },
                        Gate::crx {
                            theta: crate::angle::Angle::from_f64(theta).unwrap(),
                            control,
                            target,
                        },
                        Gate::cry {
                            theta: crate::angle::Angle::from_f64(theta).unwrap(),
                            control,
                            target,
                        },
                        Gate::crz {
                            theta: crate::angle::Angle::from_f64(theta).unwrap(),
                            control,
                            target,
                        },
                    ] {
                        let input = template.replacing_gates(vec![
                            Gate::h(control),
                            Gate::h(target),
                            Gate::cnot { control, target },
                            gate.clone(),
                            Gate::t(control),
                        ]);
                        if matches!(
                            gate,
                            Gate::cp { .. }
                                | Gate::crx { .. }
                                | Gate::cry { .. }
                                | Gate::crz { .. }
                        ) && theta.abs() == f64::from_bits(1)
                        {
                            assert!(DecomposeRotations { epsilon }.try_run(&input).is_err());
                            continue;
                        }
                        let output = DecomposeRotations { epsilon }.try_run(&input).unwrap();
                        assert!(
                            output
                                .gates
                                .iter()
                                .all(|gate| gate.kind().has_exact_matrix())
                        );
                        assert_eq!(
                            (input.num_qubits, input.num_cbits),
                            (output.num_qubits, output.num_cbits)
                        );
                        let error = spectral_error(&input, &output);
                        assert!(
                            error.is_finite() && error <= epsilon + 1e-12,
                            "gate={gate:?}, error={error}, epsilon={epsilon}"
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert_eq!(checked, 240);
    }

    #[test]
    fn controlled_rotations_and_phases_synthesize_with_per_gate_error_bound() {
        let epsilon = 1e-3;
        for theta in [-2.0 * PI, -0.73, 0.0, PI / 2.0, 2.0 * PI, 7.31] {
            for (control, target) in [(0, 1), (1, 0)] {
                for gate in [
                    Gate::p_f64(theta, target).unwrap(),
                    Gate::cp {
                        lambda: crate::angle::Angle::from_f64(theta).unwrap(),
                        control,
                        target,
                    },
                    Gate::crx {
                        theta: crate::angle::Angle::from_f64(theta).unwrap(),
                        control,
                        target,
                    },
                    Gate::cry {
                        theta: crate::angle::Angle::from_f64(theta).unwrap(),
                        control,
                        target,
                    },
                    Gate::crz {
                        theta: crate::angle::Angle::from_f64(theta).unwrap(),
                        control,
                        target,
                    },
                ] {
                    let input = {
                        let mut constructed = Circuit::with_cbits(2, 0);
                        constructed.gates =
                            vec![Gate::h(control), Gate::h(target), gate, Gate::t(control)];
                        constructed
                    };
                    let output = DecomposeRotations { epsilon }.run(&input);
                    assert!(output.gates.iter().all(|g| matches!(
                        g,
                        Gate::h(_)
                            | Gate::x(_)
                            | Gate::s(_)
                            | Gate::sdg(_)
                            | Gate::t(_)
                            | Gate::tdg(_)
                            | Gate::cnot { .. }
                    )));
                    assert!(
                        circuits_equiv(&input, &output, epsilon),
                        "theta={theta}, control={control}, input={input}"
                    );
                }
            }
        }
    }

    #[test]
    fn rotations_all_axes_angles_and_wires_preserve_semantics() {
        let epsilon = 1e-3;
        for theta in [-2.0 * PI, -PI / 2.0, -0.73, 0.0, PI / 5.0, 2.0 * PI, 7.31] {
            for q in 0..2 {
                for gate in [
                    Gate::rz_f64(theta, q).unwrap(),
                    Gate::rx_f64(theta, q).unwrap(),
                    Gate::ry_f64(theta, q).unwrap(),
                ] {
                    let mut input = Circuit::new(2);
                    input.gates = vec![
                        Gate::h(0),
                        Gate::cnot {
                            control: 0,
                            target: 1,
                        },
                        gate,
                        Gate::t(1),
                    ];
                    let output = DecomposeRotations { epsilon }.run(&input);
                    assert!(!output.gates.iter().any(|g| matches!(
                        g,
                        Gate::rz(..) | Gate::rx(..) | Gate::ry(..) | Gate::p(..)
                    )));
                    assert!(output.gates.iter().all(|g| matches!(
                        g,
                        Gate::h(_)
                            | Gate::s(_)
                            | Gate::sdg(_)
                            | Gate::t(_)
                            | Gate::x(_)
                            | Gate::cnot { .. }
                    )));
                    assert!(
                        circuits_equiv(&input, &output, 2.0 * epsilon),
                        "theta={theta}, q={q}, input={input}"
                    );
                }
            }
        }
    }
}
