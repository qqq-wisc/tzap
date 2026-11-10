//! Experimental passes retained for paper evaluation and reproduction.
//!
//! These passes are callable through the library and by their hidden
//! `--passes` names, but are deliberately not part of the advertised CLI pass
//! list or any default optimization pipeline.

mod phase_fold_dense;
mod phase_fold_sparse;

pub use phase_fold_dense::{bitvec_dimensions, phase_fold_dense, PhaseFoldDense};
pub use phase_fold_sparse::{phase_fold_sparse, PhaseFoldSparse};
pub use PhaseFoldDense as PhaseFoldBitVecExact;
pub use PhaseFoldSparse as PhaseFoldGlobalExpr;

use crate::circuit::{Gate, Qubit};

/// Gates for which these paper-era linear parity analyses have no transfer
/// rule. Each pass conservatively starts a new analysis region after one.
fn is_region_barrier(gate: &Gate) -> bool {
    matches!(
        gate,
        Gate::p(..)
            | Gate::y(_)
            | Gate::sx(_)
            | Gate::rx(..)
            | Gate::ry(..)
            | Gate::swap(..)
            | Gate::cy { .. }
            | Gate::cp { .. }
            | Gate::crx { .. }
            | Gate::cry { .. }
            | Gate::crz { .. }
            | Gate::ch { .. }
            | Gate::cswap { .. }
    )
}

#[derive(Clone)]
struct LivePhase {
    angle: crate::angle::Angle,
    qubit: Qubit,
    current_idx: usize,
    current_sign: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::{Circuit, Gate};
    use crate::unitary::circuits_equiv;

    fn paper_style_fixture() -> Circuit {
        let mut circuit = Circuit::new(3);
        circuit.apply(Gate::t(0));
        circuit.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        circuit.apply(Gate::tdg(1));
        circuit.apply(Gate::h(2));
        circuit.apply(Gate::cnot {
            control: 2,
            target: 0,
        });
        circuit.apply(Gate::t(0));
        circuit.apply(Gate::x(0));
        circuit.apply(Gate::tdg(0));
        circuit
    }

    #[test]
    fn sparse_and_dense_emit_the_same_paper_style_circuit() {
        let circuit = paper_style_fixture();
        let sparse = phase_fold_sparse(&circuit);
        let dense = phase_fold_dense(&circuit);
        assert_eq!(sparse.gates(), dense.gates());
        assert!(circuits_equiv(&circuit, &sparse, 1e-10));
    }

    #[test]
    fn sparse_and_dense_preserve_register_counts() {
        let mut circuit = Circuit::with_cbits(2, 3);
        circuit.apply(Gate::t(0));
        circuit.apply(Gate::measure { qubit: 0, cbit: 2 });
        circuit.apply(Gate::tdg(0));
        for output in [phase_fold_sparse(&circuit), phase_fold_dense(&circuit)] {
            assert_eq!(output.num_qubits, 2);
            assert_eq!(output.num_cbits, 3);
        }
    }

    #[test]
    fn native_gate_barriers_keep_the_variants_in_sync() {
        let mut circuit = paper_style_fixture();
        circuit.apply(Gate::swap(0, 1));
        circuit.apply(Gate::t(0));
        circuit.apply(Gate::t(0));
        let sparse = phase_fold_sparse(&circuit);
        let dense = phase_fold_dense(&circuit);
        assert_eq!(sparse.gates(), dense.gates());
        assert!(circuits_equiv(&circuit, &sparse, 1e-10));
    }
}
