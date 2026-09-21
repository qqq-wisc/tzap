use std::fmt;

use crate::circuit::Qubit;

/// Errors returned by [`crate::super_opt::SuperOpt::new`] and
/// [`crate::super_opt::SuperOpt::run`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SuperOptError {
    /// `window_gates` was zero.
    ZeroWindowGates,
    /// The requested [`crate::super_opt::MurmConfig`] is unusable.
    InvalidMurmConfig { reason: String },
    /// A gate in the input circuit references a qubit outside its declared range.
    InvalidQubit {
        gate_index: usize,
        qubit: Qubit,
        num_qubits: usize,
    },
    /// A window's dense unitary matrix would be too large to construct.
    MatrixTooLarge { num_qubits: usize },
}

impl fmt::Display for SuperOptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroWindowGates => write!(f, "window_gates must be greater than zero"),
            Self::InvalidMurmConfig { reason } => {
                write!(f, "invalid SuperOpt MURM config: {reason}")
            }
            Self::InvalidQubit {
                gate_index,
                qubit,
                num_qubits,
            } => write!(
                f,
                "gate {gate_index} references qubit {qubit}, but the circuit has {num_qubits} qubits"
            ),
            Self::MatrixTooLarge { num_qubits } => {
                write!(f, "a dense matrix for {num_qubits} qubits is too large")
            }
        }
    }
}

impl std::error::Error for SuperOptError {}
