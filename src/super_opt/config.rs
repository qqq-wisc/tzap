use super::BASE_GATE_SET;
use crate::circuit::GateSet;

/// Bounds and synthesis basis for a minimal unitary representative map.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MurmConfig {
    /// Maximum distinct qubits in a MURM representative.
    pub max_qubits: usize,
    /// Maximum gates in a MURM representative. A representative only replaces a
    /// window strictly larger than itself, so this never needs to exceed
    /// `window_gates - 1` for a given [`crate::super_opt::SuperOpt`].
    pub max_gates: usize,
    /// Enumeration stops independently at this many distinct unitaries per width.
    pub max_entries_per_qubit: usize,
    /// Exact gate basis the MURM may use for representatives.
    pub basis: GateSet,
}

impl Default for MurmConfig {
    fn default() -> Self {
        Self::new(3, 8, 200_000)
    }
}

impl MurmConfig {
    /// See the field docs for the meaning of each bound.
    pub const fn new(max_qubits: usize, max_gates: usize, max_entries_per_qubit: usize) -> Self {
        Self {
            max_qubits,
            max_gates,
            max_entries_per_qubit,
            basis: BASE_GATE_SET,
        }
    }

    pub const fn with_basis(mut self, basis: GateSet) -> Self {
        self.basis = basis;
        self
    }
}
