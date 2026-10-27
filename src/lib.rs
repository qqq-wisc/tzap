#![doc = include_str!("../docs/API.md")]

pub mod angle;
mod angle_export;
mod angle_expr;
pub mod angle_stats;
pub mod symbolic_angle;

pub mod cancel;
pub mod circuit;
pub mod cnot_min;
pub mod decompose;
pub mod optimize;
pub mod pass;
pub mod pbc;
pub mod phase_fold_pauli;
pub mod phase_fold_rand;
pub mod qasm;
pub mod super_opt;
mod synthesis;

#[cfg(feature = "python")]
mod python;

#[cfg(test)]
mod bench;
#[cfg(test)]
mod semantics;
#[cfg(test)]
mod unitary;

#[cfg(test)]
mod native_gate_tests;
