#![doc = include_str!("../API.md")]

pub mod cancel;
pub mod circuit;
pub mod cnot_min;
pub mod decompose;
pub mod optimize;
pub mod pass;
pub mod pbc;
pub mod phase_fold_rand;
pub mod qasm;
pub mod super_opt;

#[cfg(feature = "python")]
mod python;

#[cfg(test)]
mod bench;
#[cfg(test)]
mod semantics;
#[cfg(test)]
mod unitary;
