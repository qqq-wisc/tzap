use std::fmt;
use std::ops::Mul;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::circuit::Qubit;

use super::PbcError;

mod materialize;
pub(crate) use materialize::Factors;

/// Exact phase multiplying a Pauli expression, distinct from a rotation angle.
/// Variants are ordered by their power of i.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Phase {
    #[default]
    One,
    I,
    MinusOne,
    MinusI,
}

impl Mul for Phase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        use Phase::*;
        [One, I, MinusOne, MinusI][(self as usize + rhs as usize) % 4]
    }
}

impl fmt::Display for Phase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::One => "+",
            Self::I => "+i",
            Self::MinusOne => "-",
            Self::MinusI => "-i",
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pauli {
    I,
    X,
    Y,
    Z,
}

impl Pauli {
    /// Single-qubit product `self * other` as a phase and a Pauli.
    pub(crate) fn times(self, other: Self) -> (Phase, Self) {
        use Pauli::*;
        match (self, other) {
            (I, p) | (p, I) => (Phase::One, p),
            (X, X) | (Y, Y) | (Z, Z) => (Phase::One, I),
            (X, Y) => (Phase::I, Z),
            (Y, X) => (Phase::MinusI, Z),
            (Y, Z) => (Phase::I, X),
            (Z, Y) => (Phase::MinusI, X),
            (Z, X) => (Phase::I, Y),
            (X, Z) => (Phase::MinusI, Y),
        }
    }
}

impl fmt::Display for Pauli {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::I => "I",
            Self::X => "X",
            Self::Y => "Y",
            Self::Z => "Z",
        })
    }
}

/// Immutable reference into one circuit's arena. Equality is structural, not
/// algebraic. References from different circuits cannot be mixed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PauliRef {
    pub(crate) owner: u64,
    pub(crate) node: usize,
    pub(crate) phase: Phase,
}

impl PauliRef {
    pub fn node_index(self) -> usize {
        self.node
    }
    pub fn phase(self) -> Phase {
        self.phase
    }
    pub fn scaled(self, phase: Phase) -> Self {
        Self {
            phase: self.phase * phase,
            ..self
        }
    }
}

/// A Pauli reference known to denote a Hermitian operator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PauliAxis(pub(crate) PauliRef);

impl PauliAxis {
    pub fn as_ref(self) -> PauliRef {
        self.0
    }
    pub fn negated(self) -> Self {
        Self(self.0.scaled(Phase::MinusOne))
    }
}

/// Products are ordered and refer only to older nodes. The flat arena preserves
/// sharing and avoids recursive copying or destruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PauliNode {
    Identity,
    Single { qubit: Qubit, pauli: Pauli },
    Product(PauliRef, PauliRef),
}

#[derive(Debug)]
pub(crate) struct PauliArena {
    pub owner: u64,
    pub nodes: Vec<PauliNode>,
}

impl PauliArena {
    pub fn new() -> Self {
        static NEXT_OWNER: AtomicU64 = AtomicU64::new(0);
        let owner = NEXT_OWNER
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_add(1))
            .expect("PBC arena identifiers exhausted");
        Self {
            owner,
            nodes: vec![PauliNode::Identity],
        }
    }

    pub fn check(&self, axis: PauliRef) -> Result<(), PbcError> {
        if axis.owner != self.owner || axis.node >= self.nodes.len() {
            return Err(PbcError::ForeignPauli);
        }
        Ok(())
    }

    pub fn reference(&self, node: usize) -> PauliRef {
        PauliRef {
            owner: self.owner,
            node,
            phase: Phase::One,
        }
    }

    pub fn push(&mut self, node: PauliNode) -> PauliRef {
        let reference = self.reference(self.nodes.len());
        self.nodes.push(node);
        reference
    }
}

/// Explicit dense Pauli factors, indexed by qubit. Materialization is intended
/// for inspection of small circuits, not part of linear-time conversion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExpandedPauli {
    pub phase: Phase,
    pub factors: Vec<Pauli>,
}
