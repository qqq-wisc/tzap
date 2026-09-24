use super::{Pauli, PauliArena, PauliNode, PauliRef, Phase};
use crate::circuit::Gate;

/// Complete Clifford action, up to global phase. A row is C† Xq C or C† Zq C
/// for the postponed output Clifford C. Pauli expressions share the circuit's
/// arena, so conversion never expands a row into a dense qubit array.
#[derive(Debug)]
pub struct CliffordFrame {
    pub(crate) x: Vec<PauliRef>,
    pub(crate) z: Vec<PauliRef>,
}

impl CliffordFrame {
    pub(crate) fn identity(arena: &mut PauliArena, n: usize) -> Self {
        let mut x = Vec::with_capacity(n);
        let mut z = Vec::with_capacity(n);
        for q in 0..n {
            let q = u32::try_from(q).expect("PBC qubit count validated by caller");
            x.push(arena.push(PauliNode::Single {
                qubit: q,
                pauli: Pauli::X,
            }));
            z.push(arena.push(PauliNode::Single {
                qubit: q,
                pauli: Pauli::Z,
            }));
        }
        Self { x, z }
    }

    pub fn x(&self, q: u32) -> PauliRef {
        self.x[q as usize]
    }
    pub fn z(&self, q: u32) -> PauliRef {
        self.z[q as usize]
    }
    pub fn len(&self) -> usize {
        self.x.len()
    }
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
    pub(crate) fn roots(&self) -> impl Iterator<Item = PauliRef> + '_ {
        self.x.iter().copied().chain(self.z.iter().copied())
    }

    pub(crate) fn append(&mut self, arena: &mut PauliArena, gate: &Gate) {
        let product = |arena: &mut PauliArena, a, b| arena.push(PauliNode::Product(a, b));
        match *gate {
            Gate::h(q) => {
                std::mem::swap(&mut self.x[q as usize], &mut self.z[q as usize]);
            }
            Gate::x(q) => self.z[q as usize] = self.z[q as usize].scaled(Phase::MinusOne),
            Gate::z(q) => self.x[q as usize] = self.x[q as usize].scaled(Phase::MinusOne),
            Gate::s(q) | Gate::sdg(q) => {
                let i = q as usize;
                let phase = if matches!(gate, Gate::s(_)) {
                    Phase::MinusI
                } else {
                    Phase::I
                };
                self.x[i] = product(arena, self.x[i], self.z[i]).scaled(phase);
            }
            Gate::cnot { control, target } => {
                let c = control as usize;
                let t = target as usize;
                self.x[c] = product(arena, self.x[c], self.x[t]);
                self.z[t] = product(arena, self.z[c], self.z[t]);
            }
            Gate::cz { control, target } => {
                let c = control as usize;
                let t = target as usize;
                self.x[c] = product(arena, self.x[c], self.z[t]);
                self.x[t] = product(arena, self.z[c], self.x[t]);
            }
            _ => unreachable!("validated output Clifford"),
        }
    }
}
