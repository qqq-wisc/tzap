//! Exact phase folding with sparse symbolic parities.
//!
//! A parity is the sorted list of variable IDs that occur in it, plus a
//! constant bit. This is the `TZAP-SYMB` representation used in the POPL 2027
//! evaluation.

use std::collections::BTreeMap;

use crate::angle::{Angle, AngleContext};
use crate::circuit::{Circuit, Gate, Qubit};
use crate::pass::{Pass, PassState};

use super::{is_region_barrier, LivePhase};

#[derive(Clone, Debug)]
struct SparseParity {
    vars: Vec<usize>,
    negated: bool,
}

impl SparseParity {
    fn fresh(variable: usize) -> Self {
        Self {
            vars: vec![variable],
            negated: false,
        }
    }

    fn zero() -> Self {
        Self {
            vars: Vec::new(),
            negated: false,
        }
    }

    fn xor_assign(&mut self, other: &Self) {
        let mut result = Vec::with_capacity(self.vars.len() + other.vars.len());
        let (mut left, mut right) = (0, 0);
        while left < self.vars.len() && right < other.vars.len() {
            match self.vars[left].cmp(&other.vars[right]) {
                std::cmp::Ordering::Less => {
                    result.push(self.vars[left]);
                    left += 1;
                }
                std::cmp::Ordering::Greater => {
                    result.push(other.vars[right]);
                    right += 1;
                }
                std::cmp::Ordering::Equal => {
                    left += 1;
                    right += 1;
                }
            }
        }
        result.extend_from_slice(&self.vars[left..]);
        result.extend_from_slice(&other.vars[right..]);
        self.vars = result;
        self.negated ^= other.negated;
    }

    fn is_constant(&self) -> bool {
        self.vars.is_empty()
    }
}

/// Experimental exact phase folding with sparse symbolic parity expressions.
pub struct PhaseFoldSparse;

impl Pass for PhaseFoldSparse {
    fn name(&self) -> &str {
        "Phase folding (experimental sparse exact)"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        phase_fold_sparse(circuit)
    }

    fn run_with_state(&self, circuit: &Circuit, state: &mut PassState) -> Circuit {
        phase_fold_sparse_with_context(circuit, &mut state.angles)
    }
}

/// Free-function form of [`PhaseFoldSparse`].
pub fn phase_fold_sparse(circuit: &Circuit) -> Circuit {
    phase_fold_sparse_with_context(circuit, &mut PassState::default().angles)
}

fn phase_fold_sparse_with_context(circuit: &Circuit, ctx: &mut AngleContext) -> Circuit {
    let mut next_variable = 0;
    let mut qubits: Vec<SparseParity> = (0..circuit.num_qubits)
        .map(|_| {
            let parity = SparseParity::fresh(next_variable);
            next_variable += 1;
            parity
        })
        .collect();
    let mut live = Vec::new();
    let mut groups: BTreeMap<Vec<usize>, usize> = BTreeMap::new();
    let mut skip = vec![false; circuit.gates.len()];
    let mut emit_at: Vec<Option<(Qubit, Angle)>> = vec![None; circuit.gates.len()];

    for (idx, gate) in circuit.gates.iter().enumerate() {
        match gate {
            Gate::t(q) => record_angle(
                &qubits,
                *q,
                Angle::turns(1),
                idx,
                &mut live,
                &mut groups,
                &mut skip,
                ctx,
            ),
            Gate::tdg(q) => record_angle(
                &qubits,
                *q,
                Angle::turns(7),
                idx,
                &mut live,
                &mut groups,
                &mut skip,
                ctx,
            ),
            Gate::s(q) => record_angle(
                &qubits,
                *q,
                Angle::turns(2),
                idx,
                &mut live,
                &mut groups,
                &mut skip,
                ctx,
            ),
            Gate::sdg(q) => record_angle(
                &qubits,
                *q,
                Angle::turns(6),
                idx,
                &mut live,
                &mut groups,
                &mut skip,
                ctx,
            ),
            Gate::z(q) => record_angle(
                &qubits,
                *q,
                Angle::turns(4),
                idx,
                &mut live,
                &mut groups,
                &mut skip,
                ctx,
            ),
            Gate::rz(angle, q) => record_angle(
                &qubits,
                *q,
                angle.clone(),
                idx,
                &mut live,
                &mut groups,
                &mut skip,
                ctx,
            ),
            Gate::h(q) => {
                qubits[*q as usize] = SparseParity::fresh(next_variable);
                next_variable += 1;
            }
            Gate::x(q) => qubits[*q as usize].negated ^= true,
            Gate::cnot { control, target } => {
                let control = qubits[*control as usize].clone();
                qubits[*target as usize].xor_assign(&control);
            }
            Gate::cz { .. } | Gate::ccz { .. } | Gate::measure { .. } => {}
            Gate::ccx { target, .. } => {
                qubits[*target as usize] = SparseParity::fresh(next_variable);
                next_variable += 1;
            }
            Gate::reset(q) => qubits[*q as usize] = SparseParity::zero(),
            gate if is_region_barrier(gate) => {
                groups.clear();
                for parity in &mut qubits {
                    *parity = SparseParity::fresh(next_variable);
                    next_variable += 1;
                }
            }
            _ => unreachable!("all gate variants are covered"),
        }
    }

    for phase in &live {
        if phase.angle.projective_identity() {
            skip[phase.current_idx] = true;
        } else {
            emit_at[phase.current_idx] = Some((phase.qubit, phase.angle.clone()));
        }
    }

    let mut output = circuit.empty_like();
    for (idx, gate) in circuit.gates.clone().into_iter().enumerate() {
        if skip[idx] {
            continue;
        }
        if let Some((qubit, angle)) = &emit_at[idx] {
            angle.emit(&mut output, *qubit);
        } else {
            output.apply(gate);
        }
    }
    output
}

#[allow(clippy::too_many_arguments)]
fn record_angle(
    qubits: &[SparseParity],
    qubit: Qubit,
    angle: Angle,
    idx: usize,
    live: &mut Vec<LivePhase>,
    groups: &mut BTreeMap<Vec<usize>, usize>,
    skip: &mut [bool],
    ctx: &mut AngleContext,
) {
    let parity = &qubits[qubit as usize];
    if parity.is_constant() {
        skip[idx] = true;
        return;
    }
    let complement = parity.negated;
    if let Some(&group) = groups.get(&parity.vars) {
        let sign = if live[group].current_sign == complement {
            1
        } else {
            -1
        };
        if let Ok(merged) = live[group].angle.merge_phases_with(ctx, sign, &angle) {
            skip[live[group].current_idx] = true;
            live[group].angle = merged;
            live[group].current_idx = idx;
            live[group].qubit = qubit;
            live[group].current_sign = complement;
            return;
        }
    }
    let group = live.len();
    live.push(LivePhase {
        angle,
        qubit,
        current_idx: idx,
        current_sign: complement,
    });
    groups.insert(parity.vars.clone(), group);
}
