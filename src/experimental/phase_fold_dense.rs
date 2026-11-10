//! Exact phase folding with dense symbolic parities.
//!
//! A parity has one bit per possible symbolic variable, packed into `u64`
//! words, plus a constant bit. This is the dense representation evaluated for
//! the POPL 2027 reviewer response.

use rustc_hash::FxHashMap;

use crate::angle::{Angle, AngleContext};
use crate::circuit::{Circuit, Gate, Qubit};
use crate::pass::{Pass, PassState};

use super::{is_region_barrier, LivePhase};

#[derive(Clone, Debug)]
struct DenseParity {
    coefficients: Box<[u64]>,
    negated: bool,
}

impl DenseParity {
    fn fresh(words: usize, variable: usize) -> Self {
        let mut coefficients = vec![0; words].into_boxed_slice();
        coefficients[variable / 64] |= 1u64 << (variable % 64);
        Self {
            coefficients,
            negated: false,
        }
    }

    fn zero(words: usize) -> Self {
        Self {
            coefficients: vec![0; words].into_boxed_slice(),
            negated: false,
        }
    }

    fn xor_assign(&mut self, other: &Self) {
        for (target, control) in self.coefficients.iter_mut().zip(other.coefficients.iter()) {
            *target ^= *control;
        }
        self.negated ^= other.negated;
    }

    fn is_constant(&self) -> bool {
        self.coefficients.iter().all(|word| *word == 0)
    }
}

/// Experimental exact phase folding with dense symbolic parity bit vectors.
pub struct PhaseFoldDense;

impl Pass for PhaseFoldDense {
    fn name(&self) -> &str {
        "Phase folding (experimental dense exact)"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        phase_fold_dense(circuit)
    }

    fn run_with_state(&self, circuit: &Circuit, state: &mut PassState) -> Circuit {
        phase_fold_dense_with_context(circuit, &mut state.angles)
    }
}

/// Return `(variables, u64_words)` for the dense exact representation.
pub fn bitvec_dimensions(circuit: &Circuit) -> (usize, usize) {
    let mut variables = circuit.num_qubits;
    for gate in &circuit.gates {
        if matches!(gate, Gate::h(_) | Gate::ccx { .. }) {
            variables += 1;
        } else if is_region_barrier(gate) {
            variables += circuit.num_qubits;
        }
    }
    (variables, variables.div_ceil(64))
}

/// Free-function form of [`PhaseFoldDense`].
pub fn phase_fold_dense(circuit: &Circuit) -> Circuit {
    phase_fold_dense_with_context(circuit, &mut PassState::default().angles)
}

fn phase_fold_dense_with_context(circuit: &Circuit, ctx: &mut AngleContext) -> Circuit {
    let (variables, words) = bitvec_dimensions(circuit);
    let mut next_variable = 0;
    let mut qubits: Vec<DenseParity> = (0..circuit.num_qubits)
        .map(|_| {
            let parity = DenseParity::fresh(words, next_variable);
            next_variable += 1;
            parity
        })
        .collect();
    let mut live = Vec::new();
    let mut groups: FxHashMap<Box<[u64]>, usize> = FxHashMap::default();
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
                qubits[*q as usize] = DenseParity::fresh(words, next_variable);
                next_variable += 1;
            }
            Gate::x(q) => qubits[*q as usize].negated ^= true,
            Gate::cnot { control, target } => {
                let control = qubits[*control as usize].clone();
                qubits[*target as usize].xor_assign(&control);
            }
            Gate::cz { .. } | Gate::ccz { .. } | Gate::measure { .. } => {}
            Gate::ccx { target, .. } => {
                qubits[*target as usize] = DenseParity::fresh(words, next_variable);
                next_variable += 1;
            }
            Gate::reset(q) => qubits[*q as usize] = DenseParity::zero(words),
            gate if is_region_barrier(gate) => {
                groups.clear();
                for parity in &mut qubits {
                    *parity = DenseParity::fresh(words, next_variable);
                    next_variable += 1;
                }
            }
            _ => unreachable!("all gate variants are covered"),
        }
    }
    debug_assert_eq!(next_variable, variables);

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
    qubits: &[DenseParity],
    qubit: Qubit,
    angle: Angle,
    idx: usize,
    live: &mut Vec<LivePhase>,
    groups: &mut FxHashMap<Box<[u64]>, usize>,
    skip: &mut [bool],
    ctx: &mut AngleContext,
) {
    let parity = &qubits[qubit as usize];
    if parity.is_constant() {
        skip[idx] = true;
        return;
    }
    let complement = parity.negated;
    if let Some(&group) = groups.get(parity.coefficients.as_ref()) {
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
    groups.insert(parity.coefficients.clone(), group);
}
