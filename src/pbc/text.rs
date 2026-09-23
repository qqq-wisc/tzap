//! Sparse, line-based PBC export. Materialization is not a linear-time operation.
use super::{Pauli, PbcCircuit, PbcError, PbcOp, Phase};
use crate::circuit::Gate;
use std::fmt::Write;

#[cfg(test)]
mod tests;

/// Whether exported quantum outputs remain observable.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OutputSemantics {
    /// Preserve all quantum outputs, lowering remaining Cliffords to rotations.
    #[default]
    Quantum,
    /// Discard quantum outputs; retain only the classical result distribution.
    Classical,
}

#[derive(Clone, Copy, Debug)]
pub struct TextOptions {
    pub output_semantics: OutputSemantics,
    /// Bound on expanded arena nodes times qubit count (at least one per node).
    pub max_expansion_cells: usize,
}

impl Default for TextOptions {
    fn default() -> Self {
        Self {
            output_semantics: OutputSemantics::Quantum,
            max_expansion_cells: 16_000_000,
        }
    }
}

impl PbcCircuit {
    /// Export with quantum outputs preserved. Use `to_text_with` to explicitly
    /// discard them. Registers are c0..c(N-1); omitted Pauli factors are identity.
    /// Angles are signed integer multiples of pi/8, with exp(-i angle P).
    /// Empty factor lists represent identity. No measurement outcomes are sampled.
    pub fn to_text(&self) -> Result<String, PbcError> {
        self.to_text_with(TextOptions::default())
    }

    /// Materializes shared axes once, within the supplied cell budget. The text
    /// format supports measurements into classical registers, not hidden
    /// outcomes or feed-forward. Reject unsupported operations instead of losing them.
    pub fn to_text_with(&self, options: TextOptions) -> Result<String, PbcError> {
        for (index, op) in self.operations.iter().enumerate() {
            match op {
                PbcOp::Rotate { .. } => (),
                PbcOp::Measure {
                    target: Some(_), ..
                } => (),
                _ => return Err(PbcError::UnsupportedTextOperation { index }),
            }
        }
        let values = match self.operations.iter().map(|op| op.axis().0.node).max() {
            Some(last) => self
                .arena
                .expand(self.num_qubits, last, options.max_expansion_cells)?,
            None => Vec::new(),
        };
        let mut text = format!("qubits {}\nregisters {}\n", self.num_qubits, self.num_cbits);
        for op in &self.operations {
            let reference = op.axis().as_ref();
            let value = &values[reference.node];
            let sign = match value.phase.times(reference.phase) {
                Phase::One => "1",
                Phase::MinusOne => "-1",
                _ => return Err(PbcError::NonHermitianAxis),
            };
            match op {
                PbcOp::Rotate { angle, .. } => {
                    let k = i16::from(angle.eighths());
                    let k = if k > 4 { k - 8 } else { k };
                    write!(text, "r {k} {sign}").unwrap();
                }
                PbcOp::Measure { .. } => write!(text, "m {sign}").unwrap(),
                _ => unreachable!("validated above"),
            }
            for (q, p) in value.factors.iter().enumerate() {
                if *p != Pauli::I {
                    write!(text, " {p}{q}").unwrap();
                }
            }
            if let PbcOp::Measure {
                target: Some(c), ..
            } = op
            {
                write!(text, " -> c{c}").unwrap();
            }
            text.push('\n');
        }
        if options.output_semantics == OutputSemantics::Quantum && !self.output_cliffords.is_empty()
        {
            for gate in &self.output_cliffords {
                write_clifford(&mut text, gate);
            }
        }
        Ok(text)
    }
}

/// Constant-size exact identities, up to global phase. In time order:
/// H = Rz(pi/4) Rx(pi/4) Rz(pi/4);
/// controlled-P = R_Zc(pi/4) R_Pt(pi/4) R_ZcPt(-pi/4).
fn write_clifford(text: &mut String, gate: &Gate) {
    use Pauli::{X, Z};
    let mut rotation = |k: i8, factors: &[(u32, Pauli)]| {
        write!(text, "r {k} 1").unwrap();
        for (q, p) in factors {
            write!(text, " {p}{q}").unwrap();
        }
        text.push('\n');
    };
    match *gate {
        Gate::h(q) => {
            rotation(2, &[(q, Z)]);
            rotation(2, &[(q, X)]);
            rotation(2, &[(q, Z)]);
        }
        Gate::x(q) => rotation(4, &[(q, X)]),
        Gate::z(q) => rotation(4, &[(q, Z)]),
        Gate::s(q) => rotation(2, &[(q, Z)]),
        Gate::sdg(q) => rotation(-2, &[(q, Z)]),
        Gate::cnot { control, target } | Gate::cz { control, target } => {
            let p = if matches!(gate, Gate::cnot { .. }) {
                X
            } else {
                Z
            };
            rotation(2, &[(control, Z)]);
            rotation(2, &[(target, p)]);
            let mut factors = [(control, Z), (target, p)];
            factors.sort_unstable_by_key(|&(q, _)| q);
            rotation(-2, &factors);
        }
        _ => unreachable!("suffix construction validates Cliffords"),
    }
}
