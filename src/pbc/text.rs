//! Sparse, line-based PBC export. Materialization is not a linear-time operation.
use super::{PbcCircuit, PbcError, PbcOp, Phase};
use crate::circuit::qubit_operands;
use std::fmt::Write;

#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug)]
pub struct TextOptions {
    /// Sparse expansion work budget: roots, reachable nodes, and factors
    /// copied, multiplied, or emitted. Unrelated nodes and identity wires cost nothing.
    pub max_expansion_cells: usize,
}

impl Default for TextOptions {
    fn default() -> Self {
        Self {
            max_expansion_cells: 16_000_000,
        }
    }
}

impl PbcCircuit {
    /// Export with all quantum and classical outputs preserved, including the
    /// named Clifford suffix. See `docs/pbc.md` for the format.
    pub fn to_text(&self) -> Result<String, PbcError> {
        self.to_text_with(TextOptions::default())
    }

    /// Materializes shared axes once, within the supplied cell budget. The text
    /// format supports measurements into classical registers, not hidden
    /// outcomes or feed-forward. Reject unsupported operations instead of losing them.
    pub fn to_text_with(&self, options: TextOptions) -> Result<String, PbcError> {
        let exportable = |op: &PbcOp| {
            matches!(
                op,
                PbcOp::Rotate { .. }
                    | PbcOp::Measure {
                        target: Some(_),
                        ..
                    }
            )
        };
        if let Some(index) = self.operations.iter().position(|op| !exportable(op)) {
            return Err(PbcError::UnsupportedTextOperation { index });
        }
        let mut text = format!("qubits {}\nregisters {}\n", self.num_qubits, self.num_cbits);
        self.visit_axes(options.max_expansion_cells, |op, phase, factors| {
            let sign = match phase {
                Phase::One => "1",
                Phase::MinusOne => "-1",
                Phase::I | Phase::MinusI => return Err(PbcError::NonHermitianAxis),
            };
            match op {
                PbcOp::Rotate { angle, .. } => write!(text, "r {} {sign}", angle.signed_eighths()),
                _ => write!(text, "m {sign}"),
            }
            .unwrap();
            for (q, p) in factors {
                write!(text, " {p}{q}").unwrap();
            }
            if let PbcOp::Measure {
                target: Some(c), ..
            } = op
            {
                write!(text, " -> c{c}").unwrap();
            }
            text.push('\n');
            Ok(())
        })?;
        // Named gates form a trailing block, executed after all r/m records.
        for gate in &self.output_cliffords {
            let (n, qs) = qubit_operands(gate);
            text.push_str(gate.kind().name());
            for q in &qs[..n] {
                write!(text, " {q}").unwrap();
            }
            text.push('\n');
        }
        Ok(text)
    }
}
