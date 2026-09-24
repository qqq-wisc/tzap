//! Sparse, line-based PBC export. Materialization is not a linear-time operation.
use super::{PbcCircuit, PbcError, PbcOp, Phase};
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
    /// output Clifford frame. See `docs/pbc.md` for the format.
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
        let used = self.visit_axes(options.max_expansion_cells, |op, phase, factors| {
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
        self.visit_output_frame(
            options.max_expansion_cells - used,
            |is_z, q, phase, factors| {
                let expected = if is_z {
                    super::Pauli::Z
                } else {
                    super::Pauli::X
                };
                if phase == Phase::One && factors.len() == 1 && factors.get(&q) == Some(&expected) {
                    return Ok(());
                }
                let sign = match phase {
                    Phase::One => "1",
                    Phase::MinusOne => "-1",
                    _ => return Err(PbcError::NonHermitianAxis),
                };
                write!(text, "frame {}{q} {sign}", if is_z { 'Z' } else { 'X' }).unwrap();
                let mut sorted: Vec<_> = factors.iter().collect();
                sorted.sort_unstable_by_key(|&(q, _)| q);
                for (qubit, pauli) in sorted {
                    write!(text, " {pauli}{qubit}").unwrap();
                }
                text.push('\n');
                Ok(())
            },
        )?;
        Ok(text)
    }
}
