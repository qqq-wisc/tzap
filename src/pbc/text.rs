//! Sparse, line-based PBC export. Materialization is not a linear-time operation.
use super::{PbcCircuit, PbcError, PbcOp, Phase};
use crate::circuit::Gate;
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
    /// named Clifford suffix. Registers are c0..c(N-1); omitted factors are identity.
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
        let mut text = format!("qubits {}\nregisters {}\n", self.num_qubits, self.num_cbits);
        let roots: Vec<_> = self
            .operations
            .iter()
            .map(|op| op.axis().as_ref())
            .collect();
        let mut operations = self.operations.iter();
        self.arena
            .materialize(&roots, options.max_expansion_cells, |reference, value| {
                let op = operations.next().unwrap();
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
                for (q, p) in value.sorted_factors() {
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
        for gate in &self.output_cliffords {
            write_clifford(&mut text, gate);
        }
        Ok(text)
    }
}

/// Named gates form a trailing block, executed after all rotation/measurement records.
fn write_clifford(text: &mut String, gate: &Gate) {
    match *gate {
        Gate::h(q) => writeln!(text, "h {q}").unwrap(),
        Gate::x(q) => writeln!(text, "x {q}").unwrap(),
        Gate::z(q) => writeln!(text, "z {q}").unwrap(),
        Gate::s(q) => writeln!(text, "s {q}").unwrap(),
        Gate::sdg(q) => writeln!(text, "sdg {q}").unwrap(),
        Gate::cnot { control, target } => writeln!(text, "cx {control} {target}").unwrap(),
        Gate::cz { control, target } => writeln!(text, "cz {control} {target}").unwrap(),
        _ => unreachable!("suffix construction validates Cliffords"),
    }
}
