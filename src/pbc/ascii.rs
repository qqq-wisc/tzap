use super::{PbcCircuit, PbcError, PbcOp};

/// Explicit limits for the diagnostic renderer, which expands shared axes.
#[derive(Clone, Copy, Debug)]
pub struct AsciiOptions {
    pub max_qubits: usize,
    pub max_operations: usize,
    /// Work limit for sparse dependency evaluation (not dense qubit cells).
    pub max_expansion_cells: usize,
}

impl Default for AsciiOptions {
    fn default() -> Self {
        Self {
            max_qubits: 32,
            max_operations: 128,
            max_expansion_cells: 1_000_000,
        }
    }
}

struct Column {
    label: String,
    wires: Vec<String>,
}

impl PbcCircuit {
    /// Render operations left to right, then changed output-frame rows.
    /// Operation columns are joint operations, not independent single-qubit
    /// gates. Frame columns show generator images, not executable gates.
    /// Signs in headings apply to the entire Pauli product.
    pub fn to_ascii(&self) -> Result<String, PbcError> {
        self.to_ascii_with(AsciiOptions::default())
    }

    /// Drawing is bounded, explicit materialization, not part of conversion's
    /// linear-time contract. Oversized circuits return an error, never truncate.
    pub fn to_ascii_with(&self, options: AsciiOptions) -> Result<String, PbcError> {
        if self.num_qubits > options.max_qubits || self.operations.len() > options.max_operations {
            return Err(PbcError::DrawingLimit);
        }
        // One heading plus wires and intervening connector rows (2Q for Q>0).
        let row_count = self
            .num_qubits
            .checked_mul(2)
            .ok_or(PbcError::DrawingLimit)?
            .max(1);
        let mut columns = Vec::new();
        let used = self.visit_axes(options.max_expansion_cells, |op, phase, factors| {
            let identity = if factors.is_empty() { "I" } else { "" };
            let axis = format!("{phase}{identity}");
            let label = match op {
                PbcOp::Rotate { angle, .. } => format!("R({angle},{axis})"),
                PbcOp::Measure {
                    outcome,
                    target: Some(c),
                    ..
                } => format!("M({axis})->{outcome}/c{c}"),
                PbcOp::Measure { outcome, .. } => format!("M({axis})->{outcome}"),
                PbcOp::ConditionalRotate { angle, if_one, .. } => {
                    format!("R({angle},{axis}) if {if_one}=1")
                }
            };
            let mut wires = vec![String::new(); self.num_qubits];
            for (&q, p) in factors {
                wires[q as usize] = format!("[{p}]");
            }
            columns.push(Column { label, wires });
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
                if phase == super::Phase::One
                    && factors.len() == 1
                    && factors.get(&q) == Some(&expected)
                {
                    return Ok(());
                }
                let mut wires = vec![String::new(); self.num_qubits];
                for (&qubit, pauli) in factors {
                    wires[qubit as usize] = format!("[{pauli}]");
                }
                columns.push(Column {
                    label: format!("{}{}->{phase}", if is_z { 'Z' } else { 'X' }, q),
                    wires,
                });
                Ok(())
            },
        )?;
        if columns.len() > options.max_operations {
            return Err(PbcError::DrawingLimit);
        }
        if columns.is_empty() && self.num_qubits == 0 {
            return Ok("(empty PBC circuit; 0 qubits)\n".into());
        }
        let margin = format!("q{}: ", self.num_qubits.saturating_sub(1)).len();
        let mut rows = vec![" ".repeat(margin); row_count];
        for q in 0..self.num_qubits {
            rows[1 + 2 * q] = format!("{:>width$}", format!("q{q}: "), width = margin);
        }
        for column in columns {
            let width = column.label.len().max(5) + 2;
            rows[0].push_str(&format!("{:^width$}", column.label));
            let first = column.wires.iter().position(|s| !s.is_empty());
            let last = column.wires.iter().rposition(|s| !s.is_empty());
            for q in 0..self.num_qubits {
                let connected = first.is_some_and(|a| a <= q) && last.is_some_and(|b| q <= b);
                let symbol = if column.wires[q].is_empty() && connected {
                    "|"
                } else {
                    &column.wires[q]
                };
                rows[1 + 2 * q].push_str(&format!("{symbol:-^width$}"));
                if q + 1 < self.num_qubits {
                    let connector = if connected && last.is_some_and(|b| q < b) {
                        "|"
                    } else {
                        ""
                    };
                    rows[2 + 2 * q].push_str(&format!("{connector:^width$}"));
                }
            }
        }
        Ok(rows
            .iter()
            .map(|s| s.trim_end())
            .collect::<Vec<_>>()
            .join("\n")
            + "\n")
    }
}
