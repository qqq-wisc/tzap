use super::{PbcCircuit, PbcError, PbcOp};
use crate::circuit::{Gate, GateKind, qubit_operands};

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

fn suffix_column(gate: &Gate, num_qubits: usize) -> Column {
    let mut wires = vec![String::new(); num_qubits];
    let label = match *gate {
        Gate::cnot { control, target } => {
            wires[control as usize] = "@".into();
            wires[target as usize] = "[X]".into();
            "CX".into()
        }
        Gate::cz { control, target } => {
            wires[control as usize] = "@".into();
            wires[target as usize] = "@".into();
            "CZ".into()
        }
        _ => {
            let label = match gate.kind() {
                GateKind::Sdg => "Sdg".into(),
                kind => kind.name().to_uppercase(),
            };
            let (_, [q, ..]) = qubit_operands(gate);
            wires[q as usize] = format!("[{label}]");
            label
        }
    };
    Column { label, wires }
}

impl PbcCircuit {
    /// Render operations left to right. Each column is one joint operation,
    /// not independent single-qubit gates. Signs in headings apply to the
    /// entire Pauli product. Classical destinations appear beside outcome IDs.
    pub fn to_ascii(&self) -> Result<String, PbcError> {
        self.to_ascii_with(AsciiOptions::default())
    }

    /// Drawing is bounded, explicit materialization, not part of conversion's
    /// linear-time contract. Oversized circuits return an error, never truncate.
    pub fn to_ascii_with(&self, options: AsciiOptions) -> Result<String, PbcError> {
        if self.num_qubits > options.max_qubits
            || self
                .operations
                .len()
                .saturating_add(self.output_cliffords.len())
                > options.max_operations
        {
            return Err(PbcError::DrawingLimit);
        }
        // One heading plus wires and intervening connector rows (2Q for Q>0).
        let row_count = self
            .num_qubits
            .checked_mul(2)
            .ok_or(PbcError::DrawingLimit)?
            .max(1);
        let mut columns = Vec::new();
        self.visit_axes(options.max_expansion_cells, |op, phase, factors| {
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
        if !self.output_cliffords.is_empty() {
            columns.push(Column {
                label: "suffix".into(),
                wires: vec!["|".into(); self.num_qubits],
            });
        }
        columns.extend(
            self.output_cliffords
                .iter()
                .map(|gate| suffix_column(gate, self.num_qubits)),
        );
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
