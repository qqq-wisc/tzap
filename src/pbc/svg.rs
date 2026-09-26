//! SVG drawings of PBC circuits in the style of Litinski, "A Game of Surface
//! Codes" (Quantum 3, 128, 2019), Figs. 4 and 6.
//!
//! Each operation is one box spanning its qubits, with a large italic Pauli
//! letter on every wire it covers (𝟙 where it acts trivially inside the span).
//! A tab on the box's right edge shows the angle (π/8, π/4, 3π/8, π/2), with a
//! white strip carrying "−" when the axis is negative. π/8-type rotations are
//! green, π/4 orange, π/2 grey, and measurements blue with a meter tab. Wires
//! are labelled |q₁⟩, |q₂⟩, … (1-indexed, as in the paper).
//!
//! By default, as in the paper's figures, each operation is drawn as early as
//! its commutation allows: one column past the latest earlier operation it
//! anticommutes with, in the first column where its rows are free. Only
//! commuting operations change position, so the drawing denotes the same
//! circuit. The output Clifford is left out, as the paper does.

use std::fmt::Write;

use super::pauli::Factors;
use super::{Pauli, PbcCircuit, PbcError, PbcOp, Phase};

/// Colors sampled from the paper's figures.
const GREEN: &str = "#E3FFA1";
const ORANGE: &str = "#F5BD70";
const BLUE: &str = "#70B3F5";
const GREY: &str = "#DFDFDF";

/// Layout, in SVG user units.
const ROW: f64 = 50.0; // distance between wires
const HALF: f64 = 18.0; // box extent above and below a wire
const BOX_W: f64 = 38.0;
const TAB_W: f64 = 15.0;
const TAB_H: f64 = 30.0;
const GAP: f64 = 12.0; // between columns
const MARGIN_LEFT: f64 = 58.0; // room for the |q⟩ labels
const MARGIN: f64 = 16.0;
const COLUMN: f64 = BOX_W + TAB_W + GAP;

const FONT: &str = "'Latin Modern Math','CMU Serif','STIX Two Math','Cambria Math',\
                    'Times New Roman',serif";

/// Options for [`PbcCircuit::to_svg_with`].
#[derive(Clone, Copy, Debug)]
pub struct SvgOptions {
    /// Draw at most this many operations; the wires of a longer circuit end
    /// in dots, as in the paper's Fig. 6.
    pub max_operations: usize,
    /// Let operations on disjoint qubit ranges share a column.
    pub pack: bool,
    /// Move each operation as early as its commutation allows, as the paper
    /// draws its circuits: one column past the latest earlier operation it
    /// anticommutes with. Only commuting operations change order, so the
    /// drawing denotes the same circuit. Implies packing.
    pub reorder: bool,
    /// Every box spans all qubits (the paper's Fig. 6), instead of just its
    /// support's range (Fig. 4).
    pub full_height: bool,
    /// Draw the output Clifford, when not the identity, as a final box "C".
    /// Off by default: the paper leaves it out.
    pub show_frame: bool,
    /// Label measurements with their classical register and conditional
    /// rotations with their condition (the paper shows neither).
    pub show_targets: bool,
    /// Sparse materialization budget, as for text export.
    pub max_expansion_cells: usize,
}

impl Default for SvgOptions {
    fn default() -> Self {
        Self {
            max_operations: 400,
            pack: true,
            reorder: true,
            full_height: false,
            show_frame: false,
            show_targets: false,
            max_expansion_cells: 16_000_000,
        }
    }
}

/// One drawn box.
struct Item {
    /// Rows covered, inclusive.
    top: usize,
    bottom: usize,
    /// Letter per covered row (index `row - top`).
    letters: Vec<&'static str>,
    /// Non-identity factors, for commutation checks.
    factors: Vec<(u32, Pauli)>,
    kind: Kind,
}

impl Item {
    /// Whether the two operations must keep their order: their axes
    /// anticommute, or one of them is conditional or the frame.
    fn ordered_with(&self, other: &Item) -> bool {
        let barrier = |item: &Item| {
            matches!(
                item.kind,
                Kind::Frame
                    | Kind::Rotation {
                        condition: Some(_),
                        ..
                    }
            )
        };
        if barrier(self) || barrier(other) {
            return true;
        }
        let mut differing = 0;
        let (mut i, mut j) = (0, 0);
        while i < self.factors.len() && j < other.factors.len() {
            let ((qa, pa), (qb, pb)) = (self.factors[i], other.factors[j]);
            match qa.cmp(&qb) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    differing += usize::from(pa != pb);
                    i += 1;
                    j += 1;
                }
            }
        }
        differing % 2 == 1
    }
}

enum Kind {
    /// Angle in units of π/8 (1–4) and whether the axis is negative.
    Rotation {
        eighths: u8,
        negative: bool,
        condition: Option<String>,
    },
    Measure {
        negative: bool,
        target: Option<String>,
    },
    Frame,
}

impl PbcCircuit {
    /// An SVG drawing of the circuit with default options.
    pub fn to_svg(&self) -> Result<String, PbcError> {
        self.to_svg_with(SvgOptions::default())
    }

    /// An SVG drawing of the circuit. Only the first `max_operations`
    /// operations' axes are materialized.
    pub fn to_svg_with(&self, options: SvgOptions) -> Result<String, PbcError> {
        let n = self.num_qubits;
        let shown = self.operations.len().min(options.max_operations);
        let truncated = shown < self.operations.len();
        let mut items = Vec::with_capacity(shown + 1);
        let roots: Vec<_> = self.operations[..shown]
            .iter()
            .map(|op| op.axis().0)
            .collect();
        let used = self
            .arena
            .materialize(
                &roots,
                options.max_expansion_cells,
                |index, phase, factors| {
                    let negative = match phase {
                        Phase::One => false,
                        Phase::MinusOne => true,
                        Phase::I | Phase::MinusI => return Err(PbcError::NonHermitianAxis),
                    };
                    let (top, bottom, letters) = span(factors, n, options.full_height);
                    let kind = match self.operations[index] {
                        PbcOp::Rotate { angle, .. } => {
                            rotation(angle.signed_eighths(), negative, None)
                        }
                        PbcOp::ConditionalRotate { angle, if_one, .. } => rotation(
                            angle.signed_eighths(),
                            negative,
                            Some(format!("if {if_one}")),
                        ),
                        PbcOp::Measure { target, .. } => Kind::Measure {
                            negative,
                            target: target.map(|c| format!("c{c}")),
                        },
                    };
                    // A zero rotation is the identity: nothing to draw.
                    if !matches!(kind, Kind::Rotation { eighths: 0, .. }) {
                        items.push(Item {
                            top,
                            bottom,
                            letters,
                            factors: factors.iter().map(|(&q, &p)| (q, p)).collect(),
                            kind,
                        });
                    }
                    Ok(())
                },
            )?
            .work;
        if options.show_frame && !truncated {
            // The output Clifford acts on the qubits whose images changed.
            let mut rows: Vec<usize> = Vec::new();
            self.visit_output_frame(
                options.max_expansion_cells.saturating_sub(used),
                |is_z, q, phase, factors| {
                    let same = if is_z { Pauli::Z } else { Pauli::X };
                    if !(phase == Phase::One
                        && factors.len() == 1
                        && factors.get(&q) == Some(&same))
                    {
                        rows.push(q as usize);
                        rows.extend(factors.keys().map(|&k| k as usize));
                    }
                    Ok(())
                },
            )?;
            if let (Some(&top), Some(&bottom)) = (rows.iter().min(), rows.iter().max()) {
                let (top, bottom) = if options.full_height {
                    (0, n.saturating_sub(1))
                } else {
                    (top, bottom)
                };
                items.push(Item {
                    top,
                    bottom,
                    letters: Vec::new(),
                    factors: Vec::new(),
                    kind: Kind::Frame,
                });
            }
        }
        Ok(render(n, &items, truncated, &options))
    }
}

fn rotation(k: i8, negative: bool, condition: Option<String>) -> Kind {
    // exp(-i k π/8 (-P)) = exp(-i (-k) π/8 P): fold the angle's sign into the
    // axis sign, so the tab shows |k|. k = 4 is sign-independent up to phase.
    let eighths = k.unsigned_abs();
    let negative = eighths != 4 && (negative != (k < 0));
    Kind::Rotation {
        eighths,
        negative,
        condition,
    }
}

/// Rows covered by an axis and the letter on each.
fn span(factors: &Factors, n: usize, full_height: bool) -> (usize, usize, Vec<&'static str>) {
    let (mut top, mut bottom) = match (factors.keys().next(), factors.keys().next_back()) {
        (Some(&a), Some(&b)) => (a as usize, b as usize),
        // Identity axis: span every wire (a global phase, drawn for honesty).
        _ => (0, n.saturating_sub(1)),
    };
    if full_height {
        (top, bottom) = (0, n.saturating_sub(1));
    }
    let letters = (top..=bottom)
        .map(|row| match factors.get(&(row as u32)) {
            Some(Pauli::X) => "X",
            Some(Pauli::Y) => "Y",
            Some(Pauli::Z) => "Z",
            _ => "𝟙",
        })
        .collect();
    (top, bottom, letters)
}

fn row_y(row: usize) -> f64 {
    MARGIN + 12.0 + HALF + row as f64 * ROW
}

fn render(n: usize, items: &[Item], truncated: bool, options: &SvgOptions) -> String {
    let columns = if options.reorder {
        reordered_columns(n, items)
    } else {
        program_columns(n, items, options.pack)
    };
    draw(n, items, &columns, truncated, options.show_targets)
}

/// Columns in the paper's style: each item goes to the first column, at or
/// after one past the latest earlier item it must stay ordered with, whose
/// rows it spans are free.
fn reordered_columns(n: usize, items: &[Item]) -> Vec<usize> {
    let mut columns: Vec<usize> = Vec::with_capacity(items.len());
    // occupied[c][row]: whether column c already draws something on row.
    let mut occupied: Vec<Vec<bool>> = Vec::new();
    for (index, item) in items.iter().enumerate() {
        let earliest = items[..index]
            .iter()
            .zip(&columns)
            .filter(|(earlier, _)| earlier.ordered_with(item))
            .map(|(_, &c)| c + 1)
            .max()
            .unwrap_or(0);
        let mut column = earliest;
        while column < occupied.len() && occupied[column][item.top..=item.bottom].iter().any(|&o| o)
        {
            column += 1;
        }
        if column == occupied.len() {
            occupied.push(vec![false; n.max(1)]);
        }
        for cell in &mut occupied[column][item.top..=item.bottom] {
            *cell = true;
        }
        columns.push(column);
    }
    columns
}

/// Columns in program order: with packing, one past the latest column whose
/// rows overlap; otherwise one column per item.
fn program_columns(n: usize, items: &[Item], pack: bool) -> Vec<usize> {
    // Column assignment: one past the latest column whose occupied rows
    // overlap (with packing), or strictly sequential.
    let mut next_free = vec![0usize; n.max(1)];
    let mut columns = Vec::with_capacity(items.len());
    let mut sequential = 0;
    for item in items {
        let column = if matches!(item.kind, Kind::Frame) {
            // The output Clifford comes after everything.
            sequential
        } else if pack {
            (item.top..=item.bottom)
                .map(|r| next_free[r])
                .max()
                .unwrap_or(0)
        } else {
            sequential
        };
        for free in &mut next_free[item.top..=item.bottom] {
            *free = column + 1;
        }
        sequential = sequential.max(column + 1);
        columns.push(column);
    }
    columns
}

fn draw(
    n: usize,
    items: &[Item],
    columns: &[usize],
    truncated: bool,
    show_targets: bool,
) -> String {
    let width_columns = columns.iter().map(|c| c + 1).max().unwrap_or(0);
    let wires_end = MARGIN_LEFT + width_columns as f64 * COLUMN + GAP;
    let width = wires_end + if truncated { 44.0 } else { 0.0 } + MARGIN;
    let height = row_y(n.saturating_sub(1)) + HALF + if show_targets { 22.0 } else { 4.0 } + MARGIN;

    let mut svg = String::new();
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0}" height="{height:.0}" viewBox="0 0 {width:.1} {height:.1}" font-family="{FONT}">"#
    )
    .unwrap();
    writeln!(svg, r#"<rect width="100%" height="100%" fill="white"/>"#).unwrap();

    // Wires and their labels.
    for q in 0..n {
        let y = row_y(q);
        writeln!(
            svg,
            r#"<text x="{:.1}" y="{:.1}" font-size="22" text-anchor="end" dominant-baseline="central">|<tspan font-style="italic">q</tspan><tspan font-size="14" dy="5">{}</tspan><tspan dy="-5">⟩</tspan></text>"#,
            MARGIN_LEFT - 6.0,
            y,
            q + 1
        )
        .unwrap();
        writeln!(
            svg,
            r#"<line x1="{:.1}" y1="{y:.1}" x2="{wires_end:.1}" y2="{y:.1}" stroke="black" stroke-width="1.3"/>"#,
            MARGIN_LEFT
        )
        .unwrap();
        if truncated {
            writeln!(
                svg,
                r#"<line x1="{:.1}" y1="{y:.1}" x2="{:.1}" y2="{y:.1}" stroke="black" stroke-width="1.6" stroke-dasharray="1.6 4" stroke-linecap="round"/>"#,
                wires_end + 4.0,
                wires_end + 40.0
            )
            .unwrap();
        }
    }

    for (item, &column) in items.iter().zip(columns) {
        let x = MARGIN_LEFT + GAP + column as f64 * COLUMN;
        let top = row_y(item.top) - HALF;
        let bottom = row_y(item.bottom) + HALF;
        let middle = (top + bottom) / 2.0;
        let (fill, dashed) = match &item.kind {
            Kind::Rotation {
                eighths, condition, ..
            } => (
                match eighths {
                    2 => ORANGE,
                    4 => GREY,
                    _ => GREEN,
                },
                condition.is_some(),
            ),
            Kind::Measure { .. } => (BLUE, false),
            Kind::Frame => ("white", false),
        };
        let dash = if dashed {
            r#" stroke-dasharray="4 2""#
        } else {
            ""
        };

        // The tab is drawn first so the box's edge overlaps it.
        match &item.kind {
            Kind::Rotation {
                eighths, negative, ..
            } => tab(&mut svg, x + BOX_W, middle, fill, *negative, Some(*eighths)),
            Kind::Measure { negative, .. } => {
                tab(&mut svg, x + BOX_W, middle, fill, *negative, None)
            }
            Kind::Frame => {}
        }
        writeln!(
            svg,
            r#"<rect x="{x:.1}" y="{top:.1}" width="{BOX_W}" height="{:.1}" fill="{fill}" stroke="black" stroke-width="1.3"{dash}/>"#,
            bottom - top
        )
        .unwrap();
        match &item.kind {
            Kind::Frame => {
                writeln!(
                    svg,
                    r#"<text x="{:.1}" y="{middle:.1}" font-size="26" text-anchor="middle" dominant-baseline="central" font-style="italic">C</text>"#,
                    x + BOX_W / 2.0
                )
                .unwrap();
            }
            _ => {
                for (i, letter) in item.letters.iter().enumerate() {
                    let style = if *letter == "𝟙" {
                        ""
                    } else {
                        r#" font-style="italic""#
                    };
                    writeln!(
                        svg,
                        r#"<text x="{:.1}" y="{:.1}" font-size="26" text-anchor="middle" dominant-baseline="central"{style}>{letter}</text>"#,
                        x + BOX_W / 2.0,
                        row_y(item.top + i)
                    )
                    .unwrap();
                }
            }
        }
        // Small annotations below the box: measurement target or condition.
        let note = match &item.kind {
            _ if !show_targets => None,
            Kind::Measure { target, .. } => target.as_deref(),
            Kind::Rotation { condition, .. } => condition.as_deref(),
            Kind::Frame => None,
        };
        if let Some(note) = note {
            writeln!(
                svg,
                r#"<text x="{:.1}" y="{:.1}" font-size="11" text-anchor="middle" font-style="italic">{note}</text>"#,
                x + (BOX_W + TAB_W) / 2.0,
                bottom + 13.0
            )
            .unwrap();
        }
    }
    svg.push_str("</svg>\n");
    svg
}

/// The tab on a box's right edge: an angle fraction, or a meter for a
/// measurement; a white strip with "−" on top when negative.
fn tab(svg: &mut String, x: f64, middle: f64, fill: &str, negative: bool, eighths: Option<u8>) {
    let (top, h) = (middle - TAB_H / 2.0, TAB_H);
    // Measurement tabs are rounder, as in the paper.
    let r = if eighths.is_some() { 4.5 } else { 7.0 };
    // Rounded on the right, square where it meets the box.
    let path = format!(
        "M{x:.1},{top:.1} h{:.1} a{r},{r} 0 0 1 {r},{r} v{:.1} a{r},{r} 0 0 1 -{r},{r} h-{:.1} z",
        TAB_W - r,
        h - 2.0 * r,
        TAB_W - r
    );
    writeln!(
        svg,
        r#"<path d="{path}" fill="{fill}" stroke="black" stroke-width="1.1"/>"#
    )
    .unwrap();
    let mut content_top = top;
    if negative {
        let strip = 8.0;
        let strip_path = format!(
            "M{x:.1},{top:.1} h{:.1} a{r},{r} 0 0 1 {r},{r} v{:.1} h-{TAB_W} z",
            TAB_W - r,
            strip - r
        );
        writeln!(
            svg,
            r#"<path d="{strip_path}" fill="white" stroke="black" stroke-width="1.1"/>"#
        )
        .unwrap();
        writeln!(
            svg,
            r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="black" stroke-width="1.4"/>"#,
            x + 4.5,
            top + strip / 2.0,
            x + TAB_W - 4.5,
            top + strip / 2.0
        )
        .unwrap();
        content_top = top + strip;
    }
    let cx = x + TAB_W / 2.0;
    let cy = (content_top + top + h) / 2.0;
    match eighths {
        Some(k) => {
            // A stacked fraction: numerator over a rule over denominator.
            let (numerator, denominator) = match k {
                1 => ("π", "8"),
                2 => ("π", "4"),
                3 => ("3π", "8"),
                _ => ("π", "2"),
            };
            let size = if negative { 8.5 } else { 9.5 };
            writeln!(
                svg,
                r#"<text x="{cx:.1}" y="{:.1}" font-size="{size}" text-anchor="middle">{numerator}</text><line x1="{:.1}" y1="{cy:.1}" x2="{:.1}" y2="{cy:.1}" stroke="black" stroke-width="0.7"/><text x="{cx:.1}" y="{:.1}" font-size="{size}" text-anchor="middle">{denominator}</text>"#,
                cy - 2.0,
                cx - 4.5,
                cx + 4.5,
                cy + size - 1.0
            )
            .unwrap();
        }
        None => {
            // A meter as in the paper: a ")" arc crossed by a needle running
            // from upper left to lower right.
            writeln!(
                svg,
                r#"<path d="M{:.1},{:.1} a8,8 0 0 1 0,14" fill="none" stroke="black" stroke-width="1.2"/><line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="black" stroke-width="1.2"/>"#,
                cx - 3.5,
                cy - 7.0,
                cx - 3.5,
                cy - 1.5,
                cx + 6.0,
                cy + 4.0
            )
            .unwrap();
        }
    }
}

#[cfg(test)]
mod tests;
