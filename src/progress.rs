//! Terminal progress-bar rendering: number formatting, colored bars, live
//! redrawn progress boxes, and the final result banner. Pure formatting over
//! plain numbers — no dependency on `Circuit`, passes, or CLI options.
//!
//! Everything that writes to the terminal hangs off [`Ui`], which owns the
//! decision of whether color and in-place redraws are allowed at all (see
//! `ui.rs`). Without a live terminal the progress boxes are not drawn at all —
//! a redrawn box is meaningless in a pipe or a log file — so such a run stays
//! silent between its parse line and its final result, which is plain text
//! and always printed.

use std::io::{self, Write};

use crate::ui::Ui;

pub(crate) fn fmt_num<N: std::fmt::Display>(n: N) -> String {
    let s = n.to_string();
    let is_negative = s.starts_with('-');
    let num_part = if is_negative { &s[1..] } else { &s[..] };
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    if is_negative {
        result.push('-');
    }

    let rem = num_part.len() % 3;
    for (i, c) in num_part.chars().enumerate() {
        if i > 0 && i % 3 == rem {
            result.push(',');
        }
        result.push(c);
    }
    result
}

/// Format a byte count with the largest unit that keeps it above 1 — a
/// fixed `MB` rendered every small benchmark as a flat "0.0 MB", which reads
/// like an empty file.
pub(crate) fn fmt_size(bytes: u64) -> String {
    const KB: f64 = 1024.0;
    const MB: f64 = KB * 1024.0;
    let size = bytes as f64;
    if size >= MB {
        format!("{:.1} MB", size / MB)
    } else if size >= KB {
        format!("{:.1} KB", size / KB)
    } else {
        format!("{bytes} B")
    }
}

/// Percentage reduction from `before` to `after` (0.0 when `before` is 0).
fn pct(before: usize, after: usize) -> f64 {
    if before > 0 {
        (before as f64 - after as f64) / before as f64 * 100.0
    } else {
        0.0
    }
}

/// An arrow matching a reduction's direction, and its magnitude. A circuit
/// that grew — `--decompose-rz` expanding a rotation into Clifford+T, say —
/// used to render as "↓-5300.0%".
fn arrow(reduction: f64) -> char {
    if reduction < 0.0 { '↑' } else { '↓' }
}

fn magnitude(reduction: f64) -> String {
    format!("{:.1}", reduction.abs())
}

fn format_result_trailing(
    reduction: f64,
    before: usize,
    after: usize,
    pct_width: usize,
    count_width: usize,
) -> String {
    let arrow = arrow(reduction);
    let reduction_str = magnitude(reduction);
    let before_str = fmt_num(before);
    let after_str = fmt_num(after);
    format!(
        "{arrow}{reduction_str:>pct_width$}% · \
         {before_str:>count_width$} → {after_str:>count_width$}"
    )
}

/// Print the closing result banner. Assumes whatever ran just before it
/// (a progress box's erasure, or an "info" line like "Converged after N rounds")
/// already left exactly one blank line behind — this prints no leading
/// blank of its own.
impl Ui {
    /// Print the closing result banner (see the module docs for the box's
    /// shape). Suppressed by `--quiet`, which leaves a run reporting nothing but
    /// what was explicitly asked for.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn print_result(
        &self,
        in_gates: usize,
        out_gates: usize,
        in_2q: usize,
        out_2q: usize,
        in_depth: usize,
        out_depth: usize,
        in_t: usize,
        out_t: usize,
        in_rz: usize,
        out_rz: usize,
        secs: f64,
    ) {
        // Same box/bar rendering as the live progress boxes, but each row's
        // trailing text keeps both endpoints ("before → after") rather than just
        // the current count, since this box is printed once and never redrawn.
        let mut metrics = vec![
            ("Gates", in_gates, out_gates, GATES_BAR_COLOR),
            ("2q gates", in_2q, out_2q, TWO_QUBIT_BAR_COLOR),
            ("T/Tdg", in_t, out_t, T_BAR_COLOR),
        ];
        if in_rz > 0 || out_rz > 0 {
            metrics.push(("Rz", in_rz, out_rz, RZ_BAR_COLOR));
        }
        metrics.push(("Depth", in_depth, out_depth, DEPTH_BAR_COLOR));

        // One shared width across every row (not just each row's own before/after
        // pair), so the "→" arrows line up regardless of how much smaller a
        // metric like Rz's counts are than Gates' or Depth's.
        let width = metrics
            .iter()
            .flat_map(|&(_, before, after, _)| [before, after])
            .map(|n| fmt_num(n).chars().count())
            .max()
            .unwrap_or(0);
        let pct_width = metrics
            .iter()
            .map(|&(_, before, after, _)| magnitude(pct(before, after)).chars().count())
            .max()
            .unwrap_or(0);

        let rows: Vec<_> = metrics
            .into_iter()
            .map(|(label, before, after, color)| {
                let reduction = pct(before, after);
                (
                    label,
                    self.render_bar(reduction / 100.0, BAR_WIDTH, color),
                    format_result_trailing(reduction, before, after, pct_width, width),
                )
            })
            .collect();

        // Lead with the headline number, so the one figure most readers want is in
        // the title rather than only in the Gates row. Spelled out in words rather
        // than as the rows' ↓ arrow: "fewer"/"more" carries the direction, and a
        // circuit can grow (`--decompose-rz` expanding rotations into Clifford+T).
        let gates_reduction = pct(in_gates, out_gates);
        let direction = if gates_reduction < 0.0 {
            "more"
        } else {
            "fewer"
        };
        let title = format!(
            "Final result · {}% {direction} gates · {secs:.3}s",
            magnitude(gates_reduction)
        );
        for line in progress_box(&title, &rows) {
            self.info(&line);
        }
    }
}

/// Width, in characters, of a progress bar's fill/track region.
// 22 is 31.25% shorter than the previous 32-character bars (the nearest
// whole-character width to a 30% reduction).
const BAR_WIDTH: usize = 22;
/// Width of a progress box row's label field, with one column of padding after
/// the longest label ("2q gates").
const LABEL_WIDTH: usize = 9;

/// Green fill, used for the map-reduce chunk-completion bar.
const CHUNK_BAR_COLOR: &str = "\x1b[32m";
/// Cyan fill, used for the gate-count reduction bar.
const GATES_BAR_COLOR: &str = "\x1b[36m";
/// Yellow fill, used for the two-qubit-gate reduction bar.
const TWO_QUBIT_BAR_COLOR: &str = "\x1b[33m";
/// Magenta fill, used for the T-count reduction bar.
const T_BAR_COLOR: &str = "\x1b[35m";
/// Red fill, used for the Rz-count reduction bar.
const RZ_BAR_COLOR: &str = "\x1b[31m";
/// Blue fill, used for the depth reduction bar.
const DEPTH_BAR_COLOR: &str = "\x1b[34m";

impl Ui {
    /// Render a thin bar — heavy `color` fill, a partial tip glyph at the exact
    /// boundary, dim light-line track for the remainder — in the style of
    /// indicatif's `{bar:.color/dim}` with `━╸─` progress chars. `fraction` is
    /// clamped to `[0, 1]`.
    ///
    /// With color off the glyphs alone carry the reading (`━━━╸────`), so the
    /// bar still shows a fill level in a log file; only the escapes drop out.
    /// That also keeps its *visible* width identical either way, which
    /// [`progress_box`] relies on to stay rectangular.
    fn render_bar(&self, fraction: f64, width: usize, color: &'static str) -> String {
        let fraction = fraction.clamp(0.0, 1.0);
        let exact = fraction * width as f64;
        let full = (exact.floor() as usize).min(width);
        let has_tip = full < width && exact > full as f64;
        let empty = width - full - usize::from(has_tip);

        let mut bar = String::with_capacity(width + 16);
        bar.push_str(self.sgr(color));
        for _ in 0..full {
            bar.push('━');
        }
        if has_tip {
            bar.push('╸');
        }
        bar.push_str(self.reset());
        bar.push_str(self.sgr("\x1b[2m"));
        for _ in 0..empty {
            bar.push('─');
        }
        bar.push_str(self.reset());
        bar
    }
}

/// Number of lines a progress box with `num_rows` bar rows occupies (a top
/// and bottom border, plus one line per row).
pub(crate) fn box_lines(num_rows: usize) -> usize {
    num_rows + 2
}

/// Build the lines of a live progress box: a top border with `title`
/// embedded, one line per `(label, colored_bar, trailing)` row, and a bottom
/// border. Every line is padded to equal *visible* width — the ANSI escapes
/// inside `bar` don't count — so the box grows to fit large counts and stays
/// rectangular as values change. Indented two spaces to
/// line up with the rest of tzap's output (e.g. "  Parsing ...").
fn progress_box(title: &str, rows: &[(&str, String, String)]) -> Vec<String> {
    let row_width = |trailing: &str| LABEL_WIDTH + BAR_WIDTH + 1 + trailing.chars().count();
    let title_segment = format!("─ {title} ");
    let content_width = rows
        .iter()
        .map(|(_, _, trailing)| row_width(trailing))
        .max()
        .unwrap_or(0)
        .max(title_segment.chars().count());
    let inner_width = content_width + 2;

    let dashes = inner_width.saturating_sub(title_segment.chars().count());
    let mut lines = vec![format!("  ┌{title_segment}{}┐", "─".repeat(dashes))];
    for (label, bar, trailing) in rows {
        let pad = inner_width - (row_width(trailing) + 2);
        lines.push(format!(
            "  │ {label:<LABEL_WIDTH$}{bar} {trailing}{} │",
            " ".repeat(pad)
        ));
    }
    lines.push(format!("  └{}┘", "─".repeat(inner_width)));
    lines
}

impl Ui {
    /// Reserve `n` blank lines for a live-redrawn progress block and leave the
    /// cursor at its top-left. Pair with a later [`Ui::end_progress_block`] once
    /// the block's final frame has been drawn. A no-op without a live terminal,
    /// where nothing will be redrawn into the reserved space.
    pub(crate) fn start_progress_block(&self, n: usize) {
        if !self.live() {
            return;
        }
        eprint!("{}\x1b[{n}A", "\n".repeat(n));
        let _ = io::stderr().flush();
    }

    /// Erase a live progress block of `n` lines entirely — every line cleared,
    /// cursor returned to the block's top-left — instead of leaving its last
    /// frame on screen. Called once optimization finishes, so the box
    /// disappears rather than lingering under the closing result banner.
    pub(crate) fn end_progress_block(&self, n: usize) {
        if !self.live() {
            return;
        }
        for i in 0..n {
            eprint!("\x1b[2K");
            if i + 1 < n {
                eprintln!();
            } else if n > 1 {
                eprint!("\x1b[{}A", n - 1);
            }
        }
        let _ = io::stderr().flush();
    }

    /// Show one frame of a progress box: clear and reprint each line, then
    /// return the cursor to the block's top-left for the next frame. Must be
    /// bracketed by [`Ui::start_progress_block`] / [`Ui::end_progress_block`]
    /// with a matching line count.
    ///
    /// A no-op without a live terminal. A frame is only meaningful as
    /// something the *next* frame overwrites; appended one after another to a
    /// pipe or a log file, the same frames are just noise between the parse
    /// line and the result.
    fn draw_progress(&self, title: &str, rows: &[(&str, String, String)]) {
        if !self.live() {
            return;
        }
        self.redraw_progress_block(&progress_box(title, rows));
    }

    fn redraw_progress_block(&self, lines: &[String]) {
        let mut out = String::new();
        for (i, line) in lines.iter().enumerate() {
            out.push_str("\r\x1b[2K");
            out.push_str(line);
            if i + 1 < lines.len() {
                out.push('\n');
            }
        }
        if lines.len() > 1 {
            out.push_str(&format!("\x1b[{}A\r", lines.len() - 1));
        } else {
            out.push('\r');
        }
        eprint!("{out}");
        let _ = io::stderr().flush();
    }
    /// Redraw a live "% reduction so far" progress box under `title`: a
    /// gate, two-qubit, depth, and T-count reduction bars (reduction relative to
    /// the corresponding baselines at the start of this run), each in its own
    /// color. Shared by the fixpoint driver (title carries the
    /// iteration number) and the plain pipeline driver (no iteration — it
    /// doesn't loop). Must be bracketed by `start_progress_block(box_lines(4))`
    /// / `end_progress_block(box_lines(4))`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn update_reduction_progress(
        &self,
        title: &str,
        gates: usize,
        two_qubit: usize,
        circuit_depth: usize,
        t_count: usize,
        baseline_gates: usize,
        baseline_two_qubit: usize,
        baseline_depth: usize,
        baseline_t: usize,
        rz_count: usize,
        baseline_rz: usize,
    ) {
        let gates_pct = pct(baseline_gates, gates);
        let two_qubit_pct = pct(baseline_two_qubit, two_qubit);
        let depth_pct = pct(baseline_depth, circuit_depth);
        let t_pct = pct(baseline_t, t_count);
        let gates_str = fmt_num(gates);
        let gates_width = fmt_num(baseline_gates).chars().count();
        let two_qubit_str = fmt_num(two_qubit);
        let two_qubit_width = fmt_num(baseline_two_qubit).chars().count();
        let depth_str = fmt_num(circuit_depth);
        let depth_width = fmt_num(baseline_depth).chars().count();
        let t_str = fmt_num(t_count);
        let t_width = fmt_num(baseline_t).chars().count();
        let mut rows = vec![
            (
                "Gates",
                self.render_bar(gates_pct / 100.0, BAR_WIDTH, GATES_BAR_COLOR),
                format!("{gates_pct:>5.1}% · {gates_str:<gates_width$}"),
            ),
            (
                "2q gates",
                self.render_bar(two_qubit_pct / 100.0, BAR_WIDTH, TWO_QUBIT_BAR_COLOR),
                format!("{two_qubit_pct:>5.1}% · {two_qubit_str:<two_qubit_width$}"),
            ),
            (
                "T/Tdg",
                self.render_bar(t_pct / 100.0, BAR_WIDTH, T_BAR_COLOR),
                format!("{t_pct:>5.1}% · {t_str:<t_width$}"),
            ),
            (
                "Depth",
                self.render_bar(depth_pct / 100.0, BAR_WIDTH, DEPTH_BAR_COLOR),
                format!("{depth_pct:>5.1}% · {depth_str:<depth_width$}"),
            ),
        ];
        if baseline_rz > 0 {
            let rz_pct = pct(baseline_rz, rz_count);
            let rz_str = fmt_num(rz_count);
            let rz_width = fmt_num(baseline_rz).chars().count();
            rows.insert(
                3,
                (
                    "Rz",
                    self.render_bar(rz_pct / 100.0, BAR_WIDTH, RZ_BAR_COLOR),
                    format!("{rz_pct:>5.1}% · {rz_str:<rz_width$}"),
                ),
            );
        }
        self.draw_progress(title, &rows);
    }

    /// Redraw the live fixpoint progress box — [`update_reduction_progress`]
    /// with the current iteration number in the title.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn update_fixpoint_progress(
        &self,
        iteration: usize,
        gates: usize,
        two_qubit: usize,
        circuit_depth: usize,
        t_count: usize,
        baseline_gates: usize,
        baseline_two_qubit: usize,
        baseline_depth: usize,
        baseline_t: usize,
        rz_count: usize,
        baseline_rz: usize,
    ) {
        self.update_reduction_progress(
            &format!("Iteration {iteration} — % reduction so far"),
            gates,
            two_qubit,
            circuit_depth,
            t_count,
            baseline_gates,
            baseline_two_qubit,
            baseline_depth,
            baseline_t,
            rz_count,
            baseline_rz,
        );
    }

    /// Redraw the live parallel map-reduce progress box: how many chunks have
    /// finished, and the whole circuit's gate/T reduction achieved so far.
    /// Finished chunks contribute their optimized metrics while chunks still
    /// pending contribute their original metrics.
    ///
    /// No depth row, unlike the sequential box: a parallel run only ever measures
    /// chunks on their own, and chunk depths don't add up to the circuit's (see
    /// `Metrics::adjusted`) — the final result banner still reports it. Must be
    /// bracketed by `start_progress_block(box_lines(4))` /
    /// `end_progress_block(box_lines(4))`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn update_chunk_progress(
        &self,
        done: usize,
        total: usize,
        baseline_gates: usize,
        current_gates: usize,
        baseline_2q: usize,
        current_2q: usize,
        baseline_t: usize,
        current_t: usize,
        baseline_rz: usize,
        current_rz: usize,
    ) {
        let chunk_fraction = if total > 0 {
            done as f64 / total as f64
        } else {
            1.0
        };
        let chunk_pct = chunk_fraction * 100.0;
        let gates_pct = pct(baseline_gates, current_gates);
        let two_qubit_pct = pct(baseline_2q, current_2q);
        let t_pct = pct(baseline_t, current_t);
        let done_str = fmt_num(done);
        let done_width = fmt_num(total).chars().count();
        let total_str = fmt_num(total);
        let gates_str = fmt_num(current_gates);
        let gates_width = fmt_num(baseline_gates).chars().count();
        let two_qubit_str = fmt_num(current_2q);
        let two_qubit_width = fmt_num(baseline_2q).chars().count();
        let t_str = fmt_num(current_t);
        let t_width = fmt_num(baseline_t).chars().count();
        let mut rows = vec![
            (
                "Chunks",
                self.render_bar(chunk_fraction, BAR_WIDTH, CHUNK_BAR_COLOR),
                format!("{chunk_pct:>5.1}% · {done_str:<done_width$}/{total_str}"),
            ),
            (
                "Gates",
                self.render_bar(gates_pct / 100.0, BAR_WIDTH, GATES_BAR_COLOR),
                format!("{gates_pct:>5.1}% · {gates_str:<gates_width$}"),
            ),
            (
                "2q gates",
                self.render_bar(two_qubit_pct / 100.0, BAR_WIDTH, TWO_QUBIT_BAR_COLOR),
                format!("{two_qubit_pct:>5.1}% · {two_qubit_str:<two_qubit_width$}"),
            ),
            (
                "T/Tdg",
                self.render_bar(t_pct / 100.0, BAR_WIDTH, T_BAR_COLOR),
                format!("{t_pct:>5.1}% · {t_str:<t_width$}"),
            ),
        ];
        if baseline_rz > 0 {
            let rz_pct = pct(baseline_rz, current_rz);
            let rz_str = fmt_num(current_rz);
            let rz_width = fmt_num(baseline_rz).chars().count();
            rows.push((
                "Rz",
                self.render_bar(rz_pct / 100.0, BAR_WIDTH, RZ_BAR_COLOR),
                format!("{rz_pct:>5.1}% · {rz_str:<rz_width$}"),
            ));
        }
        self.draw_progress("Parallel optimization — % reduction so far", &rows);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A circuit that grew flips the arrow rather than negating the number —
    /// "↓-5300.0%" was previously rendered for a 3 → 162 gate expansion.
    #[test]
    fn growth_flips_the_arrow_instead_of_negating_the_percentage() {
        let trailing = format_result_trailing(pct(3, 162), 3, 162, 6, 3);
        assert!(trailing.starts_with("↑5300.0%"), "got: {trailing}");
        assert!(!trailing.contains('-'), "got: {trailing}");

        let shrunk = format_result_trailing(pct(162, 3), 162, 3, 4, 3);
        assert!(shrunk.starts_with("↓98.1%"), "got: {shrunk}");
    }

    #[test]
    fn fmt_size_picks_a_unit_that_keeps_the_number_visible() {
        assert_eq!(fmt_size(0), "0 B");
        assert_eq!(fmt_size(871), "871 B");
        assert_eq!(fmt_size(1024), "1.0 KB");
        assert_eq!(fmt_size(5_817_556), "5.5 MB");
    }

    #[test]
    fn progress_box_grows_for_large_counts() {
        let rows = [(
            "2q gates",
            Ui::plain().render_bar(1.0, BAR_WIDTH, TWO_QUBIT_BAR_COLOR),
            "↓0.0% · 1,234,567,890,123".to_string(),
        )];
        let lines = progress_box("Large counts", &rows);
        let width = lines[0].chars().count();

        // A plain `Ui` emits no escapes, so the row's rendered width *is* its
        // visible width and can be compared to the borders' directly.
        assert_eq!(lines[1].chars().count(), width);
        assert_eq!(lines[2].chars().count(), width);
        assert!(width > 60);
        assert!(lines[1].contains("1,234,567,890,123"));
    }

    #[test]
    fn progress_bars_use_the_compact_width() {
        let bar = Ui::plain().render_bar(0.5, BAR_WIDTH, GATES_BAR_COLOR);
        assert_eq!(BAR_WIDTH, 22);
        assert_eq!(bar.chars().count(), BAR_WIDTH);
        assert_eq!(bar, "━━━━━━━━━━━───────────");
    }

    #[test]
    fn final_result_aligns_single_digit_percentages() {
        let reductions = [40.5_f64, 0.0, 8.2];
        let pct_width = reductions
            .iter()
            .map(|reduction| format!("{reduction:.1}").chars().count())
            .max()
            .unwrap();
        let trailing: Vec<_> = reductions
            .iter()
            .map(|&reduction| format_result_trailing(reduction, 1, 1, pct_width, 1))
            .collect();

        let separator_columns: Vec<_> = trailing
            .iter()
            .map(|text| text.find(" · ").unwrap())
            .collect();
        assert!(separator_columns.windows(2).all(|pair| pair[0] == pair[1]));
        assert_eq!(trailing[1], "↓ 0.0% · 1 → 1");
        assert_eq!(trailing[2], "↓ 8.2% · 1 → 1");
    }
}
