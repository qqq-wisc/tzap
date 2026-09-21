//! The `--json` report: one machine-readable object on stdout describing what
//! a run did.
//!
//! Hand-rolled rather than serde-derived, deliberately. The schema is small,
//! fixed, and written in exactly one place, so the whole of serde is a poor
//! trade against this crate's build time — and every key here is a literal
//! ASCII identifier, leaving only a handful of string *values* (paths, pass
//! names) that need escaping at all.
//!
//! The object is a stable contract: keys are added over time, but an existing
//! key keeps its name, nesting, and meaning, so a script reading
//! `.metrics.output.t` keeps working across releases. Every duration is
//! `seconds` as a JSON number; every count is an integer; a value that
//! doesn't apply to this run is `null` rather than absent.

use tzap::circuit::GateSet;
use tzap::optimize::{Level, Metrics, Options, PassName, Report};

/// A value that renders itself as JSON. Only the shapes this report needs.
enum Value {
    Null,
    Bool(bool),
    Int(usize),
    /// Rendered with `{:?}`, which is the shortest representation that
    /// round-trips back to the same `f64` — and is valid JSON for every
    /// finite value. Non-finite floats can't occur here (durations and
    /// percentages are both computed from finite counts), but are mapped to
    /// `null` rather than emitted as the JSON-invalid `NaN`/`inf`.
    Float(f64),
    Str(String),
    Array(Vec<Value>),
    Object(Vec<(&'static str, Value)>),
}

impl Value {
    fn str(s: impl Into<String>) -> Value {
        Value::Str(s.into())
    }

    /// `Some` mapped through `f`, or `Value::Null` — the report's
    /// "doesn't apply to this run" case.
    fn some<T>(value: Option<T>, f: impl FnOnce(T) -> Value) -> Value {
        value.map_or(Value::Null, f)
    }

    fn write(&self, out: &mut String, indent: usize) {
        let pad = |out: &mut String, depth: usize| {
            out.push('\n');
            for _ in 0..depth {
                out.push_str("  ");
            }
        };
        match self {
            Value::Null => out.push_str("null"),
            Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Value::Int(n) => out.push_str(&n.to_string()),
            Value::Float(x) if x.is_finite() => out.push_str(&format!("{x:?}")),
            Value::Float(_) => out.push_str("null"),
            Value::Str(s) => write_string(out, s),
            Value::Array(items) if items.is_empty() => out.push_str("[]"),
            Value::Array(items) => {
                out.push('[');
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    pad(out, indent + 1);
                    item.write(out, indent + 1);
                }
                pad(out, indent);
                out.push(']');
            }
            Value::Object(fields) if fields.is_empty() => out.push_str("{}"),
            Value::Object(fields) => {
                out.push('{');
                for (i, (key, value)) in fields.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    pad(out, indent + 1);
                    write_string(out, key);
                    out.push_str(": ");
                    value.write(out, indent + 1);
                }
                pad(out, indent);
                out.push('}');
            }
        }
    }

    /// This value as a newline-terminated JSON document.
    fn document(&self) -> String {
        let mut out = String::new();
        self.write(&mut out, 0);
        out.push('\n');
        out
    }
}

/// Write `s` as a JSON string literal. Escapes the two mandatory characters,
/// the short forms for the common control characters, and `\u00xx` for the
/// rest of the C0 range — a file path really can contain a newline or a
/// quote, and unescaped either one would produce output no parser accepts.
fn write_string(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

/// One completed whole-circuit pass, as `--json` reports it.
pub(crate) struct PassRecord {
    pub(crate) stage: Option<usize>,
    pub(crate) name: String,
    pub(crate) input_gates: usize,
    pub(crate) output_gates: usize,
    pub(crate) seconds: f64,
}

/// One SuperOpt MURM load or build.
pub(crate) struct MurmRecord {
    pub(crate) stage: Option<usize>,
    pub(crate) cached: bool,
    pub(crate) basis: GateSet,
    pub(crate) seconds: f64,
}

/// How a fixpoint pipeline terminated.
pub(crate) struct FixpointRecord {
    pub(crate) stage: Option<usize>,
    pub(crate) rounds: usize,
    pub(crate) converged: bool,
}

/// Everything a run recorded that isn't already in its [`Report`]. Filled in
/// by the CLI's observer as events arrive.
#[derive(Default)]
pub(crate) struct Recording {
    pub(crate) stages: Vec<String>,
    pub(crate) passes: Vec<PassRecord>,
    pub(crate) murms: Vec<MurmRecord>,
    pub(crate) fixpoints: Vec<FixpointRecord>,
}

impl Recording {
    pub(crate) fn current_stage(&self) -> Option<usize> {
        self.stages.len().checked_sub(1)
    }
}

/// What the CLI knows about a run that the optimizer doesn't: where the
/// circuit came from and where it went.
pub(crate) struct RunInfo<'a> {
    /// The input as named on the command line, or `None` for stdin.
    pub(crate) input_path: Option<&'a str>,
    pub(crate) input_bytes: Option<u64>,
    pub(crate) input_qubits: usize,
    pub(crate) input_gate_set: GateSet,
    pub(crate) parse_seconds: f64,
    /// The output destination as named on the command line: a path, `"-"`
    /// for stdout, or `None` when the circuit was discarded.
    pub(crate) output_path: Option<&'a str>,
    pub(crate) output_gate_set: GateSet,
    pub(crate) seconds: f64,
}

fn metrics_value(m: Metrics) -> Value {
    Value::Object(vec![
        ("gates", Value::Int(m.gates)),
        ("two_qubit", Value::Int(m.two_qubit)),
        ("t", Value::Int(m.t)),
        ("rz", Value::Int(m.rz)),
        ("depth", Value::Int(m.depth)),
    ])
}

fn gate_set_value(gates: GateSet) -> Value {
    Value::Array(gates.names().map(Value::str).collect())
}

/// Percentage reduction from `before` to `after`, matching the human
/// banner's arithmetic (negative when a circuit grew, 0 when there was
/// nothing to reduce).
fn reduction(before: usize, after: usize) -> Value {
    if before == 0 {
        return Value::Float(0.0);
    }
    Value::Float((before as f64 - after as f64) / before as f64 * 100.0)
}

fn options_value(options: &Options) -> Value {
    let level = match options.level {
        Level::O1 => "O1",
        Level::O2 => "O2",
        Level::O3 => "O3",
        Level::Osuper => "Osuper",
    };
    let (qubits, window_gates, murm_entries) = options.superopt.resolved(options.level);
    Value::Object(vec![
        ("level", Value::str(level)),
        (
            "passes",
            Value::some(options.passes.as_ref(), |passes| {
                Value::Array(
                    passes
                        .iter()
                        .map(|pass| Value::str(pass_name(*pass)))
                        .collect(),
                )
            }),
        ),
        ("fixpoint", Value::Bool(options.fixpoint)),
        ("decompose_rz", Value::Bool(options.decompose_rz)),
        ("decompose_cz", Value::Bool(options.decompose_cz)),
        ("decompose_ccx", Value::Bool(options.decompose_ccx)),
        ("rz_epsilon", Value::Float(options.rz_epsilon)),
        ("parallel", Value::Bool(options.parallel)),
        (
            "superopt",
            // The *effective* bounds, presets resolved — the numbers this run
            // actually used, not the sparse set of overrides it was given.
            Value::Object(vec![
                ("qubits", Value::Int(qubits)),
                ("window_gates", Value::Int(window_gates)),
                ("murm_entries", Value::Int(murm_entries)),
                ("gates", Value::str(options.superopt_gates.as_argument())),
            ]),
        ),
    ])
}

/// A pass's canonical `--passes` spelling, so the names in a JSON report are
/// the same strings that can be fed back in on the command line.
fn pass_name(pass: PassName) -> &'static str {
    PassName::ALL
        .iter()
        .find(|(_, candidate, _)| *candidate == pass)
        .map(|(name, _, _)| *name)
        .unwrap_or("unknown")
}

/// Render a completed run as a JSON object, newline-terminated.
pub(crate) fn render(
    run: &RunInfo<'_>,
    options: &Options,
    report: &Report,
    recording: &Recording,
) -> String {
    Value::Object(vec![
        ("tzap", Value::str(env!("CARGO_PKG_VERSION"))),
        (
            "input",
            Value::Object(vec![
                ("path", Value::some(run.input_path, Value::str)),
                ("stdin", Value::Bool(run.input_path.is_none())),
                (
                    "bytes",
                    Value::some(run.input_bytes, |bytes| Value::Int(bytes as usize)),
                ),
                ("qubits", Value::Int(run.input_qubits)),
                ("gate_set", gate_set_value(run.input_gate_set)),
                ("parse_seconds", Value::Float(run.parse_seconds)),
            ]),
        ),
        (
            "output",
            Value::Object(vec![
                ("path", Value::some(run.output_path, Value::str)),
                ("stdout", Value::Bool(run.output_path == Some("-"))),
                ("gate_set", gate_set_value(run.output_gate_set)),
            ]),
        ),
        ("options", options_value(options)),
        (
            "metrics",
            Value::Object(vec![
                ("input", metrics_value(report.input)),
                ("baseline", metrics_value(report.baseline)),
                ("output", metrics_value(report.output)),
            ]),
        ),
        (
            // Against `baseline`, not `input` — currently the two are equal,
            // preserving the original circuit as the stable comparison point
            // across the native and post-decomposition stages.
            "reduction_percent",
            Value::Object(vec![
                (
                    "gates",
                    reduction(report.baseline.gates, report.output.gates),
                ),
                (
                    "two_qubit",
                    reduction(report.baseline.two_qubit, report.output.two_qubit),
                ),
                ("t", reduction(report.baseline.t, report.output.t)),
                ("rz", reduction(report.baseline.rz, report.output.rz)),
                (
                    "depth",
                    reduction(report.baseline.depth, report.output.depth),
                ),
            ]),
        ),
        (
            "stages",
            Value::Array(recording.stages.iter().cloned().map(Value::str).collect()),
        ),
        (
            "passes",
            Value::Array(
                recording
                    .passes
                    .iter()
                    .map(|pass| {
                        Value::Object(vec![
                            ("stage", Value::some(pass.stage, Value::Int)),
                            ("name", Value::str(pass.name.clone())),
                            ("input_gates", Value::Int(pass.input_gates)),
                            ("output_gates", Value::Int(pass.output_gates)),
                            ("seconds", Value::Float(pass.seconds)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            // Deprecated compatibility view: before stage-aware reporting,
            // only the last SuperOpt load was retained under `table`.
            "table",
            Value::some(recording.murms.last(), |murm| {
                Value::Object(vec![
                    ("cached", Value::Bool(murm.cached)),
                    ("seconds", Value::Float(murm.seconds)),
                ])
            }),
        ),
        (
            "murms",
            Value::Array(
                recording
                    .murms
                    .iter()
                    .map(|murm| {
                        Value::Object(vec![
                            ("stage", Value::some(murm.stage, Value::Int)),
                            ("cached", Value::Bool(murm.cached)),
                            ("basis", gate_set_value(murm.basis)),
                            ("seconds", Value::Float(murm.seconds)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            // Deprecated compatibility view: preserve the old last-fixpoint
            // object while `fixpoints` carries every stage-aware record.
            "fixpoint",
            Value::some(recording.fixpoints.last(), |fixpoint| {
                Value::Object(vec![
                    ("rounds", Value::Int(fixpoint.rounds)),
                    ("converged", Value::Bool(fixpoint.converged)),
                ])
            }),
        ),
        (
            "fixpoints",
            Value::Array(
                recording
                    .fixpoints
                    .iter()
                    .map(|fixpoint| {
                        Value::Object(vec![
                            ("stage", Value::some(fixpoint.stage, Value::Int)),
                            ("rounds", Value::Int(fixpoint.rounds)),
                            ("converged", Value::Bool(fixpoint.converged)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "cache_dir",
            Value::some(tzap::super_opt::cache_dir(), |dir| {
                Value::str(dir.display().to_string())
            }),
        ),
        ("seconds", Value::Float(run.seconds)),
    ])
    .document()
}

/// Render the `--cache-info` listing as JSON.
pub(crate) fn render_cache_info(entries: &[tzap::super_opt::CacheEntry]) -> String {
    let entries_value = || {
        Value::Array(
            entries
                .iter()
                .map(|entry| {
                    Value::Object(vec![
                        ("path", Value::str(entry.path.display().to_string())),
                        ("bytes", Value::Int(entry.bytes as usize)),
                    ])
                })
                .collect(),
        )
    };
    Value::Object(vec![
        ("tzap", Value::str(env!("CARGO_PKG_VERSION"))),
        (
            "cache_dir",
            Value::some(tzap::super_opt::cache_dir(), |dir| {
                Value::str(dir.display().to_string())
            }),
        ),
        // `tables` is the pre-MURM name and remains as a compatibility alias.
        ("tables", entries_value()),
        ("murms", entries_value()),
        (
            "total_bytes",
            Value::Int(entries.iter().map(|entry| entry.bytes as usize).sum()),
        ),
    ])
    .document()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strings_escape_everything_a_path_can_contain() {
        let mut out = String::new();
        write_string(&mut out, "a\"b\\c\nd\te\u{1}f");
        assert_eq!(out, "\"a\\\"b\\\\c\\nd\\te\\u0001f\"");
    }

    #[test]
    fn empty_containers_render_compactly() {
        let mut out = String::new();
        Value::Object(vec![
            ("a", Value::Array(vec![])),
            ("b", Value::Object(vec![])),
        ])
        .write(&mut out, 0);
        assert_eq!(out, "{\n  \"a\": [],\n  \"b\": {}\n}");
    }

    #[test]
    fn non_finite_floats_render_as_null_rather_than_invalid_json() {
        let mut out = String::new();
        Value::Float(f64::NAN).write(&mut out, 0);
        assert_eq!(out, "null");
        out.clear();
        Value::Float(f64::INFINITY).write(&mut out, 0);
        assert_eq!(out, "null");
    }

    /// A circuit that grew reports a negative reduction rather than clamping
    /// to zero, matching the human banner's `↑` row.
    #[test]
    fn reduction_is_signed_and_zero_safe() {
        let rendered = |v: Value| {
            let mut out = String::new();
            v.write(&mut out, 0);
            out
        };
        assert_eq!(rendered(reduction(10, 5)), "50.0");
        assert_eq!(rendered(reduction(5, 10)), "-100.0");
        assert_eq!(rendered(reduction(0, 0)), "0.0");
    }

    /// Every pass has a canonical name, and it is the `--passes` spelling.
    #[test]
    fn every_pass_name_round_trips_through_the_cli_spelling() {
        for (name, pass, _) in PassName::ALL {
            assert_eq!(pass_name(pass), name);
            assert_eq!(PassName::parse(pass_name(pass)), Some(pass));
        }
    }
}
