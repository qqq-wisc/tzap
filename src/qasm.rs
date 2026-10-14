//! OpenQASM 2.0 parser and serializer.
//!
//! Supports native Clifford+T, phase, axis/controlled rotations, SWAP/CSWAP,
//! CX/CZ/CCX/CCZ, U/u2/u3/id lowering, measurement, reset, and registers.
//! Custom `gate` definitions are fully inlined; unitary calls support register
//! broadcasting. Parameters use exact symbolic scalars with scoped formals.
//! Classical conditionals, opaque gates, barriers in gate bodies, and includes
//! other than `qelib1.inc` are rejected. Top-level barriers are ignored.

mod definitions;
mod expr;
#[cfg(test)]
mod gate_tests;
mod inline;
mod scan;

use crate::circuit::{CBit, Circuit, Gate, GateKind, GateSet, Qubit};
use smallvec::SmallVec;

/// Frontend resource limits for custom calls, native broadcasts, and expressions.
/// Ordinary flat circuits are not capped by `expanded_gates`.
#[derive(Clone, Copy, Debug)]
pub struct ParseLimits {
    /// Maximum nesting of custom calls (expanded with an explicit stack).
    pub call_depth: usize,
    /// Maximum expression nesting, at most 256.
    pub expression_depth: usize,
    /// Total native operations emitted by custom calls or native broadcasts.
    pub expanded_gates: usize,
    /// Work units for call visits, bound operands/parameters, and expression
    /// evaluation during implicit expansion. Cached constants cost one unit;
    /// other expressions cost one per opcode.
    pub expansion_work: usize,
}
impl Default for ParseLimits {
    fn default() -> Self {
        Self {
            call_depth: 256,
            expression_depth: 256,
            expanded_gates: 10_000_000,
            expansion_work: 20_000_000,
        }
    }
}

type Operands = SmallVec<[Qubit; 3]>;

struct Evaluation<'a> {
    limits: ParseLimits,
    angles: &'a mut AngleParser,
}

/// Parse and fully inline the supported OpenQASM 2 subset.
pub fn parse(qasm: &str) -> Result<Circuit, String> {
    parse_with_limits(qasm, ParseLimits::default())
}

/// Parse using explicit expansion limits, without changing the circuit model.
pub fn parse_with_limits(qasm: &str, limits: ParseLimits) -> Result<Circuit, String> {
    parse_with(qasm, limits, false)
}

/// Parse finite binary64 SDK literals without interpreting decimals as exact
/// mathematical rationals. Custom bodies may bind those literals structurally.
pub fn parse_numeric_sdk(qasm: &str) -> Result<Circuit, String> {
    parse_with(qasm, ParseLimits::default(), true)
}

fn parse_with(qasm: &str, limits: ParseLimits, numeric_literals: bool) -> Result<Circuit, String> {
    let mut angles = AngleParser {
        numeric_literals,
        context: crate::angle::AngleContext::with_limits(crate::angle::Limits {
            max_depth: limits
                .expression_depth
                .min(crate::angle::Limits::default().max_depth),
            ..crate::angle::Limits::default()
        }),
        ..AngleParser::default()
    };
    if limits.expression_depth > 256 {
        return Err("expression nesting limit cannot exceed 256".into());
    }
    let mut scanner = scan::Scanner::new(qasm);
    let mut registry = definitions::Registry::default();
    let mut budget = inline::Budget::default();
    let mut registers = Vec::new();
    let mut cregisters = Vec::new();
    let mut num_qubits = 0usize;
    let mut num_cbits = 0usize;
    let mut gates = Vec::new();
    let mut seen_gate = false;
    while let Some(statement) = scanner.next()? {
        let line_num = statement.line;
        let text = statement.text.trim();
        if text.is_empty() && matches!(statement.end, b';' | 0) {
            continue;
        }
        let (keyword, rest) = scan::identifier(text);
        let rest = rest.trim();
        if keyword == "gate" {
            if statement.end != b'{' {
                return Err(format!("line {line_num}: gate definition missing '{{'"));
            }
            registry.define(
                rest,
                line_num,
                &mut scanner,
                &registers,
                &cregisters,
                Evaluation {
                    limits,
                    angles: &mut angles,
                },
            )?;
            continue;
        }
        if matches!(statement.end, b'{' | b'}') {
            return Err(format!("line {line_num}: unexpected brace"));
        }
        match keyword {
            "OPENQASM" => {
                if rest != "2.0" {
                    return Err(format!(
                        "line {line_num}: unsupported OPENQASM version '{rest}'"
                    ));
                }
            }
            "include" => {
                if rest != "\"qelib1.inc\"" {
                    return Err(format!(
                        "line {line_num}: unsupported include (expected include \"qelib1.inc\")"
                    ));
                }
            }
            "barrier" => (), // Existing top-level behavior; rejected in definitions.
            "qreg" | "creg" => {
                if seen_gate {
                    return Err(format!("line {line_num}: {keyword} declaration after gate"));
                }
                let open = rest
                    .find('[')
                    .ok_or_else(|| format!("line {line_num}: invalid {keyword} declaration"))?;
                let close = rest[open + 1..]
                    .find(']')
                    .map(|i| open + 1 + i)
                    .ok_or_else(|| format!("line {line_num}: invalid {keyword} declaration"))?;
                let name = rest[..open].trim();
                if !scan::valid_name(name) || !rest[close + 1..].trim().is_empty() {
                    return Err(format!("line {line_num}: invalid {keyword} declaration"));
                }
                if registry.contains(name)
                    || (if keyword == "qreg" {
                        &registers
                    } else {
                        &cregisters
                    })
                    .iter()
                    .any(|(n, _, _)| n == name)
                {
                    return Err(format!(
                        "line {line_num}: duplicate or conflicting register name '{name}'"
                    ));
                }
                let text = rest[open + 1..close].trim();
                let size = scan::nninteger(text)
                    .ok_or_else(|| format!("line {line_num}: bad {keyword} size '{text}'"))?;
                let (lookup, count, noun) = if keyword == "qreg" {
                    (&mut registers, &mut num_qubits, "qubits")
                } else {
                    (&mut cregisters, &mut num_cbits, "classical bits")
                };
                let offset = *count;
                *count = count
                    .checked_add(size)
                    .filter(|n| *n <= u32::MAX as usize)
                    .ok_or_else(|| {
                        format!(
                            "line {line_num}: circuit declares more {noun} than the {} supported",
                            u32::MAX
                        )
                    })?;
                lookup.push((name.to_owned(), offset, size));
            }
            "measure" => {
                seen_gate = true;
                for (qubit, cbit) in parse_measure(rest, &registers, &cregisters, line_num)? {
                    gates.push(Gate::measure { qubit, cbit });
                }
            }
            "reset" => {
                seen_gate = true;
                for q in expand_qubit_operand(rest, &registers, line_num)? {
                    gates.push(Gate::reset(q));
                }
            }
            _ => {
                seen_gate = true;
                if inline::native_indexed(keyword, rest, line_num, &registers, &mut gates)? {
                    continue;
                }
                inline::call(
                    text,
                    line_num,
                    &registry,
                    &registers,
                    Evaluation {
                        limits,
                        angles: &mut angles,
                    },
                    &mut budget,
                    &mut gates,
                )?;
            }
        }
    }
    let mut circuit = Circuit::with_cbits(num_qubits, num_cbits);
    circuit.gates = gates;
    Ok(circuit)
}

/// Serialize a circuit to OpenQASM 2.0.
pub fn serialize(circuit: &Circuit) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    writeln!(s, "OPENQASM 2.0;").unwrap();
    writeln!(s, "include \"qelib1.inc\";").unwrap();
    writeln!(s, "qreg q[{}];", circuit.num_qubits).unwrap();
    if circuit.num_cbits > 0 {
        writeln!(s, "creg c[{}];", circuit.num_cbits).unwrap();
    }
    for gate in &circuit.gates {
        write_gate(&mut s, gate);
    }
    s
}

pub(crate) fn serialize_numeric(circuit: &Circuit) -> String {
    let exact = serialize(circuit);
    let mut out = String::with_capacity(exact.len());
    for line in exact.lines() {
        if let Some(start) = line.find('(') {
            let end =
                start + 1 + find_matching_paren(&line[start + 1..]).expect("serialized parameter");
            let a = crate::angle_expr::parse(&line[start + 1..end], 0).expect("validated scalar");
            let value = a
                .to_f64_exact()
                .expect("numeric export checked exact representability");
            out.push_str(&line[..start + 1]);
            out.push_str(&crate::angle::real_token(value));
            out.push_str(&line[end..]);
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }
    out
}

/// Call `f` with each statement `serialize` writes for `gate`.
///
/// Gates absent from the original qelib1.inc are expanded exactly here, so
/// metrics and gate sets computed through this function describe the emitted
/// program rather than tzap's native gate list.
pub fn for_each_emitted(gate: &Gate, mut f: impl FnMut(&Gate)) {
    match gate {
        Gate::sx(q) => {
            for term in [Gate::h(*q), Gate::s(*q), Gate::h(*q)] {
                f(&term);
            }
        }
        Gate::swap(a, b) => {
            for (control, target) in [(*a, *b), (*b, *a), (*a, *b)] {
                f(&Gate::cnot { control, target });
            }
        }
        Gate::crx {
            theta,
            control,
            target,
        }
        | Gate::cry {
            theta,
            control,
            target,
        } => {
            let cry = matches!(gate, Gate::cry { .. });
            if cry {
                f(&Gate::sdg(*target));
            }
            f(&Gate::h(*target));
            f(&Gate::crz {
                theta: theta.clone(),
                control: *control,
                target: *target,
            });
            f(&Gate::h(*target));
            if cry {
                f(&Gate::s(*target));
            }
        }
        Gate::cswap {
            control,
            first,
            second,
        } => {
            // CSWAP is absent from the original qelib1.inc. Preserve roles and
            // every branch phase with this exact expansion.
            let cx = Gate::cnot {
                control: *first,
                target: *second,
            };
            f(&cx);
            f(&Gate::ccx {
                control1: *control,
                control2: *second,
                target: *first,
            });
            f(&cx);
        }
        _ => f(gate),
    }
}

/// The statement name `serialize` writes for a kind `for_each_emitted` keeps.
pub const fn emitted_name(kind: GateKind) -> &'static str {
    match kind {
        GateKind::P => "u1",
        GateKind::Cp => "cu1",
        kind => kind.name(),
    }
}

/// Statement names in `serialize(circuit)`, in `GateKind` order.
pub fn emitted_gate_names(circuit: &Circuit) -> Vec<&'static str> {
    let mut set = GateSet::EMPTY;
    for gate in &circuit.gates {
        for_each_emitted(gate, |g| set.insert(g.kind()));
    }
    set.iter().map(emitted_name).collect()
}

fn write_gate(s: &mut String, gate: &Gate) {
    for_each_emitted(gate, |g| write_statement(s, g));
}

fn write_statement(s: &mut String, gate: &Gate) {
    use std::fmt::Write;
    match gate {
        Gate::sx(_) | Gate::swap(..) | Gate::crx { .. } | Gate::cry { .. } | Gate::cswap { .. } => {
            unreachable!("for_each_emitted expands {}", gate.kind().name())
        }
        Gate::rz(theta, _) | Gate::rx(theta, _) | Gate::ry(theta, _) | Gate::p(theta, _) => {
            write!(s, "{}({theta}) ", emitted_name(gate.kind())).unwrap();
        }
        Gate::cp { lambda, .. } => {
            write!(s, "{}({lambda}) ", emitted_name(gate.kind())).unwrap();
        }
        Gate::crz { theta, .. } => {
            write!(s, "crz({theta}) ").unwrap();
        }
        Gate::x(_)
        | Gate::h(_)
        | Gate::s(_)
        | Gate::sdg(_)
        | Gate::z(_)
        | Gate::t(_)
        | Gate::tdg(_)
        | Gate::y(_)
        | Gate::reset(_)
        | Gate::measure { .. }
        | Gate::cnot { .. }
        | Gate::cz { .. }
        | Gate::cy { .. }
        | Gate::ch { .. }
        | Gate::ccx { .. }
        | Gate::ccz { .. } => {
            s.push_str(gate.kind().name());
            s.push(' ');
        }
    }
    let (arity, qubits) = crate::circuit::qubit_operands(gate);
    for (index, &qubit) in qubits[..arity].iter().enumerate() {
        if index != 0 {
            s.push(',');
        }
        write!(s, "q[{qubit}]").unwrap();
    }
    if let Gate::measure { cbit, .. } = *gate {
        s.push_str(" -> ");
        write!(s, "c[{cbit}]").unwrap();
    }
    s.push_str(";\n");
}

#[derive(Default)]
struct AngleParser {
    numeric_literals: bool,
    context: crate::angle::AngleContext,
    cache: rustc_hash::FxHashMap<String, crate::angle::Angle>,
}
fn parse_angle_with(
    parser: &mut AngleParser,
    s: &str,
    line: usize,
) -> Result<crate::angle::Angle, String> {
    let source = s.trim();
    if parser.numeric_literals {
        return source
            .parse::<f64>()
            .map_err(|_| format!("line {line}: SDK parameter must be a numeric literal"))
            .and_then(|v| {
                crate::angle::Angle::from_f64(v).map_err(|e| format!("line {line}: {e}"))
            });
    }
    parse_math_angle_with(parser, source, line)
}

fn parse_math_angle_with(
    parser: &mut AngleParser,
    source: &str,
    line: usize,
) -> Result<crate::angle::Angle, String> {
    let source = source.trim();
    if let Some(a) = parser.cache.get(source) {
        return Ok(a.clone());
    }
    let mut result = parser.context.parse(source);
    if matches!(
        result,
        Err(crate::angle::AngleError::Limit(
            crate::symbolic_angle::Resource::Nodes
        ))
    ) {
        // Rotate the bounded weak lookup table, not the owned scalar storage.
        // Per-expression limits are still enforced on the retry.
        parser.context.clear_interner();
        result = parser.context.parse(source);
    }
    let a = result.map_err(|e| format!("line {line}: {e}"))?;
    if parser.cache.len() < 512 && source.len() <= 128 {
        parser.cache.insert(source.to_string(), a.clone());
    }
    Ok(a)
}

/// Find the position of the closing `)` that matches depth 0,
/// accounting for nested parentheses. Returns the index into `s`.
fn find_matching_paren(s: &str) -> Option<usize> {
    let mut depth = 0usize;
    for (i, c) in s.char_indices() {
        match c {
            '(' => depth += 1,
            ')' if depth == 0 => return Some(i),
            ')' => depth -= 1,
            _ => {}
        }
    }
    None
}

/// Parse a measure statement body (without the `measure ` prefix) and return one or more
/// (qubit, cbit) pairs. Supports both the indexed form `measure q[i] -> c[j];` and the
/// register-broadcast form `measure q -> c;` (per OpenQASM 2.0 §3.4, broadcast requires
/// same-size registers).
fn parse_measure(
    s: &str,
    registers: &[(String, usize, usize)],
    cregisters: &[(String, usize, usize)],
    line_num: usize,
) -> Result<Vec<(Qubit, CBit)>, String> {
    let arrow = s
        .find("->")
        .ok_or_else(|| format!("line {line_num}: measure missing '->' (got '{s}')"))?;
    let q_part = s[..arrow].trim();
    let c_part = s[arrow + 2..].trim();
    let qs = expand_qubit_operand(q_part, registers, line_num)?;
    let cs = expand_cbit_operand(c_part, cregisters, line_num)?;
    if qs.len() != cs.len() {
        return Err(format!(
            "line {line_num}: measure operand size mismatch ({} qubits, {} cbits)",
            qs.len(),
            cs.len()
        ));
    }
    Ok(qs.into_iter().zip(cs).collect())
}

/// Expand a single qubit-operand (either `name[i]` or bare `name`) into a list of qubit
/// indices. A bare register name expands to every qubit in the register, in order.
fn expand_qubit_operand(
    s: &str,
    registers: &[(String, usize, usize)],
    line_num: usize,
) -> Result<Vec<Qubit>, String> {
    let s = s.trim().trim_end_matches(';').trim();
    if s.contains('[') {
        // A bare register expands to arbitrarily many qubits, so this stays a
        // `Vec`; only the per-gate subscript path is worth keeping inline.
        return resolve_qubits(s, registers, line_num).map(Operands::into_vec);
    }
    let (_, offset, size) = registers
        .iter()
        .find(|(n, _, _)| n == s)
        .ok_or_else(|| format!("line {line_num}: unknown register '{s}'"))?;
    Ok((*offset..*offset + *size).map(|q| q as Qubit).collect())
}

fn expand_cbit_operand(
    s: &str,
    cregisters: &[(String, usize, usize)],
    line_num: usize,
) -> Result<Vec<CBit>, String> {
    let s = s.trim().trim_end_matches(';').trim();
    if s.contains('[') {
        return resolve_cbits(s, cregisters, line_num);
    }
    let (_, offset, size) = cregisters
        .iter()
        .find(|(n, _, _)| n == s)
        .ok_or_else(|| format!("line {line_num}: unknown classical register '{s}'"))?;
    Ok((*offset..*offset + *size).map(|c| c as CBit).collect())
}

fn resolve_cbits(
    s: &str,
    registers: &[(String, usize, usize)],
    line: usize,
) -> Result<Vec<CBit>, String> {
    definitions::comma_list(s)
        .map(|part| match inline::operand(part, registers, line, true)? {
            inline::Operand::Indexed(q) => Ok(q),
            _ => Err(format!("line {line}: expected indexed cbit operand")),
        })
        .collect()
}

fn resolve_qubits(
    s: &str,
    registers: &[(String, usize, usize)],
    line: usize,
) -> Result<Operands, String> {
    definitions::comma_list(s)
        .map(
            |part| match inline::operand(part, registers, line, false)? {
                inline::Operand::Indexed(q) => Ok(q),
                _ => Err(format!("line {line}: expected indexed qubit operand")),
            },
        )
        .collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn new_gates_reject_bad_angles_operands_and_arity() {
        for statement in [
            "cp(NaN) q[0],q[1]",
            "crx(inf) q[0],q[1]",
            "cry(0.3 q[0],q[1]",
            "crz(pi) q[3],q[0]",
            "cy q[0]",
            "ch q[0],q[1],q[2]",
            "cswap q[0],q[1]",
            "p(0.3) q[0],q[1]",
        ] {
            assert!(
                parse(&format!("qreg q[3]; {statement};")).is_err(),
                "{statement}"
            );
        }
    }

    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn z_from_qasm() {
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nz q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 1);
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(&c.gates[0], Gate::z(0)));
    }

    #[test]
    fn sdg_from_qasm() {
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nsdg q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 1);
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(&c.gates[0], Gate::sdg(0)));
    }

    #[test]
    fn z_qasm_roundtrip() {
        let mut c = Circuit::new(2);
        c.apply(Gate::z(0));
        c.apply(Gate::z(1));
        let qasm = serialize(&c);
        let c2 = parse(&qasm).unwrap();
        assert_eq!(c2.gates.len(), 2);
        assert!(matches!(&c2.gates[0], Gate::z(0)));
        assert!(matches!(&c2.gates[1], Gate::z(1)));
    }

    #[test]
    fn sdg_qasm_roundtrip() {
        let mut c = Circuit::new(2);
        c.apply(Gate::sdg(0));
        c.apply(Gate::sdg(1));
        let qasm = serialize(&c);
        let c2 = parse(&qasm).unwrap();
        assert_eq!(c2.gates.len(), 2);
        assert!(matches!(&c2.gates[0], Gate::sdg(0)));
        assert!(matches!(&c2.gates[1], Gate::sdg(1)));
    }

    #[test]
    fn mixed_gates_qasm_roundtrip() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(0));
        c.apply(Gate::z(0));
        c.apply(Gate::sdg(1));
        c.apply(Gate::s(2));
        c.apply(Gate::t(0));
        c.apply(Gate::tdg(1));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        c.apply(Gate::x(2));
        let qasm = serialize(&c);
        let c2 = parse(&qasm).unwrap();
        assert_eq!(c2.gates.len(), 8);
        assert!(matches!(&c2.gates[1], Gate::z(0)));
        assert!(matches!(&c2.gates[2], Gate::sdg(1)));
    }

    #[test]
    fn rz_roundtrip() {
        let mut c = Circuit::new(1);
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        let qasm = serialize(&c);
        let c2 = parse(&qasm).unwrap();
        assert_eq!(c2.gates.len(), 1);
        if let Gate::rz(theta, 0) = &c2.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz gate");
        }
    }

    #[test]
    fn z_to_qasm() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        let qasm = serialize(&c);
        assert!(qasm.contains("z q[0];"));
    }

    #[test]
    fn sdg_to_qasm() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        let qasm = serialize(&c);
        assert!(qasm.contains("sdg q[0];"));
    }

    // --- comment parsing tests ---

    #[test]
    fn semicolons_in_line_comments_never_execute_gates() {
        for body in [
            "// ignored; x q[0];\nh q[0];",
            "h q[0]; // ignored; x q[0];",
            "// ignored; x q[0];\nh q[0]; // ignored; t q[0];",
            "// unicode λ; x q[0];\nh q[0];",
        ] {
            let c = parse(&format!("OPENQASM 2.0;\nqreg q[1];\n{body}")).unwrap();
            assert_eq!(c.gates, vec![Gate::h(0)], "{body}");
        }
    }

    #[test]
    fn block_comment_markers_inside_line_comments_do_not_erase_gates() {
        for body in [
            "// /*\nh q[0];\n// */",
            "// /* unclosed\nh q[0];",
            "// */ /* ; x q[0];\nh q[0];",
        ] {
            let c = parse(&format!("OPENQASM 2.0;\nqreg q[1];\n{body}")).unwrap();
            assert_eq!(c.gates, vec![Gate::h(0)], "{body}");
        }
    }

    #[test]
    fn line_comment_markers_inside_block_comments_do_not_hide_the_end() {
        let c = parse("OPENQASM 2.0;\nqreg q[1];\n/* // ; x q[0]; */ h q[0];").unwrap();
        assert_eq!(c.gates, vec![Gate::h(0)]);
    }

    #[test]
    fn comment_markers_inside_strings_are_preserved() {
        let source = "include \"directory//file/*name*/.inc\"; // /* ignored; x q[0];\nh q[0];";
        let mut scanner = scan::Scanner::new(source);
        assert_eq!(
            scanner.next().unwrap().unwrap().text,
            "include \"directory//file/*name*/.inc\""
        );
        assert_eq!(scanner.next().unwrap().unwrap().text, "h q[0]");
        assert!(matches!(
            scan::Scanner::new("include \"directory//file.inc\";")
                .next()
                .unwrap()
                .unwrap()
                .text,
            std::borrow::Cow::Borrowed(_)
        ));
    }

    #[test]
    fn block_comments_separate_tokens_instead_of_joining_them() {
        for gate in ["r/* comment */z(1) q[0];", "rz(1/* comment */2) q[0];"] {
            let source = format!("OPENQASM 2.0;\nqreg q[1];\n{gate}");
            assert!(parse(&source).is_err(), "joined tokens in {gate}");
        }
        let c = parse("OPENQASM 2.0;\nqreg q[1];\nh/* comment */q[0];").unwrap();
        assert_eq!(c.gates, vec![Gate::h(0)]);
    }

    #[test]
    fn mixed_comments_preserve_error_line_numbers() {
        let source =
            "OPENQASM 2.0;\r\nqreg q[1];\r\n// /* ; x q[0];\r\n/* skip\nthis */\nrz(1/0) q[0];\n";
        let err = parse(source).unwrap_err();
        assert!(err.contains("line 6"), "{err}");
        assert!(err.contains("finite"), "{err}");
    }

    #[test]
    fn line_comment_only() {
        let qasm =
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n// just a comment\nqreg q[1];\nh q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
    }

    #[test]
    fn inline_line_comment() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nh q[0]; // apply hadamard\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(&c.gates[0], Gate::h(0)));
    }

    #[test]
    fn block_comment_single_line() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\n/* comment */ h q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
    }

    #[test]
    fn block_comment_multiline() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\n/* this is\na multi-line\ncomment */\nh q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
    }

    #[test]
    fn block_comment_inline() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nh /* surprise */ q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
    }

    #[test]
    fn block_comment_between_gates() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\nh q[0];\n/* between */\ncx q[0],q[1];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 2);
    }

    #[test]
    fn multiple_block_comments() {
        let qasm = "OPENQASM 2.0;\n/* a */ qreg q[1]; /* b */\n/* c */ h q[0]; /* d */\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
    }

    #[test]
    fn block_comment_spanning_gate() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\nh q[0];\n/* cx q[0],q[1]; */\nt q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 2);
        assert!(matches!(&c.gates[0], Gate::h(0)));
        assert!(matches!(&c.gates[1], Gate::t(0)));
    }

    #[test]
    fn block_and_line_comments_mixed() {
        let qasm = "\
OPENQASM 2.0;
qreg q[2];
// line comment
h q[0]; // inline
/* block */ cx q[0],q[1];
/* multi
   line */
t q[0];
";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 3);
    }

    #[test]
    fn unclosed_block_comment_is_rejected() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nh q[0];\n/* unclosed\nt q[0];\n";
        assert!(
            parse(qasm)
                .unwrap_err()
                .contains("unterminated block comment")
        );
    }

    #[test]
    fn empty_block_comment() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\n/**/ h q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
    }

    #[test]
    fn comment_only_file() {
        let qasm = "OPENQASM 2.0;\n// nothing here\n/* also nothing */\nqreg q[1];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 0);
    }

    #[test]
    fn line_comment_at_end_no_newline() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nh q[0]; // trailing";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
    }

    #[test]
    fn block_comment_preserves_line_numbers() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\n/* skip\nthis\n*/\nh q[0];\nfoo q[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(
            err.contains("line 7"),
            "expected line 7 in error, got: {err}"
        );
    }

    #[test]
    fn unsupported_gate_error() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nunknown(0.5,0.2,0.1) q[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 3"));
        assert!(err.contains("unsupported"));
        assert!(err.contains("unknown"));
    }

    // --- pi expression tests ---

    #[test]
    fn nonfinite_angles_and_intermediate_results_are_rejected() {
        for expression in ["1/0", "0/0", "-1/0", "1/(1/0)", "0*(1/0)"] {
            let source = format!("OPENQASM 2.0;\nqreg q[1];\nrz({expression}) q[0];\n");
            let err = parse(&source).unwrap_err();
            assert!(err.contains("line 3"), "{expression}: {err}");
            assert!(err.contains("finite"), "{expression}: {err}");
        }
    }

    #[test]
    fn finite_extreme_angles_and_valid_division_are_accepted() {
        for (expression, _expected) in [
            ("0", 0.0),
            ("1e308", 1e308),
            ("-1e308", -1e308),
            ("1e-300", 1e-300),
            ("1e308/1e308", 1.0),
            ("0/1", 0.0),
        ] {
            let source = format!("OPENQASM 2.0;\nqreg q[1];\nrz({expression}) q[0];\n");
            assert_eq!(
                parse(&source).unwrap().gates,
                vec![Gate::rz(
                    crate::angle::AngleContext::default()
                        .parse(expression)
                        .unwrap(),
                    0
                )]
            );
        }
    }

    #[test]
    fn multi_qubit_gates_require_distinct_operands() {
        for gate in ["cx", "cz"] {
            for a in 0..3 {
                for b in 0..3 {
                    let source = format!("OPENQASM 2.0;\nqreg q[3];\n{gate} q[{a}],q[{b}];\n");
                    let result = parse(&source);
                    if a == b {
                        let err = result.unwrap_err();
                        assert!(
                            err.contains("line 3")
                                && err.contains(gate)
                                && err.contains("distinct"),
                            "{err}"
                        );
                    } else {
                        assert!(result.is_ok(), "{source}");
                    }
                }
            }
        }
        for gate in ["ccx", "ccz"] {
            for a in 0..3 {
                for b in 0..3 {
                    for c in 0..3 {
                        let source =
                            format!("OPENQASM 2.0;\nqreg q[3];\n{gate} q[{a}],q[{b}],q[{c}];\n");
                        let result = parse(&source);
                        if a == b || a == c || b == c {
                            let err = result.unwrap_err();
                            assert!(
                                err.contains("line 3")
                                    && err.contains(gate)
                                    && err.contains("distinct"),
                                "{err}"
                            );
                        } else {
                            assert!(result.is_ok(), "{source}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn rz_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_pi_over_4() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_2_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(2*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 2.0 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_3_pi_over_4() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(3*pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 3.0 * PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_neg_pi_over_2() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(-pi/2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - (-PI / 2.0)).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_neg_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(-pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - (-PI)).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_plain_float() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(0.123456789) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 0.123456789).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_20_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(20*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 20.0 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_pi_times_2() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi*2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 2.0 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_spaces_around_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz( pi / 4 ) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_spaces_coeff() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz( 3 * pi / 4 ) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 3.0 * PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_neg_with_spaces() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(- pi / 4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - (-PI / 4.0)).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_pi_times_0_5() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(0.5*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 0.5 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_nested_parens() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz((pi/4)) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_complex_expr() {
        // 2*2*(pi/2)*(3/4*pi) = 4 * (pi/2) * (0.75*pi) = 4 * 0.5*pi * 0.75*pi = 3*pi^2/2
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(2*2*(pi/2)*(3/4*pi)) q[0];\n";
        let c = parse(qasm).unwrap();
        let expected = 2.0 * 2.0 * (PI / 2.0) * (3.0 / 4.0 * PI);
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - expected).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_addition() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi/4 + pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI / 2.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_subtraction() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi - pi/2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI / 2.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_double_neg() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(--pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_scientific_notation() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(1e-3) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - 1e-3).abs() < 1e-15);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn creg_parsed() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nh q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 1);
        assert_eq!(c.num_cbits, 1);
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(&c.gates[0], Gate::h(0)));
    }

    #[test]
    fn creg_after_gate_error() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nh q[0];\ncreg c[1];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 4"));
        assert!(err.contains("creg declaration after gate"));
    }

    // --- measurement and reset parsing ---

    #[test]
    fn measure_basic() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nmeasure q[0] -> c[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 1);
        assert_eq!(c.num_cbits, 1);
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(&c.gates[0], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(c.has_measurement());
    }

    #[test]
    fn measure_offset() {
        let qasm = "OPENQASM 2.0;\nqreg q[3];\ncreg c[3];\nmeasure q[2] -> c[1];\n";
        let c = parse(qasm).unwrap();
        assert!(matches!(&c.gates[0], Gate::measure { qubit: 2, cbit: 1 }));
    }

    #[test]
    fn measure_multi_register() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[2];\ncreg x[1];\ncreg y[2];\n\
                    measure b[1] -> y[0];\n";
        let c = parse(qasm).unwrap();
        // a[0]=0, b[0]=1, b[1]=2; x[0]=0, y[0]=1, y[1]=2
        assert!(matches!(&c.gates[0], Gate::measure { qubit: 2, cbit: 1 }));
    }

    #[test]
    fn measure_unknown_qreg() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nmeasure nope[0] -> c[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown register"));
    }

    #[test]
    fn measure_unknown_creg() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nmeasure q[0] -> nope[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown classical register"));
    }

    #[test]
    fn measure_cbit_out_of_range() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nmeasure q[0] -> c[5];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("out of range"));
    }

    #[test]
    fn measure_qubit_out_of_range() {
        // Symmetric to cbit out-of-range — the qubit side is also checked.
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[5];\nmeasure q[7] -> c[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 4"));
        assert!(err.contains("out of range"));
        assert!(err.contains("'q'"));
    }

    #[test]
    fn reset_qubit_out_of_range() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\nreset q[9];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 3"));
        assert!(err.contains("out of range"));
        assert!(err.contains("'q'"));
    }

    #[test]
    fn measure_size_check_reports_both_sizes() {
        // The size-mismatch error should include both sizes so the user can see which side is wrong.
        let qasm = "OPENQASM 2.0;\nqreg q[4];\ncreg c[2];\nmeasure q -> c;\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("size mismatch"));
        assert!(err.contains('4'));
        assert!(err.contains('2'));
    }

    #[test]
    fn programmatic_measure_skips_size_check() {
        // Construction-time validation is NOT performed: the Circuit API trusts callers.
        // (Parser-time validation is enforced; this is the internal-construction contract.)
        let mut c = Circuit::with_cbits(1, 1);
        c.apply(Gate::measure {
            qubit: 99,
            cbit: 99,
        });
        assert_eq!(c.gates.len(), 1);
        assert!(c.has_measurement());
        // But serialize → parse will fail to round-trip because the indices reference
        // bits outside the declared `qreg q[1]` / `creg c[1]`.
        let qasm = serialize(&c);
        let err = parse(&qasm).unwrap_err();
        assert!(err.contains("out of range"));
    }

    #[test]
    fn programmatic_reset_skips_size_check() {
        let mut c = Circuit::new(1);
        c.apply(Gate::reset(42));
        assert_eq!(c.gates.len(), 1);
        assert!(c.has_measurement());
        let qasm = serialize(&c);
        let err = parse(&qasm).unwrap_err();
        assert!(err.contains("out of range"));
    }

    #[test]
    fn measure_broadcast_size_mismatch_line_number() {
        // Size check is at parse time and reports the correct line number.
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[2];\nh q[0];\nh q[1];\nmeasure q -> c[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 6"));
        assert!(err.contains("size mismatch"));
    }

    #[test]
    fn measure_semicolon_separated_on_one_line() {
        // Multiple statements per line should each be parsed independently.
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[2];\nh q[0]; measure q -> c;\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 3);
        assert!(matches!(&parsed.gates[0], Gate::h(0)));
        assert!(matches!(
            &parsed.gates[1],
            Gate::measure { qubit: 0, cbit: 0 }
        ));
        assert!(matches!(
            &parsed.gates[2],
            Gate::measure { qubit: 1, cbit: 1 }
        ));
    }

    #[test]
    fn reset_semicolon_separated_on_one_line() {
        let qasm = "OPENQASM 2.0;\nqreg q[3];\nh q[0]; reset q;\n";
        let parsed = parse(qasm).unwrap();
        // H + 3 resets
        assert_eq!(parsed.gates.len(), 4);
        assert!(matches!(&parsed.gates[0], Gate::h(0)));
        for (i, g) in parsed.gates[1..].iter().enumerate() {
            assert!(matches!(g, Gate::reset(q) if *q as usize == i));
        }
    }

    #[test]
    fn measure_broadcast_with_inline_comment() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[2];\nmeasure q -> c; // sink\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 2);
    }

    #[test]
    fn reset_broadcast_with_extra_whitespace() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\nreset    q   ;\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 2);
        assert!(matches!(&parsed.gates[0], Gate::reset(0)));
        assert!(matches!(&parsed.gates[1], Gate::reset(1)));
    }

    #[test]
    fn measure_broadcast_with_extra_whitespace() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[2];\nmeasure   q   ->   c  ;\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 2);
        assert!(matches!(
            &parsed.gates[0],
            Gate::measure { qubit: 0, cbit: 0 }
        ));
        assert!(matches!(
            &parsed.gates[1],
            Gate::measure { qubit: 1, cbit: 1 }
        ));
    }

    #[test]
    fn measure_to_register_bracketed_lhs_size_one() {
        // q[0] resolves to size 1; c is a creg of size 1 → broadcast matches.
        let qasm = "OPENQASM 2.0;\nqreg q[3];\ncreg c[1];\nmeasure q[2] -> c;\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 1);
        assert!(matches!(
            &parsed.gates[0],
            Gate::measure { qubit: 2, cbit: 0 }
        ));
    }

    #[test]
    fn measure_lhs_broadcast_rhs_indexed_size_one_qreg() {
        // q is size 1 → expands to [q0]; c[1] is one cbit. Both size 1.
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[3];\nmeasure q -> c[1];\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 1);
        assert!(matches!(
            &parsed.gates[0],
            Gate::measure { qubit: 0, cbit: 1 }
        ));
    }

    #[test]
    fn measure_missing_arrow() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nmeasure q[0] c[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("missing '->'"));
    }

    #[test]
    fn reset_basic() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nreset q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 1);
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(&c.gates[0], Gate::reset(0)));
        assert!(c.has_measurement());
    }

    #[test]
    fn reset_offset_register() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[2];\nreset b[1];\n";
        let c = parse(qasm).unwrap();
        assert!(matches!(&c.gates[0], Gate::reset(2)));
    }

    #[test]
    fn reset_unknown_register() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nreset nope[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown register"));
    }

    #[test]
    fn measure_reset_roundtrip() {
        let mut c = Circuit::with_cbits(2, 2);
        c.apply(Gate::h(0));
        c.apply(Gate::reset(1));
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        c.apply(Gate::measure { qubit: 1, cbit: 1 });
        let qasm = serialize(&c);
        assert!(qasm.contains("creg c[2];"));
        assert!(qasm.contains("reset q[1];"));
        assert!(qasm.contains("measure q[0] -> c[0];"));
        assert!(qasm.contains("measure q[1] -> c[1];"));
        let c2 = parse(&qasm).unwrap();
        assert_eq!(c2.num_qubits, 2);
        assert_eq!(c2.num_cbits, 2);
        assert_eq!(c2.gates.len(), 4);
        assert!(matches!(&c2.gates[1], Gate::reset(1)));
        assert!(matches!(&c2.gates[2], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(matches!(&c2.gates[3], Gate::measure { qubit: 1, cbit: 1 }));
        assert!(c2.has_measurement());
    }

    #[test]
    fn serialize_omits_creg_when_zero() {
        let mut c = Circuit::new(1);
        c.apply(Gate::h(0));
        let qasm = serialize(&c);
        assert!(!qasm.contains("creg"));
    }

    #[test]
    fn measure_mixed_circuit_qasm() {
        let qasm = "\
OPENQASM 2.0;
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 2);
        assert_eq!(c.num_cbits, 2);
        assert_eq!(c.gates.len(), 4);
        assert!(matches!(&c.gates[0], Gate::h(0)));
        assert!(matches!(
            &c.gates[1],
            Gate::cnot {
                control: 0,
                target: 1
            }
        ));
        assert!(matches!(&c.gates[2], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(matches!(&c.gates[3], Gate::measure { qubit: 1, cbit: 1 }));
    }

    #[test]
    fn reset_with_other_gates() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nh q[0];\nreset q[0];\nh q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 3);
        assert!(matches!(&c.gates[1], Gate::reset(0)));
    }

    // --- multiple creg tests ---

    #[test]
    fn two_cregs_offsets() {
        // a[1] at q-offset 0, b[2] at q-offset 1; x[2] at c-offset 0, y[3] at c-offset 2.
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[2];\ncreg x[2];\ncreg y[3];\n\
                    measure a[0] -> x[0];\n\
                    measure a[0] -> x[1];\n\
                    measure b[0] -> y[0];\n\
                    measure b[1] -> y[2];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 3);
        assert_eq!(c.num_cbits, 5);
        assert_eq!(c.gates.len(), 4);
        // x[0]→0, x[1]→1, y[0]→2, y[2]→4
        assert!(matches!(&c.gates[0], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(matches!(&c.gates[1], Gate::measure { qubit: 0, cbit: 1 }));
        assert!(matches!(&c.gates[2], Gate::measure { qubit: 1, cbit: 2 }));
        assert!(matches!(&c.gates[3], Gate::measure { qubit: 2, cbit: 4 }));
    }

    #[test]
    fn three_cregs_full_offsets() {
        let qasm = "OPENQASM 2.0;\nqreg q[3];\n\
                    creg first[1];\ncreg second[2];\ncreg third[1];\n\
                    measure q[0] -> first[0];\n\
                    measure q[1] -> second[0];\n\
                    measure q[1] -> second[1];\n\
                    measure q[2] -> third[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_cbits, 4);
        // first[0]→0, second[0]→1, second[1]→2, third[0]→3
        assert!(matches!(&c.gates[0], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(matches!(&c.gates[1], Gate::measure { qubit: 1, cbit: 1 }));
        assert!(matches!(&c.gates[2], Gate::measure { qubit: 1, cbit: 2 }));
        assert!(matches!(&c.gates[3], Gate::measure { qubit: 2, cbit: 3 }));
    }

    #[test]
    fn multi_creg_out_of_range_uses_correct_size() {
        // x has size 1; x[1] is out of range, but if the parser were using the
        // total cbit count (3) it would wrongly accept.
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg x[1];\ncreg y[2];\n\
                    measure q[0] -> x[1];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("out of range"));
        assert!(err.contains("'x'"));
    }

    #[test]
    fn multi_creg_unknown_register_after_known() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg x[1];\ncreg y[1];\n\
                    measure q[0] -> z[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown classical register"));
        assert!(err.contains("'z'"));
    }

    #[test]
    fn multi_creg_roundtrip_collapses_to_single_creg() {
        // Serializer emits a single `creg c[N]` with N = total cbits. A parse → serialize
        // cycle on a circuit with separately-named cregs reads them as a single contiguous
        // register; the cbit indices in the gates stay the same.
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg x[1];\ncreg y[2];\n\
                    measure q[0] -> x[0];\n\
                    measure q[1] -> y[1];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_cbits, 3);
        assert!(matches!(&c.gates[0], Gate::measure { qubit: 0, cbit: 0 }));
        // y[1] → cbit offset 1 + 1 = 2
        assert!(matches!(&c.gates[1], Gate::measure { qubit: 1, cbit: 2 }));

        let qasm2 = serialize(&c);
        assert!(qasm2.contains("creg c[3];"));
        assert!(qasm2.contains("measure q[0] -> c[0];"));
        assert!(qasm2.contains("measure q[1] -> c[2];"));

        let c2 = parse(&qasm2).unwrap();
        assert_eq!(c2.num_qubits, c.num_qubits);
        assert_eq!(c2.num_cbits, c.num_cbits);
        assert_eq!(c2.gates.len(), c.gates.len());
        assert!(matches!(&c2.gates[0], Gate::measure { qubit: 0, cbit: 0 }));
        assert!(matches!(&c2.gates[1], Gate::measure { qubit: 1, cbit: 2 }));
    }

    #[test]
    fn multi_creg_same_name_qreg_distinct() {
        // A qreg and a creg may share a name (different namespaces); parser must pick
        // the right one for each side of a measure.
        let qasm = "OPENQASM 2.0;\nqreg c[2];\ncreg c[2];\n\
                    measure c[0] -> c[1];\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.num_qubits, 2);
        assert_eq!(parsed.num_cbits, 2);
        assert!(matches!(
            &parsed.gates[0],
            Gate::measure { qubit: 0, cbit: 1 }
        ));
    }

    // --- broadcast (whole-register) measure and reset (OpenQASM 2.0 §3.4) ---

    #[test]
    fn measure_broadcast_whole_register() {
        // `measure q -> c;` expands to one measure per qubit in same-size registers.
        let qasm = "OPENQASM 2.0;\nqreg q[3];\ncreg c[3];\nmeasure q -> c;\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.num_qubits, 3);
        assert_eq!(parsed.num_cbits, 3);
        assert_eq!(parsed.gates.len(), 3);
        assert!(matches!(
            &parsed.gates[0],
            Gate::measure { qubit: 0, cbit: 0 }
        ));
        assert!(matches!(
            &parsed.gates[1],
            Gate::measure { qubit: 1, cbit: 1 }
        ));
        assert!(matches!(
            &parsed.gates[2],
            Gate::measure { qubit: 2, cbit: 2 }
        ));
        assert!(parsed.has_measurement());
    }

    #[test]
    fn measure_broadcast_with_offsets() {
        // Multiple qregs and cregs; broadcast picks the named register only.
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[2];\ncreg x[1];\ncreg y[2];\n\
                    measure b -> y;\n";
        let parsed = parse(qasm).unwrap();
        // b's qubits are 1, 2; y's cbits are 1, 2.
        assert_eq!(parsed.gates.len(), 2);
        assert!(matches!(
            &parsed.gates[0],
            Gate::measure { qubit: 1, cbit: 1 }
        ));
        assert!(matches!(
            &parsed.gates[1],
            Gate::measure { qubit: 2, cbit: 2 }
        ));
    }

    #[test]
    fn measure_broadcast_size_mismatch_errors() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[3];\nmeasure q -> c;\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("size mismatch"));
    }

    #[test]
    fn measure_broadcast_unknown_qreg() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[2];\nmeasure nope -> c;\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown register"));
        assert!(err.contains("'nope'"));
    }

    #[test]
    fn measure_broadcast_unknown_creg() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\ncreg c[2];\nmeasure q -> nope;\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown classical register"));
        assert!(err.contains("'nope'"));
    }

    #[test]
    fn measure_mixed_indexed_lhs_broadcast_rhs_size_mismatch() {
        // `measure q[0] -> c;` — LHS has size 1, RHS broadcasts to 2; mismatch.
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[2];\nmeasure q[0] -> c;\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("size mismatch"));
    }

    #[test]
    fn measure_mixed_size_one_register_matches_indexed() {
        // `measure q -> c[0];` where q is size 1 — both sides have size 1, allowed.
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[2];\nmeasure q -> c[0];\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 1);
        assert!(matches!(
            &parsed.gates[0],
            Gate::measure { qubit: 0, cbit: 0 }
        ));
    }

    #[test]
    fn reset_broadcast_whole_register() {
        let qasm = "OPENQASM 2.0;\nqreg q[3];\nreset q;\n";
        let parsed = parse(qasm).unwrap();
        assert_eq!(parsed.gates.len(), 3);
        assert!(matches!(&parsed.gates[0], Gate::reset(0)));
        assert!(matches!(&parsed.gates[1], Gate::reset(1)));
        assert!(matches!(&parsed.gates[2], Gate::reset(2)));
        assert!(parsed.has_measurement());
    }

    #[test]
    fn reset_broadcast_with_offset() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[2];\nreset b;\n";
        let parsed = parse(qasm).unwrap();
        // b's qubits are 1, 2.
        assert_eq!(parsed.gates.len(), 2);
        assert!(matches!(&parsed.gates[0], Gate::reset(1)));
        assert!(matches!(&parsed.gates[1], Gate::reset(2)));
    }

    #[test]
    fn reset_broadcast_unknown_register() {
        let qasm = "OPENQASM 2.0;\nqreg q[2];\nreset nope;\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown register"));
        assert!(err.contains("'nope'"));
    }

    // --- multiple quantum register tests ---

    #[test]
    fn two_registers() {
        let qasm = "OPENQASM 2.0;\nqreg a[2];\nqreg b[3];\nh a[0];\nh a[1];\nt b[0];\nt b[2];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 5);
        assert_eq!(c.gates.len(), 4);
        // a[0] -> 0, a[1] -> 1, b[0] -> 2, b[2] -> 4
        assert!(matches!(&c.gates[0], Gate::h(0)));
        assert!(matches!(&c.gates[1], Gate::h(1)));
        assert!(matches!(&c.gates[2], Gate::t(2)));
        assert!(matches!(&c.gates[3], Gate::t(4)));
    }

    #[test]
    fn multi_register_cnot() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[1];\ncx a[0],b[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 2);
        assert!(matches!(
            &c.gates[0],
            Gate::cnot {
                control: 0,
                target: 1
            }
        ));
    }

    #[test]
    fn qreg_after_gate_error() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nh a[0];\nqreg b[1];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 4"));
        assert!(err.contains("qreg declaration after gate"));
    }

    #[test]
    fn unknown_register_error() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nh b[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown register"));
    }

    #[test]
    fn register_index_out_of_range() {
        let qasm = "OPENQASM 2.0;\nqreg a[2];\nh a[5];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("out of range"));
    }

    #[test]
    fn three_registers_offsets() {
        let qasm = "OPENQASM 2.0;\nqreg x[3];\nqreg y[2];\nqreg z[1];\n\
                     h x[0];\nh x[2];\nh y[0];\nh y[1];\nh z[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 6);
        assert_eq!(c.gates.len(), 5);
        // x[0]->0, x[2]->2, y[0]->3, y[1]->4, z[0]->5
        assert!(matches!(&c.gates[0], Gate::h(0)));
        assert!(matches!(&c.gates[1], Gate::h(2)));
        assert!(matches!(&c.gates[2], Gate::h(3)));
        assert!(matches!(&c.gates[3], Gate::h(4)));
        assert!(matches!(&c.gates[4], Gate::h(5)));
    }

    #[test]
    fn multi_register_ccx() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[1];\nqreg c[1];\nccx a[0],b[0],c[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 3);
        assert!(matches!(
            &c.gates[0],
            Gate::ccx {
                control1: 0,
                control2: 1,
                target: 2
            }
        ));
    }

    #[test]
    fn ccz_stays_native() {
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nccz q[2],q[0],q[1];\n";
        let c = parse(qasm).unwrap();

        assert_eq!(c.gates.len(), 1);
        assert!(matches!(
            c.gates[0],
            Gate::ccz {
                control1: 2,
                control2: 0,
                target: 1
            }
        ));
        assert!(!c.has_toffoli());
        assert!(c.has_ccz());

        let serialized = serialize(&c);
        assert!(serialized.contains("ccz q[2],q[0],q[1];"));
        let reparsed = parse(&serialized).unwrap();
        assert!(matches!(
            reparsed.gates.as_slice(),
            [Gate::ccz {
                control1: 2,
                control2: 0,
                target: 1
            }]
        ));
    }

    #[test]
    fn multi_register_ccz() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[1];\nqreg c[1];\nccz a[0],b[0],c[0];\n";
        let c = parse(qasm).unwrap();

        assert_eq!(c.num_qubits, 3);
        assert!(matches!(
            c.gates[0],
            Gate::ccz {
                control1: 0,
                control2: 1,
                target: 2
            }
        ));
    }

    #[test]
    fn ccz_rejects_wrong_operand_count() {
        for (operands, count) in [("q[0],q[1]", 2), ("q[0],q[1],q[2],q[3]", 4)] {
            let qasm = format!("OPENQASM 2.0;\nqreg q[4];\nccz {operands};\n");
            let err = parse(&qasm).unwrap_err();
            assert!(
                err.contains(&format!("ccz expects 3 qubit operands, got {count}")),
                "{err}"
            );
        }
    }

    #[test]
    fn serialize_programmatic_ccz() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccz {
            control1: 1,
            control2: 2,
            target: 0,
        });

        assert!(serialize(&c).ends_with("ccz q[1],q[2],q[0];\n"));
    }

    #[test]
    fn multiple_ccz_with_comments_round_trip() {
        let qasm = "OPENQASM 2.0;\nqreg q[4];\nccz q[0],q[1],q[2]; // first\n/* second */ ccz q[3],q[2],q[1];\n";
        let parsed = parse(qasm).unwrap();

        assert_eq!(parsed.gates.len(), 2);
        assert!(parsed.has_ccz());
        assert!(matches!(
            parsed.gates[1],
            Gate::ccz {
                control1: 3,
                control2: 2,
                target: 1
            }
        ));

        let reparsed = parse(&serialize(&parsed)).unwrap();
        assert_eq!(reparsed.gates.len(), 2);
        assert!(
            reparsed
                .gates
                .iter()
                .all(|gate| matches!(gate, Gate::ccz { .. }))
        );
    }

    #[test]
    fn multi_register_cz() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[1];\ncz a[0],b[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 2);
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(
            &c.gates[0],
            Gate::cz {
                control: 0,
                target: 1
            }
        ));
    }

    #[test]
    fn cz_round_trip_stays_native() {
        let input = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\ncz q[2],q[0];\n";
        let parsed = parse(input).unwrap();
        assert_eq!(parsed.gates.len(), 1);
        assert!(matches!(
            parsed.gates[0],
            Gate::cz {
                control: 2,
                target: 0
            }
        ));

        let serialized = serialize(&parsed);
        assert!(serialized.contains("cz q[2],q[0];"));
        assert!(!serialized.lines().any(|line| line.starts_with("cx ")));
        assert!(!serialized.lines().any(|line| line.starts_with("h ")));

        let reparsed = parse(&serialized).unwrap();
        assert_eq!(reparsed.gates.len(), 1);
        assert!(matches!(
            reparsed.gates[0],
            Gate::cz {
                control: 2,
                target: 0
            }
        ));
    }

    #[test]
    fn serialize_programmatic_cz() {
        let mut c = Circuit::new(2);
        c.apply(Gate::cz {
            control: 1,
            target: 0,
        });
        let qasm = serialize(&c);
        assert!(qasm.ends_with("cz q[1],q[0];\n"));
    }

    #[test]
    fn cz_rejects_wrong_operand_count() {
        for (operands, count) in [("q[0]", 1), ("q[0],q[1],q[2]", 3)] {
            let qasm = format!("OPENQASM 2.0;\nqreg q[3];\ncz {operands};\n");
            let err = parse(&qasm).unwrap_err();
            assert!(
                err.contains(&format!("cz expects 2 qubit operands, got {count}")),
                "{err}"
            );
        }
    }

    #[test]
    fn multiple_cz_gates_on_one_line_stay_native() {
        let qasm = "OPENQASM 2.0;\nqreg q[3];\ncz q[0],q[1]; cz q[2],q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 2);
        assert!(matches!(
            c.gates[0],
            Gate::cz {
                control: 0,
                target: 1
            }
        ));
        assert!(matches!(
            c.gates[1],
            Gate::cz {
                control: 2,
                target: 0
            }
        ));
        let serialized = serialize(&c);
        assert_eq!(
            serialized
                .lines()
                .filter(|line| line.starts_with("cz "))
                .count(),
            2
        );
    }

    #[test]
    fn cz_with_comments_and_mixed_gates_round_trips() {
        let qasm = "OPENQASM 2.0;\nqreg q[3];\nh q[0];\ncz q[0],q[2]; // native\nt q[2];\n/* keep this one too */ cz q[1],q[2];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 4);
        assert!(matches!(
            c.gates[1],
            Gate::cz {
                control: 0,
                target: 2
            }
        ));
        assert!(matches!(
            c.gates[3],
            Gate::cz {
                control: 1,
                target: 2
            }
        ));
        let reparsed = parse(&serialize(&c)).unwrap();
        assert_eq!(reparsed.gates.len(), 4);
        assert!(matches!(
            reparsed.gates[1],
            Gate::cz {
                control: 0,
                target: 2
            }
        ));
        assert!(matches!(
            reparsed.gates[3],
            Gate::cz {
                control: 1,
                target: 2
            }
        ));
    }

    #[test]
    fn multi_register_rz() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[2];\nrz(pi/4) b[1];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 3);
        // b[1] -> offset 1 + index 1 = 2
        if let Gate::rz(theta, 2) = &c.gates[0] {
            assert!((theta.to_f64_lossy().unwrap() - PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz on qubit 2, got {:?}", c.gates[0]);
        }
    }

    #[test]
    fn single_qubit_registers() {
        // Common pattern: many size-1 registers (like qrisp output)
        let qasm = "OPENQASM 2.0;\nqreg r0[1];\nqreg r1[1];\nqreg r2[1];\nqreg r3[1];\n\
                     cx r0[0],r3[0];\nt r2[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 4);
        // r0[0]->0, r3[0]->3, r2[0]->2
        assert!(matches!(
            &c.gates[0],
            Gate::cnot {
                control: 0,
                target: 3
            }
        ));
        assert!(matches!(&c.gates[1], Gate::t(2)));
    }

    #[test]
    fn multi_register_all_single_qubit_gates() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[1];\n\
                     x a[0];\ns b[0];\nsdg a[0];\nz b[0];\ntdg a[0];\nt b[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 2);
        assert!(matches!(&c.gates[0], Gate::x(0)));
        assert!(matches!(&c.gates[1], Gate::s(1)));
        assert!(matches!(&c.gates[2], Gate::sdg(0)));
        assert!(matches!(&c.gates[3], Gate::z(1)));
        assert!(matches!(&c.gates[4], Gate::tdg(0)));
        assert!(matches!(&c.gates[5], Gate::t(1)));
    }

    #[test]
    fn qreg_after_gate_on_same_line_error() {
        // semicolon-separated: gate then qreg on one line
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nh a[0]; qreg b[1];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("qreg declaration after gate"));
    }

    #[test]
    fn register_index_exactly_at_boundary() {
        // a[2] is out of range for size-2 register (valid: 0, 1)
        let qasm = "OPENQASM 2.0;\nqreg a[2];\nh a[2];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("out of range"));
    }

    #[test]
    fn register_index_max_valid() {
        let qasm = "OPENQASM 2.0;\nqreg a[3];\nh a[2];\n";
        let c = parse(qasm).unwrap();
        assert!(matches!(&c.gates[0], Gate::h(2)));
    }

    #[test]
    fn unknown_register_in_cnot() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\ncx a[0],nosuch[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("unknown register"));
        assert!(err.contains("nosuch"));
    }
}

/// Explicit numeric export for consumers without symbolic pi parameters.
/// Fraction conversion is approximate and its error is not certified.
#[derive(Debug)]
pub struct NumericExport {
    pub qasm: String,
    pub uncertified_conversions: usize,
}

pub fn serialize_numeric_lossy(circuit: &Circuit) -> Result<NumericExport, String> {
    crate::angle_export::lossy(circuit)
}

/// Numeric-parameter QASM for SDKs, equivalent up to overall global phase.
/// Named gates preserve mathematical quarter turns separately from floats.
/// Unsupported symbolic angles fail; this function never approximates them.
pub fn serialize_sdk(circuit: &Circuit) -> Result<String, String> {
    crate::angle_export::sdk(circuit)
}
