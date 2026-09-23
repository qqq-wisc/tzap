//! OpenQASM 2.0 parser and serializer.
//!
//! Supports `h`, `x`, `z`, `s`, `sdg`, `t`, `tdg`, `rz`, `cx`, `ccx`, `ccz`,
//! `cz`, `measure`, `reset`, `qreg`, and `creg`. Classical conditionals
//! (`if`), custom gate definitions (`gate`), barriers, and `include` files
//! other than `qelib1.inc` (which is ignored) are not supported.

use std::borrow::Cow;

use smallvec::SmallVec;

use crate::circuit::{CBit, Circuit, Gate, Qubit};

/// Qubit operands of one gate line. An arity above three only arises from
/// malformed input, which `require_arity` rejects; `SmallVec` keeps the
/// common case off the heap, where a `Vec` cost one allocation per gate
/// parsed.
type Operands = SmallVec<[Qubit; 3]>;

/// Byte offsets of the `[` and `]` of a register subscript such as `q[12]`.
///
/// A single scan, rather than `str::find('[')` and `str::find(']')`: operands
/// are a handful of bytes, so constructing two searchers costs more than the
/// scan. Requiring `]` *after* `[` also drops a latent panic on input like
/// `q]x[0`, where the two independent searches produced a reversed range.
fn subscript(part: &str) -> Option<(usize, usize)> {
    let bytes = part.as_bytes();
    let open = bytes.iter().position(|&b| b == b'[')?;
    let close = bytes[open + 1..].iter().position(|&b| b == b']')? + open + 1;
    Some((open, close))
}

/// Offset of a `//` line comment. `str::find` on a two-byte needle builds a
/// substring searcher for every line of the file.
fn line_comment(line: &str) -> Option<usize> {
    let bytes = line.as_bytes();
    (1..bytes.len())
        .find(|&i| bytes[i] == b'/' && bytes[i - 1] == b'/')
        .map(|i| i - 1)
}

/// Parse a circuit from OpenQASM 2.0 source. See the [module docs](self) for
/// the supported subset. Unrecognized lines produce an `Err`.
pub fn parse(qasm: &str) -> Result<Circuit, String> {
    let qasm = strip_block_comments(qasm);
    let mut registers: Vec<(String, usize, usize)> = Vec::new(); // (name, offset, size)
    let mut cregisters: Vec<(String, usize, usize)> = Vec::new();
    let mut num_qubits: usize = 0;
    let mut num_cbits: usize = 0;
    let mut gates = Vec::new();
    let mut seen_gate = false;
    for (line_num, raw_line) in qasm.lines().enumerate() {
        let line_num = line_num + 1;
        for line in raw_line
            .split(';')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
        {
            // strip inline comments
            let line = match line_comment(line) {
                Some(pos) => line[..pos].trim(),
                None => line,
            };
            // Group the statement prefixes by first byte. Every supported
            // statement is distinguished by its first byte, so a `t q[0];` line
            // tests one prefix instead of walking fifteen failing
            // `strip_prefix` calls to reach its arm — on a file whose every
            // line is a gate, that chain was most of the dispatch cost. Order
            // within each group is the original chain's order, which matters
            // wherever one prefix is a prefix of another.
            let Some(&first) = line.as_bytes().first() else {
                continue;
            };
            let candidates: &[&str] = match first {
                b'/' | b'O' | b'i' | b'b' => {
                    if line.starts_with("//")
                        || line.starts_with("OPENQASM")
                        || line.starts_with("include")
                        || line.starts_with("barrier")
                    {
                        continue;
                    }
                    &[]
                }
                b'q' => &["qreg"],
                b'c' => &["creg", "cx ", "ccz ", "ccx ", "cz "],
                b'm' => &["measure "],
                b'r' => &["reset ", "rz("],
                b'h' => &["h "],
                b'x' => &["x "],
                b's' => &["s ", "sdg "],
                b't' => &["tdg ", "t "],
                b'z' => &["z "],
                _ => &[],
            };
            let Some((keyword, rest)) = candidates
                .iter()
                .find_map(|&keyword| line.strip_prefix(keyword).map(|rest| (keyword, rest)))
            else {
                return Err(format!("line {line_num}: unsupported: {line}"));
            };
            match keyword {
                "qreg" => {
                    if seen_gate {
                        return Err(format!("line {line_num}: qreg declaration after gate"));
                    }
                    // parse "qreg name[size]"
                    let rest = rest.trim();
                    if let (Some(bracket), Some(end)) = (rest.find('['), rest.find(']')) {
                        let name = rest[..bracket].trim().to_string();
                        let size: usize = rest[bracket + 1..end]
                            .parse()
                            .map_err(|e| format!("line {line_num}: bad qreg size: {e}"))?;
                        registers.push((name, num_qubits, size));
                        num_qubits += size;
                        // Qubit indices are `u32` (see `circuit::Qubit`), so
                        // this is where a declaration too wide to index has to
                        // be rejected -- every later `as Qubit` is a cast on
                        // an index this bound has already cleared.
                        if num_qubits > Qubit::MAX as usize {
                            return Err(format!(
                                "line {line_num}: circuit declares {num_qubits} qubits, more than the {} supported",
                                Qubit::MAX
                            ));
                        }
                    }
                }
                "creg" => {
                    if seen_gate {
                        return Err(format!("line {line_num}: creg declaration after gate"));
                    }
                    let rest = rest.trim();
                    if let (Some(bracket), Some(end)) = (rest.find('['), rest.find(']')) {
                        let name = rest[..bracket].trim().to_string();
                        let size: usize = rest[bracket + 1..end]
                            .parse()
                            .map_err(|e| format!("line {line_num}: bad creg size: {e}"))?;
                        cregisters.push((name, num_cbits, size));
                        num_cbits += size;
                        // See the qreg case: `CBit` is `u32` too.
                        if num_cbits > CBit::MAX as usize {
                            return Err(format!(
                                "line {line_num}: circuit declares {num_cbits} classical bits, more than the {} supported",
                                CBit::MAX
                            ));
                        }
                    }
                }
                "measure " => {
                    seen_gate = true;
                    for (qubit, cbit) in parse_measure(rest, &registers, &cregisters, line_num)? {
                        gates.push(Gate::measure { qubit, cbit });
                    }
                }
                "reset " => {
                    seen_gate = true;
                    for q in expand_qubit_operand(rest, &registers, line_num)? {
                        gates.push(Gate::reset(q));
                    }
                }
                "cx " => {
                    seen_gate = true;
                    let qubits = resolve_qubits(rest, &registers, line_num)?;
                    require_arity("cx", &qubits, 2, line_num)?;
                    gates.push(Gate::cnot {
                        control: qubits[0],
                        target: qubits[1],
                    });
                }
                "ccz " => {
                    seen_gate = true;
                    let qubits = resolve_qubits(rest, &registers, line_num)?;
                    require_arity("ccz", &qubits, 3, line_num)?;
                    gates.push(Gate::ccz {
                        control1: qubits[0],
                        control2: qubits[1],
                        target: qubits[2],
                    });
                }
                "ccx " => {
                    seen_gate = true;
                    let qubits = resolve_qubits(rest, &registers, line_num)?;
                    require_arity("ccx", &qubits, 3, line_num)?;
                    gates.push(Gate::ccx {
                        control1: qubits[0],
                        control2: qubits[1],
                        target: qubits[2],
                    });
                }
                "cz " => {
                    seen_gate = true;
                    let qubits = resolve_qubits(rest, &registers, line_num)?;
                    require_arity("cz", &qubits, 2, line_num)?;
                    gates.push(Gate::cz {
                        control: qubits[0],
                        target: qubits[1],
                    });
                }
                "h " => {
                    seen_gate = true;
                    gates.push(Gate::h(resolve_single_qubit(
                        "h", rest, &registers, line_num,
                    )?));
                }
                "x " => {
                    seen_gate = true;
                    gates.push(Gate::x(resolve_single_qubit(
                        "x", rest, &registers, line_num,
                    )?));
                }
                "s " => {
                    seen_gate = true;
                    gates.push(Gate::s(resolve_single_qubit(
                        "s", rest, &registers, line_num,
                    )?));
                }
                "tdg " => {
                    seen_gate = true;
                    gates.push(Gate::tdg(resolve_single_qubit(
                        "tdg", rest, &registers, line_num,
                    )?));
                }
                "z " => {
                    seen_gate = true;
                    gates.push(Gate::z(resolve_single_qubit(
                        "z", rest, &registers, line_num,
                    )?));
                }
                "sdg " => {
                    seen_gate = true;
                    gates.push(Gate::sdg(resolve_single_qubit(
                        "sdg", rest, &registers, line_num,
                    )?));
                }
                "t " => {
                    seen_gate = true;
                    gates.push(Gate::t(resolve_single_qubit(
                        "t", rest, &registers, line_num,
                    )?));
                }
                "rz(" => {
                    seen_gate = true;
                    let paren_end = find_matching_paren(rest).ok_or_else(|| {
                        format!("line {line_num}: rz missing closing ')': {line}")
                    })?;
                    let theta = parse_angle(&rest[..paren_end], line_num)?;
                    let qubit =
                        resolve_single_qubit("rz", &rest[paren_end + 1..], &registers, line_num)?;
                    gates.push(Gate::rz(theta, qubit));
                }
                _ => unreachable!("keyword came from the candidate list"),
            }
        }
    }
    let mut c = Circuit::with_cbits(num_qubits, num_cbits);
    for g in gates {
        c.apply(g);
    }
    Ok(c)
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

fn write_gate(s: &mut String, gate: &Gate) {
    use std::fmt::Write;
    match gate {
        Gate::x(q) => writeln!(s, "x q[{q}];"),
        Gate::h(q) => writeln!(s, "h q[{q}];"),
        Gate::s(q) => writeln!(s, "s q[{q}];"),
        Gate::sdg(q) => writeln!(s, "sdg q[{q}];"),
        Gate::z(q) => writeln!(s, "z q[{q}];"),
        Gate::t(q) => writeln!(s, "t q[{q}];"),
        Gate::tdg(q) => writeln!(s, "tdg q[{q}];"),
        Gate::rz(theta, q) => writeln!(s, "rz({theta}) q[{q}];"),
        Gate::cnot { control, target } => writeln!(s, "cx q[{control}],q[{target}];"),
        Gate::cz { control, target } => writeln!(s, "cz q[{control}],q[{target}];"),
        Gate::ccx {
            control1,
            control2,
            target,
        } => {
            writeln!(s, "ccx q[{control1}],q[{control2}],q[{target}];")
        }
        Gate::ccz {
            control1,
            control2,
            target,
        } => {
            writeln!(s, "ccz q[{control1}],q[{control2}],q[{target}];")
        }
        Gate::measure { qubit, cbit } => writeln!(s, "measure q[{qubit}] -> c[{cbit}];"),
        Gate::reset(q) => writeln!(s, "reset q[{q}];"),
    }
    .unwrap();
}

fn require_arity(
    gate: &str,
    qubits: &[Qubit],
    expected: usize,
    line_num: usize,
) -> Result<(), String> {
    if qubits.len() == expected {
        Ok(())
    } else {
        let noun = if expected == 1 { "operand" } else { "operands" };
        Err(format!(
            "line {line_num}: {gate} expects {expected} qubit {noun}, got {}",
            qubits.len()
        ))
    }
}

/// Resolve a single-qubit gate's operand, rejecting anything but exactly one
/// qubit. Regression guard: `resolve_qubits(..)[0]` used to index straight
/// into the result, panicking on empty operands instead of returning a
/// parse error.
fn resolve_single_qubit(
    gate: &str,
    s: &str,
    registers: &[(String, usize, usize)],
    line_num: usize,
) -> Result<Qubit, String> {
    let qubits = resolve_qubits(s, registers, line_num)?;
    require_arity(gate, &qubits, 1, line_num)?;
    Ok(qubits[0])
}

/// Borrows when the source has no block comment at all — the overwhelmingly
/// common case, where copying the whole file (tens of megabytes on the larger
/// benchmarks) bought nothing.
fn strip_block_comments(s: &str) -> Cow<'_, str> {
    if !s.contains("/*") {
        return Cow::Borrowed(s);
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(start) = rest.find("/*") {
        out.push_str(&rest[..start]);
        match rest[start + 2..].find("*/") {
            Some(end) => {
                // preserve newlines so line numbers stay correct
                for c in rest[start..start + 2 + end + 2].chars() {
                    if c == '\n' {
                        out.push('\n');
                    }
                }
                rest = &rest[start + 2 + end + 2..];
            }
            None => {
                // unclosed block comment — treat rest as comment
                for c in rest[start..].chars() {
                    if c == '\n' {
                        out.push('\n');
                    }
                }
                return Cow::Owned(out);
            }
        }
    }
    out.push_str(rest);
    Cow::Owned(out)
}

/// Parse an angle expression with full arithmetic support.
/// Handles: numbers, `pi`, `+`, `-`, `*`, `/`, unary `-`, and parentheses.
fn parse_angle(s: &str, line_num: usize) -> Result<f64, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err(format!("line {line_num}: empty angle expression"));
    }
    let tokens = tokenize_angle(s).map_err(|e| format!("line {line_num}: {e}"))?;
    let mut pos = 0;
    let val = parse_expr(&tokens, &mut pos).map_err(|e| format!("line {line_num}: {e}"))?;
    if pos != tokens.len() {
        return Err(format!(
            "line {line_num}: unexpected token in angle expression"
        ));
    }
    Ok(val)
}

#[derive(Debug, Clone)]
enum Token {
    Num(f64),
    Pi,
    Plus,
    Minus,
    Star,
    Slash,
    LParen,
    RParen,
}

fn tokenize_angle(s: &str) -> Result<Vec<Token>, String> {
    let mut tokens = Vec::new();
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b' ' | b'\t' => i += 1,
            b'+' => {
                tokens.push(Token::Plus);
                i += 1;
            }
            b'-' => {
                tokens.push(Token::Minus);
                i += 1;
            }
            b'*' => {
                tokens.push(Token::Star);
                i += 1;
            }
            b'/' => {
                tokens.push(Token::Slash);
                i += 1;
            }
            b'(' => {
                tokens.push(Token::LParen);
                i += 1;
            }
            b')' => {
                tokens.push(Token::RParen);
                i += 1;
            }
            b'p' if s[i..].starts_with("pi")
                && (i + 2 >= bytes.len() || !bytes[i + 2].is_ascii_alphanumeric()) =>
            {
                tokens.push(Token::Pi);
                i += 2;
            }
            b'0'..=b'9' | b'.' => {
                let start = i;
                while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
                    i += 1;
                }
                // handle scientific notation e.g. 1e-10
                if i < bytes.len() && (bytes[i] == b'e' || bytes[i] == b'E') {
                    i += 1;
                    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
                        i += 1;
                    }
                    while i < bytes.len() && bytes[i].is_ascii_digit() {
                        i += 1;
                    }
                }
                let num: f64 = s[start..i]
                    .parse()
                    .map_err(|e| format!("bad number: {e}"))?;
                tokens.push(Token::Num(num));
            }
            _ => {
                return Err(format!(
                    "unexpected character '{}' in angle expression",
                    s[i..].chars().next().unwrap()
                ));
            }
        }
    }
    Ok(tokens)
}

// Recursive descent: expr = term (('+' | '-') term)*
fn parse_expr(tokens: &[Token], pos: &mut usize) -> Result<f64, String> {
    let mut val = parse_term(tokens, pos)?;
    while *pos < tokens.len() {
        match tokens[*pos] {
            Token::Plus => {
                *pos += 1;
                val += parse_term(tokens, pos)?;
            }
            Token::Minus => {
                *pos += 1;
                val -= parse_term(tokens, pos)?;
            }
            _ => break,
        }
    }
    Ok(val)
}

// term = unary (('*' | '/') unary)*
fn parse_term(tokens: &[Token], pos: &mut usize) -> Result<f64, String> {
    let mut val = parse_unary(tokens, pos)?;
    while *pos < tokens.len() {
        match tokens[*pos] {
            Token::Star => {
                *pos += 1;
                val *= parse_unary(tokens, pos)?;
            }
            Token::Slash => {
                *pos += 1;
                val /= parse_unary(tokens, pos)?;
            }
            _ => break,
        }
    }
    Ok(val)
}

// unary = '-' unary | atom
fn parse_unary(tokens: &[Token], pos: &mut usize) -> Result<f64, String> {
    if *pos < tokens.len() && matches!(tokens[*pos], Token::Minus) {
        *pos += 1;
        return Ok(-parse_unary(tokens, pos)?);
    }
    parse_atom(tokens, pos)
}

// atom = Num | Pi | '(' expr ')'
fn parse_atom(tokens: &[Token], pos: &mut usize) -> Result<f64, String> {
    if *pos >= tokens.len() {
        return Err("unexpected end of angle expression".to_string());
    }
    match &tokens[*pos] {
        Token::Num(n) => {
            let v = *n;
            *pos += 1;
            Ok(v)
        }
        Token::Pi => {
            *pos += 1;
            Ok(std::f64::consts::PI)
        }
        Token::LParen => {
            *pos += 1;
            let val = parse_expr(tokens, pos)?;
            if *pos >= tokens.len() {
                return Err("unclosed parenthesis in angle expression".to_string());
            }
            if let Token::RParen = tokens[*pos] {
                *pos += 1;
                Ok(val)
            } else {
                Err("expected ')' in angle expression".to_string())
            }
        }
        _ => Err("unexpected token in angle expression".to_string()),
    }
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
    cregisters: &[(String, usize, usize)],
    line_num: usize,
) -> Result<Vec<CBit>, String> {
    let mut result = Vec::new();
    for part in s.split(',') {
        let part = part.trim().trim_end_matches(';');
        if let (Some(bracket), Some(end)) = (part.find('['), part.find(']')) {
            let name = part[..bracket].trim();
            let idx: usize = part[bracket + 1..end]
                .parse()
                .map_err(|e| format!("line {line_num}: bad cbit index: {e}"))?;
            let (_, offset, size) = cregisters
                .iter()
                .find(|(n, _, _)| n == name)
                .ok_or_else(|| format!("line {line_num}: unknown classical register '{name}'"))?;
            if idx >= *size {
                return Err(format!(
                    "line {line_num}: index {idx} out of range for classical register '{name}' (size {size})"
                ));
            }
            result.push((offset + idx) as _);
        }
    }
    Ok(result)
}

fn resolve_qubits(
    s: &str,
    registers: &[(String, usize, usize)],
    line_num: usize,
) -> Result<Operands, String> {
    let mut result = Operands::new();
    for part in s.split(',') {
        let part = part.trim().trim_end_matches(';');
        if let Some((bracket, end)) = subscript(part) {
            let name = part[..bracket].trim();
            let idx: usize = part[bracket + 1..end]
                .parse()
                .map_err(|e| format!("line {line_num}: bad qubit index: {e}"))?;
            let (_, offset, size) = registers
                .iter()
                .find(|(n, _, _)| n == name)
                .ok_or_else(|| format!("line {line_num}: unknown register '{name}'"))?;
            if idx >= *size {
                return Err(format!(
                    "line {line_num}: index {idx} out of range for register '{name}' (size {size})"
                ));
            }
            result.push((offset + idx) as _);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
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
        c.apply(Gate::rz(PI / 4.0, 0));
        let qasm = serialize(&c);
        let c2 = parse(&qasm).unwrap();
        assert_eq!(c2.gates.len(), 1);
        if let Gate::rz(theta, 0) = &c2.gates[0] {
            assert!((theta - PI / 4.0).abs() < 1e-10);
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
    fn unclosed_block_comment_ignores_rest() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nh q[0];\n/* unclosed\nt q[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.gates.len(), 1);
        assert!(matches!(&c.gates[0], Gate::h(0)));
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
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nry(0.5) q[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 3"));
        assert!(err.contains("unsupported"));
        assert!(err.contains("ry"));
    }

    // --- pi expression tests ---

    #[test]
    fn rz_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_pi_over_4() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_2_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(2*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 2.0 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_3_pi_over_4() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(3*pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 3.0 * PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_neg_pi_over_2() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(-pi/2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - (-PI / 2.0)).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_neg_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(-pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - (-PI)).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_plain_float() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(0.123456789) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 0.123456789).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_20_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(20*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 20.0 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_pi_times_2() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi*2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 2.0 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_spaces_around_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz( pi / 4 ) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_spaces_coeff() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz( 3 * pi / 4 ) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 3.0 * PI / 4.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_neg_with_spaces() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(- pi / 4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - (-PI / 4.0)).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_pi_times_0_5() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(0.5*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 0.5 * PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_nested_parens() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz((pi/4)) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 4.0).abs() < 1e-10);
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
            assert!((theta - expected).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_addition() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi/4 + pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 2.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_subtraction() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi - pi/2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 2.0).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_double_neg() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(--pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI).abs() < 1e-10);
        } else {
            panic!("expected rz");
        }
    }

    #[test]
    fn rz_scientific_notation() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(1e-3) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 1e-3).abs() < 1e-15);
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
            assert!((theta - PI / 4.0).abs() < 1e-10);
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
