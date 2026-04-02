//! OpenQASM 2.0 parser and serializer.

use std::io::BufRead;

use crate::circuit::{Circuit, Gate};

pub fn parse(qasm: &str) -> Result<Circuit, String> {
    let qasm = strip_block_comments(qasm);
    let mut registers: Vec<(String, usize, usize)> = Vec::new(); // (name, offset, size)
    let mut num_qubits: usize = 0;
    let mut gates = Vec::new();
    let mut seen_gate = false;
    for (line_num, raw_line) in qasm.lines().enumerate() {
        let line_num = line_num + 1;
        for line in raw_line.split(';').map(|s| s.trim()).filter(|s| !s.is_empty()) {
        // strip inline comments
        let line = match line.find("//") {
            Some(pos) => line[..pos].trim(),
            None => line,
        };
        if line.is_empty()
            || line.starts_with("//")
            || line.starts_with("OPENQASM")
            || line.starts_with("include")
            || line.starts_with("barrier")
        {
            continue;
        }
        if line.starts_with("qreg") {
            if seen_gate {
                return Err(format!("line {line_num}: qreg declaration after gate"));
            }
            // parse "qreg name[size]"
            let rest = line[4..].trim();
            if let (Some(bracket), Some(end)) = (rest.find('['), rest.find(']')) {
                let name = rest[..bracket].trim().to_string();
                let size: usize = rest[bracket + 1..end].parse()
                    .map_err(|e| format!("line {line_num}: bad qreg size: {e}"))?;
                registers.push((name, num_qubits, size));
                num_qubits += size;
            }
        } else if let Some(rest) = line.strip_prefix("cx ") {
            seen_gate = true;
            let qubits = resolve_qubits(rest, &registers, line_num)?;
            gates.push(Gate::cnot { control: qubits[0], target: qubits[1] });
        } else if let Some(rest) = line.strip_prefix("ccx ") {
            seen_gate = true;
            let qubits = resolve_qubits(rest, &registers, line_num)?;
            gates.push(Gate::ccx { control1: qubits[0], control2: qubits[1], target: qubits[2] });
        } else if let Some(rest) = line.strip_prefix("cz ") {
            seen_gate = true;
            let qubits = resolve_qubits(rest, &registers, line_num)?;
            gates.push(Gate::h(qubits[1]));
            gates.push(Gate::cnot { control: qubits[0], target: qubits[1] });
            gates.push(Gate::h(qubits[1]));
        } else if let Some(rest) = line.strip_prefix("h ") {
            seen_gate = true;
            gates.push(Gate::h(resolve_qubits(rest, &registers, line_num)?[0]));
        } else if let Some(rest) = line.strip_prefix("x ") {
            seen_gate = true;
            gates.push(Gate::x(resolve_qubits(rest, &registers, line_num)?[0]));
        } else if let Some(rest) = line.strip_prefix("s ") {
            seen_gate = true;
            gates.push(Gate::s(resolve_qubits(rest, &registers, line_num)?[0]));
        } else if let Some(rest) = line.strip_prefix("tdg ") {
            seen_gate = true;
            gates.push(Gate::tdg(resolve_qubits(rest, &registers, line_num)?[0]));
        } else if let Some(rest) = line.strip_prefix("z ") {
            seen_gate = true;
            gates.push(Gate::z(resolve_qubits(rest, &registers, line_num)?[0]));
        } else if let Some(rest) = line.strip_prefix("sdg ") {
            seen_gate = true;
            gates.push(Gate::sdg(resolve_qubits(rest, &registers, line_num)?[0]));
        } else if let Some(rest) = line.strip_prefix("t ") {
            seen_gate = true;
            gates.push(Gate::t(resolve_qubits(rest, &registers, line_num)?[0]));
        } else if let Some(rest) = line.strip_prefix("rz(") {
            seen_gate = true;
            if let Some(paren_end) = find_matching_paren(rest) {
                let theta = parse_angle(&rest[..paren_end], line_num)?;
                let qubits = resolve_qubits(&rest[paren_end + 1..], &registers, line_num)?;
                gates.push(Gate::rz(theta, qubits[0]));
            }
        } else {
            return Err(format!("line {line_num}: unsupported: {line}"));
        }
        }
    }
    let mut c = Circuit::new(num_qubits);
    for g in gates {
        c.apply(g);
    }
    Ok(c)
}

pub fn serialize(circuit: &Circuit) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    writeln!(s, "OPENQASM 2.0;").unwrap();
    writeln!(s, "include \"qelib1.inc\";").unwrap();
    writeln!(s, "qreg q[{}];", circuit.num_qubits).unwrap();
    for gate in &circuit.gates {
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
            Gate::ccx { control1, control2, target } => {
                writeln!(s, "ccx q[{control1}],q[{control2}],q[{target}];")
            }
        }.unwrap();
    }
    s
}

/// Serialize just the gate lines (no header) for streaming output.
pub fn serialize_gates(gates: &[Gate]) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    for gate in gates {
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
            Gate::ccx { control1, control2, target } => {
                writeln!(s, "ccx q[{control1}],q[{control2}],q[{target}];")
            }
        }.unwrap();
    }
    s
}

fn parse_gate_stmt(
    line: &str,
    line_num: usize,
    registers: &[(String, usize, usize)],
    gates: &mut Vec<Gate>,
) -> Result<(), String> {
    if let Some(rest) = line.strip_prefix("cx ") {
        let qubits = resolve_qubits(rest, registers, line_num)?;
        gates.push(Gate::cnot { control: qubits[0], target: qubits[1] });
    } else if let Some(rest) = line.strip_prefix("ccx ") {
        let qubits = resolve_qubits(rest, registers, line_num)?;
        gates.push(Gate::ccx { control1: qubits[0], control2: qubits[1], target: qubits[2] });
    } else if let Some(rest) = line.strip_prefix("cz ") {
        let qubits = resolve_qubits(rest, registers, line_num)?;
        gates.push(Gate::h(qubits[1]));
        gates.push(Gate::cnot { control: qubits[0], target: qubits[1] });
        gates.push(Gate::h(qubits[1]));
    } else if let Some(rest) = line.strip_prefix("h ") {
        gates.push(Gate::h(resolve_qubits(rest, registers, line_num)?[0]));
    } else if let Some(rest) = line.strip_prefix("x ") {
        gates.push(Gate::x(resolve_qubits(rest, registers, line_num)?[0]));
    } else if let Some(rest) = line.strip_prefix("s ") {
        gates.push(Gate::s(resolve_qubits(rest, registers, line_num)?[0]));
    } else if let Some(rest) = line.strip_prefix("tdg ") {
        gates.push(Gate::tdg(resolve_qubits(rest, registers, line_num)?[0]));
    } else if let Some(rest) = line.strip_prefix("z ") {
        gates.push(Gate::z(resolve_qubits(rest, registers, line_num)?[0]));
    } else if let Some(rest) = line.strip_prefix("sdg ") {
        gates.push(Gate::sdg(resolve_qubits(rest, registers, line_num)?[0]));
    } else if let Some(rest) = line.strip_prefix("t ") {
        gates.push(Gate::t(resolve_qubits(rest, registers, line_num)?[0]));
    } else if let Some(rest) = line.strip_prefix("rz(") {
        if let Some(paren_end) = find_matching_paren(rest) {
            let theta = parse_angle(&rest[..paren_end], line_num)?;
            let qubits = resolve_qubits(&rest[paren_end + 1..], registers, line_num)?;
            gates.push(Gate::rz(theta, qubits[0]));
        }
    } else {
        return Err(format!("line {line_num}: unsupported: {line}"));
    }
    Ok(())
}

/// Streaming QASM reader that yields gates in batches without loading the whole file into memory.
pub struct StreamingReader<R: BufRead> {
    reader: R,
    pub num_qubits: usize,
    registers: Vec<(String, usize, usize)>,
    line_num: usize,
    in_block_comment: bool,
    done: bool,
    leftover: Vec<Gate>,
}

impl<R: BufRead> StreamingReader<R> {
    /// Create a streaming reader, consuming the QASM header (OPENQASM, include, qreg lines).
    pub fn new(mut reader: R) -> Result<Self, String> {
        let mut num_qubits = 0usize;
        let mut registers: Vec<(String, usize, usize)> = Vec::new();
        let mut line_num = 0;
        let mut in_block_comment = false;
        let mut line_buf = String::new();

        loop {
            line_buf.clear();
            let bytes = reader.read_line(&mut line_buf)
                .map_err(|e| format!("I/O error: {e}"))?;
            if bytes == 0 {
                if registers.is_empty() {
                    return Err("no qreg declaration found".into());
                }
                return Ok(StreamingReader {
                    reader, num_qubits, registers, line_num,
                    in_block_comment, done: true, leftover: Vec::new(),
                });
            }
            line_num += 1;

            let line = Self::strip_block_comment_line(&line_buf, &mut in_block_comment);
            let line = line.trim().to_string();
            if line.is_empty() { continue; }

            let mut leftover = Vec::new();
            let mut hit_gate = false;

            for stmt in line.split(';').map(|s| s.trim()).filter(|s| !s.is_empty()) {
                let stmt = match stmt.find("//") {
                    Some(pos) => stmt[..pos].trim(),
                    None => stmt,
                };
                if stmt.is_empty()
                    || stmt.starts_with("//")
                    || stmt.starts_with("OPENQASM")
                    || stmt.starts_with("include")
                    || stmt.starts_with("barrier")
                {
                    continue;
                }
                if !hit_gate && stmt.starts_with("qreg") {
                    let rest = stmt[4..].trim();
                    if let (Some(bracket), Some(end)) = (rest.find('['), rest.find(']')) {
                        let name = rest[..bracket].trim().to_string();
                        let size: usize = rest[bracket + 1..end].parse()
                            .map_err(|e| format!("line {line_num}: bad qreg size: {e}"))?;
                        registers.push((name, num_qubits, size));
                        num_qubits += size;
                    }
                    continue;
                }
                hit_gate = true;
                parse_gate_stmt(stmt, line_num, &registers, &mut leftover)?;
            }

            if hit_gate {
                if registers.is_empty() {
                    return Err("no qreg declaration found".into());
                }
                return Ok(StreamingReader {
                    reader, num_qubits, registers, line_num,
                    in_block_comment, done: false, leftover,
                });
            }
        }
    }

    /// Read the next batch of up to `batch_size` gates. Returns `None` at EOF.
    pub fn next_batch(&mut self, batch_size: usize) -> Result<Option<Vec<Gate>>, String> {
        if self.done {
            return Ok(None);
        }
        let mut gates = if !self.leftover.is_empty() {
            std::mem::take(&mut self.leftover)
        } else {
            Vec::with_capacity(batch_size.min(1_000_000))
        };
        let mut line_buf = String::new();

        while gates.len() < batch_size {
            line_buf.clear();
            let bytes = self.reader.read_line(&mut line_buf)
                .map_err(|e| format!("I/O error: {e}"))?;
            if bytes == 0 {
                self.done = true;
                break;
            }
            self.line_num += 1;
            let line_num = self.line_num;

            let line = Self::strip_block_comment_line(&line_buf, &mut self.in_block_comment);
            let line = line.trim().to_string();
            if line.is_empty() { continue; }

            for stmt in line.split(';').map(|s| s.trim()).filter(|s| !s.is_empty()) {
                let stmt = match stmt.find("//") {
                    Some(pos) => stmt[..pos].trim(),
                    None => stmt,
                };
                if stmt.is_empty()
                    || stmt.starts_with("//")
                    || stmt.starts_with("OPENQASM")
                    || stmt.starts_with("include")
                    || stmt.starts_with("barrier")
                    || stmt.starts_with("qreg")
                {
                    continue;
                }
                parse_gate_stmt(stmt, line_num, &self.registers, &mut gates)?;
            }
        }

        if gates.is_empty() { Ok(None) } else { Ok(Some(gates)) }
    }

    /// Strip block comments from a single line, tracking cross-line comment state.
    fn strip_block_comment_line(line: &str, in_comment: &mut bool) -> String {
        let mut out = String::with_capacity(line.len());
        let bytes = line.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            if *in_comment {
                if i + 1 < bytes.len() && bytes[i] == b'*' && bytes[i + 1] == b'/' {
                    *in_comment = false;
                    i += 2;
                } else {
                    i += 1;
                }
            } else if i + 1 < bytes.len() && bytes[i] == b'/' && bytes[i + 1] == b'*' {
                *in_comment = true;
                i += 2;
            } else {
                out.push(bytes[i] as char);
                i += 1;
            }
        }
        out
    }
}

fn strip_block_comments(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(start) = rest.find("/*") {
        out.push_str(&rest[..start]);
        match rest[start + 2..].find("*/") {
            Some(end) => {
                // preserve newlines so line numbers stay correct
                for c in rest[start..start + 2 + end + 2].chars() {
                    if c == '\n' { out.push('\n'); }
                }
                rest = &rest[start + 2 + end + 2..];
            }
            None => {
                // unclosed block comment — treat rest as comment
                for c in rest[start..].chars() {
                    if c == '\n' { out.push('\n'); }
                }
                return out;
            }
        }
    }
    out.push_str(rest);
    out
}

/// Parse an angle expression with full arithmetic support.
/// Handles: numbers, `pi`, `+`, `-`, `*`, `/`, unary `-`, and parentheses.
fn parse_angle(s: &str, line_num: usize) -> Result<f64, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err(format!("line {line_num}: empty angle expression"));
    }
    let tokens = tokenize_angle(s)
        .map_err(|e| format!("line {line_num}: {e}"))?;
    let mut pos = 0;
    let val = parse_expr(&tokens, &mut pos)
        .map_err(|e| format!("line {line_num}: {e}"))?;
    if pos != tokens.len() {
        return Err(format!("line {line_num}: unexpected token in angle expression"));
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
            b'+' => { tokens.push(Token::Plus); i += 1; }
            b'-' => { tokens.push(Token::Minus); i += 1; }
            b'*' => { tokens.push(Token::Star); i += 1; }
            b'/' => { tokens.push(Token::Slash); i += 1; }
            b'(' => { tokens.push(Token::LParen); i += 1; }
            b')' => { tokens.push(Token::RParen); i += 1; }
            b'p' if s[i..].starts_with("pi") && (i + 2 >= bytes.len() || !bytes[i + 2].is_ascii_alphanumeric()) => {
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
                let num: f64 = s[start..i].parse()
                    .map_err(|e| format!("bad number: {e}"))?;
                tokens.push(Token::Num(num));
            }
            _ => return Err(format!("unexpected character '{}' in angle expression", s[i..].chars().next().unwrap())),
        }
    }
    Ok(tokens)
}

// Recursive descent: expr = term (('+' | '-') term)*
fn parse_expr(tokens: &[Token], pos: &mut usize) -> Result<f64, String> {
    let mut val = parse_term(tokens, pos)?;
    while *pos < tokens.len() {
        match tokens[*pos] {
            Token::Plus => { *pos += 1; val += parse_term(tokens, pos)?; }
            Token::Minus => { *pos += 1; val -= parse_term(tokens, pos)?; }
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
            Token::Star => { *pos += 1; val *= parse_unary(tokens, pos)?; }
            Token::Slash => { *pos += 1; val /= parse_unary(tokens, pos)?; }
            _ => break,
        }
    }
    Ok(val)
}

// unary = '-' unary | atom
fn parse_unary(tokens: &[Token], pos: &mut usize) -> Result<f64, String> {
    if *pos < tokens.len() {
        if let Token::Minus = tokens[*pos] {
            *pos += 1;
            return Ok(-parse_unary(tokens, pos)?);
        }
    }
    parse_atom(tokens, pos)
}

// atom = Num | Pi | '(' expr ')'
fn parse_atom(tokens: &[Token], pos: &mut usize) -> Result<f64, String> {
    if *pos >= tokens.len() {
        return Err("unexpected end of angle expression".to_string());
    }
    match &tokens[*pos] {
        Token::Num(n) => { let v = *n; *pos += 1; Ok(v) }
        Token::Pi => { *pos += 1; Ok(std::f64::consts::PI) }
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

fn resolve_qubits(
    s: &str,
    registers: &[(String, usize, usize)],
    line_num: usize,
) -> Result<Vec<usize>, String> {
    let mut result = Vec::new();
    for part in s.split(',') {
        let part = part.trim().trim_end_matches(';');
        if let (Some(bracket), Some(end)) = (part.find('['), part.find(']')) {
            let name = part[..bracket].trim();
            let idx: usize = part[bracket + 1..end].parse()
                .map_err(|e| format!("line {line_num}: bad qubit index: {e}"))?;
            let (_, offset, size) = registers.iter()
                .find(|(n, _, _)| n == name)
                .ok_or_else(|| format!("line {line_num}: unknown register '{name}'"))?;
            if idx >= *size {
                return Err(format!(
                    "line {line_num}: index {idx} out of range for register '{name}' (size {size})"
                ));
            }
            result.push(offset + idx);
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
        c.apply(Gate::cnot { control: 0, target: 1 });
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
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n// just a comment\nqreg q[1];\nh q[0];\n";
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
        assert!(err.contains("line 7"), "expected line 7 in error, got: {err}");
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
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_pi_over_4() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 4.0).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_2_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(2*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 2.0 * PI).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_3_pi_over_4() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(3*pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 3.0 * PI / 4.0).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_neg_pi_over_2() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(-pi/2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - (-PI / 2.0)).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_neg_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(-pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - (-PI)).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_plain_float() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(0.785398163) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 0.785398163).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_20_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(20*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 20.0 * PI).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_pi_times_2() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi*2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 2.0 * PI).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_spaces_around_pi() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz( pi / 4 ) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 4.0).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_spaces_coeff() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz( 3 * pi / 4 ) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 3.0 * PI / 4.0).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_neg_with_spaces() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(- pi / 4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - (-PI / 4.0)).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_pi_times_0_5() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(0.5*pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 0.5 * PI).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_nested_parens() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz((pi/4)) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 4.0).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_complex_expr() {
        // 2*2*(pi/2)*(3/4*pi) = 4 * (pi/2) * (0.75*pi) = 4 * 0.5*pi * 0.75*pi = 3*pi^2/2
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(2*2*(pi/2)*(3/4*pi)) q[0];\n";
        let c = parse(qasm).unwrap();
        let expected = 2.0 * 2.0 * (PI / 2.0) * (3.0 / 4.0 * PI);
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - expected).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_addition() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi/4 + pi/4) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 2.0).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_subtraction() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(pi - pi/2) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI / 2.0).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_double_neg() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(--pi) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - PI).abs() < 1e-10);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn rz_scientific_notation() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\nrz(1e-3) q[0];\n";
        let c = parse(qasm).unwrap();
        if let Gate::rz(theta, 0) = &c.gates[0] {
            assert!((theta - 1e-3).abs() < 1e-15);
        } else { panic!("expected rz"); }
    }

    #[test]
    fn creg_error() {
        let qasm = "OPENQASM 2.0;\nqreg q[1];\ncreg c[1];\nh q[0];\n";
        let err = parse(qasm).unwrap_err();
        assert!(err.contains("line 3"));
        assert!(err.contains("creg"));
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
        assert!(matches!(&c.gates[0], Gate::cnot { control: 0, target: 1 }));
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
        assert!(matches!(&c.gates[0], Gate::ccx { control1: 0, control2: 1, target: 2 }));
    }

    #[test]
    fn multi_register_cz() {
        let qasm = "OPENQASM 2.0;\nqreg a[1];\nqreg b[1];\ncz a[0],b[0];\n";
        let c = parse(qasm).unwrap();
        assert_eq!(c.num_qubits, 2);
        // cz decomposes to h, cnot, h
        assert_eq!(c.gates.len(), 3);
        assert!(matches!(&c.gates[0], Gate::h(1)));
        assert!(matches!(&c.gates[1], Gate::cnot { control: 0, target: 1 }));
        assert!(matches!(&c.gates[2], Gate::h(1)));
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
        assert!(matches!(&c.gates[0], Gate::cnot { control: 0, target: 3 }));
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
