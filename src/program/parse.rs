//! The OpenQASM 3 subset of Amy and Lunderville's control-flow benchmarks.

use std::collections::HashMap;
use std::fmt;

use super::{Program, Stmt};
use crate::circuit::{Circuit, Gate, Qubit};
use crate::decompose::DecomposeToffoli;
use crate::pass::Pass;

#[derive(Debug)]
pub struct ParseError(pub String);

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

type Result<T> = std::result::Result<T, ParseError>;

fn err<T>(msg: impl Into<String>) -> Result<T> {
    Err(ParseError(msg.into()))
}

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Ident(String),
    Num(f64),
    Str(String),
    Sym(&'static str),
}

const SYMS: [&str; 17] = [
    "->", "(", ")", "[", "]", "{", "}", ";", ",", ":", "+", "-", "*", "/", "=", "!", "<",
];

fn strip_comments(src: &str) -> String {
    let mut out = String::with_capacity(src.len());
    let mut rest = src;
    while !rest.is_empty() {
        if let Some(r) = rest.strip_prefix("//") {
            rest = r.find('\n').map_or("", |i| &r[i..]);
        } else if let Some(r) = rest.strip_prefix("/*") {
            rest = r.find("*/").map_or("", |i| &r[i + 2..]);
        } else {
            let c = rest.chars().next().unwrap();
            out.push(c);
            rest = &rest[c.len_utf8()..];
        }
    }
    out
}

fn tokenize(src: &str) -> Result<Vec<Tok>> {
    let src = strip_comments(src);
    let b = src.as_bytes();
    let mut toks = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let c = b[i] as char;
        if c.is_whitespace() {
            i += 1;
        } else if c.is_ascii_alphabetic() || c == '_' {
            let start = i;
            while i < b.len() && ((b[i] as char).is_ascii_alphanumeric() || b[i] == b'_') {
                i += 1;
            }
            toks.push(Tok::Ident(src[start..i].to_string()));
        } else if c.is_ascii_digit() || (c == '.' && i + 1 < b.len() && (b[i + 1] as char).is_ascii_digit()) {
            let start = i;
            while i < b.len() && ((b[i] as char).is_ascii_digit() || b[i] == b'.' || b[i] == b'e') {
                i += 1;
            }
            let n = src[start..i].parse().map_err(|_| ParseError(format!("bad number {}", &src[start..i])))?;
            toks.push(Tok::Num(n));
        } else if c == '"' {
            let end = src[i + 1..].find('"').ok_or_else(|| ParseError("unterminated string".into()))?;
            toks.push(Tok::Str(src[i + 1..i + 1 + end].to_string()));
            i += end + 2;
        } else if let Some(s) = SYMS.iter().find(|s| src[i..].starts_with(**s)) {
            toks.push(Tok::Sym(s));
            i += s.len();
        } else {
            // Other symbols only occur inside skipped conditions (e.g. `>`, `&`).
            toks.push(Tok::Sym("?"));
            i += 1;
        }
    }
    Ok(toks)
}

/// A user-defined gate: parameter names, qubit argument names, and body.
struct GateDef {
    params: Vec<String>,
    args: Vec<String>,
    body: Vec<Tok>,
}

/// How a qubit name resolves inside a gate body or the main program.
#[derive(Clone)]
enum Binding {
    /// A register: its first qubit and size.
    Reg(Qubit, u32),
    /// A gate's formal qubit argument.
    Qubit(Qubit),
}

struct Parser<'a> {
    toks: &'a [Tok],
    pos: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos)
    }

    fn next(&mut self) -> Result<Tok> {
        let t = self.toks.get(self.pos).cloned().ok_or_else(|| ParseError("unexpected end".into()))?;
        self.pos += 1;
        Ok(t)
    }

    fn eat(&mut self, s: &str) -> bool {
        if matches!(self.peek(), Some(Tok::Sym(x)) if *x == s) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn expect(&mut self, s: &str) -> Result<()> {
        if self.eat(s) {
            Ok(())
        } else {
            err(format!("expected `{s}`, found {:?}", self.peek()))
        }
    }

    fn ident(&mut self) -> Result<String> {
        match self.next()? {
            Tok::Ident(s) => Ok(s),
            t => err(format!("expected a name, found {t:?}")),
        }
    }

    /// Skip a balanced `( ... )`, `[ ... ]` or `{ ... }` group starting here.
    fn skip_group(&mut self) -> Result<Vec<Tok>> {
        let open = match self.next()? {
            Tok::Sym(s) => s,
            t => return err(format!("expected a bracket, found {t:?}")),
        };
        let close = match open {
            "(" => ")",
            "[" => "]",
            "{" => "}",
            _ => return err(format!("expected a bracket, found {open}")),
        };
        let start = self.pos;
        let mut depth = 1;
        while depth > 0 {
            match self.next()? {
                Tok::Sym(s) if s == open => depth += 1,
                Tok::Sym(s) if s == close => depth -= 1,
                _ => {}
            }
        }
        Ok(self.toks[start..self.pos - 1].to_vec())
    }

    /// Skip to just past the next `;` at bracket depth 0.
    fn skip_statement(&mut self) -> Result<()> {
        let mut depth = 0i32;
        loop {
            match self.next()? {
                Tok::Sym("(" | "[" | "{") => depth += 1,
                Tok::Sym(")" | "]" | "}") => depth -= 1,
                Tok::Sym(";") if depth == 0 => return Ok(()),
                _ => {}
            }
        }
    }
}

/// Integer expressions over `for` variables: `i`, `i+2`, `3*i-1`.
fn eval_int(toks: &[Tok], env: &HashMap<String, i64>) -> Result<i64> {
    fn term(t: &[Tok], i: &mut usize, env: &HashMap<String, i64>) -> Result<i64> {
        let mut v = atom(t, i, env)?;
        while *i < t.len() && t[*i] == Tok::Sym("*") {
            *i += 1;
            v *= atom(t, i, env)?;
        }
        Ok(v)
    }
    fn atom(t: &[Tok], i: &mut usize, env: &HashMap<String, i64>) -> Result<i64> {
        let tok = t.get(*i).ok_or_else(|| ParseError("empty expression".into()))?;
        *i += 1;
        match tok {
            Tok::Num(n) => Ok(*n as i64),
            Tok::Ident(s) => env.get(s).copied().ok_or_else(|| ParseError(format!("unknown variable {s}"))),
            Tok::Sym("-") => Ok(-atom(t, i, env)?),
            other => err(format!("bad expression token {other:?}")),
        }
    }
    let mut i = 0;
    let mut v = term(toks, &mut i, env)?;
    while i < toks.len() {
        match &toks[i] {
            Tok::Sym("+") => {
                i += 1;
                v += term(toks, &mut i, env)?;
            }
            Tok::Sym("-") => {
                i += 1;
                v -= term(toks, &mut i, env)?;
            }
            other => return err(format!("bad expression token {other:?}")),
        }
    }
    Ok(v)
}

/// Real expressions for rotation angles: numbers, `pi`, `+ - * /`.
fn eval_real(toks: &[Tok], params: &HashMap<String, f64>) -> Result<f64> {
    fn sum(t: &[Tok], i: &mut usize, p: &HashMap<String, f64>) -> Result<f64> {
        let mut v = prod(t, i, p)?;
        while let Some(Tok::Sym(s @ ("+" | "-"))) = t.get(*i) {
            *i += 1;
            let r = prod(t, i, p)?;
            v = if *s == "+" { v + r } else { v - r };
        }
        Ok(v)
    }
    fn prod(t: &[Tok], i: &mut usize, p: &HashMap<String, f64>) -> Result<f64> {
        let mut v = atom(t, i, p)?;
        while let Some(Tok::Sym(s @ ("*" | "/"))) = t.get(*i) {
            *i += 1;
            let r = atom(t, i, p)?;
            v = if *s == "*" { v * r } else { v / r };
        }
        Ok(v)
    }
    fn atom(t: &[Tok], i: &mut usize, p: &HashMap<String, f64>) -> Result<f64> {
        let tok = t.get(*i).cloned().ok_or_else(|| ParseError("empty angle".into()))?;
        *i += 1;
        match tok {
            Tok::Num(n) => Ok(n),
            Tok::Ident(s) if s == "pi" || s == "π" => Ok(std::f64::consts::PI),
            Tok::Ident(s) => p.get(&s).copied().ok_or_else(|| ParseError(format!("unknown parameter {s}"))),
            Tok::Sym("-") => Ok(-atom(t, i, p)?),
            Tok::Sym("(") => {
                let v = sum(t, i, p)?;
                if t.get(*i) != Some(&Tok::Sym(")")) {
                    return err("expected `)` in angle");
                }
                *i += 1;
                Ok(v)
            }
            other => err(format!("bad angle token {other:?}")),
        }
    }
    let mut i = 0;
    sum(toks, &mut i, params)
}

struct Ctx {
    num_qubits: u32,
    regs: HashMap<String, Binding>,
    defs: HashMap<String, GateDef>,
}

/// Lower a Toffoli-family gate into Clifford+T, as tzap's DecomposeToffoli.
fn lower(g: Gate, n: usize) -> Vec<Stmt> {
    let mut c = Circuit::with_cbits(n, 0);
    c.apply(g);
    DecomposeToffoli.run(&c).gates.into_iter().map(Stmt::Gate).collect()
}

impl Ctx {
    /// The qubits an operand names, under `scope` (gate arguments) and `env`.
    fn operand(&self, p: &mut Parser, scope: &HashMap<String, Binding>, env: &HashMap<String, i64>) -> Result<Vec<Qubit>> {
        let name = p.ident()?;
        let bind = scope
            .get(&name)
            .or_else(|| self.regs.get(&name))
            .cloned()
            .ok_or_else(|| ParseError(format!("unknown qubit {name}")))?;
        let index = if matches!(p.peek(), Some(Tok::Sym("["))) { Some(p.skip_group()?) } else { None };
        match (bind, index) {
            (Binding::Qubit(q), None) => Ok(vec![q]),
            (Binding::Reg(start, size), None) => Ok((start..start + size).collect()),
            (Binding::Reg(start, size), Some(ix)) => {
                if let Some(colon) = ix.iter().position(|t| *t == Tok::Sym(":")) {
                    let a = eval_int(&ix[..colon], env)?;
                    let b = eval_int(&ix[colon + 1..], env)?;
                    Ok((a..=b).map(|k| start + k as u32).collect())
                } else {
                    let k = eval_int(&ix, env)?;
                    if k < 0 || k as u32 >= size {
                        return err(format!("{name}[{k}] out of range"));
                    }
                    Ok(vec![start + k as u32])
                }
            }
            (Binding::Qubit(_), Some(_)) => err(format!("cannot index qubit {name}")),
        }
    }

    fn block(&self, p: &mut Parser, scope: &HashMap<String, Binding>, params: &HashMap<String, f64>, env: &HashMap<String, i64>) -> Result<Stmt> {
        if matches!(p.peek(), Some(Tok::Sym("{"))) {
            let body = p.skip_group()?;
            let mut inner = Parser { toks: &body, pos: 0 };
            let mut stmts = Vec::new();
            while inner.peek().is_some() {
                stmts.push(self.stmt(&mut inner, scope, params, env)?);
            }
            Ok(Stmt::Seq(stmts))
        } else {
            self.stmt(p, scope, params, env)
        }
    }

    fn stmt(&self, p: &mut Parser, scope: &HashMap<String, Binding>, params: &HashMap<String, f64>, env: &HashMap<String, i64>) -> Result<Stmt> {
        let word = match p.peek() {
            Some(Tok::Ident(w)) => w.clone(),
            Some(Tok::Sym("{")) => return self.block(p, scope, params, env),
            other => return err(format!("unexpected {other:?}")),
        };
        match word.as_str() {
            "while" => {
                p.next()?;
                p.skip_group()?;
                Ok(Stmt::While(Box::new(self.block(p, scope, params, env)?)))
            }
            "if" => {
                p.next()?;
                p.skip_group()?;
                let a = self.block(p, scope, params, env)?;
                let b = if matches!(p.peek(), Some(Tok::Ident(w)) if w == "else") {
                    p.next()?;
                    self.block(p, scope, params, env)?
                } else {
                    Stmt::Seq(Vec::new())
                };
                Ok(Stmt::If(Box::new(a), Box::new(b)))
            }
            "for" => {
                // for uint i in [a:b] body
                p.next()?;
                let mut var = p.ident()?;
                if var == "uint" || var == "int" {
                    var = p.ident()?;
                }
                if p.ident()? != "in" {
                    return err("expected `in`");
                }
                let range = p.skip_group()?;
                let colon = range.iter().position(|t| *t == Tok::Sym(":")).ok_or_else(|| ParseError("expected a range".into()))?;
                let a = eval_int(&range[..colon], env)?;
                let b = eval_int(&range[colon + 1..], env)?;
                let body_start = p.pos;
                let mut stmts = Vec::new();
                for k in a..=b {
                    p.pos = body_start;
                    let mut env2 = env.clone();
                    env2.insert(var.clone(), k);
                    stmts.push(self.block(p, scope, params, &env2)?);
                }
                if a > b {
                    self.block(p, scope, params, env)?;
                }
                Ok(Stmt::Seq(stmts))
            }
            "reset" => {
                p.next()?;
                let qs = self.operand(p, scope, env)?;
                p.expect(";")?;
                Ok(Stmt::Seq(qs.into_iter().map(Stmt::Reset).collect()))
            }
            "measure" => {
                p.next()?;
                let qs = self.operand(p, scope, env)?;
                p.skip_statement()?;
                Ok(Stmt::Seq(qs.into_iter().map(Stmt::Measure).collect()))
            }
            _ => {
                // `c = measure q;`
                if matches!(p.toks.get(p.pos + 1), Some(Tok::Sym("=" | "["))) {
                    let save = p.pos;
                    p.skip_statement()?;
                    let stmt = &p.toks[save..p.pos];
                    if let Some(m) = stmt.iter().position(|t| *t == Tok::Ident("measure".into())) {
                        let mut sub = Parser { toks: &stmt[m + 1..], pos: 0 };
                        let qs = self.operand(&mut sub, scope, env)?;
                        return Ok(Stmt::Seq(qs.into_iter().map(Stmt::Measure).collect()));
                    }
                    return Ok(Stmt::Seq(Vec::new()));
                }
                self.call(p, scope, params, env)
            }
        }
    }

    fn call(&self, p: &mut Parser, scope: &HashMap<String, Binding>, params: &HashMap<String, f64>, env: &HashMap<String, i64>) -> Result<Stmt> {
        let name = p.ident()?;
        let angles: Vec<f64> = if matches!(p.peek(), Some(Tok::Sym("("))) {
            let group = p.skip_group()?;
            group
                .split(|t| *t == Tok::Sym(","))
                .map(|e| eval_real(e, params))
                .collect::<Result<_>>()?
        } else {
            Vec::new()
        };
        let mut args = Vec::new();
        loop {
            let qs = self.operand(p, scope, env)?;
            args.push(qs);
            if !p.eat(",") {
                break;
            }
        }
        p.expect(";")?;
        // Broadcast register arguments (e.g. `h q;` on a register).
        let width = args.iter().map(Vec::len).max().unwrap_or(1);
        let mut out = Vec::new();
        for k in 0..width {
            let qs: Vec<Qubit> = args.iter().map(|a| if a.len() == 1 { a[0] } else { a[k] }).collect();
            out.push(self.apply(&name, &angles, &qs)?);
        }
        Ok(Stmt::Seq(out))
    }

    fn apply(&self, name: &str, angles: &[f64], qs: &[Qubit]) -> Result<Stmt> {
        if let Some(def) = self.defs.get(name) {
            if def.args.len() != qs.len() || def.params.len() != angles.len() {
                return err(format!("wrong arity calling {name}"));
            }
            let scope: HashMap<String, Binding> =
                def.args.iter().cloned().zip(qs.iter().map(|&q| Binding::Qubit(q))).collect();
            let params: HashMap<String, f64> = def.params.iter().cloned().zip(angles.iter().copied()).collect();
            let mut p = Parser { toks: &def.body, pos: 0 };
            let mut stmts = Vec::new();
            while p.peek().is_some() {
                stmts.push(self.stmt(&mut p, &scope, &params, &HashMap::new())?);
            }
            return Ok(Stmt::Seq(stmts));
        }
        let n = self.num_qubits as usize;
        let one = |f: fn(Qubit) -> Gate| -> Result<Stmt> {
            match qs {
                [q] => Ok(Stmt::Gate(f(*q))),
                _ => err(format!("{name} takes one qubit")),
            }
        };
        match (name, qs) {
            ("h", _) => one(Gate::h),
            ("x", _) => one(Gate::x),
            ("z", _) => one(Gate::z),
            ("s", _) => one(Gate::s),
            ("sdg", _) => one(Gate::sdg),
            ("t", _) => one(Gate::t),
            ("tdg", _) => one(Gate::tdg),
            ("id", _) => Ok(Stmt::Seq(Vec::new())),
            ("rz" | "p" | "u1", [q]) => Ok(Stmt::Gate(Gate::rz(angles[0], *q))),
            ("cx" | "CX" | "cnot", [c, t]) => Ok(Stmt::Gate(Gate::cnot { control: *c, target: *t })),
            ("cz", [c, t]) => Ok(Stmt::Gate(Gate::cz { control: *c, target: *t })),
            ("swap", [a, b]) => Ok(Stmt::Seq(vec![
                Stmt::Gate(Gate::cnot { control: *a, target: *b }),
                Stmt::Gate(Gate::cnot { control: *b, target: *a }),
                Stmt::Gate(Gate::cnot { control: *a, target: *b }),
            ])),
            ("ccx", [a, b, c]) => Ok(Stmt::Seq(lower(Gate::ccx { control1: *a, control2: *b, target: *c }, n))),
            ("ccz", [a, b, c]) => Ok(Stmt::Seq(lower(Gate::ccz { control1: *a, control2: *b, target: *c }, n))),
            _ => err(format!("unsupported gate {name} on {} qubits", qs.len())),
        }
    }
}

/// Flatten nested sequences.
pub(super) fn flatten(s: Stmt) -> Stmt {
    match s {
        Stmt::Seq(xs) => {
            let mut out = Vec::new();
            for x in xs {
                match flatten(x) {
                    Stmt::Seq(ys) => out.extend(ys),
                    y => out.push(y),
                }
            }
            if out.len() == 1 { out.pop().unwrap() } else { Stmt::Seq(out) }
        }
        Stmt::If(a, b) => Stmt::If(Box::new(flatten(*a)), Box::new(flatten(*b))),
        Stmt::While(b) => Stmt::While(Box::new(flatten(*b))),
        other => other,
    }
}

/// Parse a program. Declarations and gate definitions may appear anywhere
/// before use; everything else becomes the program body.
pub fn parse(src: &str) -> Result<Program> {
    let toks = tokenize(src)?;
    let mut ctx = Ctx { num_qubits: 0, regs: HashMap::new(), defs: HashMap::new() };
    let mut p = Parser { toks: &toks, pos: 0 };
    let mut body = Vec::new();
    let empty_scope = HashMap::new();
    let empty_env = HashMap::new();
    let empty_params = HashMap::new();
    while let Some(tok) = p.peek().cloned() {
        match tok {
            Tok::Ident(w) if w == "OPENQASM" || w == "include" || w == "bit" || w == "creg" || w == "input" || w == "output" => {
                p.skip_statement()?;
            }
            Tok::Ident(w) if w == "qubit" || w == "qreg" => {
                p.next()?;
                let (name, size) = if matches!(p.peek(), Some(Tok::Sym("["))) {
                    let n = eval_int(&p.skip_group()?, &empty_env)?;
                    (p.ident()?, n as u32)
                } else {
                    let name = p.ident()?;
                    let n = if matches!(p.peek(), Some(Tok::Sym("["))) { eval_int(&p.skip_group()?, &empty_env)? as u32 } else { 1 };
                    (name, n)
                };
                p.expect(";")?;
                ctx.regs.insert(name, Binding::Reg(ctx.num_qubits, size));
                ctx.num_qubits += size;
            }
            Tok::Ident(w) if w == "gate" => {
                p.next()?;
                let name = p.ident()?;
                let params = if matches!(p.peek(), Some(Tok::Sym("("))) {
                    p.skip_group()?
                        .into_iter()
                        .filter_map(|t| if let Tok::Ident(s) = t { Some(s) } else { None })
                        .collect()
                } else {
                    Vec::new()
                };
                let mut args = vec![p.ident()?];
                while p.eat(",") {
                    args.push(p.ident()?);
                }
                let body = p.skip_group()?;
                ctx.defs.insert(name, GateDef { params, args, body });
            }
            _ => body.push(ctx.stmt(&mut p, &empty_scope, &empty_params, &empty_env)?),
        }
    }
    Ok(Program { num_qubits: ctx.num_qubits as usize, body: flatten(Stmt::Seq(body)) })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_loops_branches_and_gate_definitions() {
        let src = r#"
            include "stdgates.inc";
            qubit a; qubit[2] r; bit[2] f = "11";
            gate g p, q { cx p, q; t q; }
            t a;
            while (int[2](f) != 0) { reset r[0]; g a, r[1]; measure r[0:1] -> f[0:1]; }
            if (true) { cx a, r[0]; } else { x a; }
            for uint i in [0:1] { h r[i]; }
            tdg a;
        "#;
        let prog = parse(src).unwrap();
        assert_eq!(prog.num_qubits, 3);
        assert_eq!(prog.t_count(), 3);
        let Stmt::Seq(top) = &prog.body else { panic!() };
        assert!(matches!(top[1], Stmt::While(_)));
        assert!(matches!(top[2], Stmt::If(_, _)));
    }
}
