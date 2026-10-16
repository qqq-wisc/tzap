use super::{
    AngleParser, Evaluation, Operands, ParseLimits,
    expr::Expr,
    scan::{Scanner, identifier, valid_name},
};
use crate::{angle::Angle, circuit::Gate};
use smallvec::SmallVec;
use std::collections::HashMap;

pub(super) type Values = SmallVec<[Angle; 3]>;

#[derive(Clone, Copy, Debug)]
pub(super) enum Native {
    H,
    X,
    Z,
    S,
    Sdg,
    T,
    Tdg,
    Rz,
    Cx,
    Cz,
    Ccx,
    Ccz,
    U,
    P,
    Y,
    Sx,
    Rx,
    Ry,
    Swap,
    Cy,
    Ch,
    Cp,
    Crx,
    Cry,
    Crz,
    Cswap,
    U2,
    Id,
}

impl Native {
    pub fn lookup(name: &str) -> Option<Self> {
        Some(match name {
            "h" => Self::H,
            "x" => Self::X,
            "z" => Self::Z,
            "s" => Self::S,
            "sdg" => Self::Sdg,
            "t" => Self::T,
            "tdg" => Self::Tdg,
            "rz" => Self::Rz,
            "cx" | "CX" => Self::Cx,
            "cz" => Self::Cz,
            "ccx" => Self::Ccx,
            "ccz" => Self::Ccz,
            "U" | "u3" => Self::U,
            "p" | "u1" => Self::P,
            "y" => Self::Y,
            "sx" => Self::Sx,
            "rx" => Self::Rx,
            "ry" => Self::Ry,
            "swap" => Self::Swap,
            "cy" => Self::Cy,
            "ch" => Self::Ch,
            "cp" | "cu1" => Self::Cp,
            "crx" => Self::Crx,
            "cry" => Self::Cry,
            "crz" => Self::Crz,
            "cswap" => Self::Cswap,
            "u2" => Self::U2,
            "id" => Self::Id,
            _ => return None,
        })
    }
    /// A tzap shortcut that neither the language nor `qelib1.inc` declares,
    /// so a program may define the name itself (as Qiskit does for `ccz`).
    pub fn definable(self) -> bool {
        matches!(self, Self::Ccz)
    }
    pub fn signature(self) -> (usize, usize) {
        match self {
            Self::Rz | Self::P | Self::Rx | Self::Ry => (1, 1),
            Self::Cp | Self::Crx | Self::Cry | Self::Crz => (1, 2),
            Self::U => (3, 1),
            Self::U2 => (2, 1),
            Self::Cx | Self::Cz | Self::Swap | Self::Cy | Self::Ch => (0, 2),
            Self::Ccx | Self::Ccz | Self::Cswap => (0, 3),
            _ => (0, 1),
        }
    }
    fn emitted(self) -> usize {
        match self {
            Self::Id => 0,
            Self::U | Self::U2 => 3,
            _ => 1,
        }
    }
    pub fn emit(self, p: &[Angle], q: &[u32], gates: &mut Vec<Gate>) {
        let gate = match self {
            Self::H => Gate::h(q[0]),
            Self::X => Gate::x(q[0]),
            Self::Z => Gate::z(q[0]),
            Self::S => Gate::s(q[0]),
            Self::Sdg => Gate::sdg(q[0]),
            Self::T => Gate::t(q[0]),
            Self::Tdg => Gate::tdg(q[0]),
            Self::Rz => Gate::rz(p[0].clone(), q[0]),
            Self::P => Gate::p(p[0].clone(), q[0]),
            Self::Rx => Gate::rx(p[0].clone(), q[0]),
            Self::Ry => Gate::ry(p[0].clone(), q[0]),
            Self::Y => Gate::y(q[0]),
            Self::Sx => Gate::sx(q[0]),
            Self::Swap => Gate::swap(q[0], q[1]),
            Self::Cy => Gate::cy {
                control: q[0],
                target: q[1],
            },
            Self::Ch => Gate::ch {
                control: q[0],
                target: q[1],
            },
            Self::Cp => Gate::cp {
                lambda: p[0].clone(),
                control: q[0],
                target: q[1],
            },
            Self::Crx => Gate::crx {
                theta: p[0].clone(),
                control: q[0],
                target: q[1],
            },
            Self::Cry => Gate::cry {
                theta: p[0].clone(),
                control: q[0],
                target: q[1],
            },
            Self::Crz => Gate::crz {
                theta: p[0].clone(),
                control: q[0],
                target: q[1],
            },
            Self::Cswap => Gate::cswap {
                control: q[0],
                first: q[1],
                second: q[2],
            },
            Self::Cx => Gate::cnot {
                control: q[0],
                target: q[1],
            },
            Self::Cz => Gate::cz {
                control: q[0],
                target: q[1],
            },
            Self::Ccx => Gate::ccx {
                control1: q[0],
                control2: q[1],
                target: q[2],
            },
            Self::Ccz => Gate::ccz {
                control1: q[0],
                control2: q[1],
                target: q[2],
            },
            Self::Id => return,
            Self::U | Self::U2 => {
                let (theta, phi, lambda) = if matches!(self, Self::U2) {
                    (
                        Angle::pi_fraction(1, 2).expect("pi/2"),
                        p[0].clone(),
                        p[1].clone(),
                    )
                } else {
                    (p[0].clone(), p[1].clone(), p[2].clone())
                };
                gates.extend([
                    Gate::rz(lambda, q[0]),
                    Gate::ry(theta, q[0]),
                    Gate::rz(phi, q[0]),
                ]);
                return;
            }
        };
        gates.push(gate);
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) enum Callable {
    Native(Native),
    Custom(usize),
}

pub(super) struct BodyCall {
    pub callee: Callable,
    pub parameters: Vec<Expr>,
    pub qubits: SmallVec<[usize; 3]>,
    pub line: usize,
}

#[derive(Clone, Copy, Default)]
pub(super) struct Cost {
    pub gates: usize,
    pub work: usize,
    pub depth: usize,
}

pub(super) struct Definition {
    pub name: String,
    pub parameters: Vec<String>,
    pub qubits: Vec<String>,
    pub body: Vec<BodyCall>,
    pub line: usize,
    pub cost: Cost,
}

#[derive(Default)]
pub(super) struct Registry {
    pub definitions: Vec<Definition>,
    /// Every defined name. A definition shadows a definable native shortcut.
    names: HashMap<String, Callable>,
}

impl Registry {
    pub fn contains(&self, name: &str) -> bool {
        self.names.contains_key(name)
    }
    pub fn lookup(&self, name: &str, line: usize) -> Result<Callable, String> {
        if let Some(&callee) = self.names.get(name) {
            return Ok(callee);
        }
        Native::lookup(name).map(Callable::Native).ok_or_else(|| {
            format!("line {line}: unsupported: {name} (unknown gate or forward reference)")
        })
    }
    pub fn signature(&self, callee: Callable) -> (usize, usize) {
        match callee {
            Callable::Native(n) => n.signature(),
            Callable::Custom(i) => {
                let d = &self.definitions[i];
                (d.parameters.len(), d.qubits.len())
            }
        }
    }
    pub fn cost(&self, callee: Callable) -> Cost {
        match callee {
            Callable::Native(n) => Cost {
                gates: n.emitted(),
                work: 1 + n.signature().0 + n.signature().1,
                depth: 0,
            },
            Callable::Custom(i) => self.definitions[i].cost,
        }
    }
    pub fn define(
        &mut self,
        header: &str,
        line: usize,
        scanner: &mut Scanner<'_>,
        registers: &[(String, usize, usize)],
        cregisters: &[(String, usize, usize)],
        evaluation: Evaluation<'_>,
    ) -> Result<(), String> {
        let Evaluation { limits, angles } = evaluation;
        let (name, param_text, operand_text) = call_parts(header, line)?;
        if !valid_name(name) {
            return Err(format!("line {line}: invalid gate name '{name}'"));
        }
        let native = Native::lookup(name);
        if native.is_some_and(|native| !native.definable())
            || self.contains(name)
            || registers
                .iter()
                .chain(cregisters)
                .any(|(n, _, _)| n == name)
        {
            return Err(format!(
                "line {line}: duplicate or conflicting gate name '{name}'"
            ));
        }
        let parameters = names(param_text.unwrap_or(""), line)?;
        let qubits = names(operand_text, line)?;
        if qubits.is_empty() {
            return Err(format!(
                "line {line}: gate requires at least one qubit formal"
            ));
        }
        if parameters.iter().any(|p| qubits.contains(p)) {
            return Err(format!(
                "line {line}: overlapping parameter and qubit names"
            ));
        }
        let mut body = Vec::new();
        let mut cost = Cost {
            gates: 0,
            work: 1usize
                .saturating_add(parameters.len())
                .saturating_add(qubits.len()),
            depth: 1,
        };
        loop {
            let statement = scanner
                .next()?
                .ok_or_else(|| format!("line {line}: gate '{name}' missing closing '}}'"))?;
            if statement.end == b'}' {
                break;
            }
            if statement.end != b';' {
                return Err(format!(
                    "line {}: nested definitions are unsupported",
                    statement.line
                ));
            }
            let text = statement.text.trim();
            if text.is_empty() {
                continue;
            }
            let (callee_name, params, operands) = call_parts(text, statement.line)?;
            if callee_name == name {
                return Err(format!(
                    "line {}: recursive gate call '{name}'",
                    statement.line
                ));
            }
            if matches!(
                callee_name,
                "measure"
                    | "reset"
                    | "if"
                    | "barrier"
                    | "qreg"
                    | "creg"
                    | "gate"
                    | "include"
                    | "opaque"
            ) {
                return Err(format!(
                    "line {}: unsupported '{callee_name}' inside gate '{name}'",
                    statement.line
                ));
            }
            let callee = self.lookup(callee_name, statement.line)?;
            let expressions =
                expressions(params, Some(&parameters), limits, statement.line, angles)?;
            let mut args = SmallVec::<[usize; 3]>::new();
            for operand in comma_list(operands) {
                let i = qubits.iter().position(|q| q == operand).ok_or_else(|| {
                    format!(
                        "line {}: unknown or indexed qubit formal '{operand}' in gate '{name}'",
                        statement.line
                    )
                })?;
                args.push(i);
            }
            check_signature(
                callee_name,
                expressions.len(),
                args.len(),
                self.signature(callee),
                statement.line,
            )?;
            distinct(&args, callee_name, statement.line)?;
            let nested = self.cost(callee);
            cost.gates = cost.gates.saturating_add(nested.gates);
            cost.work = cost.work.saturating_add(nested.work);
            for expression in &expressions {
                cost.work = cost.work.saturating_add(expression.work());
            }
            cost.depth = cost.depth.max(nested.depth.saturating_add(1));
            body.push(BodyCall {
                callee,
                parameters: expressions,
                qubits: args,
                line: statement.line,
            });
        }
        // Keep the native gate only when the definition provably denotes it;
        // any other body is inlined like an ordinary custom gate.
        if let Some(native @ Native::Ccz) = native
            && parameters.is_empty()
            && qubits.len() == 3
            && is_ccz(&body)
        {
            self.names.insert(name.to_owned(), Callable::Native(native));
            return Ok(());
        }
        self.names
            .insert(name.to_owned(), Callable::Custom(self.definitions.len()));
        self.definitions.push(Definition {
            name: name.to_owned(),
            parameters,
            qubits,
            body,
            line,
            cost,
        });
        Ok(())
    }
}

/// Whether a three-qubit body is `h x; ccx y,z,x; h x` for some ordering of
/// its formals, which is CCZ: CCZ is symmetric in its three qubits.
fn is_ccz(body: &[BodyCall]) -> bool {
    let [first, middle, last] = body else {
        return false;
    };
    let hadamard = |call: &BodyCall| {
        matches!(call.callee, Callable::Native(Native::H)) && call.parameters.is_empty()
    };
    hadamard(first)
        && hadamard(last)
        && first.qubits == last.qubits
        && matches!(middle.callee, Callable::Native(Native::Ccx))
        && middle.parameters.is_empty()
        && middle.qubits.get(2) == first.qubits.first()
}

pub(super) fn comma_list(s: &str) -> impl Iterator<Item = &str> {
    s.split(',')
        .map(str::trim)
        .filter(move |_| !s.trim().is_empty())
}

fn names(s: &str, line: usize) -> Result<Vec<String>, String> {
    let mut result = Vec::new();
    for name in comma_list(s) {
        if !valid_name(name) {
            return Err(format!("line {line}: invalid formal name '{name}'"));
        }
        if result.iter().any(|n| n == name) {
            return Err(format!("line {line}: duplicate formal name '{name}'"));
        }
        result.push(name.to_owned());
    }
    Ok(result)
}

/// Return the name, optional parameter text, and remaining operands.
pub(super) fn call_parts(s: &str, line: usize) -> Result<(&str, Option<&str>, &str), String> {
    let (name, rest) = identifier(s.trim());
    if name.is_empty() {
        return Err(format!("line {line}: unsupported: {s}"));
    }
    let rest = rest.trim_start();
    if let Some(rest) = rest.strip_prefix('(') {
        let mut depth = 0usize;
        for (i, c) in rest.bytes().enumerate() {
            match c {
                b'(' => depth += 1,
                b')' if depth == 0 => return Ok((name, Some(&rest[..i]), rest[i + 1..].trim())),
                b')' => depth -= 1,
                _ => (),
            }
        }
        Err(format!("line {line}: {name} missing closing ')'"))
    } else {
        Ok((name, None, rest.trim()))
    }
}

pub(super) fn expressions(
    s: Option<&str>,
    parameters: Option<&[String]>,
    limits: ParseLimits,
    line: usize,
    angles: &mut AngleParser,
) -> Result<Vec<Expr>, String> {
    let Some(s) = s else {
        return Ok(Vec::new());
    };
    if s.trim().is_empty() {
        return Ok(Vec::new());
    }
    let mut result = Vec::new();
    let mut depth = 0usize;
    let mut start = 0;
    for (i, c) in s.bytes().enumerate() {
        match c {
            b'(' => depth += 1,
            b')' => {
                depth = depth
                    .checked_sub(1)
                    .ok_or_else(|| format!("line {line}: unexpected ')'"))?
            }
            b',' if depth == 0 => {
                result.push(Expr::parse(
                    &s[start..i],
                    parameters,
                    limits.expression_depth,
                    angles,
                    line,
                )?);
                start = i + 1;
            }
            _ => (),
        }
    }
    result.push(Expr::parse(
        &s[start..],
        parameters,
        limits.expression_depth,
        angles,
        line,
    )?);
    Ok(result)
}

pub(super) fn check_signature(
    name: &str,
    params: usize,
    qubits: usize,
    expected: (usize, usize),
    line: usize,
) -> Result<(), String> {
    if params != expected.0 {
        let message = if name == "rz" && params == 0 {
            "empty angle expression: expected a number, pi, or function; ".to_owned()
        } else {
            String::new()
        };
        return Err(format!(
            "line {line}: {message}{name} expects {} parameters, got {params}",
            expected.0
        ));
    }
    if qubits != expected.1 {
        let noun = if expected.1 == 1 {
            "operand"
        } else {
            "operands"
        };
        return Err(format!(
            "line {line}: {name} expects {} qubit {noun}, got {qubits}",
            expected.1
        ));
    }
    Ok(())
}

pub(super) fn distinct<T: PartialEq>(args: &[T], name: &str, line: usize) -> Result<(), String> {
    for (i, arg) in args.iter().enumerate() {
        if args[..i].contains(arg) {
            return Err(format!(
                "line {line}: {name} requires distinct qubit operands"
            ));
        }
    }
    Ok(())
}

pub(super) fn evaluate(
    expressions: &[Expr],
    values: &[Angle],
    line: usize,
    angles: &mut AngleParser,
) -> Result<Values, String> {
    expressions
        .iter()
        .map(|e| {
            e.evaluate(values, angles)
                .map_err(|e| format!("line {line}: {e}"))
        })
        .collect()
}

struct Frame {
    definition: usize,
    parameters: Values,
    qubits: Operands,
    next: usize,
    line: usize,
}

/// Iterative traversal keeps deeply nested user definitions off the Rust stack.
pub(super) fn expand(
    registry: &Registry,
    id: usize,
    parameters: Values,
    qubits: Operands,
    line: usize,
    gates: &mut Vec<Gate>,
    angles: &mut AngleParser,
) -> Result<(), String> {
    let mut stack = vec![Frame {
        definition: id,
        parameters,
        qubits,
        next: 0,
        line,
    }];
    while let Some(frame) = stack.last_mut() {
        let definition = &registry.definitions[frame.definition];
        let Some(call) = definition.body.get(frame.next) else {
            stack.pop();
            continue;
        };
        frame.next += 1;
        let p = match evaluate(&call.parameters, &frame.parameters, call.line, angles) {
            Ok(p) => p,
            Err(error) => {
                let mut context = error;
                for frame in stack.iter().rev().take(8) {
                    let d = &registry.definitions[frame.definition];
                    context.push_str(&format!(
                        "\n  called '{}' at line {} (defined at line {})",
                        d.name, frame.line, d.line
                    ));
                }
                return Err(context);
            }
        };
        let q: Operands = call.qubits.iter().map(|i| frame.qubits[*i]).collect();
        match call.callee {
            Callable::Native(n) => n.emit(&p, &q, gates),
            Callable::Custom(id) => stack.push(Frame {
                definition: id,
                parameters: p,
                qubits: q,
                next: 0,
                line: call.line,
            }),
        }
    }
    Ok(())
}
