//! Resolved postfix templates bind formal parameters to exact scalar angles.
use super::{AngleParser, parse_angle_with, parse_math_angle_with, scan::identifier};
use crate::angle::{Angle, AngleContext, AngleError};
use crate::symbolic_angle::{Function, Resource};

#[derive(Clone, Debug)]
enum Op {
    Value(Angle),
    Parameter(usize),
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Neg,
    Function(Function),
}

#[derive(Clone, Debug)]
pub(super) struct Expr {
    code: Vec<Op>,
    constant: Option<Angle>,
}

fn scalar(
    parser: &mut AngleParser,
    mut action: impl FnMut(&mut AngleContext) -> Result<Angle, AngleError>,
) -> Result<Angle, String> {
    let mut result = action(&mut parser.context);
    if matches!(result, Err(AngleError::Limit(Resource::Nodes))) {
        parser.context.clear_interner();
        result = action(&mut parser.context);
    }
    result.map_err(|e| e.to_string())
}
fn binary(parser: &mut AngleParser, op: &Op, a: &Angle, b: &Angle) -> Result<Angle, String> {
    scalar(parser, |ctx| match op {
        Op::Add => ctx.add(a, b),
        Op::Sub => ctx.sub(a, b),
        Op::Mul => ctx.mul(a, b),
        Op::Div => ctx.div(a, b),
        Op::Pow => ctx.power(a, b),
        _ => unreachable!(),
    })
}

impl Expr {
    pub fn work(&self) -> usize {
        if self.constant.is_some() {
            1
        } else {
            self.code.len()
        }
    }
    pub fn parse(
        source: &str,
        parameters: Option<&[String]>,
        max_depth: usize,
        angles: &mut AngleParser,
        line: usize,
    ) -> Result<Self, String> {
        // Only top-level actuals use the SDK's binary64 literal policy.
        let Some(parameters) = parameters else {
            return Ok(Self {
                code: Vec::new(),
                constant: Some(parse_angle_with(angles, source, line)?),
            });
        };
        // A body without formals still contains mathematical expressions.
        if parameters.is_empty() {
            return Ok(Self {
                code: Vec::new(),
                constant: Some(parse_math_angle_with(angles, source, line)?),
            });
        }
        if source.len() > angles.context.limits().max_source_bytes {
            return Err(format!(
                "line {line}: symbolic angle resource limit: Source"
            ));
        }
        let mut parser = Parser {
            source,
            pos: 0,
            parameters,
            max_depth,
            angles,
            line,
            code: Vec::new(),
        };
        let constant = parser
            .expression(0, 0)
            .map_err(|e| format!("line {line}: {e}"))?;
        if parser.peek().is_some() {
            return Err(format!("line {line}: unexpected token in angle expression"));
        }
        Ok(Self {
            code: parser.code,
            constant,
        })
    }
    pub fn evaluate(
        &self,
        parameters: &[Angle],
        angles: &mut AngleParser,
    ) -> Result<Angle, String> {
        if let Some(value) = &self.constant {
            return Ok(value.clone());
        }
        let mut stack = smallvec::SmallVec::<[Angle; 8]>::new();
        for op in &self.code {
            match op {
                Op::Value(v) => stack.push(v.clone()),
                Op::Parameter(i) => stack.push(parameters[*i].clone()),
                Op::Neg => {
                    let v = stack.last_mut().unwrap();
                    *v = scalar(angles, |ctx| ctx.neg(v))?;
                }
                Op::Function(fun) => {
                    let v = stack.last_mut().unwrap();
                    *v = scalar(angles, |ctx| ctx.function(*fun, v))?;
                }
                _ => {
                    let b = stack.pop().unwrap();
                    let a = stack.last_mut().unwrap();
                    *a = binary(angles, op, a, &b)?;
                }
            }
        }
        Ok(stack.pop().unwrap())
    }
}

struct Parser<'a, 'c> {
    source: &'a str,
    pos: usize,
    parameters: &'a [String],
    max_depth: usize,
    angles: &'c mut AngleParser,
    line: usize,
    code: Vec<Op>,
}
impl Parser<'_, '_> {
    fn literal(&mut self, source: &str) -> Result<Angle, String> {
        parse_math_angle_with(self.angles, source, self.line).map_err(|e| {
            e.strip_prefix(&format!("line {}: ", self.line))
                .unwrap_or(&e)
                .to_owned()
        })
    }
    fn peek(&mut self) -> Option<u8> {
        while self
            .source
            .as_bytes()
            .get(self.pos)
            .is_some_and(u8::is_ascii_whitespace)
        {
            self.pos += 1;
        }
        self.source.as_bytes().get(self.pos).copied()
    }
    fn expression(&mut self, binding: u8, depth: usize) -> Result<Option<Angle>, String> {
        if depth >= self.max_depth {
            return Err("expression nesting limit exceeded".into());
        }
        let mut value = match self.peek().ok_or("expected a number, pi, or function")? {
            b'-' => {
                self.pos += 1;
                let value = self
                    .expression(5, depth + 1)?
                    .map(|a| scalar(self.angles, |ctx| ctx.neg(&a)))
                    .transpose()?;
                self.code.push(Op::Neg);
                value
            }
            b'+' => {
                self.pos += 1;
                self.expression(5, depth + 1)?
            }
            b'(' => {
                self.pos += 1;
                let value = self.expression(0, depth + 1)?;
                if self.peek() != Some(b')') {
                    return Err("expected closing parenthesis".into());
                }
                self.pos += 1;
                value
            }
            b'0'..=b'9' | b'.' => {
                let start = self.pos;
                let bytes = self.source.as_bytes();
                while self.pos < bytes.len()
                    && (bytes[self.pos].is_ascii_digit() || bytes[self.pos] == b'.')
                {
                    self.pos += 1;
                }
                if self.pos < bytes.len() && matches!(bytes[self.pos], b'e' | b'E') {
                    self.pos += 1;
                    if self.pos < bytes.len() && matches!(bytes[self.pos], b'+' | b'-') {
                        self.pos += 1;
                    }
                    while self.pos < bytes.len() && bytes[self.pos].is_ascii_digit() {
                        self.pos += 1;
                    }
                }
                let a = self.literal(&self.source[start..self.pos])?;
                self.code.push(Op::Value(a.clone()));
                Some(a)
            }
            c if c.is_ascii_alphabetic() => {
                let (name, _) = identifier(&self.source[self.pos..]);
                self.pos += name.len();
                if let Some(i) = self.parameters.iter().position(|n| n == name) {
                    self.code.push(Op::Parameter(i));
                    None
                } else if name == "pi" {
                    let a = self.literal("pi")?;
                    self.code.push(Op::Value(a.clone()));
                    Some(a)
                } else {
                    let fun = match name {
                        "sin" => Function::Sin,
                        "cos" => Function::Cos,
                        "tan" => Function::Tan,
                        "exp" => Function::Exp,
                        "ln" => Function::Ln,
                        "sqrt" => Function::Sqrt,
                        _ => return Err(format!("unknown parameter '{name}' in angle expression")),
                    };
                    if self.peek() != Some(b'(') {
                        return Err("expected function argument".into());
                    }
                    self.pos += 1;
                    let value = self.expression(0, depth + 1)?;
                    if self.peek() != Some(b')') {
                        return Err("expected closing function parenthesis".into());
                    }
                    self.pos += 1;
                    let value = value
                        .map(|a| scalar(self.angles, |ctx| ctx.function(fun, &a)))
                        .transpose()?;
                    self.code.push(Op::Function(fun));
                    value
                }
            }
            _ => return Err("expected a number, pi, or function".into()),
        };
        while let Some(op) = self.peek() {
            let (left, right, op) = match op {
                b'+' => (1, 2, Op::Add),
                b'-' => (1, 2, Op::Sub),
                b'*' => (3, 4, Op::Mul),
                b'/' => (3, 4, Op::Div),
                b'^' => (7, 6, Op::Pow),
                _ => break,
            };
            if left < binding {
                break;
            }
            self.pos += 1;
            let rhs = self.expression(right, depth + 1)?;
            if matches!(op, Op::Div) && rhs.as_ref().is_some_and(|a| a == &Angle::ZERO) {
                return Err("angle expression must be finite: division by zero".into());
            }
            value = match (value, rhs) {
                (Some(a), Some(b)) => Some(binary(self.angles, &op, &a, &b)?),
                _ => None,
            };
            self.code.push(op);
        }
        Ok(value)
    }
}
