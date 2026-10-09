//! Bounded contextual pi recognition, with ordinary binary64 evaluation for
//! pi-free input. Syntax is retained when affine recognition exceeds limits.
use crate::angle::{Angle, AngleError, FiniteF64, PiFraction, exact_float_add};
#[derive(Clone, Debug)]
enum Expr<'a> {
    Number(&'a str),
    Pi,
    Neg(Box<Expr<'a>>),
    Binary(u8, Box<Expr<'a>>, Box<Expr<'a>>),
}
struct Parser<'a> {
    source: &'a [u8],
    pos: usize,
    nodes: usize,
    depth: usize,
}
impl<'a> Parser<'a> {
    fn ws(&mut self) {
        while self
            .source
            .get(self.pos)
            .is_some_and(u8::is_ascii_whitespace)
        {
            self.pos += 1;
        }
    }
    fn peek(&mut self) -> Option<u8> {
        self.ws();
        self.source.get(self.pos).copied()
    }
    fn expression(&mut self, min: u8) -> Result<Expr<'a>, String> {
        self.depth += 1;
        self.nodes += 1;
        if self.depth > 128 || self.nodes > 256 {
            return Err("angle expression exceeds the parser complexity limit".into());
        }
        let mut left = match self.peek() {
            Some(b'-') => {
                self.pos += 1;
                Expr::Neg(Box::new(self.expression(3)?))
            }
            Some(b'(') => {
                self.pos += 1;
                let e = self.expression(0)?;
                if self.peek() != Some(b')') {
                    return Err("unclosed parenthesis in angle expression".into());
                }
                self.pos += 1;
                e
            }
            Some(b'p') if self.source[self.pos..].starts_with(b"pi") => {
                self.pos += 2;
                Expr::Pi
            }
            Some(b'0'..=b'9' | b'.') => {
                let start = self.pos;
                while self
                    .source
                    .get(self.pos)
                    .is_some_and(|b| b.is_ascii_digit() || *b == b'.')
                {
                    self.pos += 1;
                }
                if matches!(self.source.get(self.pos), Some(b'e' | b'E')) {
                    self.pos += 1;
                    if matches!(self.source.get(self.pos), Some(b'+' | b'-')) {
                        self.pos += 1;
                    }
                    while self.source.get(self.pos).is_some_and(u8::is_ascii_digit) {
                        self.pos += 1;
                    }
                }
                let s = std::str::from_utf8(&self.source[start..self.pos]).unwrap();
                s.parse::<f64>().map_err(|e| format!("bad number: {e}"))?;
                Expr::Number(s)
            }
            None => return Err("unexpected end of angle expression".into()),
            Some(c) => {
                return Err(format!(
                    "unexpected character '{}' in angle expression",
                    c as char
                ));
            }
        };
        loop {
            let Some(op) = self.peek() else { break };
            let precedence = match op {
                b'+' | b'-' => 1,
                b'*' | b'/' => 2,
                _ => break,
            };
            if precedence < min {
                break;
            }
            self.pos += 1;
            let right = self.expression(precedence + 1)?;
            self.nodes += 1;
            if self.nodes > 256 {
                return Err("angle expression exceeds the parser complexity limit".into());
            }
            left = Expr::Binary(op, Box::new(left), Box::new(right));
        }
        self.depth -= 1;
        Ok(left)
    }
}
impl Expr<'_> {
    fn has_pi(&self) -> bool {
        match self {
            Self::Pi => true,
            Self::Neg(e) => e.has_pi(),
            Self::Binary(_, a, b) => a.has_pi() || b.has_pi(),
            _ => false,
        }
    }
    fn numeric(&self) -> Result<f64, AngleError> {
        let v = match self {
            Self::Number(s) => s.parse::<f64>().unwrap(),
            Self::Pi => std::f64::consts::PI,
            Self::Neg(e) => -e.numeric()?,
            Self::Binary(op, a, b) => {
                let a = a.numeric()?;
                let b = b.numeric()?;
                match op {
                    b'+' => a + b,
                    b'-' => a - b,
                    b'*' => a * b,
                    b'/' => {
                        if b == 0.0 {
                            return Err(AngleError::DivisionByZero);
                        }
                        a / b
                    }
                    _ => unreachable!(),
                }
            }
        };
        FiniteF64::new(v).map(FiniteF64::get)
    }
    fn scalar(&self) -> Result<PiFraction, AngleError> {
        match self {
            Self::Number(s) => PiFraction::decimal(s),
            // A signed literal can represent i64::MIN even though its
            // unsigned magnitude cannot be stored as a positive coefficient.
            Self::Neg(e) => match &**e {
                Self::Number(s) => PiFraction::decimal(&format!("-{s}")),
                _ => e.scalar()?.checked_neg(),
            },
            Self::Binary(op, a, b) => {
                let a = a.scalar()?;
                let b = b.scalar()?;
                match op {
                    b'+' => a.checked_add(b),
                    b'-' => a.checked_add(b.checked_neg()?),
                    b'*' => a.checked_mul(b),
                    b'/' => a.checked_div(b),
                    _ => unreachable!(),
                }
            }
            Self::Pi => unreachable!(),
        }
    }
    fn affine(&self) -> Result<(PiFraction, FiniteF64), RecognizeError> {
        let zero = FiniteF64::new(0.0).unwrap();
        if !self.has_pi() {
            return Ok((PiFraction::ZERO, FiniteF64::new(self.numeric()?)?));
        }
        match self {
            Self::Pi => Ok((PiFraction::new(1, 1)?, zero)),
            Self::Neg(e) => {
                let (p, r) = e.affine()?;
                Ok((p.checked_neg()?, FiniteF64::new(-r.get())?))
            }
            Self::Binary(op, a, b) => match op {
                b'+' | b'-' => {
                    let (ap, ar) = a.affine()?;
                    let (bp, br) = b.affine()?;
                    let (bp, br) = if *op == b'-' {
                        (bp.checked_neg()?, FiniteF64::new(-br.get())?)
                    } else {
                        (bp, br)
                    };
                    Ok((ap.checked_add(bp)?, exact_float_add(ar, br)?))
                }
                b'*' | b'/' => {
                    if *op == b'*' && !a.has_pi() {
                        return Self::scale(b, a, false);
                    }
                    if !b.has_pi() {
                        return Self::scale(a, b, *op == b'/');
                    }
                    // A division by a proven exact zero must not fall back.
                    let (p, r) = b.affine()?;
                    if *op == b'/' && p.numerator() == 0 && r.get() == 0.0 {
                        return Err(AngleError::DivisionByZero.into());
                    }
                    Err(RecognizeError::NonAffine)
                }
                _ => unreachable!(),
            },
            _ => unreachable!(),
        }
    }
    fn scale(
        value: &Expr,
        scalar: &Expr,
        divide: bool,
    ) -> Result<(PiFraction, FiniteF64), RecognizeError> {
        // Validate pi-free arithmetic even if its rational value would cancel.
        scalar.numeric()?;
        let (p, r) = value.affine()?;
        let s = scalar.scalar()?;
        if divide && s.numerator() == 0 {
            return Err(AngleError::DivisionByZero.into());
        }
        if r.get() != 0.0 {
            return Err(RecognizeError::NonAffine);
        }
        Ok((
            if divide {
                p.checked_div(s)?
            } else {
                p.checked_mul(s)?
            },
            r,
        ))
    }
}
enum RecognizeError {
    Arithmetic(AngleError),
    NonAffine,
}
impl From<AngleError> for RecognizeError {
    fn from(e: AngleError) -> Self {
        Self::Arithmetic(e)
    }
}
pub(crate) fn parse(source: &str, line: usize) -> Result<Angle, String> {
    if source.trim().is_empty() {
        return Err(format!("line {line}: empty angle expression"));
    }
    let mut parser = Parser {
        source: source.as_bytes(),
        pos: 0,
        nodes: 0,
        depth: 0,
    };
    let expr = parser
        .expression(0)
        .map_err(|e| format!("line {line}: {e}"))?;
    if let Some(c) = parser.peek() {
        if !c.is_ascii_digit() && !b"p.+-*/()".contains(&c) {
            return Err(format!(
                "line {line}: unexpected character '{}' in angle expression",
                source[parser.pos..].chars().next().unwrap()
            ));
        }
        return Err(format!("line {line}: unexpected token in angle expression"));
    }
    let error = |e: AngleError| format!("line {line}: {e}");
    if !expr.has_pi() {
        return Angle::from_f64(expr.numeric().map_err(error)?).map_err(error);
    }
    match expr.affine() {
        Ok((p, r)) => Ok(Angle::affine(p, r)),
        Err(RecognizeError::Arithmetic(
            AngleError::RepresentationLimit | AngleError::RoundedAddition,
        )) => {
            expr.numeric().map_err(error)?;
            Ok(Angle::preserved(source.trim().to_string(), line))
        }
        Err(RecognizeError::Arithmetic(e)) => Err(error(e)),
        Err(RecognizeError::NonAffine) => {
            Angle::numeric_fallback(expr.numeric().map_err(error)?).map_err(error)
        }
    }
}
