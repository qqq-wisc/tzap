//! Parser for the scalar API only. It is independent of the QASM frontend.
use super::*;

struct Parser<'a, 'c> {
    source: &'a str,
    pos: usize,
    depth: usize,
    ctx: &'c mut AngleContext,
    budget: Budget,
}

impl AngleContext {
    pub fn parse(&mut self, source: &str) -> Result<Angle, Error> {
        let limits = self.limits();
        if source.len() > limits.max_source_bytes {
            return Err(Error::Limit(Resource::Source));
        }
        let mut budget = Budget::new(limits);
        budget.charge(source.len())?;
        let mut parser = Parser {
            source,
            pos: 0,
            depth: 0,
            ctx: self,
            budget,
        };
        let value = parser.expression(0)?;
        if parser.peek().is_some() {
            return Err(parser.error("unexpected token"));
        }
        parser.budget.value(&value)?;
        Ok(value)
    }
}

impl Parser<'_, '_> {
    fn error(&self, message: &'static str) -> Error {
        Error::Parse {
            offset: self.pos,
            message,
        }
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
    fn expect(&mut self, byte: u8, message: &'static str) -> Result<(), Error> {
        if self.peek() != Some(byte) {
            return Err(self.error(message));
        }
        self.pos += 1;
        Ok(())
    }
    fn expression(&mut self, min_binding: u8) -> Result<Angle, Error> {
        // Display spends at most four nested parser frames per expression
        // level (a sum's parentheses, `+` and `*` operands), so this frame
        // bound accepts every printed scalar. The value's own depth is still
        // checked against `max_depth`.
        self.depth += 1;
        if self.depth > self.budget.limits.max_depth.saturating_mul(4) {
            return Err(Error::Limit(Resource::Depth));
        }
        self.budget.charge(1)?;
        let mut lhs = match self.peek() {
            Some(b'-') => {
                self.pos += 1;
                let a = self.expression(5)?;
                self.ctx
                    .scale_inner(&a, &Rational::integer(-1), &mut self.budget)?
            }
            Some(b'+') => {
                self.pos += 1;
                self.expression(5)?
            }
            Some(b'(') => {
                self.pos += 1;
                let a = self.expression(0)?;
                self.expect(b')', "expected closing parenthesis")?;
                a
            }
            Some(b'0'..=b'9' | b'.') => self.number()?,
            Some(b'a'..=b'z' | b'A'..=b'Z' | b'_') => self.identifier()?,
            _ => return Err(self.error("expected a number, pi, or function")),
        };
        while let Some(op) = self.peek() {
            let (left, right) = match op {
                b'+' | b'-' => (1, 2),
                b'*' | b'/' => (3, 4),
                b'^' => (7, 6),
                _ => break,
            };
            if left < min_binding {
                break;
            }
            self.pos += 1;
            let rhs = self.expression(right)?;
            lhs = match op {
                b'+' => self.ctx.add_inner(&lhs, &rhs, &mut self.budget)?,
                b'-' => {
                    let rhs =
                        self.ctx
                            .scale_inner(&rhs, &Rational::integer(-1), &mut self.budget)?;
                    self.ctx.add_inner(&lhs, &rhs, &mut self.budget)?
                }
                b'*' => self.ctx.mul_inner(&lhs, &rhs, &mut self.budget)?,
                b'/' => self.ctx.div_inner(&lhs, &rhs, &mut self.budget)?,
                b'^' => self.ctx.power_inner(&lhs, &rhs, &mut self.budget)?,
                _ => unreachable!(),
            };
        }
        self.depth -= 1;
        Ok(lhs)
    }
    fn number(&mut self) -> Result<Angle, Error> {
        let start = self.pos;
        while self
            .source
            .as_bytes()
            .get(self.pos)
            .is_some_and(u8::is_ascii_digit)
        {
            self.pos += 1;
        }
        if self.source.as_bytes().get(self.pos) == Some(&b'.') {
            self.pos += 1;
            while self
                .source
                .as_bytes()
                .get(self.pos)
                .is_some_and(u8::is_ascii_digit)
            {
                self.pos += 1;
            }
        }
        if matches!(self.source.as_bytes().get(self.pos), Some(b'e' | b'E')) {
            self.pos += 1;
            if matches!(self.source.as_bytes().get(self.pos), Some(b'+' | b'-')) {
                self.pos += 1;
            }
            let exponent_start = self.pos;
            while self
                .source
                .as_bytes()
                .get(self.pos)
                .is_some_and(u8::is_ascii_digit)
            {
                self.pos += 1;
            }
            if exponent_start == self.pos {
                return Err(self.error("expected decimal exponent digits"));
            }
        }
        self.ctx
            .decimal_inner(&self.source[start..self.pos], &mut self.budget)
    }
    fn identifier(&mut self) -> Result<Angle, Error> {
        let start = self.pos;
        while self
            .source
            .as_bytes()
            .get(self.pos)
            .is_some_and(|b| b.is_ascii_alphanumeric() || *b == b'_')
        {
            self.pos += 1;
        }
        let name = &self.source[start..self.pos];
        if name == "pi" {
            return Ok(Angle::PI);
        }
        let fun = match name {
            "sin" => Function::Sin,
            "cos" => Function::Cos,
            "tan" => Function::Tan,
            "exp" => Function::Exp,
            "ln" => Function::Ln,
            "sqrt" => Function::Sqrt,
            _ => return Err(self.error("unsupported function or constant")),
        };
        self.expect(b'(', "expected function argument")?;
        let value = self.expression(0)?;
        self.expect(b')', "expected closing function parenthesis")?;
        self.ctx.function_inner(fun, &value, &mut self.budget)
    }
}
