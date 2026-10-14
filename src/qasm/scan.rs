//! Streaming statement boundaries. Most statements borrow their source text.
use std::borrow::Cow;

pub(super) struct Statement<'a> {
    pub text: Cow<'a, str>,
    pub line: usize,
    pub end: u8,
}

pub(super) struct Scanner<'a> {
    source: &'a str,
    pos: usize,
    line: usize,
}

impl<'a> Scanner<'a> {
    pub fn new(source: &'a str) -> Self {
        Self {
            source,
            pos: 0,
            line: 1,
        }
    }

    fn comment(&mut self) -> Result<bool, String> {
        let rest = &self.source[self.pos..];
        let length = if rest.starts_with("//") {
            rest.find('\n').unwrap_or(rest.len())
        } else if let Some(rest) = rest.strip_prefix("/*") {
            rest.find("*/")
                .map(|n| n + 4)
                .ok_or_else(|| format!("line {}: unterminated block comment", self.line))?
        } else {
            return Ok(false);
        };
        self.line += rest.as_bytes()[..length]
            .iter()
            .filter(|&&c| c == b'\n')
            .count();
        self.pos += length;
        Ok(true)
    }

    pub fn next(&mut self) -> Result<Option<Statement<'a>>, String> {
        let bytes = self.source.as_bytes();
        loop {
            while self.pos < bytes.len() && bytes[self.pos].is_ascii_whitespace() {
                self.line += usize::from(bytes[self.pos] == b'\n');
                self.pos += 1;
            }
            if !self.comment()? {
                break;
            }
        }
        if self.pos == bytes.len() {
            return Ok(None);
        }
        let start = self.pos;
        let line = self.line;
        // Flat inputs overwhelmingly have short statements without strings or
        // comments. Scan their boundaries in bulk, keeping the slow scanner for
        // quoted delimiters and comments embedded in statements.
        if let Some(length) = memchr::memchr3(b';', b'{', b'}', &bytes[start..]) {
            let end = start + length;
            if memchr::memchr2(b'/', b'"', &bytes[start..end]).is_none() {
                let text = &self.source[start..end];
                let c = bytes[end];
                self.line += memchr::memchr_iter(b'\n', &bytes[start..end]).count();
                self.pos = end + 1;
                if c == b'}' && !text.trim().is_empty() {
                    return Err(format!("line {line}: missing ';' before '}}'"));
                }
                return Ok(Some(Statement {
                    text: Cow::Borrowed(text),
                    line,
                    end: c,
                }));
            }
        }
        let mut copied = start;
        let mut owned: Option<String> = None;
        let mut quoted = false;
        while self.pos < bytes.len() {
            let c = bytes[self.pos];
            if !quoted && c == b'/' {
                let before = self.pos;
                if self.comment()? {
                    let out = owned.get_or_insert_with(String::new);
                    out.push_str(&self.source[copied..before]);
                    out.push(' ');
                    // Keep expression source line accounting through comments.
                    for &b in &bytes[before..self.pos] {
                        if b == b'\n' {
                            out.push('\n');
                        }
                    }
                    copied = self.pos;
                    continue;
                }
            }
            if c == b'"' {
                quoted = !quoted;
            }
            if quoted && c == b'\\' && self.pos + 1 < bytes.len() {
                self.pos += 1;
            } else if !quoted && matches!(c, b';' | b'{' | b'}') {
                let end = self.pos;
                self.pos += 1;
                let text = match owned {
                    Some(mut s) => {
                        s.push_str(&self.source[copied..end]);
                        Cow::Owned(s)
                    }
                    None => Cow::Borrowed(&self.source[start..end]),
                };
                if c == b'}' && !text.trim().is_empty() {
                    return Err(format!("line {line}: missing ';' before '}}'"));
                }
                return Ok(Some(Statement { text, line, end: c }));
            }
            self.line += usize::from(bytes[self.pos] == b'\n');
            self.pos += 1;
        }
        if quoted {
            return Err(format!("line {line}: unterminated string"));
        }
        let text = match owned {
            Some(mut s) => {
                s.push_str(&self.source[copied..]);
                Cow::Owned(s)
            }
            None => Cow::Borrowed(&self.source[start..]),
        };
        if !text.trim().is_empty() {
            return Err(format!("line {line}: missing ';' or closing '}}'"));
        }
        Ok(Some(Statement { text, line, end: 0 }))
    }
}

pub(super) fn identifier(s: &str) -> (&str, &str) {
    let end = s
        .bytes()
        .take_while(|c| c.is_ascii_alphanumeric() || *c == b'_')
        .count();
    (&s[..end], &s[end..])
}

/// OpenQASM 2 `nninteger`: ASCII digits only, so no sign is accepted.
pub(super) fn nninteger(s: &str) -> Option<usize> {
    let s = s.trim();
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse().ok()
}

pub(super) fn valid_name(s: &str) -> bool {
    s.as_bytes().first().is_some_and(u8::is_ascii_lowercase)
        && s.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_')
        && !matches!(
            s,
            "gate"
                | "opaque"
                | "qreg"
                | "creg"
                | "measure"
                | "reset"
                | "include"
                | "barrier"
                | "if"
                | "pi"
                | "sin"
                | "cos"
                | "tan"
                | "exp"
                | "ln"
                | "sqrt"
        )
}
