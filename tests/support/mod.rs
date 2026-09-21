//! Shared harness for the CLI-surface tests: a runner that pins the
//! environment, escape-sequence assertions, and a small JSON parser so the
//! `--json` report can be checked by value rather than by substring.
//!
//! Included with `#[path]` rather than as its own test target, so each test
//! binary gets its own copy and Cargo doesn't try to run this file's
//! (nonexistent) tests.
#![allow(dead_code)]

use std::io::Write;
use std::path::Path;
use std::process::{Child, Command, Output, Stdio};

/// Every environment variable that can change where tzap looks for its cache.
/// Cleared on every invocation so a developer's shell — an `XDG_CACHE_HOME`
/// pointing somewhere unexpected — can't change what these tests observe.
/// Tests that are *about* one of these set it back explicitly.
///
/// `LOCALAPPDATA` and `USERPROFILE` are on the list because they are the
/// tail of the same resolution order, and unlike the others they are set by
/// default on a Windows machine: without clearing them, the tests that pin
/// down what happens with *no* cache location would find one.
///
/// Nothing about tzap's *output* is environment-driven: styling follows from
/// whether the stream is a terminal and nothing else, which
/// `cli_streams::the_environment_does_not_change_what_is_printed` pins down
/// by setting the usual color conventions and checking they do nothing.
const PINNED_ENV: [&str; 4] = [
    "TZAP_CACHE_DIR",
    "XDG_CACHE_HOME",
    "LOCALAPPDATA",
    "USERPROFILE",
];

/// A tzap invocation under construction.
pub struct Tzap {
    command: Command,
    stdin: Option<String>,
}

impl Tzap {
    pub fn new(args: &[&str]) -> Tzap {
        let mut command = Command::new(env!("CARGO_BIN_EXE_tzap"));
        for name in PINNED_ENV {
            command.env_remove(name);
        }
        command.args(args);
        Tzap {
            command,
            stdin: None,
        }
    }

    pub fn env(mut self, name: &str, value: impl AsRef<std::ffi::OsStr>) -> Tzap {
        self.command.env(name, value);
        self
    }

    pub fn env_remove(mut self, name: &str) -> Tzap {
        self.command.env_remove(name);
        self
    }

    /// Feed `input` on stdin (for the `-` input path).
    pub fn stdin(mut self, input: &str) -> Tzap {
        self.stdin = Some(input.to_string());
        self
    }

    pub fn run(mut self) -> Run {
        let output = match self.stdin {
            None => self
                .command
                .stdin(Stdio::null())
                .output()
                .expect("failed to run tzap"),
            Some(input) => {
                let mut child = self
                    .command
                    .stdin(Stdio::piped())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .spawn()
                    .expect("failed to spawn tzap");
                child
                    .stdin
                    .as_mut()
                    .expect("piped stdin")
                    .write_all(input.as_bytes())
                    .expect("failed to write tzap's stdin");
                child.wait_with_output().expect("failed to run tzap")
            }
        };
        Run::from(output)
    }

    /// Start a command without waiting, for tests of cross-process races.
    pub fn spawn(mut self) -> Child {
        assert!(
            self.stdin.is_none(),
            "spawn does not support piped test input"
        );
        self.command
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("failed to spawn tzap")
    }
}

/// A finished invocation, with its streams decoded.
pub struct Run {
    pub status: std::process::ExitStatus,
    pub stdout: String,
    pub stderr: String,
}

impl From<Output> for Run {
    fn from(output: Output) -> Run {
        Run {
            status: output.status,
            stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
            stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        }
    }
}

impl Run {
    pub fn ok(self, context: &str) -> Run {
        assert!(
            self.status.success(),
            "{context}: expected success, got {:?}\nstderr:\n{}",
            self.status.code(),
            self.stderr
        );
        self
    }

    pub fn failed(self, context: &str) -> Run {
        assert!(
            !self.status.success(),
            "{context}: expected a non-zero exit\nstdout:\n{}\nstderr:\n{}",
            self.stdout,
            self.stderr
        );
        assert!(
            self.stdout.trim().is_empty(),
            "{context}: a failing run must leave stdout empty, got:\n{}",
            self.stdout
        );
        assert!(
            self.stderr.to_lowercase().contains("error"),
            "{context}: expected an error message on stderr, got:\n{}",
            self.stderr
        );
        self
    }

    /// Both streams, for the assertions that hold over the whole run.
    pub fn both(&self) -> String {
        format!("{}{}", self.stdout, self.stderr)
    }
}

/// Run tzap with `args`.
pub fn tzap(args: &[&str]) -> Run {
    Tzap::new(args).run()
}

/// The byte every ANSI escape sequence starts with — color and cursor motion
/// alike.
pub const ESC: char = '\x1b';

/// Assert that `text` contains no escape sequence and no bare carriage
/// return: the two things that make output unreadable when it lands in a
/// pipe, a log file, or a CI transcript.
pub fn assert_plain(text: &str, context: &str) {
    assert!(
        !text.contains(ESC),
        "{context}: expected no ANSI escapes, found one at byte {}:\n{}",
        text.find(ESC).unwrap(),
        text.escape_debug()
    );
    assert!(
        !text.contains('\r'),
        "{context}: expected no carriage returns (in-place redraw), got:\n{}",
        text.escape_debug()
    );
}

/// A QASM circuit's gate lines, header and blank lines dropped.
/// `text` with every timing figure (`0.123s`) replaced by a placeholder, so
/// two runs of the same command can be compared on everything *except* how
/// long they took. Timings are the only part of tzap's output that
/// legitimately differs between identical runs, and a millisecond of
/// scheduling noise is enough to make one print `0.000s` and the next
/// `0.001s`.
pub fn without_timings(text: &str) -> String {
    let bytes = text.as_bytes();
    let digits_at = |from: usize, count: usize| {
        bytes
            .get(from..from + count)
            .is_some_and(|run| run.iter().all(|byte| byte.is_ascii_digit()))
    };
    let mut out = String::with_capacity(text.len());
    let mut i = 0;
    while i < bytes.len() {
        let whole = bytes[i..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count();
        // `<digits>.<3 digits>s` — the shape every timing tzap prints has,
        // and one nothing else in its output shares.
        let point = i + whole;
        if whole > 0
            && bytes.get(point) == Some(&b'.')
            && digits_at(point + 1, 3)
            && bytes.get(point + 4) == Some(&b's')
        {
            out.push_str("N.NNNs");
            i = point + 5;
            continue;
        }
        if whole > 0 {
            out.push_str(&text[i..point]);
            i = point;
            continue;
        }
        let character = text[i..].chars().next().expect("i is a char boundary");
        out.push(character);
        i += character.len_utf8();
    }
    out
}

pub fn gate_lines(qasm: &str) -> Vec<String> {
    qasm.lines()
        .map(str::trim)
        .filter(|line| {
            !line.is_empty()
                && !line.starts_with("OPENQASM")
                && !line.starts_with("include")
                && !line.starts_with("qreg")
                && !line.starts_with("creg")
        })
        .map(str::to_string)
        .collect()
}

/// Assert `text` is a well-formed QASM 2.0 circuit, and return its gate
/// lines.
pub fn assert_valid_qasm(text: &str, context: &str) -> Vec<String> {
    assert!(
        text.starts_with("OPENQASM 2.0;"),
        "{context}: expected a QASM 2.0 header, got:\n{text}"
    );
    assert!(
        text.contains("qreg q["),
        "{context}: expected a qubit register, got:\n{text}"
    );
    assert_plain(text, &format!("{context}: QASM output"));
    gate_lines(text)
}

pub fn read(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
}

// ---------------------------------------------------------------------------
// A JSON parser, so `--json` can be asserted by value.
//
// tzap has no serde dependency (see `src/json.rs` on why), and the point of
// these tests is that the *output* is real JSON — checking it with a parser
// written independently of the writer is exactly the check that matters.
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    /// Parse `text` as one complete JSON document, panicking with context on
    /// anything malformed — including trailing content after the value, which
    /// is how a stream with output spliced into it would show up.
    pub fn parse(text: &str) -> Json {
        let chars: Vec<char> = text.chars().collect();
        let mut i = 0;
        let value = parse_value(&chars, &mut i, text);
        skip_whitespace(&chars, &mut i);
        assert!(
            i == chars.len(),
            "trailing content after the JSON document at char {i}: {:?}",
            text.chars().skip(i).take(60).collect::<String>()
        );
        value
    }

    pub fn opt(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(fields) => fields
                .iter()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value),
            _ => None,
        }
    }

    /// The value at `key`, panicking if this isn't an object or the key is
    /// absent — a missing key is a schema regression, not a `None`.
    pub fn get(&self, key: &str) -> &Json {
        self.opt(key)
            .unwrap_or_else(|| panic!("no key {key:?} in {self:?}"))
    }

    /// The value at a `/`-separated path, for the nested report.
    pub fn at(&self, path: &str) -> &Json {
        path.split('/').fold(self, |value, key| value.get(key))
    }

    pub fn arr(&self) -> &[Json] {
        match self {
            Json::Arr(items) => items,
            _ => panic!("not an array: {self:?}"),
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            Json::Str(s) => s,
            _ => panic!("not a string: {self:?}"),
        }
    }

    pub fn as_f64(&self) -> f64 {
        match self {
            Json::Num(n) => *n,
            _ => panic!("not a number: {self:?}"),
        }
    }

    pub fn as_usize(&self) -> usize {
        let n = self.as_f64();
        assert!(
            n >= 0.0 && n.fract() == 0.0,
            "not a non-negative integer: {n}"
        );
        n as usize
    }

    pub fn as_bool(&self) -> bool {
        match self {
            Json::Bool(b) => *b,
            _ => panic!("not a boolean: {self:?}"),
        }
    }

    pub fn is_null(&self) -> bool {
        *self == Json::Null
    }

    /// This array's strings, for the `options.passes` list.
    pub fn strings(&self) -> Vec<&str> {
        self.arr().iter().map(Json::as_str).collect()
    }

    pub fn keys(&self) -> Vec<&str> {
        match self {
            Json::Obj(fields) => fields.iter().map(|(name, _)| name.as_str()).collect(),
            _ => panic!("not an object: {self:?}"),
        }
    }
}

fn skip_whitespace(chars: &[char], i: &mut usize) {
    while matches!(chars.get(*i), Some(' ' | '\t' | '\n' | '\r')) {
        *i += 1;
    }
}

fn expect(chars: &[char], i: &mut usize, c: char, text: &str) {
    assert_eq!(
        chars.get(*i).copied(),
        Some(c),
        "expected {c:?} at char {i} of:\n{text}"
    );
    *i += 1;
}

fn parse_value(chars: &[char], i: &mut usize, text: &str) -> Json {
    skip_whitespace(chars, i);
    match chars.get(*i).copied() {
        Some('{') => {
            *i += 1;
            let mut fields = Vec::new();
            skip_whitespace(chars, i);
            if chars.get(*i) == Some(&'}') {
                *i += 1;
                return Json::Obj(fields);
            }
            loop {
                skip_whitespace(chars, i);
                let key = parse_string(chars, i, text);
                skip_whitespace(chars, i);
                expect(chars, i, ':', text);
                let value = parse_value(chars, i, text);
                assert!(
                    !fields.iter().any(|(name, _): &(String, Json)| *name == key),
                    "duplicate key {key:?} in the JSON object"
                );
                fields.push((key, value));
                skip_whitespace(chars, i);
                match chars.get(*i).copied() {
                    Some(',') => *i += 1,
                    Some('}') => {
                        *i += 1;
                        return Json::Obj(fields);
                    }
                    other => panic!("expected ',' or '}}' at char {i}, got {other:?}"),
                }
            }
        }
        Some('[') => {
            *i += 1;
            let mut items = Vec::new();
            skip_whitespace(chars, i);
            if chars.get(*i) == Some(&']') {
                *i += 1;
                return Json::Arr(items);
            }
            loop {
                items.push(parse_value(chars, i, text));
                skip_whitespace(chars, i);
                match chars.get(*i).copied() {
                    Some(',') => *i += 1,
                    Some(']') => {
                        *i += 1;
                        return Json::Arr(items);
                    }
                    other => panic!("expected ',' or ']' at char {i}, got {other:?}"),
                }
            }
        }
        Some('"') => Json::Str(parse_string(chars, i, text)),
        Some('t') => {
            consume_literal(chars, i, "true", text);
            Json::Bool(true)
        }
        Some('f') => {
            consume_literal(chars, i, "false", text);
            Json::Bool(false)
        }
        Some('n') => {
            consume_literal(chars, i, "null", text);
            Json::Null
        }
        Some(c) if c == '-' || c.is_ascii_digit() => parse_number(chars, i, text),
        other => panic!("unexpected {other:?} at char {i} of:\n{text}"),
    }
}

fn consume_literal(chars: &[char], i: &mut usize, literal: &str, text: &str) {
    for expected in literal.chars() {
        assert_eq!(
            chars.get(*i).copied(),
            Some(expected),
            "expected {literal:?} at char {i} of:\n{text}"
        );
        *i += 1;
    }
}

fn parse_string(chars: &[char], i: &mut usize, text: &str) -> String {
    expect(chars, i, '"', text);
    let mut out = String::new();
    loop {
        let c = chars
            .get(*i)
            .copied()
            .unwrap_or_else(|| panic!("unterminated string in:\n{text}"));
        *i += 1;
        match c {
            '"' => return out,
            '\\' => {
                let escape = chars
                    .get(*i)
                    .copied()
                    .unwrap_or_else(|| panic!("truncated escape in:\n{text}"));
                *i += 1;
                out.push(match escape {
                    '"' => '"',
                    '\\' => '\\',
                    '/' => '/',
                    'n' => '\n',
                    'r' => '\r',
                    't' => '\t',
                    'b' => '\u{8}',
                    'f' => '\u{c}',
                    'u' => {
                        let hex: String = chars[*i..*i + 4].iter().collect();
                        *i += 4;
                        let code = u32::from_str_radix(&hex, 16)
                            .unwrap_or_else(|_| panic!("bad \\u escape {hex:?}"));
                        char::from_u32(code)
                            .unwrap_or_else(|| panic!("bad code point {code} in \\u escape"))
                    }
                    other => panic!("invalid escape \\{other} in:\n{text}"),
                });
            }
            // A literal control character in a JSON string is invalid, and
            // is exactly what an unescaped path would produce.
            c if (c as u32) < 0x20 => {
                panic!("unescaped control character {:?} in a JSON string", c)
            }
            c => out.push(c),
        }
    }
}

fn parse_number(chars: &[char], i: &mut usize, text: &str) -> Json {
    let start = *i;
    if chars.get(*i) == Some(&'-') {
        *i += 1;
    }
    while matches!(chars.get(*i), Some(c) if c.is_ascii_digit()) {
        *i += 1;
    }
    if chars.get(*i) == Some(&'.') {
        *i += 1;
        while matches!(chars.get(*i), Some(c) if c.is_ascii_digit()) {
            *i += 1;
        }
    }
    if matches!(chars.get(*i), Some('e' | 'E')) {
        *i += 1;
        if matches!(chars.get(*i), Some('+' | '-')) {
            *i += 1;
        }
        while matches!(chars.get(*i), Some(c) if c.is_ascii_digit()) {
            *i += 1;
        }
    }
    let literal: String = chars[start..*i].iter().collect();
    Json::Num(
        literal
            .parse()
            .unwrap_or_else(|_| panic!("bad number {literal:?} in:\n{text}")),
    )
}
