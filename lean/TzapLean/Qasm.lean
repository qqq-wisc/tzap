import TzapLean.Circuit

/-!
# OpenQASM 2.0

A port of `src/qasm.rs`: parser and serializer for the subset tzap reads — `h`, `x`, `z`, `s`,
`sdg`, `t`, `tdg`, `cx`, `ccx`, `ccz`, `cz`, `measure`, `reset`, `qreg`, and `creg`. Classical
conditionals (`if`), custom gate definitions (`gate`), and `include` files other than
`qelib1.inc` (which is ignored) are not supported.

**`rz` is rejected here**, with an error naming it, where Rust parses an angle expression into
an `f64`. This port is exact — angles are rationals in units of `π` — and the pass that would
consume an arbitrary `rz`, gridsynth, is not ported either, so accepting one would mean
carrying a gate nothing downstream can lower. Rejecting it says so at the door rather than
somewhere further in.

Nothing here is verified. It is a front end: it turns bytes into a `Circuit`, and every
theorem in this development is about what happens after that.
-/

namespace TzapLean

namespace Qasm

/-! ## Small string helpers -/

/-- Drop leading and trailing ASCII whitespace. -/
def trim (s : String) : String := s.trimAscii.toString

/-- Split on a character, dropping empty pieces after trimming. -/
def splitTrim (s : String) (c : Char) : List String :=
  (s.splitOn c.toString).map trim |>.filter (· ≠ "")

/-- Strip a trailing `;` if present. -/
def dropSemi (s : String) : String :=
  if s.endsWith ";" then (s.dropEnd 1).toString else s

/-- Offset of a `//` line comment. -/
def lineComment (line : String) : Option Nat :=
  let cs := line.toList
  let rec go (i : Nat) : List Char → Option Nat
    | [] => none
    | _ :: [] => none
    | a :: b :: rest => if a = '/' && b = '/' then some i else go (i + 1) (b :: rest)
  go 0 cs

/-- Remove `/* … */` comments, preserving newlines so line numbers stay correct. -/
partial def stripBlockComments (s : String) : String :=
  if !(s.splitOn "/*").length.beq 1 then
    let parts := s.splitOn "/*"
    match parts with
    | [] => s
    | head :: rest =>
        rest.foldl
          (fun acc chunk =>
            match chunk.splitOn "*/" with
            | [] => acc
            | [only] =>
                -- unclosed: the rest is comment, keep only its newlines
                acc ++ String.ofList (only.toList.filter (· = '\n'))
            | commented :: tailParts =>
                acc ++ String.ofList (commented.toList.filter (· = '\n')) ++
                  String.intercalate "*/" tailParts)
          head
  else s

/-- The `[` and `]` of a register subscript such as `q[12]`, requiring `]` after `[`. -/
def subscript (part : String) : Option (Nat × Nat) :=
  match part.toList.findIdx? (· = '[') with
  | none => none
  | some open_ =>
      match (part.toList.drop (open_ + 1)).findIdx? (· = ']') with
      | none => none
      | some rel => some (open_, rel + open_ + 1)

/-- The substring between two byte offsets, exclusive. -/
def slice (s : String) (a b : Nat) : String :=
  String.ofList ((s.toList.drop a).take (b - a))

/-- Parse a decimal `Nat`, rejecting anything else. -/
def parseNat (s : String) : Option Nat :=
  let s := trim s
  if s ≠ "" && s.toList.all Char.isDigit then s.toNat? else none

/-! ## Registers -/

/-- A declared register: name, offset into the flat index space, and size. -/
structure Reg where
  /-- The register's name in the source. -/
  name : String
  /-- Where its wires start in the flat numbering. -/
  offset : Nat
  /-- How many wires it has. -/
  size : Nat
deriving Repr, Inhabited

/-- Find a register by name. -/
def findReg (regs : List Reg) (name : String) : Option Reg :=
  regs.find? (·.name == name)

/-- Resolve a comma-separated operand list of subscripted registers to flat indices. -/
def resolveIdx (kind : String) (s : String) (regs : List Reg) (lineNum : Nat) :
    Except String (List Nat) :=
  (splitTrim s ',').foldlM (init := []) fun acc part =>
    let part := trim (dropSemi part)
    match subscript part with
    | none => .ok acc
    | some (bracket, close) =>
        let name := trim (slice part 0 bracket)
        match parseNat (slice part (bracket + 1) close) with
        | none => .error s!"line {lineNum}: bad {kind} index"
        | some i =>
            match findReg regs name with
            | none => .error s!"line {lineNum}: unknown register '{name}'"
            | some r =>
                if i ≥ r.size then
                  .error
                    s!"line {lineNum}: index {i} out of range for register '{name}' (size {r.size})"
                else .ok (acc ++ [r.offset + i])

/-- Rust's `require_arity`. -/
def requireArity (gate : String) (qs : List Nat) (expected lineNum : Nat) :
    Except String Unit :=
  if qs.length = expected then .ok ()
  else
    let noun := if expected == 1 then "operand" else "operands"
    .error s!"line {lineNum}: {gate} expects {expected} qubit {noun}, got {qs.length}"

/-- Resolve exactly one qubit operand. -/
def resolveSingle (gate s : String) (regs : List Reg) (lineNum : Nat) : Except String Nat := do
  let qs ← resolveIdx "qubit" s regs lineNum
  requireArity gate qs 1 lineNum
  match qs with
  | [q] => pure q
  | _ => .error s!"line {lineNum}: {gate} expects 1 qubit operand"

/-- Expand `name[i]` to one index, or a bare `name` to the whole register. -/
def expandOperand (kind : String) (s : String) (regs : List Reg) (lineNum : Nat) :
    Except String (List Nat) :=
  let s := trim (dropSemi s)
  if s.toList.contains '[' then resolveIdx kind s regs lineNum
  else
    match findReg regs s with
    | none => .error s!"line {lineNum}: unknown register '{s}'"
    | some r => .ok ((List.range r.size).map (r.offset + ·))

/-! ## Parsing -/

/-- Parser state threaded through the lines. -/
structure St where
  /-- Quantum registers declared so far. -/
  regs : List Reg := []
  /-- Classical registers declared so far. -/
  cregs : List Reg := []
  /-- Total qubits. -/
  numQubits : Nat := 0
  /-- Total classical bits. -/
  numCbits : Nat := 0
  /-- Gates so far, **most recent first** — `parse` reverses once at the end. Appending with
  `gates ++ [g]` copies the list per gate, which made parsing quadratic: gf2^32's 17,661 gates
  took 3.1 s to parse and gf2^256's 1.1 M were hopeless. -/
  revGates : List Gate := []
  /-- Whether a gate has been seen — declarations may not follow one. -/
  seenGate : Bool := false

/-- Declare a register: `name[size]`. -/
def parseReg (rest : String) (lineNum : Nat) : Except String (Option (String × Nat)) :=
  let rest := trim rest
  match subscript rest with
  | none => .ok none
  | some (bracket, close) =>
      let name := trim (slice rest 0 bracket)
      match parseNat (slice rest (bracket + 1) close) with
      | none => .error s!"line {lineNum}: bad register size"
      | some size => .ok (some (name, size))

/-- Parse a `measure` body into `(qubit, cbit)` pairs, supporting both the indexed form
`measure q[i] -> c[j];` and the register-broadcast form `measure q -> c;`. -/
def parseMeasure (s : String) (st : St) (lineNum : Nat) : Except String (List (Nat × Nat)) :=
  match s.splitOn "->" with
  | [qPart, cPart] => do
      let qs ← expandOperand "qubit" qPart st.regs lineNum
      let cs ← expandOperand "cbit" cPart st.cregs lineNum
      if qs.length ≠ cs.length then
        .error
          s!"line {lineNum}: measure operand size mismatch ({qs.length} qubits, {cs.length} cbits)"
      else .ok (qs.zip cs)
  | _ => .error s!"line {lineNum}: measure missing '->' (got '{trim s}')"

/-- The statement keywords, grouped by first byte exactly as Rust groups them. Order within a
group matters wherever one keyword is a prefix of another. -/
def candidates (first : Char) : List String :=
  match first with
  | 'q' => ["qreg"]
  | 'c' => ["creg", "cx ", "ccz ", "ccx ", "cz "]
  | 'm' => ["measure "]
  | 'r' => ["reset ", "rz("]
  | 'h' => ["h "]
  | 'x' => ["x "]
  | 's' => ["s ", "sdg "]
  | 't' => ["tdg ", "t "]
  | 'z' => ["z "]
  | _ => []

/-- Handle one statement. -/
def parseStatement (line : String) (st : St) (lineNum : Nat) : Except String St := do
  let first := line.front
  -- ignored statements
  if line.startsWith "//" || line.startsWith "OPENQASM" || line.startsWith "include"
      || line.startsWith "barrier" then
    return st
  let some (kw, rest) :=
      (candidates first).findSome? fun kw =>
        if line.startsWith kw then some (kw, (line.drop kw.length).toString) else none
    | .error s!"line {lineNum}: unsupported: {line}"
  let one (mk : Qubit → Gate) (name : String) : Except String St := do
    let q ← resolveSingle name rest st.regs lineNum
    return { st with revGates := mk q :: st.revGates, seenGate := true }
  let many (name : String) (arity : Nat) (mk : List Nat → Gate) : Except String St := do
    let qs ← resolveIdx "qubit" rest st.regs lineNum
    requireArity name qs arity lineNum
    return { st with revGates := mk qs :: st.revGates, seenGate := true }
  match kw with
  | "qreg" =>
      if st.seenGate then .error s!"line {lineNum}: qreg declaration after gate"
      else
        match ← parseReg rest lineNum with
        | none => return st
        | some (name, size) =>
            return { st with regs := st.regs ++ [⟨name, st.numQubits, size⟩],
                             numQubits := st.numQubits + size }
  | "creg" =>
      if st.seenGate then .error s!"line {lineNum}: creg declaration after gate"
      else
        match ← parseReg rest lineNum with
        | none => return st
        | some (name, size) =>
            return { st with cregs := st.cregs ++ [⟨name, st.numCbits, size⟩],
                             numCbits := st.numCbits + size }
  | "measure " =>
      let pairs ← parseMeasure rest st lineNum
      return { st with
                 revGates := (pairs.map fun (q, c) => Gate.measure q c).reverse ++ st.revGates,
                 seenGate := true }
  | "reset " =>
      let qs ← expandOperand "qubit" rest st.regs lineNum
      return { st with revGates := (qs.map Gate.reset).reverse ++ st.revGates,
                       seenGate := true }
  | "cx " => many "cx" 2 fun qs => .cnot (qs.getD 0 0) (qs.getD 1 0)
  | "cz " => many "cz" 2 fun qs => .cz (qs.getD 0 0) (qs.getD 1 0)
  | "ccx " => many "ccx" 3 fun qs => .ccx (qs.getD 0 0) (qs.getD 1 0) (qs.getD 2 0)
  | "ccz " => many "ccz" 3 fun qs => .ccz (qs.getD 0 0) (qs.getD 1 0) (qs.getD 2 0)
  | "h " => one .h "h"
  | "x " => one .x "x"
  | "s " => one .s "s"
  | "sdg " => one .sdg "sdg"
  | "z " => one .z "z"
  | "t " => one .t "t"
  | "tdg " => one .tdg "tdg"
  | "rz(" =>
      .error
        s!"line {lineNum}: rz is not supported by this build — angles here are exact rationals \
           in units of π, and Rz synthesis (gridsynth) is not ported"
  | _ => .error s!"line {lineNum}: unsupported: {line}"

/-- Read the statements of an OpenQASM 2.0 source into a circuit.

This is the parser proper; `parse` is this followed by `validate`, and everything outside
this module should call `parse`. -/
def parseRaw (source : String) : Except String Circuit := do
  let source := stripBlockComments source
  let mut st : St := {}
  let mut lineNum := 0
  for rawLine in source.splitOn "\n" do
    lineNum := lineNum + 1
    for piece in splitTrim rawLine ';' do
      let line :=
        match lineComment piece with
        | some pos => trim (slice piece 0 pos)
        | none => piece
      if line ≠ "" then
        st ← parseStatement line st lineNum
  return Circuit.ofGates st.numQubits st.numCbits st.revGates.reverse

/-! ## Validating

The parser checks each operand against its register's size, but nothing checks that a
multi-qubit gate's operands are *distinct*: `cx q[0],q[0]` parses. That circuit cannot be
packaged as the `Circuit.Checked` input a pass requires. The restriction is not idle —
`cnot q q` denotes `b ↦ b[q := 0]`, which is idempotent rather than
self-inverse, so `CancelGates` deleting a pair of them is a genuinely unsound rewrite. (The
Rust front end has the same gap; there it is a latent bug rather than a broken proof.)

So the front end rejects such a circuit rather than passing it on, and `parse_wf` below
turns that check into the proof used to construct a checked circuit. The range check comes
along for the ride: it is implied by the per-register bound in `resolveIdx`, but checking it
here is cheaper than threading a proof through the parser's state monad, and it is what
`parse_wellFormed` needs so the back end can only ever emit subscripts that are in range.
-/

/-- Name the gate that failed validation, for the error message. -/
def badGate (c : Circuit) : String :=
  match c.gates.find? (fun g => !decide g.Wf) with
  | some g =>
      s!"gate '{g}' repeats a qubit operand — a multi-qubit gate needs distinct qubits"
  | none =>
      match c.gates.find? (fun g => !decide (g.InRange c.numQubits c.numCbits)) with
      | some g => s!"gate '{g}' addresses a register slot that was never declared"
      | none => "circuit failed validation"

/-- Accept a parsed circuit only if every gate has distinct operands and is in range,
rebuilding it so that its cached `has*` flags are honest too. -/
def validate (c : Circuit) : Except String Circuit :=
  if ∀ g ∈ c.gates, g.Wf ∧ g.InRange c.numQubits c.numCbits then
    .ok (c.withGates c.gates)
  else
    .error (badGate c)

/-- Validation returns the circuit it was given, up to the rebuilt flags. -/
theorem validate_eq {c c' : Circuit} (h : validate c = .ok c') : c' = c.withGates c.gates := by
  unfold validate at h
  split at h
  · exact (Except.ok.injEq _ _ ▸ h).symm ▸ rfl
  · exact absurd h (by simp)

/-- **A validated circuit has distinct multi-qubit operands** — the invariant carried by
`Circuit.Checked`. -/
theorem validate_wf {c c' : Circuit} (h : validate c = .ok c') : c'.Wf := by
  unfold validate at h
  split at h
  · rename_i hall
    obtain rfl : c' = c.withGates c.gates := (Except.ok.injEq _ _ ▸ h).symm ▸ rfl
    exact fun g hg => (hall g hg).1
  · exact absurd h (by simp)

/-- **A validated circuit is well-formed**: every operand is a slot that was declared. -/
theorem validate_wellFormed {c c' : Circuit} (h : validate c = .ok c') : c'.WellFormed := by
  unfold validate at h
  split at h
  · rename_i hall
    obtain rfl : c' = c.withGates c.gates := (Except.ok.injEq _ _ ▸ h).symm ▸ rfl
    exact fun g hg => (hall g hg).2
  · exact absurd h (by simp)

/-- **A validated circuit's flags are honest**, since validation rebuilds them. -/
theorem validate_flagsOk {c c' : Circuit} (h : validate c = .ok c') : c'.FlagsOk := by
  rw [validate_eq h]; exact Circuit.flagsOk_withGates _ _

/-- Parse a circuit from OpenQASM 2.0 source, rejecting anything the pipeline may not
assume. -/
def parse (source : String) : Except String Circuit := do
  let c ← parseRaw source
  validate c

theorem parse_eq_validate {source : String} {c : Circuit} (h : parseRaw source = .ok c) :
    parse source = validate c := by
  simp only [parse, h, bind, Except.bind]

/-- **Everything the parser promises the rest of the compiler.** `Circuit.Wf` permits entry
into the checked optimizer; `WellFormed` and `FlagsOk` support checked serialization. -/
theorem parse_valid {source : String} {c : Circuit} (h : parse source = .ok c) :
    c.Wf ∧ c.WellFormed ∧ c.FlagsOk := by
  unfold parse at h
  simp only [bind, Except.bind] at h
  cases hraw : parseRaw source with
  | error e => rw [hraw] at h; exact absurd h (by simp)
  | ok c₀ =>
      rw [hraw] at h
      exact ⟨validate_wf h, validate_wellFormed h, validate_flagsOk h⟩

theorem parse_wf {source : String} {c : Circuit} (h : parse source = .ok c) : c.Wf :=
  (parse_valid h).1

theorem parse_wellFormed {source : String} {c : Circuit} (h : parse source = .ok c) :
    c.WellFormed := (parse_valid h).2.1

/-! ## Serializing -/

/-- One gate, in OpenQASM syntax. -/
def gateLine : Gate → String
  | .x q => s!"x q[{q}];"
  | .h q => s!"h q[{q}];"
  | .s q => s!"s q[{q}];"
  | .sdg q => s!"sdg q[{q}];"
  | .z q => s!"z q[{q}];"
  | .t q => s!"t q[{q}];"
  | .tdg q => s!"tdg q[{q}];"
  | .rz θ q => s!"rz({θ}*pi) q[{q}];"
  | .cnot c t => s!"cx q[{c}],q[{t}];"
  | .cz c t => s!"cz q[{c}],q[{t}];"
  | .ccx a b t => s!"ccx q[{a}],q[{b}],q[{t}];"
  | .ccz a b t => s!"ccz q[{a}],q[{b}],q[{t}];"
  | .measure q c => s!"measure q[{q}] -> c[{c}];"
  | .reset q => s!"reset q[{q}];"

/-- Serialize a circuit to OpenQASM 2.0. -/
def serialize (c : Circuit) : String :=
  let header :=
    "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n" ++ s!"qreg q[{c.numQubits}];\n" ++
      (if c.numCbits > 0 then s!"creg c[{c.numCbits}];\n" else "")
  c.gates.foldl (fun acc g => acc ++ gateLine g ++ "\n") header

/-- Serialize only when the parser confirms that the emitted text reconstructs the same
circuit.  This executable check makes the output boundary fail closed, including for syntax
such as `rz` that this build can print but intentionally cannot parse. -/
def serializeChecked (c : Circuit) : Except String String :=
  let source := serialize c
  match parse source with
  | .error e => .error s!"internal QASM serialization check failed: {e}"
  | .ok c' =>
      if c' = c.withGates c.gates then .ok source
      else .error "internal QASM serialization check failed: round trip changed the circuit"

/-- A successful checked serialization parses to exactly the circuit supplied, modulo the
intentional rebuilding of cached flags. -/
theorem serializeChecked_sound {c : Circuit} {source : String}
    (h : serializeChecked c = .ok source) :
    source = serialize c ∧ parse source = .ok (c.withGates c.gates) := by
  cases hp : parse (serialize c) with
  | error e => simp [serializeChecked, hp] at h
  | ok c' =>
      by_cases heq : c' = c.withGates c.gates
      · have hs : source = serialize c := by
          have : serialize c = source := by simpa [serializeChecked, hp, heq] using h
          exact this.symm
        subst source
        exact ⟨rfl, by simpa [heq] using hp⟩
      · simp [serializeChecked, hp, heq] at h

end Qasm

end TzapLean
