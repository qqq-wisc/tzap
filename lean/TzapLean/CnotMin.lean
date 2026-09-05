import TzapLean.Pass

/-!
# The `CnotMin` Pass

A port of `src/cnot_min.rs`: CNOT minimization by phase-polynomial resynthesis (Amy,
Azimzadeh and Mosca; Feynman's `-cnotmin`).

Between the gates it cannot interpret, a circuit is *CNOT-dihedral*: every gate either
permutes computational-basis states (`cnot`, `x`) or applies a diagonal phase (`z`, `s`,
`t`, `rz`, `cz`). Such a block is determined by two things and nothing else:

* the **linear map** it applies — each qubit ends up holding some XOR of the values the
  block's qubits held on entry, plus a constant; and
* the **phase polynomial** — a set of (parity, angle) pairs.

So a block can be discarded and re-synthesized from that pair. This file has three parts:

1. **The analysis** (`BlockState`, `feedGate`) — the forward scan of `Chunk::feed`, which
   recovers the pair. Its soundness against the density-matrix semantics is proved in
   `TzapLean/CnotMinProof.lean`.
2. **The synthesis** (`graySynth`, `linearSynth`, `synthesize`) — Gray-code phase synthesis
   followed by a Gaussian-elimination fix-up of the linear map, ported directly.
3. **The pass** — chunking, and the accept/reject decision.

## Why the resynthesis is checked rather than proved

`graySynth` is a heuristic search: it picks the column that splits the parity set most
evenly, and its correctness argument is a nontrivial invariant about Gray-code ordering.
Rather than prove that, the pass **re-analyses its own output**: a synthesized block
replaces the original only when running the (proved-sound) analysis over it returns the
*same* linear map and phase polynomial. Correctness therefore rests on the analysis alone,
and a bug in the synthesis heuristic could only cost optimization, never soundness. The
check is linear in the block size and, on every test ported from Rust, always passes — so
the pass's output matches the Rust pass gate for gate.

The proof-facing analysis and synthesis below keep the readable `Nat` parity specification.
Production chunking uses a certifying fast path: two `UInt64`s represent the same 128-bit
parity, phase terms live in a hash map, Gray-code worklists are arrays, and a Rust-style
`SynthBudget` abandons a candidate as soon as its gate or two-qubit count cannot win.  None of
those representations enter the correctness theorem: `Chunk.flush` re-analyses both the
original and the candidate with the proved `BlockState` analysis before accepting it.
-/

namespace TzapLean

/-- A parity over a block's starting values, as a bitmask over local qubit indices. -/
abbrev Parity := Nat

/-- The state of one CNOT-dihedral block: what each local qubit now holds, and the phase
polynomial accumulated so far. -/
structure BlockState where
  /-- Local index → the parity of starting values that qubit now holds. -/
  parity : List Parity
  /-- Local index → whether an odd number of `x` gates has been applied on top. -/
  consts : List Bool
  /-- Rotation angle (in units of `π`) per parity. -/
  terms : List (Parity × ℚ)
  deriving Repr, DecidableEq, Inhabited

namespace BlockState

/-- The block that has done nothing yet to `n` qubits: qubit `i` holds its own value. -/
def initial (n : Nat) : BlockState where
  parity := (List.range n).map (2 ^ ·)
  consts := List.replicate n false
  terms := []

/-- Add `s` to the angle recorded for parity `p`, appending the parity if it is new. -/
def insertTerm : List (Parity × ℚ) → Parity → ℚ → List (Parity × ℚ)
  | [], p, s => [(p, s)]
  | (q, a) :: ts, p, s => if q == p then (q, a + s) :: ts else (q, a) :: insertTerm ts p s

/-- Fold a rotation of `angle` on parity `p` (complemented when `k`) into the phase map.
A rotation on the complement of `p` is one on `p` with the opposite angle, plus a global
phase — which is dropped, as everything here compares circuits up to global phase. -/
def addTerm (st : BlockState) (p : Parity) (k : Bool) (angle : ℚ) : BlockState :=
  { st with terms := insertTerm st.terms p (if k then -angle else angle) }

/-- An angle that is a multiple of `2π`, i.e. no rotation at all. -/
def angleIsZero (a : ℚ) : Bool := a - 2 * ⌊a / 2⌋ == 0

/-- `a` reduced to `[0, 2)`. -/
def angleMod (a : ℚ) : ℚ := a - 2 * ⌊a / 2⌋

/-- Canonical form for comparing two blocks: drop zero rotations, reduce angles mod `2π`,
and sort by parity. -/
def normalize (st : BlockState) : BlockState where
  parity := st.parity
  consts := st.consts
  terms :=
    ((st.terms.filter fun t => !angleIsZero t.2).map fun t => (t.1, angleMod t.2)).mergeSort
      fun a b => a.1 ≤ b.1

end BlockState

/-! ## Fast, untrusted construction state

`BlockState` is deliberately list-backed because its small recursive definitions make the
semantic proof direct.  It is also the state used by `analyzeGates`, the certificate checker
at the boundary of this pass.

The live chunk does not need to carry that proof representation: `Chunk.flush` re-analyzes
both the original and proposed gate lists with `analyzeGates` before accepting a rewrite.
Keep its phase polynomial in a hash map instead, matching Rust's `FxHashMap` organization.
A bug here can only propose a bad synthesis, which the proved checker rejects. -/

/-- A native 128-bit parity for the untrusted search path.  Lean's proof-facing `Parity` is
an unbounded `Nat`; using it in synthesis routes every bit test through GMP. -/
structure FastParity where
  lo : UInt64 := 0
  hi : UInt64 := 0
deriving BEq, Hashable, Inhabited

namespace FastParity

def basis (i : Nat) : FastParity :=
  if i < 64 then ⟨1 <<< i.toUInt64, 0⟩ else ⟨0, 1 <<< (i - 64).toUInt64⟩

def xor (a b : FastParity) : FastParity := ⟨a.lo ^^^ b.lo, a.hi ^^^ b.hi⟩

def testBit (p : FastParity) (i : Nat) : Bool :=
  if i < 64 then ((p.lo >>> i.toUInt64) &&& 1) == 1
  else ((p.hi >>> (i - 64).toUInt64) &&& 1) == 1

/-- Numeric order, matching the `Nat` parity order used by the reference synthesis. -/
def le (a b : FastParity) : Bool :=
  a.hi < b.hi || (a.hi == b.hi && a.lo ≤ b.lo)

end FastParity

/-- Search state for a live chunk.  This is converted only through synthesized gates; it is
never consumed by a correctness theorem. -/
structure FastBlockState where
  parity : Array FastParity
  consts : Array Bool
  terms : Std.HashMap FastParity ℚ
deriving Inhabited

namespace FastBlockState

def initial (n : Nat) : FastBlockState where
  parity := (Array.range n).map FastParity.basis
  consts := Array.replicate n false
  terms := ∅

/-- Add one rotation in expected constant time. -/
def addTerm (st : FastBlockState) (p : FastParity) (k : Bool) (angle : ℚ) : FastBlockState :=
  let delta := if k then -angle else angle
  let total := (st.terms.get? p).getD 0 + delta
  { st with terms := st.terms.insert p total }

end FastBlockState

/-- The rotation a diagonal single-qubit gate applies, in units of `π`. -/
def rotAngle : Gate → Option (ℚ × Qubit)
  | .t q => some (1/4, q)
  | .tdg q => some (-1/4, q)
  | .s q => some (1/2, q)
  | .sdg q => some (-1/2, q)
  | .z q => some (1, q)
  | .rz θ q => some (θ, q)
  | _ => none

/-- The wire `q`'s local index in the block, if it has one. -/
def localIdx (qs : List Qubit) (q : Qubit) : Option Nat := qs.findIdx? (· == q)

/-- Feed one gate into a block over the fixed local qubit list `qs`. `none` when the gate
leaves the CNOT-dihedral fragment or touches a qubit the block does not cover. -/
def feedGate (qs : List Qubit) (st : BlockState) (g : Gate) : Option BlockState :=
  let idx : Qubit → Option Nat := localIdx qs
  match rotAngle g with
  | some (θ, q) =>
      match idx q with
      | some i =>
          match st.parity[i]?, st.consts[i]? with
          | some p, some k => some (st.addTerm p k θ)
          | _, _ => none
      | none => none
  | none =>
      match g with
      | .x q =>
          match idx q with
          | some i => some { st with consts := st.consts.set i (!st.consts[i]!) }
          | none => none
      | .cnot c t =>
          if c == t then none else
          match idx c, idx t with
          | some ci, some ti =>
              some { st with
                parity := st.parity.set ti (st.parity[ti]! ^^^ st.parity[ci]!)
                consts := st.consts.set ti (st.consts[ti]! != st.consts[ci]!) }
          | _, _ => none
      | .cz c t =>
          if c == t then none else
          match idx c, idx t with
          | some ci, some ti =>
              let pc := st.parity[ci]!; let pt := st.parity[ti]!
              let kc := st.consts[ci]!; let kt := st.consts[ti]!
              some (((st.addTerm pc kc (1/2)).addTerm pt kt (1/2)).addTerm
                (pc ^^^ pt) (kc != kt) (-1/2))
          | _, _ => none
      | _ => none

/-- Untrusted, hash-backed counterpart of `feedGate` used while constructing a chunk.
`Chunk.flush` validates its resulting synthesis with the proved `analyzeGates`. -/
def feedGateFast (qs : List Qubit) (st : FastBlockState) (g : Gate) : Option FastBlockState :=
  let idx : Qubit → Option Nat := localIdx qs
  match rotAngle g with
  | some (θ, q) =>
      match idx q with
      | some i =>
          match st.parity[i]?, st.consts[i]? with
          | some p, some k => some (st.addTerm p k θ)
          | _, _ => none
      | none => none
  | none =>
      match g with
      | .x q =>
          match idx q with
          | some i => some { st with consts := st.consts.set! i (!st.consts[i]!) }
          | none => none
      | .cnot c t =>
          if c == t then none else
          match idx c, idx t with
          | some ci, some ti =>
              some { st with
                parity := st.parity.set! ti (st.parity[ti]!.xor st.parity[ci]!)
                consts := st.consts.set! ti (st.consts[ti]! != st.consts[ci]!) }
          | _, _ => none
      | .cz c t =>
          if c == t then none else
          match idx c, idx t with
          | some ci, some ti =>
              let pc := st.parity[ci]!; let pt := st.parity[ti]!
              let kc := st.consts[ci]!; let kt := st.consts[ti]!
              some (((st.addTerm pc kc (1/2)).addTerm pt kt (1/2)).addTerm
                (pc.xor pt) (kc != kt) (-1/2))
          | _, _ => none
      | _ => none

/-- Analyse a gate list as one block over the local qubits `qs`. -/
def analyzeGates (qs : List Qubit) (gs : List Gate) : Option BlockState :=
  gs.foldl (fun st g => st.bind fun st => feedGate qs st g) (some (BlockState.initial qs.length))

/-! ## Synthesis -/

/-- Emit `angle · π` on `qubit` as Clifford+T gates when it is a multiple of `π/4`, and as
an `rz` otherwise — the Rust `emit_rotation`. -/
def emitRotation (q : Qubit) (a : ℚ) : List Gate :=
  let n := BlockState.angleMod a
  if n == 0 then []
  else
    let k := 4 * n
    if k.den == 1 then
      match ((k.num % 8 + 8) % 8).toNat with
      | 0 => []
      | 1 => [.t q]
      | 2 => [.s q]
      | 3 => [.s q, .t q]
      | 4 => [.z q]
      | 5 => [.z q, .t q]
      | 6 => [.sdg q]
      | _ => [.tdg q]
    else [.rz n q]

/-- One node of the Gray-code recursion. -/
structure Pt where
  candidates : List Nat
  target : Option Nat
  pending : Option Nat
  vectors : List (Parity × ℚ)

/-- Pick the column whose 0/1 split of `vectors` leaves the largest side — the Rust
`best_split`.

The fold carries the incumbent's score rather than recomputing it. `score` is a full pass
over `vectors`, and the original spelling evaluated it on both sides of the comparison, so
the incumbent's was recomputed once per candidate. -/
def bestSplit (candidates : List Nat) (vectors : List (Parity × ℚ)) (fallback : Nat) :
    List (Parity × ℚ) × Nat × List Nat × List (Parity × ℚ) :=
  let score := fun c =>
    let ones := (vectors.filter fun v => v.1.testBit c).length
    max ones (vectors.length - ones)
  let best :=
    (candidates.foldl (fun (acc, sa) c =>
      let sc := score c
      if sc ≥ sa then (c, sc) else (acc, sa)) (fallback, score fallback)).1
  let ones := vectors.filter fun v => v.1.testBit best
  let zeros := vectors.filter fun v => !v.1.testBit best
  (zeros, best, candidates.filter (· != best), ones)

/-- Gray-code phase synthesis: emit CNOTs walking the qubits through every parity that
carries a rotation, applying each rotation as its parity comes up. Returns the gates and the
linear map the emitted CNOTs left behind.

`state` and `qs` are arrays: this loop reads and updates them once per emitted gate, and on a
128-wire block that is the pass's inner loop. `acc` is built in reverse for the same reason —
`acc ++ [gate]` copies the whole accumulator, so emitting `g` gates cost `O(g²)`. `graySynth`
puts the order back. None of this is proof-carrying (see the header), so the representation
is free to be whatever runs fastest. -/
def graySynthLoop : Nat → Array Qubit → List Pt → Array Parity → List Gate →
    List Gate × Array Parity
  | 0, _, _, state, acc => (acc, state)
  | _ + 1, _, [], state, acc => (acc, state)
  | fuel + 1, qs, node :: stack, state, acc =>
      if node.vectors.isEmpty then graySynthLoop fuel qs stack state acc
      else
        match node.target, node.pending with
        | some t, some p =>
            let gate := Gate.cnot qs[p]! qs[t]!
            let state' := state.set! t (state[t]! ^^^ state[p]!)
            let stack' := stack.map fun other =>
              { other with vectors := other.vectors.map fun v =>
                  if v.1.testBit t then (v.1 ^^^ 2 ^ p, v.2) else v }
            graySynthLoop fuel qs ({ node with pending := none } :: stack') state' (gate :: acc)
        | _, _ =>
            match node.candidates with
            | [] =>
                match node.target, node.vectors with
                | some t, [(_, a)] =>
                    graySynthLoop fuel qs stack state ((emitRotation qs[t]! a).reverse ++ acc)
                | _, _ => graySynthLoop fuel qs stack state acc
            | first :: _ =>
                let (zeros, col, rest, ones) := bestSplit node.candidates node.vectors first
                let zeroNode : Pt :=
                  { candidates := rest, target := node.target, pending := none, vectors := zeros }
                let oneNode : Pt :=
                  match node.target with
                  | some t => { candidates := rest, target := some t, pending := some col,
                                vectors := ones }
                  | none => { candidates := rest, target := some col, pending := none,
                              vectors := ones }
                graySynthLoop fuel qs (zeroNode :: oneNode :: stack) state acc

/-- `graySynth` at a fuel bound that always suffices: the recursion tree has at most one
leaf per phase and depth at most `n`, and each node costs a constant number of steps. -/
def graySynth (n : Nat) (phases : List (Parity × ℚ)) (state : Array Parity) (qs : Array Qubit) :
    List Gate × Array Parity :=
  if phases.isEmpty then ([], state)
  else
    let (acc, state') :=
      graySynthLoop (4 * (phases.length + 1) * (n + 1) + 8) qs
        [{ candidates := List.range n, target := none, pending := none, vectors := phases }]
        state []
    (acc.reverse, state')

/-! ## Bounded synthesis

The unbounded functions above are useful executable specifications and direct test targets.
Production chunking uses the budget below.  Gate count and two-qubit count only increase as
synthesis proceeds, so exceeding either count of the original block proves that no
continuation can pass `Chunk.accepts`.  Returning `none` at that point is therefore exactly
the same optimizer decision as completing and rejecting the larger circuit. -/

structure SynthBudget where
  rev : List Gate := []
  count : Nat := 0
  twoQ : Nat := 0
  maxCount : Nat
  maxTwoQ : Nat

namespace SynthBudget

def push (b : SynthBudget) (g : Gate) : Option SynthBudget :=
  let count := b.count + 1
  let twoQ := b.twoQ + match g with | .cnot .. | .cz .. => 1 | _ => 0
  if count > b.maxCount || twoQ > b.maxTwoQ then none
  else some { b with rev := g :: b.rev, count, twoQ }

def pushAll (b : SynthBudget) (gs : List Gate) : Option SynthBudget :=
  gs.foldl (fun acc g => acc.bind (·.push g)) (some b)

def admitsTwoQ (b : SynthBudget) (extra : Nat) : Bool :=
  b.twoQ + extra ≤ b.maxTwoQ && b.count + extra ≤ b.maxCount

end SynthBudget

/-- Gray-code node for the native parity representation. -/
structure FastPt where
  candidates : Array Nat
  target : Option Nat
  pending : Option Nat
  vectors : Array (FastParity × ℚ)

def bestSplitFast (candidates : Array Nat) (vectors : Array (FastParity × ℚ))
    (fallback : Nat) : Array (FastParity × ℚ) × Nat × Array Nat × Array (FastParity × ℚ) :=
  Id.run do
    let mut best := fallback
    let mut bestScore := 0
    for c in candidates do
      let mut ones := 0
      for v in vectors do
        if v.1.testBit c then ones := ones + 1
      let score := max ones (vectors.size - ones)
      if score ≥ bestScore then
        best := c
        bestScore := score
    let mut zeros := #[]
    let mut ones := #[]
    for v in vectors do
      if v.1.testBit best then ones := ones.push v else zeros := zeros.push v
    let rest := candidates.filter (· != best)
    return (zeros, best, rest, ones)

/-- Budgeted native-parity Gray synthesis; `none` means the candidate has already lost. -/
def graySynthLoopBudget : Nat → Array Qubit → List FastPt → Array FastParity → SynthBudget →
    Option (SynthBudget × Array FastParity)
  | 0, _, _, state, budget => some (budget, state)
  | _ + 1, _, [], state, budget => some (budget, state)
  | fuel + 1, qs, node :: stack, state, budget =>
      if node.vectors.isEmpty then graySynthLoopBudget fuel qs stack state budget
      else
        match node.target, node.pending with
        | some t, some p =>
            match budget.push (Gate.cnot qs[p]! qs[t]!) with
            | none => none
            | some budget' =>
                let state' := state.set! t (state[t]!.xor state[p]!)
                let stack' := stack.map fun other =>
                  { other with vectors := other.vectors.map fun v =>
                      if v.1.testBit t then (v.1.xor (FastParity.basis p), v.2) else v }
                graySynthLoopBudget fuel qs ({ node with pending := none } :: stack') state' budget'
        | _, _ =>
            if node.candidates.isEmpty then
                match node.target with
                | some t =>
                    if node.vectors.size == 1 then
                      let a := node.vectors[0]!.2
                    match budget.pushAll (emitRotation qs[t]! a) with
                    | none => none
                    | some budget' => graySynthLoopBudget fuel qs stack state budget'
                    else graySynthLoopBudget fuel qs stack state budget
                | none => graySynthLoopBudget fuel qs stack state budget
            else
                let first := node.candidates[0]!
                let (zeros, col, rest, ones) := bestSplitFast node.candidates node.vectors first
                let zeroNode : FastPt :=
                  { candidates := rest, target := node.target, pending := none, vectors := zeros }
                let oneNode : FastPt :=
                  match node.target with
                  | some t => { candidates := rest, target := some t, pending := some col,
                                vectors := ones }
                  | none => { candidates := rest, target := some col, pending := none,
                              vectors := ones }
                graySynthLoopBudget fuel qs (zeroNode :: oneNode :: stack) state budget

def graySynthBudget (n : Nat) (phases : List (FastParity × ℚ)) (state : Array FastParity)
    (qs : Array Qubit) (budget : SynthBudget) : Option (SynthBudget × Array FastParity) :=
  if phases.isEmpty then some (budget, state)
  else
    graySynthLoopBudget (4 * (phases.length + 1) * (n + 1) + 8) qs
      [{ candidates := Array.range n, target := none, pending := none, vectors := phases.toArray }]
      state budget

/-- XOR together the rows of `m` selected by the set bits of `row`. -/
def rowTimesMatrix (row : Parity) (m : Array Parity) : Parity := Id.run do
  let mut acc : Parity := 0
  let mut i := 0
  for r in m do
    if row.testBit i then acc := acc ^^^ r
    i := i + 1
  return acc

/-- Invert an `n × n` bit matrix by Gauss-Jordan, or `none` if singular.

Gauss-Jordan is `O(n³)` row operations, and on a `List` every one of those pays another
`O(n)` to reach the row and `O(n)` to copy the spine — `O(n⁴)` overall, which at the 128-wire
chunk cap is the difference between a pass that finishes and one that does not. `Array.set!`
updates in place while the reference is unique, so the elimination costs what it says. -/
def invert (rows : Array Parity) (n : Nat) : Option (Array Parity) := Id.run do
  let mut a := rows
  let mut inv : Array Parity := (Array.range n).map (2 ^ ·)
  for col in List.range n do
    match (List.range n).find? (fun r => r ≥ col && a[r]!.testBit col) with
    | none => return none
    | some pivot =>
        let ac := a[col]!; let ap := a[pivot]!
        a := (a.set! col ap).set! pivot ac
        let ic := inv[col]!; let ip := inv[pivot]!
        inv := (inv.set! col ip).set! pivot ic
        for r in List.range n do
          if r != col && a[r]!.testBit col then
            a := a.set! r (a[r]! ^^^ a[col]!)
            inv := inv.set! r (inv[r]! ^^^ inv[col]!)
  return some inv

/-- Emit CNOTs taking the qubits from parities `frm` to parities `to`, by Gauss-Jordan on
`to · frm⁻¹` — the Rust `linear_synth`.

`ops` is consed rather than appended: the loop emits up to `n²/2` of them, and building that
with `ops ++ [op]` copies the list every time. Consing produces exactly the order the old
trailing `.reverse` did, so the emitted gates are unchanged. -/
def linearSynth (frm tgt : Array Parity) (qs : Array Qubit) : Option (List Gate) := Id.run do
  let n := frm.size
  if frm == tgt then return some []
  match invert frm n with
  | none => return none
  | some inverse =>
      let mut m : Array Parity := tgt.map (fun row => rowTimesMatrix row inverse)
      let mut ops : List (Nat × Nat) := []
      for col in List.range n do
        if !(m[col]!.testBit col) then
          match (List.range n).find? (fun r => r > col && m[r]!.testBit col) with
          | none => return none
          | some r =>
              m := m.set! col (m[col]! ^^^ m[r]!)
              ops := (r, col) :: ops
        for r in List.range n do
          if r != col && m[r]!.testBit col then
            m := m.set! r (m[r]! ^^^ m[col]!)
            ops := (col, r) :: ops
      return some (ops.map fun (c, t) => Gate.cnot qs[c]! qs[t]!)

def rowTimesMatrixFast (row : FastParity) (m : Array FastParity) : FastParity := Id.run do
  let mut acc : FastParity := {}
  let mut i := 0
  for r in m do
    if row.testBit i then acc := acc.xor r
    i := i + 1
  return acc

def invertFast (rows : Array FastParity) (n : Nat) : Option (Array FastParity) := Id.run do
  let mut a := rows
  let mut inv : Array FastParity := (Array.range n).map FastParity.basis
  for col in List.range n do
    match (List.range n).find? (fun r => r ≥ col && a[r]!.testBit col) with
    | none => return none
    | some pivot =>
        let ac := a[col]!; let ap := a[pivot]!
        a := (a.set! col ap).set! pivot ac
        let ic := inv[col]!; let ip := inv[pivot]!
        inv := (inv.set! col ip).set! pivot ic
        for r in List.range n do
          if r != col && a[r]!.testBit col then
            a := a.set! r (a[r]!.xor a[col]!)
            inv := inv.set! r (inv[r]!.xor inv[col]!)
  return some inv

/-- Budgeted linear fix-up.  Each pending row operation becomes exactly one CNOT, so the
budget can reject a losing block before elimination finishes. -/
def linearSynthBudget (frm tgt : Array FastParity) (qs : Array Qubit)
    (budget : SynthBudget) : Option SynthBudget := Id.run do
  let n := frm.size
  if frm == tgt then return some budget
  match invertFast frm n with
  | none => return none
  | some inverse =>
      let mut m : Array FastParity := tgt.map (fun row => rowTimesMatrixFast row inverse)
      let mut ops : List (Nat × Nat) := []
      let mut count := 0
      for col in List.range n do
        if !(m[col]!.testBit col) then
          match (List.range n).find? (fun r => r > col && m[r]!.testBit col) with
          | none => return none
          | some r =>
              m := m.set! col (m[col]!.xor m[r]!)
              ops := (r, col) :: ops
              count := count + 1
              if !budget.admitsTwoQ count then return none
        for r in List.range n do
          if r != col && m[r]!.testBit col then
            m := m.set! r (m[r]!.xor m[col]!)
            ops := (col, r) :: ops
            count := count + 1
            if !budget.admitsTwoQ count then return none
      return ops.foldl (fun acc (c, t) => acc.bind (·.push (Gate.cnot qs[c]! qs[t]!)))
        (some budget)

/-- Synthesize a block: Gray-code phase synthesis, the linear fix-up, then the `x` gates the
affine part calls for. The block's wire list and linear map cross into array form here, once,
and stay there for the rest of the synthesis. -/
def synthesize (qs : List Qubit) (st : BlockState) : Option (List Gate) :=
  let n := qs.length
  let qa := qs.toArray
  let phases :=
    ((st.terms.filter fun t => !BlockState.angleIsZero t.2).map fun t =>
      (t.1, BlockState.angleMod t.2)).mergeSort fun a b => a.1 ≤ b.1
  let state : Array Parity := (Array.range n).map (2 ^ ·)
  let (gs₁, state') := graySynth n phases state qa
  match linearSynth state' st.parity.toArray qa with
  | none => none
  | some gs₂ =>
      let xs := (st.consts.zipIdx).filterMap fun (k, i) => if k = true then some (Gate.x qa[i]!) else none
      some (gs₁ ++ gs₂ ++ xs)

/-- Synthesize from the live hash-backed state.  Sorting by parity makes the result identical
to the list-backed route regardless of hash-table iteration order. -/
def synthesizeFast (qs : List Qubit) (st : FastBlockState) (maxTwoQ maxCount : Nat) :
    Option (List Gate) :=
  let n := qs.length
  let qa := qs.toArray
  let phases :=
    ((st.terms.toList.filter fun t => !BlockState.angleIsZero t.2).map fun t =>
      (t.1, BlockState.angleMod t.2)).mergeSort fun a b => FastParity.le a.1 b.1
  let state : Array FastParity := (Array.range n).map FastParity.basis
  let budget : SynthBudget := { maxCount, maxTwoQ }
  match graySynthBudget n phases state qa budget with
  | none => none
  | some (budget, state') =>
      match linearSynthBudget state' st.parity qa budget with
      | none => none
      | some budget =>
          let xs := (st.consts.toList.zipIdx).filterMap fun (k, i) =>
            if k = true then some (Gate.x qa[i]!) else none
          match budget.pushAll xs with
          | none => none
          | some budget => some budget.rev.reverse


/-! ## Chunking -/

/-- One CNOT-dihedral block under construction, with the gates it was built from. -/
structure Chunk where
  /-- The circuit's qubit count: an operand at or above it is not interpretable, since the
  block's model of a wire is a wire the circuit actually has. -/
  numQubits : Nat
  maxQubits : Nat
  maxTerms : Nat
  qubits : List Qubit := []
  state : FastBlockState := FastBlockState.initial 0
  original : List Gate := []

namespace Chunk

def empty (numQubits maxQubits maxTerms : Nat) : Chunk := { numQubits, maxQubits, maxTerms }

def reset (ch : Chunk) : Chunk := empty ch.numQubits ch.maxQubits ch.maxTerms

/-- Give `q` a local index, appending it to the block holding its own starting value. -/
def extend (ch : Chunk) (q : Qubit) : Chunk :=
  if ch.qubits.contains q then ch
  else
    { ch with
      qubits := ch.qubits ++ [q]
      state :=
        { ch.state with
          parity := ch.state.parity.push (FastParity.basis ch.qubits.length)
          consts := ch.state.consts.push false } }

/-- The operands a gate would bring into the block, or `none` when it leaves the fragment. -/
def interpretable : Gate → Option (List Qubit)
  | .t q | .tdg q | .s q | .sdg q | .z q | .rz _ q | .x q => some [q]
  | .cnot c t | .cz c t => if c == t then none else some [c, t]
  | _ => none

/-- Whether `ops` fit under the qubit cap without evicting anything. -/
def admits (ch : Chunk) (ops : List Qubit) : Bool :=
  let fresh := (ops.eraseDups.filter fun q => !ch.qubits.contains q).length
  ch.qubits.length + fresh ≤ ch.maxQubits

/-- Number of two-qubit gates. -/
def count2q (gs : List Gate) : Nat :=
  gs.countP fun g => match g with | .cnot .. | .cz .. => true | _ => false

/-- The Rust acceptance test: no more gates overall, and lexicographically fewer
`(two-qubit, total)`. -/
def accepts (orig synth : List Gate) : Bool :=
  (synth.length ≤ orig.length) &&
    ((count2q synth < count2q orig) ||
      ((count2q synth == count2q orig) && (synth.length < orig.length)))

/-- Emit the block: the synthesized form when it is strictly better *and* re-analyses to the
same linear map and phase polynomial as the block it replaces, otherwise the gates the block
was built from.

Both sides of that comparison are re-derived with `analyzeGates` rather than taken from the
running state, so the certificate stands on its own: `flush_correct` needs nothing about how
the chunk was built. -/
def flush (ch : Chunk) : List Gate :=
  if ch.original.isEmpty then []
  else
    match synthesizeFast ch.qubits ch.state (count2q ch.original) ch.original.length with
    | none => ch.original
    | some synth =>
        match analyzeGates ch.qubits ch.original, analyzeGates ch.qubits synth with
        | some d₁, some d₂ =>
            if Chunk.accepts ch.original synth && (d₁.normalize == d₂.normalize) then synth
            else ch.original
        | _, _ => ch.original

/-- Whether a rotation's parity already has a term, or there is room for a new one — the
Rust term cap, checked before the block is mutated. -/
def capOk (ch : Chunk) (g : Gate) : Bool :=
  match rotAngle g with
  | some (_, q) =>
      match localIdx ch.qubits q with
      | some i =>
          match ch.state.parity[i]? with
          | some p => (ch.state.terms.get? p).isSome || ch.state.terms.size < ch.maxTerms
          | none => false
      | none => false
  | none => true

/-- Absorb `g` into a chunk that already covers its wires. -/
def feedInto (ch : Chunk) (g : Gate) : Option (List Gate × Chunk) :=
  if !ch.capOk g then none
  else
    match feedGateFast ch.qubits ch.state g with
    | some st => some ([], { ch with state := st, original := ch.original ++ [g] })
    | none => none

/-- Absorb `gate` if the block can take it; `none` when the block must be flushed first. -/
def feedTry (ch : Chunk) (g : Gate) : Option (List Gate × Chunk) :=
  match interpretable g with
  | none => none
  | some ops =>
      if !(ops.all fun q => decide (q < ch.numQubits)) then none
      else if !ch.admits ops then none
      else feedInto (ops.foldl Chunk.extend ch) g

/-- Feed one gate: absorb it, or flush the block and emit the gate verbatim. -/
def feed (ch : Chunk) (g : Gate) : List Gate × Chunk :=
  match ch.feedTry g with
  | some res => res
  | none =>
      if ch.qubits.isEmpty && ch.original.isEmpty then ([g], ch)
      else
        match ch.reset.feedTry g with
        | some (out, ch') => (ch.flush ++ out, ch')
        | none => (ch.flush ++ [g], ch.reset)

end Chunk

/-- Feed a gate list through a chunk, emitting as blocks close. -/
def runGates : Chunk → List Gate → List Gate × Chunk
  | ch, [] => ([], ch)
  | ch, g :: gs =>
      let (emitted, ch') := ch.feed g
      let (rest, ch'') := runGates ch' gs
      (emitted ++ rest, ch'')

/-- The gate-level `CnotMin` transformation. -/
def cnotMinGates (numQubits maxQubits maxTerms : Nat) (gs : List Gate) : List Gate :=
  let (out, ch) := runGates (Chunk.empty numQubits maxQubits maxTerms) gs
  out ++ ch.flush

/-- The Rust `MAX_CHUNK_QUBITS`, scaled down: `Parity` here is an unbounded `Nat`, so the cap
is purely the speed control the Rust comment describes. -/
def maxChunkQubits : Nat := 128

/-- The Rust `MAX_CHUNK_TERMS`. -/
def maxChunkTerms : Nat := 512

/-- `cnot_min` at the default bounds. -/
def cnotMin (numQubits : Nat) (gs : List Gate) : List Gate :=
  cnotMinGates numQubits maxChunkQubits maxChunkTerms gs

end TzapLean
