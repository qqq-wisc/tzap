import TzapLean.ExactMat

/-!
# The Bounded Clifford+T Synthesis Table

`SuperOpt` asks a question — "is there a shorter circuit with this unitary?" — that has the
same answer every time it is asked of the same unitary. So the answer is computed once, for
every unitary reachable within the bounds, and stored: this is `src/super_opt/table.rs`.

The table is built breadth-first, one gate at a time. Because layers are visited in
gate-count order, **the first circuit to reach a unitary is a smallest one**, so a table hit
*is* the synthesis answer — no search at lookup time. Two prunes keep the frontier small
without losing any unitary: a child never follows its parent's inverse, and among
qubit-disjoint neighbours only the canonically ordered interleaving is expanded.

The enumeration needs a *key*: a matrix, canonicalized so that any two representatives of the
same operator hash alike. Normalize the `√2` denominator, rotate to the canonical global
phase, flatten the coefficients, then hash them to the 64-bit fingerprint the table stores.
As in Rust, the pass independently re-verifies a candidate before rewriting, so a collision
can cost an optimization but cannot make the output wrong. The table is a source of
candidates and never load-bearing for correctness.

Circuits are stored prefix-shared, as in `synthesis_arena.rs`: BFS only ever extends a
circuit by one gate, so the entries form a tree and each node records just its last gate and
its parent.
-/

namespace TzapLean

/-! ## Canonical keys -/

namespace Cyc

/-- Whether a coefficient tuple is zero. -/
def isZero (x : Cyc) : Bool := x.a == 0 && x.b == 0 && x.c == 0 && x.d == 0

/-- The coefficients, in order. -/
def toList (x : Cyc) : List Int := [x.a, x.b, x.c, x.d]

end Cyc

/-- Lexicographic comparison on coefficient tuples. -/
def lexLt : List Int → List Int → Bool
  | [], [] => false
  | [], _ :: _ => true
  | _ :: _, [] => false
  | a :: as, b :: bs => if a < b then true else if b < a then false else lexLt as bs

/-- Multiply by `ω^p` in one step, as `src/super_opt/matrix.rs` does.

`Cyc.timesOmega` recurses `p` times, allocating a tuple at each step; in the table builder it
runs once per matrix entry and was four fifths of the enumeration's cost. Reducing `p` mod 8
and permuting once is the same function — `timesOmegaFast_eq` checks that exhaustively — and
turns `O(p)` allocations into one. -/
def Cyc.timesOmegaFast (x : Cyc) (p : Nat) : Cyc :=
  match p % 8 with
  | 0 => ⟨x.a, x.b, x.c, x.d⟩
  | 1 => ⟨-x.d, x.a, x.b, x.c⟩
  | 2 => ⟨-x.c, -x.d, x.a, x.b⟩
  | 3 => ⟨-x.b, -x.c, -x.d, x.a⟩
  | 4 => ⟨-x.a, -x.b, -x.c, -x.d⟩
  | 5 => ⟨x.d, -x.a, -x.b, -x.c⟩
  | 6 => ⟨x.c, x.d, -x.a, -x.b⟩
  | _ => ⟨x.b, x.c, x.d, -x.a⟩

/-- Lexicographic comparison of two coefficient tuples, without building lists. -/
def lexLtCyc (x y : Cyc) : Bool :=
  if x.a != y.a then decide (x.a < y.a)
  else if x.b != y.b then decide (x.b < y.b)
  else if x.c != y.c then decide (x.c < y.c)
  else decide (x.d < y.d)

/-- The power of `ω` making the first nonzero entry of a flat array lexicographically least. -/
def canonPhaseFlat (a : Array Cyc) : Nat :=
  match a.find? (fun x => !x.isZero) with
  | none => 0
  | some pivot =>
      (List.range 8).foldl
        (fun best p => if lexLtCyc (pivot.timesOmegaFast p) (pivot.timesOmegaFast best) then p else best) 0

/-- Mix one integer into a running FNV-style hash. -/
def mixHash (h : UInt64) (x : Int) : UInt64 :=
  let v : UInt64 := (((x % 4294967296) + 4294967296).toNat).toUInt64
  let h := (h ^^^ v) * 0x100000001B3
  (h ^^^ (h >>> 29)) * 0xBF58476D1CE4E5B9

/-- The hash of an already-normalized flat array. -/
def fingerprintFlat (den : Nat) (arr : Array Cyc) : UInt64 :=
  let p := canonPhaseFlat arr
  arr.foldl
    (fun h x =>
      let y := x.timesOmegaFast p
      mixHash (mixHash (mixHash (mixHash h y.a) y.b) y.c) y.d)
    (mixHash 0xcbf29ce484222325 (den : Int))

namespace ExactMat

/-- A basis state's row number: wire `0` is the most significant bit, as in Rust. Used only
to index the builder's flat arrays — nothing downstream depends on the numbering.

A loop rather than a fold over `List.finRange n`, which allocated a list on every call: this
runs twice per matrix entry, so it showed up. -/
def idx {n : Nat} (b : Basis n) : Nat := Id.run do
  let mut acc := 0
  for h : j in [0:n] do
    acc := 2 * acc + (if b ⟨j, h.2.1⟩ then 1 else 0)
  return acc

/-- Basis states of width `n`, materialized once per traversal rather than per entry. -/
def basisArray (n : Nat) : Array (Basis n) := (basisList n).toArray

/-- The entries in row-major `idx` order.

A matrix here is a *function*, so an entry re-walks every gate applied so far. That is fine
for the handful of matrices the pass builds and hopeless for the hundreds of thousands the
table enumerates, where each is one gate deeper than the last. Flattening once per candidate
keeps the chain one level deep.

Everything from here to `fingerprint` is used only by the table builder, which is unverified:
a wrong answer here costs table hits, never correctness, since `accepts` re-derives every
candidate's matrix from scratch before a rewrite is taken. -/
def flatten {n : Nat} (M : ExactMat n) : Array Cyc := Id.run do
  let bs := basisArray n
  let dim := 2 ^ n
  let mut a : Array Cyc := Array.replicate (dim * dim) 0
  for r in bs do
    for c in bs do
      a := a.set! (idx r * dim + idx c) (M.get r c)
  return a

/-- The matrix a flat array denotes. -/
def ofFlat (n : Nat) (den : Nat) (a : Array Cyc) : ExactMat n :=
  let dim := 2 ^ n
  { den, get := fun o i => a[idx o * dim + idx i]! }

/-- Strip common factors of `√2` from a flat array. -/
def normalizeFlat (den : Nat) (a : Array Cyc) : Nat × Array Cyc := Id.run do
  let mut d := den
  let mut arr := a
  for _ in [0:den] do
    if d == 0 || !arr.all (fun x => x.divisibleBySqrt2) then break
    arr := arr.map (·.divSqrt2)
    d := d - 1
  return (d, arr)

/-- **The table key**: denominator and phase both canonicalized, coefficients flattened. Two
gate lists with the same key denote the same operator up to global phase. This is the readable
statement of what `fingerprint` hashes; the builder never materializes it. -/
def key {n : Nat} (M : ExactMat n) : List Int :=
  let (d, arr) := normalizeFlat M.den (flatten M)
  let p := canonPhaseFlat arr
  (d : Int) :: (arr.toList.flatMap fun x => (x.timesOmegaFast p).toList)

/-- A 64-bit hash of the canonical key — what the table is indexed by, as in Rust.

Safe for the same reason Rust's is: a hit is only ever a *candidate*, and `accepts` recomputes
the replacement's matrix and compares it exactly before any rewrite is taken. A collision
costs a missed optimization, never a wrong one. -/
def fingerprint {n : Nat} (M : ExactMat n) : UInt64 :=
  let (d, arr) := normalizeFlat M.den (flatten M)
  fingerprintFlat d arr

/-- Flatten, normalize, hash, and hand back the flattened matrix — one traversal serving both
the table key and the next frontier, where the builder previously made two. -/
def keyedReify {n : Nat} (M : ExactMat n) : UInt64 × ExactMat n :=
  let (d, arr) := normalizeFlat M.den (flatten M)
  (fingerprintFlat d arr, ofFlat n d arr)

end ExactMat

/-! ## The library gate set

Deliberately **not** every gate the pass can read. `ccx` and `cz` are excluded so
superoptimization never *introduces* them: a Toffoli costs about seven `T` once the pipeline
lowers it, and `cz` would leave the `H`/`X`/`Z`/`S`/`T`/`CX` emission basis. Windows
*containing* those gates are still matched and simplified — their unitaries come from the
input — but such gates are never emitted. -/

/-- A gate the table may emit. -/
inductive LibGate where
  /-- Pauli `X`. -/
  | x (q : Nat)
  /-- Hadamard. -/
  | h (q : Nat)
  /-- Phase `S`. -/
  | s (q : Nat)
  /-- Inverse phase. -/
  | sdg (q : Nat)
  /-- Pauli `Z`. -/
  | z (q : Nat)
  /-- `T`. -/
  | t (q : Nat)
  /-- Inverse `T`. -/
  | tdg (q : Nat)
  /-- Controlled `X`. -/
  | cnot (control target : Nat)
deriving DecidableEq, Repr, Inhabited, Ord, Hashable

namespace LibGate

/-- The circuit gate a library gate stands for. -/
def toGate : LibGate → Gate
  | .x q => .x q
  | .h q => .h q
  | .s q => .s q
  | .sdg q => .sdg q
  | .z q => .z q
  | .t q => .t q
  | .tdg q => .tdg q
  | .cnot c tgt => .cnot c tgt

/-- The wires a library gate touches. -/
def qubits : LibGate → List Nat
  | .x q | .h q | .s q | .sdg q | .z q | .t q | .tdg q => [q]
  | .cnot c tgt => [c, tgt]

/-- Whether two library gates share no wire. -/
def isDisjoint (a b : LibGate) : Bool := a.qubits.all fun q => !b.qubits.contains q

/-- Whether `b` undoes `a` — the first prune: a child never follows its parent's inverse,
since the product would revisit the grandparent's unitary. -/
def isInverseOf : LibGate → LibGate → Bool
  | .s q, .sdg r | .sdg q, .s r | .t q, .tdg r | .tdg q, .t r => q == r
  | .x q, .x r | .h q, .h r | .z q, .z r => q == r
  | .cnot c₁ t₁, .cnot c₂ t₂ => c₁ == c₂ && t₁ == t₂
  | _, _ => false

end LibGate

/-- Every library gate on `k` wires: seven one-wire gates per wire, then every ordered pair
of distinct wires as a `cnot`. -/
def libGates (k : Nat) : List LibGate :=
  (List.range k).flatMap
      (fun q => [.x q, .h q, .s q, .sdg q, .z q, .t q, .tdg q]) ++
    (List.range k).flatMap fun c =>
      (List.range k).filterMap fun tgt => if c == tgt then none else some (.cnot c tgt)

/-! ## The builder's matrices

`ExactMat` is indexed by `Basis n`, which is right for the proofs and wrong for enumerating a
table: every entry access goes through a function and an index computation. The builder
therefore works on a flat array with the row operations of `src/super_opt/matrix.rs`, exactly
as Rust does — row `r` encodes the basis state with wire `0` most significant.

None of this is verified. A wrong answer here costs table hits, never correctness: `accepts`
re-derives every candidate's matrix through `ExactMat` and compares it exactly before a
rewrite is taken. -/

/-- A `2ⁿ × 2ⁿ` matrix as a flat row-major array over one `√2^den`. -/
structure FlatMat where
  /-- Wire count. -/
  n : Nat
  /-- Shared denominator exponent. -/
  den : Nat
  /-- Entries, row-major. -/
  data : Array Cyc
deriving Inhabited

namespace FlatMat

/-- Side length. -/
def dim (M : FlatMat) : Nat := 2 ^ M.n

/-- The row bit for wire `q`, with wire `0` most significant. -/
def bit (n q : Nat) : Nat := 1 <<< (n - 1 - q)

/-- The identity. -/
def id (n : Nat) : FlatMat := Id.run do
  let d := 2 ^ n
  let mut a : Array Cyc := Array.replicate (d * d) 0
  for i in [0:d] do
    a := a.set! (i * d + i) 1
  return { n, den := 0, data := a }

/-- Scale every row whose bits cover `mask` by `ω^p`. -/
def phaseMask (M : FlatMat) (mask p : Nat) : FlatMat := Id.run do
  let d := M.dim
  let mut a := M.data
  for r in [0:d] do
    if r &&& mask == mask then
      for c in [0:d] do
        a := a.set! (r * d + c) ((a[r * d + c]!).timesOmega p)
  return { M with data := a }

/-- Flip `tgt` in every row whose bits cover `ctrl` — `x`, `cnot`, `ccx`. -/
def xMask (M : FlatMat) (ctrl tgt : Nat) : FlatMat := Id.run do
  let d := M.dim
  let mut a := M.data
  for r in [0:d] do
    if r &&& ctrl == ctrl && r &&& tgt == 0 then
      let r' := r ||| tgt
      for c in [0:d] do
        let u := a[r * d + c]!
        let v := a[r' * d + c]!
        a := (a.set! (r * d + c) v).set! (r' * d + c) u
  return { M with data := a }

/-- Add and subtract row pairs across `q`, and bump the denominator. -/
def hadamard (M : FlatMat) (q : Nat) : FlatMat := Id.run do
  let d := M.dim
  let b := bit M.n q
  let mut a := M.data
  for r in [0:d] do
    if r &&& b == 0 then
      let r' := r ||| b
      for c in [0:d] do
        let u := a[r * d + c]!
        let v := a[r' * d + c]!
        a := (a.set! (r * d + c) (u + v)).set! (r' * d + c) (u - v)
  return { M with den := M.den + 1, data := a }

/-- One library gate, applied on the left. -/
def applyLib (M : FlatMat) : LibGate → FlatMat
  | .x q => M.xMask 0 (bit M.n q)
  | .h q => M.hadamard q
  | .s q => M.phaseMask (bit M.n q) 2
  | .sdg q => M.phaseMask (bit M.n q) 6
  | .z q => M.phaseMask (bit M.n q) 4
  | .t q => M.phaseMask (bit M.n q) 1
  | .tdg q => M.phaseMask (bit M.n q) 7
  | .cnot c tgt => M.xMask (bit M.n c) (bit M.n tgt)

/-- One circuit gate, applied on the left — the same actions as `ExactMat.applyGate`, on the
flat representation. `none` for the gates `SuperOpt` treats as window barriers. -/
def applyGate (M : FlatMat) : Gate → Option FlatMat
  | .x q => some (M.xMask 0 (bit M.n q))
  | .h q => some (M.hadamard q)
  | .s q => some (M.phaseMask (bit M.n q) 2)
  | .sdg q => some (M.phaseMask (bit M.n q) 6)
  | .z q => some (M.phaseMask (bit M.n q) 4)
  | .t q => some (M.phaseMask (bit M.n q) 1)
  | .tdg q => some (M.phaseMask (bit M.n q) 7)
  | .cnot c tgt => some (M.xMask (bit M.n c) (bit M.n tgt))
  | .ccx a b tgt => some (M.xMask (bit M.n a ||| bit M.n b) (bit M.n tgt))
  | .cz c tgt => some (M.phaseMask (bit M.n c ||| bit M.n tgt) 4)
  | .ccz a b tgt => some (M.phaseMask (bit M.n a ||| bit M.n b ||| bit M.n tgt) 4)
  | _ => none

/-- The matrix of a gate list on `k` wires. -/
def ofGates (k : Nat) (gs : List Gate) : Option FlatMat :=
  gs.foldl (fun acc g => acc.bind fun M => M.applyGate g) (some (FlatMat.id k))

/-- Strip common factors of `√2`. -/
def normalize (M : FlatMat) : FlatMat := Id.run do
  let mut d := M.den
  let mut arr := M.data
  for _ in [0:M.den] do
    if d == 0 || !arr.all (fun x => x.divisibleBySqrt2) then break
    arr := arr.map (·.divSqrt2)
    d := d - 1
  return { M with den := d, data := arr }

end FlatMat

/-! ## The table -/

/-- How far the table is built. -/
structure SuperOptTableConfig where
  /-- Widest table width, in wires. -/
  maxQubits : Nat := 2
  /-- Deepest circuit the enumeration reaches. -/
  maxGates : Nat := 4
  /-- Cap on stored unitaries per width, so a build always terminates promptly. -/
  maxEntriesPerQubit : Nat := 200000
deriving Repr, DecidableEq, Hashable

/-- One stored circuit: its last gate plus the node holding the rest. The root — the empty
circuit, i.e. the identity — has neither. -/
structure CircuitNode where
  /-- The node holding everything but the last gate. -/
  parent : Nat
  /-- The last gate. -/
  gate : LibGate
deriving Repr, Inhabited

/-- One width of the table, stored as a prefix-sharing arena. -/
structure WidthTable where
  /-- Fingerprint of each stored unitary, mapped to its node. -/
  keys : Std.HashMap UInt64 Nat := ∅
  /-- The arena; node `0` is the root. -/
  nodes : Array (Option CircuitNode) := #[none]
  /-- Whether the entry cap stopped the build. -/
  saturated : Bool := false
  /-- The last depth completed. -/
  depth : Nat := 0
deriving Inhabited

namespace WidthTable

/-- Recover a stored circuit by walking to the root. -/
def circuitOf (w : WidthTable) (node : Nat) : List LibGate :=
  go w.nodes.size node []
where
  /-- Walk up the parent chain, accumulating gates front-first. -/
  go : Nat → Nat → List LibGate → List LibGate
    | 0, _, acc => acc
    | fuel + 1, i, acc =>
        match w.nodes[i]? with
        | some (some nd) => go fuel nd.parent (nd.gate :: acc)
        | _ => acc

/-- How many unitaries this width stores. -/
def size (w : WidthTable) : Nat := w.keys.size

end WidthTable

/-- Build one width breadth-first. Layers are visited in gate-count order, so the first
circuit reaching a unitary is a smallest one. -/
def buildWidth (k : Nat) (cfg : SuperOptTableConfig) : WidthTable := Id.run do
  let idM := FlatMat.id k
  let gates := libGates k
  let mut tbl : WidthTable :=
    { keys := (∅ : Std.HashMap UInt64 Nat).insert (fingerprintFlat idM.den idM.data) 0,
      nodes := #[none] }
  let mut frontier : Array (Nat × FlatMat) := #[(0, idM)]
  let mut stop := false
  for depth in [1 : cfg.maxGates + 1] do
    if stop then break
    let mut accepted : Array (Nat × FlatMat) := #[]
    for (parent, base) in frontier do
      if stop then break
      let last : Option LibGate := (tbl.nodes[parent]?.join).map (·.gate)
      for g in gates do
        if stop then break
        -- the two prunes
        let pruned : Bool :=
          match last with
          | some l => l.isInverseOf g || (l.isDisjoint g && compare g l == .lt)
          | none => false
        if pruned then continue
        let child := (base.applyLib g).normalize
        let ky := fingerprintFlat child.den child.data
        if tbl.keys.contains ky then pure ()
        else if tbl.size ≥ cfg.maxEntriesPerQubit then
          tbl := { tbl with saturated := true }
          stop := true
        else
          let node := tbl.nodes.size
          tbl := { tbl with keys := tbl.keys.insert ky node,
                            nodes := tbl.nodes.push (some ⟨parent, g⟩) }
          accepted := accepted.push (node, child)
    tbl := { tbl with depth := depth }
    if accepted.isEmpty then break
    frontier := accepted
  return tbl

/-- The synthesis table: one `WidthTable` per width, indexed by wire count. -/
structure SynthTable where
  /-- Widths `0 … maxQubits`; index `0` is unused. -/
  widths : Array WidthTable
deriving Inhabited

/-- Build the table for a configuration.

Width `0` is built like the rest even though nothing ever queries it. There are no gates on
zero wires, so it costs a single arena node — and it keeps every `WidthTable` self-consistent
(a node for every key, a key for every node), which is what makes serializing and reading a
table back the identity. Filling it with `default` instead left one node with no key, and the
round trip did not reproduce it. -/
def buildTable (cfg : SuperOptTableConfig) : SynthTable :=
  { widths := (Array.range (cfg.maxQubits + 1)).map fun k => buildWidth k cfg }

/-- Whether the table holds anything for this fingerprint — the cheap pre-filter `SuperOpt`
uses to decide whether the verified lookup is worth running at all. Unverified in both
directions: a false negative costs an optimization, and a false positive costs one wasted
verified lookup. -/
def SynthTable.mayHold (tbl : SynthTable) (k : Nat) (M : FlatMat) : Bool :=
  match tbl.widths[k]? with
  | none => false
  | some w =>
      let N := M.normalize
      (w.keys.get? (fingerprintFlat N.den N.data)).isSome

/-- Does the table hold a circuit for this unitary *strictly shorter* than `len`?

The question `SuperOpt` actually needs answered before paying for a verified lookup. Asking
only whether the table holds the unitary at all is useless: it holds essentially every short
Clifford+T circuit, so the answer is almost always yes. -/
def SynthTable.hasShorter (tbl : SynthTable) (k : Nat) (M : FlatMat) (len : Nat) : Bool :=
  match tbl.widths[k]? with
  | none => false
  | some w =>
      let N := M.normalize
      match w.keys.get? (fingerprintFlat N.den N.data) with
      | none => false
      | some node => (w.circuitOf node).length < len

/-- Look a unitary up from the flat matrix the scan carries, rather than from an `ExactMat`.

The same lookup as `synthesize`, on the representation the search actually holds — the scan
never builds an `ExactMat`, because the checker that vets its answer does that once per
*selected* rewrite instead of once per window examined. Unverified, like every other part of
the search. -/
def SynthTable.synthesizeFlat (tbl : SynthTable) (k : Nat) (M : FlatMat) : Option (List Gate) :=
  match tbl.widths[k]? with
  | none => none
  | some w =>
      let N := M.normalize
      match w.keys.get? (fingerprintFlat N.den N.data) with
      | none => none
      | some node => some ((w.circuitOf node).map LibGate.toGate)

/-- Look a unitary up. A hit is the shortest circuit the enumeration found for it; the
caller still re-verifies before rewriting. -/
def SynthTable.synthesize (tbl : SynthTable) (k : Nat) (M : ExactMat k) : Option (List Gate) :=
  match tbl.widths[k]? with
  | none => none
  | some w =>
      match w.keys.get? M.fingerprint with
      | none => none
      | some node => some ((w.circuitOf node).map LibGate.toGate)

end TzapLean
