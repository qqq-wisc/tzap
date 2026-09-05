import TzapLean.ExactMat
import TzapLean.Locality
import TzapLean.SynthTable
import TzapLean.Rewrite

/-!
# `SuperOpt`: the Algorithm

A forward scan that carves out small *windows* — causally connected groups of gates on a few
wires — computes each window's exact matrix, and replaces it whenever a shorter gate list has
the same matrix. Nothing is matched syntactically: any identity expressible in the search
space is found without ever being written down.

Candidates come from the precomputed synthesis table (`SynthTable.lean`): a window's matrix
is canonicalized to a key and looked up, and a hit *is* the shortest circuit the enumeration
found for that unitary. The table is unverified — its BFS, its prunes, its key, its reified
matrices are all outside the proof — because every candidate is re-verified by `accepts`
before it is taken. A wrong table costs optimization, never correctness.

This is `src/super_opt/mod.rs`: one forward scan with every window open at once, greedy
selection over the whole circuit, and one splice at the end.

Windows are subsequences, not slices — gates in between that share no wire with the window are
left where they are — and the gates a window *does* claim are scattered among them. So a
window is held as a set of gate indices, and a selection as a *tagging* of the gate list. The
soundness of splicing them all in at once, however they interleave, is `applyAll_correct` in
`TzapLean/Rewrite.lean`.

Nothing in this file is verified. It proposes; `TzapLean/SuperOptProof.lean` checks.
-/

namespace TzapLean

/-- Limits on window growth and on how far the search looks. -/
structure SuperOptConfig where
  /-- Widest window, in wires. -/
  maxQubits : Nat := 2
  /-- Longest window, in gates. -/
  maxWindow : Nat := 6
deriving Repr

/-- Gates a window may contain: `rz` is outside the exact Clifford+T domain, and
`measure`/`reset` are not unitary. -/
def isWindowGate (g : Gate) : Bool :=
  match g with
  | .rz _ _ | .measure _ _ | .reset _ => false
  | _ => true

/-! ## The per-qubit index

Closing a window over a newly bridged wire means finding the earlier gates on that wire, and
the scan keeps `gatesByQubit[q]` — the indices of the gates seen so far on wire `q`, ascending
— so they are looked up rather than scanned for. It is Rust's `gates_by_qubit`, built as the
scan goes so that a window can never bridge across a gate the scan has not reached.

A wire's list is sorted, so the gates in a window's span are a contiguous slice of it, found
by binary search. -/

/-- Index of the first entry of a sorted array greater than or equal to `x`. -/
def lowerBoundIdx (a : Array Nat) (x : Nat) : Nat :=
  go a.size 0 a.size
where
  /-- Binary search, with the interval width as fuel. -/
  go : Nat → Nat → Nat → Nat
  | 0, lo, _ => lo
  | fuel + 1, lo, hi =>
      if lo < hi then
        let mid := (lo + hi) / 2
        if a[mid]! < x then go fuel (mid + 1) hi else go fuel lo mid
      else lo

/-! ## Windows

A window is the connected closure of its anchor: the gates that reach the anchor through
shared wires, everything else disjoint from its wires and so commuting past it. A gate that
brings in a *new* wire pulls in the earlier gates on that wire, retroactively — Rust's
`expand_component_closure`.

The scan holds a window as a set of *gate indices* rather than a list, because that is what
lets every window stay open at once: Rust keeps one closed component per anchor, extends every
window a gate touches, and lets the first one offered with a shorter replacement claim its
gates. Selecting the window that *completes* earliest, rather than the one that *starts*
earliest, is the whole difference — and it is only expressible if the windows are all alive
together. -/

/-- A window still growing: its anchor, the wires it covers, and the gates it holds, by index
and ascending. -/
structure LiveWin where
  /-- The gate the window started from. -/
  anchor : Nat
  /-- The wires it covers. -/
  support : List Qubit
  /-- The gates it holds, ascending. -/
  members : List Nat
deriving Inhabited

/-- Insert into an ascending list, without duplicates. -/
def insertSorted (x : Nat) : List Nat → List Nat
  | [] => [x]
  | y :: ys => if x < y then x :: y :: ys else if x == y then y :: ys else y :: insertSorted x ys

/-- Admit a wire into a support, recording it as pending if it is new. -/
def admit (sup : List Qubit) (pending : List Qubit) (q : Qubit) : List Qubit × List Qubit :=
  if sup.contains q then (sup, pending) else (sup ++ [q], q :: pending)

/-- The gates on wire `q` with index in `[lo, hi]`, ascending. -/
def onWireIn (gbq : Array (Array Nat)) (q lo hi : Nat) : List Nat :=
  match gbq[q]? with
  | none => []
  | some a => Id.run do
      let mut out : List Nat := []
      for d in [lowerBoundIdx a lo : a.size] do
        let j := a[d]!
        if j > hi then break
        out := j :: out
      return out.reverse

/-- One wire's worth of closure: absorb every gate on it between the anchor and the current
gate, admitting their wires in turn. -/
def closeWire (gs : Array Gate) (gbq : Array (Array Nat)) (maxQ maxG anchor i q : Nat)
    (sup : List Qubit) (mem : List Nat) (pending : List Qubit) :
    Option (List Qubit × List Nat × List Qubit) :=
  (onWireIn gbq q anchor i).foldl
    (fun acc j =>
      match acc with
      | none => none
      | some (sup, mem, pending) =>
          if mem.contains j then some (sup, mem, pending)
          else if !isWindowGate gs[j]! then none
          else
            let mem := insertSorted j mem
            if mem.length > maxG then none
            else
              let (sup, pending) :=
                gs[j]!.qubitsOf.foldl (fun (sp : List Qubit × List Qubit) r =>
                  admit sp.1 sp.2 r) (sup, pending)
              if sup.length > maxQ then none else some (sup, mem, pending))
    (some (sup, mem, pending))

/-- Close the window over its wires, one pending wire at a time. -/
def closeLoop (gs : Array Gate) (gbq : Array (Array Nat)) (maxQ maxG anchor i : Nat) :
    Nat → List Qubit → List Nat → List Qubit → Option (List Qubit × List Nat)
  | 0, _, _, _ => none
  | _, sup, mem, [] => some (sup, mem)
  | fuel + 1, sup, mem, q :: pending =>
      match closeWire gs gbq maxQ maxG anchor i q sup mem pending with
      | none => none
      | some (sup, mem, pending) => closeLoop gs gbq maxQ maxG anchor i fuel sup mem pending

/-- **Extend a window with gate `i` and re-close it.** `none` when the result would exceed the
bounds, or would have to absorb a gate a window may not hold. -/
def expandClosure (gs : Array Gate) (gbq : Array (Array Nat)) (maxQ maxG : Nat)
    (w : LiveWin) (i : Nat) : Option (List Qubit × List Nat) :=
  let mem := insertSorted i w.members
  if mem.length > maxG then none
  else
    let (sup, pending) :=
      gs[i]!.qubitsOf.foldl (fun (sp : List Qubit × List Qubit) r => admit sp.1 sp.2 r)
        (w.support, ([] : List Qubit))
    if sup.length > maxQ then none
    else closeLoop gs gbq maxQ maxG w.anchor i (maxG * maxQ + maxQ + 2) sup mem pending

/-! ## Proposing and verifying a replacement -/

/-- Whether a candidate is usable *and* really has the window's matrix, up to global phase.
This is the check the correctness proof consumes; the search that proposes candidates is
unverified, so everything the proof needs is re-established here. -/
def accepts {k : Nat} (target : ExactMat k) (cand : List Gate) : Bool :=
  cand.all (fun g => g.qubitsOf.all (fun q => q < k) && decide g.Wf) &&
    (match ExactMat.matrixOf k cand with
     | none => false
     | some N => (ExactMat.phaseMatch target N.normalize).isSome)

/-- Rename a local circuit back to the window's physical wires. -/
def globalizeGate (S : List Qubit) : Gate → Gate := mapQubits (fun i => S.getD i 0)

/-- **Verify a proposed replacement.** Everything the correctness proof needs is
re-established here, against the gates as they actually are: the wires are distinct and belong
to the register, the gates the rewrite claims live on them and have distinct operands, the
replacement is strictly shorter, and its matrix is the window's up to a global phase.

The search that proposed the candidate is unverified — it is a table lookup on a flat matrix,
with no `ExactMat` anywhere — and this is what makes that safe. It runs once per *selected*
rewrite, where the old arrangement built an exact matrix for every window that got past a
filter. -/
def checkRewrite (n : Nat) (S : List Qubit) (members cand : List Gate) : Bool :=
  decide S.Nodup && S.all (fun q => q < n) &&
    members.all (fun g => g.qubitsOf.all (fun q => S.contains q) && decide g.Wf) &&
    decide (cand.length < members.length) &&
    (match ExactMat.matrixOf S.length (localizeGates S members) with
     | none => false
     | some M => accepts M.normalize cand)

/-! ## The scan

One forward pass, with every window open at once — `SuperOpt::run`.

For each gate, in order: the windows registered on its wires are collected oldest-anchor
first; each is extended and re-closed; the first with a strictly shorter replacement claims
its gates, and anything overlapping a claimed gate is refused afterwards. The gate then
anchors a window of its own. Selecting the window that *completes* earliest, rather than the
one that *starts* earliest, is what this arrangement buys, and it is the reason every window
has to be alive together.

**None of this is verified.** What comes out is a *tagging* — each gate labelled with the
rewrite that claims it. `checkRewrite` vets each selected rewrite, then `onSuppB` and `sepB`
vet the whole tagging before anything is spliced. A bug here costs an optimization, or a
rewrite that is refused, never a wrong circuit. -/

/-- Whether a gate may anchor a window. -/
def canAnchor (cfg : SuperOptConfig) (n : Nat) (g : Gate) : Bool :=
  isWindowGate g && g.qubitsOf.length ≤ cfg.maxQubits &&
    g.qubitsOf.all (fun q => q < n) && decide g.Wf

/-! ## Canonical window-shape cache

Two windows with the same gate kinds and support-local operands have the same matrix,
regardless of their physical wires. Cache the synthesis answer for that exact shape, including
`none`, so a recurring shape does not rebuild and fingerprint its flat matrix. This is the
first layer of Rust's `MatrixStore`; transition states and persistence across scans are kept
separate so their effects can be measured independently. -/

/-- A compact, exact encoding of a window's support-local gate sequence. -/
structure ShapeKey where
  /-- Width is explicit so keys remain unambiguous even if a future scan carries idle wires. -/
  width : Nat
  /-- One 16-bit code per gate. -/
  gates : List UInt16
deriving BEq, Hashable

/-- Encode one support-local gate. Four bits per operand supports widths up to 16; larger
experimental configurations conservatively bypass the cache. -/
def compactGateCode : Gate → Option UInt16
  | .x q => one 0 q
  | .h q => one 1 q
  | .s q => one 2 q
  | .sdg q => one 3 q
  | .z q => one 4 q
  | .t q => one 5 q
  | .tdg q => one 6 q
  | .cnot c t => two 7 c t
  | .cz c t => two 8 c t
  | .ccx a b t => three 9 a b t
  | .ccz a b t => three 10 a b t
  | .rz _ _ | .measure _ _ | .reset _ => none
where
  one (tag q : Nat) : Option UInt16 :=
    if q < 16 then some (UInt16.ofNat (tag * 4096 + q)) else none
  two (tag a b : Nat) : Option UInt16 :=
    if a < 16 && b < 16 then some (UInt16.ofNat (tag * 4096 + a * 16 + b)) else none
  three (tag a b c : Nat) : Option UInt16 :=
    if a < 16 && b < 16 && c < 16 then
      some (UInt16.ofNat (tag * 4096 + a * 256 + b * 16 + c))
    else none

/-- Canonical key for physical gates `mem`, localized to `support`. Builds no localized gate
list, leaving that allocation and matrix construction to the cache-miss path. -/
def compactShapeKey (support : List Qubit) (gs : Array Gate) (mem : List Nat) : Option ShapeKey :=
  go mem []
where
  go : List Nat → List UInt16 → Option ShapeKey
    | [], acc => some ⟨support.length, acc.reverse⟩
    | i :: is, acc =>
        match compactGateCode (localizeGate support gs[i]!) with
        | none => none
        | some code => go is (code :: acc)

/-- The scan's state. -/
structure Scan where
  /-- Every window ever anchored, by id — ids ascend with anchor position. -/
  wins : Array LiveWin
  /-- Whether a window is still growing. -/
  alive : Array Bool
  /-- Wire → the windows registered on it. -/
  byQubit : Array (List Nat)
  /-- Wire → the gates seen on it so far, ascending. -/
  gbq : Array (Array Nat)
  /-- Whether a gate is claimed by a selected rewrite. -/
  claimed : Array Bool
  /-- Gate → the rewrite claiming it. -/
  tags : Array Claim
  /-- Rewrite → its wires. -/
  supports : Array (List Qubit)
  /-- Rewrite → its replacement, in support-local form. -/
  cands : Array (List Gate)
  /-- Canonical window shape → the table's answer; inner `none` is a cached miss. -/
  shapeCache : Std.HashMap ShapeKey (Option (List Gate))
  /-- Cache hits in this scan, retained for tests and profiling. -/
  shapeHits : Nat
  /-- Newly resolved canonical shapes in this scan. -/
  shapeMisses : Nat

/-- The state a scan starts from. -/
def Scan.initial (n : Nat) (count : Nat) : Scan where
  wins := #[]
  alive := #[]
  byQubit := Array.replicate n []
  gbq := Array.replicate n #[]
  claimed := Array.replicate count false
  tags := Array.replicate count none
  supports := #[]
  cands := #[]
  shapeCache := ∅
  shapeHits := 0
  shapeMisses := 0

/-- Retire a window. -/
def Scan.retire (st : Scan) (wid : Nat) : Scan := { st with alive := st.alive.set! wid false }

/-- Claim a window's gates for a new rewrite. -/
def Scan.select (st : Scan) (mem : List Nat) (sup : List Qubit) (cand : List Gate) : Scan :=
  let id := st.supports.size
  { st with
    claimed := mem.foldl (fun a j => a.set! j true) st.claimed
    tags := mem.foldl (fun a j => a.set! j (some id)) st.tags
    supports := st.supports.push sup
    cands := st.cands.push cand }

/-- Offer one cached or freshly synthesized answer to greedy selection. -/
def Scan.offer (st : Scan) (mem : List Nat) (sup : List Qubit)
    (answer : Option (List Gate)) : Scan × Bool :=
  match answer with
  | none => (st, false)
  | some cand =>
      if cand.length < mem.length then (st.select mem sup cand, true) else (st, false)

/-- Offer a window to the table: select it when the table holds something strictly shorter and
nothing it claims is claimed already. -/
def Scan.consider (st : Scan) (tbl : SynthTable) (gs : Array Gate) (sup : List Qubit)
    (mem : List Nat) : Scan × Bool :=
  if mem.any (st.claimed[·]!) then (st, false)
  else
    match compactShapeKey sup gs mem with
    | some key =>
        match st.shapeCache.get? key with
        | some answer =>
            let st := { st with shapeHits := st.shapeHits + 1 }
            st.offer mem sup answer
        | none =>
            let members := mem.map (gs[·]!)
            let answer :=
              (FlatMat.ofGates sup.length (localizeGates sup members)).bind
                (tbl.synthesizeFlat sup.length)
            let st := { st with shapeCache := st.shapeCache.insert key answer,
                                shapeMisses := st.shapeMisses + 1 }
            st.offer mem sup answer
    | none =>
        let members := mem.map (gs[·]!)
        let answer :=
          (FlatMat.ofGates sup.length (localizeGates sup members)).bind
            (tbl.synthesizeFlat sup.length)
        st.offer mem sup answer

/-- Offer one live window the current gate. -/
def Scan.step (cfg : SuperOptConfig) (tbl : SynthTable) (gs : Array Gate) (st : Scan)
    (i wid : Nat) : Scan :=
  let w := st.wins[wid]!
  if st.claimed[i]! || w.members.any (st.claimed[·]!) then st.retire wid
  else
    match expandClosure gs st.gbq cfg.maxQubits cfg.maxWindow w i with
    | none => st.retire wid
    | some (sup, mem) =>
        let (st, selected) := st.consider tbl gs sup mem
        if selected || mem.length ≥ cfg.maxWindow then st.retire wid
        else
          let added := sup.filter (fun q => !w.support.contains q)
          let st := { st with wins := st.wins.set! wid ⟨w.anchor, sup, mem⟩ }
          { st with byQubit := added.foldl (fun a q => a.modify q (wid :: ·)) st.byQubit }

/-- Everything one gate does to the scan. -/
def Scan.gate (cfg : SuperOptConfig) (tbl : SynthTable) (n : Nat) (gs : Array Gate)
    (st : Scan) (i : Nat) : Scan := Id.run do
  let g := gs[i]!
  let qs := g.qubitsOf
  -- the windows this gate touches, oldest anchor first; walking the lists compacts them
  let mut st := st
  let mut touched : List Nat := []
  for q in qs do
    if q < st.byQubit.size then
      let live := st.byQubit[q]!.filter (st.alive[·]!)
      st := { st with byQubit := st.byQubit.set! q live }
      for wid in live do
        touched := insertSorted wid touched
  -- the gate joins its wires' history, so no later window can bridge across it unseen
  for q in qs do
    if q < st.gbq.size then st := { st with gbq := st.gbq.modify q (·.push i) }
  if !isWindowGate g then
    -- `rz`, `measure` and `reset` kill every window on their wires
    for wid in touched do st := st.retire wid
    return st
  for wid in touched do
    st := Scan.step cfg tbl gs st i wid
  if canAnchor cfg n g && !st.claimed[i]! then
    let wid := st.wins.size
    st := { st with wins := st.wins.push ⟨i, qs, [i]⟩, alive := st.alive.push true }
    for q in qs do
      if q < st.byQubit.size then st := { st with byQubit := st.byQubit.modify q (wid :: ·) }
  return st

/-- **Propose a set of rewrites**: the tagging, and each rewrite's wires and replacement. -/
def proposeRewrites (cfg : SuperOptConfig) (tbl : SynthTable) (n : Nat) (gs : Array Gate) :
    Scan := Id.run do
  let mut st := Scan.initial n gs.size
  for i in [0 : gs.size] do
    st := Scan.gate cfg tbl n gs st i
  return st

/-! ## Checking a proposal, and taking it

The scan is unverified, so what it produces is a *proposal*, and every part of it is checked
before a gate moves. A rewrite whose check fails is untagged — the others still stand — and if
the two whole-tagging conditions then fail, nothing is rewritten at all.

In practice none of these fire. They are what makes the search's freedom safe: it can be
changed, tuned, or replaced wholesale without touching a proof. -/

/-- Each rewrite's replacement on the circuit's own wires. -/
def Scan.repl (st : Scan) (w : Nat) : List Gate :=
  (st.cands[w]?.getD []).map (globalizeGate (st.supports[w]?.getD []))

/-- Each rewrite's wires. -/
def Scan.supp (st : Scan) (w : Nat) : List Qubit := st.supports[w]?.getD []

/-- The tagging the scan produced, as gates. -/
def Scan.tagged (st : Scan) (gs : List Gate) : List Tagged :=
  gs.zipIdx.map fun (g, i) => (g, st.tags[i]?.getD none)

/-- Whether a rewrite's replacement checks out against the gates it claims.  `members` is a
batched `groupedClaim` lookup, so all rewrites share one traversal of the tagging. -/
def Scan.keep (st : Scan) (n : Nat) (members : Nat → List Gate) (w : Nat) : Bool :=
  checkRewrite n (st.supp w) (members w) (st.cands[w]?.getD [])

/-- Untag every rewrite whose replacement does not check out. -/
def Scan.vetted (st : Scan) (n : Nat) (xs : List Tagged) : List Tagged :=
  let groups := groupClaimsAux (claimBound xs) xs
  let members := fun w => groups[w]?.getD []
  vettedBy (st.keep n members) xs

/-- Peephole superoptimization of a gate list over `n` wires: propose, check, splice.

One forward scan, as `SuperOpt::run` is one forward scan — repetition is the driver's job. -/
def superOptGates (cfg : SuperOptConfig) (tbl : SynthTable) (n : Nat) (gs : List Gate) :
    List Gate :=
  let arr := gs.toArray
  let st := proposeRewrites cfg tbl n arr
  let xs := st.vetted n (st.tagged gs)
  if onSuppB st.supp xs && sepAllB st.supp xs then applyAllLinear st.repl [] xs else gs

/-- Peephole superoptimization of a circuit. -/
def superOpt (cfg : SuperOptConfig) (tbl : SynthTable) (c : Circuit) : Circuit :=
  c.withGates (superOptGates cfg tbl c.numQubits c.gates)

end TzapLean
