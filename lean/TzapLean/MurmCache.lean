import TzapLean.Murm

/-!
# Persisting the MURM

A port of the on-disk cache in `src/super_opt/table.rs`. Building the MURM is the one
expensive thing `SuperOpt` does, and the result depends on nothing but the bounds — so it is
written once to `~/.tzap-lean/superopt-tables/` and read back by later runs.

The layout follows Rust's: a magic number, a format version, the bounds the MURM was built
for, then one fixed-width record per arena node — fingerprint, parent, and last gate. The
fingerprint-to-node map is rebuilt on read rather than stored, since nodes and fingerprints
are one to one.

Every read is validated against the bounds it is being asked for, and **any failure — missing
file, truncated write, version bump, bounds mismatch — falls back to rebuilding.** A cache is
never trusted into a wrong answer; the worst a bad file can do is waste the read. Writes go to
a sibling temporary file and are renamed into place, so a reader never sees a half-written
MURM.

The version is in the file *name* as well as its header, so a bump cannot collide with old
files: they are simply never looked up again.
-/

namespace TzapLean

namespace MurmCache

/-- Format identifier. -/
def magic : List UInt8 := [0x4D, 0x55, 0x52, 0x4D]  -- "MURM"

/-- Bump whenever the byte layout below changes. -/
def formatVersion : UInt32 := 2

/-! ## Little-endian primitives -/

/-- Append a `UInt32`, little-endian. -/
def putU32 (acc : ByteArray) (v : UInt32) : ByteArray :=
  ((List.range 4).foldl (fun a i => a.push ((v >>> (8 * i.toUInt32)).toUInt8)) acc)

/-- Append a `UInt64`, little-endian. -/
def putU64 (acc : ByteArray) (v : UInt64) : ByteArray :=
  ((List.range 8).foldl (fun a i => a.push ((v >>> (8 * i.toUInt64)).toUInt8)) acc)

/-- Read a `UInt32` at `i`, little-endian. -/
def getU32 (b : ByteArray) (i : Nat) : Option UInt32 :=
  if i + 4 ≤ b.size then
    some ((List.range 4).foldl (fun a j => a ||| (b[i + j]!.toUInt32 <<< (8 * j.toUInt32))) 0)
  else none

/-- Read a `UInt64` at `i`, little-endian. -/
def getU64 (b : ByteArray) (i : Nat) : Option UInt64 :=
  if i + 8 ≤ b.size then
    some ((List.range 8).foldl (fun a j => a ||| (b[i + j]!.toUInt64 <<< (8 * j.toUInt64))) 0)
  else none

/-! ## Library gates as bytes

Four bytes: a tag plus up to three operands. `noGateTag` marks the
root, which has neither a gate nor a parent. -/

/-- Sentinel tag for the root node. -/
def noGateTag : UInt8 := 255

/-- A library gate's four bytes. -/
def gateBytes : LibGate → UInt8 × UInt8 × UInt8 × UInt8
  | .x q => (0, q.toUInt8, 0, 0)
  | .h q => (1, q.toUInt8, 0, 0)
  | .s q => (2, q.toUInt8, 0, 0)
  | .sdg q => (3, q.toUInt8, 0, 0)
  | .z q => (4, q.toUInt8, 0, 0)
  | .t q => (5, q.toUInt8, 0, 0)
  | .tdg q => (6, q.toUInt8, 0, 0)
  | .cnot c tgt => (7, c.toUInt8, tgt.toUInt8, 0)
  | .cz a b => (8, a.toUInt8, b.toUInt8, 0)
  | .ccx a b t => (9, a.toUInt8, b.toUInt8, t.toUInt8)
  | .ccz a b c => (10, a.toUInt8, b.toUInt8, c.toUInt8)

/-- The gate three bytes encode, or `none` for an unknown tag. -/
def gateOfBytes (tag a b c : UInt8) : Option LibGate :=
  match tag with
  | 0 => some (.x a.toNat)
  | 1 => some (.h a.toNat)
  | 2 => some (.s a.toNat)
  | 3 => some (.sdg a.toNat)
  | 4 => some (.z a.toNat)
  | 5 => some (.t a.toNat)
  | 6 => some (.tdg a.toNat)
  | 7 => some (.cnot a.toNat b.toNat)
  | 8 => some (.cz a.toNat b.toNat)
  | 9 => some (.ccx a.toNat b.toNat c.toNat)
  | 10 => some (.ccz a.toNat b.toNat c.toNat)
  | _ => none

/-! ## Writing -/

/-- One width: node count, then a record per node. -/
def putWidth (acc : ByteArray) (w : WidthMurm) : ByteArray :=
  -- fingerprint per node, recovered from the map (nodes and fingerprints are one to one)
  let fpOf : Array UInt64 := Id.run do
    let mut a : Array UInt64 := Array.replicate w.nodes.size 0
    for (fp, node) in w.keys.toList do
      if node < a.size then a := a.set! node fp
    return a
  let acc := putU64 acc w.nodes.size.toUInt64
  let acc := (List.range w.nodes.size).foldl (init := acc) fun a i =>
    let a := putU64 a (fpOf[i]!)
    match w.nodes[i]! with
    | none => (putU32 a 0xFFFFFFFF).push noGateTag |>.push 0 |>.push 0 |>.push 0
    | some nd =>
        let (t, x, y, z) := gateBytes nd.gate
        ((putU32 a nd.parent.toUInt32).push t).push x |>.push y |>.push z
  let acc := acc.push (if w.saturated then 1 else 0)
  putU32 acc w.depth.toUInt32

/-- Small deterministic checksum. It is corruption detection, not a trust boundary: every MURM
candidate is still checked by exact matrices before use. -/
def checksum (bytes : ByteArray) : UInt64 :=
  bytes.foldl (fun hash byte => (hash ^^^ byte.toUInt64) * 0x100000001B3) 0xcbf29ce484222325

/-- Serialize the checksummed payload of a whole MURM. -/
def serializePayload (cfg : MurmConfig) (murm : Murm) : ByteArray :=
  let acc := magic.foldl (fun a b => a.push b) ByteArray.empty
  let acc := putU32 acc formatVersion
  let acc := putU32 acc cfg.maxQubits.toUInt32
  let acc := putU32 acc cfg.maxGates.toUInt32
  let acc := putU64 acc cfg.maxEntriesPerQubit.toUInt64
  let acc := putU32 acc cfg.basis.mask.toUInt32
  let acc := putU32 acc murm.widths.size.toUInt32
  murm.widths.foldl (fun a w => putWidth a w) acc

/-- Serialize a whole MURM, including a checksum over its complete payload. -/
def serialize (cfg : MurmConfig) (murm : Murm) : ByteArray :=
  let payload := serializePayload cfg murm
  putU64 payload (checksum payload)

/-! ## Reading -/

/-- Check magic, version, and bounds. Returns the offset of the body. -/
def readHeader (b : ByteArray) (cfg : MurmConfig) : Option Nat := do
  if cfg.basis.isEmpty || !cfg.basis.isSubset supportedGateSet then none
  else if b.size < 28 then none
  else if (List.range 4).any (fun i => b[i]! != magic[i]!) then none
  else
    let ver ← getU32 b 4
    if ver != formatVersion then none
    else
      let q ← getU32 b 8
      let g ← getU32 b 12
      let e ← getU64 b 16
      let basis ← getU32 b 24
      if q != cfg.maxQubits.toUInt32 || g != cfg.maxGates.toUInt32
          || e != cfg.maxEntriesPerQubit.toUInt64 || basis != cfg.basis.mask.toUInt32 then none
      else some 28

/-- A decoded library gate is usable only when it belongs to this MURM and its operands form
a valid gate at this width. -/
def validGate (cfg : MurmConfig) (width : Nat) (g : LibGate) : Bool :=
  cfg.basis.contains g.kind && g.qubits.all (· < width) && decide g.toGate.Wf

/-- Read one width, returning it and the offset just past it. Parent links must point strictly
backward, and only node zero may be the root. -/
def getWidth (b : ByteArray) (i width : Nat) (cfg : MurmConfig) : Option (WidthMurm × Nat) := do
  let count ← getU64 b i
  let count := count.toNat
  let mut off := i + 8
  let mut nodes : Array (Option CircuitNode) := Array.emptyWithCapacity count
  let mut keys : Std.HashMap UInt64 Nat := ∅
  for idx in [0:count] do
    let some fp := getU64 b off | none
    let some parent := getU32 b (off + 8) | none
    if off + 16 > b.size then none
    let tag := b[off + 12]!
    let x := b[off + 13]!
    let y := b[off + 14]!
    let z := b[off + 15]!
    if tag == noGateTag then
      if idx != 0 || parent != 0xFFFFFFFF then none
      nodes := nodes.push none
    else
      let some g := gateOfBytes tag x y z | none
      if idx == 0 || parent.toNat ≥ idx || !validGate cfg width g then none
      nodes := nodes.push (some ⟨parent.toNat, g⟩)
    keys := keys.insert fp idx
    off := off + 16
  if off + 5 > b.size then none
  let saturated := b[off]! != 0
  let some depth := getU32 b (off + 1) | none
  some ({ keys, nodes, saturated, depth := depth.toNat }, off + 5)

/-- Deserialize a MURM, or `none` if the bytes do not match `cfg` exactly. -/
def deserialize (cfg : MurmConfig) (b : ByteArray) : Option Murm := do
  if b.size < 8 then none
  let payloadSize := b.size - 8
  let payload := b.extract 0 payloadSize
  let expected ← getU64 b payloadSize
  if checksum payload != expected then none
  let start ← readHeader payload cfg
  let widthCount ← getU32 payload start
  if widthCount.toNat != cfg.maxQubits + 1 then none
  let mut off := start + 4
  let mut widths : Array WidthMurm := Array.emptyWithCapacity widthCount.toNat
  for width in [0:widthCount.toNat] do
    let some (w, off') := getWidth payload off width cfg | none
    widths := widths.push w
    off := off'
  if off == payload.size then some { widths } else none

/-! ## The cache directory -/

/-- `~/.tzap-lean/murms`, if `HOME` is set. -/
def cacheDir : IO (Option System.FilePath) := do
  match ← IO.getEnv "HOME" with
  | none => return none
  | some home => return some ((System.FilePath.mk home) / ".tzap-lean" / "murms")

/-- One file per distinct bounds. The version is in the name as well as the header, so a bump
cannot collide with an old file. -/
def cachePath (cfg : MurmConfig) : IO (Option System.FilePath) := do
  match ← cacheDir with
  | none => return none
  | some dir =>
      return some (dir /
        s!"q{cfg.maxQubits}_g{cfg.maxGates}_e{cfg.maxEntriesPerQubit}_b{cfg.basis.mask}.v{formatVersion}.bin")

/-- Whether a usable cache file already exists — informational only; `load` re-validates. -/
def isCached (cfg : MurmConfig) : IO Bool := do
  match ← cachePath cfg with
  | none => return false
  | some p => p.pathExists

/-- Try to read a cached MURM. Any failure means "no usable cache", never a wrong MURM. -/
def load (cfg : MurmConfig) : IO (Option Murm) := do
  match ← cachePath cfg with
  | none => return none
  | some p =>
      try
        if !(← p.pathExists) then return none
        let bytes ← IO.FS.readBinFile p
        return deserialize cfg bytes
      catch _ => return none

/-- Write a MURM for later runs. Failure is ignored: a cache that cannot be written is a
missed speedup, not an error. Written to a temporary sibling and renamed, so a concurrent
reader never sees a partial file. -/
def store (cfg : MurmConfig) (murm : Murm) : IO Unit := do
  match ← cachePath cfg with
  | none => pure ()
  | some p =>
      try
        match ← cacheDir with
        | none => pure ()
        | some dir =>
            IO.FS.createDirAll dir
            let tmp := p.withExtension "tmp"
            IO.FS.writeBinFile tmp (serialize cfg murm)
            IO.FS.rename tmp p
      catch _ => pure ()

/-- The MURM for `cfg`, from disk if it is there and valid, otherwise built and then stored.
The `Bool` says whether it came from the cache. -/
def loadOrBuild (cfg : MurmConfig) : IO (Murm × Bool) := do
  match ← load cfg with
  | some murm => return (murm, true)
  | none =>
      let murm := buildMurm cfg
      store cfg murm
      return (murm, false)

end MurmCache

end TzapLean
