import TzapLean.Qasm
import TzapLean.Optimize

/-!
# QASM and CLI tests

The parser half of `src/qasm.rs`'s suite, plus the option plumbing. Every error message is
checked against the Rust text it is ported from, since a parser's error messages are most of
its user interface.
-/

namespace TzapLean

open Qasm

/-! ## Driver metrics -/

def metricsSample : RawCircuit := RawCircuit.ofGates 3 0
  [Gate.h 0, Gate.h 1, Gate.cnot 0 2, Gate.t 1, Gate.rz (1/4) 2, Gate.cz 1 2]

#guard Metrics.of metricsSample ==
  { gates := 6, twoQubit := 2, depth := 4, t := 1, rz := 1 }

/-- Parse and render compactly: qubit count, cbit count, gates. -/
def render (r : Except String RawCircuit) : String :=
  match r with
  | .error e => s!"ERR {e}"
  | .ok c =>
      s!"{c.numQubits}q {c.numCbits}c [{String.intercalate ", " (c.gates.map Gate.toString)}]"

/-! ## Gates -/

#guard render (parse "qreg q[1];\nz q[0];") == "1q 0c [z q0]"
#guard render (parse "qreg q[1];\nsdg q[0];") == "1q 0c [sdg q0]"
#guard render (parse "qreg q[1];\ntdg q[0];\nt q[0];") == "1q 0c [tdg q0, t q0]"
#guard render (parse "qreg q[1];\nsdg q[0];\ns q[0];") == "1q 0c [sdg q0, s q0]"
#guard render (parse "qreg q[3];\nh q[0];\ncx q[0],q[1];\nccx q[0],q[1],q[2];\ncz q[1],q[2];")
  == "3q 0c [h q0, cnot q0, q1, ccx q0, q1, q2, cz q1, q2]"
#guard render (parse "qreg q[3];\nccz q[0],q[1],q[2];") == "3q 0c [ccz q0, q1, q2]"

/-! ## Registers -/

-- Two registers share one flat index space, in declaration order.
#guard render (parse "qreg a[2];\nqreg b[2];\ncx a[1],b[0];") == "4q 0c [cnot q1, q2]"
-- `measure q[i] -> c[j];` and the register-broadcast form `measure q -> c;`.
#guard render (parse "qreg q[2];\ncreg c[2];\nmeasure q[0] -> c[0];")
  == "2q 2c [measure q0 -> c0]"
#guard render (parse "qreg q[2];\ncreg c[2];\nmeasure q -> c;")
  == "2q 2c [measure q0 -> c0, measure q1 -> c1]"
#guard render (parse "qreg q[3];\nreset q;") == "3q 0c [reset q0, reset q1, reset q2]"

/-! ## Lexical -/

#guard render (parse "qreg q[1];\nh q[0]; // trailing") == "1q 0c [h q0]"
#guard render (parse "qreg q[1];\n/* skip\nme */\nh q[0];") == "1q 0c [h q0]"
#guard render (parse "qreg q[1]; h q[0]; t q[0];") == "1q 0c [h q0, t q0]"
#guard render (parse "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nbarrier q;\nh q[0];")
  == "1q 0c [h q0]"

/-! ## Errors

Each message is the Rust one, with `rz` the single deliberate difference. -/

#guard render (parse "qreg q[1];\nfoo q[0];") == "ERR line 2: unsupported: foo q[0]"
#guard render (parse "qreg q[1];\nh q[0];\nqreg r[1];") == "ERR line 3: qreg declaration after gate"
#guard render (parse "qreg q[1];\nh q[5];")
  == "ERR line 2: index 5 out of range for register 'q' (size 1)"
#guard render (parse "qreg q[1];\nh r[0];") == "ERR line 2: unknown register 'r'"
#guard render (parse "qreg q[2];\ncx q[0];") == "ERR line 2: cx expects 2 qubit operands, got 1"
#guard render (parse "qreg q[2];\ncreg c[1];\nmeasure q -> c;")
  == "ERR line 3: measure operand size mismatch (2 qubits, 1 cbits)"
#guard render (parse "qreg q[1];\ncreg c[1];\nmeasure q[0] c[0];")
  == "ERR line 3: measure missing '->' (got 'q[0] c[0]')"
-- `rz` is rejected rather than silently dropped.
#guard (parse "qreg q[1];\nrz(pi/4) q[0];").isOk == false

/-! ## Serializing -/

/-- A circuit, its QASM, and the circuit that QASM parses back to. -/
def roundTrip (c : RawCircuit) : Bool :=
  match parse (serialize c) with
  | .error _ => false
  | .ok c' => c'.gates == c.gates && c'.numQubits == c.numQubits && c'.numCbits == c.numCbits

#guard roundTrip (RawCircuit.ofGates 3 1
  [.h 0, .cnot 0 1, .t 1, .ccx 0 1 2, .cz 0 1, .measure 0 0, .reset 1])

/-! ## Validation

`cx q[0],q[0]` parses but cannot be packaged as the checked circuit a pass requires. Without
`RawCircuit.Wf`, `CancelGates` would delete the pair as self-inverse — which
`cnot q q` is not. The front end rejects it rather than passing it on. -/

/-- Whether the parser accepts a source. -/
def parseAccepts (src : String) : Bool := (parse src).toOption.isSome

#guard !parseAccepts "OPENQASM 2.0;\nqreg q[1];\ncx q[0],q[0];\n"
#guard !parseAccepts "OPENQASM 2.0;\nqreg q[2];\ncz q[1],q[1];\n"
#guard !parseAccepts "OPENQASM 2.0;\nqreg q[3];\nccx q[0],q[1],q[1];\n"
#guard !parseAccepts "OPENQASM 2.0;\nqreg q[3];\nccz q[2],q[2],q[0];\n"
#guard parseAccepts "OPENQASM 2.0;\nqreg q[2];\ncx q[0],q[1];\n"

-- The offending gate is named in the message.
#guard (match parse "OPENQASM 2.0;\nqreg q[1];\ncx q[0],q[0];\n" with
        | .error e => (e.splitOn "distinct").length == 2
        | .ok _ => false)
#guard roundTrip (RawCircuit.ofGates 3 0 [.x 0, .z 1, .s 2, .sdg 0, .tdg 1, .ccz 0 1 2])
#guard serialize (RawCircuit.ofGates 1 0 [.h 0]) ==
  "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0];\n"

/-! ## Options -/

#guard (PassName.parse "CancelGates").isSome
#guard (PassName.parse "SuperOpt").isSome
#guard (PassName.parse "DecomposeRz").isNone
#guard (PassName.parse "nonsense").isNone
#guard PassName.all.length == 6
#guard tagBits == 128
-- Phase folding is probabilistically bounded; the other executable passes are unconditional.
#guard (PassName.all.filter (fun p => !p.2.1.verified)).map (·.1) == ["PhaseFoldRand"]
-- `O1` is the only level that skips the MURM.
#guard !Level.O1.usesSuperOpt
#guard Level.O2.usesSuperOpt && Level.O3.usesSuperOpt && !Level.O1.usesSuperOpt
-- `O2` is the bounded tier.
#guard Level.O2.maxRounds == some 2
#guard Level.O3.maxRounds == none
-- `CnotMin` leads the superoptimizing sweep.
#guard Level.O3.pipeline.head? == some PassName.CnotMin
#guard Level.O1.pipeline == [PassName.CancelGates, PassName.PhaseFoldRand]

/-! ## Checked serialization -/

#guard (match Qasm.serializeChecked (RawCircuit.ofGates 2 0 [.h 0, .cnot 0 1, .t 1]) with
  | .ok _ => true
  | .error _ => false)

-- This build deliberately cannot parse `rz`, so the checked writer fails closed.
#guard (match Qasm.serializeChecked (RawCircuit.ofGates 1 0 [.rz (1/3) 0]) with
  | .ok _ => false
  | .error _ => true)

/-! ## Metrics and formatting -/

#guard (Metrics.of (RawCircuit.ofGates 2 0 [.t 0, .tdg 1, .cnot 0 1, .h 0])).t == 2
#guard (Metrics.of (RawCircuit.ofGates 2 0 [.t 0, .tdg 1, .cnot 0 1, .h 0])).twoQubit == 1
#guard (Metrics.of (RawCircuit.ofGates 2 0 [.t 0, .tdg 1, .cnot 0 1, .h 0])).gates == 4
#guard fmtNum 0 == "0"
#guard fmtNum 999 == "999"
#guard fmtNum 1000 == "1,000"
#guard fmtNum 1234567 == "1,234,567"
#guard fmtPct 100 25 == "75.0"
#guard fmtPct 0 0 == "0.0"
#guard fmtPct 3 3 == "0.0"

/-! ## SuperOpt synthesis-basis selection -/

#guard SuperOptGates.parse "auto" == .ok .auto
#guard SuperOptGates.parse "base" == .ok .base
#guard SuperOptGates.parse "h,x,cz,ccx,ccz" ==
  .ok (.explicit (GateSet.ofKinds [.h, .x, .cz, .ccx, .ccz]))
#guard SuperOptGates.parse "h,h,cz" ==
  .ok (.explicit (GateSet.ofKinds [.h, .cz]))
#guard (match SuperOptGates.parse "rz" with | .error _ => true | .ok _ => false)
#guard (match SuperOptGates.parse "" with | .error _ => true | .ok _ => false)
#guard (.auto : SuperOptGates).effective (GateSet.ofKinds [.rz, .cz, .ccz]) ==
  baseGateSet.union (GateSet.ofKinds [.cz, .ccz])
#guard (.base : SuperOptGates).effective optionalGateSet == baseGateSet

/-! ## Native-gate staging matrix

The eight optional-native input gate sets are crossed with all four decomposition-flag
combinations.  This checks both lowering and the basis seen by the post-decomposition MURM. -/

def optionalCircuit (mask : Nat) : RawCircuit :=
  RawCircuit.ofGates 3 0 <|
    (if mask.testBit 0 then [.cz 0 1] else []) ++
    (if mask.testBit 1 then [.ccx 0 1 2] else []) ++
    (if mask.testBit 2 then [.ccz 0 1 2] else [])

def lowerRequested (c : RawCircuit) (decomposeCcx decomposeCz : Bool) : RawCircuit :=
  let c := if decomposeCcx then runPassRaw DecomposeToffoli c else c
  if decomposeCz then runPassRaw DecomposeCz c else c

def stagingCaseOk (mask : Nat) (decomposeCcx decomposeCz : Bool) : Bool :=
  let input := optionalCircuit mask
  let options : Options := { decomposeCcx, decomposeCz, superoptGates := .auto }
  let lowered := lowerRequested input decomposeCcx decomposeCz
  let (_, postCfg) := resolveBounds options lowered.gateSet options.decomposedGateSet
  (!decomposeCcx ||
      (!lowered.gateSet.contains .ccx && !lowered.gateSet.contains .ccz &&
       !postCfg.basis.contains .ccx && !postCfg.basis.contains .ccz)) &&
    (!decomposeCz || (!lowered.gateSet.contains .cz && !postCfg.basis.contains .cz)) &&
    (decomposeCcx ||
      (lowered.gateSet.contains .ccx == input.gateSet.contains .ccx &&
       lowered.gateSet.contains .ccz == input.gateSet.contains .ccz)) &&
    (decomposeCz || lowered.gateSet.contains .cz == input.gateSet.contains .cz)

-- `8 gate subsets × 2 CCX choices × 2 CZ choices = 32` staging configurations.
#guard (List.range 8).all fun mask =>
  [false, true].all fun decomposeCcx =>
    [false, true].all fun decomposeCz =>
      stagingCaseOk mask decomposeCcx decomposeCz

-- Auto observes the circuit entering each stage, rather than retaining the input basis.
#guard
  let input := optionalCircuit 7
  let lowered := lowerRequested input true false
  let options : Options := { decomposeCcx := true }
  let (_, before) := resolveBounds options input.gateSet
  let (_, after) := resolveBounds options lowered.gateSet options.decomposedGateSet
  before.basis == baseGateSet.union optionalGateSet &&
    after.basis == baseGateSet.insert .cz

-- Explicit mode remains exact, except that a requested lowering cannot be synthesized back.
#guard
  let requested := GateSet.ofKinds [.h, .cz, .ccx, .ccz]
  let options : Options :=
    { decomposeCcx := true, decomposeCz := true,
      superoptGates := .explicit requested }
  let (_, after) := resolveBounds options (lowerRequested (optionalCircuit 7) true true).gateSet
    options.decomposedGateSet
  after.basis == GateSet.singleton .h

/-! ## Complete basis-selection matrices -/

def optionalSetFromMask (mask : Nat) : GateSet :=
  GateSet.ofKinds <| [.cz, .ccx, .ccz].zipIdx.filterMap fun (kind, i) =>
    if mask.testBit i then some kind else none

/-- The ten documented synthesis modes: auto, base, and base plus each optional-native mask. -/
def matrixBasisMode (mode : Nat) : SuperOptGates :=
  if mode == 0 then .auto
  else if mode == 1 then .base
  else .explicit (baseGateSet.union (optionalSetFromMask (mode - 2)))

def basisSelectionCaseOk (inputMask mode : Nat) : Bool :=
  let inputSet := (optionalCircuit inputMask).gateSet
  let actual := (matrixBasisMode mode).effective inputSet
  let expected :=
    if mode == 0 then baseGateSet.union (optionalSetFromMask inputMask)
    else if mode == 1 then baseGateSet
    else baseGateSet.union (optionalSetFromMask (mode - 2))
  actual == expected

-- `8 input subsets × 10 basis modes = 80` exact basis-selection cases.
#guard (List.range 8).all fun inputMask =>
  (List.range 10).all fun mode => basisSelectionCaseOk inputMask mode

def denseGateSet (kinds : List GateKind) (mask : Nat) : GateSet :=
  GateSet.ofKinds <| kinds.zipIdx.filterMap fun (kind, i) =>
    if mask.testBit i then some kind else none

def explicitBasisRoundTrip (mask : Nat) : Bool :=
  let kinds := supportedGateSet.kinds
  let gates := denseGateSet kinds mask
  let text := String.intercalate "," (gates.kinds.map GateKind.qasmName)
  SuperOptGates.parse text == .ok (.explicit gates)

-- Every one of the 2,047 nonempty subsets of the eleven synthesis gates parses exactly.
#guard (List.range (2 ^ supportedGateSet.kinds.length - 1)).all fun i =>
  explicitBasisRoundTrip (i + 1)

def matrixLevel : Nat → Level
  | 0 => .O1
  | 1 => .O2
  | _ => .O3

def crossFeatureCaseOk (inputMask flags seed level : Nat) : Bool :=
  let decomposeCcx := flags.testBit 0
  let decomposeCz := flags.testBit 1
  let repeated :=
    (List.range (seed % 4 + 1)).flatMap fun _ => (optionalCircuit inputMask).gates
  let input := RawCircuit.ofGates 3 0 repeated
  let options : Options :=
    { level := matrixLevel level, decomposeCcx, decomposeCz,
      superoptGates := matrixBasisMode ((seed + inputMask + flags) % 10) }
  let lowered := lowerRequested input decomposeCcx decomposeCz
  let cancelled := runPassRaw CancelGates lowered
  let (_, postCfg) := resolveBounds options cancelled.gateSet options.decomposedGateSet
  cancelled.gates.all Gate.Wf && postCfg.basis.isSubset supportedGateSet &&
    (options.decomposedGateSet.inter postCfg.basis).isEmpty &&
    (!decomposeCcx ||
      (!cancelled.gateSet.contains .ccx && !cancelled.gateSet.contains .ccz)) &&
    (!decomposeCz || !cancelled.gateSet.contains .cz)

-- `8 inputs × 4 decomposition choices × 8 seeds × 3 levels = 768` cross-feature cases.
#guard (List.range 8).all fun inputMask =>
  (List.range 4).all fun flags =>
    (List.range 8).all fun seed =>
      (List.range 3).all fun level => crossFeatureCaseOk inputMask flags seed level

/-! ## The MURM cache

Round-trip through the on-disk format, and the rejections that make a bad file harmless. -/

/-- A small table to serialize. -/
def cacheCfg : MurmConfig := { maxQubits := 2, maxGates := 2, maxEntriesPerQubit := 500 }

/-- Its bytes. -/
def cacheBytes : ByteArray := MurmCache.serialize cacheCfg (buildMurm cacheCfg)

/-- Reading back what was written gives the same widths, entries and circuits. -/
def cacheRoundTrips : Bool :=
  let orig := buildMurm cacheCfg
  match MurmCache.deserialize cacheCfg cacheBytes with
  | none => false
  | some back =>
      back.widths.size == orig.widths.size &&
        (List.range orig.widths.size).all fun k =>
          let a := orig.widths[k]!
          let b := back.widths[k]!
          a.size == b.size && a.nodes.size == b.nodes.size &&
            a.saturated == b.saturated && a.depth == b.depth &&
            -- every stored circuit survives the trip
            (List.range a.nodes.size).all fun i => a.circuitOf i == b.circuitOf i

#guard cacheRoundTrips

-- A cache built for other bounds is rejected, not misread.
#guard (MurmCache.deserialize { cacheCfg with maxQubits := 3 } cacheBytes).isNone
#guard (MurmCache.deserialize { cacheCfg with maxGates := 3 } cacheBytes).isNone
#guard (MurmCache.deserialize { cacheCfg with maxEntriesPerQubit := 501 } cacheBytes).isNone
#guard (MurmCache.deserialize { cacheCfg with basis := GateSet.singleton .cz } cacheBytes).isNone
-- Garbage, a truncated write, and a wrong magic are all rejected.
#guard (MurmCache.deserialize cacheCfg (ByteArray.mk #[71, 65, 82, 66])).isNone
#guard (MurmCache.deserialize cacheCfg (cacheBytes.extract 0 100)).isNone
#guard (MurmCache.deserialize cacheCfg ByteArray.empty).isNone
#guard (MurmCache.deserialize cacheCfg (cacheBytes.push 0)).isNone
#guard (MurmCache.deserialize cacheCfg
  (cacheBytes.set! 10 (cacheBytes[10]! ^^^ 1))).isNone

-- A synthesized lookup survives the round trip: the cached table answers as the built one does.
#guard
  (let orig := buildMurm cacheCfg
   match MurmCache.deserialize cacheCfg cacheBytes with
   | none => false
   | some back =>
       [[Gate.h 0, Gate.h 0], [Gate.t 0, Gate.t 0], [Gate.x 0], [Gate.cnot 0 1]].all
         fun gs =>
           match ExactMat.matrixOf 2 gs with
           | none => false
           | some M => orig.synthesize 2 M.normalize == back.synthesize 2 M.normalize)

end TzapLean
