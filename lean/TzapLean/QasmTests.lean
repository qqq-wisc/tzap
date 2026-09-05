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

def metricsSample : Circuit := Circuit.ofGates 3 0
  [Gate.h 0, Gate.h 1, Gate.cnot 0 2, Gate.t 1, Gate.rz (1/4) 2, Gate.cz 1 2]

#guard Metrics.of metricsSample ==
  { gates := 6, twoQubit := 2, depth := 4, t := 1, rz := 1 }

/-- Parse and render compactly: qubit count, cbit count, gates. -/
def render (r : Except String Circuit) : String :=
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
def roundTrip (c : Circuit) : Bool :=
  match parse (serialize c) with
  | .error _ => false
  | .ok c' => c'.gates == c.gates && c'.numQubits == c.numQubits && c'.numCbits == c.numCbits

#guard roundTrip (Circuit.ofGates 3 1
  [.h 0, .cnot 0 1, .t 1, .ccx 0 1 2, .cz 0 1, .measure 0 0, .reset 1])

/-! ## Validation

`cx q[0],q[0]` parses but is not a circuit any pass may be handed: `Pass.correct` assumes
`Circuit.Wf`, and without it `CancelGates` would delete the pair as self-inverse — which
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
#guard roundTrip (Circuit.ofGates 3 0 [.x 0, .z 1, .s 2, .sdg 0, .tdg 1, .ccz 0 1 2])
#guard serialize (Circuit.ofGates 1 0 [.h 0]) ==
  "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0];\n"

/-! ## Options -/

#guard (PassName.parse "CancelGates").isSome
#guard (PassName.parse "SuperOpt").isSome
#guard (PassName.parse "DecomposeRz").isNone
#guard (PassName.parse "nonsense").isNone
#guard PassName.all.length == 4
-- Every executable pass carries an unconditional proof.
#guard (PassName.all.filter (fun p => !p.2.1.verified)).map (·.1) == []
-- `O1` is the only level that skips the table.
#guard !Level.O1.usesSuperOpt
#guard Level.O2.usesSuperOpt && Level.O3.usesSuperOpt && !Level.O1.usesSuperOpt
-- `O2` is the bounded tier.
#guard Level.O2.maxRounds == some 2
#guard Level.O3.maxRounds == none
-- `CnotMin` leads the superoptimizing sweep.
#guard Level.O3.pipeline.head? == some PassName.CnotMin
#guard Level.O1.pipeline == [PassName.CancelGates, PassName.PhaseFoldRand]

/-! ## Checked serialization -/

#guard (match Qasm.serializeChecked (Circuit.ofGates 2 0 [.h 0, .cnot 0 1, .t 1]) with
  | .ok _ => true
  | .error _ => false)

-- This build deliberately cannot parse `rz`, so the checked writer fails closed.
#guard (match Qasm.serializeChecked (Circuit.ofGates 1 0 [.rz (1/3) 0]) with
  | .ok _ => false
  | .error _ => true)

/-! ## Metrics and formatting -/

#guard (Metrics.of (Circuit.ofGates 2 0 [.t 0, .tdg 1, .cnot 0 1, .h 0])).t == 2
#guard (Metrics.of (Circuit.ofGates 2 0 [.t 0, .tdg 1, .cnot 0 1, .h 0])).twoQubit == 1
#guard (Metrics.of (Circuit.ofGates 2 0 [.t 0, .tdg 1, .cnot 0 1, .h 0])).gates == 4
#guard fmtNum 0 == "0"
#guard fmtNum 999 == "999"
#guard fmtNum 1000 == "1,000"
#guard fmtNum 1234567 == "1,234,567"
#guard fmtPct 100 25 == "75.0"
#guard fmtPct 0 0 == "0.0"
#guard fmtPct 3 3 == "0.0"

/-! ## The table cache

Round-trip through the on-disk format, and the rejections that make a bad file harmless. -/

/-- A small table to serialize. -/
def cacheCfg : SuperOptTableConfig := { maxQubits := 2, maxGates := 2, maxEntriesPerQubit := 500 }

/-- Its bytes. -/
def cacheBytes : ByteArray := TableCache.serialize cacheCfg (buildTable cacheCfg)

/-- Reading back what was written gives the same widths, entries and circuits. -/
def cacheRoundTrips : Bool :=
  let orig := buildTable cacheCfg
  match TableCache.deserialize cacheCfg cacheBytes with
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
#guard (TableCache.deserialize { cacheCfg with maxQubits := 3 } cacheBytes).isNone
#guard (TableCache.deserialize { cacheCfg with maxGates := 3 } cacheBytes).isNone
#guard (TableCache.deserialize { cacheCfg with maxEntriesPerQubit := 501 } cacheBytes).isNone
-- Garbage, a truncated write, and a wrong magic are all rejected.
#guard (TableCache.deserialize cacheCfg (ByteArray.mk #[71, 65, 82, 66])).isNone
#guard (TableCache.deserialize cacheCfg (cacheBytes.extract 0 100)).isNone
#guard (TableCache.deserialize cacheCfg ByteArray.empty).isNone

-- A synthesized lookup survives the round trip: the cached table answers as the built one does.
#guard
  (let orig := buildTable cacheCfg
   match TableCache.deserialize cacheCfg cacheBytes with
   | none => false
   | some back =>
       [[Gate.h 0, Gate.h 0], [Gate.t 0, Gate.t 0], [Gate.x 0], [Gate.cnot 0 1]].all
         fun gs =>
           match ExactMat.matrixOf 2 gs with
           | none => false
           | some M => orig.synthesize 2 M.normalize == back.synthesize 2 M.normalize)

end TzapLean
