import PauliFold.Pauli
import Mathlib.Data.Finsupp.Basic
import Mathlib.Data.ZMod.Defs
import Mathlib.Algebra.BigOperators.Pi

/-!
# The symbolic phase-folding abstract domain

The exact counterpart of `PhaseFoldRand` (Section 10.1 of the notes). Path variables are
natural numbers: `0, …, n-1` are the inputs, and each H introduces a fresh one. A wire
carries an *affine* Boolean form, a linear part (a finitely supported `ℕ → ZMod 2`) plus a
constant.

The transformers are written as additions, so that sums over wires are easy to track:

* X on `q` adds the constant 1 to wire `q`;
* CX from `c` to `t` adds wire `c` to wire `t`;
* H on `q` replaces wire `q` by the fresh variable (adds `fresh - w q`);
* Z, S, S†, T, T†, and CZ are diagonal and change nothing.

The state also records the *label*: the affine form of the candidate rotation, fixed at the
start of the segment. Its concretization is the operator the label describes: if the
label's linear part is the sum of the wires in a set `Q`, the transported operator is the
Z-string on `Q`, with the sign given by the constants.
-/

namespace PauliFold

/-- Linear Boolean forms over the path variables. -/
abbrev Lin := ℕ →₀ ZMod 2

/-- Affine Boolean forms: a linear part and a constant. -/
abbrev Affine := Lin × ZMod 2

/-- The phase-folding abstract state. -/
structure PhaseAbs (n : Nat) where
  wires : Fin n → Affine
  fresh : ℕ
  label : Affine

namespace PhaseAbs

variable {n : Nat}

/-- The abstract transformer of a gate. -/
noncomputable def step : CTGate n → PhaseAbs n → PhaseAbs n
  | .x q, a => { a with wires := a.wires + Pi.single q (0, 1) }
  | .cx ctl tgt _, a => { a with wires := a.wires + Pi.single tgt (a.wires ctl) }
  | .h q, a =>
      { a with
        wires := a.wires + Pi.single q ((Finsupp.single a.fresh 1, 0) - a.wires q)
        fresh := a.fresh + 1 }
  | .z _, a | .s _, a | .sdg _, a | .t _, a | .tdg _, a | .cz _ _ _, a => a

/-- The abstract transformer of a gate list. -/
noncomputable def run (gs : List (CTGate n)) (a : PhaseAbs n) : PhaseAbs n :=
  gs.foldl (fun a g => step g a) a

/-- The sum of the wires in `Q`. -/
noncomputable def wireSum (a : PhaseAbs n) (Q : Finset (Fin n)) : Affine := ∑ i ∈ Q, a.wires i

/-- `Q` expresses the label: its wires sum to the label's linear part. -/
def Expresses (a : PhaseAbs n) (Q : Finset (Fin n)) : Prop :=
  (a.wireSum Q).1 = a.label.1

/-- The sign of the Z-string on `Q`: the constants of the wires and of the label. -/
noncomputable def sign (a : PhaseAbs n) (Q : Finset (Fin n)) : Bool :=
  decide ((a.wireSum Q).2 + a.label.2 = 1)

end PhaseAbs

/-- The Z-string on `Q` with sign `neg`. -/
def zString {n : Nat} (Q : Finset (Fin n)) (neg : Bool) : SPauli n :=
  ⟨neg, fun i => if i ∈ Q then .Z else .I⟩

namespace PhaseAbs

variable {n : Nat}

open Classical in
/-- Concretization: if some wire set expresses the label, the operators `±Z_Q` for every
such `Q`; otherwise every operator. -/
noncomputable def gamma (a : PhaseAbs n) : Set (Op n) :=
  if ∃ Q, a.Expresses Q then
    {M | ∃ Q, a.Expresses Q ∧ M = (zString Q (a.sign Q)).toMatrix}
  else Set.univ

/-- The start of a segment whose first rotation is on qubit `q`: each wire carries its
input variable, and the label is wire `q`'s value. -/
noncomputable def init (n : Nat) (q : Fin n) : PhaseAbs n where
  wires := fun i => (Finsupp.single i.val 1, 0)
  fresh := n
  label := (Finsupp.single q.val 1, 0)

end PhaseAbs

end PauliFold
