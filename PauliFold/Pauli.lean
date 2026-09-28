import PauliFold.Circuit

/-!
# The Pauli-fold abstract domain

The abstraction used by `PauliFoldRand` (Section 8.1 of the notes): a candidate rotation's
transported axis is either a single *signed Pauli string* or unknown (`⊤`).

* Cliffords act exactly, as signed permutations of Pauli strings. The tables below are the
  forward conjugations `U P U†`; they were generated numerically from the gate matrices.
* A T or T† on qubit `q` rotates about `Z_q`. It leaves an axis that commutes with `Z_q`
  unchanged, and sends one that anticommutes (letter X or Y at `q`) to `⊤`.
-/

namespace PauliFold

open TzapLean

/-- A single-qubit Pauli letter. -/
inductive Letter where
  | I | X | Y | Z
  deriving DecidableEq, Repr

/-- A signed Pauli string: `neg` is the sign (`true` for `-1`), `letter q` the letter at
qubit `q`. -/
structure SPauli (n : Nat) where
  neg : Bool
  letter : Fin n → Letter

namespace Letter

/-- Matrix of a letter, rows as outputs and columns as inputs (the convention of
`TzapLean.embed1`). -/
def mat : Letter → Bool → Bool → ℂ
  | .I => fun o i => if o = i then 1 else 0
  | .X => fun o i => if o = !i then 1 else 0
  | .Y => fun o i => if o = !i then (if o then Complex.I else -Complex.I) else 0
  | .Z => fun o i => if o = i then (if o then -1 else 1) else 0

/-! Forward conjugation tables `U L U† = ± L'`: the sign flag is `true` for `-`. -/

def conjH : Letter → Bool × Letter
  | .I => (false, .I) | .X => (false, .Z) | .Y => (true, .Y) | .Z => (false, .X)

def conjS : Letter → Bool × Letter
  | .I => (false, .I) | .X => (false, .Y) | .Y => (true, .X) | .Z => (false, .Z)

def conjSdg : Letter → Bool × Letter
  | .I => (false, .I) | .X => (true, .Y) | .Y => (false, .X) | .Z => (false, .Z)

def conjX : Letter → Bool × Letter
  | .I => (false, .I) | .X => (false, .X) | .Y => (true, .Y) | .Z => (true, .Z)

def conjZ : Letter → Bool × Letter
  | .I => (false, .I) | .X => (true, .X) | .Y => (true, .Y) | .Z => (false, .Z)

/-- CX with the control's letter first. -/
def conjCX : Letter → Letter → Bool × Letter × Letter
  | .I, .I => (false, .I, .I) | .I, .X => (false, .I, .X)
  | .I, .Y => (false, .Z, .Y) | .I, .Z => (false, .Z, .Z)
  | .X, .I => (false, .X, .X) | .X, .X => (false, .X, .I)
  | .X, .Y => (false, .Y, .Z) | .X, .Z => (true, .Y, .Y)
  | .Y, .I => (false, .Y, .X) | .Y, .X => (false, .Y, .I)
  | .Y, .Y => (true, .X, .Z) | .Y, .Z => (false, .X, .Y)
  | .Z, .I => (false, .Z, .I) | .Z, .X => (false, .Z, .X)
  | .Z, .Y => (false, .I, .Y) | .Z, .Z => (false, .I, .Z)

/-- CZ with the control's letter first. -/
def conjCZ : Letter → Letter → Bool × Letter × Letter
  | .I, .I => (false, .I, .I) | .I, .X => (false, .Z, .X)
  | .I, .Y => (false, .Z, .Y) | .I, .Z => (false, .I, .Z)
  | .X, .I => (false, .X, .Z) | .X, .X => (false, .Y, .Y)
  | .X, .Y => (true, .Y, .X) | .X, .Z => (false, .X, .I)
  | .Y, .I => (false, .Y, .Z) | .Y, .X => (true, .X, .Y)
  | .Y, .Y => (false, .X, .X) | .Y, .Z => (false, .Y, .I)
  | .Z, .I => (false, .Z, .I) | .Z, .X => (false, .I, .X)
  | .Z, .Y => (false, .I, .Y) | .Z, .Z => (false, .Z, .Z)

/-- Whether the letter anticommutes with `Z`. -/
def anticommutesZ : Letter → Bool
  | .X | .Y => true
  | .I | .Z => false

end Letter

namespace SPauli

variable {n : Nat}

/-- The sign `±1` as a complex number. -/
def sgn (b : Bool) : ℂ := if b then -1 else 1

/-- The operator `±P` as a matrix: the tensor product of its letters, entry by entry. -/
noncomputable def toMatrix (P : SPauli n) : Op n :=
  fun o i => sgn P.neg * ∏ r, (P.letter r).mat (o r) (i r)

/-- Apply a single-qubit conjugation table at qubit `q`. -/
def apply1 (f : Letter → Bool × Letter) (q : Fin n) (P : SPauli n) : SPauli n :=
  ⟨xor P.neg (f (P.letter q)).1, Function.update P.letter q (f (P.letter q)).2⟩

/-- Apply a two-qubit conjugation table at qubits `c` and `t`. -/
def apply2 (f : Letter → Letter → Bool × Letter × Letter) (c t : Fin n) (P : SPauli n) :
    SPauli n :=
  let r := f (P.letter c) (P.letter t)
  ⟨xor P.neg r.1, Function.update (Function.update P.letter c r.2.1) t r.2.2⟩

/-- The exact action of a Clifford gate; the rotations T and T† are handled by the abstract
transformer, since they do not map Pauli strings to Pauli strings. -/
def conjClifford : CTGate n → SPauli n → SPauli n
  | .h q => apply1 Letter.conjH q
  | .s q => apply1 Letter.conjS q
  | .sdg q => apply1 Letter.conjSdg q
  | .x q => apply1 Letter.conjX q
  | .z q => apply1 Letter.conjZ q
  | .cx ctl tgt _ => apply2 Letter.conjCX ctl tgt
  | .cz ctl tgt _ => apply2 Letter.conjCZ ctl tgt
  | .t _ | .tdg _ => id

end SPauli

/-- The Pauli-fold abstract domain: one signed Pauli string, or `⊤`. -/
inductive PauliAbs (n : Nat) where
  | one (P : SPauli n)
  | top

namespace PauliAbs

variable {n : Nat}

/-- The abstract transformer of a gate. -/
def step : CTGate n → PauliAbs n → PauliAbs n
  | _, .top => .top
  | .t q, .one P | .tdg q, .one P => if (P.letter q).anticommutesZ then .top else .one P
  | g, .one P => .one (P.conjClifford g)

/-- The abstract transformer of a gate list. -/
def run (gs : List (CTGate n)) (a : PauliAbs n) : PauliAbs n :=
  gs.foldl (fun a g => step g a) a

/-- Concretization: a single signed Pauli denotes its matrix; `⊤` denotes every operator. -/
noncomputable def gamma : PauliAbs n → Set (Op n)
  | .one P => {P.toMatrix}
  | .top => Set.univ

end PauliAbs

end PauliFold
