/-!
# Packed arithmetic for the PhaseFoldRand fingerprint field

This module is the executable half of Rust's `GF(2^128)` representation.  Values are packed
into the low 128 bits of a `Nat`; addition is XOR and multiplication is carryless polynomial
multiplication reduced by

`z^128 + z^7 + z^2 + z + 1`.

The reduction tail is therefore `0x87`. The modulus is proved irreducible in
`GF128Certificate.lean`, and `GF128Field.lean` constructs the abstract quotient field.
`GF128Bridge.lean` proves that these operations represent addition and multiplication in the
quotient field, then transfers the nonlinear collision bound to the concrete sample format.
-/

namespace TzapLean

namespace GF128

/-- Number of bits in one fingerprint value. -/
def width : Nat := 128

/-- Low-bit mask for the packed representation. -/
def mask : Nat := 2 ^ width - 1

/-- The non-leading coefficients of `z^128 + z^7 + z^2 + z + 1`. -/
def reduction : Nat := 0x87

/-- Canonicalize a natural number to its low 128 bits. -/
def normalize (x : Nat) : Nat := x &&& mask

/-- Field addition in characteristic two. -/
def add (x y : Nat) : Nat := normalize x ^^^ normalize y

/-- Multiply by `z`, reducing the leading term when bit 127 was set. -/
def xtime (x : Nat) : Nat :=
  let x := normalize x
  normalize (x <<< 1) ^^^ (if x.testBit 127 then reduction else 0)

/-- Multiply-and-reduce loop, exposed for its inductive refinement proof. -/
def mulAux : Nat → Nat → Nat → Nat → Nat
  | 0, _, _, acc => normalize acc
  | fuel + 1, a, b, acc =>
      mulAux fuel (xtime a) (b >>> 1) (if b.testBit 0 then acc ^^^ a else acc)

/-- Carryless multiply-and-reduce. Exactly 128 iterations are used, matching the fixed
representation width and ignoring bits outside the canonical input representation. -/
def mul (x y : Nat) : Nat :=
  mulAux width (normalize x) (normalize y) 0

@[simp] theorem normalize_zero : normalize 0 = 0 := by native_decide
@[simp] theorem add_zero (x : Nat) : add x 0 = normalize x := by simp [add]

end GF128

/-! ## Degree metadata and the conservative cutoff

The degree is an upper bound for the formal polynomial represented by `value`, not part of
the field value itself.  When a CCX update would exceed `2^32`, the caller must replace the
target with a fresh opaque fingerprint instead of multiplying. -/

/-- A packed random-interpretation value and an upper bound on its formal total degree. -/
structure Fingerprint where
  value : Nat
  degree : Nat
deriving Repr, Inhabited, DecidableEq

namespace Fingerprint

/-- Maximum degree retained by nonlinear phase folding. -/
def maxDegree : Nat := 2 ^ 32

def zero : Fingerprint := ⟨0, 0⟩
def one : Fingerprint := ⟨1, 0⟩
def fresh (value : Nat) : Fingerprint := ⟨GF128.normalize value, 1⟩

/-- XOR preserves the larger of the two degree bounds. -/
def add (x y : Fingerprint) : Fingerprint :=
  ⟨GF128.add x.value y.value, max x.degree y.degree⟩

/-- The prospective product degree.  Callers inspect this before constructing the product. -/
def productDegree (x y : Fingerprint) : Nat := x.degree + y.degree

/-- Multiply fingerprints when the resulting degree remains within the proof budget. -/
def mul? (x y : Fingerprint) : Option Fingerprint :=
  let degree := productDegree x y
  if degree ≤ maxDegree then some ⟨GF128.mul x.value y.value, degree⟩ else none

/-- Nonlinear CCX update.  A fresh draw is used precisely when the product would exceed the
cutoff; otherwise the target is XORed with the controls' product. -/
def ccx (target left right : Fingerprint) (freshValue : Nat) : Fingerprint :=
  match mul? left right with
  | some product => add target product
  | none => fresh freshValue

end Fingerprint

end TzapLean
