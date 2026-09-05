import TzapLean.GateAlgebra

/-!
# The Cyclotomic Ring `ℤ[ω]`, `ω = e^{iπ/4}`

`SuperOpt` decides whether two circuit fragments are the same unitary, so it cannot use
floating point: two Clifford+T circuits are equal or they are not, and a rounding error either
way is a wrong answer. Every Clifford+T matrix entry is

```
    (a + b·ω + c·ω² + d·ω³) / √2 ᵏ ,      ω = e^{iπ/4},  ω⁴ = -1
```

with integer `a, b, c, d` and one denominator exponent shared by the whole matrix — the ring
`ℤ[1/√2, i]`. This file is that ring: exact arithmetic, the eight global phases `ω^p`, and
division by `√2 = ω - ω³`, each with its interpretation in `ℂ`.

This is `src/super_opt/matrix.rs`'s `Cyclotomic`, with two differences. Rust packs the
coefficients into four `i8`s and reports overflow, which cannot happen here; and where Rust
*checks* divisibility by `√2` with a parity test, `divSqrt2_spec` proves that test right.
-/

namespace TzapLean

/-- `ω = e^{iπ/4}`, the primitive eighth root of unity Clifford+T lives in. -/
noncomputable def ω : ℂ := ep (1/4)

@[simp] theorem ω_pow_two : ω ^ 2 = Complex.I := by
  rw [ω, pow_two, ← ep_add, show (1/4 : ℚ) + 1/4 = 1/2 by norm_num, ep_half]

@[simp] theorem ω_pow_four : ω ^ 4 = -1 := by
  rw [show (4 : Nat) = 2 * 2 from rfl, pow_mul, ω_pow_two]
  simp

theorem ω_ne_zero : ω ≠ 0 := ep_ne_zero _

/-- `√2 = ω - ω³`: the reason a `√2` denominator stays inside the ring. -/
theorem sqrt_two_eq : ((Real.sqrt 2 : ℝ) : ℂ) = ω - ω ^ 3 := by
  have h3 : ω ^ 3 = ω * Complex.I := by
    rw [show (3 : Nat) = 2 + 1 from rfl, pow_succ, ω_pow_two]
    ring
  rw [h3, ω, ep_quarter]
  push_cast
  ring_nf
  rw [Complex.I_sq]
  ring

/-- `a + b·ω + c·ω² + d·ω³`, with `ω⁴ = -1`. -/
structure Cyc where
  /-- Coefficient of `1`. -/
  a : ℤ
  /-- Coefficient of `ω`. -/
  b : ℤ
  /-- Coefficient of `ω²`. -/
  c : ℤ
  /-- Coefficient of `ω³`. -/
  d : ℤ
deriving DecidableEq, Repr, Inhabited

namespace Cyc

@[ext] theorem ext {x y : Cyc} (ha : x.a = y.a) (hb : x.b = y.b) (hc : x.c = y.c)
    (hd : x.d = y.d) : x = y := by
  cases x; cases y; simp_all

instance : Zero Cyc := ⟨⟨0, 0, 0, 0⟩⟩
instance : One Cyc := ⟨⟨1, 0, 0, 0⟩⟩
instance : Add Cyc := ⟨fun x y => ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩⟩
instance : Sub Cyc := ⟨fun x y => ⟨x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d⟩⟩
instance : Neg Cyc := ⟨fun x => ⟨-x.a, -x.b, -x.c, -x.d⟩⟩

/-- Multiplication in `ℤ[ω]`, reducing `ω⁴` to `-1`. -/
instance : Mul Cyc := ⟨fun x y =>
  ⟨x.a * y.a - x.b * y.d - x.c * y.c - x.d * y.b,
   x.a * y.b + x.b * y.a - x.c * y.d - x.d * y.c,
   x.a * y.c + x.b * y.b + x.c * y.a - x.d * y.d,
   x.a * y.d + x.b * y.c + x.c * y.b + x.d * y.a⟩⟩

@[simp] theorem zero_a : (0 : Cyc).a = 0 := rfl
@[simp] theorem zero_b : (0 : Cyc).b = 0 := rfl
@[simp] theorem zero_c : (0 : Cyc).c = 0 := rfl
@[simp] theorem zero_d : (0 : Cyc).d = 0 := rfl
@[simp] theorem one_a : (1 : Cyc).a = 1 := rfl
@[simp] theorem one_b : (1 : Cyc).b = 0 := rfl
@[simp] theorem one_c : (1 : Cyc).c = 0 := rfl
@[simp] theorem one_d : (1 : Cyc).d = 0 := rfl
@[simp] theorem add_a (x y : Cyc) : (x + y).a = x.a + y.a := rfl
@[simp] theorem add_b (x y : Cyc) : (x + y).b = x.b + y.b := rfl
@[simp] theorem add_c (x y : Cyc) : (x + y).c = x.c + y.c := rfl
@[simp] theorem add_d (x y : Cyc) : (x + y).d = x.d + y.d := rfl
@[simp] theorem sub_a (x y : Cyc) : (x - y).a = x.a - y.a := rfl
@[simp] theorem sub_b (x y : Cyc) : (x - y).b = x.b - y.b := rfl
@[simp] theorem sub_c (x y : Cyc) : (x - y).c = x.c - y.c := rfl
@[simp] theorem sub_d (x y : Cyc) : (x - y).d = x.d - y.d := rfl

/-- Multiply by `ω`: the coefficients rotate, and `ω³·ω = -1`. -/
def mulOmega (x : Cyc) : Cyc := ⟨-x.d, x.a, x.b, x.c⟩

/-- Multiply by `ω^p`, `p` taken mod 8 — the eight Clifford+T global phases. -/
def timesOmega (x : Cyc) : Nat → Cyc
  | 0 => x
  | p + 1 => (x.timesOmega p).mulOmega

/-- Whether the numerator is divisible by `√2 = ω - ω³`. Rust's parity test, verbatim. -/
def divisibleBySqrt2 (x : Cyc) : Bool := (x.a + x.c) % 2 == 0 && (x.b + x.d) % 2 == 0

/-- Divide a numerator that `divisibleBySqrt2` accepts by `√2`. -/
def divSqrt2 (x : Cyc) : Cyc :=
  ⟨(x.b - x.d) / 2, (x.a + x.c) / 2, (x.b + x.d) / 2, (x.c - x.a) / 2⟩

/-! ## Interpretation in `ℂ` -/

/-- The complex number a coefficient tuple denotes. -/
noncomputable def interp (x : Cyc) : ℂ :=
  (x.a : ℂ) + (x.b : ℂ) * ω + (x.c : ℂ) * ω ^ 2 + (x.d : ℂ) * ω ^ 3

@[simp] theorem interp_zero : interp 0 = 0 := by simp [interp]

@[simp] theorem interp_one : interp 1 = 1 := by simp [interp]

@[simp] theorem interp_add (x y : Cyc) : interp (x + y) = interp x + interp y := by
  simp only [interp, add_a, add_b, add_c, add_d]
  push_cast
  ring

@[simp] theorem interp_sub (x y : Cyc) : interp (x - y) = interp x - interp y := by
  simp only [interp, sub_a, sub_b, sub_c, sub_d]
  push_cast
  ring

@[simp] theorem interp_mulOmega (x : Cyc) : interp x.mulOmega = ω * interp x := by
  simp only [interp, mulOmega]
  push_cast
  have h4 : ω ^ 4 = -1 := ω_pow_four
  calc (-(x.d : ℂ)) + (x.a : ℂ) * ω + (x.b : ℂ) * ω ^ 2 + (x.c : ℂ) * ω ^ 3
      = (x.a : ℂ) * ω + (x.b : ℂ) * ω ^ 2 + (x.c : ℂ) * ω ^ 3 + (x.d : ℂ) * (-1) := by ring
    _ = (x.a : ℂ) * ω + (x.b : ℂ) * ω ^ 2 + (x.c : ℂ) * ω ^ 3 + (x.d : ℂ) * ω ^ 4 := by
        rw [h4]
    _ = ω * ((x.a : ℂ) + (x.b : ℂ) * ω + (x.c : ℂ) * ω ^ 2 + (x.d : ℂ) * ω ^ 3) := by ring

@[simp] theorem interp_timesOmega (x : Cyc) (p : Nat) :
    interp (x.timesOmega p) = ω ^ p * interp x := by
  induction p with
  | zero => simp [timesOmega]
  | succ p ih => rw [timesOmega, interp_mulOmega, ih, pow_succ]; ring

/-- **The division test is right.** When the parity check accepts, dividing by `√2` and
multiplying back returns the original number — so removing a common `√2` from a matrix is
exact, never approximate. -/
theorem divSqrt2_spec {x : Cyc} (h : x.divisibleBySqrt2 = true) :
    ((Real.sqrt 2 : ℝ) : ℂ) * interp x.divSqrt2 = interp x := by
  simp only [divisibleBySqrt2, Bool.and_eq_true, beq_iff_eq] at h
  obtain ⟨hac, hbd⟩ := h
  obtain ⟨u, hu⟩ : (2 : ℤ) ∣ x.a + x.c := Int.dvd_of_emod_eq_zero hac
  obtain ⟨v, hv⟩ : (2 : ℤ) ∣ x.b + x.d := Int.dvd_of_emod_eq_zero hbd
  -- The four halved coefficients, as integers.
  have hb_d : x.b - x.d = 2 * (v - x.d) := by omega
  have hc_a : x.c - x.a = 2 * (u - x.a) := by omega
  have e1 : (x.b - x.d) / 2 = v - x.d := by omega
  have e2 : (x.a + x.c) / 2 = u := by omega
  have e3 : (x.b + x.d) / 2 = v := by omega
  have e4 : (x.c - x.a) / 2 = u - x.a := by omega
  have hxa : x.a = 2 * u - x.c := by omega
  have hxb : x.b = 2 * v - x.d := by omega
  simp only [interp, divSqrt2, e1, e2, e3, e4, sqrt_two_eq]
  push_cast
  rw [hxa, hxb]
  have h4 : ω ^ 4 = -1 := ω_pow_four
  have h5 : ω ^ 5 = -ω := by
    rw [show (5 : Nat) = 4 + 1 from rfl, pow_succ, h4]; ring
  have h6 : ω ^ 6 = -ω ^ 2 := by
    rw [show (6 : Nat) = 4 + 2 from rfl, pow_add, h4]; ring
  push_cast
  ring_nf
  rw [h4, h5, h6]
  ring

end Cyc

end TzapLean
