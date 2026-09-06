import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Data.Matrix.Basic
import TzapLean.Circuit

/-!
# Density-Matrix Semantics

Circuits in `TzapLean/RawCircuit.lean` contain `measure` and `reset`, so they are not unitary
maps: they are *quantum channels* that also write classical bits. This file gives them an
exact (no floating point) semantics in three layers.

## Layer 1 — gate matrices

Every unitary gate denotes a `2ⁿ × 2ⁿ` complex matrix `gateUnitary n g`, indexed directly
by computational basis states `Basis n = Fin n → Bool`. Rows are outputs and columns are
inputs, so a gate acts by left multiplication, matching `src/unitary.rs`.

Single-qubit gates go through `embed1`, which lifts a `2 × 2` matrix onto one wire; the
diagonal ones are phases via `ep θ = e^{iπθ}`. The angle convention is the Rust one
exactly:

```
rz θ = diag(e^{-iπθ/2}, e^{iπθ/2})     t = diag(1, e^{iπ/4})     s = diag(1, i)
```

so `t` and `rz (1/4)` are equal only up to the global phase `e^{-iπ/8}` (no
equivalence-up-to-phase relation is defined yet).
`measure` and `reset` denote the identity matrix here; they are not unitary and their real
semantics is in layer 3.

An operand `q ≥ n` (out of range — possible because gate indices are `Nat`, as in Rust)
denotes the identity, so every syntactic circuit denotes something. On `WellFormed`
circuits this case never arises.

## Layer 2 — classical-quantum states

`measure` writes a classical bit, so the state is not a bare density matrix but a
*cq-state*: `CQState n m = Memory m → Density n`, one subnormalized density matrix per
classical-memory value `w : Fin m → Bool`. The trace of `ρ w` is the probability of
memory `w`, and `∑ w, tr (ρ w) = 1` for a normalized state (`IsNormalized`).

## Layer 3 — the channel semantics

`step` gives each gate its action on cq-states:

* a unitary gate `g` conjugates, memory-wise: `ρ w ↦ U ρ w U†`;
* `reset q` applies the Kraus pair `K_b = |0⟩⟨b|_q`: `ρ w ↦ ∑_b K_b ρ w K_b†`;
* `measure q c` projects and *steers the classical memory*: the outcome bit is written to
  `c`, so the new block at memory `w` collects both old blocks that differ at `c`:
  `ρ' w = P_{w c} (ρ w[c←0] + ρ w[c←1]) P_{w c}`. With `c` out of range the outcome is
  discarded instead: `ρ' w = ∑_b P_b (ρ w) P_b`.

`denote` runs a gate list head-first, and `RawCircuit.denote` is the semantics of a circuit at
its own qubit/cbit counts.

The two layers agree where they must: on a measurement-free circuit the channel semantics
is conjugation by the unitary (`denote_eq_conj_unitary`), which is the fact every
optimization pass in the Rust implementation relies on.
-/

namespace TzapLean

open Complex Matrix

noncomputable section

/-! ## Basis states -/

/-- Computational basis states of `n` qubits: one Boolean per wire. -/
abbrev Basis (n : Nat) := Fin n → Bool

namespace Basis

/-- Read qubit `q`; out-of-range reads give `false`. -/
def get {n : Nat} (b : Basis n) (q : Qubit) : Bool :=
  if h : q < n then b ⟨q, h⟩ else false

/-- Write `v` to qubit `q`; out-of-range writes do nothing. -/
def set {n : Nat} (b : Basis n) (q : Qubit) (v : Bool) : Basis n :=
  fun r => if (r : Nat) = q then v else b r

/-- Negate qubit `q` — the action of `x` on basis states. -/
def flip {n : Nat} (b : Basis n) (q : Qubit) : Basis n :=
  b.set q (!b.get q)

@[simp] theorem get_set_same {n : Nat} (b : Basis n) (q : Qubit) (v : Bool) (h : q < n) :
    (b.set q v).get q = v := by simp [get, set, h]

@[simp] theorem get_set_ne {n : Nat} (b : Basis n) (q r : Qubit) (v : Bool) (h : r ≠ q) :
    (b.set q v).get r = b.get r := by
  by_cases hr : r < n <;> simp [get, set, hr, h]

@[simp] theorem set_out_of_range {n : Nat} (b : Basis n) (q : Qubit) (v : Bool)
    (h : ¬ q < n) : b.set q v = b := by
  funext r
  have : (r : Nat) ≠ q := by
    intro hr; exact h (hr ▸ r.isLt)
  simp [set, this]

end Basis

/-! ## Complex phases at rational multiples of π -/

/-- `ep θ = e^{iπθ}`, the phase of a rotation by `θ · π`. Rotation angles in
`Gate.rz` are rational multiples of `π`, so every phase in the semantics is `ep` of a
rational. -/
def ep (θ : ℚ) : ℂ := Complex.exp ((Real.pi * (θ : ℝ) : ℝ) * Complex.I)

@[simp] theorem ep_zero : ep 0 = 1 := by simp [ep]

theorem ep_add (θ φ : ℚ) : ep (θ + φ) = ep θ * ep φ := by
  simp only [ep, ← Complex.exp_add]
  congr 1
  push_cast
  ring

theorem ep_ne_zero (θ : ℚ) : ep θ ≠ 0 := Complex.exp_ne_zero _

@[simp] theorem norm_ep (θ : ℚ) : ‖ep θ‖ = 1 := by
  simp only [ep]
  simpa using Complex.norm_exp_ofReal_mul_I (Real.pi * (θ : ℝ))

theorem star_ep (θ : ℚ) : star (ep θ) = ep (-θ) := by
  simp only [ep, RCLike.star_def, ← Complex.exp_conj]
  congr 1
  push_cast
  simp [Complex.conj_I, Complex.conj_ofReal]

theorem ep_mul_star (θ : ℚ) : ep θ * star (ep θ) = 1 := by
  rw [star_ep, ← ep_add]; simp

/-! ## Layer 1: gate matrices -/

/-- A `2ⁿ × 2ⁿ` complex matrix over basis states: rows are outputs, columns inputs. Also
the type of density matrices. -/
abbrev Density (n : Nat) := Matrix (Basis n) (Basis n) ℂ

/-- Lift a one-qubit matrix `M` (indexed `M output input`) onto wire `q` of an `n`-qubit
register: it acts as `M` on `q` and as the identity elsewhere. An out-of-range `q` gives
the identity. -/
def embed1 (n : Nat) (M : Bool → Bool → ℂ) (q : Qubit) : Density n :=
  if h : q < n then
    fun out inp =>
      if ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
      then M (out ⟨q, h⟩) (inp ⟨q, h⟩)
      else 0
  else 1

/-- `embed1`'s entries, on a wire the register has.

The definition branches on `q < n` at the level of a whole `Density n`, and `Matrix` is no
longer transparent enough for `rw [dif_pos h]` to see through that — the motive fails to
typecheck at reducible transparency. Rewriting entry by entry with this instead is both what
the proofs want and independent of how transparent `Matrix` is. -/
theorem embed1_apply_of_lt {n : Nat} (M : Bool → Bool → ℂ) {q : Qubit} (h : q < n)
    (out inp : Basis n) :
    embed1 n M q out inp =
      if ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
      then M (out ⟨q, h⟩) (inp ⟨q, h⟩) else 0 := by
  simp [embed1, h]

/-- An out-of-range wire embeds as the identity. -/
theorem embed1_eq_one {n : Nat} (M : Bool → Bool → ℂ) {q : Qubit} (h : ¬ q < n) :
    embed1 n M q = 1 := by
  simp [embed1, h]

/-- The one-qubit matrix `diag(a, b)`. -/
def diag2 (a b : ℂ) : Bool → Bool → ℂ :=
  fun out inp => if out = inp then (if inp then b else a) else 0

/-- The one-qubit `X` matrix. -/
def x2 : Bool → Bool → ℂ := fun out inp => if out = !inp then 1 else 0

/-- The one-qubit Hadamard matrix `(1/√2) · [[1, 1], [1, -1]]`. -/
def h2 : Bool → Bool → ℂ :=
  fun out inp => (if out && inp then -1 else 1) / (Real.sqrt 2 : ℂ)

/-- The permutation matrix of a map on basis states: entry `(out, inp)` is `1` when
`out = σ inp`. Used for `x`, `cnot`, and `ccx`. -/
def permMatrix {n : Nat} (σ : Basis n → Basis n) : Density n :=
  fun out inp => if out = σ inp then 1 else 0

/-- The diagonal matrix that multiplies basis state `b` by `f b`. Used for the phase gates
`z`, `s`, `sdg`, `t`, `tdg`, `rz`, `cz`, and `ccz`. -/
def phaseMatrix {n : Nat} (f : Basis n → ℂ) : Density n :=
  fun out inp => if out = inp then f inp else 0

/-- The matrix of a gate, in the same convention as `src/unitary.rs`:

```
x = [[0,1],[1,0]]     h = (1/√2)[[1,1],[1,-1]]     z = diag(1,-1)
s = diag(1, i)        sdg = diag(1,-i)             t = diag(1, e^{iπ/4})
tdg = diag(1, e^{-iπ/4})                           rz θ = diag(e^{-iπθ/2}, e^{iπθ/2})
```

with `cnot`/`ccx` permutations and `cz`/`ccz` diagonal `±1`. `measure` and `reset` are not
unitary; they denote the identity here and get their real semantics in `step`. -/
def gateUnitary (n : Nat) : Gate → Density n
  | .x q => embed1 n x2 q
  | .h q => embed1 n h2 q
  | .s q => embed1 n (diag2 1 (ep (1/2))) q
  | .sdg q => embed1 n (diag2 1 (ep (-1/2))) q
  | .z q => embed1 n (diag2 1 (ep 1)) q
  | .t q => embed1 n (diag2 1 (ep (1/4))) q
  | .tdg q => embed1 n (diag2 1 (ep (-1/4))) q
  | .rz θ q => embed1 n (diag2 (ep (-θ/2)) (ep (θ/2))) q
  | .cnot c t => permMatrix fun b => b.set t (b.get t != b.get c)
  | .ccx c₁ c₂ t => permMatrix fun b => b.set t (b.get t != (b.get c₁ && b.get c₂))
  | .cz c t => phaseMatrix fun b => if b.get c && b.get t then -1 else 1
  | .ccz c₁ c₂ t => phaseMatrix fun b => if b.get c₁ && b.get c₂ && b.get t then -1 else 1
  | .measure _ _ => 1
  | .reset _ => 1

/-- Unitary semantics of a gate list: gates execute head-first, so each new gate is applied
on the *left* of the accumulated matrix. Only meaningful for measurement-free lists; see
`denote_eq_conj_unitary`. -/
def unitary (n : Nat) : List Gate → Density n
  | [] => 1
  | g :: gs => unitary n gs * gateUnitary n g

@[simp] theorem unitary_nil (n : Nat) : unitary n [] = 1 := rfl

@[simp] theorem unitary_cons (n : Nat) (g : Gate) (gs : List Gate) :
    unitary n (g :: gs) = unitary n gs * gateUnitary n g := rfl

theorem unitary_append (n : Nat) (gs hs : List Gate) :
    unitary n (gs ++ hs) = unitary n hs * unitary n gs := by
  induction gs with
  | nil => simp
  | cons g gs ih => simp [ih, Matrix.mul_assoc]

/-- The unitary matrix of a whole circuit, at its own qubit count. -/
def RawCircuit.unitary (c : RawCircuit) : Density c.numQubits :=
  _root_.TzapLean.unitary c.numQubits c.gates

/-! ## Layer 2: classical-quantum states -/

/-- A classical memory: one bit per classical wire. -/
abbrev Memory (m : Nat) := Fin m → Bool

/-- Write `v` to classical bit `c`; out-of-range writes do nothing. -/
def Memory.write {m : Nat} (w : Memory m) (c : CBit) (v : Bool) : Memory m :=
  fun i => if (i : Nat) = c then v else w i

/-- Read classical bit `c`; out-of-range reads give `false`. -/
def Memory.read {m : Nat} (w : Memory m) (c : CBit) : Bool :=
  if h : c < m then w ⟨c, h⟩ else false

/-- Flip classical bit `c` — an involution on memories, and the pairing that makes
`measure` trace preserving. -/
def Memory.flipAt {m : Nat} (c : CBit) (w : Memory m) : Memory m :=
  w.write c (!w.read c)

@[simp] theorem Memory.read_write_same {m : Nat} (w : Memory m) (c : CBit) (v : Bool)
    (h : c < m) : (w.write c v).read c = v := by simp [read, write, h]

@[simp] theorem Memory.write_write {m : Nat} (w : Memory m) (c : CBit) (v v' : Bool) :
    (w.write c v).write c v' = w.write c v' := by
  funext i; by_cases h : (i : Nat) = c <;> simp [write, h]

theorem Memory.write_read_self {m : Nat} (w : Memory m) (c : CBit) (h : c < m) :
    w.write c (w.read c) = w := by
  funext i
  by_cases hi : (i : Nat) = c
  · have : i = ⟨c, h⟩ := Fin.ext hi
    subst this; simp [write, read, h]
  · simp [write, hi]

@[simp] theorem Memory.write_flipAt {m : Nat} (w : Memory m) (c : CBit) (v : Bool) :
    (w.flipAt c).write c v = w.write c v := by simp [flipAt]

theorem Memory.flipAt_flipAt {m : Nat} (w : Memory m) (c : CBit) (h : c < m) :
    (w.flipAt c).flipAt c = w := by
  simp [flipAt, h, Memory.write_read_self w c h]

/-- Flipping one classical bit, as an equivalence — used to re-index sums over memories. -/
def Memory.flipEquiv {m : Nat} (c : CBit) (h : c < m) : Memory m ≃ Memory m where
  toFun := Memory.flipAt c
  invFun := Memory.flipAt c
  left_inv w := Memory.flipAt_flipAt w c h
  right_inv w := Memory.flipAt_flipAt w c h

/-- A classical-quantum state over `n` qubits and `m` classical bits: a subnormalized
density matrix for each value of the classical memory. `tr (ρ w)` is the probability of
observing memory `w`. -/
abbrev CQState (n m : Nat) := Memory m → Density n

/-- Total trace of a cq-state: `1` for a physical state. -/
def CQState.totalTrace {n m : Nat} (ρ : CQState n m) : ℂ :=
  ∑ w : Memory m, Matrix.trace (ρ w)

/-- A cq-state is normalized when its blocks' traces sum to one. -/
def CQState.IsNormalized {n m : Nat} (ρ : CQState n m) : Prop := ρ.totalTrace = 1

/-- The cq-state that is the pure density matrix `ρ` with all classical bits `0`. -/
def CQState.ofDensity {n m : Nat} (ρ : Density n) : CQState n m :=
  fun w => if w = (fun _ => false) then ρ else 0

/-- The all-zeros input state `|0…0⟩⟨0…0|` with cleared classical memory. -/
def CQState.zero (n m : Nat) : CQState n m :=
  CQState.ofDensity (fun out inp =>
    if out = (fun _ => false) ∧ inp = (fun _ => false) then 1 else 0)

/-! ## Layer 3: the channel semantics -/

/-- Conjugation `ρ ↦ U ρ U†`, the action of a unitary on a density matrix. -/
def conj {n : Nat} (U ρ : Density n) : Density n := U * ρ * Uᴴ

/-- The projector onto `qubit q = b`. `P_false + P_true = 1`, and `P_b² = P_b = P_bᴴ`. -/
def proj (n : Nat) (q : Qubit) (b : Bool) : Density n :=
  phaseMatrix fun s => if s.get q = b then 1 else 0

/-- The reset Kraus operator `K_b = |0⟩⟨b|_q`: it keeps the `q = b` branch and rewrites
that qubit to `0`. `∑_b K_bᴴ K_b = 1`, so `reset` is trace preserving. -/
def resetKraus (n : Nat) (q : Qubit) (b : Bool) : Density n :=
  fun out inp => if inp.get q = b ∧ out = inp.set q false then 1 else 0

/-- The action of a single gate on a cq-state — the heart of this file.

* Unitary gates conjugate each memory block by the gate matrix.
* `reset q` applies the Kraus pair `K_b = |0⟩⟨b|_q` and sums, discarding the outcome.
* `measure q c` projects onto each outcome and records it in classical bit `c`: the block
  at memory `w` is built from the two old blocks that agree with `w` away from `c`, each
  projected onto the outcome `w` records. Overwriting `c` is what makes this a sum. If `c`
  is out of range the outcome is discarded, `ρ w ↦ ∑_b P_b (ρ w) P_b`. -/
def step {n m : Nat} (g : Gate) (ρ : CQState n m) : CQState n m :=
  match g with
  | .measure q c =>
      if c < m then
        fun w => conj (proj n q (w.read c)) (ρ (w.write c false) + ρ (w.write c true))
      else
        fun w => ∑ b : Bool, conj (proj n q b) (ρ w)
  | .reset q =>
      fun w => ∑ b : Bool, (resetKraus n q b) * ρ w * (resetKraus n q b)ᴴ
  | g => fun w => conj (gateUnitary n g) (ρ w)

/-- Semantics of a gate list: gates execute head-first. -/
def denote {n m : Nat} : List Gate → CQState n m → CQState n m
  | [], ρ => ρ
  | g :: gs, ρ => denote gs (step g ρ)

/-- Semantics of a circuit at its own qubit and classical-bit counts. -/
def RawCircuit.denote (c : RawCircuit) :
    CQState c.numQubits c.numCbits → CQState c.numQubits c.numCbits :=
  _root_.TzapLean.denote c.gates

/-- The output distribution over classical memories: `outcome ρ w` is the probability that
the classical memory ends up holding `w`. -/
def CQState.outcome {n m : Nat} (ρ : CQState n m) (w : Memory m) : ℂ := Matrix.trace (ρ w)

@[simp] theorem denote_nil {n m : Nat} (ρ : CQState n m) : denote [] ρ = ρ := rfl

@[simp] theorem denote_cons {n m : Nat} (g : Gate) (gs : List Gate) (ρ : CQState n m) :
    denote (g :: gs) ρ = denote gs (step g ρ) := rfl

theorem denote_append {n m : Nat} (gs hs : List Gate) (ρ : CQState n m) :
    denote (gs ++ hs) ρ = denote hs (denote gs ρ) := by
  induction gs generalizing ρ with
  | nil => simp
  | cons g gs ih => simp [ih]

/-! ## The two layers agree on measurement-free circuits -/

/-- On a gate that is not `measure`/`reset`, `step` is conjugation by the gate matrix. -/
theorem step_of_isUnitary {n m : Nat} {g : Gate} (h : g.isUnitary = true)
    (ρ : CQState n m) : step g ρ = fun w => conj (gateUnitary n g) (ρ w) := by
  cases g <;> simp_all [step, Gate.isUnitary, Gate.isMeasurement]

theorem conj_one {n : Nat} (ρ : Density n) : conj 1 ρ = ρ := by
  simp [conj]

theorem conj_mul {n : Nat} (U V ρ : Density n) : conj (U * V) ρ = conj U (conj V ρ) := by
  simp [conj, Matrix.mul_assoc, Matrix.conjTranspose_mul]

/-- **Main bridge theorem.** A circuit without `measure` or `reset` acts on every block of
a cq-state by conjugation with its unitary matrix — the fact that licenses reasoning about
measurement-free circuits purely in terms of `unitary`, as the Rust optimizer does. -/
theorem denote_eq_conj_unitary {n m : Nat} (gs : List Gate)
    (h : ∀ g ∈ gs, g.isUnitary = true) (ρ : CQState n m) :
    denote gs ρ = fun w => conj (unitary n gs) (ρ w) := by
  induction gs generalizing ρ with
  | nil => funext w; simp [conj_one]
  | cons g gs ih =>
      have hg : g.isUnitary = true := h g (by simp)
      have hgs : ∀ g' ∈ gs, g'.isUnitary = true := fun g' hg' => h g' (by simp [hg'])
      funext w
      rw [denote_cons, step_of_isUnitary hg, ih hgs]
      simp [unitary_cons, conj_mul]

/-- The circuit-level form of `denote_eq_conj_unitary`, phrased with the `hasMeasurement`
flag that the Rust implementation maintains. -/
theorem RawCircuit.denote_eq_conj_unitary (c : RawCircuit) (n m : Nat)
    (hn : c.numQubits = n) (hm : c.numCbits = m)
    (h : ∀ g ∈ c.gates, g.isMeasurement = false) (ρ : CQState n m) :
    _root_.TzapLean.denote c.gates ρ =
      fun w => conj (_root_.TzapLean.unitary n c.gates) (ρ w) := by
  subst hn; subst hm
  exact _root_.TzapLean.denote_eq_conj_unitary c.gates
    (fun g hg => by simp [Gate.isUnitary, h g hg]) ρ

/-! ## Kraus completeness: `measure` and `reset` are trace preserving -/

theorem proj_conjTranspose (n : Nat) (q : Qubit) (b : Bool) :
    (proj n q b)ᴴ = proj n q b := by
  funext out inp
  by_cases h : out = inp
  · subst h; simp [proj, phaseMatrix, Matrix.conjTranspose_apply]
  · have h' : inp ≠ out := fun hh => h hh.symm
    simp [proj, phaseMatrix, Matrix.conjTranspose_apply, h, h']

theorem proj_mul_self (n : Nat) (q : Qubit) (b : Bool) :
    proj n q b * proj n q b = proj n q b := by
  funext out inp
  simp only [proj, phaseMatrix, Matrix.mul_apply]
  rw [Finset.sum_eq_single out] <;> simp +contextual [eq_comm]

/-- `P_false + P_true = 1`: the measurement projectors are complete, so `measure` is trace
preserving. -/
theorem proj_sum (n : Nat) (q : Qubit) : ∑ b : Bool, proj n q b = (1 : Density n) := by
  funext out inp
  simp only [Matrix.sum_apply, Fintype.sum_bool, proj, phaseMatrix, Matrix.one_apply]
  by_cases h : out = inp
  · subst h; cases hb : out.get q <;> simp
  · simp [h]

/-- Two basis states that agree at qubit `q` and agree after zeroing `q` are equal. -/
theorem Basis.ext_of_get_set {n : Nat} {out inp : Basis n} (q : Qubit)
    (hget : out.get q = inp.get q) (hset : out.set q false = inp.set q false) : out = inp := by
  funext r
  by_cases h : (r : Nat) = q
  · have hq : q < n := h ▸ r.isLt
    have hr : r = ⟨q, hq⟩ := Fin.ext h
    subst hr
    simpa [Basis.get, hq] using hget
  · simpa [Basis.set, h] using congrFun hset r

/-- Entrywise form of `K_bᴴ K_b`: it is the projector onto "`q` reads `b`, and the two
states agree once `q` is zeroed". -/
theorem resetKraus_conjTranspose_mul (n : Nat) (q : Qubit) (b : Bool) (out inp : Basis n) :
    ((resetKraus n q b)ᴴ * resetKraus n q b) out inp =
      if out.get q = b ∧ inp.get q = b ∧ out.set q false = inp.set q false then 1 else 0 := by
  simp only [Matrix.mul_apply, Matrix.conjTranspose_apply, resetKraus]
  have key : ∀ k : Basis n,
      star (if out.get q = b ∧ k = out.set q false then 1 else 0 : ℂ) *
          (if inp.get q = b ∧ k = inp.set q false then 1 else 0)
        = if k = out.set q false then
            (if out.get q = b ∧ inp.get q = b ∧ out.set q false = inp.set q false then 1 else 0)
          else 0 := by
    intro k
    by_cases hk : k = out.set q false
    · subst hk
      by_cases h1 : out.get q = b <;> by_cases h2 : inp.get q = b <;>
        by_cases h3 : out.set q false = inp.set q false <;> simp [h1, h2, h3]
    · simp [hk]
  rw [Finset.sum_congr rfl (fun k _ => key k), Finset.sum_ite_eq' Finset.univ]
  simp

/-- `∑_b K_bᴴ K_b = 1` for the reset Kraus operators, so `reset` is trace preserving. -/
theorem resetKraus_sum (n : Nat) (q : Qubit) :
    ∑ b : Bool, (resetKraus n q b)ᴴ * resetKraus n q b = (1 : Density n) := by
  funext out inp
  simp only [Matrix.sum_apply, Fintype.sum_bool, resetKraus_conjTranspose_mul,
    Matrix.one_apply]
  by_cases h : out = inp
  · subst h; cases hb : out.get q <;> simp
  · have hne : ¬ (out.get q = inp.get q ∧ out.set q false = inp.set q false) := by
      rintro ⟨h1, h2⟩
      exact h (Basis.ext_of_get_set q h1 h2)
    cases hb1 : out.get q <;> cases hb2 : inp.get q <;> simp_all


/-! ## Trace preservation: `measure` and `reset` are channels

A cq-state's total trace is the total probability. Every gate preserves it, so a
normalized state stays normalized: unitary gates by cyclicity of the trace, `measure` and
`reset` by the Kraus completeness relations above. -/

/-- `tr (P ρ P) = tr (P ρ)`, since `P` is a Hermitian idempotent. -/
theorem trace_conj_proj (n : Nat) (q : Qubit) (b : Bool) (ρ : Density n) :
    Matrix.trace (conj (proj n q b) ρ) = Matrix.trace (proj n q b * ρ) := by
  simp only [conj, proj_conjTranspose]
  rw [Matrix.trace_mul_cycle, proj_mul_self]

/-- Summing a measurement's two outcome branches recovers the input trace. -/
theorem trace_sum_proj (n : Nat) (q : Qubit) (ρ : Density n) :
    ∑ b : Bool, Matrix.trace (conj (proj n q b) ρ) = Matrix.trace ρ := by
  simp only [trace_conj_proj]
  rw [← Matrix.trace_sum, ← Finset.sum_mul, proj_sum, Matrix.one_mul]

/-- A unitary gate preserves the total trace. -/
theorem totalTrace_step_of_isUnitary {n m : Nat} {g : Gate} (h : g.isUnitary = true)
    (ρ : CQState n m) (hU : (gateUnitary n g)ᴴ * gateUnitary n g = 1) :
    (step g ρ).totalTrace = ρ.totalTrace := by
  simp only [CQState.totalTrace, step_of_isUnitary h, conj]
  refine Finset.sum_congr rfl fun w _ => ?_
  rw [Matrix.trace_mul_cycle, hU, Matrix.one_mul]

/-- `reset` preserves the total trace: it is the Kraus channel of `resetKraus`, and
`∑ b, K_bᴴ K_b = 1`. -/
theorem totalTrace_step_reset {n m : Nat} (q : Qubit) (ρ : CQState n m) :
    (step (.reset q) ρ).totalTrace = ρ.totalTrace := by
  simp only [CQState.totalTrace, step]
  refine Finset.sum_congr rfl fun w _ => ?_
  rw [Matrix.trace_sum]
  have : ∀ b : Bool,
      Matrix.trace (resetKraus n q b * ρ w * (resetKraus n q b)ᴴ) =
        Matrix.trace ((resetKraus n q b)ᴴ * resetKraus n q b * ρ w) := by
    intro b; rw [Matrix.trace_mul_cycle, Matrix.mul_assoc]
  rw [Finset.sum_congr rfl (fun b _ => this b), ← Matrix.trace_sum, ← Finset.sum_mul,
    resetKraus_sum, Matrix.one_mul]

/-- `measure` preserves the total trace. Unlike the other gates it moves probability
*between* classical-memory blocks, so the proof pairs each memory `w` with the one that
differs at the written bit (`Memory.flipEquiv`): the two branches of the measurement sum
to the trace of the block pair, and each block pair covers each memory twice. -/
theorem totalTrace_step_measure {n m : Nat} (q : Qubit) (c : CBit) (ρ : CQState n m) :
    (step (.measure q c) ρ).totalTrace = ρ.totalTrace := by
  by_cases hc : c < m
  · simp only [CQState.totalTrace, step, hc, if_pos]
    set f : Memory m → ℂ := fun w =>
      Matrix.trace (conj (proj n q (w.read c)) (ρ (w.write c false) + ρ (w.write c true)))
      with hf
    -- the block pair seen from `w` and from its flip is the same
    have hX : ∀ w : Memory m,
        ρ ((w.flipAt c).write c false) + ρ ((w.flipAt c).write c true) =
          ρ (w.write c false) + ρ (w.write c true) := by
      intro w; simp
    -- the two measurement branches of a block pair sum to its trace
    have hpair : ∀ w : Memory m,
        f w + f (w.flipAt c) =
          Matrix.trace (ρ (w.write c false)) + Matrix.trace (ρ (w.write c true)) := by
      intro w
      have hread : (w.flipAt c).read c = !(w.read c) := by
        simp [Memory.flipAt, hc]
      rw [hf]
      simp only [hread, hX w]
      have hsum := trace_sum_proj n q (ρ (w.write c false) + ρ (w.write c true))
      rw [Fintype.sum_bool, Matrix.trace_add] at hsum
      cases hb : w.read c
      · simp only [Bool.not_false]
        exact (add_comm _ _).trans hsum
      · simp only [Bool.not_true]
        exact hsum
    -- each block pair, summed over all memories, covers every block twice
    have hsplit : ∀ w : Memory m,
        Matrix.trace (ρ (w.write c false)) + Matrix.trace (ρ (w.write c true)) =
          Matrix.trace (ρ w) + Matrix.trace (ρ (w.flipAt c)) := by
      intro w
      have hw : w.write c (w.read c) = w := Memory.write_read_self w c hc
      cases hb : w.read c
      · have h0 : w.write c false = w := by rw [← hb]; exact hw
        have h1 : w.write c true = w.flipAt c := by simp [Memory.flipAt, hb]
        rw [h0, h1]
      · have h1 : w.write c true = w := by rw [← hb]; exact hw
        have h0 : w.write c false = w.flipAt c := by simp [Memory.flipAt, hb]
        rw [h0, h1, add_comm]
    have hflip : ∑ w : Memory m, f (w.flipAt c) = ∑ w : Memory m, f w :=
      Equiv.sum_comp (Memory.flipEquiv c hc) f
    have hflip' : ∑ w : Memory m, Matrix.trace (ρ (w.flipAt c)) =
        ∑ w : Memory m, Matrix.trace (ρ w) :=
      Equiv.sum_comp (Memory.flipEquiv c hc) (fun w => Matrix.trace (ρ w))
    have key : (2 : ℂ) * ∑ w : Memory m, f w = 2 * ∑ w : Memory m, Matrix.trace (ρ w) := by
      have left : ∑ w : Memory m, (f w + f (w.flipAt c)) = 2 * ∑ w : Memory m, f w := by
        rw [Finset.sum_add_distrib, hflip]; ring
      have right : ∑ w : Memory m,
          (Matrix.trace (ρ w) + Matrix.trace (ρ (w.flipAt c))) =
            2 * ∑ w : Memory m, Matrix.trace (ρ w) := by
        rw [Finset.sum_add_distrib, hflip']; ring
      rw [← left, ← right]
      exact Finset.sum_congr rfl fun w _ => by rw [hpair w, hsplit w]
    exact mul_left_cancel₀ two_ne_zero key
  · simp only [CQState.totalTrace, step, hc, if_neg, not_false_iff]
    refine Finset.sum_congr rfl fun w _ => ?_
    rw [Matrix.trace_sum, trace_sum_proj]

end
end TzapLean
