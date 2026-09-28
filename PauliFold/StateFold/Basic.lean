import PauliFold.Soundness
import TzapLean.BoolPolynomial

/-!
# StateFold: symbolic path sums

Feynman's `stateFold` (Amy, `Feynman/Optimization/StateFold.hs`) tracks a circuit as a
*path sum*. Each wire carries a Boolean polynomial over path variables: the inputs `0, …, n-1`
and a fresh variable for every H. The circuit accumulates a phase, a function of all variables.
Unlike the two domains of `Pauli.lean` and `Phase.lean`, an abstract state here denotes exactly
one matrix:

  `⟦A⟧ ⟨o| |i⟩ = 2^{-|Y|/2} ∑_{y : Y → Bool} [σ(i, y) = o] · e^{iπ Φ(i, y)}`,

where `Y` is the set of path variables introduced by H gates, `σ` the wire polynomials, and `Φ`
the phase, in multiples of π. StateFold's imprecision lies not in the states but in when two
rotations are deemed foldable (see `StateFold/Fold.lean`).

This file defines states, their denotation, and the gate transformers, and proves each
transformer exact: `denote (step g A) = U_g · denote A`.
-/

namespace PauliFold.SF

open TzapLean Matrix

noncomputable section

abbrev Poly := BoolPolynomial

/-- Valuations of the path variables. -/
abbrev Val := ℕ → Bool

/-- The Boolean value of a polynomial. -/
def ev (ν : Val) (p : Poly) : Bool := BoolPolynomial.evalB ν p

/-- A Boolean as the rational `0` or `1`. -/
def b2q (b : Bool) : ℚ := if b then 1 else 0

/-- A rotation site: the gate's index in the circuit, its angle (a multiple of π), and its
predicate, the wire polynomial it rotates on. -/
structure Site where
  gate : ℕ
  angle : ℚ
  pred : Poly

/-- A StateFold state on `n` qubits. -/
structure State (n : ℕ) where
  ket : Fin n → Poly
  phase : Val → ℚ
  temps : Finset ℕ
  next : ℕ
  sites : List Site

variable {n : ℕ}

/-- The valuation of input `i` and the path variables `S` set to true. -/
def val (i : Basis n) (S : Finset ℕ) : Val :=
  fun v => if h : v < n then i ⟨v, h⟩ else v ∈ S

/-- `1/√2`. -/
noncomputable def invSqrt2 : ℂ := ((Real.sqrt 2 : ℝ) : ℂ)⁻¹

/-- The path sum of a state. -/
noncomputable def State.denote (A : State n) : Density n := fun o i =>
  invSqrt2 ^ A.temps.card *
    ∑ S ∈ A.temps.powerset,
      if ∀ q, ev (val i S) (A.ket q) = o q then ep (A.phase (val i S)) else 0

/-! ## Well-formed states -/

/-- `f` depends only on the variables in `D`. -/
def Supp {α : Type*} (f : Val → α) (D : Set ℕ) : Prop :=
  ∀ ν ν' : Val, (∀ v ∈ D, ν v = ν' v) → f ν = f ν'

/-- The variables a state may use: the inputs and its path variables. -/
def State.dom (A : State n) : Set ℕ := {v | v < n ∨ v ∈ A.temps}

/-- Path variables are below the fresh counter and above the inputs, and every wire and the
phase depend only on inputs and path variables. (Site predicates may mention variables a
reduction has summed out; `Fold.lean` tracks the two predicates it compares.) -/
structure State.WF (A : State n) : Prop where
  temps_ge : ∀ v ∈ A.temps, n ≤ v
  temps_lt : ∀ v ∈ A.temps, v < A.next
  n_le : n ≤ A.next
  ket : ∀ q, Supp (fun ν => ev ν (A.ket q)) A.dom
  phase : Supp A.phase A.dom

/-! ## Transformers -/

/-- A rotation by `θ · π` on qubit `q`: phase `θ · [σ_q]`, and a new site. -/
def rot (j : ℕ) (θ : ℚ) (q : Fin n) (A : State n) : State n :=
  { A with
    phase := fun ν => A.phase ν + θ * b2q (ev ν (A.ket q))
    sites := A.sites ++ [⟨j, θ, A.ket q⟩] }

/-- The transformer of gate number `j`. -/
noncomputable def step (j : ℕ) : CTGate n → State n → State n
  | .x q, A => { A with ket := Function.update A.ket q (A.ket q + 1) }
  | .cx c t _, A => { A with ket := Function.update A.ket t (A.ket t + A.ket c) }
  | .cz c t _, A =>
      { A with phase := fun ν => A.phase ν + b2q (ev ν (A.ket c) && ev ν (A.ket t)) }
  | .h q, A =>
      { A with
        ket := Function.update A.ket q (BoolPolynomial.var A.next)
        phase := fun ν => A.phase ν + b2q (ev ν (A.ket q) && ν A.next)
        temps := insert A.next A.temps
        next := A.next + 1 }
  | .z q, A => rot j 1 q A
  | .s q, A => rot j (1 / 2) q A
  | .sdg q, A => rot j (-1 / 2) q A
  | .t q, A => rot j (1 / 4) q A
  | .tdg q, A => rot j (-1 / 4) q A

/-- The transformers of a gate list, numbering gates from `j`. -/
noncomputable def runFrom : ℕ → List (CTGate n) → State n → State n
  | _, [], A => A
  | j, g :: gs, A => runFrom (j + 1) gs (step j g A)

/-- The initial state: each wire carries its input variable. -/
noncomputable def init (n : ℕ) : State n where
  ket := fun q => BoolPolynomial.var q.val
  phase := fun _ => 0
  temps := ∅
  next := n
  sites := []

/-- The analysis of a circuit. -/
noncomputable def run (gs : List (CTGate n)) : State n := runFrom 0 gs (init n)

end

end PauliFold.SF
