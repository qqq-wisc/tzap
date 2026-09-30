# pauli-fold (Lean)

A Lean 4 formalization comparing three abstract domains of Clifford+T circuits: the Pauli
domain of `PauliFoldRand`, symbolic phase folding, and a restricted model of Feynman's
StateFold (with a restricted model of Amy and Lunderville's Strengthen). See **Scope of the
StateFold model** below for what is and is not modeled.

| domain | used by | what it tracks |
|---|---|---|
| **Pauli** | `PauliFoldRand` | how the circuit conjugates Pauli strings: one signed Pauli, or `⊤` |
| **Phase** | `PhaseFoldRand`; the symbolic variant of *Linear time T gate optimization* | affine Boolean forms on the wires |
| **StateFold `d`** (+ **Strengthen**) | a restricted model of Feynman `-statefold d` and Amy and Lunderville, POPL 2025 | an exact path sum, reduced by HH (substitution degree `≤ d`) and ω; only `Z_q ↦ ±Z_{q'}` facts are kept; Strengthen adds witnessed constraints on the final state |

**Scope of the StateFold model.** `StateFold` and `Strengthen` here are *restricted models* of
Amy and Lunderville's analyses, not faithful formalizations of their full domain:

- **Only Z-to-Z facts.** `sfFacts` and `sfFactsS` keep only facts `Z_q ↦ ±Z_{q'}`. The paper
  tracks general affine (§4) and polynomial (§5) transition relations, including relations
  among several wires. For example, `1 ∈ sfGamma d [cx 0 1]` at every `d`, while
  `affGamma [cx 0 1]` rejects the identity: it misses `x'₁ = x₀ ⊕ x₁`.
- **Witnesses at the final state only.** `ProvesS` takes all witnesses on the trace's last
  state, and no constraint may mention a witness variable. Algorithm 2 accumulates
  constraints across every rewrite state and reduces phase predicates modulo the
  accumulated ideal.
- **`d` bounds only HH substitutions.** `HH.deg` checks `R.totalDegree ≤ d`, but `Witness`
  has no degree bound. So `ProvesS 1` may use non-linear constraints, and it is not a
  degree-1 Strengthen.

**What transfers to the full domain.** The full analysis derives at least these facts, so
its concretization is contained in `sfGamma`. Therefore:

- *StateFold proves more than Pauli* (`h; t; tdg; h`, `midE`): this direction transfers.
- *Pauli proves more than StateFold* (the Clifford segment): this direction is now also
  proved against the paper's own path-sum abstraction `α` (`Alpha.lean`, below).

Strengthen's accumulated constraints (Algorithm 2) are still modeled only in restricted form.

## Main results: the domains as abstractions of circuits (`Domains.lean`)

**Concrete semantics.** A circuit `gs` denotes its unitary, `circ gs`.

**Abstractions.** Each domain abstracts a circuit by the set of *Pauli facts* its analysis
derives. A fact `(P, P′)` says that `U P U† = P′`. The concretization is the set of
unitaries satisfying every fact:

```lean
def Holds (U : Density n) (f : Fact n) : Prop := conj U f.1.toMatrix = f.2.toMatrix
def gammaF (F : Set (Fact n)) : Set (Density n) := {U | ∀ f ∈ F, Holds U f}
```

| domain | facts about `gs` |
|---|---|
| Pauli (`pauliFacts`) | `(P, P′)` whenever the Pauli transformer maps `P` to `P′` rather than `⊤` |
| Phase (`phaseFacts`) | `(Z_q, ±Z_Q)` whenever the wires in `Q` sum to the input variable `x_q` |
| StateFold (`sfFacts d`) | `(Z_q, ±Z_{q′})` whenever a degree-`d` trace proves the wire predicate of `q′` equal or complementary to `x_q` |
| Strengthen (`sfFactsS d`) | the same, with witnessed constraints |

A domain is more precise on a circuit when its set of unitaries is smaller.

```lean
theorem pauli_sound (gs) : circ gs ∈ pauliGamma gs      -- also phase_sound, sf_sound, sfS_sound
```

> **In words.** Every domain is sound: the circuit's own unitary satisfies every fact it
> derives.

```lean
theorem pauli_refines_phase (gs) : pauliGamma gs ⊆ phaseGamma gs
theorem pauli_refines_phase_strict : pauliGamma hh ⊂ phaseGamma hh
```

> **In words.** Every phase fact is a Pauli fact, so the Pauli domain is at least as precise
> on every circuit. On `h; h` it is strictly more precise: Pauli proves `Z ↦ Z`, and phase
> proves nothing.

```lean
theorem pauli_stateFold_incomparable (d) (hd : 1 ≤ d) :
    (∃ n gs, ¬ sfGamma d gs ⊆ pauliGamma gs) ∧ (∃ n gs, ¬ pauliGamma gs ⊆ sfGamma d gs)
theorem midE_sf_not_pauli (d) (hd : 2 ≤ d) : ¬ pauliGamma midE ⊆ sfGamma d midE
theorem pauli_strengthen_incomparable (d) (hd : 1 ≤ d) : …   -- the same for Strengthen
```

> **In words.** Neither Pauli nor StateFold refines the other, at any degree `d ≥ 1`. Each
> direction is witnessed by a Z-to-Z fact that one domain proves and the other's facts do
> not imply:
>
> - **Pauli proves more** on `cx q0,q1; h q0; h q1; cx q1,q0; h q1`. It proves
>   `Z₁ ↦ Z₁`. StateFold proves no fact at any degree, so `X₁` satisfies all its facts but
>   sends `Z₁` to `−Z₁` (`pauli_fact_not_sf`).
> - **StateFold proves more** on `h q0; t q0; tdg q0; h q0`. It proves `Z₀ ↦ Z₀`. The
>   circuit starts with `h; t`, so every Pauli that Pauli can track has `I` or `X` at `q0`.
>   Then `U · X₀` satisfies every Pauli fact, but sends `Z₀` to `−Z₀` (`sf_fact_not_pauli`).
> - **The 21-gate relative-Toffoli circuit** `midE` is the same for `d ≥ 2`. On it, the
>   PauliFoldRand pass also misses the fold.

```lean
theorem strengthen_refines_stateFold (d gs) : sfGammaS d gs ⊆ sfGamma d gs
theorem stateFold_anti (hd : d ≤ d') (gs) : sfGamma d' gs ⊆ sfGamma d gs
```

> **In words.** Strengthen refines StateFold, and raising the degree only adds facts.

In summary:

```
Pauli ⊊ Phase     Pauli ⋈ StateFold d     Pauli ⋈ Strengthen d     Strengthen ⊆ StateFold
```

Every theorem depends only on the standard axioms `propext`, `Classical.choice`, and
`Quot.sound`. The concrete semantics is that of [`../lean`](../lean) (`tzap-lean`):
operators are its `Density n` matrices, and gates denote `TzapLean.gateUnitary`.

## The tableau abstraction (`Tableau.lean`)

This is the abstraction the Rust `PauliFoldRand` pass computes, formalized for whole
circuits. It says nothing about merging rotations; it only approximates the circuit's
semantics.

**The state.** A Clifford tableau and a list of axes:

- The rows are `x[q] = C† X_q C` and `z[q] = C† Z_q C`, where `C` is the Clifford part so far.
- There is one input-frame axis per T or T† gate.

Pauli strings carry a phase `i^k` (`PP`), so that products of strings are strings. The
pull-back `back P = C† P C` of any Pauli is computed as a product of rows.

**The transformers.**

- A Clifford gate `g` updates the rows to `back (g† G g)`, using the proved forward table of
  `g⁻¹`. It leaves the axes alone.
- A T or T† on `q` leaves the rows alone and records the current `z[q]` as an axis.

```lean
theorem PP.toMatrix_mul (P Q : PP n) : (P.mul Q).toMatrix = P.toMatrix * Q.toMatrix
theorem Tab.back_ok … : (T.back P).toMatrix = Cᴴ * P.toMatrix * C
theorem Tab.tableau_sound (gs) (P)
    (h : ∀ A ∈ ((init n).run gs).axes, PP.Comm (((init n).run gs).back P) A) :
    conj (circ gs) (((init n).run gs).back P).toMatrix = P.toMatrix
```

> **In words.**
>
> - **When a fact is proved.** After any Clifford+T circuit, take the tableau's pull-back
>   `B` of a Pauli `P`. If `B` commutes with every recorded T axis, the circuit maps `B` to
>   `P` exactly: `U B U† = P`.
> - **Why T is sound.** A T is `α·I + β·Z_q`, so pulled back to the input frame it is
>   `α·I + β·z[q]`. It therefore commutes with every Pauli that commutes with `z[q]`.
> - **The invariant.** `U = C · V`, where the rows are correct for `C` and `V` commutes with
>   every Pauli that commutes with all the axes.
>
> `Tab.gamma gs` is the resulting concretization, and `Tab.circ_mem_gamma` states that the
> circuit's unitary lies in it.

**Containment** (`TableauContainment.lean`):

```lean
theorem tab_gamma_subset_pauli (gs) : Tab.gamma gs ⊆ pauliGamma gs
theorem tableau_refines_phase (gs) : Tab.gamma gs ⊆ phaseGamma gs
theorem tableau_refines_phase_strict : Tab.gamma hh ⊂ phaseGamma hh
```

> **In words.**
>
> - **Every Pauli-domain fact is a tableau fact.** The proof runs both analyses side by
>   side from each input Pauli: while the forward run is at `P`, the rows represent `C`,
>   `C B C† = P`, and `B` commutes with every recorded axis. So the tableau refines the
>   Pauli domain, and through it phase folding, strictly on `h; h`.
> - **`PP.comm_iff`:** the combinatorial commutation test is exactly commutation of the
>   matrices.

## The affine domain as a transition relation (`Affine.lean`)

This is phase folding's domain in the whole-circuit form of Amy and Lunderville's affine
relation analysis.

**The state and transformers** (`AffAbs`). Each wire carries an affine form over the
inputs and one fresh variable per H. X adds 1, CX adds the control's wire to the target's,
H introduces a fresh variable, and diagonal gates change nothing. It runs once, with no
candidate.

**The relations** (`AffAbs.Rel`). Whenever the wires in `S` sum to `x_T ⊕ c` with no path
variables, every run of the circuit satisfies

  `(output parity on S) = (input parity on T) ⊕ c`.

**The concretization is a transition relation:**

```lean
def affGamma (gs) : Set (Density n) :=
  {U | ∀ x x', U x' x ≠ 0 → ∀ r ∈ ((AffAbs.init n).run gs).Rel,
        parity r.1 x' = parity r.2.1 x + r.2.2}

theorem aff_sound (gs) : circ gs ∈ affGamma gs
theorem pauli_refines_aff (gs) : pauliGamma gs ⊆ affGamma gs
theorem tableau_refines_aff (gs) : Tab.gamma gs ⊆ affGamma gs
theorem tableau_refines_aff_strict : Tab.gamma hh ⊂ affGamma hh
```

> **In words.** A relation `(S, T, c)` is the Pauli fact `Z_T ↦ (-1)^c Z_S`
> (`rel_pauliFact`). This generalizes the containment invariant to candidates on any input
> parity `T`. Conversely, a unitary with that fact satisfies `Z_S U = ± U Z_T`, which forces
> the parity constraint on every nonzero entry `⟨x'|U|x⟩` (`support_of_zfact`). Unitarity
> itself comes from the fact `I ↦ I`, which every run proves. So the tableau, the Pauli
> domain, and the affine relation form a chain of refinements, strict on `h; h`.

## Joins, branches and loops (`Join.lean`)

**Partial tableaux.** For joins, the tableau may forget rows: `PTab` has rows in
`Option (PP n)`, where `none` means unknown. The domain is `PAbs n = Option (PTab n)`, with
`none` as `⊤`.

**The join.**

- *Rows:* a row stays known when both branches agree on it, and becomes unknown otherwise.
- *Axes:* the union of both branches' axes. Axes are in the input coordinates of the common
  starting point, so a fact survives only if it commutes with the rotations of both
  branches.

**Transformers.**

- A Clifford gate multiplies rows, and a product involving an unknown row is unknown.
- A T whose row `z[q]` is unknown gives `⊤`.

**Facts.** A fact needs its pull-back to be defined, meaning every row it uses is known.

**Loops.** Iterate `I ← I ⊔ body♯(I)` until `body♯(I) ⊑ I`, checked decidably; after
`2n + 3` tries the loop gives `⊤`.

```lean
theorem absRun_sound (hU : Exec p U) : Sound U₀ a → Sound (U * U₀) (absRun p a)
theorem program_sound (hU : Exec p U) : U ∈ gammaA (absRun p (some (PTab.init n)))
theorem loopSimple_Z (hU : Exec loopSimple U) : conj U Z_a = Z_a     -- t a; while ★ {cx a,b}; tdg a
theorem loopH_Z      (hU : Exec loopH U)      : conj U Z_b = Z_b     -- t b; while ★ {h a};    tdg b
example : absRun loopSwap (some (PTab.init 2)) = none                -- swap forgets all rows
```

> **In words.** A partial state is sound for `U` when some *total* tableau completes it: the
> completion agrees on every known row, its axes are among the state's axes, and it
> satisfies the tableau invariant.
>
> - **Transformers.** Every partial step is matched by the total step on the completion
>   (`stepC_extends`, `back_extends`), so soundness reuses `Tab.inv_step`.
> - **The join.** The join only forgets rows and adds axes, so a completion of either branch
>   completes the join.
> - **Programs.** So every execution of a program, over all branch choices and loop
>   iterations, lies in the concretization.
> - **Examples, checked by `decide`.** `loop-simple` and `loop-h` keep the facts that let
>   their T and T† merge. `loop-swap` goes to `⊤`: its body changes all four rows.

## A bounded disjunctive domain (`Disjunctive.lean`)

An abstract value is a list of partial tableaux (or `⊤`), meaning "the program is described
by one of them".

- **Concretization:** the union of the disjuncts', so a fact must hold in every disjunct.
- **Normalization:** disjuncts with equal rows merge, unioning their axes, so loops with
  rotations in the body stabilize. Above `bound = 8` disjuncts the list collapses with the
  per-row join.
- **Transformers, join and loops:** transformers act per disjunct, the join is
  concatenation followed by normalization, and loops iterate as before.

```lean
theorem absRunD_sound (hU : Exec p U) : SoundD U₀ a → SoundD (U * U₀) (absRunD p a)
theorem programD_sound (hU : Exec p U) : U ∈ gammaD (absRunD p (some [PTab.init n]))
theorem loopSwap_ZZ (hU : Exec loopSwapSeg U) : conj U (Z_a Z_b) = Z_a Z_b    -- while ★ {swap a,b}
```

> **In words.** `while ★ { swap a,b }` has the two disjuncts identity and swap. Both pull
> `Z_aZ_b` back to itself, so `loop-swap`'s T and T† merge, as in Amy and Lunderville. The
> per-row join cannot show this. The cycle `while ★ { swap a,b; swap b,c }` correctly does
> not preserve `Z_b`.

**Not yet proved:**

- The converse `pauliGamma ⊆ Tab.gamma`, which would make the two concretizations equal.
- The StateFold incomparability results for the tableau.
- CCX, CCZ, measurement and reset, which the Rust pass adds as further blockers.

## The fold theorem (`FoldRule.lean`)

The rewrite that `PauliFoldRand` performs, justified by one fact of a segment's tableau
abstraction. This is Theorem 9.1 of `docs/tableau-abstraction.tex`.

Let `m` be the segment between two rotations and `σ = (init n).run m`. The hypotheses are:

- (a) `σ.back Z_{q'} = ± Z_q`: the segment pulls `Z_{q'}` back to `Z_q`, up to sign;
- (b) `Z_q` commutes with every axis of `σ`.

`rot θ P = diag(1, e^{iπθ})` about `P`, so `T = rot (1/4)` and `S = rot (1/2)` exactly.

```lean
theorem segment_Z : circ m * Z_q = i^{-k} • Z_{q'} * circ m          -- from tableau_sound
theorem fold_pos : circ m * rot θ Z_q = rot θ Z_{q'} * circ m                 -- sign +
theorem fold_neg : circ m * rot θ Z_q = (ep θ • rot (-θ) Z_{q'}) * circ m     -- sign −
theorem fold_merge_pos :
    K * rot θ' Z_{q'} * circ m * rot θ Z_q * L = K * rot (θ + θ') Z_{q'} * circ m * L
theorem fold_t_t     : circ (ℓ ++ t q :: m ++ t q' :: k) = circ (ℓ ++ m ++ s q' :: k)
theorem fold_t_t_neg : circ (ℓ ++ t q :: m ++ t q' :: k) = ep (1/4) • circ (ℓ ++ m ++ k)
```

> **In words.** If the segment's abstract state has the fact `(±Z_q, Z_{q'})`, a rotation on
> `q` before the segment can be moved to `q'` after it. There its angle adds to the next
> rotation's, or is subtracted when the sign is `−`. Only this one fact about `m` is
> used. An example checks the worked example's first merge (`cx; t; cx` between two `T`s)
> by `decide`.

**Not yet proved:** that the scan's check, "same axis and every live event in between
commutes", establishes (a) and (b) for the current segment (Lemma 9.3). The version for
programs with branches and loops (Corollary 9.2) is also paper-only.

## Amy and Lunderville's `α` (`Alpha.lean`)

This is the paper's path-sum abstraction (Definition 17), `α = ∃Y. ⟨X' ⊕ f(X, Y)⟩`, and
Algorithm 1, which rewrites the path sum with HH (substitution degree `≤ d`) and ω. Over `F₂`
the variety of `α` is the image of the wire map (their Proposition 12), so its
concretization is:

```lean
def SF.alphaGamma (A : State n) : Set (Density n) :=
  {U | ∀ i o, U o i ≠ 0 → ∃ S ⊆ A.temps, ∀ q, ev (val i S) (A.ket q) = o q}

def alphaGammaC (d) (gs) : Set (Density n) :=          -- every rewrite endpoint
  {U | U * Uᴴ = 1 ∧ ∀ C, SF.Reduces d (SF.run gs) C → U ∈ SF.alphaGamma C}

theorem alpha_sound : circ gs ∈ alphaGammaC d gs                      -- Prop. 18
theorem SF.alphaGamma_mono : Reduces d A C → A.WF → alphaGamma C ⊆ alphaGamma A   -- Prop. 21
theorem alpha_refines_sf : alphaGammaC d gs ⊆ sfGamma d gs
theorem alpha_cx_rejects_id : (1 : Density 2) ∉ alphaGammaC d [cx 0 1]
theorem alpha_cexMid : U * Uᴴ = 1 → U ∈ alphaGammaC d cexMid          -- α is ⊤ here
theorem alpha_pauli_incomparable (hd : 1 ≤ d) :
    (∃ n gs, ¬ alphaGammaC d gs ⊆ pauliGamma gs) ∧ (∃ n gs, ¬ pauliGamma gs ⊆ alphaGammaC d gs)
theorem midE_alpha_not_pauli (hd : 2 ≤ d) : ¬ pauliGamma midE ⊆ alphaGammaC d midE
```

> **In words.**
>
> - `α` keeps whole transition relations, so it rejects the identity on a CNOT, which the
>   Z-fact model accepts (the auditor's regression).
> - Rewriting only makes it more precise, and it refines the Z-fact model.
> - Against the Pauli domain, both separations now hold for `α` itself:
>   - On `cx q0,q1; h q0; h q1; cx q1,q0; h q1`, every path variable is in a wire, so no
>     rule applies. The wire map `(y₂ ⊕ y₃, y₄)` is onto, so `α` is `⊤`, and it admits `X₁`,
>     which breaks `Z₁ ↦ Z₁`.
>   - On `h; t; tdg; h` and `midE` (for `d ≥ 2`), `α` proves `Z₀ ↦ Z₀`, and the Pauli domain
>     allows `U · X₀`.
>
> The `cexMid` separation is a limit of the *rewrite rules*, not of the domain. `α(U)`
> computed exactly would contain `x'₁ = x₁`. A change of variables `y₂ := u ⊕ y₃` would
> expose an HH step, and so does Feynman's `-cxcz` preprocessing.

`alphaGammaC` intersects over *all* rewrite endpoints, so it is at least as precise as any
run of Algorithm 1. The Strengthen ideal inside `α` (Algorithm 2) is not modeled here.

The sections below are the supporting development. They study each domain one input Pauli
at a time (the per-candidate view), and they prove that StateFold's folding rewrite is
correct.

---

## 1. Pauli folding refines phase folding

The concrete transformer of a gate is conjugation, `A ↦ U A U†` (`CTGate.concrete`), and
`concreteRun gs` applies it gate by gate. A candidate on `q` starts from `zString {q} false`,
the Pauli string `+Z_q`.

```lean
theorem containment (q : Fin n) (gs : List (CTGate n)) :
    (PauliAbs.run gs (.one (zString {q} false))).gamma ⊆
      (PhaseAbs.run gs (PhaseAbs.init n q)).gamma
```

> **In words.** For any Clifford+T segment, every operator the Pauli domain allows as the
> candidate's axis is also allowed by phase folding. So phase folding never folds a pair
> that Pauli folding cannot.

```lean
theorem containment_strict :
    ¬ (PhaseAbs.run hh (PhaseAbs.init 1 0)).gamma ⊆
      (PauliAbs.run hh (.one (zString {0} false))).gamma
```

> **In words.** On the one-qubit segment `h; h`, the inclusion is strict. Phase folding
> introduces a fresh variable at each H and loses the axis. The Pauli domain maps `Z` to `X`
> and back to `Z`.

## 2. Soundness of the Pauli and phase domains

```lean
theorem step_sound (g : CTGate n) (p : PauliAbs n) (A : Op n) (hA : A ∈ p.gamma) :
    g.concrete A ∈ (PauliAbs.step g p).gamma

theorem candidate_sound (q : Fin n) (gs : List (CTGate n)) :
    concreteRun gs (zString {q} false).toMatrix ∈
        (PauliAbs.run gs (.one (zString {q} false))).gamma ∧
      concreteRun gs (zString {q} false).toMatrix ∈
        (PhaseAbs.run gs (PhaseAbs.init n q)).gamma
```

> **In words.** Each abstract Pauli transformer over-approximates conjugation by the gate's
> real unitary. Every Clifford conjugation table is proved as a matrix identity, and T and
> T† either fix the Pauli (on `I` or `Z`) or go to `⊤`. The candidate's true axis therefore
> lies in the Pauli concretization after any segment. Through `containment`, it also lies in
> the phase-folding concretization.

## 3. StateFold: exact path sums and sound folding

A StateFold state (`SF.State`) consists of:

- a Boolean polynomial on each wire;
- a phase function `Φ` in multiples of π;
- the set of path variables introduced by H gates.

It denotes the path sum

```
⟦A⟧⟨o|i⟩ = 2^{-|Y|/2} ∑_{y : Y → Bool} [wires(i, y) = o] · e^{iπ Φ(i, y)}.
```

```lean
theorem denote_run (gs : List (CTGate n)) :
    (run gs).denote = unitary n (gs.map CTGate.toGate)
```

> **In words.** The StateFold analysis of a circuit denotes exactly the circuit's unitary.
> Unlike the Pauli and phase domains, the analysis itself loses no information.

The two reduction rules eliminate path variables:

- **HH:** `∑_y (-1)^{y(z ⊕ R)} = 2 [z = R]`. It substitutes `z := R`, and is allowed only
  when `deg R ≤ d`.
- **ω:** it sums out a variable `y` whose phase contribution is `y · (1/2 + R)`.

```lean
theorem denote_hh  (h : HH d A B y z R Φ₀) (hA : A.WF) : B.denote = A.denote
theorem denote_omega (h : Omega A B y R Φ₀) (hA : A.WF) : B.denote = A.denote
```

> **In words.** Both rules preserve the denotation. The degree bound `d` is only a side
> condition, so soundness holds at every `d`.

A `Trace d A f g C f' g'` is a sequence of HH and ω steps from `A` to `C`. It carries two
rotation predicates `f` and `g` to `f'` and `g'`, and never sums out a variable that either
predicate mentions.

```lean
theorem fold_sound {d : ℕ} (pre mid post : List (CTGate n)) (gk gl : CTGate n)
    {θk θl : ℚ} {qk ql : Fin n} (hk : rotAngle gk = some (θk, qk))
    (hθ : θk = 1 / 4 ∨ θk = -1 / 4) (hl : rotAngle gl = some (θl, ql)) (s : Bool)
    {C : State n} {f' g' : Poly}
    (htr : Trace d (run (pre ++ gk :: mid ++ gl :: post)) (predAfter pre qk)
      (predAfter₂ pre mid ql) C f' g')
    (heq : ∀ ν, ev ν g' = (s != ev ν f')) :
    unitary n ((pre ++ mid ++ gl :: tGate ((if s then -1 else 1) * θk) ql :: post).map
        CTGate.toGate) =
      ep (if s then -θk else 0) •
        unitary n ((pre ++ gk :: mid ++ gl :: post).map CTGate.toGate)
```

> **In words.** Let `gₖ` be a T or T† and `gₗ` a later rotation.
>
> Suppose some degree-`d` reduction sequence makes their predicates equal (`s = false`) or
> complementary (`s = true`). Then StateFold's rewrite is correct:
>
> - delete `gₖ`;
> - put its angle, negated if `s`, as a T or T† right after `gₗ`.
>
> The new circuit has the same unitary up to the global phase `e^{-iπθₖ}` (when `s`).
>
> The proof uses a *shifting lemma*. The move changes the phase by `-θₖ[fₖ] + σθₖ[fₗ]`. That
> term rides along the trace unchanged, and at its end it is a constant.

```lean
theorem Trace.mono (hd : d ≤ d') (htr : Trace d A f g C f' g') : Trace d' A f g C f' g'
```

> **In words.** Raising the degree bound never loses a fold.

## 4. StateFold as an abstract domain, and incomparability with Pauli

A StateFold state is exact, so the set of operators consistent with it is just the true
axis. StateFold's imprecision lies elsewhere: in which facts its reductions can *prove*.
Its concretization is therefore defined through those facts.

- **A fact.** `Proves d (run mid) x_{qₖ} q s` holds when a degree-`d` trace makes the
  candidate's predicate `x_{qₖ}` equal (`s = false`) or complementary (`s = true`) to the
  wire predicate of qubit `q`.
- **Concretization.** `sfCandidate d mid qk` is the set of operators `O` such that
  `O = (-1)^s Z_q` for every proved fact. If nothing is proved, it is every operator.
- **For comparison.** `pauliCandidate mid qk` is the Pauli concretization
  `(PauliAbs.run mid (.one (zString {qk} false))).gamma`.

```lean
theorem proves_sound (hp : Proves d (run mid) (BoolPolynomial.var qk.val) q s) :
    axis mid qk = (zString {q} s).toMatrix
```

> **In words.** A fact StateFold proves is true: the candidate's axis after the segment
> really is `(-1)^s Z_q`.
>
> The proof starts from the shifting lemma, which gives `Z_q · U = (-1)^s U · Z_{qₖ}`. Then
> `U` is unitary, and conjugation is `concreteRun`.

```lean
theorem axis_mem_both (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) :
    axis mid qk ∈ sfCandidate d mid qk ∩ pauliCandidate mid qk
```

> **In words.** Both domains are sound, so their concretizations always share the true axis.

```lean
theorem gamma_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      sfCandidate d mid qk ⊂ pauliCandidate mid qk) ∧
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      pauliCandidate mid qk ⊂ sfCandidate d mid qk)
```

> **In words.** At every degree `d ≥ 1`, neither domain refines the other.
>
> - **StateFold is strictly more precise** on `h q0; t q0; tdg q0; h q0`. StateFold's phase
>   is a function, so the inner T and T† cancel exactly. One HH step then proves that the
>   axis is `Z₀`. PauliFold sees the inner T anticommute with `X₀` and goes to `⊤`.
> - **PauliFold is strictly more precise** on `cx q0,q1; h q0; h q1; cx q1,q0; h q1`.
>   PauliFold carries `+Z₁` through the segment exactly. StateFold ends with the wires
>   `y₂ ⊕ y₃` and `y₄`, which mention every path variable. Neither rule may sum out a
>   variable that a wire mentions, so no reduction applies and nothing is proved, at any
>   degree.

**A separation that survives the tools.** `h; t; tdg; h` separates the *domains*, but the
PauliFoldRand *pass* still clears it: it first merges the adjacent inner T and T†, which
share an axis. `DegreeTwo.lean` gives a 21-gate, 11-T circuit on which the tools differ as
well (`circuits/statefold-beats-paulifold.qasm`):

```
t q0; G; tdg q0; G; t q0
G = h q0; tdg q0; cx q2,q0; t q0; cx q1,q0; tdg q0; cx q2,q0; t q0; h q0
```

`G` is a relative-phase Toffoli with target `q0` and controls `q1`, `q2`.

```lean
theorem midE_proves : Proves 2 (run midE) (BoolPolynomial.var 0) 0 false

theorem midE_ssubset (d : ℕ) (hd : 2 ≤ d) : sfCandidate d midE 0 ⊂ pauliCandidate midE 0
```

> **In words.**
>
> - **What the circuit does.** The first `G` computes `x₀ ⊕ x₁x₂` into `q0`. The middle T†
>   acts on that nonlinear value, and the second `G` uncomputes it.
> - **What StateFold proves.** At degree 2, StateFold proves that the outer Ts' axis `Z₀`
>   returns unchanged, so they merge into an S. It takes two HH steps: `y₄ := x₀ ⊕ x₁x₂`
>   (degree 2), then `y₆ := x₀`. The first works because the four Ts of `G` add up to
>   `π · y x₁ x₂` modulo 2.
> - **What PauliFold sees.** PauliFold goes to `⊤` at the first T† inside `G`, where the axis
>   is `X₀`.
>
> The tools agree:
>
> - Feynman `-statefold 2`: 11 → 9 Ts.
> - `-statefold 1`: 11 → 11.
> - tzap's PauliFoldRand + CancelGates fixpoint: 11 → 11.
> - tzap `-O3`: 11 → 11.
>
> **Minimality.**
>
> - No deletion of one or two gates, and no gate replacement, keeps the separation.
> - The construction needs 11 Ts:
>   - producing `π · y x₁ x₂` takes four Ts per block (the parities 1, x₁, x₂, x₁ ⊕ x₂);
>   - a second block has to uncompute the first;
>   - the middle rotation must be a T, since the pass sees through a Clifford and would pair
>     the two blocks.

```lean
theorem sfCandidate_anti (hd : d ≤ d') (mid : List (CTGate n)) (qk : Fin n) :
    sfCandidate d' mid qk ⊆ sfCandidate d mid qk
```

> **In words.** A higher degree proves more facts, so its concretization only shrinks.

The same incomparability holds for the folding decision itself, on the full circuit
`t q1; cx q0,q1; h q0; h q1; cx q1,q0; h q1; tdg q1` (`Counterexample.lean`):

```lean
theorem cex_pauli :
    PauliAbs.run cexMid (.one (zString {1} false)) = .one (zString {1} false)

theorem cex_not_foldable (d : ℕ) (s : Bool) {C : State 2} {f' g' : Poly}
    (htr : Trace d (run cex) (predAfter ([] : List (CTGate 2)) 1) (predAfter₂ [] cexMid 1) C f' g') :
    ¬ ∀ ν, ev ν g' = (s != ev ν f')
```

> **In words.** The Pauli domain carries the T's axis unchanged to the T†, so the pair folds
> to the identity. No StateFold trace, at any degree, makes their predicates equal or
> complementary.

## 5. Strengthen: folding modulo witnessed constraints

Amy and Lunderville's Algorithm 2 (*Strengthen*) adds constraints witnessed by
interference. Suppose a hidden variable `y` appears in no wire, and the phase is `Φ + y·P`
with `Φ` and `P` free of `y`. Then

```
∑_y (-1)^{y P} = 2 [P = 0],
```

so every path with `P = 1` cancels, even when `P = 0` cannot be solved as an HH
substitution. Predicates are then compared only on the paths where every witnessed `P`
vanishes.

Over `F₂` with the field equations `x² = x`, vanishing on this set is ideal membership. So,
*for the constraints it collects*, this models reduction modulo a Gröbner basis exactly. It
collects fewer constraints than Algorithm 2: only witnesses on the final state, avoiding
every witness variable. The degree bound `d` does not apply to them.

`ProvesS` is `Proves` extended with a list of witnesses on the trace's final state. The
two predicates need only agree where every witnessed constraint vanishes.
`sfCandidateS` is the resulting concretization.

```lean
theorem denote_witnesses (hA : A.WF) (ws : List (ℕ × Poly))
    (hws : ∀ p ∈ ws, Witness A W p.1 p.2) (hδ : Avoids δ W) (hδ' : Avoids δ' W)
    (hag : ∀ ν, (∀ p ∈ ws, ev ν p.2 = false) → ep (δ ν) = ep (δ' ν)) :
    (A.shift δ).denote = (A.shift δ').denote
```

> **In words.** A path sum cannot distinguish two phase changes that agree wherever the
> witnessed constraints hold, provided neither depends on a witness variable.
>
> The proof pairs the `y = 0` and `y = 1` paths, which gives the factor `2 [P = 0]`. It
> handles one witness at a time.

```lean
theorem provesS_sound (hp : ProvesS d (run mid) (BoolPolynomial.var qk.val) q s) :
    axis mid qk = (zString {q} s).toMatrix

theorem sfCandidateS_subset (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) :
    sfCandidateS d mid qk ⊆ sfCandidate d mid qk

theorem gammaS_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      sfCandidateS d mid qk ⊂ pauliCandidate mid qk) ∧
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      pauliCandidate mid qk ⊂ sfCandidateS d mid qk)
```

> **In words.**
>
> - Facts proved with Strengthen are true.
> - Strengthen refines plain StateFold: it proves at least as much, so its concretization is
>   smaller.
> - Even so, it remains incomparable with PauliFold, on the same two segments as in §4.
>   Every path variable of the second segment appears in a wire, so there are no witnesses
>   there either.
>
> `sfCandidateS_anti` and `axis_mem_bothS` carry over monotonicity in `d` and the shared
> true axis.

## Modeling notes

- **Fidelity to Feynman.** Feynman's `-statefold d` on circuits
  (`Feynman/Optimization/StateFold.hs`) uses exactly two reduction rules, HH and ω. It
  applies degree-1 substitutions first, then substitutions of degree `≤ d`.
  - Its ideal-based reduction for circuits (`matchI`, `reduceAll`) is dead or commented out
    in the current code.
  - In the May 2024 version it ran only for unbounded `d`.
  - `Trace` allows any order of HH and ω steps, so it covers Feynman's fixed strategy.
  - The rewrite rule (E), which drops a variable that appears nowhere, never changes a
    predicate, so it is omitted.
- **Fidelity to the paper.** Strengthen collects constraints after every rewrite. Here
  witnesses are taken on the final state of the HH/ω trace, and no constraint may mention a
  witness variable.
  - Carrying intermediate constraints soundly through later substitutions `z := R` would
    also need the equations `z = R`, and a constraint that mentions `z` then breaks the
    pairing argument.
  - The paper's Example 23 is still covered: its `y₄` witness, which folds `ℓ` with `ℓ''`,
    is on the final path sum.
- **Candidates.** In §4 and §5 a candidate starts at the beginning of `mid` with predicate
  `x_{qₖ}`, matching the Pauli and phase domains. `fold_sound` covers arbitrary circuit
  context, `pre` and `post`.
- **Gate set.** H, X, Z, S, S†, T, T†, CX, CZ. There is no CCZ, so strictness of Strengthen
  over plain StateFold is not formalized: it needs a phase term `y·z·w`.

## Proof notes

- **Containment** (Section 10 of
  [`docs/pauli-axis-abstract-interpretation.md`](../docs/pauli-axis-abstract-interpretation.md)).
  The invariant `Inv` has two parts:
  - *freshness:* every variable in the wires and the label is below the fresh counter;
  - *agreement:* every wire set that expresses the label's linear part has, as its signed
    Z-string, exactly the Pauli state.

  Each gate preserves it:
  - diagonal gates fix Z-strings and change no wire;
  - X flips a constant and the matching sign;
  - CX changes the expressing wire set as it maps the Z-string;
  - H's fresh variable forces an expressing wire set to avoid the H's qubit.
- **Pauli soundness.** Three layers:
  1. conjugation lemmas for `embed1`, permutation matrices (CX), and phase matrices (CZ);
  2. finite letter identities: 2 × 2 cases for H, X, Z, S, S†, and all 256 two-letter
     cases for CX and CZ;
  3. assembly by induction over the gate list.
- **StateFold exactness.** Each gate's transformer is proved exact on the path-sum normal
  form. The H case introduces a variable and uses `invSqrt2 ^ card · ∑`.

## Files

| file | contents |
|---|---|
| `PauliFold/Circuit.lean` | Clifford+T gates `CTGate n`; operators `Op n`; the concrete transformer `A ↦ U A U†` |
| `PauliFold/Pauli.lean` | signed Pauli strings, Clifford conjugation tables, the domain `PauliAbs`, its transformer and `gamma` |
| `PauliFold/Phase.lean` | affine forms over path variables, the domain `PhaseAbs`, its transformer and `gamma` |
| `PauliFold/Containment.lean` | the invariant, `containment`, `containment_strict` |
| `PauliFold/Soundness.lean` | matrix-level soundness of the Pauli transformer; `candidate_sound` |
| `PauliFold/StateFold/Basic.lean` | path-sum states, their denotation, the gate transformers |
| `PauliFold/StateFold/Exact.lean` | exactness of every transformer; `denote_run` |
| `PauliFold/StateFold/Reduce.lean` | the HH and ω rules and their soundness |
| `PauliFold/StateFold/Fold.lean` | traces, the shifting lemma, `fold_sound`, `Trace.mono` |
| `PauliFold/StateFold/Counterexample.lean` | the 7-gate fold that StateFold misses at every degree |
| `PauliFold/StateFold/Incomparable.lean` | StateFold's concretization, `proves_sound`, `gamma_incomparable` |
| `PauliFold/StateFold/DegreeTwo.lean` | a 21-gate, 11-T circuit where degree-2 StateFold folds and PauliFoldRand does not, end to end |
| `PauliFold/Tableau.lean` | the tableau abstraction of the Rust pass: phased Pauli strings, their product, row updates, the invariant, and soundness |
| `PauliFold/TableauContainment.lean` | the tableau refines the Pauli and phase domains (strictly) |
| `PauliFold/Affine.lean` | the affine domain as a transition relation; soundness; the tableau and Pauli domains refine it |
| `PauliFold/Join.lean` | partial tableaux (unknown rows), the per-row join, programs with branches and loops, soundness via completions |
| `PauliFold/Disjunctive.lean` | the bounded disjunctive domain; soundness; the swap loop |
| `PauliFold/Alpha.lean` | Amy and Lunderville's `α` (Definition 17) and Algorithm 1; soundness, monotonicity, incomparability with Pauli |
| `PauliFold/FoldRule.lean` | the fold theorem: a segment fact licenses sliding and merging a rotation |
| `PauliFold/Domains.lean` | the domains as abstractions of circuits: Pauli facts, soundness, refinement and incomparability |
| `PauliFold/StateFold/Strengthen.lean` | witnesses, `provesS_sound`, `sfCandidateS_subset`, `gammaS_incomparable` |

The three circuits are in `circuits/`, as OpenQASM:

- `paulifold-beats-statefold.qasm`;
- `statefold-beats-paulifold.qasm`;
- `domain-only-h-t-tdg-h.qasm`.

## Building

The project shares `../lean`'s Mathlib checkout. After building `../lean` once:

```sh
mkdir -p .lake && ln -sfn ../../lean/.lake/packages .lake/packages
lake build
```
