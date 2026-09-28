# pauli-fold (Lean)

A Lean 4 formalization comparing three abstract domains for **rotation folding**:
merging a T or T† rotation into a later rotation on the same Pauli axis.

| domain | used by | abstract value for a candidate rotation |
|---|---|---|
| **Phase** (symbolic phase folding) | `PhaseFoldRand`; the symbolic variant of *Linear time T gate optimization* | affine Boolean forms on the wires, plus the candidate's label |
| **Pauli** (Pauli-axis folding) | `PauliFoldRand` | one signed Pauli string, or `⊤` |
| **StateFold `d`** (+ **Strengthen**) | Feynman `-statefold d`; Amy and Lunderville, POPL 2025 | an exact path sum, reduced by the HH rule at degree `≤ d` and the ω rule (Strengthen adds constraints witnessed by interference) |

All three analyses answer the same question. A candidate rotation acts on qubit `qₖ` at the
start of a Clifford+T segment `mid`. What is its axis `Z_{qₖ}` after the segment? The
concrete answer is the operator

```
axis mid qₖ  =  U Z_{qₖ} U†        (U = the unitary of mid)
```

Each domain's concretization `γ` is a set of operators that must contain this axis.
Smaller is more precise.

The main results:

```
            Pauli  ⊊  Phase                   (containment, containment_strict)

  StateFold d (+ Strengthen)  ⋈  Pauli         (gamma_incomparable, gammaS_incomparable)
           Strengthen  ⊆  StateFold d          (sfCandidateS_subset)
        StateFold d'   ⊆  StateFold d, d ≤ d'  (sfCandidate_anti)
```

Here `⋈` means that the two concretizations always overlap, since both contain the true
axis, and that neither is contained in the other in general. Every theorem depends only on
the standard axioms `propext`, `Classical.choice`, and `Quot.sound`.

The concrete semantics is that of [`../lean`](../lean) (`tzap-lean`): operators are its
`Density n` matrices, and gates denote `TzapLean.gateUnitary`.

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

Over `F₂` with the field equations `x² = x`, vanishing on this set is ideal membership. So
this models reduction modulo a Gröbner basis exactly. Feynman's implementation, which uses
non-canonical multivariate division, proves no more.

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
