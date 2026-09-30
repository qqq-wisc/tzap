# Amy–Lunderville control-flow programs

The OpenQASM 3 programs with loops, branches, resets, and measurements from Amy and
Lunderville, *Linear and non-linear relational analyses for quantum program optimization*
(POPL 2025, [arXiv:2410.23493](https://arxiv.org/abs/2410.23493)).

- **Source.** Every file except `rus-ccz.qasm` is copied verbatim from the Feynman
  repository's `benchmarks/qasm3`, the paper's artifact (checked out at `4c4c488`).
- **`rus-ccz.qasm`** reconstructs Figure 2b: the RUS loop with each Toffoli written as a CCZ
  conjugated by Hadamards.
- **Non-determinism.** Loops written `while (true)` and branches written `if (true)` stand
  for the paper's `while ★` and `if ★`. The analyses treat them as non-deterministic.
- **Supported input.** tzap reads OpenQASM 2 only, so these programs run through Feynman's
  QASM 3 frontend: `feynopt -qasm3 <passes> file.qasm`.

## Programs

In the pseudocode, `;` is sequencing, `while ★ { … }` is a non-deterministic loop, and
`if ★ { … }` a non-deterministic branch.

| file | paper | program | what it tests |
|---|---|---|---|
| `rus.qasm` | Fig. 2a, Table 1 "RUS" | `h ψ; t ψ; while ab≠00 { reset a,b; h a; h b; ccx a,b,ψ; s ψ; ccx a,b,ψ; z ψ; h a; h b; measure a,b }; tdg ψ; h ψ` | repeat-until-success. The body is diagonal on `ψ`, so the `t` before the loop cancels the `tdg` after it. |
| `rus-ccz.qasm` | Fig. 2b | the same loop, with `ccx = h·ccz·h` on `ψ` | the same invariant, which now needs interference: `z' = z` follows only from the path sum of `h·ccz·h` |
| `loop-swap.qasm` | Fig. 3, Table 1 "Loop-swap" | `cx a,b; t b; cx a,b; while ★ { swap a,b }; cx a,b; tdg b; cx a,b` | the T gates act on `a ⊕ b`, and the swap loop preserves `x' ⊕ y' = x ⊕ y` |
| `loop-nonlinear.qasm` | Fig. 11, Table 1 "Loop-nonlinear" | `reset c; x a; ccx a,b,c; t c; ccx a,b,c; x a; while ★ { cx a,b }; x a; ccx a,b,c; tdg c; ccx a,b,c; x a`, with `ccx` and `ccz` defined as gates over Clifford+T | the non-linear invariant `(1⊕x₀')x₁' = (1⊕x₀)x₁` |
| `loop-simple.qasm` | Table 1 "Loop-simple" | `t a; while ★ { cx a,b }; tdg a` | the control of a CNOT is invariant |
| `loop-h.qasm` | Table 1 "Loop-h" | `t b; while ★ { h a }; tdg b` | a loop on another qubit leaves `b` invariant |
| `loop-nested.qasm` | Table 1 "Loop-nested" | `reset a; t b; while ★ { t a; while ★ { x b } }; tdg b` | nested loops. The inner `x` flips `b` an unknown number of times, so the outer T and T† cannot merge. |
| `loop-null.qasm` | Table 1 "Loop-null" | `t b; reset a; while ★ { t b; t a }; tdg b` | phase gates inside a loop. `t a` acts on a qubit reset to 0 before the loop. |
| `if-simple.qasm` | Table 1 "If-simple" | `t a; if ★ { cx a,b }; tdg a` | a branch that preserves the control |
| `reset-simple.qasm` | Table 1 "Reset-simple" | `t a; reset a; t a` | a T after a reset acts on 0 |
| `grover.qasm` | Table 1 "Grover" | 64-bit Grover search. The ~10⁹ iterations are modeled as `while (true)`, and the oracle is a Toffoli ladder over 64 ancillas. | scale: 129 qubits |
| `loop-block.qasm` | not in the paper | `t a; h a; while ★ { t a }; h a; tdg a` | a loop of T gates inside an H frame |
| `loop-cycle.qasm` | not in the paper | `t b; while ★ { swap a,b; swap b,c }; tdg b` | a 3-cycle permutation, with no affine invariant on `b` |

## T-counts

The paper columns are from Table 1:

- PF_Aff is the affine relational analysis.
- PF_Pol is the polynomial-ideal analysis with Strengthen.
- The invariant is the loop summary PF_Pol computes.

The Feynman columns are from `feynopt -qasm3` at `4c4c488`:

- `-phasefold` runs the affine analysis.
- `-statefold d` runs the path-sum analysis, with HH substitutions up to degree `d`
  (unbounded for `0`).

| program | n | T in | PF_Aff (paper) | PF_Pol (paper) | invariant (paper) | `-phasefold` | `-statefold 2` | `-statefold 0` |
|---|---|---|---|---|---|---|---|---|
| RUS | 3 | 16 | 10 | 8 | ⟨z′ ⊕ z⟩ | 10 | 8 | 8 |
| RUS, Fig. 2b | 3 | 16 | – | – | – | 10 | 8 | 8 |
| Loop-swap | 2 | 2 | 0 | 0 | ⟨x′⊕y′⊕x⊕y, x′⊕xy⊕xx′⊕yx′⟩ | 0 | 0 | 0 |
| Loop-nonlinear | 3 | 30 | 18 | 0 | ⟨x′⊕x, z′⊕z, y′⊕y⊕xy⊕xy′⟩ | 18 | 12 | 12 |
| Loop-simple | 2 | 2 | 0 | 0 | ⟨x′⊕x, y⊕y′⊕xy⊕xy′⟩ | 0 | 0 | 0 |
| Loop-h | 2 | 2 | 0 | 0 | ⟨y′⊕y⟩ | 0 | 0 | 0 |
| Loop-nested | 2 | 3 | 2 | 2 | ⟨x′⊕x⟩, ⟨x′⊕x⟩ | 2 | 3 | 3 |
| Loop-null | 2 | 4 | 2 | 2 | ⟨x′⊕x, y′⊕y⟩ | 1 | 2 | 2 |
| If-simple | 2 | 2 | 0 | 0 | – | 0 | 0 | 0 |
| Reset-simple | 2 | 2 | 1 | 1 | – | 1 | 1 | 1 |
| Grover | 129 | 1736 × 10⁹ | 1472 × 10⁹ | timeout | – | not run | not run | not run |
| Loop-block | 2 | 3 | – | – | – | 3 | 3 | 3 |
| Loop-cycle | 3 | 2 | – | – | – | 2 | 2 | 2 |

Today's Feynman differs from the paper's numbers in three places:

- **Loop-nonlinear:** 12 Ts with `-statefold`, where the paper reports 0 for PF_Pol.
- **Loop-nested:** 3 with `-statefold`, versus 2 with `-phasefold`.
- **Loop-null:** 1 with `-phasefold`, where the paper reports 2.

`-statefold d` is the tool's current path-sum analysis over programs (`stateAnalysispp`),
not necessarily the exact PF_Pol configuration used for Table 1. The T-count of `rus.qasm`
is 16 because Feynman inlines each `ccx` into a 7-T CCZ.
