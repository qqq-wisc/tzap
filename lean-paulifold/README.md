# Pauli Abstractions for Quantum Circuits: Lean

`PauliFold/Paper/` mechanizes the paper, in namespace `PauliFold.Paper`. `lake build` checks
it. There is no `sorry` or `native_decide`; the results depend only on `propext`,
`Classical.choice` and `Quot.sound`.

Conventions:

- Qubits are 0-indexed: the paper's qubits 1 and 2 are Lean `0` and `1`.
- `Γ.x q` and `Γ.z q` are the left-hand sides of the constraints `Γ(X_q) ↦ X_q` and
  `Γ(Z_q) ↦ Z_q`; `none` is `⋆`. Axes are lists.
- `T` is `R_Z(π/4)`.
- The polynomial abstraction (`alphaGammaC`) is the most precise result of any sequence of
  degree-`d` HH and ω rewrites, over unitaries.

| paper | Lean |
|---|---|
| Defs 3.1–3.7 | `PP`, `Paulis`, `Takes`, `Cliffords`, `rzMat`, `Circuit`, `Circuit.sem`, `Circuit.Equiv` |
| Defs 4.1, 4.2 | `Models`, `ARel`, `Dom` |
| Def 4.4, Lemma 4.5 | `Cons.wp`, `Cons.wp_toMatrix` |
| Def 4.6 | `ARel.F`, `ARel.gamma` |
| Defs 4.8, 4.10, 4.11, 4.13 | `CliffordGate.cons`, `Circuit.alpha`, `Cons.comp`, `ARel.comp` |
| Lemma 4.14, Thm 4.15 | `lemma_4_14_*`, `Circuit.soundness` |
| Thm 5.1, Alg 1, Thm 5.2 | `phase_folding`, `phaseFold`, `phaseFold_correct` |
| Def 6.1, Lemma 6.2 | `alphaMeas`, `alphaReset`, `lemma_6_2` |
| Defs 6.4–6.6, Lemma 6.7 | `Dom`, `ARel.gamma`, `Dom.join`, `lemma_6_7` |
| Defs 6.9–6.11, Lemma 6.12, Thm 6.14 | `Program`, `Program.psem`, `Program.alpha`, `lemma_6_12`, `theorem_6_14` |
| §8.4, bounded disjunction | `Disj.program_sound_impl` |
| Defs 7.1–7.3, Thm 7.4 | `LinState`, `supp`, `gammaL`, `pauli_refines_linear`, `pauli_refines_linear_strict` |
| §7.2, incomparability | `pauli_polynomial_incomparable` |

Examples are `example_<section>_<number>`; the claims of Section 2 are in `Paper/Claims.lean`
and `Paper/Extras.lean`. The rest of `PauliFold/` is an earlier development that `Paper/` builds
on.
