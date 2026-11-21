# Pauli Abstractions for Quantum Circuits: Lean

`PauliFold/Paper/` mechanizes the paper, in namespace `PauliFold.Paper`. `lake build` checks
it. There is no `sorry` or `native_decide`; the results depend only on `propext`,
`Classical.choice` and `Quot.sound`.

Conventions:

- Qubits are 0-indexed: the paper's qubits 1 and 2 are Lean `0` and `1`.
- A `Cons` `Γ` is a set of transition constraints, possibly partial: `Γ.x q` and `Γ.z q` are the
  left-hand sides of the constraints `Γ(X_q) ↦ X_q` and `Γ(Z_q) ↦ Z_q`; `none` means the
  constraint is dropped. Section 3 of the paper uses full tableaux, which are the case
  `Cons.AllKnown` (`Cons.AllKnown.comp`, `Cons.wp_isSome_of_allKnown`, `Circuit.alpha_spec`).
- The paper writes composition left to right, `σ₁ ⨟ σ₂`; Lean writes the same operation right
  to left, `σ₂ ⊚ σ₁` (`Cons.comp Γ₂ Γ₁`, `ARel.comp σ₂ σ₁`), with the same meaning.
- Commutation constraints (axes) are lists; the concretization depends only on their elements.
- `T` is `R_Z(π/4)`, which differs from `diag(1, e^{iπ/4})` by a global phase.
- The polynomial abstraction (`alphaGammaC`) is the most precise result of any sequence of
  degree-`d` HH and ω rewrites, over unitaries.
- Docstrings in the Lean files cite an earlier numbering of the paper; the table below gives
  the current numbering.

| paper | Lean |
|---|---|
| §3.1, Defs 3.1–3.5 | `PP`, `Paulis`, `Takes`, `PP.comm_or_anticomm`, `Cliffords`, `CliffordGate`, `rzMat`, `Circuit`, `Circuit.sem`, `EquivP`, `Circuit.Equiv` |
| Def 3.6 | `Cons`, `Cons.AllKnown`, `Models` |
| Def 3.7, Ex 3.8 | `ARel`, `Cons.Sat`; `Examples.example_4_3`, `example_4_3_model` |
| Def 3.9, Lemma 3.10 | `Cons.wp`, `Cons.wp_isSome_of_allKnown`; `Cons.wp_toMatrix`, `Cons.wp_conj`, `Cons.wp_mul` |
| Def 3.11, Ex 3.12 | `ARel.F`, `ARel.gamma`, `ARel.mem_F_iff`; `Examples.example_4_7*` |
| Def 3.13 (Fig. 2), Ex 3.14 | `CliffordGate.cons`, `Cons.id`, `Circuit.alpha`, `Cons.comp`, `ARel.comp`, `Models.comp`; `Examples.example_4_9_H`, `example_4_9_CX`, `example_4_12`, `example_4_12_T` |
| Fig. 3 | `lemma_4_14_H`, `_S`, `_X`, `_Z`, `_CX`, `_rz` |
| Thm 3.15 | `Circuit.soundness`, `Circuit.takes_of_tracked` |
| Thm 4.1, Alg. 1, Thm 4.2 | `phase_folding`, `phaseFold`, `phaseFold_correct` |
| Fig. 1 | `Figure1.cuts`, `Figure1.cut2_table`, `Figure1.mid_eq_H`, `Figure1.figure1` |
| Def 5.1, Lemma 5.2, Ex 5.3 | `alphaMeas`, `alphaReset`, `lemma_6_2`; `Examples.example_6_3*` |
| Def 5.4, §5.2 (dropped constraints, ⊤), Lemma 5.5, Ex 5.6 | `Cons.join`, `ARel.join`, `Dom.join`, `Dom`, `comp_eq_top`, `Dom.top_join`, `Dom.join_top`; `lemma_6_7`; `Examples.example_6_8*` |
| Defs 5.7–5.9, Lemma 5.10, Ex 5.11, Thm 5.12 | `Program`, `Program.psem`, `Program.alpha`, `lemma_6_12`, `Examples.example_6_13*`, `theorem_6_14` |
| Cor 5.13 (phase folding for programs) | `program_folding`, `psem_seq_congr`, `psem_choice_congr`, `psem_loop_congr`, `loopSimple_folds` (in `Paper/ProgramFolding.lean`) |
| §7.3, bounded disjunction | `Disj.program_sound_impl`, `program_folding_disj` |
| Defs 6.1–6.3, Thm 6.4 | `LinState`, `supp`, `gammaL`, `pauli_refines_linear`, `pauli_refines_linear_strict` |
| §6.2, incomparability | `pauli_polynomial_incomparable`, `polynomial_more_precise`, `polynomial_more_precise_deg2`, `pauli_more_precise` |

Examples are `example_<section>_<number>` in the earlier numbering; the claims of Section 2 are
in `Paper/Claims.lean` and `Paper/Extras.lean`. The rest of `PauliFold/` is an earlier
development that `Paper/` builds on.
