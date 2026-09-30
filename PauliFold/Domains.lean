import PauliFold.StateFold.DegreeTwo
import PauliFold.Tableau

/-!
# Comparing the domains as abstractions of circuits

The earlier files compare the analyses on one candidate rotation's axis. This file states
the comparisons as abstractions of whole circuits, independent of folding.

**Concrete semantics.** A circuit `gs` denotes its unitary `circ gs`.

**Facts.** A *Pauli fact* `(P, P')` holds for a unitary `U` when `U P U† = P'`, as matrices.
Each domain abstracts a circuit by the set of facts its analysis derives. Its
concretization is the set of unitaries satisfying all of them:

  `γ(F) = { U | ∀ (P, P') ∈ F, U P U† = P' }`.

More facts mean a smaller set, and so a more precise abstraction.

| domain | facts about `gs` |
|---|---|
| Pauli | `(P, P')` whenever the Pauli transformer maps `P` to `P'` (not `⊤`) |
| Phase | `(Z_q, ±Z_Q)` whenever the wires in `Q` sum to the input variable `x_q` |
| StateFold `d` | `(Z_q, ±Z_{q'})` whenever a degree-`d` trace proves the wire predicate of `q'` equal or complementary to `x_q` |
| Strengthen `d` | the same, with witnessed constraints |

Results:

- **Soundness.** Every domain contains the circuit's own unitary.
- **Pauli ⊆ Phase.** `pauli_refines_phase`, strict by `pauli_refines_phase_strict`.
- **Pauli and StateFold are incomparable** (`pauli_stateFold_incomparable`), for every `d ≥ 1`.
  In each direction a *Z-to-Z* fact that one domain proves is not implied by the other:
  - `cx q0,q1; h q0; h q1; cx q1,q0; h q1`: Pauli proves `Z₁ ↦ Z₁`, and StateFold proves
    nothing at any degree.
  - `h q0; t q0; tdg q0; h q0`: StateFold proves `Z₀ ↦ Z₀`, and the Pauli facts allow `−Z₀`.
  - The 21-gate relative-Toffoli circuit `midE`: the same, for `d ≥ 2`.
- **Strengthen ⊆ StateFold.** `strengthen_refines_stateFold`. Strengthen is still
  incomparable with Pauli.
- **Monotone in the degree.** `stateFold_anti`: a higher degree only shrinks the set.

**Scope.** `sfFacts` and `sfFactsS` are restricted models of Amy and Lunderville's analyses:
they keep only `Z_q ↦ ±Z_{q'}` facts, not general affine or polynomial transition
relations. For example, `1 ∈ sfGamma d [cx 0 1]`, while `affGamma [cx 0 1]` excludes it. The
"StateFold proves more" direction of the incomparability transfers to the full domain,
because its concretization is smaller. The "Pauli proves more" direction is proved against the paper's
path-sum abstraction `α` in `Alpha.lean` (`alpha_pauli_incomparable`).
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Facts and their concretization -/

/-- A Pauli fact: a Pauli string and its claimed image under conjugation. -/
abbrev Fact (n : ℕ) := SPauli n × SPauli n

/-- `U P U† = P'`. -/
def Holds (U : Density n) (f : Fact n) : Prop := conj U f.1.toMatrix = f.2.toMatrix

/-- The unitaries satisfying every fact in `F`. -/
def gammaF (F : Set (Fact n)) : Set (Density n) := {U | ∀ f ∈ F, Holds U f}

/-- More facts, fewer unitaries. -/
theorem gammaF_anti {F F' : Set (Fact n)} (h : F ⊆ F') : gammaF F' ⊆ gammaF F :=
  fun _ hU f hf => hU f (h hf)

/-! ## The four domains -/

/-- The Pauli domain's facts: every Pauli string it tracks through `gs` without reaching `⊤`. -/
def pauliFacts (gs : List (CTGate n)) : Set (Fact n) :=
  {f | PauliAbs.run gs (.one f.1) = .one f.2}

/-- The phase domain's facts: `Z_q ↦ ±Z_Q` whenever the wires in `Q` sum to `x_q`. -/
def phaseFacts (gs : List (CTGate n)) : Set (Fact n) :=
  {f | ∃ q Q, f.1 = zString {q} false ∧
    (PhaseAbs.run gs (PhaseAbs.init n q)).Expresses Q ∧
    f.2 = zString Q ((PhaseAbs.run gs (PhaseAbs.init n q)).sign Q)}

/-- StateFold's facts: `Z_q ↦ (-1)^s Z_{q'}` whenever a degree-`d` trace proves it. -/
def sfFacts (d : ℕ) (gs : List (CTGate n)) : Set (Fact n) :=
  {f | ∃ q q' s, f.1 = zString {q} false ∧ f.2 = zString {q'} s ∧
    SF.Proves d (SF.run gs) (BoolPolynomial.var q.val) q' s}

/-- Strengthen's facts: the same, with witnessed constraints. -/
def sfFactsS (d : ℕ) (gs : List (CTGate n)) : Set (Fact n) :=
  {f | ∃ q q' s, f.1 = zString {q} false ∧ f.2 = zString {q'} s ∧
    SF.ProvesS d (SF.run gs) (BoolPolynomial.var q.val) q' s}

abbrev pauliGamma (gs : List (CTGate n)) : Set (Density n) := gammaF (pauliFacts gs)
abbrev phaseGamma (gs : List (CTGate n)) : Set (Density n) := gammaF (phaseFacts gs)
abbrev sfGamma (d : ℕ) (gs : List (CTGate n)) : Set (Density n) := gammaF (sfFacts d gs)
abbrev sfGammaS (d : ℕ) (gs : List (CTGate n)) : Set (Density n) := gammaF (sfFactsS d gs)

/-! ## Soundness -/

/-- The Pauli domain is sound: the circuit's unitary satisfies every Pauli fact. -/
theorem pauli_sound (gs : List (CTGate n)) : circ gs ∈ pauliGamma gs := fun f hf => by
  have h := run_sound gs (.one f.1) f.1.toMatrix (Set.mem_singleton _)
  rw [show PauliAbs.run gs (.one f.1) = .one f.2 from hf] at h
  rw [Holds, ← SF.concreteRun_eq]
  exact h

/-- Every phase fact is a Pauli fact: the containment invariant says that each wire set
expressing the label is, as a Z-string, the Pauli state. -/
theorem phaseFacts_subset (gs : List (CTGate n)) : phaseFacts gs ⊆ pauliFacts gs := by
  rintro ⟨P, P'⟩ ⟨q, Q, rfl, hQ, rfl⟩
  exact (inv_run gs (inv_init q)).2 Q hQ

/-- **Pauli refines phase:** on every circuit, the Pauli concretization is contained in the
phase one. -/
theorem pauli_refines_phase (gs : List (CTGate n)) : pauliGamma gs ⊆ phaseGamma gs :=
  gammaF_anti (phaseFacts_subset gs)

/-- The phase domain is sound. -/
theorem phase_sound (gs : List (CTGate n)) : circ gs ∈ phaseGamma gs :=
  pauli_refines_phase gs (pauli_sound gs)

/-- StateFold is sound. -/
theorem sf_sound (d : ℕ) (gs : List (CTGate n)) : circ gs ∈ sfGamma d gs := by
  rintro ⟨P, P'⟩ ⟨q, q', s, rfl, rfl, hp⟩
  have h := SF.proves_sound hp
  unfold SF.axis at h
  rw [Holds, ← SF.concreteRun_eq]
  exact h

/-- Strengthen is sound. -/
theorem sfS_sound (d : ℕ) (gs : List (CTGate n)) : circ gs ∈ sfGammaS d gs := by
  rintro ⟨P, P'⟩ ⟨q, q', s, rfl, rfl, hp⟩
  have h := SF.provesS_sound hp
  unfold SF.axis at h
  rw [Holds, ← SF.concreteRun_eq]
  exact h

/-- A higher degree proves more facts, so its concretization is smaller. -/
theorem stateFold_anti {d d' : ℕ} (hd : d ≤ d') (gs : List (CTGate n)) :
    sfGamma d' gs ⊆ sfGamma d gs :=
  gammaF_anti fun _ ⟨q, q', s, h1, h2, C, f', g', htr, heq⟩ =>
    ⟨q, q', s, h1, h2, C, f', g', htr.mono hd, heq⟩

/-- **Strengthen refines StateFold.** -/
theorem strengthen_refines_stateFold (d : ℕ) (gs : List (CTGate n)) :
    sfGammaS d gs ⊆ sfGamma d gs :=
  gammaF_anti fun _ ⟨q, q', s, h1, h2, hp⟩ => ⟨q, q', s, h1, h2, hp.toS⟩

/-! ## Matrix lemmas -/

theorem conj_mul (A B M : Density n) : conj (A * B) M = conj A (conj B M) := by
  unfold conj
  rw [Matrix.conjTranspose_mul]
  simp only [Matrix.mul_assoc]

theorem conj_neg (U M : Density n) : conj U (-M) = -conj U M := by
  unfold conj
  rw [Matrix.mul_neg, Matrix.neg_mul]

theorem toMatrix_flip (P : SPauli n) : (⟨!P.neg, P.letter⟩ : SPauli n).toMatrix = -P.toMatrix := by
  funext o i
  simp only [SPauli.toMatrix, Matrix.neg_apply]
  cases P.neg <;> simp [SPauli.sgn]

/-- `-Z_Q ≠ Z_Q`: the two differ at the all-zero entry. -/
theorem zString_sign_ne (Q : Finset (Fin n)) : (zString Q true).toMatrix ≠ (zString Q false).toMatrix := by
  have hp : ∀ b, (zString Q b).toMatrix (fun _ => false) (fun _ => false) = SPauli.sgn b := by
    intro b
    unfold SPauli.toMatrix
    rw [Finset.prod_eq_one fun r _ => by simp only [zString]; split_ifs <;> rfl, mul_one]
    rfl
  intro h
  have h' := congrFun (congrFun h fun _ => false) fun _ => false
  rw [hp, hp] at h'
  simp [SPauli.sgn] at h'
  norm_num at h'

/-- X on `q` fixes a Pauli with `I` or `X` at `q`. -/
theorem conj_x_fix (q : Fin n) (P : SPauli n) (hP : P.letter q = .I ∨ P.letter q = .X) :
    conj (gateUnitary n (CTGate.x q).toGate) P.toMatrix = P.toMatrix := by
  have h := step_sound (.x q) (.one P) P.toMatrix (Set.mem_singleton _)
  simp only [PauliAbs.step, SPauli.conjClifford, PauliAbs.gamma, Set.mem_singleton_iff] at h
  rw [show SPauli.apply1 Letter.conjX q P = P from
    apply1_of_fix (by rcases hP with h | h <;> rw [h] <;> rfl)] at h
  exact h

/-- X on `q` flips `Z_q`. -/
theorem conj_x_z (q : Fin n) :
    conj (gateUnitary n (CTGate.x q).toGate) (zString {q} false).toMatrix =
      -(zString {q} false).toMatrix := by
  have h := step_sound (.x q) (.one (zString {q} false)) _ (Set.mem_singleton _)
  simp only [PauliAbs.step, SPauli.conjClifford, PauliAbs.gamma, Set.mem_singleton_iff] at h
  rw [conjX_zString] at h
  have e : zString {q} (xor false (decide (q ∈ ({q} : Finset (Fin n))))) =
      ⟨!(zString {q} false).neg, (zString {q} false).letter⟩ := by simp [zString]
  change (CTGate.x q).concrete _ = _
  rw [h, e, toMatrix_flip]

/-! ## Incomparability of Pauli and StateFold -/

/-- If a circuit starts with `h q` and then a T or T† on `q`, every Pauli the Pauli domain
tracks through it has `I` or `X` at `q`. -/
theorem pauliFacts_letter {q : Fin n} {g : CTGate n} (hg : g = .t q ∨ g = .tdg q)
    (rest : List (CTGate n)) {f : Fact n} (hf : f ∈ pauliFacts (.h q :: g :: rest)) :
    f.1.letter q = .I ∨ f.1.letter q = .X := by
  have htop : ∀ gs : List (CTGate n), PauliAbs.run gs .top = .top := fun gs => by
    induction gs with
    | nil => rfl
    | cons g gs ih => exact ih
  by_contra hne
  have hZY : f.1.letter q = .Y ∨ f.1.letter q = .Z := by
    cases h : f.1.letter q <;> simp_all
  have hstep : PauliAbs.step g (PauliAbs.step (.h q) (.one f.1)) = .top := by
    rcases hg with rfl | rfl <;>
      simp only [PauliAbs.step, SPauli.conjClifford, SPauli.apply1, Function.update_self] <;>
      rcases hZY with h | h <;> simp [h, Letter.conjH, Letter.anticommutesZ]
  have := hf
  simp only [pauliFacts, Set.mem_ofPred_eq, PauliAbs.run, List.foldl] at this
  rw [hstep] at this
  rw [show List.foldl (fun a g => PauliAbs.step g a) PauliAbs.top rest = PauliAbs.top from
    htop rest] at this
  cases this

/-- On a circuit starting with `h q; t q` or `h q; tdg q`, the circuit followed by X on `q`
satisfies every Pauli fact but sends `Z_q` where the circuit does, negated. -/
theorem pauli_allows_flip {q : Fin n} {g : CTGate n} (hg : g = .t q ∨ g = .tdg q)
    (rest : List (CTGate n)) :
    circ (.h q :: g :: rest) * gateUnitary n (CTGate.x q).toGate ∈ pauliGamma (.h q :: g :: rest) :=
  fun f hf => by
    rw [Holds, conj_mul, conj_x_fix q f.1 (pauliFacts_letter hg rest hf)]
    exact pauli_sound _ f hf

/-- A StateFold fact `Z_q ↦ Z_q` on such a circuit is not implied by the Pauli facts. -/
theorem sf_fact_not_pauli {d : ℕ} {q : Fin n} {g : CTGate n} (hg : g = .t q ∨ g = .tdg q)
    (rest : List (CTGate n))
    (hp : SF.Proves d (SF.run (.h q :: g :: rest)) (BoolPolynomial.var q.val) q false) :
    (zString {q} false, zString {q} false) ∈ sfFacts d (.h q :: g :: rest) ∧
    ∃ U ∈ pauliGamma (.h q :: g :: rest), ¬ Holds U (zString {q} false, zString {q} false) := by
  have hmem : (zString {q} false, zString {q} false) ∈ sfFacts d (.h q :: g :: rest) :=
    ⟨q, q, false, rfl, rfl, hp⟩
  refine ⟨hmem, _, pauli_allows_flip hg rest, fun h => ?_⟩
  have htrue := sf_sound d (.h q :: g :: rest) _ hmem
  simp only [Holds] at h htrue
  rw [conj_mul, conj_x_z, conj_neg, htrue] at h
  exact zString_sign_ne {q} (by rw [show zString {q} true =
    ⟨!(zString {q} false).neg, (zString {q} false).letter⟩ from rfl, toMatrix_flip]; exact h)

/-- On `cx q0,q1; h q0; h q1; cx q1,q0; h q1`, the Pauli fact `Z₁ ↦ Z₁` is not implied by
StateFold's facts, since StateFold proves none at any degree. -/
theorem x_not_holds (q : Fin n) :
    ¬ Holds (gateUnitary n (CTGate.x q).toGate) (zString {q} false, zString {q} false) := by
  intro h
  simp only [Holds] at h
  rw [conj_x_z] at h
  exact zString_sign_ne {q} (by rw [show zString {q} true =
    ⟨!(zString {q} false).neg, (zString {q} false).letter⟩ from rfl, toMatrix_flip]; exact h)

theorem pauli_fact_not_sf (d : ℕ) :
    ((zString {(1 : Fin 2)} false, zString {(1 : Fin 2)} false) : Fact 2) ∈ pauliFacts SF.cexMid ∧
    ∃ U ∈ sfGamma d SF.cexMid, ¬ Holds U (zString {(1 : Fin 2)} false, zString {(1 : Fin 2)} false) := by
  refine ⟨SF.cex_pauli, gateUnitary 2 (CTGate.x 1).toGate, ?_, x_not_holds 1⟩
  rintro _ ⟨q, q', s, -, -, hp⟩
  exact absurd hp (SF.mid_not_proves_var d q.val (by omega) (by omega) q' s)

theorem midA_eq : SF.midA = .h 0 :: .t 0 :: [.tdg 0, .h 0] := rfl

theorem midE_eq : SF.midE = .h 0 :: .tdg 0 :: (SF.midE.drop 2) := rfl

/-- **Pauli and StateFold are incomparable abstractions of circuits,** at every degree
`d ≥ 1`. On `cexMid`, Pauli proves a Z-fact StateFold cannot; on `h; t; tdg; h`, StateFold
proves a Z-fact Pauli's facts do not imply. -/
theorem pauli_stateFold_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (gs : List (CTGate n)), ¬ sfGamma d gs ⊆ pauliGamma gs) ∧
    (∃ (n : ℕ) (gs : List (CTGate n)), ¬ pauliGamma gs ⊆ sfGamma d gs) := by
  refine ⟨⟨2, SF.cexMid, fun h => ?_⟩, ⟨1, SF.midA, fun h => ?_⟩⟩
  · obtain ⟨hf, U, hU, hn⟩ := pauli_fact_not_sf d
    exact hn (h hU _ hf)
  · rw [midA_eq] at h
    obtain ⟨hf, U, hU, hn⟩ := sf_fact_not_pauli (d := d) (Or.inl rfl) [.tdg 0, .h 0]
      (by rw [← midA_eq]; exact SF.midA_proves d hd)
    exact hn (h hU _ hf)

/-- The same on the 21-gate relative-Toffoli circuit, where the PauliFoldRand pass also
misses the fold: for `d ≥ 2`, StateFold proves `Z₀ ↦ Z₀` and the Pauli facts do not imply
it. -/
theorem midE_sf_not_pauli (d : ℕ) (hd : 2 ≤ d) : ¬ pauliGamma SF.midE ⊆ sfGamma d SF.midE := by
  intro h
  rw [midE_eq] at h
  obtain ⟨hf, U, hU, hn⟩ := sf_fact_not_pauli (d := d) (Or.inr rfl) (SF.midE.drop 2) (by
    rw [← midE_eq]
    obtain ⟨C, f', g', htr, heq⟩ := SF.midE_proves
    exact ⟨C, f', g', htr.mono hd, heq⟩)
  exact hn (h hU _ hf)

/-- Strengthen is also incomparable with Pauli. -/
theorem pauli_strengthen_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (gs : List (CTGate n)), ¬ sfGammaS d gs ⊆ pauliGamma gs) ∧
    (∃ (n : ℕ) (gs : List (CTGate n)), ¬ pauliGamma gs ⊆ sfGammaS d gs) := by
  refine ⟨⟨2, SF.cexMid, fun h => ?_⟩, ⟨1, SF.midA, fun h => ?_⟩⟩
  · have hU : gateUnitary 2 (CTGate.x (1 : Fin 2)).toGate ∈ sfGammaS d SF.cexMid := by
      rintro _ ⟨q, q', s, -, -, hp⟩
      exact absurd hp (SF.mid_not_provesS_var d q.val (by omega) (by omega) q' s)
    exact x_not_holds 1 (h hU _ (pauli_fact_not_sf d).1)
  · have hp := SF.midA_proves d hd
    rw [midA_eq] at h hp
    obtain ⟨-, U, hU, hn⟩ := sf_fact_not_pauli (d := d) (Or.inl rfl) [.tdg 0, .h 0] hp
    exact hn (h hU _ ⟨0, 0, false, rfl, rfl, hp.toS⟩)

/-! ## Strictness of Pauli over phase -/

/-- On `h; h`, Pauli proves `Z ↦ Z` and the phase domain proves nothing. -/
theorem pauli_refines_phase_strict : pauliGamma hh ⊂ phaseGamma hh := by
  refine ⟨pauli_refines_phase hh, fun h => ?_⟩
  -- H itself satisfies every phase fact (there are none) but sends Z to X.
  have hnone : phaseFacts hh = ∅ := by
    ext ⟨P, P'⟩
    simp only [phaseFacts, Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false]
    rintro ⟨q, Q, -, hQ, -⟩
    have hq : q = 0 := Subsingleton.elim _ _
    subst hq
    have := phase_hh_gamma
    rw [PhaseAbs.gamma, if_pos ⟨Q, hQ⟩] at this
    have h0 : (0 : Op 1) ∈ ({M | ∃ Q, (PhaseAbs.run hh (PhaseAbs.init 1 0)).Expresses Q ∧
        M = (zString Q ((PhaseAbs.run hh (PhaseAbs.init 1 0)).sign Q)).toMatrix} : Set (Op 1)) :=
      this ▸ Set.mem_univ _
    obtain ⟨Q', -, hQ'⟩ := h0
    exact SF.zString_toMatrix_ne_zero _ _ hQ'.symm
  have hH : circ [CTGate.h (0 : Fin 1)] ∈ phaseGamma hh := by
    rw [phaseGamma, hnone]; exact fun _ hf => absurd hf (Set.notMem_empty _)
  have hZ : (zString {(0 : Fin 1)} false, zString {0} false) ∈ pauliFacts hh := by
    simp only [pauliFacts, Set.mem_ofPred_eq, PauliAbs.run, hh, List.foldl, PauliAbs.step,
      SPauli.conjClifford]
    congr 1
    simp only [zString, SPauli.apply1]
    congr 1
    funext i; fin_cases i; rfl
  have hX : Holds (circ [CTGate.h (0 : Fin 1)])
      (zString {(0 : Fin 1)} false, ⟨false, fun _ => .X⟩) :=
    pauli_sound [.h 0] _ (by
      simp only [pauliFacts, Set.mem_ofPred_eq, PauliAbs.run, List.foldl, PauliAbs.step,
        SPauli.conjClifford]
      congr 1
      simp only [zString, SPauli.apply1]
      congr 1
      funext i; fin_cases i; rfl)
  have hZZ := h hH _ hZ
  simp only [Holds] at hX hZZ
  rw [hX] at hZZ
  have := congrFun (congrFun hZZ fun _ => false) fun _ => false
  simp [SPauli.toMatrix, zString, Letter.mat, SPauli.sgn] at this

end

end PauliFold
