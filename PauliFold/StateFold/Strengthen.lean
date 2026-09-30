import PauliFold.StateFold.Incomparable

/-!
# Strengthen: folding modulo witnessed constraints

Amy and Lunderville (POPL 2025, Algorithm 2) strengthen StateFold with constraints witnessed
by interference. Suppose a hidden variable `y` appears in no wire, and the phase is
`Φ + y · P` with `Φ` and `P` free of `y`. Then

  `∑_y (-1)^{y P} = 2 [P = 0]`,

so every path with `P = 1` cancels, even when `P = 0` cannot be solved as a substitution
`z := R` (which would be an HH step). Predicates may then be compared modulo the constraints:
two rotations fold if their predicates agree on the paths where every witnessed `P`
vanishes. With the field equations `x² = x`, ideal membership over `F₂` is exactly vanishing
on this set (Hilbert's Nullstellensatz). So, for the constraints collected here, this
models reduction modulo a Gröbner basis exactly. It is a *restricted* model of Algorithm 2,
which collects different constraints (see below), so it does not show that Feynman proves
no more. The degree bound `d` limits HH substitutions only; witnesses have no degree bound,
so `ProvesS 1` may use non-linear constraints.

**Modeling choice.** Witnesses are taken on the state at the end of the HH/ω trace, and
no constraint may mention a witness variable. Algorithm 2 collects constraints after every
rewrite. Carrying such constraints soundly through later substitutions `z := R` needs the
equations `z = R` as well, and a constraint that mentions `z` then breaks the pairing
argument. The paper's Example 23, where the `y₄` witness justifies folding `ℓ` with `ℓ''`,
is a witness on the final path sum.

Results:

- `denote_witnesses`: a state is insensitive to phase changes off the constrained set.
- `provesS_sound`: facts proved with Strengthen are true.
- `sfCandidateS_subset`: Strengthen refines plain StateFold.
- `gammaS_incomparable`: Strengthen and PauliFold are still incomparable, on the same two
  segments.
-/

namespace PauliFold.SF

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-- `f` does not depend on the variables in `W`. -/
def Avoids {α : Type*} (f : Val → α) (W : Set ℕ) : Prop :=
  ∀ ν v b, v ∈ W → f (Function.update ν v b) = f ν

/-- A witness on `A`: the hidden variable `y ∈ W`, in no wire, enters the phase as `y · P`,
where `P` avoids every witness variable `W`. -/
structure Witness (A : State n) (W : Set ℕ) (y : ℕ) (P : Poly) : Prop where
  hy : y ∈ A.temps
  hyW : y ∈ W
  ket : ∀ q, Avoids (fun ν => ev ν (A.ket q)) {y}
  suppP : Avoids (fun ν => ev ν P) W
  split : ∃ Φ : Val → ℚ, Avoids Φ {y} ∧ ∀ ν, ep (A.phase ν) = ep (Φ ν + b2q (ν y && ev ν P))

/-- **One witness.** Two phase shifts that avoid the witness variables and agree where
`P = 0` give the same denotation. -/
theorem denote_witness {A : State n} {W : Set ℕ} {y : ℕ} {P : Poly} (h : Witness A W y P)
    (hA : A.WF) {δ δ' : Val → ℚ} (hδ : Avoids δ W) (hδ' : Avoids δ' W)
    (hag : ∀ ν, ev ν P = false → ep (δ ν) = ep (δ' ν)) :
    (A.shift δ).denote = (A.shift δ').denote := by
  obtain ⟨hy, hyW, hket, hP, Φ, hΦ, hsplit⟩ := h
  funext o i
  set U := A.temps.erase y
  have hyU : y ∉ U := Finset.notMem_erase y _
  have hT : A.temps = insert y U := (Finset.insert_erase hy).symm
  have hyn := hA.temps_ge y hy
  rw [State.denote_eq, State.denote_eq]
  show invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset, (A.shift δ).term o (val i S) =
    invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset, (A.shift δ').term o (val i S)
  rw [hT, sum_split _ y hyU hyn, sum_split _ y hyU hyn]
  congr 1
  refine Finset.sum_congr rfl fun S hS => ?_
  have hyS : y ∉ S := fun h => hyU (Finset.mem_powerset.mp hS h)
  set ν := val i S
  have hνy : ν y = false := val_not_mem i S y hyn hyS
  have hk : ∀ q, ev (Function.update ν y true) (A.ket q) = ev ν (A.ket q) :=
    fun q => hket q ν y true rfl
  have hP1 : ev (Function.update ν y true) P = ev ν P := hP ν y true hyW
  have hΦ1 : Φ (Function.update ν y true) = Φ ν := hΦ ν y true rfl
  have hδ1 : δ (Function.update ν y true) = δ ν := hδ ν y true hyW
  have hδ1' : δ' (Function.update ν y true) = δ' ν := hδ' ν y true hyW
  unfold State.term State.shift
  simp only [hk]
  split_ifs
  · rw [ep_add, ep_add, ep_add, ep_add, hsplit, hsplit, hΦ1, hP1, hνy, Function.update_self,
      hδ1, hδ1', Bool.false_and, Bool.true_and, b2q_false, add_zero]
    have e1 : ep (Φ ν) * ep (δ ν) + ep (Φ ν) * ep (b2q (ev ν P)) * ep (δ ν) =
        ep (Φ ν) * ep (δ ν) * (1 + ep (b2q (ev ν P))) := by ring
    have e2 : ep (Φ ν) * ep (δ' ν) + ep (Φ ν) * ep (b2q (ev ν P)) * ep (δ' ν) =
        ep (Φ ν) * ep (δ' ν) * (1 + ep (b2q (ev ν P))) := by ring
    rw [ep_add (Φ ν) (b2q _), e1, e2, hh_sum]
    cases hPν : ev ν P
    · rw [hag ν hPν]
    · simp
  · simp

/-- **Several witnesses.** Two phase shifts that avoid the witness variables and agree where
every witnessed constraint vanishes give the same denotation. -/
theorem denote_witnesses {A : State n} {W : Set ℕ} (hA : A.WF) :
    ∀ (ws : List (ℕ × Poly)), (∀ p ∈ ws, Witness A W p.1 p.2) →
      ∀ {δ δ' : Val → ℚ}, Avoids δ W → Avoids δ' W →
      (∀ ν, (∀ p ∈ ws, ev ν p.2 = false) → ep (δ ν) = ep (δ' ν)) →
      (A.shift δ).denote = (A.shift δ').denote
  | [], _, δ, δ', _, _, hag =>
    denote_congr_ep _ _ rfl rfl fun ν => by
      show ep (A.phase ν + δ ν) = ep (A.phase ν + δ' ν)
      rw [ep_add, ep_add, hag ν fun p hp => absurd hp (List.not_mem_nil)]
  | (y, P) :: ws, hws, δ, δ', hδ, hδ', hag => by
    have hw := hws (y, P) (List.mem_cons_self ..)
    -- Switch to `δ'` where `P = 1`, which the witness on `y` cannot see.
    let δ₁ : Val → ℚ := fun ν => if ev ν P then δ' ν else δ ν
    have hδ₁ : Avoids δ₁ W := fun ν v b hv => by
      simp only [δ₁, hw.suppP ν v b hv, hδ ν v b hv, hδ' ν v b hv]
    have s1 : (A.shift δ).denote = (A.shift δ₁).denote :=
      denote_witness hw hA hδ hδ₁ fun ν hPν => by simp [δ₁, hPν]
    have s2 : (A.shift δ₁).denote = (A.shift δ').denote :=
      denote_witnesses hA ws (fun p hp => hws p (List.mem_cons_of_mem _ hp)) hδ₁ hδ'
        fun ν hν => by
          by_cases hPν : ev ν P
          · simp [δ₁, hPν]
          · simp only [δ₁, hPν]
            exact hag ν fun p hp => by
              rcases List.mem_cons.1 hp with rfl | hp
              · simpa using hPν
              · exact hν p hp
    exact s1.trans s2

/-- The end of a trace is well formed. -/
theorem Trace.wf {d : ℕ} {A C : State n} {f g f' g' : Poly} (htr : Trace d A f g C f' g')
    (hA : A.WF) : C.WF := by
  induction htr with
  | refl => exact hA
  | hh h _ _ _ ih => exact ih (wf_hh h hA)
  | omega h _ _ _ ih => exact ih (wf_omega h hA)

/-! ## The Strengthen concretization -/

/-- A degree-`d` trace, followed by witnesses on its final state, proves that the predicate
`f` equals (`s = false`) or complements (`s = true`) the wire predicate of `q` wherever every
witnessed constraint vanishes. -/
def ProvesS (d : ℕ) (A : State n) (f : Poly) (q : Fin n) (s : Bool) : Prop :=
  ∃ (C : State n) (f' g' : Poly) (W : Set ℕ) (ws : List (ℕ × Poly)),
    Trace d A f (A.ket q) C f' g' ∧ (∀ p ∈ ws, Witness C W p.1 p.2) ∧
    Avoids (fun ν => ev ν f') W ∧ Avoids (fun ν => ev ν g') W ∧
    ∀ ν, (∀ p ∈ ws, ev ν p.2 = false) → ev ν g' = (s != ev ν f')

/-- A plain StateFold fact is a Strengthen fact, with no witnesses. -/
theorem Proves.toS {d : ℕ} {A : State n} {f : Poly} {q : Fin n} {s : Bool}
    (h : Proves d A f q s) : ProvesS d A f q s := by
  obtain ⟨C, f', g', htr, heq⟩ := h
  exact ⟨C, f', g', ∅, [], htr, fun _ hp => absurd hp List.not_mem_nil,
    fun _ _ _ hv => absurd hv (Set.notMem_empty _), fun _ _ _ hv => absurd hv (Set.notMem_empty _),
    fun ν _ => heq ν⟩

/-- **Strengthen facts are true.** -/
theorem provesS_sound {d : ℕ} {mid : List (CTGate n)} {qk q : Fin n} {s : Bool}
    (hp : ProvesS d (run mid) (BoolPolynomial.var qk.val) q s) :
    axis mid qk = (zString {q} s).toMatrix := by
  obtain ⟨C, f', g', W, ws, htr, hws, hf', hg', heq⟩ := hp
  refine axis_of_trace htr ?_
  have hC : C.WF := htr.wf (wf_run mid)
  have hterm : ∀ a b : ℚ, Avoids (term2 a b f' g') W := fun a b ν v c hv => by
    simp only [term2, hf' ν v c hv, hg' ν v c hv]
  have hconst : Avoids (fun ν => term2 1 0 f' g' ν + if s then 1 else 0) W :=
    fun ν v c hv => by simp only [hterm 1 0 ν v c hv]
  rw [sgn_eq_ep, ← denote_shift_const]
  calc (C.shift (term2 0 1 f' g')).denote
      = (C.shift fun ν => term2 1 0 f' g' ν + if s then 1 else 0).denote :=
        denote_witnesses hC ws hws (hterm 0 1) hconst fun ν hν => by
          simp only [term2, heq ν hν]
          rw [← add_zero (0 * _ + _), ← add_zero (1 * _ + _ + _)]
          simpa using sign_key s 0 (ev ν f')
    _ = ((C.shift (term2 1 0 f' g')).shift fun _ => if s then 1 else 0).denote :=
        denote_congr' _ _ rfl rfl fun ν => by simp only [State.shift]; ring

/-- The Strengthen concretization: the operators consistent with every Strengthen fact. -/
def sfCandidateS (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) : Set (Op n) :=
  {O | ∀ q s, ProvesS d (run mid) (BoolPolynomial.var qk.val) q s →
    O = (zString {q} s).toMatrix}

theorem axis_mem_sfCandidateS (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) :
    axis mid qk ∈ sfCandidateS d mid qk :=
  fun _ _ hp => provesS_sound hp

/-- **Strengthen refines plain StateFold.** -/
theorem sfCandidateS_subset (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) :
    sfCandidateS d mid qk ⊆ sfCandidate d mid qk :=
  fun _ hO q s hp => hO q s hp.toS

theorem sfCandidateS_anti {d d' : ℕ} (hd : d ≤ d') (mid : List (CTGate n)) (qk : Fin n) :
    sfCandidateS d' mid qk ⊆ sfCandidateS d mid qk :=
  fun _ hO q s ⟨C, f', g', W, ws, htr, hws, hf', hg', heq⟩ =>
    hO q s ⟨C, f', g', W, ws, htr.mono hd, hws, hf', hg', heq⟩

/-- Both domains are sound, so they always overlap. -/
theorem axis_mem_bothS (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) :
    axis mid qk ∈ sfCandidateS d mid qk ∩ pauliCandidate mid qk :=
  ⟨axis_mem_sfCandidateS d mid qk, (candidate_sound qk mid).1⟩

/-! ## Incomparability with Strengthen -/

/-- On `h; t; tdg; h`, Strengthen is strictly more precise than PauliFold. -/
theorem midA_ssubsetS (d : ℕ) (hd : 1 ≤ d) :
    sfCandidateS d midA 0 ⊂ pauliCandidate midA 0 :=
  lt_of_le_of_lt (sfCandidateS_subset d midA 0) (midA_ssubset d hd)

/-- Every trace from a state whose wires mention all its path variables is empty. -/
theorem trace_refl_of_exposed {d : ℕ} {A C : State n} {f g f' g' : Poly}
    (hexp : ∀ y ∈ A.temps, ∃ q, ¬ Supp (fun ν => ev ν (A.ket q)) (A.dom \ {y}))
    (htr : Trace d A f g C f' g') : C = A ∧ f' = f ∧ g' = g := by
  cases htr with
  | refl => exact ⟨rfl, rfl, rfl⟩
  | hh h _ _ _ =>
    obtain ⟨q, hq⟩ := hexp _ h.hy
    exact absurd (h.ket q) hq
  | omega h _ _ _ =>
    obtain ⟨q, hq⟩ := hexp _ h.hy
    exact absurd (h.ket q) hq

/-- On `cx q0,q1; h q0; h q1; cx q1,q0; h q1`, no witness exists either, since every path
variable is in a wire: Strengthen proves nothing, at any degree. -/
theorem mid_not_provesS_var (d k : ℕ) (hk2 : k ≠ 2) (hk4 : k ≠ 4) (q : Fin 2) (s : Bool) :
    ¬ ProvesS d (run cexMid) (BoolPolynomial.var k) q s := by
  rintro ⟨C, f', g', W, ws, htr, hws, -, -, heq⟩
  have hexp : ∀ y ∈ (run cexMid).temps, ∃ q,
      ¬ Supp (fun ν => ev ν ((run cexMid).ket q)) ((run cexMid).dom \ {y}) := by
    intro y hy
    rcases (mid_temps y).1 hy with rfl | rfl | rfl
    · exact ⟨0, not_supp_of_flip _ (by simp [mid_ket0])⟩
    · exact ⟨0, not_supp_of_flip _ (by simp [mid_ket0])⟩
    · exact ⟨1, not_supp_of_flip _ (by simp [mid_ket1])⟩
  obtain ⟨rfl, rfl, rfl⟩ := trace_refl_of_exposed hexp htr
  -- No witness: each path variable flips a wire.
  have hnone : ∀ p ∈ ws, False := by
    intro p hp
    have hw := hws p hp
    have flip : ∀ q, ev (Function.update (fun _ => false) p.1 true) ((run cexMid).ket q) =
        ev (fun _ => false) ((run cexMid).ket q) := fun q => hw.ket q _ _ _ rfl
    rcases (mid_temps p.1).1 hw.hy with h | h | h
    · have := flip 0; rw [mid_ket0, mid_ket0, h] at this; simp at this
    · have := flip 0; rw [mid_ket0, mid_ket0, h] at this; simp at this
    · have := flip 1; rw [mid_ket1, mid_ket1, h] at this; simp at this
  exact mid_not_proves_var d k hk2 hk4 q s ⟨_, _, _, .refl _ _ _,
    fun ν => heq ν fun p hp => (hnone p hp).elim⟩

theorem mid_not_provesS (d : ℕ) (q : Fin 2) (s : Bool) :
    ¬ ProvesS d (run cexMid) (BoolPolynomial.var 1) q s :=
  mid_not_provesS_var d 1 (by norm_num) (by norm_num) q s

/-- On this segment, PauliFold is strictly more precise than Strengthen. -/
theorem mid_ssubsetS (d : ℕ) : pauliCandidate cexMid 1 ⊂ sfCandidateS d cexMid 1 := by
  have hS : sfCandidateS d cexMid 1 = Set.univ :=
    Set.eq_univ_of_forall fun _ q s hp => absurd hp (mid_not_provesS d q s)
  rw [hS]
  refine Set.ssubset_univ_iff.2 fun h => ?_
  have h0 : (0 : Op 2) ∈ pauliCandidate cexMid 1 := h ▸ Set.mem_univ _
  rw [pauliCandidate, cex_pauli] at h0
  exact zString_toMatrix_ne_zero _ _ h0.symm

/-- **Strengthen and PauliFold are incomparable**, at every degree `d ≥ 1`. -/
theorem gammaS_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      sfCandidateS d mid qk ⊂ pauliCandidate mid qk) ∧
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      pauliCandidate mid qk ⊂ sfCandidateS d mid qk) :=
  ⟨⟨1, midA, 0, midA_ssubsetS d hd⟩, ⟨2, cexMid, 1, mid_ssubsetS d⟩⟩

end

end PauliFold.SF
