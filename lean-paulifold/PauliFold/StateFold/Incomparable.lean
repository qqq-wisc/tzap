import PauliFold.StateFold.Counterexample
import TzapLean.SemanticsCheck

/-!
# StateFold and PauliFold are incomparable abstract domains

Both analyses track the axis of a candidate rotation on qubit `qₖ` through a Clifford+T
segment `mid`. The concrete value is the transported axis `U Z_{qₖ} U†`, where `U` is the
segment's unitary.

- **PauliFold** tracks one signed Pauli string, or `⊤`. Its concretization is
  `PauliAbs.gamma`.
- **StateFold** tracks the exact path sum of `mid`, together with the candidate's predicate
  `x_{qₖ}`. The path sum alone determines `U Z U†`, so an exact semantic concretization would
  make StateFold trivially exact. What the analysis *knows* are the facts its reductions
  prove: a degree-`d` trace makes `x_{qₖ}` equal or complementary to the wire predicate of
  qubit `q`. By the shifting lemma, such a fact means that the axis is `±Z_q`
  (`proves_sound`). The concretization `sfGamma d` is the set of operators consistent with
  every fact the analysis proves.

The results:

- **Both are sound** (`axis_mem_both`). The transported axis lies in both concretizations,
  so they always overlap.
- **Neither refines the other** (`gamma_incomparable`).
  - On `h; t; tdg; h`, StateFold proves the axis is `Z₀`, while PauliFold is `⊤`.
  - On `cx q0,q1; h q0; h q1; cx q1,q0; h q1`, PauliFold has the axis `Z₁`, while
    StateFold proves nothing at any degree.
-/

namespace PauliFold.SF

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Matrix facts -/

theorem toGate_wf : ∀ g : CTGate n, g.toGate.Wf
  | .cx _ _ hne => show _ ≠ _ from Fin.val_ne_of_ne hne
  | .cz _ _ hne => show _ ≠ _ from Fin.val_ne_of_ne hne
  | .h _ | .x _ | .z _ | .s _ | .sdg _ | .t _ | .tdg _ => trivial

/-- The unitary of a Clifford+T list is unitary. -/
theorem unitary_mul_conjTranspose (gs : List (CTGate n)) :
    unitary n (gs.map CTGate.toGate) * (unitary n (gs.map CTGate.toGate))ᴴ = 1 := by
  refine mul_eq_one_comm.mp ?_
  induction gs with
  | nil => simp
  | cons g gs ih =>
    simp only [List.map_cons, unitary_cons, Matrix.conjTranspose_mul]
    rw [Matrix.mul_assoc, ← Matrix.mul_assoc _ (unitary n _), ih, Matrix.one_mul,
      gateUnitary_unitary n _ (toGate_wf g)]

/-- The concrete transformer of a segment is conjugation by its unitary. -/
theorem concreteRun_eq (gs : List (CTGate n)) (A : Op n) :
    concreteRun gs A = conj (unitary n (gs.map CTGate.toGate)) A := by
  induction gs generalizing A with
  | nil => simp [concreteRun, conj]
  | cons g gs ih =>
    show concreteRun gs (g.concrete A) = _
    rw [ih]
    simp [CTGate.concrete, conj, Matrix.mul_assoc, Matrix.conjTranspose_mul]

/-- `Z_q` as a gate matrix. -/
abbrev zMat (q : Fin n) : Op n := embed1 n (diag2 1 (ep 1)) q.val

/-- The signed Pauli string `(-1)^s Z_q` is the matrix `(-1)^s · Z_q`. -/
theorem zString_toMatrix (q : Fin n) (s : Bool) :
    (zString {q} s).toMatrix = SPauli.sgn s • zMat q := by
  funext o i
  show _ = SPauli.sgn s * embed1 n (diag2 1 (ep 1)) q.val o i
  rw [embed1_apply_of_lt _ q.isLt]
  unfold SPauli.toMatrix
  rw [prod_split _ q]
  simp only [zString, Finset.mem_singleton, if_true, Fin.eta]
  have hrest : ∏ r ∈ Finset.univ.erase q, (if r = q then Letter.Z else Letter.I).mat (o r) (i r) =
      if ∀ r : Fin n, (r : ℕ) ≠ q.val → o r = i r then 1 else 0 := by
    split_ifs with h
    · refine Finset.prod_eq_one fun r hr => ?_
      have hr' : r ≠ q := Finset.ne_of_mem_erase hr
      simp [hr', Letter.mat, h r (Fin.val_ne_of_ne hr')]
    · push Not at h
      obtain ⟨r, hr, hne⟩ := h
      have hr' : r ≠ q := fun h => hr (congrArg Fin.val h)
      exact Finset.prod_eq_zero (Finset.mem_erase.2 ⟨hr', Finset.mem_univ r⟩)
        (by simp [hr', Letter.mat, hne])
  rw [hrest]
  split_ifs <;> cases o q <;> cases i q <;> simp [Letter.mat, diag2]

theorem zString_toMatrix_ne_zero (Q : Finset (Fin n)) (s : Bool) : (zString Q s).toMatrix ≠ 0 := by
  intro h
  have h' := congrFun (congrFun h fun _ => false) fun _ => false
  unfold SPauli.toMatrix at h'
  rw [Finset.prod_eq_one fun r _ => by
    simp only [zString]; split_ifs <;> rfl] at h'
  cases s <;> simp [SPauli.sgn, zString] at h'

/-- The denotation depends on the phase only through `e^{iπΦ}`. -/
theorem denote_congr_ep (A B : State n) (hk : A.ket = B.ket) (ht : A.temps = B.temps)
    (hp : ∀ ν, ep (A.phase ν) = ep (B.phase ν)) : A.denote = B.denote := by
  funext o i
  simp only [State.denote_eq, State.term, hk, ht, hp]

theorem ep_add_two' (θ : ℚ) : ep (θ + 2) = ep θ := by
  rw [show θ + 2 = θ + 1 + 1 by ring, ep_add, ep_add, ep_one]; ring

/-! ## The StateFold concretization -/

/-- A degree-`d` trace proves that the predicate `f` equals (`s = false`) or complements
(`s = true`) the wire predicate of qubit `q`. -/
def Proves (d : ℕ) (A : State n) (f : Poly) (q : Fin n) (s : Bool) : Prop :=
  ∃ (C : State n) (f' g' : Poly), Trace d A f (A.ket q) C f' g' ∧ ∀ ν, ev ν g' = (s != ev ν f')

/-- The StateFold concretization of the state `A` with candidate predicate `f`: the operators
consistent with every fact that the degree-`d` traces prove. -/
def sfGamma (d : ℕ) (A : State n) (f : Poly) : Set (Op n) :=
  {O | ∀ q s, Proves d A f q s → O = (zString {q} s).toMatrix}

/-- StateFold's abstract value for a candidate on `qₖ` after the segment `mid`. -/
def sfCandidate (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) : Set (Op n) :=
  sfGamma d (run mid) (BoolPolynomial.var qk.val)

/-- PauliFold's abstract value for a candidate on `qₖ` after the segment `mid`. -/
def pauliCandidate (mid : List (CTGate n)) (qk : Fin n) : Set (Op n) :=
  (PauliAbs.run mid (.one (zString {qk} false))).gamma

/-- The concrete value: the candidate's axis transported through `mid`. -/
def axis (mid : List (CTGate n)) (qk : Fin n) : Op n :=
  concreteRun mid (zString {qk} false).toMatrix

/-- A trace whose end makes the two shifts differ by `(-1)^s` pins the transported axis to
`(-1)^s Z_q`. -/
theorem axis_of_trace {d : ℕ} {mid : List (CTGate n)} {qk q : Fin n} {s : Bool}
    {C : State n} {f' g' : Poly}
    (htr : Trace d (run mid) (BoolPolynomial.var qk.val) ((run mid).ket q) C f' g')
    (hend : (C.shift (term2 0 1 f' g')).denote =
      SPauli.sgn s • (C.shift (term2 1 0 f' g')).denote) :
    axis mid qk = (zString {q} s).toMatrix := by
  set A := run mid
  set U := unitary n (mid.map CTGate.toGate)
  set f : Poly := BoolPolynomial.var qk.val
  have hA : A.WF := wf_run mid
  have hAU : A.denote = U := denote_run mid
  -- Adding `[x_{qₖ}]` to the phase is a Z on `qₖ` before the segment.
  have h1 : (A.shift (term2 1 0 f (A.ket q))).denote = U * zMat qk := by
    have r := ((rel_rot 0 1 qk (init n)).runFrom 0 1 mid).denote
    have hz : runFrom 1 mid (rot 0 1 qk (init n)) = run (.z qk :: mid) := rfl
    rw [hz, denote_run] at r
    rw [show U * zMat qk = unitary n ((CTGate.z qk :: mid).map CTGate.toGate) from rfl, r]
    exact denote_congr' _ _ rfl rfl fun ν => by
      simp only [State.shift, term2, zero_mul, add_zero]; rfl
  -- Adding `[σ_q]` is a Z on `q` after it.
  have h2 : (A.shift (term2 0 1 f (A.ket q))).denote = zMat q * U := by
    rw [← hAU, ← denote_rot 0 1 q A]
    exact denote_congr' _ _ rfl rfl fun ν => by simp [State.shift, rot, term2]
  -- Both shifts ride along the trace; at its end they differ by `(-1)^s`.
  have hf : Supp (fun ν => ev ν f) A.dom := fun ν ν' h => by
    simp only [f, ev_var]; exact h _ (Or.inl qk.isLt)
  have t1 := trace_shift htr hA hf (hA.ket q) 1 0
  have t2 := trace_shift htr hA hf (hA.ket q) 0 1
  have hmain : zMat q * U = SPauli.sgn s • (U * zMat qk) := by
    rw [← h2, ← h1, t1, t2, hend]
  -- Conclude `U Z_{qₖ} U† = (-1)^s Z_q`.
  have hsq : SPauli.sgn s * SPauli.sgn s = 1 := by cases s <;> simp [SPauli.sgn]
  have hUZ : U * zMat qk = SPauli.sgn s • (zMat q * U) := by
    rw [hmain, smul_smul, hsq, one_smul]
  unfold axis
  rw [concreteRun_eq, zString_toMatrix, zString_toMatrix, SPauli.sgn, if_neg Bool.false_ne_true,
    one_smul]
  unfold conj
  rw [hUZ, Matrix.smul_mul, Matrix.mul_assoc, unitary_mul_conjTranspose, Matrix.mul_one]

/-- `[s ⊕ x] ≡ [x] + s` modulo 2. -/
theorem sign_key (s : Bool) (P : ℚ) (x : Bool) : ep (P + (0 * b2q x + 1 * b2q (s != x))) =
    ep (P + (1 * b2q x + 0 * b2q (s != x)) + if s then 1 else 0) := by
  cases s <;> cases x <;> simp [b2q]
  rw [show P + 1 + 1 = P + 2 by ring, ep_add_two']

theorem sgn_eq_ep (s : Bool) : SPauli.sgn s = ep (if s then 1 else 0) := by
  cases s <;> simp [SPauli.sgn]

/-- **A proved fact is true.** If a trace makes `x_{qₖ}` equal or complementary to the wire
predicate of `q` after `mid`, the transported axis is `(-1)^s Z_q`. -/
theorem proves_sound {d : ℕ} {mid : List (CTGate n)} {qk q : Fin n} {s : Bool}
    (hp : Proves d (run mid) (BoolPolynomial.var qk.val) q s) :
    axis mid qk = (zString {q} s).toMatrix := by
  obtain ⟨C, f', g', htr, heq⟩ := hp
  refine axis_of_trace htr ?_
  rw [sgn_eq_ep, ← denote_shift_const]
  exact denote_congr_ep _ _ rfl rfl fun ν => by
    simp only [State.shift, term2, heq]; exact sign_key _ _ _

/-- The StateFold concretization is sound. -/
theorem axis_mem_sfCandidate (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) :
    axis mid qk ∈ sfCandidate d mid qk :=
  fun _ _ hp => proves_sound hp

/-- **Both domains are sound, so they always overlap:** the transported axis lies in both. -/
theorem axis_mem_both (d : ℕ) (mid : List (CTGate n)) (qk : Fin n) :
    axis mid qk ∈ sfCandidate d mid qk ∩ pauliCandidate mid qk :=
  ⟨axis_mem_sfCandidate d mid qk, (candidate_sound qk mid).1⟩

/-- A higher degree proves more facts, so its concretization is smaller. -/
theorem sfCandidate_anti {d d' : ℕ} (hd : d ≤ d') (mid : List (CTGate n)) (qk : Fin n) :
    sfCandidate d' mid qk ⊆ sfCandidate d mid qk := fun _ hO q s ⟨C, f', g', htr, heq⟩ =>
  hO q s ⟨C, f', g', htr.mono hd, heq⟩

/-! ## `h; t; tdg; h`: StateFold is more precise -/

/-- The segment `h q0; t q0; tdg q0; h q0`. -/
def midA : List (CTGate 1) := [.h 0, .t 0, .tdg 0, .h 0]

theorem midA_temps (y : ℕ) : y ∈ (run midA).temps ↔ y = 1 ∨ y = 2 := by
  simp [run, midA, runFrom, step, rot, init]
  omega

theorem midA_ket (ν : Val) : ev ν ((run midA).ket 0) = ν 2 := by
  simp [run, midA, runFrom, step, rot, init]

theorem midA_phase (ν : Val) : (run midA).phase ν = b2q (ν 0 && ν 1) + b2q (ν 1 && ν 2) := by
  simp [run, midA, runFrom, step, rot, init]
  ring

/-- One HH step (`y₁` summed out, `y₂ := x₀`) proves that the axis is `Z₀`. -/
theorem midA_proves (d : ℕ) (hd : 1 ≤ d) :
    Proves d (run midA) (BoolPolynomial.var 0) 0 false := by
  set A := run midA
  have h0 : (0 : ℕ) ∈ A.dom \ {1} := ⟨Or.inl (by norm_num), by simp⟩
  have h2 : (2 : ℕ) ∈ A.dom \ {1} := ⟨Or.inr ((midA_temps 2).2 (by simp)), by simp⟩
  have hh : HH d A (hhResult A 1 2 (BoolPolynomial.var 0) fun _ => 0) 1 2
      (BoolPolynomial.var 0) fun _ => 0 :=
    { hy := (midA_temps 1).2 (by simp)
      hz := (midA_temps 2).2 (by simp)
      hyz := by decide
      ket := fun q ν ν' h => by
        obtain rfl : q = 0 := Subsingleton.elim _ _
        show ev ν ((run midA).ket 0) = ev ν' ((run midA).ket 0)
        rw [midA_ket, midA_ket]
        exact h 2 h2
      deg := by
        show (MvPolynomial.X 0 : Poly).totalDegree ≤ d
        rw [MvPolynomial.totalDegree_X]; exact hd
      suppR := fun ν ν' h => by
        simp only [ev_var]
        exact h 0 ⟨Or.inl (by norm_num), by simp⟩
      supp₀ := fun _ _ _ => rfl
      split := fun ν => by
        rw [midA_phase]
        simp only [ev_var]
        cases ν 0 <;> cases ν 1 <;> cases ν 2 <;> simp [b2q] <;> norm_num
      result := rfl }
  refine ⟨_, _, _, .hh hh ?_ ?_ (.refl _ _ _), fun ν => ?_⟩
  · intro ν ν' h
    simp only [ev_var]; exact h 0 h0
  · intro ν ν' h
    show ev ν ((run midA).ket 0) = ev ν' ((run midA).ket 0)
    rw [midA_ket, midA_ket]; exact h 2 h2
  · rw [ev_subst, ev_subst, midA_ket, ev_var]
    simp

/-- The Pauli domain reaches `⊤` at the inner T: after the first H the axis is `X₀`. -/
theorem midA_pauli : PauliAbs.run midA (.one (zString {0} false)) = .top := by
  simp [PauliAbs.run, midA, PauliAbs.step, SPauli.conjClifford, SPauli.apply1, zString,
    Letter.conjH, Letter.anticommutesZ]

/-- On `h; t; tdg; h`, StateFold's concretization is strictly smaller than PauliFold's, at
every degree `d ≥ 1`. -/
theorem midA_ssubset (d : ℕ) (hd : 1 ≤ d) :
    sfCandidate d midA 0 ⊂ pauliCandidate midA 0 := by
  have hP : pauliCandidate midA 0 = Set.univ := by
    rw [pauliCandidate, midA_pauli]; rfl
  rw [hP]
  refine Set.ssubset_univ_iff.2 fun h => ?_
  have h0 : (0 : Op 1) ∈ sfCandidate d midA 0 := h ▸ Set.mem_univ _
  exact zString_toMatrix_ne_zero _ _ (h0 0 false (midA_proves d hd)).symm

/-! ## `cx q0,q1; h q0; h q1; cx q1,q0; h q1`: PauliFold is more precise -/

theorem mid_temps (y : ℕ) : y ∈ (run cexMid).temps ↔ y = 2 ∨ y = 3 ∨ y = 4 := by
  simp [run, cexMid, runFrom, step, init]
  omega

theorem mid_ket0 (ν : Val) : ev ν ((run cexMid).ket 0) = (ν 2 != ν 3) := by
  simp [run, cexMid, runFrom, step, init]

theorem mid_ket1 (ν : Val) : ev ν ((run cexMid).ket 1) = ν 4 := by
  simp [run, cexMid, runFrom, step, init]

/-- Every path variable appears in a wire, so no trace applies, and neither wire predicate is
equal or complementary to an input variable `x_k`: StateFold proves nothing, at any degree. -/
theorem mid_not_proves_var (d k : ℕ) (hk2 : k ≠ 2) (hk4 : k ≠ 4) (q : Fin 2) (s : Bool) :
    ¬ Proves d (run cexMid) (BoolPolynomial.var k) q s := by
  rintro ⟨C, f', g', htr, heq⟩
  have hexp : ∀ y ∈ (run cexMid).temps, ∃ q,
      ¬ Supp (fun ν => ev ν ((run cexMid).ket q)) ((run cexMid).dom \ {y}) := by
    intro y hy
    rcases (mid_temps y).1 hy with rfl | rfl | rfl
    · exact ⟨0, not_supp_of_flip _ (by simp [mid_ket0])⟩
    · exact ⟨0, not_supp_of_flip _ (by simp [mid_ket0])⟩
    · exact ⟨1, not_supp_of_flip _ (by simp [mid_ket1])⟩
  obtain ⟨rfl, rfl⟩ := trace_eq_of_exposed hexp htr
  have h₁ := heq fun _ => false
  fin_cases q
  · have h₂ := heq (Function.update (fun _ => false) 2 true)
    simp only [Fin.zero_eta, mid_ket0, ev_var, Function.update_of_ne hk2] at h₁ h₂
    cases s <;> simp at h₁ h₂
  · have h₂ := heq (Function.update (fun _ => false) 4 true)
    simp only [Fin.mk_one, mid_ket1, ev_var, Function.update_of_ne hk4] at h₁ h₂
    cases s <;> simp at h₁ h₂

theorem mid_not_proves (d : ℕ) (q : Fin 2) (s : Bool) :
    ¬ Proves d (run cexMid) (BoolPolynomial.var 1) q s :=
  mid_not_proves_var d 1 (by norm_num) (by norm_num) q s

/-- On this segment, PauliFold's concretization is strictly smaller than StateFold's, at every
degree. -/
theorem mid_ssubset (d : ℕ) : pauliCandidate cexMid 1 ⊂ sfCandidate d cexMid 1 := by
  have hS : sfCandidate d cexMid 1 = Set.univ :=
    Set.eq_univ_of_forall fun _ q s hp => absurd hp (mid_not_proves d q s)
  rw [hS]
  refine Set.ssubset_univ_iff.2 fun h => ?_
  have h0 : (0 : Op 2) ∈ pauliCandidate cexMid 1 := h ▸ Set.mem_univ _
  rw [pauliCandidate, cex_pauli] at h0
  exact zString_toMatrix_ne_zero _ _ h0.symm

/-! ## Incomparability -/

/-- **StateFold and PauliFold are incomparable.** Neither concretization is always contained
in the other: at every degree `d ≥ 1` there are segments on which each domain is strictly
more precise than the other. -/
theorem gamma_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      sfCandidate d mid qk ⊂ pauliCandidate mid qk) ∧
    (∃ (n : ℕ) (mid : List (CTGate n)) (qk : Fin n),
      pauliCandidate mid qk ⊂ sfCandidate d mid qk) :=
  ⟨⟨1, midA, 0, midA_ssubset d hd⟩, ⟨2, cexMid, 1, mid_ssubset d⟩⟩

end

end PauliFold.SF
