import TzapLean.SuperOpt

/-!
# `SuperOpt`: correctness

The scan proposes a whole *set* of rewrites and then every part of it is verified, so nothing
about the search enters the proof — not the order it considers windows in, not the closure it
computes, not the table it looks up. Three checks stand between a proposal and a rewritten
circuit:

* `checkRewrite` — this file. A candidate that passes `accepts` denotes the same operator as
  the gates it replaces, up to global phase, hence is equivalent to them. This is where
  `Locality.equivalent_of_local_smul` earns its keep: the matrices compared are the window's
  own, on its own wires. It runs once per *selected* rewrite, not once per window examined.
* `sepB` and `onSuppB` — `TzapLean/Rewrite.lean`. The two conditions under which a set of
  scattered, interleaved rewrites may all be spliced in at once (`applyAll_correct`).

A rewrite that fails the first is untagged and the rest still stand; if the last two fail,
nothing is rewritten. In practice none of them fire.
-/

namespace TzapLean

open ExactMat

/-! ## Localizing and globalizing preserve well-formedness -/

theorem Wf_mapQubits {f : Qubit → Qubit} {g : Gate} (hwf : g.Wf)
    (hinj : ∀ q ∈ g.qubitsOf, ∀ r ∈ g.qubitsOf, f q = f r → q = r) : (mapQubits f g).Wf := by
  cases g <;>
    simp_all [mapQubits, Gate.Wf, Gate.qubitsOf]

theorem localIdxD_inj {S : List Qubit} {q r : Qubit} (hq : q ∈ S) (hr : r ∈ S)
    (h : localIdxD S q = localIdxD S r) : q = r := by
  obtain ⟨i, hi, hdi⟩ := exists_localIdx hq
  obtain ⟨j, hj, hdj⟩ := exists_localIdx hr
  rw [hdi, hdj] at h
  subst h
  exact (localIdx_getD hi).symm.trans (localIdx_getD hj)

theorem getD_inj {S : List Qubit} (hnd : S.Nodup) {i j : Nat} (hi : i < S.length)
    (hj : j < S.length) (h : S.getD i 0 = S.getD j 0) : i = j := by
  have h1 := localIdx_getD_self hnd hi
  have h2 := localIdx_getD_self hnd hj
  rw [h, h2] at h1
  exact (Option.some.injEq _ _ ▸ h1).symm

/-! ## What `accepts` establishes -/

theorem applyGate_isUnitary {n : Nat} {g : Gate} {M M' : ExactMat n}
    (h : applyGate g M = some M') : g.isUnitary = true := by
  cases g <;> simp_all [applyGate, Gate.isUnitary, Gate.isMeasurement]

theorem matrixOfFrom_isUnitary {n : Nat} : ∀ (gs : List Gate) {M M' : ExactMat n},
    matrixOfFrom M gs = some M' → ∀ g ∈ gs, g.isUnitary = true := by
  intro gs
  induction gs with
  | nil => intro M M' _ g hg; simp at hg
  | cons a as ih =>
      intro M M' h g hg
      rw [matrixOfFrom] at h
      cases ha : applyGate a M with
      | none => rw [ha] at h; exact absurd h (by simp)
      | some M₁ =>
          rw [ha, Option.bind_some] at h
          rcases List.mem_cons.1 hg with rfl | hg
          · exact applyGate_isUnitary ha
          · exact ih h g hg

theorem matrixOf_isUnitary {n : Nat} {gs : List Gate} {M : ExactMat n}
    (h : matrixOf n gs = some M) : ∀ g ∈ gs, g.isUnitary = true :=
  matrixOfFrom_isUnitary gs h

theorem accepts_spec {k : Nat} {target : ExactMat k} {cand : List Gate}
    (h : accepts target cand = true) :
    (∀ g ∈ cand, ∀ q ∈ g.qubitsOf, q < k) ∧ (∀ g ∈ cand, g.Wf) ∧
      ∃ N p, matrixOf k cand = some N ∧ phaseMatch target N.normalize = some p := by
  rw [accepts, Bool.and_eq_true] at h
  obtain ⟨hall, hmat⟩ := h
  rw [List.all_eq_true] at hall
  refine ⟨fun g hg q hq => ?_, fun g hg => ?_, ?_⟩
  · have := hall g hg
    rw [Bool.and_eq_true, List.all_eq_true] at this
    exact of_decide_eq_true (this.1 q hq)
  · have := hall g hg
    rw [Bool.and_eq_true] at this
    exact of_decide_eq_true this.2
  · cases hN : matrixOf k cand with
    | none => rw [hN] at hmat; exact absurd hmat (by simp)
    | some N =>
        rw [hN] at hmat
        simp only at hmat
        cases hp : phaseMatch target N.normalize with
        | none => rw [hp] at hmat; exact absurd hmat (by simp)
        | some p => exact ⟨N, p, rfl, hp⟩

/-! ## A verified replacement is equivalent to the gates it replaces

The search that proposed the replacement is unverified — a table lookup on a flat matrix,
with no `ExactMat` anywhere. `checkRewrite` is what makes that safe, and this is what
`checkRewrite` buys: the two gate lists denote the same operator, so splicing one in for the
other is meaning-preserving.

It runs once per *selected* rewrite. The arrangement it replaces built an exact matrix for
every window that got past a filter, which is why the filter had to exist. -/

theorem checkRewrite_correct {n m : Nat} {S : List Qubit} {members cand : List Gate}
    (h : checkRewrite n S members cand = true) :
    Equivalent n m (cand.map (globalizeGate S)) members ∧
      (∀ g ∈ cand.map (globalizeGate S), g.Wf) ∧
      (∀ g ∈ cand.map (globalizeGate S), ∀ q ∈ g.qubitsOf, q ∈ S) ∧
      (∀ g ∈ cand.map (globalizeGate S), g.isUnitary = true) ∧
      (∀ q ∈ S, q < n) := by
  rw [checkRewrite] at h
  simp only [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at h
  obtain ⟨⟨⟨⟨hnd, hrange⟩, hmem⟩, -⟩, hmat⟩ := h
  have hrange' : ∀ q ∈ S, q < n := fun q hq => hrange q hq
  -- the members live on `S` and are well-formed
  have hsub : ∀ g ∈ members, ∀ q ∈ g.qubitsOf, q ∈ S := by
    intro g hg q hq
    have hc := (hmem g hg).1 q hq
    simpa [List.contains_eq_mem] using hc
  have hwf : ∀ g ∈ members, g.Wf := fun g hg => (hmem g hg).2
  set k := S.length with hk
  cases hM : ExactMat.matrixOf k (localizeGates S members) with
  | none => rw [hM] at hmat; exact absurd hmat (by simp)
  | some M =>
      rw [hM] at hmat
      obtain ⟨hqb, hwfc, N, p, hN, hp⟩ := accepts_spec hmat
      -- the localized members are well-formed
      have hwfl : ∀ g ∈ localizeGates S members, g.Wf := by
        intro g hg
        rcases List.mem_map.1 hg with ⟨g', hg', rfl⟩
        refine Wf_mapQubits (hwf g' hg') ?_
        intro q hq r hr hqr
        exact localIdxD_inj (hsub g' hg' q hq) (hsub g' hg' r hr) hqr
      -- globalizing the candidate and localizing it again is the identity on it
      have hround : ∀ g ∈ cand, localizeGate S (globalizeGate S g) = g := by
        intro g hg
        refine mapQubits_comp ?_
        intro q hq
        exact localIdxD_eq (localIdx_getD_self hnd (hqb g hg q hq))
      have hlocal : localizeGates S (cand.map (globalizeGate S)) = cand := by
        rw [localizeGates, List.map_map]
        refine Eq.trans (List.map_congr_left ?_) (List.map_id cand)
        intro g hg
        exact hround g hg
      have hreplsub : ∀ g ∈ cand.map (globalizeGate S), ∀ q ∈ g.qubitsOf, q ∈ S := by
        intro g hg q hq
        rcases List.mem_map.1 hg with ⟨g', hg', rfl⟩
        rw [globalizeGate, qubitsOf_mapQubits] at hq
        rcases List.mem_map.1 hq with ⟨i, hi, rfl⟩
        exact getD_mem (hqb g' hg' i hi)
      have hreplwf : ∀ g ∈ cand.map (globalizeGate S), g.Wf := by
        intro g hg
        rcases List.mem_map.1 hg with ⟨g', hg', rfl⟩
        refine Wf_mapQubits (hwfc g' hg') ?_
        intro q hq r hr hqr
        exact getD_inj hnd (hqb g' hg' q hq) (hqb g' hg' r hr) hqr
      have hcandu : ∀ g ∈ cand, g.isUnitary = true := matrixOf_isUnitary hN
      have hreplu : ∀ g ∈ cand.map (globalizeGate S), g.isUnitary = true := by
        intro g hg
        rcases List.mem_map.1 hg with ⟨g', hg', rfl⟩
        rw [globalizeGate, isUnitary_mapQubits]
        exact hcandu g' hg'
      have hmemu : ∀ g ∈ members, g.isUnitary = true := by
        intro g hg
        have := matrixOf_isUnitary hM (localizeGate S g) (List.mem_map.2 ⟨g, hg, rfl⟩)
        rwa [localizeGate, isUnitary_mapQubits] at this
      have h1 : N.interp = unitary k cand := matrixOf_sound hwfc hN
      have h2 : M.interp = unitary k (localizeGates S members) := matrixOf_sound hwfl hM
      have h3 : (N.normalize).interp = ω ^ p • (M.normalize).interp := phaseMatch_sound hp
      rw [interp_normalize, interp_normalize, h1, h2] at h3
      refine ⟨?_, hreplwf, hreplsub, hreplu, hrange'⟩
      exact equivalent_of_local_smul hnd hrange' hsub hreplsub hmemu hreplu (ω ^ p)
        (omega_pow_unit p) (by rw [hlocal]; exact h3)

/-! ## The pass

Everything above is about one rewrite. `applyAll_correct` in `TzapLean/Rewrite.lean` is about
a whole set of them, and the two conditions it needs — `Sep` and `OnSupp` — are decided by
`sepB` and `onSuppB`. So the pass is: propose (unverified), vet each rewrite, check the two
conditions, splice. If a check fails the proposal is dropped and the circuit comes back
unchanged; nothing about how the scan works enters the proof. -/

/-- The tagging carries the gates it was built from. -/
@[simp] theorem untag_tagged (st : Scan) (gs : List Gate) : untag (st.tagged gs) = gs := by
  simp only [Scan.tagged, untag, List.map_map]
  refine Eq.trans (List.map_congr_left ?_) (List.zipIdx_map_fst 0 gs)
  intro p _
  rfl

/-- Every rewrite that survives vetting replaces the gates it claims by an equivalent list. -/
theorem vetted_repl {n m : Nat} (st : Scan) (xs : List Tagged) (w : Nat)
    (hw : ∃ p ∈ st.vetted n xs, p.2 = some w) :
    Equivalent n m (st.repl w) (claimedBy w (st.vetted n xs)) ∧
      (∀ g ∈ st.repl w, g.Wf) ∧ (∀ g ∈ st.repl w, ∀ q ∈ g.qubitsOf, q ∈ st.supp w) ∧
      (∀ g ∈ st.repl w, g.isUnitary = true) ∧ (∀ q ∈ st.supp w, q < n) := by
  let groups := groupClaimsAux (claimBound xs) xs
  let members : Nat → List Gate := fun w => groups[w]?.getD []
  have hkeep : st.keep n members w = true := keep_of_mem_vettedBy hw
  have hcl : claimedBy w (st.vetted n xs) = claimedBy w xs := by
    rw [Scan.vetted, claimedBy_vettedBy, if_pos hkeep]
  rw [hcl]
  have hm : members w = claimedBy w xs := by
    simpa [members, groups, groupedClaim] using groupedClaim_eq_claimedBy xs w
  rw [← hm]
  exact checkRewrite_correct (m := m) hkeep

/-- **Superoptimization preserves meaning.** -/
theorem superOptGates_correct {n m : Nat} (cfg : SuperOptConfig) (tbl : SynthTable)
    (gs : List Gate) : Equivalent n m (superOptGates cfg tbl n gs) gs := by
  rw [superOptGates]
  split
  · rename_i hchecks
    simp only [Bool.and_eq_true] at hchecks
    set st := proposeRewrites cfg tbl n gs.toArray with hst
    set xs := st.vetted n (st.tagged gs) with hxs
    have hgates : untag xs = gs := by rw [hxs, Scan.vetted, untag_vettedBy, untag_tagged]
    have hmain := applyAll_correct (n := n) (m := m) st.supp st.repl xs
      (sepAllB_sound st.supp xs hchecks.2)
      (onSuppB_sound hchecks.1) (fun w hw => (vetted_repl st (st.tagged gs) w hw).1)
    rw [applyAllLinear_initial]
    rwa [hgates] at hmain
  · exact Equivalent.refl n m gs

/-- The two halves of the structural invariant, from the same check. -/
theorem superOptGates_pred {n m : Nat} {P : Gate → Prop} (cfg : SuperOptConfig)
    (tbl : SynthTable) (gs : List Gate) (hin : ∀ g ∈ gs, P g)
    (hrepl : ∀ (S : List Qubit) (r : List Gate), (∀ g ∈ r, g.Wf) →
      (∀ g ∈ r, ∀ q ∈ g.qubitsOf, q ∈ S) → (∀ g ∈ r, g.isUnitary = true) →
      (∀ q ∈ S, q < n) → ∀ g ∈ r, P g) :
    ∀ g ∈ superOptGates cfg tbl n gs, P g := by
  rw [superOptGates]
  split
  · set st := proposeRewrites cfg tbl n gs.toArray with hst
    set xs := st.vetted n (st.tagged gs) with hxs
    have hgates : untag xs = gs := by rw [hxs, Scan.vetted, untag_vettedBy, untag_tagged]
    rw [applyAllLinear_initial]
    refine applyAll_pred st.repl xs (hgates ▸ hin) (fun w hw => ?_)
    obtain ⟨-, hwf, hsub, hu, hrange⟩ := vetted_repl (m := m) st (st.tagged gs) w hw
    exact hrepl (st.supp w) (st.repl w) hwf hsub hu hrange
  · exact hin

/-- **Superoptimization preserves well-formedness.** -/
theorem superOptGates_wf {n : Nat} (cfg : SuperOptConfig) (tbl : SynthTable) (gs : List Gate)
    (hwf : ∀ g ∈ gs, g.Wf) : ∀ g ∈ superOptGates cfg tbl n gs, g.Wf :=
  superOptGates_pred (m := 0) cfg tbl gs hwf (fun _ _ h _ _ _ => h)

/-- **Superoptimization keeps every operand in range.** -/
theorem superOptGates_inRange {n m : Nat} (cfg : SuperOptConfig) (tbl : SynthTable)
    (gs : List Gate) (hin : ∀ g ∈ gs, g.InRange n m) :
    ∀ g ∈ superOptGates cfg tbl n gs, g.InRange n m :=
  superOptGates_pred (m := m) cfg tbl gs hin (fun S r _ hsub hu hrange g hg =>
    ⟨fun q hq => hrange q (hsub g hg q hq), fun b hb => by
      rw [Gate.cbitsOf_eq_nil_of_isUnitary (hu g hg)] at hb
      exact absurd hb (by simp)⟩)

/-- `SuperOpt`, as a `Pass`: every rewrite it takes is decided by exact matrix comparison, and
the set of them by `Sep` and `OnSupp`, so the obligations are discharged by checks the pass
already runs. -/
def SuperOpt (cfg : SuperOptConfig) (tbl : SynthTable) : Pass where
  name := "Superoptimization"
  run := fun c => ⟨superOpt cfg tbl c.raw, c.numQubits_eq, c.numCbits_eq,
    superOptGates_wf cfg tbl c.raw.gates c.wf⟩
  correct := by
    intro n m c
    rcases c with ⟨c, rfl, rfl, hc⟩
    exact superOptGates_correct cfg tbl c.gates

@[simp] theorem SuperOpt_run (cfg : SuperOptConfig) (tbl : SynthTable)
    (c : Circuit n m) : ((SuperOpt cfg tbl).run c).raw = superOpt cfg tbl c.raw := rfl

end TzapLean
