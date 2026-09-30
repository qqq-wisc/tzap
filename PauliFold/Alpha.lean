import PauliFold.Affine

/-!
# Amy and Lunderville's path-sum abstraction `α`

Definition 17 of Amy and Lunderville (POPL 2025) abstracts a path sum
`|x⟩ ↦ ∑_y Φ(x, y) |f(x, y)⟩` by the transition ideal `α = ∃Y. ⟨X' ⊕ f(X, Y)⟩`. Over `F₂`
with the field equations (their Proposition 12), the variety of this elimination ideal is
exactly the image of the wire map, `{(x, f(x, y)) | y}`. We therefore define its
concretization directly (`alphaGamma`): the matrices whose nonzero transitions `⟨o|U|i⟩` are
hit by some path of the state.

Their Algorithm 1 computes a path sum and rewrites it with the rules of their Figure 12, here
HH (substitution degree `≤ d`) and ω. `Reduces d` is that rewriting relation, and
`alphaGammaC d gs` keeps the unitaries in `alphaGamma C` for *every* reachable `C`. It is the
most precise result any run of Algorithm 1 can reach.

Results:

- `denote_mem_alphaGamma`, `alpha_sound`: `α` is sound (their Proposition 18).
- `alphaGamma_mono`: rewriting only increases precision (their Proposition 21).
- `alpha_refines_sf`: `α` refines the restricted Z-fact model `sfGamma` of `Domains.lean`, so
  every fact of that model is also a consequence of `α`.
- `alpha_cexMid`: on `cx q0,q1; h q0; h q1; cx q1,q0; h q1`, no rule applies and the wire
  map is onto, so `α` is `⊤` (every unitary).
- `alpha_pauli_incomparable`, `midE_alpha_not_pauli`: the Pauli domain and `α` are
  incomparable at every degree `d ≥ 1`, now against the paper's own abstraction.
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

namespace SF

variable {n : ℕ}

/-! ## The abstraction -/

/-- Some path of `A` starting at input `i` ends at output `o`. -/
def Reach (A : State n) (i o : Basis n) : Prop :=
  ∃ S ⊆ A.temps, ∀ q, ev (val i S) (A.ket q) = o q

/-- **γ(α(A))** (Amy and Lunderville, Definition 17 with Proposition 12): the matrices whose
nonzero transitions lie in the image of the wire map. -/
def alphaGamma (A : State n) : Set (Density n) :=
  {U | ∀ i o, U o i ≠ 0 → Reach A i o}

/-- A path sum's own operator is in its abstraction: a nonzero amplitude needs a path. -/
theorem denote_mem_alphaGamma (A : State n) : A.denote ∈ alphaGamma A := by
  intro i o h
  rw [State.denote_eq] at h
  obtain ⟨S, hS, hne⟩ := Finset.exists_ne_zero_of_sum_ne_zero (right_ne_zero_of_mul h)
  refine ⟨S, Finset.mem_powerset.mp hS, ?_⟩
  unfold State.term at hne
  split_ifs at hne with hc
  · exact hc
  · exact absurd rfl hne

/-! ## Rewriting (Algorithm 1) -/

/-- One rewrite: an HH step at substitution degree `≤ d`, or an ω step. -/
inductive Step (d : ℕ) : State n → State n → Prop
  | hh {A B : State n} {y z : ℕ} {R : Poly} {Φ₀ : Val → ℚ} : HH d A B y z R Φ₀ → Step d A B
  | omega {A B : State n} {y : ℕ} {R : Poly} {Φ₀ : Val → ℚ} : Omega A B y R Φ₀ → Step d A B

/-- Rewriting: any number of steps. -/
def Reduces (d : ℕ) : State n → State n → Prop := Relation.ReflTransGen (Step d)

theorem Step.wf_denote {d : ℕ} {A B : State n} (h : Step d A B) (hA : A.WF) :
    B.WF ∧ B.denote = A.denote := by
  cases h with
  | hh h => exact ⟨wf_hh h hA, denote_hh h hA⟩
  | omega h => exact ⟨wf_omega h hA, denote_omega h hA⟩

theorem Reduces.wf_denote {d : ℕ} {A C : State n} (h : Reduces d A C) (hA : A.WF) :
    C.WF ∧ C.denote = A.denote := by
  induction h with
  | refl => exact ⟨hA, rfl⟩
  | tail _ hst ih =>
    obtain ⟨hB, hBd⟩ := ih
    obtain ⟨hC, hCd⟩ := hst.wf_denote hB
    exact ⟨hC, hCd.trans hBd⟩

/-- A step only shrinks the set of reachable transitions. -/
theorem Step.reach {d : ℕ} {A B : State n} (h : Step d A B) (hA : A.WF) {i o : Basis n}
    (hr : Reach B i o) : Reach A i o := by
  cases h with
  | @hh y z R Φ₀ h =>
    obtain ⟨hy, hz, -, -, -, -, -, -, rfl⟩ := h
    obtain ⟨S, hS, hk⟩ := hr
    have hzn := hA.temps_ge z hz
    have hS' : S ⊆ (A.temps.erase y).erase z := hS
    have hzS : z ∉ S := fun hm => by simpa using hS' hm
    have hSA : S ⊆ A.temps :=
      hS'.trans ((Finset.erase_subset _ _).trans (Finset.erase_subset _ _))
    set r := ev (val i S) R
    have hval : val i (if r then insert z S else S) = Function.update (val i S) z r := by
      cases hr' : r
      · simp only [Bool.false_eq_true, if_false]
        rw [← val_not_mem i S z hzn hzS, Function.update_eq_self]
      · simp only [if_true]
        exact val_insert i S z hzn
    refine ⟨if r then insert z S else S, ?_, fun q => ?_⟩
    · split_ifs
      · exact Finset.insert_subset hz hSA
      · exact hSA
    · rw [hval, ← hk q]
      show _ = ev (val i S) (subst z R (A.ket q))
      rw [ev_subst]
  | @omega y R Φ₀ h =>
    obtain ⟨-, -, -, -, -, rfl⟩ := h
    obtain ⟨S, hS, hk⟩ := hr
    exact ⟨S, hS.trans (Finset.erase_subset _ _), hk⟩

/-- **Rewriting only increases precision** (their Proposition 21). -/
theorem alphaGamma_mono {d : ℕ} {A C : State n} (h : Reduces d A C) (hA : A.WF) :
    alphaGamma C ⊆ alphaGamma A := by
  induction h with
  | refl => exact le_rfl
  | tail hAB hst ih =>
    exact fun U hU i o hne => ih (fun i o hne => hst.reach (Reduces.wf_denote hAB hA).1 (hU i o hne))
      i o hne

/-! ## Traces are rewrites -/

theorem Trace.reduces {d : ℕ} {A C : State n} {f g f' g' : Poly}
    (htr : Trace d A f g C f' g') : Reduces d A C := by
  induction htr with
  | refl => exact Relation.ReflTransGen.refl
  | hh h _ _ _ ih => exact Relation.ReflTransGen.head (Step.hh h) ih
  | omega h _ _ _ ih => exact Relation.ReflTransGen.head (Step.omega h) ih

/-- A trace carrying a wire predicate ends with that wire's predicate. -/
theorem Trace.ket {d : ℕ} {A C : State n} {f g f' g' : Poly}
    (htr : Trace d A f g C f' g') (q : Fin n) (hg : g = A.ket q) : g' = C.ket q := by
  induction htr with
  | refl => exact hg
  | hh h _ _ _ ih =>
    obtain ⟨-, -, -, -, -, -, -, -, rfl⟩ := h
    exact ih (by rw [hg]; rfl)
  | omega h _ _ _ ih =>
    obtain ⟨-, -, -, -, -, rfl⟩ := h
    exact ih hg

/-- A trace carrying an input variable `x_k` leaves it evaluating to `x_k`. -/
theorem Trace.input {d : ℕ} {A C : State n} {f g f' g' : Poly}
    (htr : Trace d A f g C f' g') (hA : A.WF) {k : ℕ} (hk : k < n)
    (hf : ∀ ν, ev ν f = ν k) : ∀ ν, ev ν f' = ν k := by
  induction htr with
  | refl => exact hf
  | @hh A B C f g f' g' y z R Φ₀ h _ _ _ ih =>
    refine ih (wf_hh h hA) fun ν => ?_
    have hzk : k ≠ z := fun he => by
      have := hA.temps_ge _ h.hz; omega
    rw [ev_subst, hf, Function.update_of_ne hzk]
  | omega h _ _ _ ih => exact ih (wf_omega h hA) hf

end SF

/-! ## Circuits -/

variable {n : ℕ}

/-- **Algorithm 1's concretization at degree `d`:** the unitaries that lie in `γ(α(C))` for
every state `C` reachable by rewriting the circuit's path sum. -/
def alphaGammaC (d : ℕ) (gs : List (CTGate n)) : Set (Density n) :=
  {U | U * Uᴴ = 1 ∧ ∀ C, SF.Reduces d (SF.run gs) C → U ∈ SF.alphaGamma C}

/-- **`α` is sound** (their Proposition 18). -/
theorem alpha_sound (d : ℕ) (gs : List (CTGate n)) : circ gs ∈ alphaGammaC d gs := by
  refine ⟨SF.unitary_mul_conjTranspose gs, fun C hC => ?_⟩
  obtain ⟨-, hden⟩ := hC.wf_denote (SF.wf_run gs)
  have := SF.denote_mem_alphaGamma C
  rwa [hden, SF.denote_run] at this

/-- A unitary whose transitions all satisfy `x'_{q'} = x_q ⊕ s` maps `Z_q` to `(-1)^s Z_{q'}`. -/
theorem holds_of_support {U : Density n} (hU : U * Uᴴ = 1) {q q' : Fin n} {s : Bool}
    (h : ∀ i o, U o i ≠ 0 → o q' = (s != i q)) :
    Holds U (zString {q} false, zString {q'} s) := by
  have hc : U * (zString {q} false).toMatrix = (zString {q'} s).toMatrix * U := by
    ext o i
    rw [Matrix.mul_apply, Matrix.mul_apply,
      Finset.sum_eq_single i (fun k _ hk => by rw [zString_apply, if_neg hk, mul_zero])
        (by simp),
      Finset.sum_eq_single o (fun k _ hk => by rw [zString_apply, if_neg (Ne.symm hk), zero_mul])
        (by simp),
      zString_apply, zString_apply, if_pos rfl, if_pos rfl]
    by_cases h0 : U o i = 0
    · simp [h0]
    · have ho := h i o h0
      simp only [parity, Finset.sum_singleton, ho]
      cases s <;> cases i q <;> simp [zsgn, SPauli.sgn]
  show conj U _ = _
  rw [conj, hc, Matrix.mul_assoc, hU, Matrix.mul_one]

/-- **`α` refines the restricted model:** every Z-fact StateFold proves is a consequence of
`α` at the end of the same trace. -/
theorem alpha_refines_sf (d : ℕ) (gs : List (CTGate n)) : alphaGammaC d gs ⊆ sfGamma d gs := by
  rintro U ⟨hU, hα⟩ ⟨P, P'⟩ ⟨q, q', s, h1, h2, C, f', g', htr, heq⟩
  simp only at h1 h2
  subst h1 h2
  refine holds_of_support hU fun i o hne => ?_
  obtain ⟨S, -, hk⟩ := hα C htr.reduces i o hne
  have hg := htr.ket q' rfl
  have hf := htr.input (SF.wf_run gs) q.isLt (fun ν => by simp [SF.ev])
  rw [← hk q', ← hg, heq, hf]
  simp [SF.val]

/-! ## Incomparability with the Pauli domain -/

/-- On `cexMid`, every path variable is in a wire, so nothing rewrites. -/
theorem SF.cexMid_reduces {d : ℕ} {C : SF.State 2} (h : SF.Reduces d (SF.run SF.cexMid) C) :
    C = SF.run SF.cexMid := by
  have hexp : ∀ y ∈ (SF.run SF.cexMid).temps, ∃ q,
      ¬ SF.Supp (fun ν => SF.ev ν ((SF.run SF.cexMid).ket q))
        ((SF.run SF.cexMid).dom \ {y}) := by
    intro y hy
    rcases (SF.mid_temps y).1 hy with rfl | rfl | rfl
    · exact ⟨0, SF.not_supp_of_flip _ (by simp [SF.mid_ket0])⟩
    · exact ⟨0, SF.not_supp_of_flip _ (by simp [SF.mid_ket0])⟩
    · exact ⟨1, SF.not_supp_of_flip _ (by simp [SF.mid_ket1])⟩
  rcases Relation.ReflTransGen.cases_head h with h | ⟨B, hst, -⟩
  · exact h.symm
  · exfalso
    cases hst with
    | hh h => obtain ⟨q, hq⟩ := hexp _ h.hy; exact hq (h.ket q)
    | omega h => obtain ⟨q, hq⟩ := hexp _ h.hy; exact hq (h.ket q)

/-- On `cexMid`, `α` is `⊤`: every unitary is in it. -/
theorem alpha_cexMid (d : ℕ) (U : Density 2) (hU : U * Uᴴ = 1) : U ∈ alphaGammaC d SF.cexMid := by
  refine ⟨hU, fun C hC i o _ => ?_⟩
  rw [SF.cexMid_reduces hC]
  refine ⟨(if o 0 then {3} else ∅) ∪ (if o 1 then {4} else ∅), ?_, fun q => ?_⟩
  · intro v hv
    rw [SF.mid_temps]
    split_ifs at hv <;> simp_all
  · fin_cases q
    · rw [show ((⟨0, by decide⟩ : Fin 2)) = 0 from rfl, SF.mid_ket0]
      cases o 0 <;> cases o 1 <;> simp [SF.val]
    · rw [show ((⟨1, by decide⟩ : Fin 2)) = 1 from rfl, SF.mid_ket1]
      cases o 0 <;> cases o 1 <;> simp [SF.val]

/-- Following a circuit by X on `q` gives a unitary. -/
theorem circ_x_unitary (gs : List (CTGate n)) (q : Fin n) :
    (circ gs * gateUnitary n (CTGate.x q).toGate) *
      (circ gs * gateUnitary n (CTGate.x q).toGate)ᴴ = 1 := by
  rw [Matrix.conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc _ _ (circ gs)ᴴ,
    gate_mul_adj, Matrix.one_mul, SF.unitary_mul_conjTranspose]

/-- On a circuit starting with `h q; t q` (or `tdg`) whose `Z_q ↦ Z_q` fact StateFold proves,
the Pauli domain allows a unitary that `α` excludes. -/
theorem pauli_not_sub_alpha {d : ℕ} {q : Fin n} {g : CTGate n} (hg : g = .t q ∨ g = .tdg q)
    (rest : List (CTGate n))
    (hp : SF.Proves d (SF.run (.h q :: g :: rest)) (BoolPolynomial.var q.val) q false) :
    ¬ pauliGamma (.h q :: g :: rest) ⊆ alphaGammaC d (.h q :: g :: rest) := by
  intro h
  set U := circ (.h q :: g :: rest) * gateUnitary n (CTGate.x q).toGate
  have hmem : (zString {q} false, zString {q} false) ∈ sfFacts d (.h q :: g :: rest) :=
    ⟨q, q, false, rfl, rfl, hp⟩
  have hsf := alpha_refines_sf d _ (h (pauli_allows_flip hg rest)) _ hmem
  have htrue := sf_sound d (.h q :: g :: rest) _ hmem
  simp only [Holds] at hsf htrue
  rw [conj_mul, conj_x_z, conj_neg, htrue] at hsf
  exact zString_sign_ne {q} (by rw [show zString {q} true =
    ⟨!(zString {q} false).neg, (zString {q} false).letter⟩ from rfl, toMatrix_flip]; exact hsf)

/-- **The Pauli domain and Amy and Lunderville's `α` are incomparable**, at every degree
`d ≥ 1`.
- On `cexMid`, Pauli proves `Z₁ ↦ Z₁`, while `α` is `⊤` and contains `X₁`.
- On `h; t; tdg; h`, `α` proves `Z₀ ↦ Z₀` after one HH step, and the Pauli facts allow
  `U · X₀`. -/
theorem alpha_pauli_incomparable (d : ℕ) (hd : 1 ≤ d) :
    (∃ (n : ℕ) (gs : List (CTGate n)), ¬ alphaGammaC d gs ⊆ pauliGamma gs) ∧
    (∃ (n : ℕ) (gs : List (CTGate n)), ¬ pauliGamma gs ⊆ alphaGammaC d gs) := by
  refine ⟨⟨2, SF.cexMid, fun h => ?_⟩, ⟨1, SF.midA, ?_⟩⟩
  · have hX := h (alpha_cexMid d _ (gate_mul_adj (CTGate.x (1 : Fin 2))))
    exact x_not_holds 1 (hX _ SF.cex_pauli)
  · rw [midA_eq]
    exact pauli_not_sub_alpha (Or.inl rfl) [.tdg 0, .h 0]
      (by rw [← midA_eq]; exact SF.midA_proves d hd)

/-- The same on the 21-gate relative-Toffoli circuit, for `d ≥ 2`. -/
theorem midE_alpha_not_pauli (d : ℕ) (hd : 2 ≤ d) :
    ¬ pauliGamma SF.midE ⊆ alphaGammaC d SF.midE := by
  rw [midE_eq]
  exact pauli_not_sub_alpha (Or.inr rfl) (SF.midE.drop 2) (by
    rw [← midE_eq]
    obtain ⟨C, f', g', htr, heq⟩ := SF.midE_proves
    exact ⟨C, f', g', htr.mono hd, heq⟩)

/-! ## The auditor's regression: `α` keeps multi-wire relations -/

/-- On a single CNOT, `α` rejects the identity, since it records `x'₁ = x₀ ⊕ x₁`. The
restricted Z-fact model accepts it at every degree. -/
theorem alpha_cx_rejects_id (d : ℕ) :
    (1 : Density 2) ∉ alphaGammaC d [CTGate.cx 0 1 (by decide)] := by
  rintro ⟨-, h⟩
  obtain ⟨S, -, hk⟩ := h _ Relation.ReflTransGen.refl (fun q => q == 0) (fun q => q == 0)
    (by simp)
  have := hk 1
  simp [SF.run, SF.runFrom, SF.step, SF.init, SF.ev, SF.val] at this

end

end PauliFold
