import PauliFold.Paper.Comparison
import PauliFold.Paper.Programs
import PauliFold.StateFold.Exact

/-!
# Paper: three further claims of Sections 2.2, 7.1 and 7.2

* **Section 7.1, "How the linear abstraction works"**: "On circuits built from these gates
  alone, the formula is exact", for the gates `X`, `CX` and the diagonal gates. For every
  circuit `C` over `CTGate` without `H`, `supp ⟦C⟧` is exactly the set of pairs `(s, t)` at
  which `Φ_{ρ_C}` is satisfiable (`linear_exact`, and `linear_exact_circ` for tzap's
  semantics). The single-`CX` sentence ("the gate takes `|v₁v₂⟩` to `|v₁, v₁ ⊕ v₂⟩`, and
  nowhere else") is `linear_CX_exact`. Without the hypothesis the statement fails
  (`linear_not_exact_HH`).
* **Section 7.2, "Their abstraction of a path sum"**: the unrewritten path sum `SF.run gs` of a
  Clifford+T circuit has, on each wire, exactly the polynomial of the parity `ρ_C(q)`
  (`pathSum_ket_eq_parity`), and its path variables are exactly the Hadamard variables
  (`pathSum_temps_eq`). The variable encodings coincide: in both developments the input of
  qubit `q` is the variable `q` and the `k`-th Hadamard gate introduces the variable `n + k`
  (`pathSum_next_eq`), so the identification is the identity on variables and the map
  `parPoly` from parities to polynomials. Consequently abstracting the path sum as it
  stands gives back the linear abstraction:
  `SF.alphaGamma (SF.run gs) = γ_L(ρ) = affGamma gs` (`alphaGamma_run_eq_gammaL`,
  `alphaGamma_run_eq_affGamma`).
* **Section 2.2, dropping constraints**: after `H₁; CX_{1,2}; T₂` with `⋆ ↦ X₂`, another
  `CX_{1,2}` makes the constraint on `X₁` unknown too (`overview_drop_CX`); the dropped constraints are sound
  before and after the `CX` (`overview_drop_sound`, `overview_drop_CX_sound`).

Qubits are 0-indexed: the paper's qubits 1 and 2 are the Lean qubits `0` and `1`.
-/

namespace PauliFold.Paper

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Parities as polynomials -/

/-- The polynomial `c + ∑_{v ∈ S} v` of the parity `c ⊕ ⊕_{v ∈ S} v`, over the same variables
(input `q` is the variable `q`, Hadamard variables are `n, n + 1, …`). -/
def parPoly (p : Parity) : SF.Poly :=
  MvPolynomial.C p.2 + Finsupp.linearCombination (ZMod 2) (fun v => (MvPolynomial.X v : SF.Poly)) p.1

theorem parPoly_add (p p' : Parity) : parPoly (p + p') = parPoly p + parPoly p' := by
  simp only [parPoly, Prod.snd_add, Prod.fst_add, map_add]
  abel

theorem parPoly_var (v : ℕ) : parPoly (Finsupp.single v 1, (0 : ZMod 2)) =
    BoolPolynomial.var v := by
  simp [parPoly, BoolPolynomial.var]

theorem parPoly_one : parPoly ((0 : Lin), (1 : ZMod 2)) = 1 := by
  simp [parPoly]

theorem bit_eq_tzap (b : Bool) : Paper.bit b = TzapLean.bit b := rfl

/-- Evaluating the polynomial of a parity is evaluating the parity. -/
theorem evalF₂_parPoly (ν : SF.Val) (p : Parity) :
    BoolPolynomial.evalF₂ ν (parPoly p) = Parity.eval (fun v => TzapLean.bit (ν v)) p := by
  simp only [BoolPolynomial.evalF₂, parPoly, Parity.eval, map_add, MvPolynomial.eval_C,
    Finsupp.linearCombination_apply, map_finsuppSum, MvPolynomial.smul_eval, MvPolynomial.eval_X,
    smul_eq_mul]

theorem ev_parPoly (ν : SF.Val) (p : Parity) :
    SF.ev ν (parPoly p) = decide (Parity.eval (fun v => TzapLean.bit (ν v)) p = 1) := by
  rw [SF.ev, BoolPolynomial.evalB, evalF₂_parPoly]; rfl

/-- A parity only depends on the variables it mentions. -/
theorem Parity.eval_congr {ν μ : ℕ → ZMod 2} (p : Parity) (h : ∀ v, p.1 v ≠ 0 → ν v = μ v) :
    Parity.eval ν p = Parity.eval μ p := by
  simp only [Parity.eval, Finsupp.linearCombination_apply]
  congr 1
  refine Finsupp.sum_congr fun v hv => ?_
  rw [h v (Finsupp.mem_support_iff.1 hv)]

/-! ## The unrewritten path sum is the linear abstraction (Section 7.2) -/

/-- The correspondence between a path-sum state and a linear abstract state: each wire's
polynomial is the polynomial of the wire's parity, the counters of fresh variables agree, and
the path variables are exactly the Hadamard variables `n, …, fresh - 1`. -/
structure Corr (A : SF.State n) (a : AffAbs n) : Prop where
  ket : ∀ q, A.ket q = parPoly (a.wires q)
  next : A.next = a.fresh
  temps : A.temps = Finset.Ico n A.next
  n_le : n ≤ A.next

theorem corr_init (n : ℕ) : Corr (SF.init n) (AffAbs.init n) where
  ket q := (parPoly_var q.val).symm
  next := rfl
  temps := (Finset.Ico_self n).symm
  n_le := le_rfl

theorem corr_step (j : ℕ) (g : CTGate n) {A : SF.State n} {a : AffAbs n} (h : Corr A a) :
    Corr (SF.step j g A) (AffAbs.step g a) := by
  cases g with
  | x q =>
    refine ⟨fun r => ?_, h.next, h.temps, h.n_le⟩
    show Function.update A.ket q (A.ket q + 1) r = _
    rw [LinState.step_x]
    by_cases hr : r = q
    · subst hr; simp [h.ket, parPoly_add, parPoly_one]
    · simp [hr, h.ket]
  | cx c t hne =>
    refine ⟨fun r => ?_, h.next, h.temps, h.n_le⟩
    show Function.update A.ket t (A.ket t + A.ket c) r = _
    rw [LinState.step_cx]
    by_cases hr : r = t
    · subst hr; simp [h.ket, parPoly_add]
    · simp [hr, h.ket]
  | h q =>
    obtain ⟨hw, hf⟩ := LinState.step_h a q
    refine ⟨fun r => ?_, ?_, ?_, ?_⟩
    · show Function.update A.ket q (BoolPolynomial.var A.next) r = _
      rw [hw]
      by_cases hr : r = q
      · subst hr; simp [h.next, parPoly_var]
      · simp [hr, h.ket]
    · show A.next + 1 = _
      rw [hf, h.next]
    · show insert A.next A.temps = Finset.Ico n (A.next + 1)
      rw [h.temps, Nat.Ico_succ_right_eq_insert_Ico h.n_le]
    · show n ≤ A.next + 1
      exact h.n_le.trans (Nat.le_succ _)
  | z q => exact ⟨h.ket, h.next, h.temps, h.n_le⟩
  | s q => exact ⟨h.ket, h.next, h.temps, h.n_le⟩
  | sdg q => exact ⟨h.ket, h.next, h.temps, h.n_le⟩
  | t q => exact ⟨h.ket, h.next, h.temps, h.n_le⟩
  | tdg q => exact ⟨h.ket, h.next, h.temps, h.n_le⟩
  | cz c t hne => exact ⟨h.ket, h.next, h.temps, h.n_le⟩

theorem corr_runFrom : ∀ (gs : List (CTGate n)) (j : ℕ) {A : SF.State n} {a : AffAbs n},
    Corr A a → Corr (SF.runFrom j gs A) (a.run gs)
  | [], _, _, _, h => h
  | g :: gs, j, _, _, h => corr_runFrom gs (j + 1) (corr_step j g h)

theorem corr_run (gs : List (CTGate n)) : Corr (SF.run gs) ((AffAbs.init n).run gs) :=
  corr_runFrom gs 0 (corr_init n)

/-- **Paper, Section 7.2 ("Their abstraction of a path sum")**: "for a Clifford+T circuit, each
output bit `f_q` is exactly the parity `ρ_C(q)` of the linear abstraction". The output
polynomial of wire `q` of the unrewritten path sum is, syntactically, the polynomial of the
parity the linear abstraction gives wire `q`. -/
theorem pathSum_ket_eq_parity (gs : List (CTGate n)) (q : Fin n) :
    (SF.run gs).ket q = parPoly (((AffAbs.init n).run gs).wires q) :=
  (corr_run gs).ket q

/-- Paper, Section 7.2: "with the path variables as the Hadamard variables". Both developments
number variables in the same way, so the next fresh variable agrees ... -/
theorem pathSum_next_eq (gs : List (CTGate n)) :
    (SF.run gs).next = ((AffAbs.init n).run gs).fresh :=
  (corr_run gs).next

/-- ... and the path variables of the path sum are exactly the Hadamard variables
`n, …, fresh - 1` of the linear abstraction. -/
theorem pathSum_temps_eq (gs : List (CTGate n)) :
    (SF.run gs).temps = Finset.Ico n ((AffAbs.init n).run gs).fresh := by
  rw [(corr_run gs).temps, pathSum_next_eq]

/-- The same, for the circuits `CTCirc` of the paper and `ρ_C`. -/
theorem CTCirc.pathSum_ket_eq_rho (C : CTCirc n) (q : Fin n) :
    (SF.run C.toList).ket q = parPoly (C.rho.wires q) ∧
      (SF.run C.toList).temps = Finset.Ico n C.rho.fresh :=
  ⟨pathSum_ket_eq_parity _ q, pathSum_temps_eq _⟩

/-- A path of the unrewritten path sum from `s` to `t` is a solution of `Φ_ρ` at `(s, t)`. -/
theorem reach_iff_satAt (gs : List (CTGate n)) (s t : Basis n) :
    SF.Reach (SF.run gs) s t ↔ LinState.SatAt ((AffAbs.init n).run gs) s t := by
  have hc := corr_run gs
  set ρ := (AffAbs.init n).run gs
  have hbit : ∀ x : ZMod 2, TzapLean.bit (decide (x = 1)) = x := fun x => bit_unbit x
  constructor
  · rintro ⟨S, -, hk⟩
    refine ⟨fun v => TzapLean.bit (decide (v ∈ S)), fun q => ?_⟩
    have hass : assign s (fun v => TzapLean.bit (decide (v ∈ S))) =
        fun v => TzapLean.bit (SF.val s S v) := by
      funext v
      by_cases hv : v < n
      · simp [assign, SF.val, hv, bit_eq_tzap]
      · simp [assign, SF.val, hv]
    have e := hk q
    rw [hc.ket, ev_parPoly] at e
    show Paper.bit (t q) = _
    rw [hass, bit_eq_tzap, ← e, hbit]
  · rintro ⟨u, hu⟩
    let S := (Finset.Ico n ρ.fresh).filter (fun v => u v = 1)
    refine ⟨S, ?_, fun q => ?_⟩
    · rw [hc.temps, hc.next]; exact Finset.filter_subset _ _
    · rw [hc.ket, ev_parPoly]
      have hq := hu q
      simp only at hq
      rw [Parity.eval_congr (μ := assign s u) _ ?_, ← hq, bit_eq_tzap]
      · exact unbit_bit (t q)
      intro v hv
      by_cases hvn : v < n
      · simp [assign, SF.val, hvn, bit_eq_tzap]
      · have hvf : v < ρ.fresh := by
          by_contra hge
          exact hv ((LinState.fresh_ge gs).2 q v (show ρ.fresh ≤ v by omega))
        have hvS : v ∈ S ↔ u v = 1 := by
          simp only [S, Finset.mem_filter, Finset.mem_Ico]
          exact ⟨fun h => h.2, fun h => ⟨⟨by omega, hvf⟩, h⟩⟩
        simp only [assign, SF.val, hvn, dite_false]
        rcases zmod2_cases (u v) with h0 | h1
        · have : v ∉ S := fun h => by rw [hvS, h0] at h; exact zmod2_one_ne_zero h.symm
          simp [this, h0, TzapLean.bit]
        · simp [hvS.2 h1, h1, TzapLean.bit]

/-- **Paper, Section 7.2: "abstracting the path sum as it stands gives back the linear
abstraction."** Amy and Lunderville's concretization of the unrewritten path sum of a
Clifford+T circuit (`SF.alphaGamma`, their Definition 17 with Proposition 12) is the paper's
`γ_L(ρ_C)` (Definition 7.3). -/
theorem alphaGamma_run_eq_gammaL (gs : List (CTGate n)) :
    SF.alphaGamma (SF.run gs) = gammaL ((AffAbs.init n).run gs) := by
  ext K
  simp only [SF.alphaGamma, gammaL, supp, Set.mem_ofPred_eq, Prod.forall, reach_iff_satAt]

/-- The same, against the earlier `affGamma`. -/
theorem alphaGamma_run_eq_affGamma (gs : List (CTGate n)) :
    SF.alphaGamma (SF.run gs) = affGamma gs := by
  rw [alphaGamma_run_eq_gammaL, gammaL_eq_affGamma]

/-- The same, for `CTCirc` and `ρ_C`. -/
theorem CTCirc.alphaGamma_run_eq_gammaL (C : CTCirc n) :
    SF.alphaGamma (SF.run C.toList) = gammaL C.rho :=
  Paper.alphaGamma_run_eq_gammaL _

/-! ## Without Hadamard gates, the linear formula is exact (Section 7.1) -/

/-- A gate other than `H`. -/
def NotH (g : CTGate n) : Prop := ∀ q, g ≠ .h q

/-- A circuit without Hadamard gates. -/
def CTCirc.NoH (C : CTCirc n) : Prop := ∀ g ∈ C.toList, NotH g

theorem fresh_run_of_notH : ∀ (gs : List (CTGate n)) (a : AffAbs n), (∀ g ∈ gs, NotH g) →
    (a.run gs).fresh = a.fresh
  | [], _, _ => rfl
  | g :: gs, a, h => by
    have hg := h g List.mem_cons_self
    show ((AffAbs.step g a).run gs).fresh = a.fresh
    rw [fresh_run_of_notH gs _ fun g' hg' => h g' (List.mem_cons_of_mem _ hg')]
    cases g with
    | h q => exact (hg q rfl).elim
    | _ => rfl

/-- Without Hadamard gates the path sum has no path variables, so its support is the set of
transitions it reaches. -/
theorem supp_circ_iff_reach (gs : List (CTGate n)) (hgs : ∀ g ∈ gs, NotH g) (s t : Basis n) :
    (s, t) ∈ supp (circ gs) ↔ SF.Reach (SF.run gs) s t := by
  have hT : (SF.run gs).temps = ∅ := by
    rw [pathSum_temps_eq, fresh_run_of_notH gs _ hgs]
    exact Finset.Ico_self n
  show circ gs t s ≠ 0 ↔ _
  rw [circ, ← SF.denote_run, SF.State.denote_eq, hT]
  simp only [Finset.card_empty, pow_zero, one_mul, Finset.powerset_empty, Finset.sum_singleton,
    SF.State.term]
  constructor
  · intro h
    split_ifs at h with hc
    · exact ⟨∅, Finset.empty_subset _, hc⟩
    · exact absurd rfl h
  · rintro ⟨S, hS, hk⟩
    rw [hT, Finset.subset_empty] at hS
    subst hS
    rw [if_pos hk]
    exact ep_ne_zero _

/-- **Paper, Section 7.1: "On circuits built from these gates alone, the formula is exact"**, for
tzap's semantics: for a gate list without `H` (only `X`, `CX` and the diagonal gates `Z, S, S†,
T, T†, CZ`), `⟨t|⟦gs⟧|s⟩ ≠ 0` exactly when `Φ_ρ[v := s, y := t]` is satisfiable. -/
theorem linear_exact_circ (gs : List (CTGate n)) (hgs : ∀ g ∈ gs, NotH g) (s t : Basis n) :
    (s, t) ∈ supp (circ gs) ↔ LinState.SatAt ((AffAbs.init n).run gs) s t := by
  rw [supp_circ_iff_reach gs hgs, reach_iff_satAt]

/-- The support is unchanged by a global phase. -/
theorem supp_of_equivP {U U' : Density n} (h : EquivP U U') : supp U = supp U' := by
  obtain ⟨φ, rfl⟩ := h
  ext st
  simp only [supp, Set.mem_ofPred_eq, Matrix.smul_apply, smul_eq_mul]
  exact ⟨fun h hz => h (by rw [hz, mul_zero]), fun h => mul_ne_zero (Complex.exp_ne_zero _) h⟩

/-- **Paper, Section 7.1: "On circuits built from these gates alone, the formula is exact."** For
every circuit `C` over `X, CX` and the diagonal gates `Z, S, S†, T, T†, CZ` (no `H`), the
support of `⟦C⟧` is exactly the relation of the formula `Φ_{ρ_C}`:
`(s, t) ∈ supp(⟦C⟧) ↔ Φ_{ρ_C}[v_q := s_q, y_q := t_q]_q` is satisfiable. (With no Hadamard
variables, satisfiable just means true.) -/
theorem linear_exact (C : CTCirc n) (hC : C.NoH) :
    supp C.toPaper.sem = {st | C.rho.SatAt st.1 st.2} := by
  rw [← supp_of_equivP (CTCirc.sem_toPaper C)]
  ext ⟨s, t⟩
  exact linear_exact_circ C.toList hC s t

/-- With no Hadamard variables the formula needs no witness: `Φ_{ρ_C}` holds at `(s, t)` for
some value of the Hadamard variables iff it holds for all of them. -/
theorem satAt_iff_forall_of_noH (C : CTCirc n) (hC : C.NoH) (s t : Basis n) :
    C.rho.SatAt s t ↔ ∀ u, C.rho.Phi (assign s u) (fun q => Paper.bit (t q)) := by
  have hf : C.rho.fresh = n := fresh_run_of_notH _ _ hC
  have key : ∀ u u' q, Parity.eval (assign s u) (C.rho.wires q) =
      Parity.eval (assign s u') (C.rho.wires q) := by
    intro u u' q
    refine Parity.eval_congr _ fun v hv => ?_
    by_cases hvn : v < n
    · simp [assign, hvn]
    · exact absurd ((LinState.fresh_ge C.toList).2 q v (by rw [CTCirc.rho] at hf; omega)) hv
  constructor
  · rintro ⟨u, hu⟩ u' q
    rw [hu q, key]
  · intro h; exact ⟨0, h 0⟩

/-- **Paper, Section 7.1**: "This formula is exactly the transition relation of the `CX`: the gate
takes the basis state `|v₁v₂⟩` to `|v₁, v₁ ⊕ v₂⟩`, and nowhere else." (Lean qubits `0, 1`.) -/
theorem linear_CX_exact (s t : Basis 2) :
    (s, t) ∈ supp (CTCirc.gate (.cx (0 : Fin 2) 1 (by decide))).toPaper.sem ↔
      t 0 = s 0 ∧ t 1 = (s 0 != s 1) := by
  rw [linear_exact _ (by
    intro g hg q
    simp only [CTCirc.toList_gate, List.mem_singleton] at hg
    subst hg; intro h; cases h)]
  obtain ⟨h0, h1⟩ := linear_CX
  simp only [Set.mem_ofPred_eq, LinState.SatAt, LinState.Phi, Fin.forall_fin_two, h0, h1]
  simp only [Parity.eval, map_add, Finsupp.linearCombination_single, one_smul,
    zero_add]
  have e0 : assign s (0 : ℕ → ZMod 2) 0 = Paper.bit (s 0) := assign_lt s 0 0
  have e1 : assign s (0 : ℕ → ZMod 2) 1 = Paper.bit (s 1) := assign_lt s 0 1
  have hinj : ∀ a b : Bool, Paper.bit a = Paper.bit b ↔ a = b := by decide
  have hx : ∀ a b : Bool, Paper.bit (a != b) = Paper.bit a + Paper.bit b := TzapLean.bit_xor
  constructor
  · rintro ⟨u, hu0, hu1⟩
    have e0' : assign s u 0 = Paper.bit (s 0) := assign_lt s u 0
    have e1' : assign s u 1 = Paper.bit (s 1) := assign_lt s u 1
    rw [e0'] at hu0
    rw [e0', e1', add_comm, ← hx] at hu1
    exact ⟨(hinj _ _).1 hu0, (hinj _ _).1 hu1⟩
  · rintro ⟨ht0, ht1⟩
    refine ⟨0, ?_, ?_⟩
    · rw [e0, ht0]
    · rw [e0, e1, ht1, hx, add_comm]

/-- The hypothesis cannot be dropped: on `H₁; H₁` (the identity) the formula
`y₁ = u₂` relates `|0⟩` to `|1⟩`, which the circuit does not. -/
theorem linear_not_exact_HH :
    HH.rho.SatAt (fun _ => false) (fun _ => true) ∧
      ((fun _ => false, fun _ => true) : Basis 1 × Basis 1) ∉ supp HH.toPaper.sem := by
  refine ⟨?_, ?_⟩
  · refine ⟨fun _ => 1, fun q => ?_⟩
    rw [Fin.fin_one_eq_zero q, gammaL_HH.1]
    simp [Parity.eval, assign, Paper.bit]
  · show ¬ (HH.toPaper.sem (fun _ => true) (fun _ => false) ≠ 0)
    rw [HH_toPaper, Circuit.sem_seq]
    simp only [Circuit.sem, Circuit.gH_mul_self, Matrix.one_apply, not_not]
    rw [if_neg]
    intro h
    have := congrFun h 0
    simp at this

/-! ## Section 2.2: unknowns are contagious -/

namespace Examples

/-- Paper, Section 2.2, cut 2 of `H₁; CX_{1,2}; T₂` with the constraint `T₂` breaks dropped:
the tableau of `H₁; CX_{1,2}` with `⋆ ↦ X₂` (Lean `x 1 = none`):
`Z₁X₂ ↦ X₁, X₁ ↦ Z₁, ⋆ ↦ X₂, X₁Z₂ ↦ Z₂`. -/
def consDrop : Cons 2 := ⟨![some Z1X2, none], ![some (Xq 0), some X1Z2]⟩

theorem consDrop_eq : consDrop = { exCons with x := Function.update exCons.x 1 none } := by
  decide

/-- The constraints after one more `CX_{1,2}`: `⋆ ↦ X₁, X₁ ↦ Z₁, ⋆ ↦ X₂, Z₂ ↦ Z₂`. -/
def consDropCX : Cons 2 := ⟨![none, none], ![some (Xq 0), some (Zq 1)]⟩

/-- **Paper, Section 2.2**: "a later `CX_{1,2}` mixes the constraint on `X₂` into the one on `X₁`,
making it unknown too." Composing the dropped constraints of cut 2 with `α(CX_{1,2})` (Lean `CX_{0,1}`): the
`CX` has `X₁X₂ ↦ X₁`, whose weakest precondition needs the unknown constraint on `X₂`, so the
constraint on `X₁` becomes `⋆ ↦ X₁` (Lean `x 0 = none`); the constraint on `X₂` stays unknown, and
the constraints on `Z₁` and `Z₂` stay known (`X₁ ↦ Z₁`, and `Z₂ ↦ Z₂` since
`wp(Z₁Z₂) = X₁ · X₁Z₂ = Z₂`). Without the drop, the constraint on `X₁` would be `Z₁ ↦ X₁` (since
`CX_{1,2}; CX_{1,2}` is the identity). At the level of abstract relations, the composition is
`(consDropCX, ∅)`. -/
theorem overview_drop_CX :
    (gCX (0 : Fin 2) 1 (by decide)).cons.x 0 = some ((Xq 0).mul (Xq 1)) ∧
      ((gCX (0 : Fin 2) 1 (by decide)).cons.comp consDrop).x 0 = none ∧
      (gCX (0 : Fin 2) 1 (by decide)).cons.comp consDrop = consDropCX ∧
      ((gCX (0 : Fin 2) 1 (by decide)).cons.comp exCons).x 0 = some (Zq 0) ∧
      Circuit.alpha (.gate (gCX (0 : Fin 2) 1 (by decide))) ⊚ .rel ⟨consDrop, []⟩ =
        .rel ⟨consDropCX, []⟩ := by
  decide

/-- `consDrop` only forgets constraints of the tableau of `H₁; CX_{1,2}`. -/
theorem exCons_le_consDrop : exCons.Le consDrop := by decide

/-- **Paper, Section 2.2**: forgetting the broken constraint (`⋆ ↦ X₂`) is sound: `⟦H₁; CX_{1,2}; T₂⟧` is in the
concretization of the dropped constraints (no axes). A string tracked by `consDrop` has `I` or
`Z` on qubit 2 (Lean `1`), since `⋆ ↦ X₂`, so it commutes with `T₂`. -/
theorem overview_drop_sound :
    Circuit.sem exCirc ∈ (Dom.rel ⟨consDrop, []⟩ : Dom 2).gamma := by
  have hs := Circuit.soundness exCirc
  rw [example_4_3] at hs
  have hV := example_4_3_model
  intro P B hB _
  have hB' := Cons.wp_of_le exCons_le_consDrop hB
  refine hs P B hB' fun A hA => ?_
  have hA' : A = X1Z2 := List.mem_singleton.1 hA
  subst hA'
  have hZ : exCons.wp (Zq 1) = some X1Z2 := example_4_12_T.2
  rw [PP.comm_iff, Cons.wp_toMatrix hV hB', Cons.wp_toMatrix hV hZ, comm_conj_iff hV.unitary,
    ← PP.comm_iff, PP.comm_single]
  have hk := ((Cons.wp_isSome_iff consDrop P).1 (by rw [hB]; rfl)) 1
  rw [Cons.w_isSome_iff] at hk
  cases hL : P.letter 1 <;> rw [hL] at hk <;> first | rfl | simp [consDrop] at hk

/-- **Paper, Section 2.2**, continued: the composition with `α(CX_{1,2})` is sound too,
`⟦H₁; CX_{1,2}; T₂; CX_{1,2}⟧ ∈ γ(consDropCX, ∅)`. -/
theorem overview_drop_CX_sound :
    Circuit.sem (exCirc ⨾ .gate (gCX (0 : Fin 2) 1 (by decide))) ∈
      (Dom.rel ⟨consDropCX, []⟩ : Dom 2).gamma := by
  rw [← overview_drop_CX.2.2.2.2, Circuit.sem_seq]
  exact gamma_comp (σ₁ := .rel ⟨consDrop, []⟩)
    (show consDrop.Sat from ⟨_, example_4_3_model.of_le exCons_le_consDrop⟩) overview_drop_sound
    (Circuit.gate_mem_gamma _)

end Examples

end

end PauliFold.Paper
