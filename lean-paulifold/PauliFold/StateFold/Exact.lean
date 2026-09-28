import PauliFold.StateFold.Basic

/-!
# Exactness of the StateFold transformers

Every gate transformer is exact: `denote (step j g A) = U_g · denote A` for well-formed `A`.
So the state after a circuit denotes the circuit's unitary (`denote_run`).
-/

namespace PauliFold.SF

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Evaluation -/

@[simp] theorem ev_add (ν : Val) (p q : Poly) : ev ν (p + q) = (ev ν p != ev ν q) :=
  BoolPolynomial.evalB_add ν p q

@[simp] theorem ev_mul (ν : Val) (p q : Poly) : ev ν (p * q) = (ev ν p && ev ν q) :=
  BoolPolynomial.evalB_mul ν p q

@[simp] theorem ev_var (ν : Val) (i : ℕ) : ev ν (BoolPolynomial.var i) = ν i :=
  BoolPolynomial.evalB_var ν i

@[simp] theorem ev_one (ν : Val) : ev ν (1 : Poly) = true := by
  show BoolPolynomial.evalB ν 1 = true
  simpa [BoolPolynomial.const, bit] using BoolPolynomial.evalB_const ν true

theorem ep_b2q (θ : ℚ) (b : Bool) : ep (θ * b2q b) = if b then ep θ else 1 := by
  cases b <;> simp [b2q]

/-! ## Congruence -/

/-- Two states with the same path variables and phase, whose wire conditions agree for every
assignment, have the same entry. -/
theorem denote_congr (A B : State n) (ht : A.temps = B.temps) (hp : A.phase = B.phase)
    (o o' i : Basis n)
    (hc : ∀ S, (∀ q, ev (val i S) (A.ket q) = o q) ↔ (∀ q, ev (val i S) (B.ket q) = o' q)) :
    A.denote o i = B.denote o' i := by
  unfold State.denote
  rw [ht, hp]
  congr 1
  refine Finset.sum_congr rfl fun S _ => ?_
  simp only [hc S]

/-! ## X and CX -/

theorem conds_x (K : Fin n → Poly) (q : Fin n) (o : Basis n) (ν : Val) :
    (∀ r, ev ν (K r) = Function.update o q (!o q) r) ↔
      (∀ r, ev ν (Function.update K q (K q + 1) r) = o r) := by
  constructor
  · intro h r
    by_cases hr : r = q
    · subst hr
      have h1 := h r
      rw [Function.update_self] at h1
      rw [Function.update_self, ev_add, ev_one, h1]
      cases o r <;> rfl
    · have h1 := h r
      rwa [Function.update_of_ne hr] at h1 ⊢
  · intro h r
    by_cases hr : r = q
    · subst hr
      have h1 := h r
      rw [Function.update_self, ev_add, ev_one] at h1
      rw [Function.update_self, ← h1]
      cases ev ν (K r) <;> rfl
    · have h1 := h r
      rwa [Function.update_of_ne hr] at h1 ⊢

theorem embed1_x2_mul_apply (q : Fin n) (B : Density n) (o i : Basis n) :
    (embed1 n x2 q.val * B) o i = B (Function.update o q (!o q)) i := by
  rw [embed1_mul_apply, Fintype.sum_bool]
  cases hq : o q <;> simp only [x2, Bool.not_true, Bool.not_false] <;> norm_num

theorem denote_step_x (j : ℕ) (q : Fin n) (A : State n) :
    (step j (.x q) A).denote = gateUnitary n (CTGate.x q).toGate * A.denote := by
  funext o i
  show _ = (embed1 n x2 q.val * A.denote) o i
  rw [embed1_x2_mul_apply]
  exact (denote_congr A (step j (.x q) A) rfl rfl (Function.update o q (!o q)) o i
    fun S => conds_x A.ket q o (val i S)).symm

theorem conds_cx (K : Fin n → Poly) (c t : Fin n) (hne : c ≠ t) (o : Basis n) (ν : Val) :
    (∀ r, ev ν (K r) = cxPerm c t o r) ↔
      (∀ r, ev ν (Function.update K t (K t + K c) r) = o r) := by
  constructor
  · intro h r
    have hc := h c
    simp only [cxPerm, Function.update_of_ne hne] at hc
    by_cases hr : r = t
    · subst hr
      have ht := h r
      simp only [cxPerm, Function.update_self] at ht
      simp only [Function.update_self, ev_add, ht, hc]
      cases o r <;> cases o c <;> rfl
    · have := h r
      simp_all [cxPerm, Function.update_of_ne hr]
  · intro h r
    have hc := h c
    simp only [Function.update_of_ne hne] at hc
    by_cases hr : r = t
    · subst hr
      have ht := h r
      simp only [Function.update_self, ev_add, hc] at ht
      simp only [cxPerm, Function.update_self, ← ht]
      cases ev ν (K r) <;> cases o c <;> rfl
    · have := h r
      simp_all [cxPerm, Function.update_of_ne hr]

theorem denote_step_cx (j : ℕ) (c t : Fin n) (hne : c ≠ t) (A : State n) :
    (step j (.cx c t hne) A).denote = gateUnitary n (CTGate.cx c t hne).toGate * A.denote := by
  rw [gateUnitary_cx]
  funext o i
  rw [permMatrix_mul_apply _ (cxPerm_involutive c t hne)]
  exact (denote_congr A (step j (.cx c t hne) A) rfl rfl (cxPerm c t o) o i
    fun S => conds_cx A.ket c t hne o (val i S)).symm

/-! ## Diagonal gates -/

/-- Multiplying every term by a factor fixed by the wire condition. -/
theorem denote_phase_mul (A : State n) (δ : Val → ℚ) (f : Basis n → ℂ) (o i : Basis n)
    (hδ : ∀ S, (∀ q, ev (val i S) (A.ket q) = o q) → ep (δ (val i S)) = f o) :
    ({ A with phase := fun ν => A.phase ν + δ ν } : State n).denote o i =
      f o * A.denote o i := by
  show invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset,
      (if ∀ q, ev (val i S) (A.ket q) = o q then ep (A.phase (val i S) + δ (val i S)) else 0) =
    f o * (invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset,
      (if ∀ q, ev (val i S) (A.ket q) = o q then ep (A.phase (val i S)) else 0))
  rw [Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun S _ => ?_
  split_ifs with h
  · rw [ep_add, hδ S h]; ring
  · simp

theorem embed1_diag_mul_apply (θ : ℚ) (q : Fin n) (B : Density n) (o i : Basis n) :
    (embed1 n (diag2 1 (ep θ)) q.val * B) o i = (if o q then ep θ else 1) * B o i := by
  rw [embed1_mul_apply, Fintype.sum_bool]
  have hself : Function.update o q (o q) = o := Function.update_eq_self q o
  cases hq : o q <;> rw [hq] at hself <;> simp [diag2, hself]

theorem denote_rot (j : ℕ) (θ : ℚ) (q : Fin n) (A : State n) :
    (rot j θ q A).denote = embed1 n (diag2 1 (ep θ)) q.val * A.denote := by
  funext o i
  rw [embed1_diag_mul_apply]
  unfold rot
  have := denote_phase_mul A (fun ν => θ * b2q (ev ν (A.ket q))) (fun o => if o q then ep θ else 1)
    o i (fun S h => by rw [ep_b2q, h q])
  simpa [State.denote] using this

theorem denote_step_cz (j : ℕ) (c t : Fin n) (hne : c ≠ t) (A : State n) :
    (step j (.cz c t hne) A).denote = gateUnitary n (CTGate.cz c t hne).toGate * A.denote := by
  rw [gateUnitary_cz]
  funext o i
  rw [phaseMatrix_mul_apply]
  exact denote_phase_mul A _ (czPhase c t) o i fun S h => by
    rw [h c, h t]
    simp only [czPhase, b2q]
    split_ifs <;> simp [ep_one]

/-! ## Hadamard -/

theorem val_insert (i : Basis n) (S : Finset ℕ) (y : ℕ) (hy : n ≤ y) :
    val i (insert y S) = Function.update (val i S) y true := by
  funext v
  by_cases hv : v = y
  · subst hv; simp [val, show ¬ v < n by omega]
  · simp [val, Function.update, hv]

theorem val_not_mem (i : Basis n) (S : Finset ℕ) (y : ℕ) (hy : n ≤ y) (hS : y ∉ S) :
    val i S y = false := by
  simp [val, show ¬ y < n by omega, hS]

/-- A function supported on `D` ignores a variable outside `D`. -/
theorem Supp.update {α : Type*} {f : Val → α} {D : Set ℕ} (h : Supp f D) {y : ℕ} (hy : y ∉ D)
    (ν : Val) (b : Bool) : f (Function.update ν y b) = f ν :=
  h _ _ fun v hv => by
    have : v ≠ y := fun e => hy (e ▸ hv)
    simp [Function.update, this]

/-- The wire condition after an H on `q`, whose new variable `y` is set to `b`. -/
theorem conds_h (A : State n) (hA : A.WF) (q : Fin n) (ν : Val) (b : Bool) (o : Basis n)
    (hyD : A.next ∉ A.dom) :
    (∀ r, ev (Function.update ν A.next b) (Function.update A.ket q (BoolPolynomial.var A.next) r) =
        o r) ↔ (b = o q ∧ ∀ r, r ≠ q → ev ν (A.ket r) = o r) := by
  constructor
  · intro h
    refine ⟨?_, fun r hr => ?_⟩
    · have := h q; simpa using this
    · have := h r
      rw [Function.update_of_ne hr, (hA.ket r).update hyD] at this
      exact this
  · rintro ⟨hb, h⟩ r
    by_cases hr : r = q
    · subst hr; simpa using hb
    · rw [Function.update_of_ne hr, (hA.ket r).update hyD]; exact h r hr

/-- The wire condition for `o` updated at `q`. -/
theorem conds_update (K : Fin n → Poly) (q : Fin n) (ν : Val) (α : Bool) (o : Basis n) :
    (∀ r, ev ν (K r) = Function.update o q α r) ↔
      (ev ν (K q) = α ∧ ∀ r, r ≠ q → ev ν (K r) = o r) := by
  constructor
  · intro h
    refine ⟨by have := h q; simpa using this, fun r hr => ?_⟩
    have := h r; rwa [Function.update_of_ne hr] at this
  · rintro ⟨hf, h⟩ r
    by_cases hr : r = q
    · subst hr; simpa using hf
    · rw [Function.update_of_ne hr]; exact h r hr

theorem h2_eq (a b : Bool) : h2 a b = invSqrt2 * (if a && b then -1 else 1) := by
  simp only [h2, invSqrt2, div_eq_mul_inv, mul_comm]

theorem denote_step_h (j : ℕ) (q : Fin n) (A : State n) (hA : A.WF) :
    (step j (.h q) A).denote = gateUnitary n (CTGate.h q).toGate * A.denote := by
  funext o i
  show _ = (embed1 n h2 q.val * A.denote) o i
  have hyT : A.next ∉ A.temps := fun h => absurd (hA.temps_lt _ h) (lt_irrefl _)
  have hyD : A.next ∉ A.dom := by
    rintro (h | h)
    · exact absurd hA.n_le (by omega)
    · exact hyT h
  have hyn : n ≤ A.next := hA.n_le
  -- The common normal form: `K · ∑_S G S`.
  let rest : Finset ℕ → Prop := fun S => ∀ r, r ≠ q → ev (val i S) (A.ket r) = o r
  let f : Finset ℕ → Bool := fun S => ev (val i S) (A.ket q)
  let G : Finset ℕ → ℂ := fun S => by
    classical
    exact if rest S then invSqrt2 * (if o q && f S then -1 else 1) * ep (A.phase (val i S)) else 0
  have lhs : (step j (.h q) A).denote o i =
      invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset, G S := by
    simp only [step, State.denote]
    rw [Finset.card_insert_of_notMem hyT, Finset.sum_powerset_insert hyT, pow_succ,
      ← Finset.sum_add_distrib, mul_assoc, Finset.mul_sum]
    congr 1
    refine Finset.sum_congr rfl fun S hS => ?_
    have hyS : A.next ∉ S := fun h => hyT (Finset.mem_powerset.mp hS h)
    have hνy := val_not_mem i S A.next hyn hyS
    have h0 : Function.update (val i S) A.next false = val i S := by
      rw [← hνy]; exact Function.update_eq_self _ _
    have key0 : (∀ r, ev (val i S) (Function.update A.ket q (BoolPolynomial.var A.next) r) =
        o r) ↔ (false = o q ∧ rest S) := by
      have := conds_h A hA q (val i S) false o hyD
      rwa [h0] at this
    have key1 : (∀ r, ev (val i (insert A.next S))
        (Function.update A.ket q (BoolPolynomial.var A.next) r) = o r) ↔
          (true = o q ∧ rest S) := by
      rw [val_insert i S A.next hyn]
      exact conds_h A hA q (val i S) true o hyD
    have hq1 : ev (val i (insert A.next S)) (A.ket q) = f S := by
      rw [val_insert i S A.next hyn]; exact (hA.ket q).update hyD _ _
    have hp1 : A.phase (val i (insert A.next S)) = A.phase (val i S) := by
      rw [val_insert i S A.next hyn]; exact hA.phase.update hyD _ _
    have hy1 : val i (insert A.next S) A.next = true := by
      rw [val_insert i S A.next hyn]; simp
    simp only [key0, key1, hq1, hp1, hy1, hνy, Bool.and_false, Bool.and_true]
    simp only [G]
    by_cases hr : rest S
    · simp only [hr, and_true]
      cases hoq : o q <;> cases hf : f S <;> simp [b2q, ep_add, ep_one, hf] <;> ring
    · simp only [hr, and_false, if_false]; ring
  have rhs : (embed1 n h2 q.val * A.denote) o i =
      invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset, G S := by
    rw [embed1_mul_apply, Fintype.sum_bool]
    unfold State.denote
    rw [mul_left_comm, mul_left_comm (h2 (o q) false), ← mul_add, Finset.mul_sum,
      Finset.mul_sum, ← Finset.sum_add_distrib]
    congr 1
    refine Finset.sum_congr rfl fun S _ => ?_
    simp only [conds_update]
    simp only [G]
    by_cases hr : rest S
    · have hr' : ∀ r, r ≠ q → ev (val i S) (A.ket r) = o r := hr
      simp only [hr', and_true, hr, if_true]
      cases hf : ev (val i S) (A.ket q) <;> simp [f, h2_eq, hf] <;>
        exact fun x hx hne => absurd (hr' x hx) hne
    · have hr' : ¬ ∀ r, r ≠ q → ev (val i S) (A.ket r) = o r := hr
      simp only [hr', and_false, if_false, hr]; ring
  rw [lhs, rhs]

/-! ## Circuits -/

theorem denote_step (j : ℕ) (g : CTGate n) (A : State n) (hA : A.WF) :
    (step j g A).denote = gateUnitary n g.toGate * A.denote := by
  cases g with
  | x q => exact denote_step_x j q A
  | cx c t hne => exact denote_step_cx j c t hne A
  | cz c t hne => exact denote_step_cz j c t hne A
  | h q => exact denote_step_h j q A hA
  | z q => exact denote_rot j 1 q A
  | s q => exact denote_rot j (1 / 2) q A
  | sdg q => exact denote_rot j (-1 / 2) q A
  | t q => exact denote_rot j (1 / 4) q A
  | tdg q => exact denote_rot j (-1 / 4) q A

/-! ## Well-formedness is preserved -/

theorem Supp.mono {α : Type*} {f : Val → α} {D D' : Set ℕ} (h : Supp f D) (hD : D ⊆ D') :
    Supp f D' := fun ν ν' hν => h ν ν' fun v hv => hν v (hD hv)

theorem Supp.map2 {α β γ : Type*} {f : Val → α} {g : Val → β} {D : Set ℕ} (hf : Supp f D)
    (hg : Supp g D) (op : α → β → γ) : Supp (fun ν => op (f ν) (g ν)) D :=
  fun ν ν' hν => by simp only [hf ν ν' hν, hg ν ν' hν]

theorem Supp.map {α β : Type*} {f : Val → α} {D : Set ℕ} (hf : Supp f D) (op : α → β) :
    Supp (fun ν => op (f ν)) D :=
  fun ν ν' hν => by simp only [hf ν ν' hν]

theorem supp_update_ket {K : Fin n → Poly} {D : Set ℕ} (hK : ∀ q, Supp (fun ν => ev ν (K q)) D)
    (q : Fin n) (p : Poly) (hp : Supp (fun ν => ev ν p) D) :
    ∀ r, Supp (fun ν => ev ν (Function.update K q p r)) D := by
  intro r
  by_cases hr : r = q
  · subst hr; simpa using hp
  · simpa [Function.update_of_ne hr] using hK r

theorem wf_rot (j : ℕ) (θ : ℚ) (q : Fin n) (A : State n) (hA : A.WF) : (rot j θ q A).WF where
  temps_ge := hA.temps_ge
  temps_lt := hA.temps_lt
  n_le := hA.n_le
  ket := hA.ket
  phase := hA.phase.map2 ((hA.ket q).map fun b => θ * b2q b) (· + ·)

theorem wf_step (j : ℕ) (g : CTGate n) (A : State n) (hA : A.WF) : (step j g A).WF := by
  cases g with
  | z q => exact wf_rot j _ q A hA
  | s q => exact wf_rot j _ q A hA
  | sdg q => exact wf_rot j _ q A hA
  | t q => exact wf_rot j _ q A hA
  | tdg q => exact wf_rot j _ q A hA
  | x q =>
    exact { hA with
      ket := supp_update_ket hA.ket q _ (by simpa using (hA.ket q).map (!·)) }
  | cx c t hne =>
    exact { hA with
      ket := supp_update_ket hA.ket t _ (by
        simpa using (hA.ket t).map2 (hA.ket c) (fun a b => a != b)) }
  | cz c t hne =>
    exact { hA with
      phase := hA.phase.map2 ((hA.ket c).map2 (hA.ket t) fun a b => b2q (a && b)) (· + ·) }
  | h q =>
    have hsub : A.dom ⊆ (step j (.h q) A).dom := fun v hv => by
      rcases hv with h | h
      · exact Or.inl h
      · exact Or.inr (Finset.mem_insert_of_mem h)
    have hy : A.next ∈ (step j (.h q) A).dom := Or.inr (Finset.mem_insert_self _ _)
    have hvar : Supp (fun ν => ev ν (BoolPolynomial.var A.next)) (step j (.h q) A).dom :=
      fun ν ν' hν => by simp [hν _ hy]
    have hνy : Supp (fun ν : Val => ν A.next) (step j (.h q) A).dom := fun ν ν' hν => hν _ hy
    exact {
      temps_ge := by
        intro v hv
        rcases Finset.mem_insert.mp hv with h | h
        · rw [h]; exact hA.n_le
        · exact hA.temps_ge v h
      temps_lt := by
        intro v hv
        rcases Finset.mem_insert.mp hv with h | h
        · show v < A.next + 1; omega
        · show v < A.next + 1; have := hA.temps_lt v h; omega
      n_le := by show n ≤ A.next + 1; have := hA.n_le; omega
      ket := supp_update_ket (fun r => (hA.ket r).mono hsub) q _ hvar
      phase := (hA.phase.mono hsub).map2
        ((((hA.ket q).mono hsub).map2 hνy) fun a b => b2q (a && b)) (· + ·) }

/-! ## The initial state and whole circuits -/

theorem wf_init (n : ℕ) : (init n).WF where
  temps_ge := by simp [init]
  temps_lt := by simp [init]
  n_le := le_rfl
  ket := fun q ν ν' hν => by simp [init, hν q.val (Or.inl q.isLt)]
  phase := fun _ _ _ => rfl

theorem denote_init (n : ℕ) : (init n).denote = 1 := by
  funext o i
  simp only [State.denote, init, Finset.card_empty, pow_zero, one_mul, Finset.powerset_empty,
    Finset.sum_singleton, ev_var, ep_zero]
  have hv : ∀ q : Fin n, val i ∅ q.val = i q := fun q => by simp [val, q.isLt]
  simp only [hv, Matrix.one_apply]
  by_cases h : o = i
  · subst h; simp
  · have : ¬ ∀ q, i q = o q := fun hq => h (funext fun q => (hq q).symm)
    simp [this, h]

theorem denote_runFrom (j : ℕ) (gs : List (CTGate n)) (A : State n) (hA : A.WF) :
    (runFrom j gs A).denote = unitary n (gs.map CTGate.toGate) * A.denote ∧
      (runFrom j gs A).WF := by
  induction gs generalizing j A with
  | nil => exact ⟨by simp [runFrom], hA⟩
  | cons g gs ih =>
    obtain ⟨h1, h2⟩ := ih (j + 1) (step j g A) (wf_step j g A hA)
    refine ⟨?_, h2⟩
    simp only [runFrom, List.map_cons, unitary_cons]
    rw [h1, denote_step j g A hA, Matrix.mul_assoc]

/-- **The analysis is exact.** The state after a circuit denotes the circuit's unitary. -/
theorem denote_run (gs : List (CTGate n)) :
    (run gs).denote = unitary n (gs.map CTGate.toGate) := by
  have := (denote_runFrom 0 gs (init n) (wf_init n)).1
  unfold run
  rw [this, denote_init, Matrix.mul_one]

theorem wf_run (gs : List (CTGate n)) : (run gs).WF :=
  (denote_runFrom 0 gs (init n) (wf_init n)).2

end

end PauliFold.SF
