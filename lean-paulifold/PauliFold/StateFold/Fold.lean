import PauliFold.StateFold.Reduce

/-!
# StateFold's fold relation, and its soundness

`statefold d` folds two rotations when, after the reductions, their predicates are equal (or
complementary). We certify such a fold by a *trace*: a sequence of HH and ω reductions that
carries the two predicates along, substituting as it goes, and that never sums out a variable
either predicate depends on.

**Theorem (`fold_sound`).** Let a circuit be `pre ++ gₖ :: mid ++ gₗ :: post` with `gₖ` a T or
T† and `gₗ` any rotation. If a trace from the circuit's analysis ends with the two predicates
equal (complementary), then deleting `gₖ` and adding its angle (negated) after `gₗ` gives the
same unitary, up to global phase.

The proof moves the angle inside the phase polynomial. Deleting `gₖ` and adding the rotation
after `gₗ` changes the analysis only by the term `-θ [fₖ] ± θ [fₗ]` in the phase. Every reduction
of the trace applies equally with that term added, because the term avoids the eliminated
variable, and carries it along with the substitution. At the end the term is a constant,
because the predicates agree.
-/

namespace PauliFold.SF

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Shifted phases -/

/-- Add `δ` to the phase. -/
def State.shift (A : State n) (δ : Val → ℚ) : State n :=
  { A with phase := fun ν => A.phase ν + δ ν }

theorem State.dom_shift (A : State n) (δ : Val → ℚ) : (A.shift δ).dom = A.dom := rfl

theorem denote_shift_const (A : State n) (c : ℚ) :
    (A.shift fun _ => c).denote = ep c • A.denote := by
  funext o i
  rw [Matrix.smul_apply, State.denote_eq, State.denote_eq, smul_eq_mul]
  show invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset, (A.shift fun _ => c).term o (val i S) =
    ep c * (invSqrt2 ^ A.temps.card * ∑ S ∈ A.temps.powerset, A.term o (val i S))
  rw [Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun S _ => ?_
  unfold State.term State.shift
  split_ifs <;> simp only [ep_add, mul_zero] <;> ring

/-- The denotation depends only on the ket, the path variables, and the phase. -/
theorem denote_congr' (A B : State n) (hk : A.ket = B.ket) (ht : A.temps = B.temps)
    (hp : ∀ ν, A.phase ν = B.phase ν) : A.denote = B.denote := by
  funext o i
  simp only [State.denote_eq, State.term, hk, ht, hp]

/-! ## Well-formedness of reduction results -/

theorem Supp.subst {A : State n} {y z : ℕ} {R p : Poly}
    (hp : Supp (fun ν => ev ν p) (A.dom \ {y}))
    (hR : Supp (fun ν => ev ν R) (A.dom \ {y, z})) :
    Supp (fun ν => ev ν (subst z R p)) {v | v < n ∨ v ∈ (A.temps.erase y).erase z} := by
  intro ν ν' hν
  simp only [ev_subst]
  have hRν : ev ν R = ev ν' R := hR ν ν' fun v hv => hν v (by
    rcases hv with ⟨hv | hv, hne⟩
    · exact Or.inl hv
    · simp only [Set.mem_insert_iff, Set.mem_singleton_iff, not_or] at hne
      exact Or.inr (Finset.mem_erase.mpr ⟨hne.2, Finset.mem_erase.mpr ⟨hne.1, hv⟩⟩))
  apply hp
  intro v hv
  by_cases hvz : v = z
  · subst hvz; simp [hRν]
  · rw [Function.update_of_ne hvz, Function.update_of_ne hvz]
    rcases hv with ⟨hv | hv, hne⟩
    · exact hν v (Or.inl hv)
    · exact hν v (Or.inr (Finset.mem_erase.mpr ⟨hvz, Finset.mem_erase.mpr ⟨hne, hv⟩⟩))

theorem Supp.substFun {A : State n} {y z : ℕ} {R : Poly} {Φ : Val → ℚ}
    (hp : Supp Φ (A.dom \ {y})) (hR : Supp (fun ν => ev ν R) (A.dom \ {y, z})) :
    Supp (fun ν => Φ (Function.update ν z (ev ν R)))
      {v | v < n ∨ v ∈ (A.temps.erase y).erase z} := by
  intro ν ν' hν
  have hRν : ev ν R = ev ν' R := hR ν ν' fun v hv => hν v (by
    rcases hv with ⟨hv | hv, hne⟩
    · exact Or.inl hv
    · simp only [Set.mem_insert_iff, Set.mem_singleton_iff, not_or] at hne
      exact Or.inr (Finset.mem_erase.mpr ⟨hne.2, Finset.mem_erase.mpr ⟨hne.1, hv⟩⟩))
  apply hp
  intro v hv
  by_cases hvz : v = z
  · subst hvz; simp [hRν]
  · rw [Function.update_of_ne hvz, Function.update_of_ne hvz]
    rcases hv with ⟨hv | hv, hne⟩
    · exact hν v (Or.inl hv)
    · exact hν v (Or.inr (Finset.mem_erase.mpr ⟨hvz, Finset.mem_erase.mpr ⟨hne, hv⟩⟩))

theorem wf_hh {d : ℕ} {A B : State n} {y z : ℕ} {R : Poly} {Φ₀ : Val → ℚ}
    (h : HH d A B y z R Φ₀) (hA : A.WF) : B.WF := by
  obtain ⟨-, -, -, hket, -, hR, h₀, -, rfl⟩ := h
  exact {
    temps_ge := fun v hv => hA.temps_ge v (Finset.mem_of_mem_erase (Finset.mem_of_mem_erase hv))
    temps_lt := fun v hv => hA.temps_lt v (Finset.mem_of_mem_erase (Finset.mem_of_mem_erase hv))
    n_le := hA.n_le
    ket := fun q => (hket q).subst hR
    phase := h₀.substFun hR }

/-- A function supported off `y` is supported on the domain without `y`. -/
theorem Supp.erase {α : Type*} {A : State n} {y : ℕ} {f : Val → α} (h : Supp f (A.dom \ {y})) :
    Supp f {v | v < n ∨ v ∈ A.temps.erase y} := fun ν ν' hν => h ν ν' fun v hv => by
  rcases hv with ⟨hv | hv, hne⟩
  · exact hν v (Or.inl hv)
  · exact hν v (Or.inr (Finset.mem_erase.mpr ⟨hne, hv⟩))

theorem wf_omega {A B : State n} {y : ℕ} {R : Poly} {Φ₀ : Val → ℚ}
    (h : Omega A B y R Φ₀) (hA : A.WF) : B.WF := by
  obtain ⟨-, hket, hR, h₀, -, rfl⟩ := h
  exact {
    temps_ge := fun v hv => hA.temps_ge v (Finset.mem_of_mem_erase hv)
    temps_lt := fun v hv => hA.temps_lt v (Finset.mem_of_mem_erase hv)
    n_le := hA.n_le
    ket := fun q => (hket q).erase
    phase := (h₀.map2 (hR.map fun b => b2q b) fun a b => a + 1 / 4 - 1 / 2 * b).erase }

/-! ## Traces -/

/-- A reduction sequence at degree `d` from `A` to `C`, carrying two predicates `f`, `g` to
`f'`, `g'`. No step sums out a variable either predicate depends on. -/
inductive Trace (d : ℕ) : State n → Poly → Poly → State n → Poly → Poly → Prop
  | refl (A : State n) (f g : Poly) : Trace d A f g A f g
  | hh {A B C : State n} {f g f' g' : Poly} {y z : ℕ} {R : Poly} {Φ₀ : Val → ℚ} :
      HH d A B y z R Φ₀ →
      Supp (fun ν => ev ν f) (A.dom \ {y}) → Supp (fun ν => ev ν g) (A.dom \ {y}) →
      Trace d B (subst z R f) (subst z R g) C f' g' → Trace d A f g C f' g'
  | omega {A B C : State n} {f g f' g' : Poly} {y : ℕ} {R : Poly} {Φ₀ : Val → ℚ} :
      Omega A B y R Φ₀ →
      Supp (fun ν => ev ν f) (A.dom \ {y}) → Supp (fun ν => ev ν g) (A.dom \ {y}) →
      Trace d B f g C f' g' → Trace d A f g C f' g'

/-- The phase term `a [f] + b [g]`. -/
def term2 (a b : ℚ) (f g : Poly) (ν : Val) : ℚ := a * b2q (ev ν f) + b * b2q (ev ν g)

theorem wf_shift {A : State n} (hA : A.WF) {δ : Val → ℚ} (hδ : Supp δ A.dom) :
    (A.shift δ).WF :=
  { hA with phase := hA.phase.map2 hδ (· + ·) }

theorem Supp.term2 {D : Set ℕ} {f g : Poly} (hf : Supp (fun ν => ev ν f) D)
    (hg : Supp (fun ν => ev ν g) D) (a b : ℚ) : Supp (term2 a b f g) D :=
  (hf.map fun x => a * b2q x).map2 (hg.map fun x => b * b2q x) (· + ·)

theorem Supp.diff_sub {α : Type*} {A : State n} {y : ℕ} {f : Val → α}
    (h : Supp f (A.dom \ {y})) : Supp f A.dom :=
  h.mono Set.sdiff_subset

/-- **The shifting lemma.** A trace still applies with `a [f] + b [g]` added to the phase, and
carries the term along. -/
theorem trace_shift {d : ℕ} {A C : State n} {f g f' g' : Poly} (htr : Trace d A f g C f' g')
    (hA : A.WF) (hf : Supp (fun ν => ev ν f) A.dom) (hg : Supp (fun ν => ev ν g) A.dom)
    (a b : ℚ) :
    (A.shift (term2 a b f g)).denote = (C.shift (term2 a b f' g')).denote := by
  induction htr with
  | refl => rfl
  | @hh A B C f g f' g' y z R Φ₀ h hfy hgy _ ih =>
    have hB := wf_hh h hA
    -- The shifted HH step.
    have h' : HH d (A.shift (term2 a b f g)) (B.shift (term2 a b (subst z R f) (subst z R g)))
        y z R (fun ν => Φ₀ ν + term2 a b f g ν) := by
      obtain ⟨hy, hz, hyz, hket, hdeg, hR, h₀, hsplit, rfl⟩ := h
      refine ⟨hy, hz, hyz, hket, hdeg, hR, h₀.map2 (hfy.term2 hgy a b) (· + ·), ?_, ?_⟩
      · intro ν
        show ep (A.phase ν + term2 a b f g ν) = _
        rw [ep_add, hsplit, ← ep_add]
        congr 1; ring
      · show State.shift (hhResult A y z R Φ₀) _ = hhResult (A.shift _) y z R _
        unfold State.shift hhResult
        congr 1
        funext ν
        simp only [term2, ev_subst]
    have hB' : (B.shift (term2 a b (subst z R f) (subst z R g))).WF := by
      obtain ⟨-, -, -, -, -, hR, -, -, rfl⟩ := h
      exact wf_shift hB ((hfy.subst hR).term2 (hgy.subst hR) a b)
    rw [← denote_hh h' (wf_shift hA ((hfy.term2 hgy a b).diff_sub))]
    obtain ⟨-, -, -, -, -, hR, -, -, rfl⟩ := h
    exact ih hB (hfy.subst hR) (hgy.subst hR)
  | @omega A B C f g f' g' y R Φ₀ h hfy hgy _ ih =>
    have hB := wf_omega h hA
    have h' : Omega (A.shift (term2 a b f g)) (B.shift (term2 a b f g)) y R
        (fun ν => Φ₀ ν + term2 a b f g ν) := by
      obtain ⟨hy, hket, hR, h₀, hsplit, rfl⟩ := h
      refine ⟨hy, hket, hR, h₀.map2 (hfy.term2 hgy a b) (· + ·), ?_, ?_⟩
      · intro ν
        show ep (A.phase ν + term2 a b f g ν) = _
        rw [ep_add, hsplit, ← ep_add]
        congr 1; ring
      · show State.shift (omegaResult A y R Φ₀) _ = omegaResult (A.shift _) y R _
        unfold State.shift omegaResult
        congr 1
        funext ν
        ring
    rw [← denote_omega h' (wf_shift hA ((hfy.term2 hgy a b).diff_sub))]
    obtain ⟨-, -, -, -, -, rfl⟩ := h
    exact ih hB hfy.erase hgy.erase

/-! ## Runs that differ by a phase term -/

/-- `B` is `A` with `δ` added to the phase (sites aside). -/
structure Rel (δ : Val → ℚ) (A B : State n) : Prop where
  ket : B.ket = A.ket
  temps : B.temps = A.temps
  next : B.next = A.next
  phase : ∀ ν, B.phase ν = A.phase ν + δ ν

theorem Rel.denote {δ : Val → ℚ} {A B : State n} (h : Rel δ A B) :
    B.denote = (A.shift δ).denote :=
  denote_congr' B (A.shift δ) h.ket h.temps h.phase

theorem Rel.step {δ : Val → ℚ} {A B : State n} (h : Rel δ A B) (j j' : ℕ) (g : CTGate n) :
    Rel δ (step j g A) (step j' g B) := by
  obtain ⟨hk, ht, hn, hp⟩ := h
  cases g with
  | x q =>
    exact ⟨by show Function.update B.ket q (B.ket q + 1) = Function.update A.ket q (A.ket q + 1); rw [hk], ht, hn, hp⟩
  | cx c t hne =>
    exact ⟨by show Function.update B.ket t (B.ket t + B.ket c) = Function.update A.ket t (A.ket t + A.ket c); rw [hk], ht, hn, hp⟩
  | cz c t hne =>
    exact ⟨hk, ht, hn, fun ν => by
      show B.phase ν + b2q (ev ν (B.ket c) && ev ν (B.ket t)) =
        A.phase ν + b2q (ev ν (A.ket c) && ev ν (A.ket t)) + δ ν
      rw [hp, hk]; ring⟩
  | h q =>
    exact ⟨by show Function.update B.ket q (BoolPolynomial.var B.next) = Function.update A.ket q (BoolPolynomial.var A.next); rw [hk, hn],
      by show insert B.next B.temps = insert A.next A.temps; rw [ht, hn],
      by show B.next + 1 = A.next + 1; rw [hn],
      fun ν => by
        show B.phase ν + b2q (ev ν (B.ket q) && ν B.next) =
          A.phase ν + b2q (ev ν (A.ket q) && ν A.next) + δ ν
        rw [hp, hk, hn]; ring⟩
  | z q | s q | sdg q | t q | tdg q =>
    exact ⟨hk, ht, hn, fun ν => by
      show B.phase ν + _ * b2q (ev ν (B.ket q)) = A.phase ν + _ * b2q (ev ν (A.ket q)) + δ ν
      rw [hp, hk]; ring⟩

theorem Rel.runFrom {δ : Val → ℚ} {A B : State n} (h : Rel δ A B) (j j' : ℕ)
    (gs : List (CTGate n)) : Rel δ (runFrom j gs A) (runFrom j' gs B) := by
  induction gs generalizing A B j j' with
  | nil => exact h
  | cons g gs ih => exact ih (h.step j j' g) (j + 1) (j' + 1)

theorem Rel.trans {δ₁ δ₂ : Val → ℚ} {A B C : State n} (h₁ : Rel δ₁ A B) (h₂ : Rel δ₂ B C) :
    Rel (fun ν => δ₁ ν + δ₂ ν) A C :=
  ⟨h₂.ket.trans h₁.ket, h₂.temps.trans h₁.temps, h₂.next.trans h₁.next,
    fun ν => by rw [h₂.phase, h₁.phase]; ring⟩

theorem runFrom_append (j : ℕ) (xs ys : List (CTGate n)) (A : State n) :
    runFrom j (xs ++ ys) A = runFrom (j + xs.length) ys (runFrom j xs A) := by
  induction xs generalizing j A with
  | nil => simp [runFrom]
  | cons g xs ih =>
    simp only [List.cons_append, runFrom, List.length_cons]
    rw [ih]
    congr 1
    omega

theorem temps_step_mono (j : ℕ) (g : CTGate n) (A : State n) : A.temps ⊆ (step j g A).temps := by
  cases g <;> simp [step, rot, Finset.subset_insert]

theorem temps_runFrom_mono (j : ℕ) (gs : List (CTGate n)) (A : State n) :
    A.temps ⊆ (runFrom j gs A).temps := by
  induction gs generalizing j A with
  | nil => exact le_rfl
  | cons g gs ih => exact (temps_step_mono j g A).trans (ih (j + 1) _)

theorem dom_mono {A B : State n} (h : A.temps ⊆ B.temps) : A.dom ⊆ B.dom := fun v hv => by
  rcases hv with hv | hv
  · exact Or.inl hv
  · exact Or.inr (h hv)

/-! ## Rotation gates -/

/-- The angle and qubit of a rotation gate. -/
def rotAngle : CTGate n → Option (ℚ × Fin n)
  | .t q => some (1 / 4, q)
  | .tdg q => some (-1 / 4, q)
  | .s q => some (1 / 2, q)
  | .sdg q => some (-1 / 2, q)
  | .z q => some (1, q)
  | _ => none

theorem step_of_rotAngle {g : CTGate n} {θ : ℚ} {q : Fin n} (h : rotAngle g = some (θ, q))
    (j : ℕ) (A : State n) : step j g A = rot j θ q A := by
  cases g <;> simp only [rotAngle, Option.some.injEq, Prod.mk.injEq, reduceCtorEq] at h <;>
    obtain ⟨rfl, rfl⟩ := h <;> rfl

/-- T for `θ = 1/4` and T† otherwise. -/
def tGate (θ : ℚ) (q : Fin n) : CTGate n := if θ = 1 / 4 then .t q else .tdg q

theorem step_tGate {θ : ℚ} (hθ : θ = 1 / 4 ∨ θ = -1 / 4) (j : ℕ) (q : Fin n) (A : State n) :
    step j (tGate θ q) A = rot j θ q A := by
  rcases hθ with rfl | rfl
  · simp [tGate, step]
  · simp only [tGate, show (-1 / 4 : ℚ) ≠ 1 / 4 by norm_num, if_false]; rfl

theorem rel_rot (j : ℕ) (θ : ℚ) (q : Fin n) (A : State n) :
    Rel (fun ν => θ * b2q (ev ν (A.ket q))) A (rot j θ q A) :=
  ⟨rfl, rfl, rfl, fun _ => rfl⟩

/-! ## Soundness of folding -/

/-- The predicate of the rotation at the end of `pre`, on qubit `q`. -/
def predAfter (pre : List (CTGate n)) (q : Fin n) : Poly := (runFrom 0 pre (init n)).ket q

/-- The predicate of the rotation after `pre ++ [gₖ] ++ mid`, on qubit `q` (a rotation does not
change the ket, so `gₖ` can be dropped). -/
def predAfter₂ (pre mid : List (CTGate n)) (q : Fin n) : Poly :=
  (runFrom pre.length mid (runFrom 0 pre (init n))).ket q

theorem sign_mul_quarter {θ : ℚ} (hθ : θ = 1 / 4 ∨ θ = -1 / 4) (s : Bool) :
    (if s then -1 else 1) * θ = 1 / 4 ∨ (if s then -1 else 1) * θ = -1 / 4 := by
  rcases hθ with rfl | rfl <;> cases s <;> norm_num

/-- **Folding is sound.** If a trace from the analysis of `pre ++ gₖ :: mid ++ gₗ :: post`
carries the two rotations' predicates to equal (`s = false`) or complementary (`s = true`)
ones, then deleting the T or T† `gₖ` and adding its angle, negated when `s`, as a T or T†
right after `gₗ` gives the same unitary up to the global phase `e^{-iπθₖ}` (when `s`). -/
theorem fold_sound {d : ℕ} (pre mid post : List (CTGate n)) (gk gl : CTGate n)
    {θk θl : ℚ} {qk ql : Fin n} (hk : rotAngle gk = some (θk, qk))
    (hθ : θk = 1 / 4 ∨ θk = -1 / 4) (hl : rotAngle gl = some (θl, ql)) (s : Bool)
    {C : State n} {f' g' : Poly}
    (htr : Trace d (run (pre ++ gk :: mid ++ gl :: post)) (predAfter pre qk)
      (predAfter₂ pre mid ql) C f' g')
    (heq : ∀ ν, ev ν g' = (s != ev ν f')) :
    unitary n ((pre ++ mid ++ gl :: tGate ((if s then -1 else 1) * θk) ql :: post).map
        CTGate.toGate) =
      ep (if s then -θk else 0) •
        unitary n ((pre ++ gk :: mid ++ gl :: post).map CTGate.toGate) := by
  set σ : ℚ := if s then -1 else 1
  set P := runFrom 0 pre (init n)
  have hP : P.WF := (denote_runFrom 0 pre (init n) (wf_init n)).2
  set fk := predAfter pre qk
  set fl := predAfter₂ pre mid ql
  -- The two analyses, gate by gate.
  set O₁ := runFrom (pre.length + 1) mid (rot pre.length θk qk P)
  set W₁ := runFrom pre.length mid P
  have hO : run (pre ++ gk :: mid ++ gl :: post) =
      runFrom ((pre ++ gk :: mid).length + 1) post
        (rot (pre ++ gk :: mid).length θl ql O₁) := by
    unfold run
    rw [runFrom_append, runFrom_append 0 pre (gk :: mid)]
    simp only [runFrom, step_of_rotAngle hl, step_of_rotAngle hk, zero_add]
    rfl
  have hW : runFrom 0 (pre ++ mid ++ gl :: tGate (σ * θk) ql :: post) (init n) =
      runFrom ((pre ++ mid).length + 2) post
        (rot ((pre ++ mid).length + 1) (σ * θk) ql (rot (pre ++ mid).length θl ql W₁)) := by
    rw [runFrom_append, runFrom_append 0 pre mid]
    have hσ : σ * θk = 1 / 4 ∨ σ * θk = -1 / 4 := sign_mul_quarter hθ s
    simp only [runFrom, step_of_rotAngle hl, step_tGate hσ, zero_add]
    rfl
  -- The rewritten analysis is the original with `-θₖ [fₖ] + σ θₖ [fₗ]` added.
  have r₀ : Rel (fun ν => -θk * b2q (ev ν fk)) (rot pre.length θk qk P) P :=
    ⟨rfl, rfl, rfl, fun ν => by simp [rot, fk, predAfter, P]⟩
  have r₁ := (r₀.runFrom (pre.length + 1) pre.length mid).step (pre ++ gk :: mid).length
    (pre ++ mid).length gl
  rw [step_of_rotAngle hl, step_of_rotAngle hl] at r₁
  have hket : (rot (pre ++ mid).length θl ql W₁).ket ql = fl := rfl
  have r₂ := r₁.trans (rel_rot ((pre ++ mid).length + 1) (σ * θk) ql
    (rot (pre ++ mid).length θl ql W₁))
  rw [hket] at r₂
  have r₃ := r₂.runFrom ((pre ++ gk :: mid).length + 1) ((pre ++ mid).length + 2) post
  rw [← hO, ← hW] at r₃
  -- Denotations.
  have hrun := denote_run (pre ++ gk :: mid ++ gl :: post)
  have hrun' : (runFrom 0 (pre ++ mid ++ gl :: tGate (σ * θk) ql :: post) (init n)).denote =
      unitary n ((pre ++ mid ++ gl :: tGate (σ * θk) ql :: post).map CTGate.toGate) :=
    denote_run _
  rw [← hrun', ← hrun, r₃.denote]
  -- The predicates are supported on the final domain.
  have hfin := wf_run (pre ++ gk :: mid ++ gl :: post)
  have hPO : P.temps ⊆ O₁.temps :=
    temps_runFrom_mono (pre.length + 1) mid (rot pre.length θk qk P)
  have hOF : O₁.temps ⊆ (run (pre ++ gk :: mid ++ gl :: post)).temps := by
    rw [hO]; exact temps_runFrom_mono _ post (rot _ θl ql O₁)
  have hdomP : P.dom ⊆ (run (pre ++ gk :: mid ++ gl :: post)).dom := dom_mono (hPO.trans hOF)
  have hdomW : W₁.dom ⊆ (run (pre ++ gk :: mid ++ gl :: post)).dom := by
    apply dom_mono
    rw [show W₁.temps = O₁.temps from (r₀.runFrom (pre.length + 1) pre.length mid).temps]
    exact hOF
  have hWF₁ : W₁.WF := (denote_runFrom pre.length mid P hP).2
  have hfk : Supp (fun ν => ev ν fk) (run (pre ++ gk :: mid ++ gl :: post)).dom :=
    (hP.ket qk).mono hdomP
  have hfl : Supp (fun ν => ev ν fl) (run (pre ++ gk :: mid ++ gl :: post)).dom :=
    (hWF₁.ket ql).mono hdomW
  -- Shift through the trace.
  have hshift := trace_shift htr hfin hfk hfl (-θk) (σ * θk)
  have hzero := trace_shift htr hfin hfk hfl 0 0
  have hδ : (run (pre ++ gk :: mid ++ gl :: post)).shift
      (fun ν => -θk * b2q (ev ν fk) + σ * θk * b2q (ev ν fl)) =
      (run (pre ++ gk :: mid ++ gl :: post)).shift (term2 (-θk) (σ * θk) fk fl) := rfl
  rw [hδ, hshift]
  have hconst : (C.shift (term2 (-θk) (σ * θk) f' g')).denote =
      (C.shift fun _ => if s then -θk else 0).denote := by
    refine denote_congr' (C.shift (term2 (-θk) (σ * θk) f' g'))
      (C.shift fun _ => if s then -θk else 0) rfl rfl fun ν => ?_
    simp only [State.shift, term2, heq, σ]
    cases s <;> cases ev ν f' <;> simp [b2q] <;> ring
  rw [hconst, denote_shift_const]
  congr 1
  have h0 : (run (pre ++ gk :: mid ++ gl :: post)).denote =
      ((run (pre ++ gk :: mid ++ gl :: post)).shift (term2 0 0 fk fl)).denote :=
    denote_congr' _ _ rfl rfl fun ν => by simp [State.shift, term2]
  rw [h0, hzero]
  exact denote_congr' _ _ rfl rfl fun ν => by simp [State.shift, term2]


/-- A trace at degree `d` is a trace at every larger degree. -/
theorem Trace.mono {d d' : ℕ} (hd : d ≤ d') {A C : State n} {f g f' g' : Poly}
    (htr : Trace d A f g C f' g') : Trace d' A f g C f' g' := by
  induction htr with
  | refl A f g => exact .refl A f g
  | hh h hf hg _ ih => exact .hh { h with deg := h.deg.trans hd } hf hg ih
  | omega h hf hg _ ih => exact .omega h hf hg ih

end

end PauliFold.SF
