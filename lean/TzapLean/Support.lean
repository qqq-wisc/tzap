import TzapLean.Semantics

/-!
# Local Support and Commutation

The engine behind every rewrite in an optimization pass: *operators acting on disjoint sets
of wires commute*. Rather than prove that for each of the 12·12 pairs of gate matrices, we
prove it once for an abstract notion of locality.

`SupportedOn U S` says the matrix `U` acts only on the wires in `S`:

* `offdiag` — `U` cannot change a wire outside `S`, so entries relating basis states that
  differ off `S` are zero;
* `local'` — the entries depend only on the wires in `S`, so two index pairs that agree on
  `S` (and are diagonal off `S`) get the same entry.

`SupportedOn.mul_comm` then gives `U * V = V * U` whenever the supports are disjoint, and
`gateUnitary_supportedOn` places every gate matrix at its own operand set.

Sets of wires are `Qubit → Bool` — decidable, and directly what `Gate.qubitsOf` provides
through `Gate.support`.
-/

namespace TzapLean

open Matrix

noncomputable section

/-- A set of wires. -/
abbrev Wires := Qubit → Bool

/-- The wires a gate acts on, as a `Wires` predicate — the operands `Gate.qubitsOf`
reports (see `Gate.support_eq_qubitsOf`). For `measure` this is the qubit only; the
classical bit it writes is separate. -/
def Gate.support : Gate → Wires
  | .x q | .h q | .s q | .sdg q | .z q | .t q | .tdg q | .rz _ q | .reset q => fun r => r == q
  | .measure q _ => fun r => r == q
  | .cnot c tgt | .cz c tgt => fun r => r == c || r == tgt
  | .ccx c₁ c₂ tgt | .ccz c₁ c₂ tgt => fun r => r == c₁ || r == c₂ || r == tgt

theorem Gate.support_iff (g : Gate) (q : Qubit) : g.support q = true ↔ q ∈ g.qubitsOf := by
  cases g <;> simp [Gate.support, Gate.qubitsOf, or_assoc]

/-- Two wire sets are disjoint. -/
def Wires.Disjoint (S T : Wires) : Prop := ∀ q, S q = true → T q = false

theorem Wires.Disjoint.symm {S T : Wires} (h : S.Disjoint T) : T.Disjoint S := by
  intro q hq
  by_contra hs
  simp only [Bool.not_eq_false] at hs
  exact absurd (h q hs) (by simp [hq])

/-- `U` acts only on the wires in `S`. -/
structure SupportedOn {n : Nat} (U : Density n) (S : Wires) : Prop where
  /-- Off `S` the operator cannot move a wire: such entries vanish. -/
  offdiag : ∀ out inp : Basis n, ∀ r : Fin n, S r = false → out r ≠ inp r → U out inp = 0
  /-- On `S` the entry depends only on the wires of `S`. -/
  local' : ∀ out inp out' inp' : Basis n,
      (∀ r : Fin n, S r = true → out r = out' r) →
      (∀ r : Fin n, S r = true → inp r = inp' r) →
      (∀ r : Fin n, S r = false → out r = inp r) →
      (∀ r : Fin n, S r = false → out' r = inp' r) →
      U out inp = U out' inp'

namespace SupportedOn

/-- The conjugate transpose of a local operator is local on the same wires. -/
theorem conjTranspose {n : Nat} {U : Density n} {S : Wires} (h : SupportedOn U S) :
    SupportedOn Uᴴ S where
  offdiag out inp r hr hne := by
    simp [Matrix.conjTranspose_apply, h.offdiag inp out r hr (Ne.symm hne)]
  local' out inp out' inp' h1 h2 h3 h4 := by
    simp only [Matrix.conjTranspose_apply]
    congr 1
    exact h.local' inp out inp' out' h2 h1 (fun r hr => (h3 r hr).symm) fun r hr => (h4 r hr).symm

/-- **Disjoint supports commute.** The proof pins the single intermediate basis state that
can contribute to each product: `U` fixes the wires outside `S` and `V` those outside `T`,
so with `S` and `T` disjoint the summation variable is forced, and the two forced choices
give the same pair of entries in the opposite order. -/
theorem mul_comm {n : Nat} {U V : Density n} {S T : Wires}
    (hU : SupportedOn U S) (hV : SupportedOn V T) (hST : S.Disjoint T) :
    U * V = V * U := by
  have hTS : T.Disjoint S := hST.symm
  funext out inp
  by_cases hdiff : ∃ r : Fin n, S r = false ∧ T r = false ∧ out r ≠ inp r
  · -- Both products vanish: every intermediate state disagrees with `out` or with `inp`
    -- at a wire neither operator can move.
    obtain ⟨r, hSr, hTr, hne⟩ := hdiff
    have vanish : ∀ (A B : Density n), SupportedOn A S → SupportedOn B T →
        (A * B) out inp = 0 := by
      intro A B hA hB
      simp only [Matrix.mul_apply]
      refine Finset.sum_eq_zero fun k _ => ?_
      by_cases hk : out r = k r
      · have : k r ≠ inp r := hk ▸ hne
        simp [hB.offdiag k inp r hTr this]
      · simp [hA.offdiag out k r hSr hk]
    have h1 := vanish U V hU hV
    have h2 : (V * U) out inp = 0 := by
      simp only [Matrix.mul_apply]
      refine Finset.sum_eq_zero fun k _ => ?_
      by_cases hk : out r = k r
      · have : k r ≠ inp r := hk ▸ hne
        simp [hU.offdiag k inp r hSr this]
      · simp [hV.offdiag out k r hTr hk]
    rw [h1, h2]
  · push Not at hdiff
    have hagree : ∀ r : Fin n, S r = false → T r = false → out r = inp r := by
      intro r hs ht
      by_contra hne
      exact hne (hdiff r hs ht)
    -- The forced intermediate states.
    set k₁ : Basis n := fun r => if S r then inp r else out r with hk₁
    set k₂ : Basis n := fun r => if T r then inp r else out r with hk₂
    have sum₁ : (U * V) out inp = U out k₁ * V k₁ inp := by
      simp only [Matrix.mul_apply]
      refine Finset.sum_eq_single k₁ (fun k _ hk => ?_) (by simp)
      obtain ⟨r, hr⟩ : ∃ r : Fin n, k r ≠ k₁ r := by
        by_contra hcon
        push Not at hcon
        exact hk (funext hcon)
      by_cases hS : S r = true
      · have hT : T r = false := hST r hS
        have : k r ≠ inp r := by simpa [hk₁, hS] using hr
        simp [hV.offdiag k inp r hT this]
      · simp only [Bool.not_eq_true] at hS
        have : out r ≠ k r := by
          intro hc
          exact hr (by simpa [hk₁, hS] using hc.symm)
        simp [hU.offdiag out k r hS this]
    have sum₂ : (V * U) out inp = V out k₂ * U k₂ inp := by
      simp only [Matrix.mul_apply]
      refine Finset.sum_eq_single k₂ (fun k _ hk => ?_) (by simp)
      obtain ⟨r, hr⟩ : ∃ r : Fin n, k r ≠ k₂ r := by
        by_contra hcon
        push Not at hcon
        exact hk (funext hcon)
      by_cases hT : T r = true
      · have hS : S r = false := hTS r hT
        have : k r ≠ inp r := by simpa [hk₂, hT] using hr
        simp [hU.offdiag k inp r hS this]
      · simp only [Bool.not_eq_true] at hT
        have : out r ≠ k r := by
          intro hc
          exact hr (by simpa [hk₂, hT] using hc.symm)
        simp [hV.offdiag out k r hT this]
    have hUeq : U out k₁ = U k₂ inp := by
      refine hU.local' out k₁ k₂ inp ?_ ?_ ?_ ?_
      · intro r hr; simp [hk₂, hST r hr]
      · intro r hr; simp [hk₁, hr]
      · intro r hr; simp [hk₁, hr]
      · intro r hr
        by_cases hT : T r = true
        · simp [hk₂, hT]
        · simp only [Bool.not_eq_true] at hT
          simp [hk₂, hT, hagree r hr hT]
    have hVeq : V k₁ inp = V out k₂ := by
      refine hV.local' k₁ inp out k₂ ?_ ?_ ?_ ?_
      · intro r hr; simp [hk₁, hTS r hr]
      · intro r hr; simp [hk₂, hr]
      · intro r hr
        by_cases hS : S r = true
        · simp [hk₁, hS]
        · simp only [Bool.not_eq_true] at hS
          simp [hk₁, hS, hagree r hS hr]
      · intro r hr; simp [hk₂, hr]
    rw [sum₁, sum₂, hUeq, hVeq, _root_.mul_comm]

end SupportedOn

/-! ## The three matrix builders are local -/

theorem supportedOn_permMatrix {n : Nat} (σ : Basis n → Basis n) (S : Wires)
    (hoff : ∀ (b : Basis n) (r : Fin n), S r = false → σ b r = b r)
    (hloc : ∀ (b b' : Basis n) (r : Fin n), S r = true →
      (∀ s : Fin n, S s = true → b s = b' s) → σ b r = σ b' r) :
    SupportedOn (permMatrix σ) S where
  offdiag out inp r hr hne := by
    simp only [permMatrix, ite_eq_right_iff]
    intro h
    exact absurd ((h ▸ hoff inp r hr : out r = inp r)) hne
  local' out inp out' inp' h1 h2 h3 h4 := by
    simp only [permMatrix]
    have hiff : out = σ inp ↔ out' = σ inp' := by
      constructor
      · intro h
        funext r
        by_cases hr : S r = true
        · rw [← h1 r hr, h, hloc inp inp' r hr h2]
        · simp only [Bool.not_eq_true] at hr
          rw [hoff inp' r hr]
          exact h4 r hr
      · intro h
        funext r
        by_cases hr : S r = true
        · rw [h1 r hr, h, hloc inp' inp r hr fun s hs => (h2 s hs).symm]
        · simp only [Bool.not_eq_true] at hr
          rw [hoff inp r hr]
          exact h3 r hr
    by_cases h : out = σ inp
    · rw [if_pos h, if_pos (hiff.mp h)]
    · rw [if_neg h, if_neg fun hc => h (hiff.mpr hc)]

theorem supportedOn_phaseMatrix {n : Nat} (f : Basis n → ℂ) (S : Wires)
    (hloc : ∀ b b' : Basis n, (∀ s : Fin n, S s = true → b s = b' s) → f b = f b') :
    SupportedOn (phaseMatrix f) S where
  offdiag out inp r hr hne := by
    simp only [phaseMatrix, ite_eq_right_iff]
    intro h; exact absurd (congrFun h r) hne
  local' out inp out' inp' h1 h2 h3 h4 := by
    simp only [phaseMatrix]
    by_cases h : out = inp
    · have heq : out' = inp' := by
        funext r
        by_cases hr : S r = true
        · rw [← h1 r hr, h, ← h2 r hr]
        · simp only [Bool.not_eq_true] at hr
          exact h4 r hr
      simp only [h, heq]
      exact hloc inp inp' h2
    · have hne : out' ≠ inp' := by
        intro hc
        refine h (funext fun r => ?_)
        by_cases hr : S r = true
        · rw [h1 r hr, hc, ← h2 r hr]
        · simp only [Bool.not_eq_true] at hr
          exact h3 r hr
      simp [h, hne]

/-! ## Reading a wire through a `Wires` agreement -/

theorem Basis.get_congr {n : Nat} {S : Wires} {b b' : Basis n}
    (h : ∀ s : Fin n, S s = true → b s = b' s) (q : Qubit) (hq : S q = true) :
    b.get q = b'.get q := by
  unfold Basis.get
  split
  · exact h _ hq
  · rfl

/-- The identity is local on any wire set. -/
theorem phaseMatrix_one {n : Nat} : phaseMatrix (fun _ : Basis n => (1 : ℂ)) = 1 := by
  funext out inp
  simp [phaseMatrix, Matrix.one_apply]

theorem supportedOn_one {n : Nat} (S : Wires) : SupportedOn (1 : Density n) S := by
  rw [← phaseMatrix_one]
  exact supportedOn_phaseMatrix _ S fun _ _ _ => rfl

/-- A single-wire embedding is local on that wire.

The `dite` in `embed1` is applied pointwise rather than split in the goal: `SupportedOn` is a
structure whose fields quantify over basis states, so rewriting under it once per entry keeps
the two obligations in the shape their proofs expect. -/
theorem supportedOn_embed1 {n : Nat} (M : Bool → Bool → ℂ) (q : Qubit) :
    SupportedOn (embed1 n M q) (fun r => r == q) := by
  by_cases h : q < n
  · have he : ∀ out inp : Basis n, embed1 n M q out inp =
        if ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
        then M (out ⟨q, h⟩) (inp ⟨q, h⟩) else 0 := by
      intro out inp
      simp [embed1, h]
    refine ⟨?_, ?_⟩
    · intro out inp r hr hne
      have hrq : (r : Nat) ≠ q := by simpa using hr
      rw [he]
      simp only [ite_eq_right_iff]
      intro hall
      exact absurd (hall r hrq) hne
    · intro out inp out' inp' h1 h2 h3 h4
      have hq : ((⟨q, h⟩ : Fin n) : Nat) == q := by simp
      have hout : out ⟨q, h⟩ = out' ⟨q, h⟩ := h1 ⟨q, h⟩ hq
      have hinp : inp ⟨q, h⟩ = inp' ⟨q, h⟩ := h2 ⟨q, h⟩ hq
      have hc1 : (∀ r : Fin n, (r : Nat) ≠ q → out r = inp r) := by
        intro r hr; exact h3 r (by simpa using hr)
      have hc2 : (∀ r : Fin n, (r : Nat) ≠ q → out' r = inp' r) := by
        intro r hr; exact h4 r (by simpa using hr)
      rw [he, he, if_pos hc1, if_pos hc2, hout, hinp]
  · have hone : embed1 n M q = 1 := by simp [embed1, h]
    rw [hone]
    exact supportedOn_one _

/-! ## Every gate matrix is local on its operands -/

/-- The matrix of a gate is local on the wires the gate touches. -/
theorem gateUnitary_supportedOn (n : Nat) (g : Gate) :
    SupportedOn (gateUnitary n g) g.support := by
  cases g with
  | x q => exact supportedOn_embed1 _ q
  | h q => exact supportedOn_embed1 _ q
  | s q => exact supportedOn_embed1 _ q
  | sdg q => exact supportedOn_embed1 _ q
  | z q => exact supportedOn_embed1 _ q
  | t q => exact supportedOn_embed1 _ q
  | tdg q => exact supportedOn_embed1 _ q
  | rz θ q => exact supportedOn_embed1 _ q
  | measure q c => exact supportedOn_one _
  | reset q => exact supportedOn_one _
  | cnot c tgt =>
      refine supportedOn_permMatrix _ _ (fun b r hr => ?_) (fun b b' r hr hagree => ?_)
      · have : (r : Nat) ≠ tgt := by
          intro hc; simp [Gate.support, hc] at hr
        simp [Basis.set, this]
      · have hc : b.get c = b'.get c := Basis.get_congr hagree c (by simp [Gate.support])
        have ht : b.get tgt = b'.get tgt := Basis.get_congr hagree tgt (by simp [Gate.support])
        by_cases hrt : (r : Nat) = tgt
        · simp [Basis.set, hrt, hc, ht]
        · simp [Basis.set, hrt, hagree r hr]
  | ccx c₁ c₂ tgt =>
      refine supportedOn_permMatrix _ _ (fun b r hr => ?_) (fun b b' r hr hagree => ?_)
      · have : (r : Nat) ≠ tgt := by
          intro hc; simp [Gate.support, hc] at hr
        simp [Basis.set, this]
      · have h1 : b.get c₁ = b'.get c₁ := Basis.get_congr hagree c₁ (by simp [Gate.support])
        have h2 : b.get c₂ = b'.get c₂ := Basis.get_congr hagree c₂ (by simp [Gate.support])
        have ht : b.get tgt = b'.get tgt := Basis.get_congr hagree tgt (by simp [Gate.support])
        by_cases hrt : (r : Nat) = tgt
        · simp [Basis.set, hrt, h1, h2, ht]
        · simp [Basis.set, hrt, hagree r hr]
  | cz c tgt =>
      refine supportedOn_phaseMatrix _ _ (fun b b' hagree => ?_)
      have hc : b.get c = b'.get c := Basis.get_congr hagree c (by simp [Gate.support])
      have ht : b.get tgt = b'.get tgt := Basis.get_congr hagree tgt (by simp [Gate.support])
      simp [hc, ht]
  | ccz c₁ c₂ tgt =>
      refine supportedOn_phaseMatrix _ _ (fun b b' hagree => ?_)
      have h1 : b.get c₁ = b'.get c₁ := Basis.get_congr hagree c₁ (by simp [Gate.support])
      have h2 : b.get c₂ = b'.get c₂ := Basis.get_congr hagree c₂ (by simp [Gate.support])
      have ht : b.get tgt = b'.get tgt := Basis.get_congr hagree tgt (by simp [Gate.support])
      simp [h1, h2, ht]

/-- The measurement projector on wire `q` is local on `q`. -/
theorem proj_supportedOn (n : Nat) (q : Qubit) (b : Bool) :
    SupportedOn (proj n q b) (fun r => r == q) := by
  refine supportedOn_phaseMatrix _ _ (fun s s' hagree => ?_)
  have : s.get q = s'.get q := Basis.get_congr (S := fun r => r == q) hagree q (by simp)
  simp [this]

/-- A reset Kraus operator on wire `q` is local on `q`. It is not a permutation or a
phase, so locality is checked directly from its entries. -/
theorem resetKraus_supportedOn (n : Nat) (q : Qubit) (b : Bool) :
    SupportedOn (resetKraus n q b) (fun r => r == q) where
  offdiag out inp r hr hne := by
    have hrq : (r : Nat) ≠ q := by simpa using hr
    simp only [resetKraus, ite_eq_right_iff, and_imp]
    intro _ hout
    exact absurd (by rw [hout]; simp [Basis.set, hrq]) hne
  local' out inp out' inp' h1 h2 h3 h4 := by
    have hget : inp.get q = inp'.get q := by
      unfold Basis.get; split
      · exact h2 _ (by simp)
      · rfl
    have hiff : out = inp.set q false ↔ out' = inp'.set q false := by
      constructor
      · intro h
        funext r
        by_cases hrq : (r : Nat) = q
        · have hout : out r = false := by rw [h]; simp [Basis.set, hrq]
          rw [← h1 r (by simp [hrq]), hout]
          simp [Basis.set, hrq]
        · have hr4 := h4 r (by simpa using hrq)
          simp [Basis.set, hrq, hr4]
      · intro h
        funext r
        by_cases hrq : (r : Nat) = q
        · have hout : out' r = false := by rw [h]; simp [Basis.set, hrq]
          rw [h1 r (by simp [hrq]), hout]
          simp [Basis.set, hrq]
        · have hr3 := h3 r (by simpa using hrq)
          simp [Basis.set, hrq, hr3]
    simp only [resetKraus, hget]
    by_cases hb : inp'.get q = b
    · by_cases ho : out = inp.set q false
      · rw [if_pos ⟨hb, ho⟩, if_pos ⟨hb, hiff.mp ho⟩]
      · rw [if_neg (by tauto), if_neg (by intro hc; exact ho (hiff.mpr hc.2))]
    · rw [if_neg (by tauto), if_neg (by tauto)]

end
end TzapLean
