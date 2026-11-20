import PauliFold.TableauContainment

/-!
# The affine domain as a transition relation

This is phase folding's domain in the whole-circuit form of Amy and Lunderville's affine
relation analysis. It runs once over the circuit, with no candidate, and its concretization
is a relation between input and output basis states.

**The state** (`AffAbs`). Each wire carries an affine form over path variables:

- the inputs `x_0, …, x_{n-1}`;
- one fresh variable for each H.

The transformers are the phase domain's. X adds 1, CX adds the control's wire to the
target's, H replaces a wire by a fresh variable, and diagonal gates change nothing.

**Relations** (`AffAbs.Rel`). When the wires in `S` sum to `x_T ⊕ c`, with no path
variables, every run of the circuit satisfies

  `(output parity on S) = (input parity on T) ⊕ c`.

**Concretization** (`affGamma`). The matrices whose transitions respect every relation:

  `affGamma gs = { U | ⟨x'|U|x⟩ ≠ 0 → ∀ (S, T, c) ∈ Rel, S·x' = T·x + c }`.

**The bridge to Pauli facts.** A relation `(S, T, c)` gives the Pauli fact
`Z_T ↦ (-1)^c Z_S` (`rel_pauliFact`). Conversely, a unitary with that fact has
`Z_S U = ± U Z_T`, which forces the parity constraint on its support
(`support_of_zfact`). Hence:

- `aff_sound`: the circuit's unitary is in `affGamma`.
- `pauli_refines_aff`, and with it `tableau_refines_aff`.
- Both are strict on `h; h` (`pauli_refines_aff_strict`, `tableau_refines_aff_strict`).
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## The domain -/

/-- Affine wire values and the fresh counter. -/
structure AffAbs (n : ℕ) where
  wires : Fin n → Affine
  fresh : ℕ

namespace AffAbs

/-- Each wire carries its input variable. -/
def init (n : ℕ) : AffAbs n := ⟨fun i => (Finsupp.single i.val 1, 0), n⟩

/-- The transformer of a gate: the phase domain's, without a label. -/
def step : CTGate n → AffAbs n → AffAbs n
  | .x q, a => { a with wires := a.wires + Pi.single q (0, 1) }
  | .cx ctl tgt _, a => { a with wires := a.wires + Pi.single tgt (a.wires ctl) }
  | .h q, a =>
      { wires := a.wires + Pi.single q ((Finsupp.single a.fresh 1, 0) - a.wires q)
        fresh := a.fresh + 1 }
  | .z _, a | .s _, a | .sdg _, a | .t _, a | .tdg _, a | .cz _ _ _, a => a

def run (gs : List (CTGate n)) (a : AffAbs n) : AffAbs n := gs.foldl (fun a g => step g a) a

/-- The input parity `∑_{i ∈ T} x_i` as a linear form. -/
def inputs (T : Finset (Fin n)) : Lin := ∑ i ∈ T, Finsupp.single i.val 1

/-- The derived relations: the wires in `S` sum to `x_T ⊕ c`. -/
def Rel (a : AffAbs n) : Set (Finset (Fin n) × Finset (Fin n) × ZMod 2) :=
  {r | ∑ i ∈ r.1, a.wires i = (inputs r.2.1, r.2.2)}

/-- Attach a label, as a phase-domain state. -/
def withLabel (a : AffAbs n) (ℓ : Affine) : PhaseAbs n := ⟨a.wires, a.fresh, ℓ⟩

theorem phase_step (g : CTGate n) (a : AffAbs n) (ℓ : Affine) :
    PhaseAbs.step g (a.withLabel ℓ) = (step g a).withLabel ℓ := by
  cases g <;> rfl

theorem phase_run (gs : List (CTGate n)) (a : AffAbs n) (ℓ : Affine) :
    PhaseAbs.run gs (a.withLabel ℓ) = (run gs a).withLabel ℓ := by
  induction gs generalizing a with
  | nil => rfl
  | cons g gs ih =>
    show PhaseAbs.run gs (PhaseAbs.step g (a.withLabel ℓ)) = _
    rw [phase_step, ih]; rfl

end AffAbs

/-- The parity of the bits of `x` in `S`. -/
def parity (S : Finset (Fin n)) (x : Basis n) : ZMod 2 := ∑ i ∈ S, if x i then 1 else 0

/-- **The affine concretization:** matrices whose nonzero transitions respect every derived
relation. -/
def affGamma (gs : List (CTGate n)) : Set (Density n) :=
  {U | ∀ x x' : Basis n, U x' x ≠ 0 →
    ∀ r ∈ ((AffAbs.init n).run gs).Rel, parity r.1 x' = parity r.2.1 x + r.2.2}

/-! ## Relations give Pauli facts -/

/-- A candidate whose label is the input parity on `T`. -/
def initT (n : ℕ) (T : Finset (Fin n)) : PhaseAbs n :=
  (AffAbs.init n).withLabel (AffAbs.inputs T, 0)

theorem inv_initT (T : Finset (Fin n)) : Inv (.one (zString T false)) (initT n T) := by
  have hfst : ∀ Q : Finset (Fin n), (initT n T).wireSum Q = (AffAbs.inputs Q, 0) := by
    intro Q
    refine Prod.ext ?_ ?_ <;>
      simp [PhaseAbs.wireSum, initT, AffAbs.withLabel, AffAbs.init, AffAbs.inputs,
        Prod.fst_sum, Prod.snd_sum]
  have hfresh : ∀ (Q : Finset (Fin n)) v, n ≤ v → AffAbs.inputs Q v = 0 := by
    intro Q v hv
    simp only [AffAbs.inputs, Finsupp.finsetSum_apply, Finsupp.single_apply]
    exact Finset.sum_eq_zero fun i _ => if_neg (by have := i.isLt; omega)
  refine ⟨⟨fun i v hv => ?_, fun v hv => hfresh T v hv⟩, fun Q hQ => ?_⟩
  · change n ≤ v at hv
    change (Finsupp.single i.val (1 : ZMod 2)) v = 0
    rw [Finsupp.single_apply, if_neg]
    have := i.isLt
    omega
  · have hQeq : Q = T := by
      ext j
      have h1 := congrArg (fun f : Lin => f j.val) hQ
      simp only [PhaseAbs.Expresses, hfst] at h1
      change AffAbs.inputs Q j.val = AffAbs.inputs T j.val at h1
      simp only [AffAbs.inputs, sum_singles_apply] at h1
      by_cases hQ' : j ∈ Q <;> by_cases hT' : j ∈ T <;> simp_all
    subst hQeq
    have hsign : (initT n Q).sign Q = false := by
      unfold PhaseAbs.sign
      rw [hfst]
      simp [initT, AffAbs.withLabel]
    rw [hsign]

/-- A relation `(S, T, c)` gives the Pauli fact `Z_T ↦ (-1)^c Z_S`. -/
theorem rel_pauliFact (gs : List (CTGate n)) {S T : Finset (Fin n)} {c : ZMod 2}
    (hr : (S, T, c) ∈ ((AffAbs.init n).run gs).Rel) :
    (zString T false, zString S (decide (c = 1))) ∈ pauliFacts gs := by
  have hinv := inv_run gs (inv_initT T)
  rw [initT, AffAbs.phase_run] at hinv
  have hsum : ((AffAbs.run gs (AffAbs.init n)).withLabel (AffAbs.inputs T, 0)).wireSum S =
      (AffAbs.inputs T, c) := hr
  have hexp : ((AffAbs.run gs (AffAbs.init n)).withLabel (AffAbs.inputs T, 0)).Expresses S := by
    unfold PhaseAbs.Expresses; rw [hsum]; rfl
  have hsign : ((AffAbs.run gs (AffAbs.init n)).withLabel (AffAbs.inputs T, 0)).sign S =
      decide (c = 1) := by
    unfold PhaseAbs.sign; rw [hsum]; simp [AffAbs.withLabel]
  have := hinv.2 S hexp
  rw [hsign] at this
  exact this

/-! ## Pauli facts give support constraints -/

/-- `±1` from a bit of `ZMod 2`. -/
def zsgn (a : ZMod 2) : ℂ := if a = 1 then -1 else 1

theorem zmod2_cases : ∀ a : ZMod 2, a = 0 ∨ a = 1 := by decide

theorem zmod2_one_add_one : (1 : ZMod 2) + 1 = 0 := by decide

theorem zmod2_one_ne_zero : (1 : ZMod 2) ≠ 0 := by decide

theorem zsgn_add (a b : ZMod 2) : zsgn (a + b) = zsgn a * zsgn b := by
  rcases zmod2_cases a with rfl | rfl <;> rcases zmod2_cases b with rfl | rfl <;>
    simp [zsgn, zmod2_one_add_one, zmod2_one_ne_zero]

theorem prod_zsgn (S : Finset (Fin n)) (o : Basis n) :
    ∏ r ∈ S, (if o r then (-1 : ℂ) else 1) = zsgn (parity S o) := by
  classical
  induction S using Finset.induction_on with
  | empty => simp [zsgn, parity]
  | insert a S ha ih =>
    rw [Finset.prod_insert ha, ih]
    simp only [parity]
    rw [Finset.sum_insert ha, zsgn_add]
    congr 1
    cases o a <;> simp [zsgn, zmod2_one_ne_zero]

/-- A Z-string is diagonal, with the parity signs. -/
theorem zString_apply (Q : Finset (Fin n)) (b : Bool) (o i : Basis n) :
    (zString Q b).toMatrix o i = if o = i then SPauli.sgn b * zsgn (parity Q o) else 0 := by
  unfold SPauli.toMatrix
  split_ifs with h
  · subst h
    have hfac : ∀ r, (if r ∈ Q then Letter.Z else Letter.I).mat (o r) (o r) =
        if r ∈ Q then (if o r then (-1 : ℂ) else 1) else 1 := by
      intro r
      by_cases hr : r ∈ Q <;> simp [hr, Letter.mat]
    simp only [zString]
    rw [Finset.prod_congr rfl (fun r _ => hfac r), Finset.prod_ite_mem, Finset.univ_inter,
      prod_zsgn]
  · obtain ⟨r, hr⟩ : ∃ r, o r ≠ i r := by
      by_contra hc; push Not at hc; exact h (funext hc)
    rw [Finset.prod_eq_zero (Finset.mem_univ r)]
    · simp
    · simp only [zString]; split_ifs <;> simp [Letter.mat, hr]

theorem id_toMatrix : (⟨false, fun _ => .I⟩ : SPauli n).toMatrix = 1 := by
  rw [← PP.toMatrix_ofS]
  exact PP.toMatrix_one

/-- The identity is tracked through every circuit. -/
theorem pauli_run_id (gs : List (CTGate n)) :
    PauliAbs.run gs (.one ⟨false, fun _ => .I⟩) = .one ⟨false, fun _ => .I⟩ := by
  have hstep : ∀ g : CTGate n,
      PauliAbs.step g (.one ⟨false, fun _ => .I⟩) = .one ⟨false, fun _ => .I⟩ := by
    intro g
    cases g <;>
      simp [PauliAbs.step, SPauli.conjClifford, SPauli.apply1, SPauli.apply2, Letter.conjH,
        Letter.conjS, Letter.conjSdg, Letter.conjX, Letter.conjZ, Letter.conjCX, Letter.conjCZ,
        Letter.anticommutesZ] <;>
      funext i <;> simp [Function.update_apply]
  induction gs with
  | nil => rfl
  | cons g gs ih =>
    show PauliAbs.run gs (PauliAbs.step g _) = _
    rw [hstep, ih]

/-- **A Z-fact forces the parity constraint on the support.** If `U` satisfies `I ↦ I`
(so it is unitary) and `Z_T ↦ (-1)^s Z_S`, then every nonzero entry `⟨x'|U|x⟩` has
`S·x' = T·x + s`. -/
theorem support_of_zfact {U : Density n}
    (hI : Holds U (⟨false, fun _ => .I⟩, ⟨false, fun _ => .I⟩))
    {S T : Finset (Fin n)} {s : Bool} (h : Holds U (zString T false, zString S s))
    (x x' : Basis n) (hne : U x' x ≠ 0) :
    parity S x' = parity T x + (if s then 1 else 0) := by
  have hUU : U * Uᴴ = 1 := by
    have := hI; simp only [Holds, id_toMatrix, conj, Matrix.mul_one] at this; exact this
  have hUU' : Uᴴ * U = 1 := mul_eq_one_comm.mp hUU
  -- `U Z_T = Z_S^s U`.
  have hcomm : U * (zString T false).toMatrix = (zString S s).toMatrix * U := by
    have := h; simp only [Holds, conj] at this
    rw [← this, Matrix.mul_assoc, Matrix.mul_assoc, hUU', Matrix.mul_one]
  have he := congrFun (congrFun hcomm x') x
  simp only [Matrix.mul_apply, zString_apply] at he
  rw [Finset.sum_eq_single x (fun k _ hk => by rw [if_neg hk, mul_zero]) (by simp),
    Finset.sum_eq_single x' (fun k _ hk => by rw [if_neg (Ne.symm hk), zero_mul]) (by simp),
    if_pos rfl, if_pos rfl] at he
  have key : U x' x * zsgn (parity T x) =
      U x' x * ((if s then -1 else 1) * zsgn (parity S x')) := by
    have e := he
    simp only [SPauli.sgn, Bool.false_eq_true, if_false, one_mul] at e
    rw [e]; ring
  have hsg := mul_left_cancel₀ hne key
  rcases zmod2_cases (parity T x) with ha | ha <;> rcases zmod2_cases (parity S x') with hb | hb <;>
    cases s <;> simp [ha, hb, zsgn, zmod2_one_add_one, zmod2_one_ne_zero] at hsg ⊢ <;>
    norm_num at hsg

/-! ## Containment -/

/-- **Pauli refines affine:** every relation is implied by the Pauli facts. -/
theorem pauli_refines_aff (gs : List (CTGate n)) : pauliGamma gs ⊆ affGamma gs := by
  intro U hU x x' hne ⟨S, T, c⟩ hr
  have hI := hU (⟨false, fun _ => .I⟩, ⟨false, fun _ => .I⟩) (pauli_run_id gs)
  have hf := hU _ (rel_pauliFact gs hr)
  have := support_of_zfact hI hf x x' hne
  rw [this]
  congr 1
  fin_cases c <;> rfl

/-- The affine domain is sound. -/
theorem aff_sound (gs : List (CTGate n)) : circ gs ∈ affGamma gs :=
  pauli_refines_aff gs (pauli_sound gs)

/-- **The tableau refines affine.** -/
theorem tableau_refines_aff (gs : List (CTGate n)) : Tab.gamma gs ⊆ affGamma gs :=
  (tab_gamma_subset_pauli gs).trans (pauli_refines_aff gs)

/-- On `h; h`, the only relation is the trivial one. -/
theorem hh_rel {r : Finset (Fin 1) × Finset (Fin 1) × ZMod 2}
    (hr : r ∈ ((AffAbs.init 1).run hh).Rel) : r.1 = ∅ ∧ r.2.1 = ∅ ∧ r.2.2 = 0 := by
  obtain ⟨S, T, c⟩ := r
  have hw : ((AffAbs.init 1).run hh).wires 0 = (Finsupp.single 2 1, 0) := by
    simp [AffAbs.run, hh, AffAbs.step, AffAbs.init]
  have hS : S = ∅ ∨ S = {0} := by
    rcases Finset.eq_empty_or_nonempty S with h | ⟨i, hi⟩
    · exact Or.inl h
    · right; ext j; simp [Fin.fin_one_eq_zero j, Fin.fin_one_eq_zero i ▸ hi]
  have hT : T = ∅ ∨ T = {0} := by
    rcases Finset.eq_empty_or_nonempty T with h | ⟨i, hi⟩
    · exact Or.inl h
    · right; ext j; simp [Fin.fin_one_eq_zero j, Fin.fin_one_eq_zero i ▸ hi]
  simp only [AffAbs.Rel, Set.mem_ofPred_eq] at hr
  rcases hS with rfl | rfl <;> rcases hT with rfl | rfl
  · simp [AffAbs.inputs] at hr
    have := congrArg Prod.snd hr
    exact ⟨rfl, rfl, by simpa using this.symm⟩
  · have := congrArg (fun p : Affine => p.1 0) hr
    simp [AffAbs.inputs] at this
  · have := congrArg (fun p : Affine => p.1 2) hr
    simp [hw, AffAbs.inputs] at this
  · have := congrArg (fun p : Affine => p.1 2) hr
    simp [hw, AffAbs.inputs] at this

/-- On `h; h`, H satisfies every relation but is excluded by the Pauli fact `Z ↦ Z`. -/
theorem pauli_refines_aff_strict : pauliGamma hh ⊂ affGamma hh := by
  refine ⟨pauli_refines_aff hh, fun h => ?_⟩
  have hH : circ [CTGate.h (0 : Fin 1)] ∈ affGamma hh := by
    intro x x' _ r hr
    obtain ⟨h1, h2, h3⟩ := hh_rel hr
    rw [h1, h2, h3]; simp [parity]
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

theorem tableau_refines_aff_strict : Tab.gamma hh ⊂ affGamma hh :=
  ⟨tableau_refines_aff hh, fun h =>
    pauli_refines_aff_strict.2 (h.trans (tab_gamma_subset_pauli hh))⟩

end

end PauliFold
