import PauliFold.Phase

/-!
# Pauli folding is at least as precise as phase folding

**Theorem (`containment`).** Start a candidate rotation on qubit `q` and run both abstract
analyses over any Clifford+T segment. At the end, the Pauli-fold concretization is contained
in the phase-fold concretization.

The proof is an induction over the segment with the invariant `Inv`:

* *freshness*: every variable in the wires and the label is below the fresh counter, and
* *agreement*: whenever a wire set `Q` expresses the label, the Pauli state is exactly the
  signed Z-string on `Q`.

Agreement gives the inclusion at once: if some `Q` expresses the label, the Pauli
concretization is the single operator `±Z_Q`, which is in the phase concretization; if none
does, the phase concretization is everything.

The gate cases follow Lemmas 10.7 and 10.8 of the notes. Diagonal gates fix Z-strings and
change no wire. X flips a constant and the matching sign. CX changes which wire set expresses
the label exactly as it changes the Z-string. H cannot remove the label from all wire sets and
bring it back, because the fresh variable it introduces appears in no other wire and not in
the label; so any wire set expressing the label afterwards avoids the H's qubit and already
expressed it before.
-/

namespace PauliFold

open PhaseAbs

variable {n : Nat}

/-! ## Characteristic two -/

theorem zmod2_add_self (x : ZMod 2) : x + x = 0 := by
  fin_cases x <;> rfl

theorem affine_add_self (a : Affine) : a + a = 0 :=
  Prod.ext (Finsupp.ext fun v => by simp [zmod2_add_self]) (by simp [zmod2_add_self])

theorem affine_neg (a : Affine) : -a = a := by
  rw [neg_eq_iff_add_eq_zero, affine_add_self]

/-! ## Sums of wires -/

theorem sum_add_single (w : Fin n → Affine) (p : Fin n) (δ : Affine) (Q : Finset (Fin n)) :
    ∑ i ∈ Q, (w + Pi.single p δ : Fin n → Affine) i = (∑ i ∈ Q, w i) + if p ∈ Q then δ else 0 := by
  simp only [Pi.add_apply, Finset.sum_add_distrib, Finset.sum_pi_single']

/-- Toggle membership of `c`. -/
def toggle (Q : Finset (Fin n)) (c : Fin n) : Finset (Fin n) :=
  if c ∈ Q then Q.erase c else insert c Q

theorem sum_toggle (f : Fin n → Affine) (Q : Finset (Fin n)) (c : Fin n) :
    ∑ i ∈ toggle Q c, f i = (∑ i ∈ Q, f i) + f c := by
  unfold toggle
  split_ifs with h
  · rw [← Finset.add_sum_erase Q f h, add_comm (f c), add_assoc, affine_add_self, add_zero]
  · rw [Finset.sum_insert h, add_comm]

theorem mem_toggle {Q : Finset (Fin n)} {c i : Fin n} :
    i ∈ toggle Q c ↔ (i = c ↔ i ∉ Q) := by
  unfold toggle
  split_ifs with h <;> by_cases hi : i = c <;> simp_all

theorem toggle_toggle (Q : Finset (Fin n)) (c : Fin n) : toggle (toggle Q c) c = Q := by
  ext i
  simp only [mem_toggle]
  by_cases hi : i = c <;> simp [hi]

/-! ## The Pauli transformers on Z-strings -/

theorem letter_zString (Q : Finset (Fin n)) (b : Bool) (i : Fin n) :
    (zString Q b).letter i = if i ∈ Q then .Z else .I := rfl

theorem zString_letter_cases (Q : Finset (Fin n)) (b : Bool) (i : Fin n) :
    (zString Q b).letter i = .I ∨ (zString Q b).letter i = .Z := by
  simp only [letter_zString]; split_ifs <;> simp

/-- A single-qubit table that fixes the letter at `q` leaves the string unchanged. -/
theorem apply1_of_fix {f : Letter → Bool × Letter} {q : Fin n} {P : SPauli n}
    (h : f (P.letter q) = (false, P.letter q)) : SPauli.apply1 f q P = P := by
  cases P with
  | mk neg letter =>
    simp only [SPauli.apply1] at h ⊢
    rw [h]
    simp

/-- A two-qubit table that fixes the letters at `c` and `t` leaves the string unchanged. -/
theorem apply2_of_fix {f : Letter → Letter → Bool × Letter × Letter} {c t : Fin n}
    {P : SPauli n} (h : f (P.letter c) (P.letter t) = (false, P.letter c, P.letter t)) :
    SPauli.apply2 f c t P = P := by
  cases P with
  | mk neg letter =>
    simp only [SPauli.apply2] at h ⊢
    rw [h]
    simp

/-- A single-qubit table fixing `I` and `Z` fixes every Z-string. -/
theorem apply1_zString_fix {f : Letter → Bool × Letter} (hI : f .I = (false, .I))
    (hZ : f .Z = (false, .Z)) (q : Fin n) (Q : Finset (Fin n)) (b : Bool) :
    SPauli.apply1 f q (zString Q b) = zString Q b :=
  apply1_of_fix (by rcases zString_letter_cases Q b q with h | h <;> rw [h] <;> assumption)

/-- X on `q` flips the sign of a Z-string exactly when `q ∈ Q`. -/
theorem conjX_zString (q : Fin n) (Q : Finset (Fin n)) (b : Bool) :
    SPauli.apply1 Letter.conjX q (zString Q b) = zString Q (xor b (decide (q ∈ Q))) := by
  unfold SPauli.apply1 zString
  by_cases hq : q ∈ Q <;> simp only [hq, if_true, if_false, Letter.conjX, decide_true,
    decide_false, Bool.xor_false] <;> congr 1 <;> funext i <;> by_cases hi : i = q <;>
    simp_all [Function.update]

/-- CZ fixes every Z-string. -/
theorem conjCZ_zString (c t : Fin n) (Q : Finset (Fin n)) (b : Bool) :
    SPauli.apply2 Letter.conjCZ c t (zString Q b) = zString Q b :=
  apply2_of_fix (by
    rcases zString_letter_cases Q b c with hc | hc <;>
      rcases zString_letter_cases Q b t with ht | ht <;> rw [hc, ht] <;> rfl)

/-- CX moves a Z-string exactly as it moves the wire set expressing a label:
`Z_t ↦ Z_c Z_t` toggles `c` when `t` is in the set. -/
theorem conjCX_zString (c t : Fin n) (hne : c ≠ t) (R : Finset (Fin n)) (b : Bool) :
    SPauli.apply2 Letter.conjCX c t (zString R b) =
      zString (if t ∈ R then toggle R c else R) b := by
  unfold SPauli.apply2
  simp only [letter_zString]
  by_cases hc : c ∈ R <;> by_cases ht : t ∈ R <;>
    simp only [hc, ht, if_true, if_false, Letter.conjCX, Bool.xor_false] <;>
    congr 1 <;> funext i <;>
    by_cases hit : i = t <;> by_cases hic : i = c <;>
    simp_all [Function.update, letter_zString, mem_toggle, Ne.symm hne]

/-! ## The invariant -/

/-- Every variable in the wires and the label is below the fresh counter. -/
def Fresh (a : PhaseAbs n) : Prop :=
  (∀ i v, a.fresh ≤ v → (a.wires i).1 v = 0) ∧ ∀ v, a.fresh ≤ v → a.label.1 v = 0

/-- The analyses agree: every wire set expressing the label is the Pauli state's Z-string. -/
def Inv (p : PauliAbs n) (a : PhaseAbs n) : Prop :=
  Fresh a ∧ ∀ Q, a.Expresses Q → p = .one (zString Q (a.sign Q))

theorem sign_add_const (a : PhaseAbs n) (Q : Finset (Fin n)) (δ : ZMod 2) (b : Bool)
    (hb : b = decide ((a.wireSum Q).2 + a.label.2 = 1)) :
    decide ((a.wireSum Q).2 + δ + a.label.2 = 1) = xor b (decide (δ = 1)) := by
  subst hb
  generalize (a.wireSum Q).2 = s
  generalize a.label.2 = l
  fin_cases s <;> fin_cases l <;> fin_cases δ <;> decide

theorem inv_step (g : CTGate n) {p : PauliAbs n} {a : PhaseAbs n} (h : Inv p a) :
    Inv (PauliAbs.step g p) (PhaseAbs.step g a) := by
  obtain ⟨⟨hw, hl⟩, hm⟩ := h
  cases g with
  | z q =>
    refine ⟨⟨hw, hl⟩, fun Q hQ => ?_⟩
    rw [hm Q hQ]
    exact congrArg PauliAbs.one (apply1_zString_fix (f := Letter.conjZ) rfl rfl q Q _)
  | s q =>
    refine ⟨⟨hw, hl⟩, fun Q hQ => ?_⟩
    rw [hm Q hQ]
    exact congrArg PauliAbs.one (apply1_zString_fix (f := Letter.conjS) rfl rfl q Q _)
  | sdg q =>
    refine ⟨⟨hw, hl⟩, fun Q hQ => ?_⟩
    rw [hm Q hQ]
    exact congrArg PauliAbs.one (apply1_zString_fix (f := Letter.conjSdg) rfl rfl q Q _)
  | cz c t hne =>
    refine ⟨⟨hw, hl⟩, fun Q hQ => ?_⟩
    rw [hm Q hQ]
    exact congrArg PauliAbs.one (conjCZ_zString c t Q _)
  | t q | tdg q =>
    refine ⟨⟨hw, hl⟩, fun Q hQ => ?_⟩
    rw [hm Q hQ]
    rcases zString_letter_cases Q (a.sign Q) q with h | h <;>
      simp [PauliAbs.step, h, Letter.anticommutesZ] <;> rfl
  | x q =>
    refine ⟨⟨fun i v hv => ?_, hl⟩, fun Q hQ => ?_⟩
    · simp only [PhaseAbs.step, Pi.add_apply, Prod.fst_add, Finsupp.add_apply, hw i v hv,
        zero_add]
      by_cases hi : i = q
      · subst hi; simp
      · simp [hi]
    · have hsum : (PhaseAbs.step (.x q) a).wireSum Q =
          a.wireSum Q + if q ∈ Q then ((0 : Lin), (1 : ZMod 2)) else 0 :=
        sum_add_single a.wires q (0, 1) Q
      have hQ' : a.Expresses Q := by
        unfold Expresses at hQ ⊢
        rw [hsum] at hQ
        split_ifs at hQ <;> simpa [PhaseAbs.step] using hQ
      rw [hm Q hQ']
      show PauliAbs.one (SPauli.apply1 Letter.conjX q _) = _
      rw [conjX_zString]
      congr 2
      have hsnd : ((PhaseAbs.step (.x q) a).wireSum Q).2 =
          (a.wireSum Q).2 + if q ∈ Q then 1 else 0 := by
        rw [hsum]; split_ifs <;> simp
      unfold PhaseAbs.sign
      rw [hsnd, show (PhaseAbs.step (.x q) a).label = a.label from rfl,
        sign_add_const a Q _ _ rfl]
      by_cases hq : q ∈ Q <;> simp [hq]
  | h q =>
    have hsum : ∀ Q, (PhaseAbs.step (.h q) a).wireSum Q =
        a.wireSum Q + if q ∈ Q then ((Finsupp.single a.fresh 1, 0) - a.wires q) else 0 :=
      fun Q => sum_add_single a.wires q _ Q
    refine ⟨⟨fun i v hv => ?_, fun v hv => hl v (by simp [PhaseAbs.step] at hv; omega)⟩,
      fun Q hQ => ?_⟩
    · have hv' : a.fresh ≤ v := by simp [PhaseAbs.step] at hv; omega
      have hne : a.fresh ≠ v := by simp [PhaseAbs.step] at hv; omega
      simp only [PhaseAbs.step, Pi.add_apply, Prod.fst_add, Finsupp.add_apply, hw i v hv']
      by_cases hi : i = q
      · subst hi; simp [hne, hw i v hv']
      · simp [hi]
    · -- The fresh variable shows that `q ∉ Q`.
      have hqQ : q ∉ Q := by
        intro hq
        have h0 := hQ
        unfold Expresses at h0
        rw [hsum Q, if_pos hq] at h0
        have h1 := congrArg (fun f : Lin => f a.fresh) h0
        simp only [Prod.fst_add, Prod.fst_sub, Finsupp.add_apply, Finsupp.sub_apply,
          Finsupp.single_eq_same] at h1
        have hs : (a.wireSum Q).1 a.fresh = 0 := by
          simp [PhaseAbs.wireSum, Prod.fst_sum, Finsupp.finsetSum_apply, hw _ _ le_rfl]
        change _ = a.label.1 a.fresh at h1
        rw [hs, hw q _ le_rfl, hl _ le_rfl] at h1
        simp at h1
      have hsame : (PhaseAbs.step (.h q) a).wireSum Q = a.wireSum Q := by
        rw [hsum, if_neg hqQ, add_zero]
      have hQ' : a.Expresses Q := by
        unfold Expresses at hQ ⊢; rw [hsame] at hQ; simpa [PhaseAbs.step] using hQ
      rw [hm Q hQ']
      show PauliAbs.one (SPauli.apply1 Letter.conjH q _) = _
      rw [apply1_of_fix (by simp [letter_zString, hqQ, Letter.conjH])]
      congr 2
      unfold PhaseAbs.sign
      rw [hsame]
      rfl
  | cx c t hne =>
    have hsum : ∀ Q, (PhaseAbs.step (.cx c t hne) a).wireSum Q =
        a.wireSum Q + if t ∈ Q then a.wires c else 0 :=
      fun Q => sum_add_single a.wires t _ Q
    refine ⟨⟨fun i v hv => ?_, hl⟩, fun Q hQ => ?_⟩
    · simp only [PhaseAbs.step, Pi.add_apply, Prod.fst_add, Finsupp.add_apply, hw i v hv,
        zero_add]
      by_cases hi : i = t
      · subst hi; simp [hw c v hv]
      · simp [hi]
    · set R := if t ∈ Q then toggle Q c else Q with hR
      have hsR : a.wireSum R = (PhaseAbs.step (.cx c t hne) a).wireSum Q := by
        rw [hsum, hR]
        split_ifs with ht
        · exact sum_toggle a.wires Q c
        · rw [add_zero]
      have hQ' : a.Expresses R := by
        unfold Expresses at hQ ⊢; rw [hsR]; simpa [PhaseAbs.step] using hQ
      rw [hm R hQ']
      show PauliAbs.one (SPauli.apply2 Letter.conjCX c t _) = _
      rw [conjCX_zString _ _ hne]
      have htR : t ∈ R ↔ t ∈ Q := by
        rw [hR]; split_ifs
        · simp [mem_toggle, Ne.symm hne]
        · rfl
      have hback : (if t ∈ R then toggle R c else R) = Q := by
        by_cases ht : t ∈ Q
        · rw [if_pos (htR.mpr ht), hR, if_pos ht, toggle_toggle]
        · rw [if_neg (fun h => ht (htR.mp h)), hR, if_neg ht]
      rw [hback]
      congr 2
      unfold PhaseAbs.sign
      rw [hsR]
      rfl

theorem inv_run (gs : List (CTGate n)) {p : PauliAbs n} {a : PhaseAbs n} (h : Inv p a) :
    Inv (PauliAbs.run gs p) (PhaseAbs.run gs a) := by
  induction gs generalizing p a with
  | nil => exact h
  | cons g gs ih => exact ih (inv_step g h)

/-! ## The start of a segment -/

theorem sum_singles_apply (Q : Finset (Fin n)) (j : Fin n) :
    (∑ i ∈ Q, (Finsupp.single i.val (1 : ZMod 2))) j.val = if j ∈ Q then 1 else 0 := by
  rw [Finsupp.finsetSum_apply]
  simp [Finsupp.single_apply, Fin.val_inj]

theorem inv_init (q : Fin n) : Inv (.one (zString {q} false)) (PhaseAbs.init n q) := by
  have hfst : ∀ Q : Finset (Fin n), (PhaseAbs.init n q).wireSum Q =
      (∑ i ∈ Q, Finsupp.single i.val (1 : ZMod 2), 0) := by
    intro Q
    refine Prod.ext ?_ ?_ <;> simp [PhaseAbs.wireSum, PhaseAbs.init, Prod.fst_sum, Prod.snd_sum]
  refine ⟨⟨fun i v hv => ?_, fun v hv => ?_⟩, fun Q hQ => ?_⟩
  · change n ≤ v at hv
    change (Finsupp.single i.val (1 : ZMod 2)) v = 0
    rw [Finsupp.single_apply, if_neg]
    have := i.isLt
    omega
  · change n ≤ v at hv
    change (Finsupp.single q.val (1 : ZMod 2)) v = 0
    rw [Finsupp.single_apply, if_neg]
    have := q.isLt
    omega
  · have hQeq : Q = {q} := by
      ext j
      have h1 := congrArg (fun f : Lin => f j.val) hQ
      simp only [hfst, sum_singles_apply] at h1
      change (if j ∈ Q then 1 else 0) = (Finsupp.single q.val (1 : ZMod 2)) j.val at h1
      simp only [Finsupp.single_apply, Fin.val_inj] at h1
      have key : j ∈ Q ↔ q = j := by
        by_cases hj : j ∈ Q <;> by_cases hjq : q = j <;> simp_all
      rw [Finset.mem_singleton, key, eq_comm]
    subst hQeq
    have hsign : (PhaseAbs.init n q).sign {q} = false := by
      unfold PhaseAbs.sign
      rw [hfst]
      simp [PhaseAbs.init]
    rw [hsign]

/-! ## Containment -/

theorem gamma_subset_of_inv {p : PauliAbs n} {a : PhaseAbs n} (h : Inv p a) :
    p.gamma ⊆ a.gamma := by
  by_cases hex : ∃ Q, a.Expresses Q
  · obtain ⟨Q, hQ⟩ := hex
    rw [h.2 Q hQ]
    intro M hM
    simp only [PauliAbs.gamma, Set.mem_singleton_iff] at hM
    rw [PhaseAbs.gamma, if_pos ⟨Q, hQ⟩]
    exact ⟨Q, hQ, hM⟩
  · intro M _
    rw [PhaseAbs.gamma, if_neg hex]
    trivial

/-- **Containment.** For a candidate rotation on qubit `q` and any Clifford+T segment, the
Pauli-fold concretization after the segment is contained in the phase-fold one. -/
theorem containment (q : Fin n) (gs : List (CTGate n)) :
    (PauliAbs.run gs (.one (zString {q} false))).gamma ⊆
      (PhaseAbs.run gs (PhaseAbs.init n q)).gamma :=
  gamma_subset_of_inv (inv_run gs (inv_init q))

/-! ## Strictness -/

/-- The single-qubit segment `h; h` returns the axis to `Z`, so the Pauli analysis keeps
`One(Z)`; phase folding gives each H a fresh variable and loses the label. -/
def hh : List (CTGate 1) := [.h 0, .h 0]

theorem phase_hh_gamma : (PhaseAbs.run hh (PhaseAbs.init 1 0)).gamma = Set.univ := by
  rw [PhaseAbs.gamma, if_neg]
  rintro ⟨Q, hQ⟩
  -- The label is the input variable 0; the wire now carries the variable 2.
  have h1 := congrArg (fun f : Lin => f 0) hQ
  have hwire : (PhaseAbs.run hh (PhaseAbs.init 1 0)).wires 0 = (Finsupp.single 2 1, 0) := by
    simp [PhaseAbs.run, hh, PhaseAbs.step, PhaseAbs.init]
  have hsum : ((PhaseAbs.run hh (PhaseAbs.init 1 0)).wireSum Q).1 0 = 0 := by
    have hQ1 : Q = ∅ ∨ Q = {0} := by
      rcases Finset.eq_empty_or_nonempty Q with h | ⟨i, hi⟩
      · exact Or.inl h
      · right; ext j; simp [Fin.fin_one_eq_zero j, Fin.fin_one_eq_zero i ▸ hi]
    rcases hQ1 with rfl | rfl <;> simp [PhaseAbs.wireSum, hwire]
  change ((PhaseAbs.run hh (PhaseAbs.init 1 0)).wireSum Q).1 0 = _ at h1
  rw [hsum] at h1
  simp [PhaseAbs.run, hh, PhaseAbs.step, PhaseAbs.init] at h1

/-- **Strictness.** On `h; h`, the phase-fold concretization is everything while the
Pauli-fold one is a single operator, so the inclusion of `containment` is strict. -/
theorem containment_strict :
    ¬ (PhaseAbs.run hh (PhaseAbs.init 1 0)).gamma ⊆
      (PauliAbs.run hh (.one (zString {0} false))).gamma := by
  rw [phase_hh_gamma]
  intro h
  cases hP : PauliAbs.run hh (.one (zString {0} false)) with
  | top => simp [PauliAbs.run, hh, PauliAbs.step, SPauli.conjClifford] at hP
  | one P =>
    rw [hP] at h
    have h0 := h (Set.mem_univ (0 : Op 1))
    have h1 := h (Set.mem_univ (1 : Op 1))
    simp only [PauliAbs.gamma, Set.mem_singleton_iff] at h0 h1
    exact zero_ne_one (h0.trans h1.symm)

end PauliFold
