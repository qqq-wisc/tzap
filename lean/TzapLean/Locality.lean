import TzapLean.CnotMinProof
import TzapLean.Merge

/-!
# Locality: a Circuit on a Few Wires Acts as Itself ⊗ Identity

`SuperOpt` decides a rewrite by comparing the matrices of two *windows* — gate lists touching
a handful of wires — and a window's matrix must be computed on its own wires, not on the
whole register, or the check would be exponential in the circuit's width rather than the
window's.

This file supplies the theorem that makes that legitimate. For a list `S` of distinct
in-range wires,

```
pad S A out inp = if out and inp agree off S then A (restrict S out) (restrict S inp) else 0
```

is `A` acting on `S` and as the identity elsewhere, and `unitary_pad` says a gate list living
on `S` is exactly the padding of its localized self:

```lean
unitary n gs = pad S (unitary S.length (localizeGates S gs))
```

Padding respects products and scalars, so an equality of local matrices — up to global
phase, which is what `SuperOpt` checks — lifts to the full register.
-/

namespace TzapLean

variable {n : Nat}

/-! ## Restriction and extension of basis states -/

/-- A basis state's values on the window's wires. -/
def restrict (S : List Qubit) (b : Basis n) : Basis S.length := fun i => b.get (S.getD i 0)

/-- `b`, overwritten on the window's wires by `x`. -/
def extendB (S : List Qubit) (x : Basis S.length) (b : Basis n) : Basis n := fun r =>
  match h : localIdx S (r : Nat) with
  | some i => x ⟨i, localIdx_lt h⟩
  | none => b r

/-- Wires the window does not cover. -/
def OffS (S : List Qubit) (out inp : Basis n) : Prop :=
  ∀ r : Fin n, localIdx S (r : Nat) = none → out r = inp r

instance (S : List Qubit) (out inp : Basis n) : Decidable (OffS S out inp) := by
  unfold OffS; infer_instance

theorem localIdx_eq_none_iff {S : List Qubit} {q : Qubit} :
    localIdx S q = none ↔ q ∉ S := by
  rw [localIdx, List.findIdx?_eq_none_iff]
  constructor
  · intro h hmem
    have := h q hmem
    simp at this
  · intro h x hx
    rcases hxq : (x == q) with _ | _
    · rfl
    · have hxq' : x = q := by simpa using hxq
      subst hxq'
      exact absurd hx h

theorem getD_mem {S : List Qubit} {i : Nat} (h : i < S.length) : S.getD i 0 ∈ S := by
  rw [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem h]
  simpa using List.getElem_mem h

/-- The window's own wires are found at their own index. -/
theorem localIdx_getD_self {S : List Qubit} (hnd : S.Nodup) :
    ∀ {i : Nat}, i < S.length → localIdx S (S.getD i 0) = some i := by
  induction S with
  | nil => intro i hi; simp at hi
  | cons a as ih =>
      intro i hi
      cases i with
      | zero => simp [localIdx, List.findIdx?_cons]
      | succ j =>
          have hj : j < as.length := by simpa using hi
          have hne : ¬ (a == as.getD j 0) = true := by
            simp only [beq_iff_eq]
            intro hc
            exact (List.nodup_cons.1 hnd).1 (hc ▸ getD_mem hj)
          simp only [localIdx, List.findIdx?_cons, if_neg hne, List.getD_cons_succ]
          rw [show as.findIdx? (· == as.getD j 0) = some j from ih (List.nodup_cons.1 hnd).2 hj]
          rfl

theorem restrict_extendB {S : List Qubit} (hnd : S.Nodup) (hrange : ∀ q ∈ S, q < n)
    (x : Basis S.length) (b : Basis n) : restrict S (extendB S x b) = x := by
  funext i
  have hlt : (i : Nat) < S.length := i.isLt
  have hidx : localIdx S (S.getD (i : Nat) 0) = some (i : Nat) := localIdx_getD_self hnd hlt
  have hq : S.getD (i : Nat) 0 < n := hrange _ (getD_mem hlt)
  simp only [restrict, Basis.get, dif_pos hq, extendB]
  split
  · rename_i j hj
    have hji : j = (i : Nat) := by
      rw [hidx] at hj
      exact (Option.some.injEq _ _ ▸ hj).symm
    subst hji
    rfl
  · rename_i hj
    rw [hidx] at hj
    exact absurd hj (by simp)

theorem extendB_restrict {S : List Qubit} (hrange : ∀ q ∈ S, q < n) {b b' : Basis n}
    (h : OffS S b' b) : extendB S (restrict S b) b' = b := by
  funext r
  simp only [extendB]
  split
  · rename_i i hi
    have hgetD : S.getD i 0 = (r : Nat) := localIdx_getD hi
    simp only [restrict, hgetD, Basis.get, dif_pos r.isLt]
  · rename_i hi
    exact h r hi

/-- Agreeing off the window and on it is agreeing. -/
theorem basis_ext_of_agree {S : List Qubit} {out inp : Basis n} (hoff : OffS S out inp)
    (hon : restrict S out = restrict S inp) : out = inp := by
  funext r
  rcases hi : localIdx S (r : Nat) with _ | i
  · exact hoff r hi
  · have hgetD : S.getD i 0 = (r : Nat) := localIdx_getD hi
    have := congrFun hon ⟨i, localIdx_lt hi⟩
    simp only [restrict, hgetD, Basis.get, dif_pos r.isLt] at this
    simpa using this

/-! ## Padding -/

/-- `A` on the window's wires, the identity elsewhere. -/
noncomputable def pad (S : List Qubit) (A : Density S.length) : Density n :=
  fun out inp => if OffS S out inp then A (restrict S out) (restrict S inp) else 0

@[simp] theorem pad_one (S : List Qubit) : pad (n := n) S 1 = 1 := by
  funext out inp
  simp only [pad, Matrix.one_apply]
  by_cases h : OffS S out inp
  · rw [if_pos h]
    by_cases hr : restrict S out = restrict S inp
    · rw [if_pos hr, if_pos (basis_ext_of_agree h hr)]
    · rw [if_neg hr, if_neg]
      intro hc
      exact hr (by rw [hc])
  · rw [if_neg h, if_neg]
    intro hc
    exact h (by rw [hc]; intro r _; rfl)

theorem pad_smul (S : List Qubit) (c : ℂ) (A : Density S.length) :
    pad (n := n) S (c • A) = c • pad S A := by
  funext out inp
  by_cases h : OffS S out inp <;> simp [pad, h]

/-- Summing over the whole register collapses to a sum over the window, when the summand
vanishes off the window's fiber. -/
theorem sum_fiber {S : List Qubit} (hnd : S.Nodup) (hrange : ∀ q ∈ S, q < n)
    (out : Basis n) (f : Basis n → ℂ) (hzero : ∀ m : Basis n, ¬ OffS S out m → f m = 0) :
    ∑ m : Basis n, f m = ∑ x : Basis S.length, f (extendB S x out) := by
  classical
  have hsub : ∑ m : Basis n, f m =
      ∑ m ∈ Finset.univ.filter (fun m => OffS S out m), f m := by
    refine (Finset.sum_subset (Finset.filter_subset _ _) ?_).symm
    intro m _ hm
    exact hzero m (by simpa using hm)
  rw [hsub]
  refine Finset.sum_bij' (fun m _ => restrict S m) (fun x _ => extendB S x out) ?_ ?_ ?_ ?_ ?_
  · intro m hm
    exact Finset.mem_univ _
  · intro x _
    refine Finset.mem_filter.2 ⟨Finset.mem_univ _, ?_⟩
    intro r hr
    simp only [extendB]
    split
    · rename_i i hi; rw [hi] at hr; exact absurd hr (by simp)
    · rfl
  · intro m hm
    have hoff : OffS S out m := by simpa using (Finset.mem_filter.1 hm).2
    exact extendB_restrict hrange hoff
  · intro x _
    exact restrict_extendB hnd hrange x out
  · intro m hm
    have hoff : OffS S out m := by simpa using (Finset.mem_filter.1 hm).2
    rw [extendB_restrict hrange hoff]

theorem pad_mul {S : List Qubit} (hnd : S.Nodup) (hrange : ∀ q ∈ S, q < n)
    (A B : Density S.length) : pad (n := n) S A * pad S B = pad S (A * B) := by
  funext out inp
  by_cases hoff : OffS S out inp
  · have hRHS : pad (n := n) S (A * B) out inp =
        ∑ x : Basis S.length, A (restrict S out) x * B x (restrict S inp) := by
      simp only [pad, if_pos hoff, Matrix.mul_apply]
    rw [Matrix.mul_apply, hRHS,
      sum_fiber hnd hrange out (fun m => pad S A out m * pad S B m inp)
        (by intro m hm; simp only [pad, if_neg hm, zero_mul])]
    refine Finset.sum_congr rfl fun x _ => ?_
    have hoffx : OffS S out (extendB S x out) := by
      intro r hr
      simp only [extendB]
      split
      · rename_i i hi; rw [hi] at hr; exact absurd hr (by simp)
      · rfl
    have hoffx' : OffS S (extendB S x out) inp := by
      intro r hr
      rw [← hoffx r hr]
      exact hoff r hr
    simp only [pad, if_pos hoffx, if_pos hoffx', restrict_extendB hnd hrange x out]
  · simp only [pad, if_neg hoff, Matrix.mul_apply]
    refine Finset.sum_eq_zero fun m _ => ?_
    by_cases h1 : OffS S out m
    · have h2 : ¬ OffS S m inp := by
        intro hc
        exact hoff (fun r hr => (h1 r hr).trans (hc r hr))
      simp only [pad, if_neg h2, mul_zero]
    · simp only [pad, if_neg h1, zero_mul]

/-! ## Restriction meets the gate matrices -/

/-- A window wire's local index reads back the same bit. -/
theorem restrict_get {S : List Qubit} {q : Qubit} {i : Nat} (hi : localIdx S q = some i)
    (b : Basis n) : (restrict S b).get i = b.get q := by
  rw [Basis.get, dif_pos (localIdx_lt hi), restrict, localIdx_getD hi]

/-- A local index reads the wire it names. -/
theorem restrict_apply {S : List Qubit} {q : Qubit} {i : Nat} (hi : localIdx S q = some i)
    (hq : q < n) (b : Basis n) (hik : i < S.length) : restrict S b ⟨i, hik⟩ = b ⟨q, hq⟩ := by
  simp only [restrict, localIdx_getD hi, Basis.get, dif_pos hq]

/-- Writing a window wire is writing its local index. -/
theorem restrict_set {S : List Qubit} (hnd : S.Nodup) {q : Qubit} {i : Nat}
    (hi : localIdx S q = some i) (hq : q < n) (b : Basis n) (v : Bool) :
    restrict S (b.set q v) = (restrict S b).set i v := by
  funext j
  show (b.set q v).get (S.getD (j : Nat) 0)
      = if (j : Nat) = i then v else b.get (S.getD (j : Nat) 0)
  by_cases hj : (j : Nat) = i
  · rw [if_pos hj, show S.getD (j : Nat) 0 = q from by rw [hj]; exact localIdx_getD hi,
      Basis.get_set_same _ _ _ hq]
  · have hne : S.getD (j : Nat) 0 ≠ q := by
      intro hc
      have := localIdx_getD_self hnd j.isLt
      rw [hc, hi] at this
      exact hj (by simpa using this.symm)
    rw [if_neg hj, Basis.get_set_ne _ _ _ _ hne]

/-- A gate's operands, read through the window. -/
def localIdxD (S : List Qubit) (q : Qubit) : Qubit := (localIdx S q).getD 0

theorem localIdxD_eq {S : List Qubit} {q : Qubit} {i : Nat} (hi : localIdx S q = some i) :
    localIdxD S q = i := by rw [localIdxD, hi]; rfl

theorem all_restrict {S : List Qubit} (qs : List Qubit) (hqs : ∀ q ∈ qs, q ∈ S) (b : Basis n) :
    qs.all (fun q => (restrict S b).get (localIdxD S q)) = qs.all (fun q => b.get q) := by
  induction qs with
  | nil => rfl
  | cons c cs ih =>
      have hc : c ∈ S := hqs c (by simp)
      obtain ⟨i, hi⟩ : ∃ i, localIdx S c = some i := by
        rcases hl : localIdx S c with _ | i
        · exact absurd (localIdx_eq_none_iff.1 hl) (fun h => h hc)
        · exact ⟨i, rfl⟩
      simp only [List.all_cons, localIdxD_eq hi, restrict_get hi b,
        ih (fun q hq => hqs q (by simp [hq]))]

/-! ## One gate at a time -/

/-- A one-wire gate localizes. -/
theorem embed1_pad {S : List Qubit} (hnd : S.Nodup) (hrange : ∀ q ∈ S, q < n)
    (M : Bool → Bool → ℂ) {q : Qubit} {i : Nat} (hi : localIdx S q = some i) :
    embed1 n M q = pad S (embed1 S.length M i) := by
  have hq : q < n := hrange q (localIdx_getD hi ▸ getD_mem (localIdx_lt hi))
  have hik : i < S.length := localIdx_lt hi
  funext out inp
  simp only [embed1_apply_of_lt _ hq, embed1_apply_of_lt _ hik, pad]
  by_cases hall : ∀ r : Fin n, (r : Nat) ≠ q → out r = inp r
  · have hoff : OffS S out inp := by
      intro r hr
      refine hall r ?_
      intro hc
      rw [hc, hi] at hr
      exact absurd hr (by simp)
    have hlocal : ∀ j : Fin S.length, (j : Nat) ≠ i → restrict S out j = restrict S inp j := by
      intro j hj
      have hne : S.getD (j : Nat) 0 ≠ q := by
        intro hc
        have := localIdx_getD_self hnd j.isLt
        rw [hc, hi] at this
        exact hj (by simpa using this.symm)
      have hqn : S.getD (j : Nat) 0 < n := hrange _ (getD_mem j.isLt)
      simp only [restrict, Basis.get, dif_pos hqn]
      exact hall ⟨_, hqn⟩ hne
    rw [if_pos hall, if_pos hoff, if_pos hlocal]
    rw [restrict_apply hi hq out hik, restrict_apply hi hq inp hik]
  · rw [if_neg hall]
    by_cases hoff : OffS S out inp
    · rw [if_pos hoff, if_neg]
      intro hlocal
      refine hall fun r hr => ?_
      rcases hl : localIdx S (r : Nat) with _ | j
      · exact hoff r hl
      · have hjne : j ≠ i := by
          intro hc
          subst hc
          exact hr ((localIdx_getD hl).symm.trans (localIdx_getD hi))
        have := hlocal ⟨j, localIdx_lt hl⟩ hjne
        simp only [restrict, localIdx_getD hl, Basis.get, dif_pos r.isLt] at this
        exact this
    · rw [if_neg hoff]

/-- A diagonal gate whose phase depends only on window wires localizes. -/
theorem phaseMatrix_pad {S : List Qubit} {f : Basis n → ℂ} {g : Basis S.length → ℂ}
    (h : ∀ b : Basis n, f b = g (restrict S b)) :
    phaseMatrix f = pad S (phaseMatrix g) := by
  funext out inp
  simp only [phaseMatrix, pad]
  by_cases hoff : OffS S out inp
  · rw [if_pos hoff]
    by_cases heq : out = inp
    · subst heq
      rw [if_pos rfl, if_pos rfl, h]
    · rw [if_neg heq, if_neg]
      intro hc
      exact heq (basis_ext_of_agree hoff hc)
  · rw [if_neg hoff, if_neg]
    intro hc
    exact hoff (by rw [hc]; intro r _; rfl)

/-- A controlled-`X` on window wires localizes. -/
theorem permMatrix_pad {S : List Qubit} (hnd : S.Nodup) (hrange : ∀ q ∈ S, q < n)
    {ctrls : List Qubit} {tgt : Qubit} {i : Nat} (hctrls : ∀ q ∈ ctrls, q ∈ S)
    (htgt : localIdx S tgt = some i) (htgtn : tgt < n) :
    permMatrix (fun b : Basis n => b.set tgt (b.get tgt != ctrls.all (fun q => b.get q)))
      = pad S (permMatrix (fun x : Basis S.length =>
          x.set i (x.get i != ctrls.all (fun q => x.get (localIdxD S q))))) := by
  have hcommute : ∀ b : Basis n,
      restrict S (b.set tgt (b.get tgt != ctrls.all (fun q => b.get q)))
        = (restrict S b).set i ((restrict S b).get i !=
            ctrls.all (fun q => (restrict S b).get (localIdxD S q))) := by
    intro b
    rw [restrict_set hnd htgt htgtn, restrict_get htgt b, all_restrict ctrls hctrls b]
  have hoffset : ∀ b : Basis n, OffS S b (b.set tgt (b.get tgt != ctrls.all (fun q => b.get q))) := by
    intro b r hr
    have hne : (r : Nat) ≠ tgt := by
      intro hc
      rw [hc, htgt] at hr
      exact absurd hr (by simp)
    show b r = (if (r : Nat) = tgt then _ else b r)
    rw [if_neg hne]
  funext out inp
  simp only [permMatrix, pad]
  by_cases hoff : OffS S out inp
  · rw [if_pos hoff]
    by_cases heq : out = inp.set tgt (inp.get tgt != ctrls.all (fun q => inp.get q))
    · rw [if_pos heq, if_pos]
      rw [heq, hcommute inp]
    · rw [if_neg heq, if_neg]
      intro hc
      refine heq (basis_ext_of_agree (S := S) ?_ ?_)
      · intro r hr
        exact (hoff r hr).trans (hoffset inp r hr)
      · rw [hc, hcommute inp]
  · rw [if_neg hoff, if_neg]
    intro hc
    refine hoff fun r hr => ?_
    rw [hc]
    exact (hoffset inp r hr).symm

/-! ## Localizing a gate list -/

/-- Rename a gate's wires. -/
def mapQubits (f : Qubit → Qubit) : Gate → Gate
  | .x q => .x (f q)
  | .h q => .h (f q)
  | .s q => .s (f q)
  | .sdg q => .sdg (f q)
  | .z q => .z (f q)
  | .t q => .t (f q)
  | .tdg q => .tdg (f q)
  | .rz θ q => .rz θ (f q)
  | .cnot c t => .cnot (f c) (f t)
  | .cz c t => .cz (f c) (f t)
  | .ccx a b t => .ccx (f a) (f b) (f t)
  | .ccz a b t => .ccz (f a) (f b) (f t)
  | .measure q c => .measure (f q) c
  | .reset q => .reset (f q)

/-- Renaming twice by inverse maps is renaming not at all. -/
theorem mapQubits_comp {f g : Qubit → Qubit} {x : Gate} (h : ∀ q ∈ x.qubitsOf, f (g q) = q) :
    mapQubits f (mapQubits g x) = x := by
  cases x <;> simp_all [mapQubits, Gate.qubitsOf]

@[simp] theorem qubitsOf_mapQubits (f : Qubit → Qubit) (x : Gate) :
    (mapQubits f x).qubitsOf = x.qubitsOf.map f := by
  cases x <;> rfl

@[simp] theorem isUnitary_mapQubits (f : Qubit → Qubit) (x : Gate) :
    (mapQubits f x).isUnitary = x.isUnitary := by
  cases x <;> rfl

/-- Rename a gate's wires to their positions in the window. -/
def localizeGate (S : List Qubit) : Gate → Gate := mapQubits (localIdxD S)

/-- Rename a gate list's wires to their positions in the window. -/
def localizeGates (S : List Qubit) (gs : List Gate) : List Gate := gs.map (localizeGate S)

@[simp] theorem localizeGates_nil (S : List Qubit) : localizeGates S [] = [] := rfl

@[simp] theorem localizeGates_cons (S : List Qubit) (g : Gate) (gs : List Gate) :
    localizeGates S (g :: gs) = localizeGate S g :: localizeGates S gs := rfl

theorem exists_localIdx {S : List Qubit} {q : Qubit} (h : q ∈ S) :
    ∃ i, localIdx S q = some i ∧ localIdxD S q = i := by
  rcases hl : localIdx S q with _ | i
  · exact absurd (localIdx_eq_none_iff.1 hl) (fun hc => hc h)
  · exact ⟨i, rfl, localIdxD_eq hl⟩

/-- **One gate localizes.** -/
theorem gateUnitary_pad {S : List Qubit} (hnd : S.Nodup) (hrange : ∀ q ∈ S, q < n)
    {g : Gate} (hsub : ∀ q ∈ g.qubitsOf, q ∈ S) :
    gateUnitary n g = pad S (gateUnitary S.length (localizeGate S g)) := by
  have one : ∀ q ∈ g.qubitsOf, q < n := fun q hq => hrange q (hsub q hq)
  cases g with
  | x q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange x2 hi]
  | h q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange h2 hi]
  | s q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange _ hi]
  | sdg q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange _ hi]
  | z q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange _ hi]
  | t q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange _ hi]
  | tdg q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange _ hi]
  | rz θ q =>
      obtain ⟨i, hi, hd⟩ := exists_localIdx (hsub q (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hd, embed1_pad hnd hrange _ hi]
  | cnot c t =>
      obtain ⟨j, hj, hdj⟩ := exists_localIdx (hsub c (by simp [Gate.qubitsOf]))
      obtain ⟨i, hi, hdi⟩ := exists_localIdx (hsub t (by simp [Gate.qubitsOf]))
      have htn : t < n := one t (by simp [Gate.qubitsOf])
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hdi, hdj]
      have hleft : (permMatrix fun b : Basis n => b.set t (b.get t != b.get c))
          = permMatrix (fun b : Basis n => b.set t (b.get t != [c].all (fun q => b.get q))) := by
        congr 1; funext b; simp
      rw [hleft, permMatrix_pad hnd hrange (ctrls := [c]) (fun q hq => by
        rw [List.mem_singleton.1 hq]; exact hsub c (by simp [Gate.qubitsOf])) hi htn]
      congr 2
      funext x
      simp [localIdxD_eq hj]
  | ccx a b t =>
      obtain ⟨ja, hja, hda⟩ := exists_localIdx (hsub a (by simp [Gate.qubitsOf]))
      obtain ⟨jb, hjb, hdb⟩ := exists_localIdx (hsub b (by simp [Gate.qubitsOf]))
      obtain ⟨i, hi, hdi⟩ := exists_localIdx (hsub t (by simp [Gate.qubitsOf]))
      have htn : t < n := one t (by simp [Gate.qubitsOf])
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hdi, hda, hdb]
      have hleft : (permMatrix fun x : Basis n => x.set t (x.get t != (x.get a && x.get b)))
          = permMatrix (fun x : Basis n =>
              x.set t (x.get t != [a, b].all (fun q => x.get q))) := by
        congr 1; funext x; simp
      rw [hleft, permMatrix_pad hnd hrange (ctrls := [a, b]) (fun q hq => by
        rcases List.mem_cons.1 hq with rfl | hq
        · exact hsub q (by simp [Gate.qubitsOf])
        · rw [List.mem_singleton.1 hq]; exact hsub b (by simp [Gate.qubitsOf])) hi htn]
      congr 2
      funext x
      simp [localIdxD_eq hja, localIdxD_eq hjb]
  | cz c t =>
      obtain ⟨j, hj, hdj⟩ := exists_localIdx (hsub c (by simp [Gate.qubitsOf]))
      obtain ⟨i, hi, hdi⟩ := exists_localIdx (hsub t (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hdi, hdj]
      refine phaseMatrix_pad (fun b => ?_)
      rw [restrict_get hj b, restrict_get hi b]
  | ccz a b t =>
      obtain ⟨ja, hja, hda⟩ := exists_localIdx (hsub a (by simp [Gate.qubitsOf]))
      obtain ⟨jb, hjb, hdb⟩ := exists_localIdx (hsub b (by simp [Gate.qubitsOf]))
      obtain ⟨i, hi, hdi⟩ := exists_localIdx (hsub t (by simp [Gate.qubitsOf]))
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [hdi, hda, hdb]
      refine phaseMatrix_pad (fun b => ?_)
      rw [restrict_get hja b, restrict_get hjb b, restrict_get hi b]
  | measure q c =>
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [pad_one]
  | reset q =>
      simp only [localizeGate, mapQubits, gateUnitary]
      rw [pad_one]

/-- **A gate list living on `S` is the padding of its localized self.** -/
theorem unitary_pad {S : List Qubit} (hnd : S.Nodup) (hrange : ∀ q ∈ S, q < n) :
    ∀ (gs : List Gate), (∀ g ∈ gs, ∀ q ∈ g.qubitsOf, q ∈ S) →
      unitary n gs = pad S (unitary S.length (localizeGates S gs)) := by
  intro gs
  induction gs with
  | nil => intro _; rw [localizeGates_nil, unitary_nil, unitary_nil, pad_one]
  | cons g gs ih =>
      intro hsub
      rw [unitary_cons, ih (fun x hx => hsub x (by simp [hx])),
        gateUnitary_pad hnd hrange (hsub g (by simp)), pad_mul hnd hrange,
        localizeGates_cons, unitary_cons]

/-- **A local matrix identity, up to global phase, lifts to the whole register.** This is
what lets `SuperOpt` decide a rewrite on the window's own wires. -/
theorem equivalent_of_local_smul {m : Nat} {S : List Qubit} (hnd : S.Nodup)
    (hrange : ∀ q ∈ S, q < n) {gs hs : List Gate}
    (hgl : ∀ g ∈ gs, ∀ q ∈ g.qubitsOf, q ∈ S) (hhl : ∀ g ∈ hs, ∀ q ∈ g.qubitsOf, q ∈ S)
    (hgu : ∀ g ∈ gs, g.isUnitary = true) (hhu : ∀ g ∈ hs, g.isUnitary = true)
    (c : ℂ) (hc : c * star c = 1)
    (h : unitary S.length (localizeGates S hs) = c • unitary S.length (localizeGates S gs)) :
    Equivalent n m hs gs := by
  refine equivalent_of_unitary_smul hhu hgu c hc ?_
  rw [unitary_pad hnd hrange hs hhl, unitary_pad hnd hrange gs hgl, h, pad_smul]

end TzapLean
