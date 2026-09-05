import TzapLean.Equivalence

/-!
# Applying a Set of Rewrites at Once

Rust's superoptimizer keeps every live window open as it scans, and selects rewrites greedily
over the *whole* circuit: the first window offered with a strictly shorter replacement claims
its gates, anything overlapping a claimed gate is refused afterwards, and one final pass
splices every selected replacement in at its window's anchor (`RewriteSet::apply`). A window
is a subsequence, not a run, so the gates it claims are scattered — the ones in between belong
to other windows, or to nothing, and stay where they are.

That is a different shape of theorem from "one rewrite, on a prefix of what is left", which is
what a scan that commits its first hit and moves on needs. This file is that theorem.

The trick is to stop talking about indices. A selection is a **tagging**: every gate carries
the rewrite that claims it, or nothing. `applyAll` then emits a rewrite's replacement at the
first gate it claims, drops the others, and copies untagged gates through — a plain structural
recursion, which is what makes `applyAll_correct` a list induction instead of a permutation
argument over positions.

Two conditions make it sound, and both are checked rather than assumed (see
`TzapLean/SuperOptProof.lean`):

* `OnSupp` — a claimed gate lives on its rewrite's wires, and is unitary;
* `Sep` — a gate that some *later* gate's rewrite does not claim misses that rewrite's wires.

`Sep` is what lets a rewrite's scattered gates be gathered at its first one: everything they
cross on the way is disjoint from them. Both survive taking sublists, which is exactly what
the recursion does to the list.
-/

namespace TzapLean

/-- Which rewrite claims a gate, if any. -/
abbrev Claim := Option Nat

/-- A gate together with the rewrite that claims it. -/
abbrev Tagged := Gate × Claim

/-- The gates, forgetting the tags. -/
def untag (xs : List Tagged) : List Gate := xs.map (·.1)

/-- The gates rewrite `w` claims, in order. -/
def claimedBy (w : Nat) (xs : List Tagged) : List Gate :=
  xs.filterMap fun p => if p.2 = some w then some p.1 else none

/-! ### Gathering every claim in one pass

`claimedBy` is the proof-friendly specification for one rewrite.  Calling it once per
rewrite is quadratic in a rewrite-heavy circuit, so the executable checker groups every tag
once and indexes the result. -/

/-- One past the largest rewrite id appearing in a tagging. -/
def claimBound : List Tagged → Nat
  | [] => 0
  | (_, none) :: xs => claimBound xs
  | (_, some w) :: xs => max (w + 1) (claimBound xs)

/-- Group claims in circuit order, with a fixed output bound. -/
def groupClaimsAux (count : Nat) : List Tagged → Array (List Gate)
  | [] => Array.replicate count []
  | (_, none) :: xs => groupClaimsAux count xs
  | (g, some w) :: xs =>
      let groups := groupClaimsAux count xs
      if w < count then groups.set! w (g :: groups[w]!) else groups

/-- The indexed form of `claimedBy`; unlike its specification, all ids are gathered by one
shared traversal when the surrounding `let` is evaluated. -/
def groupedClaim (xs : List Tagged) (w : Nat) : List Gate :=
  (groupClaimsAux (claimBound xs) xs)[w]?.getD []

/-- What is left of `xs` once rewrite `w`'s gates are taken out. -/
def unclaimed (w : Nat) (xs : List Tagged) : List Tagged :=
  xs.filter fun p => p.2 ≠ some w

@[simp] theorem untag_nil : untag [] = [] := rfl
@[simp] theorem untag_cons (p : Tagged) (xs : List Tagged) :
    untag (p :: xs) = p.1 :: untag xs := rfl
@[simp] theorem claimedBy_nil (w : Nat) : claimedBy w [] = [] := rfl
@[simp] theorem unclaimed_nil (w : Nat) : unclaimed w [] = [] := rfl

theorem claimedBy_cons_self (w : Nat) (g : Gate) (xs : List Tagged) :
    claimedBy w ((g, some w) :: xs) = g :: claimedBy w xs := by
  simp [claimedBy]

theorem claimedBy_cons_other {w : Nat} {t : Claim} (g : Gate) (xs : List Tagged)
    (h : t ≠ some w) : claimedBy w ((g, t) :: xs) = claimedBy w xs := by
  simp [claimedBy, h]

@[simp] theorem groupClaimsAux_size (count : Nat) : ∀ xs,
    (groupClaimsAux count xs).size = count := by
  intro xs
  induction xs with
  | nil => simp [groupClaimsAux]
  | cons p xs ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none => simpa [groupClaimsAux] using ih
      | some v =>
          rw [groupClaimsAux]
          split
          · exact (Array.size_set! _ _ _).trans ih
          · exact ih

theorem groupClaimsAux_get (count : Nat) {w : Nat} (hw : w < count) : ∀ xs,
    (groupClaimsAux count xs)[w]! = claimedBy w xs := by
  intro xs
  induction xs with
  | nil => simp [groupClaimsAux, claimedBy, hw]
  | cons p xs ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none => simpa [groupClaimsAux, claimedBy] using ih
      | some v =>
          by_cases hv : v < count
          · by_cases hvw : v = w
            · subst v
              have hwsize : w < (groupClaimsAux count xs).size := by
                simpa [groupClaimsAux_size] using hw
              rw [groupClaimsAux, if_pos hv,
                Array.getElem!_set!_self _ _ _ hwsize, ih]
              simp [claimedBy]
            · rw [groupClaimsAux, if_pos hv,
                Array.getElem!_set!_ne _ _ _ _ hvw, ih]
              simp [claimedBy, hvw]
          · have hvw : v ≠ w := by rintro rfl; exact hv hw
            rw [groupClaimsAux, if_neg hv, ih]
            simp [claimedBy, hvw]

theorem claimedBy_eq_nil_of_bound_le {w : Nat} : ∀ xs,
    claimBound xs ≤ w → claimedBy w xs = [] := by
  intro xs
  induction xs with
  | nil => simp
  | cons p xs ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none => simpa [claimBound, claimedBy] using ih
      | some v =>
          intro h
          have hrest : claimBound xs ≤ w := le_trans (Nat.le_max_right _ _) h
          have hvlt : v < w := by
            have : v + 1 ≤ w := le_trans (Nat.le_max_left _ _) h
            omega
          have hv : some v ≠ some w := by simp [Nat.ne_of_lt hvlt]
          rw [claimedBy_cons_other (w := w) g xs hv]
          exact ih hrest

/-- The batched implementation is exactly the one-rewrite specification. -/
theorem groupedClaim_eq_claimedBy (xs : List Tagged) (w : Nat) :
    groupedClaim xs w = claimedBy w xs := by
  by_cases hw : w < claimBound xs
  · have hwsize : w < (groupClaimsAux (claimBound xs) xs).size := by
      simpa [groupClaimsAux_size] using hw
    rw [groupedClaim, Array.getElem?_eq_getElem hwsize]
    simp only [Option.getD_some]
    have hget := groupClaimsAux_get (claimBound xs) hw xs
    rw [getElem!_pos (groupClaimsAux (claimBound xs) xs) w hwsize] at hget
    exact hget
  · have hsize : (groupClaimsAux (claimBound xs) xs).size = claimBound xs :=
      groupClaimsAux_size _ _
    have hout : w ≥ (groupClaimsAux (claimBound xs) xs).size := by omega
    rw [groupedClaim, Array.getElem?_eq_none hout]
    simp [claimedBy_eq_nil_of_bound_le (xs := xs) (w := w) (by omega)]

theorem unclaimed_cons_self (w : Nat) (g : Gate) (xs : List Tagged) :
    unclaimed w ((g, some w) :: xs) = unclaimed w xs := by
  simp [unclaimed]

theorem unclaimed_cons_other {w : Nat} {t : Claim} (g : Gate) (xs : List Tagged)
    (h : t ≠ some w) : unclaimed w ((g, t) :: xs) = (g, t) :: unclaimed w xs := by
  simp [unclaimed, h]

theorem unclaimed_length_le (w : Nat) (xs : List Tagged) : (unclaimed w xs).length ≤ xs.length :=
  List.length_filter_le _ _

theorem mem_unclaimed {w : Nat} {xs : List Tagged} {p : Tagged} (h : p ∈ unclaimed w xs) :
    p ∈ xs := List.mem_of_mem_filter h

theorem mem_claimedBy {w : Nat} {xs : List Tagged} {g : Gate} (h : g ∈ claimedBy w xs) :
    (g, some w) ∈ xs := by
  rcases List.mem_filterMap.1 h with ⟨⟨y, t⟩, hy, hyg⟩
  by_cases ht : t = some w
  · subst ht
    simp only [if_true, Option.some.injEq] at hyg
    exact hyg ▸ hy
  · simp [ht] at hyg

/-! ## The two side conditions -/

/-- A claimed gate lives on its rewrite's wires, and is unitary. -/
def OnSupp (supp : Nat → List Qubit) (xs : List Tagged) : Prop :=
  ∀ p ∈ xs, ∀ w, p.2 = some w → (∀ q ∈ p.1.qubitsOf, q ∈ supp w) ∧ p.1.isUnitary = true

/-- **The separation condition, for one rewrite.** A gate that rewrite `w` does not claim, but
which one of `w`'s gates follows, misses `w`'s wires.

Note where this is used: always on the *tail after a `w`-claim*, so what it says is "everything
between two of `w`'s gates". A gate before `w`'s first claim is never constrained — it is never
crossed, because `w`'s gates only ever move left as far as the first of them. -/
def SepOne (S : List Qubit) (w : Nat) : List Tagged → Prop
  | [] => True
  | p :: rest =>
      (p.2 ≠ some w → (∃ r ∈ rest, r.2 = some w) → ∀ q ∈ p.1.qubitsOf, q ∉ S) ∧
      SepOne S w rest

/-- **The separation condition.** After each claimed gate, `SepOne` for the rewrite that
claimed it. -/
def Sep (supp : Nat → List Qubit) : List Tagged → Prop
  | [] => True
  | p :: rest => (∀ w, p.2 = some w → SepOne (supp w) w rest) ∧ Sep supp rest

theorem Sep.tail {supp : Nat → List Qubit} {p : Tagged} {rest : List Tagged}
    (h : Sep supp (p :: rest)) : Sep supp rest := h.2

theorem OnSupp.tail {supp : Nat → List Qubit} {p : Tagged} {rest : List Tagged}
    (h : OnSupp supp (p :: rest)) : OnSupp supp rest :=
  fun q hq => h q (List.mem_cons_of_mem _ hq)

theorem OnSupp.unclaimed {supp : Nat → List Qubit} {xs : List Tagged} (h : OnSupp supp xs)
    (w : Nat) : OnSupp supp (TzapLean.unclaimed w xs) :=
  fun p hp => h p (mem_unclaimed hp)

/-- Both conditions only ever say less about a shorter list, which is what lets the recursion
drop gates and keep them. -/
theorem SepOne.filter {S : List Qubit} {w : Nat} (q : Tagged → Bool) :
    ∀ {xs : List Tagged}, SepOne S w xs → SepOne S w (xs.filter q)
  | [], _ => trivial
  | p :: rest, h => by
      by_cases hq : q p
      · rw [List.filter_cons_of_pos hq]
        refine ⟨fun hne hex => h.1 hne ?_, SepOne.filter q h.2⟩
        obtain ⟨r, hr, hrw⟩ := hex
        exact ⟨r, List.mem_of_mem_filter hr, hrw⟩
      · rw [List.filter_cons_of_neg (by simpa using hq)]
        exact SepOne.filter q h.2

theorem Sep.filter {supp : Nat → List Qubit} (q : Tagged → Bool) :
    ∀ {xs : List Tagged}, Sep supp xs → Sep supp (xs.filter q)
  | [], _ => trivial
  | p :: rest, h => by
      by_cases hq : q p
      · rw [List.filter_cons_of_pos hq]
        exact ⟨fun w hw => SepOne.filter q (h.1 w hw), Sep.filter q h.2⟩
      · rw [List.filter_cons_of_neg (by simpa using hq)]
        exact Sep.filter q h.2

theorem Sep.unclaimed {supp : Nat → List Qubit} {xs : List Tagged} (h : Sep supp xs) (w : Nat) :
    Sep supp (TzapLean.unclaimed w xs) := Sep.filter _ h

/-- Two gates on disjoint operand lists have disjoint supports. -/
theorem disjoint_of_notMem {g g' : Gate} (h : ∀ q ∈ g.qubitsOf, q ∉ g'.qubitsOf) :
    Wires.Disjoint g.support g'.support := by
  intro q hq
  rcases hq' : g'.support q with _ | _
  · rfl
  · exact absurd ((Gate.support_iff g' q).1 (by simp [hq'])) (h q ((Gate.support_iff g q).1 hq))

/-! ## Gathering one rewrite's gates -/

/-- **A rewrite's gates may be gathered at the first of them.** Everything they cross on the
way is a gate the rewrite does not claim but which a claimed gate follows — and `Sep` says
those miss the rewrite's wires. -/
theorem gather_equiv {n m : Nat} (supp : Nat → List Qubit) (w : Nat) :
    ∀ xs : List Tagged, SepOne (supp w) w xs → OnSupp supp xs →
      Equivalent n m (claimedBy w xs ++ untag (unclaimed w xs)) (untag xs) := by
  intro xs
  induction xs with
  | nil => intro _ _; exact Equivalent.refl _ _ _
  | cons p rest ih =>
      intro hsep hon
      obtain ⟨g, t⟩ := p
      have ihr := ih hsep.2 hon.tail
      by_cases ht : t = some w
      · subst ht
        rw [claimedBy_cons_self, unclaimed_cons_self, untag_cons, List.cons_append]
        exact Equivalent.append_left [g] ihr
      · rw [claimedBy_cons_other _ _ ht, unclaimed_cons_other _ _ ht, untag_cons, untag_cons]
        have hcomm : Equivalent n m (claimedBy w rest ++ [g]) ([g] ++ claimedBy w rest) := by
          have hmem : ∀ x ∈ claimedBy w rest, (x, some w) ∈ rest :=
            fun x hx => mem_claimedBy hx
          have honm : ∀ x ∈ claimedBy w rest,
              (∀ q ∈ x.qubitsOf, q ∈ supp w) ∧ x.isUnitary = true := fun x hx =>
            hon.tail _ (hmem x hx) w rfl
          refine Equivalent.append_comm _ _ (fun x hx => (honm x hx).2) (fun x hx y hy => ?_)
          rw [List.mem_singleton.1 hy]
          refine disjoint_of_notMem ?_
          intro q hq hq'
          exact hsep.1 ht ⟨(x, some w), hmem x hx, rfl⟩ q hq' ((honm x hx).1 q hq)
        have step₁ : Equivalent n m (claimedBy w rest ++ g :: untag (unclaimed w rest))
            (g :: (claimedBy w rest ++ untag (unclaimed w rest))) := by
          have h := Equivalent.append_right (n := n) (m := m) (untag (unclaimed w rest)) hcomm
          simpa using h
        exact step₁.trans (Equivalent.append_left [g] ihr)

/-! ## Applying every rewrite at once -/

@[simp] theorem claimedBy_unclaimed_self (w : Nat) (xs : List Tagged) :
    claimedBy w (unclaimed w xs) = [] := by
  induction xs with
  | nil => rfl
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      by_cases ht : t = some w
      · rw [ht, unclaimed_cons_self]; exact ih
      · rw [unclaimed_cons_other _ _ ht, claimedBy_cons_other _ _ ht]; exact ih

theorem claimedBy_unclaimed_other {v w : Nat} (h : v ≠ w) (xs : List Tagged) :
    claimedBy v (unclaimed w xs) = claimedBy v xs := by
  induction xs with
  | nil => rfl
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      by_cases ht : t = some w
      · rw [ht, unclaimed_cons_self, ih,
          claimedBy_cons_other _ _ (by simp [Ne.symm h])]
      · rw [unclaimed_cons_other _ _ ht]
        by_cases htv : t = some v
        · rw [htv, claimedBy_cons_self, claimedBy_cons_self, ih]
        · rw [claimedBy_cons_other _ _ htv, claimedBy_cons_other _ _ htv, ih]

/-- Splice every rewrite in: its replacement is emitted at the first gate it claims, the
gates it claims elsewhere are dropped, and everything else is copied through. This is Rust's
`RewriteSet::apply`, with the tagging standing in for its `claimed`/`anchored` arrays. -/
def applyAll (repl : Nat → List Gate) : List Tagged → List Gate
  | [] => []
  | (g, none) :: rest => g :: applyAll repl rest
  | (_, some w) :: rest => repl w ++ applyAll repl (unclaimed w rest)
termination_by xs => xs.length
decreasing_by
  · simp
  · exact Nat.lt_succ_of_le (unclaimed_length_le w rest)

/-! The definition above mirrors its correctness proof, but filtering the whole suffix after
every first claim is quadratic in rewrite-heavy circuits.  The executable form remembers
which replacements were emitted and traverses the tagging once. -/

def dropStarted (started : List Nat) : List Tagged → List Tagged
  | [] => []
  | (g, none) :: rest => (g, none) :: dropStarted started rest
  | (g, some w) :: rest =>
      if started.contains w then dropStarted started rest
      else (g, some w) :: dropStarted started rest

def applyAllLinear (repl : Nat → List Gate) : List Nat → List Tagged → List Gate
  | _, [] => []
  | started, (g, none) :: rest => g :: applyAllLinear repl started rest
  | started, (_, some w) :: rest =>
      if started.contains w then applyAllLinear repl started rest
      else repl w ++ applyAllLinear repl (w :: started) rest

theorem dropStarted_cons (w : Nat) (started : List Nat) : ∀ xs,
    dropStarted (w :: started) xs = unclaimed w (dropStarted started xs) := by
  intro xs
  induction xs with
  | nil => simp [dropStarted, unclaimed]
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none => simp [dropStarted, unclaimed, ih]
      | some v =>
          by_cases hvw : v = w
          · subst v
            by_cases hws : w ∈ started <;>
              simp [dropStarted, unclaimed, ih, hws, List.contains_eq_mem]
          · by_cases hvs : v ∈ started <;>
              simp [dropStarted, unclaimed, ih, hvw, hvs, List.contains_eq_mem]

theorem applyAllLinear_eq (repl : Nat → List Gate) : ∀ started xs,
    applyAllLinear repl started xs = applyAll repl (dropStarted started xs) := by
  intro started xs
  induction xs generalizing started with
  | nil => simp [applyAllLinear, applyAll, dropStarted]
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none => simp [applyAllLinear, applyAll, dropStarted, ih]
      | some w =>
          by_cases hw : w ∈ started
          · simp [applyAllLinear, dropStarted, hw, ih, List.contains_eq_mem]
          · simp [applyAllLinear, applyAll, dropStarted, hw, ih, dropStarted_cons,
              List.contains_eq_mem]

@[simp] theorem dropStarted_nil (xs : List Tagged) : dropStarted [] xs = xs := by
  induction xs with
  | nil => simp [dropStarted]
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      cases t <;> simp [dropStarted, ih, List.contains_eq_mem]

theorem applyAllLinear_initial (repl : Nat → List Gate) (xs : List Tagged) :
    applyAllLinear repl [] xs = applyAll repl xs := by
  rw [applyAllLinear_eq, dropStarted_nil]

/-- Every gate of the result is either one of the input's or one supplied by a replacement
that actually claims something. -/
theorem mem_applyAll (repl : Nat → List Gate) :
    ∀ (xs : List Tagged) (g : Gate), g ∈ applyAll repl xs →
      g ∈ untag xs ∨ ∃ w, (∃ p ∈ xs, p.2 = some w) ∧ g ∈ repl w := by
  intro xs
  induction xs using applyAll.induct with
  | case1 => intro g hg; rw [applyAll] at hg; exact absurd hg (by simp)
  | case2 x rest ih =>
      intro g hg
      rw [applyAll] at hg
      rcases List.mem_cons.1 hg with rfl | hg
      · exact Or.inl (by simp)
      · rcases ih g hg with h | ⟨w, hw, h⟩
        · exact Or.inl (by simp only [untag_cons, List.mem_cons]; exact Or.inr h)
        · exact Or.inr ⟨w, ⟨hw.choose, List.mem_cons_of_mem _ hw.choose_spec.1,
            hw.choose_spec.2⟩, h⟩
  | case3 x w rest ih =>
      intro g hg
      rw [applyAll] at hg
      rcases List.mem_append.1 hg with h | h
      · exact Or.inr ⟨w, ⟨_, List.mem_cons_self, rfl⟩, h⟩
      · rcases ih g h with h' | ⟨v, hv, h'⟩
        · refine Or.inl ?_
          simp only [untag, List.mem_map] at h' ⊢
          obtain ⟨q, hq, rfl⟩ := h'
          exact ⟨q, List.mem_cons_of_mem _ (mem_unclaimed hq), rfl⟩
        · exact Or.inr ⟨v, ⟨hv.choose, List.mem_cons_of_mem _ (mem_unclaimed hv.choose_spec.1),
            hv.choose_spec.2⟩, h'⟩

/-- **A set of rewrites, spliced in at once, preserves the circuit.**

Under the two conditions the checker establishes — every claimed gate on its rewrite's wires
and unitary, and every gate a later rewrite does not claim missing that rewrite's wires — the
whole splice is meaning-preserving, however the claimed gates interleave. -/
theorem applyAll_correct {n m : Nat} (supp : Nat → List Qubit) (repl : Nat → List Gate) :
    ∀ xs : List Tagged, Sep supp xs → OnSupp supp xs →
      (∀ w, (∃ p ∈ xs, p.2 = some w) → Equivalent n m (repl w) (claimedBy w xs)) →
      Equivalent n m (applyAll repl xs) (untag xs) := by
  intro xs
  induction xs using applyAll.induct with
  | case1 => intro _ _ _; rw [applyAll]; exact Equivalent.refl _ _ _
  | case2 g rest ih =>
      intro hsep hon hrepl
      rw [applyAll, untag_cons]
      refine Equivalent.append_left [g] (ih hsep.2 hon.tail fun w hw => ?_)
      have hx := hrepl w ⟨_, List.mem_cons_of_mem _ hw.choose_spec.1, hw.choose_spec.2⟩
      rwa [claimedBy_cons_other _ _ (by simp)] at hx
  | case3 g w rest ih =>
      intro hsep hon hrepl
      have hgw : Equivalent n m (repl w) (g :: claimedBy w rest) := by
        have hx := hrepl w ⟨(g, some w), by simp, rfl⟩
        rwa [claimedBy_cons_self] at hx
      have ihr := ih ((hsep.2).unclaimed w) (hon.tail.unclaimed w) (fun v hv => ?_)
      · rw [applyAll, untag_cons]
        have step₁ : Equivalent n m (repl w ++ applyAll repl (unclaimed w rest))
            ((g :: claimedBy w rest) ++ applyAll repl (unclaimed w rest)) :=
          Equivalent.append_right _ hgw
        have step₂ : Equivalent n m ((g :: claimedBy w rest) ++ applyAll repl (unclaimed w rest))
            ((g :: claimedBy w rest) ++ untag (unclaimed w rest)) :=
          Equivalent.append_left _ ihr
        have step₃ : Equivalent n m (g :: (claimedBy w rest ++ untag (unclaimed w rest)))
            (g :: untag rest) :=
          Equivalent.append_left [g] (gather_equiv supp w rest (hsep.1 w rfl) hon.tail)
        exact step₁.trans (step₂.trans (by simpa using step₃))
      · -- a rewrite still present in what is left is not `w`, and claims the same gates
        obtain ⟨q, hq, hqv⟩ := hv
        have hvw : v ≠ w := by
          rintro rfl
          have := (List.of_mem_filter hq)
          simp only [ne_eq, decide_not, Bool.not_eq_true', decide_eq_false_iff_not] at this
          exact this hqv
        rw [claimedBy_unclaimed_other hvw]
        have hx := hrepl v ⟨q, List.mem_cons_of_mem _ (mem_unclaimed hq), hqv⟩
        rwa [claimedBy_cons_other _ _ (by simp [Ne.symm hvw])] at hx

/-- Any property of gates that the input and every replacement have, the result has: it is
built from nothing else. Both halves of the structural invariant go through this. -/
theorem applyAll_pred {P : Gate → Prop} (repl : Nat → List Gate) (xs : List Tagged)
    (hin : ∀ g ∈ untag xs, P g)
    (hr : ∀ w, (∃ p ∈ xs, p.2 = some w) → ∀ g ∈ repl w, P g) :
    ∀ g ∈ applyAll repl xs, P g := by
  intro g hg
  rcases mem_applyAll repl xs g hg with h | ⟨w, hw, h⟩
  · exact hin g h
  · exact hr w hw g h

/-! ## Deciding the two conditions

The scan that proposes a tagging is unverified, so both conditions are *checked* before the
splice is taken — this is the same certifying arrangement the window search already used, one
level up. Nothing about how a proposal was found is trusted; only that it passes these.

`sepOneB` is one right-to-left pass, carrying whether the rewrite claims anything further
right, so the whole check is linear per rewrite rather than quadratic. `sepB` runs it once per
rewrite rather than once per claimed gate: a later claim's tail is a suffix of the first
claim's, and `SepOne` only ever says less about a suffix. -/

/-- `SepOne`, decided: `.1` is the verdict, `.2` whether `w` claims anything here. -/
def sepOneAux (S : List Qubit) (w : Nat) : List Tagged → Bool × Bool
  | [] => (true, false)
  | (g, t) :: rest =>
      let r := sepOneAux S w rest
      if t = some w then (r.1, true)
      else (r.1 && (!r.2 || g.qubitsOf.all fun q => !S.contains q), r.2)

def sepOneB (S : List Qubit) (w : Nat) (xs : List Tagged) : Bool := (sepOneAux S w xs).1

theorem sepOneAux_seen (S : List Qubit) (w : Nat) :
    ∀ xs : List Tagged, (sepOneAux S w xs).2 = true ↔ ∃ r ∈ xs, r.2 = some w := by
  intro xs
  induction xs with
  | nil => simp [sepOneAux]
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      by_cases ht : t = some w
      · simp [sepOneAux, ht]
      · simp only [sepOneAux, if_neg ht]
        rw [ih]
        constructor
        · rintro ⟨r, hr, hrw⟩; exact ⟨r, List.mem_cons_of_mem _ hr, hrw⟩
        · rintro ⟨r, hr, hrw⟩
          rcases List.mem_cons.1 hr with rfl | hr
          · exact absurd hrw ht
          · exact ⟨r, hr, hrw⟩

theorem sepOneB_sound (S : List Qubit) (w : Nat) :
    ∀ xs : List Tagged, sepOneB S w xs = true → SepOne S w xs := by
  intro xs
  induction xs with
  | nil => intro _; trivial
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      intro h
      by_cases ht : t = some w
      · refine ⟨fun hne _ => absurd ht hne, ih ?_⟩
        simpa [sepOneB, sepOneAux, ht] using h
      · simp only [sepOneB, sepOneAux, if_neg ht, Bool.and_eq_true] at h
        refine ⟨fun _ hex q hq hmem => ?_, ih (by simpa [sepOneB] using h.1)⟩
        have hseen : (sepOneAux S w rest).2 = true := (sepOneAux_seen S w rest).2 hex
        simp only [hseen, Bool.not_true, Bool.false_or, List.all_eq_true] at h
        have hnc := h.2 q hq
        simp only [Bool.not_eq_true'] at hnc
        rw [List.contains_eq_mem, decide_eq_false_iff_not] at hnc
        exact hnc hmem

/-! `sepOneB` has a small proof-friendly specification, but its right-to-left traversal must
visit the whole suffix before it knows where the last claim is.  The executable whole-set
checker already knows each rewrite's member count.  `sepUntilB` consumes that count while
walking left-to-right and stops immediately after the last claim. -/

def sepUntilB (S : List Qubit) (w : Nat) : Nat → List Tagged → Bool
  | 0, _ => true
  | _ + 1, [] => false
  | k + 1, (g, t) :: rest =>
      if t = some w then sepUntilB S w k rest
      else (g.qubitsOf.all fun q => !S.contains q) && sepUntilB S w (k + 1) rest

theorem sepOne_of_no_claims (S : List Qubit) (w : Nat) : ∀ xs,
    claimedBy w xs = [] → SepOne S w xs := by
  intro xs
  induction xs with
  | nil => simp [SepOne]
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      intro h
      by_cases ht : t = some w
      · subst t
        simp [claimedBy] at h
      · have hr : claimedBy w rest = [] := by simpa [claimedBy, ht] using h
        refine ⟨?_, ih hr⟩
        intro _ hex
        obtain ⟨r, hrmem, hrw⟩ := hex
        have hm : r.1 ∈ claimedBy w rest := by
          apply List.mem_filterMap.2
          exact ⟨r, hrmem, by simp [hrw]⟩
        simp [hr] at hm

theorem sepUntilB_sound (S : List Qubit) (w : Nat) : ∀ k xs,
    (claimedBy w xs).length = k → sepUntilB S w k xs = true → SepOne S w xs := by
  intro k xs
  induction xs generalizing k with
  | nil => simp [SepOne]
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      by_cases ht : t = some w
      · subst t
        intro hlen hcheck
        cases k with
        | zero => simp [claimedBy] at hlen
        | succ k =>
            have hrest : (claimedBy w rest).length = k := by
              simpa [claimedBy] using hlen
            refine ⟨fun hne _ => absurd rfl hne, ih k hrest ?_⟩
            simpa [sepUntilB] using hcheck
      · intro hlen hcheck
        have hrest : (claimedBy w rest).length = k := by
          simpa [claimedBy, ht] using hlen
        cases k with
        | zero =>
            have hn : claimedBy w rest = [] := List.eq_nil_of_length_eq_zero hrest
            refine ⟨?_, sepOne_of_no_claims S w rest hn⟩
            intro _ hex
            obtain ⟨r, hr, hrw⟩ := hex
            have hm : r.1 ∈ claimedBy w rest := by
              apply List.mem_filterMap.2
              exact ⟨r, hr, by simp [hrw]⟩
            simp [hn] at hm
        | succ k =>
            simp only [sepUntilB, if_neg ht, Bool.and_eq_true] at hcheck
            refine ⟨?_, ih (k + 1) hrest hcheck.2⟩
            intro _ _ q hq hmem
            have hall := (List.all_eq_true.mp hcheck.1) q hq
            simp only [Bool.not_eq_true'] at hall
            rw [List.contains_eq_mem, decide_eq_false_iff_not] at hall
            exact hall hmem

/-- Count-bounded separation checker.  `members` is the batched claim lookup used by vetting,
so each first claim can stop checking at its last claim instead of traversing the circuit. -/
def sepBoundedB (supp : Nat → List Qubit) (members : Nat → List Gate) :
    List Nat → List Tagged → Bool
  | _, [] => true
  | started, (_, none) :: rest => sepBoundedB supp members started rest
  | started, (_, some w) :: rest =>
      (started.contains w || sepUntilB (supp w) w ((members w).length - 1) rest) &&
        sepBoundedB supp members (w :: started) rest

theorem sepBoundedB_sound (supp : Nat → List Qubit) (members : Nat → List Gate) :
    ∀ started xs,
      (∀ w ∈ started, SepOne (supp w) w xs) →
      (∀ w, w ∉ started → members w = claimedBy w xs) →
      sepBoundedB supp members started xs = true → Sep supp xs := by
  intro started xs
  induction xs generalizing started with
  | nil => intro _ _ _; trivial
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none =>
          intro hst hmem h
          refine ⟨fun w hw => ?_, ih started (fun w hw => (hst w hw).2) ?_
            (by simpa [sepBoundedB] using h)⟩
          · simp at hw
          · intro w hwn
            simpa [claimedBy] using hmem w hwn
      | some v =>
          intro hst hmem h
          simp only [sepBoundedB, Bool.and_eq_true, Bool.or_eq_true] at h
          have hvsep : SepOne (supp v) v rest := by
            by_cases hvstarted : started.contains v
            · exact (hst v (by simpa using hvstarted)).2
            · have hvnot : v ∉ started := by simpa using hvstarted
              have hvcheck : sepUntilB (supp v) v ((members v).length - 1) rest = true := by
                rcases h.1 with hv | hvcheck
                · exact absurd (by simpa using hv) hvnot
                · exact hvcheck
              have hm := hmem v hvnot
              have hl : (members v).length = (claimedBy v rest).length + 1 := by
                simpa [claimedBy] using congrArg List.length hm
              have hlen : (claimedBy v rest).length = (members v).length - 1 := by
                omega
              exact sepUntilB_sound _ _ _ _ hlen hvcheck
          refine ⟨?_, ih (v :: started) ?_ ?_ h.2⟩
          · intro w hw
            have hwv : w = v := (Option.some.inj hw).symm
            subst w
            exact hvsep
          · intro w hw
            simp only [List.mem_cons] at hw
            rcases hw with rfl | hw
            · exact hvsep
            · exact (hst w hw).2
          · intro w hw
            have hwv : w ≠ v := by
              intro heq; subst w; exact hw (by simp)
            have hwn : w ∉ started := fun hws => hw (by simp [hws])
            have hx := hmem w hwn
            rwa [claimedBy_cons_other (w := w) g rest (by simpa using Ne.symm hwv)] at hx

/-- Batched executable separation check for a complete tagging. -/
def sepAllB (supp : Nat → List Qubit) (xs : List Tagged) : Bool :=
  let groups := groupClaimsAux (claimBound xs) xs
  let members : Nat → List Gate := fun w => groups[w]?.getD []
  sepBoundedB supp members [] xs

theorem sepAllB_sound (supp : Nat → List Qubit) (xs : List Tagged)
    (h : sepAllB supp xs = true) : Sep supp xs := by
  let groups := groupClaimsAux (claimBound xs) xs
  let members : Nat → List Gate := fun w => groups[w]?.getD []
  apply sepBoundedB_sound supp members [] xs (by simp) ?_ ?_
  · intro w _
    simpa [members, groups, groupedClaim] using groupedClaim_eq_claimedBy xs w
  · simpa [sepAllB, members, groups] using h

/-- `Sep`, decided — `started` records the rewrites already checked. -/
def sepB (supp : Nat → List Qubit) : List Nat → List Tagged → Bool
  | _, [] => true
  | started, (_, none) :: rest => sepB supp started rest
  | started, (_, some w) :: rest =>
      (started.contains w || sepOneB (supp w) w rest) && sepB supp (w :: started) rest

theorem sepB_sound (supp : Nat → List Qubit) :
    ∀ (xs : List Tagged) (started : List Nat),
      (∀ w ∈ started, SepOne (supp w) w xs) → sepB supp started xs = true → Sep supp xs := by
  intro xs
  induction xs with
  | nil => intro _ _ _; trivial
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      intro started hst h
      cases t with
      | none =>
          refine ⟨fun w hw => absurd hw (by simp), ?_⟩
          exact ih started (fun w hw => (hst w hw).2) (by simpa [sepB] using h)
      | some w =>
          simp only [sepB, Bool.and_eq_true, Bool.or_eq_true] at h
          have hw : SepOne (supp w) w rest := by
            rcases h.1 with hc | hs
            · exact (hst w (by simpa using hc)).2
            · exact sepOneB_sound _ _ _ hs
          refine ⟨fun v hv => ?_, ih (w :: started) (fun v hv => ?_) h.2⟩
          · exact (Option.some.injEq w v ▸ hv) ▸ hw
          · rcases List.mem_cons.1 hv with rfl | hv
            · exact hw
            · exact (hst v hv).2

/-- `OnSupp`, decided. -/
def onSuppB (supp : Nat → List Qubit) (xs : List Tagged) : Bool :=
  xs.all fun p =>
    match p.2 with
    | none => true
    | some w => p.1.qubitsOf.all (fun q => (supp w).contains q) && p.1.isUnitary

theorem onSuppB_sound {supp : Nat → List Qubit} {xs : List Tagged}
    (h : onSuppB supp xs = true) : OnSupp supp xs := by
  intro p hp w hw
  have := (List.all_eq_true.1 h) p hp
  rw [hw] at this
  simp only [Bool.and_eq_true, List.all_eq_true] at this
  exact ⟨fun q hq => List.mem_of_elem_eq_true (this.1 q hq), this.2⟩

/-! ## Dropping the rewrites that do not check out

A proposal is vetted one rewrite at a time: those that fail are untagged, and the rest still
stand. Untagging only ever *removes* obligations, so neither `Sep` nor `OnSupp` has to be
rechecked because of it. -/

/-- Keep only the rewrites `keep` accepts. -/
def vettedBy (keep : Nat → Bool) (xs : List Tagged) : List Tagged :=
  xs.map fun p =>
    match p.2 with
    | none => p
    | some w => if keep w then p else (p.1, none)

@[simp] theorem vettedBy_nil (keep : Nat → Bool) : vettedBy keep [] = [] := rfl

theorem vettedBy_cons_none (keep : Nat → Bool) (g : Gate) (rest : List Tagged) :
    vettedBy keep ((g, none) :: rest) = (g, none) :: vettedBy keep rest := rfl

theorem vettedBy_cons_keep {keep : Nat → Bool} {w : Nat} (h : keep w = true) (g : Gate)
    (rest : List Tagged) :
    vettedBy keep ((g, some w) :: rest) = (g, some w) :: vettedBy keep rest := by
  simp [vettedBy, h]

theorem vettedBy_cons_drop {keep : Nat → Bool} {w : Nat} (h : keep w = false) (g : Gate)
    (rest : List Tagged) :
    vettedBy keep ((g, some w) :: rest) = (g, none) :: vettedBy keep rest := by
  simp [vettedBy, h]

@[simp] theorem untag_vettedBy (keep : Nat → Bool) : ∀ xs : List Tagged,
    untag (vettedBy keep xs) = untag xs := by
  intro xs
  induction xs with
  | nil => rfl
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none => rw [vettedBy_cons_none, untag_cons, untag_cons, ih]
      | some w =>
          by_cases hk : keep w
          · rw [vettedBy_cons_keep hk, untag_cons, untag_cons, ih]
          · rw [vettedBy_cons_drop (by simpa using hk), untag_cons, untag_cons, ih]

theorem claimedBy_vettedBy (keep : Nat → Bool) (w : Nat) : ∀ xs : List Tagged,
    claimedBy w (vettedBy keep xs) = if keep w then claimedBy w xs else [] := by
  intro xs
  induction xs with
  | nil => simp only [vettedBy_nil, claimedBy_nil]; split <;> rfl
  | cons p rest ih =>
      obtain ⟨g, t⟩ := p
      cases t with
      | none =>
          rw [vettedBy_cons_none, claimedBy_cons_other _ _ (by simp),
            claimedBy_cons_other _ _ (by simp), ih]
      | some v =>
          by_cases hk : keep v
          · rw [vettedBy_cons_keep hk]
            by_cases hvw : v = w
            · subst hvw
              rw [claimedBy_cons_self, claimedBy_cons_self, ih, if_pos hk]
              simp [hk]
            · rw [claimedBy_cons_other _ _ (by simpa using hvw),
                claimedBy_cons_other _ _ (by simpa using hvw), ih]
          · rw [vettedBy_cons_drop (by simpa using hk), claimedBy_cons_other _ _ (by simp), ih]
            by_cases hvw : v = w
            · subst hvw
              rw [if_neg (by simpa using hk), if_neg (by simpa using hk)]
            · rw [claimedBy_cons_other _ _ (by simpa using hvw)]

/-- A rewrite that survives vetting is one `keep` accepted. -/
theorem keep_of_mem_vettedBy {keep : Nat → Bool} {w : Nat} {xs : List Tagged}
    (h : ∃ p ∈ vettedBy keep xs, p.2 = some w) : keep w = true := by
  obtain ⟨p, hp, hpw⟩ := h
  rcases List.mem_map.1 hp with ⟨⟨g, t⟩, -, rfl⟩
  cases t with
  | none => exact absurd hpw (by simp)
  | some v =>
      by_cases hk : keep v
      · simp only [hk, if_true, Option.some.injEq] at hpw
        exact hpw ▸ hk
      · simp [hk] at hpw

end TzapLean
