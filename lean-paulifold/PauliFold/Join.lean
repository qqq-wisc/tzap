import PauliFold.Tableau

/-!
# Branches and loops: tableaux with unknown rows

To join the analyses of two branches, the tableau is allowed to *forget* rows.

**The state** (`PTab`). A partial tableau: each row is known (`some`) or unknown (`none`),
together with the axes. The whole domain is `PAbs n = Option (PTab n)`, where `none` is `⊤`.

**Join.** A row stays known when both branches agree on it, and becomes unknown otherwise.
The axes are the union of both branches' axes. All axes are in the input coordinates of the
common starting point, so a fact survives the join only if it commutes with the rotations
of both branches.

**Transformers.**

- A Clifford gate multiplies rows as in `Tab.step`; a product involving an unknown row is
  unknown.
- A T or T† on `q` records `z[q]` as an axis. If `z[q]` is unknown, the state goes to `⊤`.

**Facts.** `β(P)` is defined only when every row `P` needs is known; facts come only from
such `P`.

**Soundness via completions.** A partial state is sound for `U` when some *total* tableau
extends it (agrees on every known row, and has axes among its axes) and satisfies
`Tab.Inv U`. Every partial step is matched by the total step on the completion, so
soundness reduces to `Tab.inv_step`.

**Programs.** Gates, sequencing, nondeterministic choice and loops. A loop iterates
`I ← I ⊔ body♯(I)` until `body♯(I) ⊑ I`, which it checks decidably; after `2n + 3` tries it
gives `⊤`.

Main results: `absRun_sound` and `program_sound`, with examples checked by `decide`.
-/

namespace PauliFold

open TzapLean Matrix

noncomputable section

variable {n : ℕ}

/-! ## Partial tableaux -/

/-- A tableau whose rows may be unknown. -/
structure PTab (n : ℕ) where
  rx : Fin n → Option (PP n)
  rz : Fin n → Option (PP n)
  axes : List (PP n)

instance : DecidableEq (PTab n) := fun a b =>
  decidable_of_iff (a.rx = b.rx ∧ a.rz = b.rz ∧ a.axes = b.axes) (by cases a; cases b; simp)

/-- The product of two possibly unknown strings: unknown if either is. -/
def mulO (a b : Option (PP n)) : Option (PP n) := a.bind fun x => b.map fun y => x.mul y

namespace PTab

def init (n : ℕ) : PTab n := ⟨fun q => some (PP.single q .X), fun q => some (PP.single q .Z), []⟩

/-- The row for a single-qubit letter, if known. -/
def rowOf (T : PTab n) (r : Fin n) : Letter → Option (PP n)
  | .I => some PP.one
  | .X => T.rx r
  | .Z => T.rz r
  | .Y => (mulO (T.rx r) (T.rz r)).map (·.scale 1)

def backL (T : PTab n) (l : List (Fin n)) (P : PP n) : Option (PP n) :=
  l.foldr (fun r acc => mulO (T.rowOf r (P.letter r)) acc) (some PP.one)

/-- The pull-back `β(P)`, defined when every row `P` needs is known. -/
def back (T : PTab n) (P : PP n) : Option (PP n) :=
  (T.backL (List.finRange n) P).map (·.scale P.phase)

/-- The Clifford row updates, with unknown rows propagating. -/
def stepC (T : PTab n) : CTGate n → PTab n
  | .h q => { T with rx := Function.update T.rx q (T.rz q), rz := Function.update T.rz q (T.rx q) }
  | .x q => { T with rz := Function.update T.rz q ((T.rz q).map (·.scale 2)) }
  | .z q => { T with rx := Function.update T.rx q ((T.rx q).map (·.scale 2)) }
  | .s q => { T with rx := Function.update T.rx q ((mulO (T.rx q) (T.rz q)).map (·.scale 3)) }
  | .sdg q => { T with rx := Function.update T.rx q ((mulO (T.rx q) (T.rz q)).map (·.scale 1)) }
  | .cx c t _ => { T with
      rx := Function.update T.rx c (mulO (T.rx c) (T.rx t))
      rz := Function.update T.rz t (mulO (T.rz c) (T.rz t)) }
  | .cz c t _ => { T with
      rx := Function.update (Function.update T.rx c (mulO (T.rx c) (T.rz t))) t
        (mulO (T.rz c) (T.rx t)) }
  | .t _ | .tdg _ => T

/-- The transformer of a gate. A rotation whose axis row is unknown gives `⊤`. -/
def stepP (g : CTGate n) (T : PTab n) : Option (PTab n) :=
  match g with
  | .t q => (T.rz q).map fun A => { T with axes := A :: T.axes }
  | .tdg q => (T.rz q).map fun A => { T with axes := A :: T.axes }
  | g => some (T.stepC g)

/-! ## Completions -/

/-- The total tableau `T'` agrees with `T` on every known row. -/
def Extends (T' : Tab n) (T : PTab n) : Prop :=
  ∀ r A, (T.rx r = some A → T'.rx r = A) ∧ (T.rz r = some A → T'.rz r = A)

theorem mulO_eq_some {a b : Option (PP n)} {A : PP n} (h : mulO a b = some A) :
    ∃ x y, a = some x ∧ b = some y ∧ A = x.mul y := by
  cases a <;> cases b <;> simp_all [mulO]

theorem rowOf_extends {T' : Tab n} {T : PTab n} (hE : Extends T' T) (r : Fin n) (L : Letter)
    {A : PP n} (h : T.rowOf r L = some A) : T'.rowOf r L = A := by
  cases L with
  | I => simp_all [rowOf, Tab.rowOf]
  | X => exact (hE r A).1 h
  | Z => exact (hE r A).2 h
  | Y =>
    simp only [rowOf, Option.map_eq_some_iff] at h
    obtain ⟨B, hB, rfl⟩ := h
    obtain ⟨x, y, hx, hy, rfl⟩ := mulO_eq_some hB
    simp [Tab.rowOf, (hE r x).1 hx, (hE r y).2 hy]

theorem backL_extends {T' : Tab n} {T : PTab n} (hE : Extends T' T) (P : PP n) :
    ∀ (l : List (Fin n)) {B : PP n}, T.backL l P = some B → T'.backL l P = B
  | [], B, h => by simp_all [backL, Tab.backL]
  | r :: l, B, h => by
    simp only [backL, List.foldr_cons] at h
    obtain ⟨x, y, hx, hy, rfl⟩ := mulO_eq_some h
    have ih := backL_extends hE P l (B := y) hy
    simp only [Tab.backL, List.foldr_cons] at ih ⊢
    rw [rowOf_extends hE r _ hx, ih]

/-- **Known pull-backs are the completion's pull-backs.** -/
theorem back_extends {T' : Tab n} {T : PTab n} (hE : Extends T' T) {P B : PP n}
    (h : T.back P = some B) : T'.back P = B := by
  simp only [back, Option.map_eq_some_iff] at h
  obtain ⟨B', hB', rfl⟩ := h
  simp [Tab.back, backL_extends hE P _ hB']

/-- **The Clifford updates agree with the total ones on known rows.** -/
theorem stepC_extends {T' : Tab n} {T : PTab n} (hE : Extends T' T) (g : CTGate n) :
    Extends (T'.step g) (T.stepC g) := by
  intro r A
  have hx := fun s B => (hE s B).1
  have hz := fun s B => (hE s B).2
  cases g with
  | t q => exact hE r A
  | tdg q => exact hE r A
  | h q =>
    constructor <;> intro h <;> simp only [stepC, Tab.step, Function.update_apply] at h ⊢ <;>
      split_ifs at h ⊢ with hr <;> first | exact hz _ A h | exact hx _ A h
  | x q =>
    refine ⟨fun h => by simpa [stepC, Tab.step] using hx r A h, fun h => ?_⟩
    simp only [stepC, Tab.step, Function.update_apply] at h ⊢
    split_ifs at h ⊢ with hr
    · subst hr
      simp only [Option.map_eq_some_iff] at h
      obtain ⟨B, hB, rfl⟩ := h
      rw [hz _ _ hB]
    · exact hz r A h
  | z q =>
    refine ⟨fun h => ?_, fun h => by simpa [stepC, Tab.step] using hz r A h⟩
    simp only [stepC, Tab.step, Function.update_apply] at h ⊢
    split_ifs at h ⊢ with hr
    · subst hr
      simp only [Option.map_eq_some_iff] at h
      obtain ⟨B, hB, rfl⟩ := h
      rw [hx _ _ hB]
    · exact hx r A h
  | s q =>
    refine ⟨fun h => ?_, fun h => by simpa [stepC, Tab.step] using hz r A h⟩
    simp only [stepC, Tab.step, Function.update_apply] at h ⊢
    split_ifs at h ⊢ with hr
    · subst hr
      simp only [Option.map_eq_some_iff] at h
      obtain ⟨B, hB, rfl⟩ := h
      obtain ⟨u, v, hu, hv, rfl⟩ := mulO_eq_some hB
      rw [hx _ _ hu, hz _ _ hv]
    · exact hx r A h
  | sdg q =>
    refine ⟨fun h => ?_, fun h => by simpa [stepC, Tab.step] using hz r A h⟩
    simp only [stepC, Tab.step, Function.update_apply] at h ⊢
    split_ifs at h ⊢ with hr
    · subst hr
      simp only [Option.map_eq_some_iff] at h
      obtain ⟨B, hB, rfl⟩ := h
      obtain ⟨u, v, hu, hv, rfl⟩ := mulO_eq_some hB
      rw [hx _ _ hu, hz _ _ hv]
    · exact hx r A h
  | cx c t hne =>
    constructor <;> intro h <;> simp only [stepC, Tab.step, Function.update_apply] at h ⊢
    · split_ifs at h ⊢ with hr
      · obtain ⟨u, v, hu, hv, rfl⟩ := mulO_eq_some h
        subst hr; rw [hx _ _ hu, hx _ _ hv]
      · exact hx r A h
    · split_ifs at h ⊢ with hr
      · obtain ⟨u, v, hu, hv, rfl⟩ := mulO_eq_some h
        subst hr; rw [hz _ _ hu, hz _ _ hv]
      · exact hz r A h
  | cz c t hne =>
    constructor <;> intro h <;> simp only [stepC, Tab.step, Function.update_apply] at h ⊢
    · split_ifs at h ⊢ with h1 h2
      · obtain ⟨u, v, hu, hv, rfl⟩ := mulO_eq_some h
        rw [hz _ _ hu, hx _ _ hv]
      · obtain ⟨u, v, hu, hv, rfl⟩ := mulO_eq_some h
        rw [hx _ _ hu, hz _ _ hv]
      · exact hx r A h
    · exact hz r A h

end PTab

/-! ## The domain with `⊤` -/

abbrev PAbs (n : ℕ) := Option (PTab n)

/-- Soundness for `U`: `⊤`, or a completion satisfying the tableau invariant. -/
def Sound (U : Density n) : PAbs n → Prop
  | none => True
  | some T => ∃ T' : Tab n, PTab.Extends T' T ∧ (∀ A ∈ T'.axes, A ∈ T.axes) ∧ Tab.Inv U T'

/-- The facts of a partial state, and its concretization. -/
def gammaP (T : PTab n) : Set (Density n) :=
  {U | ∀ P B : PP n, T.back P = some B → (∀ A ∈ T.axes, PP.Comm B A) →
    conj U B.toMatrix = P.toMatrix}

def gammaA : PAbs n → Set (Density n)
  | none => Set.univ
  | some T => gammaP T

/-- The tableau invariant gives the facts, for a total tableau. -/
theorem fact_of_inv {U : Density n} {T : Tab n} (h : Tab.Inv U T) (P : PP n)
    (hP : ∀ A ∈ T.axes, PP.Comm (T.back P) A) : conj U (T.back P).toMatrix = P.toMatrix := by
  obtain ⟨C, V, hU, hCC, hVV, hrows, hcomm⟩ := h
  have hC : Cᴴ * C = 1 := mul_eq_one_comm.mp hCC
  have hB := hcomm _ hP
  rw [hU, conj, Matrix.conjTranspose_mul]
  simp only [Matrix.mul_assoc]
  rw [← Matrix.mul_assoc V, hB, Matrix.mul_assoc, ← Matrix.mul_assoc V Vᴴ, hVV, Matrix.one_mul,
    Tab.back_ok hC hCC hrows]
  simp only [Matrix.mul_assoc]
  rw [hCC, Matrix.mul_one, ← Matrix.mul_assoc, hCC, Matrix.one_mul]

theorem gammaA_of_sound {U : Density n} {a : PAbs n} (h : Sound U a) : U ∈ gammaA a := by
  cases a with
  | none => trivial
  | some T =>
    obtain ⟨T', hE, hax, hinv⟩ := h
    intro P B hb hc
    have hb' := PTab.back_extends hE hb
    have := fact_of_inv hinv P (fun A hA => by rw [hb']; exact hc A (hax A hA))
    rwa [hb'] at this

/-- The gate transformer on `PAbs`. -/
def stepA (g : CTGate n) (a : PAbs n) : PAbs n := a.bind (PTab.stepP g)

theorem stepA_sound {U : Density n} {a : PAbs n} (h : Sound U a) (g : CTGate n) :
    Sound (gateUnitary n g.toGate * U) (stepA g a) := by
  cases a with
  | none => trivial
  | some T =>
    obtain ⟨T', hE, hax, hinv⟩ := h
    have hinv' := Tab.inv_step hinv g
    by_cases hT : ∃ q, g = .t q ∨ g = .tdg q
    · obtain ⟨q, hq⟩ := hT
      have hstep : T'.step g = { T' with axes := T'.rz q :: T'.axes } := by
        rcases hq with rfl | rfl <;> rfl
      cases hz : T.rz q with
      | none =>
        show Sound _ (PTab.stepP g T)
        rcases hq with rfl | rfl <;> simp [PTab.stepP, hz, Sound]
      | some A =>
        have hA : T'.rz q = A := ((hE q A).2 hz)
        have hres : PTab.stepP g T = some { T with axes := A :: T.axes } := by
          rcases hq with rfl | rfl <;> simp [PTab.stepP, hz]
        show Sound _ (PTab.stepP g T)
        rw [hres]
        refine ⟨T'.step g, fun r B => by rw [hstep]; exact hE r B, fun B hB => ?_, hinv'⟩
        rw [hstep] at hB
        rcases List.mem_cons.1 hB with rfl | hB
        · rw [hA]; exact List.mem_cons_self ..
        · exact List.mem_cons_of_mem _ (hax B hB)
    · push Not at hT
      have hcl : g.IsClifford := by
        cases g with
        | t q => exact ((hT q).1 rfl).elim
        | tdg q => exact ((hT q).2 rfl).elim
        | _ => trivial
      have hres : PTab.stepP g T = some (T.stepC g) := by
        cases g <;> first | rfl | exact hcl.elim
      show Sound _ (PTab.stepP g T)
      rw [hres]
      have haxC : (T.stepC g).axes = T.axes := by
        cases g <;> first | rfl | exact hcl.elim
      refine ⟨T'.step g, PTab.stepC_extends hE g, fun B hB => ?_, hinv'⟩
      rw [Tab.axes_step_clifford hcl] at hB
      rw [haxC]; exact hax B hB

/-! ## Order and join -/

/-- `X ⊑ Y`: every known row of `Y` is the same row of `X`, and `X`'s axes are among `Y`'s. -/
def LeP (X Y : PTab n) : Prop :=
  (∀ r, Y.rx r = none ∨ Y.rx r = X.rx r) ∧ (∀ r, Y.rz r = none ∨ Y.rz r = X.rz r) ∧
    ∀ A ∈ X.axes, A ∈ Y.axes

instance (X Y : PTab n) : Decidable (LeP X Y) := by unfold LeP; infer_instance

def LeA : PAbs n → PAbs n → Prop
  | _, none => True
  | none, some _ => False
  | some X, some Y => LeP X Y

instance (a b : PAbs n) : Decidable (LeA a b) := by
  cases a <;> cases b <;> unfold LeA <;> infer_instance

theorem sound_of_le {U : Density n} {a b : PAbs n} (h : Sound U a) (hle : LeA a b) :
    Sound U b := by
  cases b with
  | none => trivial
  | some Y =>
    cases a with
    | none => exact hle.elim
    | some X =>
      obtain ⟨hx, hz, hax⟩ := hle
      obtain ⟨T', hE, hax', hinv⟩ := h
      refine ⟨T', fun r A => ⟨fun h => ?_, fun h => ?_⟩, fun A hA => hax A (hax' A hA), hinv⟩
      · rcases hx r with h' | h'
        · rw [h'] at h; cases h
        · exact (hE r A).1 (h'.symm.trans h)
      · rcases hz r with h' | h'
        · rw [h'] at h; cases h
        · exact (hE r A).2 (h'.symm.trans h)

/-- A row survives the join when both branches agree on it. -/
def rowJ (a b : Option (PP n)) : Option (PP n) := if a = b then a else none

/-- **The join**: agreeing rows are kept, others become unknown; the axes are unioned. -/
def joinP (X Y : PTab n) : PTab n :=
  ⟨fun r => rowJ (X.rx r) (Y.rx r), fun r => rowJ (X.rz r) (Y.rz r), X.axes ++ Y.axes⟩

def join : PAbs n → PAbs n → PAbs n
  | some X, some Y => some (joinP X Y)
  | _, _ => none

theorem le_join_left (a b : PAbs n) : LeA a (join a b) := by
  cases a with
  | none => cases b <;> trivial
  | some X =>
    cases b with
    | none => trivial
    | some Y =>
      refine ⟨fun r => ?_, fun r => ?_, fun A hA => List.mem_append_left _ hA⟩ <;>
        simp only [joinP, rowJ] <;> split_ifs <;> simp

theorem le_join_right (a b : PAbs n) : LeA b (join a b) := by
  cases a with
  | none => cases b <;> trivial
  | some X =>
    cases b with
    | none => trivial
    | some Y =>
      refine ⟨fun r => ?_, fun r => ?_, fun A hA => List.mem_append_right _ hA⟩ <;>
        simp only [joinP, rowJ] <;> split_ifs with h <;> simp [h]

theorem join_sound_left {U : Density n} {a : PAbs n} (h : Sound U a) (b : PAbs n) :
    Sound U (join a b) := sound_of_le h (le_join_left a b)

theorem join_sound_right {U : Density n} (a : PAbs n) {b : PAbs n} (h : Sound U b) :
    Sound U (join a b) := sound_of_le h (le_join_right a b)

/-! ## Programs -/

inductive Prog (n : ℕ) where
  | gate (g : CTGate n)
  | seq (p q : Prog n)
  | choice (p q : Prog n)
  | loop (p : Prog n)

/-- The collecting semantics: `Exec p U` when `U` is one of the unitaries `p` implements. -/
inductive Exec : Prog n → Density n → Prop
  | gate (g : CTGate n) : Exec (.gate g) (gateUnitary n g.toGate)
  | seq {p q : Prog n} {U V : Density n} : Exec p U → Exec q V → Exec (.seq p q) (V * U)
  | choiceL {p q : Prog n} {U : Density n} : Exec p U → Exec (.choice p q) U
  | choiceR {p q : Prog n} {U : Density n} : Exec q U → Exec (.choice p q) U
  | loopZero {p : Prog n} : Exec (.loop p) 1
  | loopStep {p : Prog n} {U V : Density n} : Exec (.loop p) U → Exec p V → Exec (.loop p) (V * U)

/-- Loop iteration: `I ← I ⊔ f(I)` until `f(I) ⊑ I`, with a bound on the number of tries. -/
def loopFix (f : PAbs n → PAbs n) : ℕ → PAbs n → PAbs n
  | 0, _ => none
  | k + 1, a =>
    if LeA (f (join a (f a))) (join a (f a)) then join a (f a) else loopFix f k (join a (f a))

theorem loopFix_spec (f : PAbs n → PAbs n) :
    ∀ (k : ℕ) (a : PAbs n), (∀ U, Sound U a → Sound U (loopFix f k a)) ∧
      LeA (f (loopFix f k a)) (loopFix f k a)
  | 0, a => ⟨fun _ _ => trivial, by cases f none <;> trivial⟩
  | k + 1, a => by
    simp only [loopFix]
    split_ifs with hle
    · exact ⟨fun U h => join_sound_left h _, hle⟩
    · obtain ⟨h1, h2⟩ := loopFix_spec f k (join a (f a))
      exact ⟨fun U h => h1 U (join_sound_left h _), h2⟩

/-- The abstract semantics. -/
def absRun : Prog n → PAbs n → PAbs n
  | .gate g, a => stepA g a
  | .seq p q, a => absRun q (absRun p a)
  | .choice p q, a => join (absRun p a) (absRun q a)
  | .loop p, a => loopFix (absRun p) (2 * n + 3) a

/-- **Soundness of the abstract semantics.** -/
theorem absRun_sound {p : Prog n} {U : Density n} (hU : Exec p U) :
    ∀ {a : PAbs n} {U₀ : Density n}, Sound U₀ a → Sound (U * U₀) (absRun p a) := by
  induction hU with
  | gate g => intro a U₀ h; exact stepA_sound h g
  | seq _ _ ihp ihq =>
    intro a U₀ h
    rw [Matrix.mul_assoc]
    exact ihq (ihp h)
  | choiceL _ ih => intro a U₀ h; exact join_sound_left (ih h) _
  | choiceR _ ih => intro a U₀ h; exact join_sound_right _ (ih h)
  | @loopZero p =>
    intro a U₀ h
    rw [Matrix.one_mul]
    exact (loopFix_spec (absRun p) _ a).1 _ h
  | @loopStep p U V _ _ ihloop ihbody =>
    intro a U₀ h
    have hI : Sound (U * U₀) (loopFix (absRun p) (2 * n + 3) a) := ihloop h
    rw [Matrix.mul_assoc]
    exact sound_of_le (ihbody hI) (loopFix_spec (absRun p) _ a).2

theorem sound_init : Sound (1 : Density n) (some (PTab.init n)) :=
  ⟨Tab.init n, fun r A => ⟨fun h => by cases h; rfl, fun h => by cases h; rfl⟩,
    fun _ h => absurd h List.not_mem_nil, Tab.inv_init⟩

/-- **Programs are soundly abstracted**: every unitary a program can implement, over all
branch choices and loop iterations, lies in the concretization of its abstract run. -/
theorem program_sound {p : Prog n} {U : Density n} (hU : Exec p U) :
    U ∈ gammaA (absRun p (some (PTab.init n))) := by
  have := absRun_sound hU sound_init
  rw [Matrix.mul_one] at this
  exact gammaA_of_sound this

/-- A fact read off the abstract result holds for every execution. -/
theorem fact_of_run {p : Prog n} {U : Density n} (hU : Exec p U) {T : PTab n}
    (hT : absRun p (some (PTab.init n)) = some T) {P B : PP n} (hb : T.back P = some B)
    (hc : ∀ A ∈ T.axes, PP.Comm B A) : conj U B.toMatrix = P.toMatrix := by
  have h := program_sound hU
  rw [hT] at h
  exact h P B hb hc

/-! ## Examples: the Amy–Lunderville loops -/

/-- `loop-simple`: `t a; while ★ { cx a,b }; tdg a`. -/
def loopSimple : Prog 2 :=
  .seq (.gate (.t 0)) (.seq (.loop (.gate (.cx 0 1 (by decide)))) (.gate (.tdg 0)))

theorem loopSimple_some : (absRun loopSimple (some (PTab.init 2))).isSome := by decide

/-- The loop forgets the rows it changes (`x_a`, `z_b`) but keeps `z_a`, so every execution
of `loop-simple` preserves `Z_a`, whatever the number of iterations. -/
theorem loopSimple_Z {U : Density 2} (hU : Exec loopSimple U) :
    conj U (PP.single (0 : Fin 2) .Z).toMatrix = (PP.single (0 : Fin 2) .Z).toMatrix :=
  fact_of_run hU (Option.some_get loopSimple_some).symm (by decide) (by decide)

/-- `loop-h`: `t b; while ★ { h a }; tdg b`. The loop forgets `a`'s rows and keeps `b`'s. -/
def loopH : Prog 2 :=
  .seq (.gate (.t 1)) (.seq (.loop (.gate (.h 0))) (.gate (.tdg 1)))

theorem loopH_some : (absRun loopH (some (PTab.init 2))).isSome := by decide

theorem loopH_Z {U : Density 2} (hU : Exec loopH U) :
    conj U (PP.single (1 : Fin 2) .Z).toMatrix = (PP.single (1 : Fin 2) .Z).toMatrix :=
  fact_of_run hU (Option.some_get loopH_some).symm (by decide) (by decide)

/-- `loop-swap` with the swap as three CNOTs: the loop forgets all four rows, and the
following `tdg b` needs the unknown `z_b`, so the result is `⊤`. -/
def loopSwap : Prog 2 :=
  .seq (.seq (.gate (.cx 0 1 (by decide))) (.seq (.gate (.t 1)) (.gate (.cx 0 1 (by decide)))))
    (.seq (.loop (.seq (.gate (.cx 0 1 (by decide)))
      (.seq (.gate (.cx 1 0 (by decide))) (.gate (.cx 0 1 (by decide))))))
      (.seq (.gate (.cx 0 1 (by decide))) (.seq (.gate (.tdg 1)) (.gate (.cx 0 1 (by decide))))))

example : absRun loopSwap (some (PTab.init 2)) = none := by decide

end

end PauliFold
