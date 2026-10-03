//! State folding: an experimental port of Feynman's `-statefold <d>`
//! (`Feynman.Optimization.StateFold`, Amy 2024), for comparing it with
//! PhaseFoldPauli without depending on Feynman's implementation.
//!
//! The pass symbolically executes the circuit as a path sum. Each qubit holds
//! a Boolean (F2, multilinear) polynomial over input variables `x_i` and path
//! variables `y_j`, one per Hadamard. Every phase gate is recorded twice: as a
//! *term*, keyed by its Boolean polynomial, which remembers the gates that
//! contributed to it, and in the pseudo-Boolean *phase polynomial* of the whole
//! path sum. After the circuit, path-sum reductions ([HH] and [omega]) are
//! applied while they match. An [HH] reduction substitutes a path variable by a
//! polynomial of degree at most `d` (unbounded when `d` is `None`), and term
//! keys that become equal merge. Each merged term is then re-synthesized at its
//! first gate and the other gates are deleted.
//!
//! The port follows Feynman's order exactly (candidate variables in ascending
//! order, solutions in its monomial order, linear reductions first), so it
//! makes the same choices. That includes its rule for a substitution
//! `x := p`: a term's angle is negated when `p` has constant 1, rather than
//! when the term's new key does. The two differ only when `x` occurs in the key
//! solely inside products; randomized unitary checks have not found a case
//! where Feynman's rule changes the circuit.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use smallvec::SmallVec;

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::f64::consts::PI;

use crate::circuit::{Circuit, Gate};
use crate::pass::Pass;

/// State folding up to degree `degree` (`None`: unbounded).
pub struct StateFold {
    pub degree: Option<usize>,
}

impl Pass for StateFold {
    fn name(&self) -> &str {
        "State folding"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        state_fold(circuit, self.degree)
    }
}

// ---------------------------------------------------------------------------
// Angles: Feynman's `Angle`, a dyadic multiple of pi mod 2, or a float.

#[derive(Clone, Copy, Debug)]
enum Angle {
    /// `num / 2^k` times pi, with `0 <= num < 2^(k+1)` and `num` odd unless
    /// `k == 0`.
    D { num: i64, k: u32 },
    /// Radians.
    C(f64),
}

impl PartialEq for Angle {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Angle::D { num: a, k: i }, Angle::D { num: b, k: j }) => a == b && i == j,
            (Angle::C(a), Angle::C(b)) => a == b,
            _ => false,
        }
    }
}

impl Angle {
    const ZERO: Angle = Angle::D { num: 0, k: 0 };
    const PI: Angle = Angle::D { num: 1, k: 0 };

    fn dyadic(num: i64, k: u32) -> Angle {
        let modulus = 1i64 << (k + 1);
        let (mut num, mut k) = (num.rem_euclid(modulus), k);
        while k > 0 && num % 2 == 0 {
            num /= 2;
            k -= 1;
        }
        Angle::D { num, k }
    }

    fn radians(self) -> f64 {
        match self {
            Angle::D { num, k } => PI * num as f64 / (1u64 << k) as f64,
            Angle::C(t) => t,
        }
    }

    fn add(self, other: Angle) -> Angle {
        match (self, other) {
            (Angle::D { num: a, k: i }, Angle::D { num: b, k: j }) => {
                let k = i.max(j);
                Angle::dyadic((a << (k - i)) + (b << (k - j)), k)
            }
            (a, b) => Angle::C(a.radians() + b.radians()),
        }
    }

    fn neg(self) -> Angle {
        match self {
            Angle::D { num, k } => Angle::dyadic(-num, k),
            Angle::C(t) => Angle::C(-t),
        }
    }

    /// `i` times the angle (Feynman's `power`).
    fn times(self, i: i64) -> Angle {
        match self {
            Angle::D { num, k } => Angle::dyadic(num.wrapping_mul(i), k),
            Angle::C(t) => Angle::C(i as f64 * t),
        }
    }

    fn is_zero(self) -> bool {
        self == Angle::ZERO
    }
}

// ---------------------------------------------------------------------------
// Variables, monomials and polynomials, in Feynman's orders.

/// A variable: input `x_i` is `i`, path variable `y_j` is `TEMP | j`, so the
/// numeric order is Feynman's `Ord Var` (inputs first, each by index).
type Var = u32;
const TEMP: Var = 1 << 31;

fn is_temp(v: Var) -> bool {
    v & TEMP != 0
}

/// A multilinear monomial: its variables, ascending. Ordered by Feynman's
/// `lexdegOrd`: graded reverse lexicographic on the path variables first,
/// then on the input variables.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Default)]
struct Mono(SmallVec<[Var; 4]>, u16);

impl Mono {
    /// The monomial of an ascending variable list. The second field caches
    /// the number of input variables, which the order compares on.
    fn new(vs: SmallVec<[Var; 4]>) -> Mono {
        let inits = vs.partition_point(|&v| !is_temp(v)) as u16;
        Mono(vs, inits)
    }

    fn one() -> Mono {
        Mono(SmallVec::new(), 0)
    }

    fn var(v: Var) -> Mono {
        Mono(smallvec::smallvec![v], u16::from(!is_temp(v)))
    }

    fn is_one(&self) -> bool {
        self.0.is_empty()
    }

    fn contains(&self, v: Var) -> bool {
        self.0.binary_search(&v).is_ok()
    }

    fn without(&self, v: Var) -> Mono {
        Mono::new(self.0.iter().copied().filter(|&w| w != v).collect())
    }

    /// The product: the union of the variable sets.
    fn mul(&self, other: &Mono) -> Mono {
        let (a, b) = (&self.0, &other.0);
        let mut out = SmallVec::with_capacity(a.len() + b.len());
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].cmp(&b[j]) {
                Ordering::Less => {
                    out.push(a[i]);
                    i += 1;
                }
                Ordering::Greater => {
                    out.push(b[j]);
                    j += 1;
                }
                Ordering::Equal => {
                    out.push(a[i]);
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&a[i..]);
        out.extend_from_slice(&b[j..]);
        Mono::new(out)
    }

    /// (path variables, input variables).
    fn split(&self) -> (&[Var], &[Var]) {
        let k = self.1 as usize;
        (&self.0[k..], &self.0[..k])
    }
}

fn grevlex(a: &[Var], b: &[Var]) -> Ordering {
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

impl Ord for Mono {
    fn cmp(&self, other: &Self) -> Ordering {
        let (at, ai) = self.split();
        let (bt, bi) = other.split();
        grevlex(at, bt).then_with(|| grevlex(ai, bi))
    }
}

impl PartialOrd for Mono {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A Boolean polynomial over F2, split by degree: the constant, the linear
/// variables (ascending), and the monomials of degree at least 2 (sorted,
/// distinct). Term keys and states are almost always affine, and an affine
/// substitution is then one merge of two integer lists.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Default)]
struct BoolPoly {
    one: bool,
    lin: Vec<Var>,
    high: Vec<Mono>,
}

/// The symmetric difference of two sorted, distinct lists.
fn xor_sorted<T: Ord + Clone>(a: &[T], b: &[T]) -> Vec<T> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Less => {
                out.push(a[i].clone());
                i += 1;
            }
            Ordering::Greater => {
                out.push(b[j].clone());
                j += 1;
            }
            Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Sort a list and cancel equal pairs (a sum over F2).
fn cancel_pairs<T: Ord>(mut ms: Vec<T>) -> Vec<T> {
    ms.sort_unstable();
    let mut out: Vec<T> = Vec::with_capacity(ms.len());
    for m in ms {
        if out.last() == Some(&m) {
            out.pop();
        } else {
            out.push(m);
        }
    }
    out
}

impl BoolPoly {
    fn var(v: Var) -> BoolPoly {
        BoolPoly { one: false, lin: vec![v], high: Vec::new() }
    }

    fn one() -> BoolPoly {
        BoolPoly { one: true, ..BoolPoly::default() }
    }

    fn is_affine(&self) -> bool {
        self.high.is_empty()
    }

    /// The polynomial of a list of monomials, summed over F2.
    fn from_list(ms: Vec<Mono>) -> BoolPoly {
        let mut out = BoolPoly::default();
        let mut lin = Vec::new();
        let mut high = Vec::new();
        for m in ms {
            match m.0.len() {
                0 => out.one = !out.one,
                1 => lin.push(m.0[0]),
                _ => high.push(m),
            }
        }
        out.lin = cancel_pairs(lin);
        out.high = cancel_pairs(high);
        out
    }

    /// The monomials, in Feynman's order.
    fn monos(&self) -> Vec<Mono> {
        let mut ms: Vec<Mono> = Vec::with_capacity(self.len());
        if self.one {
            ms.push(Mono::one());
        }
        ms.extend(self.lin.iter().map(|&v| Mono::var(v)));
        ms.extend(self.high.iter().cloned());
        if !self.high.is_empty() {
            ms.sort_unstable();
        }
        ms
    }

    fn len(&self) -> usize {
        usize::from(self.one) + self.lin.len() + self.high.len()
    }

    fn add(&self, other: &BoolPoly) -> BoolPoly {
        BoolPoly {
            one: self.one ^ other.one,
            lin: xor_sorted(&self.lin, &other.lin),
            high: xor_sorted(&self.high, &other.high),
        }
    }

    fn mul(&self, other: &BoolPoly) -> BoolPoly {
        let (a, b) = (self.monos(), other.monos());
        let mut ms = Vec::with_capacity(a.len() * b.len());
        for x in &a {
            for y in &b {
                ms.push(x.mul(y));
            }
        }
        BoolPoly::from_list(ms)
    }

    /// `m * self` for a monomial `m`.
    fn mul_mono(&self, m: &Mono) -> BoolPoly {
        if m.is_one() {
            return self.clone();
        }
        BoolPoly::from_list(self.monos().iter().map(|a| a.mul(m)).collect())
    }

    fn constant(&self) -> bool {
        self.one
    }

    fn drop_constant(mut self) -> BoolPoly {
        self.one = false;
        self
    }

    fn contains(&self, v: Var) -> bool {
        self.lin.binary_search(&v).is_ok() || self.high.iter().any(|m| m.contains(v))
    }

    /// The variables, ascending, without repeats.
    fn vars(&self) -> Vec<Var> {
        if self.high.is_empty() {
            return self.lin.clone();
        }
        let mut vs: Vec<Var> = self.lin.clone();
        vs.extend(self.high.iter().flat_map(|m| m.0.iter().copied()));
        vs.sort_unstable();
        vs.dedup();
        vs
    }

    /// `self := self[v := p]`, in place when both are affine.
    fn subst_in_place(&mut self, v: Var, p: &BoolPoly) {
        if self.is_affine() && p.is_affine() {
            if let Ok(i) = self.lin.binary_search(&v) {
                self.lin.remove(i);
                self.one ^= p.one;
                for &u in &p.lin {
                    match self.lin.binary_search(&u) {
                        Ok(j) => {
                            self.lin.remove(j);
                        }
                        Err(j) => self.lin.insert(j, u),
                    }
                }
            }
            return;
        }
        *self = self.subst(v, p);
    }

    /// `self[v := p]`.
    fn subst(&self, v: Var, p: &BoolPoly) -> BoolPoly {
        let in_lin = self.lin.binary_search(&v);
        let in_high = self.high.iter().any(|m| m.contains(v));
        if !in_high {
            let Ok(i) = in_lin else { return self.clone() };
            if p.is_affine() {
                // v appears once, linearly: replace it by p.
                let mut lin = self.lin.clone();
                lin.remove(i);
                return BoolPoly {
                    one: self.one ^ p.one,
                    lin: xor_sorted(&lin, &p.lin),
                    high: xor_sorted(&self.high, &p.high),
                };
            }
        }
        let mut kept = Vec::with_capacity(self.len());
        let mut products = Vec::new();
        let pm = p.monos();
        for m in self.monos() {
            if m.contains(v) {
                let rest = m.without(v);
                products.extend(pm.iter().map(|t| t.mul(&rest)));
            } else {
                kept.push(m);
            }
        }
        BoolPoly::from_list(kept).add(&BoolPoly::from_list(products))
    }

    /// The first solution of `p = 0` that Feynman's `solveForX` lists and
    /// `ok` accepts: in monomial order, a variable `u` whose linear term is its
    /// only occurrence, with `p + u` and that polynomial's degree.
    fn first_solution(&self, ok: impl Fn(Var, isize) -> bool) -> Option<(Var, BoolPoly)> {
        let mut in_high: HashSet<Var> = HashSet::default();
        let mut top = if self.one { 0isize } else { -1 };
        for m in &self.high {
            in_high.extend(m.0.iter().copied());
            top = top.max(m.0.len() as isize);
        }
        // The degree of `p` without one linear term.
        let rest_degree = if !self.high.is_empty() {
            top
        } else if self.lin.len() > 1 {
            1
        } else {
            top
        };
        // Linear terms in Feynman's order: inputs, then path variables, each
        // ascending; that is the numeric order.
        let i = self
            .lin
            .iter()
            .position(|&u| !in_high.contains(&u) && ok(u, rest_degree))?;
        let mut rest = self.clone();
        let u = rest.lin.remove(i);
        Some((u, rest))
    }
}

/// `a` distributed over a Boolean polynomial into a pseudo-Boolean one:
/// `a (m_1 + ... + m_k)` is the sum, over nonempty subsets `S`, of
/// `(-2)^(|S|-1) a` times the product of `S`. This is Feynman's recursion
/// `a(m + xs) = a m + a xs - 2a m xs`, unrolled; a subset whose coefficient
/// is zero (mod 2 pi) ends the branch, as there.
fn distribute_into(a: Angle, p: &BoolPoly, pp: &mut PhasePoly) {
    fn go(a: Angle, ms: &[Mono], start: usize, prod: &Mono, pp: &mut PhasePoly) {
        if a.is_zero() {
            return;
        }
        let next = a.times(-2);
        for i in start..ms.len() {
            let p = prod.mul(&ms[i]);
            if i + 1 < ms.len() {
                go(next, ms, i + 1, &p, pp);
            }
            pp.add_term(p, a);
        }
    }
    go(a, &p.monos(), 0, &Mono::one(), pp);
}

/// The phase polynomial. Each monomial present is stored once, in a slab, by
/// id; `by_var` lists, for each variable, the ids of monomials containing it
/// (entries of removed monomials are skipped and dropped lazily), and
/// `count` how many present monomials contain it. `dirty` collects the
/// variables whose monomials changed since it was last drained.
#[derive(Default)]
struct PhasePoly {
    ids: HashMap<Mono, u32>,
    slab: Vec<(Mono, Angle, bool)>,
    free: Vec<u32>,
    by_var: HashMap<Var, Vec<u32>>,
    count: HashMap<Var, u32>,
    dirty: HashSet<Var>,
}

impl PhasePoly {
    fn contains_var(&self, v: Var) -> bool {
        self.count.contains_key(&v)
    }

    fn vars(&self) -> Vec<Var> {
        self.count.keys().copied().collect()
    }

    fn add_term(&mut self, m: Mono, a: Angle) {
        use std::collections::hash_map::Entry;
        self.dirty.extend(m.0.iter().copied());
        match self.ids.entry(m) {
            Entry::Occupied(e) => {
                let id = *e.get();
                let slot = &mut self.slab[id as usize];
                slot.1 = slot.1.add(a);
                if slot.1.is_zero() {
                    e.remove();
                    self.release(id);
                }
            }
            Entry::Vacant(e) => {
                if a.is_zero() {
                    return;
                }
                let m = e.key().clone();
                let id = match self.free.pop() {
                    Some(id) => {
                        self.slab[id as usize] = (m.clone(), a, true);
                        id
                    }
                    None => {
                        self.slab.push((m.clone(), a, true));
                        (self.slab.len() - 1) as u32
                    }
                };
                e.insert(id);
                for &v in &m.0 {
                    self.by_var.entry(v).or_default().push(id);
                    *self.count.entry(v).or_default() += 1;
                }
            }
        }
    }

    /// Free slot `id` (already removed from `ids`).
    fn release(&mut self, id: u32) {
        let slot = &mut self.slab[id as usize];
        slot.2 = false;
        for &v in &slot.0 .0 {
            let c = self.count.get_mut(&v).expect("counted");
            *c -= 1;
            if *c == 0 {
                self.count.remove(&v);
                self.by_var.remove(&v);
            }
        }
        // The slot is reused only after `by_var` lists that may still name it
        // are rebuilt; see `live`.
        self.free.push(id);
    }

    /// The live monomials containing `v`, compacting its list. A reused slot
    /// is recognized by its monomial no longer containing `v`.
    fn live(&mut self, v: Var) -> Vec<u32> {
        let Some(list) = self.by_var.get_mut(&v) else { return Vec::new() };
        let slab = &self.slab;
        list.retain(|&id| {
            let (m, _, alive) = &slab[id as usize];
            *alive && m.contains(v)
        });
        list.sort_unstable();
        list.dedup();
        list.clone()
    }

    /// Remove and return the terms containing `v`.
    fn take_var(&mut self, v: Var) -> Vec<(Mono, Angle)> {
        let ids = self.live(v);
        let mut out = Vec::with_capacity(ids.len());
        for id in ids {
            let (m, a, _) = self.slab[id as usize].clone();
            self.ids.remove(&m);
            self.dirty.extend(m.0.iter().copied());
            self.release(id);
            out.push((m, a));
        }
        out
    }

    fn subst(&mut self, v: Var, p: &BoolPoly) {
        for (m, a) in self.take_var(v) {
            let q = p.mul_mono(&m.without(v));
            distribute_into(a, &q, self);
        }
    }

    /// `toBooleanPoly (quotVar v pp)`, with `constant` added to the constant
    /// term first: the quotient as a Boolean polynomial, if every coefficient
    /// is pi.
    fn boolean_quotient(&mut self, v: Var, constant: Option<Angle>) -> Option<BoolPoly> {
        let ids = self.live(v);
        if ids.is_empty() {
            return None;
        }
        let mut constant_sum = constant;
        for &id in &ids {
            let (m, a, _) = &self.slab[id as usize];
            if m.0.len() == 1 && constant.is_some() {
                constant_sum = Some(a.add(constant.unwrap()));
            } else if *a != Angle::PI {
                return None;
            }
        }
        if let Some(c) = constant_sum {
            if !(c.is_zero() || c == Angle::PI) {
                return None;
            }
        }
        let mut out: Vec<Mono> = ids
            .iter()
            .map(|&id| &self.slab[id as usize].0)
            .filter(|m| !(m.0.len() == 1 && constant.is_some()))
            .map(|m| m.without(v))
            .collect();
        if constant_sum.is_some_and(|c| c == Angle::PI) {
            out.push(Mono::one());
        }
        Some(BoolPoly::from_list(out))
    }
}

// ---------------------------------------------------------------------------
// The analysis.

/// A mergeable phase: its key, the gates `(location, parity)` that
/// contributed to it, and its total angle. A term merged into another is dead.
struct Term {
    key: BoolPoly,
    hash: u64,
    locs: Vec<(usize, bool)>,
    angle: Option<Angle>,
    alive: bool,
}

struct Ctx {
    temps: u32,
    ket: Vec<BoolPoly>,
    /// For each variable, how many qubits' states contain it.
    in_ket: HashMap<Var, usize>,
    terms: Vec<Term>,
    /// The live terms with each key hash.
    term_of_key: HashMap<u64, SmallVec<[usize; 1]>>,
    /// Terms whose key contained the variable when they were listed; stale
    /// entries are skipped on use.
    terms_by_var: HashMap<Var, Vec<usize>>,
    pp: PhasePoly,
}

impl Ctx {
    fn new(n: usize) -> Self {
        Ctx {
            temps: 0,
            ket: (0..n as u32).map(BoolPoly::var).collect(),
            in_ket: (0..n as u32).map(|i| (i, 1)).collect(),
            terms: Vec::new(),
            term_of_key: HashMap::default(),
            terms_by_var: HashMap::default(),
            pp: PhasePoly::default(),
        }
    }

    fn fresh(&mut self) -> Var {
        self.temps += 1;
        TEMP | (self.temps - 1)
    }

    /// Set qubit `i`'s state, keeping `in_ket` current. A variable that
    /// enters or leaves the state changes whether it may be reduced.
    fn set_ket(&mut self, i: usize, p: BoolPoly) {
        let old = self.ket[i].vars();
        let new = p.vars();
        for &v in &old {
            if new.binary_search(&v).is_err() {
                let c = self.in_ket.get_mut(&v).expect("counted");
                *c -= 1;
                if *c == 0 {
                    self.in_ket.remove(&v);
                    self.pp.dirty.insert(v);
                }
            }
        }
        for &v in &new {
            if old.binary_search(&v).is_err() {
                *self.in_ket.entry(v).or_default() += 1;
                self.pp.dirty.insert(v);
            }
        }
        self.ket[i] = p;
    }

    fn key_hash(key: &BoolPoly) -> u64 {
        use std::hash::BuildHasher;
        rustc_hash::FxBuildHasher.hash_one(key)
    }

    /// The live term with key `key`, if any.
    fn find(&self, key: &BoolPoly, hash: u64) -> Option<usize> {
        self.term_of_key
            .get(&hash)?
            .iter()
            .copied()
            .find(|&id| self.terms[id].key == *key)
    }

    /// Give term `id` (keyless) the key `key`, merging it into the term
    /// already there. A new key is indexed under `new_vars`: every variable
    /// of `key` that the term is not already listed under.
    fn rekey(&mut self, id: usize, key: BoolPoly, new_vars: &[Var]) {
        let hash = Self::key_hash(&key);
        if let Some(other) = self.find(&key, hash) {
            let locs = std::mem::take(&mut self.terms[id].locs);
            let angle = self.terms[id].angle.take();
            self.terms[id].alive = false;
            let t = &mut self.terms[other];
            t.locs.extend(locs);
            t.angle = match (t.angle, angle) {
                (Some(a), Some(b)) => Some(a.add(b)),
                (a, b) => a.or(b),
            };
            return;
        }
        for &v in new_vars {
            self.terms_by_var.entry(v).or_default().push(id);
        }
        self.term_of_key.entry(hash).or_default().push(id);
        self.terms[id].key = key;
        self.terms[id].hash = hash;
    }

    /// Forget term `id`'s key in the key map.
    fn unkey(&mut self, id: usize) {
        let hash = self.terms[id].hash;
        let ids = self.term_of_key.get_mut(&hash).expect("keyed term");
        ids.retain(|&mut i| i != id);
        if ids.is_empty() {
            self.term_of_key.remove(&hash);
        }
    }

    fn add_term(&mut self, theta: Angle, loc: usize, bexp: &BoolPoly) {
        let parity = bexp.constant();
        let theta2 = if parity { theta.neg() } else { theta };
        let key = bexp.clone().drop_constant();
        if let Some(id) = self.find(&key, Self::key_hash(&key)) {
            let t = &mut self.terms[id];
            t.locs.push((loc, parity));
            t.angle = Some(t.angle.map_or(theta2, |a| theta2.add(a)));
        } else {
            let id = self.terms.len();
            self.terms.push(Term {
                key: BoolPoly::default(),
                hash: 0,
                locs: vec![(loc, parity)],
                angle: Some(theta2),
                alive: true,
            });
            let vars = key.vars();
            self.rekey(id, key, &vars);
        }
        distribute_into(theta, bexp, &mut self.pp);
    }

    fn gate(&mut self, g: &Gate, loc: usize) {
        let q = |q: u32| q as usize;
        match *g {
            Gate::t(v) => self.add_term(Angle::dyadic(1, 2), loc, &self.ket[q(v)].clone()),
            Gate::tdg(v) => self.add_term(Angle::dyadic(7, 2), loc, &self.ket[q(v)].clone()),
            Gate::s(v) => self.add_term(Angle::dyadic(1, 1), loc, &self.ket[q(v)].clone()),
            Gate::sdg(v) => self.add_term(Angle::dyadic(3, 1), loc, &self.ket[q(v)].clone()),
            Gate::z(v) => self.add_term(Angle::PI, loc, &self.ket[q(v)].clone()),
            Gate::rz(theta, v) => self.add_term(Angle::C(theta), loc, &self.ket[q(v)].clone()),
            Gate::cnot { control, target } => {
                let p = self.ket[q(control)].add(&self.ket[q(target)]);
                self.set_ket(q(target), p);
            }
            Gate::cz { control, target } => {
                let p = self.ket[q(control)].mul(&self.ket[q(target)]);
                distribute_into(Angle::PI, &p, &mut self.pp);
            }
            Gate::x(v) => {
                let p = BoolPoly::one().add(&self.ket[q(v)]);
                self.set_ket(q(v), p);
            }
            Gate::h(v) => {
                let y = BoolPoly::var(self.fresh());
                let p = self.ket[q(v)].mul(&y);
                distribute_into(Angle::PI, &p, &mut self.pp);
                self.set_ket(q(v), y);
            }
            // Feynman's catch-all: the operands become fresh, unreducible
            // path variables.
            _ => {
                let (len, operands) = crate::circuit::qubit_operands(g);
                for &v in &operands[..len] {
                    let y = self.fresh();
                    self.pp.add_term(Mono::var(y), Angle::C(2f64.sqrt()));
                    self.set_ket(q(v), BoolPoly::var(y));
                }
            }
        }
    }

    /// Substitute `x := p` in the phase polynomial, the term keys and the
    /// state.
    fn subst(&mut self, x: Var, p: &BoolPoly) {
        self.pp.subst(x, p);
        // Feynman's rule: negate when the substituted polynomial has constant 1.
        let flip = p.constant();
        let mut ids = self.terms_by_var.remove(&x).unwrap_or_default();
        ids.sort_unstable();
        ids.dedup();
        let p_vars: Vec<Var> = p.vars();
        let mut renamed = Vec::with_capacity(ids.len());
        for id in ids {
            if !self.terms[id].alive || !self.terms[id].key.contains(x) {
                continue;
            }
            self.unkey(id);
            let t = &mut self.terms[id];
            if flip {
                for (_, par) in &mut t.locs {
                    *par = !*par;
                }
                t.angle = t.angle.map(Angle::neg);
            }
            let mut key = std::mem::take(&mut t.key);
            key.subst_in_place(x, p);
            renamed.push((id, key.drop_constant()));
        }
        // Rekeyed after all keys are removed, so that two terms merging here
        // find each other however they are ordered.
        for (id, key) in renamed {
            self.rekey(id, key, &p_vars);
        }
        for i in 0..self.ket.len() {
            if self.ket[i].contains(x) {
                let new = self.ket[i].subst(x, p);
                self.set_ket(i, new);
            }
        }
    }

    /// Whether `x` may be reduced: a path variable in the phase polynomial
    /// and not in the state.
    fn is_candidate(&self, x: Var) -> bool {
        is_temp(x) && self.pp.contains_var(x) && !self.in_ket.contains_key(&x)
    }

    /// Feynman's `applyReductions`: [HH] while one matches, else [omega],
    /// always at the smallest variable that matches. The matches of each
    /// variable are cached, and only variables whose terms or state
    /// membership changed are looked at again.
    fn reduce(&mut self, cutoff: Option<usize>) {
        let mut hh: BTreeMap<Var, (Var, BoolPoly)> = BTreeMap::new();
        let mut omega: BTreeMap<Var, BoolPoly> = BTreeMap::new();
        let mut todo: Vec<Var> = self.pp.vars();
        self.pp.dirty.clear();
        let ok = |u: Var, d: isize| is_temp(u) && cutoff.is_none_or(|c| d <= c as isize);
        loop {
            for x in todo.drain(..) {
                hh.remove(&x);
                omega.remove(&x);
                if !self.is_candidate(x) {
                    continue;
                }
                if let Some(q) = self.pp.boolean_quotient(x, None) {
                    if let Some(sol) = q.first_solution(ok) {
                        hh.insert(x, sol);
                    }
                }
                if let Some(q) = self.pp.boolean_quotient(x, Some(Angle::dyadic(3, 1))) {
                    omega.insert(x, q);
                }
            }
            if let Some((x, (y, sub))) = hh.pop_first() {
                self.pp.take_var(x);
                self.subst(y, &sub);
            } else if let Some((x, q)) = omega.pop_first() {
                self.pp.take_var(x);
                self.pp.add_term(Mono::one(), Angle::dyadic(1, 2));
                distribute_into(Angle::dyadic(3, 1), &q, &mut self.pp);
            } else {
                return;
            }
            todo.extend(self.pp.dirty.drain());
        }
    }
}

/// The new angle of each phase gate location (absent: unchanged).
fn analyze(circuit: &Circuit, degree: Option<usize>) -> HashMap<usize, Angle> {
    let mut ctx = Ctx::new(circuit.num_qubits);
    for (loc, g) in circuit.gates.iter().enumerate() {
        ctx.gate(g, loc);
    }
    ctx.reduce(Some(1));
    match degree {
        Some(1) => {}
        Some(d) => ctx.reduce(Some(d)),
        None => ctx.reduce(None),
    }
    let mut out = HashMap::default();
    for term in ctx.terms.iter().filter(|t| t.alive) {
        let Some(angle) = term.angle else { continue };
        let &(first, parity) = term.locs.iter().min().expect("a term has a location");
        let angle = if parity { angle.neg() } else { angle };
        for &(l, _) in &term.locs {
            out.insert(l, if l == first { angle } else { Angle::ZERO });
        }
    }
    out
}

/// Feynman's `synthesizePhase`.
fn synthesize(q: u32, angle: Angle, out: &mut Vec<Gate>) {
    match angle {
        Angle::C(t) => out.push(Gate::rz(t, q)),
        Angle::D { num: 0, .. } => {}
        Angle::D { k: 0, .. } => out.push(Gate::z(q)),
        Angle::D { num, k: 1 } => out.push(if num % 4 == 1 { Gate::s(q) } else { Gate::sdg(q) }),
        Angle::D { num, k: 2 } => match num % 8 {
            1 => out.push(Gate::t(q)),
            3 => out.extend([Gate::tdg(q), Gate::z(q)]),
            5 => out.extend([Gate::t(q), Gate::z(q)]),
            _ => out.push(Gate::tdg(q)),
        },
        a => out.push(Gate::rz(a.radians(), q)),
    }
}

pub fn state_fold(circuit: &Circuit, degree: Option<usize>) -> Circuit {
    if circuit
        .gates
        .iter()
        .any(|g| matches!(g, Gate::measure { .. } | Gate::reset(_)))
    {
        return circuit.clone();
    }
    let angles = analyze(circuit, degree);
    let mut gates = Vec::with_capacity(circuit.gates.len());
    for (loc, g) in circuit.gates.iter().enumerate() {
        match (angles.get(&loc), g) {
            (
                Some(&a),
                Gate::t(q) | Gate::tdg(q) | Gate::s(q) | Gate::sdg(q) | Gate::z(q) | Gate::rz(_, q),
            ) => synthesize(*q, a, &mut gates),
            _ => gates.push(g.clone()),
        }
    }
    Circuit {
        gates,
        ..circuit.clone()
    }
}

#[cfg(test)]
mod tests_support {
    use super::*;

    pub struct Rng(pub u64);

    impl Rng {
        pub fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        pub fn up_to(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    pub fn random_circuit(rng: &mut Rng, n: usize, len: usize) -> Circuit {
        let mut c = Circuit::with_cbits(n, 0);
        for _ in 0..len {
            let q = rng.up_to(n) as u32;
            let r = ((q as usize + 1 + rng.up_to(n - 1)) % n) as u32;
            c.apply(match rng.up_to(10) {
                0 | 1 => Gate::h(q),
                2 | 3 => Gate::t(q),
                4 => Gate::tdg(q),
                5 => Gate::s(q),
                6 => Gate::x(q),
                7 => Gate::rz(0.3 * (1 + rng.up_to(3)) as f64, q),
                8 => Gate::cz { control: q, target: r },
                _ => Gate::cnot { control: q, target: r },
            });
        }
        c
    }

}

#[cfg(test)]
mod tests {
    use super::tests_support::*;
    use super::*;
    use crate::unitary::circuits_equiv;

    /// Every degree preserves the unitary up to
    /// global phase, and removes T gates on some circuits.
    #[test]
    fn random_circuits_preserve_unitary() {
        let mut rng = Rng(0x5747_e_f01d_0001);
        let mut folded = 0;
        for case in 0..400 {
            let n = 2 + case % 3;
            let len = 8 + rng.up_to(30);
            let c = random_circuit(&mut rng, n, len);
            for degree in [Some(1), Some(2), Some(3), None] {
                let out = state_fold(&c, degree);
                assert!(circuits_equiv(&c, &out, 1e-8), "case {case}, degree {degree:?}");
                folded += (out.gates.len() < c.gates.len()) as usize;
            }
        }
        assert!(folded > 100, "too few folds: {folded}");
    }

    #[test]
    fn t_gates_merge_across_a_hadamard_pair() {
        // T q0; H q1; H q1; T q0 -> S q0 (the HH pair cancels).
        let mut c = Circuit::with_cbits(2, 0);
        for g in [Gate::t(0), Gate::h(1), Gate::h(1), Gate::t(0)] {
            c.apply(g);
        }
        let out = state_fold(&c, Some(1));
        assert!(circuits_equiv(&c, &out, 1e-9));
        assert_eq!(out.gates.iter().filter(|g| matches!(g, Gate::t(_))).count(), 0);
    }
}


#[cfg(test)]
mod oracle {
    use super::tests_support::*;
    use super::*;

    /// The optimized implementation rewrites exactly as the first port did.
    #[test]
    fn matches_the_reference_port() {
        let mut rng = Rng(0x0dd5_eed5_0001);
        for case in 0..2000 {
            let n = 2 + case % 5;
            let len = 8 + rng.up_to(60);
            let c = random_circuit(&mut rng, n, len);
            for degree in [Some(1), Some(2), Some(3), None] {
                let a = state_fold(&c, degree);
                let b = crate::state_fold_ref::state_fold(&c, degree);
                assert_eq!(format!("{:?}", a.gates), format!("{:?}", b.gates), "case {case}, degree {degree:?}");
            }
        }
    }
}
