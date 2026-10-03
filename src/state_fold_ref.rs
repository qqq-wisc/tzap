//! The first, straightforward port of state folding, kept as a test oracle for
//! the optimized implementation in `state_fold.rs`: both must rewrite every
//! circuit identically.
#![allow(dead_code)]

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::f64::consts::PI;

use crate::circuit::{Circuit, Gate};


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

/// `Init` before `Temp`, each by index, as in Feynman's `Ord Var`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Var {
    Init(u32),
    Temp(u32),
}

impl Var {
    fn is_temp(self) -> bool {
        matches!(self, Var::Temp(_))
    }
}

/// A multilinear monomial: its variables, ascending. Ordered by Feynman's
/// `lexdegOrd`: graded reverse lexicographic on the path variables first,
/// then on the input variables.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Mono(Vec<Var>);

impl Mono {
    const ONE: Mono = Mono(Vec::new());

    fn var(v: Var) -> Mono {
        Mono(vec![v])
    }

    fn contains(&self, v: Var) -> bool {
        self.0.binary_search(&v).is_ok()
    }

    fn without(&self, v: Var) -> Mono {
        Mono(self.0.iter().copied().filter(|&w| w != v).collect())
    }

    /// The product: the union of the variable sets.
    fn mul(&self, other: &Mono) -> Mono {
        let mut out = Vec::with_capacity(self.0.len() + other.0.len());
        let (mut i, mut j) = (0, 0);
        while i < self.0.len() || j < other.0.len() {
            let next = match (self.0.get(i), other.0.get(j)) {
                (Some(&a), Some(&b)) if a == b => {
                    i += 1;
                    j += 1;
                    a
                }
                (Some(&a), Some(&b)) if a < b => {
                    i += 1;
                    a
                }
                (Some(_), Some(&b)) => {
                    j += 1;
                    b
                }
                (Some(&a), None) => {
                    i += 1;
                    a
                }
                (None, Some(&b)) => {
                    j += 1;
                    b
                }
                (None, None) => unreachable!(),
            };
            out.push(next);
        }
        Mono(out)
    }

    fn split(&self) -> (&[Var], &[Var]) {
        let first_temp = self.0.iter().position(|v| v.is_temp()).unwrap_or(self.0.len());
        (&self.0[first_temp..], &self.0[..first_temp])
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

/// A Boolean polynomial: a set of monomials, summed over F2.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
struct BoolPoly(BTreeSet<Mono>);

impl BoolPoly {
    fn var(v: Var) -> BoolPoly {
        BoolPoly(BTreeSet::from([Mono::var(v)]))
    }

    fn one() -> BoolPoly {
        BoolPoly(BTreeSet::from([Mono::ONE]))
    }

    fn toggle(&mut self, m: Mono) {
        if !self.0.remove(&m) {
            self.0.insert(m);
        }
    }

    fn add(&self, other: &BoolPoly) -> BoolPoly {
        let mut out = self.clone();
        for m in &other.0 {
            out.toggle(m.clone());
        }
        out
    }

    fn mul(&self, other: &BoolPoly) -> BoolPoly {
        let mut out = BoolPoly::default();
        for a in &self.0 {
            for b in &other.0 {
                out.toggle(a.mul(b));
            }
        }
        out
    }

    fn constant(&self) -> bool {
        self.0.contains(&Mono::ONE)
    }

    fn drop_constant(&self) -> BoolPoly {
        let mut out = self.clone();
        out.0.remove(&Mono::ONE);
        out
    }

    fn contains(&self, v: Var) -> bool {
        self.0.iter().any(|m| m.contains(v))
    }

    fn vars(&self) -> impl Iterator<Item = Var> + '_ {
        self.0.iter().flat_map(|m| m.0.iter().copied())
    }

    /// `self[v := p]`.
    fn subst(&self, v: Var, p: &BoolPoly) -> BoolPoly {
        let mut out = BoolPoly::default();
        for m in &self.0 {
            if m.contains(v) {
                let rest = BoolPoly(BTreeSet::from([m.without(v)]));
                for t in rest.mul(p).0 {
                    out.toggle(t);
                }
            } else {
                out.toggle(m.clone());
            }
        }
        out
    }

    /// The first solution of `p = 0` that Feynman's `solveForX` lists and
    /// `ok` accepts: in monomial order, a variable `u` whose linear term is its
    /// only occurrence, with `p + u` and that polynomial's degree.
    fn first_solution(&self, ok: impl Fn(Var, isize) -> bool) -> Option<(Var, BoolPoly)> {
        let mut occurrences: HashMap<Var, usize> = HashMap::default();
        let (mut top, mut top_count, mut second) = (-1isize, 0usize, -1isize);
        for m in &self.0 {
            for &v in &m.0 {
                *occurrences.entry(v).or_default() += 1;
            }
            let d = m.0.len() as isize;
            if d > top {
                second = top;
                top = d;
                top_count = 1;
            } else if d == top {
                top_count += 1;
            } else if d > second {
                second = d;
            }
        }
        // The degree of `p` without one linear term.
        let rest_degree = if top == 1 && top_count == 1 { second } else { top };
        let m = self
            .0
            .iter()
            .filter(|m| m.0.len() == 1)
            .find(|m| occurrences[&m.0[0]] == 1 && ok(m.0[0], rest_degree))?;
        let mut rest = self.clone();
        rest.0.remove(m);
        Some((m.0[0], rest))
    }
}

/// `a` distributed over a Boolean polynomial into a pseudo-Boolean one:
/// `a (m_1 + ... + m_k)` is the sum, over nonempty subsets `S`, of
/// `(-2)^(|S|-1) a` times the product of `S`. This is Feynman's recursion
/// `a(m + xs) = a m + a xs - 2a m xs`, unrolled; a subset whose coefficient
/// is zero (mod 2 pi) ends the branch, as there.
fn distribute(a: Angle, p: &BoolPoly) -> Vec<(Mono, Angle)> {
    let ms: Vec<&Mono> = p.0.iter().collect();
    let mut out = Vec::new();
    distribute_from(a, &ms, 0, &Mono::ONE, &mut out);
    out
}

fn distribute_from(a: Angle, ms: &[&Mono], start: usize, prod: &Mono, out: &mut Vec<(Mono, Angle)>) {
    if a.is_zero() {
        return;
    }
    let next = a.times(-2);
    for i in start..ms.len() {
        let p = prod.mul(ms[i]);
        if i + 1 < ms.len() {
            distribute_from(next, ms, i + 1, &p, out);
        }
        out.push((p, a));
    }
}

/// The phase polynomial, with an index from variables to monomials and the
/// set of variables whose terms changed since it was last drained.
#[derive(Default)]
struct PhasePoly {
    terms: HashMap<Mono, Angle>,
    by_var: HashMap<Var, HashSet<Mono>>,
    dirty: HashSet<Var>,
}

impl PhasePoly {
    fn add_term(&mut self, m: Mono, a: Angle) {
        self.dirty.extend(m.0.iter().copied());
        let sum = match self.terms.get(&m) {
            Some(&b) => b.add(a),
            None => a,
        };
        if sum.is_zero() {
            if self.terms.remove(&m).is_some() {
                self.unindex(&m);
            }
        } else if self.terms.insert(m.clone(), sum).is_none() {
            for &v in &m.0 {
                self.by_var.entry(v).or_default().insert(m.clone());
            }
        }
    }

    fn unindex(&mut self, m: &Mono) {
        for v in &m.0 {
            if let Some(set) = self.by_var.get_mut(v) {
                set.remove(m);
                if set.is_empty() {
                    self.by_var.remove(v);
                }
            }
        }
    }

    fn add(&mut self, terms: Vec<(Mono, Angle)>) {
        for (m, a) in terms {
            self.add_term(m, a);
        }
    }

    /// Remove and return the terms containing `v`.
    fn take_var(&mut self, v: Var) -> Vec<(Mono, Angle)> {
        let Some(ms) = self.by_var.remove(&v) else { return Vec::new() };
        let mut out = Vec::with_capacity(ms.len());
        for m in ms {
            let a = self.terms.remove(&m).expect("indexed term");
            self.unindex(&m);
            self.dirty.extend(m.0.iter().copied());
            out.push((m, a));
        }
        out
    }

    /// Feynman's `elimVar`: drop the terms containing `v`.
    fn remove_var(&mut self, v: Var) {
        self.take_var(v);
    }

    fn subst(&mut self, v: Var, p: &BoolPoly) {
        for (m, a) in self.take_var(v) {
            let rest = BoolPoly(BTreeSet::from([m.without(v)]));
            self.add(distribute(a, &rest.mul(p)));
        }
    }

    /// `toBooleanPoly (quotVar v pp)`, with `constant` added to the constant
    /// term first: the quotient as a Boolean polynomial, if every coefficient
    /// is pi.
    fn boolean_quotient(&self, v: Var, constant: Option<Angle>) -> Option<BoolPoly> {
        let ms = self.by_var.get(&v)?;
        let unit = Mono::var(v);
        let mut has_unit = false;
        for m in ms {
            let a = self.terms[m];
            if *m == unit && constant.is_some() {
                has_unit = true;
                let sum = a.add(constant.unwrap());
                if !(sum.is_zero() || sum == Angle::PI) {
                    return None;
                }
            } else if a != Angle::PI {
                return None;
            }
        }
        if let Some(c) = constant {
            if !has_unit && c != Angle::PI {
                return None;
            }
        }
        let mut out = BoolPoly::default();
        for m in ms {
            let q = m.without(v);
            if q == Mono::ONE && constant.is_some() {
                if self.terms[m].add(constant.unwrap()).is_zero() {
                    continue;
                }
            }
            out.0.insert(q);
        }
        if let Some(c) = constant {
            if !has_unit && !c.is_zero() {
                out.0.insert(Mono::ONE);
            }
        }
        Some(out)
    }
}

// ---------------------------------------------------------------------------
// The analysis.

/// A mergeable phase: the gates `(location, parity)` that contributed to it,
/// and its total angle.
#[derive(Clone, Default)]
struct Term {
    locs: BTreeSet<(usize, bool)>,
    angle: Option<Angle>,
}

struct Ctx {
    temps: u32,
    ket: Vec<BoolPoly>,
    /// For each variable, how many qubits' states contain it.
    in_ket: HashMap<Var, usize>,
    terms: HashMap<BoolPoly, Term>,
    terms_by_var: HashMap<Var, HashSet<BoolPoly>>,
    pp: PhasePoly,
}

impl Ctx {
    fn new(n: usize) -> Self {
        Ctx {
            temps: 0,
            ket: (0..n as u32).map(|i| BoolPoly::var(Var::Init(i))).collect(),
            in_ket: (0..n as u32).map(|i| (Var::Init(i), 1)).collect(),
            terms: HashMap::default(),
            terms_by_var: HashMap::default(),
            pp: PhasePoly::default(),
        }
    }

    fn fresh(&mut self) -> Var {
        self.temps += 1;
        Var::Temp(self.temps - 1)
    }

    /// Set qubit `i`'s state, keeping `in_ket` current. A variable that
    /// enters or leaves the state changes whether it may be reduced.
    fn set_ket(&mut self, i: usize, p: BoolPoly) {
        let old: BTreeSet<Var> = self.ket[i].vars().collect();
        let new: BTreeSet<Var> = p.vars().collect();
        for &v in old.difference(&new) {
            let c = self.in_ket.get_mut(&v).expect("counted");
            *c -= 1;
            if *c == 0 {
                self.in_ket.remove(&v);
                self.pp.dirty.insert(v);
            }
        }
        for &v in new.difference(&old) {
            *self.in_ket.entry(v).or_default() += 1;
            self.pp.dirty.insert(v);
        }
        self.ket[i] = p;
    }

    /// Merge `term` into the term with key `key`.
    fn merge_term(&mut self, key: BoolPoly, term: Term) {
        if !self.terms.contains_key(&key) {
            for v in key.vars().collect::<BTreeSet<_>>() {
                self.terms_by_var.entry(v).or_default().insert(key.clone());
            }
        }
        let entry = self.terms.entry(key).or_default();
        entry.locs.extend(term.locs);
        entry.angle = match (entry.angle, term.angle) {
            (Some(a), Some(b)) => Some(a.add(b)),
            (a, b) => a.or(b),
        };
    }

    fn add_term(&mut self, theta: Angle, loc: usize, bexp: &BoolPoly) {
        let parity = bexp.constant();
        let theta2 = if parity { theta.neg() } else { theta };
        let term = Term {
            locs: BTreeSet::from([(loc, parity)]),
            angle: Some(theta2),
        };
        self.merge_term(bexp.drop_constant(), term);
        self.pp.add(distribute(theta, bexp));
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
                self.pp.add(distribute(Angle::PI, &p));
            }
            Gate::x(v) => {
                let p = BoolPoly::one().add(&self.ket[q(v)]);
                self.set_ket(q(v), p);
            }
            Gate::h(v) => {
                let y = BoolPoly::var(self.fresh());
                let p = self.ket[q(v)].mul(&y);
                self.pp.add(distribute(Angle::PI, &p));
                self.set_ket(q(v), y);
            }
            // Feynman's catch-all: the operands become fresh, unreducible
            // path variables.
            _ => {
                let (len, operands) = crate::circuit::qubit_operands(g);
                for &v in &operands[..len] {
                    let y = self.fresh();
                    self.pp.add(vec![(Mono::var(y), Angle::C(2f64.sqrt()))]);
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
        let moved = self.terms_by_var.remove(&x).unwrap_or_default();
        let mut renamed = Vec::with_capacity(moved.len());
        for key in moved {
            let mut term = self.terms.remove(&key).expect("indexed term");
            for v in key.vars().collect::<BTreeSet<_>>() {
                if let Some(set) = self.terms_by_var.get_mut(&v) {
                    set.remove(&key);
                    if set.is_empty() {
                        self.terms_by_var.remove(&v);
                    }
                }
            }
            if flip {
                term.locs = term.locs.into_iter().map(|(l, par)| (l, !par)).collect();
                term.angle = term.angle.map(Angle::neg);
            }
            renamed.push((key.subst(x, p).drop_constant(), term));
        }
        for (key, term) in renamed {
            self.merge_term(key, term);
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
        x.is_temp() && self.pp.by_var.contains_key(&x) && !self.in_ket.contains_key(&x)
    }

    /// Feynman's `applyReductions`: [HH] while one matches, else [omega],
    /// always at the smallest variable that matches. The matches of each
    /// variable are cached, and only variables whose terms or state
    /// membership changed are looked at again.
    fn reduce(&mut self, cutoff: Option<usize>) {
        let mut hh: BTreeMap<Var, (Var, BoolPoly)> = BTreeMap::new();
        let mut omega: BTreeMap<Var, BoolPoly> = BTreeMap::new();
        let mut todo: Vec<Var> = self.pp.by_var.keys().copied().collect();
        self.pp.dirty.clear();
        loop {
            for x in todo.drain(..) {
                hh.remove(&x);
                omega.remove(&x);
                if !self.is_candidate(x) {
                    continue;
                }
                if let Some(q) = self.pp.boolean_quotient(x, None) {
                    let ok = |u: Var, d: isize| u.is_temp() && cutoff.map_or(true, |c| d <= c as isize);
                    if let Some(sol) = q.first_solution(ok) {
                        hh.insert(x, sol);
                    }
                }
                if let Some(q) = self.pp.boolean_quotient(x, Some(Angle::dyadic(3, 1))) {
                    omega.insert(x, q);
                }
            }
            if let Some((x, (y, sub))) = hh.pop_first() {
                self.pp.remove_var(x);
                self.subst(y, &sub);
            } else if let Some((x, q)) = omega.pop_first() {
                self.pp.remove_var(x);
                self.pp.add(vec![(Mono::ONE, Angle::dyadic(1, 2))]);
                self.pp.add(distribute(Angle::dyadic(3, 1), &q));
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
    for term in ctx.terms.values() {
        let Some(angle) = term.angle else { continue };
        let &(first, parity) = term.locs.first().expect("a term has a location");
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

