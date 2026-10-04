//! StateFold over programs: a port of Feynman's `stateFoldpp`
//! (`-qasm3 -statefold <d>`), the state folding of loops and branches.
//!
//! Each branch or loop body is analyzed from a fresh state, as a block. Its
//! effect is summarized as a transition ideal over the pre-state `X` and
//! post-state `X'` (a reduced Gröbner basis over F2, with path variables
//! eliminated); a loop's summary is the Kleene closure of its body's. The
//! summary then fast-forwards the enclosing state: the post-state becomes new
//! path variables constrained by the summary, and the body's phase polynomial,
//! reduced modulo the summary and scaled by an irrational constant (so it never
//! enables a reduction), is added. Phase terms inside a body never merge with
//! terms outside it (they become *orphans*). The port keeps Feynman's orders of
//! iteration, so it makes the same choices.

use super::*;

/// The `Prime` variables (post-states in summaries), ordered between inputs
/// and path variables, as in Feynman.
pub(super) const PRIME: Var = 1 << 30;

fn is_prime(v: Var) -> bool {
    v & PRIME != 0 && !is_temp(v)
}

// ---------------------------------------------------------------------------
// Gröbner bases over F2[X]/(x^2 - x), as in Feynman's `Groebner` module.

/// Feynman's `lexdegOrd`, with `elim` choosing the eliminated variables.
fn cmp_with(a: &Mono, b: &Mono, elim: &dyn Fn(Var) -> bool) -> Ordering {
    let split = |m: &Mono| -> (Vec<Var>, Vec<Var>) { m.0.iter().partition(|&&v| elim(v)) };
    let (ae, ar) = split(a);
    let (be, br) = split(b);
    grevlex(&ae, &be).then_with(|| grevlex(&ar, &br))
}

/// A polynomial as an ascending list of monomials under some order.
type GPoly = Vec<Mono>;

struct Order<'a>(&'a dyn Fn(Var) -> bool);

impl Order<'_> {
    fn sort(&self, mut ms: Vec<Mono>) -> GPoly {
        ms.sort_by(|a, b| cmp_with(a, b, self.0));
        let mut out: Vec<Mono> = Vec::with_capacity(ms.len());
        for m in ms {
            if out.last() == Some(&m) {
                out.pop();
            } else {
                out.push(m);
            }
        }
        out
    }

    fn add(&self, p: &GPoly, q: &GPoly) -> GPoly {
        self.sort(p.iter().chain(q).cloned().collect())
    }

    fn mul_mono(&self, p: &GPoly, m: &Mono) -> GPoly {
        self.sort(p.iter().map(|t| t.mul(m)).collect())
    }

    fn mul(&self, p: &GPoly, q: &GPoly) -> GPoly {
        self.sort(p.iter().flat_map(|a| q.iter().map(move |b| a.mul(b))).collect())
    }

    /// `leadReduce`.
    fn lead_reduce(&self, f: &GPoly, g: &GPoly) -> GPoly {
        let (Some(lf), Some(lg)) = (f.last(), g.last()) else { return f.clone() };
        if divides(lg, lf) {
            self.add(f, &self.mul_mono(g, &quot(lf, lg)))
        } else {
            f.clone()
        }
    }

    /// `mvd`: multivariate division, keeping irreducible leading terms.
    fn mvd(&self, f: &GPoly, basis: &[GPoly]) -> GPoly {
        let mut f = f.clone();
        let mut rest: Vec<Mono> = Vec::new();
        while !f.is_empty() {
            let f2 = basis.iter().fold(f.clone(), |acc, g| self.lead_reduce(&acc, g));
            if f2 == f {
                rest.push(f.pop().expect("nonzero"));
            } else {
                f = f2;
            }
        }
        self.sort(rest)
    }

    fn s_poly(&self, p: &GPoly, q: &GPoly) -> GPoly {
        let (m, n) = (p.last().unwrap(), q.last().unwrap());
        let lc = m.mul(n);
        self.add(&self.mul_mono(p, &quot(&lc, m)), &self.mul_mono(q, &quot(&lc, n)))
    }

    fn s_polys(&self, p: &GPoly, xs: &[GPoly]) -> Vec<GPoly> {
        let Some(lm) = p.last() else { return Vec::new() };
        // The implicit field equations x^2 = x: v * p for each v of LM(p).
        let mut out: Vec<GPoly> = lm.0.iter().map(|&v| self.mul_mono(p, &Mono::var(v))).collect();
        for q in xs {
            if let Some(lq) = q.last() {
                if !lm.0.iter().all(|v| !lq.contains(*v)) {
                    out.push(self.s_poly(p, q));
                }
            }
        }
        out
    }

    fn add_to_basis(&self, xs: &[GPoly], p: GPoly) -> Vec<GPoly> {
        let mut todo: std::collections::VecDeque<GPoly> = self.s_polys(&p, xs).into();
        let mut basis: Vec<GPoly> = xs.to_vec();
        basis.push(p);
        while let Some(s) = todo.pop_front() {
            let r = self.mvd(&s, &basis);
            if !r.is_empty() {
                todo.extend(self.s_polys(&r, &basis));
                basis.push(r);
            }
        }
        basis
    }

    fn reduce_basis(&self, basis: Vec<GPoly>) -> Vec<GPoly> {
        let mut done: Vec<GPoly> = Vec::new();
        for i in 0..basis.len() {
            let others: Vec<GPoly> = done.iter().chain(&basis[i + 1..]).cloned().collect();
            let r = self.mvd(&basis[i], &others);
            if !r.is_empty() {
                done.insert(0, r);
            }
        }
        done
    }

    /// `rbuchberger`: a reduced Gröbner basis.
    fn rbuchberger(&self, polys: Vec<GPoly>) -> Vec<GPoly> {
        polys.into_iter().fold(Vec::new(), |x, p| {
            let p = self.sort(p);
            self.reduce_basis(self.add_to_basis(&x, p))
        })
    }
}

fn divides(m: &Mono, n: &Mono) -> bool {
    m.0.iter().all(|v| n.contains(*v))
}

fn quot(n: &Mono, m: &Mono) -> Mono {
    Mono::new(n.0.iter().copied().filter(|v| !m.contains(*v)).collect())
}

const DEFAULT: Order<'static> = Order(&is_temp);

fn to_g(p: &BoolPoly) -> GPoly {
    DEFAULT.sort(p.monos())
}

fn from_g(p: &GPoly) -> BoolPoly {
    BoolPoly::from_list(p.clone())
}

fn rename(p: &GPoly, f: &dyn Fn(Var) -> Var) -> GPoly {
    DEFAULT.sort(
        p.iter()
            .map(|m| {
                let mut vs: SmallVec<[Var; 4]> = m.0.iter().map(|&v| f(v)).collect();
                vs.sort_unstable();
                vs.dedup();
                Mono::new(vs)
            })
            .collect(),
    )
}

fn has_temp(p: &GPoly) -> bool {
    p.iter().any(|m| m.0.iter().any(|&v| is_temp(v)))
}

/// `eliminateAll`: the basis polynomials free of path variables.
fn eliminate_all(ideal: Vec<GPoly>) -> Vec<GPoly> {
    DEFAULT.rbuchberger(ideal).into_iter().filter(|p| !has_temp(p)).collect()
}

/// `eliminateVars`: project `elim` out, under an order eliminating it.
fn eliminate_vars(elim: &[Var], ideal: Vec<GPoly>) -> Vec<GPoly> {
    let pred = |v: Var| elim.contains(&v);
    let order = Order(&pred);
    order
        .rbuchberger(ideal)
        .into_iter()
        .filter(|p| !p.iter().any(|m| m.0.iter().any(|&v| pred(v))))
        .map(|p| DEFAULT.sort(p))
        .collect()
}

fn ideal_plus(i: &[GPoly], j: &[GPoly]) -> Vec<GPoly> {
    DEFAULT.rbuchberger(i.iter().chain(j).cloned().collect())
}

/// `join`: the ideal product (the union of varieties).
fn join(i: &[GPoly], j: &[GPoly]) -> Vec<GPoly> {
    DEFAULT.rbuchberger(i.iter().flat_map(|p| j.iter().map(move |q| DEFAULT.mul(p, q))).collect())
}

fn compose(i: &[GPoly], j: &[GPoly]) -> Vec<GPoly> {
    // I's post-state and J's pre-state become the same path variables.
    let primes_to_temps = |v: Var| if is_prime(v) { TEMP | (v & !PRIME) } else { v };
    let inits_to_temps = |v: Var| if !is_prime(v) && !is_temp(v) { TEMP | v } else { v };
    let ishift: Vec<GPoly> = i.iter().map(|p| rename(p, &primes_to_temps)).collect();
    let jshift: Vec<GPoly> = j.iter().map(|p| rename(p, &inits_to_temps)).collect();
    eliminate_all(ideal_plus(&ishift, &jshift))
}

fn star(f: &[GPoly]) -> Vec<GPoly> {
    let mut i = f.to_vec();
    loop {
        let next = join(&i, &compose(&i, f));
        let as_set = |v: &[GPoly]| v.iter().cloned().collect::<std::collections::HashSet<_>>();
        if as_set(&i) == as_set(&next) {
            return i;
        }
        i = next;
    }
}

// ---------------------------------------------------------------------------
// The analysis of programs.

/// A program statement with each gate's location.
pub(crate) enum PStmt {
    Gate(usize, Gate),
    Reset(u32),
    Measure(u32),
    Seq(Vec<PStmt>),
    If(Box<PStmt>, Box<PStmt>),
    While(Box<PStmt>),
}

/// The angle `√2 · a` that marks a summarized phase as opaque.
fn scale_sqrt2(a: Angle) -> Angle {
    Angle::C(2f64.sqrt() * a.radians())
}

/// `mvdInPP`: each monomial of a phase polynomial, reduced modulo `ideal`
/// (with the coefficient distributed over the result).
fn mvd_in_pp(pp: &PhasePoly, ideal: &[GPoly], scale: bool, out: &mut PhasePoly) {
    let mut terms: Vec<(Mono, Angle)> = pp.terms();
    terms.sort_by(|a, b| a.0.cmp(&b.0));
    for (m, a) in terms {
        let r = DEFAULT.mvd(&vec![m], ideal);
        let a = if scale { scale_sqrt2(a) } else { a };
        distribute_into(a, &from_g(&r), out);
    }
}

impl Ctx {
    fn summary_ideal(&self) -> Vec<GPoly> {
        DEFAULT.rbuchberger(
            self.ket
                .iter()
                .enumerate()
                .map(|(i, p)| to_g(&p.add(&BoolPoly::var(PRIME | i as u32))))
                .collect(),
        )
    }

    fn take_terms_as_orphans(&mut self, other: &mut Ctx) {
        self.orphans.append(&mut other.orphans);
        for t in other.terms.iter().filter(|t| t.alive) {
            self.orphans.push((t.locs.clone(), t.angle));
        }
    }

    fn loop_summary(&mut self, mut body: Ctx) -> (PhasePoly, Vec<GPoly>) {
        self.take_terms_as_orphans(&mut body);
        let ideal = body.summary_ideal();
        let mut pp = PhasePoly::default();
        mvd_in_pp(&body.pp, &ideal, true, &mut pp);
        (pp, star(&eliminate_all(ideal)))
    }

    fn branch_summary(&mut self, mut a: Ctx, mut b: Ctx) -> (PhasePoly, Vec<GPoly>) {
        self.take_terms_as_orphans(&mut a);
        self.take_terms_as_orphans(&mut b);
        let (ia, ib) = (a.summary_ideal(), b.summary_ideal());
        let mut pp = PhasePoly::default();
        mvd_in_pp(&a.pp, &ia, true, &mut pp);
        mvd_in_pp(&b.pp, &ib, true, &mut pp);
        (pp, join(&eliminate_all(ia), &eliminate_all(ib)))
    }

    /// `fastForward`: apply a summary to the current state.
    fn fast_forward(&mut self, (poly, summary): (PhasePoly, Vec<GPoly>)) {
        let n = self.ket.len() as u32;
        let t = self.temps;
        let pre = |i: u32| TEMP | (i + n + t);
        let post = |i: u32| TEMP | (i + t);
        let trans = DEFAULT.rbuchberger(
            self.ket.iter().enumerate().map(|(i, p)| to_g(&p.add(&BoolPoly::var(pre(i as u32))))).collect(),
        );
        let shift = |v: Var| {
            if is_prime(v) {
                post(v & !PRIME)
            } else if !is_temp(v) {
                pre(v)
            } else {
                v
            }
        };
        let summary: Vec<GPoly> = summary.iter().map(|p| rename(p, &shift)).collect();
        let ideal = ideal_plus(&trans, &summary);
        let mut shifted = PhasePoly::default();
        for (m, a) in poly.terms() {
            let vs: SmallVec<[Var; 4]> = {
                let mut vs: SmallVec<[Var; 4]> = m.0.iter().map(|&v| shift(v)).collect();
                vs.sort_unstable();
                vs.dedup();
                vs
            };
            shifted.add_term(Mono::new(vs), a);
        }
        let mut reduced = PhasePoly::default();
        mvd_in_pp(&shifted, &ideal, false, &mut reduced);
        let evars: Vec<Var> = (0..n).map(pre).collect();
        let trans2 = eliminate_vars(&evars, ideal.clone());
        let ideal2 = ideal_plus(&ideal, &trans2);
        for i in 0..n {
            let k = DEFAULT.mvd(&vec![Mono::var(post(i))], &trans2);
            self.set_ket(i as usize, from_g(&k));
        }
        self.temps = t + n;
        for (m, a) in reduced.terms() {
            self.pp.add_term(m, a);
        }
        self.ideal = ideal2;
    }

    fn apply_stmt(&mut self, d: Option<usize>, s: &PStmt) {
        match s {
            PStmt::Gate(loc, g) => self.gate(g, *loc),
            PStmt::Reset(q) => self.set_ket(*q as usize, BoolPoly::default()),
            PStmt::Measure(_) => {}
            PStmt::Seq(xs) => xs.iter().for_each(|x| self.apply_stmt(d, x)),
            PStmt::If(a, b) => {
                let n = self.ket.len();
                let ca = Ctx::process_block(n, d, a);
                let cb = Ctx::process_block(n, d, b);
                let summary = self.branch_summary(ca, cb);
                self.fast_forward(summary);
            }
            PStmt::While(b) => {
                let body = Ctx::process_block(self.ket.len(), d, b);
                let summary = self.loop_summary(body);
                self.fast_forward(summary);
            }
        }
    }

    /// `processBlock`: run a block from a fresh state, then reduce.
    fn process_block(n: usize, d: Option<usize>, s: &PStmt) -> Ctx {
        let mut ctx = Ctx::new(n);
        ctx.run_block(d, s);
        ctx
    }

    fn run_block(&mut self, d: Option<usize>, s: &PStmt) {
        self.apply_stmt(d, s);
        self.reduce(Some(1));
        match d {
            Some(1) => {}
            Some(k) => self.reduce(Some(k)),
            None => self.reduce(None),
        }
        // `reduceTerms`: term keys modulo the accumulated ideal.
        if !self.ideal.is_empty() {
            let ideal = self.ideal.clone();
            let ids: Vec<usize> = (0..self.terms.len()).filter(|&i| self.terms[i].alive).collect();
            for id in ids {
                if !self.terms[id].alive {
                    continue;
                }
                let key = from_g(&DEFAULT.mvd(&to_g(&self.terms[id].key), &ideal));
                if key != self.terms[id].key {
                    self.unkey(id);
                    self.terms[id].key = BoolPoly::default();
                    let vars = key.vars();
                    self.rekey(id, key, &vars);
                }
            }
        }
    }
}

/// The new angle of each gate location, as Feynman's `stateAnalysispp`.
fn analyze(n: usize, d: Option<usize>, s: &PStmt) -> HashMap<usize, Angle> {
    let ctx = Ctx::process_block(n, d, s);
    let mut out = HashMap::default();
    let mut assign = |locs: &[(usize, bool)], angle: Option<Angle>, nonzero: bool| {
        let Some(angle) = angle else { return };
        let &(first, parity) = locs.iter().min().expect("a term has a location");
        let angle = if parity { angle.neg() } else { angle };
        for &(l, _) in locs {
            out.insert(l, if nonzero && l == first { angle } else { Angle::ZERO });
        }
    };
    for t in ctx.terms.iter().filter(|t| t.alive) {
        assign(&t.locs, t.angle, t.key != BoolPoly::default());
    }
    for (locs, angle) in &ctx.orphans {
        assign(locs, *angle, true);
    }
    out
}

/// Re-synthesize every phase gate the analysis assigned an angle.
fn rewrite(s: &PStmt, angles: &HashMap<usize, Angle>) -> Vec<crate::program::Stmt> {
    use crate::program::Stmt;
    match s {
        PStmt::Gate(loc, g) => match (angles.get(loc), g) {
            (
                Some(&a),
                Gate::t(q) | Gate::tdg(q) | Gate::s(q) | Gate::sdg(q) | Gate::z(q) | Gate::rz(_, q),
            ) => {
                let mut gates = Vec::new();
                synthesize(*q, a, &mut gates);
                gates.into_iter().map(Stmt::Gate).collect()
            }
            _ => vec![Stmt::Gate(g.clone())],
        },
        PStmt::Reset(q) => vec![Stmt::Reset(*q)],
        PStmt::Measure(q) => vec![Stmt::Measure(*q)],
        PStmt::Seq(xs) => xs.iter().flat_map(|x| rewrite(x, angles)).collect(),
        PStmt::If(a, b) => vec![Stmt::If(Box::new(Stmt::Seq(rewrite(a, angles))), Box::new(Stmt::Seq(rewrite(b, angles))))],
        PStmt::While(b) => vec![Stmt::While(Box::new(Stmt::Seq(rewrite(b, angles))))],
    }
}

/// Feynman's `stateFoldpp`: the program, with its phase gates re-synthesized.
pub(crate) fn fold(n: usize, d: Option<usize>, s: &PStmt) -> Vec<crate::program::Stmt> {
    rewrite(s, &analyze(n, d, s))
}
