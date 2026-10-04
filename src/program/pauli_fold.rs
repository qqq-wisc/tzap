//! PhaseFoldPauli over programs: the abstraction of programs and the phase
//! folding algorithm of the paper, extended to branches, loops, measurement
//! and reset.
//!
//! An abstract relation is a set of constraints `x_q = P`, `z_q = Q` (each
//! possibly unknown, written ⋆) and a set of axes, or ⊤. Cliffords update the
//! constraints by the parallel assignments of Figure 1; a rotation, a
//! measurement or a reset adds axes; a branch joins its two sides; and a loop
//! is the fixpoint σ_{i+1} = σ_i ⊔ (α(body) ∘ σ_i) from the identity.
//!
//! Folding runs Algorithm 1 on each block, left to right. A branch or loop in
//! the block is a segment: its relation updates the constraints and its axes,
//! pulled back to the block's input, block the folds that cross it. The body of
//! each branch and loop is folded on its own, which merges rotations within one
//! run of the body.

use super::{Program, Stmt};
use crate::circuit::{Gate, Qubit};
use std::f64::consts::FRAC_PI_4;

// ---------------------------------------------------------------------------
// Pauli strings: i^phase · Π_q X_q^{x_q} Z_q^{z_q}.

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Pauli {
    x: Vec<u64>,
    z: Vec<u64>,
    /// The power of i, mod 4.
    phase: u8,
}

impl Pauli {
    fn identity(n: usize) -> Pauli {
        let w = n.div_ceil(64).max(1);
        Pauli { x: vec![0; w], z: vec![0; w], phase: 0 }
    }

    fn single(n: usize, q: usize, x: bool, z: bool) -> Pauli {
        let mut p = Pauli::identity(n);
        if x {
            p.x[q / 64] |= 1 << (q % 64);
        }
        if z {
            p.z[q / 64] |= 1 << (q % 64);
        }
        p
    }

    fn bit(v: &[u64], q: usize) -> bool {
        v[q / 64] >> (q % 64) & 1 == 1
    }

    /// `self · other`.
    fn mul(&self, other: &Pauli) -> Pauli {
        // Moving each Z of `self` past an X of `other` on the same qubit
        // contributes a factor -1 (i^2).
        let mut flips = 0u32;
        for (sz, ox) in self.z.iter().zip(&other.x) {
            flips += (sz & ox).count_ones();
        }
        Pauli {
            x: self.x.iter().zip(&other.x).map(|(a, b)| a ^ b).collect(),
            z: self.z.iter().zip(&other.z).map(|(a, b)| a ^ b).collect(),
            phase: ((self.phase as u32 + other.phase as u32 + 2 * flips) % 4) as u8,
        }
    }

    fn scale(mut self, k: u8) -> Pauli {
        self.phase = (self.phase + k) % 4;
        self
    }

    fn anticommutes(&self, other: &Pauli) -> bool {
        let mut n = 0u32;
        for i in 0..self.x.len() {
            n += (self.x[i] & other.z[i]).count_ones() + (self.z[i] & other.x[i]).count_ones();
        }
        n % 2 == 1
    }

    /// Equal up to sign (the phase is 0 or 2 apart).
    fn same_axis(&self, other: &Pauli) -> Option<bool> {
        (self.x == other.x && self.z == other.z && (self.phase + 4 - other.phase) % 2 == 0)
            .then(|| self.phase != other.phase)
    }
}

// ---------------------------------------------------------------------------
// Abstract relations.

/// An axis: a known Pauli string, or one whose value is unknown (it blocks
/// every fold across it).
#[derive(Clone, Debug, PartialEq)]
enum Axis {
    Known(Pauli),
    Unknown,
}

impl Axis {
    fn blocks(&self, p: &Pauli) -> bool {
        match self {
            Axis::Known(a) => a.anticommutes(p),
            Axis::Unknown => true,
        }
    }
}

/// `(φ, axes)`: `cons[2q]` is `x_q`, `cons[2q+1]` is `z_q` (`None`: ⋆).
#[derive(Clone, Debug, PartialEq)]
struct Rel {
    cons: Vec<Option<Pauli>>,
    axes: Vec<Axis>,
}

impl Rel {
    fn identity(n: usize) -> Rel {
        Rel {
            cons: (0..2 * n).map(|v| Some(Pauli::single(n, v / 2, v % 2 == 0, v % 2 == 1))).collect(),
            axes: Vec::new(),
        }
    }

    fn n(&self) -> usize {
        self.cons.len() / 2
    }

    /// `wp_φ(P)`, or `None` if it needs an unknown constraint.
    fn wp(&self, p: &Pauli) -> Option<Pauli> {
        let n = self.n();
        let mut out = Pauli::identity(n).scale(p.phase);
        for q in 0..n {
            if Pauli::bit(&p.x, q) {
                out = out.mul(self.cons[2 * q].as_ref()?);
            }
            if Pauli::bit(&p.z, q) {
                out = out.mul(self.cons[2 * q + 1].as_ref()?);
            }
        }
        Some(out)
    }

    fn wp_axis(&self, a: &Axis) -> Axis {
        match a {
            Axis::Known(p) => self.wp(p).map_or(Axis::Unknown, Axis::Known),
            Axis::Unknown => Axis::Unknown,
        }
    }

    fn push_axis(&mut self, a: Axis) {
        if !self.axes.contains(&a) {
            self.axes.push(a);
        }
    }

    /// `other ∘ self`: `self` first.
    fn then(&self, other: &Rel) -> Rel {
        let mut out = Rel {
            cons: other.cons.iter().map(|c| c.as_ref().and_then(|p| self.wp(p))).collect(),
            axes: self.axes.clone(),
        };
        for a in &other.axes {
            out.push_axis(self.wp_axis(a));
        }
        out
    }

    fn join(&self, other: &Rel) -> Rel {
        let mut out = Rel {
            cons: self.cons.iter().zip(&other.cons).map(|(a, b)| if a == b { a.clone() } else { None }).collect(),
            axes: self.axes.clone(),
        };
        for a in &other.axes {
            out.push_axis(a.clone());
        }
        out
    }

    /// Apply a Clifford gate (Figure 1), as `α(G) ∘ self`.
    fn clifford(&mut self, g: &Gate) {
        let mul = |a: &Option<Pauli>, b: &Option<Pauli>| match (a, b) {
            (Some(a), Some(b)) => Some(a.mul(b)),
            _ => None,
        };
        let neg = |a: &Option<Pauli>| a.clone().map(|p| p.scale(2));
        match *g {
            Gate::h(q) => {
                let (x, z) = x_z(q);
                self.cons.swap(x, z);
            }
            Gate::s(q) => {
                // S†XS = -Y = -i X Z.
                let (x, z) = x_z(q);
                self.cons[x] = mul(&self.cons[x], &self.cons[z]).map(|p| p.scale(3));
            }
            Gate::sdg(q) => {
                // S X S† = Y = i X Z.
                let (x, z) = x_z(q);
                self.cons[x] = mul(&self.cons[x], &self.cons[z]).map(|p| p.scale(1));
            }
            Gate::x(q) => {
                let (_, z) = x_z(q);
                self.cons[z] = neg(&self.cons[z]);
            }
            Gate::z(q) => {
                let (x, _) = x_z(q);
                self.cons[x] = neg(&self.cons[x]);
            }
            Gate::cnot { control, target } => {
                let (xc, zc) = x_z(control);
                let (xt, zt) = x_z(target);
                self.cons[xc] = mul(&self.cons[xc], &self.cons[xt]);
                self.cons[zt] = mul(&self.cons[zc], &self.cons[zt]);
            }
            Gate::cz { control, target } => {
                let (xc, zc) = x_z(control);
                let (xt, zt) = x_z(target);
                self.cons[xc] = mul(&self.cons[xc], &self.cons[zt]);
                self.cons[xt] = mul(&self.cons[xt], &self.cons[zc]);
            }
            _ => unreachable!("not a Clifford gate: {g:?}"),
        }
    }
}

fn x_z(q: Qubit) -> (usize, usize) {
    (2 * q as usize, 2 * q as usize + 1)
}

// ---------------------------------------------------------------------------
// Programs with mergeable rotations.

/// An angle `eighths · π/4 + residual`.
#[derive(Clone, Copy, Debug)]
struct Angle {
    eighths: i64,
    residual: f64,
}

impl Angle {
    fn add(self, other: Angle, sign: i64) -> Angle {
        Angle {
            eighths: (self.eighths + sign * other.eighths).rem_euclid(8),
            residual: self.residual + sign as f64 * other.residual,
        }
    }

    fn emit(self, q: Qubit, out: &mut Vec<Stmt>) {
        if self.residual != 0.0 {
            out.push(Stmt::Gate(Gate::rz(self.eighths as f64 * FRAC_PI_4 + self.residual, q)));
            return;
        }
        let gates: &[fn(Qubit) -> Gate] = match self.eighths {
            0 => &[],
            1 => &[Gate::t],
            2 => &[Gate::s],
            3 => &[Gate::s, Gate::t],
            4 => &[Gate::z],
            5 => &[Gate::z, Gate::t],
            6 => &[Gate::sdg],
            _ => &[Gate::tdg],
        };
        out.extend(gates.iter().map(|g| Stmt::Gate(g(q))));
    }
}

/// A statement with its rotations explicit, so folds can change their angles.
enum Node {
    Clifford(Gate),
    Rot(Qubit, Angle),
    Reset(Qubit),
    Measure(Qubit),
    If(Vec<Node>, Vec<Node>),
    While(Vec<Node>),
}

fn rotation_angle(g: &Gate) -> Option<(Qubit, Angle)> {
    match *g {
        Gate::t(q) => Some((q, Angle { eighths: 1, residual: 0.0 })),
        Gate::tdg(q) => Some((q, Angle { eighths: 7, residual: 0.0 })),
        Gate::rz(theta, q) => Some((q, Angle { eighths: 0, residual: theta })),
        _ => None,
    }
}

fn to_nodes(s: &Stmt, out: &mut Vec<Node>) {
    match s {
        Stmt::Gate(g) => out.push(match rotation_angle(g) {
            Some((q, a)) => Node::Rot(q, a),
            None => Node::Clifford(g.clone()),
        }),
        Stmt::Reset(q) => out.push(Node::Reset(*q)),
        Stmt::Measure(q) => out.push(Node::Measure(*q)),
        Stmt::Seq(xs) => xs.iter().for_each(|x| to_nodes(x, out)),
        Stmt::If(a, b) => {
            let (mut na, mut nb) = (Vec::new(), Vec::new());
            to_nodes(a, &mut na);
            to_nodes(b, &mut nb);
            out.push(Node::If(na, nb));
        }
        Stmt::While(b) => {
            let mut nb = Vec::new();
            to_nodes(b, &mut nb);
            out.push(Node::While(nb));
        }
    }
}

fn to_stmt(nodes: &[Node]) -> Stmt {
    let mut out = Vec::new();
    for node in nodes {
        match node {
            Node::Clifford(g) => out.push(Stmt::Gate(g.clone())),
            Node::Rot(q, a) => a.emit(*q, &mut out),
            Node::Reset(q) => out.push(Stmt::Reset(*q)),
            Node::Measure(q) => out.push(Stmt::Measure(*q)),
            Node::If(a, b) => out.push(Stmt::If(Box::new(to_stmt(a)), Box::new(to_stmt(b)))),
            Node::While(b) => out.push(Stmt::While(Box::new(to_stmt(b)))),
        }
    }
    Stmt::Seq(out)
}

/// How precisely to analyze: the bound on disjuncts (1 is the plain domain
/// with joins), and whether to track eigenstate facts from resets and
/// measurements.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub disjuncts: usize,
    pub zero_facts: bool,
}

impl Options {
    pub const PLAIN: Options = Options { disjuncts: 1, zero_facts: false };
}

impl Rel {
    /// Equal as relations: the same constraints and the same set of axes.
    fn same(&self, other: &Rel) -> bool {
        self.cons == other.cons
            && self.axes.len() == other.axes.len()
            && self.axes.iter().all(|a| other.axes.contains(a))
    }
}

/// Add `r` to a set of relations, unless an equal one is there.
fn insert(set: &mut Vec<Rel>, r: Rel) -> bool {
    if set.iter().any(|s| s.same(&r)) {
        return false;
    }
    set.push(r);
    true
}

/// The join of a set of relations.
fn join_all(set: &[Rel]) -> Rel {
    set[1..].iter().fold(set[0].clone(), |acc, r| acc.join(r))
}

/// Keep a set within the bound by joining it into one relation.
fn bound(set: Vec<Rel>, k: usize) -> Vec<Rel> {
    if set.len() > k { vec![join_all(&set)] } else { set }
}

/// The disjunctive abstraction of a block: a set of relations, one of which
/// holds on each path.
fn alpha_set(n: usize, nodes: &[Node], k: usize) -> Vec<Rel> {
    let mut set = vec![Rel::identity(n)];
    for node in nodes {
        set = bound(step_set(n, set, node, k), k);
    }
    set
}

/// `α(node) ∘ set`, disjunct by disjunct.
fn step_set(n: usize, set: Vec<Rel>, node: &Node, k: usize) -> Vec<Rel> {
    let after: Vec<Rel> = match node {
        Node::If(a, b) => {
            let mut sides = alpha_set(n, a, k);
            for r in alpha_set(n, b, k) {
                insert(&mut sides, r);
            }
            compose_sets(&set, &sides)
        }
        Node::While(b) => compose_sets(&set, &loop_set(n, &alpha_set(n, b, k), k)),
        Node::Clifford(g) => set
            .into_iter()
            .map(|mut r| {
                r.clifford(g);
                r
            })
            .collect(),
        _ => {
            let own = own_rel(n, node);
            set.iter().map(|r| r.then(&own)).collect()
        }
    };
    let mut out = Vec::new();
    for r in after {
        insert(&mut out, r);
    }
    out
}

fn compose_sets(set: &[Rel], segment: &[Rel]) -> Vec<Rel> {
    set.iter().flat_map(|r| segment.iter().map(move |s| r.then(s))).collect()
}

/// The relation of a single gate, rotation, measurement or reset.
fn own_rel(n: usize, node: &Node) -> Rel {
    let mut rel = Rel::identity(n);
    match node {
        Node::Clifford(g) => rel.clifford(g),
        Node::Rot(q, _) | Node::Measure(q) => rel.push_axis(Axis::Known(Pauli::single(n, *q as usize, false, true))),
        Node::Reset(q) => {
            rel.push_axis(Axis::Known(Pauli::single(n, *q as usize, true, false)));
            rel.push_axis(Axis::Known(Pauli::single(n, *q as usize, false, true)));
        }
        Node::If(..) | Node::While(..) => unreachable!("compound"),
    }
    rel
}

/// `α(while body)`: the least set containing the identity and closed under
/// running the body once more; past the bound, the join's fixpoint.
fn loop_set(n: usize, body: &[Rel], k: usize) -> Vec<Rel> {
    let mut set = vec![Rel::identity(n)];
    loop {
        let mut changed = false;
        for s in compose_sets(&set, body) {
            changed |= insert(&mut set, s);
        }
        if set.len() > k {
            let b = join_all(body);
            let mut r = join_all(&set);
            loop {
                let next = r.join(&r.then(&b));
                if next.same(&r) {
                    return vec![r];
                }
                r = next;
            }
        }
        if !changed {
            return set;
        }
    }
}

// ---------------------------------------------------------------------------
// Eigenstate facts: an unsigned stabilizer group of the current state.

/// Pauli strings `P` (up to sign) with `P|ψ⟩ = ±|ψ⟩` on every path, kept as
/// a basis of their span; `x` and `z` are over the current qubits.
#[derive(Clone, Debug, PartialEq)]
struct Facts {
    basis: Vec<(Vec<u64>, Vec<u64>)>,
    words: usize,
}

impl Facts {
    fn new(n: usize) -> Facts {
        Facts { basis: Vec::new(), words: n.div_ceil(64).max(1) }
    }

    fn reduce(&self, mut x: Vec<u64>, mut z: Vec<u64>) -> (Vec<u64>, Vec<u64>) {
        // Each basis element has a distinct pivot bit; eliminate it.
        for (bx, bz) in &self.basis {
            if let Some(p) = pivot(bx, bz) {
                if get(&x, &z, p, self.words) {
                    xor(&mut x, bx);
                    xor(&mut z, bz);
                }
            }
        }
        (x, z)
    }

    fn contains(&self, x: &[u64], z: &[u64]) -> bool {
        let (x, z) = self.reduce(x.to_vec(), z.to_vec());
        x.iter().all(|w| *w == 0) && z.iter().all(|w| *w == 0)
    }

    fn add(&mut self, x: Vec<u64>, z: Vec<u64>) {
        let (x, z) = self.reduce(x, z);
        if let Some(p) = pivot(&x, &z) {
            // Keep pivots unique: clear p from the existing elements.
            for (bx, bz) in &mut self.basis {
                if get(bx, bz, p, self.words) {
                    xor(bx, &x);
                    xor(bz, &z);
                }
            }
            self.basis.push((x, z));
        }
    }

    /// Keep the subgroup commuting with `(ax, az)`.
    fn commute_with(&mut self, ax: &[u64], az: &[u64]) {
        let anti = |x: &[u64], z: &[u64]| {
            let mut c = 0u32;
            for i in 0..x.len() {
                c += (x[i] & az[i]).count_ones() + (z[i] & ax[i]).count_ones();
            }
            c % 2 == 1
        };
        let elems: Vec<_> = std::mem::take(&mut self.basis);
        let mut first: Option<(Vec<u64>, Vec<u64>)> = None;
        let mut keep = Vec::new();
        for (x, z) in elems {
            if anti(&x, &z) {
                match &first {
                    None => first = Some((x, z)),
                    Some((fx, fz)) => {
                        let (mut x, mut z) = (x, z);
                        xor(&mut x, fx);
                        xor(&mut z, fz);
                        keep.push((x, z));
                    }
                }
            } else {
                keep.push((x, z));
            }
        }
        for (x, z) in keep {
            self.add(x, z);
        }
    }

    /// The facts in both (a subset of the intersection of the spans).
    fn meet(&self, other: &Facts) -> Facts {
        let mut out = Facts { basis: Vec::new(), words: self.words };
        for (x, z) in &self.basis {
            if other.contains(x, z) {
                out.add(x.clone(), z.clone());
            }
        }
        out
    }

    fn single(&self, q: usize, x: bool, z: bool) -> (Vec<u64>, Vec<u64>) {
        let p = Pauli::single(self.words * 64, q, x, z);
        (p.x[..self.words].to_vec(), p.z[..self.words].to_vec())
    }

    /// Conjugate every fact by a Clifford gate (signs are not tracked).
    fn clifford(&mut self, g: &Gate) {
        let elems: Vec<_> = std::mem::take(&mut self.basis);
        for (mut x, mut z) in elems {
            let bit = |v: &[u64], q: u32| v[q as usize / 64] >> (q % 64) & 1 == 1;
            let flip = |v: &mut [u64], q: u32| v[q as usize / 64] ^= 1 << (q % 64);
            match *g {
                Gate::h(q) => {
                    if bit(&x, q) != bit(&z, q) {
                        flip(&mut x, q);
                        flip(&mut z, q);
                    }
                }
                Gate::s(q) | Gate::sdg(q) => {
                    if bit(&x, q) {
                        flip(&mut z, q);
                    }
                }
                Gate::x(_) | Gate::z(_) => {}
                Gate::cnot { control, target } => {
                    if bit(&x, control) {
                        flip(&mut x, target);
                    }
                    if bit(&z, target) {
                        flip(&mut z, control);
                    }
                }
                Gate::cz { control, target } => {
                    let (xc, xt) = (bit(&x, control), bit(&x, target));
                    if xc {
                        flip(&mut z, target);
                    }
                    if xt {
                        flip(&mut z, control);
                    }
                }
                _ => {}
            }
            self.add(x, z);
        }
    }

    /// The facts after `node`, given the facts before it.
    fn step(&self, node: &Node) -> Facts {
        let mut f = self.clone();
        match node {
            Node::Clifford(g) => f.clifford(g),
            Node::Rot(q, _) => {
                let (x, z) = f.single(*q as usize, false, true);
                f.commute_with(&x, &z);
            }
            Node::Measure(q) | Node::Reset(q) => {
                // Afterwards the qubit is a Z eigenstate (|0⟩ after a reset).
                let (x, z) = f.single(*q as usize, false, true);
                f.commute_with(&x, &z);
                f.add(x, z);
            }
            Node::If(a, b) => f = f.block(a).meet(&f.block(b)),
            Node::While(b) => f = f.loop_invariant(b),
        }
        f
    }

    fn block(&self, nodes: &[Node]) -> Facts {
        nodes.iter().fold(self.clone(), |f, node| f.step(node))
    }

    /// The facts at a loop's head: what holds on entry and after every run.
    fn loop_invariant(&self, body: &[Node]) -> Facts {
        let mut f = self.clone();
        loop {
            let next = f.meet(&f.block(body));
            if next.basis.len() == f.basis.len() {
                return f;
            }
            f = next;
        }
    }
}

fn pivot(x: &[u64], z: &[u64]) -> Option<usize> {
    let w = x.len();
    for i in 0..w {
        if x[i] != 0 {
            return Some(i * 64 + x[i].trailing_zeros() as usize);
        }
    }
    for i in 0..w {
        if z[i] != 0 {
            return Some(w * 64 + i * 64 + z[i].trailing_zeros() as usize);
        }
    }
    None
}

fn get(x: &[u64], z: &[u64], p: usize, w: usize) -> bool {
    let (v, p) = if p < w * 64 { (x, p) } else { (z, p - w * 64) };
    v[p / 64] >> (p % 64) & 1 == 1
}

fn xor(a: &mut [u64], b: &[u64]) {
    for (x, y) in a.iter_mut().zip(b) {
        *x ^= y;
    }
}

// ---------------------------------------------------------------------------
// Folding.

/// One disjunct of a block's state: its relation, the axes recorded so far
/// (each with the number of rotations recorded before it), and each recorded
/// rotation's axis (`None`: unknown, so it cannot be folded into).
#[derive(Clone)]
struct Disjunct {
    rel: Rel,
    events: Vec<(usize, Axis)>,
    rot_axes: Vec<Option<Pauli>>,
}

impl Disjunct {
    fn same(&self, other: &Disjunct) -> bool {
        self.rel.same(&other.rel) && self.events == other.events && self.rot_axes == other.rot_axes
    }

    /// Whether rotation `j`'s axis commutes with every axis recorded after it.
    fn free_after(&self, j: usize, a: &Pauli) -> bool {
        self.events.iter().all(|(e, b)| *e <= j || !b.blocks(a))
    }
}

fn collapse(ds: Vec<Disjunct>, k: usize) -> Vec<Disjunct> {
    let mut out: Vec<Disjunct> = Vec::new();
    for d in ds {
        if !out.iter().any(|o| o.same(&d)) {
            out.push(d);
        }
    }
    if out.len() <= k {
        return out;
    }
    let rels: Vec<Rel> = out.iter().map(|d| d.rel.clone()).collect();
    let mut events: Vec<(usize, Axis)> = Vec::new();
    for d in &out {
        for e in &d.events {
            if !events.contains(e) {
                events.push(e.clone());
            }
        }
    }
    let rot_axes = (0..out[0].rot_axes.len())
        .map(|j| {
            let a = &out[0].rot_axes[j];
            out.iter().all(|d| &d.rot_axes[j] == a).then(|| a.clone()).flatten()
        })
        .collect();
    vec![Disjunct { rel: join_all(&rels), events, rot_axes }]
}

/// Algorithm 1 on one block, recursing into branches and loops. Returns the
/// eigenstate facts after the block.
fn fold_block(n: usize, nodes: &mut [Node], opts: Options, facts_in: &Facts) -> Facts {
    let k = opts.disjuncts;
    let mut ds = vec![Disjunct { rel: Rel::identity(n), events: Vec::new(), rot_axes: Vec::new() }];
    // The node index of each recorded rotation.
    let mut rots: Vec<usize> = Vec::new();
    let mut facts = facts_in.clone();
    for i in 0..nodes.len() {
        match &mut nodes[i] {
            Node::If(a, b) => {
                let fa = fold_block(n, a, opts, &facts);
                let fb = fold_block(n, b, opts, &facts);
                if opts.zero_facts {
                    facts = fa.meet(&fb);
                }
            }
            Node::While(b) => {
                let head = facts.loop_invariant(b);
                fold_block(n, b, opts, &head);
                if opts.zero_facts {
                    facts = head;
                }
            }
            _ => {}
        }
        if let Node::Rot(q, angle) = nodes[i] {
            let (zx, zz) = facts.single(q as usize, false, true);
            if opts.zero_facts && facts.contains(&zx, &zz) {
                // Acting on an eigenstate: a global phase on every path.
                nodes[i] = Node::Rot(q, Angle { eighths: 0, residual: 0.0 });
                continue;
            }
            if opts.zero_facts {
                facts.commute_with(&zx, &zz);
            }
            let zq = Pauli::single(n, q as usize, false, true);
            let axes: Vec<Option<Pauli>> = ds.iter().map(|d| d.rel.wp(&zq)).collect();
            // The latest rotation it folds into in every disjunct, with one sign.
            let hit = (0..rots.len()).rev().find_map(|j| {
                let mut sign = None;
                for (d, p) in ds.iter().zip(&axes) {
                    let (p, a) = (p.as_ref()?, d.rot_axes[j].as_ref()?);
                    let s = p.same_axis(a)?;
                    if sign.is_some_and(|t| t != s) || !d.free_after(j, a) {
                        return None;
                    }
                    sign = Some(s);
                }
                sign.map(|s| (j, s))
            });
            if let Some((j, negated)) = hit {
                if let Node::Rot(_, prior) = &mut nodes[rots[j]] {
                    *prior = prior.add(angle, if negated { -1 } else { 1 });
                }
                nodes[i] = Node::Rot(q, Angle { eighths: 0, residual: 0.0 });
                continue;
            }
            rots.push(i);
            for (d, p) in ds.iter_mut().zip(axes) {
                d.events.push((rots.len(), p.clone().map_or(Axis::Unknown, Axis::Known)));
                d.rot_axes.push(p);
            }
            continue;
        }
        if opts.zero_facts {
            facts = match &nodes[i] {
                Node::If(..) | Node::While(..) => facts,
                node => facts.step(node),
            };
        }
        if let Node::Clifford(g) = &nodes[i] {
            for d in &mut ds {
                d.rel.clifford(g);
            }
            continue;
        }
        // The statement's relations; their axes, pulled back, block the folds
        // across it.
        let segment: Vec<Rel> = match &nodes[i] {
            Node::If(..) | Node::While(..) => alpha_set(n, std::slice::from_ref(&nodes[i]), k),
            node => vec![own_rel(n, node)],
        };
        let epoch = rots.len();
        let next: Vec<Disjunct> = ds
            .iter()
            .flat_map(|d| {
                segment.iter().map(move |s| {
                    let mut e = d.events.clone();
                    for a in &s.axes {
                        e.push((epoch, d.rel.wp_axis(a)));
                    }
                    Disjunct { rel: d.rel.then(s), events: e, rot_axes: d.rot_axes.clone() }
                })
            })
            .collect();
        ds = collapse(next, k);
    }
    facts
}

/// PhaseFoldPauli over a program, with the plain domain.
pub fn fold(prog: &Program) -> Program {
    fold_with(prog, Options::PLAIN)
}

/// PhaseFoldPauli over a program.
pub fn fold_with(prog: &Program, opts: Options) -> Program {
    let mut nodes = Vec::new();
    to_nodes(&prog.body, &mut nodes);
    fold_block(prog.num_qubits, &mut nodes, opts, &Facts::new(prog.num_qubits));
    Program { num_qubits: prog.num_qubits, body: super::parse::flatten(to_stmt(&nodes)) }
}

#[cfg(test)]
mod tests {
    use super::super::testing::*;
    use super::*;

    const MODES: [Options; 4] = [
        Options { disjuncts: 1, zero_facts: false },
        Options { disjuncts: 100, zero_facts: false },
        Options { disjuncts: 1, zero_facts: true },
        Options { disjuncts: 100, zero_facts: true },
    ];

    /// Every path of the folded program has the operator of the matching path
    /// of the original, up to a phase: loops up to twice, every branch and
    /// outcome.
    #[test]
    fn random_programs_keep_every_path() {
        let mut rng = 0x9a11_5eed_u64;
        let mut changed = 0;
        for case in 0..1500 {
            let n = 2 + (case % 2) as u32;
            let body = random_stmt(&mut rng, n, 5 + case % 5, 2);
            let prog = Program { num_qubits: n as usize, body: super::super::parse::flatten(body) };
            if count(&prog.body, 2) > 2000 {
                continue;
            }
            for opts in MODES {
                let out = fold_with(&prog, opts);
                changed += (out.t_count() < prog.t_count()) as usize;
                let (pa, pb) = (paths(&prog.body, 2), paths(&out.body, 2));
                assert_eq!(pa.len(), pb.len(), "case {case}: path structure changed");
                for (i, (a, b)) in pa.iter().zip(&pb).enumerate().take(400) {
                    let (ma, mb) = (operator(n as usize, a), operator(n as usize, b));
                    assert!(equiv(&ma, &mb, 1e-9), "case {case}, {opts:?}, path {i}:\n{}\n=>\n{}", prog.to_qasm3(), out.to_qasm3());
                }
            }
        }
        assert!(changed > 30, "too few folds: {changed}");
    }

    /// `t q; C; S; C⁻¹; t q` with random Cliffords `C` and a random
    /// measurement, reset, branch or loop `S`: the two rotations share an
    /// axis, so whether they fold is decided by `S` alone.
    #[test]
    fn sandwiches_keep_every_path() {
        let mut rng = 0x5a4d_0001_u64;
        let mut next = |m: u64| -> u64 {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng % m
        };
        let (mut folded, mut kept) = (0, 0);
        for case in 0..2000 {
            let n = 2 + (case % 2) as u32;
            let q = next(n as u64) as u32;
            let mut cliffords = Vec::new();
            for _ in 0..1 + next(4) {
                let a = next(n as u64) as u32;
                let b = (a + 1 + next(n as u64 - 1) as u32) % n;
                cliffords.push(match next(4) {
                    0 => Gate::h(a),
                    1 => Gate::s(a),
                    2 => Gate::x(a),
                    _ => Gate::cnot { control: a, target: b },
                });
            }
            let inverse: Vec<Gate> = cliffords
                .iter()
                .rev()
                .map(|g| match *g {
                    Gate::s(a) => Gate::sdg(a),
                    ref g => g.clone(),
                })
                .collect();
            let b = next(n as u64) as u32;
            let inner = random_stmt(&mut (next(1 << 30) | 1), n, 1 + next(2) as usize, 0);
            let middle = match next(4) {
                0 => Stmt::Measure(b),
                1 => Stmt::Reset(b),
                2 => Stmt::If(Box::new(inner), Box::new(Stmt::Seq(Vec::new()))),
                _ => Stmt::While(Box::new(inner)),
            };
            let mut body = vec![Stmt::Gate(Gate::t(q))];
            body.extend(cliffords.into_iter().map(Stmt::Gate));
            body.push(middle);
            body.extend(inverse.into_iter().map(Stmt::Gate));
            body.push(Stmt::Gate(Gate::t(q)));
            let prog = Program { num_qubits: n as usize, body: super::super::parse::flatten(Stmt::Seq(body)) };
            if count(&prog.body, 2) > 2000 {
                continue;
            }
            for opts in MODES {
                let out = fold_with(&prog, opts);
                if opts.disjuncts == 1 && !opts.zero_facts {
                    if out.t_count() < prog.t_count() { folded += 1 } else { kept += 1 }
                }
                for (a, b) in paths(&prog.body, 2).iter().zip(&paths(&out.body, 2)) {
                    assert!(
                        equiv(&operator(n as usize, a), &operator(n as usize, b), 1e-9),
                        "case {case}, {opts:?}:\n{}\n=>\n{}",
                        prog.to_qasm3(),
                        out.to_qasm3()
                    );
                }
            }
        }
        assert!(folded > 200 && kept > 200, "folded {folded}, kept {kept}");
    }

    /// Folds that would cross an anticommuting measurement, reset, branch or
    /// loop are refused, and the path check catches them.
    #[test]
    fn blockers_refuse_folds() {
        use super::super::parse;
        for src in [
            "qubit a; t a; h a; measure a; h a; t a;",
            "qubit a; t a; h a; reset a; h a; t a;",
            "qubit a; qubit b; t a; while (true) { h a; } t a;",
            "qubit a; qubit b; t a; if (true) { h a; t a; h a; } t a;",
        ] {
            let prog = parse(src).unwrap();
            let out = fold(&prog);
            assert_eq!(out.t_count(), prog.t_count(), "{src}");
            let (pa, pb) = (paths(&prog.body, 2), paths(&out.body, 2));
            for (a, b) in pa.iter().zip(&pb) {
                assert!(equiv(&operator(prog.num_qubits, a), &operator(prog.num_qubits, b), 1e-9));
            }
        }
    }

    /// The path check itself rejects an unsound fold.
    #[test]
    fn path_check_catches_a_bad_fold() {
        use super::super::parse;
        let bad = parse("qubit a; h a; measure a; h a; s a;").unwrap();
        let orig = parse("qubit a; t a; h a; measure a; h a; t a;").unwrap();
        let (pa, pb) = (paths(&orig.body, 2), paths(&bad.body, 2));
        assert!(pa.iter().zip(&pb).any(|(a, b)| !equiv(&operator(1, a), &operator(1, b), 1e-9)));
    }
}
