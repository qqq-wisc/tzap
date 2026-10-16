use super::*;
use num_bigint::BigInt;
use num_traits::One;
use rustc_hash::FxHashMap;
use std::sync::Weak;

#[derive(Clone, Copy, Debug, Default)]
pub struct InternerStats {
    pub hits: usize,
    pub misses: usize,
    pub pruned: usize,
    pub retained_entries: usize,
}

/// Local weak interning. Returned angles own their nodes and outlive this
/// context. Every public operation has a fresh, bounded work budget.
pub struct AngleContext {
    limits: Limits,
    remaining_run_work: Option<usize>,
    pub(super) nodes: FxHashMap<u64, Vec<Weak<Node>>>,
    pub(super) stats: InternerStats,
}
impl Default for AngleContext {
    fn default() -> Self {
        Self::with_limits(Limits::default())
    }
}
impl AngleContext {
    pub fn with_limits(limits: Limits) -> Self {
        Self {
            limits,
            remaining_run_work: None,
            nodes: FxHashMap::default(),
            stats: InternerStats::default(),
        }
    }
    pub(crate) fn with_run_budget(mut self) -> Self {
        self.remaining_run_work = Some(self.limits.max_work);
        self
    }
    /// Merge phase rotations in an explicit projective circuit scope.
    /// The scalar add/sub APIs retain literal coefficients and periods.
    #[inline]
    pub(crate) fn merge_phases(&mut self, a: &Angle, sign: i8, b: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a, b])?;
        if let (Some(a), Some(b)) = (
            a.quarter_turns(RotationPeriod::TwoPi),
            b.quarter_turns(RotationPeriod::TwoPi),
        ) {
            let out =
                Angle::phase_quarters((sign as i16 * a as i16 + b as i16).rem_euclid(8) as u8);
            budget.value(&out)?;
            return Ok(out);
        }
        let a = if sign < 0 {
            self.scale_inner(a, &Rational::integer(-1), &mut budget)?
        } else {
            a.clone()
        };
        let value = self.add_inner(&a, b, &mut budget)?;
        self.rotation_view(&value, RotationPeriod::TwoPi)
    }
    pub fn limits(&self) -> Limits {
        self.limits
    }
    pub fn stats(&self) -> InternerStats {
        self.stats
    }
    /// Remove expired weak entries without affecting any live angle.
    pub fn prune(&mut self) -> usize {
        let mut removed = 0;
        self.nodes.retain(|_, bucket| {
            let before = bucket.len();
            bucket.retain(|node| node.strong_count() != 0);
            removed += before - bucket.len();
            !bucket.is_empty()
        });
        self.stats.retained_entries -= removed;
        self.stats.pruned = self.stats.pruned.saturating_add(removed);
        removed
    }
    /// Drop the lookup table. Existing values remain valid and independent.
    pub fn clear_interner(&mut self) {
        self.nodes.clear();
        self.stats.retained_entries = 0;
    }

    #[inline]
    pub fn rational(&mut self, n: i64, d: i64) -> Result<Angle, Error> {
        let mut budget = Budget::new(self.limits);
        let r = Rational::new(n, d)?;
        self.number(r, false, &mut budget)
    }
    #[inline]
    pub fn pi_fraction(&mut self, n: i64, d: i64) -> Result<Angle, Error> {
        let mut budget = Budget::new(self.limits);
        let r = Rational::new(n, d)?;
        self.number(r, true, &mut budget)
    }
    pub fn big_rational(&mut self, n: BigInt, d: BigInt) -> Result<Angle, Error> {
        let mut budget = Budget::new(self.limits);
        let r = Rational::from_big(n, d, &mut budget)?;
        self.number(r, false, &mut budget)
    }
    pub fn big_pi_fraction(&mut self, n: BigInt, d: BigInt) -> Result<Angle, Error> {
        let mut budget = Budget::new(self.limits);
        let r = Rational::from_big(n, d, &mut budget)?;
        self.number(r, true, &mut budget)
    }
    #[inline]
    pub fn from_f64(&mut self, value: f64) -> Result<Angle, Error> {
        let value = Angle::from_f64(value)?;
        Budget::new(self.limits).value(&value)?;
        Ok(value)
    }
    pub fn decimal(&mut self, source: &str) -> Result<Angle, Error> {
        self.decimal_inner(source, &mut Budget::new(self.limits))
    }
    #[inline]
    pub fn add(&mut self, a: &Angle, b: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a, b])?;
        self.add_inner(a, b, &mut budget)
    }
    pub fn sub(&mut self, a: &Angle, b: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a, b])?;
        let b = self.scale_inner(b, &Rational::integer(-1), &mut budget)?;
        self.add_inner(a, &b, &mut budget)
    }
    #[inline]
    pub fn neg(&mut self, a: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a])?;
        self.scale_inner(a, &Rational::integer(-1), &mut budget)
    }
    #[inline]
    pub fn scale(&mut self, a: &Angle, coefficient: &Rational) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a])?;
        budget.bits(coefficient.bits())?;
        self.scale_inner(a, coefficient, &mut budget)
    }
    #[inline]
    pub fn half(&mut self, a: &Angle) -> Result<Angle, Error> {
        self.scale(a, &Rational::canonical(1, 2))
    }
    pub fn mul(&mut self, a: &Angle, b: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a, b])?;
        self.mul_inner(a, b, &mut budget)
    }
    pub fn div(&mut self, a: &Angle, b: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a, b])?;
        self.div_inner(a, b, &mut budget)
    }
    pub fn pow(&mut self, a: &Angle, exponent: i64) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a])?;
        self.pow_inner(a, exponent, &mut budget)
    }
    pub fn power(&mut self, a: &Angle, b: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a, b])?;
        self.power_inner(a, b, &mut budget)
    }
    pub fn function(&mut self, fun: Function, a: &Angle) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a])?;
        self.function_inner(fun, a, &mut budget)
    }
    pub fn sin(&mut self, a: &Angle) -> Result<Angle, Error> {
        self.function(Function::Sin, a)
    }
    pub fn cos(&mut self, a: &Angle) -> Result<Angle, Error> {
        self.function(Function::Cos, a)
    }
    pub fn tan(&mut self, a: &Angle) -> Result<Angle, Error> {
        self.function(Function::Tan, a)
    }
    pub fn exp(&mut self, a: &Angle) -> Result<Angle, Error> {
        self.function(Function::Exp, a)
    }
    pub fn ln(&mut self, a: &Angle) -> Result<Angle, Error> {
        self.function(Function::Ln, a)
    }
    pub fn sqrt(&mut self, a: &Angle) -> Result<Angle, Error> {
        self.function(Function::Sqrt, a)
    }
    pub fn prove_equal(&mut self, a: &Angle, b: &Angle) -> Result<Proof, Error> {
        let mut budget = self.inputs(&[a, b])?;
        if a == b {
            return Ok(Proof::Proven);
        }
        let neg = self.scale_inner(b, &Rational::integer(-1), &mut budget)?;
        Ok(self.add_inner(a, &neg, &mut budget)?.prove_zero())
    }
    /// Explicit reduction of a pure pi coefficient in a chosen gate view.
    /// General expressions are left intact; their periodicity is not guessed.
    pub fn rotation_view(&mut self, a: &Angle, period: RotationPeriod) -> Result<Angle, Error> {
        let mut budget = self.inputs(&[a])?;
        let Some(r) = a.pi_coefficient() else {
            return Ok(a.clone());
        };
        let r = if period == RotationPeriod::FourPi {
            let half = r.div(&Rational::integer(2), &mut budget)?;
            half.reduced_pi(&mut budget)?
                .mul(&Rational::integer(2), &mut budget)?
        } else {
            r.reduced_pi(&mut budget)?
        };
        self.number(r, true, &mut budget)
    }
    pub(super) fn charge_run(&mut self, values: &[&Angle]) -> Result<(), Error> {
        if let Some(remaining) = &mut self.remaining_run_work {
            let mut cost = 0usize;
            for value in values {
                let m = value.metadata();
                if m.nodes > 1 || m.bits > 128 {
                    cost = cost
                        .saturating_add(m.nodes.saturating_mul(m.bits.div_ceil(64).max(1).pow(2)));
                }
            }
            *remaining = remaining
                .checked_sub(cost)
                .ok_or(Error::Limit(Resource::Work))?;
        }
        Ok(())
    }
    #[inline]
    fn inputs(&mut self, values: &[&Angle]) -> Result<Budget, Error> {
        self.charge_run(values)?;
        let mut budget = Budget::new(self.limits);
        for value in values {
            budget.value(value)?;
        }
        Ok(budget)
    }
    pub(super) fn number(
        &mut self,
        r: Rational,
        pi: bool,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        budget.charge(1)?;
        budget.bits(r.bits())?;
        if budget.limits.max_nodes == 0 {
            return Err(Error::Limit(Resource::Nodes));
        }
        if budget.limits.max_depth == 0 {
            return Err(Error::Limit(Resource::Depth));
        }
        if r.is_zero() {
            return Ok(Angle::ZERO);
        }
        if let Some((n, d)) = r.as_small()
            && let (Ok(n), Ok(d)) = (i32::try_from(n), u32::try_from(d))
        {
            return Ok(Angle(if pi {
                Value::Pi { n, d }
            } else {
                Value::Rational { n, d }
            }));
        }
        if !pi && let Some(value) = dyadic::from_rational(&r) {
            return Ok(value);
        }
        let sign = Sign::from_integer(r.signum());
        self.intern(
            if pi {
                NodeKind::Pi(r)
            } else {
                NodeKind::Rational(r)
            },
            sign,
            budget,
        )
    }
    fn intern(&mut self, kind: NodeKind, sign: Sign, budget: &mut Budget) -> Result<Angle, Error> {
        let mut meta = Metadata {
            nodes: 1,
            depth: 1,
            bits: 1,
            sign,
        };
        let mut child = |a: &Angle| -> Result<(), Error> {
            let m = a.metadata();
            meta.nodes = meta
                .nodes
                .checked_add(m.nodes)
                .ok_or(Error::Limit(Resource::Nodes))?;
            meta.depth = meta.depth.max(
                m.depth
                    .checked_add(1)
                    .ok_or(Error::Limit(Resource::Depth))?,
            );
            meta.bits = meta.bits.max(m.bits);
            Ok(())
        };
        match &kind {
            NodeKind::Rational(r) | NodeKind::Pi(r) => meta.bits = r.bits(),
            NodeKind::Sum(terms) => {
                if terms.len() > budget.limits.max_terms {
                    return Err(Error::Limit(Resource::Terms));
                }
                for term in terms {
                    child(&term.atom)?;
                }
                for term in terms {
                    meta.bits = meta.bits.max(term.coefficient.bits());
                }
            }
            NodeKind::Product {
                coefficient,
                factors,
            } => {
                if factors.len() > budget.limits.max_terms {
                    return Err(Error::Limit(Resource::Terms));
                }
                for a in factors {
                    child(a)?;
                }
                meta.bits = meta.bits.max(coefficient.bits());
            }
            NodeKind::Quotient(a, b) => {
                child(a)?;
                child(b)?;
            }
            NodeKind::Power(a, e) => {
                child(a)?;
                meta.bits = meta
                    .bits
                    .max((64 - e.unsigned_abs().leading_zeros()) as usize);
            }
            NodeKind::Function(_, a) => child(a)?,
        }
        if meta.depth > budget.limits.max_depth {
            return Err(Error::Limit(Resource::Depth));
        }
        if meta.nodes > budget.limits.max_nodes {
            return Err(Error::Limit(Resource::Nodes));
        }
        budget.bits(meta.bits)?;
        budget.charge(meta.nodes)?;
        let hash = node_hash(&kind);
        if let Some(value) = self.lookup(&kind, hash, budget)? {
            return Ok(value);
        }
        if self.stats.retained_entries >= self.limits.max_nodes {
            budget.charge(self.stats.retained_entries)?;
            self.prune();
            if self.stats.retained_entries >= self.limits.max_nodes {
                return Err(Error::Limit(Resource::Nodes));
            }
        }
        let node = Arc::new(Node { kind, hash, meta });
        self.nodes
            .entry(hash)
            .or_default()
            .push(Arc::downgrade(&node));
        self.stats.misses = self.stats.misses.saturating_add(1);
        self.stats.retained_entries += 1;
        Ok(Angle(Value::Shared(node)))
    }

    fn lookup(
        &mut self,
        kind: &NodeKind,
        hash: u64,
        budget: &mut Budget,
    ) -> Result<Option<Angle>, Error> {
        if let Some(bucket) = self.nodes.get_mut(&hash) {
            let mut index = 0;
            while index < bucket.len() {
                budget.charge(1)?;
                if let Some(node) = bucket[index].upgrade() {
                    budget.charge(node.meta.nodes)?;
                    if &node.kind == kind {
                        self.stats.hits = self.stats.hits.saturating_add(1);
                        return Ok(Some(Angle(Value::Shared(node))));
                    }
                    index += 1;
                } else {
                    bucket.swap_remove(index);
                    self.stats.retained_entries -= 1;
                    self.stats.pruned = self.stats.pruned.saturating_add(1);
                }
            }
        }
        Ok(None)
    }

    #[inline]
    pub(super) fn add_inner(
        &mut self,
        a: &Angle,
        b: &Angle,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        budget.charge(1)?;
        if a.prove_zero() == Proof::Proven {
            return Ok(b.clone());
        }
        if b.prove_zero() == Proof::Proven {
            return Ok(a.clone());
        }
        // Equal or power-of-two denominators cover quarter turns and fused rotations.
        // Reduce using machine integers before creating a Rational or a node.
        let inline = match (&a.0, &b.0) {
            (Value::Pi { n: an, d: ad }, Value::Pi { n: bn, d: bd })
                if ad == bd || (ad.is_power_of_two() && bd.is_power_of_two()) =>
            {
                Some((*an, *bn, *ad, *bd, true))
            }
            (Value::Rational { n: an, d: ad }, Value::Rational { n: bn, d: bd })
                if ad == bd || (ad.is_power_of_two() && bd.is_power_of_two()) =>
            {
                Some((*an, *bn, *ad, *bd, false))
            }
            _ => None,
        };
        if let Some((an, bn, ad, bd, pi)) = inline {
            use num_integer::Integer;
            let d = ad.max(bd);
            let n = if ad == bd {
                an as i64 + bn as i64
            } else {
                ((an as i64) << (d.trailing_zeros() - ad.trailing_zeros()))
                    + ((bn as i64) << (d.trailing_zeros() - bd.trailing_zeros()))
            };
            let gcd = if d.is_power_of_two() {
                1_u64 << n.unsigned_abs().trailing_zeros().min(d.trailing_zeros())
            } else {
                n.unsigned_abs().gcd(&(d as u64))
            };
            if let Ok(n) = i32::try_from(n / gcd as i64) {
                let d = d / gcd as u32;
                let result = if n == 0 {
                    Angle::ZERO
                } else {
                    Angle(if pi {
                        Value::Pi { n, d }
                    } else {
                        Value::Rational { n, d }
                    })
                };
                budget.charge(1)?;
                budget.value(&result)?;
                return Ok(result);
            }
        }
        if let (Value::Dyadic(a), Value::Dyadic(b)) = (&a.0, &b.0)
            && let Some(result) = dyadic::add(*a, *b)
        {
            budget.charge(1)?;
            budget.value(&result)?;
            return Ok(result);
        }
        // Recognize a single negative term directly, avoiding temporary vectors.
        for (term, atom) in [(a, b), (b, a)] {
            if let Some(NodeKind::Sum(terms)) = term.kind()
                && terms.len() == 1
                && terms[0].coefficient.as_small() == Some((-1, 1))
            {
                budget.charge(atom.expression_nodes())?;
                if terms[0].atom == *atom {
                    return Ok(Angle::ZERO);
                }
            }
        }
        if let (Some(a), Some(b)) = (a.as_rational(), b.as_rational()) {
            return self.number(a.add(&b, budget)?, false, budget);
        }
        if let (Some(a), Some(b)) = (a.pi_coefficient(), b.pi_coefficient()) {
            return self.number(a.add(&b, budget)?, true, budget);
        }
        let mut terms = Vec::new();
        self.terms(a, &Rational::ONE, &mut terms, budget)?;
        self.terms(b, &Rational::ONE, &mut terms, budget)?;
        self.sum(terms, budget)
    }
    pub(super) fn terms(
        &mut self,
        a: &Angle,
        scale: &Rational,
        out: &mut Vec<Term>,
        budget: &mut Budget,
    ) -> Result<(), Error> {
        budget.charge(a.expression_nodes())?;
        if let Some(r) = a.as_rational() {
            out.push(Term {
                coefficient: r.mul(scale, budget)?,
                atom: Angle::ONE,
            });
        } else if let Some(r) = a.pi_coefficient() {
            out.push(Term {
                coefficient: r.mul(scale, budget)?,
                atom: Angle::PI,
            });
        } else if let Some(NodeKind::Sum(terms)) = a.kind() {
            if out.len().saturating_add(terms.len()) > budget.limits.max_terms {
                return Err(Error::Limit(Resource::Terms));
            }
            for term in terms {
                out.push(Term {
                    coefficient: term.coefficient.mul(scale, budget)?,
                    atom: term.atom.clone(),
                });
            }
        } else if let Some(NodeKind::Product {
            coefficient,
            factors,
        }) = a.kind()
        {
            let atom = if coefficient.is_one() {
                a.clone()
            } else {
                self.intern(
                    NodeKind::Product {
                        coefficient: Rational::ONE,
                        factors: factors.clone(),
                    },
                    product_sign(&Rational::ONE, factors),
                    budget,
                )?
            };
            out.push(Term {
                coefficient: coefficient.mul(scale, budget)?,
                atom,
            });
        } else {
            out.push(Term {
                coefficient: scale.clone(),
                atom: a.clone(),
            });
        }
        if out.len() > budget.limits.max_terms {
            return Err(Error::Limit(Resource::Terms));
        }
        Ok(())
    }
    pub(super) fn sum(
        &mut self,
        mut terms: Vec<Term>,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        self.sort_budget(
            terms.len(),
            terms.iter().map(|t| t.atom.expression_nodes()).sum(),
            budget,
        )?;
        terms.sort_unstable_by(|a, b| a.atom.cmp(&b.atom));
        let mut merged: Vec<Term> = Vec::with_capacity(terms.len());
        for term in terms {
            if let Some(last) = merged.last_mut().filter(|last| last.atom == term.atom) {
                last.coefficient = last.coefficient.add(&term.coefficient, budget)?;
            } else {
                merged.push(term);
            }
        }
        merged.retain(|term| !term.coefficient.is_zero());
        if merged.is_empty() {
            return Ok(Angle::ZERO);
        }
        // c*sin(x)^2 + c*cos(x)^2 = c, wherever the pair occurs in the sum.
        let squares: Vec<usize> = (0..merged.len())
            .filter(|&i| trig_square(&merged[i].atom).is_some())
            .collect();
        if squares.len() >= 2 {
            budget.charge(squares.len() * squares.len())?;
            let pair = squares.iter().enumerate().find_map(|(k, &i)| {
                let (a, x) = trig_square(&merged[i].atom)?;
                squares[k + 1..]
                    .iter()
                    .copied()
                    .find(|&j| {
                        merged[i].coefficient == merged[j].coefficient
                            && trig_square(&merged[j].atom).is_some_and(|(b, y)| a != b && x == y)
                    })
                    .map(|j| (i, j))
            });
            if let Some((i, j)) = pair {
                // i < j, so removing j first leaves i in place.
                merged.swap_remove(j);
                let Term { coefficient, .. } = merged.swap_remove(i);
                merged.push(Term {
                    coefficient,
                    atom: Angle::ONE,
                });
                return self.sum(merged, budget);
            }
        }
        if merged.len() == 1 {
            let term = merged.pop().unwrap();
            if term.atom == Angle::ONE {
                return self.number(term.coefficient, false, budget);
            }
            if term.atom == Angle::PI {
                return self.number(term.coefficient, true, budget);
            }
            if term.coefficient.is_one() {
                return Ok(term.atom);
            }
            if matches!(term.atom.kind(), Some(NodeKind::Product { .. })) {
                return self.scale_inner(&term.atom, &term.coefficient, budget);
            }
            merged.push(term);
        }
        let sign = sum_sign(&merged);
        self.intern(NodeKind::Sum(merged.into_boxed_slice()), sign, budget)
    }
    fn sort_budget(&self, length: usize, nodes: usize, budget: &mut Budget) -> Result<(), Error> {
        let log = (usize::BITS - length.max(1).leading_zeros()) as usize;
        budget.charge(nodes.checked_mul(log).ok_or(Error::Limit(Resource::Work))?)
    }
    pub(super) fn scale_inner(
        &mut self,
        a: &Angle,
        r: &Rational,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        budget.charge(1)?;
        if r.is_zero() || a.prove_zero() == Proof::Proven {
            return Ok(Angle::ZERO);
        }
        if r.is_one() {
            return Ok(a.clone());
        }
        if let Value::Dyadic(bits) = &a.0
            && let Some((n, d)) = r.as_small()
            && n.unsigned_abs().is_power_of_two()
            && d.is_power_of_two()
            && let Some(result) = dyadic::scale(
                *bits,
                n.unsigned_abs().trailing_zeros() as i32 - d.trailing_zeros() as i32,
                n < 0,
            )
        {
            budget.charge(1)?;
            budget.value(&result)?;
            return Ok(result);
        }
        if let Some(v) = a.as_rational() {
            return self.number(v.mul(r, budget)?, false, budget);
        }
        if let Some(v) = a.pi_coefficient() {
            return self.number(v.mul(r, budget)?, true, budget);
        }
        if let Some(NodeKind::Product {
            coefficient,
            factors,
        }) = a.kind()
        {
            let coefficient = coefficient.mul(r, budget)?;
            let sign = product_sign(&coefficient, factors);
            return self.intern(
                NodeKind::Product {
                    coefficient,
                    factors: factors.clone(),
                },
                sign,
                budget,
            );
        }
        let mut terms = Vec::new();
        self.terms(a, r, &mut terms, budget)?;
        self.sum(terms, budget)
    }
    pub(super) fn mul_inner(
        &mut self,
        a: &Angle,
        b: &Angle,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        if let Some(r) = a.as_rational() {
            return self.scale_inner(b, &r, budget);
        }
        if let Some(r) = b.as_rational() {
            return self.scale_inner(a, &r, budget);
        }
        let mut coefficient = Rational::ONE;
        let mut factors = Vec::new();
        for value in [a, b] {
            self.push_factor(value, &mut coefficient, &mut factors, budget)?;
        }
        // Simplifying a merged power (sqrt(x)^2 = x) can expose a coefficient
        // or a base shared with another factor, so regroup until stable. Each
        // regrouping removes a square root, which bounds the iterations.
        let mut combined = loop {
            if factors.len() > budget.limits.max_terms {
                return Err(Error::Limit(Resource::Terms));
            }
            self.sort_budget(
                factors.len(),
                factors.iter().map(Angle::expression_nodes).sum(),
                budget,
            )?;
            factors.sort_unstable_by(|a, b| power_factor(a).0.cmp(power_factor(b).0));
            let mut combined = Vec::with_capacity(factors.len());
            let mut regroup = false;
            let mut index = 0;
            while index < factors.len() {
                let (base, exponent) = power_factor(&factors[index]);
                let mut exponent = exponent as i128;
                let mut end = index + 1;
                while end < factors.len() && power_factor(&factors[end]).0 == base {
                    exponent += power_factor(&factors[end]).1 as i128;
                    end += 1;
                }
                // Keep valid factors when their combined exponent exceeds i64;
                // normalization must not reject a representable exact expression.
                if end - index == 1 {
                    combined.push(factors[index].clone());
                } else if let Ok(exponent) = i64::try_from(exponent) {
                    let factor = self.pow_inner(base, exponent, budget)?;
                    if factor.as_rational().is_none() && power_factor(&factor).0 == base {
                        combined.push(factor);
                    } else {
                        regroup |= factor.as_rational().is_none();
                        self.push_factor(&factor, &mut coefficient, &mut combined, budget)?;
                    }
                } else {
                    combined.extend(factors[index..end].iter().cloned());
                }
                index = end;
            }
            if !regroup {
                break combined;
            }
            factors = combined;
        };
        if combined.is_empty() {
            return self.number(coefficient, false, budget);
        }
        if combined.len() == 1 {
            return self.scale_inner(&combined[0], &coefficient, budget);
        }
        self.sort_budget(
            combined.len(),
            combined.iter().map(Angle::expression_nodes).sum(),
            budget,
        )?;
        combined.sort_unstable();
        let sign = product_sign(&coefficient, &combined);
        self.intern(
            NodeKind::Product {
                coefficient,
                factors: combined.into_boxed_slice(),
            },
            sign,
            budget,
        )
    }
    /// Add a multiplicand to a product in normal form: rational and pi
    /// coefficients, nested products, single scaled terms, and the content of
    /// a sum are all moved into the product's coefficient.
    fn push_factor(
        &mut self,
        value: &Angle,
        coefficient: &mut Rational,
        factors: &mut Vec<Angle>,
        budget: &mut Budget,
    ) -> Result<(), Error> {
        budget.charge(1)?;
        if let Some(r) = value.as_rational() {
            *coefficient = coefficient.mul(&r, budget)?;
        } else if let Some(r) = value.pi_coefficient() {
            *coefficient = coefficient.mul(&r, budget)?;
            factors.push(Angle::PI);
        } else if let Some((content, primitive)) = self.primitive_sum(value, budget)? {
            *coefficient = coefficient.mul(&content, budget)?;
            factors.push(primitive);
        } else {
            match value.kind() {
                Some(NodeKind::Product {
                    coefficient: c,
                    factors: f,
                }) => {
                    *coefficient = coefficient.mul(c, budget)?;
                    factors.extend(f.iter().cloned());
                }
                Some(NodeKind::Sum(terms)) if terms.len() == 1 => {
                    *coefficient = coefficient.mul(&terms[0].coefficient, budget)?;
                    self.push_factor(&terms[0].atom, coefficient, factors, budget)?;
                }
                _ => factors.push(value.clone()),
            }
        }
        Ok(())
    }
    /// Split a sum of two or more terms into its content and a primitive sum
    /// whose leading term has coefficient one, unless it already is one.
    /// Sums stay distributed on their own (`2*(a+b)` is `2a+2b`), but as a
    /// factor or a power's base a sum is primitive, so that `2*(a+b)*c` and
    /// `(a+b)*(2*c)` share one representation.
    fn primitive_sum(
        &mut self,
        a: &Angle,
        budget: &mut Budget,
    ) -> Result<Option<(Rational, Angle)>, Error> {
        let Some(NodeKind::Sum(terms)) = a.kind() else {
            return Ok(None);
        };
        if terms.len() < 2 || terms[0].coefficient.is_one() {
            return Ok(None);
        }
        let content = terms[0].coefficient.clone();
        let inverse = Rational::ONE.div(&content, budget)?;
        let primitive = self.scale_inner(a, &inverse, budget)?;
        Ok(Some((content, primitive)))
    }
    pub(super) fn div_inner(
        &mut self,
        a: &Angle,
        b: &Angle,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        match b.prove_zero() {
            Proof::Proven => return Err(Error::DivisionByZero),
            Proof::Unknown => return Err(Error::DomainUnproved),
            _ => {}
        }
        if a == b {
            return Ok(Angle::ONE);
        }
        if let Some(r) = b.as_rational() {
            let r = Rational::ONE.div(&r, budget)?;
            return self.scale_inner(a, &r, budget);
        }
        if a.prove_zero() == Proof::Proven {
            return Ok(Angle::ZERO);
        }
        if let (Some(a), Some(b)) = (a.pi_coefficient(), b.pi_coefficient()) {
            return self.number(a.div(&b, budget)?, false, budget);
        }
        let sign = sign_product(a.metadata().sign, b.metadata().sign);
        self.intern(NodeKind::Quotient(a.clone(), b.clone()), sign, budget)
    }
    pub(super) fn pow_inner(
        &mut self,
        a: &Angle,
        e: i64,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        if e <= 0 {
            match a.prove_zero() {
                Proof::Proven => {
                    return Err(if e == 0 {
                        Error::ZeroPowerZero
                    } else {
                        Error::DivisionByZero
                    });
                }
                Proof::Unknown => return Err(Error::DomainUnproved),
                _ => {}
            }
        }
        if e == 0 {
            return Ok(Angle::ONE);
        }
        if e == 1 {
            return Ok(a.clone());
        }
        if let Some(r) = a.as_rational() {
            return self.number(r.pow(e, budget)?, false, budget);
        }
        if let Some(NodeKind::Function(Function::Sqrt, arg)) = a.kind()
            && e % 2 == 0
        {
            return self.pow_inner(arg, e / 2, budget);
        }
        // Distribute over rational coefficients: (c*x)^e = c^e * x^e.
        if let Some(r) = a.pi_coefficient()
            && !r.is_one()
        {
            let coefficient = r.pow(e, budget)?;
            let power = self.pow_inner(&Angle::PI, e, budget)?;
            return self.scale_inner(&power, &coefficient, budget);
        }
        if let Some(NodeKind::Sum(terms)) = a.kind()
            && terms.len() == 1
        {
            let coefficient = terms[0].coefficient.pow(e, budget)?;
            let power = self.pow_inner(&terms[0].atom, e, budget)?;
            return self.scale_inner(&power, &coefficient, budget);
        }
        if let Some((content, primitive)) = self.primitive_sum(a, budget)? {
            let coefficient = content.pow(e, budget)?;
            let power = self.pow_inner(&primitive, e, budget)?;
            return self.scale_inner(&power, &coefficient, budget);
        }
        // A proven nonzero product has proven nonzero factors, so negative
        // exponents distribute as well.
        if let Some(NodeKind::Product {
            coefficient,
            factors,
        }) = a.kind()
        {
            let coefficient = coefficient.pow(e, budget)?;
            let mut power = Angle::ONE;
            for factor in factors {
                let factor = self.pow_inner(factor, e, budget)?;
                power = self.mul_inner(&power, &factor, budget)?;
            }
            return self.scale_inner(&power, &coefficient, budget);
        }
        if let Some(NodeKind::Power(base, old)) = a.kind()
            && let Some(e) = old.checked_mul(e)
        {
            return self.pow_inner(base, e, budget);
        }
        budget.charge(e.unsigned_abs().ilog2() as usize + 1)?;
        let sign = if e & 1 != 0 {
            a.metadata().sign
        } else {
            match a.prove_zero() {
                Proof::Disproven => Sign::Positive,
                Proof::Proven => Sign::Zero,
                Proof::Unknown => Sign::NonNegative,
            }
        };
        self.intern(NodeKind::Power(a.clone(), e), sign, budget)
    }

    /// Real-valued `a^b`. Integer exponents stay exact powers; other exponents
    /// use the principal value `exp(b·ln a)`, real only for `a > 0` (and
    /// `0^b = 0` for proven-positive `b`).
    pub(super) fn power_inner(
        &mut self,
        a: &Angle,
        b: &Angle,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        let exponent = b.as_rational();
        if let Some(r) = exponent.as_ref().filter(|r| r.is_integer()) {
            let (e, _) = r.as_small().ok_or(Error::Limit(Resource::Work))?;
            return self.pow_inner(a, e, budget);
        }
        match a.prove_positive() {
            Proof::Proven => {}
            Proof::Unknown => return Err(Error::DomainUnproved),
            Proof::Disproven => {
                return match (a.prove_zero(), b.prove_positive()) {
                    (Proof::Proven, Proof::Proven) => Ok(Angle::ZERO),
                    (Proof::Proven, Proof::Unknown) => Err(Error::DomainUnproved),
                    _ => Err(Error::NonRealPower),
                };
            }
        }
        if exponent.is_some_and(|r| r == Rational::canonical(1, 2)) {
            return self.function_inner(Function::Sqrt, a, budget);
        }
        let ln = self.function_inner(Function::Ln, a, budget)?;
        let product = self.mul_inner(b, &ln, budget)?;
        self.function_inner(Function::Exp, &product, budget)
    }

    pub(super) fn function_inner(
        &mut self,
        fun: Function,
        a: &Angle,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        budget.charge(1)?;
        let key = NodeKind::Function(fun, a.clone());
        if let Some(value) = self.lookup(&key, node_hash(&key), budget)? {
            return Ok(value);
        }
        if fun == Function::Ln {
            match a.prove_positive() {
                Proof::Disproven => return Err(Error::InvalidDomain(fun)),
                Proof::Unknown => return Err(Error::DomainUnproved),
                _ => {}
            }
            if a == &Angle::ONE {
                return Ok(Angle::ZERO);
            }
            if let Some(NodeKind::Function(Function::Exp, arg)) = a.kind() {
                return Ok(arg.clone());
            }
        }
        if fun == Function::Sqrt {
            match a.prove_nonnegative() {
                Proof::Disproven => return Err(Error::InvalidDomain(fun)),
                Proof::Unknown => return Err(Error::DomainUnproved),
                _ => {}
            }
            if let Some(r) = a.as_rational()
                && let Some(root) = r.perfect_sqrt(budget)?
            {
                return self.number(root, false, budget);
            }
            if let Some(NodeKind::Power(base, 2)) = a.kind() {
                match base.metadata().sign {
                    Sign::Positive | Sign::NonNegative | Sign::Zero => return Ok(base.clone()),
                    Sign::Negative | Sign::NonPositive => {
                        return self.scale_inner(base, &Rational::integer(-1), budget);
                    }
                    _ => {}
                }
            }
        }
        if fun == Function::Exp {
            if a.prove_zero() == Proof::Proven {
                return Ok(Angle::ONE);
            }
            if let Some(NodeKind::Function(Function::Ln, arg)) = a.kind() {
                return Ok(arg.clone());
            }
        }
        if matches!(fun, Function::Sin | Function::Cos | Function::Tan) {
            if a.prove_zero() == Proof::Proven {
                return Ok(if fun == Function::Cos {
                    Angle::ONE
                } else {
                    Angle::ZERO
                });
            }
            if let Some(k) = a.quarter_turns(RotationPeriod::TwoPi) {
                let sign = match fun {
                    Function::Sin => match k {
                        0 | 4 => return Ok(Angle::ZERO),
                        2 => return Ok(Angle::ONE),
                        6 => return self.number(Rational::integer(-1), false, budget),
                        1 | 3 => 1,
                        _ => -1,
                    },
                    Function::Cos => match k {
                        2 | 6 => return Ok(Angle::ZERO),
                        0 => return Ok(Angle::ONE),
                        4 => return self.number(Rational::integer(-1), false, budget),
                        1 | 7 => 1,
                        _ => -1,
                    },
                    Function::Tan => {
                        return match k {
                            2 | 6 => Err(Error::InvalidDomain(fun)),
                            0 | 4 => Ok(Angle::ZERO),
                            1 | 5 => Ok(Angle::ONE),
                            _ => self.number(Rational::integer(-1), false, budget),
                        };
                    }
                    _ => unreachable!(),
                };
                let two = self.number(Rational::integer(2), false, budget)?;
                let root = self.function_inner(Function::Sqrt, &two, budget)?;
                return self.scale_inner(&root, &Rational::new(sign, 2)?, budget);
            }
            if let Some(q) = a.pi_coefficient() {
                let reduced = if fun == Function::Tan {
                    q.mul(&Rational::integer(2), budget)?
                        .reduced_pi(budget)?
                        .div(&Rational::integer(2), budget)?
                } else {
                    q.reduced_pi(budget)?
                };
                if q != reduced {
                    let a = self.number(reduced, true, budget)?;
                    return self.function_inner(fun, &a, budget);
                }
            }
            if let Some(positive) = self.strip_negative(a, budget)? {
                let value = self.function_inner(fun, &positive, budget)?;
                return if fun == Function::Cos {
                    Ok(value)
                } else {
                    self.scale_inner(&value, &Rational::integer(-1), budget)
                };
            }
            if fun == Function::Tan {
                let cosine = self.function_inner(Function::Cos, a, budget)?;
                match cosine.prove_zero() {
                    Proof::Proven => return Err(Error::InvalidDomain(fun)),
                    Proof::Unknown => return Err(Error::DomainUnproved),
                    _ => {}
                }
            }
        }
        let sign = match fun {
            Function::Exp => Sign::Positive,
            Function::Sqrt => {
                if a.prove_zero() == Proof::Disproven {
                    Sign::Positive
                } else {
                    Sign::NonNegative
                }
            }
            Function::Ln => match a.as_rational() {
                Some(r) => Sign::from_integer(match r.cmp(&Rational::ONE) {
                    Ordering::Less => -1,
                    Ordering::Equal => 0,
                    Ordering::Greater => 1,
                }),
                _ => Sign::Unknown,
            },
            Function::Sin | Function::Cos => self.trig_sign(fun, a, budget)?,
            Function::Tan => sign_product(
                self.trig_sign(Function::Sin, a, budget)?,
                self.trig_sign(Function::Cos, a, budget)?,
            ),
        };
        self.intern(NodeKind::Function(fun, a.clone()), sign, budget)
    }
    fn strip_negative(&mut self, a: &Angle, budget: &mut Budget) -> Result<Option<Angle>, Error> {
        let negative = if let Some(r) = a.as_rational() {
            r.signum() < 0
        } else if let Some(r) = a.pi_coefficient() {
            r.signum() < 0
        } else {
            match a.kind() {
                Some(NodeKind::Sum(terms)) => terms.iter().all(|t| t.coefficient.signum() < 0),
                Some(NodeKind::Product { coefficient, .. }) => coefficient.signum() < 0,
                _ => false,
            }
        };
        if negative {
            Ok(Some(self.scale_inner(a, &Rational::integer(-1), budget)?))
        } else {
            Ok(None)
        }
    }
    fn trig_sign(&self, fun: Function, a: &Angle, budget: &mut Budget) -> Result<Sign, Error> {
        if let Some(q) = a.pi_coefficient() {
            let q = q.reduced_pi(budget)?;
            if fun == Function::Sin {
                if q.is_zero() || q == Rational::integer(-1) {
                    return Ok(Sign::Zero);
                }
                return Ok(Sign::from_integer(q.signum()));
            }
            let q = if q.signum() < 0 { q.neg(budget)? } else { q };
            return Ok(match q.cmp(&Rational::new(1, 2)?) {
                Ordering::Less => Sign::Positive,
                Ordering::Equal => Sign::Zero,
                Ordering::Greater => Sign::Negative,
            });
        }
        if let Some(r) = a.as_rational() {
            // Proven analytic facts: sin(x)>0 for 0<x<3<pi;
            // cos(x)>=1-x*x/2>=1/2 for |x|<=1.
            if fun == Function::Sin && r > Rational::integer(-3) && r < Rational::integer(3) {
                return Ok(Sign::from_integer(r.signum()));
            }
            if fun == Function::Cos && r >= Rational::integer(-1) && r <= Rational::ONE {
                return Ok(Sign::Positive);
            }
        }
        Ok(Sign::Unknown)
    }
    pub(super) fn decimal_inner(
        &mut self,
        source: &str,
        budget: &mut Budget,
    ) -> Result<Angle, Error> {
        if source.len() > budget.limits.max_source_bytes {
            return Err(Error::Limit(Resource::Source));
        }
        budget.charge(source.len())?;
        let invalid = || Error::Parse {
            offset: 0,
            message: "invalid decimal literal",
        };
        let (mantissa, exponent) = source.split_once(['e', 'E']).unwrap_or((source, "0"));
        let exponent: i64 = exponent.parse().map_err(|_| invalid())?;
        let (negative, mantissa) = if let Some(s) = mantissa.strip_prefix('-') {
            (true, s)
        } else {
            (false, mantissa.strip_prefix('+').unwrap_or(mantissa))
        };
        let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        if whole.is_empty() && fraction.is_empty() {
            return Err(invalid());
        }
        if !whole
            .bytes()
            .chain(fraction.bytes())
            .all(|b| b.is_ascii_digit())
        {
            return Err(invalid());
        }
        let significant_digits = whole
            .bytes()
            .chain(fraction.bytes())
            .skip_while(|b| *b == b'0')
            .count();
        if significant_digits == 0 {
            return self.number(Rational::ZERO, false, budget);
        }
        let scale = (fraction.len() as i64)
            .checked_sub(exponent)
            .ok_or(Error::Limit(Resource::Bits))?;
        let scale_digits =
            usize::try_from(scale.unsigned_abs()).map_err(|_| Error::Limit(Resource::Bits))?;
        // The value is a D-digit integer times 10^-scale, so its magnitude
        // lies in [10^(D-1-scale), 10^(D-scale)). A reduced rational needs at
        // least log2 of the larger of |value| and 1/|value| bits; reject only
        // when that lower bound already exceeds the limit (3.321 < log2 10).
        // The exact reduced coefficient is checked once it is built, so a
        // spelling with more digits than its value needs is not rejected.
        let decades = if scale <= 0 {
            significant_digits - 1 + scale_digits
        } else if significant_digits > scale_digits {
            significant_digits - 1 - scale_digits
        } else {
            scale_digits - significant_digits
        };
        budget.bits(decades.checked_mul(3321).ok_or(Error::Limit(Resource::Bits))? / 1000)?;
        // Work follows the spelled-out numerator or power of ten, whichever
        // is longer (3.322 > log2 10).
        let spelled = if scale <= 0 {
            significant_digits + scale_digits
        } else {
            significant_digits.max(scale_digits)
        };
        let estimate = spelled
            .checked_mul(3322)
            .ok_or(Error::Limit(Resource::Bits))?
            / 1000
            + 1;
        budget.arithmetic(estimate)?;
        if significant_digits <= 19 && scale.unsigned_abs() <= 19 {
            let n = whole
                .bytes()
                .chain(fraction.bytes())
                .try_fold(0_i128, |n, b| {
                    n.checked_mul(10)?.checked_add((b - b'0') as i128)
                });
            if let Some(n) = n {
                let n = if negative { -n } else { n };
                let power = 10_i128.pow(scale.unsigned_abs() as u32);
                let pair = if scale >= 0 {
                    Some((n, power))
                } else {
                    n.checked_mul(power).map(|n| (n, 1))
                };
                if let Some((n, d)) = pair {
                    let r = Rational::wide(n, d, budget)?;
                    return self.number(r, false, budget);
                }
            }
        }
        let digits = format!("{whole}{fraction}");
        let digits = digits.trim_start_matches('0');
        let mut n = BigInt::parse_bytes(digits.as_bytes(), 10).ok_or_else(invalid)?;
        if negative {
            n = -n;
        }
        let magnitude =
            u32::try_from(scale.unsigned_abs()).map_err(|_| Error::Limit(Resource::Bits))?;
        let power = BigInt::from(10).pow(magnitude);
        let (n, d) = if scale >= 0 {
            (n, power)
        } else {
            (n * power, BigInt::one())
        };
        let r = Rational::from_big(n, d, budget)?;
        self.number(r, false, budget)
    }
}

fn power_factor(a: &Angle) -> (&Angle, i64) {
    match a.kind() {
        Some(NodeKind::Power(base, e)) => (base, *e),
        _ => (a, 1),
    }
}

fn trig_square(a: &Angle) -> Option<(Function, &Angle)> {
    let NodeKind::Power(base, 2) = a.kind()? else {
        return None;
    };
    let NodeKind::Function(fun, arg) = base.kind()? else {
        return None;
    };
    matches!(fun, Function::Sin | Function::Cos).then_some((*fun, arg))
}

fn sign_product(a: Sign, b: Sign) -> Sign {
    if a == Sign::Zero || b == Sign::Zero {
        return Sign::Zero;
    }
    match (a, b) {
        (Sign::Positive, x) | (x, Sign::Positive) => x,
        (Sign::Negative, x) | (x, Sign::Negative) => x.neg(),
        (Sign::NonNegative, Sign::NonNegative) | (Sign::NonPositive, Sign::NonPositive) => {
            Sign::NonNegative
        }
        (Sign::NonNegative, Sign::NonPositive) | (Sign::NonPositive, Sign::NonNegative) => {
            Sign::NonPositive
        }
        _ => Sign::Unknown,
    }
}
fn product_sign(coefficient: &Rational, factors: &[Angle]) -> Sign {
    factors
        .iter()
        .fold(Sign::from_integer(coefficient.signum()), |s, a| {
            sign_product(s, a.metadata().sign)
        })
}
fn sum_sign(terms: &[Term]) -> Sign {
    let mut positive = false;
    let mut negative = false;
    let mut strict_positive = false;
    let mut strict_negative = false;
    for term in terms {
        let sign = sign_product(
            Sign::from_integer(term.coefficient.signum()),
            term.atom.metadata().sign,
        );
        match sign {
            Sign::Unknown => return Sign::Unknown,
            Sign::Positive => {
                positive = true;
                strict_positive = true;
            }
            Sign::Negative => {
                negative = true;
                strict_negative = true;
            }
            Sign::NonNegative => positive = true,
            Sign::NonPositive => negative = true,
            Sign::Zero => {}
        }
    }
    match (positive, negative, strict_positive, strict_negative) {
        (true, false, true, _) => Sign::Positive,
        (false, true, _, true) => Sign::Negative,
        (true, false, _, _) => Sign::NonNegative,
        (false, true, _, _) => Sign::NonPositive,
        (false, false, _, _) => Sign::Zero,
        _ => Sign::Unknown,
    }
}
