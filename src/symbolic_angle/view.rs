//! Borrowed scalar structure for enclosing evaluation and serialization.
//! Does not expose interner identity, node keys, or fingerprints.
use super::*;
use std::borrow::Cow;

pub(crate) enum ExpressionRef<'a> {
    Rational(Cow<'a, Rational>),
    Pi(Cow<'a, Rational>),
    Sum(SumTerms<'a>),
    Product {
        coefficient: &'a Rational,
        factors: &'a [Angle],
    },
    Power(&'a Angle, i64),
    Function(Function, &'a Angle),
}
pub(crate) struct SumTerms<'a>(std::slice::Iter<'a, Term>);
impl<'a> Iterator for SumTerms<'a> {
    type Item = (&'a Rational, &'a Angle);
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|term| (&term.coefficient, &term.atom))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}
impl ExactSizeIterator for SumTerms<'_> {}

impl Angle {
    pub(crate) fn expression(&self) -> ExpressionRef<'_> {
        match &self.0 {
            Value::Rational { n, d } => {
                ExpressionRef::Rational(Cow::Owned(Rational::canonical(*n as i64, *d as u64)))
            }
            Value::Pi { n, d } => {
                ExpressionRef::Pi(Cow::Owned(Rational::canonical(*n as i64, *d as u64)))
            }
            Value::Dyadic(bits) => ExpressionRef::Rational(Cow::Owned(Rational::from_bits(*bits))),
            Value::Shared(node) => match &node.kind {
                NodeKind::Rational(r) => ExpressionRef::Rational(Cow::Borrowed(r)),
                NodeKind::Pi(r) => ExpressionRef::Pi(Cow::Borrowed(r)),
                NodeKind::Sum(terms) => ExpressionRef::Sum(SumTerms(terms.iter())),
                NodeKind::Product {
                    coefficient,
                    factors,
                } => ExpressionRef::Product {
                    coefficient,
                    factors,
                },
                NodeKind::Power(a, e) => ExpressionRef::Power(a, *e),
                NodeKind::Function(fun, a) => ExpressionRef::Function(*fun, a),
            },
        }
    }
}
