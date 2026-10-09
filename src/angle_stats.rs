//! Per-run arithmetic diagnostics, propagated explicitly into rayon workers.
use crate::{
    angle::AngleError,
    circuit::{Circuit, Gate},
    optimize::{Options, PassName},
};
use std::{
    cell::RefCell,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NumericalReport {
    pub preserved_expressions: usize,
    pub numerical_fallbacks: usize,
    pub skipped_rounded_folds: usize,
    pub skipped_nonfinite_folds: usize,
    pub skipped_coefficient_limit_folds: usize,
    pub uncertified_syntheses: usize,
    pub randomized_matching: bool,
}
#[derive(Default)]
pub(crate) struct Counters {
    rounded: AtomicUsize,
    nonfinite: AtomicUsize,
    coefficient: AtomicUsize,
    synthesis: AtomicUsize,
}
thread_local! {static ACTIVE:RefCell<Option<Arc<Counters>>>=const{RefCell::new(None)};}
pub(crate) fn current() -> Option<Arc<Counters>> {
    ACTIVE.with(|a| a.borrow().clone())
}
pub(crate) struct Scope(Option<Arc<Counters>>);
impl Scope {
    pub(crate) fn install(counters: Option<Arc<Counters>>) -> Self {
        Self(ACTIVE.with(|a| a.replace(counters)))
    }
}
impl Drop for Scope {
    fn drop(&mut self) {
        ACTIVE.with(|a| a.replace(self.0.take()));
    }
}
pub(crate) fn failure(error: AngleError) {
    ACTIVE.with(|a| {
        if let Some(c) = a.borrow().as_ref() {
            let target = match error {
                AngleError::RoundedAddition => &c.rounded,
                AngleError::NonFinite => &c.nonfinite,
                AngleError::RepresentationLimit => &c.coefficient,
                _ => return,
            };
            target.fetch_add(1, Ordering::Relaxed);
        }
    });
}
pub(crate) fn synthesis() {
    ACTIVE.with(|a| {
        if let Some(c) = a.borrow().as_ref() {
            c.synthesis.fetch_add(1, Ordering::Relaxed);
        }
    });
}
impl NumericalReport {
    pub fn input(circuit: &Circuit, options: &Options) -> Self {
        let mut result = Self {
            randomized_matching: options
                .passes
                .as_ref()
                .map_or(true, |passes| passes.contains(&PassName::PhaseFoldRand)),
            ..Self::default()
        };
        for gate in &circuit.gates {
            if let Gate::rz(a, _) = gate {
                result.preserved_expressions += usize::from(a.is_preserved());
                result.numerical_fallbacks += usize::from(a.used_numeric_fallback());
            }
        }
        result
    }
    pub(crate) fn collected(circuit: &Circuit, options: &Options, c: &Counters) -> Self {
        Self {
            skipped_rounded_folds: c.rounded.load(Ordering::Relaxed),
            skipped_nonfinite_folds: c.nonfinite.load(Ordering::Relaxed),
            skipped_coefficient_limit_folds: c.coefficient.load(Ordering::Relaxed),
            uncertified_syntheses: c.synthesis.load(Ordering::Relaxed),
            ..Self::input(circuit, options)
        }
    }
}
