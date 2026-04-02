use indicatif::ProgressBar;

use crate::circuit::{Circuit, Gate};

pub trait Pass: Sync {
    fn name(&self) -> &str;
    fn run(&self, circuit: &Circuit) -> Circuit {
        self.run_with_progress(circuit, &ProgressBar::hidden())
    }
    fn run_with_progress(&self, circuit: &Circuit, pb: &ProgressBar) -> Circuit;
}

pub struct PassResult {
    pub circuit: Circuit,
    pub t_after_first: usize,
    pub gates_after_first: usize,
}

pub fn run_passes(circuit: &Circuit, passes: &[&dyn Pass]) -> PassResult {
    let mut c = circuit.clone();
    let mut t_after_first = 0;
    let mut gates_after_first = 0;
    for (i, p) in passes.iter().enumerate() {
        c = p.run(&c);
        if i == 0 {
            t_after_first = count_t(&c);
            gates_after_first = c.gates.len();
        }
    }
    PassResult { circuit: c, t_after_first, gates_after_first }
}

pub fn count_t(c: &Circuit) -> usize {
    c.gates.iter().filter(|g| matches!(g, Gate::t(_) | Gate::tdg(_))).count()
}

pub fn count_rz(c: &Circuit) -> usize {
    c.gates.iter().filter(|g| matches!(g, Gate::rz(..))).count()
}

/// Runs a sequence of passes N times.
pub struct Repeat<'a> {
    pub passes: &'a [&'a dyn Pass],
    pub times: usize,
}

impl Pass for Repeat<'_> {
    fn name(&self) -> &str {
        "repeat"
    }

    fn run_with_progress(&self, circuit: &Circuit, _pb: &ProgressBar) -> Circuit {
        let mut c = circuit.clone();
        for i in 0..self.times {
            eprintln!("  round {}/{}", i + 1, self.times);
            c = run_passes(&c, self.passes).circuit;
        }
        c
    }
}
