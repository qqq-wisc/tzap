//! StateFold over programs: see `crate::state_fold::program`.

use super::{Program, Stmt};
use crate::state_fold::program::{fold as fold_located, PStmt};

fn locate(s: &Stmt, next: &mut usize) -> PStmt {
    match s {
        Stmt::Gate(g) => {
            *next += 1;
            PStmt::Gate(*next - 1, g.clone())
        }
        Stmt::Reset(q) => PStmt::Reset(*q),
        Stmt::Measure(q, bit) => PStmt::Measure(*q, bit.clone()),
        Stmt::Seq(xs) => PStmt::Seq(xs.iter().map(|x| locate(x, next)).collect()),
        Stmt::If(cond, a, b) => PStmt::If(cond.clone(), Box::new(locate(a, next)), Box::new(locate(b, next))),
        Stmt::While(cond, b) => PStmt::While(cond.clone(), Box::new(locate(b, next))),
    }
}

/// Feynman's `-qasm3 -statefold <d>` (`None`: unbounded degree).
pub fn fold(prog: &Program, degree: Option<usize>) -> Program {
    let located = locate(&prog.body, &mut 0);
    let body = super::parse::flatten(Stmt::Seq(fold_located(prog.num_qubits, degree, &located)));
    Program { num_qubits: prog.num_qubits, decls: prog.decls.clone(), body }
}
