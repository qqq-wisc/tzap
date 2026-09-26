//! Machine-checked proof of `docs/pbc-weight-example.md`: every program in
//! that walkthrough is compared, with the exact-arithmetic oracle in
//! `crate::semantics`, against the input circuit's unitary (up to global
//! phase), and every quoted output, T count, and weight is pinned.

use super::*;
use crate::optimize::{Level, Options};
use crate::semantics::circuit_unitary;
use crate::semantics::test_support::assert_equivalent;

fn cx(control: u32, target: u32) -> Gate {
    Gate::cnot { control, target }
}

/// `cx q0,q1; t q1; t q1; tdg q1; cx q0,q1; h q0; cx q1,q0; t q0`.
fn input() -> Circuit {
    Circuit {
        num_qubits: 2,
        num_cbits: 0,
        gates: vec![
            cx(0, 1),
            Gate::t(1),
            Gate::t(1),
            Gate::tdg(1),
            cx(0, 1),
            Gate::h(0),
            cx(1, 0),
            Gate::t(0),
        ],
    }
}

/// The output Clifford C from converting the input: its Clifford gates, in
/// order.
const C: fn() -> Vec<Gate> = || vec![cx(0, 1), cx(0, 1), Gate::h(0), cx(1, 0)];

/// A row of the walk table as a program: rotations (processed, then F, then
/// the rest), followed by the output Clifford C.
fn row(layers: &[(&str, i64)]) -> PbcCircuit {
    let mut c = rotations(2, layers);
    for gate in C() {
        c.push_output_clifford(gate).unwrap();
    }
    c
}

fn weight(c: &PbcCircuit) -> usize {
    let text = c.to_text().unwrap();
    text.lines()
        .filter(|l| l.starts_with("r "))
        .map(|l| l.split_whitespace().count() - 3)
        .sum()
}

#[test]
fn step_1_conversion_is_exact() {
    let pbc = to_pbc(&input()).unwrap();
    assert_equivalent(&input(), &pbc);
    assert_eq!(
        pbc.to_text().unwrap(),
        "qubits 2\nregisters 0\nr 1 1 Z0 Z1\nr 1 1 Z0 Z1\nr -1 1 Z0 Z1\nr 1 1 X0 Z1\n\
         f X0 1 Z0\nf X1 1 Z0 X1\nf Z0 1 X0 Z1\n"
    );
    assert_eq!((pbc.t_count(), weight(&pbc)), (4, 8));
}

#[test]
fn pipeline_a_tzap_o1_then_conversion_is_exact() {
    let options = Options {
        level: Level::O1,
        ..Options::default()
    };
    let (optimized, _) = crate::optimize::optimize(&input(), &options).unwrap();
    assert_eq!(
        optimized.gates,
        vec![
            cx(0, 1),
            Gate::t(1),
            cx(0, 1),
            Gate::h(0),
            cx(1, 0),
            Gate::t(0)
        ]
    );
    // tzap's own rewrite preserves the unitary...
    let limits = Limits::default();
    assert!(
        circuit_unitary(&optimized, limits)
            .unwrap()
            .equivalent_up_to_global_phase(&circuit_unitary(&input(), limits).unwrap())
    );
    // ...and so does converting its output.
    let pbc = to_pbc(&optimized).unwrap();
    assert_equivalent(&input(), &pbc);
    assert_eq!(
        pbc.to_text().unwrap(),
        "qubits 2\nregisters 0\nr 1 1 Z0 Z1\nr 1 1 X0 Z1\n\
         f X0 1 Z0\nf X1 1 Z0 X1\nf Z0 1 X0 Z1\n"
    );
    assert_eq!((pbc.t_count(), weight(&pbc)), (2, 4));
}

#[test]
fn pipeline_b_conversion_then_the_pass_is_exact() {
    let mut pbc = to_pbc(&input()).unwrap();
    let stats = pbc.optimize_rotations(OptimizeOptions::default()).unwrap();
    assert_equivalent(&input(), &pbc);
    assert_eq!(
        pbc.to_text().unwrap(),
        "qubits 2\nregisters 0\nr -1 1 Z0 Z1\nr -1 1 Y0\n\
         f X0 1 Z0\nf X1 -1 Y1\nf Z0 -1 Y0\n"
    );
    assert_eq!((pbc.t_count(), weight(&pbc)), (2, 3));
    assert_eq!((stats.merges, stats.cliffords_to_frame), (1, 1));
}

/// Every row of the walk table, as "processed ; F ; not yet processed ; C",
/// is the input's unitary. F = S(Z0Z1) is written as the rotation `r 2` about
/// Z0Z1; T(-Y0) as `r 1` about -Y0.
#[test]
fn every_row_of_the_walk_table_is_exact() {
    let rows: [&[(&str, i64)]; 5] = [
        // start: F is the identity; everything still to process.
        &[("ZZ", 1), ("ZZ", 1), ("ZZ", -1), ("XZ", 1)],
        // step 1: T(Z0Z1) processed.
        &[("ZZ", 1), ("ZZ", 1), ("ZZ", -1), ("XZ", 1)],
        // step 2: merged into F = S(Z0Z1).
        &[("ZZ", 2), ("ZZ", -1), ("XZ", 1)],
        // step 3: T†(Z0Z1) processed; F moved past it, unchanged.
        &[("ZZ", -1), ("ZZ", 2), ("XZ", 1)],
        // step 4: T(X0Z1) processed; F moved past it, axis now -Y0.
        &[("ZZ", -1), ("-YI", 1), ("ZZ", 2)],
    ];
    for layers in rows {
        assert_equivalent(&input(), &row(layers));
    }
    // Row "end": F merged into C. Rebuild C' = S(Z0Z1) then C by the gates
    // (S on Z0Z1 is CX; S on q1; CX), and check it matches the output frame
    // quoted in the doc.
    let mut end = rotations(2, &[("ZZ", -1), ("-YI", 1)]);
    for gate in [cx(0, 1), Gate::s(1), cx(0, 1)].into_iter().chain(C()) {
        end.push_output_clifford(gate).unwrap();
    }
    assert_equivalent(&input(), &end);
    assert_eq!(
        end.to_text().unwrap(),
        "qubits 2\nregisters 0\nr -1 1 Z0 Z1\nr 1 -1 Y0\n\
         f X0 1 Z0\nf X1 -1 Y1\nf Z0 -1 Y0\n"
    );
}

/// The swapped axis in step 4 and the frame table: moving S(Z0Z1) past each
/// Pauli as claimed. Checked on unitaries: S Q' = Q S, i.e. running S then
/// a rotation about Q equals running a rotation about Q' then S.
#[test]
fn moving_s_past_each_axis_is_exact() {
    let limits = Limits::default();
    let unitary = |layers: &[(&str, i64)]| pbc_unitary(&rotations(2, layers), limits).unwrap();
    for (q, q_after) in [
        ("ZZ", "ZZ"),  // T†(Z0Z1): commutes, unchanged
        ("XZ", "-YI"), // T(X0Z1) -> T(-Y0)
        ("ZI", "ZI"),  // frame X0: Z0 unchanged
        ("ZX", "-IY"), // frame X1: Z0X1 -> -Y1
        ("IZ", "IZ"),  // frame Z1: unchanged
    ] {
        let before = unitary(&[("ZZ", 2), (q, 1)]);
        let after = unitary(&[(q_after, 1), ("ZZ", 2)]);
        assert!(
            before.equivalent_up_to_global_phase(&after),
            "{q} -> {q_after}"
        );
    }
}
