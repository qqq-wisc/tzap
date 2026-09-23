//! Shared assertions for the independent gate and PBC reference interpreters.
use super::{Limits, Matrix, circuit_unitary, pbc_unitary, scalar::Scalar};
use crate::{circuit::Circuit, pbc::PbcCircuit};

pub(crate) fn assert_equivalent(input: &Circuit, pbc: &PbcCircuit) {
    let actual = pbc_unitary(pbc, Limits::default()).expect("PBC unitary within test limits");
    let expected =
        circuit_unitary(input, Limits::default()).expect("gate unitary within test limits");
    assert!(
        actual.equivalent_up_to_global_phase(&expected),
        "PBC mismatch for input:\n{input}"
    );
    check_numeric(&actual, input);
    assert_eq!(actual.adjoint().mul(&actual), Matrix::identity(actual.dim));
}

/// Compare to the existing floating-point interpreter using a largest-magnitude
/// pivot to align phase. Exact equality above remains the primary oracle.
pub(crate) fn check_numeric(actual: &Matrix, circuit: &Circuit) {
    let expected = crate::unitary::circuit_unitary(circuit);
    let (r, c) = (0..actual.dim)
        .flat_map(|r| (0..actual.dim).map(move |c| (r, c)))
        .max_by(|&(r, c), &(s, t)| {
            expected[r][c]
                .norm_sq()
                .total_cmp(&expected[s][t].norm_sq())
        })
        .unwrap();
    let a = approximate(actual.get(r, c));
    let b = expected[r][c].components();
    assert!((a.0 * a.0 + a.1 * a.1 - b.0 * b.0 - b.1 * b.1).abs() < 1e-10);
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    for (r, row) in expected.iter().enumerate() {
        for (c, value) in row.iter().enumerate() {
            let left = mul(approximate(actual.get(r, c)), b);
            let right = mul(value.components(), a);
            assert!(
                (left.0 - right.0).abs() < 1e-10 && (left.1 - right.1).abs() < 1e-10,
                "numeric mismatch at ({r}, {c}): {left:?} != {right:?}\n{circuit}"
            );
        }
    }
}

pub(crate) fn approximate(s: &Scalar) -> (f64, f64) {
    let c: Vec<f64> =
        s.0.iter()
            .map(|r| {
                r.numer().to_string().parse::<f64>().unwrap()
                    / r.denom().to_string().parse::<f64>().unwrap()
            })
            .collect();
    (c[0] + 2.0_f64.sqrt() * c[1], c[2] + 2.0_f64.sqrt() * c[3])
}
