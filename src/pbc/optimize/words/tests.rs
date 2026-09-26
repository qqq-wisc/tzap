use super::*;
use crate::pbc::{PauliAngle, PbcCircuit, Phase};
use rand::{Rng, SeedableRng, rngs::StdRng};

const PAULIS: [Pauli; 4] = [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z];

/// Pack a list of single-qubit factors (qubit i gets factors[i]).
fn packed(factors: &[Pauli], l: usize) -> Vec<u64> {
    let mut w = vec![0u64; 2 * l];
    for (q, p) in factors.iter().enumerate() {
        let bit = 1u64 << (q % 64);
        if matches!(p, Pauli::X | Pauli::Y) {
            w[q / 64] |= bit;
        }
        if matches!(p, Pauli::Z | Pauli::Y) {
            w[l + q / 64] |= bit;
        }
    }
    w
}

/// Product by the single-qubit table: (phase exponent, factors).
fn table_product(a: &[Pauli], b: &[Pauli]) -> (u8, Vec<Pauli>) {
    let mut phase = Phase::One;
    let factors = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| {
            let (p, f) = x.times(y);
            phase = phase * p;
            f
        })
        .collect();
    (phase as u8, factors)
}

#[test]
fn kernels_match_the_pauli_table_on_random_multiword_words() {
    let mut rng = StdRng::seed_from_u64(7);
    for n in [1usize, 2, 3, 63, 64, 65, 130] {
        let l = n.div_ceil(64);
        for _ in 0..200 {
            let a: Vec<_> = (0..n).map(|_| PAULIS[rng.gen_range(0..4)]).collect();
            let b: Vec<_> = (0..n).map(|_| PAULIS[rng.gen_range(0..4)]).collect();
            let (e, s) = table_product(&a, &b);
            let mut out = vec![0; 2 * l];
            assert_eq!(product(&packed(&a, l), &packed(&b, l), &mut out, l), e);
            assert_eq!(out, packed(&s, l));
            // Anticommuting exactly when the product phase is imaginary.
            assert_eq!(anticommutes(&packed(&a, l), &packed(&b, l), l), e % 2 == 1);
        }
    }
}

#[test]
fn exhaustive_two_qubit_products_including_xx_yy() {
    for a0 in PAULIS {
        for a1 in PAULIS {
            for b0 in PAULIS {
                for b1 in PAULIS {
                    let (a, b) = ([a0, a1], [b0, b1]);
                    let (e, s) = table_product(&a, &b);
                    let mut out = vec![0; 2];
                    assert_eq!(product(&packed(&a, 1), &packed(&b, 1), &mut out, 1), e);
                    assert_eq!(out, packed(&s, 1));
                }
            }
        }
    }
    // (XX)(YY) = -ZZ.
    let mut out = vec![0; 2];
    let e = product(
        &packed(&[Pauli::X, Pauli::X], 1),
        &packed(&[Pauli::Y, Pauli::Y], 1),
        &mut out,
        1,
    );
    assert_eq!((e, out), (2, packed(&[Pauli::Z, Pauli::Z], 1)));
}

#[test]
fn packing_matches_expand_including_signs_inside_the_dag() {
    let mut rng = StdRng::seed_from_u64(11);
    for n in [1usize, 3, 70] {
        let l = n.div_ceil(64);
        let mut c = PbcCircuit::new(n, 0);
        let mut pool = vec![c.identity().as_ref()];
        for _ in 0..60 {
            let r = if rng.gen_bool(0.4) {
                let q = rng.gen_range(0..n as u32);
                c.single(q, PAULIS[rng.gen_range(1..4)]).unwrap().as_ref()
            } else {
                let a = pool[rng.gen_range(0..pool.len())];
                let b = pool[rng.gen_range(0..pool.len())];
                c.product(a, b).unwrap()
            };
            pool.push(r);
            // Rotate about every Hermitian expression, fixing imaginary ones.
            let phase = c.expand(r, 1 << 20).unwrap().phase;
            let axis = match phase {
                Phase::I | Phase::MinusI => r.scaled(Phase::I),
                _ => r,
            };
            let axis = c.hermitian_axis(axis, 1 << 20).unwrap();
            c.rotate(axis, PauliAngle::new(rng.gen_range(0..8)))
                .unwrap();
        }
        let roots: Vec<_> = c.operations().iter().map(|op| op.axis().as_ref()).collect();
        let (axes, records) = pack(&c.arena, n, &roots, usize::MAX).unwrap();
        for (op, (id, sign)) in c.operations().iter().zip(records) {
            let expanded = c.expand(op.axis().as_ref(), 1 << 20).unwrap();
            assert_eq!(axes.get(id), packed(&expanded.factors, l));
            let expected_sign = if expanded.phase == Phase::One { 1 } else { -1 };
            assert_eq!(sign, expected_sign);
            // The canonical handle evaluates to +W.
            let handle = axes.handles[id as usize].unwrap();
            let canonical = c.expand(handle, 1 << 20).unwrap();
            assert_eq!(canonical.phase, Phase::One);
            assert_eq!(canonical.factors, expanded.factors);
        }
    }
}

#[test]
fn packing_respects_the_storage_budget() {
    let mut c = PbcCircuit::new(2, 0);
    let z = c.z(0).unwrap();
    c.rotate(z, PauliAngle::new(1)).unwrap();
    assert_eq!(
        pack(&c.arena, 2, &[z.as_ref()], 1).err(),
        Some(PbcError::ExpansionLimit)
    );
}
