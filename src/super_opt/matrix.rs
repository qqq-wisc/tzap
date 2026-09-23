//! Dense exact unitary matrices for Clifford+T window semantics.
//!
//! Each entry is `(a + b*omega + c*omega^2 + d*omega^3) / sqrt(2)^k`, where
//! `omega = exp(i*pi/4)`, `omega^4 = -1`, and one denominator exponent is
//! shared by the whole matrix. Numerators use four symmetric `i8` coefficients
//! (`-127..=127`); arithmetic widens to `i16`, and an unrepresentable window or
//! MURM child is conservatively skipped.
//!
//! Hadamards add and subtract row pairs and increment `k`; phase gates multiply
//! rows by powers of `omega`; controlled-X gates permute rows. Common factors of
//! `sqrt(2)` are removed exactly after each Hadamard. Equality and deterministic
//! 64-bit fingerprints canonicalize the eight possible Clifford+T global
//! phases, and every fingerprint hit is confirmed by exact matrix comparison.
//! Rz gates are window barriers and never reach this code.

use crate::circuit::{Gate, Qubit};

use super::SuperOptError;

#[cfg(test)]
pub(super) const IDENTITY_TOLERANCE: f64 = 1e-10;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct CoefficientOverflow;

fn narrow_coefficient(value: i16) -> Option<i8> {
    let value = i8::try_from(value).ok()?;
    (value != i8::MIN).then_some(value)
}

/// `a + b*omega + c*omega^2 + d*omega^3`, with `omega^4 = -1`.
///
/// The symmetric range `-127..=127` is maintained so multiplying by any power
/// of `omega` (a signed permutation) remains representable. Arithmetic widens
/// to `i16` and reports overflow before narrowing back to `i8`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Cyclotomic {
    coefficients: [i8; 4],
}

impl Cyclotomic {
    const ZERO: Self = Self {
        coefficients: [0; 4],
    };
    const ONE: Self = Self {
        coefficients: [1, 0, 0, 0],
    };

    fn is_zero(self) -> bool {
        self == Self::ZERO
    }

    fn checked_add(self, rhs: Self) -> Option<Self> {
        self.checked_zip(rhs, i16::checked_add)
    }

    fn checked_sub(self, rhs: Self) -> Option<Self> {
        self.checked_zip(rhs, i16::checked_sub)
    }

    fn checked_zip(self, rhs: Self, op: fn(i16, i16) -> Option<i16>) -> Option<Self> {
        let mut coefficients = [0; 4];
        for (index, output) in coefficients.iter_mut().enumerate() {
            let value = op(
                i16::from(self.coefficients[index]),
                i16::from(rhs.coefficients[index]),
            )?;
            *output = narrow_coefficient(value)?;
        }
        Some(Self { coefficients })
    }

    /// Multiply by `omega^power` using `omega^4 = -1`.
    fn times_omega(self, power: u8) -> Self {
        let [a, b, c, d] = self.coefficients;
        debug_assert!(self.coefficients.iter().all(|&value| value != i8::MIN));
        let neg = |value: i8| -value;
        let coefficients = match power & 7 {
            0 => [a, b, c, d],
            1 => [neg(d), a, b, c],
            2 => [neg(c), neg(d), a, b],
            3 => [neg(b), neg(c), neg(d), a],
            4 => [neg(a), neg(b), neg(c), neg(d)],
            5 => [d, neg(a), neg(b), neg(c)],
            6 => [c, d, neg(a), neg(b)],
            7 => [b, c, d, neg(a)],
            _ => unreachable!(),
        };
        Self { coefficients }
    }

    /// Whether this numerator is divisible by `sqrt(2) = omega - omega^3`.
    fn divisible_by_sqrt_2(self) -> bool {
        let [a, b, c, d] = self.coefficients;
        (a & 1) == (c & 1) && (b & 1) == (d & 1)
    }

    /// Divide an exactly divisible numerator by `sqrt(2)`.
    fn divide_by_sqrt_2(self) -> Option<Self> {
        debug_assert!(self.divisible_by_sqrt_2());
        let [a, b, c, d] = self.coefficients.map(i16::from);
        let values = [(b - d) / 2, (a + c) / 2, (b + d) / 2, (c - a) / 2];
        let mut coefficients = [0; 4];
        for (output, value) in coefficients.iter_mut().zip(values) {
            *output = narrow_coefficient(value)?;
        }
        Some(Self { coefficients })
    }

    #[cfg(test)]
    fn to_cartesian(self, denominator_exponent: u16) -> (f64, f64) {
        let [a, b, c, d] = self.coefficients.map(f64::from);
        let root_half = std::f64::consts::FRAC_1_SQRT_2;
        let denominator = root_half.powi(i32::from(denominator_exponent));
        (
            (a + (b - d) * root_half) * denominator,
            (c + (b + d) * root_half) * denominator,
        )
    }
}

/// A dense `2^n` by `2^n` exact unitary matrix in row-major order.
///
/// Qubit zero is the most significant basis-state bit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnitaryMatrix {
    num_qubits: usize,
    dim: usize,
    denominator_exponent: u16,
    data: Vec<Cyclotomic>,
}

impl UnitaryMatrix {
    pub(super) fn identity(num_qubits: usize) -> Result<Self, SuperOptError> {
        if num_qubits >= usize::BITS as usize / 2 {
            return Err(SuperOptError::MatrixTooLarge { num_qubits });
        }
        let dim = 1 << num_qubits;
        let mut data = vec![Cyclotomic::ZERO; dim * dim];
        for i in 0..dim {
            data[i * dim + i] = Cyclotomic::ONE;
        }
        Ok(Self {
            num_qubits,
            dim,
            denominator_exponent: 0,
            data,
        })
    }

    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    #[cfg(test)]
    pub(super) fn entry_coefficients(&self, row: usize, column: usize) -> ([i8; 4], u16) {
        (
            self.data[row * self.dim + column].coefficients,
            self.denominator_exponent,
        )
    }

    pub(super) fn equivalent_up_to_global_phase(&self, other: &Self) -> bool {
        if self.num_qubits != other.num_qubits
            || self.denominator_exponent != other.denominator_exponent
        {
            return false;
        }
        let left_phase = canonical_phase_power(self);
        let right_phase = canonical_phase_power(other);
        self.data
            .iter()
            .zip(&other.data)
            .all(|(&left, &right)| left.times_omega(left_phase) == right.times_omega(right_phase))
    }

    pub(super) fn copy_from(&mut self, source: &Self) {
        self.num_qubits = source.num_qubits;
        self.dim = source.dim;
        self.denominator_exponent = source.denominator_exponent;
        if self.data.len() == source.data.len() {
            self.data.copy_from_slice(&source.data);
        } else {
            self.data.clone_from(&source.data);
        }
    }

    fn row_pair_mut(&mut self, row0: usize, row1: usize) -> (&mut [Cyclotomic], &mut [Cyclotomic]) {
        let dim = self.dim;
        let (head, tail) = self.data.split_at_mut(row1 * dim);
        (&mut head[row0 * dim..(row0 + 1) * dim], &mut tail[..dim])
    }

    fn apply_h_left(&mut self, bit: usize) -> Result<(), CoefficientOverflow> {
        self.denominator_exponent = self
            .denominator_exponent
            .checked_add(1)
            .ok_or(CoefficientOverflow)?;
        for row0 in 0..self.dim {
            if row0 & bit != 0 {
                continue;
            }
            let (top, bottom) = self.row_pair_mut(row0, row0 | bit);
            for (a, b) in top.iter_mut().zip(bottom) {
                let (x, y) = (*a, *b);
                let sum = x.checked_add(y).ok_or(CoefficientOverflow)?;
                let difference = x.checked_sub(y).ok_or(CoefficientOverflow)?;
                (*a, *b) = (sum, difference);
            }
        }
        self.normalize_denominator()?;
        Ok(())
    }

    fn normalize_denominator(&mut self) -> Result<(), CoefficientOverflow> {
        while self.denominator_exponent > 0
            && self.data.iter().all(|entry| entry.divisible_by_sqrt_2())
        {
            for entry in &mut self.data {
                *entry = entry.divide_by_sqrt_2().ok_or(CoefficientOverflow)?;
            }
            self.denominator_exponent -= 1;
        }
        Ok(())
    }

    /// Multiply target-one rows by a power of `omega`.
    fn apply_phase_left(&mut self, bit: usize, power: u8) {
        for row in 0..self.dim {
            if row & bit == 0 {
                continue;
            }
            for value in &mut self.data[row * self.dim..(row + 1) * self.dim] {
                *value = value.times_omega(power);
            }
        }
    }

    fn apply_controlled_x_left(&mut self, control_mask: usize, target_bit: usize) {
        for row0 in 0..self.dim {
            if row0 & control_mask == control_mask && row0 & target_bit == 0 {
                let (top, bottom) = self.row_pair_mut(row0, row0 | target_bit);
                top.swap_with_slice(bottom);
            }
        }
    }

    fn apply_phase_flip_left(&mut self, mask: usize) {
        for row in 0..self.dim {
            if row & mask == mask {
                for value in &mut self.data[row * self.dim..(row + 1) * self.dim] {
                    *value = value.times_omega(4);
                }
            }
        }
    }

    pub(super) fn apply_gate_left(
        &mut self,
        gate: &Gate,
        support: &[Qubit],
    ) -> Result<(), CoefficientOverflow> {
        let bit = |q: &Qubit| {
            let local = support.binary_search(q).expect("gate qubit is in support");
            qubit_bit(self.num_qubits, local)
        };
        match gate {
            Gate::x(q) => self.apply_controlled_x_left(0, bit(q)),
            Gate::h(q) => return self.apply_h_left(bit(q)),
            Gate::s(q) => self.apply_phase_left(bit(q), 2),
            Gate::sdg(q) => self.apply_phase_left(bit(q), 6),
            Gate::z(q) => self.apply_phase_flip_left(bit(q)),
            Gate::t(q) => self.apply_phase_left(bit(q), 1),
            Gate::tdg(q) => self.apply_phase_left(bit(q), 7),
            Gate::rz(..) => unreachable!("Rz windows are rejected before matrix construction"),
            Gate::cnot { control, target } => {
                self.apply_controlled_x_left(bit(control), bit(target))
            }
            Gate::cz { control, target } => self.apply_phase_flip_left(bit(control) | bit(target)),
            Gate::ccx {
                control1,
                control2,
                target,
            } => self.apply_controlled_x_left(bit(control1) | bit(control2), bit(target)),
            Gate::ccz {
                control1,
                control2,
                target,
            } => self.apply_phase_flip_left(bit(control1) | bit(control2) | bit(target)),
            Gate::measure { .. } | Gate::reset(_) => {
                unreachable!("measurement and reset are window barriers")
            }
        }
        Ok(())
    }
}

/// Pick the power of `omega` that makes the first nonzero entry's exact
/// coefficient tuple lexicographically smallest. This canonicalizes the
/// eight possible Clifford+T global phases without division or rounding.
fn canonical_phase_power(matrix: &UnitaryMatrix) -> u8 {
    let Some(&pivot) = matrix.data.iter().find(|entry| !entry.is_zero()) else {
        return 0;
    };
    (0..8)
        .min_by_key(|&power| pivot.times_omega(power).coefficients)
        .expect("the phase candidate range is nonempty")
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub(super) struct UnitaryFingerprint(u64);

impl UnitaryFingerprint {
    pub(super) fn to_bits(self) -> u64 {
        self.0
    }

    pub(super) fn from_bits(bits: u64) -> Self {
        Self(bits)
    }
}

/// A deterministic 64-bit hash of the exact phase-canonical matrix.
/// Hash collisions are still confirmed by exact matrix comparison before a
/// rewrite is accepted.
pub(super) fn unitary_fingerprint(matrix: &UnitaryMatrix) -> UnitaryFingerprint {
    let phase = canonical_phase_power(matrix);
    let mut first = 0xcbf2_9ce4_8422_2325u64;
    let mut second = 0x9e37_79b9_7f4a_7c15u64;
    let mut mix = |word: u64| {
        first ^= word;
        first = first.wrapping_mul(0x0000_0100_0000_01b3);
        second ^= word.wrapping_add(0x517c_c1b7_2722_0a95);
        second = second.rotate_left(27).wrapping_mul(0x94d0_49bb_1331_11eb);
    };
    mix(matrix.num_qubits as u64);
    mix(u64::from(matrix.denominator_exponent));
    for &entry in &matrix.data {
        for coefficient in entry.times_omega(phase).coefficients {
            mix(i64::from(coefficient) as u64);
        }
    }
    // Fold the independently mixed streams, then avalanche so every output
    // bit depends on both halves. Keeping both streams preserves the previous
    // mixer's distribution while halving the stored and in-memory key width.
    let mut fingerprint = first ^ second.rotate_left(32);
    fingerprint ^= fingerprint >> 30;
    fingerprint = fingerprint.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    fingerprint ^= fingerprint >> 27;
    fingerprint = fingerprint.wrapping_mul(0x94d0_49bb_1331_11eb);
    fingerprint ^= fingerprint >> 31;
    UnitaryFingerprint(fingerprint)
}

fn qubit_bit(num_qubits: usize, position: usize) -> usize {
    1usize << (num_qubits - 1 - position)
}

#[cfg(test)]
mod exact_tests {
    use std::collections::HashMap;
    use std::f64::consts::{FRAC_PI_4, TAU};

    use super::*;
    use crate::circuit::Circuit;

    const FLOAT_TOLERANCE: f64 = 2e-10;

    fn exact_matrix(circuit: &Circuit) -> UnitaryMatrix {
        let support: Vec<Qubit> = (0..circuit.num_qubits as Qubit).collect();
        let mut matrix = UnitaryMatrix::identity(circuit.num_qubits).unwrap();
        for gate in &circuit.gates {
            matrix.apply_gate_left(gate, &support).unwrap();
            assert_normalized(&matrix);
        }
        matrix
    }

    fn assert_normalized(matrix: &UnitaryMatrix) {
        assert!(
            matrix.denominator_exponent == 0
                || !matrix.data.iter().all(|entry| entry.divisible_by_sqrt_2()),
            "matrix retained a common sqrt(2) factor at denominator exponent {}",
            matrix.denominator_exponent
        );
    }

    fn assert_exactly_equal(left: &UnitaryMatrix, right: &UnitaryMatrix) {
        assert_eq!(left.num_qubits, right.num_qubits);
        assert_eq!(left.denominator_exponent, right.denominator_exponent);
        assert_eq!(left.data, right.data);
    }

    fn assert_matches_floating_oracle(circuit: &Circuit) {
        let exact = exact_matrix(circuit);
        let floating = crate::unitary::circuit_unitary(circuit);
        let dim = 1 << circuit.num_qubits;
        for (row, floating_row) in floating.iter().enumerate().take(dim) {
            for (column, &expected) in floating_row.iter().enumerate().take(dim) {
                let actual =
                    exact.data[row * dim + column].to_cartesian(exact.denominator_exponent);
                let (expected_re, expected_im) = expected.components();
                let error = (actual.0 - expected_re).hypot(actual.1 - expected_im);
                assert!(
                    error <= FLOAT_TOLERANCE,
                    "entry ({row}, {column}) differs after {:?}: exact={actual:?}, floating=({expected_re}, {expected_im}), error={error}",
                    circuit.gates
                );
            }
        }
    }

    fn circuit(num_qubits: usize, gates: &[Gate]) -> Circuit {
        let mut circuit = Circuit::new(num_qubits);
        for gate in gates {
            circuit.apply(gate.clone());
        }
        circuit
    }

    #[test]
    fn omega_powers_match_eighth_roots_of_unity() {
        for power in 0..8 {
            let exact = Cyclotomic::ONE.times_omega(power).to_cartesian(0);
            let angle = f64::from(power) * FRAC_PI_4;
            assert!((exact.0 - angle.cos()).abs() < 1e-15, "power {power}");
            assert!((exact.1 - angle.sin()).abs() < 1e-15, "power {power}");
        }
    }

    #[test]
    fn omega_multiplication_obeys_the_cyclic_group_law() {
        let samples = [
            Cyclotomic::ZERO,
            Cyclotomic::ONE,
            Cyclotomic {
                coefficients: [2, -3, 5, -7],
            },
            Cyclotomic {
                coefficients: [-11, 13, -17, 19],
            },
        ];
        for value in samples {
            assert_eq!(value.times_omega(8), value);
            for left in 0..8 {
                for right in 0..8 {
                    assert_eq!(
                        value.times_omega(left).times_omega(right),
                        value.times_omega((left + right) & 7)
                    );
                }
            }
        }
    }

    #[test]
    fn sqrt_two_division_round_trips_many_coefficients() {
        for a in -2..=2 {
            for b in -2..=2 {
                for c in -2..=2 {
                    for d in -2..=2 {
                        let value = Cyclotomic {
                            coefficients: [a, b, c, d],
                        };
                        // sqrt(2) = omega - omega^3.
                        let multiplied = value
                            .times_omega(1)
                            .checked_sub(value.times_omega(3))
                            .unwrap();
                        assert!(multiplied.divisible_by_sqrt_2());
                        assert_eq!(multiplied.divide_by_sqrt_2().unwrap(), value);
                    }
                }
            }
        }
    }

    #[test]
    fn sqrt_two_divisibility_is_exactly_the_parity_condition() {
        for a in -3..=3 {
            for b in -3..=3 {
                for c in -3..=3 {
                    for d in -3..=3 {
                        let value = Cyclotomic {
                            coefficients: [a, b, c, d],
                        };
                        assert_eq!(
                            value.divisible_by_sqrt_2(),
                            (a & 1) == (c & 1) && (b & 1) == (d & 1)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn coefficient_bound_rejects_both_unrepresentable_endpoints() {
        let positive_limit = Cyclotomic {
            coefficients: [127, 0, 0, 0],
        };
        let negative_limit = Cyclotomic {
            coefficients: [-127, 0, 0, 0],
        };
        assert_eq!(positive_limit.checked_add(Cyclotomic::ONE), None);
        assert_eq!(negative_limit.checked_sub(Cyclotomic::ONE), None);
    }

    #[test]
    fn exact_single_qubit_gate_relations_hold() {
        let identity = exact_matrix(&Circuit::new(1));
        let cases: &[(&[Gate], &[Gate])] = &[
            (&[Gate::h(0), Gate::h(0)], &[]),
            (&[Gate::x(0), Gate::x(0)], &[]),
            (&[Gate::z(0), Gate::z(0)], &[]),
            (&[Gate::s(0), Gate::sdg(0)], &[]),
            (&[Gate::t(0), Gate::tdg(0)], &[]),
            (&[Gate::t(0), Gate::t(0)], &[Gate::s(0)]),
            (&[Gate::s(0), Gate::s(0)], &[Gate::z(0)]),
            (
                &[Gate::t(0), Gate::t(0), Gate::t(0), Gate::t(0)],
                &[Gate::z(0)],
            ),
        ];
        for &(left, right) in cases {
            assert_exactly_equal(
                &exact_matrix(&circuit(1, left)),
                &exact_matrix(&circuit(1, right)),
            );
        }
        let t_eight = circuit(1, &vec![Gate::t(0); 8]);
        assert_exactly_equal(&exact_matrix(&t_eight), &identity);
    }

    #[test]
    fn exact_controlled_gate_relations_hold() {
        for gates in [
            vec![
                Gate::cnot {
                    control: 0,
                    target: 1,
                };
                2
            ],
            vec![
                Gate::cz {
                    control: 0,
                    target: 1,
                };
                2
            ],
        ] {
            assert_exactly_equal(
                &exact_matrix(&circuit(2, &gates)),
                &exact_matrix(&Circuit::new(2)),
            );
        }
        for gates in [
            vec![
                Gate::ccx {
                    control1: 0,
                    control2: 1,
                    target: 2,
                };
                2
            ],
            vec![
                Gate::ccz {
                    control1: 0,
                    control2: 1,
                    target: 2,
                };
                2
            ],
        ] {
            assert_exactly_equal(
                &exact_matrix(&circuit(3, &gates)),
                &exact_matrix(&Circuit::new(3)),
            );
        }
    }

    #[test]
    fn denominator_normalization_handles_nested_hadamards() {
        for num_qubits in 1..=5 {
            let mut gates = Vec::new();
            for qubit in 0..num_qubits as Qubit {
                gates.push(Gate::h(qubit));
            }
            for qubit in (0..num_qubits as Qubit).rev() {
                gates.push(Gate::h(qubit));
            }
            let matrix = exact_matrix(&circuit(num_qubits, &gates));
            assert_eq!(matrix.denominator_exponent, 0);
            assert_exactly_equal(&matrix, &exact_matrix(&Circuit::new(num_qubits)));
        }
    }

    #[test]
    fn all_eighth_root_global_phases_have_one_fingerprint() {
        let identity = exact_matrix(&Circuit::new(1));
        let identity_fingerprint = unitary_fingerprint(&identity);
        let omega_identity = [Gate::x(0), Gate::t(0), Gate::x(0), Gate::t(0)];
        let mut gates = Vec::new();
        for power in 0..8 {
            let phased = exact_matrix(&circuit(1, &gates));
            assert!(
                identity.equivalent_up_to_global_phase(&phased),
                "power {power}"
            );
            assert_eq!(unitary_fingerprint(&phased), identity_fingerprint);
            gates.extend_from_slice(&omega_identity);
        }
    }

    #[test]
    fn distinct_one_gate_projective_unitaries_have_distinct_fingerprints() {
        let gates = [
            Gate::x(0),
            Gate::h(0),
            Gate::s(0),
            Gate::sdg(0),
            Gate::z(0),
            Gate::t(0),
            Gate::tdg(0),
        ];
        let mut seen = HashMap::new();
        for gate in gates {
            let matrix = exact_matrix(&circuit(1, std::slice::from_ref(&gate)));
            assert!(
                seen.insert(unitary_fingerprint(&matrix), gate).is_none(),
                "distinct one-gate matrices shared a fingerprint"
            );
        }
    }

    #[test]
    fn fingerprint_occupies_eight_bytes() {
        assert_eq!(std::mem::size_of::<UnitaryFingerprint>(), 8);
    }

    #[test]
    fn exhaustive_single_qubit_circuits_match_floating_oracle() {
        fn visit(
            circuit: &mut Circuit,
            gates: &[Gate],
            remaining: usize,
            fingerprints: &mut HashMap<UnitaryFingerprint, UnitaryMatrix>,
            checked: &mut usize,
        ) {
            assert_matches_floating_oracle(circuit);
            let matrix = exact_matrix(circuit);
            let fingerprint = unitary_fingerprint(&matrix);
            if let Some(previous) = fingerprints.get(&fingerprint) {
                assert!(matrix.equivalent_up_to_global_phase(previous));
            } else {
                fingerprints.insert(fingerprint, matrix);
            }
            *checked += 1;
            if remaining == 0 {
                return;
            }
            for gate in gates {
                circuit.apply(gate.clone());
                visit(circuit, gates, remaining - 1, fingerprints, checked);
                circuit.gates.pop();
            }
        }

        let gates = [
            Gate::x(0),
            Gate::h(0),
            Gate::s(0),
            Gate::sdg(0),
            Gate::z(0),
            Gate::t(0),
            Gate::tdg(0),
        ];
        let mut checked = 0;
        visit(
            &mut Circuit::new(1),
            &gates,
            5,
            &mut HashMap::new(),
            &mut checked,
        );
        assert_eq!(checked, 19_608);
    }

    #[test]
    fn exhaustive_two_qubit_circuits_match_floating_oracle() {
        fn visit(circuit: &mut Circuit, gates: &[Gate], remaining: usize, checked: &mut usize) {
            assert_matches_floating_oracle(circuit);
            *checked += 1;
            if remaining == 0 {
                return;
            }
            for gate in gates {
                circuit.apply(gate.clone());
                visit(circuit, gates, remaining - 1, checked);
                circuit.gates.pop();
            }
        }

        let gates = [
            Gate::h(0),
            Gate::h(1),
            Gate::t(0),
            Gate::t(1),
            Gate::s(0),
            Gate::s(1),
            Gate::x(0),
            Gate::x(1),
            Gate::cnot {
                control: 0,
                target: 1,
            },
            Gate::cnot {
                control: 1,
                target: 0,
            },
            Gate::cz {
                control: 0,
                target: 1,
            },
        ];
        let mut checked = 0;
        visit(&mut Circuit::new(2), &gates, 4, &mut checked);
        assert_eq!(checked, 16_105);
    }

    struct TestRng(u64);

    impl TestRng {
        fn next(&mut self, upper: usize) -> usize {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0 as usize % upper
        }

        /// [`TestRng::next`] as a qubit index.
        fn qubit(&mut self, upper: usize) -> Qubit {
            self.next(upper) as Qubit
        }
    }

    fn random_distinct_qubits(rng: &mut TestRng, num_qubits: usize, count: usize) -> Vec<Qubit> {
        let mut qubits = Vec::with_capacity(count);
        while qubits.len() < count {
            let qubit = rng.qubit(num_qubits);
            if !qubits.contains(&qubit) {
                qubits.push(qubit);
            }
        }
        qubits
    }

    fn random_gate(rng: &mut TestRng, num_qubits: usize) -> Gate {
        let choices = if num_qubits >= 3 {
            11
        } else if num_qubits == 2 {
            9
        } else {
            7
        };
        match rng.next(choices) {
            0 => Gate::x(rng.qubit(num_qubits)),
            1 => Gate::h(rng.qubit(num_qubits)),
            2 => Gate::s(rng.qubit(num_qubits)),
            3 => Gate::sdg(rng.qubit(num_qubits)),
            4 => Gate::z(rng.qubit(num_qubits)),
            5 => Gate::t(rng.qubit(num_qubits)),
            6 => Gate::tdg(rng.qubit(num_qubits)),
            7 => {
                let q = random_distinct_qubits(rng, num_qubits, 2);
                Gate::cnot {
                    control: q[0],
                    target: q[1],
                }
            }
            8 => {
                let q = random_distinct_qubits(rng, num_qubits, 2);
                Gate::cz {
                    control: q[0],
                    target: q[1],
                }
            }
            9 => {
                let q = random_distinct_qubits(rng, num_qubits, 3);
                Gate::ccx {
                    control1: q[0],
                    control2: q[1],
                    target: q[2],
                }
            }
            10 => {
                let q = random_distinct_qubits(rng, num_qubits, 3);
                Gate::ccz {
                    control1: q[0],
                    control2: q[1],
                    target: q[2],
                }
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn randomized_multi_qubit_circuits_match_floating_oracle() {
        let mut rng = TestRng(0xc1c1_0700_1c5e_5eed);
        for case in 0..400 {
            let num_qubits = 1 + rng.next(5);
            let gate_count = 1 + rng.next(60);
            let mut circuit = Circuit::new(num_qubits);
            for _ in 0..gate_count {
                circuit.apply(random_gate(&mut rng, num_qubits));
            }
            assert_matches_floating_oracle(&circuit);
            assert_normalized(&exact_matrix(&circuit));
            assert!(
                circuit
                    .gates
                    .iter()
                    .all(|gate| !matches!(gate, Gate::rz(..))),
                "case {case} unexpectedly generated Rz"
            );
        }
    }

    #[test]
    fn inserting_inverse_pairs_preserves_exact_matrix_and_fingerprint() {
        let mut rng = TestRng(0x1d3e_1717_1e5e_5eed);
        for _ in 0..250 {
            let num_qubits = 2 + rng.next(3);
            let mut base = Circuit::new(num_qubits);
            for _ in 0..20 {
                base.apply(random_gate(&mut rng, num_qubits));
            }
            let mut decorated = base.clone();
            for _ in 0..10 {
                let qubit = rng.qubit(num_qubits);
                match rng.next(6) {
                    0 => {
                        decorated.apply(Gate::h(qubit));
                        decorated.apply(Gate::h(qubit));
                    }
                    1 => {
                        decorated.apply(Gate::x(qubit));
                        decorated.apply(Gate::x(qubit));
                    }
                    2 => {
                        decorated.apply(Gate::s(qubit));
                        decorated.apply(Gate::sdg(qubit));
                    }
                    3 => {
                        decorated.apply(Gate::t(qubit));
                        decorated.apply(Gate::tdg(qubit));
                    }
                    4 => {
                        let q = random_distinct_qubits(&mut rng, num_qubits, 2);
                        let gate = Gate::cnot {
                            control: q[0],
                            target: q[1],
                        };
                        decorated.apply(gate.clone());
                        decorated.apply(gate);
                    }
                    5 => {
                        let q = random_distinct_qubits(&mut rng, num_qubits, 2);
                        let gate = Gate::cz {
                            control: q[0],
                            target: q[1],
                        };
                        decorated.apply(gate.clone());
                        decorated.apply(gate);
                    }
                    _ => unreachable!(),
                }
            }
            let base = exact_matrix(&base);
            let decorated = exact_matrix(&decorated);
            assert_exactly_equal(&base, &decorated);
            assert_eq!(unitary_fingerprint(&base), unitary_fingerprint(&decorated));
        }
    }

    #[test]
    fn copy_from_preserves_every_exact_field() {
        let source = exact_matrix(&circuit(
            3,
            &[
                Gate::h(0),
                Gate::t(1),
                Gate::cnot {
                    control: 0,
                    target: 2,
                },
                Gate::sdg(2),
            ],
        ));
        let mut destination = UnitaryMatrix::identity(1).unwrap();
        destination.copy_from(&source);
        assert_exactly_equal(&source, &destination);
    }

    #[test]
    #[should_panic(expected = "Rz windows are rejected before matrix construction")]
    fn direct_rz_matrix_application_is_rejected() {
        let mut matrix = UnitaryMatrix::identity(1).unwrap();
        matrix
            .apply_gate_left(&Gate::rz(TAU / 7.0, 0), &[0])
            .unwrap();
    }
}
