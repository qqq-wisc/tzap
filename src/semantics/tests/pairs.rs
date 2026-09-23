//! Hand-written normal-circuit/PBC pairs, checked against both gate interpreters.
//! These fixtures do not use or implement automatic conversion.
use super::*;

fn cx(control: u32, target: u32) -> Gate {
    Gate::cnot { control, target }
}

// Arguments: test name, qubit count, normal gates, (Pauli word, pi/8 units),
// and optionally the Clifford suffix. Both interpreters must agree up to phase.
macro_rules! pair {
    ($name:ident, $n:expr, $gates:expr, $rotations:expr) => {
        pair!($name, $n, $gates, $rotations, std::iter::empty::<Gate>());
    };
    ($name:ident, $n:expr, $gates:expr, $rotations:expr, $suffix:expr) => {
        #[test]
        fn $name() {
            let mut pbc = rotations($n, $rotations);
            for gate in $suffix {
                pbc.push_output_clifford(gate).unwrap();
            }
            check(&pbc, $gates);
        }
    };
}

// Single-qubit identities, sign changes, and noncommuting rotation sequences.
pair!(sdg, 1, vec![Gate::sdg(0)], &[("Z", -2)]);
pair!(z, 1, vec![Gate::z(0)], &[("Z", 4)]);
pair!(
    x_then_z_is_y_up_to_phase,
    1,
    vec![Gate::x(0), Gate::z(0)],
    &[("Y", 4)]
);
pair!(
    x_t_x_negates_axis,
    1,
    vec![Gate::x(0), Gate::t(0), Gate::x(0)],
    &[("Z", -1)]
);
pair!(
    z_h_t_h_z_negates_x,
    1,
    vec![Gate::z(0), Gate::h(0), Gate::t(0), Gate::h(0), Gate::z(0)],
    &[("X", -1)]
);
pair!(
    s_h_t_h_sdg_is_negative_y,
    1,
    vec![Gate::s(0), Gate::h(0), Gate::t(0), Gate::h(0), Gate::sdg(0)],
    &[("Y", -1)]
);
pair!(
    two_t_gates_merge,
    1,
    vec![Gate::t(0), Gate::t(0)],
    &[("Z", 2)]
);
pair!(
    t_h_t_h_keeps_order,
    1,
    vec![Gate::t(0), Gate::h(0), Gate::t(0), Gate::h(0)],
    &[("Z", 1), ("X", 1)]
);
pair!(
    h_t_h_t_keeps_order,
    1,
    vec![Gate::h(0), Gate::t(0), Gate::h(0), Gate::t(0)],
    &[("X", 1), ("Z", 1)]
);
pair!(
    h_t_s_h_tdg_merges_in_x_basis,
    1,
    vec![Gate::h(0), Gate::t(0), Gate::s(0), Gate::h(0), Gate::tdg(0)],
    &[("X", 3), ("Z", -1)]
);

// Clifford conjugation and parity gadgets on two qubits.
pair!(
    cx_control_phase_stays_local,
    2,
    vec![cx(0, 1), Gate::t(0), cx(0, 1)],
    &[("ZI", 1)]
);
pair!(
    reversed_cx_target_phase,
    2,
    vec![cx(1, 0), Gate::tdg(0), cx(1, 0)],
    &[("ZZ", -1)]
);
pair!(
    cx_spreads_control_x,
    2,
    vec![cx(0, 1), Gate::h(0), Gate::t(0), Gate::h(0), cx(0, 1)],
    &[("XX", 1)]
);
pair!(
    cx_spreads_control_y,
    2,
    vec![
        cx(0, 1),
        Gate::sdg(0),
        Gate::h(0),
        Gate::t(0),
        Gate::h(0),
        Gate::s(0),
        cx(0, 1)
    ],
    &[("YX", 1)]
);
pair!(
    cz_spreads_x_to_xz,
    2,
    vec![
        Gate::cz {
            control: 0,
            target: 1
        },
        Gate::h(0),
        Gate::t(0),
        Gate::h(0),
        Gate::cz {
            control: 0,
            target: 1
        }
    ],
    &[("XZ", 1)]
);
pair!(
    cz_spreads_other_x_to_zx,
    2,
    vec![
        Gate::cz {
            control: 0,
            target: 1
        },
        Gate::h(1),
        Gate::tdg(1),
        Gate::h(1),
        Gate::cz {
            control: 0,
            target: 1
        }
    ],
    &[("ZX", -1)]
);
pair!(
    cx_conjugation_of_yy_has_negative_sign,
    2,
    vec![
        cx(0, 1),
        Gate::sdg(0),
        Gate::sdg(1),
        Gate::h(0),
        Gate::h(1),
        cx(0, 1),
        Gate::t(1),
        cx(0, 1),
        Gate::h(1),
        Gate::h(0),
        Gate::s(1),
        Gate::s(0),
        cx(0, 1)
    ],
    &[("XZ", -1)]
);
pair!(
    parity_phase_then_local_x_rotation,
    2,
    vec![
        cx(0, 1),
        Gate::t(1),
        cx(0, 1),
        Gate::h(0),
        Gate::tdg(0),
        Gate::h(0)
    ],
    &[("ZZ", 1), ("XI", -1)]
);

// Larger supports and idle wires test tensor ordering as well as angles.
pair!(
    three_qubit_ladder_parity,
    3,
    vec![cx(0, 1), cx(1, 2), Gate::t(2), cx(1, 2), cx(0, 1)],
    &[("ZZZ", 1)]
);
pair!(
    three_qubit_fanin_inverse_parity,
    3,
    vec![cx(0, 2), cx(1, 2), Gate::tdg(2), cx(1, 2), cx(0, 2)],
    &[("ZZZ", -1)]
);
pair!(
    three_qubit_x_parity,
    3,
    vec![
        Gate::h(0),
        Gate::h(1),
        Gate::h(2),
        cx(0, 2),
        cx(1, 2),
        Gate::t(2),
        cx(1, 2),
        cx(0, 2),
        Gate::h(2),
        Gate::h(1),
        Gate::h(0)
    ],
    &[("XXX", 1)]
);
pair!(
    three_qubit_mixed_xyz_parity,
    3,
    vec![
        Gate::h(0),
        Gate::sdg(1),
        Gate::h(1),
        cx(0, 2),
        cx(1, 2),
        Gate::t(2),
        cx(1, 2),
        cx(0, 2),
        Gate::h(1),
        Gate::s(1),
        Gate::h(0)
    ],
    &[("XYZ", 1)]
);
pair!(
    nonadjacent_parity_with_two_idle_wires,
    4,
    vec![cx(0, 3), Gate::s(3), cx(0, 3)],
    &[("ZIIZ", 2)]
);
pair!(
    ccx_with_target_on_first_wire,
    3,
    vec![Gate::ccx {
        control1: 1,
        control2: 2,
        target: 0
    }],
    &[
        ("XII", 1),
        ("IZI", 1),
        ("IIZ", 1),
        ("XZI", -1),
        ("XIZ", -1),
        ("IZZ", -1),
        ("XZZ", 1)
    ]
);

// Nonempty output frames: moving gates to the suffix changes rotation axes.
pair!(
    x_then_t_retains_x_frame,
    1,
    vec![Gate::x(0), Gate::t(0)],
    &[("Z", -1)],
    [Gate::x(0)]
);
pair!(
    h_t_s_retains_h_frame,
    1,
    vec![Gate::h(0), Gate::t(0), Gate::s(0)],
    &[("X", 3)],
    [Gate::h(0)]
);
pair!(
    h_t_h_t_h_retains_h_frame,
    1,
    vec![Gate::h(0), Gate::t(0), Gate::h(0), Gate::t(0), Gate::h(0)],
    &[("X", 1), ("Z", 1)],
    [Gate::h(0)]
);
pair!(
    bell_preparation_then_t_retains_entangling_frame,
    2,
    vec![Gate::h(0), cx(0, 1), Gate::t(1)],
    &[("XZ", 1)],
    [Gate::h(0), cx(0, 1)]
);
pair!(
    ghz_preparation_then_t_retains_entangling_frame,
    3,
    vec![Gate::h(0), cx(0, 1), cx(1, 2), Gate::t(2)],
    &[("XZZ", 1)],
    [Gate::h(0), cx(0, 1), cx(1, 2)]
);
pair!(
    cx_with_two_different_phase_locations,
    2,
    vec![cx(0, 1), Gate::t(0), Gate::tdg(1)],
    &[("ZI", 1), ("ZZ", -1)],
    [cx(0, 1)]
);
