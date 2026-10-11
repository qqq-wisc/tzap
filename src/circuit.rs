//! Circuit representation: gates, qubits, and display.

use std::fmt;

/// Index of a qubit within a [`Circuit`].
pub type Qubit = u32;
/// Index of a classical bit within a [`Circuit`].
pub type CBit = u32;

/// A supported gate kind, ordered canonically for user-facing gate-set output.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum GateKind {
    H,
    X,
    Z,
    S,
    Sdg,
    T,
    Tdg,
    Rz,
    Cx,
    Cz,
    Ccx,
    Ccz,
    Measure,
    Reset,
    P,
    Y,
    Sx,
    Rx,
    Ry,
    Swap,
    Cy,
    Cp,
    Crx,
    Cry,
    Crz,
    Ch,
    Cswap,
}

impl GateKind {
    pub const ALL: [GateKind; 27] = [
        Self::H,
        Self::X,
        Self::Z,
        Self::S,
        Self::Sdg,
        Self::T,
        Self::Tdg,
        Self::Rz,
        Self::Cx,
        Self::Cz,
        Self::Ccx,
        Self::Ccz,
        Self::Measure,
        Self::Reset,
        Self::P,
        Self::Y,
        Self::Sx,
        Self::Rx,
        Self::Ry,
        Self::Swap,
        Self::Cy,
        Self::Cp,
        Self::Crx,
        Self::Cry,
        Self::Crz,
        Self::Ch,
        Self::Cswap,
    ];

    pub const fn name(self) -> &'static str {
        match self {
            Self::H => "h",
            Self::X => "x",
            Self::Z => "z",
            Self::S => "s",
            Self::Sdg => "sdg",
            Self::T => "t",
            Self::Tdg => "tdg",
            Self::Rz => "rz",
            Self::Cx => "cx",
            Self::Cz => "cz",
            Self::Ccx => "ccx",
            Self::Ccz => "ccz",
            Self::Measure => "measure",
            Self::Reset => "reset",
            Self::P => "p",
            Self::Y => "y",
            Self::Sx => "sx",
            Self::Rx => "rx",
            Self::Ry => "ry",
            Self::Swap => "swap",
            Self::Cy => "cy",
            Self::Cp => "cp",
            Self::Crx => "crx",
            Self::Cry => "cry",
            Self::Crz => "crz",
            Self::Ch => "ch",
            Self::Cswap => "cswap",
        }
    }

    pub fn parse(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|kind| kind.name() == name)
    }
}

/// Compact set of [`GateKind`] values with stable canonical iteration.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct GateSet(u32);

impl GateSet {
    pub const EMPTY: Self = Self(0);

    /// Construct a set from already-validated bits in a constant context.
    pub(crate) const fn from_bits_const(bits: u32) -> Self {
        Self(bits)
    }

    pub const fn singleton(kind: GateKind) -> Self {
        Self(1 << kind as u8)
    }

    pub fn from_kinds(kinds: impl IntoIterator<Item = GateKind>) -> Self {
        let mut set = Self::EMPTY;
        for kind in kinds {
            set.insert(kind);
        }
        set
    }

    pub fn insert(&mut self, kind: GateKind) {
        self.0 |= 1 << kind as u8;
    }

    pub const fn contains(self, kind: GateKind) -> bool {
        self.0 & (1 << kind as u8) != 0
    }

    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    pub const fn intersection(self, other: Self) -> Self {
        Self(self.0 & other.0)
    }

    pub const fn difference(self, other: Self) -> Self {
        Self(self.0 & !other.0)
    }

    pub const fn is_subset(self, other: Self) -> bool {
        self.0 & !other.0 == 0
    }

    pub const fn bits(self) -> u32 {
        self.0
    }

    pub fn from_bits(bits: u32) -> Option<Self> {
        (bits & !((1 << GateKind::ALL.len()) - 1) == 0).then_some(Self(bits))
    }

    pub fn iter(self) -> impl Iterator<Item = GateKind> {
        GateKind::ALL
            .into_iter()
            .filter(move |kind| self.contains(*kind))
    }

    pub fn names(self) -> impl Iterator<Item = &'static str> {
        self.iter().map(GateKind::name)
    }
}

impl fmt::Display for GateSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{{}}}", self.names().collect::<Vec<_>>().join(", "))
    }
}

/// A single quantum (or classical `measure`/`reset`) operation.
///
/// Variants are lowercase to mirror their QASM gate names.
#[allow(non_camel_case_types)]
#[derive(Clone, Debug)]
pub enum Gate {
    x(Qubit),
    h(Qubit),
    s(Qubit),
    sdg(Qubit),
    z(Qubit),
    t(Qubit),
    tdg(Qubit),
    rz(crate::angle::Angle, Qubit),
    cnot {
        control: Qubit,
        target: Qubit,
    },
    cz {
        control: Qubit,
        target: Qubit,
    },
    ccx {
        control1: Qubit,
        control2: Qubit,
        target: Qubit,
    },
    ccz {
        control1: Qubit,
        control2: Qubit,
        target: Qubit,
    },
    measure {
        qubit: Qubit,
        cbit: CBit,
    },
    reset(Qubit),
    /// Phase gate: diag(1, exp(i * lambda)). Angles are in radians.
    p(crate::angle::Angle, Qubit),
    /// Pauli Y: [[0, -i], [i, 0]].
    y(Qubit),
    /// Principal square root of X: H S H.
    sx(Qubit),
    /// X rotation: exp(-i * theta * X / 2), with theta in radians.
    rx(crate::angle::Angle, Qubit),
    /// Y rotation: exp(-i * theta * Y / 2), with theta in radians.
    ry(crate::angle::Angle, Qubit),
    /// Exchange the states of two distinct qubits.
    swap(Qubit, Qubit),
    /// Apply Y to target when control is 1.
    cy {
        control: Qubit,
        target: Qubit,
    },
    /// Controlled phase: diag(1, 1, 1, exp(i * lambda)), in radians.
    cp {
        lambda: crate::angle::Angle,
        control: Qubit,
        target: Qubit,
    },
    /// Controlled exp(-i * theta * X / 2), with theta in radians.
    crx {
        theta: crate::angle::Angle,
        control: Qubit,
        target: Qubit,
    },
    /// Controlled exp(-i * theta * Y / 2), with theta in radians.
    cry {
        theta: crate::angle::Angle,
        control: Qubit,
        target: Qubit,
    },
    /// Controlled exp(-i * theta * Z / 2), with theta in radians.
    crz {
        theta: crate::angle::Angle,
        control: Qubit,
        target: Qubit,
    },
    /// Apply H to target when control is 1.
    ch {
        control: Qubit,
        target: Qubit,
    },
    /// Exchange first and second when control is 1. All operands are distinct.
    cswap {
        control: Qubit,
        first: Qubit,
        second: Qubit,
    },
}

// Controlled parameter equality must retain the relative phase of a 2π turn.
impl PartialEq for Gate {
    fn eq(&self, rhs: &Self) -> bool {
        match (self, rhs) {
            (Self::x(a), Self::x(b)) => a == b,
            (Self::h(a), Self::h(b)) => a == b,
            (Self::s(a), Self::s(b)) => a == b,
            (Self::sdg(a), Self::sdg(b)) => a == b,
            (Self::z(a), Self::z(b)) => a == b,
            (Self::t(a), Self::t(b)) => a == b,
            (Self::tdg(a), Self::tdg(b)) => a == b,
            (Self::y(a), Self::y(b)) => a == b,
            (Self::sx(a), Self::sx(b)) => a == b,
            (Self::reset(a), Self::reset(b)) => a == b,
            (Self::rz(a, q), Self::rz(b, r)) => q == r && a == b,
            (Self::p(a, q), Self::p(b, r)) => q == r && a == b,
            (Self::rx(a, q), Self::rx(b, r)) => q == r && a == b,
            (Self::ry(a, q), Self::ry(b, r)) => q == r && a == b,
            (
                Self::cnot {
                    control: a,
                    target: b,
                },
                Self::cnot {
                    control: c,
                    target: d,
                },
            ) => a == c && b == d,
            (
                Self::cz {
                    control: a,
                    target: b,
                },
                Self::cz {
                    control: c,
                    target: d,
                },
            ) => a == c && b == d,
            (
                Self::cy {
                    control: a,
                    target: b,
                },
                Self::cy {
                    control: c,
                    target: d,
                },
            ) => a == c && b == d,
            (
                Self::ch {
                    control: a,
                    target: b,
                },
                Self::ch {
                    control: c,
                    target: d,
                },
            ) => a == c && b == d,
            (
                Self::ccx {
                    control1: a,
                    control2: b,
                    target: c,
                },
                Self::ccx {
                    control1: d,
                    control2: e,
                    target: f,
                },
            ) => a == d && b == e && c == f,
            (
                Self::ccz {
                    control1: a,
                    control2: b,
                    target: c,
                },
                Self::ccz {
                    control1: d,
                    control2: e,
                    target: f,
                },
            ) => a == d && b == e && c == f,
            (
                Self::cp {
                    lambda: a,
                    control: b,
                    target: c,
                },
                Self::cp {
                    lambda: d,
                    control: e,
                    target: f,
                },
            ) => a.literal_eq(d) && b == e && c == f,
            (
                Self::crx {
                    theta: a,
                    control: b,
                    target: c,
                },
                Self::crx {
                    theta: d,
                    control: e,
                    target: f,
                },
            ) => a.literal_eq(d) && b == e && c == f,
            (
                Self::cry {
                    theta: a,
                    control: b,
                    target: c,
                },
                Self::cry {
                    theta: d,
                    control: e,
                    target: f,
                },
            ) => a.literal_eq(d) && b == e && c == f,
            (
                Self::crz {
                    theta: a,
                    control: b,
                    target: c,
                },
                Self::crz {
                    theta: d,
                    control: e,
                    target: f,
                },
            ) => a.literal_eq(d) && b == e && c == f,
            (Self::swap(a, b), Self::swap(c, d)) => a == c && b == d,
            (
                Self::cswap {
                    control: a,
                    first: b,
                    second: c,
                },
                Self::cswap {
                    control: d,
                    first: e,
                    second: f,
                },
            ) => a == d && b == e && c == f,
            (Self::measure { qubit: a, cbit: b }, Self::measure { qubit: c, cbit: d }) => {
                a == c && b == d
            }
            _ => false,
        }
    }
}
/// An ordered sequence of [`Gate`]s over a fixed number of qubits.
///
/// Gate metadata is derived from `gates` rather than cached beside it. This
/// keeps the public `gates` vector safe to mutate without silently leaving
/// optimizer policy decisions out of sync with the actual circuit.
#[derive(Clone, Debug)]
pub struct Circuit {
    pub num_qubits: usize,
    pub num_cbits: usize,
    pub gates: Vec<Gate>,
}

impl Circuit {
    /// An empty circuit over `num_qubits` qubits and no classical bits.
    pub fn new(num_qubits: usize) -> Self {
        Circuit {
            num_qubits,
            num_cbits: 0,
            gates: Vec::new(),
        }
    }

    /// An empty circuit over `num_qubits` qubits and `num_cbits` classical
    /// bits. Use this instead of [`Circuit::new`] when the circuit contains
    /// `measure` gates.
    pub fn with_cbits(num_qubits: usize, num_cbits: usize) -> Self {
        Circuit {
            num_qubits,
            num_cbits,
            gates: Vec::new(),
        }
    }

    pub(crate) fn empty_like(&self) -> Self {
        Self::with_cbits(self.num_qubits, self.num_cbits)
    }
    #[cfg(test)]
    pub(crate) fn replacing_gates(&self, gates: Vec<Gate>) -> Self {
        Self {
            gates,
            ..self.empty_like()
        }
    }
    /// Append `gate`.
    pub fn apply(&mut self, gate: Gate) {
        self.gates.push(gate);
    }

    /// Gate kinds currently present, derived from the circuit's actual gates.
    pub fn gate_set(&self) -> GateSet {
        GateSet::from_kinds(self.gates.iter().map(Gate::kind))
    }

    pub fn has_toffoli(&self) -> bool {
        self.gates
            .iter()
            .any(|gate| matches!(gate, Gate::ccx { .. }))
    }

    pub fn has_ccz(&self) -> bool {
        self.gates
            .iter()
            .any(|gate| matches!(gate, Gate::ccz { .. }))
    }

    pub fn has_measurement(&self) -> bool {
        self.gates
            .iter()
            .any(|gate| matches!(gate, Gate::measure { .. } | Gate::reset(_)))
    }

    /// Serialize to OpenQASM 2.0. See [`crate::qasm`] for the supported subset.
    pub fn to_qasm(&self) -> String {
        crate::qasm::serialize(self)
    }

    /// Parse from OpenQASM 2.0. See [`crate::qasm`] for the supported subset.
    pub fn from_qasm(qasm: &str) -> Result<Self, String> {
        crate::qasm::parse(qasm)
    }
}

impl Gate {
    /// Borrow the checked parameter of a phase or rotation gate.
    pub fn angle(&self) -> Option<&crate::angle::Angle> {
        match self {
            Self::rz(a, _) | Self::p(a, _) | Self::rx(a, _) | Self::ry(a, _) => Some(a),
            Self::cp { lambda: a, .. }
            | Self::crx { theta: a, .. }
            | Self::cry { theta: a, .. }
            | Self::crz { theta: a, .. } => Some(a),
            _ => None,
        }
    }

    pub fn p_f64(theta: f64, q: Qubit) -> Result<Self, crate::angle::AngleError> {
        Ok(Self::p(crate::angle::Angle::from_f64(theta)?, q))
    }
    pub fn rx_f64(theta: f64, q: Qubit) -> Result<Self, crate::angle::AngleError> {
        Ok(Self::rx(crate::angle::Angle::from_f64(theta)?, q))
    }
    pub fn ry_f64(theta: f64, q: Qubit) -> Result<Self, crate::angle::AngleError> {
        Ok(Self::ry(crate::angle::Angle::from_f64(theta)?, q))
    }
    pub fn rz_f64(value: f64, q: Qubit) -> Result<Self, crate::angle::AngleError> {
        Ok(Self::rz(crate::angle::Angle::from_f64(value)?, q))
    }

    /// This gate's kind, independent of its operands.
    pub fn kind(&self) -> GateKind {
        GateKind::of(self)
    }

    /// The same gate with every qubit operand sent through `f`. Classical
    /// bits are untouched.
    pub fn map_qubits(&self, mut f: impl FnMut(Qubit) -> Qubit) -> Gate {
        match self {
            Gate::x(q) => Gate::x(f(*q)),
            Gate::h(q) => Gate::h(f(*q)),
            Gate::s(q) => Gate::s(f(*q)),
            Gate::sdg(q) => Gate::sdg(f(*q)),
            Gate::z(q) => Gate::z(f(*q)),
            Gate::t(q) => Gate::t(f(*q)),
            Gate::tdg(q) => Gate::tdg(f(*q)),
            Gate::rz(theta, q) => Gate::rz(theta.clone(), f(*q)),
            Gate::p(theta, q) => Gate::p(theta.clone(), f(*q)),
            Gate::y(q) => Gate::y(f(*q)),
            Gate::sx(q) => Gate::sx(f(*q)),
            Gate::rx(theta, q) => Gate::rx(theta.clone(), f(*q)),
            Gate::ry(theta, q) => Gate::ry(theta.clone(), f(*q)),
            Gate::swap(a, b) => Gate::swap(f(*a), f(*b)),
            Gate::cy { control, target } => Gate::cy {
                control: f(*control),
                target: f(*target),
            },
            Gate::cp {
                lambda,
                control,
                target,
            } => Gate::cp {
                lambda: lambda.clone(),
                control: f(*control),
                target: f(*target),
            },
            Gate::crx {
                theta,
                control,
                target,
            } => Gate::crx {
                theta: theta.clone(),
                control: f(*control),
                target: f(*target),
            },
            Gate::cry {
                theta,
                control,
                target,
            } => Gate::cry {
                theta: theta.clone(),
                control: f(*control),
                target: f(*target),
            },
            Gate::crz {
                theta,
                control,
                target,
            } => Gate::crz {
                theta: theta.clone(),
                control: f(*control),
                target: f(*target),
            },
            Gate::ch { control, target } => Gate::ch {
                control: f(*control),
                target: f(*target),
            },
            Gate::cswap {
                control,
                first,
                second,
            } => Gate::cswap {
                control: f(*control),
                first: f(*first),
                second: f(*second),
            },
            Gate::cnot { control, target } => Gate::cnot {
                control: f(*control),
                target: f(*target),
            },
            Gate::cz { control, target } => Gate::cz {
                control: f(*control),
                target: f(*target),
            },
            Gate::ccx {
                control1,
                control2,
                target,
            } => Gate::ccx {
                control1: f(*control1),
                control2: f(*control2),
                target: f(*target),
            },
            Gate::ccz {
                control1,
                control2,
                target,
            } => Gate::ccz {
                control1: f(*control1),
                control2: f(*control2),
                target: f(*target),
            },
            Gate::measure { qubit, cbit } => Gate::measure {
                qubit: f(*qubit),
                cbit: *cbit,
            },
            Gate::reset(q) => Gate::reset(f(*q)),
        }
    }
}

/// Canonical operand order for a symmetric two-qubit gate.
pub(crate) fn canonical_pair<T: Copy + Ord>(a: T, b: T) -> (T, T) {
    if a <= b { (a, b) } else { (b, a) }
}

/// Canonical operand order for a Toffoli: sorted controls and fixed target.
pub(crate) fn canonical_ccx<T: Copy + Ord>(control1: T, control2: T, target: T) -> (T, T, T) {
    let (control1, control2) = canonical_pair(control1, control2);
    (control1, control2, target)
}

/// Canonical operand order for a fully symmetric three-qubit gate.
pub(crate) fn canonical_triple<T: Copy + Ord>(a: T, b: T, c: T) -> (T, T, T) {
    let mut operands = [a, b, c];
    operands.sort_unstable();
    (operands[0], operands[1], operands[2])
}

impl GateKind {
    /// Whether the exact Clifford+T matrix engine can interpret this kind.
    /// Parametric rotations remain boundaries even at special numeric angles.
    pub(crate) const fn has_exact_matrix(self) -> bool {
        !matches!(
            self,
            Self::Rz
                | Self::P
                | Self::Rx
                | Self::Ry
                | Self::Cp
                | Self::Crx
                | Self::Cry
                | Self::Crz
                | Self::Measure
                | Self::Reset
        )
    }

    pub const fn of(gate: &Gate) -> Self {
        match gate {
            Gate::h(_) => Self::H,
            Gate::x(_) => Self::X,
            Gate::z(_) => Self::Z,
            Gate::s(_) => Self::S,
            Gate::sdg(_) => Self::Sdg,
            Gate::t(_) => Self::T,
            Gate::tdg(_) => Self::Tdg,
            Gate::rz(..) => Self::Rz,
            Gate::cnot { .. } => Self::Cx,
            Gate::cz { .. } => Self::Cz,
            Gate::ccx { .. } => Self::Ccx,
            Gate::ccz { .. } => Self::Ccz,
            Gate::measure { .. } => Self::Measure,
            Gate::reset(_) => Self::Reset,
            Gate::p(..) => Self::P,
            Gate::y(_) => Self::Y,
            Gate::sx(_) => Self::Sx,
            Gate::rx(..) => Self::Rx,
            Gate::ry(..) => Self::Ry,
            Gate::swap(..) => Self::Swap,
            Gate::cy { .. } => Self::Cy,
            Gate::cp { .. } => Self::Cp,
            Gate::crx { .. } => Self::Crx,
            Gate::cry { .. } => Self::Cry,
            Gate::crz { .. } => Self::Crz,
            Gate::ch { .. } => Self::Ch,
            Gate::cswap { .. } => Self::Cswap,
        }
    }
}

impl fmt::Display for Gate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Gate::x(q) => write!(f, "x q{q}"),
            Gate::h(q) => write!(f, "h q{q}"),
            Gate::s(q) => write!(f, "s q{q}"),
            Gate::sdg(q) => write!(f, "sdg q{q}"),
            Gate::z(q) => write!(f, "z q{q}"),
            Gate::t(q) => write!(f, "t q{q}"),
            Gate::tdg(q) => write!(f, "tdg q{q}"),
            Gate::rz(theta, q) => write!(f, "rz({theta:.4}) q{q}"),
            Gate::p(theta, q) => write!(f, "p({theta:.4}) q{q}"),
            Gate::y(q) => write!(f, "y q{q}"),
            Gate::sx(q) => write!(f, "sx q{q}"),
            Gate::rx(theta, q) => write!(f, "rx({theta:.4}) q{q}"),
            Gate::ry(theta, q) => write!(f, "ry({theta:.4}) q{q}"),
            Gate::swap(a, b) => write!(f, "swap q{a}, q{b}"),
            Gate::cy { control, target } => write!(f, "cy q{control}, q{target}"),
            Gate::cp {
                lambda,
                control,
                target,
            } => write!(f, "cp({lambda:.4}) q{control}, q{target}"),
            Gate::crx {
                theta,
                control,
                target,
            } => write!(f, "crx({theta:.4}) q{control}, q{target}"),
            Gate::cry {
                theta,
                control,
                target,
            } => write!(f, "cry({theta:.4}) q{control}, q{target}"),
            Gate::crz {
                theta,
                control,
                target,
            } => write!(f, "crz({theta:.4}) q{control}, q{target}"),
            Gate::ch { control, target } => write!(f, "ch q{control}, q{target}"),
            Gate::cswap {
                control,
                first,
                second,
            } => write!(f, "cswap q{control}, q{first}, q{second}"),
            Gate::cnot { control, target } => write!(f, "cnot q{control}, q{target}"),
            Gate::cz { control, target } => write!(f, "cz q{control}, q{target}"),
            Gate::ccx {
                control1,
                control2,
                target,
            } => {
                write!(f, "ccx q{control1}, q{control2}, q{target}")
            }
            Gate::ccz {
                control1,
                control2,
                target,
            } => {
                write!(f, "ccz q{control1}, q{control2}, q{target}")
            }
            Gate::measure { qubit, cbit } => write!(f, "measure q{qubit} -> c{cbit}"),
            Gate::reset(q) => write!(f, "reset q{q}"),
        }
    }
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "Circuit ({} qubits, {} gates):",
            self.num_qubits,
            self.gates.len()
        )?;
        for (i, gate) in self.gates.iter().enumerate() {
            writeln!(f, "  {i}: {gate}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod gate_set_tests {
    use super::*;
    use crate::super_opt::{BASE_GATE_SET, SUPPORTED_GATE_SET};

    #[test]
    fn gate_set_uses_the_canonical_qasm_order() {
        let mut circuit = Circuit::with_cbits(3, 1);
        circuit.apply(Gate::reset(2));
        circuit.apply(Gate::ccz {
            control1: 0,
            control2: 1,
            target: 2,
        });
        circuit.apply(Gate::h(0));
        circuit.apply(Gate::rz_f64(0.3, 1).unwrap());
        circuit.apply(Gate::measure { qubit: 0, cbit: 0 });
        assert_eq!(
            circuit.gate_set().to_string(),
            "{h, rz, ccz, measure, reset}"
        );
    }

    #[test]
    fn superopt_gate_sets_are_exact() {
        assert_eq!(BASE_GATE_SET.to_string(), "{h, x, z, s, sdg, t, tdg, cx}");
        assert_eq!(
            SUPPORTED_GATE_SET.to_string(),
            "{h, x, z, s, sdg, t, tdg, cx, cz, ccx, ccz}"
        );
    }

    #[test]
    fn metadata_tracks_direct_public_gate_mutation() {
        let mut circuit = Circuit::new(3);
        circuit.gates.push(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        assert_eq!(circuit.gate_set(), GateSet::singleton(GateKind::Ccx));
        assert!(circuit.has_toffoli());

        circuit.gates.clear();
        circuit.gates.push(Gate::reset(1));
        assert_eq!(circuit.gate_set(), GateSet::singleton(GateKind::Reset));
        assert!(!circuit.has_toffoli());
        assert!(circuit.has_measurement());
    }
}

/// The qubits a gate acts on, as `(count, operands)` — the first `count`
/// entries of the array are the operands, in the same order [`qubits_of`]
/// reports them.
///
/// This is the form to reach for on a hot path: no gate acts on more than
/// three qubits, so the operands ride in registers and a per-gate scan of a
/// million-gate circuit allocates nothing. [`qubits_of`]'s `Vec` costs a
/// malloc/free per gate, which dominated `pass::depth` before this existed.
pub fn qubit_operands(gate: &Gate) -> (usize, [Qubit; 3]) {
    match gate {
        Gate::x(q)
        | Gate::h(q)
        | Gate::s(q)
        | Gate::sdg(q)
        | Gate::z(q)
        | Gate::t(q)
        | Gate::tdg(q)
        | Gate::rz(_, q)
        | Gate::p(_, q)
        | Gate::y(q)
        | Gate::sx(q)
        | Gate::rx(_, q)
        | Gate::ry(_, q)
        | Gate::reset(q) => (1, [*q, 0, 0]),
        Gate::swap(control, target)
        | Gate::cnot { control, target }
        | Gate::cz { control, target }
        | Gate::cy { control, target }
        | Gate::ch { control, target }
        | Gate::cp {
            control, target, ..
        }
        | Gate::crx {
            control, target, ..
        }
        | Gate::cry {
            control, target, ..
        }
        | Gate::crz {
            control, target, ..
        } => (2, [*control, *target, 0]),
        Gate::ccx {
            control1,
            control2,
            target,
        }
        | Gate::ccz {
            control1,
            control2,
            target,
        } => (3, [*control1, *control2, *target]),
        Gate::cswap {
            control,
            first,
            second,
        } => (3, [*control, *first, *second]),
        Gate::measure { qubit, .. } => (1, [*qubit, 0, 0]),
    }
}

/// Return the qubits a gate acts on. See [`qubit_operands`] for an
/// allocation-free equivalent.
pub fn qubits_of(gate: &Gate) -> Vec<Qubit> {
    let (count, operands) = qubit_operands(gate);
    operands[..count].to_vec()
}

/// Remap a gate's qubits through a lookup table: qubit i becomes its index in `qubits`.
/// Classical bits are not remapped.
pub fn remap_gate(gate: &Gate, qubits: &[Qubit]) -> Gate {
    gate.map_qubits(|q| {
        qubits
            .iter()
            .position(|&x| x == q)
            .expect("gate operand is in the remap table") as Qubit
    })
}

/// Build a compact circuit with qubits remapped to 0..n.
pub fn remap_subcircuit(gates: &[Gate], qubits: &[Qubit]) -> Circuit {
    let n = qubits.len();
    let mut c = Circuit::new(n);
    for g in gates {
        c.apply(remap_gate(g, qubits));
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bell_pair() {
        let mut c = Circuit::new(2);
        c.apply(Gate::h(0));
        c.apply(Gate::cnot {
            control: 0,
            target: 1,
        });
        assert_eq!(c.gates.len(), 2);
        let s = format!("{c}");
        assert!(s.contains("h q0"));
        assert!(s.contains("cnot q0, q1"));
        println!("{c}");
    }

    #[test]
    fn ghz_state() {
        let n = 4;
        let mut c = Circuit::new(n);
        c.apply(Gate::h(0));
        for i in 0..n as Qubit - 1 {
            c.apply(Gate::cnot {
                control: i,
                target: i + 1,
            });
        }
        assert_eq!(c.gates.len(), 4);
        println!("{c}");
    }

    #[test]
    fn t_gate_decomposition_of_rz() {
        let mut c = Circuit::new(1);
        c.apply(Gate::t(0));
        c.apply(Gate::s(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        let s = format!("{c}");
        assert!(s.contains("t q0"));
        assert!(s.contains("s q0"));
        assert!(s.contains("rz((1*pi/4)) q0"));
        println!("{c}");
    }

    #[test]
    fn ccx_gate() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(2));
        c.apply(Gate::ccx {
            control1: 0,
            control2: 1,
            target: 2,
        });
        c.apply(Gate::h(2));
        assert_eq!(c.gates.len(), 3);
        let s = format!("{c}");
        assert!(s.contains("ccx q0, q1, q2"));
        println!("{c}");
    }

    #[test]
    fn ccz_gate_display_and_metadata() {
        let mut c = Circuit::new(3);
        c.apply(Gate::ccz {
            control1: 2,
            control2: 0,
            target: 1,
        });

        assert!(format!("{c}").contains("ccz q2, q0, q1"));
        assert!(!c.has_toffoli());
        assert!(c.has_ccz());
        assert_eq!(qubits_of(&c.gates[0]), vec![2, 0, 1]);
    }

    #[test]
    fn ccz_remap() {
        let gate = Gate::ccz {
            control1: 8,
            control2: 2,
            target: 5,
        };

        assert!(matches!(
            remap_gate(&gate, &[2, 5, 8]),
            Gate::ccz {
                control1: 2,
                control2: 0,
                target: 1
            }
        ));
    }

    #[test]
    fn qft_3qubit() {
        let mut c = Circuit::new(3);
        c.apply(Gate::h(0));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 2.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::cnot {
            control: 1,
            target: 0,
        });
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 4.0", 1).unwrap(),
            0,
        ));
        c.apply(Gate::cnot {
            control: 2,
            target: 0,
        });
        c.apply(Gate::h(1));
        c.apply(Gate::rz(
            crate::angle_expr::parse("pi / 2.0", 1).unwrap(),
            1,
        ));
        c.apply(Gate::cnot {
            control: 2,
            target: 1,
        });
        c.apply(Gate::h(2));
        assert_eq!(c.num_qubits, 3);
        assert_eq!(c.gates.len(), 9);
        println!("{c}");
    }

    #[test]
    fn z_gate_display() {
        let mut c = Circuit::new(1);
        c.apply(Gate::z(0));
        let s = format!("{c}");
        assert!(s.contains("z q0"));
    }

    #[test]
    fn sdg_gate_display() {
        let mut c = Circuit::new(1);
        c.apply(Gate::sdg(0));
        let s = format!("{c}");
        assert!(s.contains("sdg q0"));
    }

    #[test]
    fn cz_gate_display_and_metadata() {
        let mut c = Circuit::new(3);
        c.apply(Gate::cz {
            control: 2,
            target: 0,
        });
        assert!(format!("{c}").contains("cz q2, q0"));
        assert!(!c.has_toffoli());
        assert!(!c.has_measurement());
    }

    #[test]
    fn cz_qubits_of_preserves_operand_order() {
        let g = Gate::cz {
            control: 4,
            target: 1,
        };
        assert_eq!(qubits_of(&g), vec![4, 1]);
    }

    #[test]
    fn cz_remap() {
        let g = Gate::cz {
            control: 7,
            target: 3,
        };
        let remapped = remap_gate(&g, &[3, 7]);
        assert!(matches!(
            remapped,
            Gate::cz {
                control: 1,
                target: 0
            }
        ));
    }

    #[test]
    fn cz_remap_subcircuit() {
        let gates = vec![
            Gate::t(8),
            Gate::cz {
                control: 8,
                target: 2,
            },
        ];
        let remapped = remap_subcircuit(&gates, &[2, 8]);
        assert_eq!(remapped.num_qubits, 2);
        assert!(matches!(remapped.gates[0], Gate::t(1)));
        assert!(matches!(
            remapped.gates[1],
            Gate::cz {
                control: 1,
                target: 0
            }
        ));
    }

    #[test]
    fn measure_gate_display() {
        let mut c = Circuit::with_cbits(1, 1);
        c.apply(Gate::measure { qubit: 0, cbit: 0 });
        let s = format!("{c}");
        assert!(s.contains("measure q0 -> c0"));
        assert!(c.has_measurement());
    }

    #[test]
    fn reset_gate_display() {
        let mut c = Circuit::new(1);
        c.apply(Gate::reset(0));
        let s = format!("{c}");
        assert!(s.contains("reset q0"));
        assert!(c.has_measurement());
    }

    #[test]
    fn measure_qubits_of() {
        let g = Gate::measure { qubit: 3, cbit: 7 };
        assert_eq!(qubits_of(&g), vec![3]);
    }

    #[test]
    fn reset_qubits_of() {
        let g = Gate::reset(2);
        assert_eq!(qubits_of(&g), vec![2]);
    }

    #[test]
    fn measure_remap() {
        let g = Gate::measure { qubit: 5, cbit: 2 };
        let remapped = remap_gate(&g, &[5]);
        match remapped {
            Gate::measure { qubit: 0, cbit: 2 } => {}
            _ => panic!("expected measure q0 -> c2, got {remapped:?}"),
        }
    }

    #[test]
    fn reset_remap() {
        let g = Gate::reset(5);
        let remapped = remap_gate(&g, &[5]);
        assert!(matches!(remapped, Gate::reset(0)));
    }

    #[test]
    fn with_cbits_default_fields() {
        let c = Circuit::with_cbits(2, 3);
        assert_eq!(c.num_qubits, 2);
        assert_eq!(c.num_cbits, 3);
        assert!(!c.has_measurement());
        assert!(!c.has_toffoli());
        assert!(!c.has_ccz());
        assert_eq!(c.gates.len(), 0);
    }

    #[test]
    fn has_measurement_flag_set_by_reset_alone() {
        // reset has no cbits but still counts as measurement
        let mut c = Circuit::new(1);
        assert!(!c.has_measurement());
        c.apply(Gate::reset(0));
        assert!(c.has_measurement());
    }
}
