//! Independent formal phase oracle: function subexpressions are indeterminates.
//! This proves preservation for every value of those functions, without asking
//! the production normalizer or a floating evaluator to establish equivalence.
use super::{Angle, AngleContext, Limits, Rational, view::ExpressionRef};
use crate::{
    cancel::CancelGates,
    circuit::{Circuit, Gate},
    cnot_min::CnotMin,
    pass::{Pass, PassState},
    phase_fold_pauli::PhaseFoldPauli,
    phase_fold_rand::PhaseFoldRand,
};
use num_rational::BigRational as R;
use num_traits::Zero;
use std::collections::BTreeMap;
type Phase = BTreeMap<String, R>;
fn scalar(a: &Angle, scale: &R, out: &mut Phase) {
    match a.expression() {
        ExpressionRef::Rational(r) => {
            *out.entry("1".into()).or_default() += scale * r.to_big_rational()
        }
        ExpressionRef::Pi(r) => *out.entry("pi".into()).or_default() += scale * r.to_big_rational(),
        ExpressionRef::Sum(terms) => {
            for (c, a) in terms {
                scalar(a, &(scale * c.to_big_rational()), out);
            }
        }
        ExpressionRef::Product {
            coefficient,
            factors,
        } => {
            let scale = scale * coefficient.to_big_rational();
            if factors.len() == 1 {
                scalar(&factors[0], &scale, out);
            } else {
                let key = format!(
                    "product({})",
                    factors
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join("*")
                );
                *out.entry(key).or_default() += scale;
            }
        }
        _ => *out.entry(a.to_string()).or_default() += scale,
    }
}
pub(super) fn signatures(c: &Circuit) -> Vec<(usize, Phase)> {
    let mut result = Vec::new();
    for input in 0..1usize << c.num_qubits {
        let mut bits = input;
        let mut phase = Phase::new();
        for g in &c.gates {
            let named = match g {
                Gate::t(q) => Some((*q, 1)),
                Gate::tdg(q) => Some((*q, -1)),
                Gate::s(q) => Some((*q, 2)),
                Gate::sdg(q) => Some((*q, -2)),
                Gate::z(q) => Some((*q, 4)),
                _ => None,
            };
            if let Some((q, k)) = named {
                if bits >> q & 1 != 0 {
                    *phase.entry("pi".into()).or_default() += R::new(k.into(), 4.into());
                }
                continue;
            }
            match g {
                Gate::x(q) => bits ^= 1 << q,
                Gate::cnot { control, target } => {
                    if bits >> control & 1 != 0 {
                        bits ^= 1 << target
                    }
                }
                Gate::cz { control, target } => {
                    if bits >> control & 1 != 0 && bits >> target & 1 != 0 {
                        *phase.entry("pi".into()).or_default() += R::from_integer(1.into());
                    }
                }
                Gate::rz(a, q) | Gate::p(a, q) => {
                    if bits >> q & 1 != 0 {
                        scalar(a, &R::from_integer(1.into()), &mut phase)
                    }
                }
                _ => panic!("unexpected monomial gate {g:?}"),
            }
        }
        result.push((bits, phase));
    }
    let origin = result[0].1.clone();
    for (_, phase) in &mut result {
        for (atom, c) in &origin {
            *phase.entry(atom.clone()).or_default() -= c;
        }
        if let Some(p) = phase.get_mut("pi") {
            *p = R::new(p.numer() % (p.denom() * 2), p.denom().clone());
            if *p < R::zero() {
                *p += R::from_integer(2.into());
            }
        }
        phase.retain(|_, c| !c.is_zero());
    }
    result
}
pub(super) fn passes() -> Vec<Box<dyn Pass>> {
    vec![
        Box::new(CancelGates),
        Box::new(PhaseFoldRand),
        Box::new(PhaseFoldPauli),
        Box::new(CnotMin::default()),
    ]
}
#[test]
fn general_expressions_fold_across_complemented_and_parity_axes() {
    for body in [
        "rz(sin(1)+sqrt(2)/7) q[0]; rz(cos(2)) q[0]; rz(-sin(1)-sqrt(2)/7) q[0];",
        "rz(exp(1)+ln(2)) q[0]; x q[0]; rz(exp(1)+ln(2)) q[0]; x q[0];",
        "cx q[0],q[1]; rz(sin(1)/3) q[1]; cx q[0],q[1]; cx q[0],q[1]; rz(2*sin(1)/3) q[1]; cx q[0],q[1];",
        "rz(tan(1)/3+pi/7) q[0]; rz(-tan(1)/3+pi/7) q[0];",
        "rz(sin(1)) q[0]; rz(sin(1)) q[0]; rz(-2*sin(1)) q[0];",
    ] {
        let input = crate::qasm::parse(&format!("qreg q[2]; {body}")).unwrap();
        let reference = signatures(&input);
        for p in passes() {
            let output = p.run(&input);
            assert_eq!(signatures(&output), reference, "{}: {body}", p.name());
            assert_eq!(
                signatures(&crate::qasm::parse(&output.to_qasm()).unwrap()),
                reference
            );
        }
    }
}
#[test]
fn limit_failures_preserve_original_gates_and_survive_context_drop() {
    let input = crate::qasm::parse("qreg q[1]; rz(sin(1)) q[0]; rz(cos(2)) q[0];").unwrap();
    for p in passes() {
        let mut state = PassState::with_limits(Limits {
            max_terms: 1,
            ..Limits::default()
        });
        let output = p.run_with_state(&input, &mut state);
        drop(state);
        assert_eq!(output.gates, input.gates, "{}", p.name());
        assert_eq!(signatures(&output), signatures(&input));
        let mut state = PassState::with_limits(Limits {
            max_work: 0,
            ..Limits::default()
        });
        assert_eq!(
            p.run_with_state(&input, &mut state).gates,
            input.gates,
            "{} work",
            p.name()
        );
    }
}
#[test]
fn every_gate_family_cancels_symbolic_inverses_with_its_own_period() {
    for family in ["p", "rx", "ry", "rz", "cp", "crx", "cry", "crz"] {
        let operands = if family.starts_with('c') {
            "q[0],q[1]"
        } else {
            "q[0]"
        };
        let source = format!(
            "qreg q[2]; {family}(sin(1)+sqrt(2)) {operands}; {family}(-sin(1)-sqrt(2)) {operands};"
        );
        let input = crate::qasm::parse(&source).unwrap();
        assert!(CancelGates.run(&input).gates.is_empty(), "{family}");
        for (period, cancels) in [(2, !family.starts_with("cr")), (4, true)] {
            let input = crate::qasm::parse(&format!(
                "qreg q[2]; {family}(pi) {operands}; {family}({}*pi) {operands};",
                period - 1
            ))
            .unwrap();
            assert_eq!(
                CancelGates.run(&input).gates.is_empty(),
                cancels,
                "{family} {period}pi"
            );
        }
        let input = crate::qasm::parse(&format!(
            "qreg q[2]; {family}(sin(1)) {operands}; {family}(cos(2)) {operands};"
        ))
        .unwrap();
        assert_eq!(CancelGates.run(&input).gates, input.gates);
    }
}
#[test]
fn sparse_general_accumulation_retains_exact_coefficients() {
    let mut ctx = AngleContext::default();
    let atom = ctx.parse("sin(1)").unwrap();
    let mut input = Circuit::new(1);
    for _ in 0..1000 {
        input.apply(Gate::rz(atom.clone(), 0));
    }
    let negative = ctx.scale(&atom, &Rational::integer(-999)).unwrap();
    input.apply(Gate::rz(negative, 0));
    let reference = signatures(&input);
    for p in passes().into_iter().skip(1) {
        let output = p.run(&input);
        assert_eq!(
            output.gates,
            vec![Gate::rz(atom.clone(), 0)],
            "{}",
            p.name()
        );
        assert_eq!(signatures(&output), reference);
    }
}
#[test]
fn full_symbolic_expressions_are_certified_at_synthesis_boundary() {
    for family in ["rz", "crz"] {
        let operands = if family == "crz" { "q[0],q[1]" } else { "q[0]" };
        let input = crate::qasm::parse(&format!(
            "qreg q[2]; {family}(sin(1)+cos(1)+tan(1)/7+exp(1)/9+ln(2)/11+sqrt(2)/13) {operands};"
        ))
        .unwrap();
        let output = crate::decompose::DecomposeRotations { epsilon: 1e-4 }
            .try_run(&input)
            .unwrap();
        assert!(crate::unitary::circuits_equiv(&input, &output, 1e-4));
        assert!(output.gates.iter().all(|g| g.angle().is_none()));
    }
}

#[test]
fn seeded_general_expression_streams_preserve_formal_phases() {
    let mut ctx = AngleContext::default();
    let atoms: Vec<_> = [
        "sin(1)",
        "cos(2)/3",
        "sqrt(2)/7",
        "exp(1)/11",
        "ln(2)/13",
        "tan(1)/17",
        "pi/7",
        "0.1",
    ]
    .into_iter()
    .map(|s| ctx.parse(s).unwrap())
    .collect();
    for seed in 0..64u64 {
        let mut random = seed + 1;
        let mut input = Circuit::new(3);
        for _ in 0..48 {
            random = random
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let q = ((random >> 32) % 3) as u32;
            let gate = match random % 5 {
                0 => Gate::x(q),
                1 => Gate::cnot {
                    control: q,
                    target: (q + 1) % 3,
                },
                2 => Gate::t(q),
                _ => {
                    let a = &atoms[((random >> 16) as usize) % atoms.len()];
                    let a = if random & 256 != 0 {
                        ctx.neg(a).unwrap()
                    } else {
                        a.clone()
                    };
                    Gate::rz(a, q)
                }
            };
            input.apply(gate);
        }
        let reference = signatures(&input);
        for p in passes() {
            assert_eq!(
                signatures(&p.run(&input)),
                reference,
                "seed={seed} {}",
                p.name()
            );
        }
    }
}
#[test]
fn normalized_inverses_and_pythagorean_sums_fold_to_identity() {
    let inverse =
        crate::qasm::parse("qreg q[1]; rz(2*pi*sin(1)) q[0]; rz(-2*pi*sin(1)) q[0];").unwrap();
    for p in passes() {
        assert!(p.run(&inverse).gates.is_empty(), "{}", p.name());
    }
    let pythagorean =
        crate::qasm::parse("qreg q[1]; rz(-1) q[0]; rz(sin(1)^2) q[0]; rz(cos(1)^2) q[0];")
            .unwrap();
    for p in passes().into_iter().skip(1) {
        assert!(p.run(&pythagorean).gates.is_empty(), "{}", p.name());
    }
}
#[test]
fn printed_circuits_reparse_within_the_expression_depth_limit() {
    let mut angle = "1".to_string();
    for _ in 0..16 {
        angle = format!("sin(1+{angle})");
    }
    let input = crate::qasm::parse(&format!("qreg q[1]; rz({angle}) q[0];")).unwrap();
    let output = crate::qasm::parse(&input.to_qasm()).unwrap();
    assert_eq!(output.gates, input.gates);
}
#[test]
fn differently_scaled_sum_products_cancel() {
    let inverse = crate::qasm::parse(
        "qreg q[1]; rz(2*(sin(1)+cos(1))*sin(2)) q[0]; rz(-(sin(1)+cos(1))*(2*sin(2))) q[0];",
    )
    .unwrap();
    for p in passes() {
        assert!(p.run(&inverse).gates.is_empty(), "{}", p.name());
    }
}
