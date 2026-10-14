use super::*;
use std::f64::consts::PI;
fn angle(source: &str) -> crate::angle::Angle {
    crate::angle_expr::parse(source, 1).unwrap()
}

fn gates(source: &str) -> Vec<Gate> {
    parse(source).unwrap().gates
}
fn fails(source: &str, expected: &str) {
    let error = parse(source).unwrap_err();
    assert!(
        error.contains(expected),
        "expected {expected:?}, got {error:?}"
    );
}

#[test]
fn nested_parameters_are_scoped_and_not_periodically_reduced() {
    let source = "gate half(theta) a,b { rz(theta/2) b; cx a,b; }
        gate outer(phi,theta) a,b { half(phi-theta) b,a; half(theta) a,b; }
        qreg unused[2]; qreg data[2]; outer(3*pi,pi) data[0],data[1];";
    assert_eq!(
        gates(source),
        vec![
            Gate::rz(angle("pi"), 2),
            Gate::cnot {
                control: 3,
                target: 2
            },
            Gate::rz(angle("pi/2"), 3),
            Gate::cnot {
                control: 2,
                target: 3
            }
        ]
    );
}

#[test]
fn arbitrary_custom_arity_empty_calls_and_flat_output() {
    let c = parse(
        "gate empty() a,b,c,d {} gate four a,b,c,d { ccx a,b,c; ccz b,c,d; }
        qreg q[4]; empty q[0],q[1],q[2],q[3]; four() q[3],q[2],q[1],q[0];",
    )
    .unwrap();
    assert_eq!(
        c.gates,
        vec![
            Gate::ccx {
                control1: 3,
                control2: 2,
                target: 1
            },
            Gate::ccz {
                control1: 2,
                control2: 1,
                target: 0
            }
        ]
    );
    let output = serialize(&c);
    assert!(!output.contains("gate "));
    assert_eq!(parse(&output).unwrap().gates, c.gates);
    fails(
        "qreg q[1]; gate empty a {} empty q[0]; qreg later[1];",
        "declaration after gate",
    );
}

#[test]
fn broadcasting_preserves_lane_order_and_register_offsets() {
    assert_eq!(
        gates(
            "gate pair a,b { h a; cx a,b; } qreg pad[1]; qreg a[2]; qreg b[2];
        pair a,b; pair a,pad[0];"
        ),
        vec![
            Gate::h(1),
            Gate::cnot {
                control: 1,
                target: 3
            },
            Gate::h(2),
            Gate::cnot {
                control: 2,
                target: 4
            },
            Gate::h(1),
            Gate::cnot {
                control: 1,
                target: 0
            },
            Gate::h(2),
            Gate::cnot {
                control: 2,
                target: 0
            }
        ]
    );
    assert_eq!(
        gates("qreg a[2]; qreg b[2]; h a; CX a,b; U(0,0,pi) b;").len(),
        10
    );
    fails(
        "gate pair a,b {} qreg a[1]; qreg b[2]; pair a,b;",
        "size mismatch",
    );
    fails("gate pair a,b {} qreg a[2]; pair a,a[0];", "distinct");
    // A scalar may overlap a later lane's register qubit, provided it is a
    // different argument in each instantiated call.
    fails("qreg q[2]; CX q[0],q;", "distinct");
    assert!(gates("gate pair a,b {} qreg a[0]; qreg b[0]; pair a,b;").is_empty());
    fails(
        "gate pair(theta) a {} qreg q[0]; pair(1/0) q;",
        "division by zero",
    );
    fails(
        "gate pair a,b {} qreg a[0]; qreg b[1]; pair a,b;",
        "size mismatch",
    );
}

#[test]
fn definitions_validate_unused_bodies_and_signatures() {
    for (source, message) in [
        ("gate bad a { h a[0]; }", "indexed qubit formal"),
        ("qreg global[1]; gate bad a { h global; }", "qubit formal"),
        ("gate bad(theta) a { rz(phi) a; }", "unknown parameter"),
        ("gate bad a { rz(1/0) a; }", "division by zero"),
        ("gate bad(theta) a { rz(theta*(1e999999)) a; }", "limit"),
        ("gate bad a { CX a,a; }", "distinct"),
        ("gate bad a { CX a; }", "expects 2 qubit"),
        ("gate bad a { U(0,0) a; }", "expects 3 parameters"),
        (
            "gate bad a { later a; } gate later a {}",
            "forward reference",
        ),
        ("gate bad a { bad a; }", "recursive"),
        ("gate bad a { gate nested b {} }", "nested definitions"),
        ("gate bad a { barrier a; }", "inside gate"),
        ("gate bad a { reset a; }", "inside gate"),
        ("gate bad a { measure a -> c; }", "inside gate"),
        ("gate bad a { if(c==1) h a; }", "inside gate"),
        ("gate bad a { qreg q[1]; }", "inside gate"),
        ("gate bad() {}", "at least one qubit"),
        ("gate bad(a) a {}", "overlapping"),
        ("gate bad(a,a) q {}", "duplicate formal"),
        ("gate bad a,a {}", "duplicate formal"),
        ("gate bad(pi) a {}", "invalid formal"),
        ("gate bad a { h a }", "missing ';'"),
        ("gate bad a { h a;", "missing closing"),
        ("gate bad a;", "missing '{'"),
        ("gate bad a {} gate bad a {}", "conflicting gate"),
        ("gate h a {}", "conflicting gate"),
        ("qreg bad[1]; gate bad a {}", "conflicting gate"),
        ("gate bad a {} creg bad[1];", "conflicting register"),
    ] {
        fails(source, message);
    }
}

#[test]
fn call_errors_are_diagnostic_and_never_partially_return_a_circuit() {
    for (source, message) in [
        (
            "gate g(theta) a {} qreg q[1]; g q[0];",
            "expects 1 parameters",
        ),
        ("gate g a {} qreg q[1]; g(1) q[0];", "expects 0 parameters"),
        ("gate g a {} qreg q[1]; g q[0],q[0];", "expects 1 qubit"),
        ("gate g a,b {} qreg q[1]; g q[0],q[0];", "distinct"),
        ("gate g a {} qreg q[1]; g q[1];", "out of range"),
        (
            "gate g a {} qreg q[1]; g q[0]tail;",
            "invalid qubit operand",
        ),
        ("gate g a {} qreg q[1]; g missing[0];", "unknown register"),
        ("gate g a {} qreg q[1]; g q[-1];", "bad qubit index"),
        (
            "gate g(x) a {} qreg q[1]; g(pi_foo) q[0];",
            "unsupported function or constant",
        ),
        ("qreg q[1]; rz(1e999999) q[0];", "limit"),
        ("qreg q[1]; rz(0/0) q[0];", "division by zero"),
        ("include \"other.inc\";", "unsupported include"),
        ("include \"qelib1.inc;h q[0]\";", "unsupported include"),
        ("include \"unterminated", "unterminated string"),
        ("OPENQASM 3.0;", "unsupported OPENQASM"),
        ("gatekeeper g a {}", "unexpected brace"),
        ("qreg q[1]; if(c==1) h q[0];", "unsupported: if"),
        ("opaque g a;", "unsupported: opaque"),
    ] {
        fails(source, message);
    }
    let err = parse(
        "gate inner(theta) a { rz(1/theta) a; }
        gate outer(x) b { inner(x) b; }
        qreg q[1]; outer(0) q[0];",
    )
    .unwrap_err();
    assert!(err.contains("division by zero"));
    assert!(err.contains("called 'inner'"));
    assert!(err.contains("called 'outer'"));
    assert!(err.contains("defined at line 1"));
}

#[test]
fn comments_separate_tokens_and_track_body_lines() {
    let c = parse(
        "// { ; /* ignored\ngate /* comment */ g(theta) a {
        rz( // ; } ignored
            theta /* embedded */ / 2) a; // more ;
        h a;
    }qreg q[1]; g(2e-1) q[0];",
    )
    .unwrap();
    assert_eq!(c.gates, vec![Gate::rz(angle("0.1"), 0), Gate::h(0)]);
    fails(
        "gate g(theta) a { rz(theta// hi\n theta) a; }",
        "unexpected token",
    );
    let err = parse("/*\ncomment\n*/gate g a {\n  h wrong;\n}").unwrap_err();
    assert!(err.starts_with("line 4:"), "{err}");
    fails("gate g a { /* unfinished", "unterminated block comment");
}

#[test]
fn expansion_limits_include_depth_total_work_and_empty_chains() {
    let limits = ParseLimits {
        expanded_gates: 2,
        expansion_work: 10,
        call_depth: 2,
        ..ParseLimits::default()
    };
    let input = "gate g a { h a; h a; } qreg q[2]; g q;";
    assert!(
        parse_with_limits(input, limits)
            .unwrap_err()
            .contains("expanded gate limit")
    );
    let input = "gate g a { h a; } qreg q[1]; g q[0]; g q[0]; g q[0];";
    assert!(
        parse_with_limits(input, limits)
            .unwrap_err()
            .contains("expanded gate limit")
    );
    let input = "gate g a {} gate gg a { g a; } gate ggg a { gg a; } qreg q[1]; ggg q[0];";
    assert!(
        parse_with_limits(input, limits)
            .unwrap_err()
            .contains("call depth")
    );
    let mut empty = "gate empty a {}\n".to_owned();
    let mut previous = "empty".to_owned();
    for n in 0..60 {
        let name = format!("empty{n}");
        empty.push_str(&format!(
            "gate {name} a {{ {previous} a; {previous} a; }}\n"
        ));
        previous = name;
    }
    empty.push_str(&format!("qreg q[1]; {previous} q[0];"));
    assert!(parse(&empty).unwrap_err().contains("work limit"));
    // Unused huge definitions do not consume the execution budget.
    assert!(parse(empty.split("qreg").next().unwrap()).is_ok());
    // Saturation must not admit exponential work when callers raise limits.
    let overflow = empty.replace(
        "qreg q[1]; empty59 q[0];",
        "gate beyond a { empty59 a; empty59 a; }
        gate further a { beyond a; beyond a; } gate overflow a { further a; further a; }
        gate overflowed a { overflow a; overflow a; } qreg q[1]; overflowed q[0];",
    );
    let generous = ParseLimits {
        expansion_work: usize::MAX,
        ..ParseLimits::default()
    };
    assert!(
        parse_with_limits(&overflow, generous)
            .unwrap_err()
            .contains("overflow")
    );
    let unsafe_depth = ParseLimits {
        expression_depth: usize::MAX,
        ..ParseLimits::default()
    };
    assert!(parse_with_limits("qreg q[1];", unsafe_depth).is_err());
    let deep = format!(
        "qreg q[1]; rz({}pi{}) q[0];",
        "(".repeat(300),
        ")".repeat(300)
    );
    fails(&deep, "Depth");
    let unary = format!("qreg q[1]; rz({}pi) q[0];", "-".repeat(300));
    fails(&unary, "Depth");
    // No custom operation cap on ordinary flat input.
    assert!(parse_with_limits("qreg q[1]; h q[0]; h q[0]; h q[0];", limits).is_ok());
}

#[test]
fn native_broadcasts_obey_output_and_work_limits() {
    let limits = ParseLimits {
        expanded_gates: 1,
        expansion_work: 1,
        ..Default::default()
    };
    for source in [
        "qreg q[2]; h q;",
        "qreg q[2]; U(0,0,0) q;",
        "qreg q[2]; u2(0,0) q;",
        "qreg q[2]; u3(0,0,0) q;",
        "qreg q[2]; rz(0) q;",
        "qreg a[2]; qreg b[2]; CX a,b;",
        "qreg q[3]; CX q[0],q;",
    ] {
        let error = parse_with_limits(source, limits).unwrap_err();
        assert!(error.contains("expanded gate limit"), "{source}: {error}");
    }
    // Identity broadcasts emit no gates but still visit every lane.
    let error = parse_with_limits("qreg q[1000]; id q;", limits).unwrap_err();
    assert!(error.contains("expansion work limit"));

    // Reject these tiny inputs before reservation or billions of loop visits.
    fails("qreg q[4294967295]; U(0,0,0) q;", "expanded gate limit");
    fails("qreg q[4294967295]; id q;", "expansion work limit");
}

#[test]
fn implicit_expansion_budgets_are_shared_and_allow_exact_boundaries() {
    let source = "gate g a { h a; } qreg q[2]; h q; g q[0];";
    let limits = ParseLimits {
        expanded_gates: 3,
        expansion_work: 8,
        ..Default::default()
    };
    assert_eq!(parse_with_limits(source, limits).unwrap().gates.len(), 3);
    let short_output = ParseLimits {
        expanded_gates: 2,
        ..limits
    };
    assert!(
        parse_with_limits(source, short_output)
            .unwrap_err()
            .contains("expanded gate limit")
    );
    let short_work = ParseLimits {
        expansion_work: 7,
        ..limits
    };
    assert!(
        parse_with_limits(source, short_work)
            .unwrap_err()
            .contains("expansion work limit")
    );
    let singles = "qreg q[1]; h q; h q;";
    assert!(parse_with_limits(singles, short_output).is_ok());
    assert!(
        parse_with_limits(
            singles,
            ParseLimits {
                expanded_gates: 1,
                ..limits
            }
        )
        .is_err()
    );

    // Explicit indexed native statements remain uncapped, including aliases
    // that lower to multiple gates and calls with optional empty parentheses.
    let zero = ParseLimits {
        expanded_gates: 0,
        expansion_work: 0,
        ..limits
    };
    let flat = parse_with_limits("qreg q[1]; h q[0]; h() q[0]; U(0,0,0) q[0];", zero).unwrap();
    assert_eq!(flat.gates.len(), 5);
    assert!(
        parse_with_limits("qreg q[0]; h q; U(0,0,0) q; id q;", zero)
            .unwrap()
            .gates
            .is_empty()
    );
}

#[test]
fn large_empty_gate_broadcasts_charge_operand_binding() {
    let arity = 1000;
    let formals = (0..arity)
        .map(|i| format!("a{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let actuals = (0..arity - 1)
        .map(|i| format!("q[{i}]"))
        .chain(std::iter::once("r".to_owned()))
        .collect::<Vec<_>>()
        .join(",");
    let source = |lanes| {
        format!(
            "gate g {formals} {{}} qreg q[{}]; qreg r[{lanes}]; g {actuals};",
            arity - 1
        )
    };
    // The audit reproducer formerly fit a budget of one frame visit per lane.
    let limits = ParseLimits {
        expanded_gates: 0,
        expansion_work: 20_000,
        ..Default::default()
    };
    assert!(
        parse_with_limits(&source(20_000), limits)
            .unwrap_err()
            .contains("expansion work limit")
    );
    // Its dangerous default-budget version must also fail before visiting lanes.
    fails(&source(20_000_000), "expansion work limit");
    // Keep arbitrary arity supported with exact accounting: one frame + arity
    // bound qubits for each lane, even though the gate emits nothing.
    let exact = ParseLimits {
        expanded_gates: 0,
        expansion_work: 2 * (arity + 1),
        ..Default::default()
    };
    assert!(
        parse_with_limits(&source(2), exact)
            .unwrap()
            .gates
            .is_empty()
    );
    assert!(
        parse_with_limits(
            &source(2),
            ParseLimits {
                expansion_work: exact.expansion_work - 1,
                ..exact
            }
        )
        .is_err()
    );
}

#[test]
fn descriptor_alias_checks_match_explicit_lane_checks() {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(219);
    for _ in 0..1000 {
        let width = rng.gen_range(0..6);
        let count = rng.gen_range(1..9);
        let formals = (0..count)
            .map(|i| format!("a{i}"))
            .collect::<Vec<_>>()
            .join(",");
        let mut operands = vec!["q".to_owned()];
        let mut addresses = vec![(0, true)];
        for _ in 1..count {
            match rng.gen_range(0..3) {
                0 => {
                    operands.push("q".to_owned());
                    addresses.push((0, true));
                }
                1 => {
                    operands.push("r".to_owned());
                    addresses.push((width, true));
                }
                _ => {
                    let q = rng.gen_range(0..2 * width + 1);
                    let (name, index) = if q < width {
                        ("q", q)
                    } else if q < 2 * width {
                        ("r", q - width)
                    } else {
                        ("s", 0)
                    };
                    operands.push(format!("{name}[{index}]"));
                    addresses.push((q, false));
                }
            }
        }
        let aliased = (0..width).any(|lane| {
            let qubits: Vec<_> = addresses
                .iter()
                .map(|&(q, register)| q + if register { lane } else { 0 })
                .collect();
            qubits
                .iter()
                .enumerate()
                .any(|(i, q)| qubits[..i].contains(q))
        });
        let source = format!(
            "gate g {formals} {{}} qreg q[{width}]; qreg r[{width}]; qreg s[1]; g {};",
            operands.join(",")
        );
        let parsed = parse(&source);
        assert_eq!(parsed.is_err(), aliased, "{source}: {parsed:?}");
        if let Err(error) = parsed {
            assert!(error.contains("distinct"), "{error}");
        }
    }
    // Detect an alias at the very last lane, before reservation/expansion.
    fails("gate g a,b {} qreg q[1000000]; g q,q[999999];", "distinct");
    // Equal-width registers have different qubits at the same lane; fixed
    // operands just beyond a register's end must also remain valid.
    assert!(parse("gate g a,b {} qreg q[2]; qreg r[2]; g q,r; g q,r[0];").is_ok());
}

#[test]
fn nested_work_budget_counts_parameter_expression_execution() {
    // inner: frame + parameter + qubit = 3.
    // outer: frame + parameter + qubit + inner + five expression opcodes = 11.
    let source = "gate inner(y) a {} gate outer(x) a { inner(x+x+x) a; } qreg q[1]; outer(1) q[0];";
    let limits = ParseLimits {
        expanded_gates: 0,
        expansion_work: 11,
        ..Default::default()
    };
    assert!(parse_with_limits(source, limits).unwrap().gates.is_empty());
    assert!(
        parse_with_limits(
            source,
            ParseLimits {
                expansion_work: 10,
                ..limits
            }
        )
        .unwrap_err()
        .contains("expansion work limit")
    );
    // Cached constants count as one step regardless of their source size.
    let constant = source.replace("x+x+x", "1+2+3");
    assert!(
        parse_with_limits(
            &constant,
            ParseLimits {
                expansion_work: 7,
                ..limits
            }
        )
        .is_ok()
    );
}

#[test]
fn malformed_source_fuzz_never_panics() {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(9173);
    let alphabet = [
        'g', 'a', 't', 'e', 'q', 'r', 'z', '(', ')', '[', ']', '{', '}', ';', '/', '*', '\n', ' ',
        ',', '"', '\\', '0', '1', '.', '-', 'π', '🦀',
    ];
    for _ in 0..2000 {
        let source: String = (0..rng.gen_range(0..200))
            .map(|_| alphabet[rng.gen_range(0..alphabet.len())])
            .collect();
        let _ = parse(&source);
    }
}

#[test]
fn random_nested_binding_matches_independent_flat_construction() {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(712);
    for _ in 0..100 {
        let theta = rng.gen_range(-20.0..20.0);
        let phi = rng.gen_range(-20.0..20.0);
        let source = format!(
            "gate inner(x,y) a,b {{ rz(x/2) b; cx a,b; rz(-y/3) a; }}
            gate outer(y,x) a,b {{ inner(x-y,y+x) b,a; }}
            qreg pad[1]; qreg q[2]; outer({theta},{phi}) q[0],q[1];"
        );
        let c = parse(&source).unwrap();
        assert_eq!(
            c.gates,
            vec![
                Gate::rz(angle(&format!("({phi}-({theta}))/2")), 1),
                Gate::cnot {
                    control: 2,
                    target: 1
                },
                Gate::rz(angle(&format!("-(({theta})+({phi}))/3")), 2)
            ]
        );
    }
}

#[test]
fn u_lowering_matches_independent_analytic_matrix() {
    for (theta, phi, lambda) in [
        (0.0, 0.0, 0.0),
        (PI, 0.0, PI),
        (-0.37, 1.4, -2.3),
        (8.0 * PI, -4.0 * PI, PI / 3.0),
        (0.2, 0.3, 0.4),
    ] {
        let c = parse(&format!(
            "gate g(t,p,l) a {{ U(t,p,l) a; }} qreg q[1]; g({theta},{phi},{lambda}) q[0];"
        ))
        .unwrap();
        let actual = crate::unitary::circuit_unitary(&c);
        let angles = [
            -(phi + lambda) / 2.0,
            -(phi - lambda) / 2.0,
            (phi - lambda) / 2.0,
            (phi + lambda) / 2.0,
        ];
        let scales = [
            (theta / 2.0).cos(),
            -(theta / 2.0).sin(),
            (theta / 2.0).sin(),
            (theta / 2.0).cos(),
        ];
        for i in 0..4 {
            let (re, im) = actual[i / 2][i % 2].components();
            assert!((re - scales[i] * angles[i].cos()).abs() < 1e-12);
            assert!((im - scales[i] * angles[i].sin()).abs() < 1e-12);
        }
    }
    assert_eq!(
        gates("qreg q[1]; u2(0,pi) q[0];"),
        gates("qreg q[1]; U(pi/2,0,pi) q[0];")
    );
    assert_eq!(
        gates("qreg q[1]; u3(0.1,0.2,0.3) q[0];"),
        gates("qreg q[1]; U(0.1,0.2,0.3) q[0];")
    );
    assert_eq!(
        gates("qreg q[1]; u1(pi) q[0]; id q[0];"),
        vec![Gate::p(angle("pi"), 0)]
    );
}

#[test]
fn custom_controlled_rotation_retains_relative_phase() {
    let c = parse(
        "gate controlled(theta) a,b { rz(theta/2) b; cx a,b; rz(-theta/2) b; cx a,b; }
        qreg q[2]; controlled(2*pi) q[0],q[1];",
    )
    .unwrap();
    let expected = parse("qreg q[2]; z q[0];").unwrap();
    assert!(crate::unitary::circuits_equiv(&c, &expected, 1e-12));
    assert!(!crate::unitary::circuits_equiv(&c, &Circuit::new(2), 1e-12));
}

#[test]
fn inlining_and_optimization_preserve_exact_measurement_channels() {
    use crate::semantics::channel::{ChannelLimits, circuit_channel};
    let source = "gate pair a,b { h a; cx a,b; t b; tdg b; }
        gate twice a,b { pair a,b; pair a,b; }
        qreg q[2]; creg c[1]; twice q[0],q[1]; measure q[1] -> c[0];";
    let c = parse(source).unwrap();
    let reference = parse(
        "qreg q[2]; creg c[1]; h q[0]; cx q[0],q[1]; t q[1]; tdg q[1];
        h q[0]; cx q[0],q[1]; t q[1]; tdg q[1]; measure q[1] -> c[0];",
    )
    .unwrap();
    let options = crate::optimize::Options {
        level: crate::optimize::Level::O1,
        ..Default::default()
    };
    let (optimized, _) = crate::optimize::optimize(&c, &options).unwrap();
    for initial in [false, true] {
        let expected = circuit_channel(&reference, &[initial], ChannelLimits::default()).unwrap();
        for actual in [&c, &optimized] {
            let actual = circuit_channel(actual, &[initial], ChannelLimits::default()).unwrap();
            assert_eq!(expected.compare(&actual), Ok(()));
        }
    }
}

#[test]
fn surrounding_measurement_and_reset_survive_inlining() {
    assert_eq!(
        gates(
            "gate pulse a { x a; } qreg q[1]; creg c[1]; pulse q[0];
        measure q -> c; reset q; pulse q[0];"
        ),
        vec![
            Gate::x(0),
            Gate::measure { qubit: 0, cbit: 0 },
            Gate::reset(0),
            Gate::x(0)
        ]
    );
}

#[test]
fn exact_symbolic_bindings_preserve_main_scalar_semantics() {
    let error = parse("gate bad a { rz(1/0) a; }").unwrap_err();
    assert_eq!(error.matches("line 1:").count(), 1, "{error}");
    let definitions = "gate inner(t) a { rz(t) a; }
        gate outer(x) a { inner(-x^2+sqrt(4)+sin(pi/2)) a; inner(2^3^2) a; }
        qreg q[1]; outer(0.1+0.2) q[0];";
    let c = parse(definitions).unwrap();
    assert_eq!(
        c.gates,
        vec![Gate::rz(angle("-(3/10)^2+3"), 0), Gate::rz(angle("512"), 0)]
    );
    let c = parse("gate g(theta) a { rz(theta/2) a; } qreg q[1]; g(pi/7) q[0];").unwrap();
    assert_eq!(c.gates, vec![Gate::rz(angle("pi/14"), 0)]);
    let expression = "gate g(theta) a { rz(sin(theta)+sqrt(2)/3) a; } qreg q[1]; g(pi/7) q[0];";
    let c = parse(expression).unwrap();
    assert_eq!(c.gates, vec![Gate::rz(angle("sin(pi/7)+sqrt(2)/3"), 0)]);
    assert_eq!(parse(&serialize(&c)).unwrap().gates, c.gates);
    fails(
        "gate g(x) a { rz(ln(x)) a; } qreg q[1]; g(-1) q[0];",
        "real domain",
    );
}

#[test]
fn custom_bodies_keep_all_current_native_gates_and_literal_controls() {
    let body = "p(t) a; y a; sx a; rx(t) a; ry(t) a; swap a,b; cy a,b; ch a,b;
        cp(t) a,b; crx(t) a,b; cry(t) a,b; crz(t) a,b; cswap a,b,c;";
    let custom = parse(&format!(
        "gate g(t) a,b,c {{ {body} }} qreg q[3]; g(2*pi) q[0],q[1],q[2];"
    ))
    .unwrap();
    let expected = parse(
        "qreg q[3]; p(2*pi) q[0]; y q[0]; sx q[0]; rx(2*pi) q[0]; ry(2*pi) q[0];
        swap q[0],q[1]; cy q[0],q[1]; ch q[0],q[1]; cp(2*pi) q[0],q[1];
        crx(2*pi) q[0],q[1]; cry(2*pi) q[0],q[1]; crz(2*pi) q[0],q[1]; cswap q[0],q[1],q[2];",
    )
    .unwrap();
    assert_eq!(custom.gates, expected.gates);
    assert!(matches!(&custom.gates[11], Gate::crz { theta, .. } if theta == &angle("2*pi")));
}

#[test]
fn sdk_literal_policy_survives_inlining() {
    let source = "qreg q[1]; rz(0.1) q[0];";
    let sdk = parse_numeric_sdk(source).unwrap();
    assert_eq!(
        sdk.gates,
        vec![Gate::rz(crate::angle::Angle::from_f64(0.1).unwrap(), 0)]
    );
    assert_ne!(sdk.gates, parse(source).unwrap().gates);
    assert!(parse_numeric_sdk("qreg q[1]; rz(0.1+0.2) q[0];").is_err());
    let custom = parse_numeric_sdk("gate g(t) a { rz(t/2) a; } qreg q[1]; g(0.1) q[0];").unwrap();
    let expected = crate::angle::AngleContext::default()
        .half(&crate::angle::Angle::from_f64(0.1).unwrap())
        .unwrap();
    assert_eq!(custom.gates, vec![Gate::rz(expected, 0)]);
}

#[test]
fn sdk_custom_bodies_use_mathematical_constants_and_bound_binary64_values() {
    let source = "gate inner(t) a { rz(t+pi/2) a; rz(t+0.1) a; }
        gate constants a { rz(sin(pi/2)+0.1) a; }
        gate outer(t) a { inner(t/2) a; constants a; }
        qreg q[1]; outer(0.1) q[0]; rz(0.1) q[0];";
    let custom = parse_numeric_sdk(source).unwrap();
    let mut context = crate::angle::AngleContext::default();
    let input = crate::angle::Angle::from_f64(0.1).unwrap();
    let half = context.half(&input).unwrap();
    let shifted_pi = context.add(&half, &angle("pi/2")).unwrap();
    let shifted_decimal = context.add(&half, &angle("1/10")).unwrap();
    assert_eq!(
        custom.gates,
        vec![
            Gate::rz(shifted_pi, 0),
            Gate::rz(shifted_decimal, 0),
            Gate::rz(angle("11/10"), 0),
            Gate::rz(input, 0),
        ]
    );
    // Mathematical body constants must not bypass the SDK input policy via
    // the shared scalar cache, even when the same source occurs at top level.
    for source in ["pi", "0.1+0.2"] {
        let error = parse_numeric_sdk(&format!(
            "gate constants a {{ rz({source}) a; }} qreg q[1]; rz({source}) q[0];"
        ))
        .unwrap_err();
        assert!(error.contains("SDK parameter must be a numeric literal"));
    }
}

#[test]
fn qasm2_lexical_rules_reject_signed_integers_reserved_names_and_truncation() {
    for (source, message) in [
        ("qreg q[2]; h q[+1];", "bad qubit index '+1'"),
        ("qreg a[1]; qreg q[2]; h q[+1];", "bad qubit index '+1'"),
        (
            "qreg q[1]; creg c[2]; measure q[0] -> c[+1];",
            "bad cbit index '+1'",
        ),
        ("qreg q[+2];", "bad qreg size '+2'"),
        ("creg c[+2];", "bad creg size '+2'"),
        ("gate g(sin) a { rz(sin) a; }", "invalid formal name 'sin'"),
        ("gate g a,sqrt { CX a,sqrt; }", "invalid formal name 'sqrt'"),
        ("gate exp a { h a; }", "invalid gate name 'exp'"),
        ("qreg ln[1];", "invalid qreg declaration"),
        // A file cut off mid-statement must not broadcast over the register.
        ("qreg q[4]; h q", "missing ';'"),
        ("qreg q[4]; h q[3]", "missing ';'"),
    ] {
        fails(source, message);
    }
    assert_eq!(gates("qreg q[2]; h q[ 1 ];"), vec![Gate::h(1)]);
    assert_eq!(gates("qreg q[4]; h q[3]; // trailing\n"), vec![Gate::h(3)]);
    assert_eq!(
        gates("gate g(sine) a { rz(sine) a; } qreg q[1]; g(pi/4) q[0];"),
        gates("qreg q[1]; rz(pi/4) q[0];")
    );
}

#[test]
fn non_integer_exponents_bind_exactly_in_gate_bodies() {
    let expected = gates("qreg q[1]; rz(sqrt(2)) q[0]; rz(exp(ln(3)/4)) q[0];");
    assert_eq!(
        gates("qreg q[1]; rz(2^0.5) q[0]; rz(3^0.25) q[0];"),
        expected
    );
    assert_eq!(
        gates("gate g(a,b) x { rz(a^0.5) x; rz(b^(1/4)) x; } qreg q[1]; g(2,3) q[0];"),
        expected
    );
    fails(
        "gate g(a) x { rz(a^0.5) x; } qreg q[1]; g(-1) q[0];",
        "non-integer power of a non-positive base",
    );
    fails("qreg q[1]; rz((-2)^0.5) q[0];", "non-integer power");
}

#[test]
fn reports_describe_the_emitted_statements_for_every_native_gate() {
    let source = "qreg q[3]; creg c[1];
        h q[0]; x q[1]; y q[2]; z q[0]; s q[1]; sdg q[2]; t q[0]; tdg q[1]; sx q[2];
        rz(0.3) q[0]; rx(pi/4) q[1]; ry(0.2) q[2]; p(pi/8) q[0];
        cx q[0],q[1]; cz q[1],q[2]; cy q[2],q[0]; ch q[0],q[2]; swap q[1],q[0];
        cp(0.4) q[0],q[1]; crx(pi/2) q[1],q[2]; cry(0.5) q[2],q[0]; crz(pi) q[0],q[2];
        ccx q[0],q[1],q[2]; ccz q[2],q[1],q[0]; cswap q[1],q[2],q[0];
        measure q[0] -> c[0]; reset q[1];";
    let circuit = parse(source).unwrap();
    let emitted = serialize(&circuit);
    let reparsed = parse(&emitted).unwrap();
    assert_eq!(
        crate::optimize::Metrics::of(&circuit),
        crate::optimize::Metrics::of(&reparsed)
    );
    assert_eq!(
        crate::optimize::Metrics::of(&circuit).gates,
        emitted.lines().filter(|l| l.contains("q[")).count() - 1
    );
    let mut names: Vec<_> = emitted
        .lines()
        .skip(4)
        .map(|l| scan::identifier(l).0)
        .collect();
    names.sort_unstable();
    names.dedup();
    let mut reported = emitted_gate_names(&circuit);
    reported.sort_unstable();
    assert_eq!(reported, names);
    for absent in ["p", "cp", "sx", "swap", "crx", "cry", "cswap"] {
        assert!(!reported.contains(&absent), "{absent}");
    }
}
