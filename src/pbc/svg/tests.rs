use super::*;
use crate::circuit::{Circuit, Gate};
use crate::pbc::{PauliAngle, to_pbc};

#[test]
fn draws_boxes_letters_tabs_and_frame() {
    let circuit = Circuit {
        num_qubits: 3,
        num_cbits: 1,
        gates: vec![
            Gate::h(0),
            Gate::cnot {
                control: 0,
                target: 2,
            },
            Gate::t(2),
            Gate::tdg(1),
            Gate::measure { qubit: 2, cbit: 0 },
        ],
    };
    let pbc = to_pbc(&circuit).unwrap();
    let svg = pbc.to_svg().unwrap();
    assert!(svg.starts_with("<svg") && svg.trim_end().ends_with("</svg>"));
    // Wire labels are 1-indexed kets.
    for q in 1..=3 {
        assert!(svg.contains(&format!(r#"<tspan font-size="14" dy="5">{q}</tspan>"#)));
    }
    // T on q2 after CX from an H-rotated q0: axis X0 Z2 spans rows 0-2, with 𝟙 on row 1.
    assert!(svg.contains(">𝟙</text>"));
    assert!(svg.contains(GREEN) && svg.contains(BLUE));
    // T† is a negative π/8 rotation: a minus strip on its tab.
    assert!(svg.contains(r#"fill="white" stroke="black" stroke-width="1.1""#));
    // The output Clifford is left out by default, as in the paper, and drawn
    // as a final box C on request.
    assert!(!svg.contains(">C</text>"));
    let with_frame = pbc
        .to_svg_with(SvgOptions {
            show_frame: true,
            ..SvgOptions::default()
        })
        .unwrap();
    assert!(with_frame.contains(">C</text>"));
}

#[test]
fn angle_tabs_cover_every_angle_and_signs_fold() {
    let mut c = PbcCircuit::new(1, 0);
    let z = c.z(0).unwrap();
    for k in [1, 2, 3, 4, -1] {
        c.rotate(z, PauliAngle::new(k)).unwrap();
    }
    let svg = c.to_svg().unwrap();
    for (num, den) in [("π", "8"), ("π", "4"), ("3π", "8"), ("π", "2")] {
        assert!(svg.contains(&format!(">{num}</text>")) && svg.contains(&format!(">{den}</text>")));
    }
    assert!(svg.contains(ORANGE) && svg.contains(GREY));
    // Only the k = -1 rotation (and k = 3, shown as 3π/8 positive) carry
    // signs: exactly one minus strip.
    assert_eq!(
        svg.matches(r#"fill="white" stroke="black" stroke-width="1.1""#)
            .count(),
        1
    );
}

#[test]
fn disjoint_operations_share_a_column_and_long_circuits_truncate() {
    let mut c = PbcCircuit::new(2, 0);
    let z0 = c.z(0).unwrap();
    let z1 = c.z(1).unwrap();
    c.rotate(z0, PauliAngle::new(1)).unwrap();
    c.rotate(z1, PauliAngle::new(1)).unwrap();
    let packed = c.to_svg().unwrap();
    let unpacked = c
        .to_svg_with(SvgOptions {
            pack: false,
            reorder: false,
            ..SvgOptions::default()
        })
        .unwrap();
    let width = |svg: &str| -> f64 {
        let start = svg.find("width=\"").unwrap() + 7;
        svg[start..].split('"').next().unwrap().parse().unwrap()
    };
    assert!(width(&packed) < width(&unpacked));
    let truncated = c
        .to_svg_with(SvgOptions {
            max_operations: 1,
            ..SvgOptions::default()
        })
        .unwrap();
    assert!(truncated.contains("stroke-dasharray=\"1.6 4\""));
}

/// x coordinates of the colored (operation) boxes, in drawing order.
fn box_xs(svg: &str) -> Vec<f64> {
    svg.lines()
        .filter(|l| l.starts_with("<rect x=") && !l.contains("fill=\"white\""))
        .map(|l| {
            let start = l.find("x=\"").unwrap() + 3;
            l[start..].split('"').next().unwrap().parse().unwrap()
        })
        .collect()
}

/// The paper's Fig. 4 (bottom right): Z, XY and -Y all commute with each
/// other and with ZZZY; the paper draws the first three in one column and
/// ZZZY, whose span overlaps them, in the next. Reordering reproduces that
/// even though ZZZY comes before XY in program order.
#[test]
fn reorder_matches_the_papers_figure_4_layout() {
    let mut c = PbcCircuit::new(4, 0);
    let axis = |c: &mut PbcCircuit, factors: &[(u32, Pauli)]| {
        let mut p = c.identity().as_ref();
        for &(q, f) in factors {
            let single = c.single(q, f).unwrap();
            p = c.product(p, single.as_ref()).unwrap();
        }
        c.hermitian_axis(p, 1000).unwrap()
    };
    use Pauli::{X, Y, Z};
    let z0 = axis(&mut c, &[(0, Z)]);
    let y3 = axis(&mut c, &[(3, Y)]);
    let zzzy = axis(&mut c, &[(0, Z), (1, Z), (2, Z), (3, Y)]);
    let xy = axis(&mut c, &[(1, X), (2, Y)]);
    c.rotate(z0, PauliAngle::new(1)).unwrap();
    c.rotate(y3, PauliAngle::new(-1)).unwrap();
    c.rotate(zzzy, PauliAngle::new(-1)).unwrap();
    c.rotate(xy, PauliAngle::new(1)).unwrap();
    let xs = box_xs(&c.to_svg().unwrap());
    // Drawing order is program order: Z0, Y3, ZZZY, XY.
    assert_eq!(xs[0], xs[1]);
    assert_eq!(xs[0], xs[3], "XY moves into the first column");
    assert!(xs[2] > xs[0], "ZZZY needs its own column");
    // In program order, XY comes after ZZZY.
    let program = c
        .to_svg_with(SvgOptions {
            reorder: false,
            ..SvgOptions::default()
        })
        .unwrap();
    let xs = box_xs(&program);
    assert!(xs[3] > xs[2]);
}

/// Anticommuting operations keep their order even when reordering.
#[test]
fn reorder_keeps_anticommuting_operations_in_order() {
    let mut c = PbcCircuit::new(2, 0);
    let z0 = c.z(0).unwrap();
    let x0 = c.x(0).unwrap();
    let z1 = c.z(1).unwrap();
    c.rotate(z0, PauliAngle::new(1)).unwrap();
    c.rotate(x0, PauliAngle::new(1)).unwrap();
    c.rotate(z1, PauliAngle::new(1)).unwrap();
    let xs = box_xs(&c.to_svg().unwrap());
    assert!(xs[1] > xs[0], "X0 after Z0");
    assert_eq!(
        xs[2], xs[0],
        "Z1 commutes with both and joins the first column"
    );
}
