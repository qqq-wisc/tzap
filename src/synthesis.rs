//! Certified synthesis of accepted angles. The dependency proposes a word;
//! independent integer interval arithmetic decides whether it is acceptable.
mod exact_word;
pub(crate) mod interval;
use crate::angle::Angle;
use crate::circuit::Gate;
use dashu_float::{FBig, round::mode::HalfEven};
use dashu_int::IBig;
use interval::{Complex, Context, Interval, Matrix};
use num_bigint::BigInt;
use num_integer::Integer;
use std::{
    collections::VecDeque,
    sync::{Mutex, OnceLock},
};

const MAX_BITS: usize = 2048;
const MAX_WORD: usize = 8192;
const MAX_ATTEMPTS: usize = 3;
const MIN_EPSILON: f64 = 1e-80;
/// Temporary exact scaling; never round or store a half in the public Angle IR.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Scale {
    Whole,
    Half,
    NegativeHalf,
}
type Cache = VecDeque<(Angle, Scale, u64, Vec<char>)>;
static VERIFIED_WORDS: OnceLock<Mutex<Cache>> = OnceLock::new();

#[derive(Debug)]
pub(crate) enum Failure {
    Limit(&'static str),
    Verification(&'static str),
    Interval(&'static str),
    Unsupported(&'static str),
}
impl From<&'static str> for Failure {
    fn from(message: &'static str) -> Self {
        Self::Unsupported(message)
    }
}
impl std::fmt::Display for Failure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let message = match self {
            Self::Limit(m) | Self::Verification(m) | Self::Interval(m) | Self::Unsupported(m) => m,
        };
        f.write_str(message)
    }
}

struct Target {
    ctx: Context,
    pi: Interval,
    reduced: Interval,
    matrix: Matrix,
}
#[cfg(test)]
fn prepare(angle: &Angle, bits: usize) -> Result<Target, Failure> {
    prepare_scaled_with_period(angle, bits, 2, Scale::Whole)
}
fn prepare_with_period(angle: &Angle, bits: usize, period: i64) -> Result<Target, Failure> {
    prepare_scaled_with_period(angle, bits, period, Scale::Whole)
}
fn prepare_scaled_with_period(
    angle: &Angle,
    bits: usize,
    period: i64,
    scale: Scale,
) -> Result<Target, Failure> {
    let ctx = Context::new(bits);
    let pi = ctx.pi();
    let theta = crate::symbolic_angle::evaluate::enclose(angle, &ctx).map_err(Failure::Interval)?;
    let theta = match scale {
        Scale::Whole => theta,
        Scale::Half => theta.div_int(2),
        Scale::NegativeHalf => theta.div_int(2).neg(),
    };
    // Any integer quotient is valid. Enclose the subtraction with mathematical
    // pi; only require the resulting interval to be small enough for the kernel.
    let tau = pi.scale(&period.into());
    let turns = (theta.middle() + tau.middle() / 2u32).div_floor(&tau.middle());
    let reduced = theta.sub(&tau.scale(&turns));
    let (sin, cos) = ctx
        .sin_cos(&reduced.div_int(2))
        .map_err(Failure::Interval)?;
    let mut matrix = interval::identity(&ctx);
    matrix[0][0] = Complex {
        re: cos.clone(),
        im: sin.neg(),
    };
    matrix[1][1] = Complex { re: cos, im: sin };
    Ok(Target {
        ctx,
        pi,
        reduced,
        matrix,
    })
}

type TwoQubitMatrix = [[Complex; 4]; 4];
fn two_identity(ctx: &Context) -> TwoQubitMatrix {
    std::array::from_fn(|i| std::array::from_fn(|j| Complex::real(ctx, i64::from(i == j))))
}
/// Verify the complete controlled operator, including its inactive subspace.
pub(crate) fn verify_controlled(input: &Gate, gates: &[Gate], epsilon: f64) -> Result<(), Failure> {
    let (angle, control, target) = match input {
        Gate::cp {
            lambda,
            control,
            target,
        } => (lambda, *control, *target),
        Gate::crx {
            theta,
            control,
            target,
        }
        | Gate::cry {
            theta,
            control,
            target,
        }
        | Gate::crz {
            theta,
            control,
            target,
        } => (theta, *control, *target),
        _ => return Err("expected a controlled rotation".into()),
    };
    let exact = controlled_word(gates, control, target)?;
    let accuracy = (-epsilon.min(1.0).log2()).ceil() as usize;
    let mut bits = initial_bits(angle, accuracy);
    loop {
        match verify_controlled_at_precision(input, angle, &exact, epsilon, bits) {
            Ok(ControlledVerification::Verified) => return Ok(()),
            Ok(ControlledVerification::Rejected) => {
                return Err(Failure::Verification(
                    "complete controlled operator failed certification",
                ));
            }
            Ok(ControlledVerification::Inconclusive) if bits < MAX_BITS => {}
            Err(Failure::Interval(_)) if bits < MAX_BITS => {}
            Ok(ControlledVerification::Inconclusive) => {
                return Err(Failure::Verification(
                    "complete controlled operator failed certification at the precision limit",
                ));
            }
            Err(error) => return Err(error),
        }
        bits = (bits * 2).min(MAX_BITS);
    }
}

fn controlled_word(
    gates: &[Gate],
    control: u32,
    target: u32,
) -> Result<exact_word::ExactWord<4>, Failure> {
    let mut exact = exact_word::ExactWord::<4>::identity();
    for gate in gates {
        if let Gate::cnot {
            control: c,
            target: t,
        } = gate
        {
            let c = if *c == control {
                0
            } else if *c == target {
                1
            } else {
                return Err("invalid expanded control".into());
            };
            let t = if *t == control {
                0
            } else if *t == target {
                1
            } else {
                return Err("invalid expanded target".into());
            };
            exact.cnot(c, t);
            continue;
        }
        let (word, qubit): (&[char], u32) = match gate {
            Gate::h(q) => (&['H'], *q),
            Gate::x(q) => (&['X'], *q),
            Gate::z(q) => (&['S', 'S'], *q),
            Gate::s(q) => (&['S'], *q),
            Gate::sdg(q) => (&['S', 'S', 'S'], *q),
            Gate::t(q) => (&['T'], *q),
            Gate::tdg(q) => (&['T', 'T', 'T', 'T', 'T', 'T', 'T'], *q),
            _ => return Err("unsupported expanded synthesis instruction".into()),
        };
        let q = if qubit == control {
            0
        } else if qubit == target {
            1
        } else {
            return Err("invalid expanded qubit".into());
        };
        for instruction in word {
            exact.apply(*instruction, q).map_err(Failure::Interval)?;
        }
    }
    Ok(exact)
}

enum ControlledVerification {
    Verified,
    Rejected,
    Inconclusive,
}

fn verify_controlled_at_precision(
    input: &Gate,
    angle: &Angle,
    exact: &exact_word::ExactWord<4>,
    epsilon: f64,
    bits: usize,
) -> Result<ControlledVerification, Failure> {
    let prepared = prepare_with_period(angle, bits, 4)?;
    let ctx = &prepared.ctx;
    let mut expected = two_identity(ctx);
    // Local bit 0 is the control, local bit 1 the target, independently of
    // physical operand ordering. CR* must retain its literal 4*pi period.
    let cos = prepared.matrix[0][0].re.clone();
    let sin = prepared.matrix[1][1].im.clone();
    match input {
        Gate::cp { .. } => {
            let (sin, cos) = ctx
                .sin_cos(&prepared.reduced.div_int(2))
                .map_err(Failure::Interval)?;
            // Square exp(i theta/2); avoid a larger trig kernel argument.
            let half = Complex { re: cos, im: sin };
            expected[3][3] = half.mul(&half, ctx);
        }
        Gate::crz { .. } => {
            expected[1][1] = prepared.matrix[0][0].clone();
            expected[3][3] = prepared.matrix[1][1].clone();
        }
        Gate::crx { .. } => {
            expected[1][1].re = cos.clone();
            expected[3][3].re = cos;
            let off = Complex {
                re: ctx.integer(0),
                im: sin.neg(),
            };
            expected[1][3] = off.clone();
            expected[3][1] = off;
        }
        Gate::cry { .. } => {
            expected[1][1].re = cos.clone();
            expected[3][3].re = cos;
            expected[1][3] = Complex {
                re: sin.neg(),
                im: ctx.integer(0),
            };
            expected[3][1] = Complex {
                re: sin,
                im: ctx.integer(0),
            };
        }
        _ => unreachable!(),
    }
    let actual = exact.intervals(ctx);
    // Normalize by the phase of the inactive-control |00> amplitude. Its
    // exact magnitude is nonzero when the lower interval bound proves it.
    let anchor = &actual[0][0];
    let magnitude = ctx
        .mul(&anchor.re, &anchor.re)
        .add(&ctx.mul(&anchor.im, &anchor.im));
    let reciprocal = ctx.reciprocal_sqrt(&magnitude).map_err(Failure::Interval)?;
    let phase = Complex {
        re: ctx.mul(&anchor.re, &reciprocal),
        im: ctx.mul(&anchor.im.neg(), &reciprocal),
    };
    let mut squared_lower = BigInt::from(0);
    let mut squared_upper = BigInt::from(0);
    for i in 0..4 {
        for j in 0..4 {
            let (lower, upper) = actual[i][j]
                .mul(&phase, ctx)
                .difference_square_bounds(&expected[i][j]);
            squared_lower += lower;
            squared_upper += upper;
        }
    }
    let budget = ctx.float(epsilon);
    if squared_upper <= &budget.lo * &budget.lo {
        Ok(ControlledVerification::Verified)
    } else if squared_lower > &budget.hi * &budget.hi {
        // This phase-normalized Frobenius certificate cannot pass even with
        // tighter enclosures. In particular, a nonidentity rotation need not
        // exhaust the precision budget during the initial empty-word probe.
        Ok(ControlledVerification::Rejected)
    } else {
        Ok(ControlledVerification::Inconclusive)
    }
}
/// First fixed-point precision to try. Rationals enclose exactly at any
/// precision, so only an angle's magnitude (which scales pi's error during
/// period reduction) needs extra bits, not the length of its denominator.
/// Symbolic expressions fall back to their coefficient size. Callers widen on
/// an inconclusive enclosure, up to `MAX_BITS`.
fn initial_bits(angle: &Angle, accuracy: usize) -> usize {
    let magnitude = match angle.as_rational().or_else(|| angle.pi_coefficient()) {
        Some(r) => {
            let r = r.to_big_rational();
            r.numer().bits().saturating_sub(r.denom().bits()) as usize + 2
        }
        None => angle.coefficient_bits(),
    };
    magnitude
        .max(64)
        .saturating_add(accuracy + 128)
        .min(MAX_BITS)
}
fn dyadic(integer: &BigInt, bits: usize, precision: usize) -> FBig<HalfEven> {
    // Integer spelling is exact; never reinterpret an f64 decimal as the target.
    FBig::from_parts(
        IBig::from_str_radix(&integer.to_string(), 10).unwrap(),
        -(bits as isize),
    )
    .with_precision(precision)
    .value()
}

/// Returned words have independently verified projective spectral error <= eps.
pub(crate) fn synthesize(angle: &Angle, epsilon: f64) -> Result<Vec<char>, Failure> {
    synthesize_scaled(angle, Scale::Whole, epsilon)
}
pub(crate) fn synthesize_scaled(
    angle: &Angle,
    scale: Scale,
    epsilon: f64,
) -> Result<Vec<char>, Failure> {
    if epsilon < MIN_EPSILON {
        return Err(Failure::Limit(
            "certification precision limit: epsilon must be at least 1e-80",
        ));
    }
    let cache = VERIFIED_WORDS.get_or_init(|| Mutex::new(VecDeque::new()));
    if let Some((_, _, _, word)) = cache.lock().unwrap_or_else(|e| e.into_inner()).iter().find(
        |(prior, prior_scale, bound, _)| {
            prior.literal_eq(angle) && *prior_scale == scale && *bound == epsilon.to_bits()
        },
    ) {
        return Ok(word.clone());
    }
    let accuracy_bits = (-epsilon.min(1.0).log2()).ceil() as usize;
    let mut bits = initial_bits(angle, accuracy_bits);
    let target = loop {
        match prepare_scaled_with_period(angle, bits, 2, scale) {
            Ok(target)
                if target.reduced.hi.clone() - &target.reduced.lo
                    <= target.ctx.float(epsilon / 16.0).lo =>
            {
                break target;
            }
            _ if bits < MAX_BITS => bits = (bits * 2).min(MAX_BITS),
            _ => {
                return Err(Failure::Interval(
                    "symbolic target could not be enclosed within the precision budget",
                ));
            }
        }
    };
    let budget = target.ctx.float(epsilon.min(1.0));
    // Tiny rotations may be proved close to I without calling the synthesizer.
    if interval::verifies(&target.matrix, &[], &target.pi, &budget, &target.ctx)
        .map_err(Failure::Interval)?
    {
        return Ok(Vec::new());
    }
    // Initialize the dependency's lazy constants at ample precision. A caller
    // may already have initialized them; the independent verifier still protects
    // against stale constants and every dependency rounding error.
    rsgridsynth::common::set_prec_bits(MAX_BITS);
    let _ = &*rsgridsynth::common::PI;
    let _ = &*rsgridsynth::common::TAU;
    for attempt in 0..MAX_ATTEMPTS {
        let requested = epsilon.min(1.0) / (4u64 << (attempt * 2)) as f64;
        let mut config = rsgridsynth::config::config_from_theta_epsilon(
            0.0,
            requested,
            attempt as u64,
            false,
            true,
        );
        let precision = rsgridsynth::common::get_prec_bits().max(accuracy_bits + 128);
        rsgridsynth::common::set_prec_bits(precision);
        config.theta = dyadic(&target.reduced.middle(), bits, precision);
        let exact_budget = Context::new(1074).float(requested);
        config.epsilon = dyadic(&exact_budget.lo, 1074, precision);
        let result = rsgridsynth::gridsynth::gridsynth_gates(&mut config);
        let word: Vec<_> = result.gates.chars().take(MAX_WORD + 1).collect();
        if word.len() > MAX_WORD {
            return Err(Failure::Limit("synthesis output limit reached"));
        }
        if interval::verifies(&target.matrix, &word, &target.pi, &budget, &target.ctx)
            .map_err(Failure::Interval)?
        {
            // Reuse only a previously verified word at the identical literal
            // angle and budget. The cache is bounded independently of inputs.
            let mut cache = cache.lock().unwrap_or_else(|e| e.into_inner());
            if cache.len() == 64 {
                cache.pop_front();
            }
            cache.push_back((angle.clone(), scale, epsilon.to_bits(), word.clone()));
            return Ok(word);
        }
    }
    Err(Failure::Verification(
        "candidate verification failed after 3 attempts",
    ))
}

#[cfg(test)]
mod tests;
