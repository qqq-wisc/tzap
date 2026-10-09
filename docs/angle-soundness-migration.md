# Angle soundness implementation and migration

This branch implements the accepted-angle policy in
[the design](angle-soundness-design.md), starting from `main` at `755c935`.
The release version remains unchanged; this is a Rust source compatibility
change to schedule for the next API release.

## Inputs and public API

`Gate::rz` takes an immutable checked `Angle`. Replace numeric construction
`Gate::rz(theta, q)` with `Gate::rz_f64(theta, q)?`, or construct an angle with
`Angle::from_f64(theta)?`. Use `Angle::pi_fraction(n, d)?` for an intentional
rational multiple of mathematical pi. Angle and coefficient fields are private;
NaN, infinity, zero denominators, and overflowing coefficients cannot enter
through public constructors.

`PiFraction` uses gcd normalization and checked widened arithmetic.
`Angle::components` exposes read-only normalized rational and finite-float
components. Equality and hashing compare structural rotation values modulo
exact `2*pi`, with signed zero canonicalized. They never compare with a
numerical tolerance. A floating `PI/4` remains numeric, independently of the
bits matching Rust's approximation to pi.

The guarantee begins at the accepted angle. Pi-free source expressions retain
binary64 evaluation order. Scalar expressions building pi coefficients use
bounded rational evaluation. Thus `0.1+0.2` has its usual rounded binary64
value, whereas `0.1*pi` means exactly `pi/10`.

Supported non-affine pi expressions, such as `pi*pi/4`, use an explicitly
recorded numeric input fallback. Powers, calculator functions, formal parameters,
and custom gate bodies remain outside the supported parser subset. Unknown
syntax is rejected. Expressions can span lines; lexical comment handling and
source line diagnostics are preserved. The parser bounds expression complexity
(128 recursive levels and 256 syntax nodes) to avoid unbounded recursive AST
construction, evaluation, or destruction.

Coefficient-limit and uncertifiable affine residual expressions retain their
validated source spelling and line. They are not folded, classified, or silently
approximated. The CLI prints an aggregate warning; the Python wrapper emits a
`RuntimeWarning`. Quiet CLI mode suppresses the human warning while JSON retains
the count. Invalid divisions and nonfinite numerical intermediates remain errors.

## Arithmetic and pass behavior

A checked TwoSum operation accepts a single-float merge only if its exact
correction is zero and every required intermediate is finite. Independent
test-only big-rational decoding of binary64 values verifies accepted results,
including subnormals, cancellation, large values, and seeded random pairs.
Production angle arithmetic has no big-integer or general symbolic dependency.

Failed phase merges leave the earlier group and current gate intact.
`PhaseFoldRand` starts a fresh group; `PhaseFoldPauli` commits deletion and frame
changes only after the merge and emitted cost succeed; `CnotMin` abandons the
whole chunk and copies its saved input if any phase arithmetic fails. Numeric
residuals are never wrapped using floating pi. Emission keeps pi contributions
and residuals separate and counts the resulting operations before committing
folds. Named Clifford+T phases use a modulo-eight path without per-T allocation.

`CancelGates`, Pauli folding, the SuperOpt pass adapter, and PBC conversion lower
certified quarter-turn Rz exactly. Other fractions and numeric Rz remain outside
SuperOpt's exact matrix domain and PBC's discrete fragment. The raw
`SuperOpt::run` analyzer retains input instruction indices: callers can explicitly
call `angle::lower_exact_rotations` before analysis. Its `Pass` adapter performs
that preprocessing automatically. MURM gate encodings and cache format stay the
same.

Toffoli and CZ decomposition copy unrelated angles unchanged. Exact channel
reference tests now cover certified quarter-turn Rz, measurement, and reset,
including arbitrary quantum inputs represented by Choi matrices. Floating
matrix comparisons remain supplemental checks.

## Output and approximation

QASM serialization prints exact coefficients using integer arithmetic and pi.
A mixed angle serializes as two adjacent commuting Rz gates. Finite numeric
values retain binary64 round-trip spelling with real-token lexical adjustments.
Preserved expressions retain source spelling.

`qasm::serialize_numeric_lossy` is an explicit approximate export for numeric-only
consumers. It lowers certified quarter-turns without approximation, counts
uncertified general fraction conversions, and rejects preserved expressions.
`Angle::to_f64_lossy` similarly makes a numeric conversion explicit.

Existing explicit synthesis requests count as consent to approximation; no new
opt-in flag is required. `DecomposeRz::try_run` validates epsilon, lowers exact
quarter-turns first, and returns `SynthesisError` for unavailable conversions.
The optimization driver propagates that error. The infallible `Pass` adapter
returns its original circuit on failure; use `try_run` or the driver when synthesis
is required. Calls remain serialized around gridsynth's global precision state.

Gridsynth's current `config_from_theta_epsilon` converts a floating target through
its decimal spelling and computes its working precision from epsilon. This
implementation does not certify that decimal target conversion, mathematical-pi
conversion, or the synthesis result against the accepted angle. Reports therefore
state that total conversion and whole-circuit error are uncertified. Epsilon is a
per-rotation parameter for the numeric synthesis target. The human CLI summary
also states this scope whenever synthesis occurs; quiet mode retains only the
structured report. Bounds greater than one use one internally to avoid unsigned
underflow in the dependency's precision calculation. Certified interval
conversion, verified synthesis error, and whole-circuit budgeting remain the
separate feature specified by the design.

The numeric Qiskit and PennyLane input adapters continue to supply finite numeric
angles. Pi fractions are never inferred from their floats. Mixed folds emit named
quarter-turn gates plus numeric residuals, so those adapters do not introduce
implicit symbolic-to-numeric conversion. Existing symbolic/trainable-input
rejections remain in place.

## Diagnostics

Rust `Report::numerical`, Python `OptimizationReport.numerical`, and CLI JSON
`.numerical` expose:

- preserved input expressions and numerical input fallbacks;
- skipped folds due to rounding, nonfinite arithmetic, or coefficient limits;
- whether the configured pipeline uses unchecked randomized matching;
- the number of uncertified synthesis calls.

Counters are per optimization run and propagated to parallel workers, rather than
shared globally between unrelated runs. Their values can differ when chunk
boundaries change optimization opportunities. The input policy is
`accepted_angles`; numerical folding adds no rounding. Randomized Boolean-function
matching in `PhaseFoldRand` still has its separate collision probability. A forced
collision regression documents this distinction. Pauli folding confirms sketch
hits with exact axes and commutation checks. An explicit deterministic pipeline
can omit `PhaseFoldRand`.

## Storage and validation

`Angle` occupies 16 bytes and `Gate` occupies 24 bytes on the tested 64-bit host,
compared with 16-byte gates on the baseline. Finite numeric values are inline;
fractions, mixed values, and preserved expressions use shared immutable payloads.
There is no per-numeric-Rz allocation in the angle representation and no per-T
allocation in phase accumulation. The larger gate stride affects scratch buffers
and circuit clones; performance measurements below report its practical effect.

Regression tests cover the audit's tiny-angle, near-quarter, cancellation,
huge-float, and overflow examples. They also check coefficient bounds, structural
hashes, symbolic round trips, exact PBC eligibility, synthesis errors, signed
rollback, parallel diagnostics, and measurement/reset channels. Tests that intend
mathematical pi now supply explicit pi provenance; separate numeric tests check
that floating pi is not classified.

### Local performance comparison

On the development macOS host, alternating three release CLI runs per version
with `-O1 -q --json` gave the following medians. Time includes process startup,
parsing, optimization, metrics, and JSON output; peak RSS comes from
`/usr/bin/time -l`. The baseline snapshot has the same production sources as
`main` at `755c935`. These are indicative measurements on one host, especially
for the small Rz case where startup dominates.

| Workload | Main time | New time | Main RSS | New RSS | Main/new output gates |
|---|---:|---:|---:|---:|---:|
| Feynman `hwb12` | 117.2 ms | 122.1 ms | 86.81 MiB | 105.41 MiB | 346,646 / 346,646 |
| Cobble Rz `laplacian-filter` | 4.2 ms | 4.3 ms | 2.78 MiB | 2.91 MiB | 1,263 / 1,290 |
| QFT `qft_q040_d68921` | 164.0 ms | 167.4 ms | 161.00 MiB | 191.59 MiB | 657,332 / 657,332 |

T counts remain 85,611, 172, and 310,907 respectively. The larger workloads
show a 2–4% time increase and a 19–21% peak-memory increase. The 24-byte gate
layout retains inline numeric payloads in exchange for that memory cost.
The sampled Rz run reports 35 skipped rounded folds; the other two workloads
report no skipped arithmetic. No workload reports coefficient-limit or nonfinite
folds, input fallback, preserved expressions, or synthesis.
The Rz workload's higher output count reflects conservative numerical folding
and separate emission of exact pi contributions and numeric residuals; recovering
those gates by adding numerical tolerances would violate the accepted contract.

### Validation performed

The full `cargo test --offline` suite passes, including 875 library tests,
12 angle regression tests, CLI integration suites, and API doctests.
`PBC_FUZZ_CASES=4 cargo test --offline --release --lib -- --ignored` passes all
23 extended tests, including benchmark rewrite audits, the 5,120-case native
configuration matrix, and exact PBC channel fuzz samples. The final signed
coefficient boundary change additionally passes the 12 angle regressions and
all 121 parser tests.

All 334 Python tests pass locally. Ruff lint and formatting, Rust formatting,
native Python stub checking, and diff whitespace checks pass. These are local
results; the cross-platform CI matrix has not run on this branch.
