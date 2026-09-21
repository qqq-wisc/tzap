# Native CZ, CCX, and CCZ Optimization Plan

This plan adds native `cz`, `ccx`, and `ccz` support throughout tzap while
keeping all decomposition options explicit. Lean changes are out of scope for
this branch.

## 1. Circuit gate-set detection and reporting

Add a gate-set profile derived from `Circuit::gates`. Derivation keeps the
profile correct even when Rust callers mutate the public gate vector directly.

After loading a circuit, print a stable, QASM-named set directly beneath the
existing parsed-circuit metrics:

```text
Parsed hamiltonian-simulation-n13.qasm (22.2 MB) in 0.255s
  ├─ 20 qubits · 1,955,569 gates
  └─ Circuit gates: {h, t, tdg, cx, ccx, ccz}
```

- Include every gate kind actually present, including `rz`, `measure`, and
  `reset`.
- Always list gates in this fixed canonical order rather than alphabetically
  or by first appearance:

  ```text
  h, x, z, s, sdg, t, tdg, rz, cx, cz, ccx, ccz, measure, reset
  ```

- Render the gate-set line in the normal foreground color, not the grey/dim
  style used for secondary metadata.
- Keep it grouped with the existing parsed-circuit message rather than
  printing a separate top-level status message.
- Print through the informational UI so `--quiet` suppresses it and QASM on
  stdout remains clean.
- Include the input gate set as an array in JSON reports.
- Recompute the profile after optimization and decomposition.
- Derive profiles for the whole circuit before parallel chunking, so all
  chunks make the same pass and MURM-selection decisions.

## 2. PhaseFoldRand support

Extend fingerprints to model nonlinear `ccx` behavior:

```text
X(q):           fp[q] += 1
CX(c, t):       fp[t] += fp[c]
CCX(a, b, t):   fp[t] += fp[a] * fp[b]
H(q):           fp[q] = fresh()
reset(q):       fp[q] = 0
CZ/CCZ:         no wire-state change
```

- Preserve native `cz` and `ccz` gates while allowing ordinary rotations to
  fold across them.
- Do not initially collect or re-emit the phase contributed by `ccz` itself.
- Track polynomial degree bounds.
- Conservatively revert to an opaque fresh target when nonlinear degree makes
  the finite-field collision guarantee ineffective.
- Add focused nonlinear-folding, diagonal-pass-through, and randomized
  equivalence tests.

## 3. CancelGates support

Extend the existing combined cancellation pass:

- Recognize `ccx` with swapped controls as the same gate.
- Retain symmetric matching for `cz` and fully symmetric matching for `ccz`.
- Add commuting lookahead cancellation for `ccx` and `ccz`.
- Add conservative commutation rules and explicit blocker tests.
- Enable native-gate lookahead logic according to the detected circuit gate
  set.
- Keep adjacent cancellation, Hadamard reduction, and lookahead cancellation
  in one combined fixpoint.

## 4. SuperOpt and MURM gate selection

Rename the SuperOpt synthesis table to **MURM**, short for **minimal unitary
representative map**. A MURM maps each exactly represented unitary to a minimal
representative circuit in the configured synthesis basis. SuperOpt remains the
name of the optimization pass; SuperOpt queries a MURM when replacing a
window.

The ordinary SuperOpt synthesis library remains available:

```text
h, x, z, s, sdg, t, tdg, cx
```

The only automatically optional generators are:

```text
cz, ccx, ccz
```

By default, each optional generator is added when it occurs in that
optimization stage's circuit.

Symmetric gates are enumerated canonically to avoid duplicate generators:

- Generate `cz(a, b)` only for `a < b`.
- Generate `ccx(a, b, target)` only for `a < b`, with the target distinct.
- Generate `ccz(a, b, c)` only for `a < b < c`.

This affects MURM construction only. Input circuits may use any valid operand
ordering.

Rename the corresponding implementation and reporting concepts, including:

- `UnitaryCircuitTable` to `Murm`.
- `SuperOptTableConfig` to `MurmConfig`.
- `shared_synthesis_table` and table-cache helpers to MURM-named equivalents.
- `src/super_opt/table.rs` to `src/super_opt/murm.rs`.
- Observer/JSON records that specifically describe loading or building the
  synthesis map.
- Cache directory, filenames, headers, diagnostics, documentation, and tests.
- Table-specific tuning names such as `table_entries` to `murm_entries`, while
  keeping `SuperOpt` terminology for the pass and its window configuration.

Also update:

- `LibraryGate` conversion, arity, inverse, disjointness, and serialization.
- The serialized gate record to support three operands.
- `MurmConfig` with the effective gate-library mask.
- In-memory MURM-cache keys, filenames, and cache headers.
- The cache format version.

The pre- and post-decomposition stages may intentionally use different cached
MURMs.

Whenever a MURM is loaded or built, report the exact effective
library beneath the existing status line:

```text
Loaded MURM in 0.021s
  └─ Synthesis basis: {h, x, z, s, sdg, t, tdg, cx, ccx}
```

Use the same canonical gate order as the circuit-gate display, and render the
basis line in the normal foreground color rather than grey/dim.

## 5. `--superopt-gates`

Add a SuperOpt synthesis-basis option with automatic, base-only, and explicit
forms:

```bash
tzap input.qasm --superopt-gates auto
tzap input.qasm --superopt-gates base
tzap input.qasm --superopt-gates h,x,z,s,sdg,t,tdg,cx,cz
```

Semantics:

- `auto` is the default. It selects the base SuperOpt library plus any of
  `{cz, ccx, ccz}` present in that optimization stage's circuit.
- `base` selects exactly `h,x,z,s,sdg,t,tdg,cx`.
- A comma-separated list selects exactly the listed synthesis basis. An
  explicitly listed optional gate may therefore be emitted even when it was
  absent from the stage's input circuit.
- The option affects only gates SuperOpt may emit in replacements.
- It does not change parsing or reject input gates.
- It does not constrain CancelGates, PhaseFoldRand, or decomposition output.
- Internally, represent the setting as:

  ```text
  SuperOptGates::Auto
  SuperOptGates::Base
  SuperOptGates::Explicit(GateSet)
  ```

Additional behavior:

- Use lowercase QASM gate names.
- Reject unsupported names and an empty list with clear diagnostics.
- Define and test the treatment of duplicate names.
- Permit the flag with optimization levels and explicit pipelines containing
  SuperOpt.
- Propagate it through Rust, Python, Qiskit, PennyLane, JSON, help text, and
  documentation.
- Include the effective synthesis basis in the MURM cache identity.
- A distinct explicit basis requires a distinct MURM: filtering a MURM built
  from a larger basis is not sufficient because that map may
  retain a shortest forbidden representative instead of a longer permitted
  one.

## 6. Optimize, decompose, optimize

All decomposition options are opt-in middle stages:

```text
native input
  -> optimize
  -> requested decompositions
  -> optimize again
  -> output
```

Flags:

- `--decompose-ccx`: lower both `ccx` and `ccz`.
- `--decompose-cz`: lower `cz`.
- `--decompose-rz`: synthesize `rz`.
- No `ccx`, `ccz`, `cz`, or `rz` decomposition runs by default.

Apply requested decompositions in this fixed order:

1. CCX/CCZ
2. CZ
3. Rz

Run the second optimization stage only when a requested decomposition actually
changed the circuit. The pre- and post-decomposition stages independently
compute their gate profiles and effective SuperOpt libraries.

For example, with `--superopt-gates auto --decompose-ccx --decompose-cz`:

```text
Input circuit gates:
  {h, t, rz, cx, cz, ccx}

Native-stage synthesis basis:
  {h, x, z, s, sdg, t, tdg, cx, cz, ccx}

Post-decomposition circuit gates:
  {h, t, tdg, rz, cx}

Post-decomposition synthesis basis:
  {h, x, z, s, sdg, t, tdg, cx}
```

The two optimization stages therefore load different MURMs when their
effective synthesis bases differ.

Explicit `--passes` pipelines remain exactly user-ordered. Decomposition flags
remain incompatible with `--passes`; users place the corresponding
decomposition passes directly into an explicit pipeline.

## 7. Reporting and APIs

Propagate `decompose_ccx = false`, `superopt_gates`, and the staged pipeline
through:

- Rust `Options`.
- CLI parsing and help.
- JSON options, gate-set fields, pass records, and MURM load/build records.
- Python native bindings and type stubs.
- Qiskit.
- PennyLane.

Progress output should distinguish:

```text
Optimizing input circuit
Decomposing CCX/CCZ
Decomposing CZ
Decomposing Rz
Optimizing decomposed circuit
```

JSON should retain the original input metrics and record both optimization
stages and each decomposition pass separately.

## 8. README and API documentation

Update the README tagline to:

> A super fast, Rust-based optimizer for large Clifford+T/Rz circuits.

Document:

- Native `ccx`, `ccz`, `cz`, and `rz` behavior.
- All decomposition flags being opt-in.
- The optimize, decompose, optimize ordering.
- `--decompose-ccx` covering both `ccx` and `ccz`.
- Automatic optional SuperOpt generators.
- `--superopt-gates` modes, syntax, scope, and examples.

Apply corresponding updates to Rust, Qiskit, PennyLane, and Python
documentation.

## 9. Exhaustive tests

Use reduced SuperOpt bounds for combinational tests so they remain practical.

### Normal-CI integration matrix

Do not run the literal 5,120-case Cartesian product in normal CI. Most of its
cases multiply orthogonal choices and add runtime without improving failure
localization. Instead, use these structured suites:

- **64 pipeline cases:** eight input subsets of `{cz, ccx, ccz}` times eight
  combinations of `{decompose_ccx, decompose_cz, decompose_rz}`.
- **80 SuperOpt-selection cases:** eight native-gate input subsets times ten
  basis modes: `auto`, `base`, and the eight explicit base-plus-subset
  combinations over `{cz, ccx, ccz}`.
- **32 optimization-level cases:** four levels times eight representative
  circuits covering all native-gate subsets.
- **24 parallel cases:** eight native-gate subsets times `auto`, `base`, and a
  representative explicit basis.
- **Approximately 20–40 targeted cross-feature regressions** for interactions
  that do not fit cleanly into the factored suites.

This yields approximately 220–250 meaningful end-to-end cases in normal CI.
Each decomposition fixture also contains representative base gates and an
`rz`, so every requested decomposition has something to act on.

Across these suites, check:

- Circuit equivalence.
- Expected gates removed by decomposition.
- Effective SuperOpt generator set.
- Explicit synthesis-basis compliance.
- Sequential/parallel agreement.
- Valid QASM round-trip.
- Correct input/output gate-set reporting.

### Extended full matrix

Keep the complete 5,120-case Cartesian product as an ignored or scheduled
extended test:

- Eight input subsets of `{cz, ccx, ccz}`.
- Ten SuperOpt-basis modes.
- Eight decomposition combinations.
- Four optimization levels.
- Sequential and parallel execution.

This suite is for nightly/manual validation rather than every commit.

### Full explicit-basis logic

SuperOpt has eleven candidate emitted gate kinds:

```text
h, x, z, s, sdg, t, tdg, cx, cz, ccx, ccz
```

Test all 2,048 basis masks at the pure configuration-selection level: 2,047
non-empty explicit bases and the empty basis as an error case.
Build actual MURMs for:

- All eight native-gate masks.
- Every single-gate exclusion.
- Representative mixed base-gate bases.
- Empty and invalid explicit bases as error cases.

### Pass-specific tests

Add focused tests for:

- Nonlinear CCX phase folding.
- CZ/CCZ transparency.
- CCX control symmetry.
- CCZ operand symmetry.
- Positive and negative commutation cases.
- Native SuperOpt replacements.
- Gate serialization and MURM-cache invalidation.
- Pre/post-decomposition cleanup.
- CLI, JSON, Python, Qiskit, and PennyLane propagation.

## 10. Fuzzing

Expand the circuit generators so every optional gate subset is exercised.

For each `{cz, ccx, ccz}` subset:

- Generate circuits containing every selected gate and no excluded optional
  gate.
- Exercise every decomposition combination.
- Randomize the complete explicit `--superopt-gates` basis.
- Exercise sequential and parallel modes.
- Check equivalence with the existing equivalence infrastructure.
- Check SuperOpt output-basis compliance and decomposition postconditions.
- Check CancelGates idempotence.
- Save and report reproducible seeds for failures.

Add dedicated fuzz campaigns for:

- PhaseFoldRand with nonlinear CCX dependencies.
- CancelGates commutation boundaries.
- SuperOpt native-gate MURMs.
- Full optimize, decompose, optimize pipelines.

Run a bounded, deterministic fuzz budget in normal CI. Put longer campaigns
that sample the entire cross-product of gate sets, synthesis bases,
decomposition combinations, optimization levels, and execution modes in the
extended/nightly suite.

## 11. Markdown test report

Create `docs/native-gate-test-report.md` during implementation. It will
contain:

- A table of every test added, with direct source links.
- The property each test verifies.
- The exhaustive matrix dimensions and result counts.
- Fuzz seeds, case counts, and commands.
- Cold/warm MURM-cache results for every native-gate mask.
- Example OpenQASM circuits for nonlinear CCX phase folding, CZ cancellation,
  CCX cancellation with swapped controls, CCZ cancellation, native SuperOpt
  replacement, explicit SuperOpt output bases, combined decomposition flags, and
  pre/post-decomposition optimization.
- Expected and actual output gate sets.
- Final Rust, Python, QCEC, CLI, and fuzzing results.

The implementation handoff will link both this report and the specific test
files.
