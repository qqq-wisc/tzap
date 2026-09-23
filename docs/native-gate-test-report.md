# Native CZ/CCX/CCZ Test Report

This report covers the native controlled-gate work on
`codex/native-controlled-gates`. Lean was intentionally left unchanged.

## Coverage summary

| Suite | Cases | CI policy | What it covers |
|---|---:|---|---|
| Staged pipeline matrix | 64 | Normal | 8 input subsets × 8 decomposition subsets |
| SuperOpt basis selection | 80 | Normal | 8 input subsets × `auto`, `base`, and 8 explicit native subsets |
| Optimization levels | 32 | Normal | 4 levels × 8 input subsets |
| Parallel execution | 24 | Normal | 8 input subsets × 3 basis modes |
| Deterministic cross-feature fuzz | 512 | Normal | Eight random samples for every input/decomposition subset, cycling basis and execution modes |
| Explicit basis parser | 2,048 | Normal | All 2,047 non-empty masks plus the empty-basis error |
| Native MURM cache round trips | 8 | Normal | Cold build and warm read for every optional-native mask |
| Full Cartesian matrix | 5,120 | Ignored/manual | 8 inputs × 10 bases × 8 decompositions × 4 levels × 2 execution modes |

The normal-CI matrix contains 712 end-to-end pipeline/fuzz cases before the
focused pass, CLI, JSON, cache, and Python regressions are counted.

## Tests added

### Circuit profiles and reporting

Source: [`src/circuit.rs`](../src/circuit.rs),
[`tests/cli.rs`](../tests/cli.rs), and
[`tests/cli_json.rs`](../tests/cli_json.rs).

| Test | Property |
|---|---|
| `gate_set_uses_the_canonical_qasm_order` | A derived `GateSet` includes all 14 supported gate kinds in fixed display order. |
| `superopt_gate_sets_are_exact` | Base, optional, and supported SuperOpt masks contain exactly the documented gates. |
| `metadata_tracks_direct_public_gate_mutation` | Gate profiles and legacy presence queries are derived from the actual gate vector, so direct Rust API mutation cannot leave stale optimizer metadata. |
| `parsed_circuit_gate_set_uses_canonical_order` | The parse banner prints `Circuit gates: {...}` with canonical ordering. |
| `murm_reports_auto_base_and_explicit_synthesis_bases` | MURM status reports the exact effective basis for all three CLI modes. |
| `stages_and_superopt_gate_mode_are_reported` | JSON records the native/decomposition/post-decomposition stages and canonicalizes duplicate explicit gates. |

### PhaseFoldRand

Source: [`src/phase_fold_rand.rs`](../src/phase_fold_rand.rs).

| Test | Property |
|---|---|
| `nonlinear_ccx_fingerprint_recovers_after_inverse_pair` | The GF(2^128) nonlinear target fingerprint changes across CCX and returns after its inverse, allowing phase cancellation without removing the native CCX gates. |
| `gf128_multiplication_obeys_field_identities` | The carry-less field implementation satisfies zero, identity, commutativity, distributivity, and the configured reduction polynomial. |
| `nonlinear_shared_control_dependencies_remain_equivalent` | CCX controls that share algebraic variables after CNOT propagation remain sound through phase folding. |
| `extreme_polynomial_degree_conservatively_stops_phase_folding` | A synthetic Fibonacci-degree CCX circuit crosses the `2^32` safety limit and proves that PhaseFoldRand forgets the target instead of making a probabilistically weak phase merge. |
| Existing CZ/CCZ transparency tests | Rotations fold across every diagonal operand while CZ/CCZ gate order and multiplicity are preserved. |
| 512-case deterministic fuzz suite | Random nonlinear CCX/CCZ dependencies remain equivalent through the full staged driver. |

### CancelGates

Source: [`src/cancel.rs`](../src/cancel.rs).

| Test | Property |
|---|---|
| `ccx_swapped_controls_cancel_across_commuting_control_phase` | CCX controls are symmetric and a control phase commutes through the cancellation. |
| `ccx_lookahead_is_blocked_by_a_target_phase` | A phase on the CCX target is a sound blocker. |
| `ccz_lookahead_cancels_across_diagonal_gates` | Fully symmetric CCZ pairs cancel through diagonal gates. |
| `every_positive_ccx_and_ccz_commutation_rule_is_unitary_sound` | Every gate shape accepted by the new CCX/CCZ positive commutation predicates is checked by exact unitary comparison. |
| Existing operand/blocker tests | All CCX/CCZ operands, measurement/reset boundaries, disjoint gates, and idempotence remain covered. |

### SuperOpt and MURMs

Source: [`src/super_opt/tests.rs`](../src/super_opt/tests.rs).

| Test | Property |
|---|---|
| `native_library_generators_are_canonical_and_complete` | CZ uses `a < b`, CCX uses ordered controls with a distinct target, and CCZ uses `a < b < c`; no symmetric duplicates are generated. |
| `native_library_gates_round_trip_the_four_byte_encoding` | All native generators survive the new tag-plus-three-operands cache encoding. |
| `compact_keys_canonicalize_symmetric_native_operands` | Reordered CZ, CCX-control, and CCZ operands share one normalized matrix-cache key. |
| `native_murms_synthesize_native_representatives` | Singleton CZ, CCX, and CCZ bases synthesize their native representative. |
| `explicit_cz_murm_rewrites_h_cx_h_to_native_cz` | An explicit CZ basis may emit CZ even when the input stage did not contain it. |
| `every_native_mask_builds_and_round_trips_its_own_murm` | All 8 base-plus-native masks build, persist, reload, and retain identical width/depth data. |
| `actual_murms_cover_every_single_gate_exclusion_and_mixed_basis` | Actual MURMs build for every one-gate exclusion and representative sparse mixed bases. |
| `disk_read_rejects_a_mismatched_config` | A cache created with one basis cannot be loaded for another basis. |
| `disk_read_rejects_structurally_corrupt_bodies_before_using_their_lengths` | Corrupt width/table counts, forward parents, and out-of-basis or out-of-width gates are rejected before unsafe allocation or traversal. |
| `disk_read_rejects_checksum_mismatches_and_trailing_bytes` | Structurally plausible body changes and appended data are rejected by the versioned full-body checksum and exact EOF check. |
| Renamed MURM cache/build tests | Process sharing, deterministic parallel builds, invalid config errors, saturation, header/body corruption, structural bounds, checksums, trailing bytes, and crate-version invalidation use MURM terminology. |

### Driver matrices and fuzzing

Source: [`src/optimize.rs`](../src/optimize.rs).

| Test | Cases | Property |
|---|---:|---|
| `ccx_decomposition_is_opt_in` | 2 | Native default and explicit CCX/CCZ lowering. |
| `auto_basis_is_recomputed_after_decomposition` | 1 | The observed native MURM is base+CCX and the observed post-decomposition MURM is base-only. |
| `requested_decomposition_constrains_an_explicit_post_basis` | 3 | Exact singleton CZ/CCX/CCZ bases cannot synthesize a requested decomposition back into the output. |
| `absent_requested_gate_is_still_forbidden_from_the_post_basis` | 1 | A requested gate family stays excluded even when another decomposition is what triggers the second optimizer stage. |
| `explicit_pipeline_recomputes_auto_at_every_superopt_boundary` | 6 | Three ordered-pipeline shapes, sequential and parallel, resolve `auto` from the circuit at each SuperOpt occurrence. |
| `every_explicit_superopt_basis_mask_round_trips` | 2,048 | Every supported explicit mask round-trips; empty and unsupported masks fail; duplicates deduplicate. |
| `auto_superopt_basis_covers_all_native_gate_subsets` | 8 | `auto` adds exactly the optional gates in the whole-stage profile; `base` never does. |
| `all_native_and_decomposition_subsets_follow_the_staged_contract` | 64 | Equivalence, QASM round trip, and decomposition postconditions for every input/decomposition subset. |
| `superopt_selection_matrix_covers_80_end_to_end_cases` | 80 | Effective MURM basis, equivalence, QASM validity, and no unrequested optional emissions. |
| `optimization_level_matrix_covers_32_end_to_end_cases` | 32 | O1/O2/O3/Osuper preserve every optional-native input subset. |
| `parallel_superopt_matrix_covers_24_end_to_end_cases` | 24 | Sequential and parallel results are mutually equivalent for auto/base/explicit modes. |
| `parallel_fixpoint_reports_the_maximum_chunk_round` | 1 | Parallel chunk telemetry collapses to one stage record with the maximum round count and convergence only when every chunk converged. |
| `bounded_native_pipeline_fuzz_covers_512_cross_feature_cases` | 512 | Eight deterministic random circuits cover every native/decomposition subset, cycling basis modes and sequential/parallel execution, plus QASM round trip, CancelGates idempotence, and every requested CCX/CCZ/CZ/Rz output postcondition. |
| `extended_full_native_configuration_matrix_has_5120_cases` | 5,120 | Ignored exhaustive Cartesian product with equivalence and every requested decomposition postcondition. |

The bounded fuzz seed is
`0x4e41544956450000 | (sample << 12) | (native_mask << 8) | decomposition_mask`;
every failure
prints the complete seed. Run it with:

```bash
cargo test --lib bounded_native_pipeline_fuzz_covers_512_cross_feature_cases -- --nocapture
```

Run the extended matrix manually with:

```bash
cargo test --release --lib extended_full_native_configuration_matrix_has_5120_cases -- --ignored --nocapture
```

### CLI, JSON, Python, Qiskit, and PennyLane

Sources: [`tests/cli.rs`](../tests/cli.rs),
[`tests/cli_errors.rs`](../tests/cli_errors.rs),
[`tests/cli_json.rs`](../tests/cli_json.rs), and
[`tests/python`](../tests/python).

| Test group | Property |
|---|---|
| CLI native/decomposition regressions | Native CCX/CCZ/CZ remain by default; `--decompose-ccx`, `--decompose-cz`, and `--decompose-rz` lower only when requested. |
| `explicit_superopt_basis_cannot_reintroduce_requested_decompositions` | Reproduces the formerly failing exact CZ/CCX/CCZ basis combinations through the CLI and verifies the written QASM. |
| `superopt_gates_rejects_missing_empty_and_unsupported_values` | Missing, empty, Rz, and unknown explicit bases have clear errors. |
| `passes_conflicts_with_decomposition_flags` | All three decomposition flags remain incompatible with explicit `--passes`. |
| JSON schema/stage/MURM tests | Input/output gate sets and effective MURM bases are machine-readable; pass, MURM, and every fixpoint record carry their owning stage index; parallel stages report the maximum chunk round; deprecated `table`, `fixpoint`, and cache-info `tables` views remain compatible. |
| `concurrent_cold_writers_leave_one_valid_warm_murm` | Two processes may build the same cold MURM concurrently without sharing a temporary file, leaking it, or corrupting the warm result. |
| `a_structurally_valid_cache_with_body_corruption_is_rebuilt` | A valid-header cache with a damaged body is rejected, rebuilt deterministically, and usable by the next warm process. |
| Python binding tests | Native defaults, all decomposition flags, all basis modes, invalid bases, and original-input baseline semantics reach Rust. |
| Qiskit tests | Default native preservation and opt-in CCX/CCZ lowering round-trip through Qiskit. |
| PennyLane tests | Default native preservation and opt-in Toffoli/CCZ lowering round-trip through PennyLane. |
| Python QCEC matrix | 8 deterministic random programs × 10 optimizer configs × 2 adapters = 160 external equivalence checks per Python version. |

## Example circuits

Each snippet is a gate body; prepend the standard OpenQASM 2 header and the
shown three-qubit register where needed.

### Nonlinear CCX phase folding

```qasm
qreg q[3];
t q[2];
ccx q[0],q[1],q[2];
ccx q[0],q[1],q[2];
tdg q[2];
```

`PhaseFoldRand` removes the cancelling `t`/`tdg` pair and preserves both CCX
gates. Input gate set: `{t, tdg, ccx}`. Output gate set: `{ccx}`.

### CZ cancellation through a diagonal

```qasm
qreg q[2];
cz q[0],q[1];
t q[0];
cz q[1],q[0];
```

`CancelGates` recognizes reversed CZ operands and emits only `t q[0]`.

### CCX cancellation with swapped controls

```qasm
qreg q[3];
ccx q[0],q[1],q[2];
s q[0];
ccx q[1],q[0],q[2];
```

The CCX pair cancels through the control phase; the output gate set is `{s}`.

### CCZ cancellation through diagonal gates

```qasm
qreg q[3];
ccz q[0],q[1],q[2];
rz(pi/7) q[0];
t q[2];
ccz q[2],q[0],q[1];
```

The fully symmetric CCZ pair cancels; the output gate set is `{t, rz}` in
canonical order.

### Native SuperOpt replacement and explicit output basis

```qasm
qreg q[2];
h q[1];
cx q[0],q[1];
h q[1];
```

With `--superopt-gates cz`, the one-gate explicit MURM representative is
`cz q[0],q[1]`. Input gate set: `{h, cx}`. Output gate set: `{cz}`. With
`base`, SuperOpt cannot introduce CZ.

### Combined decomposition flags

```qasm
qreg q[3];
ccx q[0],q[1],q[2];
ccz q[0],q[1],q[2];
cz q[0],q[1];
rz(pi/7) q[2];
```

Running with `--decompose-ccx --decompose-cz --decompose-rz
--superopt-gates base` uses the fixed CCX/CCZ → CZ → Rz order. The final gate
set is a subset of `{h, x, z, s, sdg, t, tdg, cx}`.

### Native and post-decomposition cleanup

```qasm
qreg q[3];
h q[0];
h q[0];
ccx q[0],q[1],q[2];
rz(pi/7) q[2];
```

The input stage first removes `h; h`. Requested decomposition then lowers
CCX and/or Rz, and the post-decomposition stage cleans the generated sequence.
JSON records both optimization stages separately.

## Verification results

The following are required before handoff; the final branch was checked with
the commands shown here:

| Command | Result |
|---|---|
| `cargo test --lib --quiet` | 635 passed, 18 ignored |
| `cargo test --release --lib extended_full_native_configuration_matrix_has_5120_cases -- --ignored --nocapture` | All 5,120 configurations passed |
| `cargo test --test cli --test cli_cache --test cli_errors --test cli_json --test cli_streams --quiet` | 207 passed |
| `cargo check --features python` | Passed |
| `cargo clippy --all-targets --features python -- -D warnings` | Passed |
| `cargo test --doc --quiet` | 2 passed, 10 ignored |
| `uv run ruff check python tests/python pyproject.toml` | Passed |
| `uv run ruff format --check python tests/python` | Passed |
| `uv run --no-sync pytest tests/python -m 'not qcec' -q` | 169 passed |
| `uv run --no-sync pytest tests/python/test_qcec.py -m qcec -q` | 161 passed: 160 adapter equivalence proofs plus the negative checker control |
