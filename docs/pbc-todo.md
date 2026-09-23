# PBC branch: next steps

Branch: `codex/clifford-to-pbc`

## 1. PBC visualizer

- Build a visualizer for `.pbc` files with publication-quality, Litinski-style
  circuit diagrams.
- Show multi-qubit Pauli rotations and measurements as joint operations, with
  clear angles, signs, and classical-register destinations.
- Distinguish Clifford rotations, non-Clifford rotations, and measurements.
- Support readable layouts for larger circuits and SVG/PDF export.

## 2. Large-circuit experiments

- Run conversion and export on large benchmark circuits, including CCX/CCZ-heavy
  circuits and circuits with terminal measurements.
- Measure optimization, compressed PBC conversion, and text export separately:
  runtime, peak memory, output size, and expansion-budget failures.
- Record rotation counts, Pauli-weight distributions, and measurement counts.
- Check scaling with gate count and qubit count independently; keep benchmark
  commands and configurations reproducible.

## 3. Optimization opportunities

- Merge rotations about equal Pauli axes and cancel inverse rotations.
- Explore commuting rotations to expose merges and reduce non-Clifford depth.
- Remove rotations that cannot affect the requested measurement results.
- For partial readout, explore discarding measured quantum outputs while
  preserving the unmeasured outputs; avoid unnecessary frame-restoring rotations.
- Investigate bounded-Pauli-weight conversion and its representation/runtime
  trade-offs.
- Validate each transformation with exact small-circuit semantics and fuzzing;
  measure its cost separately from the linear-time converter.
