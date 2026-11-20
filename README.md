# tzap, `pauli-eval` branch

This branch holds the code for the paper *Pauli Abstractions for Quantum Circuits*:

- **PauliAbs**, the Pauli-abstraction phase-folding pass (`PhaseFoldPauli` in tzap);
- our Rust port of Feynman's **state folding** (`StateFold1`, `StateFold2`, `StateFoldInf`);
- phase folding over **programs with control flow** (the `tzap-prog` binary);
- the **Lean development** that mechanizes the paper (`lean-paulifold/`).

For tzap in general, see the README on the `main` branch.

## Build

```bash
cargo build --release
```

This builds `target/release/tzap` (circuits) and `target/release/tzap-prog` (programs).

## Circuits (RQ1)

Each configuration runs one pass alone on a Clifford+T(+Rz) OpenQASM circuit:

```bash
tzap input.qasm --passes <PASS> --no-parallel -q --json -o output.qasm
```

| Paper name | `<PASS>` |
|---|---|
| PauliAbs | `PhaseFoldPauli` |
| StateFold(1) | `StateFold1` |
| StateFold(2) | `StateFold2` |
| StateFold(∞) | `StateFoldInf` |
| Linear | `PhaseFoldRand` |

`--json` prints the gate counts and the pass time. For example:

```bash
target/release/tzap benchmarks/feynman/hwb12.qasm --passes PhaseFoldPauli --no-parallel -q --json -o out.qasm
```

Benchmarks with Toffoli gates must be expanded into Clifford+T first:

```bash
tzap input.qasm --passes DecomposeToffoli -q -o expanded.qasm
```

### PauliAbs knobs

`PhaseFoldPauli` reads three environment variables, used for the ablations:

| Variable | Effect |
|---|---|
| `TZAP_PFP_BUDGET=<n>` | cap each fold attempt at `n` steps (default: unbounded, as in the paper) |
| `TZAP_PFP_INDEX=0` | turn off the fingerprint index (search all earlier axes) |
| `TZAP_PFP_SUPPORT=0` | turn off support lists for the commutation check |

## Programs with control flow (RQ2)

```bash
tzap-prog <CONFIG> program.qasm [--print]
```

It prints the qubit count, the static T count before and after, and the time as JSON; `--print`
also prints the optimized program.

| Paper name | `<CONFIG>` |
|---|---|
| PauliAbs (as reported) | `pauli-full` |
| PauliAbs, plain | `pauli` |
| StateFold(1), (2), (∞) | `statefold1`, `statefold2`, `statefold0` |

`pauli-full` keeps up to 100 abstract relations before joining and uses eigenstate facts after
resets and measurements; `pauli-disj` and `pauli-zero` turn on only one of the two, and
`pauli-full<k>` keeps up to `k` relations. The programs of Amy and Lunderville are in
`benchmarks/amy-lunderville/`, for example:

```bash
target/release/tzap-prog pauli-full benchmarks/amy-lunderville/loop-simple.qasm
```

## Feynman baseline

The paper also runs Feynman's own state folding, with inputs treated as unknown states:

```bash
feynopt -purecircuit -statefold <d> input.qasm   # d = 1, 2, or 0 for unbounded
```

## Lean

```bash
cd lean-paulifold
lake build
```

See `lean-paulifold/README.md` for how the Lean names map to the paper.
