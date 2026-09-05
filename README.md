# ⚡️ tzap

[![CI](https://github.com/qqq-wisc/tzap/actions/workflows/ci.yml/badge.svg)](https://github.com/qqq-wisc/tzap/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/tzap-opt.svg)](https://crates.io/crates/tzap-opt)
[![PyPI](https://img.shields.io/pypi/v/tzap.svg)](https://pypi.org/project/tzap/)
![Rust](https://img.shields.io/badge/Rust-000000?logo=rust&logoColor=white)
![Lean 4](https://img.shields.io/badge/Lean_4-black?logo=lean&logoColor=white)
![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue)
[![arXiv](https://img.shields.io/badge/arXiv-2605.13929-b31b1b.svg)](https://arxiv.org/abs/2605.13929)

[**Installation**](#installation) · [**Using tzap**](#running-tzap) &nbsp;**|**&nbsp;  [Qiskit integration](https://github.com/qqq-wisc/tzap/blob/main/docs/qiskit.md) · [PennyLane integration](https://github.com/qqq-wisc/tzap/blob/main/docs/pennylane.md)

A super fast, Rust-based optimizer for large Clifford+T circuits.
- tzap is state-of-the-art in *speed*, *scalability*, and *gate-count reduction*.
- tzap **minimizes T-count** with a new linear-time phase folding algorithm, based on [this paper](https://arxiv.org/abs/2605.13929).
- tzap implements a new and fast **superoptimization** pass.
- A **formally verified** Lean port of tzap is available in [`lean`](lean/).

tzap is **multiple orders of magnitude** faster than other optimizers&mdash;and **linearly** **scales** to **millions** of gates!
Here's a runtime comparison to two powerful optimizers on increasingly larger circuits.

<p align="center">
  <img src="https://raw.githubusercontent.com/qqq-wisc/tzap/main/assets/comparison.svg"
       alt="Runtime comparison of tzap, VOQC, and QuiZX on GF multipliers"
       style="width: 90%; height: auto;">
</p>

## Installation

You can use tzap as a command-line utility or a library.

### Install the binary

These options install the standalone native `tzap` executable.

**Homebrew** (macOS/Linux):

```bash
brew install qqq-wisc/tap/tzap
```

**Prebuilt release binary** (macOS/Linux):

```bash
curl -LsSf https://github.com/qqq-wisc/tzap/releases/latest/download/tzap-opt-installer.sh | sh
```

You can also build and install tzap from [crates.io](https://crates.io/crates/tzap-opt) (`cargo install tzap-opt`) or build from source (`cargo install --path .`).
You can also use tzap through the Rust API; see the [Rust API documentation](https://github.com/qqq-wisc/tzap/blob/main/API.md).

### Integrations with Qiskit and PennyLane
You can also use tzap as a Python library and apply it as a
Qiskit optimization pass or PennyLane transform. See the
[Qiskit API guide](https://github.com/qqq-wisc/tzap/blob/main/docs/qiskit.md) or
[PennyLane API guide](https://github.com/qqq-wisc/tzap/blob/main/docs/pennylane.md) for framework-specific setup.



## Running tzap

The standard command-line workflow is described below.

**Optimize a circuit**

```bash
tzap input.qasm -o output.qasm
```

For example, using a benchmark in this repo:

```console
$ tzap benchmarks/feynman/hwb12.qasm -o optimized.qasm
⚡️ tzap v0.6.0
  Parsed benchmarks/feynman/hwb12.qasm (5.5 MB) in 0.080s
	└─ 20 qubits · 514,412 gates
  Loaded superoptimizer table in 0.021s

  Converged after 6 rounds

  ┌─ Final result · 43.7% fewer gates · 1.595s ──────────────────────────┐
  │ Gates    ━━━━━━━━━━━━━╸────────────────── ↓43.7% · 514,412 → 289,484 │
  │ 2q gates ━━━━━╸────────────────────────── ↓18.7% · 191,803 → 155,914 │
  │ T/Tdg    ━━━━━━━━━━━━━━━╸──────────────── ↓49.9% · 171,465 →  85,897 │
  │ Depth    ━━━━━━━╸──────────────────────── ↓24.3% · 274,781 → 207,940 │
  └──────────────────────────────────────────────────────────────────────┘
  wrote optimized.qasm
```

**Optimization levels**

| Level | Description |
|---|---|
| `-O1` | phase folding + basic gate cancellation. Fastest; captures most of the T-gate reduction. |
| `-O2` | Adds superoptimization to `-O1`. |
| **`-O3`** | **Default.** Repeats `-O2` until reaching a fixpoint.  |
| `-Osuper` | Like `-O3`, but with more superoptimization power (slower on first use). |

```bash
tzap benchmarks/feynman/hwb12.qasm -O1 -o optimized.qasm
```

**Decompose Rz into Clifford+T**

Use `--decompose-rz` when the target backend only accepts Clifford+T; tzap uses [gridsynth](https://crates.io/crates/rsgridsynth). `--epsilon` trades approximation accuracy for circuit size (default `1e-10`; larger is coarser).

```bash
tzap input.qasm -o output.qasm --decompose-rz --epsilon 1e-6
```

Use `--decompose-cz` to decompose CZ gates into `H`+`CX`+`H` before the
optimization pipeline.

**Custom pipeline**

`--passes` runs an explicit, ordered sequence of passes in place of the default pipeline.

```bash
tzap input.qasm -o output.qasm --passes CancelGates,PhaseFoldRand
tzap input.qasm -o output.qasm --passes DecomposeCz,CancelGates,PhaseFoldRand
```

## Circuit support

tzap supports a subset of OpenQASM 2.0:

- **Gates:** `h`, `x`, `z`, `s`, `sdg`, `t`, `tdg`, `rz`, `cx`, `ccx`, `ccz`, `cz`, `measure`, `reset`
- **Declarations:** `qreg`, `creg`
- **Not supported:** classical conditionals (`if`), custom gate definitions (`gate`), barriers, `include` files (besides `qelib1.inc`, which is ignored)
- Unrecognized lines produce an error

Toffoli (`ccx`) and doubly controlled-Z (`ccz`) are auto-decomposed into Clifford+T. Controlled-Z (`cz`) is kept native so phase folding and cancellation can operate through it; use `--decompose-cz` for `H`+`CX` output. `Rz` is left as-is unless you pass `--decompose-rz`.

## Correctness

1. **Fuzzing and equivalence verification** on small random circuits and benchmark circuits.
2. **Lean port:** the core optimizer is implemented and proven sound in Lean 4 — see [`lean/`](lean/).

## Citation

If you use tzap in your research, please cite:

```bibtex
@misc{albarghouthi2026tzap,
      title={Linear-Time T-Gate Optimization via Random Abstraction}, 
      author={Aws Albarghouthi},
      year={2026},
      eprint={2605.13929},
      archivePrefix={arXiv},
      primaryClass={cs.PL},
      url={https://arxiv.org/abs/2605.13929}, 
}
```
