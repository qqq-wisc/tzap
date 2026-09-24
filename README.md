# ⚡️ tzap

[![CI](https://github.com/qqq-wisc/tzap/actions/workflows/ci.yml/badge.svg)](https://github.com/qqq-wisc/tzap/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/tzap-opt.svg)](https://crates.io/crates/tzap-opt)
[![PyPI](https://img.shields.io/pypi/v/tzap.svg)](https://pypi.org/project/tzap/)
![Rust](https://img.shields.io/badge/Rust-000000?logo=rust&logoColor=white)
![Lean 4](https://img.shields.io/badge/Lean_4-black?logo=lean&logoColor=white)
![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue)
[![arXiv](https://img.shields.io/badge/arXiv-2605.13929-b31b1b.svg)](https://arxiv.org/abs/2605.13929)

[**Installation**](#installation) · [**Using tzap**](#running-tzap) &nbsp;**|**&nbsp;  [Qiskit integration](https://github.com/qqq-wisc/tzap/blob/main/docs/qiskit.md) · [PennyLane integration](https://github.com/qqq-wisc/tzap/blob/main/docs/pennylane.md)

A super fast, Rust-based optimizer for large Clifford+T/Rz circuits.
- tzap is state-of-the-art in *speed*, *scalability*, and *gate-count reduction*.
- tzap **minimizes T-count** with a new linear-time phase folding algorithm, based on [this paper](https://arxiv.org/abs/2605.13929).
- tzap implements a new and fast **superoptimization** pass, based on [this paper](https://ia.cr/2026/2115).
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
⚡️ tzap v0.6.1
  Parsed benchmarks/feynman/hwb12.qasm (5.5 MB) in 0.079s
	├─ 20 qubits · 514,412 gates
	└─ Circuit gates: {h, x, t, tdg, cx}
  Optimizing input circuit
  Loaded MURM in 0.039s
	└─ Synthesis basis: {h, x, z, s, sdg, t, tdg, cx}

  Converged after 6 rounds

  ┌─ Final result · 44.6% fewer gates · 1.379s ────────────────┐
  │ Gates    ━━━━━━━━━╸──────────── ↓44.6% · 514,412 → 284,848 │
  │ 2q gates ━━━━╸───────────────── ↓22.1% · 191,803 → 149,500 │
  │ T/Tdg    ━━━━━━━━━━╸─────────── ↓49.9% · 171,465 →  85,889 │
  │ Depth    ━━━━━━╸─────────────── ↓28.4% · 274,781 → 196,865 │
  └────────────────────────────────────────────────────────────┘
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

**Optional decomposition**

CCX, CCZ, CZ, and Rz stay native by default. To decompose them, use:

- `--decompose-ccx` to decompose CCX and CCZ
- `--decompose-cz` to decompose CZ into CX+H
- `--decompose-rz` to decompose Rz via gridsynth

**PBC output**

See the [PBC exchange format](docs/pbc.md) for syntax and measurement examples.

```bash
tzap input.qasm --to-pbc -o output.pbc
```

Conversion runs last, after optimization and requested decompositions. Inputs
must have no resets and only terminal measurements; use `--decompose-rz` for Rz.
Export preserves all quantum and classical outputs, including post-measurement
states. The entire remaining Clifford suffix is retained as named gates,
whether none, some, or all qubits are measured.

```text
qubits 2
registers 1
r 1 1 X0 Z1
m -1 Z1 -> c0
h 0
```

`r <k> <sign> <factors>` rotates by `k*pi/8` using `exp(-i*k*pi/8*P)`.
`m <sign> <factors> -> cN` measures the signed Pauli product: +1 gives bit 0,
-1 gives bit 1. Signs are `1` or `-1`; omitted factors are identity (an empty
list is the identity). Register writes may overwrite earlier values. Named
Clifford gates follow the `r` and `m` instructions to restore quantum outputs.
`-o -` writes PBC to stdout. JSON metrics describe the optimized gate circuit
before conversion.
Export materializes Pauli strings with a default budget of 16 million sparse-work units;
unlike compressed conversion, expanded output is not guaranteed linear in size.

## Circuit support

tzap supports a subset of OpenQASM 2.0:

- **Gates:** `h`, `x`, `z`, `s`, `sdg`, `t`, `tdg`, `rz`, `cx`, `ccx`, `ccz`, `cz`, `measure`, `reset`
- **Declarations:** `qreg`, `creg`
- **Not supported:** classical conditionals (`if`), custom gate definitions (`gate`), barriers, `include` files (besides `qelib1.inc`, which is ignored)
- Unrecognized lines produce an error

## Correctness

1. **Fuzzing and equivalence verification** on small random circuits and benchmark circuits.
2. **Lean port:** the core optimizer is implemented and proven sound in Lean 4 — see [`lean/`](lean/).

## Citation

If you use tzap in your research, please cite the following papers:

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

```bibtex
@misc{cryptoeprint:2026/2115,
      author = {Aws Albarghouthi},
      title = {Fast Quantum-Circuit Superoptimization},
      howpublished = {Cryptology {ePrint} Archive, Paper 2026/2115},
      year = {2026},
      url = {https://eprint.iacr.org/2026/2115}
}
```
