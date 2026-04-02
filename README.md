# ⚡️ tzap

[![CI](https://github.com/qqq-wisc/tzap/actions/workflows/ci.yml/badge.svg)](https://github.com/qqq-wisc/tzap/actions/workflows/ci.yml)
![Rust](https://img.shields.io/badge/Rust-000000?logo=rust&logoColor=white)
![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue)

A fast, Rust-based T-gate optimizer for large Clifford+T circuits. tzap (pronounced *T-zap*) applies phase folding that is O(n) in the number of gates.

It takes OpenQASM 2.0 circuits as input, optimizes them, and outputs optimized OpenQASM 2.0. The supported gate set is: `h`, `x`, `z`, `s`, `sdg`, `t`, `tdg`, `rz`, `cx`, `ccx`, `cz`.

**Gate handling:**

- **Toffoli (`ccx`)** gates are automatically decomposed into Clifford+T before optimization.
- **Rz** gates are left as-is by default. Pass `--decompose-rz` to decompose them into Clifford+T via [gridsynth](https://crates.io/crates/rsgridsynth).

## Usage

CLI usage or API usage (see [API.md](API.md)).

```
tzap input.qasm                           # optimize, print stats only
tzap input.qasm -o output.qasm            # write optimized circuit to file
tzap input.qasm -o output.qasm --cancel   # enable gate cancellation pass
tzap input.qasm -o output.qasm --decompose-rz  # decompose Rz via gridsynth
```

Output is only written when `-o` is given.

### Example

```
$ tzap qasm/barenco_tof_5.qasm

⚡️ tzap
  Parsing qasm/barenco_tof_5.qasm (0.0 MB)
	└─ 9 qubits · 50 gates · 0 T/Tdg · 0.000s

  Toffoli decomposition
	└─ 218 gates · 84 T · 0.000s
  Phase folding
	└─ 146 gates · 40 T · 0.000s

  ⚡️ Result
	├─ Gates  218 → 146 (↓33.0%)
	├─ T/Tdg  84 → 40 (↓52.4%)
	└─ Time   0.000s
```

## Limitations

tzap supports a subset of OpenQASM 2.0:

- **Supported gates:** `h`, `x`, `z`, `s`, `sdg`, `t`, `tdg`, `rz`, `cx`, `ccx`, `cz`
- **Not supported:** classical registers (`creg`), measurement (`measure`), conditionals (`if`), custom gate definitions (`gate`), barriers, and `include` files (besides `qelib1.inc`, which is ignored)
- Unrecognized lines will produce an error

## Building

```
cargo build --release
```
