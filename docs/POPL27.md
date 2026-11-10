# POPL 2027 commands

The paper's main benchmark suites are under `benchmarks/`:

- `benchmarks/feynman/` contains the Feynman-suite circuits accepted by the current parser. The paper's runtime plots use the circuits with at least 2,000 gates. The frozen POPL artifact also retains `cycle_17_3` and `mod_adder_1048576`; they were removed from the current tree because they contain self-CNOT instructions.
- `benchmarks/cobble-t/` contains the synthesized Clifford+T Cobble circuits used for the reported comparisons.
- `benchmarks/cobble-rz/` contains the same Cobble programs before Rz synthesis; these are source circuits, not the inputs used for the cross-tool Clifford+T comparison.

`benchmarks/qft/` is an additional corpus and was not used in the POPL paper.

The RQ2 Cobble and Tower scaling instances contain up to hundreds of millions of gates and are too large to include in the repository; they were generated for the experiment. The synthetic RQ3 HWB-like family is generated as follows; the paper uses 16 qubits and five gates per stage:

Generate an HWB-like scaling circuit (16 qubits, 40,000 stages, 200,000 gates):

```sh
cargo run --release --bin gen_hwb_like -- 16 40000 input.qasm
```

Full `tzap` pipeline used in the paper:

```sh
tzap input.qasm --no-parallel --passes CancelGates,PhaseFoldRand -o output.qasm
```

Randomized phase folding only (`TZAP-PF`):

```sh
tzap input.qasm --no-parallel --passes PhaseFoldRand -o output.qasm
```

Sparse exact phase folding (`TZAP-SYMB`):

```sh
tzap input.qasm --no-parallel --passes PhaseFoldGlobalExpr -o output.qasm
```

Dense exact phase folding:

```sh
tzap input.qasm --no-parallel --passes PhaseFoldBitVecExact -o output.qasm
```
