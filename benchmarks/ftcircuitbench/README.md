# FTCircuitBench Clifford+T suite

Clifford+T circuits synthesized from the inputs of
[FTCircuitBench](https://github.com/pnnl/FTCircuitBench)
([arXiv:2601.03185](https://arxiv.org/abs/2601.03185)), a benchmark suite for
fault-tolerant compilation from PNNL, using FTCircuitBench's own Gridsynth
pipeline at precision 5.

| family | circuits | qubits | gates | T/Tdg |
| --- | --- | --- | --- | --- |
| `adder` | 4 | 4–64 | 1,577 | 624 |
| `hamiltonians` | 12 | 9–100 | 5,738,548 | 2,070,200 |
| `hamiltonians_5trotter` | 36 | 9–200 | 18,242,741 | 6,747,130 |
| `hhl` | 3 | 4–12 | 333,991 | 128,332 |
| `qft` | 4 | 4–63 | 404,370 | 156,614 |
| `qpe` | 1 | 8–8 | 1,212,538 | 466,229 |
| `qsvt` | 3 | 6–8 | 1,301,314 | 506,632 |
| **total** | 63 | | 27,235,079 | 10,075,761 |

Families: ripple-carry adders; Hamiltonian simulation (Ising, Heisenberg and
Fermi-Hubbard models on 1D, 2D and triangular lattices), both full-length
(`hamiltonians`) and five Trotter steps (`hamiltonians_5trotter`); HHL; QFT;
QPE on a Hubbard model; and QSVT on a banded circulant matrix.

## How it was built

`generate.py` (in this directory) reproduces the suite from an FTCircuitBench
checkout. For each input under FTCircuitBench's `qasm/`, it:

1. **Selects by size:** inputs with more than 60,000 arbitrary rotations
   (`rz`/`u1`/`u3`/`cu1`/`rx`/`ry`) are skipped, because synthesis turns each
   rotation into about 25 T gates. This keeps 63 of FTCircuitBench's 95
   inputs.
2. **Synthesizes:** calls `transpile_to_gridsynth_clifford_t` with
   `gridsynth_precision=5`, the same call FTCircuitBench's `gs` pipeline
   makes, and writes the result as OpenQASM 2.0.
3. **Records:** writes `manifest.csv` with each circuit's source file,
   input rotation count, qubits, gates, T count, `cx` count and synthesis
   time, headed by the FTCircuitBench commit (116df00).

Synthesis of all 63 circuits takes about 90 seconds on 8 cores. As a check,
`qft/qft_4q.qasm` has 459 T, matching the value FTCircuitBench's README
gives for this pipeline and precision.

Differences from FTCircuitBench's raw inputs, both inherited from its
pipeline: final measurements are removed, and `reset` gates are dropped. The
HHL and QPE inputs reset ancilla qubits mid-circuit (1–4 resets each), so
their Clifford+T versions here are unitary circuits that no longer match the
source algorithm at those points. They are the same circuits FTCircuitBench
itself benchmarks. The skipped inputs are the larger full-length Hamiltonian circuits,
QPE on H2 and the larger Hubbard QPE, HHL at 21 qubits, and QSVT at 9–13
qubits.

**Regenerating is not bit-for-bit reproducible.** FTCircuitBench's Gridsynth
pipeline (nwqec) is nondeterministic. Two runs on the same input give the same
T and `cx` counts but differ slightly in single-qubit Clifford gates, e.g.
48,370 vs 48,226 gates for `qft_18q` with 18,769 T both times. Expect about 1%
variation in total gate count; T and `cx` counts are stable.

To rebuild, see the instructions at the top of `generate.py`.
