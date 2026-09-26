# FTCircuitBench Clifford+T suite, large tier

The remaining 32 of FTCircuitBench's 95 inputs: those with more than 60,000
arbitrary rotations, which `../ftcircuitbench` skips. Same source, same
pipeline (FTCircuitBench's Gridsynth at precision 5), same caveats; see
[`../ftcircuitbench/README.md`](../ftcircuitbench/README.md).

| family | circuits | qubits | gates | T/Tdg |
| --- | --- | --- | --- | --- |
| `hamiltonians` | 21 | 9–128 | 51,068,156 | 15,820,480 |
| `hhl` | 1 | 21–21 | 2,949,323 | 1,133,342 |
| `qpe` | 5 | 10–12 | 101,414,705 | 38,995,149 |
| `qsvt` | 5 | 9–13 | 54,541,232 | 21,283,980 |
| **total** | 32 | | 209,973,416 | 77,232,951 |

Circuits range from 478k to 24.6M gates. The largest are QPE on H2 (12
qubits, about 24M gates and 9.2M T each) and QSVT at 13 qubits.

## How it was built

With the same `generate.py`, in its `--streaming` mode:

```bash
uv run --no-sync --project ~/git/FTCircuitBench python benchmarks/ftcircuitbench/generate.py \
    --ftcb ~/git/FTCircuitBench --streaming --min-rotations 60001 \
    --max-rotations 100000000 --jobs 2 --out benchmarks/ftcircuitbench-large
```

The normal mode holds each synthesized circuit as a Qiskit circuit, about
1.6 KB of memory per gate, which cannot fit these circuits. `--streaming`
makes the same calls into FTCircuitBench's pipeline (input preparation, the
Rz intermediate, and nwqec's C++ Gridsynth with epsilon = 10⁻⁵), then writes
nwqec's QASM output directly. FTCircuitBench's final basis step is a
pass-through for that output, and the generator checks every gate against
that basis.

**Validation:** on nine circuits spanning every family, the streaming and
normal modes gave identical T and `cx` counts. They differ only in
single-qubit Clifford gates, which FTCircuitBench's pipeline also varies from
run to run.

**Resources:** peak memory is about 8 GB per worker, on the input side
(Qiskit's conversion to the Rz intermediate), hence `--jobs 2` on a 24 GB
machine. Synthesis took about 90 minutes of CPU time, and the output is
1.7 GB, so this tier should be regenerated rather than committed.
