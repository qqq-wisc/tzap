"""Build the FTCircuitBench Clifford+T suite.

Synthesizes FTCircuitBench's input circuits (https://github.com/pnnl/FTCircuitBench,
arXiv:2601.03185) into Clifford+T with its own Gridsynth pipeline, and writes
one OpenQASM 2.0 file per circuit plus a manifest.

Run inside FTCircuitBench's environment. One input (qpe_Hubbard_8q) is
OpenQASM 3, which needs the optional `qiskit_qasm3_import` package; `--no-sync`
keeps uv from removing it again:

    cd ~/git/FTCircuitBench && uv sync && uv pip install qiskit_qasm3_import
    uv run --no-sync --project ~/git/FTCircuitBench python benchmarks/ftcircuitbench/generate.py \
        --ftcb ~/git/FTCircuitBench [--precision 5] [--max-rotations 60000] [--jobs 8]

For circuits too large to hold in Qiskit after synthesis, add --streaming
(see ../ftcircuitbench-large/README.md).

Circuits whose inputs have more than --max-rotations arbitrary rotations
(rz/u1/u3/cu1/rx/ry lines) are skipped, which bounds the synthesized size.
"""

import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROTATION_PREFIXES = ("rz", "u3", "u1", "cu1", "rx", "ry")


def rotations(path: Path) -> int:
    with open(path) as handle:
        return sum(1 for line in handle if line.lstrip().startswith(ROTATION_PREFIXES))


def synthesize(source: str, target: str, precision: int) -> dict:
    """Synthesize one circuit; runs in a worker process."""
    from ftcircuitbench.parser import load_qasm_circuit
    from ftcircuitbench.transpilers import transpile_to_gridsynth_clifford_t
    from qiskit.qasm2 import dump

    start = time.time()
    circuit = load_qasm_circuit(source)
    clifford_t = transpile_to_gridsynth_clifford_t(
        circuit, gridsynth_precision=precision, return_intermediate=False, prefer_cpp=True
    )
    elapsed = time.time() - start
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(target, "w") as handle:
        dump(clifford_t, handle)
    counts = clifford_t.count_ops()
    return {
        "qubits": clifford_t.num_qubits,
        "gates": sum(v for k, v in counts.items() if k not in ("barrier", "measure")),
        "t": counts.get("t", 0) + counts.get("tdg", 0),
        "cx": counts.get("cx", 0),
        "synthesis_s": round(elapsed, 2),
    }


def synthesize_streaming(source: str, target: str, precision: int) -> dict:
    """The same pipeline as `synthesize`, for circuits too large to hold as a
    Qiskit circuit after synthesis (about 1.6 KB of memory per gate).

    It makes FTCircuitBench's nwqec path's calls directly: `prepare_input`,
    `to_intermediate_rz`, then nwqec's C++ Gridsynth with epsilon =
    10^-precision and keep_ccx=False. nwqec's output is already in the
    Clifford+T basis, so FTCircuitBench's final `enforce_pbc_basis` is a
    pass-through; here nwqec's QASM text is checked against that basis and
    written directly instead of being re-parsed into Qiskit.
    """
    import tempfile

    import nwqec
    from ftcircuitbench.transpilers._basis import (
        PBC_COMPATIBLE_CLIFFORD_T_BASIS,
        is_clifford_t_basis,
        prepare_input,
        to_intermediate_rz,
    )
    from ftcircuitbench.transpilers.nwqec_ct import _strip_non_semantic
    from qiskit.qasm2 import dumps

    start = time.time()
    circuit = prepare_input(source, is_file=True, remove_final_measurements=True)
    if is_clifford_t_basis(circuit):
        text = dumps(circuit)
    else:
        intermediate = _strip_non_semantic(to_intermediate_rz(circuit))
        del circuit
        with tempfile.NamedTemporaryFile("w", suffix=".qasm", delete=False) as tmp:
            tmp.write(dumps(intermediate))
        del intermediate
        try:
            synthesized = nwqec.to_clifford_t(
                nwqec.load_qasm(tmp.name), keep_ccx=False, epsilon=10.0 ** (-precision)
            )
        finally:
            os.remove(tmp.name)
        text = synthesized.to_qasm()
        del synthesized
    elapsed = time.time() - start

    allowed = set(PBC_COMPATIBLE_CLIFFORD_T_BASIS) | {"measure", "barrier"}
    counts: dict = {}
    qubits = 0
    for line in text.splitlines():
        word = line.split(" ", 1)[0].split("(", 1)[0]
        if word in ("OPENQASM", "include", "creg", "") or word.startswith("//"):
            continue
        if word == "qreg":
            qubits += int(line[line.index("[") + 1 : line.index("]")])
            continue
        if word not in allowed:
            raise ValueError(f"unexpected gate {word!r} in synthesized output")
        counts[word] = counts.get(word, 0) + 1
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(target, "w") as handle:
        handle.write(text)
    return {
        "qubits": qubits,
        "gates": sum(v for k, v in counts.items() if k not in ("barrier", "measure")),
        "t": counts.get("t", 0) + counts.get("tdg", 0),
        "cx": counts.get("cx", 0),
        "synthesis_s": round(elapsed, 2),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ftcb", type=Path, required=True, help="FTCircuitBench checkout")
    parser.add_argument("--precision", type=int, default=5, help="Gridsynth precision (digits)")
    parser.add_argument("--min-rotations", type=int, default=0)
    parser.add_argument("--max-rotations", type=int, default=60000)
    parser.add_argument("--only", nargs="*", help="restrict to these input stems")
    parser.add_argument("--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--streaming", action="store_true",
                        help="write nwqec's output directly (for circuits too large for Qiskit)")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()

    commit = subprocess.run(
        ["git", "-C", str(args.ftcb), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    inputs = sorted((args.ftcb / "qasm").glob("*/*.qasm"))
    selected = [(p, rotations(p)) for p in inputs]
    selected = [(p, r) for p, r in selected if args.min_rotations <= r <= args.max_rotations]
    if args.only:
        selected = [(p, r) for p, r in selected if p.stem in args.only]
    print(f"{len(selected)} of {len(inputs)} circuits selected "
          f"({args.min_rotations}-{args.max_rotations} rotations), FTCircuitBench {commit}",
          file=sys.stderr)

    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {}
        for path, count in selected:
            family = path.parent.name
            target = args.out / family / path.name
            run = synthesize_streaming if args.streaming else synthesize
            futures[pool.submit(run, str(path), str(target), args.precision)] = (
                family, path, count, target)
        for future in as_completed(futures):
            family, path, count, target = futures[future]
            try:
                result = future.result()
            except Exception as error:  # report and continue with the rest
                print(f"FAILED {path}: {error}", file=sys.stderr)
                continue
            rows.append({
                "family": family,
                "circuit": path.stem,
                "file": str(target.relative_to(args.out)),
                "source": str(path.relative_to(args.ftcb)),
                "input_rotations": count,
                **result,
            })
            print(f"{family}/{path.stem}: {result['gates']:,} gates, {result['t']:,} T "
                  f"({result['synthesis_s']}s)", file=sys.stderr)

    rows.sort(key=lambda r: (r["family"], r["input_rotations"]))
    manifest = args.out / "manifest.csv"
    with open(manifest, "w", newline="") as handle:
        handle.write(f"# FTCircuitBench {commit}, gridsynth precision {args.precision}, "
                     f"{args.min_rotations}-{args.max_rotations} input rotations\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} circuits and {manifest}", file=sys.stderr)


if __name__ == "__main__":
    main()
