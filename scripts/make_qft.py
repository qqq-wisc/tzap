# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "qiskit>=1.0",
#     "qiskit-aer",
# ]
# ///

import sys
from math import log2, pi
import qiskit.synthesis.qft
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import UnitarySynthesis, BasisTranslator
from qiskit.circuit.equivalence_library import SessionEquivalenceLibrary

NAM = ["cx", "h", "rz", "x"]
CLIFFORD_T = ["cx", "h", "t", "tdg", "x"]


def precompute_cp_decomps(unique_ks, gridsynth=False):
    """Decompose CP(2π/2^k) for each k into Clifford+T using a tiny 2-qubit circuit.

    Returns a dict mapping k -> list of (gate_name, [local_qubit_indices])
    where local qubit 0 = control, 1 = target.
    """
    decomps = {}
    sel = SessionEquivalenceLibrary
    gridsynth_pm = PassManager(
        passes=[
            BasisTranslator(target_basis=NAM, equivalence_library=sel),
            UnitarySynthesis(
                method="gridsynth",
                plugin_config={"epsilon": 1e-12},
                basis_gates=CLIFFORD_T,
                synth_gates=["rz"],
            ),
        ]
    )
    for k in sorted(unique_ks):
        qc = QuantumCircuit(2)
        qc.cp(2 * pi / (2**k), 0, 1)
        if gridsynth:
            tc = gridsynth_pm.run(qc)
        else:
            tc = transpile(qc, basis_gates=CLIFFORD_T, optimization_level=0)
        gates = []
        for inst in tc.data:
            local_qubits = [tc.find_bit(q).index for q in inst.qubits]
            gates.append((inst.operation.name, local_qubits))
        decomps[k] = gates
        print(f"  CP(2π/2^{k}): {len(gates)} gates", file=sys.stderr)
    return decomps


def make_qft_streaming(n, output_path, gridsynth=False):
    threshold = int(log2(n) + 2)
    approx_degree = (n - 1) - threshold
    print(
        f"n={n}, threshold={threshold}, approx_degree={approx_degree}", file=sys.stderr
    )

    print("Building logical QFT circuit...", file=sys.stderr)
    logical = qiskit.synthesis.qft.synth_qft_full(
        num_qubits=n, do_swaps=False, approximation_degree=approx_degree
    )

    unique_ks = set()
    for inst in logical.data:
        if inst.operation.name == "cp":
            angle = float(inst.operation.params[0])
            k = int(round(log2(2 * pi / angle)))
            unique_ks.add(k)

    print(f"Unique CP angles (k values): {sorted(unique_ks)}", file=sys.stderr)
    print(
        "Pre-computing Clifford+T decompositions (2-qubit circuits)...", file=sys.stderr
    )
    cp_decomps = precompute_cp_decomps(unique_ks, gridsynth=gridsynth)

    print(f"Writing QASM to {output_path} ...", file=sys.stderr)
    total = 0
    with open(output_path, "w") as f:
        f.write("OPENQASM 2.0;\n")
        f.write('include "qelib1.inc";\n')
        f.write(f"qreg q[{n}];\n")

        for inst in logical.data:
            name = inst.operation.name
            qubits = [logical.find_bit(q).index for q in inst.qubits]

            if name == "h":
                f.write(f"h q[{qubits[0]}];\n")
                total += 1
            elif name == "cp":
                angle = float(inst.operation.params[0])
                k = int(round(log2(2 * pi / angle)))
                control, target = qubits[0], qubits[1]
                local_to_global = {0: control, 1: target}
                for gname, local_qs in cp_decomps[k]:
                    global_qs = ", ".join(f"q[{local_to_global[q]}]" for q in local_qs)
                    f.write(f"{gname} {global_qs};\n")
                    total += 1
            else:
                qs = ", ".join(f"q[{q}]" for q in qubits)
                f.write(f"{name} {qs};\n")
                total += 1

    print(f"Done. Total gates written: {total}", file=sys.stderr)


if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Generate a QFT circuit in QASM")
    parser.add_argument("n", type=int, nargs="?", default=3, help="Number of qubits")
    parser.add_argument(
        "--gridsynth",
        action="store_true",
        help="Use gridsynth for RZ decomposition (default: use Qiskit transpiler)",
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file path (default: qft_<n>_<method>.qasm in current directory)",
    )
    args = parser.parse_args()

    suffix = "gridsynth" if args.gridsynth else "cliffordt"
    output_path = args.output or os.path.join(os.getcwd(), f"qft_{args.n}_{suffix}.qasm")
    make_qft_streaming(args.n, output_path, gridsynth=args.gridsynth)
