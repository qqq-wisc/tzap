"""Console entry point shipped in the Python wheel."""

from __future__ import annotations

import argparse
import sys

from ._core import TzapError, optimize_qasm
from ._native import __version__


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tzap", description="fast Clifford+T/Rz circuit optimizer"
    )
    parser.add_argument("input", help="input OpenQASM 2.0 file")
    parser.add_argument("output_positional", nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("-o", "--output", help="write optimized QASM to this file")
    levels = parser.add_mutually_exclusive_group()
    levels.add_argument("-O1", action="store_const", const="O1", dest="level")
    levels.add_argument("-O2", action="store_const", const="O2", dest="level")
    levels.add_argument("-O3", action="store_const", const="O3", dest="level")
    levels.add_argument("-Osuper", action="store_const", const="Osuper", dest="level")
    parser.set_defaults(level="O3")
    parser.add_argument("--passes", help="comma-separated explicit pass pipeline")
    parser.add_argument("--fixpoint", action="store_true")
    parser.add_argument("--decompose-rz", action="store_true")
    parser.add_argument("--decompose-cz", action="store_true")
    parser.add_argument("--decompose-ccx", action="store_true")
    parser.add_argument("--superopt-gates", default="auto")
    parser.add_argument("--epsilon", type=float, default=1e-10)
    parser.add_argument("--parallel", action="store_true")
    parser.add_argument(
        "-v", "--version", action="version", version="tzap " + __version__
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    output = args.output or args.output_positional
    try:
        with open(args.input, "r", encoding="utf-8") as source:
            qasm = source.read()
        passes = None
        if args.passes is not None:
            passes = [name.strip() for name in args.passes.split(",") if name.strip()]
        result = optimize_qasm(
            qasm,
            level=args.level,
            passes=passes,
            fixpoint=args.fixpoint,
            decompose_rz=args.decompose_rz,
            decompose_cz=args.decompose_cz,
            decompose_ccx=args.decompose_ccx,
            superopt_gates=args.superopt_gates,
            rz_epsilon=args.epsilon,
            parallel=args.parallel,
        )
        if output:
            with open(output, "w", encoding="utf-8") as destination:
                destination.write(result.qasm)
    except (OSError, TzapError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    report = result.report
    print(
        f"{report.baseline.gates} -> {report.output.gates} gates; {report.baseline.t} -> {report.output.t} T/Tdg",
        file=sys.stderr,
    )
    if output:
        print(f"wrote {output}", file=sys.stderr)
    return 0
