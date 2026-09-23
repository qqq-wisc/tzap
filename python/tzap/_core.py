"""Python interface to tzap's native optimizer."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from . import _native

TzapError = _native.TzapError
QasmError = _native.QasmError
OptimizationError = _native.OptimizationError


@dataclass(frozen=True)
class Metrics:
    """Circuit metrics collected in one traversal by tzap."""

    gates: int
    two_qubit: int
    depth: int
    t: int
    rz: int


@dataclass(frozen=True)
class OptimizationReport:
    """Input/baseline metrics and metrics after the staged optimization."""

    input: Metrics
    baseline: Metrics
    output: Metrics


@dataclass(frozen=True)
class OptimizationResult:
    """Optimized OpenQASM and the metrics for that optimization run."""

    qasm: str
    report: OptimizationReport


def _metrics(raw: tuple[int, int, int, int, int]) -> Metrics:
    return Metrics(*raw)


def optimize_qasm(
    qasm: str,
    *,
    level: str = "O3",
    passes: Iterable[str] | None = None,
    fixpoint: bool = False,
    decompose_rz: bool = False,
    decompose_cz: bool = False,
    decompose_ccx: bool = False,
    rz_epsilon: float = 1e-10,
    parallel: bool = False,
    superopt_qubits: int | None = None,
    superopt_window_gates: int | None = None,
    superopt_murm_entries: int | None = None,
    superopt_gates: str = "auto",
) -> OptimizationResult:
    """Optimize an OpenQASM 2 program with tzap.

    ``level`` accepts ``"O1"``, ``"O2"``, ``"O3"`` (the default), or
    ``"Osuper"``. Supplying ``passes`` replaces the level's default pipeline.
    Pass names are the same as the CLI's ``--passes`` names.

    Native ``ccx``, ``ccz``, ``cz``, and ``rz`` gates are preserved unless
    their matching ``decompose_*`` option is enabled. ``superopt_gates`` is
    ``"auto"``, ``"base"``, or an explicit comma-separated synthesis basis.

    The CPU-heavy optimizer releases Python's GIL while it runs.
    """

    pass_list: Sequence[str] | None
    pass_list = None if passes is None else tuple(passes)
    optimized, raw_report = _native._optimize_qasm(
        qasm,
        level=level,
        passes=pass_list,
        fixpoint=fixpoint,
        decompose_rz=decompose_rz,
        decompose_cz=decompose_cz,
        decompose_ccx=decompose_ccx,
        rz_epsilon=rz_epsilon,
        parallel=parallel,
        superopt_qubits=superopt_qubits,
        superopt_window_gates=superopt_window_gates,
        superopt_murm_entries=superopt_murm_entries,
        superopt_gates=superopt_gates,
    )
    report = OptimizationReport(
        input=_metrics(raw_report[0]),
        baseline=_metrics(raw_report[1]),
        output=_metrics(raw_report[2]),
    )
    return OptimizationResult(optimized, report)
