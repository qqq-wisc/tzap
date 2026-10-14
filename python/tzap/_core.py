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
class NumericalReport:
    """No added rounding in folding; input interpretation and synthesis differ."""

    preserved_expressions: int
    numerical_fallbacks: int
    skipped_rounded_folds: int
    skipped_nonfinite_folds: int
    skipped_coefficient_limit_folds: int
    uncertified_syntheses: int
    randomized_matching: bool
    certified_syntheses: int = 0
    failed_certifications: int = 0
    certified_error_bound: float | None = 0.0
    input_policy: str = "exact_symbolic"
    folding: str = "no_added_rounding"
    synthesis_error_scope: str = (
        "accepted_angle_per_input_gate; summed_bound_covers_synthesis_only"
    )


@dataclass(frozen=True)
class OptimizationReport:
    """Input/baseline metrics and metrics after the staged optimization."""

    input: Metrics
    baseline: Metrics
    output: Metrics
    numerical: NumericalReport


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
    decompose_rotations: bool = False,
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

    Native ``ccx``, ``ccz``, ``cz``, phase gates, and rotations are preserved unless
    their matching ``decompose_*`` option is enabled. ``decompose_rotations``
    lowers P, RZ/RX/RY, CP, CRZ/CRX/CRY, and CH; ``decompose_rz`` is a compatibility alias.
    ``decompose_ccx`` also lowers CSWAP. ``superopt_gates`` is
    ``"auto"``, ``"base"``, or an explicit comma-separated synthesis basis.

    The returned OpenQASM uses only original ``qelib1.inc`` names: P and CP are
    written as ``u1``/``cu1``, and SX, SWAP, CSWAP, CRX, and CRY as exact
    expansions (CRX/CRY around ``crz``). Report metrics count those statements.

    The CPU-heavy optimizer releases Python's GIL while it runs.
    """

    pass_list: Sequence[str] | None
    pass_list = None if passes is None else tuple(passes)
    optimized, raw_report = _native._optimize_qasm(
        qasm,
        level=level,
        passes=pass_list,
        fixpoint=fixpoint,
        decompose_rotations=decompose_rotations,
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
        numerical=NumericalReport(*raw_report[3]),
    )
    return OptimizationResult(optimized, report)
