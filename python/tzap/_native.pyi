class TzapError(Exception): ...
class QasmError(TzapError): ...
class OptimizationError(TzapError): ...

__version__: str
__all__ = [
    "OptimizationError",
    "QasmError",
    "TzapError",
    "__version__",
    "_optimize_qasm",
]

def _optimize_qasm(
    qasm: str,
    *,
    level: str = "O3",
    passes: list[str] | None = None,
    fixpoint: bool = False,
    decompose_rz: bool = False,
    decompose_cz: bool = False,
    decompose_ccx: bool = False,
    rz_epsilon: float = ...,
    parallel: bool = False,
    superopt_qubits: int | None = None,
    superopt_window_gates: int | None = None,
    superopt_murm_entries: int | None = None,
    superopt_gates: str = "auto",
) -> tuple[
    str,
    tuple[
        tuple[int, int, int, int, int],
        tuple[int, int, int, int, int],
        tuple[int, int, int, int, int],
    ],
]: ...
