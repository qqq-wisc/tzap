#!/usr/bin/env uv run
# /// script
# requires-python = ">=3.12"
# dependencies = ["pyzx"]
# ///
"""Compare tzap and PyZX T-gate optimization on QASM files."""

import argparse
import signal
import subprocess
import time
from pathlib import Path

import pyzx

QASM_DIR = Path(__file__).resolve().parent.parent / "qasm"
TZAP_BIN = Path(__file__).resolve().parent.parent / "target" / "release" / "tzap"


def count_gates(path: Path) -> int:
    return sum(1 for line in path.read_text().splitlines() if line.strip().endswith(";"))


def run_tzap(path: Path, timeout: float) -> tuple[int, float]:
    start = time.perf_counter()
    try:
        result = subprocess.run(
            [str(TZAP_BIN), str(path)],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return -2, time.perf_counter() - start
    elapsed = time.perf_counter() - start
    if result.returncode != 0:
        return -1, elapsed
    t_count = 0
    for line in result.stdout.splitlines():
        line = line.strip()
        if line.startswith("t ") or line.startswith("tdg "):
            t_count += 1
    return t_count, elapsed


class TimeoutError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise TimeoutError()


def run_pyzx(path: Path, timeout: float) -> tuple[int, float]:
    old = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(int(timeout) + 1)
    start = time.perf_counter()
    try:
        c = pyzx.Circuit.from_qasm_file(str(path))
        c_opt = pyzx.optimize.full_reduce(c)
        elapsed = time.perf_counter() - start
        return c_opt.tcount(), elapsed
    except TimeoutError:
        return -2, time.perf_counter() - start
    except Exception:
        return -1, time.perf_counter() - start
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def fmt_result(t: int) -> str:
    if t == -2:
        return "T/O"
    if t < 0:
        return "ERR"
    return str(t)


def main():
    parser = argparse.ArgumentParser(description="Compare tzap and PyZX")
    parser.add_argument("--max-gates", type=int, default=1000, help="max gate count (default: 1000)")
    parser.add_argument("--timeout", type=float, default=60, help="per-circuit timeout in seconds (default: 60)")
    args = parser.parse_args()

    if not TZAP_BIN.exists():
        print("Build tzap first: cargo build --release")
        raise SystemExit(1)

    files = sorted(QASM_DIR.glob("*.qasm"))
    files = [(f, count_gates(f)) for f in files]
    files = [(f, g) for f, g in files if g <= args.max_gates]
    files.sort(key=lambda x: x[1])

    if not files:
        print(f"No QASM files with <= {args.max_gates} gates")
        raise SystemExit(1)

    name_w = max(len(f.stem) for f, _ in files)
    print(f"{'Circuit':<{name_w}}  {'Gates':>5}  {'tzap T':>6}  {'PyZX T':>6}  {'tzap':>7}  {'PyZX':>7}")
    print("─" * (name_w + 42))

    for path, gates in files:
        tzap_t, tzap_time = run_tzap(path, args.timeout)
        pyzx_t, pyzx_time = run_pyzx(path, args.timeout)

        print(
            f"{path.stem:<{name_w}}  {gates:>5}  {fmt_result(tzap_t):>6}  {fmt_result(pyzx_t):>6}"
            f"  {tzap_time:>6.3f}s  {pyzx_time:>6.3f}s"
        )


if __name__ == "__main__":
    main()
