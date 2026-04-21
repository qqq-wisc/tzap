#!/usr/bin/env python3
"""
Compare T-count reduction and runtime across tzap, QuiZX, PyZX, and Feynman.

Usage:
    python3 scripts/compare.py [qasm_dir]   (default: qasm/)

Only runs on files up to --max-kb (default 2000 KB).
"""

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

# ── locate tzap ──────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT   = os.path.dirname(SCRIPT_DIR)

TZAP = os.path.join(REPO_ROOT, "target", "release", "tzap")
if not os.path.isfile(TZAP):
    TZAP = os.path.join(REPO_ROOT, "target", "debug", "tzap")

# ── runner: tzap ─────────────────────────────────────────────────────────────

# Matches tzap's Result lines: "T/Tdg  21 → 15 (↓28.6%)" and "Gates  146 → 99 (↓32.2%)"
_TZAP_RESULT_T_RE = re.compile(r"T/Tdg\s+([\d,]+)\s+→\s+([\d,]+)")
_TZAP_RESULT_G_RE = re.compile(r"Gates\s+([\d,]+)\s+→\s+([\d,]+)")
# Matches tzap's parsing line: "└─ 7 qubits · 34 gates · 0 T/Tdg · 0.000s"
_TZAP_PARSE_RE = re.compile(r"qubits\s*·\s*([\d,]+)\s*gates\s*·\s*[\d,]+\s*T/Tdg\s*·\s*([\d.]+)s")

def run_tzap(qasm_path: str, timeout: float = 60.0) -> dict:
    t0 = time.perf_counter()
    try:
        r = subprocess.run(
            [TZAP, qasm_path, "-o", "/dev/null", "--cancel"],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"error": f"timeout >{timeout:g}s"}
    elapsed = time.perf_counter() - t0
    if r.returncode != 0:
        return {"error": r.stderr.strip()}
    t_in = t_after = None
    g_in = g_after = None
    parse_s = None
    for line in r.stderr.splitlines():
        m = _TZAP_RESULT_T_RE.search(line)
        if m:
            t_in    = int(m.group(1).replace(",", ""))
            t_after = int(m.group(2).replace(",", ""))
        mg = _TZAP_RESULT_G_RE.search(line)
        if mg:
            g_in    = int(mg.group(1).replace(",", ""))
            g_after = int(mg.group(2).replace(",", ""))
        p = _TZAP_PARSE_RE.search(line)
        if p and parse_s is None:
            if g_in is None:
                g_in = int(p.group(1).replace(",", ""))
            parse_s = float(p.group(2))
    return {"t_in": t_in, "t_after": t_after,
            "g_in": g_in, "g_after": g_after,
            "time_s": elapsed, "parse_s": parse_s}

# ── runner: quizx ────────────────────────────────────────────────────────────

QUIZX = "quizx"

# quizx outputs rz gates with angles like "0.25*pi", "-0.75*pi", etc.
# T-count = number of rz gates whose angle is an odd multiple of pi/4.
_QUIZX_RZ_RE = re.compile(r"rz\(\s*(-?\d*\.?\d+)\s*\*\s*pi\s*\)")

def _count_quizx_t(qasm_text: str) -> int:
    count = 0
    for m in _QUIZX_RZ_RE.finditer(qasm_text):
        val = float(m.group(1))  # fraction of pi, e.g. 0.25 or -0.75
        # Normalise to [0, 2) and check if it's an odd multiple of 0.25
        val_mod = val % 2.0
        # Multiples of 0.25 that are NOT multiples of 0.5 are non-Clifford (T-type)
        quarter_steps = round(val_mod / 0.25)
        if quarter_steps % 2 != 0:  # odd => non-Clifford
            count += 1
    return count

def run_quizx(qasm_path: str, timeout: float = 60.0) -> dict:
    try:
        t0 = time.perf_counter()
        r = subprocess.run(
            [QUIZX, "opt", qasm_path],
            capture_output=True, text=True, timeout=timeout,
        )
        elapsed = time.perf_counter() - t0
        if r.returncode != 0:
            return {"error": r.stderr.strip() or "non-zero exit"}
        t_count = _count_quizx_t(r.stdout)
        return {"t_after": t_count, "time_s": elapsed}
    except subprocess.TimeoutExpired:
        return {"error": f"timeout >{timeout:g}s"}
    except FileNotFoundError:
        return {"error": "quizx not found"}
    except Exception as e:
        return {"error": str(e)}

# ── runner: pyzx ─────────────────────────────────────────────────────────────

def run_pyzx(qasm_path: str, timeout: float = 60.0) -> dict:
    # Run in a subprocess so we can enforce a hard timeout.
    script = (
        "import pyzx, sys\n"
        f"c = pyzx.Circuit.load({qasm_path!r})\n"
        "g = c.to_graph()\n"
        "pyzx.full_reduce(g)\n"
        "c2 = pyzx.extract_circuit(g).to_basic_gates()\n"
        "c2 = pyzx.basic_optimization(c2)\n"
        "print(c2.tcount())\n"
    )
    try:
        t0 = time.perf_counter()
        r = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True, text=True, timeout=timeout,
        )
        elapsed = time.perf_counter() - t0
        if r.returncode != 0:
            return {"error": r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "non-zero exit"}
        return {"t_after": int(r.stdout.strip()), "time_s": elapsed}
    except subprocess.TimeoutExpired:
        return {"error": f"timeout >{timeout:g}s"}
    except Exception as e:
        return {"error": str(e)}

# ── runner: voqc ─────────────────────────────────────────────────────────────

VOQC_DEFAULT = os.path.expanduser("~/git/voqc/VOQC/_build/default/voqc.exe")
VOQC = VOQC_DEFAULT if os.path.isfile(VOQC_DEFAULT) else (shutil.which("voqc") or VOQC_DEFAULT)

_VOQC_ORIG_RE  = re.compile(r"Original:.*?\bT\s+(\d+)")
_VOQC_FINAL_RE = re.compile(r"Final:.*?\bT\s+(\d+)")

def run_voqc(qasm_path: str, timeout: float = 60.0) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".qasm", delete=False) as tmp:
        out_path = tmp.name
    try:
        t0 = time.perf_counter()
        r = subprocess.run(
            [VOQC, "-i", qasm_path, "-o", out_path],
            capture_output=True, text=True, timeout=timeout,
        )
        elapsed = time.perf_counter() - t0
        if r.returncode != 0:
            return {"error": (r.stderr.strip().splitlines() or ["non-zero exit"])[-1]}
        m_orig  = _VOQC_ORIG_RE.search(r.stdout)
        m_final = _VOQC_FINAL_RE.search(r.stdout)
        if not m_final:
            return {"error": "no Final line"}
        return {
            "t_in":    int(m_orig.group(1)) if m_orig else None,
            "t_after": int(m_final.group(1)),
            "time_s":  elapsed,
        }
    except subprocess.TimeoutExpired:
        return {"error": f"timeout >{timeout:g}s"}
    except FileNotFoundError:
        return {"error": "voqc not found"}
    except Exception as e:
        return {"error": str(e)}
    finally:
        if os.path.exists(out_path):
            os.unlink(out_path)

# ── runner: feynman (feynopt) ────────────────────────────────────────────────

FEYN_DEFAULT = (
    "/Users/aws/git/feynman/dist-newstyle/build/aarch64-osx/"
    "ghc-9.14.1/Feynman-0.1.0.0/x/feynopt/build/feynopt/feynopt"
)
FEYNMAN = FEYN_DEFAULT if os.path.isfile(FEYN_DEFAULT) else (shutil.which("feynopt") or FEYN_DEFAULT)

# feynopt prints "//   T: <n>" under a "// Result" block and "// Original" block.
_FEYN_RESULT_RE = re.compile(r"// Result.*?(?=^OPENQASM|\Z)", re.DOTALL | re.MULTILINE)
_FEYN_ORIG_RE = re.compile(r"// Original.*?// Result", re.DOTALL)
_FEYN_T_RE = re.compile(r"^//\s+T:\s+(\d+)", re.MULTILINE)

def run_feynman(qasm_path: str, timeout: float = 60.0) -> dict:
    try:
        t0 = time.perf_counter()
        r = subprocess.run(
            [FEYNMAN, "-toCliffordT", "-O2", qasm_path],
            capture_output=True, text=True, timeout=timeout,
        )
        elapsed = time.perf_counter() - t0
        if r.returncode != 0:
            return {"error": (r.stderr.strip().splitlines() or ["non-zero exit"])[-1]}
        m = _FEYN_RESULT_RE.search(r.stdout)
        if not m:
            return {"error": "no Result block"}
        t_match = _FEYN_T_RE.search(m.group(0))
        t_after = int(t_match.group(1)) if t_match else 0
        t_in = None
        om = _FEYN_ORIG_RE.search(r.stdout)
        if om:
            tm = _FEYN_T_RE.search(om.group(0))
            if tm:
                t_in = int(tm.group(1))
        return {"t_in": t_in, "t_after": t_after, "time_s": elapsed}
    except subprocess.TimeoutExpired:
        return {"error": f"timeout >{timeout:g}s"}
    except FileNotFoundError:
        return {"error": "feynopt not found"}
    except Exception as e:
        return {"error": str(e)}

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qasm_dir", nargs="?",
                        default=os.path.join(REPO_ROOT, "qasm"))
    parser.add_argument("--max-kb", type=float, default=2000.0,
                        help="Skip files larger than this many KB (default: 2000)")
    parser.add_argument("--filter", default="*",
                        help="Glob pattern to match circuit names (default: *)")
    parser.add_argument("--tools", default="tzap,quizx,feynman,voqc",
                        help="Comma-separated list of tools to run "
                             "(choices: tzap, quizx, pyzx, feynman, voqc; "
                             "default: tzap,quizx,feynman,voqc)")
    parser.add_argument("--timeout", type=float, default=60.0,
                        help="Per-tool timeout in seconds (default: 60)")
    parser.add_argument("--csv", default="compare_results.csv",
                        help="CSV output path (default: compare_results.csv)")
    args = parser.parse_args()

    all_runners = {
        "tzap":    ("tzap",    run_tzap),
        "quizx":   ("quizx",   run_quizx),
        "pyzx":    ("pyzx",    run_pyzx),
        "feynman": ("feynman", run_feynman),
        "voqc":    ("voqc",    run_voqc),
    }
    selected = [t.strip() for t in args.tools.split(",") if t.strip()]
    for t in selected:
        if t not in all_runners:
            print(f"Unknown tool: {t}. Choices: {list(all_runners)}", file=sys.stderr)
            sys.exit(2)
    runners = [all_runners[t] for t in selected]

    files = sorted(
        p for p in glob.glob(os.path.join(args.qasm_dir, "*.qasm"))
        if os.path.getsize(p) <= args.max_kb * 1024
        and glob.fnmatch.fnmatch(os.path.basename(p), args.filter + ".qasm")
    )
    if not files:
        print(f"No .qasm files ≤ {args.max_kb} KB in {args.qasm_dir}")
        sys.exit(1)

    W = 22
    C = 18
    TW = 12
    header = f"\n{'Circuit':<{W}}  {'T_in':>{TW}}  {'Gates_in':>{TW}}  {'Gates_out':>{TW}}"
    for name, _ in runners:
        header += f"  {name+' (T, %red)':<{C}}    ms"
        if name == "tzap":
            header += "  parse"
    print(header)
    print("─" * (W + 7 + 2*(TW+2) + len(runners) * (C + 10)))

    def fmt(r, t_in, tool=None):
        if "error" in r:
            err = r["error"][:10]
            cell = f"{'ERR:'+err:>12} {'':>5}  "
        else:
            t = r["t_after"]
            ms = r["time_s"] * 1000
            c = f"{t:,} ({(t_in-t)/t_in*100:.0f}%)" if t_in else f"{t:,}"
            cell = f"{c:>{C}} {ms:>5.0f}ms"
        if tool == "tzap":
            parse_s = r.get("parse_s")
            parse_str = f" {parse_s*1000:>5.0f}ms" if parse_s is not None else f" {'':>7}"
            cell += parse_str
        return cell

    # If tzap is not among runners, still count T_in from source qasm
    _T_IN_RE = re.compile(r"^\s*(t|tdg)\b", re.IGNORECASE | re.MULTILINE)

    def count_t_in(path: str) -> int:
        with open(path) as fh:
            return len(_T_IN_RE.findall(fh.read()))

    csv_path = os.path.abspath(args.csv)
    csv_fh = open(csv_path, "w", newline="")
    csv_cols = ["circuit", "t_in", "gates_in", "gates_after"]
    for tname, _ in runners:
        csv_cols += [f"{tname}_t_after", f"{tname}_time_ms"]
        if tname == "tzap":
            csv_cols.append("tzap_parse_ms")
    csv_w = csv.writer(csv_fh)
    csv_w.writerow(csv_cols)

    for path in files:
        name = os.path.splitext(os.path.basename(path))[0]

        results = {}
        for tname, fn in runners:
            results[tname] = fn(path, timeout=args.timeout)

        t_in = None
        for tname in ("tzap", "feynman", "voqc"):
            if tname in results and results[tname].get("t_in") is not None:
                t_in = results[tname]["t_in"]
                break
        if t_in is None:
            t_in = count_t_in(path)

        g_in = g_after = None
        tzap_res = results.get("tzap") or {}
        if tzap_res.get("g_in") is not None:
            g_in = tzap_res["g_in"]
        if tzap_res.get("g_after") is not None:
            g_after = tzap_res["g_after"]

        g_in_s    = f"{g_in:,}"    if g_in    is not None else "-"
        g_after_s = f"{g_after:,}" if g_after is not None else "-"
        row = f"{name:<{W}}  {t_in:>{TW},}  {g_in_s:>{TW}}  {g_after_s:>{TW}}"
        for tname, _ in runners:
            row += "  " + fmt(results[tname], t_in, tool=tname)
        print(row, flush=True)

        csv_row = [name, t_in, g_in, g_after]
        for tname, _ in runners:
            r = results[tname]
            if "error" in r:
                csv_row += [f"ERR:{r['error']}", ""]
            else:
                csv_row += [r.get("t_after"), f"{r['time_s']*1000:.1f}"]
            if tname == "tzap":
                ps = r.get("parse_s")
                csv_row.append(f"{ps*1000:.1f}" if ps is not None else "")
        csv_w.writerow(csv_row)
        csv_fh.flush()

    csv_fh.close()
    print(f"\nResults written to {csv_path}")

if __name__ == "__main__":
    main()
