import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Raw data
data = [
    ("gf2^4_mult",   3,      4),
    ("gf2^5_mult",   3,      5),
    ("gf2^6_mult",   7,     11),
    ("gf2^7_mult",   3,     10),
    ("gf2^8_mult",   3,     16),
    ("gf2^9_mult",   3,     21),
    ("gf2^10_mult",  8,     54),
    ("gf2^16_mult", 20,    297),
    ("gf2^32_mult",  5,  13093),
    ("gf2^64_mult", 10, 920557),
    ("gf2^128_mult",56,   None),  # None = timeout
]

# Sort by exponent in name (2^4, 2^5, ..., 2^128)
import re
data.sort(key=lambda r: int(re.search(r'\^(\d+)', r[0]).group(1)))

circuits  = [r[0] for r in data]
tzap_ms   = [r[1] for r in data]
quizx_ms  = [r[2] for r in data]

TIMEOUT_PLACEHOLDER = 1_500_000

# --- palette (light theme) --------------------------------------------------
TZAP_COL    = "#1971C2"   # deep blue
QUIZX_COL   = "#E03131"   # deep red
TIMEOUT_COL = "#F08C00"   # amber
BG          = "#FFFFFF"
PANEL_BG    = "#F8F9FA"
GRID_COL    = "#DEE2E6"
TEXT_COL    = "#212529"
LABEL_COL   = "#495057"

fig, ax = plt.subplots(figsize=(11, 6.5))
fig.patch.set_facecolor(BG)
ax.set_facecolor(PANEL_BG)

n = len(circuits)
y = np.arange(n)
bar_h = 0.36

# --- quizx bars -------------------------------------------------------------
quizx_vals = [v if v is not None else TIMEOUT_PLACEHOLDER for v in quizx_ms]
is_timeout = [v is None for v in quizx_ms]

bars_quizx = ax.barh(y + bar_h / 2, quizx_vals, bar_h,
                     color=QUIZX_COL, alpha=0.80, zorder=3, label="quizx",
                     linewidth=0)

def fmt_time(ms):
    if ms >= 1000:
        s = ms / 1000
        return f"{s:.1f}s" if s < 100 else f"{s:.0f}s"
    return f"{ms} ms"

def fmt_speedup(tzap, quizx):
    if quizx is None:
        return None
    ratio = quizx / tzap
    if ratio >= 10:
        return f"{ratio:.0f}×"
    return f"{ratio:.1f}×"

for bar, tmo in zip(bars_quizx, is_timeout):
    if tmo:
        bar.set_color(TIMEOUT_COL)
        bar.set_alpha(0.65)
        bar.set_hatch("////")
        bar.set_edgecolor("#CCA300")
        ax.text(TIMEOUT_PLACEHOLDER * 1.05,
                bar.get_y() + bar.get_height() / 2,
                "timeout", va="center", ha="left", fontsize=8.5,
                color=TIMEOUT_COL, fontweight="bold")

# --- tzap bars --------------------------------------------------------------
bars_tzap = ax.barh(y - bar_h / 2, tzap_ms, bar_h,
                    color=TZAP_COL, alpha=0.85, zorder=3, label="tzap",
                    linewidth=0)

# --- value labels -----------------------------------------------------------
for bar, val in zip(bars_tzap, tzap_ms):
    ax.text(val * 1.08, bar.get_y() + bar.get_height() / 2,
            fmt_time(val), va="center", ha="left", fontsize=8,
            color=TZAP_COL, fontweight="semibold")

for bar, val, tmo in zip(bars_quizx, quizx_ms, is_timeout):
    if not tmo:
        ax.text(val * 1.08, bar.get_y() + bar.get_height() / 2,
                fmt_time(val), va="center", ha="left", fontsize=8,
                color=QUIZX_COL, fontweight="semibold")

# --- speedup labels on quizx bars (or at right edge for timeout) ------------
for bar, tzap, quizx, tmo in zip(bars_quizx, tzap_ms, quizx_ms, is_timeout):
    label = fmt_speedup(tzap, quizx)
    if label is None:
        continue
    bar_cx = (bar.get_x() + bar.get_width()) / 2 if not tmo else TIMEOUT_PLACEHOLDER / 2
    bar_cy = bar.get_y() + bar.get_height() / 2
    ax.text(bar_cx, bar_cy, label,
            va="center", ha="center", fontsize=7.5, fontweight="bold",
            color="white", alpha=0.92)

# --- axes -------------------------------------------------------------------
ax.set_xscale("log")
ax.set_xlim(1, 9_000_000)

ax.set_yticks(y)
ax.set_yticklabels(circuits, fontsize=10.5, color=TEXT_COL, fontfamily="monospace")

ax.set_xlabel("Time (ms, log scale)", color=LABEL_COL, fontsize=10)
ax.tick_params(colors=LABEL_COL, which="both", labelsize=9)
ax.tick_params(axis="x", which="minor", length=3)

for spine in ["top", "right", "left"]:
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_color(GRID_COL)

ax.xaxis.grid(True, which="both", color=GRID_COL, linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
ax.invert_yaxis()

# --- legend -----------------------------------------------------------------
timeout_patch = mpatches.Patch(facecolor=TIMEOUT_COL, hatch="////",
                                edgecolor="#CCA300", alpha=0.65,
                                label="quizx (timeout)")
tzap_patch  = mpatches.Patch(color=TZAP_COL,  alpha=0.85, label="tzap")
quizx_patch = mpatches.Patch(color=QUIZX_COL, alpha=0.80, label="quizx")
ax.legend(handles=[tzap_patch, quizx_patch, timeout_patch],
          facecolor=BG, edgecolor=GRID_COL, labelcolor=TEXT_COL,
          fontsize=9.5, loc="upper right", framealpha=1)

ax.set_title("T-gate reduction: tzap vs quizx — runtime",
             color=TEXT_COL, fontsize=13, pad=14, fontweight="bold")

plt.tight_layout()
plt.savefig("scripts/comparison.png", dpi=160, bbox_inches="tight",
            facecolor=fig.get_facecolor())
print("Saved to scripts/comparison.png")
