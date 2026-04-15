import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import re

# --- data -------------------------------------------------------------------
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
    ("gf2^128_mult",56,   None),
]
data.sort(key=lambda r: int(re.search(r'\^(\d+)', r[0]).group(1)))

circuits = [r[0] for r in data]
tzap_ms  = [r[1] for r in data]
quizx_ms = [r[2] for r in data]

TIMEOUT_PLACEHOLDER = 1_500_000

def fmt_time(ms):
    if ms is None:
        return "timeout"
    if ms >= 1000:
        s = ms / 1000
        return f"{s:.0f} s" if s >= 10 else f"{s:.1f} s"
    return f"{ms} ms"

def fmt_speedup(tzap, quizx):
    if quizx is None:
        return None
    r = quizx / tzap
    return f"{r:.0f}×" if r >= 10 else f"{r:.1f}×"

# --- style ------------------------------------------------------------------
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         9,
    "axes.linewidth":    0.7,
    "xtick.major.width": 0.7,
    "xtick.minor.width": 0.4,
})

DARK  = "#1a1a1a"
MID   = "#888888"
LIGHT = "#cccccc"
WHITE = "#ffffff"
GRID  = "#e8e8e8"

fig, ax = plt.subplots(figsize=(7.5, 4.0))
fig.patch.set_facecolor(WHITE)
ax.set_facecolor(WHITE)

n     = len(circuits)
y     = np.arange(n)
bar_h = 0.34

is_timeout = [v is None for v in quizx_ms]
quizx_draw = [v if v is not None else TIMEOUT_PLACEHOLDER for v in quizx_ms]

# quizx bars
bars_q = ax.barh(y + bar_h / 2, quizx_draw, bar_h,
                 color=LIGHT, edgecolor=MID, linewidth=0.5, zorder=3)
for bar, tmo in zip(bars_q, is_timeout):
    if tmo:
        bar.set_hatch("////")

# tzap bars
bars_t = ax.barh(y - bar_h / 2, tzap_ms, bar_h,
                 color=DARK, linewidth=0, zorder=3)

# --- axes (set before annotating so transforms are stable) ------------------
ax.set_xscale("log")
ax.set_xlim(1, 4e7)          # extra right margin for labels
ax.invert_yaxis()

ax.set_yticks(y)
ax.set_yticklabels(circuits, fontsize=8.5, fontfamily="monospace", color=DARK)
ax.set_xlabel("Time (ms, log scale)", fontsize=9, color=DARK)
ax.tick_params(axis="x", colors=DARK, which="both", labelsize=8)
ax.tick_params(axis="y", length=0)

for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_color(DARK)

ax.xaxis.grid(True, which="major", color=GRID, linewidth=0.7, zorder=0)
ax.xaxis.grid(True, which="minor", color=GRID, linewidth=0.35, zorder=0)
ax.set_axisbelow(True)

# --- value labels (just outside each bar) -----------------------------------
for bar, val in zip(bars_t, tzap_ms):
    ax.text(val * 1.15, bar.get_y() + bar.get_height() / 2,
            fmt_time(val), va="center", ha="left", fontsize=7, color=DARK)

for bar, val, tmo in zip(bars_q, quizx_ms, is_timeout):
    x = (TIMEOUT_PLACEHOLDER if tmo else val) * 1.15
    ax.text(x, bar.get_y() + bar.get_height() / 2,
            fmt_time(val), va="center", ha="left", fontsize=7,
            color=MID, style="italic" if tmo else "normal")

# --- speedup column: right-aligned at a fixed x position -------------------
X_SPEEDUP = 1.2e7   # fixed x for speedup column (in data coords)
ax.axvline(X_SPEEDUP * 0.72, color=GRID, linewidth=0.6, zorder=0)   # subtle separator

for i, (tzap, quizx) in enumerate(zip(tzap_ms, quizx_ms)):
    label = fmt_speedup(tzap, quizx)
    if label is None:
        label = "—"
    ax.text(X_SPEEDUP, i, label,
            va="center", ha="center", fontsize=7.5,
            fontweight="bold", color=DARK)

# column header
ax.text(X_SPEEDUP, -0.85, "speedup",
        va="center", ha="center", fontsize=7.5, color=MID,
        style="italic")

# --- legend -----------------------------------------------------------------
timeout_patch = mpatches.Patch(facecolor=LIGHT, hatch="////",
                                edgecolor=MID, linewidth=0.5,
                                label="quizx (timeout)")
tzap_patch  = mpatches.Patch(facecolor=DARK,  linewidth=0,    label="tzap")
quizx_patch = mpatches.Patch(facecolor=LIGHT, edgecolor=MID,
                              linewidth=0.5,                   label="quizx")
ax.legend(handles=[tzap_patch, quizx_patch, timeout_patch],
          facecolor=WHITE, edgecolor=GRID, labelcolor=DARK,
          fontsize=8, loc="upper right", framealpha=1,
          handlelength=1.2, handleheight=0.9,
          bbox_to_anchor=(0.62, 1.0))

ax.set_title("tzap vs. quizx: runtime on GF(2$^k$) multiplier circuits",
             fontsize=10, color=DARK, pad=10)

plt.tight_layout(pad=0.6)
plt.savefig("scripts/comparison.png", dpi=200, bbox_inches="tight",
            facecolor=WHITE)
print("Saved to scripts/comparison.png")
