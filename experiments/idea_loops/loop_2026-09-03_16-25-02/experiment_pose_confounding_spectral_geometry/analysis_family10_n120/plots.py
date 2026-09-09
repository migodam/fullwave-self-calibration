#!/usr/bin/env python
"""Small supplementary figures from the stored Family10 n120 JSON."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "results" / "family10_online_slam_toy_n120.json"
OUT = Path(__file__).resolve().parent / "fig_diagnostics_family10_n120.png"

data = json.loads(JSON_PATH.read_text())
mc = data["monte_carlo"]


def rows(mode: str, kind: str) -> list[dict]:
    return [t[kind] for t in mc[mode]["trial_rows"]]


fig, axes = plt.subplots(2, 2, figsize=(12.5, 10))

for ax, mode in zip(axes[0], ("born", "full_wave")):
    rk = np.array([r["c_error_l2"] for r in rows(mode, "known")])
    rf = np.array([r["c_error_l2"] for r in rows(mode, "free")])
    okk = np.array([r["success"] for r in rows(mode, "known")])
    okf = np.array([r["success"] for r in rows(mode, "free")])
    e2k = rk**2
    e2f = rf**2
    ok = okk & okf
    med_f = np.median(e2f[okf])
    out = okf & (e2f > 20.0 * med_f)

    ax.loglog(
        e2k[ok & ~out], e2f[ok & ~out], "o", ms=4, alpha=0.7,
        label="successful, non-outlier",
    )
    ax.loglog(
        e2k[out], e2f[out], "o", ms=7, color="tab:red",
        label="free outlier (>20x median)",
    )
    fail = okk & ~okf
    ax.loglog(
        e2k[fail], e2f[fail], "x", ms=9, color="black",
        label="free non-success",
    )
    lim = [1e3, 1e8]
    ax.plot(lim, lim, ":", color="gray", lw=1)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"known-pose $\|e\|^2$")
    ax.set_ylabel(r"free-pose $\|e\|^2$")
    ax.set_title(f"{mode}: free vs known map error")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)

# Sorted free squared errors + 20x median line + cumulative sum share.
ax = axes[1][0]
for mode, color in (("born", "tab:blue"), ("full_wave", "tab:orange")):
    okf = np.array([r["success"] for r in rows(mode, "free")])
    e2f = np.sort(np.array([r["c_error_l2"] for r in rows(mode, "free")])[okf] ** 2)
    xs = np.arange(1, e2f.size + 1)
    ax.semilogy(xs, e2f, "-o", ms=3, color=color, label=mode)
    med = np.median(e2f)
    ax.axhline(20.0 * med, color=color, ls="--", lw=0.8)
ax.set_xlabel("trial rank (successful, sorted)")
ax.set_ylabel(r"$\|e\|^2$")
ax.set_title("sorted free map errors (dashed = 20x median)")
ax.grid(alpha=0.3)
ax.legend(fontsize=7)

# Cumulative share of sum-of-squares for born free.
ax = axes[1][1]
okf = np.array([r["success"] for r in rows("born", "free")])
e2f = np.sort(np.array([r["c_error_l2"] for r in rows("born", "free")])[okf] ** 2)[::-1]
cum = np.cumsum(e2f) / e2f.sum()
ax.plot(np.arange(1, e2f.size + 1), cum, "-o", ms=3, color="tab:blue")
top14 = np.where(np.arange(1, e2f.size + 1) == 14)[0][0]
ax.axvline(14, color="gray", ls=":")
ax.axhline(cum[13], color="gray", ls=":")
ax.annotate(
    f"top 14 trials: {100*cum[13]:.1f}% of sum|e|^2",
    xy=(14, cum[13]), xytext=(20, 0.55),
    arrowprops=dict(arrowstyle="->", lw=0.8),
)
ax.set_xlabel("number of largest born-free trials included")
ax.set_ylabel("cumulative share of sum||e||^2")
ax.set_title("born free concentration of squared error")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUT, dpi=160)
print(f"[written] {OUT}")
