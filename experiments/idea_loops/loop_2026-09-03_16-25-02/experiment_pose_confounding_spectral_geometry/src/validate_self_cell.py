"""Numerical radial-quadrature validation of the corrected self-cell integral.

Run from the experiment root:
    .venv/bin/python src/validate_self_cell.py

Purpose (PARENT_CORRECTIONS.md, workshop rule 16): independently verify

    I_self = (i/4) * 2 pi * int_0^a r H_0^{(1)}(k_b r) dr
           = (i pi a / (2 k_b)) * H_1^{(1)}(k_b a) - 1/k_b^2,

with equal-area disk radius a = h/sqrt(pi), by Gauss-Legendre radial
quadrature over the disk, for several grid sizes N (h = 1/N).  The output
records the analytic value, the quadrature value, their difference, and the
behaviour of |I_self| under refinement.  The invalidated rule (upper-bound
term only) is retained in the output for comparison and clearly labelled
INVALID so its saturation near 1/k_b^2 can be audited.

Outputs:
    results/self_cell_quadrature_validation.json
    figures/self_cell_quadrature_validation.png
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.special import hankel1

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import helmholtz as hh  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


KB = 2.0 * np.pi
NS = [8, 12, 16, 24, 32, 40, 64, 128]
QUAD_NODES = 2048


def wrong_self_cell_invalid(k_b: float, h: float) -> complex:
    """INVALIDATED rule: (i/4) 2 pi (a/k_b) H_1(k_b a), missing -1/k_b^2."""
    k_b = float(k_b)
    a = float(h) / np.sqrt(np.pi)
    return (1j / 4.0) * (2.0 * np.pi) * (a / k_b) * hankel1(1, k_b * a)


def disk_quadrature(k_b: float, h: float, nodes: int = QUAD_NODES) -> complex:
    """2 pi int_0^a r (i/4) H_0^{(1)}(k_b r) dr on Gauss-Legendre nodes."""
    a = float(h) / np.sqrt(np.pi)
    x, w = np.polynomial.legendre.leggauss(nodes)
    r = 0.5 * a * (x + 1.0)                # [-1, 1] -> [0, a]
    integrand = (1j / 4.0) * (2.0 * np.pi) * r * hankel1(0, k_b * r)
    return complex(np.dot(w, integrand) * (0.5 * a))


def c2l(z: complex) -> list[float]:
    return [float(z.real), float(z.imag)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    t0 = time.perf_counter()
    started_utc = datetime.now(timezone.utc)

    rows = []
    for N in NS:
        h = 1.0 / N
        a = h / np.sqrt(np.pi)
        analytic = hh.self_cell_green(KB, h)
        numeric = disk_quadrature(KB, h)
        wrong = wrong_self_cell_invalid(KB, h)
        abs_err = abs(analytic - numeric)
        rows.append(
            {
                "N": N,
                "h": h,
                "a": a,
                "k_b_a": KB * a,
                "I_self_analytic_re_im": c2l(analytic),
                "I_self_quadrature_re_im": c2l(numeric),
                "quadrature_abs_error": abs_err,
                "quadrature_rel_error": float(
                    abs_err / max(abs(analytic), np.finfo(float).eps)
                ),
                "I_self_invalidated_wrong_re_im": c2l(wrong),
                "abs_I_self": abs(analytic),
                "abs_I_self_invalidated_wrong": abs(wrong),
                "G_D_diag_kb2_I_self_re_im": c2l(KB**2 * analytic),
                "G_D_diag_kb2_I_self_abs": abs(KB**2 * analytic),
                "G_D_diag_kb2_invalidated_wrong_abs": abs(KB**2 * wrong),
            }
        )

    mags = [r["abs_I_self"] for r in rows]
    quad_errs = [r["quadrature_abs_error"] for r in rows]
    tends_to_zero = bool(
        len(mags) >= 2
        and mags[-1] < mags[0] / 50.0
        and all(f < c for c, f in zip(mags, mags[1:]))
        and max(quad_errs) < 1e-11
    )
    summary = {
        "formula": (
            "I_self = (i*pi*a/(2*k_b))*H_1^(1)(k_b*a) - 1/k_b^2, "
            "a = h/sqrt(pi); [G_D]_nn = k_b^2*I_self"
        ),
        "invalidated_rule": (
            "(i/4)*2*pi*(a/k_b)*H_1^(1)(k_b*a)  [omits the -1/k_b^2 "
            "lower-endpoint term; labelled INVALID]"
        ),
        "N_values": NS,
        "quadrature": {
            "method": "Gauss-Legendre radial quadrature on [0, a]",
            "nodes_per_integral": QUAD_NODES,
        },
        "tends_to_zero_under_refinement": tends_to_zero,
        "max_quadrature_abs_error": float(max(quad_errs)),
        "max_abs_I_self_over_ladder": float(max(mags)),
        "first_abs_I_self": float(mags[0]),
        "last_abs_I_self": float(mags[-1]),
        "first_to_last_ratio": float(mags[0] / mags[-1]),
        "k_b": KB,
        "k_b_inv_squared": 1.0 / KB**2,
        "note": (
            "Analytic quadrature agreement is reported raw (no tolerance "
            "threshold applied). The tends-to-zero flag combines monotone "
            "|I_self| decrease over the ladder, a >=50x reduction, and "
            "quadrature consistency <1e-11; it is a finite-ladder "
            "diagnostic, not a continuum-limit proof."
        ),
    }

    runtime_s = time.perf_counter() - t0
    src_sha = sha256(_HERE / "helmholtz.py")
    script_sha = sha256(_HERE / "validate_self_cell.py")
    metadata = {
        "generated_utc": started_utc.isoformat(),
        "command": ".venv/bin/python src/validate_self_cell.py",
        "working_directory": str(_ROOT),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "matplotlib": matplotlib.__version__,
        "runtime_seconds": runtime_s,
        "seeds": None,
        "source_sha256": {
            "src/helmholtz.py": src_sha,
            "src/validate_self_cell.py": script_sha,
        },
        "parent_correction_gate": (
            "context/PARENT_CORRECTIONS.md and workshop rule 16 (2026-09-03)"
        ),
    }
    results = {
        "schema": "self_cell_quadrature_validation",
        "rows": rows,
        "summary": summary,
        "metadata": metadata,
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    results_path = results_dir / "self_cell_quadrature_validation.json"
    results_path.write_text(json.dumps(results, indent=2) + "\n")

    # Figure: magnitude under refinement, wrong rule kept for contrast.
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    ax.loglog(NS, mags, "o-", color="#1f77b4", lw=1.4, ms=5,
              label=r"corrected $|I_{\rm self}|$ (analytic)")
    ax.loglog(NS, [r["quadrature_abs_error"] for r in rows], "x", ms=3,
              color="#7f7f7f", label="quadrature |analytic - numeric|")
    ax.loglog(
        NS,
        [r["abs_I_self_invalidated_wrong"] for r in rows],
        "s--",
        color="#d62728",
        lw=1.2,
        ms=4,
        label="invalidated rule (no $-1/k_b^2$)",
    )
    ax.axhline(1.0 / KB**2, color="#d62728", lw=0.8, ls=":", alpha=0.6)
    ax.annotate(
        r"$1/k_b^2 \approx 2.53\times10^{-2}$ (wrong limit)",
        xy=(NS[-2], 1.0 / KB**2),
        xytext=(NS[-4], 3.5e-2),
        fontsize=8,
        color="#d62728",
    )
    ax.annotate(
        "corrected $I_{\\rm self}\\to 0$ as $h\\to 0$",
        xy=(NS[-2], mags[-2]),
        xytext=(NS[-4], 4e-4),
        fontsize=8,
        color="#1f77b4",
    )
    ax.set_xlabel("grid size N (h = 1/N)")
    ax.set_ylabel("|self-cell disk integral|")
    ax.set_title(
        "Self-cell validation: complete equal-area disk integral\n"
        rf"$g=(i/4)H_0^{{(1)}}(k_b r)$, $k_b=2\pi$, "
        r"$a=h/\sqrt{\pi}$, Gauss-Legendre quadrature"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    figure_path = figures_dir / "self_cell_quadrature_validation.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    print("===== SELF-CELL QUADRATURE VALIDATION =====")
    print("max quadrature abs error over N:", summary["max_quadrature_abs_error"])
    for r in rows:
        print(
            "  N=%3d  |I_self|=%.6e  |wrong|=%.6e  |G_D_diag|=%.6e"
            % (r["N"], r["abs_I_self"], r["abs_I_self_invalidated_wrong"],
               r["G_D_diag_kb2_I_self_abs"])
        )
    print("tends to zero under refinement:", tends_to_zero)
    print("results ->", results_path)
    print("figure  ->", figure_path)


if __name__ == "__main__":
    main()
