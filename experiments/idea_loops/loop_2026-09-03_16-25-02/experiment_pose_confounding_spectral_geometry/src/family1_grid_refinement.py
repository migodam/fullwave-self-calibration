"""Forward-model grid-refinement diagnostic for the Family-1 pilot scene.

Run from the experiment root:
    .venv/bin/python src/family1_grid_refinement.py

Purpose (workshop rule 15 and PLAN_PART_3_TO_7.md, "one refinement
diagnostic"): compare the same physical forward-data configuration -- same
continuous contrast function chi0, same poses, same receivers, same k_b -- at
increasing CPU-feasible grid resolutions N in {16, 24, 32, 40} on Apple
Silicon CPU.  The grid samples the continuous scene at cell centers and
applies the corrected equal-area disk self-cell integral at each resolution.

Terminology guard: this is a *forward-model refinement diagnostic*, not a
continuum-convergence proof.  The discrete model is not a nested restriction
of a single continuum problem, so shrinking successive differences is only
evidence that the implemented discretization responds consistently to h
refinement.  No continuum-transfer claim is made.  N=32 and N=40 are both
completed here (not deferred) and the runtimes are recorded.

Outputs:
    results/family1_grid_refinement.json
    figures/family1_grid_refinement.png
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

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import helmholtz as hh  # noqa: E402
import family1_pilot as fp  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


NS = [16, 24, 32, 40]


def c2l(z: complex) -> list[float]:
    return [float(z.real), float(z.imag)]


def complex_vec_to_list(v: np.ndarray) -> list[list[float]]:
    return [[float(z.real), float(z.imag)] for z in v]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    t_start = time.perf_counter()
    started_utc = datetime.now(timezone.utc)
    cfg = dict(fp.CONFIG)
    poses = fp.build_poses(cfg)
    rx_offsets = np.asarray(cfg["rx_offsets"], dtype=float)
    tx_offset = np.asarray(cfg["tx_offset"], dtype=float)

    per_grid = []
    F_vectors = {}
    for N in NS:
        t_n = time.perf_counter()
        points, h = hh.make_grid(N)
        chi0 = fp.make_chi0(points, cfg)
        ops = hh.build_operators(
            chi0, poses, rx_offsets, tx_offset, N, cfg["k_b"]
        )
        F = hh.forward_measurements(
            chi0, poses, rx_offsets, tx_offset, N, cfg["k_b"]
        )
        runtime_s = time.perf_counter() - t_n
        F_vectors[N] = F
        per_grid.append(
            {
                "N": N,
                "h": 1.0 / N,
                "a": h / np.sqrt(np.pi),
                "I_self_re_im": c2l(hh.self_cell_green(cfg["k_b"], h)),
                "abs_I_self": abs(hh.self_cell_green(cfg["k_b"], h)),
                "G_D_diag_kb2_I_self_abs": abs(
                    cfg["k_b"] ** 2 * hh.self_cell_green(cfg["k_b"], h)
                ),
                "sigma_min_M": ops["sigma_min_M"],
                "M_norm": ops["M_norm"],
                "sigma_min_over_M_norm": ops["sigma_min_over_M_norm"],
                "max_state_residual": ops["max_state_residual"],
                "max_abs_E_tot": ops["max_abs_E_tot"],
                "norm_F": float(np.linalg.norm(F, ord=2)),
                "F_re_im": complex_vec_to_list(F),
                "runtime_seconds": runtime_s,
            }
        )

    # Successive and finest-referenced comparisons on the fixed-dimension
    # physical data vector (T*n_rx = 24 complex entries at shared receivers).
    comparisons = []
    finest = NS[-1]
    F_finest = F_vectors[finest]
    for i, N in enumerate(NS):
        F = F_vectors[N]
        row = {
            "N": N,
            "rel_diff_to_finest": float(
                np.linalg.norm(F - F_finest, ord=2)
                / np.linalg.norm(F_finest, ord=2)
            ),
        }
        if i > 0:
            F_prev = F_vectors[NS[i - 1]]
            row["rel_diff_to_previous"] = float(
                np.linalg.norm(F - F_prev, ord=2) / np.linalg.norm(F, ord=2)
            )
        else:
            row["rel_diff_to_previous"] = None
        comparisons.append(row)

    summary = {
        "discipline": (
            "Forward-model grid-refinement diagnostic only. The grids sample "
            "one continuous two-Gaussian contrast scene with shared "
            "poses/receivers/k_b; shrinking differences under h refinement "
            "are NOT asserted to be continuum convergence."
        ),
        "deferred_grids": [],
        "completed_grids": NS,
        "status_32_and_40": "both completed (not deferred)",
        "self_cell_formula": hh.SELF_CELL_FORMULA,
        "self_cell_formula_version": hh.SELF_CELL_FORMULA_VERSION,
        "n_finest_reference": finest,
        "physical_data_dimension": int(len(F_finest)),
    }

    runtime_total = time.perf_counter() - t_start
    metadata = {
        "generated_utc": started_utc.isoformat(),
        "command": ".venv/bin/python src/family1_grid_refinement.py",
        "working_directory": str(_ROOT),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "matplotlib": matplotlib.__version__,
        "runtime_seconds_total": runtime_total,
        "seeds": None,
        "tolerances": {
            "fd_or_quadrature_tolerance": "none applied; raw numbers recorded",
            "lu_check_finite": False,
        },
        "source_sha256": {
            "src/helmholtz.py": sha256(_HERE / "helmholtz.py"),
            "src/family1_pilot.py": sha256(_HERE / "family1_pilot.py"),
            "src/family1_grid_refinement.py": sha256(
                _HERE / "family1_grid_refinement.py"
            ),
        },
        "parent_correction_gate": (
            "context/PARENT_CORRECTIONS.md and workshop rule 16 (2026-09-03)"
        ),
    }
    results = {
        "schema": "family1_grid_refinement",
        "config": {
            "k_b": float(cfg["k_b"]),
            "poses": [p.tolist() for p in poses],
            "rx_offsets": rx_offsets.tolist(),
            "tx_offset": tx_offset.tolist(),
            "chi_blobs": cfg["chi_blobs"],
            "grid_N": NS,
        },
        "per_grid": per_grid,
        "comparisons": comparisons,
        "summary": summary,
        "metadata": metadata,
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    results_path = results_dir / "family1_grid_refinement.json"
    results_path.write_text(json.dumps(results, indent=2) + "\n")

    # Figure: refinement response of the physical data vector and the
    # corrected domain-propagator diagonal magnitude.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8))
    prev_vals = [c["rel_diff_to_previous"] for c in comparisons[1:]]
    fin_vals = [c["rel_diff_to_finest"] for c in comparisons]
    ax1.loglog(NS[1:], prev_vals, "o-", color="#1f77b4", ms=5, lw=1.3,
               label=r"$||F_N - F_{N-1}||_2/||F_N||_2$")
    ax1.loglog(NS, fin_vals, "s--", color="#2ca02c", ms=5, lw=1.2,
               label=r"$||F_N - F_{40}||_2/||F_{40}||_2$")
    ax1.set_xlabel("grid size N (h = 1/N)")
    ax1.set_ylabel("relative difference")
    ax1.set_title("Forward data under grid refinement\n(discrete-model diagnostic)")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=8)
    diag_abs = [g["G_D_diag_kb2_I_self_abs"] for g in per_grid]
    ax2.loglog(NS, diag_abs, "o-", color="#d62728", ms=5, lw=1.3,
               label=r"$|k_b^2 I_{\rm self}|$ (corrected)")
    ax2.set_xlabel("grid size N (h = 1/N)")
    ax2.set_ylabel(r"$|\,[G_D]_{nn}\,|$")
    ax2.set_title("Corrected self-cell diagonal under refinement\n"
                  "(NOT a continuum-convergence claim)")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    figure_path = figures_dir / "family1_grid_refinement.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    print("===== FAMILY 1 GRID REFINEMENT DIAGNOSTIC =====")
    for r in comparisons:
        print("  N=%3d  rel-to-finetest=%.4e  rel-to-previous=%s"
              % (r["N"], r["rel_diff_to_finest"],
                 "%.4e" % r["rel_diff_to_previous"]
                 if r["rel_diff_to_previous"] is not None else "-"))
    print("completed grids:", NS)
    print("status N=32/40: completed (not deferred)")
    print("results ->", results_path)
    print("figure  ->", figure_path)


if __name__ == "__main__":
    main()
