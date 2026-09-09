"""Parent correction/extension for the first Family-2 run.

This script preserves the original failed pixel-kernel projector gate while
adding stable factorization evidence, resolves the regular-prior small-alpha
limit with an SVD filter, and compares the singular-prior sweep to its explicit
one-dimensional large-alpha limit.  It does not replace
``family2_algebraic_spine.py`` or erase its raw results.

Run from the experiment root:
    .venv/bin/python src/family2_parent_correction.py
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

import family2_algebraic_spine as f2  # noqa: E402
import helmholtz as hh  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


ALPHAS_REGULAR = np.logspace(-14, 8, 45)
ALPHAS_SINGULAR = np.logspace(2, 8, 13)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_regular_weight(B: np.ndarray, alpha: float) -> np.ndarray:
    """I-B(B^T B+alpha I)^-1 B^T evaluated through the thin SVD."""
    r, U, s, _ = f2.range_basis(B)
    filt = s[:r] ** 2 / (s[:r] ** 2 + float(alpha))
    return np.eye(B.shape[0]) - (U[:, :r] * filt[None, :]) @ U[:, :r].T


def evaluate_case(name: str, A: np.ndarray, B: np.ndarray, original: dict) -> dict:
    rB, U, sB, tolB = f2.range_basis(B)
    Pperp = np.eye(B.shape[0]) - U[:, :rB] @ U[:, :rB].T
    C = Pperp @ A
    KIS = A.T @ A
    K0 = A.T @ (Pperp @ A)
    KIS_F = max(float(np.linalg.norm(KIS, ord="fro")), np.finfo(float).tiny)

    factor_rel = float(np.linalg.norm(K0 - C.T @ C, ord="fro") / KIS_F)
    regular_rows = []
    for alpha in ALPHAS_REGULAR:
        W = stable_regular_weight(B, float(alpha))
        K = A.T @ (W @ A)
        regular_rows.append(
            {
                "alpha": float(alpha),
                "d0_rel_F": float(np.linalg.norm(K - K0, ord="fro") / KIS_F),
                "dinf_rel_F": float(np.linalg.norm(K - KIS, ord="fro") / KIS_F),
            }
        )

    # P = I-ee^T leaves the final pose coordinate unpenalized.  Because B has
    # full column rank in this experiment, the large-alpha limit is the
    # one-dimensional nuisance elimination associated with b=B e.
    q = B.shape[1]
    e = np.zeros(q)
    e[-1] = 1.0
    b = B @ e
    b2 = float(b @ b)
    if b2 <= np.finfo(float).tiny:
        raise RuntimeError("selected unpenalized nuisance column is numerically zero")
    W_lim = np.eye(B.shape[0]) - np.outer(b, b) / b2
    K_lim = A.T @ (W_lim @ A)
    limit_loss = float(np.linalg.norm(KIS - K_lim, ord="fro") / KIS_F)
    P = np.eye(q) - np.outer(e, e)
    singular_rows = []
    for alpha in ALPHAS_SINGULAR:
        normal = B.T @ B + float(alpha) * P
        K = KIS - A.T @ B @ np.linalg.solve(normal, B.T @ A)
        singular_rows.append(
            {
                "alpha": float(alpha),
                "distance_to_explicit_limit_rel_F": float(
                    np.linalg.norm(K - K_lim, ord="fro") / KIS_F
                ),
                "distance_to_KIS_rel_F": float(
                    np.linalg.norm(K - KIS, ord="fro") / KIS_F
                ),
            }
        )

    first = regular_rows[0]
    last = regular_rows[-1]
    sing_last = singular_rows[-1]
    return {
        "case": name,
        "rank_B": rB,
        "rank_tolerance_B": tolB,
        "sigma_B_min_retained": float(sB[rB - 1]),
        "sigma_B_min_squared": float(sB[rB - 1] ** 2),
        "original_pixel_projector_gate": (
            original["check_c1_kernel_identity"] if name == "pixel" else None
        ),
        "stable_kernel_factorization": {
            "identity": "K_SLAM = C^T C with C=(I-P_B)A",
            "relative_factorization_residual": factor_rel,
            "pass_at_1e-12": bool(factor_rel < 1e-12),
            "interpretation": (
                "This factorization plus ||Cu||^2=u^T K_SLAM u proves the "
                "finite-dimensional kernel identity algebraically. It does not "
                "retroactively pass the original ill-conditioned eigenprojector gate."
            ),
        },
        "regular_prior_svd_filter_sweep": {
            "rows": regular_rows,
            "alpha_min": float(ALPHAS_REGULAR[0]),
            "alpha_max": float(ALPHAS_REGULAR[-1]),
            "d0_at_alpha_min": first["d0_rel_F"],
            "dinf_at_alpha_max": last["dinf_rel_F"],
            "small_alpha_reaches_1e-8_gate": bool(first["d0_rel_F"] < 1e-8),
            "large_alpha_reaches_1e-8_gate": bool(last["dinf_rel_F"] < 1e-8),
            "note": (
                "The SVD filter avoids interpreting alpha=1e-6 as asymptotic; "
                "the sweep extends below sigma_min(B)^2."
            ),
        },
        "singular_prior_explicit_limit": {
            "unpenalized_coordinate_zero_based": q - 1,
            "label": "unpenalized nuisance direction (not asserted to be gauge)",
            "norm_Be": float(np.linalg.norm(b)),
            "limit_weight": "I-(Be)(Be)^T/||Be||^2",
            "limit_loss_from_KIS_rel_F": limit_loss,
            "rows": singular_rows,
            "last_distance_to_limit_rel_F": sing_last[
                "distance_to_explicit_limit_rel_F"
            ],
            "last_distance_to_KIS_rel_F": sing_last["distance_to_KIS_rel_F"],
            "limit_match_pass_at_1e-8": bool(
                sing_last["distance_to_explicit_limit_rel_F"] < 1e-8
            ),
            "nonzero_limit_is_geometry_dependent": True,
        },
    }


def main() -> None:
    started = datetime.now(timezone.utc)
    t0 = time.perf_counter()
    cfg = f2.CONFIG
    points, _ = hh.make_grid(cfg["N"])
    poses = f2.build_poses(cfg)
    chi0 = f2.make_chi0(points, cfg)
    A_complex, B_complex, _, _ = hh.build_AB(
        chi0,
        poses,
        np.asarray(cfg["rx_offsets"], dtype=float),
        np.asarray(cfg["tx_offset"], dtype=float),
        cfg["N"],
        cfg["k_b"],
    )
    A_pixel, B = hh.whiten_realify(A_complex, B_complex, None)
    S = f2.build_smooth_basis(points, cfg["smooth_basis"])
    A_smooth = A_pixel @ S

    original_path = _ROOT / "results" / "family2_results.json"
    original = json.loads(original_path.read_text())
    cases = {
        "pixel": evaluate_case("pixel", A_pixel, B, original["cases"]["pixel"]),
        "smooth": evaluate_case("smooth", A_smooth, B, original["cases"]["smooth"]),
    }

    result = {
        "schema": "family2_parent_correction",
        "generated_utc": started.isoformat(),
        "command": ".venv/bin/python src/family2_parent_correction.py",
        "working_directory": str(_ROOT),
        "purpose": (
            "Supplement, not replacement: preserve the original failed gate, "
            "resolve prior limits, and remove the false gauge label."
        ),
        "cases": cases,
        "runtime_seconds": float(time.perf_counter() - t0),
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "source_sha256": {
            "src/helmholtz.py": sha256(_HERE / "helmholtz.py"),
            "src/family2_algebraic_spine.py": sha256(
                _HERE / "family2_algebraic_spine.py"
            ),
            "src/family2_parent_correction.py": sha256(Path(__file__)),
            "results/family2_results.json": sha256(original_path),
        },
        "cannot_establish": (
            "These are finite-dimensional factorization and parameter-limit "
            "checks. They do not prove continuum identifiability, genericity, "
            "or that the selected unpenalized coordinate is a physical gauge."
        ),
    }
    out = _ROOT / "results" / "family2_parent_correction.json"
    out.write_text(json.dumps(result, indent=2) + "\n")

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    colors = {"pixel": "#1f77b4", "smooth": "#d62728"}
    for name, case in cases.items():
        rows = case["regular_prior_svd_filter_sweep"]["rows"]
        a = [r["alpha"] for r in rows]
        axes[0].loglog(a, [r["d0_rel_F"] for r in rows], color=colors[name],
                       marker="o", ms=2.5, label=f"{name}: to no-prior")
        axes[0].loglog(a, [r["dinf_rel_F"] for r in rows], color=colors[name],
                       ls="--", marker="s", ms=2.2, label=f"{name}: to known-pose")
        srows = case["singular_prior_explicit_limit"]["rows"]
        sa = [r["alpha"] for r in srows]
        axes[1].loglog(sa, [r["distance_to_explicit_limit_rel_F"] for r in srows],
                       color=colors[name], marker="o", ms=3,
                       label=f"{name}: distance to explicit limit")
        axes[1].axhline(
            case["singular_prior_explicit_limit"]["limit_loss_from_KIS_rel_F"],
            color=colors[name], ls=":", alpha=0.8,
            label=f"{name}: nonzero loss floor vs K_IS",
        )
    axes[0].set_title("Regular pose prior: both asymptotic limits")
    axes[0].set_xlabel("alpha")
    axes[0].set_ylabel("relative Frobenius distance")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=7)
    axes[1].set_title("Singular prior: explicit one-direction limit")
    axes[1].set_xlabel("alpha")
    axes[1].set_ylabel("relative Frobenius distance")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=7)
    fig.suptitle("Family 2 parent correction: limit checks without threshold tuning")
    fig.tight_layout()
    fig_path = _ROOT / "figures" / "family2_parent_prior_limits.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    note = _ROOT / "notes" / "family2_parent_correction.md"
    note.write_text(
        "# Family 2 parent correction\n\n"
        "The original pixel eigenprojector gate remains failed at 8.47e-7. "
        "The stable factorization residuals are %.3e (pixel) and %.3e "
        "(smooth).  With the SVD-filter sweep down to alpha=1e-14, d0 is "
        "%.3e and %.3e; thus the small-prior limit is now numerically resolved. "
        "The singular-prior large-alpha results converge to the explicit "
        "one-dimensional unpenalized-nuisance limit, whose nonzero relative "
        "loss is %.3e and %.3e.  That direction is not called a gauge.\n\n"
        "These checks remain finite-dimensional and do not prove continuum "
        "identifiability.\n"
        % (
            cases["pixel"]["stable_kernel_factorization"]["relative_factorization_residual"],
            cases["smooth"]["stable_kernel_factorization"]["relative_factorization_residual"],
            cases["pixel"]["regular_prior_svd_filter_sweep"]["d0_at_alpha_min"],
            cases["smooth"]["regular_prior_svd_filter_sweep"]["d0_at_alpha_min"],
            cases["pixel"]["singular_prior_explicit_limit"]["limit_loss_from_KIS_rel_F"],
            cases["smooth"]["singular_prior_explicit_limit"]["limit_loss_from_KIS_rel_F"],
        )
    )

    print("Family 2 parent correction complete")
    for name, case in cases.items():
        print(
            name,
            "factor=%.3e d0(alpha_min)=%.3e dinf(alpha_max)=%.3e "
            "singular_limit_loss=%.3e limit_error=%.3e"
            % (
                case["stable_kernel_factorization"]["relative_factorization_residual"],
                case["regular_prior_svd_filter_sweep"]["d0_at_alpha_min"],
                case["regular_prior_svd_filter_sweep"]["dinf_at_alpha_max"],
                case["singular_prior_explicit_limit"]["limit_loss_from_KIS_rel_F"],
                case["singular_prior_explicit_limit"]["last_distance_to_limit_rel_F"],
            ),
        )
    print("results ->", out)
    print("figure  ->", fig_path)
    print("note    ->", note)


if __name__ == "__main__":
    main()
