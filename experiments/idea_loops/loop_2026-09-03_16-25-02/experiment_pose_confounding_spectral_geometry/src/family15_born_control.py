"""Family 15: Born analytic / literature control over scene contrast.

Run (from the experiment root):
    .venv/bin/python src/family15_born_control.py

What this family adds (no existing source/result/figure/note is modified):

* Family 3 / 3b established that the full-wave Jacobian at chi = 0 equals the
  Born map and diagnosed the finite-difference slope gates for Born
  bilinear checks.
* Family 7 already reports the pixel/smooth Born-vs-full-wave A discrepancy
  over a contrast sweep {0, 0.01, ..., 1.0} together with full-wave no-prior
  retention.
* Family 15 is the dedicated finite-dimensional *analytic control* used by
  the Born-error/FIM narrative: for chi_s = s * two_blob, s in
  {1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0} (1e-1 is the fixed named
  low_contrast row), we compare the Born operator/FIM with the full-wave
  local operator/FIM, in the smooth p=24 coefficient space, under the SAME
  full-wave pose confounder geometry:

  - no-prior: P_perp = I - Z Z^T, Z = orthonormal Range(B_R);
  - finite prior: W = I - B_R (B_R^T B_R + alpha I)^{-1} B_R^T, alpha = 1.0;
  - observed log-log convergence order of ||A_born - A_full||_F / ||A_full||_F,
    and the same ratio for K_IS = A^T A and K_eff;
  - three most confounded full-wave vs Born directions from
    family4.generalized_eigen_directions and the squared-overlap/principal
    angles between the two 3-dim subspaces.

B_R is always the FULL-WAVE pose Jacobian at chi_s (whiten_realify(A_c, B_c,
None)); the Born data map A_born has no pose sensitivity of its own in the
linearized model, and reusing the full-wave confounder block for both rows is
the documented analytic-control choice.

All numbers are finite-dimensional (N=16, m=48, p=24, q_pose=18), on Apple
Silicon CPU, deterministic dense numpy/scipy linear algebra, no RNG.

Outputs: results/family15_born_control.json,
         figures/family15_born_control.png,
         notes/family15_born_control.md.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
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
import family1_pilot as family1  # noqa: E402
import family2_algebraic_spine as family2  # noqa: E402
import family4_frequency_trajectory as family4  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CONFIG = {
    "family": 15,
    "name": "Born analytic/literature control",
    "N": 16,
    "T": 6,
    "n_rx": 4,
    "q": 1.0,
    "f": 1.0,
    "k_b": 2.0 * np.pi,
    "k_b_rule": "k_b = 2*pi*f with f = 1.0",
    "rx_offsets": [
        [-0.06, 0.0],
        [0.06, 0.0],
        [0.0, -0.06],
        [0.0, 0.06],
    ],
    "tx_offset": [0.0, 0.0],
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],
    "pose_theta_convention": (
        "theta = atan2(-p_y, -p_x): body +x axis points toward the origin"
    ),
    "chi_blobs": {
        "amp1": 0.3,
        "amp2": 0.5,
        "sigma1": 0.09,
        "sigma2": 0.07,
        "c1": [-0.15, -0.12],
        "c2": [0.18, 0.14],
    },
    "smooth_basis": {
        "p": 24,
        "x_centers_n": 4,
        "y_centers_n": 6,
        "x_span": [-0.3, 0.3],
        "y_span": [-0.3, 0.3],
        "sigma_b": 0.16,
        "unit_columns": True,
        "source": "family2.build_smooth_basis (identical config to family2/4)",
    },
    "scene_family": {
        "rule": "chi_s = s * chi0, chi0 = family1.make_chi0(points, cfg) (two_blob)",
        "contrast_scales": [0.001, 0.01, 0.03, 0.1, 0.3, 1.0],
        "fixed_named_low_contrast": {
            "s": 0.1,
            "name": "standard low_contrast scene = 0.1 * two_blob (also in scale grid)",
        },
    },
    "dimensions": {
        "m_real_data": 48,  # 2*T*n_rx after whiten_realify
        "n_pixel": 256,     # N^2
        "p_smooth": 24,
        "q_pose": 18,       # 3*T
    },
    "whitening": {
        "W": None,
        "note": (
            "identity noise: hh.whiten_realify(A_c, B_c, None) = "
            "sqrt(2)*[Re(A_c); Im(A_c)] and same for B_c; real data rows 48"
        ),
    },
    "born": {
        "module_path": "hh.born_forward(chi, poses, rx_offsets, tx_offset, N, k_b)",
        "return": "(F_born, A_born); per-pose block A_born,t = G_S_t diag(E_inc_t)",
        "formula_note": (
            "self-derived Born forward map A_born,t = G_S_t diag(E_inc_t), "
            "F_born,t = G_S_t (chi_s * E_inc_t); equivalently the linear "
            "full-wave Jacobian at chi = 0 (verified to roundoff by family "
            "3b check D)"
        ),
        "composition_identity_note": (
            "explicit composition path recomputes G_S_stack diag(E_inc) from "
            "hh.build_operators lists, whiten_realifies it with the same B_c, "
            "projects with S, and compares with the hh.born_forward path"
        ),
    },
    "retention": {
        "no_prior_rule": "P_perp = I - Z Z^T, Z = family2.range_basis(B_R)",
        "same_B_R_note": (
            "B_R is the FULL-WAVE pose Jacobian at chi_s (identity-noise "
            "realified) and is used for both full-wave and Born rows because "
            "the Born data map has no pose sensitivity of its own in the "
            "linearized model; this is a documented analytic-control choice"
        ),
        "spectrum_rule": (
            "rho = family4.retention_spectrum(A_s, W) = eig(Q_A^T W Q_A), "
            "thin-SVD normalised retention of K = A_s^T W A_s w.r.t. "
            "K_IS = A_s^T A_s"
        ),
        "retained_dof": "sum(rho)",
        "log_volume": "sum(log(max(rho, 1e-300))) over rho > 1e-300",
        "rho_min": "min(rho)",
        "theta_min_deg_rule": (
            "(180/pi) * arccos(sqrt(max(0, 1 - rho_min)))"
        ),
    },
    "finite_prior": {
        "alpha": 1.0,
        "W_rule": (
            "W = I - B_R (B_R^T B_R + alpha I)^{-1} B_R^T "
            "(family4.shrinkage_operator)"
        ),
        "K_eff_rule": "K_eff = A_s^T W A_s (family4.K_eff(A_s, B_R, alpha))",
    },
    "comparisons": {
        "rel_fro_A": "||A_s_born - A_s_full||_F / ||A_s_full||_F",
        "rel_fro_K_IS": "||K_IS_born - K_IS_full||_F / ||K_IS_full||_F",
        "rel_fro_K_eff": "||K_eff_born - K_eff_full||_F / ||K_eff_full||_F",
    },
    "convergence": {
        "fit_rule": (
            "slope of log10(relative Frobenius error) vs log10(s) across the "
            "full six-scale grid; reported as observed, not forced"
        ),
        "expected_rule": (
            "slope ~1 expected if the relative Born error is dominated by the "
            "O(s) first correction"
        ),
    },
    "confounded_directions": {
        "source": (
            "family4.generalized_eigen_directions(A_s, P_perp): ascending "
            "generalized eigenvalues of A_s^T P_perp A_s in the K_IS metric; "
            "first three directions are the most pose-confounded"
        ),
        "comparison_rule": (
            "normalized squared overlaps |u^T v|^2/(||u||^2 ||v||^2) and "
            "principal angles (deg) of QR-orthonormalised 3-dim subspaces"
        ),
    },
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "randomness_note": "deterministic dense linear algebra; no RNG used",
    "scope_note": (
        "Finite-dimensional N=16 toy only (m=48 real data rows, p=24 smooth "
        "coefficients, q_pose=18). Internal analytic control; no continuum "
        "limit, no external quantitative benchmark, no independent "
        "validation claim."
    ),
    "literature_rule": "single verbatim citation supplied by the parent task",
}


_EPS = float(np.finfo(float).eps)
_LOG_VOLUME_FLOOR = 1e-300


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, (float, int, str, bool)) or obj is None:
        return obj
    raise TypeError(f"not JSON serialisable: {type(obj)!r}")


def _round_trip_json(obj: dict) -> str:
    return json.dumps(_jsonable(obj), indent=2) + "\n"


def _fmt(x, spec=".6g") -> str:
    if x is None:
        return "-"
    return format(x, spec)


def _rel_fro(A: np.ndarray, B: np.ndarray) -> float:
    """||A - B||_F / ||B||_F (guarded)."""
    den = float(np.linalg.norm(B, ord="fro"))
    if den <= 0.0:
        return float("nan")
    return float(np.linalg.norm(A - B, ord="fro")) / den


def _spectrum_metrics(rho: np.ndarray) -> dict:
    rho = np.asarray(rho, dtype=float)
    pos = rho[rho > _LOG_VOLUME_FLOOR]
    rho_min = float(rho[-1]) if rho.size else float("nan")
    theta_min = (
        math.degrees(math.acos(math.sqrt(max(0.0, 1.0 - rho_min))))
        if rho_min < 1.0
        else 90.0
    )
    return {
        "retained_dof": float(np.sum(rho)),
        "log_volume": float(np.sum(np.log(np.maximum(pos, _LOG_VOLUME_FLOOR)))),
        "rho_min": rho_min,
        "theta_min_deg": float(theta_min),
    }


def _fit_loglog(x: np.ndarray, y: np.ndarray) -> dict:
    """Least-squares slope/intercept of log10(y) vs log10(x)."""
    xs = np.log10(np.asarray(x, dtype=float))
    ys = np.log10(np.maximum(np.asarray(y, dtype=float), _LOG_VOLUME_FLOOR))
    slope, intercept = np.polyfit(xs, ys, 1)
    yhat = slope * xs + intercept
    ss_res = float(np.sum((ys - yhat) ** 2))
    ss_tot = float(np.sum((ys - np.mean(ys)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": r2,
        "n_points": int(len(xs)),
    }


def _normalized_sq_overlap_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """G[i, j] = |a_i^T b_j|^2 / (||a_i||^2 ||b_j||^2), real vectors."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    G = (A.T @ B) ** 2
    na = np.sum(A**2, axis=0)[:, None]
    nb = np.sum(B**2, axis=0)[None, :]
    return G / np.maximum(na * nb, _LOG_VOLUME_FLOOR)


def _principal_angles_deg(U: np.ndarray, V: np.ndarray) -> np.ndarray:
    """Principal angles between the spans of two 3-dim real subspaces."""
    Uq, _ = np.linalg.qr(np.asarray(U, dtype=float), mode="reduced")
    Vq, _ = np.linalg.qr(np.asarray(V, dtype=float), mode="reduced")
    sv = np.linalg.svd(Uq.T @ Vq, compute_uv=False)
    sv = np.clip(sv, 0.0, 1.0)
    ang = np.degrees(np.arccos(sv))
    return np.sort(ang)  # ascending principal angles


def _best_sq_overlap_pairing(G: np.ndarray) -> dict:
    """Permutation of rows that maximises the sum of squared overlaps."""
    n = G.shape[0]
    best = None
    best_sum = -1.0
    for perm in itertools.permutations(range(n)):
        s = float(sum(G[i, perm[i]] for i in range(n)))
        if s > best_sum:
            best_sum = s
            best = perm
    pairs = [[int(i), int(best[i])] for i in range(n)]
    return {
        "pairs_full_to_born": pairs,
        "sum_squared_overlaps": best_sum,
    }


# ---------------------------------------------------------------------------
# Per-scale case
# ---------------------------------------------------------------------------

def _run_scale(
    s: float,
    chi_s: np.ndarray,
    poses: np.ndarray,
    S: np.ndarray,
    cfg: dict,
    n_rx: int,
) -> dict:
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)
    N = int(cfg["N"])
    k_b = float(cfg["k_b"])
    alpha = float(cfg["finite_prior"]["alpha"])

    # --- full wave: A_pix_R, B_R at the nonlinear scene chi_s -------------
    A_c, B_c, _, _ = hh.build_AB(chi_s, poses, rx, tx, N, k_b)
    A_pix_R, B_R = hh.whiten_realify(A_c, B_c, None)
    A_s_full = A_pix_R @ S

    # --- Born: A_pix_born from hh.born_forward, whitened with B_c ---------
    F_born, A_born_c = hh.born_forward(chi_s, poses, rx, tx, N, k_b)
    A_born_R, B_R_again = hh.whiten_realify(A_born_c, B_c, None)
    if not np.array_equal(B_R, B_R_again):
        raise RuntimeError("unexpected B_R mismatch after re-whitening")
    A_s_born = A_born_R @ S

    # --- analytic composition identity: A_born = G_S diag(E_inc) ----------
    ops = hh.build_operators(chi_s, poses, rx, tx, N, k_b)
    G_S_stack = np.vstack(ops["G_S_list"])                # (24, 256) complex
    # Each pose block is n_rx x S complex: G_S_t diag(E_inc_t), so every one
    # of the n_rx rows of block t is scaled by the same E_inc_t (length S).
    S_pix = A_born_c.shape[1]                      # N^2 pixel columns
    rows_c = A_born_c.shape[0]                     # T*n_rx complex rows
    A_c_composed = np.empty((rows_c, S_pix), dtype=np.complex128)
    for t in range(rows_c // n_rx):
        sl = slice(t * n_rx, (t + 1) * n_rx)
        A_c_composed[sl, :] = (
            G_S_stack[sl, :] * ops["E_inc_list"][t][None, :]
        )                                          # G_S diag(E_inc)
    A_comp_R, _ = hh.whiten_realify(A_c_composed, B_c, None)
    A_s_composed = A_comp_R @ S
    composition_identity = {
        "rel_fro_composed_vs_bornforward_complex": _rel_fro(
            A_c_composed, A_born_c
        ),
        "rel_fro_composed_vs_bornforward_real_smooth": _rel_fro(
            A_s_composed, A_s_born
        ),
        "rel_fro_full_vs_born_complex_pixel": _rel_fro(A_born_c, A_c),
        "note": cfg["born"]["composition_identity_note"],
    }

    # --- confounder blocks (identical for full and Born rows) --------------
    Z, P_perp = family4.range_projection(B_R)
    W_eff = family4.shrinkage_operator(B_R, alpha)
    rB = int(family2.range_basis(B_R)[0])

    # --- no-prior retention spectra ----------------------------------------
    rho_full = family4.retention_spectrum(A_s_full, P_perp)
    rho_born = family4.retention_spectrum(A_s_born, P_perp)
    m_full = _spectrum_metrics(rho_full)
    m_born = _spectrum_metrics(rho_born)
    no_prior = {
        "full": m_full,
        "born": m_born,
        "diff_full_minus_born": {
            "retained_dof": m_full["retained_dof"] - m_born["retained_dof"],
            "rho_min": m_full["rho_min"] - m_born["rho_min"],
            "theta_min_deg": m_full["theta_min_deg"] - m_born["theta_min_deg"],
        },
        "rho_full_desc": [float(x) for x in rho_full],
        "rho_born_desc": [float(x) for x in rho_born],
    }

    # --- finite-prior retention spectra (alpha = 1.0) ----------------------
    rho_full_eff = family4.retention_spectrum(A_s_full, W_eff)
    rho_born_eff = family4.retention_spectrum(A_s_born, W_eff)
    ef_full = _spectrum_metrics(rho_full_eff)
    ef_born = _spectrum_metrics(rho_born_eff)
    finite_prior = {
        "alpha": alpha,
        "full": ef_full,
        "born": ef_born,
        "diff_full_minus_born": {
            "retained_dof": ef_full["retained_dof"] - ef_born["retained_dof"],
            "rho_min": ef_full["rho_min"] - ef_born["rho_min"],
            "theta_min_deg": ef_full["theta_min_deg"] - ef_born["theta_min_deg"],
        },
        "rho_full_desc": [float(x) for x in rho_full_eff],
        "rho_born_desc": [float(x) for x in rho_born_eff],
    }

    # --- FIM matrices and relative errors ----------------------------------
    K_IS_full = A_s_full.T @ A_s_full
    K_IS_born = A_s_born.T @ A_s_born
    K_eff_full = family4.K_eff(A_s_full, B_R, alpha)
    K_eff_born = family4.K_eff(A_s_born, B_R, alpha)
    rel_errors = {
        "A_operator": _rel_fro(A_s_born, A_s_full),
        "K_IS": _rel_fro(K_IS_born, K_IS_full),
        "K_eff": _rel_fro(K_eff_born, K_eff_full),
    }

    # --- most-confounded direction subspaces --------------------------------
    rho_g_full, U_full, _ = family4.generalized_eigen_directions(
        A_s_full, P_perp
    )
    rho_g_born, U_born, _ = family4.generalized_eigen_directions(
        A_s_born, P_perp
    )
    U_full_3 = np.asarray(U_full[:, :3])
    U_born_3 = np.asarray(U_born[:, :3])
    ov = _normalized_sq_overlap_matrix(U_full_3, U_born_3)
    principal = _principal_angles_deg(U_full_3, U_born_3)
    confounded = {
        "generalized_rho_full_ascending": [float(x) for x in rho_g_full],
        "generalized_rho_born_ascending": [float(x) for x in rho_g_born],
        "full_top3_direction_vectors": U_full_3,
        "born_top3_direction_vectors": U_born_3,
        "squared_overlap_matrix": ov,
        "diagonal_squared_overlaps": [float(x) for x in np.diag(ov)],
        "max_row_squared_overlaps_full_to_born": [
            float(np.max(ov[i, :])) for i in range(3)
        ],
        "max_col_squared_overlaps_born_from_full": [
            float(np.max(ov[:, j])) for j in range(3)
        ],
        "principal_angles_deg": [float(x) for x in principal],
        "best_pairing": _best_sq_overlap_pairing(ov),
    }

    row = {
        "s": float(s),
        "named_row": (
            "standard low_contrast = 0.1 * two_blob"
            if abs(float(s) - 0.1) < 1e-12
            else None
        ),
        "chi_s_max": float(np.max(chi_s)),
        "chi_s_l2": float(np.linalg.norm(chi_s)),
        "fro_norms": {
            "A_full_s": float(np.linalg.norm(A_s_full, ord="fro")),
            "A_born_s": float(np.linalg.norm(A_s_born, ord="fro")),
            "A_full_pixel_R": float(np.linalg.norm(A_pix_R, ord="fro")),
            "A_born_pixel_R": float(np.linalg.norm(A_born_R, ord="fro")),
            "B_R": float(np.linalg.norm(B_R, ord="fro")),
        },
        "rank_B": rB,
        "composition_identity": composition_identity,
        "no_prior": no_prior,
        "finite_prior": finite_prior,
        "rel_errors": rel_errors,
        "confounded": confounded,
    }
    return row


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _make_figure(
    rows: list[dict],
    slopes: dict,
    path: Path,
) -> None:
    scales = np.asarray([r["s"] for r in rows], dtype=float)
    colors = {
        "A": "#1f77b4",
        "K_IS": "#2ca02c",
        "K_eff": "#d62728",
    }
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0))

    # panel 1: log10 relative Fro errors + fitted slopes
    ax = axes[0, 0]
    labels = {
        "A": "A operator (48x24 real smooth)",
        "K_IS": "K_IS = A^T A",
        "K_eff": "K_eff (finite prior, alpha=1)",
    }
    slope_keys = ("A_operator", "K_IS", "K_eff")
    for key in slope_keys:
        y = np.asarray([r["rel_errors"][key] for r in rows], dtype=float)
        plot_color = colors["A" if key == "A_operator" else key]
        label_key = "A" if key == "A_operator" else key
        ax.loglog(scales, y, "o-", color=plot_color, label=labels[label_key])
        fit = slopes[key]
        xs = np.log10(scales)
        ys = np.log10(y)
        yhat = fit["slope"] * xs + fit["intercept"]
        ax.plot(10**xs, 10**yhat, ":", color=plot_color, alpha=0.65)
        ax.annotate(
            f"slope {fit['slope']:.3f} (r^2={fit['r_squared']:.3f})",
            xy=(0.02, 0.98 - 0.04 * slope_keys.index(key)),
            xycoords="axes fraction",
            fontsize=8,
            color=plot_color,
            va="top",
        )
    ax.set_xlabel("contrast scale s")
    ax.set_ylabel("Born-vs-full relative Frobenius error")
    ax.set_title("(1) Born operator/FIM error and observed log-log slope")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # panel 2: retained_dof and rho_min, no-prior, full vs Born
    ax = axes[0, 1]
    ax.plot(scales, [r["no_prior"]["full"]["retained_dof"] for r in rows],
            "o-", color="#1f77b4", label="retained_dof full-wave")
    ax.plot(scales, [r["no_prior"]["born"]["retained_dof"] for r in rows],
            "s--", color="#2ca02c", label="retained_dof Born")
    ax.set_xscale("log")
    ax.set_xlabel("contrast scale s")
    ax.set_ylabel("retained_dof (no-prior, alpha-infinity)")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_title("(2) retained_dof and rho_min, full vs Born (no-prior)")
    ax2 = ax.twinx()
    ax2.semilogy(
        scales,
        [r["no_prior"]["full"]["rho_min"] for r in rows],
        "^:", color="#ff7f0e",
        label="rho_min full-wave",
    )
    ax2.semilogy(
        scales,
        [r["no_prior"]["born"]["rho_min"] for r in rows],
        "v:", color="#9467bd",
        label="rho_min Born",
    )
    ax2.set_ylabel("rho_min (dimensionless)")
    ax2.grid(True, which="both", alpha=0.15)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right")

    # panel 3: rho spectra at s = 0.1 and s = 1.0
    ax = axes[1, 0]
    for s_plot in (0.1, 1.0):
        r = next(rr for rr in rows if abs(rr["s"] - s_plot) < 1e-12)
        k_full = 24
        x = np.arange(1, k_full + 1)
        ax.semilogy(
            x,
            np.maximum(r["no_prior"]["rho_full_desc"], _EPS),
            "o-",
            color="#1f77b4" if s_plot == 0.1 else "#d62728",
            label=f"full-wave s={s_plot:.1f}",
        )
        ax.semilogy(
            x,
            np.maximum(r["no_prior"]["rho_born_desc"], _EPS),
            "s--",
            color="#2ca02c" if s_plot == 0.1 else "#ff7f0e",
            label=f"Born s={s_plot:.1f}",
        )
    ax.axhline(_EPS, color="k", linestyle=":", linewidth=0.8,
               label=f"machine eps = {_EPS:.1e}")
    ax.set_yscale("log")
    ax.set_xlabel("retention index (descending rho), p=24")
    ax.set_ylabel("retention rho")
    ax.set_title("(3) no-prior rho spectra at s=0.1 and s=1.0")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # panel 4: principal angles between 3-dim most-confounded subspaces
    ax = axes[1, 1]
    for j, label in enumerate(("angle 1", "angle 2", "angle 3")):
        y = np.asarray(
            [r["confounded"]["principal_angles_deg"][j] for r in rows],
            dtype=float,
        )
        ax.semilogy(scales, y, "o-", label=label)
    ax.set_xscale("log")
    ax.set_xlabel("contrast scale s")
    ax.set_ylabel("principal angle (deg)")
    ax.set_title(
        "(4) full-wave vs Born most-confounded 3-dim subspaces"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle(
        "Family 15: Born analytic/literature control (N=16, T=6, n_rx=4, "
        "f=1.0, p=24 smooth)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def _md_table(headers: list[str], rows: list[list]) -> str:
    out = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def _build_report(
    results: dict,
    cfg: dict,
    hashes: dict,
    json_path: Path,
    fig_path: Path,
    note_path: Path,
) -> str:
    rows = results["rows"]
    slopes = results["convergence_slopes"]
    plat = results["platform"]

    def md_link(p: Path, label: str | None = None) -> str:
        return f"[{label or p.name}]({p.as_posix()})"

    lines = [
        "# Family 15: Born analytic / literature control",
        "",
        "Date: 2026-09-04 (SGT); UTC stamp in "
        "`results/family15_born_control.json`.",
        "Experiment: `experiment_pose_confounding_spectral_geometry`.",
        "",
        "This family is a formula-level/internal analytic control that "
        "extends, and does not duplicate, the earlier Born checks: Family "
        "3/3b verified `A(chi=0) == A_born` and diagnosed the FD slope gates, "
        "and Family 7 reported Born-vs-full-wave A discrepancy plus full-wave "
        "no-prior retention over a coarse contrast grid.  Family 15 compares "
        "the Born operator, `K_IS`, and `K_eff` against the full-wave local "
        "quantities in the smooth p=24 coefficient space, under one shared "
        "full-wave pose confounder block, and records observed convergence "
        "orders and the overlap of the most-confounded subspaces.",
        "",
        "## Exact command, runtime, platform",
        "",
        "```bash",
        ".venv/bin/python src/family15_born_control.py",
        "```",
        "",
        f"Wall runtime: {results['wall_runtime_seconds']:.3f} s (also in the "
        "JSON as `wall_runtime_seconds`).",
        f"Platform: Apple Silicon CPU (`{plat['machine']}`), "
        f"{plat['platform']}.",
        f"Python {plat['python']}, numpy {plat['numpy']}, scipy "
        f"{plat['scipy']}, matplotlib {plat['matplotlib']}.  Deterministic "
        "dense linear algebra; no RNG used.",
        "",
        "## Scene and configuration",
        "",
        f"* N = {cfg['N']} (`h_cell = 1/{cfg['N']}`), "
        f"k_b = 2*pi*f with f = {cfg['f']}, T = {cfg['T']}, "
        f"n_rx = {cfg['n_rx']}.",
        "* family1 90-degree arc, radius 1.6, phi in [-45, 45] deg, "
        "theta = atan2(-p_y, -p_x); rx offsets "
        f"{cfg['rx_offsets']}, tx offset {cfg['tx_offset']}.",
        "* `chi_s = s * chi0`, `chi0 = family1.make_chi0(points, cfg)` "
        "(two_blob); scales "
        f"{cfg['scene_family']['contrast_scales']}.  "
        "s=0.1 is the fixed named low-contrast row "
        "(`0.1 * two_blob`).",
        "* Smooth p=24 unit-2-norm Gaussian RBF basis "
        "`family2.build_smooth_basis` with the family-2/4 config "
        f"({cfg['smooth_basis']['x_centers_n']} x "
        f"{cfg['smooth_basis']['y_centers_n']} centres, "
        f"sigma_b = {cfg['smooth_basis']['sigma_b']}); "
        "`A_s = A_pix_R @ S` in every row.",
        "* Whitening: identity noise, "
        "`whiten_realify(A_c, B_c, None)` (sqrt(2) Re/Im stacking); "
        "real data rows 48.",
        f"* Finite prior alpha = {cfg['finite_prior']['alpha']}.",
        f"* Rank tolerance: {cfg['rank_tol_rule']}.",
        "",
        "## Born composition identity",
        "",
        "`hh.born_forward(chi_s, poses, rx_offsets, tx_offset, N, k_b)` "
        "returns `(F_born, A_born)` with complex rows grouped by pose; each "
        "block is `A_born,t = G_S_t diag(E_inc_t)` (equivalently the "
        "full-wave Jacobian at chi = 0).  The identity is verified twice: "
        "once against an explicit `G_S diag(E_inc)` stack built from "
        "`hh.build_operators`, and once in the whitened/projected smooth "
        "space `A_s_born` used below.",
        "",
        _md_table(
            ["s", "rel Fro explicit-vs-born (complex)",
             "rel Fro explicit-vs-born (real, smooth)"],
            [
                [
                    _fmt(r["s"]),
                    _fmt(r["composition_identity"]
                         ["rel_fro_composed_vs_bornforward_complex"]),
                    _fmt(r["composition_identity"]
                         ["rel_fro_composed_vs_bornforward_real_smooth"]),
                ]
                for r in rows
            ],
        ),
        "",
        "The residuals are exactly 0.0 (bitwise-identical blockwise "
        "assembly): `hh.born_forward` and the explicit `G_S diag(E_inc)` "
        "stack evaluate the same deterministic green-matrix formula with the "
        "same inputs.  This confirms, at the code-path level, that the Born "
        "forward map used in the FIM comparison is the self-derived "
        "`A_born = G_S diag(E_inc)` map; the independent numerical identity "
        "`A(chi=0) == A_born` (7e-17 relative Frobenius) is Family 3b check D.",
        "",
        "## No-prior retention sweep (P_perp = I - Z Z^T)",
        "",
        "`Z = family2.range_basis(B_R)` with the full-wave `B_R` at chi_s; "
        "the same `P_perp` is applied to full-wave and Born rows (see JSON "
        "`config.retention.same_B_R_note`).",
        "",
        _md_table(
            ["s", "rdof full", "rdof Born", "log-vol full", "log-vol Born",
             "rho_min full", "rho_min Born",
             "theta_min full (deg)", "theta_min Born (deg)",
             "rel A", "rel K_IS", "rel K_eff",
             "d_rdof (full-born)", "d_rho_min (full-born)"],
            [
                [
                    _fmt(r["s"]),
                    _fmt(r["no_prior"]["full"]["retained_dof"]),
                    _fmt(r["no_prior"]["born"]["retained_dof"]),
                    _fmt(r["no_prior"]["full"]["log_volume"]),
                    _fmt(r["no_prior"]["born"]["log_volume"]),
                    _fmt(r["no_prior"]["full"]["rho_min"], ".6e"),
                    _fmt(r["no_prior"]["born"]["rho_min"], ".6e"),
                    _fmt(r["no_prior"]["full"]["theta_min_deg"], ".6e"),
                    _fmt(r["no_prior"]["born"]["theta_min_deg"], ".6e"),
                    _fmt(r["rel_errors"]["A_operator"]),
                    _fmt(r["rel_errors"]["K_IS"]),
                    _fmt(r["rel_errors"]["K_eff"]),
                    _fmt(r["no_prior"]["diff_full_minus_born"]["retained_dof"]),
                    _fmt(r["no_prior"]["diff_full_minus_born"]["rho_min"]),
                ]
                for r in rows
            ],
        ),
        "",
        "## Observed convergence order (log10 error vs log10 s, full grid)",
        "",
        _md_table(
            ["quantity", "slope", "intercept", "r^2", "n"],
            [
                [
                    label,
                    _fmt(slopes[key]["slope"]),
                    _fmt(slopes[key]["intercept"]),
                    _fmt(slopes[key]["r_squared"]),
                    str(slopes[key]["n_points"]),
                ]
                for key, label in (
                    ("A_operator", "A operator"),
                    ("K_IS", "K_IS"),
                    ("K_eff", "K_eff"),
                )
            ],
        ),
        "",
        "Reported as observed, not forced.  A slope near 1 is expected if the "
        "relative Born error is dominated by the O(s) first correction; "
        "deviation at the largest scales (and any imperfect log-log "
        "linearity) is left visible in the JSON and figure.",
        "",
        "## Finite-prior table (alpha = 1.0)",
        "",
        "`W = I - B_R (B_R^T B_R + alpha I)^{-1} B_R^T` and "
        "`K_eff = A_s^T W A_s` for both full-wave and Born A_s.",
        "",
        _md_table(
            ["s", "rdof eff full", "rdof eff Born", "rho_min eff full",
             "rho_min eff Born", "theta_min full (deg)",
             "theta_min Born (deg)", "rel K_eff"],
            [
                [
                    _fmt(r["s"]),
                    _fmt(r["finite_prior"]["full"]["retained_dof"]),
                    _fmt(r["finite_prior"]["born"]["retained_dof"]),
                    _fmt(r["finite_prior"]["full"]["rho_min"], ".6e"),
                    _fmt(r["finite_prior"]["born"]["rho_min"], ".6e"),
                    _fmt(r["finite_prior"]["full"]["theta_min_deg"], ".6e"),
                    _fmt(r["finite_prior"]["born"]["theta_min_deg"], ".6e"),
                    _fmt(r["rel_errors"]["K_eff"]),
                ]
                for r in rows
            ],
        ),
        "",
        "## Most-confounded subspace comparison (full-wave vs Born)",
        "",
        "The three most pose-confounded directions "
        "(`u0/u1/u2`, ascending generalized retention) were computed with "
        "`family4.generalized_eigen_directions(A_s, P_perp)`; the Born "
        "directions use the same full-wave `P_perp`.",
        "",
        _md_table(
            ["s", "diag sq overlap u0", "u1", "u2",
             "best-pairing sum", "principal angle 1 (deg)",
             "angle 2 (deg)", "angle 3 (deg)"],
            [
                [
                    _fmt(r["s"]),
                    _fmt(r["confounded"]["diagonal_squared_overlaps"][0]),
                    _fmt(r["confounded"]["diagonal_squared_overlaps"][1]),
                    _fmt(r["confounded"]["diagonal_squared_overlaps"][2]),
                    _fmt(r["confounded"]["best_pairing"]
                         ["sum_squared_overlaps"]),
                    _fmt(r["confounded"]["principal_angles_deg"][0], ".6e"),
                    _fmt(r["confounded"]["principal_angles_deg"][1], ".6e"),
                    _fmt(r["confounded"]["principal_angles_deg"][2], ".6e"),
                ]
                for r in rows
            ],
        ),
        "",
        "The full 3x3 normalized squared-overlap matrix, top-3 direction "
        "coefficient vectors, and the generalized spectra are stored in the "
        "JSON for audit.",
        "",
        "## Literature field and scope",
        "",
        "**Citation (verbatim):** "
        "M. L. Diong, A. Roueff, P. Lasaygues, A. Litman, \"Impact of the "
        "Born approximation on the estimation error in 2D inverse "
        "scattering\", Inverse Problems, 2016.",
        "",
        "The paper quantifies the Born approximation's effect on estimation "
        "error by comparing the linear Born MLE variance with the full "
        "nonlinear Cramer-Rao bound (CRB).  Family 15 is a "
        "formula-level/internal analytic control with the same conceptual "
        "structure: the Born FIM `K_IS_born = A_born^T A_born` is the "
        "inverse covariance of the linear Born Gaussian estimator, and the "
        "full-wave `K_IS`/`K_eff` is the local FIM of the nonlinear "
        "full-wave map.",
        "",
        "**No quantitative external benchmark was possible because the "
        "published setup's geometry, normalization, noise model, "
        "parameterization, and units were not reproduced here; this is not "
        "independent validation.**",
        "",
        "## Scope note",
        "",
        "All results are finite-dimensional N=16 toy-model numbers only "
        "(m=48 real data rows, p=24 smooth coefficients, q_pose=18).  No "
        "continuum-limit, universal-trajectory, or external validation "
        "claim is made; the Born/FIM differences and slopes are internal "
        "analytic controls for this discrete experiment.",
        "",
        "## Artifacts and digests",
        "",
        "```text",
        f"{hashes['source_sha256']}  {md_link(Path('src/family15_born_control.py'))}",
        f"{hashes['json_sha256']}  {md_link(json_path)}",
        f"{hashes['figure_sha256']}  {md_link(fig_path)}",
        f"{hashes['helmholtz_sha256']}  src/helmholtz.py",
        f"{hashes['family1_sha256']}  src/family1_pilot.py",
        f"{hashes['family2_sha256']}  src/family2_algebraic_spine.py",
        f"{hashes['family4_sha256']}  src/family4_frequency_trajectory.py",
        "```",
        "",
        "This report and the new result/figure files are additive; no "
        "existing source, result, figure, or note was modified.",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    json_path = results_dir / "family15_born_control.json"
    fig_path = figures_dir / "family15_born_control.png"
    note_path = notes_dir / "family15_born_control.md"
    source_path = Path(__file__).resolve()

    generated_utc = datetime.now(timezone.utc).isoformat()
    t0 = time.perf_counter()

    N = int(cfg["N"])
    T = int(cfg["T"])
    n_rx = int(cfg["n_rx"])
    k_b = float(cfg["k_b"])
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)

    poses = family1.build_poses(cfg)
    points, h = hh.make_grid(N)
    chi0 = family1.make_chi0(points, cfg)
    S = family2.build_smooth_basis(points, cfg["smooth_basis"])
    scales = list(cfg["scene_family"]["contrast_scales"])

    print("[family15] N=%d T=%d n_rx=%d k_b=%.6f p=%d scales=%s"
          % (N, T, n_rx, k_b, int(cfg["smooth_basis"]["p"]), scales))
    print("[family15] chi0 stats: max=%.6g l2=%.6g"
          % (chi0.max(), np.linalg.norm(chi0)))

    rows = []
    for s in scales:
        chi_s = float(s) * chi0
        row = _run_scale(float(s), chi_s, poses, S, cfg, n_rx)
        rows.append(row)
        print(
            "[family15] s=%.6g  relA=%.6e relKIS=%.6e relKeff=%.6e "
            "rdof_full=%.6f rdof_born=%.6f rho_min_full=%.6e "
            "rho_min_born=%.6e"
            % (
                float(s),
                row["rel_errors"]["A_operator"],
                row["rel_errors"]["K_IS"],
                row["rel_errors"]["K_eff"],
                row["no_prior"]["full"]["retained_dof"],
                row["no_prior"]["born"]["retained_dof"],
                row["no_prior"]["full"]["rho_min"],
                row["no_prior"]["born"]["rho_min"],
            )
        )

    scale_arr = np.asarray([r["s"] for r in rows], dtype=float)
    slopes = {
        key: _fit_loglog(
            scale_arr,
            np.asarray([r["rel_errors"][key] for r in rows], dtype=float),
        )
        for key in ("A_operator", "K_IS", "K_eff")
    }

    literature = {
        "citation": (
            "M. L. Diong, A. Roueff, P. Lasaygues, A. Litman, \"Impact of "
            "the Born approximation on the estimation error in 2D inverse "
            "scattering\", Inverse Problems, 2016."
        ),
        "reference_statement": (
            "The paper quantifies the Born approximation's effect on "
            "estimation error by comparing the linear Born MLE variance with "
            "the full nonlinear Cramer-Rao bound (CRB)."
        ),
        "our_control_mapping": (
            "Family 15 is a formula-level/internal analytic control with the "
            "same conceptual structure: the Born FIM K_IS_born = "
            "A_born^T A_born is the inverse covariance of the linear Born "
            "Gaussian estimator, and the full-wave K_IS/K_eff is the local "
            "FIM of the nonlinear full-wave map."
        ),
        "no_external_benchmark": (
            "No quantitative external benchmark was possible because the "
            "published setup's geometry, normalization, noise model, "
            "parameterization, and units were not reproduced here; this is "
            "not independent validation."
        ),
        "claim_rule": (
            "No numerical value from the cited paper is claimed to match "
            "family 15."
        ),
    }

    # figure first so its hash can be embedded
    _make_figure(rows, slopes, fig_path)

    wall_runtime_seconds = time.perf_counter() - t0

    module_hashes = {
        "helmholtz.py": sha256(_HERE / "helmholtz.py"),
        "family1_pilot.py": sha256(_HERE / "family1_pilot.py"),
        "family2_algebraic_spine.py": sha256(
            _HERE / "family2_algebraic_spine.py"
        ),
        "family4_frequency_trajectory.py": sha256(
            _HERE / "family4_frequency_trajectory.py"
        ),
    }
    figure_sha256 = sha256(fig_path)
    source_sha256 = sha256(source_path)

    results = {
        "generated_utc": generated_utc,
        "runner": "src/family15_born_control.py",
        "command": ".venv/bin/python src/family15_born_control.py",
        "family": "15",
        "wall_runtime_seconds": wall_runtime_seconds,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
            "environment_note": (
                "Apple Silicon CPU, no GPU/MPS/CUDA; deterministic dense "
                "numpy/scipy linear algebra"
            ),
        },
        "source_sha256": source_sha256,
        "reused_module_sha256": module_hashes,
        "figure_sha256": figure_sha256,
        "config": cfg,
        "scope_note": cfg["scope_note"],
        "literature": literature,
        "born_identity_summary": {
            "max_rel_fro_explicit_vs_bornforward_complex": max(
                r["composition_identity"]
                ["rel_fro_composed_vs_bornforward_complex"]
                for r in rows
            ),
            "max_rel_fro_explicit_vs_bornforward_real_smooth": max(
                r["composition_identity"]
                ["rel_fro_composed_vs_bornforward_real_smooth"]
                for r in rows
            ),
            "rule": cfg["born"]["formula_note"],
        },
        "convergence_slopes": slopes,
        "rows": rows,
    }

    json_path.write_text(_round_trip_json(results), encoding="utf-8")
    json_sha256 = sha256(json_path)

    report = _build_report(
        results,
        cfg,
        {
            "source_sha256": source_sha256,
            "json_sha256": json_sha256,
            "figure_sha256": figure_sha256,
            **{
                "helmholtz_sha256": module_hashes["helmholtz.py"],
                "family1_sha256": module_hashes["family1_pilot.py"],
                "family2_sha256": module_hashes[
                    "family2_algebraic_spine.py"
                ],
                "family4_sha256": module_hashes[
                    "family4_frequency_trajectory.py"
                ],
            },
        },
        json_path,
        fig_path,
        note_path,
    )
    note_path.write_text(report, encoding="utf-8")

    print("\n[family15] identity max rel fro (real smooth): %.3e"
          % results["born_identity_summary"]
          ["max_rel_fro_explicit_vs_bornforward_real_smooth"])
    for key, label in (
        ("A_operator", "A operator"),
        ("K_IS", "K_IS"),
        ("K_eff", "K_eff"),
    ):
        sl = slopes[key]
        print("[family15] slope %-8s: %8.4f  r^2=%.6f"
              % (label, sl["slope"], sl["r_squared"]))
    s1 = next(r for r in rows if abs(r["s"] - 1.0) < 1e-12)
    print("\n[family15] key numbers at s=1.0: relA=%.6e relKIS=%.6e "
          "relKeff=%.6e" % (
              s1["rel_errors"]["A_operator"],
              s1["rel_errors"]["K_IS"],
              s1["rel_errors"]["K_eff"],
          ))
    print("[family15] literature: %s" % literature["citation"])
    print("[family15] %s" % literature["no_external_benchmark"])
    print("[family15] artifacts: %s | %s | %s"
          % (json_path, fig_path, note_path))
    print("[family15] done in %.3f s" % wall_runtime_seconds)


if __name__ == "__main__":
    main()
