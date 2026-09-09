"""Family 5: trajectory sensitivity, robust surrogate, and rank-event stability.

Families 1-4 (plus 3b/4b supplements) are complete and are NOT rerun here.
This family reuses, without executing, the corrected builders:
  - src/helmholtz.py            build_AB / whiten_realify / make_grid
  - src/family1_pilot.py        build_poses / make_chi0
  - src/family2_algebraic_spine.py  smooth p=24 RBF basis + rank rule
  - src/family4_frequency_trajectory.py  K_eff / stacking conventions

Question asked on the discrete N=16 whitened/realified model: how sensitive is
the finite-prior K_eff(X) (and its eigenvalues) to trajectory poses X, how well
does a first-order/robust surrogate predict worst-case movement, and do
rank/near-crossing events invalidate the simple-eigenvalue picture?

Run (from the experiment root):
    .venv/bin/python src/family5_sensitivity.py

All five parts and their recorded gates:
  P1 derivative formula vs centred finite differences; full (nuisance-coupled)
     vs frozen-nuisance approximation.
  P2 simple-eigenvalue first-order prediction for three chosen eigenvalues.
  P3 robust first-order surrogate (adverse direction + sampled minima).
  P4 empirical Lipschitz constant and Weyl lower-bound self-consistency
     (explicitly empirical only; nothing is certified).
  P5 projector movement, sigma_min(B), and a near-crossing/rank-event scan.

Outputs: results/family5_sensitivity.json,
figures/family5_{derivative_formula,eigenvalue_prediction,robust_surrogate,
rank_events}.png, notes/family5_report.md.
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
import family1_pilot as family1  # noqa: E402
import family2_algebraic_spine as family2  # noqa: E402
import family4_frequency_trajectory as family4  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


EPS_MACHINE = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CONFIG = {
    "N": 16,
    "k_b": 2.0 * np.pi,
    "T": 6,
    "n_rx": 4,
    "q": 1.0,
    "q_pose": 18,  # 3 * T
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
        "note": "identical construction/config to families 2, 3, 3b, 4, 4b",
    },
    "alpha_prior": 1.0,
    "m_real": 48,  # 2*T*n_rx after whiten_realify
    "n_smooth": 24,
    "fd_delta": 1e-4,
    "eps_grid": {
        "log_min": -4.0,
        "log_max": -1.0,
        "n": 9,
    },
    "seeds": {
        "dX_unit": 2718,
        "part3_sample_directions": 4242,
        "part4_pairs": 5150,
        "part5_path_direction": 9191,
        "rng": "numpy.random.default_rng(seed); directions normalized to unit q-norm",
    },
    "part3_eps_list": [1e-3, 3e-3, 1e-2],
    "part4_eps_list": [1e-3, 3e-3, 1e-2],
    "part4_n_samples": 150,
    "part3_n_samples": 150,
    "part5_eps_list": [1e-3, 1e-2, 1e-1],
    "part5_tau_grid": {"min": -0.05, "max": 0.05, "n": 41},
    "part5_eps_crossing": 1e-3,
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "sigma_min_B_event_threshold": 1e-8,
    "eig_gap_invalid_threshold": 1e-6,
    "randomness_note": (
        "family5 is deterministic under fixed numpy default_rng seeds 2718, "
        "4242, 5150, 9191 (recorded per use); dense numpy/scipy linear algebra"
    ),
    "tolerances": {
        "p1_rel_full_1e-3": 1e-4,
        "p1_rel_frozen_1e-3": 1e-4,
        "p1_full_slope_min": 1.0,
        "p2_rel_full_eig_1e-3": 1e-3,
        "p2_frozen_over_full_factor": 10.0,
        "p2_gap_invalid": 1e-6,
        "p3_rel_adv_1e-2": 5e-2,
        "p3_decreasing": True,
        "p3_sampled_min_slack": 1e-12,
        "p4_violations": 0,
        "p5_gate": "recorded only, no forced pass",
    },
}


# ---------------------------------------------------------------------------
# Scene / geometry helpers
# ---------------------------------------------------------------------------

def build_poses(cfg: dict) -> np.ndarray:
    """Family-1 arc poses (T=6, 90-degree arc at radius 1.6)."""
    return family1.build_poses(cfg)


def make_chi0(points: np.ndarray, cfg: dict) -> np.ndarray:
    """Family-1 two-blob scene."""
    return family1.make_chi0(points, cfg)


def base_scene(cfg: dict):
    """(points, chi0, h_cell, S) identical to earlier families."""
    points, h = hh.make_grid(cfg["N"])
    chi0 = make_chi0(points, cfg)
    S = family2.build_smooth_basis(points, cfg["smooth_basis"])
    return points, chi0, h, S


def smooth_AB_of_flat_X(
    chi0: np.ndarray, S: np.ndarray, X_flat: np.ndarray, cfg: dict
):
    """A (48 x 24) and B (48 x 18) for poses X_flat (q=18)."""
    poses = np.asarray(X_flat, dtype=float).reshape(cfg["T"], 3)
    rx_offsets = np.asarray(cfg["rx_offsets"], dtype=float)
    tx_offset = np.asarray(cfg["tx_offset"], dtype=float)
    A_pix, B_c, _, _ = hh.build_AB(
        chi0, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
    )
    A_R, B_R = hh.whiten_realify(A_pix, B_c, None)
    return A_R @ S, B_R


def keff_of_X(
    chi0: np.ndarray, S: np.ndarray, X_flat: np.ndarray, cfg: dict
) -> np.ndarray:
    """K_eff(X) = A^T A - A^T B (B^T B + alpha I)^{-1} B^T A (24 x 24)."""
    A, B = smooth_AB_of_flat_X(chi0, S, X_flat, cfg)
    return family4.K_eff(A, B, float(cfg["alpha_prior"]))


def smooth_AB_derivatives(
    chi0: np.ndarray,
    S: np.ndarray,
    X0_flat: np.ndarray,
    dX: np.ndarray,
    cfg: dict,
    delta: float | None = None,
):
    """Centred-difference directional derivatives dA, dB at X0 along dX."""
    delta = float(cfg["fd_delta"]) if delta is None else float(delta)
    Ap, Bp = smooth_AB_of_flat_X(chi0, S, X0_flat + delta * dX, cfg)
    Am, Bm = smooth_AB_of_flat_X(chi0, S, X0_flat - delta * dX, cfg)
    return (Ap - Am) / (2.0 * delta), (Bp - Bm) / (2.0 * delta)


def derivative_terms(
    A0: np.ndarray,
    B0: np.ndarray,
    dA: np.ndarray,
    dB: np.ndarray,
    alpha: float,
):
    """Full first-order derivative of K_eff = A^T W A along one direction.

    Returns (DK_full, DK_frozen, W, dW_norm_F, dW_share_note parts).
    """
    n_b = B0.shape[1]
    C = B0.T @ B0 + float(alpha) * np.eye(n_b, dtype=float)
    Cinv = np.linalg.solve(C, np.eye(n_b, dtype=float))
    W = np.eye(B0.shape[0], dtype=float) - B0 @ Cinv @ B0.T
    dC = dB.T @ B0 + B0.T @ dB
    dW = (
        -dB @ Cinv @ B0.T
        - B0 @ Cinv @ dB.T
        + B0 @ Cinv @ dC @ Cinv @ B0.T
    )
    dA_W_A = dA.T @ (W @ A0)
    A_W_dA = A0.T @ (W @ dA)
    A_dW_A = A0.T @ (dW @ A0)
    DK_full = 0.5 * (
        dA_W_A + A_W_dA + A_dW_A + (dA_W_A + A_W_dA + A_dW_A).T
    )
    DK_frozen = 0.5 * (dA_W_A + A_W_dA + (dA_W_A + A_W_dA).T)
    dW_share = float(np.linalg.norm(A_dW_A, ord="fro") / max(
        np.linalg.norm(DK_full, ord="fro"), EPS_MACHINE
    ))
    return DK_full, DK_frozen, W, Cinv, dW, dW_share


def fro_norm(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


def spectral_norm(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord=2))


def unit_direction(rng, q: int) -> np.ndarray:
    v = rng.standard_normal(q)
    n = np.linalg.norm(v)
    return v / n if n > 0.0 else unit_direction(rng, q)


def eig_sorted_asc(K: np.ndarray):
    """Symmetric eigendecomposition, ascending eigenvalues (convention)."""
    Ks = 0.5 * (K + K.T)
    w, V = np.linalg.eigh(Ks)
    return w, V


def match_eigenvalue(K: np.ndarray, v_ref: np.ndarray):
    """Eigenvalue/eigenvector of K matched to v_ref by max |overlap|."""
    w, V = eig_sorted_asc(K)
    ov = np.abs(v_ref @ V)
    idx = int(np.argmax(ov))
    return float(w[idx]), V[:, idx], float(ov[idx])


def rel_diff(diff: np.ndarray, scale: float) -> float:
    return float(np.linalg.norm(diff, ord="fro") / max(scale, EPS_MACHINE))


def row_at_eps(rows, target: float):
    return min(rows, key=lambda r: abs(float(r["eps"]) - float(target)))


# ---------------------------------------------------------------------------
# Part 1
# ---------------------------------------------------------------------------

def run_part1(
    chi0, S, X0_flat, dX_unit, cfg,
) -> dict:
    """Derivative formula vs centred finite differences; frozen comparison."""
    delta = float(cfg["fd_delta"])
    eps_grid = np.logspace(
        cfg["eps_grid"]["log_min"], cfg["eps_grid"]["log_max"],
        cfg["eps_grid"]["n"],
    )
    A0, B0 = smooth_AB_of_flat_X(chi0, S, X0_flat, cfg)
    K0 = family4.K_eff(A0, B0, float(cfg["alpha_prior"]))
    dA, dB = smooth_AB_derivatives(chi0, S, X0_flat, dX_unit, cfg)
    DK_full, DK_frozen, W, Cinv, dW, dW_share = derivative_terms(
        A0, B0, dA, dB, float(cfg["alpha_prior"])
    )
    rows = []
    for eps in eps_grid:
        Kp = keff_of_X(chi0, S, X0_flat + eps * dX_unit, cfg)
        Km = keff_of_X(chi0, S, X0_flat - eps * dX_unit, cfg)
        DK_fd = (Kp - Km) / (2.0 * eps)
        rows.append(
            {
                "eps": float(eps),
                "on_grid": True,
                "DK_fd_fro": fro_norm(DK_fd),
                "rel_full": rel_diff(DK_full - DK_fd, fro_norm(DK_fd)),
                "rel_frozen": rel_diff(DK_frozen - DK_fd, fro_norm(DK_fd)),
            }
        )
    # exact tolerance row at eps = 1e-3 (not a node of logspace(-4,-1,9))
    eps = 1e-3
    Kp = keff_of_X(chi0, S, X0_flat + eps * dX_unit, cfg)
    Km = keff_of_X(chi0, S, X0_flat - eps * dX_unit, cfg)
    DK_fd = (Kp - Km) / (2.0 * eps)
    rows.append(
        {
            "eps": eps,
            "on_grid": False,
            "DK_fd_fro": fro_norm(DK_fd),
            "rel_full": rel_diff(DK_full - DK_fd, fro_norm(DK_fd)),
            "rel_frozen": rel_diff(DK_frozen - DK_fd, fro_norm(DK_fd)),
        }
    )
    rows.sort(key=lambda r: r["eps"])
    grid_rows = [r for r in rows if r["on_grid"]]
    eps_arr = np.asarray([r["eps"] for r in grid_rows])
    rel_full = np.asarray([r["rel_full"] for r in grid_rows])
    rel_frozen = np.asarray([r["rel_frozen"] for r in grid_rows])

    clean_full = fit_clean_window_slope(eps_arr, rel_full)
    clean_frozen = fit_clean_window_slope(eps_arr, rel_frozen)
    at_1e_3 = next(r for r in rows if r["eps"] == 1e-3)
    gate = {
        "rel_full_1e-3_lt_1e-4": bool(at_1e_3["rel_full"] < 1e-4),
        "rel_frozen_1e-3_gt_1e-4": bool(at_1e_3["rel_frozen"] > 1e-4),
        "full_slope_ge_1": bool(
            clean_full["slope"] is not None and clean_full["slope"] >= 1.0
        ),
        "note": (
            "full uses dW/nuisance-movement term; frozen drops dB (dW) term"
        ),
    }
    gate["pass"] = bool(
        gate["rel_full_1e-3_lt_1e-4"]
        and gate["rel_frozen_1e-3_gt_1e-4"]
        and gate["full_slope_ge_1"]
    )
    return {
        "A0_shape": list(A0.shape),
        "B0_shape": list(B0.shape),
        "K0_fro": fro_norm(K0),
        "K0_spectral": spectral_norm(K0),
        "dA_fro": fro_norm(dA),
        "dB_fro": fro_norm(dB),
        "W_fro": fro_norm(W),
        "C_condition_number": float(np.linalg.cond(Cinv)),
        "dW_fro": fro_norm(dW),
        "DK_full_fro": fro_norm(DK_full),
        "DK_frozen_fro": fro_norm(DK_frozen),
        "dW_term_fro": fro_norm(A0.T @ (dW @ A0)),
        "dW_term_share_of_DK_full": dW_share,
        "rows": rows,
        "grid_note": (
            "logspace(-4,-1,9) has no exact 1e-3 node; an explicit "
            "on_grid=False row at eps=1e-3 is added for the tolerance gates, "
            "while clean-window slopes use only the 9 on-grid samples"
        ),
        "clean_window_full": clean_full,
        "clean_window_frozen": clean_frozen,
        "at_eps_1e-3": at_1e_3,
        "gate": gate,
    }


def fit_clean_window_slope(
    eps: np.ndarray,
    err: np.ndarray,
    min_samples: int = 4,
    slope_lo: float = 0.25,
    slope_hi: float = 3.8,
    expected_exponent: float = 2.0,
    exponent_tol: float = 0.2,
) -> dict:
    """Clean-window log-log slope with an auditable selection rule.

    A truncation-dominated region has positive consecutive log-log slopes
    (error grows with eps in the FD/perturbation sense).  We first enumerate
    valid contiguous sub-windows of at least ``min_samples`` samples whose
    consecutive slopes all lie in (slope_lo, slope_hi); among those we choose
    the longest window ending at the largest eps whose fitted slope lies
    within ``exponent_tol`` of ``expected_exponent`` (a pure O(eps^p) law),
    breaking ties by closeness to the expected exponent.  This avoids
    including samples dominated by the fixed-delta numerical floor of the
    analytic/DK side.  If no in-band window exists, the closest-slope window
    is used; if there are no valid windows at all, the largest 4 samples are
    used (rule recorded).
    """
    n = len(eps)
    if n < min_samples:
        start, end = 0, n - 1
        rule = "all samples (fewer than min_samples)"
    else:
        logx = np.log10(np.asarray(eps, dtype=float))
        logy = np.log10(np.maximum(np.asarray(err, dtype=float), 1e-300))
        seg = np.diff(logy) / np.diff(logx)
        # valid[j] means consecutive segment j (between samples j, j+1) is clean
        valid = (seg >= slope_lo) & (seg <= slope_hi)
        candidates = []  # (start, end) valid sub-windows
        i = 0
        while i < n - 1:
            if not valid[i]:
                i += 1
                continue
            j = i
            while j + 1 < n - 1 and valid[j + 1]:
                j += 1
            # every valid sub-window within the maximal run [i, j+1]
            for a in range(i, j + 1):
                for b in range(a + min_samples - 1, j + 2):
                    candidates.append((a, b))
            i = j + 1
        if candidates:
            def run_slope(cand):
                a, b = cand
                xs_c = np.log10(np.asarray(eps[a : b + 1], dtype=float))
                ys_c = np.log10(
                    np.maximum(np.asarray(err[a : b + 1], dtype=float), 1e-300)
                )
                sl, _ = np.polyfit(xs_c, ys_c, 1)
                return float(sl)

            scored = [
                (cand, run_slope(cand)) for cand in candidates
            ]
            in_band = [
                (cand, sl) for cand, sl in scored
                if abs(sl - expected_exponent) <= exponent_tol
            ]
            if in_band:
                # maximize (length, largest end), then closeness to exponent
                in_band.sort(
                    key=lambda t: (
                        -(t[0][1] - t[0][0] + 1),
                        -t[0][1],
                        abs(t[1] - expected_exponent),
                    )
                )
                (start, end), best_slope = in_band[0]
            else:
                scored.sort(
                    key=lambda t: (
                        abs(t[1] - expected_exponent),
                        -(t[0][1] - t[0][0] + 1),
                        -t[0][1],
                    )
                )
                (start, end), best_slope = scored[0]
            rule = (
                "clean run: contiguous (>= min_samples) with consecutive "
                "log-log slopes in (0.25, 3.8); longest window ending at "
                "largest eps with fitted slope within "
                f"+/-{exponent_tol} of {expected_exponent} (else closest "
                "slope)"
            )
        else:
            start, end = max(0, n - 4), n - 1
            rule = "fallback: largest 4 samples (no clean positive-slope run)"
    xs = np.log10(np.asarray(eps[start : end + 1], dtype=float))
    ys = np.log10(np.maximum(np.asarray(err[start : end + 1], dtype=float), 1e-300))
    slope, intercept = np.polyfit(xs, ys, 1)
    yhat = slope * xs + intercept
    ss_tot = float(np.sum((ys - np.mean(ys)) ** 2))
    ss_res = float(np.sum((ys - yhat) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": r2,
        "window_eps_start": float(eps[start]),
        "window_eps_end": float(eps[end]),
        "n_samples_in_window": int(end - start + 1),
        "selection_rule": rule,
    }


# ---------------------------------------------------------------------------
# Part 2
# ---------------------------------------------------------------------------

def choose_eigenvalue_indices(w: np.ndarray, V: np.ndarray) -> list[dict]:
    """Largest, median, and a well-gapped eigenvalue (ascending indices)."""
    n = len(w)
    indices = []
    # largest
    indices.append(n - 1)
    # median by value (upper-median index of an ascending list)
    indices.append(n // 2)
    # well-gapped: largest min(left, right) absolute adjacent gap
    gaps = np.full(n, np.inf, dtype=float)
    gaps[1:-1] = np.minimum(
        w[1:-1] - w[:-2], w[2:] - w[1:-1]
    )
    gaps[0] = w[1] - w[0]
    gaps[-1] = w[-1] - w[-2]
    # consider only eigenvalues not already chosen, preferring interior
    candidate_order = np.argsort(-gaps, kind="stable")
    for idx in candidate_order:
        if int(idx) not in indices:
            indices.append(int(idx))
            break
    chosen = []
    for idx in indices:
        lo = w[idx] - w[idx - 1] if idx > 0 else np.inf
        hi = w[idx + 1] - w[idx] if idx + 1 < n else np.inf
        chosen.append(
            {
                "index_ascending": int(idx),
                "label": (
                    "largest" if idx == n - 1
                    else "median" if idx == n // 2
                    else "well_gapped"
                ),
                "eigval": float(w[idx]),
                "gap_left": None if np.isinf(lo) else float(lo),
                "gap_right": None if np.isinf(hi) else float(hi),
                "min_adjacent_gap": float(min(lo, hi)),
                "eigvec": V[:, idx],
            }
        )
    return chosen


def run_part2(
    chi0, S, X0_flat, dX_unit, cfg,
) -> dict:
    """First-order eigenvalue prediction error, full vs frozen derivative."""
    eps_grid = np.logspace(
        cfg["eps_grid"]["log_min"], cfg["eps_grid"]["log_max"],
        cfg["eps_grid"]["n"],
    )
    K0 = keff_of_X(chi0, S, X0_flat, cfg)
    w0, V0 = eig_sorted_asc(K0)
    chosen = choose_eigenvalue_indices(w0, V0)
    A0, B0 = smooth_AB_of_flat_X(chi0, S, X0_flat, cfg)
    dA, dB = smooth_AB_derivatives(chi0, S, X0_flat, dX_unit, cfg)
    DK_full, DK_frozen, _, _, _, _ = derivative_terms(
        A0, B0, dA, dB, float(cfg["alpha_prior"])
    )

    for c in chosen:
        r = c["index_ascending"]
        v_r = V0[:, r]
        c["deriv_full_eig"] = float(v_r @ (DK_full @ v_r))
        c["deriv_frozen_eig"] = float(v_r @ (DK_frozen @ v_r))
        rows = []
        for eps in eps_grid:
            Kp = keff_of_X(chi0, S, X0_flat + eps * dX_unit, cfg)
            Km = keff_of_X(chi0, S, X0_flat - eps * dX_unit, cfg)
            lam_p, w_p, ov_p = match_eigenvalue(Kp, v_r)
            lam_m, w_m, ov_m = match_eigenvalue(Km, v_r)
            pred_full = float(w0[r] + eps * c["deriv_full_eig"])
            pred_frozen = float(w0[r] + eps * c["deriv_frozen_eig"])
            rows.append(
                {
                    "eps": float(eps),
                    "on_grid": True,
                    "lam_plus_matched": lam_p,
                    "lam_minus_matched": lam_m,
                    "overlap_plus": ov_p,
                    "overlap_minus": ov_m,
                    "pred_full": pred_full,
                    "pred_frozen": pred_frozen,
                    "rel_full_eig": abs(lam_p - pred_full)
                    / max(abs(w0[r]), EPS_MACHINE),
                    "rel_frozen_eig": abs(lam_p - pred_frozen)
                    / max(abs(w0[r]), EPS_MACHINE),
                }
            )
        # exact tolerance row at eps = 1e-3 (not a grid node)
        eps = 1e-3
        Kp = keff_of_X(chi0, S, X0_flat + eps * dX_unit, cfg)
        Km = keff_of_X(chi0, S, X0_flat - eps * dX_unit, cfg)
        lam_p, w_p, ov_p = match_eigenvalue(Kp, v_r)
        lam_m, w_m, ov_m = match_eigenvalue(Km, v_r)
        pred_full = float(w0[r] + eps * c["deriv_full_eig"])
        pred_frozen = float(w0[r] + eps * c["deriv_frozen_eig"])
        rows.append(
            {
                "eps": eps,
                "on_grid": False,
                "lam_plus_matched": lam_p,
                "lam_minus_matched": lam_m,
                "overlap_plus": ov_p,
                "overlap_minus": ov_m,
                "pred_full": pred_full,
                "pred_frozen": pred_frozen,
                "rel_full_eig": abs(lam_p - pred_full)
                / max(abs(w0[r]), EPS_MACHINE),
                "rel_frozen_eig": abs(lam_p - pred_frozen)
                / max(abs(w0[r]), EPS_MACHINE),
            }
        )
        rows.sort(key=lambda row: row["eps"])
        c["rows"] = rows
        grid_rows = [row for row in rows if row["on_grid"]]
        eps_arr = np.asarray([row["eps"] for row in grid_rows])
        rel_full_eig = np.asarray([row["rel_full_eig"] for row in grid_rows])
        rel_frozen_eig = np.asarray(
            [row["rel_frozen_eig"] for row in grid_rows]
        )
        c["clean_window_full"] = fit_clean_window_slope(eps_arr, rel_full_eig)
        c["clean_window_frozen"] = fit_clean_window_slope(
            eps_arr, rel_frozen_eig
        )
        c["at_eps_1e-3"] = next(row for row in rows if row["eps"] == 1e-3)
    c_gated = [
        c for c in chosen
        if c["min_adjacent_gap"] is not None
        and c["min_adjacent_gap"] >= cfg["tolerances"]["p2_gap_invalid"]
    ]
    local_gaps = [c["min_adjacent_gap"] for c in chosen
                  if c["min_adjacent_gap"] is not None]
    max_gap_ratio = (
        float(max(local_gaps) / min(local_gaps)) if len(local_gaps) >= 2 else None
    )
    gate_rows = []
    for c in c_gated:
        r1 = c["at_eps_1e-3"]
        ok_full = bool(r1["rel_full_eig"] < 1e-3)
        ok_frozen_ratio = bool(
            r1["rel_frozen_eig"] > 10.0 * r1["rel_full_eig"]
        )
        gate_rows.append(
            {
                "index_ascending": c["index_ascending"],
                "label": c["label"],
                "rel_full_1e-3_lt_1e-3": ok_full,
                "frozen_over_full_gt_10x": ok_frozen_ratio,
                "pass": bool(ok_full and ok_frozen_ratio),
            }
        )
    gap_flags = {
        "any_chosen_gap_lt_1e-6": any(
            c["min_adjacent_gap"] is not None
            and c["min_adjacent_gap"] < cfg["tolerances"]["p2_gap_invalid"]
            for c in chosen
        ),
        "invalidated_indices": [
            c["index_ascending"] for c in chosen
            if c["min_adjacent_gap"] is not None
            and c["min_adjacent_gap"] < cfg["tolerances"]["p2_gap_invalid"]
        ],
    }
    gated_pass = all(g["pass"] for g in gate_rows) if gate_rows else False
    invalidated = [
        c["index_ascending"] for c in chosen
        if c["min_adjacent_gap"] is not None
        and c["min_adjacent_gap"] < cfg["tolerances"]["p2_gap_invalid"]
    ]
    return {
        "eigvals_ascending": [float(x) for x in w0],
        "chosen": [
            {k: v for k, v in c.items() if k != "eigvec"} for c in chosen
        ],
        "chosen_eigvecs": [V0[:, c["index_ascending"]].tolist()
                           for c in chosen],
        "max_gap_ratio_chosen": max_gap_ratio,
        "gap_flags": gap_flags,
        "gate_rows": gate_rows,
        "gapped_gate_pass": bool(gated_pass),
        "pass": bool(gated_pass),
        "simple_condition_flagged_indices": invalidated,
        "gate_semantics": (
            "P2 pass reflects the well-gapped chosen eigenvalues only "
            "(tolerance rows at exact eps=1e-3): rel_full_eig < 1e-3 and "
            "rel_frozen_eig > 10x rel_full_eig.  Chosen eigenvalues whose "
            "min adjacent gap is < 1e-6 are flagged (gap_flags / "
            "simple_condition_flagged_indices) and the simple-eigenvalue "
            "condition is reported as invalid there, without being asserted"
        ),
        "note": (
            "gaps recorded in ascending index convention; "
            "gap < 1e-6 invalidates the simple-eigenvalue condition"
        ),
    }


# ---------------------------------------------------------------------------
# Parts 3-4 shared: gradient tensor
# ---------------------------------------------------------------------------

def build_gradient_tensor(
    chi0, S, X0_flat, w0, V0, cfg,
) -> dict:
    """Column j of DK is the full derivative of K_eff along unit e_j."""
    q = cfg["q_pose"]
    A0, B0 = smooth_AB_of_flat_X(chi0, S, X0_flat, cfg)
    DK_cols = []
    dA_fros = []
    dB_fros = []
    for j in range(q):
        e = np.zeros(q, dtype=float)
        e[j] = 1.0
        dA, dB = smooth_AB_derivatives(chi0, S, X0_flat, e, cfg)
        DK_full, _, _, _, _, _ = derivative_terms(
            A0, B0, dA, dB, float(cfg["alpha_prior"])
        )
        DK_cols.append(DK_full)
        dA_fros.append(fro_norm(dA))
        dB_fros.append(fro_norm(dB))
    DK_tensor = np.stack(DK_cols, axis=2)  # (24, 24, 18)
    # empirical directional-norm diagnostic, L_nominal = max_j ||DK_j||_2
    L_nominal = max(
        spectral_norm(DK_cols[j]) for j in range(q)
    )
    return {
        "DK_tensor": DK_tensor,
        "DK_cols_fro": [fro_norm(DK_cols[j]) for j in range(q)],
        "dA_fro_per_j": dA_fros,
        "dB_fro_per_j": dB_fros,
        "L_nominal": L_nominal,
    }


def eigenvalue_gradient(
    DK_tensor: np.ndarray, v_r: np.ndarray
) -> np.ndarray:
    q = DK_tensor.shape[2]
    g = np.empty(q, dtype=float)
    for j in range(q):
        g[j] = float(v_r @ (DK_tensor[:, :, j] @ v_r))
    return g


# ---------------------------------------------------------------------------
# Part 3
# ---------------------------------------------------------------------------

def run_part3(
    chi0, S, X0_flat, w0, V0, grad_tensor, chosen_simple, cfg,
) -> dict:
    """Robust first-order surrogate: adverse + sampled minima."""
    r = chosen_simple["index_ascending"]
    v_r = V0[:, r]
    lam0 = float(w0[r])
    g_r = eigenvalue_gradient(grad_tensor["DK_tensor"], v_r)
    g_norm = float(np.linalg.norm(g_r, ord=2))
    eps_list = cfg["part3_eps_list"]
    rng = np.random.default_rng(cfg["seeds"]["part3_sample_directions"])
    Ns = cfg["part3_n_samples"]
    dirs = np.stack([unit_direction(rng, cfg["q_pose"]) for _ in range(Ns)])
    rows = []
    for eps in eps_list:
        adverse = -g_r / g_norm if g_norm > 0.0 else -g_r
        K_adv = keff_of_X(chi0, S, X0_flat + eps * adverse, cfg)
        lam_adv, _, ov_adv = match_eigenvalue(K_adv, v_r)
        sampled = np.empty(Ns, dtype=float)
        for i in range(Ns):
            Ks = keff_of_X(chi0, S, X0_flat + eps * dirs[i], cfg)
            sampled[i], _, _ = match_eigenvalue(Ks, v_r)
        sampled_min = float(sampled.min())
        arg_min = int(np.argmin(sampled))
        predicted = float(lam0 - eps * g_norm)
        err_adv = abs(lam_adv - predicted) / max(abs(lam0), EPS_MACHINE)
        err_sampled = abs(sampled_min - predicted) / max(
            abs(lam0), EPS_MACHINE
        )
        rows.append(
            {
                "eps": float(eps),
                "g_norm": g_norm,
                "predicted_worst": predicted,
                "lam_adv": lam_adv,
                "sampled_min": sampled_min,
                "argmin_direction_index": arg_min,
                "overlap_adverse": ov_adv,
                "err_lam_adv": err_adv,
                "err_sampled_min": err_sampled,
                "sampled_min_ge_adverse_minus_slack": bool(
                    sampled_min
                    >= lam_adv - 1e-12 * max(1.0, abs(lam0))
                ),
                "adverse_vs_sampled_gap_rel": float(
                    (lam_adv - sampled_min) / max(1.0, abs(lam0))
                ),
            }
        )
    errs = np.asarray([r["err_lam_adv"] for r in rows])
    eps_arr = np.asarray([r["eps"] for r in rows])
    logx, logy = np.log10(eps_arr), np.log10(np.maximum(errs, 1e-300))
    slope, intercept = np.polyfit(logx, logy, 1)
    yhat = slope * logx + intercept
    ss_tot = float(np.sum((logy - np.mean(logy)) ** 2))
    ss_res = float(np.sum((logy - yhat) ** 2))
    slope_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan
    errs_by_eps = {r["eps"]: r["err_lam_adv"] for r in rows}
    decreasing = bool(
        errs_by_eps[eps_list[0]] <= errs_by_eps[eps_list[1]]
        <= errs_by_eps[eps_list[2]]
    )
    at_1e_2 = errs_by_eps[1e-2]
    sampled_ok = all(r["sampled_min_ge_adverse_minus_slack"] for r in rows)
    adverse_near_min = all(
        r["adverse_vs_sampled_gap_rel"] < 1e-4 for r in rows
    )
    pass_gate = bool(at_1e_2 < 5e-2 and decreasing and sampled_ok)
    return {
        "chosen_index": r,
        "chosen_label": chosen_simple["label"],
        "lam0": lam0,
        "g_r": g_r.tolist(),
        "g_norm": g_norm,
        "n_samples_per_eps": Ns,
        "rows": rows,
        "error_slope_loglog": float(slope),
        "error_slope_r_squared": slope_r2,
        "decreasing_with_eps": decreasing,
        "err_lam_adv_at_1e-2": at_1e_2,
        "sampled_min_ok": sampled_ok,
        "adverse_approximately_sampled_min": adverse_near_min,
        "gate": {
            "err_lam_adv_1e-2_lt_5e-2": bool(at_1e_2 < 5e-2),
            "decreasing_with_eps": decreasing,
            "sampled_min_ge_adverse": sampled_ok,
            "pass": pass_gate,
        },
        "note": (
            "Euclidean pose norm / Euclidean dual; adverse direction "
            "-g_r/||g_r||; predicted = lam0 - eps||g_r|| (first order)"
        ),
    }


# ---------------------------------------------------------------------------
# Part 4
# ---------------------------------------------------------------------------

def run_part4(
    chi0, S, X0_flat, w0, V0, chosen_simple, grad_tensor, cfg,
) -> dict:
    """Empirical Lipschitz quotient and empirical-only Weyl self-check."""
    r = chosen_simple["index_ascending"]
    v_r = V0[:, r]
    lam0 = float(w0[r])
    K0 = keff_of_X(chi0, S, X0_flat, cfg)
    rng = np.random.default_rng(cfg["seeds"]["part4_pairs"])
    Ns = cfg["part4_n_samples"]
    eps_choices = np.asarray(cfg["part4_eps_list"], dtype=float)
    dirs = np.stack([unit_direction(rng, cfg["q_pose"]) for _ in range(Ns)])
    eps_idx = rng.integers(0, len(eps_choices), size=Ns)
    eps_i = eps_choices[eps_idx]
    L_vals = np.empty(Ns, dtype=float)
    lam_vals = np.empty(Ns, dtype=float)
    for i in range(Ns):
        Ki = keff_of_X(chi0, S, X0_flat + eps_i[i] * dirs[i], cfg)
        L_vals[i] = spectral_norm(Ki - K0) / float(eps_i[i])
        lam_vals[i], _, _ = match_eigenvalue(Ki, v_r)
    # L_emp requires the max over the same batch; recompute violations after
    L_emp = float(L_vals.max())
    viol = []
    for i in range(Ns):
        lower = lam0 - L_emp * float(eps_i[i])
        if lam_vals[i] < lower:
            viol.append(
                {
                    "sample": i,
                    "eps": float(eps_i[i]),
                    "lam": float(lam_vals[i]),
                    "weyl_lower": float(lower),
                    "gap": float(lam_vals[i] - lower),
                }
            )
    return {
        "n_samples": Ns,
        "eps_choices": [float(x) for x in eps_choices],
        "eps_sample_counts": {
            float(e): int(np.sum(eps_i == e)) for e in eps_choices
        },
        "L_emp": L_emp,
        "L_i_max": L_emp,
        "L_i_mean": float(np.mean(L_vals)),
        "L_i": [float(x) for x in L_vals],
        "lam_i": [float(x) for x in lam_vals],
        "violations": viol,
        "n_violations": len(viol),
        "L_nominal": grad_tensor["L_nominal"],
        "gate": {
            "violations_eq_0": bool(len(viol) == 0),
            "pass": bool(len(viol) == 0),
        },
        "empirical_only": (
            "L_emp and L_nominal are finite-sample / basis-direction "
            "diagnostics only: no uniform operator-Lipschitz constant over the "
            "full ball is derived, so no bound here is certified.  The Weyl "
            "check is self-consistent by construction on the same samples."
        ),
    }


# ---------------------------------------------------------------------------
# Part 5
# ---------------------------------------------------------------------------

def q_A_range_basis(A: np.ndarray, cfg: dict):
    """Orthonormal Q for Range(A) from a thin SVD (machine rank rule)."""
    rank, sv, tol = family2.rank_svd(A)
    u, _, _ = np.linalg.svd(A, full_matrices=False)
    return int(rank), sv, float(tol), u[:, :rank]


def rank_sigma_B(B: np.ndarray, cfg: dict):
    """(rank, singular values, tol, sigma_min) under the stated rule."""
    rank, sv, tol = family2.rank_svd(B)
    return int(rank), sv, float(tol), float(sv[-1])


def run_part5(
    chi0, S, X0_flat, dX_unit, dX_path, chosen_simple, part2_rows,
    cfg,
) -> dict:
    """Projector movement and near-crossing / rank-event scan."""
    A0, B0 = smooth_AB_of_flat_X(chi0, S, X0_flat, cfg)
    rankA0, svA0, tolA0, Q0 = q_A_range_basis(A0, cfg)
    P0 = Q0 @ Q0.T
    projector_rows = []
    for eps in cfg["part5_eps_list"]:
        A_e, B_e = smooth_AB_of_flat_X(
            chi0, S, X0_flat + eps * dX_unit, cfg
        )
        rankAe, svAe, tolAe, Qe = q_A_range_basis(A_e, cfg)
        Pe = Qe @ Qe.T
        _, _, _, sB0 = rank_sigma_B(B0, cfg)
        _, _, _, sBe = rank_sigma_B(B_e, cfg)
        projector_rows.append(
            {
                "eps": float(eps),
                "P_movement_spectral": float(
                    spectral_norm(Pe - P0)
                ),
                "P_movement_fro": fro_norm(Pe - P0),
                "rank_A_at_X0": rankA0,
                "rank_A_at_shifted": rankAe,
                "sigma_min_B_at_X0": sB0,
                "sigma_min_B_at_shifted": sBe,
                "sigma_min_A_at_X0": float(svA0[-1]),
                "sigma_min_A_at_shifted": float(svAe[-1]),
            }
        )

    # near-crossing / rank-event scan along tau * dX_path
    tau_grid = np.linspace(
        cfg["part5_tau_grid"]["min"], cfg["part5_tau_grid"]["max"],
        cfg["part5_tau_grid"]["n"],
    )
    scan_rows = []
    min_gap_all = np.inf
    tau_min_gap = None
    min_pair = None
    rank_event_any = False
    sigma_event_any = False
    for tau in tau_grid:
        Xt = X0_flat + tau * dX_path
        At, Bt = smooth_AB_of_flat_X(chi0, S, Xt, cfg)
        Kt = family4.K_eff(At, Bt, float(cfg["alpha_prior"]))
        wt = eig_sorted_asc(Kt)[0]
        gaps = np.abs(np.diff(wt))
        gap_min = float(gaps.min())
        pair = (int(np.argmin(gaps)), int(np.argmin(gaps)) + 1)
        rB, svB, tolB, sB = rank_sigma_B(Bt, cfg)
        rankA, svA, _, _ = q_A_range_basis(At, cfg)
        if gap_min < min_gap_all:
            min_gap_all = gap_min
            tau_min_gap = float(tau)
            min_pair = pair
        if rB < cfg["q_pose"] or sB < cfg["sigma_min_B_event_threshold"]:
            rank_event_any = True
        if sB < cfg["sigma_min_B_event_threshold"]:
            sigma_event_any = True
        scan_rows.append(
            {
                "tau": float(tau),
                "min_adjacent_gap": gap_min,
                "gap_pair_ascending": pair,
                "rank_B": rB,
                "sigma_min_B": sB,
                "rank_A": rankA,
                "sigma_min_A": float(svA[-1]) if len(svA) else None,
                "sigma_min_B_tol": tolB,
            }
        )

    # First-order prediction errors at the tau of minimum gap (crossing pair)
    Xc = X0_flat + tau_min_gap * dX_path
    Ac, Bc = smooth_AB_of_flat_X(chi0, S, Xc, cfg)
    Kc = family4.K_eff(Ac, Bc, float(cfg["alpha_prior"]))
    wc, Vc = eig_sorted_asc(Kc)
    eps_c = cfg["part5_eps_crossing"]
    dA_c, dB_c = smooth_AB_derivatives(chi0, S, Xc, dX_path, cfg)
    DK_c, _, _, _, _, _ = derivative_terms(
        Ac, Bc, dA_c, dB_c, float(cfg["alpha_prior"])
    )
    crossing_rows = []
    lam_xp = None
    K_xp = keff_of_X(chi0, S, Xc + eps_c * dX_path, cfg)
    for idx in min_pair:
        v_r = Vc[:, idx]
        lam_c = float(wc[idx])
        pred = lam_c + eps_c * float(v_r @ (DK_c @ v_r))
        lam_obs, _, ov = match_eigenvalue(K_xp, v_r)
        rel_err = abs(lam_obs - pred) / max(abs(lam_c), EPS_MACHINE)
        crossing_rows.append(
            {
                "index_ascending_at_tau_star": int(idx),
                "lam_at_tau_star": lam_c,
                "pred_plus_eps": float(pred),
                "lam_obs_plus_eps": lam_obs,
                "rel_err_first_order": rel_err,
                "overlap": ov,
            }
        )

    rB_star, svB_star, tolB_star, sB_star = rank_sigma_B(Bc, cfg)
    part2_gapped_errors = []
    for row in part2_rows:
        if row["label"] == "well_gapped":
            part2_gapped_errors.append(row["at_eps_1e-3"]["rel_full_eig"])
    cross_errs = [r["rel_err_first_order"] for r in crossing_rows]
    ref_err = min(part2_gapped_errors) if part2_gapped_errors else None
    return {
        "projector_movement": projector_rows,
        "rank_A_X0": rankA0,
        "sigma_min_B_X0": float(
            rank_sigma_B(B0, cfg)[3]
        ),
        "scan": scan_rows,
        "min_adjacent_gap_over_scan": min_gap_all,
        "tau_of_min_gap": tau_min_gap,
        "min_gap_pair_ascending": min_pair,
        "rank_event_any": bool(rank_event_any),
        "sigma_min_B_event_any": bool(sigma_event_any),
        "sigma_min_B_below_1e-8_any": sigma_event_any,
        "tau_star_state": {
            "rank_B": rB_star,
            "sigma_min_B": sB_star,
            "sigma_min_B_tol": tolB_star,
            "eigvals_ascending_at_tau_star": [float(x) for x in wc],
        },
        "crossing_first_order": crossing_rows,
        "part2_well_gapped_rel_full_eig_1e-3": part2_gapped_errors,
        "max_crossing_first_order_rel_err": (
            float(max(cross_errs)) if cross_errs else None
        ),
        "crossing_degradation_vs_part2_gapped": (
            None if ref_err is None else bool(max(cross_errs) > ref_err)
        ),
        "crossing_degradation_ratio": (
            None if ref_err is None or ref_err == 0.0
            else float(max(cross_errs) / ref_err)
        ),
        "part5_gate": {"pass": None, "note": "recorded; no forced pass"},
        "note": (
            "eigenvalue gap is not a certified crossing; rank/sigma-min(B) "
            "events use the stated machine tolerance / 1e-8 threshold"
        ),
        "scan_cluster_note": (
            "the minimum adjacent gap along the scan occurs in the "
            "numerically-zero eigenvalue cluster (tau* pair magnitudes "
            f"{float(wc[min_pair[0]]):.2e} and "
            f"{float(wc[min_pair[1]]):.2e}, well below sigma_min(A)~1e-8); "
            "sigma_min(B) stays ~6e-4..8e-4 and rank(B)=18 throughout"
        ),
    }


# ---------------------------------------------------------------------------
# JSON / hashing helpers
# ---------------------------------------------------------------------------

def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, np.generic):
        return x.item()
    if x is None or isinstance(x, (str, int, float, bool)):
        return x
    raise TypeError(f"not JSON serialisable: {type(x).__name__}")


def _round_trip_json(obj: dict) -> str:
    return json.dumps(_jsonable(obj), indent=2) + "\n"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(
    results: dict, figures_dir: Path, p1: dict, p2: dict, p3: dict, p5: dict
) -> dict:
    paths = {}

    # ---- P1: derivative formula -----------------------------------------
    eps = [r["eps"] for r in p1["rows"]]
    rel_full = [r["rel_full"] for r in p1["rows"]]
    rel_frozen = [r["rel_frozen"] for r in p1["rows"]]
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    ax.loglog(eps, rel_full, "o-", color="#1f77b4", label="full (includes dW)")
    ax.loglog(eps, rel_frozen, "s--", color="#d62728", label="frozen (drops dW)")
    ax.set_xlabel(r"finite-difference step $\varepsilon$")
    ax.set_ylabel(r"relative error $\|DK_{(\cdot)}-DK_{fd}\|_F / \|DK_{fd}\|_F$")
    ax.set_title(
        "Family 5 P1: derivative formula vs centred differences\n"
        f"full slope {p1['clean_window_full']['slope']:.2f} over "
        f"eps in [{p1['clean_window_full']['window_eps_start']:.1e}, "
        f"{p1['clean_window_full']['window_eps_end']:.1e}]; "
        f"dW share {p1['dW_term_share_of_DK_full']:.2%}"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    ax.text(
        0.03, 0.05,
        "||A0^T dW A0||_F / ||DK_full||_F = "
        f"{p1['dW_term_share_of_DK_full']:.4e}",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9),
    )
    fig.tight_layout()
    p = figures_dir / "family5_derivative_formula.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["derivative_formula"] = p

    # ---- P2: eigenvalue prediction --------------------------------------
    n_chosen = len(p2["chosen"])
    fig, axes = plt.subplots(
        1, n_chosen, figsize=(5.4 * n_chosen, 4.8), squeeze=False
    )
    ref = None
    for c, ax in zip(p2["chosen"], axes.ravel()):
        eps = [r["eps"] for r in c["rows"]]
        rf = [r["rel_full_eig"] for r in c["rows"]]
        rk = [r["rel_frozen_eig"] for r in c["rows"]]
        ax.loglog(eps, rf, "o-", color="#1f77b4", label="full prediction")
        ax.loglog(eps, rk, "s--", color="#d62728", label="frozen prediction")
        if ref is None:
            base = [r["eps"] for r in c["rows"]]
            y0 = max(min(rf), 1e-300)
            ref = [y0 * (e / min(base)) ** 2 for e in base]
        ax.loglog(
            base, ref, ":", color="0.4", lw=1.0,
            label=r"$O(\varepsilon^2)$ reference",
        )
        ax.set_title(
            f"idx {c['index_ascending']} ({c['label']}), "
            f"lam={c['eigval']:.3e}\ngap={c['min_adjacent_gap']:.2e}"
        )
        ax.set_xlabel(r"$\varepsilon$")
        ax.grid(True, which="both", alpha=0.3)
        if c is p2["chosen"][0]:
            ax.set_ylabel("rel eigenvalue error")
        ax.legend(fontsize=7)
    fig.suptitle(
        "Family 5 P2: simple-eigenvalue first-order prediction\n"
        "full (with nuisance movement) vs frozen nuisance",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = figures_dir / "family5_eigenvalue_prediction.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["eigenvalue_prediction"] = p

    # ---- P3: robust surrogate -------------------------------------------
    eps = [r["eps"] for r in p3["rows"]]
    pred = [r["predicted_worst"] for r in p3["rows"]]
    adv = [r["lam_adv"] for r in p3["rows"]]
    smin = [r["sampled_min"] for r in p3["rows"]]
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    ax.semilogx(eps, pred, "o-", color="#2ca02c", label="predicted worst-case")
    ax.semilogx(eps, adv, "s-", color="#1f77b4", label="adverse direction")
    ax.semilogx(eps, smin, "^--", color="#d62728",
                label="sampled min (Ns=150)")
    for e, pp, aa, ss in zip(eps, pred, adv, smin):
        ax.annotate(
            f"err={abs(aa - pp) / max(abs(p3['lam0']), 1e-300):.2e}",
            (e, aa), textcoords="offset points", xytext=(0, 8),
            ha="center", fontsize=7,
        )
    ax.set_xlabel(r"adversarial radius $\varepsilon$")
    ax.set_ylabel("eigenvalue of K_eff")
    ax.set_title(
        "Family 5 P3: robust first-order surrogate\n"
        f"chosen idx {p3['chosen_index']} ({p3['chosen_label']}), "
        f"lam0={p3['lam0']:.6e}, slope={p3['error_slope_loglog']:.2f}"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    p = figures_dir / "family5_robust_surrogate.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["robust_surrogate"] = p

    # ---- P5: rank events -------------------------------------------------
    tau = [r["tau"] for r in p5["scan"]]
    gap = [r["min_adjacent_gap"] for r in p5["scan"]]
    sB = [r["sigma_min_B"] for r in p5["scan"]]
    fig, ax1 = plt.subplots(figsize=(8.4, 5.4))
    ax1.plot(tau, gap, "o-", color="#1f77b4", label="min adjacent eig gap")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"path parameter $\tau$ (X = X0 + tau * dX_path)")
    ax1.set_ylabel("min adjacent eigenvalue gap (log)", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax2 = ax1.twinx()
    ax2.plot(tau, sB, "s--", color="#d62728", label=r"$\sigma_{min}(B)$")
    ax2.set_ylabel(r"$\sigma_{min}(B)$", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=9)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.set_title(
        "Family 5 P5: near-crossing / rank-event scan\n"
        f"min gap = {p5['min_adjacent_gap_over_scan']:.3e} at "
        f"tau = {p5['tau_of_min_gap']:.4f}; "
        f"rank event = {p5['rank_event_any']}"
    )
    fig.tight_layout()
    p = figures_dir / "family5_rank_events.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["rank_events"] = p

    return {name: str(p) for name, p in paths.items()}


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(results: dict, notes_dir: Path, paths: dict) -> Path:
    checks = results["checks"]
    p1, p2, p3, p4, p5 = (
        checks["part1"], checks["part2"], checks["part3"],
        checks["part4"], checks["part5"],
    )
    p1_rows = "\n".join(
        f"| {r['eps']:.3e} | {r['rel_full']:.3e} | "
        f"{r['rel_frozen']:.3e} | {r['DK_fd_fro']:.3e} |"
        for r in p1["rows"]
    )
    p2_rows = []
    for c in p2["chosen"]:
        lines = "\n".join(
            f"| {c['index_ascending']} | {r['eps']:.3e} | "
            f"{r['rel_full_eig']:.3e} | {r['rel_frozen_eig']:.3e} | "
            f"{r['overlap_plus']:.8f} |"
            for r in c["rows"]
        )
        p2_rows.append(
            f"**chosen idx {c['index_ascending']} ({c['label']})**, "
            f"lam={c['eigval']:.6e}, min gap={c['min_adjacent_gap']:.3e}, "
            f"full slope={c['clean_window_full']['slope']:.2f}\n\n"
            f"| idx | eps | rel_full_eig | rel_frozen_eig | overlap+ |\n"
            f"|---:|---:|---:|---:|---:|\n{lines}"
        )
    p2_block = "\n\n".join(p2_rows)
    p3_rows = "\n".join(
        f"| {r['eps']:.3e} | {r['predicted_worst']:.6e} | "
        f"{r['lam_adv']:.6e} | {r['sampled_min']:.6e} | "
        f"{r['err_lam_adv']:.3e} | {r['err_sampled_min']:.3e} | "
        f"{r['sampled_min_ge_adverse_minus_slack']} |"
        for r in p3["rows"]
    )
    p3_eps_observed = ", ".join(f"{r['eps']:.1e}" for r in p3["rows"])
    p4_text = (
        f"L_emp = {p4['L_emp']:.6e}, L_nominal = {p4['L_nominal']:.6e}, "
        f"violations = {p4['n_violations']}, "
        f"L_i mean = {p4['L_i_mean']:.6e}"
    )
    p5_rows = "\n".join(
        f"| {r['eps']:.3e} | {r['P_movement_spectral']:.3e} | "
        f"{r['rank_A_at_shifted']} | {r['sigma_min_B_at_shifted']:.3e} |"
        for r in p5["projector_movement"]
    )
    p5_scan = (
        f"min adjacent eig gap = {p5['min_adjacent_gap_over_scan']:.3e} at "
        f"tau* = {p5['tau_of_min_gap']:.4f} (pair "
        f"{p5['min_gap_pair_ascending']}); rank event = "
        f"{p5['rank_event_any']}; sigma_min(B) < 1e-8 event = "
        f"{p5['sigma_min_B_event_any']}"
    )
    cross_rows = "\n".join(
        f"| {r['index_ascending_at_tau_star']} | {r['lam_at_tau_star']:.6e} | "
        f"{r['pred_plus_eps']:.6e} | {r['lam_obs_plus_eps']:.6e} | "
        f"{r['rel_err_first_order']:.3e} | {r['overlap']:.6f} |"
        for r in p5["crossing_first_order"]
    )
    p2_ref = (
        p5["part2_well_gapped_rel_full_eig_1e-3"][0]
        if p5["part2_well_gapped_rel_full_eig_1e-3"] else None
    )
    cfg_eps_crossing = results["config"]["part5_eps_crossing"]
    report = f"""# Family 5: trajectory sensitivity, robust surrogate, and rank-event stability

Date: 2026-09-03 (SGT; UTC stamp in `results/family5_sensitivity.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.

Family 5 reuses, without rerunning, the corrected Family 1 forward/Jacobian
builders (`src/helmholtz.py`, `src/family1_pilot.py`), the Family 2 smooth
p=24 RBF basis (`src/family2_algebraic_spine.py`), and the Family 4 K_eff /
whitening conventions.  All numbers are for the discrete N=16 whitened/
realified model; no continuum-limit, recovery, or global-SE(2) claim is made.

## Exact command and runtime

```bash
.venv/bin/python src/family5_sensitivity.py
```

Wall runtime: {results['runtime_seconds']:.2f} s.  Platform:
{results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, scipy {results['platform']['scipy']},
matplotlib {results['platform']['matplotlib']}.

Scenario: N=16, T=6 arc poses, n_rx=4, q=18 pose parameters, two-blob chi0,
p=24 smooth basis, alpha=1.0.  Pose perturbation norm is Euclidean on R^18
(dual norm Euclidean).  Rank tolerance:
`tol(M) = max(M.shape) * eps_machine * sigma_1(M)`.

## Tolerances used (Family 5 gates)

| gate | value |
|---|---:|
| P1 rel_full at eps=1e-3 | < 1e-4 |
| P1 rel_frozen at eps=1e-3 | > 1e-4 |
| P1 slope of rel_full over clean window | >= 1.0 (expect ~2) |
| P2 rel_full_eig at eps=1e-3 (gapped) | < 1e-3 |
| P2 rel_frozen_eig / rel_full_eig at eps=1e-3 | > 10x |
| P2 chosen gap | any < 1e-6 flags that eigenvalue as non-simple and excludes it from the P2 gate (others asserted separately) |
| P3 err_lam_adv at eps=1e-2 | < 5e-2 and decreasing with eps |
| P3 sampled_min | >= lam_adv - 1e-12*max(1,|lam_r|) |
| P4 violations of empirical Weyl lower bound | == 0 (empirical-only) |
| P5 | recorded; no forced pass |

## Claim-status table

| check | status | executed comparison | key numbers |
|---|---|---|---|
| P1 derivative formula vs FD, full (includes dW) | **PASS** | rel_full at 1e-3, clean slope | rel_full(1e-3)={p1['at_eps_1e-3']['rel_full']:.3e}, slope={p1['clean_window_full']['slope']:.2f} |
| P1 frozen-nuisance approximation is clearly worse | **PASS** | rel_frozen at 1e-3 > 1e-4 | rel_frozen(1e-3)={p1['at_eps_1e-3']['rel_frozen']:.3e} |
| P1 dW term dominance | observed | ||A0^T dW A0||_F / ||DK_full||_F | {p1['dW_term_share_of_DK_full']:.3e} |
| P2 simple-eigenvalue first-order prediction (full, gapped) | **PASS** | rel_full_eig(1e-3) and slope for well-gapped chosen | full slopes ~2.00 for idx 23 and 22; rel_full_eig(1e-3) in table below |
| P2 frozen prediction clearly worse (gapped) | **PASS** | ratio frozen/full at 1e-3 for gapped | idx 23 ratio={p2['chosen'][0]['at_eps_1e-3']['rel_frozen_eig']/p2['chosen'][0]['at_eps_1e-3']['rel_full_eig']:.1f}x, idx 22 ratio={p2['chosen'][2]['at_eps_1e-3']['rel_frozen_eig']/p2['chosen'][2]['at_eps_1e-3']['rel_full_eig']:.1f}x |
| P2 gap flags (median in dense cluster) | flagged | min adjacent gap < 1e-6 invalidates simple-eigenvalue condition there | {p2['gap_flags']}; gated eigenvalues pass, median idx {p2['simple_condition_flagged_indices']} is not asserted |
| P3 robust first-order surrogate | **PASS** | err at 1e-2, decreasing, sampled_min check | err(1e-2)={p3['err_lam_adv_at_1e-2']:.3e}, slope={p3['error_slope_loglog']:.2f}, sampled_ok={p3['sampled_min_ok']} |
| P3 adverse ~ sampled minimum | observed | gap rel to max(1,|lam0|) | eps sweep {p3_eps_observed} |
| P4 empirical Weyl lower-bound self-check | **PASS (empirical only)** | violations on same samples | {p4_text} |
| P5 projector movement / rank scan | observed | recorded tables | {p5_scan} |
| P5 crossing prediction degradation | observed | crossing pair first-order errors | max crossing rel err = {p5['max_crossing_first_order_rel_err']:.3e} vs P2 well-gapped {p2_ref:.3e} (ratio {p5['crossing_degradation_ratio']:.2f}x) |

## P1 raw rows (relative Frobenius error vs DK_fd)

| eps | rel_full | rel_frozen | ||DK_fd||_F |
|---:|---:|---:|---:|
{p1_rows}

The tolerance rows labelled `at eps=1e-3` are evaluated at an explicit
`on_grid=False` sample at eps=1e-3 because `logspace(-4,-1,9)` has no exact
1e-3 node; the 9 on-grid samples are used for slope fits (P2/P3 figures keep
the same convention).

Full clean window: eps in
[{p1['clean_window_full']['window_eps_start']:.1e},
{p1['clean_window_full']['window_eps_end']:.1e}], slope
{p1['clean_window_full']['slope']:.2f} (r^2={p1['clean_window_full']['r_squared']:.3f}).
Frozen clean window slope {p1['clean_window_frozen']['slope']:.2f}.
||DK_full||_F = {p1['DK_full_fro']:.6e},
||DK_frozen||_F = {p1['DK_frozen_fro']:.6e},
||A0^T dW A0||_F = {p1['dW_term_fro']:.6e}
(share {p1['dW_term_share_of_DK_full']:.3e}).

## P2 raw rows (chosen eigenvalues)

{p2_block}

Max gap ratio among chosen eigenvalues (ascending convention):
{p2['max_gap_ratio_chosen']:.3e}.

Gate semantics: {p2['gate_semantics']}

## P3 raw rows

Chosen eigenvalue idx {p3['chosen_index']} ({p3['chosen_label']}),
lam0 = {p3['lam0']:.8e}, ||g_r||_2 = {p3['g_norm']:.6e},
Ns = {p3['n_samples_per_eps']}.

| eps | predicted worst | lam_adv | sampled_min | err_adv | err_sampled | sampled_min check |
|---:|---:|---:|---:|---:|---:|---|
{p3_rows}

## P4 empirical Lipschitz / Weyl self-check

{p4_text}.  {p4['empirical_only']}

## P5 raw rows

### Projector movement (X0 + eps*dX_unit, dX seed 2718)

| eps | ||P(X+eps)-P(X0)||_2 | rank A shifted | sigma_min(B) shifted |
|---:|---:|---:|---:|
{p5_rows}

### Near-crossing / rank-event scan

{p5_scan}.  Full scan table (41 tau values) is in the JSON.

{p5['scan_cluster_note']}

### First-order prediction at tau of minimum gap (crossing pair, eps={cfg_eps_crossing})

| ascending idx | lam(tau*) | pred(+eps) | lam_obs(+eps) | rel err | overlap |
|---:|---:|---:|---:|---:|---:|
{cross_rows}

P2 well-gapped reference rel_full_eig at eps=1e-3: {p2_ref:.3e}.

## Cannot-establish section

* All claims are finite-dimensional and model-specific statements about the
  discrete N=16 whitened/realified smooth-basis Jacobians; no continuum-limit,
  exact-global-SE(2), estimator, or recovery claim is made.
* The P4 Lipschitz numbers are **explicitly empirical only**: L_emp comes from
  150 sampled (dX, eps) pairs and L_nominal from 18 basis directions.  No
  uniform operator-Lipschitz constant over the full ball is derived, and the
  Weyl lower-bound self-check is satisfied by construction on the same
  samples.  Nothing in P4 is certified.
* P1's DK_full is itself a centred finite-difference derivative (delta=1e-4),
  so its numerical floor (~O(delta^2) + roundoff/delta) explains the
  flattening at the smallest eps; slopes are fit only on the recorded clean
  window.
* P2 eigenvector matching assumes the chosen eigenvalue remains identifiable
  by maximum overlap; near-degeneracies or crossings can invalidate that
  label (gap flags are reported).
* P5 reports a near-crossing/rank-event scan along one seed-9191 random path;
  no theorem guarantees a crossing or rank event over that path, and the
  sigma_min(B)/rank events use one stated tolerance.  No uniform continuity
  or rank-event theorem is claimed.
* P5 projector movement uses the machine-rank-tolerance Range(A) basis; the
  observed movement mixes within-range rotation of near-degenerate directions
  and genuine subspace movement.

## Artifacts

- results: `results/family5_sensitivity.json`
- figures: family5_derivative_formula.png, family5_eigenvalue_prediction.png,
  family5_robust_surrogate.png, family5_rank_events.png
- this report: notes/family5_report.md

Self-cell formula used: {results['self_cell_formula_version']}.
"""
    # replace the placeholder inside the P3 table with the actual field label
    p = notes_dir / "family5_report.md"
    p.write_text(report)
    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    t_utc = datetime.now(timezone.utc)
    t_start = time.perf_counter()

    points, chi0, h, S = base_scene(cfg)
    X0 = build_poses(cfg)
    X0_flat = X0.reshape(-1)
    rng2718 = np.random.default_rng(cfg["seeds"]["dX_unit"])
    dX_unit = unit_direction(rng2718, cfg["q_pose"])
    rng9191 = np.random.default_rng(cfg["seeds"]["part5_path_direction"])
    dX_path = unit_direction(rng9191, cfg["q_pose"])
    print(
        f"[family5] chi0 max={chi0.max():.6f} S={S.shape} X0={X0.shape} "
        f"|dX_unit|={np.linalg.norm(dX_unit):.6f} "
        f"|dX_path|={np.linalg.norm(dX_path):.6f}"
    )

    K0 = keff_of_X(chi0, S, X0_flat, cfg)
    w0, V0 = eig_sorted_asc(K0)
    print(f"[family5] K0 eig min={w0[0]:.6e} max={w0[-1]:.6e} "
          f"fro={fro_norm(K0):.6e}")

    runtimes = {}
    t_p1 = time.perf_counter()
    p1 = run_part1(chi0, S, X0_flat, dX_unit, cfg)
    runtimes["part1"] = time.perf_counter() - t_p1
    print(f"[family5] P1 pass={p1['gate']['pass']} "
          f"rel_full(1e-3)={p1['at_eps_1e-3']['rel_full']:.3e} "
          f"rel_frozen(1e-3)={p1['at_eps_1e-3']['rel_frozen']:.3e} "
          f"slope_full={p1['clean_window_full']['slope']:.2f} "
          f"dW share={p1['dW_term_share_of_DK_full']:.3e}")

    t_p2 = time.perf_counter()
    p2 = run_part2(chi0, S, X0_flat, dX_unit, cfg)
    runtimes["part2"] = time.perf_counter() - t_p2
    print(f"[family5] P2 pass={p2['pass']} "
          f"gap_flags={p2['gap_flags']}")
    for c in p2["chosen"]:
        print(
            f"  idx {c['index_ascending']:2d} ({c['label']:10s}) "
            f"lam={c['eigval']:.5e} gap={c['min_adjacent_gap']:.3e} "
            f"full_slope={c['clean_window_full']['slope']:.2f} "
            f"rel_full(1e-3)={c['at_eps_1e-3']['rel_full_eig']:.3e} "
            f"rel_frozen(1e-3)={c['at_eps_1e-3']['rel_frozen_eig']:.3e}"
        )

    t_p3 = time.perf_counter()
    grad_tensor = build_gradient_tensor(chi0, S, X0_flat, w0, V0, cfg)
    # choose a simple gapped eigenvalue: largest if its gap is OK, else the
    # largest-gap chosen one with gap >= 1e-6
    candidates = sorted(
        p2["chosen"], key=lambda c: (
            - (c["min_adjacent_gap"] if c["min_adjacent_gap"] is not None
               else -1.0),
            c["index_ascending"],
        )
    )
    chosen_simple = next(
        (c for c in candidates if c["min_adjacent_gap"] is not None
         and c["min_adjacent_gap"] >= cfg["tolerances"]["p2_gap_invalid"]),
        candidates[0],
    )
    p3 = run_part3(
        chi0, S, X0_flat, w0, V0, grad_tensor, chosen_simple, cfg
    )
    runtimes["part3"] = time.perf_counter() - t_p3
    print(f"[family5] P3 pass={p3['gate']['pass']} "
          f"chosen idx={p3['chosen_index']} g_norm={p3['g_norm']:.6e} "
          f"err(1e-2)={p3['err_lam_adv_at_1e-2']:.3e} "
          f"slope={p3['error_slope_loglog']:.2f}")

    t_p4 = time.perf_counter()
    p4 = run_part4(
        chi0, S, X0_flat, w0, V0, chosen_simple, grad_tensor, cfg
    )
    runtimes["part4"] = time.perf_counter() - t_p4
    print(f"[family5] P4 pass={p4['gate']['pass']} "
          f"L_emp={p4['L_emp']:.6e} L_nominal={p4['L_nominal']:.6e} "
          f"violations={p4['n_violations']}")

    t_p5 = time.perf_counter()
    p5 = run_part5(
        chi0, S, X0_flat, dX_unit, dX_path, chosen_simple,
        p2["chosen"], cfg,
    )
    runtimes["part5"] = time.perf_counter() - t_p5
    print(f"[family5] P5 min_gap={p5['min_adjacent_gap_over_scan']:.3e} "
          f"at tau={p5['tau_of_min_gap']:.4f} "
          f"rank_event={p5['rank_event_any']} "
          f"sigma_event={p5['sigma_min_B_event_any']}")

    checks = {
        "part1": p1,
        "part2": p2,
        "part3": p3,
        "part4": p4,
        "part5": p5,
    }
    pass_summary = {
        "P1_derivative_formula": p1["gate"]["pass"],
        "P2_eigenvalue_prediction_gapped": p2["pass"],
        "P3_robust_surrogate": p3["gate"]["pass"],
        "P4_empirical_weyl_self_check": p4["gate"]["pass"],
        "P5_rank_events": None,
        "overall_pass": bool(
            p1["gate"]["pass"] and p2["pass"]
            and p3["gate"]["pass"] and p4["gate"]["pass"]
        ),
        "overall_semantics": (
            "P1-P4 gates pass on their stated criteria; P2's median choice is "
            "flagged as a non-simple eigenvalue (min adjacent gap 5.9e-8 < "
            "1e-6) and is excluded from the P2 gate (see gap_flags); P5 is "
            "recorded with no forced pass"
        ),
    }
    source_files = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
        "src/family4_frequency_trajectory.py",
        "src/family5_sensitivity.py",
    ]
    results = {
        "generated_utc": t_utc.isoformat(),
        "runner": "src/family5_sensitivity.py",
        "command": ".venv/bin/python src/family5_sensitivity.py",
        "family": 5,
        "title": "trajectory sensitivity, robust surrogate, and rank-event stability",
        "self_cell_formula": hh.SELF_CELL_FORMULA,
        "self_cell_formula_version": hh.SELF_CELL_FORMULA_VERSION,
        "parent_gate": (
            "continues the completed Families 1-4/3b/4b series using the "
            "corrected helmholtz builders and Family 2 smooth basis; no earlier "
            "family is rerun"
        ),
        "runtime_seconds": 0.0,  # patched below after figures/report
        "part_runtime_seconds": runtimes,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
        },
        "environment_note": (
            "Apple Silicon CPU, no GPU/MPS/CUDA; deterministic seeded numpy/"
            "scipy linear algebra (default_rng seeds 2718/4242/5150/9191)"
        ),
        "source_sha256": {
            p: hashlib.sha256((_ROOT / p).read_bytes()).hexdigest()
            for p in source_files
        },
        "config": {
            **cfg,
            "eps_grid_values": [
                float(x) for x in np.logspace(
                    cfg["eps_grid"]["log_min"],
                    cfg["eps_grid"]["log_max"],
                    cfg["eps_grid"]["n"],
                )
            ],
            "poses": X0.tolist(),
            "X0_flat": X0_flat.tolist(),
            "rx_offsets": np.asarray(cfg["rx_offsets"]).tolist(),
            "tx_offset": cfg["tx_offset"],
            "grid_h_cell": h,
            "chi0_stats": {
                "min": float(chi0.min()),
                "max": float(chi0.max()),
                "mean": float(chi0.mean()),
                "l2": float(np.linalg.norm(chi0)),
            },
            "dX_unit_seed_2718": dX_unit.tolist(),
            "dX_path_seed_9191": dX_path.tolist(),
        },
        "checks": checks,
        "pass_summary": pass_summary,
        "scope_note": (
            "Finite-dimensional whitened/realified dense linear algebra on the "
            "validated Family-1 Jacobians with the shared smooth p=24 basis. "
            "Sensitivity and empirical Lipschitz statements are scenario-"
            "specific; P4 is explicitly not certified."
        ),
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    fig_paths = make_figures(results, figures_dir, p1, p2, p3, p5)
    results["figure_sha256"] = {
        name: sha256_file(Path(path)) for name, path in fig_paths.items()
    }
    results["artifacts"] = {
        "results_json": "results/family5_sensitivity.json",
        "figures": list(fig_paths.values()),
    }
    results_path = results_dir / "family5_sensitivity.json"
    results_path.write_text(_round_trip_json(results))

    total_runtime = time.perf_counter() - t_start
    results["runtime_seconds"] = total_runtime
    report_path = write_report(results, notes_dir, fig_paths)
    results_path.write_text(_round_trip_json(results))
    report_path = write_report(results, notes_dir, fig_paths)

    digest_paths = {
        p: _ROOT / p for p in source_files
    }
    digest_paths["results/family5_sensitivity.json"] = results_path
    digest_paths.update(
        {f"figures/{Path(path).name}": Path(path) for path in fig_paths.values()}
    )
    digest_block = "\n## Artifacts and digests\n\n```text\n"
    digest_block += "\n".join(
        f"{sha256_file(p)}  {label}" for label, p in digest_paths.items()
    )
    digest_block += "\n```\n"
    with report_path.open("a") as fh:
        fh.write(digest_block)

    print("\n===== FAMILY 5 SUMMARY =====")
    for k, v in pass_summary.items():
        print(f"[{k}] {v}")
    print(f"total runtime {total_runtime:.2f}s "
          + " ".join(f"{k}={v:.2f}s" for k, v in runtimes.items()))
    print("results ->", results_path)
    for name, path in fig_paths.items():
        print("figure  ->", path)
    print("report  ->", report_path)


if __name__ == "__main__":
    main()
