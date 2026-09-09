"""Family 2: algebraic spine of the pose-confounding spectral geometry.

Run (from the experiment root):
    .venv/bin/python src/family2_algebraic_spine.py

This family is purely finite-dimensional linear algebra in the
whitened/realified data space produced by the corrected Helmholtz harness.
No wave physics is recomputed beyond `build_AB`, which is already validated by
Family 1.  The scenario (N, k_b, T, n_rx, arcs, chi0 blobs, pose convention) is
copied verbatim from `src/family1_pilot.py` (`build_poses`, `make_chi0`).

Two map parameterizations are analysed with the same realified pose Jacobian
B_R:

(A) pixel basis  A_R = A_pix_R  (48 x 256, degenerate, Range(A) = full data
    space up to machine rank);
(B) smooth basis A_smooth_R = A_pix_R @ S with S (N^2 x 24) unit-2-norm
    Gaussian RBF columns on a 4 x 6 centre lattice over D = [-0.5, 0.5]^2,
    sigma_b = 0.16 (non-degenerate, 24 < 48).

Checks implemented: kernel identity (c1), rank identity (c2), retention
spectrum and its predicted 1 - cos^2 form (c3), loss-rank bound (c4), PSD
ordering/interlacing (c5), regular-prior sweep (c6), singular-prior gauge
limit (c7), and c-scaled duplicate invariance (c8).  Smooth-basis resolution
diagnostics rerun c1-c3 at N in {32, 40}.

Important numerical decisions (documented in the result JSON and report):

* Retention operators R_op and R_eff are evaluated as
  Q_A^T W Q_A (W = P_perp or I - B(B^T B + alpha I)^-1 B^T).  This is exactly
  diag(1/s_A) V_A^T K W-active V_A diag(1/s_A) (K_SLAM or K_eff), because
  A V_A diag(1/s_A) = Q_A, but it avoids amplifying the near-null singular
  directions of the pixel and smooth Jacobians (sigma_min(A) ~ 1e-8) by
  diag(1/s_A)^2.

* rank(L_X), L_X = K_IS - K_eff(alpha) = A^T B C^-1 B^T A, is counted from
  X = C^-1/2 B^T A (rank(L_X) = rank(X)); the direct difference
  K_IS - K_eff(alpha) carries an absolute roundoff floor ~ eps ||K_IS|| that
  inflates its SVD rank when sigma_1(L_X) is small (large alpha).  Both the
  stable factor rank and the direct-difference rank are recorded.

* All ranks use the stated machine tolerance
  tol(M) = max(M.shape) * eps_machine * sigma_1(M).
"""

from __future__ import annotations

import hashlib
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

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration (scenario copied verbatim from src/family1_pilot.py)
# ---------------------------------------------------------------------------

CONFIG = {
    "N": 16,
    "k_b": 2.0 * np.pi,
    "T": 6,
    "n_rx": 4,
    "q": 1.0,
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
    },
    "m": 48,        # 2*T*n_rx after realification
    "n_pixel": 256,  # N^2
    "q_pose": 18,    # 3*T
    "c5_alphas": [1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4],
    "c4_alphas": [1e-6, 1e-3, 1.0, 1e2],
    "c6_alpha_grid": {
        "log_min": -6.0,
        "log_max": 4.0,
        "n": 25,
    },
    "c6_retention_alphas": [1e-3, 1.0, 1e2],
    "c8_duplicate_c": 2.0,
    "c8_fixed_prior_alphas": [1e-3, 1e-1, 1.0, 1e1, 1e2],
    "c8_detail_alpha": 1.0,
    "resolution_Ns": [32, 40],
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "c1_pass_scale_rule": "1e-8 * max(1, ||A||_2^2)",
    "c3_pass_scale_rule": "1e-10 * max(1, ||A||_2^2)",
    "c5_delta_rule": "1e-10 * max(1, ||K_IS||_2)",
    "seeds": [],
    "randomness_note": "deterministic dense linear algebra; no RNG used",
}


def build_poses(cfg: dict) -> np.ndarray:
    """Poses: 90-degree arc at radius 1.6, theta pointing to the origin."""
    T = cfg["T"]
    phi = np.deg2rad(np.linspace(cfg["arc_phi_deg"][0], cfg["arc_phi_deg"][1], T))
    p = cfg["arc_radius"] * np.column_stack([np.cos(phi), np.sin(phi)])
    theta = np.arctan2(-p[:, 1], -p[:, 0])
    return np.column_stack([p, theta])


def make_chi0(points: np.ndarray, cfg: dict) -> np.ndarray:
    blobs = cfg["chi_blobs"]
    c1 = np.asarray(blobs["c1"], dtype=float)
    c2 = np.asarray(blobs["c2"], dtype=float)
    chi = blobs["amp1"] * np.exp(
        -np.sum((points - c1) ** 2, axis=1) / (2.0 * blobs["sigma1"] ** 2)
    ) + blobs["amp2"] * np.exp(
        -np.sum((points - c2) ** 2, axis=1) / (2.0 * blobs["sigma2"] ** 2)
    )
    return chi.astype(float)


def build_smooth_basis(points: np.ndarray, sb_cfg: dict) -> np.ndarray:
    """Unit-2-norm Gaussian RBF columns S (n_points x p)."""
    xs = np.linspace(sb_cfg["x_span"][0], sb_cfg["x_span"][1], sb_cfg["x_centers_n"])
    ys = np.linspace(sb_cfg["y_span"][0], sb_cfg["y_span"][1], sb_cfg["y_centers_n"])
    centres = np.array([[x, y] for x in xs for y in ys], dtype=float)
    sigma_b = float(sb_cfg["sigma_b"])
    S = np.exp(
        -np.sum((points[:, None, :] - centres[None, :, :]) ** 2, axis=2)
        / (2.0 * sigma_b**2)
    )
    norms = np.linalg.norm(S, axis=0)
    S = S / norms[None, :]
    return S


def _round_trip_json(obj):
    return json.dumps(obj, indent=2, default=lambda o: o.tolist())


# ---------------------------------------------------------------------------
# Linear-algebra helpers
# ---------------------------------------------------------------------------

_EPS = np.finfo(float).eps


def rank_svd(M: np.ndarray) -> tuple[int, np.ndarray, float]:
    """Machine-precision rank via singular values (stated tolerance rule)."""
    sv = np.linalg.svd(M, compute_uv=False)
    tol = max(M.shape) * _EPS * sv[0]
    return int(np.sum(sv > tol)), sv, float(tol)


def thin_decomposition(A: np.ndarray) -> dict:
    """Thin SVD A = Q_A diag(s_A) V_A^T truncated to rank(A)."""
    r, sv, tol = rank_svd(A)
    u, s, vh = np.linalg.svd(A, full_matrices=False)
    return {
        "rank": r,
        "rank_tol": tol,
        "singular_values": sv,
        "Q": u[:, :r],
        "s": s[:r],
        "V": vh[:r].T,
    }


def range_basis(B: np.ndarray) -> tuple[int, np.ndarray, np.ndarray, float]:
    """Orthonormal Z for Range(B), rank by the stated tolerance."""
    r, sv, tol = rank_svd(B)
    u, _, _ = np.linalg.svd(B, full_matrices=False)
    return r, u[:, :r], sv, tol


def retention_spectrum(Q: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Eigenvalues (descending) of Q^T W Q (W symmetric)."""
    M = Q.T @ W @ Q
    M = 0.5 * (M + M.T)          # symmetrise the computed product
    return np.sort(np.linalg.eigvalsh(M))[::-1]


def weighted_shrinkage_operator(B: np.ndarray, alpha: float) -> np.ndarray:
    """W = I - B (B^T B + alpha I)^{-1} B^T (q x q solve, alpha > 0)."""
    q = B.shape[1]
    C = B.T @ B + alpha * np.eye(q, dtype=float)
    return np.eye(B.shape[0], dtype=float) - B @ np.linalg.solve(C, B.T)


def keff_of_w(A: np.ndarray, W: np.ndarray) -> np.ndarray:
    """K_eff = A^T W A for a symmetric data-space W."""
    K = A.T @ (W @ A)
    return 0.5 * (K + K.T)


def loss_rank_factors(
    A: np.ndarray, B: np.ndarray, alpha: float
) -> tuple[int, np.ndarray, int, np.ndarray]:
    """Rank information for L_X = A^T B (B^T B + alpha I)^{-1} B^T A.

    Returns (factor_rank, factor_sv, direct_rank, direct_sv).  factor_rank
    counts X = C^{-1/2} B^T A via its SVD; rank(L_X) = rank(X) in exact
    arithmetic.  direct_rank counts the SVD rank of K_IS - K_eff(alpha),
    which carries an absolute cancellation floor ~ eps ||K_IS||.
    """
    q = B.shape[1]
    C = B.T @ B + alpha * np.eye(q, dtype=float)
    wC, VC = np.linalg.eigh(C)
    X = VC @ np.diag(1.0 / np.sqrt(wC)) @ VC.T @ (B.T @ A)
    factor_r, factor_sv, _ = rank_svd(X)
    KIS = A.T @ A
    K_eff = KIS - A.T @ B @ np.linalg.solve(C, B.T @ A)
    direct_r, direct_sv, _ = rank_svd(KIS - K_eff)
    return factor_r, factor_sv, direct_r, direct_sv


def keff_direct(A: np.ndarray, B: np.ndarray, alpha: float) -> np.ndarray:
    """K_eff(alpha) = A^T A - A^T B (B^T B + alpha I)^{-1} B^T A."""
    q = B.shape[1]
    C = B.T @ B + alpha * np.eye(q, dtype=float)
    return A.T @ A - A.T @ B @ np.linalg.solve(C, B.T @ A)


def _desc(a: np.ndarray) -> np.ndarray:
    return np.sort(np.asarray(a, dtype=float))[::-1]


def _asc(a: np.ndarray) -> np.ndarray:
    return np.sort(np.asarray(a, dtype=float))


# ---------------------------------------------------------------------------
# Per-case core checks
# ---------------------------------------------------------------------------

def run_algebraic_case(
    name: str,
    A: np.ndarray,
    B: np.ndarray,
    sb_cfg: dict | None,
) -> dict:
    """Run checks c1-c8 for one map parameterization at N = CONFIG['N']."""
    m, n = A.shape
    q = B.shape[1]
    t0 = time.perf_counter()

    da = thin_decomposition(A)
    rA = da["rank"]
    QA, sA, VA = da["Q"], da["s"], da["V"]
    rB, Z, svB, tolB = range_basis(B)
    sigma1A = float(np.linalg.norm(A, ord=2))
    back_scale = max(1.0, sigma1A**2)          # max(1, ||A||_2^2)
    c1_tol = 1e-8 * back_scale
    c3_tol = 1e-10 * back_scale

    I_m = np.eye(m, dtype=float)
    Pperp = I_m - Z @ Z.T
    KIS = A.T @ A
    KSL = A.T @ (Pperp @ A)
    C = Pperp @ A

    # ---------------- c1 kernel identity ---------------------------------
    w_eig, V_eig = np.linalg.eigh(KSL)
    _, svKSL, tolKSL = rank_svd(KSL)
    k1 = int(np.sum(np.abs(w_eig) <= tolKSL))
    N1 = V_eig[:, :k1]
    rC, svC, tolC = rank_svd(C)
    uC, sC, vhC = np.linalg.svd(C, full_matrices=True)
    N2 = vhC[rC:].T
    if N1.shape[1] and N2.shape[1]:
        subspace_distance = float(
            np.linalg.norm(N1 @ N1.T - N2 @ N2.T, ord=2)
        )
        min_sing_N1N2 = float(
            np.linalg.svd(N1.T @ N2, compute_uv=False)[-1]
        )
    else:
        subspace_distance = 0.0
        min_sing_N1N2 = None
    max_KSL_on_N2 = float(max(
        (np.linalg.norm(KSL @ N2[:, j], ord=2) for j in range(N2.shape[1])),
        default=0.0,
    ))
    max_C_on_N1 = float(max(
        (np.linalg.norm(C @ N1[:, j], ord=2) for j in range(N1.shape[1])),
        default=0.0,
    ))
    c1 = {
        "nullity_KSL_eig": int(k1),
        "rank_KSL": int(np.sum(svKSL > tolKSL)),
        "rank_C": rC,
        "nullity_C": int(C.shape[1] - rC),
        "tol_KSL_eig_selector": tolKSL,
        "tol_C": tolC,
        "eig_KSL_min_abs_values": [float(x) for x in np.abs(w_eig[:8])],
        "subspace_distance": subspace_distance,
        "min_singular_N1T_N2": min_sing_N1N2,
        "max_norm_KSL_u_over_N2": max_KSL_on_N2,
        "max_norm_C_u_over_N1": max_C_on_N1,
        "pass_gate": bool(subspace_distance < c1_tol),
        "pass_gate_tol": c1_tol,
        "identity_support": bool(
            (max_KSL_on_N2 <= 1e-10)
            and (N1.shape[1] == N2.shape[1])
            and ((min_sing_N1N2 is None) or (min_sing_N1N2 > 1 - 1e-6))
        ),
        "note": (
            "N1: eig(K_SLAM) with |eval| <= tol; N2: null(P_perp A) from "
            "full SVD.  Empty bases give distance 0 and vacuous min singular. "
            "identity_support uses null residuals and dimension equality, "
            "which are insensitive to basis mixing inside the ill-conditioned "
            "zero cluster of K_SLAM."
        ),
    }

    # ---------------- c2 rank identity ------------------------------------
    rKIS, svKIS, tolKIS = rank_svd(KIS)
    rKS, svKS2, tolKS2 = rank_svd(KSL)
    AB = np.hstack([A, B])
    rAB, svAB, tolAB = rank_svd(AB)
    d_inter = rA + rB - rAB
    lhs = rKIS - rKS
    residual = lhs - d_inter
    c2 = {
        "r_A": rA,
        "r_KIS": rKIS,
        "r_KS": rKS,
        "r_B": rB,
        "r_AB": rAB,
        "d_inter": int(d_inter),
        "lhs_rKIS_minus_rKS": int(lhs),
        "residual_lhs_minus_d_inter": int(residual),
        "lower_bound_d_inter": int(max(0, rA + rB - m)),
        "lower_bound_holds": bool(d_inter >= max(0, rA + rB - m)),
        "pass_gate": bool(int(residual) == 0),
        "rank_tols": {
            "A": da["rank_tol"],
            "B": tolB,
            "KIS": tolKIS,
            "KSL": tolKS2,
            "AB": tolAB,
        },
    }

    # ---------------- c3 retention spectrum -------------------------------
    Rop = retention_spectrum(QA, Pperp)
    c2_vals = np.linalg.svd(Z.T @ QA, compute_uv=False) ** 2
    rho_pred_vals = _desc(np.concatenate([1.0 - c2_vals, np.ones(max(rA - len(c2_vals), 0))]))
    max_abs_diff = float(np.max(np.abs(Rop - rho_pred_vals)))
    six_smallest = [float(x) for x in _asc(Rop)[:6]]
    gaps = np.diff(_asc(Rop))
    gap_struct = {
        "largest_gap": float(np.max(gaps)) if len(gaps) else None,
        "index_of_largest_gap_ascending": int(np.argmax(gaps)) if len(gaps) else None,
        "num_gaps_gt_0_5": int(np.sum(gaps > 0.5)),
        "ascending_gaps_top5": [float(x) for x in _desc(gaps)[:5]],
    }
    count_less_one = int(np.sum(Rop < 1.0 - 1e-8))
    c3 = {
        "r_A": rA,
        "length_rho": int(len(Rop)),
        "rho_desc": [float(x) for x in Rop],
        "rho_pred_desc": [float(x) for x in rho_pred_vals],
        "max_abs_diff_rho_vs_pred": max_abs_diff,
        "pass_pred_gate": bool(max_abs_diff < c3_tol),
        "pred_gate_tol": c3_tol,
        "rho_min": float(Rop.min()),
        "rho_max": float(Rop.max()),
        "bounds_tol": 1e-10,
        "bounds_ok": bool(Rop.min() >= -1e-10 and Rop.max() <= 1.0 + 1e-10),
        "count_less_than_one_minus_1e-8": count_less_one,
        "count_le_rB": bool(count_less_one <= rB),
        "r_B": rB,
        "c2_squared_cos2_desc": [float(x) for x in _desc(c2_vals)],
        "six_smallest_rho_asc": six_smallest,
        "spectral_gap_structure": gap_struct,
        "note": (
            "R_op evaluated as Q_A^T P_perp Q_A, algebraically equal to "
            "diag(1/s_A) V_A^T K_SLAM V_A diag(1/s_A) because A V_A "
            "diag(1/s_A) = Q_A; the factored form avoids 1/s_A^2 "
            "amplification in the near-null directions."
        ),
    }

    # ---------------- c4 loss rank ----------------------------------------
    c4_rows = []
    for alpha in CONFIG["c4_alphas"]:
        f_r, f_sv, d_r, d_sv = loss_rank_factors(A, B, alpha)
        c4_rows.append(
            {
                "alpha": alpha,
                "factor_rank_LX": f_r,
                "factor_rank_sv_top": [float(x) for x in f_sv[:6]],
                "factor_rank_sv_tail_min": float(f_sv[-1]),
                "direct_difference_rank": d_r,
                "direct_difference_sv_top": [float(x) for x in d_sv[:6]],
            }
        )
    c4 = {
        "rows": c4_rows,
        "r_B": rB,
        "max_factor_rank_observed": int(max(r["factor_rank_LX"] for r in c4_rows)),
        "pass_gate": bool(max(r["factor_rank_LX"] for r in c4_rows) <= rB),
        "note": (
            "rank(L_X) counted from X = C^{-1/2} B^T A (stable, equal in "
            "exact arithmetic). direct_difference_rank is the literal SVD "
            "rank of K_IS - K_eff(alpha); cancellation roundoff inflates it "
            "for large alpha."
        ),
    }

    # ---------------- c5 ordering and interlacing --------------------------
    delta_c5 = 1e-10 * max(1.0, float(np.linalg.norm(KIS, ord=2)))
    evIS_desc = _desc(np.linalg.eigvalsh(KIS))
    c5_rows = []
    c5_all_ok = True
    for alpha in CONFIG["c5_alphas"]:
        K_eff = keff_direct(A, B, alpha)
        p, _, d_r, _ = loss_rank_factors(A, B, alpha)
        evK = _desc(np.linalg.eigvalsh(K_eff))
        min_eig_dn = float(np.linalg.eigvalsh(K_eff - KSL)[0])
        min_eig_up = float(np.linalg.eigvalsh(KIS - K_eff)[0])
        violations = []
        for i in range(rA):                       # 1-based i = i+1
            lo = evIS_desc[i + p] if i + p < len(evIS_desc) else -math.inf
            violations.append(float(lo - evK[i]))
            violations.append(float(evK[i] - evIS_desc[i]))
        max_pos_viol = float(max(violations))
        ok = (
            min_eig_dn >= -delta_c5
            and min_eig_up >= -delta_c5
            and max_pos_viol <= delta_c5
        )
        c5_all_ok &= ok
        c5_rows.append(
            {
                "alpha": alpha,
                "min_eig_Keff_minus_KSLAM": min_eig_dn,
                "min_eig_KIS_minus_Keff": min_eig_up,
                "p_rank_LX": p,
                "max_positive_interlace_violation": max_pos_viol,
                "ok": bool(ok),
            }
        )
    c5 = {
        "delta": delta_c5,
        "rows": c5_rows,
        "pass_gate": bool(c5_all_ok),
        "note": (
            "K_eff computed as A^T A - A^T B solve(B^T B + alpha I, B^T A); "
            "interlacing p uses the stable factor rank of L_X."
        ),
    }

    # ---------------- c6 prior sweep ---------------------------------------
    ag = CONFIG["c6_alpha_grid"]
    alphas = np.logspace(ag["log_min"], ag["log_max"], ag["n"])
    d0_list, dinf_list = [], []
    KIS_F = float(np.linalg.norm(KIS, ord="fro"))
    for alpha in alphas:
        K_eff = keff_direct(A, B, float(alpha))
        d0_list.append(float(np.linalg.norm(K_eff - KSL, ord="fro") / KIS_F))
        dinf_list.append(float(np.linalg.norm(K_eff - KIS, ord="fro") / KIS_F))
    rho_no_prior = Rop
    c6_ret = []
    for alpha in CONFIG["c6_retention_alphas"]:
        W = weighted_shrinkage_operator(B, alpha)
        rhoX = retention_spectrum(QA, W)
        diff_vs_no_prior = float(np.max(np.abs(rhoX - rho_no_prior)))
        c6_ret.append(
            {
                "alpha": alpha,
                "rho_X_min": float(rhoX.min()),
                "rho_X_max": float(rhoX.max()),
                "all_in_0_1": bool(rhoX.min() >= -1e-10 and rhoX.max() <= 1.0 + 1e-10),
                "max_abs_diff_sorted_rhoX_vs_rho_no_prior": diff_vs_no_prior,
                "rho_X_desc": [float(x) for x in rhoX],
            }
        )
    c6 = {
        "alpha_grid": [float(x) for x in alphas],
        "d0_rel_F": d0_list,
        "dinf_rel_F": dinf_list,
        "sweep_rows": [
            {
                "alpha": float(a),
                "d0_rel_F": float(d0),
                "dinf_rel_F": float(dinf),
            }
            for a, d0, dinf in zip(alphas, d0_list, dinf_list)
        ],
        "d0_at_1e-6": d0_list[0],
        "d0_min_over_grid": float(min(d0_list)),
        "dinf_at_1e4": dinf_list[-1],
        "dinf_min_over_grid": float(min(dinf_list)),
        "finite_prior_retention": c6_ret,
        "note": (
            "d0 = ||K_eff-K_SLAM||_F/||K_IS||_F; dinf = "
            "||K_eff-K_IS||_F/||K_IS||_F.  Left-grid d0 is not small for "
            "this geometry: sigma_min(B_R) ~ 6.7e-4 has sigma^2 ~ 4.5e-7 < "
            "alpha = 1e-6, so the weakest pose direction is still weakly "
            "removed at the grid's left edge (convergence to K_SLAM is "
            "monotone as alpha -> 0)."
        ),
    }

    # ---------------- c7 singular prior ------------------------------------
    P_sing = np.eye(q, dtype=float)
    P_sing[q - 1, q - 1] = 0.0
    d0s_list, dinfs_list = [], []
    for alpha in alphas:
        C_s = B.T @ B + alpha * P_sing
        K_s = KIS - A.T @ B @ np.linalg.pinv(C_s) @ B.T @ A
        d0s_list.append(float(np.linalg.norm(K_s - KSL, ord="fro") / KIS_F))
        dinfs_list.append(float(np.linalg.norm(K_s - KIS, ord="fro") / KIS_F))
    dinf_sing_1e4 = dinfs_list[-1]
    dinf_reg_1e4 = dinf_list[-1]
    c7 = {
        "singular_prior_P_diag": P_sing.tolist(),
        "unpenalized_pose_coordinate": int(q - 1),
        "alpha_grid": [float(x) for x in alphas],
        "d0_sing_rel_F": d0s_list,
        "dinf_sing_rel_F": dinfs_list,
        "sweep_rows": [
            {
                "alpha": float(a),
                "d0_sing_rel_F": float(d0),
                "dinf_sing_rel_F": float(dinf),
            }
            for a, d0, dinf in zip(alphas, d0s_list, dinfs_list)
        ],
        "d0_sing_at_1e-6": d0s_list[0],
        "dinf_sing_at_1e4": dinf_sing_1e4,
        "dinf_regular_at_1e4": dinf_reg_1e4,
        "ratio_singular_to_regular_dinf_1e4": float(dinf_sing_1e4 / max(dinf_reg_1e4, 1e-300)),
        "pass_nonconvergence_gate": bool(
            dinf_sing_1e4 > max(1e-6, 10.0 * dinf_reg_1e4)
        ),
        "example_magnitude_gate_gt_0_1": bool(dinf_sing_1e4 > 0.1),
        "note": (
            "The unpenalized coordinate is pose-column q-1 = 18 (last pose, "
            "theta).  dinf_sing(alpha=1e4) stays nonzero (non-convergence to "
            "K_IS), but its magnitude ~1e-4..5e-4 is below the illustrative "
            "0.1 level because the residual gauge direction has modest data "
            "sensitivity in this geometry; the gap over the regular-prior "
            "dinf at the same alpha is two to three orders."
        ),
    }

    # ---------------- c8 duplicate invariance ------------------------------
    c = CONFIG["c8_duplicate_c"]
    scale_dup = 1.0 + c**2
    Adup = np.vstack([A, c * A])
    Bdup = np.vstack([B, c * B])
    KIS2 = Adup.T @ Adup
    rB2, Z2, _, _ = range_basis(Bdup)
    P2 = np.eye(Adup.shape[0], dtype=float) - Z2 @ Z2.T
    KSL2 = Adup.T @ (P2 @ Adup)
    dA2 = thin_decomposition(Adup)
    Q2 = dA2["Q"]
    rho2 = retention_spectrum(Q2, P2)

    alpha_detail = CONFIG["c8_detail_alpha"]
    W1 = weighted_shrinkage_operator(B, alpha_detail)
    W2 = weighted_shrinkage_operator(Bdup, alpha_detail)
    rhoX1 = retention_spectrum(QA, W1)
    rhoX2 = retention_spectrum(Q2, W2)
    W2j = np.eye(Adup.shape[0]) - Bdup @ np.linalg.solve(
        Bdup.T @ Bdup + scale_dup * alpha_detail * np.eye(q, dtype=float),
        Bdup.T,
    )
    Keff1 = keff_of_w(A, W1)
    Keff2j = keff_of_w(Adup, W2j)
    rhoX1j = retention_spectrum(QA, W1)
    rhoX2j = retention_spectrum(Q2, W2j)

    fixed_prior_effect = []
    for alpha_f in CONFIG["c8_fixed_prior_alphas"]:
        Wf1 = weighted_shrinkage_operator(B, alpha_f)
        Wf2 = weighted_shrinkage_operator(Bdup, alpha_f)
        r1 = retention_spectrum(QA, Wf1)
        r2 = retention_spectrum(Q2, Wf2)
        fixed_prior_effect.append(
            {
                "alpha": alpha_f,
                "max_abs_diff_rhoX2_fixed_prior": float(np.max(np.abs(r2 - r1))),
            }
        )
    no_prior_effect = float(np.max(np.abs(rho2 - Rop)))
    c8 = {
        "duplicate_factor_c": c,
        "scale_1_plus_c2": scale_dup,
        "rel_KIS2_vs_scaled_KIS": float(
            np.linalg.norm(KIS2 - scale_dup * KIS, ord="fro") / KIS_F
        ),
        "rel_KSL2_vs_scaled_KSL": float(
            np.linalg.norm(KSL2 - scale_dup * KSL, ord="fro")
            / max(np.linalg.norm(KSL, ord="fro"), 1e-300)
        ),
        "no_prior_max_abs_rho2_minus_rho": no_prior_effect,
        "no_prior_pass_1e-10": bool(no_prior_effect < 1e-10),
        "fixed_prior_alpha_detail": alpha_detail,
        "fixed_prior_max_abs_rhoX2_minus_rhoX": float(
            np.max(np.abs(rhoX2 - rhoX1))
        ),
        "fixed_prior_effect_nonzero": bool(
            np.max(np.abs(rhoX2 - rhoX1)) > 1e-8
        ),
        "joint_scaled_prior_rel_K_eff": float(
            np.linalg.norm(Keff2j - scale_dup * Keff1, ord="fro")
            / max(np.linalg.norm(Keff1, ord="fro"), 1e-300)
        ),
        "joint_scaled_prior_max_abs_rhoX2j_minus_rhoX": float(
            np.max(np.abs(rhoX2j - rhoX1j))
        ),
        "joint_scaled_pass_1e-10": bool(
            np.linalg.norm(Keff2j - scale_dup * Keff1, ord="fro")
            / max(np.linalg.norm(Keff1, ord="fro"), 1e-300)
            < 1e-10
            and np.max(np.abs(rhoX2j - rhoX1j)) < 1e-10
        ),
        "fixed_prior_effect_curve": fixed_prior_effect,
        "note": (
            "Adup = [A; c A], Bdup = [B; c B] (duplicate data blocks, not a "
            "rescaled single copy).  K_IS2 = Adup^T Adup = (1+c^2) K_IS and "
            "likewise K_SLAM2 in exact arithmetic; retention spectra are "
            "invariant without a prior or under the jointly scaled prior "
            "(J_X2 = (1+c^2) alpha I), while a fixed finite prior shifts the "
            "data-to-prior weighting."
        ),
    }

    # Summary metadata
    sigmaA = da["singular_values"]
    checks_seconds = time.perf_counter() - t0
    return {
        "name": name,
        "n_param": int(n),
        "m_data": int(m),
        "q_pose": int(q),
        "p_smooth": int(sb_cfg["p"]) if sb_cfg else None,
        "ranks_summary": {
            "r_A": rA,
            "r_B": rB,
            "sigma_A_max": float(sigmaA[0]),
            "sigma_A_min": float(sigmaA[-1]),
            "sigma_B_max": float(svB[0]),
            "sigma_B_min": float(svB[-1]),
            "cond_A": float(sigmaA[0] / max(sigmaA[-1], 1e-300)),
        },
        "checks_seconds": checks_seconds,
        "check_c1_kernel_identity": c1,
        "check_c2_rank_identity": c2,
        "check_c3_retention_spectrum": c3,
        "check_c4_loss_rank": c4,
        "check_c5_ordering_interlacing": c5,
        "check_c6_prior_sweep": c6,
        "check_c7_singular_prior": c7,
        "check_c8_duplicate_invariance": c8,
    }


def run_resolution_smooth(N: int) -> dict:
    """Smooth-basis core checks (c1-c3 max residuals only) at N > 16."""
    t_build = time.perf_counter()
    cfg = dict(CONFIG)
    cfg["N"] = N
    poses = build_poses(cfg)
    points, h = hh.make_grid(N)
    chi0 = make_chi0(points, cfg)
    A, B, _, _ = hh.build_AB(
        chi0, poses, np.asarray(cfg["rx_offsets"], dtype=float),
        np.asarray(cfg["tx_offset"], dtype=float), N, cfg["k_b"],
    )
    A_R, B_R = hh.whiten_realify(A, B, None)
    S = build_smooth_basis(points, CONFIG["smooth_basis"])
    A_s = A_R @ S
    build_seconds = time.perf_counter() - t_build

    t_check = time.perf_counter()
    da = thin_decomposition(A_s)
    rA = da["rank"]
    QA = da["Q"]
    rB, Z, _, _ = range_basis(B_R)
    Pperp = np.eye(A_s.shape[0]) - Z @ Z.T
    KSL = A_s.T @ (Pperp @ A_s)
    C = Pperp @ A_s

    # c1 summary
    w_eig, V_eig = np.linalg.eigh(KSL)
    _, svKSL, tolKSL = rank_svd(KSL)
    k1 = int(np.sum(np.abs(w_eig) <= tolKSL))
    N1 = V_eig[:, :k1]
    rC, _, _ = rank_svd(C)
    _, _, vhC = np.linalg.svd(C, full_matrices=True)
    N2 = vhC[rC:].T
    c1_max_residual = {
        "subspace_distance": float(
            np.linalg.norm(N1 @ N1.T - N2 @ N2.T, ord=2)
        )
        if N1.shape[1] and N2.shape[1] else 0.0,
        "max_norm_KSL_u_over_N2": float(
            max((np.linalg.norm(KSL @ N2[:, j], ord=2)
                 for j in range(N2.shape[1])), default=0.0)
        ),
        "max_norm_C_u_over_N1": float(
            max((np.linalg.norm(C @ N1[:, j], ord=2)
                 for j in range(N1.shape[1])), default=0.0)
        ),
    }
    # c2 summary
    rKIS, _, _ = rank_svd(A_s.T @ A_s)
    rKS, _, _ = rank_svd(KSL)
    rAB, _, _ = rank_svd(np.hstack([A_s, B_R]))
    d_inter = rA + rB - rAB
    c2_max_residual = {
        "r_A": rA, "r_KIS": rKIS, "r_KS": rKS, "r_B": rB, "r_AB": rAB,
        "d_inter": int(d_inter),
        "residual_lhs_minus_d_inter": int(rKIS - rKS - d_inter),
    }
    # c3 summary
    rho = retention_spectrum(QA, Pperp)
    c2_vals = np.linalg.svd(Z.T @ QA, compute_uv=False) ** 2
    rho_pred = _desc(np.concatenate([1.0 - c2_vals, np.ones(max(rA - len(c2_vals), 0))]))
    c3_max_residual = {
        "max_abs_diff_rho_vs_pred": float(np.max(np.abs(rho - rho_pred))),
        "count_less_than_one": int(np.sum(rho < 1.0 - 1e-8)),
        "r_B": rB,
    }
    check_seconds = time.perf_counter() - t_check
    return {
        "N": N,
        "grid_h": float(h),
        "build_seconds": build_seconds,
        "check_seconds": check_seconds,
        "r_B": rB,
        "c1_max_residual": c1_max_residual,
        "c2_max_residual": c2_max_residual,
        "c3_max_residual": c3_max_residual,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    t0 = datetime.now(timezone.utc)
    t_start = time.perf_counter()

    points, h = hh.make_grid(cfg["N"])
    poses = build_poses(cfg)
    chi0 = make_chi0(points, cfg)
    rx_offsets = np.asarray(cfg["rx_offsets"], dtype=float)
    tx_offset = np.asarray(cfg["tx_offset"], dtype=float)

    t_build = time.perf_counter()
    A, B, F_full, _ = hh.build_AB(
        chi0, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
    )
    build_seconds = time.perf_counter() - t_build
    A_pix_R, B_R = hh.whiten_realify(A, B, None)
    S = build_smooth_basis(points, cfg["smooth_basis"])
    A_smooth_R = A_pix_R @ S
    print(f"[family2] build N={cfg['N']}: A_pix_R {A_pix_R.shape}, "
          f"B_R {B_R.shape}, S {S.shape}, A_smooth_R {A_smooth_R.shape}")

    case_pixel = run_algebraic_case("pixel", A_pix_R, B_R, None)
    case_smooth = run_algebraic_case("smooth", A_smooth_R, B_R, cfg["smooth_basis"])

    resolutions = []
    for N_res in cfg["resolution_Ns"]:
        t_r = time.perf_counter()
        res = run_resolution_smooth(N_res)
        res["total_seconds"] = time.perf_counter() - t_r
        resolutions.append(res)
        print(f"[family2] resolution N={N_res}: "
              f"build={res['build_seconds']:.2f}s check={res['check_seconds']:.2f}s "
              f"c3 residual={res['c3_max_residual']['max_abs_diff_rho_vs_pred']:.3e}")

    total_runtime = time.perf_counter() - t_start
    source_files = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
    ]
    results = {
        "generated_utc": t0.isoformat(),
        "runner": "src/family2_algebraic_spine.py",
        "command": ".venv/bin/python src/family2_algebraic_spine.py",
        "family": 2,
        "self_cell_formula": hh.SELF_CELL_FORMULA,
        "self_cell_formula_version": hh.SELF_CELL_FORMULA_VERSION,
        "parent_correction_gate": (
            "context/PARENT_CORRECTIONS.md and workshop rule 16 (2026-09-03); "
            "this file only consumes the corrected helmholtz build_AB"
        ),
        "runtime_seconds": total_runtime,
        "build_seconds_N16": build_seconds,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
        },
        "environment_note": (
            "Apple Silicon CPU, no GPU/MPS/CUDA; deterministic dense "
            "numpy/scipy linear algebra"
        ),
        "source_sha256": {
            p: hashlib.sha256((_ROOT / p).read_bytes()).hexdigest()
            for p in source_files
        },
        "config": {
            **cfg,
            "k_b": float(cfg["k_b"]),
            "poses": [p.tolist() for p in poses],
            "rx_offsets": rx_offsets.tolist(),
            "tx_offset": tx_offset.tolist(),
            "grid_h_cell": h,
            "chi0_stats": {
                "min": float(chi0.min()),
                "max": float(chi0.max()),
                "mean": float(chi0.mean()),
                "l2": float(np.linalg.norm(chi0)),
            },
        },
        "dimensions_common": {
            "m_real_data": int(A_pix_R.shape[0]),
            "q_pose": int(B_R.shape[1]),
            "S_N2": int(A_pix_R.shape[1]),
            "p_smooth": int(S.shape[1]),
            "T": cfg["T"],
            "n_rx": cfg["n_rx"],
        },
        "cases": {"pixel": case_pixel, "smooth": case_smooth},
        "resolution_smooth_c1c2c3": resolutions,
        "linear_algebra_scope_note": (
            "Finite-dimensional whitened/realified dense linear algebra on the "
            "validated Family-1 Jacobians only.  These checks establish "
            "identities for the specific discrete (A_R, B_R); they are not "
            "continuum-limit or inverse-problem recovery claims."
        ),
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)
    results_path = results_dir / "family2_results.json"
    results_path.write_text(_round_trip_json(results) + "\n")

    # ------------------------- figures --------------------------------------
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    styles = {"pixel": ("#1f77b4", "o"), "smooth": ("#d62728", "s")}
    for cname, col, mk in [("pixel", "#1f77b4", "o"), ("smooth", "#d62728", "s")]:
        cs = case_pixel if cname == "pixel" else case_smooth
        rho = np.asarray(cs["check_c3_retention_spectrum"]["rho_desc"])
        pred = np.asarray(cs["check_c3_retention_spectrum"]["rho_pred_desc"])
        x = np.arange(1, len(rho) + 1)
        ax.plot(x, rho, mk + "-", color=col, ms=4.5, lw=1.1,
                label=f"{cname} rho")
        ax.plot(x, pred, "--", color=col, lw=0.9, alpha=0.75,
                label=f"{cname} predicted (1-cos^2, padded)")
    ax.axvline(B_R.shape[1] + 0.5, color="k", ls=":", lw=1.2,
               label=f"rank B = {B_R.shape[1]}")
    ax.set_xlabel("sorted index (descending retention)")
    ax.set_ylabel("retention eigenvalue rho")
    ax.set_title("Family 2: retention spectra, pixel vs smooth basis\n"
                 f"(N={cfg['N']}, m=48, q=18; corrected Helmholtz harness)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig_path1 = figures_dir / "family2_retention_spectra.png"
    fig.savefig(fig_path1, dpi=200)
    plt.close(fig)

    # Prior sweep: d0 and dinf vs alpha (regular prior), plus singular dinf.
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    for cname, col, mk in [("pixel", "#1f77b4", "o"), ("smooth", "#d62728", "s")]:
        cs = case_pixel if cname == "pixel" else case_smooth
        al = np.asarray(cs["check_c6_prior_sweep"]["alpha_grid"])
        ax.loglog(al, cs["check_c6_prior_sweep"]["d0_rel_F"],
                  mk + "-", color=col, ms=4.0, lw=1.1,
                  label=f"{cname} d0(alpha)=||K_eff-K_SLAM||/||K_IS||")
        ax.loglog(al, cs["check_c6_prior_sweep"]["dinf_rel_F"],
                  mk + "--", color=col, ms=4.0, lw=1.1,
                  label=f"{cname} dinf(alpha)=||K_eff-K_IS||/||K_IS||")
        dinf_s = cs["check_c7_singular_prior"]["dinf_sing_at_1e4"]
        ax.axhline(dinf_s, color=col, ls=":", lw=1.0)
    ax.set_xlabel(r"regular-prior strength $\alpha$ (log)")
    ax.set_ylabel("relative Frobenius distance")
    ax.set_title("Family 2: prior sweep d0/dinf and singular-prior limit\n"
                 "(horizontal dotted = singular-prior dinf at alpha=1e4)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig_path2 = figures_dir / "family2_prior_sweep.png"
    fig.savefig(fig_path2, dpi=200)
    plt.close(fig)

    # Duplicate prior effect.
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    alphas_f = cfg["c8_fixed_prior_alphas"]
    for cname, col, mk in [("pixel", "#1f77b4", "o"), ("smooth", "#d62728", "s")]:
        cs = case_pixel if cname == "pixel" else case_smooth
        curve = cs["check_c8_duplicate_invariance"]["fixed_prior_effect_curve"]
        ys = [r["max_abs_diff_rhoX2_fixed_prior"] for r in curve]
        ax.loglog(alphas_f, ys, mk + "-", color=col, ms=5.0, lw=1.2,
                  label=f"{cname}: fixed prior J_X=alpha I (duplicated data)")
    ax.axhline(1e-15, color="k", ls=":", lw=1.0,
               label="no-prior invariance level (~1e-15)")
    ax.set_xlabel(r"fixed prior strength $\alpha$")
    ax.set_ylabel(r"max |rho_X(duplicate) - rho_X(single)|")
    ax.set_title("Family 2: c-scaled duplicate data - fixed-prior effect\n"
                 "(no prior and jointly scaled prior are invariant)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig_path3 = figures_dir / "family2_duplicate_prior_effect.png"
    fig.savefig(fig_path3, dpi=200)
    plt.close(fig)

    # ------------------------- console summary -----------------------------
    def gate(cs, check):
        return cs[check].get("pass_gate")

    print("\n===== FAMILY 2 ALGEBRAIC SPINE SUMMARY =====")
    for cname, cs in [("pixel", case_pixel), ("smooth", case_smooth)]:
        print(f"[{cname}] ranks rA={cs['check_c2_rank_identity']['r_A']} "
              f"rB={cs['check_c2_rank_identity']['r_B']} "
              f"rKS={cs['check_c2_rank_identity']['r_KS']} "
              f"d_inter={cs['check_c2_rank_identity']['d_inter']} "
              f"c2 residual={cs['check_c2_rank_identity']['residual_lhs_minus_d_inter']}")
        print(f"[{cname}] c1 gate pass={gate(cs, 'check_c1_kernel_identity')} "
              f"dist={cs['check_c1_kernel_identity']['subspace_distance']:.3e}")
        print(f"[{cname}] c3 gate pass={cs['check_c3_retention_spectrum']['pass_pred_gate']} "
              f"maxdiff={cs['check_c3_retention_spectrum']['max_abs_diff_rho_vs_pred']:.3e} "
              f"count<1={cs['check_c3_retention_spectrum']['count_less_than_one_minus_1e-8']}")
        print(f"[{cname}] c4 gate pass={gate(cs, 'check_c4_loss_rank')}")
        print(f"[{cname}] c5 gate pass={gate(cs, 'check_c5_ordering_interlacing')}")
        print(f"[{cname}] c6 d0@1e-6={cs['check_c6_prior_sweep']['d0_at_1e-6']:.4e} "
              f"dinf@1e4={cs['check_c6_prior_sweep']['dinf_at_1e4']:.3e}")
        print(f"[{cname}] c7 dinf_sing@1e4={cs['check_c7_singular_prior']['dinf_sing_at_1e4']:.4e} "
              f"pass_nonconvergence={cs['check_c7_singular_prior']['pass_nonconvergence_gate']}")
        c8d = cs["check_c8_duplicate_invariance"]
        print(f"[{cname}] c8 relKIS={c8d['rel_KIS2_vs_scaled_KIS']:.2e} "
              f"relKSL={c8d['rel_KSL2_vs_scaled_KSL']:.2e} "
              f"no-prior rho diff={c8d['no_prior_max_abs_rho2_minus_rho']:.2e} "
              f"fixed-prior effect={c8d['fixed_prior_max_abs_rhoX2_minus_rhoX']:.4e} "
              f"joint-pass={c8d['joint_scaled_pass_1e-10']}")
    print(f"total runtime {total_runtime:.2f}s; "
          f"resolution runs {[(r['N'], r['total_seconds']) for r in resolutions]}")
    print("results ->", results_path)
    print("figure  ->", fig_path1)
    print("figure  ->", fig_path2)
    print("figure  ->", fig_path3)


if __name__ == "__main__":
    main()
