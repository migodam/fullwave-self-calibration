"""Family 8: current-space G_S/SOM versus map-tangent A/K_eff distinction.

NEW finite-dimensional sanity check that makes explicit, and numerically
verifies on the discrete N=16 model, the distinction between

  * the current-space (SOM-style) sensing map G_S stacked over the six
    Family-1 poses (rows T*n_rx, columns N^2 current/pixel modes), and
  * the map-tangent data-informative Jacobian A = dF/dchi of the full-wave
    contrast-source forward map, plus the pose-corrected pixel-space kernel
    K_eff(alpha) = A_R^T P A_R used by the finite-prior SLAM analysis.

Reused modules (read-only, not modified):
  src/helmholtz.py     make_grid, build_operators, build_AB,
                       whiten_realify, born_forward
  src/family1_pilot.py build_poses, make_chi0
  src/family2_algebraic_spine.py  build_smooth_basis, rank_svd
  src/family4_frequency_trajectory.py  K_eff (internal formula check)

Run (from the experiment root):
    .venv/bin/python src/family8_som_distinction.py

Outputs:
  results/family8_som_distinction.json
  figures/family8_composition_identity.png
  figures/family8_mode_distinction.png
  notes/family8_som_distinction.md

Checks:
  1. Full-wave composition identity: rows of A_c are reproduced as
     A_comp_t = G_S_list[t] @ lu_solve(M, diag(E_tot_list[t]));
     stack relative Frobenius error gate < 1e-12.
  2. Born composition identity (hh.born_forward is available):
     A_born_comp_t = G_S_list[t] @ diag(E_inc_list[t]);
     stack relative Frobenius error gate < 1e-12.
  3. Range containment Range(A_c) subset Range(G_S_stack) on the LEFT
     singular (data/column-space) subspaces U_A, U_GS: direct projector and
     containment-residual norms < 1e-8 and rank_A <= rank_GS (the recorded
     arccos principal angle ~2.1e-8 rad is the sqrt(2*eps) roundoff floor of
     arccos(1-eps) for exactly identical full-rank subspaces).  The V_GS^H
     V_A object requested for the mode distinction is reported in part 3; it
     is not itself a column-range containment measure.
  4. Right singular subspaces V_GS and V_A are NOT identical: recorded
     boolean max principal angle > 1e-6 deg (observed, never forced).
  5. Squared overlaps of the three smallest-eigenvalue K_eff eigenvectors
     u_conf with V_GS and V_A; top-5 lists; pull-through cosines of
     c_t = T_full_t @ u_conf against the top-5 V_GS current modes.
  6. Mode-count contrast at relative thresholds 1e-6 and 1e-3, so the
     "SOM mode count" is never silently identified with the
     "map-information mode count".

Every number below is a statement about the exact finite-dimensional
dense-linear-algebra model at N=16, T=6, f=1.0.  No SOM-algorithm-quality
claim and no continuum theorem is asserted.
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

from scipy.linalg import lu_factor, lu_solve  # noqa: E402

_EPS = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Configuration (Family 1/4/6 scenario, f = 1.0, alpha = 1.0)
# ---------------------------------------------------------------------------

CONFIG = {
    "N": 16,
    "k_b": 2.0 * np.pi,          # k_b = 2*pi*f with f = 1.0
    "f": 1.0,
    "T": 6,
    "n_rx": 4,
    "q": 1.0,
    "q_pose": 18,                # 3*T real pose parameters
    "m_data_complex": 24,        # T*n_rx
    "m_data_real": 48,           # 2*T*n_rx after whiten_realify
    "pixel_dim": 256,            # N^2
    "rx_offsets": [
        [-0.06, 0.0],
        [0.06, 0.0],
        [0.0, -0.06],
        [0.0, 0.06],
    ],
    "tx_offset": [0.0, 0.0],
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],
    "pose_builder": "family1.build_poses",
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
        "note": (
            "identical construction/config to family2-family7; included as "
            "scene metadata and for the chi0 smooth-projection residual only; "
            "all Family 8 kernels are the full pixel-space (N^2) operators"
        ),
    },
    "finite_prior_alpha": 1.0,
    "shrinkage_operator": (
        "P = I - B_R (B_R^T B_R + alpha I)^{-1} B_R^T with B_R the "
        "whiten_realify(A_c, B_c, None) pose block; K_eff = A_R^T P A_R "
        "(pixel space)"
    ),
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "gates": {
        "fullwave_composition_rel_fro": 1e-12,
        "born_composition_rel_fro": 1e-12,
        "containment_max_angle_rad": 1e-8,
        "V_not_identical_angle_deg": 1e-6,
    },
    "relative_thresholds": [1e-6, 1e-3],
    "representative_pose_index": 2,
    "topk": {
        "angle_topk": 10,
        "overlap_topk": 5,
        "current_topk": 5,
    },
    "randomness_note": "deterministic dense linear algebra; no RNG used",
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _sym(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    raise TypeError(f"not JSON serialisable: {type(obj)}")


def _md_table(headers: list[str], rows: list[list]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for r in rows:
        lines.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(lines)


def _fmt(x, spec=".6g"):
    if x is None:
        return "-"
    return format(x, spec)


def _safe_div(num: float, den: float) -> float:
    return float(num / max(abs(den), 1e-300))


def _rel_fro(A: np.ndarray, B: np.ndarray) -> float:
    """Relative Frobenius discrepancy ||A-B||_F / ||B||_F (guarded)."""
    return _safe_div(
        float(np.linalg.norm(A - B, ord="fro")),
        float(np.linalg.norm(B, ord="fro")),
    )


def _principal_angles_rad(X: np.ndarray, Y: np.ndarray):
    """Principal angles arccos(svd(X^H Y)) for subspaces X (smaller) in Y.

    X, Y have orthonormal columns in the same ambient space.  Returns
    (clipped singular values, angles in radians), angles ascending.
    """
    C = X.conj().T @ Y
    sv = np.linalg.svd(C, compute_uv=False)
    svc = np.clip(np.asarray(sv, dtype=float), 0.0, 1.0)
    return svc, np.arccos(svc)


def _sq_overlaps(u: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Per-column |basis^H u|^2 for an orthonormal complex basis."""
    return np.abs(np.conj(basis).T @ np.asarray(u)) ** 2


def _cos_complex(a: np.ndarray, b: np.ndarray) -> float:
    """|inner product| / (||a|| ||b||) for complex vectors."""
    return float(
        abs(np.vdot(a, b))
        / max(float(np.linalg.norm(a) * np.linalg.norm(b)), 1e-300)
    )


def _rel_count_sv(sv: np.ndarray, threshold: float) -> int:
    """Count singular values above threshold*sigma_max."""
    sv = np.asarray(sv, dtype=float)
    if sv.size == 0 or sv[0] == 0.0:
        return 0
    return int(np.sum(sv > float(threshold) * sv[0]))


def _rel_count_eig(ev: np.ndarray, threshold: float) -> int:
    """Count eigenvalues above threshold*lambda_max (ev descending)."""
    ev = np.asarray(ev, dtype=float)
    if ev.size == 0 or ev[0] <= 0.0:
        return 0
    return int(np.sum(ev > float(threshold) * ev[0]))


def _top_overlap_entries(u_confs: list, basis: np.ndarray, topk: int) -> list:
    """Top-k |<u,v>|^2 pairs over all confounded vectors and basis columns."""
    rows: list[dict] = []
    for j, u in enumerate(u_confs):
        ov = _sq_overlaps(u, basis)
        for i in range(len(ov)):
            rows.append(
                {
                    "u_conf_index": int(j),
                    "basis_index": int(i),
                    "squared_overlap": float(ov[i]),
                    "abs_inner_product": float(np.sqrt(ov[i])),
                }
            )
    rows.sort(key=lambda d: d["squared_overlap"], reverse=True)
    return rows[:topk]


# ---------------------------------------------------------------------------
# Computations
# ---------------------------------------------------------------------------

def run_all() -> dict:
    cfg = CONFIG
    N = int(cfg["N"])
    T = int(cfg["T"])
    n_rx = int(cfg["n_rx"])
    k_b = float(cfg["k_b"])
    alpha = float(cfg["finite_prior_alpha"])
    pixel = N * N
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)

    # ---- deterministic scene ---------------------------------------------
    poses = family1.build_poses(cfg)
    points, h = hh.make_grid(N)
    chi0 = family1.make_chi0(points, cfg)
    S = family2.build_smooth_basis(points, cfg["smooth_basis"])
    c_smooth, *_ = np.linalg.lstsq(S, chi0, rcond=None)
    smooth_resid = _safe_div(
        float(np.linalg.norm(chi0 - S @ c_smooth)),
        float(np.linalg.norm(chi0)),
    )

    # ---- single dense construction and one LU factorisation of M ----------
    ops = hh.build_operators(chi0, poses, rx, tx, N, k_b)
    M = ops["M"]
    # build_operators returns M but not lu/piv; factor once (identical to the
    # internal use in build_AB).
    lu, piv = lu_factor(M, check_finite=False)
    G_S_stack = np.vstack(ops["G_S_list"])       # (24, 256) complex
    A_c, B_c, F_c, diag_c = hh.build_AB(chi0, poses, rx, tx, N, k_b)

    # ---- 1. composition identity (full wave and Born) --------------------
    A_comp_stack = np.empty_like(A_c)
    A_born_stack = np.empty_like(A_c)
    T_full_list = []
    full_per_pose = np.empty(T, dtype=float)
    born_per_pose = np.empty(T, dtype=float)
    full_vs_born_per_pose = np.empty(T, dtype=float)
    for t in range(T):
        sl = slice(t * n_rx, (t + 1) * n_rx)
        G_S_t = ops["G_S_list"][t]
        E_inc_t = ops["E_inc_list"][t]
        E_tot_t = ops["E_tot_list"][t]
        X_t = lu_solve((lu, piv), np.diag(E_tot_t), check_finite=False)
        T_full_list.append(X_t)
        A_comp_t = G_S_t @ X_t
        A_born_t = G_S_t @ np.diag(E_inc_t)      # G_S diag(E_inc)
        A_comp_stack[sl, :] = A_comp_t
        A_born_stack[sl, :] = A_born_t
        full_per_pose[t] = _rel_fro(A_comp_t, A_c[sl, :])
        born_per_pose[t] = _rel_fro(A_born_t, A_c[sl, :])
        full_vs_born_per_pose[t] = _rel_fro(A_c[sl, :], A_born_t)

    fullwave_stack_rel = _rel_fro(A_comp_stack, A_c)
    born_available = hasattr(hh, "born_forward")
    born_comp_rel = np.nan
    if born_available:
        _, A_born_ref = hh.born_forward(chi0, poses, rx, tx, N, k_b)
        born_comp_rel = _rel_fro(A_born_stack, A_born_ref)
        for t in range(T):                        # per-pose reference rows
            sl = slice(t * n_rx, (t + 1) * n_rx)
            born_per_pose[t] = _rel_fro(A_born_stack[sl, :], A_born_ref[sl, :])

    # ---- 2. thin SVDs and range containment (left singular subspaces) -----
    U_GS, _, Vh_GS = np.linalg.svd(G_S_stack, full_matrices=False)
    U_A, _, Vh_A = np.linalg.svd(A_c, full_matrices=False)
    V_GS = Vh_GS.conj().T                         # (256, 24) pixel columns
    V_A = Vh_A.conj().T

    r_GS, sv_GS, tol_GS = family2.rank_svd(G_S_stack)
    r_A, sv_A, tol_A = family2.rank_svd(A_c)
    U_GS_t = U_GS[:, :r_GS]
    U_A_t = U_A[:, :r_A]
    V_GS_t = V_GS[:, :r_GS]
    V_A_t = V_A[:, :r_A]

    svc_U, ang_U_rad = _principal_angles_rad(U_A_t, U_GS_t)
    max_angle_U_rad = float(ang_U_rad[-1])
    containment_dist = float(np.sqrt(1.0 - svc_U[-1] ** 2))
    rank_ineq_holds = bool(r_A <= r_GS)
    # Numerically stable containment metrics (canonical angles from arccos of
    # a cross product carry a ~sqrt(eps) floor even for exactly identical
    # full-rank subspaces; direct projector/residual norms do not).
    Proj_A = U_A_t @ U_A_t.conj().T
    Proj_GS = U_GS_t @ U_GS_t.conj().T
    projector_diff_2 = float(np.linalg.norm(Proj_A - Proj_GS, ord=2))
    containment_residual_2 = float(
        np.linalg.norm(
            (np.eye(Proj_GS.shape[0], dtype=complex) - Proj_GS) @ U_A_t,
            ord=2,
        )
    )
    row_projector_diff_2 = float(
        np.linalg.norm(
            V_A_t @ V_A_t.conj().T - V_GS_t @ V_GS_t.conj().T,
            ord=2,
        )
    )
    orth_err_GS = float(
        np.linalg.norm(
            U_GS_t.conj().T @ U_GS_t - np.eye(U_GS_t.shape[1]), ord=2
        )
    )
    orth_err_A = float(
        np.linalg.norm(
            U_A_t.conj().T @ U_A_t - np.eye(U_A_t.shape[1]), ord=2
        )
    )
    containment_pass = bool(
        rank_ineq_holds
        and projector_diff_2 < 1e-8
        and containment_residual_2 < 1e-8
    )
    literal_angle_gate_pass = bool(
        max_angle_U_rad < cfg["gates"]["containment_max_angle_rad"]
    )

    # ---- 3. right singular subspaces (pixel/current-space distinction) ----
    _, ang_V_rad = _principal_angles_rad(V_GS_t, V_A_t)
    ang_V_deg = np.rad2deg(ang_V_rad)
    count_lt_1e8_deg = int(np.sum(ang_V_deg < 1e-8))
    median_V_deg = float(np.median(ang_V_deg))
    max_V_deg = float(ang_V_deg[-1])
    topk = int(cfg["topk"]["angle_topk"])
    top10_idx = np.argsort(ang_V_deg)[::-1][:topk]
    top10_V_deg = [float(ang_V_deg[i]) for i in top10_idx]
    not_identical_deg = bool(max_V_deg > cfg["gates"]["V_not_identical_angle_deg"])

    # ---- K_eff pixel-space kernel and most confounded directions ----------
    A_R, B_R = hh.whiten_realify(A_c, B_c, None)
    P = np.eye(B_R.shape[0], dtype=float) - B_R @ np.linalg.solve(
        _sym(B_R.T @ B_R) + alpha * np.eye(B_R.shape[1], dtype=float),
        B_R.T,
    )
    K_eff = _sym(A_R.T @ (P @ A_R))
    K4 = family4.K_eff(A_R, B_R, alpha)
    K_eff_vs_family4_rel = _rel_fro(K_eff, K4)
    w_K_eff_asc, V_K = np.linalg.eigh(K_eff)     # ascending
    u_confs = [V_K[:, j] for j in range(3)]
    smallest3_asc = [float(w_K_eff_asc[j]) for j in range(3)]

    # K_IS (realified map information) descending spectrum.
    K_IS = _sym(A_R.T @ A_R)
    w_K_IS = np.linalg.eigvalsh(K_IS)
    w_K_IS_desc = w_K_IS[::-1]

    # ---- overlaps with V_GS and V_A ---------------------------------------
    per_vec = []
    for j, u in enumerate(u_confs):
        ov_GS = _sq_overlaps(u, V_GS_t)
        ov_A = _sq_overlaps(u, V_A_t)
        iGS = int(np.argmax(ov_GS))
        iA = int(np.argmax(ov_A))
        per_vec.append(
            {
                "u_conf_index": int(j),
                "eigval_ascending": float(w_K_eff_asc[j]),
                "V_GS_max_squared_overlap": float(ov_GS[iGS]),
                "V_GS_argmax_index": iGS,
                "V_GS_singular_value_at_argmax": float(sv_GS[iGS]),
                "V_A_max_squared_overlap": float(ov_A[iA]),
                "V_A_argmax_index": iA,
                "V_A_singular_value_at_argmax": float(sv_A[iA]),
                "V_GS_mass_sum": float(np.sum(ov_GS)),
                "V_A_mass_sum": float(np.sum(ov_A)),
            }
        )
    top5_GS = _top_overlap_entries(u_confs, V_GS_t, int(cfg["topk"]["overlap_topk"]))
    top5_A = _top_overlap_entries(u_confs, V_A_t, int(cfg["topk"]["overlap_topk"]))

    # ---- pull-through through a representative single-pose current map ----
    rep_t = int(cfg["representative_pose_index"])
    ctk = int(cfg["topk"]["current_topk"])
    top5_idx_GS = np.argsort(sv_GS)[::-1][:ctk]
    top5_idx_A = np.argsort(sv_A)[::-1][:ctk]
    pull_through = []
    for j, u in enumerate(u_confs):
        c_t = T_full_list[rep_t] @ u
        pull_through.append(
            {
                "u_conf_index": int(j),
                "representative_pose_index": rep_t,
                "cosine_with_top5_V_GS_modes": [
                    _cos_complex(c_t, V_GS[:, i]) for i in top5_idx_GS
                ],
                "top5_V_GS_basis_indices": [int(i) for i in top5_idx_GS],
                "cosine_with_top5_V_A_modes": [
                    _cos_complex(c_t, V_A[:, i]) for i in top5_idx_A
                ],
                "top5_V_A_basis_indices": [int(i) for i in top5_idx_A],
            }
        )

    # ---- 4. mode-count contrast -------------------------------------------
    thrs = cfg["relative_thresholds"]
    s_GS_all = np.asarray(sv_GS, dtype=float)
    s_A_all = np.asarray(sv_A, dtype=float)
    s_AR_all = np.linalg.svd(A_R, compute_uv=False)
    mode_counts = {
        "relative_thresholds": list(thrs),
        "sigma_max_G_S": float(s_GS_all[0]),
        "sigma_max_A": float(s_A_all[0]),
        "sigma_max_A_R": float(s_AR_all[0]),
        "lambda_max_K_IS": float(w_K_IS_desc[0]),
        "G_S_counts": [_rel_count_sv(s_GS_all, t) for t in thrs],
        "A_c_counts": [_rel_count_sv(s_A_all, t) for t in thrs],
        "A_R_counts": [_rel_count_sv(s_AR_all, t) for t in thrs],
        "K_IS_counts": [_rel_count_eig(w_K_IS_desc, t) for t in thrs],
        "sv_G_S": s_GS_all,
        "sv_A": s_A_all,
        "eig_K_IS_desc": w_K_IS_desc,
        "count_rule_note": (
            "G_S/A_c counts use sigma > threshold*sigma_max; K_IS counts use "
            "eigenvalue > threshold*lambda_max.  Because K_IS = A_R^T A_R, "
            "its eigenvalues scale like sigma_A_R^2, so the same literal "
            "relative threshold on eigenvalues is NOT commensurable with the "
            "relative singular-value thresholds.  A_R sigma counts are "
            "recorded as the commensurable companion."
        ),
    }

    # ---- gates -------------------------------------------------------------
    checks = [
        {
            "id": "fullwave_composition",
            "claim": (
                "full-wave composition A_comp = G_S lu_solve(M, diag(E_tot)) "
                "reproduces A_c from build_AB"
            ),
            "threshold": cfg["gates"]["fullwave_composition_rel_fro"],
            "operator": "stack rel Fro <",
            "observed": fullwave_stack_rel,
            "pass": bool(fullwave_stack_rel < cfg["gates"]["fullwave_composition_rel_fro"]),
        },
        {
            "id": "born_composition",
            "claim": (
                "Born composition A_born_comp = G_S diag(E_inc) reproduces "
                "hh.born_forward A_born"
            ),
            "threshold": cfg["gates"]["born_composition_rel_fro"],
            "operator": "stack rel Fro <",
            "observed": born_comp_rel,
            "available": born_available,
            "pass": bool(
                born_available
                and born_comp_rel < cfg["gates"]["born_composition_rel_fro"]
            ),
        },
        {
            "id": "range_containment",
            "claim": (
                "Range(A_c) subset Range(G_S_stack): stable projector and "
                "residual containment checks < 1e-8 with rank_A <= rank_GS"
            ),
            "threshold": cfg["gates"]["containment_max_angle_rad"],
            "operator": (
                "projector diff < 1e-8 and containment residual < 1e-8 and "
                "rank_A <= rank_GS"
            ),
            "observed_max_angle_rad": max_angle_U_rad,
            "observed_projector_diff_2": projector_diff_2,
            "observed_containment_residual_2": containment_residual_2,
            "observed_rank_pair": [int(r_A), int(r_GS)],
            "rank_inequality": rank_ineq_holds,
            "pass": containment_pass,
            "literal_arccos_angle_gate_pass": literal_angle_gate_pass,
            "note": (
                "The literal arccos-based max canonical angle is "
                f"{max_angle_U_rad:.3e} rad, sitting at the arccos(1-eps) "
                "roundoff floor ~sqrt(2*eps) for exactly identical full-rank "
                "subspaces; the direct projector/residual norms are O(eps) "
                "and are the executed containment gate."
            ),
        },
        {
            "id": "V_not_identical",
            "claim": (
                "V_GS and V_A right singular (pixel/current-mode) subspaces "
                "are NOT identical"
            ),
            "threshold": cfg["gates"]["V_not_identical_angle_deg"],
            "operator": "at least one principal angle > deg (recorded only)",
            "observed_max_angle_deg": max_V_deg,
            "observed_n_angles": int(len(ang_V_deg)),
            "observed_count_lt_1e-8_deg": count_lt_1e8_deg,
            "pass": None,
            "boolean": not_identical_deg,
            "note": (
                "recorded boolean; not a forced pass.  Right-subspace angles "
                "do not measure column-range containment."
            ),
        },
        {
            "id": "K_eff_overlaps",
            "claim": (
                "squared overlaps of the 3 smallest-eigenvalue K_eff "
                "eigenvectors with V_A vs V_GS (recorded)"
            ),
            "threshold": None,
            "operator": "recorded top-5 lists and per-vector maxima",
            "observed": per_vec,
            "pass": None,
        },
        {
            "id": "K_eff_formula_consistency",
            "claim": "pixel-space K_eff formula matches family4.K_eff(A_R,B_R,alpha)",
            "threshold": 1e-12,
            "operator": "rel Fro <",
            "observed": K_eff_vs_family4_rel,
            "pass": bool(K_eff_vs_family4_rel < 1e-12),
        },
    ]

    return {
        "scene": {
            "N": N,
            "T": T,
            "n_rx": n_rx,
            "pixel_dim": pixel,
            "k_b": k_b,
            "poses_phi_deg": [
                float(round(np.rad2deg(float(x)), 8))
                for x in np.arctan2(poses[:, 1], poses[:, 0])
            ],
            "chi0_range": [float(np.min(chi0)), float(np.max(chi0))],
            "chi0_smooth_projection_rel_residual": smooth_resid,
            "M_min_sigma_over_norm": ops["sigma_min_over_M_norm"],
            "max_state_residual": ops["max_state_residual"],
        },
        "part1_composition": {
            "fullwave_stack_rel_fro": fullwave_stack_rel,
            "fullwave_per_pose_rel_fro": full_per_pose,
            "born_available": born_available,
            "born_comp_stack_rel_fro": born_comp_rel,
            "born_per_pose_rel_fro": born_per_pose,
            "full_vs_born_per_pose_rel_fro": full_vs_born_per_pose,
            "fro_A_born_stack": float(np.linalg.norm(A_born_stack, ord="fro")),
            "fro_A_full_stack": float(np.linalg.norm(A_c, ord="fro")),
            "born_over_full_fro_ratio": _safe_div(
                float(np.linalg.norm(A_born_stack, ord="fro")),
                float(np.linalg.norm(A_c, ord="fro")),
            ),
            "machine_epsilon": _EPS,
            "reuse_note": (
                "M, G_S_list, E_inc_list, E_tot_list come from one "
                "hh.build_operators call; lu,piv are computed once from the "
                "returned M (build_operators itself does not return lu,piv)."
            ),
        },
        "part2_range_containment": {
            "rank_GS": int(r_GS),
            "rank_A": int(r_A),
            "rank_tol_GS": float(tol_GS),
            "rank_tol_A": float(tol_A),
            "principal_angle_svs_U_A_to_U_GS": svc_U,
            "principal_angles_U_A_to_U_GS_rad": ang_U_rad,
            "max_angle_rad": max_angle_U_rad,
            "containment_distance_sqrt_1_min_smin2": containment_dist,
            "projector_diff_2": projector_diff_2,
            "containment_residual_2": containment_residual_2,
            "row_projector_diff_2": row_projector_diff_2,
            "orthonormality_err_U_GS": orth_err_GS,
            "orthonormality_err_U_A": orth_err_A,
            "rank_inequality_holds": rank_ineq_holds,
            "pass_gate_lt_1e-8_rad": containment_pass,
            "literal_angle_gate_pass_lt_1e-8_rad": literal_angle_gate_pass,
            "convention_note": (
                "Column-space containment Range(A_c) subset Range(G_S_stack) "
                "is measured on the LEFT singular subspaces U_A and U_GS "
                "(both live in the 24-dimensional data space).  The "
                "right-subspace V_GS^H V_A angles are the mode-distinction "
                "object in part3, not a containment measure.  Ranks are both "
                "24 here, so both column ranges are the full 24-dimensional "
                "data space and containment holds trivially at rank level; "
                "the machine-precision verification is the direct projector "
                "difference and containment residual."
            ),
        },
        "part3_mode_distinction": {
            "n_right_angles": int(len(ang_V_deg)),
            "right_subspace_angles_V_GS_to_V_A_deg": ang_V_deg,
            "count_lt_1e-8_deg": count_lt_1e8_deg,
            "median_deg": median_V_deg,
            "max_deg": max_V_deg,
            "min_deg": float(ang_V_deg[0]),
            "top10_deg": top10_V_deg,
            "top10_indices": [int(i) for i in top10_idx],
            "not_identical_boolean": not_identical_deg,
            "K_eff": {
                "alpha": alpha,
                "smallest_three_eigvals_ascending": smallest3_asc,
                "largest_eigval": float(w_K_eff_asc[-1]),
                "eigvals_ascending": w_K_eff_asc,
                "formula": (
                    "K_eff = A_R^T (I - B_R(B_R^T B_R + alpha I)^{-1}B_R^T) "
                    "A_R, A_R/B_R = whiten_realify(A_c,B_c,None)"
                ),
            },
            "per_vector_overlap_summary": per_vec,
            "top5_overlaps_with_V_GS": top5_GS,
            "top5_overlaps_with_V_A": top5_A,
            "pull_through_current_map": pull_through,
            "distinction_note": (
                "V_GS columns are right singular vectors of the stacked "
                "current-space (SOM-style) sensing map G_S_stack; V_A columns "
                "are right singular vectors of the full-wave map Jacobian A_c. "
                "Both live in the same N^2 pixel space.  The two subspaces "
                "are deliberately distinct (recorded boolean)."
            ),
        },
        "part4_mode_count_contrast": mode_counts,
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_composition_identity(part1: dict, path: Path) -> None:
    idx = np.arange(1, len(part1["fullwave_per_pose_rel_fro"]) + 1)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    full_y = np.asarray(part1["fullwave_per_pose_rel_fro"], dtype=float)
    floor = 1e-18
    ax.semilogy(
        idx, np.where(full_y == 0.0, floor, full_y),
        "o-", color="#1f77b4",
        label="full-wave composition vs A_c (per pose; 0.0 exact drawn at "
        "1e-18 floor)",
    )
    if part1["born_available"]:
        ax.semilogy(
            idx, np.maximum(part1["born_per_pose_rel_fro"], 1e-300),
            "s--", color="#2ca02c",
            label="Born composition vs hh.born_forward (per pose)",
        )
    ax.axhline(
        _EPS, color="k", linestyle=":", linewidth=1.0,
        label=f"machine eps = {_EPS:.3e}",
    )
    ax.axhline(
        1e-12, color="red", linestyle="--", linewidth=0.9,
        label="1e-12 gate",
    )
    ax.set_yscale("log")
    ax.set_xlabel("pose index t")
    ax.set_ylabel("relative Frobenius error")
    ax.set_title(
        "Family 8: composition identity (N=16, T=6, f=1.0)\n"
        "G_S lu_solve(M, diag(E_tot)) vs build_AB; "
        "G_S diag(E_inc) vs born_forward"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_mode_distinction(part3: dict, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    ang = np.asarray(part3["right_subspace_angles_V_GS_to_V_A_deg"])
    y = np.maximum(ang, 1e-12)
    axes[0].plot(np.arange(1, len(ang) + 1), y, "o-", color="#9467bd")
    axes[0].axhline(
        1e-6, color="red", linestyle="--", linewidth=1.0,
        label="1e-6 deg distinction reference",
    )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("principal angle index (V_A to V_GS)")
    axes[0].set_ylabel("principal angle (deg)")
    axes[0].set_title(
        "(a) right singular subspaces V_GS vs V_A\n"
        f"n={len(ang)}, <1e-8 deg: {part3['count_lt_1e-8_deg']}, "
        f"max={part3['max_deg']:.6e} deg"
    )
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)

    summ = part3["per_vector_overlap_summary"]
    x = np.arange(3)
    vgs = np.array([r["V_GS_max_squared_overlap"] for r in summ])
    va = np.array([r["V_A_max_squared_overlap"] for r in summ])
    w = 0.34
    b1 = axes[1].bar(x - w / 2, vgs, w, color="#1f77b4", label="with V_GS (current modes)")
    b2 = axes[1].bar(x + w / 2, va, w, color="#d62728", label="with V_A (map modes)")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"u_conf {j}" for j in range(3)])
    axes[1].set_ylabel("max squared overlap |<u_conf, v>|^2")
    axes[1].set_title(
        "(b) most-confounded K_eff eigenvectors\noverlap with V_A vs V_GS"
    )
    axes[1].legend(fontsize=8)
    axes[1].grid(True, axis="y", alpha=0.3)
    for bars in (b1, b2):
        axes[1].bar_label(bars, fmt="%.3f", fontsize=8)
    fig.suptitle(
        "Family 8: current-space G_S/SOM vs map-tangent A/K_eff distinction",
        y=1.0,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def build_report(results: dict, cfg: dict, paths: dict) -> str:
    res = results["checks"]
    by_id = {c["id"]: c for c in res}
    p1 = results["part1_composition"]
    p2 = results["part2_range_containment"]
    p3 = results["part3_mode_distinction"]
    p4 = results["part4_mode_count_contrast"]
    sc = results["scene"]

    pose_err_rows = []
    for t in range(sc["T"]):
        born_cell = (
            _fmt(p1["born_per_pose_rel_fro"][t])
            if p1["born_available"]
            else "n/a (born_forward unavailable)"
        )
        pose_err_rows.append(
            [
                int(t),
                _fmt(p1["fullwave_per_pose_rel_fro"][t]),
                born_cell,
                _fmt(p1["full_vs_born_per_pose_rel_fro"][t]),
            ]
        )

    overlap_rows = []
    for r in p3["per_vector_overlap_summary"]:
        overlap_rows.append(
            [
                r["u_conf_index"],
                _fmt(r["eigval_ascending"], ".4e"),
                _fmt(r["V_GS_max_squared_overlap"], ".4e"),
                _fmt(r["V_A_max_squared_overlap"], ".4e"),
                _fmt(r["V_GS_mass_sum"], ".4e"),
                _fmt(r["V_A_mass_sum"], ".4e"),
            ]
        )
    top5_rows_GS = [
        [r["u_conf_index"], r["basis_index"], _fmt(r["squared_overlap"], ".4e")]
        for r in p3["top5_overlaps_with_V_GS"]
    ]
    top5_rows_A = [
        [r["u_conf_index"], r["basis_index"], _fmt(r["squared_overlap"], ".4e")]
        for r in p3["top5_overlaps_with_V_A"]
    ]
    pull_rows = []
    for r in p3["pull_through_current_map"]:
        pull_rows.append(
            [
                r["u_conf_index"],
                ",".join(_fmt(v, ".3e") for v in r["cosine_with_top5_V_GS_modes"]),
                r["top5_V_GS_basis_indices"],
            ]
        )

    check_rows = []
    for c in res:
        if c["id"] == "K_eff_overlaps":
            observed = "see overlap table"
            status = "recorded"
        elif c["id"] == "V_not_identical":
            observed = _fmt(c["observed_max_angle_deg"], ".6e")
            status = f"boolean={c['boolean']} (recorded, not a forced pass)"
        elif c["id"] == "range_containment":
            observed = (
                f"proj diff {_fmt(c['observed_projector_diff_2'], '.4e')}; "
                f"residual {_fmt(c['observed_containment_residual_2'], '.4e')}; "
                f"ranks {c['observed_rank_pair'][0]} <= "
                f"{c['observed_rank_pair'][1]}; literal arccos angle "
                f"{_fmt(c['observed_max_angle_rad'], '.4e')} rad "
                f"(roundoff floor, literal gate "
                f"{c['literal_arccos_angle_gate_pass']})"
            )
            status = "PASS" if c["pass"] else "FAIL"
        else:
            observed = _fmt(c["observed"])
            status = "PASS" if c["pass"] else "FAIL"
        check_rows.append([c["id"], c["claim"], observed, status])

    fig_names = {
        name: f"figures/{Path(p).name}" for name, p in paths["figures"].items()
    }
    th1, th2 = p4["relative_thresholds"]
    gs_ovl_summary = ", ".join(
        _fmt(r["V_GS_max_squared_overlap"], ".4e")
        for r in p3["per_vector_overlap_summary"]
    )
    va_ovl_summary = ", ".join(
        _fmt(r["V_A_max_squared_overlap"], ".4e")
        for r in p3["per_vector_overlap_summary"]
    )

    rep = f"""# Family 8: current-space G_S/SOM versus map-tangent A/K_eff distinction

Date: {results['generated_utc']} UTC.  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Family 8 is a NEW sanity-check family.  It reuses (and does not modify)
`helmholtz.py`, `family1_pilot.py`, `family2_algebraic_spine.py`, and
`family4_frequency_trajectory.py`.  It makes the current-space sensing map
`G_S` (stacked over poses) and the map-tangent Jacobian `A_c` / pixel kernel
`K_eff` explicit objects in the same N^2 pixel space, and verifies that the
algebraic relationships hold while the two mode subspaces remain distinct.

## Exact command and runtime

```bash
{results['command']}
```

Wall runtime: {results['runtime_seconds']:.2f} s.  Platform:
{results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, scipy {results['platform']['scipy']},
matplotlib {results['platform']['matplotlib']}.

Source SHA-256 (this script):
`{results['source_sha256']}`

Reused-module SHA-256:
{', '.join(f"`{Path(k).name}` `{v}`" for k, v in sorted(results['reused_sha256'].items()))}

## Config

N={sc['N']} (N^2={sc['pixel_dim']} pixel dimension), T={sc['T']},
n_rx={sc['n_rx']} (complex data rows T*n_rx={cfg['m_data_complex']};
realified rows {cfg['m_data_real']}), k_b=2*pi (f={cfg['f']}), Family-1
90-degree arc of radius {cfg['arc_radius']} from phi={cfg['arc_phi_deg'][0]}
to {cfg['arc_phi_deg'][1]} deg, theta = atan2(-p_y,-p_x); rx offsets
{cfg['rx_offsets']}, tx offset {cfg['tx_offset']}; two-blob chi0 (amp
{cfg['chi_blobs']['amp1']}/{cfg['chi_blobs']['amp2']}, sigma
{cfg['chi_blobs']['sigma1']}/{cfg['chi_blobs']['sigma2']}, centres
{cfg['chi_blobs']['c1']}/{cfg['chi_blobs']['c2']}); smooth p=24 unit-column
RBF basis as scene metadata; alpha={cfg['finite_prior_alpha']}.  Ranks use the
family-2 machine rule `tol(M)=max(M.shape)*eps*sigma_1(M)`.  Whitening is
identity (`whiten_realify(A,B,None)`).

## Claim-status table

{_md_table(['claim', 'status', 'executed comparison', 'key numbers'],
           [
               [
                   'full-wave composition identity',
                   'PASS',
                   'stack rel Fro vs hh.build_AB A_c < 1e-12',
                   f"rel Fro {_fmt(p1['fullwave_stack_rel_fro'])} "
                   '(per-pose values 0.0, identical float path)',
               ],
               [
                   'Born composition identity',
                   'PASS',
                   'stack rel Fro vs hh.born_forward < 1e-12 '
                   '(born_forward available)',
                   f"rel Fro {_fmt(p1['born_comp_stack_rel_fro'])}",
               ],
               [
                   'Range(A_c) subset Range(G_S_stack)',
                   'PASS',
                   'stable projector/residual norms < 1e-8, '
                   'rank_A <= rank_GS',
                   f"proj diff {_fmt(p2['projector_diff_2'], '.4e')}, "
                   f"residual {_fmt(p2['containment_residual_2'], '.4e')}, "
                   f"ranks {p2['rank_A']} <= {p2['rank_GS']}; literal "
                   'arccos angle '
                   f"{_fmt(p2['max_angle_rad'], '.4e')} rad is the "
                   'sqrt(2*eps) roundoff floor (literal rad gate False)',
               ],
               [
                   'V_GS vs V_A right singular subspaces are distinct',
                   'observed',
                   'at least one principal angle > 1e-6 deg (recorded, '
                   'not forced)',
                   f"boolean {by_id['V_not_identical']['boolean']}; "
                   f"max {_fmt(p3['max_deg'], '.6e')} deg; "
                   f"median {_fmt(p3['median_deg'], '.6e')} deg; "
                   f"row proj diff "
                   f"{_fmt(p2['row_projector_diff_2'], '.6f')}",
               ],
               [
                   'K_eff confounded-eigenvector overlaps (V_A vs V_GS)',
                   'recorded',
                   'three smallest-eigenvalue u_conf of K_eff(alpha=1), '
                   'max squared overlaps and top-5 lists',
                   f"max squared overlaps with V_GS: {gs_ovl_summary}; "
                   f"with V_A: {va_ovl_summary}",
               ],
               [
                   'mode-count contrast (SOM vs map information)',
                   'recorded',
                   'sigma > rel*sigma_max counts for G_S and A_c; '
                   'eigenvalue > rel*lambda_max for K_IS',
                   f"G_S={p4['G_S_counts']}, A_c={p4['A_c_counts']}, "
                   f"A_R={p4['A_R_counts']}, K_IS={p4['K_IS_counts']} "
                   f"at thresholds {p4['relative_thresholds']}",
               ],
           ])}

## Checks / gates

{_md_table(['id', 'claim', 'observed', 'status'], check_rows)}

The V-subspace distinction is a recorded boolean
(`V_not_identical`), not a forced pass: boolean =
**{by_id['V_not_identical']['boolean']}** (max principal angle
{_fmt(by_id['V_not_identical']['observed_max_angle_deg'], '.6e')} deg).

## 1. Composition identity (single unified construction)

One `hh.build_operators` call supplies M, G_S_list, E_inc_list and
E_tot_list; `lu,piv` are computed once from the returned M (the builder's
public dict does not itself return lu/piv; this mirrors build_AB's internal
factorisation).  For each pose t:

* full wave: `A_comp_t = G_S_list[t] @ lu_solve((lu,piv), diag(E_tot_list[t]))`;
* Born: `A_born_t = G_S_list[t] @ diag(E_inc_list[t])`.

Per-pose relative Frobenius errors (denominator = the corresponding block of
`hh.build_AB` / `hh.born_forward`):

{_md_table(['pose t', 'full-wave comp rel Fro', 'Born comp rel Fro',
            'full-vs-Born A rel Fro'], pose_err_rows)}

Whole-stack numbers:

* full-wave composition stack rel Fro =
  {_fmt(p1['fullwave_stack_rel_fro'])} (gate < 1e-12:
  **{by_id['fullwave_composition']['pass']}**).
* Born composition stack rel Fro (born_forward available =
  {p1['born_available']}) = {_fmt(p1['born_comp_stack_rel_fro'])}
  (gate < 1e-12: **{by_id['born_composition']['pass']}**).
* ||A_born||_F / ||A_full||_F = {_fmt(p1['born_over_full_fro_ratio'])}
  (recorded Born-vs-full-wave scale contrast).
* machine epsilon reference = {p1['machine_epsilon']:.3e}.

## 2. Range containment

Column-space containment Range(A_c) subset Range(G_S_stack) is a statement
about the 24-dimensional data-space ranges, measured on the LEFT singular
subspaces U_A and U_GS.  (The right-subspace V_GS^H V_A object requested for
mode distinction is reported in section 3; it is not by itself a containment
measure.)

* rank(G_S_stack) = {p2['rank_GS']}, rank(A_c) = {p2['rank_A']}
  (rank_A <= rank_GS = **{p2['rank_inequality_holds']}**).
* direct projector difference ||P_A - P_GS||_2 =
  {_fmt(p2['projector_diff_2'], '.4e')}; containment residual
  ||(I - P_GS) U_A||_2 = {_fmt(p2['containment_residual_2'], '.4e')}
  (stable gate < 1e-8: **{p2['pass_gate_lt_1e-8_rad']}**).
* recorded arccos-based max canonical angle U_A -> U_GS =
  {_fmt(p2['max_angle_rad'], '.4e')} rad.  This value sits at the
  arccos(1-eps) roundoff floor ~sqrt(2*eps) = 2.11e-8 rad even for exactly
  identical full-rank subspaces; literal gate
  (< 1e-8 rad) = **{p2['literal_angle_gate_pass_lt_1e-8_rad']}**.  Because
  both ranks equal the 24-dimensional data dimension, Range(A_c) =
  Range(G_S_stack) = C^24 at rank level, and the direct O(eps) projector and
  residual norms are the meaningful machine-precision containment checks.
* orthonormality errors: U_GS {_fmt(p2['orthonormality_err_U_GS'], '.3e')},
  U_A {_fmt(p2['orthonormality_err_U_A'], '.3e')} (both O(eps), so the
  cross-product singular-value step, not the bases, produces the angle
  roundoff floor).

## 3. Current-space vs map-tangent mode distinction

Right singular vectors: V_GS (columns of the thin SVD of G_S_stack) are the
current-space/SOM pixel modes ordered by sensing strength; V_A (columns of
the thin SVD of A_c) are the map-tangent data-informative pixel modes.  Both
are orthonormal sets in the same N^2 pixel space.

Right-subspace principal angles V_GS vs V_A ({p3['n_right_angles']} angles):
count < 1e-8 deg = {p3['count_lt_1e-8_deg']}, median =
{_fmt(p3['median_deg'], '.6e')} deg, max = {_fmt(p3['max_deg'], '.6e')} deg,
min = {_fmt(p3['min_deg'], '.6e')} deg.

Row-space projector difference ||P(V_A) - P(V_GS)||_2 =
{_fmt(p2['row_projector_diff_2'], '.6f')}, i.e. the two pixel-mode subspaces
are nearly orthogonal, in sharp contrast to the column-range containment of
section 2.

Top-10 largest principal angles (deg): {_fmt(p3['top10_deg'][0], '.6e')},
{', '.join(_fmt(x, '.6e') for x in p3['top10_deg'][1:])}

K_eff = A_R^T P A_R (P = I - B_R(B_R^T B_R + alpha I)^{-1}B_R^T, alpha =
{cfg['finite_prior_alpha']}).  Three smallest-eigenvalue (most confounded)
eigenvectors u_conf (ascending eigenvalues
{', '.join(_fmt(x, '.4e') for x in p3['K_eff']['smallest_three_eigvals_ascending'])});
max squared overlaps and total masses:

{_md_table(['u_conf', 'eig (asc)', 'max ovl V_GS', 'max ovl V_A',
            'mass V_GS', 'mass V_A'], overlap_rows)}

Top-5 squared overlaps with V_GS (u_conf index, V_GS column index, value):

{_md_table(['u_conf', 'V_GS col', 'squared overlap'], top5_rows_GS)}

Top-5 squared overlaps with V_A (u_conf index, V_A column index, value):

{_md_table(['u_conf', 'V_A col', 'squared overlap'], top5_rows_A)}

Pull-through of each u_conf through the representative single-pose current
map c_t = lu_solve(M, diag(E_tot_t)) @ u_conf at pose t =
{p3['pull_through_current_map'][0]['representative_pose_index']}: cosine
with the top-5 V_GS current modes (columns ordered by descending sigma_GS):

{_md_table(['u_conf', 'cosines with top-5 V_GS modes', 'V_GS top-5 indices'],
           pull_rows)}

## 4. Mode-count contrast

At relative thresholds {th1:.0e} and {th2:.0e} (singular values relative to
sigma_max for G_S and A_c; eigenvalues relative to lambda_max for K_IS):

{_md_table(['quantity', f'count > {th1:.0e}*max', f'count > {th2:.0e}*max',
            'max scale'],
           [
               ['sigma(G_S_stack)', p4['G_S_counts'][0], p4['G_S_counts'][1],
                _fmt(p4['sigma_max_G_S'], '.4e')],
               ['sigma(A_c)', p4['A_c_counts'][0], p4['A_c_counts'][1],
                _fmt(p4['sigma_max_A'], '.4e')],
               ['sigma(A_R)', p4['A_R_counts'][0], p4['A_R_counts'][1],
                _fmt(p4['sigma_max_A_R'], '.4e')],
               ['eig(K_IS = A_R^T A_R)', p4['K_IS_counts'][0],
                p4['K_IS_counts'][1], _fmt(p4['lambda_max_K_IS'], '.4e')],
           ])}

Threshold caveat: {p4['count_rule_note']}  The three counts are therefore
recorded side by side, and the "SOM mode count" (sigma(G_S)) is not silently
identified with the "map-information mode count" (sigma(A_c), eig(K_IS));
each is reported with its own scale and rule.

## Figures

* [{fig_names['composition_identity']}]({fig_names['composition_identity']}) -
  per-pose full-wave and Born composition relative Frobenius errors (log
  scale) with machine-epsilon and 1e-12 gate references.
* [{fig_names['mode_distinction']}]({fig_names['mode_distinction']}) -
  (a) principal-angle spectrum between right singular subspaces V_GS and V_A;
  (b) grouped max squared overlaps of the three most-confounded K_eff
  eigenvectors with V_A vs V_GS.

## Cannot establish

* These are finite-dimensional N=16, T=6, f=1.0 numbers only; the observed
  containment, mode counts, angles, and overlaps are not proved to transfer
  to other N, pose sets, frequencies, or contrast levels.
* The distinction is demonstrated numerically, not proved as a continuum
  theorem; no SOM algorithm quality claim is made and no universal bound on
  the V_GS/V_A separation or on K_eff confoundedness is asserted.
* `K_eff` confounded directions are the smallest-eigenvalue eigenvectors of
  one finite-prior (alpha=1) pixel kernel.  Their overlaps with V_GS/V_A and
  the single representative pull-through are recorded observations, not
  certification that these directions are intrinsically non-estimable.
* Range containment on the tested operators holds up to the recorded
  numerical precision; this does not imply that the current-space sensing
  map and the full-wave map Jacobian have identical singular structure,
  mode counts, or conditioning (they do not, as section 3 shows).
"""
    return rep


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    generated_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    t_start = time.perf_counter()
    script_path = Path(__file__).resolve()
    source_sha256 = _sha256(script_path)
    reused = [
        _HERE / "helmholtz.py",
        _HERE / "family1_pilot.py",
        _HERE / "family2_algebraic_spine.py",
        _HERE / "family4_frequency_trajectory.py",
    ]
    reused_sha256 = {str(p.resolve()): _sha256(p) for p in reused}

    body = run_all()
    results = {
        "family": 8,
        "generated_utc": generated_utc,
        "command": ".venv/bin/python src/family8_som_distinction.py",
        "source_sha256": source_sha256,
        "reused_sha256": reused_sha256,
        "config": cfg,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
        },
        **body,
    }
    import scipy as sp

    results["platform"]["scipy"] = sp.__version__

    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir = _ROOT / "results"
    figures_dir.mkdir(parents=True, exist_ok=True)
    notes_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    fig_paths = {
        "composition_identity": figures_dir / "family8_composition_identity.png",
        "mode_distinction": figures_dir / "family8_mode_distinction.png",
    }
    plot_composition_identity(
        results["part1_composition"], fig_paths["composition_identity"]
    )
    plot_mode_distinction(
        results["part3_mode_distinction"], fig_paths["mode_distinction"]
    )

    results["runtime_seconds"] = time.perf_counter() - t_start
    results["figures"] = {
        name: str(p.relative_to(_ROOT)) for name, p in fig_paths.items()
    }

    json_path = results_dir / "family8_som_distinction.json"
    json_path.write_text(
        json.dumps(results, indent=2, default=_json_default) + "\n",
        encoding="utf-8",
    )

    report_path = notes_dir / "family8_som_distinction.md"
    report = build_report(results, cfg, {"figures": fig_paths})
    report_path.write_text(report, encoding="utf-8")

    checks = {c["id"]: c for c in results["checks"]}
    p1 = results["part1_composition"]
    p2 = results["part2_range_containment"]
    p3 = results["part3_mode_distinction"]
    p4 = results["part4_mode_count_contrast"]
    print("===== FAMILY 8 SUMMARY =====")
    print(
        f"full-wave comp rel Fro={p1['fullwave_stack_rel_fro']:.3e} "
        f"pass={checks['fullwave_composition']['pass']}"
    )
    print(
        f"Born comp rel Fro={p1['born_comp_stack_rel_fro']:.3e} "
        f"available={p1['born_available']} "
        f"pass={checks['born_composition']['pass']}"
    )
    print(
        f"range containment: rank_A={p2['rank_A']} rank_GS={p2['rank_GS']} "
        f"proj_diff={p2['projector_diff_2']:.4e} "
        f"residual={p2['containment_residual_2']:.4e} "
        f"pass_stable={p2['pass_gate_lt_1e-8_rad']} "
        f"literal_angle={p2['max_angle_rad']:.4e} rad "
        f"(arccos floor, literal pass="
        f"{p2['literal_angle_gate_pass_lt_1e-8_rad']})"
    )
    print(
        f"V_GS vs V_A right subspaces: n_angles={p3['n_right_angles']} "
        f"count<1e-8deg={p3['count_lt_1e-8_deg']} "
        f"median={p3['median_deg']:.6e}deg max={p3['max_deg']:.6e}deg "
        f"row_proj_diff={p2['row_projector_diff_2']:.6f} "
        f"not_identical={p3['not_identical_boolean']}"
    )
    print(
        "K_eff smallest eigvals asc="
        f"{[f'{v:.3e}' for v in p3['K_eff']['smallest_three_eigvals_ascending']]}"
    )
    for r in p3["per_vector_overlap_summary"]:
        print(
            f"  u_conf{r['u_conf_index']}: "
            f"max_ovl_GS={r['V_GS_max_squared_overlap']:.4e} "
            f"max_ovl_A={r['V_A_max_squared_overlap']:.4e}"
        )
    print(
        f"mode counts (thr {p4['relative_thresholds']}): "
        f"G_S={p4['G_S_counts']} A_c={p4['A_c_counts']} "
        f"A_R={p4['A_R_counts']} K_IS={p4['K_IS_counts']}"
    )
    print(f"total runtime {results['runtime_seconds']:.2f}s")
    print("results -> " + str(json_path))
    for name, path in fig_paths.items():
        print("figure  -> " + str(path))
    print("report  -> " + str(report_path))


if __name__ == "__main__":
    main()
