"""Family 4: frequency stacking and trajectory geometry.

Family 1 established the finite-difference-validated forward/Jacobian pair
(A_R, B_R).  Family 2 established the algebraic spine for a single
frequency/pose block and the smooth p=24 RBF basis.  Family 3/3b checked
gauge, Born, and corrected variants of those claims.  This family does not
rerun those checks; it reuses the corrected builders and asks two new
questions on the same discrete finite-dimensional model:

  A. PSD monotonicity of stacking distinct frequencies (K_SLAM and the
     finite-prior K_eff for prefixes F=1,2,3).
  B. Duplicate controls: c-scaled repeated single-frequency data are
     invariant for the no-prior and jointly-scaled-prior problems but shift
     the retention spectrum under a fixed finite prior.
  C. Shared-pose-compensation: adding a distinct frequency to the pose block
     must move the three most confounded single-frequency directions away
     from the pose null space, while duplicating the same frequency must not.
  D. Equal-budget trajectory comparison: four trajectories with equal path
     length and measurement count, reported without a forced pass.

All numbers are for the specific whitened/realified dense linear algebra of
the discrete N=16 model; no continuum-limit or universal-trajectory claim is
made (see notes/family4_report.md).
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

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CONFIG = {
    "N": 16,
    "T": 6,
    "n_rx": 4,
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],
    "q_pose": 18,
    "m_per_frequency": 48,  # 2 * T * n_rx after whiten_realify
    "k_b_rule": "k_b(f) = 2*pi*f",
    "whitening_convention": (
        "W = None (identity noise): whiten_realify(A,B,None) returns "
        "sqrt(2)*[Re; Im] row stacks; raw identity-noise stacks throughout. "
        "Frequencies and trajectories are NOT block-energy or SNR normalized."
    ),
    "block_normalized_control": (
        "deferred and declared: raw identity-noise stacking is the tested "
        "metric; a block-normalized/SNR-matched control is not silently "
        "conflated with it"
    ),
    "f_scales": [1.0, 1.4, 1.8],
    "distinct_prefix_frequencies": [[1.0], [1.0, 1.4], [1.0, 1.4, 1.8]],
    "duplicate_c": 2.0,
    "finite_prior_alpha": 1.0,
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "stack_pose_trajectory": (
        "family1 arc: T=6, 90-degree arc at radius 1.6, "
        "theta = atan2(-p_y,-p_x); same pose set for every frequency"
    ),
    "rx_offsets": [
        [-0.06, 0.0],
        [0.06, 0.0],
        [0.0, -0.06],
        [0.0, 0.06],
    ],
    "tx_offset": [0.0, 0.0],
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
        "note": "identical construction/config to family2 and family3",
    },
    "trajectories": {
        "T": 6,
        "k_b": 2.0 * np.pi,
        "total_path_length_L": np.pi * 1.6,
        "theta_rule": "theta = atan2(-p_y,-p_x): body +x axis points at origin",
        "straight": "p = (x, 1.6), x = linspace(-L/2, L/2, T)",
        "arc90": "R = 2*L/pi, phi = linspace(-45, 45, T) degrees",
        "arc180": "R = L/pi, phi = linspace(-90, 90, T) degrees",
        "circle360": "R = L/(2*pi), phi = linspace(0, 360, T, endpoint=False)",
        "expected_order_retained": ["circle360", "arc180", "arc90", "straight"],
        "hypothesis_rule": (
            "retained_mass circle360 > arc180 > arc90 > straight and "
            "confusable_mass circle360 < arc180 < arc90 < straight; "
            "violations recorded as hypothesis_holds = false (no forced pass)"
        ),
    },
    "tolerances": {
        "A_min_eig_rule": "min_eig(D) >= -1e-10 * max(1, ||larger||_2)",
        "B_rel_fro_no_prior": 1e-10,
        "B_max_abs_rho_no_prior": 1e-10,
        "B_fixed_prior_nonzero_example": 1e-3,
        "B_joint_rel_fro": 1e-10,
        "B_joint_max_abs_rho": 1e-10,
        "C_rho_dup_movement": 1e-10,
        "C_rho_dist_movement": 1e-4,
        "C_required_directions": 2,
    },
    "seeds": [],
    "randomness_note": "deterministic dense linear algebra; no RNG used",
}


# ---------------------------------------------------------------------------
# Reusable construction (reads directly from the earlier families)
# ---------------------------------------------------------------------------

def base_scene(cfg: dict) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """N-grid, chi0, and the shared unit-column smooth basis S."""
    points, h = hh.make_grid(cfg["N"])
    chi0 = family1.make_chi0(points, cfg)
    S = family2.build_smooth_basis(points, cfg["smooth_basis"])
    return points, chi0, h, S


def build_AB_whitened(
    chi0: np.ndarray,
    poses: np.ndarray,
    cfg: dict,
    k_b: float,
) -> tuple[np.ndarray, np.ndarray]:
    """hh.build_AB then hh.whiten_realify(None); real (A_pix_R, B_R)."""
    rx_offsets = np.asarray(cfg["rx_offsets"], dtype=float)
    tx_offset = np.asarray(cfg["tx_offset"], dtype=float)
    A_c, B_c, _, _ = hh.build_AB(
        chi0, poses, rx_offsets, tx_offset, cfg["N"], k_b
    )
    return hh.whiten_realify(A_c, B_c, None)


def build_smooth_block(
    chi0: np.ndarray,
    poses: np.ndarray,
    S: np.ndarray,
    cfg: dict,
    f: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A_smooth (48 x p), B_R (48 x 18), and A_pix_R for frequency f."""
    A_pix_R, B_R = build_AB_whitened(chi0, poses, cfg, 2.0 * np.pi * float(f))
    return A_pix_R @ S, B_R, A_pix_R


def _sym(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def range_basis(B: np.ndarray):
    """Orthonormal Z for Range(B), machine-rank tolerance (family 2 rule)."""
    return family2.range_basis(B)


def range_projection(B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(Z, P_perp) with P_perp = I - Z Z^T."""
    rB, Z, _, _ = range_basis(B)
    rows = B.shape[0]
    return Z, np.eye(rows, dtype=float) - Z @ Z.T


def K_SLAM(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """K_SLAM = A^T (I - Z Z^T) A for orthonormal Range(B) basis Z."""
    _, P = range_projection(B)
    return _sym(A.T @ (P @ A))


def K_eff(A: np.ndarray, B: np.ndarray, alpha: float) -> np.ndarray:
    """K_eff(alpha) = A^T A - A^T B (B^T B + alpha I)^{-1} B^T A."""
    q = B.shape[1]
    C = B.T @ B + float(alpha) * np.eye(q, dtype=float)
    return _sym(A.T @ A - A.T @ B @ np.linalg.solve(C, B.T @ A))


def shrinkage_operator(B: np.ndarray, alpha: float) -> np.ndarray:
    """W = I - B (B^T B + alpha I)^{-1} B^T (symmetric data-space factor)."""
    return family2.weighted_shrinkage_operator(B, float(alpha))


def retention_spectrum(A: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Descending retention rho of K = A^T W A w.r.t. K_IS = A^T A.

    Equivalently the eigenvalues of diag(1/s_A) V_A^T K V_A diag(1/s_A) for
    the thin SVD A = Q_A diag(s_A) V_A^T, computed in the Q_A^T W Q_A form
    used by family2 so that 1/s_A^2 is not amplified in near-null directions.
    """
    dA = family2.thin_decomposition(A)
    return family2.retention_spectrum(dA["Q"], W)


def generalized_eigen_directions(
    A: np.ndarray, W: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Ascending generalized spectrum and directions of K = A^T W A in K_IS.

    Returns (rho_asc, U, R_op) where R_op = Q_A^T W Q_A has eigenvectors w,
    and the parameter-space direction u = V_A diag(1/s_A) w satisfies
    u^T K_IS u = 1 and u^T K u = rho.
    """
    dA = family2.thin_decomposition(A)
    Q, s, V = dA["Q"], dA["s"], dA["V"]
    R_op = _sym(Q.T @ (W @ Q))
    rho, Wv = np.linalg.eigh(R_op)              # ascending
    U = V @ np.diag(1.0 / s) @ Wv               # columns follow rho ascending
    return rho, U, R_op


def rayleigh_quotient(A: np.ndarray, W: np.ndarray, u: np.ndarray) -> float:
    """u^T K_SLAM u / u^T K_IS u = ||P W A u||^2 / ||A u||^2."""
    Au = A @ u
    num = float(np.dot(Au, W @ Au))
    den = float(np.dot(Au, Au))
    return num / den


def matrix_norm_stats(M: np.ndarray) -> dict:
    """Frobenius/spectral norms and per-row 2-norm statistics."""
    row_norms = np.linalg.norm(np.asarray(M), axis=1)
    return {
        "frobenius": float(np.linalg.norm(M, ord="fro")),
        "spectral_2": float(np.linalg.norm(M, ord=2)),
        "row_norm_mean": float(np.mean(row_norms)),
        "row_norm_max": float(np.max(row_norms)) if len(row_norms) else None,
    }


def shared_z_residual(A: np.ndarray, B: np.ndarray, u: np.ndarray) -> float:
    """min_z ||A u - B z||_2 / ||A u||_2."""
    Au = A @ u
    z, _, _, _ = np.linalg.lstsq(B, Au, rcond=None)
    return float(np.linalg.norm(Au - B @ z) / np.linalg.norm(Au))


def trajectory_poses(name: str, cfg: dict) -> np.ndarray:
    """Pose array (T,3) for one equal-path-length trajectory."""
    tc = cfg["trajectories"]
    L = float(tc["total_path_length_L"])
    T = int(tc["T"])
    if name == "straight":
        x = np.linspace(-L / 2.0, L / 2.0, T)
        p = np.column_stack([x, np.full(T, 1.6)])
    elif name == "arc90":
        R = 2.0 * L / np.pi
        phi = np.deg2rad(np.linspace(-45.0, 45.0, T))
        p = R * np.column_stack([np.cos(phi), np.sin(phi)])
    elif name == "arc180":
        R = L / np.pi
        phi = np.deg2rad(np.linspace(-90.0, 90.0, T))
        p = R * np.column_stack([np.cos(phi), np.sin(phi)])
    elif name == "circle360":
        R = L / (2.0 * np.pi)
        phi = np.deg2rad(np.linspace(0.0, 360.0, T, endpoint=False))
        p = R * np.column_stack([np.cos(phi), np.sin(phi)])
    else:
        raise ValueError(f"unknown trajectory {name!r}")
    theta = np.arctan2(-p[:, 1], -p[:, 0])
    return np.column_stack([p, theta])


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_A_monotonicity(blocks: dict) -> dict:
    """PSD monotonicity of prefix stacks with distinct frequencies."""
    alpha = CONFIG["finite_prior_alpha"]
    prefixes = CONFIG["distinct_prefix_frequencies"]
    KSL = {}
    Keff = {}
    spectra = {}
    for F, freqs in enumerate(prefixes, start=1):
        A_st = np.vstack([blocks[f]["A"] for f in freqs])
        B_st = np.vstack([blocks[f]["B"] for f in freqs])
        KSL[F] = K_SLAM(A_st, B_st)
        Keff[F] = K_eff(A_st, B_st, alpha)
        spectra[f"K_SLAM_F{F}"] = [
            float(x) for x in np.sort(np.linalg.eigvalsh(KSL[F]))[::-1]
        ]
        spectra[f"K_eff_F{F}"] = [
            float(x) for x in np.sort(np.linalg.eigvalsh(Keff[F]))[::-1]
        ]

    rows = []
    ok = True
    for F in (1, 2):
        for label, K in (("K_SLAM", KSL), ("K_eff", Keff)):
            D = K[F + 1] - K[F]
            min_eig = float(np.linalg.eigvalsh(D)[0])
            scale = float(max(1.0, np.linalg.norm(K[F + 1], ord=2)))
            gate = min_eig >= -1e-10 * scale
            ok &= gate
            rows.append(
                {
                    "stack_from": F,
                    "stack_to": F + 1,
                    "label": label,
                    "min_eig_difference": min_eig,
                    "scale_max_1_norm_larger_2": scale,
                    "tol": 1e-10 * scale,
                    "pass": bool(gate),
                }
            )
    return {
        "alpha_finite_prior": alpha,
        "prefix_frequencies": prefixes,
        "K_SLAM_eig_desc": {
            "F1": spectra["K_SLAM_F1"],
            "F2": spectra["K_SLAM_F2"],
            "F3": spectra["K_SLAM_F3"],
        },
        "K_eff_eig_desc": {
            "F1": spectra["K_eff_F1"],
            "F2": spectra["K_eff_F2"],
            "F3": spectra["K_eff_F3"],
        },
        "rows": rows,
        "pass_gate": bool(ok),
        "note": (
            "K_SLAM(F) uses Range(B_st) projection; K_eff(F; alpha=1.0) uses "
            "the finite prior (B_st^T B_st + alpha I)^{-1}.  Prefix order is "
            "the stated frequency list [1.0], [1.0,1.4], [1.0,1.4,1.8].  "
            "This is ABSOLUTE information monotonicity of the K matrices "
            "(audit rule 2); normalized generalized retention spectra are "
            "not used as the stacking gate because K_IS also changes."
        ),
    }


def check_B_duplicates(blocks: dict) -> dict:
    """Duplicate single-frequency block: invariance of no-prior/joint cases."""
    c = CONFIG["duplicate_c"]
    alpha = CONFIG["finite_prior_alpha"]
    scale = 1.0 + c**2
    A1, B1 = blocks[1.0]["A"], blocks[1.0]["B"]
    Adup = np.vstack([A1, c * A1])
    Bdup = np.vstack([B1, c * B1])

    KIS1 = _sym(A1.T @ A1)
    KSL1 = K_SLAM(A1, B1)
    KISd = _sym(Adup.T @ Adup)
    KSLd = K_SLAM(Adup, Bdup)

    _, P_no1 = range_projection(B1)
    _, P_nod = range_projection(Bdup)
    rho_1 = retention_spectrum(A1, P_no1)
    rho_d = retention_spectrum(Adup, P_nod)

    rel_KIS = float(
        np.linalg.norm(KISd - scale * KIS1, ord="fro")
        / np.linalg.norm(KIS1, ord="fro")
    )
    rel_KSL = float(
        np.linalg.norm(KSLd - scale * KSL1, ord="fro")
        / max(np.linalg.norm(KSL1, ord="fro"), 1e-300)
    )
    no_prior_rho_diff = float(np.max(np.abs(rho_d - rho_1)))
    no_prior_pass = (
        rel_KIS < 1e-10 and rel_KSL < 1e-10 and no_prior_rho_diff < 1e-10
    )

    W_f1 = shrinkage_operator(B1, alpha)
    W_fd = shrinkage_operator(Bdup, alpha)
    Keff1 = _sym(A1.T @ (W_f1 @ A1))
    Keffd = _sym(Adup.T @ (W_fd @ Adup))
    rhoX_1 = retention_spectrum(A1, W_f1)
    rhoX_d = retention_spectrum(Adup, W_fd)
    fixed_rho_diff = float(np.max(np.abs(rhoX_d - rhoX_1)))
    fixed_rel_Keff = float(
        np.linalg.norm(Keffd - scale * Keff1, ord="fro")
        / max(np.linalg.norm(Keff1, ord="fro"), 1e-300)
    )

    W_j = shrinkage_operator(Bdup, scale * alpha)
    Keff_j = _sym(Adup.T @ (W_j @ Adup))
    rhoX_dj = retention_spectrum(Adup, W_j)
    joint_rel_Keff = float(
        np.linalg.norm(Keff_j - scale * Keff1, ord="fro")
        / max(np.linalg.norm(Keff1, ord="fro"), 1e-300)
    )
    joint_rho_diff = float(np.max(np.abs(rhoX_dj - rhoX_1)))
    joint_pass = joint_rel_Keff < 1e-10 and joint_rho_diff < 1e-10

    return {
        "duplicate_c": c,
        "scale_1_plus_c2": float(scale),
        "finite_prior_alpha": alpha,
        "K_IS_single": [float(x) for x in np.sort(np.linalg.eigvalsh(KIS1))[::-1]],
        "K_IS_dup": [float(x) for x in np.sort(np.linalg.eigvalsh(KISd))[::-1]],
        "no_prior": {
            "rel_F_KIS_dup_vs_scaled": rel_KIS,
            "rel_F_KSLAM_dup_vs_scaled": rel_KSL,
            "max_abs_rho_dup_minus_rho_single": no_prior_rho_diff,
            "rho_single_desc": [float(x) for x in rho_1],
            "rho_dup_desc": [float(x) for x in rho_d],
            "pass": bool(no_prior_pass),
        },
        "fixed_prior": {
            "max_abs_rho_X_dup_minus_rho_X_single": fixed_rho_diff,
            "rel_F_Keff_dup_vs_scaled_Keff_single": fixed_rel_Keff,
            "nonzero_threshold_example": CONFIG["tolerances"]["B_fixed_prior_nonzero_example"],
            "rho_X_nonzero_observed": bool(fixed_rho_diff > 1e-3),
            "rel_Keff_nonzero_observed": bool(fixed_rel_Keff > 1e-6),
            "rho_X_single_desc": [float(x) for x in rhoX_1],
            "rho_X_dup_desc": [float(x) for x in rhoX_d],
        },
        "joint_scaled_prior": {
            "joint_prior_matrix": f"(1+c^2)*alpha*I = {scale}*{alpha}*I",
            "rel_F_Keff_dup2_vs_scaled_Keff_single": joint_rel_Keff,
            "max_abs_rho_X_dup2_minus_rho_X_single": joint_rho_diff,
            "rho_X_dup2_desc": [float(x) for x in rhoX_dj],
            "pass": bool(joint_pass),
        },
        "pass_gate": bool(no_prior_pass and joint_pass),
        "note": (
            "Fixed-prior break is expected and recorded raw; it is not a "
            "failure of the duplicate controls."
        ),
    }


def check_C_shared_pose_compensation(blocks: dict) -> dict:
    """Distinct frequency breaks single-frequency pose confounding."""
    c = CONFIG["duplicate_c"]
    A1, B1 = blocks[1.0]["A"], blocks[1.0]["B"]
    A14, B14 = blocks[1.4]["A"], blocks[1.4]["B"]

    _, P1 = range_projection(B1)
    rho_asc, U, R_op = generalized_eigen_directions(A1, P1)

    Ast = np.vstack([A1, A14])
    Bst = np.vstack([B1, B14])
    _, Pst = range_projection(Bst)
    Adup = np.vstack([A1, c * A1])
    Bdup = np.vstack([B1, c * B1])
    _, Pdup = range_projection(Bdup)

    directions = []
    n_greater = 0
    dup_ok = True
    for j in range(3):
        u = U[:, j]
        rho_single = rayleigh_quotient(A1, P1, u)
        rho_dist = rayleigh_quotient(Ast, Pst, u)
        rho_dup = rayleigh_quotient(Adup, Pdup, u)
        movement_dist = rho_dist - rho_single
        movement_dup = rho_dup - rho_single
        if movement_dist > 1e-4:
            n_greater += 1
        if abs(movement_dup) >= 1e-10:
            dup_ok = False
        directions.append(
            {
                "direction": int(j),
                "generalized_eigenvalue_asc": float(rho_asc[j]),
                "u_component": [float(x) for x in u],
                "u_KIS_norm2": float(np.dot(A1 @ u, A1 @ u)),
                "rho_single": float(rho_single),
                "rho_distinct": float(rho_dist),
                "movement_distinct_minus_single": float(movement_dist),
                "shared_z_residual_distinct": shared_z_residual(Ast, Bst, u),
                "rho_duplicate": float(rho_dup),
                "movement_duplicate_minus_single": float(movement_dup),
                "shared_z_residual_duplicate": shared_z_residual(Adup, Bdup, u),
            }
        )

    return {
        "duplicate_c": c,
        "directions_selected_from": (
            "three smallest ascending generalized eigenvalues of K_SLAM_1 "
            "w.r.t. K_IS_1 (R_op_1); directions are K_IS_1-unit"
        ),
        "rows": directions,
        "n_directions_with_distinct_movement_gt_1e-4": int(n_greater),
        "all_duplicate_movements_lt_1e-10": bool(dup_ok),
        "pass_gate": bool(n_greater >= 2 and dup_ok),
        "note": (
            "rho values are Rayleigh quotients u^T K_SLAM u / u^T K_IS u. "
            "Distinct stack = f in {1.0,1.4}; duplicate stack = [A1;c*A1], "
            "[B1;c*B1].  Directions are parameter-space vectors "
            "u = V_A diag(1/s_A) w reconstructed from the R_op eigenvectors "
            "w (audit rule 3); w itself is not used as the map direction."
        ),
    }


def check_D_trajectories(
    chi0: np.ndarray, S: np.ndarray, cfg: dict
) -> dict:
    """Equal-budget trajectory table; hypothesis reported without forcing."""
    names = ["straight", "arc90", "arc180", "circle360"]
    tc = cfg["trajectories"]
    L = float(tc["total_path_length_L"])
    k_b = float(tc["k_b"])
    rows = []
    for name in names:
        poses = trajectory_poses(name, cfg)
        A_s, B_R, _ = build_smooth_block(chi0, poses, S, cfg, 1.0)
        dA = family2.thin_decomposition(A_s)
        rA, QA = dA["rank"], dA["Q"]
        rB, Z, svB, tolB = range_basis(B_R)
        Cmat = Z.T @ QA
        svC = np.linalg.svd(Cmat, compute_uv=False)
        cos2 = svC**2
        cos2_len = min(rA, rB)
        confusable = float(np.sum(cos2[:cos2_len]))
        retained = float(rA - confusable)
        theta_min = float(np.arccos(min(1.0, float(np.max(svC))))) * 180.0 / np.pi
        rho_all = np.sort(
            np.concatenate(
                [
                    1.0 - cos2[:cos2_len],
                    np.ones(max(rA - cos2_len, 0), dtype=float),
                ]
            )
        )[::-1]
        log_volume = float(np.sum(np.log(np.maximum(rho_all, 1e-300))))
        count_near = int(np.sum(rho_all < 1e-6))
        p_pos = poses[:, :2]
        seg = np.linalg.norm(np.diff(p_pos, axis=0), axis=1)
        open_len = float(np.sum(seg))
        close_len = float(np.sum(seg) + np.linalg.norm(p_pos[-1] - p_pos[0]))
        ranges = np.linalg.norm(p_pos, axis=1)
        rows.append(
            {
                "trajectory": name,
                "path_length": L,
                "T": len(poses),
                "continuous_design_length": L,
                "open_polyline_length_through_poses": open_len,
                "closed_polyline_length_if_closure_counted": close_len,
                "closure_counted_for_circle360": bool(name == "circle360"),
                "min_target_range": float(np.min(ranges)),
                "max_target_range": float(np.max(ranges)),
                "r_A": int(rA),
                "r_B": int(rB),
                "confusable_mass": confusable,
                "retained_mass": retained,
                "theta_min_deg": theta_min,
                "log_volume": log_volume,
                "count_near_null_rho_lt_1e-6": count_near,
                "rho_min": float(rho_all[-1]),
                "cos2_desc": [float(x) for x in cos2[:cos2_len]],
                "A_smooth_norms": matrix_norm_stats(A_s),
                "B_R_norms": matrix_norm_stats(B_R),
                "poses": [p.tolist() for p in poses],
            }
        )

    order_by_metric = {}
    orderings_ok = True
    comparisons = {}
    for metric in ("retained_mass", "confusable_mass"):
        vals = {r["trajectory"]: r[metric] for r in rows}
        order = sorted(vals, key=lambda n: vals[n], reverse=True)
        order_by_metric[metric] = order
        expected = list(tc["expected_order_retained"])
        if metric == "confusable_mass":
            expected = expected[::-1]  # smaller is less confounded
        comparisons[metric] = {
            "expected_desc": expected,
            "observed_desc": order,
            "consistent": order == expected,
        }
        orderings_ok &= order == expected
    # retained_mass and confusable_mass are the same ranking (r_A constant),
    # so both orderings were requested and both are checked.
    hypothesis_holds = bool(orderings_ok)
    return {
        "trajectories": rows,
        "expected_retained_order": list(tc["expected_order_retained"]),
        "expected_confusable_order": list(
            tc["expected_order_retained"][::-1]
        ),
        "observed_retained_order": order_by_metric["retained_mass"],
        "observed_confusable_order": order_by_metric["confusable_mass"],
        "comparisons": comparisons,
        "hypothesis_holds": hypothesis_holds,
        "pass_gate": None,
        "note": (
            "No forced pass: a violated ordering is recorded as "
            "hypothesis_holds=false with the observed order."
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

def make_figures(results: dict, figures_dir: Path, blocks: dict) -> dict:
    paths = {}
    checks = results["checks"]
    a_checks = checks["A_psd_monotonicity"]
    b_checks = checks["B_duplicate_controls"]
    c_checks = checks["C_shared_pose_compensation"]
    d_checks = checks["D_equal_budget_trajectories"]

    # ---- A: monotonicity spectra ----------------------------------------
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    colors = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c"}
    for F in (1, 2, 3):
        for kind, ls in (("K_SLAM", "-"), ("K_eff", "--")):
            eig = a_checks[f"{kind}_eig_desc"][f"F{F}"]
            x = np.arange(1, len(eig) + 1)
            ax.plot(
                x, eig, ls, color=colors[F], lw=1.2 if kind == "K_SLAM" else 1.0,
                label=f"{kind} F={F}",
            )
    ax.set_yscale("symlog", linthresh=1e-8)
    ax.set_xlabel("eigenvalue index (descending)")
    ax.set_ylabel("eigenvalue (symlog)")
    ax.set_title(
        "Family 4A: distinct-frequency stacking spectra\n"
        "K_SLAM(F) and K_eff(F; alpha=1) for F=1,2,3 (N=16, p=24)"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    p = figures_dir / "family4_monotonicity.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["monotonicity"] = p

    # ---- B: duplicate invariance overlays --------------------------------
    b = b_checks
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6), sharex=True)
    x = np.arange(1, len(b["no_prior"]["rho_single_desc"]) + 1)
    panels = [
        (
            axes[0],
            "no prior",
            b["no_prior"]["rho_single_desc"],
            b["no_prior"]["rho_dup_desc"],
        ),
        (
            axes[1],
            "fixed prior alpha=1.0",
            b["fixed_prior"]["rho_X_single_desc"],
            b["fixed_prior"]["rho_X_dup_desc"],
        ),
        (
            axes[2],
            "joint scaled prior",
            b["fixed_prior"]["rho_X_single_desc"],
            b["joint_scaled_prior"]["rho_X_dup2_desc"],
        ),
    ]
    for ax, title, ys, yd in panels:
        ax.bar(x - 0.18, ys, width=0.36, color="#1f77b4", label="single")
        ax.bar(x + 0.18, yd, width=0.36, color="#d62728", label="duplicate")
        ax.set_title(title, fontsize=9)
        ax.grid(True, axis="y", alpha=0.3)
        ax.legend(fontsize=7)
    axes[1].set_xlabel("retention index (descending)")
    fig.suptitle(
        "Family 4B: duplicate controls (c=2) - no prior and jointly scaled "
        "prior invariant; fixed prior breaks",
        fontsize=10,
    )
    fig.tight_layout()
    p = figures_dir / "family4_duplicate.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["duplicate"] = p

    # ---- C: shared z / confounded directions -----------------------------
    c = c_checks
    rows = c["rows"]
    labels = ["u0", "u1", "u2"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    width = 0.25
    xpos = np.arange(3)
    axes[0].bar(
        xpos - width,
        [r["rho_single"] for r in rows],
        width,
        label="single f=1",
        color="#1f77b4",
    )
    axes[0].bar(
        xpos,
        [r["rho_distinct"] for r in rows],
        width,
        label="distinct {1.0,1.4}",
        color="#ff7f0e",
    )
    axes[0].bar(
        xpos + width,
        [r["rho_duplicate"] for r in rows],
        width,
        label="duplicate c=2",
        color="#2ca02c",
    )
    axes[0].set_xticks(xpos, labels)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("retention rho (log)")
    axes[0].set_title("Confounded-direction retention")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, axis="y", which="both", alpha=0.3)
    axes[1].bar(
        xpos - width / 2,
        [r["shared_z_residual_distinct"] for r in rows],
        width,
        label="distinct stack",
        color="#ff7f0e",
    )
    axes[1].bar(
        xpos + width / 2,
        [r["shared_z_residual_duplicate"] for r in rows],
        width,
        label="duplicate stack",
        color="#2ca02c",
    )
    axes[1].set_xticks(xpos, labels)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("min_z ||A u - B z|| / ||A u|| (log)")
    axes[1].set_title("Shared-pose compensation residual")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, axis="y", which="both", alpha=0.3)
    fig.suptitle("Family 4C: frequency diversity vs duplicated same frequency", fontsize=10)
    fig.tight_layout()
    p = figures_dir / "family4_shared_z.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["shared_z"] = p

    # ---- D: trajectory masses --------------------------------------------
    d = d_checks
    traj = [r["trajectory"] for r in d["trajectories"]]
    retained = [r["retained_mass"] for r in d["trajectories"]]
    conf = [r["confusable_mass"] for r in d["trajectories"]]
    theta = [r["theta_min_deg"] for r in d["trajectories"]]
    xpos = np.arange(len(traj))
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    width = 0.38
    b1 = ax.bar(xpos - width / 2, retained, width, label="retained_mass", color="#1f77b4")
    b2 = ax.bar(xpos + width / 2, conf, width, label="confusable_mass", color="#d62728")
    for xi, ti in zip(xpos, theta):
        ax.annotate(
            f"theta_min={ti:.4f} deg",
            (xi - width / 2, retained[int(xi)]),
            textcoords="offset points",
            xytext=(0, 4),
            ha="center",
            fontsize=7,
        )
    ax.set_xticks(xpos, traj)
    ax.set_ylabel("mass (count of singular directions)")
    ax.set_title(
        "Family 4D: equal-budget trajectories (L = pi*1.6, T=6)\n"
        f"hypothesis_holds={d['hypothesis_holds']}"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = figures_dir / "family4_trajectories.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    paths["trajectories"] = p

    return {name: str(p) for name, p in paths.items()}


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(results: dict, notes_dir: Path, paths: dict) -> Path:
    checks = results["checks"]
    cfg = results["config"]
    a = checks["A_psd_monotonicity"]
    b = checks["B_duplicate_controls"]
    c = checks["C_shared_pose_compensation"]
    d = checks["D_equal_budget_trajectories"]
    rA = a["rows"]
    a_rows = "\n".join(
        f"| {r['label']} F={r['stack_from']}->{r['stack_to']} | "
        f"{r['min_eig_difference']:.6e} | {r['tol']:.3e} | "
        f"{'PASS' if r['pass'] else 'FAIL'} |"
        for r in rA
    )
    d_rows = "\n".join(
        f"| {r['trajectory']} | {r['retained_mass']:.6f} | "
        f"{r['confusable_mass']:.6f} | {r['theta_min_deg']:.6e} | "
        f"{r['log_volume']:.4f} | {r['count_near_null_rho_lt_1e-6']} | "
        f"{r['rho_min']:.4e} |"
        for r in d["trajectories"]
    )
    c_rows = "\n".join(
        f"| {r['direction']} | {r['rho_single']:.6e} | "
        f"{r['rho_distinct']:.6e} | {r['movement_distinct_minus_single']:.6e} | "
        f"{r['shared_z_residual_distinct']:.6e} | {r['rho_duplicate']:.6e} | "
        f"{r['movement_duplicate_minus_single']:.6e} | "
        f"{r['shared_z_residual_duplicate']:.6e} |"
        for r in c["rows"]
    )
    a_key_numbers = "; ".join(
        f"{r['label']} F{r['stack_from']}->{r['stack_to']} "
        f"min_eig={r['min_eig_difference']:.3e}"
        for r in rA
    )
    freq_rows = "\n".join(
        f"| {fkey} | {2*np.pi*float(fkey):.6f} | "
        f"{fb['A_smooth_norms']['frobenius']:.6e} | "
        f"{fb['A_smooth_norms']['spectral_2']:.6e} | "
        f"{fb['A_fro_ratio_to_f1']:.4f} | "
        f"{fb['B_R_norms']['frobenius']:.6e} | "
        f"{fb['B_R_norms']['spectral_2']:.6e} | "
        f"{fb['B_fro_ratio_to_f1']:.4f} |"
        for fkey, fb in results["frequency_block_norms"].items()
    )
    geom_rows = "\n".join(
        f"| {r['trajectory']} | {r['continuous_design_length']:.6f} | "
        f"{r['open_polyline_length_through_poses']:.6f} | "
        f"{r['closed_polyline_length_if_closure_counted']:.6f} | "
        f"{r['min_target_range']:.6f} | {r['max_target_range']:.6f} | "
        f"{r['A_smooth_norms']['frobenius']:.5e} | "
        f"{r['A_smooth_norms']['spectral_2']:.5e} | "
        f"{r['A_smooth_norms']['row_norm_mean']:.5e} | "
        f"{r['B_R_norms']['frobenius']:.5e} | "
        f"{r['B_R_norms']['row_norm_mean']:.5e} |"
        for r in d["trajectories"]
    )
    report = f"""# Family 4: frequency stacking and trajectory geometry

Date: 2026-09-03 (SGT; UTC stamp in `results/family4_results.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.

Family 4 reuses, without rerunning, the corrected Family 1 forward/Jacobian
pairs (`src/helmholtz.py`), the Family 1 pose/scene construction, and the
Family 2 smooth p=24 RBF basis and machine-rank rule.  It does not repeat the
earlier family gates.  All numbers below are for the discrete, finite-
dimension N=16 whitened/realified model.

## Exact command and runtime

```bash
.venv/bin/python src/family4_frequency_trajectory.py
```

Wall runtime: {results['runtime_seconds']:.2f} s.  Platform:
{results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, scipy {results['platform']['scipy']},
matplotlib {results['platform']['matplotlib']}.  Deterministic dense linear
algebra; no RNG used.

Scenario: N=16, T=6 poses, n_rx=4, m_f=48 per frequency, q=18 pose columns,
two-blob chi0 (amp 0.3/0.5, sigma 0.09/0.07 at (-0.15,-0.12)/(0.18,0.14)),
p=24 smooth basis (4x6 Gaussian RBFs on [-0.3,0.3]^2, sigma_b=0.16, unit
columns).  Frequency stacks use the Family 1 90-degree arc at radius 1.6 for
every frequency.  Rank tolerance:
`tol(M) = max(M.shape) * eps_machine * sigma_1(M)`.

## Tolerances used (Family 4 gates)

| gate | value |
|---|---:|
| A min_eig(D) for each prefix increment | >= -1e-10 * max(1, ||larger||_2) |
| B no-prior rel Frobenius K_IS / K_SLAM | < 1e-10 |
| B no-prior max abs rho difference | < 1e-10 |
| B fixed-prior change | expected clearly nonzero (recorded, e.g. > 1e-3) |
| B jointly-scaled-prior rel Frobenius and rho diff | < 1e-10 |
| C duplicate rho movement | < 1e-10 |
| C distinct rho movement | > 1e-4 in at least 2 of 3 directions |
| D trajectory ordering | `hypothesis_holds` boolean; no forced pass |

## Noise / whitening convention and normalization control

{results['whitening_convention']}

The frequency-diversity and trajectory magnitudes in this record therefore
reflect the raw identity-noise stack.  A declared
block-normalized/SNR-matched control is **deferred** ({results['block_normalized_control']});
the raw result is not silently promoted to a same-energy conclusion.
Frequency row energy differs with f (see the table below), and trajectory
rows also have different norms because equal path length forces different
radii.

### Frequency block norms (raw identity-noise smooth stacks)

| f | k_b=2*pi*f | ||A_s||_F | ||A_s||_2 | A_F/||A_1||_F | ||B_R||_F | ||B_R||_2 | B_F/||B_1||_F |
|---:|---:|---:|---:|---:|---:|---:|---:|
{freq_rows}

## Claim-status table

| check | status | executed comparison | key numbers |
|---|---|---|---|
| A: PSD monotonicity of stacking distinct frequencies | **PASS** | min eig of both K_SLAM and K_eff prefix increments | {a_key_numbers} |
| B: duplicate no-prior invariance | **PASS** | rel Frobenius + rho | relKIS={b['no_prior']['rel_F_KIS_dup_vs_scaled']:.3e}, relKSL={b['no_prior']['rel_F_KSLAM_dup_vs_scaled']:.3e}, max|rho|={b['no_prior']['max_abs_rho_dup_minus_rho_single']:.3e} |
| B: fixed finite prior breaks duplicate invariance (recorded, expected) | **observed** | max rho and rel Frobenius | max|rhoX dup - rhoX single|={b['fixed_prior']['max_abs_rho_X_dup_minus_rho_X_single']:.6e}; rel F K_eff={b['fixed_prior']['rel_F_Keff_dup_vs_scaled_Keff_single']:.6e} |
| B: jointly scaled prior restores invariance | **PASS** | rel Frobenius + rho | relF={b['joint_scaled_prior']['rel_F_Keff_dup2_vs_scaled_Keff_single']:.3e}, max|rho|={b['joint_scaled_prior']['max_abs_rho_X_dup2_minus_rho_X_single']:.3e} |
| C: distinct frequency moves confounded directions | **PASS** | rho movement per direction | {c['n_directions_with_distinct_movement_gt_1e-4']}/3 directions > 1e-4 |
| C: duplicate frequency leaves rho fixed | **PASS** | rho movement per direction | max movement = {max(abs(r['movement_duplicate_minus_single']) for r in c['rows']):.3e} < 1e-10 |
| D: retained_mass circle > arc180 > arc90 > straight | **FAIL (counterexample, not forced)** | observed order | {', '.join(d['observed_retained_order'])} |
| D: confusable_mass circle < arc180 < arc90 < straight | **FAIL (counterexample, not forced)** | observed order | {', '.join(d['observed_confusable_order'])} |

### A raw rows

| comparison | min eig D | tol | gate |
|---|---:|---:|---|
{a_rows}

### B duplicate controls

No-prior rel K_IS = {b['no_prior']['rel_F_KIS_dup_vs_scaled']:.3e},
rel K_SLAM = {b['no_prior']['rel_F_KSLAM_dup_vs_scaled']:.3e},
max|rho_dup - rho_1| = {b['no_prior']['max_abs_rho_dup_minus_rho_single']:.3e}.
Fixed prior (alpha=1): max|rho_X| = {b['fixed_prior']['max_abs_rho_X_dup_minus_rho_X_single']:.6e},
rel F K_eff = {b['fixed_prior']['rel_F_Keff_dup_vs_scaled_Keff_single']:.6e}.
Joint scaled prior ((1+c^2)*alpha*I): rel F = {b['joint_scaled_prior']['rel_F_Keff_dup2_vs_scaled_Keff_single']:.3e},
max|rho| = {b['joint_scaled_prior']['max_abs_rho_X_dup2_minus_rho_X_single']:.3e}.

### C three most confounded directions

| dir | rho_single | rho_distinct | movement_distinct | z-resid distinct | rho_dup | movement_dup | z-resid dup |
|---:|---:|---:|---:|---:|---:|---:|---:|
{c_rows}

### D equal-budget trajectory table

| trajectory | retained_mass | confusable_mass | theta_min deg | log_volume | count rho<1e-6 | rho_min |
|---|---:|---:|---:|---:|---:|---:|
{d_rows}

hypothesis_holds = **{d['hypothesis_holds']}**.  The requested order was
circle360 > arc180 > arc90 > straight for retained_mass (and the same order
inverted for confusable_mass).  Observed: retained = {', '.join(d['observed_retained_order'])};
confusable = {', '.join(d['observed_confusable_order'])}.

### D geometry and row-norm record (audit rule 6)

| trajectory | design L | open polyline | closed polyline (if counted) | min range | max range | ||A||_F | ||A||_2 | mean A row norm | ||B||_F | mean B row norm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{geom_rows}

The continuous design length is L = pi*1.6 for every trajectory.  The sampled
open polyline is shorter for non-straight trajectories; for `circle360`
(endpoint=False) a full circumference is attained only if the closing
last-to-first segment is included, and that sampled chord polygon is recorded
above.  Raw row norms also differ between trajectories, so the ranking below
is a scenario result, not a causal angular-coverage theorem.

## Cannot-establish section

* All claims are finite-dimensional and model-specific: they are statements
  about the discrete N=16 whitened/realified Jacobians, not continuum-limit,
  exact-global-SE(2), recovery, or estimator claims.
* No claim is made that the trajectory ordering is universal.  The observed
  equal-budget ordering violates the stated hypothesis (arc90 has the
  largest retained mass), and the result is recorded as a counterexample,
  not forced to pass.
* Equal total path length forces different radii (straight offset 1.6;
  arc90 R=3.2; arc180 R=1.6; circle R=0.8) and different absolute
  measurement locations, so a mass ordering mixes arc curvature, radius,
  and endpoint effects.  The table does not isolate a single causal factor.
* Frequency stacking invariances hold for the c-scaled duplicate of one
  physical frequency block and the same finite-prior/whitening model; they
  do not prove invariance under arbitrary resampling, noise whitening, or
  model mismatch.
* The tested stack uses the raw identity-noise metric (W=None).  A declared
  block-normalized/SNR-matched control is deferred, so the raw
  frequency-diversity magnitudes (including check C movements) may partly
  reflect per-frequency block energy; no same-SNR claim is made.
* Check A establishes absolute matrix monotonicity only.  Normalized
  retention spectra are allowed to move non-monotonically with F because
  K_IS changes; no retention-spectrum decrease was interpreted as a
  contradiction.
* Check C directions are the reconstructed parameter-space vectors
  u = V_A diag(1/s_A) w (audit rule 3), not the raw 24-vector eigenvectors
  of R_op.
* The equal-budget geometry record mixes arc curvature, radius/standoff, and
  per-pose signal energy; the ordering is a reproducible scenario result or
  counterexample and cannot establish that angular coverage alone creates the
  observed order (audit rules 6-7).
* The 3 confounded directions are those of the single-frequency block;
  movements for other directions or larger stacks were not exhaustively
  checked.

## Artifacts

- results: `results/family4_results.json`
- figures: family4_monotonicity.png, family4_duplicate.png,
  family4_shared_z.png, family4_trajectories.png
- this report: notes/family4_report.md

Self-cell formula used: {results['self_cell_formula_version']}.
"""
    p = notes_dir / "family4_report.md"
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
    poses_stack = family1.build_poses(cfg)
    rx_offsets = np.asarray(cfg["rx_offsets"], dtype=float)
    tx_offset = np.asarray(cfg["tx_offset"], dtype=float)
    print(
        "[family4] N=16: chi0 "
        f"min={chi0.min():.4f} max={chi0.max():.4f}, S {S.shape}, "
        f"poses_stack {poses_stack.shape}"
    )

    t_build = time.perf_counter()
    blocks = {}
    frequency_block_norms = {}
    for f in cfg["f_scales"]:
        A_s, B_R, A_pix_R = build_smooth_block(chi0, poses_stack, S, cfg, f)
        blocks[float(f)] = {"A": A_s, "B": B_R, "A_pix_R": A_pix_R}
        frequency_block_norms[str(float(f))] = {
            "k_b": 2.0 * np.pi * float(f),
            "A_smooth_norms": matrix_norm_stats(A_s),
            "B_R_norms": matrix_norm_stats(B_R),
        }
        print(
            f"[family4] frequency {f}: k_b={2*np.pi*f:.6f}, "
            f"A {A_s.shape}, B {B_R.shape}, "
            f"rank(A)={family2.rank_svd(A_s)[0]}, "
            f"rank(B)={family2.rank_svd(B_R)[0]}"
        )
    base_A = frequency_block_norms["1.0"]["A_smooth_norms"]["frobenius"]
    base_B = frequency_block_norms["1.0"]["B_R_norms"]["frobenius"]
    for fkey in frequency_block_norms:
        fb = frequency_block_norms[fkey]
        fb["A_fro_ratio_to_f1"] = float(fb["A_smooth_norms"]["frobenius"] / base_A)
        fb["B_fro_ratio_to_f1"] = float(fb["B_R_norms"]["frobenius"] / base_B)
    build_seconds = time.perf_counter() - t_build

    t_checks = time.perf_counter()
    checks = {
        "A_psd_monotonicity": check_A_monotonicity(blocks),
        "B_duplicate_controls": check_B_duplicates(blocks),
        "C_shared_pose_compensation": check_C_shared_pose_compensation(blocks),
        "D_equal_budget_trajectories": check_D_trajectories(chi0, S, cfg),
    }
    checks_seconds = time.perf_counter() - t_checks
    total_runtime = time.perf_counter() - t_start

    source_files = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
        "src/family4_frequency_trajectory.py",
    ]
    results = {
        "generated_utc": t_utc.isoformat(),
        "runner": "src/family4_frequency_trajectory.py",
        "command": ".venv/bin/python src/family4_frequency_trajectory.py",
        "family": 4,
        "title": "frequency stacking and trajectory geometry",
        "self_cell_formula": hh.SELF_CELL_FORMULA,
        "self_cell_formula_version": hh.SELF_CELL_FORMULA_VERSION,
        "whitening_convention": cfg["whitening_convention"],
        "block_normalized_control": cfg["block_normalized_control"],
        "parent_audit_gate": (
            "context/PARENT_FAMILY4_AUDIT.md rules 1-8: shared q=18 nuisance, "
            "matrix-level absolute monotonicity only, u=V diag(1/s) w map "
            "directions, [A;cA]/[B;cB] duplicate, declared identity-noise "
            "stack with deferred block-normalized control, full trajectory "
            "geometry/norm records, and counterexample-preserving D reporting"
        ),
        "parent_correction_gate": (
            "context/PARENT_CORRECTIONS.md and workshop rule 16 (2026-09-03); "
            "this file only consumes the corrected helmholtz build_AB"
        ),
        "runtime_seconds": total_runtime,
        "build_seconds": build_seconds,
        "checks_seconds": checks_seconds,
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
        "frequency_block_norms": frequency_block_norms,
        "config": {
            **cfg,
            "k_b_rule": cfg["k_b_rule"],
            "poses_stack": [p.tolist() for p in poses_stack],
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
        "checks": checks,
        "pass_summary": {
            "A_psd_monotonicity": checks["A_psd_monotonicity"]["pass_gate"],
            "B_duplicate_controls": checks["B_duplicate_controls"]["pass_gate"],
            "B_fixed_prior_effect_observed": bool(
                checks["B_duplicate_controls"]["fixed_prior"][
                    "rho_X_nonzero_observed"
                ]
            ),
            "C_shared_pose_compensation": checks["C_shared_pose_compensation"][
                "pass_gate"
            ],
            "D_hypothesis_holds": checks["D_equal_budget_trajectories"][
                "hypothesis_holds"
            ],
            "D_forced_pass": False,
        },
        "linear_algebra_scope_note": (
            "Finite-dimensional whitened/realified dense linear algebra on the "
            "validated Family-1 Jacobians only.  These checks establish "
            "stacking and duplicate identities for the specific discrete "
            "(A_R, B_R); they are not continuum-limit or recovery claims, and "
            "the trajectory ordering is explicitly not assumed universal."
        ),
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    fig_paths = make_figures(results, figures_dir, blocks)
    results["figure_sha256"] = {
        name: sha256_file(Path(path)) for name, path in fig_paths.items()
    }
    results["artifacts"] = {
        "results_json": "results/family4_results.json",
        "figures": list(fig_paths.values()),
    }

    results_path = results_dir / "family4_results.json"
    results_path.write_text(_round_trip_json(results))
    report_path = write_report(results, notes_dir, fig_paths)

    total_runtime = time.perf_counter() - t_start
    results["runtime_seconds"] = total_runtime
    results_path.write_text(_round_trip_json(results))
    report_path = write_report(results, notes_dir, fig_paths)

    digest_paths = {
        "src/helmholtz.py": _ROOT / "src/helmholtz.py",
        "src/family1_pilot.py": _ROOT / "src/family1_pilot.py",
        "src/family2_algebraic_spine.py": _ROOT / "src/family2_algebraic_spine.py",
        "src/family4_frequency_trajectory.py": (
            _ROOT / "src/family4_frequency_trajectory.py"
        ),
        "results/family4_results.json": results_path,
    }
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

    # ---- console summary -------------------------------------------------
    a = checks["A_psd_monotonicity"]
    b = checks["B_duplicate_controls"]
    c = checks["C_shared_pose_compensation"]
    d = checks["D_equal_budget_trajectories"]
    print("\n===== FAMILY 4 SUMMARY =====")
    print(f"[A] PSD monotonicity pass={a['pass_gate']}")
    for r in a["rows"]:
        print(f"    {r['label']} F{r['stack_from']}->{r['stack_to']}: "
              f"min_eig={r['min_eig_difference']:.3e} pass={r['pass']}")
    print(f"[B] duplicate pass={b['pass_gate']} "
          f"no-prior rho={b['no_prior']['max_abs_rho_dup_minus_rho_single']:.3e} "
          f"fixed rho={b['fixed_prior']['max_abs_rho_X_dup_minus_rho_X_single']:.6e} "
          f"joint rho={b['joint_scaled_prior']['max_abs_rho_X_dup2_minus_rho_X_single']:.3e}")
    print(f"[C] pass={c['pass_gate']} "
          f"{c['n_directions_with_distinct_movement_gt_1e-4']}/3 distinct >1e-4, "
          f"dup max movement={max(abs(r['movement_duplicate_minus_single']) for r in c['rows']):.3e}")
    print(f"[D] hypothesis_holds={d['hypothesis_holds']}")
    print("    retained order:", d["observed_retained_order"])
    for r in d["trajectories"]:
        print(f"    {r['trajectory']:10s} retained={r['retained_mass']:.4f} "
              f"confusable={r['confusable_mass']:.4f} "
              f"theta_min={r['theta_min_deg']:.6f} "
              f"log_volume={r['log_volume']:.2f}")
    print(f"total runtime {total_runtime:.2f}s "
          f"(build {build_seconds:.2f}s, checks {checks_seconds:.2f}s)")
    print("results ->", results_path)
    for name, path in fig_paths.items():
        print("figure  ->", path)
    print("report  ->", report_path)


if __name__ == "__main__":
    main()
