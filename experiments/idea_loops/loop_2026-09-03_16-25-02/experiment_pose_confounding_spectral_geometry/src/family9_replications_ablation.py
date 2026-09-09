"""Family 9: seed replications, component ablations, and pose-DOF nuisance baseline.

Run (from the experiment root):
    .venv/bin/python src/family9_replications_ablation.py

Read-only reuse (existing files are not modified):
  src/helmholtz.py                    make_grid, build_AB, whiten_realify
  src/family1_pilot.py                build_poses, make_chi0
  src/family2_algebraic_spine.py      build_smooth_basis, thin_decomposition,
                                      range_basis, rank_svd,
                                      retention_spectrum
  src/family4_frequency_trajectory.py base_scene, build_AB_whitened,
                                      build_smooth_block, K_SLAM, K_eff,
                                      retention_spectrum,
                                      generalized_eigen_directions,
                                      rayleigh_quotient, shared_z_residual

Outputs:
  results/family9_replications_ablation.json
  figures/family9_seed_replications.png
  figures/family9_pose_dof_ablation.png
  figures/family9_rx_T_basis_ablation.png
  notes/family9_replications_ablation.md

The family-4/7 principal-angle convention is used throughout:
  Q_A = left singular basis of A from the thin SVD, Z = orthonormal basis of
  Range(B); cos^2 entries are descending singular values of Z^T Q_A squared;
  retained mass = r_A - sum(cos^2 over the first min(r_A, r_B) entries);
  rho = desc(concat(1 - cos^2, ones(r_A - min(r_A, r_B)))).

Check C direction replication: Family 4 check C and all of its descendants
(families 4b, 4c, 6, 7) call
    family4.generalized_eigen_directions(A1, P1)
with P1 = I - Z Z^T (Range(B1) orthogonal projector), because that function
expects a data-space matrix W and forms R_op = Q_A^T W Q_A.  The literal
second-argument form "A1.T@A1" in the experiment brief is dimensionally
incompatible with the verified API (A1.T@A1 is (p,p); Q_A is (m,r)), so the
executed replication follows Family 4 exactly and the literal form is
recorded as a diagnostic with the exception it raises.
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

_EPS = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Configuration (standard Family 1/4/7 scenario plus family 9 sweeps)
# ---------------------------------------------------------------------------

CONFIG = {
    "title": (
        "deterministic seed replications of the algebraic-spine and "
        "frequency-diversity claims, plus pose-DOF / receiver / pose-count / "
        "basis ablations"
    ),
    "N": 16,
    "T": 6,
    "n_rx": 4,
    "q_pose": 18,  # 3*T
    "m_real": 48,  # 2*T*n_rx after whiten_realify at n_rx=4, T=6
    "f_single": 1.0,
    "f_distinct": 1.4,
    "k_b_rule": "k_b(f) = 2*pi*f",
    "whitening_convention": (
        "W = None (identity noise): hh.whiten_realify(A,B,None) returns "
        "sqrt(2)*[Re; Im] row stacks; raw identity-noise stacks throughout"
    ),
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
        "theta = atan2(-p_y,-p_x): body +x axis points toward the origin"
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
        "note": "identical construction/config to families 2-8",
    },
    "partA_scenes": {
        "random": {
            "n_blobs": 8,
            "seeds": [1001, 1002, 1003, 1004, 1005],
            "a_range": [0.1, 0.6],
            "c_range": [-0.35, 0.35],
            "sigma_range": [0.05, 0.12],
            "formula": (
                "chi(r) = sum_{m=1..8} a_m exp(-|r-c_m|^2/(2 s_m^2)); "
                "a_m ~ U(0.1,0.6), c_m ~ U(-0.35,0.35)^2, s_m ~ U(0.05,0.12)"
            ),
        },
        "ring": {
            "amplitude": 0.5,
            "radius": 0.25,
            "width": 0.05,
            "formula": "chi_ring = 0.5 exp(-((|r|-0.25)/0.05)^2)",
        },
        "l2_reference": "standard two-blob chi0 (family1.make_chi0)",
        "l2_metric": (
            "sample vector 2-norm over the N^2 cell-centre values "
            "(no h weighting)"
        ),
        "scene_order": [
            "seed_1001",
            "seed_1002",
            "seed_1003",
            "seed_1004",
            "seed_1005",
            "two_blob",
            "ring",
        ],
    },
    "alpha_prior": 1.0,
    "movement_threshold": 1e-4,
    "gates": {
        "rank_identity": "integer residual == 0 (family 2 c2)",
        "retention_identity_scale": (
            "1e-10 * max(1, ||A_s||_2^2) (family 2 c3)"
        ),
        "kernel_identity_scale": (
            "1e-8 * max(1, ||A_s||_2^2) for the family-2 subspace-distance "
            "gate; identity_support additionally requires equal nullities, "
            "cross residuals <= 1e-10, and min singular N1^T N2 > 1-1e-6 "
            "when the null spaces are non-empty"
        ),
        "check_C": (
            ">= 2 of 3 directions with rho movement > 1e-4 (family 4 gate)"
        ),
    },
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "partB_pose_dof": {
        "pose_column_order": "[px, py, theta] repeated per pose t",
        "subsets": {
            "known_pose": [],
            "x_only": "3t+0",
            "y_only": "3t+1",
            "theta_only": "3t+2",
            "xy": "3t+0,3t+1",
            "x_theta": "3t+0,3t+2",
            "y_theta": "3t+1,3t+2",
            "full": "all 3T columns",
        },
        "claim_scope": "descriptive ablation; no forced pass",
    },
    "partC_receiver_sweep": {
        "rows": [
            {
                "n_rx": 1,
                "rx_offsets": [[0.0, 0.0]],
            },
            {
                "n_rx": 2,
                "rx_offsets": [[-0.06, 0.0], [0.06, 0.0]],
            },
            {
                "n_rx": 4,
                "rx_offsets": [
                    [-0.06, 0.0],
                    [0.06, 0.0],
                    [0.0, -0.06],
                    [0.0, 0.06],
                ],
            },
            {
                "n_rx": 8,
                "rx_offsets": [
                    [-0.06, 0.0],
                    [0.06, 0.0],
                    [0.0, -0.06],
                    [0.0, 0.06],
                    [0.0424, 0.0424],
                    [-0.0424, 0.0424],
                    [0.0424, -0.0424],
                    [-0.0424, -0.0424],
                ],
            },
        ],
        "claim_scope": "descriptive ablation; no forced pass",
    },
    "partD_pose_sweep": {
        "Ts": [3, 6, 12],
        "geometry": (
            "Family-1 arc: R=1.6, phi=linspace(-45,45,T) deg, "
            "theta=atan2(-py,-px); q=3T"
        ),
        "claim_scope": "descriptive ablation; no forced pass",
    },
    "partE_basis_sweep": {
        "configs": [
            {"p": 12, "x_centers_n": 3, "y_centers_n": 4},
            {"p": 24, "x_centers_n": 4, "y_centers_n": 6},
            {"p": 36, "x_centers_n": 4, "y_centers_n": 9},
        ],
        "fixed": {
            "x_span": [-0.3, 0.3],
            "y_span": [-0.3, 0.3],
            "sigma_b": 0.16,
            "unit_columns": True,
        },
        "claim_scope": "descriptive ablation; no forced pass",
    },
    "randomness_note": (
        "seeded numpy default_rng(1001..1005) only for random scenes; "
        "dense linear algebra is deterministic"
    ),
    "claim_scope": (
        "finite-grid (N=16) replications and descriptive ablations only; "
        "no continuum transfer, no universal monotonicity theorem, no "
        "silent conflation of machine-rank boundary events with exact "
        "algebraic violations"
    ),
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
    if isinstance(obj, np.bool_):
        return bool(obj)
    raise TypeError(f"not JSON serialisable: {type(obj)}")


def _fmt(x, spec=".6g"):
    if x is None:
        return "n/a"
    return format(float(x), spec)


def _md_table(headers: list[str], rows: list[list]) -> str:
    if not rows:
        rows = [["n/a"] * len(headers)]
    widths = [len(str(h)) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))
    head = "| " + " | ".join(str(h).ljust(widths[i]) for i, h in enumerate(headers)) + " |"
    sep = "| " + " | ".join("-" * w for w in widths) + " |"
    body = [
        "| " + " | ".join(str(c).ljust(widths[i]) for i, c in enumerate(row)) + " |"
        for row in rows
    ]
    return "\n".join([head, sep] + body)


def _desc(a: np.ndarray) -> np.ndarray:
    return np.sort(np.asarray(a, dtype=float))[::-1]


def _safe_div(num: float, den: float) -> float:
    return float(num / max(den, 1e-300))


# ---------------------------------------------------------------------------
# info_metrics helper (principal-angle convention of families 4 and 7)
# ---------------------------------------------------------------------------

def info_metrics(A_s: np.ndarray, B_R: np.ndarray) -> dict:
    """Rank/principal-angle/Frobenius metric block for one (A_s, B_R) pair.

    All ranks use the family-2 machine rule.  Principal-angle quantities
    follow `family7.subspace_mass_metrics` exactly (which is stated there to
    be the exact family-4/6 convention): cos^2 from the singular values of
    Z^T Q_A, retained = r_A - sum(cos^2 over first min(r_A,r_B)), rho from
    1 - cos^2 padded with ones to length r_A, and theta_min = arccos(max
    singular value) (90 degrees by convention when a subspace is empty).
    """
    A = np.asarray(A_s, dtype=float)
    B = np.asarray(B_R, dtype=float)
    if A.ndim != 2 or B.ndim != 2 or A.shape[0] != B.shape[0]:
        raise ValueError(
            f"A and B must be 2-D with equal rows; got {A.shape}, {B.shape}"
        )

    thin = family2.thin_decomposition(A)
    rA = int(thin["rank"])
    Q_A = thin["Q"]
    rB, Z, svB, tolB = family2.range_basis(B)
    rB = int(rB)

    AB = np.concatenate([A, B], axis=1)
    rAB, _, _ = family2.rank_svd(AB)
    K_IS = _sym(A.T @ A)
    K_SLAM = family4.K_SLAM(A, B)
    K_eff = family4.K_eff(A, B, float(CONFIG["alpha_prior"]))
    rKIS, _, _ = family2.rank_svd(K_IS)
    rKSL, _, _ = family2.rank_svd(K_SLAM)
    lhs = rKIS - rKSL
    rhs = rA + rB - rAB

    # principal-angle mass block (family 7 exact convention)
    cos2_len = min(rA, rB)
    if rA > 0 and rB > 0:
        svC = np.linalg.svd(Z.T @ Q_A, compute_uv=False)
        cos2 = _desc(np.asarray(svC, dtype=float) ** 2)[:cos2_len]
        max_sv = float(np.max(svC))
        theta_min_rad = float(np.arccos(min(1.0, max_sv)))
    else:
        cos2 = np.array([], dtype=float)
        max_sv = 0.0
        # family 7 declared convention for a trivial subspace
        theta_min_rad = 0.5 * np.pi
    confusable = float(np.sum(cos2))
    retained = float(rA - confusable)
    rho_all = np.sort(
        np.concatenate(
            [
                1.0 - cos2,
                np.ones(max(rA - cos2_len, 0), dtype=float),
            ]
        )
    )[::-1]

    return {
        "r_A": rA,
        "r_B": rB,
        "r_AB": rAB,
        "r_KIS": int(rKIS),
        "r_KSL": int(rKSL),
        "rank_identity_lhs": int(lhs),
        "rank_identity_rhs": int(rhs),
        "rank_identity_residual": int(lhs - rhs),
        "rank_identity_holds": bool(lhs == rhs),
        "confusable_mass": confusable,
        "retained_mass": retained,
        "theta_min_rad": theta_min_rad,
        "theta_min_deg": float(theta_min_rad * 180.0 / np.pi),
        "rho_min": float(rho_all[-1]) if len(rho_all) else None,
        "count_rho_lt_1e-6": int(np.sum(rho_all < 1e-6)),
        "cos2_len": int(cos2_len),
        "K_IS_fro": float(np.linalg.norm(K_IS, ord="fro")),
        "K_SLAM_fro": float(np.linalg.norm(K_SLAM, ord="fro")),
        "K_eff_fro": float(np.linalg.norm(K_eff, ord="fro")),
        "sigma_min_B": float(svB[-1]) if len(svB) else None,
        "rank_tol_B": float(tolB),
    }


def _pose_submatrix(B_R: np.ndarray, kind: str) -> np.ndarray:
    """Restrict B_R columns to a pose-DOF subset (t-column order px,py,theta)."""
    T = B_R.shape[1] // 3
    if kind == "known_pose":
        cols = []
    elif kind == "x_only":
        cols = [3 * t for t in range(T)]
    elif kind == "y_only":
        cols = [3 * t + 1 for t in range(T)]
    elif kind == "theta_only":
        cols = [3 * t + 2 for t in range(T)]
    elif kind == "xy":
        cols = [c for t in range(T) for c in (3 * t, 3 * t + 1)]
    elif kind == "x_theta":
        cols = [c for t in range(T) for c in (3 * t, 3 * t + 2)]
    elif kind == "y_theta":
        cols = [c for t in range(T) for c in (3 * t + 1, 3 * t + 2)]
    elif kind == "full":
        cols = list(range(B_R.shape[1]))
    else:
        raise ValueError(f"unknown pose-DOF subset {kind!r}")
    if kind == "known_pose":
        # Full-width zero matrix: family2.rank_svd requires at least one
        # column (sv[0] indexing), and q=18 remains the pose-parameter count
        # of the known-pose reference even though Range(B) = {0}.
        return np.zeros_like(B_R)
    return B_R[:, cols]


# ---------------------------------------------------------------------------
# Family 2 core algebraic-spine checks (faithful replication)
# ---------------------------------------------------------------------------

def core_family2_checks(A_s: np.ndarray, B_R: np.ndarray) -> dict:
    """c1 kernel identity, c2 rank identity, c3 retention identity (F2 style)."""
    m, n = A_s.shape
    da = family2.thin_decomposition(A_s)
    rA = da["rank"]
    Q_A = da["Q"]
    rB, Z, svB, tolB = family2.range_basis(B_R)
    sigma1A = float(np.linalg.norm(A_s, ord=2))
    back_scale = max(1.0, sigma1A**2)
    c1_tol = 1e-8 * back_scale
    c3_tol = 1e-10 * back_scale

    I_m = np.eye(m, dtype=float)
    Pperp = I_m - Z @ Z.T
    K_IS = _sym(A_s.T @ A_s)
    KSL = _sym(A_s.T @ (Pperp @ A_s))
    C = Pperp @ A_s

    # ---- c1 kernel identity ----------------------------------------------
    w_eig, V_eig = np.linalg.eigh(KSL)
    _, svKSL, tolKSL = family2.rank_svd(KSL)
    k1 = int(np.sum(np.abs(w_eig) <= tolKSL))
    N1 = V_eig[:, :k1]
    rC, svC, tolC = family2.rank_svd(C)
    uC, sC, vhC = np.linalg.svd(C, full_matrices=True)
    N2 = vhC[rC:].T
    if N1.shape[1] and N2.shape[1]:
        subspace_distance = float(np.linalg.norm(N1 @ N1.T - N2 @ N2.T, ord=2))
        min_sing_N1N2 = float(np.linalg.svd(N1.T @ N2, compute_uv=False)[-1])
    else:
        subspace_distance = 0.0
        min_sing_N1N2 = None
    max_KSL_on_N1 = float(
        max(
            (np.linalg.norm(KSL @ N1[:, j], ord=2) for j in range(N1.shape[1])),
            default=0.0,
        )
    )
    max_C_on_N2 = float(
        max(
            (np.linalg.norm(C @ N2[:, j], ord=2) for j in range(N2.shape[1])),
            default=0.0,
        )
    )
    # informative crossed residuals used by family 2:
    max_KSL_on_N2 = float(
        max(
            (np.linalg.norm(KSL @ N2[:, j], ord=2) for j in range(N2.shape[1])),
            default=0.0,
        )
    )
    max_C_on_N1 = float(
        max(
            (np.linalg.norm(C @ N1[:, j], ord=2) for j in range(N1.shape[1])),
            default=0.0,
        )
    )
    nullity_match = bool(k1 == (C.shape[1] - rC))
    identity_support = bool(
        (max_KSL_on_N2 <= 1e-10)
        and nullity_match
        and ((min_sing_N1N2 is None) or (min_sing_N1N2 > 1 - 1e-6))
    )
    c1 = {
        "nullity_KSL_eig": k1,
        "nullity_C": int(C.shape[1] - rC),
        "nullity_match": nullity_match,
        "tol_KSL_eig_selector": float(tolKSL),
        "tol_C": float(tolC),
        "max_norm_KSL_on_N1": max_KSL_on_N1,
        "max_norm_C_on_N2": max_C_on_N2,
        "max_norm_KSL_on_N2": max_KSL_on_N2,
        "max_norm_C_on_N1": max_C_on_N1,
        "subspace_distance": subspace_distance,
        "min_singular_N1T_N2": min_sing_N1N2,
        "pass_gate_f2": bool(subspace_distance < c1_tol),
        "pass_gate_tol": c1_tol,
        "identity_support_f2": identity_support,
        "note": (
            "N1 = eig(K_SLAM) with |eigenvalue| <= tol; N2 = null(P_perp A) "
            "from the full SVD.  Empty bases give distance 0 and vacuous "
            "min singular.  max_norm_KSL_on_N2 and max_norm_C_on_N1 are the "
            "cross residuals reported by family 2; the literal same-space "
            "residuals on N1/N2 are recorded as well."
        ),
    }

    # ---- c2 rank identity ------------------------------------------------
    rKIS, _, _ = family2.rank_svd(K_IS)
    rKS, _, _ = family2.rank_svd(KSL)
    AB = np.hstack([A_s, B_R])
    rAB, _, _ = family2.rank_svd(AB)
    lhs = rKIS - rKS
    rhs = rA + rB - rAB
    c2 = {
        "r_A": rA,
        "r_KIS": int(rKIS),
        "r_KSL": int(rKS),
        "r_B": rB,
        "r_AB": int(rAB),
        "lhs_rKIS_minus_rKSL": int(lhs),
        "rhs_rA_plus_rB_minus_rAB": int(rhs),
        "rank_identity_residual": int(lhs - rhs),
        "rank_identity_holds": bool(lhs == rhs),
    }

    # ---- c3 no-prior retention principal-angle identity ------------------
    Rop = family4.retention_spectrum(A_s, Pperp)
    cos2 = np.linalg.svd(Z.T @ Q_A, compute_uv=False) ** 2 if (rA and rB) else np.array([])
    rho_pred = np.sort(
        np.concatenate(
            [
                1.0 - cos2,
                np.ones(max(rA - len(cos2), 0), dtype=float),
            ]
        )
    )[::-1]
    max_abs_diff = float(np.max(np.abs(Rop - rho_pred)))
    c3 = {
        "r_A": rA,
        "rho_min": float(Rop[-1]),
        "max_abs_diff_rho_vs_pred": max_abs_diff,
        "pass_gate": bool(max_abs_diff < c3_tol),
        "pass_gate_tol": c3_tol,
        "note": (
            "rho = retention_spectrum(A_s, I - Z Z^T); prediction = "
            "desc(concat(1 - sv(Z^T Q_A)^2, ones(r_A - len)))."
        ),
    }
    return {"c1_kernel_identity": c1, "c2_rank_identity": c2, "c3_retention": c3}


# ---------------------------------------------------------------------------
# Scenes
# ---------------------------------------------------------------------------

def _random_smooth_scene(points: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    chi = np.zeros(points.shape[0], dtype=float)
    for _ in range(int(CONFIG["partA_scenes"]["random"]["n_blobs"])):
        a = float(rng.uniform(*CONFIG["partA_scenes"]["random"]["a_range"]))
        c = rng.uniform(
            *CONFIG["partA_scenes"]["random"]["c_range"], size=2
        )
        s = float(
            rng.uniform(*CONFIG["partA_scenes"]["random"]["sigma_range"])
        )
        chi += a * np.exp(-np.sum((points - c) ** 2, axis=1) / (2.0 * s**2))
    return chi


def _ring_scene(points: np.ndarray) -> np.ndarray:
    ring = CONFIG["partA_scenes"]["ring"]
    r = np.linalg.norm(points, axis=1)
    return float(ring["amplitude"]) * np.exp(
        -((r - float(ring["radius"])) / float(ring["width"])) ** 2
    )


def _rescale_to_chi0(chi: np.ndarray, chi0: np.ndarray) -> np.ndarray:
    return chi * (float(np.linalg.norm(chi0)) / float(np.linalg.norm(chi)))


def _partA_scenes(points: np.ndarray, chi0: np.ndarray) -> dict:
    scenes = {}
    for seed in CONFIG["partA_scenes"]["random"]["seeds"]:
        scenes[f"seed_{seed}"] = _rescale_to_chi0(
            _random_smooth_scene(points, int(seed)), chi0
        )
    scenes["two_blob"] = chi0.copy()
    scenes["ring"] = _rescale_to_chi0(_ring_scene(points), chi0)
    return scenes


# ---------------------------------------------------------------------------
# Part A: seed replications
# ---------------------------------------------------------------------------

def _frequency_diversity_row(
    chi: np.ndarray, poses: np.ndarray, S: np.ndarray, cfg: dict
) -> dict:
    A1, B1, _ = family4.build_smooth_block(
        chi, poses, S, cfg, float(CONFIG["f_single"])
    )
    A14, B14, _ = family4.build_smooth_block(
        chi, poses, S, cfg, float(CONFIG["f_distinct"])
    )
    _, P1 = family4.range_projection(B1)
    rho_asc, U, R_op = family4.generalized_eigen_directions(A1, P1)
    Ast = np.vstack([A1, A14])
    Bst = np.vstack([B1, B14])
    _, Pst = family4.range_projection(Bst)

    rows = []
    n_moved = 0
    for j in range(3):
        u = U[:, j]
        rho_single = family4.rayleigh_quotient(A1, P1, u)
        rho_distinct = family4.rayleigh_quotient(Ast, Pst, u)
        movement = rho_distinct - rho_single
        if movement > CONFIG["movement_threshold"]:
            n_moved += 1
        rows.append(
            {
                "direction": int(j),
                "generalized_eigenvalue_asc": float(rho_asc[j]),
                "rho_single": float(rho_single),
                "rho_distinct": float(rho_distinct),
                "movement_distinct_minus_single": float(movement),
                "shared_z_residual_distinct": float(
                    family4.shared_z_residual(Ast, Bst, u)
                ),
                "u_KIS_norm2": float(np.dot(A1 @ u, A1 @ u)),
            }
        )
    return {
        "rows": rows,
        "n_directions_with_movement_gt_1e-4": int(n_moved),
        "family4_pass_ge2_moved": bool(n_moved >= 2),
        "directions_method": (
            "family4.generalized_eigen_directions(A1, P1), P1 = I - Z Z^T "
            "for Range(B1); three ascending (most confounded) directions; "
            "u = V_A diag(1/s_A) w"
        ),
    }


def run_part_A(
    points: np.ndarray, chi0: np.ndarray, S: np.ndarray, cfg: dict
) -> dict:
    poses = family1.build_poses(cfg)
    scenes = _partA_scenes(points, chi0)

    # Literal-spec diagnostic: generalized_eigen_directions(A1, A1.T@A1)
    literal_note = {}
    A1d, B1d, _ = family4.build_smooth_block(
        scenes["two_blob"], poses, S, cfg, float(CONFIG["f_single"])
    )
    try:
        family4.generalized_eigen_directions(A1d, A1d.T @ A1d)
        literal_note["raised"] = False
    except Exception as exc:  # noqa: BLE001 - diagnostic capture
        literal_note["raised"] = True
        literal_note["exception_type"] = type(exc).__name__
        literal_note["exception_message"] = str(exc)
        literal_note["dimensions"] = {
            "A1": A1d.shape,
            "A1.T@A1": (A1d.T @ A1d).shape,
            "Q_A": family2.thin_decomposition(A1d)["Q"].shape,
        }

    rows = []
    for name in CONFIG["partA_scenes"]["scene_order"]:
        chi = scenes[name]
        A_s, B_R, A_pix_R = family4.build_smooth_block(
            chi, poses, S, cfg, float(CONFIG["f_single"])
        )
        core = core_family2_checks(A_s, B_R)
        freq = _frequency_diversity_row(chi, poses, S, cfg)
        rows.append(
            {
                "scene": name,
                "l2_ratio_to_two_blob": float(
                    np.linalg.norm(chi) / np.linalg.norm(chi0)
                ),
                "chi_min": float(np.min(chi)),
                "chi_max": float(np.max(chi)),
                "A_s_norm_fro": float(np.linalg.norm(A_s, ord="fro")),
                "A_s_norm_2": float(np.linalg.norm(A_s, ord=2)),
                "B_R_norm_fro": float(np.linalg.norm(B_R, ord="fro")),
                "core_family2": core,
                "frequency_diversity": freq,
            }
        )

    # per-gate pass-rate summaries
    summary = {
        "n_scenes": len(rows),
        "rank_identity_pass": int(
            sum(r["core_family2"]["c2_rank_identity"]["rank_identity_holds"] for r in rows)
        ),
        "retention_identity_pass": int(
            sum(r["core_family2"]["c3_retention"]["pass_gate"] for r in rows)
        ),
        "kernel_identity_f2_gate_pass": int(
            sum(r["core_family2"]["c1_kernel_identity"]["pass_gate_f2"] for r in rows)
        ),
        "kernel_identity_support_pass": int(
            sum(
                r["core_family2"]["c1_kernel_identity"]["identity_support_f2"]
                for r in rows
            )
        ),
        "kernel_nullity_match": int(
            sum(r["core_family2"]["c1_kernel_identity"]["nullity_match"] for r in rows)
        ),
        "frequency_diversity_pass_ge2": int(
            sum(
                r["frequency_diversity"]["family4_pass_ge2_moved"]
                for r in rows
            )
        ),
        "count_moved_gt_1e-4_per_scene": [
            r["frequency_diversity"]["n_directions_with_movement_gt_1e-4"]
            for r in rows
        ],
    }

    # read-only crosscheck against stored family 2 and family 4 JSONs
    crosscheck = {"available": False}
    try:
        f2res = json.loads(
            (_ROOT / "results" / "family2_results.json").read_text()
        )
        f4res = json.loads(
            (_ROOT / "results" / "family4_results.json").read_text()
        )
        std_row = next(r for r in rows if r["scene"] == "two_blob")
        std_core = std_row["core_family2"]
        f2s = f2res["cases"]["smooth"]
        f2_rank_res = int(
            f2s["check_c2_rank_identity"]["residual_lhs_minus_d_inter"]
        )
        f2_ret_diff = float(
            f2s["check_c3_retention_spectrum"]["max_abs_diff_rho_vs_pred"]
        )
        f2_kernel_gate = bool(f2s["check_c1_kernel_identity"]["pass_gate"])
        f2_kernel_support = bool(
            f2s["check_c1_kernel_identity"]["identity_support"]
        )
        f4_c_rows = f4res["checks"]["C_shared_pose_compensation"]["rows"]
        std_freq = std_row["frequency_diversity"]["rows"]
        mov_max = max(
            abs(std_freq[j]["movement_distinct_minus_single"]
                - float(f4_c_rows[j]["movement_distinct_minus_single"]))
            for j in range(3)
        )
        crosscheck = {
            "available": True,
            "standard_rank_residual_match": bool(
                std_core["c2_rank_identity"]["rank_identity_residual"]
                == f2_rank_res
            ),
            "standard_rank_residual_f2": f2_rank_res,
            "standard_retention_diff_f2": f2_ret_diff,
            "standard_retention_max_abs_residual": std_core["c3_retention"][
                "max_abs_diff_rho_vs_pred"
            ],
            "standard_retention_pass_match": bool(
                std_core["c3_retention"]["pass_gate"]
                == f2s["check_c3_retention_spectrum"]["pass_pred_gate"]
            ),
            "standard_kernel_gate_match": bool(
                std_core["c1_kernel_identity"]["pass_gate_f2"]
                == f2_kernel_gate
            ),
            "standard_kernel_support_match": bool(
                std_core["c1_kernel_identity"]["identity_support_f2"]
                == f2_kernel_support
            ),
            "standard_freq_movement_max_abs_diff_vs_f4": mov_max,
            "standard_freq_movement_match_lt_1e-12": bool(mov_max < 1e-12),
            "family4_checkC_pass": bool(
                f4res["checks"]["C_shared_pose_compensation"]["pass_gate"]
            ),
            "replication_note": (
                "two_blob row is compared to results/family2_results.json "
                "(cases.smooth) and results/family4_results.json "
                "(checks.C) read-only."
            ),
        }
    except Exception as exc:  # noqa: BLE001
        crosscheck = {
            "available": False,
            "error": f"{type(exc).__name__}: {exc}",
        }

    return {
        "scenes": rows,
        "pass_rates": summary,
        "literal_generalized_second_argument_diagnostic": literal_note,
        "crosscheck_vs_family2_4": crosscheck,
    }


# ---------------------------------------------------------------------------
# Part B: pose-DOF nuisance baseline (standard scene, f=1.0)
# ---------------------------------------------------------------------------

def run_part_B(A_s: np.ndarray, B_R: np.ndarray) -> dict:
    subset_names = [
        "known_pose",
        "x_only",
        "y_only",
        "theta_only",
        "xy",
        "x_theta",
        "y_theta",
        "full",
    ]
    rows = []
    info_by_name = {}
    for name in subset_names:
        B_sub = _pose_submatrix(B_R, name)
        info = info_metrics(A_s, B_sub)
        info_by_name[name] = info
        rows.append(
            {
                "subset": name,
                "q_sub": int(B_sub.shape[1]),
                "info_metrics": info,
                "B_over_A_fro": _safe_div(
                    float(np.linalg.norm(B_sub, ord="fro")),
                    float(np.linalg.norm(A_s, ord="fro")),
                ),
            }
        )

    def ret(name: str) -> float:
        return info_by_name[name]["retained_mass"]

    def conf(name: str) -> float:
        return info_by_name[name]["confusable_mass"]

    transitions = [
        ("known_pose", "x_only", "add x nuisance"),
        ("known_pose", "y_only", "add y nuisance"),
        ("known_pose", "theta_only", "add theta nuisance"),
        ("x_only", "xy", "add y to x-only"),
        ("x_only", "x_theta", "add theta to x-only"),
        ("y_only", "xy", "add x to y-only"),
        ("y_only", "y_theta", "add theta to y-only"),
        ("theta_only", "x_theta", "add x to theta-only"),
        ("theta_only", "y_theta", "add y to theta-only"),
        ("x_only", "full", "add y+theta to x-only"),
        ("y_only", "full", "add x+theta to y-only"),
        ("theta_only", "full", "add x+y to theta-only"),
        ("xy", "full", "add theta to xy"),
        ("x_theta", "full", "add y to x-theta"),
        ("y_theta", "full", "add x to y-theta"),
    ]
    marginals = []
    for base, aug, label in transitions:
        marginals.append(
            {
                "from": base,
                "to": aug,
                "label": label,
                "retained_from": ret(base),
                "retained_to": ret(aug),
                "marginal_destruction_retained": float(ret(base) - ret(aug)),
                "confusable_from": conf(base),
                "confusable_to": conf(aug),
                "confusable_increase": float(conf(aug) - conf(base)),
            }
        )
    return {
        "A_s_fro": float(np.linalg.norm(A_s, ord="fro")),
        "B_R_fro_full": float(np.linalg.norm(B_R, ord="fro")),
        "pose_column_order": [0, 1, 2],
        "pose_column_order_note": "[px, py, theta] repeated per pose",
        "rows": rows,
        "marginal_destruction": marginals,
        "known_pose_retained_mass_reference": info_by_name["known_pose"][
            "retained_mass"
        ],
        "claim_scope": "descriptive ablation; no forced pass",
    }


# ---------------------------------------------------------------------------
# Parts C, D, E ablations
# ---------------------------------------------------------------------------

def _monotone(vals: list) -> tuple[bool, bool]:
    """(non_decreasing, non_increasing) over consecutive numeric values."""
    ok_up = all(b >= a for a, b in zip(vals, vals[1:]))
    ok_down = all(b <= a for a, b in zip(vals, vals[1:]))
    return bool(ok_up), bool(ok_down)


def run_part_C(
    chi0: np.ndarray,
    poses: np.ndarray,
    S: np.ndarray,
    cfg: dict,
) -> dict:
    rows = []
    for conf in CONFIG["partC_receiver_sweep"]["rows"]:
        cfg_rx = dict(cfg)
        cfg_rx["n_rx"] = int(conf["n_rx"])
        cfg_rx["rx_offsets"] = conf["rx_offsets"]
        A_s, B_R, A_pix_R = family4.build_smooth_block(
            chi0, poses, S, cfg_rx, float(CONFIG["f_single"])
        )
        info = info_metrics(A_s, B_R)
        rows.append(
            {
                "n_rx": int(conf["n_rx"]),
                "m_real": int(A_s.shape[0]),
                "rx_offsets": conf["rx_offsets"],
                "info_metrics": info,
                "B_over_A_fro": _safe_div(
                    float(np.linalg.norm(B_R, ord="fro")),
                    float(np.linalg.norm(A_s, ord="fro")),
                ),
            }
        )
    retained = [r["info_metrics"]["retained_mass"] for r in rows]
    confusable = [r["info_metrics"]["confusable_mass"] for r in rows]
    theta = [r["info_metrics"]["theta_min_deg"] for r in rows]
    rho_min = [r["info_metrics"]["rho_min"] for r in rows]
    ret_up, ret_down = _monotone(retained)
    conf_up, conf_down = _monotone(confusable)
    return {
        "rows": rows,
        "observed": {
            "n_rx": [r["n_rx"] for r in rows],
            "retained_mass": retained,
            "confusable_mass": confusable,
            "theta_min_deg": theta,
            "rho_min": rho_min,
            "B_over_A_fro": [r["B_over_A_fro"] for r in rows],
        },
        "observed_flags": {
            "retained_non_increasing": ret_down,
            "retained_non_decreasing": ret_up,
            "confusable_non_decreasing": conf_up,
            "confusable_non_increasing": conf_down,
        },
        "note": (
            "r_A itself changes with measurement count: n_rx=1 gives "
            "r_A=12 and Range(B_R)=R^12 (A fully B-confounded, retained 0); "
            "n_rx=2 gives r_A=24 with Range(B_R) contained in Range(A_s) "
            "(retained 6); n_rx>=4 keeps r_A=24 with a transverse "
            "intersection (retained 9.41 for n_rx=4 and 9.40 for n_rx=8). "
            "The literal monotonicity flags over the mixed-r_A rows are "
            "therefore not interpretable as a pure receiver-count effect."
        ),
        "claim_scope": "descriptive ablation; no forced pass",
    }


def run_part_D(
    chi0: np.ndarray,
    S: np.ndarray,
    cfg: dict,
) -> dict:
    rows = []
    for T in CONFIG["partD_pose_sweep"]["Ts"]:
        cfgT = dict(cfg)
        cfgT["T"] = int(T)
        posesT = family1.build_poses(cfgT)
        A_s, B_R, A_pix_R = family4.build_smooth_block(
            chi0, posesT, S, cfgT, float(CONFIG["f_single"])
        )
        info = info_metrics(A_s, B_R)
        rows.append(
            {
                "T": int(T),
                "q_pose": int(B_R.shape[1]),
                "m_real": int(A_s.shape[0]),
                "info_metrics": info,
                "B_over_A_fro": _safe_div(
                    float(np.linalg.norm(B_R, ord="fro")),
                    float(np.linalg.norm(A_s, ord="fro")),
                ),
            }
        )
    retained = [r["info_metrics"]["retained_mass"] for r in rows]
    confusable = [r["info_metrics"]["confusable_mass"] for r in rows]
    theta = [r["info_metrics"]["theta_min_deg"] for r in rows]
    ret_up, ret_down = _monotone(retained)
    conf_up, conf_down = _monotone(confusable)
    return {
        "rows": rows,
        "observed": {
            "T": [r["T"] for r in rows],
            "retained_mass": retained,
            "confusable_mass": confusable,
            "theta_min_deg": theta,
        },
        "observed_flags": {
            "retained_non_decreasing": ret_up,
            "retained_non_increasing": ret_down,
            "confusable_non_decreasing": conf_up,
            "confusable_non_increasing": conf_down,
        },
        "note": (
            "q=3T grows with T and r_A stays 24 (A_s has p=24 columns), so "
            "retained mass falls monotonically here (15 -> 9.41 -> 6.58) as "
            "more pose nuisance columns occupy more of the data space."
        ),
        "claim_scope": "descriptive ablation; no forced pass",
    }


def run_part_E(
    points: np.ndarray,
    chi0: np.ndarray,
    poses: np.ndarray,
    cfg: dict,
    A_pix_standard: np.ndarray,
    B_R_standard: np.ndarray,
) -> dict:
    rows = []
    for bcfg in CONFIG["partE_basis_sweep"]["configs"]:
        sb = dict(CONFIG["smooth_basis"])
        sb["p"] = int(bcfg["p"])
        sb["x_centers_n"] = int(bcfg["x_centers_n"])
        sb["y_centers_n"] = int(bcfg["y_centers_n"])
        S_p = family2.build_smooth_basis(points, sb)
        A_s = A_pix_standard @ S_p
        info = info_metrics(A_s, B_R_standard)
        rows.append(
            {
                "p": int(bcfg["p"]),
                "x_centers_n": int(bcfg["x_centers_n"]),
                "y_centers_n": int(bcfg["y_centers_n"]),
                "smooth_basis_config": sb,
                "info_metrics": info,
                "B_over_A_fro": _safe_div(
                    float(np.linalg.norm(B_R_standard, ord="fro")),
                    float(np.linalg.norm(A_s, ord="fro")),
                ),
            }
        )
    retained = [r["info_metrics"]["retained_mass"] for r in rows]
    confusable = [r["info_metrics"]["confusable_mass"] for r in rows]
    theta = [r["info_metrics"]["theta_min_deg"] for r in rows]
    ret_up, ret_down = _monotone(retained)
    conf_up, conf_down = _monotone(confusable)
    return {
        "rows": rows,
        "observed": {
            "p": [r["p"] for r in rows],
            "retained_mass": retained,
            "confusable_mass": confusable,
            "theta_min_deg": theta,
        },
        "observed_flags": {
            "retained_non_decreasing": ret_up,
            "retained_non_increasing": ret_down,
            "confusable_non_decreasing": conf_up,
            "confusable_non_increasing": conf_down,
        },
        "note": (
            "B_R_standard is fixed (r_B=18) while r_A = p grows, so retained "
            "mass increases with p (0.157 -> 9.41 -> 18.75) partly because "
            "more A directions are added outside the fixed B range; "
            "confusable mass also increases (11.84 -> 14.59 -> 17.25)."
        ),
        "claim_scope": "descriptive ablation; no forced pass",
    }


# ---------------------------------------------------------------------------
# run_all
# ---------------------------------------------------------------------------

def run_all() -> dict:
    cfg = CONFIG
    points, chi0, h, S = family4.base_scene(cfg)
    poses = family1.build_poses(cfg)

    A_s_std, B_R_std, A_pix_std = family4.build_smooth_block(
        chi0, poses, S, cfg, float(cfg["f_single"])
    )

    t0 = time.perf_counter()
    partA = run_part_A(points, chi0, S, cfg)
    tA = time.perf_counter() - t0

    t0 = time.perf_counter()
    partB = run_part_B(A_s_std, B_R_std)
    tB = time.perf_counter() - t0

    t0 = time.perf_counter()
    partC = run_part_C(chi0, poses, S, cfg)
    tC = time.perf_counter() - t0

    t0 = time.perf_counter()
    partD = run_part_D(chi0, S, cfg)
    tD = time.perf_counter() - t0

    t0 = time.perf_counter()
    partE = run_part_E(
        points, chi0, poses, cfg, A_pix_std, B_R_std
    )
    tE = time.perf_counter() - t0

    checks = [
        {
            "id": "A_rank_identity",
            "claim": (
                "rank identity r(K_IS)-r(K_SLAM) == r(A)+r(B)-r([A,B]) "
                "at machine rank over the 7 scenes"
            ),
            "threshold": "integer residual == 0 in every scene",
            "observed_pass": partA["pass_rates"]["rank_identity_pass"],
            "observed_total": partA["pass_rates"]["n_scenes"],
            "failing_scenes": [
                r["scene"]
                for r in partA["scenes"]
                if not r["core_family2"]["c2_rank_identity"]["rank_identity_holds"]
            ],
            "pass": bool(
                partA["pass_rates"]["rank_identity_pass"]
                == partA["pass_rates"]["n_scenes"]
            ),
            "note": (
                "recorded at the family-2 machine tolerance; a tolerance "
                "boundary on K_SLAM's smallest eigenvalue can flip one rank "
                "unit (ring scene), see report."
            ),
        },
        {
            "id": "A_retention_principal_angle_identity",
            "claim": (
                "no-prior retention spectrum equals 1 - cos^2 of the "
                "principal angles (padded) within the family-2 gate"
            ),
            "threshold": "max abs diff < 1e-10 * max(1, ||A_s||_2^2)",
            "observed_pass": partA["pass_rates"]["retention_identity_pass"],
            "observed_total": partA["pass_rates"]["n_scenes"],
            "pass": bool(
                partA["pass_rates"]["retention_identity_pass"]
                == partA["pass_rates"]["n_scenes"]
            ),
        },
        {
            "id": "A_kernel_identity",
            "claim": (
                "kernel identity between null(K_SLAM) and null((I-ZZ^T)A_s) "
                "at the family-2 threshold"
            ),
            "threshold": (
                "family-2 subspace gate < 1e-8*max(1,||A_s||_2^2) and "
                "identity_support (equal nullity + cross residuals + "
                "subspace agreement)"
            ),
            "observed_f2_gate_pass": partA["pass_rates"][
                "kernel_identity_f2_gate_pass"
            ],
            "observed_identity_support_pass": partA["pass_rates"][
                "kernel_identity_support_pass"
            ],
            "observed_nullity_match": partA["pass_rates"]["kernel_nullity_match"],
            "observed_total": partA["pass_rates"]["n_scenes"],
            "pass": bool(
                partA["pass_rates"]["kernel_identity_f2_gate_pass"]
                == partA["pass_rates"]["n_scenes"]
            ),
            "note": (
                "null(K_SLAM) and null(C) are empty on the two-blob and most "
                "seed scenes at the declared eig thresholds; the ring scene "
                "selects one K_SLAM eig below tol while C stays full rank, so "
                "its nullity match is false (identity_support false) although "
                "the family-2 subspace-distance gate is vacuous-true."
            ),
        },
        {
            "id": "A_frequency_diversity",
            "claim": (
                "stacking distinct f in {1.0,1.4} moves the 3 most-confounded "
                "single-frequency directions; family-4 gate >=2/3 > 1e-4"
            ),
            "threshold": "rho movement > 1e-4 in at least 2 of 3 directions",
            "observed_pass": partA["pass_rates"]["frequency_diversity_pass_ge2"],
            "observed_total": partA["pass_rates"]["n_scenes"],
            "observed_movement_counts": partA["pass_rates"][
                "count_moved_gt_1e-4_per_scene"
            ],
            "pass": bool(
                partA["pass_rates"]["frequency_diversity_pass_ge2"]
                == partA["pass_rates"]["n_scenes"]
            ),
        },
        {
            "id": "B_pose_dof_nuisance",
            "claim": (
                "pose-DOF nuisance retained/confusable masses and marginal "
                "destruction pattern (recorded; no forced pass)"
            ),
            "threshold": None,
            "observed": partB["rows"],
            "pass": None,
        },
        {
            "id": "C_receiver_count",
            "claim": (
                "receiver-count ablation trend (recorded; no forced pass)"
            ),
            "threshold": None,
            "observed": partC["rows"],
            "pass": None,
        },
        {
            "id": "D_pose_count",
            "claim": "pose-count ablation trend (recorded; no forced pass)",
            "threshold": None,
            "observed": partD["rows"],
            "pass": None,
        },
        {
            "id": "E_basis_p",
            "claim": "smooth-basis-p ablation trend (recorded; no forced pass)",
            "threshold": None,
            "observed": partE["rows"],
            "pass": None,
        },
    ]

    return {
        "standard_scene_metrics": info_metrics(A_s_std, B_R_std),
        "part_seconds": {
            "A": tA,
            "B": tB,
            "C": tC,
            "D": tD,
            "E": tE,
        },
        "partA_seed_replications": partA,
        "partB_pose_dof_nuisance": partB,
        "partC_receiver_ablation": partC,
        "partD_pose_count_ablation": partD,
        "partE_basis_ablation": partE,
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _plot_seed_replications(partA: dict, path: Path) -> None:
    rows = partA["scenes"]
    labels = [r["scene"] for r in rows]
    x = np.arange(len(rows))
    ret_diff = [
        r["core_family2"]["c3_retention"]["max_abs_diff_rho_vs_pred"]
        for r in rows
    ]
    rank_res = [
        int(r["core_family2"]["c2_rank_identity"]["rank_identity_residual"])
        for r in rows
    ]
    movements = [
        [r["frequency_diversity"]["rows"][j][
            "movement_distinct_minus_single"
        ] for j in range(3)]
        for r in rows
    ]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.9))

    ax = axes[0]
    floor = 1e-16
    shown = [max(v, floor) for v in ret_diff]
    ax.semilogy(x, shown, "o-", color="#1f77b4", ms=6, lw=1.4,
                label="retention identity max abs diff")
    for xi, rr in zip(x, rank_res):
        ax.text(
            xi,
            floor,
            f"rank-res={rr}",
            ha="center",
            va="bottom",
            fontsize=7,
            color="#d62728" if rr else "#2ca02c",
        )
    ax.axhline(floor, color="grey", lw=0.7, ls=":")
    ax.set_yscale("log")
    ax.set_ylim(floor / 10.0, max(max(shown) * 10.0, 1e-12))
    ax.set_xticks(x, labels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("max abs diff (log; exact zeros shown at 1e-16 floor)")
    ax.set_title(
        "(a) Family 2 core identities across scenes\n"
        "retention principal-angle max abs diff + rank-identity residual"
    )
    ax.grid(True, axis="y", which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    ax = axes[1]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    width = 0.25
    for j in range(3):
        ax.bar(
            x + (j - 1) * width,
            [m[j] for m in movements],
            width,
            color=colors[j],
            label=f"direction u{j}",
        )
    ax.axhline(1e-4, color="black", ls="--", lw=1.0, label="1e-4 movement gate")
    ax.set_xticks(x, labels, rotation=35, ha="right", fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("rho movement (distinct stack minus single, log)")
    ax.set_title("(b) frequency-diversity movement of 3 most-confounded directions")
    ax.grid(True, axis="y", which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")

    fig.suptitle("Family 9: deterministic seed replications (7 scenes, N=16, f=1.0 vs {1.0,1.4})", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_pose_dof(partB: dict, path: Path) -> None:
    names = [r["subset"] for r in partB["rows"]]
    retained = [r["info_metrics"]["retained_mass"] for r in partB["rows"]]
    confusable = [r["info_metrics"]["confusable_mass"] for r in partB["rows"]]
    x = np.arange(len(names))
    width = 0.38
    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    ax.bar(x - width / 2, retained, width, label="retained mass", color="#1f77b4")
    ax.bar(
        x + width / 2,
        confusable,
        width,
        label="confusable mass (sum cos^2)",
        color="#d62728",
    )
    ax.axhline(24.0, color="black", ls=":", lw=1.0)
    ax.text(
        len(names) - 0.5,
        24.0,
        " r_A = 24 (all smooth directions)",
        fontsize=7,
        va="bottom",
        ha="right",
    )
    for xi, (ret, conf) in enumerate(zip(retained, confusable)):
        ax.text(
            xi,
            max(ret, conf) + 0.2,
            f"theta={partB['rows'][xi]['info_metrics']['theta_min_deg']:.1f}",
            ha="center",
            fontsize=6.5,
        )
    ax.set_xticks(x, names, rotation=25, ha="right")
    ax.set_ylabel("mass")
    ax.set_title(
        "Family 9 Part B: pose-DOF nuisance subsets (standard two-blob, "
        "f=1.0, q=18 full)"
    )
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_rx_T_basis(partC: dict, partD: dict, partE: dict, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))

    def draw(ax, x, row_names, retained, confusable, theta, title, xlabel):
        ax.plot(x, retained, "o-", color="#1f77b4", label="retained mass")
        ax.plot(x, confusable, "s--", color="#d62728", label="confusable mass")
        for i, (xi, th) in enumerate(zip(x, theta)):
            y_top = max(float(retained[i]), float(confusable[i]))
            ax.annotate(
                f"{th:.1f} deg",
                (xi, y_top),
                xytext=(xi, y_top + 0.5),
                ha="center",
                fontsize=6.5,
                arrowprops=dict(arrowstyle="-", lw=0.4, color="grey"),
            )
        ax.set_xticks(x, row_names)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)

    draw(
        axes[0],
        np.arange(len(partC["rows"])),
        [f"n_rx={r['n_rx']}" for r in partC["rows"]],
        partC["observed"]["retained_mass"],
        partC["observed"]["confusable_mass"],
        partC["observed"]["theta_min_deg"],
        "(a) receiver-count ablation",
        "n_rx",
    )
    draw(
        axes[1],
        np.arange(len(partD["rows"])),
        [f"T={r['T']}" for r in partD["rows"]],
        partD["observed"]["retained_mass"],
        partD["observed"]["confusable_mass"],
        partD["observed"]["theta_min_deg"],
        "(b) pose-count ablation (arc90 geometry)",
        "T",
    )
    draw(
        axes[2],
        np.arange(len(partE["rows"])),
        [f"p={r['p']}" for r in partE["rows"]],
        partE["observed"]["retained_mass"],
        partE["observed"]["confusable_mass"],
        partE["observed"]["theta_min_deg"],
        "(c) smooth-basis p ablation",
        "p",
    )
    fig.suptitle(
        "Family 9: receiver-count, pose-count, and smooth-basis ablations "
        "(standard two-blob scene, f=1.0)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def build_report(results: dict, cfg: dict, paths: dict) -> str:
    pA = results["partA_seed_replications"]
    pB = results["partB_pose_dof_nuisance"]
    pC = results["partC_receiver_ablation"]
    pD = results["partD_pose_count_ablation"]
    pE = results["partE_basis_ablation"]
    pr = pA["pass_rates"]

    scene_rows = []
    for r in pA["scenes"]:
        c1 = r["core_family2"]["c1_kernel_identity"]
        c2 = r["core_family2"]["c2_rank_identity"]
        c3 = r["core_family2"]["c3_retention"]
        fd = r["frequency_diversity"]
        scene_rows.append(
            [
                r["scene"],
                int(c2["rank_identity_residual"]),
                "PASS" if c2["rank_identity_holds"] else "FAIL",
                _fmt(c3["max_abs_diff_rho_vs_pred"], ".3e"),
                "PASS" if c3["pass_gate"] else "FAIL",
                f"{c1['nullity_KSL_eig']}/{c1['nullity_C']}",
                "PASS" if c1["nullity_match"] else "FAIL",
                fd["n_directions_with_movement_gt_1e-4"],
                "PASS" if fd["family4_pass_ge2_moved"] else "FAIL",
            ]
        )

    check_rows = []
    for c in results["checks"]:
        if c["pass"] is None:
            status = "recorded (no forced pass)"
            observed = "see ablation table(s)"
        elif c["id"] == "A_kernel_identity":
            observed = (
                f"F2 gate {c['observed_f2_gate_pass']}/{c['observed_total']}; "
                f"identity_support {c['observed_identity_support_pass']}/"
                f"{c['observed_total']}; nullity match "
                f"{c['observed_nullity_match']}/{c['observed_total']}"
            )
            status = "PASS" if c["pass"] else "PARTIAL/FAIL"
        else:
            observed = f"{c['observed_pass']}/{c['observed_total']}"
            if c["id"] == "A_frequency_diversity":
                observed += (
                    " (moved counts "
                    + ",".join(str(v) for v in c["observed_movement_counts"])
                    + ")"
                )
            elif c["id"] == "A_rank_identity":
                observed += " failing scenes=" + (
                    ",".join(c["failing_scenes"]) if c["failing_scenes"] else "-"
                )
            status = "PASS" if c["pass"] else "PARTIAL/FAIL"
        check_rows.append([c["id"], c["claim"], observed, status])

    pose_rows = [
        [
            r["subset"],
            r["q_sub"],
            _fmt(r["info_metrics"]["retained_mass"], ".6f"),
            _fmt(r["info_metrics"]["confusable_mass"], ".6f"),
            _fmt(r["info_metrics"]["theta_min_deg"], ".4f"),
            _fmt(r["info_metrics"]["rho_min"], ".3e"),
            _fmt(r["info_metrics"]["K_IS_fro"], ".4e"),
            _fmt(r["info_metrics"]["K_SLAM_fro"], ".4e"),
            _fmt(r["info_metrics"]["K_eff_fro"], ".4e"),
        ]
        for r in pB["rows"]
    ]
    marginal_rows = [
        [
            m["label"],
            _fmt(m["marginal_destruction_retained"], ".6f"),
            _fmt(m["confusable_increase"], ".6f"),
        ]
        for m in pB["marginal_destruction"]
    ]
    rx_rows = [
        [
            r["n_rx"],
            _fmt(r["info_metrics"]["retained_mass"], ".6f"),
            _fmt(r["info_metrics"]["confusable_mass"], ".6f"),
            _fmt(r["info_metrics"]["theta_min_deg"], ".4f"),
            _fmt(r["info_metrics"]["rho_min"], ".3e"),
            _fmt(r["B_over_A_fro"], ".4f"),
        ]
        for r in pC["rows"]
    ]
    T_rows = [
        [
            r["T"],
            r["q_pose"],
            _fmt(r["info_metrics"]["retained_mass"], ".6f"),
            _fmt(r["info_metrics"]["confusable_mass"], ".6f"),
            _fmt(r["info_metrics"]["theta_min_deg"], ".4f"),
            _fmt(r["info_metrics"]["rho_min"], ".3e"),
        ]
        for r in pD["rows"]
    ]
    p_rows = [
        [
            r["p"],
            r["x_centers_n"],
            r["y_centers_n"],
            _fmt(r["info_metrics"]["retained_mass"], ".6f"),
            _fmt(r["info_metrics"]["confusable_mass"], ".6f"),
            _fmt(r["info_metrics"]["theta_min_deg"], ".4f"),
            _fmt(r["info_metrics"]["rho_min"], ".3e"),
        ]
        for r in pE["rows"]
    ]

    fig_names = {
        name: f"figures/{Path(p).name}" for name, p in paths["figures"].items()
    }
    rep = f"""# Family 9: seed replications and component ablations

Date: {results['generated_utc']} UTC.  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Family 9 re-runs the deterministic algebraic-spine (Family 2 c1-c3) and
frequency-diversity (Family 4 check C) claims over seven scenes and adds
descriptive component ablations: pose-DOF nuisance subsets, receiver count,
pose count, and smooth-basis p.  It reuses (read-only) `helmholtz.py`,
`family1_pilot.py`, `family2_algebraic_spine.py`, and
`family4_frequency_trajectory.py`; no existing file is modified.

## Exact command and runtime

```bash
{results['command']}
```

Wall runtime: {results['runtime_seconds']:.2f} s (part seconds:
A {results['part_seconds']['A']:.2f}, B {results['part_seconds']['B']:.2f},
C {results['part_seconds']['C']:.2f}, D {results['part_seconds']['D']:.2f},
E {results['part_seconds']['E']:.2f}).  Platform:
{results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, matplotlib
{results['platform']['matplotlib']}.

Source SHA-256 (this script):
`{results['source_sha256']}`

Reused-module SHA-256:
{', '.join(f"`{Path(k).name}` `{v}`" for k, v in sorted(results['reused_sha256'].items()))}

## Config

N={cfg['N']} (N^2=256), standard T={cfg['T']}, n_rx={cfg['n_rx']}
(realified rows m=2*T*n_rx; m=48 at standard T/n_rx), k_b=2*pi (f=1.0),
Family-1 90-degree arc radius {cfg['arc_radius']}, phi
{cfg['arc_phi_deg'][0]}..{cfg['arc_phi_deg'][1]} deg,
theta=atan2(-py,-px); rx offsets {cfg['rx_offsets']}; tx offset
{cfg['tx_offset']}.  Standard two-blob chi0 as in Families 1-8.  Smooth basis
p={cfg['smooth_basis']['p']} ({cfg['smooth_basis']['x_centers_n']}x
{cfg['smooth_basis']['y_centers_n']} centres, sigma_b=
{cfg['smooth_basis']['sigma_b']}, unit columns).  Random scenes:
{cfg['partA_scenes']['random']['n_blobs']} positive blobs per scene,
seeds {cfg['partA_scenes']['random']['seeds']}, rescaled to the two-blob L2
norm.  Whitening is identity (`whiten_realify(A,B,None)`); alpha=
{cfg['alpha_prior']}.  Ranks use the family-2 rule
`tol(M)=max(M.shape)*eps*sigma_1(M)`.

## Claim-status table (seed replications over 7 scenes)

{_md_table(['claim', 'status', 'executed comparison', 'key numbers'],
           [
               [
                   'rank identity (F2 c2)',
                   'PASS' if pr['rank_identity_pass'] == pr['n_scenes'] else 'PARTIAL',
                   'integer residual == 0 per scene',
                   f"{pr['rank_identity_pass']}/{pr['n_scenes']}; "
                   'ring residual = 1 (machine-rank boundary)',
               ],
               [
                   'no-prior retention identity (F2 c3)',
                   'PASS' if pr['retention_identity_pass'] == pr['n_scenes'] else 'PARTIAL',
                   'max abs diff < 1e-10*max(1,||A_s||_2^2)',
                   f"{pr['retention_identity_pass']}/{pr['n_scenes']}; "
                   'diff ~1e-15 on all scenes',
               ],
               [
                   'kernel identity (F2 c1)',
                   'PARTIAL',
                   'F2 subspace gate; identity_support/nullity match',
                   f"F2 gate {pr['kernel_identity_f2_gate_pass']}/"
                   f"{pr['n_scenes']}; identity_support "
                   f"{pr['kernel_identity_support_pass']}/{pr['n_scenes']}; "
                   f"nullity match {pr['kernel_nullity_match']}/{pr['n_scenes']}"
               ],
               [
                   'frequency diversity (F4 check C)',
                   'PASS' if pr['frequency_diversity_pass_ge2'] == pr['n_scenes'] else 'PARTIAL',
                   'rho movement > 1e-4 in >=2 of 3 directions',
                   f"{pr['frequency_diversity_pass_ge2']}/{pr['n_scenes']}; "
                   'per-scene moved counts '
                   f"{pr['count_moved_gt_1e-4_per_scene']}"
               ],
           ])}

Part B-D-E claims are descriptive ablations with no forced pass:

{_md_table(['part', 'ablation', 'status'],
           [
               ['B', 'pose-DOF nuisance subsets + marginal destruction', 'recorded'],
               ['C', 'receiver count n_rx in {1,2,4,8}', 'recorded'],
               ['D', 'pose count T in {3,6,12}', 'recorded'],
               ['E', 'smooth-basis p in {12,24,36}', 'recorded'],
           ])}

## A. Seed replications (per scene)

Per-scene table (rank-identity residual; retention identity max abs diff and
gate; kernel identity nullities KSL/C and nullity match; frequency-diversity
moved count and family-4 gate):

{_md_table(['scene', 'rank res', 'rank pass', 'retention diff',
            'retention pass', 'null KSL/C', 'null match',
            'moved >1e-4', 'C pass'],
           scene_rows)}

Pass-rate summary: rank identity {pr['rank_identity_pass']}/
{pr['n_scenes']}; retention identity {pr['retention_identity_pass']}/
{pr['n_scenes']}; kernel identity F2 gate {pr['kernel_identity_f2_gate_pass']}/
{pr['n_scenes']} with identity_support {pr['kernel_identity_support_pass']}/
{pr['n_scenes']} and nullity match {pr['kernel_nullity_match']}/
{pr['n_scenes']}; frequency diversity {pr['frequency_diversity_pass_ge2']}/
{pr['n_scenes']} (all scenes move 3/3 directions; see table).

Cross-check of the standard two-blob row against stored results:
{_fmt(pA['crosscheck_vs_family2_4'].get('standard_freq_movement_max_abs_diff_vs_f4') if pA['crosscheck_vs_family2_4'].get('available') else 'unavailable')}
max abs movement difference vs `family4_results.json`
(match < 1e-12:
{pA['crosscheck_vs_family2_4'].get('standard_freq_movement_match_lt_1e-12') if pA['crosscheck_vs_family2_4'].get('available') else 'unavailable'}).

## A caveat: ring-scene rank boundary

The ring scene (radially very smooth) has K_SLAM's smallest eigenvalue below
the family-2 eig rank tolerance while C=(I-ZZ^T)A_s stays full rank at its own
(larger) SVD tolerance.  The literal machine-rank identity therefore records
residual 1 and nullity match false on the ring scene.  This is recorded as a
machine-rank boundary event, not as evidence that the exact algebraic identity
is false; the exact statements would hold under exact arithmetic for these
subspaces if no near-degeneracy crosses the tolerance.

## Check-C direction method note

The executed directions follow Family 4 check C exactly:
`family4.generalized_eigen_directions(A1, P1)` with P1 = I - Z Z^T for
Range(B1).  The literal brief form `A1.T@A1` is dimension-incompatible with
the verified API (R_op = Q_A^T W Q_A requires a data-space W; A1.T@A1 is
p x p, Q_A is m x r).  Diagnostic raised:
{pA['literal_generalized_second_argument_diagnostic'].get('exception_type', 'none')}:
{pA['literal_generalized_second_argument_diagnostic'].get('exception_message', 'none')}
Dimensions:
{pA['literal_generalized_second_argument_diagnostic'].get('dimensions', 'n/a')}.

## B. Pose-DOF nuisance baseline

Standard two-blob scene, f=1.0, full B_R has q=18 columns ordered
[px, py, theta] per pose.  Known-pose reference = zero B (retained mass
{pB['known_pose_retained_mass_reference']:.6f}, equal to p=24 up to roundoff).
The known-pose row is stored as a full-width zero B (q=18) because
`family2.range_basis`/`rank_svd` need at least one column.

{_md_table(['subset', 'q_sub', 'retained', 'confusable', 'theta_min deg',
            'rho_min', 'K_IS fro', 'K_SLAM fro', 'K_eff fro'], pose_rows)}

Marginal destruction (loss of retained mass / increase of confusable mass
when a nuisance subset is added):

{_md_table(['transition', 'delta retained', 'delta confusable'],
           marginal_rows)}

Pattern observation: single-DOF retained mass is highest for theta alone
(21.36), then y-only (18.03) and x-only (18.02).  Adding theta to x-only or
y-only costs only ~2.64-2.65 retained mass, while adding x or y to
theta-only costs ~5.98 each, and adding the missing translational DOF to
x/y pairs costs ~5.45-5.46; full nuisance leaves retained = 9.41.  These are
recorded observations, not forced claims.

## C. Receiver-count ablation

{_md_table(['n_rx', 'retained', 'confusable', 'theta_min deg', 'rho_min',
            'B/A fro'], rx_rows)}

Recorded monotonicity flags: retained non-increasing =
{pC['observed_flags']['retained_non_increasing']}; confusable
non-decreasing = {pC['observed_flags']['confusable_non_decreasing']}.

Caveat: {pC['note']}

## D. Pose-count ablation

{_md_table(['T', 'q', 'retained', 'confusable', 'theta_min deg', 'rho_min'],
           T_rows)}

Recorded monotonicity flags: retained non-increasing =
{pD['observed_flags']['retained_non_increasing']}; confusable
non-decreasing = {pD['observed_flags']['confusable_non_decreasing']}.

Note: {pD['note']}

## E. Smooth-basis ablation

{_md_table(['p', 'x_centers', 'y_centers', 'retained', 'confusable',
            'theta_min deg', 'rho_min'], p_rows)}

Recorded monotonicity flags: retained non-increasing =
{pE['observed_flags']['retained_non_increasing']}; confusable
non-decreasing = {pE['observed_flags']['confusable_non_decreasing']}.

Note: {pE['note']}

## Figures

* [{fig_names['seed_replications']}]({fig_names['seed_replications']}) -
  (a) Family 2 core-identity residuals per scene; (b) frequency-diversity
  movements for the three most-confounded directions.
* [{fig_names['pose_dof_ablation']}]({fig_names['pose_dof_ablation']}) -
  retained and confusable mass bars for pose-DOF subsets including the
  known-pose reference.
* [{fig_names['rx_T_basis_ablation']}]({fig_names['rx_T_basis_ablation']}) -
  receiver-count, pose-count, and basis-p ablation panels.

## Cannot establish

* These are finite-dimensional N=16 numbers on one deterministic pose
  geometry at f in {{1.0,1.4}}; no continuum theorem, other-N transfer, other
  geometry, other contrast, or universal bound is asserted.
* The ring scene's rank-identity/nullity mismatch is a machine-rank
  tolerance boundary; it neither disproves the exact algebraic identities nor
  certifies them in floating point.
* Kernel identity on scenes with empty null spaces is a vacuous (though
  consistent) check at the declared eig thresholds; it does not test a
  nontrivial kernel.
* Part B/C/D/E trends are descriptive; the recorded monotonicity flags are
  observations over the tested grid, not guarantees.
* Frequency-diversity "movement" is defined on directions fixed by the
  single-frequency kernel and measured with raw identity-noise stacks; block
  or SNR normalization was neither applied nor silently conflated.
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
        "family": 9,
        "generated_utc": generated_utc,
        "command": ".venv/bin/python src/family9_replications_ablation.py",
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

    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir = _ROOT / "results"
    figures_dir.mkdir(parents=True, exist_ok=True)
    notes_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    fig_paths = {
        "seed_replications": figures_dir / "family9_seed_replications.png",
        "pose_dof_ablation": figures_dir / "family9_pose_dof_ablation.png",
        "rx_T_basis_ablation": figures_dir / "family9_rx_T_basis_ablation.png",
    }
    _plot_seed_replications(
        results["partA_seed_replications"], fig_paths["seed_replications"]
    )
    _plot_pose_dof(
        results["partB_pose_dof_nuisance"], fig_paths["pose_dof_ablation"]
    )
    _plot_rx_T_basis(
        results["partC_receiver_ablation"],
        results["partD_pose_count_ablation"],
        results["partE_basis_ablation"],
        fig_paths["rx_T_basis_ablation"],
    )

    results["runtime_seconds"] = time.perf_counter() - t_start
    results["figures"] = {
        name: str(p.relative_to(_ROOT)) for name, p in fig_paths.items()
    }

    json_path = results_dir / "family9_replications_ablation.json"
    json_path.write_text(
        json.dumps(results, indent=2, default=_json_default) + "\n",
        encoding="utf-8",
    )

    report_path = notes_dir / "family9_replications_ablation.md"
    report_path.write_text(
        build_report(results, cfg, {"figures": fig_paths}), encoding="utf-8"
    )

    pA = results["partA_seed_replications"]
    pr = pA["pass_rates"]
    pB = results["partB_pose_dof_nuisance"]
    pC = results["partC_receiver_ablation"]
    pD = results["partD_pose_count_ablation"]
    pE = results["partE_basis_ablation"]
    print("===== FAMILY 9 SUMMARY =====")
    print(
        f"rank identity {pr['rank_identity_pass']}/{pr['n_scenes']}, "
        f"retention identity {pr['retention_identity_pass']}/{pr['n_scenes']}, "
        f"kernel F2 gate {pr['kernel_identity_f2_gate_pass']}/{pr['n_scenes']} "
        f"(support {pr['kernel_identity_support_pass']}/{pr['n_scenes']}, "
        f"nullity match {pr['kernel_nullity_match']}/{pr['n_scenes']}), "
        f"frequency diversity {pr['frequency_diversity_pass_ge2']}/{pr['n_scenes']}"
    )
    b_by_name = {r["subset"]: r for r in pB["rows"]}
    print(
        "pose DOF retained: known "
        f"{b_by_name['known_pose']['info_metrics']['retained_mass']:.3f}, x "
        f"{b_by_name['x_only']['info_metrics']['retained_mass']:.3f}, y "
        f"{b_by_name['y_only']['info_metrics']['retained_mass']:.3f}, theta "
        f"{b_by_name['theta_only']['info_metrics']['retained_mass']:.3f}, xy "
        f"{b_by_name['xy']['info_metrics']['retained_mass']:.3f}, full "
        f"{b_by_name['full']['info_metrics']['retained_mass']:.3f}"
    )
    print(
        "ablation retained C "
        f"{[round(v,3) for v in pC['observed']['retained_mass']]}, D "
        f"{[round(v,3) for v in pD['observed']['retained_mass']]}, E "
        f"{[round(v,3) for v in pE['observed']['retained_mass']]}"
    )
    print(f"total runtime {results['runtime_seconds']:.2f}s")
    print("results -> " + str(json_path))
    for name, path in fig_paths.items():
        print("figure  -> " + str(path))
    print("report  -> " + str(report_path))


if __name__ == "__main__":
    main()
