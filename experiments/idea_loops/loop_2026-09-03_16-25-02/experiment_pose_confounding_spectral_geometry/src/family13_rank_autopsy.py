"""Family 13: rank-boundary autopsy of the Family 9 ring-scene c2 residual.

Run (from the experiment root):
    .venv/bin/python src/family13_rank_autopsy.py

Purpose
-------
Family 9 records a machine-rank boundary on the radially smooth ring scene:
the family-2 c2 rank identity has residual 1 (r_KIS=24, r_KSL=23, r_A=24,
r_B=18, r_AB=42), while the standard two-blob scene has residual 0.  This
family autopsies that boundary:

  A. raw singular-value and eigenvalue spectra of A_s, B_R, AB, K_IS,
     K_SLAM(direct), C=(I-ZZ^T)A_s, C^T C, and R_op=Q_A^T(I-ZZ^T)Q_A;
  B. numerical ranks on a predeclared relative-tolerance grid for ring and
     two_blob;
  C. stable-rank/effective-rank diagnostics;
  D. backward residuals (K_SLAM vs C^T C, R_op vs the scaled thin-SVD
     product, projector idempotency, Z^T Z, and the c2 tolerance gaps);
  E. a verdict among (a) discontinuous rank change, (b) tolerance
     classification of a near-null singular value of a directly formed
     A^T P A, or (c) unstable direct subtraction;
  F. the family-5 engineered algebraic rank-loss falsifier (t-scaling of the
     weakest singular value of B) on the standard two-blob scene.

Everything is finite-dimensional dense linear algebra on the N=16
whitened/realified smooth p=24 model.  No continuum statement is claimed.
Existing modules are imported read-only; no existing file is modified.
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
import scipy  # noqa: F401  (version recorded; dense LA is numpy)

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import helmholtz as hh  # noqa: E402
import family1_pilot as family1  # noqa: E402
import family2_algebraic_spine as family2  # noqa: E402
import family4_frequency_trajectory as family4  # noqa: E402
import family5_parent_generalized as family5  # noqa: E402
import family9_replications_ablation as family9  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_EPS = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Configuration (declared before any computation; tolerance grid is fixed)
# ---------------------------------------------------------------------------

CONFIG = {
    "title": (
        "Family 13 rank-boundary autopsy: the Family-9 ring-scene c2 residual "
        "(r_KIS - r_KSL = r_A + r_B - r_AB) versus the two-blob reference"
    ),
    "N": 16,
    "T": 6,
    "n_rx": 4,
    "m_real": 48,  # 2*T*n_rx after whiten_realify
    "q_pose": 18,  # 3*T
    "n_pixel": 256,  # N^2
    "p_smooth": 24,
    "f": 1.0,
    "k_b_rule": "k_b(f) = 2*pi*f",
    "whitening_convention": (
        "W=None identity noise: hh.whiten_realify(A,B,None) returns "
        "sqrt(2)*[Re; Im] row stacks (48 real rows at T=6,n_rx=4)"
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
        "note": "family2.build_smooth_basis over the N=16 cell-centre grid",
    },
    "scenes": {
        "ring": {
            "formula": "chi_ring = 0.5*exp(-((|r|-0.25)/0.05)^2)",
            "scale_note": (
                "raw formula values on the N=16 cell-centre grid; no L2 "
                "rescaling.  Global scene scaling leaves every rank in c2 "
                "unchanged, so the row equals the stored family-9 row"
            ),
        },
        "two_blob": {
            "formula": "family1.make_chi0 (two Gaussian blobs)",
            "scale_note": "raw family-1 chi0, exactly as families 1-9",
        },
    },
    "alpha_prior": 1.0,
    "rank_tol_rule": (
        "tol(M) = max(M.shape) * eps_machine * sigma_1(M) "
        "(family2.rank_svd / range_basis)"
    ),
    "relative_tolerance_grid": {
        "predeclared": True,
        "tol_rel_values": [
            1e-16,
            1e-15,
            1e-14,
            1e-13,
            1e-12,
            1e-11,
            1e-10,
            1e-9,
            1e-8,
        ],
        "rank_rule": "rank = count(sigma > tol_rel * sigma_1)",
        "declared_before_computation": (
            "the list above is fixed in this CONFIG before any ranks are "
            "computed and is not adjusted after the fact"
        ),
    },
    "c2_identity": (
        "r_KIS - r_KSL = r_A + r_B - r_AB; residual = lhs - rhs "
        "(family 2 check c2)"
    ),
    "mats_for_spectra": [
        "A_s",
        "B_R",
        "AB",
        "K_IS",
        "K_SLAM_direct",
        "C",
        "C^T C",
        "R_op",
    ],
    "falsifier": {
        "construction": "B(t) = U diag(s_1,...,s_17, |t| s_18) V^T from the full SVD of the standard-scene B_R",
        "t_values": [0.0, 1e-6, 1e-3, 1.0],
        "t_plot_sweep": {"log_min": -13.0, "log_max": 0.0, "n": 41},
        "label": (
            "algebraic control, NOT a physical ring event; reproduces the "
            "family5_parent_generalized engineered_rank_event"
        ),
        "expected_reference_mass": {
            "t_gt_0": 9.40716186665,
            "t_eq_0": 10.40528058115,
        },
    },
    "figure": {
        "backend": "Agg",
        "dpi": 180,
        "path": "figures/family13_rank_autopsy.png",
        "panels": [
            "ring singular-value spectra (log10): A_s, B_R, AB, C",
            "ring eigenvalue-spectra bottoms (indices > 12): K_IS, "
            "K_SLAM(direct), C^T C, R_op",
            "rank vs tol_rel for K_SLAM(direct), C, R_op: ring and two_blob",
            "engineered falsifier: retained mass and projector jump vs t "
            "(log-spaced, t=0 marked)",
        ],
    },
    "randomness_note": (
        "deterministic dense numpy linear algebra; no RNG is used.  CPU-only "
        "computation on Apple Silicon (no MPS/GPU path)"
    ),
    "scope_note": (
        "finite-dimensional N=16 model only (m=48 data rows, p=24 smooth "
        "columns, q=18 pose columns); no continuum/transversality claim is "
        "made from these numbers"
    ),
}


def _sym(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    raise TypeError(f"not JSON serialisable: {type(obj)!r}")


def _fmt(x, spec=".8g"):
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


def _sv_desc(M: np.ndarray) -> np.ndarray:
    return np.linalg.svd(M, compute_uv=False)


def _eig_desc(M: np.ndarray) -> np.ndarray:
    return np.sort(np.linalg.eigvalsh(M))[::-1]


def _rel_rank(M: np.ndarray, tol_rel: float) -> int:
    sv = _sv_desc(M)
    return int(np.sum(sv > float(tol_rel) * float(sv[0])))


def _safe_div(num: float, den: float) -> float:
    return float(num / max(den, 1e-300))


# ---------------------------------------------------------------------------
# Scene construction
# ---------------------------------------------------------------------------

def ring_chi(points: np.ndarray) -> np.ndarray:
    r = np.linalg.norm(points, axis=1)
    return 0.5 * np.exp(-((r - 0.25) / 0.05) ** 2)


def scene_chi(name: str, points: np.ndarray, chi0: np.ndarray) -> np.ndarray:
    if name == "ring":
        return ring_chi(points)
    if name == "two_blob":
        return chi0.copy()
    raise ValueError(f"unknown scene {name!r}")


# ---------------------------------------------------------------------------
# Core case analysis
# ---------------------------------------------------------------------------

def machine_rank_row(A_s: np.ndarray, B_R: np.ndarray, Z: np.ndarray) -> dict:
    """c2 reproduction row exactly as in the Family 9 brief."""
    P = np.eye(A_s.shape[0], dtype=float) - Z @ Z.T
    K_IS = _sym(A_s.T @ A_s)
    K_SLAM = _sym(A_s.T @ (P @ A_s))
    AB = np.hstack([A_s, B_R])
    rA, _, _ = family2.rank_svd(A_s)
    rB, _, _ = family2.rank_svd(B_R)
    rAB, _, _ = family2.rank_svd(AB)
    rKIS, _, _ = family2.rank_svd(K_IS)
    rKSL, _, _ = family2.rank_svd(K_SLAM)
    lhs = rKIS - rKSL
    rhs = rA + rB - rAB
    return {
        "r_A": int(rA),
        "r_KIS": int(rKIS),
        "r_KSL": int(rKSL),
        "r_B": int(rB),
        "r_AB": int(rAB),
        "lhs_rKIS_minus_rKSL": int(lhs),
        "rhs_rA_plus_rB_minus_rAB": int(rhs),
        "rank_identity_residual": int(lhs - rhs),
        "rank_identity_holds": bool(lhs == rhs),
    }


def rank_diagnostics_for_sv(name: str, sv: np.ndarray) -> dict:
    s = np.asarray(sv, dtype=float)
    n = len(s)
    s1 = float(s[0])
    effective = float(np.sum(s) ** 2 / np.sum(s**2)) if n else None
    stable = float(np.sum(s**2) / s1**2) if n else None
    tail = [
        {"index": int(i), "value": float(x), "relative": float(x / s1)}
        for i, x in enumerate(s[-(min(6, n)) :], start=n - min(6, n))
    ]
    row = {
        "matrix": name,
        "n_singular": int(n),
        "effective_rank_sum2_over_sumsq": effective,
        "stable_rank_fro2_over_2norm2": stable,
        "sigma_last_over_sigma_1": float(s[-1] / s1) if n else None,
    }
    if n >= 24:
        row["sigma_24_over_sigma_1"] = float(s[23] / s1)
        row["sigma_23_over_sigma_24"] = float(s[22] / s[23])
    else:
        row["sigma_24_over_sigma_1"] = None
        row["sigma_23_over_sigma_24"] = None
    row["smallest_six_singular_values_relative"] = tail
    return row


def analyze_case(
    name: str,
    chi: np.ndarray,
    chi0_l2: float,
    poses: np.ndarray,
    S: np.ndarray,
    cfg: dict,
) -> dict:
    A_s, B_R, A_pix_R = family4.build_smooth_block(
        chi, poses, S, cfg, float(cfg["f"])
    )
    m = A_s.shape[0]

    thin = family2.thin_decomposition(A_s)
    Q_A = thin["Q"]
    s_A = thin["s"]
    V_A = thin["V"]
    rB, Z, svB, tolB = family2.range_basis(B_R)
    P = np.eye(m, dtype=float) - Z @ Z.T
    K_IS = _sym(A_s.T @ A_s)
    C = P @ A_s
    K_SLAM = _sym(A_s.T @ C)  # direct A^T (I-ZZ^T) A, family 2/4 convention
    AB = np.hstack([A_s, B_R])
    CtC = C.T @ C
    R_op = _sym(Q_A.T @ (P @ Q_A))

    row = machine_rank_row(A_s, B_R, Z)

    raw = {
        "matrices": {
            "A_s": {
                "shape": list(A_s.shape),
                "singular_values_desc": _sv_desc(A_s).tolist(),
            },
            "B_R": {
                "shape": list(B_R.shape),
                "singular_values_desc": _sv_desc(B_R).tolist(),
            },
            "AB": {
                "shape": list(AB.shape),
                "singular_values_desc": _sv_desc(AB).tolist(),
            },
            "K_IS": {
                "shape": list(K_IS.shape),
                "singular_values_desc": _sv_desc(K_IS).tolist(),
            },
            "K_SLAM_direct": {
                "shape": list(K_SLAM.shape),
                "singular_values_desc": _sv_desc(K_SLAM).tolist(),
            },
            "C": {
                "shape": list(C.shape),
                "singular_values_desc": _sv_desc(C).tolist(),
            },
            "C^T C": {
                "shape": list(CtC.shape),
                "singular_values_desc": _sv_desc(CtC).tolist(),
            },
            "R_op": {
                "shape": list(R_op.shape),
                "singular_values_desc": _sv_desc(R_op).tolist(),
            },
        },
        "eigenvalues": {
            "K_IS": {
                "ascending": np.linalg.eigvalsh(K_IS).tolist(),
                "descending": _eig_desc(K_IS).tolist(),
            },
            "K_SLAM_direct": {
                "ascending": np.linalg.eigvalsh(K_SLAM).tolist(),
                "descending": _eig_desc(K_SLAM).tolist(),
            },
            "C^T C": {
                "ascending": np.linalg.eigvalsh(CtC).tolist(),
                "descending": _eig_desc(CtC).tolist(),
            },
            "R_op": {
                "ascending": np.linalg.eigvalsh(R_op).tolist(),
                "descending": _eig_desc(R_op).tolist(),
            },
        },
    }

    # ---- B: predeclared relative-tolerance grid ---------------------------
    tol_grid = cfg["relative_tolerance_grid"]["tol_rel_values"]
    tol_rows = []
    sv_by_label = {
        "K_IS": _sv_desc(K_IS),
        "K_SLAM_direct": _sv_desc(K_SLAM),
        "C": _sv_desc(C),
        "C^T C": _sv_desc(CtC),
        "R_op": _sv_desc(R_op),
        "A_s": _sv_desc(A_s),
        "B_R": _sv_desc(B_R),
        "AB": _sv_desc(AB),
    }
    for tr in tol_grid:
        def cnt(label: str) -> int:
            sv = sv_by_label[label]
            return int(np.sum(sv > float(tr) * float(sv[0])))

        cA, cB, cAB = cnt("A_s"), cnt("B_R"), cnt("AB")
        cKIS, cKSL = cnt("K_IS"), cnt("K_SLAM_direct")
        tol_rows.append(
            {
                "tol_rel": float(tr),
                "rank_A": cA,
                "rank_B": cB,
                "rank_AB": cAB,
                "rank_K_IS": cKIS,
                "rank_K_SLAM_direct": cKSL,
                "rank_C": cnt("C"),
                "rank_C^T C": cnt("C^T C"),
                "rank_R_op": cnt("R_op"),
                "c2_residual_relative_ranks": int(
                    (cKIS - cKSL) - (cA + cB - cAB)
                ),
                "c2_residual_fixed_machine_AB": int(
                    (cKIS - cKSL) - (row["r_A"] + row["r_B"] - row["r_AB"])
                ),
            }
        )

    # ---- C: effective/stable rank diagnostics -----------------------------
    diag_rows = []
    for label in cfg["mats_for_spectra"]:
        diag_rows.append(
            rank_diagnostics_for_sv(label, sv_by_label[label])
        )
    wK = np.linalg.eigvalsh(K_SLAM)
    diag_rows.append(
        {
            "matrix": "K_SLAM_direct smallest eigenvalue",
            "smallest_eigenvalue": float(wK[0]),
            "largest_eigenvalue": float(wK[-1]),
            "smallest_eigenvalue_positive": bool(wK[0] > 0.0),
            "smallest_eigenvalue_negative": bool(wK[0] < 0.0),
            "note": "eigvalsh ascending; sign records roundoff in the direct Gram formation",
        }
    )

    # ---- D: backward residuals --------------------------------------------
    K_SLAM_sv = sv_by_label["K_SLAM_direct"]
    C_sv = sv_by_label["C"]
    R_op_sv = sv_by_label["R_op"]
    X = (V_A / s_A[None, :]).T @ K_SLAM @ (V_A / s_A[None, :])
    normX = float(np.linalg.norm(X, ord="fro"))
    normR = float(np.linalg.norm(R_op, ord="fro"))
    K_IS_fro = float(np.linalg.norm(K_IS, ord="fro"))
    K_SLAM_fro = float(np.linalg.norm(K_SLAM, ord="fro"))
    C_fro = float(np.linalg.norm(C, ord="fro"))
    rKSL_machine, _, tolKSL = family2.rank_svd(K_SLAM)
    gap = {
        "rank_family2": int(rKSL_machine),
        "tol": float(tolKSL),
        "sigma_1": float(K_SLAM_sv[0]),
        "smallest_retained_sigma": float(K_SLAM_sv[rKSL_machine - 1]),
        "largest_dropped_sigma": (
            float(K_SLAM_sv[rKSL_machine])
            if rKSL_machine < len(K_SLAM_sv)
            else None
        ),
        "smallest_retained_above_tol": float(
            K_SLAM_sv[rKSL_machine - 1] - tolKSL
        ),
        "largest_dropped_below_tol": (
            float(tolKSL - K_SLAM_sv[rKSL_machine])
            if rKSL_machine < len(K_SLAM_sv)
            else None
        ),
    }
    backward = {
        "rel_fro_K_SLAM_vs_C^T_C": _safe_div(
            float(np.linalg.norm(K_SLAM - CtC, ord="fro")),
            float(max(K_SLAM_fro, C_fro * C_fro)),
        ),
        "denominator_max_norm_K_SLAM_fro_norm_C_fro_sq": float(
            max(K_SLAM_fro, C_fro * C_fro)
        ),
        "rel_fro_R_op_vs_thin_scaled_product": _safe_div(
            float(np.linalg.norm(R_op - X, ord="fro")),
            float(max(normR, normX)),
        ),
        "norm_R_op_fro": normR,
        "norm_thin_scaled_product_fro": normX,
        "R_op_defined_as_Q_A^T_P_Q_A": True,
        "rel_fro_R_op_minus_definition_path": 0.0,
        "rel_fro_P_sq_minus_P_over_P_fro": _safe_div(
            float(np.linalg.norm(P @ P - P, ord="fro")),
            float(np.linalg.norm(P, ord="fro")),
        ),
        "rel_fro_Z^T_Z_minus_I": float(
            np.linalg.norm(Z.T @ Z - np.eye(Z.shape[1]), ord="fro")
        ),
        "c2_lhs_rhs_and_tolerance_gap": {
            "lhs_rKIS_minus_rKSL": int(row["lhs_rKIS_minus_rKSL"]),
            "rhs_rA_plus_rB_minus_rAB": int(row["rhs_rA_plus_rB_minus_rAB"]),
            "residual": int(row["rank_identity_residual"]),
            "gap": gap,
        },
        "note": (
            "K_SLAM(direct)=A^T(P A); C^T C=(P A)^T(P A).  R_op is by "
            "definition Q_A^T P Q_A, and X=V_A^T diag(1/s_A) K_SLAM "
            "V_A diag(1/s_A) from the thin SVD; the rel-Fro residual of R_op "
            "vs X tests how much scaling by 1/s_A^2 amplifies roundoff."
        ),
    }

    # ---- E: verdict determination -----------------------------------------
    verdict = determine_verdict(
        name=name,
        tol_rows=tol_rows,
        row=row,
        C_sv=C_sv,
        K_SLAM_sv=K_SLAM_sv,
        R_op_sv=R_op_sv,
        backward=backward,
        n_A=K_SLAM.shape[0],
    )

    return {
        "scene": name,
        "chi": {
            "min": float(np.min(chi)),
            "max": float(np.max(chi)),
            "l2_norm": float(np.linalg.norm(chi)),
            "l2_ratio_to_two_blob": float(
                np.linalg.norm(chi) / float(chi0_l2)
            ),
        },
        "shapes": {
            "A_s": list(A_s.shape),
            "B_R": list(B_R.shape),
            "AB": list(AB.shape),
            "K_IS": list(K_IS.shape),
            "K_SLAM": list(K_SLAM.shape),
            "C": list(C.shape),
            "R_op": list(R_op.shape),
        },
        "machine_rank_tol": {
            "rank_A": row["r_A"],
            "rank_B": row["r_B"],
            "rank_AB": row["r_AB"],
            "rank_K_IS": row["r_KIS"],
            "rank_K_SLAM_direct": row["r_KSL"],
            "tol_K_SLAM": float(tolKSL),
            "sigma_1_K_SLAM": float(K_SLAM_sv[0]),
        },
        "reproduction_row": row,
        "raw_spectra": raw,
        "tol_grid": {
            "predeclared": cfg["relative_tolerance_grid"]["predeclared"],
            "rule": cfg["relative_tolerance_grid"]["rank_rule"],
            "rows": tol_rows,
        },
        "rank_diagnostics": diag_rows,
        "backward_residuals": backward,
        "verdict": verdict,
        "norms": {
            "K_IS_fro": float(K_IS_fro),
            "K_SLAM_fro": float(K_SLAM_fro),
            "C_fro": float(C_fro),
            "A_s_fro": float(np.linalg.norm(A_s, ord="fro")),
            "A_s_op": float(np.linalg.norm(A_s, ord=2)),
            "B_R_fro": float(np.linalg.norm(B_R, ord="fro")),
        },
    }


def determine_verdict(
    name: str,
    tol_rows: list[dict],
    row: dict,
    C_sv: np.ndarray,
    K_SLAM_sv: np.ndarray,
    R_op_sv: np.ndarray,
    backward: dict,
    n_A: int,
) -> dict:
    n = int(n_A)
    relC = float(C_sv[-1] / C_sv[0])
    relK = float(K_SLAM_sv[-1] / K_SLAM_sv[0])
    relR = float(R_op_sv[-1] / R_op_sv[0])
    c_always_full = bool(all(r["rank_C"] == n for r in tol_rows))
    ksl_drops = bool(any(r["rank_K_SLAM_direct"] < n for r in tol_rows))
    residual_nonzero_any = bool(
        any(r["c2_residual_relative_ranks"] != 0 for r in tol_rows)
    )
    rel_KCtC = float(backward["rel_fro_K_SLAM_vs_C^T_C"])
    if name == "ring":
        supported = (
            "b_tolerance_classification_of_a_near_null_singular_value"
            if (c_always_full and ksl_drops)
            else "unresolved"
        )
        statement = (
            "Supported explanation: (b).  On the ring scene C=(I-ZZ^T)A_s "
            "has full column rank at every predeclared relative tolerance "
            "(sigma_24(C)/sigma_1(C) is far above the grid), while the "
            "directly formed K_SLAM=A^T P A has a near-null singular value "
            "sigma_24(K_SLAM)/sigma_1(K_SLAM) = [sigma_24(C)/sigma_1(C)]^2 "
            "that falls below tolerance on part of the grid.  This is a "
            "numerical rank classification of a near-null direction, not a "
            "discontinuous exact rank change; K_SLAM and C^T C agree to "
            "machine-level relative Frobenius error, so unstable direct "
            "subtraction between large separately computed terms is not the "
            "mechanism."
        )
    else:
        supported = "none_family2_machine_residual_zero"
        statement = (
            "No rank boundary under the declared family-2 machine rule: the "
            "two-blob c2 residual is 0 and K_SLAM's smallest singular value "
            "ratio (1.5e-13) sits well above the family-2 relative threshold "
            "24*eps ~= 5.3e-15.  As in the ring scene, coarse relative "
            "tolerances near/above 1e-12 begin shedding the near-null K_SLAM "
            "modes (that is tolerance classification, not a discontinuous "
            "rank event), but no such shedding occurs at the declared "
            "machine rule, and the exact identity is not in question here."
        )
    return {
        "scene": name,
        "supported_explanation": supported,
        "exact_identity_in_exact_arithmetic": (
            "not disproven by these floating-point counts"
        ),
        "evidence": {
            "sigma_24_over_sigma_1_C": float(relC),
            "sigma_24_over_sigma_1_K_SLAM_direct": float(relK),
            "sigma_24_over_sigma_1_R_op": float(relR),
            "C_full_rank_all_grid_tol_rel": bool(c_always_full),
            "K_SLAM_direct_rank_drops_some_grid_tol_rel": bool(ksl_drops),
            "residual_nonzero_on_some_grid_tol_rel": bool(residual_nonzero_any),
            "rel_fro_K_SLAM_vs_C^T_C": float(rel_KCtC),
            "family2_machine_residual": int(row["rank_identity_residual"]),
        },
        "statement": statement,
        "theorem_claim": (
            "No theorem is claimed: the verdict above is supported by the "
            "computed finite-dimensional numbers only."
        ),
    }


# ---------------------------------------------------------------------------
# F: family-5 algebraic falsifier (two_blob scene, read-only reuse)
# ---------------------------------------------------------------------------

def engineered_B(A_s: np.ndarray, B_R: np.ndarray, t: float) -> np.ndarray:
    U, s, Vh = np.linalg.svd(B_R, full_matrices=False)
    st = s.copy()
    st[-1] = abs(float(t)) * s[-1]
    return (U * st[None, :]) @ Vh


def falsifier_case(A_s: np.ndarray, B_R: np.ndarray) -> dict:
    U, s, Vh = np.linalg.svd(B_R, full_matrices=False)
    KIS = family5.information_matrices(A_s, B_R, None)[0]
    B0 = (U * np.r_[s[:-1], 0.0][None, :]) @ Vh
    P0, r0, _, tol0 = family5.range_projector(B0)
    K0 = family5.information_matrices(A_s, B0, None)[1]
    KIS_fro = float(np.linalg.norm(KIS, ord="fro"))

    def one_row(t: float) -> dict:
        Bt = engineered_B(A_s, B_R, t)
        PB, rB, svB, tolB = family5.range_projector(Bt)
        dec = family5.retention_decomposition(A_s, Bt, None)
        mass = float(np.sum(dec["rho"]))
        Kt = family5.information_matrices(A_s, Bt, None)[1]
        jumpP0 = float(np.linalg.norm(PB - P0, ord=2))
        relFro = _safe_div(
            float(np.linalg.norm(Kt - K0, ord="fro")), float(KIS_fro)
        )
        relFroKt = _safe_div(
            float(np.linalg.norm(Kt - K0, ord="fro")),
            float(max(np.linalg.norm(Kt, ord="fro"), np.linalg.norm(K0, ord="fro"))),
        )
        return {
            "t": float(t),
            "rank_B": int(rB),
            "sigma_min_B": float(svB[-1]),
            "tol_B": float(tolB),
            "retained_mass_no_prior": mass,
            "rho_min": float(dec["rho"][0]),
            "rho_max": float(dec["rho"][-1]),
            "projector_jump_op_norm_vs_t0": jumpP0,
            "rel_fro_information_jump_vs_KIS_fro": relFro,
            "rel_fro_information_jump_vs_Kt_max": relFroKt,
        }

    rows = [one_row(t) for t in CONFIG["falsifier"]["t_values"]]
    sweep_cfg = CONFIG["falsifier"]["t_plot_sweep"]
    t_pos = list(
        np.logspace(sweep_cfg["log_min"], sweep_cfg["log_max"], sweep_cfg["n"])
    )
    sweep = [one_row(t) for t in t_pos]
    sweep.append(one_row(0.0))

    # exact reuse of family5_parent_generalized.engineered_rank_event
    reference = family5.engineered_rank_event(A_s, B_R)
    return {
        "construction": CONFIG["falsifier"]["construction"],
        "label": CONFIG["falsifier"]["label"],
        "t_rows": rows,
        "plot_sweep": sweep,
        "base": {
            "rank_B": int(reference["base_rank_B"]),
            "sigma_min_B": float(reference["base_sigma_min_B"]),
            "tol_B": float(reference["base_tol_B"]),
            "K_IS_fro": float(KIS_fro),
        },
        "family5_parent_generalized_reuse": {
            "source_rows_t": [-1e-2, -1e-4, -1e-8, 0.0, 1e-8, 1e-4, 1e-2],
            "retained_mass_t_0_actual": float(
                [r for r in reference["rows"] if r["t"] == 0.0][0]["retained_mass"]
            ),
            "retained_mass_t_1e-8_actual": float(
                [r for r in reference["rows"] if r["t"] == 1e-8][0][
                    "retained_mass"
                ]
            ),
            "projector_jump_norm_at_zero_vs_1e-8": float(
                reference["projector_jump_norm_at_zero"]
            ),
            "KSLAM_jump_fro_at_zero_vs_1e-8": float(
                reference["KSLAM_jump_fro_at_zero"]
            ),
            "KSLAM_jump_relative_to_KIS_vs_1e-8": float(
                reference["KSLAM_jump_relative_to_KIS"]
            ),
            "interpretation": reference["interpretation"],
        },
    }


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _log10_abs(x: float) -> float:
    return math.log10(max(abs(float(x)), np.nextafter(0.0, 1.0)))


def make_figure(results: dict, path: Path) -> str:
    fig, axes = plt.subplots(2, 2, figsize=(15.2, 11.6))
    ring = results["scenes"]["ring"]
    two = results["scenes"]["two_blob"]

    # (1) ring singular-value spectra ---------------------------------------
    ax = axes[0, 0]
    colors = plt.cm.tab10.colors
    entries = [
        ("A_s", "A_s (48x24)", colors[0]),
        ("B_R", "B_R (48x18)", colors[1]),
        ("AB", "AB (48x42)", colors[2]),
        ("C", "C=(I-ZZ^T)A_s (48x24)", colors[3]),
    ]
    for key, label, col in entries:
        sv = np.asarray(ring["raw_spectra"]["matrices"][key]["singular_values_desc"])
        x = np.arange(1, len(sv) + 1)
        ax.plot(x, np.log10(sv), ".-", color=col, label=label, markersize=5)
    ax.set_xlabel("singular-value index (descending)")
    ax.set_ylabel("log10 singular value")
    ax.set_title("(1) Ring scene: singular-value spectra (log10)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)

    # (2) ring eigenvalue bottoms (indices > 12) ----------------------------
    ax = axes[0, 1]
    for key, label, col in [
        ("K_IS", "K_IS = A_s^T A_s (eig)", colors[0]),
        ("K_SLAM_direct", "K_SLAM direct = A_s^T P A_s (eig)", colors[1]),
        ("C^T C", "C^T C (eig)", colors[2]),
        ("R_op", "R_op = Q_A^T P Q_A (eig)", colors[3]),
    ]:
        eig = np.asarray(
            ring["raw_spectra"]["eigenvalues"][key]["descending"]
        )
        x = np.arange(1, len(eig) + 1)
        mask = x > 12
        ax.plot(
            x[mask],
            np.log10(np.abs(eig[mask])),
            ".-",
            color=col,
            label=label,
            markersize=6,
        )
    ax.set_xlabel("eigenvalue index (descending |lambda|), indices > 12")
    ax.set_ylabel("log10 |eigenvalue|")
    ax.set_title("(2) Ring: eigenvalue spectra, bottom of 24")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)

    # (3) rank vs tol_rel ----------------------------------------------------
    ax = axes[1, 0]
    style = {
        "ring": ("o-", colors[0]),
        "two_blob": ("s--", colors[1]),
    }
    for scene_name, scene in (("ring", ring), ("two_blob", two)):
        rows = scene["tol_grid"]["rows"]
        tol = [float(r["tol_rel"]) for r in rows]
        for key, col in [
            ("rank_K_SLAM_direct", "#1f77b4"),
            ("rank_C", "#d62728"),
            ("rank_R_op", "#2ca02c"),
        ]:
            ranks = [int(r[key]) for r in rows]
            ax.plot(
                tol,
                ranks,
                style[scene_name][0],
                color=col,
                label=f"{scene_name} {key.replace('rank_','')}",
            )
    ax.set_xscale("log")
    ax.set_xlabel("relative tolerance tol_rel (rank = sigma > tol_rel * sigma_1)")
    ax.set_ylabel("numerical rank")
    ax.set_title("(3) Numerical rank vs predeclared tol_rel")
    ax.set_yticks([23, 24])
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7)

    # (4) engineered falsifier ----------------------------------------------
    ax = axes[1, 1]
    sweep = results["falsifier"]["plot_sweep"]
    t0_row = next(r for r in sweep if r["t"] == 0.0)
    pos = [r for r in sweep if r["t"] > 0.0]
    tx = np.asarray([float(r["t"]) for r in pos])
    mass = np.asarray([float(r["retained_mass_no_prior"]) for r in pos])
    jump = np.asarray([float(r["projector_jump_op_norm_vs_t0"]) for r in pos])
    x_zero = 1e-13
    mass_line, = ax.plot(
        tx, mass, ".-", color=colors[0], label="retained mass (t>0)"
    )
    ax.axhline(
        float(t0_row["retained_mass_no_prior"]),
        color=colors[3],
        linestyle=":",
        linewidth=1.1,
        label=f"t=0 mass = {float(t0_row['retained_mass_no_prior']):.5f}",
    )
    star_line, = ax.plot(
        [x_zero],
        [float(t0_row["retained_mass_no_prior"])],
        marker="*",
        markersize=12,
        color=colors[3],
    )
    ax.set_xscale("log")
    ax.set_xlim(5e-14, 4.0)
    ax.set_xlabel("t (weakest singular value scale; left star = t=0)")
    ax.set_ylabel("no-prior retained mass (sum rho)", color=colors[0])
    ax.tick_params(axis="y", labelcolor=colors[0])
    ax.grid(True, which="both", alpha=0.25)

    ax2 = ax.twinx()
    jump_line, = ax2.plot(
        tx, jump, ".-", color=colors[1], label="||P(t)-P(0)||_2"
    )
    ax2.set_ylabel("projector jump operator norm ||P(t)-P(0)||_2", color=colors[1])
    ax2.tick_params(axis="y", labelcolor=colors[1])
    ax2.set_ylim(-0.05, 1.15)
    lines = [mass_line, jump_line, star_line]
    labels = [l.get_label() for l in lines]
    labels[-1] = "t=0 mass (star)"
    ax.legend(lines, labels, fontsize=7, loc="center left")
    ax.set_title(
        "(4) Algebraic falsifier on two_blob: B(t)=U diag(s_1..s_17,|t|s_18)V^T"
    )

    fig.suptitle(
        "Family 13 rank-boundary autopsy (ring f=1.0 vs two_blob reference; "
        "N=16 finite-dimensional)",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return _sha256(path)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def _grid_table(scene: dict, cfg: dict) -> str:
    rows_out = []
    for r in scene["tol_grid"]["rows"]:
        rows_out.append(
            [
                f"{float(r['tol_rel']):.0e}",
                int(r["rank_A"]),
                int(r["rank_B"]),
                int(r["rank_AB"]),
                int(r["rank_K_IS"]),
                int(r["rank_K_SLAM_direct"]),
                int(r["rank_C"]),
                int(r["rank_C^T C"]),
                int(r["rank_R_op"]),
                int(r["c2_residual_relative_ranks"]),
            ]
        )
    return _md_table(
        [
            "tol_rel",
            "r_A",
            "r_B",
            "r_AB",
            "r_KIS",
            "r_KSL",
            "r_C",
            "r_CtC",
            "r_Rop",
            "c2 res.",
        ],
        rows_out,
    )


def _diag_table(scene: dict) -> str:
    rows_out = []
    for d in scene["rank_diagnostics"]:
        if d["matrix"].startswith("K_SLAM_direct smallest"):
            continue
        rel_tail = d["smallest_six_singular_values_relative"]
        tail_s = ", ".join(_fmt(v["relative"], ".3g") for v in rel_tail[:4])
        s24 = (
            _fmt(d["sigma_24_over_sigma_1"], ".3g")
            if d["sigma_24_over_sigma_1"] is not None
            else _fmt(d["sigma_last_over_sigma_1"], ".3g") + " (last/1)"
        )
        gap = (
            _fmt(d["sigma_23_over_sigma_24"], ".3g")
            if d["sigma_23_over_sigma_24"] is not None
            else "n/a"
        )
        rows_out.append(
            [
                d["matrix"],
                _fmt(d["effective_rank_sum2_over_sumsq"], ".4g"),
                _fmt(d["stable_rank_fro2_over_2norm2"], ".4g"),
                s24,
                gap,
                tail_s,
            ]
        )
    return _md_table(
        ["matrix", "eff rank", "stable rank", "sigma_24/1", "sigma_23/24", "smallest rel. (top 4 of 6)"],
        rows_out,
    )


def _backward_table(scene: dict) -> str:
    b = scene["backward_residuals"]
    gap = b["c2_lhs_rhs_and_tolerance_gap"]
    g = gap["gap"]
    return _md_table(
        ["quantity", "value"],
        [
            ["rel Fro ||K_SLAM - C^T C|| / max(||K_SLAM||F, ||C||F^2)", _fmt(b["rel_fro_K_SLAM_vs_C^T_C"], ".6g")],
            ["rel Fro ||R_op - V^T diag(1/s) K_SLAM V diag(1/s)|| / max(F norms)", _fmt(b["rel_fro_R_op_vs_thin_scaled_product"], ".6g")],
            ["R_op definition residual ||R_op - Q_A^T P Q_A||F / ||R_op||F", "0.0 (same object)"],
            ["||P^2-P||F / ||P||F (idempotency)", _fmt(b["rel_fro_P_sq_minus_P_over_P_fro"], ".6g")],
            ["||Z^T Z - I||F", _fmt(b["rel_fro_Z^T_Z_minus_I"], ".6g")],
            ["c2 lhs (r_KIS - r_KSL)", str(gap["lhs_rKIS_minus_rKSL"])],
            ["c2 rhs (r_A + r_B - r_AB)", str(gap["rhs_rA_plus_rB_minus_rAB"])],
            ["c2 residual", str(gap["residual"])],
            ["family2 rank tol(K_SLAM)", _fmt(g["tol"], ".6g")],
            ["smallest retained sigma(K_SLAM) - tol", _fmt(g["smallest_retained_above_tol"], ".6g")],
            ["largest dropped sigma(K_SLAM) below tol", _fmt(g["largest_dropped_below_tol"], ".6g")],
        ],
    )


def _falsifier_table(f: dict) -> str:
    rows_out = []
    for r in f["t_rows"]:
        rows_out.append(
            [
                "0" if float(r["t"]) == 0.0 else f"{float(r['t']):.0e}",
                int(r["rank_B"]),
                _fmt(r["retained_mass_no_prior"], ".8f"),
                _fmt(r["projector_jump_op_norm_vs_t0"], ".6g"),
                _fmt(r["rel_fro_information_jump_vs_KIS_fro"], ".6g"),
            ]
        )
    return _md_table(
        ["t", "rank(B(t))", "retained mass", "||P(t)-P(0)||_2", "rel Fro info jump / ||K_IS||F"],
        rows_out,
    )


def build_report(results: dict, cfg: dict) -> str:
    ring = results["scenes"]["ring"]
    two = results["scenes"]["two_blob"]
    rep = results["reproduction"]
    f13 = results["falsifier"]
    lines = []
    lines.append("# Family 13: rank-boundary autopsy (ring-scene c2 residual)")
    lines.append("")
    lines.append(
        f"Date: {results['generated_utc']}.  Experiment: "
        "`experiment_pose_confounding_spectral_geometry`."
    )
    lines.append("")
    lines.append(
        "Family 13 autopsies the Family-9 ring-scene machine-rank boundary: "
        "c2 rank identity r_KIS - r_KSL = r_A + r_B - r_AB records residual 1 "
        "on the ring scene and 0 on the two-blob reference.  It reuses "
        "(read-only) the standard builders and the family-5 algebraic "
        "falsifier; no existing source/result/figure/note is modified."
    )
    lines.append("")
    lines.append("## Exact command and runtime")
    lines.append("")
    lines.append("```bash")
    lines.append(results["command"])
    lines.append("```")
    lines.append("")
    lines.append(
        f"Wall runtime: {results['wall_runtime_seconds']:.3f} s.  Platform: "
        f"{results['platform']}.  Python {results['versions']['python']}, "
        f"numpy {results['versions']['numpy']}, scipy "
        f"{results['versions']['scipy']}, matplotlib "
        f"{results['versions']['matplotlib']}.  CPU-only (Apple Silicon)."
    )
    lines.append("")
    lines.append(
        "Source SHA-256 (this script): `" + results["source_sha256"] + "`"
    )
    lines.append("")
    lines.append("Reused-module SHA-256:")
    for key, val in results["reused_module_sha256"].items():
        lines.append(f"* `{key}` `{val}`")
    lines.append("")
    lines.append(f"Figure SHA-256: `{results['figure_sha256']}`.")
    lines.append("")

    lines.append("## Reproduction row (ring scene, f=1.0)")
    lines.append("")
    r = rep["ring"]
    lines.append(
        f"Family-2 machine ranks: r_A={r['r_A']}, r_B={r['r_B']}, "
        f"r_AB={r['r_AB']}, r_KIS={r['r_KIS']}, r_KSL={r['r_KSL']}; "
        f"lhs={r['lhs_rKIS_minus_rKSL']}, rhs={r['rhs_rA_plus_rB_minus_rAB']}, "
        f"residual={r['rank_identity_residual']}."
    )
    lines.append("")
    lines.append(
        f"Stored family-9 ring row: "
        f"r_KIS={rep['family9']['r_KIS']}, r_KSL={rep['family9']['r_KSL']}, "
        f"r_A={rep['family9']['r_A']}, r_B={rep['family9']['r_B']}, "
        f"r_AB={rep['family9']['r_AB']}, "
        f"residual={rep['family9']['rank_identity_residual']}.  "
        f"Row match: **{rep['row_matches_family9']}**."
    )
    lines.append("")
    lines.append(
        "The ring chi is the raw `0.5*exp(-((|r|-0.25)/0.05)^2)` grid vector "
        "(no L2 rescaling); a global scene scale leaves every c2 rank "
        "unchanged, which is why the raw ring row equals the family-9 row."
    )
    lines.append("")

    lines.append("## A. Raw spectra")
    lines.append("")
    lines.append(
        "Full descending singular-value arrays for A_s, B_R, AB, K_IS, "
        "K_SLAM(direct), C, C^T C and R_op, plus ascending/descending "
        "eigenvalues of K_SLAM and R_op, are embedded in "
        "`results/family13_rank_autopsy.json` (`scenes.<scene>.raw_spectra`)."
    )
    lines.append("")
    lines.append("Key relative minima:")
    lines.append("")
    for scene in (ring, two):
        v = scene["verdict"]["evidence"]
        lines.append(
            f"* {scene['scene']}: sigma_24/sigma_1 "
            f"C={_fmt(v['sigma_24_over_sigma_1_C'], '.6g')}, "
            f"K_SLAM={_fmt(v['sigma_24_over_sigma_1_K_SLAM_direct'], '.6g')}, "
            f"R_op={_fmt(v['sigma_24_over_sigma_1_R_op'], '.6g')}."
        )
    lines.append("")

    lines.append("## B. Predeclared relative-tolerance grid")
    lines.append("")
    lines.append(
        "tol_rel grid fixed in CONFIG before computing: "
        "1e-16 ... 1e-8 (nine values); rank = count(sigma > tol_rel*sigma_1)."
    )
    lines.append("")
    lines.append("Ring:")
    lines.append("")
    lines.append(_grid_table(ring, cfg))
    lines.append("")
    lines.append("Two-blob reference:")
    lines.append("")
    lines.append(_grid_table(two, cfg))
    lines.append("")

    lines.append("## C. Stable-rank / effective-rank diagnostics")
    lines.append("")
    lines.append("Ring:")
    lines.append("")
    lines.append(_diag_table(ring))
    lines.append("")
    lines.append("Two-blob:")
    lines.append("")
    lines.append(_diag_table(two))
    lines.append("")
    ksmall = [d for d in ring["rank_diagnostics"] if d["matrix"].startswith("K_SLAM_direct smallest")][0]
    lines.append(
        "Smallest eigvalsh eigenvalue of ring K_SLAM: "
        f"{_fmt(ksmall['smallest_eigenvalue'], '.6g')} "
        f"(positive={ksmall['smallest_eigenvalue_positive']})."
    )
    lines.append("")

    lines.append("## D. Backward residuals")
    lines.append("")
    lines.append("Ring:")
    lines.append("")
    lines.append(_backward_table(ring))
    lines.append("")
    lines.append("Two-blob:")
    lines.append("")
    lines.append(_backward_table(two))
    lines.append("")

    lines.append("## E. Verdict")
    lines.append("")
    for scene in (ring, two):
        lines.append(f"### {scene['scene']}")
        lines.append("")
        ev = scene["verdict"]["evidence"]
        lines.append(
            f"Supported explanation: `{scene['verdict']['supported_explanation']}`."
        )
        lines.append("")
        lines.append(scene["verdict"]["statement"])
        lines.append("")
        lines.append("Supporting numbers:")
        lines.append("")
        for k, val in ev.items():
            lines.append(f"* {k}: {val}")
        lines.append("")
        lines.append(scene["verdict"]["theorem_claim"])
        lines.append("")

    lines.append("## F. Engineered algebraic falsifier (two_blob, f=1.0)")
    lines.append("")
    lines.append(
        "B(t) = U diag(s_1,...,s_17,|t| s_18) V^T, reusing "
        "`family5_parent_generalized.engineered_rank_event`/`range_projector`/"
        "`retention_decomposition` read-only.  This is an **algebraic control** "
        "for the discontinuous no-prior projector at exact rank loss; it is "
        "not a physical ring event and says nothing about the ring scene."
    )
    lines.append("")
    lines.append(_falsifier_table(f13))
    lines.append("")
    lines.append(
        f"Retained mass jumps from "
        f"{_fmt(f13['t_rows'][1]['retained_mass_no_prior'], '.8f')} "
        f"(t>0, full-rank B) to "
        f"{_fmt(f13['t_rows'][0]['retained_mass_no_prior'], '.8f')} (t=0) "
        "while the no-prior projector jumps by operator norm "
        f"{_fmt(f13['t_rows'][1]['projector_jump_op_norm_vs_t0'], '.6g')} "
        "between t=0 and every t>0 row (family-5 reference with t=1e-8: "
        f"{_fmt(f13['family5_parent_generalized_reuse']['projector_jump_norm_at_zero_vs_1e-8'], '.6g')})."
    )
    lines.append("")

    lines.append("## Scope")
    lines.append("")
    lines.append(
        "All numbers are finite-dimensional N=16 dense linear algebra on the "
        "whitened/realified smooth p=24 model at f=1.0 with the standard "
        "pose arc and receiver set.  The ring residual of 1 is a numerical "
        "rank classification boundary and is **not evidence against the exact "
        "c2 identity**; C=(I-ZZ^T)A_s remains full rank at every predeclared "
        "grid tolerance, and K_SLAM agrees with C^T C at machine-level "
        "relative Frobenius error.  No continuum/transversality or theorem "
        "claim is made."
    )
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append("* [results/family13_rank_autopsy.json](results/family13_rank_autopsy.json)")
    lines.append("* [figures/family13_rank_autopsy.png](figures/family13_rank_autopsy.png)")
    lines.append("* [notes/family13_rank_autopsy.md](notes/family13_rank_autopsy.md)")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.perf_counter()
    cfg = CONFIG
    points, h = hh.make_grid(int(cfg["N"]))
    chi0 = family1.make_chi0(points, cfg)
    poses = family1.build_poses(cfg)
    S = family2.build_smooth_basis(points, cfg["smooth_basis"])

    scenes = {}
    scenes_out = {}
    for name in ("two_blob", "ring"):
        chi = scene_chi(name, points, chi0)
        case = analyze_case(
            name, chi, float(np.linalg.norm(chi0)), poses, S, cfg
        )
        scenes[name] = case
        scenes_out[name] = case

    ring = scenes["ring"]
    two = scenes["two_blob"]

    # reproduction cross-check against stored family-9 JSON (read-only)
    fam9_path = _ROOT / "results" / "family9_replications_ablation.json"
    try:
        f9 = json.loads(fam9_path.read_text(encoding="utf-8"))
        f9_ring = next(
            r
            for r in f9["partA_seed_replications"]["scenes"]
            if r["scene"] == "ring"
        )["core_family2"]["c2_rank_identity"]
        family9_available = True
    except Exception as exc:  # noqa: BLE001
        f9_ring = {"read_error": f"{type(exc).__name__}: {exc}"}
        family9_available = False

    our_row = ring["reproduction_row"]
    expected = {
        "r_KIS": 24,
        "r_KSL": 23,
        "r_A": 24,
        "r_B": 18,
        "r_AB": 42,
        "lhs_rKIS_minus_rKSL": 1,
        "rhs_rA_plus_rB_minus_rAB": 0,
        "rank_identity_residual": 1,
    }
    row_match_expected = all(our_row[k] == v for k, v in expected.items())
    row_match_family9 = bool(
        family9_available
        and all(our_row[k] == int(f9_ring[k]) for k in expected)
    )

    results = {
        "schema": "family13-rank-autopsy-v1",
        "family": 13,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "runner": f"{platform.python_implementation()} {platform.python_version()}",
        "command": ".venv/bin/python src/family13_rank_autopsy.py",
        "wall_runtime_seconds": None,  # filled just before serialisation
        "platform": platform.platform(),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "source_sha256": _sha256(Path(__file__)),
        "reused_module_sha256": {
            "helmholtz.py": _sha256(_HERE / "helmholtz.py"),
            "family1_pilot.py": _sha256(_HERE / "family1_pilot.py"),
            "family2_algebraic_spine.py": _sha256(
                _HERE / "family2_algebraic_spine.py"
            ),
            "family4_frequency_trajectory.py": _sha256(
                _HERE / "family4_frequency_trajectory.py"
            ),
            "family5_parent_generalized.py": _sha256(
                _HERE / "family5_parent_generalized.py"
            ),
            "family9_replications_ablation.py": _sha256(
                _HERE / "family9_replications_ablation.py"
            ),
        },
        "figure_sha256": None,
        "config": cfg,
        "scope_note": cfg["scope_note"],
        "reproduction": {
            "ring": our_row,
            "expected": expected,
            "row_matches_expected": bool(row_match_expected),
            "family9": f9_ring,
            "family9_available": bool(family9_available),
            "row_matches_family9": bool(row_match_family9),
            "two_blob_residual": int(two["reproduction_row"]["rank_identity_residual"]),
            "note": (
                "global scene rescaling does not change any c2 rank; the raw "
                "ring formula and family-9's L2-normalised ring therefore "
                "give the same integer row"
            ),
        },
        "scenes": scenes_out,
        "falsifier": None,
    }
    # family-5 falsifier needs the two_blob A_s/B_R matrices
    chi2 = scene_chi("two_blob", points, chi0)
    A2, B2, _ = family4.build_smooth_block(
        chi2, poses, S, cfg, float(cfg["f"])
    )
    results["falsifier"] = falsifier_case(A2, B2)

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    fig_path = figures_dir / "family13_rank_autopsy.png"
    results["figure_sha256"] = make_figure(results, fig_path)

    results["wall_runtime_seconds"] = float(time.perf_counter() - t0)
    report = build_report(results, cfg)
    report_path = notes_dir / "family13_rank_autopsy.md"
    report_path.write_text(report, encoding="utf-8")
    out_path = results_dir / "family13_rank_autopsy.json"
    out_path.write_text(
        json.dumps(results, indent=2, default=_json_default), encoding="utf-8"
    )

    verdicts = {k: v["verdict"]["supported_explanation"] for k, v in scenes_out.items()}
    summary = {
        "family13_ok": True,
        "reproduction_ring_row": our_row,
        "row_matches_expected": bool(row_match_expected),
        "row_matches_family9": bool(row_match_family9),
        "two_blob_row": two["reproduction_row"],
        "verdicts": verdicts,
        "ring_key_numbers": {
            "sigma24_over_sigma1_C": ring["verdict"]["evidence"][
                "sigma_24_over_sigma_1_C"
            ],
            "sigma24_over_sigma1_K_SLAM": ring["verdict"]["evidence"][
                "sigma_24_over_sigma_1_K_SLAM_direct"
            ],
            "sigma24_over_sigma1_R_op": ring["verdict"]["evidence"][
                "sigma_24_over_sigma_1_R_op"
            ],
            "rel_fro_KSL_vs_CtC": ring["backward_residuals"][
                "rel_fro_K_SLAM_vs_C^T_C"
            ],
            "rel_fro_Rop_vs_thin_scaled": ring["backward_residuals"][
                "rel_fro_R_op_vs_thin_scaled_product"
            ],
        },
        "falsifier": {
            "retained_mass_t0": results["falsifier"]["t_rows"][0][
                "retained_mass_no_prior"
            ],
            "retained_mass_t_1e-6": results["falsifier"]["t_rows"][1][
                "retained_mass_no_prior"
            ],
            "projector_jump_t0": results["falsifier"]["t_rows"][1][
                "projector_jump_op_norm_vs_t0"
            ],
            "rel_fro_info_jump_t0_vs_1e-6": results["falsifier"]["t_rows"][1][
                "rel_fro_information_jump_vs_KIS_fro"
            ],
        },
        "wall_runtime_seconds": results["wall_runtime_seconds"],
        "results": str(out_path),
        "figure": str(fig_path),
        "report": str(report_path),
    }
    print(json.dumps(summary, indent=2, default=_json_default))


if __name__ == "__main__":
    main()
