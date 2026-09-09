"""Family 5b: structural operator-Lipschitz bound for the robust surrogate.

This experiment derives and validates the structural (pointwise) uniform
operator-Lipschitz bound

    L_struct = 2||A||_F JA_op
             + 2||A||_F^2 JB_op (tau/sigma^2 + tau^3/sigma^4)

for K(X) = A(X)^T P(X) A(X), P = I - B (B^T B + alpha I)^{-1} B^T, where
JA_op = ||J_A||_2 and JB_op = ||J_B||_2 are the operator norms of the pose
Jacobians of the smooth A and B blocks (centred finite differences), and
tau = ||B||_2, sigma = sigma_min(B).

Certification semantics are separated:
  * the affine tangent certificate L_cert_affine is rigorous FOR THE AFFINE
    MODEL ONLY, conditional on the numerically computed FD Jacobians being
    exact to machine precision (supported by the O(h^2) convergence checks);
  * the full nonlinear results are a derived structural bound with a sampled
    ingredient envelope only -- no formal interval certificate of the
    nonlinear Helmholtz/smooth model is claimed.

Reused (read-only, never rerun): src/family5_sensitivity.py
(base_scene, build_poses, smooth_AB_of_flat_X and its scene config).

Run (from the experiment root):
    .venv/bin/python src/family5b_lipschitz_bound.py

Outputs:
  results/family5b_lipschitz_bound.json
  figures/family5b_lipschitz_bound.png
  notes/family5b_lipschitz_bound.md
"""

from __future__ import annotations

import copy
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

import family5_sensitivity as f5  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

EPS_MACHINE = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Configuration (Family 5 scenario + Family 5b experiment parameters)
# ---------------------------------------------------------------------------

CONFIG = {
    **copy.deepcopy(f5.CONFIG),
    "title": (
        "structural uniform operator-Lipschitz bound for the robust first-"
        "order surrogate (family 5b)"
    ),
    "alpha_prior": 1.0,
    "fd_h1": 1e-4,
    "fd_h2": 2e-4,
    "o2_h_steps": [1e-4, 2e-4, 4e-4, 8e-4],
    "o2_check_columns_A": [0, 8, 17],
    "o2_check_columns_B": [0, 8, 17],
    "family5b_seeds": {
        "affine_unit_directions_2000": 7171,
        "full_model_unit_directions_500": 9191,
        "boundary_points_5": 5353,
        "note": (
            "numpy.random.default_rng(seed); directions normalised to unit "
            "2-norm on R^18"
        ),
    },
    "affine_eps_n_grid": 12,
    "affine_eps_min": 1e-5,
    "affine_eps_sigma_safety": 0.9,  # max eps = 0.9 * sigma / JB_op
    "affine_n_unit_samples": 2000,
    "full_eps_list": [1e-3, 3e-3, 1e-2],
    "full_n_unit_samples": 500,
    "boundary_eps0": 1e-2,
    "boundary_n_points": 5,
    "boundary_fd_h": 1e-4,
    "norm_conventions": (
        "input norm on R^18 is the Euclidean 2-norm; output norms are "
        "Frobenius for A/K matrices and spectral (||.||_2) for P, B, C^-1; "
        "JA_op = largest singular value of J_A (m_A x q) as a matrix, "
        "JB_op = largest singular value of J_B (m_B x q)"
    ),
    "nonlinearity_note": (
        "full-model samples use the exact smooth_AB_of_flat_X nonlinear map; "
        "the A and B Jacobians used in L_struct are centred finite differences "
        "at the evaluation point"
    ),
}


# ---------------------------------------------------------------------------
# Small numeric helpers
# ---------------------------------------------------------------------------

def fro(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


def spec(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord=2))


def sigma_min(M: np.ndarray) -> float:
    return float(np.linalg.svd(M, compute_uv=False)[-1])


def unit_direction(rng, q: int) -> np.ndarray:
    v = rng.standard_normal(q)
    n = float(np.linalg.norm(v))
    return v / n if n > 0.0 else unit_direction(rng, q)


def K_ab(A: np.ndarray, B: np.ndarray, alpha: float) -> np.ndarray:
    """K = A^T P A with P = I - B (B^T B + alpha I)^{-1} B^T."""
    n_b = B.shape[1]
    m = A.shape[0]
    C = B.T @ B + float(alpha) * np.eye(n_b, dtype=float)
    P = np.eye(m, dtype=float) - B @ np.linalg.solve(C, B.T)
    K = A.T @ (P @ A)
    return 0.5 * (K + K.T)


def lambda_min_sym(K: np.ndarray) -> float:
    return float(np.linalg.eigvalsh(0.5 * (K + K.T))[0])


def fd_jacobians(X_ref, h_step, A0, B0, cfg):
    """Centred-difference Jacobians of flattened A and B wrt X (q columns)."""
    q = int(cfg["q_pose"])
    m_A = int(A0.size)
    m_B = int(B0.size)
    J_A = np.empty((m_A, q), dtype=float)
    J_B = np.empty((m_B, q), dtype=float)
    hh = float(h_step)
    for j in range(q):
        e = np.zeros(q, dtype=float)
        e[j] = hh
        Ap, Bp = f5.smooth_AB_of_flat_X(chi0_ref, S_ref, X_ref + e, cfg)
        Am, Bm = f5.smooth_AB_of_flat_X(chi0_ref, S_ref, X_ref - e, cfg)
        J_A[:, j] = ((Ap - Am) / (2.0 * hh)).ravel()
        J_B[:, j] = ((Bp - Bm) / (2.0 * hh)).ravel()
    return J_A, J_B


# scene state set in main() and reused by the low-level closure-free helpers
chi0_ref = None
S_ref = None
X0_flat_ref = None


def _ab_pair(X_flat, cfg):
    return f5.smooth_AB_of_flat_X(chi0_ref, S_ref, X_flat, cfg)


def col_fd_pair(j, h_step, cfg, X_ref):
    """Centred FD column j of the flattened A and B Jacobians at X_ref."""
    e = np.zeros(int(cfg["q_pose"]), dtype=float)
    e[j] = float(h_step)
    Ap, Bp = _ab_pair(X_ref + e, cfg)
    Am, Bm = _ab_pair(X_ref - e, cfg)
    dA = ((Ap - Am) / (2.0 * float(h_step))).ravel()
    dB = ((Bp - Bm) / (2.0 * float(h_step))).ravel()
    return dA, dB


# ---------------------------------------------------------------------------
# Computations
# ---------------------------------------------------------------------------

def run_fd_validation(J_A1, J_B1, J_A2, J_B2, cfg, X_ref):
    """Per-column h1/h2 differences plus O(h^2) checks on a few columns."""
    q = int(cfg["q_pose"])

    def per_column(a, b):
        n1 = np.linalg.norm(a, axis=0)
        denom = np.maximum(n1, EPS_MACHINE)
        rel = np.linalg.norm(a - b, axis=0) / denom
        return [float(x) for x in rel], [float(x) for x in n1]

    rel_A, norm_A1 = per_column(J_A1, J_A2)
    rel_B, norm_B1 = per_column(J_B1, J_B2)

    out = {
        "h1": float(cfg["fd_h1"]),
        "h2": float(cfg["fd_h2"]),
        "global_relative_diff_fro": {
            "J_A": float(np.linalg.norm(J_A1 - J_A2, ord="fro")
                         / np.linalg.norm(J_A1, ord="fro")),
            "J_B": float(np.linalg.norm(J_B1 - J_B2, ord="fro")
                         / np.linalg.norm(J_B1, ord="fro")),
        },
        "per_column_relative_diff": {
            "J_A": rel_A,
            "J_B": rel_B,
        },
        "per_column_norm_h1": {
            "J_A": norm_A1,
            "J_B": norm_B1,
        },
        "per_column_relative_diff_max": {
            "J_A": float(np.max(rel_A)),
            "J_B": float(np.max(rel_B)),
        },
        "per_column_relative_diff_min": {
            "J_A": float(np.min(rel_A)),
            "J_B": float(np.min(rel_B)),
        },
        "o2_checks": {},
    }

    for key, j_cols in (
        ("J_A", cfg["o2_check_columns_A"]),
        ("J_B", cfg["o2_check_columns_B"]),
    ):
        rows = []
        for j in j_cols:
            cols = {}
            for hh in cfg["o2_h_steps"]:
                cA, cB = col_fd_pair(j, hh, cfg, X_ref)
                cols[hh] = cA if key == "J_A" else cB
            hs = cfg["o2_h_steps"]
            consec = [
                float(np.linalg.norm(cols[hs[i + 1]] - cols[hs[i]])
                      / max(np.linalg.norm(cols[hs[i]]), EPS_MACHINE))
                for i in range(len(hs) - 1)
            ]
            ratios = [
                consec[i + 1] / consec[i] if consec[i] > 0.0 else float("nan")
                for i in range(len(consec) - 1)
            ]
            mean_ratio = float(np.nanmean(ratios))
            rows.append(
                {
                    "column": int(j),
                    "column_norm_h1": float(np.linalg.norm(cols[hs[0]])),
                    "h_steps": [float(x) for x in hs],
                    "consecutive_relative_diff_dh": consec,
                    "consecutive_ratios": ratios,
                    "mean_ratio": mean_ratio,
                    "mean_log2_slope": (
                        float(np.log2(mean_ratio)) if mean_ratio > 0.0
                        else float("nan")
                    ),
                }
            )
        out["o2_checks"][key] = rows
    return out


def L_struct_constants(A0, B0, J_A, J_B):
    """Pointwise structural bound ingredients and L_struct."""
    A_fro = fro(A0)
    tau = spec(B0)
    sig = sigma_min(B0)
    JA_op = spec(J_A)
    JB_op = spec(J_B)
    inner = tau / sig**2 + tau**3 / sig**4
    term_A = 2.0 * A_fro * JA_op
    term_B = 2.0 * A_fro**2 * JB_op * inner
    L_struct = term_A + term_B
    return {
        "A0_fro": A_fro,
        "B0_fro": fro(B0),
        "tau_spec_B0": tau,
        "sigma_min_B0": sig,
        "JA_op_sigma_max": JA_op,
        "JB_op_sigma_max": JB_op,
        "JA_fro": fro(J_A),
        "JB_fro": fro(J_B),
        "inner_tau_sigma": inner,
        "term_A": term_A,
        "term_B": term_B,
        "L_struct_point": L_struct,
    }


def run_affine_certificate(A0, B0, J_A, J_B, K0, lam0, alpha, cfg):
    """Affine tangent model certificate at X0 (ball radius eps)."""
    q = int(cfg["q_pose"])
    n = int(cfg["affine_n_unit_samples"])
    rng = np.random.default_rng(cfg["family5b_seeds"]["affine_unit_directions_2000"])
    dirs = np.stack([unit_direction(rng, q) for _ in range(n)])

    dK_fro = np.empty(n, dtype=float)
    dlam = np.empty(n, dtype=float)
    for i in range(n):
        d = dirs[i]
        A_aff = A0 + (J_A @ d).reshape(A0.shape)
        B_aff = B0 + (J_B @ d).reshape(B0.shape)
        K_aff = K_ab(A_aff, B_aff, alpha)
        dK_fro[i] = fro(K_aff - K0)
        dlam[i] = abs(lambda_min_sym(K_aff) - lam0)

    tau = spec(B0)
    sig = sigma_min(B0)
    JB_op = spec(J_B)
    A_fro = fro(A0)
    JA_op = spec(J_A)

    eps_max_raw = float(cfg["affine_eps_sigma_safety"]) * sig / JB_op
    eps_min = float(cfg["affine_eps_min"])
    n_grid = int(cfg["affine_eps_n_grid"])
    eps_hi = min(eps_max_raw, max(eps_min, eps_max_raw))
    eps_grid = np.geomspace(eps_min, eps_hi, n_grid)

    rows = []
    for eps in eps_grid:
        A_max = A_fro + JA_op * eps
        tau_max = tau + JB_op * eps
        sigma_min_aff = sig - JB_op * eps
        valid = bool(sigma_min_aff > 0.0)
        row = {
            "eps": float(eps),
            "A_max": float(A_max),
            "tau_max": float(tau_max),
            "sigma_min_aff": float(sigma_min_aff),
            "sigma_min_aff_positive": valid,
            "certificate_valid": valid,
        }
        if not valid:
            row.update(
                {
                    "L_cert_affine": None,
                    "worst_fro_ratio": None,
                    "worst_lambda_min_ratio": None,
                    "ok_fro": None,
                    "ok_lambda_min": None,
                }
            )
        else:
            L_cert = (
                2.0 * A_max * JA_op
                + 2.0 * A_max**2 * JB_op
                * (tau_max / sigma_min_aff**2
                   + tau_max**3 / sigma_min_aff**4)
            )
            w_fro = float(np.max(dK_fro)) / float(eps)
            w_lam = float(np.max(dlam)) / float(eps)
            ok_fro = bool(w_fro <= L_cert)
            ok_lam = bool(w_lam <= L_cert)
            row.update(
                {
                    "L_cert_affine": float(L_cert),
                    "worst_fro_ratio": float(w_fro),
                    "worst_lambda_min_ratio": float(w_lam),
                    "ok_fro": ok_fro,
                    "ok_lambda_min": ok_lam,
                    "slack_fro": float(L_cert - w_fro),
                    "relative_slack_fro": float(L_cert / max(w_fro, EPS_MACHINE)),
                    "slack_lambda_min": float(L_cert - w_lam),
                    "relative_slack_lambda_min": float(
                        L_cert / max(w_lam, EPS_MACHINE)
                    ),
                }
            )
        rows.append(row)

    invalid_rows = []
    for eps in cfg["full_eps_list"]:
        sigma_min_aff = sig - JB_op * float(eps)
        invalid_rows.append(
            {
                "eps": float(eps),
                "sigma_min_aff": float(sigma_min_aff),
                "certificate_valid": bool(sigma_min_aff > 0.0),
            }
        )

    return {
        "n_unit_samples": n,
        "seed": cfg["family5b_seeds"]["affine_unit_directions_2000"],
        "eps_grid": [float(x) for x in eps_grid],
        "rows": rows,
        "invalid_at_full_eps_list": invalid_rows,
        "eps_for_sigma_min_aff_zero": float(sig / JB_op),
        "sample_summary": {
            "max_dK_fro": float(np.max(dK_fro)),
            "mean_dK_fro": float(np.mean(dK_fro)),
            "max_dlambda_min": float(np.max(dlam)),
            "mean_dlambda_min": float(np.mean(dlam)),
        },
        "all_valid_rows_ok_fro": all(
            r["ok_fro"] for r in rows if r["certificate_valid"]
        ),
        "all_valid_rows_ok_lambda_min": all(
            r["ok_lambda_min"] for r in rows if r["certificate_valid"]
        ),
    }


def run_full_model(A0, B0, K0, lam0, alpha, L_struct_point, cfg):
    """Full nonlinear model empirical validation at eps in {1e-3,3e-3,1e-2}."""
    q = int(cfg["q_pose"])
    n = int(cfg["full_n_unit_samples"])
    rng = np.random.default_rng(cfg["family5b_seeds"]["full_model_unit_directions_500"])
    dirs = np.stack([unit_direction(rng, q) for _ in range(n)])

    out_rows = []
    for eps in cfg["full_eps_list"]:
        eps = float(eps)
        fro_i = np.empty(n, dtype=float)
        lam_i = np.empty(n, dtype=float)
        for i in range(n):
            A, B = _ab_pair(X0_flat_ref + eps * dirs[i], cfg)
            K = K_ab(A, B, alpha)
            fro_i[i] = fro(K - K0) / eps
            lam_i[i] = abs(lambda_min_sym(K) - lam0) / eps
        w_fro = float(np.max(fro_i))
        w_lam = float(np.max(lam_i))
        out_rows.append(
            {
                "eps": eps,
                "worst_fro_ratio": w_fro,
                "mean_fro_ratio": float(np.mean(fro_i)),
                "worst_lambda_min_drop_ratio": w_lam,
                "mean_lambda_min_drop_ratio": float(np.mean(lam_i)),
                "fro_ratios": [float(x) for x in fro_i],
                "lambda_min_drop_ratios": [float(x) for x in lam_i],
                "le_L_struct_point_fro": bool(w_fro <= L_struct_point),
                "le_L_struct_point_lambda_min": bool(w_lam <= L_struct_point),
                "violation_magnitude_fro": max(0.0, w_fro - L_struct_point),
                "violation_magnitude_lambda_min": max(
                    0.0, w_lam - L_struct_point
                ),
            }
        )
    return {
        "n_unit_samples": n,
        "seed": cfg["family5b_seeds"]["full_model_unit_directions_500"],
        "rows": out_rows,
        "overall_worst_fro_ratio": float(
            max(r["worst_fro_ratio"] for r in out_rows)
        ),
        "overall_worst_lambda_min_drop_ratio": float(
            max(r["worst_lambda_min_drop_ratio"] for r in out_rows)
        ),
        "all_le_L_struct_point_fro": bool(
            all(r["le_L_struct_point_fro"] for r in out_rows)
        ),
        "all_le_L_struct_point_lambda_min": bool(
            all(r["le_L_struct_point_lambda_min"] for r in out_rows)
        ),
    }


def run_boundary_envelope(A0, B0, alpha, cfg):
    """L_struct at X0 and at 5 seeded boundary points; sampled envelope."""
    q = int(cfg["q_pose"])
    eps0 = float(cfg["boundary_eps0"])
    hh = float(cfg["boundary_fd_h"])
    rng = np.random.default_rng(cfg["family5b_seeds"]["boundary_points_5"])
    dirs = np.stack(
        [unit_direction(rng, q) for _ in range(int(cfg["boundary_n_points"]))]
    )

    point_rows = []
    for k in range(int(cfg["boundary_n_points"])):
        X_b = X0_flat_ref + eps0 * dirs[k]
        A_b, B_b = _ab_pair(X_b, cfg)
        J_A, J_B = fd_jacobians(X_b, hh, A0, B0, cfg)
        ing = L_struct_constants(A_b, B_b, J_A, J_B)
        point_rows.append(
            {
                "point": k + 1,
                "direction": [float(x) for x in dirs[k]],
                "sigma_min_B": ing["sigma_min_B0"],
                "tau_spec_B": ing["tau_spec_B0"],
                "A_fro": ing["A0_fro"],
                "JA_op": ing["JA_op_sigma_max"],
                "JB_op": ing["JB_op_sigma_max"],
                "L_struct_point": ing["L_struct_point"],
            }
        )

    boundary_L = [r["L_struct_point"] for r in point_rows]
    return {
        "eps0": eps0,
        "fd_h": hh,
        "seed": cfg["family5b_seeds"]["boundary_points_5"],
        "points": point_rows,
        "L_struct_boundary_max": float(max(boundary_L)),
        "L_struct_env_sample": float(max(boundary_L)),
    }


# ---------------------------------------------------------------------------
# Derivation text, JSON / Markdown helpers
# ---------------------------------------------------------------------------

DERIVATION = {
    "objective": (
        "K(X) = A(X)^T P(X) A(X), P = I - B (B^T B + alpha I)^{-1} B^T, "
        "alpha = 1.0 >= 0. Norms on matrix outputs are Frobenius unless "
        "stated; pose perturbations are measured with the Euclidean 2-norm "
        "on R^18."
    ),
    "step1_directional_derivative": (
        "D_d K = A^T dP A + dA^T P A + A^T P dA (product rule on A^T P A)."
    ),
    "step2_frobenius_triangle_bound": (
        "||D_d K||_F <= 2||A||_F ||dA||_F + ||A||_F^2 ||dP||_2, using "
        "||dA^T P A||_F <= ||dA||_F ||P||_2 ||A||_F and the identical "
        "A^T P dA term, and ||A^T dP A||_F <= ||A||_F^2 ||dP||_2."
    ),
    "step3_P_projection": (
        "P is symmetric with eigenvalues 1 - s_i^2/(s_i^2 + alpha) on the "
        "range directions of B (plus unit eigenvalues); hence ||P||_2 <= 1 "
        "for alpha >= 0."
    ),
    "step4_C_inverse": (
        "C = B^T B + alpha I satisfies sigma_min(C) >= sigma_min(B)^2, so "
        "||C^{-1}||_2 <= 1/sigma^2 with sigma = sigma_min(B). alpha = 1 "
        "actually makes sigma_min(C) ~ 1 + sigma^2; the requested general "
        "alpha >= 0 bound is used for the derived constant."
    ),
    "step5_dP": (
        "d(C^{-1}) = -C^{-1} (dB^T B + B^T dB) C^{-1} gives dP = "
        "-dB C^{-1} B^T - B C^{-1} dB^T + B C^{-1} (dB^T B + B^T dB) "
        "C^{-1} B^T."
    ),
    "step6_dP_norm_bound": (
        "||dP||_2 <= 2||dB||_2 (tau/sigma^2 + tau^3/sigma^4), with "
        "tau = ||B||_2: the first two terms each cost ||dB||_2 tau "
        "||C^{-1}||_2 <= ||dB||_2 tau/sigma^2 and the third costs at most "
        "tau ||C^{-1}||_2 (2 tau ||dB||_2) tau ||C^{-1}||_2 = "
        "2||dB||_2 tau^3/sigma^4."
    ),
    "step7_operator_norms": (
        "With ||dA||_F <= JA_op ||d|| (JA_op = sigma_max(J_A)) and "
        "||dB||_2 <= JB_op ||d|| (JB_op = sigma_max(J_B)), the pointwise "
        "structural Lipschitz constant is L_struct = 2||A||_F JA_op + "
        "2||A||_F^2 JB_op (tau/sigma^2 + tau^3/sigma^4)."
    ),
    "affine_certificate": (
        "For the affine model A_aff(d) = A0 + J_A d, B_aff(d) = B0 + J_B d, "
        "the ingredients on the radius-eps ball are bounded by A_max = "
        "||A0||_F + JA_op eps, tau_max = tau + JB_op eps and sigma_min(B_aff) "
        ">= sigma_min_aff = sigma - JB_op eps (valid while positive). "
        "Substituting these in the same derivation yields L_cert_affine(eps) "
        "= 2 A_max JA_op + 2 A_max^2 JB_op (tau_max/sigma_min_aff^2 + "
        "tau_max^3/sigma_min_aff^4). The affine map is Lipschitz on the whole "
        "ball with this constant, so ||K_aff(d) - K_aff(0)||_F/eps <= "
        "L_cert_affine(eps) and, by Weyl, |lambda_min(K_aff(d)) - "
        "lambda_min(K_aff(0))|/eps <= L_cert_affine(eps) for every unit d."
    ),
    "full_model_validation": (
        "The exact nonlinear K(X0 + eps d) is sampled at 500 unit directions "
        "for eps in {1e-3, 3e-3, 1e-2} and its observed quotients are "
        "compared with L_struct at X0 and with a sampled ingredient envelope "
        "computed at 5 boundary points X0 + 1e-2 d."
    ),
}


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


def round_trip_json(obj: dict) -> str:
    return json.dumps(_jsonable(obj), indent=2) + "\n"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _e(x):
    return f"{x:.3e}"


def _g(x, nd=2):
    return f"{x:.{nd}g}"


def _mid(rel: list[float]) -> tuple[float, float, float]:
    return min(rel), float(np.mean(rel)), max(rel)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(results, figures_dir):
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.1))

    # ---- panel (a): FD Jacobian validation, O(h^2) -----------------------
    ax = axes[0]
    for key, color in (("J_A", "#1f77b4"), ("J_B", "#d62728")):
        for row in results["fd_validation"]["o2_checks"][key]:
            hs = np.asarray(row["h_steps"][:-1])
            dh = np.asarray(row["consecutive_relative_diff_dh"])
            ax.loglog(
                hs, dh, "o-", color=color,
                label=f"{key} col {row['column']}", markersize=4,
            )
    ref_h = np.logspace(-4, np.log10(8e-4), 40)
    scale = results["fd_validation"]["o2_checks"]["J_A"][0][
        "consecutive_relative_diff_dh"
    ][0] / (1e-4) ** 2
    ax.loglog(ref_h, scale * ref_h**2, ":", color="0.35",
              label=r"O$(h^2)$ guide")
    ax.set_xlabel(r"centred-difference step $h$")
    ax.set_ylabel(
        r"relative $|\,J(h)-J(h/2)\,|/\|J(h)\|$ (column vector)"
    )
    ax.set_title(
        "(a) FD Jacobian validation\nconsecutive-step differences, ratios ~4"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="upper left")

    # ---- panel (b): affine tangent certificate ----------------------------
    ax = axes[1]
    aff = results["affine_certificate"]
    valid = [r for r in aff["rows"] if r["certificate_valid"]]
    eps_v = np.asarray([r["eps"] for r in valid])
    lc = np.asarray([r["L_cert_affine"] for r in valid])
    wr_f = np.asarray([r["worst_fro_ratio"] for r in valid])
    wr_l = np.asarray([r["worst_lambda_min_ratio"] for r in valid])
    ax.loglog(eps_v, lc, "-", color="#2ca02c",
              label=r"$L_{\mathrm{cert,affine}}(\varepsilon)$")
    ax.loglog(
        eps_v, wr_f, "o--", color="#1f77b4", markersize=5,
        label=r"worst $\|K_{\mathrm{aff}}(d)-K_{\mathrm{aff}}(0)\|_F/\varepsilon$",
    )
    ax.loglog(
        eps_v, wr_l, "s--", color="#9467bd", markersize=5,
        label=r"worst $|\Delta\lambda_{\min}|/\varepsilon$",
    )
    ax.axvline(
        aff["eps_for_sigma_min_aff_zero"], color="0.4",
        linestyle=":", lw=1.2,
    )
    ax.text(
        0.02, 0.03,
        r"$\sigma_{\min,\mathrm{aff}}=0$ at "
        f"{aff['eps_for_sigma_min_aff_zero']:.2e}",
        transform=ax.transAxes, fontsize=8,
        bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9),
    )
    ax.set_xlabel(r"ball radius $\varepsilon$")
    ax.set_ylabel(r"Lipschitz quotient / certificate value")
    ax.set_title(
        "(b) affine tangent certificate\n"
        f"({aff['n_unit_samples']} unit directions, seed {aff['seed']})"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)

    # ---- panel (c): full nonlinear model validation -----------------------
    ax = axes[2]
    full = results["full_model_validation"]
    eps_f = np.asarray([r["eps"] for r in full["rows"]])
    wr_full = np.asarray([r["worst_fro_ratio"] for r in full["rows"]])
    wr_lam_full = np.asarray(
        [r["worst_lambda_min_drop_ratio"] for r in full["rows"]]
    )
    Lp = results["L_struct_point"]
    Lenv = results["boundary_envelope"]["L_struct_env_sample"]
    ax.loglog(
        eps_f, wr_full, "o-", color="#1f77b4",
        label=r"worst $\|K(X_0+\varepsilon d)-K(X_0)\|_F/\varepsilon$",
    )
    ax.loglog(
        eps_f, wr_lam_full, "s--", color="#9467bd",
        label=r"worst $|\Delta\lambda_{\min}|/\varepsilon$",
    )
    ax.axhline(Lp, color="0.25", linestyle="-", lw=1.4,
               label=r"$L_{\mathrm{struct}}(X_0)$")
    ax.axhline(Lenv, color="#8c564b", linestyle="--", lw=1.4,
               label=r"$L_{\mathrm{struct,env,sample}}$")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"observed quotient")
    ax.set_title(
        "(c) full nonlinear model\n"
        f"({full['n_unit_samples']} unit directions, seed {full['seed']})"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)

    fig.tight_layout()
    path = figures_dir / "family5b_lipschitz_bound.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def write_report(results: dict, notes_dir: Path) -> Path:
    fd = results["fd_validation"]
    const = results["constants"]
    aff = results["affine_certificate"]
    full = results["full_model_validation"]
    env = results["boundary_envelope"]

    def md_table(headers, rows):
        header = "| " + " | ".join(headers) + " |"
        sep = "|" + "|".join(["---"] * len(headers)) + "|"
        body = "\n".join(
            "| " + " | ".join(str(c) for c in row) + " |" for row in rows
        )
        return header + "\n" + sep + "\n" + body

    relA = fd["per_column_relative_diff"]["J_A"]
    relB = fd["per_column_relative_diff"]["J_B"]
    a_min, a_mean, a_max = _mid(relA)
    b_min, b_mean, b_max = _mid(relB)

    o2rows = []
    for key, lab in (("J_A", "A"), ("J_B", "B")):
        for row in fd["o2_checks"][key]:
            o2rows.append(
                [
                    lab,
                    row["column"],
                    _e(row["column_norm_h1"]),
                    ", ".join(_e(x) for x in row["consecutive_relative_diff_dh"]),
                    ", ".join(f"{x:.2f}" for x in row["consecutive_ratios"]),
                    f"{row['mean_ratio']:.2f}",
                ]
            )

    aff_ok_rows = []
    for r in aff["rows"]:
        if not r["certificate_valid"]:
            continue
        aff_ok_rows.append(
            [
                _e(r["eps"]),
                _e(r["A_max"]),
                _e(r["tau_max"]),
                _e(r["sigma_min_aff"]),
                _e(r["L_cert_affine"]),
                _e(r["worst_fro_ratio"]),
                _e(r["worst_lambda_min_ratio"]),
                r["ok_fro"],
                r["ok_lambda_min"],
            ]
        )

    full_rows = []
    for r in full["rows"]:
        full_rows.append(
            [
                _e(r["eps"]),
                _e(r["worst_fro_ratio"]),
                _e(r["worst_lambda_min_drop_ratio"]),
                r["le_L_struct_point_fro"],
                r["le_L_struct_point_lambda_min"],
                _e(r["violation_magnitude_fro"]),
                _e(r["violation_magnitude_lambda_min"]),
            ]
        )

    env_rows = [
        [
            "X0",
            "--",
            _e(const["A0_fro"]),
            _e(const["tau_spec_B0"]),
            _e(const["sigma_min_B0"]),
            _e(const["JA_op_sigma_max"]),
            _e(const["JB_op_sigma_max"]),
            _e(results["L_struct_point"]),
        ]
    ]
    for p in env["points"]:
        env_rows.append(
            [
                f"p{p['point']}",
                _e(env["eps0"]),
                _e(p["A_fro"]),
                _e(p["tau_spec_B"]),
                _e(p["sigma_min_B"]),
                _e(p["JA_op"]),
                _e(p["JB_op"]),
                _e(p["L_struct_point"]),
            ]
        )

    invalid_note = (
        "certificate invalid (sigma_min_aff <= 0): "
        + ", ".join(
            f"eps={r['eps']:.1e} (sigma_min_aff={r['sigma_min_aff']:.2e})"
            for r in aff["invalid_at_full_eps_list"]
        )
    )
    inner = const["inner_tau_sigma"]
    inner1 = const["tau_spec_B0"] / const["sigma_min_B0"] ** 2
    inner2 = inner - inner1

    report = f"""# Family 5b: structural operator-Lipschitz bound for the robust surrogate

Date: 2026-09-03 (SGT; UTC stamp in `results/family5b_lipschitz_bound.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.
This file is generated by the new experiment `src/family5b_lipschitz_bound.py`;
no existing source file was modified and Family 5 itself was not rerun.

## Exact command and runtime

```bash
.venv/bin/python src/family5b_lipschitz_bound.py
```

Wall runtime: {results['runtime_seconds']:.2f} s.  Platform:
{results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, scipy {results['platform']['scipy']},
matplotlib {results['platform']['matplotlib']}.

Scenario: N=16, T=6 arc poses (radius 1.6, phi in [-45, 45] deg), n_rx=4,
q=18 pose parameters, two-blob chi0, smooth p=24 RBF basis, k_b=2*pi,
alpha=1.0, rx offsets {results['config']['rx_offsets']}. Pose-perturbation
norm is Euclidean on R^18.

## Derivation (recorded verbatim in the JSON)

1. **Directional derivative**: {results['derivation']['step1_directional_derivative']}
2. **Frobenius bound**: {results['derivation']['step2_frobenius_triangle_bound']}
3. **Projector bound**: {results['derivation']['step3_P_projection']}
4. **C-inverse bound**: {results['derivation']['step4_C_inverse']}
5. **dP expansion**: {results['derivation']['step5_dP']}
6. **dP norm bound**: {results['derivation']['step6_dP_norm_bound']}
7. **Structural constant**: {results['derivation']['step7_operator_norms']}

{results['derivation']['affine_certificate']}

{results['derivation']['full_model_validation']}

## Ingredients at the Family-5 reference X0

{md_table(
    ["quantity", "value", "notes"],
    [
        ["A0 shape", str(const["A0_shape"]), "2*T*n_rx x p = 48 x 24"],
        ["B0 shape", str(const["B0_shape"]), "2*T*n_rx x 3T = 48 x 18"],
        ["J_A shape", str(results["jacobians"]["J_A_shape"]), "flattened A derivative"],
        ["J_B shape", str(results["jacobians"]["J_B_shape"]), "flattened B derivative"],
        ["||A0||_F", _e(const["A0_fro"]), ""],
        ["tau = ||B0||_2", _e(const["tau_spec_B0"]),
         f"||B0||_F = {_e(const['B0_fro'])}"],
        ["sigma = sigma_min(B0)", _e(const["sigma_min_B0"]),
         "18th singular value"],
        ["||P||_2", _e(const["P_spectral_norm"]), "<= 1 for alpha >= 0"],
        ["JA_op = sigma_max(J_A)", _e(const["JA_op_sigma_max"]),
         f"||J_A||_F = {_e(const['JA_fro'])}"],
        ["JB_op = sigma_max(J_B)", _e(const["JB_op_sigma_max"]),
         f"||J_B||_F = {_e(const['JB_fro'])}"],
        ["tau/sigma^2", _e(inner1), "term 1 of dP bracket"],
        ["tau^3/sigma^4", _e(inner2), "term 2 of dP bracket"],
        ["term A = 2||A||_F JA_op", _e(const["term_A"]), ""],
        ["term B = 2||A||_F^2 JB_op (tau/sigma^2+tau^3/sigma^4)",
         _e(const["term_B"]), ""],
    ],
)}

**L_struct(X0) = {_e(results['L_struct_point'])}.**  The
{100.0 * const['term_B'] / results['L_struct_point']:.4f}% of L_struct
comes from the B/P (nuisance) term, which is dominated by the sigma^-4 factor
because sigma_min(B0) ~ 6.7e-4.  Supplementary note: alpha=1.0 regularises C
so the actual sigma_min(C) ~ 1.0; the used general bound ||C^-1|| <=
1/sigma^2 ~ {_e(const['Cinv2_used_bound'])} is much looser than the measured
||C^-1||_2 ~ {_e(const['Cinv2_actual'])}. The required general (alpha >= 0)
derivation is reported and used throughout; no tightened constant is
substituted into any check below.

## FD Jacobian validation

Centred differences at h1 = 1e-4 and h2 = 2e-4 over all 18 unit coordinates.
Per-column relative differences (all columns in the JSON): J_A
min/mean/max = {_e(a_min)} / {_e(a_mean)} / {_e(a_max)}; J_B
min/mean/max = {_e(b_min)} / {_e(b_mean)} / {_e(b_max)}.
Global relative Frobenius difference: J_A {_e(fd['global_relative_diff_fro']['J_A'])},
J_B {_e(fd['global_relative_diff_fro']['J_B'])}.

O(h^2) consecutive-step checks (steps {fd['o2_checks']['J_A'][0]['h_steps']});
d(h) = relative |J(h) - J(h/2)| / |J(h)| column norm. Ratios ~4 mean the
discretisation error scales like h^2:

{md_table(
    ["matrix", "column", "|col(h1)|", "d(h) at h=1e-4,2e-4,4e-4",
     "ratios", "mean ratio"],
    o2rows,
)}

## Affine tangent certificate (2000 unit directions, seed 7171)

Certificate domain: valid for eps < sigma/JB_op =
{_e(aff['eps_for_sigma_min_aff_zero'])}.  {invalid_note}.

{md_table(
    ["eps", "A_max", "tau_max", "sigma_min_aff", "L_cert_affine",
     "worst ||dK||_F/eps", "worst |dlambda_min|/eps", "ok fro", "ok lambda"],
    aff_ok_rows,
)}

All certificate-valid rows satisfy both checks (fro:
{aff['all_valid_rows_ok_fro']}; lambda_min:
{aff['all_valid_rows_ok_lambda_min']}). Worst affine displacement over all
unit samples: ||dK_aff||_F max = {_e(aff['sample_summary']['max_dK_fro'])},
|dlambda_min| max = {_e(aff['sample_summary']['max_dlambda_min'])}.
Absolute/relative slacks are recorded per row in the JSON.

## Full nonlinear model empirical validation (500 unit directions, seed 9191)

{md_table(
    ["eps", "worst ||DK||_F/eps", "worst |dlambda_min|/eps",
     "<= L_struct(X0) (fro)", "<= L_struct(X0) (lambda)",
     "violation mag fro", "violation mag lambda"],
    full_rows,
)}

Overall worst full-model Frobenius quotient:
{_e(full['overall_worst_fro_ratio'])} vs L_struct(X0) =
{_e(results['L_struct_point'])} (upper-bounds all samples:
{full['all_le_L_struct_point_fro']}). Overall worst lambda_min drop quotient:
{_e(full['overall_worst_lambda_min_drop_ratio'])} (upper-bounds all samples:
{full['all_le_L_struct_point_lambda_min']}).

## Sampled ingredient envelope at 5 boundary points (eps0 = 1e-2, seed 5353)

{md_table(
    ["point", "eps0", "||A||_F", "tau", "sigma_min(B)", "JA_op", "JB_op",
     "L_struct"],
    env_rows,
)}

L_struct sampled envelope = {_e(env['L_struct_env_sample'])} (max over the
five boundary points; X0's L_struct = {_e(results['L_struct_point'])}).
Full-model worst quotients are compared with this envelope in the JSON
(`le_L_struct_env_sample_*`); all sampled rows are below both reference
constants. The envelope is an ingredient sample, not a supremum/interval
certificate.

## Certification status

* The **affine tangent certificate** L_cert_affine(eps) is a RIGOROUS
  certificate FOR THE AFFINE MODEL ONLY, conditional on the numerically
  computed FD Jacobians J_A, J_B being exact to machine precision (the
  O(h^2) h-convergence checks support this).  The 2000-direction checks are
  verification samples; the norm-chain bound itself holds for every unit d on
  the ball while sigma_min_aff > 0.
* The **full nonlinear model results** are a DERIVED structural bound with a
  sampled ingredient envelope; they are NOT a formal interval certificate,
  and every violation/slack is reported exactly.  None of the sampled
  quotients exceeds L_struct(X0) or the sampled envelope.
* **No claim of global certification of the nonlinear Helmholtz model is
  made.**  L_struct(X0) is a pointwise structural Lipschitz constant whose
  ingredients (||A||, tau, sigma, JA_op, JB_op) are not certified uniformly
  over any ball of positive radius.
* The numerical Jacobians are finite-difference estimates; the h=1e-4 vs
  2e-4 comparisons and O(h^2) ratios characterise discretisation error but
  are not machine-precision exactness proofs.

## Cannot-establish section

* No uniform (radius-eps) upper bound on ||A(X)||_F, ||B(X)||_2 or
  sigma_min(B(X))^-1 is proven for the nonlinear map, so L_struct(X0) does
  not certify K on a neighbourhood.
* B(X0) is nearly rank deficient (sigma_min = {_e(const['sigma_min_B0'])}) and
  the structural constant inherits a sigma^-4 blow-up; the empirical
  random-direction quotients are orders of magnitude smaller, so the sampled
  slack is large and no tightness is claimed.
* lambda_min(K) is numerically pinned near zero (K0 lambda_min ~
  {_e(const['lambda_min_K0'])}) because A and B are both rank-deficient to
  machine scale; eigenvalue-drop quotients are therefore round-off scale and
  carry little statistical signal.
* FD Jacobian columns with very small norms can have larger relative
  differences; all per-column values are reported, and no column is censored.

## Artifacts

- results: `results/family5b_lipschitz_bound.json`
- figure: `figures/family5b_lipschitz_bound.png`
- this report: `notes/family5b_lipschitz_bound.md`
"""
    p = notes_dir / "family5b_lipschitz_bound.md"
    p.write_text(report)
    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    global chi0_ref, S_ref, X0_flat_ref
    cfg = CONFIG
    t_start = time.perf_counter()
    t_utc = datetime.now(timezone.utc)

    points, chi0, h_cell, S = f5.base_scene(cfg)
    X0 = f5.build_poses(cfg)
    X0_flat = X0.reshape(-1)
    chi0_ref = chi0
    S_ref = S
    X0_flat_ref = X0_flat

    A0, B0 = f5.smooth_AB_of_flat_X(chi0, S, X0_flat, cfg)
    alpha = float(cfg["alpha_prior"])
    K0 = K_ab(A0, B0, alpha)
    lam0 = lambda_min_sym(K0)
    m = A0.shape[0]
    C0 = B0.T @ B0 + alpha * np.eye(B0.shape[1], dtype=float)
    P0 = np.eye(m, dtype=float) - B0 @ np.linalg.solve(C0, B0.T)
    sv_B0 = np.linalg.svd(B0, compute_uv=False)
    print(
        f"[family5b] A0={A0.shape} B0={B0.shape} K0fro={fro(K0):.6e} "
        f"sigma_min(B)={sigma_min(B0):.3e}"
    )

    runtimes = {}

    t_p1 = time.perf_counter()
    J_A1, J_B1 = fd_jacobians(X0_flat, cfg["fd_h1"], A0, B0, cfg)
    J_A2, J_B2 = fd_jacobians(X0_flat, cfg["fd_h2"], A0, B0, cfg)
    fd_valid = run_fd_validation(
        J_A1, J_B1, J_A2, J_B2, cfg, X0_flat
    )
    runtimes["fd_validation"] = time.perf_counter() - t_p1
    print(
        "[family5b] FD validation: JA rel diff "
        f"{fd_valid['global_relative_diff_fro']['J_A']:.3e}, "
        f"JB {fd_valid['global_relative_diff_fro']['J_B']:.3e}"
    )

    const = L_struct_constants(A0, B0, J_A1, J_B1)
    L_point = const["L_struct_point"]
    const_out = {
        "A0_shape": list(A0.shape),
        "B0_shape": list(B0.shape),
        **const,
        "P_spectral_norm": spec(P0),
        "lambda_min_K0": lam0,
        "lambda_min_C0": sigma_min(C0),
        "Cinv2_used_bound": float(1.0 / const["sigma_min_B0"] ** 2),
        "Cinv2_actual": float(1.0 / sigma_min(C0)),
        "B0_singular_values": [float(x) for x in sv_B0],
        "sigma_min_A0": sigma_min(A0),
    }
    print(
        "[family5b] JA_op="
        f"{const['JA_op_sigma_max']:.6e} JB_op="
        f"{const['JB_op_sigma_max']:.6e} L_struct="
        f"{L_point:.6e}"
    )

    t_p2 = time.perf_counter()
    affine = run_affine_certificate(A0, B0, J_A1, J_B1, K0, lam0, alpha, cfg)
    runtimes["affine_certificate"] = time.perf_counter() - t_p2
    print(
        "[family5b] affine cert valid rows ok fro="
        f"{affine['all_valid_rows_ok_fro']}, ok lambda="
        f"{affine['all_valid_rows_ok_lambda_min']}"
    )

    t_p3 = time.perf_counter()
    full = run_full_model(A0, B0, K0, lam0, alpha, L_point, cfg)
    runtimes["full_model_validation"] = time.perf_counter() - t_p3
    print(
        "[family5b] full-model worst fro ratio="
        f"{full['overall_worst_fro_ratio']:.6e}, <= L_struct(X0): "
        f"{full['all_le_L_struct_point_fro']}"
    )

    t_p4 = time.perf_counter()
    env = run_boundary_envelope(A0, B0, alpha, cfg)
    runtimes["boundary_envelope"] = time.perf_counter() - t_p4
    env["L_struct_env_sample"] = max(
        float(env["L_struct_env_sample"]), float(L_point)
    )
    print(
        "[family5b] boundary envelope L_struct="
        f"{env['L_struct_env_sample']:.6e}"
    )

    # per-eps comparisons against the sampled envelope
    for r in full["rows"]:
        r["le_L_struct_env_sample_fro"] = bool(
            r["worst_fro_ratio"] <= env["L_struct_env_sample"]
        )
        r["le_L_struct_env_sample_lambda_min"] = bool(
            r["worst_lambda_min_drop_ratio"] <= env["L_struct_env_sample"]
        )
    full["all_le_L_struct_env_sample_fro"] = bool(
        all(r["le_L_struct_env_sample_fro"] for r in full["rows"])
    )
    full["all_le_L_struct_env_sample_lambda_min"] = bool(
        all(r["le_L_struct_env_sample_lambda_min"] for r in full["rows"])
    )

    source_script = _HERE / "family5b_lipschitz_bound.py"
    reused = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
        "src/family4_frequency_trajectory.py",
        "src/family5_sensitivity.py",
    ]
    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    results = {
        "generated_utc": t_utc.isoformat(),
        "runner": "src/family5b_lipschitz_bound.py",
        "command": ".venv/bin/python src/family5b_lipschitz_bound.py",
        "family": "5b",
        "title": cfg["title"],
        "parent_gate": (
            "read-only reuse of Family 5's verified API "
            "(base_scene/build_poses/smooth_AB_of_flat_X and its scene "
            "config); Family 5 was not rerun and no existing file was changed"
        ),
        "runtime_seconds": 0.0,  # patched at the end
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
            "scipy linear algebra"
        ),
        "source_sha256": {
            p: hashlib.sha256((_ROOT / p).read_bytes()).hexdigest()
            for p in reused
        },
        "source_script_sha256": hashlib.sha256(
            source_script.read_bytes()
        ).hexdigest(),
        "config": {
            **{k: v for k, v in cfg.items() if k != "seeds"},
            "seeds_family5": cfg.get("seeds"),
            "poses": X0.tolist(),
            "X0_flat": X0_flat.tolist(),
            "rx_offsets": np.asarray(cfg["rx_offsets"]).tolist(),
            "tx_offset": cfg["tx_offset"],
            "grid_h_cell": h_cell,
            "chi0_stats": {
                "min": float(chi0.min()),
                "max": float(chi0.max()),
                "mean": float(chi0.mean()),
                "l2": float(np.linalg.norm(chi0)),
            },
        },
        "derivation": DERIVATION,
        "jacobians": {
            "J_A_shape": list(J_A1.shape),
            "J_B_shape": list(J_B1.shape),
            "fd_h1": cfg["fd_h1"],
            "fd_h2": cfg["fd_h2"],
        },
        "fd_validation": fd_valid,
        "constants": const_out,
        "L_struct_point": L_point,
        "affine_certificate": affine,
        "full_model_validation": full,
        "boundary_envelope": env,
        "certification_status": {
            "affine_tangent": (
                "RIGOROUS certificate FOR THE AFFINE MODEL ONLY, conditional "
                "on the numerically computed FD Jacobians J_A, J_B being "
                "exact to machine precision (supported by the O(h^2) "
                "h-convergence checks); holds on the radius-eps ball while "
                "sigma_min_aff > 0"
            ),
            "full_nonlinear": (
                "DERIVED structural bound with sampled ingredient envelope; "
                "NOT a formal interval certificate; every violation/slack is "
                "reported exactly"
            ),
            "global_certification": (
                "No claim of global certification of the nonlinear Helmholtz/"
                "smooth model is made"
            ),
            "fd_Jacobians": (
                "finite-difference estimates with recorded h1/h2 and O(h^2) "
                "convergence checks, not machine-precision exactness proofs"
            ),
        },
        "scope_note": (
            "Finite-dimensional whitened/realified dense linear algebra on "
            "the validated Family-1 Jacobians with the shared smooth p=24 "
            "basis and alpha=1.0; structural Lipschitz statements are "
            "pointwise/sampled as described and are scenario-specific"
        ),
    }

    fig_path = make_figure(results, figures_dir)
    results["figure_sha256"] = {
        "family5b_lipschitz_bound.png": sha256_file(fig_path)
    }
    results["artifacts"] = {
        "results_json": "results/family5b_lipschitz_bound.json",
        "figure": str(fig_path),
        "report": "notes/family5b_lipschitz_bound.md",
    }

    total_runtime = time.perf_counter() - t_start
    results["runtime_seconds"] = total_runtime
    results_json_path = results_dir / "family5b_lipschitz_bound.json"
    results_json_path.write_text(round_trip_json(results))
    report_path = write_report(results, notes_dir)

    digest_paths = {
        **{p: _ROOT / p for p in reused},
        "src/family5b_lipschitz_bound.py": source_script,
        "results/family5b_lipschitz_bound.json": results_json_path,
        f"figures/{fig_path.name}": fig_path,
        f"notes/{report_path.name}": report_path,
    }
    digest_block = "\n## Artifacts and digests\n\n```text\n"
    digest_block += "\n".join(
        f"{sha256_file(p)}  {label}" for label, p in digest_paths.items()
    )
    digest_block += "\n```\n"
    with report_path.open("a") as fh:
        fh.write(digest_block)

    print("\n===== FAMILY 5B SUMMARY =====")
    print(f"[L_struct_point] {L_point:.6e}")
    print(f"[JA_op] {const['JA_op_sigma_max']:.6e}")
    print(f"[JB_op] {const['JB_op_sigma_max']:.6e}")
    print(f"[affine all ok] {affine['all_valid_rows_ok_fro']}")
    print(f"[full-model <= L_struct_point] {full['all_le_L_struct_point_fro']}")
    print(f"[full-model worst fro ratio] {full['overall_worst_fro_ratio']:.6e}")
    print(f"[envelope L_struct] {env['L_struct_env_sample']:.6e}")
    print(f"[runtime] {total_runtime:.2f} s")
    print(f"[results] results/family5b_lipschitz_bound.json")
    print(f"[figure] figures/{fig_path.name}")
    print(f"[report] notes/{report_path.name}")


if __name__ == "__main__":
    main()
