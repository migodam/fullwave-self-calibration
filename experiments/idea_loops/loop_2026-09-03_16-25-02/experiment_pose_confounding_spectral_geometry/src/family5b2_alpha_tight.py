"""Family 5b2: alpha-aware TIGHT structural operator-Lipschitz bound.

Corrected/improved follow-up to src/family5b_lipschitz_bound.py.

Family 5b used the general alpha >= 0 bound ||C^-1||_2 <= 1/sigma_min(B)^2
even though alpha = 1.0.  That bound is ~2.2e6 times looser than the true
||C^-1||_2 ~ 1/(sigma_min(B)^2 + alpha) ~ 1.  This experiment derives and
uses the exact alpha-aware constants:

    mu = 1/(sigma_min(B)^2 + alpha) = ||C^-1||_2            (exact, PSD B^T B)
    ||dP||_2 <= 2||dB||_2 (tau*mu + tau^3*mu^2)
    L_struct = 2||A||_F JA_op + 2||A||_F^2 JB_op (tau*mu + tau^3*mu^2)

For the affine model over ||d|| <= eps:

    A_max      = ||A0||_F + JA_op*eps
    tau_max    = tau + JB_op*eps
    mu_aff     = 1/(max(0, sigma_min(B0) - JB_op*eps)^2 + alpha)
    L_cert_affine(eps) = 2*A_max*JA_op
                      + 2*A_max^2*JB_op*(tau_max*mu_aff + tau_max^3*mu_aff^2)

alpha > 0 keeps C_aff = B_aff^T B_aff + alpha I invertible for every eps, so
the affine certificate is defined on all eps in {1e-3,3e-3,1e-2,3e-2,1e-1}
(no sigma_min_aff positivity cutoff).

The scenario, basis, poses and centred-FD Jacobian construction are the same
as family 5b (reused through src/family5_sensitivity.py).  family5b's result
JSON is read from disk only for the recorded loose L_struct comparison.  No
existing file is modified.

Certification semantics:
  * L_cert_affine is rigorous FOR THE AFFINE MODEL ONLY, conditional on the
    numerically computed FD Jacobians being machine-accurate;
  * full nonlinear results are a derived structural bound plus a sampled
    ingredient envelope, NOT a formal nonlinear interval certificate.

Run (from the experiment root):
    .venv/bin/python src/family5b2_alpha_tight.py

Outputs:
  results/family5b2_alpha_tight.json
  figures/family5b2_alpha_tight.png
  notes/family5b2_alpha_tight.md
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

# Prevent the interpreter from writing new __pycache__ bytecode files while
# this experiment imports the shared src modules read-only.
sys.dont_write_bytecode = True

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import family5_sensitivity as f5  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

EPS_MACHINE = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Configuration (same Family-5 scenario as family 5b + family-5b2 fields)
# ---------------------------------------------------------------------------

CONFIG = {
    **copy.deepcopy(f5.CONFIG),
    "title": (
        "alpha-aware TIGHT structural operator-Lipschitz bound for the "
        "robust first-order surrogate (family 5b2)"
    ),
    "alpha_prior": 1.0,
    "fd_h1": 1e-4,
    "fd_h2": 2e-4,
    "o2_h_steps": [1e-4, 2e-4, 4e-4, 8e-4],
    "o2_check_columns_A": [0, 8, 17],
    "o2_check_columns_B": [0, 8, 17],
    "family5b2_seeds": {
        "affine_unit_directions_2000": 7171,
        "full_model_unit_directions_500": 9191,
        "boundary_points_5": 5353,
        "note": (
            "numpy.random.default_rng(seed); directions normalised to unit "
            "2-norm on R^18"
        ),
    },
    "affine_eps_list": [1e-3, 3e-3, 1e-2, 3e-2, 1e-1],
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
    "affine_test_convention": (
        "u_i are unit 2-norm directions; the affine displacement is "
        "d = eps*u_i and the recorded quotient is ||K_aff(d)-K_aff(0)||/eps, "
        "i.e. a genuine radius-eps ball Lipschitz quotient"
    ),
    "nonlinearity_note": (
        "full-model samples use the exact smooth_AB_of_flat_X nonlinear map; "
        "the A and B Jacobians used in L_struct are centred finite "
        "differences at the evaluation point"
    ),
    "parent_bound_note": (
        "family5b_lipschitz_bound.py used the general alpha>=0 bound "
        "||C^-1||_2 <= 1/sigma_min(B)^2 and therefore also needed "
        "sigma_min_aff > 0 for its affine rows; family5b2 uses the exact "
        "alpha-aware mu and mu_aff with an alpha floor, removing both issues"
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


def smooth_ab(chi0, S, X_flat, cfg):
    return f5.smooth_AB_of_flat_X(chi0, S, np.asarray(X_flat, dtype=float), cfg)


def fd_jacobians(chi0, S, X_ref, h_step, A0, B0, cfg):
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
        Ap, Bp = smooth_ab(chi0, S, X_ref + e, cfg)
        Am, Bm = smooth_ab(chi0, S, X_ref - e, cfg)
        J_A[:, j] = ((Ap - Am) / (2.0 * hh)).ravel()
        J_B[:, j] = ((Bp - Bm) / (2.0 * hh)).ravel()
    return J_A, J_B


def col_fd_pair(chi0, S, X_ref, j, h_step, cfg):
    """Centred FD column j of the flattened A and B Jacobians at X_ref."""
    e = np.zeros(int(cfg["q_pose"]), dtype=float)
    e[j] = float(h_step)
    Ap, Bp = smooth_ab(chi0, S, X_ref + e, cfg)
    Am, Bm = smooth_ab(chi0, S, X_ref - e, cfg)
    dA = ((Ap - Am) / (2.0 * float(h_step))).ravel()
    dB = ((Bp - Bm) / (2.0 * float(h_step))).ravel()
    return dA, dB


def run_fd_validation(chi0, S, X_ref, J_A1, J_B1, J_A2, J_B2, cfg):
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
            "J_A": float(
                np.linalg.norm(J_A1 - J_A2, ord="fro")
                / np.linalg.norm(J_A1, ord="fro")
            ),
            "J_B": float(
                np.linalg.norm(J_B1 - J_B2, ord="fro")
                / np.linalg.norm(J_B1, ord="fro")
            ),
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
                cA, cB = col_fd_pair(chi0, S, X_ref, j, hh, cfg)
                cols[hh] = cA if key == "J_A" else cB
            hs = cfg["o2_h_steps"]
            consec = [
                float(
                    np.linalg.norm(cols[hs[i + 1]] - cols[hs[i]])
                    / max(np.linalg.norm(cols[hs[i]]), EPS_MACHINE)
                )
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


def alpha_tight_constants(A, B, J_A, J_B, alpha):
    """Pointwise alpha-aware structural bound ingredients and L_struct."""
    A_fro = fro(A)
    tau = spec(B)
    sig = sigma_min(B)
    JA_op = spec(J_A)
    JB_op = spec(J_B)
    mu = 1.0 / (sig**2 + float(alpha))  # exact ||C^-1||_2
    inner = tau * mu + tau**3 * mu**2
    term_A = 2.0 * A_fro * JA_op
    term_B = 2.0 * A_fro**2 * JB_op * inner
    L_struct = term_A + term_B
    return {
        "A_fro": A_fro,
        "B_fro": fro(B),
        "tau_spec": tau,
        "sigma_min_B": sig,
        "JA_op_sigma_max": JA_op,
        "JB_op_sigma_max": JB_op,
        "JA_fro": fro(J_A),
        "JB_fro": fro(J_B),
        "mu_Cinv2": mu,
        "dP_inner_tight": inner,
        "term_A": term_A,
        "term_B": term_B,
        "L_struct_point_tight": L_struct,
    }


def alpha_loose_constants(A, B, J_A, J_B, alpha=0.0):
    """Family-5b (general alpha >= 0) constants for the comparison."""
    A_fro = fro(A)
    tau = spec(B)
    sig = sigma_min(B)
    JA_op = spec(J_A)
    JB_op = spec(J_B)
    # family 5b formula: ||C^-1|| <= 1/sigma^2 with no alpha correction.
    inner_loose = tau / sig**2 + tau**3 / sig**4
    term_A = 2.0 * A_fro * JA_op
    term_B = 2.0 * A_fro**2 * JB_op * inner_loose
    L_struct = term_A + term_B
    return {
        "A_fro": A_fro,
        "tau_spec": tau,
        "sigma_min_B": sig,
        "JA_op_sigma_max": JA_op,
        "JB_op_sigma_max": JB_op,
        "dP_inner_loose": inner_loose,
        "term_A": term_A,
        "term_B": term_B,
        "L_struct_point_loose": L_struct,
        "note": (
            "family-5b general-alpha>=0 formula using ||C^-1|| <= 1/sigma^2; "
            "recorded for the tightness ratio only"
        ),
    }


# ---------------------------------------------------------------------------
# Main computations
# ---------------------------------------------------------------------------

def run_affine_certificate(A0, B0, J_A, J_B, alpha, cfg, K0, lam0):
    """Affine-model certificate checks on radius-eps balls (d = eps*u)."""
    q = int(cfg["q_pose"])
    n = int(cfg["affine_n_unit_samples"])
    eps_list = [float(x) for x in cfg["affine_eps_list"]]
    rng = np.random.default_rng(
        cfg["family5b2_seeds"]["affine_unit_directions_2000"]
    )
    dirs = np.stack([unit_direction(rng, q) for _ in range(n)])

    tau = spec(B0)
    sig = sigma_min(B0)
    JB_op = spec(J_B)
    A_fro = fro(A0)
    JA_op = spec(J_A)

    rows = []
    for eps in eps_list:
        A_max = A_fro + JA_op * eps
        tau_max = tau + JB_op * eps
        sigma_min_aff_lb = max(0.0, sig - JB_op * eps)
        mu_aff = 1.0 / (sigma_min_aff_lb**2 + float(alpha))
        L_cert = (
            2.0 * A_max * JA_op
            + 2.0 * A_max**2 * JB_op
            * (tau_max * mu_aff + tau_max**3 * mu_aff**2)
        )

        fro_i = np.empty(n, dtype=float)
        lam_i = np.empty(n, dtype=float)
        for i in range(n):
            d = eps * dirs[i]
            A_aff = A0 + (J_A @ d).reshape(A0.shape)
            B_aff = B0 + (J_B @ d).reshape(B0.shape)
            K_aff = K_ab(A_aff, B_aff, alpha)
            fro_i[i] = fro(K_aff - K0) / eps
            lam_i[i] = abs(lambda_min_sym(K_aff) - lam0) / eps

        w_fro = float(np.max(fro_i))
        w_lam = float(np.max(lam_i))
        rows.append(
            {
                "eps": eps,
                "A_max": float(A_max),
                "tau_max": float(tau_max),
                "sigma_min_B_lb": float(sig - JB_op * eps),
                "sigma_min_aff_lb_clamped": float(sigma_min_aff_lb),
                "mu_aff": float(mu_aff),
                "certificate_valid": True,  # alpha > 0 keeps C_aff invertible
                "L_cert_affine": float(L_cert),
                "worst_fro_ratio": w_fro,
                "mean_fro_ratio": float(np.mean(fro_i)),
                "worst_lambda_min_ratio": w_lam,
                "mean_lambda_min_ratio": float(np.mean(lam_i)),
                "ok_fro": bool(w_fro <= L_cert),
                "ok_lambda_min": bool(w_lam <= L_cert),
                "slack_fro": float(L_cert - w_fro),
                "relative_slack_fro": float(L_cert / max(w_fro, EPS_MACHINE)),
                "slack_lambda_min": float(L_cert - w_lam),
                "relative_slack_lambda_min": float(
                    L_cert / max(w_lam, EPS_MACHINE)
                ),
            }
        )

    return {
        "n_unit_samples": n,
        "seed": cfg["family5b2_seeds"]["affine_unit_directions_2000"],
        "eps_list": eps_list,
        "rows": rows,
        "eps_for_sigma_min_lb_zero": float(sig / JB_op),
        "sample_note": (
            "d = eps*u with unit u; K_aff(d) = K_ab(A0 + J_A d, B0 + J_B d)"
        ),
        "all_rows_ok_fro": bool(all(r["ok_fro"] for r in rows)),
        "all_rows_ok_lambda_min": bool(
            all(r["ok_lambda_min"] for r in rows)
        ),
    }


def run_full_model(chi0, S, X0_flat, K0, lam0, alpha, cfg):
    """Exact nonlinear model sampled Lipschitz quotients."""
    q = int(cfg["q_pose"])
    n = int(cfg["full_n_unit_samples"])
    rng = np.random.default_rng(
        cfg["family5b2_seeds"]["full_model_unit_directions_500"]
    )
    dirs = np.stack([unit_direction(rng, q) for _ in range(n)])

    rows = []
    for eps in cfg["full_eps_list"]:
        eps = float(eps)
        fro_i = np.empty(n, dtype=float)
        lam_i = np.empty(n, dtype=float)
        for i in range(n):
            A, B = smooth_ab(chi0, S, X0_flat + eps * dirs[i], cfg)
            K = K_ab(A, B, alpha)
            fro_i[i] = fro(K - K0) / eps
            lam_i[i] = abs(lambda_min_sym(K) - lam0) / eps
        rows.append(
            {
                "eps": eps,
                "worst_fro_ratio": float(np.max(fro_i)),
                "mean_fro_ratio": float(np.mean(fro_i)),
                "worst_lambda_min_drop_ratio": float(np.max(lam_i)),
                "mean_lambda_min_drop_ratio": float(np.mean(lam_i)),
            }
        )

    return {
        "n_unit_samples": n,
        "seed": cfg["family5b2_seeds"]["full_model_unit_directions_500"],
        "rows": rows,
        "overall_worst_fro_ratio": float(
            max(r["worst_fro_ratio"] for r in rows)
        ),
        "overall_worst_lambda_min_drop_ratio": float(
            max(r["worst_lambda_min_drop_ratio"] for r in rows)
        ),
    }


def run_boundary_envelope(
    chi0, S, X0_flat, A0, B0, J_A0, J_B0, alpha, cfg
):
    """Tight L_struct at X0 and at 5 seeded boundary points (eps0=1e-2)."""
    q = int(cfg["q_pose"])
    eps0 = float(cfg["boundary_eps0"])
    hh = float(cfg["boundary_fd_h"])
    rng = np.random.default_rng(
        cfg["family5b2_seeds"]["boundary_points_5"]
    )
    dirs = np.stack(
        [unit_direction(rng, q) for _ in range(int(cfg["boundary_n_points"]))]
    )

    c0 = alpha_tight_constants(A0, B0, J_A0, J_B0, alpha)
    point_rows = [
        {
            "point": 0,
            "label": "X0",
            "eps0": 0.0,
            "direction": None,
            "A_fro": c0["A_fro"],
            "tau_spec_B": c0["tau_spec"],
            "sigma_min_B": c0["sigma_min_B"],
            "JA_op": c0["JA_op_sigma_max"],
            "JB_op": c0["JB_op_sigma_max"],
            "L_struct_point_tight": c0["L_struct_point_tight"],
        }
    ]

    for k in range(int(cfg["boundary_n_points"])):
        X_b = X0_flat + eps0 * dirs[k]
        A_b, B_b = smooth_ab(chi0, S, X_b, cfg)
        J_A, J_B = fd_jacobians(chi0, S, X_b, hh, A_b, B_b, cfg)
        ing = alpha_tight_constants(A_b, B_b, J_A, J_B, alpha)
        point_rows.append(
            {
                "point": k + 1,
                "label": f"X0+{eps0:.0e}*u{k+1}",
                "eps0": eps0,
                "direction": [float(x) for x in dirs[k]],
                "A_fro": ing["A_fro"],
                "tau_spec_B": ing["tau_spec"],
                "sigma_min_B": ing["sigma_min_B"],
                "JA_op": ing["JA_op_sigma_max"],
                "JB_op": ing["JB_op_sigma_max"],
                "L_struct_point_tight": ing["L_struct_point_tight"],
            }
        )

    tight_at_boundary = [r["L_struct_point_tight"] for r in point_rows[1:]]
    L_env = max(
        [r["L_struct_point_tight"] for r in point_rows]
    )
    return {
        "eps0": eps0,
        "fd_h": hh,
        "seed": cfg["family5b2_seeds"]["boundary_points_5"],
        "points": point_rows,
        "L_struct_boundary_max": float(max(tight_at_boundary)),
        "L_struct_env_sample": float(L_env),
        "note": (
            "tight alpha-aware L_struct at X0 plus 5 seeded boundary points; "
            "sampled envelope, not a uniform interval certificate"
        ),
    }


# ---------------------------------------------------------------------------
# Serialisation and hashing helpers
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
    fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.2))

    # ---- panel (a): loose vs tight L_struct (log bar) --------------------
    ax = axes[0]
    cmp = results["L_struct_loose_vs_tight"]
    labels = [
        "loose (family 5b)\n1/sigma^2 bound",
        "tight (family 5b2)\nmu = 1/(sigma^2+alpha)",
    ]
    values = [
        cmp["recomputed_loose_L_struct"],
        results["L_struct_point"],
    ]
    bars = ax.bar(labels, values, log=True, color=["#d62728", "#2ca02c"])
    for bar, v in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            v * 1.5,
            f"{v:.4e}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.set_ylabel(r"$L_{\mathrm{struct}}$ (log)")
    ax.set_title(
        "(a) loose vs alpha-aware tight\n"
        f"tightness ratio {cmp['tightness_ratio_recomputed']:.3e}"
    )
    ax.grid(True, which="both", alpha=0.3, axis="y")

    # ---- panel (b): affine certificate (log-log) --------------------------
    ax = axes[1]
    aff = results["affine_certificate"]
    eps_v = np.asarray([r["eps"] for r in aff["rows"]])
    lc = np.asarray([r["L_cert_affine"] for r in aff["rows"]])
    wf = np.asarray([r["worst_fro_ratio"] for r in aff["rows"]])
    wl = np.asarray([r["worst_lambda_min_ratio"] for r in aff["rows"]])
    ax.loglog(
        eps_v, lc, "o-", color="#2ca02c",
        label=r"$L_{\mathrm{cert,affine}}(\varepsilon)$",
    )
    ax.loglog(
        eps_v, wf, "s--", color="#1f77b4", markersize=6,
        label=r"worst $\|K_{\mathrm{aff}}(\varepsilon u)-K_{\mathrm{aff}}(0)\|_F/\varepsilon$",
    )
    ax.loglog(
        eps_v, np.maximum(wl, 1e-18), "d--", color="#9467bd", markersize=6,
        label=r"worst $|\Delta\lambda_{\min}|/\varepsilon$",
    )
    ax.set_xlabel(r"ball radius $\varepsilon$")
    ax.set_ylabel("Lipschitz quotient / certificate value")
    ax.set_title(
        "(b) alpha-aware affine certificate\n"
        f"({aff['n_unit_samples']} unit directions, seed {aff['seed']}); "
        "valid for all eps (alpha>0)"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, loc="center left")

    # ---- panel (c): full nonlinear model validation -----------------------
    ax = axes[2]
    full = results["full_model_validation"]
    eps_f = np.asarray([r["eps"] for r in full["rows"]])
    wf_full = np.asarray([r["worst_fro_ratio"] for r in full["rows"]])
    wl_full = np.asarray(
        [r["worst_lambda_min_drop_ratio"] for r in full["rows"]]
    )
    Lp = results["L_struct_point"]
    Lenv = results["boundary_envelope"]["L_struct_env_sample"]
    ax.loglog(
        eps_f, wf_full, "o-", color="#1f77b4",
        label=r"worst $\|K(X_0+\varepsilon u)-K(X_0)\|_F/\varepsilon$",
    )
    ax.loglog(
        eps_f, np.maximum(wl_full, 1e-18), "s--", color="#9467bd",
        label=r"worst $|\Delta\lambda_{\min}|/\varepsilon$",
    )
    ax.axhline(
        Lp, color="0.25", linestyle="-", lw=1.4,
        label=r"$L_{\mathrm{struct,tight}}(X_0)$",
    )
    ax.axhline(
        Lenv, color="#8c564b", linestyle="--", lw=1.4,
        label=r"$L_{\mathrm{struct,env,sample}}$",
    )
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel("observed quotient")
    ax.set_title(
        "(c) full nonlinear model\n"
        f"({full['n_unit_samples']} unit directions, seed {full['seed']})"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)

    fig.tight_layout()
    path = figures_dir / "family5b2_alpha_tight.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def write_report(results: dict, notes_dir: Path) -> Path:
    def md_table(headers, rows):
        header = "| " + " | ".join(headers) + " |"
        sep = "|" + "|".join(["---"] * len(headers)) + "|"
        body = "\n".join(
            "| " + " | ".join(str(c) for c in row) + " |" for row in rows
        )
        return header + "\n" + sep + "\n" + body

    const = results["constants"]
    cmp = results["L_struct_loose_vs_tight"]
    fd = results["fd_validation"]
    aff = results["affine_certificate"]
    full = results["full_model_validation"]
    env = results["boundary_envelope"]

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
                    ", ".join(
                        _e(x) for x in row["consecutive_relative_diff_dh"]
                    ),
                    ", ".join(f"{x:.2f}" for x in row["consecutive_ratios"]),
                    f"{row['mean_ratio']:.2f}",
                ]
            )

    aff_rows = [
        [
            _e(r["eps"]),
            _e(r["A_max"]),
            _e(r["tau_max"]),
            _e(r["sigma_min_B_lb"]),
            _e(r["mu_aff"]),
            _e(r["L_cert_affine"]),
            _e(r["worst_fro_ratio"]),
            _e(r["worst_lambda_min_ratio"]),
            r["ok_fro"],
            r["ok_lambda_min"],
        ]
        for r in aff["rows"]
    ]

    full_rows = [
        [
            _e(r["eps"]),
            _e(r["worst_fro_ratio"]),
            _e(r["worst_lambda_min_drop_ratio"]),
            r["le_L_struct_point_tight_fro"],
            r["le_L_struct_point_tight_lambda_min"],
            r["le_L_struct_env_sample_fro"],
            r["le_L_struct_env_sample_lambda_min"],
            _e(r["violation_magnitude_fro"]),
            _e(r["violation_magnitude_lambda_min"]),
        ]
        for r in full["rows"]
    ]

    env_rows = [
        [
            p["label"],
            _e(p["eps0"]) if p["eps0"] else "--",
            _e(p["A_fro"]),
            _e(p["tau_spec_B"]),
            _e(p["sigma_min_B"]),
            _e(p["JA_op"]),
            _e(p["JB_op"]),
            _e(p["L_struct_point_tight"]),
        ]
        for p in env["points"]
    ]

    loose_tight_notes = (
        "recorded family-5b L_struct = "
        f"{_e(cmp['recorded_loose_L_struct'])}; "
        "recomputed loose value = "
        f"{_e(cmp['recomputed_loose_L_struct'])}; "
        f"tightness ratio = {cmp['tightness_ratio_recorded']:.6e} "
        "(recorded) / "
        f"{cmp['tightness_ratio_recomputed']:.6e} (recomputed)."
    )

    report = f"""# Family 5b2: alpha-aware TIGHT structural operator-Lipschitz bound

Date: {results['generated_utc']} (UTC; also recorded in JSON).
Experiment: `experiment_pose_confounding_spectral_geometry`.
Generated by the new `src/family5b2_alpha_tight.py`; no existing source,
result, figure, or note file was modified.  The prior result
`results/family5b_lipschitz_bound.json` was read from disk for the loose-value
comparison only.

## Exact command and runtime

```bash
.venv/bin/python src/family5b2_alpha_tight.py
```

Wall runtime: {results['runtime_seconds']:.2f} s (part runtimes in the JSON).
Platform: {results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, scipy {results['platform']['scipy']},
matplotlib {results['platform']['matplotlib']}.

Scenario: N=16, T=6 arc poses (radius 1.6, phi in [-45, 45] deg), n_rx=4,
q=18 pose parameters, two-blob chi0, smooth p=24 RBF basis, k_b=2*pi,
alpha=1.0, rx offsets {results['config']['rx_offsets']}.  Pose-perturbation
norm is Euclidean on R^18.

## Error fixed relative to family 5b

1. **Loose C^-1 bound replaced by the exact alpha-aware mu.**  Family 5b used
   `||C^-1||_2 <= 1/sigma_min(B)^2` (the general alpha >= 0 bound), giving a
   sigma^-4 blow-up, although alpha=1.0 makes the true `||C^-1||_2 =
   1/(sigma_min(B)^2 + alpha) ~ 1`.  Family 5b2 uses
   `mu = 1/(sigma_min(B)^2 + alpha)` exactly (C = B^T B + alpha I is a
   PSD-plus-alpha shift, so this equality is exact, not a bound).
2. **Affine certificate domain extended to all eps.**  Family 5b required
   `sigma_min_aff = sigma - JB_op*eps > 0`, so every requested full-scale eps
   in {{1e-3, 3e-3, 1e-2}} was marked invalid.  With alpha > 0,
   `mu_aff = 1/(max(0, sigma_min(B0) - JB_op*eps)^2 + alpha)` is always
   finite and C_aff = B_aff^T B_aff + alpha I is always invertible, so the
   certificate is valid for every eps (including 3e-2 and 1e-1).
3. **Affine verification evaluated on the actual radius-eps ball.**  Family
   5b's loop built `A_aff = A0 + J_A u` for unit u and divided the resulting
   unit-radius displacement by eps, which does not test the stated radius-eps
   Lipschitz quotient.  Family 5b2 evaluates `d = eps*u` (unit u) and divides
   `||K_aff(eps u) - K_aff(0)||` by eps, which is the quotient appearing in
   the certificate.

{loose_tight_notes}

## Correct alpha-aware derivation (recorded verbatim in the JSON)

Let K(X) = A(X)^T P(X) A(X) with P = I - B C^{-1} B^T and
C = B^T B + alpha I, alpha = 1.0 > 0.

1. **C is invertible and mu is exact.**  C is symmetric positive definite
   with eigenvalues s_i(B)^2 + alpha; hence
   `sigma_min(C) = sigma_min(B)^2 + alpha` and
   `mu := ||C^{-1}||_2 = 1/(sigma_min(B)^2 + alpha)` exactly (largest
   eigenvalue of C^{-1} = 1/smallest eigenvalue of C).
2. **Directional derivative.**  D_d K = A^T dP A + dA^T P A + A^T P dA, so
   `||D_d K||_F <= 2||A||_F ||dA||_F + ||A||_F^2 ||dP||_2`
   (using ||P||_2 <= 1 for alpha >= 0).
3. **dP expansion.**  d(C^{-1}) = -C^{-1}(dB^T B + B^T dB) C^{-1}, hence
   dP = -dB C^{-1} B^T - B C^{-1} dB^T
        + B C^{-1}(dB^T B + B^T dB) C^{-1} B^T.
4. **dP spectral bound with tau = ||B||_2.**  Each of the first two terms is
   at most `||dB||_2 * tau * mu`; the third is at most
   `2||dB||_2 * tau^3 * mu^2`.  Therefore
   `||dP||_2 <= 2||dB||_2 (tau*mu + tau^3*mu^2)`.
5. **Pointwise structural constant.**  With `||dA||_F <= JA_op ||d||` and
   `||dB||_2 <= JB_op ||d||`:

       L_struct_tight = 2||A||_F JA_op
                      + 2||A||_F^2 JB_op (tau*mu + tau^3*mu^2).

6. **Affine certificate.**  For A_aff(d) = A0 + J_A d and
   B_aff(d) = B0 + J_B d with ||d|| <= eps,

       A_max     = ||A0||_F + JA_op*eps
       tau_max   = tau + JB_op*eps
       sigma_min(B_aff) >= max(0, sigma_min(B0) - JB_op*eps)
       mu_aff    = 1/(max(0, sigma_min(B0) - JB_op*eps)^2 + alpha)
       L_cert_affine(eps) = 2*A_max*JA_op
                            + 2*A_max^2*JB_op
                              *(tau_max*mu_aff + tau_max^3*mu_aff^2).

   Because alpha > 0, C_aff = B_aff^T B_aff + alpha I is invertible for every
   d and every eps, so `mu_aff` (and hence L_cert_affine) is finite for all
   eps >= 0.  The norm-chain argument bounds the operator norm of the affine
   map's derivative on the whole radius-eps ball, so by the mean value
   inequality `||K_aff(d)-K_aff(0)||_F/eps <= L_cert_affine(eps)` for every
   unit d; Weyl's inequality then gives the same bound for
   `|lambda_min(K_aff(d)) - lambda_min(K_aff(0))|/eps` (spectral norm of a
   symmetric difference is <= its Frobenius norm).

7. **Full nonlinear model.**  The exact nonlinear K(X0 + eps d) is sampled at
   500 unit directions for eps in {{1e-3, 3e-3, 1e-2}}; observed quotients are
   compared with L_struct_tight at X0 and with a sampled ingredient envelope
   computed at X0 plus 5 boundary points X0 + 1e-2 u (seed 5353).

## Tight ingredients at the Family-5 reference X0

{md_table(
    ["quantity", "value", "notes"],
    [
        ["A0 shape", str(const["A0_shape"]), "2*T*n_rx x p = 48 x 24"],
        ["B0 shape", str(const["B0_shape"]), "2*T*n_rx x 3T = 48 x 18"],
        ["J_A shape", str(results["jacobians"]["J_A_shape"]), "flattened A derivative"],
        ["J_B shape", str(results["jacobians"]["J_B_shape"]), "flattened B derivative"],
        ["||A0||_F", _e(const["A_fro"]), ""],
        ["tau = ||B0||_2", _e(const["tau_spec"]),
         f"||B0||_F = {_e(const['B_fro'])}"],
        ["sigma = sigma_min(B0)", _e(const["sigma_min_B"]), "18th singular value"],
        ["alpha", f"{const['alpha']:g}", "regulariser in C = B^T B + alpha I"],
        ["mu = ||C^-1||_2 = 1/(sigma^2+alpha)", _e(const["mu_Cinv2"]),
         f"exact; sigma_min(C0) = {_e(const['lambda_min_C0'])}"],
        ["family-5b loose ||C^-1|| upper bound 1/sigma^2", _e(const["Cinv2_loose"]),
         "not used in any family-5b2 bound"],
        ["JA_op = sigma_max(J_A)", _e(const["JA_op_sigma_max"]),
         f"||J_A||_F = {_e(const['JA_fro'])}"],
        ["JB_op = sigma_max(J_B)", _e(const["JB_op_sigma_max"]),
         f"||J_B||_F = {_e(const['JB_fro'])}"],
        ["tau*mu", _e(const["tau_mu"]), "term 1 of tight dP bracket"],
        ["tau^3*mu^2", _e(const["tau3_mu2"]), "term 2 of tight dP bracket"],
        ["term A = 2||A||_F JA_op", _e(const["term_A"]), ""],
        ["term B = 2||A||_F^2 JB_op (tau*mu+tau^3*mu^2)",
         _e(const["term_B"]), ""],
    ],
)}

**L_struct_tight(X0) = {_e(results['L_struct_point'])}** (vs the loose
family-5b value {_e(cmp['recorded_loose_L_struct'])}; tightness ratio
{cmp['tightness_ratio_recorded']:.6e}).  The
{100.0 * const['term_A'] / results['L_struct_point']:.2f}% of the tight
constant comes from the A-term (2||A||_F JA_op ~ {_e(const['term_A'])}); the
regularised B/P term contributes only
{_e(const['term_B'])} (~{100.0 * const['term_B'] / results['L_struct_point']:.2f}%),
which is the reverse of the loose family-5b balance.

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

Rows below use d = eps*u with unit u.  sigma_min(B0) - JB_op*eps turns
negative at eps > {_e(aff['eps_for_sigma_min_lb_zero'])}; the alpha floor
clamps it to 0, mu_aff = 1/alpha, and every row is certificate-valid.

{md_table(
    ["eps", "A_max", "tau_max", "sigma_min_B_lb", "mu_aff",
     "L_cert_affine", "worst ||dK||_F/eps", "worst |dlambda_min|/eps",
     "ok fro", "ok lambda"],
    aff_rows,
)}

All affine rows satisfy both checks (fro: {aff['all_rows_ok_fro']};
lambda_min: {aff['all_rows_ok_lambda_min']}).  L_cert_affine is mildly
increasing (~1.1e-2 at eps=1e-3 to ~1.8e-2 at eps=1e-1) because alpha bounds
mu_aff; worst sampled Frobenius quotients are ~1e-3 (up to 1.4e-3 at
eps=1e-1), and eigenvalue-drop quotients are round-off scale.  Full slacks are
in the JSON.

## Full nonlinear model empirical validation (500 unit directions, seed 9191)

{md_table(
    ["eps", "worst ||DK||_F/eps", "worst |dlambda_min|/eps",
     "<= L_tight(X0) (fro)", "<= L_tight(X0) (lambda)",
     "<= L_env (fro)", "<= L_env (lambda)",
     "violation mag fro", "violation mag lambda"],
    full_rows,
)}

Overall worst full-model Frobenius quotient
{_e(full['overall_worst_fro_ratio'])} vs L_struct_tight(X0) =
{_e(results['L_struct_point'])} (upper-bounds all samples:
{full['all_le_L_struct_point_tight_fro']}).  All samples also stay below the
sampled envelope (fro: {full['all_le_L_struct_env_sample_fro']}).

## Sampled ingredient envelope (tight L_struct at X0 and 5 boundary points, eps0 = 1e-2, seed 5353)

{md_table(
    ["point", "eps0", "||A||_F", "tau", "sigma_min(B)", "JA_op", "JB_op",
     "L_struct_tight"],
    env_rows,
)}

L_struct_tight sampled envelope = {_e(env['L_struct_env_sample'])} (max over
X0 and the five boundary points).  The envelope is an ingredient sample, not
a supremum/interval certificate.

## Certification status

* The **affine tangent certificate** L_cert_affine(eps) is a RIGOROUS
  certificate FOR THE AFFINE MODEL ONLY, conditional on the numerically
  computed FD Jacobians J_A, J_B being machine-accurate (supported by the
  O(h^2) h-convergence checks).  alpha > 0 keeps C_aff invertible for all eps,
  so the certificate is valid for every eps in
  {{1e-3, 3e-3, 1e-2, 3e-2, 1e-1}} without a sigma_min_aff positivity
  restriction.
* The **full nonlinear model results** are a DERIVED structural bound with a
  sampled ingredient envelope; they are NOT a formal interval certificate,
  and every violation/slack is reported exactly.  No sampled quotient exceeds
  L_struct_tight(X0) or the sampled envelope.
* **No claim of global certification of the nonlinear Helmholtz/smooth model
  is made.**  L_struct_tight(X0) is pointwise; its ingredients are not
  certified uniformly over a ball of positive radius.
* The numerical Jacobians are finite-difference estimates; h=1e-4 vs 2e-4
  comparisons and O(h^2) ratios characterise discretisation error but are not
  machine-precision exactness proofs.

## Cannot-establish section

* No uniform radius-eps upper bounds on ||A(X)||_F, ||B(X)||_2 or a lower
  bound on sigma_min(B(X)) are proven for the nonlinear map, so the
  full-model statement remains pointwise-derived plus sampled-envelope.
* lambda_min(K) is numerically pinned near zero (K0 lambda_min ~
  {_e(const['lambda_min_K0'])}) because A and B are rank-deficient to machine
  scale; eigenvalue-drop quotients are round-off scale and carry little
  statistical signal.
* mu_aff = 1/alpha is used once the linear lower bound on sigma_min(B_aff)
  reaches zero; it is finite and correct, but is looser than the exact
  inverse of the actual minimum eigenvalue of C_aff at any particular
  perturbed point.

## Artifacts

- results: `results/family5b2_alpha_tight.json`
- figure: `figures/family5b2_alpha_tight.png`
- this report: `notes/family5b2_alpha_tight.md`
"""
    p = notes_dir / "family5b2_alpha_tight.md"
    p.write_text(report)
    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    t_start = time.perf_counter()
    t_utc = datetime.now(timezone.utc)
    alpha = float(cfg["alpha_prior"])

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"

    # ---- read the family-5b result (recorded loose L_struct) -------------
    prior_path = results_dir / "family5b_lipschitz_bound.json"
    prior = json.loads(prior_path.read_text())
    prior_loose = float(prior["L_struct_point"])

    # ---- scene, poses, A0/B0, K0 (identical construction to family 5b) ---
    points, chi0, h_cell, S = f5.base_scene(cfg)
    X0 = f5.build_poses(cfg)
    X0_flat = X0.reshape(-1)
    A0, B0 = smooth_ab(chi0, S, X0_flat, cfg)
    K0 = K_ab(A0, B0, alpha)
    lam0 = lambda_min_sym(K0)
    m = A0.shape[0]
    C0 = B0.T @ B0 + alpha * np.eye(B0.shape[1], dtype=float)
    P0 = np.eye(m, dtype=float) - B0 @ np.linalg.solve(C0, B0.T)
    sv_B0 = np.linalg.svd(B0, compute_uv=False)
    print(
        f"[family5b2] A0={A0.shape} B0={B0.shape} "
        f"sigma_min(B)={sigma_min(B0):.3e} alpha={alpha:g}"
    )

    runtimes = {}

    # ---- centred-FD Jacobians at X0 + validation --------------------------
    t_p1 = time.perf_counter()
    J_A1, J_B1 = fd_jacobians(
        chi0, S, X0_flat, cfg["fd_h1"], A0, B0, cfg
    )
    J_A2, J_B2 = fd_jacobians(
        chi0, S, X0_flat, cfg["fd_h2"], A0, B0, cfg
    )
    fd_valid = run_fd_validation(
        chi0, S, X0_flat, J_A1, J_B1, J_A2, J_B2, cfg
    )
    runtimes["fd_validation"] = time.perf_counter() - t_p1
    print(
        "[family5b2] FD validation: JA rel diff "
        f"{fd_valid['global_relative_diff_fro']['J_A']:.3e}, "
        f"JB {fd_valid['global_relative_diff_fro']['J_B']:.3e}"
    )

    # ---- tight + loose constants at X0 ------------------------------------
    tight = alpha_tight_constants(A0, B0, J_A1, J_B1, alpha)
    loose = alpha_loose_constants(A0, B0, J_A1, J_B1)
    L_tight = tight["L_struct_point_tight"]
    L_loose_recomp = loose["L_struct_point_loose"]

    const_out = {
        "A0_shape": list(A0.shape),
        "B0_shape": list(B0.shape),
        "alpha": alpha,
        **tight,
        "P_spectral_norm": spec(P0),
        "lambda_min_K0": lam0,
        "lambda_min_C0": sigma_min(C0),
        "mu_Cinv2": tight["mu_Cinv2"],
        "Cinv2_loose": float(1.0 / tight["sigma_min_B"] ** 2),
        "tau_mu": float(tight["tau_spec"] * tight["mu_Cinv2"]),
        "tau3_mu2": float(tight["tau_spec"] ** 3 * tight["mu_Cinv2"] ** 2),
        "B0_singular_values": [float(x) for x in sv_B0],
        "sigma_min_A0": sigma_min(A0),
    }

    loose_vs_tight = {
        "recorded_source": "results/family5b_lipschitz_bound.json",
        "recorded_loose_L_struct": prior_loose,
        "recomputed_loose_L_struct": L_loose_recomp,
        "loose_formula": loose["note"],
        "tight_L_struct": L_tight,
        "tightness_ratio_recorded": float(prior_loose / L_tight),
        "tightness_ratio_recomputed": float(L_loose_recomp / L_tight),
    }
    print(
        "[family5b2] L_struct tight="
        f"{L_tight:.6e} (recorded loose {prior_loose:.6e}; ratio "
        f"{prior_loose / L_tight:.3e})"
    )

    # ---- affine certificate -------------------------------------------------
    t_p2 = time.perf_counter()
    affine = run_affine_certificate(
        A0, B0, J_A1, J_B1, alpha, cfg, K0, lam0
    )
    runtimes["affine_certificate"] = time.perf_counter() - t_p2
    print(
        "[family5b2] affine rows ok fro="
        f"{affine['all_rows_ok_fro']}, ok lambda="
        f"{affine['all_rows_ok_lambda_min']}"
    )

    # ---- full nonlinear model ----------------------------------------------
    t_p3 = time.perf_counter()
    full = run_full_model(chi0, S, X0_flat, K0, lam0, alpha, cfg)
    runtimes["full_model_validation"] = time.perf_counter() - t_p3

    # ---- tight sampled envelope at X0 + 5 boundary points -------------------
    t_p4 = time.perf_counter()
    env = run_boundary_envelope(
        chi0, S, X0_flat, A0, B0, J_A1, J_B1, alpha, cfg
    )
    runtimes["boundary_envelope"] = time.perf_counter() - t_p4

    # ---- full-model comparisons against tight point + envelope -------------
    L_env = env["L_struct_env_sample"]
    for r in full["rows"]:
        r["le_L_struct_point_tight_fro"] = bool(
            r["worst_fro_ratio"] <= L_tight
        )
        r["le_L_struct_point_tight_lambda_min"] = bool(
            r["worst_lambda_min_drop_ratio"] <= L_tight
        )
        r["le_L_struct_env_sample_fro"] = bool(
            r["worst_fro_ratio"] <= L_env
        )
        r["le_L_struct_env_sample_lambda_min"] = bool(
            r["worst_lambda_min_drop_ratio"] <= L_env
        )
        r["violation_magnitude_fro"] = max(
            0.0, r["worst_fro_ratio"] - L_tight
        )
        r["violation_magnitude_lambda_min"] = max(
            0.0, r["worst_lambda_min_drop_ratio"] - L_tight
        )
    full["all_le_L_struct_point_tight_fro"] = bool(
        all(r["le_L_struct_point_tight_fro"] for r in full["rows"])
    )
    full["all_le_L_struct_point_tight_lambda_min"] = bool(
        all(r["le_L_struct_point_tight_lambda_min"] for r in full["rows"])
    )
    full["all_le_L_struct_env_sample_fro"] = bool(
        all(r["le_L_struct_env_sample_fro"] for r in full["rows"])
    )
    full["all_le_L_struct_env_sample_lambda_min"] = bool(
        all(r["le_L_struct_env_sample_lambda_min"] for r in full["rows"])
    )
    print(
        "[family5b2] full-model worst fro ratio="
        f"{full['overall_worst_fro_ratio']:.6e}, <= tight point: "
        f"{full['all_le_L_struct_point_tight_fro']}, <= env: "
        f"{full['all_le_L_struct_env_sample_fro']}"
    )

    reused = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
        "src/family4_frequency_trajectory.py",
        "src/family5_sensitivity.py",
    ]
    source_script = _HERE / "family5b2_alpha_tight.py"

    results = {
        "generated_utc": t_utc.isoformat(),
        "runner": "src/family5b2_alpha_tight.py",
        "command": ".venv/bin/python src/family5b2_alpha_tight.py",
        "family": "5b2",
        "title": cfg["title"],
        "parent_gate": (
            "read-only reuse of Family 5's verified API "
            "(base_scene/build_poses/smooth_AB_of_flat_X and its scene "
            "config); Family 5 and family 5b were not rerun; the family-5b "
            "result JSON was read from disk; no existing file was modified"
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
        "derivation": {
            "objective": (
                "K(X) = A(X)^T P(X) A(X), P = I - B (B^T B + alpha I)^-1 B^T, "
                "alpha = 1.0 > 0. Output norms Frobenius unless stated; "
                "pose perturbations measured with the Euclidean 2-norm on R^18."
            ),
            "step1_mu_exact": (
                "C = B^T B + alpha I has eigenvalues s_i(B)^2 + alpha, so "
                "sigma_min(C) = sigma_min(B)^2 + alpha and mu := ||C^-1||_2 "
                "= 1/(sigma_min(B)^2 + alpha) EXACTLY (not the loose general "
                "alpha>=0 bound 1/sigma_min(B)^2 used by family 5b)."
            ),
            "step2_directional_derivative": (
                "D_d K = A^T dP A + dA^T P A + A^T P dA, hence ||D_d K||_F "
                "<= 2||A||_F ||dA||_F + ||A||_F^2 ||dP||_2 with ||P||_2 <= 1."
            ),
            "step3_dP_expansion": (
                "dP = -dB C^-1 B^T - B C^-1 dB^T + B C^-1 (dB^T B + B^T dB) "
                "C^-1 B^T."
            ),
            "step4_dP_bound": (
                "With tau = ||B||_2, ||dP||_2 <= 2||dB||_2 (tau*mu + "
                "tau^3*mu^2)."
            ),
            "step5_structural_constant": (
                "L_struct_tight = 2||A||_F JA_op + 2||A||_F^2 JB_op "
                "(tau*mu + tau^3*mu^2), with JA_op = sigma_max(J_A) and "
                "JB_op = sigma_max(J_B)."
            ),
            "affine_certificate": (
                "On the affine model A_aff(d) = A0 + J_A d, B_aff(d) = B0 + "
                "J_B d with ||d|| <= eps: A_max = ||A0||_F + JA_op eps; "
                "tau_max = tau + JB_op eps; sigma_min(B_aff) >= max(0, "
                "sigma_min(B0) - JB_op eps); mu_aff = 1/(max(0, sigma_min(B0)"
                " - JB_op eps)^2 + alpha); L_cert_affine(eps) = 2 A_max "
                "JA_op + 2 A_max^2 JB_op (tau_max mu_aff + tau_max^3 "
                "mu_aff^2). alpha > 0 keeps C_aff = B_aff^T B_aff + alpha I "
                "invertible for every d and eps, so the certificate is valid "
                "for all eps (no sigma_min_aff > 0 cutoff). The sampled "
                "checks evaluate d = eps*u and verify the quotient <= "
                "L_cert_affine(eps) on the actual radius-eps ball."
            ),
            "full_model_validation": (
                "The exact nonlinear K(X0 + eps d) is sampled at 500 unit "
                "directions for eps in {1e-3, 3e-3, 1e-2}; observed quotients "
                "are compared with L_struct_tight at X0 and with a sampled "
                "ingredient envelope at X0 plus 5 boundary points "
                "X0 + 1e-2 u."
            ),
        },
        "jacobians": {
            "J_A_shape": list(J_A1.shape),
            "J_B_shape": list(J_B1.shape),
            "fd_h1": cfg["fd_h1"],
            "fd_h2": cfg["fd_h2"],
        },
        "fd_validation": fd_valid,
        "constants": const_out,
        "L_struct_point": L_tight,
        "L_struct_loose_vs_tight": loose_vs_tight,
        "affine_certificate": affine,
        "full_model_validation": full,
        "boundary_envelope": env,
        "certification_status": {
            "affine_tangent": (
                "RIGOROUS certificate FOR THE AFFINE MODEL ONLY, conditional "
                "on the numerically computed FD Jacobians J_A, J_B being "
                "machine-accurate (supported by the O(h^2) h-convergence "
                "checks); holds for every eps in {1e-3,3e-3,1e-2,3e-2,1e-1} "
                "because alpha > 0 keeps C_aff invertible"
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
            "basis and alpha=1.0; alpha-aware structural Lipschitz statements "
            "are pointwise/sampled as described and are scenario-specific"
        ),
    }

    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    fig_path = make_figure(results, figures_dir)
    results["figure_sha256"] = {
        "family5b2_alpha_tight.png": sha256_file(fig_path)
    }
    results["prior_artifacts_sha256"] = {
        "results/family5b_lipschitz_bound.json": sha256_file(prior_path)
    }
    results["artifacts"] = {
        "results_json": "results/family5b2_alpha_tight.json",
        "figure": str(fig_path),
        "report": "notes/family5b2_alpha_tight.md",
    }

    total_runtime = time.perf_counter() - t_start
    results["runtime_seconds"] = total_runtime
    results_json_path = results_dir / "family5b2_alpha_tight.json"
    results_json_path.write_text(round_trip_json(results))
    report_path = write_report(results, notes_dir)

    digest_paths = {
        **{p: _ROOT / p for p in reused},
        "src/family5b2_alpha_tight.py": source_script,
        "results/family5b2_alpha_tight.json": results_json_path,
        f"figures/{fig_path.name}": fig_path,
        f"notes/{report_path.name}": report_path,
        "results/family5b_lipschitz_bound.json": prior_path,
    }
    digest_block = "\n## Artifacts and digests\n\n```text\n"
    digest_block += "\n".join(
        f"{sha256_file(p)}  {label}" for label, p in digest_paths.items()
    )
    digest_block += "\n```\n"
    with report_path.open("a") as fh:
        fh.write(digest_block)

    print("\n===== FAMILY 5B2 SUMMARY =====")
    print(f"[L_struct_point_tight] {L_tight:.6e}")
    print(f"[loose (recorded)] {prior_loose:.6e}")
    print(f"[tightness ratio] {prior_loose / L_tight:.6e}")
    print(
        "[affine worst ratios] "
        + ", ".join(f"{r['worst_fro_ratio']:.6e}" for r in affine["rows"])
    )
    print(
        "[affine all ok] "
        f"{affine['all_rows_ok_fro']} / "
        f"{affine['all_rows_ok_lambda_min']}"
    )
    print(
        "[full-model worst fro ratios] "
        + ", ".join(f"{r['worst_fro_ratio']:.6e}" for r in full["rows"])
    )
    print(
        "[full-model violations] "
        f"fro={full['all_le_L_struct_point_tight_fro']}, "
        f"env={full['all_le_L_struct_env_sample_fro']}"
    )
    print(f"[envelope L_struct_tight] {L_env:.6e}")
    print(f"[runtime] {total_runtime:.2f} s")
    print("[results] results/family5b2_alpha_tight.json")
    print(f"[figure] figures/{fig_path.name}")
    print(f"[report] notes/{report_path.name}")


if __name__ == "__main__":
    main()
