"""Seed scene/data generation, shared material-only initialisation and
evaluation metrics.  Never touches final seeds except under explicit final
mode (enforced by run_e4)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import minimize

from common import (
    BudgetExhausted,
    CostLedger,
    SQRT2,
    lever_metric_error,
    path_for_seed,
    select_rows,
)


STRATA_ORDER = (
    ("weak", "full"),
    ("weak", "limited"),
    ("strong", "full"),
    ("strong", "limited"),
)
STRENGTH_RANGES = {"weak": (0.15, 0.5), "strong": (0.7, 1.4)}


def stratum_for_seed(seed: int) -> tuple[str, str]:
    return STRATA_ORDER[(int(seed) - 1) % len(STRATA_ORDER)]


def nominal_sigma(model32, alpha_ref: float = 0.5, snr_db: float = 30.0) -> float:
    alpha = np.full(model32.n_alpha, alpha_ref)
    fw = model32.forward(alpha, np.zeros(3), jacobian=False)
    power = float(np.mean(np.abs(fw["total"]) ** 2))
    return float(np.sqrt(power / 10.0 ** (snr_db / 10.0)))


def generate_seed_scene(
    seed: int,
    model32: Any,
    sigma: float,
    rng_seed: int | None = None,
    allow_final: bool = False,
) -> dict[str, Any]:
    path_for_seed(seed, allow_final=allow_final)
    strength, aperture = stratum_for_seed(seed)
    lo, hi = STRENGTH_RANGES[strength]
    rng = np.random.default_rng(rng_seed if rng_seed is not None else int(seed))
    alpha_true = rng.uniform(lo, hi, size=model32.n_alpha)
    x_true = np.zeros(3)
    fw = model32.forward(alpha_true, x_true, jacobian=False)
    y_true = fw["total"]
    noise = (sigma / SQRT2) * (
        rng.standard_normal(y_true.shape)
        + 1j * rng.standard_normal(y_true.shape)
    )
    y = y_true + noise
    noise_val = (sigma / SQRT2) * (
        rng.standard_normal(y_true.shape)
        + 1j * rng.standard_normal(y_true.shape)
    )
    y_val = y_true + noise_val
    return {
        "seed": int(seed),
        "stratum_strength": strength,
        "aperture": aperture,
        "alpha_true": alpha_true.tolist(),
        "x_true": x_true.tolist(),
        "y": y.tolist(),
        "y_val": y_val.tolist(),
        "sigma": float(sigma),
    }


def init_offsets(seed: int, lambda_min: float = 0.5, lever: float = 1.5):
    """Twelve x0 error combinations from A2 (half split lever metric)."""
    radii = lambda_min * np.array([1.0 / 8.0, 1.0 / 2.0, 1.0])
    angles = np.array([0.0, np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0])
    sign = 1.0 if int(seed) % 2 == 0 else -1.0
    out = []
    for r in radii:
        trans = r / SQRT2
        dtheta = sign * r / (SQRT2 * lever)
        for a in angles:
            out.append(
                {
                    "radius_m": float(r),
                    "angle": float(a),
                    "x0": [
                        float(trans * np.cos(a)),
                        float(trans * np.sin(a)),
                        float(dtheta),
                    ],
                }
            )
    return out


def initial_material_fit(
    model: Any,
    y: np.ndarray,
    sigma: float,
    x0: np.ndarray,
    init_max_units: float,
) -> dict[str, Any]:
    """Shared bounded low-frequency material-only fit at the nominal pose."""
    from objective import gaussian_gradient, prior_terms

    ledger = CostLedger()
    ledger.calibration = getattr(model, '_calibration', None)
    alpha0 = np.full(model.n_alpha, 0.5)
    freq_ids = (0,)

    def fg(z: np.ndarray):
        ratio = ledger.calibration.matvec_ratio() if ledger.calibration else 0.0
        if ledger.units + 12 + 60*ratio > init_max_units + 1e-9:
            raise BudgetExhausted("initialisation budget cap")
        pt = prior_terms(z, {"prior_strength": 0.0})
        out = gaussian_gradient(
            model, z, np.asarray(x0, dtype=float), y, sigma, freq_ids, ledger, pt
        )
        g = out["grad_alpha"]
        return float(out["loss"]), g

    res = None
    fallback = False
    try:
        res = minimize(
            fg,
            alpha0,
            method="L-BFGS-B",
            jac=True,
            bounds=[(0.03, 2.5)] * model.n_alpha,
            options={"maxiter": 4, "maxls": 10, "ftol": 1e-6, "gtol": 1e-5},
        )
    except BudgetExhausted:
        fallback = True
    if fallback or res is None or not np.all(np.isfinite(res.x)):
        # predeclared fixed nonzero nominal alpha=0.5 design change
        return {
            "alpha0": alpha0.tolist(),
            "fallback": True,
            "ledger": ledger.snapshot(),
            "units": ledger.units,
            "reason": "init_budget_or_numerical_fallback_to_nominal_alpha_0_5",
        }
    return {
        "alpha0": np.clip(res.x, 0.03, 2.5).tolist(),
        "fallback": False,
        "ledger": ledger.snapshot(),
        "units": ledger.units,
        "reason": None,
    }


def evaluate_estimate(
    model16: Any,
    model32: Any,
    alpha_est: np.ndarray,
    x_est: np.ndarray,
    y: np.ndarray,
    y_val: np.ndarray,
    sigma: float,
    alpha_true: np.ndarray,
    x_true: np.ndarray,
    Q_res: np.ndarray,
    margins: dict[str, float],
    eval_ledger: CostLedger | None = None,
) -> dict[str, Any]:
    """Evaluation-only physical metrics with declared discretisation
    allowance (never feeds estimator decisions)."""
    from objective import coherent_loss

    if eval_ledger is None:
        eval_ledger = CostLedger()
    fw16 = model16.forward(alpha_est, x_est, jacobian=False)
    eval_ledger.charge_model_work(fw16["work"])
    fw32 = model32.forward(alpha_est, x_est, jacobian=False)
    eval_ledger.charge_model_work(fw32["work"])
    r_fit = realify_arr(fw16["total"] - y, sigma)
    loss = 0.5 * float(np.sum(r_fit**2))
    disp = realify_arr(fw16["total"] - fw32["total"], sigma)
    r_val = realify_arr(fw16["total"] - y_val, sigma)
    # whitened validation consistency under the declared discretisation
    # allowance: r_val = disp + noise, so r_val-disp should be unit real noise.
    r_noise = r_val - disp
    val_chi2 = float(np.dot(r_noise, r_noise))
    dof = int(r_noise.size)
    from scipy.stats import chi2

    val_ok = bool(val_chi2 <= chi2.ppf(0.999, df=dof))
    pose_err = lever_metric_error(x_est, x_true)
    a = np.asarray(alpha_est) - np.asarray(alpha_true)
    task_map = np.sqrt(float(np.dot(a, Q_res.T @ Q_res @ a)) / 9.0)
    full_map = float(np.linalg.norm(a)) / np.sqrt(9.0)
    pose_ok = bool(pose_err <= float(margins["d_x_pose_m"]))
    map_ok = bool(task_map <= float(margins["d_chi_task"]))
    return {
        "loss_exact_final": float(loss),
        "residual_fit_norm": float(np.linalg.norm(r_fit)),
        "pose_error_lever_m": float(pose_err),
        "pose_ok": pose_ok,
        "task_map_rmse": float(task_map),
        "full_map_rmse": float(full_map),
        "map_ok": map_ok,
        "validation_chi2": val_chi2,
        "validation_dof": dof,
        "validation_ok": val_ok,
        "discretisation_allowance_norm": float(np.linalg.norm(disp)),
        "eval_ledger": eval_ledger.snapshot(),
    }


def realify_arr(z: np.ndarray, sigma: float) -> np.ndarray:
    z = np.asarray(z)
    return (SQRT2 / sigma) * np.concatenate([np.real(z), np.imag(z)])
