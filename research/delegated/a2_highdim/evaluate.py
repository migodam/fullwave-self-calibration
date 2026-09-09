"""Evaluation-only metrics for a2_highdim estimates.

Truth is used only here and for explicitly labelled oracle runs.  Nothing in
this module feeds an estimator decision.  All forward solves performed here
are kept in a separate evaluation ledger, never in the method work ledger.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.stats import chi2, circstd

from common import (
    CostLedger,
    SQRT2,
    complex_from_json,
    lever_metric_error,
    make_model,
)

ROWS_PER_FREQ = 72


def _wrap_phase(z: np.ndarray) -> np.ndarray:
    return (np.asarray(z) + np.pi) % (2.0 * np.pi) - np.pi


def evaluate_estimate(
    seed: int,
    method: str,
    alpha_est: np.ndarray,
    x_est: np.ndarray,
    scene: dict[str, Any],
    mismatch_kind: str | None = None,
    truth_spatial: np.ndarray | None = None,
) -> tuple[dict[str, Any], np.ndarray | None, dict[str, Any]]:
    """Returns (metrics dict, masked wrapped phase residual array, eval
    ledger snapshot).  ``scene`` provides truth arrays and observations."""
    model_inv = make_model(20, "gaussian49")
    model_data = make_model(32, "gaussian49")
    alpha_est = np.asarray(alpha_est, dtype=float).reshape(-1)
    x_est = np.asarray(x_est, dtype=float).reshape(3)
    sigma = float(scene["sigma"])
    y_true = complex_from_json(scene["y_true"])
    y = complex_from_json(scene["y"])
    y_val = complex_from_json(scene["y_val"])
    x_true = np.asarray(scene["x_true"], dtype=float)
    alpha_true = scene.get("alpha_true")
    if alpha_true is not None and not isinstance(alpha_true, np.ndarray):
        alpha_true = np.asarray(alpha_true, dtype=float)

    eval_ledger = CostLedger()

    # Spatial maps on the N32 data grid.
    Phi32 = model_data.basis
    chi_est = Phi32 @ alpha_est
    if truth_spatial is not None:
        chi_true = np.asarray(truth_spatial, dtype=float).reshape(-1)
        coeff_rmse = None
    elif alpha_true is not None:
        chi_true = Phi32 @ alpha_true
        coeff_rmse = float(np.sqrt(np.mean((alpha_est - alpha_true) ** 2)))
    else:
        raise ValueError("evaluate_estimate needs truth spatial or alpha truth")
    spatial_rmse = float(np.sqrt(np.mean((chi_est - chi_true) ** 2)))
    spatial_nrmse = float(
        np.sqrt(np.mean((chi_est - chi_true) ** 2))
        / max(float(np.sqrt(np.mean(chi_true**2))), 1e-12)
    )
    pose_err = lever_metric_error(x_est, x_true)

    # Full N32 forward predictions (evaluation only).
    t0 = time.perf_counter()
    fw32 = model_data.forward(alpha_est, x_est, jacobian=False)
    eval_ledger.charge_model_work(fw32["work"])
    pred32 = fw32["total"]
    fw20 = model_inv.forward(alpha_est, x_est, jacobian=False)
    eval_ledger.charge_model_work(fw20["work"])
    eval_ledger.add_wall(time.perf_counter() - t0)

    fit_res = realify_vec(pred32 - y, sigma)
    val_res = realify_vec(pred32 - y_val, sigma)
    fit_chi2 = float(np.dot(fit_res, fit_res))
    val_chi2 = float(np.dot(val_res, val_res))
    dof = int(fit_res.size)
    fit_p = float(chi2.sf(fit_chi2, df=dof))
    val_p = float(chi2.sf(val_chi2, df=dof))

    # Wrapped forward-phase residual against the clean truth total field
    # (never the measured phase as ground truth), masked to channels with
    # truth total amplitude > 3 sigma.  Mask is evaluation-only.
    amp_mask = np.abs(y_true) > 3.0 * sigma
    mask_count = int(np.count_nonzero(amp_mask))
    phase_err = None
    phase_metrics: dict[str, Any] = {
        "mask_count": mask_count,
        "mask_definition": "truth_total_amplitude_gt_3sigma",
        "masked_fraction": float(mask_count / y_true.size),
    }
    if mask_count > 0:
        phase_err = _wrap_phase(
            np.angle(pred32[amp_mask]) - np.angle(y_true[amp_mask])
        )
        phase_metrics.update(
            {
                "mean_abs_rad": float(np.mean(np.abs(phase_err))),
                "median_abs_rad": float(np.median(np.abs(phase_err))),
                "rms_rad": float(np.sqrt(np.mean(phase_err**2))),
                "circular_std_rad": float(circstd(phase_err, high=np.pi, low=-np.pi)),
                "p95_abs_rad": float(np.percentile(np.abs(phase_err), 95)),
            }
        )

    metrics: dict[str, Any] = {
        "spatial_map_rmse": spatial_rmse,
        "spatial_map_nrmse": spatial_nrmse,
        "coefficient_rmse": coeff_rmse,
        "pose_error_lever_m": pose_err,
        "pose_error_components_m_rad": x_est.tolist(),
        "fit_whitened_residual_norm": float(np.linalg.norm(fit_res)),
        "fit_chi2": fit_chi2,
        "fit_dof": dof,
        "fit_chi2_p": fit_p,
        "validation_whitened_residual_norm": float(np.linalg.norm(val_res)),
        "validation_chi2": val_chi2,
        "validation_dof": dof,
        "validation_chi2_p": val_p,
        "discretisation_n20_vs_n32_est_norm": float(
            np.linalg.norm(realify_vec(fw20["total"] - pred32, sigma))
        ),
        "phase": phase_metrics,
    }
    if mismatch_kind is not None:
        fitted_ok = bool(fit_p > 0.01)
        val_ok = bool(val_p > 0.05)
        # A joint estimator may silently absorb mismatch into pose: if it
        # reports an in-sample-calibrated residual while pose drifts well
        # beyond the in-basis joint scale (~0.02 m), flag false assurance.
        pose_drift = bool(
            method == "coherent_joint" and pose_err > 0.05
        )
        map_drift = bool(spatial_rmse > 0.10)
        metrics["false_assurance_potential"] = bool(
            fitted_ok and (not val_ok or pose_drift or map_drift)
        )
        metrics["false_assurance_note"] = (
            "descriptive only: an in-sample calibrated residual combined "
            "with a poor validation/map outcome can hide model mismatch; "
            "no calibration guarantee is claimed"
        )
    return metrics, phase_err, eval_ledger.snapshot()


def realify_vec(z: np.ndarray, sigma: float) -> np.ndarray:
    z = np.asarray(z)
    return (SQRT2 / float(sigma)) * np.concatenate([np.real(z), np.imag(z)])


def truth_spatial_map_for_seed(seed: int) -> tuple[np.ndarray, np.ndarray]:
    """N32 data-grid points and true spatial chi map for an in-basis scene."""
    from common import make_scene

    scene = make_scene(seed)
    model_data = make_model(32, "gaussian49")
    alpha_true = np.asarray(scene["alpha_true"], dtype=float)
    return model_data.basis @ alpha_true, model_data.basis
