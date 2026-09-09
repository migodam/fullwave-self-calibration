"""Development joint material/pose solvers for the A3 ROM study.

Methods implemented here share one exact physical acceptance controller:
every trial is evaluated with the full A2 forward model and only an exact
loss decrease is accepted.  Direct baselines are

- ``direct_gn``: analytic-Jacobian damped Gauss-Newton (QR normal step,
  exact trial acceptance);
- ``direct_adjoint``: exact-adjoint L-BFGS-B baseline.

ROM methods use the frozen-chart fixed LS proxy.  The chart is frozen for
every derivative/line-search trial, rebuilt between accepted iterates with
full cost accounting, and reported separately from reduced solves.  All ROM
methods share the same trial-acceptance controller and damping schedule.
"""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from scipy.linalg import qr, solve_triangular
from scipy.optimize import minimize

from common import (
    FREQ_IDS,
    POSE_SCALE,
    WorkLedger,
    exact_loss,
    lever_metric_error,
    params_from_scaled,
    realify,
    scaled_from_params,
    select_rows,
    z_bounds,
    complex_from_json,
)
from basis import build_chart
from wave import fixed_chart_evaluate

DAMPING_SCHEDULE = (1e-14, 1e-9, 1e-6, 1e-3, 1e-1, 10.0, 1e3, 1e5)
ARMJO_C1 = 1e-4

CONT_STAGES = ((0,), (0, 1), (0, 1, 3))
CONT_MAXITER = (8, 8, 14)
ALL_STAGES = ((0, 1, 3),)
ALL_MAXITER = (24,)

# Predeclared development chart budgets (grid16).  These come from the
# offline rank audit on development scenes only; they are not final-seed
# tuning and every chart's actual rank is recorded.
CHART_RANKS: dict[str, dict[int, int]] = {
    "sensing": {0: 36, 1: 36, 3: 36},
    "twofold": {0: 96, 1: 128, 3: 192},
    "generic_task": {0: 64, 1: 96, 3: 128},
    "sensing_task": {0: 96, 1: 96, 3: 128},
    "twofold_task": {0: 96, 1: 96, 3: 128},
    "block_krylov": {0: 96, 1: 192, 3: 224},
}

ROM_METHODS = tuple(CHART_RANKS)


class BudgetExhausted(RuntimeError):
    pass


def _z_clip(z: np.ndarray) -> np.ndarray:
    bounds = z_bounds()
    out = z.copy()
    for i, (lo, hi) in enumerate(bounds):
        out[i] = min(max(out[i], lo), hi)
    return out


def _loss_jacobian_from_forward(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: WorkLedger,
) -> tuple[float, np.ndarray, np.ndarray, dict[str, Any]]:
    fw = model.forward(alpha, x, list(freq_ids), jacobian=True)
    ledger.charge_model_work(fw["work"])
    y_sel = select_rows(y, freq_ids)
    res_complex = fw["total"] - y_sel
    loss = exact_loss(res_complex, sigma)
    J = np.hstack([fw["A"], fw["B"] / POSE_SCALE[None, :]])
    J = realify(J, sigma)
    r = realify(res_complex, sigma)
    return loss, J, r, fw


def _exact_loss_forward(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: WorkLedger,
) -> tuple[float, dict[str, Any]]:
    fw = model.forward(alpha, x, list(freq_ids), jacobian=False)
    ledger.charge_model_work(fw["work"])
    ledger.charge_acceptance(1)
    y_sel = select_rows(y, freq_ids)
    return exact_loss(fw["total"] - y_sel, sigma), fw


def _adjoint_loss_gradient(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: WorkLedger,
) -> tuple[float, np.ndarray]:
    y_sel = select_rows(y, freq_ids)
    out = model.adjoint_gradient(alpha, x, y_sel, sigma, list(freq_ids))
    ledger.charge_model_work(out["work"])
    g = np.concatenate(
        [
            np.asarray(out["grad_alpha"]),
            np.asarray(out["grad_x"]) / POSE_SCALE,
        ]
    )
    return float(out["loss"]), g


def _gn_trial_step(
    J: np.ndarray,
    r: np.ndarray,
    lam: float,
    ledger: WorkLedger | None = None,
) -> tuple[np.ndarray, float]:
    """Damped GN step by QR of the stacked system (no normal inverse)."""
    m, d = J.shape
    A = np.vstack([J, np.sqrt(max(lam, 1e-300)) * np.eye(d)])
    b = np.concatenate([-r, np.zeros(d)])
    Q, R = qr(A, mode="economic", check_finite=False)
    if ledger is not None:
        ledger.charge_reduced(qr=1, dims=A.shape)
    step = solve_triangular(R, Q.T @ b, lower=False, check_finite=False)
    return step, 0.0


def accept_trial(Lcur: float, Lnew: float, step_norm: float) -> bool:
    if not (np.isfinite(Lnew) and np.isfinite(Lcur)):
        return False
    if step_norm < 1e-12:
        return False
    threshold = 1e-8 * max(1.0, float(abs(Lcur)))
    return float(Lnew) <= float(Lcur) - threshold


def _reduce_method_outcome_counts(events: list[dict[str, Any]]) -> dict[str, Any]:
    by_freq: dict[int, int] = {}
    n_reduced = 0
    n_rejected = 0
    n_fallback = 0
    for ev in events:
        if ev.get("accepted_reduced"):
            n_reduced += 1
            for fi in ev.get("freq_ids", []):
                by_freq[int(fi)] = by_freq.get(int(fi), 0) + 1
        if ev.get("rejected_reduced"):
            n_rejected += 1
        if ev.get("fallback_direct"):
            n_fallback += 1
    return {
        "accepted_reduced_total": n_reduced,
        "accepted_reduced_by_freq": by_freq,
        "rejected_reduced_trials": n_rejected,
        "fallback_direct_steps": n_fallback,
    }


def _run_direct_gn(
    model: Any,
    y: np.ndarray,
    sigma: float,
    alpha0: np.ndarray,
    x0: np.ndarray,
    freq_stages: tuple[tuple[int, ...], ...],
    stage_maxiter: tuple[int, ...],
    ledger: WorkLedger,
    budget_rhs: int,
) -> dict[str, Any]:
    alpha = alpha0.copy()
    x = x0.copy()
    events: list[dict[str, Any]] = []
    status = "schedule_complete"
    t_start = time.perf_counter()

    def est_rhs(n_freq: int, jac: bool) -> int:
        return n_freq * (78 if jac else 6)

    for si, stage_freqs in enumerate(freq_stages):
        stage_loss, _ = _exact_loss_forward(
            model, alpha, x, y, sigma, stage_freqs, ledger
        )
        for _it in range(stage_maxiter[si]):
            if ledger.full_rhs_columns + est_rhs(len(stage_freqs), True) > budget_rhs:
                status = "budget_stopped"
                break
            Lcur, J, r, _ = _loss_jacobian_from_forward(
                model, alpha, x, y, sigma, stage_freqs, ledger
            )
            event: dict[str, Any] = {"stage": si, "freq_ids": list(stage_freqs)}
            accepted = False
            for damp in DAMPING_SCHEDULE:
                step, _ = _gn_trial_step(J, r, damp, ledger)
                trial = _z_clip(scaled_from_params(alpha, x) + step)
                trial_alpha, trial_x = params_from_scaled(trial)
                if np.max(np.abs(trial - scaled_from_params(alpha, x))) < 1e-13:
                    continue
                if ledger.full_rhs_columns + est_rhs(len(stage_freqs), False) > budget_rhs:
                    status = "budget_stopped"
                    break
                Lnew, _ = _exact_loss_forward(
                    model,
                    trial_alpha,
                    trial_x,
                    y,
                    sigma,
                    stage_freqs,
                    ledger,
                )
                if accept_trial(Lcur, Lnew, float(np.linalg.norm(step))):
                    alpha, x = trial_alpha, trial_x
                    accepted = True
                    event.update(
                        {"damping": damp, "loss": Lcur, "loss_trial": Lnew}
                    )
                    events.append(event)
                    break
            if status != "schedule_complete":
                break
            if not accepted:
                status = "no_improving_trial"
                event["no_accept"] = True
                events.append(event)
                break
            if ledger.full_rhs_columns > budget_rhs:
                status = "budget_stopped"
                break
        if status != "schedule_complete":
            break
    return {
        "status": status,
        "alpha": alpha,
        "x": x,
        "events": events,
        "wall_seconds": time.perf_counter() - t_start,
    }


def _run_direct_adjoint(
    model: Any,
    y: np.ndarray,
    sigma: float,
    alpha0: np.ndarray,
    x0: np.ndarray,
    freq_stages: tuple[tuple[int, ...], ...],
    stage_maxiter: tuple[int, ...],
    ledger: WorkLedger,
    budget_rhs: int,
) -> dict[str, Any]:
    alpha = alpha0.copy()
    x = x0.copy()
    z = scaled_from_params(alpha, x)
    accepted = z.copy()
    status = "schedule_complete"
    t_start = time.perf_counter()

    def fg(zc: np.ndarray) -> tuple[float, np.ndarray]:
        key = bytes(np.asarray(zc, dtype=float).tobytes())
        if key in cache:
            return cache[key]
        a, xx = params_from_scaled(zc)
        loss, g = _adjoint_loss_gradient(
            model, a, xx, y, sigma, stage_freqs, ledger
        )
        if ledger.full_rhs_columns > budget_rhs:
            raise BudgetExhausted("adjoint RHS budget exceeded")
        cache[key] = (loss, g)
        return cache[key]

    def callback(zk: np.ndarray) -> None:
        nonlocal accepted
        accepted = np.asarray(zk, dtype=float).copy()

    for si, stage_freqs in enumerate(freq_stages):
        cache: dict[bytes, tuple[float, np.ndarray]] = {}
        try:
            res = minimize(
                fg,
                z,
                method="L-BFGS-B",
                jac=True,
                bounds=z_bounds(),
                options={
                    "maxiter": stage_maxiter[si],
                    "maxls": 24,
                    "ftol": 1e-12,
                    "gtol": 1e-7,
                },
                callback=callback,
            )
            if bool(getattr(res, "success", False)) and np.all(np.isfinite(res.x)):
                accepted = np.asarray(res.x).copy()
            z = accepted.copy()
        except BudgetExhausted:
            z = accepted.copy()
            status = "budget_stopped"
            break
    alpha, x = params_from_scaled(accepted)
    return {
        "status": status,
        "alpha": alpha,
        "x": x,
        "events": [],
        "wall_seconds": time.perf_counter() - t_start,
    }


def _run_rom(
    model: Any,
    y: np.ndarray,
    sigma: float,
    alpha0: np.ndarray,
    x0: np.ndarray,
    chart_method: str,
    freq_stages: tuple[tuple[int, ...], ...],
    stage_maxiter: tuple[int, ...],
    ledger: WorkLedger,
    budget_rhs: int,
) -> dict[str, Any]:
    alpha = alpha0.copy()
    x = x0.copy()
    events: list[dict[str, Any]] = []
    chart_builds = 0
    chart_ranks_seen: dict[int, list[int]] = {}
    status = "schedule_complete"
    t_start = time.perf_counter()
    chart: dict[int, np.ndarray] | None = None

    def est_rhs(n_freq: int, jac: bool) -> int:
        return n_freq * (78 if jac else 6)

    for si, stage_freqs in enumerate(freq_stages):
        # Exact loss at the accepted stage-start point.
        Lcur, _ = _exact_loss_forward(
            model, alpha, x, y, sigma, stage_freqs, ledger
        )
        for _it in range(stage_maxiter[si]):
            if ledger.full_rhs_columns + est_rhs(len(stage_freqs), True) > budget_rhs:
                status = "budget_stopped"
                break
            # Rebuild the frozen chart at the current accepted physical
            # point.  This cost is charged to the run every accepted iterate.
            chart = build_chart(
                model,
                alpha,
                x,
                stage_freqs,
                chart_method,
                256,
                y,
                ledger,
                include_tangent_rhs=True,
                rank_by_freq=CHART_RANKS[chart_method],
            )
            chart_builds += 1
            for fi in stage_freqs:
                if fi in chart:
                    chart_ranks_seen.setdefault(int(fi), []).append(chart[fi].shape[1])
            if not chart:
                status = "empty_chart"
                break

            red = fixed_chart_evaluate(
                model,
                alpha,
                x,
                chart,
                y,
                sigma,
                stage_freqs,
                jacobian=True,
                ledger=ledger,
            )
            y_sel = select_rows(y, stage_freqs)
            res = red["total"] - y_sel
            r = realify(res, sigma)
            A = red["A"]
            B = red["B"]
            J = np.hstack([A, B / POSE_SCALE[None, :]])
            J = realify(J, sigma)
            gdot = None
            accepted = False
            event: dict[str, Any] = {
                "stage": si,
                "freq_ids": list(stage_freqs),
                "proxy_loss": red["loss"],
                "max_state_res_rel": red["max_state_res_rel"],
            }
            for damp in DAMPING_SCHEDULE:
                step, _ = _gn_trial_step(J, r, damp, ledger)
                trial = _z_clip(scaled_from_params(alpha, x) + step)
                trial_alpha, trial_x = params_from_scaled(trial)
                if np.max(np.abs(trial - scaled_from_params(alpha, x))) < 1e-13:
                    continue
                if ledger.full_rhs_columns + est_rhs(len(stage_freqs), False) > budget_rhs:
                    status = "budget_stopped"
                    break
                Lnew, _ = _exact_loss_forward(
                    model,
                    trial_alpha,
                    trial_x,
                    y,
                    sigma,
                    stage_freqs,
                    ledger,
                )
                if accept_trial(Lcur, Lnew, float(np.linalg.norm(step))):
                    alpha, x = trial_alpha, trial_x
                    Lcur = Lnew
                    accepted = True
                    event.update(
                        {
                            "accepted_reduced": True,
                            "damping": damp,
                            "loss_cur": Lcur,
                            "loss_trial": Lnew,
                            "reduced_step_norm": float(np.linalg.norm(step)),
                            "ranks": [chart[fi].shape[1] for fi in stage_freqs if fi in chart],
                        }
                    )
                    events.append(event)
                    break
                event["rejected_reduced"] = True
                event.setdefault("damping_attempts", []).append(damp)
            if status != "schedule_complete":
                break
            if accepted:
                continue

            # Direct fallback at the frozen/rejected chart: one exact analytic
            # GN step with the same acceptance controller.
            Lcur, Jfull, rfull, _ = _loss_jacobian_from_forward(
                model, alpha, x, y, sigma, stage_freqs, ledger
            )
            fallback_ok = False
            for damp in DAMPING_SCHEDULE:
                step, _ = _gn_trial_step(Jfull, rfull, damp, ledger)
                trial = _z_clip(scaled_from_params(alpha, x) + step)
                trial_alpha, trial_x = params_from_scaled(trial)
                if np.max(np.abs(trial - scaled_from_params(alpha, x))) < 1e-13:
                    continue
                if ledger.full_rhs_columns + est_rhs(len(stage_freqs), False) > budget_rhs:
                    status = "budget_stopped"
                    break
                Lnew, _ = _exact_loss_forward(
                    model,
                    trial_alpha,
                    trial_x,
                    y,
                    sigma,
                    stage_freqs,
                    ledger,
                )
                if accept_trial(Lcur, Lnew, float(np.linalg.norm(step))):
                    alpha, x = trial_alpha, trial_x
                    Lcur = Lnew
                    fallback_ok = True
                    event.update(
                        {
                            "fallback_direct": True,
                            "damping": damp,
                            "loss_cur": Lcur,
                            "loss_trial": Lnew,
                        }
                    )
                    events.append(event)
                    break
            if not fallback_ok:
                status = "no_improving_trial_and_fallback"
                event["fallback_failed"] = True
                events.append(event)
                break
        if status != "schedule_complete":
            break

    counts = _reduce_method_outcome_counts(events)
    return {
        "status": status,
        "alpha": alpha,
        "x": x,
        "events": events,
        "chart_builds": chart_builds,
        "chart_ranks_seen": {str(k): v for k, v in chart_ranks_seen.items()},
        "wall_seconds": time.perf_counter() - t_start,
        **counts,
    }


def run_solver(
    model: Any,
    scene: dict[str, Any],
    method: str,
    continuation: bool,
    seed: int,
    budget_rhs: int = 5000,
) -> dict[str, Any]:
    """Run one development method under the shared controller."""
    if method not in ("direct_gn", "direct_adjoint", *ROM_METHODS):
        raise ValueError(f"unknown method {method!r}")
    y = complex_from_json(scene["y"])
    sigma = float(scene["sigma"])
    alpha_true = np.asarray(scene["alpha_true"], dtype=float)
    x_true = np.asarray(scene["x_true"], dtype=float)
    alpha0 = np.full(model.n_alpha, 0.08)
    from common import pose_init

    x0 = pose_init(seed)
    freq_stages = CONT_STAGES if continuation else ALL_STAGES
    maxiters = CONT_MAXITER if continuation else ALL_MAXITER
    ledger = WorkLedger()
    t_start = time.perf_counter()
    if method == "direct_gn":
        outcome = _run_direct_gn(
            model,
            y,
            sigma,
            alpha0,
            x0,
            freq_stages,
            maxiters,
            ledger,
            budget_rhs,
        )
    elif method == "direct_adjoint":
        outcome = _run_direct_adjoint(
            model,
            y,
            sigma,
            alpha0,
            x0,
            freq_stages,
            maxiters,
            ledger,
            budget_rhs,
        )
    else:
        outcome = _run_rom(
            model,
            y,
            sigma,
            alpha0,
            x0,
            method,
            freq_stages,
            maxiters,
            ledger,
            budget_rhs,
        )
    alpha_end = outcome["alpha"]
    x_end = outcome["x"]

    # Final exact audit on the full three-frequency development data.
    loss_end, fw_end = _exact_loss_forward(
        model, alpha_end, x_end, y, sigma, FREQ_IDS, ledger
    )
    loss_init, _ = _exact_loss_forward(
        model, alpha0, x0, y, sigma, FREQ_IDS, ledger
    )
    alpha_rmse = float(
        np.sqrt(np.mean((np.asarray(alpha_end) - np.asarray(alpha_true)) ** 2))
    )
    pose_err = lever_metric_error(np.asarray(x_end), np.asarray(x_true))
    counts = (
        _reduce_method_outcome_counts(outcome.get("events", []))
        if method not in ("direct_gn", "direct_adjoint")
        else {
            "accepted_reduced_total": 0,
            "accepted_reduced_by_freq": {},
            "rejected_reduced_trials": 0,
            "fallback_direct_steps": 0,
        }
    )
    return {
        "seed": int(seed),
        "method": method,
        "continuation": bool(continuation),
        "status": outcome["status"],
        "final_loss": loss_end,
        "initial_loss_all_freq": loss_init,
        "alpha_rmse": alpha_rmse,
        "pose_metric_error_m": pose_err,
        "x_end": list(map(float, np.asarray(x_end))),
        "alpha_end": list(map(float, np.asarray(alpha_end))),
        "chart_builds": outcome.get("chart_builds", 0),
        "chart_ranks_seen": outcome.get("chart_ranks_seen", {}),
        "reduced": counts,
        "events": outcome.get("events", []),
        "ledger": ledger.snapshot(),
        "wall_seconds": time.perf_counter() - t_start,
    }
