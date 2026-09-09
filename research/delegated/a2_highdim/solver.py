"""Budgeted exact-adjoint L-BFGS-B driver for the a2_highdim comparison.

Every method uses the same low-to-high cumulative frequency schedule and the
same frozen L-BFGS-B settings/boxes.  Pose coordinates are optimised in
scaled coordinates z=x*POSE_SCALE and the gradient is transformed with
grad_z = grad_x / POSE_SCALE (the corrected chain rule; the legacy A2
solver used grad_x*POSE_SCALE and is deliberately not reused).

The last ACCEPTED iterate is preserved at any budget/time stop, never an
unaccepted trial.  Work is charged per full-wave RHS solve column
(forward states, material/pose solves and adjoints), with factorisations,
operator products and wall seconds reported separately.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.optimize import minimize

from common import (
    CostLedger,
    N_ALPHA,
    POSE_SCALE,
    jsonable,
    params_from_scaled,
    scaled_from_params,
    z_bounds,
)
from objective import BudgetExhausted, _fg_with_fixed_x

METHODS = ("coherent_fixedpose", "coherent_joint", "intensity_joint", "oracle")


@dataclass
class RunOutcome:
    method: str
    seed: int
    radius_label: str
    alpha: np.ndarray
    x: np.ndarray
    ledger: CostLedger
    status: str = ""
    stages_visited: list[int] = field(default_factory=list)
    stage_events: list[dict[str, Any]] = field(default_factory=list)
    task_errors: list[str] = field(default_factory=list)
    is_oracle: bool = False
    all_stages_visited: bool = False
    wall_seconds: float = 0.0


def _update_if_accepted(
    tracker: dict[str, Any], new_z: np.ndarray, success: bool
) -> None:
    """Callback updates only accepted major iterates."""
    tracker["current"] = np.asarray(new_z, dtype=float).copy()
    tracker["callback_calls"] += 1


def run_method(
    model: Any,
    y: np.ndarray,
    sigma: float,
    method: str,
    x0: np.ndarray,
    settings: dict[str, Any],
    seed: int,
    radius_label: str,
    x_true: np.ndarray,
    wall_cap_seconds: float = 900.0,
) -> RunOutcome:
    """Run one method under the frozen schedule/budget; returns the last
    accepted estimate and full work ledger."""
    if method not in METHODS:
        raise ValueError(f"unknown method {method!r}")
    if method == "oracle":
        x_fixed = np.asarray(x_true, dtype=float).copy()
    elif method in ("coherent_fixedpose",):
        x_fixed = np.asarray(x0, dtype=float).copy()
    else:
        x_fixed = None

    is_oracle = method == "oracle"
    free_x = x_fixed is None
    alpha0 = np.full(N_ALPHA, float(settings["alpha_init"]))
    x_cur = np.asarray(x0, dtype=float).copy()
    if free_x:
        z_cur = scaled_from_params(alpha0, x_cur)
    else:
        z_cur = alpha0.copy()
    accepted = z_cur.copy()
    tracker = {"current": accepted.copy(), "callback_calls": 0}
    ledger = CostLedger()
    budget = int(settings["budget_rhs_columns"])
    maxls = int(settings["maxls"])
    ftol = float(settings["ftol"])
    gtol = float(settings["gtol"])
    stage_maxiter = [int(v) for v in settings["stage_maxiter"]]
    stage_freqs = [tuple(int(v) for v in s) for s in settings["stage_freqs"]]

    outcome = RunOutcome(
        method=method,
        seed=int(seed),
        radius_label=str(radius_label),
        alpha=alpha0.copy(),
        x=x_cur.copy(),
        ledger=ledger,
        is_oracle=is_oracle,
    )
    t_start = time.perf_counter()

    if free_x:
        bounds = z_bounds()
    else:
        from common import ALPHA_BOX

        bounds = [tuple(ALPHA_BOX)] * N_ALPHA

    for si, freq_ids in enumerate(stage_freqs):
        if len(stage_maxiter) <= si:
            outcome.task_errors.append(
                "stage_maxiter shorter than stage_freqs in frozen settings"
            )
            break
        if ledger.units >= budget:
            break
        z0 = accepted.copy()
        cache: dict[bytes, tuple[float, np.ndarray]] = {}

        # A budgeted gradient objective in scaled z coordinates.
        fg, _ = _fg_with_fixed_x(
            model,
            y,
            sigma,
            freq_ids,
            method,
            ledger,
            budget,
            x_fixed,
        )

        # Memoising wrapper (the physics core is deterministic; identical
        # repeated line-search points must not be double charged).
        def fg_memo(z: np.ndarray) -> tuple[float, np.ndarray]:
            z = np.asarray(z, dtype=float).reshape(-1)
            key = bytes(z.tobytes())
            if key in cache:
                return cache[key]
            val = fg(z)
            cache[key] = val
            return val

        def cb(zk: np.ndarray) -> None:
            _update_if_accepted(tracker, zk, True)

        event: dict[str, Any] = {
            "stage": si,
            "freq_ids": list(freq_ids),
            "units_before": ledger.units,
            "wall_before": time.perf_counter() - t_start,
        }
        try:
            res = minimize(
                fg_memo,
                z0,
                method="L-BFGS-B",
                jac=True,
                bounds=bounds,
                options={
                    "maxiter": stage_maxiter[si],
                    "maxls": maxls,
                    "ftol": ftol,
                    "gtol": gtol,
                    "maxfun": 100000,
                },
                callback=cb,
            )
            event.update(
                {
                    "res_success": bool(getattr(res, "success", False)),
                    "res_message": str(getattr(res, "message", "")),
                    "nit": int(getattr(res, "nit", -1)),
                    "nfev": int(getattr(res, "nfev", -1)),
                }
            )
            if bool(getattr(res, "success", False)) and np.all(
                np.isfinite(res.x)
            ):
                accepted = np.asarray(res.x, dtype=float).copy()
            elif tracker["callback_calls"] > 0:
                accepted = tracker["current"].copy()
        except BudgetExhausted as exc:
            event["budget_exhausted"] = str(exc)
            accepted = tracker["current"].copy()
            outcome.status = "budget_stopped"
            outcome.stage_events.append(event)
            outcome.stages_visited.append(si)
            break
        except Exception as exc:  # task error, never fabricate convergence
            outcome.task_errors.append(f"stage {si}: {type(exc).__name__}: {exc}")
            outcome.status = "task_error"
            accepted = tracker["current"].copy()
            outcome.stage_events.append(event)
            outcome.stages_visited.append(si)
            break
        outcome.stages_visited.append(si)
        outcome.stage_events.append(event)
        if time.perf_counter() - t_start > wall_cap_seconds:
            outcome.status = "time_stopped"
            break

    if not outcome.status:
        outcome.status = "schedule_complete"
    if free_x:
        alpha_end, x_end = params_from_scaled(accepted)
    else:
        alpha_end = accepted.copy()
        x_end = x_cur if x_fixed is None else x_fixed.copy()
    outcome.alpha = alpha_end.copy()
    outcome.x = x_end.copy()
    outcome.all_stages_visited = (
        len(outcome.stages_visited) >= len(stage_freqs)
        and outcome.status not in ("task_error",)
    )
    outcome.wall_seconds = time.perf_counter() - t_start
    if outcome.status == "task_error":
        outcome.task_errors.append(
            "task errors reported separately; no primary success is claimed"
        )
    return outcome
