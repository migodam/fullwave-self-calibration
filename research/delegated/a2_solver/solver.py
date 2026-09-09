"""E4 numerical solver driver (direct, prasc, fixed_rank, phaseless,
direct_control) under one charged-work ledger and stage budget policy."""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from scipy.optimize import minimize

from common import (
    BudgetExhausted,
    CostLedger,
    FREQ_STAGES,
    PRASC_RANK_LADDER,
    SQRT2,
    STAGE_FRACTIONS,
    default_settings,
    pose_scale,
    select_rows,
    settings_hash,
)
from objective import (
    coherent_loss,
    forward_predict,
    gaussian_gradient,
    phaseless_gradient,
    phaseless_loss,
    prior_terms,
    realify,
)
from reduced import build_sensing_basis, reduced_evaluate


METHODS = ("direct", "prasc", "fixed_rank", "phaseless", "direct_control")


class StageLimit(Exception):
    """Current cumulative frequency-stage budget is spent (not terminal)."""


def _scaled_from_params(alpha: np.ndarray, x: np.ndarray, lam: float) -> np.ndarray:
    return np.concatenate(
        [np.asarray(alpha, dtype=float), np.asarray(x, dtype=float) * pose_scale(lam)]
    )


def _params_from_scaled(z: np.ndarray, lam: float) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    return z[:9].copy(), z[9:] / pose_scale(lam)


def _bounds(lam: float):
    sx = pose_scale(lam)
    return [(0.03, 2.5)] * 9 + [(-0.7 * float(s), 0.7 * float(s)) for s in sx]


class Tracker:
    def __init__(self, z0: np.ndarray) -> None:
        self.current = np.asarray(z0, dtype=float).copy()
        self.loss_cache: dict[bytes, float] = {}

    def callback(self, zk: np.ndarray) -> None:
        self.current = np.asarray(zk, dtype=float).copy()

    def remember(self, z: np.ndarray, loss: float) -> None:
        self.loss_cache[bytes(np.asarray(z, dtype=float).tobytes())] = float(loss)

    def cached(self, z: np.ndarray) -> float | None:
        return self.loss_cache.get(bytes(np.asarray(z, dtype=float).tobytes()))


@dataclass
class RunState:
    method: str
    alpha: np.ndarray
    x: np.ndarray
    ledger: CostLedger
    settings: dict[str, Any]
    budget: int
    lam: float = 0.5
    mode: str = "exact"  # exact | reduced
    rank: int = 0
    U: np.ndarray | None = None
    basis_x: np.ndarray | None = None
    incomplete_stages: list[int] = field(default_factory=list)
    rank_trajectory: list[dict[str, Any]] = field(default_factory=list)
    cert_history: list[dict[str, Any]] = field(default_factory=list)
    fallback_events: list[dict[str, Any]] = field(default_factory=list)
    exact_checks: int = 0
    last_exact_loss: float | None = None
    last_exact_alpha: np.ndarray | None = None
    last_exact_x: np.ndarray | None = None
    reduced_moves: int = 0
    failure: str | None = None

    def z(self) -> np.ndarray:
        return _scaled_from_params(self.alpha, self.x, self.lam)

    def store_exact(self, loss: float) -> None:
        self.last_exact_loss = float(loss)
        self.last_exact_alpha = self.alpha.copy()
        self.last_exact_x = self.x.copy()

    def restore_exact(self) -> None:
        if self.last_exact_alpha is not None:
            self.alpha = self.last_exact_alpha.copy()
            self.x = self.last_exact_x.copy()


class StageGuard:
    """Raises before an evaluation that would exceed the current cap."""

    def __init__(self, run: RunState, stage_idx: int) -> None:
        self.run = run
        self.stage_idx = stage_idx
        self.total_cap = STAGE_FRACTIONS[-1] * float(run.budget)
        self.cap = STAGE_FRACTIONS[stage_idx] * float(run.budget)

    def _check_value(self, projected: float) -> None:
        if projected > self.total_cap + 1e-9:
            raise BudgetExhausted(
                f"total work cap {self.total_cap:.3f} would be exceeded "
                f"({projected:.3f})"
            )
        if projected > self.cap + 1e-9:
            raise StageLimit(
                f"stage {self.stage_idx} cumulative cap {self.cap:.3f} "
                f"would be exceeded ({projected:.3f})"
            )

    def check(self, estimated_units: float, what: str) -> None:
        self._check_value(self.run.ledger.units + float(estimated_units))

    def raise_if_over(self, what: str) -> None:
        self._check_value(self.run.ledger.units)


def make_exact_fg(
    model: Any,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    run: RunState,
    guard: StageGuard,
    phaseless: bool,
    tracker: Tracker | None,
):
    lam = run.lam
    scale = pose_scale(lam)
    f_est = 12.0 * len(freq_ids)

    def fg(z: np.ndarray):
        z = np.asarray(z, dtype=float)
        guard.check(f_est, "exact combined evaluation")
        alpha, x = _params_from_scaled(z, lam)
        pt = prior_terms(alpha, run.settings)
        if phaseless:
            out = phaseless_gradient(
                model, alpha, x, y, sigma, freq_ids, run.ledger, pt
            )
        else:
            out = gaussian_gradient(
                model, alpha, x, y, sigma, freq_ids, run.ledger, pt
            )
        guard.raise_if_over("exact combined evaluation")
        if tracker is not None:
            tracker.remember(z, float(out["loss"]))
        g = np.concatenate(
            [out["grad_alpha"], np.asarray(out["grad_x"]) * scale]
        )
        return float(out["loss"]), g

    return fg


def make_exact_f_only(
    model: Any,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    run: RunState,
    guard: StageGuard,
    phaseless: bool,
):
    lam = run.lam

    def f(z: np.ndarray):
        z = np.asarray(z, dtype=float)
        guard.check(6.0 * len(freq_ids), "exact objective-only evaluation")
        alpha, x = _params_from_scaled(z, lam)
        pt = prior_terms(alpha, run.settings)
        if phaseless:
            out = phaseless_loss(
                model, alpha, x, y, sigma, freq_ids, run.ledger, pt
            )
        else:
            out = forward_predict(
                model, alpha, x, y, sigma, freq_ids, run.ledger, pt
            )
        guard.raise_if_over("exact objective-only evaluation")
        return float(out["loss"])

    return f


def _run_lbfgsb_step(
    fg: Callable[[np.ndarray], tuple[float, np.ndarray]],
    z0: np.ndarray,
    bounds: list[tuple[float, float]],
    settings: dict[str, Any],
    maxiter: int,
    tracker: Tracker,
) -> Any:
    return minimize(
        fg,
        z0,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={
            "maxiter": int(maxiter),
            "maxls": int(settings.get("maxls", 10)),
            "ftol": float(settings.get("ftol", 1e-8)),
            "gtol": float(settings.get("gtol", 1e-6)),
            "maxfun": 100000,
        },
        callback=tracker.callback,
    )


def exact_checkpoint(
    model: Any,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    run: RunState,
    guard: StageGuard,
    phaseless: bool,
) -> float:
    """Exact physical objective at the current iterate (charged)."""
    f = make_exact_f_only(
        model, y, sigma, freq_ids, run, guard, phaseless
    )
    loss = f(run.z())
    run.exact_checks += 1
    run.store_exact(loss)
    return float(loss)


def run_exact_stage(
    model: Any,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    stage_idx: int,
    run: RunState,
    settings: dict[str, Any],
    phaseless: bool = False,
) -> None:
    """L-BFGS-B exact stage.  control_rounds=True restarts curvature at every
    accepted iterate (the prasc/direct_control chart cadence)."""
    guard = StageGuard(run, stage_idx)
    lam = run.lam
    bounds = _bounds(lam)
    tracker = Tracker(run.z())
    control = bool(settings.get("control_rounds", True)) and (
        run.method in ("direct_control", "prasc")
    )
    while True:
        z0 = run.z()
        fg = make_exact_fg(
            model, y, sigma, freq_ids, run, guard, phaseless, tracker
        )
        try:
            res = _run_lbfgsb_step(
                fg, z0, bounds, settings, 1 if control else 3000, tracker
            )
        except StageLimit:
            run.alpha, run.x = _params_from_scaled(tracker.current, lam)
            run.incomplete_stages.append(stage_idx)
            return
        except BudgetExhausted:
            run.alpha, run.x = _params_from_scaled(tracker.current, lam)
            raise
        z_new = np.asarray(res.x, dtype=float)
        moved = float(np.linalg.norm(z_new - z0)) > 1e-12 * max(
            1.0, float(np.linalg.norm(z0))
        )
        loss_new = tracker.cached(z_new)
        run.alpha, run.x = _params_from_scaled(z_new, lam)
        if loss_new is not None:
            run.store_exact(loss_new)
        if not moved or int(getattr(res, "nit", 0)) == 0:
            break
        if not control and getattr(res, "success", False):
            break


def _reduced_objective_gradient(
    model: Any,
    y_sel: np.ndarray,
    sigma: float,
    z: np.ndarray,
    out: dict[str, Any],
    lam: float,
    settings: dict[str, Any],
) -> tuple[float, np.ndarray]:
    alpha, _ = _params_from_scaled(z, lam)
    pt = prior_terms(alpha, settings)
    total = out["total"]
    loss = coherent_loss(total - y_sel, sigma)
    A = out["A"]
    B = out["B"]
    AB = np.hstack([A, B])
    res_re = realify(total - y_sel, sigma)
    JR = (SQRT2 / sigma) * np.vstack([np.real(AB), np.imag(AB)])
    g_phys = JR.T @ res_re
    if pt is not None:
        loss += float(pt[0])
        g_phys[:9] += pt[1]
    scale = pose_scale(lam)
    g = np.concatenate([g_phys[:9], g_phys[9:] * scale])
    return float(loss), g


def build_chart(
    model: Any,
    run: RunState,
    freq_ids: tuple[int, ...],
    rank: int,
    settings: dict[str, Any],
) -> dict[str, Any]:
    basis = build_sensing_basis(
        model,
        run.x,
        freq_ids,
        rank,
        run.ledger,
        tsom_append=bool(settings.get("tsom_append", False)),
        tsom_target_fraction=float(settings.get("tsom_target_fraction", 0.0)),
    )
    run.U = basis["U"]
    run.basis_x = run.x.copy()
    run.rank = int(basis["actual_rank"])
    run.rank_trajectory.append(
        {
            "x": run.x.tolist(),
            "rank": run.rank,
            "ldet": basis["ldet"],
            "tsom_rank": basis["tsom_rank"],
            "sensing_target": rank,
        }
    )
    return basis


def _reduced_trial_units(
    ledger: CostLedger, freq_ids: tuple[int, ...], rank: int
) -> float:
    cal = ledger.calibration
    if cal is None:
        return 0.0
    F = len(freq_ids)
    rr = cal.reduced_solve_ratio(rank)
    qr = cal.reduced_qr_ratio(rank)
    return rr * 78.0 * F + qr * F


def run_reduced_stage(
    model: Any,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    stage_idx: int,
    run: RunState,
    settings: dict[str, Any],
) -> None:
    guard = StageGuard(run, stage_idx)
    lam = run.lam
    y_sel = select_rows(y, freq_ids)
    fixed = run.method == "fixed_rank"
    if fixed:
        run.rank = int(settings.get("fixed_rank") or 8)
    if run.rank <= 0:
        run.rank = 8
    ladder = list(PRASC_RANK_LADDER)
    max_n = model.n_cells

    while True:
        if run.mode == "exact":
            run_exact_stage(
                model, y, sigma, freq_ids, stage_idx, run, settings, False
            )
            return
        chart_rank = min(int(run.rank), max_n)
        if chart_rank >= max_n:
            run.mode = "exact"
            run.fallback_events.append(
                {
                    "stage": stage_idx,
                    "reason": "rank_reached_n_direct",
                    "rank": chart_rank,
                }
            )
            continue
        build_chart(model, run, freq_ids, chart_rank, settings)
        r_actual = int(run.U.shape[1]) if run.U is not None else chart_rank
        est = _reduced_trial_units(run.ledger, freq_ids, r_actual)
        guard.check(est, "reduced chart evaluation")
        out0 = reduced_evaluate(
            model,
            run.alpha,
            run.x,
            run.U,
            freq_ids,
            run.ledger,
            need_jacobian=True,
        )
        guard.raise_if_over("reduced chart evaluation")
        z0 = run.z()
        loss0, grad0 = _reduced_objective_gradient(
            model, y_sel, sigma, z0, out0, lam, settings
        )
        cert = None
        if not fixed:
            from cert import reduced_certification

            cert = reduced_certification(
                model, run.alpha, y_sel, sigma, freq_ids, out0, settings
            )
        run.cert_history.append(
            {
                "stage": stage_idx,
                "rank": r_actual,
                "certification": cert,
                "loss": loss0,
            }
        )
        if cert is not None and not cert["pass"]:
            nxt = next((r for r in ladder if r > chart_rank), None)
            if nxt is not None and nxt < max_n:
                nxt_est = _reduced_trial_units(run.ledger, freq_ids, nxt)
                if guard.total_cap < run.ledger.units + nxt_est + 1e-9:
                    # cannot afford a higher certified rank under the total
                    # work cap; conservative certificate forces exact fallback
                    run.mode = "exact"
                    run.fallback_events.append(
                        {
                            "stage": stage_idx,
                            "reason": "certificate_rejected_no_rank_budget",
                            "rank_from": chart_rank,
                            "rank_to": nxt,
                        }
                    )
                    continue
                run.rank = nxt
                run.fallback_events.append(
                    {
                        "stage": stage_idx,
                        "reason": "certificate_rejected_promote",
                        "rank_from": chart_rank,
                        "rank_to": nxt,
                    }
                )
                continue
            if settings.get("exact_fallback", True):
                run.mode = "exact"
                run.fallback_events.append(
                    {
                        "stage": stage_idx,
                        "reason": "certificate_unavailable_or_failed",
                        "rank": chart_rank,
                        "cert": cert,
                    }
                )
                continue
            run.failure = "no_certified_reduced_state"
            return

        # One L-BFGS-B iteration with the basis frozen through the line
        # search.  A fresh call resets quasi-Newton curvature after the
        # accepted chart event.
        bounds = _bounds(lam)
        tracker = Tracker(z0)

        def fg(z: np.ndarray):
            z = np.asarray(z, dtype=float)
            if np.array_equal(z, z0):
                return loss0, grad0
            guard.check(
                _reduced_trial_units(run.ledger, freq_ids, r_actual),
                "reduced trial evaluation",
            )
            alpha, x = _params_from_scaled(z, lam)
            out = reduced_evaluate(
                model,
                alpha,
                x,
                run.U,
                freq_ids,
                run.ledger,
                need_jacobian=True,
            )
            guard.raise_if_over("reduced trial evaluation")
            loss, g = _reduced_objective_gradient(
                model, y_sel, sigma, z, out, lam, settings
            )
            tracker.remember(z, loss)
            return loss, g

        try:
            res = _run_lbfgsb_step(fg, z0, bounds, settings, 1, tracker)
        except StageLimit:
            run.alpha, run.x = _params_from_scaled(tracker.current, lam)
            run.incomplete_stages.append(stage_idx)
            return
        except BudgetExhausted:
            run.alpha, run.x = _params_from_scaled(tracker.current, lam)
            raise
        z_new = np.asarray(res.x, dtype=float)
        moved = float(np.linalg.norm(z_new - z0)) > 1e-12 * max(
            1.0, float(np.linalg.norm(z0))
        )
        loss_new = tracker.cached(z_new)
        if loss_new is None or not moved:
            if fixed:
                return
            # reduced stage converged: exact physical checkpoint
            prev_loss = run.last_exact_loss
            prev_alpha = run.alpha.copy() if run.last_exact_alpha is not None else None
            prev_x = run.x.copy() if run.last_exact_x is not None else None
            try:
                exact_checkpoint(
                    model, y, sigma, freq_ids, run, guard, False
                )
            except StageLimit:
                run.incomplete_stages.append(stage_idx)
                return
            if (
                run.reduced_moves
                and prev_loss is not None
                and prev_alpha is not None
                and run.last_exact_loss is not None
                and run.last_exact_loss > prev_loss + 1e-12
            ):
                # The certified reduced path did not improve the exact
                # physical objective: reject the candidate and use exact
                # solves for the remainder of this stage.
                run.fallback_events.append(
                    {
                        "stage": stage_idx,
                        "reason": "reduced_exact_checkpoint_rejected",
                        "loss_before": prev_loss,
                        "loss_after": run.last_exact_loss,
                    }
                )
                run.alpha = prev_alpha
                run.x = prev_x
                run.last_exact_loss = prev_loss
                run.last_exact_alpha = prev_alpha
                run.last_exact_x = prev_x
                run.mode = "exact"
            return
        alpha_new, x_new = _params_from_scaled(z_new, lam)
        run.alpha = alpha_new
        run.x = x_new
        run.reduced_moves += 1


def solve(
    model: Any,
    y: np.ndarray,
    sigma: float,
    alpha0: np.ndarray,
    x0: np.ndarray,
    method: str,
    budget: int,
    settings: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Callable E4 solver API (see TASK.md artifacts/API)."""
    method = str(method)
    if method not in METHODS:
        raise ValueError(f"unknown method {method!r}")
    if budget not in (200, 800):
        raise ValueError(f"budget must be 200 or 800, got {budget!r}")
    merged = default_settings()
    if settings:
        merged.update(settings)
    settings = merged
    if method == "fixed_rank" and settings.get("fixed_rank") is None:
        raise ValueError("fixed_rank requires settings['fixed_rank']")

    ledger = CostLedger()
    if hasattr(model, "_calibration") and getattr(model, "_calibration", None) is not None:
        ledger.calibration = model._calibration
    run = RunState(
        method=method,
        alpha=np.asarray(alpha0, dtype=float).copy(),
        x=np.asarray(x0, dtype=float).copy(),
        ledger=ledger,
        settings=settings,
        budget=int(budget),
        mode="reduced" if method in ("prasc", "fixed_rank") else "exact",
    )
    init_units = int(round(float(settings.get("init_units", 0.0))))
    if init_units:
        ledger.charge_rhs(init_units, "initialisation")

    t_start = time.perf_counter()
    try:
        for stage_idx, freq_ids in enumerate(FREQ_STAGES):
            try:
                if method in ("prasc", "fixed_rank"):
                    run_reduced_stage(
                        model, y, sigma, freq_ids, stage_idx, run, settings
                    )
                else:
                    run_exact_stage(
                        model,
                        y,
                        sigma,
                        freq_ids,
                        stage_idx,
                        run,
                        settings,
                        method == "phaseless",
                    )
            except StageLimit:
                run.incomplete_stages.append(stage_idx)
                continue
    except BudgetExhausted as exc:
        run.failure = f"budget_exhausted: {exc}"
    ledger.add_wall(time.perf_counter() - t_start)

    return {
        "method": method,
        "alpha_est": run.alpha.tolist(),
        "x_est": run.x.tolist(),
        "failure": run.failure,
        "status": (
            run.failure
            if run.failure
            else (
                "incomplete_stage_budget"
                if run.incomplete_stages
                else "completed_within_budget"
            )
        ),
        "mode": run.mode,
        "ledger": ledger.snapshot(),
        "settings": settings,
        "config_hash": settings_hash(settings),
        "rank": run.rank,
        "rank_trajectory": run.rank_trajectory,
        "cert_history": run.cert_history,
        "fallback_events": run.fallback_events,
        "incomplete_stages": run.incomplete_stages,
        "exact_checks": run.exact_checks,
        "reduced_moves": run.reduced_moves,
    }
