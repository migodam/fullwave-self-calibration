"""Exact physical objectives for the a2_highdim comparison.

Coherent methods minimise the proper-complex Gaussian negative loglikelihood
0.5*||realify(total - y, sigma)||^2 through the physics core's analytic
adjoint gradient.  The matched-intensity joint method minimises the induced
Rice / noncentral-chi-square magnitude NLL of the same parent complex noise
with an exact adjoint gradient against the returned full-wave state
matrices.  No additive Gaussian intensity surrogate is used.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import i0e, i1e

from common import SQRT2, CostLedger, N_ALPHA, POSE_SCALE, params_from_scaled

ROWS_PER_FREQ = 72  # 3 poses x 2 illuminations x 12 receivers


def select_freq_rows(y: np.ndarray, freq_ids: tuple[int, ...]) -> np.ndarray:
    y = np.asarray(y, dtype=np.complex128).reshape(-1)
    ids = sorted(set(int(f) for f in freq_ids))
    parts = [y[f * ROWS_PER_FREQ : (f + 1) * ROWS_PER_FREQ] for f in ids]
    return np.concatenate(parts)


def realify(z: np.ndarray, sigma: float) -> np.ndarray:
    return (SQRT2 / float(sigma)) * np.concatenate([np.real(z), np.imag(z)])


def coherent_loss(res: np.ndarray, sigma: float) -> float:
    return 0.5 * float(np.sum(np.abs(realify(res, sigma)) ** 2))


# ---------------------------------------------------------------------------
# Coherent exact adjoint gradient
# ---------------------------------------------------------------------------


def coherent_value_gradient(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
) -> dict[str, Any]:
    y_sel = select_freq_rows(y, freq_ids)
    out = model.adjoint_gradient(
        alpha, x, y_sel, float(sigma), list(freq_ids)
    )
    ledger.charge_model_work(out["work"])
    return {
        "loss": float(out["loss"]),
        "grad_alpha": np.asarray(out["grad_alpha"], dtype=float).copy(),
        "grad_x": np.asarray(out["grad_x"], dtype=float).copy(),
        "work": out["work"],
    }


def coherent_forward(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
) -> dict[str, Any]:
    fw = model.forward(alpha, x, list(freq_ids), jacobian=False)
    ledger.charge_model_work(fw["work"])
    total = fw["total"]
    loss = coherent_loss(total - select_freq_rows(y, freq_ids), sigma)
    return {"loss": loss, "total": total, "work": fw["work"]}


# ---------------------------------------------------------------------------
# Induced Rice / noncentral-chi-square magnitude NLL
# ---------------------------------------------------------------------------


def rice_nll_weights(
    total: np.ndarray, y: np.ndarray, sigma: float
) -> tuple[float, np.ndarray]:
    """Rice magnitude NLL and complex gradient weights for one data row.

    The parent complex noise has per-real-component variance sigma^2/2
    (physics convention).  For observed magnitude rho=|y| and predicted mean
    mu=total, dropping the rho-independent constants:

        NLL = (rho - |mu|)^2/sigma^2 - log i0e(2 rho |mu| / sigma^2)
              + log(sigma^2/2).

    Returns the sum over rows and weights b with
    dL/dtheta = Re(J^H b), J = d total/dtheta.
    """
    total = np.asarray(total, dtype=np.complex128).reshape(-1)
    y = np.asarray(y, dtype=np.complex128).reshape(-1)
    sigma = float(sigma)
    rho = np.abs(y)
    mu_abs = np.abs(total)
    small = mu_abs <= 1e-30
    mu_safe = np.maximum(mu_abs, 1e-30)
    xarg = 2.0 * rho * mu_safe / sigma**2
    e0 = i0e(xarg)
    e1 = i1e(xarg)
    ratio = np.divide(e1, e0, out=np.zeros_like(e0), where=e0 > 0)
    nll = (rho - mu_abs) ** 2 / sigma**2 - np.log(
        np.maximum(e0, np.finfo(float).tiny)
    ) + np.log(sigma**2 / 2.0)
    scale = (2.0 / sigma**2) * (mu_abs - rho * ratio)
    unit = np.divide(total, mu_safe, out=np.zeros_like(total), where=~small)
    weights = scale * unit
    weights[small] = 0.0
    return float(np.sum(nll)), weights


def intensity_value_gradient(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
) -> dict[str, Any]:
    """Exact adjoint gradient of the induced Rice NLL against the returned
    physical state matrices.  The physics core does not expose its internal
    LU for arbitrary adjoint RHS, so the per-frequency LU factorisations are
    repeated and charged as factorisations (never hidden)."""
    fw = model.forward(alpha, x, list(freq_ids), jacobian=False)
    ledger.charge_model_work(fw["work"])
    y_sel = select_freq_rows(y, freq_ids)
    nll, weights = rice_nll_weights(fw["total"], y_sel, sigma)

    states = fw["states"]
    lu_by_freq: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for st in states:
        fi = int(st["freq_idx"])
        if fi not in lu_by_freq:
            lu_by_freq[fi] = lu_factor(st["M"], check_finite=False)
            ledger.charge(factorizations=1)

    grad_alpha = np.zeros(model.n_alpha, dtype=np.complex128)
    grad_x = np.zeros(3, dtype=np.complex128)
    for st in states:
        lu_piv = lu_by_freq[int(st["freq_idx"])]
        w = weights[st["row_start"] : st["row_stop"]]
        S_mat = st["S"]
        v = S_mat.conj().T @ w
        ledger.charge(operator_products=1)
        lam = lu_solve(lu_piv, v, trans=2, check_finite=False)
        ledger.charge(rhs=1, kind="adjoint")

        Tmat = st["Tmat"]
        Etot = st["E_total"]
        grad_alpha += Tmat.conj().T @ (np.conj(Etot) * lam)
        ledger.charge(operator_products=1)

        C = st["direct_derivative_dx"]  # (3, R)
        bx = np.einsum("la,a->l", C.conj(), w)
        ledger.charge(operator_products=1)
        j = st["j"]
        for ell in range(3):
            u = st["dS_dx"][ell].conj().T @ w
            bx[ell] += np.conj(j) @ u
            ledger.charge(operator_products=1)
        bx += np.conj(st["db_dx"]).T @ lam
        ledger.charge(operator_products=1)
        grad_x += bx

    return {
        "loss": float(nll),
        "grad_alpha": np.real(grad_alpha),
        "grad_x": np.real(grad_x),
        "work": fw["work"],
    }


def intensity_forward(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
) -> dict[str, Any]:
    fw = model.forward(alpha, x, list(freq_ids), jacobian=False)
    ledger.charge_model_work(fw["work"])
    y_sel = select_freq_rows(y, freq_ids)
    nll, _ = rice_nll_weights(fw["total"], y_sel, sigma)
    return {"loss": float(nll), "total": fw["total"], "work": fw["work"]}


# ---------------------------------------------------------------------------
# Scaled-coordinate L-BFGS-B objective wrappers
# ---------------------------------------------------------------------------


class BudgetExhausted(RuntimeError):
    pass


def _fg_with_fixed_x(
    model: Any,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    method: str,
    ledger: CostLedger,
    budget: int,
    x_fixed: np.ndarray,
) -> tuple[
    Callable[[np.ndarray], tuple[float, np.ndarray]],
    Callable[[np.ndarray], float],
]:
    rhs_est = 12 * len(freq_ids)
    cache: dict[bytes, tuple[float, np.ndarray]] = {}
    if method in ("coherent_fixedpose", "coherent_joint", "oracle"):
        grad_fn = coherent_value_gradient
    elif method in ("intensity_fixedpose", "intensity_joint"):
        grad_fn = intensity_value_gradient
    else:
        raise ValueError(f"unknown objective method {method!r}")

    def fg(z: np.ndarray) -> tuple[float, np.ndarray]:
        z = np.asarray(z, dtype=float).reshape(-1)
        key = bytes(z.tobytes())
        if key in cache:
            return cache[key]
        if ledger.units + rhs_est > int(budget):
            raise BudgetExhausted(
                f"RHS budget {budget} would be exceeded "
                f"({ledger.units}+{rhs_est})"
            )
        if x_fixed is not None:
            alpha = z.copy()
            x = x_fixed
        else:
            alpha, x = params_from_scaled(z)
        out = grad_fn(model, alpha, x, y, sigma, freq_ids, ledger)
        if x_fixed is not None:
            g = np.asarray(out["grad_alpha"], dtype=float).copy()
        else:
            g = np.concatenate(
                [
                    np.asarray(out["grad_alpha"], dtype=float),
                    np.asarray(out["grad_x"], dtype=float) / POSE_SCALE,
                ]
            )
        cache[key] = (float(out["loss"]), g)
        return cache[key]

    def f_only(z: np.ndarray) -> float:
        return fg(np.asarray(z, dtype=float))[0]

    return fg, f_only
