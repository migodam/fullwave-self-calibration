"""Physical objectives: exact coherent Gaussian and matched intensity NLL.

The phaseless objective uses the exact induced noncentral-chi-square /
Rice intensity law of the proper complex parent and an adjoint gradient
chain to the total field.  No additive Gaussian intensity surrogate is
used.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import i0e, i1e

from common import SQRT2, CostLedger, select_rows


def coherent_loss(res_complex: np.ndarray, sigma: float) -> float:
    return 0.5 * float(np.sum(np.abs(realify(res_complex, sigma)) ** 2))


def realify(z: np.ndarray, sigma: float) -> np.ndarray:
    return (SQRT2 / sigma) * np.concatenate([np.real(z), np.imag(z)])


def gaussian_gradient(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
    prior_terms=None,
):
    """Exact physical Gaussian objective + analytic adjoint gradient."""
    out = model.adjoint_gradient(
        alpha, x, select_rows(y, freq_ids), float(sigma), list(freq_ids)
    )
    ledger.charge_model_work(out["work"])
    loss = float(out["loss"])
    g_a = np.asarray(out["grad_alpha"], dtype=float).copy()
    g_x = np.asarray(out["grad_x"], dtype=float).copy()
    if prior_terms is not None:
        loss += float(prior_terms[0])
        g_a = g_a + prior_terms[1]
    return {
        "loss": loss,
        "grad_alpha": g_a,
        "grad_x": g_x,
        "work": out["work"],
        "residual": None,
    }


def forward_predict(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
    prior_terms=None,
):
    fw = model.forward(alpha, x, list(freq_ids), jacobian=False)
    ledger.charge_model_work(fw["work"])
    total = fw["total"]
    res = total - select_rows(y, freq_ids)
    loss = coherent_loss(res, sigma)
    if prior_terms is not None:
        loss += float(prior_terms[0])
    return {
        "loss": loss,
        "total": total,
        "residual_complex": res,
        "work": fw["work"],
    }


def intensity_nll_and_weights(total: np.ndarray, y: np.ndarray, sigma: float):
    """Matched intensity NLL and real gradient weights per complex row.

    y has the same complex rows as total.  For each row the observed parent
    intensity is r=|y|^2 and the predicted mean is mu=total.  The NLL is
    (sqrt(r)-|mu|)^2/sigma^2 - log i0e(2|mu|sqrt(r)/sigma^2) + log(sigma^2).
    weights b = g_r + i g_i satisfy grad_theta = Re(J^H b).
    """
    y = np.asarray(y, dtype=np.complex128).reshape(-1)
    total = np.asarray(total, dtype=np.complex128).reshape(-1)
    sigma = float(sigma)
    r = np.abs(y) ** 2
    sr = np.sqrt(np.maximum(r, 0.0))
    mu = total
    L = np.abs(mu)
    small = L <= 1e-30
    Ls = np.maximum(L, 1e-30)
    x = 2.0 * Ls * sr / sigma**2
    e0 = i0e(x)
    e1 = i1e(x)
    ratio = np.divide(e1, e0, out=np.zeros_like(e0), where=e0 > 0)
    nll = (sr - L) ** 2 / sigma**2 - np.log(np.maximum(e0, np.finfo(float).tiny)) + np.log(
        sigma**2
    )
    scale = (2.0 / sigma**2) * (L - sr * ratio)
    unit = np.divide(mu, Ls, out=np.zeros_like(mu), where=~small)
    weights = scale * unit
    weights[small] = 0.0
    return float(np.sum(nll)), weights


def phaseless_loss(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
    prior_terms=None,
):
    fw = forward_predict(model, alpha, x, y, sigma, freq_ids, ledger)
    nll, _ = intensity_nll_and_weights(fw["total"], select_rows(y, freq_ids), sigma)
    if prior_terms is not None:
        nll += float(prior_terms[0])
    return {"loss": nll, "total": fw["total"], "work": fw["work"]}


def phaseless_gradient(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
    prior_terms=None,
):
    """Adjoint gradient of the intensity NLL against returned state matrices."""
    fw = model.forward(alpha, x, list(freq_ids), jacobian=False)
    ledger.charge_model_work(fw["work"])
    y_sel = select_rows(y, freq_ids)
    nll, weights = intensity_nll_and_weights(fw["total"], y_sel, sigma)

    states = fw["states"]
    # The physical core does not expose its internal LU for arbitrary
    # adjoint right-hand sides, so the solver re-factors and charges the
    # factorisations honestly.
    lu_by_freq: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for st in states:
        fi = int(st["freq_idx"])
        if fi not in lu_by_freq:
            lu_by_freq[fi] = lu_factor(st["M"], check_finite=False)
            ledger.charge_factorizations(1)

    g_a = np.zeros(model.n_alpha, dtype=np.complex128)
    g_x = np.zeros(3, dtype=np.complex128)
    for st in states:
        lu_piv = lu_by_freq[int(st["freq_idx"])]
        w = weights[st["row_start"] : st["row_stop"]]
        S_mat = st["S"]
        v = S_mat.conj().T @ w
        ledger.charge_operator_products(1)
        lam = lu_solve(lu_piv, v, trans=2, check_finite=False)
        ledger.charge_rhs(1, "adjoint")
        Tmat = st["Tmat"]
        Etot = st["E_total"]
        g_a += Tmat.conj().T @ (np.conj(Etot) * lam)
        ledger.charge_operator_products(1)
        C = st["direct_derivative_dx"]
        b_x = np.einsum("la,a->l", C.conj(), w)
        ledger.charge_operator_products(1)
        j = st["j"]
        for ell in range(3):
            u = st["dS_dx"][ell].conj().T @ w
            b_x[ell] += np.conj(j) @ u
            ledger.charge_operator_products(1)
        b_x += np.conj(st["db_dx"]).T @ lam
        ledger.charge_operator_products(1)
        g_x += b_x

    grad_alpha = np.real(g_a)
    grad_x = np.real(g_x)
    if prior_terms is not None:
        nll += float(prior_terms[0])
        grad_alpha = grad_alpha + prior_terms[1]
    return {
        "loss": nll,
        "grad_alpha": grad_alpha,
        "grad_x": grad_x,
        "work": fw["work"],
    }


def prior_terms(
    alpha: np.ndarray, settings: dict[str, Any]
) -> tuple[float, np.ndarray] | None:
    strength = float(settings.get("prior_strength", 0.0))
    if strength <= 0.0:
        return None
    center = np.asarray(settings.get("prior_center", [0.5] * 9), dtype=float)
    d = np.asarray(alpha, dtype=float) - center
    return 0.5 * strength * float(np.dot(d, d)), strength * d
