"""Conservative acceptance certificates and diagnostics for reduced states."""
from __future__ import annotations

from typing import Any

import numpy as np

from common import SQRT2


def _realify(z: np.ndarray, sigma: float) -> np.ndarray:
    return (SQRT2 / sigma) * np.concatenate([np.real(z), np.imag(z)])


def data_error_bound(
    model: Any, alpha: np.ndarray, reduced_out: dict[str, Any], sigma: float
) -> tuple[float, dict[str, Any]]:
    """Conservative realified data error bound of the reduced prediction."""
    from certificates import state_bound
    from reduced import passive_certificate

    total_bound = 0.0
    per_state: list[float] = []
    per_freq: dict[int, float] = {}
    for rec in reduced_out["states"]:
        fi = int(rec["freq_idx"])
        cert = passive_certificate(model, alpha, fi)
        b = state_bound(cert, rec["state_residual"], rec["S"], float(sigma))
        total_bound += b*b
        per_state.append(float(b))
        per_freq[fi] = per_freq.get(fi, 0.0) + float(b)
    realified = SQRT2 * np.sqrt(total_bound)
    return float(realified), {
        "per_state_complex": per_state,
        "per_freq_complex": {str(k): v for k, v in per_freq.items()},
        "available": all(
            passive_certificate(model, alpha, int(rec["freq_idx"])).available
            for rec in reduced_out["states"]
        ),
    }


def gradient_error_summary(
    model: Any,
    alpha: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    reduced_out: dict[str, Any],
) -> dict[str, Any]:
    """Certified bound on ||g_exact - g_reduced|| for the physical objective."""
    from certificates import derivative_bound, gradient_error_bound
    from reduced import domain_operator, passive_certificate

    q = model.n_alpha
    total = reduced_out["total"]
    res_re = _realify(total - y, sigma)
    A = reduced_out.get("A")
    B = reduced_out.get("B")
    if A is None or B is None:
        raise ValueError("gradient_error_summary requires reduced A/B blocks")
    AB = np.hstack([A, B])
    Jreal = (SQRT2 / sigma) * np.vstack([np.real(AB), np.imag(AB)])
    gred = Jreal.T @ res_re
    # per-column complex data derivative error bounds
    col_bounds = np.zeros(q + 3, dtype=float)
    per_state = []
    for rec in reduced_out["states"]:
        fi = int(rec["freq_idx"])
        cert = passive_certificate(model, alpha, fi)
        D = domain_operator(model, fi)
        T = rec["Tmat"]
        z0 = rec["state_residual"]
        S = rec["S"]
        der = rec["derivative_residual"]
        row = np.zeros(q + 3)
        for a in range(q):
            # M_a = -diag(T_a) D, b_a = diag(T_a) E_inc
            Mv = -T[:, a][:, None] * D
            row[a] = derivative_bound(
                cert,
                z0,
                der[:, a],
                Mv,
                S,
                None,
                float(sigma),
            )
        for ell in range(3):
            Sv = rec["dS_dx"][ell]
            row[q + ell] = derivative_bound(
                cert, z0, der[:, q + ell], np.zeros_like(D), S, Sv, float(sigma)
            )
        col_bounds += row**2
        per_state.append(row)
    col_bounds_real = SQRT2 * np.sqrt(col_bounds)
    # In optimizer coordinates z=(alpha, x*scale), derivative columns scale
    # reciprocally. The error certificate must use that same norm.
    from common import pose_scale
    scaling = np.r_[np.ones(q),1/pose_scale(0.5)]
    col_bounds_real *= scaling
    Jreal *= scaling[None,:]
    gred = Jreal.T @ res_re
    jac_eb = float(np.linalg.norm(col_bounds_real))
    eb_y, _ = data_error_bound(model, alpha, reduced_out, sigma)
    bound = float(
        gradient_error_bound(
            res_re,
            Jreal,
            eb_y,
            jac_eb,
        )
    )
    gn = float(np.linalg.norm(gred))
    return {
        "gradient_error_bound": bound,
        "gradient_norm": gn,
        "relative_bound": float(bound / max(gn, 1e-15)),
        "jacobian_error_fro": jac_eb,
        "data_error_bound": eb_y,
        "column_bounds_real": col_bounds_real.tolist(),
    }


def reduced_certification(
    model: Any,
    alpha: np.ndarray,
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...],
    reduced_out: dict[str, Any],
    settings: dict[str, Any],
) -> dict[str, Any]:
    """Practical initial candidate acceptance with conservative bounds."""
    from reduced import passive_certificate

    total = reduced_out["total"]
    res_re = _realify(total - y, sigma)
    nr_real = int(res_re.size)
    eb = data_error_bound(model, alpha, reduced_out, sigma)
    eb_y = eb[0]
    grad = gradient_error_summary(
        model, alpha, y, sigma, freq_ids, reduced_out
    )
    avail = all(
        passive_certificate(model, alpha, int(rec["freq_idx"])).available
        for rec in reduced_out["states"]
    )
    data_ok = bool(
        avail
        and eb_y
        <= float(settings.get("cert_residual_factor", 0.05)) * np.sqrt(nr_real)
    )
    rel = grad["relative_bound"]
    grad_ok = bool(
        avail
        and rel <= float(settings.get("cert_gradient_rel", 0.1))
    )
    info = {
        "available": avail,
        "data_error_bound": eb_y,
        "data_ok": data_ok,
        "gradient_error_bound": grad["gradient_error_bound"],
        "gradient_norm": grad["gradient_norm"],
        "relative_gradient_bound": rel,
        "gradient_ok": grad_ok,
        "reduced_residual_realified": float(np.linalg.norm(res_re)),
        "threshold_real_measurements": float(np.sqrt(nr_real)),
        "pass": bool(data_ok and grad_ok),
    }
    return info


def objective_interval(
    model: Any,
    alpha: np.ndarray,
    y: np.ndarray,
    sigma: float,
    reduced_out: dict[str, Any],
    prior_value: float = 0.0,
) -> tuple[float, float]:
    from certificates import objective_interval

    res_re = _realify(reduced_out["total"] - y, sigma)
    r = float(np.linalg.norm(res_re))
    eb, _ = data_error_bound(model, alpha, reduced_out, sigma)
    lo, hi = objective_interval(res_re, eb, 0.0)
    return float(lo) + prior_value, float(hi) + prior_value


def covariance_diagnostic(A: np.ndarray, B: np.ndarray, sigma: float) -> dict[str, Any]:
    """Surrogate pose/map covariance gate. Never a coverage certificate."""
    if A is None or B is None:
        return {"available": False}
    AB = np.hstack([A, B])
    JR = (SQRT2 / sigma) * np.vstack([np.real(AB), np.imag(AB)])
    gram = np.real(JR.conj().T @ JR)
    try:
        cov = np.linalg.pinv(gram, rcond=1e-12)
        min_eig_pose = float(np.linalg.eigvalsh(cov[9:, 9:])[0])
        min_eig_joint = float(np.linalg.eigvalsh(cov)[0])
        rank = int(np.linalg.matrix_rank(gram, tol=1e-9 * np.max(gram)))
        ok = bool(min_eig_pose > 0 and np.isfinite(min_eig_pose) and rank == gram.shape[0])
    except np.linalg.LinAlgError:
        min_eig_pose = float("nan")
        min_eig_joint = float("nan")
        rank = 0
        ok = False
    return {
        "available": True,
        "surrogate_only": True,
        "pose_min_eigenvalue": min_eig_pose,
        "joint_min_eigenvalue": min_eig_joint,
        "gram_rank": rank,
        "surrogate_pass": ok,
    }
