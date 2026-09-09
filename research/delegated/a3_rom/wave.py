"""Low-level full-wave fixed-chart ROM algebra for A3 development.

The physical operators are reused read-only from ``a2_physics/physics.py``.
No exact full-state solve is used while a chart is built or while the
fixed-chart proxy is differentiated.  The proxy state solves

    c = argmin_z || M U z - b ||_2^2          (C = M U),

use an economic QR of C.  The fixed-chart coefficient derivative keeps the
residual term

    (C*C) c_v = C*(b_v - M_v U c) + (M_v U)* r_s,   r_s = b - C c,

and is solved through the QR factors (two triangular solves per RHS column),
never through an explicit normal inverse.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.linalg import qr, solve_triangular

from common import A2, FREQ_IDS, ROWS_PER_FREQ, WorkLedger, realify, select_rows


def tau(model: Any, fi: int) -> float:
    omega = float(model.ks[int(fi)]) * A2.C0
    return float(model.config.loss_coefficient) / (omega * A2.EPS0)


def chi_and_Tmat(
    model: Any, alpha: np.ndarray, fi: int
) -> tuple[np.ndarray, np.ndarray]:
    loss = 1.0 + 1j * tau(model, fi)
    u = model.basis @ np.asarray(alpha, dtype=float)
    return u * loss, model.basis * loss


def state_geometry(
    model: Any, geom: dict[str, np.ndarray], fi: int, p: int, u: int
) -> dict[str, np.ndarray]:
    """Pose/illumination geometry (no current solve, no truth)."""
    from physics import (
        green_grad_first,
        green_grad_source,
        green_matrix,
    )

    k = float(model.ks[int(fi)])
    h2 = (model.h**2) * (k**2)
    rx_world = geom["rx_world"][p]
    tx_world = geom["tx_world"][p]
    drx_x = geom["drx_x"][p]
    dtx_x = geom["dtx_x"][p]

    S = h2 * green_matrix(model.points, rx_world, k)
    grad_rx = green_grad_first(model.points, rx_world, k)
    dS_dx = [
        h2 * np.einsum("asd,ad->as", grad_rx, drx_x[l]) for l in range(3)
    ]
    tx_u = tx_world[u]
    E_inc = green_matrix(tx_u[None, :], model.points, k)[:, 0]
    direct = green_matrix(tx_u[None, :], rx_world, k)[:, 0]
    grad_src = green_grad_source(model.points, tx_u, k)
    db_geom_dx = np.einsum("nd,ld->nl", grad_src, dtx_x[:, u, :])
    grad_at_rx = green_grad_first(tx_u[None, :], rx_world, k)[:, 0, :]
    rel_dx = drx_x - dtx_x[:, u, None, :]
    direct_dx = np.einsum("ad,lad->la", grad_at_rx, rel_dx)
    return {
        "S": S,
        "E_inc": E_inc,
        "direct": direct,
        "dS_dx": dS_dx,
        "db_geom_dx": db_geom_dx,
        "direct_derivative_dx": direct_dx,
        "rx_world": rx_world,
        "tx_world": tx_world,
    }


def domain_operator(model: Any, fi: int) -> np.ndarray:
    return model._domain_operator(int(fi))


def solve_normal_via_qr(
    Rtri: np.ndarray,
    rhs: np.ndarray,
    ledger: WorkLedger,
    n_rhs: int,
) -> np.ndarray:
    """Solve ``R^H R c = rhs`` without forming the normal inverse.

    Each of the ``n_rhs`` right-hand sides costs two triangular solves.
    """
    y = solve_triangular(
        Rtri.conj().T, rhs, lower=True, check_finite=False
    )
    c = solve_triangular(Rtri, y, lower=False, check_finite=False)
    ledger.charge_reduced(tri_rhs=2 * int(n_rhs))
    return c


def fixed_chart_evaluate(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    chart: dict[int, np.ndarray],
    y: np.ndarray,
    sigma: float,
    freq_ids: tuple[int, ...] = FREQ_IDS,
    jacobian: bool = True,
    ledger: WorkLedger | None = None,
) -> dict[str, Any]:
    """Reduced forward and fixed-chart Jacobian on frozen charts.

    ``chart`` maps each model frequency index to an orthonormal complex
    ``n x r`` current chart.  Returns approximate states/outputs plus, when
    requested, the exact derivative of the fixed LS proxy (including the
    residual term).  Work is appended to ``ledger``; no full M^{-1} state
    solve is performed here.
    """
    ledger = ledger or WorkLedger()
    alpha = np.asarray(alpha, dtype=float)
    x = np.asarray(x, dtype=float)
    n = model.n_cells
    geom = model._effective_geometry(x)
    q = model.n_alpha
    total_parts: list[np.ndarray] = []
    state_records: list[dict[str, Any]] = []
    row = 0
    ids = sorted(set(int(f) for f in freq_ids))
    y_sel = select_rows(y, ids)
    y_arr = y_sel
    per_freq: dict[int, dict[str, Any]] = {}

    for fi in ids:
        U = chart[int(fi)]
        r = int(U.shape[1])
        chi, Tmat = chi_and_Tmat(model, alpha, fi)
        D = domain_operator(model, fi)
        t0 = time.perf_counter()
        DU = D @ U
        ledger.charge_reduced(products=1, product_columns=int(r))
        C = U - chi[:, None] * DU
        t0 = time.perf_counter()
        Qred, Rtri = qr(C, mode="economic", check_finite=False)
        ledger.charge_reduced(qr=1, wall=time.perf_counter() - t0, dims=(n, r))
        QH = Qred.conj().T
        freq_rank_ok = True
        diag_abs = np.abs(np.diag(Rtri))
        if diag_abs.size:
            tol = 1e-11 * float(diag_abs[0]) if diag_abs[0] > 0 else 1e-15
            if float(np.min(diag_abs)) <= tol:
                freq_rank_ok = False
        per_freq[fi] = {"r": r, "rank_ok": bool(freq_rank_ok)}

        for p in range(model.n_pose):
            for u in range(model.n_tx):
                g = state_geometry(model, geom, fi, p, u)
                S = g["S"]
                b = chi * g["E_inc"]
                t0 = time.perf_counter()
                rhs_c = QH @ b
                c = solve_triangular(Rtri, rhs_c, lower=False, check_finite=False)
                ledger.charge_reduced(
                    tri_rhs=1,
                    products=1,
                    basis=1,
                    wall=time.perf_counter() - t0,
                )
                j_tilde = U @ c
                Cc = C @ c
                state_res = b - Cc
                Etot = g["E_inc"] + DU @ c
                pred = g["direct"] + S @ j_tilde
                ledger.charge_reduced(basis=3, products=1, product_columns=2)
                state_res_norm = float(np.linalg.norm(state_res))
                b_norm = float(np.linalg.norm(b))
                rec: dict[str, Any] = {
                    "freq_idx": fi,
                    "pose_idx": p,
                    "illum_idx": u,
                    "k": float(model.ks[fi]),
                    "U": U,
                    "rank": r,
                    "j_tilde": j_tilde,
                    "c": c,
                    "E_total_approx": Etot,
                    "state_residual": state_res,
                    "state_residual_abs": state_res_norm,
                    "state_residual_rel": (
                        state_res_norm / b_norm if b_norm > 0 else np.nan
                    ),
                    "prediction": pred,
                    "direct": g["direct"],
                    "S": S,
                    "chi": chi,
                    "Tmat": Tmat,
                    "db_geom_dx": g["db_geom_dx"],
                    "dS_dx": g["dS_dx"],
                    "direct_derivative_dx": g["direct_derivative_dx"],
                    "DU": DU,
                    "C": C,
                    "Qred": Qred,
                    "Rtri": Rtri,
                    "row_start": row,
                    "row_stop": row + model.n_rx,
                }

                if jacobian:
                    # alpha derivative RHS: C^H(T*E) - (DU)^H(conj(T)*r_s)
                    rhs_a_all = (
                        C.conj().T @ (Tmat * Etot[:, None])
                        - DU.conj().T @ (np.conj(Tmat) * state_res[:, None])
                    )
                    # pose derivative RHS: C^H(chi * db_geom_dx)
                    rhs_p_all = C.conj().T @ (chi[:, None] * g["db_geom_dx"])
                    ledger.charge_reduced(products=3, basis=3)
                    c_a = solve_normal_via_qr(Rtri, rhs_a_all, ledger, q)
                    c_x = solve_normal_via_qr(Rtri, rhs_p_all, ledger, 3)
                    j_a = U @ c_a
                    j_x = U @ c_x
                    A_block = S @ j_a
                    B_block = S @ j_x
                    ledger.charge_reduced(basis=2, products=2, product_columns=q + 3)
                    for ell in range(3):
                        B_block[:, ell] += g["dS_dx"][ell] @ j_tilde
                        ledger.charge_reduced(products=1, product_columns=1)
                    B_block += g["direct_derivative_dx"].T
                    rec["A_block"] = A_block
                    rec["B_block"] = B_block
                    rec["coef_alpha"] = c_a
                    rec["coef_pose"] = c_x

                state_records.append(rec)
                total_parts.append(pred)
                row += model.n_rx

    if not total_parts:
        raise RuntimeError("no reduced states produced")
    total = np.concatenate(total_parts)
    if y_arr.size != total.size:
        raise ValueError(
            f"y length {y_arr.size} != reduced output length {total.size}"
        )
    residual = total - y_arr
    loss = 0.5 * float(np.sum(np.abs(realify(residual, sigma)) ** 2))
    result: dict[str, Any] = {
        "total": total,
        "loss": loss,
        "states": state_records,
        "per_freq": per_freq,
        "max_state_res_abs": float(
            max(r["state_residual_abs"] for r in state_records)
        ),
        "max_state_res_rel": float(
            max(r["state_residual_rel"] for r in state_records)
        ),
        "q": q,
        "n": n,
    }
    if jacobian:
        m = total.size
        A = np.empty((m, q), dtype=np.complex128)
        B = np.empty((m, 3), dtype=np.complex128)
        for rec in state_records:
            sl = slice(rec["row_start"], rec["row_stop"])
            A[sl, :] = rec["A_block"]
            B[sl, :] = rec["B_block"]
        result["A"] = A
        result["B"] = B
        result["J_complex"] = np.hstack([A, B])
    return result


def scaled_reduced_jacobian(
    A: np.ndarray, B: np.ndarray, sigma: float
) -> np.ndarray:
    """Realified fixed-chart data Jacobian in scaled coordinates z."""
    from common import POSE_SCALE

    J = np.hstack([A, B / POSE_SCALE[None, :]])
    return realify(J, sigma)


def reduced_loss_gradient(
    out: dict[str, Any], y: np.ndarray, sigma: float, freq_ids: tuple[int, ...]
) -> np.ndarray:
    """Gradient (scaled coordinates) of the fixed-chart proxy loss."""
    from common import POSE_SCALE, select_rows

    A = out["A"]
    B = out["B"]
    J = np.hstack([A, B / POSE_SCALE[None, :]])
    res = out["total"] - select_rows(y, freq_ids)
    r = realify(res, sigma)
    Jr = realify(J, sigma)
    return Jr.T @ r
