"""Certified SOM-informed reduced physical state solver for A2 E4."""
from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.linalg import qr, solve_triangular, svd

from common import (
    EPS0,
    CostLedger,
)


def tau_at(model: Any, fi: int) -> float:
    k = float(model.ks[fi])
    omega = k * 299792458.0
    return float(model.config.loss_coefficient) / (omega * EPS0)


def chi_and_t(model: Any, alpha: np.ndarray, fi: int) -> tuple[np.ndarray, np.ndarray]:
    tau = tau_at(model, fi)
    loss = 1.0 + 1j * tau
    u = model.basis @ alpha
    return u * loss, model.basis * loss


def domain_operator(model: Any, fi: int) -> np.ndarray:
    return model._domain_operator(int(fi))


def _state_geometry(model: Any, geom: dict[str, np.ndarray], fi: int, p: int, u: int):
    """Green/sensing geometry fields for one state (no current solves)."""
    from physics import (
        green_grad_first,
        green_grad_source,
        green_matrix,
    )

    k = float(model.ks[fi])
    rx_world = geom["rx_world"][p]
    tx_world = geom["tx_world"][p]
    drx_x = geom["drx_x"][p]
    dtx_x = geom["dtx_x"][p]
    S = (k**2) * (model.h**2) * green_matrix(model.points, rx_world, k)
    grad_rx = green_grad_first(model.points, rx_world, k)
    dS_dx = [
        (k**2)
        * (model.h**2)
        * np.einsum("asd,ad->as", grad_rx, drx_x[l])
        for l in range(3)
    ]
    tx_u = tx_world[u]
    E_inc = green_matrix(tx_u[None, :], model.points, k)[:, 0]
    direct = green_matrix(tx_u[None, :], rx_world, k)[:, 0]
    grad_src = green_grad_source(model.points, tx_u, k)
    db_dx = np.empty((model.n_cells, 3), dtype=np.complex128)
    for l in range(3):
        db_dx[:, l] = np.einsum(
            "nd,d->n", grad_src, dtx_x[l, u, :]
        )
    grad_at_rx = green_grad_first(tx_u[None, :], rx_world, k)[:, 0, :]
    rel_dx = drx_x - dtx_x[:, u, None, :]
    direct_dx = np.einsum("ad,lad->la", grad_at_rx, rel_dx)
    return {
        "S": S,
        "E_inc": E_inc,
        "direct": direct,
        "dS_dx": dS_dx,
        "db_dx": db_dx,
        "direct_derivative_dx": direct_dx,
        "rx_world": rx_world,
        "tx_world": tx_world,
    }


def build_sensing_basis(
    model: Any,
    x: np.ndarray,
    freq_ids: tuple[int, ...],
    rank: int,
    ledger: CostLedger,
    tsom_append: bool = False,
    tsom_target_fraction: float = 0.0,
) -> dict[str, Any]:
    """Leading right singular vectors of the stacked per-pose sensing rows.

    Optional sequential TSOM block P_(S-)V_(D,+) is appended with RRQR and
    the actual retained rank recorded; it is never called an exact
    intersection.
    """
    from physics import green_matrix

    n = model.n_cells
    rank = max(1, min(int(rank), n))
    geom = model._effective_geometry(x)
    rows = []
    for fi in sorted(freq_ids):
        k = float(model.ks[fi])
        for p in range(3):
            rx_world = geom["rx_world"][p]
            S = (k**2) * (model.h**2) * green_matrix(model.points, rx_world, k)
            rows.append(S)
    S_stack = np.vstack(rows)
    t0 = time.perf_counter()
    _, s, Vh = svd(S_stack, full_matrices=False, check_finite=False)
    ledger.charge_svd(1)
    ledger.add_wall(time.perf_counter() - t0)
    V = Vh.conj().T
    if s.size == 0:
        raise RuntimeError("empty sensing SVD")
    s0 = float(np.max(s))
    tol = max(1e-10 * s0, 1e-14)
    ldet = int(np.count_nonzero(s > tol))
    U = V[:, :rank].copy()
    actual_rank = U.shape[1]
    tsom_rank = 0
    if tsom_append and rank < n:
        target = max(1, int(round(tsom_target_fraction * rank)))
        extra: list[np.ndarray] = []
        P = np.eye(n, dtype=np.complex128) - U @ U.conj().T
        for fi in sorted(freq_ids):
            D = domain_operator(model, fi)
            t0 = time.perf_counter()
            _, _, VhD = svd(D, full_matrices=False, check_finite=False)
            ledger.charge_svd(1)
            ledger.add_wall(time.perf_counter() - t0)
            VD = VhD.conj().T[:, : min(target, n)]
            extra.append(P @ VD)
        if extra:
            W = np.hstack(extra)
            t0 = time.perf_counter()
            Qp, Rp, perm = qr(W, mode="economic", pivoting=True, check_finite=False)
            ledger.charge_rrqr(1)
            ledger.add_wall(time.perf_counter() - t0)
            diag = np.abs(np.diag(Rp))
            d0 = float(diag[0]) if diag.size else 0.0
            keep = int(np.count_nonzero(diag > max(1e-10 * d0, 1e-13)))
            keep = min(keep, n - rank)
            if keep > 0:
                U = np.hstack([U, Qp[:, :keep]])
                tsom_rank = keep
                actual_rank += keep
    # enforce complex-orthonormal columns
    t0 = time.perf_counter()
    U, _ = qr(U, mode="economic", check_finite=False)
    ledger.charge_rrqr(1)
    ledger.add_wall(time.perf_counter() - t0)
    if U.shape[1] > rank + tsom_rank:
        U = U[:, : rank + tsom_rank]
    return {
        "U": U,
        "actual_rank": int(U.shape[1]),
        "ldet": ldet,
        "singular_values": s[: min(len(s), 96)].tolist(),
        "tsom_rank": tsom_rank,
        "sensing_target": rank,
    }


def reduced_evaluate(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    U: np.ndarray,
    freq_ids: tuple[int, ...],
    ledger: CostLedger,
    need_jacobian: bool,
) -> dict[str, Any]:
    """Approximate physical states/sensitivities on the frozen basis U.

    All state coefficients minimise ||M U c - b|| in the full state space
    (physical state residual, not measurement fit).  The derivative solve
    is the same reduced tangent equation and is reported as an inexact
    physical Jacobian covered by the parent residual certificate.
    """
    n = model.n_cells
    r = int(U.shape[1])
    geom = model._effective_geometry(x)
    total_parts = []
    block_rows = 0
    state_records = []
    q = model.n_alpha
    # each state row contributes 12 receivers
    rows_per_state = 12

    for fi in sorted(freq_ids):
        k = float(model.ks[fi])
        chi, Tmat = chi_and_t(model, alpha, fi)
        D = domain_operator(model, fi)
        DU = D @ U
        ledger.charge_basis_products(r)
        R = U - chi[:, None] * DU
        t0 = time.perf_counter()
        Qred, Rred = qr(R, mode="economic", check_finite=False)
        ledger.charge_reduced(factorizations=1, rank=r)
        ledger.add_wall(time.perf_counter() - t0)
        QH = Qred.conj().T
        for p in range(3):
            for u in range(2):
                g = _state_geometry(model, geom, fi, p, u)
                S = g["S"]
                b = chi * g["E_inc"]
                rhs_c = QH @ b
                t0 = time.perf_counter()
                c = solve_triangular(Rred, rhs_c, lower=False, check_finite=False)
                ledger.charge_reduced(solves=1, rank=r)
                ledger.add_wall(time.perf_counter() - t0)
                jt = U @ c
                ledger.charge_basis_products(1)
                Etot = g["E_inc"] + DU @ c
                ledger.charge_basis_products(1)
                z_state = R @ c - b
                ledger.charge_basis_products(1)
                pred = g["direct"] + S @ jt
                ledger.charge_operator_products(1)
                rec = {
                    "freq_idx": fi,
                    "pose_idx": p,
                    "illum_idx": u,
                    "S": S,
                    "chi": chi,
                    "Tmat": Tmat,
                    "E_inc": g["E_inc"],
                    "E_total_approx": Etot,
                    "j_tilde": jt,
                    "c": c,
                    "b": b,
                    "state_residual": z_state,
                    "direct": g["direct"],
                    "prediction": pred,
                    "dS_dx": g["dS_dx"],
                    "db_dx": g["db_dx"],
                    "direct_derivative_dx": g["direct_derivative_dx"],
                    "R": R,
                    "Qred": Qred,
                    "Rred": Rred,
                    "row_start": block_rows,
                    "row_stop": block_rows + rows_per_state,
                }
                if need_jacobian:
                    # alpha tangent RHS diag(Etot) T
                    rhs_a = Etot[:, None] * Tmat
                    rhs_b = g["db_dx"] * chi[:, None]
                    all_rhs = np.hstack([rhs_a, rhs_b])  # n x (q+3)
                    t0 = time.perf_counter()
                    coef = solve_triangular(
                        Rred, QH @ all_rhs, lower=False, check_finite=False
                    )
                    ledger.charge_reduced(solves=int(all_rhs.shape[1]), rank=r)
                    ledger.add_wall(time.perf_counter() - t0)
                    coeff_a = coef[:, :q]
                    coeff_b = coef[:, q:]
                    t_a = U @ coeff_a
                    t_b = U @ coeff_b
                    ledger.charge_basis_products(q+3)
                    Ared = S @ t_a
                    Bred = S @ t_b
                    ledger.charge_operator_products(q+3)
                    for ell in range(3):
                        Bred[:, ell] += g["dS_dx"][ell] @ rec["j_tilde"]
                        ledger.charge_operator_products(1)
                    Bred += g["direct_derivative_dx"].T
                    der_res = R @ coef - all_rhs
                    ledger.charge_basis_products(q+3)
                    rec["A_block"] = Ared
                    rec["B_block"] = Bred
                    rec["derivative_residual"] = der_res
                    rec["alpha_coef"] = coeff_a
                    rec["pose_coef"] = coeff_b
                state_records.append(rec)
                block_rows += rows_per_state

    total = np.concatenate([r["prediction"] for r in state_records])
    out: dict[str, Any] = {"total": total, "states": state_records, "rank": r}
    if need_jacobian:
        m = total.size
        A = np.empty((m, q), dtype=np.complex128)
        B = np.empty((m, 3), dtype=np.complex128)
        for rec in state_records:
            sl = slice(rec["row_start"], rec["row_stop"])
            A[sl, :] = rec["A_block"]
            B[sl, :] = rec["B_block"]
        out["A"] = A
        out["B"] = B
    return out


def passive_certificate(model: Any, alpha: np.ndarray, fi: int) -> Any:
    """Parent passive-medium certificate from certificates.py."""
    from certificates import passive_bound

    u = model.basis @ alpha
    return passive_bound(u, tau_at(model, fi), float(model.ks[fi]), float(model.h))
