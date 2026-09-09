"""Complex current-space chart seeds/enrichment for the A3 ROM study.

Definitions used here are deliberately explicit:

- ``sensing``: leading right singular vectors of the stacked receiver-to-cell
  sensing rows ``S`` (three poses x 12 receivers).
- ``twofold``: sequential projection ``[V_S, (I-P_S) V_D]`` where ``V_D`` are
  leading right singular vectors of the domain Green/VIE operator ``D``.
  This is the implementable Chen form; it is never called an exact
  intersection.
- ``generic``: the physical state right-hand-side columns ``b`` as seed.
- ``generic_task`` / ``sensing_task`` / ``twofold_task``: the same seed plus
  a fixed budget of residual/primal/adjoint/tangent-RHS enrichment rounds.
- ``block_krylov``: block Krylov columns ``[B, K B, K^2 B, ...]`` with
  ``K = diag(chi) D`` (the natural Born/Neumann expansion operator), where
  ``B`` contains the physical right-hand sides for all poses/illuminations.

No exact full state solve is used in any chart construction.
"""

from __future__ import annotations

import time
from typing import Any, Iterable

import numpy as np
from scipy.linalg import qr, svd
from scipy.linalg import solve_triangular

from common import WorkLedger
from wave import (
    chi_and_Tmat,
    domain_operator,
    state_geometry,
)


def _orthonormalize(
    cols: np.ndarray,
    ledger: WorkLedger,
    max_rank: int | None = None,
    tol_rel: float = 1e-12,
) -> np.ndarray:
    if cols.size == 0:
        return np.zeros((cols.shape[0], 0), dtype=np.complex128)
    W = np.asarray(cols, dtype=np.complex128)
    if W.ndim == 1:
        W = W[:, None]
    max_rank = int(max_rank) if max_rank is not None else W.shape[1]
    max_rank = min(max_rank, W.shape[1])
    t0 = time.perf_counter()
    Q, R, perm = qr(W, mode="economic", pivoting=True, check_finite=False)
    ledger.charge_reduced(
        qr=1, wall=time.perf_counter() - t0, dims=(W.shape[0], W.shape[1])
    )
    diag = np.abs(np.diag(R))
    if diag.size == 0:
        return np.zeros((W.shape[0], 0), dtype=np.complex128)
    ref = max(float(diag[0]), 1e-300)
    keep = int(np.count_nonzero(diag > tol_rel * ref))
    keep = min(keep, max_rank)
    return Q[:, :keep].copy()


def _append_orth(
    U: np.ndarray,
    W: np.ndarray,
    ledger: WorkLedger,
    max_rank: int,
) -> np.ndarray:
    """Append columns of W to orthonormal U in order (stable MGS)."""
    U = np.asarray(U, dtype=np.complex128)
    W = np.asarray(W, dtype=np.complex128)
    if W.ndim == 1:
        W = W[:, None]
    out = U.copy()
    added = 0
    for j in range(W.shape[1]):
        if out.shape[1] >= int(max_rank):
            break
        v = W[:, j].copy()
        if out.shape[1]:
            for _pass in range(2):  # re-orthogonalise for stability
                v = v - out @ (out.conj().T @ v)
        nrm = float(np.linalg.norm(v))
        if nrm > 1e-13 * max(1.0, float(np.linalg.norm(W[:, j]))):
            out = np.hstack([out, (v / nrm)[:, None]])
            added += 1
    ledger.charge_reduced(products=int(W.shape[1]), basis=added)
    return out[:, : int(max_rank)].copy()


def _stack_sensing_rows(model: Any, geom: dict[str, np.ndarray], fi: int) -> np.ndarray:
    from physics import green_matrix

    k = float(model.ks[int(fi)])
    h2 = (model.h**2) * (k**2)
    rows = [
        h2
        * green_matrix(model.points, geom["rx_world"][p], k)
        for p in range(model.n_pose)
    ]
    return np.vstack(rows)


def _b_matrix(model: Any, chi: np.ndarray, geom: dict[str, np.ndarray], fi: int) -> np.ndarray:
    cols = []
    for p in range(model.n_pose):
        for u in range(model.n_tx):
            g = state_geometry(model, geom, fi, p, u)
            cols.append(chi * g["E_inc"])
    return np.column_stack(cols)


def sensing_chart(
    model: Any,
    geom: dict[str, np.ndarray],
    fi: int,
    r_target: int,
    ledger: WorkLedger,
) -> np.ndarray:
    """Leading right singular vectors of stacked sensing rows."""
    S_stack = _stack_sensing_rows(model, geom, int(fi))
    t0 = time.perf_counter()
    _, s, Vh = svd(S_stack, full_matrices=False, check_finite=False)
    ledger.charge_reduced(
        svd=1, wall=time.perf_counter() - t0, dims=S_stack.shape
    )
    tol = max(1e-10 * float(s[0]), 1e-14) if s.size else 0.0
    eff = int(np.count_nonzero(s > tol)) if s.size else 0
    k = min(max(0, int(r_target)), eff)
    return Vh.conj().T[:, :k].copy()


def twofold_chart(
    model: Any,
    alpha: np.ndarray,
    geom: dict[str, np.ndarray],
    fi: int,
    r_target: int,
    ledger: WorkLedger,
    vs: np.ndarray | None = None,
    vd: np.ndarray | None = None,
) -> np.ndarray:
    """Sequential Twofold ``[V_S, (I-P_S)V_D]`` (not exact intersection)."""
    r_target = max(1, min(int(r_target), model.n_cells))
    if vs is None or vs.shape[1] == 0:
        vs = sensing_chart(model, geom, fi, r_target, ledger)
    P = np.eye(model.n_cells, dtype=np.complex128) - vs @ vs.conj().T
    if vd is None or vd.shape[1] == 0:
        D = domain_operator(model, fi)
        t0 = time.perf_counter()
        _, _, Vh = svd(D, full_matrices=False, check_finite=False)
        ledger.charge_reduced(
            svd=1, wall=time.perf_counter() - t0, dims=(D.shape[0], D.shape[1])
        )
        vd = Vh.conj().T
    if vd.ndim != 2 or vd.shape[0] != model.n_cells:
        raise ValueError("vd must be an n x n_vd right-singular block")
    need = r_target - vs.shape[1]
    if need > 0:
        projected = P @ vd[:, : min(need, vd.shape[1])]
        return _append_orth(vs, projected, ledger, max_rank=r_target)
    return vs[:, :r_target].copy()


def _reduced_states_on_U(
    model: Any,
    alpha: np.ndarray,
    geom: dict[str, np.ndarray],
    fi: int,
    U: np.ndarray,
    ledger: WorkLedger,
) -> list[dict[str, np.ndarray]]:
    """One-round reduced states on U used only for enrichment directions."""
    chi, Tmat = chi_and_Tmat(model, alpha, fi)
    D = domain_operator(model, fi)
    DU = D @ U
    ledger.charge_reduced(products=1, product_columns=int(U.shape[1]))
    C = U - chi[:, None] * DU
    t0 = time.perf_counter()
    Q, Rtri = qr(C, mode="economic", check_finite=False)
    ledger.charge_reduced(
        qr=1, wall=time.perf_counter() - t0, dims=(model.n_cells, U.shape[1])
    )
    QH = Q.conj().T
    out = []
    for p in range(model.n_pose):
        for u in range(model.n_tx):
            g = state_geometry(model, geom, fi, p, u)
            b = chi * g["E_inc"]
            c = solve_triangular(
                Rtri, QH @ b, lower=False, check_finite=False
            )
            ledger.charge_reduced(tri_rhs=1, products=1, basis=1)
            j = U @ c
            state_res = b - C @ c
            Etot = g["E_inc"] + DU @ c
            pred = g["direct"] + g["S"] @ j
            ledger.charge_reduced(basis=3, products=1, product_columns=1)
            out.append(
                {
                    "g": g,
                    "b": b,
                    "j": j,
                    "state_res": state_res,
                    "Etot": Etot,
                    "pred": pred,
                    "chi": chi,
                    "Tmat": Tmat,
                }
            )
    return out


def _task_directions(
    model: Any,
    recs: list[dict[str, np.ndarray]],
    y_sel_rows: list[np.ndarray],
    include_tangent_rhs: bool,
) -> list[np.ndarray]:
    """Residual, output-adjoint and (optionally) tangent-RHS directions."""
    directions: list[np.ndarray] = []
    for idx, rec in enumerate(recs):
        sr = rec["state_res"]
        norm = float(np.linalg.norm(sr))
        if norm > 1e-12:
            directions.append(sr / norm)
        g = rec["g"]
        w = rec["pred"] - y_sel_rows[idx]
        adj = g["S"].conj().T @ w
        a_norm = float(np.linalg.norm(adj))
        if a_norm > 1e-12:
            directions.append(adj / a_norm)
        if include_tangent_rhs:
            # map tangent RHS T diag(E) -- physical X columns, not a solve.
            Wmap = rec["Tmat"] * rec["Etot"][:, None]
            Uo, _, _ = svd(Wmap, full_matrices=False, check_finite=False)
            directions.append(np.ascontiguousarray(Uo[:, 0]))
            if Wmap.shape[1] > 1:
                directions.append(np.ascontiguousarray(Uo[:, 1]))
            for ell in range(3):
                bp = rec["chi"] * g["db_geom_dx"][:, ell]
                nrm = float(np.linalg.norm(bp))
                if nrm > 1e-12:
                    directions.append(bp / nrm)
    return directions


def task_enriched_chart(
    model: Any,
    alpha: np.ndarray,
    geom: dict[str, np.ndarray],
    y_freq: np.ndarray,
    fi: int,
    r_target: int,
    ledger: WorkLedger,
    seed_kind: str = "generic",
    max_rounds: int = 14,
    include_tangent_rhs: bool = False,
    vd_by_freq: dict[int, np.ndarray] | None = None,
) -> np.ndarray:
    """Seed + residual/primal/adjoint enrichment chart (no full solves)."""
    chi, _ = chi_and_Tmat(model, alpha, fi)
    n = model.n_cells
    if seed_kind == "sensing":
        U = sensing_chart(model, geom, fi, r_target, ledger)
    elif seed_kind == "twofold":
        U = twofold_chart(
            model,
            alpha,
            geom,
            fi,
            min(r_target, 64),
            ledger,
            vd=(vd_by_freq or {}).get(int(fi)),
        )
    elif seed_kind == "generic":
        B = _b_matrix(model, chi, geom, fi)
        U = _orthonormalize(B, ledger, max_rank=r_target)
    else:
        raise ValueError(f"unknown seed_kind {seed_kind!r}")
    if U.shape[1] == 0 or U.shape[1] >= r_target:
        if U.shape[1] > r_target:
            U = U[:, :r_target]
        return U.copy()
    y_rows = list(np.asarray(y_freq).reshape(-1).reshape(-1, model.n_rx))
    for _round in range(int(max_rounds)):
        recs = _reduced_states_on_U(model, alpha, geom, fi, U, ledger)
        directions = _task_directions(
            model,
            recs,
            y_rows,
            include_tangent_rhs=include_tangent_rhs,
        )
        if not directions:
            break
        U_new = _append_orth(
            U, np.column_stack(directions), ledger, max_rank=r_target
        )
        grew = U_new.shape[1] > U.shape[1]
        U = U_new
        if not grew or U.shape[1] >= r_target:
            break
    return U.copy()


def block_krylov_chart(
    model: Any,
    alpha: np.ndarray,
    geom: dict[str, np.ndarray],
    fi: int,
    r_target: int,
    ledger: WorkLedger,
) -> np.ndarray:
    """Block Krylov columns of K = diag(chi)D applied to RHS block B."""
    chi, _ = chi_and_Tmat(model, alpha, fi)
    B = _b_matrix(model, chi, geom, fi)
    U = _orthonormalize(B, ledger, max_rank=r_target)
    D = domain_operator(model, fi)
    rounds = 0
    while U.shape[1] < r_target and rounds < 80:
        t0 = time.perf_counter()
        DU = D @ U
        ledger.charge_reduced(
            products=1,
            product_columns=int(U.shape[1]),
            wall=time.perf_counter() - t0,
        )
        W = chi[:, None] * DU
        Unew = _append_orth(U, W, ledger, max_rank=r_target)
        if Unew.shape[1] <= U.shape[1]:
            break
        U = Unew
        rounds += 1
    return U.copy()


def build_chart(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    freq_ids: Iterable[int],
    method: str,
    rank: int,
    y_full: np.ndarray,
    ledger: WorkLedger,
    rank_by_freq: dict[int, int] | None = None,
    include_tangent_rhs: bool = False,
    vd_by_freq: dict[int, np.ndarray] | None = None,
) -> dict[int, np.ndarray]:
    """Build one chart per requested model frequency index."""
    alpha = np.asarray(alpha, dtype=float)
    x = np.asarray(x, dtype=float)
    geom = model._effective_geometry(x)
    charts: dict[int, np.ndarray] = {}
    ids = sorted(set(int(f) for f in freq_ids))
    y_all = np.asarray(y_full, dtype=np.complex128).reshape(-1)
    for fi in ids:
        r = int((rank_by_freq or {}).get(fi, rank))
        if r < 1:
            continue
        y_freq = y_all[int(fi) * 72 : (int(fi) + 1) * 72]
        if method in ("sensing",):
            U = sensing_chart(model, geom, fi, r, ledger)
        elif method in ("twofold",):
            U = twofold_chart(
                model,
                alpha,
                geom,
                fi,
                r,
                ledger,
                vd=(vd_by_freq or {}).get(int(fi)),
            )
        elif method in ("generic", "generic_task"):
            seed = "generic"
            U = task_enriched_chart(
                model,
                alpha,
                geom,
                y_freq,
                fi,
                r,
                ledger,
                seed_kind=seed,
                include_tangent_rhs=include_tangent_rhs,
                vd_by_freq=vd_by_freq,
            )
        elif method in ("sensing_task",):
            U = task_enriched_chart(
                model,
                alpha,
                geom,
                y_freq,
                fi,
                r,
                ledger,
                seed_kind="sensing",
                include_tangent_rhs=include_tangent_rhs,
                vd_by_freq=vd_by_freq,
            )
        elif method in ("twofold_task",):
            U = task_enriched_chart(
                model,
                alpha,
                geom,
                y_freq,
                fi,
                r,
                ledger,
                seed_kind="twofold",
                include_tangent_rhs=include_tangent_rhs,
                vd_by_freq=vd_by_freq,
            )
        elif method in ("block_krylov",):
            U = block_krylov_chart(model, alpha, geom, fi, r, ledger)
        else:
            raise ValueError(f"unknown chart method {method!r}")
        if U.shape[1] > 0:
            charts[int(fi)] = U
    return charts
