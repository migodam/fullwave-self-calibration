"""Audited cost accounting and chart builders for the A3 chart-reuse pass.

The A3 ``a3_rom`` modules are imported read-only.  Their ``basis._task_directions``
performs an uncharged map-RHS SVD; this module mirrors the generic/twofold-task
enrichment loop with an explicit, dimensioned charge for every SVD/QR.  It is
not an alternative chart algebra: it calls the same A3 seed/enrichment
primitives so the resulting columns are identical to ``a3_rom.basis``.

Everything here is development evidence only.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.linalg import svd

_HERE = Path(__file__).resolve().parent
_ROM = _HERE.parent / "a3_rom"
if str(_ROM) not in sys.path:
    sys.path.insert(0, str(_ROM))

from common import WorkLedger  # noqa: E402
from basis import (  # noqa: E402
    _append_orth,
    _b_matrix,
    _orthonormalize,
    _reduced_states_on_U,
    sensing_chart,
    twofold_chart,
)


class AuditedLedger(WorkLedger):
    """WorkLedger that additionally records every SVD/QR/LU event.

    Existing A3 helper functions call ``charge_reduced``, so this subclass
    sees their QR/SVD calls and dimensions.  Full model calls report LU counts
    through ``charge_model_work``.  The extra keyword ``phase`` is used by this
    pass's own chart builder and is ignored by existing helpers.
    """

    def __init__(self, domain_n: int | None = None) -> None:
        super().__init__()
        self.domain_n = domain_n
        self.factorization_events: list[dict[str, Any]] = []
        self.rhs_solve_seconds = 0.0
        self.factorization_seconds = 0.0
        self.model_wall_seconds = 0.0
        self.setup_wall_seconds = 0.0

    def _record(
        self,
        kind: str,
        dims: list[int] | tuple[int, ...] | None,
        wall: float,
        phase: str,
    ) -> None:
        if dims is not None:
            dims = [int(v) for v in dims]
        elif kind == "lu" and self.domain_n is not None:
            dims = [int(self.domain_n), int(self.domain_n)]
        self.factorization_events.append(
            {
                "kind": kind,
                "dims": dims,
                "wall_seconds": float(wall),
                "phase": str(phase),
            }
        )

    def charge_reduced(
        self,
        tri_rhs: int = 0,
        qr: int = 0,
        svd: int = 0,
        products: int = 0,
        basis: int = 0,
        product_columns: int = 0,
        wall: float = 0.0,
        dims: list[int] | tuple[int, ...] | None = None,
        phase: str = "reduced",
    ) -> None:
        super().charge_reduced(
            tri_rhs=tri_rhs,
            qr=qr,
            svd=svd,
            products=products,
            basis=basis,
            product_columns=product_columns,
            wall=wall,
            dims=dims,
        )
        if qr:
            self._record("qr", dims, wall, phase)
        if svd:
            self._record("svd", dims, wall, phase)

    def charge_model_work(self, work: dict[str, Any]) -> None:
        n_lu = int(work.get("factorizations", 0))
        for _ in range(n_lu):
            self._record("lu", None, 0.0, "full_model")
        super().charge_model_work(work)
        self.rhs_solve_seconds += float(work.get("rhs_solve_seconds", 0.0))
        self.factorization_seconds += float(work.get("factorization_seconds", 0.0))
        self.model_wall_seconds += float(work.get("wall_seconds", 0.0))

    def charge_setup(self, wall: float) -> None:
        self.setup_wall_seconds += float(wall)
        self.wall_seconds += float(wall)

    def snapshot(self) -> dict[str, Any]:
        out = super().snapshot()
        out["rhs_solve_seconds"] = self.rhs_solve_seconds
        out["factorization_seconds"] = self.factorization_seconds
        out["model_wall_seconds"] = self.model_wall_seconds
        out["setup_wall_seconds"] = self.setup_wall_seconds
        out["factorization_events"] = self.factorization_events
        return out


def _task_directions_charged(
    model: Any,
    recs: list[dict[str, np.ndarray]],
    y_sel_rows: list[np.ndarray],
    include_tangent_rhs: bool,
    ledger: AuditedLedger,
    freq_idx: int,
) -> list[np.ndarray]:
    """Audited mirror of ``a3_rom.basis._task_directions``.

    The only numerical change is that the map-RHS SVD (previously uncharged)
    is timed and charged explicitly with its matrix dimensions.
    """
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
            Wmap = rec["Tmat"] * rec["Etot"][:, None]
            t0 = time.perf_counter()
            Uo, _, _ = svd(Wmap, full_matrices=False, check_finite=False)
            ledger.charge_reduced(
                svd=1,
                wall=time.perf_counter() - t0,
                dims=(Wmap.shape[0], Wmap.shape[1]),
                phase=f"enrich_map_rhs_freq{int(freq_idx)}",
            )
            directions.append(np.ascontiguousarray(Uo[:, 0]))
            if Wmap.shape[1] > 1:
                directions.append(np.ascontiguousarray(Uo[:, 1]))
            for ell in range(3):
                bp = rec["chi"] * g["db_geom_dx"][:, ell]
                nrm = float(np.linalg.norm(bp))
                if nrm > 1e-12:
                    directions.append(bp / nrm)
    return directions


def task_enriched_chart_charged(
    model: Any,
    alpha: np.ndarray,
    geom: dict[str, np.ndarray],
    y_freq: np.ndarray,
    fi: int,
    r_target: int,
    ledger: AuditedLedger,
    seed_kind: str = "generic",
    max_rounds: int = 14,
    include_tangent_rhs: bool = False,
    vd: np.ndarray | None = None,
) -> np.ndarray:
    """Audited generic/twofold-task enrichment chart (identical columns)."""
    from basis import chi_and_Tmat

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
            vd=vd,
        )
    elif seed_kind == "generic":
        B = _b_matrix(model, chi, geom, fi)
        U = _orthonormalize(B, ledger, max_rank=r_target)
    else:
        raise ValueError(f"unknown seed_kind {seed_kind!r}")
    if U.shape[1] == 0 or U.shape[1] >= r_target:
        if U.shape[1] > r_target:
            U = U[:, : r_target]
        return U.copy()
    y_rows = list(np.asarray(y_freq).reshape(-1).reshape(-1, model.n_rx))
    for _round in range(int(max_rounds)):
        recs = _reduced_states_on_U(model, alpha, geom, fi, U, ledger)
        directions = _task_directions_charged(
            model,
            recs,
            y_rows,
            include_tangent_rhs=include_tangent_rhs,
            ledger=ledger,
            freq_idx=fi,
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


def build_chart_charged(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    freq_ids: Iterable[int],
    method: str,
    rank: int,
    y_full: np.ndarray,
    ledger: AuditedLedger,
    rank_by_freq: dict[int, int] | None = None,
    include_tangent_rhs: bool = True,
    vd_by_freq: dict[int, np.ndarray] | None = None,
) -> dict[int, np.ndarray]:
    """Chart build with identical rules to ``a3_rom.basis.build_chart``.

    ``method`` is restricted to ``generic_task``/``twofold_task`` (the two ROM
    variants compared in this pass); unsupported names fail loudly.
    """
    if method not in ("generic_task", "twofold_task"):
        raise ValueError(f"chart reuse pass only builds {method!r} charts")
    alpha = np.asarray(alpha, dtype=float)
    x = np.asarray(x, dtype=float)
    geom = model._effective_geometry(x)
    charts: dict[int, np.ndarray] = {}
    ids = sorted(set(int(f) for f in freq_ids))
    y_all = np.asarray(y_full, dtype=np.complex128).reshape(-1)
    seed_kind = "generic" if method == "generic_task" else "twofold"
    for fi in ids:
        r = int((rank_by_freq or {}).get(fi, rank))
        if r < 1:
            charts[int(fi)] = np.zeros((model.n_cells, 0), dtype=np.complex128)
            continue
        y_freq = y_all[int(fi) * 72 : (int(fi) + 1) * 72]
        U = task_enriched_chart_charged(
            model,
            alpha,
            geom,
            y_freq,
            fi,
            r,
            ledger,
            seed_kind=seed_kind,
            include_tangent_rhs=include_tangent_rhs,
            vd=(vd_by_freq or {}).get(int(fi)),
        )
        if U.shape[1] > 0:
            charts[int(fi)] = U
        else:
            charts[int(fi)] = np.zeros((model.n_cells, 0), dtype=np.complex128)
    return charts
