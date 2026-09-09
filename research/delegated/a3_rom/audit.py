"""Offline exact-state/derivative audits for the A3 fixed-chart study.

Exact state solves and exact analytic derivatives from the A2 core are used
only as OFFLINE audits.  They are never part of online chart construction.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from common import WorkLedger


def exact_audit_forward(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    freq_ids: tuple[int, ...],
    ledger: WorkLedger | None = None,
) -> dict[str, Any]:
    ledger = ledger or WorkLedger()
    fw = model.forward(
        np.asarray(alpha, dtype=float),
        np.asarray(x, dtype=float),
        list(freq_ids),
        jacobian=True,
    )
    ledger.charge_model_work(fw["work"])
    return {
        "total": fw["total"],
        "states": fw["states"],
        "A": fw["A"],
        "B": fw["B"],
        "work": fw["work"],
    }


def compare_to_exact(
    model: Any,
    alpha: np.ndarray,
    x: np.ndarray,
    reduced: dict[str, Any],
    freq_ids: tuple[int, ...],
) -> dict[str, Any]:
    """Relative errors of a fixed-chart reduced evaluation against direct
    exact states/tangents at the same physical parameters."""
    fw = exact_audit_forward(model, alpha, x, tuple(freq_ids))
    exact_states = {
        (int(s["freq_idx"]), int(s["pose_idx"]), int(s["illum_idx"])): s
        for s in fw["states"]
    }
    state_rel: list[float] = []
    for rec in reduced["states"]:
        es = exact_states[
            (int(rec["freq_idx"]), int(rec["pose_idx"]), int(rec["illum_idx"]))
        ]
        num = float(np.linalg.norm(rec["j_tilde"] - es["j"]))
        den = float(np.linalg.norm(es["j"]))
        state_rel.append(num / den if den > 0 else float("nan"))
    out_rel = float(
        np.linalg.norm(reduced["total"] - fw["total"])
        / max(np.linalg.norm(fw["total"]), 1e-300)
    )
    map_rel = float(
        np.linalg.norm(reduced["A"] - fw["A"], "fro")
        / max(np.linalg.norm(fw["A"], "fro"), 1e-300)
    )
    pose_rel = float(
        np.linalg.norm(reduced["B"] - fw["B"], "fro")
        / max(np.linalg.norm(fw["B"], "fro"), 1e-300)
    )
    return {
        "state_rel_max": float(max(state_rel)) if state_rel else float("nan"),
        "state_rel_list": state_rel,
        "output_rel": out_rel,
        "map_jac_rel": map_rel,
        "pose_jac_rel": pose_rel,
    }

