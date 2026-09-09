#!/usr/bin/env python
"""Fixed-chart rank-error study on the A3 development grid16 model.

This is an OFFLINE audit: direct exact states and exact derivatives from the
A2 core are used only to measure chart error after the fact.  They are never
used inside chart construction.

Run:
    experiments/.../.venv/bin/python research/delegated/a3_rom/rank_study.py
"""

from __future__ import annotations

import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.linalg import svd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common import (  # noqa: E402
    FREQ_IDS,
    WorkLedger,
    complex_from_json,
    config_json,
    guard_seed,
    jsonable,
    make_alpha_true,
    make_model,
    make_scene,
    pose_init,
    write_json,
)
from basis import build_chart  # noqa: E402
from wave import domain_operator, fixed_chart_evaluate  # noqa: E402

RANK_LADDER = (16, 24, 32, 48, 64, 96, 128, 160, 192, 256)
STATE_TOL = 1e-3
OUTPUT_TOL = 1e-3
MAP_JAC_TOL = 1e-2
POSE_JAC_TOL = 1e-2
METHODS = (
    "sensing",
    "twofold",
    "generic_task",
    "block_krylov",
    "sensing_task",
    "twofold_task",
)


def _err_vs_exact(reduced: dict, exact: dict) -> dict:
    es = {
        (int(s["freq_idx"]), int(s["pose_idx"]), int(s["illum_idx"])): s
        for s in exact["states"]
    }
    state_vals = []
    for rec in reduced["states"]:
        key = (int(rec["freq_idx"]), int(rec["pose_idx"]), int(rec["illum_idx"]))
        s = es[key]
        den = float(np.linalg.norm(s["j"]))
        state_vals.append(float(np.linalg.norm(rec["j_tilde"] - s["j"])) / den)
    return {
        "state_rel": float(max(state_vals)),
        "output_rel": float(
            np.linalg.norm(reduced["total"] - exact["total"])
            / max(np.linalg.norm(exact["total"]), 1e-300)
        ),
        "map_jac_rel": float(
            np.linalg.norm(reduced["A"] - exact["A"], "fro")
            / max(np.linalg.norm(exact["A"], "fro"), 1e-300)
        ),
        "pose_jac_rel": float(
            np.linalg.norm(reduced["B"] - exact["B"], "fro")
            / max(np.linalg.norm(exact["B"], "fro"), 1e-300)
        ),
    }


def passed(err: dict) -> bool:
    return (
        err["state_rel"] <= STATE_TOL
        and err["output_rel"] <= OUTPUT_TOL
        and err["map_jac_rel"] <= MAP_JAC_TOL
        and err["pose_jac_rel"] <= POSE_JAC_TOL
    )


def main() -> None:
    model = make_model(n=16)
    scenes = {s: make_scene(s) for s in (4101, 4102, 4103)}
    y = {s: complex_from_json(scenes[s]["y"]) for s in scenes}

    # D is frequency-fixed, so its right singular block is cached once per
    # frequency and its reuse count is disclosed (not charged to every chart).
    vd_ledger = WorkLedger()
    vd_by_freq: dict[int, np.ndarray] = {}
    for fi in FREQ_IDS:
        D = domain_operator(model, fi)
        _, _, Vh = svd(D, full_matrices=False, check_finite=False)
        vd_ledger.charge_reduced(svd=1, dims=(D.shape[0], D.shape[1]))
        vd_by_freq[fi] = Vh.conj().T

    exact_by_seed: dict[int, dict[int, dict]] = {}
    exact_ledger = WorkLedger()
    for seed in scenes:
        exact_by_seed[seed] = {}
        alpha = make_alpha_true(seed)
        x = pose_init(seed)
        for fi in FREQ_IDS:
            fw = model.forward(alpha, x, [fi], jacobian=True)
            exact_ledger.charge_model_work(fw["work"])
            exact_by_seed[seed][fi] = {
                "alpha": alpha,
                "x": x,
                "total": fw["total"],
                "states": fw["states"],
                "A": fw["A"],
                "B": fw["B"],
            }

    rows: list[dict] = []
    first_pass: dict[str, dict[str, dict]] = {}
    total_t0 = time.perf_counter()
    chart_ledger = WorkLedger()
    for seed, exfreq in exact_by_seed.items():
        alpha = exfreq[FREQ_IDS[0]]["alpha"]
        x = exfreq[FREQ_IDS[0]]["x"]
        y_full = y[seed]
        for fi in FREQ_IDS:
            exact = exfreq[fi]
            for method in METHODS:
                if method not in first_pass:
                    first_pass[method] = {}
                for r in RANK_LADDER:
                    charts = build_chart(
                        model,
                        alpha,
                        x,
                        (fi,),
                        method,
                        r,
                        y_full,
                        chart_ledger,
                        include_tangent_rhs=True,
                        vd_by_freq=vd_by_freq,
                    )
                    if fi not in charts or charts[fi].shape[1] == 0:
                        continue
                    actual_rank = int(charts[fi].shape[1])
                    red = fixed_chart_evaluate(
                        model,
                        alpha,
                        x,
                        charts,
                        y_full,
                        1.0,
                        (fi,),
                        jacobian=True,
                    )
                    err = _err_vs_exact(red, exact)
                    row = {
                        "seed": seed,
                        "alpha_sample": list(map(float, alpha)),
                        "x_sample": list(map(float, x)),
                        "frequency_k": float(model.ks[fi]),
                        "freq_idx": int(fi),
                        "method": method,
                        "target_rank": int(r),
                        "actual_rank": actual_rank,
                        **{f"err_{k}": v for k, v in err.items()},
                        "passed_all": bool(passed(err)),
                        "chart_build_wall": float(chart_ledger.wall_seconds),
                    }
                    rows.append(row)
                    fp = first_pass[method].setdefault(
                        f"k={model.ks[fi]:g},seed={seed}", None
                    )
                    if fp is None and passed(err):
                        first_pass[method][f"k={model.ks[fi]:g},seed={seed}"] = {
                            "target_rank": int(r),
                            "actual_rank": actual_rank,
                            "err": err,
                        }

    # Remove the placeholder None entries for samples that never passed.
    for method in list(first_pass):
        first_pass[method] = {
            k: v for k, v in first_pass[method].items() if v is not None
        }

    # Summaries: aggregate first-pass ranks per method/frequency.
    summary: dict[str, dict[str, dict]] = {}
    for method in METHODS:
        summary[method] = {}
        for key, val in first_pass.get(method, {}).items():
            freq = key.split(",")[0]
            summary[method].setdefault(freq, {})[key] = val

    payload = {
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "COMPLETE",
        "config": config_json(n=16),
        "thresholds": {
            "state_rel": STATE_TOL,
            "output_rel": OUTPUT_TOL,
            "map_jac_rel": MAP_JAC_TOL,
            "pose_jac_rel": POSE_JAC_TOL,
        },
        "rank_ladder": list(RANK_LADDER),
        "methods": list(METHODS),
        "samples": sorted(scenes),
        "rows": jsonable(rows),
        "first_pass": jsonable(first_pass),
        "work": {
            "chart_total_wall": chart_ledger.wall_seconds,
            "chart_qr": chart_ledger.qr_factorizations,
            "chart_svd": chart_ledger.svd_factorizations,
            "chart_products": chart_ledger.full_operator_products,
            "exact_audit_full_rhs": exact_ledger.full_rhs_columns,
            "exact_audit_lu": exact_ledger.lu_factorizations,
            "vd_svd_setup": vd_ledger.svd_factorizations,
        },
        "wall_seconds": time.perf_counter() - total_t0,
    }
    write_json(HERE / "rank_study_raw.json", payload)

    # Markdown summary
    lines = [
        "# A3 fixed-chart rank-error study (development, grid16)",
        "",
        f"Timestamp: {payload['timestamp']}",
        f"Wall: {payload['wall_seconds']:.2f} s. "
        f"Three material samples (seeds 4101-4103), k=3,6,12 rad/m.",
        "",
        "Thresholds: state/output rel error <= 1e-3; map/pose Jacobian rel "
        "Frobenius error <= 1e-2. First pass is the lowest ladder rank that "
        "passes all four on that sample/frequency. Actual rank is shown when "
        "it differs from the target (sensing is capped by its 36-row SVD).",
        "",
        "| Method | k=3 first-pass ranks | k=6 | k=12 |",
        "|---|---|---|---|",
    ]
    for method in METHODS:
        cells = []
        for freq in ("k=3", "k=6", "k=12"):
            vals = summary.get(method, {}).get(freq, {})
            ranks = [v["target_rank"] for v in vals.values()] if vals else []
            if ranks:
                cells.append(" / ".join(str(v) for v in ranks))
            else:
                cells.append("未通过至256")
        lines.append(f"| {method} | {cells[0]} | {cells[1]} | {cells[2]} |")
    lines += [
        "",
        "## Notes and limitations",
        "",
        "- Chart construction never uses exact states/derivatives; exact "
        "forward/Jacobian solves are offline audits only.",
        "- The Twofold seed is the implementable sequential projection "
        "``[Vs,(I-Ps)Vd]``, not an exact intersection.",
        "- ``sensing`` is capped by the numerical rank of the 36-row sensing "
        "matrix; the sensing/twofold rows labelled *_task include the same "
        "task-residual/tangent-RHS enrichment budget as ``generic_task``.",
        "- D right singular blocks are frequency-fixed and cached once; the "
        "single setup SVD per frequency is reported separately, not charged "
        "to each chart.",
        "- Raw rows are in ``rank_study_raw.json``. This is a development "
        "rank-error audit, not a final sweep and not SOM falsification.",
    ]
    (HERE / "rank_study_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("A3 rank-error study complete")
    for method in METHODS:
        r3 = [v["target_rank"] for v in summary.get(method, {}).get("k=3", {}).values()]
        r6 = [v["target_rank"] for v in summary.get(method, {}).get("k=6", {}).values()]
        r12 = [v["target_rank"] for v in summary.get(method, {}).get("k=12", {}).values()]
        print(
            f"- {method}: k3={r3 or ['fail']}, k6={r6 or ['fail']}, "
            f"k12={r12 or ['fail']}"
        )


if __name__ == "__main__":
    main()
