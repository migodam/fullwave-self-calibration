#!/usr/bin/env python
"""Development nonlinear joint material/pose runs for the A3 ROM study.

Only DEVELOPMENT seeds 4101-4104 are used.  Exact physical acceptance is
shared by direct and ROM methods; ROM reduced-step acceptance is recorded per
frequency.  No large final sweep and no SOM-falsification claim is made.

Run:
    experiments/.../.venv/bin/python research/delegated/a3_rom/run_nonlinear.py
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common import (  # noqa: E402
    DEVELOPMENT_SEEDS,
    complex_from_json,
    config_json,
    make_model,
    make_scene,
    write_json,
)
from solver import run_solver  # noqa: E402

METHODS = (
    "direct_gn",
    "direct_adjoint",
    "sensing",
    "twofold",
    "generic_task",
    "sensing_task",
    "twofold_task",
    "block_krylov",
)
# Continuation/no-continuation comparison uses the same three-frequency data.
CONTINUATIONS = {
    "direct_gn": (True, False),
    "direct_adjoint": (True, False),
    "generic_task": (True, False),
    "sensing_task": (True, False),
    "twofold_task": (True, False),
    "block_krylov": (True, False),
    "sensing": (True,),
    "twofold": (True,),
}
BUDGET_RHS = 8000


def main() -> None:
    model = make_model(n=16)
    scenes = {seed: make_scene(seed) for seed in DEVELOPMENT_SEEDS}
    records: list[dict] = []
    t0 = time.perf_counter()
    runs = 0
    with (HERE / "records_nonlinear.jsonl").open("w", encoding="utf-8") as fh:
        for seed in DEVELOPMENT_SEEDS:
            scene = scenes[seed]
            for method in METHODS:
                for continuation in CONTINUATIONS[method]:
                    wall_before = time.perf_counter()
                    rec = run_solver(
                        model,
                        scene,
                        method,
                        continuation,
                        seed,
                        budget_rhs=BUDGET_RHS,
                    )
                    rec["timestamp"] = datetime.now().astimezone().isoformat(
                        timespec="seconds"
                    )
                    rec["scene_seed"] = seed
                    rec["wall_total_seconds"] = time.perf_counter() - wall_before
                    fh.write(json.dumps(rec, sort_keys=True) + "\n")
                    fh.flush()
                    records.append(rec)
                    runs += 1
                    print(
                        f"[{runs:02d}] seed={seed} method={method} "
                        f"cont={int(continuation)} status={rec['status']} "
                        f"loss={rec['final_loss']:.2f} "
                        f"pose={rec['pose_metric_error_m']:.4f} "
                        f"red={rec['reduced']['accepted_reduced_total']} "
                        f"wall={rec['wall_seconds']:.2f}s",
                        flush=True,
                    )

    summary = summarize(records)
    payload = {
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "DEVELOPMENT_COMPLETE",
        "config": config_json(n=16),
        "budget_rhs": BUDGET_RHS,
        "seeds": list(DEVELOPMENT_SEEDS),
        "methods": list(METHODS),
        "summary": summary,
        "total_wall_seconds": time.perf_counter() - t0,
        "runs": runs,
        "notes": [
            "Development seeds 4101-4104 only; no E4 final seeds.",
            "Chart ranks are predeclared in solver.CHART_RANKS from the "
            "offline rank audit; exact state/derivatives never enter chart "
            "construction.",
            "A no_improving/fallback stop at a low exact loss is an exact "
            "controller termination, not a ROM shortcut; accepted-reduced "
            "counts are separate.",
            "These runs do not claim a wall-time speedup or SOM falsification.",
        ],
    }
    write_json(HERE / "nonlinear_summary.json", payload)
    write_markdown(summary, payload)
    print(f"Nonlinear development runs complete: {runs} runs")


def summarize(records: list[dict]) -> dict:
    out: dict = {}
    for rec in records:
        key = f"{rec['method']}::cont={int(rec['continuation'])}"
        item = out.setdefault(
            key,
            {
                "loss": [],
                "pose": [],
                "alpha_rmse": [],
                "accepted_reduced_total": [],
                "accepted_by_freq": {},
                "fallback": [],
                "status": [],
                "wall": [],
                "full_rhs": [],
            },
        )
        item["loss"].append(rec["final_loss"])
        item["pose"].append(rec["pose_metric_error_m"])
        item["alpha_rmse"].append(rec["alpha_rmse"])
        item["accepted_reduced_total"].append(
            rec["reduced"]["accepted_reduced_total"]
        )
        for fi, count in rec["reduced"]["accepted_reduced_by_freq"].items():
            item["accepted_by_freq"].setdefault(str(fi), []).append(int(count))
        item["fallback"].append(rec["reduced"]["fallback_direct_steps"])
        item["status"].append(rec["status"])
        item["wall"].append(rec["wall_seconds"])
        item["full_rhs"].append(rec["ledger"]["full_rhs_columns"])

    def mean(v):
        return float(sum(v) / len(v)) if v else None

    agg = {}
    for key, item in out.items():
        agg[key] = {
            "n": len(item["loss"]),
            "final_loss_mean": mean(item["loss"]),
            "pose_metric_mean": mean(item["pose"]),
            "alpha_rmse_mean": mean(item["alpha_rmse"]),
            "accepted_reduced_total_mean": mean(item["accepted_reduced_total"]),
            "accepted_by_freq_mean": {
                k: mean(v) for k, v in item["accepted_by_freq"].items()
            },
            "fallback_direct_mean": mean(item["fallback"]),
            "full_rhs_mean": mean(item["full_rhs"]),
            "wall_mean": mean(item["wall"]),
            "statuses": sorted(set(item["status"])),
        }
    return agg


def write_markdown(summary: dict, payload: dict) -> None:
    lines = [
        "# A3 development nonlinear ROM runs (grid16)",
        "",
        f"Timestamp: {payload['timestamp']}; total wall "
        f"{payload['total_wall_seconds']:.1f} s; {payload['runs']} runs.",
        "",
        "All runs use identical final multifrequency data (k=3,6,12). "
        "``cont=1`` is the low-to-high stage schedule; ``cont=0`` is the same "
        "data added at once. Reduced acceptance counts are exact-objective "
        "accepted ROM steps, grouped by every frequency active in that stage.",
        "",
        "| method / continuation | n | final loss mean | pose m mean | "
        "alpha RMSE mean | accepted reduced mean | fallback mean | full RHS mean | wall s mean | statuses |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for key in sorted(summary):
        s = summary[key]
        by_freq = s["accepted_by_freq_mean"]
        red = f"{s['accepted_reduced_total_mean']:.1f}"
        if by_freq:
            red += " (" + ", ".join(
                f"k{int(f):g}:{v:.1f}" for f, v in sorted(by_freq.items())
            ) + ")"
        lines.append(
            f"| {key} | {s['n']} | {s['final_loss_mean']:.2f} | "
            f"{s['pose_metric_mean']:.4f} | {s['alpha_rmse_mean']:.4f} | "
            f"{red} | {s['fallback_direct_mean']:.1f} | "
            f"{s['full_rhs_mean']:.0f} | {s['wall_mean']:.2f} | "
            f"{','.join(s['statuses'])} |"
        )
    lines += [
        "",
        "## Interpretation boundaries",
        "",
        "- A ``no_improving_trial_and_fallback`` stop at low exact loss is "
        "the exact acceptance controller refusing further numerical "
        "improvement; it is not a solver failure by itself.",
        "- Accepted ROM steps in the single all-frequency schedule directly "
        "include k=12 (model frequency 3); continuation accepted-step counts "
        "are split by active stage frequency.",
        "- Chart rebuild is charged to every accepted iterate; full/reduced "
        "work are not merged into a speedup claim.",
        "- These are development outcomes only (seeds 4101-4104). No final "
        "sweep, no SOM-falsification claim.",
    ]
    (HERE / "nonlinear_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
