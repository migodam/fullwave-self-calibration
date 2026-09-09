"""Runner for the a2_highdim exploratory 2D imaging comparison.

Usage (single-thread CPU):
    python runner.py --mode smoke
    python runner.py --mode tune            # writes frozen.json when done
    python runner.py --mode test            # requires frozen.json
    python runner.py --mode mismatch        # requires frozen.json
    python runner.py --mode summarize       # records.jsonl + figures -> summary
    python runner.py --mode all             # tune, test, mismatch, summarize
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

# Enforce single-thread CPU before NumPy/BLAS initialisation in this process.
for _var in (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_var, "1")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (  # noqa: E402
    DIR,
    MISMATCH_SEED,
    TEST_SEEDS,
    TUNING_SEEDS,
    complex_from_json,
    default_settings,
    guard_seed,
    jsonable,
    load_frozen,
    make_model,
    make_scene,
    mismatch_init_x0,
    nominal_sigma,
    pose_error_radius,
    settings_hash,
    write_frozen,
)
from evaluate import evaluate_estimate  # noqa: E402
from mismatch import make_mismatch_scene, mismatch_scene_arrays  # noqa: E402
from solver import METHODS, run_method  # noqa: E402

RECORDS_PATH = DIR / "records.jsonl"


def _record_line(rec: dict[str, Any]) -> str:
    return json.dumps(jsonable(rec), sort_keys=True) + "\n"


def _save_record(rec: dict[str, Any]) -> None:
    with RECORDS_PATH.open("a", encoding="utf-8") as fh:
        fh.write(_record_line(rec))


def _run_one(
    phase: str,
    settings: dict[str, Any],
    method: str,
    seed: int,
    radius_index: int | None = None,
    mismatch_kind: str | None = None,
    wall_cap_seconds: float = 900.0,
) -> dict[str, Any]:
    guard_seed(seed)
    model_inv = make_model(20, "gaussian49")
    if mismatch_kind is None:
        scene = make_scene(seed)
        offset = pose_error_radius(seed, int(radius_index))
        x0 = np.asarray(offset["x0"], dtype=float)
        radius_label = offset["radius_label"]
    else:
        scene = make_mismatch_scene(mismatch_kind)
        x0 = mismatch_init_x0()
        radius_label = "0.125lambda"
    sigma = float(scene["sigma"])
    y = complex_from_json(scene["y"])
    x_true = np.asarray(scene["x_true"], dtype=float)

    outcome = run_method(
        model_inv,
        y,
        sigma,
        method,
        x0,
        settings,
        seed,
        radius_label,
        x_true,
        wall_cap_seconds=wall_cap_seconds,
    )
    rec: dict[str, Any] = {
        "phase": phase,
        "seed": int(seed),
        "radius_index": radius_index,
        "radius_label": radius_label,
        "mismatch_kind": mismatch_kind,
        "method": method,
        "is_oracle": bool(outcome.is_oracle),
        "settings_hash": settings_hash(settings),
        "x_init": x0.tolist(),
        "alpha_init": [float(settings["alpha_init"])] * 49,
        "status": outcome.status,
        "all_stages_visited": bool(outcome.all_stages_visited),
        "stages_visited": outcome.stages_visited,
        "stage_events": outcome.stage_events,
        "task_errors": outcome.task_errors,
        "alpha_est": outcome.alpha.tolist(),
        "x_est": outcome.x.tolist(),
        "work": outcome.ledger.snapshot(),
        "run_wall_seconds": outcome.wall_seconds,
    }
    if mismatch_kind is None:
        truth_spatial = None
    else:
        arrs = mismatch_scene_arrays(scene)
        truth_spatial = arrs["truth_spatial"]
    metrics, phase_err, eval_snapshot = evaluate_estimate(
        seed,
        method,
        outcome.alpha,
        outcome.x,
        scene,
        mismatch_kind=mismatch_kind,
        truth_spatial=truth_spatial,
    )
    rec["metrics"] = metrics
    rec["evaluation_ledger"] = eval_snapshot
    if phase_err is not None:
        rec["phase_residual_masked_rad"] = [float(v) for v in phase_err]
    else:
        rec["phase_residual_masked_rad"] = []
    return rec


def _write_batch(records: list[dict[str, Any]]) -> None:
    with RECORDS_PATH.open("a", encoding="utf-8") as fh:
        for rec in records:
            fh.write(_record_line(rec))


def _smoke() -> None:
    """Basic pipeline runs (all four stages but one iteration per stage)."""
    print("[smoke] basic pipeline checks on seed 71 ...")
    settings = default_settings()
    settings["stage_maxiter"] = [1, 1, 1, 1]
    records = []
    for method in ("coherent_fixedpose", "coherent_joint", "intensity_joint"):
        records.append(
            _run_one("smoke", settings, method, TUNING_SEEDS[0], 0)
        )
    _write_batch(records)
    for r in records:
        print(
            f"[smoke] {r['method']:<20} status={r['status']:<18} "
            f"units={r['work']['units_rhs_columns']:>5} "
            f"map_rmse={r['metrics']['spatial_map_rmse']:.4f}"
        )


def _tune() -> None:
    print("[tune] seeds 71-72, radius .125 lambda, all four methods ...")
    settings = default_settings()
    records = []
    for seed in TUNING_SEEDS:
        for method in METHODS:
            records.append(_run_one("tune", settings, method, seed, 0))
    _write_batch(records)
    walls = [r["run_wall_seconds"] for r in records]
    print(
        f"[tune] finished {len(records)} runs; per-run wall s "
        f"min/median/max = {min(walls):.2f}/{np.median(walls):.2f}/"
        f"{max(walls):.2f}"
    )
    frozen = write_frozen(settings)
    print(f"[tune] frozen settings -> {DIR / 'frozen.json'}")
    print(f"[tune] sha256 = {frozen['sha256']}")


def _test() -> None:
    settings = load_frozen()
    print("[test] seeds 801-806, radii .125 and .5 lambda, four methods ...")
    count = 0
    for seed in TEST_SEEDS:
        for radius_index in (0, 1):
            for method in METHODS:
                _save_record(_run_one("test", settings, method, seed, radius_index))
                count += 1
                print(
                    f"[test] {count:>3}/48 seed={seed} radius={radius_index} "
                    f"method={method}"
                )
    print("[test] 48 runs recorded")


def _mismatch() -> None:
    settings = load_frozen()
    kinds = ("clock_phase", "receiver_coupling", "outside_basis")
    print("[mismatch] seed 881, three mismatch controls x two coherent methods")
    count = 0
    for kind in kinds:
        for method in ("coherent_fixedpose", "coherent_joint"):
            _save_record(
                _run_one(
                    "mismatch",
                    settings,
                    method,
                    MISMATCH_SEED,
                    mismatch_kind=kind,
                )
            )
            count += 1
            print(f"[mismatch] {count}/6 kind={kind} method={method}")


def _summarize() -> None:
    from summarize import build_all

    build_all()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("smoke", "tune", "test", "mismatch",
                                       "summarize", "all"))
    args = ap.parse_args()
    if args.mode in ("smoke", "tune", "all"):
        if not RECORDS_PATH.exists():
            RECORDS_PATH.write_text("", encoding="utf-8")
    if args.mode in ("smoke", "all"):
        _smoke()
    if args.mode in ("tune", "all"):
        _tune()
    if args.mode in ("test", "all"):
        _test()
    if args.mode in ("mismatch", "all"):
        _mismatch()
    if args.mode in ("summarize", "all"):
        _summarize()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
