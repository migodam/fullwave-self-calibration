#!/usr/bin/env python3
"""A2 E4 runner: tune | final | summarize.

Strict boundary: final mode requires an explicit --frozen-config path and
never accepts implicit tuning; tuning never touches seeds 1001-1020.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOLVER_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SOLVER_DIR))
sys.path.insert(0, str(ROOT / "research/delegated/a2_physics"))
sys.path.insert(0, str(ROOT / "research/trispace_self_calibration/a2_research"))

from common import (  # noqa: E402
    Calibration,
    CostLedger,
    FINAL_SEED_RANGE,
    TUNING_SEED_RANGE,
    default_settings,
    frozen_timing_path,
    jsonable,
    read_json,
    settings_hash,
    write_json,
)
from calibrate import run_microbenchmark  # noqa: E402
from physics import Config, Model  # noqa: E402


def make_models():
    m16f = Model(Config(N=16, aperture="full"))
    m16l = Model(Config(N=16, aperture="limited"))
    m32f = Model(Config(N=32, aperture="full"))
    m32l = Model(Config(N=32, aperture="limited"))
    return {
        ("full", 16): m16f,
        ("limited", 16): m16l,
        ("full", 32): m32f,
        ("limited", 32): m32l,
    }


def candidates_for(method: str, index: int, full_grid: bool = False) -> list[dict]:
    from common import FIXED_RANKS

    maxls_list = [10, 20]
    ftol_list = [1e-8, 1e-10]
    opt = []
    for ml in maxls_list:
        for ft in ftol_list:
            opt.append({"maxls": ml, "ftol": ft, "gtol": 1e-6})
    if method == "fixed_rank" and not full_grid:
        # equal candidate count per physical method: cycle ranks over the
        # same four optimizer settings (all three ranks are exercised).
        return [
            {**o, "fixed_rank": int(FIXED_RANKS[(index + j) % len(FIXED_RANKS)])}
            for j, o in enumerate(opt)
        ]
    if method == "fixed_rank":
        return [
            {**o, "fixed_rank": int(r)}
            for o in opt
            for r in FIXED_RANKS
        ]
    return list(opt)


def parse_seeds(spec: str | None, mode: str) -> list[int]:
    if spec is None:
        if mode == "final":
            return list(FINAL_SEED_RANGE)
        return list(TUNING_SEED_RANGE)
    out = []
    for part in str(spec).split(","):
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    if mode != "final":
        bad = [s for s in out if s in FINAL_SEED_RANGE]
        if bad:
            raise SystemExit(f"refusing final seed access in non-final mode: {bad}")
    if mode == "final":
        bad = [s for s in out if s not in FINAL_SEED_RANGE]
        if bad:
            raise SystemExit(f"final seeds must be in 1001-1020, got {bad}")
    return out


def init_indices(spec: str | None, seeds: list[int], mode: str) -> list[int]:
    if mode == "final" and spec is None:
        return list(range(12))
    if spec is None or spec == "by_seed":
        return "by_seed"
    return [int(x) for x in str(spec).split(",")]


def init_indices_for_seed(
    spec: str | None, seed: int, mode: str
) -> list[int]:
    if mode == "final" and spec is None:
        return list(range(12))
    if spec is None or spec == "by_seed":
        # tuning cost control: one representative init per seed
        return [int(seed) % 12]
    return [int(x) for x in str(spec).split(",")]


def prepare_frozen(args, calibration, reference, tuned_settings, sigma_by_ap):
    frozen_dir = SOLVER_DIR / "frozen"
    frozen_dir.mkdir(exist_ok=True)
    name = f"frozen_{args.tag or datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
    path = frozen_dir / name
    payload = {
        "boundary": "final_mode_only_explicit_frozen_config",
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "budget": args.budget,
        "methods": args.methods,
        "tuned_settings_by_method": tuned_settings,
        "calibration": calibration.to_dict(),
        "reference_design": reference,
        "sigma_by_aperture": sigma_by_ap,
        "fractions": [0.12, 0.28, 0.5, 1.0],
        "frequency_stages": [[0], [0, 1], [0, 1, 2], [0, 1, 2, 3]],
        "init_protocol": "12_offsets_half_squared_lever_metric",
    }
    write_json(path, payload)
    return path


def parse_methods(spec: str | None) -> list[str]:
    if spec:
        out = [m.strip() for m in spec.split(",")]
    else:
        out = ["direct", "prasc", "fixed_rank", "phaseless", "direct_control"]
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("tune", "final", "summarize"), required=True)
    p.add_argument("--budget", type=int, choices=(200, 800), required=True)
    p.add_argument("--methods", default=None, help="comma-separated method names")
    p.add_argument("--seeds", default=None)
    p.add_argument("--inits", default=None)
    p.add_argument("--frozen-config", default=None)
    p.add_argument("--tag", default=None)
    p.add_argument("--grid", choices=("equal4", "exhaustive"), default="equal4")
    p.add_argument("--records", default=None)
    p.add_argument("--init-max-units", type=float, default=32.0)
    args = p.parse_args(argv)
    if args.mode == "final" and not args.frozen_config:
        p.error("--mode final requires an explicit --frozen-config path")

    models = make_models()
    methods = parse_methods(args.methods)
    seeds = parse_seeds(args.seeds, args.mode)

    if args.mode == "summarize":
        records_path = Path(args.records or (SOLVER_DIR / "records_tune.jsonl"))
        summarize(records_path, out_dir=SOLVER_DIR, methods=methods)
        return 0

    if args.mode == "tune":
        return run_tune(args, models, methods, seeds)
    return run_final(args, models, methods, seeds)


def run_tune(args, models, methods, seeds):
    cal_path = SOLVER_DIR / "calibration_tune.json"
    cal = (
        run_microbenchmark(
            models[("full", 16)], reps=4, out_path=cal_path
        )
        if not cal_path.exists()
        else Calibration_from_path(cal_path)
    )
    write_json(cal_path, cal.to_dict())

    # reference design and noise are declared before tuning data generation
    sigma_by_ap = {}
    from generate import nominal_sigma

    for ap in ("full", "limited"):
        sigma_by_ap[ap] = nominal_sigma(models[(ap, 32)])
    from reference import predeclare_reference

    ref = predeclare_reference(
        models[("full", 16)],
        models[("limited", 16)],
        models[("full", 32)],
        models[("limited", 32)],
        sigma_by_ap,
        SOLVER_DIR / "reference_design.json",
    )
    write_json(
        SOLVER_DIR / "noise_sigma_tune.json",
        {"sigma_by_aperture": sigma_by_ap, "snr_db": 30.0, "reference": "nominal_alpha_0.5_total_field_N32"},
    )

    _ = init_indices(args.inits, seeds, "tune")
    records_path = SOLVER_DIR / (
        f"records_tune_{args.tag}.jsonl"
        if args.tag
        else "records_tune.jsonl"
    )
    full_grid = args.grid == "exhaustive"
    # fixed settings for all non-optimizer choices
    base = default_settings()
    base["init_max_units"] = args.init_max_units
    with records_path.open("a", encoding="utf-8") as fh:
        for seed in seeds:
            inits = init_indices_for_seed(args.inits, seed, "tune")
            scene, init_info = prepare_scene_impl(args, models, seed, inits)
            for idx, info in zip(inits, init_info):
                for method in methods:
                    cands = candidates_for(method, int(seed) % 10, full_grid)
                    for cand in cands:
                        st = {**base, **cand}
                        rec = run_one(
                            args, models, seed, scene, info, method, st, ref
                        )
                        fh.write(json.dumps(jsonable(rec), separators=(",", ":")) + "\n")
                        fh.flush()
    print(f"tune records -> {records_path}")
    summarize(records_path, out_dir=SOLVER_DIR, methods=methods)
    best = choose_tuned(records_path, methods)
    frozen_path = prepare_frozen(
        args, cal, ref, best, sigma_by_ap
    )
    print(f"frozen config -> {frozen_path}")
    return 0


def Calibration_from_path(path):
    from common import Calibration

    return Calibration.from_dict(read_json(path))


def prepare_scene_impl(args, models, seed, inits):
    """Scene object shared across the twelve initialisation runs of a seed."""
    from generate import generate_seed_scene, init_offsets, stratum_for_seed

    aperture = stratum_for_seed(seed)[1]
    sigma_by_ap = read_json(SOLVER_DIR / "noise_sigma_tune.json")["sigma_by_aperture"]
    sigma = float(sigma_by_ap[aperture])
    m32 = models[(aperture, 32)]
    scene = generate_seed_scene(seed, m32, sigma)
    offsets = init_offsets(seed)
    infos = []
    for ii in inits:
        off = offsets[ii]
        infos.append({"index": ii, **off})
    return scene, infos


def run_one(
    args, models, seed, scene, init_info, method, settings, ref, frozen=None
):
    from generate import (
        evaluate_estimate,
        initial_material_fit,
    )
    from solver import solve

    aperture = scene["aperture"]
    m16 = models[(aperture, 16)]
    m32 = models[(aperture, 32)]
    y = np.asarray(scene["y"], dtype=np.complex128)
    y_val = np.asarray(scene["y_val"], dtype=np.complex128)
    sigma = float(scene["sigma"])
    x0 = np.asarray(init_info["x0"], dtype=float)
    y_freq0 = y[:72]
    fit = initial_material_fit(
        m16,
        y_freq0,
        sigma,
        x0,
        float(settings.get("init_max_units", 32.0)),
    )
    st = {**settings}
    st["init_units"] = fit["units"]
    alpha0 = np.asarray(fit["alpha0"], dtype=float)
    if frozen is not None:
        m16._calibration = Calibration.from_dict(frozen["calibration"])
    else:
        m16._calibration = Calibration_from_path(
            SOLVER_DIR / "calibration_tune.json"
        )
    t_run = time.perf_counter()
    res = None
    err = None
    try:
        res = solve(m16, y, sigma, alpha0, x0, method, args.budget, st)
    except Exception as exc:  # robustness: a failing run is a record, not a crash
        err = f"{type(exc).__name__}: {exc}"
    wall = time.perf_counter() - t_run
    if res is None:
        res = {
            "method": method,
            "alpha_est": alpha0.tolist(),
            "x_est": x0.tolist(),
            "failure": err or "solver_exception",
            "status": err or "solver_exception",
            "ledger": {"units": fit["units"], "wall_seconds": wall},
            "rank_trajectory": [],
            "cert_history": [],
            "fallback_events": [],
            "incomplete_stages": [],
            "mode": "unknown",
            "reduced_moves": 0,
            "exact_checks": 0,
        }
    margins = ref[aperture]
    Q = np.asarray(margins["Q_res"], dtype=float)
    eval_res = evaluate_estimate(
        m16,
        m32,
        np.asarray(res["alpha_est"]),
        np.asarray(res["x_est"]),
        y,
        y_val,
        sigma,
        np.asarray(scene["alpha_true"]),
        np.asarray(scene["x_true"]),
        Q,
        margins,
    )
    # calibration flags: the covariance surrogate is diagnostic only and is
    # computed here from the exact final Jacobian for the run record.
    from cert import covariance_diagnostic

    cov_eval_ledger = CostLedger()
    fwj = m16.forward(np.asarray(res["alpha_est"]), np.asarray(res["x_est"]), jacobian=True)
    cov_eval_ledger.charge_model_work(fwj["work"])
    cov_diag = covariance_diagnostic(fwj["A"], fwj["B"], sigma)
    cov_pass = bool(cov_diag.get("available") and cov_diag.get("surrogate_pass"))
    failure = res.get("failure")
    status = res.get("status") or (
        f"solver_exception" if failure else "completed_within_budget"
    )
    cov_ok = cov_pass if cov_pass is not None else False
    success = bool(
        failure is None
        and status == "completed_within_budget"
        and eval_res["pose_ok"]
        and eval_res["map_ok"]
        and eval_res["validation_ok"]
    )
    false_calibrated = bool(cov_ok and not success)
    rec = {
        "seed": seed,
        "stratum_strength": scene["stratum_strength"],
        "aperture": aperture,
        "init_index": init_info["index"],
        "init_radius_m": init_info["radius_m"],
        "init_angle": init_info["angle"],
        "x0": x0.tolist(),
        "method": method,
        "config_hash": res.get("config_hash"),
        "budget": args.budget,
        "alpha_truth": scene["alpha_true"],
        "x_truth": scene["x_true"],
        "alpha0": alpha0.tolist(),
        "init_fallback": fit["fallback"],
        "init_reason": fit["reason"],
        "alpha_est": res["alpha_est"],
        "x_est": res["x_est"],
        "failure": failure,
        "status": status,
        "success": success,
        "mode": res.get("mode"),
        "rank": res.get("rank"),
        "r_num": (res.get("rank_trajectory") or [{}])[-1].get("rank", res.get("rank")),
        "L_det": (res.get("rank_trajectory") or [{}])[-1].get("ldet"),
        "r_free": 0,
        "reduced_moves": res.get("reduced_moves"),
        "exact_checks": res.get("exact_checks"),
        "covariance_surrogate_pass": cov_pass,
        "covariance_diagnostic": cov_diag,
        "false_calibrated": false_calibrated,
        "ledger": res.get("ledger"),
        "rank_trajectory": res.get("rank_trajectory", []),
        "cert_history": res.get("cert_history", []),
        "fallback_events": res.get("fallback_events", []),
        "incomplete_stages": res.get("incomplete_stages", []),
        "settings": st,
        "metrics": eval_res,
        "wall_seconds": wall,
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    return rec


def summarize(records_path: Path, out_dir: Path, methods: list[str] | None):
    """Method-level summary; adds paired seed-cluster bootstrap when final
    records (seeds 1001-1020) are present."""
    rows = []
    with Path(records_path).open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    if not rows:
        print("no records")
        return
    methods = methods or sorted({r["method"] for r in rows})
    out = {"records": len(rows), "methods": methods}
    table = []
    for m in methods:
        rs = [r for r in rows if r["method"] == m]
        succ = [r for r in rs if r.get("success")]
        table.append(
            {
                "method": m,
                "runs": len(rs),
                "success_rate": len(succ) / max(1, len(rs)),
                "mean_pose_rmse": float(np.mean([r["metrics"]["pose_error_lever_m"] for r in rs])),
                "mean_task_map_rmse": float(np.mean([r["metrics"]["task_map_rmse"] for r in rs])),
                "mean_full_map_rmse": float(np.mean([r["metrics"]["full_map_rmse"] for r in rs])),
                "mean_wall_s": float(np.mean([r["wall_seconds"] for r in rs])),
                "failure_runs": sum(1 for r in rs if r.get("failure")),
            }
        )
    out["summary_table"] = table
    write_json(out_dir / "summary.json", out)
    if any(r.get("seed", 0) in FINAL_SEED_RANGE for r in rows):
        final_stats = final_statistics(rows)
        write_json(out_dir / "final_statistics.json", final_stats)
        lines.append("")
        lines.append("## Final-seed paired statistics")
        for row in final_stats.get("tables", []):
            lines.append("")
            lines.append(row)
    lines = [
        "# A2 E4 tuning summary",
        "",
        f"- records: {len(rows)}",
        "",
        "| method | runs | success | pose RMSE m | task map RMSE | full map RMSE | wall s | failures |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for t in table:
        lines.append(
            f"| {t['method']} | {t['runs']} | {t['success_rate']:.3f} | "
            f"{t['mean_pose_rmse']:.4g} | {t['mean_task_map_rmse']:.4g} | "
            f"{t['mean_full_map_rmse']:.4g} | {t['mean_wall_s']:.3g} | {t['failure_runs']} |"
        )
    (out_dir / "tuning_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def final_statistics(rows: list[dict]) -> dict[str, Any]:
    """Seed-cluster paired bootstrap for the prasc-vs-direct compound gate."""
    import numpy as np

    prasc = [r for r in rows if r["method"] == "prasc"]
    direct = [r for r in rows if r["method"] == "direct"]
    seeds = sorted({int(r["seed"]) for r in rows if int(r["seed"]) in FINAL_SEED_RANGE})
    # pair by seed and init index
    by_key = {
        (int(r["seed"]), int(r["init_index"]), r["method"]): r
        for r in rows
        if r["method"] in ("prasc", "direct")
    }
    paired_keys = [
        k
        for k in by_key
        if (k[0], k[1], "prasc") in by_key and (k[0], k[1], "direct") in by_key
    ]
    stats = {
        "clusters": len(seeds),
        "paired_runs": len(paired_keys),
        "tables": [],
        "notes": "seed-cluster bootstrap over 10000 draws",
    }
    if not paired_keys:
        return stats
    diff_success = np.empty(len(seeds), dtype=float)
    diff_pose = np.empty(len(seeds), dtype=float)
    diff_map = np.empty(len(seeds), dtype=float)
    for si, seed in enumerate(seeds):
        ks = [k for k in paired_keys if k[0] == seed]
        ds = [
            int(by_key[(seed, ii, "prasc")]["success"])
            - int(by_key[(seed, ii, "direct")]["success"])
            for (_, ii, _) in ks
        ]
        po = [
            by_key[(seed, ii, "prasc")]["metrics"]["pose_error_lever_m"]
            - by_key[(seed, ii, "direct")]["metrics"]["pose_error_lever_m"]
            for (_, ii, _) in ks
        ]
        ma = [
            by_key[(seed, ii, "prasc")]["metrics"]["task_map_rmse"]
            - by_key[(seed, ii, "direct")]["metrics"]["task_map_rmse"]
            for (_, ii, _) in ks
        ]
        diff_success[si] = float(np.mean(ds))
        diff_pose[si] = float(np.mean(po))
        diff_map[si] = float(np.mean(ma))
    rng = np.random.default_rng(20260905)
    n = len(seeds)
    draws = 10000
    idx = rng.integers(0, n, size=(draws, n))
    boot_s = diff_success[idx].mean(axis=1)
    boot_p = diff_pose[idx].mean(axis=1)
    boot_m = diff_map[idx].mean(axis=1)
    stats["success_diff_mean"] = float(diff_success.mean())
    stats["pose_rmse_diff_mean_m"] = float(diff_pose.mean())
    stats["map_rmse_diff_mean"] = float(diff_map.mean())
    stats["success_diff_95_ci"] = list(np.percentile(boot_s, [2.5, 97.5]))
    stats["pose_diff_95_ci_m"] = list(np.percentile(boot_p, [2.5, 97.5]))
    stats["map_diff_95_ci"] = list(np.percentile(boot_m, [2.5, 97.5]))
    stats["tables"] = [
        "mean(success prasc - direct) = "
        f"{diff_success.mean():.4f} (95% CI {np.percentile(boot_s,[2.5,97.5]).tolist()})",
        "mean(pose RMSE prasc - direct, m) = "
        f"{diff_pose.mean():.5f} (95% CI {np.percentile(boot_p,[2.5,97.5]).tolist()})",
        "mean(task-map RMSE prasc - direct) = "
        f"{diff_map.mean():.5f} (95% CI {np.percentile(boot_m,[2.5,97.5]).tolist()})",
    ]
    return stats


def choose_tuned(records_path: Path, methods: list[str]) -> dict[str, dict]:
    """Pick one locked setting per physical method from tuning records only."""
    rows = []
    with Path(records_path).open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    best: dict[str, dict] = {}
    for m in methods:
        groups: dict[str, list] = {}
        for r in rows:
            if r["method"] != m:
                continue
            key = json.dumps(r["settings"], sort_keys=True, default=str)
            groups.setdefault(key, []).append(r)
        if not groups:
            continue
        ranked = []
        for key, rs in groups.items():
            sr = np.mean([1.0 if x.get("success") else 0.0 for x in rs])
            loss = np.mean([x["metrics"]["loss_exact_final"] for x in rs])
            pose = np.mean([x["metrics"]["pose_error_lever_m"] for x in rs])
            ranked.append(((-sr, loss, pose), key, rs[0]["settings"]))
        ranked.sort(key=lambda t: t[0])
        best[m] = ranked[0][2]
    return best


def run_final(args, models, methods, seeds):
    frozen = read_json(Path(args.frozen_config))
    if int(frozen.get("budget")) != args.budget:
        raise SystemExit("frozen budget does not match --budget")
    if set(frozen["methods"]) != set(methods):
        raise SystemExit("frozen methods do not match --methods")
    write_json(
        frozen_timing_path(Path(args.frozen_config).parent),
        frozen["calibration"],
    )
    sigma_by_ap = frozen["sigma_by_aperture"]
    write_json(
        SOLVER_DIR / "noise_sigma_final.json",
        {"sigma_by_aperture": sigma_by_ap, "locked": True},
    )
    settings_by_method = {
        m: {**default_settings(), **frozen["tuned_settings_by_method"].get(m, {})}
        for m in methods
    }
    records_path = SOLVER_DIR / (
        f"records_final_{args.tag}.jsonl" if args.tag else "records_final.jsonl"
    )
    # per-run scenes are generated once per seed
    scene_by_seed = {}
    init_info_by_seed = {}
    for seed in seeds:
        from generate import generate_seed_scene, init_offsets, stratum_for_seed

        inits = init_indices_for_seed(args.inits, seed, "final")
        ap = stratum_for_seed(seed)[1]
        m32 = models[(ap, 32)]
        scene_by_seed[seed] = generate_seed_scene(
            seed, m32, float(sigma_by_ap[ap]), allow_final=True
        )
        offsets = init_offsets(seed)
        init_info_by_seed[seed] = [{"index": i, **offsets[i]} for i in inits]
    ref = frozen["reference_design"]
    with records_path.open("a", encoding="utf-8") as fh:
        for seed in seeds:
            scene = scene_by_seed[seed]
            for info in init_info_by_seed[seed]:
                for m in methods:
                    rec = run_one(
                        args,
                        models,
                        seed,
                        scene,
                        info,
                        m,
                        dict(settings_by_method[m]),
                        ref,
                        frozen,
                    )
                    fh.write(json.dumps(jsonable(rec), separators=(",", ":")) + "\n")
                    fh.flush()
    print(f"final records -> {records_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
