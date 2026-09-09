"""Aggregate summaries and Markdown for the a2_highdim worker."""

from __future__ import annotations

import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from common import DIR
from figures import make_all_figures

RECORDS_PATH = DIR / "records.jsonl"
METHOD_LABELS = {
    "coherent_fixedpose": "coherent_fixedpose",
    "coherent_joint": "coherent_joint",
    "intensity_joint": "intensity_joint",
    "oracle": "oracle",
}


def load_records() -> list[dict[str, Any]]:
    recs = []
    if not RECORDS_PATH.exists():
        return recs
    with RECORDS_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                recs.append(json.loads(line))
    return recs


def _test_map(test: list[dict[str, Any]]) -> dict[tuple[int, int, str], dict[str, Any]]:
    return {
        (int(r["seed"]), int(r["radius_index"]), str(r["method"])): r
        for r in test
    }


def _cluster_bootstrap_diff(
    test: list[dict[str, Any]],
    method_a: str,
    method_b: str,
    metric_key: str,
    b: int = 10_000,
    seed: int = 20260906,
) -> dict[str, Any]:
    mapped = _test_map(test)
    seeds = sorted({int(r["seed"]) for r in test})
    radii = sorted({int(r["radius_index"]) for r in test})
    diffs = []
    for s in seeds:
        for ri in radii:
            ra = mapped[(s, ri, method_a)]
            rb = mapped[(s, ri, method_b)]
            diffs.append(
                float(ra["metrics"][metric_key])
                - float(rb["metrics"][metric_key])
            )
    diffs = np.asarray(diffs, dtype=float)
    if len(diffs) < 2:
        return {"n": int(len(diffs)), "mean_diff": None, "ci": None}
    mean_a = float(np.mean([
        mapped[(s, ri, method_a)]["metrics"][metric_key]
        for s in seeds for ri in radii
    ]))
    mean_b = float(np.mean([
        mapped[(s, ri, method_b)]["metrics"][metric_key]
        for s in seeds for ri in radii
    ]))
    rng = np.random.default_rng(seed)
    resample = np.empty(b, dtype=float)
    for i in range(b):
        chosen = rng.choice(seeds, size=len(seeds), replace=True)
        vals_a = [
            mapped[(s, ri, method_a)]["metrics"][metric_key]
            for s in chosen for ri in radii
        ]
        vals_b = [
            mapped[(s, ri, method_b)]["metrics"][metric_key]
            for s in chosen for ri in radii
        ]
        resample[i] = float(np.mean(vals_a) - np.mean(vals_b))
    lo, hi = np.percentile(resample, [2.5, 97.5])
    return {
        "n": int(len(diffs)),
        "mean_a": mean_a,
        "mean_b": mean_b,
        "mean_diff_a_minus_b": float(np.mean(diffs)),
        "cluster_bootstrap_95ci": [float(lo), float(hi)],
        "bootstrap_n": b,
        "clusters": len(seeds),
        "note": "descriptive 6-seed cluster bootstrap; exploratory, not "
                "confirmatory",
    }


def _aggregate_stats(test: list[dict[str, Any]], metric_key: str) -> dict[str, Any]:
    out = {}
    for method in sorted({r["method"] for r in test}):
        vals = np.asarray([
            r["metrics"][metric_key] for r in test if r["method"] == method
            and r["metrics"].get(metric_key) is not None
        ], dtype=float)
        out[method] = {
            "n": int(vals.size),
            "mean": float(np.mean(vals)) if vals.size else None,
            "median": float(np.median(vals)) if vals.size else None,
            "std": float(np.std(vals, ddof=1)) if vals.size > 1 else None,
            "min": float(np.min(vals)) if vals.size else None,
            "max": float(np.max(vals)) if vals.size else None,
        }
    return out


def _phase_stats(test: list[dict[str, Any]]) -> dict[str, Any]:
    out = {}
    for method in sorted({r["method"] for r in test}):
        samples = np.asarray([
            v for r in test if r["method"] == method
            for v in r.get("phase_residual_masked_rad", [])
        ], dtype=float)
        if samples.size == 0:
            out[method] = {"n": 0}
            continue
        from scipy.stats import circstd

        out[method] = {
            "masked_channel_samples": int(samples.size),
            "mean_abs_rad": float(np.mean(np.abs(samples))),
            "median_abs_rad": float(np.median(np.abs(samples))),
            "rms_rad": float(np.sqrt(np.mean(samples**2))),
            "circular_std_rad": float(
                circstd(samples, high=np.pi, low=-np.pi)
            ),
        }
    return out


def _summary_markdown(summary: dict[str, Any]) -> str:
    lines = []
    lines.append("# a2_highdim exploratory summary")
    lines.append("")
    lines.append(f"- Generated: {summary['generated_at']}")
    lines.append(f"- Frozen settings sha256: {summary['settings_sha256']}")
    lines.append("")
    lines.append("## Scope and epistemic status")
    lines.append("")
    lines.append(summary["scope_statement"])
    lines.append("")
    lines.append("## Run counts and statuses")
    lines.append("")
    lines.append("| phase | method | n | statuses |")
    lines.append("|---|---|---|---|")
    for row in summary["counts_table"]:
        lines.append(
            f"| {row['phase']} | {row['method']} | {row['n']} | "
            f"{row['statuses']} |"
        )
    lines.append("")
    lines.append("## Unconditional test metrics (seeds 801-806, radii "
                 ".125/.5 lambda, n=12 per method)")
    lines.append("")
    lines.append("### Spatial map RMSE")
    lines.append("")
    lines.append("| method | n | mean | median | std | min | max |")
    lines.append("|---|---|---|---|---|---|---|")
    for method, st in summary["spatial_rmse"].items():
        lines.append(_row(method, st))
    lines.append("")
    lines.append("### Coefficient RMSE (49 real coefficients)")
    lines.append("")
    lines.append("| method | n | mean | median | std | min | max |")
    lines.append("|---|---|---|---|---|---|---|")
    for method, st in summary["coeff_rmse"].items():
        lines.append(_row(method, st))
    lines.append("")
    lines.append("### Pose lever-arm error (m)")
    lines.append("")
    lines.append("| method | n | mean | median | std | min | max |")
    lines.append("|---|---|---|---|---|---|---|")
    for method, st in summary["pose_error"].items():
        lines.append(_row(method, st))
    lines.append("")
    lines.append("### Masked forward-phase residual (vs clean truth total; "
                 "|truth|>3 sigma, evaluation-only)")
    lines.append("")
    lines.append("| method | samples | mean abs (rad) | median abs (rad) | "
                 "rms (rad) | circular std (rad) |")
    lines.append("|---|---|---|---|---|---|")
    for method, st in summary["phase_stats"].items():
        if "mean_abs_rad" not in st:
            lines.append(f"| {method} | 0 | - | - | - | - |")
            continue
        lines.append(
            f"| {method} | {st['masked_channel_samples']} | "
            f"{st['mean_abs_rad']:.4f} | {st['median_abs_rad']:.4f} | "
            f"{st['rms_rad']:.4f} | {st['circular_std_rad']:.4f} |"
        )
    lines.append("")
    lines.append("## Paired differences, 6-seed cluster bootstrap "
                 "(descriptive 95% CIs)")
    lines.append("")
    lines.append("| metric | comparison | mean difference | 95% CI |")
    lines.append("|---|---|---|---|")
    for row in summary["paired_differences"]:
        ci = row["ci"]
        if ci is None:
            ci_txt = "-"
        else:
            ci_txt = f"[{ci[0]:.4f}, {ci[1]:.4f}]"
        lines.append(
            f"| {row['metric']} | {row['comparison']} | "
            f"{row['mean_diff']:.4f} | {ci_txt} |"
        )
    lines.append("")
    lines.append("Bootstrap convention: resample the six seed clusters "
                 "(both pose-error radii travel together), 10,000 draws; "
                 "percentile 2.5-97.5. Exploratory descriptive intervals, "
                 "not a confirmatory gate.")
    lines.append("")
    lines.append("## Mismatch controls (seed 881)")
    lines.append("")
    lines.append("| kind | method | spatial map RMSE | pose err (m) | "
                 "fit chi2 p | val chi2 p | false-assurance flag |")
    lines.append("|---|---|---|---|---|---|---|")
    for row in summary["mismatch_rows"]:
        lines.append(
            f"| {row['kind']} | {row['method']} | {row['map_rmse']:.4f} | "
            f"{row['pose_err']:.4f} | {row['fit_p']:.3e} | "
            f"{row['val_p']:.3e} | {row['false_assurance']} |"
        )
    lines.append("")
    lines.append(summary["mismatch_statement"])
    lines.append("")
    lines.append("## Figures")
    lines.append("")
    lines.append("- `figures/reconstruction_map.png`")
    lines.append("- `figures/error_distributions.png`")
    lines.append("- `figures/phase_residuals.png`")
    lines.append("- `figures/mismatch_residuals.png`")
    lines.append("")
    lines.append("## Work/wall accounting")
    lines.append("")
    lines.append(summary["work_statement"])
    lines.append("")
    lines.append("| phase | n runs | median solver wall (s) | total solver wall (s) | "
                 "median eval wall (s) | total eval wall (s) |")
    lines.append("|---|---|---|---|---|---|")
    for row in summary["wall_table"]:
        lines.append(
            f"| {row['phase']} | {row['n']} | {row['median_solver_s']:.2f} | "
            f"{row['total_solver_s']:.1f} | {row['median_eval_s']:.2f} | "
            f"{row['total_eval_s']:.1f} |"
        )
    lines.append("")
    lines.append("## Exact commands")
    lines.append("")
    lines.append("```bash")
    lines.append(summary["commands"])
    lines.append("```")
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    for lim in summary["limitations"]:
        lines.append(f"- {lim}")
    return "\n".join(lines)


def _row(method: str, st: dict[str, Any]) -> str:
    fmt = lambda v: "None" if v is None else f"{v:.4f}"
    return (
        f"| {method} | {st['n']} | {fmt(st['mean'])} | {fmt(st['median'])} | "
        f"{fmt(st['std'])} | {fmt(st['min'])} | {fmt(st['max'])} |"
    )


def build_all() -> None:
    records = load_records()
    test = [r for r in records if r.get("phase") == "test"]
    tune = [r for r in records if r.get("phase") == "tune"]
    smoke = [r for r in records if r.get("phase") == "smoke"]
    mismatch = [r for r in records if r.get("phase") == "mismatch"]
    counts_rows = []
    for phase, recs in (("tune", tune), ("test", test), ("mismatch", mismatch),
                        ("smoke", smoke)):
        counts = Counter((r["method"], r["status"]) for r in recs)
        methods = sorted({r["method"] for r in recs})
        for method in methods:
            statuses = counts.get((method, "schedule_complete"), 0)
            budget = counts.get((method, "budget_stopped"), 0)
            errs = counts.get((method, "task_error"), 0)
            time_stops = counts.get((method, "time_stopped"), 0)
            counts_rows.append(
                {
                    "phase": phase,
                    "method": method,
                    "n": int(sum(
                        v for (m, _), v in counts.items() if m == method
                    )),
                    "statuses": (
                        f"complete={statuses}, budget_stopped={budget}, "
                        f"task_error={errs}, time_stopped={time_stops}"
                    ),
                }
            )

    paired = []
    comparisons = [
        ("coherent_joint", "coherent_fixedpose"),
        ("intensity_joint", "coherent_fixedpose"),
        ("intensity_joint", "coherent_joint"),
        ("oracle", "coherent_joint"),
    ]
    for a, b in comparisons:
        for key, label in (
            ("spatial_map_rmse", "spatial map RMSE"),
            ("coefficient_rmse", "coefficient RMSE"),
            ("pose_error_lever_m", "pose lever-arm error"),
        ):
            if key == "coefficient_rmse" and (a == "oracle" or b == "oracle"):
                # oracle has the same coefficient metric; keep it.
                pass
            row = _cluster_bootstrap_diff(test, a, b, key)
            if row.get("mean_diff_a_minus_b") is None:
                continue
            ci = row.get("cluster_bootstrap_95ci")
            paired.append({
                "metric": label,
                "comparison": f"{a} minus {b}",
                "mean_diff": row["mean_diff_a_minus_b"],
                "ci": ci,
                "cluster_info": {
                    "n": row["n"],
                    "clusters": row["clusters"],
                    "bootstrap_n": row["bootstrap_n"],
                },
            })

    mismatch_rows = []
    for r in mismatch:
        kind = r.get("mismatch_kind", "")
        mm = r["metrics"]
        mismatch_rows.append({
            "kind": kind,
            "method": r["method"],
            "map_rmse": float(mm["spatial_map_rmse"]),
            "pose_err": float(mm["pose_error_lever_m"]),
            "fit_p": float(mm["fit_chi2_p"]),
            "val_p": float(mm["validation_chi2_p"]),
            "false_assurance": bool(mm.get("false_assurance_potential", False)),
        })

    frozen_path = DIR / "frozen.json"
    settings_sha = "-"
    if frozen_path.exists():
        settings_sha = json.loads(frozen_path.read_text(encoding="utf-8"))["sha256"]

    wall_table = []
    for phase in ("tune", "test", "mismatch", "smoke"):
        recs = [r for r in records if r.get("phase") == phase]
        if not recs:
            continue
        sw = np.asarray([r.get("run_wall_seconds", 0.0) for r in recs])
        ew = np.asarray([
            r.get("evaluation_ledger", {}).get("wall_seconds", 0.0)
            for r in recs
        ])
        wall_table.append({
            "phase": phase,
            "n": len(recs),
            "median_solver_s": float(np.median(sw)),
            "total_solver_s": float(np.sum(sw)),
            "median_eval_s": float(np.median(ew)) if len(ew) else 0.0,
            "total_eval_s": float(np.sum(ew)),
        })

    commands = (
        "VENV=experiments/idea_loops/loop_2026-09-04_02-58-31/"
        "experiment_geometry_lifted_trispace_som/.venv/bin/python\n"
        "# bounded implementation checks\n"
        "OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV "
        "research/delegated/a2_highdim/test_highdim.py\n"
        "# basic one-iteration pipeline smoke (seed 71, pre-freeze settings)\n"
        "OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV "
        "research/delegated/a2_highdim/runner.py --mode smoke\n"
        "# tuning seeds 71-72, then frozen.json written\n"
        "OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV "
        "research/delegated/a2_highdim/runner.py --mode tune\n"
        "# 48 final exploratory runs, seeds 801-806, frozen settings only\n"
        "OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV "
        "research/delegated/a2_highdim/runner.py --mode test\n"
        "# 6 controlled mismatch runs, seed 881\n"
        "OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV "
        "research/delegated/a2_highdim/runner.py --mode mismatch\n"
        "# aggregate summary and figures (reproducible, no reruns)\n"
        "OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV "
        "research/delegated/a2_highdim/runner.py --mode summarize"
    )
    summary: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "settings_sha256": settings_sha,
        "scope_statement": (
            "EXPLORATORY 2D imaging feasibility only; not the A2 E4 primary "
            "endpoint and not publication or production acceptance. "
            "Reconstructions are model-conditional within the 49-Gaussian "
            "basis and the scalar full-wave physics core. The known-true-pose "
            "method is an explicit oracle upper reference and is not a "
            "practical baseline."
        ),
        "counts_table": counts_rows,
        "spatial_rmse": _aggregate_stats(test, "spatial_map_rmse"),
        "coeff_rmse": _aggregate_stats(test, "coefficient_rmse"),
        "pose_error": _aggregate_stats(test, "pose_error_lever_m"),
        "phase_stats": _phase_stats(test),
        "paired_differences": paired,
        "mismatch_rows": mismatch_rows,
        "mismatch_statement": (
            "These controlled mismatches test the coherent geometry-only "
            "model against an unmodelled clock phase, receiver coupling and "
            "an outside-49-basis inclusion. A fitted residual that looks "
            "calibrated under model mismatch is flagged as false-assurance "
            "potential; no geometry proof and no calibration guarantee is "
            "claimed. The fixed-pose result merely fixes geometry and does "
            "not remove mismatch; target artifacts must not be read as "
            "validated localisation."
        ),
        "work_statement": (
            "Per-run RHS solve columns (forward states, adjoints and any "
            "initialisation) are reported in work.units_rhs_columns with "
            "kind counts; LU factorisations, operator products and wall "
            "seconds are separate raw ledger fields and are never converted "
            "into solve-equivalent units. Evaluation-only N32/N20 forward "
            "solves are in evaluation_ledger and never feed estimator "
            "decisions. This report draws no work-advantage conclusion."
        ),
        "wall_table": wall_table,
        "commands": commands,
        "limitations": [
            "Exploratory single aperture (full), one noise level (30 dB) and "
            "six test scenes with two pose-error radii; no final E4 "
            "seeds 1001-1020 were used.",
            "The inverse N20 grid required a runtime-only extension of the "
            "physics Config allowed-N tuple; physics.py itself was not "
            "changed and no formula was altered.",
            "The 49-Gaussian coefficient basis is strongly overlapping "
            "(Gram condition ~2.7e4 at N20), so coefficient RMSE is not an "
            "interpretable spatial measure; spatial N32 map RMSE is the "
            "imaging metric.",
            "Cluster-bootstrap intervals are descriptive only and cannot "
            "establish equivalence or superiority.",
            "The matched-intensity comparison uses the induced Rice "
            "magnitude NLL and exactly the parent |y|; it is not a "
            "direct-intensity sensor model.",
            "Fixed-pose and joint methods share the frozen schedule and "
            "identical settings; oracle shares them too but uses x_true, so "
            "its advantage is not a calibration claim.",
            "Mismatch controls are misspecification tests of the coherent "
            "geometry-only model; an augmented nuisance estimator was not "
            "implemented.",
        ],
    }
    (DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (DIR / "summary.md").write_text(
        _summary_markdown(summary), encoding="utf-8"
    )
    print("[summarize] wrote summary.json and summary.md")
    make_all_figures(records)
