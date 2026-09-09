#!/usr/bin/env python3
"""Compute requested 95% CI statistics from existing A2 result JSONs (read-only).

Reads the five result files listed in the task and writes
results/ci_stats.json and results/ci_stats.txt. It does not rerun any
experiment and does not modify the source JSON files.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def _t_ppf_975(df: float) -> float:
    """Two-sided 97.5th percentile of Student t, scipy first, then manual."""
    try:
        from scipy import stats  # type: ignore

        return float(stats.t.ppf(0.975, df))
    except Exception:
        pass

    # Manual fallback: bisection against a composite-Simpson integral of the
    # Student t density.  The 2.5% right tail lies well below 20 for df >= 1.
    def pdf(x: float) -> float:
        v = float(df)
        c = math.gamma((v + 1.0) / 2.0) / (
            math.sqrt(v * math.pi) * math.gamma(v / 2.0)
        )
        return c * (1.0 + x * x / v) ** (-(v + 1.0) / 2.0)

    def cdf(x: float) -> float:
        a = -20.0
        b = min(max(x, a + 1e-6), 20.0)
        n = 40_000
        h = (b - a) / n
        total = pdf(a) + pdf(b)
        for i in range(1, n):
            total += 4.0 * pdf(a + i * h) if i % 2 == 1 else 2.0 * pdf(a + i * h)
        return (h / 3.0) * total

    lo, hi = 0.0, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < 0.975:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sf4(value: float | int | None):
    """Round a scalar to four significant figures (0 is unchanged)."""
    if value is None:
        return None
    value = float(value)
    if value == 0.0 or not math.isfinite(value):
        return value
    return float(f"{value:.4g}")


def fmt4(value: float | int | None) -> str:
    if value is None:
        return "None"
    return f"{value:.4g}"


def summary(vals, with_ci: bool = True) -> dict:
    """n, mean, sample sd (ddof=1), min/max and 95% t-CI of the mean."""
    vals = [float(v) for v in vals if v is not None]
    n = len(vals)
    mean = sum(vals) / n if n else float("nan")
    if n > 1:
        sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (n - 1))
    else:
        sd = None
    out = {"n": n}
    if n:
        out["mean"] = sf4(mean)
        out["min"] = sf4(min(vals))
        out["max"] = sf4(max(vals))
    if n > 1 and with_ci:
        tcrit = _t_ppf_975(n - 1)
        half = tcrit * sd / math.sqrt(n)
        out["sd"] = sf4(sd)
        out["ci_mean"] = [sf4(mean - half), sf4(mean + half)]
    elif sd is not None:
        out["sd"] = sf4(sd)
    return out


def ci_of_shift(vals, shift: float) -> dict:
    """Statistics of (val - shift) plus the shift source metadata."""
    diffs = [float(v) - shift for v in vals if v is not None]
    s = summary(diffs, with_ci=True)
    s["analytic"] = sf4(shift)
    return s


def collect_missing(missing, label, found):
    if not found:
        missing.append(label)


def load(name: str) -> dict:
    with open(RESULTS / name, encoding="utf-8") as fh:
        return json.load(fh)


def main() -> int:
    missing: list[str] = []
    stats_out: dict = {}

    # ---------------------------------------------------------------- E1
    e1 = load("e1_results.json")
    e1_models = e1.get("models", {})
    e1_out = {}
    skipped = []
    for model, model_dict in e1_models.items():
        per_seed = model_dict.get("per_seed") or []
        if not per_seed:
            skipped.append(model)
            continue
        model_stats = {}
        for emp_key, analytic_key, metric_name in (
            ("j_coh_emp", "j_coh_analytic", "j_coh"),
            ("j_ph_emp", "j_ph_analytic", "j_ph"),
        ):
            vals = [rec.get(emp_key) for rec in per_seed if emp_key in rec]
            analytic = model_dict.get(analytic_key)
            if not vals or analytic is None:
                collect_missing(
                    missing,
                    f"e1.{model}.{emp_key} or {analytic_key}",
                    bool(vals and analytic is not None),
                )
                continue
            entry = summary(vals, with_ci=True)
            entry["analytic"] = sf4(analytic)
            entry["analytic_key"] = analytic_key
            entry["empirical_minus_analytic"] = ci_of_shift(vals, analytic)
            entry["empirical_minus_analytic"]["analytic_key"] = analytic_key
            entry["source_key"] = emp_key
            model_stats[metric_name] = entry
        e1_out[model] = model_stats
    stats_out["e1"] = {
        "source_file": "results/e1_results.json",
        "skipped_models_empty_per_seed": skipped,
        "models": e1_out,
    }

    # ---------------------------------------------------------------- E2
    e2 = load("e2_results.json")
    e2_records = e2.get("results", {}).get("random_tangents", {}).get("records", [])
    if not e2_records:
        collect_missing(missing, "e2.results.random_tangents.records", bool(e2_records))
    e2_out = {}
    for field, want in (
        ("scaled_residual", ("n", "mean", "sd", "min", "max", "ci_mean")),
        ("loss_residual_fro", ("n", "min", "mean", "max")),
    ):
        vals = [rec.get(field) for rec in e2_records if field in rec]
        collect_missing(missing, f"e2.record.{field}", len(vals) == len(e2_records))
        if not vals:
            continue
        e2_out[field] = summary(vals, with_ci=True)
        e2_out[field]["source_path"] = "results.random_tangents.records[].%s" % field
    stats_out["e2"] = {
        "source_file": "results/e2_results.json",
        "source_path": "results.random_tangents.records",
        "n_records": len(e2_records),
        **e2_out,
    }

    # ---------------------------------------------------------------- E3
    e3 = load("e3_results.json")
    e3_records = e3.get("results", {}).get("theorem7", {}).get("mc", {}).get("records", [])
    if not e3_records:
        collect_missing(missing, "e3.results.theorem7.mc.records", bool(e3_records))
    actual_e3 = {
        "variance_rel": "variance_rel_error",
        "bias_rel": "bias_rel_error",
        "sample_cov_vs_Jx_inv_fro_rel": "sample_cov_vs_Jx_inv_fro_rel",
    }
    e3_out = {}
    for output_key, actual_key in actual_e3.items():
        vals = [rec.get(actual_key) for rec in e3_records if actual_key in rec]
        collect_missing(
            missing, f"e3.record.{actual_key}", len(vals) == len(e3_records)
        )
        if not vals:
            continue
        entry = summary(vals, with_ci=True)
        entry["source_key"] = actual_key
        entry["source_path"] = "results.theorem7.mc.records[].%s" % actual_key
        e3_out[output_key] = entry
    stats_out["e3"] = {
        "source_file": "results/e3_results.json",
        "source_path": "results.theorem7.mc.records",
        "n_records": len(e3_records),
        **e3_out,
    }

    # ---------------------------------------------------------------- E5
    e5 = load("e5_results.json")
    e5_records = (
        e5.get("results", {})
        .get("rank_acquisition_budget", {})
        .get("seeds", {})
        .get("records", [])
    )
    if not e5_records:
        collect_missing(
            missing,
            "e5.results.rank_acquisition_budget.seeds.records",
            bool(e5_records),
        )
    e5_out = {}
    id_vals = [rec.get("identity_fro_residual") for rec in e5_records]
    collect_missing(
        missing, "e5.record.identity_fro_residual", len(id_vals) == len(e5_records)
    )
    if e5_records:
        entry = summary(id_vals, with_ci=True)
        entry["count_le_1e-12"] = int(sum(v is not None and v <= 1e-12 for v in id_vals))
        entry["source_path"] = (
            "results.rank_acquisition_budget.seeds.records[].identity_fro_residual"
        )
        e5_out["identity_fro_residual"] = entry

        flag_key = "loewner_I_ge_L"
        flags = [rec.get(flag_key) for rec in e5_records]
        collect_missing(missing, "e5.record.loewner_I_ge_L", all(isinstance(f, bool) for f in flags))
        if all(isinstance(f, bool) for f in flags):
            e5_out["loewner_flag"] = {
                "source_key": flag_key,
                "source_path": (
                    "results.rank_acquisition_budget.seeds.records[].%s" % flag_key
                ),
                "meaning": "flag is the Loewner condition I_acq >= L_rank",
                "true": int(sum(flags)),
                "false": int(len(flags) - sum(flags)),
                "n": int(len(flags)),
            }
    stats_out["e5"] = {
        "source_file": "results/e5_results.json",
        "note": (
            "The 24 requested budget records are stored under "
            "results.rank_acquisition_budget.seeds.records."
        ),
        "source_path": "results.rank_acquisition_budget.seeds.records",
        "n_records": len(e5_records),
        **e5_out,
    }

    # ------------------------------------------------- physical optional
    phys = load("physical_optional_results.json").get("results", {})
    phase1_records = phys.get("phase1", {}).get("per_seed_records", [])
    if len(phase1_records) != 12:
        collect_missing(
            missing, "physical.phase1.per_seed_records (expected 12)", False
        )
    seed_maxima = []
    for rec in phase1_records:
        vals = [p.get("backward_scaled_residual") for p in rec.get("t5_pairs", [])]
        collect_missing(
            missing,
            "physical.phase1.record.t5_pairs[].backward_scaled_residual",
            bool(vals),
        )
        if vals:
            seed_maxima.append(max(vals))
    agg1 = phys.get("phase1", {}).get("aggregate", {})
    phase1_out = {}
    if seed_maxima:
        phase1_out["per_seed_max_t5a_backward_scaled_residual"] = {
            "source_path": (
                "results.phase1.per_seed_records[].t5_pairs[]"
                ".backward_scaled_residual -> max per seed"
            ),
            "n": len(seed_maxima),
            "min": sf4(min(seed_maxima)),
            "mean": sf4(sum(seed_maxima) / len(seed_maxima)),
            "max": sf4(max(seed_maxima)),
            "seed_maxima": [sf4(v) for v in seed_maxima],
        }
    aggregate_vals = {
        "max_t5a_backward_scaled_residual": agg1.get("max_t5a_backward_scaled_residual"),
        "max_theorem3_residual_max_abs": agg1.get("max_theorem3_residual_max_abs"),
    }
    for key, value in aggregate_vals.items():
        if value is None:
            collect_missing(missing, f"physical.phase1.aggregate.{key}", False)
    phase1_out["aggregate"] = {
        "source_path": "results.phase1.aggregate",
        "max_t5a_backward_scaled_residual": sf4(aggregate_vals["max_t5a_backward_scaled_residual"]),
        "max_theorem3_residual_max_abs": sf4(aggregate_vals["max_theorem3_residual_max_abs"]),
    }
    stats_out["physical_optional"] = {
        "source_file": "results/physical_optional_results.json",
        "phase1": phase1_out,
    }

    phase2_records = phys.get("phase2", {}).get("per_seed_records", [])
    if len(phase2_records) != 12:
        collect_missing(missing, "physical.phase2.per_seed_records (expected 12)", False)
    min_eigs = [
        rec.get("stack_vs_sum", {}).get("min_eig_J_stack_minus_sum")
        for rec in phase2_records
    ]
    collect_missing(
        missing,
        "physical.phase2.record.stack_vs_sum.min_eig_J_stack_minus_sum",
        all(v is not None for v in min_eigs) and len(min_eigs) == 12,
    )
    agg2 = phys.get("phase2", {}).get("aggregate", {})
    phase2_out = {}
    if all(v is not None for v in min_eigs) and min_eigs:
        vals = [float(v) for v in min_eigs]
        phase2_out["per_seed_min_eig_J_stack_minus_sum"] = {
            "source_path": (
                "results.phase2.per_seed_records[].stack_vs_sum"
                ".min_eig_J_stack_minus_sum"
            ),
            "n": len(vals),
            "min": sf4(min(vals)),
            "mean": sf4(sum(vals) / len(vals)),
            "max": sf4(max(vals)),
            "per_seed": [sf4(v) for v in vals],
        }
    agg2_out = {}
    for key in ("max_innovation_rel_residual", "max_budget_identity_rel_residual"):
        value = agg2.get(key)
        if value is None:
            collect_missing(missing, f"physical.phase2.aggregate.{key}", False)
        else:
            agg2_out[key] = sf4(value)
    phase2_out["aggregate"] = {
        "source_path": "results.phase2.aggregate",
        **agg2_out,
    }
    stats_out["physical_optional"]["phase2"] = phase2_out

    stats_out["meta"] = {
        "ci_description": (
            "95% two-sided Student-t confidence interval for the mean; "
            "sample sd uses ddof=1"
        ),
        "significant_figures": 4,
        "generated_from_existing_results": True,
    }
    stats_out["missing_keys"] = missing

    json_path = RESULTS / "ci_stats.json"
    json_path.write_text(
        json.dumps(stats_out, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    # ------------------------------------------------------- text layout
    lines: list[str] = []
    lines.append("A2 PRASC theory-validation 95% confidence statistics")
    lines.append("Method: two-sided Student-t CI for the mean; sample sd ddof=1; 4 significant figures.")
    lines.append("")

    for model, model_stats in e1_out.items():
        for metric, entry in model_stats.items():
            emp = entry
            d = emp["empirical_minus_analytic"]
            lines.append(
                "E1 %s %s (n=%d): mean=%.4g sd=%.4g ci95_mean=[%.4g, %.4g] "
                "ci95(emp-%.4g)=[%.4g, %.4g] analytic_key=%s"
                % (
                    model,
                    metric,
                    emp["n"],
                    emp["mean"],
                    emp["sd"],
                    emp["ci_mean"][0],
                    emp["ci_mean"][1],
                    emp["analytic"],
                    d["ci_mean"][0],
                    d["ci_mean"][1],
                    emp["analytic_key"],
                )
            )
    lines.append("")
    lines.append(
        "E2 scaled_residual (n=%d): mean=%.4g sd=%.4g min=%.4g max=%.4g "
        "ci95_mean=[%.4g, %.4g]"
        % (
            stats_out["e2"]["scaled_residual"]["n"],
            stats_out["e2"]["scaled_residual"]["mean"],
            stats_out["e2"]["scaled_residual"]["sd"],
            stats_out["e2"]["scaled_residual"]["min"],
            stats_out["e2"]["scaled_residual"]["max"],
            stats_out["e2"]["scaled_residual"]["ci_mean"][0],
            stats_out["e2"]["scaled_residual"]["ci_mean"][1],
        )
    )
    lines.append(
        "E2 loss_residual_fro (n=%d): min=%.4g mean=%.4g max=%.4g"
        % (
            stats_out["e2"]["loss_residual_fro"]["n"],
            stats_out["e2"]["loss_residual_fro"]["min"],
            stats_out["e2"]["loss_residual_fro"]["mean"],
            stats_out["e2"]["loss_residual_fro"]["max"],
        )
    )
    lines.append("")
    for output_key, entry in e3_out.items():
        lines.append(
            "E3 %s (actual key %s, n=%d): mean=%.4g sd=%.4g min=%.4g max=%.4g "
            "ci95_mean=[%.4g, %.4g]"
            % (
                output_key,
                entry["source_key"],
                entry["n"],
                entry["mean"],
                entry["sd"],
                entry["min"],
                entry["max"],
                entry["ci_mean"][0],
                entry["ci_mean"][1],
            )
        )
    lines.append("")
    id_entry = stats_out["e5"]["identity_fro_residual"]
    lines.append(
        "E5 identity_fro_residual (n=%d): mean=%.4g sd=%.4g min=%.4g max=%.4g "
        "ci95_mean=[%.4g, %.4g] count<=1e-12=%d"
        % (
            id_entry["n"],
            id_entry["mean"],
            id_entry["sd"],
            id_entry["min"],
            id_entry["max"],
            id_entry["ci_mean"][0],
            id_entry["ci_mean"][1],
            id_entry["count_le_1e-12"],
        )
    )
    flag_entry = stats_out["e5"]["loewner_flag"]
    lines.append(
        "E5 loewner_I_ge_L flag (I_acq>=L_rank; n=%d): true=%d false=%d"
        % (flag_entry["n"], flag_entry["true"], flag_entry["false"])
    )
    lines.append("")
    p1 = stats_out["physical_optional"]["phase1"]
    ps1 = p1["per_seed_max_t5a_backward_scaled_residual"]
    lines.append(
        "Physical phase1 per-seed max T5a backward_scaled_residual (n=%d): "
        "min=%.4g mean=%.4g max=%.4g"
        % (ps1["n"], ps1["min"], ps1["mean"], ps1["max"])
    )
    lines.append(
        "Physical phase1 aggregate: max_t5a_backward_scaled_residual=%.4g "
        "max_theorem3_residual_max_abs=%.4g"
        % (
            p1["aggregate"]["max_t5a_backward_scaled_residual"],
            p1["aggregate"]["max_theorem3_residual_max_abs"],
        )
    )
    p2 = stats_out["physical_optional"]["phase2"]
    ps2 = p2["per_seed_min_eig_J_stack_minus_sum"]
    lines.append(
        "Physical phase2 per-seed min_eig_J_stack_minus_sum (n=%d): "
        "min=%.4g mean=%.4g max=%.4g"
        % (ps2["n"], ps2["min"], ps2["mean"], ps2["max"])
    )
    lines.append(
        "Physical phase2 aggregate: max_innovation_rel_residual=%.4g "
        "max_budget_identity_rel_residual=%.4g"
        % (
            p2["aggregate"]["max_innovation_rel_residual"],
            p2["aggregate"]["max_budget_identity_rel_residual"],
        )
    )

    txt_path = RESULTS / "ci_stats.txt"
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # -------------------------------------------------- short confirmation
    if missing:
        print("Keys not found:", "; ".join(missing))
    else:
        print("Confirmation: all requested keys were found in the result JSONs.")
    print("Wrote", json_path)
    print("Wrote", txt_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
