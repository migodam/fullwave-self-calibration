#!/usr/bin/env python
"""Run experiment E1 and write all required results/figures.

Usage (from the experiment working directory):

    python scripts/run_e1.py
"""

from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import e1

SEEDS = list(range(101, 111))
N_SAMPLES = 2000
N_BINS = 40

THETA = {
    "phase_only": 1.0,
    "nuisance_phase": 1.0,
    "scale_family": 0.0,
    "total_field_ref": np.pi / 2.0,
    "mismatch_intensity_noise": np.pi / 2.0,
}


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    return obj


def _mean_over_seeds(records, key):
    vals = [r[key] for r in records if r.get(key) is not None]
    return float(np.mean(vals)) if vals else None


def _max_over_seeds(records, key):
    vals = [r[key] for r in records if r.get(key) is not None]
    return float(np.max(vals)) if vals else None


def main():
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    analytic = {}
    per_seed = {}
    for model in e1.CONTROL_MODELS:
        theta = THETA[model]
        analytic[model] = e1.analytic_values(model, theta)
        # Attach full quadrature diagnostics where the phaseless law is used.
        if model == "total_field_ref":
            analytic[model]["j_ph_detail"] = e1.phaseless_fisher_detail(model, theta)
        recs = [
            e1.mc_score_check(model, theta, seed=seed, n_samples=N_SAMPLES,
                              n_bins=N_BINS)
            for seed in SEEDS
        ]
        per_seed[model] = recs

    # Mismatch control (analytic Gaussian-sensor Fisher only; the exact
    # intensity likelihood is not defined for this deliberately mismatched
    # sensor, so no matched-score identity / no exact J_ph).
    mismatch_model = "mismatch_intensity_noise"
    analytic[mismatch_model] = e1.analytic_values(mismatch_model, THETA[mismatch_model])
    per_seed[mismatch_model] = []

    # ---- Assemble the structured JSON -------------------------------------
    doc = {
        "title": "E1 matched coherent vs phaseless information contraction "
                 "(exact induced intensity likelihoods)",
        "date": "2026-09-05",
        "package": "a2val",
        "settings": {
            "seeds": SEEDS,
            "n_samples_per_seed": N_SAMPLES,
            "n_bins": N_BINS,
            "contraction_models": list(e1.CONTROL_MODELS),
            "mismatch_control": mismatch_model,
            "quantile_bounds": "ncx2 ppf(1e-12) .. ppf(1-1e-12) with full-domain cross check",
        },
        "models": {},
    }
    for model in list(e1.CONTROL_MODELS) + [mismatch_model]:
        recs = per_seed[model]
        m = dict(analytic[model])
        m["per_seed"] = _jsonable(recs)
        m["j_ph_emp_mean_over_seeds"] = _mean_over_seeds(recs, "j_ph_emp")
        m["j_ph_emp_se_mean_over_seeds"] = _mean_over_seeds(recs, "j_ph_emp_se")
        m["j_coh_emp_mean_over_seeds"] = _mean_over_seeds(recs, "j_coh_emp")
        m["j_coh_emp_se_mean_over_seeds"] = _mean_over_seeds(recs, "j_coh_emp_se")
        m["max_bin_weighted_rmse_over_seeds"] = _max_over_seeds(recs, "bin_weighted_rmse")
        m["max_bin_mean_diff_over_seeds"] = _max_over_seeds(recs, "max_bin_mean_diff")
        m["max_std_bin_diff_over_seeds"] = _max_over_seeds(recs, "max_std_bin_diff")
        m["identity_within_4se_all_seeds"] = all(
            r["identity_within_4_per_bin_se"] for r in recs
        ) if recs else None
        m["score_expectation_ok_all_seeds"] = all(
            r["score_expectation_ok"] for r in recs
        ) if recs else None
        m["mean_mean_coh_over_seeds"] = _mean_over_seeds(recs, "mean_coh")
        m["mean_mean_ph_over_seeds"] = _mean_over_seeds(recs, "mean_ph")
        if model == "nuisance_phase":
            m["j_eff_coh_emp_mean_over_seeds"] = _mean_over_seeds(recs, "j_eff_coh_emp")
            m["j_coh_raw_emp_mean_over_seeds"] = _mean_over_seeds(recs, "j_coh_raw_emp")
        doc["models"][model] = _jsonable(m)

    # ---- Pass/fail summary ------------------------------------------------
    summary = []
    for model in e1.CONTROL_MODELS:
        a = analytic[model]
        recs = per_seed[model]
        emp_coh = doc["models"][model]["j_coh_emp_mean_over_seeds"]
        emp_ph = doc["models"][model]["j_ph_emp_mean_over_seeds"]
        max_rmse = doc["models"][model]["max_bin_weighted_rmse_over_seeds"]
        if model == "phase_only":
            passed = (a["j_ph_analytic"] < a["j_coh_analytic"]) and all(
                r["mc_contraction_ok"] for r in recs
            )
            note = "strict contraction: J_ph=0 < J_coh=2"
        elif model == "nuisance_phase":
            passed = (a["j_eff_coh_analytic"] == a["j_ph_analytic"] == 0.0)
            note = "efficient coherent and phaseless target information both zero"
        elif model == "scale_family":
            passed = (abs(a["j_ph_analytic"] - a["j_coh_analytic"]) < 1e-6) and all(
                r["mc_contraction_ok"] for r in recs
            )
            note = "positive equality control: intensity sufficient for scale score"
        else:  # total_field_ref
            passed = (0.0 < a["j_ph_analytic"] < a["j_coh_analytic"]) and all(
                r["mc_contraction_ok"] for r in recs
            )
            note = "strict contraction with known reference (0 < J_ph < J_coh)"
        summary.append({
            "model": model,
            "theta": a["theta"],
            "j_coh_analytic": a["j_coh_analytic"],
            "j_ph_analytic": a["j_ph_analytic"],
            "j_coh_emp_mean": emp_coh,
            "j_ph_emp_mean": emp_ph,
            "max_bin_weighted_rmse": max_rmse,
            "contraction_pass": bool(passed),
            "note": note,
        })
    mm = analytic[mismatch_model]
    summary.append({
        "model": mismatch_model,
        "theta": mm["theta"],
        "j_coh_analytic": mm["j_coh_analytic"],
        "j_ph_analytic": None,
        "j_mismatch_analytic": mm["j_mismatch_analytic"],
        "j_coh_emp_mean": None,
        "j_ph_emp_mean": None,
        "max_bin_weighted_rmse": None,
        "contraction_pass": False,
        "note": "deliberate mismatch control, NOT covered by the contraction theorem",
    })
    doc["summary"] = summary
    doc["summary_verdict"] = {
        "all_control_contraction_pass": all(s["contraction_pass"] for s in summary if s["model"] != mismatch_model),
        "all_control_score_expectation_ok": all(
            doc["models"][model]["score_expectation_ok_all_seeds"]
            for model in e1.CONTROL_MODELS
        ),
        "all_control_identity_within_4se": all(
            doc["models"][model]["identity_within_4se_all_seeds"]
            for model in e1.CONTROL_MODELS
        ),
        "mismatch_control_violates_ordering_by_design": True,
    }
    with open("results/e1_results.json", "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")

    # ---- CSV summary --------------------------------------------------------
    fieldnames = [
        "model", "j_coh_analytic", "j_ph_analytic", "j_coh_emp_mean",
        "j_ph_emp_mean", "max_bin_weighted_rmse", "contraction_pass",
        "j_mismatch_analytic", "note",
    ]
    with open("results/e1_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary:
            out = {k: row.get(k) for k in fieldnames}
            if out["j_ph_analytic"] is None:
                out["j_ph_analytic"] = ""
            if out["j_coh_emp_mean"] is None:
                out["j_coh_emp_mean"] = ""
            if out["j_ph_emp_mean"] is None:
                out["j_ph_emp_mean"] = ""
            if out["max_bin_weighted_rmse"] is None:
                out["max_bin_weighted_rmse"] = ""
            if "j_mismatch_analytic" not in out or out["j_mismatch_analytic"] is None:
                out["j_mismatch_analytic"] = ""
            writer.writerow(out)

    # ---- Figures -------------------------------------------------------------
    rep_seed = SEEDS[0]
    model4 = "total_field_ref"
    theta4 = THETA[model4]

    # --- e1_fisher_comparison.png ------------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(13.5, 5.0), gridspec_kw={"width_ratios": [3.7, 1.3]}
    )
    colors = {
        "coh_a": "#1f4e79",
        "ph_a": "#b03a2e",
        "coh_e": "#6da8d4",
        "ph_e": "#e5a199",
    }
    models_left = ["phase_only", "nuisance_phase", "scale_family", "total_field_ref"]
    x = np.arange(len(models_left))
    w = 0.19
    for i, model in enumerate(models_left):
        a = analytic[model]
        emp_coh = doc["models"][model]["j_coh_emp_mean_over_seeds"]
        emp_ph = doc["models"][model]["j_ph_emp_mean_over_seeds"]
        ax1.bar(i - 1.5 * w, a["j_coh_analytic"], w, color=colors["coh_a"], label="J_coh analytic" if i == 0 else None)
        ax1.bar(i - 0.5 * w, a["j_ph_analytic"], w, color=colors["ph_a"], label="J_ph analytic" if i == 0 else None)
        if emp_coh is not None:
            ax1.bar(i + 0.5 * w, emp_coh, w, color=colors["coh_e"], label="J_coh empirical (mean seeds)" if i == 0 else None)
        if emp_ph is not None:
            ax1.bar(i + 1.5 * w, emp_ph, w, color=colors["ph_e"], label="J_ph empirical (mean seeds)" if i == 0 else None)
    ax1.set_xticks(x)
    ax1.set_xticklabels(["phase_only", "nuisance_phase\n(efficient coh.)", "scale_family", "total_field_ref"])
    ax1.set_ylabel("Fisher information")
    ax1.set_title("Contraction models 1-4 (J_ph <= J_coh)")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.set_ylim(0, 2.6)

    mm = analytic[mismatch_model]
    labels = ["J_coh\nanalytic", "J_mismatch\n(direct Gaussian sensor)"]
    vals = [mm["j_coh_analytic"], mm["j_mismatch_analytic"]]
    bars = ax2.bar([0, 1], vals, 0.5, color=[colors["coh_a"], "#e09a00"])
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{v:.3g}",
                 ha="center", va="bottom", fontsize=9)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_ylabel("Fisher information")
    ax2.set_title("Model 5: mismatch control\n(ordering violation by design)")
    ax2.set_ylim(0, 2250)
    fig.suptitle("E1: coherent vs exact-intensity Fisher; mismatch sensor is not an exact intensity likelihood")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig("figures/e1_fisher_comparison.png", dpi=160)
    plt.close(fig)

    # Rerun one seed of model 4 with detailed bin arrays for the plots.
    rng = np.random.default_rng(rep_seed)
    sigma2 = e1.model_sigma2(model4, theta4)
    mu = e1.model_mu(model4, theta4)
    dmu = e1.model_dmu(model4, theta4)
    y = e1.draw_complex_gaussian(rng, complex(mu), sigma2, N_SAMPLES)
    I = np.abs(y) ** 2
    s_coh = np.asarray(e1.coherent_score(y, mu, dmu, sigma2, u=0.0))
    s_ph = np.asarray(e1.phaseless_score(model4, I, theta4))
    bin_res = e1._binned_check(I, s_coh, s_ph, n_bins=N_BINS)

    # --- e1_conditional_expectation.png ------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.plot([-1.2, 1.2], [-1.2, 1.2], color="grey", lw=1, ls="--", label="identity y=x")
    sc = ax.scatter(bin_res["bin_mean_s_ph"], bin_res["bin_mean_s_coh"],
                    c=bin_res["bin_centers"], cmap="viridis", s=42, zorder=3)
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("bin mean s_ph(I)")
    ax.set_ylabel("bin mean s_coh(Y)")
    ax.set_title(f"E1 model 4 (total_field_ref), seed {rep_seed}: "
                 f"conditional-expectation identity\nslope={bin_res['regression_slope']:.3f}, "
                 f"R2={bin_res['regression_r2']:.4f}, bin RMSE={bin_res['bin_weighted_rmse']:.4f}")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("bin center I")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig("figures/e1_conditional_expectation.png", dpi=160)
    plt.close(fig)

    # --- e1_score_scatter.png -----------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(I, s_coh, s=2, alpha=0.16, color="#1f4e79", label="per-sample s_coh(Y)")
    ax.plot(bin_res["bin_centers"], bin_res["bin_mean_s_coh"], "-o", color="#b03a2e",
            lw=1.8, ms=4, label="bin mean s_coh(Y)")
    ax.plot(bin_res["bin_centers"], bin_res["bin_mean_s_ph"], "s", color="#e09a00",
            ms=4, alpha=0.85, label="bin mean s_ph(I)")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("intensity I = |Y|^2")
    ax.set_ylabel("score")
    ax.set_title(f"E1 model 4 (total_field_ref), seed {rep_seed}: coherent score vs I\n"
                 f"(bin means of s_coh and s_ph overlap per the matched identity)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig("figures/e1_score_scatter.png", dpi=160)
    plt.close(fig)

    # ---- Print compact report ------------------------------------------------
    print("E1 complete. Analytic + empirical results:")
    for row in doc["summary"]:
        print(row)
    print("\nSaved: results/e1_results.json, results/e1_summary.csv")
    print("Saved: figures/e1_fisher_comparison.png, "
          "figures/e1_conditional_expectation.png, figures/e1_score_scatter.png")


if __name__ == "__main__":
    main()
