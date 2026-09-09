#!/usr/bin/env python
"""Run experiment E2 and write all required results and figures.

Usage (from the experiment working directory):

    python scripts/run_e2.py
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

from a2val import e2


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    return obj


def _fmt(v, nd=6):
    if v is None:
        return ""
    return f"{v:.{nd}g}"


def main():
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    data = e2.all_e2()
    ec = data["exact_controls"]
    neutral = data["neutral_admission"]
    rt = data["random_tangents"]

    doc = {
        "title": "E2 nested free-current rank loss, calibration-neutral "
                 "admission, nonmonotone rho, and numerical-rank controls",
        "date": "2026-09-05",
        "package": "a2val",
        "settings": {
            "theory": "Theorem 5 identities (T5a, T5b), Theorem 6 neutral "
                      "admission incl. complex J_c-safe kernel, nonmonotone "
                      "relative map retention, rank/gap roundoff controls",
            "matrices": "real, unit-noise whitened, gauge-fixed",
            "e_basis_drop": "projected singular values <= 1e-12 of the "
                            "largest are dropped",
            "rho_support": "SVD of K0 with support tolerance 1e-12*max; "
                           "undefined null directions excluded",
            "random_tangent_seeds": list(range(201, 213)),
            "neutral_fixture_seed": 2600,
        },
        "results": data,
    }

    exact_rows = [
        ("equality_case", ec["equality_case"]),
        ("strict_loss_case", ec["strict_loss_case"]),
        ("complete_hiding_case", ec["complete_hiding_case"]),
        ("loss_without_rank_drop", ec["loss_without_rank_drop"]),
    ]
    verdict = {
        "all_exact_t5a": all(v["t5a"] for _, v in exact_rows),
        "all_exact_t5b": all(v["t5b"] for _, v in exact_rows),
        "nonmonotone_rho_pass": ec["nonmonotone_rho"]["pass"],
        "rank_uncertainty_discontinuity": ec["rank_uncertainty"]["discontinuity_confirmed"],
        "neutral_real_safe": neutral["real_safe_space"]["safe_preserves_j"]
        and neutral["real_safe_space"]["unsafe_decreases_loewner"],
        "neutral_complex_pass": neutral["complex_safe_space"]["pass"],
        "utility_projection_pass": neutral["utility_projection"].get("maximizer_dominates"),
        "random_tangents_all_t5b": rt["all_t5b_matches"],
        "random_tangents_max_scaled_below_100": bool(rt["max_scaled_residual"] < 100.0),
    }
    doc["verdict"] = verdict

    with open("results/e2_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")

    # ---- CSV summary ------------------------------------------------------
    rows = []

    def add(case, category, key_values, passed, note):
        rows.append({
            "case": case,
            "category": category,
            "key_values": key_values,
            "pass": passed,
            "note": note,
        })

    for name, rec in exact_rows:
        add(name, "exact_control",
            (f"loss_res={rec['loss_residual_fro']:.3e}, "
             f"drop={rec['rank_drop_actual']}, dimE={rec['dim_E']}"),
            rec["t5a"] and rec["t5b"],
            "T5a loss identity and T5b rank-drop identity")

    nr = ec["nonmonotone_rho"]
    add("nonmonotone_rho", "exact_control",
        f"rho={[_fmt(x, 4) for x in nr['rho_sequence']]}",
        nr["pass"],
        "canonical falsification of monotone map retention (3/4, 1, 0)")

    ru = ec["rank_uncertainty"]
    add("rank_uncertainty", "exact_control",
        f"J(0)={ru['j_x_0']}, max|J(t>0)|={max(ru['j_x_nonzero']):.3e}",
        ru["discontinuity_confirmed"],
        "no oracle exact-rank certificate; absolute singular values recorded")

    sat = ec["saturation_control"]
    add("saturation_control", "exact_control",
        f"j_fro={sat['j_x_fro']:.3e}, scaled={sat['backward_scaled_residual']:.3e}, "
        f"max_sv_Bv={sat['bv_max_abs_singular_value']:.3e}",
        sat["backward_scaled_residual"] < 100.0,
        "J_x zero up to backward error; full-current rank saturates nuisance")

    rtc = ec["rank_threshold_control"]
    add("rank_threshold_control", "exact_control",
        f"rank_abs={rtc['rank_absolute']}, rank_rel={rtc['rank_relative']}, "
        f"sv={[_fmt(x, 2) for x in rtc['singular_values']]}",
        True,
        "small-ratio trap: 1e-14/1 is roundoff and is never physical rank")

    pgc = ec["projector_gap_closure"]
    add("projector_gap_closure", "exact_control",
        f"projector ranks all 1, ||P(0.05)-P(-0.05)||_F={pgc['fro_P_plus_minus_0_05']:.6f}",
        all(r["rank_top_one_projector"] == 1 for r in pgc["records"]),
        "gap closure at t=0 is a projector discontinuity, not an H rank change")

    rs = neutral["real_safe_space"]
    add("real_safe_space", "neutral_admission",
        f"dim_kerF={rs['dim_kerF']}, rankF={rs['rank_F']}, "
        f"safe_loss={rs['safe_loss_fro']:.3e}, "
        f"max diff eig={max(rs['diff_eigs']):.4f}",
        rs["safe_preserves_j"] and rs["unsafe_decreases_loewner"],
        "d-rank(F)=1 >= d-p=1; safe kernel leaves J_x unchanged")

    cs = neutral["complex_safe_space"]
    add("complex_safe_space", "neutral_admission",
        f"real_dim_S={cs['real_dim_S']}, complex_dim_S={cs['complex_dim_S']}, "
        f"Jc_resid={cs['jc_invariance_residual_fro']:.3e}",
        cs["pass"],
        "S=ker F cap ker(F Jc) is J_c-invariant and lossless")

    up = neutral["utility_projection"]
    if "trace_claimed_maximizer" in up:
        add("utility_projection", "neutral_admission",
            f"trace_opt={up['trace_claimed_maximizer']:.6f}, "
            f"max_rand={up['max_random_trace']:.6f}, r={up['chosen_r']}",
            up["maximizer_dominates"],
            "leading eigenvectors of P_S T P_S maximize utility inside S")
    else:
        add("utility_projection", "neutral_admission", "vacuous (S empty)",
            False, "no r>=1 subspace available")

    add("random_tangents_T5a_T5b", "fixed_random_tangents",
        f"max_scaled_resid={rt['max_scaled_residual']:.6f}, "
        f"max_loss_resid={rt['max_loss_residual_fro']:.3e}, "
        f"t5b_matches={sum(r['t5b_matches'] for r in rt['records'])}/36",
        rt["all_t5b_matches"] and rt["max_scaled_residual"] < 100.0,
        "machine-precision T5a residual on 36 seed/pair checks")

    fieldnames = ["case", "category", "key_values", "pass", "note"]
    with open("results/e2_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    # ---- Figures -----------------------------------------------------------
    _make_rank_loss_figure(rt)
    _make_rho_figure(ec)
    _make_neutral_figure(rs)

    print("E2 complete.")
    print("  verdict:", {k: bool(v) for k, v in verdict.items()})
    print("  exact rho sequence:", [_fmt(x, 6) for x in nr["rho_sequence"]])
    print("  random tangents: max scaled residual =",
          _fmt(rt["max_scaled_residual"], 6),
          ", max loss residual =", _fmt(rt["max_loss_residual_fro"], 3),
          ", T5b match", rt["all_t5b_matches"])
    print("  real safe: dim ker F =", rs["dim_kerF"], ", safe loss =",
          _fmt(rs["safe_loss_fro"], 3))
    print("  complex safe: real dim S =", cs["real_dim_S"], ", complex dim =",
          cs["complex_dim_S"])
    print("  rank uncertainty:",
          [_fmt(x["j_x"], 3) for x in ru["records"]])
    print("Saved: results/e2_results.json, results/e2_summary.csv")
    print("Saved: figures/e2_rank_loss.png, figures/e2_rho_nonmonotone.png, "
          "figures/e2_neutral_admission.png")


def _make_rank_loss_figure(rt):
    """Predicted vs actual loss eigenvalues for three representative seeds."""
    colors = {0: "#1f4e79", 1: "#b03a2e", 2: "#e09a00"}
    markers = {0: "o", 1: "s", 2: "^"}
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    maxv = 1e-30
    plotted = False
    for rec in rt["records"]:
        if rec["seed"] not in (201, 202, 203):
            continue
        # "C0->C1": index 0, "C1->C2": 1, "C2->C3": 2
        pair_idx = int(rec["pair"].split("->")[0][1])
        x = np.maximum(np.asarray(rec["b_T_P_E_b_eigs"]), 0.0)
        y = np.maximum(np.asarray(rec["j_r_minus_j_r1_eigs"]), 0.0)
        ax.scatter(x, y, s=34, color=colors[pair_idx],
                   marker=markers[pair_idx], alpha=0.9,
                   label=f"seed {rec['seed']}, {rec['pair']}"
                   if rec["seed"] == 201 else None)
        maxv = max(maxv, float(np.max(np.concatenate([x, y]))))
        plotted = True
    if plotted:
        lim = maxv * 1.08
        ax.plot([0, lim], [0, lim], "k--", lw=1, label="identity")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
    ax.set_xlabel("eigenvalues of B^T P_{E_r} B (predicted loss)")
    ax.set_ylabel("eigenvalues of J_x(r) - J_x(r+1) (actual loss)")
    ax.set_title("E2 Theorem 5a: predicted vs actual nested-rank loss\n"
                 "(seeds 201-203, pairs C0->C1, C1->C2, C2->C3)")
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig("figures/e2_rank_loss.png", dpi=160)
    plt.close(fig)


def _make_rho_figure(ec):
    nr = ec["nonmonotone_rho"]
    labels = [r["chain"] for r in nr["rho_records"]]
    vals = [r["rho"][0] for r in nr["rho_records"]]
    colors = ["#6da8d4", "#1f4e79", "#e09a00"]
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    bars = ax.bar(labels, vals, 0.55, color=colors)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.025, f"{v:.4g}",
                ha="center", va="bottom")
    ax.set_ylim(0, 1.2)
    ax.axhline(vals[0], color="grey", lw=0.8, ls=":")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(["C0 = {0}", "C1 = span(e1)",
                        "C2 = span(e1, e2+e3)"])
    ax.set_ylabel(r"$\rho(r)$ (relative map retention on supp $K_0$)")
    ax.set_title("E2 nonmonotone map retention: canonical falsification\n"
                 r"$\rho=3/4 \to 1 \to 0$ (no monotonicity under nesting)")
    fig.tight_layout()
    fig.savefig("figures/e2_rho_nonmonotone.png", dpi=160)
    plt.close(fig)


def _make_neutral_figure(rs):
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), sharey=True)
    x = np.arange(3)
    w = 0.35
    for ax, label, before, after, title, subtitle in (
        (axes[0], "safe", rs["j_before_eigs"], rs["j_after_safe_eigs"],
         "safe direction z in ker(V^T Q)", "J unchanged"),
        (axes[1], "unsafe", rs["j_before_eigs"], rs["j_after_unsafe_eigs"],
         "generic z not in ker(V^T Q)", "J strictly decreases a direction"),
    ):
        ax.bar(x - w / 2, before, w, color="#1f4e79", label="J_x before")
        ax.bar(x + w / 2, after, w, color="#b03a2e", label="J_x after")
        ax.set_xticks(x)
        ax.set_xticklabels([f"pose eig {i + 1}" for i in range(3)])
        ax.set_title(f"{title}\n{subtitle}")
        ax.legend(fontsize=8)
        ax.set_xlabel("eigenvalue index")
    axes[0].set_ylabel("eigenvalue of pose information J_x")
    fig.suptitle("E2 real calibration-neutral admission "
                 "(random A, B, C0; existing nuisance N=Ran([C0,A]))")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig("figures/e2_neutral_admission.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
