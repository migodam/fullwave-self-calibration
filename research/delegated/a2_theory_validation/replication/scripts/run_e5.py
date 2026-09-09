#!/usr/bin/env python
"""Run experiment E5 and write all required results and figures.

Usage (from the experiment working directory):

    python scripts/run_e5.py
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

from a2val import e5


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

    data = e5.all_e5()
    sa = data["shared_map_algebra"]
    rb = data["rank_acquisition_budget"]
    pol = data["policies"]
    gauge = data["gauge_symmetry"]

    min_pair = sa["minimal_pair"]
    rand_inn = sa["random_innovation"]
    rand_sweep = sa["random_innovation_sweep"]
    singular = sa["singular_fallback"]
    comp = sa["compensation_criterion"]
    svs = sa["stack_vs_sum"]
    sat = sa["saturated_new_frame"]
    budget = rb["seeds"]
    crafted = rb["crafted_cases"]
    neutral = rb["neutral_control"]
    nonsub = pol["nonsubmodularity"]
    policies = pol["acquisition_policies"]
    neutral_ld = pol["neutral_logdet"]
    gauge_null = gauge["global_gauge_null"]
    anchor = gauge["anchor_restores_absolute"]
    source = gauge["source_stabilizer_note"]

    verdict = {
        "minimal_pair": bool(min_pair["pass"]),
        "random_innovation": bool(rand_inn["pass"]),
        "random_innovation_sweep": bool(rand_sweep["pass"]),
        "singular_fallback": bool(singular["pass"]),
        "compensation_criterion": bool(comp["pass"]),
        "stack_vs_sum": bool(svs["pass"]),
        "saturated_new_frame": bool(sat["pass"]),
        "rank_budget_seeds": bool(budget["pass"]),
        "rank_budget_crafted_compensate": bool(crafted["compensate_ok"]),
        "rank_budget_crafted_noncompensate": bool(crafted["noncompensate_ok"]),
        "rank_budget_neutral_control": bool(neutral["pass"]),
        "nonsubmodularity": bool(nonsub["pass"]),
        "acquisition_policies_run": bool(policies["pass"]),
        "acquisition_policies_greedy_optimal": bool(
            policies["greedy_is_optimal"]),
        "gauge_null_remains_null": bool(gauge_null["pass"]),
        "anchor_restores_absolute": bool(anchor["pass"]),
        "source_stabilizer": source["status"],
    }

    doc = {
        "title": ("E5 shared-map acquisition innovation, rank-acquisition "
                  "information budget, non-submodularity and policies, and "
                  "gauge/symmetry controls"),
        "date": "2026-09-05",
        "package": "a2val",
        "settings": {
            "theory": ("Theorem 9 shared-map innovation and compensation, "
                       "Corollary 10 rank-acquisition budget, pose-Schur "
                       "logdet non-submodularity, greedy/pair/exhaustive/"
                       "random policies"),
            "matrices": "real, unit-noise whitened, gauge-fixed",
            "matrices_note": ("Each frame is current-cleaned before use "
                              "(a_l,b_l); only full-column-rank map models "
                              "use the Theorem 9 inverse formula, otherwise "
                              "the singular fallback/variational stack is "
                              "used."),
            "random_innovation_seed": rand_inn["seed"],
            "random_innovation_sweep_seeds": list(range(2501, 2521)),
            "singular_fallback_seed": singular["seed"],
            "compensation_seed": comp["seed"],
            "stack_vs_sum_seed": svs["seed"],
            "saturated_seed": sat["seed"],
            "rank_budget_seeds": list(range(401, 413)),
            "rank_budget_crafted_seed": crafted["seed"],
            "neutral_control_seed": neutral["seed"],
            "policy_seed": policies["seed"],
            "policy_objective": policies["objective_note"],
            "policy_random_seeds": list(range(401, 413)),
            "loewner_tolerance": 1e-10,
            "identity_tolerance_fro": 1e-12,
        },
        "results": data,
        "verdict": verdict,
    }

    with open("results/e5_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")

    rows = []

    def add(case, category, key_values, passed, note):
        rows.append({
            "case": case,
            "category": category,
            "key_values": key_values,
            "pass": passed,
            "note": note,
        })

    add("minimal_pair_contrast", "shared_map_algebra",
        f"J1={_fmt(float(np.linalg.norm(min_pair['J1'])), 2)}, "
        f"J2={_fmt(float(np.linalg.norm(min_pair['J2_contrast'])), 2)}, "
        f"J_stack={_fmt(float(min_pair['J_stack_contrast'].flat[0]), 6)}, "
        f"sum_J={_fmt(float(min_pair['sum_J_contrast'].flat[0]), 2)}, "
        f"J_stack>sum={min_pair['J_stack_gt_sum']}",
        min_pair["pass"],
        "per-frame pose info zero; shared-map stack has information two")
    add("minimal_pair_duplicate", "shared_map_algebra",
        f"J_stack_dup={_fmt(min_pair['J_stack_duplicate_fro'], 2)}",
        min_pair["pass"],
        "duplicated b2=b1 gives zero stacked information")
    add("random_innovation", "shared_map_algebra",
        f"fro_res={_fmt(rand_inn['formula_direct_fro_residual'], 2)}, "
        f"max_res={_fmt(rand_inn['formula_direct_max_residual'], 2)}, "
        f"V_eq_fro={_fmt(rand_inn['equality_V_fro'], 2)}, "
        f"ker_dist={_fmt(rand_inn['ker_subspace_distance'], 2)}",
        rand_inn["pass"],
        "Theorem 9 equals direct stack; V=0 no-innovation; ker identity")
    add("random_innovation_sweep", "shared_map_algebra",
        f"worst_fro={_fmt(rand_sweep['worst_fro_residual'], 2)}",
        rand_sweep["pass"],
        "20 random fixtures; formula-direct Frobenius residual")
    add("singular_fallback", "shared_map_algebra",
        f"rank_G={singular['rank_G']}, G_min_eig={_fmt(singular['G_min_eig'], 2)}, "
        f"J_direct_min_eig={_fmt(singular['J_new_psd'], 4)}, "
        f"pinv_wrong_diff={_fmt(singular['naive_pinv_wrong_fro_difference'], 4)}",
        singular["pass"],
        "singular G: naive pinv formula invalid; variational/direct PSD match")
    add("compensation_visible", "shared_map_algebra",
        f"rank_D={comp['visible']['rank_D']}, "
        f"J_stack_min_eig={_fmt(comp['visible']['min_eig_J_stack'], 4)}",
        comp["pass"],
        "distinct compensators: every frame hides its pose, stack sees it")
    add("compensation_hidden", "shared_map_algebra",
        f"rank_D={comp['hidden']['rank_D']}, "
        f"J_stack_min_eig={_fmt(comp['hidden']['min_eig_J_stack'], 2)}",
        comp["pass"],
        "equal compensators: stack stays hidden")
    add("stack_vs_sum_random", "shared_map_algebra",
        f"min_eig(J_stack-sum)="
        f"{_fmt(svs['random_three_frame']['min_eig_stack_minus_sum'], 4)}",
        svs["pass"],
        "J_stack >= sum J_l (Loewner) on a random 3-frame model")
    add("stack_vs_sum_equality", "shared_map_algebra",
        f"gap_fro={_fmt(svs['equality_nontrivial_orthogonal_residual']['fro_gap'], 2)}, "
        f"J_stack_min_eig="
        f"{_fmt(svs['equality_nontrivial_orthogonal_residual']['J_stack_min_eig'], 4)}",
        svs["pass"],
        "common per-frame optimal map direction gives exact equality "
        "(nonzero J)")
    add("saturated_new_frame", "shared_map_algebra",
        f"a_clean={_fmt(sat['a_clean_fro'], 2)}, "
        f"b_clean={_fmt(sat['b_clean_fro'], 2)}, "
        f"J_change={_fmt(sat['J_change_fro'], 2)}",
        sat["pass"],
        "frame current spanning its data rows yields zero cleaned tangent")

    add("rank_budget_overall", "rank_acquisition_budget",
        f"worst_identity={_fmt(budget['worst_identity_fro_residual'], 2)}",
        budget["pass"],
        "seeds 401-412, J_final-J_original = I_acq - L_rank exact")
    loewner_ok = sum(r["loewner_I_ge_L"] for r in budget["records"])
    add("rank_budget_loewner_counts", "rank_acquisition_budget",
        f"loewner_holds={loewner_ok}/{len(budget['records'])}",
        budget["pass"],
        "count of seeds/rank pairs satisfying I_acq >= L_rank")
    c1 = [r for r in budget["records"] if r["label"] == "C0->C1"]
    c2 = [r for r in budget["records"] if r["label"] == "C1->C2"]
    add("rank_budget_C0_C1", "rank_acquisition_budget",
        f"worst_identity="
        f"{_fmt(max(r['identity_fro_residual'] for r in c1), 2)}, "
        f"loewner={sum(r['loewner_I_ge_L'] for r in c1)}/{len(c1)}",
        all(r["pass"] for r in c1),
        "empty -> 2-column current enlargement")
    add("rank_budget_C1_C2", "rank_acquisition_budget",
        f"worst_identity="
        f"{_fmt(max(r['identity_fro_residual'] for r in c2), 2)}, "
        f"loewner={sum(r['loewner_I_ge_L'] for r in c2)}/{len(c2)}",
        all(r["pass"] for r in c2),
        "2-column -> 4-column current enlargement")
    add("rank_budget_crafted_compensate", "rank_acquisition_budget",
        f"J_final-J_orig_min_eig="
        f"{_fmt(crafted['compensate']['J_final_minus_original_min_eig'], 2)}, "
        f"I-L_min_eig="
        f"{_fmt(crafted['compensate']['I_minus_L_min_eig'], 2)}",
        crafted["compensate_ok"],
        "b = aH + W^{-1/2}U_E^T B makes I_acq=L_rank exactly")
    add("rank_budget_crafted_noncompensate", "rank_acquisition_budget",
        f"I-L_min_eig="
        f"{_fmt(crafted['noncompensate']['I_minus_L_min_eig'], 4)}, "
        f"J_final-J_orig_min_eig="
        f"{_fmt(crafted['noncompensate']['J_final_minus_original_min_eig'], 4)}",
        crafted["noncompensate_ok"],
        "half-strength acquisition: I_acq < L_rank, J_final < J_original")
    add("rank_budget_neutral_control", "rank_acquisition_budget",
        f"L_max_eig={_fmt(neutral['L_rank_max_eig'], 2)}, "
        f"I_min_eig={_fmt(neutral['I_acq_min_eig'], 2)}, "
        f"delta_min_eig="
        f"{_fmt(neutral['min_eig_J_final_minus_J_original'], 2)}",
        neutral["pass"],
        "L_rank=0 neutral enlargement: PSD innovation cannot degrade J")

    for r in nonsub["records"]:
        add(f"nonsubmodular_lambda_{r['Lambda_x']}", "non_submodularity",
            f"g1+g2={_fmt(r['g1_plus_g2'], 6)}, "
            f"g12+g0={_fmt(r['g12_plus_g_empty'], 6)}, "
            f"gap={_fmt(r['submodularity_gap'], 6)}",
            r["submodularity_violated"],
            "g(1)+g(2) < g(1,2)+g(empty) violates submodularity")
    add("acquisition_policies", "policies",
        f"greedy={_fmt(policies['greedy']['logdet'], 6)}, "
        f"pair={_fmt(policies['pair_lookahead']['logdet'], 6)}, "
        f"exhaustive={_fmt(policies['exhaustive']['best_logdet'], 6)}, "
        f"random_best={_fmt(policies['random']['best_of_12']['logdet'], 6)}, "
        f"greedy_gap={_fmt(policies['greedy_optimality_gap'], 6)}",
        True,
        "K=3 over 8 candidates on the fixed old model; greedy is "
        f"{'optimal' if policies['greedy_is_optimal'] else 'not optimal'}")
    add("policy_evaluation_counts", "policies",
        f"rank1_greedy={policies['rank1_innovation_counts']['greedy']}, "
        f"pair={policies['rank1_innovation_counts']['pair_lookahead']}, "
        f"exhaustive={policies['rank1_innovation_counts']['exhaustive']}, "
        f"random={policies['rank1_innovation_counts']['random_best_of_12']}",
        True,
        "matrix evaluations charged as rank-1 Schur updates; no forward/adjoint")
    add("neutral_logdet_lambda0", "policies",
        f"full_rank={neutral_ld['full_rank_triples']}/"
        f"{neutral_ld['total_triples']}, "
        f"best0={_fmt(neutral_ld['best_lambda0']['logdet'], 6)}, "
        f"best0.1={_fmt(neutral_ld['best_lambda_0.1']['logdet'], 6)}, "
        f"same_best={neutral_ld['same_best_subset']}",
        True,
        "Lambda_x=0 pose-Schur logdet on full-rank sets vs regularized value")

    for r in gauge_null["records"]:
        add(f"gauge_null_L{r['L']}", "gauge_symmetry",
            f"min_eig_J={_fmt(r['min_eig_J'], 2)}, "
            f"null_dim={r['tangent_null_dim']}",
            r["min_eig_J"] <= 1e-10 and r["tangent_null_dim"] == 1,
            "global rigid gauge (1,-1) survives additional relative frames")
    for r in anchor["records"]:
        add(f"anchor_restores_L{r['L']}", "gauge_symmetry",
            f"min_eig_J={_fmt(r['min_eig_J'], 4)}, "
            f"null_dim={r['tangent_null_dim']}",
            r["min_eig_J"] > 1e-10 and r["tangent_null_dim"] == 0,
            "absolute anchor row (a=0,b=1) removes the gauge")
    add("source_stabilizer", "gauge_symmetry",
        "not_yet_run", None,
        source["reason"])

    fieldnames = ["case", "category", "key_values", "pass", "note"]
    with open("results/e5_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    _make_innovation_figure(rand_sweep, singular)
    _make_budget_figure(budget)
    _make_policies_figure(policies)
    _make_nonsubmodular_figure(nonsub)

    # Final key numbers
    print("E5 all cases ran.")
    print("  verdict:", verdict)
    print(f"  minimal pair: J1={float(np.linalg.norm(min_pair['J1'])):g}, "
          f"J2={float(np.linalg.norm(min_pair['J2_contrast'])):g}, "
          f"J_stack={float(min_pair['J_stack_contrast'].flat[0]):g}, "
          f"sum J_l={float(min_pair['sum_J_contrast'].flat[0]):g}, "
          f"duplicate J_stack={min_pair['J_stack_duplicate_fro']:g}")
    print(f"  max innovation residual (single fixture): "
          f"{rand_inn['formula_direct_max_residual']:.3e}; "
          f"sweep worst Frobenius: {rand_sweep['worst_fro_residual']:.3e}")
    print(f"  singular fallback: naive-pinv wrong by "
          f"{singular['naive_pinv_wrong_fro_difference']:.6g} Frobenius")
    print(f"  budget identity worst Frobenius residual: "
          f"{budget['worst_identity_fro_residual']:.3e}")
    print(f"  non-submodularity gaps: "
          + "; ".join(f"lambda={r['Lambda_x']}: {r['submodularity_gap']:.6f}"
                      for r in nonsub["records"]))
    print(f"  policies: greedy={policies['greedy']['logdet']:.6f}, "
          f"pair={policies['pair_lookahead']['logdet']:.6f}, "
          f"exhaustive={policies['exhaustive']['best_logdet']:.6f}, "
          f"random-best-of-12={policies['random']['best_of_12']['logdet']:.6f}; "
          f"greedy-optimality gap={policies['greedy_optimality_gap']:.6g}")
    print(f"  gauge null min eig L=2: "
          f"{gauge_null['records'][0]['min_eig_J']:.3e}, "
          f"L=4: {gauge_null['records'][1]['min_eig_J']:.3e}; "
          f"anchor min eig L=2: {anchor['records'][0]['min_eig_J']:.6g}")
    print("Saved: results/e5_results.json, results/e5_summary.csv")
    print("Saved: figures/e5_innovation.png, figures/e5_budget.png, "
          "figures/e5_policies.png, figures/e5_nonsubmodular.png")


def _make_innovation_figure(sweep, singular):
    recs = sweep["records"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for r in recs:
        Jf = np.asarray(r["J_formula"])
        Jd = np.asarray(r["J_direct"])
        axes[0].scatter(Jd.ravel(), Jf.ravel(), s=13, alpha=0.65,
                        color="#1f4e79", edgecolors="none")
    lim = axes[0].get_xlim()
    lo, hi = min(lim[0], 0.0), max(lim[1], 1e-12)
    axes[0].plot([lo, hi], [lo, hi], ls="--", color="#b03a2e", lw=1.2)
    axes[0].set_xlim(lo, hi)
    axes[0].set_ylim(lo, hi)
    axes[0].set_xlabel("direct stacked J_new entries")
    axes[0].set_ylabel("Theorem 9 J_new entries")
    axes[0].set_title("J_new formula vs direct stack\n(20 random fixtures)")
    axes[0].grid(alpha=0.25)

    seeds = [r["seed"] for r in recs]
    res = [r["fro_residual"] for r in recs]
    axes[1].scatter(seeds, res, s=26, color="#1f4e79", label="random fixtures")
    axes[1].axhline(1e-12, color="#b03a2e", ls="--",
                    label="1e-12 acceptance line")
    axes[1].set_yscale("log")
    axes[1].set_ylim(min(min(res), 1e-15) * 0.5, max(max(res), 1e-10) * 5)
    axes[1].set_xlabel("seed")
    axes[1].set_ylabel("Frobenius residual")
    axes[1].set_title("Innovation residual (log scale)")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.25, which="both")
    fig.suptitle("E5 Theorem 9: exact shared-map innovation identity", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("figures/e5_innovation.png", dpi=160)
    plt.close(fig)


def _make_budget_figure(budget):
    reps = [r for r in budget["records"] if r["seed"] in range(401, 405)]
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    colors = ["#1f4e79", "#2e75b6", "#e09a00", "#b03a2e"]
    for r in reps:
        xs = np.asarray(r["I_minus_L_eigs"])
        ys = np.asarray(r["delta_eigs"])
        ax.scatter(xs, ys, s=30, color=colors[(r["seed"] - 401) % 4],
                   label=f"seed {r['seed']} ({r['label']})", alpha=0.9)
    lo = min(min(r["I_minus_L_eigs"].min() for r in reps),
             min(r["delta_eigs"].min() for r in reps))
    hi = max(max(r["I_minus_L_eigs"].max() for r in reps),
             max(r["delta_eigs"].max() for r in reps))
    ax.plot([lo - 0.2, hi + 0.2], [lo - 0.2, hi + 0.2], ls="--",
            color="#7f7f7f", lw=1.0, label="exact equality")
    ax.axhline(0, color="#c8c8c8", lw=0.8)
    ax.axvline(0, color="#c8c8c8", lw=0.8)
    ax.text(0.97, 0.97, "$I_{\\mathrm{acq}} > L_{\\mathrm{rank}}$\n"
            "(pose information improves)",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            color="#1f4e79")
    ax.text(0.03, 0.03, "$I_{\\mathrm{acq}} < L_{\\mathrm{rank}}$\n"
            "(pose information degrades)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=9,
            color="#b03a2e")
    ax.set_xlabel("eigenvalues of $I_{\\mathrm{acq}}-L_{\\mathrm{rank}}$")
    ax.set_ylabel("eigenvalues of $J_{\\mathrm{final}}-J_{\\mathrm{original}}$")
    ax.set_title("E5 Corollary 10: rank-acquisition budget identity\n"
                 "(4 representative seeds, both rank steps; Loewner sign "
                 "matches exactly)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig("figures/e5_budget.png", dpi=160)
    plt.close(fig)


def _make_policies_figure(policies):
    det = [policies["greedy"]["logdet"],
           policies["pair_lookahead"]["logdet"],
           policies["exhaustive"]["best_logdet"]]
    random_vals = [r["logdet"] for r in policies["random"]["rows"]]
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    x = np.arange(4)
    ax.scatter(x[:3], det, s=70, color="#1f4e79", zorder=3,
               label="policy final logdet")
    bp = ax.boxplot(random_vals, positions=[3], widths=0.45,
                    showfliers=False, patch_artist=True)
    bp["boxes"][0].set_facecolor("#dce6f1")
    bp["medians"][0].set_color("#1f4e79")
    rng = np.random.default_rng(5)
    jitter = rng.uniform(-0.12, 0.12, size=len(random_vals))
    ax.scatter(3 + jitter, random_vals, s=24, color="#b03a2e", zorder=3,
               alpha=0.8, label="random selection (seeds 401-412)")
    ax.axhline(policies["exhaustive"]["best_logdet"], color="#2e75b6",
               ls=":", lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(["greedy", "pair look-ahead", "exhaustive",
                        "random\n(12 subsets)"])
    ax.set_ylabel("final logdet "
                  r"$\log\det(J_{old\cup S}+\Lambda_x)$")
    ax.set_title("E5 acquisition policies (K=3 of 8 candidates); greedy "
                 "gap = %.4g" % policies["greedy_optimality_gap"])
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig("figures/e5_policies.png", dpi=160)
    plt.close(fig)


def _make_nonsubmodular_figure(nonsub):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    labels = ["g(empty)", "g(1)", "g(2)", "g(1)+g(2)", "g(1,2)"]
    for ax, r in zip(axes, nonsub["records"]):
        vals = [r["g_empty"], r["g1"], r["g2"], r["g1_plus_g2"],
                r["g12"]]
        colors = ["#8c8c8c", "#8c8c8c", "#8c8c8c", "#e09a00", "#1f4e79"]
        ax.bar(np.arange(5), vals, color=colors, width=0.68)
        ax.set_xticks(np.arange(5))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("value")
        ax.set_title(r"$\Lambda_x=%.0g$" % r["Lambda_x"])
        for xi, v in zip(np.arange(5), vals):
            ax.annotate(f"{v:.3f}", (xi, v), textcoords="offset points",
                        xytext=(0, 4), ha="center", fontsize=7)
        ax.grid(alpha=0.25, axis="y")
    axes[0].set_title("Non-submodularity gap = "
                      "%.4f" % nonsub["records"][0]["submodularity_gap"])
    axes[1].set_title("Non-submodularity gap = "
                      "%.4f" % nonsub["records"][1]["submodularity_gap"])
    fig.suptitle("E5 pose-Schur logdet: g(1)+g(2) < g(1,2)+g(empty) "
                 "violates submodularity", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig("figures/e5_nonsubmodular.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
