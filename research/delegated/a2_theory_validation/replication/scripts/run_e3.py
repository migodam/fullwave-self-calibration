#!/usr/bin/env python
"""Run experiment E3 and write all required results and figures.

Usage (from the experiment working directory):

    python scripts/run_e3.py
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
import matplotlib.ticker as mticker

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import e3


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

    data = e3.all_e3()
    dual = data["dual_spectra"]
    prior = data["pose_prior"]
    thm7 = data["theorem7"]
    caveat = data["adaptive_rank_caveat"]
    conf = data["confounding"]

    verdict = {
        "exact_small_dual_spectra": dual["exact_small"]["all_pass"],
        "random_nuisance_theorem3": dual["random_nuisance"]["pass"],
        "rho_undefined_control": dual["rho_undefined_control"]["pass"],
        "scaling_control": dual["scaling_control"]["pass"],
        "b_zero_control": dual["b_zero_control"]["pass"],
        "pose_prior_extension": prior["pass"],
        "theorem7_mc": thm7["pass"],
        "adaptive_rank_caveat_demonstrated": caveat["pass"],
        "confounding_sweep": conf["sweep_pass"],
    }
    doc = {
        "title": ("E3 dual tangent spectra, pose-prior extension, bias-aware "
                  "calibration risk, and relative-only gate failures"),
        "date": "2026-09-05",
        "package": "a2val",
        "settings": {
            "theory": ("Theorem 3 normalized dual spectra, Lambda_x pose-prior "
                       "Loewner/variational/augmented-angle identities, "
                       "Theorem 7 bias-aware risk, adaptive-rank caveat, "
                       "Section-10 confounding model"),
            "matrices": "real, unit-noise whitened, gauge-fixed",
            "support_tol": "K0/B0 eigen-directions below 1e-12*max are "
                           "reported undefined, never assigned 0 or 1",
            "random_fixture_seed_A2": 2150,
            "theorem7_fixture_seed": 2200,
            "theorem7_mc_seeds": list(range(301, 311)),
            "theorem7_mc_draws_per_seed": 1000,
            "adaptive_rank_seed": 2300,
            "adaptive_rank_draws": 500,
        },
        "results": data,
        "verdict": verdict,
    }

    with open("results/e3_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")

    rows = []

    def add(case, category, key_values, passed, note):
        rows.append({
            "case": case,
            "category": category,
            "key_values": key_values,
            "pass": bool(passed),
            "note": note,
        })

    es = dual["exact_small"]["records"]
    for label in ("pi/4", "pi/2"):
        r = es[label]
        add(f"exact_small_{label}", "dual_spectra",
            f"map={[_fmt(x, 4) for x in r['map_spectrum']]}, "
            f"pose={[_fmt(x, 4) for x in r['pose_spectrum']]}, "
            f"res={_fmt(max(r['map_residual'], r['pose_residual']), 2)}",
            r["pass"],
            "Theorem 3 spectra and multiplicities at known principal angles")
    rn = dual["random_nuisance"]
    add("random_nuisance", "dual_spectra",
        f"max_res={_fmt(rn['max_spectral_residual'], 2)}, "
        f"raw_Ke={[_fmt(x, 4) for x in rn['raw_K_e_eigs']]}, "
        f"norm_map={[_fmt(x, 4) for x in rn['map_spectrum_normalized']]}",
        rn["pass"],
        "Theorem 3 spectral identity on random nuisance fixture")
    ru = dual["rho_undefined_control"]
    add("rho_undefined_control", "dual_spectra",
        f"supp_map={ru['support_map']}, supp_pose={ru['support_pose']}, "
        f"undef_map={ru['undefined_map']}, undef_pose={ru['undefined_pose']}",
        ru["pass"],
        "null K0/B0 directions excluded, never rho 0 or 1")
    sc = dual["scaling_control"]
    add("scaling_control", "dual_spectra",
        "; ".join(
            f"s={r['scale']} rho_delta={_fmt(max(r['map_abs_delta_vs_unscaled'], r['pose_abs_delta_vs_unscaled']), 2)}, "
            f"J_s2={_fmt(r['J_x_fro_ratio_to_s2'], 3)}"
            for r in sc["records"]
        ),
        sc["pass"],
        "canonical rho invariant, J_x scales as s^2")
    bz = dual["b_zero_control"]
    add("b_zero_control", "dual_spectra",
        f"map_retention={bz['map_spectrum']}, J_x_sv={bz['J_x_singular_values']}",
        bz["pass"],
        "relative-only gate failure: map retention maximal with zero pose info")

    sand = prior["sandwich"]
    add("pose_prior_Loewner", "pose_prior_extension",
        f"min_eig(K_eL-K_e)={_fmt(sand['min_eig_K_eL-K_e'], 2)}, "
        f"min_eig(K0-K_eL)={_fmt(sand['min_eig_K0-K_eL'], 2)}",
        sand["K_e<=K_eL"] and sand["K_eL<=K0"],
        "K_e <= K_eL <= K0")
    add("pose_prior_monotone", "pose_prior_extension",
        "; ".join(f"{k}={_fmt(v, 2)}" for k, v in
                  prior["monotonicity"]["min_eig_differences"].items()),
        prior["monotonicity"]["pass"],
        "K_eL monotone nondecreasing in Lambda_x")
    add("pose_prior_variational", "pose_prior_extension",
        f"max_res={_fmt(prior['variational']['max_residual'], 2)}",
        prior["variational"]["pass"],
        "v^T K_eL v = min_h (||A_c v - B_c h||^2 + ||Lambda^1/2 h||^2)")
    aa = prior["augmented_principal_angles"]
    add("pose_prior_augmented_angles", "pose_prior_extension",
        f"aug_res={_fmt(aa['augmented_residual'], 2)}, "
        f"closed_form_res={_fmt(aa['closed_form_vs_augmented_residual'], 2)}",
        aa["augmented_residual"] <= 1e-10,
        "augmented A/B canonical correlations reproduce whitened K_eL")
    lz = prior["lambda_zero_recovery"]
    add("pose_prior_lambda_zero", "pose_prior_extension",
        f"fro_res={_fmt(lz['fro_residual_vs_K_e'], 2)}",
        lz["recovers_K_e"],
        "Lambda_x=0 recovers K_e")

    an = thm7["analytic"]
    add("theorem7_analytic", "theorem7",
        f"var_tr={_fmt(an['variance_trace'], 6)}, "
        f"bias2={_fmt(an['worst_bias_squared'], 6)}, "
        f"total={_fmt(an['analytic_total_risk'], 6)}",
        True,
        "analytic risk decomposition")
    mc = thm7["mc"]
    add("theorem7_mc_variance", "theorem7",
        f"max_rel={_fmt(mc['max_variance_rel_error'], 2)}, "
        f"within_3se={mc['all_variance_within_3se']}",
        mc["all_variance_within_3se"],
        "per-seed MC variance trace vs tr(M_x J_x^{-1})")
    add("theorem7_mc_bias", "theorem7",
        f"max_rel={_fmt(mc['max_bias_rel_error'], 2)}, "
        f"within_3se={mc['all_bias_within_3se']}",
        mc["all_bias_within_3se"],
        "per-seed MC total risk at worst SVD bias direction")
    add("theorem7_cov_check", "theorem7",
        f"cov_rel=[{_fmt(mc['min_sample_cov_rel_error'], 2)}, "
        f"{_fmt(mc['max_sample_cov_rel_error'], 2)}]",
        True,
        "sample covariance of B_v^+ eps vs J_x^{-1}, Frobenius relative")
    db = thm7["deterministic_bound"]
    add("theorem7_deterministic_bound", "theorem7",
        f"sigma_min={_fmt(db['sigma_min_B_v_Mx_invhalf'], 6)}, "
        f"bound={_fmt(db['bound'], 6)}, "
        f"worst_err={_fmt(db['worst_analytic_Mx_error'], 6)}",
        db["bound_satisfied"],
        "worst analytic M_x error <= (delta+beta_r)/sigma_min")
    add("adaptive_rank_caveat", "post_selection_caveat",
        f"sel_N1={caveat['n_selected_N1']}/{caveat['n_draws']}, "
        f"pool_vs_majority={_fmt(caveat['fro_diff_pooled_vs_majority_inverse_info'], 3)}",
        caveat["pass"],
        "pooled covariance of selected estimator != selected-model inverse info")

    cr = conf["records"]
    for row in cr:
        add(f"confounding_eps_{row['eps_param']}", "confounding",
            f"full={_fmt(row['full_ls_risk'], 6)}, "
            f"trunc={_fmt(row['truncated_risk'], 6)}, "
            f"better_trunc={row['truncation_better']}",
            row["truncation_better"] == (row["eps_param"] < 1.0),
            "full vs truncated risk")
    slb = conf["small_residual_large_bias"]
    add("confounding_small_residual_large_bias", "confounding",
        f"full={_fmt(slb['full_ls_risk_exact'], 8)}, "
        f"trunc={_fmt(slb['truncated_mse_exact'], 8)}, "
        f"trunc_bias={slb['truncation_bias']}, "
        f"residual_y2={_fmt(slb['residual_y2_under_truncation'], 3)}, "
        f"resid/pose_err={_fmt(slb['residual_y2_relative_to_pose_error'], 3)}",
        True,
        "residual-only gate passes while pose is wrong by c*")

    fieldnames = ["case", "category", "key_values", "pass", "note"]
    with open("results/e3_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    _make_dual_figure(dual)
    _make_prior_figure(prior)
    _make_confounding_figure(conf)
    _make_thm7_figure(thm7)

    print("E3 complete.")
    print("  verdict:", verdict)
    rn_max = dual["random_nuisance"]["max_spectral_residual"]
    print(f"  random Theorem-3 max spectral residual: {rn_max:.3e}")
    print(f"  pose-prior all pass: {prior['pass']}")
    print("  Theorem 7 analytic: variance =", _fmt(an["variance_trace"], 6),
          ", worst bias^2 =", _fmt(an["worst_bias_squared"], 6),
          ", total =", _fmt(an["analytic_total_risk"], 6))
    print("  Theorem 7 MC: max variance rel err =", _fmt(mc["max_variance_rel_error"], 2),
          ", max total-risk rel err =", _fmt(mc["max_bias_rel_error"], 2))
    print("  confounding crossover eps_param =", _fmt(conf["crossover_eps_param"], 6))
    print("  map retention vs pose bias:", conf["map_observation"]["map_information_retention"],
          "vs", conf["map_observation"]["pose_truncation_bias"])
    print("Saved: results/e3_results.json, results/e3_summary.csv")
    print("Saved: figures/e3_dual_spectra.png, figures/e3_prior_extension.png, "
          "figures/e3_bias_risk.png, figures/e3_thm7_mc.png")


def _make_dual_figure(dual):
    rn = dual["random_nuisance"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    map_x = np.arange(1, len(rn["map_spectrum_normalized"]) + 1)
    pose_x = np.arange(1, len(rn["pose_spectrum_normalized"]) + 1)
    canon_map = rn["expected_map_normalized"]
    canon_pose = rn["expected_pose_normalized"]
    axes[0].bar(map_x - 0.15, rn["map_spectrum_normalized"], 0.28,
                color="#1f4e79", label="normalized map spectrum")
    axes[0].scatter(map_x, rn["raw_K_e_eigs"], marker="x", s=45,
                    color="#b03a2e", label="raw K_e eigenvalues", zorder=3)
    axes[0].scatter(map_x + 0.15, canon_map, marker="o", s=38,
                    color="#e09a00", label="1-c^2 (+ units)", zorder=3)
    axes[0].set_xticks(map_x)
    axes[0].set_title("Map spectrum (a=%d): normalized vs raw K_e"
                      % len(map_x))
    axes[1].bar(pose_x - 0.15, rn["pose_spectrum_normalized"], 0.28,
                color="#1f4e79", label="normalized pose spectrum")
    axes[1].scatter(pose_x, rn["raw_J_x_eigs"], marker="x", s=45,
                    color="#b03a2e", label="raw J_x eigenvalues", zorder=3)
    axes[1].scatter(pose_x + 0.15, canon_pose, marker="o", s=38,
                    color="#e09a00", label="1-c^2 (+ units)", zorder=3)
    axes[1].set_xticks(pose_x)
    axes[1].set_title("Pose spectrum (b=%d): normalized vs raw J_x"
                      % len(pose_x))
    for ax in axes:
        ax.set_xlabel("sorted eigenvalue index")
        ax.set_ylabel("value")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25)
    fig.suptitle("E3 Theorem 3: whitened dual spectra match canonical "
                 "1-c^2+ones; raw spectra carry absolute scale",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("figures/e3_dual_spectra.png", dpi=160)
    plt.close(fig)


def _make_prior_figure(prior):
    sand = prior["sandwich"]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    data = {
        "K_e": np.sort(sand["K_e_eigs"]),
        "K_eL(0.5L)": np.sort(prior["monotonicity"]["records"][0]["eigenvalues"]),
        "K_eL(1.0L)": np.sort(prior["monotonicity"]["records"][1]["eigenvalues"]),
        "K_eL(2.0L)": np.sort(prior["monotonicity"]["records"][2]["eigenvalues"]),
        "K0": np.sort(sand["K0_eigs"]),
    }
    colors = ["#8c8c8c", "#7aa6d0", "#1f4e79", "#12314d", "#e09a00"]
    x = np.arange(1, len(next(iter(data.values()))) + 1)
    for (label, vals), color in zip(data.items(), colors):
        ax.plot(x, vals, marker="o", ls="-", ms=4.5, color=color,
                label=label)
    ax.set_xticks(x)
    ax.set_xlabel("sorted eigenvalue index")
    ax.set_ylabel("eigenvalue")
    ax.set_title("E3 pose prior: Loewner ordering\n"
                 r"$K_e \leq K_{eL}(0.5\Lambda)\leq K_{eL}(\Lambda) "
                 r"\leq K_{eL}(2\Lambda)\leq K_0$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig("figures/e3_prior_extension.png", dpi=160)
    plt.close(fig)


def _make_confounding_figure(conf):
    eps = np.asarray([r["eps_param"] for r in conf["records"]])
    full = np.asarray([r["full_ls_risk"] for r in conf["records"]])
    trunc = np.asarray([r["truncated_risk"] for r in conf["records"]])
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(eps, full, marker="o", color="#b03a2e", label="full LS risk")
    ax.plot(eps, trunc, marker="s", color="#1f4e79", label="truncated MSE")
    ax.axvline(conf["crossover_eps_param"], color="grey", ls=":",
               label=f"crossover e={conf['crossover_eps_param']:g}")
    ax.set_xscale("log")
    ax.set_xticks(eps)
    ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    ax.set_xlabel(r"$\epsilon_{\mathrm{param}}$ (nuisance coefficient scale)")
    ax.set_ylabel(r"risk / MSE of $h$ ($\sigma^2=1$, $c^*=1$)")
    ax.set_title("E3 confounding: truncation is better only for "
                 r"$\epsilon_{\mathrm{param}}<1$")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig("figures/e3_bias_risk.png", dpi=160)
    plt.close(fig)


def _make_thm7_figure(thm7):
    mc = thm7["mc"]
    recs = mc["records"]
    seeds = [r["seed"] for r in recs]
    est = [r["variance_estimate"] for r in recs]
    err = [3.0 * r["variance_se"] for r in recs]
    ana = thm7["analytic"]["variance_trace"]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.errorbar(seeds, est, yerr=err, fmt="o", ms=5, capsize=3,
                color="#1f4e79", label="MC estimate +/- 3 SE")
    ax.axhline(ana, color="#b03a2e", ls="--",
               label=f"analytic tr($M_x J_x^{{-1}}$) = {ana:.4f}")
    ax.set_xlabel("MC seed")
    ax.set_ylabel("variance trace estimate")
    ax.set_title("E3 Theorem 7: per-seed MC variance trace vs analytic value")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig("figures/e3_thm7_mc.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
