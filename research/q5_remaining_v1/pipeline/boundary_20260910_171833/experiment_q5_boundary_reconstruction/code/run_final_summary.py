#!/usr/bin/env python3
"""PART I -- FINAL aggregation for the bounded Q5 boundary reconstruction.

Consumes the authoritative artifacts produced by Parts A-H (never overwrites any
existing results JSON), emits:
    results/FINAL_SUMMARY.json   (machine-readable, one row per supplied constant)
    results/FINAL_REPORT.md      (<= 100 lines, human-readable)

No installs, single thread, no figures.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import sys
import time

import mpmath as mp
import numpy as np
import scipy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WALL0 = time.time()


def jload(rel):
    with open(os.path.join(ROOT, rel)) as fh:
        return json.load(fh)


def sha256(rel):
    p = os.path.join(ROOT, rel)
    if not os.path.exists(p):
        return None
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def newest(pat):
    hits = [p for p in glob.glob(os.path.join(ROOT, pat))
            if not p.endswith((".superseded", ".bad", ".tmp"))]
    hits.sort()
    return os.path.relpath(hits[-1], ROOT) if hits else None


def flat(o, p=""):
    out = {}
    if isinstance(o, dict):
        for k, v in o.items():
            out.update(flat(v, p + "/" + str(k)))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out.update(flat(v, p + "/%d" % i))
    else:
        out[p] = o
    return out


def main():
    cov_p = newest("results/coverage_audit_*.json")
    cex_p = newest("results/counterexample_*.json")
    cov = jload(cov_p)
    cex = jload(cex_p)
    mie_post_p = "results/mie_checks_20260910T101602Z.json"
    mie_pre_p = "results/mie_checks_20260910T092815Z.json"
    mie_post, mie_pre = jload(mie_post_p), jload(mie_pre_p)
    refgls = jload("results/ref_gls_20260910T100641Z.json")
    fixedborn = jload("results/fixed_loss_born_20260910T100202Z.json")
    interval = jload("results/interval_20260910T093408Z.json")
    qprime = jload("results/qprime_variants_20260910T094847Z.json")
    finrisk = jload("results/finite_risk_20260910T095306Z.json")
    signconv = jload("results/mie_sign_convention_20260910T100913Z.json")

    table = cov["constant_table"]
    counts = cov["counts"]

    # ---- Part G: independent pre/post fix diff on the three-way grid ------------
    fa, fb = flat(mie_pre), flat(mie_post)
    keys = set(fa) | set(fb)
    diff = sorted(k for k in keys if fa.get(k) != fb.get(k))
    meta_diff = [k for k in diff if k.startswith("/meta")]
    num_diff = [k for k in diff if not k.startswith("/meta")]
    agreement_keys = [
        "C_q_routes/max_abs_a_minus_b", "C_q_routes/max_rel_a_minus_b",
        "C_q_routes/max_rel_a_minus_c_consistent_pairing",
        "C_q_routes/max_rel_a_minus_c_prompt_literal",
        "C_q_routes/max_rel_a_minus_b_x_ge_0p15",
        "C_q_routes/max_rel_a_minus_c_x_ge_0p15",
    ]
    agree = {}
    for k in agreement_keys:
        pre = fa.get("/" + k)
        post = fb.get("/" + k)
        d = None
        if isinstance(pre, (int, float)) and isinstance(post, (int, float)):
            d = 0.0 if pre == post else abs(post - pre) / max(abs(pre), 1e-300)
        agree[k] = {"pre_fix": pre, "post_fix": post, "rel_change": d}
    grid_unchanged = all(v["rel_change"] == 0.0 for v in agree.values())
    other_secs = {s: sorted([k for k in num_diff if k.startswith("/" + s)])
                  for s in ("A_closed_forms", "B_derivative_fd", "C_q_routes",
                            "D_optical_and_rayleigh", "E_qprime_hprime_fd",
                            "F_supplied_constants")}

    # ---- Part H: counterexample summary ----------------------------------------
    h1 = cex["H1_admissible_lambda_interval"]
    h0 = cex["H0_theory_instance_raw"]
    h2 = cex["H2_no_rescale_counterexamples"]
    counterexamples = [{
        "q": r["q"], "q_prime": r["q_prime"],
        "abs_g": r["abs_g"], "abs_g_prime": r["abs_g_prime"],
        "both_gains_in_annulus": r["both_gains_in_annulus"],
        "eps_world1": r["eps_world1"], "eps_world2": r["eps_world2"],
        "eps_min_margin": min(r["material_margins"]["world1"]["min_margin_over_all"],
                              r["material_margins"]["world2"]["min_margin_over_all"]),
        "max_rel_data_residual": r["max_rel_residual"],
        "residual_B_identity": r["rel_residual_B_identity"],
        "note": r["note"],
    } for r in h2["instances"]]

    # ---- environment ------------------------------------------------------------
    prov_dir = "/Users/migodam/.cache/uv/archive-v0/VjSSQm31390egyfy"
    prov_list = "logs/uv_cache_VjSSQm31390egyfy.sha256"
    environment = {
        "interpreter": sys.executable,
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "mpmath_version": mp.__version__,
        "threads": 1,
        "mpmath_provenance": {
            "pythonpath_dir": prov_dir,
            "sha256_list_path": prov_list,
            "sha256_list_sha256": sha256(prov_list),
            "note": "mpmath 1.3.0 imported only via PYTHONPATH=%s; the sha256 list is the "
                    "frozen archive manifest captured at run start." % prov_dir,
        },
    }

    # ---- artifacts --------------------------------------------------------------
    art_patterns = [
        cov_p, cex_p, mie_post_p, mie_pre_p,
        "results/fixed_loss_born_20260910T100202Z.json",
        "results/ref_gls_20260910T100641Z.json",
        "results/mie_sign_convention_20260910T100913Z.json",
        "results/qprime_variants_20260910T094847Z.json",
        "results/finite_risk_20260910T095306Z.json",
        "results/interval_20260910T093408Z.json",
        "results/FINDINGS_increment2.md", "results/FINDINGS_increment3.md",
        "results/DERIVATIONS.md", "results/EXTRACTION_REPORT.md",
        "logs/coverage_audit.log", "logs/mie_fix.log", "logs/fixed_loss_born.log",
        "logs/mie_checks.log", "logs/mie_checks_prevrun_20260910T092815Z.log",
        "logs/qprime_variants.log", "logs/finite_risk.log", "logs/interval.log",
        "logs/ref_gls.log", "logs/mie_sign_convention.log",
        prov_list, "inputs/SOURCE_HASHES.json",
        "code/mie_reactance.py", "code/run_coverage_audit.py", "code/run_mie_fix_check.py",
        "code/run_counterexample.py", "code/run_final_summary.py",
    ]
    artifacts = []
    for p in art_patterns:
        if p is None:
            continue
        s = sha256(p)
        artifacts.append({"path": p, "sha256": s,
                          "bytes": os.path.getsize(os.path.join(ROOT, p)) if s else None})

    # ---- caveats ----------------------------------------------------------------
    qp_rows = {r["name"]: r for r in table}
    def rel_of(name):
        r = qp_rows.get(name)
        return None if r is None else r["rel_diff"]
    nuisance = fixedborn["C6_lower_bound"].get("extra_geometry_nuisance_probe", {})
    annulus = refgls["D2_logdet_and_constraint"]["constrained_gain_annulus"]
    logdet = refgls["D2_logdet_and_constraint"]["logdet_vs_gls"]
    c2 = fixedborn["C2_exact_reconstruction"]

    failures = [
        {
            "id": "qprime_lower_amplitude_derivative_lower_not_bit_reproduced",
            "detail": "supplied qprime_lower / amplitude_derivative_lower differ from every one "
                      "of the 7 interval arrangements tried (V1..V7); the supplied value lies "
                      "inside the V4/V6 envelope, so both are conservative lower bounds on the "
                      "same true infimum at the right interval endpoint.",
            "rel_diff_qprime_lower_I1": rel_of("monotonicity_checks[0].qprime_lower"),
            "rel_diff_qprime_lower_I2": rel_of("monotonicity_checks[1].qprime_lower"),
            "rel_diff_amplitude_lower_I1": rel_of("monotonicity_checks[0].amplitude_derivative_lower"),
            "rel_diff_amplitude_lower_I2": rel_of("monotonicity_checks[1].amplitude_derivative_lower"),
            "qprime_upper_rel_diff_I1": rel_of("monotonicity_checks[0].qprime_upper"),
            "qprime_upper_rel_diff_I2": rel_of("monotonicity_checks[1].qprime_upper"),
            "tightest_variant": "V4 fully expanded: +2.423% (I1), +3.755% (I2) vs supplied",
        },
        {
            "id": "section3_gain_instance_inadmissible",
            "detail": "A5 sec.3 quotes q=30, q'=45 with g=1; the implied partner gain has "
                      "|g'| = 0.6668722947 < 0.75, outside gain_allowed_modulus [0.75,1.25], "
                      "so that exact instance is NOT admissible under the declared gain domain. "
                      "It is repaired by a shared positive rescale lambda in [1.1246531096707642, 1.25].",
            "abs_g_prime": h0["abs_g_prime"],
            "gain_allowed_modulus": h0["gain_allowed_modulus"],
        },
        {
            "id": "logdet_vs_gls_distinction",
            "detail": "adding log det(sigma^2 I + sigma_r^2 f f^*) to the concentrated GLS "
                      "objective is a different (marginal) model, not the same profiling target; "
                      "it moves theta-hat by up to %.3e." % logdet["max_abs_theta_diff"],
            "max_abs_theta_diff": logdet["max_abs_theta_diff"],
            "median_abs_theta_diff": logdet["median_abs_theta_diff"],
        },
        {
            "id": "annulus_constrained_breakdown",
            "detail": "the unconstrained reference GLS gain |g*| = %.6f lies outside the declared "
                      "annulus; the constrained optimum sits on the boundary and the objective "
                      "jumps from %.4f to %.4f." % (annulus["abs_g_star"],
                                                    annulus["L_unconstrained_formula"],
                                                    annulus["L_constrained"]),
            "abs_g_star": annulus["abs_g_star"],
            "L_unconstrained": annulus["L_unconstrained_formula"],
            "L_constrained": annulus["L_constrained"],
        },
        {
            "id": "geometry_nuisance_breaks_sigma_min_Avis_bound",
            "detail": "the quoted sigma_min(A_vis) >= |g| sigma_min(B) sqrt(lambda_min(H)) bound "
                      "carries no geometry/position nuisance term; one labelled extra nuisance "
                      "direction drops the measured ratio to %.6f, i.e. the bound can fail."
                      % nuisance.get("ratio_with_extra", float("nan")),
            "ratio_with_extra": nuisance.get("ratio_with_extra"),
        },
        {
            "id": "ill_conditioned_LS_solver_artifact",
            "detail": "with cond(B)=1e4 the naive unweighted nonlinear least-squares re-fit loses "
                      "precision (only a few of 20 starts reach the global optimum) while the "
                      "explicit inverse stays exact; this is a scaling/solver artifact of the "
                      "check, not a failure of the sec.3 map.",
            "cond_B": 1e4,
        },
        {
            "id": "exact_reactance_identity_error_unidentified",
            "detail": "no (identity, precision) pair among (a)-(d) at float64 / dps=15,20,25,30,35,50 "
                      "landed within 20% of 8.879228496839643e-20 over the four world points; the "
                      "nearest is identity (a) in float64 at the two world-2 points (2.168409e-19, "
                      "2.44x). A dps=17 probe on a guessed 36-point grid reaches 0.90x but dps=17 is "
                      "outside the required set and the grid is a guess, so the constant stays "
                      "NOT_REPRODUCED.",
            "target": cov["exact_reactance_identity_error_scan"]["target"],
            "near_misses": cov["exact_reactance_identity_error_scan"]["near_misses"],
        },
        {
            "id": "superseded_and_quarantined_artifacts",
            "detail": "results/interval_20260910T093248Z.json is truncated JSON and was quarantined "
                      "as PARTIAL_interval_20260910T093248Z.json.bad; Part D ref_gls T100317Z "
                      "(Woodbury numerator bug) and Part E mie_sign_convention T100743Z "
                      "(wrong chi1 sign) are superseded and must not be cited; "
                      "counterexample T101904Z and T102008Z are two identical seeded runs.",
        },
        {
            "id": "part_B_q_mie_offset_diagnostics_are_pre_fix",
            "detail": "the frozen Part B JSON (finite_risk T095306Z) records two *diagnostic* "
                      "leaves q_mie_xi_minus_vs_q_rel (1.9999947 / 1.9999976) and "
                      "q_mie_xi_plus_vs_q_rel (2.6664e-3 / 1.8048e-3) produced by the OLD "
                      "wrong-arrangement q_mie(); after the Part G fix the same diagnostics read "
                      "~5.3e-12. No scientific number in Part B (scale c, delta_mu identity, "
                      "sigma factor, D, p_*, error bounds) depends on them, and a re-run was "
                      "deliberately NOT issued so the frozen artifact stays the pre-fix record.",
            "pre_fix_values": {"q_mie_xi_minus_vs_q_rel": [1.9999946677220815, 1.9999975570309303],
                               "q_mie_xi_plus_vs_q_rel": [0.0026664010227751404,
                                                          0.001804796416290883]},
        },
    ]

    hypothesis = (
        "Every numeric constant supplied with the round is either (i) reproduced by an "
        "independent computation to <= 1e-9 relative, (ii) identified as a provenance/array "
        "field with no independent numeric target, or (iii) explicitly reported as not "
        "reproduced with the reason; and the A5 sec.3 parallel u/ell counterexample is "
        "admissible under the round's declared gain annulus and material intervals."
    )
    overall_verdict = (
        "PARTIAL (bounded Class M known-geometry modal check only). %d supplied constants were "
        "compared: %d REPRODUCED, %d REPRODUCED_WEAKER, %d NOT_REPRODUCED, %d INCOMPARABLE. "
        "The A5 sec.3 raw instance (q=30, q'=45, g=1) is INADMISSIBLE (|g'| = %.9f < 0.75) but "
        "is repairable by a shared rescale lambda in [%s, %s]; three further no-rescale "
        "counterexamples are admissible and reproduce the data to <= 1.7e-16. The q_mie() "
        "wrong-arrangement default was replaced by the branch-correct arrangement, with the "
        "three-way cross-implementation grid unchanged (max relative change 0)."
        % (counts["constants_compared"], counts["REPRODUCED"], counts["REPRODUCED_WEAKER"],
           counts["NOT_REPRODUCED"], counts["INCOMPARABLE"], h0["abs_g_prime"],
           h1["lambda_min"], h1["lambda_max"])
    )

    summary = {
        "part": "I",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "hypothesis": hypothesis,
        "overall_verdict": overall_verdict,
        "counts": {
            "leaves_in_inputs_modal_boundary_json": cov["leaf_enumeration"]["total_leaves"],
            "leaves_by_category": cov["leaf_enumeration"]["counts_by_category"],
            "constants_compared": counts["constants_compared"],
            "REPRODUCED": counts["REPRODUCED"],
            "REPRODUCED_WEAKER": counts["REPRODUCED_WEAKER"],
            "NOT_REPRODUCED": counts["NOT_REPRODUCED"],
            "INCOMPARABLE": counts["INCOMPARABLE"],
            "not_reproduced_names": counts["not_reproduced_names"],
            "reproduced_weaker_names": counts["reproduced_weaker_names"],
            "incomparable_names": counts["incomparable_names"],
        },
        "constant_table": [{
            "name": r["name"], "supplied": r.get("supplied"),
            "computed": r.get("computed"), "rel_diff": r.get("rel_diff"),
            "verdict": r.get("verdict"), "abs_diff": r.get("abs_diff"),
            "method": r.get("method"),
        } for r in table],
        "part_F_coverage_audit": {
            "source_json": cov_p,
            "leaf_enumeration_counts": cov["leaf_enumeration"]["counts_by_category"],
            "relative_data_difference": {
                "supplied": 0.0003036059493932306,
                "winner": cov["relative_data_difference_winner"],
                "winner_value": cov["relative_data_difference_candidates"][
                    cov["relative_data_difference_winner"]],
                "candidates": cov["relative_data_difference_candidates"],
            },
            "exact_reactance_identity_error": {
                "identified": cov["exact_reactance_identity_error_scan"]["identified"],
                "target": cov["exact_reactance_identity_error_scan"]["target"],
                "matches_within_20pct": cov["exact_reactance_identity_error_scan"][
                    "matches_within_20pct_at_required_precisions"],
                "near_misses": cov["exact_reactance_identity_error_scan"]["near_misses"],
            },
        },
        "part_G_mie_fix": {
            "pre_fix_json": mie_pre_p, "post_fix_json": mie_post_p,
            "differing_leaves_across_all_sections": len(diff),
            "differing_leaves_in_numerical_sections": len(num_diff),
            "differing_meta_leaves": meta_diff,
            "differing_leaves_by_section": other_secs,
            "three_way_agreement_pre_vs_post": agree,
            "three_way_grid_unchanged": bool(grid_unchanged),
            "convention": ("outgoing branch h1^(1)(z) = psi1 - i*chi1 = z*h1^(1)(z) "
                           "(xi_sign=-1) is the h^(1) branch tied to e^{-i omega t}; on it "
                           "Im(a1) < 0 (Re a1 = +|a1|^2 for lossless eps) and the exact "
                           "conversion to the A5 sec.8 reactance is q = +i*a1/(1-a1). The module "
                           "default xi_sign=+1 is the conjugated ingoing h1^(2) branch, "
                           "Im(a1) > 0, Re a1 = -|a1|^2, with q = -i*a1/(1-a1)."),
            "legacy_wrong_arrangement": "-1j*a1/(1+a1) (prompt literal) differs from q by rel "
                                        "2.67e-3 on xi_sign=+1; retained only as "
                                        "q_mie_legacy_prompt_literal() and the named form raises "
                                        "ValueError.",
        },
        "part_H_counterexample": {
            "source_json": cex_p,
            "raw_theory_instance": {
                "q": h0["q"], "q_prime": h0["q_prime"], "abs_g": h0["abs_g"],
                "abs_g_prime": h0["abs_g_prime"],
                "gain_allowed_modulus": h0["gain_allowed_modulus"],
                "both_gains_in_annulus": h0["both_gains_in_annulus"],
                "max_rel_data_residual": h0["max_rel_residual"],
                "eps_world1": h0["eps_world1"], "eps_world2": h0["eps_world2"],
            },
            "admissible_lambda_interval": {
                "lambda_min": h1["lambda_min"], "lambda_max": h1["lambda_max"],
                "lambda_min_float": h1["lambda_min_float"],
                "world1_constraint": h1["constraint_gain_world1_lambda_range"],
                "world2_constraint": h1["constraint_gain_world2_lambda_range"],
                "lambda_1p2_check": h1["lambda_equals_1p2_check"],
            },
            "no_rescale_counterexamples": counterexamples,
            "all_admissible": h2["all_admissible"],
            "all_eps_inside_material_domain": h2["all_eps_inside_material_domain"],
            "all_delta_mu_zero_to_1e15": h2["all_delta_mu_zero_to_1e15"],
            "verdicts": cex["H4_verdicts"],
        },
        "other_parts_numbers": {
            "A_qprime_variants": qprime.get("summary"),
            "B_finite_risk": finrisk.get("summary"),
            "C_fixed_loss_born": {
                "det_identity_max_rel_err": {
                    k: fixedborn["C2_exact_reconstruction"][k].get("max_det_identity_rel_err")
                    for k in ("m6_well_conditioned", "m12_well_conditioned",
                              "m6_ill_conditioned_1e4", "m12_ill_conditioned_1e4")},
                "C4": fixedborn["C4_parallel_counterexample"]["same_data_rel_residual"],
                "C5_slopes": {k: v.get("slope_log_adv_du_vs_log_absD")
                              for k, v in fixedborn["C5_near_parallel"].items()},
                "C6_lambda_min_H_max_abs_err": fixedborn["C6_lower_bound"]["lambda_min_H"]["max_abs_err"],
                "C6_ratio_B_identity": [fixedborn["C6_lower_bound"]["sigma_min_Avis_bound"]["B_identity_m2"]["min_ratio"],
                                        fixedborn["C6_lower_bound"]["sigma_min_Avis_bound"]["B_identity_m2"]["max_ratio"]],
                "C6_extra_geometry_nuisance_ratio": nuisance.get("ratio_with_extra"),
            },
            "D_ref_gls": {
                "logdet_max_abs_theta_diff": logdet["max_abs_theta_diff"],
                "abs_g_star": annulus["abs_g_star"],
                "L_unconstrained": annulus["L_unconstrained_formula"],
                "L_constrained": annulus["L_constrained"],
                "max_rel_err_full_derivative": refgls["D4_residual_derivative"]["max_rel_err_full_derivative"],
                "max_rel_err_naive_drop_term": refgls["D4_residual_derivative"]["max_rel_err_naive_drop_term"],
            },
            "E_sign_convention": signconv["PART_E_sign_convention"],
            "interval": {"m_i_reproduction": interval.get("m_i_reproduction"),
                         "comparison_with_supplied": interval.get("comparison_with_supplied")},
        },
        "artifacts": artifacts,
        "environment": environment,
        "declared_limits": [
            "this is the restricted Class M known-geometry modal check only",
            "not joint two-sphere recovery",
            "not Class C",
            "no hardware",
            "no novelty claim",
            "interval verification is machine-assisted and not independently re-implemented by a second party",
        ],
        "failures_and_caveats": failures,
    }

    out_p = os.path.join(ROOT, "results", "FINAL_SUMMARY.json")
    tmp = out_p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    os.replace(tmp, out_p)
    print("WROTE", os.path.relpath(out_p, ROOT))

    # ------------------------------------------------------------------ report
    n_rep = counts["REPRODUCED"]
    lines = []
    lines.append("# FINAL REPORT -- bounded Q5 boundary reconstruction (Parts A-I)")
    lines.append("")
    lines.append("cwd `%s`" % ROOT)
    lines.append("interpreter `%s`; numpy %s, scipy %s, mpmath %s (PYTHONPATH `%s`)."
                 % (sys.executable, np.__version__, scipy.__version__, mp.__version__, prov_dir))
    lines.append("Machine-readable table: `results/FINAL_SUMMARY.json` (%d constant rows)."
                 % len(table))
    lines.append("")
    lines.append("## Hypothesis")
    lines.append(hypothesis)
    lines.append("")
    lines.append("## Overall verdict")
    lines.append(overall_verdict)
    lines.append("")
    lines.append("## Coverage (Part F)")
    lines.append("- `inputs/modal_boundary.json` has %d leaves: category i=%d, ii=%d, iii=%d, iv=%d."
                 % (cov["leaf_enumeration"]["total_leaves"],
                    cov["leaf_enumeration"]["counts_by_category"]["i"],
                    cov["leaf_enumeration"]["counts_by_category"]["ii"],
                    cov["leaf_enumeration"]["counts_by_category"]["iii"],
                    cov["leaf_enumeration"]["counts_by_category"]["iv"]))
    lines.append("- Constants compared: %d; REPRODUCED %d, REPRODUCED_WEAKER %d, NOT_REPRODUCED %d, INCOMPARABLE %d."
                 % (counts["constants_compared"], counts["REPRODUCED"],
                    counts["REPRODUCED_WEAKER"], counts["NOT_REPRODUCED"], counts["INCOMPARABLE"]))
    lines.append("- NOT_REPRODUCED: %s." % ", ".join("`%s`" % n for n in counts["not_reproduced_names"]))
    lines.append("- REPRODUCED_WEAKER: %s." % ", ".join("`%s`" % n for n in counts["reproduced_weaker_names"]))
    lines.append("- INCOMPARABLE: %s." % ", ".join("`%s`" % n for n in counts["incomparable_names"]))
    wr = cov["relative_data_difference_winner"]
    lines.append("- `relative_data_difference` = 3.036059493932306e-04 REPRODUCED by candidate "
                 "`%s` = %s at dps=50 (rel 1.69e-15); the other four normalizations do not match."
                 % (wr, cov["relative_data_difference_candidates"][wr]))
    lines.append("- `exact_reactance_identity_error` = 8.879228496839643e-20 NOT IDENTIFIED; "
                 "nearest probe 2.168409e-19 (2.44x) identity (a) float64 world-2.")
    lines.append("")
    lines.append("## Mie helper fix (Part G)")
    lines.append("- `code/mie_reactance.py` sha256 post-fix `%s` (pre-fix `%s`)."
                 % (sha256("code/mie_reactance.py"), "6bdea5e475d5b4176396b59ebc5d4da94b190885dec381bc12075920b463b62e"))
    lines.append("- Three-way grid (route a closed form / b scipy / c Mie a1): pre-fix vs post-fix "
                 "differing leaves in any numerical section = %d (of %d total; the other %d are "
                 "/meta timestamps + the script hash); max relative change in every agreement "
                 "figure = 0." % (len(num_diff), len(diff), len(meta_diff)))
    lines.append("  - max rel |q_a-q_b| = %.6e ; max rel |q_a-q_c| (consistent pairing) = %.6e ; "
                 "prompt-literal = %.6e (unchanged, the legacy value is retained only as a record)."
                 % (agree["C_q_routes/max_rel_a_minus_b"]["post_fix"],
                    agree["C_q_routes/max_rel_a_minus_c_consistent_pairing"]["post_fix"],
                    agree["C_q_routes/max_rel_a_minus_c_prompt_literal"]["post_fix"]))
    lines.append("- Route c in that grid is generated by `q_from_a1_conversion(..., xi_sign=-1, "
                 "form='i_over_1minus')`, i.e. it never called `q_mie()`, which is why the fix "
                 "cannot move it; the only pre-fix `q_mie()` numbers left in the frozen artifacts "
                 "are the two Part B offset diagnostics (1.99999 and 2.67e-3/1.80e-3).")
    lines.append("- Final convention: outgoing `h1^(1)(z) = psi1 - i*chi1 = z*h1^(1)(z)` "
                 "(`xi_sign=-1`) pairs with `e^{-i omega t}`, `Im(a1) < 0`, `Re a1 = +|a1|^2`, "
                 "exact conversion `q = +i*a1/(1-a1)`; the module default `xi_sign=+1` is the "
                 "conjugated ingoing `h1^(2)` branch with `q = -i*a1/(1-a1)`.")
    lines.append("")
    lines.append("## Admissible parallel counterexample (Part H)")
    lines.append("- Raw A5 sec.3 instance q=30, q'=45, g=1: eps worlds (1.9,2.5) and (2.35,3.25) "
                 "inside I_1=[1.5,4], I_2=[2,5] with margins (0.400/2.100, 0.500/2.500) and "
                 "(0.850/1.650, 0.350/1.750); data residual %.3e over 9 random complex B + B=I; "
                 "BUT |g'| = %.10f < 0.75 -> INADMISSIBLE under `gain_allowed_modulus` [0.75,1.25]."
                 % (h0["max_rel_residual"], h0["abs_g_prime"]))
    lines.append("- Repair by a shared positive rescale g -> lambda g, g' -> lambda g': "
                 "lambda in [%s, %s] (lambda_min = 0.75/|g'| exactly); lambda = 1.2 gives "
                 "|g| = 1.2, |g'| = %.10f with B=I residual %.3e."
                 % (h1["lambda_min"], h1["lambda_max"],
                    h1["lambda_equals_1p2_check"]["abs_g_prime"],
                    h1["lambda_equals_1p2_check"]["rel_residual_B_identity"]))
    for r in counterexamples:
        def fmt_eps(vals, doms):
            return ", ".join("%.4g [margin %.3g/%.3g]"
                             % (v, v - doms[i][0], doms[i][1] - v)
                             for i, v in enumerate(vals))
        lines.append("- No-rescale admissible: q=%g, q'=%g, |g'|=%.12f, both gains in annulus %s; "
                     "eps1=(%s) in [1.5,4]x[2,5], eps2=(%s); max data residual %.3e."
                     % (r["q"], r["q_prime"], r["abs_g_prime"], r["both_gains_in_annulus"],
                        fmt_eps(r["eps_world1"], [[1.5, 4.0], [2.0, 5.0]]),
                        fmt_eps(r["eps_world2"], [[1.5, 4.0], [2.0, 5.0]]),
                        r["max_rel_data_residual"]))
    lines.append("- All admissible instances satisfy delta_mu = 0 to rel <= 1e-15 "
                 "(`all_delta_mu_zero_to_1e15` = %s) and lie strictly inside the material domain."
                 % h2["all_delta_mu_zero_to_1e15"])
    lines.append("")
    lines.append("## Other parts, exact numbers")
    lines.append("- Part A q' sweep: no variant reproduces the supplied lower bounds; tightest "
                 "rigorous V4 = +2.423% (I1) / +3.755% (I2) vs supplied; infimum at the right "
                 "endpoint (eps=4.0, eps=5.0).")
    lines.append("- Part C: det M = -|g|^2 det[u,ell] verified (max rel 1.56e-13 / 1.82e-13); exact "
                 "reconstruction max rel 1.56e-13 (m6) / 1.82e-13 (m12); near-parallel slopes "
                 "-1.0052 / -0.9978 / -0.9986; lambda_min(H) max abs err %.3e."
                 % fixedborn["C6_lower_bound"]["lambda_min_H"]["max_abs_err"])
    lines.append("- Part D: scalar GLS identity rel 1.34e-16; log-det objective moves theta by up "
                 "to %.3e; |g*| = %.4f outside the annulus (L %.1f -> %.1f); full residual "
                 "derivative rel %.1e vs mpmath dps=60, naive drop of the term errs %.4f."
                 % (logdet["max_abs_theta_diff"], annulus["abs_g_star"],
                    annulus["L_unconstrained_formula"], annulus["L_constrained"],
                    refgls["D4_residual_derivative"]["max_rel_err_full_derivative"],
                    refgls["D4_residual_derivative"]["max_rel_err_naive_drop_term"]))
    lines.append("- Part E: textbook a1 (h1^(1)) = 1.7774330813186167e-06 - 0.0013332028810538396i; "
                 "`+i*a1/(1-a1)` reproduces q to rel 2.1e-15, `-i*a1/(1+a1)` errs rel 2.67e-3.")
    lines.append("")
    lines.append("## Declared limits")
    for d in summary["declared_limits"]:
        lines.append("- %s." % d)
    lines.append("")
    lines.append("## Failures and caveats")
    for f in failures:
        lines.append("- `%s`: %s" % (f["id"], " ".join(f["detail"].split())))
    lines.append("")
    lines.append("## Artifacts (path -- sha256)")
    for a in artifacts:
        lines.append("- `%s` -- `%s`" % (a["path"], (a["sha256"] or "")[:32]))
    rep_p = os.path.join(ROOT, "results", "FINAL_REPORT.md")
    tmp = rep_p + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    os.replace(tmp, rep_p)
    print("WROTE results/FINAL_REPORT.md lines=%d wall=%.2fs" % (len(lines), time.time() - WALL0))


if __name__ == "__main__":
    main()
