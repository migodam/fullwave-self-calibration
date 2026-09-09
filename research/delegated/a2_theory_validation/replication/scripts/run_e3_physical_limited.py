#!/usr/bin/env python
"""E3 physical extension, FULL vs LIMITED aperture, with and without pose
priors (N=8 smoke, seeds 301..310, m=1000 Theorem-7 MC draws).

Usage (from the experiment working directory):

    python scripts/run_e3_physical_limited.py

The shared physical core is imported read-only from
``research/delegated/a2_physics``; nothing under that directory is modified.
All current-space directions are DECLARED candidate left singular directions
of the realified full data tangent (no oracle rank claim).
"""

from __future__ import annotations

import csv
import json
import os
import sys
import traceback

import numpy as np
import scipy.stats as st
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PHYSICS_DIR = (
    "/Volumes/migodam's-external-brain/Research/Inv_SLAM/"
    "research/delegated/a2_physics"
)
if PHYSICS_DIR not in sys.path:
    sys.path.insert(0, PHYSICS_DIR)
import physics  # noqa: E402  (read-only shared core)


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import e2  # noqa: E402
from a2val import e3  # noqa: E402
from a2val import physical_runner as pr  # noqa: E402


EPS = float(np.finfo(float).eps)
SEEDS = list(range(301, 311))
APERTURES = ("full", "limited")
N_GRID = 8
M_DRAWS = 1000
DF = 9
T_CRIT = float(st.t.ppf(0.975, DF))


def _sym(X):
    X = np.asarray(X, dtype=float)
    return 0.5 * (X + X.T)


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


def _fmt(v, nd=8):
    if v is None:
        return ""
    return f"{v:.{nd}g}"


def _reference(seed: int):
    """Deterministic non-oracle reference (clamped as specified)."""
    rng = np.random.default_rng(seed)
    alpha = 0.2 + 0.05 * rng.standard_normal(9)
    alpha = np.clip(alpha, -0.5, 0.5)
    x = 0.02 * rng.standard_normal(3)
    x = np.clip(x, -0.05, 0.05)
    return alpha, x, rng


def _stats(values):
    values = [float(v) for v in values]
    n = len(values)
    if n == 0:
        return {
            "n": 0,
            "mean": None,
            "sd": None,
            "ci_low": None,
            "ci_high": None,
            "ci_half": None,
        }
    arr = np.asarray(values, dtype=float)
    mean = float(np.mean(arr))
    sd = float(np.std(arr, ddof=1))
    half = T_CRIT * sd / np.sqrt(float(n))
    return {
        "n": int(n),
        "mean": mean,
        "sd": sd,
        "ci_low": mean - half,
        "ci_high": mean + half,
        "ci_half": half,
    }


def _pair_stats(diffs):
    diffs = [float(v) for v in diffs]
    if not diffs:
        return {
            "n_pairs": 0,
            "mean": None,
            "sd": None,
            "ci_low": None,
            "ci_high": None,
            "ci_half": None,
            "ci_strictly_positive": False,
        }
    s = _stats(diffs)
    s["n_pairs"] = s.pop("n")
    s["ci_strictly_positive"] = bool(s["ci_low"] > 0.0)
    return s


def _build_theorem7_record(R, Jx_min_eig, seed, rng) -> dict:
    """Theorem-7 analytic/MC block following the declared-bias convention."""
    n = int(R.shape[0])
    full_rank = bool(Jx_min_eig > 1e-10)
    if not full_rank:
        return {
            "full_rank_B_v": False,
            "min_eig_J_x": float(Jx_min_eig),
            "bias_direction_source": None,
            "note": (
                "min_eig(J_x) <= 1e-10, so the Theorem-7 inverse risk is not "
                "evaluated (rank-deficient visible residual B_v); no exception "
                "raised."
            ),
        }

    R = np.asarray(R, dtype=float)
    u, s, _ = np.linalg.svd(R, full_matrices=False)
    bias_source = "first_left_singular_of_residual"
    if s.size == 0 or s[0] <= max(1e-14, EPS * float(np.linalg.norm(R, "fro"))):
        u0 = rng.standard_normal(n)
        u0 = u0 / max(float(np.linalg.norm(u0)), np.finfo(float).tiny)
        bias_source = "random_unit_fallback_residual_zero"
    else:
        u0 = u[:, 0]
    D_r = np.asarray(u0, dtype=float).reshape(-1, 1)
    fx = {
        "B_v": np.asarray(R, dtype=float),
        "D_r": D_r,
        "M_x": np.eye(3, dtype=float),
        "n": int(n),
    }
    ana = e3.theorem7_analytic(fx)
    mc = e3.theorem7_mc_seed(fx, int(seed), m=int(M_DRAWS))
    return {
        "full_rank_B_v": True,
        "min_eig_J_x": float(Jx_min_eig),
        "bias_direction_source": bias_source,
        "D_r_norm": float(np.linalg.norm(D_r)),
        "analytic": {
            "variance_trace": float(ana["variance_trace"]),
            "worst_bias_squared": float(ana["worst_bias_squared"]),
            "analytic_total_risk": float(ana["analytic_total_risk"]),
        },
        "mc": _jsonable(mc),
    }


def _run_one(model, seed: int, aperture: str) -> dict:
    alpha, x, rng = _reference(int(seed))
    out = model.forward(alpha, x, jacobian=True)
    A_r = physics.realify_jacobian(out["A"])
    B_r = physics.realify_jacobian(out["B"])

    T = np.hstack([A_r, B_r])
    U, sT, _ = np.linalg.svd(T, full_matrices=False)
    C1 = U[:, :4]

    # Theorem 3 spectral identity at the declared candidate space.
    thm3 = pr.theorem3_physical(A_r, B_r, C1)

    # Residual visible pose tangent and its information Gram.
    R = e2.Bv(A_r, B_r, C1)
    Jx = _sym(R.T @ R)
    jx_eig = np.sort(np.linalg.eigvalsh(Jx))
    pose_info_trace = float(np.sum(jx_eig))
    pose_info_min_eig = float(jx_eig[0]) if jx_eig.size else 0.0

    # Pose-prior Loewner sandwich checks.
    priors = {}
    for lam, key in ((0.0, "lambda0"), (0.1, "lambda01"), (1.0, "lambda1")):
        pr_rec = pr.prior_sandwich_physical(
            A_r, B_r, C1, lam * np.eye(3, dtype=float)
        )
        priors[key] = {
            "lambda_x": lam,
            "loewner_sandwich_ok": bool(pr_rec["loewner_sandwich_ok"]),
            "min_eig_K_eL_minus_K_e": float(
                pr_rec["min_eig_K_eL_minus_K_e"]
            ),
            "min_eig_K0_minus_K_eL": float(pr_rec["min_eig_K0_minus_K_eL"]),
            "K_e_eigenvalues": pr_rec["K_e_eigenvalues"],
            "K_eL_eigenvalues": pr_rec["K_eL_eigenvalues"],
            "K0_eigenvalues": pr_rec["K0_eigenvalues"],
        }

    thm7 = _build_theorem7_record(R, pose_info_min_eig, seed, rng)

    return {
        "seed": int(seed),
        "aperture": aperture,
        "work_wall_seconds": float(out["work"]["wall_seconds"]),
        "reference": {
            "alpha": alpha.tolist(),
            "x": x.tolist(),
            "note": (
                "alpha=0.2+0.05*N(0,1) clipped to [-0.5,0.5]; "
                "x=0.02*N(0,1) clipped to [-0.05,0.05]; "
                "rng=np.random.default_rng(seed)"
            ),
        },
        "declared_C1": {
            "width": 4,
            "source": (
                "first 4 left singular vectors of realified full data tangent "
                "T=hstack([A_r,B_r]); declared candidate current directions, "
                "NO oracle rank"
            ),
            "T_singular_values_top": np.sort(sT)[::-1][:8].tolist(),
        },
        "Jx": {
            "eigenvalues_ascending": jx_eig.tolist(),
            "pose_info_trace": pose_info_trace,
            "pose_info_min_eig": pose_info_min_eig,
            "source": (
                "R=(I-P_[C1,A])B via a2val.e2.Bv; J_x=R^T R; eigenvalues "
                "ascending via np.linalg.eigvalsh"
            ),
        },
        "theorem3": {
            "a": thm3["a"],
            "b": thm3["b"],
            "canonical_c": thm3["canonical_c"],
            "rank_Ccan": thm3["rank_Ccan"],
            "map_normalized_spectrum": thm3["map_normalized_spectrum"],
            "pose_normalized_spectrum": thm3["pose_normalized_spectrum"],
            "map_expected_spectrum": thm3["map_expected_spectrum"],
            "pose_expected_spectrum": thm3["pose_expected_spectrum"],
            "map_support_dim": thm3["map_support_dim"],
            "pose_support_dim": thm3["pose_support_dim"],
            "map_residual_max_abs": float(thm3["map_residual_max_abs"]),
            "pose_residual_max_abs": float(thm3["pose_residual_max_abs"]),
            "max_residual_abs": float(
                max(thm3["map_residual_max_abs"], thm3["pose_residual_max_abs"])
            ),
            "raw_K0_eigenvalues": thm3["raw_K0_eigenvalues"],
            "raw_K_e_eigenvalues": thm3["raw_K_e_eigenvalues"],
            "raw_J_x_eigenvalues": thm3["raw_J_x_eigenvalues"],
        },
        "priors": priors,
        "theorem7": thm7,
    }


def _by_aperture(records):
    by = {ap: {} for ap in APERTURES}
    for rec in records:
        by[rec["aperture"]][int(rec["seed"])] = rec
    return by


def _map_retention(rec):
    ev = np.asarray(rec["theorem3"]["raw_K_e_eigenvalues"], dtype=float)
    return float(np.sum(ev)), float(np.min(ev))


def _theorem3_max_res(rec):
    return float(rec["theorem3"]["max_residual_abs"])


def _mc_block(rec):
    th7 = rec["theorem7"]
    if not th7.get("full_rank_B_v") or "mc" not in th7:
        return None
    return th7["mc"]


def _metric_series(records, name):
    out = []
    for rec in records:
        if name == "pose_info_trace":
            v = rec["Jx"]["pose_info_trace"]
        elif name == "pose_info_min_eig":
            v = rec["Jx"]["pose_info_min_eig"]
        elif name == "map_retention_trace":
            v, _ = _map_retention(rec)
        elif name == "map_retention_min_eig":
            _, v = _map_retention(rec)
        elif name == "theorem3_max_residual":
            v = _theorem3_max_res(rec)
        elif name in (
            "theorem7_variance_rel_error",
            "theorem7_bias_rel_error",
            "sample_cov_vs_Jx_inv_fro_rel",
        ):
            mc = _mc_block(rec)
            if mc is None:
                continue
            if name == "theorem7_variance_rel_error":
                v = mc["variance_rel_error"]
            elif name == "theorem7_bias_rel_error":
                v = mc["bias_rel_error"]
            else:
                v = mc["sample_cov_vs_Jx_inv_fro_rel"]
        else:
            raise ValueError(f"unknown metric {name}")
        out.append((rec["seed"], float(v)))
    return out


def _per_aperture_metric_stats(records_by_ap, name):
    out = {}
    for ap in APERTURES:
        series = _metric_series(records_by_ap[ap].values(), name)
        if name in (
            "theorem7_variance_rel_error",
            "theorem7_bias_rel_error",
            "sample_cov_vs_Jx_inv_fro_rel",
        ):
            stats = _stats([v for _, v in series]) if series else _stats([])
            stats["n_available"] = len(series)
        else:
            stats = _stats([v for _, v in series])
            stats["n_available"] = len(series)
        stats["seeds"] = [int(s) for s, _ in series]
        out[ap] = stats
    return out


def _paired_metric_stats(records_by_ap, name):
    pairs = []
    for seed in SEEDS:
        f = records_by_ap["full"].get(seed)
        l = records_by_ap["limited"].get(seed)
        if f is None or l is None:
            continue
        if name in (
            "theorem7_variance_rel_error",
            "theorem7_bias_rel_error",
            "sample_cov_vs_Jx_inv_fro_rel",
        ):
            mf = _mc_block(f)
            ml = _mc_block(l)
            if mf is None or ml is None:
                continue
            if name == "theorem7_variance_rel_error":
                pairs.append((seed, mf["variance_rel_error"]
                              - ml["variance_rel_error"]))
            elif name == "theorem7_bias_rel_error":
                pairs.append((seed, mf["bias_rel_error"]
                              - ml["bias_rel_error"]))
            else:
                pairs.append((seed, mf["sample_cov_vs_Jx_inv_fro_rel"]
                              - ml["sample_cov_vs_Jx_inv_fro_rel"]))
            continue
        sf = _metric_series([f], name)
        sl = _metric_series([l], name)
        if not sf or not sl:
            continue
        pairs.append((seed, float(sf[0][1]) - float(sl[0][1])))
    stats = _pair_stats([v for _, v in pairs])
    stats["seeds"] = [int(s) for s, _ in pairs]
    return stats


def _make_figure(records_by_ap) -> None:
    left_metrics = ["pose_info_trace", "map_retention_trace"]
    left_means = {ap: {} for ap in APERTURES}
    left_errs = {ap: {} for ap in APERTURES}
    for ap in APERTURES:
        recs = list(records_by_ap[ap].values())
        for name in left_metrics:
            vals = [float(v) for _, v in _metric_series(recs, name)]
            s = _stats(vals)
            left_means[ap][name] = s["mean"]
            left_errs[ap][name] = s["ci_half"]

    # Lambda_x = 0.1 prior sandwich matrix traces (K_e <= K_eL <= K0).
    matrix_names = ["K_e", "K_eL", "K0"]
    right_means = {ap: {mn: 0.0 for mn in matrix_names} for ap in APERTURES}
    right_errs = {ap: {mn: 0.0 for mn in matrix_names} for ap in APERTURES}
    for ap in APERTURES:
        recs = list(records_by_ap[ap].values())
        for mn in matrix_names:
            vals = []
            for rec in recs:
                eig = rec["priors"]["lambda01"][f"{mn}_eigenvalues"]
                vals.append(float(np.sum(np.asarray(eig, dtype=float))))
            s = _stats(vals)
            right_means[ap][mn] = s["mean"]
            right_errs[ap][mn] = s["ci_half"]

    colors = {"full": "#1f4e79", "limited": "#b03a2e"}
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4))

    ax = axes[0]
    x = np.arange(len(left_metrics))
    w = 0.34
    for j, ap in enumerate(APERTURES):
        means = [left_means[ap][m] for m in left_metrics]
        errs = [left_errs[ap][m] for m in left_metrics]
        ax.bar(
            x + (j - 0.5) * w,
            means,
            width=w,
            yerr=errs,
            color=colors[ap],
            alpha=0.9,
            capsize=4,
            label=f"{ap} aperture",
            error_kw={"elinewidth": 1.0, "ecolor": "black"},
        )
    ax.set_xticks(x)
    ax.set_xticklabels(
        [r"pose info trace  tr($J_x$)",
         r"map retention trace  tr($K_e$)"]
    )
    ax.set_ylabel("mean across seeds 301-310")
    ax.set_title(
        "Full vs limited aperture: pose information and map retention\n"
        "bars = mean, whiskers = 95% Student-t CI (df=9)"
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, axis="y")

    ax = axes[1]
    x = np.arange(len(matrix_names))
    for j, ap in enumerate(APERTURES):
        means = [right_means[ap][m] for m in matrix_names]
        errs = [right_errs[ap][m] for m in matrix_names]
        ax.bar(
            x + (j - 0.5) * w,
            means,
            width=w,
            yerr=errs,
            color=colors[ap],
            alpha=0.9,
            capsize=4,
            label=f"{ap} aperture",
            error_kw={"elinewidth": 1.0, "ecolor": "black"},
        )
    ax.set_xticks(x)
    ax.set_xticklabels([r"$K_e$", r"$K_{eL}(0.1 I_3)$", r"$K_0$"])
    ax.set_ylabel("mean matrix trace across seeds 301-310")
    ax.set_title(
        "Pose-prior Loewner sandwich at Lambda_x = 0.1 I3\n"
        "K_e <= K_eL <= K0 in Loewner order (trace demonstration); "
        "whiskers = 95% CI"
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, axis="y")

    fig.tight_layout()
    fig.savefig("figures/e3_physical_limited_prior.png", dpi=160)
    plt.close(fig)


def _write_summary_csv(records_by_ap) -> None:
    fieldnames = [
        "seed",
        "full_pose_info_trace",
        "limited_pose_info_trace",
        "full_map_retention_trace",
        "limited_map_retention_trace",
        "full_theorem3_max_res",
        "limited_theorem3_max_res",
        "full_prior01_ok",
        "limited_prior01_ok",
        "full_mc_var_rel",
        "limited_mc_var_rel",
        "full_mc_bias_rel",
        "limited_mc_bias_rel",
    ]
    with open("results/e3_physical_limited_summary.csv", "w",
              newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for seed in SEEDS:
            f = records_by_ap["full"][seed]
            l = records_by_ap["limited"][seed]
            f_tr, _ = _map_retention(f)
            l_tr, _ = _map_retention(l)
            fmc = _mc_block(f)
            lmc = _mc_block(l)
            writer.writerow({
                "seed": int(seed),
                "full_pose_info_trace": _fmt(f["Jx"]["pose_info_trace"]),
                "limited_pose_info_trace": _fmt(l["Jx"]["pose_info_trace"]),
                "full_map_retention_trace": _fmt(f_tr),
                "limited_map_retention_trace": _fmt(l_tr),
                "full_theorem3_max_res": _fmt(_theorem3_max_res(f)),
                "limited_theorem3_max_res": _fmt(_theorem3_max_res(l)),
                "full_prior01_ok": bool(
                    f["priors"]["lambda01"]["loewner_sandwich_ok"]
                ),
                "limited_prior01_ok": bool(
                    l["priors"]["lambda01"]["loewner_sandwich_ok"]
                ),
                "full_mc_var_rel": (
                    _fmt(fmc["variance_rel_error"]) if fmc else ""
                ),
                "limited_mc_var_rel": (
                    _fmt(lmc["variance_rel_error"]) if lmc else ""
                ),
                "full_mc_bias_rel": (
                    _fmt(fmc["bias_rel_error"]) if fmc else ""
                ),
                "limited_mc_bias_rel": (
                    _fmt(lmc["bias_rel_error"]) if lmc else ""
                ),
            })


def main() -> int:
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    models = {
        ap: physics.Model(physics.Config(N=N_GRID, aperture=ap))
        for ap in APERTURES
    }
    records = []
    for seed in SEEDS:
        for ap in APERTURES:
            rec = _run_one(models[ap], seed, ap)
            records.append(rec)
            print(
                f"seed {seed} aperture={ap:8s} "
                f"pose_tr={rec['Jx']['pose_info_trace']:.6g} "
                f"map_tr={sum(rec['theorem3']['raw_K_e_eigenvalues']):.6g} "
                f"thm3_max_res={rec['theorem3']['max_residual_abs']:.3g} "
                f"thm7_rank={rec['theorem7']['full_rank_B_v']}"
            )

    if len(records) != len(SEEDS) * len(APERTURES):
        raise RuntimeError("not all expected records were collected")
    records_by_ap = _by_aperture(records)

    always_metrics = [
        "pose_info_trace",
        "pose_info_min_eig",
        "map_retention_trace",
        "map_retention_min_eig",
        "theorem3_max_residual",
    ]
    available_metrics = [
        "theorem7_variance_rel_error",
        "theorem7_bias_rel_error",
        "sample_cov_vs_Jx_inv_fro_rel",
    ]
    summary = {
        "per_aperture_stats": {},
        "paired_full_minus_limited": {},
    }
    for name in always_metrics + available_metrics:
        summary["per_aperture_stats"][name] = _per_aperture_metric_stats(
            records_by_ap, name
        )
        summary["paired_full_minus_limited"][name] = _paired_metric_stats(
            records_by_ap, name
        )

    pose_pairs = []
    map_pairs = []
    limited_reduces_pose_per_seed = {}
    limited_reduces_map_per_seed = {}
    for seed in SEEDS:
        f = records_by_ap["full"][seed]
        l = records_by_ap["limited"][seed]
        ft, _ = _map_retention(f)
        lt, _ = _map_retention(l)
        fp = float(f["Jx"]["pose_info_trace"])
        lp = float(l["Jx"]["pose_info_trace"])
        pose_pairs.append(fp - lp)
        map_pairs.append(ft - lt)
        limited_reduces_pose_per_seed[int(seed)] = bool(fp > lp)
        limited_reduces_map_per_seed[int(seed)] = bool(ft > lt)

    pose_reduce_stats = _pair_stats(pose_pairs)
    map_reduce_stats = _pair_stats(map_pairs)

    counts = {}
    for ap in APERTURES:
        recs = list(records_by_ap[ap].values())
        counts[ap] = {
            "n_records": len(recs),
            "loewner_sandwich_ok": {
                "lambda0": int(
                    sum(r["priors"]["lambda0"]["loewner_sandwich_ok"]
                        for r in recs)
                ),
                "lambda01": int(
                    sum(r["priors"]["lambda01"]["loewner_sandwich_ok"]
                        for r in recs)
                ),
                "lambda1": int(
                    sum(r["priors"]["lambda1"]["loewner_sandwich_ok"]
                        for r in recs)
                ),
            },
            "theorem7_available": int(
                sum(1 for r in recs if r["theorem7"]["full_rank_B_v"])
            ),
            "variance_within_3se": int(
                sum(
                    1
                    for r in recs
                    if _mc_block(r) is not None
                    and _mc_block(r)["variance_within_3se"]
                )
            ),
            "bias_within_3se": int(
                sum(
                    1
                    for r in recs
                    if _mc_block(r) is not None
                    and _mc_block(r)["bias_within_3se"]
                )
            ),
            "both_within_3se": int(
                sum(
                    1
                    for r in recs
                    if _mc_block(r) is not None
                    and _mc_block(r)["variance_within_3se"]
                    and _mc_block(r)["bias_within_3se"]
                )
            ),
        }

    all_thm3 = all(
        float(r["theorem3"]["max_residual_abs"]) <= 1e-8 for r in records
    )
    all_prior = all(
        r["priors"]["lambda01"]["loewner_sandwich_ok"]
        and r["priors"]["lambda1"]["loewner_sandwich_ok"]
        for r in records
    )
    mc_records = [r for r in records if _mc_block(r) is not None]
    all_mc = bool(
        mc_records
        and all(_mc_block(r)["variance_within_3se"] for r in mc_records)
        and all(_mc_block(r)["bias_within_3se"] for r in mc_records)
    ) if mc_records else True

    pass_flags = {
        "all_theorem3_residuals_pass": bool(all_thm3),
        "all_prior_sandwiches_pass": bool(all_prior),
        "all_mc_within_3se": bool(all_mc),
        "limited_reduces_pose_info": bool(pose_reduce_stats["ci_low"] > 0.0),
        "limited_reduces_map_retention": bool(
            map_reduce_stats["ci_low"] > 0.0
        ),
        "hard_pass_criteria_note": (
            "only theorem3 residuals (<=1e-8), prior Loewner sandwiches, and "
            "Theorem-7 MC within 3 SE are hard pass criteria; the "
            "limited-reduces booleans are reported observations."
        ),
    }
    hard_pass = bool(
        pass_flags["all_theorem3_residuals_pass"]
        and pass_flags["all_prior_sandwiches_pass"]
        and pass_flags["all_mc_within_3se"]
    )

    summary["counts_by_aperture"] = counts
    summary["limited_reduces"] = {
        "pose_info_trace": {
            "paired_full_minus_limited_stats": pose_reduce_stats,
            "per_seed": limited_reduces_pose_per_seed,
        },
        "map_retention_trace": {
            "paired_full_minus_limited_stats": map_reduce_stats,
            "per_seed": limited_reduces_map_per_seed,
        },
    }
    summary["pass_flags"] = pass_flags

    interpreter = sys.executable
    script = os.path.join("scripts", os.path.basename(__file__))
    doc = {
        "title": (
            "E3 physical extension: FULL vs LIMITED aperture, with and "
            "without pose priors (N=8 smoke, seeds 301-310)"
        ),
        "date": "2026-09-05",
        "package": "a2val",
        "settings": {
            "physics_module": {
                "name": physics.__name__,
                "file": physics.__file__,
            },
            "config": f"physics.Config(N={N_GRID}, aperture in "
                      "('full','limited'))",
            "seeds": SEEDS,
            "apertures": APERTURES,
            "theorem7_mc_draws": int(M_DRAWS),
            "ci": "95% Student-t, df=9 (scipy.stats.t.ppf(0.975,9))",
            "reference": (
                "rng=default_rng(seed); alpha=0.2+0.05*N clipped to "
                "[-0.5,0.5]; x=0.02*N clipped to [-0.05,0.05]"
            ),
            "declared_C1": (
                "first 4 left singular vectors of T=hstack([A_r,B_r]); "
                "no oracle rank"
            ),
            "realification": (
                "physics.realify_jacobian: sqrt(2)*[Re(J); Im(J)] "
                "(unit-noise whitened)"
            ),
            "theorem7_gate": (
                "evaluated only when min_eig(J_x)=min_eig(R^T R) > 1e-10; "
                "otherwise recorded full_rank_B_v=False without exception"
            ),
            "bias_direction_convention": (
                "D_r = first left singular vector of R when "
                "S[0] > max(1e-14, eps*||R||_F), else random unit fallback "
                "from the per-seed rng"
            ),
            "read_only_boundary": (
                "nothing under research/delegated/a2_physics/ is modified"
            ),
            "pass_thresholds": {
                "theorem3_max_residual": 1e-8,
                "prior_loewner_relative": 1e-8,
                "theorem7_mc_3se": True,
            },
        },
        "records": records,
        "summary": summary,
        "verdict": {
            "run_completed_ok": True,
            "hard_pass": hard_pass,
            **pass_flags,
        },
    }

    with open("results/e3_physical_limited_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")

    _write_summary_csv(records_by_ap)
    _make_figure(records_by_ap)

    command = [f"cd {ROOT}", f"{interpreter} {script}"]
    with open("results/e3_physical_limited_commands.txt", "w") as fh:
        fh.write("Exact command run (E3 physical limited aperture)\n")
        fh.write("\n".join(command))
        fh.write("\n")

    print("\nSummary (mean across 10 seeds; 95% CI paired full-limited):")
    pt = summary["paired_full_minus_limited"]["pose_info_trace"]
    mt = summary["paired_full_minus_limited"]["map_retention_trace"]
    print(f"  pose_info_trace  diff: {_fmt(pt['mean'])} "
          f"[{_fmt(pt['ci_low'])}, {_fmt(pt['ci_high'])}]")
    print(f"  map_retention_trace diff: {_fmt(mt['mean'])} "
          f"[{_fmt(mt['ci_low'])}, {_fmt(mt['ci_high'])}]")
    print("  pass flags:", pass_flags)
    print("  hard_pass:", hard_pass)
    print("Saved: results/e3_physical_limited_results.json, "
          "results/e3_physical_limited_summary.csv, "
          "figures/e3_physical_limited_prior.png, "
          "results/e3_physical_limited_commands.txt")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:  # noqa: BLE001 - surface the failure for fixing
        traceback.print_exc()
        raise SystemExit(1)
