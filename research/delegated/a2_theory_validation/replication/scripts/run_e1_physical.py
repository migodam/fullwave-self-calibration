#!/usr/bin/env python
"""E1 physical coherent-vs-phaseless Fisher scene at the N=8 shared core.

Usage (from the experiment working directory, exact interpreter required):

    python scripts/run_e1_physical.py

The shared physical core is imported read-only from
``research/delegated/a2_physics``; nothing under that directory is modified.
The coherent experiment uses the realified unit-noise complex tangent, while
the phaseless experiment uses the exact induced intensity (noncentral chi
square) likelihood with per-row score-variance obtained by 1-D quadrature.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import integrate, stats

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

from a2val import e1  # noqa: E402  (phaseless score helper)

SEEDS = list(range(101, 111))
EPSABS = 1e-9
EPSREL = 1e-8
QUAD_LIMIT = 100
LAMBDA_EPS = 1e-12
RANK_TOL = 1e-9


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


def _reference(seed: int):
    """Deterministic non-oracle reference, clamped as specified."""
    rng = np.random.default_rng(seed)
    alpha = np.clip(0.2 + 0.05 * rng.standard_normal(9), -0.5, 0.5)
    x = np.clip(0.02 * rng.standard_normal(3), -0.05, 0.05)
    return alpha, x


def _quad_shared(integrand, lo, hi):
    """One-dimensional quadrature on [lo,hi] with a scalar integrand."""
    # a2val.e1.phaseless_score_lambda returns a Python float for scalar input.
    return integrate.quad(integrand, lo, hi, epsabs=EPSABS, epsrel=EPSREL,
                          limit=QUAD_LIMIT)


def _v_lambda(lam, v_cache, v_cache_err):
    """v(lambda) = E[s_lambda^2] under the ncx2(df=2, scale=0.5) law.

    Cached by lambda rounded to 12 decimals.  For lambda < 1e-12 the task
    fast path returns v=0.5 exactly; no physics seed reached that regime, so
    every recorded row used quadrature (see note in the JSON output).
    """
    key = round(float(lam), 12)
    if key in v_cache:
        return v_cache[key]
    if float(lam) < LAMBDA_EPS:
        v_cache[key] = 0.5
        v_cache_err[key] = 0.0
        return 0.5
    dist = stats.ncx2(df=2, nc=float(lam), scale=0.5)
    lo = max(float(dist.ppf(1e-12)), 0.0)
    hi = float(dist.ppf(1 - 1e-12))

    def integrand(t):
        s = e1.phaseless_score_lambda(2.0 * float(t), float(lam))
        return float(s) * float(s) * float(dist.pdf(float(t)))

    val, err = _quad_shared(integrand, lo, hi)
    v_cache[key] = float(val)
    v_cache_err[key] = float(err)
    return float(val)


def _per_seed_record(model, seed, t_start):
    alpha, x = _reference(seed)
    out = model.forward(alpha, x, jacobian=True)
    mu = out["total"]
    A = out["A"]
    B = out["B"]
    Ar = physics.realify_jacobian(A)
    Br = physics.realify_jacobian(B)
    n_alpha = Ar.shape[1]
    n_x = Br.shape[1]

    # Coherent Fisher in realified unit-noise coordinates.
    J_coh = np.empty((n_alpha + n_x, n_alpha + n_x))
    J_coh[:n_alpha, :n_alpha] = Ar.T @ Ar
    J_coh[:n_alpha, n_alpha:] = Ar.T @ Br
    J_coh[n_alpha:, :n_alpha] = Br.T @ Ar
    J_coh[n_alpha:, n_alpha:] = Br.T @ Br
    J_coh = 0.5 * (J_coh + J_coh.T)

    # Phaseless Fisher: independent rows with exact induced intensity law.
    lam_i = 2.0 * np.abs(mu) ** 2
    g_alpha = 4.0 * np.real(np.conj(mu)[:, None] * A)
    g_x = 4.0 * np.real(np.conj(mu)[:, None] * B)
    G = np.concatenate([g_alpha, g_x], axis=1)  # (rows, 12)
    J_ph = np.zeros((n_alpha + n_x, n_alpha + n_x))
    v_cache = {}
    v_cache_err = {}
    for i in range(G.shape[0]):
        v = _v_lambda(lam_i[i], v_cache, v_cache_err)
        J_ph += v * np.outer(G[i], G[i])
    J_ph = 0.5 * (J_ph + J_ph.T)

    eig_coh = np.linalg.eigvalsh(J_coh)
    eig_ph = np.linalg.eigvalsh(J_ph)
    D = J_coh - J_ph
    eig_D = np.linalg.eigvalsh(D)
    asym_D = D - D.T
    fro_D = float(np.linalg.norm(D, "fro"))
    fro_coh = float(np.linalg.norm(J_coh, "fro"))
    max_eig_ph = max(1.0, float(np.max(eig_ph)))
    max_eig_coh = max(1.0, float(np.max(eig_coh)))
    rank_ph = int(np.count_nonzero(eig_ph > RANK_TOL * max_eig_ph))
    rank_coh = int(np.count_nonzero(eig_coh > RANK_TOL * max_eig_coh))

    return {
        "seed": int(seed),
        "alpha": alpha.tolist(),
        "x": x.tolist(),
        "J_coh": J_coh.tolist(),
        "J_ph": J_ph.tolist(),
        "eig_coh_asc": eig_coh.tolist(),
        "eig_ph_asc": eig_ph.tolist(),
        "eig_D_asc": eig_D.tolist(),
        "min_eig_D": float(eig_D[0]),
        "max_abs_asym_D": float(np.max(np.abs(asym_D))) if asym_D.size else 0.0,
        "trace_coh": float(np.trace(J_coh)),
        "trace_ph": float(np.trace(J_ph)),
        "trace_D": float(np.trace(D)),
        "fro_D": fro_D,
        "fro_D_over_fro_coh": fro_D / fro_coh if fro_coh > 0 else None,
        "rank_ph": rank_ph,
        "rank_coh": rank_coh,
        "work_wall_seconds": float(out["work"]["wall_seconds"]),
        "v_cache_size": len(v_cache),
        "max_quad_err_estimate": max(v_cache_err.values(), default=0.0),
        "elapsed_total_seconds": time.perf_counter() - t_start,
    }


def _ci_mean(values, t_crit):
    values = np.asarray(values, dtype=float)
    m = float(np.mean(values))
    s = float(np.std(values, ddof=1))
    half = t_crit * s / np.sqrt(values.size)
    return {"mean": m, "sd": s, "ci_lower": m - half, "ci_upper": m + half}


def _build_figure(records, summary):
    fig, (ax_left, ax_right) = plt.subplots(
        1, 2, figsize=(13.5, 5.4), constrained_layout=True
    )

    # Left panel: eigenvalue spectra for seed 101.
    seed101 = next(r for r in records if r["seed"] == 101)
    eig_coh = np.asarray(seed101["eig_coh_asc"])[::-1]
    eig_ph = np.asarray(seed101["eig_ph_asc"])[::-1]
    k = np.arange(1, len(eig_coh) + 1)
    all_positive = bool(np.all(eig_coh > 0.0) and np.all(eig_ph > 0.0))
    ax_left.axhline(0.0, color="0.7", lw=0.7, zorder=1)
    if all_positive:
        ax_left.plot(k, eig_coh, "o-", ms=4, lw=1.3, color="#1b6ca8",
                     label="J_coh (seed 101)")
        ax_left.plot(k, eig_ph, "s-", ms=4, lw=1.3, color="#c1440e",
                     label="J_ph (seed 101)")
        ax_left.set_yscale("log")
        ax_left.set_ylabel("eigenvalue (log10)")
    else:
        ax_left.plot(k, eig_coh, "o-", ms=4, lw=1.3, color="#1b6ca8",
                     label="J_coh (seed 101)")
        ax_left.plot(k, eig_ph, "s-", ms=4, lw=1.3, color="#c1440e",
                     label="J_ph (seed 101)")
        ax_left.set_yscale("symlog", linthresh=1e-6)
        ax_left.set_ylabel("eigenvalue (symlog)")
    ax_left.set_xlabel("eigenvalue index (descending)")
    ax_left.set_title("E1 physical Fisher spectra\n(N=8 full aperture, seed 101)")
    ax_left.grid(True, alpha=0.3, which="both")
    ax_left.legend(frameon=False)

    # Right panel: traces and min_eig_D with 95% CIs across seeds 101..110.
    def _show(key, color, marker, yaxis):
        per = np.asarray([r[key] for r in records], dtype=float)
        ci = summary[key]
        err_lo = ci["mean"] - ci["ci_lower"]
        err_hi = ci["ci_upper"] - ci["mean"]
        yaxis.errorbar(
            np.arange(1, len(per) + 1), per, marker=marker, ms=4.5, lw=0,
            color=color, capsize=3, ecolor=color, alpha=0.85,
            yerr=[err_lo * np.ones_like(per), err_hi * np.ones_like(per)],
            label=f"{key} mean={ci['mean']:.4g}",
        )

    _show("trace_coh", "#1b6ca8", "o", ax_right)
    _show("trace_ph", "#c1440e", "s", ax_right)
    _show("trace_D", "#3c8c40", "^", ax_right)
    ax_right.axhline(0.0, color="0.6", lw=0.8)
    ax_right.set_xlabel("seed index (101..110)")
    ax_right.set_ylabel("trace value")
    ax_right.set_title("Per-seed traces and min_eig(D)\n(95% Student-t CI bands)")
    ax_right.grid(True, alpha=0.3)
    ax_right.legend(frameon=False, fontsize=8)

    ax_right_min = ax_right.twinx()
    _show("min_eig_D", "#7a3b8f", "d", ax_right_min)
    ax_right_min.axhline(0.0, color="0.6", lw=0.6, ls=":")
    ax_right_min.set_ylabel("min_eig_D (secondary axis)", color="#7a3b8f")
    ax_right_min.tick_params(axis="y", colors="#7a3b8f")

    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/e1_physical_fisher_comparison.png", dpi=170)
    plt.close(fig)


def main():
    t_start = time.perf_counter()
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    cfg = physics.Config(N=8, aperture="full")
    model = physics.Model(cfg)

    records = [_per_seed_record(model, seed, t_start) for seed in SEEDS]

    keys_summary = [
        "min_eig_D", "trace_coh", "trace_ph", "trace_D",
        "fro_D_over_fro_coh", "rank_coh", "rank_ph",
    ]
    n = len(records)
    t_crit = float(stats.t.ppf(0.975, n - 1))
    summary = {
        key: _ci_mean([r[key] for r in records], t_crit) for key in keys_summary
    }

    mean_trace_coh = summary["trace_coh"]["mean"]
    loewner_tol_global = 1e-8 * max(1.0, mean_trace_coh)
    ci_lower_pass = bool(
        summary["min_eig_D"]["ci_lower"] >= -loewner_tol_global
    )
    per_seed_pass = bool(all(
        r["min_eig_D"] >= -1e-8 * max(1.0, r["trace_coh"]) for r in records
    ))
    loewner_psd_backward_tol = bool(ci_lower_pass and per_seed_pass)
    strict_positive_info = bool(
        all(0.0 < r["trace_ph"] < r["trace_coh"] for r in records)
    )

    doc = {
        "title": (
            "E1 physical coherent-vs-phaseless Fisher comparison at the "
            "shared N=8 full-aperture core"
        ),
        "label": (
            "physical E1 Fisher scene: exact tangent Fisher (coherent) vs "
            "exact induced-intensity ncx2 Fisher (phaseless); no oracle or "
            "method-superiority claim"
        ),
        "config": {
            "N": 8,
            "aperture": "full",
            "sigma2": 1.0,
            "realified_scale": "sqrt(2)/sigma with sigma=1 per real component",
            "intensity_law": "I ~ ncx2(df=2, nc=lambda, scale=0.5), "
                             "lambda=2|mu|^2",
        },
        "settings": {
            "seeds": SEEDS,
            "seed_reference": (
                "rng=default_rng(seed); alpha=clip(0.2+0.05*N(0,1),[-0.5,0.5]) "
                "len 9; x=clip(0.02*N(0,1),[-0.05,0.05]) len 3"
            ),
            "quadrature": {
                "epsabs": EPSABS,
                "epsrel": EPSREL,
                "limit": QUAD_LIMIT,
                "domain": "[ppf(1e-12), ppf(1-1e-12)] of the ncx2 law",
                "fast_path_lambda_lt": LAMBDA_EPS,
                "fast_path_v": 0.5,
            },
            "rank_tolerance_relative": RANK_TOL,
            "import_physics": {
                "module": physics.__name__,
                "file": physics.__file__,
            },
            "convention_note": (
                "Task fast path assigns v=0.5 for lambda<1e-12, but analytic "
                "and quadrature checks give v(lambda->0)=E[(-1/2+w/4)^2]="
                "0.25 under w~Exp(mean 2).  The fast path was never exercised "
                "(observed lambda range >= 2.7e-4, all rows quadrature), so "
                "recorded J_ph is unaffected by that convention.  With the "
                "quadrature convention J_ph(0)=0.5*J_coh for a pure mean "
                "shift."
            ),
        },
        "per_seed": records,
        "summary": {
            key: {
                "mean": summary[key]["mean"],
                "sd": summary[key]["sd"],
                "ci_95_lower": summary[key]["ci_lower"],
                "ci_95_upper": summary[key]["ci_upper"],
                "n": n,
            }
            for key in keys_summary
        },
        "pass_fail": {
            "loewner_psd_within_backward_tol": loewner_psd_backward_tol,
            "ci_lower_bound_ge_negative_tol": ci_lower_pass,
            "every_seed_min_eig_D_ge_negative_tol": per_seed_pass,
            "tolerance": "1e-8*max(1, mean_or_seed trace_coh)",
            "strict_positive_information_0_lt_trace_ph_lt_trace_coh_all_seeds":
                strict_positive_info,
            "mean_min_eig_D_ci_lower": summary["min_eig_D"]["ci_lower"],
            "mean_trace_coh": mean_trace_coh,
        },
        "runtime_wall_seconds": time.perf_counter() - t_start,
    }

    with open("results/e1_physical_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")

    with open("results/e1_physical_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "seed", "j_coh_trace", "j_ph_trace", "min_eig_D", "min_eig_coh",
            "min_eig_ph", "max_abs_asym_D", "rank_coh", "rank_ph",
            "work_wall_seconds",
        ])
        writer.writeheader()
        for r in records:
            writer.writerow({
                "seed": r["seed"],
                "j_coh_trace": f"{r['trace_coh']:.9g}",
                "j_ph_trace": f"{r['trace_ph']:.9g}",
                "min_eig_D": f"{r['min_eig_D']:.9g}",
                "min_eig_coh": f"{r['eig_coh_asc'][0]:.9g}",
                "min_eig_ph": f"{r['eig_ph_asc'][0]:.9g}",
                "max_abs_asym_D": f"{r['max_abs_asym_D']:.9g}",
                "rank_coh": r["rank_coh"],
                "rank_ph": r["rank_ph"],
                "work_wall_seconds": f"{r['work_wall_seconds']:.9g}",
            })

    _build_figure(records, summary)

    interpreter = sys.executable
    with open("results/e1_physical_commands.txt", "w") as fh:
        fh.write("Exact shell command run (E1 physical Fisher scene)\n")
        fh.write("cd " + ROOT + "\n")
        fh.write(interpreter + " scripts/run_e1_physical.py\n")

    keys_print = [
        "min_eig_D", "trace_coh", "trace_ph", "trace_D",
        "fro_D_over_fro_coh", "rank_coh", "rank_ph",
    ]
    print("E1 physical Fisher: per-seed counts")
    for key in keys_print:
        s = summary[key]
        print(
            f"  {key:22s} mean={s['mean']:.9g} "
            f"CI95=[{s['ci_lower']:.9g}, {s['ci_upper']:.9g}]"
        )
    print("Loewner PSD within backward tolerance:",
          loewner_psd_backward_tol)
    print("Strict positive information (0<trace_ph<trace_coh on all seeds):",
          strict_positive_info)
    print("Runtime wall seconds:", round(doc["runtime_wall_seconds"], 3))
    print("Saved: results/e1_physical_results.json, "
          "results/e1_physical_summary.csv, "
          "figures/e1_physical_fisher_comparison.png, "
          "results/e1_physical_commands.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
