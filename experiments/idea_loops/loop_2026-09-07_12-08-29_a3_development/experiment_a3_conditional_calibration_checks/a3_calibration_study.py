#!/usr/bin/env python3
"""A3 extended finite-noise calibration and forward-model verification study.

Self-contained companion to the frozen deterministic A3 checks.  This script
only READS checks.json (frozen baseline audit); it never modifies
regenerate_checks.py, checks.json, or SUMMARY.md.

Parts
-----
0.  frozen baseline audit of checks.json counts
1.  B3 finite-noise branch-bound calibration (Gaussian + heavy-tailed)
2.  sharpened B3 reused-validation (3E) control
3.  B2 finite-noise ratio sufficiency / projector statistics
4.  forward-model B4 dual-residual verification on two toy linearized SLAM
    systems (2D pose graph and 3D bundle adjustment)
5.  near-singular forward-model B4 control (small residual alone is not an
    output-error certificate; the certified gamma bound still holds)
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import scipy
from scipy.stats import chi2 as scipy_chi2

MPL_OK = True
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover - environment dependent
    MPL_OK = False

HERE = os.path.dirname(os.path.abspath(__file__))
CHECKS_PATH = os.path.join(HERE, "checks.json")


# ---------------------------------------------------------------------------
# generic helpers
# ---------------------------------------------------------------------------
def _norm(x):
    return float(np.linalg.norm(x))


def _clean_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _clean_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def _wilson_95(p, n, z=1.96):
    """Wilson 95% interval for a binomial proportion (clipped)."""
    p = float(p)
    n = float(n)
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denom
    half = (z / denom) * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n))
    lo = max(0.0, centre - half)
    hi = min(1.0, centre + half)
    return float(lo), float(hi)


# ---------------------------------------------------------------------------
# B2 helpers copied from the frozen deterministic baseline
# ---------------------------------------------------------------------------
def cross_ratio_matrix(H):
    """Return (C, flag); r,t>1 in python indices (>=1)."""
    R, T = H.shape
    out = np.empty((R - 1, T - 1), dtype=complex)
    for r in range(1, R):
        for t in range(1, T):
            den1 = H[r, 0]
            den2 = H[0, t]
            if abs(den1) < 1e-9 or abs(den2) < 1e-9:
                return None, "near_zero_denominator"
            out[r - 1, t - 1] = (H[r, t] * H[0, 0]) / (den1 * den2)
    return out, None


def _random_gain_vector(length, rng):
    mag = 0.5 + rng.random_sample(length)
    ph = 2.0 * np.pi * rng.random_sample(length)
    return mag * np.exp(1j * ph)


def _log_cross_ratio_operator(R, T):
    """Complex K x (R*T) linearized log cross-ratio operator."""
    K = (R - 1) * (T - 1)
    N = R * T
    L = np.zeros((K, N), dtype=complex)
    k = 0
    for r in range(1, R):
        for t in range(1, T):
            L[k, r * T + t] += 1.0
            L[k, 0] += 1.0
            L[k, r * T + 0] -= 1.0
            L[k, 0 * T + t] -= 1.0
            k += 1
    return L


def _gain_subspace(R, T):
    """Complex N x (R+T-1) basis of {E: E_rt = u_r + v_t}."""
    N = R * T
    D = R + T - 1
    G = np.zeros((N, D), dtype=complex)
    d = 0
    for r in range(R):
        E = np.zeros((R, T), dtype=complex)
        E[r, :] = 1.0
        G[:, d] = E.ravel(order="C")
        d += 1
    for t in range(1, T):
        E = np.zeros((R, T), dtype=complex)
        E[:, t] = 1.0
        G[:, d] = E.ravel(order="C")
        d += 1
    return G


def _realify_matrix(Lc):
    """2K x 2N real matrix realizing a complex K x N linear map."""
    K, N = Lc.shape
    LR = np.zeros((2 * K, 2 * N))
    LR[:K, :N] = Lc.real
    LR[:K, N:] = -Lc.imag
    LR[K:, :N] = Lc.imag
    LR[K:, N:] = Lc.real
    return LR


def _realify_complex_subspace(Gc):
    """Realification (2N x 2D) of a complex-linear subspace with basis Gc."""
    N, D = Gc.shape
    cols = []
    for d in range(D):
        g = Gc[:, d]
        cols.append(np.concatenate([g.real, g.imag]))
        cols.append(np.concatenate([-g.imag, g.real]))
    B = np.column_stack(cols)
    Q, _ = np.linalg.qr(B)
    sv_all = np.linalg.svd(B, compute_uv=False)
    thresh = 1e-10 * max(1.0, float(sv_all[0]))
    rank = int(np.sum(sv_all > thresh))
    return Q[:, :rank], rank


def _orth_projector(A):
    """P = A (A^T A)^{-1} A^T for full-column-rank real A."""
    return A @ np.linalg.solve(A.T @ A, A.T)


def _realify_vec(v):
    return np.concatenate([np.asarray(v).real.ravel(),
                           np.asarray(v).imag.ravel()])


def _lift_right(LR, v):
    """Right-inverse lift of a realified quotient vector into 2N-space."""
    return LR.T @ np.linalg.solve(LR @ LR.T, v)


# ---------------------------------------------------------------------------
# PART 0: frozen baseline audit
# ---------------------------------------------------------------------------
def part0_frozen_audit():
    with open(CHECKS_PATH, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    counts = data["counts"]
    expected = {"pass": 16, "expected_fail": 4, "fail": 0, "total": 20}
    ok = bool(counts == expected)
    if not ok:
        raise AssertionError(
            f"frozen baseline counts mismatch: {counts} != {expected}")
    return {
        "frozen_baseline_ok": ok,
        "frozen_counts": dict(counts),
        "checks_path": CHECKS_PATH,
        "cases": len(data["cases"]),
    }


# ---------------------------------------------------------------------------
# PART 1: B3 finite-noise branch-bound calibration
# ---------------------------------------------------------------------------
def _branch_bound(K, mus, i_star, beta, sigma):
    """Plug-in Gaussian branch bound, per-coordinate noise variance sigma^2."""
    d_ij = np.linalg.norm(mus - mus[i_star], axis=1)
    terms = []
    for j in range(K):
        if j == i_star:
            continue
        arg = max(d_ij[j] - 2.0 * beta, 0.0)
        terms.append(np.exp(-(arg ** 2) / (8.0 * sigma * sigma)))
    return float(min(1.0, sum(terms)))


def _b3_geometries():
    m = 4
    geo = {}
    geo["geo_d4"] = {
        "name": "geo_d4", "d": 4.0, "K": 3, "m": m, "i_star": 0,
        "competitor": 1,
        "mus": np.array([[0, 0, 0, 0], [4, 0, 0, 0], [0, 4, 0, 0]],
                        dtype=float),
    }
    geo["geo_d6"] = {
        "name": "geo_d6", "d": 6.0, "K": 3, "m": m, "i_star": 0,
        "competitor": 1,
        "mus": np.array([[0, 0, 0, 0], [6, 0, 0, 0], [0, 6, 0, 0]],
                        dtype=float),
    }
    geo["geo_d8"] = {
        "name": "geo_d8", "d": 8.0, "K": 3, "m": m, "i_star": 0,
        "competitor": 1,
        "mus": np.array([[0, 0, 0, 0], [8, 0, 0, 0], [0, 8, 0, 0]],
                        dtype=float),
    }
    return geo


def _b3_configs():
    """Deterministic list of (noise_model, geometry_name, beta, sigma)."""
    configs = []
    # Gaussian
    for beta in (0.25, 0.5, 1.0):
        for sigma in (0.5, 1.0, 2.0):
            configs.append(("gaussian", "geo_d4", beta, sigma))
    for sigma in (0.5, 1.0, 2.0):
        configs.append(("gaussian", "geo_d6", 1.0, sigma))
        configs.append(("gaussian", "geo_d8", 1.0, sigma))
    # heavy-tailed Student-t(3), scaled to per-coordinate variance sigma^2
    for beta in (0.25, 0.5):
        for sigma in (1.0, 2.0):
            configs.append(("heavy_t3", "geo_d4", beta, sigma))
    for sigma in (1.0, 2.0):
        configs.append(("heavy_t3", "geo_d6", 1.0, sigma))
        configs.append(("heavy_t3", "geo_d8", 1.0, sigma))
    return configs


def _run_b3_seed(seed, geom, beta, sigma, noise_model, reps_per_seed):
    rng = np.random.RandomState(seed)
    m = geom["m"]
    mus = geom["mus"]
    i_star = geom["i_star"]
    j = geom["competitor"]
    e = mus[j] - mus[i_star]
    e = beta * e / float(np.linalg.norm(e))
    if noise_model == "gaussian":
        noise = sigma * rng.standard_normal((m, reps_per_seed))
    else:
        noise = (sigma / np.sqrt(3.0)) * rng.standard_t(3, size=(m, reps_per_seed))
    y = (mus[i_star] + e)[:, None] + noise          # (m, reps)
    y = y.T                                         # (reps, m)
    dists = np.linalg.norm(y[:, None, :] - mus[None, :, :], axis=2)
    sel = np.argmin(dists, axis=1)
    return int(np.sum(sel != i_star))


def part1_b3_calibration():
    seeds = list(range(9101, 9151))
    reps_per_seed = 2000
    n_seed_boot = 2000
    geos = _b3_geometries()
    rows = []
    for noise_model, geo_name, beta, sigma in _b3_configs():
        geom = geos[geo_name]
        mus = geom["mus"]
        i_star = geom["i_star"]
        bound_s1 = _branch_bound(geom["K"], mus, i_star, beta, sigma=1.0)
        bound_s = _branch_bound(geom["K"], mus, i_star, beta, sigma=sigma)
        errors = [
            _run_b3_seed(s, geom, beta, sigma, noise_model, reps_per_seed)
            for s in seeds
        ]
        total_err = int(sum(errors))
        total_draws = int(len(seeds) * reps_per_seed)
        p_hat = total_err / total_draws
        rates = [err / reps_per_seed for err in errors]
        rates_arr = np.asarray(rates, dtype=float)
        wil_lo, wil_hi = _wilson_95(p_hat, total_draws)
        boot_rng = np.random.RandomState(
            int(1000 + int(sigma * 100) + beta * 10000))
        boot = np.mean(boot_rng.choice(rates_arr, size=(n_seed_boot, len(rates_arr)),
                                       replace=True), axis=1)
        boot_lo, boot_hi = (float(x) for x in
                            np.percentile(boot, [2.5, 97.5]))
        rows.append({
            "noise_model": noise_model,
            "geometry": geo_name,
            "d": geom["d"],
            "K": geom["K"],
            "m": geom["m"],
            "beta": float(beta),
            "sigma": float(sigma),
            "i_star": int(i_star),
            "bound_sigma1": float(bound_s1),
            "bound_sigma": float(bound_s),
            "p_hat": float(p_hat),
            "wilson_lo": wil_lo,
            "wilson_hi": wil_hi,
            "boot_lo": boot_lo,
            "boot_hi": boot_hi,
            "per_seed_rates": rates,
            "covered_sigma1": bool(p_hat <= bound_s1 + 0.03),
            "covered_sigma": bool(p_hat <= bound_s + 0.03),
            "strong_violation_sigma": bool(boot_lo > bound_s),
            "strong_violation_sigma1": bool(boot_lo > bound_s1),
            "total_draws": total_draws,
        })
    gauss = [r for r in rows if r["noise_model"] == "gaussian"]
    heavy = [r for r in rows if r["noise_model"] == "heavy_t3"]
    summary = {
        "gaussian_configs_total": len(gauss),
        "gaussian_covered_sigma_count": sum(r["covered_sigma"] for r in gauss),
        "gaussian_covered_sigma1_count": sum(r["covered_sigma1"] for r in gauss),
        "gaussian_uncovered_sigma": [r["geometry"] + "/beta=" + str(r["beta"])
                                     + "/sigma=" + str(r["sigma"])
                                     for r in gauss if not r["covered_sigma"]],
        "gaussian_uncovered_sigma1": [r["geometry"] + "/beta=" + str(r["beta"])
                                      + "/sigma=" + str(r["sigma"])
                                      for r in gauss if not r["covered_sigma1"]],
        "heavy_configs_total": len(heavy),
        "heavy_strong_violation_vs_matched_sigma_count":
            sum(r["strong_violation_sigma"] for r in heavy),
        "heavy_strong_violation_vs_sigma1_count":
            sum(r["strong_violation_sigma1"] for r in heavy),
        "heavy_uncovered_sigma1": [r["geometry"] + "/beta=" + str(r["beta"])
                                   + "/sigma=" + str(r["sigma"])
                                   for r in heavy if not r["covered_sigma1"]],
        "heavy_uncovered_sigma": [r["geometry"] + "/beta=" + str(r["beta"])
                                  + "/sigma=" + str(r["sigma"])
                                  for r in heavy if not r["covered_sigma"]],
    }
    return {"configs": rows, "summary": summary}


# ---------------------------------------------------------------------------
# PART 2: sharpened B3 reused-validation (3E) control
# ---------------------------------------------------------------------------
def part2_sharpened_3e():
    m, beta = 4, 0.1
    seeds = list(range(9151, 9201))
    reps_per_seed = 2000
    n_lt1 = 0
    n_lt09 = 0
    n_lt05 = 0
    total = 0
    plugin_sum = 0.0
    bound_min = 1.0
    bound_max = 0.0
    per_seed_lt1 = []
    per_seed_avg = []
    for seed in seeds:
        rng = np.random.RandomState(seed)
        y = rng.standard_normal((m, reps_per_seed))           # m x reps
        d_refit = np.linalg.norm(y, axis=0)                   # reps
        arg = np.clip(d_refit - 2.0 * beta, 0.0, None)
        B = np.minimum(1.0, np.exp(-(arg ** 2) / 8.0))
        n_lt1 += int(np.sum(B < 1.0))
        n_lt09 += int(np.sum(B < 0.9))
        n_lt05 += int(np.sum(B < 0.5))
        plugin_sum += float(np.sum(B))
        bound_min = min(bound_min, float(np.min(B)))
        bound_max = max(bound_max, float(np.max(B)))
        total += reps_per_seed
        per_seed_lt1.append(float(np.mean(B < 1.0)))
        per_seed_avg.append(float(np.mean(B)))
    out = {
        "m": m,
        "beta": float(beta),
        "seeds": [int(min(seeds)), int(max(seeds))],
        "reps_per_seed": reps_per_seed,
        "total_reps": int(total),
        "empirical_error_by_construction": 1.0,
        "n_bound_lt_1": int(n_lt1),
        "n_bound_lt_0p9": int(n_lt09),
        "n_bound_lt_0p5": int(n_lt05),
        "fraction_bound_lt_1": float(n_lt1 / total),
        "avg_plugin_bound": float(plugin_sum / total),
        "min_plugin_bound": float(bound_min),
        "max_plugin_bound": float(bound_max),
        "per_seed_bound_lt_1_fraction": per_seed_lt1,
        "per_seed_avg_plugin": per_seed_avg,
        "violation_count": int(n_lt1),
    }
    return out


# ---------------------------------------------------------------------------
# PART 3: B2 finite-noise ratio sufficiency
# ---------------------------------------------------------------------------
def part3_b2_finite_noise():
    R, T = 5, 4
    K = (R - 1) * (T - 1)
    N = R * T
    Lc = _log_cross_ratio_operator(R, T)
    Gc = _gain_subspace(R, T)
    LR = _realify_matrix(Lc)
    GR_raw, rankG = _realify_complex_subspace(Gc)

    # frozen 2D covariance machinery
    rng = np.random.RandomState(4204)
    Qr, _ = np.linalg.qr(rng.standard_normal((2 * N, 2 * N)))
    eig = np.linspace(1.0, 5.0, 2 * N)
    Sigma = (Qr * eig) @ Qr.T
    S12 = (Qr * np.sqrt(eig)) @ Qr.T
    W = (Qr * (1.0 / np.sqrt(eig))) @ Qr.T
    Lw = LR @ S12
    Gw = W @ GR_raw
    PL = Lw.T @ np.linalg.solve(Lw @ Lw.T, Lw)
    PG = np.eye(2 * N) - _orth_projector(Gw)
    fro_diff = float(np.linalg.norm(PL - PG, "fro"))
    assert fro_diff < 1e-9, f"algebraic projector identity failed: {fro_diff}"
    gain_ann = float(np.max(np.abs(PL @ Gw)))
    assert gain_ann < 1e-9

    # finite-noise Gaussian calibration: 50 seeds x 2000 = 100000 draws
    cal_seeds = list(range(9001, 9051))
    reps = 2000
    ntot = len(cal_seeds) * reps
    Y = np.empty((2 * N, ntot))
    col = 0
    for seed in cal_seeds:
        rr = np.random.RandomState(seed)
        z = rr.standard_normal((2 * N, reps))
        eps = S12 @ z
        Y[:, col:col + reps] = PL @ eps
        col += reps
    ybar = Y.mean(axis=1)
    Cov_emp = (Y @ Y.T) / (ntot - 1.0) \
        - (ntot / (ntot - 1.0)) * np.outer(ybar, ybar)
    Cov_theory = PL @ Sigma @ PL
    cov_rel_err = float(np.linalg.norm(Cov_emp - Cov_theory, "fro")
                        / np.linalg.norm(Cov_theory, "fro"))
    evals, evecs = np.linalg.eigh(Cov_theory)
    rank = int(np.sum(evals > 1e-9 * float(np.trace(Cov_theory))))
    assert rank == 24, f"expected rank 24, got {rank}"
    idx = np.argsort(evals)[::-1][:rank]
    V = evecs[:, idx]
    Ppinv = V @ np.diag(1.0 / evals[idx]) @ V.T
    Qv = np.sum(Y * (Ppinv @ Y), axis=0)
    chi95 = float(scipy_chi2.ppf(0.95, df=rank))
    chi_df = {"mean": float(np.mean(Qv)), "std": float(np.std(Qv, ddof=1)),
              "p95": float(np.percentile(Qv, 95.0)),
              "coverage_chi2_95": float(np.mean(Qv <= chi95)),
              "theoretical_mean": float(rank), "theoretical_var": float(2 * rank),
              "chi2_95_quantile": chi95}

    # exact nonlinear cross-ratio gain invariance with finite random H
    n_inv = 0
    max_inv = 0.0
    any_den_flag_inv = False
    n_inv2 = 0
    max_inv2 = 0.0
    any_den_flag_inv2 = False
    for seed in range(9201, 9251):
        rr = np.random.RandomState(seed)
        for _ in range(200):
            H0 = _random_gain_vector(R * T, rr).reshape(R, T)
            a = _random_gain_vector(R, rr)
            b = _random_gain_vector(T, rr)
            C0, f0 = cross_ratio_matrix(H0)
            Hp = np.diag(a) @ H0 @ np.diag(b)
            C1, f1 = cross_ratio_matrix(Hp)
            any_den_flag_inv = any_den_flag_inv or bool(f0 or f1)
            max_inv = max(max_inv, float(np.max(np.abs(C1 - C0))))
            n_inv += 1
            # multiplicative log-noise version
            En = 0.05 * (rr.standard_normal((R, T))
                         + 1j * rr.standard_normal((R, T)))
            Hn = H0 * np.exp(En)
            Hnp = np.diag(a) @ Hn @ np.diag(b)
            Cn, fn = cross_ratio_matrix(Hn)
            Cn2, fn2 = cross_ratio_matrix(Hnp)
            any_den_flag_inv2 = any_den_flag_inv2 or bool(fn or fn2)
            max_inv2 = max(max_inv2, float(np.max(np.abs(Cn2 - Cn))))
            n_inv2 += 1

    # linearization consistency for three noise scales
    lin_rows = []
    for eps_noise in (0.02, 0.05, 0.1):
        rr = np.random.RandomState(9301)
        rel_sum = 0.0
        pl_rel_sum = 0.0
        n_draw = 10000
        for _ in range(n_draw):
            H0 = _random_gain_vector(R * T, rr).reshape(R, T)
            E = eps_noise * (rr.standard_normal((R, T))
                             + 1j * rr.standard_normal((R, T)))
            H = H0 * np.exp(E)
            C0, f0 = cross_ratio_matrix(H0)
            C, f1 = cross_ratio_matrix(H)
            if f0 or f1:
                raise AssertionError("unexpected near-zero denominator")
            delta = np.log(C / C0).ravel(order="C")
            lin_pred = Lc @ E.ravel(order="C")
            rel_sum += float(np.linalg.norm(delta - lin_pred)
                             / max(float(np.linalg.norm(lin_pred)), 1e-12))
            vd = _realify_vec(delta)
            vl = _realify_vec(lin_pred)
            pd = PL @ _lift_right(LR, vd)
            pll = PL @ _lift_right(LR, vl)
            pl_rel_sum += float(np.linalg.norm(pd - pll)
                                / max(float(np.linalg.norm(pll)), 1e-12))
        lin_rows.append({
            "eps_noise": float(eps_noise),
            "n_draws": int(n_draw),
            "rel_lin_err_mean": float(rel_sum / n_draw),
            "realified_pl_rel_err_mean": float(pl_rel_sum / n_draw),
        })
    lin_note = (
        "realified_pl_rel_err_mean: quotient differences are lifted back to "
        "the 2N log-H space with the right-inverse of LR, then compared "
        "through the 2N x 2N whitened projector PL; per-draw relative error, "
        "then averaged."
    )

    payload = {
        "algebraic_identity": {
            "R": R, "T": T, "K": K, "N": N, "real_dim": 2 * N,
            "fro_diff_PL_vs_PG": fro_diff,
            "max_gain_annihilation_err": gain_ann,
            "gain_subspace_real_dim": int(rankG),
            "cond_Sigma": float(np.linalg.cond(Sigma)),
            "rank_of_Cov_theory": int(rank),
        },
        "finite_noise_gaussian": {
            "seeds": [int(min(cal_seeds)), int(max(cal_seeds))],
            "reps_per_seed": reps,
            "total_draws": int(ntot),
            "cov_rel_err_fro": cov_rel_err,
            "rank": int(rank),
            "quadratic_form": chi_df,
            "chi2_df": int(rank),
            "max_gain_annihilation_finite_noise":
                float(np.max(np.abs(PL @ Gw))),
        },
        "exact_finite_noise_gain_invariance": {
            "n_draws": int(n_inv),
            "max_abs_invariance_err": float(max_inv),
            "any_near_zero_denominator": bool(any_den_flag_inv),
        },
        "log_noise_gain_invariance": {
            "log_noise_scale": 0.05,
            "n_draws": int(n_inv2),
            "max_abs_invariance_err": float(max_inv2),
            "any_near_zero_denominator": bool(any_den_flag_inv2),
        },
        "linearization_consistency": {"note": lin_note, "rows": lin_rows},
    }
    return payload, Qv


# ---------------------------------------------------------------------------
# PART 4: forward-model B4 verification
# ---------------------------------------------------------------------------
def _cg(A, rhs, x0, maxiter, tol=1e-30):
    """Plain conjugate gradient for symmetric positive-definite A."""
    x = np.array(x0, dtype=float).copy()
    r = rhs - A @ x
    p = r.copy()
    rsold = float(r @ r)
    iters = 0
    for it in range(1, maxiter + 1):
        Ap = A @ p
        denom = float(p @ Ap)
        if denom <= 0.0:
            break
        alpha = rsold / denom
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = float(r @ r)
        iters = it
        if rsnew <= tol or rsnew < 1e-24 * max(1.0, float(rhs @ rhs)):
            break
        p = r + (rsnew / rsold) * p
        rsold = rsnew
    return x, iters


def _bound_verification(M, b, w_chart, w_ls, Lrow, qval, chart_r, gamma):
    """Chart + dual-residual + ridge bookkeeping for one sweep point."""
    n = M.shape[0]
    Lt1 = (Lrow.T * qval).ravel()
    r = b - M @ w_chart
    # exact adjoint
    p_exact = np.linalg.solve(M.T, Lt1)
    target_ls = float((Lrow @ w_ls)[0])
    corrected = float((Lrow @ w_chart)[0] + p_exact @ r)
    identity_rel = abs(corrected - target_ls) / max(abs(target_ls), 1e-12)
    # p_tilde = 0 bound
    r_d0 = Lt1 - M.T @ np.zeros(n)
    bound_p0 = gamma * _norm(r_d0) * _norm(r)
    actual_p0 = abs(float((Lrow @ (w_ls - w_chart))[0]))
    margin_p0 = float(bound_p0 - actual_p0)
    # CG-approximate adjoint (maxiter 8)
    p_cg, cg_iters = _cg(M, Lt1, np.zeros(n), maxiter=8)
    r_d_cg = Lt1 - M.T @ p_cg
    bound_cg = gamma * _norm(r_d_cg) * _norm(r)
    actual_cg = abs(float((Lrow @ (w_ls - w_chart))[0]) - p_cg @ r)
    return {
        "identity_rel": float(identity_rel),
        "bound_p0": float(bound_p0),
        "actual_p0": float(actual_p0),
        "bound_margin_p0": float(margin_p0),
        "bound_p0_holds": bool(actual_p0 <= bound_p0 + 1e-8),
        "cg_iters_used": int(cg_iters),
        "bound_cg": float(bound_cg),
        "actual_cg": float(actual_cg),
        "bound_cg_holds": bool(actual_cg <= bound_cg + 1e-8),
    }


def _run_sweep(M_list, b_list, labels, Lrow, qval, r_chart, extra):
    rows = []
    for lab, M, b in zip(labels, M_list, b_list):
        n = M.shape[0]
        sv = np.linalg.svd(M, compute_uv=False)
        sigma_min = float(np.min(sv))
        gamma = 1.0 / sigma_min
        cond_M = float(sigma_min and np.max(sv) / sigma_min)
        w_ls = np.linalg.solve(M, b)
        evals, evecs = np.linalg.eigh(M)
        U = evecs[:, -r_chart:]
        w_chart = U @ np.linalg.solve(U.T @ M @ U, U.T @ b)
        lam_r = 0.05 * float(np.mean(np.diag(M)))
        w_ridge = np.linalg.solve(M + lam_r * np.eye(n), b)
        ls_err = abs(float((Lrow @ w_ls)[0]))
        chart_err = abs(float((Lrow @ w_chart)[0]))
        p_ex_vec = np.linalg.solve(M.T, (Lrow.T * qval).ravel())
        corr = float((Lrow @ w_chart)[0]
                     + p_ex_vec @ (b - M @ w_chart))
        a3_err = abs(corr)
        ridge_err = abs(float((Lrow @ w_ridge)[0]))
        bv = _bound_verification(M, b, w_chart, w_ls, Lrow, qval,
                                 r_chart, gamma)
        row = dict(extra(lab))
        row.update({
            "cond_M": cond_M, "gamma": gamma,
            "ls_err": ls_err, "chart_err": chart_err,
            "a3_corrected_err": a3_err, "ridge_err": ridge_err,
            "n": int(n), "r": int(r_chart),
        })
        row.update(bv)
        rows.append(row)
    return rows


def _rot2(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s], [s, c]])


def _se2_jac_blocks(i, j, poses):
    """Jacobian of predicted relative pose wrt pose perturbations (per spec)."""
    th_i = poses[i, 2]
    p_i = poses[i, :2]
    p_j = poses[j, :2]
    Ri = _rot2(th_i)
    Rj = _rot2(th_j := poses[j, 2])
    Rrel = Ri.T @ Rj
    cd, sd = float(Rrel[0, 0]), float(Rrel[1, 0])
    p_rel = Ri.T @ (p_j - p_i)
    dxr, dyr = float(p_rel[0]), float(p_rel[1])
    Ji = np.array([[-cd, -sd, dyr],
                   [sd, -cd, -dxr],
                   [0.0, 0.0, -1.0]])
    Jj = np.array([[cd, sd, 0.0],
                   [-sd, cd, 0.0],
                   [0.0, 0.0, 1.0]])
    return Ji, Jj, th_j - th_i, p_rel


def _pose_graph_systems():
    """Build the five alpha_lc pose-graph systems shared by PARTs 4 and 5."""
    n_poses = 20
    ang = 2.0 * np.pi * np.arange(n_poses) / n_poses
    radius = 5.0
    poses = np.column_stack([radius * np.cos(ang),
                             radius * np.sin(ang),
                             ang + np.pi / 2.0])
    sigma_t = 0.02
    sigma_r = 0.01
    odo_pairs = [(i - 1, i) for i in range(1, n_poses)]
    cand = [(i, j) for i in range(n_poses) for j in range(i + 2, n_poses)
            if min(j - i, n_poses - (j - i)) >= 3]
    rs_pairs = np.random.RandomState(8121)
    perm = rs_pairs.permutation(len(cand))
    lc_pairs = [cand[k] for k in perm[:8]]
    noise_rs = np.random.RandomState(8122)

    def rel_from_pairs(pairs):
        Js, rs = [], []
        for i, j in pairs:
            Ji, Jj, threl, p_rel = _se2_jac_blocks(i, j, poses)
            pred = np.array([p_rel[0], p_rel[1], threl])
            meas = pred + np.array([
                noise_rs.standard_normal() * sigma_t,
                noise_rs.standard_normal() * sigma_t,
                noise_rs.standard_normal() * sigma_r,
            ])
            Jr = np.zeros((3, 3 * n_poses))
            Jr[:, 3 * i:3 * i + 3] = Ji
            Jr[:, 3 * j:3 * j + 3] = Jj
            Js.append(Jr)
            rs.append(meas - pred)
        return np.vstack(Js), np.concatenate(rs)

    J_odo_full, r_odo = rel_from_pairs(odo_pairs)
    J_lc_full, r_lc = rel_from_pairs(lc_pairs)
    J_odo = J_odo_full[:, 3:]
    J_lc = J_lc_full[:, 3:]
    n = J_odo.shape[1]
    assert n == 57
    target_col = 3 * 19 + 0 - 3
    assert target_col == 54
    Lrow = np.zeros((1, n))
    Lrow[0, target_col] = 1.0
    qval = 1.0
    Ms, bs, labels = [], [], []
    for alpha_lc in (1.0, 0.1, 0.01, 0.001, 0.0001):
        M = J_odo.T @ J_odo + alpha_lc * (J_lc.T @ J_lc)
        b = J_odo.T @ r_odo + alpha_lc * (J_lc.T @ r_lc)
        Ms.append(M)
        bs.append(b)
        labels.append(alpha_lc)
    meta = {
        "description": "linearized 2D SE(2) pose graph, pose-0 anchor dropped",
        "n_poses": n_poses, "n": n, "r_chart": 10,
        "n_odometry": len(odo_pairs), "n_loop_closure": len(lc_pairs),
        "loop_closure_pairs": [list(p) for p in lc_pairs],
        "pair_seed": 8121, "noise_seed": 8122,
        "target": "x-coordinate of pose 19 (state col 54)",
    }
    return {"meta": meta, "Ms": Ms, "bs": bs, "labels": labels,
            "Lrow": Lrow, "qval": qval}


def part4a_pose_graph():
    s = _pose_graph_systems()
    sweep = _run_sweep(s["Ms"], s["bs"], s["labels"], s["Lrow"], s["qval"],
                       10, lambda a: {"alpha_lc": float(a)})
    meta = dict(s["meta"])
    meta["sweep"] = sweep
    meta["summary_flags"] = {
        "all_identity_rel_lt_1e-8":
            bool(max(r["identity_rel"] for r in sweep) < 1e-8),
        "all_bound_p0_holds": bool(all(r["bound_p0_holds"]
                                       for r in sweep)),
        "all_bound_cg_holds": bool(all(r["bound_cg_holds"]
                                       for r in sweep)),
    }
    return meta


def _lookat_camera(phi):
    t = 8.0 * np.array([np.cos(phi), np.sin(phi), 0.0])
    f = -t / _norm(t)
    up = np.array([0.0, 0.0, 1.0])
    x = np.cross(up, f)
    x = x / _norm(x)
    y = np.cross(f, x)
    R = np.column_stack([x, y, f])
    return t, R


def _bundle_systems():
    """Build the five lambda_p bundle systems shared by PARTs 4 and 5."""
    n_cams = 5
    m_points = 12
    rng = np.random.RandomState(8123)
    f = 1.0
    sigma_pix = 0.02
    cams = [_lookat_camera(2.0 * np.pi * c / n_cams) for c in range(n_cams)]
    pts = 2.0 * (rng.random_sample((m_points, 3)) - 0.5) * 2.0
    obs = []
    for c in range(n_cams):
        t_c, R_c = cams[c]
        RcT = R_c.T
        for mm in range(m_points):
            x = RcT @ (pts[mm] - t_c)
            if x[2] <= 1e-9:
                continue
            if abs(x[0] / x[2]) >= 2.0 or abs(x[1] / x[2]) >= 2.0:
                continue
            u_t = f * x[0] / x[2]
            v_t = f * x[1] / x[2]
            u_o = u_t + rng.standard_normal() * sigma_pix
            v_o = v_t + rng.standard_normal() * sigma_pix
            obs.append((c, mm, u_o, v_o, u_t, v_t))
    # unknown columns: camera translations c=1..4 (12), then points (36)
    n_cam_unknown = (n_cams - 1) * 3
    n = n_cam_unknown + m_points * 3
    assert n == 48
    J = np.zeros((2 * len(obs), n))
    res = np.empty(2 * len(obs))
    for row0, (c, mm, u_o, v_o, u_t, v_t) in enumerate(obs):
        t_c, R_c = cams[c]
        RcT = R_c.T
        x = RcT @ (pts[mm] - t_c)
        x0, x1, x2 = float(x[0]), float(x[1]), float(x[2])
        du = np.array([1.0 / x2, 0.0, -x0 / (x2 * x2)]) * f
        dv = np.array([0.0, 1.0 / x2, -x1 / (x2 * x2)]) * f
        du_dp = du @ RcT
        dv_dp = dv @ RcT
        p_cols = np.arange(n_cam_unknown + 3 * mm,
                           n_cam_unknown + 3 * mm + 3)
        J[2 * row0, p_cols] = du_dp
        J[2 * row0 + 1, p_cols] = dv_dp
        if c >= 1:
            c_cols = np.arange(3 * (c - 1), 3 * (c - 1) + 3)
            J[2 * row0, c_cols] = -du_dp
            J[2 * row0 + 1, c_cols] = -dv_dp
        res[2 * row0] = u_o - u_t
        res[2 * row0 + 1] = v_o - v_t
    target_col = n_cam_unknown + 0 * 3 + 2
    assert target_col == 14
    Lrow = np.zeros((1, n))
    Lrow[0, target_col] = 1.0
    qval = 1.0
    Ms, bs, labels = [], [], []
    for lam_p in (1.0, 0.1, 0.01, 0.001, 0.0001):
        M = J.T @ J + lam_p * np.eye(n)
        b = J.T @ res
        Ms.append(M)
        bs.append(b)
        labels.append(lam_p)
    meta = {
        "description": "toy linearized 3D bundle adjustment, camera 0 anchored",
        "n_cams": n_cams, "m_points": m_points, "n": n, "r_chart": 10,
        "n_observations": len(obs),
        "target": "z-coordinate of point 0 (state col 14)",
    }
    return {"meta": meta, "Ms": Ms, "bs": bs, "labels": labels,
            "Lrow": Lrow, "qval": qval}


def part4b_bundle():
    s = _bundle_systems()
    sweep = _run_sweep(s["Ms"], s["bs"], s["labels"], s["Lrow"], s["qval"],
                       10, lambda a: {"lambda_p": float(a)})
    meta = dict(s["meta"])
    meta["sweep"] = sweep
    meta["summary_flags"] = {
        "all_identity_rel_lt_1e-8":
            bool(max(r["identity_rel"] for r in sweep) < 1e-8),
        "all_bound_p0_holds": bool(all(r["bound_p0_holds"]
                                       for r in sweep)),
        "all_bound_cg_holds": bool(all(r["bound_cg_holds"]
                                       for r in sweep)),
    }
    return meta


# ---------------------------------------------------------------------------
# PART 5: near-singular forward-model B4 control
# ---------------------------------------------------------------------------
def _near_singular_control_row(task, base_label, M, b, eps_out=0.5):
    """Small-residual control: p_tilde=0, L-row on the near-singular v_min."""
    n = M.shape[0]
    U_svd, s, Vh = np.linalg.svd(M)
    sigma_min = float(s[-1])
    sigma_max = float(s[0])
    cond_M = float(sigma_max / sigma_min)
    v_min = Vh[-1, :]
    u_min = U_svd[:, -1]
    gamma = 1.0 / sigma_min

    w_ls = np.linalg.solve(M, b)
    w_tilde = w_ls - eps_out * v_min
    Lrow = v_min.reshape(1, -1)          # row vector targeting v_min
    qval = 1.0
    p_tilde = np.zeros(n)

    r = b - M @ w_tilde
    residual_norm = float(np.linalg.norm(r))
    r_d = Lrow.T - (M.T @ p_tilde)[:, None]   # equals v_min as a column
    bound = gamma * float(np.linalg.norm(r_d)) * residual_norm
    actual = abs(float(qval * (Lrow @ (w_ls - w_tilde))[0]
                       - float(p_tilde @ r)))

    p_exact = np.linalg.solve(M.T, Lrow.T)

    # Dual identity on the NONZERO correction functional:
    #   q^T L (w_ls - w_tilde) = eps_out * q^T L v_min
    #                           = p_exact^T r   (exact adjoint correction)
    correction_true = qval * float((Lrow @ (w_ls - w_tilde))[0])
    correction_est = float(np.squeeze(p_exact) @ r)
    correction_abs = abs(correction_est - correction_true)
    correction_rel = correction_abs / max(abs(correction_true), 1e-12)

    # Full-functional identity, kept for the audit record:
    #   q^T L w_tilde + p_exact^T r == q^T L w_ls
    full_corrected = (qval * float((Lrow @ w_tilde)[0])
                      + float(np.squeeze(p_exact) @ r))
    full_target_ls = qval * float((Lrow @ w_ls)[0])
    full_identity_abs = abs(full_corrected - full_target_ls)
    full_identity_rel = full_identity_abs / max(abs(full_target_ls), 1e-12)
    identity_ok = bool(correction_rel < 1e-8 and full_identity_abs < 1e-9)

    return {
        "task": task,
        "base_row": base_label,
        "n": int(n),
        "eps_out": float(eps_out),
        "cond_M": cond_M,
        "sigma_min": sigma_min,
        "gamma": gamma,
        "residual_norm": residual_norm,
        "output_error": float(actual),
        "small_residual_certifies":
            bool(actual <= residual_norm + 1e-12),
        "certified_bound": bound,
        "bound_holds": bool(actual <= bound + 1e-8),
        "correction_true": correction_true,
        "correction_est": correction_est,
        "correction_abs": correction_abs,
        "correction_rel": correction_rel,
        "full_identity_abs": full_identity_abs,
        "full_identity_rel": full_identity_rel,
        "target_ls": full_target_ls,
        "identity_ok": identity_ok,
    }


def part5_near_singular():
    pg = _pose_graph_systems()
    ba = _bundle_systems()
    # Base rows fixed by PART 4 sweep order:
    # pose graph alpha_lc=1.0 is index 0; bundle lambda_p=0.001 is index 3.
    pose_graph = _near_singular_control_row(
        "pose_graph", "alpha_lc=1.0 (PART 4 sweep index 0)",
        pg["Ms"][0], pg["bs"][0])
    bundle = _near_singular_control_row(
        "bundle", "lambda_p=0.001 (PART 4 sweep index 3)",
        ba["Ms"][3], ba["bs"][3])
    return {
        "pose_graph": pose_graph,
        "bundle": bundle,
        "summary_flags": {
            "all_bound_holds": bool(pose_graph["bound_holds"]
                                    and bundle["bound_holds"]),
            "all_small_residual_controls_fail":
                bool(not pose_graph["small_residual_certifies"]
                     and not bundle["small_residual_certifies"]),
            "all_identity_ok": bool(pose_graph["identity_ok"]
                                    and bundle["identity_ok"]),
            "max_correction_rel": float(max(pose_graph["correction_rel"],
                                            bundle["correction_rel"])),
            "max_full_identity_abs": float(max(
                pose_graph["full_identity_abs"],
                bundle["full_identity_abs"])),
        },
    }


def part4_forward_model():
    return {
        "pose_graph": part4a_pose_graph(),
        "bundle": part4b_bundle(),
    }


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def _plot_b3(b3, path):
    g = [r for r in b3["configs"] if r["noise_model"] == "gaussian"]
    h = [r for r in b3["configs"] if r["noise_model"] == "heavy_t3"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axes[0]
    for r in g:
        ax.errorbar(r["bound_sigma"], r["p_hat"],
                    yerr=[[r["p_hat"] - r["boot_lo"]],
                          [r["boot_hi"] - r["p_hat"]]],
                    fmt="o", ms=5, capsize=2,
                    label=(f"{r['geometry']} b={r['beta']} s={r['sigma']}"))
    lims = [0.0, max(max(r["bound_sigma"], r["boot_hi"])
                     for r in g) * 1.08]
    ax.plot(lims, lims, "k--", lw=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("bound_sigma")
    ax.set_ylabel("p_hat (error rate)")
    ax.set_title("B3 Gaussian calibration (finite noise)")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)
    ax2 = axes[1]
    for r in h:
        uc = "o" if r["covered_sigma1"] else "v"
        ax2.plot(r["bound_sigma1"], r["p_hat"], uc, ms=7,
                 label=(f"{r['geometry']} b={r['beta']} s={r['sigma']}"
                        + ("" if r["covered_sigma1"] else "  (under)")))
    lims2 = [0.0, max(max(r["bound_sigma1"], r["p_hat"]) for r in h) * 1.1]
    ax2.plot(lims2, lims2, "k--", lw=1)
    ax2.set_xlim(lims2)
    ax2.set_ylim(lims2)
    ax2.set_xlabel("bound_sigma1 (Gaussian claim)")
    ax2.set_ylabel("p_hat (error rate)")
    ax2.set_title("B3 heavy-tailed t(3) vs Gaussian bound")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _plot_b2(q_out, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    df = q_out["rank"]
    qf_in = q_out["quadratic_form"]
    Qv = np.asarray(q_out.get("_Q_sample", []))
    if Qv.size == 0:
        ax.text(0.5, 0.5, "Q sample unavailable", ha="center")
    else:
        ax.hist(Qv, bins=90, density=True, alpha=0.6, label="empirical Q")
        xs = np.linspace(0, max(np.percentile(Qv, 99.9), df + 4.0 * np.sqrt(2 * df)),
                         400)
        ax.plot(xs, scipy_chi2.pdf(xs, df=df), "r-", lw=1.8,
                label=f"chi2(df={df}) pdf")
        q95 = scipy_chi2.ppf(0.95, df=df)
        ax.axvline(q95, color="k", ls="--", lw=1.2,
                   label=f"chi2 95%={q95:.2f}")
        ax.set_title(f"Q distribution: mean={qf_in['mean']:.2f} "
                     f"coverage={qf_in['coverage_chi2_95']*100:.1f}%")
    ax.set_xlabel("Q = y^T Cov^+ y")
    ax.set_ylabel("density")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _plot_b4(b4, path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, key, xkey in ((axes[0], "pose_graph", "alpha_lc"),
                          (axes[1], "bundle", "lambda_p")):
        rows = b4[key]["sweep"]
        xs = [r["cond_M"] for r in rows]
        for lbl, kk in (("LS", "ls_err"), ("chart(r=10)", "chart_err"),
                        ("A3 corrected", "a3_corrected_err"),
                        ("ridge", "ridge_err")):
            ax.semilogx(xs, [r[kk] for r in rows], "o-", lw=1.4,
                        label=lbl)
        for i, r in enumerate(rows):
            ax.annotate(f"id={r['identity_rel']:.1e}\n"
                        f"cg={str(r['bound_cg_holds'])[0]}",
                        (xs[i], r["a3_corrected_err"]),
                        textcoords="offset points", xytext=(4, -12),
                        fontsize=6, color="0.25")
        ax.set_xlabel(f"cond(M)   ({xkey} decreasing left->right)")
        ax.set_ylabel("|q^T L (w - 0)| (truth = 0)")
        ax.set_title(f"B4 {key}: chart vs corrected vs LS vs ridge")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _plot_b5(b5, path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
    for ax, key in ((axes[0], "pose_graph"), (axes[1], "bundle")):
        row = b5[key]
        names = ["residual\nnorm", "output\nerror", "certified\nbound"]
        vals = [row["residual_norm"], row["output_error"],
                row["certified_bound"]]
        colors = ["tab:blue", "tab:red", "tab:green"]
        bars = ax.bar(names, vals, color=colors, alpha=0.75)
        ax.set_yscale("log")
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2.0, v,
                    f"{v:.3e}", ha="center", va="bottom", fontsize=9)
        ax.set_title(f"{key}: cond(M)={row['cond_M']:.4g}, "
                     f"sigma_min={row['sigma_min']:.3e}")
        ax.set_ylabel("value (log scale)")
        ax.grid(alpha=0.3, which="both", axis="y")
        ax.set_ylim(bottom=max(1e-12, min(vals) * 0.3))
        ax.text(0.03, 0.03,
                f"identity (correction functional):\n"
                f"correction_rel={row['correction_rel']:.2e}, "
                f"full_identity_abs={row['full_identity_abs']:.2e}, "
                f"identity_ok={row['identity_ok']}",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", ec="0.6",
                          alpha=0.85))
    fig.suptitle("B4 near-singular control: residual << output error, but "
                 "certified bound == output error (gamma bound holds); "
                 "identity verified on the nonzero correction functional",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _save_plots(results):
    saved = []
    if not MPL_OK:
        return saved
    try:
        p1 = os.path.join(HERE, "b3_calibration.png")
        _plot_b3(results["b3_calibration"], p1)
        saved.append(os.path.basename(p1))
    except Exception as exc:  # pragma: no cover
        print(f"  [warn] b3 plot failed: {exc}")
    try:
        p2 = os.path.join(HERE, "b2_ratio_sufficiency.png")
        payload = results["b2_finite_noise"]["finite_noise_gaussian"]
        payload["_Q_sample"] = results["_q_sample"]
        _plot_b2(payload, p2)
        saved.append(os.path.basename(p2))
        payload.pop("_Q_sample", None)
    except Exception as exc:  # pragma: no cover
        print(f"  [warn] b2 plot failed: {exc}")
    try:
        p3 = os.path.join(HERE, "b4_forward_model.png")
        _plot_b4(results["b4_forward_model"], p3)
        saved.append(os.path.basename(p3))
    except Exception as exc:  # pragma: no cover
        print(f"  [warn] b4 plot failed: {exc}")
    try:
        p5 = os.path.join(HERE, "b4_near_singular.png")
        _plot_b5(results["b4_near_singular"], p5)
        saved.append(os.path.basename(p5))
    except Exception as exc:  # pragma: no cover
        print(f"  [warn] b4_near_singular plot failed: {exc}")
    return saved


# ---------------------------------------------------------------------------
# summary writer
# ---------------------------------------------------------------------------
def _fmt_table(rows, key_label, err_keys):
    head = ("| " + key_label + " | cond(M) | LS err | chart err | A3 corr | "
            "ridge err | identity rel | bound_cg holds |")
    sep = "|---|---|---|---|---|---|---|---|"
    lines = [head, sep]
    for r in rows:
        lines.append(
            f"| {r[key_label]:.6g} | {r['cond_M']:.3g} "
            f"| {r['ls_err']:.3g} | {r['chart_err']:.3g} "
            f"| {r['a3_corrected_err']:.3g} | {r['ridge_err']:.3g} "
            f"| {r['identity_rel']:.2e} | {r['bound_cg_holds']} |")
    return "\n".join(lines)


def write_summary(results):
    b3s = results["b3_calibration"]["summary"]
    b2 = results["b2_finite_noise"]
    qf = b2["finite_noise_gaussian"]["quadratic_form"]
    sh = results["b3_reused_validation_sharpened"]
    pg = results["b4_forward_model"]["pose_graph"]
    ba = results["b4_forward_model"]["bundle"]
    b5 = results["b4_near_singular"]
    gauss_bad = len(b3s["gaussian_uncovered_sigma1"])
    heavy_bad = len(b3s["heavy_uncovered_sigma1"])
    a3_note = (
        "A3 dual-residual correction equals the full-LS functional "
        "(identity), so it provides no endpoint-accuracy gain over LS; it "
        "removes the subspace/chart error. Chart-only and ridge are "
        "regularized estimators and can be smaller or larger than LS "
        "depending on the setting."
    )
    lines = [
        "# A3 extended calibration study - results",
        "",
        f"Generated by {results['generation']['generator']} "
        f"(numpy {results['generation']['numpy']}, "
        f"scipy {results['generation']['scipy']}, "
        f"cpu_threads={results['generation']['cpu_threads']}).",
        "",
        "## PART 0 frozen baseline audit",
        "",
        f"- counts = {results['frozen_baseline']['frozen_counts']}; "
        f"frozen_baseline_ok = {results['frozen_baseline']['frozen_baseline_ok']}.",
        "",
        "## PART 1 B3 finite-noise branch-bound calibration",
        "",
        f"- Gaussian configs: {b3s['gaussian_covered_sigma_count']}/"
        f"{b3s['gaussian_configs_total']} covered at matched sigma, "
        f"{b3s['gaussian_covered_sigma1_count']}/"
        f"{b3s['gaussian_configs_total']} covered vs sigma=1 bound.",
        f"- Gaussian not covered vs sigma1 ({gauss_bad}): "
        f"{b3s['gaussian_uncovered_sigma1']}.",
        f"- Heavy-tail configs: {b3s['heavy_configs_total']} total; "
        f"strong bootstrap violations vs matched sigma: "
        f"{b3s['heavy_strong_violation_vs_matched_sigma_count']}, "
        f"vs sigma=1 bound: "
        f"{b3s['heavy_strong_violation_vs_sigma1_count']}.",
        f"- Heavy-tail not covered vs sigma1 ({heavy_bad}): "
        f"{b3s['heavy_uncovered_sigma1']}.",
        "",
        "## PART 2 sharpened B3 reused-validation (3E) control",
        "",
        f"- {sh['total_reps']} reps, empirical error = "
        f"{sh['empirical_error_by_construction']} by construction.",
        f"- Bound < 1 in {sh['violation_count']} reps "
        f"(fraction {sh['fraction_bound_lt_1']:.4f}); avg plug-in bound "
        f"{sh['avg_plugin_bound']:.4f} "
        f"(min {sh['min_plugin_bound']:.4f}, max {sh['max_plugin_bound']:.4f}).",
        "",
        "## PART 3 B2 finite-noise ratio sufficiency",
        "",
        f"- Projector identity fro_diff = {b2['algebraic_identity']['fro_diff_PL_vs_PG']:.3e}; "
        f"rank(Cov_theory) = {b2['finite_noise_gaussian']['rank']} "
        f"(rank of {b2['algebraic_identity']['gain_subspace_real_dim']} gain "
        f"directions annihilated).",
        f"- Finite-noise covariance rel Frobenius err = "
        f"{b2['finite_noise_gaussian']['cov_rel_err_fro']:.3e} "
        f"over {b2['finite_noise_gaussian']['total_draws']} draws.",
        f"- Q: empirical mean {qf['mean']:.3f} (theory {qf['theoretical_mean']}), "
        f"std {qf['std']:.3f} (theory sqrt(48)={np.sqrt(48):.3f}), "
        f"95th pct {qf['p95']:.3f} vs chi2 95% {qf['chi2_95_quantile']:.3f}; "
        f"coverage {qf['coverage_chi2_95']*100:.2f}%.",
        f"- Nonlinear gain invariance exact max err = "
        f"{b2['exact_finite_noise_gain_invariance']['max_abs_invariance_err']:.3e}; "
        f"multiplicative-noise version = "
        f"{b2['log_noise_gain_invariance']['max_abs_invariance_err']:.3e}.",
        f"- Log-H exact linearity (not an approximation): for multiplicative "
        f"log-noise H = H0*exp(E), log C(H) - log C(H0) = Lc*vec(E) exactly "
        f"(no higher-order terms), because the log cross-ratio map is exactly "
        f"linear in log-H space. The ~1e-15 mean rel errs below are therefore "
        f"exact-identity roundoff, not linearization error: "
        + "; ".join(
            f"eps={r['eps_noise']}: {r['rel_lin_err_mean']:.3e} "
            f"(realified PL {r['realified_pl_rel_err_mean']:.3e})"
            for r in b2["linearization_consistency"]["rows"])
        + ".",
        "",
        "## PART 4 forward-model B4 verification",
        "",
        "### Pose graph (2D SE(2), anchor 0 dropped)",
        "",
        f"- identity max rel = {max(r['identity_rel'] for r in pg['sweep']):.2e}; "
        f"p0 bound holds all = {pg['summary_flags']['all_bound_p0_holds']}; "
        f"CG bound holds all = {pg['summary_flags']['all_bound_cg_holds']}.",
        "",
        _fmt_table(pg["sweep"], "alpha_lc", ()),
        "",
        a3_note,
        "",
        "### Bundle adjustment (3D pinhole, camera 0 anchored)",
        "",
        f"- identity max rel = {max(r['identity_rel'] for r in ba['sweep']):.2e}; "
        f"p0 bound holds all = {ba['summary_flags']['all_bound_p0_holds']}; "
        f"CG bound holds all = {ba['summary_flags']['all_bound_cg_holds']}.",
        "",
        _fmt_table(ba["sweep"], "lambda_p", ()),
        "",
        a3_note,
        "",
        "## PART 5 near-singular forward-model B4 control",
        "",
        f"- L-row = v_min (right singular vector of the smallest singular "
        f"value), q = 1, p_tilde = 0, eps_out = 0.5; base systems: pose graph "
        f"alpha_lc=1.0 and bundle lambda_p=0.001.",
        f"- summary_flags = {b5['summary_flags']}.",
        "",
        "| task | cond(M) | sigma_min | gamma | residual_norm | output_error "
        "| certified_bound | bound_holds | correction_rel | "
        "full_identity_abs | identity_ok |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for key in ("pose_graph", "bundle"):
        r = b5[key]
        lines.append(
            f"| {r['task']} | {r['cond_M']:.3g} | {r['sigma_min']:.3e} "
            f"| {r['gamma']:.3e} | {r['residual_norm']:.3e} "
            f"| {r['output_error']:.3e} | {r['certified_bound']:.3e} "
            f"| {r['bound_holds']} | {r['correction_rel']:.2e} "
            f"| {r['full_identity_abs']:.2e} | {r['identity_ok']} |")
    lines += [
        "",
        "Expected-fail note: residual_norm << output_error in both tasks, so "
        "a small residual alone is not an endpoint-output certificate "
        "(small_residual_certifies = False is the expected control outcome). "
        "The certified gamma bound holds with certified_bound == output_error, "
        "confirming the verified bound, not the residual size, is the "
        "certificate.",
        "",
    ]
    id_parts = []
    for key in ("pose_graph", "bundle"):
        r = b5[key]
        id_parts.append(
            f"{r['task']}: correction_rel={r['correction_rel']:.3e}, "
            f"full_identity_abs={r['full_identity_abs']:.3e}, "
            f"identity_ok={r['identity_ok']}")
    lines += [
        f"- Adjoint-identity note: the dual-residual identity is verified on "
        f"the NONZERO correction functional "
        f"q^T L(w_ls - w_tilde) = eps_out = 0.5, not on q^T L w_ls, whose "
        f"bundle value is machine-zero because b is orthogonal to the damped "
        f"null direction. correction_rel compares p_exact^T r to "
        f"q^T L(w_ls - w_tilde) (target 0.5). The full-functional absolute "
        f"identity |q^T L w_tilde + p_exact^T r - q^T L w_ls| is also "
        f"recorded. Results: {'/'.join(id_parts)}.",
        "",
        "## Overall flags",
        "",
        "```json",
        json.dumps(results["overall_flags"], indent=2),
        "```",
        "",
    ]
    path = os.path.join(HERE, "calibration_summary.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return path


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    t_start = time.time()
    print("PART 0: frozen baseline audit")
    frozen = part0_frozen_audit()
    print(f"  frozen baseline ok={frozen['frozen_baseline_ok']} "
          f"counts={frozen['frozen_counts']}")

    print("PART 1: B3 finite-noise calibration "
          "(50 seeds x 2000 draws/config)")
    b3 = part1_b3_calibration()
    b3s = b3["summary"]
    print(f"  gaussian covered(matched sigma)="
          f"{b3s['gaussian_covered_sigma_count']}/"
          f"{b3s['gaussian_configs_total']}; "
          f"heavy strong violations(vs sigma1)="
          f"{b3s['heavy_strong_violation_vs_sigma1_count']}/"
          f"{b3s['heavy_configs_total']}")

    print("PART 2: sharpened B3 reused-validation control")
    sharp = part2_sharpened_3e()
    print(f"  violation_count={sharp['violation_count']} "
          f"avg_plugin_bound={sharp['avg_plugin_bound']:.4f}")

    print("PART 3: B2 finite-noise ratio sufficiency")
    b2, q_sample = part3_b2_finite_noise()
    qf = b2["finite_noise_gaussian"]["quadratic_form"]
    print(f"  fro_diff={b2['algebraic_identity']['fro_diff_PL_vs_PG']:.3e} "
          f"cov_rel_err={b2['finite_noise_gaussian']['cov_rel_err_fro']:.3e} "
          f"Q mean={qf['mean']:.3f} coverage={qf['coverage_chi2_95']:.4f}")

    print("PART 4: forward-model B4 verification")
    b4 = part4_forward_model()
    for key in ("pose_graph", "bundle"):
        row0 = b4[key]["sweep"][0]
        flags = b4[key]["summary_flags"]
        print(f"  {key}: identity_rel={row0['identity_rel']:.2e} "
              f"flags={flags}")

    print("PART 5: near-singular forward-model B4 control")
    b5 = part5_near_singular()
    for key in ("pose_graph", "bundle"):
        r = b5[key]
        print(f"  {key}: cond={r['cond_M']:.4g} "
              f"sigma_min={r['sigma_min']:.3e} "
              f"gamma={r['gamma']:.4g} "
              f"residual_norm={r['residual_norm']:.3e} "
              f"output_error={r['output_error']:.3e} "
              f"bound={r['certified_bound']:.3e} "
              f"bound_holds={r['bound_holds']} "
              f"small_residual_certifies={r['small_residual_certifies']} "
              f"correction_rel={r['correction_rel']:.2e} "
              f"full_identity_abs={r['full_identity_abs']:.2e} "
              f"identity_ok={r['identity_ok']}")

    results = {
        "generation": {
            "generator": os.path.basename(__file__),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpu_threads": 1,
        },
        "frozen_baseline": frozen,
        "b3_calibration": b3,
        "b3_reused_validation_sharpened": sharp,
        "b2_finite_noise": b2,
        "b4_forward_model": b4,
        "b4_near_singular": b5,
    }

    pgf = b4["pose_graph"]["summary_flags"]
    baf = b4["bundle"]["summary_flags"]
    b5f = b5["summary_flags"]
    overall_flags = {
        "frozen_baseline_ok": bool(frozen["frozen_baseline_ok"]),
        "b2_algebraic_identity_ok":
            bool(b2["algebraic_identity"]["fro_diff_PL_vs_PG"] < 1e-9),
        "b2_rank_ok":
            bool(b2["finite_noise_gaussian"]["rank"] == 24),
        "b2_chi2_coverage_in_0p90_0p99":
            bool(0.90 <= qf["coverage_chi2_95"] <= 0.99),
        "b3_gaussian_all_covered_matched_sigma":
            bool(b3s["gaussian_covered_sigma_count"]
                 == b3s["gaussian_configs_total"]),
        "b3_gaussian_all_covered_sigma1":
            bool(b3s["gaussian_covered_sigma1_count"]
                 == b3s["gaussian_configs_total"]),
        "b3_heavy_violations_are_findings":
            True,
        "b3_heavy_strong_violation_vs_sigma1_count":
            int(b3s["heavy_strong_violation_vs_sigma1_count"]),
        "b3_sharpened_3e_violation_count": int(sharp["violation_count"]),
        "b4_pose_graph_identity_ok": bool(pgf["all_identity_rel_lt_1e-8"]),
        "b4_pose_graph_bound_p0_ok": bool(pgf["all_bound_p0_holds"]),
        "b4_pose_graph_bound_cg_ok": bool(pgf["all_bound_cg_holds"]),
        "b4_bundle_identity_ok": bool(baf["all_identity_rel_lt_1e-8"]),
        "b4_bundle_bound_p0_ok": bool(baf["all_bound_p0_holds"]),
        "b4_bundle_bound_cg_ok": bool(baf["all_bound_cg_holds"]),
        "b4_near_singular_bound_holds": bool(b5f["all_bound_holds"]),
        "b4_near_singular_small_residual_controls_fail":
            bool(b5f["all_small_residual_controls_fail"]),
        "b4_near_singular_identity_ok":
            bool(b5f["all_identity_ok"]),
        "plots_saved": None,
    }
    results["overall_flags"] = overall_flags
    results["_q_sample"] = q_sample.tolist()

    plots = _save_plots(results)
    results.pop("_q_sample", None)
    overall_flags["plots_saved"] = plots

    json_path = os.path.join(HERE, "calibration_results.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(_clean_jsonable(results), fh, indent=2)
        fh.write("\n")
    md_path = write_summary(results)

    fails = [
        key for key, val in overall_flags.items()
        if isinstance(val, bool) and key.startswith(("b2_", "b4_"))
        and not val
    ]
    print("\nfinal audit")
    print(f"  frozen_baseline_ok={overall_flags['frozen_baseline_ok']}")
    print(f"  b2 algebraic/rank/coverage flags ok (see json): "
          f"identity={overall_flags['b2_algebraic_identity_ok']}, "
          f"rank={overall_flags['b2_rank_ok']}, "
          f"chi2_coverage={overall_flags['b2_chi2_coverage_in_0p90_0p99']}")
    print(f"  b3 gaussian covered matched sigma="
          f"{overall_flags['b3_gaussian_all_covered_matched_sigma']}, "
          f"covered sigma1="
          f"{overall_flags['b3_gaussian_all_covered_sigma1']}; "
          f"heavy strong violations vs sigma1="
          f"{overall_flags['b3_heavy_strong_violation_vs_sigma1_count']}")
    print(f"  b3 sharpened 3E violations="
          f"{overall_flags['b3_sharpened_3e_violation_count']} "
          f"(expected: all reps)")
    print(f"  b4 pose_graph bound_cg_ok="
          f"{overall_flags['b4_pose_graph_bound_cg_ok']} "
          f"bundle bound_cg_ok={overall_flags['b4_bundle_bound_cg_ok']}")
    print(f"  b4 near-singular: bound_holds="
          f"{overall_flags['b4_near_singular_bound_holds']}, "
          f"small_residual_controls_fail="
          f"{overall_flags['b4_near_singular_small_residual_controls_fail']}, "
          f"identity_ok={overall_flags['b4_near_singular_identity_ok']}")
    print(f"  non-bool fail flags = {len(fails)}; "
          f"plots = {plots or 'skipped'}")
    audit = "PASS" if not fails else "FAIL"
    print(f"AUDIT: {audit}  (elapsed {time.time() - t_start:.1f}s)")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
