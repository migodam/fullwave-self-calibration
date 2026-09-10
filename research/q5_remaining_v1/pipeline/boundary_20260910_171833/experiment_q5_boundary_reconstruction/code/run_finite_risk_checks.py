#!/usr/bin/env python3
"""PART B: independent reproduction of the A5 sections 5-7 finite-risk constants.

Everything here is recomputed from the exact Mie reactance q(eps, x)
(code/mie_reactance.py) and the supplied world values in inputs/modal_boundary.json.
Nothing is tuned to match the supplied constants; where a constant does not
reproduce, the exact discrepancy is reported.

Usage:
    PYTHONPATH=<uv archive with mpmath> <a3 venv python> code/run_finite_risk_checks.py
"""

from __future__ import annotations

import json
import os
import platform
import shlex
import sys
import time
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import mie_reactance as M               # noqa: E402
import run_qprime_variant_sweep as A    # noqa: E402  (Sink / Tee / dec30 / json_safe)

SEED = 20260910
# my own, more conservative interval lower bounds of h' (results/interval_20260910T093408Z.json)
MY_AMPLITUDE_DERIVATIVE_LOWER = [4.2929719577e-04, 1.2943650641e-04]


def h_of(eps, x):
    """|t(eps,x)| = q / sqrt(1 + q^2), numpy vectorised."""
    q = M.q_closed(eps, x)
    return q / np.sqrt(1.0 + q * q)


def make_h_inverse(x, eps_lo, eps_hi, ngrid=400001):
    """Monotone table inversion of h(.,x) on [eps_lo, eps_hi] + Newton refinement."""
    grid = np.linspace(eps_lo, eps_hi, ngrid)
    hgrid = h_of(grid, x)
    if not np.all(np.diff(hgrid) > 0):
        raise SystemExit("h is not strictly increasing on the material interval")

    def inv(target):
        e = np.interp(target, hgrid, grid)
        for _ in range(4):
            e = e - (h_of(e, x) - target) / M.hprime_closed(e, x)
        return e

    return inv, float(hgrid[0]), float(hgrid[-1])


def main():
    sys.stdout = A.Tee(os.path.join(ROOT, "logs", "finite_risk.log"))
    sys.stderr = sys.stdout
    t0 = time.time()

    from mpmath import mp, mpf, mpc, sqrt, log, erfc, erfinv, pi

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_path = os.path.join(ROOT, "results", f"finite_risk_{stamp}.json")
    sink = A.Sink(out_path)
    mp.dps = 50

    print("=" * 100)
    print("PART B  independent reproduction of A5 sections 5-7 finite-risk constants")
    print("started_utc :", datetime.now(timezone.utc).isoformat())
    print("interpreter :", sys.executable)
    print("platform    :", platform.platform())
    print("numpy       :", np.__version__)
    print("rng seed    :", SEED, "(numpy default_rng; child seeds SEED/SEED+1)")
    print("=" * 100)

    hashes = {p: A.sha256_file(os.path.join(ROOT, p)) for p in [
        "code/mie_reactance.py", "code/run_finite_risk_checks.py",
        "inputs/modal_boundary.json"]}
    for k, v in hashes.items():
        print(f"  sha256 {v}  {k}")

    J = json.load(open(os.path.join(ROOT, "inputs", "modal_boundary.json")))
    sink.put("meta", {
        "task": "PART B finite-risk constant checks",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command_line": " ".join(shlex.quote(a) for a in sys.argv),
        "interpreter": sys.executable, "numpy": np.__version__,
        "platform": platform.platform(), "cwd": os.getcwd(),
        "file_sha256": hashes, "rng_seed": SEED,
        "mpmath_dps_for_exact_steps": 50,
        "my_amplitude_derivative_lower": MY_AMPLITUDE_DERIVATIVE_LOWER,
    })

    eps1 = [mpf(repr(v)) for v in J["world1_eps"]]
    eps2 = [mpf(repr(v)) for v in J["world2_eps"]]
    x_i = [mpf(repr(J["monotonicity_checks"][i]["x"])) for i in range(2)]
    g1 = mpf(repr(J["gain1"]))
    g2 = mpf(repr(J["gain2"]))
    a_min = mpf(repr(J["gain_allowed_modulus"][0]))
    sigma = [mpf(repr(v)) for v in J["complex_noise_sigma"]]
    sigma_a = mpf(repr(J["reference_scalar_real_standard_deviation"]))
    cdps = lambda phi: 0.5 * erfc(-phi / sqrt(2))          # Phi
    b_a = sigma_a * sqrt(2) * erfinv(2 * mpf("0.9995") - 1)  # sigma_a * Phi^{-1}(0.9995)

    print(f"\n  Phi^{-1}(0.9995) = {mp.nstr(sqrt(2)*erfinv(2*mpf('0.9995')-1), 20)}")
    print(f"  b_a = sigma_a * Phi^-1(0.9995) = {mp.nstr(b_a, 20)}"
          f"   (supplied reference_error_bound = {J['reference_error_bound']!r})")
    print(f"  a_min = {mp.nstr(a_min, 20)}")

    # ---------------------------------------------------------------- step 1: scale c
    print("\n" + "-" * 100)
    print("STEP 1  shared scale c_i = q(world2_eps_i, x_i) / q(world1_eps_i, x_i)")
    q1 = [M.mp_q(eps1[i], x_i[i], 50) for i in range(2)]
    q2 = [M.mp_q(eps2[i], x_i[i], 50) for i in range(2)]
    c_i = [q2[i] / q1[i] for i in range(2)]
    step1 = {"q_world1": [A.dec30(v) for v in q1], "q_world2": [A.dec30(v) for v in q2]}
    for i in range(2):
        print(f"  ch{i+1} x={float(x_i[i])}: q1={mp.nstr(q1[i], 20)}  q2={mp.nstr(q2[i], 20)}"
              f"  c={mp.nstr(c_i[i], 22)}  c-1.25={mp.nstr(c_i[i]-mpf('1.25'), 4)}")
        print(f"      supplied world1_q={J['world1_q'][i]!r}  world2_q={J['world2_q'][i]!r}")
    print(f"  gain1={float(g1)} gain2={float(g2)}  gain1/c={mp.nstr((g1/(c_i[0]+c_i[1])*2), 20)}"
          f"  gain2-(gain1/c_1)={mp.nstr(g2-g1/c_i[0], 4)}")
    step1["c"] = [A.dec30(v) for v in c_i]
    step1["c_minus_1p25"] = [A.dec30(v - mpf("1.25")) for v in c_i]
    step1["gain2_equals_gain1_over_c_rel"] = [float(abs(g2 - g1 / v) / g2) for v in c_i]
    sink.put("step1_scale", step1)

    # ------------------------------------------------- step 2: section 5 identity
    print("\n" + "-" * 100)
    print("STEP 2  section 5 channel-difference identity, two ways")
    I = mpc(0, 1)
    t_of_q = lambda q: I * q / (1 - I * q)
    step2 = {}
    for i in range(2):
        c = c_i[i]
        d_direct = g2 * t_of_q(c * q1[i]) - g1 * t_of_q(q1[i])
        d_closed = -(g1 * (c - 1) * q1[i] ** 2) / ((1 - I * c * q1[i]) * (1 - I * q1[i]))
        # Cross-check the reactance against the Riccati-Bessel a1 route.
        # Module convention (see logs/mie_checks.log block [C]): the standard a1 uses
        # xi = psi - i*chi, i.e. xi_sign=-1, and then A5's t = -a1 exactly.  The other
        # branch xi_sign=+1 returns -conj(t), so q_mie(...,+1) is a *different*
        # (conjugated-convention) quantity, q/(1+2iq), not q.
        e1f, x1f = float(eps1[i]), float(x_i[i])
        a1_std = complex(M.mie_a1(e1f, x1f, -1))      # standard a1, xi = psi - i*chi
        a1_plus = complex(M.mie_a1(e1f, x1f, +1))     # other branch, xi = psi + i*chi
        t_alt = complex(t_of_q(q1[i]))
        q_num = complex(q1[i])
        q_route_a1 = 1j * a1_std / (1.0 - a1_std)     # correct q from the standard a1
        qm_std = complex(M.q_mie(e1f, x1f, -1))       # module wrapper, prompt-literal conversion
        qm_plus = complex(M.q_mie(e1f, x1f, +1))
        step2[f"ch{i+1}"] = {
            "c": A.dec30(c),
            "delta_mu_direct_re": A.dec30(mp.re(d_direct)),
            "delta_mu_direct_im": A.dec30(mp.im(d_direct)),
            "abs_direct": A.dec30(abs(d_direct)),
            "abs_closed": A.dec30(abs(d_closed)),
            "abs_disc": A.dec30(abs(d_direct - d_closed)),
            "rel_disc": float(abs(d_direct - d_closed) / abs(d_direct)),
            "abs_rel_formula": A.dec30(abs(g2 * t_of_q(c * q1[i]) - g1 * t_of_q(q1[i])) / abs(g1 * t_of_q(q1[i]))),
            "bound_c_minus_1_times_q": A.dec30(abs(c - 1) * abs(q1[i])),
            "t_from_q_vs_neg_a1_xi_minus_abs": repr(abs(t_alt + a1_std)),
            "q_from_a1_route_vs_q_rel": repr(abs(q_route_a1 - q_num) / abs(q_num)),
            "q_mie_xi_minus_vs_q_rel": repr(abs(qm_std - q_num) / abs(q_num)),
            "a1_xi_plus_vs_neg_conj_t_abs": repr(abs(a1_plus + t_alt.conjugate())),
            "q_mie_xi_plus_vs_q_rel": repr(abs(qm_plus - q_num) / abs(q_num)),
            "note_q_mie_wrapper": ("code/mie_reactance.q_mie uses the prompt-literal "
                                   "conversion -i*a1/(1+a1); it is NOT equivalent to q "
                                   "for either xi sign (see logs/mie_checks.log block [C])"),
        }
        s = step2[f"ch{i+1}"]
        print(f"  ch{i+1}: |delta_mu| direct={s['abs_direct']}  closed={s['abs_closed']}"
              f"  abs disc={s['abs_disc']}  rel disc={s['rel_disc']:.3e}")
        print(f"      delta_mu = {s['delta_mu_direct_re']} + i*{s['delta_mu_direct_im']}")
        print(f"      |delta|/|g t| = {s['abs_rel_formula']}  <= |c-1|q = {s['bound_c_minus_1_times_q']}")
        print(f"      |t - (-a1(xi=-1))| = {s['t_from_q_vs_neg_a1_xi_minus_abs']}"
              f"   |(i*a1/(1-a1))(xi=-1) - q|/|q| = {s['q_from_a1_route_vs_q_rel']}")
        print(f"      |a1(xi=+1) - (-conj t)| = {s['a1_xi_plus_vs_neg_conj_t_abs']}"
              f"   wrapper q_mie offsets: xi=-1 {s['q_mie_xi_minus_vs_q_rel']},"
              f" xi=+1 {s['q_mie_xi_plus_vs_q_rel']}")
    sink.put("step2_identity", step2)
    delta = [g2 * t_of_q(c_i[i] * q1[i]) - g1 * t_of_q(q1[i]) for i in range(2)]

    # ------------------------------------------------- step 3: sigma = 0.2% of |t1|
    print("\n" + "-" * 100)
    print("STEP 3  sigma_i = 0.002 * |g1 * t(eps1_i, x_i)| ?")
    step3 = {}
    for i in range(2):
        amp = abs(g1 * t_of_q(q1[i]))
        factor = sigma[i] / amp
        step3[f"ch{i+1}"] = {"abs_g1_t1": A.dec30(amp), "sigma_supplied": A.dec30(sigma[i]),
                             "implied_factor": A.dec30(factor),
                             "rel_diff_vs_0p002": float(abs(factor - mpf("0.002")) / mpf("0.002"))}
        print(f"  ch{i+1}: |g1 t1| = {mp.nstr(amp, 20)}  sigma = {mp.nstr(sigma[i], 20)}"
              f"  implied factor = {mp.nstr(factor, 20)}  (rel vs 0.002: {step3[f'ch{i+1}']['rel_diff_vs_0p002']:.3e})")
    sink.put("step3_sigma_factor", step3)

    # ------------------------------------------- step 4: D, D^2 and p_*
    print("\n" + "-" * 100)
    print("STEP 4  whitened distance D and equal-prior optimal error p_*")
    D2 = sum(abs(delta[i]) ** 2 / sigma[i] ** 2 for i in range(2))
    D = sqrt(D2)
    p_star = 0.5 * erfc(D / 2)
    sup_D = mpf(repr(J["complex_whitened_distance"]))
    sup_p = mpf(repr(J["no_reference_equal_prior_optimal_error"]))
    step4 = {"D2": A.dec30(D2), "D": A.dec30(D),
             "D_supplied": repr(J["complex_whitened_distance"]),
             "D_abs_diff": A.dec30(D - sup_D),
             "D_rel_diff": float(abs(D - sup_D) / sup_D),
             "p_star": A.dec30(p_star), "p_star_supplied": repr(J["no_reference_equal_prior_optimal_error"]),
             "p_star_abs_diff": A.dec30(p_star - sup_p),
             "p_star_rel_diff": float(abs(p_star - sup_p) / sup_p)}
    for k, v in step4.items():
        print(f"  {k:22s} = {v}")
    sink.put("step4_distance", step4)

    # --------------------------------------- step 5: Monte-Carlo theorem 4
    print("\n" + "-" * 100)
    print("STEP 5  independent Monte-Carlo check of Theorem 4 (equal-prior optimal test)")
    N = 4_000_000
    mu0 = np.array([complex(g1 * t_of_q(q1[i])) for i in range(2)], dtype=np.complex128)
    mu1 = np.array([complex(g2 * t_of_q(c_i[i] * q1[i])) for i in range(2)], dtype=np.complex128)
    dl = mu1 - mu0
    sig = np.array([float(sigma[i]) for i in range(2)])
    inv_s2 = 1.0 / sig ** 2
    tau = float(sum(abs(dl[i]) ** 2 * inv_s2[i] for i in range(2)))
    rng = np.random.default_rng(SEED)
    chunk = 500_000
    err_cnt = [0, 0]
    means = [np.zeros(2, dtype=np.complex128), np.zeros(2, dtype=np.complex128)]
    for j, mu in ((0, mu0), (1, mu1)):
        done = 0
        while done < N:
            m = min(chunk, N - done)
            n = (rng.standard_normal((m, 2)) + 1j * rng.standard_normal((m, 2))) * (sig / np.sqrt(2.0))
            y = mu + n
            L = 2.0 * np.real(((y - mu0) * np.conj(dl) * inv_s2).sum(axis=1))
            dec = L > tau
            err_cnt[j] += int(np.count_nonzero(dec != (j == 1)))
            means[j] += y.sum(axis=0)
            done += m
    means = [mm / N for mm in means]
    emp = [(err_cnt[0] + err_cnt[1]) / (2.0 * N), err_cnt[0] / N, err_cnt[1] / N]
    ci = 1.96 * np.sqrt(emp[0] * (1 - emp[0]) / (2.0 * N))
    D_hat = float(np.sqrt(sum(abs(means[1][i] - means[0][i]) ** 2 * inv_s2[i] for i in range(2))))
    step5 = {"N_per_world": N, "tau_D2": repr(tau), "emp_mean_error": emp[0],
             "emp_error_world0": emp[1], "emp_error_world1": emp[2],
             "emp_err_counts": err_cnt, "binomial_ci95_halfwidth": ci,
             "ci95": [emp[0] - ci, emp[0] + ci],
             "theoretical_Phi_minus_D_over_sqrt2": float(p_star),
             "theoretical_in_ci": bool(emp[0] - ci <= float(p_star) <= emp[0] + ci),
             "empirical_whitened_distance_from_sample_means": D_hat,
             "D_rel_diff_from_sample_means": float(abs(D_hat - float(D)) / float(D))}
    for k, v in step5.items():
        print(f"  {k:44s} = {v}")
    sink.put("step5_monte_carlo_theorem4", step5)

    # ------------------------------------------- step 6: section 7 error bounds
    print("\n" + "-" * 100)
    print("STEP 6  section 7 uniform material error bounds")
    b_i = [sigma[i] * sqrt(log(mpf(1000))) for i in range(2)]
    H_i = [mpf(repr(J["monotonicity_checks"][i]["q_upper"])) for i in range(2)]
    m_i_sup = [mpf(repr(J["monotonicity_checks"][i]["amplitude_derivative_lower"])) for i in range(2)]
    m_i_my = [mpf(repr(v)) for v in MY_AMPLITUDE_DERIVATIVE_LOWER]

    def bound(m, H, b):
        return (b + H * b_a) / ((a_min - b_a) * m)

    bounds_sup = [bound(m_i_sup[i], H_i[i], b_i[i]) for i in range(2)]
    bounds_my = [bound(m_i_my[i], H_i[i], b_i[i]) for i in range(2)]
    sup_bounds = [mpf(repr(v)) for v in J["uniform_eps_error_bounds"]]
    fail_p = 2 * (mpf(1) / 1000) + 2 * (mpf(1) - cdps(mpf("3.2905267314919255")))
    step6 = {"b_i": [A.dec30(v) for v in b_i], "b_a": A.dec30(b_a), "a_min": A.dec30(a_min),
             "H_i_from_supplied_q_upper": [A.dec30(v) for v in H_i],
             "m_i_supplied": [A.dec30(v) for v in m_i_sup],
             "m_i_mine": [A.dec30(v) for v in m_i_my],
             "bounds_with_supplied_m": [A.dec30(v) for v in bounds_sup],
             "bounds_supplied": [repr(v) for v in J["uniform_eps_error_bounds"]],
             "bounds_with_supplied_m_rel_diff_vs_supplied":
                 [float(abs(bounds_sup[i] - sup_bounds[i]) / sup_bounds[i]) for i in range(2)],
             "bounds_with_my_m": [A.dec30(v) for v in bounds_my],
             "bounds_growth_factor_my_over_supplied_m":
                 [float(bounds_my[i] / bounds_sup[i]) for i in range(2)],
             "bounds_with_my_m_rel_growth_pct":
                 [100.0 * float(bounds_my[i] / bounds_sup[i] - 1) for i in range(2)],
             "bounds_with_my_m_below_0p10": [bool(v < mpf("0.1")) for v in bounds_my],
             "per_channel_failure_prob": "exp(-b_i^2/sigma_i^2) = 1e-3 each",
             "reference_failure_prob": "2*(1-Phi(3.2905267315)) = 1e-3",
             "union_bound_failure_prob": A.dec30(fail_p),
             "implied_coverage_lower": A.dec30(1 - fail_p),
             "supplied_coverage_lower": repr(J["simultaneous_coverage_probability_lower"])}
    for k, v in step6.items():
        print(f"  {k:46s} = {v}")
    sink.put("step6_section7_bounds", step6)

    # ------------------------------------- step 7: Monte-Carlo section 7 estimator
    print("\n" + "-" * 100)
    print("STEP 7  independent Monte-Carlo check of the section 7 estimator claim")
    NM = 2_000_000
    rng7 = np.random.default_rng(SEED + 1)
    grid = [(1.5, 4.0), (2.0, 5.0)]
    inv = []
    for i in range(2):
        f, hlo, hhi = make_h_inverse(float(x_i[i]), grid[i][0], grid[i][1])
        inv.append((f, hlo, hhi))
    eps_pt = [eps1, eps2]
    gains = [float(g1), float(g2)]
    step7 = {"N_per_world": NM, "seed": SEED + 1, "channels": {}}
    for i in range(2):
        f, hlo, hhi = inv[i]
        bi = float(b_i[i])
        si = float(sigma[i])
        ch = {"x": float(x_i[i]), "material_interval": list(grid[i]),
              "clip_range_h": [hlo, hhi], "bound": float(bounds_my[i]),
              "b_i": bi, "b_a": float(b_a), "sigma_i": si,
              "worlds": {}}
        for j in range(2):
            g = gains[j]
            e = float(eps_pt[j][i])
            t = complex(g * t_of_q(M.mp_q(eps_pt[j][i], x_i[i], 50)))
            n_good = 0
            n_viol = 0
            max_err_good = 0.0
            max_viol = 0.0
            max_err_all = 0.0
            done = 0
            while done < NM:
                m = min(chunk, NM - done)
                xi = rng7.standard_normal(m) * float(sigma_a)
                eta = (rng7.standard_normal(m) + 1j * rng7.standard_normal(m)) * (si / np.sqrt(2.0))
                yy = t + eta
                zz = g + xi
                ratio = np.abs(yy) / zz
                ep = f(np.clip(ratio, hlo, hhi))
                err = np.abs(ep - e)
                good = (np.abs(eta) <= bi) & (np.abs(xi) <= float(b_a))
                n_good += int(np.count_nonzero(good))
                eg = err[good]
                if eg.size:
                    max_err_good = max(max_err_good, float(eg.max()))
                    bad = eg[eg > float(bounds_my[i])]
                    if bad.size:
                        n_viol += int(bad.size)
                        max_viol = max(max_viol, float(bad.max()))
                max_err_all = max(max_err_all, float(err.max()))
                done += m
            gf = n_good / NM
            gci = 1.96 * float(np.sqrt(gf * (1.0 - gf) / NM))
            ch["worlds"][f"world{j}"] = {
                "gain": g, "eps": e, "t": [t.real, t.imag],
                "good_fraction": gf,
                "good_fraction_ci95": [gf - gci, gf + gci],
                "per_channel_theory_0p998001_in_ci": bool(gf - gci <= 0.998001 <= gf + gci),
                "joint_theory_0p997_below_good_fraction": bool(gf > 0.997),
                # one eta for this channel + one shared xi -> (1-1e-3)^2 = 0.998001;
                # the supplied >=0.997 figure is the joint bound over both channels
                # and the reference (three events, union bound 3e-3).
                "theoretical_good_lower_per_channel": 0.998001,
                "supplied_joint_good_lower": 0.997,
                "max_err_on_good_event": max_err_good,
                "bound": float(bounds_my[i]),
                "bound_satisfied_on_good_event": bool(max_err_good <= float(bounds_my[i])),
                "violations_count": n_viol, "max_violation": max_viol,
                "unconditional_max_err": max_err_all}
            w = ch["worlds"][f"world{j}"]
            print(f"  ch{i+1} world{j} (eps={e:.6f}, g={g}): good frac={w['good_fraction']:.6f}"
                  f"  max|err| on good={w['max_err_on_good_event']:.6e}  bound={w['bound']:.6e}"
                  f"  ok={w['bound_satisfied_on_good_event']}  viol={n_viol} (max {max_viol:.3e})"
                  f"  uncond max err={max_err_all:.6e}")
        step7["channels"][f"ch{i+1}"] = ch
    sink.put("step7_monte_carlo_estimator", step7)

    # ------------------------------------------------------------------ summary
    print("\n" + "=" * 100)
    print("SUMMARY")
    summary = {
        "c_i": [float(v) for v in c_i], "c_i_minus_1.25": [float(v - mpf("1.25")) for v in c_i],
        "abs_delta_mu": [float(abs(d)) for d in delta],
        "delta_mu_closed_form_rel_disc": [step2[f"ch{i+1}"]["rel_disc"] for i in range(2)],
        "sigma_implied_factor": [float(sigma[i] / abs(g1 * t_of_q(q1[i]))) for i in range(2)],
        "D": float(D), "D_supplied": J["complex_whitened_distance"],
        "D_rel_diff": float(abs(D - sup_D) / sup_D),
        "p_star": float(p_star), "p_star_supplied": J["no_reference_equal_prior_optimal_error"],
        "p_star_rel_diff": float(abs(p_star - sup_p) / sup_p),
        "mc_mean_error": emp[0], "mc_ci95": [emp[0] - ci, emp[0] + ci],
        "bounds_with_supplied_m": [float(v) for v in bounds_sup],
        "bounds_with_my_m": [float(v) for v in bounds_my],
        "bounds_supplied": J["uniform_eps_error_bounds"],
        "good_fractions": [[step7["channels"][f"ch{i+1}"]["worlds"][f"world{j}"]["good_fraction"]
                            for j in range(2)] for i in range(2)],
        "violations_total": sum(step7["channels"][f"ch{i+1}"]["worlds"][f"world{j}"]["violations_count"]
                                for i in range(2) for j in range(2)),
        "wall_seconds": time.time() - t0,
    }
    sink.put("summary", summary)
    for k, v in summary.items():
        print(f"  {k} = {v}")
    print("wrote", os.path.relpath(out_path, ROOT))
    print("=" * 100)


if __name__ == "__main__":
    main()
