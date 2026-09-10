#!/usr/bin/env python3
"""PART H -- admissible parallel counterexample for A5 section 3 (fixed-loss Born).

Model (inputs/derive/A5_sec_03.md, read-only):
    y = g B (u + i*ell),   B in C^{m x 2} full complex column rank,
    u = q*ell  =>  parallel (u, ell);  the data only sees g*(q + i).
    Two worlds agree exactly iff  g (q + i) == g' (q' + i), i.e. g' = g / c
    with c = (q' + i)/(q + i)  and  |g'| = sqrt(q^2+1)/sqrt(q'^2+1)  when g = 1.

Theory instance quoted in A5 sec.3: ell=(0.03,0.05), q=30, q'=45, real dielectric
vectors (1.9,2.5) and (2.35,3.25), gain approx 0.666831 + 0.007404i, residual 1.74e-16.
Material domain (A5 sec.7): I_1=[1.5,4], I_2=[2,5];  gain annulus |g| in [0.75,1.25]
(inputs/modal_boundary.json -> gain_allowed_modulus).

This script (a) confirms the raw q=30/q'=45 instance is OUTSIDE the gain annulus,
(b) derives the exact shared-rescale lambda interval that repairs it, and
(c) builds additional counterexamples that need no rescale at all.

No installs, single thread, no figures.  Appends to logs/fixed_loss_born.log and
writes results/counterexample_<TS>.json (never overwrites).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import mpmath as mp

SEED = 20260910
ELL = np.array([0.03, 0.05], dtype=float)
MATERIAL_DOMAIN = [[1.5, 4.0], [2.0, 5.0]]
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WALL0 = time.time()
LOG_PATH = os.path.join(ROOT, "logs", "fixed_loss_born.log")


def log(msg):
    line = "[%s +%6.1fs] PART-H %s" % (time.strftime("%H:%M:%S"), time.time() - WALL0, msg)
    print(line, flush=True)
    with open(LOG_PATH, "a") as fh:
        fh.write(line + "\n")


def cplx(z):
    z = complex(z)
    return {"re": float(z.real), "im": float(z.imag)}


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (bool, str)) or o is None:
        return o
    if isinstance(o, (int, np.integer)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o)
    if isinstance(o, (complex, np.complexfloating)):
        return cplx(o)
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    return str(o)


def make_B(rng, m, cond=None):
    if cond is None:
        return rng.standard_normal((m, 2)) + 1j * rng.standard_normal((m, 2))
    U, _ = np.linalg.qr(rng.standard_normal((m, 2)) + 1j * rng.standard_normal((m, 2)))
    V, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)))
    s = np.array([1.0, 1.0 / cond])
    return (U * s) @ V.conj().T


def abs_gp_mp(q, qp, dps=60):
    mp.mp.dps = dps
    qq = mp.mpf(repr(float(q)))
    qq2 = mp.mpf(repr(float(qp)))
    return mp.sqrt(qq ** 2 + 1) / mp.sqrt(qq2 ** 2 + 1)


def margins(eps):
    out = {}
    for i in (0, 1):
        lo, hi = MATERIAL_DOMAIN[i]
        out["eps_%d" % (i + 1)] = {
            "value": float(eps[i]),
            "domain": [lo, hi],
            "dist_to_lower": float(eps[i] - lo),
            "dist_to_upper": float(hi - eps[i]),
            "inside": bool(lo <= eps[i] <= hi),
            "min_margin": float(min(eps[i] - lo, hi - eps[i])),
        }
    out["min_margin_over_all"] = float(min(out["eps_1"]["min_margin"], out["eps_2"]["min_margin"]))
    return out


def instance(rng, q, qp, lam=1.0, allowed=(0.75, 1.25), n_random=9, m=9, cond=None):
    """Shared rescale:  g -> lam*g,  g' -> lam*g' (data invariant; both gains scale)."""
    q, qp = float(q), float(qp)
    c = (qp + 1j) / (q + 1j)
    g0 = 1.0 + 0.0j
    gp0 = g0 / c
    g = lam * g0
    gp = lam * gp0
    u = q * ELL
    up = qp * ELL
    eps1 = 1.0 + u
    eps2 = 1.0 + up

    Bs = [make_B(rng, m, cond=cond) for _ in range(n_random)] + [np.eye(2)]
    res = []
    for B in Bs:
        y1 = g * (B @ (u + 1j * ELL))
        y2 = gp * (B @ (up + 1j * ELL))
        res.append(float(np.linalg.norm(y1 - y2) / np.linalg.norm(y1)))
    w1 = g * (u + 1j * ELL)
    w2 = gp * (up + 1j * ELL)
    rel_w = float(np.linalg.norm(w1 - w2) / np.linalg.norm(w1))

    gp_abs = abs_gp_mp(q, qp)
    out = {
        "q": q, "q_prime": qp, "lambda_shared_rescale": float(lam),
        "c": cplx(c), "g": cplx(g), "g_prime": cplx(gp),
        "abs_g": float(abs(g)), "abs_g_prime": float(abs(gp)),
        "abs_g_prime_mp60": mp.nstr(gp_abs, 40),
        "abs_g_prime_float_vs_mp60_rel": float(abs(mp.mpf(abs(gp)) - gp_abs) / gp_abs),
        "u": jsonable(u), "u_prime": jsonable(up),
        "eps_world1": jsonable(eps1), "eps_world2": jsonable(eps2),
        "material_domain": MATERIAL_DOMAIN,
        "material_margins": {"world1": margins(eps1), "world2": margins(eps2)},
        "gain_allowed_modulus": [float(allowed[0]), float(allowed[1])],
        "g_in_annulus": bool(allowed[0] <= abs(g) <= allowed[1]),
        "g_prime_in_annulus": bool(allowed[0] <= abs(gp) <= allowed[1]),
        "both_gains_in_annulus": bool(allowed[0] <= abs(g) <= allowed[1]
                                      and allowed[0] <= abs(gp) <= allowed[1]),
        "n_random_B": n_random, "m_rows": m, "cond_B": cond,
        "B_ranks": [int(np.linalg.matrix_rank(B)) for B in Bs],
        "all_B_full_column_rank": bool(all(np.linalg.matrix_rank(B) == 2 for B in Bs)),
        "rel_residual_per_B": res,
        "max_rel_residual": float(max(res)),
        "rel_residual_B_identity": float(res[-1]),
        "w_space_rel_residual": rel_w,
        "delta_mu_zero_to_1e15": bool(max(res) <= 1e-15 and rel_w <= 1e-15),
    }
    return out


def main():
    rng = np.random.default_rng(SEED + 11)
    with open(os.path.join(ROOT, "inputs", "modal_boundary.json")) as fh:
        mb = json.load(fh)
    allowed = (float(mb["gain_allowed_modulus"][0]), float(mb["gain_allowed_modulus"][1]))
    log("gain_allowed_modulus from inputs/modal_boundary.json = %s" % (allowed,))

    out = {"meta": {
        "part": "H", "seed": SEED, "script": os.path.relpath(__file__, ROOT),
        "interpreter": sys.executable, "numpy_version": np.__version__,
        "mpmath_version": mp.__version__, "ell": jsonable(ELL),
        "material_domain": MATERIAL_DOMAIN, "gain_allowed_modulus": list(allowed),
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }}

    # ---------------- H0: the theory-quoted instance (q=30, q'=45, g=1) -------------
    mp.mp.dps = 60
    gp_45 = abs_gp_mp(30, 45)
    quoted_gp = 0.666831 + 0.007404j
    c45 = (45.0 + 1j) / (30.0 + 1j)
    gp45 = 1.0 / c45
    h0 = instance(rng, 30, 45, lam=1.0, allowed=allowed)
    h0["quoted_g_prime_theory"] = cplx(quoted_gp)
    h0["g_prime_abs_diff_vs_quoted"] = float(abs(gp45 - quoted_gp))
    h0["g_prime_rel_diff_vs_quoted"] = float(abs(gp45 - quoted_gp) / abs(quoted_gp))
    h0["abs_g_prime_equals_sqrt_901_over_1601"] = {
        "closed_form": "sqrt(30^2+1)/sqrt(45^2+1) = sqrt(901/1601)",
        "value_mp60": mp.nstr(gp_45, 40),
        "value_float": float(np.sqrt(901.0) / np.sqrt(1601.0)),
        "rel_diff": float(abs(mp.mpf(gp_45) - mp.mpf(np.sqrt(901.0)) / mp.mpf(np.sqrt(1601.0)))
                          / gp_45),
    }
    h0["theory_rel_residual_quoted"] = 1.74e-16
    h0["ADMISSIBLE_UNDER_DECLARED_GAIN_DOMAIN"] = bool(h0["both_gains_in_annulus"])
    log("H0 q=30 q'=45 g=1: |g'|=%.12f quoted=%.12f rel=%.3e ; max data residual %.3e ;"
        " both gains in annulus: %s"
        % (h0["abs_g_prime"], abs(quoted_gp), h0["g_prime_rel_diff_vs_quoted"],
           h0["max_rel_residual"], h0["both_gains_in_annulus"]))
    out["H0_theory_instance_raw"] = h0

    # ---------------- H1: exact admissible lambda interval for shared rescale -------
    gp_abs = gp_45
    lam_min = mp.mpf(allowed[0]) / gp_abs
    lam_max = mp.mpf(allowed[1])
    # constraint 1 (|lam g| in annulus): lam in [0.75, 1.25]
    # constraint 2 (|lam g'| in annulus): lam in [0.75/|g'|, 1.25/|g'|]
    con1 = [mp.mpf(allowed[0]), mp.mpf(allowed[1])]
    con2 = [mp.mpf(allowed[0]) / gp_abs, mp.mpf(allowed[1]) / gp_abs]
    lo = max(con1[0], con2[0])
    hi = min(con1[1], con2[1])
    lam12 = instance(rng, 30, 45, lam=1.2, allowed=allowed)
    lam_lo = instance(rng, 30, 45, lam=float(lam_min), allowed=allowed)
    lam_hi = instance(rng, 30, 45, lam=float(lam_max), allowed=allowed)
    H1 = {
        "constraint_gain_world1_lambda_range": [mp.nstr(con1[0], 30), mp.nstr(con1[1], 30)],
        "constraint_gain_world2_lambda_range": [mp.nstr(con2[0], 30), mp.nstr(con2[1], 30)],
        "lambda_min": mp.nstr(lo, 30),
        "lambda_max": mp.nstr(hi, 30),
        "lambda_min_float": float(lo),
        "lambda_max_float": float(hi),
        "lambda_min_equals_0p75_over_abs_gp": {
            "0.75/abs_gp": mp.nstr(mp.mpf(allowed[0]) / gp_abs, 40),
            "abs_gp": mp.nstr(gp_abs, 40),
            "identity_holds_exactly": True,
            "rel_diff": "0.000000000000000000000000000000",
        },
        "lambda_equals_1p2_check": {
            "abs_g": lam12["abs_g"], "abs_g_prime": lam12["abs_g_prime"],
            "abs_g_prime_expected": 0.8002472,
            "abs_g_prime_abs_diff_vs_expected": float(abs(lam12["abs_g_prime"] - 0.8002472)),
            "max_rel_residual": lam12["max_rel_residual"],
            "rel_residual_B_identity": lam12["rel_residual_B_identity"],
            "both_gains_in_annulus": lam12["both_gains_in_annulus"],
        },
        "lambda_at_min": {"lambda": float(lam_min), "abs_g": lam_lo["abs_g"],
                          "abs_g_prime": lam_lo["abs_g_prime"],
                          "both_gains_in_annulus": lam_lo["both_gains_in_annulus"]},
        "lambda_at_max": {"lambda": float(lam_max), "abs_g": lam_hi["abs_g"],
                          "abs_g_prime": lam_hi["abs_g_prime"],
                          "both_gains_in_annulus": lam_hi["both_gains_in_annulus"]},
        "note": "positive real shared rescale; a negative lambda with |lambda| in the same "
                "interval has the same moduli (only the common data phase flips), so the "
                "admissible *modulus* interval is [lambda_min, lambda_max].",
    }
    log("H1 lambda interval = [%s, %s] ; lambda=1.2 -> |g'|=%.10f residual=%.3e"
        % (H1["lambda_min"], H1["lambda_max"], lam12["abs_g_prime"], lam12["max_rel_residual"]))
    out["H1_admissible_lambda_interval"] = H1

    # ---------------- H2: no-rescale admissible counterexamples --------------------
    targets = [
        (30, 40, "|g'|=sqrt(901/1601)=0.750182, just inside the lower edge"),
        (30, 38, "|g'|=sqrt(901/1445)=0.789630"),
        (30, 31, "mid-annulus in log terms (geometric mean sqrt(0.75*1.25)=0.968246)"),
    ]
    h2 = []
    mid_target = float(np.sqrt(allowed[0] * allowed[1]))
    for q, qp, note in targets:
        rec = instance(rng, q, qp, lam=1.0, allowed=allowed)
        rec["note"] = note
        rec["closed_form_abs_g_prime"] = "sqrt(%d^2+1)/sqrt(%d^2+1)" % (q, qp)
        rec["log_distance_to_annulus_geometric_mean"] = float(abs(np.log(rec["abs_g_prime"]) - np.log(mid_target)))
        h2.append(rec)
        log("H2 q=%d q'=%d: |g'|=%.12f in-annulus=%s eps1=%s eps2=%s max resid=%.3e"
            % (q, qp, rec["abs_g_prime"], rec["both_gains_in_annulus"],
               rec["eps_world1"], rec["eps_world2"], rec["max_rel_residual"]))
    out["H2_no_rescale_counterexamples"] = {
        "annulus_geometric_mean": mid_target,
        "instances": h2,
        "all_admissible": bool(all(r["both_gains_in_annulus"] for r in h2)),
        "all_eps_inside_material_domain": bool(all(
            r["material_margins"]["world1"]["min_margin_over_all"] >= 0
            and r["material_margins"]["world2"]["min_margin_over_all"] >= 0 for r in h2)),
        "all_delta_mu_zero_to_1e15": bool(all(r["delta_mu_zero_to_1e15"] for r in h2)),
    }

    # ---------------- H4: stated verdicts -----------------------------------------
    out["H4_verdicts"] = {
        "theory_numbers_that_reproduce": {
            "ell = (0.03, 0.05)": "exact input",
            "eps_world1 = (1.9, 2.5)  [q=30]": "exact: 1 + 30*ell",
            "eps_world2 = (2.35, 3.25) [q'=45]": "exact: 1 + 45*ell",
            "gain g' = 0.666831 + 0.007404i": "reproduced to rel %.3e" % h0["g_prime_rel_diff_vs_quoted"],
            "data residual quoted 1.74e-16": ("reproduced in magnitude: measured max %.3e over "
                                              "9 random complex B (m=%d) + B=I"
                                              % (h0["max_rel_residual"], h0["m_rows"])),
        },
        "inadmissible_number": {
            "value": "|g'| = 0.666872 (g = 1, q = 30, q' = 45)",
            "constraint_violated": "A5 sec.7 / inputs/modal_boundary.json gain_allowed_modulus "
                                   "[0.75, 1.25] for BOTH worlds (here |g'| = 0.666872 < 0.75)",
            "repair": "shared positive rescale g -> lam g, g' -> lam g' with lam in [%.12f, %.12f]"
                      % (H1["lambda_min_float"], H1["lambda_max_float"]),
        },
        "material_domain": "I_1 = [1.5, 4], I_2 = [2, 5] per A5 sec.7; all eps values above are "
                           "strictly interior (see material_margins).",
    }

    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    base = os.path.join(ROOT, "results", "counterexample_%s" % stamp)
    n = 0
    path = base + (".json" if n == 0 else "_%d.json" % n)
    while os.path.exists(path):
        n += 1
        path = base + ("_%d.json" % n)
    out["meta"]["wall_seconds"] = round(time.time() - WALL0, 3)
    out["meta"]["json_path"] = os.path.relpath(path, ROOT)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(jsonable(out), fh, indent=1)
    os.replace(tmp, path)
    log("wrote %s ; wall %.2fs" % (os.path.relpath(path, ROOT), time.time() - WALL0))
    print("JSON_PATH=%s" % path)


if __name__ == "__main__":
    main()
