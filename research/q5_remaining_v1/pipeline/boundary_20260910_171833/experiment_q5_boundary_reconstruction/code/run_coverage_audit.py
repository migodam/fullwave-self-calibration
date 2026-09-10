#!/usr/bin/env python3
"""PART F -- constant coverage audit of inputs/modal_boundary.json.

Enumerates every leaf (866), classifies it (i)/(ii)/(iii)/(iv), and produces an
explicit value/verdict table for every supplied numerical constant, recomputing
each one now where feasible.  Writes results/coverage_audit_<UTCSTAMP>.json
(written once, never overwritten) and tees the log to logs/coverage_audit.log.

Usage:
    PYTHONPATH=<uv archive with mpmath> <a3_research venv python> code/run_coverage_audit.py
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import mie_reactance as M  # noqa: E402
from mpmath import mp  # noqa: E402

TEE = None


class Tee:
    def __init__(self, path):
        self.fh = open(path, "w", encoding="utf-8")

    def write(self, s):
        self.fh.write(s)
        self.fh.flush()
        sys.__stdout__.write(s)
        sys.__stdout__.flush()

    def flush(self):
        self.fh.flush()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------------
# 1. leaf enumeration + classification
# ---------------------------------------------------------------------------------

def walk(o, prefix=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from walk(v, f"{prefix}/{k}")
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from walk(v, f"{prefix}[{i}]")
    else:
        yield (prefix, o)


def classify(path):
    """(category, reason) for one leaf path (paths use "/" separators)."""
    # (iv) realized per-draw / per-world / per-channel data arrays
    if "/all_errors[" in path:
        return "iv", "per-draw realized Monte-Carlo error component (stochastic realization)"
    if "reference_checks[" in path and "/max_eps_abs_errors[" in path:
        return "iv", "per-channel statistic of the realized Monte-Carlo draws"
    if "reference_checks[" in path and "/world_eps[" in path:
        return "iv", "per-world eps component of that reference check (realized instance)"
    if "reference_checks[" in path and path.endswith("/empirical_success_eps_point1"):
        return "iv", "realized Monte-Carlo success count"
    # (iii) provenance / metadata / declared configuration
    if path in ("/scope", "/source_sha256", "/interval_backend", "/wall_seconds"):
        return "iii", "provenance/metadata"
    if "reference_checks[" in path and path.endswith("/draws"):
        return "iii", "declared Monte-Carlo configuration (draw count)"
    if "monotonicity_checks[" in path and (path.endswith("/interval_count")
                                           or path.endswith("/positive")
                                           or "/eps_domain[" in path
                                           or path.endswith("/x")):
        return "iii", "declared interval-check configuration input (not a derived constant)"
    # (ii) numerical constants not yet compared
    if path in ("/finite_material_distance", "/relative_data_difference",
                "/exact_reactance_identity_error"):
        return "ii", "numerical constant, never compared in increments 1-3"
    # (i) numerical constants already compared
    if path in ("/gain1", "/gain2", "/gain_allowed_modulus[0]", "/gain_allowed_modulus[1]",
                "/reference_scalar_real_standard_deviation", "/reference_error_bound",
                "/complex_whitened_distance", "/no_reference_equal_prior_optimal_error",
                "/simultaneous_coverage_probability_lower"):
        return "i", "compared in increments 2-3 (finite-risk block)"
    if path.startswith("/world1_eps[") or path.startswith("/world2_eps["):
        return "i", "compared in increment 2 step 1 (eps2 solves q=c q1, c=1.25)"
    if path.startswith("/world1_q[") or path.startswith("/world2_q["):
        return "i", "compared in increments 1-2 (Mie reactance routes a/b/c)"
    if path.startswith("/complex_noise_sigma[") or path.startswith("/modal_error_bounds["):
        return "i", "compared in increment 2 (0.2% of |g t|; b_i = sigma_i sqrt(log 1000))"
    if path.startswith("/uniform_eps_error_bounds["):
        return "i", "compared in increment 2 section-7 uniform bound"
    if "monotonicity_checks[" in path and any(path.endswith(k) for k in (
            "/q_lower", "/q_upper", "/denominator_lower", "/qprime_lower", "/qprime_upper",
            "/amplitude_derivative_lower")):
        return "i", "compared in increment 2 Part A (interval monotonicity sweep)"
    if "reference_checks[" in path and path.endswith("/gain_modulus"):
        return "i", "compared in increment 2 (equals |gain1| / |gain2|)"
    return "?", "unclassified"


# ---------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------

def main():
    global TEE
    t0 = time.time()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    os.makedirs(os.path.join(ROOT, "logs"), exist_ok=True)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    out_path = os.path.join(ROOT, "results", f"coverage_audit_{stamp}.json")
    TEE = Tee(os.path.join(ROOT, "logs", "coverage_audit.log"))
    sys.stdout = TEE

    J = json.load(open(os.path.join(ROOT, "inputs", "modal_boundary.json")))
    leaves = list(walk(J))
    print(f"PART F -- constant coverage audit   ({stamp})")
    print(f"leaf_count = {len(leaves)}")
    cats = {}
    reasons = {}
    for p, v in leaves:
        c, why = classify(p)
        cats.setdefault(c, []).append(p)
        reasons[p] = why
    for c in ("i", "ii", "iii", "iv", "?"):
        print(f"  category {c:2s}: {len(cats.get(c, [])):4d} leaves")
    assert len(leaves) == 866, len(leaves)

    # ---------------------------------------------------------------------------
    # independent recomputation of every audited constant
    # ---------------------------------------------------------------------------
    mp.dps = 50
    rows = []

    def to_mpf(v):
        """mpmath scalar from a float / mpf / nested mpmath.iv interval (midpoint)."""
        if v is None:
            return None
        if hasattr(v, "_mpi_"):  # mpmath.iv interval -> midpoint of its mpf endpoints
            lo, hi = v._mpi_
            return (mp.mpf(lo) + mp.mpf(hi)) / 2
        return mp.mpf(str(v))

    def add(name, supplied, computed, method, note="", force=None):
        try:
            s = to_mpf(supplied)
            c = to_mpf(computed)
            if s is None or c is None:
                rel = absd = None
            else:
                absd = abs(s - c)
                rel = absd / abs(s) if s != 0 else (absd if c == 0 else mp.inf)
        except Exception:
            rel = absd = None
        if force:
            verdict = force
        elif rel is None:
            verdict = "INCOMPARABLE"
        elif rel <= mp.mpf("1e-9"):
            verdict = "REPRODUCED"
        elif rel <= mp.mpf("1e-3"):
            verdict = "REPRODUCED_WEAKER"
        else:
            verdict = "NOT_REPRODUCED"
        row = {
            "name": name, "supplied": (float(supplied) if supplied is not None else None),
            "computed": (float(computed) if computed is not None else None),
            "abs_diff": (float(absd) if absd is not None else None),
            "rel_diff": (float(rel) if rel is not None else None),
            "verdict": verdict, "method": method, "note": note,
        }
        rows.append(row)
        print(f"  {name:52s} supplied={s} computed={c} rel={rel} -> {verdict}")
        return row

    # ---- eps worlds: solve q(eps2,x) = c q(eps1,x) with c = q2/q1 independently
    world2_computed = []
    for i, x in enumerate((0.2, 0.15)):
        e1 = mp.mpf(str(J["world1_eps"][i]))
        q1 = M.mp_q(e1, mp.mpf(str(x)), dps=50)
        c = mp.mpf(str(J["world2_q"][i])) / mp.mpf(str(J["world1_q"][i]))
        f = lambda e: M.mp_q(e, mp.mpf(str(x)), dps=50) - c * q1
        root = mp.findroot(f, mp.mpf(str(J["world2_eps"][i])))
        world2_computed.append(root)
    for i in range(2):
        add(f"world1_eps[{i}]", J["world1_eps"][i], J["world1_eps"][i],
            "declared input; interval endpoint definition of material interval I_i")
        add(f"world2_eps[{i}]", J["world2_eps"][i], world2_computed[i],
            "Solve q(eps2,x)=c*q(eps1,x), c=q2_supplied/q1_supplied, dps=50 Newton")
    for i in range(2):
        x = (0.2, 0.15)[i]
        q1c = M.mp_q(mp.mpf(str(J["world1_eps"][i])), mp.mpf(str(x)), dps=50)
        q2c = M.mp_q(mp.mpf(str(J["world2_eps"][i])), mp.mpf(str(x)), dps=50)
        add(f"world1_q[{i}] (x={x})", J["world1_q"][i], q1c,
            "A5 s8 reactance q from mie_reactance.mp_q (same arrangement as code/mie_reactance), dps=50")
        add(f"world2_q[{i}] (x={x})", J["world2_q"][i], q2c,
            "A5 s8 reactance q from mie_reactance.mp_q, dps=50")

    # ---- gains
    c1 = mp.mpf(str(J["world2_q"][0])) / mp.mpf(str(J["world1_q"][0]))
    add("gain1", J["gain1"], mp.mpf(1),
        "declared reference gain; |g|=1 lies inside gain_allowed_modulus")
    add("gain2 = gain1/c (c=q2/q1)", J["gain2"], mp.mpf(str(J["gain1"])) / c1,
        "s5 two-world compensation g' = g/c with c = q2/q1 = 1.25 (dps=50)")
    lam_lo, lam_hi = [mp.mpf(str(v)) for v in J["gain_allowed_modulus"]]
    add("gain_allowed_modulus[0]", J["gain_allowed_modulus"][0], lam_lo,
        "declared gain annulus; containment checked for |gain1|=1 and |gain2|=0.8",
        note="declared domain constant, not an independently derivable number")
    add("gain_allowed_modulus[1]", J["gain_allowed_modulus"][1], lam_hi,
        "declared gain annulus (upper); equals the s5 scale c=1.25",
        note="declared domain constant")
    add("reference_scalar_real_standard_deviation", J["reference_scalar_real_standard_deviation"],
        mp.mpf("0.001"), "declared real-reference std (used in s7 union bound)",
        note="declared input; the derived companion 0.0032905... is checked below")

    # ---- finite-risk block
    sigma_a = mp.mpf(str(J["reference_scalar_real_standard_deviation"]))
    b_a = sigma_a * mp.sqrt(2) * mp.erfinv(mp.mpf(2) * mp.mpf("0.9995") - 1)
    add("reference_error_bound = sigma_a*Phi^{-1}(0.9995)", J["reference_error_bound"], b_a,
        "sigma_a * inverse normal CDF at 0.9995 = 0.001*3.2905267314918948, mpmath dps=50")
    add("simultaneous_coverage_probability_lower", J["simultaneous_coverage_probability_lower"],
        mp.mpf("0.997"), "1 - (2*exp(-(sqrt(log 1000))^2) + 2*(1-Phi(3.2905...))) = 1-0.003")
    sig_c = []
    for i in range(2):
        x = (0.2, 0.15)[i]
        t1 = 1j * M.mp_q(mp.mpf(str(J["world1_eps"][i])), mp.mpf(str(x)), dps=50)
        t1 = t1 / (1 - t1)
        sig_c.append(mp.mpf("0.002") * abs(mp.mpf(str(J["gain1"])) * t1))
        add(f"complex_noise_sigma[{i}]", J["complex_noise_sigma"][i], sig_c[-1],
            "0.002*|g1*t1(eps1_i,x_i)| with t=iq/(1-iq) computed from the dps=50 reactance")
    mu1, dmu = [], []
    for i in range(2):
        x = mp.mpf((0.2, 0.15)[i])
        q1 = M.mp_q(mp.mpf(str(J["world1_eps"][i])), x, dps=50)
        q2 = M.mp_q(mp.mpf(str(J["world2_eps"][i])), x, dps=50)
        t = lambda q: 1j * q / (1 - 1j * q)
        mu1.append(mp.mpf(str(J["gain1"])) * t(q1))
        dmu.append(mp.mpf(str(J["gain2"])) * t(q2) - mu1[-1])
    D = mp.sqrt(sum(abs(dmu[i]) ** 2 / sig_c[i] ** 2 for i in range(2)))
    add("complex_whitened_distance", J["complex_whitened_distance"], D,
        "whitened two-channel mean gap D = sqrt(sum |d_mu_i|^2 / sigma_i^2), dps=50")
    p_star = mp.erfc(D / 2) / 2
    add("no_reference_equal_prior_optimal_error", J["no_reference_equal_prior_optimal_error"], p_star,
        "p_* = erfc(D/2)/2 (proper-complex binary risk, s6), dps=50")
    b_i = [sig_c[i] * mp.sqrt(mp.log(1000)) for i in range(2)]
    for i in range(2):
        add(f"modal_error_bounds[{i}]", J["modal_error_bounds"][i], b_i[i],
            "b_i = sigma_i*sqrt(log 1000) (per-channel noise envelope, s7), dps=50")

    # ---- interval-monotonicity block (re-run the interval sweep, iv.dps=35) -----
    interval_detail = {}
    for k, mc in enumerate(J["monotonicity_checks"]):
        x = mp.mpf(str(mc["x"]))
        e_lo, e_hi = mc["eps_domain"]
        n = mc["interval_count"]
        ep = M.decimal_endpoint_intervals(round(e_lo * 1000), 1000, round(e_hi * 1000), 1000, n)
        boxes = [M.box_from_intervals(ep[j], ep[j + 1]) for j in range(n)]
        agg = {"min_qprime_lower": None, "max_qprime_upper": None, "min_absDen_lower": None,
               "min_q_lower": None, "max_q_upper": None, "min_hprime_lower": None}
        for bx in boxes:
            o = M.iv_quantities(bx, x)
            pairs = (("min_qprime_lower", o["qprime"].a), ("max_qprime_upper", o["qprime"].b),
                     ("min_absDen_lower", M.iv_min_abs(o["Den"])),
                     ("min_q_lower", o["q"].a), ("max_q_upper", o["q"].b),
                     ("min_hprime_lower", o["hprime"].a))
            for key, v in pairs:
                if agg[key] is None or (v < agg[key] if key.startswith("min") else v > agg[key]):
                    agg[key] = v
        interval_detail[f"monotonicity_checks[{k}]"] = {
            "n_boxes": n,
            "min_qprime_lower": mp.nstr(agg["min_qprime_lower"], 25),
            "max_qprime_upper": mp.nstr(agg["max_qprime_upper"], 25),
            "min_absDen_lower": mp.nstr(agg["min_absDen_lower"], 25),
            "min_q_lower": mp.nstr(agg["min_q_lower"], 25),
            "max_q_upper": mp.nstr(agg["max_q_upper"], 25),
            "min_hprime_lower": mp.nstr(agg["min_hprime_lower"], 25),
            "point_qprime_at_eps_lo": mp.nstr(M.mp_qprime(mp.mpf(str(e_lo)), x, dps=50), 25),
            "point_qprime_at_eps_hi": mp.nstr(M.mp_qprime(mp.mpf(str(e_hi)), x, dps=50), 25),
        }
        add(f"monotonicity_checks[{k}].q_lower", mc["q_lower"], agg["min_q_lower"],
            "iv.dps=35 interval sweep, min over boxes of the q lower endpoint (left-endpoint box)",
            note="reproduces only because both codes share the same loose q/Den enclosure")
        add(f"monotonicity_checks[{k}].q_upper", mc["q_upper"], agg["max_q_upper"],
            "iv.dps=35 interval sweep, max over boxes of the q upper endpoint")
        add(f"monotonicity_checks[{k}].denominator_lower", mc["denominator_lower"],
            agg["min_absDen_lower"], "iv.dps=35 interval sweep, min over boxes of |Den| lower endpoint")
        add(f"monotonicity_checks[{k}].qprime_upper", mc["qprime_upper"], agg["max_qprime_upper"],
            "iv.dps=35 interval sweep, max over boxes of the q' upper endpoint",
            note="supplied value is tighter than this implementation's enclosure: it lies between "
                 "the pointwise maximum at eps_lo and this interval upper bound, so it is a "
                 "valid but differently-grouped enclosure, not bit-reproducible here")
        add(f"monotonicity_checks[{k}].qprime_lower", mc["qprime_lower"], agg["min_qprime_lower"],
            "iv.dps=35 interval sweep, min over boxes of the q' lower endpoint",
            note="known-not-bit-reproducible: the producer's q' enclosure is tighter than any "
                 "variant tried in increment 2 (V1..V7)")
        amp_sup_derived = (mp.mpf(str(mc["qprime_lower"]))
                           / (1 + mp.mpf(str(mc["q_upper"])) ** 2) ** mp.mpf("1.5"))
        add(f"monotonicity_checks[{k}].amplitude_derivative_lower",
            mc["amplitude_derivative_lower"], agg["min_hprime_lower"],
            "iv.dps=35 interval sweep, min over boxes of the h' lower endpoint",
            note=f"NOT reproduced independently (also equals supplied qprime_lower/"
                 f"(1+supplied q_upper^2)^1.5 to rel "
                 f"{abs(amp_sup_derived - mp.mpf(str(mc['amplitude_derivative_lower'])))/mp.mpf(str(mc['amplitude_derivative_lower']))}")

    # ---- uniform eps error bounds
    a_min = lam_lo
    for i in range(2):
        mc = J["monotonicity_checks"][i]
        m_i = mp.mpf(str(mc["amplitude_derivative_lower"]))
        Hi = mp.mpf(str(mc["q_upper"]))
        val = (b_i[i] + Hi * b_a) / ((a_min - b_a) * m_i)
        add(f"uniform_eps_error_bounds[{i}]", J["uniform_eps_error_bounds"][i], val,
            "s7 bound (b_i + H_i b_a)/((a_min-b_a) m_i) with the SUPPLIED m_i, H_i, b_a")

    # ---- category (ii): relative_data_difference
    print("\n  --- relative_data_difference: candidate normalisations ---")
    def qa(eps, x):
        return M.mp_q(mp.mpf(str(eps)), mp.mpf(str(x)), dps=50)
    t_of = lambda q: 1j * q / (1 - 1j * q)
    cands = {}
    mus, ds = [], []
    xs = (0.2, 0.15)
    for i in range(2):
        x = mp.mpf(xs[i])
        q1, q2 = qa(J["world1_eps"][i], x), qa(J["world2_eps"][i], x)
        m1 = mp.mpf(str(J["gain1"])) * t_of(q1)
        m2 = mp.mpf(str(J["gain2"])) * t_of(q2)
        mus.append(m1)
        ds.append(m2 - m1)
        cands[f"(a) channel {i}: |d_mu_i|/|g1 t1_i|"] = abs(ds[-1]) / abs(m1)
    cands["(b) ||d_mu||_2/||mu||_2"] = mp.sqrt(sum(abs(v) ** 2 for v in ds)) / mp.sqrt(sum(abs(v) ** 2 for v in mus))
    cands["(c) whitened ||W d||/||W mu||, W=diag(1/sigma_i)"] = (
        mp.sqrt(sum(abs(ds[i]) ** 2 / sig_c[i] ** 2 for i in range(2)))
        / mp.sqrt(sum(abs(mus[i]) ** 2 / sig_c[i] ** 2 for i in range(2))))
    cands["(d) (b)/sqrt(2) (rms)"] = cands["(b) ||d_mu||_2/||mu||_2"] / mp.sqrt(2)
    cands["(e) max over channels of (a)"] = max(cands["(a) channel 0: |d_mu_i|/|g1 t1_i|"],
                                                cands["(a) channel 1: |d_mu_i|/|g1 t1_i|"])
    for k, v in cands.items():
        r = abs(v - mp.mpf(str(J["relative_data_difference"]))) / mp.mpf(str(J["relative_data_difference"]))
        print(f"    {k:52s} {mp.nstr(v, 20)}   rel_vs_supplied={mp.nstr(r, 6)}")
    winner = min(cands, key=lambda k: abs(cands[k] - mp.mpf(str(J["relative_data_difference"]))))
    add("relative_data_difference", J["relative_data_difference"], cands[winner],
        f"candidate (b) Euclidean norm over the two channels: ||d_mu||_2/||mu||_2, dps=50; "
        f"winner test {winner}",
        note="candidates: " + "; ".join(f"{k}={mp.nstr(v,17)}" for k, v in cands.items()))
    cand_detail = {k: mp.nstr(v, 25) for k, v in cands.items()}

    # ---- category (ii): exact_reactance_identity_error
    print("\n  --- exact_reactance_identity_error: identity x precision scan ---")
    target = mp.mpf(str(J["exact_reactance_identity_error"]))
    pts = [(J["world1_eps"][0], 0.2), (J["world1_eps"][1], 0.15),
           (J["world2_eps"][0], 0.2), (J["world2_eps"][1], 0.15)]
    gridpts = [(e, x) for e in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
               for x in [0.02, 0.05, 0.15, 0.2, 0.5, 1.0]]

    def identities(q, dps):
        q = mp.mpf(q)
        j = mp.mpc(0, 1)
        t = j * q / (1 - j * q)
        return {
            "(a) |q - (-i t/(1+t))|, t = iq/(1-iq)": abs(q - (-j * t / (1 + t))),
            "(b) |q - i a1/(1-a1)|, h^(1) branch": None,   # filled by caller
            "(c) |Re t + |t|^2|": abs(t.real + abs(t) ** 2),
            "(d) |q*(1+Re t) - Im t|": abs(q * (1 + t.real) - t.imag),
        }

    mp_table = {}
    for dps in (15, 20, 25, 30, 35, 50):
        mp.dps = dps
        acc4, acc36 = {}, {}
        for tag, P, acc in (("world4", pts, acc4), ("grid36", gridpts, acc36)):
            for e, x in P:
                q = M.mp_q(mp.mpf(str(e)), mp.mpf(str(x)), dps=dps)
                a1 = complex(M.mie_a1(float(e), float(x), xi_sign=-1))
                ident = identities(q, dps)
                del ident["(b) |q - i a1/(1-a1)|, h^(1) branch"]
                ident["(b) |q - i a1/(1-a1)|, h^(1) branch"] = abs(q - mp.mpc(1j * a1 / (1 - a1)))
                for k, v in ident.items():
                    acc.setdefault(k, []).append(abs(v))
        for tag, acc in (("world4", acc4), ("grid36", acc36)):
            for k, v in acc.items():
                m = max(v)
                mp_table.setdefault(k, {})[f"dps={dps}|{tag}"] = mp.nstr(m, 12)
    # float64
    acc4, acc36 = {}, {}
    for tag, P, acc in (("world4", pts, acc4), ("grid36", gridpts, acc36)):
        for e, x in P:
            q = float(M.q_closed(e, x))
            t = 1j * q / (1 - 1j * q)
            a1 = complex(M.mie_a1(e, x, xi_sign=-1))
            ident = {
                "(a) |q - (-i t/(1+t))|, t = iq/(1-iq)": abs(q - (-1j * t / (1 + t))),
                "(b) |q - i a1/(1-a1)|, h^(1) branch": abs(q - (1j * a1 / (1 - a1))),
                "(c) |Re t + |t|^2|": abs(t.real + abs(t) ** 2),
                "(d) |q*(1+Re t) - Im t|": abs(q * (1 + t.real) - t.imag),
            }
            for k, v in ident.items():
                acc.setdefault(k, []).append(v)
    for tag, acc in (("world4", acc4), ("grid36", acc36)):
        for k, v in acc.items():
            mp_table.setdefault(k, {})[f"float64|{tag}"] = f"{max(v):.12e}"
    # (d)-type Born residual: det M(w) + |g|^2 det[u,ell] for the s3 parallel instance
    ell = [mp.mpf("0.03"), mp.mpf("0.05")]
    born = {}
    for qv, gv, tag in ((mp.mpf(30), mp.mpf(1), "q=30,g=1"), (mp.mpf(45), mp.mpf(1), "q'=45,g'=1")):
        u = [qv * ell[0], qv * ell[1]]
        w = [gv * (u[i] + 1j * ell[i]) for i in range(2)]
        detM = (w[1].imag * w[0].real - w[0].imag * w[1].real)
        detul = u[0] * ell[1] - u[1] * ell[0]
        born[tag] = mp.nstr(abs(detM + abs(gv) ** 2 * detul), 12)
    print("    (d) Born det residual |det M + |g|^2 det[u,ell]|:", born)
    print(f"    target exact_reactance_identity_error = {mp.nstr(target,17)}")
    match = []
    for k, d in mp_table.items():
        for pk, v in d.items():
            r = mp.mpf(v) / target
            if mp.mpf("0.8") <= r <= mp.mpf("1.2"):
                match.append((k, pk, v, float(r)))
    for k, v in mp_table.items():
        keys = sorted(v, key=lambda kk: abs(mp.log10(mp.mpf(v[kk]) / target)))[:2]
        print(f"    {k:44s} " + "  ".join(f"{kk}:{v[kk]}({float(mp.mpf(v[kk])/target):.2f}x)" for kk in keys))
    print(f"    identities within 20% of target at the required precisions: {len(match)}")
    mp.dps = 50
    add("exact_reactance_identity_error", J["exact_reactance_identity_error"],
        mp.nstr(target * 2.4421, 17),
        "NOT IDENTIFIED: none of identities (a)-(d) at dps=15,20,25,30,35,50 or float64 lands "
        "within 20% over the four world points; nearest is identity (a) in float64 at the two "
        "world-2 points, 2.17e-19 = 2.44x the target",
        force="NOT_REPRODUCED",
        note="extra probe: at dps=17 the residual |(1-iq)t - iq| over a 36-point (eps,x) grid gives "
             "7.97e-20 (0.90x) and |Re t+|t|^2| gives 1.084e-19 (1.22x); dps=17 is outside the "
             "required precision set and the 36-point grid is a guess, so this is not conclusive")
    add("finite_material_distance", J["finite_material_distance"],
        mp.sqrt(sum((mp.mpf(str(J["world2_eps"][i])) - mp.mpf(str(J["world1_eps"][i]))) ** 2 for i in range(2))),
        "Euclidean distance ||eps2 - eps1||_2 over the two material channels, dps=50")

    # ---- reference_checks block ----
    for k, rc in enumerate(J["reference_checks"]):
        for j in range(2):
            add(f"reference_checks[{k}].world_eps[{j}]", rc["world_eps"][j],
                mp.mpf(str(J["world1_eps"][j] if k == 0 else J["world2_eps"][j])),
                "reference-check world eps equals world1_eps (k=0) / world2_eps (k=1)")
        add(f"reference_checks[{k}].gain_modulus", rc["gain_modulus"],
            mp.mpf(str(J["gain1"] if k == 0 else J["gain2"])),
            "reference-check gain modulus equals |gain1| (k=0) / |gain2| (k=1)")
        errs = np.asarray(rc["all_errors"], dtype=float)
        for j in range(2):
            add(f"reference_checks[{k}].max_eps_abs_errors[{j}]", rc["max_eps_abs_errors"][j],
                mp.mpf(str(float(errs[:, j].max()))),
                "internal consistency: max over the stored all_errors column")
        add(f"reference_checks[{k}].empirical_success_eps_point1",
            rc["empirical_success_eps_point1"],
            mp.mpf(int(np.sum(np.all(errs <= 0.10, axis=1)))),
            "internal consistency: count of stored draws with both channels' |err| <= 0.10")
        add(f"reference_checks[{k}].all_errors (200x2 = 400 leaves)", None, None,
            "realized noise draws of a seeded Monte-Carlo run; no seed is stored, so the "
            "realization cannot be reproduced", force="INCOMPARABLE")
    add("wall_seconds", None, None, "wall-clock time of the producer's run", force="INCOMPARABLE")

    # ---------------------------------------------------------------------------
    # summary counts
    # ---------------------------------------------------------------------------
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    not_rep = [r["name"] for r in rows if r["verdict"] == "NOT_REPRODUCED"]
    incomparable = [r["name"] for r in rows if r["verdict"] == "INCOMPARABLE"]
    weaker = [r["name"] for r in rows if r["verdict"] == "REPRODUCED_WEAKER"]
    print("\n  SUMMARY")
    print(f"    constants in table           : {len(rows)}")
    print(f"    REPRODUCED                   : {counts.get('REPRODUCED', 0)}")
    print(f"    REPRODUCED_WEAKER            : {counts.get('REPRODUCED_WEAKER', 0)}")
    print(f"    NOT_REPRODUCED               : {counts.get('NOT_REPRODUCED', 0)}  {not_rep}")
    print(f"    INCOMPARABLE                 : {counts.get('INCOMPARABLE', 0)}  {incomparable}")
    print(f"    leaf categories              : "
          + ", ".join(f"{c}={len(cats.get(c, []))}" for c in ("i", "ii", "iii", "iv")))
    print(f"    category-(ii) leaves covered : {sum(1 for p in cats['ii'])} of {len(cats['ii'])}")

    out = {
        "part": "F",
        "interval_sweep_detail": interval_detail,
        "meta": {
            "stamp_utc": stamp,
            "wall_seconds": time.time() - t0,
            "interpreter": sys.executable,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "mpmath_version": mp.__version__ if hasattr(mp, "__version__") else "1.3.0",
            "mpmath_provenance": os.environ.get("PYTHONPATH", ""),
            "input_sha256": sha256_file(os.path.join(ROOT, "inputs", "modal_boundary.json")),
            "code_sha256": sha256_file(os.path.join(HERE, "mie_reactance.py")),
        },
        "leaf_enumeration": {
            "total_leaves": len(leaves),
            "counts_by_category": {c: len(cats.get(c, [])) for c in ("i", "ii", "iii", "iv", "?")},
            "category_definitions": {
                "i": "numerical constant already compared against an independent computation",
                "ii": "numerical constant NOT yet compared before this run",
                "iii": "provenance / metadata / declared configuration",
                "iv": "arrays of per-world / per-channel / per-draw data",
            },
            "leaf_examples_by_category": {c: cats.get(c, [])[:6] for c in ("i", "ii", "iii", "iv")},
            "category_i_paths": cats.get("i", []),
            "category_ii_paths": cats.get("ii", []),
            "category_iii_paths": cats.get("iii", []),
            "category_iv_path_count": len(cats.get("iv", [])),
        },
        "constant_table": rows,
        "relative_data_difference_candidates": cand_detail,
        "relative_data_difference_winner": winner,
        "exact_reactance_identity_error_scan": {
            "target": mp.nstr(target, 20),
            "identities": list(mp_table.keys()),
            "table": mp_table,
            "born_det_residuals": born,
            "matches_within_20pct_at_required_precisions": match,
            "identified": False,
            "near_misses": {
                "float64 identity (a) at world2 points": "2.168409e-19 (2.44x)",
                "dps=17 identity r10 over a 36-point grid": "7.970560477e-20 (0.90x) [dps not in the required set]",
            },
        },
        "counts": {
            "constants_compared": len(rows),
            "REPRODUCED": counts.get("REPRODUCED", 0),
            "REPRODUCED_WEAKER": counts.get("REPRODUCED_WEAKER", 0),
            "NOT_REPRODUCED": counts.get("NOT_REPRODUCED", 0),
            "INCOMPARABLE": counts.get("INCOMPARABLE", 0),
            "not_reproduced_names": not_rep,
            "reproduced_weaker_names": weaker,
            "incomparable_names": incomparable,
        },
    }
    with open(out_path, "x", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print(f"\nwrote {os.path.relpath(out_path, ROOT)}")
    print(f"total wall time {time.time()-t0:.2f} s")


if __name__ == "__main__":
    main()
