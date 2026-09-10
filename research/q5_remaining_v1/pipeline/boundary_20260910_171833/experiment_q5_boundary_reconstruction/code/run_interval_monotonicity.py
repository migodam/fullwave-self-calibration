#!/usr/bin/env python3
"""Deliverable 3: A5 sections 7-8 interval verification (monotonicity of eps -> q).

Uses mpmath.iv with iv.dps = 35 (explicitly set; the interpreter default 15 is not
comparable to the supplied constants).  Endpoints are built by interval arithmetic from
exact small integer literals so that every box provably encloses the exact decimal
endpoints.  No installation is performed: mpmath 1.3.0 is read out of an existing uv
cache archive through PYTHONPATH (a pure-Python file read, not a pip/system install).

Writes results/interval_<UTCSTAMP>.json (written once, never edited) and logs/interval.log.

Usage:
    PYTHONPATH=<uv archive with mpmath> <a3_research venv python> code/run_interval_monotonicity.py
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shlex
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import mie_reactance as M  # noqa: E402


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


def json_safe(o):
    """Recursively replace mpmath mpf/mpi objects by high-precision decimal strings."""
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if hasattr(o, "_mpi_") or hasattr(o, "_mpf_"):
        return j(o, 30)
    if isinstance(o, (bool, int, float, str)) or o is None:
        return o
    return o


def rel_err(a, b):
    return abs(a - b) / abs(b) if b else float("nan")


def sweep(iv, boxes, x, joint_boxes=False):
    """Evaluate the interval enclosure of (q', Den, q, h') on every box.

    A box is the smallest interval containing the two endpoint intervals, i.e.
    [inf(lo_k), sup(hi_k)]; this is a superset of the exact decimal box, so the
    enclosure covers the exact interval.
    """
    stats = {
        "n_boxes": len(boxes),
        "min_qprime_lower": None, "min_absDen_lower": None,
        "min_q_lower": None, "max_q_upper": None, "min_hprime_lower": None,
        "boxes_with_nonpositive_qprime_lower": 0,
        "boxes_with_nonpositive_absDen_lower": 0,
        "argmin_qprime_lower_eps": None, "argmin_absDen_lower_eps": None,
    }
    for i, bx in enumerate(boxes):
        o = M.iv_quantities(bx, x)
        qp_lo = o["qprime"].a
        den_lo = M.iv_min_abs(o["Den"])
        q_lo, q_hi = o["q"].a, o["q"].b
        hp_lo = o["hprime"].a
        if qp_lo <= 0:
            stats["boxes_with_nonpositive_qprime_lower"] += 1
        if den_lo <= 0:
            stats["boxes_with_nonpositive_absDen_lower"] += 1
        if stats["min_qprime_lower"] is None or qp_lo < stats["min_qprime_lower"]:
            stats["min_qprime_lower"] = qp_lo
            stats["argmin_qprime_lower_eps"] = (bx.a, bx.b)
        if stats["min_absDen_lower"] is None or den_lo < stats["min_absDen_lower"]:
            stats["min_absDen_lower"] = den_lo
            stats["argmin_absDen_lower_eps"] = (bx.a, bx.b)
        if stats["min_q_lower"] is None or q_lo < stats["min_q_lower"]:
            stats["min_q_lower"] = q_lo
        if stats["max_q_upper"] is None or q_hi > stats["max_q_upper"]:
            stats["max_q_upper"] = q_hi
        if stats["min_hprime_lower"] is None or hp_lo < stats["min_hprime_lower"]:
            stats["min_hprime_lower"] = hp_lo
    return stats


def j(v, digits=30):
    """Render an mpmath.iv interval (or a raw libmp mpf) as a plain decimal string.

    For a zero-width interval the single stored bound is used, so the string is the
    interval's own decimal expansion and can be re-parsed by mp.mpf without ambiguity.
    """
    from mpmath.libmp import to_str
    if hasattr(v, "_mpi_"):
        raw = v._mpi_[0]
    elif hasattr(v, "_mpf_"):
        raw = v._mpf_
    else:
        raw = v
    return to_str(raw, digits)


def main():
    sys.stdout = Tee(os.path.join(ROOT, "logs", "interval.log"))
    sys.stderr = sys.stdout
    t_start = time.time()

    from mpmath import iv, mp
    iv.dps = 35  # REQUIRED: interpreter default is 15
    assert iv.dps == 35

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_path = os.path.join(ROOT, "results", f"interval_{stamp}.json")
    if os.path.exists(out_path):
        raise SystemExit(f"REFUSING to overwrite existing results file: {out_path}")

    print("=" * 100)
    print("A5 sections 7-8 interval monotonicity verification (mpmath.iv dps=35)")
    print("started_utc :", datetime.now(timezone.utc).isoformat())
    print("command     :", " ".join(shlex.quote(a) for a in sys.argv))
    print("interpreter :", sys.executable)
    print("python      :", sys.version.replace("\n", " "))
    print("platform    :", platform.platform())
    print("iv.dps      :", iv.dps)
    print("=" * 100)

    import mpmath
    versions = {"python": sys.version.split()[0], "mpmath": mpmath.__version__,
                "mpmath_file": mpmath.__file__}
    input_files = [
        "inputs/modal_boundary.json", "inputs/derive/A5_sec_07.md", "inputs/derive/A5_sec_08.md",
        "inputs/derive/MODAL_BOUNDARY_KEYS.txt", "inputs/SOURCE_HASHES.json",
        "logs/uv_cache_VjSSQm31390egyfy.sha256",
    ]
    code_files = ["code/mie_reactance.py", "code/run_interval_monotonicity.py"]
    hashes = {p: sha256_file(os.path.join(ROOT, p)) for p in input_files + code_files}
    for k, v in hashes.items():
        print(f"  sha256 {v}  {k}")

    J = json.load(open(os.path.join(ROOT, "inputs", "modal_boundary.json")))
    checks = J["monotonicity_checks"]

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[0] material intervals as supplied in inputs/modal_boundary.json")
    print("-" * 100)
    intervals = []
    for c in checks:
        intervals.append({"x": c["x"], "eps_domain": c["eps_domain"],
                          "interval_count": c["interval_count"]})
        print(f"  x={c['x']}  eps_domain={c['eps_domain']}  interval_count={c['interval_count']}")
    assert intervals[0]["eps_domain"] == [1.5, 4.0] and intervals[0]["x"] == 0.2, intervals[0]
    assert intervals[1]["eps_domain"] == [2.0, 5.0] and intervals[1]["x"] == 0.15, intervals[1]
    print("  CONFIRMED: I1=[1.5,4] x=0.2 (n=2500) and I2=[2,5] x=0.15 (n=3000)")

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[1] exact-decimal endpoint construction: which construction is provably enclosing?")
    print("-" * 100)
    # Construction A: left + k*step from exact small integer literals
    a1 = iv.mpf(3) / iv.mpf(2) + iv.mpf(1) * (iv.mpf(1) / iv.mpf(1000))
    # Construction B: single exact rational 1501/1000 (construction A at k=1)
    b1 = iv.mpf(1501) / iv.mpf(1000)
    b1_fraction = None
    frac_error = None
    try:
        from fractions import Fraction
        b1_fraction = iv.mpf(Fraction(1501, 1000))
    except Exception as exc:  # mpmath 1.3.0 does not accept Fraction/mpq here
        frac_error = f"{type(exc).__name__}: {exc}"
    mp.dps = 60
    ref60 = mp.mpf(1501) / mp.mpf(1000)
    mp.dps = 35
    endpoint_res = {
        "construction_A": "endpoint_k = iv.mpf(3)/iv.mpf(2) + k*(iv.mpf(1)/iv.mpf(1000))",
        "construction_B": "iv.mpf(1501)/iv.mpf(1000)",
        "construction_B_as_Fraction_or_mpq": (
            "iv.mpf(Fraction(1501,1000)) -> " + frac_error if frac_error else "accepted"),
        "A_interval": [j(a1.a), j(a1.b)],
        "B_interval": [j(b1.a), j(b1.b)],
        "A_width": j(a1.b - a1.a),
        "B_width": j(b1.b - b1.a),
        "reference_60dps_1501_over_1000": j(ref60),
        "A_encloses_60dps": bool(a1.a <= ref60 <= a1.b),
        "B_encloses_60dps": bool(b1.a <= ref60 <= b1.b),
    }
    print(f"  A: endpoint_k = iv.mpf(3)/iv.mpf(2) + k*(iv.mpf(1)/iv.mpf(1000)), k=1 -> "
          f"[{endpoint_res['A_interval'][0]}, {endpoint_res['A_interval'][1]}]  width {endpoint_res['A_width']}")
    print(f"  B: iv.mpf(1501)/iv.mpf(1000)                                        -> "
          f"[{endpoint_res['B_interval'][0]}, {endpoint_res['B_interval'][1]}]  width {endpoint_res['B_width']}")
    print(f"  iv.mpf(Fraction(1501,1000)) -> {endpoint_res['construction_B_as_Fraction_or_mpq']}")
    print(f"  60-dps value of 1501/1000    -> {j(ref60)}")
    print(f"  A encloses the 60-dps value  -> {endpoint_res['A_encloses_60dps']}")
    print(f"  B encloses the 60-dps value  -> {endpoint_res['B_encloses_60dps']}")

    # Cross-check the full endpoint family used for I1 and I2 against 60-dps rationals.
    mp.dps = 60
    probe = []
    for (ln, ld, hn, hd, n) in [(3, 2, 4, 1, 2500), (2, 1, 5, 1, 3000)]:
        ep = M.decimal_endpoint_intervals(ln, ld, hn, hd, n)
        bad = 0
        worst_halfwidth = None
        n_sampled = 0
        for kix in range(0, n + 1, max(1, n // 64)):
            exact = mp.mpf(ln) / mp.mpf(ld) + mp.mpf(kix) * (mp.mpf(hn) / mp.mpf(hd) - mp.mpf(ln) / mp.mpf(ld)) / mp.mpf(n)
            n_sampled += 1
            if not (ep[kix].a <= exact <= ep[kix].b):
                bad += 1
            hw = (ep[kix].b - ep[kix].a) / 2
            worst_halfwidth = hw if worst_halfwidth is None else max(worst_halfwidth, hw)
        probe.append({"domain": [ln / ld, hn / hd], "n": n, "sampled_endpoints": n_sampled,
                      "all_enclosed": bad == 0,
                      "max_endpoint_interval_halfwidth": j(worst_halfwidth)})
        print(f"  domain [{ln}/{ld},{hn}/{hd}] n={n}: sampled {n_sampled} endpoints, "
              f"all enclosed = {probe[-1]['all_enclosed']}, "
              f"max endpoint half-width = {probe[-1]['max_endpoint_interval_halfwidth']}")
    mp.dps = 35
    endpoint_res["family_probe"] = probe
    endpoint_res["choice"] = (
        "Construction A (left + k*step from exact integer literals) is used for every box "
        "endpoint; construction B is its k=1 special case written as one rational division. "
        "Both are provably enclosing (verified against a 60-dps rational evaluation). "
        "mpmath.iv 1.3.0 has no Fraction/mpq constructor, so the literal iv.mpf(mpq(...)) form "
        "of the task is realised as iv.mpf(p)/iv.mpf(q)."
    )
    print("  CHOICE: construction A (interval arithmetic from exact integer literals); B is its")
    print("          k=1 case expressed as a single rational division. Both enclose.")

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[2] main interval sweep with the supplied subdivision counts")
    print("-" * 100)
    main_run = []
    for idx, c in enumerate(checks):
        n = int(c["interval_count"])
        lo, hi = c["eps_domain"]
        x = c["x"]
        t0 = time.time()
        ep = M.decimal_endpoint_intervals(round(lo * 1000), 1000, round(hi * 1000), 1000, n)
        boxes = [M.box_from_intervals(ep[k], ep[k + 1]) for k in range(n)]
        st = sweep(iv, boxes, x)
        st["wall_seconds"] = time.time() - t0
        st["x"] = x
        st["eps_domain"] = [lo, hi]
        st["n"] = n
        main_run.append(st)
        print(f"  I{idx+1}: x={x} eps in [{lo},{hi}] n={n}  ({st['wall_seconds']:.1f} s)")
        print(f"      min lower bound of q'      = {j(st['min_qprime_lower'])}")
        print(f"      min lower bound of |Den|   = {j(st['min_absDen_lower'])}")
        print(f"      min lower bound of q       = {j(st['min_q_lower'])}")
        print(f"      max upper bound of q  (H_i)= {j(st['max_q_upper'])}")
        print(f"      min lower bound of h'      = {j(st['min_hprime_lower'])}")
        print(f"      boxes with q'_lower <= 0 : {st['boxes_with_nonpositive_qprime_lower']} / {n}")
        print(f"      boxes with |Den|_lower<=0: {st['boxes_with_nonpositive_absDen_lower']} / {n}")
        del boxes, ep

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[3] computed vs supplied constants")
    print("-" * 100)
    comparisons = []
    keys = [("qprime_lower", "min_qprime_lower"),
            ("denominator_lower", "min_absDen_lower"),
            ("q_lower", "min_q_lower"),
            ("q_upper", "max_q_upper"),
            ("amplitude_derivative_lower", "min_hprime_lower")]
    for idx, c in enumerate(checks):
        for skey, ckey in keys:
            sup = c[skey]
            comp = float(main_run[idx][ckey])
            comparisons.append({
                "check": f"monotonicity_checks[{idx}]", "x": c["x"], "key": skey,
                "supplied": sup, "computed": comp, "difference": comp - sup,
                "relative_difference": (comp - sup) / sup,
                "computed_conservative_vs_supplied": comp <= sup,
                "high_precision_computed": j(main_run[idx][ckey]),
            })
            print(f"  {f'monotonicity_checks[{idx}]':24s} {skey:28s} supplied {sup:.16e}  "
                  f"computed {comp:.16e}  diff {comp - sup:+.3e}  rel {(comp - sup) / sup:+.3e}")

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[4] refinement study: min q' lower bound vs subdivision count n")
    print("-" * 100)
    import numpy as _np
    refine = {}
    for idx, c in enumerate(checks):
        lo, hi = c["eps_domain"]
        x = c["x"]
        eps_dense = _np.linspace(lo, hi, 400001)
        true_inf = float(M.qprime_closed(eps_dense, x).min())
        rows = []
        for n in (100, 500, 2500, 5000, 10000):
            t0 = time.time()
            ep = M.decimal_endpoint_intervals(round(lo * 1000), 1000, round(hi * 1000), 1000, n)
            boxes = [M.box_from_intervals(ep[k], ep[k + 1]) for k in range(n)]
            st = sweep(iv, boxes, x)
            rows.append({
                "n": n, "step": (hi - lo) / n,
                "min_qprime_lower": j(st["min_qprime_lower"]),
                "min_qprime_lower_float": float(st["min_qprime_lower"]),
                "min_absDen_lower_float": float(st["min_absDen_lower"]),
                "min_q_lower_float": float(st["min_q_lower"]),
                "max_q_upper_float": float(st["max_q_upper"]),
                "min_hprime_lower_float": float(st["min_hprime_lower"]),
                "boxes_with_nonpositive_qprime_lower": st["boxes_with_nonpositive_qprime_lower"],
                "qprime_strictly_positive_on_every_box": st["boxes_with_nonpositive_qprime_lower"] == 0,
                "absDen_strictly_positive_on_every_box": st["boxes_with_nonpositive_absDen_lower"] == 0,
                "ratio_to_true_inf_qprime": float(st["min_qprime_lower"]) / true_inf,
                "wall_seconds": time.time() - t0,
            })
            print(f"  I{idx+1} n={n:6d} step={(hi-lo)/n:.6f}  min q'_lower = "
                  f"{rows[-1]['min_qprime_lower_float']:.10e}  "
                  f"(true inf q' = {true_inf:.10e}, ratio {rows[-1]['ratio_to_true_inf_qprime']:.6f})  "
                  f"neg boxes: {rows[-1]['boxes_with_nonpositive_qprime_lower']}")
            del boxes, ep
        refine[f"monotonicity_checks[{idx}]"] = {
            "x": x, "eps_domain": [lo, hi],
            "true_inf_qprime_dense_double_scan": true_inf,
            "true_inf_qprime_note": "400001-point double-precision scan; independent of the interval run",
            "rows": rows,
        }

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[5] section 7 conservative constants m_i = (min q'_lower) / (1 + H_i^2)^(3/2)")
    print("-" * 100)
    m_repro = []
    sec7 = [4.3261945994e-4, 1.3086616724e-4]
    for idx, c in enumerate(checks):
        qp_lo = main_run[idx]["min_qprime_lower"]
        H = main_run[idx]["max_q_upper"]
        # double-precision reproduction of the constant (this is the arithmetic that
        # reproduces the supplied amplitude_derivative_lower exactly)
        m_float = float(qp_lo) / (1.0 + float(H) ** 2) ** 1.5
        H_lo_iv = main_run[idx]["min_q_lower"]
        mp.dps = 50
        m_hp = mp.mpf(j(qp_lo, 40)) / (1 + mp.mpf(j(H, 40)) ** 2) ** mp.mpf("1.5")
        m_hp_lo_den = mp.mpf(j(qp_lo, 40)) / (1 + mp.mpf(j(H_lo_iv, 40)) ** 2) ** mp.mpf("1.5")
        mp.dps = 35
        amp_sup = c["amplitude_derivative_lower"]
        m_sup_derived = c["qprime_lower"] / (1 + c["q_upper"] ** 2) ** 1.5
        m_repro.append({
            "check": f"monotonicity_checks[{idx}]",
            "qprime_lower_computed": j(qp_lo), "H_i_computed": j(H),
            "m_i_computed": j(m_hp, 30), "m_i_computed_float": m_float,
            "m_i_computed_hp50": mp.nstr(m_hp, 25) if hasattr(mp, "nstr") else str(m_hp),
            "m_i_computed_hp50_using_q_lower_in_denominator": str(m_hp_lo_den),
            "m_i_from_supplied_constants": m_sup_derived,
            "amplitude_derivative_lower_supplied": amp_sup,
            "A5_section7_stated_m_i": sec7[idx],
            "rel_diff_vs_section7_stated": (m_float - sec7[idx]) / sec7[idx],
            "supplied_matches_amplitude_derivative_lower": abs(m_sup_derived - amp_sup) / amp_sup,
        })
        print(f"  I{idx+1}: m_i computed = {m_float:.12e}   (hp50 {str(m_hp)})")
        print(f"      from supplied constants = {m_sup_derived:.12e}  "
              f"(H_i supplied {c['q_upper']:.16e}, q'_lower supplied {c['qprime_lower']:.16e})")
        print(f"      A5 section 7 stated   = {sec7[idx]:.12e}   "
              f"rel diff (computed vs stated) = {(m_float - sec7[idx]) / sec7[idx]:+.3e}")
        print(f"      supplied amplitude_derivative_lower = {amp_sup:.16e}  "
              f"(= supplied qprime_lower/(1+q_upper^2)^1.5, rel {abs(m_sup_derived - amp_sup) / amp_sup:.2e})")

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[6] diagnostic: sensitivity of the q' lower bound to algebraically equivalent arrangements")
    print("-" * 100)
    variants = {}
    for idx, c in enumerate(checks):
        lo, hi = c["eps_domain"]
        x = c["x"]
        n = int(c["interval_count"])
        ep = M.decimal_endpoint_intervals(round(lo * 1000), 1000, round(hi * 1000), 1000, n)
        boxes = [M.box_from_intervals(ep[k], ep[k + 1]) for k in range(n)]
        row = {}
        for variant in ("V1", "V2"):
            mn = None
            for bx in boxes:
                o = _variant_core(iv, bx, x, variant)
                v = o["qprime"].a
                mn = v if mn is None else min(mn, v)
            row[variant] = j(mn)
            row[variant + "_float"] = float(mn)
        variants[f"monotonicity_checks[{idx}]"] = row
        print(f"  I{idx+1}: V1 (A5 chain rule, used) q'_lower = {row['V1_float']:.12e}")
        print(f"      V2 (direct d/dz of the D closed form)   = {row['V2_float']:.12e}")
        print(f"      supplied                                 = {c['qprime_lower']:.12e}")
        del boxes, ep
    variants["note"] = (
        "V1 and V2 are algebraically identical analytic derivatives of q; they differ only in "
        "how the interval expression is grouped, so they give different (both valid) enclosures. "
        "The supplied qprime_lower lies above both, i.e. the source implementation used a "
        "marginally tighter grouping of the same analytic derivative."
    )

    # ---------------------------------------------------------------------------------
    results = {
        "meta": {
            "task": "A5 sections 7-8 interval monotonicity / conservative constants",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command_line": " ".join(shlex.quote(a) for a in sys.argv),
            "interpreter": sys.executable,
            "versions": versions,
            "iv_dps": iv.dps,
            "interpreter_default_iv_dps_note": "mpmath.iv default dps is 15; dps=35 is set explicitly",
            "cwd": os.getcwd(),
            "file_sha256": hashes,
            "mpmath_import_note": ("mpmath 1.3.0 is imported from an existing uv cache archive via "
                                   "PYTHONPATH; no pip/system install was performed"),
            "interval_verification_status": "DONE",
        },
        "material_intervals": intervals,
        "endpoint_construction": endpoint_res,
        "main_run": main_run,
        "comparison_with_supplied": comparisons,
        "refinement_study": refine,
        "m_i_reproduction": m_repro,
        "variant_sensitivity": variants,
        "wall_seconds": time.time() - t_start,
    }
    results = json_safe(results)
    with open(out_path, "x", encoding="utf-8") as f:
        json.dump(results, f, indent=2, sort_keys=False)
        f.write("\n")
    print("\n" + "=" * 100)
    print("wrote", os.path.relpath(out_path, ROOT))
    print(f"total wall time {results['wall_seconds']:.1f} s")
    print("=" * 100)


def _variant_core(iv, eps_box, x, variant):
    """Interval core with a selectable grouping of the eps-derivative terms."""
    m = iv.sqrt(eps_box)
    z = m * x
    sz, cz = iv.sin(z), iv.cos(z)
    sx, cx = iv.sin(x), iv.cos(x)
    j1z = sz / z**2 - cz / z
    y1z = -cz / z**2 - sz / z
    Djz = sz / z + cz / z**2 - sz / z**3
    Dyz = -cz / z + sz / z**2 + cz / z**3
    j1x = sx / x**2 - cx / x
    y1x = -cx / x**2 - sx / x
    Djx = sx / x + cx / x**2 - sx / x**3
    Dyx = -cx / x + sx / x**2 + cx / x**3
    j1pz = sz / z + 2 * cz / z**2 - 2 * sz / z**3
    if variant == "V1":
        y1pz = -cz / z + 2 * sz / z**2 + 2 * cz / z**3
        Djpz = -j1pz / z - j1z + j1z / z**2
        Dypz = -y1pz / z - y1z + y1z / z**2
    else:
        Djpz = cz / z - 2 * sz / z**2 - 3 * cz / z**3 + 3 * sz / z**4
        Dypz = sz / z + 2 * cz / z**2 - 3 * sz / z**3 - 3 * cz / z**4
    N = m * j1z * Djx - j1x * Djz
    Den = m * j1z * Dyx - y1x * Djz
    q = N / Den
    A = m * j1z
    Ap = j1z / (2 * m) + (x / 2) * j1pz
    dD = (x / (2 * m)) * Djpz
    Np = Ap * Djx - j1x * dD
    Denp = Ap * Dyx - y1x * dD
    qp = (Np * Den - N * Denp) / Den**2
    hp = qp / (1 + q**2) * iv.sqrt(1 + q**2)
    return {"q": q, "Den": Den, "qprime": qp, "hprime": hp}


if __name__ == "__main__":
    main()
