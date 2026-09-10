#!/usr/bin/env python3
"""PART A: arrangement-dependence sweep of the interval lower bound of dq/deps.

Question: the supplied qprime_lower (inputs/modal_boundary.json -> monotonicity_checks[i]
.qprime_lower) is *tighter* than the value produced by code/mie_reactance.py at iv.dps=35.
Is it reproducible by SOME mathematically equivalent interval arrangement, and how much
does the answer depend on the arrangement?

Also: rigorously locate inf q' on both intervals (dense mp.dps=50 scan) and verify
q'' < 0 on the whole interval by interval arithmetic (so inf q' is at the right endpoint).

Usage:
    PYTHONPATH=<uv archive with mpmath> <a3 venv python> code/run_qprime_variant_sweep.py
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

import mie_reactance as M          # noqa: E402
import qprime_variants_lib as L    # noqa: E402

SUPPLIED_QPRIME_LOWER = [4.3262426032187571e-04, 1.3086650394735435e-04]
INTERVALS = [(3, 2, 4, 1, 2500, 0.2), (2, 1, 5, 1, 3000, 0.15)]  # lo,hi,n,x


class Sink:
    """Results file created on first write, then updated in place (partial progress)."""

    def __init__(self, path):
        self.path = path
        if os.path.exists(path):
            raise SystemExit(f"REFUSING to overwrite existing results file: {path}")
        self.data = {}
        self.flush()

    def put(self, key, value):
        self.data[key] = value
        self.flush()
        return value

    def flush(self):
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(json_safe(self.data), f, indent=2, sort_keys=False)
            f.write("\n")
        os.replace(tmp, self.path)


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


def dec30(v):
    """Render an mpmath mpf / iv value as a plain decimal string (30 digits).

    For an interval this renders BOTH endpoints explicitly; mpmath's own ``v.a``
    would silently give the lower bound of the lower endpoint.
    """
    from mpmath import mp
    from mpmath.libmp import to_str
    with mp.extradps(25):
        if hasattr(v, "_mpi_"):
            return "[%s, %s]" % (to_str(v._mpi_[0], 30), to_str(v._mpi_[1], 30))
        if hasattr(v, "_mpf_"):
            return to_str(v._mpf_, 30)
        if isinstance(v, float):
            return repr(v)
        return str(v)


def json_safe(o):
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if hasattr(o, "_mpi_") or hasattr(o, "_mpf_"):
        return dec30(o)
    if isinstance(o, (bool, int, float, str)) or o is None:
        return o
    return o


def rel_diff(a, b):
    return float((a - b) / b) if b else float("nan")


def boxes_for(lo_num, lo_den, hi_num, hi_den, n):
    ep = M.decimal_endpoint_intervals(lo_num, lo_den, hi_num, hi_den, n)
    return [M.box_from_intervals(ep[k], ep[k + 1]) for k in range(n)]


def sweep_variant(iv, boxes, x, variant):
    best = None
    best_box = None
    for bx in boxes:
        o = L.evaluate(bx, x, iv, variant)
        lo = L.iv_bounds(o["qprime"])[0]
        if best is None or lo < best:
            best, best_box = lo, L.iv_bounds(bx)
    return best, best_box


def main():
    sys.stdout = Tee(os.path.join(ROOT, "logs", "qprime_variants.log"))
    sys.stderr = sys.stdout
    t0 = time.time()

    from mpmath import iv, mp
    import mpmath

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_path = os.path.join(ROOT, "results", f"qprime_variants_{stamp}.json")
    sink = Sink(out_path)

    print("=" * 100)
    print("PART A  arrangement-dependence sweep of the interval lower bound of dq/deps")
    print("started_utc :", datetime.now(timezone.utc).isoformat())
    print("interpreter :", sys.executable)
    print("platform    :", platform.platform())
    print("mpmath      :", mpmath.__version__)
    print("=" * 100)

    hashes = {p: sha256_file(os.path.join(ROOT, p)) for p in [
        "code/mie_reactance.py", "code/qprime_variants_lib.py",
        "code/run_qprime_variant_sweep.py", "inputs/modal_boundary.json",
        "logs/uv_cache_VjSSQm31390egyfy.sha256"]}
    for k, v in hashes.items():
        print(f"  sha256 {v}  {k}")
    print("  quarantined: results/interval_20260910T093248Z.json -> "
          "results/PARTIAL_interval_20260910T093248Z.json.bad (truncated JSON, content unchanged)")

    J = json.load(open(os.path.join(ROOT, "inputs", "modal_boundary.json")))
    checks = J["monotonicity_checks"]
    for i, c in enumerate(checks):
        assert abs(c["qprime_lower"] - SUPPLIED_QPRIME_LOWER[i]) / SUPPLIED_QPRIME_LOWER[i] < 1e-9
    print("  supplied qprime_lower:", checks[0]["qprime_lower"], checks[1]["qprime_lower"])

    sink.put("meta", {
        "task": "PART A qprime variant sweep",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command_line": " ".join(shlex.quote(a) for a in sys.argv),
        "interpreter": sys.executable, "mpmath": mpmath.__version__,
        "platform": platform.platform(), "cwd": os.getcwd(),
        "file_sha256": hashes,
        "supplied_qprime_lower": SUPPLIED_QPRIME_LOWER,
        "interval_definition": [
            {"eps_domain": [1.5, 4.0], "x": 0.2, "n_boxes": 2500, "index": 0},
            {"eps_domain": [2.0, 5.0], "x": 0.15, "n_boxes": 3000, "index": 1}],
        "quarantined_file": ("results/interval_20260910T093248Z.json -> "
                             "results/PARTIAL_interval_20260910T093248Z.json.bad"),
        "rigor_note": ("every reported lower bound is the lower endpoint of an interval "
                       "enclosure computed on boxes that provably contain the exact decimal "
                       "sub-intervals; the boxes themselves are built at iv.dps=35 and are "
                       "reused unchanged when the arithmetic precision is varied"),
    })

    # ---------------------------------------------------------------- variant sweep
    iv.dps = 35
    results = {}
    boxes_by_i = {}
    for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
        boxes_by_i[idx] = boxes_for(ln, ld, hn, hd, n)
        print(f"\n  I{idx+1}: eps in [{ln}/{ld},{hn}/{hd}] x={x} n={n} boxes")
        results[f"I{idx+1}"] = {}

    # V1..V4 at dps=35
    for variant in ("V1", "V2", "V3", "V4"):
        for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
            iv.dps = 35
            val, bx = sweep_variant(iv, boxes_by_i[idx], x, variant)
            results[f"I{idx+1}"][variant] = {
                "min_lower": float(val), "min_lower_str": dec30(val),
                "rel_diff_vs_supplied": rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx]),
                "matches_supplied_within_1e-9": bool(
                    abs(rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx])) < 1e-9),
                "argmin_box": [dec30(bx[0]), dec30(bx[1])],
                "iv_dps": 35,
            }
            print(f"    {variant}: q'_lower = {float(val):.12e}  "
                  f"rel vs supplied = {rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx]):+.3e}")

    # V5: V1 at dps=25 and dps=50 (same boxes)
    for ic_dps in (25, 50):
        for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
            iv.dps = ic_dps
            val, bx = sweep_variant(iv, boxes_by_i[idx], x, "V1")
            results[f"I{idx+1}"][f"V5_dps{ic_dps}"] = {
                "min_lower": float(val), "min_lower_str": dec30(val),
                "rel_diff_vs_supplied": rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx]),
                "matches_supplied_within_1e-9": bool(
                    abs(rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx])) < 1e-9),
                "argmin_box": [dec30(bx[0]), dec30(bx[1])], "iv_dps": ic_dps}
            print(f"    V5@dps={ic_dps}: q'_lower = {float(val):.12e}  "
                  f"rel = {rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx]):+.3e}")
    iv.dps = 35
    sink.put("variant_min_lower_bounds", results)

    # V7: centered form eps = mid + delta
    v7 = {}
    for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
        iv.dps = 35
        best = None
        for bx in boxes_by_i[idx]:
            blo, bhi = L.iv_bounds(bx)
            mid = (blo + bhi) / 2
            delta = bx - iv.mpf([mid, mid])
            o = L.evaluate(iv.mpf([mid, mid]) + delta, x, iv, "V1")
            lo = L.iv_bounds(o["qprime"])[0]
            best = lo if best is None else min(best, lo)
        v7[f"I{idx+1}"] = {"min_lower": float(best), "min_lower_str": dec30(best),
                           "rel_diff_vs_supplied": rel_diff(float(best), SUPPLIED_QPRIME_LOWER[idx]),
                           "matches_supplied_within_1e-9": bool(
                               abs(rel_diff(float(best), SUPPLIED_QPRIME_LOWER[idx])) < 1e-9),
                           "iv_dps": 35}
        print(f"    V7 (centered form): q'_lower = {float(best):.12e}  "
              f"rel = {rel_diff(float(best), SUPPLIED_QPRIME_LOWER[idx]):+.3e}")
    sink.put("V7_centered_form", v7)

    # V6: midpoint Taylor with an interval bound on |q''|
    v6 = {f"I{i+1}": {} for i in range(2)}
    kappas = [("V6_halfwidth", 0.5), ("V6_0p5h", 0.5), ("V6_1p0h", 1.0)]
    for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
        best = {k[0]: None for k in kappas}
        best["V6_mid_reference"] = None
        for bx in boxes_by_i[idx]:
            iv.dps = 35
            qpp_box = L.qpp(bx, x, iv)
            qpp_abs_max = L.abs_interval_upper(qpp_box)
            blo, bhi = L.iv_bounds(bx)
            half = (bhi - blo) / 2
            iv.dps = 50
            mid = (blo + bhi) / 2
            qmid = L.evaluate(iv.mpf([mid, mid]), x, iv, "V1")["qprime"]
            qlo = L.iv_bounds(qmid)[0]
            for name, kap in kappas:
                lb = L.iv_bounds(qmid - iv.mpf(kap) * 2 * half * qpp_abs_max)[0]
                best[name] = lb if best[name] is None else min(best[name], lb)
            if best["V6_mid_reference"] is None or qlo < best["V6_mid_reference"]:
                best["V6_mid_reference"] = qlo
        for name in list(best):
            val = best[name]
            v6[f"I{idx+1}"][name] = {
                "min_lower": float(val), "min_lower_str": dec30(val),
                "rel_diff_vs_supplied": rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx]),
                "matches_supplied_within_1e-9": bool(
                    abs(rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx])) < 1e-9),
                "note": ("midpoint value only; NOT a rigorous lower bound" if
                         name == "V6_mid_reference" else
                         f"q'_mid - {name.split('_')[-1]}*|q''|_max (width = box width)")}
            print(f"    I{idx+1} {name}: {float(val):.12e}  "
                  f"rel vs supplied = {rel_diff(float(val), SUPPLIED_QPRIME_LOWER[idx]):+.3e}")
    iv.dps = 35
    sink.put("V6_midpoint_taylor", v6)

    # ------------------------------------------------------- range across variants
    rng = {}
    for idx in (0, 1):
        vals = {}
        for name, dd in results[f"I{idx+1}"].items():
            vals[name] = dd["min_lower"]
        vals["V7"] = v7[f"I{idx+1}"]["min_lower"]
        for name, dd in v6[f"I{idx+1}"].items():
            vals[name] = dd["min_lower"]
        rng[f"I{idx+1}"] = {
            "values": {k: float(v) for k, v in vals.items()},
            "min": float(min(vals.values())), "max": float(max(vals.values())),
            "range_max_minus_min": float(max(vals.values()) - min(vals.values())),
            "relative_range": float((max(vals.values()) - min(vals.values())) / min(vals.values())),
            "supplied": SUPPLIED_QPRIME_LOWER[idx],
            "best_matching_variants": sorted(
                [k for k, v in vals.items()
                 if abs(rel_diff(float(v), SUPPLIED_QPRIME_LOWER[idx])) < 1e-9]),
        }
        print(f"\n  I{idx+1} variant range: min={rng[f'I{idx+1}']['min']:.12e} "
              f"max={rng[f'I{idx+1}']['max']:.12e} "
              f"range={rng[f'I{idx+1}']['range_max_minus_min']:.3e} "
              f"supplied={SUPPLIED_QPRIME_LOWER[idx]:.12e}")
        print(f"      variants matching supplied within 1e-9 rel: "
              f"{rng[f'I{idx+1}']['best_matching_variants'] or 'NONE'}")
    sink.put("variant_range", rng)

    # ------------------------------------------------------------ dense point scan
    mp.dps = 50
    dense = {}
    for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
        NPT = 20001
        lo = mp.mpf(ln) / mp.mpf(ld)
        hi = mp.mpf(hn) / mp.mpf(hd)
        step = (hi - lo) / (NPT - 1)
        best_v, best_e = None, None
        vals_end = {}
        for k in range(NPT):
            e = lo + k * step
            v = M.mp_qprime(e, x, 50)
            if best_v is None or v < best_v:
                best_v, best_e = v, e
            if k in (0, NPT - 1):
                vals_end[k] = v
        dense[f"I{idx+1}"] = {
            "n_points": NPT, "dps": 50,
            "argmin_eps": dec30(best_e), "argmin_eps_float": float(best_e),
            "min_qprime": dec30(best_v), "min_qprime_float": float(best_v),
            "qprime_at_left_endpoint": dec30(vals_end[0]),
            "qprime_at_right_endpoint": dec30(vals_end[NPT - 1]),
            "argmin_is_right_endpoint": bool(abs(best_e - hi) < abs(step) * 1e-9),
            "argmin_within_one_grid_step_of_right_endpoint": bool(hi - best_e < step),
        }
        print(f"  dense scan I{idx+1}: n={NPT} argmin_eps={float(best_e):.15g} "
              f"min q'={float(best_v):.12e}  (left {float(vals_end[0]):.12e}, "
              f"right {float(vals_end[NPT-1]):.12e})")
        print(f"      infimum at right endpoint: "
              f"{dense[f'I{idx+1}']['argmin_within_one_grid_step_of_right_endpoint']}")
    sink.put("dense_scan_mp_dps50", dense)

    # ------------------------------------------------- q'' sign check (dps=35, intervals)
    qpp_check = {}
    iv.dps = 35
    for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
        worst_upper, worst_box, all_neg, argmax_up, argmin_lo = None, None, True, None, None
        for bx in boxes_by_i[idx]:
            up, lo = L.iv_bounds(L.qpp(bx, x, iv))[::-1]
            if up >= 0:
                all_neg = False
            if worst_upper is None or up > worst_upper:
                worst_upper, worst_box = up, L.iv_bounds(bx)
            argmax_up = up if argmax_up is None else max(argmax_up, up)
            argmin_lo = lo if argmin_lo is None else min(argmin_lo, lo)
        qpp_check[f"I{idx+1}"] = {
            "n_boxes": len(boxes_by_i[idx]), "iv_dps": 35,
            "max_qpp_upper": dec30(worst_upper), "max_qpp_upper_float": float(worst_upper),
            "argmax_box": [dec30(worst_box[0]), dec30(worst_box[1])],
            "min_qpp_lower_float": float(argmin_lo),
            "qpp_strictly_negative_on_whole_interval": bool(all_neg),
            "conclusion": ("q'' < 0 everywhere on the interval => q' is strictly "
                           "decreasing in eps => inf q' is attained at the RIGHT endpoint"
                           if all_neg else "NOT PROVEN: some box has q''_upper >= 0"),
        }
        print(f"  q'' sign I{idx+1}: max q''_upper = {float(worst_upper):+.6e}  "
              f"strictly negative everywhere: {all_neg}")
        print(f"      => {qpp_check[f'I{idx+1}']['conclusion']}")
    sink.put("qpp_sign_check_dps35", qpp_check)

    # ------------------------------------------- q'' point verification vs finite diff
    mp.dps = 50
    pts = []
    for idx, (ln, ld, hn, hd, n, x) in enumerate(INTERVALS):
        lo, hi = float(ln) / float(ld), float(hn) / float(hd)
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            pts.append((idx, lo + frac * (hi - lo), x))
    ver = []
    iv.dps = 35
    for idx, e, x in pts:
        em = mp.mpf(repr(e))
        a = None
        # analytic scalar q'' by the same chain rule (mp at dps=50)
        m = mp.sqrt(em); xx = mp.mpf(x); z = m * xx
        sz, cz = mp.sin(z), mp.cos(z)
        sx, cx = mp.sin(xx), mp.cos(xx)
        j1z = sz / z**2 - cz / z
        Djz = sz / z + cz / z**2 - sz / z**3
        j1pz = sz / z + 2 * cz / z**2 - 2 * sz / z**3
        Djpz = cz / z - 2 * sz / z**2 - 3 * cz / z**3 + 3 * sz / z**4
        Djppz = -sz / z - 3 * cz / z**2 + 7 * sz / z**3 + 12 * cz / z**4 - 12 * sz / z**5
        j1ppz = -(2 / z) * j1pz - (1 - 2 / z**2) * j1z
        j1x = sx / xx**2 - cx / xx
        y1x = -cx / xx**2 - sx / xx
        Djx = sx / xx + cx / xx**2 - sx / xx**3
        Dyx = -cx / xx + sx / xx**2 + cx / xx**3
        u = 1 / (2 * m); zp = xx * u; zpp = -xx / (4 * m**3)
        N = m * j1z * Djx - j1x * Djz
        Den = m * j1z * Dyx - y1x * Djz
        Ap = u * j1z + (xx / 2) * j1pz
        App = (-1 / (4 * m**3)) * j1z + xx * u**2 * j1pz + (xx**2 / 2) * u * j1ppz
        d2D = Djppz * zp**2 + Djpz * zpp
        Np = Ap * Djx - j1x * Djpz * zp
        Npp = App * Djx - j1x * d2D
        Denp = Ap * Dyx - y1x * Djpz * zp
        Denpp = App * Dyx - y1x * d2D
        F = Np * Den - N * Denp
        a = ((Npp * Den - N * Denpp) * Den - 2 * F * Denp) / Den**3
        # Richardson-extrapolated central first difference of q' estimates q''.
        # The step dependence is reported because a single unlucky step (h=1e-3)
        # can look ~1e-15 off even though the analytic value is exact to ~1e-29.
        fds = {}
        for hs, hh in (("1e-2", mp.mpf("1e-2")), ("1e-3", mp.mpf("1e-3")),
                       ("1e-4", mp.mpf("1e-4"))):
            d1 = (M.mp_qprime(em + hh, x, 60) - M.mp_qprime(em - hh, x, 60)) / (2 * hh)
            d2 = (M.mp_qprime(em + hh / 2, x, 60) - M.mp_qprime(em - hh / 2, x, 60)) / hh
            fds[hs] = (4 * d2 - d1) / 3
        rich = min(fds.values(), key=lambda v: abs(v - a))
        # interval containment: enclosing box at dps=35
        iv.dps = 35
        bx = iv.mpf([repr(e), repr(e)])
        q2i = L.qpp(bx, x, iv)
        i_lo, i_hi = L.iv_bounds(q2i)
        ver.append({
            "interval_index": idx + 1, "eps": repr(e), "x": x,
            "qpp_analytic_dps50": dec30(a),
            "qpp_richardson_fd_by_step": {k: dec30(v) for k, v in fds.items()},
            "qpp_richardson_fd_best": dec30(rich),
            "rel_diff_analytic_vs_best_fd": float(abs(a - rich) / abs(rich)),
            "rel_diff_by_step": {k: float(abs(a - v) / abs(a)) for k, v in fds.items()},
            "qpp_interval_dps35": [dec30(i_lo), dec30(i_hi)],
            "interval_contains_analytic": bool(i_lo <= a <= i_hi),
            "interval_width": dec30(i_hi - i_lo),
            "analytic_vs_interval_rel": float(abs(a - i_lo) / abs(a))})
    print("\n  q'' point verification (analytic vs Richardson FD of q', and interval containment):")
    print("    max rel diff analytic-vs-best-FD : %.3e" %
          max(v["rel_diff_analytic_vs_best_fd"] for v in ver))
    print("    rel diff by FD step: %s" % ver[0]["rel_diff_by_step"])
    print("    all intervals contain analytic value: %s" %
          all(v["interval_contains_analytic"] for v in ver))
    sink.put("qpp_point_verification", ver)

    # ------------------------------------------------------------------- summary
    summary = {
        "supplied_qprime_lower": SUPPLIED_QPRIME_LOWER,
        "V1_min_lower": [results["I1"]["V1"]["min_lower"], results["I2"]["V1"]["min_lower"]],
        "variant_range_I1": rng["I1"]["range_max_minus_min"],
        "variant_range_I2": rng["I2"]["range_max_minus_min"],
        "variants_matching_supplied_I1": rng["I1"]["best_matching_variants"],
        "variants_matching_supplied_I2": rng["I2"]["best_matching_variants"],
        "dense_argmin_eps": [dense["I1"]["argmin_eps_float"], dense["I2"]["argmin_eps_float"]],
        "dense_min_qprime": [dense["I1"]["min_qprime_float"], dense["I2"]["min_qprime_float"]],
        "qpp_strictly_negative": [qpp_check["I1"]["qpp_strictly_negative_on_whole_interval"],
                                  qpp_check["I2"]["qpp_strictly_negative_on_whole_interval"]],
        "wall_seconds": time.time() - t0,
    }
    sink.put("summary", summary)
    print("\n" + "=" * 100)
    print("SUMMARY")
    for k, v in summary.items():
        print(f"  {k} = {v}")
    print("wrote", os.path.relpath(out_path, ROOT))
    print("=" * 100)


if __name__ == "__main__":
    main()
