#!/usr/bin/env python3
"""Deliverable 2: run every Deliverable-1 check on small grids.

Writes results/mie_checks_<UTCSTAMP>.json (written once, never edited) and tees the
full log to logs/mie_checks.log.

Usage:
    PYTHONPATH=<uv archive with mpmath> <a3_research venv python> code/run_mie_checks.py
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import mie_reactance as M  # noqa: E402

TEE_STREAM = None


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


def hash_many(rel_paths):
    out = {}
    for rp in rel_paths:
        p = os.path.join(ROOT, rp)
        out[rp] = sha256_file(p) if os.path.exists(p) else None
    return out


def rel_err(a, b):
    return abs(a - b) / abs(b) if abs(b) > 0 else abs(a - b)


def main():
    global TEE_STREAM
    os.makedirs(os.path.join(ROOT, "logs"), exist_ok=True)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    TEE_STREAM = Tee(os.path.join(ROOT, "logs", "mie_checks.log"))
    sys.stdout = TEE_STREAM
    sys.stderr = TEE_STREAM

    t_start = time.time()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_path = os.path.join(ROOT, "results", f"mie_checks_{stamp}.json")
    if os.path.exists(out_path):
        raise SystemExit(f"REFUSING to overwrite existing results file: {out_path}")

    print("=" * 100)
    print("A5 Section 8 reactance q -- Deliverable 1 cross-implementation checks")
    print("started_utc :", datetime.now(timezone.utc).isoformat())
    print("command     :", " ".join(shlex.quote(a) for a in sys.argv))
    print("interpreter :", sys.executable)
    print("python      :", sys.version.replace("\n", " "))
    print("platform    :", platform.platform())
    print("cwd         :", os.getcwd())
    print("=" * 100)

    import scipy
    versions = {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__}
    try:
        import mpmath
        versions["mpmath"] = mpmath.__version__
        versions["mpmath_file"] = mpmath.__file__
    except ImportError as exc:  # pragma: no cover
        versions["mpmath"] = f"NOT IMPORTABLE ({exc})"
        versions["mpmath_file"] = None
    print("versions    :", versions)

    input_files = [
        "inputs/WORKER_TASK.md", "inputs/A5_THEORY.md", "inputs/modal_boundary.json",
        "inputs/SOURCE_HASHES.json", "inputs/derive/A5_sec_05.md", "inputs/derive/A5_sec_07.md",
        "inputs/derive/A5_sec_08.md", "inputs/derive/MODAL_BOUNDARY_KEYS.txt",
        "logs/uv_cache_VjSSQm31390egyfy.sha256",
    ]
    code_files = ["code/mie_reactance.py", "code/run_mie_checks.py"]
    hashes = {"inputs": hash_many(input_files), "code": hash_many(code_files)}
    for k, v in hashes["inputs"].items():
        print(f"  sha256(inputs) {v}  {k}")
    for k, v in hashes["code"].items():
        print(f"  sha256(code)   {v}  {k}")

    results = {
        "meta": {
            "task": "A5 section 8 exact spherical electric-dipole reactance; Deliverable 1 checks",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command_line": " ".join(shlex.quote(a) for a in sys.argv),
            "interpreter": sys.executable,
            "versions": versions,
            "cwd": os.getcwd(),
            "file_sha256": hashes,
        }
    }

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[A] closed-form trigonometric expressions vs their definitions / scipy")
    print("-" * 100)
    zc = np.array([1e-3, 0.02, 0.05, 0.1, 0.2, 0.2449, 0.5, 1.0, 2.0, 3.5, 5.0, 8.0])
    a_res = {
        "z_grid": zc.tolist(),
        "j1_closed_vs_scipy_maxabs": float(np.max(np.abs(M.j1_closed(zc) - M.special.spherical_jn(1, zc)))),
        "y1_closed_vs_scipy_maxabs": float(np.max(np.abs(M.y1_closed(zc) - M.special.spherical_yn(1, zc)))),
        "j1p_closed_vs_scipy_derivative_kwarg_maxabs": float(
            np.max(np.abs(M.j1p_closed(zc) - M.special.spherical_jn(1, zc, derivative=True)))),
        "y1p_closed_vs_scipy_derivative_kwarg_maxabs": float(
            np.max(np.abs(M.y1p_closed(zc) - M.special.spherical_yn(1, zc, derivative=True)))),
        "j1p_vs_j0_minus_2j1_over_z_maxabs": float(np.max(np.abs(
            M.special.spherical_jn(0, zc) - 2 * M.special.spherical_jn(1, zc) / zc
            - M.special.spherical_jn(1, zc, derivative=True)))),
        "Dj_closed_vs_definition_maxabs": float(np.max(np.abs(M.Dj_closed(zc) - M.Dj_def(zc)))),
        "Dy_closed_vs_definition_maxabs": float(np.max(np.abs(M.Dy_closed(zc) - M.Dy_def(zc)))),
        "j1_closed_vs_scipy_maxrel": float(np.max(np.abs(M.j1_closed(zc) - M.special.spherical_jn(1, zc))
                                                   / np.abs(M.special.spherical_jn(1, zc)))),
        "Dj_closed_vs_definition_maxrel": float(np.max(np.abs(M.Dj_closed(zc) - M.Dj_def(zc)) / np.abs(M.Dj_def(zc)))),
        "Dy_closed_vs_definition_maxrel": float(np.max(np.abs(M.Dy_closed(zc) - M.Dy_def(zc)) / np.abs(M.Dy_def(zc)))),
    }
    for k, v in a_res.items():
        if k != "z_grid":
            print(f"  {k:52s} {v:.6e}")
    results["A_closed_forms"] = a_res

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[B] D_j'(z) and D_y'(z) analytic forms vs mpmath 50-dps central differences (h=1e-20)")
    print("-" * 100)
    from mpmath import mp
    mp.dps = 50
    h = mp.mpf("1e-20")
    zfd = [mp.mpf(s) for s in ["0.02", "0.05", "0.1", "0.15", "0.2", "0.2449", "0.3", "0.5",
                               "0.75", "1.0", "1.5", "2.0", "3.0", "5.0", "8.0", "12.0"]]
    b_rows = []
    for z in zfd:
        mp.dps = 50
        d = M.mp_derivatives(z, dps=50)
        mp.dps = 60
        dp = M.mp_derivatives(z + h, dps=60)
        dm = M.mp_derivatives(z - h, dps=60)
        mp.dps = 50
        fd_Djp = (dp["Dj"] - dm["Dj"]) / (2 * h)
        fd_Dyp = (dp["Dy"] - dm["Dy"]) / (2 * h)
        b_rows.append({
            "z": str(z),
            "Dj_prime_analytic": str(d["Djp"]),
            "Dj_prime_fd": str(fd_Djp),
            "Dj_prime_rel_err": str(abs(fd_Djp - d["Djp"]) / abs(d["Djp"])),
            "Dy_prime_analytic": str(d["Dyp"]),
            "Dy_prime_fd": str(fd_Dyp),
            "Dy_prime_rel_err": str(abs(fd_Dyp - d["Dyp"]) / abs(d["Dyp"])),
        })
    b_res = {
        "fd_step": "1e-20",
        "mpmath_dps": 50,
        "max_rel_err_Dj_prime": str(max(mp.mpf(r["Dj_prime_rel_err"]) for r in b_rows)),
        "max_rel_err_Dy_prime": str(max(mp.mpf(r["Dy_prime_rel_err"]) for r in b_rows)),
        "rows": b_rows,
    }
    print(f"  max |rel err| D_j'(z) : {mp.nstr(mp.mpf(b_res['max_rel_err_Dj_prime']), 6)}")
    print(f"  max |rel err| D_y'(z) : {mp.nstr(mp.mpf(b_res['max_rel_err_Dy_prime']), 6)}")
    results["B_derivative_fd"] = b_res

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[C] q via three independent routes on the required grid")
    print("-" * 100)
    eps_grid = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
    x_grid = [0.02, 0.05, 0.15, 0.2, 0.5, 1.0]
    rows = []
    max_ab_abs = max_ab_rel = 0.0
    for eps in eps_grid:
        for x in x_grid:
            qa = float(M.q_closed(eps, x))
            qb = float(M.q_scipy(eps, x))
            variants = {f"{s:+d}|{f}": complex(v) for (s, f), v in M.q_route_c_matrix(eps, x).items()}
            qc = variants["-1|i_over_1minus"].real  # consistent pairing (see DERIVATIONS.md)
            qc_prompt = variants["-1|minus_i_over_1plus"]
            rows.append({
                "eps": eps, "x": x,
                "q_a_closed": qa, "q_b_scipy": qb,
                "q_c_a1_i_over_1minus_xi_minus": qc,
                "q_c_prompt_literal_xi_minus": [qc_prompt.real, qc_prompt.imag],
                "abs_a_minus_b": abs(qa - qb),
                "rel_a_minus_b": rel_err(qa, qb),
                "rel_a_minus_c": rel_err(qc, qa),
                "rel_a_minus_c_prompt": rel_err(abs(qc_prompt), qa),
            })
            max_ab_abs = max(max_ab_abs, abs(qa - qb))
            max_ab_rel = max(max_ab_rel, rel_err(qa, qb))
    rel_ac = max(r["rel_a_minus_c"] for r in rows)
    rel_ac_prompt = max(r["rel_a_minus_c_prompt"] for r in rows)
    rows_big_x = [r for r in rows if r["x"] >= 0.15]
    ab_big_rel = max(r["rel_a_minus_b"] for r in rows_big_x)
    ac_big_rel = max(r["rel_a_minus_c"] for r in rows_big_x)
    w_ab = max(rows, key=lambda r: r["rel_a_minus_b"])
    w_ac = max(rows, key=lambda r: r["rel_a_minus_c"])
    w_prompt = max(rows, key=lambda r: r["rel_a_minus_c_prompt"])
    c_res = {
        "eps_grid": eps_grid, "x_grid": x_grid, "n_points": len(rows),
        "max_abs_a_minus_b": max_ab_abs, "max_rel_a_minus_b": max_ab_rel,
        "max_rel_a_minus_c_consistent_pairing": rel_ac,
        "max_rel_a_minus_c_prompt_literal": rel_ac_prompt,
        "max_rel_a_minus_b_x_ge_0p15": ab_big_rel,
        "max_rel_a_minus_c_x_ge_0p15": ac_big_rel,
        "worst_a_minus_b": {"eps": w_ab["eps"], "x": w_ab["x"], "rel": w_ab["rel_a_minus_b"],
                            "q_a": w_ab["q_a_closed"], "q_b": w_ab["q_b_scipy"]},
        "worst_a_minus_c": {"eps": w_ac["eps"], "x": w_ac["x"], "rel": w_ac["rel_a_minus_c"],
                            "q_a": w_ac["q_a_closed"], "q_c": w_ac["q_c_a1_i_over_1minus_xi_minus"]},
        "rows": rows,
    }
    print(f"  grid: {len(rows)} points")
    print(f"  max |q_a - q_b|          abs {max_ab_abs:.3e}   rel {max_ab_rel:.3e}")
    print(f"  max  rel |q_a - q_c|     (xi=-1, q=i*a1/(1-a1))        {rel_ac:.3e}")
    print(f"  max  rel |q_a - q_c|     (prompt literal xi=-1, -i*a1/(1+a1)) {rel_ac_prompt:.3e}")
    print(f"  worst |q_a-q_b| at eps={w_ab['eps']} x={w_ab['x']} (q_a={w_ab['q_a_closed']:.6e})")
    print(f"  worst |q_a-q_c| at eps={w_ac['eps']} x={w_ac['x']} (q_a={w_ac['q_a_closed']:.6e}, "
          f"q_c={w_ac['q_c_a1_i_over_1minus_xi_minus']:.6e})")
    print(f"  restricted to x>=0.15: max rel |q_a-q_b| {ab_big_rel:.3e} , max rel |q_a-q_c| {ac_big_rel:.3e}")
    print(f"  prompt-literal worst point eps={w_prompt['eps']} x={w_prompt['x']}: "
          f"q_c={w_prompt['q_c_prompt_literal_xi_minus']}")
    results["C_q_routes"] = c_res

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[D] optical theorem and Rayleigh limit")
    print("-" * 100)
    d_res = {"optical_rows": []}
    worst_std = worst_t = 0.0
    for eps in eps_grid:
        for x in x_grid:
            a1 = complex(M.mie_a1(eps, x, xi_sign=-1))     # standard Mie coefficient, xi = psi - i*chi
            t = -a1                                        # A5 section 5 T-matrix element
            e_std = abs(a1.real - abs(a1) ** 2)
            e_t = abs(t.real + abs(t) ** 2)
            worst_std = max(worst_std, e_std)
            worst_t = max(worst_t, e_t)
            d_res["optical_rows"].append({"eps": eps, "x": x, "a1": [a1.real, a1.imag],
                                          "Re_a1_minus_abs_a1_sq": a1.real - abs(a1) ** 2,
                                          "Re_t_plus_abs_t_sq": t.real + abs(t) ** 2})
    d_res["max_abs_Re_a1_minus_abs2"] = worst_std
    d_res["max_abs_Re_t_plus_abs2"] = worst_t
    print(f"  max |Re(a1) - |a1|^2|          (standard a1, xi=psi-i*chi) : {worst_std:.6e}")
    print(f"  max |Re(t)  + |t|^2|           (A5 t = -a1)                 : {worst_t:.6e}")

    ray = {"rows": [], "note": "high-precision ratio uses mpmath at 50 dps; the double-precision "
           "route loses digits for x <= 1e-5 because q ~ x^3 is a cancellation-limited difference."}
    xseq = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
    for eps in eps_grid:
        A = (2.0 / 3.0) * (eps - 1.0) / (eps + 2.0)
        r = []
        for x in xseq:
            mp.dps = 50
            q_mp = M.mp_q(mp.mpf(eps), mp.mpf(x), dps=50)
            t_mp = 1j * q_mp / (1 - 1j * q_mp)
            q_dbl = float(M.q_closed(eps, x))
            r.append({"x": x,
                      "q_over_q_rayleigh_mp50": str(q_mp / (mp.mpf(A) * mp.mpf(x) ** 3)),
                      "Im_t_over_rayleigh_tE_mp50": str(t_mp.imag / (mp.mpf(A) * mp.mpf(x) ** 3)),
                      "Re_t_over_minus_rayleigh_tE_sq_mp50": str(
                          t_mp.real / (-(mp.mpf(A) * mp.mpf(x) ** 3) ** 2)) if A != 0 else None,
                      "q_over_q_rayleigh_double": q_dbl / (A * x ** 3) if A != 0 else None})
        ray["rows"].append({"eps": eps, "A": A, "samples": r})
        last = r[-1]
        print(f"  eps={eps}: x=1e-7 mp50  q/q_Rayleigh={float(last['q_over_q_rayleigh_mp50']):.12f}"
              f"   Im t/Rayleigh t_E={float(last['Im_t_over_rayleigh_tE_mp50']):.12f}"
              f"   (double at x=1e-7: {last['q_over_q_rayleigh_double']:.6f})")
    d_res["rayleigh"] = ray
    results["D_optical_and_rayleigh"] = d_res

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[E] q' and h' analytic derivatives vs mpmath 50-dps central differences in eps (h=1e-20)")
    print("-" * 100)
    e_rows = []
    worst_qp = worst_hp = 0.0
    for eps in eps_grid:
        for x in x_grid:
            mp.dps = 50
            vals = M.mp_quantities(mp.mpf(eps), mp.mpf(x), dps=50)
            mp.dps = 60
            vp = M.mp_quantities(mp.mpf(eps) + h, mp.mpf(x), dps=60)
            vm = M.mp_quantities(mp.mpf(eps) - h, mp.mpf(x), dps=60)
            mp.dps = 50
            fd_qp = (vp["q"] - vm["q"]) / (2 * h)
            fd_hp = ((vp["q"] / mp.sqrt(1 + vp["q"] ** 2))
                     - (vm["q"] / mp.sqrt(1 + vm["q"] ** 2))) / (2 * h)
            rqp = abs(fd_qp - vals["qprime"]) / abs(vals["qprime"])
            rhp = abs(fd_hp - vals["hprime"]) / abs(vals["hprime"])
            worst_qp = max(worst_qp, rqp)
            worst_hp = max(worst_hp, rhp)
            e_rows.append({"eps": eps, "x": x,
                           "q": str(vals["q"]), "qprime": str(vals["qprime"]),
                           "Den": str(vals["Den"]), "hprime": str(vals["hprime"]),
                           "qprime_fd_rel_err": str(rqp), "hprime_fd_rel_err": str(rhp),
                           "q_double_vs_mp50_rel_err": str(
                               abs(float(M.q_closed(eps, x)) - float(vals["q"])) / abs(float(vals["q"])))})
    e_res = {
        "fd_step": "1e-20", "mpmath_dps": 50,
        "max_rel_err_qprime": str(worst_qp), "max_rel_err_hprime": str(worst_hp),
        "max_rel_err_q_double_vs_mp50": str(max(mp.mpf(r["q_double_vs_mp50_rel_err"]) for r in e_rows)),
        "rows": e_rows,
    }
    print(f"  max |rel err| q'(eps,x)  : {mp.nstr(worst_qp, 6)}")
    print(f"  max |rel err| h'(eps,x)  : {mp.nstr(worst_hp, 6)}")
    print(f"  max rel err q double vs mp50 : {mp.nstr(mp.mpf(e_res['max_rel_err_q_double_vs_mp50']), 6)}")
    results["E_qprime_hprime_fd"] = e_res

    # ---------------------------------------------------------------------------------
    print("\n" + "-" * 100)
    print("[F] reproducibility of the supplied JSON constants (world q + exact-reactance residual)")
    print("-" * 100)
    J = json.load(open(os.path.join(ROOT, "inputs", "modal_boundary.json")))
    f_res = {"world_checks": []}
    for i, (eps, x) in enumerate([(J["world1_eps"][0], 0.2), (J["world1_eps"][1], 0.15),
                                  (J["world2_eps"][0], 0.2), (J["world2_eps"][1], 0.15)]):
        q = float(M.q_closed(eps, x))
        ref = (J["world1_q"] + J["world2_q"])[i]
        f_res["world_checks"].append({"eps": eps, "x": x, "q_computed": q, "q_supplied": ref,
                                      "abs_diff": abs(q - ref), "rel_diff": rel_err(q, ref)})
        print(f"  world q[{i}] eps={eps:.12f} x={x}: computed {q:.16e} supplied {ref:.16e} "
              f"rel {rel_err(q, ref):.3e}")
    cands = {}
    for i, (eps, x) in enumerate([(2.0, 0.2), (3.0, 0.15), (2.359740580371098, 0.2),
                                  (3.989319415017825, 0.15)]):
        qa = float(M.q_closed(eps, x))
        qc = complex(M.q_from_a1_conversion(eps, x, xi_sign=-1, form="i_over_1minus"))
        cands[f"|q_a-q_c| eps={eps} x={x}"] = abs(qa - qc.real)
    cands["max_over_world_points"] = max(cands.values())
    f_res["exact_reactance_identity_error_supplied"] = J["exact_reactance_identity_error"]
    f_res["identity_residual_candidates"] = {k: float(v) for k, v in cands.items()}
    print(f"  supplied exact_reactance_identity_error = {J['exact_reactance_identity_error']:.6e}")
    print(f"  |q_a - q_c| over world points           = {cands['max_over_world_points']:.6e}")
    results["F_supplied_constants"] = f_res

    # ---------------------------------------------------------------------------------
    results["meta"]["wall_seconds"] = time.time() - t_start
    with open(out_path, "x", encoding="utf-8") as f:
        json.dump(results, f, indent=2, sort_keys=False)
        f.write("\n")
    print("\n" + "=" * 100)
    print("wrote", os.path.relpath(out_path, ROOT))
    print(f"total wall time {results['meta']['wall_seconds']:.2f} s")
    print("=" * 100)
    TEE_STREAM.flush()


if __name__ == "__main__":
    main()
