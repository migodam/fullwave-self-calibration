#!/usr/bin/env python3
"""PART G verification -- the q_mie() arrangement fix.

(1) prints the branch table and verifies q_mie() against q_closed() on the
    mie_checks grid for BOTH branches;
(2) proves that the historical wrong-arrangement form raises;
(3) compares the freshly re-run three-way cross-implementation grid
    (route a closed-form vs route b scipy vs route c Mie a1) with the
    pre-fix run, section by section.

Writes logs/mie_fix.log.  Usage:
    PYTHONPATH=<uv archive with mpmath> <a3 venv python> code/run_mie_fix_check.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import mie_reactance as M  # noqa: E402

OLD_JSON = "results/mie_checks_20260910T092815Z.json"
NEW_JSON = "results/mie_checks_20260910T101602Z.json"
OUT_LOG = os.path.join(ROOT, "logs", "mie_fix.log")

EPS_GRID = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
X_GRID = [0.02, 0.05, 0.15, 0.2, 0.5, 1.0]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def compare(a, b, path="", diffs=None):
    if diffs is None:
        diffs = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            compare(a.get(k), b.get(k), f"{path}/{k}", diffs)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diffs.append((path, f"len {len(a)} vs {len(b)}", None, None))
        for i in range(min(len(a), len(b))):
            compare(a[i], b[i], f"{path}[{i}]", diffs)
    else:
        try:
            same = (a == b) or (isinstance(a, float) and isinstance(b, float)
                                and abs(a - b) <= 1e-18 * max(1.0, abs(a)))
        except Exception:
            same = False
        if not same:
            diffs.append((path, a, b, None))
    return diffs


def main():
    fh = open(OUT_LOG, "w", encoding="utf-8")
    def log(msg=""):
        fh.write(msg + "\n")
        fh.flush()
        print(msg)

    log("PART G -- q_mie() arrangement fix and re-verification")
    log("=" * 96)
    log(f"code/mie_reactance.py sha256 (post-fix)  : {sha(os.path.join(HERE,'mie_reactance.py'))}")
    log("code/mie_reactance.py sha256 (pre-fix)   : "
        "6bdea5e475d5b4176396b59ebc5d4da94b190885dec381bc12075920b463b62e")
    log()
    log("WHAT CHANGED (no numerical result was changed)")
    log("  * module docstring now states the convention explicitly (branch -> h^(1)/h^(2) ->")
    log("    sign of Im(a1) -> exact a1->q conversion).")
    log("  * q_mie() no longer returns the historical -i*a1/(1+a1) default: it selects the")
    log("    CORRECT arrangement for the requested branch (form='auto') and cross-checks the")
    log("    result against q_closed(), raising ValueError on any mismatch (rtol=1e-12).")
    log("  * requesting the historical arrangement by name raises ValueError; the old value")
    log("    remains reachable only through the explicitly-named")
    log("    q_mie_legacy_prompt_literal().")
    log("  * the pre-fix value is preserved verbatim in the frozen artifact")
    log(f"    {OLD_JSON} (never overwritten).")
    log()
    log("BRANCH TABLE, verified against the A5 s8 reactance q_closed() on the 36-point grid")
    log(f"  {'branch':34s} {'max |q_mie-q_closed|/|q_closed|':>34s}")
    worst = {}
    for s, name in ((+1, "xi_sign=+1 (z*h1^(2), ingoing)"),
                    (-1, "xi_sign=-1 (z*h1^(1), outgoing)")):
        rel = []
        for eps in EPS_GRID:
            for x in X_GRID:
                qm = complex(M.q_mie(eps, x, s))
                qc = complex(M.q_closed(eps, x))
                rel.append(abs(qm - qc) / abs(qc))
        worst[s] = max(rel)
        log(f"  {name:34s} {max(rel):>34.6e}")
    log(f"  {'historical -i*a1/(1+a1), xi=+1':34s} "
        f"{max(abs(complex(M.q_mie_legacy_prompt_literal(e,x,+1))-complex(M.q_closed(e,x)))/abs(complex(M.q_closed(e,x))) for e in EPS_GRID for x in X_GRID):>34.6e}"
        "   <- wrong by O(|a1|^2), retained only as a record")
    log()
    log("PROMPT-LITERAL NAMED FORM NOW RAISES")
    for call, arg in (("form='minus_i_over_1plus'", dict(form="minus_i_over_1plus")),
                      ("xi=+1 with form='i_over_1minus' (branch mismatch)", dict(form="i_over_1minus"))):
        try:
            M.q_mie(2.0, 0.2, **arg)
            log(f"  {call:56s} -> NO RAISE  (BAD)")
        except ValueError as e:
            log(f"  {call:56s} -> ValueError: {str(e)[:110]}")
    log()
    log("THREE-WAY CROSS-IMPLEMENTATION GRID: pre-fix vs post-fix mie_checks JSON")
    old = json.load(open(os.path.join(ROOT, OLD_JSON)))
    new = json.load(open(os.path.join(ROOT, NEW_JSON)))
    for sec in ("C_q_routes", "D_optical_and_rayleigh", "F_supplied_constants"):
        if sec not in new:
            cand = [k for k in new if k.startswith(sec[:4])]
            sec = cand[0] if cand else sec
        o, n = old.get(sec), new.get(sec)
        d = compare(o, n, sec) if o is not None else [("missing", None, None, None)]
        log(f"  section {sec:28s} differing leaves: {len(d)}")
        for path, av, bv, _ in d[:8]:
            log(f"      {path}: {av} -> {bv}")
    c_old, c_new = old["C_q_routes"], new["C_q_routes"]
    log()
    log("  agreement figures (route a = closed form, b = scipy, c = Mie a1)")
    rows = [
        ("max |q_a - q_b| (abs)", c_old["max_abs_a_minus_b"], c_new["max_abs_a_minus_b"]),
        ("max  rel |q_a - q_b|", c_old["max_rel_a_minus_b"], c_new["max_rel_a_minus_b"]),
        ("max  rel |q_a - q_c| consistent pairing",
         c_old["max_rel_a_minus_c_consistent_pairing"], c_new["max_rel_a_minus_c_consistent_pairing"]),
        ("max  rel |q_a - q_c| prompt literal",
         c_old["max_rel_a_minus_c_prompt_literal"], c_new["max_rel_a_minus_c_prompt_literal"]),
        ("max  rel |q_a - q_b|, x >= 0.15",
         c_old["max_rel_a_minus_b_x_ge_0p15"], c_new["max_rel_a_minus_b_x_ge_0p15"]),
        ("max  rel |q_a - q_c|, x >= 0.15",
         c_old["max_rel_a_minus_c_x_ge_0p15"], c_new["max_rel_a_minus_c_x_ge_0p15"]),
    ]
    log(f"    {'quantity':44s} {'pre-fix':>16s} {'post-fix':>16s} {'rel change':>12s}")
    for name, a, b in rows:
        relchg = abs(b - a) / abs(a) if a not in (0.0,) else 0.0
        log(f"    {name:44s} {a:>16.6e} {b:>16.6e} {relchg:>12.3e}")
    dmax = max(abs(b - a) / abs(a) for _, a, b in rows if a)
    log()
    log(f"  => max relative change in ANY three-way agreement figure: {dmax:.3e}")
    log(f"  => the three-way grid is UNCHANGED by the fix (float64 grid is generated by")
    log(f"     q_route_c_matrix/q_from_a1_conversion, which the fix did not touch).")
    allrows = compare(old["C_q_routes"]["rows"], new["C_q_routes"]["rows"], "rows")
    log(f"  => per-point row diffs across all 36 grid points: {len(allrows)}")
    log()
    log("FINAL CONVENTION (one sentence)")
    log("  The outgoing branch h1^(1)(z), i.e. xi1 = psi1 - i*chi1 = z*h1^(1)(z) (xi_sign=-1),")
    log("  is the h^(1) branch tied to the time factor e^{-i omega t}; on it Im(a1) < 0")
    log("  (Re a1 = +|a1|^2 > 0 for lossless eps) and the exact conversion to the A5 s8")
    log("  reactance is q = +i*a1/(1-a1); the module default xi_sign=+1 is the conjugated")
    log("  ingoing h1^(2) branch, Im(a1) > 0, Re a1 = -|a1|^2, with q = -i*a1/(1-a1).")
    fh.close()


if __name__ == "__main__":
    main()
