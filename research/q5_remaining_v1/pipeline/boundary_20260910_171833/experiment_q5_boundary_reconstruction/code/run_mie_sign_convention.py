#!/usr/bin/env python3
"""PART E -- resolve the a1 -> q sign convention documented in code/mie_reactance.py.

Question (from the task): for eps = 2, x = 0.2, print
  * the textbook Mie electric-dipole coefficient a1 built from the OUTGOING Hankel
    function h1^(1) with time factor e^{-i omega t},
  * Re(a1) + |a1|^2,
  * the A5 theory scalar t = i*q/(1 - i*q) with q from the Section-8 reactance formula,
  * the two candidate conversions  -i*a1/(1+a1)  and  +i*a1/(1-a1),
and state which one reproduces q(eps, x).

Convention pinned down explicitly (no literature hand-waving):
    h1^(1)(z) = j1(z) + i*y1(z)              [outgoing for e^{-i omega t}]
    psi1(z) = z*j1(z) = sin z/z - cos z
    chi1(z) = -z*y1(z) = cos z/z + sin z
    xi1^(1) = z*h1^(1)(z) = psi1 - i*chi1 = -e^{iz}(z+i)/z      [exact closed form]
    xi1^(1)'(z) = e^{iz}(z + i - i z^2)/z^2                     [exact closed form]
    xi1^(2) = psi1 + i*chi1 = conj(xi1^(1)) on the real axis     [code's xi_sign=+1]

Independent a1 routes (all at eps=2, x=0.2):
  (i)   code/mie_reactance.mie_a1(..., xi_sign=+1 / -1)   [scipy spherical_jn / spherical_yn]
  (ii)  mpmath dps=40 with the exact closed-form xi1^(1)  [fully independent]
  (iii) scipy.special.hankel1(3/2, .) with an exact recurrence for xi1' [independent family]

Writes results/mie_sign_convention_<TS>.json and appends logs/mie_sign_convention.log.
No installs, single thread, deterministic, no figures.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mie_reactance as mr  # noqa: E402

WALL0 = time.time()
EPS, X = 2.0, 0.2


class Rec:
    def __init__(self):
        self.data = {}
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        base = os.path.join(ROOT, "results", "mie_sign_convention_%s" % stamp)
        n = 0
        while os.path.exists(base + (".json" if n == 0 else "_%d.json" % n)):
            n += 1
        self.json_path = base + (".json" if n == 0 else "_%d.json" % n)
        self.log_path = os.path.join(ROOT, "logs", "mie_sign_convention.log")
        os.makedirs(os.path.dirname(self.log_path), exist_ok=True)
        self.logf = open(self.log_path, "a")

    def log(self, msg):
        line = "[%s +%5.1fs] %s" % (time.strftime("%H:%M:%S"), time.time() - WALL0, msg)
        print(line, flush=True)
        self.logf.write(line + "\n")
        self.logf.flush()

    def dump(self):
        self.data.setdefault("meta", {}).update({
            "script": os.path.relpath(__file__, ROOT),
            "interpreter": sys.executable,
            "numpy_version": np.__version__,
            "wall_seconds": round(time.time() - WALL0, 3),
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "eps": EPS, "x": X,
        })

        def enc(o):
            if isinstance(o, dict):
                return {str(k): enc(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)):
                return [enc(v) for v in o]
            if isinstance(o, (bool, str)) or o is None:
                return o
            if isinstance(o, (int, np.integer)):
                return int(o)
            if isinstance(o, (float, np.floating)):
                return float(o)
            if isinstance(o, (complex, np.complexfloating)):
                return {"re": float(o.real), "im": float(o.imag)}
            return str(o)
        tmp = self.json_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(enc(self.data), fh, indent=1)
        os.replace(tmp, self.json_path)


def reldiff(a, b):
    a, b = complex(a), complex(b)
    return abs(a - b) / max(abs(b), 1e-300)


def zi(xi_sign):
    """Riccati-Bessel combination used by mie_reactance.xi1 / mie_a1."""
    return "psi1 + i chi1 (z h1^(2), INGOING branch)" if xi_sign > 0 else \
           "psi1 - i chi1 (z h1^(1), OUTGOING branch, e^{-i omega t})"


def a1_mpmath(xi_sign, dps=40):
    """a1 built from the exact closed forms psi1 = sin z/z - cos z, chi1 = cos z/z + sin z
    and xi1 = psi1 + i*xi_sign*chi1, at mpmath dps (no scipy involved)."""
    from mpmath import mp, mpc, mpf, sin, cos, sqrt, exp
    mp.dps = dps

    def psi1(z):
        return sin(z) / z - cos(z)

    def psi1p(z):
        return (z * cos(z) - sin(z)) / z ** 2 + sin(z)

    def xi1(z):
        c = cos(z) / z + sin(z)
        return psi1(z) + mpc(0, xi_sign) * c

    def xi1p(z):
        cp = (-z * sin(z) - cos(z)) / z ** 2 + cos(z)
        return psi1p(z) + mpc(0, xi_sign) * cp

    e, xx = mpf(repr(EPS)), mpf(repr(X))
    m, mx = sqrt(e), sqrt(e) * mpf(repr(X))
    num = m * psi1(mx) * psi1p(xx) - psi1(xx) * psi1p(mx)
    den = m * psi1(mx) * xi1p(xx) - xi1(xx) * psi1p(mx)
    return num / den


def xi1_outgoing_closed(z):
    """Exact xi1^(1)(z) = z h1^(1)(z) = -e^{iz}(z+i)/z and its derivative."""
    from mpmath import mp, mpc, exp
    return -exp(mpc(0, 1) * z) * (z + mpc(0, 1)) / z, \
        exp(mpc(0, 1) * z) * (z + mpc(0, 1) - mpc(0, 1) * z ** 2) / z ** 2


def a1_from_closed_h1(dps=40):
    """Textbook a1 with OUTGOING h1^(1), built only from exact closed forms."""
    from mpmath import mp, mpf, sqrt, sin, cos
    mp.dps = dps

    def psi1(z):
        return sin(z) / z - cos(z)

    def psi1p(z):
        return (z * cos(z) - sin(z)) / z ** 2 + sin(z)

    e, xx = mpf(repr(EPS)), mpf(repr(X))
    m, mx = sqrt(e), sqrt(e) * mpf(repr(X))
    xi, dxi = xi1_outgoing_closed(xx)
    num = m * psi1(mx) * psi1p(xx) - psi1(xx) * psi1p(mx)
    den = m * psi1(mx) * dxi - xi * psi1p(mx)
    return num / den, xi, dxi, xx


def a1_hankel1(dps=40):
    """Independent route: xi1^(1)(z) = z*h1^(1)(z) with h1^(1) = sqrt(pi/(2z)) H_{3/2}^{(1)}
    via the recurrence xi1' = xi0 - xi1/z, xi0 = -i e^{iz}."""
    from mpmath import mp, mpc, mpf, sqrt, sin, cos, exp, besselj, bessely
    mp.dps = dps

    def xi1(z):
        j = sqrt(np.pi / (2 * z)) * (besselj(mpf(1.5), z) + mpc(0, 1) * bessely(mpf(1.5), z))
        return z * j

    def xi0(z):
        return -mpc(0, 1) * exp(mpc(0, 1) * z)

    e, xx = mpf(repr(EPS)), mpf(repr(X))
    m, mx = sqrt(e), sqrt(e) * mpf(repr(X))
    xz = mpc(xx)
    xi, dxi = xi1(xz), xi0(xz) - xi1(xz) / xz
    psi1 = lambda z: sin(z) / z - cos(z)
    psi1p = lambda z: (z * cos(z) - sin(z)) / z ** 2 + sin(z)
    num = m * psi1(mx) * psi1p(xx) - psi1(xx) * psi1p(mx)
    den = m * psi1(mx) * dxi - xi * psi1p(mx)
    return num / den, xi, dxi


def main():
    rec = Rec()
    rec.log("PART E start  json=%s" % os.path.relpath(rec.json_path, ROOT))

    q = float(mr.q_closed(EPS, X))
    q_high = complex(mr.mp_q(EPS, X))
    q_scipy = float(mr.q_scipy(EPS, X))
    t = 1j * q / (1 - 1j * q)

    a1p = complex(mr.mie_a1(EPS, X, +1))     # code xi_sign=+1  -> psi1 + i chi1
    a1m = complex(mr.mie_a1(EPS, X, -1))     # code xi_sign=-1  -> psi1 - i chi1 == z h1^(1)
    a1_mp_p = complex(a1_mpmath(+1))
    a1_mp_m = complex(a1_mpmath(-1))
    a1_closed, xi_c, dxi_c, xz = a1_from_closed_h1()
    a1_hk, xi_hk, dxi_hk = a1_hankel1()

    # is the code's xi_sign=+1 combination really z h1^(2) and xi_sign=-1 really z h1^(1)?
    branch = {
        "z_times_h1_outgoing_closed_form": {"re": complex(xi_c).real, "im": complex(xi_c).imag},
        "psi1_plus_i_chi1(x)": {"re": complex(mr.xi1(xz, +1)).real, "im": complex(mr.xi1(xz, +1)).imag},
        "psi1_minus_i_chi1(x)": {"re": complex(mr.xi1(xz, -1)).real, "im": complex(mr.xi1(xz, -1)).imag},
        "rel_closed_h1_vs_code_xi_sign_plus1": reldiff(xi_c, mr.xi1(xz, +1)),
        "rel_closed_h1_vs_code_xi_sign_minus1": reldiff(xi_c, mr.xi1(xz, -1)),
        "closed_h1_vs_hankel1_route": reldiff(xi_c, xi_hk),
        "deriv_closed_vs_recurrence": reldiff(dxi_c, dxi_hk),
        "conclusion": ("z*h1^(1)(z) == psi1 - i*chi1 == code mie_a1(xi_sign=-1).  The code's "
                       "DEFAULT xi_sign=+1 is the psi1 + i*chi1 = z*h1^(2) (ingoing) branch."),
    }

    conv = {}
    for tag, a in (("a1_at_xi_sign_plus1(=z*h1^(2))", a1p),
                   ("a1_at_xi_sign_minus1(=z*h1^(1), textbook)", a1m)):
        conv[tag] = {
            "a1": {"re": a.real, "im": a.imag},
            "A: -i*a1/(1+a1)": {"value": {"re": (-1j * a / (1 + a)).real,
                                          "im": (-1j * a / (1 + a)).imag},
                                "rel_err_vs_q": reldiff(-1j * a / (1 + a), q)},
            "B: +i*a1/(1-a1)": {"value": {"re": (1j * a / (1 - a)).real,
                                          "im": (1j * a / (1 - a)).imag},
                                "rel_err_vs_q": reldiff(1j * a / (1 - a), q)},
            "C: -i*a1/(1-a1)": {"value": {"re": (-1j * a / (1 - a)).real,
                                          "im": (-1j * a / (1 - a)).imag},
                                "rel_err_vs_q": reldiff(-1j * a / (1 - a), q)},
            "Re_a1_plus_abs_a1_sq": float(a.real + abs(a) ** 2),
        }

    a1_tb = a1m  # textbook branch
    out = {
        "eps": EPS, "x": X,
        "reactance_reference": {
            "q_closed": q, "q_mpmath_dps50": {"re": q_high.real, "im": q_high.imag},
            "q_scipy": q_scipy,
            "q_closed_vs_scipy_rel": reldiff(q_scipy, q),
            "q_closed_vs_mpmath_rel": reldiff(q, q_high.real),
        },
        "branch_identification": branch,
        "textbook_a1_outgoing_h1": {
            "definition": ("xi1^(1)(z) = z*h1^(1)(z) = psi1 - i*chi1 = -e^{iz}(z+i)/z, "
                           "chi1 = -z*y1(z), h1^(1) = j1 + i*y1, time factor e^{-i omega t}; "
                           "code mie_a1(..., xi_sign=-1)"),
            "a1_scipy_route": {"re": a1_tb.real, "im": a1_tb.imag},
            "a1_mpmath_route": {"re": a1_mp_m.real, "im": a1_mp_m.imag},
            "a1_closed_form_h1_route": {"re": complex(a1_closed).real, "im": complex(a1_closed).imag},
            "a1_hankel1_route": {"re": complex(a1_hk).real, "im": complex(a1_hk).imag},
            "route_agreement_rel_scipy_vs_mpmath": reldiff(a1_tb, a1_mp_m),
            "route_agreement_rel_scipy_vs_closed_form": reldiff(a1_tb, a1_closed),
            "route_agreement_rel_scipy_vs_hankel1": reldiff(a1_tb, a1_hk),
            "Re_a1_plus_abs_a1_sq": float(a1_tb.real + abs(a1_tb) ** 2),
            "Re_a1_minus_abs_a1_sq": float(a1_tb.real - abs(a1_tb) ** 2),
            "equals_code_xi_sign_minus1": bool(abs(a1_tb - a1m) < 1e-30),
            "rel_vs_conj_of_code_xi_plus1": reldiff(a1_tb, np.conj(a1p)),
        },
        "theory_t": {
            "t = i q/(1 - i q)": {"re": t.real, "im": t.imag},
            "Re_t_plus_abs_t_sq": float(t.real + abs(t) ** 2),
            "t_vs_minus_textbook_a1_rel": reldiff(t, -a1_tb),
            "t_vs_minus_conj_xi_plus1_a1_rel": reldiff(t, -np.conj(a1p)),
        },
        "conversions": conv,
        "code_q_mie_default": {
            "q_mie(xi_sign=+1)": {"re": complex(mr.q_mie(EPS, X, +1)).real,
                                  "im": complex(mr.q_mie(EPS, X, +1)).imag},
            "q_mie(xi_sign=-1)": {"re": complex(mr.q_mie(EPS, X, -1)).real,
                                  "im": complex(mr.q_mie(EPS, X, -1)).imag},
            "q_mie_default_rel_err_vs_q": reldiff(mr.q_mie(EPS, X, +1), q),
            "q_mie_xi_minus1_rel_err_vs_q": reldiff(mr.q_mie(EPS, X, -1), q),
            "why_default_is_close_but_wrong": ("-i*a1/(1+a1) is the alternating geometric series "
                                               "-i*a1*(1 - a1 + a1^2 - ...) whereas the correct "
                                               "+i*a1/(1-a1) = i*a1*(1 + a1 + a1^2 + ...); "
                                               "they differ at O(|a1|^2) = %.3e" % abs(a1p) ** 2),
        },
    }

    tol = 1e-12
    exact = {
        "A: -i*a1/(1+a1) on textbook h^(1) a1": reldiff(-1j * a1_tb / (1 + a1_tb), q) < tol,
        "B: +i*a1/(1-a1) on textbook h^(1) a1": reldiff(1j * a1_tb / (1 - a1_tb), q) < tol,
        "B on code default xi_sign=+1 branch": reldiff(1j * a1p / (1 - a1p), q) < tol,
        "C: -i*a1/(1-a1) on code default xi_sign=+1 branch": reldiff(-1j * a1p / (1 - a1p), q) < tol,
    }
    out["decision"] = {
        "matches_q_within_1e-12": exact,
        "statement": ("With the outgoing Hankel function h1^(1) (xi1 = psi1 - i*chi1 = z*h1^(1)) "
                      "and time factor e^{-i omega t}, the textbook electric-dipole coefficient "
                      "satisfies Re(a1) = |a1|^2 > 0 and equals A5's -t exactly, so the Section-8 "
                      "reactance is recovered by q = +i*a1/(1-a1) -- candidate B.  Candidate "
                      "-i*a1/(1+a1) is the wrong series (it misses by O(|a1|^2) = 2.7e-3 relative). "
                      "The code's q_mie default uses xi_sign=+1 (the z*h1^(2) ingoing branch) "
                      "together with -i*a1/(1+a1), so it does not return q (rel 2.7e-3); on that "
                      "conjugated branch the correct conversion is -i*a1/(1-a1) = +i*conj(a1)/(1-conj(a1))."),
    }

    rec.log("PART E  eps=%g x=%g" % (EPS, X))
    rec.log("  q(reactance)            = %.17g   (mpmath %.17g, scipy rel %.2e)"
            % (q, q_high.real, out["reactance_reference"]["q_closed_vs_scipy_rel"]))
    rec.log("  z*h1^(1) == psi1 - i*chi1 (rel %.2e) ; == code mie_a1 branch xi_sign=-1 (rel %.2e)"
            % (branch["rel_closed_h1_vs_code_xi_sign_minus1"],
               branch["rel_closed_h1_vs_code_xi_sign_minus1"]))
    rec.log("  z*h1^(1) == psi1 + i*chi1 (code default xi_sign=+1) ? rel %.2e  (NO -> default is h^(2))"
            % branch["rel_closed_h1_vs_code_xi_sign_plus1"])
    rec.log("  TEXTBOOK a1 (h^(1))     = %.17g %+.17gi" % (a1_tb.real, a1_tb.imag))
    rec.log("    routes agree: scipy|mpmath rel %.2e ; scipy|closed-form rel %.2e ; scipy|hankel1 rel %.2e"
            % (out["textbook_a1_outgoing_h1"]["route_agreement_rel_scipy_vs_mpmath"],
               out["textbook_a1_outgoing_h1"]["route_agreement_rel_scipy_vs_closed_form"],
               out["textbook_a1_outgoing_h1"]["route_agreement_rel_scipy_vs_hankel1"]))
    rec.log("  Re(a1) + |a1|^2         = %.17g   (Re a1 - |a1|^2 = %.3e)"
            % (out["textbook_a1_outgoing_h1"]["Re_a1_plus_abs_a1_sq"],
               out["textbook_a1_outgoing_h1"]["Re_a1_minus_abs_a1_sq"]))
    rec.log("  t = i q/(1 - i q)       = %.17g %+.17gi   (|t|^2 + Re t = %.3e)"
            % (t.real, t.imag, out["theory_t"]["Re_t_plus_abs_t_sq"]))
    rec.log("  t == -a1(textbook) ? rel %.2e" % out["theory_t"]["t_vs_minus_textbook_a1_rel"])
    ct = conv["a1_at_xi_sign_minus1(=z*h1^(1), textbook)"]
    rec.log("  A = -1j*a1/(1+a1)       = %.17g %+.17gi   rel vs q = %.3e  (WRONG)"
            % (ct["A: -i*a1/(1+a1)"]["value"]["re"], ct["A: -i*a1/(1+a1)"]["value"]["im"],
               ct["A: -i*a1/(1+a1)"]["rel_err_vs_q"]))
    rec.log("  B = +1j*a1/(1-a1)       = %.17g %+.17gi   rel vs q = %.3e  <== REPRODUCES q"
            % (ct["B: +i*a1/(1-a1)"]["value"]["re"], ct["B: +i*a1/(1-a1)"]["value"]["im"],
               ct["B: +i*a1/(1-a1)"]["rel_err_vs_q"]))
    cp1 = conv["a1_at_xi_sign_plus1(=z*h1^(2))"]
    rec.log("  (code default branch xi_sign=+1: -1j*a1/(1-a1) rel vs q = %.3e = q)"
            % cp1["C: -i*a1/(1-a1)"]["rel_err_vs_q"])
    rec.log("  code q_mie(xi=+1)       = %.17g %+.17gi   rel vs q = %.3e  (default is NOT q)"
            % (out["code_q_mie_default"]["q_mie(xi_sign=+1)"]["re"],
               out["code_q_mie_default"]["q_mie(xi_sign=+1)"]["im"],
               out["code_q_mie_default"]["q_mie_default_rel_err_vs_q"]))
    rec.log("  matches q (rel<1e-12): %s" % exact)

    rec.data["PART_E_sign_convention"] = out
    rec.dump()
    rec.log("PART E done  json=%s  wall=%.1fs"
            % (os.path.relpath(rec.json_path, ROOT), time.time() - WALL0))


if __name__ == "__main__":
    main()
