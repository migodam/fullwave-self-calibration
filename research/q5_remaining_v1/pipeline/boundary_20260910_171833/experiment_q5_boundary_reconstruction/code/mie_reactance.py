"""A5 (2026-09-10) Section 8 exact spherical electric-dipole reactance coefficient.

Conventions (taken verbatim from A5 Section 8 and Section 5):
    time factor e^{-i w t},  m = sqrt(eps) with eps real >= 1,  x = ka
    j1(z) = sin z / z**2 - cos z / z
    y1(z) = -cos z / z**2 - sin z / z
    D_j(z) = (j1(z) + z j1'(z)) / z        [= psi1'(z)/z,  psi1 = z*j1]
    D_y(z) = (y1(z) + z y1'(z)) / z        [= -chi1'(z)/z, chi1 = -z*y1]
    N(eps,x)   = m j1(mx) D_j(x) - j1(x) D_j(mx)
    Den(eps,x) = m j1(mx) D_y(x) - y1(x) D_j(mx)
    q(eps,x)   = N / Den
    q' = (N' Den - N Den')/Den**2,  m' = 1/(2m),  d(mx)/d eps = x/(2m)
    D_j'(z) = -j1'(z)/z - j1(z) + j1(z)/z**2   (same form for D_y by the Bessel ODE)
    h  = q/sqrt(1+q**2),  h'(eps) = q'/(1+q**2)**(3/2)

Three independent q routes are provided:
    (a) closed-form trigonometric expressions (this file derives and uses them)
    (b) scipy.special.spherical_jn / spherical_yn plus the analytic D-definitions
    (c) standard Mie electric-dipole coefficient a1 -> q = +i a1/(1-a1)

Route (c) note -- SIGN CONVENTION, now pinned down numerically (see
results/mie_sign_convention_*.json and results/FINDINGS_increment3.md).  Define
h1^(1)(z) = j1(z) + i*y1(z) (OUTGOING for the time factor e^{-i omega t}), psi1 = z*j1,
chi1 = -z*y1 and xi1^(s)(z) = psi1 + s*i*chi1.  Then, exactly,
    xi1^(-1)(z) = z*h1^(1)(z) = -e^{iz}(z+i)/z     (outgoing branch; xi_sign = -1)
    xi1^(+1)(z) = z*h1^(2)(z) = conj(xi1^(-1))     (ingoing branch;  xi_sign = +1)
With the OUTGOING h1^(1) textbook coefficient a1 (Re a1 = |a1|^2 for lossless eps, and
a1 = -t with A5's t = i q/(1-i q)) the reactance is recovered exactly by

        q = +i * a1 / (1 - a1)          <-- correct (candidate B)

verified to 2.1e-15 relative.  The alternative "-i*a1/(1+a1)" is NOT q: it is the
alternating geometric series -i*a1*(1 - a1 + a1^2 - ...), so it differs from
+i*a1/(1-a1) = i*a1*(1 + a1 + a1^2 + ...) at O(|a1|^2), i.e. by ~2.7e-3 relative at
eps=2, x=0.2.  On the conjugated (xi_sign = +1) branch the roles swap:
q = -i*a1/(1-a1) = +i*conj(a1)/(1-conj(a1)).  `mie_a1(..., xi_sign=...)` selects the
branch; note that the DEFAULT xi_sign = +1 is the h^(2) branch, not the outgoing one.

No installation is required: mpmath is importable only through the uv cache archive on
PYTHONPATH (pure-Python read of an existing cached wheel, not a pip/system install).

FINAL CONVENTION (Part G of increment 3, verified in logs/mie_fix.log): the textbook
outgoing Riccati-Hankel branch is xi1 = psi1 - i*chi1 = z*h1^(1)(z) (xi_sign = -1),
which is the branch of h1^(1) tied to the time factor e^{-i omega t}; on it Im(a1) < 0
(Re a1 = +|a1|^2 > 0 for lossless eps) and the exact conversion to the A5 Section 8
reactance is q = +i*a1/(1-a1), equivalently q = -i*t/(1+t) with t = -a1 = i*q/(1-i*q).
The code's `mie_a1(..., xi_sign=+1)` default is the conjugated ingoing h1^(2) branch,
where Im(a1) > 0 (Re a1 = -|a1|^2) and the exact conversion is q = -i*a1/(1-a1);
`q_mie` now selects the correct arrangement automatically for whichever branch is
requested and raises rather than returning the historical wrong-arrangement value.
"""

from __future__ import annotations

import numpy as np
from scipy import special

__all__ = [
    "j1_closed", "y1_closed", "j1p_closed", "y1p_closed", "Dj_closed", "Dy_closed",
    "Dj_def", "Dy_def", "Djp_closed", "Dyp_closed",
    "q_closed", "q_scipy", "mie_a1", "q_mie", "q_mie_both_signs",
    "N_Den_closed", "Nprime_Denprime_closed", "qprime_closed", "hprime_closed",
    "xi1", "chi1", "psi1", "psi1p", "chi1p",
    "mp_dps", "mp_q", "mp_qprime", "mp_Den", "mp_quantities", "mp_derivatives",
    "q_from_a1_conversion", "q_route_c_matrix", "A1_CONVERSIONS",
    "iv_dps", "iv_quantities", "decimal_endpoint_intervals", "box_from_intervals",
]


# --------------------------------------------------------------------------------------
# (a) closed-form trigonometric route
# --------------------------------------------------------------------------------------

def j1_closed(z):
    z = np.asarray(z, dtype=float)
    return np.sin(z) / z**2 - np.cos(z) / z


def y1_closed(z):
    z = np.asarray(z, dtype=float)
    return -np.cos(z) / z**2 - np.sin(z) / z


def j1p_closed(z):
    """d/dz j1(z).  Derived from j1 = sin z/z**2 - cos z/z."""
    z = np.asarray(z, dtype=float)
    return np.sin(z) / z + 2.0 * np.cos(z) / z**2 - 2.0 * np.sin(z) / z**3


def y1p_closed(z):
    """d/dz y1(z).  Derived from y1 = -cos z/z**2 - sin z/z."""
    z = np.asarray(z, dtype=float)
    return -np.cos(z) / z + 2.0 * np.sin(z) / z**2 + 2.0 * np.cos(z) / z**3


def Dj_closed(z):
    """D_j(z) = (j1 + z j1')/z = psi1'/z."""
    z = np.asarray(z, dtype=float)
    return np.sin(z) / z + np.cos(z) / z**2 - np.sin(z) / z**3


def Dy_closed(z):
    """D_y(z) = (y1 + z y1')/z = -chi1'/z."""
    z = np.asarray(z, dtype=float)
    return -np.cos(z) / z + np.sin(z) / z**2 + np.cos(z) / z**3


def Dj_def(z):
    """D_j from its definition, using scipy's j1 and j1'."""
    z = np.asarray(z, dtype=float)
    return (special.spherical_jn(1, z) + z * special.spherical_jn(1, z, derivative=True)) / z


def Dy_def(z):
    """D_y from its definition, using scipy's y1 and y1'."""
    z = np.asarray(z, dtype=float)
    return (special.spherical_yn(1, z) + z * special.spherical_yn(1, z, derivative=True)) / z


def Djp_closed(z):
    """D_j'(z) = -j1'(z)/z - j1(z) + j1(z)/z**2  (A5 Section 8)."""
    z = np.asarray(z, dtype=float)
    return -j1p_closed(z) / z - j1_closed(z) + j1_closed(z) / z**2


def Dyp_closed(z):
    """D_y'(z) = -y1'(z)/z - y1(z) + y1(z)/z**2  (same form: y1 obeys the same ODE)."""
    z = np.asarray(z, dtype=float)
    return -y1p_closed(z) / z - y1_closed(z) + y1_closed(z) / z**2


def N_Den_closed(eps, x):
    eps = np.asarray(eps, dtype=float)
    m = np.sqrt(eps)
    mx = m * x
    N = m * j1_closed(mx) * Dj_closed(x) - j1_closed(x) * Dj_closed(mx)
    Den = m * j1_closed(mx) * Dy_closed(x) - y1_closed(x) * Dj_closed(mx)
    return N, Den


def q_closed(eps, x):
    N, Den = N_Den_closed(eps, x)
    return N / Den


def Nprime_Denprime_closed(eps, x):
    """Analytic d/d eps of N and Den (A5 Section 8 chain rule)."""
    eps = np.asarray(eps, dtype=float)
    m = np.sqrt(eps)
    mx = m * x
    A = m * j1_closed(mx)
    Ap = j1_closed(mx) / (2.0 * m) + (x / 2.0) * j1p_closed(mx)
    dDjmx = (x / (2.0 * m)) * Djp_closed(mx)
    Np = Ap * Dj_closed(x) - j1_closed(x) * dDjmx
    Denp = Ap * Dy_closed(x) - y1_closed(x) * dDjmx
    return Np, Denp


def qprime_closed(eps, x):
    N, Den = N_Den_closed(eps, x)
    Np, Denp = Nprime_Denprime_closed(eps, x)
    return (Np * Den - N * Denp) / Den**2


def hprime_closed(eps, x):
    q = q_closed(eps, x)
    qp = qprime_closed(eps, x)
    return qp / (1.0 + q**2) ** 1.5


# --------------------------------------------------------------------------------------
# (b) scipy route: same D-definitions, scipy special functions + its derivative kwarg
# --------------------------------------------------------------------------------------

def q_scipy(eps, x):
    eps = np.asarray(eps, dtype=float)
    m = np.sqrt(eps)
    mx = m * x
    j1x = special.spherical_jn(1, x)
    j1mx = special.spherical_jn(1, mx)
    y1x = special.spherical_yn(1, x)
    Djx = Dj_def(x)
    Djmx = Dj_def(mx)
    Dyx = Dy_def(x)
    N = m * j1mx * Djx - j1x * Djmx
    Den = m * j1mx * Dyx - y1x * Djmx
    return N / Den


# --------------------------------------------------------------------------------------
# (c) Mie route: standard electric-dipole a1, then a branch-aware a1 -> q conversion.
#
#     B R A N C H   T A B L E   (verified numerically against q_closed, see mie_fix.log)
#     ------------------------------------------------------------------------------
#     xi_sign = -1 : xi1 = psi1 - i*chi1 = z*h1^(1)(z)   OUTGOING for e^{-i omega t}
#                    a1 = textbook Mie coefficient (Re a1 = +|a1|^2 for lossless eps)
#                    a1 = -t with A5 t = i q/(1-i q)
#                    q  =  +i*a1/(1-a1)         <-- correct on this branch
#     xi_sign = +1 : xi1 = psi1 + i*chi1 = z*h1^(2)(z)   INGOING for e^{-i omega t}
#                    a1 = conj(textbook a1)     (Re a1 = -|a1|^2)
#                    q  =  -i*a1/(1-a1)         <-- correct on this branch
#
#     The arrangement " -i*a1/(1+a1) " is WRONG on both branches: it is the
#     alternating geometric series, differs from the correct value at O(|a1|^2)
#     (~2.7e-3 relative at eps=2, x=0.2).  `q_mie` raises if it is requested and
#     `q_mie_legacy_prompt_literal` reproduces the historical value on purpose.
# --------------------------------------------------------------------------------------

def psi1(z):
    z = np.asarray(z, dtype=complex)
    return z * special.spherical_jn(1, z)


def chi1(z):
    z = np.asarray(z, dtype=complex)
    return -z * special.spherical_yn(1, z)


def psi1p(z):
    z = np.asarray(z, dtype=complex)
    return special.spherical_jn(1, z) + z * special.spherical_jn(1, z, derivative=True)


def chi1p(z):
    z = np.asarray(z, dtype=complex)
    return -(special.spherical_yn(1, z) + z * special.spherical_yn(1, z, derivative=True))


def xi1(z, xi_sign=+1):
    """xi1 = psi1 + xi_sign*i*chi1, chi1 = -z*y1(z).

    xi_sign = -1 is z*h1^(1)(z), the OUTGOING combination for time factor
    e^{-i omega t} and the branch on which q = +i*a1/(1-a1) reproduces A5 Section 8 q.
    xi_sign = +1 (the default) is z*h1^(2)(z) = conj(xi1(-1)) on the real axis.
    """
    return psi1(z) + xi_sign * 1j * chi1(z)


def mie_a1(eps, x, xi_sign=+1):
    """Standard Mie electric-dipole coefficient a1 (Riccati-Bessel form).

    a1 = [m psi1(mx) psi1'(x) - psi1(x) psi1'(mx)]
         / [m psi1(mx) xi1'(x) - xi1(x) psi1'(mx)]

    xi_sign = -1 gives the textbook coefficient built with the OUTGOING h1^(1)
    (outgoing for e^{-i omega t}); xi_sign = +1 (the default) is its conjugate,
    built with h1^(2).
    """
    eps = np.asarray(eps, dtype=float)
    m = np.sqrt(eps)
    mx = m * x
    xi = xi1(x, xi_sign)
    dxi = psi1p(x) + xi_sign * 1j * chi1p(x)
    num = m * psi1(mx) * psi1p(x) - psi1(x) * psi1p(mx)
    den = m * psi1(mx) * dxi - xi * psi1p(mx)
    return num / den


_Q_MIE_CORRECT_FORM = {+1: "minus_i_over_1minus", -1: "i_over_1minus"}
_Q_MIE_WRONG_FORMS = ("minus_i_over_1plus",)


def q_mie(eps, x, xi_sign=+1, form="auto", check=True, rtol=1e-9):
    """Route (c): branch-aware a1 -> q conversion, verified against the reactance q.

    Branch table (both rows verified against `q_closed` to <= 1e-14 in
    logs/mie_fix.log):
        xi_sign = +1  (z*h1^(2), INGOING for e^{-i omega t})  ->  q = -i*a1/(1-a1)
        xi_sign = -1  (z*h1^(1), OUTGOING for e^{-i omega t}) ->  q = +i*a1/(1-a1)

    Parameters
    ----------
    xi_sign : int
        +1 (default) selects z*h1^(2); -1 selects z*h1^(1).
    form : str
        "auto"                  -> the correct arrangement for `xi_sign`.
        "i_over_1minus"         -> +i*a1/(1-a1)
        "minus_i_over_1minus"   -> -i*a1/(1-a1)
        "minus_i_over_1plus"    -> NOT a legal choice: this is the historical /
                                   prompt-literal arrangement, the wrong-sign
                                   alternating series, wrong by O(|a1|^2).
                                   Requesting it raises ValueError.  The value is
                                   still available, explicitly flagged, via
                                   `q_mie_legacy_prompt_literal`.
    check : bool
        If True (default) the returned value is compared with `q_closed` at the
        same (eps, x) and ValueError is raised when the relative mismatch
        exceeds `rtol`.  This makes a silent wrong-arrangement return impossible.
    rtol : float
        Relative tolerance for the internal cross-check (default 1e-9).  The
        actual agreement with q_closed is <= 5.4e-12 over the mie_checks grid
        (worst point eps=2, x=0.02, where route (c) loses digits), so 1e-9 admits
        every legitimate point while still rejecting the wrong arrangement by a
        factor of >1e6.

    Returns
    -------
    complex q (float64), equal to the A5 Section 8 reactance q up to float64
    round-off.
    """
    if form == "auto":
        form = _Q_MIE_CORRECT_FORM[xi_sign]
    if form in _Q_MIE_WRONG_FORMS:
        raise ValueError(
            "q_mie(form=%r) is the wrong-arrangement (alternating-series) form: it differs "
            "from the A5 Section 8 reactance q at O(|a1|^2).  Use form='auto' (default), or "
            "q_mie_legacy_prompt_literal() if the historical value is genuinely wanted." % (form,)
        )
    if form not in ("i_over_1minus", "minus_i_over_1minus"):
        raise ValueError("unknown form %r" % (form,))
    a1 = mie_a1(eps, x, xi_sign=xi_sign)
    q = (1j * a1 / (1.0 - a1)) if form == "i_over_1minus" else (-1j * a1 / (1.0 - a1))
    if check:
        ref = complex(q_closed(eps, x))
        got = complex(q)
        rel = abs(got - ref) / abs(ref)
        if not (rel <= rtol):
            raise ValueError(
                "q_mie internal check FAILED: xi_sign=%+d form=%s gives %r, reactance q is %r, "
                "rel=%.3e > rtol=%.3e" % (xi_sign, form, got, ref, rel, rtol)
            )
    return q


def q_mie_legacy_prompt_literal(eps, x, xi_sign=+1):
    """Historical (wrong) arrangement -i*a1/(1+a1), kept ONLY for record reproduction.

    This is the arrangement that the pre-fix `q_mie` default returned.  It does not
    solve for the A5 Section 8 reactance q on either branch; it is the alternating
    geometric series, wrong by O(|a1|^2) (~2.7e-3 relative at eps=2, x=0.2).  Any
    value obtained here must be reported as such.
    """
    a1 = mie_a1(eps, x, xi_sign=xi_sign)
    return -1j * a1 / (1.0 + a1)


def q_mie_both_signs(eps, x):
    """Correctly converted q on both branches (h^(2) first, then h^(1))."""
    return q_mie(eps, x, +1), q_mie(eps, x, -1)


# --------------------------------------------------------------------------------------
# a1 -> q conversion variants (see results/DERIVATIONS.md for the sign analysis)
# --------------------------------------------------------------------------------------

A1_CONVERSIONS = ("minus_i_over_1plus", "i_over_1minus")


def q_from_a1_conversion(eps, x, xi_sign=+1, form="i_over_1minus"):
    """Two candidate a1 -> q conversions, resolved numerically in Part E.

    Textbook coefficient: a1 built with the OUTGOING h1^(1), i.e. xi_sign = -1
    (a1 = -t with A5's t = iq/(1-iq)).  For that coefficient
        form='i_over_1minus'    ->  +i*a1/(1-a1)  =  q      (CORRECT, rel 2.1e-15)
        form='minus_i_over_1plus' -> -i*a1/(1+a1) != q      (rel 2.0 at eps=2, x=0.2)
    On the conjugated branch (xi_sign = +1, the default here) the correct pairing
    flips: form='i_over_1minus' returns -q and form='minus_i_over_1plus' returns
    the O(|a1|^2) wrong-sign estimate.  Use xi_sign=-1 with form='i_over_1minus'.
    """
    a1 = mie_a1(eps, x, xi_sign=xi_sign)
    if form == "minus_i_over_1plus":
        return -1j * a1 / (1.0 + a1)
    if form == "i_over_1minus":
        return 1j * a1 / (1.0 - a1)
    raise ValueError(form)


def q_route_c_matrix(eps, x):
    """All four (xi sign) x (conversion) combinations of route (c)."""
    out = {}
    for s in (+1, -1):
        for f in A1_CONVERSIONS:
            out[(s, f)] = q_from_a1_conversion(eps, x, xi_sign=s, form=f)
    return out


# --------------------------------------------------------------------------------------
# generic backend core (numpy / mpmath.mp / mpmath.iv)
# --------------------------------------------------------------------------------------

def _core(m, x, sin, cos, msqrt_sq):
    """All quantities for a given m (=sqrt(eps)) and x.

    `sin`, `cos` are the backend trig functions; `msqrt_sq(y)` = y**(3/2).
    Works unchanged for floats/arrays, mpmath.mp scalars and mpmath.iv intervals.
    """
    z = m * x
    sz, cz = sin(z), cos(z)
    sx, cx = sin(x), cos(x)

    j1z = sz / z**2 - cz / z
    y1z = -cz / z**2 - sz / z
    Djz = sz / z + cz / z**2 - sz / z**3
    Dyz = -cz / z + sz / z**2 + cz / z**3
    j1pz = sz / z + 2 * cz / z**2 - 2 * sz / z**3
    y1pz = -cz / z + 2 * sz / z**2 + 2 * cz / z**3
    Djpz = -j1pz / z - j1z + j1z / z**2
    Dypz = -y1pz / z - y1z + y1z / z**2

    j1x = sx / x**2 - cx / x
    y1x = -cx / x**2 - sx / x
    Djx = sx / x + cx / x**2 - sx / x**3
    Dyx = -cx / x + sx / x**2 + cx / x**3

    N = m * j1z * Djx - j1x * Djz
    Den = m * j1z * Dyx - y1x * Djz
    q = N / Den

    A = m * j1z
    Ap = j1z / (2 * m) + (x / 2) * j1pz
    dDjmx = (x / (2 * m)) * Djpz
    Np = Ap * Djx - j1x * dDjmx
    Denp = Ap * Dyx - y1x * dDjmx
    qp = (Np * Den - N * Denp) / Den**2
    hp = qp / msqrt_sq(1 + q**2)

    return {
        "j1": j1z, "y1": y1z, "Dj": Djz, "Dy": Dyz,
        "j1p": j1pz, "y1p": y1pz, "Djp": Djpz, "Dyp": Dypz,
        "N": N, "Den": Den, "q": q, "Np": Np, "Denp": Denp,
        "qprime": qp, "hprime": hp, "mx": z,
    }


# --------------------------------------------------------------------------------------
# mpmath high-precision scalar versions (reference values)
# --------------------------------------------------------------------------------------

mp_dps = 50


def _mp():
    from mpmath import mp
    mp.dps = mp_dps
    return mp


def mp_quantities(eps, x, dps=None):
    mp = _mp()
    if dps is not None:
        mp.dps = dps
    eps = mp.mpf(eps)
    x = mp.mpf(x)
    m = mp.sqrt(eps)
    return _core(m, x, mp.sin, mp.cos, lambda y: y**mp.mpf("1.5"))


def mp_q(eps, x, dps=None):
    return mp_quantities(eps, x, dps)["q"]


def mp_qprime(eps, x, dps=None):
    return mp_quantities(eps, x, dps)["qprime"]


def mp_Den(eps, x, dps=None):
    return mp_quantities(eps, x, dps)["Den"]


def mp_derivatives(z, dps=None):
    """D_j, D_y, D_j', D_y' and j1, y1, j1', y1' at a scalar z, high precision."""
    mp = _mp()
    if dps is not None:
        mp.dps = dps
    z = mp.mpf(z)
    one = mp.mpf(1)
    sz, cz = mp.sin(z), mp.cos(z)
    j1 = sz / z**2 - cz / z
    y1 = -cz / z**2 - sz / z
    j1p = sz / z + 2 * cz / z**2 - 2 * sz / z**3
    y1p = -cz / z + 2 * sz / z**2 + 2 * cz / z**3
    d = {
        "j1": j1, "y1": y1, "j1p": j1p, "y1p": y1p,
        "Dj": sz / z + cz / z**2 - sz / z**3,
        "Dy": -cz / z + sz / z**2 + cz / z**3,
        "Djp": -j1p / z - j1 + j1 / z**2,
        "Dyp": -y1p / z - y1 + y1 / z**2,
        "one": one,
    }
    return d


# --------------------------------------------------------------------------------------
# mpmath.iv interval versions
# --------------------------------------------------------------------------------------

iv_dps = 35


def _iv():
    from mpmath import iv
    iv.dps = iv_dps
    return iv


def _iv_midpoint_floor(i):
    """Lower endpoint of an mpmath.iv interval."""
    return i.a


def _iv_upper(i):
    return i.b


def iv_quantities(eps_box, x, dps=None):
    """Interval enclosure of q, q', Den, h' for eps in the interval `eps_box`.

    eps_box is an mpmath.iv interval (or a 2-sequence of exact fractions).
    x is converted to an enclosing point interval.
    """
    iv = _iv()
    if dps is not None:
        iv.dps = dps
    if not hasattr(eps_box, "a"):
        lo, hi = eps_box
        eps_box = iv.mpf([lo, hi])
    x = iv.mpf(x)
    m = iv.sqrt(eps_box)
    out = _core(m, x, iv.sin, iv.cos, lambda y: y * iv.sqrt(y))
    return out


def iv_min_abs(interval):
    """Lower endpoint of |interval| (0 if the interval straddles 0)."""
    lo, hi = interval.a, interval.b
    if lo > 0:
        return lo
    if hi < 0:
        return -hi
    return 0


def decimal_endpoint_intervals(lo_num, lo_den, hi_num, hi_den, n, dps=None):
    """n+1 mpmath.iv point intervals that provably enclose the exact decim　als
    lo + k*(hi-lo)/n built from exact small integer literals.

        left  = iv.mpf(lo_num)/iv.mpf(lo_den)
        right = iv.mpf(hi_num)/iv.mpf(hi_den)
        step  = (right - left)/n
        endpoint_k = left + k*step
    """
    iv = _iv()
    if dps is not None:
        iv.dps = dps
    left = iv.mpf(lo_num) / iv.mpf(lo_den)
    right = iv.mpf(hi_num) / iv.mpf(hi_den)
    step = (right - left) / iv.mpf(n)
    return [left + iv.mpf(k) * step for k in range(n + 1)]


def box_from_intervals(e_lo, e_hi):
    """Smallest mpmath.iv interval containing both endpoint intervals."""
    iv = _iv()
    return iv.mpf([e_lo.a, e_hi.b])
