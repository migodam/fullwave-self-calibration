"""Interval variants of dq/deps and the analytic second derivative q'' (A5 section 8).

Everything here is written against an injected `iv` module object (mpmath.iv) so the
caller controls iv.dps.  All expressions are mathematically identical; they differ only
in the expression DAG, hence in how much the interval enclosure is widened by
dependency.  See results/DERIVATIONS.md for the derivations.

Naming of the variants follows the task description:
    V1  A5 section 8 chain rule, Djp(z) = -j1'(z)/z - j1(z) + j1(z)/z**2
    V2  same chain rule, Djp(z) from the directly differentiated closed form
    V3  log-derivative rearrangement: q * (N'/N - Den'/Den)
    V4  fully expanded:  W * ( Dj(z)*P - (x/2)*j1(z)*Djp(z) ) / Den**2,
        with W = j1(x)*Dy(x) - y1(x)*Dj(x) a *constant* interval and
        P = j1(z)/(2 sqrt(eps)) + (x/2) j1'(z).  This form removes the N/Den
        interval products entirely from the numerator.
"""

from __future__ import annotations


def _as_iv(x, iv):
    """Coerce x to an mpmath.iv interval.

    IMPORTANT: x must never stay a Python float inside the interval expressions.
    A float would make sub-expressions such as ``x**2 / 2`` evaluate in IEEE double
    arithmetic *outside* the interval machinery, which silently injects a ~1e-16
    relative perturbation into the computed formula (measured: 2.2e-15 relative on
    q'' at (1.5, 0.2)) and breaks containment of the enclosure.
    """
    if hasattr(x, "_mpi_"):
        return x
    return iv.mpf(x)


def prims(m, x, iv):
    """Closed-form primitives at z = m*x and at the fixed argument x."""
    x = _as_iv(x, iv)
    z = m * x
    sz, cz = iv.sin(z), iv.cos(z)
    sx, cx = iv.sin(x), iv.cos(x)
    one = iv.mpf(1)

    j1z = sz / z**2 - cz / z
    Djz = sz / z + cz / z**2 - sz / z**3
    j1pz = sz / z + 2 * cz / z**2 - 2 * sz / z**3

    # two algebraically identical expressions for d/dz D_j(z)
    Djp_chain = -j1pz / z - j1z + j1z / z**2
    Djp_closed = cz / z - 2 * sz / z**2 - 3 * cz / z**3 + 3 * sz / z**4
    # d2/dz2 D_j(z), obtained by direct differentiation of Djp_closed
    Djppz = -sz / z - 3 * cz / z**2 + 7 * sz / z**3 + 12 * cz / z**4 - 12 * sz / z**5
    # j1''(z) from the spherical Bessel ODE (n = 1)
    j1ppz = -(2 / z) * j1pz - (one - 2 / z**2) * j1z

    j1x = sx / x**2 - cx / x
    y1x = -cx / x**2 - sx / x
    Djx = sx / x + cx / x**2 - sx / x**3
    Dyx = -cx / x + sx / x**2 + cx / x**3

    return dict(z=z, j1z=j1z, Djz=Djz, j1pz=j1pz, Djp_chain=Djp_chain,
                Djp_closed=Djp_closed, Djppz=Djppz, j1ppz=j1ppz,
                j1x=j1x, y1x=y1x, Djx=Djx, Dyx=Dyx)


def evaluate(eps_box, x, iv, variant="V1"):
    """Interval enclosure of q, Den and q' on eps_box for the requested variant."""
    x = _as_iv(x, iv)
    m = iv.sqrt(eps_box)
    P = prims(m, x, iv)
    Djpz = P["Djp_closed"] if variant == "V2" else P["Djp_chain"]

    N = m * P["j1z"] * P["Djx"] - P["j1x"] * P["Djz"]
    Den = m * P["j1z"] * P["Dyx"] - P["y1x"] * P["Djz"]
    q = N / Den

    Ap = P["j1z"] / (2 * m) + (x / 2) * P["j1pz"]
    dD = (x / (2 * m)) * Djpz
    Np = Ap * P["Djx"] - P["j1x"] * dD
    Denp = Ap * P["Dyx"] - P["y1x"] * dD

    if variant == "V3":
        qp = q * (Np / N - Denp / Den)
    elif variant == "V4":
        W = P["j1x"] * P["Dyx"] - P["y1x"] * P["Djx"]
        qp = W * (P["Djz"] * Ap - (x / 2) * P["j1z"] * Djpz) / Den**2
    else:
        qp = (Np * Den - N * Denp) / Den**2

    return dict(q=q, Den=Den, qprime=qp, N=N, Denp=Denp, Np=Np, Ap=Ap, P=P)


def qpp(eps_box, x, iv, Djpz=None):
    """Interval enclosure of d2q/deps2 by applying the section 8 chain rule twice.

        s = sqrt(eps),  z = s*x,  u = 1/(2s),  z' = x*u,  z'' = -x/(4 s^3)
        A = s*j1(z),    A' = u*j1(z) + (x/2) j1'(z)
        A'' = -j1(z)/(4 s^3) + x*u^2 j1'(z) + (x^2/2) u j1''(z)
        N''' terms use  d2/de2[D_j(z)] = D_j''(z) z'^2 + D_j'(z) z''.
        q'  = F/Den^2 with F = N' Den - N Den'
        q'' = [ (N'' Den - N Den'') Den - 2 F Den' ] / Den^3
    """
    x = _as_iv(x, iv)
    m = iv.sqrt(eps_box)
    P = prims(m, x, iv)
    D1 = P["Djp_chain"] if Djpz is None else Djpz

    u = 1 / (2 * m)
    zp = x * u
    zpp = -x / (4 * m**3)

    N = m * P["j1z"] * P["Djx"] - P["j1x"] * P["Djz"]
    Den = m * P["j1z"] * P["Dyx"] - P["y1x"] * P["Djz"]
    Ap = P["j1z"] / (2 * m) + (x / 2) * P["j1pz"]
    App = (-1 / (4 * m**3)) * P["j1z"] + x * u**2 * P["j1pz"] + (x**2 / 2) * u * P["j1ppz"]

    d2D = P["Djppz"] * zp**2 + D1 * zpp
    Np = Ap * P["Djx"] - P["j1x"] * D1 * zp
    Denp = Ap * P["Dyx"] - P["y1x"] * D1 * zp
    Npp = App * P["Djx"] - P["j1x"] * d2D
    Denpp = App * P["Dyx"] - P["y1x"] * d2D

    F = Np * Den - N * Denp
    return ((Npp * Den - N * Denpp) * Den - 2 * F * Denp) / Den**3


def iv_bounds(v):
    """(lower, upper) endpoints of an mpmath.iv value as mp.mpf objects.

    NOTE: for mpmath intervals the attributes ``v.a`` / ``v.b`` are themselves
    intervals, so reading them directly (or rendering them) silently returns the
    *lower* bound of each endpoint.  Always go through this helper instead.
    """
    from mpmath import mp
    from mpmath.libmp import to_str
    a, b = v._mpi_[0], v._mpi_[1]
    with mp.extradps(25):
        return mp.mpf(to_str(a, 45)), mp.mpf(to_str(b, 45))


def abs_interval_upper(i):
    """Upper bound of |i| for an interval that does not straddle zero."""
    lo, hi = iv_bounds(i)
    return max(abs(lo), abs(hi))
