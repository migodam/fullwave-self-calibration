#!/usr/bin/env python3
"""PART D -- A5 section 9 (reference GLS == joint gain elimination) checks.

Theory read directly from inputs/derive/A5_sec_09.md:
    L(g) = ||y - g f||^2/sigma^2 + |z - g|^2/sigma_r^2
    min_g L(g) = (y - z f)^* (sigma^2 I + sigma_r^2 f f^*)^(-1) (y - z f)
    g_* = (f^* y/sigma^2 + z/sigma_r^2) / (||f||^2/sigma^2 + 1/sigma_r^2)
    no log-determinant term; constrained |g| in [0.75,1.25] needs the constrained min.
    v(theta) = (f(theta)/sigma, 1/sigma_r), b = (y/sigma, z/sigma_r),
    g = v^*b/(v^*v), r = b - v g,
    dr = -Q_v (dv) g - v ((dv)^* r)/(v^* v),  Q_v = I - v v^*/(v^* v).

Q_v is the orthogonal projector onto the orthogonal complement of v: the derivation
    dg = [(dv)^*r - v^*(dv) g]/(v^*v),  dr = -(dv) g - v dg
gives -(dv)g + v v^*(dv) g/(v^*v) = -(I - v v^*/(v^*v))(dv) g, i.e. Q_v = I - vv^*/(v^*v).
The alternative "2I - vv^*/(v^*v)" is NOT a projector and is not used.

No installs, single thread, fixed seed, incremental output.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2, kstest

SEED = 20260910
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WALL0 = time.time()
ANNULUS = (0.75, 1.25)


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


class Rec:
    def __init__(self):
        self.data = {}
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        base = os.path.join(ROOT, "results", "ref_gls_%s" % stamp)
        n = 0
        while os.path.exists(base + (".json" if n == 0 else "_%d.json" % n)):
            n += 1
        self.json_path = base + (".json" if n == 0 else "_%d.json" % n)
        self.log_path = os.path.join(ROOT, "logs", "ref_gls.log")
        os.makedirs(os.path.dirname(self.log_path), exist_ok=True)
        self.logf = open(self.log_path, "a")

    def log(self, msg):
        line = "[%s +%6.1fs] %s" % (time.strftime("%H:%M:%S"), time.time() - WALL0, msg)
        print(line, flush=True)
        self.logf.write(line + "\n")
        self.logf.flush()

    def put(self, key, val):
        self.data[key] = val

    def dump(self):
        self.data.setdefault("meta", {})
        self.data["meta"].update(
            {
                "seed": SEED,
                "script": os.path.relpath(__file__, ROOT),
                "interpreter": sys.executable,
                "numpy_version": np.__version__,
                "wall_seconds": round(time.time() - WALL0, 3),
                "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
        tmp = self.json_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(jsonable(self.data), fh, indent=1)
        os.replace(tmp, self.json_path)


rec = None


# ------------------------------------------------------------------ core identities
def L_of(g, y, f, z, sigma, sigma_r):
    return float(np.linalg.norm(y - g * f) ** 2 / sigma ** 2
                 + abs(z - g) ** 2 / sigma_r ** 2)


def g_star_of(y, f, z, sigma, sigma_r):
    return (np.vdot(f, y) / sigma ** 2 + z / sigma_r ** 2) / (
        np.vdot(f, f).real / sigma ** 2 + 1.0 / sigma_r ** 2)


def gls_rhs(y, f, z, sigma, sigma_r, mode="woodbury"):
    u = y - z * f
    nf2 = float(np.vdot(f, f).real)
    if mode == "woodbury":
        return float((np.vdot(u, u).real - sigma_r ** 2 * abs(np.vdot(f, u)) ** 2
                      / (sigma ** 2 + sigma_r ** 2 * nf2)) / sigma ** 2)
    m = f.shape[0]
    A = sigma ** 2 * np.eye(m) + sigma_r ** 2 * np.outer(f, f.conj())
    return float(np.vdot(u, np.linalg.solve(A, u)).real)


def augmented_min(y, f, z, sigma, sigma_r):
    v = np.concatenate([f / sigma, [1.0 / sigma_r]])
    b = np.concatenate([y / sigma, [z / sigma_r]])
    gs = np.vdot(v, b) / np.vdot(v, v).real
    return float(np.linalg.norm(b - gs * v) ** 2), gs


# ------------------------------------------------------------------ nonlinear field
def make_field(rng, K=5, p=1, m=None, om_scale=1.0):
    m = m or K
    C = rng.standard_normal((K, m)) + 1j * rng.standard_normal((K, m))
    W = rng.uniform(-om_scale, om_scale, size=(K, p))
    return C, W


def field_F(C, W, thetas):
    """F[i, :] = sum_k C[k, :] exp(i W[k] . thetas[i]);  thetas: (n, p)."""
    E = np.exp(1j * (thetas @ W.T))
    return E @ C


def field_J(C, W, theta):
    """df/dtheta: (m, p) at a single theta (p-vector)."""
    K, m = C.shape
    E = np.exp(1j * (W @ theta))
    return (C * (1j * E)[:, None]).T @ W


def J1_theta(C, W, theta, y, z, sigma, sigma_r):
    f = (np.exp(1j * (W @ theta))[:, None] * C).sum(axis=0)
    return gls_rhs(y, f, z, sigma, sigma_r), f


def J1_grid(C, W, thetas, y, z, sigma, sigma_r):
    F = field_F(C, W, thetas)
    u = y[None, :] - z * F
    nf2 = (F.conj() * F).real.sum(axis=1)
    # Woodbury numerator is f^*(y - z f); the earlier y-only version was a bug
    # (it dropped the -z||f||^2 part and made the grid minimise a different objective).
    fu = (F.conj() * u).sum(axis=1)
    J = (np.abs(u) ** 2).sum(axis=1) / sigma ** 2 - sigma_r ** 2 * np.abs(fu) ** 2 / (
        sigma ** 2 + sigma_r ** 2 * nf2) / sigma ** 2
    return J


def Jplug_grid(C, W, thetas, y, z, sigma):
    F = field_F(C, W, thetas)
    return (np.abs(y[None, :] - z * F) ** 2).sum(axis=1) / sigma ** 2


def profiled_theta(obj, C, W, y, z, sigma, sigma_r, lo, hi, n_grid=1201):
    thetas = np.linspace(lo, hi, n_grid)[:, None]
    grid = (J1_grid(C, W, thetas, y, z, sigma, sigma_r) if obj is J1_theta
            else Jplug_grid(C, W, thetas, y, z, sigma))
    i = int(np.argmin(grid))
    a = thetas[max(i - 1, 0), 0]
    b = thetas[min(i + 1, n_grid - 1), 0]
    res = minimize_scalar(lambda t: (J1_theta(C, W, np.array([t]), y, z, sigma, sigma_r)[0]
                                     if obj is J1_theta
                                     else float(Jplug_grid(C, W, np.array([[t]]), y, z, sigma)[0])),
                          bounds=(a, b), method="bounded", options={"xatol": 1e-13})
    return float(res.x), float(res.fun)


# ============================================================ D1: scalar identity
def part_d1():
    rng = np.random.default_rng(SEED)
    sigma, sigma_r = 0.3, 0.05
    m = 6
    f = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    y = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    out = {}
    for ztag, z in (("real_z", 0.8), ("complex_z", 0.7 - 0.45j)):
        gs = g_star_of(y, f, z, sigma, sigma_r)
        L_gs = L_of(gs, y, f, z, sigma, sigma_r)
        # (a) dense grid + local refine
        gr = max(abs(gs), 1.0)
        xs = np.linspace(gs.real - 4 * gr, gs.real + 4 * gr, 401)
        ys = np.linspace(gs.imag - 4 * gr, gs.imag + 4 * gr, 401)
        gx, gy = np.meshgrid(xs, ys)
        G = gx + 1j * gy
        vals = (np.abs(y[None, None, :] - G[..., None] * f[None, None, :]) ** 2).sum(-1) / sigma ** 2 \
            + np.abs(z - G) ** 2 / sigma_r ** 2
        g0 = G[np.unravel_index(int(np.argmin(vals)), G.shape)]
        ref = minimize(lambda p: L_of(p[0] + 1j * p[1], y, f, z, sigma, sigma_r),
                       np.array([g0.real, g0.imag]), method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-15, "fatol": 1e-20})
        L_dense = L_of(ref.x[0] + 1j * ref.x[1], y, f, z, sigma, sigma_r)
        rhs_w = gls_rhs(y, f, z, sigma, sigma_r, "woodbury")
        rhs_d = gls_rhs(y, f, z, sigma, sigma_r, "direct")
        L_aug, g_aug = augmented_min(y, f, z, sigma, sigma_r)
        refv = max(abs(rhs_w), 1e-300)
        out[ztag] = {
            "z": cplx(z), "sigma": sigma, "sigma_r": sigma_r, "m": m,
            "g_star_formula": cplx(gs),
            "g_star_dense_argmin": cplx(ref.x[0] + 1j * ref.x[1]),
            "g_star_dense_vs_formula_rel": float(abs((ref.x[0] + 1j * ref.x[1]) - gs) / abs(gs)),
            "L_at_g_star_formula": L_gs,
            "L_dense_min": L_dense,
            "RHS_woodbury": rhs_w,
            "RHS_direct_solve": rhs_d,
            "augmented_min_value": L_aug,
            "augmented_g_star": cplx(g_aug),
            "augmented_g_vs_formula_rel": float(abs(g_aug - gs) / abs(gs)),
            "rel_discrepancy_L_vs_RHS": float(abs(L_gs - rhs_w) / refv),
            "rel_discrepancy_L_vs_augmented": float(abs(L_gs - L_aug) / refv),
            "rel_discrepancy_woodbury_vs_direct": float(abs(rhs_w - rhs_d) / refv),
            "rel_discrepancy_dense_vs_formula_value": float(abs(L_dense - L_gs) / refv),
        }
        rec.log("D1 %-10s g*=%s L(g*)=%.12e RHS(w)=%.12e RHS(d)=%.12e aug=%.12e "
                "rel(L-RHS)=%.2e rel(w-d)=%.2e" % (ztag, np.format_float_scientific(abs(gs), 6),
                                                   L_gs, rhs_w, rhs_d, L_aug,
                                                   abs(L_gs - rhs_w) / refv,
                                                   abs(rhs_w - rhs_d) / refv))
    rec.put("D1_scalar_identity", out)
    rec.dump()


# ============================================================ D2: log-det + constraint
def part_d2():
    rng = np.random.default_rng(SEED + 1)
    out = {}
    # ---- log-determinant marginal-likelihood variant vs concentrated GLS
    sigma = 0.5
    trials = []
    for _ in range(5):
        C, W = make_field(rng, K=4, p=1, m=5, om_scale=1.3)
        theta0 = 0.4
        f0 = (np.exp(1j * (W @ np.array([theta0])))[:, None] * C).sum(axis=0)
        g = 1.0 * np.exp(1j * 0.2)
        sigma_r = 0.9   # deliberately strong reference noise so the log-det term bites
        y = g * f0 + (sigma / np.sqrt(2)) * (rng.standard_normal(5) + 1j * rng.standard_normal(5))
        z = g + (sigma_r / np.sqrt(2)) * (rng.standard_normal() + 1j * rng.standard_normal())

        def Jgls(t):
            return J1_theta(C, W, np.array([t]), y, z, sigma, sigma_r)[0]

        def Jld(t):
            f = (np.exp(1j * (W @ np.array([t])))[:, None] * C).sum(axis=0)
            nf2 = float(np.vdot(f, f).real)
            m = f.shape[0]
            return (J1_theta(C, W, np.array([t]), y, z, sigma, sigma_r)[0]
                    + m * np.log(sigma ** 2) + np.log(1 + sigma_r ** 2 * nf2 / sigma ** 2))

        r1 = minimize_scalar(Jgls, bounds=(-2.0, 3.0), method="bounded",
                             options={"xatol": 1e-13})
        r2 = minimize_scalar(Jld, bounds=(-2.0, 3.0), method="bounded",
                             options={"xatol": 1e-13})
        trials.append({"theta_hat_GLS": float(r1.x), "theta_hat_logdet": float(r2.x),
                       "abs_diff": float(abs(r1.x - r2.x)),
                       "rel_diff": float(abs(r1.x - r2.x) / max(abs(r1.x), 1e-12)),
                       "Jgls_at_gls": float(r1.fun), "Jld_at_logdet": float(r2.fun),
                       "Jgls_at_logdet": float(Jgls(r2.x)), "Jld_at_gls": float(Jld(r1.x))})
    out["logdet_vs_gls"] = {
        "sigma": sigma, "sigma_r": 0.9, "n_trials": len(trials),
        "max_abs_theta_diff": max(t["abs_diff"] for t in trials),
        "median_abs_theta_diff": float(np.median([t["abs_diff"] for t in trials])),
        "note": ("log-det objective = concentrated GLS + log det(sigma^2 I + sigma_r^2 f f^*); "
                 "this is a different (marginal) model, not the same profiling target."),
        "trials": trials,
    }
    rec.log("D2 log-det: max |theta_gls - theta_ld| = %.6e, median %.6e"
            % (out["logdet_vs_gls"]["max_abs_theta_diff"],
               out["logdet_vs_gls"]["median_abs_theta_diff"]))

    # ---- constrained |g| in [0.75, 1.25]
    sigma_r = 0.05
    f = rng.standard_normal(6) + 1j * rng.standard_normal(6)
    z = 3.0 + 0.0j
    y = 0.9 * f + (sigma / np.sqrt(2)) * (rng.standard_normal(6) + 1j * rng.standard_normal(6))
    gs = g_star_of(y, f, z, sigma, sigma_r)
    L_unc = gls_rhs(y, f, z, sigma, sigma_r)
    res = minimize(lambda p: L_of(p[0] + 1j * p[1], y, f, z, sigma, sigma_r), [1.0, 0.0],
                   constraints=[{"type": "ineq",
                                 "fun": lambda p: ANNULUS[1] ** 2 - p[0] ** 2 - p[1] ** 2},
                                {"type": "ineq",
                                 "fun": lambda p: p[0] ** 2 + p[1] ** 2 - ANNULUS[0] ** 2}],
                   method="SLSQP", options={"maxiter": 3000, "ftol": 1e-16})
    gc = res.x[0] + 1j * res.x[1]
    out["constrained_gain_annulus"] = {
        "annulus": list(ANNULUS), "z": cplx(z), "sigma": sigma, "sigma_r": sigma_r,
        "g_star_unconstrained": cplx(gs), "abs_g_star": float(abs(gs)),
        "unconstrained_outside_annulus": bool(not (ANNULUS[0] <= abs(gs) <= ANNULUS[1])),
        "L_unconstrained_formula": L_unc,
        "L_at_unconstrained_g_star": L_of(gs, y, f, z, sigma, sigma_r),
        "constrained_g_hat": cplx(gc), "abs_constrained_g": float(abs(gc)),
        "L_constrained": float(L_of(gc, y, f, z, sigma, sigma_r)),
        "constrained_minus_unconstrained": float(L_of(gc, y, f, z, sigma, sigma_r) - L_unc),
        "relative_difference": float(abs(L_of(gc, y, f, z, sigma, sigma_r) - L_unc) / L_unc),
        "slsqp_success": bool(res.success),
    }
    rec.log("D2 constrained: |g*|=%.4f outside annulus=%s ; L_unc=%.6e L_con=%.6e rel=%.4e"
            % (abs(gs), out["constrained_gain_annulus"]["unconstrained_outside_annulus"],
               L_unc, out["constrained_gain_annulus"]["L_constrained"],
               out["constrained_gain_annulus"]["relative_difference"]))
    rec.put("D2_logdet_and_constraint", out)
    rec.dump()


# ============================================================ D3: non-linear profiling
def part_d3(n_draws=5000):
    rng = np.random.default_rng(SEED + 2)
    sigma, sigma_r = 0.4, 0.35
    C, W = make_field(rng, K=5, p=1, m=5, om_scale=1.1)
    theta0 = 0.6
    out = {"sigma": sigma, "sigma_r": sigma_r, "m": int(C.shape[1]), "K": int(C.shape[0]),
           "theta_true": theta0}

    f0 = (np.exp(1j * (W @ np.array([theta0])))[:, None] * C).sum(axis=0)
    g0 = 1.0 * np.exp(1j * 0.9)
    y = g0 * f0 + (sigma / np.sqrt(2)) * (rng.standard_normal(5) + 1j * rng.standard_normal(5))
    z = g0 + (sigma_r / np.sqrt(2)) * (rng.standard_normal() + 1j * rng.standard_normal())

    t_gls, v_gls = profiled_theta(J1_theta, C, W, y, z, sigma, sigma_r, -2.0, 3.0)
    t_plug, v_plug = profiled_theta(lambda *a: None, C, W, y, z, sigma, sigma_r, -2.0, 3.0)
    # joint VarPro over (theta, Re g, Im g) with the reference
    Fg = field_F(C, W, np.array([[t_gls]]))

    def J2(p):
        f = (np.exp(1j * (W @ p[:1]))[:, None] * C).sum(axis=0)
        g = p[1] + 1j * p[2]
        return L_of(g, y, f, z, sigma, sigma_r)

    def J3(p):
        f = (np.exp(1j * (W @ p[:1]))[:, None] * C).sum(axis=0)
        g = p[1] + 1j * p[2]
        v = np.concatenate([f / sigma, [1.0 / sigma_r]])
        b = np.concatenate([y / sigma, [z / sigma_r]])
        return float(np.linalg.norm(b - g * v) ** 2)

    p0 = np.array([t_gls, 0.9, 0.4])
    r2 = minimize(J2, p0, method="Nelder-Mead",
                  options={"maxiter": 20000, "xatol": 1e-14, "fatol": 1e-22})
    r3 = minimize(J3, p0, method="Nelder-Mead",
                  options={"maxiter": 20000, "xatol": 1e-14, "fatol": 1e-22})
    out["single_draw_profiling"] = {
        "theta_profiled_reference": t_gls, "J1_min": v_gls,
        "theta_joint_varpro": float(r2.x[0]), "J2_min": float(r2.fun),
        "theta_augmented_gls": float(r3.x[0]), "J3_min": float(r3.fun),
        "theta_plugin": t_plug, "Jplug_min": v_plug,
        "rel_diff_profiled_vs_varpro": float(abs(t_gls - r2.x[0]) / abs(t_gls)),
        "rel_diff_profiled_vs_augmented": float(abs(t_gls - r3.x[0]) / abs(t_gls)),
        "abs_diff_profiled_vs_varpro": float(abs(t_gls - r2.x[0])),
        "abs_diff_profiled_vs_augmented": float(abs(t_gls - r3.x[0])),
        "rel_diff_varpro_vs_augmented": float(abs(r2.x[0] - r3.x[0]) / abs(r2.x[0])),
        "abs_diff_profiled_vs_plugin": float(abs(t_gls - t_plug)),
        "J1_grid_vs_J2_min_rel": float(abs(v_gls - r2.fun) / abs(v_gls)),
        "J1_min_vs_J3_min_rel": float(abs(v_gls - r3.fun) / abs(v_gls)),
        "g_star_at_theta_hat": cplx(g_star_of(y, Fg[0], z, sigma, sigma_r)),
        "g_varpro": cplx(r2.x[1] + 1j * r2.x[2]),
    }
    rec.log("D3 single draw: theta_ref=%.12f theta_varpro=%.12f theta_aug=%.12f "
            "theta_plug=%.12f ; rel(ref,varpro)=%.2e rel(ref,aug)=%.2e"
            % (t_gls, r2.x[0], r3.x[0], t_plug,
               out["single_draw_profiling"]["rel_diff_profiled_vs_varpro"],
               out["single_draw_profiling"]["rel_diff_profiled_vs_augmented"]))

    # ---- Monte Carlo: bias, null calibration and chi-square gate behaviour
    # Statistic convention: T = 2 * (concentrated weighted residual), so a proper
    # complex noise of variance sigma^2 contributes 2|.|^2/sigma^2 ~ chi^2_2 per
    # complex degree of freedom (no stray factor of 1/2).
    #   * reference GLS at the TRUE theta: r = eps_y - eps_z f has covariance
    #     sigma^2 I + sigma_r^2 f f^*, so T_true ~ chi^2_{2m} exactly
    #     (2m real data dims; the 2 real gain params are eliminated in closed form).
    #   * reference GLS profiled over theta as well: 2m+2 real data, 3 real params
    #     => T_gls ~ chi^2_{2m-1} asymptotically (exact law is a curved-family
    #     mixture; the empirical quantile is reported alongside).
    #   * plug-in baseline (z treated as an exact gain): the unmodelled
    #     sigma_r^2 |f|^2 reference variance inflates T_plugin, which is NOT
    #     chi^2_{2m-1} under the null, so a naive gate over-rejects.
    m = int(C.shape[1])
    dof_gls = 2 * m - 1
    dof_true = 2 * m
    thr_gls = float(chi2.ppf(0.95, dof_gls))
    thr_true = float(chi2.ppf(0.95, dof_true))
    th_g, th_p, Tg, Tp, Tt = [], [], [], [], []
    for _ in range(n_draws):
        ey = (sigma / np.sqrt(2)) * (rng.standard_normal(m) + 1j * rng.standard_normal(m))
        ez = (sigma_r / np.sqrt(2)) * (rng.standard_normal() + 1j * rng.standard_normal())
        y = g0 * f0 + ey
        z = g0 + ez
        tg, vg = profiled_theta(J1_theta, C, W, y, z, sigma, sigma_r, -2.0, 3.0)
        tp, vp = profiled_theta(lambda *a: None, C, W, y, z, sigma, sigma_r, -2.0, 3.0)
        th_g.append(tg); th_p.append(tp)
        Tg.append(2.0 * vg); Tp.append(2.0 * vp)
        Tt.append(2.0 * gls_rhs(y, f0, z, sigma, sigma_r, "woodbury"))
    th_g = np.array(th_g); th_p = np.array(th_p)
    Tg = np.array(Tg); Tp = np.array(Tp); Tt = np.array(Tt)
    n = len(th_g)
    emp_q95 = {"gls_profiled": float(np.quantile(Tg, 0.95)),
               "plugin_profiled": float(np.quantile(Tp, 0.95)),
               "gls_true_theta": float(np.quantile(Tt, 0.95))}
    ks_true = kstest(Tt, "chi2", args=(dof_true,))
    out["monte_carlo"] = {
        "n_draws": n, "m": m,
        "statistic": "T = 2 * concentrated weighted residual",
        "dof_fixed_true_theta": dof_true,
        "dof_profiled_gls_asymptotic": dof_gls,
        "chi2_95_threshold_dof_2m": thr_true,
        "chi2_95_threshold_dof_2m_minus_1": thr_gls,
        "gls_mean_theta": float(th_g.mean()), "gls_bias": float(th_g.mean() - theta0),
        "gls_bias_se": float(th_g.std(ddof=1) / np.sqrt(n)),
        "gls_rmse": float(np.sqrt(((th_g - theta0) ** 2).mean())),
        "plugin_mean_theta": float(th_p.mean()), "plugin_bias": float(th_p.mean() - theta0),
        "plugin_bias_se": float(th_p.std(ddof=1) / np.sqrt(n)),
        "plugin_rmse": float(np.sqrt(((th_p - theta0) ** 2).mean())),
        "rmse_ratio_plugin_over_gls": float(np.sqrt(((th_p - theta0) ** 2).mean())
                                            / np.sqrt(((th_g - theta0) ** 2).mean())),
        "gls_median_abs_error": float(np.median(np.abs(th_g - theta0))),
        "plugin_median_abs_error": float(np.median(np.abs(th_p - theta0))),
        "gls_fraction_within_0p1": float((np.abs(th_g - theta0) < 0.1).mean()),
        "plugin_fraction_within_0p1": float((np.abs(th_p - theta0) < 0.1).mean()),
        "gls_minus_plugin_mean_diff": float((th_g - th_p).mean()),
        "gls_minus_plugin_diff_se": float((th_g - th_p).std(ddof=1) / np.sqrt(n)),
        "T_true_fixed_theta": {
            "mean": float(Tt.mean()), "var": float(Tt.var(ddof=1)),
            "q95": emp_q95["gls_true_theta"],
            "theory_mean": float(dof_true), "theory_var": float(2 * dof_true),
            "ks_p_vs_chi2_2m": float(ks_true.pvalue),
        },
        "T_gls_profiled": {
            "mean": float(Tg.mean()), "var": float(Tg.var(ddof=1)),
            "q95": emp_q95["gls_profiled"],
            "theory_mean": float(dof_gls), "theory_var": float(2 * dof_gls),
            "ks_p_vs_chi2_2m_minus_1": float(kstest(Tg, "chi2", args=(dof_gls,)).pvalue),
        },
        "T_plugin_profiled": {
            "mean": float(Tp.mean()), "var": float(Tp.var(ddof=1)),
            "q95": emp_q95["plugin_profiled"],
            "theory_mean_if_gain_exact": float(dof_gls),
            "ks_p_vs_chi2_2m_minus_1": float(kstest(Tp, "chi2", args=(dof_gls,)).pvalue),
            "mean_over_dof_2m_minus_1": float(Tp.mean() / dof_gls),
        },
        "empirical_95pct_quantiles": emp_q95,
        "gate_rejection_rates": {
            "gls_vs_chi2_2m_minus_1": float((Tg > thr_gls).mean()),
            "gls_vs_chi2_2m": float((Tg > thr_true).mean()),
            "gls_vs_empirical_quantile": float((Tg > emp_q95["gls_profiled"]).mean()),
            "plugin_vs_chi2_2m_minus_1": float((Tp > thr_gls).mean()),
            "plugin_vs_chi2_2m": float((Tp > thr_true).mean()),
            "plugin_vs_empirical_quantile": float((Tp > emp_q95["plugin_profiled"]).mean()),
            "true_theta_vs_chi2_2m": float((Tt > thr_true).mean()),
        },
        "rejection_rate_difference_gls_minus_plugin": float(
            (Tg > thr_gls).mean() - (Tp > thr_gls).mean()),
        "note": ("Reference GLS dof at fixed true theta = 2m (exact); profiled over "
                 "theta = 2m-1 (asymptotic). Plug-in is not chi-square under the null "
                 "because the reference variance sigma_r^2|f|^2 is unmodelled."),
    }
    rec.log("D3 MC n=%d : GLS bias=%+.3e (se %.1e) rej[chi2_2m-1]=%.4f rej[chi2_2m]=%.4f ; "
            "plugin bias=%+.3e (se %.1e) rej[chi2_2m-1]=%.4f ; mean T true=%.3f (chi2_2m=%.1f) "
            "gls=%.3f (%.1f) plug=%.3f ; KS p (T_true~chi2_2m)=%.3f"
            % (n, out["monte_carlo"]["gls_bias"], out["monte_carlo"]["gls_bias_se"],
               out["monte_carlo"]["gate_rejection_rates"]["gls_vs_chi2_2m_minus_1"],
               out["monte_carlo"]["gate_rejection_rates"]["gls_vs_chi2_2m"],
               out["monte_carlo"]["plugin_bias"], out["monte_carlo"]["plugin_bias_se"],
               out["monte_carlo"]["gate_rejection_rates"]["plugin_vs_chi2_2m_minus_1"],
               out["monte_carlo"]["T_true_fixed_theta"]["mean"], float(dof_true),
               out["monte_carlo"]["T_gls_profiled"]["mean"], float(dof_gls),
               out["monte_carlo"]["T_plugin_profiled"]["mean"],
               out["monte_carlo"]["T_true_fixed_theta"]["ks_p_vs_chi2_2m"]))
    rec.put("D3_nonlinear_profiling", out)
    rec.dump()


# ============================================================ D4: full residual derivative
def part_d4(n_dir=20):
    from mpmath import mp
    mp.dps = 60
    rng = np.random.default_rng(SEED + 3)
    sigma, sigma_r = 0.4, 0.35
    K, p, m = 4, 3, 5
    C = rng.standard_normal((K, m)) + 1j * rng.standard_normal((K, m))
    W = rng.uniform(-1.0, 1.0, size=(K, p))
    theta0 = rng.uniform(-0.5, 0.5, size=p)
    y = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    z = 0.8 - 0.3j

    def to_mp(v):
        return [mp.mpc(repr(float(a.real)), repr(float(a.imag))) if np.iscomplexobj(np.asarray(v))
                else mp.mpf(repr(float(a))) for a in np.asarray(v).ravel()]

    Cm = [[mp.mpc(repr(float(C[i, j].real)), repr(float(C[i, j].imag))) for j in range(m)]
          for i in range(K)]
    Wm = [[mp.mpf(repr(float(W[i, j]))) for j in range(p)] for i in range(K)]
    ym = [mp.mpc(repr(float(y[j].real)), repr(float(y[j].imag))) for j in range(m)]
    zm = mp.mpc(repr(float(z.real)), repr(float(z.imag)))
    sig = mp.mpf(repr(sigma))
    sigr = mp.mpf(repr(sigma_r))

    def f_mp(th):
        out = [mp.mpc(0) for _ in range(m)]
        for i in range(K):
            arg = mp.mpc(0)
            for j in range(p):
                arg += Wm[i][j] * th[j]
            e = mp.exp(mp.mpc(0, 1) * arg)
            for j in range(m):
                out[j] += Cm[i][j] * e
        return out

    def J_mp(th):
        # (m x p) complex
        Jm = [[mp.mpc(0) for _ in range(p)] for _ in range(m)]
        for i in range(K):
            arg = mp.mpc(0)
            for j in range(p):
                arg += Wm[i][j] * th[j]
            e = mp.exp(mp.mpc(0, 1) * arg)
            for j in range(m):
                for q in range(p):
                    Jm[j][q] += Cm[i][j] * mp.mpc(0, 1) * Wm[i][q] * e
        return Jm

    def v_mp(th):
        f = f_mp(th)
        return [f[j] / sig for j in range(m)] + [mp.mpf(1) / sigr]

    def r_mp(th):
        v = v_mp(th)
        b = [ym[j] / sig for j in range(m)] + [zm / sigr]
        vv = sum((abs(v[j]) ** 2 for j in range(m + 1)), mp.mpf(0))
        vb = sum((v[j].conjugate() * b[j] for j in range(m + 1)), mp.mpc(0))
        g = vb / vv
        return [b[j] - v[j] * g for j in range(m + 1)]

    def analytic_dr(th, dvec):
        v = v_mp(th)
        b = [ym[j] / sig for j in range(m)] + [zm / sigr]
        Jm = J_mp(th)
        df = [sum((Jm[j][q] * dvec[q] for q in range(p)), mp.mpc(0)) for j in range(m)]
        dv = [df[j] / sig for j in range(m)] + [mp.mpf(0)]
        vv = sum((abs(v[j]) ** 2 for j in range(m + 1)), mp.mpf(0))
        vb = sum((v[j].conjugate() * b[j] for j in range(m + 1)), mp.mpc(0))
        g = vb / vv
        r = [b[j] - v[j] * g for j in range(m + 1)]
        dvc = sum((dv[j].conjugate() * r[j] for j in range(m + 1)), mp.mpc(0))
        # Q_v dv  with Q_v = I - v v^*/(v^*v)
        vcdv = sum((v[j].conjugate() * dv[j] for j in range(m + 1)), mp.mpc(0))
        Qdv = [dv[j] - v[j] * vcdv / vv for j in range(m + 1)]
        dr = [-Qdv[j] * g - v[j] * (dvc / vv) for j in range(m + 1)]
        dr_naive = [-Qdv[j] * g for j in range(m + 1)]
        vdotr = sum((v[j].conjugate() * r[j] for j in range(m + 1)), mp.mpc(0))
        return dr, dr_naive, vdotr

    th0 = [mp.mpf(repr(float(t))) for t in theta0]
    rows = []
    for _ in range(n_dir):
        d = rng.standard_normal(p)
        d /= np.linalg.norm(d)
        dm = [mp.mpf(repr(float(v))) for v in d]
        h = mp.mpf("1e-15")
        rp = r_mp([th0[j] + h * dm[j] for j in range(p)])
        rm_ = r_mp([th0[j] - h * dm[j] for j in range(p)])
        fd = [(rp[j] - rm_[j]) / (2 * h) for j in range(m + 1)]
        dr, dr_naive, vdotr = analytic_dr(th0, dm)
        def relerr(a, b):
            num = mp.sqrt(sum((abs(a[j] - b[j]) ** 2 for j in range(m + 1)), mp.mpf(0)))
            den = mp.sqrt(sum((abs(b[j]) ** 2 for j in range(m + 1)), mp.mpf(0)))
            return float(num / den) if den > 0 else float(num)
        rows.append({
            "direction": jsonable(d),
            "rel_err_full_derivative": relerr(dr, fd),
            "rel_err_naive_drop_term": relerr(dr_naive, fd),
            "abs_vdotr": float(abs(vdotr)),
            "norm_fd": float(mp.sqrt(sum((abs(fd[j]) ** 2 for j in range(m + 1)), mp.mpf(0)))),
            "norm_drop_term": float(mp.sqrt(sum((abs(dr[j] - dr_naive[j]) ** 2
                                                 for j in range(m + 1)), mp.mpf(0)))),
        })
    out = {
        "m": m, "p": p, "K": K, "sigma": sigma, "sigma_r": sigma_r,
        "theta0": jsonable(theta0), "mpmath_dps": 60, "fd_step": 1e-15,
        "Q_v_used": "I - v v^*/(v^* v)  (orthogonal projector onto the complement of v)",
        "max_rel_err_full_derivative": max(r["rel_err_full_derivative"] for r in rows),
        "median_rel_err_full_derivative": float(np.median([r["rel_err_full_derivative"] for r in rows])),
        "min_rel_err_full_derivative": min(r["rel_err_full_derivative"] for r in rows),
        "max_rel_err_naive_drop_term": max(r["rel_err_naive_drop_term"] for r in rows),
        "median_rel_err_naive_drop_term": float(np.median([r["rel_err_naive_drop_term"] for r in rows])),
        "max_abs_v_star_r": max(r["abs_vdotr"] for r in rows),
        "drop_term_is_negligible": bool(max(r["rel_err_naive_drop_term"] for r in rows) < 1e-3),
        "rows": rows,
    }
    rec.log("D4 full derivative: max rel err %.3e (median %.3e) ; naive drop-term max rel "
            "err %.3e ; max |v^*r| %.3e"
            % (out["max_rel_err_full_derivative"], out["median_rel_err_full_derivative"],
               out["max_rel_err_naive_drop_term"], out["max_abs_v_star_r"]))
    rec.put("D4_residual_derivative", out)
    rec.dump()


def main():
    global rec
    rec = Rec()
    rec.log("PART D start  json=%s" % os.path.relpath(rec.json_path, ROOT))
    part_d1()
    part_d2()
    part_d3()
    part_d4()
    rec.log("PART D done  json=%s  wall=%.1fs"
            % (os.path.relpath(rec.json_path, ROOT), time.time() - WALL0))
    rec.dump()


if __name__ == "__main__":
    main()
