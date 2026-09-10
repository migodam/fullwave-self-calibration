#!/usr/bin/env python3
"""PART C -- A5 section 3 (fixed-loss Born gain ambiguity) independent checks.

Theory read directly from inputs/derive/A5_sec_03.md:
    y = g B (u + i*ell),   B in C^{m x 2} of full complex column rank,
    u in R^2 the two-region real contrast, ell in R^2 known exactly, g != 0.
  * if D = det[u, ell] != 0 then (u, g) is uniquely determined by noiseless data;
    if ell == 0 a positive real scale ambiguity survives;
    if ell != 0 and u = q*ell then every allowed q' gives the same data.
  * w = B^dagger y = g (u + i ell),  h = 1/g = a + i b,
    M(w) = [Im w, Re w]  (2x2, columns Im w and Re w),  M(w)[a;b] = ell,
    det M(w) = -|g|^2 det[u, ell],  u = Re(h w),  g = 1/h.
  * H = I - (u u^T + ell ell^T)/S,  S = |u|^2 + |ell|^2,
    lambda_min(H) = (1 - sqrt(1 - 4 D^2/S^2))/2,
    sigma_min(A_vis) >= |g| sigma_min(B) sqrt(lambda_min(H)).
  * perturbation bounds (||delta w|| <= eta, mu = sigma_min(M(w)), eta < mu):
    |delta h| <= |h| eta/(mu-eta),
    ||delta u|| <= |h| eta [1 + (||w||+eta)/(mu-eta)],
    eta = ||B^dagger|| beta  for raw data error ||delta y|| <= beta.

No installs, single thread, fixed seed reported in the JSON, incremental output.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.optimize import least_squares, minimize

SEED = 20260910
ELL = np.array([0.03, 0.05], dtype=float)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WALL0 = time.time()


# ---------------------------------------------------------------------------- helpers
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
        base = os.path.join(ROOT, "results", "fixed_loss_born_%s" % stamp)
        n = 0
        while os.path.exists(base + (".json" if n == 0 else "_%d.json" % n)):
            n += 1
        self.json_path = base + (".json" if n == 0 else "_%d.json" % n)
        self.log_path = os.path.join(ROOT, "logs", "fixed_loss_born.log")
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
                "ell": jsonable(ELL),
            }
        )
        tmp = self.json_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(jsonable(self.data), fh, indent=1)
        os.replace(tmp, self.json_path)


def det2(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def M_of(w):
    return np.column_stack([w.imag, w.real])


def reconstruct(w, ell=ELL):
    M = M_of(w)
    ab = np.linalg.solve(M, ell)
    h = ab[0] + 1j * ab[1]
    return h, (h * w).real, 1.0 / h, M


def reconstruct_batch(w, ell=ELL):
    n = w.shape[0]
    M = np.empty((n, 2, 2))
    M[:, 0, 0] = w[:, 0].imag
    M[:, 1, 0] = w[:, 1].imag
    M[:, 0, 1] = w[:, 0].real
    M[:, 1, 1] = w[:, 1].real
    bb = np.tile(np.asarray(ell, dtype=float), (n, 1))[:, :, None]
    ab = np.linalg.solve(M, bb)[..., 0]
    h = ab[:, 0] + 1j * ab[:, 1]
    return h, (h[:, None] * w).real, 1.0 / h, M


def make_B(rng, m, cond=None):
    if cond is None:
        return rng.standard_normal((m, 2)) + 1j * rng.standard_normal((m, 2))
    U, _ = np.linalg.qr(rng.standard_normal((m, 2)) + 1j * rng.standard_normal((m, 2)))
    V, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)))
    s = np.array([1.0, 1.0 / cond])
    return (U * s) @ V.conj().T


def draw_u(rng, lo=0.1, hi=2.5, det_min=1e-2):
    while True:
        u = rng.uniform(lo, hi, size=2)
        if abs(det2(u, ELL)) > det_min:
            return u


def lam_min_H(u, ell=ELL):
    S = float(u @ u + ell @ ell)
    D = det2(u, ell)
    H = np.eye(2) - (np.outer(u, u) + np.outer(ell, ell)) / S
    return float(np.linalg.eigvalsh(H)[0]), S, D, H


def lam_min_H_formula(u, ell=ELL, dps=50):
    """Literal A5 formula (1 - sqrt(1 - 4 D^2/S^2))/2 evaluated in mpmath."""
    from mpmath import mp
    mp.dps = dps
    uu = [mp.mpf(repr(float(v))) for v in u]
    ee = [mp.mpf(repr(float(v))) for v in ell]
    S = uu[0] ** 2 + uu[1] ** 2 + ee[0] ** 2 + ee[1] ** 2
    D = uu[0] * ee[1] - uu[1] * ee[0]
    return float((1 - mp.sqrt(1 - 4 * D ** 2 / S ** 2)) / 2)


def mu_pred_stable(u, g, ell=ELL):
    """sigma_min(M(w)) = |g| sqrt(S lambda_min(H)); stable rearrangement.

    (S - sqrt(S^2-4D^2))/2 suffers catastrophic cancellation for small D;
    the algebraically identical 2(S t)^2... form below is used instead.
    """
    S = float(u @ u + ell @ ell)
    D = det2(u, ell)
    t = D / S
    lam = 2 * t * t / (1 + np.sqrt(1 - 4 * t * t))
    return float(abs(g) * np.sqrt(S * lam))


rec = None


# ============================================================ C1: B generation
def part_c1():
    rng = np.random.default_rng(SEED)
    out = {}
    for m in (6, 12):
        for kind, cond in (("well_conditioned", None), ("ill_conditioned_1e4", 1e4)):
            B = make_B(rng, m, cond)
            sv = np.linalg.svd(B, compute_uv=False)
            key = "m%d_%s" % (m, kind)
            out[key] = {
                "m": m,
                "requested_cond": cond,
                "singular_values": jsonable(sv),
                "cond": float(sv[0] / sv[-1]),
                "rank": int(np.linalg.matrix_rank(B)),
                "B": jsonable(B),
            }
            rec.log("C1 %-28s sigma=%s cond=%.6g rank=%d"
                    % (key, np.array2string(sv, precision=5), sv[0] / sv[-1],
                       np.linalg.matrix_rank(B)))
    rec.put("C1_B_generation", out)
    rec.dump()


# ============================================ C2: exact reconstruction + identity
def nl_fit_from(B, y, start, ell=ELL):
    def resid(p):
        r = (p[2] + 1j * p[3]) * (B @ (p[:2] + 1j * ell)) - y
        return np.concatenate([r.real, r.imag])

    sol = least_squares(resid, start, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    return sol


def part_c2(n_draws=200, n_nl=10):
    rng = np.random.default_rng(SEED + 1)
    out = {}
    for m in (6, 12):
        for kind, cond in (("well_conditioned", None), ("ill_conditioned_1e4", 1e4)):
            key = "m%d_%s" % (m, kind)
            rel_u, rel_g, det_abs, det_rel, nl_u, nl_g, nls = [], [], [], [], [], [], []
            nl_cost, nl_global = [], []
            N_START = 20
            for k in range(n_draws):
                B = make_B(rng, m, cond)
                u = draw_u(rng)
                rho = rng.uniform(0.75, 1.25)
                g = rho * np.exp(1j * rng.uniform(0, 2 * np.pi))
                y = g * (B @ (u + 1j * ELL))
                w = np.linalg.lstsq(B, y, rcond=None)[0]
                h, u_hat, g_hat, M = reconstruct(w)
                rel_u.append(float(np.linalg.norm(u_hat - u) / np.linalg.norm(u)))
                rel_g.append(float(abs(g_hat - g) / abs(g)))
                D = det2(u, ELL)
                dm = float(np.linalg.det(M))
                det_abs.append(abs(dm + abs(g) ** 2 * D))
                det_rel.append(abs(dm + abs(g) ** 2 * D) / abs(abs(g) ** 2 * D))
                if k < n_nl:
                    # multipstart: a local solver from a random start can stall in a
                    # spurious minimum (the problem is unique, not convex); report both
                    # the best-residual solution and how many starts reached it.
                    best, nglob = None, 0
                    for _ in range(N_START):
                        start = np.concatenate([rng.uniform(-1, 1, 2),
                                                rng.uniform(-1, 1, 2)])
                        sol = nl_fit_from(B, y, start)
                        if sol.cost < 1e-18:
                            nglob += 1
                        if best is None or sol.cost < best.cost:
                            best = sol
                    nl_u.append(float(np.linalg.norm(best.x[:2] - u) / np.linalg.norm(u)))
                    nl_g.append(float(abs((best.x[2] + 1j * best.x[3]) - g) / abs(g)))
                    nl_cost.append(float(best.cost))
                    nl_global.append(nglob)
                    nls.append(int(best.status))
            out[key] = {
                "n_draws": n_draws,
                "max_rel_u_err": max(rel_u),
                "median_rel_u_err": float(np.median(rel_u)),
                "max_rel_g_err": max(rel_g),
                "median_rel_g_err": float(np.median(rel_g)),
                "max_det_identity_abs_err": max(det_abs),
                "max_det_identity_rel_err": max(det_rel),
                "n_nl_fits": n_nl,
                "nl_starts_per_fit": N_START,
                "max_nl_ls_rel_u_diff": max(nl_u),
                "max_nl_ls_rel_g_diff": max(nl_g),
                "max_nl_ls_best_cost": max(nl_cost),
                "n_starts_reaching_global_per_fit": nl_global,
                "nl_ls_status": nls,
            }
            rec.log("C2 %-28s max|du|/|u|=%.3e max|dg|/|g|=%.3e det_rel=%.3e "
                    "nl_best=%.3e nl_global_starts=%s"
                    % (key, max(rel_u), max(rel_g), max(det_rel), max(nl_u), nl_global))
    rec.put("C2_exact_reconstruction", out)
    rec.dump()


# ============================================ C3: uniqueness + scale ambiguity
def part_c3(n_dir=20000):
    rng = np.random.default_rng(SEED + 2)
    out = {}
    B = make_B(rng, 8)
    u = draw_u(rng)
    g = 1.0 + 0.0j
    w = g * (B @ (u + 1j * ELL))

    def y_of(uu, gg):
        return gg * (B @ (uu + 1j * ELL))

    eps_p = 1e-3
    J = np.zeros((2 * B.shape[0], 4))
    for j, e in enumerate(np.eye(4)):
        du = e[:2] * eps_p
        dg = (e[2] + 1j * e[3]) * eps_p
        d = y_of(u + du, g + dg) - w
        J[:, j] = np.concatenate([d.real, d.imag]) / eps_p
    sv = np.linalg.svd(J, compute_uv=False)
    dirs = rng.standard_normal((n_dir, 4))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    steps = dirs * eps_p
    dys = np.array([np.linalg.norm(y_of(u + s[:2], g + s[2] + 1j * s[3]) - w)
                    for s in steps])
    out["noiseless_uniqueness"] = {
        "m": int(B.shape[0]),
        "step_norm": eps_p,
        "n_directions": n_dir,
        "min_delta_y": float(dys.min()),
        "median_delta_y": float(np.median(dys)),
        "max_delta_y": float(dys.max()),
        "jacobian_singular_values": jsonable(sv),
        "predicted_min_from_sigma_min_J": float(eps_p * sv[-1]),
        "relative_margin_min_over_median": float(dys.min() / np.median(dys)),
        "strictly_positive": bool(dys.min() > 0),
    }
    rec.log("C3 uniqueness: min||dy||=%.6e  sigma_min(J)*1e-3=%.6e (median %.6e)"
            % (dys.min(), eps_p * sv[-1], np.median(dys)))

    u0 = np.array([0.8, 1.3])
    g0 = 0.9 + 0.3j
    Bp = make_B(rng, 7)
    y0 = g0 * (Bp @ u0)
    amb = {}
    for c in (0.3, 2.0, 7.0):
        gp = g0 / c
        up = c * u0
        yp = gp * (Bp @ up)
        amb["c=%g" % c] = {
            "c": c, "g_prime": cplx(gp), "u_prime": jsonable(up),
            "rel_residual": float(np.linalg.norm(yp - y0) / np.linalg.norm(y0)),
        }
    out["ell_zero_scale_ambiguity"] = {
        "u": jsonable(u0), "g": cplx(g0), "m": int(Bp.shape[0]),
        "det_u_ell": det2(u0, np.zeros(2)),
        "cases": amb,
        "all_identical_to_roundoff": bool(all(v["rel_residual"] < 1e-14 for v in amb.values())),
    }
    rec.log("C3 ell=0 scale ambiguity: max rel residual %.3e over c=(0.3,2,7)"
            % max(v["rel_residual"] for v in amb.values()))

    Det = det2(u, ELL)
    wc = g * (u + 1j * ELL)          # compressed data B^dagger y, 2-vector
    Mw = M_of(wc)
    out["det_nondegenerate"] = {"det_u_ell": Det, "det_M_w": float(np.linalg.det(Mw)),
                                "det_M_plus_g2D": float(np.linalg.det(Mw) + abs(g) ** 2 * Det)}
    rec.put("C3_uniqueness_and_ambiguity", out)
    rec.dump()


# ============================================ C4: parallel counterexample
def part_c4():
    with open(os.path.join(ROOT, "inputs", "modal_boundary.json")) as fh:
        mb = json.load(fh)
    allowed = mb["gain_allowed_modulus"]
    q, qp = 30.0, 45.0
    u = q * ELL
    up = qp * ELL
    c = (qp + 1j) / (q + 1j)
    g, gp = 1.0 + 0.0j, 1.0 / c
    rng = np.random.default_rng(SEED + 3)
    B = make_B(rng, 9)
    y1 = g * (B @ (u + 1j * ELL))
    y2 = gp * (B @ (up + 1j * ELL))
    relres = float(np.linalg.norm(y1 - y2) / np.linalg.norm(y1))
    quoted = 0.666831 + 0.007404j
    eps1 = 1 + u
    eps2 = 1 + up
    lam = 1.2
    y3 = (lam * g) * (B @ (u + 1j * ELL))
    y4 = (lam * gp) * (B @ (up + 1j * ELL))
    out = {
        "ell": jsonable(ELL), "q": q, "q_prime": qp,
        "u": jsonable(u), "u_prime": jsonable(up),
        "c": cplx(c), "g": cplx(g), "g_prime": cplx(gp),
        "g_prime_vs_quoted": {
            "quoted": cplx(quoted),
            "abs_diff": float(abs(gp - quoted)),
            "rel_diff": float(abs(gp - quoted) / abs(quoted)),
            "matches_quoted_1e_6": bool(abs(gp - quoted) < 1e-6),
        },
        "m": int(B.shape[0]),
        "same_data_rel_residual": relres,
        "same_data_to_1e_15": bool(relres < 1e-15),
        "eps_world1": jsonable(eps1),
        "eps_world2": jsonable(eps2),
        "gain_allowed_modulus": allowed,
        "abs_g": float(abs(g)),
        "abs_g_prime": float(abs(gp)),
        "g_in_annulus": bool(allowed[0] <= abs(g) <= allowed[1]),
        "g_prime_in_annulus": bool(allowed[0] <= abs(gp) <= allowed[1]),
        "abs_g_prime_over_g": float(abs(gp) / abs(g)),
        "rescaled_pair_rel_residual": float(np.linalg.norm(y3 - y4) / np.linalg.norm(y3)),
        "shared_gain_rescale": {
            "lambda_min_needed": float(allowed[0] / abs(gp)),
            "lambda_max_allowed": float(allowed[1]),
            "example_lambda": lam,
            "then_abs_g": float(lam),
            "then_abs_g_prime": float(lam * abs(gp)),
            "then_both_in_annulus": bool(allowed[0] <= lam <= allowed[1]
                                         and allowed[0] <= lam * abs(gp) <= allowed[1]),
        },
    }
    rec.log("C4 parallel: rel residual %.3e ; |g'|=%.9f ; g' in annulus: %s"
            % (relres, abs(gp), out["g_prime_in_annulus"]))
    rec.put("C4_parallel_counterexample", out)
    rec.dump()


# ============================================ C5: near-parallel deterioration
def _err_of(X, Bdag, w, u, h):
    dw = X @ Bdag.T
    hh, uu, gg, _ = reconstruct_batch(w + dw)
    return np.linalg.norm(uu - u[None, :], axis=1), np.abs(hh - h), np.abs(gg - 1.0 / h)


def adversarial(rng, beta, Bdag, w, u, h, target, n_sample=8000, n_ref=2, maxiter=1200):
    m = Bdag.shape[1]
    X = rng.standard_normal((n_sample, m)) + 1j * rng.standard_normal((n_sample, m))
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    X *= beta
    eu, eh, eg = _err_of(X, Bdag, w, u, h)
    key = {"u": eu, "h": eh, "g": eg}[target]
    order = np.argsort(-key)
    best = float(key[order[0]])
    bx = X[order[0]].copy()

    def negobj(z):
        d = z[:m] + 1j * z[m:]
        n = np.linalg.norm(d)
        if n == 0:
            return 1e30
        _, e_h, e_g = _err_of((d / n * beta)[None, :], Bdag, w, u, h)
        return -float(e_h[0] if target == "h" else e_g[0])

    for k in range(min(n_ref, n_sample)):
        res = minimize(negobj, np.concatenate([X[order[k]].real, X[order[k]].imag]),
                       method="Nelder-Mead",
                       options={"maxiter": maxiter, "xatol": 1e-15, "fatol": 1e-22})
        if -res.fun > best:
            best = -res.fun
            d = res.x[:m] + 1j * res.x[m:]
            bx = d / np.linalg.norm(d) * beta
    return bx, best


def part_c5():
    rng = np.random.default_rng(SEED + 4)
    q = 30.0
    perp = np.array([-ELL[1], ELL[0]])
    s_grid = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
    out = {}
    for tag, m, cond in (("B_identity", 2, None), ("B_random_m6", 6, None),
                         ("B_illcond_m6_1e4", 6, 1e4)):
        B = np.eye(2, dtype=complex) if (cond is None and m == 2) else make_B(rng, m, cond)
        Bdag = np.linalg.pinv(B)
        nm = float(np.linalg.svd(Bdag, compute_uv=False)[0])
        rows = []
        s_small = s_grid[-1]
        u_s = q * ELL + s_small * perp
        mu_small = float(np.linalg.svd(M_of(u_s + 1j * ELL), compute_uv=False)[-1])
        beta = 0.1 * mu_small / nm
        for s in s_grid:
            u = q * ELL + s * perp
            g = 1.0 + 0.0j
            w = g * (u + 1j * ELL)
            M = M_of(w)
            mu = float(np.linalg.svd(M, compute_uv=False)[-1])
            D = det2(u, ELL)
            S = float(u @ u + ELL @ ELL)
            mu_pre = mu_pred_stable(u, g)
            lam, _, _, _ = lam_min_H(u)
            eta = nm * beta
            h = 1.0 / g
            bnd_h = abs(h) * eta / (mu - eta)
            bnd_u = abs(h) * eta * (1 + (np.linalg.norm(w) + eta) / (mu - eta))
            Xr = rng.standard_normal((4000, m)) + 1j * rng.standard_normal((4000, m))
            Xr /= np.linalg.norm(Xr, axis=1, keepdims=True)
            Xr *= beta
            eru, erh, erg = _err_of(Xr, Bdag, w, u, h)
            _, euh = adversarial(rng, beta, Bdag, w, u, h, "u")
            _, ehh = adversarial(rng, beta, Bdag, w, u, h, "h")
            rows.append({
                "s": s, "abs_D": abs(D), "mu": mu, "mu_pred": mu_pre,
                "mu_pred_rel_err": abs(mu - mu_pre) / abs(mu_pre),
                "lambda_min_H": lam, "beta": beta, "eta": eta,
                "eta_over_mu": eta / mu, "bound_dh": bnd_h, "bound_du": bnd_u,
                "random_mean_du": float(eru.mean()), "random_max_du": float(eru.max()),
                "random_max_dh": float(erh.max()), "random_max_dg": float(erg.max()),
                "adv_du": float(euh), "adv_dh": float(ehh),
                "ratio_random_du": float(eru.max() / bnd_u),
                "ratio_random_dh": float(erh.max() / bnd_h),
                "ratio_adv_du": float(euh / bnd_u),
                "ratio_adv_dh": float(ehh / bnd_h),
            })
            rec.log("C5 %-18s s=%.0e |D|=%.3e mu=%.3e adv_du=%.3e bnd_u=%.3e "
                    "ratio_adv=%.4f ratio_rand=%.4f"
                    % (tag, s, abs(D), mu, euh, bnd_u, euh / bnd_u, eru.max() / bnd_u))

        def slope(key):
            x = np.log(np.array([r["abs_D"] for r in rows]))
            y = np.log(np.array([r[key] for r in rows]))
            return float(np.polyfit(x, y, 1)[0])

        y_round = g * (B @ (u_s + 1j * ELL))
        floor = float(np.linalg.norm(Bdag @ y_round - g * (u_s + 1j * ELL)))
        out[tag] = {
            "m": m, "requested_cond": cond,
            "cond_B": float(np.linalg.svd(B, compute_uv=False)[0]
                            / np.linalg.svd(B, compute_uv=False)[-1]),
            "Bdag_operator_norm": nm, "beta": beta,
            "pinv_apply_roundoff_floor": floor,
            "slope_log_adv_du_vs_log_absD": slope("adv_du"),
            "slope_log_random_mean_du_vs_log_absD": slope("random_mean_du"),
            "slope_log_adv_dh_vs_log_absD": slope("adv_dh"),
            "max_ratio_random_du": max(r["ratio_random_du"] for r in rows),
            "max_ratio_random_dh": max(r["ratio_random_dh"] for r in rows),
            "max_ratio_adv_du": max(r["ratio_adv_du"] for r in rows),
            "max_ratio_adv_dh": max(r["ratio_adv_dh"] for r in rows),
            "all_ratios_le_1": bool(all(r["ratio_random_du"] <= 1 + 1e-9
                                        and r["ratio_random_dh"] <= 1 + 1e-9
                                        and r["ratio_adv_du"] <= 1 + 1e-9
                                        and r["ratio_adv_dh"] <= 1 + 1e-9 for r in rows)),
            "max_ratio_over_all_measures": max(
                max(r["ratio_random_du"], r["ratio_random_dh"],
                    r["ratio_adv_du"], r["ratio_adv_dh"]) for r in rows),
            "max_relative_overshoot_over_1": max(
                0.0, max(max(r["ratio_random_du"], r["ratio_random_dh"],
                             r["ratio_adv_du"], r["ratio_adv_dh"]) for r in rows) - 1.0),
            "ratio_tolerance_used": 1e-9,
            "max_mu_pred_rel_err": max(r["mu_pred_rel_err"] for r in rows),
            "rows": rows,
        }
        rec.log("C5 %-18s slope(adv du)=%.4f slope(rand du)=%.4f ; max ratios "
                "rand_du=%.4f rand_dh=%.4f adv_du=%.4f adv_dh=%.4f ; all<=1: %s"
                % (tag, out[tag]["slope_log_adv_du_vs_log_absD"],
                   out[tag]["slope_log_random_mean_du_vs_log_absD"],
                   out[tag]["max_ratio_random_du"], out[tag]["max_ratio_random_dh"],
                   out[tag]["max_ratio_adv_du"], out[tag]["max_ratio_adv_dh"],
                   out[tag]["all_ratios_le_1"]))
        rec.put("C5_near_parallel", out)
        rec.dump()


# ============================================ C6: exact local lower bound
def sigma_min_real(B, g, u, ell=ELL, extra_nuisance=False, sigma=1.0, rng=None):
    m = B.shape[0]
    scale = np.sqrt(2.0) / sigma
    F = B @ (u + 1j * ell)

    def rv(v, sc):
        return np.concatenate([v.real, v.imag]) * sc

    A = np.column_stack([rv(g * B[:, 0], scale), rv(g * B[:, 1], scale)])
    Ncols = [rv(F, scale), rv(1j * F, scale)]
    if extra_nuisance:
        hh = rng.standard_normal(m) + 1j * rng.standard_normal(m)
        Ncols += [rv(g * hh, scale), rv(1j * g * hh, scale)]
    N = np.column_stack(Ncols)
    Qn, _ = np.linalg.qr(N)
    Avis = (np.eye(2 * m) - Qn @ Qn.T) @ A
    return float(np.linalg.svd(Avis, compute_uv=False)[-1]), Avis


def part_c6(n_lambda=20000):
    rng = np.random.default_rng(SEED + 5)
    out = {}
    errs = []
    for _ in range(n_lambda):
        u = draw_u(rng, -2.0, 3.0)
        lam, S, D, H = lam_min_H(u)
        pred = lam_min_H_formula(u)
        errs.append(abs(lam - pred))
    out["lambda_min_H"] = {
        "n_draws": n_lambda, "max_abs_err": float(max(errs)),
        "median_abs_err": float(np.median(errs)),
        "formula": "(1 - sqrt(1 - 4 D^2/S^2))/2, S = |u|^2+|ell|^2, D = det[u,ell]",
        "evaluated_with": "mpmath dps=50 literal expression; compared to numpy eigvalsh",
    }
    rec.log("C6 lambda_min(H) formula: max abs err %.3e over %d draws" % (max(errs), n_lambda))

    blocks = {}
    for tag, m, cond in (("B_identity_m2", 2, None), ("B_random_m8", 8, None),
                         ("B_random_m12", 12, None), ("B_illcond_m8_1e4", 8, 1e4)):
        B = np.eye(2, dtype=complex) if (cond is None and m == 2) else make_B(rng, m, cond)
        smin_B = float(np.linalg.svd(B, compute_uv=False)[-1])
        rows = []
        for _ in range(30):
            u = draw_u(rng)
            g = rng.uniform(0.75, 1.25) * np.exp(1j * rng.uniform(0, 2 * np.pi))
            sm_b = float(np.linalg.svd(B, compute_uv=False)[-1])
            lam, S, D, _ = lam_min_H(u)
            scale = np.sqrt(2.0)
            lhs, _ = sigma_min_real(B, g, u)
            rhs = scale * abs(g) * sm_b * np.sqrt(lam)
            rows.append({"lhs": lhs, "rhs": rhs,
                         "ratio": (lhs / rhs if rhs > 0 else None),
                         "smin_B_complex": sm_b, "lambda_min_H": lam})
        ratios = [r["ratio"] for r in rows if r["ratio"] is not None]
        blocks[tag] = {
            "m": m, "cond_B": float(np.linalg.svd(B, compute_uv=False)[0] / smin_B),
            "n_draws": len(rows), "min_ratio": min(ratios), "max_ratio": max(ratios),
            "median_ratio": float(np.median(ratios)),
            "inequality_holds_all": bool(min(ratios) >= 1 - 1e-9),
            "max_relative_undershoot_below_1": max(0.0, 1.0 - min(ratios)),
            "ratio_tolerance_used": 1e-9,
            "tight_within_1e_6": bool(abs(np.median(ratios) - 1) < 1e-6),
            "rows_head": rows[:5],
        }
        rec.log("C6 %-18s sigma_min(Avis)/bound in [%.6f, %.6f] median %.6f holds=%s"
                % (tag, min(ratios), max(ratios), float(np.median(ratios)),
                   blocks[tag]["inequality_holds_all"]))
    out["sigma_min_Avis_bound"] = blocks

    B = make_B(rng, 8)
    u = draw_u(rng)
    g = 1.0 + 0.0j
    smin_B = float(np.linalg.svd(B, compute_uv=False)[-1])
    lam, S, D, _ = lam_min_H(u)
    lhs0, _ = sigma_min_real(B, g, u, extra_nuisance=False)
    lhs1, _ = sigma_min_real(B, g, u, extra_nuisance=True, rng=rng)
    rhs = np.sqrt(2.0) * abs(g) * smin_B * np.sqrt(lam)
    out["extra_geometry_nuisance_probe"] = {
        "note": ("N = [D_x(gF), F, iF] with a NONZERO position/geometry nuisance "
                 "direction can only shrink sigma_min(A_vis); the quoted bound carries "
                 "no such term, so it may fail. Labelled probe, not the prescribed check."),
        "sigma_min_Avis_no_extra": lhs0, "sigma_min_Avis_with_extra": lhs1,
        "bound": rhs, "ratio_without_extra": lhs0 / rhs, "ratio_with_extra": lhs1 / rhs,
        "fails_with_extra": bool(lhs1 / rhs < 1),
    }
    rec.log("C6 extra-nuisance probe: ratio %.6f (no extra) vs %.6f (with extra), fails=%s"
            % (lhs0 / rhs, lhs1 / rhs, lhs1 / rhs < 1))
    rec.put("C6_lower_bound", out)
    rec.dump()


def main():
    global rec
    rec = Rec()
    rec.log("PART C start  json=%s" % os.path.relpath(rec.json_path, ROOT))
    part_c1()
    part_c2()
    part_c3()
    part_c4()
    part_c5()
    part_c6()
    rec.log("PART C done  json=%s  wall=%.1fs"
            % (os.path.relpath(rec.json_path, ROOT), time.time() - WALL0))
    rec.dump()


if __name__ == "__main__":
    main()
