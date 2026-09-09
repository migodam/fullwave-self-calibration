"""Bounded full-wave 2D TM pilot on the 2001 Fresnel dielTM_dec8f data.

Forward model: analytic Mie series for a homogeneous lossless/lossy circular
dielectric cylinder (radius a, permittivity eps_r) illuminated by a plane wave
travelling from the emitter position towards the cylinder centre. The
view-dependent propagation phase exp(+i k0 |r_s - c|) is applied explicitly
(the time convention exp(+i w t) was confirmed against the measured incident
phase coherence); a single complex gain per frequency is profiled out
identically for every fit variant ("gain per frequency" calibration nuisance).

Units: the .exp files carry no stated field units; the per-frequency complex
gain absorbs the absolute scale/phase, so the fit is scale-free.

Self-verification of the analytic core (all asserted in self_checks()):
  1. Jacobi-Anger expansion of the incident plane wave (numerical check),
  2. E_z and dE_z/drho continuity across the cylinder boundary,
  3. net Poynting flux through an enclosing circle ~ 0 (lossless) / < 0 (lossy),
  4. PEC limit of the scattering coefficients,
  5. series truncation tail decay.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import hankel1, jv
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
C0 = 299792458.0
DE, DR = 0.720, 0.760  # emitter and receiver distances from the setup centre [m]
A_RAD = 0.015          # published cylinder radius [m]
MU0 = 4.0e-7 * np.pi
_GEOM = None


def cached_geometry():
    """Return (th_s, th_r, src, rec), computed once."""
    global _GEOM
    if _GEOM is None:
        _GEOM = geometry()
    return _GEOM


def load_data():
    """Load dielTM_dec8f.exp into structured arrays."""
    p = (
        ROOT.parent
        / "a2_literature"
        / "data"
        / "2001_iop_17_6_301"
        / "dielTM_dec8f.exp"
    )
    rows = []
    with open(p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            rows.append([float(x) for x in s.split()])
    arr = np.asarray(rows)
    view = arr[:, 0].astype(int)
    rec = arr[:, 1].astype(int)
    freq = arr[:, 2].astype(int)
    total = arr[:, 3] + 1j * arr[:, 4]
    incident = arr[:, 5] + 1j * arr[:, 6]
    return view, rec, freq, total, incident


def geometry():
    """Source and receiver positions for the file's (view, receiver) indices."""
    th_s = np.deg2rad((np.arange(1, 37) - 1) * 10.0)      # sources, 10 deg steps
    th_r = np.deg2rad((np.arange(1, 73) - 1) * 5.0)       # receivers, 5 deg steps
    src = np.stack([DE * np.cos(th_s), DE * np.sin(th_s)], axis=1)
    rec = np.stack([DR * np.cos(th_r), DR * np.sin(th_r)], axis=1)
    return th_s, th_r, src, rec


def jacobi_anger(k0: float, u: np.ndarray, nmax: int, a: float = A_RAD) -> np.ndarray:
    """Incident plane-wave expansion weights w_n around the cylinder centre.

    E_i(r) = exp(-i k0 u.(r-c)) = sum_n w_n J_n(k0 rho) exp(i n phi).
    With u at polar angle phi_u: w_n = (-i)^n exp(-i n phi_u).
    """
    phi_u = np.arctan2(u[1], u[0])
    n = np.arange(-nmax, nmax + 1)
    return ((-1j) ** n) * np.exp(-1j * n * phi_u)


def mie_coefficients(k0: float, eps: complex, a: float, nmax: int) -> np.ndarray:
    """TM scattering coefficients c_n for a dielectric circular cylinder.

    a_n = w_n * c_n.  c_n from continuity of E_z and dE_z/drho at rho = a.
    """
    k1 = k0 * np.sqrt(eps)
    n = np.arange(-nmax, nmax + 1)
    j0 = jv(n, k0 * a)
    j1 = jv(n, k1 * a)
    h0 = hankel1(n, k0 * a)
    dj0 = (jv(n - 1, k0 * a) - jv(n + 1, k0 * a)) / 2.0
    dj1 = (jv(n - 1, k1 * a) - jv(n + 1, k1 * a)) / 2.0
    dh0 = (hankel1(n - 1, k0 * a) - hankel1(n + 1, k0 * a)) / 2.0
    num = k1 * j0 * dj1 - k0 * dj0 * j1
    den = k0 * dh0 * j1 - k1 * h0 * dj1
    return num / den


def forward_scattered(
    freq_ghz: int,
    view_idx: int,
    center: tuple[float, float],
    eps: complex,
    a: float = A_RAD,
) -> np.ndarray:
    """Mie scattered field at all 72 receiver positions of one view (unit amp).

    Includes the view-dependent propagation phase exp(+i k0 |r_s - c|).
    """
    _, _, src, rec = cached_geometry()
    k0 = 2.0 * np.pi * freq_ghz * 1e9 / C0
    c = np.asarray(center, dtype=float)
    rs = src[view_idx - 1]
    u = c - rs
    dist_sc = np.linalg.norm(u)
    u = u / dist_sc
    nmax = int(np.ceil(k0 * a)) + 15
    w = jacobi_anger(k0, u, nmax, a)
    cn = mie_coefficients(k0, eps, a, nmax)
    an = w * cn * np.exp(1j * k0 * dist_sc)
    d = rec - c
    rho = np.linalg.norm(d, axis=1)
    phi = np.arctan2(d[:, 1], d[:, 0])
    n = np.arange(-nmax, nmax + 1)
    H = hankel1(n[:, None], k0 * rho[None, :])
    phase = np.exp(1j * n[:, None] * phi[None, :])
    out = np.sum(an[:, None] * H * phase, axis=0)
    return out


def forward_view_batch(
    freq_ghz: int,
    views,
    center: tuple[float, float],
    eps: complex,
    a: float = A_RAD,
) -> dict:
    """Mie scattered field per view; Mie coefficients computed once per freq."""
    _, _, src, rec = cached_geometry()
    k0 = 2.0 * np.pi * freq_ghz * 1e9 / C0
    c = np.asarray(center, dtype=float)
    nmax = int(np.ceil(k0 * a)) + 15
    cn = mie_coefficients(k0, eps, a, nmax)
    n = np.arange(-nmax, nmax + 1)
    out = {}
    d = rec - c
    rho = np.linalg.norm(d, axis=1)
    phi = np.arctan2(d[:, 1], d[:, 0])
    phase_r = np.exp(1j * n[:, None] * phi[None, :])
    H = hankel1(n[:, None], k0 * rho[None, :])
    for v in views:
        rs = src[int(v) - 1]
        u = c - rs
        dist_sc = np.linalg.norm(u)
        u = u / dist_sc
        w = jacobi_anger(k0, u, nmax, a)
        an = w * cn * np.exp(1j * k0 * dist_sc)
        out[int(v)] = np.sum(an[:, None] * H * phase_r, axis=0)
    return out


def self_checks() -> dict:
    """Verify the analytic core; raise on failure."""
    res: dict = {}
    _, _, src, rec = geometry()
    k0 = 2.0 * np.pi * 8e9 / C0
    c = np.array([0.03, 0.0])
    rs = src[0]
    u = c - rs
    u = u / np.linalg.norm(u)
    nmax = int(np.ceil(k0 * A_RAD)) + 15
    n = np.arange(-nmax, nmax + 1)

    # 1. Jacobi-Anger reconstruction on the cylinder boundary
    w = jacobi_anger(k0, u, nmax, A_RAD)
    phi_t = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False)
    rho_pts = np.stack([A_RAD * np.cos(phi_t), A_RAD * np.sin(phi_t)], axis=1) + c
    exact = np.exp(-1j * k0 * np.einsum("i,ni->n", u, rho_pts - c))
    recn = np.sum(
        w[:, None]
        * jv(n, k0 * A_RAD)[:, None]
        * np.exp(1j * n[:, None] * phi_t[None, :]),
        axis=0,
    )
    res["jacobi_anger_max_err"] = float(np.max(np.abs(recn - exact)))
    assert res["jacobi_anger_max_err"] < 1e-10, "Jacobi-Anger expansion failed"

    # 2. Boundary continuity for eps = 3.0 - 0.1j at 8 GHz
    eps = 3.0 - 0.1j
    k1 = k0 * np.sqrt(eps)
    cn = mie_coefficients(k0, eps, A_RAD, nmax)
    an = w * cn
    j0a, j1a = jv(n, k0 * A_RAD), jv(n, k1 * A_RAD)
    h0a = hankel1(n, k0 * A_RAD)
    dj0a = (jv(n - 1, k0 * A_RAD) - jv(n + 1, k0 * A_RAD)) / 2.0
    dj1a = (jv(n - 1, k1 * A_RAD) - jv(n + 1, k1 * A_RAD)) / 2.0
    dh0a = (hankel1(n - 1, k0 * A_RAD) - hankel1(n + 1, k0 * A_RAD)) / 2.0
    bn = (w * j0a + an * h0a) / j1a
    d_out = k0 * (w * dj0a + an * dh0a)
    d_in = k1 * bn * dj1a
    res["bc_derivative_max_rel_err"] = float(
        np.max(np.abs(d_out - d_in)) / np.max(np.abs(d_out))
    )
    assert res["bc_derivative_max_rel_err"] < 1e-10, "boundary continuity failed"

    # 3. Poynting flux conservation through rho = 0.2 m around the cylinder
    def flux(eps_c: complex) -> float:
        k1c = k0 * np.sqrt(eps_c)
        cnc = mie_coefficients(k0, eps_c, A_RAD, nmax)
        anc = w * cnc
        radius = 0.2
        z0 = k0 * radius
        ang = np.linspace(0.0, 2.0 * np.pi, 8001, endpoint=False)
        e_phi = np.exp(1j * n[:, None] * ang[None, :])
        j_n = jv(n, z0)
        h_n = hankel1(n, z0)
        dj_n = (jv(n - 1, z0) - jv(n + 1, z0)) / 2.0
        dh_n = (hankel1(n - 1, z0) - hankel1(n + 1, z0)) / 2.0
        et = np.sum(
            (w[:, None] * j_n[:, None] + anc[:, None] * h_n[:, None]) * e_phi, axis=0
        )
        det = np.sum(
            k0
            * (w[:, None] * dj_n[:, None] + anc[:, None] * dh_n[:, None])
            * e_phi,
            axis=0,
        )
        integrand = np.real(
            et * np.conj(det) / (1j * (2.0 * np.pi * 8e9) * MU0)
        )
        return float(np.trapezoid(integrand * radius, ang))

    f_lossless = flux(3.0 + 0j)
    f_lossy = flux(3.0 - 0.2j)
    res["lossless_net_flux"] = f_lossless
    res["lossy_net_flux"] = f_lossy
    assert f_lossy < 0, "lossy flux sign wrong (must absorb power)"
    res["flux_lossless_over_lossy"] = float(abs(f_lossless) / abs(f_lossy))
    assert (
        res["flux_lossless_over_lossy"] < 0.01
    ), "lossless flux not conserved within quadrature tolerance"

    # 3b. Exact conservation certificate: optical theorem for a lossless cylinder
    cn0 = mie_coefficients(k0, 3.0 + 0j, A_RAD, nmax)
    sig_scat = (4.0 / k0) * np.sum(np.abs(cn0) ** 2)
    sig_ext = -(4.0 / k0) * np.sum(np.real(cn0))
    res["optical_theorem_ratio"] = float(sig_scat / sig_ext)
    assert abs(res["optical_theorem_ratio"] - 1.0) < 1e-10, "optical theorem failed"

    # 4. PEC limit (approached as O(1/sqrt(eps)); check at eps = 1e8)
    c_pec = mie_coefficients(k0, 1e8 + 0j, A_RAD, nmax)
    target = -jv(n, k0 * A_RAD) / hankel1(n, k0 * A_RAD)
    res["pec_limit_max_err"] = float(np.max(np.abs(c_pec - target)))
    assert res["pec_limit_max_err"] < 1e-3, "PEC limit failed"

    # 5. Series truncation tail
    k0_1 = 2.0 * np.pi * 1e9 / C0
    nm1 = int(np.ceil(k0_1 * A_RAD)) + 15
    c_1 = mie_coefficients(k0_1, 3.0 + 0j, A_RAD, nm1)
    res["series_tail_rel_last"] = float(np.abs(c_1[-1]) / np.abs(c_1).max())
    assert res["series_tail_rel_last"] < 1e-10, "series truncation insufficient"
    return res


def rec_idx_for_view(view: int) -> np.ndarray:
    """0-based indices (1..72) of the 49 receivers recorded for a view."""
    start = 13 + 2 * (view - 1)
    idx = (start + np.arange(49)) % 72
    idx[idx == 0] = 72
    return idx - 1


def build_index(view, rec, freq, scattered):
    """Map (view, freq) -> row indices and the scattered (total - incident)."""
    key_to_row = {}
    for i in range(view.size):
        key = (int(view[i]), int(freq[i]))
        key_to_row.setdefault(key, []).append(i)
    return key_to_row, scattered


def gains_and_residual(params, eps, views, freqs, s_meas, key_to_row):
    """Profile per-frequency complex gains on `views`; return total residual."""
    gain, per_f = {}, {}
    tot_num = tot_den = 0.0
    for fq in freqs:
        batch = forward_view_batch(int(fq), views, params[:2], eps)
        m_list, s_list = [], []
        for v in views:
            rows = key_to_row.get((int(v), int(fq)))
            if rows is None:
                continue
            m_list.append(batch[int(v)][rec_idx_for_view(v)])
            s_list.append(s_meas[rows])
        Mv = np.concatenate(m_list)
        sv = np.concatenate(s_list)
        g = np.sum(np.conj(Mv) * sv) / np.sum(np.abs(Mv) ** 2)
        gain[int(fq)] = complex(g)
        num = np.sum(np.abs(g * Mv - sv) ** 2)
        den = np.sum(np.abs(sv) ** 2)
        per_f[int(fq)] = float(num / den)
        tot_num += num
        tot_den += den
    return float(tot_num / tot_den), {"gain": gain, "per_frequency": per_f}


def run_fits() -> dict:
    """Nominal-geometry vs joint low-dimensional fit on a held-out split."""
    view, rec, freq, total, incident = load_data()
    scattered = total - incident
    key_to_row, s_meas = build_index(view, rec, freq, scattered)
    train_views = np.arange(2, 37, 2)
    test_views = np.arange(1, 37, 2)
    freqs = np.arange(1, 9)

    nominal = (0.030, 0.0)
    eps_nom = 3.0 + 0j
    res_nom_tr, meta_nom_tr = gains_and_residual(
        (nominal[0], nominal[1], 3.0), eps_nom, train_views, freqs, s_meas, key_to_row
    )
    res_nom_te, meta_nom_te = gains_and_residual(
        (nominal[0], nominal[1], 3.0), eps_nom, test_views, freqs, s_meas, key_to_row
    )

    def objective(x):
        cx, cy, er = x
        return gains_and_residual(
            (cx, cy, er), er + 0j, train_views, freqs, s_meas, key_to_row
        )[0]

    bounds = [(-0.06, 0.06), (-0.06, 0.06), (1.2, 6.0)]
    starts = [
        (0.030, 0.0, 3.0),
        (0.0, 0.030, 3.0),
        (-0.030, 0.0, 3.0),
        (0.0, -0.030, 3.0),
        (0.0, 0.0, 3.0),
        (0.021, 0.021, 3.0),
        (-0.021, -0.021, 3.0),
        (0.030, 0.0, 2.4),
    ]
    best = None
    runs = []
    for st in starts:
        r = minimize(
            objective,
            np.asarray(st),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 150, "ftol": 1e-12},
        )
        runs.append(
            {
                "start": list(st),
                "x": [float(v) for v in r.x],
                "f": float(r.fun),
                "niter": int(r.nit),
                "success": bool(r.success),
            }
        )
        if best is None or r.fun < best.fun:
            best = r

    bx = best.x
    res_joint_tr, meta_joint_tr = gains_and_residual(
        (bx[0], bx[1], bx[2]), bx[2] + 0j, train_views, freqs, s_meas, key_to_row
    )
    res_joint_te, meta_joint_te = gains_and_residual(
        (bx[0], bx[1], bx[2]), bx[2] + 0j, test_views, freqs, s_meas, key_to_row
    )

    def objective4(x):
        cx, cy, er, ei = x
        return gains_and_residual(
            (cx, cy, er), er - 1j * ei, train_views, freqs, s_meas, key_to_row
        )[0]

    r4 = minimize(
        objective4,
        np.asarray([bx[0], bx[1], bx[2], 0.0]),
        method="L-BFGS-B",
        bounds=[(-0.06, 0.06), (-0.06, 0.06), (1.2, 6.0), (0.0, 0.6)],
        options={"maxiter": 150, "ftol": 1e-12},
    )
    x4 = r4.x
    res_joint4_tr, meta_joint4_tr = gains_and_residual(
        (x4[0], x4[1], x4[2]),
        x4[2] - 1j * x4[3],
        train_views,
        freqs,
        s_meas,
        key_to_row,
    )
    res_joint4_te, meta_joint4_te = gains_and_residual(
        (x4[0], x4[1], x4[2]),
        x4[2] - 1j * x4[3],
        test_views,
        freqs,
        s_meas,
        key_to_row,
    )

    return {
        "split": {
            "train_views": [int(v) for v in train_views],
            "test_views": [int(v) for v in test_views],
            "frequencies_ghz": [int(f) for f in freqs],
        },
        "nominal": {
            "params": {"cx_m": nominal[0], "cy_m": nominal[1], "eps_r": eps_nom.real},
            "anchor_note": "centre direction unpublished; +x used as gauge anchor",
            "train_residual": float(res_nom_tr),
            "test_residual": float(res_nom_te),
            "train_per_freq": meta_nom_tr["per_frequency"],
            "test_per_freq": meta_nom_te["per_frequency"],
        },
        "joint_3param": {
            "params": {"cx_m": float(bx[0]), "cy_m": float(bx[1]), "eps_r": float(bx[2])},
            "center_radius_mm": float(np.hypot(bx[0], bx[1]) * 1e3),
            "center_angle_deg": float(np.rad2deg(np.arctan2(bx[1], bx[0]))),
            "train_residual": float(res_joint_tr),
            "test_residual": float(res_joint_te),
            "train_per_freq": meta_joint_tr["per_frequency"],
            "test_per_freq": meta_joint_te["per_frequency"],
            "optimizer_runs": runs,
        },
        "joint_4param": {
            "params": {
                "cx_m": float(x4[0]),
                "cy_m": float(x4[1]),
                "eps_r": float(x4[2]),
                "loss_tan": float(x4[3] / max(x4[2], 1e-9)),
            },
            "train_residual": float(res_joint4_tr),
            "test_residual": float(res_joint4_te),
            "train_per_freq": meta_joint4_tr["per_frequency"],
            "test_per_freq": meta_joint4_te["per_frequency"],
        },
    }


def robustness_pilot() -> dict:
    """Inject a known perturbed geometry into synthetic data and refit it."""
    view, rec, freq, total, incident = load_data()
    rng = np.random.default_rng(7)
    c_inj = np.array([0.024, -0.008])
    eps_inj = 2.7 + 0j
    train_views = np.arange(2, 37, 2)
    test_views = np.arange(1, 37, 2)
    freqs = np.arange(1, 9)

    key_syn, key_test = {}, {}
    for fq in freqs:
        batch = forward_view_batch(int(fq), range(1, 37), c_inj, eps_inj)
        syn = [batch[v][rec_idx_for_view(v)] for v in range(1, 37)]
        syn = np.concatenate(syn)
        sigma = 0.05 * np.sqrt(np.mean(np.abs(syn) ** 2))
        syn = syn + sigma * (
            rng.standard_normal(syn.size) + 1j * rng.standard_normal(syn.size)
        )
        for j, v in enumerate(range(1, 37)):
            target = key_syn if v in train_views else key_test
            target[(int(v), int(fq))] = syn[j * 49 : (j + 1) * 49]

    def obj(x):
        cx, cy, er = x
        tot_num = tot_den = 0.0
        for fq in freqs:
            batch = forward_view_batch(int(fq), train_views, (cx, cy), er + 0j)
            m_list, s_list = [], []
            for v in train_views:
                m_list.append(batch[int(v)][rec_idx_for_view(v)])
                s_list.append(key_syn[(int(v), int(fq))])
            Mv, sv = np.concatenate(m_list), np.concatenate(s_list)
            g = np.sum(np.conj(Mv) * sv) / np.sum(np.abs(Mv) ** 2)
            tot_num += np.sum(np.abs(g * Mv - sv) ** 2)
            tot_den += np.sum(np.abs(sv) ** 2)
        return tot_num / tot_den

    bounds = [(-0.06, 0.06), (-0.06, 0.06), (1.2, 6.0)]
    best = None
    for st in [(0.0, 0.0, 3.0), (0.03, 0.0, 3.0), (-0.03, 0.0, 2.5), (0.0, 0.03, 3.5)]:
        r = minimize(
            obj,
            np.asarray(st),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 120, "ftol": 1e-13},
        )
        if best is None or r.fun < best.fun:
            best = r

    bx = best.x

    def eval_res(views, key):
        tot_num = tot_den = 0.0
        for fq in freqs:
            batch = forward_view_batch(int(fq), views, (bx[0], bx[1]), bx[2] + 0j)
            m_list, s_list = [], []
            for v in views:
                m_list.append(batch[int(v)][rec_idx_for_view(v)])
                s_list.append(key[(int(v), int(fq))])
            Mv, sv = np.concatenate(m_list), np.concatenate(s_list)
            g = np.sum(np.conj(Mv) * sv) / np.sum(np.abs(Mv) ** 2)
            tot_num += np.sum(np.abs(g * Mv - sv) ** 2)
            tot_den += np.sum(np.abs(sv) ** 2)
        return tot_num / tot_den

    return {
        "injected": {
            "cx_m": float(c_inj[0]),
            "cy_m": float(c_inj[1]),
            "eps_r": float(eps_inj.real),
        },
        "recovered": {
            "cx_m": float(bx[0]),
            "cy_m": float(bx[1]),
            "eps_r": float(bx[2]),
        },
        "center_error_mm": float(np.hypot(bx[0] - c_inj[0], bx[1] - c_inj[1]) * 1e3),
        "eps_error": float(bx[2] - eps_inj.real),
        "train_residual": float(eval_res(train_views, key_syn)),
        "test_residual": float(eval_res(test_views, key_test)),
        "noise_level": "5% RMS complex Gaussian per frequency",
        "caveat": "injected perturbation is synthetic, not natural unknown-array truth",
    }


def low_frequency_diagnostic() -> dict:
    """Stability of the joint fit as the frequency band is extended."""
    view, rec, freq, total, incident = load_data()
    scattered = total - incident
    key, s = build_index(view, rec, freq, scattered)
    train = np.arange(2, 37, 2)
    test = np.arange(1, 37, 2)
    out = {}
    for fmax in [2, 3, 4]:
        freqs = np.arange(1, fmax + 1)

        def obj(x):
            return gains_and_residual(
                (x[0], x[1], x[2]), x[2] + 0j, train, freqs, s, key
            )[0]

        bounds = [(-0.06, 0.06), (-0.06, 0.06), (1.2, 6.0)]
        best = None
        for st in [
            (0.03, 0.0, 3.0),
            (0.0, 0.03, 3.0),
            (-0.03, 0.0, 3.0),
            (0.0, -0.03, 3.0),
            (0.0, 0.0, 3.0),
        ]:
            r = minimize(
                obj,
                np.asarray(st),
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 150, "ftol": 1e-12},
            )
            if best is None or r.fun < best.fun:
                best = r
        bx = best.x
        rtr, _ = gains_and_residual(
            (bx[0], bx[1], bx[2]), bx[2] + 0j, train, freqs, s, key
        )
        rte, _ = gains_and_residual(
            (bx[0], bx[1], bx[2]), bx[2] + 0j, test, freqs, s, key
        )
        out[f"up_to_{fmax}GHz"] = {
            "cx_m": float(bx[0]),
            "cy_m": float(bx[1]),
            "center_radius_mm": float(np.hypot(bx[0], bx[1]) * 1e3),
            "center_angle_deg": float(np.rad2deg(np.arctan2(bx[1], bx[0]))),
            "eps_r": float(bx[2]),
            "train_residual": float(rtr),
            "test_residual": float(rte),
        }
    return out


def main() -> int:
    checks = self_checks()
    fits = run_fits()
    pilot = robustness_pilot()
    diag = low_frequency_diagnostic()
    results = {
        "self_checks": checks,
        "fits": fits,
        "robustness_pilot": pilot,
        "low_frequency_diagnostic": diag,
    }
    out = ROOT / "results_model.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, sort_keys=True, default=complex)
    print(json.dumps(results, indent=2, default=str))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
