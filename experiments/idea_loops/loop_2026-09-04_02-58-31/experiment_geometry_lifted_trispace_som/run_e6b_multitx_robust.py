"""E6b -- robustness sweep + directional phased-dipole multi-transmitter model.

Strengthens E6 (oriented multi-transmitter theta identifiability) by:
  * monopole multi-transmitter sweeps: more SNR levels and seeds for L=2/L=3
    tx_only and an L=2 co-moving (p_rx == p_tx) identifiability check;
  * a new phased-dipole transmitter model (each "source" is a tangent-oriented
    pair of phase-weighted monopoles) with finite-difference validation,
    B_t gauge/rank analysis (L=1 vs L=2), and nonlinear tx_only runs.

Monopole nonlinear solves reuse E6's run_case/make_observation exactly
(method strings "wrongpose" and "direct"; the prompt's "direct_tx" is the
E6 "direct" method at pose_mode="tx_only" -- mapping recorded in
parameters.method_alias and per-run method_alias).  The dipole nonlinear
runner mirrors E6's run_case 1:1 but calls the dipole forward/Jacobian.

Outputs: results_e6b_multitx_robust.json, plot_e6b_multitx_robust.png
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares

from geom_som_core import (
    Params,
    pixel_grid,
    green_domain_matrix,
    incident_field,
    state_operator,
    solve_current,
    total_field,
    receiver_positions,
    receiver_directions,
    data_matrix,
    receiver_data_derivative,
)
from run_e5 import basis_matrix, _transmitter_incident_derivative, _realified_residual
from run_e6_multitx import (
    make_t_locals,
    _mode_split,
    run_case,
    make_observation,
    multi_jacobians,
    LS_KWARGS,
    DIRS,
    _rot,
    _real_ranks,
    _theta_gauge_residual,
)


HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Fixed scene / experiment configuration (identical base scene to E6)
# ---------------------------------------------------------------------------
P = Params(N=16, k=12.0, M=8, R_r=1.6, R_t=2.0, phi0=0.0, s=0.12)
ALPHA_TRUE = np.array([1.5, 2.0, 0.0])
ALPHA_INIT = np.array([1.0, 1.0, 0.0])
P_TRUE = np.array([0.08, -0.06, 0.05])
P_RX_FIXED = np.zeros(3)
P_INIT_ZERO = np.zeros(3)
P_INITS = {
    "zero": np.zeros(3),
    "pert_pos": np.array([0.18, 0.04, 0.15]),
}

MONO_ARRAYS = {
    "mono_L2": {
        "angles": [0.7, 2.0],
        "pose_mode": "tx_only",
        "snr_dB": [15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
        "model": "monopole",
    },
    "mono_L3": {
        "angles": [0.7, 2.0, 1.3],
        "pose_mode": "tx_only",
        "snr_dB": [20.0, 30.0, 40.0],
        "model": "monopole",
    },
    "mono_L2_comoving": {
        "angles": [0.7, 2.0],
        "pose_mode": "co_moving",
        "snr_dB": [20.0, 30.0, 40.0],
        "model": "monopole",
    },
    "dipole_L2": {
        "angles": [0.7, 2.0],
        "pose_mode": "tx_only",
        "snr_dB": [20.0, 30.0, 40.0],
        "model": "dipole",
    },
}

DIPOLE_D = 0.08
DIPOLE_PHASE = 0.5 * np.pi
NOISE_SEEDS = list(range(8))
EPS_LIST = (1e-2, 1e-3, 1e-4, 1e-5)

METHOD_ALIAS = {
    ("tx_only", "wrongpose"): "wrongpose",
    ("tx_only", "direct"): "direct_tx",
    ("co_moving", "wrongpose"): "wrongpose",
    ("co_moving", "direct"): "direct_comoving",
}


def method_alias_of(pose_mode: str, method: str) -> str:
    return METHOD_ALIAS[(pose_mode, method)]


# ---------------------------------------------------------------------------
# Directional phased-dipole transmitter model
# ---------------------------------------------------------------------------
def _dipole_axis(orientation: float) -> np.ndarray:
    return np.array([np.cos(orientation), np.sin(orientation)])


def dipole_incident_field(
    xs: np.ndarray,
    t_center: np.ndarray,
    orientation: float,
    d: float,
    phase: float,
    k: float,
) -> np.ndarray:
    """u_inc of a phased dipole centred at t_center, axis u(orientation)."""
    u = _dipole_axis(orientation)
    t_plus = t_center + d * u
    t_minus = t_center - d * u
    w_plus = np.exp(0.5j * phase)
    w_minus = np.exp(-0.5j * phase)
    return w_plus * incident_field(xs, t_plus, k) + w_minus * incident_field(
        xs, t_minus, k
    )


def dipole_incident_derivative(
    xs: np.ndarray,
    t_center: np.ndarray,
    orientation: float,
    d: float,
    phase: float,
    k: float,
    w_t: np.ndarray,
) -> np.ndarray:
    """du_inc/dt when the whole dipole translates rigidly at velocity w_t."""
    u = _dipole_axis(orientation)
    t_plus = t_center + d * u
    t_minus = t_center - d * u
    w_plus = np.exp(0.5j * phase)
    w_minus = np.exp(-0.5j * phase)
    return w_plus * _transmitter_incident_derivative(
        xs, t_plus, w_t, k
    ) + w_minus * _transmitter_incident_derivative(xs, t_minus, w_t, k)


def _dipole_source_components(
    tl: np.ndarray,
    p_tx: np.ndarray,
    R_tx: np.ndarray,
    k: float,
    d: float = DIPOLE_D,
    phase: float = DIPOLE_PHASE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, complex, complex]:
    """Body-frame geometry of one tangent phased dipole at local centre tl.

    The two monopole offsets are defined in the transmitter body frame
    (q = tl +/- d*u(local_tangent)) and then rotated into the world frame,
    so that the theta derivative per monopole is R_tx @ R90(q).
    Returns q_plus, q_minus, t_plus, t_minus, w_plus, w_minus.
    """
    local_angle = float(np.arctan2(tl[1], tl[0]))
    tangent_local = local_angle + 0.5 * np.pi
    u = _dipole_axis(tangent_local)
    q_plus = tl + d * u
    q_minus = tl - d * u
    t_plus = R_tx @ q_plus + p_tx[:2]
    t_minus = R_tx @ q_minus + p_tx[:2]
    w_plus = np.exp(0.5j * phase)
    w_minus = np.exp(-0.5j * phase)
    return q_plus, q_minus, t_plus, t_minus, w_plus, w_minus


def _rot90(v: np.ndarray) -> np.ndarray:
    return np.array([-v[1], v[0]])


def dipole_multi_forward(
    alpha: np.ndarray,
    p_tx: np.ndarray,
    p_rx: np.ndarray,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    t_locals: np.ndarray,
    P: Params,
    want_sources: bool = False,
):
    """Stacked L-dipole data d = [d_0; ...] (complex length L*M)."""
    chi = Phi @ alpha
    y, _ = receiver_positions(p_rx, P)
    G_s = data_matrix(y, xs, h, k)
    R_tx = _rot(p_tx[2])
    ds = []
    for tl in t_locals:
        _, _, t_plus, t_minus, w_plus, w_minus = _dipole_source_components(
            tl, p_tx, R_tx, k
        )
        u_inc = w_plus * incident_field(xs, t_plus, k) + w_minus * incident_field(
            xs, t_minus, k
        )
        J = solve_current(chi, u_inc, G_D)
        ds.append(G_s @ J)
    d_all = np.concatenate(ds)
    return (d_all, ds) if want_sources else d_all


def dipole_multi_jacobians(
    alpha: np.ndarray,
    p_tx: np.ndarray,
    p_rx: np.ndarray,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    t_locals: np.ndarray,
    P: Params,
    include_receiver: bool = True,
) -> dict:
    """Analytic Jacobians for the stacked phased-dipole transmitter model.

    Transmitter columns: tx/ty move both monopoles rigidly; theta rotates the
    body frame, so each monopole's velocity is R_tx @ R90(q_mono_local).
    Receiver-side columns are identical to E6's multi_jacobians.
    """
    chi = Phi @ alpha
    A = state_operator(chi, G_D)
    lu, piv = lu_factor(A)

    y, y_local = receiver_positions(p_rx, P)
    G_s = data_matrix(y, xs, h, k)
    R_tx = _rot(p_tx[2])
    R_rx = _rot(p_rx[2])

    d_l = []
    J_alpha_l = []
    B_t_l = []
    B_r_l = []
    for tl in t_locals:
        q_plus, q_minus, t_plus, t_minus, w_plus, w_minus = (
            _dipole_source_components(tl, p_tx, R_tx, k)
        )
        u_inc = w_plus * incident_field(xs, t_plus, k) + w_minus * incident_field(
            xs, t_minus, k
        )
        J = lu_solve((lu, piv), chi * u_inc)
        u = total_field(u_inc, J, G_D)
        d_l.append(G_s @ J)

        J_chi = lu_solve((lu, piv), np.diag(u))  # Np x Np
        J_alpha_l.append(G_s @ J_chi @ Phi)

        B_t = np.empty((P.M, 3), dtype=complex)
        for c, direction in enumerate(DIRS):
            # E6 column labels are ("rx","ry","theta"); for a transmitter they
            # mean global x/y translation and body rotation respectively.
            if direction == "theta":  # per-monopole rigid-rotation velocities
                v_plus = R_tx @ _rot90(q_plus)
                v_minus = R_tx @ _rot90(q_minus)
            elif direction == "rx":
                v_plus = v_minus = np.array([1.0, 0.0])
            elif direction == "ry":
                v_plus = v_minus = np.array([0.0, 1.0])
            else:
                raise ValueError(f"unknown pose column {direction!r}")
            du_inc = w_plus * _transmitter_incident_derivative(
                xs, t_plus, v_plus, k
            ) + w_minus * _transmitter_incident_derivative(
                xs, t_minus, v_minus, k
            )
            J_p = lu_solve((lu, piv), chi * du_inc)
            B_t[:, c] = G_s @ J_p
        B_t_l.append(B_t)

        B_r = np.zeros((P.M, 3), dtype=complex)
        if include_receiver:
            for c, direction in enumerate(DIRS):
                v = receiver_directions(direction, y_local)
                if direction == "theta":
                    v = v @ R_rx.T
                B_r[:, c] = receiver_data_derivative(y, v, xs, h, k) @ J
        B_r_l.append(B_r)

    d = np.concatenate(d_l)
    J_alpha = np.vstack(J_alpha_l)  # L*M x K
    B_t = np.vstack(B_t_l)  # L*M x 3
    B_r = np.vstack(B_r_l)
    B_total = B_t + B_r
    return {
        "d": d,
        "J_alpha": J_alpha,
        "B_t": B_t,
        "B_r": B_r,
        "B_total": B_total,
        "J_alpha_real": np.vstack([J_alpha.real, J_alpha.imag]),
        "B_t_real": np.vstack([B_t.real, B_t.imag]),
        "B_r_real": np.vstack([B_r.real, B_r.imag]),
        "B_real": np.vstack([B_total.real, B_total.imag]),
    }


# ---------------------------------------------------------------------------
# Finite-difference verification of dipole_multi_jacobians
# ---------------------------------------------------------------------------
def _per_col_rel(ana: np.ndarray, fd_cols: list[np.ndarray]) -> list[float]:
    """Per-column relative errors norm(ana-fd)/norm(fd)."""
    out = []
    for c, fd in enumerate(fd_cols):
        ana_col = np.concatenate([ana[:, c].real, ana[:, c].imag])
        out.append(float(np.linalg.norm(ana_col - fd) / np.linalg.norm(fd)))
    return out


def dipole_fd_check(
    alpha: np.ndarray,
    p0: np.ndarray,
    t_locals: np.ndarray,
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    eps_list,
    pose_mode: str,
    p_rx_fixed: np.ndarray,
) -> list[dict]:
    """Centred FD vs analytic dipole B columns.

    For tx_only the receiver is fixed, so the model FD equals B_t (= B_total).
    For co_moving we check two FD conventions:
      * receiver_fixed : p_rx stays at p0's receiver -> analytic B_t
      * co_moving      : p_rx == p_tx moves -> analytic B_total = B_t + B_r
    """
    p0 = np.asarray(p0, dtype=float)
    p_tx0, p_rx0 = _mode_split(p0, pose_mode, p_rx_fixed)
    jb = dipole_multi_jacobians(
        alpha,
        p_tx0,
        p_rx0,
        xs,
        h,
        k,
        G_D,
        Phi,
        t_locals,
        P,
        include_receiver=(pose_mode == "co_moving"),
    )

    rows = []
    for eps in eps_list:
        checks = [("B_t", jb["B_t"], False)]
        if pose_mode == "co_moving":
            # receiver follows the transmitter: FD over the full pose change
            # validates B_total = B_t + B_r.
            checks.append(("B_total", jb["B_total"], True))
        for matrix_name, ana, fd_receiver_moves in checks:
            per_col = []
            fd_cols = []
            for c in range(3):
                ep = np.zeros(3)
                ep[c] = eps
                if fd_receiver_moves:
                    p_tx_plus, p_rx_plus = _mode_split(
                        p0 + ep, pose_mode, p_rx_fixed
                    )
                    p_tx_minus, p_rx_minus = _mode_split(
                        p0 - ep, pose_mode, p_rx_fixed
                    )
                else:
                    p_tx_plus, p_tx_minus = p0 + ep, p0 - ep
                    p_rx_plus = p_rx_minus = p_rx0
                d_plus = dipole_multi_forward(
                    alpha,
                    p_tx_plus,
                    p_rx_plus,
                    xs,
                    h,
                    k,
                    G_D,
                    Phi,
                    t_locals,
                    P,
                )
                d_minus = dipole_multi_forward(
                    alpha,
                    p_tx_minus,
                    p_rx_minus,
                    xs,
                    h,
                    k,
                    G_D,
                    Phi,
                    t_locals,
                    P,
                )
                fd_cols.append(
                    _realified_residual(d_plus, d_minus) / (2.0 * eps)
                )
            per_col = _per_col_rel(ana, fd_cols)
            rows.append(
                {
                    "pose_mode": pose_mode,
                    "matrix": matrix_name,
                    "fd_convention": (
                        "co_moving"
                        if fd_receiver_moves
                        else "receiver_fixed"
                    ),
                    "eps": float(eps),
                    "per_col_rel_err": per_col,
                    "max_rel_err": float(max(per_col)),
                }
            )
            print(
                f"[E6b FD dipole] mode={pose_mode:>9s} matrix={matrix_name:>8s} "
                f"fd={rows[-1]['fd_convention']:>14s} eps={eps:.0e} "
                f"max_rel_err={max(per_col):.3e}"
            )
    return rows


# ---------------------------------------------------------------------------
# Dipole rank/gauge analysis (mirrors E6's rank_gauge_analysis helpers)
# ---------------------------------------------------------------------------
def dipole_rank_gauge_analysis(
    alpha_init: np.ndarray,
    t_locals_by_L: dict,
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
) -> tuple[list[dict], list[dict]]:
    out = []
    p0 = np.zeros(3)
    for L in sorted(t_locals_by_L):
        tl = t_locals_by_L[L]
        jb = dipole_multi_jacobians(
            alpha_init,
            p0,
            p0,
            xs,
            h,
            k,
            G_D,
            Phi,
            tl,
            P,
            include_receiver=False,
        )
        B_t = jb["B_t"]
        rank_c, real_sv, normed_real_sv = _real_ranks(B_t)
        c_sv = np.linalg.svd(B_t, compute_uv=False)
        res = _theta_gauge_residual(B_t)
        res_rel = float(res / np.linalg.norm(B_t[:, 2]))
        out.append(
            {
                "L": int(L),
                "model": "dipole",
                "d": float(DIPOLE_D),
                "phase": float(DIPOLE_PHASE),
                "t_local_center_angles": [
                    float(np.arctan2(t[1], t[0])) for t in tl
                ],
                "complex_rank": rank_c,
                "complex_sv": [float(v) for v in c_sv],
                "realified_sv": [float(v) for v in real_sv],
                "realified_colnormed_sv": [float(v) for v in normed_real_sv],
                "theta_residual_abs": res,
                "theta_residual_rel": res_rel,
            }
        )
        print(
            f"[E6b gauge dipole] L={L} complex_rank={rank_c} complex_sv="
            + ",".join(f"{v:.3e}" for v in c_sv)
            + f" theta_res_abs={res:.3e} rel={res_rel:.3e}"
        )

    comov = []
    for L in sorted(t_locals_by_L):
        tl = t_locals_by_L[L]
        jb = dipole_multi_jacobians(
            alpha_init,
            p0,
            p0,
            xs,
            h,
            k,
            G_D,
            Phi,
            tl,
            P,
            include_receiver=True,
        )
        rank_c, real_sv, _ = _real_ranks(jb["B_total"])
        comov.append(
            {
                "L": int(L),
                "model": "dipole",
                "B_total_complex_rank": rank_c,
                "B_total_realified_sv": [float(v) for v in real_sv],
            }
        )
        print(
            f"[E6b gauge dipole] co-moving L={L} B_total complex_rank={rank_c} "
            "realified_sv=" + ",".join(f"{v:.3e}" for v in real_sv)
        )
    return out, comov


# ---------------------------------------------------------------------------
# Dipole observation + nonlinear runner (E6-compatible fields/conventions)
# ---------------------------------------------------------------------------
def dipole_make_observation(
    alpha_true: np.ndarray,
    p_true: np.ndarray,
    pose_mode: str,
    p_rx_fixed: np.ndarray,
    snr_db: float,
    seed: int,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    t_locals: np.ndarray,
    P: Params,
) -> tuple[np.ndarray, np.ndarray]:
    """Same noise convention as E6 make_observation, dipole forward."""
    p_tx, p_rx = _mode_split(p_true, pose_mode, p_rx_fixed)
    d_true = dipole_multi_forward(
        alpha_true, p_tx, p_rx, xs, h, k, G_D, Phi, t_locals, P
    )
    ratio = float(10.0 ** (-snr_db / 20.0))
    rng = np.random.default_rng(seed)
    raw = (
        rng.standard_normal(d_true.size)
        + 1j * rng.standard_normal(d_true.size)
    ) / np.sqrt(2.0)
    noise = raw * (ratio * np.linalg.norm(d_true) / np.linalg.norm(raw))
    return d_true + noise, d_true


def dipole_run_case(
    method: str,
    pose_mode: str,
    snr_db: float,
    seed: int,
    p_init_label: str,
    p_init: np.ndarray,
    alpha_init: np.ndarray,
    alpha_true: np.ndarray,
    p_true: np.ndarray,
    d_obs: np.ndarray,
    ref_norm: float,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    t_locals: np.ndarray,
    P: Params,
    p_rx_fixed: np.ndarray,
) -> dict:
    """One least-squares solve with the dipole model (E6 run_case fields)."""
    assert method in ("wrongpose", "direct")
    K = Phi.shape[1]

    def residual(theta: np.ndarray) -> np.ndarray:
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
        else:
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
        p_tx, p_rx = _mode_split(p, pose_mode, p_rx_fixed)
        return _realified_residual(
            dipole_multi_forward(
                alpha, p_tx, p_rx, xs, h, k, G_D, Phi, t_locals, P
            ),
            d_obs,
        )

    def jacobian(theta: np.ndarray) -> np.ndarray:
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
        else:
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
        p_tx, p_rx = _mode_split(p, pose_mode, p_rx_fixed)
        jb = dipole_multi_jacobians(
            alpha,
            p_tx,
            p_rx,
            xs,
            h,
            k,
            G_D,
            Phi,
            t_locals,
            P,
            include_receiver=(pose_mode == "co_moving"),
        )
        if method == "wrongpose":
            return jb["J_alpha_real"]
        return np.hstack([jb["J_alpha_real"], jb["B_real"]])

    if method == "wrongpose":
        x0 = np.asarray(alpha_init, dtype=float)
    else:
        x0 = np.concatenate(
            [np.asarray(alpha_init, float), np.asarray(p_init, float)]
        )

    result = least_squares(residual, x0, jac=jacobian, **LS_KWARGS)

    if method == "wrongpose":
        alpha_est = np.asarray(result.x[:K], dtype=float)
        p_est = np.zeros(3)
    else:
        alpha_est = np.asarray(result.x[:K], dtype=float)
        p_est = np.asarray(result.x[K:], dtype=float)

    b_final = residual(result.x)
    final_residual = float(np.linalg.norm(b_final) / ref_norm)
    chi_true = Phi @ alpha_true
    chi_est = Phi @ alpha_est
    map_error = float(
        np.linalg.norm(chi_est - chi_true) / np.linalg.norm(chi_true)
    )
    delta_p = p_est - p_true
    coord_err = np.abs(delta_p)
    success = bool(
        result.status > 0
        and (result.optimality < 1e-8 or result.status in (1, 2, 3, 4))
    )
    record = {
        "method": method,
        "pose_mode": pose_mode,
        "SNR_dB": float(snr_db),
        "noise_seed": int(seed),
        "p_init_label": p_init_label,
        "p_init": [float(v) for v in p_init],
        "alpha_est": [float(v) for v in alpha_est],
        "p_est": [float(v) for v in p_est],
        "pose_error": float(np.linalg.norm(delta_p)),
        "tx_err": float(coord_err[0]),
        "ty_err": float(coord_err[1]),
        "theta_err": float(coord_err[2]),
        "map_error": map_error,
        "final_residual": final_residual,
        "success": success,
        "nfev": int(result.nfev),
        "optimality": float(result.optimality),
        "n_vis": None,
        "hid_sv": None,
        "hidden_rank": None,
    }
    return record


def finalize_record(rec: dict, group: str) -> dict:
    """Add E6b group/alias fields and the requested per-run field names."""
    out = dict(rec)
    out["group"] = group
    out["seed"] = int(rec["noise_seed"])
    out["method_alias"] = method_alias_of(rec["pose_mode"], rec["method"])
    out["tx_error"] = float(rec["tx_err"])
    out["ty_error"] = float(rec["ty_err"])
    out["theta_error"] = float(rec["theta_err"])
    out["final_residual_ratio"] = float(rec["final_residual"])
    return out


def med_iqr(values: list[float]) -> tuple[float, float]:
    vals = np.asarray(values, dtype=float)
    return (
        float(np.median(vals)),
        float(np.percentile(vals, 75) - np.percentile(vals, 25)),
    )


def make_summary(records: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = defaultdict(list)
    order: list[tuple] = []
    for rec in records:
        key = (
            rec["group"],
            rec["pose_mode"],
            rec["method_alias"],
            rec["SNR_dB"],
            rec["p_init_label"],
        )
        if key not in groups:
            order.append(key)
        groups[key].append(rec)
    rows = []
    for key in order:
        sel = groups[key]

        def stat(field):
            return med_iqr([r[field] for r in sel])

        pose_med, pose_iqr = stat("pose_error")
        theta_med, theta_iqr = stat("theta_error")
        map_med, map_iqr = stat("map_error")
        rows.append(
            {
                "group": key[0],
                "pose_mode": key[1],
                "method_alias": key[2],
                "method": sel[0]["method"],
                "SNR_dB": key[3],
                "p_init_label": key[4],
                "n_seeds": len(sel),
                "n_success": int(sum(r["success"] for r in sel)),
                "median_pose_error": pose_med,
                "iqr_pose_error": pose_iqr,
                "median_theta_err": theta_med,
                "iqr_theta_err": theta_iqr,
                "median_map_error": map_med,
                "iqr_map_error": map_iqr,
                "median_final_residual": med_iqr(
                    [r["final_residual"] for r in sel]
                )[0],
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
GROUP_COLORS = {
    "mono_L2": "#1f77b4",
    "mono_L3": "#2ca02c",
    "mono_L2_comoving": "#ff7f0e",
    "dipole_L2": "#d62728",
}
GROUP_LABELS = {
    "mono_L2": "L=2 mono tx_only",
    "mono_L3": "L=3 mono tx_only",
    "mono_L2_comoving": "L=2 mono co_moving",
    "dipole_L2": "L=2 dipole tx_only",
}
SNR_COLS = [15.0, 20.0, 25.0, 30.0, 35.0, 40.0]


def _series(records: list[dict], group: str, alias: str, field: str):
    """Median/IQR series vs SNR for direct-type runs of one tx group."""
    rows = [
        r
        for r in records
        if r["group"] == group
        and r["method_alias"] == alias
        and r["p_init_label"] == "zero"
    ]
    by_snr: dict[float, list[float]] = defaultdict(list)
    for r in rows:
        by_snr[float(r["SNR_dB"])].append(float(r[field]))
    snrs = sorted(by_snr)
    med = [float(np.median(by_snr[s])) for s in snrs]
    lo = [float(np.percentile(by_snr[s], 25)) for s in snrs]
    hi = [float(np.percentile(by_snr[s], 75)) for s in snrs]
    return np.asarray(snrs), np.asarray(med), np.asarray(lo), np.asarray(hi)


def _plot_series(ax, records, groups, field, ylabel, title):
    for g in groups:
        snrs, med, lo, hi = _series(records, g, "direct_tx", field)
        c = GROUP_COLORS[g]
        (line,) = ax.plot(snrs, med, "-o", color=c, lw=1.8, ms=5,
                          label=GROUP_LABELS[g])
        lo_c = np.maximum(lo, med.min() * 1e-6)
        ax.fill_between(snrs, lo_c, hi, color=c, alpha=0.18, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)


def make_plot(
    records: list[dict],
    gauge_mono: list[dict],
    gauge_dipole: list[dict],
    png_path: Path,
) -> None:
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.24)

    # 1: median pose error vs SNR (direct tx-only groups)
    ax1 = fig.add_subplot(gs[0, 0])
    _plot_series(
        ax1,
        records,
        ("mono_L2", "mono_L3", "dipole_L2"),
        "pose_error",
        "median pose error (medians, IQR shaded)",
        "(1) tx_only direct: median pose error vs SNR",
    )

    # 2: median theta error vs SNR
    ax2 = fig.add_subplot(gs[0, 1])
    _plot_series(
        ax2,
        records,
        ("mono_L2", "mono_L3", "dipole_L2"),
        "theta_err",
        "median |theta_err| (medians, IQR shaded)",
        "(2) tx_only direct: median theta error vs SNR",
    )

    # 3: map error wrongpose vs direct_tx (L=2 mono) + direct dipole
    ax3 = fig.add_subplot(gs[1, 0])
    map_groups = [
        ("mono_L2", "wrongpose", "L=2 mono wrongpose"),
        ("mono_L2", "direct_tx", "L=2 mono direct_tx"),
        ("dipole_L2", "direct_tx", "L=2 dipole direct_tx"),
    ]
    colors3 = ["#d62728", "#1f77b4", "#9467bd"]
    for (g, alias, lab), c in zip(map_groups, colors3):
        snrs, med, lo, hi = _series(records, g, alias, "map_error")
        ax3.plot(snrs, med, "-o", color=c, lw=1.8, ms=5, label=lab)
        ax3.fill_between(
            snrs,
            np.maximum(lo, np.maximum(med.min(), 1e-10) * 1e-5),
            hi,
            color=c,
            alpha=0.18,
        )
    ax3.set_yscale("log")
    ax3.set_xlabel("SNR (dB)")
    ax3.set_ylabel("median ||chi_est-chi_true||/||chi_true||")
    ax3.set_title("(3) map error: L=2 mono wrongpose/direct + dipole direct",
                  fontsize=10)
    ax3.grid(alpha=0.3, which="both")
    ax3.legend(fontsize=7)

    # 4: tx_only vs co_moving, L=2 mono direct, pose+theta grouped bars
    ax4 = fig.add_subplot(gs[1, 1])
    snr4 = [20.0, 30.0, 40.0]
    metrics = [
        ("pose_error", "pose error", "#1f77b4"),
        ("theta_err", "theta error", "#2ca02c"),
    ]
    widths = 0.32
    for mi, (field, lab, c) in enumerate(metrics):
        base = mi * len(snr4) + np.arange(len(snr4))
        tx = np.array(
            [
                float(np.median([
                    r[field] for r in records
                    if r["group"] == "mono_L2"
                    and r["method_alias"] == "direct_tx"
                    and r["SNR_dB"] == s
                ]))
                for s in snr4
            ]
        )
        cm = np.array(
            [
                float(np.median([
                    r[field] for r in records
                    if r["group"] == "mono_L2_comoving"
                    and r["method_alias"] == "direct_comoving"
                    and r["SNR_dB"] == s
                ]))
                for s in snr4
            ]
        )
        tx = np.maximum(tx, 1e-14)
        cm = np.maximum(cm, 1e-14)
        b_tx = ax4.bar(base - widths / 2, tx, widths, color=c, alpha=0.62)
        b_cm = ax4.bar(base + widths / 2, cm, widths, color=c, alpha=0.95,
                       hatch="//")
        for b, v in zip(base, tx):
            ax4.text(b - widths / 2, v * 1.3, f"{v:.1e}", ha="center",
                     fontsize=5.5, rotation=90)
        for b, v in zip(base, cm):
            ax4.text(b + widths / 2, v * 1.3, f"{v:.1e}", ha="center",
                     fontsize=5.5, rotation=90)
    ax4.set_yscale("log")
    ax4.set_xticks(np.arange(len(snr4) * 2))
    ax4.set_xticklabels(
        [f"pose {int(s)}" for s in snr4] + [f"theta {int(s)}" for s in snr4],
        fontsize=8,
    )
    ax4.set_xlabel("grouped by metric and SNR (log y)")
    ax4.set_ylabel("median error")
    ax4.set_title("(4) L=2 mono direct: tx_only vs co_moving at SNR 20/30/40",
                  fontsize=10)
    ax4.grid(alpha=0.3, which="both", axis="y")
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor="#1f77b4", alpha=0.62, label="pose error, tx_only"),
        Patch(facecolor="#1f77b4", alpha=0.95, hatch="//",
              label="pose error, co_moving"),
        Patch(facecolor="#2ca02c", alpha=0.62, label="theta error, tx_only"),
        Patch(facecolor="#2ca02c", alpha=0.95, hatch="//",
              label="theta error, co_moving"),
    ]
    ax4.legend(handles=legend_handles, fontsize=6.5, loc="best")

    # 5: B_t realified singular values, mono + dipole models
    ax5 = fig.add_subplot(gs[2, 0])
    all_g = []
    for g in gauge_mono:
        all_g.append((f"L={g['L']} mono", g["realified_sv"]))
    for g in gauge_dipole:
        all_g.append((f"L={g['L']} dipole", g["realified_sv"]))
    xt = np.arange(len(all_g))
    gauge_colors = plt.cm.tab10(np.linspace(0.0, 1.0, len(all_g)))
    for i, (lab, svs) in enumerate(all_g):
        svs = np.maximum(svs, 1e-16)
        x = xt[i] + np.arange(len(svs)) - 0.3
        ax5.bar(x, svs, 0.25, color=gauge_colors[i], label=lab, alpha=0.9)
        for xi, v in zip(x, svs):
            ax5.text(xi, max(v, 1e-16) * 2.0, f"{v:.1e}", ha="center",
                     fontsize=4.6, rotation=90)
    ax5.set_yscale("log")
    ax5.set_xticks(xt)
    ax5.set_xticklabels([lab for lab, _ in all_g], fontsize=8)
    ax5.set_ylabel(r"$\sigma_j$ realified $B_t$ (log)")
    ax5.set_title(
        "(5) realified B_t singular values (theta-gauge break at row 3)",
        fontsize=10,
    )
    ax5.grid(alpha=0.3, which="both", axis="y")
    ax5.legend(fontsize=7, loc="upper left", ncol=2)
    ymin = float(min(np.maximum(np.asarray([v for _, sv in all_g for v in sv]), 1e-16)))
    ax5.set_ylim(ymin * 0.05, None)

    # 6: success counts heatmap (rows = group+method, cols = SNR)
    ax6 = fig.add_subplot(gs[2, 1])
    row_keys = []
    row_meta = []
    for g in ("mono_L2", "mono_L3", "mono_L2_comoving", "dipole_L2"):
        for alias in ("wrongpose", "direct_tx", "direct_comoving"):
            if alias == "direct_comoving" and g != "mono_L2_comoving":
                continue
            if alias in ("wrongpose", "direct_tx") and any(
                r["group"] == g and r["method_alias"] == alias
                for r in records
            ):
                row_meta.append((g, alias))
    mat = np.full((len(row_meta), len(SNR_COLS)), np.nan)
    for i, (g, alias) in enumerate(row_meta):
        for j, s in enumerate(SNR_COLS):
            sel = [
                r for r in records
                if r["group"] == g
                and r["method_alias"] == alias
                and float(r["SNR_dB"]) == s
            ]
            if sel:
                mat[i, j] = sum(r["success"] for r in sel)
    cmap = matplotlib.colormaps["YlGn"]
    cmap.set_bad("#eeeeee")
    masked = np.ma.masked_invalid(mat)
    im = ax6.imshow(masked, cmap=cmap, vmin=0, vmax=8, aspect="auto")
    ax6.set_xticks(np.arange(len(SNR_COLS)))
    ax6.set_xticklabels([f"{int(s)}" for s in SNR_COLS], fontsize=8)
    ax6.set_yticks(np.arange(len(row_meta)))
    ax6.set_yticklabels(
        [
            f"{GROUP_LABELS[g].replace(' tx_only', '').replace(' co_moving', '')} "
            f"{alias}"
            for g, alias in row_meta
        ],
        fontsize=7,
    )
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not np.isnan(mat[i, j]):
                ax6.text(j, i, f"{int(mat[i, j])}", ha="center", va="center",
                         fontsize=7)
    ax6.set_xlabel("SNR (dB)")
    ax6.set_title("(6) success counts out of 8 seeds", fontsize=10)
    cb = fig.colorbar(im, ax=ax6, fraction=0.046, pad=0.04)
    cb.set_label("n_success", fontsize=8)

    fig.suptitle(
        "E6b: robust multi-transmitter sweeps + directional phased-dipole "
        "transmitter model",
        y=0.995,
    )
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def _base_setup():
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    Phi = basis_matrix(xs, P)
    K = Phi.shape[1]
    assert K == 3
    return xs, h, G_D, Phi, K


def _run_monopole_sweep(
    xs, h, G_D, Phi, records: list[dict]
) -> None:
    """Step 2: monopole sweeps reusing E6 run_case/make_observation."""
    obs_cache = {}
    jobs = []
    for group, cfg in MONO_ARRAYS.items():
        if cfg["model"] != "monopole":
            continue
        tls = make_t_locals(cfg["angles"])
        for snr in cfg["snr_dB"]:
            for seed in NOISE_SEEDS:
                for method in ("wrongpose", "direct"):
                    jobs.append(
                        (group, tls, cfg["pose_mode"], method, float(snr), seed)
                    )

    def obs_key(group, mode, snr, seed):
        # angles and model are fixed per group -> group is enough to separate
        # arrays with identical (mode, SNR, seed) noise keys.
        return (group, mode, float(snr), int(seed))

    for group, tls, mode, method, snr, seed in jobs:
        key = obs_key(group, mode, snr, seed)
        if key not in obs_cache:
            obs_cache[key] = make_observation(
                ALPHA_TRUE,
                P_TRUE,
                mode,
                P_RX_FIXED,
                snr,
                seed,
                xs,
                h,
                P.k,
                G_D,
                Phi,
                tls,
                P,
            )

    print(f"\n[E6b] running {len(jobs)} monopole least-squares cases ...")
    for job_i, (group, tls, mode, method, snr, seed) in enumerate(jobs):
        d_obs, d_true = obs_cache[obs_key(group, mode, snr, seed)]
        ref_norm = float(np.linalg.norm(d_obs))
        rec = run_case(
            method,
            mode,
            snr,
            seed,
            "zero",
            P_INIT_ZERO,
            ALPHA_INIT,
            ALPHA_TRUE,
            P_TRUE,
            d_obs,
            ref_norm,
            None,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            tls,
            P,
            P_RX_FIXED,
        )
        records.append(finalize_record(rec, group))
        if (job_i + 1) % 24 == 0 or job_i == len(jobs) - 1:
            print(
                f"  [{job_i+1:03d}/{len(jobs)}] {group:>16s} "
                f"{method_alias_of(mode, method):>15s} SNR={snr:.0f} "
                f"seed={seed} pose_err={rec['pose_error']:.3e} "
                f"ok={int(rec['success'])}"
            )


def _run_dipole_nonlinear(
    xs, h, G_D, Phi, t_locals, records: list[dict]
) -> None:
    """Step 3 nonlinear tx_only runs with the L=2 dipole model."""
    mode = "tx_only"
    snrs = [20.0, 30.0, 40.0]
    obs_cache = {}
    jobs = []
    for snr in snrs:
        for seed in NOISE_SEEDS:
            for method in ("wrongpose", "direct"):
                jobs.append((method, snr, seed))
            obs_cache[(snr, seed)] = dipole_make_observation(
                ALPHA_TRUE,
                P_TRUE,
                mode,
                P_RX_FIXED,
                snr,
                seed,
                xs,
                h,
                P.k,
                G_D,
                Phi,
                t_locals,
                P,
            )

    print(
        f"\n[E6b] running {len(jobs)} dipole least-squares cases (tx_only) ..."
    )
    for job_i, (method, snr, seed) in enumerate(jobs):
        d_obs, d_true = obs_cache[(snr, seed)]
        ref_norm = float(np.linalg.norm(d_obs))
        rec = dipole_run_case(
            method,
            mode,
            snr,
            seed,
            "zero",
            P_INIT_ZERO,
            ALPHA_INIT,
            ALPHA_TRUE,
            P_TRUE,
            d_obs,
            ref_norm,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            t_locals,
            P,
            P_RX_FIXED,
        )
        records.append(finalize_record(rec, "dipole_L2"))
        if (job_i + 1) % 24 == 0 or job_i == len(jobs) - 1:
            print(
                f"  [{job_i+1:03d}/{len(jobs)}] dipole_L2 "
                f"{method_alias_of(mode, method):>15s} SNR={snr:.0f} "
                f"seed={seed} pose_err={rec['pose_error']:.3e} "
                f"ok={int(rec['success'])}"
            )


def main() -> dict:
    xs, h, G_D, Phi, K = _base_setup()

    # --- Step 3 setup: FD + gauge for the dipole model --------------------
    p_fd = np.array([0.08, -0.06, 0.05])
    alpha_fd = ALPHA_INIT
    t_fd_dipole = make_t_locals([0.7, 2.0])  # dipole centres
    dipole_fd = []
    for mode in ("tx_only", "co_moving"):
        dipole_fd += dipole_fd_check(
            alpha_fd,
            p_fd,
            t_fd_dipole,
            P,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            EPS_LIST,
            mode,
            P_RX_FIXED,
        )
    max_at_1e4 = max(
        r["max_rel_err"]
        for r in dipole_fd
        if abs(r["eps"] - 1e-4) < 1e-15
    )
    print(f"[E6b FD dipole] max rel err at eps=1e-4 = {max_at_1e4:.3e}")
    assert max_at_1e4 < 1e-5, (
        "dipole finite-difference verification FAILED -- do not run nonlinear"
    )
    for r in dipole_fd:
        r["per_col_rel_err"] = [float(v) for v in r["per_col_rel_err"]]

    # --- Step 3 gauge / rank analysis (L=1 vs L=2 dipoles) ----------------
    t_locals_dipole_by_L = {
        1: make_t_locals([0.7]),
        2: make_t_locals([0.7, 2.0]),
    }
    gauge_dipole, comov_dipole = dipole_rank_gauge_analysis(
        ALPHA_INIT, t_locals_dipole_by_L, P, xs, h, P.k, G_D, Phi
    )

    # E6 monopole gauge (loaded through helpers, used for plot panel 5).
    t_locals_mono_by_L = {
        1: make_t_locals([0.7]),
        2: make_t_locals([0.7, 2.0]),
        3: make_t_locals([0.7, 2.0, 1.3]),
    }
    gauge_mono = []
    p0 = np.zeros(3)
    for L in sorted(t_locals_mono_by_L):
        jb = multi_jacobians(
            ALPHA_INIT,
            p0,
            p0,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            t_locals_mono_by_L[L],
            P,
            include_receiver=False,
        )
        rank_c, real_sv, normed_real_sv = _real_ranks(jb["B_t"])
        c_sv = np.linalg.svd(jb["B_t"], compute_uv=False)
        res = _theta_gauge_residual(jb["B_t"])
        gauge_mono.append(
            {
                "L": int(L),
                "model": "monopole",
                "t_local_angles": [
                    float(np.arctan2(t[1], t[0]))
                    for t in t_locals_mono_by_L[L]
                ],
                "complex_rank": rank_c,
                "complex_sv": [float(v) for v in c_sv],
                "realified_sv": [float(v) for v in real_sv],
                "realified_colnormed_sv": [float(v) for v in normed_real_sv],
                "theta_residual_abs": res,
                "theta_residual_rel": float(
                    res / np.linalg.norm(jb["B_t"][:, 2])
                ),
            }
        )

    # --- Step 2 + Step 3 nonlinear records --------------------------------
    records: list[dict] = []
    _run_monopole_sweep(xs, h, G_D, Phi, records)
    _run_dipole_nonlinear(xs, h, G_D, Phi, t_fd_dipole, records)

    summary = make_summary(records)

    # --- Parameters block -------------------------------------------------
    parameters = {
        "experiment": "run_e6b_multitx_robust",
        "extends": "run_e6_multitx",
        "scene": {
            "N": P.N,
            "k": P.k,
            "M": P.M,
            "R_r": P.R_r,
            "R_t": P.R_t,
            "s": P.s,
            "K": K,
            "basis_centers": [
                list(c)
                for c in ((-0.15, 0.10), (0.20, -0.10), (0.00, 0.00))
            ],
        },
        "physics_true": {
            "alpha_true": [float(v) for v in ALPHA_TRUE],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_true": [float(v) for v in P_TRUE],
            "p_rx_fixed": [float(v) for v in P_RX_FIXED],
            "p_init_label_zero": [float(v) for v in P_INIT_ZERO],
        },
        "arrays": {
            key: {
                "model": cfg["model"],
                "angles": list(cfg["angles"]),
                "pose_mode": cfg["pose_mode"],
                "SNR_dB_values": list(cfg["snr_dB"]),
                "noise_seeds": list(NOISE_SEEDS),
            }
            for key, cfg in MONO_ARRAYS.items()
        },
        "dipole_model": {
            "kind": "tangent phased dipole pair per source",
            "d": float(DIPOLE_D),
            "phase": float(DIPOLE_PHASE),
            "weights": ["exp(+i*phase/2) plus", "exp(-i*phase/2) minus"],
            "orientation": "local angle + pi/2 (tangent to R_t circle)",
            "theta_velocity": "per-monopole R_tx @ R90(q_local) "
                              "(not centre velocity)",
        },
        "noise": {
            "convention": (
                "exact E6 make_observation: complex Gaussian raw vector "
                "rescaled to ||noise|| = 10^(-SNR/20) * ||d_true||"
            ),
            "seed_range": [int(min(NOISE_SEEDS)), int(max(NOISE_SEEDS))],
        },
        "least_squares": dict(LS_KWARGS),
        "success_gate": (
            "copied from E6 run_case: result.status > 0 and "
            "(result.optimality < 1e-8 or result.status in (1,2,3,4))"
        ),
        "method_alias": {
            "wrongpose": "wrongpose (alpha only at p=0)",
            "direct@tx_only": "direct_tx",
            "direct@co_moving": "direct_comoving",
            "note": (
                "E6 run_case method strings discovered as 'wrongpose' and "
                "'direct'; the prompt's direct_tx name maps to the E6 'direct' "
                "method with pose_mode='tx_only'"
            ),
        },
        "fd_dipole": {
            "t_local_center_angles": [0.7, 2.0],
            "p_fd": [float(v) for v in p_fd],
            "alpha_fd": [float(v) for v in alpha_fd],
            "eps_list": [float(v) for v in EPS_LIST],
            "requirement": "max per-column rel err at eps=1e-4 < 1e-5",
            "pass_at_eps_1e-4_max_rel": max_at_1e4,
        },
        "gauge_dipole": {
            "t_local_center_angles_by_L": {
                str(L): [
                    float(np.arctan2(t[1], t[0]))
                    for t in t_locals_dipole_by_L[L]
                ]
                for L in sorted(t_locals_dipole_by_L)
            },
            "analysis_at": {
                "alpha": [float(v) for v in ALPHA_INIT],
                "p": [0.0, 0.0, 0.0],
            },
        },
    }

    results = {
        "parameters": parameters,
        "fd_checks_dipole": dipole_fd,
        "gauge_partA_dipole": {
            "B_t_analysis": gauge_dipole,
            "comoving_B_total_sanity": comov_dipole,
            "monopole_B_t_reference": gauge_mono,
        },
        "runs": records,
        "summary_medians": summary,
    }

    json_path = HERE / "results_e6b_multitx_robust.json"
    png_path = HERE / "plot_e6b_multitx_robust.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2, allow_nan=False)
    make_plot(records, gauge_mono, gauge_dipole, png_path)

    with open(json_path) as fh:
        loaded = json.load(fh)
    assert len(loaded["runs"]) == len(records)
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    print("\n--- E6b text summary (median | IQR, medians over seeds) ---")
    hdr = (
        f"{'group':>16s} {'alias':>15s} {'SNR':>3s} "
        f"{'pose_err':>9s} {'theta_err':>9s} {'map_err':>9s} ok"
    )
    print(hdr)
    for s in summary:
        print(
            f"{s['group']:>16s} {s['method_alias']:>15s} "
            f"{int(s['SNR_dB']):>3d} "
            f"{s['median_pose_error']:.3e}|{s['iqr_pose_error']:.1e} "
            f"{s['median_theta_err']:.3e}|{s['iqr_theta_err']:.1e} "
            f"{s['median_map_error']:.3e}|{s['iqr_map_error']:.1e} "
            f"{s['n_success']}/{s['n_seeds']}"
        )

    n_fail = sum(not r["success"] for r in records)
    print(f"\nrecords={len(records)} failed={n_fail}")
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
