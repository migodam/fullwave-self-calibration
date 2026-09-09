"""E6 -- multi-transmitter / oriented-array theta-gauge experiment.

Test whether a multi-transmitter (or oriented-array) transmitter can break the
scalar point-transmitter theta gauge, so that transmitter-only self-calibration
identifies theta.

Physics/convention notes (identical to the existing harness):
  * J = A^{-1}(chi .* u_inc), A = I - diag(chi) G_D, u = u_inc + G_D J.
  * d_l = G_s J_l for source l; the multi data vector is the vertical stack.
  * Receiver rows rotate as v @ R.T (body tangent -> global) and transmitter
    column vectors as R @ w (body -> global), R = [[c,-s],[s,c]].
  * Realification of a complex m-vector is [Re; Im] (2m).  Jacobian real
    matrices are the same vertical stack.

Outputs: results_e6.json, plot_e6.png
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
    transmitter_direction,
    data_matrix,
    receiver_data_derivative,
)
from run_e2_e4 import (
    realify_complex_cols,
    realify_real_cols,
    norm_cols,
    proj_complement,
)
from run_e5 import basis_matrix, _transmitter_incident_derivative, _realified_residual


HERE = Path(__file__).resolve().parent
DIRS = ("rx", "ry", "theta")
COORD_NAMES = ("tx", "ty", "theta")
TOL = 1e-10
HIDDEN_SV_REL = 1e-8
METHOD_ORDER = ("wrongpose", "direct", "reduced_r4", "reduced_r6")

LS_KWARGS = dict(
    method="trf",
    x_scale="jac",
    max_nfev=1000,
    ftol=1e-10,
    xtol=1e-10,
    gtol=1e-10,
)


def _rot(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def make_t_locals(angles, R_t=2.0) -> np.ndarray:
    """Lx2 local transmitter positions from polar angles."""
    ang = np.asarray(angles, dtype=float)
    return R_t * np.stack((np.cos(ang), np.sin(ang)), axis=1)


def _mode_split(p: np.ndarray, pose_mode: str, p_rx_fixed: np.ndarray):
    """Return (p_tx, p_rx) for the requested pose mode."""
    p = np.asarray(p, dtype=float)
    if pose_mode == "tx_only":
        return p, np.asarray(p_rx_fixed, dtype=float)
    if pose_mode == "co_moving":
        return p, p
    raise ValueError(f"unknown pose_mode {pose_mode!r}")


# ---------------------------------------------------------------------------
# Multi-transmitter forward model and analytic Jacobians
# ---------------------------------------------------------------------------
def multi_forward(
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
    """Stacked L-source data d = [d_0; d_1; ...] (complex length L*M)."""
    chi = Phi @ alpha
    y, _ = receiver_positions(p_rx, P)
    G_s = data_matrix(y, xs, h, k)
    R = _rot(p_tx[2])
    ds = []
    for tl in t_locals:
        t = R @ tl + p_tx[:2]
        u_inc = incident_field(xs, t, k)
        J = solve_current(chi, u_inc, G_D)
        ds.append(G_s @ J)
    d_all = np.concatenate(ds)
    return (d_all, ds) if want_sources else d_all


def multi_jacobians(
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
    """Analytic Jacobians for the stacked multi-transmitter model.

    include_receiver=True means the receiver array pose varies along with the
    transmitter pose (co-moving calibration: p_rx == p_tx); then
    B_total = B_t + B_r.  For transmitter-only self-calibration the receiver
    array is fixed and include_receiver=False, so B_total = B_t.
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
        t = R_tx @ tl + p_tx[:2]
        u_inc = incident_field(xs, t, k)
        J = lu_solve((lu, piv), chi * u_inc)
        u = total_field(u_inc, J, G_D)
        d_l.append(G_s @ J)

        J_chi = lu_solve((lu, piv), np.diag(u))  # Np x Np
        J_alpha_l.append(G_s @ J_chi @ Phi)

        B_t = np.empty((P.M, 3), dtype=complex)
        for c, direction in enumerate(DIRS):
            w = transmitter_direction(direction, tl)
            if direction == "theta":
                w = R_tx @ w
            du_inc = _transmitter_incident_derivative(xs, t, w, k)
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
# Visible-pose subspace (generalized to stacked multi-transmitter rows)
# ---------------------------------------------------------------------------
def multi_visible_pose_subspace(
    alpha: np.ndarray,
    p: np.ndarray,
    r: int,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    t_locals: np.ndarray,
    P: Params,
    pose_mode: str,
    p_rx_fixed: np.ndarray,
) -> dict:
    """Reduced visible pose basis (3 x n_vis) at (alpha, p) for multi-tx rows."""
    L = t_locals.shape[0]
    p_tx, p_rx = _mode_split(p, pose_mode, p_rx_fixed)
    jb = multi_jacobians(
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
    y, _ = receiver_positions(p_rx, P)
    G_s = data_matrix(y, xs, h, k)
    _, _, Vh = np.linalg.svd(G_s, full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    Q = G_s @ V_r
    Qtile = np.vstack([Q] * L)
    Hn = norm_cols(
        np.hstack(
            [
                realify_complex_cols(Qtile),
                realify_real_cols(jb["J_alpha"]),
            ]
        )
    )
    P_perp = proj_complement(Hn, TOL)
    B_red = P_perp @ jb["B_real"]
    _, sb, Vhb = np.linalg.svd(B_red, full_matrices=False)
    sb = np.asarray(sb, dtype=float)
    threshold = HIDDEN_SV_REL * sb[0] if sb[0] > 0.0 else 0.0
    n_vis = int(np.count_nonzero(sb > threshold))
    hidden_rank = 3 - n_vis
    V_vis = Vhb[:n_vis].T  # 3 x n_vis
    hidden_directions = [Vhb[i].tolist() for i in range(n_vis, 3)]
    return {
        "rank": int(r),
        "n_vis": n_vis,
        "hidden_rank": hidden_rank,
        "hid_sv": [float(v) for v in sb],
        "threshold": float(threshold),
        "V_vis": V_vis,
        "hidden_directions": hidden_directions,
    }


# ---------------------------------------------------------------------------
# Finite-difference verification
# ---------------------------------------------------------------------------
def fd_check(
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
    """Centered finite differences vs analytic B columns for one pose mode."""
    p0 = np.asarray(p0, dtype=float)
    p_tx0, p_rx0 = _mode_split(p0, pose_mode, p_rx_fixed)
    jb = multi_jacobians(
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
    ana = jb["B_t"] if pose_mode == "tx_only" else jb["B_total"]

    rows = []
    for eps in eps_list:
        per_col = []
        for c in range(3):
            ep = np.zeros(3)
            ep[c] = eps
            p_plus, p_rx_plus = _mode_split(p0 + ep, pose_mode, p_rx_fixed)
            p_minus, p_rx_minus = _mode_split(p0 - ep, pose_mode, p_rx_fixed)
            d_plus = multi_forward(
                alpha, p_plus, p_rx_plus, xs, h, k, G_D, Phi, t_locals, P
            )
            d_minus = multi_forward(
                alpha, p_minus, p_rx_minus, xs, h, k, G_D, Phi, t_locals, P
            )
            fd_col = _realified_residual(d_plus, d_minus) / (2.0 * eps)
            ana_col = np.concatenate([ana[:, c].real, ana[:, c].imag])
            rel = float(np.linalg.norm(ana_col - fd_col) / np.linalg.norm(fd_col))
            per_col.append(rel)
        rows.append(
            {
                "pose_mode": pose_mode,
                "eps": float(eps),
                "per_col_rel_err": per_col,
                "max_rel_err": float(max(per_col)),
            }
        )
        print(
            f"[E6 FD] mode={pose_mode:>9s} eps={eps:.0e} "
            f"max_rel_err={max(per_col):.3e} per_col="
            + ",".join(f"{v:.2e}" for v in per_col)
        )
    return rows


# ---------------------------------------------------------------------------
# Rank / gauge analysis (Part A)
# ---------------------------------------------------------------------------
def _theta_gauge_residual(B_t: np.ndarray) -> float:
    """|theta_col - U (U^H theta_col)| with U from QR of the two translations."""
    Q, _ = np.linalg.qr(B_t[:, :2], mode="reduced")
    res = B_t[:, 2] - Q @ (Q.conj().T @ B_t[:, 2])
    return float(np.linalg.norm(res))


def _real_ranks(B: np.ndarray) -> tuple[int, np.ndarray, np.ndarray]:
    """Complex rank + realified singular values + column-normed real sv."""
    Bc = np.asarray(B, dtype=complex)
    Br = np.vstack([Bc.real, Bc.imag])
    return (
        int(np.linalg.matrix_rank(Bc)),
        np.linalg.svd(Br, compute_uv=False),
        np.linalg.svd(norm_cols(Br), compute_uv=False),
    )


def rank_gauge_analysis(
    alpha_init: np.ndarray,
    t_locals_by_L: dict,
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
) -> tuple[list[dict], list[dict]]:
    """Per-L complex rank, singular values and theta-gauge residual of B_t."""
    out = []
    p0 = np.zeros(3)
    for L in sorted(t_locals_by_L):
        tl = t_locals_by_L[L]
        jb = multi_jacobians(
            alpha_init, p0, p0, xs, h, k, G_D, Phi, tl, P, include_receiver=False
        )
        B_t = jb["B_t"]
        rank_c, real_sv, normed_real_sv = _real_ranks(B_t)
        c_sv = np.linalg.svd(B_t, compute_uv=False)
        res = _theta_gauge_residual(B_t)
        res_rel = float(res / np.linalg.norm(B_t[:, 2]))
        out.append(
            {
                "L": int(L),
                "t_local_angles": [
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
            f"[E6 gauge] L={L} complex_rank={rank_c} complex_sv="
            + ",".join(f"{v:.3e}" for v in c_sv)
            + f" theta_res_abs={res:.3e} rel={res_rel:.3e}"
        )

    # Co-moving sanity ranks (B_total = B_r + B_t) for L=1 and L=2.
    comov = []
    for L in (1, 2):
        tl = t_locals_by_L[L]
        jb = multi_jacobians(
            alpha_init, p0, p0, xs, h, k, G_D, Phi, tl, P, include_receiver=True
        )
        rank_c, real_sv, _ = _real_ranks(jb["B_total"])
        comov.append(
            {
                "L": int(L),
                "B_total_complex_rank": rank_c,
                "B_total_realified_sv": [float(v) for v in real_sv],
            }
        )
        print(
            f"[E6 gauge] co-moving L={L} B_total complex_rank={rank_c} "
            "realified_sv=" + ",".join(f"{v:.3e}" for v in real_sv)
        )
    return out, comov


# ---------------------------------------------------------------------------
# Nonlinear self-calibration (Part B)
# ---------------------------------------------------------------------------
def make_observation(
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
    """d_obs = d_true + scaled complex Gaussian noise (exact norm ratio)."""
    p_tx, p_rx = _mode_split(p_true, pose_mode, p_rx_fixed)
    d_true = multi_forward(
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


def run_case(
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
    subspace: dict | None,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    t_locals: np.ndarray,
    P: Params,
    p_rx_fixed: np.ndarray,
) -> dict:
    """One scipy least_squares solve; returns the per-run record."""
    K = Phi.shape[1]
    V_vis = subspace["V_vis"] if subspace is not None else None

    def residual(theta: np.ndarray) -> np.ndarray:
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
        elif method == "direct":
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
        else:  # reduced
            alpha = np.asarray(theta[:K], dtype=float)
            q = np.asarray(theta[K:], dtype=float)
            p = p_init + V_vis @ q
        p_tx, p_rx = _mode_split(p, pose_mode, p_rx_fixed)
        return _realified_residual(
            multi_forward(
                alpha, p_tx, p_rx, xs, h, k, G_D, Phi, t_locals, P
            ),
            d_obs,
        )

    def jacobian(theta: np.ndarray) -> np.ndarray:
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
        elif method == "direct":
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
        else:
            alpha = np.asarray(theta[:K], dtype=float)
            q = np.asarray(theta[K:], dtype=float)
            p = p_init + V_vis @ q
        p_tx, p_rx = _mode_split(p, pose_mode, p_rx_fixed)
        jb = multi_jacobians(
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
        if method == "direct":
            return np.hstack([jb["J_alpha_real"], jb["B_real"]])
        return np.hstack([jb["J_alpha_real"], jb["B_real"] @ V_vis])

    if method == "wrongpose":
        x0 = np.asarray(alpha_init, dtype=float)
    elif method == "direct":
        x0 = np.concatenate(
            [np.asarray(alpha_init, float), np.asarray(p_init, float)]
        )
    else:
        n_vis_solver = int(subspace["n_vis"])
        x0 = np.concatenate(
            [np.asarray(alpha_init, float), np.zeros(n_vis_solver)]
        )

    result = least_squares(residual, x0, jac=jacobian, **LS_KWARGS)

    if method == "wrongpose":
        alpha_est = np.asarray(result.x[:K], dtype=float)
        p_est = np.zeros(3)
    elif method == "direct":
        alpha_est = np.asarray(result.x[:K], dtype=float)
        p_est = np.asarray(result.x[K:], dtype=float)
    else:
        alpha_est = np.asarray(result.x[:K], dtype=float)
        q_est = np.asarray(result.x[K:], dtype=float)
        p_est = p_init + V_vis @ q_est

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
    }
    if method.startswith("reduced_"):
        record["n_vis"] = int(subspace["n_vis"])
        record["hid_sv"] = [float(v) for v in subspace["hid_sv"]]
        record["hidden_rank"] = int(subspace["hidden_rank"])
    else:
        record["n_vis"] = None
        record["hid_sv"] = None
        record["hidden_rank"] = None
    return record


def summarize(records: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    order = []
    for rec in records:
        key = (
            rec["pose_mode"],
            rec["method"],
            rec["SNR_dB"],
            rec["p_init_label"],
        )
        if key not in groups:
            order.append(key)
        groups[key].append(rec)

    rows = []
    for key in order:
        sel = groups[key]

        def med(f):
            return float(np.median([f(r) for r in sel]))

        def iqr(f):
            return float(
                np.percentile([f(r) for r in sel], 75)
                - np.percentile([f(r) for r in sel], 25)
            )

        rows.append(
            {
                "pose_mode": key[0],
                "method": key[1],
                "SNR_dB": key[2],
                "p_init_label": key[3],
                "n_seeds": len(sel),
                "n_success": int(sum(r["success"] for r in sel)),
                "median_pose_error": med(lambda r: r["pose_error"]),
                "iqr_pose_error": iqr(lambda r: r["pose_error"]),
                "median_theta_err": med(lambda r: r["theta_err"]),
                "median_map_error": med(lambda r: r["map_error"]),
                "iqr_map_error": iqr(lambda r: r["map_error"]),
                "median_final_residual": med(lambda r: r["final_residual"]),
                "iqr_final_residual": iqr(lambda r: r["final_residual"]),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def make_plot(
    gauge: list[dict],
    records: list[dict],
    summary: list[dict],
    png_path: Path,
) -> None:
    colors = {
        "wrongpose": "#d62728",
        "direct": "#1f77b4",
        "reduced_r4": "#2ca02c",
        "reduced_r6": "#9467bd",
    }
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.4))

    # Left: complex singular values of B_t per L (rank/gauge panel).
    ax = axes[0]
    for gi, g in enumerate(gauge):
        x = gi + np.arange(3) - 0.3
        sv = np.maximum(np.asarray(g["complex_sv"]), 1e-16)
        ax.bar(x, sv, 0.25, label=f'L={g["L"]}', alpha=0.85)
        for xi, v in zip(x, sv):
            ax.text(
                xi, v * 1.4, f"{v:.1e}", ha="center", fontsize=6.5, rotation=90
            )
    ax.set_yscale("log")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["L=1", "L=2", "L=3"])
    ax.set_ylabel(r"$\sigma_j(B_t)$ (complex)")
    ax.set_title(
        "E6 A: rank/gauge singular values of $B_t$\n"
        "(theta in span of tx/ty iff L=1)"
    )
    ax.grid(axis="y", which="both", alpha=0.3)
    ax.legend(fontsize=8)
    for gi, g in enumerate(gauge):
        ax.text(
            gi,
            2e-15 if gi == 0 else 2e-14,
            f"rank={g['complex_rank']} res={g['theta_residual_abs']:.1e}",
            ha="center",
            fontsize=7.5,
        )

    # Middle: transmitter-only pose error by method at SNR 30/20 (zero init).
    ax = axes[1]

    def sel(mode, snr):
        return [
            s
            for s in summary
            if s["pose_mode"] == mode
            and s["SNR_dB"] == snr
            and s["p_init_label"] == "zero"
        ]

    tx30 = sel("tx_only", 30.0)
    tx20 = sel("tx_only", 20.0)
    x30 = np.arange(len(METHOD_ORDER))
    m30 = {s["method"]: s["median_pose_error"] for s in tx30}
    m20 = {s["method"]: s["median_pose_error"] for s in tx20}
    all_v = []
    for xv, m, off, alpha in (
        (x30 - 0.19, m30, 0.0, 0.95),
        (x30 + 0.19, m20, 0.0, 0.6),
    ):
        for j, meth in enumerate(METHOD_ORDER):
            if meth not in m:
                continue
            v = max(m[meth], 1e-12)
            all_v.append(v)
            ax.bar(
                xv[j],
                v,
                0.34,
                color=colors[meth],
                alpha=alpha,
                label=(
                    ("SNR 30" if alpha > 0.9 else "SNR 20")
                    if j == 0
                    else None
                ),
            )
            ax.text(
                xv[j],
                v * 1.35,
                f"{v:.1e}",
                ha="center",
                fontsize=6.2,
                rotation=90,
            )
    for meth in METHOD_ORDER:
        pts = [
            r
            for r in records
            if r["pose_mode"] == "tx_only"
            and r["method"] == meth
            and r["SNR_dB"] == 30.0
            and r["p_init_label"] == "zero"
        ]
        for r in pts:
            j = METHOD_ORDER.index(meth)
            ax.scatter(
                j - 0.19 + 0.04 * r["noise_seed"],
                max(r["pose_error"], 1e-13),
                marker=".",
                color="black",
                s=12,
                alpha=0.45,
                zorder=3,
            )
    ax.set_yscale("log")
    ax.set_xticks(x30)
    ax.set_xticklabels(list(METHOD_ORDER))
    ax.set_ylabel("median pose error")
    ax.set_title(
        "E6 B: transmitter-only pose error\n(zero init; medians over seeds)"
    )
    ax.grid(axis="y", which="both", alpha=0.3)
    if all_v:
        ax.set_ylim(0.4 * min(all_v), 3.0 * max(all_v))
    ax.legend(fontsize=8)

    # Right: transmitter-only theta error at SNR 30 (zero init).
    ax = axes[2]
    vals30 = {}
    for s in tx30:
        vals30[s["method"]] = s["median_theta_err"]
    xr = np.arange(len(METHOD_ORDER))
    for j, meth in enumerate(METHOD_ORDER):
        v = max(vals30[meth], 1e-14)
        ax.bar(xr[j], v, 0.5, color=colors[meth], alpha=0.88)
        ax.text(
            xr[j],
            v * 1.5,
            f"{v:.1e}",
            ha="center",
            fontsize=7,
            rotation=90,
        )
    for r in records:
        if (
            r["pose_mode"] == "tx_only"
            and r["SNR_dB"] == 30.0
            and r["p_init_label"] == "zero"
        ):
            ax.scatter(
                METHOD_ORDER.index(r["method"]) + 0.04 * r["noise_seed"],
                max(r["theta_err"], 1e-15),
                marker=".",
                color="black",
                s=14,
                alpha=0.5,
                zorder=3,
            )
    ax.axhline(0.05, color="k", ls=":", lw=0.8)
    ax.text(3.45, 0.052, "true $\\theta=0.05$", fontsize=7, ha="right")
    ax.set_yscale("log")
    ax.set_xticks(xr)
    ax.set_xticklabels(list(METHOD_ORDER))
    ax.set_ylabel(r"median $|\theta_{\rm est}-0.05|$")
    ax.set_title(
        "E6 B: transmitter-only theta error, SNR 30\n"
        "(direct/reduced recover $\\theta$; wrongpose stays at 0)"
    )
    ax.grid(axis="y", which="both", alpha=0.3)
    ax.set_ylim(1e-14, 1.0)

    fig.suptitle(
        "E6: oriented multi-transmitter array vs "
        "scalar point-transmitter theta gauge",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def main() -> dict:
    P = Params(N=16, k=12.0, M=8, R_r=1.6, R_t=2.0, phi0=0.0, s=0.12)
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    Phi = basis_matrix(xs, P)
    K = Phi.shape[1]
    assert K == 3

    t_fd = make_t_locals([0.3, 2.0])
    p_fd = np.array([0.08, -0.06, 0.05])
    p_rx_fixed = np.zeros(3)
    alpha_fd = np.array([1.0, 1.0, 0.0])
    eps_list = (1e-2, 1e-3, 1e-4, 1e-5)

    # --- Section 2: finite-difference verification (must pass) ------------
    fd_rows = []
    for mode in ("tx_only", "co_moving"):
        fd_rows += fd_check(
            alpha_fd,
            p_fd,
            t_fd,
            P,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            eps_list,
            mode,
            p_rx_fixed,
        )
    max_at_1e4 = max(
        r["max_rel_err"] for r in fd_rows if abs(r["eps"] - 1e-4) < 1e-15
    )
    print(f"[E6 FD] max rel err at eps=1e-4 = {max_at_1e4:.3e}")
    assert max_at_1e4 < 1e-5, "finite-difference verification FAILED"
    for r in fd_rows:
        r["per_col_rel_err"] = [float(v) for v in r["per_col_rel_err"]]

    # --- Section 3: rank/gauge analysis (Part A) --------------------------
    alpha_init = np.array([1.0, 1.0, 0.0])
    t_locals_by_L = {
        1: make_t_locals([0.7]),
        2: make_t_locals([0.7, 2.0]),
        3: make_t_locals([0.7, 2.0, 1.3]),
    }
    gauge, comov_rank = rank_gauge_analysis(
        alpha_init, t_locals_by_L, P, xs, h, P.k, G_D, Phi
    )

    # --- Section 4: nonlinear self-calibration (Part B) -------------------
    alpha_true = np.array([1.5, 2.0, 0.0])
    p_true = np.array([0.08, -0.06, 0.05])
    t_b = t_fd  # L=2, t_locals from the FD block
    L_b = t_b.shape[0]

    p_inits = {
        "zero": np.zeros(3),
        "pert_pos": np.array([0.18, 0.04, 0.15]),
    }

    # Subspace cache keyed by (pose_mode, p_init_label, r).
    subspaces = {}
    for mode in ("tx_only", "co_moving"):
        for label, p_init in p_inits.items():
            for r in (4, 6):
                subspaces[(mode, label, r)] = multi_visible_pose_subspace(
                    alpha_init,
                    p_init,
                    r,
                    xs,
                    h,
                    P.k,
                    G_D,
                    Phi,
                    t_b,
                    P,
                    mode,
                    p_rx_fixed,
                )
    print("\n[E6] visible-pose subspaces at (alpha_init, p_init):")
    for key, sp in subspaces.items():
        print(
            f"  mode={key[0]:>9s} init={key[1]:>8s} r={key[2]} "
            f"n_vis={sp['n_vis']} hid_sv="
            + ",".join(f"{v:.2e}" for v in sp["hid_sv"])
        )

    jobs = []
    for mode in ("tx_only", "co_moving"):
        for method in METHOD_ORDER:
            for seed in (0, 1, 2, 3):
                jobs.append((method, mode, 30.0, "zero", seed))
    for method in ("direct", "reduced_r4"):
        for seed in (0, 1, 2, 3):
            jobs.append((method, "tx_only", 20.0, "zero", seed))
    for method in ("direct", "reduced_r4", "reduced_r6"):
        for seed in (0, 1):
            jobs.append((method, "tx_only", 30.0, "pert_pos", seed))

    # Observation cache by (pose_mode, snr_db, seed) so every method sees the
    # same noisy data vector within a case.
    obs_cache = {}
    for method, mode, snr, label, seed in jobs:
        key = (mode, snr, seed)
        if key not in obs_cache:
            obs_cache[key] = make_observation(
                alpha_true,
                p_true,
                mode,
                p_rx_fixed,
                snr,
                seed,
                xs,
                h,
                P.k,
                G_D,
                Phi,
                t_b,
                P,
            )

    records = []
    print(f"\n[E6] running {len(jobs)} least-squares cases ...")
    for job_i, (method, mode, snr, label, seed) in enumerate(jobs):
        d_obs, d_true = obs_cache[(mode, snr, seed)]
        ref_norm = float(np.linalg.norm(d_obs))
        reduced_r = (
            int(method.split("_")[1][1:])
            if method.startswith("reduced_")
            else None
        )
        subspace = (
            subspaces[(mode, label, reduced_r)]
            if reduced_r is not None
            else None
        )
        rec = run_case(
            method,
            mode,
            snr,
            seed,
            label,
            p_inits[label],
            alpha_init,
            alpha_true,
            p_true,
            d_obs,
            ref_norm,
            subspace,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            t_b,
            P,
            p_rx_fixed,
        )
        records.append(rec)
        print(
            f"  [{job_i+1:02d}/{len(jobs)}] {mode:>9s} {method:>11s} "
            f"SNR={snr:.0f} {label:>8s} seed={seed} "
            f"pose_err={rec['pose_error']:.3e} "
            f"theta_err={rec['theta_err']:.3e} "
            f"map_err={rec['map_error']:.3e} "
            f"res={rec['final_residual']:.3e} "
            f"nfev={rec['nfev']} opt={rec['optimality']:.1e} "
            f"ok={int(rec['success'])}"
        )

    summary = summarize(records)

    # --- Section 5: summary and artifacts ---------------------------------
    jsonable_subspaces = {}
    for key, sp in subspaces.items():
        jsonable_subspaces[f"{key[0]}|{key[1]}|r{key[2]}"] = {
            "pose_mode": key[0],
            "p_init_label": key[1],
            "r": key[2],
            "rank": int(sp["rank"]),
            "n_vis": int(sp["n_vis"]),
            "hidden_rank": int(sp["hidden_rank"]),
            "hid_sv": [float(v) for v in sp["hid_sv"]],
            "threshold": float(sp["threshold"]),
            "hidden_directions": [
                list(v) for v in sp["hidden_directions"]
            ],
            "V_vis_cols": [
                sp["V_vis"][:, j].tolist()
                for j in range(sp["V_vis"].shape[1])
            ],
        }

    results = {
        "experiment": "run_e6_multitx",
        "parameters": {
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
            "fd": {
                "t_locals_angles": [0.3, 2.0],
                "t_locals_radius": float(P.R_t),
                "p_fd": [float(v) for v in p_fd],
                "alpha_fd": [float(v) for v in alpha_fd],
                "eps_list": [float(v) for v in eps_list],
                "pass_at_eps_1e-4_max_rel": max_at_1e4,
            },
            "partB": {
                "t_locals_angles": [0.3, 2.0],
                "L": L_b,
                "alpha_true": [float(v) for v in alpha_true],
                "alpha_init": [float(v) for v in alpha_init],
                "p_true": [float(v) for v in p_true],
                "p_inits": {
                    k: [float(v) for v in vv] for k, vv in p_inits.items()
                },
                "SNR_dB_values": [20.0, 30.0],
                "noise_seeds_30": [0, 1, 2, 3],
                "least_squares": dict(LS_KWARGS),
            },
        },
        "finite_difference_checks": fd_rows,
        "gauge_partA": {
            "B_t_analysis": gauge,
            "comoving_B_total_sanity": comov_rank,
        },
        "visible_pose_subspaces": jsonable_subspaces,
        "records": records,
        "summary_medians": summary,
    }

    json_path = HERE / "results_e6.json"
    png_path = HERE / "plot_e6.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2, allow_nan=False)
    make_plot(gauge, records, summary, png_path)

    with open(json_path) as fh:
        loaded = json.load(fh)
    assert len(loaded["records"]) == len(jobs)
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    print("\n--- E6 text summary (medians over seeds) ---")
    hdr = (
        f"{'pose_mode':>9s} {'method':>11s} {'SNR':>3s} {'init':>8s} "
        f"{'pose_err':>9s} {'theta_err':>9s} {'map_err':>9s} "
        f"{'resid':>9s} ok"
    )
    print(hdr)
    for s in summary:
        print(
            f"{s['pose_mode']:>9s} {s['method']:>11s} "
            f"{int(s['SNR_dB']):>3d} {s['p_init_label']:>8s} "
            f"{s['median_pose_error']:.3e} {s['median_theta_err']:.3e} "
            f"{s['median_map_error']:.3e} "
            f"{s['median_final_residual']:.3e} "
            f"{s['n_success']}/{s['n_seeds']}"
        )

    print("\nArtifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
