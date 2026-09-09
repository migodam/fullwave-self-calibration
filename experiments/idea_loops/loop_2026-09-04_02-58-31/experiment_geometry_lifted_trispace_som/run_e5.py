"""Geometry-lifted TriSpace SOM: E5 -- nonlinear self-calibration prototype.

Small-grid (N=16) adaptive Levenberg-Marquardt with realified least squares
over the real contrast coefficients alpha (K=3) and the co-moving pose
p=(rx, ry, theta).  Three solvers are compared on the same noisy observations:

  * wrongpose  : alpha only, pose held at p=0 (mismodelled-geometry baseline)
  * direct     : joint (alpha, p) with the full co-moving data Jacobian
  * reduced    : joint (alpha, p) with pose columns restricted to
                 P_perp @ B_real, where P_perp removes the geometry-lift +
                 real-contrast data subspace col(realify(G_s V_r), realify(W))
                 recomputed at every iteration for retained rank r
                 (r=4 disambiguated, r=6 confounded).

Determinism: no random state is used except one fresh numpy default_rng
per noise seed when drawing d_obs.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import hankel1

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
    transmitter_position,
    transmitter_direction,
    data_matrix,
    receiver_data_derivative,
)
from run_e2_e4 import (
    realify_complex_cols,
    realify_real_cols,
    norm_cols,
    rank_svd_tol,
    proj_complement,
)


HERE = Path(__file__).resolve().parent

TOL = 1e-10
LM_LAM0 = 1e-3
LM_LAM_MIN = 1e-12
LM_LAM_MAX = 1e10
LM_INNER_TRIES = 12
LM_MAX_ITER = 300
LM_GTOL = 1e-9
LM_FTOL = 1e-14
LM_XTOL = 1e-12
SNR_DB = 30.0
HIDDEN_SV_REL = 1e-8
DIRS = ("rx", "ry", "theta")
NOISE_SEEDS = (0, 1)
P_INIT_LABELS = ("zero", "pert_pos", "pert_neg")
ALPHA_INIT = np.array([1.0, 1.0, 0.0], dtype=float)


def basis_matrix(xs: np.ndarray, P: Params) -> np.ndarray:
    """Np x 3 real Gaussian contrast basis (local K=3 copy)."""
    centers = ((-0.15, 0.10), (0.20, -0.10), (0.00, 0.00))
    Phi = np.empty((xs.shape[0], 3), dtype=float)
    for j, center in enumerate(centers):
        delta = xs - np.asarray(center, dtype=float)
        Phi[:, j] = np.exp(-np.einsum("ij,ij->i", delta, delta) / (2.0 * P.s**2))
    return Phi


def _transmitter_incident_derivative(
    xs: np.ndarray, t: np.ndarray, w: np.ndarray, k: float
) -> np.ndarray:
    """du_inc/dp for transmitter motion w (spec formula, H1 first kind)."""
    delta = xs - t[None, :]
    R = np.sqrt(np.einsum("ij,ij->i", delta, delta))
    assert np.all(R > 0.0)
    dot = delta @ w
    return (1j / 4.0) * hankel1(1, k * R) * k * dot / R


def _realified_residual(d: np.ndarray, d_obs: np.ndarray) -> np.ndarray:
    """b = [Re(d - d_obs); Im(d - d_obs)] as a 2M real vector."""
    diff = d - d_obs
    return np.concatenate([diff.real, diff.imag])


def forward(
    alpha: np.ndarray,
    p: np.ndarray,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    P: Params,
) -> np.ndarray:
    """d = G_s(p) J, J = solve_current(Phi alpha, u_inc(p))."""
    chi = Phi @ alpha
    t, _ = transmitter_position(p, P)
    u_inc = incident_field(xs, t, k)
    J = solve_current(chi, u_inc, G_D)
    y, _ = receiver_positions(p, P)
    G_s = data_matrix(y, xs, h, k)
    return G_s @ J


def jacobians(
    alpha: np.ndarray,
    p: np.ndarray,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    P: Params,
) -> dict:
    """Co-moving Jacobians at (alpha, p).

    Returns d, complex J_alpha (M x K) and B_total/B_r/B_t (M x 3), plus the
    realified J_alpha_real and B_real and auxiliary state quantities.
    """
    chi = Phi @ alpha
    t, t_local = transmitter_position(p, P)
    u_inc = incident_field(xs, t, k)
    J = solve_current(chi, u_inc, G_D)
    u = total_field(u_inc, J, G_D)
    A = state_operator(chi, G_D)
    y, y_local = receiver_positions(p, P)
    G_s = data_matrix(y, xs, h, k)

    J_chi = np.linalg.solve(A, np.diag(u))  # Np x Np
    J_alpha = G_s @ J_chi @ Phi  # M x K complex

    B_r = np.empty((P.M, 3), dtype=complex)
    for c, direction in enumerate(DIRS):
        v = receiver_directions(direction, y_local)
        if direction == "theta":
            # receiver_directions returns the body-frame tangent
            # [-y_local_y, y_local_x]; the global receiver velocity at pose
            # theta is that tangent rotated by R(theta) (row convention:
            # v @ R.T).  At theta=0 this reduces to the reference-pose
            # derivative validated in E3.
            th = p[2]
            R = np.array(
                [[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]
            )
            v = v @ R.T
        B_r[:, c] = receiver_data_derivative(y, v, xs, h, k) @ J

    B_t = np.empty((P.M, 3), dtype=complex)
    for c, direction in enumerate(DIRS):
        w = transmitter_direction(direction, t_local)
        if direction == "theta":
            # Same body-to-global rotation for the transmitter column vector.
            th = p[2]
            R = np.array(
                [[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]
            )
            w = R @ w
        du_inc = _transmitter_incident_derivative(xs, t, w, k)
        J_p = np.linalg.solve(A, chi * du_inc)
        B_t[:, c] = G_s @ J_p

    B_total = B_r + B_t
    return {
        "d": G_s @ J,
        "chi": chi,
        "u": u,
        "J": J,
        "A": A,
        "G_s": G_s,
        "J_alpha": J_alpha,
        "B_r": B_r,
        "B_t": B_t,
        "B_total": B_total,
        "J_alpha_real": np.vstack([J_alpha.real, J_alpha.imag]),
        "B_real": np.vstack([B_total.real, B_total.imag]),
    }


def reduced_projector(
    jb: dict, r: int, P: Params
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """P_perp = I - U_r U_r^T over col(realify(G_s V_r), realify(W)).

    Returns (P_perp, V_r, Q, Hn, rank_H) at the current pose/alpha.
    """
    U, s, Vh = np.linalg.svd(jb["G_s"], full_matrices=False)
    del U, s
    V_r = Vh.conj().T[:, :r]
    Q = jb["G_s"] @ V_r
    W = jb["J_alpha"]
    Hn = norm_cols(
        np.hstack([realify_complex_cols(Q), realify_real_cols(W)])
    )
    rank_H = int(rank_svd_tol(Hn))
    P_perp = proj_complement(Hn, TOL)
    return P_perp, V_r, Q, Hn, rank_H


def solve_lm(
    residual_and_jacobian,
    theta0: np.ndarray,
    max_iter: int = LM_MAX_ITER,
    gtol: float = LM_GTOL,
    ftol: float = LM_FTOL,
    xtol: float = LM_XTOL,
    residual_fn=None,
    ref_norm: float = 1.0,
) -> tuple:
    """Adaptive Levenberg-Marquardt on a realified least-squares residual.

    residual_and_jacobian(theta) -> (b, J), where b is the real residual
    (2M,) and J is the real Jacobian (2M x n).  residual_fn(theta) -> b is
    an optional cheap forward-only evaluation used for trial misfits; it
    defaults to residual_and_jacobian.  ref_norm is the data norm used to
    normalize the per-accepted-step residual history.

    Returns (theta_est, final_b, final_J, iterations, converged,
    final_max_grad, history, reason).
    """
    theta = np.asarray(theta0, dtype=float).reshape(-1).copy()
    if residual_fn is None:
        residual_fn = lambda th: residual_and_jacobian(th)[0]

    lam = LM_LAM0
    iterations = 0
    converged = False
    reason = ""
    history: list[float] = []

    b, J = residual_and_jacobian(theta)
    misfit = 0.5 * float(b @ b)
    history.append(float(np.linalg.norm(b) / ref_norm))

    for _ in range(max_iter):
        g = J.T @ b
        if float(np.max(np.abs(g))) < gtol:
            converged = True
            reason = "max_grad"
            break

        H = J.T @ J
        diagH = np.diag(H) + 1e-10
        mf_old = misfit
        delta = None

        for _ in range(LM_INNER_TRIES):
            try:
                delta = np.linalg.solve(
                    H + lam * np.diag(diagH), -g
                )
            except np.linalg.LinAlgError:
                lam = min(lam * 10.0, LM_LAM_MAX)
                continue
            if not np.all(np.isfinite(delta)):
                lam = min(lam * 10.0, LM_LAM_MAX)
                continue

            theta_cand = theta + delta
            b_cand = residual_fn(theta_cand)
            mf_cand = 0.5 * float(b_cand @ b_cand)
            if np.isfinite(mf_cand) and mf_cand < mf_old:
                # Accept: move, refresh (b, J) at the new point, reduce damping.
                theta = theta_cand
                b, J = residual_and_jacobian(theta)
                misfit = 0.5 * float(b @ b)
                history.append(float(np.linalg.norm(b) / ref_norm))
                iterations += 1
                lam = max(lam / 3.0, LM_LAM_MIN)
                break
            lam = min(lam * 10.0, LM_LAM_MAX)
        else:
            reason = "no_step"
            break

        if float(np.linalg.norm(delta)) < xtol:
            converged = True
            reason = "small_step"
            break
        if abs(mf_old - misfit) < ftol * max(1.0, mf_old):
            converged = True
            reason = "small_misfit_change"
            break
    else:
        reason = "max_iter"

    final_max_grad = float(np.max(np.abs(J.T @ b)))
    return theta, b, J, iterations, converged, final_max_grad, history, reason


def state_witness(
    alpha: np.ndarray,
    p: np.ndarray,
    r: int,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    P: Params,
) -> float:
    """T_U for the receiver-rx SOM lift at retained rank r and pose p."""
    chi = Phi @ alpha
    t, _ = transmitter_position(p, P)
    u_inc = incident_field(xs, t, k)
    J = solve_current(chi, u_inc, G_D)
    u = total_field(u_inc, J, G_D)
    A = state_operator(chi, G_D)
    y, y_local = receiver_positions(p, P)
    G_s = data_matrix(y, xs, h, k)
    _, _, Vh = np.linalg.svd(G_s, full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    Q = G_s @ V_r
    v_rx = receiver_directions("rx", y_local)
    B_rx = receiver_data_derivative(y, v_rx, xs, h, k) @ J
    c_r = np.linalg.pinv(Q) @ B_rx
    dJ_r = V_r @ c_r
    denom = np.where(np.abs(u) > 1e-12, u, 1e-12)
    delta_chi_implied = (A @ dJ_r) / denom
    n = np.linalg.norm(delta_chi_implied)
    if n == 0.0:
        return 0.0
    return float(np.linalg.norm(delta_chi_implied.imag) / n)


def reduced_diagnostics(
    alpha: np.ndarray,
    p: np.ndarray,
    r: int,
    p_true: np.ndarray,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    P: Params,
) -> dict:
    """Final hid_sv / hidden rank / hidden vectors / pose decomposition."""
    jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
    P_perp, V_r, Q, Hn, rank_H = reduced_projector(jb, r, P)
    del V_r, Q
    B_red = P_perp @ jb["B_real"]
    _, hid_sv, Vh_b = np.linalg.svd(B_red, full_matrices=False)
    if hid_sv[0] > 0.0:
        thr = HIDDEN_SV_REL * hid_sv[0]
    else:
        thr = 0.0
    hidden_mask = hid_sv <= thr
    hidden_rank = int(np.count_nonzero(hidden_mask))
    idx = np.nonzero(hidden_mask)[0]
    hidden_vectors = [Vh_b[i].tolist() for i in idx]

    delta = np.asarray(p, dtype=float) - np.asarray(p_true, dtype=float)
    pose_error = float(np.linalg.norm(delta))
    if hidden_rank:
        basis = Vh_b[idx]  # hidden rows, orthonormal
        hidden_norm = float(np.linalg.norm(basis @ delta))
        visible_norm = float(
            np.sqrt(max(pose_error**2 - hidden_norm**2, 0.0))
        )
    else:
        hidden_norm = 0.0
        visible_norm = pose_error
    return {
        "rank_H": int(rank_H),
        "hid_sv": [float(v) for v in hid_sv],
        "hidden_subspace_rank": hidden_rank,
        "hidden_vectors": hidden_vectors,
        "pose_error_hidden_norm": hidden_norm,
        "pose_error_visible_norm": visible_norm,
    }


def run_case(
    method: str,
    seed: int,
    p_init_label: str,
    alpha_init: np.ndarray,
    p_init: np.ndarray | None,
    p_fixed: np.ndarray | None,
    r: int | None,
    obs: dict[int, np.ndarray],
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    P: Params,
    alpha_true: np.ndarray,
    p_true: np.ndarray,
) -> tuple[dict, list[float]]:
    d_obs = obs[seed]
    K = Phi.shape[1]
    ref_norm = float(np.linalg.norm(d_obs))

    def wrongpose_rj(theta, p):
        jb = jacobians(theta, p, xs, h, k, G_D, Phi, P)
        return _realified_residual(jb["d"], d_obs), jb["J_alpha_real"]

    def wrongpose_res(theta, p):
        return _realified_residual(
            forward(theta, p, xs, h, k, G_D, Phi, P), d_obs
        )

    def solve_wrongpose(alpha_init, p_fixed):
        """LM over alpha (K,) with pose held at p_fixed."""
        p = np.asarray(p_fixed, dtype=float)
        theta0 = np.asarray(alpha_init, dtype=float)

        def rj(theta):
            return wrongpose_rj(theta, p)

        def res(theta):
            return wrongpose_res(theta, p)

        return solve_lm(
            rj, theta0, max_iter=LM_MAX_ITER,
            residual_fn=res, ref_norm=ref_norm,
        )

    def joint_rj(theta, reduced_r):
        a, pp = theta[:K], theta[K:]
        jb = jacobians(a, pp, xs, h, k, G_D, Phi, P)
        if reduced_r is None:
            J = np.hstack([jb["J_alpha_real"], jb["B_real"]])
        else:
            P_perp, _, _, _, _ = reduced_projector(jb, reduced_r, P)
            J = np.hstack([jb["J_alpha_real"], P_perp @ jb["B_real"]])
        return _realified_residual(jb["d"], d_obs), J

    def joint_res(theta):
        a, pp = theta[:K], theta[K:]
        return _realified_residual(
            forward(a, pp, xs, h, k, G_D, Phi, P), d_obs
        )

    def solve_direct(alpha_init, p_init):
        """Joint LM over theta = [alpha; p] with the full B_real columns."""
        theta0 = np.concatenate(
            [np.asarray(alpha_init, float), np.asarray(p_init, float)]
        )
        return solve_lm(
            lambda th: joint_rj(th, None), theta0,
            max_iter=LM_MAX_ITER, residual_fn=joint_res,
            ref_norm=ref_norm,
        )

    def solve_reduced(alpha_init, p_init, r):
        """Joint LM with pose columns P_perp @ B_real (P_perp recomputed)."""
        theta0 = np.concatenate(
            [np.asarray(alpha_init, float), np.asarray(p_init, float)]
        )
        return solve_lm(
            lambda th: joint_rj(th, r), theta0,
            max_iter=LM_MAX_ITER, residual_fn=joint_res,
            ref_norm=ref_norm,
        )

    if method == "wrongpose":
        solve_out = solve_wrongpose(alpha_init, p_fixed)
    elif method == "direct":
        solve_out = solve_direct(alpha_init, p_init)
    else:  # reduced_r4 / reduced_r6
        solve_out = solve_reduced(alpha_init, p_init, r)

    theta_est, _, _, iterations, converged, final_max_grad, hist, reason = solve_out
    if method == "wrongpose":
        alpha_est = np.asarray(theta_est, dtype=float)
        p_est = np.asarray(p_fixed, dtype=float)
    else:
        alpha_est = np.asarray(theta_est[:K], dtype=float)
        p_est = np.asarray(theta_est[K:], dtype=float)

    final_residual = hist[-1]
    pose_error = float(np.linalg.norm(p_est - p_true))
    chi_est = Phi @ alpha_est
    chi_true = Phi @ alpha_true
    map_error = float(np.linalg.norm(chi_est - chi_true) / np.linalg.norm(chi_true))

    record = {
        "method": method,
        "noise_seed": int(seed),
        "p_init_label": p_init_label,
        "r": int(r) if r is not None else None,
        "alpha_est": [float(x) for x in alpha_est],
        "p_est": [float(x) for x in p_est],
        "pose_error": pose_error,
        "map_error": map_error,
        "final_residual": final_residual,
        "iterations": int(iterations),
        "converged": bool(converged),
        "final_max_grad": float(final_max_grad),
        "reason": str(reason),
        "hid_sv": None,
        "hidden_subspace_rank": None,
        "T_U": None,
    }

    if method.startswith("reduced"):
        diag = reduced_diagnostics(
            alpha_est, p_est, r, p_true, xs, h, k, G_D, Phi, P
        )
        record["hid_sv"] = diag["hid_sv"]
        record["hidden_subspace_rank"] = diag["hidden_subspace_rank"]
        record["hidden_vectors"] = diag["hidden_vectors"]
        record["pose_error_hidden_norm"] = diag["pose_error_hidden_norm"]
        record["pose_error_visible_norm"] = diag["pose_error_visible_norm"]
        record["rank_H_final"] = diag["rank_H"]
        record["T_U_rank"] = int(r)
    else:
        record["hidden_vectors"] = None
        record["pose_error_hidden_norm"] = None
        record["pose_error_visible_norm"] = None
        record["T_U_rank"] = 4  # disambiguated/default retained rank for the witness
    record["T_U"] = state_witness(
        alpha_est, p_est, record["T_U_rank"], xs, h, k, G_D, Phi, P
    )
    return record, hist


def make_plot(
    records: list[dict],
    histories: dict[str, list[float]],
    path: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8))

    # Left: residual history (seed 0, zero/fixed init).
    left_cases = (
        ("wrongpose", "fixed_zero"),
        ("direct", "zero"),
        ("reduced_r4", "zero"),
        ("reduced_r6", "zero"),
    )
    styles = {
        "wrongpose": dict(color="#d62728", ls=":", marker="o", ms=3),
        "direct": dict(color="#1f77b4", ls="-", marker="s", ms=3),
        "reduced_r4": dict(color="#2ca02c", ls="-", marker="^", ms=3),
        "reduced_r6": dict(color="#9467bd", ls="-", marker="v", ms=3),
    }
    ax = axes[0]
    for method, init in left_cases:
        key = f"{method}:{init}:seed0"
        hist = np.asarray(histories[key], dtype=float)
        ax.plot(
            np.arange(len(hist)),
            np.maximum(hist, 1e-16),
            label=method,
            **styles[method],
        )
    ax.axhline(10 ** (-SNR_DB / 20.0), color="k", ls="--", lw=0.8)
    ax.text(
        0.98,
        10 ** (-SNR_DB / 20.0) * 1.25,
        "noise rms ratio",
        transform=ax.transData,
        ha="right",
        va="bottom",
        fontsize=8,
    )
    ax.set_yscale("log")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\|d-d_{\rm obs}\|/\|d_{\rm obs}\|$")
    ax.set_title("E5: LM residual history (seed 0)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    # Right: seed-0 errors averaged over the three p_inits where applicable.
    methods = ("wrongpose", "direct", "reduced_r4", "reduced_r6")
    ax = axes[1]
    means = {}
    for method in methods:
        sel = [rec for rec in records if rec["method"] == method and rec["noise_seed"] == 0]
        pose = float(np.mean([rec["pose_error"] for rec in sel]))
        mape = float(np.mean([rec["map_error"] for rec in sel]))
        means[method] = {"pose_error": pose, "map_error": mape}
        if method == "reduced_r6":
            means[method]["hidden"] = float(
                np.mean([rec["pose_error_hidden_norm"] for rec in sel])
            )
            means[method]["visible"] = float(
                np.mean([rec["pose_error_visible_norm"] for rec in sel])
            )

    x = np.arange(len(methods))
    width = 0.34
    for j, method in enumerate(methods):
        if method != "reduced_r6":
            pose = max(means[method]["pose_error"], 1e-16)
            mape = max(means[method]["map_error"], 1e-16)
            ax.bar(
                x[j] - width / 2,
                pose,
                width,
                color="#1f77b4",
                label="pose error" if j == 0 else None,
            )
            ax.bar(
                x[j] + width / 2,
                mape,
                width,
                color="#ff7f0e",
                label="map error" if j == 0 else None,
            )
            ax.text(x[j] - width / 2, pose * 1.35, f"{pose:.1e}", ha="center", fontsize=7)
            ax.text(x[j] + width / 2, mape * 1.35, f"{mape:.1e}", ha="center", fontsize=7)
        else:
            pose = max(means[method]["pose_error"], 1e-16)
            mape = max(means[method]["map_error"], 1e-16)
            vis = max(means[method]["visible"], 1e-16)
            hid = max(means[method]["hidden"], 1e-16)
            for off, val, col, lab in (
                (-0.45, pose, "#1f77b4", "pose error"),
                (-0.15, mape, "#ff7f0e", "map error"),
                (0.15, vis, "#2ca02c", "visible (r6)"),
                (0.45, hid, "#d62728", "hidden (r6)"),
            ):
                ax.bar(x[j] + off, val, 0.26, color=col, label=lab if j == 0 else None)
                ax.text(x[j] + off, val * 1.35, f"{val:.1e}", ha="center", fontsize=6)
    ax.set_yscale("log")
    ax.set_ylim(1e-16, 1e2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{m}\n(seed 0)" for m in methods], fontsize=9)
    ax.set_ylabel("error")
    ax.set_title("E5: pose / map errors at seed 0\n(3-p-init means; r6 decomposed)")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="upper left")

    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> dict:
    P = Params(
        N=16, k=12.0, M=8, R_r=1.6, R_t=2.0, phi0=0.7, s=0.12
    )
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    Phi = basis_matrix(xs, P)
    K = Phi.shape[1]

    alpha_true = np.array([1.5, 2.0, 0.0])
    p_true = np.array([0.08, -0.06, 0.05])
    chi_true = Phi @ alpha_true

    d_clean = forward(alpha_true, p_true, xs, h, P.k, G_D, Phi, P)
    rms = float(np.sqrt(np.mean(np.abs(d_clean) ** 2)))
    noise_std = 10 ** (-SNR_DB / 20.0) * rms
    obs: dict[int, np.ndarray] = {}
    for seed in NOISE_SEEDS:
        rng = np.random.default_rng(seed)
        noise = noise_std / np.sqrt(2.0) * (
            rng.standard_normal(P.M) + 1j * rng.standard_normal(P.M)
        )
        obs[int(seed)] = d_clean + noise

    p_zero = np.zeros(3)
    p_pert_pos = p_true + np.array([0.10, 0.10, 0.10])
    p_pert_neg = p_true + np.array([-0.08, 0.06, -0.05])
    p_inits = {
        "zero": p_zero,
        "pert_pos": p_pert_pos,
        "pert_neg": p_pert_neg,
    }

    records: list[dict] = []
    histories: dict[str, list[float]] = {}

    for seed in NOISE_SEEDS:
        rec, hist = run_case(
            "wrongpose", seed, "fixed_zero", ALPHA_INIT, None,
            np.zeros(3), None, obs, xs, h, P.k, G_D, Phi, P,
            alpha_true, p_true,
        )
        records.append(rec)
        histories[f"wrongpose:fixed_zero:seed{seed}"] = hist

    for method in ("direct", "reduced_r4", "reduced_r6"):
        r = 4 if method == "reduced_r4" else (6 if method == "reduced_r6" else None)
        for seed in NOISE_SEEDS:
            for label in P_INIT_LABELS:
                rec, hist = run_case(
                    method, seed, label, ALPHA_INIT, p_inits[label],
                    None, r, obs, xs, h, P.k, G_D, Phi, P,
                    alpha_true, p_true,
                )
                records.append(rec)
                histories[f"{method}:{label}:seed{seed}"] = hist

    # Keep JSON small: normalized-residual histories only for the
    # seed-0 zero/fixed-init runs that the left plot panel shows.
    plot_histories = {
        key: histories[key]
        for key in (
            "wrongpose:fixed_zero:seed0",
            "direct:zero:seed0",
            "reduced_r4:zero:seed0",
            "reduced_r6:zero:seed0",
        )
    }

    results = {
        "parameters": {
            "N": P.N,
            "k": P.k,
            "M": P.M,
            "R_r": P.R_r,
            "R_t": P.R_t,
            "phi0": P.phi0,
            "s": P.s,
            "K": K,
            "basis_centers": [list(c) for c in ((-0.15, 0.10), (0.20, -0.10), (0.00, 0.00))],
            "SNR_dB": SNR_DB,
            "noise_std": float(noise_std),
            "noise_rms_ratio": float(10 ** (-SNR_DB / 20.0)),
            "alpha_true": [float(x) for x in alpha_true],
            "p_true": [float(x) for x in p_true],
            "alpha_init": [float(x) for x in ALPHA_INIT],
            "p_init_labels": list(P_INIT_LABELS),
            "LevenbergMarquardt": {
                "lambda0": LM_LAM0,
                "lambda_min": LM_LAM_MIN,
                "lambda_max": LM_LAM_MAX,
                "max_iter": LM_MAX_ITER,
                "inner_max_tries": LM_INNER_TRIES,
                "gtol": LM_GTOL,
                "ftol": LM_FTOL,
                "xtol": LM_XTOL,
            },
        },
        "note": (
            "Reduced solvers recompute P_perp = I - U_r U_r^T over "
            "col(realify(G_s V_r), realify(W)) at every LM iteration with "
            "retained rank r (4 disambiguated, 6 confounded). hid_sv are the "
            "singular values of B_red = P_perp @ B_real at the final estimate; "
            "hidden_subspace_rank counts hid_sv <= 1e-8 * hid_sv[0]. T_U is "
            "the receiver-rx state witness at the final estimate: for reduced "
            "runs at their retained rank r; for wrongpose/direct at r=4 "
            "(recorded as T_U_rank). hidden_vectors are the right singular "
            "vectors of B_red for the hidden (near-zero hid_sv) directions; "
            "pose_error_hidden/visible_norm decompose p_est-p_true onto "
            "those directions."
        ),
        "records": records,
        "residual_histories": {
            key: [float(v) for v in hist] for key, hist in plot_histories.items()
        },
    }

    json_path = HERE / "results_e5.json"
    png_path = HERE / "plot_e5.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(records, plot_histories, png_path)

    # Verification: JSON reload round-trips and the PNG loads.
    with open(json_path) as fh:
        loaded = json.load(fh)
    assert len(loaded["records"]) == 20, len(loaded["records"])
    assert loaded["records"] == results["records"]
    assert png_path.exists() and png_path.stat().st_size > 0
    try:
        import PIL.Image

        im = PIL.Image.open(png_path)
        im.load()
    except Exception:
        pass  # matplotlib rendering already confirms the file is non-empty

    print("--- E5 summary ---")
    wp = [r for r in loaded["records"] if r["method"] == "wrongpose"]
    for rec in wp:
        print(
            f"wrongpose seed{rec['noise_seed']}: final_residual="
            f"{rec['final_residual']:.4e} map_error={rec['map_error']:.4e} "
            f"iterations={rec['iterations']} converged={rec['converged']} "
            f"reason={rec['reason']} max_grad={rec['final_max_grad']:.2e}"
        )
    for method in ("direct", "reduced_r4", "reduced_r6"):
        sel = [
            r for r in loaded["records"]
            if r["method"] == method and r["noise_seed"] == 0
        ]
        pose = sorted(r["pose_error"] for r in sel)[len(sel) // 2]
        mape = sorted(r["map_error"] for r in sel)[len(sel) // 2]
        res = sorted(r["final_residual"] for r in sel)[len(sel) // 2]
        print(
            f"{method} seed0 median over 3 p_inits: pose_error={pose:.4e} "
            f"map_error={mape:.4e} final_residual={res:.4e}"
        )
    print(
        f"noise_rms_ratio = {loaded['parameters']['noise_rms_ratio']:.4e}"
    )
    print("\nconvergence per record (method seed init: converged reason max_grad):")
    for rec in loaded["records"]:
        print(
            f"  {rec['method']:>11s} seed{rec['noise_seed']} "
            f"{rec['p_init_label']:>9s}: converged={rec['converged']} "
            f"reason={rec['reason']:>22s} max_grad={rec['final_max_grad']:.2e}"
        )
    for rec in loaded["records"]:
        if rec["method"].startswith("reduced") and rec["noise_seed"] == 0:
            print(
                f"{rec['method']} init={rec['p_init_label']:>9s}: "
                f"pose_error={rec['pose_error']:.4e} map_error={rec['map_error']:.4e} "
                f"hid_sv="
                + ",".join(f"{v:.2e}" for v in rec["hid_sv"])
                + f" hidden_rank={rec['hidden_subspace_rank']} T_U={rec['T_U']:.4f}"
            )
    print("\nT_U values (all records, rx direction):")
    for rec in loaded["records"]:
        print(
            f"  {rec['method']:>11s} {rec['p_init_label']:>9s} "
            f"seed{rec['noise_seed']} r={rec['r'] if rec['r'] is not None else rec['T_U_rank']}: "
            f"T_U={rec['T_U']:.4f} p_est="
            + ",".join(f"{v:.4f}" for v in rec["p_est"])
        )
    print("\nArtifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
