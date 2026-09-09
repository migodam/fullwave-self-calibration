"""E7 -- state-consistent reduced self-calibration in the confounded regime.

Single-transmitter co-moving fixture (identical to the E5 confounded case):
N=16, k=12, M=8, R_r=1.6, R_t=2.0, phi0=0.7, s=0.12, K=3 Gaussian contrast
basis, alpha_true=[1.5,2.0,0], alpha_init=[1,1,0], p_true=[0.08,-0.06,0.05],
p_init=0, retained rank r=6, SNR 30 dB.  At (alpha_init, p_init=0) the
linearised visible-pose subspace has n_vis=1 (hidden_rank=2); direct joint
GN recovers pose while reduced_r6 freezes the hidden directions.

The experiment asks whether making the *state* participate algorithmically
can lift that freeze:

  A direct     : theta=[alpha;p], full-physics data, analytic Jacobian.
  B reduced_r6 : theta=[alpha;q], p=p_init+V_vis@q (V_vis fixed r6 subspace),
                 full-physics data, analytic Jacobian.
  C state_c    : theta=[alpha;c_real;c_imag;p], data from retained current
                 G_s(p) V_r(p) c plus lam-weighted state-consistency residual
                 (V_r(p) c vs the physical current J_phys(alpha,p)), solved
                 with a built-in finite-difference Jacobian; lam sweep.
  D tu_penalty : theta=[alpha;p], full-physics data plus the scalar penalty
                 lam_T * state_witness(alpha,p,r=6) appended to the realified
                 data residual; finite-difference Jacobian; lam_T sweep.

Deterministic: only two fresh numpy default_rng draws (seeds 0,1) for noise.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from geom_som_core import (
    Params,
    pixel_grid,
    green_domain_matrix,
    data_matrix,
    solve_current,
    incident_field,
    transmitter_position,
    receiver_positions,
)
from run_e5 import (
    basis_matrix,
    forward,
    jacobians,
    state_witness,
    _realified_residual,
)
from run_e5_final import visible_pose_subspace


HERE = Path(__file__).resolve().parent

K = 3
RETAINED_R = 6
SNR_DB = 30.0
NOISE_RATIO = float(10.0 ** (-SNR_DB / 20.0))
NOISE_SEEDS = (0, 1)

ALPHA_TRUE = np.array([1.5, 2.0, 0.0], dtype=float)
ALPHA_INIT = np.array([1.0, 1.0, 0.0], dtype=float)
P_TRUE = np.array([0.08, -0.06, 0.05], dtype=float)
P_INIT = np.zeros(3)

STATE_C_LAMS = (1e-4, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0)
TU_LAMS = (1e-3, 1e-2, 1e-1, 1.0)

LS_KWARGS = dict(
    method="trf",
    x_scale="jac",
    max_nfev=2000,
    ftol=1e-10,
    xtol=1e-10,
    gtol=1e-10,
)
J_FLOOR_REL = 1e-9  # guard for ||J_phys|| -> 0 in the state residual


def realify(z: np.ndarray) -> np.ndarray:
    """realify(z) = [Re z; Im z] preserving the 2-norm."""
    return np.concatenate([z.real, z.imag])


def physical_current(
    alpha: np.ndarray, p: np.ndarray, xs, h, k, G_D, Phi, P
) -> np.ndarray:
    chi = Phi @ alpha
    t, _ = transmitter_position(p, P)
    u_inc = incident_field(xs, t, k)
    return solve_current(chi, u_inc, G_D)


def retained_basis(
    p: np.ndarray,
    r: int,
    xs,
    h,
    k,
    P,
    V_ref: np.ndarray | None = None,
) -> np.ndarray:
    """First r right singular vectors of G_s(p) as Np x r complex columns.

    A complex SVD frame is unique only up to per-column unit phases, and the
    hard r-cut can sit on a (near-)degenerate singular-value pair.  When a
    reference frame V_ref is supplied the columns are Procrustes-rotated to be
    closest to V_ref, giving a deterministic gauge that is continuous wherever
    the retained subspace itself is well defined.  The span (hence data model
    and the best c) is unchanged.
    """
    y, _ = receiver_positions(p, P)
    G_s = data_matrix(y, xs, h, k)
    _, _, Vh = np.linalg.svd(G_s, full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    if V_ref is not None:
        O = V_r.conj().T @ V_ref  # r x r overlap with the reference frame
        Uo, _, Vho = np.linalg.svd(O, full_matrices=False)
        Q = Uo @ Vho  # unitary Procrustes rotation (no conjugation needed
        # because Vho holds the right vectors for O = Uo Sigma Vho^H; the
        # optimal rotation is Uo @ (Vho^H)^H = Uo @ Vho).
        V_r = V_r @ Q
    return V_r


def run_baseline_method(
    method: str,
    seed: int,
    obs: dict[int, np.ndarray],
    subspace: dict,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
) -> dict:
    """A (direct) and B (reduced_r6) with analytic Jacobians (E5-final style)."""
    d_obs = obs[seed]
    V_vis = subspace["V_vis"]  # 3 x n_vis at (alpha_init, p_init), r=6

    def residual(theta: np.ndarray) -> np.ndarray:
        if method == "direct":
            alpha = theta[:K]
            p = theta[K:]
        else:
            alpha = theta[:K]
            q = theta[K:]
            p = P_INIT + V_vis @ q
        return _realified_residual(
            forward(alpha, p, xs, h, k, G_D, Phi, P), d_obs
        )

    def jacobian(theta: np.ndarray) -> np.ndarray:
        if method == "direct":
            alpha = theta[:K]
            p = theta[K:]
            jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
            return np.hstack([jb["J_alpha_real"], jb["B_real"]])
        alpha = theta[:K]
        q = theta[K:]
        p = P_INIT + V_vis @ q
        jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
        return np.hstack(
            [jb["J_alpha_real"], jb["B_real"] @ V_vis]
        )

    if method == "direct":
        x0 = np.concatenate([ALPHA_INIT, P_INIT])
    else:
        x0 = np.concatenate([ALPHA_INIT, np.zeros(int(subspace["n_vis"]))])

    result = least_squares(residual, x0, jac=jacobian, **LS_KWARGS)
    alpha_est = np.asarray(result.x[:K], dtype=float)
    if method == "direct":
        p_est = np.asarray(result.x[K:], dtype=float)
    else:
        q_est = np.asarray(result.x[K:], dtype=float)
        p_est = P_INIT + V_vis @ q_est
    return finalize_record(
        method=method,
        lam=None,
        seed=seed,
        alpha_est=alpha_est,
        p_est=p_est,
        result=result,
        xs=xs,
        h=h,
        k=k,
        G_D=G_D,
        Phi=Phi,
        P=P,
        obs=obs,
        subspace=subspace,
        extra=None,
    )


def make_state_c_problem(
    lam: float,
    d_obs: np.ndarray,
    j_ref: float,
    V_ref: np.ndarray,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
) -> tuple:
    """Residual for the state-consistent reduced variant C at weight lam."""

    def residual(theta: np.ndarray) -> np.ndarray:
        alpha = np.asarray(theta[:K], dtype=float)
        c = (
            np.asarray(theta[K : K + RETAINED_R], dtype=float)
            + 1j * np.asarray(theta[K + RETAINED_R : K + 2 * RETAINED_R], dtype=float)
        )
        p = np.asarray(theta[K + 2 * RETAINED_R :], dtype=float)
        chi = Phi @ alpha
        t, _ = transmitter_position(p, P)
        J_phys = solve_current(chi, incident_field(xs, t, k), G_D)
        V_r = retained_basis(p, RETAINED_R, xs, h, k, P, V_ref)
        cur = V_r @ c
        y, _ = receiver_positions(p, P)
        G_s = data_matrix(y, xs, h, k)
        r_data = realify(G_s @ cur - d_obs)
        denom = max(float(np.linalg.norm(J_phys)), J_FLOOR_REL * j_ref)
        r_state = lam * realify(cur - J_phys) / denom
        return np.concatenate([r_data, r_state])

    return residual


def run_state_c(
    lam: float,
    seed: int,
    obs: dict[int, np.ndarray],
    j_ref: float,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
    subspace: dict,
) -> dict:
    d_obs = obs[seed]
    V_ref = retained_basis(P_INIT, RETAINED_R, xs, h, k, P)
    residual = make_state_c_problem(
        lam, d_obs, j_ref, V_ref, xs, h, k, G_D, Phi, P
    )
    c0 = V_ref.conj().T @ (
        physical_current(ALPHA_INIT, P_INIT, xs, h, k, G_D, Phi, P)
    )
    x0 = np.concatenate(
        [ALPHA_INIT, c0.real, c0.imag, P_INIT]
    )
    try:
        result = least_squares(residual, x0, **LS_KWARGS)
        alpha_est = np.asarray(result.x[:K], dtype=float)
        c_est = (
            np.asarray(result.x[K : K + RETAINED_R], dtype=float)
            + 1j
            * np.asarray(result.x[K + RETAINED_R : K + 2 * RETAINED_R], dtype=float)
        )
        p_est = np.asarray(result.x[K + 2 * RETAINED_R :], dtype=float)
        # State/diagnostics at the solution.
        V_r = retained_basis(p_est, RETAINED_R, xs, h, k, P, V_ref)
        y, _ = receiver_positions(p_est, P)
        G_s = data_matrix(y, xs, h, k)
        cur = V_r @ c_est
        J_phys = physical_current(alpha_est, p_est, xs, h, k, G_D, Phi, P)
        denom = max(float(np.linalg.norm(J_phys)), J_FLOOR_REL * j_ref)
        extra = {
            "retained_data_residual": float(
                np.linalg.norm(G_s @ cur - d_obs) / np.linalg.norm(d_obs)
            ),
            "state_resid": float(np.linalg.norm(cur - J_phys) / denom),
            "c_norm": float(np.linalg.norm(c_est)),
        }
    except Exception as exc:  # noqa: BLE001 - a lam failure is a record
        return failed_record(
            method="state_c",
            lam=lam,
            seed=seed,
            alpha_est=ALPHA_INIT,
            p_est=P_INIT,
            xs=xs,
            h=h,
            k=k,
            G_D=G_D,
            Phi=Phi,
            P=P,
            obs=obs,
            subspace=subspace,
            reason=f"exception: {type(exc).__name__}: {exc}",
        )
    return finalize_record(
        method="state_c",
        lam=lam,
        seed=seed,
        alpha_est=alpha_est,
        p_est=p_est,
        result=result,
        xs=xs,
        h=h,
        k=k,
        G_D=G_D,
        Phi=Phi,
        P=P,
        obs=obs,
        subspace=subspace,
        extra=extra,
    )


def run_tu_penalty(
    lam_t: float,
    seed: int,
    obs: dict[int, np.ndarray],
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
    subspace: dict,
) -> dict:
    d_obs = obs[seed]

    def residual(theta: np.ndarray) -> np.ndarray:
        alpha = theta[:K]
        p = theta[K:]
        r_data = _realified_residual(
            forward(alpha, p, xs, h, k, G_D, Phi, P), d_obs
        )
        w = state_witness(
            alpha, p, RETAINED_R, xs, h, k, G_D, Phi, P
        )
        return np.concatenate([r_data, [lam_t * w]])

    x0 = np.concatenate([ALPHA_INIT, P_INIT])
    try:
        result = least_squares(residual, x0, **LS_KWARGS)
        alpha_est = np.asarray(result.x[:K], dtype=float)
        p_est = np.asarray(result.x[K:], dtype=float)
        extra = None
    except Exception as exc:  # noqa: BLE001
        return failed_record(
            method="tu_penalty",
            lam=lam_t,
            seed=seed,
            alpha_est=ALPHA_INIT,
            p_est=P_INIT,
            xs=xs,
            h=h,
            k=k,
            G_D=G_D,
            Phi=Phi,
            P=P,
            obs=obs,
            subspace=subspace,
            reason=f"exception: {type(exc).__name__}: {exc}",
        )
    return finalize_record(
        method="tu_penalty",
        lam=lam_t,
        seed=seed,
        alpha_est=alpha_est,
        p_est=p_est,
        result=result,
        xs=xs,
        h=h,
        k=k,
        G_D=G_D,
        Phi=Phi,
        P=P,
        obs=obs,
        subspace=subspace,
        extra=extra,
    )


def decompose_pose(
    p_est: np.ndarray, subspace: dict
) -> tuple[float, float, float]:
    delta = np.asarray(p_est, dtype=float) - P_TRUE
    V_vis = subspace["V_vis"]
    vis = V_vis @ (V_vis.T @ delta)
    return (
        float(np.linalg.norm(delta)),
        float(np.linalg.norm(vis)),
        float(np.linalg.norm(delta - vis)),
    )


def finalize_record(
    method: str,
    lam: float | None,
    seed: int,
    alpha_est: np.ndarray,
    p_est: np.ndarray,
    result,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
    obs: dict[int, np.ndarray],
    subspace: dict,
    extra: dict | None,
) -> dict:
    d_obs = obs[seed]
    ref_norm = float(np.linalg.norm(d_obs))
    pose_error, pose_vis, pose_hid = decompose_pose(p_est, subspace)
    chi_true = Phi @ ALPHA_TRUE
    chi_est = Phi @ alpha_est
    map_error = float(
        np.linalg.norm(chi_est - chi_true) / np.linalg.norm(chi_true)
    )
    full_res = float(
        np.linalg.norm(
            forward(alpha_est, p_est, xs, h, k, G_D, Phi, P) - d_obs
        )
        / ref_norm
    )
    t_u = state_witness(
        alpha_est, p_est, RETAINED_R, xs, h, k, G_D, Phi, P
    )
    rec = {
        "method": method,
        "lam": lam if lam is not None else None,
        "noise_seed": int(seed),
        "success": bool(result.success),
        "status": int(result.status),
        "optimality": float(result.optimality),
        "nfev": int(result.nfev),
        "alpha_est": [float(v) for v in alpha_est],
        "p_est": [float(v) for v in p_est],
        "pose_error": pose_error,
        "pose_error_visible": pose_vis,
        "pose_error_hidden": pose_hid,
        "map_error": map_error,
        "final_full_data_residual": full_res,
        "T_U": t_u,
    }
    p_shift = float(np.linalg.norm(np.asarray(p_est) - P_INIT))
    rec["p_shift_from_init"] = p_shift
    rec["p_frozen_at_init"] = bool(p_shift < 1e-6)
    if extra is not None:
        rec.update(extra)
    return rec


def failed_record(
    method: str,
    lam: float | None,
    seed: int,
    alpha_est: np.ndarray,
    p_est: np.ndarray,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
    obs: dict[int, np.ndarray],
    subspace: dict,
    reason: str,
) -> dict:
    pose_error, pose_vis, pose_hid = decompose_pose(p_est, subspace)
    chi_true = Phi @ ALPHA_TRUE
    map_error = float(
        np.linalg.norm(Phi @ alpha_est - chi_true) / np.linalg.norm(chi_true)
    )
    return {
        "method": method,
        "lam": lam,
        "noise_seed": int(seed),
        "success": False,
        "status": -99,
        "optimality": np.inf,
        "nfev": 0,
        "alpha_est": [float(v) for v in alpha_est],
        "p_est": [float(v) for v in p_est],
        "pose_error": pose_error,
        "pose_error_visible": pose_vis,
        "pose_error_hidden": pose_hid,
        "map_error": map_error,
        "final_full_data_residual": None,
        "T_U": None,
        "failure_reason": reason,
    }


def summarize(records: list[dict]) -> list[dict]:
    """Median over noise seeds per (method, lam)."""
    keys: dict[tuple, list[dict]] = {}
    for rec in records:
        key = (rec["method"], rec["lam"])
        keys.setdefault(key, []).append(rec)
    out = []
    for (method, lam), group in keys.items():
        ok = [r for r in group if r.get("final_full_data_residual") is not None]
        if not ok:
            out.append(
                {
                    "method": method,
                    "lam": lam,
                    "n_seeds": len(group),
                    "all_failed": True,
                }
            )
            continue
        base = {
            "method": method,
            "lam": lam,
            "n_seeds": len(group),
            "n_success": int(sum(bool(r["success"]) for r in group)),
            "n_strict": int(
                sum(
                    bool(r["success"] and r["optimality"] < 1e-8)
                    for r in group
                )
            ),
        }
        fields = (
            "pose_error",
            "pose_error_visible",
            "pose_error_hidden",
            "map_error",
            "final_full_data_residual",
            "T_U",
        )
        for f in fields:
            vals = [r[f] for r in ok if r.get(f) is not None]
            base[f + "_median"] = float(np.median(vals)) if vals else None
        if method == "state_c":
            for f in (
                "retained_data_residual",
                "state_resid",
                "c_norm",
            ):
                vals = [r[f] for r in ok if r.get(f) is not None]
                base[f + "_median"] = float(np.median(vals)) if vals else None
        out.append(base)
    return out


def structural_diagnostics(
    obs: dict[int, np.ndarray],
    j_ref: float,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
) -> dict:
    """Evidence for the finite-difference pinning of state_c at p_init=0."""
    y0, _ = receiver_positions(P_INIT, P)
    G_s0 = data_matrix(y0, xs, h, k)
    s = np.linalg.svd(G_s0, compute_uv=False)
    V_ref = retained_basis(P_INIT, RETAINED_R, xs, h, k, P)
    # Numerical Jacobian column norms of the state_c residual at x0 (seed 0).
    d_obs = obs[0]
    residual = make_state_c_problem(
        1e-2, d_obs, j_ref, V_ref, xs, h, k, G_D, Phi, P
    )
    c0 = V_ref.conj().T @ physical_current(
        ALPHA_INIT, P_INIT, xs, h, k, G_D, Phi, P
    )
    x0 = np.concatenate([ALPHA_INIT, c0.real, c0.imag, P_INIT])
    b0 = residual(x0)
    n = len(x0)
    J = np.empty((len(b0), n))
    step = 1.5e-8  # scipy default abs step for an x0 entry equal to zero
    for j in range(n):
        xp = x0.copy()
        xm = x0.copy()
        xp[j] += step
        xm[j] -= step
        J[:, j] = (residual(xp) - residual(xm)) / (2.0 * step)
    col_norm = np.linalg.norm(J, axis=0)
    p_idx = slice(K + 2 * RETAINED_R, n)
    return {
        "G_s_singular_values_at_p_init": [float(v) for v in s],
        "s6_minus_s7_at_p_init": float(s[5] - s[6]),
        "note": (
            "s6 and s7 of G_s(p_init) are degenerate to machine precision and "
            "sit exactly on the hard r=6 retained cut, so the 'first r right "
            "singular vectors' frame has a conical singularity at p=0.  "
            "Finite-difference Jacobian columns for p then measure a "
            "subspace-jump of an ill-defined cut rather than a physical "
            "derivative; TRF accepts only ~1e-11 steps and reports xtol with "
            "a large optimality (status 2-4), freezing state_c at p_init for "
            "every lambda.  Diagnostics below were computed at the seed-0 "
            "state_c starting point with lambda=1e-2 and step 1.5e-8."
        ),
        "numeric_jacobian_column_norms_at_x0": {
            "alpha_max": float(np.max(col_norm[:K])),
            "c_max": float(np.max(col_norm[K : K + 2 * RETAINED_R])),
            "p_columns": [float(v) for v in col_norm[p_idx]],
            "p_max": float(np.max(col_norm[p_idx])),
        },
    }


def make_plot(
    records: list[dict],
    summaries: list[dict],
    baselines: dict,
    png_path: Path,
) -> None:
    def meds(method: str, key_lam: str = "lam") -> dict:
        sel = [s for s in summaries if s["method"] == method]
        return {
            s[key_lam]: s
            for s in sel
            if s.get("pose_error_median") is not None
        }

    sc = meds("state_c")
    tu = meds("tu_penalty", "lam")
    direct = baselines["direct"]
    reduced = baselines["reduced_r6"]

    def x_for(lam):
        return 3e-5 if lam == 0.0 else lam

    fig, axes = plt.subplots(1, 3, figsize=(19.5, 5.8))

    # Panel 1: state_c pose error vs lam (visible/hidden split).
    ax = axes[0]
    lams = sorted(sc)
    lams_nonzero = [lam for lam in lams if lam is not None and lam != 0.0]
    xvals = [x_for(lam) for lam in lams if lam is not None]
    vals = {f: [] for f in ("pose_error", "pose_error_visible", "pose_error_hidden")}
    keep_lams = []
    for lam, s in sc.items():
        if lam is None:
            continue
        keep_lams.append(lam)
        for f in vals:
            vals[f].append(s.get(f + "_median"))
    styles = {
        "pose_error": dict(color="#1f77b4", marker="o", label="pose error"),
        "pose_error_visible": dict(
            color="#2ca02c", marker="s", label="visible (r6 decomp)"
        ),
        "pose_error_hidden": dict(
            color="#d62728", marker="^", label="hidden (r6 decomp)"
        ),
    }
    for f, sty in styles.items():
        keep = [
            (xv, yv)
            for xv, yv in zip([x_for(l) for l in keep_lams], vals[f])
            if yv is not None
        ]
        if keep:
            ax.semilogx(
                [a for a, _ in keep],
                [b for _, b in keep],
                ls="-",
                ms=5,
                **sty,
            )
    ax.axhline(
        direct["pose_error_median"],
        color="#1f77b4",
        ls="--",
        lw=1.0,
        label=f"direct {direct['pose_error_median']:.2e}",
    )
    ax.axhline(
        reduced["pose_error_median"],
        color="#9467bd",
        ls=":",
        lw=1.0,
        label=f"reduced_r6 {reduced['pose_error_median']:.2e}",
    )
    ax.axhline(
        reduced["pose_error_hidden_median"],
        color="#d62728",
        ls=":",
        lw=0.8,
        label=f"r6 hidden {reduced['pose_error_hidden_median']:.2e}",
    )
    failed_sc = [lam for lam in STATE_C_LAMS if lam not in sc]
    fail_txt = f"failed lam={failed_sc}" if failed_sc else ""
    ax.set_title(
        "E7: state_c pose error vs lam\n"
        "(medians over seeds; state_c never left p_init=0 - finite-difference "
        "Jacobians pinned by the s6=s7 degeneracy at the r=6 cut)\n"
        + fail_txt,
        fontsize=9,
    )
    ax.set_xlabel(r"state weight $\lambda$")
    ax.set_ylabel("pose error")
    ax.set_ylim(1e-4, 1.0)
    ax.set_xticks([3e-5] + lams_nonzero)
    ax.set_xticklabels(["0"] + [f"{v:.0e}" for v in lams_nonzero], rotation=45)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    # Panel 2: state_c final full-data residual vs lam.
    ax = axes[1]
    keep = [
        (x_for(lam), s["final_full_data_residual_median"])
        for lam, s in sc.items()
        if lam is not None
        and s.get("final_full_data_residual_median") is not None
    ]
    if keep:
        ax.semilogx(
            [a for a, _ in keep],
            [b for _, b in keep],
            color="#1f77b4",
            marker="o",
            ls="-",
            ms=5,
            label="state_c",
        )
    ax.axhline(
        direct["final_full_data_residual_median"],
        color="#1f77b4",
        ls="--",
        lw=1.0,
        label=f"direct {direct['final_full_data_residual_median']:.2e}",
    )
    ax.axhline(
        reduced["final_full_data_residual_median"],
        color="#9467bd",
        ls=":",
        lw=1.0,
        label=f"reduced_r6 {reduced['final_full_data_residual_median']:.2e}",
    )
    ax.set_ylim(1e-4, 1.0)
    ax.set_xticks([3e-5] + lams_nonzero)
    ax.set_xticklabels(["0"] + [f"{v:.0e}" for v in lams_nonzero], rotation=45)
    ax.set_title(
        "E7: state_c final full-physics data residual vs lam\n"
        "(||forward(alpha_est,p_est)-d_obs||/||d_obs||, medians)",
        fontsize=9,
    )
    ax.set_xlabel(r"state weight $\lambda$")
    ax.set_ylabel("normalized full-data residual")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    # Panel 3: tu_penalty pose error vs lam_T.
    ax = axes[2]
    keep = sorted(
        lam for lam, s in tu.items()
        if lam is not None and s.get("n_success", 0) > 0
    )
    omitted = sorted(
        lam for lam, s in tu.items()
        if lam is not None and s.get("n_success", 0) == 0
    )
    for f, sty in (
        ("pose_error", dict(color="#1f77b4", marker="o", label="pose error")),
        (
            "pose_error_hidden",
            dict(color="#d62728", marker="^", label="hidden (r6 decomp)"),
        ),
    ):
        pts = [
            (lam, tu[lam][f + "_median"])
            for lam in keep
            if tu[lam].get(f + "_median") is not None
        ]
        if pts:
            ax.semilogx(
                [a for a, _ in pts],
                [b for _, b in pts],
                ls="-",
                ms=5,
                **sty,
            )
    ax.axhline(
        direct["pose_error_median"],
        color="#1f77b4",
        ls="--",
        lw=1.0,
        label=f"direct {direct['pose_error_median']:.2e}",
    )
    ax.axhline(
        reduced["pose_error_median"],
        color="#9467bd",
        ls=":",
        lw=1.0,
        label=f"reduced_r6 {reduced['pose_error_median']:.2e}",
    )
    fail_txt = (
        f"omitted (status 0, max_nfev): lam_T={omitted}"
        if omitted
        else ""
    )
    ax.set_title(
        "E7: tu_penalty pose error vs lam_T\n"
        "(raw data residual + lam_T*state_witness; medians)\n"
        + fail_txt,
        fontsize=9,
    )
    ax.set_xlabel(r"witness penalty $\lambda_T$")
    ax.set_ylabel("pose error")
    ax.set_ylim(1e-4, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    fig.suptitle(
        "E7: state-consistent reduced self-calibration (r=6, confounded regime)",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> dict:
    t_start = time.time()
    P = Params(N=16, k=12.0, M=8, R_r=1.6, R_t=2.0, phi0=0.7, s=0.12)
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    Phi = basis_matrix(xs, P)
    chi_true = Phi @ ALPHA_TRUE

    subspace = visible_pose_subspace(
        ALPHA_INIT, P_INIT, RETAINED_R, xs, h, P.k, G_D, Phi, P
    )
    assert subspace["n_vis"] == 1 and subspace["hidden_rank"] == 2

    d_true = forward(ALPHA_TRUE, P_TRUE, xs, h, P.k, G_D, Phi, P)
    obs: dict[int, np.ndarray] = {}
    for seed in NOISE_SEEDS:
        rng = np.random.default_rng(seed)
        raw = (
            rng.standard_normal(P.M) + 1j * rng.standard_normal(P.M)
        ) / np.sqrt(2.0)
        noise = raw * NOISE_RATIO * np.linalg.norm(d_true) / np.linalg.norm(raw)
        obs[int(seed)] = d_true + noise

    j_ref = float(
        np.linalg.norm(
            physical_current(
                ALPHA_INIT, P_INIT, xs, h, P.k, G_D, Phi, P
            )
        )
    )

    records: list[dict] = []
    for seed in NOISE_SEEDS:
        # A and B: analytic Jacobian baseline cases.
        records.append(
            run_baseline_method(
                "direct", seed, obs, subspace, xs, h, P.k, G_D, Phi, P
            )
        )
        records.append(
            run_baseline_method(
                "reduced_r6", seed, obs, subspace, xs, h, P.k, G_D, Phi, P
            )
        )
        # C: state_consistent sweep (lam=0 included as degenerate data-only).
        for lam in (0.0,) + STATE_C_LAMS:
            rec = run_state_c(
                lam, seed, obs, j_ref, xs, h, P.k, G_D, Phi, P, subspace
            )
            records.append(rec)
            print(
                f"state_c seed={seed} lam={lam:>6}: "
                f"success={rec['success']} opt={rec['optimality']:.2e} "
                f"pose={rec['pose_error']:.4e} "
                f"vis={rec['pose_error_visible']:.4e} "
                f"hid={rec['pose_error_hidden']:.4e} "
                f"map={rec['map_error']:.4e} nfev={rec['nfev']}"
            )
        # D: T_U penalty sweep.
        for lam_t in TU_LAMS:
            rec = run_tu_penalty(
                lam_t, seed, obs, xs, h, P.k, G_D, Phi, P, subspace
            )
            records.append(rec)
            print(
                f"tu_penalty seed={seed} lam_T={lam_t:>6}: "
                f"success={rec['success']} opt={rec['optimality']:.2e} "
                f"pose={rec['pose_error']:.4e} "
                f"vis={rec['pose_error_visible']:.4e} "
                f"hid={rec['pose_error_hidden']:.4e} "
                f"T_U={rec['T_U'] if rec['T_U'] is not None else float('nan'):.4f}"
            )

    summaries = summarize(records)
    baselines = {}
    for method in ("direct", "reduced_r6"):
        sel = [s for s in summaries if s["method"] == method]
        assert len(sel) == 1
        baselines[method] = sel[0]

    subspace_json = {
        "r": subspace["rank"],
        "n_vis": int(subspace["n_vis"]),
        "hidden_rank": int(subspace["hidden_rank"]),
        "hid_sv": [float(v) for v in subspace["hid_sv"]],
        "threshold": float(subspace["threshold"]),
        "hidden_directions": subspace["hidden_directions"],
        "hidden_direction_note": (
            "rows of Vhb (right-singular 3-vectors) for hid_sv at or below "
            "1e-8*hid_sv[0]; pose errors are decomposed with V_vis columns "
            "and this row complement at (alpha_init,p_init)"
        ),
    }

    struct_diag = structural_diagnostics(
        obs, j_ref, xs, h, P.k, G_D, Phi, P
    )

    results = {
        "experiment": "run_e7_state_consistent",
        "goal": (
            "Does a state-consistent reduced variant (state_c, explicit "
            "retained current c plus lam*(V_r(p)c-J_phys) consistency) or the "
            "T_U witness as a penalty (tu_penalty) recover the hidden pose "
            "directions that reduced_r6 freezes (r=6, n_vis=1 regime)?"
        ),
        "parameters": {
            "N": P.N,
            "k": P.k,
            "M": P.M,
            "R_r": P.R_r,
            "R_t": P.R_t,
            "phi0": P.phi0,
            "s": P.s,
            "K": K,
            "basis_centers": [
                list(c)
                for c in ((-0.15, 0.10), (0.20, -0.10), (0.00, 0.00))
            ],
            "r": RETAINED_R,
            "SNR_dB": SNR_DB,
            "noise_ratio": NOISE_RATIO,
            "alpha_true": [float(v) for v in ALPHA_TRUE],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_true": [float(v) for v in P_TRUE],
            "p_init": [float(v) for v in P_INIT],
            "noise_seeds": [int(s) for s in NOISE_SEEDS],
            "state_c_lams": [float(v) for v in (0.0,) + STATE_C_LAMS],
            "tu_lam_T": [float(v) for v in TU_LAMS],
            "least_squares": dict(LS_KWARGS),
            "d_true_norm": float(np.linalg.norm(d_true)),
            "chi_true_norm": float(np.linalg.norm(chi_true)),
        },
        "visible_pose_subspace_at_p_init_r6": subspace_json,
        "notes": [
            "state_c optimizes theta=[alpha; c_real; c_imag; p], d_model="
            "G_s(p)@(V_r(p)@c), residual=[realify(d_model-d_obs); "
            "lam*realify(V_r(p)@c-J_phys(alpha,p))/||J_phys||]; J_phys is the "
            "full physical contrast current. lam=0 is the degenerate "
            "data-only retained-current baseline.",
            "state_c data residual is NOT normalized by ||d_obs||; d entries "
            "are O(1e-4), while the dimensionless state residual is O(1), so "
            "lambda values near 1e-3 already dominate the objective. Raw "
            "residual scaling is the reason final numbers must be read "
            "per-lam, not as a universal weighting.",
            "tu_penalty appends the scalar lam_T*T_U to the raw realified "
            "data residual (also unnormalized); T_U in [0,1] so the penalty "
            "dominates for lam_T >= ~1e-3 given raw data residuals.",
            "V_r columns are Procrustes-aligned to a fixed reference frame at "
            "p_init (deterministic smooth gauge where the retained subspace "
            "is well defined); spans are unchanged.",
            "success is scipy result.success; strict convergence additionally "
            "requires optimality < 1e-8.",
            "final_full_data_residual is the full-physics forward misfit of "
            "the returned (alpha,p), including for state_c (whose minimized "
            "residual is the retained-data one).",
            "state_c never left p_init=0 for any lambda in this fixture: "
            "s6=s7 of G_s(p_init) are degenerate to machine precision exactly "
            "on the hard r=6 cut, so finite-difference p-Jacobians of the "
            "retained frame are dominated by the subspace cut discontinuity "
            "and TRF xtol-stops with a large optimality.  See "
            "structural_diagnostics and the strict-convergence counts.",
        ],
        "records": records,
        "summaries_median_over_seeds": summaries,
        "baselines": baselines,
        "structural_diagnostics": struct_diag,
    }

    json_path = HERE / "results_e7.json"
    png_path = HERE / "plot_e7.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(records, summaries, baselines, png_path)

    # Verification: JSON round-trip and PNG opens.
    with open(json_path) as fh:
        loaded = json.load(fh)
    assert loaded["records"] == results["records"]
    assert len(loaded["records"]) == 2 * (
        2 + len(STATE_C_LAMS) + 1 + len(TU_LAMS)
    )
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    print("\n--- E7 summary (medians over seeds 0,1) ---")
    print(
        f"direct     : pose={baselines['direct']['pose_error_median']:.4e} "
        f"vis={baselines['direct']['pose_error_visible_median']:.4e} "
        f"hid={baselines['direct']['pose_error_hidden_median']:.4e} "
        f"map={baselines['direct']['map_error_median']:.4e} "
        f"res={baselines['direct']['final_full_data_residual_median']:.4e} "
        f"T_U={baselines['direct']['T_U_median']:.4f}"
    )
    print(
        f"reduced_r6 : pose={baselines['reduced_r6']['pose_error_median']:.4e} "
        f"vis={baselines['reduced_r6']['pose_error_visible_median']:.4e} "
        f"hid={baselines['reduced_r6']['pose_error_hidden_median']:.4e} "
        f"map={baselines['reduced_r6']['map_error_median']:.4e} "
        f"res={baselines['reduced_r6']['final_full_data_residual_median']:.4e} "
        f"T_U={baselines['reduced_r6']['T_U_median']:.4f}"
    )
    print("\nstate_c:")
    for s in summaries:
        if s["method"] != "state_c":
            continue
        lam = s["lam"]
        if s.get("all_failed"):
            print(f"  lam={lam:>6}: ALL FAILED")
            continue
        print(
            f"  lam={lam:>6}: pose={s['pose_error_median']:.4e} "
            f"vis={s['pose_error_visible_median']:.4e} "
            f"hid={s['pose_error_hidden_median']:.4e} "
            f"map={s['map_error_median']:.4e} "
            f"full_res={s['final_full_data_residual_median']:.4e} "
            f"ret_res={s['retained_data_residual_median']:.4e} "
            f"state={s['state_resid_median']:.4e} "
            f"c_norm={s['c_norm_median']:.4e} "
            f"T_U={s['T_U_median']:.4f} ok={s['n_success']}/{s['n_seeds']} "
            f"strict={s['n_strict']}/{s['n_seeds']}"
        )
    print("\ntu_penalty:")
    for s in summaries:
        if s["method"] != "tu_penalty":
            continue
        lam = s["lam"]
        if s.get("all_failed"):
            print(f"  lam_T={lam:>6}: ALL FAILED")
            continue
        print(
            f"  lam_T={lam:>6}: pose={s['pose_error_median']:.4e} "
            f"vis={s['pose_error_visible_median']:.4e} "
            f"hid={s['pose_error_hidden_median']:.4e} "
            f"map={s['map_error_median']:.4e} "
            f"res={s['final_full_data_residual_median']:.4e} "
            f"T_U={s['T_U_median']:.4f} ok={s['n_success']}/{s['n_seeds']} "
            f"strict={s['n_strict']}/{s['n_seeds']}"
        )
    sd = struct_diag
    print("\nstructural (state_c pinning evidence):")
    print(
        "  s6-s7(G_s at p_init) = "
        f"{sd['s6_minus_s7_at_p_init']:.3e} "
        "(machine-degenerate; hard r=6 cut on the pair)"
    )
    print(
        "  FD Jacobian col norms at x0: "
        f"alpha_max={sd['numeric_jacobian_column_norms_at_x0']['alpha_max']:.2e} "
        f"c_max={sd['numeric_jacobian_column_norms_at_x0']['c_max']:.2e} "
        "p=["
        + ",".join(
            f"{v:.2e}"
            for v in sd["numeric_jacobian_column_norms_at_x0"]["p_columns"]
        )
        + "]"
    )
    print(f"\nwall time {time.time()-t_start:.1f} s")
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
