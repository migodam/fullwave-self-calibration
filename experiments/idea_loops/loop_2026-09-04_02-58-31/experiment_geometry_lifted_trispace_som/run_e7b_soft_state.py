"""E7b -- soft-spectral-filter state-consistent retained-current self-calibration.

Scene is identical to E7 (single-transmitter co-moving, N=16, k=12, M=8,
R_r=1.6, R_t=2.0, phi0=0.7, s=0.12, K=3, alpha_true=[1.5,2,0],
alpha_init=[1,1,0], p_true=[0.08,-0.06,0.05], p_init=0, SNR 30 dB, seeds 0,1).

E7's hard rank-r=6 retained basis freezes at p_init because s6=s7 of
G_s(p_init) are degenerate to machine precision and lie exactly on the hard
cut, so the "first r right singular vectors" frame is conically singular.

E7b replaces the hard cut with a smooth spectral filter:

  w(s) = 0.5*(1 + tanh(a*(s - s_cut))),  transition parameters fixed once at
  p_init (s5=s[4], s8=s[7], s_cut=0.5*(s5+s8), a=4/max(s5-s8,1e-12)),

and solves the state-consistent retained-current problem (state_c_soft):

  theta = [alpha(K); c_real(M); c_imag(M); p(3)]
  d_model = G_s(p) @ (V_soft(p) @ c)
  residual = [realify(d_model-d_obs);
              lam*realify(V_soft(p)@c - J_phys(alpha,p))/||J_phys||]

Gauge note (important implementation finding): applying smooth weights to the
raw numpy SVD *individual* singular-vector frame does not by itself make the
parameterised map smooth.  The two right singular vectors of the numerically
degenerate s6=s7 pair mix by O(1) under arbitrarily small pose moves (a pure
basis-selection artifact), so raw V_soft = Vh^H diag(w) has conical FD
columns even though the weighted *subspace* is smooth.  We therefore use the
same span/operator with a smooth coordinate gauge: V_soft(p) is rotated
(Procrustes, V_soft -> V_soft @ Q, Q unitary) to be closest to the fixed
soft frame V_soft(p_init).  This only re-coordinates the free coefficient c
(d_model, the attainable state-residual set and the least-squares optimum
over (alpha,p) are unchanged), but it makes the finite-difference Jacobian
well behaved.  Both raw-frame and aligned-frame structural FD columns are
recorded for comparison.

Deterministic: two fresh numpy default_rng draws (seeds 0,1) for noise.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

# Bound BLAS threading inside worker processes (many independent solves).
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "2")

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
from run_e7_state_consistent import (
    realify,
    physical_current,
    decompose_pose,
    finalize_record,
    failed_record,
)


HERE = Path(__file__).resolve().parent

K = 3
RETAINED_R = 6          # E7 hard rank (used only for T_U and the V_vis basis)
SOFT_M = 8              # c length = full receiver count M
SNR_DB = 30.0
NOISE_RATIO = float(10.0 ** (-SNR_DB / 20.0))
NOISE_SEEDS = (0, 1)

ALPHA_TRUE = np.array([1.5, 2.0, 0.0], dtype=float)
ALPHA_INIT = np.array([1.0, 1.0, 0.0], dtype=float)
P_TRUE = np.array([0.08, -0.06, 0.05], dtype=float)
P_INIT = np.zeros(3)

STATE_C_SOFT_LAMS = (1e-4, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0)
ALL_SOFT_LAMS = (0.0,) + STATE_C_SOFT_LAMS

LS_BASE_KWARGS = dict(
    method="trf",
    x_scale="jac",
    max_nfev=1000,
    ftol=1e-10,
    xtol=1e-10,
    gtol=1e-10,
)
LS_SOFT_KWARGS = dict(
    method="trf",
    x_scale="jac",
    max_nfev=3000,
    ftol=1e-10,
    xtol=1e-10,
    gtol=1e-10,
)
J_FLOOR_REL = 1e-9


def build_fixture(seed: int | None = None) -> dict:
    """Deterministic fixture (optionally only the requested noise seed)."""
    P = Params(N=16, k=12.0, M=8, R_r=1.6, R_t=2.0, phi0=0.7, s=0.12)
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    Phi = basis_matrix(xs, P)
    chi_true = Phi @ ALPHA_TRUE
    subspace = visible_pose_subspace(
        ALPHA_INIT, P_INIT, RETAINED_R, xs, h, P.k, G_D, Phi, P
    )
    assert subspace["n_vis"] == 1 and subspace["hidden_rank"] == 2
    s_cut, a = soft_transition_at_p_init(xs, h, P.k, P)
    d_true = forward(ALPHA_TRUE, P_TRUE, xs, h, P.k, G_D, Phi, P)
    obs: dict[int, np.ndarray] = {}
    seeds = (int(seed),) if seed is not None else NOISE_SEEDS
    for s in seeds:
        rng = np.random.default_rng(s)
        raw = (
            rng.standard_normal(P.M) + 1j * rng.standard_normal(P.M)
        ) / np.sqrt(2.0)
        noise = raw * NOISE_RATIO * np.linalg.norm(d_true) / np.linalg.norm(raw)
        obs[int(s)] = d_true + noise
    j_ref = float(
        np.linalg.norm(
            physical_current(ALPHA_INIT, P_INIT, xs, h, P.k, G_D, Phi, P)
        )
    )
    return {
        "P": P,
        "xs": xs,
        "h": h,
        "G_D": G_D,
        "Phi": Phi,
        "chi_true": chi_true,
        "subspace": subspace,
        "s_cut": s_cut,
        "a": a,
        "d_true": d_true,
        "obs": obs,
        "j_ref": j_ref,
    }


def _soft_worker(job: tuple[int, float]) -> dict:
    """Module-level pool worker: one (seed, lam) state_c_soft solve."""
    seed, lam = job
    f = build_fixture(seed)
    rec = run_state_c_soft(
        lam,
        int(seed),
        f["obs"],
        f["j_ref"],
        f["s_cut"],
        f["a"],
        f["xs"],
        f["h"],
        f["P"].k,
        f["G_D"],
        f["Phi"],
        f["P"],
        f["subspace"],
    )
    return rec


def make_soft_filter(
    p: np.ndarray,
    s_cut: float,
    a: float,
    xs: np.ndarray,
    h: float,
    k: float,
    P: Params,
    V_ref: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """G_s(p), V_soft(p), s(p), w(p) for the fixed soft transition.

    V_soft = Vh^H diag(w) (Np x M), optionally Procrustes-rotated to be
    closest (as columns, under a unitary right rotation Q) to V_ref.  The
    rotation only re-coordinates the free current coefficients c; the span of
    the retained-soft model is unchanged.
    """
    y, _ = receiver_positions(p, P)
    G_s = data_matrix(y, xs, h, k)
    _, s, Vh = np.linalg.svd(G_s, full_matrices=False)
    w = 0.5 * (1.0 + np.tanh(a * (s - s_cut)))
    V_soft = Vh.conj().T @ np.diag(w)
    if V_ref is not None:
        O = V_soft.conj().T @ V_ref
        Uo, _, Vho = np.linalg.svd(O, full_matrices=False)
        # Deterministic unitary gauge aligned to the fixed p_init frame.
        # Convention Q = Uo @ Vho (same as run_e7.retained_basis).  Numerically
        # this branch is continuous across the near-degenerate singular
        # vectors at every finite-difference scale tested (1.5e-8 .. 1e-3),
        # whereas the raw frame and the alternative closest-rotation polar
        # factor remain conically singular inside the degenerate s6=s7 pair.
        V_soft = V_soft @ (Uo @ Vho)
    return G_s, V_soft, s, w


def soft_transition_at_p_init(
    xs: np.ndarray, h: float, k: float, P: Params
) -> tuple[float, float]:
    y0, _ = receiver_positions(P_INIT, P)
    s = np.linalg.svd(data_matrix(y0, xs, h, k), compute_uv=False)
    s5, s8 = float(s[4]), float(s[7])
    s_cut = 0.5 * (s5 + s8)
    a = 4.0 / max(s5 - s8, 1e-12)
    return s_cut, a


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
    """direct and reduced_r6 with analytic Jacobians (E5-final/E7 style)."""
    d_obs = obs[seed]
    V_vis = subspace["V_vis"]  # 3 x n_vis fixed at (alpha_init, p_init), r=6

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
        return np.hstack([jb["J_alpha_real"], jb["B_real"] @ V_vis])

    if method == "direct":
        x0 = np.concatenate([ALPHA_INIT, P_INIT])
    else:
        x0 = np.concatenate([ALPHA_INIT, np.zeros(int(subspace["n_vis"]))])

    result = least_squares(residual, x0, jac=jacobian, **LS_BASE_KWARGS)
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


def run_state_c_soft(
    lam: float,
    seed: int,
    obs: dict[int, np.ndarray],
    j_ref: float,
    s_cut: float,
    a: float,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
    subspace: dict,
) -> dict:
    """state_c_soft at weight lam with the fixed soft filter and gauge."""
    d_obs = obs[seed]
    _, V_ref, _, _ = make_soft_filter(
        P_INIT, s_cut, a, xs, h, k, P, V_ref=None
    )
    J0 = physical_current(ALPHA_INIT, P_INIT, xs, h, k, G_D, Phi, P)
    c0 = np.linalg.lstsq(V_ref, J0, rcond=None)[0]
    x0 = np.concatenate([ALPHA_INIT, c0.real, c0.imag, P_INIT])

    def residual(theta: np.ndarray) -> np.ndarray:
        alpha = np.asarray(theta[:K], dtype=float)
        c = (
            np.asarray(theta[K : K + SOFT_M], dtype=float)
            + 1j * np.asarray(theta[K + SOFT_M : K + 2 * SOFT_M], dtype=float)
        )
        p = np.asarray(theta[K + 2 * SOFT_M :], dtype=float)
        G_s, V_soft, _, _ = make_soft_filter(
            p, s_cut, a, xs, h, k, P, V_ref=V_ref
        )
        J_phys = physical_current(alpha, p, xs, h, k, G_D, Phi, P)
        r_data = realify(G_s @ (V_soft @ c) - d_obs)
        denom = max(float(np.linalg.norm(J_phys)), J_FLOOR_REL * j_ref)
        r_state = lam * realify(V_soft @ c - J_phys) / denom
        return np.concatenate([r_data, r_state])

    try:
        result = least_squares(residual, x0, **LS_SOFT_KWARGS)
        alpha_est = np.asarray(result.x[:K], dtype=float)
        c_est = (
            np.asarray(result.x[K : K + SOFT_M], dtype=float)
            + 1j
            * np.asarray(
                result.x[K + SOFT_M : K + 2 * SOFT_M], dtype=float
            )
        )
        p_est = np.asarray(result.x[K + 2 * SOFT_M :], dtype=float)
        G_s, V_soft, _, _ = make_soft_filter(
            p_est, s_cut, a, xs, h, k, P, V_ref=V_ref
        )
        cur = V_soft @ c_est
        J_phys = physical_current(
            alpha_est, p_est, xs, h, k, G_D, Phi, P
        )
        denom = max(
            float(np.linalg.norm(J_phys)), J_FLOOR_REL * j_ref
        )
        extra = {
            "retained_data_residual": float(
                np.linalg.norm(G_s @ cur - d_obs) / np.linalg.norm(d_obs)
            ),
            "state_resid": float(np.linalg.norm(cur - J_phys) / denom),
            "c_norm": float(np.linalg.norm(c_est)),
        }
        rec = finalize_record(
            method="state_c_soft",
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
        # Keep the caller-visible success flag aligned with status>0 too.
        rec["success"] = bool(result.success and result.status > 0)
        return rec
    except Exception as exc:  # noqa: BLE001 - a lam failure is a record
        return failed_record(
            method="state_c_soft",
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


def summarize(records: list[dict]) -> list[dict]:
    """Median over noise seeds per (method, lam)."""
    keys: dict[tuple, list[dict]] = {}
    for rec in records:
        keys.setdefault((rec["method"], rec["lam"]), []).append(rec)
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
            "n_frozen_at_init": int(
                sum(bool(r.get("p_frozen_at_init")) for r in group)
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
        if method == "state_c_soft":
            for f in ("retained_data_residual", "state_resid", "c_norm"):
                vals = [r[f] for r in ok if r.get(f) is not None]
                base[f + "_median"] = float(np.median(vals)) if vals else None
        out.append(base)
    return out


def structural_diagnostics(
    obs: dict[int, np.ndarray],
    j_ref: float,
    s_cut: float,
    a: float,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
) -> dict:
    """Raw- and aligned-frame FD columns of the state_c_soft residual at x0.

    Centered finite differences in each pose coordinate p[c] at eps=1e-4 and
    eps=1e-3, with c=c0, alpha=alpha_init, lam=1e-1 (seed-0 observations).
    Full-column norm, data-part (first 2M) norm, state-part norm and the
    ratio full/data-part are reported for both the raw spec frame
    V_soft=Vh^H diag(w) and the Procrustes gauge used by the solves.
    """
    d_obs = obs[0]
    G_s0, V_ref, s0, w0 = make_soft_filter(
        P_INIT, s_cut, a, xs, h, k, P, V_ref=None
    )
    J0 = physical_current(ALPHA_INIT, P_INIT, xs, h, k, G_D, Phi, P)
    c0 = np.linalg.lstsq(V_ref, J0, rcond=None)[0]
    x0 = np.concatenate([ALPHA_INIT, c0.real, c0.imag, P_INIT])
    lam_diag = 1e-1
    n_p = 3
    n_data = 2 * P.M

    def residual_fn(align: bool):
        def res(theta: np.ndarray) -> np.ndarray:
            alpha = np.asarray(theta[:K], dtype=float)
            c = (
                np.asarray(theta[K : K + SOFT_M], dtype=float)
                + 1j
                * np.asarray(
                    theta[K + SOFT_M : K + 2 * SOFT_M], dtype=float
                )
            )
            p = np.asarray(theta[K + 2 * SOFT_M :], dtype=float)
            if align:
                G_s, V_soft, _, _ = make_soft_filter(
                    p, s_cut, a, xs, h, k, P, V_ref=V_ref
                )
            else:
                y, _ = receiver_positions(p, P)
                G_s = data_matrix(y, xs, h, k)
                _, sv, Vh = np.linalg.svd(G_s, full_matrices=False)
                wv = 0.5 * (1.0 + np.tanh(a * (sv - s_cut)))
                V_soft = Vh.conj().T @ np.diag(wv)
            J_phys = physical_current(
                alpha, p, xs, h, k, G_D, Phi, P
            )
            r_data = realify(G_s @ (V_soft @ c) - d_obs)
            denom = max(float(np.linalg.norm(J_phys)), J_FLOOR_REL * j_ref)
            r_state = lam_diag * realify(V_soft @ c - J_phys) / denom
            return np.concatenate([r_data, r_state])

        return res

    jb0 = jacobians(ALPHA_INIT, P_INIT, xs, h, k, G_D, Phi, P)
    direct_data_pose_col_norms = [
        float(np.linalg.norm(jb0["B_real"][:, c])) for c in range(n_p)
    ]

    out = {}
    p_slice = slice(K + 2 * SOFT_M, x0.size)
    for align in (False, True):
        res = residual_fn(align)
        per_eps = {}
        for eps in (1e-4, 1e-3):
            cols = []
            for c in range(n_p):
                xp = x0.copy()
                xm = x0.copy()
                xp[p_slice.start + c] += eps
                xm[p_slice.start + c] -= eps
                col = (res(xp) - res(xm)) / (2.0 * eps)
                full = float(np.linalg.norm(col))
                data_part = float(np.linalg.norm(col[:n_data]))
                state_part = float(np.linalg.norm(col[n_data:]))
                cols.append(
                    {
                        "col_norm": full,
                        "data_part_norm": data_part,
                        "state_part_norm": state_part,
                        "ratio_col_to_data_part": (
                            full / data_part if data_part > 0 else None
                        ),
                    }
                )
            per_eps[f"eps_{eps:.0e}"] = cols
        out["aligned_frame" if align else "raw_frame"] = per_eps

    return {
        "soft_filter_at_p_init": {
            "singular_values_s": [float(v) for v in s0],
            "s5_index4": float(s0[4]),
            "s6_index5": float(s0[5]),
            "s7_index6": float(s0[6]),
            "s8_index7": float(s0[7]),
            "weights_w": [float(v) for v in w0],
            "w5_index4": float(w0[4]),
            "w6_index5": float(w0[5]),
            "w7_index6": float(w0[6]),
            "w8_index7": float(w0[7]),
            "s_cut": float(s_cut),
            "slope_a": float(a),
            "soft_eff_rank_sum_w": float(np.sum(w0)),
            "Q_soft_finite": bool(np.all(np.isfinite(G_s0 @ V_ref))),
            "Q_soft_condition_number": float(np.linalg.cond(G_s0 @ V_ref)),
        },
        "direct_analytic_B_real_column_norms_at_p_init": direct_data_pose_col_norms,
        "residual_frame_fd_columns_at_lam_1e-1": out,
        "note": (
            "Raw V_soft=Vh^H diag(w) still has conical FD columns because the "
            "numpy SVD chooses an arbitrary individual frame inside the "
            "numerically degenerate s6=s7 pair, and that frame mixes by O(1) "
            "under O(eps) pose moves.  Smooth weights remove the hard-cut "
            "cliff but not the raw singular-vector basis-selection artifact.  "
            "The Procrustes gauge (rotate V_soft(p) to closest V_soft(p_init)"
            ", re-coordinating free c only) makes the FD columns "
            "eps-converged and moderate; state_c_soft solves use that gauge.  "
            "Runs were also probed against the raw frame (same parameters)."
        ),
    }


def pose_identifiability_probe(
    obs: dict[int, np.ndarray],
    s_cut: float,
    a: float,
    xs,
    h,
    k,
    G_D,
    Phi,
    P,
    j_ref: float,
) -> dict:
    """Min-over-c objective along p=tau*p_true (evidence for pose flatness)."""
    d_obs = obs[0]
    _, V_ref, _, _ = make_soft_filter(
        P_INIT, s_cut, a, xs, h, k, P, V_ref=None
    )
    out = {}
    for lam in (1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
        rows = {}
        for alpha, alab in ((ALPHA_INIT, "alpha_init"), (ALPHA_TRUE, "alpha_true")):
            vals = {}
            for tau in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
                p = tau * P_TRUE
                G_s, V_soft, _, _ = make_soft_filter(
                    p, s_cut, a, xs, h, k, P, V_ref=V_ref
                )
                J = physical_current(alpha, p, xs, h, k, G_D, Phi, P)
                A = np.vstack([G_s @ V_soft, lam * V_soft])
                b = np.concatenate([d_obs, lam * J])
                c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
                denom = max(float(np.linalg.norm(J)), J_FLOOR_REL * j_ref)
                r = np.concatenate(
                    [
                        realify(G_s @ (V_soft @ c) - d_obs),
                        lam * realify(V_soft @ c - J) / denom,
                    ]
                )
                vals[float(tau)] = 0.5 * float(r @ r)
            rows[alab] = vals
        out[float(lam)] = rows
    return out


def make_plot(
    records: list[dict],
    summaries: list[dict],
    baselines: dict,
    png_path: Path,
) -> None:
    def meds(method: str) -> dict:
        return {
            s["lam"]: s
            for s in summaries
            if s["method"] == method and s.get("pose_error_median") is not None
        }

    sc = meds("state_c_soft")
    direct = baselines["direct"]
    reduced = baselines["reduced_r6"]
    soft_lams = [lam for lam in ALL_SOFT_LAMS if lam is not None]
    lams_nonzero = [lam for lam in soft_lams if lam != 0.0]

    def x_for(lam):
        return 3e-5 if lam == 0.0 else lam

    fig, axes = plt.subplots(1, 3, figsize=(19.5, 5.8))

    # Panel 1: pose error and hidden component vs lam.
    ax = axes[0]
    styles = {
        "pose_error": dict(color="#1f77b4", marker="o", label="pose error"),
        "pose_error_hidden": dict(
            color="#d62728", marker="^", label="hidden (r6 decomp)"
        ),
        "pose_error_visible": dict(
            color="#2ca02c", marker="s", label="visible (r6 decomp)"
        ),
    }
    for f, sty in styles.items():
        pts = [
            (x_for(lam), s[f + "_median"])
            for lam, s in sc.items()
            if s.get(f + "_median") is not None
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
    ax.axhline(
        reduced["pose_error_hidden_median"],
        color="#d62728",
        ls=":",
        lw=0.8,
        label=f"r6 hidden {reduced['pose_error_hidden_median']:.2e}",
    )
    failed = [lam for lam in soft_lams if lam not in sc]
    fail_txt = f"failed lam={failed}" if failed else ""
    ax.set_title(
        "E7b: state_c_soft pose error vs lam\n"
        "(soft spectral filter, Procrustes gauge; medians over seeds)\n"
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

    # Panel 2: final full-physics data residual vs lam.
    ax = axes[1]
    pts = [
        (x_for(lam), s["final_full_data_residual_median"])
        for lam, s in sc.items()
        if s.get("final_full_data_residual_median") is not None
    ]
    if pts:
        ax.semilogx(
            [a for a, _ in pts],
            [b for _, b in pts],
            color="#1f77b4",
            marker="o",
            ls="-",
            ms=5,
            label="state_c_soft",
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
        "E7b: state_c_soft final full-physics data residual vs lam\n"
        "(||forward(alpha_est,p_est)-d_obs||/||d_obs||, medians)",
        fontsize=9,
    )
    ax.set_xlabel(r"state weight $\lambda$")
    ax.set_ylabel("normalized full-data residual")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    # Panel 3: state consistency and retained-data residuals vs lam.
    ax = axes[2]
    for f, sty in (
        ("state_resid", dict(color="#d62728", marker="^", label="state_resid")),
        (
            "retained_data_residual",
            dict(color="#1f77b4", marker="o", label="retained data residual"),
        ),
    ):
        pts = [
            (x_for(lam), s[f + "_median"])
            for lam, s in sc.items()
            if s.get(f + "_median") is not None and s[f + "_median"] > 0
        ]
        if pts:
            ax.semilogx(
                [a for a, _ in pts],
                [b for _, b in pts],
                ls="-",
                ms=5,
                **sty,
            )
    ax.set_xticks([3e-5] + lams_nonzero)
    ax.set_xticklabels(["0"] + [f"{v:.0e}" for v in lams_nonzero], rotation=45)
    ax.set_title(
        "E7b: state_c_soft consistency residuals vs lam\n"
        "(||V_soft c - J_phys||/||J_phys|| and retained-data misfit, medians)",
        fontsize=9,
    )
    ax.set_xlabel(r"state weight $\lambda$")
    ax.set_ylabel("normalized residual")
    ax.set_ylim(1e-5, 2.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    fig.suptitle(
        "E7b: soft-spectral-filter state-consistent self-calibration "
        "(r=6/n_vis=1 regime)",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> dict:
    t_start = time.time()
    fx = build_fixture()
    P, xs, h, G_D, Phi = (
        fx["P"],
        fx["xs"],
        fx["h"],
        fx["G_D"],
        fx["Phi"],
    )
    chi_true, subspace = fx["chi_true"], fx["subspace"]
    s_cut, a = fx["s_cut"], fx["a"]
    d_true, obs, j_ref = fx["d_true"], fx["obs"], fx["j_ref"]

    records: list[dict] = []
    for seed in NOISE_SEEDS:
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

    tasks = [
        (int(seed), float(lam))
        for seed in NOISE_SEEDS
        for lam in ALL_SOFT_LAMS
    ]
    soft_records: list[dict] | None = None
    try:
        import multiprocessing as mp

        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=4) as pool:
            soft_records = pool.map(_soft_worker, tasks, chunksize=1)
        print("state_c_soft: parallel pool finished", flush=True)
    except Exception as exc:  # noqa: BLE001 - serial fallback keeps results
        print(
            f"state_c_soft: pool unavailable ({type(exc).__name__}: {exc}); "
            "running serially",
            flush=True,
        )
        soft_records = [_soft_worker(t) for t in tasks]

    for seed in NOISE_SEEDS:
        for lam in ALL_SOFT_LAMS:
            rec = soft_records[
                list(NOISE_SEEDS).index(seed) * len(ALL_SOFT_LAMS)
                + list(ALL_SOFT_LAMS).index(lam)
            ]
            records.append(rec)
            print(
                f"state_c_soft seed={seed} lam={lam:>6}: "
                f"success={rec['success']} status={rec['status']} "
                f"opt={rec['optimality']:.2e} nfev={rec['nfev']} "
                f"pose={rec['pose_error']:.4e} "
                f"vis={rec['pose_error_visible']:.4e} "
                f"hid={rec['pose_error_hidden']:.4e} "
                f"map={rec['map_error']:.4e} "
                f"frozen={rec.get('p_frozen_at_init')}", flush=True
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
        "V_vis": [list(row) for row in subspace["V_vis"].T],
    }

    struct_diag = structural_diagnostics(
        obs, j_ref, s_cut, a, xs, h, P.k, G_D, Phi, P
    )
    probe = pose_identifiability_probe(
        obs, s_cut, a, xs, h, P.k, G_D, Phi, P, j_ref
    )

    results = {
        "experiment": "run_e7b_soft_state",
        "goal": (
            "E7b: replace E7's hard rank-r=6 retained basis with a smooth "
            "spectral filter (tanh weights fixed at p_init) so the "
            "state-consistent retained-current method has a well-behaved "
            "Jacobian, then test whether state consistency can recover the "
            "hidden pose directions that reduced_r6 freezes."
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
            "hard_rank_r_reference": RETAINED_R,
            "soft_c_dim": SOFT_M,
            "SNR_dB": SNR_DB,
            "noise_ratio": NOISE_RATIO,
            "alpha_true": [float(v) for v in ALPHA_TRUE],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_true": [float(v) for v in P_TRUE],
            "p_init": [float(v) for v in P_INIT],
            "noise_seeds": [int(s) for s in NOISE_SEEDS],
            "state_c_soft_lams": [float(v) for v in ALL_SOFT_LAMS],
            "least_squares_baselines": dict(LS_BASE_KWARGS),
            "least_squares_state_c_soft": dict(LS_SOFT_KWARGS),
            "d_true_norm": float(np.linalg.norm(d_true)),
            "chi_true_norm": float(np.linalg.norm(chi_true)),
        },
        "visible_pose_subspace_at_p_init_r6": subspace_json,
        "soft_filter_transition": {
            "s_cut": float(s_cut),
            "a": float(a),
            "fixed_at_p_init": True,
            "formula": "w(s)=0.5*(1+tanh(a*(s-s_cut)))",
            "frame_gauge": (
                "V_soft(p) = Vh(p)^H diag(w(s(p))); solves rotate each "
                "V_soft(p) by the unitary Procrustes factor that makes it "
                "closest to V_soft(p_init).  Raw-frame values are recorded "
                "in structural_diagnostics."
            ),
        },
        "notes": [
            "state_c_soft optimizes theta=[alpha;c_real;c_imag;p] with "
            "d_model=G_s(p)@(V_soft(p)@c) and residual=[realify(d_model-"
            "d_obs); lam*realify(V_soft(p)@c-J_phys(alpha,p))/||J_phys||]; "
            "c is complex length M=8.",
            "Structural finding: smooth weights alone do not smooth the raw "
            "SVD individual-vector frame inside the degenerate s6=s7 pair; "
            "the raw V_soft FD p-columns remain conical.  A Procrustes gauge "
            "aligned to the fixed p_init soft frame re-coordinates only the "
            "free c and yields eps-converged moderate FD columns.",
            "The min-over-c objective is essentially flat in p along "
            "p_init->p_true (see pose_identifiability_probe): with M free "
            "soft coefficients, the retained model can fit the data at "
            "essentially every pose and the state-consistency cost changes "
            "only at ~1e-8 levels for lam<=1e-3, while for lam>=1e-2 it "
            "slightly favors p_init.  state_c_soft is therefore expected to "
            "leave p near p_init (frozen) rather than recover hidden pose.",
            "state_c_soft data residual is NOT normalized by ||d_obs||; d "
            "entries are O(1e-4) while the dimensionless state residual is "
            "O(1), so per-lam numbers must be read with that scaling in "
            "mind (same convention as E7 state_c).",
            "success = scipy result.success and status>0; n_strict counts "
            "success with optimality<1e-8.  High optimality at ftol/xtol "
            "stops is expected on the pose-flat objective.",
            "final_full_data_residual is the full-physics forward misfit of "
            "the returned (alpha,p), including for state_c_soft.",
            "Baselines (direct, reduced_r6) were recomputed once for this "
            "fixture with analytic Jacobians and max_nfev=1000.",
        ],
        "records": records,
        "summaries_median_over_seeds": summaries,
        "baselines": baselines,
        "structural_diagnostics": struct_diag,
        "pose_identifiability_probe": probe,
    }

    json_path = HERE / "results_e7b.json"
    png_path = HERE / "plot_e7b.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(records, summaries, baselines, png_path)

    with open(json_path) as fh:
        loaded = json.load(fh)
    assert loaded["records"] == results["records"]
    assert len(loaded["records"]) == 2 * (2 + len(ALL_SOFT_LAMS))
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    print("\n--- E7b summary (medians over seeds 0,1) ---", flush=True)
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
    print("state_c_soft:")
    for s in summaries:
        if s["method"] != "state_c_soft":
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
            f"strict={s['n_strict']}/{s['n_seeds']} "
            f"frozen={s['n_frozen_at_init']}/{s['n_seeds']}"
        )
    sf = struct_diag["soft_filter_at_p_init"]
    print("\nsoft filter at p_init:")
    print(
        "  s5,s6,s7,s8 =",
        f"{sf['s5_index4']:.6e}, {sf['s6_index5']:.6e}, "
        f"{sf['s7_index6']:.6e}, {sf['s8_index7']:.6e}",
    )
    print(
        "  w5,w6,w7,w8 =",
        f"{sf['w5_index4']:.6f}, {sf['w6_index5']:.6f}, "
        f"{sf['w7_index6']:.6f}, {sf['w8_index7']:.6f}",
    )
    print(
        f"  s_cut={sf['s_cut']:.6e} a={sf['slope_a']:.3f} "
        f"soft_eff_rank={sf['soft_eff_rank_sum_w']:.4f} "
        f"cond(Q_soft)={sf['Q_soft_condition_number']:.4f}"
    )
    print("structural FD p-columns at lam=1e-1 (raw / aligned, eps=1e-4,1e-3):")
    for frame in ("raw_frame", "aligned_frame"):
        parts = []
        for eps in ("eps_1e-04", "eps_1e-03"):
            norms = [
                f"{c['col_norm']:.3e}"
                for c in struct_diag["residual_frame_fd_columns_at_lam_1e-1"][
                    frame
                ][eps]
            ]
            parts.append(eps + "=[" + ",".join(norms) + "]")
        print(" ", frame, " ".join(parts))
    print(f"\nwall time {time.time()-t_start:.1f} s")
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
