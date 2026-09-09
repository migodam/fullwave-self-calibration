"""E5 final -- nonlinear self-calibration with scipy robust least squares.

Compared with run_e5.py (adaptive LM with a projected pose Jacobian
P_perp @ B_real recomputed every iteration), this version uses a clean
reduced parameterization:

  reduced(r) : theta = [alpha; q], p = p_init + V_vis @ q
    with V_vis (3 x n_vis) the *fixed* visible pose subspace computed once
    at p_init from the geometry-lift hiding condition, and with the TRUE
    co-moving pose Jacobian B_real @ V_vis (no projection of B_real).

  direct     : theta = [alpha; p],  J = [J_alpha_real, B_real]
  wrongpose  : theta = alpha only at fixed p = (0,0,0), J = J_alpha_real

All methods use scipy.optimize.least_squares(method='trf', x_scale='jac',
max_nfev=1000, ftol=xtol=gtol=1e-10).  Deterministic: fixed seeds for the
noise, no random state anywhere else.

Visible-pose subspace at (alpha_init, p_init):
  SVD G_s -> V_r (first r right singular vectors); Q = G_s @ V_r;
  Hn = norm_cols([realify(G_s V_r), realify(J_alpha)]);
  P_perp = I - U_r U_r^T over col(Hn);
  B_red = P_perp @ B_real;  SVD B_red = Ub sb Vhb^T;
  n_vis = #{sb > 1e-8*sb[0]};  V_vis = first n_vis columns of Vhb.

Pose-error decomposition uses V_vis fixed at p_init.  For reduced runs it is
the run's own r subspace; for direct and wrongpose runs it is the r=6
subspace at the same p_init (recorded as decomposition_basis="r6") so the
confounded-regime error location can be compared across methods.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from geom_som_core import Params, pixel_grid, green_domain_matrix
from run_e2_e4 import (
    realify_complex_cols,
    realify_real_cols,
    norm_cols,
    proj_complement,
)
from run_e5 import (
    basis_matrix,
    forward,
    jacobians,
    state_witness,
    _realified_residual,
)


HERE = Path(__file__).resolve().parent

SNR_DB = 30.0
NOISE_RATIO = float(10.0 ** (-SNR_DB / 20.0))  # ~0.0316
HIDDEN_SV_REL = 1e-8
TOL = 1e-10
NOISE_SEEDS = (0, 1)
P_INIT_LABELS = ("zero", "pert_pos", "pert_neg")
ALPHA_INIT = np.array([1.0, 1.0, 0.0], dtype=float)
METHOD_ORDER = ("wrongpose", "direct", "reduced_r4", "reduced_r6")

LS_KWARGS = dict(
    method="trf",
    x_scale="jac",
    max_nfev=1000,
    ftol=1e-10,
    xtol=1e-10,
    gtol=1e-10,
)


# ---------------------------------------------------------------------------
# Geometry-lift visible pose subspace (computed once at alpha_init, p_init)
# ---------------------------------------------------------------------------
def visible_pose_subspace(
    alpha: np.ndarray,
    p: np.ndarray,
    r: int,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    P: Params,
) -> dict:
    """Reduced visible pose basis from the geometry-lift hiding condition.

    Returns the 3 x n_vis basis V_vis (columns orthonormal), the three sorted
    (descending) singular values of B_red = P_perp @ B_real, and the hidden
    directions as rows of Vhb corresponding to singular values at or below
    the 1e-8 * sb[0] cutoff.  Note: Vhb rows are the right-singular 3-vectors;
    V_vis columns are the transpose of the first n_vis rows.
    """
    jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
    _, _, Vh = np.linalg.svd(jb["G_s"], full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    Q = jb["G_s"] @ V_r
    W = jb["J_alpha"]
    Hn = norm_cols(
        np.hstack([realify_complex_cols(Q), realify_real_cols(W)])
    )
    P_perp = proj_complement(Hn, TOL)
    B_red = P_perp @ jb["B_real"]
    _, sb, Vhb = np.linalg.svd(B_red, full_matrices=False)
    sb = np.asarray(sb, dtype=float)
    threshold = HIDDEN_SV_REL * sb[0] if sb[0] > 0.0 else 0.0
    n_vis = int(np.count_nonzero(sb > threshold))
    hidden_rank = 3 - n_vis
    V_vis = Vhb[:n_vis].T  # 3 x n_vis, columns are visible pose directions
    hidden_directions = [Vhb[i].tolist() for i in range(n_vis, 3)]
    return {
        "rank": int(r),
        "n_vis": n_vis,
        "hidden_rank": hidden_rank,
        "hid_sv": [float(v) for v in sb],  # sorted descending
        "threshold": float(threshold),
        "V_vis": V_vis,
        "hidden_directions": hidden_directions,  # 3-vectors stored as rows
    }


def make_record(
    method: str,
    seed: int,
    p_init_label: str,
    alpha_init: np.ndarray,
    p_init: np.ndarray,
    alpha_true: np.ndarray,
    p_true: np.ndarray,
    d_obs: np.ndarray,
    ref_norm: float,
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    reduced_r: int | None,
    subspace: dict,
    decomp_basis_rank: int,
    decomp_subspace: dict,
) -> dict:
    """Run one scipy least-squares case and return its record."""
    K = Phi.shape[1]
    V_vis = decomp_subspace["V_vis"]  # basis used for pose-error decomposition

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
            p = p_init + subspace["V_vis"] @ q
        return _realified_residual(
            forward(alpha, p, xs, h, k, G_D, Phi, P), d_obs
        )

    def jacobian(theta: np.ndarray) -> np.ndarray:
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
            jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
            return jb["J_alpha_real"]
        if method == "direct":
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
            jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
            return np.hstack([jb["J_alpha_real"], jb["B_real"]])
        # reduced: true pose Jacobian columns restricted to the fixed V_vis.
        alpha = np.asarray(theta[:K], dtype=float)
        q = np.asarray(theta[K:], dtype=float)
        p = p_init + subspace["V_vis"] @ q
        jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
        return np.hstack(
            [jb["J_alpha_real"], jb["B_real"] @ subspace["V_vis"]]
        )

    if method == "wrongpose":
        x0 = np.asarray(alpha_init, dtype=float)
        n_vis_solver = None
    elif method == "direct":
        x0 = np.concatenate(
            [np.asarray(alpha_init, float), np.asarray(p_init, float)]
        )
        n_vis_solver = None
    else:
        n_vis_solver = int(subspace["n_vis"])
        q0 = np.zeros(n_vis_solver)
        x0 = np.concatenate([np.asarray(alpha_init, float), q0])

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
        p_est = p_init + subspace["V_vis"] @ q_est

    b_final = residual(result.x)
    final_residual = float(np.linalg.norm(b_final) / ref_norm)

    chi_true = Phi @ alpha_true
    chi_est = Phi @ alpha_est
    map_error = float(
        np.linalg.norm(chi_est - chi_true) / np.linalg.norm(chi_true)
    )
    delta_p = p_est - p_true
    pose_error = float(np.linalg.norm(delta_p))
    visible_proj = V_vis @ (V_vis.T @ delta_p)
    pose_error_visible = float(np.linalg.norm(visible_proj))
    pose_error_hidden = float(np.linalg.norm(delta_p - visible_proj))

    if reduced_r is not None:
        t_u_rank = reduced_r
    else:
        t_u_rank = 4  # disambiguated/default retained rank, as in run_e5
    t_u = state_witness(
        alpha_est, p_est, t_u_rank, xs, h, k, G_D, Phi, P
    )

    record = {
        "method": method,
        "noise_seed": int(seed),
        "p_init_label": p_init_label,
        "p_init": [float(v) for v in p_init],
        "r": int(reduced_r) if reduced_r is not None else None,
        "n_vis": int(decomp_subspace["n_vis"]),
        "decomposition_basis": f"r{decomp_basis_rank}",
        "alpha_est": [float(v) for v in alpha_est],
        "p_est": [float(v) for v in p_est],
        "pose_error": pose_error,
        "pose_error_visible": pose_error_visible,
        "pose_error_hidden": pose_error_hidden,
        "map_error": map_error,
        "final_residual": final_residual,
        "success": bool(result.success),
        "nfev": int(result.nfev),
        "optimality": float(result.optimality),
        "hid_sv": [float(v) for v in decomp_subspace["hid_sv"]],
        "hidden_rank": int(decomp_subspace["hidden_rank"]),
        "T_U": t_u,
        "T_U_rank": int(t_u_rank),
    }
    return record


def finite_difference_check(
    method: str,
    alpha_init: np.ndarray,
    p_init: np.ndarray,
    d_obs: np.ndarray,
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    subspace: dict | None,
    step: float = 1e-6,
) -> tuple[str, float]:
    """Central-difference Jacobian check at the solver start point."""
    K = Phi.shape[1]

    def residual_and_jac(theta: np.ndarray):
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
            jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
            J = jb["J_alpha_real"]
        elif method == "direct":
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
            jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
            J = np.hstack([jb["J_alpha_real"], jb["B_real"]])
        else:
            alpha = np.asarray(theta[:K], dtype=float)
            q = np.asarray(theta[K:], dtype=float)
            p = p_init + subspace["V_vis"] @ q
            jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
            J = np.hstack(
                [jb["J_alpha_real"], jb["B_real"] @ subspace["V_vis"]]
            )
        b = _realified_residual(jb["d"], d_obs)
        return b, J

    def residual_only(theta: np.ndarray) -> np.ndarray:
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
        elif method == "direct":
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
        else:
            alpha = np.asarray(theta[:K], dtype=float)
            q = np.asarray(theta[K:], dtype=float)
            p = p_init + subspace["V_vis"] @ q
        return _realified_residual(
            forward(alpha, p, xs, h, k, G_D, Phi, P), d_obs
        )

    if method == "wrongpose":
        x0 = np.asarray(alpha_init, dtype=float)
    elif method == "direct":
        x0 = np.concatenate(
            [np.asarray(alpha_init, float), np.asarray(p_init, float)]
        )
    else:
        x0 = np.concatenate(
            [
                np.asarray(alpha_init, float),
                np.zeros(int(subspace["n_vis"])),
            ]
        )
    b0, J_ana = residual_and_jac(x0)
    n = J_ana.shape[1]
    J_num = np.empty_like(J_ana)
    for j in range(n):
        xp = np.array(x0, dtype=float)
        xm = np.array(x0, dtype=float)
        xp[j] += step
        xm[j] -= step
        J_num[:, j] = (residual_only(xp) - residual_only(xm)) / (2.0 * step)
    col_norm = np.linalg.norm(J_ana, axis=0)
    scale = np.where(col_norm > 1e-12, col_norm, 1.0)
    rel = np.linalg.norm(J_num - J_ana, axis=0) / scale
    del b0
    return method, float(np.max(rel))


def make_plot(
    records: list[dict],
    png_path: Path,
) -> None:
    """Two-panel E5-final figure (seed-0 zero-init view + median errors)."""
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8))
    colors = {
        "wrongpose": "#d62728",
        "direct": "#1f77b4",
        "reduced_r4": "#2ca02c",
        "reduced_r6": "#9467bd",
    }

    # ---- Left: final normalized residual vs noise floor (seed 0, zero) ----
    ax = axes[0]
    left_values = {}
    for method in METHOD_ORDER:
        for rec in records:
            if (
                rec["method"] == method
                and rec["noise_seed"] == 0
                and rec["p_init_label"] == "zero"
            ):
                left_values[method] = rec["final_residual"]
                break
    x = np.arange(len(METHOD_ORDER))
    vals = [max(left_values[m], 1e-16) for m in METHOD_ORDER]
    bars = ax.bar(
        x,
        vals,
        0.58,
        color=[colors[m] for m in METHOD_ORDER],
        alpha=0.85,
    )
    ax.axhline(NOISE_RATIO, color="k", ls="--", lw=1.0)
    ax.text(
        0.985,
        NOISE_RATIO * 1.25,
        f"noise ratio {NOISE_RATIO:.4f}",
        transform=ax.transData,
        ha="right",
        va="bottom",
        fontsize=8,
    )
    for bar, m, v in zip(bars, METHOD_ORDER, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            v * 1.4,
            f"{v:.4e}",
            ha="center",
            va="bottom",
            fontsize=7.5,
            rotation=0,
        )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(list(METHOD_ORDER), fontsize=9)
    ax.set_ylabel(r"$\|d-d_{\rm obs}\|/\|d_{\rm obs}\|$")
    ax.set_title(
        "E5 final: normalized residual (seed 0, zero init)\n"
        "trf least_squares vs noise floor"
    )
    ax.set_ylim(0.5 * min(vals + [NOISE_RATIO]), 2.0 * max(vals))
    ax.grid(axis="y", alpha=0.3)
    ax.margins(x=0.06)

    # ---- Right: seed-0 medians over inits (pose / map / r6 decomposition) --
    ax = axes[1]
    medians = {}
    for method in METHOD_ORDER:
        sel = [
            rec for rec in records
            if rec["method"] == method and rec["noise_seed"] == 0
        ]
        pose = float(np.median([rec["pose_error"] for rec in sel]))
        mape = float(np.median([rec["map_error"] for rec in sel]))
        medians[method] = {"pose": pose, "map": mape}
        if method == "reduced_r6":
            medians[method]["visible"] = float(
                np.median([rec["pose_error_visible"] for rec in sel])
            )
            medians[method]["hidden"] = float(
                np.median([rec["pose_error_hidden"] for rec in sel])
            )

    x = np.arange(len(METHOD_ORDER))
    legend_seen: set[str] = set()
    for j, method in enumerate(METHOD_ORDER):
        pose = max(medians[method]["pose"], 1e-16)
        mape = max(medians[method]["map"], 1e-16)
        if method == "reduced_r6":
            # Decomposition of the median r6 pose error into its orthogonal
            # visible/hidden parts (shown side-by-side because a log y-axis
            # cannot honestly stack orthogonal components; pose_error is the
            # root-sum-square, not the sum, of these components).
            width = 0.16
            offsets = (-0.27, -0.09, 0.09, 0.27)
            cols = (
                (pose, "#1f77b4", "median pose error"),
                (mape, "#ff7f0e", "median map error"),
            )
            vis = max(medians[method]["visible"], 1e-16)
            hid = max(medians[method]["hidden"], 1e-16)
            cols += (
                (vis, "#2ca02c", "visible (r6)"),
                (hid, "#d62728", "hidden (r6)"),
            )
        else:
            width = 0.22
            offsets = (-0.13, 0.13)
            cols = (
                (pose, "#1f77b4", "median pose error"),
                (mape, "#ff7f0e", "median map error"),
            )
        for (val, col, lab), off in zip(cols, offsets):
            first = lab not in legend_seen
            legend_seen.add(lab)
            ax.bar(
                x[j] + off,
                val,
                width,
                color=col,
                label=lab if first else None,
            )
            ax.text(
                x[j] + off,
                val * 1.55,
                f"{val:.1e}" if method == "reduced_r6" else f"{val:.2e}",
                ha="center",
                fontsize=6 if method == "reduced_r6" else 6.5,
            )

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(list(METHOD_ORDER), fontsize=9)
    ax.set_ylabel("error")
    ax.set_title(
        "E5 final: seed-0 median errors over inits\n"
        "(r6 pose error decomposed into visible/hidden)"
    )
    ax.set_ylim(1e-4, 1.0)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(
        fontsize=7,
        ncol=2,
        loc="upper left",
        bbox_to_anchor=(1.012, 1.0),
        borderaxespad=0.0,
    )

    fig.suptitle("E5 final: scipy trf nonlinear self-calibration", y=1.02)
    fig.tight_layout()
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> dict:
    P = Params(N=16, k=12.0, M=8, R_r=1.6, R_t=2.0, phi0=0.7, s=0.12)
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    Phi = basis_matrix(xs, P)
    K = Phi.shape[1]

    alpha_true = np.array([1.5, 2.0, 0.0])
    p_true = np.array([0.08, -0.06, 0.05])
    chi_true = Phi @ alpha_true

    d_clean = forward(alpha_true, p_true, xs, h, P.k, G_D, Phi, P)
    rms = float(np.sqrt(np.mean(np.abs(d_clean) ** 2)))
    noise_std = NOISE_RATIO * rms
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

    # Fixed visible-pose subspaces at (alpha_init, p_init) for r=4 and r=6.
    subspaces = {}
    for label, p_init in p_inits.items():
        subspaces[label] = {}
        for r in (4, 6):
            sp = visible_pose_subspace(
                ALPHA_INIT, p_init, r, xs, h, P.k, G_D, Phi, P
            )
            subspaces[label][f"r{r}"] = sp

    # Analytic vs central-difference Jacobian sanity checks (seed-0 data).
    for method, r_key, sp in (
        ("wrongpose", None, None),
        ("direct", None, None),
        ("reduced_r4", "r4", subspaces["zero"]["r4"]),
        ("reduced_r6", "r6", subspaces["zero"]["r6"]),
    ):
        name, rel_err = finite_difference_check(
            method,
            ALPHA_INIT,
            p_zero,
            obs[0],
            P,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            sp,
        )
        print(
            f"FD Jacobian check [{name:>10s} "
            f"{r_key if r_key else 'n/a'}] max rel err = {rel_err:.3e}"
        )
        assert rel_err < 1e-4, f"Jacobian mismatch for {method}: {rel_err:.3e}"

    records: list[dict] = []
    for seed in NOISE_SEEDS:
        ref_norm = float(np.linalg.norm(obs[seed]))
        # wrongpose: alpha only, pose fixed at (0,0,0); decomposition uses
        # the r=6 subspace at the zero p_init.
        records.append(
            make_record(
                "wrongpose",
                seed,
                "zero",
                ALPHA_INIT,
                p_zero,
                alpha_true,
                p_true,
                obs[seed],
                ref_norm,
                P,
                xs,
                h,
                P.k,
                G_D,
                Phi,
                None,
                None,
                6,
                subspaces["zero"]["r6"],
            )
        )
        for label in P_INIT_LABELS:
            p_init = p_inits[label]
            records.append(
                make_record(
                    "direct",
                    seed,
                    label,
                    ALPHA_INIT,
                    p_init,
                    alpha_true,
                    p_true,
                    obs[seed],
                    ref_norm,
                    P,
                    xs,
                    h,
                    P.k,
                    G_D,
                    Phi,
                    None,
                    None,
                    6,
                    subspaces[label]["r6"],
                )
            )
        for method, r in (("reduced_r4", 4), ("reduced_r6", 6)):
            for label in P_INIT_LABELS:
                p_init = p_inits[label]
                records.append(
                    make_record(
                        method,
                        seed,
                        label,
                        ALPHA_INIT,
                        p_init,
                        alpha_true,
                        p_true,
                        obs[seed],
                        ref_norm,
                        P,
                        xs,
                        h,
                        P.k,
                        G_D,
                        Phi,
                        r,
                        subspaces[label][f"r{r}"],
                        r,
                        subspaces[label][f"r{r}"],
                    )
                )

    # JSON-friendly copies of the subspace data (V_vis rows for readability;
    # columns are the visible pose basis vectors, hidden_directions are rows).
    subspace_json = {}
    for label in P_INIT_LABELS:
        subspace_json[label] = {}
        for r in (4, 6):
            sp = subspaces[label][f"r{r}"]
            subspace_json[label][f"r{r}"] = {
                "rank": sp["rank"],
                "n_vis": sp["n_vis"],
                "hidden_rank": sp["hidden_rank"],
                "hid_sv": sp["hid_sv"],
                "threshold": sp["threshold"],
                "hidden_directions": sp["hidden_directions"],
                "hidden_direction_note": (
                    "rows of Vhb (right-singular 3-vectors) for hid_sv at or "
                    "below 1e-8*hid_sv[0]"
                ),
            }

    results = {
        "experiment": "run_e5_final",
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
            "SNR_dB": SNR_DB,
            "noise_ratio": NOISE_RATIO,
            "noise_std": float(noise_std),
            "alpha_true": [float(v) for v in alpha_true],
            "p_true": [float(v) for v in p_true],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_init_variants": {
                label: [float(v) for v in p] for label, p in p_inits.items()
            },
            "noise_seeds": [int(s) for s in NOISE_SEEDS],
            "least_squares": {
                "method": LS_KWARGS["method"],
                "x_scale": LS_KWARGS["x_scale"],
                "max_nfev": LS_KWARGS["max_nfev"],
                "ftol": LS_KWARGS["ftol"],
                "xtol": LS_KWARGS["xtol"],
                "gtol": LS_KWARGS["gtol"],
            },
        },
        "visible_pose_subspaces_at_p_init": subspace_json,
        "note": (
            "Reduced method optimizes (alpha, q) with p = p_init + V_vis @ q "
            "and the TRUE pose Jacobian B_real @ V_vis; V_vis is fixed once "
            "at (alpha_init, p_init) from col(realify(G_s V_r), "
            "realify(J_alpha)). hid_sv are sorted descending singular values "
            "of B_red = P_perp @ B_real at p_init. n_vis = count(hid_sv > "
            "1e-8*hid_sv[0]); hidden_rank = 3 - n_vis. Records: wrongpose "
            "and direct use the r=6 V_vis decomposition basis at p_init "
            "(decomposition_basis='r6'); reduced records decompose on their "
            "own V_vis basis. pose_error_visible/hidden are the orthogonal "
            "projection norms of p_est-p_true onto V_vis and its complement "
            "(pose_error is their root-sum-square). T_U is the receiver-rx "
            "state-consistency witness at the final estimate at rank 4 for "
            "wrongpose/direct/reduced_r4 and rank 6 for reduced_r6."
        ),
        "records": records,
    }

    json_path = HERE / "results_e5_final.json"
    png_path = HERE / "plot_e5_final.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(records, png_path)

    # Verification: JSON round-trip and PNG loads.
    with open(json_path) as fh:
        loaded = json.load(fh)
    assert len(loaded["records"]) == 20, len(loaded["records"])
    assert loaded["records"] == results["records"]
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    print("\n--- E5 final summary ---")
    print("visible-pose subspaces at (alpha_init, p_init):")
    for label in P_INIT_LABELS:
        for r in (4, 6):
            sp = loaded["visible_pose_subspaces_at_p_init"][label][f"r{r}"]
            print(
                f"  {label:>9s} r{r}: n_vis={sp['n_vis']} "
                f"hidden_rank={sp['hidden_rank']} hid_sv="
                + ",".join(f"{v:.2e}" for v in sp["hid_sv"])
            )
    for method in METHOD_ORDER:
        sel = [
            rec for rec in loaded["records"]
            if rec["method"] == method and rec["noise_seed"] == 0
        ]
        pose = float(np.median([rec["pose_error"] for rec in sel]))
        mape = float(np.median([rec["map_error"] for rec in sel]))
        res = float(np.median([rec["final_residual"] for rec in sel]))
        t_u = float(np.median([rec["T_U"] for rec in sel]))
        ok = sum(int(rec["success"]) for rec in sel)
        print(
            f"{method:>11s} seed0 median: final_residual={res:.4e} "
            f"pose_error={pose:.4e} map_error={mape:.4e} T_U={t_u:.4f} "
            f"success={ok}/{len(sel)}"
        )
    for rec in loaded["records"]:
        if rec["method"] == "reduced_r6":
            print(
                f"  r6 seed{rec['noise_seed']} {rec['p_init_label']:>9s}: "
                f"pose={rec['pose_error']:.4e} vis={rec['pose_error_visible']:.4e} "
                f"hid={rec['pose_error_hidden']:.4e} map={rec['map_error']:.4e} "
                f"res={rec['final_residual']:.4e} nfev={rec['nfev']} "
                f"success={rec['success']} T_U={rec['T_U']:.4f}"
            )
    print(
        f"noise_ratio = {loaded['parameters']['noise_ratio']:.4e} "
        f"(30 dB)"
    )
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
