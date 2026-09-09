"""E9 -- bounded grid/scene robustness check of the core E5 claims.

Reuses the E5-final nonlinear machinery (single transmitter, co-moving pose
p=(tx,ty,theta); scipy trf least_squares with analytic Jacobians; fixed
visible-pose subspaces V_vis computed once at (alpha_init, p_init) from the
geometry-lift hiding condition).

Scene/grid matrix:
  * base scene N=24 and N=32 (alpha_true = [1.5, 2.0, 0.0])
  * alt  scene N=16  (same basis centers/geometry, alpha_true = [2.0, -1.0, 0.8])
Methods: wrongpose (alpha only, p=0), direct ([alpha; p]), reduced_r4,
reduced_r6; p_init_label = "zero" everywhere; 30 dB noise, seeds 0,1,2.

The core E5 claims being re-checked are:
  * reduced_r4 matches direct (lossless for identifiable pose), and
  * reduced_r6 freezes hidden directions with larger pose error whenever
    hidden_rank > 0.
This file adds larger deterministic grids (N=24, N=32) and one alternate
contrast scene to test whether those claims survive.
"""

from __future__ import annotations

import json
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
    receiver_positions,
    data_matrix,
)
from run_e5 import (
    basis_matrix,
    forward,
    jacobians,
    state_witness,
    _realified_residual,
)
from run_e5_final import visible_pose_subspace, LS_KWARGS


HERE = Path(__file__).resolve().parent

SNR_DB = 30.0
NOISE_RATIO = float(10.0 ** (-SNR_DB / 20.0))  # 10^-1.5
NOISE_SEEDS = (0, 1, 2)
SUCCESS_OPTIMALITY = 1e-7
MSG_LIMIT = 120

ALPHA_INIT = np.array([1.0, 1.0, 0.0], dtype=float)
P_TRUE = np.array([0.08, -0.06, 0.05], dtype=float)
P_ZERO = np.zeros(3)
P_PERT_POS = np.array([0.18, 0.04, 0.15], dtype=float)
P_INIT_LABEL = "zero"

METHOD_ORDER = ("wrongpose", "direct", "reduced_r4", "reduced_r6")
BASIS_CENTERS = ((-0.15, 0.10), (0.20, -0.10), (0.00, 0.00))

SCENE_GRIDS = (
    {
        "scene_id": "base",
        "N": 24,
        "alpha_true": [1.5, 2.0, 0.0],
    },
    {
        "scene_id": "base",
        "N": 32,
        "alpha_true": [1.5, 2.0, 0.0],
    },
    {
        "scene_id": "alt",
        "N": 16,
        "alpha_true": [2.0, -1.0, 0.8],
    },
)


def scene_label(scene_id: str, N: int) -> str:
    return f"{scene_id} N{N}"


# ---------------------------------------------------------------------------
# Noise (exactly the E8/E9 convention: ||noise||/||d_true|| = 10^-1.5)
# ---------------------------------------------------------------------------
def make_noise(seed: int, M: int) -> np.ndarray:
    rng = np.random.default_rng(int(seed))
    return (
        rng.standard_normal(M) + 1j * rng.standard_normal(M)
    ) / np.sqrt(2.0)


def scaled_observation(
    d_true: np.ndarray, raw: np.ndarray, snr_db: float
) -> tuple[np.ndarray, float]:
    ratio = float(10.0 ** (-snr_db / 20.0))
    scale = ratio * float(np.linalg.norm(d_true)) / float(np.linalg.norm(raw))
    return d_true + scale * raw, ratio


# ---------------------------------------------------------------------------
# Per-run solver (E5-final fixed-V_vis reduced pattern; E8 record layout)
# ---------------------------------------------------------------------------
def run_case(
    method: str,
    seed: int,
    scene_id: str,
    alpha_true: np.ndarray,
    d_obs: np.ndarray,
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    subspaces_zero: dict,
) -> dict:
    """One scipy TRF run and the requested error/state diagnostics."""
    K = Phi.shape[1]
    p_init = P_ZERO
    if method in ("wrongpose", "direct"):
        reduced_r = None
        solver_subspace = None
        decomp_subspace = subspaces_zero["r6"]
        decomp_basis = "r6"
    else:
        reduced_r = int(method[-1])
        solver_subspace = subspaces_zero[f"r{reduced_r}"]
        decomp_subspace = solver_subspace
        decomp_basis = f"r{reduced_r}"
    V_vis = decomp_subspace["V_vis"]

    def residual(theta: np.ndarray) -> np.ndarray:
        if method == "wrongpose":
            alpha = np.asarray(theta, dtype=float)
            p = np.zeros(3)
        elif method == "direct":
            alpha = np.asarray(theta[:K], dtype=float)
            p = np.asarray(theta[K:], dtype=float)
        else:
            alpha = np.asarray(theta[:K], dtype=float)
            q = np.asarray(theta[K:], dtype=float)
            p = p_init + solver_subspace["V_vis"] @ q
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
        alpha = np.asarray(theta[:K], dtype=float)
        q = np.asarray(theta[K:], dtype=float)
        p = p_init + solver_subspace["V_vis"] @ q
        jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
        return np.hstack(
            [jb["J_alpha_real"], jb["B_real"] @ solver_subspace["V_vis"]]
        )

    if method == "wrongpose":
        x0 = np.asarray(ALPHA_INIT, dtype=float)
    elif method == "direct":
        x0 = np.concatenate(
            [np.asarray(ALPHA_INIT, dtype=float), np.asarray(p_init, float)]
        )
    else:
        x0 = np.concatenate(
            [
                np.asarray(ALPHA_INIT, dtype=float),
                np.zeros(int(solver_subspace["n_vis"])),
            ]
        )

    result = least_squares(residual, x0, jac=jacobian, **LS_KWARGS)
    status = int(result.status)
    optimality = float(result.optimality)
    success = bool((status > 0) and (optimality < SUCCESS_OPTIMALITY))

    if method == "wrongpose":
        alpha_est = np.asarray(result.x[:K], dtype=float)
        p_est = np.zeros(3)
    elif method == "direct":
        alpha_est = np.asarray(result.x[:K], dtype=float)
        p_est = np.asarray(result.x[K:], dtype=float)
    else:
        alpha_est = np.asarray(result.x[:K], dtype=float)
        q_est = np.asarray(result.x[K:], dtype=float)
        p_est = p_init + solver_subspace["V_vis"] @ q_est

    ref_norm = float(np.linalg.norm(d_obs))
    final_residual = float(
        np.linalg.norm(
            _realified_residual(
                forward(alpha_est, p_est, xs, h, k, G_D, Phi, P), d_obs
            )
        )
        / ref_norm
    )
    chi_true = Phi @ alpha_true
    chi_est = Phi @ alpha_est
    map_error = float(
        np.linalg.norm(chi_est - chi_true) / np.linalg.norm(chi_true)
    )
    delta_p = np.asarray(p_est, dtype=float) - P_TRUE
    pose_error = float(np.linalg.norm(delta_p))
    visible_proj = V_vis @ (V_vis.T @ delta_p)
    pose_error_visible = float(np.linalg.norm(visible_proj))
    pose_error_hidden = float(np.linalg.norm(delta_p - visible_proj))

    if reduced_r is not None:
        t_u_rank = reduced_r
    else:
        t_u_rank = 4  # disambiguated/default retained rank, as in E5/E8
    t_u = state_witness(
        alpha_est, p_est, t_u_rank, xs, h, k, G_D, Phi, P
    )
    message = str(result.message).strip().replace("\n", " ")
    if len(message) > MSG_LIMIT:
        message = message[: MSG_LIMIT - 3] + "..."

    return {
        "scene_id": str(scene_id),
        "N": int(P.N),
        "method": method,
        "r": int(reduced_r) if reduced_r is not None else None,
        "noise_seed": int(seed),
        "p_init_label": P_INIT_LABEL,
        "p_init": [float(v) for v in p_init],
        "success": success,
        "status": status,
        "optimality": optimality,
        "nfev": int(result.nfev),
        "message": message,
        "alpha_est": [float(v) for v in alpha_est],
        "p_est": [float(v) for v in p_est],
        "pose_error": pose_error,
        "pose_error_visible": pose_error_visible,
        "pose_error_hidden": pose_error_hidden,
        "tx_err": float(abs(delta_p[0])),
        "ty_err": float(abs(delta_p[1])),
        "theta_err": float(abs(delta_p[2])),
        "map_error": map_error,
        "final_residual": final_residual,
        "n_vis": int(decomp_subspace["n_vis"]),
        "hidden_rank": int(decomp_subspace["hidden_rank"]),
        "decomposition_basis": decomp_basis,
        "hid_sv": [float(v) for v in decomp_subspace["hid_sv"]],
        "threshold": float(decomp_subspace["threshold"]),
        "T_U": t_u,
        "T_U_rank": int(t_u_rank),
    }


# ---------------------------------------------------------------------------
# Rank / threshold checks
# ---------------------------------------------------------------------------
def rank_check(
    P: Params, xs: np.ndarray, h: float
) -> dict:
    """8 singular values of G_s at p_init=zero plus the 0-indexed gaps."""
    y, _ = receiver_positions(P_ZERO, P)
    G_s = data_matrix(y, xs, h, P.k)
    s = np.linalg.svd(G_s, compute_uv=False)
    s = np.asarray(s, dtype=float)
    assert s.shape[0] == P.M, "expected 8 singular values"
    return {
        "singular_values": [float(v) for v in s],
        "gaps_s3_minus_s4": float(s[3] - s[4]),
        "gaps_s4_minus_s5": float(s[4] - s[5]),
        "gaps_s5_minus_s6": float(s[5] - s[6]),
        "gaps_s6_minus_s7": float(s[6] - s[7]),
        "p_init": [float(v) for v in P_ZERO],
    }


def subspace_json(sp: dict) -> dict:
    return {
        "rank": int(sp["rank"]),
        "p_init_label": P_INIT_LABEL,
        "n_vis": int(sp["n_vis"]),
        "hidden_rank": int(sp["hidden_rank"]),
        "hid_sv": [float(v) for v in sp["hid_sv"]],
        "threshold": float(sp["threshold"]),
        "V_vis_columns": np.asarray(sp["V_vis"]).tolist(),
        "hidden_directions": sp["hidden_directions"],
        "hidden_direction_note": (
            "rows of Vhb (right-singular 3-vectors) for hid_sv at or below "
            "1e-8*hid_sv[0]"
        ),
    }


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def percentile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    return float(np.percentile(values, q))


def build_summary(records: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = {}
    for rec in records:
        key = (rec["scene_id"], rec["N"], rec["method"])
        groups.setdefault(key, []).append(rec)
    out = []
    for (scene_id, N, method) in sorted(groups):
        sel = groups[(scene_id, N, method)]
        pose = [r["pose_error"] for r in sel]
        out.append(
            {
                "scene_id": scene_id,
                "N": int(N),
                "method": method,
                "n_seeds": len(sel),
                "n_success": sum(int(r["success"]) for r in sel),
                "median_pose_error": percentile(pose, 50),
                "iqr_pose_error": percentile(pose, 75) - percentile(pose, 25),
                "median_theta_err": percentile(
                    [r["theta_err"] for r in sel], 50
                ),
                "iqr_theta_err": percentile(
                    [r["theta_err"] for r in sel], 75
                )
                - percentile([r["theta_err"] for r in sel], 25),
                "median_map_error": percentile(
                    [r["map_error"] for r in sel], 50
                ),
                "iqr_map_error": percentile(
                    [r["map_error"] for r in sel], 75
                )
                - percentile([r["map_error"] for r in sel], 25),
                "median_final_residual": percentile(
                    [r["final_residual"] for r in sel], 50
                ),
                "iqr_final_residual": percentile(
                    [r["final_residual"] for r in sel], 75
                )
                - percentile([r["final_residual"] for r in sel], 25),
            }
        )
    return out


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
COLORS = {
    "wrongpose": "#d62728",
    "direct": "#1f77b4",
    "reduced_r4": "#2ca02c",
    "reduced_r6": "#9467bd",
}


def _grouped_median_axes(
    ax,
    labels: list[str],
    values: dict[tuple[str, int, str], float],
    title: str,
    ylabel: str,
) -> None:
    x = np.arange(len(labels))
    width = 0.19
    for j, method in enumerate(METHOD_ORDER):
        off = (j - 1.5) * width
        vals = []
        for lab in labels:
            scene_id, Ns = lab.rsplit(" ", 1)
            N = int(Ns[1:]) if Ns.startswith("N") else int(Ns)
            v = values.get((scene_id, N, method), float("nan"))
            if np.isfinite(v):
                v = max(v, 1e-16)
            vals.append(v)
        ax.bar(
            x + off,
            vals,
            width,
            color=COLORS[method],
            alpha=0.88,
            label=method,
        )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(axis="y", which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.14))


def make_plot(
    records: list[dict],
    subspace_diag: dict[tuple, dict],
    png_path: Path,
) -> None:
    group_labels = [
        scene_label("base", 24),
        scene_label("base", 32),
        scene_label("alt", 16),
    ]
    medians = {}
    for scene_id in ("base", "alt"):
        for N in (16, 24, 32):
            for method in METHOD_ORDER:
                sel = [
                    r
                    for r in records
                    if r["scene_id"] == scene_id and r["N"] == N and r["method"] == method
                ]
                if not sel:
                    continue
                medians[(scene_id, N, method)] = {
                    "pose": float(np.median([r["pose_error"] for r in sel])),
                    "theta": float(np.median([r["theta_err"] for r in sel])),
                    "res": float(np.median([r["final_residual"] for r in sel])),
                }

    fig, axes = plt.subplots(2, 2, figsize=(15.0, 10.2))
    pose_map = {
        k: v["pose"] for k, v in medians.items()
    }
    theta_map = {
        k: v["theta"] for k, v in medians.items()
    }
    res_map = {
        k: v["res"] for k, v in medians.items()
    }
    _grouped_median_axes(
        axes[0, 0],
        group_labels,
        pose_map,
        "E9 panel 1: median pose error\n"
        r"$\|p_{\rm est}-p_{\rm true}\|$ by scene-grid and method",
        r"median pose error",
    )
    _grouped_median_axes(
        axes[0, 1],
        group_labels,
        theta_map,
        "E9 panel 2: median absolute theta error",
        r"median $|\theta_{\rm est}-\theta_{\rm true}|$",
    )
    _grouped_median_axes(
        axes[1, 0],
        group_labels,
        res_map,
        "E9 panel 3: median final normalized residual",
        r"median $\|d(\hat\alpha,\hat p)-d_{\rm obs}\|/\|d_{\rm obs}\|$",
    )

    # Panel 4: n_vis / hidden_rank text table per scene-grid and retained rank.
    ax = axes[1, 1]
    ax.axis("off")
    rows = []
    for lab in group_labels:
        scene_id, Ns = lab.rsplit(" ", 1)
        N = int(Ns[1:]) if Ns.startswith("N") else int(Ns)
        d = subspace_diag[(scene_id, N)]
        rows.append(
            (
                lab,
                f"n_vis={d['zero_r4']['n_vis']}  "
                f"hidden={d['zero_r4']['hidden_rank']}",
                f"n_vis={d['zero_r6']['n_vis']}  "
                f"hidden={d['zero_r6']['hidden_rank']}",
                f"n_vis={d['pert_r6']['n_vis']}  "
                f"hidden={d['pert_r6']['hidden_rank']}",
            )
        )
    table = ax.table(
        cellText=rows,
        colLabels=(
            "scene-grid",
            "r4 @ zero\nn_vis / hidden",
            "r6 @ zero\nn_vis / hidden",
            "r6 @ pert_pos\nn_vis / hidden",
        ),
        loc="center",
        cellLoc="center",
        colWidths=(0.22, 0.26, 0.26, 0.26),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.5)
    ax.set_title(
        "E9 panel 4: visible/hidden pose subspace size at alpha_init\n"
        "(hidden = 3 - n_vis, cutoff 1e-8*hid_sv[0])",
        fontsize=10,
    )

    fig.suptitle(
        "E9: bounded grid/scene robustness of E5 claims "
        "(single-transmitter co-moving, 30 dB, seeds 0-2, p_init=zero)",
        fontsize=12.5,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> dict:
    records: list[dict] = []
    rank_checks: dict[str, dict] = {}
    subspace_diag: dict[tuple, dict] = {}
    grid_configs: dict[str, dict] = {}
    noise_diag: list[dict] = []
    dropped_n32_note = None

    for cfg in SCENE_GRIDS:
        scene_id = cfg["scene_id"]
        N = int(cfg["N"])
        alpha_true = np.asarray(cfg["alpha_true"], dtype=float)
        gkey = scene_label(scene_id, N)
        P = Params(
            N=N,
            k=12.0,
            M=8,
            R_r=1.6,
            R_t=2.0,
            phi0=0.7,
            s=0.12,
        )
        xs, h = pixel_grid(P.N)
        G_D = green_domain_matrix(xs, h, P.k)
        Phi = basis_matrix(xs, P)
        K = int(Phi.shape[1])

        d_true = forward(
            alpha_true, P_TRUE, xs, h, P.k, G_D, Phi, P
        )

        # ---- rank / threshold checks at (alpha_init, p_init=zero) -------
        rc = rank_check(P, xs, h)
        sp_zero_r4 = visible_pose_subspace(
            ALPHA_INIT, P_ZERO, 4, xs, h, P.k, G_D, Phi, P
        )
        sp_zero_r6 = visible_pose_subspace(
            ALPHA_INIT, P_ZERO, 6, xs, h, P.k, G_D, Phi, P
        )
        sp_pert_r6 = visible_pose_subspace(
            ALPHA_INIT, P_PERT_POS, 6, xs, h, P.k, G_D, Phi, P
        )
        subspace_diag[(scene_id, N)] = {
            "zero_r4": sp_zero_r4,
            "zero_r6": sp_zero_r6,
            "pert_r6": sp_pert_r6,
        }
        rank_checks[gkey] = {
            "scene_id": scene_id,
            "N": N,
            "rank_check": rc,
            "subspaces": {
                "zero_r4": subspace_json(sp_zero_r4),
                "zero_r6": subspace_json(sp_zero_r6),
                "pert_pos_r6": subspace_json(sp_pert_r6),
            },
        }
        grid_configs[gkey] = {
            "scene_id": scene_id,
            "N": N,
            "alpha_true": [float(v) for v in alpha_true],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_true": [float(v) for v in P_TRUE],
            "p_init_zero": [float(v) for v in P_ZERO],
            "p_init_pert_pos": [float(v) for v in P_PERT_POS],
            "K": K,
            "basis_centers": [list(c) for c in BASIS_CENTERS],
            "transmitter_mode": "single, co-moving pose p=(tx,ty,theta)",
        }

        # ---- noise and observations ----
        observations: dict[int, np.ndarray] = {}
        for seed in NOISE_SEEDS:
            raw = make_noise(int(seed), P.M)
            d_obs, _ = scaled_observation(d_true, raw, SNR_DB)
            observations[int(seed)] = d_obs
            noise_diag.append(
                {
                    "scene_grid": gkey,
                    "scene_id": scene_id,
                    "N": N,
                    "noise_seed": int(seed),
                    "noise_over_true_norm": float(
                        np.linalg.norm(d_obs - d_true)
                        / np.linalg.norm(d_true)
                    ),
                    "requested_ratio": NOISE_RATIO,
                }
            )

        subspaces_zero = {"r4": sp_zero_r4, "r6": sp_zero_r6}
        print(
            f"\n[{gkey}] alpha_true={[float(v) for v in alpha_true]} "
            f"n_vis: r4={sp_zero_r4['n_vis']} "
            f"r6={sp_zero_r6['n_vis']} "
            f"r6@pert={sp_pert_r6['n_vis']}"
        )
        print(
            f"[{gkey}] G_s s = "
            + ", ".join(f"{v:.4e}" for v in rc["singular_values"])
        )
        print(
            f"[{gkey}] gaps s3-s4..s6-s7 = "
            + ", ".join(
                f"{rc[key]:.3e}"
                for key in (
                    "gaps_s3_minus_s4",
                    "gaps_s4_minus_s5",
                    "gaps_s5_minus_s6",
                    "gaps_s6_minus_s7",
                )
            )
        )

        for seed in NOISE_SEEDS:
            for method in METHOD_ORDER:
                rec = run_case(
                    method,
                    int(seed),
                    scene_id,
                    alpha_true,
                    observations[int(seed)],
                    P,
                    xs,
                    h,
                    P.k,
                    G_D,
                    Phi,
                    subspaces_zero,
                )
                records.append(rec)
                print(
                    f"  [{gkey}] {method:>11s} seed={int(seed)}: "
                    f"status={rec['status']} opt={rec['optimality']:.1e} "
                    f"nfev={rec['nfev']} pose={rec['pose_error']:.3e} "
                    f"theta={rec['theta_err']:.3e} map={rec['map_error']:.3e} "
                    f"res={rec['final_residual']:.3e} "
                    f"T_U={rec['T_U']:.3f} ok={rec['success']}"
                )

    summary = build_summary(records)
    failures = [
        r
        for r in records
        if not r["success"]
    ]

    results = {
        "experiment": "run_e9_grid_scenes",
        "goal": (
            "Bounded grid/scene robustness check of the core E5 claims: "
            "reduced_r4 lossless vs direct for identifiable pose, and "
            "reduced_r6 freezes hidden directions whenever hidden_rank>0."
        ),
        "methods": list(METHOD_ORDER),
        "parameters": {
            "k": 12.0,
            "M": 8,
            "R_r": 1.6,
            "R_t": 2.0,
            "phi0": 0.7,
            "s": 0.12,
            "basis_centers": [list(c) for c in BASIS_CENTERS],
            "SNR_dB": SNR_DB,
            "noise_ratio": NOISE_RATIO,
            "noise_seeds": [int(s) for s in NOISE_SEEDS],
            "p_init_label": P_INIT_LABEL,
            "least_squares": dict(LS_KWARGS),
            "success_definition": (
                "status > 0 and optimality < 1e-7"
            ),
            "visible_hidden_cutoff": "1e-8 * hid_sv[0]",
            "noise_formula": (
                "noise = scale * (rng.standard_normal(M) + "
                "1j*rng.standard_normal(M))/sqrt(2) with "
                "||noise||/||d_true|| = 10^(-SNR_dB/20) exactly"
            ),
        },
        "scene_grids": list(grid_configs.values()),
        "rank_and_subspace_checks_per_scene_grid": rank_checks,
        "noise_diagnostics": noise_diag,
        "seed_policy": (
            "All three seeds 0,1,2 were retained for every scene-grid "
            "including N=32; N=32 analytic-Jacobian TRF runs were fast enough "
            "that no seed dropping was required."
        ),
        "records": records,
        "summary_groups": summary,
        "failures": failures,
        "counts": {
            "n_records": len(records),
            "n_success": len(records) - len(failures),
            "n_failures": len(failures),
        },
        "dropped_n32_note": dropped_n32_note,
        "note": (
            "Reduced runs optimize (alpha, q) with p = p_zero + V_vis @ q "
            "using the fixed V_vis computed at (alpha_init, p_zero); "
            "wrongpose/direct decompose their pose error on the r=6 V_vis at "
            "p_zero (decomposition_basis='r6'), while reduced runs use their "
            "own V_vis. T_U is the receiver-rx state-consistency witness at "
            "rank 4 for wrongpose/direct/reduced_r4 and rank 6 for "
            "reduced_r6. pose_error_visible/hidden are orthogonal "
            "projection norms of p_est-p_true onto V_vis and its complement. "
            "Alt scene keeps the same geometry but a different contrast "
            "coefficient vector."
        ),
    }

    json_path = HERE / "results_e9_grid_scenes.json"
    png_path = HERE / "plot_e9_grid_scenes.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(records, subspace_diag, png_path)

    # Verification: JSON round-trip, expected record count, PNG loads.
    with open(json_path) as fh:
        loaded = json.load(fh)
    expected_records = (
        len(SCENE_GRIDS) * len(NOISE_SEEDS) * len(METHOD_ORDER)
    )
    assert len(loaded["records"]) == expected_records, len(loaded["records"])
    assert loaded["records"] == results["records"]
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    print("\n--- E9 rank/threshold checks ---")
    for gkey, rk in rank_checks.items():
        rc = rk["rank_check"]
        subs = rk["subspaces"]
        print(
            f"{gkey:>9s}: svals[3..7]="
            + ", ".join(f"{v:.4e}" for v in rc["singular_values"][3:])
        )
        print(
            f"{gkey:>9s}: gaps s3-s4={rc['gaps_s3_minus_s4']:.3e} "
            f"s4-s5={rc['gaps_s4_minus_s5']:.3e} "
            f"s5-s6={rc['gaps_s5_minus_s6']:.3e} "
            f"s6-s7={rc['gaps_s6_minus_s7']:.3e}"
        )
        for key, label in (
            ("zero_r4", "r4 @ zero"),
            ("zero_r6", "r6 @ zero"),
            ("pert_pos_r6", "r6 @ pert_pos"),
        ):
            sp = subs[key]
            print(
                f"{gkey:>9s} {label:>14s}: n_vis={sp['n_vis']} "
                f"hidden_rank={sp['hidden_rank']} hid_sv="
                + ",".join(f"{v:.2e}" for v in sp["hid_sv"])
            )

    print("\n--- E9 summary (median over seeds per scene-grid/method) ---")
    header = (
        f"{'scene-grid':>12s} {'method':>11s} {'n':>2s} {'ok':>3s} | "
        f"{'pose med':>10s} {'pose IQR':>10s} {'theta med':>10s} "
        f"{'theta IQR':>10s} | {'map med':>9s} {'res med':>9s}"
    )
    print(header)
    for g in summary:
        lab = scene_label(g["scene_id"], g["N"])
        print(
            f"{lab:>12s} {g['method']:>11s} {g['n_seeds']:2d} "
            f"{g['n_success']:3d} | "
            f"{g['median_pose_error']:10.3e} {g['iqr_pose_error']:10.3e} "
            f"{g['median_theta_err']:10.3e} {g['iqr_theta_err']:10.3e} | "
            f"{g['median_map_error']:9.3e} {g['median_final_residual']:9.3e}"
        )
    print(
        f"\nTotal records: {len(records)}; successful: "
        f"{len(records) - len(failures)}; failures: {len(failures)}"
    )
    for f in failures:
        print(
            f"  FAIL {f['scene_id']} N={f['N']} {f['method']:>11s} "
            f"seed={f['noise_seed']} status={f['status']} "
            f"opt={f['optimality']:.2e} res={f['final_residual']:.2e} "
            f"pose={f['pose_error']:.2e}"
        )
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
