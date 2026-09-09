"""E10b -- completed retained-rank hiding transition for M=12 and M=16.

E10 swept retained rank r in {2,...,8} at (alpha_init, p_init) and found
hidden_rank = 0 for every r <= 8 when M=12 or M=16 (full-circle receivers,
base scene, k=12), whereas the canonical M=8 case already has hidden_rank=2
at r=6.  E10b closes that gap: it sweeps r = 2..M and runs the reduced
nonlinear solver in the newly hidden regimes.

Why the transition sits at r = M-2 (empirical rule confirmed below):
the geometry-lift hiding span is col(realify(G_s V_r), realify(J_alpha)),
a real matrix with 2r+K = 2r+3 columns inside R^(2M).  When
2r+3 = 2M-1 (r = M-2) the orthocomplement has dimension 1, so the
projected pose block B_red = P_perp @ B_real has rank <= 1: n_vis=1 and
hidden_rank=2.  For r >= M-1 the columns 2r+3 >= 2M+1 exceed the data
dimension, P_perp is numerically 0, B_red ~ 1e-18, and the relative-gap
rule (n_vis = count(hid_sv > 1e-8*hid_sv[0])) reads hidden_rank=0 again.
That trailing state is the same "empty-lift / B_red vanished" state E10
already saw for M=8 at r=7,8.

Nonlinear runs (same solver/noise conventions as E10):
  * reduced_r{r} for the first (and here only) r with hidden_rank>0,
  * a clearly flagged diagnostic reduced_r{M-1} in the empty-lift state
    right after the hidden regime (hidden_rank back to 0, B_red~0), and
  * direct ([alpha; p], full physics) for comparison.
Seeds 0,1,2, SNR 30 dB, noise scaled so ||noise||/||d_true||=10^(-30/20).

Pose-error decomposition uses the run's own V_vis_r for reduced runs and the
first hidden-r V_vis (r = M-2) for direct runs, so reduced/direct visible and
hidden components are directly comparable in the newly hidden regimes.

CPU only, deterministic, reuses E5/E9/E10 functions; does not touch existing
result files.
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
from run_e5 import (
    basis_matrix,
    forward,
    jacobians,
    state_witness,
    _realified_residual,
)
from run_e5_final import visible_pose_subspace, LS_KWARGS
from run_e10_settings_sweep import (
    BASE_N,
    BASE_P,
    ALPHA_TRUE,
    ALPHA_INIT,
    P_TRUE,
    P_INIT,
    SNR_DB,
    SUCCESS_OPTIMALITY,
    make_params,
    make_noise,
    scaled_observation,
)


HERE = Path(__file__).resolve().parent

NOISE_SEEDS = (0, 1, 2)
VANISH_SV_ABS = 1e-12  # diagnostic: B_red singular values at this absolute
# level mean P_perp ~ 0 (empty-lift state), not a genuinely visible pose.

# The two E10 settings whose r<=8 sweep showed no hidden directions.
SETTINGS = (
    dict(label="M12_k12_full", M=12, k=12.0, aperture="full"),
    dict(label="M16_k12_full", M=16, k=12.0, aperture="full"),
)


def sweep_setting(
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
) -> tuple[list[dict], dict[int, dict]]:
    """Rank-transition table for r=2..M plus the in-memory subspaces."""
    rows: list[dict] = []
    spaces: dict[int, dict] = {}
    for r in range(2, int(P.M) + 1):
        sp = visible_pose_subspace(
            ALPHA_INIT, P_INIT, int(r), xs, h, k, G_D, Phi, P
        )
        spaces[int(r)] = sp
        sv = np.asarray(sp["hid_sv"], dtype=float)
        if int(sp["hidden_rank"]) > 0:
            regime = "hidden"
        elif float(sv[0]) < VANISH_SV_ABS:
            regime = "vanished_Bred"
        else:
            regime = "visible"
        rows.append(
            {
                "r": int(r),
                "n_vis": int(sp["n_vis"]),
                "hidden_rank": int(sp["hidden_rank"]),
                "hid_sv": [float(v) for v in sv],
                "threshold": float(sp["threshold"]),
                "regime": regime,
            }
        )
    return rows, spaces


def summarize_transition(rows: list[dict]) -> dict:
    """Transition summary; generic in case more than one r is hidden."""
    hidden = [row["r"] for row in rows if int(row["hidden_rank"]) > 0]
    vanished = [row["r"] for row in rows if row["regime"] == "vanished_Bred"]
    n_vis_zero = [
        row["r"] for row in rows if int(row["n_vis"]) == 0
    ]
    summary = {
        "r_transition": int(hidden[0]) if hidden else None,
        "hidden_rs": [int(r) for r in hidden],
        "n_hidden_rs": len(hidden),
        "vanished_Bred_rs": [int(r) for r in vanished],
        "n_vis_zero_rs": [int(r) for r in n_vis_zero],
        "has_hidden_rank_3": bool(n_vis_zero),
    }
    if hidden:
        first = next(row for row in rows if int(row["r"]) == hidden[0])
        summary["hidden_rank_at_transition"] = int(first["hidden_rank"])
        summary["n_vis_at_transition"] = int(first["n_vis"])
    return summary


def run_case(
    method: str,
    r: int | None,
    seed: int,
    setting: dict,
    P: Params,
    d_obs: np.ndarray,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    solver_subspace: dict | None,
    decomp_subspace: dict,
    regime: str,
) -> dict:
    """One least-squares record (direct or reduced), E10 field conventions."""
    K = Phi.shape[1]
    V_dec = np.asarray(decomp_subspace["V_vis"], dtype=float)  # 3 x n_vis
    V_sol = (
        None
        if solver_subspace is None
        else np.asarray(solver_subspace["V_vis"], dtype=float)
    )

    def xy(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if method == "direct":
            return (
                np.asarray(theta[:K], dtype=float),
                np.asarray(theta[K:], dtype=float),
            )
        # reduced: [alpha; q], p = p_init + V_vis @ q
        alpha = np.asarray(theta[:K], dtype=float)
        q = np.asarray(theta[K:], dtype=float)
        p = P_INIT.copy()
        if V_sol is not None and V_sol.shape[1] > 0:
            p = p + V_sol @ q
        return alpha, p

    def residual(theta: np.ndarray) -> np.ndarray:
        alpha, p = xy(theta)
        return _realified_residual(
            forward(alpha, p, xs, h, k, G_D, Phi, P), d_obs
        )

    def jacobian(theta: np.ndarray) -> np.ndarray:
        alpha, p = xy(theta)
        jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
        if method == "direct":
            return np.hstack([jb["J_alpha_real"], jb["B_real"]])
        pose_cols = jb["B_real"]
        if V_sol is not None and V_sol.shape[1] > 0:
            pose_cols = pose_cols @ V_sol
        else:
            pose_cols = pose_cols[:, :0]
        return np.hstack([jb["J_alpha_real"], pose_cols])

    if method == "direct":
        x0 = np.concatenate(
            [np.asarray(ALPHA_INIT, dtype=float), np.asarray(P_INIT, dtype=float)]
        )
        n_vis_solver = None
        solver_r = None
    else:
        n_vis_solver = int(solver_subspace["n_vis"])
        solver_r = int(solver_subspace["rank"])
        x0 = np.concatenate(
            [
                np.asarray(ALPHA_INIT, dtype=float),
                np.zeros(n_vis_solver, dtype=float),
            ]
        )

    result = least_squares(residual, x0, jac=jacobian, **LS_KWARGS)
    status = int(result.status)
    optimality = float(result.optimality)
    success = bool((status > 0) and (optimality < SUCCESS_OPTIMALITY))
    njev = getattr(result, "njev", None)
    njev = int(njev) if njev is not None else None

    alpha_est, p_est = xy(np.asarray(result.x, dtype=float))
    alpha_est = np.asarray(alpha_est, dtype=float)
    p_est = np.asarray(p_est, dtype=float)

    ref_norm = float(np.linalg.norm(d_obs))
    final_residual_ratio = float(
        np.linalg.norm(
            _realified_residual(
                forward(alpha_est, p_est, xs, h, k, G_D, Phi, P), d_obs
            )
        )
        / ref_norm
    )
    chi_true = Phi @ ALPHA_TRUE
    chi_est = Phi @ alpha_est
    map_error = float(
        np.linalg.norm(chi_est - chi_true) / np.linalg.norm(chi_true)
    )
    delta_p = p_est - P_TRUE
    pose_error = float(np.linalg.norm(delta_p))
    visible_proj = V_dec @ (V_dec.T @ delta_p)
    pose_error_visible = float(np.linalg.norm(visible_proj))
    pose_error_hidden = float(np.linalg.norm(delta_p - visible_proj))

    t_u_rank = int(r) if r is not None else int(decomp_subspace["rank"])
    t_u = state_witness(
        alpha_est, p_est, t_u_rank, xs, h, k, G_D, Phi, P
    )

    return {
        "setting": str(setting["label"]),
        "M": int(P.M),
        "k": float(P.k),
        "aperture": str(setting["aperture"]),
        "method": method if method == "direct" else f"reduced_r{r}",
        "r": int(r) if r is not None else None,
        "regime": regime,
        "noise_seed": int(seed),
        "success": success,
        "status": status,
        "optimality": optimality,
        "nfev": int(result.nfev),
        "njev": njev,
        "final_residual_ratio": final_residual_ratio,
        "alpha_est": [float(v) for v in alpha_est],
        "p_est": [float(v) for v in p_est],
        "map_error": map_error,
        "pose_error": pose_error,
        "pose_error_visible": pose_error_visible,
        "pose_error_hidden": pose_error_hidden,
        "decomposition_basis": f"r{int(decomp_subspace['rank'])}",
        "n_vis": int(decomp_subspace["n_vis"]),
        "hidden_rank": int(decomp_subspace["hidden_rank"]),
        "solver_n_vis": n_vis_solver,
        "solver_hidden_rank": (
            int(solver_subspace["hidden_rank"])
            if solver_subspace is not None
            else None
        ),
        "solver_r": solver_r,
        "T_U": float(t_u),
        "T_U_rank": int(t_u_rank),
    }


def _median(records: list[dict], field: str) -> float:
    vals = np.array([float(rec[field]) for rec in records], dtype=float)
    if vals.size == 0:
        return float("nan")
    return float(np.median(vals))


def make_plot(
    results: dict,
    png_path: Path,
) -> None:
    """2x2 figure: transition, singular values, per-M pose errors."""
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 12.5))
    plt.rcParams.update({"font.size": 10})
    mcolors = {"M12_k12_full": "#1f77b4", "M16_k12_full": "#d62728"}

    # ---- (1) n_vis / hidden_rank vs r ----
    ax = axes[0, 0]
    for lab in results["per_M"]:
        rows = results["per_M"][lab]["rank_table"]
        rs = [row["r"] for row in rows]
        nv = [row["n_vis"] for row in rows]
        hk = [row["hidden_rank"] for row in rows]
        col = mcolors[lab]
        ax.plot(rs, nv, "o-", color=col, ms=5, label=f"{lab} n_vis")
        ax.plot(rs, hk, "s--", ms=4, color=col, alpha=0.65,
                label=f"{lab} hidden_rank")
        trans = results["per_M"][lab]["transition_summary"]
        if trans["r_transition"] is not None:
            rt = trans["r_transition"]
            ax.annotate(
                f"hidden at r={rt}",
                xy=(rt, 1.0),
                xytext=(rt + 0.5, 2.55),
                fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.9, color=col),
                color=col,
            )
            ax.axvline(rt, color=col, ls=":", lw=0.9, alpha=0.7)
    ax.set_xlabel("retained rank r")
    ax.set_ylabel("pose DOF")
    ax.set_ylim(-0.2, 3.2)
    ax.set_xlim(1, max(r["r"] for M in results["per_M"].values()
                       for r in M["rank_table"]) + 0.5)
    ax.set_xticks(range(2, 17, 2))
    ax.set_title(
        "(1) Visible/hidden pose DOF vs r (M=12,16 full circle)\n"
        "transition at r=M-2; r>=M-1 is empty-lift (B_red~0, hidden_rank=0)"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # ---- (2) singular values of B_red vs r (log y) ----
    ax = axes[0, 1]
    styles = ("-", "--", ":")
    for lab in results["per_M"]:
        rows = results["per_M"][lab]["rank_table"]
        rs = [row["r"] for row in rows]
        col = mcolors[lab]
        for j in range(3):
            sv = np.array([row["hid_sv"][j] for row in rows], dtype=float)
            ax.semilogy(
                rs,
                np.maximum(sv, 1e-22),
                styles[j],
                color=col,
                lw=1.4 if j == 0 else 1.0,
                label=f"{lab} sv{j+1}",
            )
    ax.set_xlabel("retained rank r")
    ax.set_ylabel(r"singular values of $B_{\rm red}=P_\perp B_{\rm real}$")
    ax.set_ylim(1e-22, 1e-1)
    ax.set_xticks(range(2, 17, 2))
    ax.set_title(
        "(2) Pose singular-value spectrum vs r\n"
        "(solid sv1; dashed sv2; dotted sv3; log y)"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2)

    # ---- (3a)/(3b) median pose errors per M ----
    for a, lab in enumerate(("M12_k12_full", "M16_k12_full")):
        ax = axes[1, a]
        M = results["per_M"][lab]["M"]
        trans = results["per_M"][lab]["transition_summary"]
        runs = results["runs"]
        supp_r = trans["vanished_Bred_rs"][0] if trans["vanished_Bred_rs"] else None
        cats = []
        rcats = {}
        for r, regime in ((trans["r_transition"], "hidden"), (supp_r, "empty-lift")):
            if r is None:
                continue
            rcats[str(r)] = (regime, r)
            cats.append((f"r{r}", regime, r))

        x = np.arange(len(cats))
        width = 0.19
        colors = {
            "red_vis": "#2ca02c",
            "red_hid": "#d62728",
            "dir_vis": "#7fbf7b",
            "dir_hid": "#fbb4ae",
        }
        legend_seen: set[str] = set()
        for xi, (_name, regime, r) in enumerate(cats):
            red = [rec for rec in runs if rec["method"] == f"reduced_r{r}"
                   and rec["setting"] == lab and rec["success"]]
            # direct record exists per seed once; decomposition basis is the
            # first hidden r (transition), as recorded in decomposition_basis.
            dirs = [rec for rec in runs if rec["method"] == "direct"
                    and rec["setting"] == lab and rec["success"]]
            values = {
                "red_vis": _median(red, "pose_error_visible"),
                "red_hid": _median(red, "pose_error_hidden"),
                "dir_vis": _median(dirs, "pose_error_visible"),
                "dir_hid": _median(dirs, "pose_error_hidden"),
            }
            offs = {"red_vis": -1.5, "red_hid": -0.5,
                    "dir_vis": 0.5, "dir_hid": 1.5}
            for key, val in values.items():
                v = max(float(val), 1e-16)
                off = offs[key] * width
                lab_s = key.replace("red_", "reduced ").replace("dir_", "direct ")
                lab_s += " visible" if key.endswith("vis") else " hidden"
                show = lab_s not in legend_seen
                if show:
                    legend_seen.add(lab_s)
                ax.bar(
                    x[xi] + off,
                    v,
                    width,
                    color=colors[key],
                    label=lab_s if show else None,
                    alpha=0.95 if key.startswith("red") else 0.75,
                    edgecolor=colors[key],
                )
            for key in ("red_vis", "red_hid"):
                ax.text(
                    x[xi] + offs[key] * width,
                    max(float(values[key]), 1e-16) * 2.0,
                    f"{values[key]:.1e}",
                    ha="center",
                    va="bottom",
                    fontsize=5.5,
                )
        ax.set_yscale("log")
        ax.set_xticks(x)
        labels = []
        for name, regime, _r in cats:
            labels.append(f"{name}\n({regime})")
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("median pose error (3 seeds)")
        ax.set_ylim(1e-16, 1e1)
        ax.set_title(
            f"(3) {lab} (M={M}) direct vs reduced\n"
            "direct decomposition basis = first hidden r (r=M-2); "
            "reduced uses its own basis"
        )
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=7, loc="upper left")

    fig.suptitle(
        "E10b: completed retained-rank hiding transition (M=12,16; k=12; "
        "full circle; 30 dB)",
        y=1.01,
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> dict:
    xs, h = pixel_grid(BASE_N)
    Phi = basis_matrix(xs, Params(N=BASE_N, s=BASE_P["s"]))
    k = 12.0
    G_D = green_domain_matrix(xs, h, k)
    K = int(Phi.shape[1])

    results: dict = {
        "experiment": "run_e10b_m12_m16_highrank",
        "parameters": {
            "N": BASE_N,
            "K": K,
            "k": k,
            "R_r": BASE_P["R_r"],
            "R_t": BASE_P["R_t"],
            "phi0": BASE_P["phi0"],
            "s": BASE_P["s"],
            "alpha_true": [float(v) for v in ALPHA_TRUE],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_true": [float(v) for v in P_TRUE],
            "p_init": [float(v) for v in P_INIT],
            "SNR_dB": SNR_DB,
            "noise_ratio": float(10.0 ** (-SNR_DB / 20.0)),
            "noise_seeds": [int(s) for s in NOISE_SEEDS],
            "success_optimality_max": SUCCESS_OPTIMALITY,
            "settings": [
                {
                    "label": s["label"],
                    "M": int(s["M"]),
                    "aperture": s["aperture"],
                    "receiver_angles": None,
                }
                for s in SETTINGS
            ],
            "least_squares": dict(LS_KWARGS),
        },
        "per_M": {},
        "runs": [],
    }

    all_runs: list[dict] = []
    per_m = {}
    for setting in SETTINGS:
        lab = setting["label"]
        P = make_params(setting)
        assert P.receiver_angles is None  # full circle
        rows, spaces = sweep_setting(P, xs, h, k, G_D, Phi)
        trans = summarize_transition(rows)

        # Nonlinear selection: first (here the only) hidden r, then the
        # post-transition empty-lift r as a flagged diagnostic.
        selected: list[tuple[int, str]] = []
        if trans["r_transition"] is not None:
            selected.append((int(trans["r_transition"]), "hidden"))
        for rr in trans["n_vis_zero_rs"]:
            if rr not in [s[0] for s in selected]:
                selected.append((int(rr), "n_vis_zero"))
        if trans["vanished_Bred_rs"]:
            supp = int(trans["vanished_Bred_rs"][0])
            if supp not in [s[0] for s in selected]:
                selected.append((supp, "empty_lift_diagnostic"))
        selected = selected[:4]

        d_true = forward(ALPHA_TRUE, P_TRUE, xs, h, k, G_D, Phi, P)
        obs = {}
        for seed in NOISE_SEEDS:
            raw = make_noise(int(seed), int(P.M))
            obs[int(seed)], _ = scaled_observation(d_true, raw)

        decomp_sp = spaces[int(trans["r_transition"])] if trans["r_transition"] else None
        for rr, regime in selected:
            solver_sp = spaces[int(rr)]
            for seed in NOISE_SEEDS:
                rec = run_case(
                    "reduced",
                    int(rr),
                    int(seed),
                    setting,
                    P,
                    obs[int(seed)],
                    xs,
                    h,
                    k,
                    G_D,
                    Phi,
                    solver_sp,
                    solver_sp,
                    regime,
                )
                all_runs.append(rec)
                print(
                    f"[E10b] {lab:>16s} {rec['method']:>12s} "
                    f"seed={int(seed)} regime={regime:>18s} "
                    f"status={rec['status']} pose={rec['pose_error']:.3e} "
                    f"vis={rec['pose_error_visible']:.3e} "
                    f"hid={rec['pose_error_hidden']:.3e} "
                    f"map={rec['map_error']:.3e}"
                )

        for seed in NOISE_SEEDS:
            rec = run_case(
                "direct",
                None,
                int(seed),
                setting,
                P,
                obs[int(seed)],
                xs,
                h,
                k,
                G_D,
                Phi,
                None,
                decomp_sp,
                "direct_full_physics",
            )
            all_runs.append(rec)
            print(
                f"[E10b] {lab:>16s} {'direct':>12s} "
                f"seed={int(seed)} regime={'direct_full_physics':>18s} "
                f"status={rec['status']} pose={rec['pose_error']:.3e} "
                f"vis={rec['pose_error_visible']:.3e} "
                f"hid={rec['pose_error_hidden']:.3e} "
                f"map={rec['map_error']:.3e}"
            )

        per_m[lab] = {
            "M": int(P.M),
            "k": float(P.k),
            "aperture": setting["aperture"],
            "rank_table": rows,
            "transition_summary": trans,
            "selected_run_ranks": [
                {"r": int(rr), "regime": regime} for rr, regime in selected
            ],
            "selected_subspaces": {
                f"r{int(rr)}": {
                    "rank": int(spaces[int(rr)]["rank"]),
                    "n_vis": int(spaces[int(rr)]["n_vis"]),
                    "hidden_rank": int(spaces[int(rr)]["hidden_rank"]),
                    "hid_sv": [float(v) for v in spaces[int(rr)]["hid_sv"]],
                    "threshold": float(spaces[int(rr)]["threshold"]),
                    "V_vis_columns": np.asarray(
                        spaces[int(rr)]["V_vis"]
                    ).tolist(),
                    "hidden_directions": spaces[int(rr)]["hidden_directions"],
                }
                for rr, _regime in selected
            },
        }

    results["per_M"] = per_m
    results["runs"] = all_runs
    results["counts"] = {
        "n_runs": len(all_runs),
        "n_converged": sum(1 for rec in all_runs if rec["success"]),
        "n_direct": sum(1 for rec in all_runs if rec["method"] == "direct"),
        "n_reduced": sum(1 for rec in all_runs if rec["method"] != "direct"),
    }
    results["note"] = (
        "Rank-transition rows sweep r=2..M at (alpha_init, p_init) with the "
        "canonical hidden_rank convention of run_e5_final/E10: n_vis = "
        "count(hid_sv > 1e-8*hid_sv[0]) on B_red = P_perp @ B_real, where "
        "P_perp projects onto the orthocomplement of col(realify(G_s V_r), "
        "realify(J_alpha)) in R^(2M). The first hidden regime occurs at "
        "r=M-2 (2r+3 = 2M-1 columns -> 1D complement -> hidden_rank=2). "
        "For r>=M-1 the realified lift span exceeds the 2M data dimension, "
        "P_perp~0, B_red~1e-18, and the relative-gap rule returns "
        "hidden_rank=0 (empty-lift state, same as M=8 at r=7,8 in E10); the "
        "reduced parameterization then covers the full pose again via a "
        "near-orthonormal V_vis. Reduced runs are selected at the first "
        "hidden r and, as a flagged diagnostic, at the first empty-lift r. "
        "Pose-error decomposition: reduced uses its own V_vis_r; direct uses "
        "the first hidden-r V_vis (r=M-2), so reduced/direct components are "
        "comparable in the hidden regime. Noise: raw=(N(0,1)+1jN(0,1))/"
        "sqrt(2), scaled so ||noise||/||d_true||=10^(-SNR_dB/20) (E10 "
        "convention). success = status>0 and optimality<1e-7. Stalls or "
        "failures are reported per record and never rewritten."
    )

    json_path = HERE / "results_e10b_m12_m16_highrank.json"
    png_path = HERE / "plot_e10b_m12_m16_highrank.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(results, png_path)

    # Verification.
    with open(json_path) as fh:
        loaded = json.load(fh)
    for lab in per_m:
        n_rows = per_m[lab]["M"] - 1
        assert len(per_m[lab]["rank_table"]) == n_rows
        ts = per_m[lab]["transition_summary"]
        assert ts["r_transition"] == per_m[lab]["M"] - 2
        assert ts["hidden_rank_at_transition"] == 2
        n_sel = len(per_m[lab]["selected_run_ranks"])
        assert n_sel == 2  # hidden r + first empty-lift r
        n_runs_setting = (1 + n_sel) * len(NOISE_SEEDS)
        setting_runs = [
            rec for rec in loaded["runs"] if rec["setting"] == lab
        ]
        assert len(setting_runs) == n_runs_setting
    assert len(loaded["runs"]) == len(results["runs"])
    assert loaded["per_M"] == results["per_M"]
    assert loaded["runs"] == results["runs"]
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    # Console summary.
    print("\n--- E10b summary ---")
    for lab in per_m:
        ts = per_m[lab]["transition_summary"]
        print(
            f"\n{lab} (M={per_m[lab]['M']}): "
            f"r_transition={ts['r_transition']} "
            f"hidden_rs={ts['hidden_rs']} "
            f"hidden_rank_at_transition={ts['hidden_rank_at_transition']} "
            f"(n_vis={ts['n_vis_at_transition']})"
        )
        print(
            f"  empty-lift (B_red~0) rs: {ts['vanished_Bred_rs']} "
            f"(hidden_rank=0 again, E10-style trailing state)"
        )
        for method_r in sorted(
            {rec["method"] for rec in results["runs"] if rec["setting"] == lab}
        ):
            recs = [
                rec for rec in results["runs"]
                if rec["setting"] == lab and rec["method"] == method_r
            ]
            ok = sum(int(rec["success"]) for rec in recs)
            reg = recs[0]["regime"] if recs else ""
            print(
                f"  {method_r:>12s} ({reg}): ok={ok}/{len(recs)} "
                f"med_pose={_median(recs, 'pose_error'):.3e} "
                f"med_vis={_median(recs, 'pose_error_visible'):.3e} "
                f"med_hid={_median(recs, 'pose_error_hidden'):.3e} "
                f"med_map={_median(recs, 'map_error'):.3e}"
            )
    c = results["counts"]
    print(
        f"\nConverged {c['n_converged']}/{c['n_runs']} runs "
        f"(direct={c['n_direct']}, reduced={c['n_reduced']})."
    )
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
