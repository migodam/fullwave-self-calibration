"""E8 -- multi-seed / multi-SNR replication of the E5 self-calibration study.

Reuses the E5-final geometry (single transmitter, co-moving pose
p=(tx,ty,theta)), scene, solver options, fixed visible-pose subspaces and the
validated analytic Jacobians.  For every (SNR_dB, noise_seed) pair the same
four methods are run at p_init='zero' (8 seeds) and an initialization
robustness subset runs direct / reduced_r4 / reduced_r6 at p_init='pert_pos',
SNR 30 dB, seeds 0..3.  Records keep per-run solver status, pose-error
decomposition, map error, residual, and the T_U state-consistency witness so
failures can be diagnosed individually rather than being hidden behind
medians.

Noise convention (per task spec): raw complex Gaussian noise with unit
expected norm is scaled so that ||noise||/||d_true|| = 10^(-SNR_dB/20)
exactly for each realization.
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


HERE = Path(__file__).resolve().parent

ALPHA_INIT = np.array([1.0, 1.0, 0.0], dtype=float)
ALPHA_TRUE = np.array([1.5, 2.0, 0.0], dtype=float)
P_TRUE = np.array([0.08, -0.06, 0.05], dtype=float)

SNR_DBS = (30.0, 20.0)
SEEDS_ALL = (0, 1, 2, 3, 4, 5, 6, 7)
SEEDS_PERT = (0, 1, 2, 3)
METHOD_ORDER = ("wrongpose", "direct", "reduced_r4", "reduced_r6")
P_INITS = {
    "zero": np.zeros(3),
    "pert_pos": np.array([0.18, 0.04, 0.15], dtype=float),
}
SUCCESS_OPTIMALITY = 1e-7
MSG_LIMIT = 120


# ---------------------------------------------------------------------------
# Noise generation
# ---------------------------------------------------------------------------
def make_noise(seed: int, M: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = (
        rng.standard_normal(M) + 1j * rng.standard_normal(M)
    ) / np.sqrt(2.0)
    return raw


def scaled_observation(
    d_true: np.ndarray, raw: np.ndarray, snr_db: float
) -> tuple[np.ndarray, float]:
    """d_obs = d_true + noise with ||noise||/||d_true|| = 10^(-SNR/20)."""
    ratio = float(10.0 ** (-snr_db / 20.0))
    scale = ratio * float(np.linalg.norm(d_true)) / float(np.linalg.norm(raw))
    noise = scale * raw
    return d_true + noise, ratio


# ---------------------------------------------------------------------------
# Per-run solver
# ---------------------------------------------------------------------------
def run_case(
    method: str,
    snr_db: float,
    seed: int,
    p_init_label: str,
    p_init: np.ndarray,
    d_obs: np.ndarray,
    P: Params,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    subspaces: dict,
) -> dict:
    """One scipy TRF run; record keeps full solver + error diagnostics."""
    K = Phi.shape[1]
    if method in ("wrongpose", "direct"):
        reduced_r = None
        solver_subspace = None
        decomp_subspace = subspaces[p_init_label]["r6"]
    else:
        reduced_r = int(method[-1])
        solver_subspace = subspaces[p_init_label][f"r{reduced_r}"]
        decomp_subspace = solver_subspace
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
    success = bool(
        (status in (1, 2, 3, 4))
        or ((status > 0) and (optimality < SUCCESS_OPTIMALITY))
    )

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
    chi_true = Phi @ ALPHA_TRUE
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
        decomp_basis = f"r{reduced_r}"
    else:
        t_u_rank = 4
        decomp_basis = "r6"
    t_u = state_witness(
        alpha_est, p_est, t_u_rank, xs, h, k, G_D, Phi, P
    )
    message = str(result.message).strip().replace("\n", " ")
    if len(message) > MSG_LIMIT:
        message = message[: MSG_LIMIT - 3] + "..."

    return {
        "method": method,
        "SNR_dB": float(snr_db),
        "noise_seed": int(seed),
        "p_init_label": p_init_label,
        "p_init": [float(v) for v in p_init],
        "r": int(reduced_r) if reduced_r is not None else None,
        "n_vis": int(decomp_subspace["n_vis"]),
        "hidden_rank": int(decomp_subspace["hidden_rank"]),
        "decomposition_basis": decomp_basis,
        "success": bool(success),
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
        "T_U": t_u,
        "T_U_rank": int(t_u_rank),
    }


# ---------------------------------------------------------------------------
# Summary + failures
# ---------------------------------------------------------------------------
def percentile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    return float(np.percentile(values, q))


def build_summary(records: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = {}
    for rec in records:
        key = (rec["SNR_dB"], rec["p_init_label"], rec["method"])
        groups.setdefault(key, []).append(rec)
    out = []
    for (snr, label, method) in sorted(groups):
        sel = groups[(snr, label, method)]
        pose = [r["pose_error"] for r in sel]
        out.append(
            {
                "SNR_dB": float(snr),
                "p_init_label": label,
                "method": method,
                "n_seeds": len(sel),
                "n_success": sum(int(r["success"]) for r in sel),
                "median_pose_error": percentile(pose, 50),
                "iqr_pose_error": percentile(pose, 75) - percentile(pose, 25),
                "min_pose_error": percentile(pose, 0),
                "max_pose_error": percentile(pose, 100),
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
                "median_T_U": percentile([r["T_U"] for r in sel], 50),
                "iqr_T_U": percentile([r["T_U"] for r in sel], 75)
                - percentile([r["T_U"] for r in sel], 25),
            }
        )
    return out


def failures_list(records: list[dict]) -> list[dict]:
    return [
        {
            "method": r["method"],
            "SNR_dB": r["SNR_dB"],
            "noise_seed": r["noise_seed"],
            "p_init_label": r["p_init_label"],
            "status": r["status"],
            "optimality": r["optimality"],
            "nfev": r["nfev"],
            "final_residual": r["final_residual"],
            "pose_error": r["pose_error"],
            "message": r["message"],
        }
        for r in records
        if not r["success"]
    ]


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
COLORS = {
    "wrongpose": "#d62728",
    "direct": "#1f77b4",
    "reduced_r4": "#2ca02c",
    "reduced_r6": "#9467bd",
}


def _strip(ax, x: float, values: list[float], color: str):
    """Jittered per-seed points plus horizontal median bar."""
    rng = np.random.default_rng(1234)
    xs = x + (rng.random(len(values)) - 0.5) * 0.22
    ax.plot(xs, values, "o", ms=4.5, color=color, alpha=0.85, zorder=3)
    median = float(np.median(values))
    ax.plot(
        [x - 0.28, x + 0.28],
        [median, median],
        "-",
        color="black",
        lw=2.2,
        zorder=4,
    )
    ax.plot([x], [median], "o", ms=6, color=color, mec="black", zorder=5)
    return median


def make_plot(records: list[dict], png_path: Path) -> None:
    fig = plt.figure(figsize=(15.0, 11.6))
    grid = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.20)

    # Panels 1 and 2: per-seed pose error (zero-init full seed set) SNR 30/20.
    for panel, snr in ((1, 30.0), (2, 20.0)):
        ax = fig.add_subplot(grid[(panel - 1) // 2, (panel - 1) % 2])
        x = np.arange(len(METHOD_ORDER))
        medians = {}
        for j, method in enumerate(METHOD_ORDER):
            sel = [
                r["pose_error"]
                for r in records
                if r["method"] == method
                and r["p_init_label"] == "zero"
                and r["SNR_dB"] == snr
            ]
            medians[method] = _strip(ax, x[j], sel, COLORS[method])
            ok = sum(
                1
                for r in records
                if r["method"] == method
                and r["p_init_label"] == "zero"
                and r["SNR_dB"] == snr
                and r["success"]
            )
            ax.text(
                x[j],
                min(sel) * 0.32,
                f"med {medians[method]:.2e}",
                ha="center",
                fontsize=7,
                color=COLORS[method],
            )
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(list(METHOD_ORDER), fontsize=9)
        ax.set_ylabel(r"pose error $\|p_{\rm est}-p_{\rm true}\|$")
        ax.set_title(
            f"E8 panel {panel}: per-seed pose error, "
            f"zero init, SNR {int(snr)} dB\n"
            "points per seed, black median bar"
        )
        ax.grid(axis="y", which="both", alpha=0.3)
        ax.margins(x=0.08)

    # Panel 3: per-seed theta error, SNR 30 and 20 grouped per method.
    ax = fig.add_subplot(grid[1, 0])
    xs = np.arange(len(METHOD_ORDER)) * 2.2
    snr_styles = {
        30.0: ("o", "#0b5394", "SNR 30 dB"),
        20.0: ("s", "#e07b00", "SNR 20 dB"),
    }
    rng = np.random.default_rng(42)
    medians = {}
    for snr, (marker, color, lab) in snr_styles.items():
        vals_by_method: dict[str, list[float]] = {}
        for method in METHOD_ORDER:
            vals_by_method[method] = [
                r["theta_err"]
                for r in records
                if r["method"] == method
                and r["p_init_label"] == "zero"
                and r["SNR_dB"] == snr
            ]
        for j, method in enumerate(METHOD_ORDER):
            xm = xs[j] + (-0.32 if snr == 30.0 else 0.32)
            values = vals_by_method[method]
            xpts = xm + (rng.random(len(values)) - 0.5) * 0.18
            ax.plot(xpts, values, marker, ms=4.5, color=color, alpha=0.85)
            med = float(np.median(values))
            medians[(snr, method)] = med
            ax.plot(
                [xm - 0.22, xm + 0.22], [med, med], "-", color="black", lw=2
            )
        ax.plot([], [], marker, color=color, label=lab)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(list(METHOD_ORDER), fontsize=9)
    ax.set_ylabel(r"$|\theta_{\rm est}-\theta_{\rm true}|$")
    ax.set_title(
        "E8 panel 3: per-seed theta error by SNR\n"
        "(zero init; black bars are medians)"
    )
    ax.grid(axis="y", which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")

    # Panel 4: residual vs pose-error diagnostic scatter, failures in red X.
    ax = fig.add_subplot(grid[1, 1])
    for method in METHOD_ORDER:
        sel = [r for r in records if r["method"] == method]
        ax.scatter(
            [r["final_residual"] for r in sel],
            [r["pose_error"] for r in sel],
            s=26,
            color=COLORS[method],
            alpha=0.75,
            label=method,
            zorder=3,
        )
    failed_labeled = False
    for r in records:
        if not r["success"]:
            ax.scatter(
                [r["final_residual"]],
                [r["pose_error"]],
                s=150,
                marker="X",
                color="red",
                linewidths=1.2,
                edgecolors="black",
                zorder=5,
                label="failed run" if not failed_labeled else None,
            )
            failed_labeled = True
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"final residual $\|d(\hat\alpha,\hat p)-d_{\rm obs}\|/\|d_{\rm obs}\|$")
    ax.set_ylabel(r"pose error")
    ax.set_title(
        "E8 panel 4: failure diagnostics\n"
        "all runs; red X marks success=False records"
    )
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, loc="upper left", ncol=2)

    fig.suptitle(
        "E8: multi-seed / multi-SNR replication "
        "(E5 single-transmitter co-moving self-calibration)",
        fontsize=13,
        y=0.995,
    )
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> dict:
    P = Params(N=16, k=12.0, M=8, R_r=1.6, R_t=2.0, phi0=0.7, s=0.12)
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    Phi = basis_matrix(xs, P)
    K = int(Phi.shape[1])
    alpha_true = np.asarray(ALPHA_TRUE)
    p_true = np.asarray(P_TRUE)
    d_true = forward(
        alpha_true, p_true, xs, h, P.k, G_D, Phi, P
    )

    # Fixed visible-pose subspaces at (alpha_init, p_init) for r in (4, 6).
    subspaces: dict[str, dict] = {}
    for label, p_init in P_INITS.items():
        subspaces[label] = {}
        for r in (4, 6):
            sp = visible_pose_subspace(
                ALPHA_INIT, p_init, r, xs, h, P.k, G_D, Phi, P
            )
            subspaces[label][f"r{r}"] = sp
    print("Visible-pose subspaces at (alpha_init, p_init):")
    for label in P_INITS:
        for r in (4, 6):
            sp = subspaces[label][f"r{r}"]
            print(
                f"  {label:>9s} r{r}: n_vis={sp['n_vis']} "
                f"hidden_rank={sp['hidden_rank']} hid_sv="
                + ",".join(f"{v:.3e}" for v in sp["hid_sv"])
            )

    raw_noise = {
        int(seed): make_noise(int(seed), P.M) for seed in SEEDS_ALL
    }
    observations: dict[tuple[float, int], np.ndarray] = {}
    noise_ratio_actual: dict[tuple[float, int], float] = {}
    for snr in SNR_DBS:
        for seed in SEEDS_ALL:
            d_obs, ratio = scaled_observation(d_true, raw_noise[seed], snr)
            observations[(snr, seed)] = d_obs
            noise_ratio_actual[(snr, seed)] = float(
                np.linalg.norm(d_obs - d_true) / np.linalg.norm(d_true)
            )

    jobs: list[tuple[str, float, int, str]] = []
    for snr in SNR_DBS:
        for seed in SEEDS_ALL:
            for method in METHOD_ORDER:
                jobs.append((method, snr, seed, "zero"))
    for seed in SEEDS_PERT:
        for method in ("direct", "reduced_r4", "reduced_r6"):
            jobs.append((method, 30.0, seed, "pert_pos"))

    records: list[dict] = []
    print(f"\nRunning {len(jobs)} least-squares cases ...")
    for method, snr, seed, label in jobs:
        rec = run_case(
            method,
            snr,
            seed,
            label,
            P_INITS[label],
            observations[(snr, seed)],
            P,
            xs,
            h,
            P.k,
            G_D,
            Phi,
            subspaces,
        )
        records.append(rec)
        print(
            f"  {method:>11s} SNR={snr:4.1f} seed={seed} {label:>9s}: "
            f"status={rec['status']} opt={rec['optimality']:.1e} "
            f"pose={rec['pose_error']:.3e} theta={rec['theta_err']:.3e} "
            f"map={rec['map_error']:.3e} res={rec['final_residual']:.3e} "
            f"T_U={rec['T_U']:.3f} ok={rec['success']}"
        )

    summary = build_summary(records)
    failures = failures_list(records)

    # JSON-friendly subspace diagnostics per (p_init, r).
    subspace_json = {}
    for label in P_INITS:
        subspace_json[label] = {}
        for r in (4, 6):
            sp = subspaces[label][f"r{r}"]
            subspace_json[label][f"r{r}"] = {
                "rank": int(sp["rank"]),
                "n_vis": int(sp["n_vis"]),
                "hidden_rank": int(sp["hidden_rank"]),
                "hid_sv": [float(v) for v in sp["hid_sv"]],
                "threshold": float(sp["threshold"]),
                "V_vis_columns": np.asarray(sp["V_vis"]).tolist(),
                "hidden_directions": sp["hidden_directions"],
                "hidden_direction_note": (
                    "rows of Vhb (right-singular 3-vectors) for hid_sv at or "
                    "below 1e-8*hid_sv[0]"
                ),
            }

    results = {
        "experiment": "run_e5_seeds",
        "scene": {
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
            "alpha_true": [float(v) for v in alpha_true],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_true": [float(v) for v in p_true],
            "p_init_variants": {
                label: [float(v) for v in p] for label, p in P_INITS.items()
            },
            "transmitter_mode": "single, co-moving pose p=(tx,ty,theta)",
        },
        "data_noise": {
            "SNR_dB_list": list(SNR_DBS),
            "seeds": list(SEEDS_ALL),
            "noise_formula": (
                "noise = scale * (rng.standard_normal(M) + "
                "1j*rng.standard_normal(M))/sqrt(2) with "
                "||noise||/||d_true|| = 10^(-SNR_dB/20) exactly"
            ),
            "actual_ratio_per_snr_seed": [
                {
                    "SNR_dB": snr,
                    "noise_seed": seed,
                    "noise_over_true_norm": noise_ratio_actual[(snr, seed)],
                }
                for snr in SNR_DBS
                for seed in SEEDS_ALL
            ],
        },
        "least_squares": dict(LS_KWARGS),
        "success_definition": (
            "(status in {1,2,3,4}) or "
            "((status > 0) and (optimality < 1e-7))"
        ),
        "visible_pose_subspaces_at_p_init": subspace_json,
        "note": (
            "Each record: solver details, pose-error decomposition on the "
            "run's own V_vis (wrongpose/direct use the r=6 subspace at the "
            "same p_init, decomposition_basis='r6'), map error, normalized "
            "final residual, and T_U = receiver-rx state-consistency witness "
            "at rank 4 (wrongpose/direct/reduced_r4) or 6 (reduced_r6). "
            "pose_error_visible/hidden are orthogonal projections of "
            "p_est-p_true onto V_vis and its complement."
        ),
        "records": records,
        "summary_groups": summary,
        "failures": failures,
        "counts": {
            "n_records": len(records),
            "n_failures": len(failures),
            "n_success": len(records) - len(failures),
        },
    }

    json_path = HERE / "results_e5_seeds.json"
    png_path = HERE / "plot_e5_seeds.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(records, png_path)

    # Verification: JSON round-trip, record count, PNG loads.
    with open(json_path) as fh:
        loaded = json.load(fh)
    assert len(loaded["records"]) == 76, len(loaded["records"])
    assert loaded["records"] == results["records"]
    assert png_path.exists() and png_path.stat().st_size > 0
    import PIL.Image

    im = PIL.Image.open(png_path)
    im.load()

    print("\n--- E8 summary (median over seeds per group) ---")
    header = (
        f"{'SNR':>5s} {'init':>9s} {'method':>11s} "
        f"{'n':>2s} {'ok':>3s} | {'pose med':>10s} {'pose IQR':>10s} "
        f"{'theta med':>10s} {'theta IQR':>10s} | "
        f"{'map med':>9s} {'res med':>9s} {'T_U med':>8s}"
    )
    print(header)
    for g in summary:
        print(
            f"{int(g['SNR_dB']):5d} {g['p_init_label']:>9s} "
            f"{g['method']:>11s} {g['n_seeds']:2d} {g['n_success']:3d} | "
            f"{g['median_pose_error']:10.3e} {g['iqr_pose_error']:10.3e} "
            f"{g['median_theta_err']:10.3e} {g['iqr_theta_err']:10.3e} | "
            f"{g['median_map_error']:9.3e} {g['median_final_residual']:9.3e} "
            f"{g['median_T_U']:8.3f}"
        )
    print(f"\nTotal records: {len(records)}; "
          f"successful: {len(records) - len(failures)}; failures: {len(failures)}")
    for f in failures:
        print(
            f"  FAIL {f['method']:>11s} SNR={f['SNR_dB']:4.1f} "
            f"seed={f['noise_seed']} {f['p_init_label']:>9s} "
            f"status={f['status']} opt={f['optimality']:.2e} "
            f"res={f['final_residual']:.2e} pose={f['pose_error']:.2e}"
        )
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
