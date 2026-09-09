"""Family 1 pilot: derivative consistency of the 2D contrast-source harness.

Run (from the experiment root):
    .venv/bin/python src/family1_pilot.py

Pilot configuration (documented convention):
  N = 16 (grid cell side h_cell = 1/16), k_b = 2*pi,
  T = 6 poses on a 90-degree arc of radius 1.6 centred on the origin,
      arc angle phi = -45 deg .. +45 deg about the +x direction;
      pose orientation theta = atan2(-p_y, -p_x)  (body +x axis points at
      the origin).
  n_rx = 4 receivers at body offsets [(-0.06, 0), (0.06, 0), (0, -0.06),
      (0, 0.06)]; tx body offset (0, 0).
  Scene chi_0 = two Gaussian blobs 0.3*exp(-|r-c1|^2/(2*0.09^2)) +
      0.5*exp(-|r-c2|^2/(2*0.07^2)), c1 = (-0.15, -0.12), c2 = (0.18, 0.14).

Checks implemented:
  (a) state-equation residual, sigma_min(M)/||M||_2, max |E_tot|;
  (b) map FD: centred differences in chi vs A_t dchi, random unit dchi,
      seeds {0,1,2}, h in logspace(-4,-1,12);
  (c) pose FD: centred differences in X vs B dX, random unit dX,
      seeds {10,11,12}, same sweep;
  (d) second-order convergence slope in the smallest-h truncation region;
  (e) Born vs full-wave discrepancy table over scale s in {0.01,...,0.2}
      plus ||D_chi G_D||_2 at s = 1;
  (f) raw config, seeds, and tolerance-free numbers saved to JSON.

Outputs: results/family1_pilot_results.json, figures/family1_fd_convergence.png.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import helmholtz as hh  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CONFIG = {
    "N": 16,
    "k_b": 2.0 * np.pi,
    "T": 6,
    "n_rx": 4,
    "q": 1.0,
    "rx_offsets": [
        [-0.06, 0.0],
        [0.06, 0.0],
        [0.0, -0.06],
        [0.0, 0.06],
    ],
    "tx_offset": [0.0, 0.0],
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],  # phi span of the 90-degree arc
    "pose_theta_convention": (
        "theta = atan2(-p_y, -p_x): body +x axis points toward the origin"
    ),
    "chi_blobs": {
        "amp1": 0.3,
        "amp2": 0.5,
        "sigma1": 0.09,
        "sigma2": 0.07,
        "c1": [-0.15, -0.12],
        "c2": [0.18, 0.14],
    },
    "fd_h_min": 1e-4,
    "fd_h_max": 1e-1,
    "fd_n_h": 12,
    "map_seeds": [0, 1, 2],
    "pose_seeds": [10, 11, 12],
    "born_scales": [0.01, 0.05, 0.1, 0.2],
    "norm_guard_eps": float(np.finfo(float).eps),
}


def build_poses(cfg: dict) -> np.ndarray:
    """Poses: 90-degree arc at radius 1.6, theta pointing to the origin."""
    T = cfg["T"]
    phi = np.deg2rad(np.linspace(cfg["arc_phi_deg"][0], cfg["arc_phi_deg"][1], T))
    p = cfg["arc_radius"] * np.column_stack([np.cos(phi), np.sin(phi)])
    theta = np.arctan2(-p[:, 1], -p[:, 0])
    return np.column_stack([p, theta])


def make_chi0(points: np.ndarray, cfg: dict) -> np.ndarray:
    blobs = cfg["chi_blobs"]
    c1 = np.asarray(blobs["c1"], dtype=float)
    c2 = np.asarray(blobs["c2"], dtype=float)
    chi = blobs["amp1"] * np.exp(
        -np.sum((points - c1) ** 2, axis=1) / (2.0 * blobs["sigma1"] ** 2)
    ) + blobs["amp2"] * np.exp(
        -np.sum((points - c2) ** 2, axis=1) / (2.0 * blobs["sigma2"] ** 2)
    )
    return chi.astype(float)


def relative_fd_error(diff: np.ndarray, ref: np.ndarray, guard: float) -> float:
    return float(
        np.linalg.norm(diff, ord=2) / max(np.linalg.norm(ref, ord=2), guard)
    )


def map_fd_max_error(
    seed: int,
    step: float,
    dchi: np.ndarray,
    chi0: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    A: np.ndarray,
    cfg: dict,
    guard: float,
) -> float:
    """Max over poses of ||D_h F_t[dchi] - A_t dchi|| / max(||A_t dchi||, eps)."""
    Fp = hh.forward_measurements(
        chi0 + step * dchi, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
    )
    Fm = hh.forward_measurements(
        chi0 - step * dchi, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
    )
    fd_t = (Fp - Fm) / (2.0 * step)
    A_ref = A @ dchi
    errs_t = []
    n_rx = cfg["n_rx"]
    for t in range(cfg["T"]):
        sl = slice(t * n_rx, (t + 1) * n_rx)
        errs_t.append(relative_fd_error(fd_t[sl] - A_ref[sl], A_ref[sl], guard))
    return float(max(errs_t))


def pose_fd_error(
    seed: int,
    step: float,
    dX: np.ndarray,
    chi0: np.ndarray,
    X0_flat: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    B: np.ndarray,
    cfg: dict,
    guard: float,
) -> float:
    """||D_h F[dX] - B dX|| / max(||B dX||, eps) over all poses."""
    Xp = (X0_flat + step * dX).reshape(cfg["T"], 3)
    Xm = (X0_flat - step * dX).reshape(cfg["T"], 3)
    Fp = hh.forward_measurements(
        chi0, Xp, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
    )
    Fm = hh.forward_measurements(
        chi0, Xm, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
    )
    fd = (Fp - Fm) / (2.0 * step)
    ref = B @ dX
    return relative_fd_error(fd - ref, ref, guard)


def estimate_convergence_slope(
    hs: np.ndarray, errs: np.ndarray
) -> dict:
    """Slope of log10(err) vs log10(h) in the O(h^2) truncation regime.

    Roundoff contamination dominates the smallest h values (error rises as
    h -> 0), so the clean second-order region starts just above the h value
    with the smallest observed error.  We fit the smallest contiguous window
    of 5 samples immediately to the right of that minimum (or from the
    minimum itself if fewer than 5 samples remain).  Full curves are saved so
    the choice is auditable.
    """
    hs = np.asarray(hs, dtype=float)
    errs = np.asarray(errs, dtype=float)
    n = len(hs)
    i_min = int(np.argmin(errs))
    start = min(i_min + 1, n - 5)
    if start < 0 or start + 4 >= n:  # pragma: no cover - fallback
        start = max(0, n - 5)
    end = min(start + 4, n - 1)
    xs = np.log10(hs[start : end + 1])
    ys = np.log10(errs[start : end + 1])
    slope, intercept = np.polyfit(xs, ys, 1)
    yhat = slope * xs + intercept
    ss_res = float(np.sum((ys - yhat) ** 2))
    ss_tot = float(np.sum((ys - np.mean(ys)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": r2,
        "window_h_start": float(hs[start]),
        "window_h_end": float(hs[end]),
        "i_min_error_index": i_min,
        "h_at_min_error": float(hs[i_min]),
        "min_error": float(errs[i_min]),
    }


def _round_trip_json(obj):
    return json.dumps(obj, indent=2, default=lambda o: o.tolist())


def main() -> None:
    cfg = CONFIG
    t0 = datetime.now(timezone.utc)

    poses = build_poses(cfg)
    rx_offsets = np.asarray(cfg["rx_offsets"], dtype=float)
    tx_offset = np.asarray(cfg["tx_offset"], dtype=float)
    X0_flat = poses.reshape(-1)
    points, h_cell = hh.make_grid(cfg["N"])
    chi0 = make_chi0(points, cfg)
    S = cfg["N"] ** 2
    h_sweep = np.logspace(
        np.log10(cfg["fd_h_min"]), np.log10(cfg["fd_h_max"]), cfg["fd_n_h"]
    )
    guard = cfg["norm_guard_eps"]
    q_poses = np.asarray(poses, dtype=float)

    print("[family1] config:", {k: v for k, v in cfg.items() if k != "chi_blobs"})
    print("[family1] chi0 stats: max=%.6g mean=%.6g L2=%.6g"
          % (chi0.max(), chi0.mean(), np.linalg.norm(chi0)))

    # ---------------- (a) state checks + base Jacobians -------------------
    ops = hh.build_operators(chi0, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"])
    A, B, F_full0, diag_AB = hh.build_AB(
        chi0, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
    )
    print("[family1] A", A.shape, "B", B.shape)

    state_check = {
        "state_residual_per_pose": ops["state_residual_list"],
        "max_state_residual": ops["max_state_residual"],
        "sigma_min_M": ops["sigma_min_M"],
        "M_min_sigma": ops["M_min_sigma"],
        "M_norm": ops["M_norm"],
        "sigma_min_over_M_norm": ops["sigma_min_over_M_norm"],
        "max_abs_E_inc": ops["max_abs_E_inc"],
        "max_abs_E_tot": ops["max_abs_E_tot"],
        "max_field": ops["max_field"],
        "max_abs_E_tot_per_pose": [
            float(np.max(np.abs(Et))) for Et in ops["E_tot_list"]
        ],
    }

    # ----------------- (b) map finite differences -------------------------
    map_errs: dict[str, list[float]] = {}
    map_dchis: dict[str, np.ndarray] = {}
    for seed in cfg["map_seeds"]:
        rng = np.random.default_rng(seed)
        dchi = rng.standard_normal(S)
        dchi /= np.linalg.norm(dchi)
        map_dchis[str(seed)] = dchi
        curve = []
        for step in h_sweep:
            curve.append(
                map_fd_max_error(
                    seed, step, dchi, chi0, q_poses, rx_offsets, tx_offset,
                    A, cfg, guard,
                )
            )
        map_errs[str(seed)] = curve

    # ----------------- (c) pose finite differences ------------------------
    pose_errs: dict[str, list[float]] = {}
    pose_dXs: dict[str, np.ndarray] = {}
    for seed in cfg["pose_seeds"]:
        rng = np.random.default_rng(seed)
        dX = rng.standard_normal(3 * cfg["T"])
        dX /= np.linalg.norm(dX)
        pose_dXs[str(seed)] = dX
        curve = []
        for step in h_sweep:
            curve.append(
                pose_fd_error(
                    seed, step, dX, chi0, X0_flat, rx_offsets, tx_offset,
                    B, cfg, guard,
                )
            )
        pose_errs[str(seed)] = curve

    # Errors at the requested reporting step h_fd = 1e-3 (not on the 12-point
    # log sweep, so they are evaluated separately and kept in the JSON too).
    map_err_at_1e3 = {
        s: map_fd_max_error(
            int(s), 1e-3, map_dchis[s], chi0, q_poses, rx_offsets, tx_offset,
            A, cfg, guard,
        )
        for s in map_dchis
    }
    pose_err_at_1e3 = {
        s: pose_fd_error(
            int(s), 1e-3, pose_dXs[s], chi0, X0_flat, rx_offsets, tx_offset,
            B, cfg, guard,
        )
        for s in pose_dXs
    }
    for s, e in map_err_at_1e3.items():
        print("[family1] map seed %s: max rel err at h=1e-3 = %.3e" % (s, e))
    for s, e in pose_err_at_1e3.items():
        print("[family1] pose seed %s: rel err at h=1e-3 = %.3e" % (s, e))

    # ----------------- (d) convergence slopes -----------------------------
    slopes = {
        "map_seed_0": estimate_convergence_slope(h_sweep, map_errs["0"]),
        "pose_seed_10": estimate_convergence_slope(h_sweep, pose_errs["10"]),
    }
    # all seeds too, for diagnostics
    slopes["map_all_seeds"] = {
        s: estimate_convergence_slope(h_sweep, map_errs[s]) for s in map_errs
    }
    slopes["pose_all_seeds"] = {
        s: estimate_convergence_slope(h_sweep, pose_errs[s]) for s in pose_errs
    }

    # ----------------- (e) Born vs full wave ------------------------------
    born_table = []
    for s in cfg["born_scales"]:
        chi_s = s * chi0
        A_s, _, F_s, _ = hh.build_AB(
            chi_s, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
        )
        F_born_s, A_born_s = hh.born_forward(
            chi_s, poses, rx_offsets, tx_offset, cfg["N"], cfg["k_b"]
        )
        norm_F_full = float(np.linalg.norm(F_s, ord=2))
        norm_F_born = float(np.linalg.norm(F_born_s, ord=2))
        norm_A_full = float(np.linalg.norm(A_s, ord="fro"))
        norm_A_born = float(np.linalg.norm(A_born_s, ord="fro"))
        born_table.append(
            {
                "s": s,
                "rel_discrepancy_F": float(
                    np.linalg.norm(F_s - F_born_s, ord=2) / norm_F_full
                    if norm_F_full > 0.0 else np.nan
                ),
                "rel_discrepancy_A_Fro": float(
                    np.linalg.norm(A_s - A_born_s, ord="fro") / norm_A_full
                    if norm_A_full > 0.0 else np.nan
                ),
                "norm_F_full": norm_F_full,
                "norm_F_born": norm_F_born,
                "norm_A_full_Fro": norm_A_full,
                "norm_A_born_Fro": norm_A_born,
            }
        )
    # Born validity indicator at s = 1
    Dchi_G_D = chi0[:, None] * ops["G_D"]
    DchiG_norm2 = float(np.linalg.svd(Dchi_G_D, compute_uv=False)[0])

    # ----------------- (f) JSON + figure -----------------------------------
    results = {
        "generated_utc": t0.isoformat(),
        "runner": "src/family1_pilot.py",
        "config": {
            **cfg,
            "k_b": float(cfg["k_b"]),
            "poses": [p.tolist() for p in q_poses],
            "rx_offsets": rx_offsets.tolist(),
            "tx_offset": tx_offset.tolist(),
            "grid_h_cell": h_cell,
            "chi0_stats": {
                "min": float(chi0.min()),
                "max": float(chi0.max()),
                "mean": float(chi0.mean()),
                "l2": float(np.linalg.norm(chi0)),
                "positive_frac": float(np.mean(chi0 > 0.0)),
            },
        },
        "state_check": state_check,
        "map_fd": {
            "h_sweep": h_sweep.tolist(),
            "seeds": cfg["map_seeds"],
            "errors_max_over_poses_per_seed": map_errs,
            "errors_at_h_1e-3_max_over_poses": {
                k: v for k, v in map_err_at_1e3.items()
            },
            "norm_guard_eps": guard,
        },
        "pose_fd": {
            "h_sweep": h_sweep.tolist(),
            "seeds": cfg["pose_seeds"],
            "errors_per_seed": pose_errs,
            "errors_at_h_1e-3": {k: v for k, v in pose_err_at_1e3.items()},
            "norm_guard_eps": guard,
        },
        "convergence_slopes": slopes,
        "born_fullwave": {
            "table": born_table,
            "Dchi_G_D_spectral_norm_s1": DchiG_norm2,
            "note": (
                "rel_discrepancy_F = ||F_full - F_born||_2 / ||F_full||_2; "
                "rel_discrepancy_A = ||A_full - A_born||_F / ||A_full||_F"
            ),
        },
        "dimensions": {
            "A": [int(x) for x in A.shape],
            "B": [int(x) for x in B.shape],
            "F_full_length": int(len(F_full0)),
            "S_N2": S,
        },
        "per_pose_diagnostics": [
            {
                "pose": d["pose"].tolist(),
                "norm_A_t": d["norm_A_t"],
                "sigma_min_A_t": d["sigma_min_A_t"],
                "norm_B_t": d["norm_B_t"],
                "max_abs_E_tot_t": float(np.max(np.abs(d["E_tot_t"]))),
            }
            for d in diag_AB
        ],
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    results_path = results_dir / "family1_pilot_results.json"
    results_path.write_text(_round_trip_json(results) + "\n")

    # Figure: log-log relative FD error vs FD step, slope-2 reference.
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    colors = ["#d62728", "#1f77b4", "#2ca02c"]
    for idx, seed in enumerate(cfg["map_seeds"]):
        ax.loglog(
            h_sweep, map_errs[str(seed)], "o-", color=colors[idx],
            lw=1.2, ms=4.0, label=f"map dchi, seed {seed}",
        )
    for idx, seed in enumerate(cfg["pose_seeds"]):
        ax.loglog(
            h_sweep, pose_errs[str(seed)], "s--", color=colors[idx],
            lw=1.2, ms=4.0, label=f"pose dX, seed {seed}",
        )
    # slope-2 reference through the geometric mean point of the map seed 0 curve
    i_ref = min(len(h_sweep) - 1, int(np.argmin(map_errs["0"])) + 1)
    h_ref = h_sweep[i_ref]
    e_ref = map_errs["0"][i_ref]
    ax.loglog(
        [h_ref, h_sweep[-1]],
        [e_ref, e_ref * (h_sweep[-1] / h_ref) ** 2],
        "k:",
        lw=1.4,
        label="slope-2 reference",
    )
    ax.set_xlabel(r"central-difference step $h_{\rm fd}$")
    ax.set_ylabel("relative FD error")
    ax.set_title(
        "Family 1 pilot: analytic Jacobians vs centred finite differences\n"
        f"(N={cfg['N']}, k_b=2pi, T={cfg['T']}, n_rx={cfg['n_rx']}, two-blob chi0)"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    figure_path = figures_dir / "family1_fd_convergence.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    # ------------------------------- summary ------------------------------
    print("\n===== FAMILY 1 PILOT SUMMARY =====")
    print("max map FD rel err at h=1e-3 : %.3e"
          % max(map_err_at_1e3.values()))
    print("max pose FD rel err at h=1e-3: %.3e"
          % max(pose_err_at_1e3.values()))
    print("slope map seed 0 : %.3f (window %g..%g, R2=%.3f)"
          % (slopes["map_seed_0"]["slope"],
             slopes["map_seed_0"]["window_h_start"],
             slopes["map_seed_0"]["window_h_end"],
             slopes["map_seed_0"]["r_squared"]))
    print("slope pose seed 10: %.3f (window %g..%g, R2=%.3f)"
          % (slopes["pose_seed_10"]["slope"],
             slopes["pose_seed_10"]["window_h_start"],
             slopes["pose_seed_10"]["window_h_end"],
             slopes["pose_seed_10"]["r_squared"]))
    print("sigma_min(M)/||M|| = %.3e   (sigma_min=%.3e, ||M||=%.3e)"
          % (ops["sigma_min_over_M_norm"], ops["sigma_min_M"], ops["M_norm"]))
    print("max state residual = %.3e, max|E_tot| = %.3e"
          % (ops["max_state_residual"], ops["max_abs_E_tot"]))
    print("Born discrepancy table (s, relF, relA_Fro):")
    for row in born_table:
        print("  s=%.2f  relF=%.3e  relA=%.3e" % (
            row["s"], row["rel_discrepancy_F"], row["rel_discrepancy_A_Fro"]))
    print("||D_chi G_D||_2 at s=1 = %.3e" % DchiG_norm2)
    print("results ->", results_path)
    print("figure  ->", figure_path)


if __name__ == "__main__":
    main()
