"""Family 4c supplementary control: same-standoff trajectory geometry with
per-pose A-block-energy normalization.

Family 4D compared equal-continuous-path-length trajectories whose radii,
standoff ranges, and per-pose signal energies changed together, so its
ordering (arc90 > straight > circle360 > arc180) is a scenario counterexample
and cannot be attributed to angular coverage alone.  This supplementary
script runs the control that Family 4's parent audit requested: keep every
curved path at one standoff (R = 1.6), include one varying-standoff straight
trajectory, and equalize every pose's total realified A-block energy before
projecting onto the smooth basis.

All matrices reuse the verified Family 4 / Family 2 pipeline:
  - helmholtz.make_grid / family1_pilot.make_chi0 / family2.build_smooth_basis
  - family4.build_AB_whitened (W = None identity noise, then
    whiten_realify -> sqrt(2)*[Re; Im])
  - family2.thin_decomposition / range_basis and family4.matrix_norm_stats
  - Family 4 check-D definitions (retained_mass, confusable_mass,
    theta_min_deg, log_volume, near-null counts, rho_min) exactly.

Per-pose row blocks in the realified (A_pix_R, B_R) stacks: whiten_realify
returns 2*T*n_rx real rows ordered [Re of all complex rows; Im of all
complex rows], and complex rows are grouped pose-by-pose.  A pose's full
realified block is therefore the two slices
    re_t = t*n_rx .. (t+1)*n_rx
    im_t = T*n_rx + re_t,
which is the pose_block_indices convention already used by
src/family4_parent_controls.py.  A single contiguous slice
slice(t*n_rx, (t+1)*n_rx) would touch only the real-half rows and would not
equalize the pose's total A-block energy, so that convention is applied here.

Run (from the experiment root):
    .venv/bin/python src/family4c_trajectory_control.py

No existing src/results/notes file is modified; outputs are new family4c
artifacts.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import helmholtz as hh  # noqa: E402
import family1_pilot as family1  # noqa: E402
import family2_algebraic_spine as family2  # noqa: E402
import family4_frequency_trajectory as f4  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


CONFIG = {
    "family": "4c",
    "title": "same-standoff per-pose-energy-normalized trajectory control",
    "N": 16,
    "T": 6,
    "n_rx": 4,
    "rx_offsets": [
        [-0.06, 0.0],
        [0.06, 0.0],
        [0.0, -0.06],
        [0.0, 0.06],
    ],
    "tx_offset": [0.0, 0.0],
    "chi_blobs": {
        "amp1": 0.3,
        "amp2": 0.5,
        "sigma1": 0.09,
        "sigma2": 0.07,
        "c1": [-0.15, -0.12],
        "c2": [0.18, 0.14],
    },
    "smooth_basis": {
        "p": 24,
        "x_centers_n": 4,
        "y_centers_n": 6,
        "x_span": [-0.3, 0.3],
        "y_span": [-0.3, 0.3],
        "sigma_b": 0.16,
        "unit_columns": True,
        "note": "identical construction/config to family2/3/4",
    },
    "f": 1.0,
    "k_b": 2.0 * np.pi,
    "arc_radius_R": 1.6,
    "straight_mid_standoff": 1.6,
    "straight_half_length": np.pi * 1.6 / 2.0,
    "theta_rule": "theta = atan2(-p_y,-p_x): body +x axis points at origin",
    "trajectories": {
        "names": [
            "straight_same_mid",
            "arc90_same_standoff",
            "arc180_same_standoff",
            "arc270_same_standoff",
            "circle360_same_standoff",
        ],
        "builders": {
            "straight_same_mid": (
                "line p=(x, 1.6), x=linspace(-L/2, L/2, T), L=pi*1.6 "
                "(varying standoff; midpoint standoff 1.6)"
            ),
            "arc90_same_standoff": (
                "R=1.6, phi=linspace(-45, 45, T) degrees"
            ),
            "arc180_same_standoff": (
                "R=1.6, phi=linspace(-90, 90, T) degrees"
            ),
            "arc270_same_standoff": (
                "R=1.6, phi=linspace(-135, 135, T) degrees"
            ),
            "circle360_same_standoff": (
                "R=1.6, phi=linspace(0, 360, T, endpoint=False) degrees"
            ),
        },
        "core_hypothesis_retained": [
            "circle360_same_standoff",
            "arc180_same_standoff",
            "arc90_same_standoff",
            "straight_same_mid",
        ],
        "core_names": [
            "straight_same_mid",
            "arc90_same_standoff",
            "arc180_same_standoff",
            "circle360_same_standoff",
        ],
        "hypothesis_rule": (
            "normalized retained_mass circle360 > arc180 > arc90 > straight "
            "and normalized confusable_mass circle360 < arc180 < arc90 < "
            "straight on {straight, arc90, arc180, circle360}; violations "
            "recorded as hypothesis_holds_for_four_core=false (no forced pass)"
        ),
    },
    "normalization_rule": (
        "for pose t, let block_t be its full realified row block (both Re and "
        "Im slices, 2*n_rx rows in A_pix_R/B_R; whiten_realify stacks "
        "sqrt(2)*[Re; Im]); n_t = ||A_pix_R[block_t]||_F; scale A_pix_R and "
        "B_R block rows by the SAME scalar 1/n_t; then A_s_norm = "
        "A_pix_R_norm @ S and B_R_norm = B_R.  Each pose's total A-block "
        "energy is normalized to 1 and standoff/signal-strength scaling is "
        "removed."
    ),
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "randomness_note": "deterministic dense linear algebra; no RNG used",
}


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, np.generic):
        return x.item()
    if x is None or isinstance(x, (str, int, float, bool)):
        return x
    raise TypeError(f"not JSON serialisable: {type(x).__name__}")


def _round_trip_json(obj: dict) -> str:
    return json.dumps(_jsonable(obj), indent=2) + "\n"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trajectory_control_poses(name: str, cfg: dict) -> np.ndarray:
    """Pose array (T,3) for one family4c trajectory builder."""
    T = int(cfg["T"])
    R = float(cfg["arc_radius_R"])
    Lhalf = float(cfg["straight_half_length"])
    if name == "straight_same_mid":
        x = np.linspace(-Lhalf, Lhalf, T)
        p = np.column_stack([x, np.full(T, float(cfg["straight_mid_standoff"]))])
    elif name == "arc90_same_standoff":
        phi = np.deg2rad(np.linspace(-45.0, 45.0, T))
        p = R * np.column_stack([np.cos(phi), np.sin(phi)])
    elif name == "arc180_same_standoff":
        phi = np.deg2rad(np.linspace(-90.0, 90.0, T))
        p = R * np.column_stack([np.cos(phi), np.sin(phi)])
    elif name == "arc270_same_standoff":
        phi = np.deg2rad(np.linspace(-135.0, 135.0, T))
        p = R * np.column_stack([np.cos(phi), np.sin(phi)])
    elif name == "circle360_same_standoff":
        phi = np.deg2rad(np.linspace(0.0, 360.0, T, endpoint=False))
        p = R * np.column_stack([np.cos(phi), np.sin(phi)])
    else:
        raise ValueError(f"unknown family4c trajectory {name!r}")
    theta = np.arctan2(-p[:, 1], -p[:, 0])
    return np.column_stack([p, theta])


def continuous_design_length(name: str, cfg: dict) -> float:
    """Nominal continuous path length of the family4c design."""
    R = float(cfg["arc_radius_R"])
    if name == "straight_same_mid":
        return 2.0 * float(cfg["straight_half_length"])
    if name == "arc90_same_standoff":
        return R * np.deg2rad(90.0)
    if name == "arc180_same_standoff":
        return R * np.deg2rad(180.0)
    if name == "arc270_same_standoff":
        return R * np.deg2rad(270.0)
    if name == "circle360_same_standoff":
        return 2.0 * np.pi * R
    raise ValueError(name)


def pose_block_indices(T: int, n_rx: int, t: int) -> np.ndarray:
    """Rows of the full realified pose block (Re rows + Im rows)."""
    m_complex = T * n_rx
    re = np.arange(t * n_rx, (t + 1) * n_rx)
    return np.concatenate([re, m_complex + re])


def check_D_metrics(A_s: np.ndarray, B_R: np.ndarray) -> dict:
    """Family 4 check-D metrics, exactly as in family4 check_D_trajectories."""
    dA = family2.thin_decomposition(A_s)
    rA, QA = dA["rank"], dA["Q"]
    rB, Z, _, _ = family2.range_basis(B_R)
    Cmat = Z.T @ QA
    svC = np.linalg.svd(Cmat, compute_uv=False)
    cos2 = svC**2
    cos2_len = min(rA, rB)
    confusable = float(np.sum(cos2[:cos2_len]))
    retained = float(rA - confusable)
    theta_min = (
        float(np.arccos(min(1.0, float(np.max(svC))))) * 180.0 / np.pi
    )
    rho_all = np.sort(
        np.concatenate(
            [
                1.0 - cos2[:cos2_len],
                np.ones(max(rA - cos2_len, 0), dtype=float),
            ]
        )
    )[::-1]
    log_volume = float(np.sum(np.log(np.maximum(rho_all, 1e-300))))
    count_near = int(np.sum(rho_all < 1e-6))
    return {
        "r_A": int(rA),
        "r_B": int(rB),
        "confusable_mass": confusable,
        "retained_mass": retained,
        "theta_min_deg": theta_min,
        "log_volume": log_volume,
        "count_near_null_rho_lt_1e-6": count_near,
        "rho_min": float(rho_all[-1]),
        "cos2_desc": [float(x) for x in cos2[:cos2_len]],
    }


def path_geometry(name: str, poses: np.ndarray, cfg: dict) -> dict:
    p = poses[:, :2]
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    open_len = float(np.sum(seg))
    close_len = float(np.sum(seg) + np.linalg.norm(p[-1] - p[0]))
    ranges = np.linalg.norm(p, axis=1)
    return {
        "continuous_design_length": float(continuous_design_length(name, cfg)),
        "open_polyline_length_through_poses": open_len,
        "closed_polyline_length_if_closure_counted": close_len,
        "closure_counted_for_circle360": bool(name == "circle360_same_standoff"),
        "min_target_range": float(np.min(ranges)),
        "max_target_range": float(np.max(ranges)),
    }


def build_family4c_block(
    chi0: np.ndarray,
    poses: np.ndarray,
    cfg: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """Raw identity-noise realified blocks A_pix_R (2Tn x N^2), B_R (2Tn x 3T)."""
    return f4.build_AB_whitened(chi0, poses, cfg, float(cfg["k_b"]))


def order_names(rows: list[dict], metric: str) -> list[str]:
    vals = {r["trajectory"]: r[metric] for r in rows}
    return sorted(vals, key=lambda n: vals[n], reverse=True)


def compare_core_hypothesis(
    rows: list[dict], cfg: dict, metric_kind: str = "retained_mass"
) -> dict:
    core = list(cfg["trajectories"]["core_names"])
    core_rows = [r for r in rows if r["trajectory"] in core]
    obs_ret = order_names(core_rows, "retained_mass")
    obs_conf_asc = sorted(
        core_rows,
        key=lambda r: r["confusable_mass"],
    )
    obs_conf = [r["trajectory"] for r in obs_conf_asc]
    exp_ret = list(cfg["trajectories"]["core_hypothesis_retained"])
    exp_conf = exp_ret[::-1]
    return {
        "core_names": core,
        "observed_retained_desc": obs_ret,
        "observed_confusable_asc": obs_conf,
        "expected_retained_desc": exp_ret,
        "expected_confusable_asc": exp_conf,
        "retained_consistent": bool(obs_ret == exp_ret),
        "confusable_consistent": bool(obs_conf == exp_conf),
        "hypothesis_holds": bool(obs_ret == exp_ret and obs_conf == exp_conf),
    }


def run_experiment() -> dict:
    cfg = CONFIG
    t_utc = datetime.now(timezone.utc)
    t_start = time.perf_counter()

    points, chi0, h, S = f4.base_scene(f4.CONFIG)
    T = int(cfg["T"])
    n_rx = int(cfg["n_rx"])

    table_raw = []
    table_norm = []
    block_norm_rows = []
    verification_rows = []
    for name in cfg["trajectories"]["names"]:
        poses = trajectory_control_poses(name, cfg)
        A_pix_R, B_R = build_family4c_block(chi0, poses, cfg)
        geo = path_geometry(name, poses, cfg)

        # ---- per-pose pre/post block Fro norms ---------------------------
        A_pre = np.linalg.norm(
            [A_pix_R[pose_block_indices(T, n_rx, t)] for t in range(T)],
            axis=(1, 2),
        )
        B_pre = np.linalg.norm(
            [B_R[pose_block_indices(T, n_rx, t)] for t in range(T)],
            axis=(1, 2),
        )
        if np.any(A_pre <= 0.0):
            raise RuntimeError(f"zero per-pose A block norm for {name}")

        scales = 1.0 / A_pre
        A_norm = A_pix_R.copy()
        B_norm = B_R.copy()
        for t in range(T):
            idx = pose_block_indices(T, n_rx, t)
            A_norm[idx] *= scales[t]
            B_norm[idx] *= scales[t]

        A_post = np.linalg.norm(
            [A_norm[pose_block_indices(T, n_rx, t)] for t in range(T)],
            axis=(1, 2),
        )
        B_post = np.linalg.norm(
            [B_norm[pose_block_indices(T, n_rx, t)] for t in range(T)],
            axis=(1, 2),
        )
        dev = np.max(np.abs(A_post - 1.0))
        verification_rows.append(
            {
                "trajectory": name,
                "max_abs_deviation_A_block_fro_from_1": float(dev),
                "pass_within_1e-12": bool(dev <= 1e-12),
                "all_A_pre_positive": bool(np.all(A_pre > 0.0)),
            }
        )
        block_norm_rows.append(
            {
                "trajectory": name,
                "per_pose_A_block_fro_pre": [float(x) for x in A_pre],
                "per_pose_B_block_fro_pre": [float(x) for x in B_pre],
                "per_pose_A_block_fro_post": [float(x) for x in A_post],
                "per_pose_B_block_fro_post": [float(x) for x in B_post],
                "per_pose_scale_1_over_A_fro": [float(x) for x in scales],
            }
        )

        # ---- raw and normalized smooth matrices --------------------------
        A_s_raw = A_pix_R @ S
        A_s_norm = A_norm @ S
        metrics_raw = check_D_metrics(A_s_raw, B_R)
        metrics_norm = check_D_metrics(A_s_norm, B_norm)

        row_raw = {
            "trajectory": name,
            **geo,
            "poses": [p.tolist() for p in poses],
            **metrics_raw,
            "A_smooth_norms": f4.matrix_norm_stats(A_s_raw),
            "B_R_norms": f4.matrix_norm_stats(B_R),
        }
        row_norm = {
            "trajectory": name,
            **geo,
            "poses": [p.tolist() for p in poses],
            **metrics_norm,
            "A_smooth_norms": f4.matrix_norm_stats(A_s_norm),
            "B_R_norms": f4.matrix_norm_stats(B_norm),
        }
        table_raw.append(row_raw)
        table_norm.append(row_norm)

    # ---- orders ----------------------------------------------------------
    raw_orders = {
        "retained_mass_desc": order_names(table_raw, "retained_mass"),
        "confusable_mass_desc": order_names(table_raw, "confusable_mass"),
    }
    norm_orders = {
        "retained_mass_desc": order_names(table_norm, "retained_mass"),
        "confusable_mass_desc": order_names(table_norm, "confusable_mass"),
    }
    norm_core = compare_core_hypothesis(table_norm, cfg)
    raw_core = compare_core_hypothesis(table_raw, cfg)

    retained_flip = raw_orders["retained_mass_desc"] != norm_orders[
        "retained_mass_desc"
    ]
    confusable_flip = raw_orders["confusable_mass_desc"] != norm_orders[
        "confusable_mass_desc"
    ]
    family4_ref = json.loads(
        (_ROOT / "results" / "family4_results.json").read_text()
    )
    # Family 4 D names map onto the family4c core by dropping suffixes.
    f4_order = [
        {
            "arc90": "arc90_same_standoff",
            "straight": "straight_same_mid",
            "circle360": "circle360_same_standoff",
            "arc180": "arc180_same_standoff",
        }[n]
        for n in family4_ref["checks"]["D_equal_budget_trajectories"][
            "observed_retained_order"
        ]
    ]

    raw_retained_by_name = {r["trajectory"]: r["retained_mass"] for r in table_raw}
    norm_retained_by_name = {
        r["trajectory"]: r["retained_mass"] for r in table_norm
    }
    counterexample_discussion = {
        "family4_raw_observed_retained_order": family4_ref["checks"][
            "D_equal_budget_trajectories"
        ]["observed_retained_order"],
        "family4c_raw_observed_core_order": raw_core["observed_retained_desc"],
        "family4c_normalized_observed_core_order": norm_core[
            "observed_retained_desc"
        ],
        "family4_raw_to_family4c_raw_core_order_changed": bool(
            f4_order != raw_core["observed_retained_desc"]
        ),
        "family4_hypothesis_holds_normalized_four_core": bool(
            norm_core["hypothesis_holds"]
        ),
        "family4_raw_counterexample_status": (
            "resolved by the normalized control"
            if norm_core["hypothesis_holds"]
            else "persists as a counterexample under the normalized control"
        ),
        "note": (
            "Family 4's raw order is a scenario counterexample under equal "
            "design length with different radii; family4c raw removes the "
            "radius/standoff difference between curved paths (all R=1.6) "
            "while retaining a varying-standoff straight path and changing "
            "path length with angular coverage.  Only the normalized four-core "
            "comparison is the requested hypothesis gate."
        ),
    }

    results = {
        "schema": "family4c_trajectory_control",
        "generated_utc": t_utc.isoformat(),
        "runner": "src/family4c_trajectory_control.py",
        "command": ".venv/bin/python src/family4c_trajectory_control.py",
        "family": "4c",
        "title": cfg["title"],
        "normalization_rule": cfg["normalization_rule"],
        "pose_block_row_rule": (
            "whiten_realify stacks 2*T*n_rx rows as [Re(A_hat) (T*n_rx); "
            "Im(A_hat) (T*n_rx)] with complex rows grouped by pose; pose t "
            "block rows = pose_block_indices(T,n_rx,t) = "
            "concat(arange(t*n_rx,(t+1)*n_rx), arange(T*n_rx+t*n_rx,"
            "T*n_rx+(t+1)*n_rx)), matching src/family4_parent_controls.py. "
            "A contiguous slice(t*n_rx,(t+1)*n_rx) alone would cover only "
            "the real-half rows and was therefore not used."
        ),
        "runtime_seconds": float(time.perf_counter() - t_start),
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
        },
        "source_sha256": {
            p: sha256_file(_ROOT / p)
            for p in (
                "src/helmholtz.py",
                "src/family1_pilot.py",
                "src/family2_algebraic_spine.py",
                "src/family4_frequency_trajectory.py",
                "src/family4c_trajectory_control.py",
            )
        },
        "reference_sha256": {
            "results/family4_results.json": sha256_file(
                _ROOT / "results" / "family4_results.json"
            ),
        },
        "config": {
            "N": cfg["N"],
            "T": cfg["T"],
            "n_rx": cfg["n_rx"],
            "f": cfg["f"],
            "k_b": float(cfg["k_b"]),
            "arc_radius_R": cfg["arc_radius_R"],
            "straight_mid_standoff": cfg["straight_mid_standoff"],
            "chi_blobs": cfg["chi_blobs"],
            "smooth_basis": cfg["smooth_basis"],
            "trajectories": cfg["trajectories"],
            "rank_tol_rule": cfg["rank_tol_rule"],
            "grid_h_cell": float(h),
            "chi0_stats": {
                "min": float(chi0.min()),
                "max": float(chi0.max()),
                "mean": float(chi0.mean()),
                "l2": float(np.linalg.norm(chi0)),
            },
        },
        "family4_raw_reference": {
            "observed_retained_order": family4_ref["checks"][
                "D_equal_budget_trajectories"
            ]["observed_retained_order"],
            "hypothesis_holds": family4_ref["checks"][
                "D_equal_budget_trajectories"
            ]["hypothesis_holds"],
            "rows": family4_ref["checks"]["D_equal_budget_trajectories"][
                "trajectories"
            ],
        },
        "checks": {
            "1_normalization_verification": {
                "rows": verification_rows,
                "max_deviation_over_all_trajectories": float(
                    max(r["max_abs_deviation_A_block_fro_from_1"]
                        for r in verification_rows)
                ),
                "pass_all_post_pose_A_block_fro_1_within_1e-12": bool(
                    all(r["pass_within_1e-12"] for r in verification_rows)
                ),
            },
            "2_hypothesis_core_normalized": {
                "expected_retained_desc": norm_core["expected_retained_desc"],
                "expected_confusable_asc": norm_core["expected_confusable_asc"],
                "observed_retained_desc": norm_core["observed_retained_desc"],
                "observed_confusable_asc": norm_core["observed_confusable_asc"],
                "retained_consistent": norm_core["retained_consistent"],
                "confusable_consistent": norm_core["confusable_consistent"],
                "hypothesis_holds_for_four_core": norm_core["hypothesis_holds"],
            },
            "2b_hypothesis_core_raw_for_context": {
                "observed_retained_desc": raw_core["observed_retained_desc"],
                "observed_confusable_asc": raw_core["observed_confusable_asc"],
                "hypothesis_holds_for_four_core": raw_core["hypothesis_holds"],
            },
            "3_order_flips_raw_vs_normalized": {
                "raw_retained_desc": raw_orders["retained_mass_desc"],
                "normalized_retained_desc": norm_orders["retained_mass_desc"],
                "retained_order_flipped": bool(retained_flip),
                "raw_confusable_desc": raw_orders["confusable_mass_desc"],
                "normalized_confusable_desc": norm_orders[
                    "confusable_mass_desc"
                ],
                "confusable_order_flipped": bool(confusable_flip),
            },
            "3b_prior_counterexample_under_family4c_controls": (
                counterexample_discussion
            ),
        },
        "tables": {
            "raw": {"rows": table_raw},
            "normalized": {"rows": table_norm},
        },
        "per_pose_block_norms": {"rows": block_norm_rows},
        "raw_retained_mass_by_trajectory": {
            k: float(v) for k, v in raw_retained_by_name.items()
        },
        "normalized_retained_mass_by_trajectory": {
            k: float(v) for k, v in norm_retained_by_name.items()
        },
        "cannot_establish": [
            (
                "All claims are finite-dimensional and model-specific (N=16, "
                "p=24 smooth basis, T=6, R=1.6, this two-blob scene); no "
                "continuum-limit, universal-trajectory, or estimator claim."
            ),
            (
                "Per-pose normalization is a declared gain/whitening control "
                "on the linearized blocks; it does not prescribe a physical "
                "noise covariance and is not a claim about optimal SNR "
                "whitening under correlated noise."
            ),
            (
                "Same standoff fixes curved-path range but changes "
                "measurement count/budget with angular coverage (continuous "
                "design length and sampled polyline lengths differ), and the "
                "straight path still varies standoff, so even the normalized "
                "five-trajectory set is not a pure angular-coverage isolation."
            ),
            (
                "No universal trajectory-ordering theorem is asserted; all "
                "ordering statements are observations on this discrete model "
                "and normalization/geometry convention."
            ),
        ],
    }
    return results


def write_figure(results: dict, figures_dir: Path) -> Path:
    raw_rows = results["tables"]["raw"]["rows"]
    norm_rows = results["tables"]["normalized"]["rows"]
    traj = [r["trajectory"] for r in norm_rows]
    labels = [t.replace("_same_standoff", "").replace("_same_mid", "") for t in traj]
    norm_ret = [r["retained_mass"] for r in norm_rows]
    norm_conf = [r["confusable_mass"] for r in norm_rows]
    raw_ret = [r["retained_mass"] for r in raw_rows]

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.6))
    xpos = np.arange(len(traj))
    width = 0.36
    axes[0].bar(
        xpos - width / 2,
        norm_ret,
        width,
        label="normalized retained_mass",
        color="#1f77b4",
    )
    axes[0].bar(
        xpos + width / 2,
        norm_conf,
        width,
        label="normalized confusable_mass",
        color="#d62728",
    )
    axes[0].set_xticks(xpos, labels)
    axes[0].set_ylabel("mass (count of singular directions)")
    axes[0].set_title(
        "Family 4c: normalized retained / confusable mass\n"
        "per-pose A-block energy normalized to 1"
    )
    axes[0].legend(fontsize=8)
    axes[0].grid(True, axis="y", alpha=0.3)

    axes[1].bar(
        xpos - width / 2,
        raw_ret,
        width,
        label="raw (same-standoff geometry)",
        color="#ff7f0e",
    )
    axes[1].bar(
        xpos + width / 2,
        norm_ret,
        width,
        label="normalized retained_mass",
        color="#1f77b4",
    )
    axes[1].set_xticks(xpos, labels)
    axes[1].set_ylabel("retained_mass")
    axes[1].set_title("Family 4c: raw vs normalized retained_mass")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, axis="y", alpha=0.3)

    holds = results["checks"]["2_hypothesis_core_normalized"][
        "hypothesis_holds_for_four_core"
    ]
    fig.suptitle(
        "Family 4c trajectory control "
        f"(R=1.6 arcs/circle; four-core hypothesis holds = {holds})",
        fontsize=11,
    )
    fig.tight_layout()
    p = figures_dir / "family4c_trajectory_control.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


def _fmt(x: float, spec: str = ".6f") -> str:
    return f"{x:{spec}}"


def write_report(results: dict, notes_dir: Path) -> Path:
    raw_rows = results["tables"]["raw"]["rows"]
    norm_rows = results["tables"]["normalized"]["rows"]
    norm_checks = results["checks"]["2_hypothesis_core_normalized"]
    raw_checks = results["checks"]["2b_hypothesis_core_raw_for_context"]
    flips = results["checks"]["3_order_flips_raw_vs_normalized"]
    cex = results["checks"]["3b_prior_counterexample_under_family4c_controls"]
    ver = results["checks"]["1_normalization_verification"]
    bnorm = results["per_pose_block_norms"]["rows"]

    def metric_table(rows, label: str) -> str:
        lines = [
            "| trajectory | retained | confusable | theta_min deg | log_volume | rho<1e-6 | rho_min | r_A | r_B | design L | open polyline | min range | max range |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in rows:
            lines.append(
                f"| {r['trajectory']} | {_fmt(r['retained_mass'], '.6f')} | "
                f"{_fmt(r['confusable_mass'], '.6f')} | "
                f"{_fmt(r['theta_min_deg'], '.6e')} | "
                f"{_fmt(r['log_volume'], '.4f')} | "
                f"{r['count_near_null_rho_lt_1e-6']} | "
                f"{_fmt(r['rho_min'], '.4e')} | {r['r_A']} | {r['r_B']} | "
                f"{_fmt(r['continuous_design_length'], '.6f')} | "
                f"{_fmt(r['open_polyline_length_through_poses'], '.6f')} | "
                f"{_fmt(r['min_target_range'], '.6f')} | "
                f"{_fmt(r['max_target_range'], '.6f')} |"
            )
        return f"\n#### {label}\n\n" + "\n".join(lines)

    norm_table = metric_table(norm_rows, "Normalized table")
    raw_table = metric_table(raw_rows, "Raw table (same-standoff geometry, un-normalized)")

    ver_rows = "\n".join(
        f"| {r['trajectory']} | {_fmt(r['max_abs_deviation_A_block_fro_from_1'], '.3e')} | "
        f"{'PASS' if r['pass_within_1e-12'] else 'FAIL'} |"
        for r in ver["rows"]
    )

    norm_block_rows = "\n".join(
        f"| {r['trajectory']} | "
        f"{', '.join(_fmt(x, '.6e') for x in r['per_pose_A_block_fro_pre'])} | "
        f"{', '.join(_fmt(x, '.6e') for x in r['per_pose_B_block_fro_pre'])} | "
        f"{', '.join(_fmt(x, '.6e') for x in r['per_pose_A_block_fro_post'])} | "
        f"{', '.join(_fmt(x, '.6e') for x in r['per_pose_B_block_fro_post'])} | "
        f"{', '.join(_fmt(x, '.6f') for x in r['per_pose_scale_1_over_A_fro'])} |"
        for r in bnorm
    )

    geom_rows = "\n".join(
        f"| {r['trajectory']} | {_fmt(r['continuous_design_length'], '.6f')} | "
        f"{_fmt(r['open_polyline_length_through_poses'], '.6f')} | "
        f"{_fmt(r['closed_polyline_length_if_closure_counted'], '.6f')} | "
        f"{r['closure_counted_for_circle360']} | "
        f"{_fmt(r['min_target_range'], '.6f')} | "
        f"{_fmt(r['max_target_range'], '.6f')} |"
        for r in norm_rows
    )

    source_lines = "\n".join(
        f"- `{p}`: `{h}`" for p, h in results["source_sha256"].items()
    )
    norm_rows_by_name = {r["trajectory"]: r for r in norm_rows}
    raw_rows_by_name = {r["trajectory"]: r for r in raw_rows}
    A_norm_rows = "\n".join(
        f"| {r['trajectory']} | "
        f"{_fmt(raw_rows_by_name[r['trajectory']]['A_smooth_norms']['frobenius'], '.6e')} | "
        f"{_fmt(r['A_smooth_norms']['frobenius'], '.6e')} | "
        f"{_fmt(raw_rows_by_name[r['trajectory']]['B_R_norms']['frobenius'], '.6e')} | "
        f"{_fmt(r['B_R_norms']['frobenius'], '.6e')} | "
        f"{_fmt(r['A_smooth_norms']['spectral_2'], '.6e')} | "
        f"{_fmt(r['B_R_norms']['spectral_2'], '.6e')} |"
        for r in norm_rows
    )

    report = f"""# Family 4c: same-standoff, per-pose-energy-normalized trajectory control

Date: 2026-09-03 (SGT; UTC stamp in
`results/family4c_trajectory_control.json`).  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Supplementary control requested by the Family 4 parent audit rule 7: hold
curved-path standoff fixed (R = 1.6) and equalize each pose's realified
A-block energy before projection, so angular coverage is compared after
removing standoff/radius and per-pose signal-strength scaling.  **No existing
src/results/notes file was modified.**

## Exact command and runtime

```bash
.venv/bin/python src/family4c_trajectory_control.py
```

Wall runtime: {results['runtime_seconds']:.3f} s.  Platform:
{results['platform']['platform']}, machine {results['platform']['machine']},
Python {results['platform']['python']}, numpy {results['platform']['numpy']},
scipy {results['platform']['scipy']}, matplotlib
{results['platform']['matplotlib']}.  Deterministic dense linear algebra; no
RNG used.

## Source SHA-256

{source_lines}

## Scenario and normalization rule

Same scene as Family 4/6: N=16, T=6, n_rx=4, receivers at the four body
offsets, tx at the body origin, two-blob chi0 (amp 0.3/0.5, sigma 0.09/0.07
at (-0.15,-0.12) and (0.18,0.14)), p=24 unit-column smooth Gaussian-RBF basis
(4x6 centres on [-0.3,0.3]^2, sigma_b=0.16), f=1.0 with k_b=2*pi.

{results['normalization_rule']}

{results['pose_block_row_rule']}

### Trajectory builders

| trajectory | builder | design length | min/max range |
|---|---|---:|---:|
"""

    traj_lines = "\n".join(
        f"| {t} | {CONFIG['trajectories']['builders'][t]} | "
        f"{_fmt(raw_rows_by_name[t]['continuous_design_length'], '.6f')} | "
        f"{_fmt(raw_rows_by_name[t]['min_target_range'], '.6f')} / "
        f"{_fmt(raw_rows_by_name[t]['max_target_range'], '.6f')} |"
        for t in CONFIG["trajectories"]["names"]
    )

    report += traj_lines + f"""

## 1. Normalization verification

For every normalized trajectory, each pose's post-normalization total
A-block Frobenius norm must be 1 within 1e-12.

| trajectory | max |1 - ||A_block||_F| | pass <= 1e-12 |
|---|---:|---|
{ver_rows}

Max deviation over all trajectories:
{ver['max_deviation_over_all_trajectories']:.3e}; pass =
{ver['pass_all_post_pose_A_block_fro_1_within_1e-12']}.

### Per-pose A/B block Frobenius norms (raw pre, normalized post) and scales

| trajectory | A_fro pre (per pose) | B_fro pre | A_fro post | B_fro post | scale = 1/A_fro |
|---|---|---|---|---|---|
{norm_block_rows}

## 2. Hypothesis check on the four core trajectories (normalized)

Expected retained order (core): circle360 > arc180 > arc90 > straight.
Expected confusable order (ascending): straight < arc90 < arc180 < circle360.

Observed retained desc: {', '.join(norm_checks['observed_retained_desc'])}.
Observed confusable asc: {', '.join(norm_checks['observed_confusable_asc'])}.
retained_consistent = {norm_checks['retained_consistent']};
confusable_consistent = {norm_checks['confusable_consistent']}.

**hypothesis_holds_for_four_core = {norm_checks['hypothesis_holds_for_four_core']}**
(not forced to pass).

Raw four-core context: retained desc = {', '.join(raw_checks['observed_retained_desc'])};
hypothesis_holds_for_four_core (raw) =
{raw_checks['hypothesis_holds_for_four_core']}.

## 3. Ordering flips and the previous counterexample

Raw full retained order: {', '.join(flips['raw_retained_desc'])}.
Normalized full retained order: {', '.join(flips['normalized_retained_desc'])}.
retained_order_flipped = {flips['retained_order_flipped']}.

Raw full confusable desc: {', '.join(flips['raw_confusable_desc'])}.
Normalized full confusable desc:
{', '.join(flips['normalized_confusable_desc'])}.
confusable_order_flipped = {flips['confusable_order_flipped']}.

Family 4 raw counterexample retained order:
{', '.join(cex['family4_raw_observed_retained_order'])}.  Family 4c raw core
order: {', '.join(cex['family4c_raw_observed_core_order'])}.  Family 4c
normalized core order:
{', '.join(cex['family4c_normalized_observed_core_order'])}.
Family-4-to-family4c-raw core change =
{cex['family4_raw_to_family4c_raw_core_order_changed']}; normalized gate =
**{cex['family4_hypothesis_holds_normalized_four_core']}**;
status: {cex['family4_raw_counterexample_status']}.

{cex['note']}

## Raw and normalized metric tables

{raw_table}

{norm_table}

## Path geometry record (same in both tables)

| trajectory | continuous design L | open polyline | closed polyline if counted | closure counted for circle | min range | max range |
|---|---:|---:|---:|:---:|---:|---:|
{geom_rows}

## Smooth A and B norm record (matrix_norm_stats)

| trajectory | ||A_s||_F raw | ||A_s_norm||_F | ||B_R||_F raw | ||B_R_norm||_F | ||A_s_norm||_2 | ||B_R_norm||_2 |
|---|---:|---:|---:|---:|---:|---:|
{A_norm_rows}

## Cannot-establish section

"""
    for i, item in enumerate(results["cannot_establish"], start=1):
        report += f"{i}. {item}\n"

    report += f"""
## Artifacts

- results: `results/family4c_trajectory_control.json`
- figure: `figures/family4c_trajectory_control.png`
- this report: `notes/family4c_trajectory_control.md`
"""
    p = notes_dir / "family4c_trajectory_control.md"
    p.write_text(report)
    return p


def main() -> None:
    results = run_experiment()
    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    fig_path = write_figure(results, figures_dir)
    results["figure_sha256"] = {
        "figures/family4c_trajectory_control.png": sha256_file(fig_path)
    }
    results["artifacts"] = {
        "results_json": "results/family4c_trajectory_control.json",
        "figure": "figures/family4c_trajectory_control.png",
        "report_md": "notes/family4c_trajectory_control.md",
    }

    results_path = results_dir / "family4c_trajectory_control.json"
    results_path.write_text(_round_trip_json(results))
    report_path = write_report(results, notes_dir)

    v = results["checks"]["1_normalization_verification"]
    h = results["checks"]["2_hypothesis_core_normalized"]
    f = results["checks"]["3_order_flips_raw_vs_normalized"]
    cex = results["checks"]["3b_prior_counterexample_under_family4c_controls"]
    print("\n===== FAMILY 4C TRAJECTORY CONTROL SUMMARY =====")
    print(
        "normalization max deviation: "
        f"{v['max_deviation_over_all_trajectories']:.3e}; pass = "
        f"{v['pass_all_post_pose_A_block_fro_1_within_1e-12']}"
    )
    print("raw retained order: " + ", ".join(f["raw_retained_desc"]))
    print("norm retained order: " + ", ".join(f["normalized_retained_desc"]))
    print(
        "four-core normalized hypothesis holds = "
        f"{h['hypothesis_holds_for_four_core']}"
    )
    print("prior counterexample status: " + cex["family4_raw_counterexample_status"])
    for r in results["tables"]["normalized"]["rows"]:
        print(
            f"  {r['trajectory']}: retained={r['retained_mass']:.6f} "
            f"confusable={r['confusable_mass']:.6f} "
            f"theta={r['theta_min_deg']:.6e}"
        )
    print("results ->", results_path)
    print("figure  ->", fig_path)
    print("report  ->", report_path)


if __name__ == "__main__":
    main()
