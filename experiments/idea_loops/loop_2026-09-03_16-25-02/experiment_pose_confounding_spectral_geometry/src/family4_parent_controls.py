"""Parent controls for Family 4 metric and geometry confounds.

Adds frequency-block SNR matching, per-pose SNR matching for trajectory
comparisons, same-radius curved paths, explicit sampled path geometry, and
conditioning diagnostics for the most-confounded map directions.  The raw
Family-4 results and its trajectory-ordering counterexample remain unchanged.

Run from the experiment root:
    .venv/bin/python src/family4_parent_controls.py
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

import family4_frequency_trajectory as f4  # noqa: E402
import helmholtz as hh  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def realify_vector(z: np.ndarray) -> np.ndarray:
    return np.sqrt(2.0) * np.concatenate([z.real, z.imag])


def build_block(
    chi: np.ndarray,
    poses: np.ndarray,
    S: np.ndarray,
    cfg: dict,
    f_scale: float = 1.0,
) -> dict:
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)
    A_c, B_c, F_c, _ = hh.build_AB(
        chi, poses, rx, tx, cfg["N"], 2 * np.pi * f_scale
    )
    A_pix, B = hh.whiten_realify(A_c, B_c, None)
    return {"A": A_pix @ S, "B": B, "F": realify_vector(F_c)}


def retention_stats(A: np.ndarray, B: np.ndarray) -> dict:
    dA = f4.family2.thin_decomposition(A)
    rA, Q = dA["rank"], dA["Q"]
    rB, Z, _, _ = f4.range_basis(B)
    cos = np.linalg.svd(Z.T @ Q, compute_uv=False)
    cos2 = cos**2
    rho = np.sort(
        np.concatenate([1.0 - cos2, np.ones(max(rA - len(cos2), 0))])
    )[::-1]
    return {
        "rank_A": rA,
        "rank_B": rB,
        "retained_mass": float(np.sum(rho)),
        "confusable_mass": float(np.sum(cos2)),
        "rho_min": float(rho[-1]),
        "theta_min_deg": float(np.degrees(np.arccos(min(1.0, cos[0])))),
        "log_volume": float(np.sum(np.log(np.maximum(rho, 1e-300)))),
        "count_rho_lt_1e-6": int(np.sum(rho < 1e-6)),
        "rho_desc": [float(x) for x in rho],
    }


def pose_block_indices(T: int, n_rx: int, t: int) -> np.ndarray:
    m_complex = T * n_rx
    re = np.arange(t * n_rx, (t + 1) * n_rx)
    im = m_complex + re
    return np.concatenate([re, im])


def snr_match_poses(block: dict, T: int, n_rx: int) -> tuple[dict, dict]:
    norms = []
    for t in range(T):
        idx = pose_block_indices(T, n_rx, t)
        norms.append(float(np.linalg.norm(block["F"][idx])))
    target = float(np.median(norms))
    scales = np.asarray([target / max(x, np.finfo(float).tiny) for x in norms])
    row_scale = np.ones_like(block["F"])
    for t, scale in enumerate(scales):
        row_scale[pose_block_indices(T, n_rx, t)] = scale
    return (
        {
            "A": row_scale[:, None] * block["A"],
            "B": row_scale[:, None] * block["B"],
            "F": row_scale * block["F"],
        },
        {
            "nominal_signal_norms": norms,
            "target_median_norm": target,
            "row_scales_by_pose": [float(x) for x in scales],
            "metric": (
                "Pose-block row scaling corresponding to noise standard "
                "deviation proportional to each nominal signal norm; this is "
                "an SNR-matched sensitivity metric, not raw identity noise."
            ),
        },
    )


def path_geometry(poses: np.ndarray, closed: bool) -> dict:
    p = poses[:, :2]
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    closure = float(np.linalg.norm(p[-1] - p[0]))
    ranges = np.linalg.norm(p, axis=1)
    return {
        "open_polyline_length": float(np.sum(seg)),
        "closure_segment_length": closure,
        "counted_path_length": float(np.sum(seg) + (closure if closed else 0.0)),
        "closed": closed,
        "min_target_range": float(np.min(ranges)),
        "max_target_range": float(np.max(ranges)),
    }


def same_radius_poses(name: str, T: int, radius: float = 1.6) -> np.ndarray:
    if name == "arc90":
        phi = np.linspace(-np.pi / 4, np.pi / 4, T)
    elif name == "arc180":
        phi = np.linspace(-np.pi / 2, np.pi / 2, T)
    elif name == "circle360":
        phi = np.linspace(0, 2 * np.pi, T, endpoint=False)
    else:
        raise ValueError(name)
    p = radius * np.column_stack([np.cos(phi), np.sin(phi)])
    theta = np.arctan2(-p[:, 1], -p[:, 0])
    return np.column_stack([p, theta])


def trajectory_controls(chi: np.ndarray, S: np.ndarray, cfg: dict) -> dict:
    T, n_rx = cfg["T"], cfg["n_rx"]
    equal_rows = []
    for name in ("straight", "arc90", "arc180", "circle360"):
        poses = f4.trajectory_poses(name, cfg)
        block = build_block(chi, poses, S, cfg)
        matched, match_meta = snr_match_poses(block, T, n_rx)
        equal_rows.append(
            {
                "trajectory": name,
                "path_geometry": path_geometry(poses, name == "circle360"),
                "identity_noise": retention_stats(block["A"], block["B"]),
                "snr_matched": retention_stats(matched["A"], matched["B"]),
                "snr_match_metadata": match_meta,
            }
        )

    radius_rows = []
    for name in ("arc90", "arc180", "circle360"):
        poses = same_radius_poses(name, T)
        block = build_block(chi, poses, S, cfg)
        matched, match_meta = snr_match_poses(block, T, n_rx)
        radius_rows.append(
            {
                "trajectory": name,
                "path_geometry": path_geometry(poses, name == "circle360"),
                "identity_noise": retention_stats(block["A"], block["B"]),
                "snr_matched": retention_stats(matched["A"], matched["B"]),
                "snr_match_metadata": match_meta,
            }
        )

    def orders(rows: list[dict], metric_name: str) -> dict:
        out = {}
        for metric in ("retained_mass", "rho_min", "log_volume"):
            out[metric] = sorted(
                [r["trajectory"] for r in rows],
                key=lambda name: next(
                    r[metric_name][metric] for r in rows if r["trajectory"] == name
                ),
                reverse=True,
            )
        return out

    return {
        "equal_continuous_length": {
            "rows": equal_rows,
            "orders_identity_noise": orders(equal_rows, "identity_noise"),
            "orders_snr_matched": orders(equal_rows, "snr_matched"),
        },
        "same_radius_curved_paths": {
            "radius": 1.6,
            "rows": radius_rows,
            "orders_identity_noise": orders(radius_rows, "identity_noise"),
            "orders_snr_matched": orders(radius_rows, "snr_matched"),
        },
        "interpretation": (
            "These controls show whether rankings persist after changing the "
            "metric or holding curved-path standoff fixed.  They still do not "
            "prove a universal trajectory ordering."
        ),
    }


def frequency_control(
    chi: np.ndarray, poses: np.ndarray, S: np.ndarray, cfg: dict
) -> dict:
    freqs = [1.0, 1.4, 1.8]
    raw = {f: build_block(chi, poses, S, cfg, f) for f in freqs}
    norms = {f: float(np.linalg.norm(raw[f]["F"])) for f in freqs}
    target = norms[1.0]
    scales = {f: target / max(norms[f], np.finfo(float).tiny) for f in freqs}
    blocks = {
        f: {
            "A": scales[f] * raw[f]["A"],
            "B": scales[f] * raw[f]["B"],
            "F": scales[f] * raw[f]["F"],
        }
        for f in freqs
    }

    A1, B1 = blocks[1.0]["A"], blocks[1.0]["B"]
    _, P1 = f4.range_projection(B1)
    rho, U, _ = f4.generalized_eigen_directions(A1, P1)
    A2 = np.vstack([blocks[1.0]["A"], blocks[1.4]["A"]])
    B2 = np.vstack([blocks[1.0]["B"], blocks[1.4]["B"]])
    _, P2 = f4.range_projection(B2)
    rows = []
    for j in range(3):
        u = U[:, j]
        rows.append(
            {
                "direction": j,
                "euclidean_norm_u": float(np.linalg.norm(u)),
                "max_abs_component_u": float(np.max(np.abs(u))),
                "rho_single": f4.rayleigh_quotient(A1, P1, u),
                "rho_two_frequency_snr_matched": f4.rayleigh_quotient(A2, P2, u),
                "shared_z_residual_two_frequency_snr_matched": f4.shared_z_residual(
                    A2, B2, u
                ),
            }
        )
    K1 = f4.K_SLAM(A1, B1)
    K2 = f4.K_SLAM(A2, B2)
    min_eig = float(np.linalg.eigvalsh(K2 - K1)[0])
    return {
        "nominal_signal_norm_by_frequency": {str(k): v for k, v in norms.items()},
        "frequency_row_scale": {str(k): v for k, v in scales.items()},
        "metric": "one scalar row scaling per frequency to match nominal ||F_f||_2",
        "min_eig_KSLAM_F2_minus_F1": min_eig,
        "rows": rows,
        "direction_conditioning_note": (
            "The K_IS-unit directions can have very large Euclidean coefficient "
            "norms when A is ill-conditioned.  Their local physical amplitude "
            "must therefore be scaled accordingly; the reported rho values are "
            "Fisher-direction diagnostics, not unit-contrast perturbations."
        ),
    }


def main() -> None:
    started = datetime.now(timezone.utc)
    t0 = time.perf_counter()
    cfg = f4.CONFIG
    _, chi, _, S = f4.base_scene(cfg)
    poses = f4.family1.build_poses(cfg)
    original_path = _ROOT / "results" / "family4_results.json"
    original = json.loads(original_path.read_text())
    frequency = frequency_control(chi, poses, S, cfg)
    trajectories = trajectory_controls(chi, S, cfg)

    result = {
        "schema": "family4_parent_controls",
        "generated_utc": started.isoformat(),
        "command": ".venv/bin/python src/family4_parent_controls.py",
        "working_directory": str(_ROOT),
        "raw_family4_hypothesis_holds": original["pass_summary"][
            "D_hypothesis_holds"
        ],
        "frequency_snr_matched": frequency,
        "trajectory_controls": trajectories,
        "runtime_seconds": float(time.perf_counter() - t0),
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "source_sha256": {
            "src/helmholtz.py": sha256(_HERE / "helmholtz.py"),
            "src/family4_frequency_trajectory.py": sha256(
                _HERE / "family4_frequency_trajectory.py"
            ),
            "src/family4_parent_controls.py": sha256(Path(__file__)),
            "results/family4_results.json": sha256(original_path),
        },
        "cannot_establish": (
            "SNR matching is a declared alternative noise metric, not evidence "
            "of the real sensor noise model.  Same-radius and equal-length "
            "finite scenarios do not prove universal or causal trajectory laws."
        ),
    }
    out = _ROOT / "results" / "family4_parent_controls.json"
    out.write_text(json.dumps(result, indent=2) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8))
    frows = frequency["rows"]
    x = np.arange(3)
    axes[0].bar(x - 0.18, [r["rho_single"] for r in frows], 0.36, label="single")
    axes[0].bar(
        x + 0.18,
        [r["rho_two_frequency_snr_matched"] for r in frows],
        0.36,
        label="1.0+1.4, SNR matched",
    )
    axes[0].set_xticks(x, [f"u{i}" for i in x])
    axes[0].set_ylabel("directional retention")
    axes[0].set_title("Frequency diversity after block SNR matching")
    axes[0].legend(fontsize=7)
    axes[0].grid(True, axis="y", alpha=0.25)

    eq = trajectories["equal_continuous_length"]["rows"]
    names = [r["trajectory"] for r in eq]
    xx = np.arange(len(names))
    axes[1].bar(
        xx - 0.18,
        [r["identity_noise"]["retained_mass"] for r in eq],
        0.36,
        label="identity noise",
    )
    axes[1].bar(
        xx + 0.18,
        [r["snr_matched"]["retained_mass"] for r in eq],
        0.36,
        label="per-pose SNR matched",
    )
    axes[1].set_xticks(xx, names, rotation=20)
    axes[1].set_ylabel("retained mass")
    axes[1].set_title("Equal design length: metric sensitivity")
    axes[1].legend(fontsize=7)
    axes[1].grid(True, axis="y", alpha=0.25)

    sr = trajectories["same_radius_curved_paths"]["rows"]
    sn = [r["trajectory"] for r in sr]
    sx = np.arange(len(sn))
    axes[2].bar(
        sx - 0.18,
        [r["identity_noise"]["retained_mass"] for r in sr],
        0.36,
        label="identity noise",
    )
    axes[2].bar(
        sx + 0.18,
        [r["snr_matched"]["retained_mass"] for r in sr],
        0.36,
        label="per-pose SNR matched",
    )
    axes[2].set_xticks(sx, sn, rotation=20)
    axes[2].set_ylabel("retained mass")
    axes[2].set_title("Same radius 1.6: curved-path control")
    axes[2].legend(fontsize=7)
    axes[2].grid(True, axis="y", alpha=0.25)
    fig.suptitle("Family 4 parent controls: metric, standoff, and path sampling")
    fig.tight_layout()
    fig_path = _ROOT / "figures" / "family4_parent_controls.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    note = _ROOT / "notes" / "family4_parent_controls.md"
    note.write_text(
        "# Family 4 parent controls\n\n"
        "The raw trajectory-ordering hypothesis remains false.  Frequency "
        "diversity was rerun after matching each frequency block's nominal "
        "signal norm; the three selected retentions become %s.  Trajectory "
        "rankings are reported under identity noise, per-pose SNR matching, "
        "equal continuous design length, and a same-radius curved-path control. "
        "Differences between these rankings demonstrate objective/metric and "
        "standoff sensitivity, not a universal ordering.\n"
        % [round(r["rho_two_frequency_snr_matched"], 6) for r in frows]
    )
    print("Family 4 parent controls complete")
    print("SNR-matched frequency directions:", [
        r["rho_two_frequency_snr_matched"] for r in frows
    ])
    print("equal-length orders:", trajectories["equal_continuous_length"][
        "orders_snr_matched"
    ])
    print("same-radius orders:", trajectories["same_radius_curved_paths"][
        "orders_snr_matched"
    ])
    print("results ->", out)
    print("figure  ->", fig_path)
    print("note    ->", note)


if __name__ == "__main__":
    main()
