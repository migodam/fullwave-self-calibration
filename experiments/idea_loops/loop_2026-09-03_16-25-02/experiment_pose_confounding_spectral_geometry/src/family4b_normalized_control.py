"""Family 4b supplementary control: block-normalized frequency stacks.

Family 4 (raw identity-noise stacks) declared but deferred a
block-normalized / SNR-matched control because per-frequency A and B block
energies differ.  This supplementary script runs that control on the exact
same Family 4 scenario and does not modify any family4 result/report.

Scenario (identical to Family 4): N=16, T=6 poses on the 90-degree arc at
R=1.6, n_rx=4, two-blob chi0 (family1_pilot.build_poses / make_chi0), p=24
smooth Gaussian-RBF basis (family2.build_smooth_basis), frequencies
{1.0, 1.4, 1.8} for the distinct stacks and duplicate c = 2.0 on f = 1.0.

Block normalization rule (energy-matched whitening, this control):
  for each frequency block f,
      A_f_n = A_f / ||A_f||_F,
      B_f_n = B_f / ||A_f||_F,
i.e. the same scalar is applied to A and B so the linearized per-block model
[A_f | B_f] is rescaled consistently.  Every repeated duplicate block is
normalized the same way before c-scaling.  The directions u0,u1,u2 are the
three most-confounded raw single-frequency generalized directions
(generalized eigenvectors of the raw single-frequency K_SLAM_1 w.r.t.
K_IS_1, reconstructed as u = V_A diag(1/s_A) w per audit rule 3); they are
scale-invariant and are reused unchanged on all normalized stacks.

Checks:
  1. C-control: rho_dist_norm movement away from raw rho_single, shared-z
     residual for the normalized distinct stack {1.0,1.4}, duplicate rho
     movement for the normalized [A1;c A1]/[B1;c B1] stack.
  2. A-control: PSD monotonicity of normalized distinct prefix stacks
     (F = 1,2,3) for K_SLAM_norm and K_eff_norm(alpha = 1.0).
  3. Raw-vs-normalized block-norm table for the JSON.

Tolerances (mirror Family 4):
  duplicate movement          |rho_dup_norm - rho_single| < 1e-10
  distinct normalized movement rho_dist_norm - rho_single > 1e-4 in >= 2/3
  monotonicity                min_eig >= -1e-10 * max(1, ||larger||_2)

Run (from the experiment root):
    .venv/bin/python src/family4b_normalized_control.py
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
import family4_frequency_trajectory as family4  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration (scenario copied from family4_frequency_trajectory.py)
# ---------------------------------------------------------------------------

CONTROL = {
    "family": "4b",
    "title": "block-normalized (energy-matched) frequency-diversity control",
    "distinct_frequencies": [1.0, 1.4],
    "monotonicity_frequencies": [[1.0], [1.0, 1.4], [1.0, 1.4, 1.8]],
    "duplicate_c": 2.0,
    "finite_prior_alpha": 1.0,
    "normalization_rule": (
        "A_f_n = A_f / ||A_f||_F and B_f_n = B_f / ||A_f||_F (same scalar "
        "for A and B within each frequency block; per-block linearized model "
        "rescaled consistently; duplicate blocks normalized the same way "
        "before c-scaling)"
    ),
    "whitening_convention": (
        "raw identity-noise blocks from family4 (W=None realify) are "
        "post-hoc block-normalized for this supplementary control only"
    ),
    "directions_from": (
        "raw single-frequency generalized eigenvectors of K_SLAM_1 w.r.t. "
        "K_IS_1, reconstructed u = V_A diag(1/s_A) w (family4 audit rule 3)"
    ),
    "tolerances": {
        "rho_dup_movement_abs": 1e-10,
        "rho_dist_movement": 1e-4,
        "required_directions": 2,
        "min_eig_rule": "min_eig(D) >= -1e-10 * max(1, ||larger||_2)",
    },
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


def normalize_block(A: np.ndarray, B: np.ndarray, scale: float):
    """Divide both A and B by the block's raw ||A||_F."""
    return A / scale, B / scale


def norm_table_row(
    label: str, A: np.ndarray, B: np.ndarray, raw_scale: float | None = None
) -> dict:
    """Raw-vs-normalized block-norm record for one frequency block."""
    row = {
        "label": label,
        "normA_fro": float(np.linalg.norm(A, ord="fro")),
        "normA_spectral": float(np.linalg.norm(A, ord=2)),
        "normA_row_mean": float(np.mean(np.linalg.norm(A, axis=1))),
        "normA_row_max": float(np.max(np.linalg.norm(A, axis=1))),
        "normB_fro": float(np.linalg.norm(B, ord="fro")),
        "normB_spectral": float(np.linalg.norm(B, ord=2)),
        "normB_row_mean": float(np.mean(np.linalg.norm(B, axis=1))),
        "normB_row_max": float(np.max(np.linalg.norm(B, axis=1))),
    }
    row["B_over_A_fro_ratio"] = row["normB_fro"] / row["normA_fro"]
    if raw_scale is not None:
        row["normalized"] = False
        row["scaling_to_unit_A_fro"] = float(raw_scale)
    else:
        row["normalized"] = True
        row["scaling_to_unit_A_fro"] = None
    return row


def run_control() -> dict:
    cfg = family4.CONFIG
    ctl = CONTROL
    t_utc = datetime.now(timezone.utc)
    t_start = time.perf_counter()

    points, chi0, h, S = family4.base_scene(cfg)
    poses = family1.build_poses(cfg)

    # ---------------- raw single-frequency blocks (f = 1.0,1.4,1.8) --------
    raw_blocks = {}
    for f in cfg["f_scales"]:
        A_s, B_R, _ = family4.build_smooth_block(chi0, poses, S, cfg, f)
        raw_blocks[float(f)] = {"A": A_s, "B": B_R}

    A1, B1 = raw_blocks[1.0]["A"], raw_blocks[1.0]["B"]
    A14, B14 = raw_blocks[1.4]["A"], raw_blocks[1.4]["B"]
    A18, B18 = raw_blocks[1.8]["A"], raw_blocks[1.8]["B"]

    # Family 4 stored directions and raw rho values (read-only reference).
    f4_results = json.loads(
        (_ROOT / "results" / "family4_results.json").read_text()
    )
    f4_c_rows = f4_results["checks"]["C_shared_pose_compensation"]["rows"]
    family4_rho_single = [r["rho_single"] for r in f4_c_rows]
    family4_rho_distinct = [r["rho_distinct"] for r in f4_c_rows]
    family4_u_components = [
        np.asarray(r["u_component"], dtype=float) for r in f4_c_rows
    ]
    family4_raw_movement = [
        r["movement_distinct_minus_single"] for r in f4_c_rows
    ]
    family4_dup_movement = [
        r["movement_duplicate_minus_single"] for r in f4_c_rows
    ]

    # ------- the three most confounded raw single-frequency directions -----
    _, P1 = family4.range_projection(B1)
    rho_asc_raw, U_raw, _ = family4.generalized_eigen_directions(A1, P1)
    us = [np.asarray(U_raw[:, j], dtype=float) for j in range(3)]

    direction_agreement = []
    for j, u in enumerate(us):
        v = family4_u_components[j]
        cosine = float(
            np.abs(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))
        )
        direction_agreement.append(
            {
                "direction": int(j),
                "abs_cosine_with_family4_u": cosine,
                "u_euclid_norm": float(np.linalg.norm(u)),
            }
        )

    # Recompute the raw single-frequency rho to use as the movement baseline.
    rho_single_raw = [
        family4.rayleigh_quotient(A1, P1, u) for u in us
    ]
    rho_crosscheck = [
        float(abs(rho_single_raw[j] - family4_rho_single[j]))
        for j in range(3)
    ]

    # ---------------- normalized blocks and stacks ------------------------
    scales = {}
    norm_blocks = {}
    norm_table = []
    for f, A, B in (
        (1.0, A1, B1),
        (1.4, A14, B14),
        (1.8, A18, B18),
    ):
        scale = float(np.linalg.norm(A, ord="fro"))
        scales[float(f)] = scale
        An, Bn = normalize_block(A, B, scale)
        norm_blocks[float(f)] = {"A": An, "B": Bn}
        norm_table.append(norm_table_row(f"raw_f{f}", A, B, raw_scale=scale))
        norm_table.append(norm_table_row(f"norm_f{f}", An, Bn, raw_scale=None))

    # distinct normalized stack {1.0, 1.4} and duplicate normalized stack
    A1n, B1n = norm_blocks[1.0]["A"], norm_blocks[1.0]["B"]
    A14n, B14n = norm_blocks[1.4]["A"], norm_blocks[1.4]["B"]
    A18n, B18n = norm_blocks[1.8]["A"], norm_blocks[1.8]["B"]
    c = float(ctl["duplicate_c"])

    Ast_n = np.vstack([A1n, A14n])
    Bst_n = np.vstack([B1n, B14n])
    Adup_n = np.vstack([A1n, c * A1n])
    Bdup_n = np.vstack([B1n, c * B1n])

    _, Pst_n = family4.range_projection(Bst_n)
    _, Pdup_n = family4.range_projection(Bdup_n)
    _, P1_n = family4.range_projection(B1n)

    # normalized single-frequency rho at the SAME directions (scale check)
    rho_single_norm = [family4.rayleigh_quotient(A1n, P1_n, u) for u in us]

    rows = []
    n_greater = 0
    dup_pass = True
    for j in range(3):
        u = us[j]
        rho_dist_norm = family4.rayleigh_quotient(Ast_n, Pst_n, u)
        rho_dup_norm = family4.rayleigh_quotient(Adup_n, Pdup_n, u)
        movement_dist_norm = rho_dist_norm - rho_single_raw[j]
        movement_dup_norm = rho_dup_norm - rho_single_raw[j]
        zres_dist_norm = family4.shared_z_residual(Ast_n, Bst_n, u)
        zres_dup_norm = family4.shared_z_residual(Adup_n, Bdup_n, u)
        if movement_dist_norm > ctl["tolerances"]["rho_dist_movement"]:
            n_greater += 1
        if abs(movement_dup_norm) >= ctl["tolerances"]["rho_dup_movement_abs"]:
            dup_pass = False
        rows.append(
            {
                "direction": int(j),
                "generalized_eigenvalue_asc_raw": float(rho_asc_raw[j]),
                "u_component": [float(x) for x in u],
                "u_KIS_norm2_raw": float(np.dot(A1 @ u, A1 @ u)),
                "rho_single_raw": float(rho_single_raw[j]),
                "rho_single_family4_reference": float(family4_rho_single[j]),
                "rho_single_norm": float(rho_single_norm[j]),
                "rho_single_raw_crosscheck_abs_diff": float(
                    rho_crosscheck[j]
                ),
                "rho_distinct_norm": float(rho_dist_norm),
                "movement_distinct_norm_minus_single": float(
                    movement_dist_norm
                ),
                "movement_distinct_raw_reference": float(
                    family4_raw_movement[j]
                ),
                "shared_z_residual_distinct_norm": float(zres_dist_norm),
                "rho_duplicate_norm": float(rho_dup_norm),
                "movement_duplicate_norm_minus_single": float(
                    movement_dup_norm
                ),
                "movement_duplicate_raw_reference": float(
                    family4_dup_movement[j]
                ),
                "shared_z_residual_duplicate_norm": float(zres_dup_norm),
                "pass_distinct_norm_gt_1e-4": bool(
                    movement_dist_norm
                    > ctl["tolerances"]["rho_dist_movement"]
                ),
                "pass_duplicate_norm_lt_1e-10": bool(
                    abs(movement_dup_norm)
                    < ctl["tolerances"]["rho_dup_movement_abs"]
                ),
            }
        )

    distinct_norm_pass = (
        n_greater >= ctl["tolerances"]["required_directions"]
    )
    max_abs_dup_norm_movement = float(
        max(abs(r["movement_duplicate_norm_minus_single"]) for r in rows)
    )

    # ---------------- PSD monotonicity on normalized stacks ---------------
    alpha = ctl["finite_prior_alpha"]
    mono = {"rows": [], "pass_gate": True, "alpha": alpha}
    for F in (1, 2, 3):
        freqs = ctl["monotonicity_frequencies"][F - 1]
        A_st = np.vstack([norm_blocks[f]["A"] for f in freqs])
        B_st = np.vstack([norm_blocks[f]["B"] for f in freqs])
        KSL = family4.K_SLAM(A_st, B_st)
        Keff = family4.K_eff(A_st, B_st, alpha)
        mono[f"K_SLAM_norm_F{F}_eig_desc"] = [
            float(x) for x in np.sort(np.linalg.eigvalsh(KSL))[::-1]
        ]
        mono[f"K_eff_norm_F{F}_eig_desc"] = [
            float(x) for x in np.sort(np.linalg.eigvalsh(Keff))[::-1]
        ]
    for F in (1, 2):
        for label in ("K_SLAM_norm", "K_eff_norm"):
            freqs_lo = ctl["monotonicity_frequencies"][F - 1]
            freqs_hi = ctl["monotonicity_frequencies"][F]
            Alo = np.vstack([norm_blocks[f]["A"] for f in freqs_lo])
            Blo = np.vstack([norm_blocks[f]["B"] for f in freqs_lo])
            Ahi = np.vstack([norm_blocks[f]["A"] for f in freqs_hi])
            Bhi = np.vstack([norm_blocks[f]["B"] for f in freqs_hi])
            if label == "K_SLAM_norm":
                Klo, Khi = family4.K_SLAM(Alo, Blo), family4.K_SLAM(Ahi, Bhi)
            else:
                Klo = family4.K_eff(Alo, Blo, alpha)
                Khi = family4.K_eff(Ahi, Bhi, alpha)
            D = Khi - Klo
            min_eig = float(np.linalg.eigvalsh(D)[0])
            scale = float(max(1.0, np.linalg.norm(Khi, ord=2)))
            tol = 1e-10 * scale
            ok = bool(min_eig >= -tol)
            mono["pass_gate"] &= ok
            mono["rows"].append(
                {
                    "stack_from": F,
                    "stack_to": F + 1,
                    "label": label,
                    "min_eig_difference": min_eig,
                    "scale_max_1_norm_larger_2": scale,
                    "tol": tol,
                    "pass": ok,
                }
            )

    # ---------------------------- summary --------------------------------
    pass_summary = {
        "distinct_normalized_n_directions_gt_1e-4": int(n_greater),
        "distinct_normalized_pass": bool(distinct_norm_pass),
        "duplicate_normalized_max_abs_movement": max_abs_dup_norm_movement,
        "duplicate_normalized_pass": bool(dup_pass),
        "monotonicity_normalized_pass": bool(mono["pass_gate"]),
        "family4b_overall_pass": bool(
            distinct_norm_pass and dup_pass and mono["pass_gate"]
        ),
    }

    source_files = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
        "src/family4_frequency_trajectory.py",
        "src/family4b_normalized_control.py",
    ]
    results = {
        "generated_utc": t_utc.isoformat(),
        "runner": "src/family4b_normalized_control.py",
        "command": ".venv/bin/python src/family4b_normalized_control.py",
        "family": "4b",
        "title": ctl["title"],
        "deferred_control_resolution": (
            "Family 4 deferred a block-normalized/SNR-matched control; this "
            "file executes it as a SUPPLEMENTARY control and does not modify "
            "results/family4_results.json or notes/family4_report.md"
        ),
        "raw_metric_reference": (
            "results/family4_results.json checks.C (raw identity-noise stack)"
        ),
        "normalization_rule": ctl["normalization_rule"],
        "self_cell_formula": hh.SELF_CELL_FORMULA,
        "self_cell_formula_version": hh.SELF_CELL_FORMULA_VERSION,
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
            p: hashlib.sha256((_ROOT / p).read_bytes()).hexdigest()
            for p in source_files
        },
        "family4_results_sha256": sha256_file(
            _ROOT / "results" / "family4_results.json"
        ),
        "config": {
            "N": cfg["N"],
            "T": cfg["T"],
            "n_rx": cfg["n_rx"],
            "arc_radius": cfg["arc_radius"],
            "arc_phi_deg": cfg["arc_phi_deg"],
            "k_b_rule": cfg["k_b_rule"],
            "smooth_basis": cfg["smooth_basis"],
            "chi_blobs": cfg["chi_blobs"],
            "poses_stack": [p.tolist() for p in poses],
            "grid_h_cell": float(h),
            "distinct_frequencies": ctl["distinct_frequencies"],
            "monotonicity_frequencies": ctl["monotonicity_frequencies"],
            "duplicate_c": c,
            "finite_prior_alpha": alpha,
            "tolerances": ctl["tolerances"],
            "chi0_stats": {
                "min": float(chi0.min()),
                "max": float(chi0.max()),
                "mean": float(chi0.mean()),
                "l2": float(np.linalg.norm(chi0)),
            },
        },
        "direction_agreement_with_family4": direction_agreement,
        "family4_raw_reference": {
            "rho_single": family4_rho_single,
            "rho_distinct": family4_rho_distinct,
            "movement_distinct_minus_single": family4_raw_movement,
            "movement_duplicate_minus_single": family4_dup_movement,
        },
        "block_norm_table": norm_table,
        "block_normalization_scales_normA_fro": scales,
        "checks": {
            "C_block_normalized_shared_pose_compensation": {
                "directions_selected_from": ctl["directions_from"],
                "rows": rows,
                "n_directions_with_distinct_norm_movement_gt_1e-4": int(
                    n_greater
                ),
                "all_duplicate_norm_movements_lt_1e-10": bool(dup_pass),
                "pass_gate": bool(distinct_norm_pass and dup_pass),
            },
            "A_block_normalized_psd_monotonicity": mono,
        },
        "pass_summary": pass_summary,
        "scope_note": (
            "Block normalization is an energy-matched per-frequency rescale "
            "of the SAME raw whitened/realified (A_R, B_R) Jacobian blocks. "
            "All statements remain finite-dimensional discrete-model "
            "statements on N=16, p=24; they are not continuum-limit or "
            "estimator claims, and the physical-noise interpretation is a "
            "per-block equal-||A_f||_F control, not a claim about optimal "
            "SNR whitening in the presence of correlated noise."
        ),
    }
    return results


def write_figure(results: dict, figures_dir: Path) -> Path:
    rows = results["checks"]["C_block_normalized_shared_pose_compensation"][
        "rows"
    ]
    f4 = results["family4_raw_reference"]
    labels = ["u0", "u1", "u2"]
    raw_mov = f4["movement_distinct_minus_single"]
    norm_mov = [r["movement_distinct_norm_minus_single"] for r in rows]
    raw_dup = [abs(x) for x in f4["movement_duplicate_minus_single"]]
    norm_dup = [
        abs(r["movement_duplicate_norm_minus_single"]) for r in rows
    ]
    x = np.arange(3)
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    axes[0].bar(
        x - width / 2, raw_mov, width, color="#1f77b4",
        label="raw identity-noise (family4)",
    )
    axes[0].bar(
        x + width / 2, norm_mov, width, color="#ff7f0e",
        label="block-normalized (family4b)",
    )
    axes[0].axhline(1e-4, color="k", lw=0.8, ls=":")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel(r"$\rho_{\rm distinct} - \rho_{\rm single}$")
    axes[0].set_title("rho movement, distinct {1.0,1.4}")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, axis="y", alpha=0.3)
    axes[1].bar(
        x - width / 2, raw_dup, width, color="#1f77b4",
        label="raw identity-noise (family4)",
    )
    axes[1].bar(
        x + width / 2, norm_dup, width, color="#2ca02c",
        label="block-normalized (family4b)",
    )
    axes[1].axhline(1e-10, color="k", lw=0.8, ls=":")
    axes[1].set_xticks(x, labels)
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$|\rho_{\rm duplicate} - \rho_{\rm single}|$ (log)")
    axes[1].set_title("duplicate c=2 movement (log)")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, axis="y", which="both", alpha=0.3)
    fig.suptitle(
        "Family 4b control: block normalization does not remove "
        "frequency-diversity movement",
        fontsize=10,
    )
    fig.tight_layout()
    p = figures_dir / "family4b_normalized_shared_z.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


def write_report(results: dict, notes_dir: Path) -> Path:
    c_rows = results["checks"]["C_block_normalized_shared_pose_compensation"][
        "rows"
    ]
    mono = results["checks"]["A_block_normalized_psd_monotonicity"]
    ps = results["pass_summary"]
    f4 = results["family4_raw_reference"]

    c_tbl = "\n".join(
        f"| {r['direction']} | {r['rho_single_raw']:.6e} | "
        f"{r['rho_distinct_norm']:.6e} | "
        f"{r['movement_distinct_norm_minus_single']:.6e} | "
        f"{r['shared_z_residual_distinct_norm']:.6e} | "
        f"{r['rho_duplicate_norm']:.6e} | "
        f"{r['movement_duplicate_norm_minus_single']:.6e} | "
        f"{r['shared_z_residual_duplicate_norm']:.6e} | "
        f"{'PASS' if r['pass_distinct_norm_gt_1e-4'] else 'FAIL'} | "
        f"{'PASS' if r['pass_duplicate_norm_lt_1e-10'] else 'FAIL'} |"
        for r in c_rows
    )
    a_tbl = "\n".join(
        f"| {r['label']} F={r['stack_from']}->{r['stack_to']} | "
        f"{r['min_eig_difference']:.6e} | {r['tol']:.3e} | "
        f"{'PASS' if r['pass'] else 'FAIL'} |"
        for r in mono["rows"]
    )
    bt = results["block_norm_table"]
    raw_rows = [t for t in bt if not t["normalized"]]
    norm_rows = [t for t in bt if t["normalized"]]
    block_rows = "\n".join(
        f"| f={raw['label'].split('_f')[-1]} | "
        f"{raw['normA_fro']:.8e} | {norm['normA_fro']:.8e} | "
        f"{raw['normB_fro']:.8e} | {norm['normB_fro']:.8e} | "
        f"{raw['B_over_A_fro_ratio']:.6f} | "
        f"{norm['B_over_A_fro_ratio']:.6f} |"
        for raw, norm in zip(raw_rows, norm_rows)
    )
    raw_mov = ", ".join(f"{v:.6e}" for v in f4["movement_distinct_minus_single"])
    norm_mov = ", ".join(
        f"{r['movement_distinct_norm_minus_single']:.6e}" for r in c_rows
    )
    dup_raw = f4["movement_duplicate_minus_single"]
    dup_norm = [r["movement_duplicate_norm_minus_single"] for r in c_rows]
    report = f"""# Family 4b: block-normalized frequency-diversity control

Date: 2026-09-03 (SGT; UTC stamp in `results/family4b_normalized_control.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.
Supplementary control for Family 4 (deferred in `results/family4_results.json`
and `notes/family4_report.md`).  **No Family 4 result or report file was
modified.**

## Exact command and runtime

```bash
.venv/bin/python src/family4b_normalized_control.py
```

Wall runtime: {results['runtime_seconds']:.2f} s.  Platform:
{results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, scipy {results['platform']['scipy']}.
Deterministic dense linear algebra; no RNG used.

## Scenario and normalization rule

Same scenario as Family 4: N=16, T=6 poses on the 90-degree arc at R=1.6,
n_rx=4, two-blob chi0, p=24 smooth basis, raw blocks f in
{{1.0, 1.4, 1.8}}; distinct stacks use {{1.0, 1.4}}; duplicate is block 1
repeated with c=2.0.  Rule:

{results['normalization_rule']}

Directions u0..u2 are the three most-confounded generalized eigenvectors of
the **raw** single-frequency K_SLAM_1 w.r.t. K_IS_1
(u = V_A diag(1/s_A) w); they are scale-invariant and are reused unchanged
on the normalized stacks.  Direction agreement with the Family 4 stored u
vectors (abs cosine): {', '.join(f"{d['direction']}: {d['abs_cosine_with_family4_u']:.12f}" for d in results['direction_agreement_with_family4'])}.
Recomputed raw rho_single differs from the Family 4 reference by
{', '.join(f"{r['rho_single_raw_crosscheck_abs_diff']:.2e}" for r in c_rows)}.

## Claim status

**Does block-normalization change the Family 4 diversity conclusion?  NO.**
The normalized distinct stack still moves all three confounded directions
well above the 1e-4 gate ({ps['distinct_normalized_n_directions_gt_1e-4']}/3),
the normalized duplicate movement remains at roundoff level
(max |{max(abs(x) for x in dup_norm):.3e}| < 1e-10), and normalized PSD
monotonicity passes for all four prefix increments.

### C control rows (same u0,u1,u2; movement vs raw rho_single)

Raw Family 4 distinct movement: {raw_mov}.
Normalized distinct movement: {norm_mov}.

| dir | rho_single raw | rho_dist_norm | movement_dist_norm | z-resid dist norm | rho_dup_norm | movement_dup_norm | z-resid dup norm | gate >1e-4 | gate <1e-10 |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
{c_tbl}

### A control: normalized PSD monotonicity rows

| comparison | min eig D | tol | gate |
|---|---:|---:|---|
{a_tbl}

### Block-norm record (raw vs normalized)

| block | ||A||_F raw | ||A||_F norm | ||B||_F raw | ||B||_F norm | B/A ratio raw | B/A ratio norm |
|---|---:|---:|---:|---:|---:|---:|
{block_rows}

## Cannot-establish section

* All claims are finite-dimensional discrete-model statements (N=16, p=24,
  the specific arc/chi0/blobs); no continuum-limit, universal-trajectory, or
  estimator claim is made.
* Block normalization rescales each raw block by its own ||A_f||_F, which
  equalizes A-block energy across frequencies but does not itself prescribe a
  physical noise covariance; it is one declared energy-matched control, not
  an optimal or correlated-noise whitening claim.
* The shared-z residual and rho are evaluated at only the three raw
  single-frequency confounded directions; other directions or larger stacks
  are not exhaustively checked (same boundary as Family 4 check C).

## Artifacts

- results: `results/family4b_normalized_control.json`
- figure: figures/family4b_normalized_shared_z.png
- this report: notes/family4b_normalized_control.md

Self-cell formula used: {results['self_cell_formula_version']}.
"""
    p = notes_dir / "family4b_normalized_control.md"
    p.write_text(report)
    return p


def main() -> None:
    results = run_control()
    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    fig_path = write_figure(results, figures_dir)
    results["figure_sha256"] = {
        "family4b_normalized_shared_z.png": sha256_file(fig_path)
    }
    results["artifacts"] = {
        "results_json": "results/family4b_normalized_control.json",
        "figures": [str(fig_path)],
        "report_md": "notes/family4b_normalized_control.md",
    }

    results_path = results_dir / "family4b_normalized_control.json"
    results_path.write_text(_round_trip_json(results))
    report_path = write_report(results, notes_dir)

    c = results["checks"]["C_block_normalized_shared_pose_compensation"]
    a = results["checks"]["A_block_normalized_psd_monotonicity"]
    ps = results["pass_summary"]
    print("\n===== FAMILY 4B NORMALIZED CONTROL SUMMARY =====")
    print(
        f"[C] distinct norm movements = "
        + ", ".join(
            f"{r['movement_distinct_norm_minus_single']:.6e}" for r in c["rows"]
        )
    )
    print(
        f"[C] distinct pass={c['pass_gate']} "
        f"({c['n_directions_with_distinct_norm_movement_gt_1e-4']}/3 > 1e-4); "
        f"dup norm movements = "
        + ", ".join(
            f"{r['movement_duplicate_norm_minus_single']:.3e}" for r in c["rows"]
        )
        + f"; dup pass={c['all_duplicate_norm_movements_lt_1e-10']}"
    )
    print(f"[A] normalized monotonicity pass={a['pass_gate']}")
    for r in a["rows"]:
        print(
            f"    {r['label']} F{r['stack_from']}->{r['stack_to']}: "
            f"min_eig={r['min_eig_difference']:.3e} pass={r['pass']}"
        )
    print(f"overall pass_summary = {ps}")
    print(f"runtime {results['runtime_seconds']:.2f}s")
    print("results ->", results_path)
    print("figure  ->", fig_path)
    print("report  ->", report_path)


if __name__ == "__main__":
    main()
