"""Scientifically corrected Family-3 follow-up.

The first automated Family-3 run is deliberately preserved.  This follow-up
tests a compactly supported scene for discrete gauge refinement, interprets a
radial scene's rotational stabilizer by absolute residuals, implements a real
pose anchor, and isolates the Born mixed term with a joint centered experiment.

Run from the experiment root:
    .venv/bin/python src/family3_parent_correction.py
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

import family3_gauge_born as f3  # noqa: E402
import helmholtz as hh  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


GRID_NS = [16, 24, 32, 40]
BORN_STEPS = np.logspace(-4, -1, 9)
ANCHOR_ALPHAS = [0.0, 1e-4, 1e-2, 1.0, 1e2, 1e6]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def realify_vector(z: np.ndarray) -> np.ndarray:
    return np.sqrt(2.0) * np.concatenate([z.real, z.imag])


def compact_two_bump(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """C-infinity interior contrast and its analytic gradient.

    Each bump is exp(1-1/(1-q)) for q=|r-c|^2/R^2<1 and zero otherwise, so
    both the contrast and all derivatives vanish at the support boundary.
    """
    chi = np.zeros(points.shape[0])
    grad = np.zeros((points.shape[0], 2))
    specs = [
        (0.30, np.array([-0.15, -0.12]), 0.16),
        (0.50, np.array([0.18, 0.14]), 0.14),
    ]
    for amp, center, radius in specs:
        dr = points - center
        q = np.sum(dr**2, axis=1) / radius**2
        mask = q < 1.0
        val = amp * np.exp(1.0 - 1.0 / (1.0 - q[mask]))
        chi[mask] += val
        grad[mask] += val[:, None] * (
            -2.0
            * dr[mask]
            / (radius**2 * (1.0 - q[mask])[:, None] ** 2)
        )
    return chi, grad


def pixel_gauge_refinement(cfg: dict) -> dict:
    rows = []
    for N in GRID_NS:
        cfg_n = dict(cfg)
        cfg_n["N"] = N
        points, _ = hh.make_grid(N)
        poses = f3.build_poses(cfg_n)
        dXs = f3.gauge_dX(poses)

        scene_rows = {}
        for label in ("truncated_gaussian", "compact_C_infinity"):
            if label == "truncated_gaussian":
                chi = f3.make_chi0(points, cfg_n)
                grad = f3.grad_chi0(points, cfg_n)
            else:
                chi, grad = compact_two_bump(points)
            setup = f3.analysis_setup(chi, cfg_n)
            dchis = f3.gauge_dchi(grad, points)
            report = f3.gauge_rows(
                setup["A_pix_R"],
                setup["B_R"],
                None,
                dchis,
                dXs,
                include_kernel_distance=False,
            )
            scene_rows[label] = report["rows"]
        rows.append({"N": N, "scenes": scene_rows})
    return {
        "rows": rows,
        "interpretation": (
            "The noncompact Gaussians are truncated by the fixed square and "
            "can retain a boundary-action residual.  The compact C-infinity "
            "scene removes that confound and shows decreasing discrete gauge "
            "residuals.  Four grids are a refinement diagnostic, not a proof."
        ),
    }


def fit_log_slope(xs: list[float], ys: list[float], first_n: int = 5) -> dict:
    x = np.log10(np.asarray(xs[:first_n]))
    y = np.log10(np.maximum(np.asarray(ys[:first_n]), np.finfo(float).tiny))
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "slope": float(slope),
        "r_squared": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0,
        "indices": list(range(first_n)),
    }


def born_joint_test(cfg: dict) -> dict:
    N = cfg["N"]
    points, _ = hh.make_grid(N)
    poses = f3.build_poses(cfg)
    X0 = poses.reshape(-1)
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)
    k_b = cfg["k_b"]

    center = np.array([0.1, 0.05])
    dchi = np.exp(-np.sum((points - center) ** 2, axis=1) / (2 * 0.1**2))
    dchi /= np.linalg.norm(dchi)
    rng = np.random.default_rng(31415)
    dX = rng.standard_normal(X0.size)
    dX /= np.linalg.norm(dX)

    def born_A_real(X: np.ndarray) -> np.ndarray:
        _, A = hh.born_forward(
            np.zeros(N * N), X.reshape(-1, 3), rx, tx, N, k_b
        )
        return np.sqrt(2.0) * np.vstack([A.real, A.imag])

    def born_F_real(chi: np.ndarray, X: np.ndarray) -> np.ndarray:
        F, _ = hh.born_forward(chi, X.reshape(-1, 3), rx, tx, N, k_b)
        return realify_vector(F)

    def full_F_real(chi: np.ndarray, X: np.ndarray) -> np.ndarray:
        F = hh.forward_measurements(chi, X.reshape(-1, 3), rx, tx, N, k_b)
        return realify_vector(F)

    A0 = born_A_real(X0)
    delta_ref = 1e-4
    D_ref = (
        -born_A_real(X0 + 2 * delta_ref * dX)
        + 8 * born_A_real(X0 + delta_ref * dX)
        - 8 * born_A_real(X0 - delta_ref * dX)
        + born_A_real(X0 - 2 * delta_ref * dX)
    ) @ dchi / (12 * delta_ref)
    D_norm = float(np.linalg.norm(D_ref))
    linear_norm = float(np.linalg.norm(A0 @ dchi))

    rows = []
    for eps in BORN_STEPS:
        A_plus = born_A_real(X0 + eps * dX)
        A_minus = born_A_real(X0 - eps * dX)
        D_center = ((A_plus - A_minus) / (2 * eps)) @ dchi

        # Because F_B is linear in chi, J_B is exactly a joint-order eps^2
        # central term.  J_FW uses the same joint scaling and approaches it in
        # the weak-contrast limit while retaining multiple scattering.
        J_born = (
            born_F_real(eps * dchi, X0 + eps * dX)
            - born_F_real(eps * dchi, X0 - eps * dX)
        ) / 2.0
        J_full = (
            full_F_real(eps * dchi, X0 + eps * dX)
            - full_F_real(eps * dchi, X0 - eps * dX)
        ) / 2.0
        coeff_born = J_born / eps**2
        coeff_full = J_full / eps**2
        rows.append(
            {
                "epsilon": float(eps),
                "central_D_error_rel": float(
                    np.linalg.norm(D_center - D_ref) / D_norm
                ),
                "joint_Born_norm": float(np.linalg.norm(J_born)),
                "joint_Born_over_eps2_Dnorm": float(
                    np.linalg.norm(J_born) / (eps**2 * D_norm)
                ),
                "joint_Born_coefficient_error_rel": float(
                    np.linalg.norm(coeff_born - D_ref) / D_norm
                ),
                "Born_linearity_identity_error_rel": float(
                    np.linalg.norm(coeff_born - D_center) / D_norm
                ),
                "joint_fullwave_coefficient_error_rel": float(
                    np.linalg.norm(coeff_full - D_ref) / D_norm
                ),
            }
        )

    epss = [r["epsilon"] for r in rows]
    return {
        "definition": (
            "J_B(eps)=[F_B(eps*dchi,X+eps*dX)-F_B(eps*dchi,X-eps*dX)]/2; "
            "J_B/eps^2 -> (D_X A_B[dX])dchi"
        ),
        "dchi_l2": float(np.linalg.norm(dchi)),
        "dX_l2": float(np.linalg.norm(dX)),
        "dX_seed": 31415,
        "reference_derivative": "five-point centered stencil, delta=1e-4",
        "reference_D_norm": D_norm,
        "linear_A0_dchi_norm": linear_norm,
        "rows": rows,
        "slopes_first_five": {
            "central_D_error": fit_log_slope(
                epss, [r["central_D_error_rel"] for r in rows]
            ),
            "joint_Born_norm": fit_log_slope(
                epss, [r["joint_Born_norm"] for r in rows]
            ),
            "fullwave_coefficient_error": fit_log_slope(
                epss, [r["joint_fullwave_coefficient_error_rel"] for r in rows]
            ),
        },
        "interpretation": (
            "The Born joint signal scales as eps^2 and its coefficient has "
            "centered O(eps^2) error.  The full-wave coefficient approaches "
            "the Born coefficient roughly linearly in contrast amplitude "
            "because multiple scattering is retained."
        ),
    }


def pose_anchor_test(cfg: dict) -> dict:
    points, _ = hh.make_grid(cfg["N"])
    chi = f3.make_chi0(points, cfg)
    setup = f3.analysis_setup(chi, cfg)
    A, B = setup["A_pix_R"], setup["B_R"]
    dchis = f3.gauge_dchi(f3.grad_chi0(points, cfg), points)
    KIS = A.T @ A
    rB, Z, _, _ = f3.range_basis(B)
    K0 = A.T @ ((np.eye(A.shape[0]) - Z[:, :rB] @ Z[:, :rB].T) @ A)
    rows = []
    for alpha in ANCHOR_ALPHAS:
        if alpha == 0:
            K = K0
        else:
            J = np.zeros((B.shape[1], B.shape[1]))
            J[:3, :3] = alpha * np.eye(3)
            K = KIS - A.T @ B @ np.linalg.solve(B.T @ B + J, B.T @ A)
        directional = {}
        for key, vec in dchis.items():
            denom = float(vec @ KIS @ vec)
            directional[key] = float(vec @ K @ vec / denom)
        rows.append({"alpha": alpha, "directional_retention": directional})
    return {
        "prior": "J_X anchors the first pose's x,y,theta coordinates",
        "rows": rows,
        "interpretation": (
            "A physical pose prior lifts the tested global map-gauge directions "
            "from near-zero retention, but the values saturate below one because "
            "the remaining poses are still nuisance variables."
        ),
    }


def main() -> None:
    started = datetime.now(timezone.utc)
    t0 = time.perf_counter()
    cfg = f3.CONFIG
    original_path = _ROOT / "results" / "family3_gauge_born.json"
    original = json.loads(original_path.read_text())

    gauge = pixel_gauge_refinement(cfg)
    born = born_joint_test(cfg)
    anchor = pose_anchor_test(cfg)
    sym_rows = original["gauge_symmetric_smooth"]["rows"]
    sym_rot = next(r for r in sym_rows if r["generator"] == "rot")
    sym_trans_scale = max(
        r["pose_norm_B_dX"] for r in sym_rows if r["generator"] in ("tx", "ty")
    )
    symmetric_stabilizer = {
        "dchi_rotation_l2": sym_rot["dchi_l2"],
        "B_dX_rotation_l2": sym_rot["pose_norm_B_dX"],
        "B_dX_rotation_over_translation_scale": float(
            sym_rot["pose_norm_B_dX"] / sym_trans_scale
        ),
        "interpretation": (
            "The radial map tangent is zero.  The small absolute B*dX value is "
            "a discrete symmetry-breaking residual for a pose-only rotational "
            "stabilizer; the original normalized ratio near one is ill-scaled "
            "because both sides should be zero."
        ),
    }

    result = {
        "schema": "family3_parent_correction",
        "generated_utc": started.isoformat(),
        "command": ".venv/bin/python src/family3_parent_correction.py",
        "working_directory": str(_ROOT),
        "original_failures_preserved": original["gates"],
        "pixel_gauge_refinement": gauge,
        "symmetric_rotation_stabilizer": symmetric_stabilizer,
        "born_joint_second_order": born,
        "pose_anchor": anchor,
        "mask_experiment_interpretation": {
            "original_rows": original["anchor_variants"]["rows"],
            "conclusion": (
                "The automatic tenfold-improvement expectation failed.  These "
                "support masks are not universal physical anchors and two made "
                "rho_min smaller; no positive ordering is claimed."
            ),
        },
        "runtime_seconds": float(time.perf_counter() - t0),
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "source_sha256": {
            "src/helmholtz.py": sha256(_HERE / "helmholtz.py"),
            "src/family3_gauge_born.py": sha256(_HERE / "family3_gauge_born.py"),
            "src/family3_parent_correction.py": sha256(Path(__file__)),
            "results/family3_gauge_born.json": sha256(original_path),
        },
        "cannot_establish": (
            "Finite grids and local Taylor tests do not prove continuum gauge "
            "theorems, global nonlinear identifiability, or universal anchor "
            "and support-mask orderings."
        ),
    }
    out = _ROOT / "results" / "family3_parent_correction.json"
    out.write_text(json.dumps(result, indent=2) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8))
    styles = {
        "tx": ("o", "#1f77b4"),
        "ty": ("s", "#d62728"),
        "rot": ("^", "#2ca02c"),
    }
    for scene, ls in (("truncated_gaussian", "--"), ("compact_C_infinity", "-")):
        for gen, (marker, color) in styles.items():
            vals = [
                next(x for x in row["scenes"][scene] if x["generator"] == gen)[
                    "gauge_residual"
                ]
                for row in gauge["rows"]
            ]
            axes[0].semilogy(
                GRID_NS, vals, marker=marker, color=color, ls=ls,
                label=f"{scene}, {gen}",
            )
    axes[0].set_title("Pixel gauge residual: boundary vs compact support")
    axes[0].set_xlabel("grid N")
    axes[0].set_ylabel("normalized residual")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=6)

    epss = [r["epsilon"] for r in born["rows"]]
    axes[1].loglog(
        epss,
        [r["joint_Born_coefficient_error_rel"] for r in born["rows"]],
        "o-", label="Born mixed coefficient error",
    )
    axes[1].loglog(
        epss,
        [r["joint_fullwave_coefficient_error_rel"] for r in born["rows"]],
        "s-", label="full-wave weak-contrast error",
    )
    axes[1].set_title("Empty-background mixed term under joint scaling")
    axes[1].set_xlabel("epsilon")
    axes[1].set_ylabel("relative coefficient error")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=7)

    for gen, (_, color) in styles.items():
        axes[2].semilogx(
            [max(r["alpha"], 1e-8) for r in anchor["rows"]],
            [r["directional_retention"][gen] for r in anchor["rows"]],
            "o-", color=color, label=gen,
        )
    axes[2].set_title("First-pose anchor lifts gauge-direction retention")
    axes[2].set_xlabel("anchor prior alpha (zero plotted at 1e-8)")
    axes[2].set_ylabel("directional retention")
    axes[2].grid(True, which="both", alpha=0.25)
    axes[2].legend(fontsize=7)
    fig.suptitle("Family 3 parent correction: gauge semantics and Born bilinear term")
    fig.tight_layout()
    fig_path = _ROOT / "figures" / "family3_parent_correction.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    note = _ROOT / "notes" / "family3_parent_correction.md"
    note.write_text(
        "# Family 3 parent correction\n\n"
        "The original failed gates remain preserved.  The compact-support pixel "
        "scene reduces translation-x gauge residual from %.3e at N=16 to %.3e "
        "at N=40; this is only a finite-grid diagnostic.  The radial scene's "
        "rotation has zero map tangent and an absolute pose residual %.3e "
        "(%.3e of the translation scale), so it is treated as a pose-only "
        "stabilizer rather than a map-retention direction.  The corrected Born "
        "joint signal has slope %.3f (expected 2), its coefficient error has "
        "slope %.3f, and the full-wave weak-contrast discrepancy has slope "
        "%.3f.  A first-pose prior lifts the tested map-gauge directional "
        "retentions; arbitrary support masks have no asserted ordering.\n"
        % (
            next(x for x in gauge["rows"][0]["scenes"]["compact_C_infinity"]
                 if x["generator"] == "tx")["gauge_residual"],
            next(x for x in gauge["rows"][-1]["scenes"]["compact_C_infinity"]
                 if x["generator"] == "tx")["gauge_residual"],
            symmetric_stabilizer["B_dX_rotation_l2"],
            symmetric_stabilizer["B_dX_rotation_over_translation_scale"],
            born["slopes_first_five"]["joint_Born_norm"]["slope"],
            born["slopes_first_five"]["central_D_error"]["slope"],
            born["slopes_first_five"]["fullwave_coefficient_error"]["slope"],
        )
    )

    print("Family 3 parent correction complete")
    print("Born slopes:", born["slopes_first_five"])
    print("symmetric stabilizer:", symmetric_stabilizer)
    print("results ->", out)
    print("figure  ->", fig_path)
    print("note    ->", note)


if __name__ == "__main__":
    main()
