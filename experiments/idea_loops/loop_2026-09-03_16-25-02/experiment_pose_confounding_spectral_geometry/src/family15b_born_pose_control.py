"""Family 15b: complete Born-pair pose-tangent control.

This additive parent correction complements Family 15.  Family 15 deliberately
used (A_born | B_full) to isolate replacement of the map tangent.  Here we
construct and finite-difference validate the actual Born pose Jacobian

    B_born = D_X [G_S(X) diag(chi) e_inc(X)]

and compare the complete pairs (A_full | B_full), (A_born | B_born), and the
mixed isolation pair (A_born | B_full).  It also checks the exact finite-
dimensional contrast-scaling law A_born(s)=A0 and B_born(s)=s B1.

Run from the experiment root:
    .venv/bin/python src/family15b_born_pose_control.py
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import helmholtz as hh  # noqa: E402
import family1_pilot as family1  # noqa: E402
import family2_algebraic_spine as family2  # noqa: E402
import family4_frequency_trajectory as family4  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


CONFIG = {
    "family": "15b",
    "N": 16,
    "T": 6,
    "n_rx": 4,
    "k_b": 2.0 * math.pi,
    "rx_offsets": [[-0.06, 0.0], [0.06, 0.0], [0.0, -0.06], [0.0, 0.06]],
    "tx_offset": [0.0, 0.0],
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],
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
    },
    "contrast_scales": [0.001, 0.01, 0.03, 0.1, 0.3, 1.0],
    "alpha": 1.0,
    "fd_validation_scale": 0.1,
    "fd_steps": [0.001, 0.0003, 0.0001],
    "rank_rule": "max(M.shape) * eps_machine * sigma_1(M)",
    "scope": (
        "finite-dimensional N=16 scalar Helmholtz toy; local tangents only; "
        "no continuum, estimator, hardware, or external-benchmark claim"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sym(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def rel_fro(X: np.ndarray, Y: np.ndarray) -> float:
    return float(np.linalg.norm(X - Y, "fro") / max(np.linalg.norm(Y, "fro"), np.finfo(float).tiny))


def rel_vec(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.linalg.norm(x - y) / max(np.linalg.norm(y), np.finfo(float).tiny))


def build_poses() -> np.ndarray:
    phi = np.deg2rad(np.linspace(CONFIG["arc_phi_deg"][0], CONFIG["arc_phi_deg"][1], CONFIG["T"]))
    p = CONFIG["arc_radius"] * np.column_stack([np.cos(phi), np.sin(phi)])
    theta = np.arctan2(-p[:, 1], -p[:, 0])
    return np.column_stack([p, theta])


def build_born_pair(
    chi: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return exact discrete Born A, pose derivative B, and data F (complex)."""
    N = CONFIG["N"]
    k_b = CONFIG["k_b"]
    points, h = hh.make_grid(N)
    geom = hh.pose_geometry(poses, rx_offsets, tx_offset)
    T = poses.shape[0]
    n_rx = rx_offsets.shape[0]
    n_pix = points.shape[0]
    A = np.empty((T * n_rx, n_pix), dtype=np.complex128)
    B = np.zeros((T * n_rx, 3 * T), dtype=np.complex128)
    F = np.empty(T * n_rx, dtype=np.complex128)

    for t in range(T):
        sl = slice(t * n_rx, (t + 1) * n_rx)
        rx_t = geom["rx_world"][t]
        tx_t = geom["tx_world"][t]
        drx_t = geom["drx_world"][t]
        dtx_t = geom["dtx_world"][t]

        G_S = (k_b**2) * (h**2) * hh.green_matrix(points, rx_t, k_b)
        E_inc = hh.green_matrix(tx_t[None, :], points, k_b)[:, 0]
        J_born = chi * E_inc
        A_t = G_S * E_inc[None, :]
        F_t = G_S @ J_born

        grad_rx = hh.green_grad_first(points, rx_t, k_b)
        grad_tx = hh.green_grad_source(points, tx_t, k_b)
        B_t = np.empty((n_rx, 3), dtype=np.complex128)
        for ell in range(3):
            dG_S = (k_b**2) * (h**2) * np.einsum(
                "and,ad->an", grad_rx, drx_t[ell]
            )
            dE_inc = np.einsum("nd,d->n", grad_tx, dtx_t[ell])
            B_t[:, ell] = dG_S @ J_born + G_S @ (chi * dE_inc)

        A[sl] = A_t
        B[sl, 3 * t : 3 * t + 3] = B_t
        F[sl] = F_t

    return A, B, F


def born_pose_fd(
    chi: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    step: float,
) -> np.ndarray:
    """Centered finite-difference Jacobian of the Born data with respect to all poses."""
    T = poses.shape[0]
    rows = T * rx_offsets.shape[0]
    B = np.empty((rows, 3 * T), dtype=np.complex128)
    flat = poses.ravel()
    for j in range(flat.size):
        plus = flat.copy()
        minus = flat.copy()
        plus[j] += step
        minus[j] -= step
        Fp, _ = hh.born_forward(
            chi, plus.reshape(T, 3), rx_offsets, tx_offset, CONFIG["N"], CONFIG["k_b"]
        )
        Fm, _ = hh.born_forward(
            chi, minus.reshape(T, 3), rx_offsets, tx_offset, CONFIG["N"], CONFIG["k_b"]
        )
        B[:, j] = (Fp - Fm) / (2.0 * step)
    return B


def data_projector(B: np.ndarray) -> tuple[int, np.ndarray, np.ndarray]:
    r, Z, sv, _ = family2.range_basis(B)
    P_B = Z @ Z.T
    return r, P_B, sv


def retention(A: np.ndarray, W: np.ndarray) -> dict:
    d = family2.thin_decomposition(A)
    rho = family2.retention_spectrum(d["Q"], W)
    return {
        "rank_A": int(d["rank"]),
        "rho_desc": rho.tolist(),
        "retained_dof": float(np.sum(rho)),
        "rho_min": float(rho[-1]),
        "rho_max": float(rho[0]),
    }


def pair_metrics(A: np.ndarray, B: np.ndarray, alpha: float) -> dict:
    rB, P_B, svB = data_projector(B)
    P_perp = np.eye(A.shape[0]) - P_B
    W = family4.shrinkage_operator(B, alpha)
    return {
        "rank_B": int(rB),
        "sigma_B_desc": svB.tolist(),
        "no_prior": retention(A, P_perp),
        "finite_prior": retention(A, W),
        "K_IS": sym(A.T @ A),
        "K_SLAM": sym(A.T @ (P_perp @ A)),
        "K_eff": sym(A.T @ (W @ A)),
        "P_B": P_B,
    }


def fit_slope(scales: np.ndarray, values: np.ndarray) -> dict:
    x = np.log10(scales)
    y = np.log10(values)
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(1.0 - ss_res / ss_tot) if ss_tot else 1.0,
    }


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items() if k not in {"K_IS", "K_SLAM", "K_eff", "P_B"}}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def main() -> None:
    start = time.perf_counter()
    points, _ = hh.make_grid(CONFIG["N"])
    poses = build_poses()
    rx = np.asarray(CONFIG["rx_offsets"], dtype=float)
    tx = np.asarray(CONFIG["tx_offset"], dtype=float)
    chi0 = family1.make_chi0(points, CONFIG)
    S = family2.build_smooth_basis(points, CONFIG["smooth_basis"])
    alpha = float(CONFIG["alpha"])
    scales = np.asarray(CONFIG["contrast_scales"], dtype=float)

    rows = []
    born_B1 = None
    born_P_ref = None
    for s in scales:
        chi = float(s) * chi0
        A_full_c, B_full_c, _, _ = hh.build_AB(
            chi, poses, rx, tx, CONFIG["N"], CONFIG["k_b"]
        )
        A_born_c, B_born_c, F_born = build_born_pair(chi, poses, rx, tx)
        F_ref, A_ref_c = hh.born_forward(
            chi, poses, rx, tx, CONFIG["N"], CONFIG["k_b"]
        )
        A_full_R, B_full_R = hh.whiten_realify(A_full_c, B_full_c, None)
        A_born_R, B_born_R = hh.whiten_realify(A_born_c, B_born_c, None)
        A_full = A_full_R @ S
        A_born = A_born_R @ S

        full = pair_metrics(A_full, B_full_R, alpha)
        born = pair_metrics(A_born, B_born_R, alpha)
        mixed = pair_metrics(A_born, B_full_R, alpha)

        B1 = B_born_R / float(s)
        if born_B1 is None:
            born_B1 = B1
            born_P_ref = born["P_B"]
        KIS_full = full["K_IS"]
        KIS_born = born["K_IS"]

        loss = KIS_born - born["K_eff"]
        cross = A_born.T @ born_B1
        loss_bound = float((s**2 / alpha) * np.linalg.norm(cross, ord=2) ** 2)
        loss_norm = float(np.linalg.norm(loss, ord=2))

        rows.append(
            {
                "s": float(s),
                "rank_B_full": full["rank_B"],
                "rank_B_born": born["rank_B"],
                "born_forward_rel_error": rel_vec(F_born, F_ref),
                "born_A_path_rel_error": rel_fro(A_born_c, A_ref_c),
                "born_A_scale_invariance_rel": rel_fro(A_born, rows[0]["_A_born"] if rows else A_born),
                "born_B_over_s_invariance_rel": rel_fro(B1, born_B1),
                "born_projector_invariance_op": float(np.linalg.norm(born["P_B"] - born_P_ref, ord=2)),
                "full_vs_born": {
                    "rel_A": rel_fro(A_born, A_full),
                    "rel_B": rel_fro(B_born_R, B_full_R),
                    "rel_K_IS": rel_fro(KIS_born, KIS_full),
                    "rel_K_eff_complete_pairs": rel_fro(born["K_eff"], full["K_eff"]),
                    "rel_K_SLAM_complete_pairs": rel_fro(born["K_SLAM"], full["K_SLAM"]),
                    "pose_projector_gap_op": float(np.linalg.norm(full["P_B"] - born["P_B"], ord=2)),
                },
                "full_pair": jsonable(full),
                "born_pair": jsonable(born),
                "mixed_pair_A_born_B_full": jsonable(mixed),
                "finite_prior_born_loss": {
                    "spectral_norm": loss_norm,
                    "O_s2_upper_bound_alpha_I": loss_bound,
                    "bound_holds": bool(loss_norm <= loss_bound * (1.0 + 1e-12) + 1e-20),
                    "ratio_bound_over_actual": float(loss_bound / max(loss_norm, np.finfo(float).tiny)),
                },
                "_A_born": A_born,
            }
        )

    # Validate the analytic Born pose Jacobian against centered FD.
    s_fd = float(CONFIG["fd_validation_scale"])
    chi_fd = s_fd * chi0
    _, B_exact_c, _ = build_born_pair(chi_fd, poses, rx, tx)
    fd_rows = []
    for hfd in CONFIG["fd_steps"]:
        B_fd_c = born_pose_fd(chi_fd, poses, rx, tx, float(hfd))
        fd_rows.append({"h": float(hfd), "relative_fro_error": rel_fro(B_fd_c, B_exact_c)})
    fd_slope = fit_slope(
        np.asarray([r["h"] for r in fd_rows]),
        np.asarray([r["relative_fro_error"] for r in fd_rows]),
    )

    for row in rows:
        row.pop("_A_born", None)

    fits = {}
    for key in ["rel_A", "rel_B", "rel_K_IS", "rel_K_eff_complete_pairs", "rel_K_SLAM_complete_pairs", "pose_projector_gap_op"]:
        vals = np.asarray([r["full_vs_born"][key] for r in rows])
        fits[key] = fit_slope(scales, vals)

    # Exact s=0 Born endpoint: B=0, so no-prior elimination equals known pose.
    A0_c, B0_c, _ = build_born_pair(np.zeros_like(chi0), poses, rx, tx)
    A0_R, B0_R = hh.whiten_realify(A0_c, B0_c, None)
    zero_pair = pair_metrics(A0_R @ S, B0_R, alpha)

    # Figure.
    fig_path = ROOT / "figures" / "family15b_born_pose_control.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0), constrained_layout=True)
    ax = axes[0, 0]
    for key, label in [
        ("rel_A", "A map tangent"),
        ("rel_B", "B pose tangent"),
        ("rel_K_IS", "K_IS"),
        ("rel_K_eff_complete_pairs", "K_eff complete pair"),
    ]:
        ax.loglog(scales, [r["full_vs_born"][key] for r in rows], "o-", label=f"{label}; slope {fits[key]['slope']:.2f}")
    ax.set_xlabel("contrast scale s")
    ax.set_ylabel("relative Frobenius error")
    ax.set_title("(1) Complete Born-pair versus full-wave tangents")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.semilogx(scales, [r["full_pair"]["no_prior"]["retained_dof"] for r in rows], "o-", label="full pair")
    ax.semilogx(scales, [r["born_pair"]["no_prior"]["retained_dof"] for r in rows], "s--", label="complete Born pair")
    ax.semilogx(scales, [r["mixed_pair_A_born_B_full"]["no_prior"]["retained_dof"] for r in rows], "^:", label="mixed A_born | B_full")
    ax.axhline(zero_pair["no_prior"]["retained_dof"], color="black", linestyle="--", linewidth=1.0, label="Born endpoint s=0 (B=0)")
    ax.set_xlabel("contrast scale s (s>0)")
    ax.set_ylabel("no-prior retained DOF")
    ax.set_title("(2) Scale-invariant Born range and singular s=0 endpoint")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    actual = [r["finite_prior_born_loss"]["spectral_norm"] for r in rows]
    bound = [r["finite_prior_born_loss"]["O_s2_upper_bound_alpha_I"] for r in rows]
    ax.loglog(scales, actual, "o-", label="actual ||K_IS-K_eff||_2")
    ax.loglog(scales, bound, "s--", label="s^2 ||A^T B1||_2^2 / alpha")
    ax.set_xlabel("contrast scale s")
    ax.set_ylabel("spectral norm")
    ax.set_title("(3) Finite-prior Born loss and O(s^2) bound")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.loglog(scales, [r["full_vs_born"]["pose_projector_gap_op"] for r in rows], "o-", label="||P_Bfull-P_Bborn||_2")
    ax.loglog(scales, [max(r["born_B_over_s_invariance_rel"], 1e-18) for r in rows], "s--", label="B_born/s invariance residual")
    ax.loglog(scales, [max(r["born_projector_invariance_op"], 1e-18) for r in rows], "^:", label="P_Bborn invariance residual")
    ax.set_xlabel("contrast scale s")
    ax.set_ylabel("operator/relative residual")
    ax.set_title("(4) Nuisance-range convergence and exact Born scaling")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle("Family 15b: complete Born pose-tangent control (finite N=16)", fontsize=14)
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    runtime = time.perf_counter() - start
    result_path = ROOT / "results" / "family15b_born_pose_control.json"
    note_path = ROOT / "notes" / "family15b_born_pose_control.md"
    result = {
        "schema": "family15b_born_pose_control.v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "runner": sys.executable,
        "command": ".venv/bin/python src/family15b_born_pose_control.py",
        "wall_runtime_seconds": runtime,
        "platform": platform.platform(),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "config": CONFIG,
        "analytic_formula": (
            "B_born,t[:,ell] = (D_X G_S,t[ell]) (chi * E_inc,t) + "
            "G_S,t (chi * D_X E_inc,t[ell])"
        ),
        "finite_dimensional_proposition": {
            "statement": (
                "For F_B(chi,X)=A0(X)chi and chi=s*chi_bar, the map tangent is "
                "A_B(s)=A0 while the pose tangent is B_B(s)=s*B1. Thus for "
                "every s!=0 the no-prior range projector and K_SLAM are scale "
                "invariant, but at s=0 B_B=0 and K_SLAM=K_IS; the endpoint is "
                "discontinuous unless Range(A0) is orthogonal to Range(B1). "
                "For J=alpha I, alpha>0, the loss is continuous and bounded by "
                "s^2 ||A0^T B1||_2^2/alpha."
            ),
            "proof": (
                "Differentiate the bilinear map A0(X)(s chi_bar): D_chi F=A0 "
                "and D_X F=s D_X[A0(X)chi_bar]=sB1. Nonzero scalar multiplication "
                "does not change a range projector. At s=0 the pose tangent is "
                "zero. For alpha I, use (s^2 B1^T B1+alpha I)^-1 <= alpha^-1 I "
                "inside the Schur loss."
            ),
            "semidefinite_prior_caveat": (
                "The O(s^2) bound requires a strictly positive prior floor on the "
                "relevant pose directions; a singular prior needs a separate range/kernel analysis."
            ),
        },
        "born_pose_fd_validation": {"scale": s_fd, "rows": fd_rows, "fit": fd_slope},
        "rows": rows,
        "convergence_fits": fits,
        "born_s_zero_endpoint": jsonable(zero_pair),
        "gates": {
            "born_pose_fd_smallest_step_rel_lt_1e_minus_6": bool(fd_rows[-1]["relative_fro_error"] < 1e-6),
            "born_pose_fd_all_rel_lt_2e_minus_5": bool(max(r["relative_fro_error"] for r in fd_rows) < 2e-5),
            "born_pose_fd_slope_near_2": bool(1.8 < fd_slope["slope"] < 2.2),
            "A_scale_invariance_max_lt_1e_minus_12": bool(max(r["born_A_scale_invariance_rel"] for r in rows) < 1e-12),
            "B_over_s_invariance_max_lt_1e_minus_12": bool(max(r["born_B_over_s_invariance_rel"] for r in rows) < 1e-12),
            "projector_invariance_max_lt_1e_minus_10": bool(max(r["born_projector_invariance_op"] for r in rows) < 1e-10),
            "finite_prior_bounds_all_hold": bool(all(r["finite_prior_born_loss"]["bound_holds"] for r in rows)),
        },
        "literature": {
            "citation": (
                "M. L. Diong, A. Roueff, P. Lasaygues, and A. Litman, "
                "Impact of the Born approximation on the estimation error in 2D inverse scattering, "
                "Inverse Problems 32(6), 065006 (2016), DOI 10.1088/0266-5611/32/6/065006."
            ),
            "boundary": (
                "Conceptual/statistical neighbor only. This run does not reproduce its geometry, "
                "parameterization, noise, estimator bias, units, or numerical values."
            ),
        },
        "scope_note": CONFIG["scope"],
        "source_sha256": sha256(Path(__file__)),
        "figure_sha256": sha256(fig_path),
    }
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    first = rows[0]
    last = rows[-1]
    note = f"""# Family 15b: complete Born pose-tangent control

This additive parent correction closes the semantic gap in Family 15. Family 15's `(A_born | B_full)` rows remain a useful **mixed-tangent isolation control**; they are not a complete Born-SLAM Fisher model. Family 15b constructs the actual analytic Born pose Jacobian and compares complete tangent pairs.

## Execution

- Command: `.venv/bin/python src/family15b_born_pose_control.py`
- Runtime: {runtime:.3f} s on Apple Silicon CPU ({platform.platform()})
- Dimensions: N=16, 48 real data rows, p=24 smooth map coefficients, q=18 pose coordinates
- Scope: {CONFIG['scope']}

## Analytic Born pose derivative

For each pose and coordinate, the implemented exact discrete formula is

`B_born,t[:,ell] = (D_X G_S,t[ell]) (chi * E_inc,t) + G_S,t (chi * D_X E_inc,t[ell])`.

At contrast scale s={s_fd:g}, centered finite differences over all 18 pose coordinates give:

| h | relative Frobenius error |
|---:|---:|
""" + "\n".join(f"| {r['h']:.1e} | {r['relative_fro_error']:.6e} |" for r in fd_rows) + f"""

The observed error-vs-step slope is {fd_slope['slope']:.4f} (R^2={fd_slope['r_squared']:.6f}), supporting the analytic derivative with the expected centered-difference second-order trend on this discrete model.

## Exact finite-dimensional contrast-scaling proposition

For the linear Born map `F_B(chi,X)=A0(X) chi`, set `chi=s chi_bar`. Then `A_B(s)=A0` and `B_B(s)=s B1`. Hence, for every nonzero s, `Range(B_B(s))=Range(B1)`: the no-prior projector and Born `K_SLAM` are independent of contrast amplitude. At s=0, `B_B(0)=0`, so `K_SLAM(0)=K_IS`. The endpoint is discontinuous unless the map tangent is orthogonal to the nonzero Born pose range.

With `J=alpha I`, alpha>0, the finite-prior loss is continuous and obeys

`||K_IS-K_eff(s)||_2 <= s^2 ||A0^T B1||_2^2 / alpha`.

This follows directly from differentiation of the bilinear Born map, invariance of a range under nonzero scalar multiplication, and `(s^2 B1^T B1+alpha I)^-1 <= alpha^-1 I`. A semidefinite prior needs a separate range/kernel analysis.

Executed checks: max `A_B(s)` scale residual {max(r['born_A_scale_invariance_rel'] for r in rows):.3e}; max `B_B(s)/s` residual {max(r['born_B_over_s_invariance_rel'] for r in rows):.3e}; max nonzero-s Born projector residual {max(r['born_projector_invariance_op'] for r in rows):.3e}. All finite-prior bound rows hold: {all(r['finite_prior_born_loss']['bound_holds'] for r in rows)}.

The nonzero-s complete-Born retained DOF is {first['born_pair']['no_prior']['retained_dof']:.8f} at s={first['s']:g} and {last['born_pair']['no_prior']['retained_dof']:.8f} at s={last['s']:g}; the s=0 endpoint is {zero_pair['no_prior']['retained_dof']:.8f}. The tiny nonzero-s variation is numerical. This is a rank-stratum singular limit, not evidence that arbitrarily weak scatterers provide finite practical pose information.

## Complete-pair versus mixed-pair result

| s | rel A | rel B | rel K_IS | rel K_eff complete | projector gap | rdof full | rdof Born complete | rdof mixed |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
""" + "\n".join(
        f"| {r['s']:.3g} | {r['full_vs_born']['rel_A']:.3e} | {r['full_vs_born']['rel_B']:.3e} | {r['full_vs_born']['rel_K_IS']:.3e} | {r['full_vs_born']['rel_K_eff_complete_pairs']:.3e} | {r['full_vs_born']['pose_projector_gap_op']:.3e} | {r['full_pair']['no_prior']['retained_dof']:.6f} | {r['born_pair']['no_prior']['retained_dof']:.6f} | {r['mixed_pair_A_born_B_full']['no_prior']['retained_dof']:.6f} |"
        for r in rows
    ) + f"""

Observed log-log slopes are A {fits['rel_A']['slope']:.3f}, B {fits['rel_B']['slope']:.3f}, K_IS {fits['rel_K_IS']['slope']:.3f}, complete-pair K_eff {fits['rel_K_eff_complete_pairs']['slope']:.3f}, and pose-projector gap {fits['pose_projector_gap_op']['slope']:.3f}. These are finite-grid observations, not continuum asymptotics.

## Scientific boundary

- Only the complete `(A_born | B_born)` rows represent the implemented Born pose-elimination model.
- `(A_born | B_full)` remains explicitly labelled mixed-tangent.
- At s=0, first-order pose information vanishes; for s>0, a no-prior projector ignores the magnitude of B and is therefore a singular statistical idealization as s tends to zero.
- The Diong et al. paper (DOI `10.1088/0266-5611/32/6/065006`) is a conceptual neighbor only. This run does not reproduce its estimator bias, geometry, noise, parameterization, units, or numerical values and is not an external benchmark.

## Artifacts

- `results/family15b_born_pose_control.json`
- `figures/family15b_born_pose_control.png`
- `src/family15b_born_pose_control.py`
"""
    note_path.write_text(note, encoding="utf-8")
    print(json.dumps({
        "status": "success",
        "runtime_seconds": runtime,
        "fd_slope": fd_slope["slope"],
        "born_nonzero_s_rdof_first_last": [first["born_pair"]["no_prior"]["retained_dof"], last["born_pair"]["no_prior"]["retained_dof"]],
        "born_s0_rdof": zero_pair["no_prior"]["retained_dof"],
        "gates": result["gates"],
        "artifacts": [str(result_path), str(fig_path), str(note_path)],
    }, indent=2))


if __name__ == "__main__":
    main()
