"""Parent audit supplement for Family 5 generalized retention derivatives.

The original Family 5 correctly differentiates ordinary eigenvalues of the
finite-prior K_eff matrix.  The paper's central normalized quantity is instead
the generalized retention eigenvalue

    K_eff(X) v(X) = rho(X) K_IS(X) v(X),
    v(X)^T K_IS(X) v(X) = 1.

This supplement verifies the missing ``-rho * dot(K_IS)`` term for both the
finite-prior and no-prior (constant-rank) cases, treats the exact rho=1 cluster
with a compressed derivative, and constructs an explicit rank-loss event at
which the no-prior projector is discontinuous.  It is a finite-dimensional
diagnostic, not a continuum or physical-trajectory theorem.

Run from the experiment root:
    .venv/bin/python src/family5_parent_generalized.py
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

import family5_sensitivity as f5  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


EPS = float(np.finfo(float).eps)


def sym(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def rank_tol(M: np.ndarray) -> tuple[int, np.ndarray, float]:
    s = np.linalg.svd(M, compute_uv=False)
    tol = max(M.shape) * EPS * (float(s[0]) if s.size else 0.0)
    return int(np.sum(s > tol)), s, tol


def range_projector(B: np.ndarray) -> tuple[np.ndarray, int, np.ndarray, float]:
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    tol = max(B.shape) * EPS * (float(s[0]) if s.size else 0.0)
    r = int(np.sum(s > tol))
    return U[:, :r] @ U[:, :r].T, r, s, tol


def information_matrices(
    A: np.ndarray, B: np.ndarray, alpha: float | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    KIS = sym(A.T @ A)
    if alpha is None:
        PB, _, _, _ = range_projector(B)
        E = np.eye(A.shape[0]) - PB
    else:
        C = B.T @ B + float(alpha) * np.eye(B.shape[1])
        E = np.eye(A.shape[0]) - B @ np.linalg.solve(C, B.T)
    K = sym(A.T @ E @ A)
    return KIS, K, sym(E)


def retention_decomposition(
    A: np.ndarray, B: np.ndarray, alpha: float | None
) -> dict:
    """Stable generalized eigendecomposition through the SVD of A.

    If A=Q Sigma V^T has full column rank, eigenpairs of Q^T E Q map to
    coefficient vectors v=V Sigma^{-1} y with v^T A^T A v=1.
    """
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    tol = max(A.shape) * EPS * float(s[0])
    r = int(np.sum(s > tol))
    if r != A.shape[1]:
        raise RuntimeError(
            f"generalized audit requires full column rank A; got {r}/{A.shape[1]}"
        )
    Q = U[:, :r]
    KIS, K, E = information_matrices(A, B, alpha)
    R = sym(Q.T @ E @ Q)
    rho, Y = np.linalg.eigh(R)
    Vcoef = Vh.T[:, :r] @ ((1.0 / s[:r])[:, None] * Y)
    D = A @ Vcoef
    # The direct coefficient-space checks below are deliberately retained as
    # conditioning diagnostics.  The scientifically relevant backward error
    # is evaluated in the factored R problem, avoiding cancellation through
    # cond(A)^2 in KIS=A^T A.
    norm_err_direct = float(np.max(np.abs(np.diag(Vcoef.T @ KIS @ Vcoef) - 1.0)))
    norm_err_factored = float(np.max(np.abs(np.diag(D.T @ D) - 1.0)))
    residuals_direct = []
    residuals_factored = []
    residuals_relative_eigenvalue = []
    Rnorm = float(np.linalg.norm(R, ord=2))
    for j in range(r):
        v = Vcoef[:, j]
        num = np.linalg.norm(K @ v - rho[j] * (KIS @ v))
        den = max(np.linalg.norm(K @ v) + abs(rho[j]) * np.linalg.norm(KIS @ v), EPS)
        residuals_direct.append(float(num / den))
        y = Y[:, j]
        rnum = np.linalg.norm(R @ y - rho[j] * y)
        # Normwise backward residual.  Scaling only by |rho| would turn an
        # O(eps) absolute residual for the near-zero modes into a misleading
        # O(1e-6) number.
        rden = max(Rnorm + abs(rho[j]), EPS)
        residuals_factored.append(float(rnum / rden))
        residuals_relative_eigenvalue.append(float(rnum / max(abs(rho[j]), EPS)))
    return {
        "rho": rho,
        "Vcoef": Vcoef,
        "data_dirs": D,
        "KIS": KIS,
        "K": K,
        "E": E,
        "rank_A": r,
        "singular_A": s,
        "tol_A": tol,
        "normalization_direct_coefficient_max_error": norm_err_direct,
        "normalization_factored_data_max_error": norm_err_factored,
        "generalized_residual_direct_coefficient_max": float(max(residuals_direct)),
        "generalized_residual_factored_max": float(max(residuals_factored)),
        "generalized_residual_relative_eigenvalue_max": float(max(residuals_relative_eigenvalue)),
    }


def match_mode(base_data_dir: np.ndarray, dec: dict) -> tuple[int, float, float]:
    overlaps = np.abs(base_data_dir @ dec["data_dirs"])
    j = int(np.argmax(overlaps))
    return j, float(dec["rho"][j]), float(overlaps[j])


def choose_simple_modes(rho: np.ndarray, n: int = 3) -> list[int]:
    gaps = np.minimum(
        np.r_[np.inf, np.diff(rho)], np.r_[np.diff(rho), np.inf]
    )
    order = np.argsort(-gaps, kind="stable")
    chosen = []
    for j in order:
        if gaps[j] <= 1e-8:
            continue
        chosen.append(int(j))
        if len(chosen) == n:
            break
    if len(chosen) < n:
        raise RuntimeError("fewer than three well-gapped retention modes")
    return chosen


def derivative_matrices(
    A: np.ndarray,
    B: np.ndarray,
    dA: np.ndarray,
    dB: np.ndarray,
    alpha: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    dKIS = sym(dA.T @ A + A.T @ dA)
    if alpha is not None:
        dK, _, E, _, dE, _ = f5.derivative_terms(A, B, dA, dB, alpha)
        return dKIS, sym(dK), sym(E), sym(dE), {
            "case": "finite_prior",
            "dE_fro": float(np.linalg.norm(dE, ord="fro")),
            "E_fro": float(np.linalg.norm(E, ord="fro")),
        }

    # Constant-rank derivative of P_B = B B^dagger and K=A^T(I-P_B)A.
    PB, rB, sB, tolB = range_projector(B)
    Bdag = np.linalg.pinv(B, rcond=tolB / max(float(sB[0]), EPS))
    IminusP = np.eye(B.shape[0]) - PB
    dPB = sym(IminusP @ dB @ Bdag + Bdag.T @ dB.T @ IminusP)
    dK = sym(
        dA.T @ IminusP @ A
        + A.T @ IminusP @ dA
        - A.T @ dPB @ A
    )
    return dKIS, dK, sym(IminusP), sym(-dPB), {
        "case": "no_prior_constant_rank",
        "rank_B": rB,
        "sigma_min_B": float(sB[rB - 1]),
        "tol_B": float(tolB),
        "dPB_fro": float(np.linalg.norm(dPB, ord="fro")),
    }


def audit_case(
    chi0: np.ndarray,
    S: np.ndarray,
    X0: np.ndarray,
    dX: np.ndarray,
    alpha: float | None,
    eps_grid: list[float],
) -> dict:
    cfg = f5.CONFIG
    A0, B0 = f5.smooth_AB_of_flat_X(chi0, S, X0, cfg)
    dA, dB = f5.smooth_AB_derivatives(chi0, S, X0, dX, cfg)
    dKIS, dK, E, dE, derivative_meta = derivative_matrices(A0, B0, dA, dB, alpha)
    base = retention_decomposition(A0, B0, alpha)
    chosen = choose_simple_modes(base["rho"], 3)

    modes = []
    for j0 in chosen:
        rho0 = float(base["rho"][j0])
        v0 = base["Vcoef"][:, j0]
        d0 = base["data_dirs"][:, j0]
        data0 = A0 @ v0
        ddata0 = dA @ v0
        drho_stable = float(
            ddata0 @ ((E - rho0 * np.eye(E.shape[0])) @ data0)
            + data0 @ ((E - rho0 * np.eye(E.shape[0])) @ ddata0)
            + data0 @ (dE @ data0)
        )
        drho_direct = float(v0 @ ((dK - rho0 * dKIS) @ v0))
        drho_wrong_missing_denominator = float(
            ddata0 @ (E @ data0)
            + data0 @ (E @ ddata0)
            + data0 @ (dE @ data0)
        )
        drho = drho_stable
        left = float(base["rho"][j0] - base["rho"][j0 - 1]) if j0 > 0 else None
        right = float(base["rho"][j0 + 1] - base["rho"][j0]) if j0 + 1 < len(base["rho"]) else None
        min_gap = min(x for x in (left, right) if x is not None)
        rows = []
        for eps in eps_grid:
            Ap, Bp = f5.smooth_AB_of_flat_X(chi0, S, X0 + eps * dX, cfg)
            Am, Bm = f5.smooth_AB_of_flat_X(chi0, S, X0 - eps * dX, cfg)
            dp = retention_decomposition(Ap, Bp, alpha)
            dm = retention_decomposition(Am, Bm, alpha)
            jp, rp, op = match_mode(d0, dp)
            jm, rm, om = match_mode(d0, dm)
            fd = (rp - rm) / (2.0 * eps)
            pred = rho0 + eps * drho
            pred_wrong = rho0 + eps * drho_wrong_missing_denominator
            rows.append(
                {
                    "eps": eps,
                    "rho_plus": rp,
                    "rho_minus": rm,
                    "matched_index_plus": jp,
                    "matched_index_minus": jm,
                    "data_overlap_plus": op,
                    "data_overlap_minus": om,
                    "central_fd_derivative": float(fd),
                    "relative_derivative_error": float(
                        abs(fd - drho) / max(abs(drho), 1e-14)
                    ),
                    "first_order_abs_error": float(abs(rp - pred)),
                    "missing_dKIS_abs_error": float(abs(rp - pred_wrong)),
                }
            )
        # Fit the truncation-dominated largest four points.  The complete rows
        # remain in JSON so roundoff-floor selection is auditable.
        xs = np.log10(np.asarray([r["eps"] for r in rows[-4:]]))
        ys = np.log10(np.maximum(
            np.asarray([r["first_order_abs_error"] for r in rows[-4:]]), 1e-300
        ))
        slope = float(np.polyfit(xs, ys, 1)[0])
        modes.append(
            {
                "index_ascending": j0,
                "rho0": rho0,
                "gap_left": left,
                "gap_right": right,
                "min_adjacent_gap": min_gap,
                "drho_formula": drho,
                "drho_direct_coefficient_formula": drho_direct,
                "stable_vs_direct_abs_difference": abs(drho_stable - drho_direct),
                "drho_wrong_if_dKIS_omitted": drho_wrong_missing_denominator,
                "first_order_error_slope_largest4": slope,
                "rows": rows,
            }
        )

    at_1e3 = [
        min(m["rows"], key=lambda row: abs(row["eps"] - 1e-3)) for m in modes
    ]
    return {
        "alpha": alpha,
        "base_rho_ascending": [float(x) for x in base["rho"]],
        "rank_A": base["rank_A"],
        "condition_A": float(base["singular_A"][0] / base["singular_A"][-1]),
        "normalization_direct_coefficient_max_error": base["normalization_direct_coefficient_max_error"],
        "normalization_factored_data_max_error": base["normalization_factored_data_max_error"],
        "generalized_residual_direct_coefficient_max": base["generalized_residual_direct_coefficient_max"],
        "generalized_residual_factored_max": base["generalized_residual_factored_max"],
        "generalized_residual_relative_eigenvalue_max": base["generalized_residual_relative_eigenvalue_max"],
        "derivative_meta": derivative_meta,
        "chosen_modes": modes,
        "gate": {
            "all_simple_gaps_gt_1e-8": bool(all(m["min_adjacent_gap"] > 1e-8 for m in modes)),
            "all_fd_derivative_relerr_at_1e-3_lt_1e-4": bool(
                all(r["relative_derivative_error"] < 1e-4 for r in at_1e3)
            ),
            "all_prediction_slopes_approx_2": bool(
                all(1.75 <= m["first_order_error_slope_largest4"] <= 2.25 for m in modes)
            ),
            "factored_generalized_residual_lt_1e10": bool(base["generalized_residual_factored_max"] < 1e-10),
        },
    }


def cluster_audit(
    chi0: np.ndarray, S: np.ndarray, X0: np.ndarray, dX: np.ndarray
) -> dict:
    """Compressed derivative for the exact no-prior rho=1 cluster."""
    cfg = f5.CONFIG
    A0, B0 = f5.smooth_AB_of_flat_X(chi0, S, X0, cfg)
    dA, dB = f5.smooth_AB_derivatives(chi0, S, X0, dX, cfg)
    dKIS, dK, E, dE, _ = derivative_matrices(A0, B0, dA, dB, None)
    dec = retention_decomposition(A0, B0, None)
    cluster = np.flatnonzero(np.abs(dec["rho"] - 1.0) < 1e-10)
    Vc = dec["Vcoef"][:, cluster]
    H_direct = sym(Vc.T @ (dK - dKIS) @ Vc)
    D = A0 @ Vc
    dD = dA @ Vc
    H = sym(
        dD.T @ ((E - np.eye(E.shape[0])) @ D)
        + D.T @ ((E - np.eye(E.shape[0])) @ dD)
        + D.T @ dE @ D
    )
    mu = np.linalg.eigvalsh(H)
    eps_rows = []
    for eps in (1e-4, 1e-3, 1e-2):
        Ap, Bp = f5.smooth_AB_of_flat_X(chi0, S, X0 + eps * dX, cfg)
        dp = retention_decomposition(Ap, Bp, None)
        top = dp["rho"][-len(cluster):]
        eps_rows.append({
            "eps": eps,
            "top_cluster_rho": [float(x) for x in top],
            "max_abs_rho_minus_one": float(np.max(np.abs(top - 1.0))),
        })
    return {
        "cluster_indices": [int(x) for x in cluster],
        "cluster_size": int(len(cluster)),
        "compressed_derivative_eigenvalues": [float(x) for x in mu],
        "compressed_derivative_spectral_norm": float(np.linalg.norm(H, ord=2)),
        "direct_coefficient_compressed_spectral_norm": float(np.linalg.norm(H_direct, ord=2)),
        "perturbed_rows": eps_rows,
        "interpretation": (
            "Repeated generalized eigenvalues require the compressed matrix; "
            "individual-vector derivatives are not invariant.  Here the six "
            "structural rho=1 modes remain at one under constant ranks."
        ),
    }


def engineered_rank_event(A: np.ndarray, B: np.ndarray) -> dict:
    """Make the weakest singular value of B cross zero at t=0."""
    U, s, Vh = np.linalg.svd(B, full_matrices=False)
    r0, _, tol0 = rank_tol(B)
    if r0 != B.shape[1]:
        raise RuntimeError("base B is not full column rank")
    rows = []
    projectors = {}
    for t in (-1e-2, -1e-4, -1e-8, 0.0, 1e-8, 1e-4, 1e-2):
        st = s.copy()
        st[-1] = abs(t) * s[-1]
        Bt = (U * st[None, :]) @ Vh
        PB, rB, svB, tolB = range_projector(Bt)
        dec = retention_decomposition(A, Bt, None)
        projectors[t] = PB
        rows.append({
            "t": t,
            "rank_B": rB,
            "sigma_min_B": float(svB[-1]),
            "tol_B": float(tolB),
            "rho_min": float(dec["rho"][0]),
            "rho_max": float(dec["rho"][-1]),
            "retained_mass": float(np.sum(dec["rho"])),
        })
    KIS, Kleft, _ = information_matrices(A, (U * np.r_[s[:-1], 1e-8 * s[-1]][None, :]) @ Vh, None)
    _, Kzero, _ = information_matrices(A, (U * np.r_[s[:-1], 0.0][None, :]) @ Vh, None)
    jump = Kzero - Kleft
    return {
        "construction": "B(t)=U diag(s1,...,s17,|t| s18) V^T",
        "base_rank_B": r0,
        "base_sigma_min_B": float(s[-1]),
        "base_tol_B": float(tol0),
        "rows": rows,
        "projector_jump_norm_at_zero": float(np.linalg.norm(projectors[0.0] - projectors[1e-8], ord=2)),
        "KSLAM_jump_fro_at_zero": float(np.linalg.norm(jump, ord="fro")),
        "KSLAM_jump_relative_to_KIS": float(
            np.linalg.norm(jump, ord="fro") / np.linalg.norm(KIS, ord="fro")
        ),
        "interpretation": (
            "For every nonzero t above the numerical rank threshold the range "
            "of B(t) is unchanged, but at t=0 it loses one dimension.  The "
            "orthogonal projector therefore jumps; no smooth constant-rank "
            "derivative formula applies at the event."
        ),
    }


def make_figure(results: dict, path: Path) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    colors = plt.cm.tab10.colors
    for case_i, key in enumerate(("finite_prior", "no_prior")):
        case = results[key]
        for i, mode in enumerate(case["chosen_modes"]):
            eps = np.asarray([r["eps"] for r in mode["rows"]])
            err = np.asarray([r["relative_derivative_error"] for r in mode["rows"]])
            axes[0].loglog(
                eps, err, "o-", color=colors[case_i * 3 + i],
                label=f"{key}, idx {mode['index_ascending']}"
            )
            pred = np.asarray([r["first_order_abs_error"] for r in mode["rows"]])
            wrong = np.asarray([r["missing_dKIS_abs_error"] for r in mode["rows"]])
            axes[1].loglog(eps, pred, "o-", color=colors[case_i * 3 + i])
            axes[1].loglog(eps, wrong, ":", color=colors[case_i * 3 + i], alpha=0.7)
    axes[0].set_title("Generalized derivative vs centered FD")
    axes[0].set_xlabel("trajectory perturbation epsilon")
    axes[0].set_ylabel("relative derivative error")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=7)
    axes[1].set_title("First-order error (solid) vs omit dKIS (dotted)")
    axes[1].set_xlabel("trajectory perturbation epsilon")
    axes[1].set_ylabel("absolute rho prediction error")
    axes[1].grid(True, which="both", alpha=0.25)

    rows = results["engineered_rank_event"]["rows"]
    x = np.arange(len(rows))
    mass = [r["retained_mass"] for r in rows]
    labels = [f"{r['t']:.0e}" if r["t"] != 0 else "0" for r in rows]
    axes[2].plot(x, mass, "o-")
    axes[2].axvline(3, color="crimson", linestyle="--", alpha=0.7)
    axes[2].set_xticks(x, labels, rotation=35)
    axes[2].set_title("Engineered rank event: projector jump at t=0")
    axes[2].set_xlabel("weakest singular-value scale t")
    axes[2].set_ylabel("sum of no-prior retention eigenvalues")
    axes[2].grid(True, alpha=0.25)
    fig.suptitle("Family 5 parent audit: generalized retention and rank-event semantics")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    t0 = time.perf_counter()
    cfg = f5.CONFIG
    _, chi0, _, S = f5.base_scene(cfg)
    X0 = f5.build_poses(cfg).reshape(-1)
    rng = np.random.default_rng(cfg["seeds"]["dX_unit"])
    dX = rng.standard_normal(cfg["q_pose"])
    dX /= np.linalg.norm(dX)
    eps_grid = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2]

    finite = audit_case(chi0, S, X0, dX, 1.0, eps_grid)
    no_prior = audit_case(chi0, S, X0, dX, None, eps_grid)
    finite["gate"]["pass"] = bool(all(finite["gate"].values()))
    no_prior["gate"]["pass"] = bool(all(no_prior["gate"].values()))
    cluster = cluster_audit(chi0, S, X0, dX)
    A0, B0 = f5.smooth_AB_of_flat_X(chi0, S, X0, cfg)
    event = engineered_rank_event(A0, B0)

    results = {
        "schema": "family5-parent-generalized-v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": ".venv/bin/python src/family5_parent_generalized.py",
        "scope": (
            "N=16 finite-dimensional whitened/realified smooth p=24 model; "
            "outer perturbation uses Euclidean norm in R^18; dA,dB use the "
            "same centered finite-difference builders as Family 5."
        ),
        "formula": {
            "generalized_simple": "dot(rho)=v^T(dot(Keff)-rho dot(KIS))v, v^T KIS v=1",
            "finite_prior": "E=I-B C^-1 B^T, C=B^T B+alpha I; differentiate E and A^T E A",
            "no_prior": "dot(PB)=(I-PB)dot(B)Bdag+Bdag^T dot(B)^T(I-PB), constant rank only",
            "cluster": "eigenvalues of Vc^T(dot(K)-rho dot(KIS))Vc",
        },
        "finite_prior": finite,
        "no_prior": no_prior,
        "no_prior_rho_one_cluster": cluster,
        "engineered_rank_event": event,
        "overall_gate": bool(finite["gate"]["pass"] and no_prior["gate"]["pass"]),
        "cannot_establish": [
            "The checks do not certify a continuum derivative or a uniform trajectory theorem.",
            "dA and dB are centered finite-difference approximations, so a numerical floor remains.",
            "The rank event is an engineered algebraic control, not an observed physical trajectory event.",
            "At repeated eigenvalues only the compressed cluster derivative is invariant; individual eigenvector derivatives are not claimed.",
        ],
        "environment": {
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
        },
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    fig_path = figures_dir / "family5_parent_generalized.png"
    results["figure_sha256"] = make_figure(results, fig_path)
    results["runtime_seconds"] = float(time.perf_counter() - t0)
    out = results_dir / "family5_parent_generalized.json"
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")

    print(json.dumps({
        "overall_gate": results["overall_gate"],
        "finite_gate": finite["gate"],
        "no_prior_gate": no_prior["gate"],
        "cluster_size": cluster["cluster_size"],
        "cluster_derivative_norm": cluster["compressed_derivative_spectral_norm"],
        "rank_event_projector_jump": event["projector_jump_norm_at_zero"],
        "rank_event_relative_information_jump": event["KSLAM_jump_relative_to_KIS"],
        "runtime_seconds": results["runtime_seconds"],
        "results": str(out),
        "figure": str(fig_path),
    }, indent=2))


if __name__ == "__main__":
    main()
