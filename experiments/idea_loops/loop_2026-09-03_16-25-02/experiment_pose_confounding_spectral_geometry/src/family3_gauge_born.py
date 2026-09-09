"""Family 3: gauge generators and Born empty-background degeneration.

Run (from the experiment root):
    .venv/bin/python src/family3_gauge_born.py

Scenario (identical to Families 1-2): N = 16, k_b = 2*pi, T = 6 poses on the
90-degree arc radius 1.6, n_rx = 4 body-fixed receivers, tx at body origin,
and chi0 = two Gaussian blobs.  The whitened/realified data space is
m = 2*T*n_rx = 48 and q = 3*T = 18.

Families 1-2 are NOT rerun.  This file consumes the validated
`helmholtz.build_AB` / `whiten_realify` / `forward_measurements` functions and
reuses the Family-2 smooth RBF basis and machine rank-tolerance conventions.

Contents:
  * smooth- and pixel-basis SE(2) gauge checks (Tx, Ty, rotation) for the
    two-blob scene and smooth-basis checks for the symmetric scene;
  * smooth-basis gauge refinement at N in {16, 32, 40};
  * retention spectra (R_op = Q_A^T (I - Z Z^T) Q_A) for the two-blob and
    symmetric scenes;
  * anchored / known-background / fixed-boundary variants (smooth basis);
  * Born empty-background degeneration: B0 = 0, K_SLAM0 = K_IS0, and the
    bilinear operator-uncertainty dominance check with h logspace(-4,-1,9).

Self-cell marker: this file consumes hh.SELF_CELL_FORMULA / VERSION (corrected
v2; the corrected harness is already used by the accepted Family 2 run).
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
    "arc_phi_deg": [-45.0, 45.0],
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
    "smooth_basis": {
        "p": 24,
        "x_centers_n": 4,
        "y_centers_n": 6,
        "x_span": [-0.3, 0.3],
        "y_span": [-0.3, 0.3],
        "sigma_b": 0.16,
    },
    "symmetric_scene": {
        "amp": 0.5,
        "sigma": 0.12,
    },
    "m": 48,          # 2*T*n_rx after realification
    "n_pixel": 256,   # N^2
    "q_pose": 18,     # 3*T
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "gauge_tol_5e-2": 5.0e-2,
    "gauge_refinement_Ns": [32, 40],
    "anchor_masks": {
        "base": "all True (unanchored)",
        "anchor": "free = NOT (x > 0.30 AND y > 0.30)  [fix upper-right corner]",
        "known_bg": "free = (|r| <= 0.35)  [fix outside disk]",
        "boundary": "free = (|x| <= 0.4375 AND |y| <= 0.4375)  [fix outer ring]",
    },
    "born": {
        "dchi_bump_center": [0.1, 0.05],
        "dchi_bump_sigma": 0.1,
        "dX_seed": 31415,
        "h_log_min": -4.0,
        "h_log_max": -1.0,
        "h_n": 9,
        "h_summary": 1.0e-3,
        "B0_gate": 1.0e-12,
        "KSLAM0_gate": 1.0e-12,
        "rel_err_1e-3_gate": 1.0e-3,
        "slope_gate_min": 1.8,
        "slope_gate_max": 2.2,
        "anchor_ratio_gate": 10.0,
    },
    "randomness_note": (
        "only RNG use is the fixed-seed 31415 pose direction in the Born "
        "bilinear check; all gauge/algebra blocks are deterministic"
    ),
}

_EPS = np.finfo(float).eps
_J90 = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=float)


# ---------------------------------------------------------------------------
# Geometry / scene helpers (copied conventions from Families 1-2)
# ---------------------------------------------------------------------------

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


def make_chi0_sym(points: np.ndarray, cfg: dict) -> np.ndarray:
    ss = cfg["symmetric_scene"]
    return (ss["amp"] * np.exp(
        -np.sum(points**2, axis=1) / (2.0 * ss["sigma"] ** 2)
    )).astype(float)


def build_smooth_basis(points: np.ndarray, sb_cfg: dict) -> np.ndarray:
    """Unit-2-norm Gaussian RBF columns S (n_points x p), Family-2 pattern."""
    xs = np.linspace(sb_cfg["x_span"][0], sb_cfg["x_span"][1], sb_cfg["x_centers_n"])
    ys = np.linspace(sb_cfg["y_span"][0], sb_cfg["y_span"][1], sb_cfg["y_centers_n"])
    centres = np.array([[x, y] for x in xs for y in ys], dtype=float)
    sigma_b = float(sb_cfg["sigma_b"])
    S = np.exp(
        -np.sum((points[:, None, :] - centres[None, :, :]) ** 2, axis=2)
        / (2.0 * sigma_b**2)
    )
    S = S / np.linalg.norm(S, axis=0)[None, :]
    return S


def grad_chi0(points: np.ndarray, cfg: dict) -> np.ndarray:
    """Analytic gradient of the two-blob chi0, shape (N^2, 2)."""
    g = np.zeros((points.shape[0], 2), dtype=float)
    blobs = cfg["chi_blobs"]
    for amp, c, sig in [
        (blobs["amp1"], np.asarray(blobs["c1"]), blobs["sigma1"]),
        (blobs["amp2"], np.asarray(blobs["c2"]), blobs["sigma2"]),
    ]:
        w = amp * np.exp(
            -np.sum((points - c) ** 2, axis=1) / (2.0 * sig**2)
        )
        g += w[:, None] * (-(points - c) / sig**2)
    return g


def grad_chi_sym(points: np.ndarray, cfg: dict) -> np.ndarray:
    """Analytic gradient of chi0_sym = 0.5 exp(-|r|^2/(2*0.12^2))."""
    ss = cfg["symmetric_scene"]
    sigma2 = ss["sigma"] ** 2
    return (
        -ss["amp"]
        * np.exp(-np.sum(points**2, axis=1) / (2.0 * sigma2))[:, None]
        * (points / sigma2)
    )


# ---------------------------------------------------------------------------
# Rank / projector helpers (Family-2 conventions copied verbatim)
# ---------------------------------------------------------------------------

def rank_svd(M: np.ndarray) -> tuple[int, np.ndarray, float]:
    sv = np.linalg.svd(M, compute_uv=False)
    tol = max(M.shape) * _EPS * sv[0]
    return int(np.sum(sv > tol)), sv, float(tol)


def thin_decomposition(A: np.ndarray) -> dict:
    r, sv, tol = rank_svd(A)
    u, s, vh = np.linalg.svd(A, full_matrices=False)
    return {
        "rank": r,
        "rank_tol": tol,
        "singular_values": sv,
        "Q": u[:, :r],
        "s": s[:r],
        "V": vh[:r].T,
    }


def range_basis(B: np.ndarray) -> tuple[int, np.ndarray, np.ndarray, float]:
    r, sv, tol = rank_svd(B)
    u, _, _ = np.linalg.svd(B, full_matrices=False)
    return r, u[:, :r], sv, tol


def retention_spectrum(Q: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Descending eigenvalues of Q^T W Q (W symmetric), Family-2 form."""
    M = Q.T @ W @ Q
    M = 0.5 * (M + M.T)
    return np.sort(np.linalg.eigvalsh(M))[::-1]


def spectrum_report(A: np.ndarray, B: np.ndarray) -> dict:
    """R_op = Q_A^T (I - Z Z^T) Q_A eigenvalues plus rank info."""
    da = thin_decomposition(A)
    QA = da["Q"]
    rB, Z, svB, tolB = range_basis(B)
    P = np.eye(A.shape[0], dtype=float) - Z @ Z.T
    rho = retention_spectrum(QA, P)
    return {
        "rank_A": da["rank"],
        "rank_A_tol": da["rank_tol"],
        "rank_B": rB,
        "rank_B_tol": tolB,
        "sigma_A_min": float(da["singular_values"][-1]),
        "sigma_B_min": float(svB[-1]) if len(svB) else None,
        "rho_desc": [float(x) for x in rho],
        "rho_min": float(rho.min()) if len(rho) else None,
        "rho_max": float(rho.max()) if len(rho) else None,
        "count_rho_lt_1_minus_1e-8": int(np.sum(rho < 1.0 - 1e-8)),
    }


# ---------------------------------------------------------------------------
# SE(2) generators
# ---------------------------------------------------------------------------

def gauge_dchi(grad: np.ndarray, points: np.ndarray) -> dict[str, np.ndarray]:
    """Pixel vectors dchi_g = - (generator field) . grad chi(r)."""
    jr = (_J90 @ points.T).T            # (N^2, 2), omega = +1
    return {
        "tx": -grad[:, 0],
        "ty": -grad[:, 1],
        "rot": -np.einsum("ni,ni->n", jr, grad),
    }


def gauge_dX(poses: np.ndarray) -> dict[str, np.ndarray]:
    T = poses.shape[0]
    tx = np.tile(np.array([1.0, 0.0, 0.0]), T)
    ty = np.tile(np.array([0.0, 1.0, 0.0]), T)
    rot = np.concatenate(
        [np.r_[(_J90 @ pose[:2]), 1.0] for pose in poses]
    )
    return {"tx": tx, "ty": ty, "rot": rot}


# ---------------------------------------------------------------------------
# Scene linearisations
# ---------------------------------------------------------------------------

def analysis_setup(
    chi: np.ndarray, cfg: dict
) -> dict:
    """A_pix_R, B_R, S, A_smooth_R for one scene at cfg['N']."""
    poses = build_poses(cfg)
    points, h = hh.make_grid(cfg["N"])
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)
    A, B, F, _ = hh.build_AB(chi, poses, rx, tx, cfg["N"], cfg["k_b"])
    A_pix_R, B_R = hh.whiten_realify(A, B, None)
    S = build_smooth_basis(points, cfg["smooth_basis"])
    A_smooth_R = A_pix_R @ S
    return {
        "poses": poses,
        "points": points,
        "grid_h": h,
        "A_pix_R": A_pix_R,
        "B_R": B_R,
        "S": S,
        "A_smooth_R": A_smooth_R,
        "F_complex": F,
    }


def fit_coefficients(S: np.ndarray | None, dchi: np.ndarray):
    """c_g = lstsq(S, dchi) for smooth basis; exact pixel vector for pixel."""
    if S is None:
        return dchi.copy(), 0.0
    c, _, _, _ = np.linalg.lstsq(S, dchi, rcond=None)
    resid = float(
        np.linalg.norm(S @ c - dchi) / max(np.linalg.norm(dchi), _EPS)
    )
    return c, resid


def gauge_rows(
    A: np.ndarray,
    B: np.ndarray,
    S: np.ndarray | None,
    dchis: dict[str, np.ndarray],
    dXs: dict[str, np.ndarray],
    include_kernel_distance: bool = True,
    generators: list[str] | None = None,
) -> dict:
    rows = []
    P = None
    KSL = None
    if generators is None:
        generators = ["tx", "ty", "rot"]
    if include_kernel_distance:
        _, Z, _, _ = range_basis(B)
        P = np.eye(A.shape[0], dtype=float) - Z @ Z.T
        KSL = A.T @ (P @ A)
    for name in generators:
        dchi = dchis[name]
        dX = dXs[name]
        c, rep = fit_coefficients(S, dchi)
        Ac = A @ c
        Bd = B @ dX
        num = float(np.linalg.norm(Ac + Bd, ord=2))
        den = float(np.linalg.norm(Ac, ord=2) + np.linalg.norm(Bd, ord=2))
        gauge = num / max(den, _EPS)
        row = {
            "generator": name,
            "dchi_l2": float(np.linalg.norm(dchi, ord=2)),
            "c_l2": float(np.linalg.norm(c, ord=2)),
            "representation_residual": rep,
            "representation_guard": _EPS,
            "map_norm_A_c": float(np.linalg.norm(Ac, ord=2)),
            "pose_norm_B_dX": float(np.linalg.norm(Bd, ord=2)),
            "gauge_numerator": num,
            "gauge_residual": gauge,
        }
        if include_kernel_distance:
            row["kernel_distance"] = float(
                np.linalg.norm(KSL @ c, ord=2)
                / max(np.linalg.norm(c, ord=2), _EPS)
            )
        rows.append(row)
    return {"rows": rows}


# ---------------------------------------------------------------------------
# Convergence-slope helper (window just above the roundoff minimum)
# ---------------------------------------------------------------------------

def fit_clean_slope(hs: np.ndarray, errs: np.ndarray, window_len: int = 4) -> dict:
    hs = np.asarray(hs, dtype=float)
    errs = np.asarray(errs, dtype=float)
    n = len(hs)
    i_min = int(np.argmin(errs))
    if n < window_len:  # pragma: no cover - config has 9 samples
        start, end = 0, n - 1
    else:
        start = min(i_min + 1, n - window_len)
        if start < 0:
            start = max(0, n - window_len)
        end = min(start + window_len - 1, n - 1)
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
        "window_indices": list(range(start, end + 1)),
        "h_at_min_error": float(hs[i_min]),
        "min_error": float(errs[i_min]),
        "note": (
            "window starts one index after the global error minimum (or the "
            "last four samples if the minimum is at the end), matching the "
            "roundoff-conscious convention used in Family 1"
        ),
    }


# ---------------------------------------------------------------------------
# Born empty-background checks
# ---------------------------------------------------------------------------

def realify_complex_vector(z: np.ndarray) -> np.ndarray:
    return np.sqrt(2.0) * np.concatenate([np.real(z), np.imag(z)])


def run_born_checks(
    cfg: dict,
    poses: np.ndarray,
    rx: np.ndarray,
    tx: np.ndarray,
    points: np.ndarray,
) -> dict:
    N = cfg["N"]
    T = cfg["T"]
    X0 = poses.reshape(-1)
    zeros = np.zeros(N * N, dtype=float)

    t0 = time.perf_counter()
    A0, B0, _, _ = hh.build_AB(zeros, poses, rx, tx, N, cfg["k_b"])
    A0_R, B0_R = hh.whiten_realify(A0, B0, None)
    build0_seconds = time.perf_counter() - t0

    rel_B0 = float(
        np.linalg.norm(B0_R, ord="fro") / max(1.0, np.linalg.norm(A0_R, ord="fro"))
    )
    K_IS0 = A0_R.T @ A0_R
    rB0, Z0, _, _ = range_basis(B0_R)
    P0 = np.eye(A0_R.shape[0], dtype=float) - Z0 @ Z0.T
    K_SLAM0 = A0_R.T @ (P0 @ A0_R)
    rel_K0 = float(
        np.linalg.norm(K_SLAM0 - K_IS0, ord="fro")
        / max(np.linalg.norm(K_IS0, ord="fro"), _EPS)
    )

    # Fixed smooth unit-norm pixel direction and fixed pose direction.
    center = np.asarray(cfg["born"]["dchi_bump_center"], dtype=float)
    sigma_bump = cfg["born"]["dchi_bump_sigma"]
    dchi = np.exp(
        -np.sum((points - center) ** 2, axis=1) / (2.0 * sigma_bump**2)
    )
    dchi = dchi / np.linalg.norm(dchi)
    rng = np.random.default_rng(cfg["born"]["dX_seed"])
    dX = rng.standard_normal(3 * T)
    dX = dX / np.linalg.norm(dX)

    def A0_real_at(Xflat: np.ndarray) -> np.ndarray:
        A_x, B_x, _, _ = hh.build_AB(
            zeros, Xflat.reshape(T, 3), rx, tx, N, cfg["k_b"]
        )
        return hh.whiten_realify(A_x, B_x, None)[0]

    def F_real_at(chi: np.ndarray, Xflat: np.ndarray) -> np.ndarray:
        Fz = hh.forward_measurements(
            chi, Xflat.reshape(T, 3), rx, tx, N, cfg["k_b"]
        )
        return realify_complex_vector(Fz)

    linear = A0_R @ dchi
    F0_R = F_real_at(dchi, X0)
    R0 = F0_R - linear

    steps = np.logspace(
        cfg["born"]["h_log_min"],
        cfg["born"]["h_log_max"],
        cfg["born"]["h_n"],
    )
    rows = []
    t_fd0 = time.perf_counter()
    for step in steps:
        A_p = A0_real_at(X0 + step * dX)
        A_m = A0_real_at(X0 - step * dX)
        BIL = ((A_p - A_m) / (2.0 * step)) @ dchi
        Fp = F_real_at(dchi, X0 + step * dX)
        Fm = F_real_at(dchi, X0 - step * dX)
        FD = (Fp - Fm) / (2.0 * step)
        rel_err = float(
            np.linalg.norm(FD - BIL, ord=2) / max(np.linalg.norm(BIL, ord=2), _EPS)
        )
        R_h = Fp - linear - step * BIL
        rows.append(
            {
                "h": float(step),
                "BIL_norm": float(np.linalg.norm(BIL, ord=2)),
                "FD_norm": float(np.linalg.norm(FD, ord=2)),
                "rel_err": rel_err,
                "remainder_norm_R_h": float(np.linalg.norm(R_h, ord=2)),
                "R_h_minus_R0_norm": float(
                    np.linalg.norm(R_h - R0, ord=2)
                ),
            }
        )
    fd_seconds = time.perf_counter() - t_fd0

    # Exact h = 1e-3 summary point (not an element of the 9-point log grid).
    hq = cfg["born"]["h_summary"]
    A_p = A0_real_at(X0 + hq * dX)
    A_m = A0_real_at(X0 - hq * dX)
    BIL_q = ((A_p - A_m) / (2.0 * hq)) @ dchi
    Fp_q = F_real_at(dchi, X0 + hq * dX)
    Fm_q = F_real_at(dchi, X0 - hq * dX)
    FD_q = (Fp_q - Fm_q) / (2.0 * hq)
    R_hq = Fp_q - linear - hq * BIL_q
    rel_err_q = float(
        np.linalg.norm(FD_q - BIL_q, ord=2)
        / max(np.linalg.norm(BIL_q, ord=2), _EPS)
    )
    summary_q = {
        "h": hq,
        "bilinear_norm": float(np.linalg.norm(BIL_q, ord=2)),
        "linear_norm_A0_dchi": float(np.linalg.norm(linear, ord=2)),
        "bilinear_over_linear": float(
            np.linalg.norm(BIL_q, ord=2) / max(np.linalg.norm(linear, ord=2), _EPS)
        ),
        "rel_err": rel_err_q,
        "remainder_norm_R_h": float(np.linalg.norm(R_hq, ord=2)),
        "R0_norm": float(np.linalg.norm(R0, ord=2)),
        "R_h_minus_R0_norm": float(np.linalg.norm(R_hq - R0, ord=2)),
    }

    hs = np.asarray([r["h"] for r in rows])
    rel_errs = np.asarray([r["rel_err"] for r in rows])
    rem_diffs = np.asarray([r["R_h_minus_R0_norm"] for r in rows])
    fit_rel_clean = fit_clean_slope(hs, rel_errs)
    fit_rem_clean = fit_clean_slope(hs, rem_diffs)

    # First-four slope for audit alongside the clean-window slope.
    def first4_slope(ys: np.ndarray) -> dict:
        xs = np.log10(hs[:4])
        yy = np.log10(ys[:4])
        s, ic = np.polyfit(xs, yy, 1)
        return {"slope": float(s), "window_indices": [0, 1, 2, 3]}

    return {
        "build0_seconds": build0_seconds,
        "fd_seconds": fd_seconds,
        "A0_norm_F": float(np.linalg.norm(A0_R, ord="fro")),
        "B0_norm_F": float(np.linalg.norm(B0_R, ord="fro")),
        "B0_rel_gate_metric": rel_B0,
        "rank_B0": rB0,
        "K_SLAM0_minus_K_IS0_rel_F": rel_K0,
        "dchi": {
            "center": center.tolist(),
            "sigma": sigma_bump,
            "l2": float(np.linalg.norm(dchi)),
            "max_abs": float(np.max(np.abs(dchi))),
        },
        "dX": {
            "seed": cfg["born"]["dX_seed"],
            "l2": float(np.linalg.norm(dX)),
            "first_entries": [float(x) for x in dX[:6]],
        },
        "linear_norm_A0_dchi": float(np.linalg.norm(linear)),
        "R0_norm": float(np.linalg.norm(R0)),
        "rows": rows,
        "summary_at_1e-3": summary_q,
        "fit_rel_err_clean_window": fit_rel_clean,
        "fit_rel_err_first4": first4_slope(rel_errs),
        "fit_RminusR0_clean_window": fit_rem_clean,
        "fit_RminusR0_first4": first4_slope(rem_diffs),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    t0 = datetime.now(timezone.utc)
    t_start = time.perf_counter()

    poses = build_poses(cfg)
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)
    points16, h16 = hh.make_grid(cfg["N"])
    chi0 = make_chi0(points16, cfg)
    chi0_sym = make_chi0_sym(points16, cfg)

    print(f"[family3] building N={cfg['N']} scenes ...")
    t_two = time.perf_counter()
    two = analysis_setup(chi0, cfg)
    build_two_seconds = time.perf_counter() - t_two
    t_sym = time.perf_counter()
    sym = analysis_setup(chi0_sym, cfg)
    build_sym_seconds = time.perf_counter() - t_sym

    grad2 = grad_chi0(points16, cfg)
    dchis_two = gauge_dchi(grad2, points16)
    dXs = gauge_dX(poses)
    g2_tx = gauge_rows(
        two["A_smooth_R"], two["B_R"], two["S"], dchis_two, dXs
    )
    g2_pix = gauge_rows(
        two["A_pix_R"], two["B_R"], None, dchis_two, dXs
    )
    print("[family3] smooth two-blob gauge residuals:",
          {r["generator"]: round(r["gauge_residual"], 6) for r in g2_tx["rows"]})
    print("[family3] pixel two-blob gauge residuals:",
          {r["generator"]: round(r["gauge_residual"], 8) for r in g2_pix["rows"]})

    # Symmetric scene (smooth basis)
    grads = grad_chi_sym(points16, cfg)
    dchis_sym = gauge_dchi(grads, points16)
    g_sym = gauge_rows(
        sym["A_smooth_R"], sym["B_R"], sym["S"], dchis_sym, dXs,
        include_kernel_distance=False,
    )
    print("[family3] symmetric smooth gauge residuals:",
          {r["generator"]: round(r["gauge_residual"], 6) for r in g_sym["rows"]})

    spectra_two = spectrum_report(two["A_smooth_R"], two["B_R"])
    spectra_sym = spectrum_report(sym["A_smooth_R"], sym["B_R"])
    print(f"[family3] rho_min smooth two-blob {spectra_two['rho_min']:.3e}; "
          f"symmetric {spectra_sym['rho_min']:.3e}")

    # Smooth-basis gauge refinement
    refinement = []
    for N_res in [cfg["N"]] + list(cfg["gauge_refinement_Ns"]):
        cfg_r = dict(cfg)
        cfg_r["N"] = N_res
        t_r = time.perf_counter()
        if N_res == cfg["N"]:
            setup_r = two
            pts_r = points16
            dchi_r = dchis_two
        else:
            pts_r, _ = hh.make_grid(N_res)
            chi_r = make_chi0(pts_r, cfg)
            setup_r = analysis_setup(chi_r, cfg_r)
            dchi_r = gauge_dchi(grad_chi0(pts_r, cfg), pts_r)
        subset = {k: dchi_r[k] for k in ["tx", "rot"]}
        dX_sub = {k: dXs[k] for k in ["tx", "rot"]}
        rr = gauge_rows(
            setup_r["A_smooth_R"], setup_r["B_R"], setup_r["S"],
            subset, dX_sub,
            generators=["tx", "rot"],
        )
        refinement.append(
            {
                "N": N_res,
                "grid_h": float(setup_r["grid_h"]),
                "rows": rr["rows"],
                "seconds": time.perf_counter() - t_r,
            }
        )
    print("[family3] refinement:", [
        (r["N"], {x["generator"]: round(x["gauge_residual"], 6)
                 for x in r["rows"]}) for r in refinement
    ])

    # Anchor / known-background / fixed-boundary variants (translation x)
    tx_dchi = dchis_two["tx"]
    tx_dX = dXs["tx"]
    x, y = points16[:, 0], points16[:, 1]
    masks = {
        "base": np.ones(points16.shape[0], dtype=bool),
        "anchor": ~((x > 0.30) & (y > 0.30)),
        "known_bg": np.linalg.norm(points16, axis=1) <= 0.35,
        "boundary": (np.abs(x) <= 0.4375) & (np.abs(y) <= 0.4375),
    }
    anchor_rows = []
    base_min = None
    for name, mask in masks.items():
        A_red = two["A_pix_R"][:, mask] @ two["S"][mask, :]
        dchi_free = tx_dchi[mask]
        c, rep = fit_coefficients(two["S"][mask, :], dchi_free)
        Ac = A_red @ c
        Bd = two["B_R"] @ tx_dX
        gauge = float(
            np.linalg.norm(Ac + Bd, ord=2)
            / max(np.linalg.norm(Ac, ord=2) + np.linalg.norm(Bd, ord=2), _EPS)
        )
        sp = spectrum_report(A_red, two["B_R"])
        if name == "base":
            base_min = sp["rho_min"]
        anchor_rows.append(
            {
                "name": name,
                "mask_definition": cfg["anchor_masks"][name],
                "n_free_pixels": int(np.sum(mask)),
                "n_fixed_pixels": int(np.sum(~mask)),
                "c_l2": float(np.linalg.norm(c)),
                "representation_residual_free": rep,
                "gauge_residual_tx": gauge,
                "spectrum": sp,
                "rho_min": sp["rho_min"],
            }
        )
    for row in anchor_rows:
        row["rho_min_ratio_vs_base"] = (
            float(row["rho_min"] / max(base_min, 0.0))
            if row["name"] != "base" else 1.0
        )
    print("[family3] anchor rho_min:", {r["name"]: r["rho_min"] for r in anchor_rows})

    # Born empty-background checks
    print("[family3] running Born empty-background checks ...")
    born = run_born_checks(cfg, poses, rx, tx, points16)
    print("[family3] born rel B0:", born["B0_rel_gate_metric"],
          " rel K:", born["K_SLAM0_minus_K_IS0_rel_F"],
          " rel_err@1e-3:", born["summary_at_1e-3"]["rel_err"])

    # ---------------- gates ------------------------------------------------
    smooth16 = {r["generator"]: r for r in g2_tx["rows"]}
    ref16 = {r["generator"]: r for r in refinement[0]["rows"]}
    ref32 = {r["generator"]: r for r in refinement[1]["rows"]}
    ref40 = {r["generator"]: r for r in refinement[2]["rows"]}
    gates = {
        "smooth_two_blob_gauge_le_5e-2_N16": bool(
            all(r["gauge_residual"] <= cfg["gauge_tol_5e-2"]
                for r in g2_tx["rows"])
        ),
        "smooth_two_blob_gauge_tol_5e-2": cfg["gauge_tol_5e-2"],
        "smooth_tx_gauge_decreases_16_to_40": bool(
            ref40["tx"]["gauge_residual"] < ref16["tx"]["gauge_residual"]
        ),
        "smooth_rot_gauge_decreases_16_to_40": bool(
            ref40["rot"]["gauge_residual"] < ref16["rot"]["gauge_residual"]
        ),
        "smooth_gauge_decreases_16_to_40_overall": bool(
            ref40["tx"]["gauge_residual"] < ref16["tx"]["gauge_residual"]
            and ref40["rot"]["gauge_residual"] < ref16["rot"]["gauge_residual"]
        ),
        "born_B0_rel_gate_lt_1e-12": bool(
            born["B0_rel_gate_metric"] < cfg["born"]["B0_gate"]
        ),
        "born_B0_gate_tol": cfg["born"]["B0_gate"],
        "born_KSLAM0_eq_KIS0_rel_lt_1e-12": bool(
            born["K_SLAM0_minus_K_IS0_rel_F"] < cfg["born"]["KSLAM0_gate"]
        ),
        "born_KSLAM0_gate_tol": cfg["born"]["KSLAM0_gate"],
        "born_rel_err_at_1e-3_lt_1e-3": bool(
            born["summary_at_1e-3"]["rel_err"]
            < cfg["born"]["rel_err_1e-3_gate"]
        ),
        "born_rel_err_1e-3_gate_tol": cfg["born"]["rel_err_1e-3_gate"],
        "born_rel_err_clean_slope_in_1.8_2.2": bool(
            cfg["born"]["slope_gate_min"]
            <= born["fit_rel_err_clean_window"]["slope"]
            <= cfg["born"]["slope_gate_max"]
        ),
        "born_slope_gate_interval": [
            cfg["born"]["slope_gate_min"],
            cfg["born"]["slope_gate_max"],
        ],
        "anchor_all_restricted_gt_10x_base": bool(
            all(
                row["rho_min_ratio_vs_base"] > cfg["born"]["anchor_ratio_gate"]
                for row in anchor_rows
                if row["name"] != "base"
            )
        ),
        "anchor_ratio_gate": cfg["born"]["anchor_ratio_gate"],
    }

    total_runtime = time.perf_counter() - t_start
    source_files = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
        "src/family3_gauge_born.py",
    ]
    results = {
        "generated_utc": t0.isoformat(),
        "runner": "src/family3_gauge_born.py",
        "command": ".venv/bin/python src/family3_gauge_born.py",
        "family": 3,
        "self_cell_formula": hh.SELF_CELL_FORMULA,
        "self_cell_formula_version": hh.SELF_CELL_FORMULA_VERSION,
        "parent_correction_gate": (
            "context/PARENT_CORRECTIONS.md / workshop rule 16; family3 "
            "consumes only the corrected helmholtz build_AB/forward_measurements"
        ),
        "runtime_seconds": total_runtime,
        "build_two_blob_seconds_N16": build_two_seconds,
        "build_symmetric_seconds_N16": build_sym_seconds,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
        },
        "environment_note": (
            "Apple Silicon CPU, no GPU/MPS/CUDA; deterministic dense "
            "numpy/scipy linear algebra"
        ),
        "source_sha256": {
            p: hashlib.sha256((_ROOT / p).read_bytes()).hexdigest()
            for p in source_files
        },
        "config": {
            **cfg,
            "k_b": float(cfg["k_b"]),
            "poses": [p.tolist() for p in poses],
            "rx_offsets": rx.tolist(),
            "tx_offset": tx.tolist(),
            "grid_h_cell_N16": h16,
            "chi0_stats": {
                "min": float(chi0.min()),
                "max": float(chi0.max()),
                "mean": float(chi0.mean()),
                "l2": float(np.linalg.norm(chi0)),
            },
            "chi0_sym_stats": {
                "min": float(chi0_sym.min()),
                "max": float(chi0_sym.max()),
                "mean": float(chi0_sym.mean()),
                "l2": float(np.linalg.norm(chi0_sym)),
            },
        },
        "dimensions_common": {
            "m_real_data": int(two["A_pix_R"].shape[0]),
            "q_pose": int(two["B_R"].shape[1]),
            "S_N2_N16": int(two["A_pix_R"].shape[1]),
            "p_smooth": int(two["S"].shape[1]),
            "T": cfg["T"],
            "n_rx": cfg["n_rx"],
        },
        "gauge_two_blob_smooth": g2_tx,
        "gauge_two_blob_pixel": g2_pix,
        "gauge_symmetric_smooth": g_sym,
        "smooth_spectrum_two_blob": spectra_two,
        "smooth_spectrum_symmetric": spectra_sym,
        "gauge_refinement_smooth": refinement,
        "anchor_variants": {
            "rows": anchor_rows,
            "generator": "translation-x",
            "dchi_l2_full": float(np.linalg.norm(tx_dchi)),
        },
        "born_empty_background": born,
        "gates": gates,
        "gauge_note": (
            "The prompt's gauge definition is implemented literally: dchi_g is "
            "the analytic pixel derivative and c_g = lstsq(S, dchi_g) for the "
            "smooth basis (pixel basis: c_g = dchi_g).  Exact cancellation "
            "would require A c_g = -B dX_g; for the smooth RBF basis this is "
            "limited by how well span(S) contains the derivative fields, and "
            "for the pixel basis by discrete forward-model equivariance."
        ),
        "born_note": (
            "hh.forward_measurements is the full-wave forward map.  For the "
            "fixed unit-L2 dchi the finite-difference derivative D_X F(dchi) "
            "contains full-wave terms beyond the Born bilinear term "
            "(D_X A(0)) dchi, so rel_err has a nonzero floor controlled by "
            "the scene amplitude and does not follow O(h^2); raw values and "
            "fits are recorded without hiding the floor."
        ),
        "linear_algebra_scope_note": (
            "Finite-dimensional, discrete, whitened/realified checks on the "
            "validated Family-1 Jacobians.  No continuum-limit or exact-global-"
            "SE(2) claim is made; a bounded grid has no exact global SE(2) "
            "equivariance."
        ),
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)
    results_path = results_dir / "family3_gauge_born.json"
    results_path.write_text(json.dumps(results, indent=2, default=lambda o: o.tolist()) + "\n")

    # --------------------------- figures -----------------------------------
    # 1. gauge residual bars
    labels = ["Tx", "Ty", "Rot"]
    groups = [
        ("smooth two-blob", [smooth16[g]["gauge_residual"] for g in ["tx", "ty", "rot"]]),
        ("pixel two-blob", [g2_pix["rows"][i]["gauge_residual"] for i in range(3)]),
        ("smooth symmetric", [g_sym["rows"][i]["gauge_residual"] for i in range(3)]),
    ]
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    xpos = np.arange(3)
    width = 0.26
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    for gi, ((name, vals), col) in enumerate(zip(groups, colors)):
        ax.bar(xpos + (gi - 1) * width, vals, width, label=name, color=col,
               log=True)
        for xx, vv in zip(xpos + (gi - 1) * width, vals):
            ax.text(xx, vv * 1.12, f"{vv:.2e}", ha="center", va="bottom",
                    fontsize=7.5, rotation=90)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("gauge residual (log)")
    ax.set_ylim(top=max(v for _, vs in groups for v in vs) * 8.0)
    ax.set_title("Family 3: SE(2) gauge residual\n"
                 "$||A\\,c_g+B\\,dX_g||/(||A\\,c_g||+||B\\,dX_g||)$")
    ax.grid(True, which="both", axis="y", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig_path1 = figures_dir / "family3_gauge_residuals.png"
    fig.savefig(fig_path1, dpi=200)
    plt.close(fig)

    # 2. gauge refinement vs N
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 5.0), sharey=False)
    Ns = [r["N"] for r in refinement]
    for gi, gen in enumerate(["tx", "rot"]):
        gr = [next(x for x in r["rows"] if x["generator"] == gen)
              for r in refinement]
        ax = axs[gi]
        ax.semilogy(Ns, [x["gauge_residual"] for x in gr], "o-", color="#1f77b4",
                    label="gauge residual")
        ax.semilogy(Ns, [x["representation_residual"] for x in gr], "s--",
                    color="#d62728", label="representation residual")
        ax.set_xlabel("N")
        ax.set_title(f"smooth translation-{gen} refinement")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
    axs[0].set_ylabel("residual (log)")
    fig.suptitle("Family 3: smooth-basis gauge residual vs N (p=24 fixed)")
    fig.tight_layout()
    fig_path2 = figures_dir / "family3_gauge_refinement.png"
    fig.savefig(fig_path2, dpi=200)
    plt.close(fig)

    # 3. anchor retention spectra (ascending order, log y)
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    colors3 = {"base": "k", "anchor": "#1f77b4", "known_bg": "#d62728",
               "boundary": "#2ca02c"}
    markers3 = {"base": "o", "anchor": "s", "known_bg": "^", "boundary": "v"}
    for row in anchor_rows:
        rho = np.sort(np.asarray(row["spectrum"]["rho_desc"]))  # ascending
        rho = rho[rho > 0.0]
        ax.semilogy(np.arange(1, len(rho) + 1), rho, markers3[row["name"]] + "-",
                    ms=4.0, lw=1.0, color=colors3[row["name"]],
                    label=(f"{row['name']} (rho_min={row['rho_min']:.2e}, "
                           f"free={row['n_free_pixels']})"))
    ax.axhline(1.0, color="grey", ls=":", lw=0.8)
    ax.set_xlabel("sorted index (ascending retention)")
    ax.set_ylabel("retention eigenvalue rho (log)")
    ax.set_title("Family 3: retention spectra, unanchored vs restricted supports\n"
                 "(smooth basis, translation-x scene generator)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig_path3 = figures_dir / "family3_anchor_rho.png"
    fig.savefig(fig_path3, dpi=200)
    plt.close(fig)

    # 4. Born bilinear check
    hs_arr = np.asarray([r["h"] for r in born["rows"]])
    rel_arr = np.asarray([r["rel_err"] for r in born["rows"]])
    rem_arr = np.asarray([r["R_h_minus_R0_norm"] for r in born["rows"]])
    hq = born["summary_at_1e-3"]["h"]
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 5.0))
    axs[0].loglog(hs_arr, rel_arr, "o-", color="#1f77b4",
                  label="rel_err(h)")
    ref_y = born["summary_at_1e-3"]["rel_err"] * (hq / hs_arr) ** 2
    axs[0].loglog(hs_arr, ref_y, "--", color="grey",
                  label=r"$O(h^2)$ through h=1e-3")
    axs[0].axhline(born["summary_at_1e-3"]["rel_err"], color="crimson",
                   ls=":", lw=1.0, label="observed floor")
    axs[0].set_xlabel("h")
    axs[0].set_ylabel("||FD_h - BIL_h|| / ||BIL_h||")
    axs[0].set_title("Born empty background: bilinear FD match\n"
                     "(full-wave D_XF floor is not O(h^2))")
    axs[0].grid(True, which="both", alpha=0.3)
    axs[0].legend(fontsize=8)

    axs[1].loglog(hs_arr, rem_arr, "o-", color="#d62728",
                  label=r"$||R(h)-R_0||$")
    ref2 = born["summary_at_1e-3"]["R_h_minus_R0_norm"] * (hq / hs_arr) ** 2
    axs[1].loglog(hs_arr, ref2, "--", color="grey",
                  label=r"$O(h^2)$ through h=1e-3")
    axs[1].set_xlabel("h")
    axs[1].set_ylabel(r"$||R(h)-R_0||$")
    axs[1].set_title("Born empty background: remainder vs h\n"
                     "(observed linear-in-h full-wave term)")
    axs[1].grid(True, which="both", alpha=0.3)
    axs[1].legend(fontsize=8)
    fig.suptitle("Family 3: Born bilinear operator-uncertainty dominance")
    fig.tight_layout()
    fig_path4 = figures_dir / "family3_born_empty_background.png"
    fig.savefig(fig_path4, dpi=200)
    plt.close(fig)

    # --------------------------- console summary ---------------------------
    print("\n===== FAMILY 3 GAUGE / BORN SUMMARY =====")
    print("gauge smooth two-blob:")
    for r in g2_tx["rows"]:
        print(f"  {r['generator']:4s} rep={r['representation_residual']:.4e} "
              f"gauge={r['gauge_residual']:.4e} "
              f"Kdist={r['kernel_distance']:.4e}")
    print("gauge pixel two-blob:")
    for r in g2_pix["rows"]:
        print(f"  {r['generator']:4s} gauge={r['gauge_residual']:.4e} "
              f"Kdist={r['kernel_distance']:.4e}")
    print("gauge symmetric smooth:")
    for r in g_sym["rows"]:
        print(f"  {r['generator']:4s} dchi={r['dchi_l2']:.4e} "
              f"rep={r['representation_residual']:.4e} "
              f"gauge={r['gauge_residual']:.4e}")
    print("anchor rho_min / ratio:")
    for r in anchor_rows:
        print(f"  {r['name']:9s} rho_min={r['rho_min']:.4e} "
              f"ratio={r['rho_min_ratio_vs_base']:.4f} "
              f"gauge_tx={r['gauge_residual_tx']:.4e}")
    print("born: relB0", born["B0_rel_gate_metric"],
          " relK", born["K_SLAM0_minus_K_IS0_rel_F"])
    print("born bilinear:", {
        "BIL/linear": born["summary_at_1e-3"]["bilinear_over_linear"],
        "rel_err": born["summary_at_1e-3"]["rel_err"],
        "R-R0": born["summary_at_1e-3"]["R_h_minus_R0_norm"],
        "slope_rel": born["fit_rel_err_clean_window"]["slope"],
        "slope_R": born["fit_RminusR0_clean_window"]["slope"],
    })
    print("gates:", {k: v for k, v in gates.items() if isinstance(v, bool)})
    print(f"total runtime {total_runtime:.2f}s")
    print("results ->", results_path)
    print("figure  ->", fig_path1)
    print("figure  ->", fig_path2)
    print("figure  ->", fig_path3)
    print("figure  ->", fig_path4)


if __name__ == "__main__":
    main()
