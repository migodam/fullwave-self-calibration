"""Family 3b: correctly-conditioned refinements of the Family 3 claims.

Run (from the experiment root):
    .venv/bin/python src/family3b_refinements.py

Family 3 recorded three checks whose specified form did not match the theory's
validity conditions:
  * the smooth RBF basis (p=24) cannot represent the derivative/generator
    fields it is asked to gauge-test (representation-limited smooth checks);
  * the Born bilinear/remainder check used an O(1) scene amplitude, so the
    full-wave finite-difference error was dominated by the chi^2 remainder
    rather than by the O(h^2) discretization whose slope was being claimed;
  * the anchored/known-background/boundary claim used the global rho_min of a
    full retention spectrum instead of the generator-specific retention.

This file does NOT modify or delete results/family3_gauge_born.json or
notes/family3_report.md.  It keeps those as the raw Family 3 record and writes
separate refined outputs:
  results/family3b_refinements.json
  figures/family3b_gauge_pixel_refinement.png
  figures/family3b_anchor_generator_retention.png
  figures/family3b_born_bilinear_refined.png
  notes/family3b_report.md

Checks
  A  pixel-basis gauge residual and kernel distance vs N in {16, 32, 40};
  B  p=27 augmented smooth basis (24 RBF columns + the three normalised
     generator fields) at N=16, side by side with the original p=24 basis;
  C  generator-specific retention rho_g for the translation-x generator with
     pixel basis and base/corner-anchor/known-disk/outer-ring masks;
  D  corrected Born empty-background bilinear check with ||dchi_small||=1e-3.

Raw Family 3 gates that failed remain failed in the record; this script only
supplies the correctly conditioned tests and states which refined checks pass.
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
import family1_pilot as f1  # noqa: E402  (build_poses, make_chi0, cfg)
import family2_algebraic_spine as f2  # noqa: E402  (basis, rank tolerances)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration (scenario verbatim from Families 1-3)
# ---------------------------------------------------------------------------

CONFIG = {
    "N": 16,
    "k_b": 2.0 * np.pi,
    "T": 6,
    "n_rx": 4,
    "q": 1.0,
    "rx_offsets": f1.CONFIG["rx_offsets"],
    "tx_offset": f1.CONFIG["tx_offset"],
    "arc_radius": f1.CONFIG["arc_radius"],
    "arc_phi_deg": f1.CONFIG["arc_phi_deg"],
    "pose_theta_convention": f1.CONFIG["pose_theta_convention"],
    "chi_blobs": f1.CONFIG["chi_blobs"],
    "smooth_basis": f2.CONFIG["smooth_basis"],
    "symmetric_scene": {
        "amp": 0.5,
        "sigma": 0.12,
    },
    "m": 48,
    "n_pixel": 256,
    "q_pose": 18,
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "pixel_gauge_Ns": [16, 32, 40],
    "gauge_tol_refined": 1.0e-3,
    "augmented_representation_tol": 1.0e-10,
    "augmented_gauge_tol": 1.0e-3,
    "anchor_masks": {
        "base": "all True (unanchored)",
        "anchor": "free = NOT (x > 0.30 AND y > 0.30)  [fix upper-right corner]",
        "known_bg": "free = (|r| <= 0.35)  [fix outside disk]",
        "boundary": "free = (|x| <= 0.4375 AND |y| <= 0.4375)  [fix outer ring]",
    },
    "anchor_generator": "tx",
    "base_rho_g_tol": 1.0e-8,
    "anchor_ratio_gate": 10.0,
    "born": {
        "dchi_bump_center": [0.1, 0.05],
        "dchi_bump_sigma": 0.1,
        "dchi_l2": 1.0e-3,
        "dX_seed": 31415,
        "h_log_min": -3.0,
        "h_log_max": -1.0,
        "h_n": 9,
        "h_clean_min": 3.0e-3,
        "h_clean_max": 1.0e-1,
        "h_ratio": 1.0e-2,
        "A0_A_born_gate": 1.0e-12,
        "born_slope_gate": [1.8, 2.2],
        "born_rel_err_1e-3_gate": 1.0e-4,
        "full_rel_err_1e-2_gate": 1.0e-2,
    },
    "randomness_note": (
        "the only RNG use is the fixed-seed 31415 pose direction in the Born "
        "bilinear check; all gauge/algebra blocks are deterministic"
    ),
}

_EPS = np.finfo(float).eps
_J90 = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=float)


# ---------------------------------------------------------------------------
# Scene / generator helpers (same conventions as Family 3)
# ---------------------------------------------------------------------------

def make_chi0_sym(points: np.ndarray, cfg: dict) -> np.ndarray:
    ss = cfg["symmetric_scene"]
    return (
        ss["amp"]
        * np.exp(-np.sum(points**2, axis=1) / (2.0 * ss["sigma"] ** 2))
    ).astype(float)


def grad_chi0(points: np.ndarray, cfg: dict) -> np.ndarray:
    """Analytic gradient of the two-blob chi0, shape (N^2, 2)."""
    g = np.zeros((points.shape[0], 2), dtype=float)
    blobs = cfg["chi_blobs"]
    for amp, c, sig in [
        (blobs["amp1"], np.asarray(blobs["c1"], dtype=float), blobs["sigma1"]),
        (blobs["amp2"], np.asarray(blobs["c2"], dtype=float), blobs["sigma2"]),
    ]:
        w = amp * np.exp(-np.sum((points - c) ** 2, axis=1) / (2.0 * sig**2))
        g += w[:, None] * (-(points - c) / sig**2)
    return g


def gauge_dchi(grad: np.ndarray, points: np.ndarray) -> dict[str, np.ndarray]:
    """Pixel vectors dchi_g = - (generator field) . grad chi0 (Family-3 form)."""
    jr = (_J90 @ points.T).T
    return {
        "tx": -grad[:, 0],
        "ty": -grad[:, 1],
        "rot": -np.einsum("ni,ni->n", jr, grad),
    }


def gauge_dX(poses: np.ndarray) -> dict[str, np.ndarray]:
    """Pose-direction increments for Tx, Ty and rotation (omega = +1)."""
    T = poses.shape[0]
    tx = np.tile(np.array([1.0, 0.0, 0.0]), T)
    ty = np.tile(np.array([0.0, 1.0, 0.0]), T)
    rot = np.concatenate(
        [np.r_[(_J90 @ pose[:2]), 1.0] for pose in poses]
    )
    return {"tx": tx, "ty": ty, "rot": rot}


def realify_complex_vector(z: np.ndarray) -> np.ndarray:
    return np.sqrt(2.0) * np.concatenate([np.real(z), np.imag(z)])


# ---------------------------------------------------------------------------
# Rank / projector helpers (Family-2 conventions imported and wrapped)
# ---------------------------------------------------------------------------

def range_orthonormal_basis(B: np.ndarray) -> tuple[int, np.ndarray, float]:
    r, U, _, tol = f2.range_basis(B)
    return r, U[:, :r], tol


def symm(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def slam_kernel(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, int, float]:
    """K_SLAM = A^T (I - Z Z^T) A with Z an orthonormal basis of Range(B)."""
    r, Z, tol = range_orthonormal_basis(B)
    P = np.eye(A.shape[0], dtype=float) - Z @ Z.T
    return symm(A.T @ (P @ A)), r, tol


# ---------------------------------------------------------------------------
# Common gauge row / convergence helpers
# ---------------------------------------------------------------------------

def gauge_row_pixel(
    A: np.ndarray,
    B: np.ndarray,
    dchi: np.ndarray,
    dX: np.ndarray,
    generator: str,
    KSL: np.ndarray | None = None,
) -> dict:
    """One pixel-basis gauge row: c = dchi (identity coefficients)."""
    c = dchi.copy()
    Ac = A @ c
    Bd = B @ dX
    num = float(np.linalg.norm(Ac + Bd, ord=2))
    den = float(np.linalg.norm(Ac, ord=2) + np.linalg.norm(Bd, ord=2))
    row = {
        "generator": generator,
        "basis": "pixel",
        "dchi_l2": float(np.linalg.norm(dchi, ord=2)),
        "c_l2": float(np.linalg.norm(c, ord=2)),
        "representation_residual": 0.0,
        "representation_guard": _EPS,
        "map_norm_A_c": float(np.linalg.norm(Ac, ord=2)),
        "pose_norm_B_dX": float(np.linalg.norm(Bd, ord=2)),
        "gauge_numerator": num,
        "gauge_residual": num / max(den, _EPS),
    }
    if KSL is not None:
        row["kernel_distance"] = float(
            np.linalg.norm(KSL @ c, ord=2) / max(np.linalg.norm(c, ord=2), _EPS)
        )
    return row


def fit_coefficients(S: np.ndarray, dchi: np.ndarray) -> tuple[np.ndarray, float]:
    c, _, _, _ = np.linalg.lstsq(S, dchi, rcond=None)
    rep = float(
        np.linalg.norm(S @ c - dchi, ord=2)
        / max(np.linalg.norm(dchi, ord=2), _EPS)
    )
    return c, rep


def smooth_gauge_rows(
    A_smooth: np.ndarray,
    B: np.ndarray,
    S: np.ndarray,
    dchis: dict[str, np.ndarray],
    dXs: dict[str, np.ndarray],
    basis_name: str,
    KSL: np.ndarray | None = None,
) -> list[dict]:
    K = KSL
    if K is None:
        K, _, _ = slam_kernel(A_smooth, B)
    rows = []
    for gen in ["tx", "ty", "rot"]:
        c, rep = fit_coefficients(S, dchis[gen])
        Ac = A_smooth @ c
        Bd = B @ dXs[gen]
        num = float(np.linalg.norm(Ac + Bd, ord=2))
        den = float(np.linalg.norm(Ac, ord=2) + np.linalg.norm(Bd, ord=2))
        row = {
            "generator": gen,
            "basis": basis_name,
            "dchi_l2": float(np.linalg.norm(dchis[gen], ord=2)),
            "c_l2": float(np.linalg.norm(c, ord=2)),
            "representation_residual": rep,
            "representation_guard": _EPS,
            "map_norm_A_c": float(np.linalg.norm(Ac, ord=2)),
            "pose_norm_B_dX": float(np.linalg.norm(Bd, ord=2)),
            "gauge_numerator": num,
            "gauge_residual": num / max(den, _EPS),
        }
        if K is not None:
            row["kernel_distance"] = float(
                np.linalg.norm(K @ c, ord=2) / max(np.linalg.norm(c, ord=2), _EPS)
            )
        rows.append(row)
    return rows


def fit_loglog_slope(
    hs: np.ndarray,
    errs: np.ndarray,
    indices: np.ndarray | list[int],
) -> dict:
    hs = np.asarray(hs, dtype=float)
    errs = np.asarray(errs, dtype=float)
    idx = np.asarray(indices, dtype=int)
    xs = np.log10(hs[idx])
    ys = np.log10(errs[idx])
    slope, intercept = np.polyfit(xs, ys, 1)
    yhat = slope * xs + intercept
    ss_res = float(np.sum((ys - yhat) ** 2))
    ss_tot = float(np.sum((ys - np.mean(ys)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": r2,
        "window_h_start": float(hs[idx[0]]),
        "window_h_end": float(hs[idx[-1]]),
        "window_indices": [int(i) for i in idx],
    }


# ---------------------------------------------------------------------------
# Check A: pixel-basis gauge residual vs N
# ---------------------------------------------------------------------------

def run_check_A(cfg: dict, poses: np.ndarray, rx: np.ndarray, tx: np.ndarray):
    Ns = list(cfg["pixel_gauge_Ns"])
    N_blocks = []
    for N in Ns:
        t0 = time.perf_counter()
        pts, h = hh.make_grid(N)
        chi0 = f1.make_chi0(pts, cfg)
        A, B, _, _ = hh.build_AB(chi0, poses, rx, tx, N, cfg["k_b"])
        A_R, B_R = hh.whiten_realify(A, B, None)
        KSL, rank_B, tolB = slam_kernel(A_R, B_R)
        grad = grad_chi0(pts, cfg)
        dchis = gauge_dchi(grad, pts)
        dXs = gauge_dX(poses)
        rows = []
        for gen in ["tx", "ty", "rot"]:
            dchi = dchis[gen]
            dX = dXs[gen]
            row = gauge_row_pixel(A_R, B_R, dchi, dX, gen, KSL)
            rows.append(row)
        N_blocks.append(
            {
                "N": int(N),
                "grid_h": float(h),
                "n_pixel": int(N * N),
                "rank_B": int(rank_B),
                "rank_B_tol": float(tolB),
                "rows": rows,
                "seconds": time.perf_counter() - t0,
            }
        )
    # Per-generator non-increasing check N=16 -> 32 -> 40 and <= tol.
    ref = {r["generator"]: r["gauge_residual"] for r in N_blocks[0]["rows"]}
    monotone = True
    all_le_tol = True
    comparisons = {}
    for gen, r16 in ref.items():
        vals = [
            next(r["gauge_residual"] for r in block["rows"] if r["generator"] == gen)
            for block in N_blocks
        ]
        gen_ok = all(v <= r16 for v in vals[1:]) and vals[0] <= cfg["gauge_tol_refined"]
        # strict pairwise non-increasing
        pair_ok = all(vals[i + 1] <= vals[i] + 1e-16 for i in range(len(vals) - 1))
        all_le = all(v <= cfg["gauge_tol_refined"] for v in vals)
        monotone = monotone and pair_ok
        all_le_tol = all_le_tol and all_le
        comparisons[gen] = {
            "values": [float(v) for v in vals],
            "vs_N16_leq": gen_ok,
            "pairwise_non_increasing": bool(pair_ok),
            "all_le_1e-3": bool(all_le),
        }
    return {
        "Ns": Ns,
        "blocks": N_blocks,
        "comparisons": comparisons,
        "pass": bool(monotone and all_le_tol),
    }


# ---------------------------------------------------------------------------
# Check B: augmented smooth basis at N=16
# ---------------------------------------------------------------------------

def run_check_B(cfg: dict, poses: np.ndarray, rx: np.ndarray, tx: np.ndarray):
    N = cfg["N"]
    pts, h = hh.make_grid(N)
    chi0 = f1.make_chi0(pts, cfg)
    A, B, _, _ = hh.build_AB(chi0, poses, rx, tx, N, cfg["k_b"])
    A_R, B_R = hh.whiten_realify(A, B, None)
    S0 = f2.build_smooth_basis(pts, cfg["smooth_basis"])

    grad = grad_chi0(pts, cfg)
    dchis = gauge_dchi(grad, pts)
    dXs = gauge_dX(poses)

    gx = grad[:, 0]
    gy = grad[:, 1]
    jr = (_J90 @ pts.T).T
    g_rot = np.einsum("ni,ni->n", jr, grad)
    gcols = {
        "gx": gx / np.linalg.norm(gx),
        "gy": gy / np.linalg.norm(gy),
        "g_rot": g_rot / np.linalg.norm(g_rot),
    }
    S_aug = np.column_stack([S0, gcols["gx"], gcols["gy"], gcols["g_rot"]])
    A_smooth0 = A_R @ S0
    A_smooth_aug = A_R @ S_aug

    KSL0, rank_B0, tolB0 = slam_kernel(A_smooth0, B_R)
    KSL_aug, rank_B_aug, tolB_aug = slam_kernel(A_smooth_aug, B_R)

    original_rows = smooth_gauge_rows(
        A_smooth0, B_R, S0, dchis, dXs, "smooth-p24", KSL0
    )
    augmented_rows = smooth_gauge_rows(
        A_smooth_aug, B_R, S_aug, dchis, dXs, "smooth-p27-aug", KSL_aug
    )

    # Compare the original p=24 recomputation with the raw Family 3 record.
    f3_path = _ROOT / "results" / "family3_gauge_born.json"
    comparison = None
    if f3_path.exists():
        f3 = json.loads(f3_path.read_text())
        f3_rows = {
            r["generator"]: r for r in f3["gauge_two_blob_smooth"]["rows"]
        }
        diffs = {}
        for row in original_rows:
            fr = f3_rows[row["generator"]]
            diffs[row["generator"]] = {
                key: float(row[key] - fr[key])
                for key in [
                    "representation_residual",
                    "map_norm_A_c",
                    "pose_norm_B_dX",
                    "gauge_numerator",
                    "gauge_residual",
                    "kernel_distance",
                ]
                if key in fr and key in row
            }
        comparison = {
            "source": "results/family3_gauge_born.json",
            "diffs": diffs,
            "max_abs_gauge_residual_diff": float(
                max(abs(d["gauge_residual"]) for d in diffs.values())
            ),
        }

    rep_ok = all(r["representation_residual"] < cfg["augmented_representation_tol"]
                 for r in augmented_rows)
    gauge_ok = all(r["gauge_residual"] <= cfg["augmented_gauge_tol"]
                   for r in augmented_rows)
    return {
        "N": N,
        "grid_h": float(h),
        "p0": int(S0.shape[1]),
        "p_aug": int(S_aug.shape[1]),
        "S0": {"columns": int(S0.shape[1]), "normed": True},
        "S_aug": {
            "columns": int(S_aug.shape[1]),
            "appended": ["gx/||gx||", "gy/||gy||", "g_rot/||g_rot||"],
            "p_lt_m": int(S_aug.shape[1]) < cfg["m"],
            "m": cfg["m"],
            "column_norms": [
                float(np.linalg.norm(S_aug[:, j]))
                for j in range(S_aug.shape[1])
            ],
        },
        "rank_B": int(rank_B0),
        "rank_B_tol": float(tolB0),
        "rank_B_aug_tol": float(tolB_aug),
        "original_p24_rows": original_rows,
        "augmented_p27_rows": augmented_rows,
        "family3_comparison": comparison,
        "pass_representation": bool(rep_ok),
        "pass_gauge": bool(gauge_ok),
    }


# ---------------------------------------------------------------------------
# Check C: generator-specific anchor retention (pixel basis, translation-x)
# ---------------------------------------------------------------------------

def run_check_C(cfg: dict, poses: np.ndarray, rx: np.ndarray, tx: np.ndarray):
    N = cfg["N"]
    pts, h = hh.make_grid(N)
    chi0 = f1.make_chi0(pts, cfg)
    A, B, _, _ = hh.build_AB(chi0, poses, rx, tx, N, cfg["k_b"])
    A_R, B_R = hh.whiten_realify(A, B, None)
    grad = grad_chi0(pts, cfg)
    dchis = gauge_dchi(grad, pts)
    dXs = gauge_dX(poses)
    dchi_tx = dchis["tx"]
    dX_tx = dXs["tx"]

    x = pts[:, 0]
    y = pts[:, 1]
    masks = {
        "base": np.ones(pts.shape[0], dtype=bool),
        "anchor": ~((x > 0.30) & (y > 0.30)),
        "known_bg": np.linalg.norm(pts, axis=1) <= 0.35,
        "boundary": (np.abs(x) <= 0.4375) & (np.abs(y) <= 0.4375),
    }
    # Z for Range(B_R) is mask-independent (the full pose Jacobian).
    _, Z, _ = range_orthonormal_basis(B_R)
    P = np.eye(A_R.shape[0], dtype=float) - Z @ Z.T

    rows = []
    for name in ["base", "anchor", "known_bg", "boundary"]:
        mask = masks[name]
        A_red = A_R[:, mask]
        c = dchi_tx[mask]
        K_IS_red = symm(A_red.T @ A_red)
        K_SLAM_red = symm(A_red.T @ (P @ A_red))
        num = float(c @ (K_SLAM_red @ c))
        den = float(c @ (K_IS_red @ c))
        rho_g = num / max(den, _EPS)
        Ac = A_red @ c
        Bd = B_R @ dX_tx
        gauge_red = float(
            np.linalg.norm(Ac + Bd, ord=2)
            / max(np.linalg.norm(Ac, ord=2) + np.linalg.norm(Bd, ord=2), _EPS)
        )
        rows.append(
            {
                "name": name,
                "mask_definition": cfg["anchor_masks"][name],
                "generator": "translation-x",
                "n_free_pixels": int(np.sum(mask)),
                "n_fixed_pixels": int(np.sum(~mask)),
                "c_l2": float(np.linalg.norm(c, ord=2)),
                "K_IS_quad": den,
                "K_SLAM_quad": num,
                "rho_g": rho_g,
                "gauge_residual_red": gauge_red,
            }
        )
    base_rho = rows[0]["rho_g"]
    for row in rows:
        row["rho_g_ratio_vs_base"] = (
            1.0 if row["name"] == "base" else float(row["rho_g"] / base_rho)
        )
    by_name = {r["name"]: r for r in rows}
    return {
        "N": N,
        "generator": "translation-x (Tx)",
        "pixel_basis": True,
        "rows": rows,
        "pass_base_tiny": bool(by_name["base"]["rho_g"] <= cfg["base_rho_g_tol"]),
        "pass_corner_anchor_ratio": bool(
            by_name["anchor"]["rho_g_ratio_vs_base"] > cfg["anchor_ratio_gate"]
        ),
    }


# ---------------------------------------------------------------------------
# Check D: corrected Born empty-background bilinear check
# ---------------------------------------------------------------------------

def run_check_D(cfg: dict, poses: np.ndarray, rx: np.ndarray, tx: np.ndarray):
    N = cfg["N"]
    T = cfg["T"]
    X0 = poses.reshape(-1)
    pts, h = hh.make_grid(N)
    zeros = np.zeros(N * N, dtype=float)
    born_cfg = cfg["born"]

    # ---- D.1 A0(full wave at chi=0) == Born map ---------------------------
    A0, B0, _, _ = hh.build_AB(zeros, poses, rx, tx, N, cfg["k_b"])
    F_born0, A_born = hh.born_forward(zeros, poses, rx, tx, N, cfg["k_b"])
    diff_F = A0 - A_born
    rel_fro_global = float(
        np.linalg.norm(diff_F, ord="fro")
        / max(np.linalg.norm(A_born, ord="fro"), _EPS)
    )
    per_pose = []
    for t in range(T):
        sl = slice(t * cfg["n_rx"], (t + 1) * cfg["n_rx"])
        per_pose.append(
            float(
                np.linalg.norm(diff_F[sl, :], ord="fro")
                / max(np.linalg.norm(A_born[sl, :], ord="fro"), _EPS)
            )
        )
    rel_fro_max = float(max(per_pose))

    A0_R, B0_R = hh.whiten_realify(A0, B0, None)
    linear = A0_R @ np.zeros(N * N)  # placeholder replaced below

    # ---- dchi_small, fixed pose direction --------------------------------
    center = np.asarray(born_cfg["dchi_bump_center"], dtype=float)
    sigma_bump = born_cfg["dchi_bump_sigma"]
    dchi_raw = np.exp(
        -np.sum((pts - center) ** 2, axis=1) / (2.0 * sigma_bump**2)
    )
    dchi = born_cfg["dchi_l2"] * dchi_raw / np.linalg.norm(dchi_raw)
    dchi_max_abs = float(np.max(np.abs(dchi)))

    rng = np.random.default_rng(born_cfg["dX_seed"])
    dX = rng.standard_normal(3 * T)
    dX = dX / np.linalg.norm(dX)
    linear = A0_R @ dchi

    def A0_R_at(Xflat: np.ndarray) -> np.ndarray:
        Ax, Bx, _, _ = hh.build_AB(
            zeros, Xflat.reshape(T, 3), rx, tx, N, cfg["k_b"]
        )
        return hh.whiten_realify(Ax, Bx, None)[0]

    def full_real(chi: np.ndarray, Xflat: np.ndarray) -> np.ndarray:
        Fz = hh.forward_measurements(
            chi, Xflat.reshape(T, 3), rx, tx, N, cfg["k_b"]
        )
        return realify_complex_vector(Fz)

    def born_real(chi: np.ndarray, Xflat: np.ndarray) -> np.ndarray:
        Fb, _ = hh.born_forward(chi, Xflat.reshape(T, 3), rx, tx, N, cfg["k_b"])
        return realify_complex_vector(Fb)

    F0_full_R = full_real(dchi, X0)
    R0 = F0_full_R - linear

    hs = np.logspace(born_cfg["h_log_min"], born_cfg["h_log_max"], born_cfg["h_n"])
    rows = []
    for step in hs:
        Xp = X0 + step * dX
        Xm = X0 - step * dX
        A_Rp = A0_R_at(Xp)
        A_Rm = A0_R_at(Xm)
        BIL = ((A_Rp - A_Rm) / (2.0 * step)) @ dchi
        Fp_full = full_real(dchi, Xp)
        Fm_full = full_real(dchi, Xm)
        FD_full = (Fp_full - Fm_full) / (2.0 * step)
        Fp_born = born_real(dchi, Xp)
        Fm_born = born_real(dchi, Xm)
        FD_born = (Fp_born - Fm_born) / (2.0 * step)
        rel_born = float(
            np.linalg.norm(FD_born - BIL, ord=2)
            / max(np.linalg.norm(BIL, ord=2), _EPS)
        )
        rel_full = float(
            np.linalg.norm(FD_full - BIL, ord=2)
            / max(np.linalg.norm(BIL, ord=2), _EPS)
        )
        R_h = Fp_full - linear - step * BIL
        rows.append(
            {
                "h": float(step),
                "BIL_norm": float(np.linalg.norm(BIL, ord=2)),
                "FD_born_norm": float(np.linalg.norm(FD_born, ord=2)),
                "FD_full_norm": float(np.linalg.norm(FD_full, ord=2)),
                "rel_err_born": rel_born,
                "rel_err_full": rel_full,
                "remainder_norm_R_h": float(np.linalg.norm(R_h, ord=2)),
                "R_h_minus_R0_norm": float(np.linalg.norm(R_h - R0, ord=2)),
            }
        )

    hs_arr = np.asarray([r["h"] for r in rows])
    rel_born_arr = np.asarray([r["rel_err_born"] for r in rows])
    rel_full_arr = np.asarray([r["rel_err_full"] for r in rows])
    rdiff_arr = np.asarray([r["R_h_minus_R0_norm"] for r in rows])

    clean_idx = np.where(
        (hs_arr >= born_cfg["h_clean_min"]) & (hs_arr <= born_cfg["h_clean_max"])
    )[0]
    fit_born = fit_loglog_slope(hs_arr, rel_born_arr, clean_idx)
    fit_full = fit_loglog_slope(hs_arr, rel_full_arr, clean_idx)
    fit_rdiff = fit_loglog_slope(hs_arr, rdiff_arr, clean_idx)

    i_floor = int(np.argmin(rel_full_arr))
    small_idx = np.where(hs_arr <= 1.0e-2)[0]
    i_floor_small = int(small_idx[np.argmin(rel_full_arr[small_idx])])
    hq = born_cfg["h_ratio"]
    iq = int(np.argmin(np.abs(hs_arr - hq)))
    BIL_q = rows[iq]["BIL_norm"]
    ratio_q = BIL_q / max(float(np.linalg.norm(linear, ord=2)), _EPS)

    rel_born_1e_3 = next(
        r["rel_err_born"] for r in rows if abs(r["h"] - 1.0e-3) < 1e-12
    )
    rel_full_1e_2 = rows[iq]["rel_err_full"]

    return {
        "N": N,
        "dchi_small": {
            "l2": float(np.linalg.norm(dchi)),
            "max_abs": dchi_max_abs,
            "center": center.tolist(),
            "sigma": sigma_bump,
            "shape": "same Gaussian bump as Family 3, rescaled to L2=1e-3",
        },
        "dX": {
            "seed": born_cfg["dX_seed"],
            "l2": float(np.linalg.norm(dX)),
            "first_entries": [float(x) for x in dX[:6]],
        },
        "A0_equals_A_born": {
            "rel_fro_global": rel_fro_global,
            "rel_fro_max_over_pose_blocks": rel_fro_max,
            "per_pose_rel_fro": per_pose,
            "A0_norm_F": float(np.linalg.norm(A0, ord="fro")),
            "A_born_norm_F": float(np.linalg.norm(A_born, ord="fro")),
        },
        "linear_norm_A0_R_dchi": float(np.linalg.norm(linear, ord=2)),
        "R0_norm": float(np.linalg.norm(R0, ord=2)),
        "rows": rows,
        "clean_window": {
            "h_min": born_cfg["h_clean_min"],
            "h_max": born_cfg["h_clean_max"],
            "indices": [int(i) for i in clean_idx],
        },
        "fits": {
            "rel_err_born": fit_born,
            "rel_err_full": fit_full,
            "R_full_minus_R0": fit_rdiff,
        },
        "full_wave_floor": {
            "value": float(rel_full_arr[i_floor]),
            "h": float(hs_arr[i_floor]),
            "small_h_value": float(rel_full_arr[i_floor_small]),
            "small_h_definition": "min over grid h <= 1e-2",
        },
        "summary_at_1e-2": {
            "h": hq,
            "BIL_norm": BIL_q,
            "linear_norm_A0_R_dchi": float(np.linalg.norm(linear, ord=2)),
            "bilinear_over_linear": ratio_q,
            "rel_err_born": rows[iq]["rel_err_born"],
            "rel_err_full": rel_full_1e_2,
            "R_h_minus_R0_norm": rows[iq]["R_h_minus_R0_norm"],
            "R0_norm": float(np.linalg.norm(R0, ord=2)),
        },
        "summary_at_1e-3": {
            "h": 1.0e-3,
            "rel_err_born": rel_born_1e_3,
        },
    }


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def plot_check_A(a_result: dict, path: Path) -> None:
    gens = ["tx", "ty", "rot"]
    labels = {"tx": "Tx", "ty": "Ty", "rot": "Rot"}
    Ns = [b["N"] for b in a_result["blocks"]]
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 5.2))
    for gen in gens:
        gauges = [
            next(r["gauge_residual"] for r in b["rows"] if r["generator"] == gen)
            for b in a_result["blocks"]
        ]
        kds = [
            next(r["kernel_distance"] for r in b["rows"] if r["generator"] == gen)
            for b in a_result["blocks"]
        ]
        axs[0].semilogy(Ns, gauges, "o-", label=labels[gen], lw=1.7)
        axs[1].semilogy(Ns, kds, "s--", label=labels[gen], lw=1.4, ms=5)
    axs[0].axhline(1e-3, color="crimson", ls=":", lw=1.2, label="tol 1e-3")
    axs[1].set_xlabel("N")
    axs[0].set_xlabel("N")
    axs[1].set_ylabel(r"$||K_{SLAM}\,d\chi_g||_2 / ||d\chi_g||_2$")
    axs[0].set_ylabel("gauge residual (log)")
    axs[0].set_title("Pixel-basis gauge residual vs N\n"
                     r"$||A_R\,d\chi_g+B_R\,dX_g||/(||A_R\,d\chi_g||+||B_R\,dX_g||)$")
    axs[1].set_title("Pixel-basis SLAM-kernel distance vs N")
    for ax in axs:
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=9)
    fig.suptitle("Family 3b (A): discrete pixel forward model is the meaningful gauge space")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_check_C(c_result: dict, path: Path) -> None:
    names = [r["name"] for r in c_result["rows"]]
    rho = [r["rho_g"] for r in c_result["rows"]]
    gauge = [r["gauge_residual_red"] for r in c_result["rows"]]
    ratios = [r["rho_g_ratio_vs_base"] for r in c_result["rows"]]
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 5.2))
    xpos = np.arange(len(names))
    colors = {"base": "0.4", "anchor": "#1f77b4", "known_bg": "#d62728",
              "boundary": "#2ca02c"}
    for i, name in enumerate(names):
        axs[0].bar(xpos[i], rho[i], 0.58, color=colors[name], log=True)
        axs[1].bar(xpos[i], gauge[i], 0.58, color=colors[name], log=True)
    for i, (rr, gg, ratio) in enumerate(zip(rho, gauge, ratios)):
        axs[0].text(xpos[i], max(rr, 1e-18) * 1.5, f"{rr:.2e}\nratio {ratio:.1f}",
                    ha="center", va="bottom", fontsize=7.5)
        axs[1].text(xpos[i], max(gg, 1e-18) * 1.5, f"{gg:.2e}",
                    ha="center", va="bottom", fontsize=7.5)
    axs[0].set_xticks(xpos)
    axs[0].set_xticklabels(names)
    axs[1].set_xticks(xpos)
    axs[1].set_xticklabels(names)
    axs[0].set_ylabel(r"generator retention $\rho_g$ (log)")
    axs[1].set_ylabel("gauge residual (log)")
    axs[0].axhline(1e-8, color="crimson", ls=":", lw=1.0, label="base tol 1e-8")
    for ax in axs:
        ax.grid(True, which="both", axis="y", alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.suptitle(
        "Family 3b (C): generator-specific retention, translation-x pixel gauge "
        "(fixing pixels raises rho_g)"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_check_D(d_result: dict, path: Path) -> None:
    hs = np.asarray([r["h"] for r in d_result["rows"]])
    rb = np.asarray([r["rel_err_born"] for r in d_result["rows"]])
    rf = np.asarray([r["rel_err_full"] for r in d_result["rows"]])
    rd = np.asarray([r["R_h_minus_R0_norm"] for r in d_result["rows"]])
    fig, axs = plt.subplots(2, 2, figsize=(11.0, 8.4))
    axs[0, 0].loglog(hs, rb, "o-", color="#1f77b4", label="rel_err_born")
    axs[0, 0].loglog(hs, rb[0] * (hs / hs[0]) ** 2, "--", color="grey",
                     label=r"$O(h^2)$ through h=1e-3")
    axs[0, 0].set_title(f"Born FD vs BIL (slope {d_result['fits']['rel_err_born']['slope']:.2f})")
    axs[0, 0].set_ylabel(r"$||FD_{born,h}-BIL_h||/||BIL_h||$")
    axs[0, 1].loglog(hs, rf, "s-", color="#d62728", label="rel_err_full")
    axs[0, 1].axhline(d_result["full_wave_floor"]["value"], color="crimson",
                      ls=":", lw=1.0,
                      label=f"floor {d_result['full_wave_floor']['value']:.2e}")
    axs[0, 1].set_title(
        f"Full-wave FD vs BIL (slope {d_result['fits']['rel_err_full']['slope']:.2f})"
    )
    axs[1, 0].loglog(hs, rd, "^--", color="#2ca02c", label=r"$||R(h)-R(0)||$")
    axs[1, 0].loglog(hs, rd[0] * (hs / hs[0]) ** 2, ":", color="grey",
                     label=r"$O(h^2)$ ref")
    axs[1, 0].set_title(
        f"Remainder drift (slope {d_result['fits']['R_full_minus_R0']['slope']:.2f})"
    )
    axs[1, 1].loglog(hs, [r["BIL_norm"] for r in d_result["rows"]], "o-",
                     color="#7f7f7f", label=r"$||BIL_h||$")
    axs[1, 1].axhline(d_result["linear_norm_A0_R_dchi"], color="k", ls=":",
                      label=r"$||A_0\,d\chi_{small}||$")
    axs[1, 1].set_title(
        f"Bilinear norm, ratio@1e-2 = "
        f"{d_result['summary_at_1e-2']['bilinear_over_linear']:.3f}"
    )
    for ax in axs.ravel():
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
    for ax in [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]:
        ax.set_xlabel("h")
    fig.suptitle(
        "Family 3b (D): corrected Born bilinear check, "
        r"$||d\chi_{small}||_2=10^{-3}$"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    t0 = datetime.now(timezone.utc)
    t_start = time.perf_counter()

    poses = f1.build_poses(cfg)
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)
    pts16, h16 = hh.make_grid(cfg["N"])
    chi0_16 = f1.make_chi0(pts16, cfg)

    print("[family3b] check A: pixel-basis gauge residual vs N ...")
    a_result = run_check_A(cfg, poses, rx, tx)
    for block in a_result["blocks"]:
        print("  N=%d:" % block["N"], {
            r["generator"]: round(r["gauge_residual"], 8) for r in block["rows"]
        })

    print("[family3b] check B: augmented smooth basis at N=16 ...")
    b_result = run_check_B(cfg, poses, rx, tx)
    print("  original p=24:", {
        r["generator"]: round(r["gauge_residual"], 6) for r in b_result["original_p24_rows"]
    })
    print("  augmented p=27:", {
        r["generator"]: round(r["gauge_residual"], 8) for r in b_result["augmented_p27_rows"]
    })

    print("[family3b] check C: generator-specific anchor retention ...")
    c_result = run_check_C(cfg, poses, rx, tx)
    print("  rho_g:", {r["name"]: r["rho_g"] for r in c_result["rows"]})
    print("  ratios:", {
        r["name"]: r["rho_g_ratio_vs_base"] for r in c_result["rows"]
    })

    print("[family3b] check D: corrected Born bilinear check ...")
    d_result = run_check_D(cfg, poses, rx, tx)
    print("  A0==A_born relF:", d_result["A0_equals_A_born"]["rel_fro_max_over_pose_blocks"])
    print("  rel_born@1e-3:", d_result["summary_at_1e-3"]["rel_err_born"],
          " rel_full@1e-2:", d_result["summary_at_1e-2"]["rel_err_full"])
    print("  slopes:", d_result["fits"]["rel_err_born"]["slope"],
          d_result["fits"]["R_full_minus_R0"]["slope"])

    # ------------------------------------------------------------------ gates
    gates = {
        "A_pixel_gauge_le_1e-3": bool(a_result["pass"]),
        "A_pixel_gauge_tol": cfg["gauge_tol_refined"],
        "B_augmented_representation_lt_1e-10": bool(
            b_result["pass_representation"]
        ),
        "B_augmented_representation_tol": cfg["augmented_representation_tol"],
        "B_augmented_gauge_le_1e-3": bool(b_result["pass_gauge"]),
        "B_augmented_gauge_tol": cfg["augmented_gauge_tol"],
        "C_base_rho_g_le_1e-8": bool(c_result["pass_base_tiny"]),
        "C_base_rho_g_tol": cfg["base_rho_g_tol"],
        "C_corner_anchor_ratio_gt_10": bool(c_result["pass_corner_anchor_ratio"]),
        "C_anchor_ratio_gate": cfg["anchor_ratio_gate"],
        "D_A0_equals_A_born_lt_1e-12": bool(
            d_result["A0_equals_A_born"]["rel_fro_max_over_pose_blocks"]
            < cfg["born"]["A0_A_born_gate"]
        ),
        "D_A0_A_born_gate": cfg["born"]["A0_A_born_gate"],
        "D_born_slope_in_1.8_2.2": bool(
            cfg["born"]["born_slope_gate"][0]
            <= d_result["fits"]["rel_err_born"]["slope"]
            <= cfg["born"]["born_slope_gate"][1]
        ),
        "D_born_slope_gate_interval": cfg["born"]["born_slope_gate"],
        "D_rel_err_born_at_1e-3_lt_1e-4": bool(
            d_result["summary_at_1e-3"]["rel_err_born"]
            < cfg["born"]["born_rel_err_1e-3_gate"]
        ),
        "D_rel_err_born_1e-3_gate": cfg["born"]["born_rel_err_1e-3_gate"],
        "D_rel_err_full_at_1e-2_lt_1e-2": bool(
            d_result["summary_at_1e-2"]["rel_err_full"]
            < cfg["born"]["full_rel_err_1e-2_gate"]
        ),
        "D_rel_err_full_1e-2_gate": cfg["born"]["full_rel_err_1e-2_gate"],
    }

    source_files = [
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
        "src/family3_gauge_born.py",
        "src/family3b_refinements.py",
    ]
    results = {
        "generated_utc": t0.isoformat(),
        "runner": "src/family3b_refinements.py",
        "command": ".venv/bin/python src/family3b_refinements.py",
        "family": "3b",
        "supplementary_to": "results/family3_gauge_born.json",
        "raw_family3_record_note": (
            "results/family3_gauge_born.json and notes/family3_report.md are "
            "the raw Family 3 record and were not modified or deleted."
        ),
        "self_cell_formula": hh.SELF_CELL_FORMULA,
        "self_cell_formula_version": hh.SELF_CELL_FORMULA_VERSION,
        "runtime_seconds": None,  # filled after plots
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
            "rx_offsets": rx.tolist(),
            "tx_offset": tx.tolist(),
            "poses": [p.tolist() for p in poses],
            "grid_h_cell_N16": h16,
            "chi0_stats": {
                "min": float(chi0_16.min()),
                "max": float(chi0_16.max()),
                "mean": float(chi0_16.mean()),
                "l2": float(np.linalg.norm(chi0_16)),
            },
        },
        "dimensions_common": {
            "m_real_data": cfg["m"],
            "q_pose": cfg["q_pose"],
            "S_N2_N16": cfg["n_pixel"],
            "T": cfg["T"],
            "n_rx": cfg["n_rx"],
        },
        "check_A_pixel_gauge_vs_N": a_result,
        "check_B_augmented_smooth_gauge": b_result,
        "check_C_anchor_generator_retention": c_result,
        "check_D_born_bilinear_corrected": d_result,
        "gates": gates,
        "notes": {
            "A_pixel_gauge_vs_N": (
                "All pixel gauge residuals are <= 1e-3, but the measured "
                "residuals are NOT non-increasing from N=16 to N=40: e.g. Tx "
                "6.613e-5 -> 7.955e-5 -> 8.135e-5.  The kernel distance does "
                "decrease.  The refined record therefore reports A as failed "
                "on the monotonicity condition while passing the absolute "
                "1e-3 magnitude condition."
            ),
            "D_born_slope_metric_degenerate": (
                "FD_born_h and BIL_h are both centred finite differences of "
                "the same linear Born map (A0 == A_born, verified to 7.1e-17 "
                "relative Frobenius), so their difference is roundoff-limited "
                "rather than O(h^2): the requested born-slope-in-[1.8,2.2] "
                "gate cannot be established from this comparison and is "
                "reported false.  rel_err_born <= 1e-4 at h=1e-3 does pass "
                "(as an exact-equality consistency check), and the meaningful "
                "full-wave floor + O(h^2) remainder drift is recorded "
                "separately (floor 2.166e-4, slope 1.986)."
            ),
        },
        "conditioning_note": (
            "A uses the pixel (discrete forward) basis; B augments the 24 "
            "smooth columns with the three normalised generator fields so the "
            "gauge generators live in the smooth map tangent space; C uses "
            "rho_g = c^T K_SLAM_red c / c^T K_IS_red c for the translation-x "
            "pixel generator rather than the global rho_min of a spectrum; "
            "D rescales the scene perturbation to ||dchi_small||_2 = 1e-3 so "
            "the chi^2 remainder is below the O(h^2) bilinear FD error."
        ),
    }

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)
    results_path = results_dir / "family3b_refinements.json"
    results_path.write_text(
        json.dumps(results, indent=2, default=lambda o: o.tolist()) + "\n"
    )

    fig_a = figures_dir / "family3b_gauge_pixel_refinement.png"
    fig_c = figures_dir / "family3b_anchor_generator_retention.png"
    fig_d = figures_dir / "family3b_born_bilinear_refined.png"
    plot_check_A(a_result, fig_a)
    plot_check_C(c_result, fig_c)
    plot_check_D(d_result, fig_d)

    total_runtime = time.perf_counter() - t_start
    results["runtime_seconds"] = total_runtime
    results_path.write_text(
        json.dumps(results, indent=2, default=lambda o: o.tolist()) + "\n"
    )

    print("\n===== FAMILY 3B REFINEMENT SUMMARY =====")
    print("A pass:", gates["A_pixel_gauge_le_1e-3"],
          "| B pass:", gates["B_augmented_representation_lt_1e-10"],
          gates["B_augmented_gauge_le_1e-3"],
          "| C pass:", gates["C_base_rho_g_le_1e-8"],
          gates["C_corner_anchor_ratio_gt_10"],
          "| D pass:", gates["D_A0_equals_A_born_lt_1e-12"],
          gates["D_born_slope_in_1.8_2.2"],
          gates["D_rel_err_born_at_1e-3_lt_1e-4"],
          gates["D_rel_err_full_at_1e-2_lt_1e-2"])
    print("results ->", results_path)
    print("figure  ->", fig_a)
    print("figure  ->", fig_c)
    print("figure  ->", fig_d)
    print(f"total runtime {total_runtime:.2f}s")


if __name__ == "__main__":
    main()
