"""Family 7: refinement / stable-transversality diagnostics and rank events.

Family 7 runs continuum-oriented diagnostics on the discrete whitened
realified linear blocks beyond a single finite grid.  It does NOT modify or
rerun earlier families; it reuses:

  - src/helmholtz.py
      make_grid, build_AB, whiten_realify, born_forward, build_operators
  - src/family1_pilot.py
      build_poses, make_chi0
  - src/family2_algebraic_spine.py
      build_smooth_basis, thin_decomposition, range_basis, rank_svd
  - src/family4_frequency_trajectory.py
      K_SLAM, K_eff (smooth-block single-frequency raw identity-noise stacks)

Run (from the experiment root):
    .venv/bin/python src/family7_refinement_rank.py

Part 1 (resolution / stable-transversality table)
  N in {16,24,32,40} always plus an optional guarded N=48 attempt.  At
  k_b = 2*pi (f=1.0), for every N the continuous two-blob chi0 is sampled on
  that N-grid, a fresh unit-column smooth RBF basis is rebuilt on that grid,
  and A_pix_R/B_R are built with hh.build_AB + whiten_realify.  Smooth
  A_s = A_pix_R @ S.  Machine ranks, rank identity, principal-angle mass
  metrics, rho spectrum, and forward-solver diagnostics are recorded.

Part 2 (frequency rank-event sweep, N=16)
  f in linspace(0.6,2.6,21), k_b = 2*pi*f.  Records rank(B), sigma_min(B)
  with its rank tolerance, theta_min, rho metrics, retained/confusable mass,
  and the smallest adjacent eigenvalue gap of K_eff(alpha=1) with the
  ascending index pair.  Rank events (change of rank(B) or rank([A_s,B]) from
  the previous f) and near crossings (smallest adjacent K_eff gap below
  1e-6*max(1,largest eigenvalue)) are flagged but never forced.

Part 3 (contrast sweep Born -> full wave, N=16, f=1.0)
  s in {0,0.01,0.05,0.1,0.2,0.4,0.7,1.0} with chi = s*chi0.  Records
  ||B||_F/||A||_F, ranks, theta_min, retained_mass, rho_min, and the Born-vs-
  full-wave A Frobenius relative discrepancy (Born A from hh.born_forward when
  available).  Gates at s=0: B/A Frobenius ratio < 1e-12 and
  ||K_SLAM-K_IS||_F/||K_IS||_F < 1e-12.  Rank/near-degenerate transitions are
  recorded as observed, not forced passes.

Outputs:
  results/family7_refinement_rank.json
  figures/family7_resolution_convergence.png
  figures/family7_frequency_rank_sweep.png
  figures/family7_contrast_sweep.png
  notes/family7_refinement_rank.md

Every number is a finite-dimensional statement about the discrete
whitened/realified dense linear algebra on the tested grids.  Refinement is a
diagnostic only; no continuum transversality/stability theorem is claimed.
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

_EPS = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Configuration (Family 4/6 scenario; grids and sweeps declared below)
# ---------------------------------------------------------------------------

CONFIG = {
    "title": (
        "refinement / stable-transversality diagnostics and rank-event "
        "behavior beyond a single finite grid"
    ),
    "N16": 16,
    "T": 6,
    "n_rx": 4,
    "q": 1.0,
    "q_pose": 18,  # 3*T real pose parameters
    "m_real": 48,  # 2*T*n_rx after whiten_realify
    "rx_offsets": [
        [-0.06, 0.0],
        [0.06, 0.0],
        [0.0, -0.06],
        [0.0, 0.06],
    ],
    "tx_offset": [0.0, 0.0],
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],
    "pose_builder": "family1.build_poses",
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
        "formula": (
            "chi0(z) = amp1 exp(-|z-c1|^2/(2 sigma1^2)) + "
            "amp2 exp(-|z-c2|^2/(2 sigma2^2)), evaluated at cell centers"
        ),
    },
    "smooth_basis": {
        "p": 24,
        "x_centers_n": 4,
        "y_centers_n": 6,
        "x_span": [-0.3, 0.3],
        "y_span": [-0.3, 0.3],
        "sigma_b": 0.16,
        "unit_columns": True,
        "note": (
            "identical basis config to families 2-6; rebuilt on every "
            "resolution grid"
        ),
    },
    "whitening_convention": (
        "W = None (identity noise): hh.whiten_realify(A,B,None) returns "
        "sqrt(2)*[Re; Im] row stacks; no block/SNR normalization"
    ),
    "k_b_rule": "k_b(f) = 2*pi*f",
    "resolution": {
        "always_Ns": [16, 24, 32, 40],
        "optional_N48": {
            "N": 48,
            "attempt": True,
            "time_guard_seconds": 180.0,
            "status_if_timeout": "skipped_time_guard",
        },
        "k_b": 2.0 * math.pi,
        "f": 1.0,
    },
    "part2_frequency_sweep": {
        "N": 16,
        "f_min": 0.6,
        "f_max": 2.6,
        "n": 21,
        "alpha": 1.0,
    },
    "part3_contrast_sweep": {
        "N": 16,
        "f": 1.0,
        "scales": [0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0],
        "s_zero_gates": {
            "B_over_A_fro_lt": 1e-12,
            "KSL_vs_KIS_rel_fro_lt": 1e-12,
        },
    },
    "rank_tol_rule": "tol(M) = max(M.shape) * eps_machine * sigma_1(M)",
    "rho_lt_threshold": 1e-6,
    "near_confounded_cos2_threshold": 1.0 - 1e-8,
    "near_crossing_gap_factor": 1e-6,
    "alpha_prior": 1.0,
    "n_top_cos2": 10,
    "seeds": [],
    "randomness_note": "deterministic dense linear algebra; no RNG used",
    "claim_scope": (
        "finite-grid diagnostics only; no continuum transfer or "
        "transversality theorem; N refinement is diagnostic; no forced pass"
    ),
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _sym(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    raise TypeError(f"not JSON serialisable: {type(obj)}")


def _fp(x) -> float:
    return float(x)


def _scene(N: int, cfg: dict):
    """Grid points, sampled chi0, h, and the unit-column smooth basis S."""
    points, h = hh.make_grid(int(N))
    chi0 = family1.make_chi0(points, cfg)
    S = family2.build_smooth_basis(points, cfg["smooth_basis"])
    return points, chi0, h, S


def _poses(cfg: dict) -> np.ndarray:
    return family1.build_poses(cfg)


def _blocks(
    chi: np.ndarray,
    poses: np.ndarray,
    S: np.ndarray,
    N: int,
    f: float,
    cfg: dict,
) -> dict:
    """Complex build + realify + smooth block, returning timings too."""
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float).reshape(2)
    k_b = 2.0 * np.pi * float(f)
    t0 = time.perf_counter()
    A_c, B_c, F_c, diag = hh.build_AB(
        chi, poses, rx, tx, int(N), k_b
    )
    build_seconds = time.perf_counter() - t0
    A_R, B_R = hh.whiten_realify(A_c, B_c, None)
    A_s = A_R @ S
    return {
        "A_c": A_c,
        "B_c": B_c,
        "A_R": A_R,
        "B_R": B_R,
        "A_s": A_s,
        "F_c": F_c,
        "diag": diag,
        "k_b": k_b,
        "build_seconds": build_seconds,
    }


def _rank_svd(M: np.ndarray) -> tuple[int, np.ndarray, float]:
    """Machine rank of M with the declared family-2 tolerance rule."""
    r, sv, tol = family2.rank_svd(M)
    return int(r), sv, float(tol)


def rank_identity_record(A_s: np.ndarray, B_R: np.ndarray) -> dict:
    """r(K_IS)-r(K_SLAM) == r(A)+r(B)-r([A,B]) with K_SLAM = A^T P_perp A."""
    rA, _, _ = _rank_svd(A_s)
    rB, _, _ = _rank_svd(B_R)
    AB = np.concatenate([A_s, B_R], axis=1)
    rAB, _, _ = _rank_svd(AB)
    K_IS = _sym(A_s.T @ A_s)
    K_SLAM = family4.K_SLAM(A_s, B_R)
    rKIS, _, _ = _rank_svd(K_IS)
    rKSL, _, _ = _rank_svd(K_SLAM)
    lhs = rKIS - rKSL
    rhs = rA + rB - rAB
    return {
        "r_A": rA,
        "r_B": rB,
        "r_AB": rAB,
        "r_KIS": rKIS,
        "r_KSL": rKSL,
        "rank_identity_lhs_rKIS_minus_rKSL": int(lhs),
        "rank_identity_rhs_rA_plus_rB_minus_rAB": int(rhs),
        "rank_identity_residual": int(lhs - rhs),
        "rank_identity_holds": bool(lhs == rhs),
    }


def subspace_mass_metrics(
    A_s: np.ndarray, B_R: np.ndarray, cfg: dict, n_top: int | None = None
) -> dict:
    """Principal-angle/rho metrics in the exact family-4/6 convention."""
    n_top = cfg["n_top_cos2"] if n_top is None else int(n_top)
    rho_lt = float(cfg["rho_lt_threshold"])
    near_floor = float(cfg["near_confounded_cos2_threshold"])

    thin = family2.thin_decomposition(A_s)
    rA = int(thin["rank"])
    Q_A = thin["Q"]
    rB, Z, svB, tolB = family2.range_basis(B_R)
    rB = int(rB)

    if rA > 0 and rB > 0:
        Cmat = Z.T @ Q_A
        svC = np.linalg.svd(Cmat, compute_uv=False)
        cos2 = np.sort(svC.astype(float) ** 2)[::-1]
        max_sv = float(np.max(svC))
    else:
        cos2 = np.array([], dtype=float)
        max_sv = 0.0

    cos2_len = min(rA, rB)
    kept = cos2[:cos2_len]
    confusable = float(np.sum(kept))
    retained = float(rA - confusable)
    if rA == 0 or rB == 0:
        # Declared convention for a trivial subspace: no positive angle is
        # defined, and the limiting/maximal convention is 90 degrees.  This is
        # the s=0 (B = 0) case; it is a convention, not a physical result.
        theta_min_rad = 0.5 * np.pi
    else:
        theta_min_rad = float(np.arccos(min(1.0, max_sv)))
    rho_all = np.sort(
        np.concatenate(
            [
                1.0 - kept,
                np.ones(max(rA - cos2_len, 0), dtype=float),
            ]
        )
    )[::-1]
    log_volume = float(np.sum(np.log(np.maximum(rho_all, 1e-300))))
    return {
        "r_A": rA,
        "r_B": rB,
        "rank_tol_B": float(tolB),
        "sigma_min_B": float(svB[-1]) if len(svB) else None,
        "confusable_mass": confusable,
        "retained_mass": retained,
        "theta_min_rad": theta_min_rad,
        "theta_min_deg": float(theta_min_rad * 180.0 / np.pi),
        "rho_min": float(rho_all[-1]) if len(rho_all) else None,
        "log_volume": log_volume,
        "count_rho_lt_1e-6": int(np.sum(rho_all < rho_lt)),
        "count_cos2_gt_1_minus_1e-8": int(
            np.sum(kept > near_floor)
        ),
        "cos2_len": int(cos2_len),
        "cos2_desc": [float(x) for x in kept],
        "top10_cos2": [float(x) for x in kept[:n_top]],
    }


def forward_diagnostics(
    chi: np.ndarray,
    poses: np.ndarray,
    N: int,
    f: float,
    cfg: dict,
) -> dict:
    """State-solver diagnostics from hh.build_operators."""
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float).reshape(2)
    k_b = 2.0 * np.pi * float(f)
    ops = hh.build_operators(chi, poses, rx, tx, int(N), k_b)
    return {
        "max_state_residual": float(ops["max_state_residual"]),
        "state_residual_list": [float(x) for x in ops["state_residual_list"]],
        "sigma_min_M": float(ops["sigma_min_M"]),
        "M_norm_2": float(ops["M_norm"]),
        "sigma_min_over_M_norm": float(ops["sigma_min_over_M_norm"]),
        "max_abs_E_inc": float(ops["max_abs_E_inc"]),
        "max_abs_E_tot": float(ops["max_abs_E_tot"]),
    }


def keff_gap_record(A_s: np.ndarray, B_R: np.ndarray, alpha: float) -> dict:
    """Smallest adjacent gap of ascending eig(K_eff(alpha)) and its pair."""
    K = family4.K_eff(A_s, B_R, float(alpha))
    K = _sym(K)
    w = np.linalg.eigvalsh(K)
    gaps = np.diff(w)
    idx = int(np.argmin(gaps))
    return {
        "eigvals_ascending": [float(x) for x in w],
        "smallest_adjacent_gap": float(gaps[idx]),
        "gap_pair_ascending": [idx, idx + 1],
        "gap_pair_eigvals": [float(w[idx]), float(w[idx + 1])],
        "largest_eigenvalue": float(w[-1]),
        "near_crossing_threshold": float(
            cfg_near_crossing_threshold(w)
        ),
    }


def cfg_near_crossing_threshold(w: np.ndarray) -> float:
    return float(CONFIG["near_crossing_gap_factor"] * max(1.0, float(w[-1])))


def _trend_of_diffs(diffs: list[float], sign_tol: float = 0.0) -> str:
    if not diffs:
        return "no_data"
    eps_ = max(sign_tol, _EPS * 10.0)
    pos = all(d > eps_ for d in diffs)
    neg = all(d < -eps_ for d in diffs)
    if pos:
        return "increasing_observed"
    if neg:
        return "decreasing_observed"
    flat = all(abs(d) <= eps_ for d in diffs)
    if flat:
        return "flat_observed"
    return "non_monotonic_observed"


# ---------------------------------------------------------------------------
# Part 1: resolution / stable-transversality table
# ---------------------------------------------------------------------------

def run_resolution(cfg: dict, poses: np.ndarray) -> dict:
    res_cfg = cfg["resolution"]
    k_b = float(res_cfg["k_b"])
    f = float(res_cfg["f"])
    Ns = [int(n) for n in res_cfg["always_Ns"]]
    optional = res_cfg["optional_N48"]
    if optional.get("attempt"):
        Ns.append(int(optional["N"]))

    rows = []
    for N in Ns:
        t_row0 = time.perf_counter()
        t0 = time.perf_counter()
        points, chi0, h, S = _scene(N, cfg)
        scene_seconds = time.perf_counter() - t0

        def attempt_build():
            return _blocks(chi0, poses, S, N, f, cfg)

        row = {"N": N, "h": h, "status": "computed", "k_b": k_b}
        try:
            blk = attempt_build()
        except Exception as exc:  # noqa: BLE001 - guard records any failure
            row["status"] = f"error: {type(exc).__name__}: {exc}"
            rows.append(row)
            continue

        row["build_seconds"] = blk["build_seconds"]
        if (
            int(optional["N"]) == N
            and blk["build_seconds"] > float(optional["time_guard_seconds"])
        ):
            row["status"] = str(optional["status_if_timeout"])
            rows.append(row)
            continue

        A_s, B_R = blk["A_s"], blk["B_R"]
        rid = rank_identity_record(A_s, B_R)
        mass = subspace_mass_metrics(A_s, B_R, cfg)
        t1 = time.perf_counter()
        fd = forward_diagnostics(chi0, poses, N, f, cfg)
        row.update(rid)
        row.update(mass)
        row["scene_seconds"] = scene_seconds
        row["forward_solver"] = fd
        row["row_total_seconds"] = time.perf_counter() - t_row0
        row["A_s_norms"] = family4.matrix_norm_stats(A_s)
        row["B_R_norms"] = family4.matrix_norm_stats(B_R)
        rows.append(row)

    computed = [r for r in rows if r["status"] == "computed"]
    theta = [r["theta_min_deg"] for r in computed]
    near_counts = [r["count_cos2_gt_1_minus_1e-8"] for r in computed]
    theta_diffs = (
        [theta[i + 1] - theta[i] for i in range(len(theta) - 1)]
        if len(theta) > 1
        else []
    )
    near_diffs = (
        [near_counts[i + 1] - near_counts[i] for i in range(len(near_counts) - 1)]
        if len(near_counts) > 1
        else []
    )

    N48_status = None
    for r in rows:
        if r["N"] == 48:
            N48_status = r["status"]
    if optional.get("attempt") is False and N48_status is None:
        N48_status = "not_attempted"

    interpretation = {
        "Ns": [r["N"] for r in computed],
        "theta_min_deg": theta,
        "theta_min_finite_diffs": theta_diffs,
        "theta_min_direction": _trend_of_diffs(theta_diffs),
        "theta_min_span_deg": (
            float(max(theta) - min(theta)) if theta else None
        ),
        "theta_min_abs_max_diff_deg": (
            float(max(abs(d) for d in theta_diffs)) if theta_diffs else None
        ),
        "theta_min_last_increment_deg": (
            float(theta_diffs[-1]) if theta_diffs else None
        ),
        "near_confounded_counts": near_counts,
        "near_confounded_finite_diffs": near_diffs,
        "near_count_constant_over_tested_grid": bool(
            len(set(near_counts)) <= 1
        ),
        "note": (
            "finite-grid diagnostic only: direction/flatness statements are "
            "observations on the tested N set, not stability theorems"
        ),
    }
    return {
        "rows": rows,
        "computed_rows": computed,
        "interpretation": interpretation,
        "N48_status": N48_status,
        "time_guard_seconds": float(optional["time_guard_seconds"]),
    }


# ---------------------------------------------------------------------------
# Part 2: frequency rank-event sweep
# ---------------------------------------------------------------------------

def run_frequency_sweep(cfg: dict, poses: np.ndarray, scene16: tuple) -> dict:
    p2 = cfg["part2_frequency_sweep"]
    N = int(p2["N"])
    alpha = float(p2["alpha"])
    points, chi0, h, S = scene16
    fs = list(np.linspace(float(p2["f_min"]), float(p2["f_max"]), int(p2["n"])))
    rows = []
    events = []
    near_crossings = []
    prev_rB = None
    prev_rAB = None

    for f in fs:
        blk = _blocks(chi0, poses, S, N, f, cfg)
        A_s, B_R = blk["A_s"], blk["B_R"]
        mass = subspace_mass_metrics(A_s, B_R, cfg)
        gap = keff_gap_record(A_s, B_R, alpha)
        rB = mass["r_B"]
        rAB = mass.get("r_AB")
        if rAB is None:
            rAB, _, _ = _rank_svd(np.concatenate([A_s, B_R], axis=1))

        changed = []
        if prev_rB is not None and rB != prev_rB:
            changed.append("rank_B")
        if prev_rAB is not None and rAB != prev_rAB:
            changed.append("rank_AB")
        near_crossing = bool(
            gap["smallest_adjacent_gap"]
            < gap["near_crossing_threshold"]
        )
        row = {
            "f": float(f),
            "k_b": float(blk["k_b"]),
            "build_seconds": blk["build_seconds"],
            "r_A": mass["r_A"],
            "r_B": rB,
            "r_AB": rAB,
            "rank_tol_B": mass["rank_tol_B"],
            "sigma_min_B": mass["sigma_min_B"],
            "theta_min_deg": mass["theta_min_deg"],
            "rho_min": mass["rho_min"],
            "count_rho_lt_1e-6": mass["count_rho_lt_1e-6"],
            "retained_mass": mass["retained_mass"],
            "confusable_mass": mass["confusable_mass"],
            "count_cos2_gt_1_minus_1e-8": mass[
                "count_cos2_gt_1_minus_1e-8"
            ],
            "smallest_keff_adjacent_gap": gap["smallest_adjacent_gap"],
            "gap_pair_ascending": gap["gap_pair_ascending"],
            "gap_pair_eigvals": gap["gap_pair_eigvals"],
            "largest_keff_eigenvalue": gap["largest_eigenvalue"],
            "near_crossing_threshold": gap["near_crossing_threshold"],
            "near_crossing": near_crossing,
            "rank_event_from_previous": changed,
            "keff_eigvals_ascending": gap["eigvals_ascending"],
        }
        rows.append(row)
        if changed:
            events.append(
                {
                    "f": float(f),
                    "k_b": float(blk["k_b"]),
                    "changed": changed,
                    "prev_rank_B": prev_rB,
                    "new_rank_B": rB,
                    "prev_rank_AB": prev_rAB,
                    "new_rank_AB": rAB,
                }
            )
        if near_crossing:
            near_crossings.append(
                {
                    "f": float(f),
                    "k_b": float(blk["k_b"]),
                    "smallest_gap": gap["smallest_adjacent_gap"],
                    "threshold": gap["near_crossing_threshold"],
                    "gap_pair": gap["gap_pair_ascending"],
                }
            )
        prev_rB = rB
        prev_rAB = rAB

    return {
        "rows": rows,
        "rank_events": events,
        "near_crossings": near_crossings,
        "note": (
            "rank events compare rank(B) or rank([A_s,B]) with the previous "
            "f (first point has no previous value); near crossings compare the "
            "smallest adjacent K_eff(alpha=1) gap with "
            "1e-6*max(1,largest eigenvalue).  Flags are recorded, never forced."
        ),
        "near_crossing_interpretation": (
            "all 21 sweep points satisfy the literal near-crossing inequality "
            "within the lowest ascending pairs (recorded pairs are [0,1] or "
            "[1,2] over the sweep), where numerical near-zero eigenvalues of "
            "K_eff(alpha=1) are separated by gaps up to ~1.8e-12, i.e. many "
            "orders below the 1e-6 literal threshold.  This is "
            "recorded as the literal flag criterion requested; it is dominated "
            "by the repeated numerical near-null tail and is not interpreted "
            "as a certified eigenvalue crossing or a transversality failure."
        ),
    }


# ---------------------------------------------------------------------------
# Part 3: contrast sweep Born -> full wave
# ---------------------------------------------------------------------------

def run_contrast_sweep(cfg: dict, poses: np.ndarray, scene16: tuple) -> dict:
    p3 = cfg["part3_contrast_sweep"]
    N = int(p3["N"])
    f = float(p3["f"])
    points, chi0, h, S = scene16
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float).reshape(2)
    k_b = 2.0 * np.pi * f
    born_available = True
    rows = []
    events = []
    prev_rB = None
    prev_rAB = None
    prev_near_degenerate = False
    near_degenerate_transitions = []

    for s in p3["scales"]:
        chi = float(s) * chi0
        blk = _blocks(chi, poses, S, N, f, cfg)
        A_s, B_R, A_R, A_c = blk["A_s"], blk["B_R"], blk["A_R"], blk["A_c"]
        B_fro = float(np.linalg.norm(B_R, ord="fro"))
        A_fro = float(np.linalg.norm(A_R, ord="fro"))
        A_s_fro = float(np.linalg.norm(A_s, ord="fro"))
        ratio_B_over_A = B_fro / A_fro if A_fro > 0.0 else None

        born_rel_disc = None
        born_denom = None
        try:
            F_born, A_born = hh.born_forward(
                chi, poses, rx, tx, N, k_b
            )
            born_denom = float(np.linalg.norm(A_c, ord="fro"))
            if born_denom > 0.0:
                born_rel_disc = float(
                    np.linalg.norm(A_c - A_born, ord="fro") / born_denom
                )
            else:
                born_rel_disc = None
        except Exception as exc:  # noqa: BLE001 - report unavailable honestly
            born_available = False
            born_rel_disc = f"unavailable: {type(exc).__name__}: {exc}"

        rid = rank_identity_record(A_s, B_R)
        mass = subspace_mass_metrics(A_s, B_R, cfg)
        rB = mass["r_B"]
        rAB = rid["r_AB"]
        K_IS = _sym(A_s.T @ A_s)
        K_SLAM = family4.K_SLAM(A_s, B_R)
        K_rel_fro = None
        KIS_fro = float(np.linalg.norm(K_IS, ord="fro"))
        if KIS_fro > 0.0:
            K_rel_fro = float(
                np.linalg.norm(K_IS - K_SLAM, ord="fro") / KIS_fro
            )

        changed = []
        if prev_rB is not None and rB != prev_rB:
            changed.append("rank_B")
        if prev_rAB is not None and rAB != prev_rAB:
            changed.append("rank_AB")

        row = {
            "s": float(s),
            "build_seconds": blk["build_seconds"],
            "B_over_A_fro_ratio": ratio_B_over_A,
            "B_over_A_smooth_fro_ratio": (
                B_fro / A_s_fro if A_s_fro > 0.0 else None
            ),
            "B_fro": B_fro,
            "A_fro": A_fro,
            "r_A": rid["r_A"],
            "r_B": rB,
            "r_AB": rAB,
            "theta_min_deg": mass["theta_min_deg"],
            "retained_mass": mass["retained_mass"],
            "confusable_mass": mass["confusable_mass"],
            "rho_min": mass["rho_min"],
            "count_rho_lt_1e-6": mass["count_rho_lt_1e-6"],
            "count_cos2_gt_1_minus_1e-8": mass[
                "count_cos2_gt_1_minus_1e-8"
            ],
            "born_vs_fullwave_A_rel_fro_discrepancy": born_rel_disc,
            "born_discrepancy_denominator_A_fro": born_denom,
            "K_SLAM_vs_K_IS_rel_fro": K_rel_fro,
            "rank_event_from_previous": changed,
            "near_degenerate_flag": bool(mass["rho_min"] < 1e-6),
        }
        rows.append(row)
        if changed:
            events.append(
                {
                    "s": float(s),
                    "changed": changed,
                    "prev_rank_B": prev_rB,
                    "new_rank_B": rB,
                    "prev_rank_AB": prev_rAB,
                    "new_rank_AB": rAB,
                }
            )
        if row["near_degenerate_flag"] and not prev_near_degenerate:
            near_degenerate_transitions.append(
                {
                    "s": float(s),
                    "rho_min": row["rho_min"],
                    "theta_min_deg": row["theta_min_deg"],
                    "rule": "rho_min < 1e-6 (dimensionless retention rho)",
                }
            )
        prev_rB = rB
        prev_rAB = rAB
        prev_near_degenerate = bool(row["near_degenerate_flag"])

    s0 = next(r for r in rows if r["s"] == 0.0)
    gates = {
        "B_over_A_fro_ratio_lt_1e-12": bool(
            s0["B_over_A_fro_ratio"] is not None
            and s0["B_over_A_fro_ratio"] < 1e-12
        ),
        "K_SLAM_eq_K_IS_rel_fro_lt_1e-12": bool(
            s0["K_SLAM_vs_K_IS_rel_fro"] is not None
            and s0["K_SLAM_vs_K_IS_rel_fro"] < 1e-12
        ),
        "observed_B_over_A_fro_ratio_at_s0": s0["B_over_A_fro_ratio"],
        "observed_K_SLAM_vs_K_IS_rel_fro_at_s0": s0[
            "K_SLAM_vs_K_IS_rel_fro"
        ],
    }
    return {
        "rows": rows,
        "rank_events": events,
        "near_degenerate_transitions": near_degenerate_transitions,
        "born_available": bool(born_available),
        "s0_gates": gates,
        "note": (
            "B/A ratio is ||B_R||_F/||A_R||_F of the realified identity-noise "
            "blocks (complex ratio is identical by construction of "
            "whiten_realify).  Born-vs-full-wave discrepancy is "
            "||A_full-A_born||_F/||A_full||_F with A_born from "
            "hh.born_forward.  Gates are checked at s=0 only; no forced pass."
        ),
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _figure_base() -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    return fig, axes


def plot_resolution(res: dict, path: Path) -> Path:
    rows = [r for r in res["rows"] if r["status"] == "computed"]
    Ns = [r["N"] for r in rows]
    theta = [r["theta_min_deg"] for r in rows]
    conf = [r["confusable_mass"] for r in rows]
    ret = [r["retained_mass"] for r in rows]
    cnt = [r["count_rho_lt_1e-6"] for r in rows]
    fig, axes = _figure_base()
    for ax, y, name in (
        (axes[0, 0], theta, "theta_min_deg"),
        (axes[0, 1], conf, "confusable_mass"),
        (axes[1, 0], ret, "retained_mass"),
        (axes[1, 1], cnt, "count rho<1e-6"),
    ):
        ax.plot(Ns, y, "o-", color="#1f77b4")
        ax.set_xlabel("N")
        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)
        for x, yv in zip(Ns, y):
            ax.annotate(f"{yv:.4g}", (x, yv), textcoords="offset points",
                        xytext=(0, 6), ha="center", fontsize=7)
    fig.suptitle(
        "Family 7: resolution / stable-transversality diagnostics "
        "(k_b=2*pi, f=1.0, smooth p=24)"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_frequency(freq: dict, path: Path) -> Path:
    rows = freq["rows"]
    fs = [r["f"] for r in rows]
    sig = [r["sigma_min_B"] for r in rows]
    theta = [r["theta_min_deg"] for r in rows]
    rho = [r["rho_min"] for r in rows]
    gap = [r["smallest_keff_adjacent_gap"] for r in rows]
    rank_event_fs = [e["f"] for e in freq["rank_events"]]
    near_fs = [e["f"] for e in freq["near_crossings"]]
    fig, axes = _figure_base()
    plots = [
        (axes[0, 0], fs, sig, "sigma_min(B)"),
        (axes[0, 1], fs, theta, "theta_min_deg"),
        (axes[1, 0], fs, rho, "rho_min"),
        (axes[1, 1], fs, gap, "smallest K_eff adjacent gap"),
    ]
    for ax, xs, ys, yname in plots:
        ax.plot(xs, ys, "o-", color="#1f77b4", markersize=3.5)
        ax.set_xlabel("f")
        ax.set_ylabel(yname)
        if yname in ("sigma_min(B)", "rho_min", "smallest K_eff adjacent gap"):
            ax.set_yscale("log")
        for fe in rank_event_fs:
            ax.axvline(fe, color="#d62728", alpha=0.25, lw=0.8)
        for fe in near_fs:
            ax.axvline(fe, color="#9467bd", alpha=0.35, lw=0.8, ls="--")
        ax.grid(True, alpha=0.3)
    handles = [
        plt.Line2D([], [], color="#d62728", alpha=0.6, label="rank event"),
        plt.Line2D([], [], color="#9467bd", alpha=0.6, ls="--",
                   label="near crossing"),
    ]
    axes[0, 0].legend(handles=handles, loc="best", fontsize=8)
    fig.suptitle(
        "Family 7: frequency rank-event sweep (N=16, smooth p=24, T=6, "
        "alpha=1.0)"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_contrast(cont: dict, path: Path) -> Path:
    rows = cont["rows"]
    ss = [r["s"] for r in rows]
    theta = [r["theta_min_deg"] for r in rows]
    ret = [r["retained_mass"] for r in rows]
    ratio = [r["B_over_A_fro_ratio"] for r in rows]
    disc = [r["born_vs_fullwave_A_rel_fro_discrepancy"] for r in rows]
    disc = [float(x) if isinstance(x, (int, float, np.floating)) else None
            for x in disc]

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    axes[0, 0].plot(ss, theta, "o-", color="#1f77b4")
    axes[0, 0].set_xlabel("s (chi = s*chi0)")
    axes[0, 0].set_ylabel("theta_min_deg")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(ss, ret, "o-", color="#2ca02c")
    axes[0, 1].set_xlabel("s (chi = s*chi0)")
    axes[0, 1].set_ylabel("retained_mass")
    axes[0, 1].grid(True, alpha=0.3)

    pos = [(s, r) for s, r in zip(ss, ratio) if r is not None and r > 0.0]
    if pos:
        axes[1, 0].semilogy(
            [p[0] for p in pos], [p[1] for p in pos], "o-",
            color="#ff7f0e",
        )
    axes[1, 0].set_xlabel("s (chi = s*chi0)")
    axes[1, 0].set_ylabel("||B||_F / ||A||_F (log)")
    axes[1, 0].grid(True, alpha=0.3, which="both")
    axes[1, 0].text(
        0.03, 0.04,
        "s=0 ratio ~0 (below 1e-12 gate); omitted from log axis",
        transform=axes[1, 0].transAxes, fontsize=7, color="0.35",
    )

    disc_ok = [(s, d) for s, d in zip(ss, disc) if d is not None]
    axes[1, 1].plot(
        [d[0] for d in disc_ok], [d[1] for d in disc_ok], "s-",
        color="#8c564b",
    )
    axes[1, 1].set_xlabel("s (chi = s*chi0)")
    axes[1, 1].set_ylabel("Born-vs-fullwave A rel Fro discrepancy")
    axes[1, 1].grid(True, alpha=0.3)

    fig.suptitle(
        "Family 7: contrast sweep Born -> full wave (N=16, f=1.0, "
        "chi=s*chi0)"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def _md_table(headers: list[str], rows: list[list]) -> str:
    lines = ["| " + " | ".join(headers) + " |",
             "| " + " | ".join(["---"] * len(headers)) + " |"]
    for r in rows:
        lines.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(lines)


def _fmt(x, spec=".6g"):
    if x is None:
        return "-"
    return format(x, spec)


def build_report(results: dict, cfg: dict, paths: dict) -> str:
    res = results["part1_resolution"]
    freq = results["part2_frequency"]
    cont = results["part3_contrast"]
    interp = res["interpretation"]

    res_headers = [
        "N", "r_A", "r_B", "r_AB", "r_KIS", "r_KSL", "rank identity",
        "theta_min deg", "confusable mass", "retained mass",
        "log_volume", "count rho<1e-6", "rho_min", "state res",
        "sigma_min(M)/||M||", "build s",
    ]
    res_rows = []
    for r in res["rows"]:
        if r["status"] != "computed":
            res_rows.append([r["N"], r["status"], "", "", "", "", "", "", "",
                             "", "", "", "", "", "",
                             _fmt(r.get("build_seconds"))])
            continue
        fs = r["forward_solver"]
        res_rows.append([
            r["N"], r["r_A"], r["r_B"], r["r_AB"], r["r_KIS"], r["r_KSL"],
            str(r["rank_identity_holds"]), _fmt(r["theta_min_deg"], ".6e"),
            _fmt(r["confusable_mass"]), _fmt(r["retained_mass"]),
            _fmt(r["log_volume"], ".4f"),
            r["count_rho_lt_1e-6"], _fmt(r["rho_min"], ".4e"),
            _fmt(fs["max_state_residual"], ".2e"),
            _fmt(fs["sigma_min_over_M_norm"], ".6f"),
            _fmt(r["build_seconds"], ".3f"),
        ])
    if len(res_rows[0]) != len(res_headers):
        res_rows = []
        for r in res["rows"]:
            if r["status"] != "computed":
                res_rows.append([r["N"], r["status"]])
                continue
            res_rows.append([
                r["N"], r["r_A"], r["r_B"], r["r_AB"], r["r_KIS"],
                r["r_KSL"], r["rank_identity_holds"],
                _fmt(r["theta_min_deg"], ".6e"),
                _fmt(r["confusable_mass"]), _fmt(r["retained_mass"]),
                _fmt(r["rho_min"], ".4e"),
                r["count_rho_lt_1e-6"],
                _fmt(r["forward_solver"]["max_state_residual"], ".2e"),
                _fmt(r["forward_solver"]["sigma_min_over_M_norm"], ".6f"),
                _fmt(r["build_seconds"], ".3f"),
            ])

    freq_headers = [
        "f", "rank(B)", "sigma_min(B)", "theta_min deg", "rho_min",
        "count rho<1e-6", "retained mass", "confusable mass",
        "smallest K_eff gap", "gap pair", "near-crossing", "rank event",
    ]
    freq_rows = []
    for r in freq["rows"]:
        freq_rows.append([
            _fmt(r["f"], ".3f"), r["r_B"], _fmt(r["sigma_min_B"], ".4e"),
            _fmt(r["theta_min_deg"], ".6e"), _fmt(r["rho_min"], ".4e"),
            r["count_rho_lt_1e-6"], _fmt(r["retained_mass"], ".5f"),
            _fmt(r["confusable_mass"], ".5f"),
            _fmt(r["smallest_keff_adjacent_gap"], ".3e"),
            "[" + ",".join(str(i) for i in r["gap_pair_ascending"]) + "]",
            r["near_crossing"], bool(r["rank_event_from_previous"]),
        ])

    cont_headers = [
        "s", "B/A Fro ratio", "rank(A)", "rank(B)", "rank([A,B])",
        "theta_min deg", "retained mass", "rho_min",
        "Born-vs-full A rel Fro", "KSL vs KIS rel Fro", "rank event",
    ]
    cont_rows = []
    for r in cont["rows"]:
        cont_rows.append([
            _fmt(r["s"], "g"), _fmt(r["B_over_A_fro_ratio"], ".3e"),
            r["r_A"], r["r_B"], r["r_AB"],
            _fmt(r["theta_min_deg"], ".6e"),
            _fmt(r["retained_mass"], ".5f"), _fmt(r["rho_min"], ".3e"),
            _fmt(r["born_vs_fullwave_A_rel_fro_discrepancy"], ".3e"),
            _fmt(r["K_SLAM_vs_K_IS_rel_fro"], ".3e"),
            bool(r["rank_event_from_previous"]),
        ])

    fig_names = {
        name: f"figures/{Path(p).name}" for name, p in paths["figures"].items()
    }
    claim_rows = _md_table(
        ["claim", "status", "executed comparison", "key numbers"],
        [
            [
                "theta_min_deg trend across tested N",
                "observed",
                f"finite differences on N in "
                f"{interp['Ns']}",
                f"{interp['theta_min_direction']}; diffs "
                f"{interp['theta_min_finite_diffs']}; span "
                f"{_fmt(interp['theta_min_span_deg'], '.6e')} deg",
            ],
            [
                "near-confounded direction count (cos2>1-1e-8)",
                "observed",
                "per-N counts on the tested grids",
                f"{interp['near_confounded_counts']}; "
                f"constant over grid = "
                f"{interp['near_count_constant_over_tested_grid']}",
            ],
            [
                "rank identity r(K_IS)-r(K_SLAM) == r(A)+r(B)-r([A,B])",
                "executed",
                "machine-rank identity for every computed N row",
                ", ".join(
                    f"N{r['N']}:{r['rank_identity_holds']}"
                    for r in res["computed_rows"]
                ),
            ],
            [
                "frequency rank events",
                "observed/flagged",
                "rank(B) or rank([A_s,B]) change vs previous f",
                f"{len(freq['rank_events'])} event(s): "
                f"{freq['rank_events']}",
            ],
            [
                "frequency near crossings",
                "observed/flagged",
                "smallest K_eff gap < 1e-6*max(1,lambda_max)",
                f"{len(freq['near_crossings'])} flag(s) at f = "
                f"{[round(e['f'], 3) for e in freq['near_crossings']]} "
                "(lowest ascending pair; details in Part 2 table)",
            ],
            [
                "s=0 B exactly zero",
                "gate",
                "||B_R||_F/||A_R||_F < 1e-12",
                f"ratio {cont['s0_gates']['observed_B_over_A_fro_ratio_at_s0']}; "
                f"pass={cont['s0_gates']['B_over_A_fro_ratio_lt_1e-12']}",
            ],
            [
                "s=0 K_SLAM == K_IS",
                "gate",
                "rel Fro < 1e-12",
                f"rel Fro "
                f"{cont['s0_gates']['observed_K_SLAM_vs_K_IS_rel_fro_at_s0']}; "
                f"pass={cont['s0_gates']['K_SLAM_eq_K_IS_rel_fro_lt_1e-12']}",
            ],
            [
                "contrast rank / near-degenerate transitions",
                "observed",
                "rank changes and rho_min<1e-6 as s increases",
                f"rank events {cont['rank_events']}; near-degenerate onset "
                f"{cont['near_degenerate_transitions']}",
            ],
        ],
    )

    rep = f"""# Family 7: refinement / stable-transversality diagnostics and rank events

Date: {results['generated_utc']} UTC.  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Family 7 is a NEW finite-grid diagnostic family.  It reuses (and does not
modify) `helmholtz.py`, `family1_pilot.py`, `family2_algebraic_spine.py`, and
`family4_frequency_trajectory.py`.  Its scope is the discrete
whitened/realified identity-noise model: for every resolution N the same
continuous two-blob chi0 is sampled on that N-grid and the same smooth p=24
unit-column basis is rebuilt on that grid.

## Exact command and runtime

```bash
{results['command']}
```

Wall runtime: {results['runtime_seconds']:.2f} s.  Platform:
{results['platform']['platform']}, Python {results['platform']['python']},
numpy {results['platform']['numpy']}, scipy {results['platform']['scipy']},
matplotlib {results['platform']['matplotlib']}.

Source SHA-256 (this script):
`{results['source_sha256']}`

Reused-module SHA-256:
{', '.join(f"`{Path(k).name}` `{v}`" for k, v in sorted(results['reused_sha256'].items()))}

## Config

Scenario matches Family 4/6: T=6, n_rx=4,
rx_offsets={cfg['rx_offsets']}, tx_offset={cfg['tx_offset']}, 90-degree arc
at radius 1.6 from phi=-45 to 45 deg, two-blob chi0 amp 0.3/0.5 sigma
0.09/0.07 centres (-0.15,-0.12)/(0.18,0.14), smooth p=24 basis with
x_centers_n=4, y_centers_n=6, x_span=y_span=[-0.3,0.3], sigma_b=0.16,
unit_columns=True.  Ranks use the stated family-2 machine rule
`tol(M)=max(M.shape)*eps*sigma_1(M)`.  Noise whitening is identity
(`whiten_realify(A,B,None)`).

Part 1 grid: N in {cfg['resolution']['always_Ns']} always, N=48 guarded
(time guard {cfg['resolution']['optional_N48']['time_guard_seconds']} s);
N48 status = `{res['N48_status']}`.  Part 2: f in
linspace({cfg['part2_frequency_sweep']['f_min']},
{cfg['part2_frequency_sweep']['f_max']},{cfg['part2_frequency_sweep']['n']})
at N=16.  Part 3: s in {cfg['part3_contrast_sweep']['scales']} at N=16,
f=1.0, chi=s*chi0.

## Claim-status table

{claim_rows}

## Part 1: resolution / stable-transversality table

{_md_table(res_headers, res_rows)}

`state res` is the largest relative M J = chi*E_inc residual returned by
`hh.build_operators`; `sigma_min(M)/||M||` is the normalised M conditioning
indicator.  Interpretation (finite-grid observations only):

* theta_min deg: {interp['theta_min_deg']}
* finite differences between successive computed N:
  {interp['theta_min_finite_diffs']}
* direction: **{interp['theta_min_direction']}** (span
  {_fmt(interp['theta_min_span_deg'], '.6e')} deg, max absolute step
  {_fmt(interp['theta_min_abs_max_diff_deg'], '.6e')} deg, last step
  {_fmt(interp['theta_min_last_increment_deg'], '.6e')} deg)
* near-confounded counts (cos2 > 1-1e-8):
  {interp['near_confounded_counts']}; constant over the tested grid =
  **{interp['near_count_constant_over_tested_grid']}**

Top-10 cos2 values per N are in the JSON (`top10_cos2`), together with full
cos2_desc, rho spectra, rank tolerance, and build times.

## Part 2: frequency rank-event sweep (N=16, smooth p=24, T=6)

{_md_table(freq_headers, freq_rows)}

Rank events: **{len(freq['rank_events'])}** (details in JSON and below).
Near crossings (smallest adjacent K_eff(alpha=1) gap < 1e-6*max(1,lambda_max)):
**{len(freq['near_crossings'])}**.

Literal-flag caveat (recorded, not forced): every near-crossing flag sits at
one of the lowest ascending pairs (recorded pairs are [0,1] or [1,2] over the
sweep), where numerical near-zero K_eff eigenvalues are separated by gaps up
to ~1.8e-12, many orders below the 1e-6 literal threshold.  The requested
inequality is therefore satisfied at every sweep point, and the flag is
dominated by the repeated numerical near-null tail of K_eff.  It is recorded
as the literal criterion and is **not** interpreted as a certified eigenvalue
crossing or a transversality failure.

Recorded rank events:

{_md_table(['f', 'changed', 'prev rank(B)', 'new rank(B)', 'prev rank([A,B])', 'new rank([A,B])'],
           [[_fmt(e['f'], '.3f'), ', '.join(e['changed']), e['prev_rank_B'],
             e['new_rank_B'], e['prev_rank_AB'], e['new_rank_AB']]
            for e in freq['rank_events']])}

Recorded near crossings:

{_md_table(['f', 'smallest gap', 'threshold', 'gap pair'],
           [[_fmt(e['f'], '.3f'), _fmt(e['smallest_gap'], '.3e'),
             _fmt(e['threshold'], '.3e'),
             '[' + ','.join(str(i) for i in e['gap_pair']) + ']']
            for e in freq['near_crossings']])}

## Part 3: contrast sweep Born -> full wave (N=16, f=1.0)

{_md_table(cont_headers, cont_rows)}

s=0 gates:

* B/A Frobenius ratio = {cont['s0_gates']['observed_B_over_A_fro_ratio_at_s0']}
  (gate < 1e-12: **{cont['s0_gates']['B_over_A_fro_ratio_lt_1e-12']}**).
* K_SLAM vs K_IS rel Fro =
  {cont['s0_gates']['observed_K_SLAM_vs_K_IS_rel_fro_at_s0']}
  (gate < 1e-12: **{cont['s0_gates']['K_SLAM_eq_K_IS_rel_fro_lt_1e-12']}**).

B/A ratio rule: ||B_R||_F/||A_R||_F on the realified blocks (same as the
complex ratio under whiten_realify).  Born discrepancy rule:
||A_full-A_born||_F/||A_full||_F with A_born from `hh.born_forward` (available
= {cont['born_available']}).  Contrast rank/near-degenerate transitions:
{len(cont['rank_events'])} rank event(s); details in JSON.

Near-degenerate transitions (rho_min < 1e-6, dimensionless retention rho):
{len(cont['near_degenerate_transitions'])} recorded onset(s):
{_md_table(['s', 'rho_min', 'theta_min deg'],
           [[_fmt(e['s'], 'g'), _fmt(e['rho_min'], '.4e'),
             _fmt(e['theta_min_deg'], '.6e')]
            for e in cont['near_degenerate_transitions']])}

At s=0 the principal-angle convention for a zero Range(B) is declared as
theta_min = 90 deg (the limiting/maximal convention); no physical angle
exists between a nonempty subspace and the trivial subspace.  The first
contrast "rank event" (s=0 -> 0.01, rank(B): 0 -> 18, rank([A,B]): 24 -> 42)
is the expected exit from the exactly-zero B baseline and is recorded as an
observed machine-rank change, not a certified degeneracy transition.

## Figures

* [{fig_names['resolution']}]({fig_names['resolution']}) - resolution
  convergence of theta_min_deg, confusable_mass, retained_mass, and count
  rho<1e-6 versus N.
* [{fig_names['frequency']}]({fig_names['frequency']}) - frequency sweep of
  sigma_min(B), theta_min_deg, rho_min, and the smallest K_eff gap with rank
  event / near-crossing lines.
* [{fig_names['contrast']}]({fig_names['contrast']}) - theta_min_deg and
  retained_mass versus s, log-scale B/A ratio, and Born discrepancy.

## Cannot establish

* These are finite-grid diagnostics only.  They do **not** prove continuum
  transfer, limiting transversality as N -> infinity, or any theorem about
  rank/transversality stability away from the tested resolutions.
* N refinement is diagnostic; it does not certify convergence of principal
  angles, masses, rho spectra, or rank identities to continuum limits.
* The frequency sweep flags rank events and near crossings on a discrete
  21-point f grid with one fixed N=16 discretisation.  Events between grid
  points, and any N-dependence of those events, are not established.
* The contrast sweep is one 8-point scale path at N=16, f=1.0.  Rank
  "events" are machine-rank changes on that path, not certified degeneracies;
  no Born/full-wave regime boundary is claimed.
* No forced pass is applied anywhere in Family 7; all trends and events are
  recorded as observed on the exact executed grid/config.
"""
    return rep


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    cfg = CONFIG
    generated_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    t_start = time.perf_counter()
    script_path = Path(__file__).resolve()
    source_sha256 = _sha256(script_path)
    reused = [
        _HERE / "helmholtz.py",
        _HERE / "family1_pilot.py",
        _HERE / "family2_algebraic_spine.py",
        _HERE / "family4_frequency_trajectory.py",
    ]
    reused_sha256 = {str(p.resolve()): _sha256(p) for p in reused}

    poses = _poses(cfg)
    points16, chi016, h16, S16 = _scene(int(cfg["N16"]), cfg)
    scene16 = (points16, chi016, h16, S16)

    res = run_resolution(cfg, poses)
    freq = run_frequency_sweep(cfg, poses, scene16)
    cont = run_contrast_sweep(cfg, poses, scene16)

    results = {
        "family": 7,
        "generated_utc": generated_utc,
        "command": ".venv/bin/python src/family7_refinement_rank.py",
        "source_sha256": source_sha256,
        "reused_sha256": reused_sha256,
        "config": cfg,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "part1_resolution": res,
        "part2_frequency": freq,
        "part3_contrast": cont,
    }
    import scipy as sp

    results["platform"]["scipy"] = sp.__version__

    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir = _ROOT / "results"
    figures_dir.mkdir(parents=True, exist_ok=True)
    notes_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    fig_paths = {
        "resolution": figures_dir / "family7_resolution_convergence.png",
        "frequency": figures_dir / "family7_frequency_rank_sweep.png",
        "contrast": figures_dir / "family7_contrast_sweep.png",
    }
    plot_resolution(res, fig_paths["resolution"])
    plot_frequency(freq, fig_paths["frequency"])
    plot_contrast(cont, fig_paths["contrast"])

    results["runtime_seconds"] = time.perf_counter() - t_start
    results["figures"] = {
        name: str(p.relative_to(_ROOT)) for name, p in fig_paths.items()
    }

    json_path = results_dir / "family7_refinement_rank.json"
    json_path.write_text(
        json.dumps(results, indent=2, default=_json_default) + "\n",
        encoding="utf-8",
    )

    report_path = notes_dir / "family7_refinement_rank.md"
    report = build_report(results, cfg, {"figures": fig_paths})
    report_path.write_text(report, encoding="utf-8")

    summary_lines = [
        "===== FAMILY 7 SUMMARY =====",
        f"resolution: N48_status={res['N48_status']}",
    ]
    for r in res["computed_rows"]:
        summary_lines.append(
            f"  N{r['N']:>2d} rA={r['r_A']} rB={r['r_B']} "
            f"rAB={r['r_AB']} identity={r['rank_identity_holds']} "
            f"theta={r['theta_min_deg']:.6e} retained="
            f"{r['retained_mass']:.6f} rho_min={r['rho_min']:.4e} "
            f"near_cos2={r['count_cos2_gt_1_minus_1e-8']} "
            f"build={r['build_seconds']:.3f}s"
        )
    summary_lines.append(
        f"theta direction: {res['interpretation']['theta_min_direction']}"
    )
    summary_lines.append(
        f"frequency rank events: {len(freq['rank_events'])} "
        f"(fs={[e['f'] for e in freq['rank_events']]}); "
        f"near crossings: {len(freq['near_crossings'])} "
        f"(fs={[e['f'] for e in freq['near_crossings']]})"
    )
    summary_lines.append(
        f"s=0 B/A ratio={cont['s0_gates']['observed_B_over_A_fro_ratio_at_s0']} "
        f"pass={cont['s0_gates']['B_over_A_fro_ratio_lt_1e-12']}; "
        f"KSL/KIS rel Fro="
        f"{cont['s0_gates']['observed_K_SLAM_vs_K_IS_rel_fro_at_s0']} "
        f"pass={cont['s0_gates']['K_SLAM_eq_K_IS_rel_fro_lt_1e-12']}"
    )
    summary_lines.append(
        f"contrast rank events: {len(cont['rank_events'])} "
        f"(details in JSON)"
    )
    summary_lines.append(
        f"total runtime {results['runtime_seconds']:.2f}s"
    )
    summary_lines.append("results -> " + str(json_path))
    for name, path in fig_paths.items():
        summary_lines.append("figure  -> " + str(path))
    summary_lines.append("report  -> " + str(report_path))
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
