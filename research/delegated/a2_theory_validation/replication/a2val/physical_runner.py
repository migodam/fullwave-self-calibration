"""Array helpers for the OPTIONAL physical second-stage tangent runner.

This module deliberately does NOT import or call the shared physical core
(``physics.py``): the script and the smoke test own that import and pass plain
real, unit-noise-whitened tangents here.

Semantics (kept explicit so no claim can leak from the algebraic checks):

* every result is a "physical tangent smoke, N=8, declared cutoffs, no oracle
  rank" record;
* the nested current spaces C0/C1/C2 are DECLARED left singular directions of
  the realified full data tangent, never a claim about physical current rank;
* numerical rank cutoffs are absolute and tied to ``||B|| * eps * n``, never a
  small relative-to-self ratio.
"""

from __future__ import annotations

import numpy as np

from a2val import common
from a2val import e2
from a2val import e3
from a2val import e5

EPS = float(np.finfo(float).eps)


def _fro(X) -> float:
    return float(np.linalg.norm(np.asarray(X, dtype=float), ord="fro"))


def _sym(X) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    return 0.5 * (X + X.T)


def declared_note() -> str:
    return "physical tangent smoke, N=8, declared cutoffs, no oracle rank"


# ---------------------------------------------------------------------------
# Numerical-rank controls (absolute, backward-error-scaled)
# ---------------------------------------------------------------------------


def backward_rank_tol(B: np.ndarray) -> float:
    """Absolute rank cutoff tied to ``10*eps*max(1,||B||_F)*n``.

    ``B`` is the original (unprojected) pose tangent of the same experiment;
    the cutoff is deliberately computed from that scale so it never becomes a
    tiny relative ratio of a shrunk residual matrix.
    """
    B = np.asarray(B, dtype=float)
    n = int(B.shape[0])
    return 10.0 * EPS * max(1.0, float(np.linalg.norm(B, ord="fro"))) * max(1, n)


def absolute_singular_values(X: np.ndarray):
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return np.zeros(0)
    return np.linalg.svd(X, compute_uv=False)


def numerical_rank_abs(X: np.ndarray, B: np.ndarray) -> int:
    """Rank of ``X`` under :func:`backward_rank_tol` (computed from ``B``)."""
    X = np.asarray(X, dtype=float)
    tol = backward_rank_tol(B)
    s = absolute_singular_values(X)
    if s.size == 0 or s[0] <= 0.0:
        return 0
    return int(np.count_nonzero(s > tol))


# ---------------------------------------------------------------------------
# E2 Theorem-5 helpers on physical tangents
# ---------------------------------------------------------------------------


def declared_E_basis(A, C_low, C_high) -> tuple[np.ndarray, dict]:
    """Robust basis of ``E = Ran([C_high,A]) cap Ran([C_low,A])^perp``.

    Mirrors ``a2val.e2.E_basis`` but replaces the purely relative projected-SVD
    drop with an absolute backward-error floor.  When the high and low nuisance
    ranges coincide up to projector roundoff (the physical C1->C2 case here),
    the naive relative floor retains machine-noise columns and would make T5a
    fail spuriously.  The floor is tied to unit-norm basis columns and ``n``:
    ``1000*eps*n``, so genuine small-angle directions whose information loss is
    below the backward-error scale are legitimately dropped.
    """
    A = np.asarray(A, dtype=float)
    C_low = np.asarray(C_low, dtype=float)
    C_high = np.asarray(C_high, dtype=float)
    n = A.shape[0]
    N_lo = np.hstack([C_low, A]) if C_low.shape[1] else A
    N_hi = np.hstack([C_high, A]) if C_high.shape[1] else A
    basis_hi = e2._col_basis(N_hi)
    if basis_hi.shape[1] == 0:
        return np.zeros((n, 0)), {"rank": 0}
    P_lo = common.orth_proj(N_lo) if N_lo.shape[1] else np.zeros((n, n))
    Y = (np.eye(n, dtype=float) - P_lo) @ basis_hi
    u, s, _ = np.linalg.svd(Y, full_matrices=False)
    floor = 1e3 * EPS * max(1.0, float(n))
    r = int(np.count_nonzero(s > floor))
    info = {
        "rank": int(r),
        "floor_absolute": float(floor),
        "raw_projected_singular_values": s.tolist(),
        "note": (
            "absolute backward-error floor on projected basis columns; "
            "relative-only 1e-12 floor of e2.E_basis is NOT used here because "
            "it retains machine-noise dimensions when the nuisance ranges "
            "coincide"
        ),
    }
    if r == 0:
        return np.zeros((n, 0)), info
    return u[:, :r].copy(), info


def t5_adjacent_pair(A, B, C_low, C_high, pair_label: str) -> dict:
    """One T5a/T5b check for an adjacent declared pair ``C_low -> C_high``.

    Reuses ``a2val.e2`` (``Jx``, ``Bv``, and the ``E_basis`` projection
    definition) for the core algebra.
    The scaled residual follows the existing E2 convention
    ``eps*max(1,||B||_F**2)*n``, and the rank cutoff is the absolute
    backward-error cutoff of :func:`backward_rank_tol`.  The ``E_r`` basis is
    obtained from :func:`declared_E_basis` (the same projection definition as
    ``e2.E_basis`` with a machine-noise floor, which is necessary for the
    physical coincidence case).
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    C_low = np.asarray(C_low, dtype=float)
    C_high = np.asarray(C_high, dtype=float)
    n = A.shape[0]

    E, E_info = declared_E_basis(A, C_low, C_high)
    P_E = E @ E.T if E.shape[1] else np.zeros((n, n), dtype=float)
    J_r = e2.Jx(A, B, C_low)
    J_r1 = e2.Jx(A, B, C_high)
    predicted_loss = B.T @ P_E @ B
    actual_loss = J_r - J_r1
    loss_residual = _fro(predicted_loss - actual_loss)
    scale = EPS * max(1.0, float(np.linalg.norm(B, ord="fro") ** 2)) * n
    scaled_residual = loss_residual / scale

    Bv_r = e2.Bv(A, B, C_low)
    Bv_r1 = e2.Bv(A, B, C_high)
    sv_r = absolute_singular_values(Bv_r)
    sv_r1 = absolute_singular_values(Bv_r1)
    tol = backward_rank_tol(B)
    rank_r = numerical_rank_abs(Bv_r, B)
    rank_r1 = numerical_rank_abs(Bv_r1, B)
    drop_actual = rank_r - rank_r1
    if E.shape[1]:
        inter = (
            rank_r
            + int(E.shape[1])
            - numerical_rank_abs(np.hstack([Bv_r, E]), B)
        )
    else:
        inter = 0

    return {
        "pair": pair_label,
        "declared_dims": {
            "C_low": int(C_low.shape[1]),
            "C_high": int(C_high.shape[1]),
            "E_dim": int(E.shape[1]),
        },
        "E_basis_floor": E_info,
        "loss_residual_fro": loss_residual,
        "backward_scaled_residual": float(scaled_residual),
        "scale": float(scale),
        "rank_bv_low": int(rank_r),
        "rank_bv_high": int(rank_r1),
        "rank_drop_actual": int(drop_actual),
        "rank_drop_predicted": int(inter),
        "t5b_matches": bool(inter == drop_actual),
        "bv_low_abs_singular_values": np.sort(sv_r)[::-1].tolist(),
        "bv_high_abs_singular_values": np.sort(sv_r1)[::-1].tolist(),
        "rank_tol_absolute": tol,
    }


# ---------------------------------------------------------------------------
# E3 helpers on physical tangents
# ---------------------------------------------------------------------------


def theorem3_physical(A, B, C, support_tol_rel: float = 1e-12) -> dict:
    """Theorem-3 canonical spectral identity for declared nuisance ``C``."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    C = np.asarray(C, dtype=float)
    grams = e3.map_pose_matrices(A, B, C)
    Ac, Bc = grams["A_c"], grams["B_c"]
    norm = e3.dual_normalized_spectra(Ac, Bc, support_tol_rel=support_tol_rel)
    info = e3.canonical_info(Ac, Bc)
    map_exp, pose_exp = e3.expected_dual_spectra(info)

    def _max_diff(vals, exp):
        if len(vals) != len(exp):
            return float("inf")
        if len(vals) == 0:
            return 0.0
        return float(np.max(np.abs(np.asarray(vals) - np.asarray(exp))))

    return {
        "a": info["a"],
        "b": info["b"],
        "canonical_c": info["c"].tolist(),
        "rank_Ccan": info["rank_Ccan"],
        "map_normalized_spectrum": norm["map"]["values"].tolist(),
        "pose_normalized_spectrum": norm["pose"]["values"].tolist(),
        "map_expected_spectrum": map_exp.tolist(),
        "pose_expected_spectrum": pose_exp.tolist(),
        "map_support_dim": norm["map"]["support_dim"],
        "pose_support_dim": norm["pose"]["support_dim"],
        "map_residual_max_abs": _max_diff(
            norm["map"]["values"], map_exp
        ),
        "pose_residual_max_abs": _max_diff(
            norm["pose"]["values"], pose_exp
        ),
        "raw_K0_eigenvalues": np.linalg.eigvalsh(_sym(grams["K0"])).tolist(),
        "raw_K_e_eigenvalues": np.linalg.eigvalsh(_sym(grams["K_e"])).tolist(),
        "raw_J_x_eigenvalues": np.linalg.eigvalsh(_sym(grams["J_x"])).tolist(),
    }


def prior_sandwich_physical(A, B, C, Lambda_x) -> dict:
    """Loewner sandwich ``K_e <= K_eL(Lambda_x) <= K0`` (pose-prior)."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    C = np.asarray(C, dtype=float)
    grams = e3.map_pose_matrices(A, B, C)
    Ac, Bc = grams["A_c"], grams["B_c"]
    K0, Ke = grams["K0"], grams["K_e"]
    K_eL = e3.K_eL(Ac, Bc, Lambda_x)
    min_low = float(np.min(np.linalg.eigvalsh(_sym(K_eL - Ke))))
    min_up = float(np.min(np.linalg.eigvalsh(_sym(K0 - K_eL))))
    return {
        "min_eig_K_eL_minus_K_e": min_low,
        "min_eig_K0_minus_K_eL": min_up,
        "loewner_sandwich_ok": bool(
            min_low >= -1e-8 * max(1.0, float(np.linalg.norm(K0, 2)))
            and min_up >= -1e-8 * max(1.0, float(np.linalg.norm(K0, 2)))
        ),
        "K_e_eigenvalues": np.linalg.eigvalsh(_sym(Ke)).tolist(),
        "K_eL_eigenvalues": np.linalg.eigvalsh(_sym(K_eL)).tolist(),
        "K0_eigenvalues": np.linalg.eigvalsh(_sym(K0)).tolist(),
    }


def theorem7_physical(A, B, C, rng: np.random.Generator) -> dict:
    """Theorem-7 analytic risk at declared nuisance ``C``.

    Bias directions are declared as the first left singular vector of the
    residual matrix ``R = (I-P_[C,A])B`` (a declared model-error direction,
    not an oracle), falling back to a fixed random unit vector from the same
    ``rng`` if ``R`` is numerically zero.  The pose error metric is declared
    as ``M_x = I_3`` (unit metric on the three pose coordinates).
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    C = np.asarray(C, dtype=float)
    n = A.shape[0]
    R = e2.Bv(A, B, C)  # residual matrix (I-P_[C,A]) B
    u, s, _ = np.linalg.svd(R, full_matrices=False)
    bias_source = "first_left_singular_of_residual"
    if s.size == 0 or s[0] <= max(1e-14, EPS * float(np.linalg.norm(R, "fro"))):
        u0 = rng.standard_normal(n)
        u0 = u0 / max(np.linalg.norm(u0), np.finfo(float).tiny)
        bias_source = "random_unit_fallback_residual_zero"
    else:
        u0 = u[:, 0]
    D_r = u0.reshape(-1, 1)
    fx = {
        "B_v": R,
        "D_r": D_r,
        "M_x": np.eye(3, dtype=float),
    }
    Jx = R.T @ R
    ev = np.linalg.eigvalsh(_sym(Jx))
    full_rank_Bv = bool(ev[0] > 1e-12 * max(1.0, float(ev[-1])))
    if not full_rank_Bv:
        return {
            "bias_direction_source": bias_source,
            "full_rank_B_v": False,
            "min_eig_J_x": float(ev[0]),
            "note": (
                "B_v is not full column rank, so the Theorem-7 inverse risk "
                "is not evaluated (rank-deficient physical smoke record)."
            ),
        }
    ana = e3.theorem7_analytic(fx)
    return {
        "bias_direction_source": bias_source,
        "full_rank_B_v": True,
        "min_eig_J_x": float(ev[0]),
        "variance_trace": ana["variance_trace"],
        "worst_bias_squared": ana["worst_bias_squared"],
        "analytic_total_risk": ana["analytic_total_risk"],
        "J_x_eigenvalues": np.linalg.eigvalsh(Jx).tolist(),
    }


# ---------------------------------------------------------------------------
# E5 helpers on physical frequency frames
# ---------------------------------------------------------------------------


def split_frequency_frames(A_r, B_r, blocks) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split realified full tangents into four per-frequency frames.

    ``A_r``/``B_r`` follow the realification layout ``[Re; Im]``, so for every
    complex row index ``r`` the two realified rows are ``r`` and ``r+half``.
    ``blocks`` supplies the complex-row intervals per frequency; within each
    frequency the complex rows are contiguous, which is asserted.
    """
    A_r = np.asarray(A_r, dtype=float)
    B_r = np.asarray(B_r, dtype=float)
    n_total = A_r.shape[0]
    half = n_total // 2
    if n_total != 2 * half or B_r.shape[0] != n_total:
        raise ValueError("realified tangents must have even row count and match")
    freq_idx = list(blocks["freq_idx"])
    row_start = list(blocks["row_start"])
    row_stop = list(blocks["row_stop"])
    frames = []
    for fi in sorted(set(freq_idx)):
        rows = []
        for f, rs, re in zip(freq_idx, row_start, row_stop):
            if f == fi:
                rows.extend(range(int(rs), int(re)))
        rows = sorted(set(rows))
        if not rows or rows != list(range(rows[0], rows[-1] + 1)):
            raise ValueError(f"frequency {fi} rows are not contiguous")
        real_rows = np.asarray(rows, dtype=int)
        idx = np.concatenate([real_rows, real_rows + half])
        frames.append((A_r[idx, :].copy(), B_r[idx, :].copy()))
    return frames


def stack_vs_sum_physical(frames) -> dict:
    """Per-frame pose info, stacked info, and the ``J_stack >= sum J_l`` check."""
    J_l = [e5._pose_info(a, b) for a, b in frames]
    A, B = e5._stack_AB(frames)
    J_stack = e5._pose_info(A, B)
    sumJ = sum(J_l)
    diff = _sym(J_stack - sumJ)
    min_eig = float(np.min(np.linalg.eigvalsh(diff)))
    scale = max(1.0, float(np.linalg.norm(B, ord="fro") ** 2))
    tol = max(1e-10, 20.0 * EPS * scale * B.shape[0])
    return {
        "n_frames": len(frames),
        "frames_shapes": [list(a.shape) for a, b in frames],
        "J_per_frame_fro": [float(np.linalg.norm(J, "fro")) for J in J_l],
        "sum_J_fro": float(np.linalg.norm(sumJ, "fro")),
        "J_stack_fro": float(np.linalg.norm(J_stack, "fro")),
        "min_eig_J_stack_minus_sum": min_eig,
        "psd_tol_absolute": tol,
        "stack_ge_sum": bool(min_eig >= -tol),
        "J_stack_eigenvalues": np.linalg.eigvalsh(_sym(J_stack)).tolist(),
    }


def sequential_innovation(frames) -> dict:
    """Sequential shared-map innovation for frames 0..3.

    Uses the Theorem-9 formula from ``a2val.e5`` whenever the accumulated old
    map Gram is positive definite (min eig > 1e-10), otherwise falls back to
    the exact variational/SVD stacked criterion and marks the step singular.
    """
    A_old, B_old = frames[0]
    records = []
    singular_count = 0
    for step, (a, b) in enumerate(frames[1:], start=1):
        J_old = e5._pose_info(A_old, B_old)
        G = A_old.T @ A_old
        g_min = e5._min_eig(G)
        full_rank = bool(g_min > 1e-10)
        if full_rank:
            th = e5._theorem9_innovation(A_old, B_old, a, b)
            I_acq = th["I_acq"]
            J_formula = _sym(J_old + I_acq)
            method = "theorem9_inverse_formula"
        else:
            var = e5._variational_innovation(A_old, B_old, a, b)
            I_acq = var["I_acq"]
            J_formula = var["J_new"]
            method = "svd_qr_variational_fallback"
            singular_count += 1
        A_new = np.vstack([A_old, a])
        B_new = np.vstack([B_old, b])
        J_direct = e5._pose_info(A_new, B_new)
        residual = _fro(J_formula - J_direct)
        records.append({
            "step": int(step),
            "method": method,
            "G_min_eig": float(g_min),
            "G_full_rank": bool(full_rank),
            "singular": bool(not full_rank),
            "formula_direct_fro_residual": residual,
            "formula_direct_rel_residual": residual
            / (1.0 + float(np.linalg.norm(J_direct, "fro"))),
            "J_old_fro": float(np.linalg.norm(J_old, "fro")),
            "I_acq_fro": float(np.linalg.norm(I_acq, "fro")),
            "J_new_fro": float(np.linalg.norm(J_direct, "fro")),
            "J_new_min_eig": e5._min_eig(J_direct),
        })
        A_old, B_old = A_new, B_new
    return {
        "records": records,
        "singular_steps": int(singular_count),
        "n_steps": int(len(records)),
    }


def rank_budget_physical(A, B, a, b, seed: int, label: str) -> dict:
    """Corollary-10 rank budget for the physical transition old -> old + frame.

    Declared (no-oracle) enlargement: ``C_r`` = first 4 and ``C_r1`` = first 8
    left singular directions of the old full data tangent ``T=[A,B]``.  Reuses
    the E5 primitives; when the enlarged map Gram is singular the innovation is
    evaluated with the exact variational/SVD stacked fallback instead of the
    inverse formula (marked ``singular_fallback``).
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = A.shape[0]
    T = np.hstack([A, B])
    U, sT, _ = np.linalg.svd(T, full_matrices=False)
    C_r = U[:, : min(4, U.shape[1])]
    C_r1 = U[:, : min(8, U.shape[1])]

    def _proj(M):
        if M.shape[1] == 0:
            return np.zeros((n, n), dtype=float)
        return common.orth_proj(M)

    I_n = np.eye(n, dtype=float)
    P_Nr = _proj(np.hstack([C_r, A]))
    P_Nr1 = _proj(np.hstack([C_r1, A]))
    J_original = _sym(B.T @ (I_n - P_Nr) @ B)
    J_after_rank = _sym(B.T @ (I_n - P_Nr1) @ B)
    L_rank = _sym(J_original - J_after_rank)

    Y = (I_n - P_Nr) @ C_r1
    U_E = e5._col_basis(Y)
    P_E = U_E @ U_E.T if U_E.shape[1] else np.zeros((n, n))
    L_rank_E = _sym(B.T @ P_E @ B)

    P_Cr1 = _proj(C_r1)
    A_enl = (I_n - P_Cr1) @ A
    B_enl = (I_n - P_Cr1) @ B
    J_after_enlarged = e5._pose_info(A_enl, B_enl)
    G_enl = A_enl.T @ A_enl
    g_min = e5._min_eig(G_enl)
    singular_fallback = bool(g_min <= 1e-10)
    if singular_fallback:
        var = e5._variational_innovation(A_enl, B_enl, a, b)
        I_acq = var["I_acq"]
        J_final_formula = var["J_new"]
        method = "svd_qr_variational_fallback"
    else:
        th = e5._theorem9_innovation(A_enl, B_enl, a, b)
        I_acq = th["I_acq"]
        J_final_formula = _sym(J_after_enlarged + I_acq)
        method = "theorem9_inverse_formula"

    Ast = np.vstack([A_enl, a])
    Bst = np.vstack([B_enl, b])
    J_final_direct = e5._pose_info(Ast, Bst)
    identity_fro = _fro((J_final_direct - J_original) - (I_acq - L_rank))
    min_eig_delta = e5._min_eig(J_final_direct - J_original)
    min_eig_IL = e5._min_eig(I_acq - L_rank)
    tol_psd = max(1e-10, 20.0 * EPS * max(1.0, float(np.linalg.norm(B, "fro") ** 2)) * n)
    return {
        "seed": int(seed),
        "label": label,
        "method": method,
        "singular_fallback": bool(singular_fallback),
        "declared_dims": {"C_r": int(C_r.shape[1]), "C_r1": int(C_r1.shape[1])},
        "old_T_singular_values_top": sT.tolist(),
        "J_original_eigenvalues": np.linalg.eigvalsh(_sym(J_original)).tolist(),
        "J_after_rank_eigenvalues": np.linalg.eigvalsh(_sym(J_after_rank)).tolist(),
        "L_rank_fro": float(np.linalg.norm(L_rank, "fro")),
        "L_rank_E_residual_fro": _fro(L_rank - L_rank_E),
        "E_dim": int(U_E.shape[1]),
        "J_after_enlarged_vs_after_rank_residual": _fro(
            J_after_enlarged - J_after_rank
        ),
        "G_enlarged_min_eig": float(g_min),
        "G_enlarged_full_rank": bool(not singular_fallback),
        "I_acq_fro": float(np.linalg.norm(I_acq, "fro")),
        "identity_fro_residual": identity_fro,
        "identity_rel_residual": identity_fro
        / (1.0 + float(np.linalg.norm(J_final_direct - J_original, "fro"))),
        "min_eig_J_final_minus_J_original": min_eig_delta,
        "min_eig_I_acq_minus_L_rank": min_eig_IL,
        "loewner_I_ge_L": bool(min_eig_IL >= -tol_psd),
        "psd_tol_absolute": tol_psd,
    }
