"""Experiment E2: nested free-current rank loss, neutral admission, and
numerical-rank controls.

All matrices are real, unit-noise whitened, gauge-fixed tangent matrices:

* ``A``  : n x q map tangent (may be empty, then ``n x 0``),
* ``B``  : n x p pose tangent,
* ``C_r``: n x r_r free-current nuisance basis, nested ``C_r subset C_{r+1}``.

The current-and-map nuisance space is ``N_r = Ran([C_r, A])`` and

    P_Nr     = orth_proj([C_r, A])
    J_x(r)   = B^T (I - P_Nr) B
    B_v(r)   = (I - P_Nr) B
    E_r      = N_{r+1} cap N_r^perp

The module verifies Theorem 5

    J_x(r) - J_x(r+1) = B^T P_{E_r} B
    rank(B_v(r)) - rank(B_v(r+1)) = dim(Ran(B_v(r)) cap E_r),

the nonmonotone relative map-retention spectrum, calibration-neutral
admission (Theorem 6, including the complex-linear ``J_c``-invariant kernel),
and the numerical-rank / projector-gap controls.
"""

from __future__ import annotations

import numpy as np

from a2val import common

EPS = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Core linear-algebra helpers
# ---------------------------------------------------------------------------


def _mat(X, name="X"):
    X = np.asarray(X, dtype=float)
    if X.ndim == 0:
        raise ValueError(f"{name} must be 1-D or 2-D")
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    return X


def _nuisance_cols(A, C):
    """Column stack ``[C | A]`` forming the nuisance basis of N_r."""
    A = _mat(A, "A")
    C = _mat(C, "C")
    if A.shape[0] != C.shape[0]:
        raise ValueError("A and C must have the same number of rows")
    if C.shape[1] == 0 and A.shape[1] == 0:
        return np.zeros((A.shape[0], 0))
    return np.hstack([C, A])


def _col_basis(X, rcond=None):
    """Orthonormal basis of the column range of ``X`` via SVD."""
    X = _mat(X, "X")
    n = X.shape[0]
    if X.shape[1] == 0:
        return np.zeros((n, 0))
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    if s.size == 0 or s[0] <= 0.0:
        return np.zeros((n, 0))
    if rcond is None:
        rcond = max(X.shape) * EPS
    r = int(np.count_nonzero(s > rcond * s[0]))
    if r == 0:
        return np.zeros((n, 0))
    return u[:, :r].copy()


def _complement_projector(A, C=None):
    """``I - orth_proj([C, A])`` as an n x n real matrix."""
    n = _mat(A).shape[0]
    cols = _nuisance_cols(A, np.zeros((n, 0)) if C is None else C)
    P = common.orth_proj(cols)
    return np.eye(n, dtype=P.dtype) - P


def Jx(A, B, C):
    """Pose information after removing the current+map nuisance: ``B^T R B``."""
    B = _mat(B, "B")
    R = _complement_projector(A, C)
    return B.T @ R @ B


def Bv(A, B, C):
    """Residual visible pose tangent ``(I - P_{[C,A]}) B``."""
    B = _mat(B, "B")
    R = _complement_projector(A, C)
    return R @ B


def E_basis(A, C_low, C_high, abs_rel=1e-12):
    """Orthonormal basis of ``E = N_high cap N_low^perp``.

    Uses the definition from the validation protocol: take an orthonormal
    basis of ``N_high``, project it onto ``N_low^perp``, and take an
    SVD-based orthonormal basis of the projected image, dropping projected
    singular values <= ``abs_rel`` times the largest projected value.
    """
    A = _mat(A, "A")
    C_low = _mat(C_low, "C_low")
    C_high = _mat(C_high, "C_high")
    n = A.shape[0]
    N_hi = _nuisance_cols(A, C_high)
    if N_hi.shape[1] == 0:
        return np.zeros((n, 0))
    basis_hi = _col_basis(N_hi)
    P_lo = common.orth_proj(_nuisance_cols(A, C_low))
    Y = (np.eye(n, dtype=float) - P_lo) @ basis_hi
    u, s, _ = np.linalg.svd(Y, full_matrices=False)
    if s.size == 0 or s[0] <= 0.0:
        return np.zeros((n, 0))
    r = int(np.count_nonzero(s > abs_rel * s[0]))
    if r == 0:
        return np.zeros((n, 0))
    return u[:, :r].copy()


def K0_map(A, C):
    """Map Gram on the complement of the current space: ``A^T (I-P_C) A``."""
    A = _mat(A, "A")
    if A.shape[1] == 0:
        return np.zeros((A.shape[0], 0))
    C = _mat(C, "C")
    P_C = common.orth_proj(C) if C.shape[1] else np.zeros((A.shape[0], A.shape[0]))
    return A.T @ (np.eye(A.shape[0]) - P_C) @ A


def K_eff_map(A, B, C):
    """Map Gram after also removing pose: ``A^T (I - P_{[C,B]}) A``."""
    A = _mat(A, "A")
    if A.shape[1] == 0:
        return np.zeros((A.shape[0], 0))
    B = _mat(B, "B")
    C = _mat(C, "C")
    P_CB = common.orth_proj(_nuisance_cols(B, C))
    return A.T @ (np.eye(A.shape[0]) - P_CB) @ A


def rho_values(A, B, C, support_tol_rel=1e-12):
    """Sorted relative map-retention values ``rho(r)`` on the support of K0.

    ``K0 = A^T (I-P_C) A``, ``K_eff = A^T (I-P_{[C,B]}) A``.  The support is
    the range of ``K0`` with singular values above ``support_tol_rel`` times
    the largest one; generalized eigenvalues are computed from
    ``K0^{-1/2} K_eff K0^{-1/2}`` restricted to that support.  Directions in
    which ``K0 = 0`` are excluded (they are undefined, not assigned 0 or 1).
    """
    A = _mat(A, "A")
    B = _mat(B, "B")
    C = _mat(C, "C")
    n = A.shape[0]
    if A.shape[1] == 0:
        return {
            "rho": [],
            "support_dim": 0,
            "q": 0,
            "k0_singular_values": [],
            "support_tol": None,
            "k0_rank_all": 0,
            "k_eff": [],
        }
    K0 = K0_map(A, C)
    Kef = K_eff_map(A, B, C)
    sK0 = np.linalg.svd(K0, compute_uv=False)
    smax = float(sK0[0]) if sK0.size else 0.0
    if smax <= 0.0:
        return {
            "rho": [],
            "support_dim": 0,
            "q": int(A.shape[1]),
            "k0_singular_values": sK0.tolist(),
            "support_tol": 0.0,
            "k0_rank_all": 0,
            "k_eff": Kef.tolist(),
        }
    tol = support_tol_rel * smax
    w, U = np.linalg.eigh(K0)
    support = w > tol
    r_sup = int(np.count_nonzero(support))
    if r_sup == 0:
        return {
            "rho": [],
            "support_dim": 0,
            "q": int(A.shape[1]),
            "k0_singular_values": sK0.tolist(),
            "support_tol": float(tol),
            "k0_rank_all": int(np.count_nonzero(sK0 > tol)),
            "k_eff": Kef.tolist(),
        }
    Ur = U[:, support]
    K0r = Ur.T @ K0 @ Ur
    Kefr = Ur.T @ Kef @ Ur
    w0, V0 = np.linalg.eigh(K0r)
    wpos = np.clip(w0, 0.0, None)
    half_inv = np.diag(1.0 / np.sqrt(wpos))
    K0inv_half = (V0 @ half_inv) @ V0.T
    M = (K0inv_half @ Kefr) @ K0inv_half
    vals = np.linalg.eigvalsh(M)
    return {
        "rho": np.sort(vals).tolist(),
        "support_dim": int(r_sup),
        "q": int(A.shape[1]),
        "k0_singular_values": sK0.tolist(),
        "support_tol": float(tol),
        "k0_rank_all": int(np.count_nonzero(sK0 > tol)),
        "k_eff": Kef.tolist(),
    }


def _fro(X):
    return float(np.linalg.norm(np.asarray(X, dtype=float), "fro"))


def _eig_diff_psd(L):
    """Sorted descending eigenvalues of the symmetric matrix ``L``."""
    w = np.linalg.eigvalsh((np.asarray(L, dtype=float) + np.asarray(L, dtype=float).T) / 2.0)
    return np.sort(w)[::-1]


def _rank(X, rcond=None):
    """Numerical rank with a backward-error-scaled cutoff.

    Exact-rank identities on projected residual matrices are only meaningful
    above the roundoff left by computing ``(I - P) B``, so the default cutoff
    is ``10*eps*max(1, sigma_max)*max(shape)`` (an absolute floor scaled to
    the matrix, not a pure sigma_i/sigma_1 ratio).  ``rcond`` overrides this
    with the classic relative convention when explicitly passed.
    """
    X = _mat(X, "X")
    s = np.linalg.svd(X, compute_uv=False)
    if s.size == 0 or s[0] <= 0.0:
        return 0
    if rcond is not None:
        return int(np.count_nonzero(s > rcond * s[0]))
    tol = 10.0 * EPS * max(1.0, float(s[0])) * max(X.shape)
    return int(np.count_nonzero(s > tol))


# ---------------------------------------------------------------------------
# Exact small controls
# ---------------------------------------------------------------------------


def _e1(A, B, C0, C1, name, expected=None):
    """Run one adjacent-pair exact control and return a full result dict."""
    n = _mat(A).shape[0]
    A0 = _mat(A)
    B0 = _mat(B)
    C0m = _mat(C0)
    C1m = _mat(C1)
    E = E_basis(A0, C0m, C1m)
    P_E = E @ E.T if E.shape[1] else np.zeros((n, n))
    J0 = Jx(A0, B0, C0m)
    J1 = Jx(A0, B0, C1m)
    pred_loss = B0.T @ P_E @ B0
    actual_loss = J0 - J1
    loss_residual = _fro(pred_loss - actual_loss)
    Bv0 = Bv(A0, B0, C0m)
    Bv1 = Bv(A0, B0, C1m)
    r0 = _rank(Bv0)
    r1 = _rank(Bv1)
    if E.shape[1]:
        inter = r0 + E.shape[1] - _rank(np.hstack([Bv0, E]))
    else:
        inter = 0
    drop_actual = r0 - r1
    drop_pred = inter
    t5a = bool(loss_residual <= 1e-12)
    t5b = bool(drop_pred == drop_actual)
    rec = {
        "case": name,
        "n": int(n),
        "A": A0.tolist(),
        "B": B0.tolist(),
        "C0": C0m.tolist(),
        "C1": C1m.tolist(),
        "E": E.tolist(),
        "j0": J0.tolist(),
        "j1": J1.tolist(),
        "actual_loss": actual_loss.tolist(),
        "predicted_loss": pred_loss.tolist(),
        "loss_residual_fro": loss_residual,
        "rank_bv0": int(r0),
        "rank_bv1": int(r1),
        "dim_E": int(E.shape[1]),
        "rank_drop_actual": int(drop_actual),
        "rank_drop_predicted": int(drop_pred),
        "intersection_dim": int(inter),
        "t5a": t5a,
        "t5b": t5b,
    }
    if expected is not None:
        rec["expected"] = expected
    return rec


def exact_controls():
    """Exact small controls (A1-A9)."""
    out = {}

    # A1 equality case: current direction orthogonal to visible pose.
    n3 = 3
    A = np.zeros((n3, 0))
    B = np.eye(n3)[:, :2]
    C0 = np.zeros((n3, 0))
    C1 = np.eye(n3)[:, 2:3]
    rec = _e1(A, B, C0, C1, "equality_case",
              expected={"j0": "I2", "j1": "I2", "loss": 0, "rank_drop": 0})
    rec["j0_expect"] = [[1.0, 0.0], [0.0, 1.0]]
    rec["j1_expect"] = [[1.0, 0.0], [0.0, 1.0]]
    out["equality_case"] = rec

    # A2 strict loss: current direction aligned with one pose direction.
    C1 = np.eye(n3)[:, 0:1]
    rec = _e1(A, B, C0, C1, "strict_loss_case",
              expected={"j0": "I2", "j1": "diag(0,1)",
                        "loss": "e1 e1^T", "rank_drop": 1})
    rec["j0_expect"] = [[1.0, 0.0], [0.0, 1.0]]
    rec["j1_expect"] = [[0.0, 0.0], [0.0, 1.0]]
    out["strict_loss_case"] = rec

    # A3 complete hiding.
    C1 = np.eye(n3)[:, :2]
    rec = _e1(A, B, C0, C1, "complete_hiding_case",
              expected={"j0": "I2", "j1": 0, "loss": "I2", "rank_drop": 2})
    rec["j0_expect"] = [[1.0, 0.0], [0.0, 1.0]]
    rec["j1_expect"] = [[0.0, 0.0], [0.0, 0.0]]
    out["complete_hiding_case"] = rec

    # A4 information loss without rank drop: B=(1,1)^T, C1=span(e1).
    n2 = 2
    A2 = np.zeros((n2, 0))
    B4 = np.array([[1.0], [1.0]])
    C0_4 = np.zeros((n2, 0))
    C1_4 = np.eye(n2)[:, 0:1]
    rec = _e1(A2, B4, C0_4, C1_4, "loss_without_rank_drop",
              expected={"j0": 2.0, "j1": 1.0, "rank_drop": 0})
    rec["j0_scalar"] = float(Jx(A2, B4, C0_4).item())
    rec["j1_scalar"] = float(Jx(A2, B4, C1_4).item())
    out["loss_without_rank_drop"] = rec

    # A5 nonmonotone rho (canonical falsification of monotone retention).
    A5 = np.array([[1.0], [1.0], [0.0]])
    B5 = np.array([[1.0], [0.0], [1.0]])
    C0_5 = np.zeros((3, 0))
    C1_5 = np.eye(3)[:, 0:1]
    C2_5 = np.column_stack([np.eye(3)[:, 0], np.array([0.0, 1.0, 1.0])])
    steps = [
        ("C0", C0_5, [0.75]),
        ("C1", C1_5, [1.0]),
        ("C2", C2_5, [0.0]),
    ]
    rho_records = []
    ok = True
    for label, Cstep, expect in steps:
        rec = rho_values(A5, B5, Cstep)
        rho = rec["rho"]
        good = len(rho) == 1 and abs(float(rho[0]) - expect[0]) <= 1e-12
        ok = ok and good
        rho_records.append({
            "chain": label,
            "C": Cstep.tolist(),
            "rho": rho,
            "expected": expect,
            "support_dim": rec["support_dim"],
            "k0_singular_values": rec["k0_singular_values"],
        })
    out["nonmonotone_rho"] = {
        "case": "nonmonotone_rho",
        "n": 3,
        "A": A5.tolist(),
        "B": B5.tolist(),
        "expected_rho": [0.75, 1.0, 0.0],
        "rho_records": rho_records,
        "rho_sequence": [r["rho"][0] for r in rho_records],
        "rho_monotone_falsified": bool(
            not (rho_records[0]["rho"] <= rho_records[1]["rho"] <= rho_records[2]["rho"]) or
            rho_records[2]["rho"][0] < rho_records[1]["rho"][0]
        ),
        "pass": ok,
    }

    # A6 rank uncertainty: N(t) = t e1, B = e1, n = 1.
    tvals = [0.0, 1e-14, 1e-10, 1e-6, 0.1]
    A6 = np.zeros((1, 0))
    B6 = np.ones((1, 1))
    uncertainty = []
    for t in tvals:
        N = np.array([[t]])
        j = float(Jx(A6, B6, N).item())
        s = float(np.linalg.svd(N, compute_uv=False)[0]) if N.shape[1] else 0.0
        # Default numpy-style relative rank (scale-invariant), plus absolute
        # context: the model's true rank at t is not oracle-certifiable.
        rel_rank = int(common.numerical_rank(N))
        uncertainty.append({
            "t": float(t),
            "j_x": j,
            "singular_value": s,
            "default_relative_rank": int(rel_rank),
        })
    j_zero = uncertainty[0]["j_x"]
    j_nonzero = [u["j_x"] for u in uncertainty[1:]]
    out["rank_uncertainty"] = {
        "case": "rank_uncertainty",
        "n": 1,
        "A": A6.tolist(),
        "B": B6.tolist(),
        "t_values": tvals,
        "records": uncertainty,
        "j_x_0": j_zero,
        "j_x_nonzero": j_nonzero,
        "discontinuity_confirmed": bool(
            abs(j_zero - 1.0) < 1e-12 and all(abs(x) < 1e-12 for x in j_nonzero)
        ),
        "numerical_rank_ambiguity": (
            "t=0 is exactly distinguishable only in exact arithmetic; for any "
            "tiny-but-nonzero t at or below roundoff of ||N|| no finite "
            "precision computation can certify the true rank without an oracle. "
            "A relative-to-itself SVD rank always reports rank 1 for a nonzero "
            "column, so record absolute singular values and the J discontinuity."
        ),
    }

    # A7 saturation control: C = full orthonormal basis of R^5.
    rng = np.random.default_rng(777)
    n5 = 5
    A7 = _scaled_orthonormal(rng, n5, 2)
    B7 = _scaled_orthonormal(rng, n5, 3)
    C7, _ = np.linalg.qr(rng.standard_normal((n5, n5)))
    R = _complement_projector(A7, C7)
    j7 = Jx(A7, B7, C7)
    Bv7 = Bv(A7, B7, C7)
    sv7 = np.linalg.svd(Bv7, compute_uv=False)
    scale = EPS * max(1.0, float(np.linalg.norm(B7, "fro") ** 2)) * n5
    out["saturation_control"] = {
        "case": "saturation_control",
        "n": n5,
        "seed": 777,
        "A": A7.tolist(),
        "B": B7.tolist(),
        "C_rank": int(common.numerical_rank(C7)),
        "j_x_fro": float(np.linalg.norm(j7, "fro")),
        "residual_R_fro": float(np.linalg.norm(R, "fro")),
        "bv_absolute_singular_values": sv7.tolist(),
        "bv_max_abs_singular_value": float(sv7[0]) if sv7.size else 0.0,
        "backward_scaled_residual": float(np.linalg.norm(j7, "fro") / scale),
        "backward_scale": float(scale),
        "note": "J_x is zero up to backward error; absolute singular values "
                "are recorded because a relative-to-self rank is meaningless "
                "for an all-roundoff residual matrix.",
    }

    # A8 rank-threshold control.
    dvals = np.array([1.0, 1e-8, 1e-14, 0.0])
    D = np.diag(dvals)
    sD = np.linalg.svd(D, compute_uv=False)
    tol_abs = 1e-12 * max(1, D.shape[0]) * float(sD[0])  # spectral norm of D
    rank_abs = int(np.count_nonzero(sD > tol_abs))
    tol_rel = 1e-12 * float(sD[0])
    rank_rel = int(np.count_nonzero(sD > tol_rel))
    ratios = [float(sD[i] / sD[0]) for i in range(sD.size) if sD[0] > 0]
    out["rank_threshold_control"] = {
        "case": "rank_threshold_control",
        "D_diag": dvals.tolist(),
        "singular_values": sD.tolist(),
        "tol_abs_formula": "1e-12 * max(1, n) * ||D||_2",
        "tol_abs": float(tol_abs),
        "rank_absolute": int(rank_abs),
        "tol_rel_formula": "sigma_i <= 1e-12 * sigma_1 (counted when sigma_i > tol)",
        "tol_rel_rank_rule": "rank = count(sigma_i > 1e-12 * sigma_1)",
        "tol_rel": float(tol_rel),
        "rank_relative": int(rank_rel),
        "ratios_to_sigma1": ratios,
        "small_ratio_trap": (
            "sigma_3/sigma_1 = 1e-14 is roundoff-level and must not be counted "
            "as physical rank; both cutoffs give rank 2.  A cutoff anchored to "
            "the previous singular value (sigma_3/sigma_2 = 1e-6 > 1e-12) or a "
            "naive eps-level absolute threshold would wrongly count sigma_3."
        ),
    }

    # A9 projector gap closure.
    gap_records = []
    t_rank = [0.5, 0.05, 1e-8, -1e-8, -0.05]
    for t in t_rank:
        eig = sorted([1.0 + t, 1.0 - t], reverse=True)
        P = np.array([[1.0, 0.0], [0.0, 0.0]]) if t > 0 else np.array(
            [[0.0, 0.0], [0.0, 1.0]]
        )
        gap_records.append({
            "t": float(t),
            "eigenvalues": eig,
            "spectral_gap": float(abs(2.0 * t)),
            "rank_top_one_projector": int(np.linalg.matrix_rank(P)),
            "P": P.tolist(),
        })
    P_pos = np.diag([1.0, 0.0])
    P_neg = np.diag([0.0, 1.0])
    d_pos_neg = float(np.linalg.norm(P_pos - P_neg, "fro"))
    out["projector_gap_closure"] = {
        "case": "projector_gap_closure",
        "H_t": "diag(1+t, 1-t)",
        "records": gap_records,
        "fro_P_plus_minus_0_05": float(np.sqrt(2.0)),
        "left_right_limits": {
            "P_plus_limit": P_pos.tolist(),
            "P_minus_limit": P_neg.tolist(),
            "fro_distance": d_pos_neg,
        },
        "distance_t0": 0.0,
        "limit_nonunique": (
            "At t=0 the eigenvalues are degenerate (I), so the top-one spectral "
            "projector has no unique limit: t->0+ gives diag(1,0) and t->0- "
            "gives diag(0,1); any unit vector defines an admissible rank-one "
            "choice.  distance(P(0),P(0)) is trivially 0 for a fixed choice."
        ),
        "rank_change_event_distinct": (
            "rank(H(t)) = 2 for every t, so closing the spectral gap at t=0 is "
            "NOT a rank-change event of H; it is a discontinuity/non-uniqueness "
            "event of the rank-one spectral projector (coordinate instability)."
        ),
    }
    return out


def _scaled_orthonormal(rng, n, k, lo=0.5, hi=2.0):
    """Random n x k matrix with orthonormal columns times sigmas in [lo, hi]."""
    if k == 0:
        return np.zeros((n, 0))
    q, _ = np.linalg.qr(rng.standard_normal((n, k)))
    s = rng.uniform(lo, hi, size=k)
    return q * s


# ---------------------------------------------------------------------------
# Neutral admission
# ---------------------------------------------------------------------------


def neutral_fixture(seed=2600):
    """Deterministic real fixture shared by the neutral-admission cases."""
    rng = np.random.default_rng(seed)
    n, q, p = 8, 2, 3
    A = _scaled_orthonormal(rng, n, q)
    B = _scaled_orthonormal(rng, n, p)
    C0 = _scaled_orthonormal(rng, n, 2, lo=0.8, hi=1.5)
    Qreal = _scaled_orthonormal(rng, n, 4, lo=0.8, hi=1.5)
    return {
        "seed": int(seed),
        "n": n,
        "q": q,
        "p": p,
        "A": A,
        "B": B,
        "C0": C0,
        "Qreal": Qreal,
    }


def _existing_nuisance_basis(A, C0):
    return _col_basis(np.hstack([C0, A]))


def _visible_pose(A, B, C0):
    N = _existing_nuisance_basis(A, C0)
    P_N = N @ N.T if N.shape[1] else np.zeros((A.shape[0], A.shape[0]))
    return (np.eye(A.shape[0]) - P_N) @ B


def _null_basis(X, rcond=1e-12):
    """Orthonormal basis of the (right) null space of ``X`` via SVD."""
    X = _mat(X)
    if X.shape[0] == 0:
        return np.eye(X.shape[1])
    u, s, vh = np.linalg.svd(X, full_matrices=True)
    if s.size == 0 or s[0] <= 0:
        return vh.T
    tol = rcond * s[0] if rcond is not None else max(X.shape) * EPS * s[0]
    r = int(np.count_nonzero(s > tol))
    return vh.T[:, r:].copy()


def _add_current(A, B, C0, extra):
    """J_x after appending one/more current columns ``extra`` to C0."""
    extra = _mat(extra)
    if extra.shape[1] == 0:
        return Jx(A, B, C0)
    return Jx(A, B, np.hstack([C0, extra]))


def real_safe_space():
    """Real calibration-neutral admission (Theorem 6, real case)."""
    fx = neutral_fixture()
    A, B, C0, Q = fx["A"], fx["B"], fx["C0"], fx["Qreal"]
    V = _visible_pose(A, B, C0)
    F = V.T @ Q
    Z_safe = _null_basis(F)
    dim_kerF = int(Z_safe.shape[1])
    rF = _rank(F, rcond=1e-12)
    d = int(Q.shape[1])
    p = int(B.shape[1])
    J_before = Jx(A, B, C0)

    safe_losses = []
    z_safe = Z_safe[:, 0:1]
    J_after_safe = _add_current(A, B, C0, Q @ z_safe)
    diff_safe = J_before - J_after_safe
    scale = EPS * max(1.0, float(np.linalg.norm(B, "fro") ** 2)) * A.shape[0]
    safe_losses.append({
        "fro_loss": float(np.linalg.norm(diff_safe, "fro")),
        "backward_scaled": float(np.linalg.norm(diff_safe, "fro") / scale),
    })

    # Generic unsafe direction: principal right singular vector of F.
    _, sF, vhF = np.linalg.svd(F, full_matrices=True)
    z_bad = vhF.T[:, 0:1]
    J_after_unsafe = _add_current(A, B, C0, Q @ z_bad)
    diff_unsafe = J_before - J_after_unsafe
    w_unsafe = _eig_diff_psd(diff_unsafe)
    loewner_ok = bool(
        float(np.min(w_unsafe)) >= -1e-10 and float(np.max(w_unsafe)) > 1e-10
    )
    J_safe_eigs = _eig_diff_psd(J_after_safe)
    J_before_eigs = _eig_diff_psd(J_before)
    J_unsafe_eigs = _eig_diff_psd(J_after_unsafe)
    return {
        "case": "real_safe_space",
        "seed": fx["seed"],
        "dims": {"n": fx["n"], "q": fx["q"], "p": fx["p"],
                 "d_candidate_cols": d},
        "existing_nuisance": "N = Ran([C0, A])  (current basis plus map tangent)",
        "nuisance_dim": int(_existing_nuisance_basis(A, C0).shape[1]),
        "F_shape": list(F.shape),
        "F_singular_values": sF.tolist(),
        "rank_F": int(rF),
        "dim_kerF": dim_kerF,
        "bound": {"d_minus_rankF": d - rF, "d_minus_p": d - p,
                  "holds": bool(d - rF >= d - p)},
        "j_before": J_before.tolist(),
        "j_before_eigs": J_before_eigs.tolist(),
        "safe_z": z_safe.tolist(),
        "j_after_safe": J_after_safe.tolist(),
        "j_after_safe_eigs": J_safe_eigs.tolist(),
        "safe_losses": safe_losses,
        "safe_loss_fro": safe_losses[0]["fro_loss"],
        "safe_preserves_j": bool(safe_losses[0]["fro_loss"] <= 1e-10),
        "unsafe_z": z_bad.tolist(),
        "unsafe_norm_Fz": float(np.linalg.norm(F @ z_bad)),
        "j_after_unsafe": J_after_unsafe.tolist(),
        "j_after_unsafe_eigs": J_unsafe_eigs.tolist(),
        "diff_eigs": w_unsafe.tolist(),
        "unsafe_decreases_loewner": loewner_ok,
        "interpretation_note": (
            "J_x always eliminates the map tangent A as well as current C0, so "
            "the existing nuisance of Theorem 6 is Ran([C0,A]); V is computed "
            "with that full nuisance.  Restricting N to C0 alone does not "
            "guarantee preservation of this J_x when A is present."
        ),
    }


def _designed_complex_candidate(A, B, C0, seed=2601, w0=None):
    """Complex candidate with a declared nonzero complex-linear safe kernel.

    The candidate ``K = Vq M + (I - Vq Vq^T) K_rand`` is genuinely complex and
    full-rank in general, but ``V^T K`` has the declared kernel ``span(w0)``
    (``Vq`` is an orthonormal basis of the visible-pose range ``V`` and
    ``V = Vq R`` with invertible ``R``).  This makes the complex-safe kernel
    nonempty so the J_c-invariance and utility checks are non-vacuous.
    """
    rng = np.random.default_rng(seed)
    V = _visible_pose(A, B, C0)
    Vq = _col_basis(V)
    n = A.shape[0]
    d = 3
    if w0 is None:
        w0 = np.array([1.0 + 0.0j, 0.0 + 1.0j, 0.0 + 0.0j])
    m1 = rng.standard_normal(3) + 1j * rng.standard_normal(3)
    m3 = rng.standard_normal(3) + 1j * rng.standard_normal(3)
    m2 = 1j * m1  # guarantees M w0 = m1 + i*m2 = 0
    M = np.column_stack([m1, m2, m3])
    K_rand = rng.standard_normal((n, d)) + 1j * rng.standard_normal((n, d))
    W = (np.eye(n) - Vq @ Vq.T) @ K_rand if Vq.shape[1] else K_rand
    K = Vq @ M + W
    Q = np.column_stack([K.real, -K.imag])
    return K, Q, w0, M


def complex_safe_space():
    """Complex calibration-neutral admission (Theorem 6, complex case)."""
    fx = neutral_fixture()
    A, B, C0 = fx["A"], fx["B"], fx["C0"]
    V = _visible_pose(A, B, C0)
    d = 3
    Jc = np.block([
        [np.zeros((d, d)), -np.eye(d)],
        [np.eye(d), np.zeros((d, d))],
    ])
    K, Q, w0, _ = _designed_complex_candidate(A, B, C0)
    F = V.T @ Q
    FJ = F @ Jc
    stack = np.vstack([F, FJ])
    Z_S = _null_basis(stack)
    dimS = int(Z_S.shape[1])
    # For numerical bookkeeping the designed complex dimension is dim/2.
    complex_dim = dimS // 2
    Z_safe_real = _null_basis(F)

    # (i) J_c invariance of S.
    if dimS:
        P_S = Z_S @ Z_S.T
        Jc_resid = np.linalg.norm((np.eye(2 * d) - P_S) @ (Jc @ Z_S), "fro")
    else:
        Jc_resid = 0.0
    # (ii) columns of Z_S in ker F (real safe directions).
    max_abs_FZS = float(np.max(np.abs(F @ Z_S))) if dimS else 0.0
    J_before = Jx(A, B, C0)
    safe_j_loss = None
    if dimS:
        J_after_S = _add_current(A, B, C0, Q @ Z_S)
        safe_j_loss = float(np.linalg.norm(J_before - J_after_S, "fro"))

    # (iii) z in ker F not in ker(F Jc): Jc-image unsafe.
    unsafe_record = {}
    z_candidates = Z_safe_real
    if z_candidates.shape[1]:
        # choose the kernel vector with largest component outside S
        if dimS:
            resid = (np.eye(2 * d) - Z_S @ Z_S.T) @ z_candidates
            nrm = np.linalg.norm(resid, axis=0)
            jj = int(np.argmax(nrm))
        else:
            jj = 0
        z_unsafe = z_candidates[:, jj:jj + 1]
        nrm_FJz = float(np.linalg.norm(F @ (Jc @ z_unsafe)))
        J_after_Jcz = _add_current(A, B, C0, Q @ (Jc @ z_unsafe))
        diff = J_before - J_after_Jcz
        w_diff = _eig_diff_psd(diff)
        unsafe_record = {
            "z_in_kerF": z_unsafe.tolist(),
            "z_in_ker_F": float(np.linalg.norm(F @ z_unsafe)) <= 1e-10,
            "norm_F_Jc_z": nrm_FJz,
            "in_ker_F_Jc": bool(nrm_FJz <= 1e-10),
            "j_after_Jc_z": J_after_Jcz.tolist(),
            "diff_eigs": w_diff.tolist(),
            "decreases_loewner": bool(
                float(np.min(w_diff)) >= -1e-10 and float(np.max(w_diff)) > 1e-10
            ),
        }
    # Generic random full-rank complex candidate control.
    rng = np.random.default_rng(2602)
    K_gen = rng.standard_normal((A.shape[0], d)) + 1j * rng.standard_normal((A.shape[0], d))
    Q_gen = np.column_stack([K_gen.real, -K_gen.imag])
    F_gen = V.T @ Q_gen
    S_gen = _null_basis(np.vstack([F_gen, F_gen @ Jc]))
    Vfull_rank = int(_rank(V))
    F_rank = int(_rank(F))
    return {
        "case": "complex_safe_space",
        "seed": fx["seed"],
        "dims": {"n": fx["n"], "q": fx["q"], "p": fx["p"], "d_complex": d},
        "candidate_note": (
            "designed genuinely-complex candidate K with declared complex null "
            "vector w0=(1,i,0) of V^T K so S is nonempty and the invariance "
            "test is non-vacuous; a generic full-rank random K gives S={0} "
            "when p=d (reported as generic_control)."
        ),
        "w0_declared": {"re": [1.0, 0.0, 0.0], "im": [0.0, 1.0, 0.0]},
        "V_rank": int(Vfull_rank),
        "F_shape": list(F.shape),
        "rank_F": int(F_rank),
        "real_dim_S": dimS,
        "complex_dim_S": complex_dim,
        "lower_bound_real": 2 * (d - fx["p"]),
        "lower_bound_note": "2*(d-p)=0 with d=p=3, so the bound is not binding",
        "jc_invariance_residual_fro": float(Jc_resid),
        "jc_invariant": bool(Jc_resid <= 1e-10),
        "max_abs_F_Z_S": max_abs_FZS,
        "z_S_in_ker_F": bool(max_abs_FZS <= 1e-10),
        "safe_loss_fro_after_S": safe_j_loss,
        "safe_complex_subspace_lossless": bool(
            safe_j_loss is not None and safe_j_loss <= 1e-10
        ),
        "unsafe": unsafe_record,
        "generic_control": {
            "real_dim_S": int(S_gen.shape[1]),
            "note": "generic complex K with d=p has empty complex-safe kernel",
        },
        "pass": bool(
            Jc_resid <= 1e-10
            and max_abs_FZS <= 1e-10
            and (safe_j_loss is None or safe_j_loss <= 1e-10)
            and bool(unsafe_record.get("decreases_loewner"))
        ),
    }


def utility_projection():
    """Best-r subspace in the safe kernel for a declared PSD utility T."""
    fx = neutral_fixture()
    A, B, C0 = fx["A"], fx["B"], fx["C0"]
    V = _visible_pose(A, B, C0)
    d = 3
    Jc = np.block([
        [np.zeros((d, d)), -np.eye(d)],
        [np.eye(d), np.zeros((d, d))],
    ])
    _, Q, _, _ = _designed_complex_candidate(A, B, C0)
    F = V.T @ Q
    Z_S = _null_basis(np.vstack([F, F @ Jc]))
    sdim = int(Z_S.shape[1])
    rng = np.random.default_rng(2700)
    G = rng.standard_normal((2 * d, 2 * d))
    T = G @ G.T
    r = 1 if sdim >= 1 else 0
    rec = {
        "case": "utility_projection",
        "coeff_dim": int(2 * d),
        "real_dim_S": int(sdim),
        "chosen_r": int(r),
        "T_psd": True,
        "T_eigenvalues": np.linalg.eigvalsh(T).tolist(),
    }
    if r == 0:
        rec["note"] = "safe kernel is empty; utility maximization vacuous"
        return rec
    P_S = Z_S @ Z_S.T
    M = (P_S @ T) @ P_S
    w, U = np.linalg.eigh(M)
    Z_opt = U[:, -r:]
    trace_opt = float(np.trace(Z_opt.T @ T @ Z_opt))
    # Compare with many random r-dimensional subspaces inside S.
    n_trials = 300
    traces = []
    max_trace = -np.inf
    for _ in range(n_trials):
        H = rng.standard_normal((sdim, r))
        Hq, _ = np.linalg.qr(H)
        Ztrial = Z_S @ Hq
        tr = float(np.trace(Ztrial.T @ T @ Ztrial))
        traces.append(tr)
        max_trace = max(max_trace, tr)
    tol = 1e-8 * max(1.0, trace_opt)
    rec.update({
        "trace_claimed_maximizer": trace_opt,
        "n_random_trials": int(n_trials),
        "max_random_trace": float(max_trace),
        "mean_random_trace": float(np.mean(traces)),
        "maximizer_dominates": bool(trace_opt + tol >= max_trace),
        "leading_eigenvalue": float(w[-1]),
    })
    return rec


# ---------------------------------------------------------------------------
# Fixed random tangents (seeds 201-212)
# ---------------------------------------------------------------------------


def _random_U(seed=900, n=12):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    return q


def random_tangent_checks(seeds=None):
    """Adjacent nested-pair checks for seeds 201-212 (36 checks)."""
    if seeds is None:
        seeds = list(range(201, 213))
    n, q, p = 12, 3, 4
    U = _random_U()
    records = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        A = _scaled_orthonormal(rng, n, q)
        B = _scaled_orthonormal(rng, n, p)
        chains = {
            "C0": np.zeros((n, 0)),
            "C1": U[:, :2],
            "C2": U[:, :4],
            "C3": U[:, :6],
        }
        for r in range(3):
            Cr = chains[f"C{r}"]
            Cr1 = chains[f"C{r + 1}"]
            E = E_basis(A, Cr, Cr1)
            P_E = E @ E.T if E.shape[1] else np.zeros((n, n))
            J_r = Jx(A, B, Cr)
            J_r1 = Jx(A, B, Cr1)
            pred = B.T @ P_E @ B
            actual = J_r - J_r1
            loss_res = float(np.linalg.norm(pred - actual, "fro"))
            scale = EPS * max(1.0, float(np.linalg.norm(B, "fro") ** 2)) * n
            scaled = loss_res / scale
            Bv_r = Bv(A, B, Cr)
            Bv_r1 = Bv(A, B, Cr1)
            r_r = _rank(Bv_r)
            r_r1 = _rank(Bv_r1)
            drop_actual = r_r - r_r1
            if E.shape[1]:
                inter = r_r + E.shape[1] - _rank(np.hstack([Bv_r, E]))
            else:
                inter = 0
            records.append({
                "seed": int(seed),
                "pair": f"C{r}->C{r+1}",
                "loss_residual_fro": loss_res,
                "scaled_residual": scaled,
                "rank_bv_r": int(r_r),
                "rank_bv_r1": int(r_r1),
                "rank_drop_predicted": int(inter),
                "rank_drop_actual": int(drop_actual),
                "t5b_matches": bool(inter == drop_actual),
                "dim_E": int(E.shape[1]),
                "j_r_minus_j_r1_eigs": _eig_diff_psd(actual).tolist(),
                "b_T_P_E_b_eigs": _eig_diff_psd(pred).tolist(),
            })
    max_scaled = max((r["scaled_residual"] for r in records), default=0.0)
    max_res = max((r["loss_residual_fro"] for r in records), default=0.0)
    all_match = all(r["t5b_matches"] for r in records)
    return {
        "dims": {"n": n, "q": q, "p": p},
        "seeds": list(seeds),
        "n_checks": len(records),
        "max_scaled_residual": float(max_scaled),
        "max_loss_residual_fro": float(max_res),
        "all_t5b_matches": bool(all_match),
        "records": records,
    }


def all_e2():
    """All E2 cases in one nested dictionary (for results/e2_results.json)."""
    return {
        "exact_controls": exact_controls(),
        "neutral_admission": {
            "real_safe_space": real_safe_space(),
            "complex_safe_space": complex_safe_space(),
            "utility_projection": utility_projection(),
        },
        "random_tangents": random_tangent_checks(),
    }
