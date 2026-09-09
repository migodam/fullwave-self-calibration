"""Experiment E3: dual tangent spectra, pose-prior extension, bias-aware
calibration risk, and relative-only gate failures.

All matrices are real, unit-noise whitened, gauge-fixed tangent matrices:

* ``A``  : n x q map tangent,
* ``B``  : n x p pose tangent,
* ``C``  : n x r current/nuisance tangent (``N = Ran(K_r)``),

with ``Pi = I - P_C``, ``A_c = Pi A``, ``B_c = Pi B``.  The central objects
are the map Gram and the pose information Gram

    K0  = A_c^T A_c,
    K_e = A_c^T (I - P_{Ran B_c}) A_c,
    J_x = B_c^T (I - P_{Ran A_c}) B_c.

The experiment checks Theorem 3 (whitened dual spectra = ``1 - c_i^2`` plus
unit multiplicities), the ``Lambda_x`` pose-prior extension of the same
spectral theorem, Theorem 7 (bias-aware risk of the projected pose estimator),
an adaptive-rank post-selection caveat, and the Section-10 confounding model
used as a relative-only gate failure.
"""

from __future__ import annotations

import numpy as np

from a2val import common

EPS = float(np.finfo(float).eps)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _mat(X, name="X"):
    X = np.asarray(X, dtype=float)
    if X.ndim == 0:
        raise ValueError(f"{name} must be 1-D or 2-D")
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    return X


def _sym(X):
    X = np.asarray(X, dtype=float)
    return 0.5 * (X + X.T)


def _rank(X, rcond=None):
    """Numerical column rank under an explicit or numpy-default cutoff."""
    X = _mat(X)
    if X.shape[1] == 0:
        return 0
    s = np.linalg.svd(X, compute_uv=False)
    if s.size == 0 or s[0] <= 0.0:
        return 0
    if rcond is None:
        return int(np.linalg.matrix_rank(X))
    return int(np.count_nonzero(s > rcond * s[0]))


def _col_basis(X, rcond=None):
    """Orthonormal basis of the column range of ``X`` (compact SVD)."""
    X = _mat(X)
    n = X.shape[0]
    if X.shape[1] == 0:
        return np.zeros((n, 0))
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    if s.size == 0 or s[0] <= 0.0:
        return np.zeros((n, 0))
    if rcond is None:
        r = int(np.linalg.matrix_rank(X))
    else:
        r = int(np.count_nonzero(s > rcond * s[0]))
    if r == 0:
        return np.zeros((n, 0))
    return u[:, :r].copy()


def _id(n):
    return np.eye(n, dtype=float)


def _eigh_sorted(X):
    return np.sort(np.linalg.eigvalsh(_sym(X)))


def map_pose_matrices(A, B, C):
    """Project out the current nuisance and return all E3 Gram matrices."""
    A = _mat(A, "A")
    B = _mat(B, "B")
    C = _mat(C, "C")
    n = A.shape[0]
    if B.shape[0] != n or (C.shape[1] and C.shape[0] != n):
        raise ValueError("A, B, C must share the leading dimension")
    Pi = _id(n) - common.orth_proj(C) if C.shape[1] else _id(n)
    Ac = Pi @ A
    Bc = Pi @ B
    P_Bc = common.orth_proj(Bc) if Bc.shape[1] else np.zeros((n, n))
    P_Ac = common.orth_proj(Ac) if Ac.shape[1] else np.zeros((n, n))
    K0 = Ac.T @ Ac
    Ke = Ac.T @ (_id(n) - P_Bc) @ Ac
    Jx = Bc.T @ (_id(n) - P_Ac) @ Bc
    return {
        "A_c": Ac,
        "B_c": Bc,
        "K0": _sym(K0),
        "K_e": _sym(Ke),
        "J_x": _sym(Jx),
        "P_C": _id(n) - Pi,
    }


def _support_whitened(S, S_eff, support_tol_rel=1e-12):
    """``S^{-1/2} S_eff S^{-1/2}`` restricted to the support of PSD ``S``.

    Eigen-directions of ``S`` below the support tolerance are returned as
    undefined (they are never assigned a spectral value 0 or 1).
    """
    S = _sym(S)
    S_eff = _sym(S_eff)
    q = int(S.shape[0])
    w, U = np.linalg.eigh(S)
    smax = float(np.max(np.abs(w))) if q else 0.0
    if smax <= 0.0:
        return {
            "values": np.zeros(0, dtype=float),
            "support_dim": 0,
            "q": q,
            "undefined_dim": q,
            "support_tol": 0.0,
        }
    tol = support_tol_rel * smax
    support = w > tol
    r = int(np.count_nonzero(support))
    if r == 0:
        return {
            "values": np.zeros(0, dtype=float),
            "support_dim": 0,
            "q": q,
            "undefined_dim": q,
            "support_tol": float(tol),
        }
    ws = np.clip(w[support], 0.0, None)
    Us = U[:, support]
    D = np.diag(1.0 / np.sqrt(ws))
    M = ((D @ Us.T) @ S_eff) @ (Us @ D)
    vals = np.linalg.eigvalsh(_sym(M))
    vals = np.sort(vals)
    # Clamp only numerical pollution below the tolerance used for zeros.
    vals = np.where(vals > -support_tol_rel, np.maximum(vals, 0.0), vals)
    return {
        "values": np.asarray(vals, dtype=float),
        "support_dim": int(r),
        "q": q,
        "undefined_dim": int(q - r),
        "support_tol": float(tol),
    }


def dual_normalized_spectra(A_c, B_c, support_tol_rel=1e-12):
    """Whitened map and pose spectra on their Gram supports."""
    K0 = _sym(A_c.T @ A_c)
    B0 = _sym(B_c.T @ B_c)
    n = A_c.shape[0]
    P_B = common.orth_proj(B_c) if B_c.shape[1] else np.zeros((n, n))
    P_A = common.orth_proj(A_c) if A_c.shape[1] else np.zeros((n, n))
    Ke = _sym(A_c.T @ (_id(n) - P_B) @ A_c)
    Jx = _sym(B_c.T @ (_id(n) - P_A) @ B_c)
    map_spec = _support_whitened(K0, Ke, support_tol_rel=support_tol_rel)
    pose_spec = _support_whitened(B0, Jx, support_tol_rel=support_tol_rel)
    return {
        "map": map_spec,
        "pose": pose_spec,
        "K0": K0,
        "B0": B0,
        "K_e": Ke,
        "J_x": Jx,
    }


def canonical_info(A_c, B_c, rcond=1e-12):
    """Canonical-correlation data of the two visible ranges."""
    QA = _col_basis(A_c, rcond=rcond)
    QB = _col_basis(B_c, rcond=rcond)
    a = int(QA.shape[1])
    b = int(QB.shape[1])
    Ccan = QA.T @ QB if a and b else np.zeros((a, b))
    c = np.clip(np.linalg.svd(Ccan, compute_uv=False), 0.0, 1.0)
    return {"QA": QA, "QB": QB, "a": a, "b": b,
            "Ccan": Ccan, "c": c, "rank_Ccan": int(_rank(Ccan, rcond=rcond))}


def expected_dual_spectra(info):
    """Theorem-3 predicted whitened map/pose spectra from singular values."""
    a, b, c = info["a"], info["b"], info["c"]
    k = min(a, b)
    paired = 1.0 - c[:k] ** 2
    map_vals = np.sort(np.concatenate(
        [paired, np.ones(max(0, a - k))])) if a else np.zeros(0)
    pose_vals = np.sort(np.concatenate(
        [paired, np.ones(max(0, b - k))])) if b else np.zeros(0)
    return map_vals, pose_vals


def _spectral_match(values, expected, tol=1e-10):
    if len(values) != len(expected):
        return np.inf
    return float(np.max(np.abs(np.asarray(values) - np.asarray(expected))))


def _intersection_dim(A_c, B_c):
    a = _rank(A_c)
    b = _rank(B_c)
    if a == 0 or b == 0:
        return 0
    stacked = np.hstack([_col_basis(A_c), _col_basis(B_c)])
    return int(a + b - _rank(stacked))


# ---------------------------------------------------------------------------
# Part A: dual tangent spectra
# ---------------------------------------------------------------------------


def exact_small_dual_spectra():
    """A1: two exactly known principal angles with C empty."""
    n = 4
    A = _id(n)[:, :2]
    records = {}
    passes = []
    for t, label in ((np.pi / 4, "pi/4"), (np.pi / 2, "pi/2")):
        v = np.cos(t) * _id(n)[:, 0] + np.sin(t) * _id(n)[:, 2]
        B = np.column_stack([v, _id(n)[:, 3]])
        C = np.zeros((n, 0))
        grams = map_pose_matrices(A, B, C)
        Ac, Bc = grams["A_c"], grams["B_c"]
        norm = dual_normalized_spectra(Ac, Bc)
        info = canonical_info(Ac, Bc)
        c = info["c"]
        expected_map, expected_pose = expected_dual_spectra(info)
        map_res = _spectral_match(norm["map"]["values"], expected_map)
        pose_res = _spectral_match(norm["pose"]["values"], expected_pose)
        zero_pred = _intersection_dim(Ac, Bc)
        zero_map = int(np.count_nonzero(norm["map"]["values"] <= 1e-9))
        zero_pose = int(np.count_nonzero(norm["pose"]["values"] <= 1e-9))
        unit_map_pred = info["a"] - info["rank_Ccan"]
        unit_pose_pred = info["b"] - info["rank_Ccan"]
        unit_map_obs = int(np.count_nonzero(norm["map"]["values"] >= 1.0 - 1e-9))
        unit_pose_obs = int(np.count_nonzero(norm["pose"]["values"] >= 1.0 - 1e-9))
        ok = bool(
            map_res <= 1e-10 and pose_res <= 1e-10
            and zero_map == zero_pred == zero_pose
            and unit_map_obs == unit_map_pred
            and unit_pose_obs == unit_pose_pred
        )
        passes.append(ok)
        records[label] = {
            "t": float(t),
            "n": n,
            "A": A.tolist(),
            "B": B.tolist(),
            "C": C.tolist(),
            "c": c.tolist(),
            "expected_c": [np.cos(t), 0.0],
            "map_spectrum": norm["map"]["values"].tolist(),
            "pose_spectrum": norm["pose"]["values"].tolist(),
            "expected_map_spectrum": expected_map.tolist(),
            "expected_pose_spectrum": expected_pose.tolist(),
            "map_residual": map_res,
            "pose_residual": pose_res,
            "zero_multiplicity_predicted": zero_pred,
            "zero_multiplicity_map": zero_map,
            "zero_multiplicity_pose": zero_pose,
            "unit_multiplicity_map_predicted": unit_map_pred,
            "unit_multiplicity_map": unit_map_obs,
            "unit_multiplicity_pose_predicted": unit_pose_pred,
            "unit_multiplicity_pose": unit_pose_obs,
            "support_map": norm["map"]["support_dim"],
            "support_pose": norm["pose"]["support_dim"],
            "undefined_map": norm["map"]["undefined_dim"],
            "undefined_pose": norm["pose"]["undefined_dim"],
            "pass": ok,
        }
    return {
        "case": "exact_small",
        "records": records,
        "all_pass": bool(all(passes)),
    }


def random_dual_spectral_fixture(seed=2150):
    """Fixed random fixture for A2/A4 and the B-part pose-prior checks."""
    rng = np.random.default_rng(seed)
    n, q, p, r_c = 10, 4, 3, 2
    A, _ = np.linalg.qr(rng.standard_normal((n, q)))
    B, _ = np.linalg.qr(rng.standard_normal((n, p)))
    C, _ = np.linalg.qr(rng.standard_normal((n, r_c)))
    # Deterministic, invertible row re-shaping so the ranges are generic but
    # never collinear with the coordinate axes.
    A = A + 0.15 * rng.standard_normal((n, q))
    B = B + 0.15 * rng.standard_normal((n, p))
    C = C + 0.15 * rng.standard_normal((n, r_c))
    return {
        "seed": int(seed),
        "n": n,
        "q": q,
        "p": p,
        "r_c": r_c,
        "A": A,
        "B": B,
        "C": C,
    }


def theorem3_fixture_check(seed=2150, rcond=1e-12):
    """A2: random nuisance check of Theorem 3 on one fixture."""
    fx = random_dual_spectral_fixture(seed)
    A, B, C = fx["A"], fx["B"], fx["C"]
    grams = map_pose_matrices(A, B, C)
    Ac, Bc = grams["A_c"], grams["B_c"]
    norm = dual_normalized_spectra(Ac, Bc)
    info = canonical_info(Ac, Bc, rcond=rcond)
    expected_map, expected_pose = expected_dual_spectra(info)
    map_res = _spectral_match(norm["map"]["values"], expected_map)
    pose_res = _spectral_match(norm["pose"]["values"], expected_pose)
    zero_pred = _intersection_dim(Ac, Bc)
    zero_obs = int(np.count_nonzero(norm["map"]["values"] <= 1e-9))
    unit_map_obs = int(np.count_nonzero(norm["map"]["values"] >= 1.0 - 1e-9))
    unit_pose_obs = int(np.count_nonzero(norm["pose"]["values"] >= 1.0 - 1e-9))
    unit_map_pred = info["a"] - info["rank_Ccan"]
    unit_pose_pred = info["b"] - info["rank_Ccan"]
    raw_Ke = np.sort(np.linalg.eigvalsh(_sym(grams["K_e"])))
    raw_Jx = np.sort(np.linalg.eigvalsh(_sym(grams["J_x"])))
    all_ok = bool(
        map_res <= 1e-10 and pose_res <= 1e-10
        and zero_obs == zero_pred
        and unit_map_obs == unit_map_pred
        and unit_pose_obs == unit_pose_pred
        and norm["map"]["undefined_dim"] == 0
        and norm["pose"]["undefined_dim"] == 0
        and _rank(Ac, rcond=rcond) == info["a"]
        and _rank(Bc, rcond=rcond) == info["b"]
    )
    return {
        "case": "random_nuisance",
        "dims": fx,
        "a": info["a"],
        "b": info["b"],
        "c": info["c"].tolist(),
        "rank_Ccan": info["rank_Ccan"],
        "canonical_1_minus_c2": (1.0 - info["c"] ** 2).tolist(),
        "map_spectrum_normalized": norm["map"]["values"].tolist(),
        "pose_spectrum_normalized": norm["pose"]["values"].tolist(),
        "expected_map_normalized": expected_map.tolist(),
        "expected_pose_normalized": expected_pose.tolist(),
        "map_residual": map_res,
        "pose_residual": pose_res,
        "raw_K_e_eigs": raw_Ke.tolist(),
        "raw_J_x_eigs": raw_Jx.tolist(),
        "zero_multiplicity_predicted": zero_pred,
        "zero_multiplicity_observed": zero_obs,
        "unit_map_multiplicity_predicted": unit_map_pred,
        "unit_map_multiplicity_observed": unit_map_obs,
        "unit_pose_multiplicity_predicted": unit_pose_pred,
        "unit_pose_multiplicity_observed": unit_pose_obs,
        "raw_spectra_differ_from_normalized": bool(
            not np.allclose(np.sort(raw_Ke), norm["map"]["values"], atol=1e-8)
            or not np.allclose(np.sort(raw_Jx), norm["pose"]["values"], atol=1e-8)
        ),
        "max_spectral_residual": max(map_res, pose_res),
        "pass": all_ok,
    }


def rho_undefined_control():
    """A3: K0 and B0 null directions are reported undefined, never 0 or 1."""
    n = 4
    C = _id(n)[:, 3:4]
    A = np.column_stack([_id(n)[:, 0], _id(n)[:, 1], _id(n)[:, 3]])
    B = np.column_stack([_id(n)[:, 2], _id(n)[:, 3]])
    grams = map_pose_matrices(A, B, C)
    Ac, Bc = grams["A_c"], grams["B_c"]
    norm = dual_normalized_spectra(Ac, Bc)
    info = canonical_info(Ac, Bc)
    expected_map, expected_pose = expected_dual_spectra(info)
    map_res = _spectral_match(norm["map"]["values"], expected_map)
    pose_res = _spectral_match(norm["pose"]["values"], expected_pose)
    k0_vals = np.linalg.svd(grams["K0"], compute_uv=False)
    b0_vals = np.linalg.svd(norm["B0"], compute_uv=False)
    pass_ = bool(
        norm["map"]["undefined_dim"] == 1
        and norm["pose"]["undefined_dim"] == 1
        and map_res <= 1e-10 and pose_res <= 1e-10
        and info["a"] == 2 and info["b"] == 1
        and all(v == 1.0 for v in norm["map"]["values"])
        and norm["pose"]["values"][0] == 1.0
    )
    return {
        "case": "rho_undefined_control",
        "n": n,
        "A": A.tolist(),
        "B": B.tolist(),
        "C": C.tolist(),
        "map_spectrum": norm["map"]["values"].tolist(),
        "pose_spectrum": norm["pose"]["values"].tolist(),
        "expected_map": expected_map.tolist(),
        "expected_pose": expected_pose.tolist(),
        "map_residual": map_res,
        "pose_residual": pose_res,
        "support_map": norm["map"]["support_dim"],
        "support_pose": norm["pose"]["support_dim"],
        "undefined_map": norm["map"]["undefined_dim"],
        "undefined_pose": norm["pose"]["undefined_dim"],
        "K0_singular_values": k0_vals.tolist(),
        "B0_singular_values": b0_vals.tolist(),
        "note": (
            "A_c has an exact null direction inside Ran(C) (third A column) "
            "and B_c has one inside Ran(C) as well. Those directions are "
            "excluded from the normalized spectra; they are not assigned rho "
            "0 or 1.  Remaining visible map modes are retained (rho=1) because "
            "Ran B_c is orthogonal to Ran A_c."
        ),
        "pass": pass_,
    }


def scaling_control(seed=2150):
    """A4: scaling A,B,C by s leaves canonical rho, squares J_x."""
    fx = random_dual_spectral_fixture(seed)
    base = theorem3_fixture_check(seed=seed)
    base_norm_map = np.asarray(base["map_spectrum_normalized"])
    base_norm_pose = np.asarray(base["pose_spectrum_normalized"])
    base_jx_fro = float(np.linalg.norm(
        map_pose_matrices(fx["A"], fx["B"], fx["C"])["J_x"], "fro"
    ))
    records = []
    ok = True
    for s in (0.1, 1.0, 10.0):
        A, B, C = s * fx["A"], s * fx["B"], s * fx["C"]
        grams = map_pose_matrices(A, B, C)
        norm = dual_normalized_spectra(grams["A_c"], grams["B_c"])
        map_d = float(np.max(np.abs(np.asarray(norm["map"]["values"]) - base_norm_map)))
        pose_d = float(np.max(np.abs(np.asarray(norm["pose"]["values"]) - base_norm_pose)))
        jx_fro = float(np.linalg.norm(grams["J_x"], "fro"))
        s2_ratio = jx_fro / (s * s * base_jx_fro)
        good = bool(map_d <= 1e-10 and pose_d <= 1e-10 and abs(s2_ratio - 1.0) <= 1e-10)
        ok = ok and good
        records.append({
            "scale": float(s),
            "map_spectrum": norm["map"]["values"].tolist(),
            "pose_spectrum": norm["pose"]["values"].tolist(),
            "map_abs_delta_vs_unscaled": map_d,
            "pose_abs_delta_vs_unscaled": pose_d,
            "J_x_fro": jx_fro,
            "J_x_fro_ratio_to_s2": s2_ratio,
            "pass": good,
        })
    return {
        "case": "scaling_control",
        "unscaled_J_x_fro": base_jx_fro,
        "records": records,
        "pass": bool(ok),
    }


def b_zero_control():
    """A5: B=0 gives full map retention but zero pose information."""
    n = 4
    A = _id(n)[:, :2]
    B = np.zeros((n, 2))
    C = np.zeros((n, 0))
    grams = map_pose_matrices(A, B, C)
    norm = dual_normalized_spectra(grams["A_c"], grams["B_c"])
    jx_abs = float(np.linalg.norm(grams["J_x"], "fro"))
    jx_sv = np.linalg.svd(grams["J_x"], compute_uv=False).tolist()
    return {
        "case": "b_zero_control",
        "n": n,
        "map_spectrum": norm["map"]["values"].tolist(),
        "pose_support_dim": norm["pose"]["support_dim"],
        "pose_spectrum": norm["pose"]["values"].tolist(),
        "pose_undefined_dim": norm["pose"]["undefined_dim"],
        "J_x_singular_values": jx_sv,
        "J_x_fro": jx_abs,
        "map_retention_all_one": bool(
            len(norm["map"]["values"]) == 2
            and np.allclose(norm["map"]["values"], 1.0, atol=1e-12)
        ),
        "relative_only_gate_failure": (
            "map retention rho=1 is maximal while pose information J_x=0: a "
            "map-only (relative retention) gate would pass despite no pose "
            "information being recoverable."
        ),
        "pass": bool(
            len(norm["map"]["values"]) == 2
            and np.allclose(norm["map"]["values"], 1.0, atol=1e-12)
            and jx_abs == 0.0
        ),
    }


def dual_spectra_all():
    return {
        "exact_small": exact_small_dual_spectra(),
        "random_nuisance": theorem3_fixture_check(),
        "rho_undefined_control": rho_undefined_control(),
        "scaling_control": scaling_control(),
        "b_zero_control": b_zero_control(),
    }


# ---------------------------------------------------------------------------
# Part B: pose-prior extension
# ---------------------------------------------------------------------------


def K_eL(A_c, B_c, Lambda_x):
    """Pose-prior regularized map Gram in the map-coefficient space."""
    Lambda_x = _mat(Lambda_x, "Lambda_x")
    M = _sym(B_c.T @ B_c + Lambda_x)
    Mpinv = np.linalg.pinv(M)
    return _sym(A_c.T @ A_c - (A_c.T @ B_c) @ Mpinv @ (B_c.T @ A_c))


def _min_eig_diff(Low, High):
    return float(np.min(np.linalg.eigvalsh(_sym(High - Low))))


def _loewner_le(Low, High, tol=1e-10):
    return bool(_min_eig_diff(Low, High) >= -tol)


def pose_prior_checks(seed=2150):
    """B1/B2: Loewner sandwich, monotonicity, variational and augmented IDs."""
    fx = random_dual_spectral_fixture(seed)
    grams = map_pose_matrices(fx["A"], fx["B"], fx["C"])
    Ac, Bc = grams["A_c"], grams["B_c"]
    K0, Ke = grams["K0"], grams["K_e"]
    p = int(Bc.shape[1])
    Lam = np.diag([10.0, 1.0, 0.1])
    K_lam = K_eL(Ac, Bc, Lam)
    # Loewner sandwich K_e <= K_eL <= K0.
    min_low = _min_eig_diff(Ke, K_lam)
    min_up = _min_eig_diff(K_lam, K0)
    sandwich_ok = _loewner_le(Ke, K_lam) and _loewner_le(K_lam, K0)
    # Monotone under Lambda scaling.
    monotonic_records = []
    mono_ok = True
    for scale in (0.5, 1.0, 2.0):
        Ls = scale * Lam
        Ks = K_eL(Ac, Bc, Ls)
        rec = {
            "scale": float(scale),
            "eigenvalues": np.linalg.eigvalsh(_sym(Ks)).tolist(),
        }
        monotonic_records.append(rec)
    diffs = {
        "min_eig_K(0.5L)-K_e": _min_eig_diff(Ke, K_eL(Ac, Bc, 0.5 * Lam)),
        "min_eig_K(1.0L)-K(0.5L)": _min_eig_diff(
            K_eL(Ac, Bc, 0.5 * Lam), K_eL(Ac, Bc, Lam)),
        "min_eig_K(2.0L)-K(1.0L)": _min_eig_diff(
            K_eL(Ac, Bc, Lam), K_eL(Ac, Bc, 2.0 * Lam)),
        "min_eig_K0-K(2.0L)": _min_eig_diff(K_eL(Ac, Bc, 2.0 * Lam), K0),
    }
    mono_ok = all(v >= -1e-10 for v in diffs.values())
    # Variational min identity for five random v directions.
    rng = np.random.default_rng(2160)
    M = _sym(Bc.T @ Bc + Lam)
    var_records = []
    var_ok = True
    for _ in range(5):
        v = rng.standard_normal(Ac.shape[1])
        h_opt = np.linalg.solve(M, Bc.T @ Ac @ v)
        min_direct = float(
            (Ac @ v - Bc @ h_opt) @ (Ac @ v - Bc @ h_opt)
            + (Lam @ h_opt) @ h_opt
        )
        quad = float(v @ (K_lam @ v))
        res = abs(min_direct - quad)
        var_ok = var_ok and res <= 1e-8
        var_records.append({
            "v": v.tolist(),
            "h_star": h_opt.tolist(),
            "variational_min": min_direct,
            "v^T K_eL v": quad,
            "residual": res,
        })
    # Augmented principal-angle identity.
    A_aug = np.vstack([Ac, np.zeros((p, Ac.shape[1]))])
    B_aug = np.vstack([Bc, np.linalg.cholesky(Lam).T])
    norm_aug = dual_normalized_spectra(A_aug, B_aug)
    info_aug = canonical_info(A_aug, B_aug)
    expected_aug, _ = expected_dual_spectra(info_aug)
    aug_res = _spectral_match(norm_aug["map"]["values"], expected_aug)
    # Compare to the closed-form K_eL spectrum directly.
    direct_spec = _support_whitened(K0, K_lam)
    direct_res = _spectral_match(direct_spec["values"], norm_aug["map"]["values"])
    # Lambda=0 recovers K_e.
    K_zero = K_eL(Ac, Bc, np.zeros((p, p)))
    recovery_res = float(np.linalg.norm(K_zero - Ke, "fro"))
    all_ok = bool(
        sandwich_ok and mono_ok and var_ok
        and aug_res <= 1e-10 and direct_res <= 1e-10
        and recovery_res <= 1e-10
    )
    return {
        "case": "pose_prior_extension",
        "Lambda_x": Lam.tolist(),
        "sandwich": {
            "K_e<=K_eL": sandwich_ok,
            "K_eL<=K0": sandwich_ok,
            "min_eig_K_eL-K_e": min_low,
            "min_eig_K0-K_eL": min_up,
            "K_e_eigs": np.linalg.eigvalsh(_sym(Ke)).tolist(),
            "K_eL_eigs": np.linalg.eigvalsh(_sym(K_lam)).tolist(),
            "K0_eigs": np.linalg.eigvalsh(_sym(K0)).tolist(),
        },
        "monotonicity": {
            "pass": mono_ok,
            "records": monotonic_records,
            "min_eig_differences": diffs,
        },
        "variational": {
            "pass": var_ok,
            "max_residual": max(r["residual"] for r in var_records),
            "records": var_records,
        },
        "augmented_principal_angles": {
            "a_aug_rank": info_aug["a"],
            "b_aug_rank": info_aug["b"],
            "c_aug": info_aug["c"].tolist(),
            "expected_spectrum": expected_aug.tolist(),
            "augmented_whitened_spectrum": norm_aug["map"]["values"].tolist(),
            "closed_form_whitened_spectrum": direct_spec["values"].tolist(),
            "augmented_residual": aug_res,
            "closed_form_vs_augmented_residual": direct_res,
        },
        "lambda_zero_recovery": {
            "fro_residual_vs_K_e": recovery_res,
            "recovers_K_e": bool(recovery_res <= 1e-10),
        },
        "pass": all_ok,
    }


# ---------------------------------------------------------------------------
# Part C: Theorem 7 bias-aware risk
# ---------------------------------------------------------------------------


def theorem7_fixture(seed=2200):
    """Random fixed design with N, B, D_r, and PSD metric M_x."""
    rng = np.random.default_rng(seed)
    n, n_n, p, n_d = 12, 3, 4, 2
    N = rng.standard_normal((n, n_n))
    B = rng.standard_normal((n, p))
    D = rng.standard_normal((n, n_d))
    G = rng.standard_normal((p, p))
    Mx = G.T @ G + 0.1 * np.eye(p)
    Bv = (_id(n) - common.orth_proj(N)) @ B
    return {
        "seed": int(seed),
        "n": n,
        "rank_N": n_n,
        "p": p,
        "rank_D": n_d,
        "N": N,
        "B": B,
        "D_r": D,
        "M_x": _sym(Mx),
        "B_v": Bv,
        "rank_B_v": int(np.linalg.matrix_rank(Bv)),
    }


def theorem7_analytic(fx):
    Bv, D, Mx = fx["B_v"], fx["D_r"], fx["M_x"]
    p = int(Bv.shape[1])
    Jx = _sym(Bv.T @ Bv)
    Jx_inv = np.linalg.inv(Jx)
    Bv_dag = np.linalg.pinv(Bv)
    w, V = np.linalg.eigh(Mx)
    Mhalf = (V * np.sqrt(w)) @ V.T
    var_tr = float(np.trace(Mx @ Jx_inv))
    MBD = Mhalf @ Bv_dag @ D
    bias2 = float(np.linalg.norm(MBD, 2) ** 2)
    return {
        "J_x": Jx,
        "J_x_inv": Jx_inv,
        "B_v_dag": Bv_dag,
        "M_x_half": Mhalf,
        "M_x_B_v_dag_D": MBD,
        "variance_trace": var_tr,
        "worst_bias_squared": bias2,
        "analytic_total_risk": var_tr + bias2,
    }


def theorem7_mc_seed(fx, seed, m=1000):
    """One Monte-Carlo seed for the fixed design."""
    ana = theorem7_analytic(fx)
    rng = np.random.default_rng(seed)
    Bv_dag, Mhalf = ana["B_v_dag"], ana["M_x_half"]
    D, z_star = fx["D_r"], None
    # Worst deterministic direction for the bias term.
    u, s, vh = np.linalg.svd(ana["M_x_B_v_dag_D"], full_matrices=False)
    z_star = vh.T[:, 0]
    # Coordinate-space bias added to h_hat: B_v^+ D_r z* (the metric root is
    # applied only when computing || . ||_{M_x}).
    worst_bias_vec = Bv_dag @ (D @ z_star)
    E = rng.standard_normal((m, fx["n"]))
    Hvar = E @ Bv_dag.T                      # m x p, B_v^+ eps_i
    TVar = np.einsum("ij,ij->i", Hvar @ Mhalf.T, Hvar @ Mhalf.T)
    var_est = float(np.mean(TVar))
    var_se = float(np.std(TVar, ddof=1) / np.sqrt(m))
    S_cov = np.cov(Hvar, rowvar=False)
    cov_rel = float(
        np.linalg.norm(S_cov - ana["J_x_inv"], "fro")
        / np.linalg.norm(ana["J_x_inv"], "fro")
    )
    # Worst-bias + noise realizations.
    Hbias = E @ Bv_dag.T + worst_bias_vec
    Tbias = np.einsum("ij,ij->i", Hbias @ Mhalf.T, Hbias @ Mhalf.T)
    bias_est = float(np.mean(Tbias))
    bias_se = float(np.std(Tbias, ddof=1) / np.sqrt(m))
    total_ana = ana["analytic_total_risk"]
    var_diff = abs(var_est - ana["variance_trace"])
    bias_diff = abs(bias_est - total_ana)
    var_pass = bool(var_diff <= 3.0 * var_se)
    bias_pass = bool(bias_diff <= 3.0 * bias_se)
    return {
        "seed": int(seed),
        "m": int(m),
        "variance_estimate": var_est,
        "variance_analytic": ana["variance_trace"],
        "variance_se": var_se,
        "variance_rel_error": abs(var_est - ana["variance_trace"])
        / ana["variance_trace"],
        "variance_within_3se": var_pass,
        "sample_cov_trace_of_Mx": float(np.trace(Mhalf @ S_cov @ Mhalf)),
        "sample_cov_vs_Jx_inv_fro_rel": cov_rel,
        "bias_estimate_total_risk": bias_est,
        "bias_analytic_total_risk": total_ana,
        "bias_se": bias_se,
        "bias_rel_error": abs(bias_est - total_ana) / total_ana,
        "bias_within_3se": bias_pass,
        "z_star": z_star.tolist(),
        "worst_bias_squared": float(np.sum(worst_bias_vec ** 2)),
    }


def theorem7_checks(seeds=None, m=1000):
    """C1-C4 under the Theorem 7 banner."""
    if seeds is None:
        seeds = list(range(301, 311))
    fx = theorem7_fixture()
    if fx["rank_B_v"] != fx["p"]:
        raise RuntimeError("random fixture lost full column rank of B_v")
    ana = theorem7_analytic(fx)
    records = [theorem7_mc_seed(fx, seed, m=m) for seed in seeds]
    # Deterministic bound.
    Mhalf = ana["M_x_half"]
    W = fx["B_v"] @ np.linalg.inv(Mhalf)
    sW = np.linalg.svd(W, compute_uv=False)
    sigma_min = float(sW[-1])
    beta_r = float(np.linalg.norm(fx["D_r"], 2))
    delta = 0.0
    bound = (delta + beta_r) / sigma_min
    worst_err = float(np.sqrt(ana["worst_bias_squared"]))
    bound_ok = bool(worst_err <= bound * (1.0 + 1e-12) + 1e-12)
    max_var_seed = max(r["variance_rel_error"] for r in records)
    max_bias_seed = max(r["bias_rel_error"] for r in records)
    all_var = all(r["variance_within_3se"] for r in records)
    all_bias = all(r["bias_within_3se"] for r in records)
    return {
        "fixture": {
            "seed": fx["seed"],
            "dims": {"n": fx["n"], "rank_N": fx["rank_N"],
                     "p": fx["p"], "rank_D": fx["rank_D"]},
            "rank_B_v": fx["rank_B_v"],
            "N": fx["N"].tolist(),
            "B": fx["B"].tolist(),
            "D_r": fx["D_r"].tolist(),
            "M_x": fx["M_x"].tolist(),
        },
        "analytic": {
            "variance_trace": ana["variance_trace"],
            "worst_bias_squared": ana["worst_bias_squared"],
            "analytic_total_risk": ana["analytic_total_risk"],
        },
        "deterministic_bound": {
            "sigma_min_B_v_Mx_invhalf": sigma_min,
            "beta_r": beta_r,
            "delta": delta,
            "bound": bound,
            "worst_analytic_Mx_error": worst_err,
            "worst_error_squared": ana["worst_bias_squared"],
            "bound_satisfied": bound_ok,
            "relation_note": (
                "||M^{1/2} B_v^+ D_r z|| <= ||D_r z||/sigma_min(B_v M^{-1/2}), "
                "so worst-error <= bound with equality when z aligns with the "
                "left singular direction of W^+ D_r and ||D_r z||=||D_r||."
            ),
        },
        "mc": {
            "m_per_seed": m,
            "seeds": list(seeds),
            "records": records,
            "all_variance_within_3se": all_var,
            "all_bias_within_3se": all_bias,
            "max_variance_rel_error": max_var_seed,
            "max_bias_rel_error": max_bias_seed,
            "min_sample_cov_rel_error": min(r["sample_cov_vs_Jx_inv_fro_rel"]
                                            for r in records),
            "max_sample_cov_rel_error": max(r["sample_cov_vs_Jx_inv_fro_rel"]
                                            for r in records),
        },
        "pass": bool(all_var and all_bias and bound_ok),
    }


def adaptive_rank_caveat(seed=2300, n_draws=500):
    """Small post-selection caveat control (not a theorem check)."""
    rng = np.random.default_rng(seed)
    n, p = 6, 2
    B = np.zeros((n, p))
    B[0, 0] = 1.0
    B[1, 1] = 1.0
    n0 = np.array([0.7, 0.0, 0.5, 0.0, 0.0, 0.0])
    n0 /= np.linalg.norm(n0)
    n1raw = np.array([0.0, 0.6, 0.2, 0.4, 0.0, 0.0])
    n1raw /= np.linalg.norm(n1raw)
    N1 = np.column_stack([n0, n1raw])
    # Penalized-residual selector: residual reduction from the extra nuisance
    # column must exceed tau (so the larger rank is chosen only for
    # "significant" reductions).
    tau = 0.5
    h_true = np.array([1.0, 0.5])
    a0 = 0.8
    eps = rng.standard_normal((n_draws, n))
    y = eps + B @ h_true + n0 * a0
    # Selection statistic is d^T y where d is the unit component of n1raw
    # orthogonal to [B, n0]; under the truth d^T eps ~ N(0,1).
    S_full = np.column_stack([B, n0])
    Qs, _ = np.linalg.qr(S_full)
    d = (_id(n) - Qs @ Qs.T) @ n1raw
    d /= np.linalg.norm(d)
    scores = (y @ d) ** 2
    sel = scores > tau
    n_sel1 = int(np.count_nonzero(sel))
    n_sel0 = n_draws - n_sel1
    Bv0 = (_id(n) - common.orth_proj(n0[:, None])) @ B
    Bv1 = (_id(n) - common.orth_proj(N1)) @ B
    J0 = _sym(Bv0.T @ Bv0)
    J1 = _sym(Bv1.T @ Bv1)
    J0_inv, J1_inv = np.linalg.inv(J0), np.linalg.inv(J1)
    rank0 = int(np.linalg.matrix_rank(Bv0))
    rank1 = int(np.linalg.matrix_rank(Bv1))

    def est(N):
        if N.shape[1]:
            P = _id(n) - common.orth_proj(N)
        else:
            P = _id(n)
        return np.linalg.pinv(P @ B) @ (P @ y.T)

    H0 = est(n0[:, None]) - h_true[:, None]
    H1 = est(N1) - h_true[:, None]
    Hsel = np.where(sel[:, None], H1.T, H0.T)
    cov_pool = np.cov(Hsel, rowvar=False)
    frac1 = n_sel1 / n_draws
    weighted_inv = (1.0 - frac1) * J0_inv + frac1 * J1_inv
    majority = 1 if n_sel1 >= n_sel0 else 0
    maj_inv = J1_inv if majority == 1 else J0_inv
    fro_diff_majority = float(np.linalg.norm(cov_pool - maj_inv, "fro"))
    fro_diff_weighted = float(np.linalg.norm(cov_pool - weighted_inv, "fro"))
    # Conditional covariances (selection is independent of estimator noise in
    # this small design, so the mismatch is the *unconditional* mixture).
    cond0 = np.cov(H0.T[~sel], rowvar=False) if n_sel0 > 1 else None
    cond1 = np.cov(H1.T[sel], rowvar=False) if n_sel1 > 1 else None
    cond0_res = (float(np.linalg.norm(cond0 - J0_inv, "fro"))
                 if cond0 is not None else None)
    cond1_res = (float(np.linalg.norm(cond1 - J1_inv, "fro"))
                 if cond1 is not None else None)
    return {
        "case": "adaptive_rank_caveat",
        "seed": int(seed),
        "n_draws": int(n_draws),
        "tau": tau,
        "n_selected_N0": n_sel0,
        "n_selected_N1": n_sel1,
        "fraction_selected_N1": frac1,
        "ranks_visible": {"N0": rank0, "N1": rank1},
        "B_v_ranks_full": bool(rank0 == p and rank1 == p),
        "J0_inv": J0_inv.tolist(),
        "J1_inv": J1_inv.tolist(),
        "pooled_empirical_covariance": cov_pool.tolist(),
        "weighted_inverse_information": weighted_inv.tolist(),
        "majority_model": f"N{majority}",
        "majority_model_inverse_information": maj_inv.tolist(),
        "fro_diff_pooled_vs_majority_inverse_info": fro_diff_majority,
        "fro_diff_pooled_vs_weighted_inverse_info": fro_diff_weighted,
        "conditional_cov_residuals": {
            "N0_vs_J0_inv_fro": cond0_res,
            "N1_vs_J1_inv_fro": cond1_res,
        },
        "caveat_note": (
            "The selector uses the same noise realization as the estimator. "
            "Even though each conditional (given model choice) law is close to "
            "the selected-model Gaussian, the pooled empirical covariance of "
            "the selected estimator is a mixture and is NOT equal to the "
            "inverse information of a single fixed selected model.  This is a "
            "post-selection caveat control, not a theorem check."
        ),
        "pass": bool(
            rank0 == p and rank1 == p
            and 0.05 < frac1 < 0.95
            and fro_diff_majority
            > 0.05 * max(np.linalg.norm(J0_inv, "fro"),
                         np.linalg.norm(J1_inv, "fro"))
        ),
    }


# ---------------------------------------------------------------------------
# Part D: confounding / bias-vs-variance model
# ---------------------------------------------------------------------------


def confounding_case(eps_param, sigma2=1.0, c_star=1.0):
    eps_param = float(eps_param)
    full_risk = sigma2 * (1.0 + eps_param ** (-2))
    trunc_risk = sigma2 + c_star * c_star
    better_trunc = bool(trunc_risk < full_risk)
    return {
        "eps_param": eps_param,
        "sigma2": sigma2,
        "c_star": c_star,
        "full_ls_risk": full_risk,
        "truncated_risk": trunc_risk,
        "truncation_better": better_trunc,
        "criterion": c_star * c_star < sigma2 / (eps_param * eps_param),
    }


def confounding_sweep():
    """D1-D4 with exact closed-form risk values."""
    eps_values = [0.2, 0.5, 1.0, 2.0, 5.0]
    rows = [confounding_case(e) for e in eps_values]
    # Crossover: sigma2*(1+e^-2) = sigma2 + c*^2  =>  e^2 = sigma2/c*^2.
    crossover = float(np.sqrt(1.0 / 1.0))
    # Map-observation example (D3).
    map_obs = {
        "y3_model": "y3 = chi + eps3, eps3 ~ N(0, sigma2)",
        "map_information_retention": 1.0,
        "pose_truncation_bias": 1.0,
        "note": (
            "An independent map observation gives unit relative map retention "
            "(map is not confounded), while the truncated pose estimator is "
            "still biased by c*=1.  A relative map-retention gate would pass "
            "despite a large pose error."
        ),
    }
    # Small residual / large bias control (D4).
    e4 = 1e-3
    s4 = 1e-6
    c4 = 1.0
    full4 = s4 * (1.0 + e4 ** (-2))
    trunc4 = s4 + c4 * c4
    small_large = {
        "sigma2": s4,
        "eps_param": e4,
        "c_star": c4,
        "full_ls_risk_exact": full4,
        "truncated_mse_exact": trunc4,
        "truncation_bias": c4,
        "residual_y2_under_truncation": abs(e4 * c4),
        "residual_y2_relative_to_pose_error": abs(e4 * c4) / c4,
        "residual_scale_relative_to_sigma": abs(e4 * c4) / np.sqrt(s4),
        "note": (
            "The fitted y2 residual under the truncated (c=0) estimator is "
            "only eps_param*c_star = 1e-3 (0.1% of the pose error c_star=1.0), "
            "so a residual-only gate passes although the pose estimate is "
            "wrong by c_star=1.0. "
            "Both full LS and truncated MSE equal 1.000001 because "
            "e=sqrt(sigma2)/c_star; the full-LS variance is e^-2*sigma2 = 1, "
            "i.e. full LS is only 'exact' relative to this inflated variance, "
            "not relative to sigma2=1e-6."
        ),
    }
    return {
        "model": "y1 = h + c + eps1, y2 = eps_param*c + eps2, eps_i~N(0,sigma2)",
        "full_ls_risk_formula": "sigma2*(1 + eps_param^-2)",
        "truncated_risk_formula": "sigma2 + c_star^2",
        "crossover_criterion": "truncation better iff c_star^2 < sigma2/eps_param^2",
        "sigma2": 1.0,
        "c_star": 1.0,
        "crossover_eps_param": crossover,
        "records": rows,
        "expected_better": [e < 1.0 for e in eps_values],
        "sweep_pass": all(
            row["truncation_better"] == (row["eps_param"] < 1.0)
            and row["full_ls_risk"] == 2.0 if row["eps_param"] == 1.0
            else True
            for row in rows
        ),
        "map_observation": map_obs,
        "small_residual_large_bias": small_large,
    }


def all_e3():
    return {
        "dual_spectra": dual_spectra_all(),
        "pose_prior": pose_prior_checks(),
        "theorem7": theorem7_checks(),
        "adaptive_rank_caveat": adaptive_rank_caveat(),
        "confounding": confounding_sweep(),
    }
