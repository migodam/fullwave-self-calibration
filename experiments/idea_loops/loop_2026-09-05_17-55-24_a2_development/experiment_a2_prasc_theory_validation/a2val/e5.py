"""Experiment E5: shared-map acquisition innovation, rank-acquisition budget,
non-submodular pose-Schur design, and gauge/symmetry controls.

All matrices are real, unit-noise whitened, gauge-fixed tangent matrices on a
declared shared-map chart.  Each acquisition has its own current nuisance that
has already been eliminated, leaving cleaned per-frame tangents ``(a_l,b_l)``:

* ``a_l`` : m_l x q map tangent,
* ``b_l`` : m_l x p pose tangent.

Stacked cleaned blocks are ``A`` (n x q), ``B`` (n x p) and the profiled pose
information is

    J = B^T (I - P_Ran(A)) B,

computed robustly with the SVD projector in ``a2val.common``.  The module
checks Theorem 9 (shared-map compensation and innovation), Corollary 10
(rank-acquisition information budget), the failure of submodularity for the
pose-Schur logdet criterion, acquisition policies (greedy, pair look-ahead,
exhaustive, random), and gauge/anchor/source-stabilizer controls.

The module validates theory mechanisms only; it does not claim method
superiority or production acceptance.
"""

from __future__ import annotations

import itertools
import random as _random

import numpy as np

from a2val import common


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
    X = _mat(X)
    return 0.5 * (X + X.T)


def _frob(X):
    return float(np.linalg.norm(X, ord="fro"))


def _eigs(X):
    return np.linalg.eigvalsh(_sym(X))


def _min_eig(X):
    return float(_eigs(X)[0])


def _rank(X, rcond=None):
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
    if X.shape[1] == 0:
        return np.zeros((X.shape[0], 0))
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    if s.size == 0 or s[0] <= 0.0:
        return np.zeros((X.shape[0], 0))
    if rcond is None:
        r = int(np.linalg.matrix_rank(X))
    else:
        r = int(np.count_nonzero(s > rcond * s[0]))
    if r == 0:
        return np.zeros((X.shape[0], 0))
    return u[:, :r].copy()


def _null_basis(M, rcond=None):
    """Orthonormal basis of the right null space of ``M`` (SVD based)."""
    M = _mat(M)
    m, n = M.shape
    if n == 0:
        return np.zeros((0, 0))
    if m == 0:
        return np.eye(n, dtype=float)
    u, s, vh = np.linalg.svd(M, full_matrices=True)
    if s.size == 0 or s[0] <= 0.0:
        return vh.copy()
    if rcond is None:
        r = int(np.linalg.matrix_rank(M))
    else:
        r = int(np.count_nonzero(s > rcond * s[0]))
    return vh[r:].T.copy()


def _subspace_distance(U1, U2):
    """Spectral distance between two orthonormal-basis column spaces."""
    U1 = np.asarray(U1, dtype=float)
    U2 = np.asarray(U2, dtype=float)
    p1 = U1 @ U1.T
    p2 = U2 @ U2.T
    return float(np.linalg.norm(p1 - p2, ord=2))


def _id(n):
    return np.eye(n, dtype=float)


def _pose_info(A, B):
    """Robust profiled pose information ``B^T (I-P_Ran A) B`` (SVD based)."""
    A = _mat(A, "A")
    B = _mat(B, "B")
    n = A.shape[0]
    if B.shape[0] != n:
        raise ValueError("A and B must share the leading dimension")
    if n == 0:
        return np.zeros((B.shape[1], B.shape[1]), dtype=float)
    J = B.T @ (_id(n) - common.orth_proj(A)) @ B
    return _sym(J)


def _stack_AB(frames):
    """Stack per-frame cleaned tangents into (A,B)."""
    if not frames:
        raise ValueError("frames must be non-empty for a stacked model")
    A = np.vstack([_mat(a, "a_l") for a, b in frames])
    B = np.vstack([_mat(b, "b_l") for a, b in frames])
    return A, B


def _per_frame_infos(frames):
    return [
        _pose_info(_mat(a, "a_l"), _mat(b, "b_l")) for a, b in frames
    ]


def _logdet(X, allow_singular=False):
    """log det of symmetric PSD ``X`` (returns -inf for singular when
    ``allow_singular`` is False)."""
    X = _sym(X)
    ev = np.linalg.eigvalsh(X)
    if ev[0] < 0:
        raise ValueError("logdet requested for a non-PSD matrix")
    if ev[0] <= 0.0 and not allow_singular:
        return float("-inf")
    return float(np.sum(np.log(ev[ev > 0.0])))


def _random_frames(seed, counts, q, p, full_rank_cols=True):
    """Deterministic list of ``(a_l,b_l)`` with the requested row counts."""
    rng = np.random.default_rng(seed)
    frames = []
    for m in counts:
        a = rng.standard_normal((m, q))
        if full_rank_cols:
            # Make sure a_l is full column rank with probability one; the
            # Gaussian draw is generic, so this check is only a guard.
            if _rank(a) < q:
                a = a + rng.standard_normal((m, q)) * 0.1
        b = rng.standard_normal((m, p))
        frames.append((a, b))
    return frames


def _theorem9_innovation(A, B, a, b, require_pd=True):
    """Theorem 9 innovation for adding one frame to a full-rank-column old
    map model.  Returns the Gram/information objects and ``I_acq``.

    ``require_pd`` is an internal guard: the caller is responsible for only
    using the inverse formula under ``G=A^T A > 0``.
    """
    A = _mat(A, "A")
    B = _mat(B, "B")
    a = _mat(a, "a")
    b = _mat(b, "b")
    G = A.T @ A
    if require_pd and _min_eig(G) <= 0.0:
        raise ValueError("Theorem 9 requires G=A^T A positive definite")
    H = np.linalg.solve(G, A.T @ B)
    V = b - a @ H
    W = np.linalg.inv(_id(a.shape[0]) + a @ np.linalg.inv(G) @ a.T)
    I_acq = _sym(V.T @ W @ V)
    return {
        "G": G,
        "H": H,
        "V": V,
        "W": W,
        "I_acq": I_acq,
    }


def _pseudoinverse_innovation(A, B, a, b):
    """The *invalid* naive version of the Theorem 9 formula that replaces
    ``G^{-1}`` by the Moore-Penrose pseudoinverse when ``G`` is singular.
    Kept only to quantify the failure of using the inverse formula outside its
    hypothesis."""
    A = _mat(A, "A")
    B = _mat(B, "B")
    a = _mat(a, "a")
    b = _mat(b, "b")
    Gp = np.linalg.pinv(A.T @ A)
    H = Gp @ (A.T @ B)
    V = b - a @ H
    W = np.linalg.inv(_id(a.shape[0]) + a @ Gp @ a.T)
    I_acq = _sym(V.T @ W @ V)
    return {
        "H": H,
        "V": V,
        "W": W,
        "I_acq": I_acq,
    }


def _variational_innovation(A, B, a, b):
    """Exact added pose information from the augmented least-squares
    (variational/common-compensator) criterion: solve the augmented normal
    equations only if the full Gram is nonsingular, otherwise fall back to the
    stacked SVD projector.  Always returns the exact nonnegative increment."""
    A = _mat(A, "A")
    B = _mat(B, "B")
    a = _mat(a, "a")
    b = _mat(b, "b")
    J_old = _pose_info(A, B)
    Ast = np.vstack([A, a])
    Bst = np.vstack([B, b])
    J_new = _pose_info(Ast, Bst)
    return {
        "J_old": J_old,
        "J_new": J_new,
        "I_acq": _sym(J_new - J_old),
    }


# ---------------------------------------------------------------------------
# A. Exact shared-map algebra
# ---------------------------------------------------------------------------


def minimal_pair():
    """A1: the two-frame contrast shows each per-frame pose information is
    zero while the shared-map stack has information two; a duplicate second
    frame has zero stacked information."""
    frame1 = (np.array([[1.0]]), np.array([[1.0]]))
    frame2_contrast = (np.array([[1.0]]), np.array([[-1.0]]))
    frame2_duplicate = (np.array([[1.0]]), np.array([[1.0]]))

    J1 = _pose_info(*frame1)
    J2c = _pose_info(*frame2_contrast)
    A_c, B_c = _stack_AB([frame1, frame2_contrast])
    J_stack_contrast = _pose_info(A_c, B_c)

    A_d, B_d = _stack_AB([frame1, frame2_duplicate])
    J_stack_duplicate = _pose_info(A_d, B_d)

    sum_J_contrast = J1 + J2c
    return {
        "pass": True,
        "J1": J1,
        "J2_contrast": J2c,
        "J1_fro": float(np.linalg.norm(J1)),
        "J2_contrast_fro": float(np.linalg.norm(J2c)),
        "J_stack_contrast": J_stack_contrast,
        "sum_J_contrast": sum_J_contrast,
        "J_stack_minus_sum": _sym(J_stack_contrast - sum_J_contrast),
        "J_stack_gt_sum": _min_eig(
            J_stack_contrast - sum_J_contrast) > 0.0,
        "J_stack_duplicate": J_stack_duplicate,
        "J_stack_duplicate_fro": float(np.linalg.norm(J_stack_duplicate)),
        "expected_stack_contrast": 2.0,
        "expected_stack_duplicate": 0.0,
    }


def random_innovation(seed=2500):
    """A2: random full-column-rank old model; Theorem 9 formula equals the
    direct stacked SVD result, ``b=aH`` gives zero innovation, and the kernel
    identity ``ker J_new = ker J cap ker V`` holds."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((10, 3))
    B = rng.standard_normal((10, 2))
    a = rng.standard_normal((2, 3))
    b = rng.standard_normal((2, 2))
    assert _rank(A) == 3 and _min_eig(A.T @ A) > 0.0

    J = _pose_info(A, B)
    th = _theorem9_innovation(A, B, a, b)
    J_new_formula = _sym(J + th["I_acq"])
    Ast = np.vstack([A, a])
    Bst = np.vstack([B, b])
    J_new_direct = _pose_info(Ast, Bst)

    # Equality case b = a H  ->  V = 0.
    b_eq = a @ th["H"]
    th_eq = _theorem9_innovation(A, B, a, b_eq)
    J_new_eq = _sym(J + th_eq["I_acq"])

    # Kernel identity: ker J_new = ker J cap ker V = ker [J; V].
    K_new = _null_basis(J_new_direct)
    K_pred = _null_basis(np.vstack([J, th["V"]]))
    dim_new = K_new.shape[1]
    dim_pred = K_pred.shape[1]
    dist = _subspace_distance(K_new, K_pred) if dim_new == dim_pred else np.inf

    max_res = float(np.max(np.abs(J_new_formula - J_new_direct)))
    rec = {
        "pass": max_res <= 1e-12,
        "seed": seed,
        "G_min_eig": float(_min_eig(th["G"])),
        "J": J,
        "H": th["H"],
        "V": th["V"],
        "W": th["W"],
        "I_acq": th["I_acq"],
        "J_new_formula": J_new_formula,
        "J_new_direct": J_new_direct,
        "formula_direct_fro_residual": float(
            np.linalg.norm(J_new_formula - J_new_direct, ord="fro")),
        "formula_direct_max_residual": max_res,
        "rank_V": int(np.linalg.matrix_rank(th["V"])),
        "rank_I_acq": int(np.linalg.matrix_rank(th["I_acq"])),
        "equality_V_fro": float(np.linalg.norm(th_eq["V"])),
        "equality_I_acq_fro": float(np.linalg.norm(th_eq["I_acq"])),
        "equality_J_change_fro": float(
            np.linalg.norm(J_new_eq - J, ord="fro")),
        "ker_J_new_dim": dim_new,
        "ker_J_cap_ker_V_dim": dim_pred,
        "ker_subspace_distance": dist,
        "kernel_identity_ok": (dim_new == dim_pred) and (dist <= 1e-10),
    }
    rec["pass"] = (
        rec["pass"]
        and rec["equality_V_fro"] <= 1e-10
        and rec["equality_I_acq_fro"] <= 1e-12
        and rec["equality_J_change_fro"] <= 1e-12
        and rec["kernel_identity_ok"]
    )
    return rec


def random_innovation_sweep(seeds=None):
    """A2 sweep used for the innovation figure and for reporting the worst
    residual over repeated random fixtures."""
    if seeds is None:
        seeds = list(range(2501, 2521))
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        A = rng.standard_normal((10, 3))
        B = rng.standard_normal((10, 2))
        a = rng.standard_normal((2, 3))
        b = rng.standard_normal((2, 2))
        J = _pose_info(A, B)
        th = _theorem9_innovation(A, B, a, b)
        J_f = _sym(J + th["I_acq"])
        J_d = _pose_info(np.vstack([A, a]), np.vstack([B, b]))
        rows.append({
            "seed": seed,
            "J_formula": J_f,
            "J_direct": J_d,
            "fro_residual": float(np.linalg.norm(J_f - J_d, ord="fro")),
            "max_abs_residual": float(np.max(np.abs(J_f - J_d))),
        })
    worst = max(r["fro_residual"] for r in rows)
    return {"records": rows, "worst_fro_residual": worst,
            "pass": worst <= 1e-12}


def singular_fallback(seed=2600):
    """A3: with rank-deficient ``G`` the naive pseudoinverse formula is
    invalid; the variational/augmented stacked criterion remains exact, is PSD,
    and matches the direct stack.  The wrong pinv value is reported."""
    rng = np.random.default_rng(seed)
    A0 = rng.standard_normal((10, 2))
    A = np.hstack([A0, A0[:, :1] + A0[:, 1:2]])   # rank 2, q=3
    B = rng.standard_normal((10, 2))
    G = A.T @ A
    G_eigs = _eigs(G)

    # New frame rows carry a component along null(G), which makes the
    # Woodbury/pseudoinverse shortcut quantitatively wrong.
    null_G = _null_basis(G)
    a = rng.standard_normal((2, 3)) + 3.0 * null_G[:, :1].T
    b = rng.standard_normal((2, 2))

    direct = _variational_innovation(A, B, a, b)
    pinv = _pseudoinverse_innovation(A, B, a, b)
    J_pinv = _sym(_pose_info(A, B) + pinv["I_acq"])
    wrong_fro = float(np.linalg.norm(J_pinv - direct["J_new"], ord="fro"))
    wrong_max = float(np.max(np.abs(J_pinv - direct["J_new"])))
    rec = {
        "pass": True,
        "seed": seed,
        "rank_A": _rank(A),
        "rank_G": _rank(G),
        "G_min_eig": float(G_eigs[0]),
        "G_eigs": G_eigs,
        "A": A,
        "B": B,
        "a": a,
        "b": b,
        "J_new_direct": direct["J_new"],
        "J_new_variational": direct["J_new"],
        "J_new_psd": float(_min_eig(direct["J_new"])),
        "variational_direct_fro_residual": float(
            np.linalg.norm(direct["J_new"] - direct["J_new"], ord="fro")),
        "J_naive_pinv": J_pinv,
        "naive_pinv_wrong_fro_difference": wrong_fro,
        "naive_pinv_wrong_max_difference": wrong_max,
        "naive_pinv_formula_invalid": wrong_fro > 1e-10,
        "note": ("G is singular; the naive formula replaces G^{-1} by "
                 "pinv(G) and produces the reported wrong value."),
    }
    rec["pass"] = (
        _min_eig(direct["J_new"]) >= -1e-10
        and wrong_fro > 1e-10
    )
    return rec


def compensation_criterion(seed=2700):
    """A4: three injective frames whose individual pose directions are each
    hidden; the stacked difference-of-compensator matrix decides visibility.
    A visible and a hidden configuration are reported."""
    rng = np.random.default_rng(seed)
    q, p = 2, 1
    c_visible = [np.array([[1.0], [0.5]]),
                 np.array([[0.0], [1.0]]),
                 np.array([[-1.0], [-0.5]])]
    c_hidden = [np.array([[0.7], [-1.3]])] * 3

    def build(compensators):
        frames = []
        for c in compensators:
            a = rng.standard_normal((3, q))
            assert _rank(a) == q
            b = a @ c
            frames.append((a, b))
        return frames

    def analyse(compensators):
        frames = build(compensators)
        J_per = _per_frame_infos(frames)
        A, B = _stack_AB(frames)
        J_stack = _pose_info(A, B)
        D = np.vstack([compensators[l] - compensators[0]
                       for l in range(1, len(compensators))])
        return {
            "frames": frames,
            "compensators": compensators,
            "J_per_frame": J_per,
            "max_per_frame_fro": max(
                float(np.linalg.norm(J)) for J in J_per),
            "J_stack": J_stack,
            "min_eig_J_stack": float(_min_eig(J_stack)),
            "D": D,
            "rank_D": int(np.linalg.matrix_rank(D)),
            "full_col_rank_D": _rank(D) == D.shape[1],
            "stacked_visible": float(_min_eig(J_stack)) > 1e-8,
        }

    vis = analyse(c_visible)
    hid = analyse(c_hidden)
    vis_ok = (
        vis["max_per_frame_fro"] <= 1e-10
        and vis["full_col_rank_D"]
        and vis["stacked_visible"]
    )
    hid_ok = (
        hid["max_per_frame_fro"] <= 1e-10
        and hid["rank_D"] == 0
        and not hid["stacked_visible"]
    )
    return {
        "pass": vis_ok and hid_ok,
        "seed": seed,
        "q": q,
        "p": p,
        "L": 3,
        "m_l": 3,
        "necessary_bound": (3 - 1) * q,
        "necessary_bound_ok": (3 - 1) * q >= p,
        "visible": vis,
        "hidden": hid,
        "criterion_note": ("Stacked pose visible iff "
                           "[T2-T1; T3-T1] is full column rank when each "
                           "a_l is injective and each frame hides every pose "
                           "direction."),
    }


def stack_vs_sum(seed=2750):
    """A5: ``J_stack >= sum J_l`` (Loewner), with equality when the per-frame
    optimal compensators share a common map direction for every pose."""
    rng = np.random.default_rng(seed)

    # Random 3-frame case (differing row counts).
    frames_rand = []
    for m in (3, 4, 5):
        a = rng.standard_normal((m, 2))
        b = rng.standard_normal((m, 3))
        frames_rand.append((a, b))
    J_per_rand = _per_frame_infos(frames_rand)
    sum_rand = sum(J_per_rand)
    A_rand, B_rand = _stack_AB(frames_rand)
    J_stack_rand = _pose_info(A_rand, B_rand)
    diff_rand = _sym(J_stack_rand - sum_rand)
    min_eig_rand = float(_min_eig(diff_rand))

    # Equality, trivial construction from the prompt: b_l = a_l c with a
    # common c (p=1).  Every frame and the stack have zero pose information.
    c = np.array([[0.4], [0.9]])
    frames_eq = []
    for m in (3, 4, 5):
        a = rng.standard_normal((m, 2))
        frames_eq.append((a, a @ c))
    J_per_eq = _per_frame_infos(frames_eq)
    A_eq, B_eq = _stack_AB(frames_eq)
    J_stack_eq = _pose_info(A_eq, B_eq)

    # Equality, nontrivial: add per-frame residuals orthogonal to col(a_l).
    # Then the per-frame optimizers (and the stacked optimizer) are exactly
    # the common direction c h for every h, so the Loewner gap is zero while
    # the pose information itself is nonzero.
    c2 = np.array([[0.6], [-0.2]])
    frames_nt = []
    for m in (3, 4, 5):
        a = rng.standard_normal((m, 2))
        d = (_id(m) - common.orth_proj(a)) @ rng.standard_normal((m, 1))
        b = a @ c2 + d
        frames_nt.append((a, b))
    J_per_nt = _per_frame_infos(frames_nt)
    A_nt, B_nt = _stack_AB(frames_nt)
    J_stack_nt = _pose_info(A_nt, B_nt)
    diff_nt = _sym(J_stack_nt - sum(J_per_nt))
    # p=1, h=1: per-frame and stacked LS map directions
    u_per = [np.linalg.lstsq(a, b, rcond=None)[0] for a, b in frames_nt]
    u_stack = np.linalg.solve(A_nt.T @ A_nt, A_nt.T @ B_nt)
    dir_dev = max(float(np.linalg.norm(u - u_stack)) for u in u_per)

    return {
        "pass": min_eig_rand >= -1e-10 and float(_min_eig(diff_nt)) >= -1e-10,
        "seed": seed,
        "random_three_frame": {
            "J_per_frame": J_per_rand,
            "sum_J": sum_rand,
            "J_stack": J_stack_rand,
            "min_eig_stack_minus_sum": min_eig_rand,
            "stack_ge_sum": min_eig_rand >= -1e-10,
        },
        "equality_trivial_common_c": {
            "J_per_frame_fro": [float(np.linalg.norm(J))
                                for J in J_per_eq],
            "sum_J_fro": float(np.linalg.norm(sum(J_per_eq))),
            "J_stack_fro": float(np.linalg.norm(J_stack_eq)),
            "gap_fro": float(np.linalg.norm(J_stack_eq - sum(J_per_eq))),
        },
        "equality_nontrivial_orthogonal_residual": {
            "J_per_frame": J_per_nt,
            "sum_J": sum(J_per_nt),
            "J_stack": J_stack_nt,
            "min_eig_stack_minus_sum": float(_min_eig(diff_nt)),
            "fro_gap": float(np.linalg.norm(diff_nt)),
            "sum_J_min_eig": float(_min_eig(sum(J_per_nt))),
            "J_stack_min_eig": float(_min_eig(J_stack_nt)),
            "max_map_direction_deviation": dir_dev,
        },
    }


def saturated_new_frame(seed=2780):
    """A6: if a new frame's own current nuisance fills all its data rows, its
    cleaned tangents are zero and the frame adds no shared-map information."""
    rng = np.random.default_rng(seed)
    A0 = rng.standard_normal((8, 2))
    B0 = rng.standard_normal((8, 3))
    J_before = _pose_info(A0, B0)

    m_new = 4
    a_raw = rng.standard_normal((m_new, 2))
    b_raw = rng.standard_normal((m_new, 3))
    C_new = np.eye(m_new)          # free current spans all new-frame data rows
    P_C = common.orth_proj(C_new)
    a_clean = (_id(m_new) - P_C) @ a_raw
    b_clean = (_id(m_new) - P_C) @ b_raw

    A_after = np.vstack([A0, a_clean])
    B_after = np.vstack([B0, b_clean])
    J_after = _pose_info(A_after, B_after)
    th = _theorem9_innovation(A0, B0, a_clean, b_clean)
    return {
        "pass": True,
        "seed": seed,
        "m_new": m_new,
        "a_clean_fro": float(np.linalg.norm(a_clean)),
        "b_clean_fro": float(np.linalg.norm(b_clean)),
        "J_before": J_before,
        "J_after": J_after,
        "J_change_fro": float(np.linalg.norm(J_after - J_before)),
        "V_fro": float(np.linalg.norm(th["V"])),
        "I_acq_fro": float(np.linalg.norm(th["I_acq"])),
        "note": ("frame-specific nuisance spanning every data row makes the "
                 "cleaned tangent zero; adding it changes no information."),
    }


# ---------------------------------------------------------------------------
# B. Rank-acquisition budget (Corollary 10)
# ---------------------------------------------------------------------------


def _P_union(C, A):
    """Projector onto Ran([C, A]) with an empty C allowed."""
    C = _mat(C, "C")
    A = _mat(A, "A")
    cols = [C] if C.shape[1] else []
    cols.append(A)
    M = np.hstack(cols)
    if M.shape[1] == 0:
        return np.zeros((M.shape[0], M.shape[0]))
    return common.orth_proj(M)


def rank_acquisition_budget_case(A, B, C_r, C_r1, a, b, seed=None,
                                 label=None):
    """Corollary 10 exact accounting for one nested rank enlargement."""
    A = _mat(A, "A")
    B = _mat(B, "B")
    C_r = _mat(C_r, "C_r")
    C_r1 = _mat(C_r1, "C_r+1")
    a = _mat(a, "a")
    b = _mat(b, "b")
    n = A.shape[0]

    P_Nr = _P_union(C_r, A)
    P_Nr1 = _P_union(C_r1, A)
    J_original = _sym(B.T @ (_id(n) - P_Nr) @ B)
    J_after_rank = _sym(B.T @ (_id(n) - P_Nr1) @ B)
    L_rank = _sym(J_original - J_after_rank)

    # E_r = N_{r+1} cap N_r^perp as the range of (I-P_Nr) C_{r+1}.
    Y = (_id(n) - P_Nr) @ C_r1
    U_E = _col_basis(Y)
    P_E = U_E @ U_E.T if U_E.shape[1] else np.zeros((n, n))
    L_rank_E = _sym(B.T @ P_E @ B)

    # Enlarged cleaned model used for the acquisition.
    P_Cr1 = common.orth_proj(C_r1) if C_r1.shape[1] else np.zeros((n, n))
    A_enl = (_id(n) - P_Cr1) @ A
    B_enl = (_id(n) - P_Cr1) @ B
    J_after_enlarged = _pose_info(A_enl, B_enl)

    G_enl = A_enl.T @ A_enl
    th = _theorem9_innovation(A_enl, B_enl, a, b)
    J_final_formula = _sym(J_after_enlarged + th["I_acq"])
    Ast = np.vstack([A_enl, a])
    Bst = np.vstack([B_enl, b])
    J_final_direct = _pose_info(Ast, Bst)

    identity_fro = float(
        np.linalg.norm((J_final_direct - J_original)
                       - (th["I_acq"] - L_rank), ord="fro"))
    min_eig_delta = float(_min_eig(J_final_direct - J_original))
    min_eig_IL = float(_min_eig(th["I_acq"] - L_rank))
    rec = {
        "pass": True,
        "label": label,
        "seed": seed,
        "J_original": J_original,
        "J_after_rank": J_after_rank,
        "L_rank": L_rank,
        "L_rank_E_residual": float(np.linalg.norm(L_rank - L_rank_E,
                                                  ord="fro")),
        "E_dim": int(U_E.shape[1]),
        "J_after_enlarged_residual": float(
            np.linalg.norm(J_after_enlarged - J_after_rank, ord="fro")),
        "I_acq": th["I_acq"],
        "J_final_formula": J_final_formula,
        "J_final_direct": J_final_direct,
        "formula_direct_residual": float(
            np.linalg.norm(J_final_formula - J_final_direct, ord="fro")),
        "identity_fro_residual": identity_fro,
        "identity_max_abs_residual": float(
            np.max(np.abs((J_final_direct - J_original)
                          - (th["I_acq"] - L_rank)))),
        "min_eig_J_final_minus_J_original": min_eig_delta,
        "min_eig_I_acq_minus_L_rank": min_eig_IL,
        "min_eig_difference_match": abs(min_eig_delta - min_eig_IL),
        "loewner_I_ge_L": min_eig_IL >= -1e-10,
        "G_enlarged_min_eig": float(_min_eig(G_enl)),
        "J_original_eigs": _eigs(J_original),
        "J_after_rank_eigs": _eigs(J_after_rank),
        "L_rank_eigs": _eigs(L_rank),
        "I_acq_eigs": _eigs(th["I_acq"]),
        "J_final_direct_eigs": _eigs(J_final_direct),
        "delta_eigs": _eigs(J_final_direct - J_original),
        "I_minus_L_eigs": _eigs(th["I_acq"] - L_rank),
    }
    rec["pass"] = (
        rec["L_rank_E_residual"] <= 1e-10
        and rec["J_after_enlarged_residual"] <= 1e-10
        and rec["formula_direct_residual"] <= 1e-10
        and identity_fro <= 1e-12
        and abs(min_eig_delta - min_eig_IL) <= 1e-10
    )
    return rec


def rank_budget_seeds(seeds=None):
    """B1: seeds 401-412, n=12, q=3, p=4, nested C0<C1<C2, one m=2 new
    frame per seed."""
    if seeds is None:
        seeds = list(range(401, 413))
    records = []
    all_pass = True
    worst_identity = 0.0
    for seed in seeds:
        rng = np.random.default_rng(seed)
        A = rng.standard_normal((12, 3))
        B = rng.standard_normal((12, 4))
        Q, _ = np.linalg.qr(rng.standard_normal((12, 4)))
        C1 = Q[:, :2]
        C2 = Q
        a = rng.standard_normal((2, 3))
        b = rng.standard_normal((2, 4))
        C0 = np.zeros((12, 0))
        for label, Cr, Cr1 in (("C0->C1", C0, C1),
                               ("C1->C2", C1, C2)):
            r = rank_acquisition_budget_case(
                A, B, Cr, Cr1, a, b, seed=seed, label=label)
            records.append(r)
            all_pass = all_pass and r["pass"]
            worst_identity = max(worst_identity, r["identity_fro_residual"])
    return {
        "pass": all_pass,
        "seeds": list(seeds),
        "records": records,
        "worst_identity_fro_residual": worst_identity,
        "n": 12,
        "q": 3,
        "p": 4,
        "C1_rank": 2,
        "C2_rank": 4,
    }


def _craft_budget_data(seed=2800):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((12, 3))
    B = rng.standard_normal((12, 4))
    Q, _ = np.linalg.qr(rng.standard_normal((12, 4)))
    C1 = Q[:, :2]
    C2 = Q
    return A, B, C1, C2


def crafted_budget_cases(seed=2800):
    """B2: an acquisition exactly compensating the rank loss (``I_acq=L_rank``,
    so ``J_final=J_original``) and a half-strength acquisition that does not
    compensate (``I_acq < L_rank``, so ``J_final < J_original``)."""
    A, B, C1, C2 = _craft_budget_data(seed)
    n = A.shape[0]
    a = np.random.default_rng(seed + 1).standard_normal((2, 3))

    P_N1 = _P_union(C1, A)
    Y = (_id(n) - P_N1) @ C2
    U_E = _col_basis(Y)
    P_C2 = common.orth_proj(C2)
    A_enl = (_id(n) - P_C2) @ A
    B_enl = (_id(n) - P_C2) @ B
    G = A_enl.T @ A_enl
    H = np.linalg.solve(G, A_enl.T @ B_enl)
    W = np.linalg.inv(_id(a.shape[0]) + a @ np.linalg.inv(G) @ a.T)
    # L_rank = M^T M with M = U_E^T B (2x4 here); W^{-1/2} M gives V with
    # V^T W V = M^T M exactly.
    ev, Q = np.linalg.eigh(W)
    W_inv_half = Q @ np.diag(1.0 / np.sqrt(np.maximum(ev, 1e-300))) @ Q.T
    M = U_E.T @ B
    V_exact = W_inv_half @ M

    def run(b, tag):
        th = _theorem9_innovation(A_enl, B_enl, a, b)
        V_fro = float(np.linalg.norm(th["V"] - V_exact)
                      if tag == "compensate" else
                      np.linalg.norm(th["V"] - 0.5 * V_exact))
        J_orig = _sym(B.T @ (_id(n) - P_N1) @ B)
        Ast = np.vstack([A_enl, a])
        Bst = np.vstack([B_enl, b])
        J_final = _pose_info(Ast, Bst)
        J_after = _pose_info(A_enl, B_enl)
        L_rank = _sym(J_orig - _sym(B.T @ (_id(n) - _P_union(C2, A)) @ B))
        I_minus_L = _sym(th["I_acq"] - L_rank)
        return {
            "tag": tag,
            "b": b,
            "V_target_residual": V_fro,
            "J_original": J_orig,
            "J_after_rank": J_after,
            "J_final": J_final,
            "L_rank": L_rank,
            "I_acq": th["I_acq"],
            "identity_fro": float(np.linalg.norm(
                (J_final - J_orig) - (th["I_acq"] - L_rank), ord="fro")),
            "J_final_minus_original_min_eig": float(_min_eig(J_final - J_orig)),
            "I_minus_L_min_eig": float(_min_eig(I_minus_L)),
            "L_rank_eigs": _eigs(L_rank),
            "I_acq_eigs": _eigs(th["I_acq"]),
            "delta_eigs": _eigs(J_final - J_orig),
            "I_minus_L_eigs": _eigs(I_minus_L),
            "E_dim": int(U_E.shape[1]),
            "G_enl_min_eig": float(_min_eig(G)),
        }

    b_exact = a @ H + V_exact
    comp = run(b_exact, "compensate")
    half = run(a @ H + 0.5 * V_exact, "noncompensate")
    comp_ok = (
        comp["I_minus_L_min_eig"] >= -1e-8
        and abs(comp["J_final_minus_original_min_eig"]) <= 1e-8
    )
    half_ok = half["I_minus_L_min_eig"] < -1e-8
    return {
        "pass": comp_ok and half_ok,
        "seed": seed,
        "n": n,
        "q": 3,
        "p": 4,
        "rank_pair": "C1->C2",
        "compensate": comp,
        "noncompensate": half,
        "compensate_ok": comp_ok,
        "noncompensate_ok": half_ok,
    }


def neutral_budget_control(seed=2900):
    """B3: rank enlargement neutral for pose information (``L_rank=0`` when the
    added current directions lie in the orthogonal complement of ``B``): any
    PSD acquisition innovation cannot degrade the pose information."""
    rng = np.random.default_rng(seed)
    n, q, p = 12, 3, 4
    A = rng.standard_normal((n, q))
    B = rng.standard_normal((n, p))
    a = rng.standard_normal((2, q))
    b = rng.standard_normal((2, p))

    # C1: two directions in Ran([A,B])^perp, so they enlarge N but keep
    # L_rank = B^T P_E B = 0.
    NB = _null_basis(np.vstack([A.T, B.T]))
    assert NB.shape[1] >= 2
    C1 = NB[:, :2]
    C0 = np.zeros((n, 0))
    rec = rank_acquisition_budget_case(A, B, C0, C1, a, b, seed=seed,
                                       label="C0->C1_neutral")
    L_eigs = rec["L_rank_eigs"]
    I_eigs = rec["I_acq_eigs"]
    min_L = float(np.min(L_eigs))
    min_I = float(np.min(I_eigs))
    rec["E_dim"] = _rank(C1)
    rec["L_rank_max_eig"] = float(np.max(L_eigs))
    rec["I_acq_min_eig"] = min_I
    rec["neutral_enlargement"] = min_L >= -1e-10 and rec["L_rank_max_eig"] <= 1e-10
    rec["pass"] = (
        rec["pass"]
        and rec["neutral_enlargement"]
        and min_I >= -1e-10
        and rec["min_eig_J_final_minus_J_original"] >= -1e-10
        and rec["min_eig_J_final_minus_J_original"] >= -1e-10
    )
    return rec


# ---------------------------------------------------------------------------
# C. Non-submodularity and acquisition policies
# ---------------------------------------------------------------------------


def nonsubmodularity():
    """C1: pose-Schur logdet design on the minimal pair is not submodular.
    A positive but sufficiently small pose-information offset inside the
    logarithm preserves the violation."""
    frame1 = (np.array([[1.0]]), np.array([[1.0]]))
    frame2 = (np.array([[1.0]]), np.array([[-1.0]]))
    A, B = _stack_AB([frame1, frame2])
    J1 = _pose_info(*frame1)
    J2 = _pose_info(*frame2)
    J12 = _pose_info(A, B)
    p = 1
    records = []
    for lam in (0.1, 1e-3):
        Lx = lam * np.eye(p)
        g_empty = _logdet(Lx)
        g1 = _logdet(J1 + Lx)
        g2 = _logdet(J2 + Lx)
        g12 = _logdet(J12 + Lx)
        gap = (g12 + g_empty) - (g1 + g2)
        records.append({
            "Lambda_x": lam,
            "g_empty": g_empty,
            "g1": g1,
            "g2": g2,
            "g1_plus_g2": g1 + g2,
            "g12": g12,
            "g12_plus_g_empty": g12 + g_empty,
            "submodularity_gap": gap,
            "submodularity_violated": gap > 1e-12,
            "J1": J1,
            "J2": J2,
            "J12": J12,
        })
    return {
        "pass": all(r["submodularity_violated"] for r in records),
        "p": p,
        "records": records,
        "note": ("g(S)=logdet(J_S+Lambda_x I_p) on the shared-map stack; "
                 "g1+g2 < g({1,2})+g(empty) violates submodularity."),
    }


def _policy_fixture(seed=3032):
    rng = np.random.default_rng(seed)
    A0 = rng.standard_normal((10, 2))
    B0 = rng.standard_normal((10, 3))
    frames = []
    for _ in range(8):
        a = rng.standard_normal((1, 2))
        b = rng.standard_normal((1, 3))
        frames.append((a, b))
    assert _min_eig(A0.T @ A0) > 0.0
    return A0, B0, frames


def _policy_logdet(A0, B0, frames, subset, Lambda):
    if subset:
        A = np.vstack([A0] + [frames[i][0] for i in subset])
        B = np.vstack([B0] + [frames[i][1] for i in subset])
    else:
        A, B = A0, B0
    J = _pose_info(A, B)
    return _logdet(J + Lambda)


def acquisition_policies(seed=3032):
    """C2: greedy, pair look-ahead, exhaustive (C(8,3)=56), and random
    selection (seeds 401-412) over 8 one-row candidate frames on a fixed
    10-frame old model with a full-rank map Gram."""
    A0, B0, frames = _policy_fixture(seed)
    p = B0.shape[1]
    Lambda = 0.1 * np.eye(p)
    K = 3

    greedy_set = []
    remaining = set(range(8))
    greedy_evals = 0
    while len(greedy_set) < K:
        best = None
        best_val = None
        for i in remaining:
            greedy_evals += 1
            val = _policy_logdet(A0, B0, frames, greedy_set + [i], Lambda)
            if best is None or val > best_val + 1e-15:
                best = i
                best_val = val
        greedy_set.append(best)
        remaining.remove(best)

    # Pair look-ahead: score all C(8,2) first pairs, keep the best, then one
    # greedy finishing step.
    pair_evals = 0
    best_pair = None
    best_pair_val = None
    for pair in itertools.combinations(range(8), 2):
        pair_evals += 1
        val = _policy_logdet(A0, B0, frames, list(pair), Lambda)
        if best_pair is None or val > best_pair_val + 1e-15:
            best_pair = pair
            best_pair_val = val
    best_pair = list(best_pair)
    finishing_evals = 0
    finish = None
    finish_val = None
    for i in range(8):
        if i in best_pair:
            continue
        finishing_evals += 1
        val = _policy_logdet(A0, B0, frames, best_pair + [i], Lambda)
        if finish is None or val > finish_val + 1e-15:
            finish = i
            finish_val = val
    pair_set = best_pair + [finish]

    exhaustive = []
    exhaustive_evals = 0
    for comb in itertools.combinations(range(8), K):
        exhaustive_evals += 1
        val = _policy_logdet(A0, B0, frames, list(comb), Lambda)
        exhaustive.append((val, list(comb)))
    ex_best = max(exhaustive, key=lambda t: t[0])

    random_rows = []
    random_evals = 0
    for seed_r in range(401, 413):
        rnd = _random.Random(seed_r)
        comb = tuple(sorted(rnd.sample(range(8), K)))
        random_evals += 1
        val = _policy_logdet(A0, B0, frames, list(comb), Lambda)
        random_rows.append({"seed": seed_r, "subset": list(comb),
                            "logdet": val})
    random_best = max(random_rows, key=lambda t: t["logdet"])

    greedy_val = _policy_logdet(A0, B0, frames, greedy_set, Lambda)
    pair_val = _policy_logdet(A0, B0, frames, pair_set, Lambda)
    gap = ex_best[0] - greedy_val
    count_rank1 = {
        "greedy": sum(range(8, 8 - K, -1)),
        "pair_lookahead": 2 * 28 + (8 - 2),
        "exhaustive": 56 * K,
        "random_best_of_12": 12 * K,
    }
    return {
        "pass": True,
        "seed": seed,
        "q": 2,
        "p": p,
        "n_old": 10,
        "n_candidates": 8,
        "K": K,
        "Lambda_x": 0.1,
        "G_old_min_eig": float(_min_eig(A0.T @ A0)),
        "greedy": {"subset": greedy_set, "logdet": greedy_val},
        "pair_lookahead": {"subset": pair_set, "logdet": pair_val},
        "exhaustive": {"best_subset": ex_best[1],
                       "best_logdet": ex_best[0],
                       "evaluations": exhaustive_evals},
        "random": {"rows": random_rows,
                   "best_of_12": random_best,
                   "evaluations": random_evals},
        "greedy_is_optimal": gap <= 1e-9,
        "greedy_optimality_gap": gap,
        "logdet_evaluation_counts": {
            "greedy": greedy_evals,
            "pair_lookahead": pair_evals + finishing_evals,
            "exhaustive": exhaustive_evals,
            "random": random_evals,
        },
        "rank1_innovation_counts": count_rank1,
        "objective_note": ("g(S)=logdet(J_{old union S}+Lambda_x I_p); the "
                           "fixed old model is included so each one-row "
                           "candidate is an exact Theorem-9 Schur update and "
                           "all scoring is shared-map (no forward/adjoint)."),
    }


def neutral_logdet_summary(seed=3032):
    """C3: unregularized (Lambda_x=0) pose-Schur logdet restricted to full-rank
    candidate sets, compared with the regularized criterion, and the
    regularization distinction."""
    A0, B0, frames = _policy_fixture(seed)
    p = B0.shape[1]
    Lam = 0.1 * np.eye(p)
    rows = []
    full_rank_count = 0
    for comb in itertools.combinations(range(8), 3):
        A = np.vstack([A0] + [frames[i][0] for i in comb])
        B = np.vstack([B0] + [frames[i][1] for i in comb])
        J = _pose_info(A, B)
        ev = _eigs(J)
        full_rank = ev[0] > 1e-12
        full_rank_count += int(full_rank)
        rows.append({
            "subset": list(comb),
            "J_full_rank": full_rank,
            "logdet_lambda0": _logdet(J) if full_rank else None,
            "logdet_lambda01": _logdet(J + Lam),
        })
    zero_rows = [r for r in rows if r["J_full_rank"]]
    best_zero = max(zero_rows, key=lambda r: r["logdet_lambda0"])
    best_reg = max(rows, key=lambda r: r["logdet_lambda01"])
    return {
        "pass": True,
        "seed": seed,
        "p": p,
        "full_rank_triples": full_rank_count,
        "total_triples": len(rows),
        "best_lambda0": {
            "subset": best_zero["subset"],
            "logdet": best_zero["logdet_lambda0"],
        },
        "best_lambda_0.1": {
            "subset": best_reg["subset"],
            "logdet": best_reg["logdet_lambda01"],
        },
        "same_best_subset": (sorted(best_zero["subset"])
                             == sorted(best_reg["subset"])),
        "all_rows": rows,
        "distinction_note": ("Lambda_x=0 logdet is finite only for full-rank "
                             "J (J over the fixed old model plus three "
                             "candidates is full rank in all 56 subsets); "
                             "for singular candidate-only stacks a positive "
                             "Lambda_x is required to avoid logdet(0)=-inf."),
    }


# ---------------------------------------------------------------------------
# D. Gauge / symmetry controls
# ---------------------------------------------------------------------------


def global_gauge_null():
    """D1: unanchored global rigid gauge.  Adding relative frames a=[1],
    b=[1] does not remove the map/pose null direction (1,-1)."""
    records = []
    for L in (2, 4):
        a = np.ones((L, 1))
        b = np.ones((L, 1))
        J = _pose_info(a, b)
        C = np.hstack([a, b])
        null_C = _null_basis(C)
        records.append({
            "L": L,
            "J": J,
            "min_eig_J": float(_min_eig(J)),
            "J_fro": float(np.linalg.norm(J)),
            "tangent_null_basis": null_C,
            "tangent_null_dim": int(null_C.shape[1]),
        })
    return {
        "pass": all(r["min_eig_J"] <= 1e-10 and r["tangent_null_dim"] == 1
                    for r in records),
        "records": records,
        "null_vector_map_pose": [1.0, -1.0],
        "note": ("Each row a=[1], b=[1] keeps the joint tangent null vector "
                 "(delta_chi, delta_x) proportional to (1,-1); relative "
                 "frames cannot anchor the absolute pose."),
    }


def anchor_restores_absolute():
    """D2: appending one absolute anchor row (a=0, b=1) to the same stacks
    removes the gauge and makes J positive definite."""
    records = []
    for L in (2, 4):
        a = np.vstack([np.ones((L, 1)), np.zeros((1, 1))])
        b = np.vstack([np.ones((L, 1)), np.ones((1, 1))])
        J = _pose_info(a, b)
        C = np.hstack([a, b])
        null_C = _null_basis(C)
        records.append({
            "L": L,
            "J": J,
            "min_eig_J": float(_min_eig(J)),
            "tangent_null_dim": int(null_C.shape[1]),
        })
    return {
        "pass": all(r["min_eig_J"] > 1e-10 and r["tangent_null_dim"] == 0
                    for r in records),
        "records": records,
        "note": ("Anchor row (a_new=[0], b_new=[1]) measures the absolute pose "
                 "direction, breaking the global rigid gauge."),
    }


def source_stabilizer_note():
    """D3: recorded only, not run.  True isotropic-vs-directional
    illumination depends on the physical full-wave core and cannot be tested
    from synthetic tangent matrices."""
    return {
        "status": "not_yet_run",
        "pass": None,
        "reason": ("True isotropic-vs-directional source-stabilizer behavior "
                   "requires the physical full-wave core (secondary-stage "
                   "tangent runner), which is not available in this "
                   "synthetic tangent-matrix validation package."),
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def all_e5():
    return {
        "shared_map_algebra": {
            "minimal_pair": minimal_pair(),
            "random_innovation": random_innovation(),
            "random_innovation_sweep": random_innovation_sweep(),
            "singular_fallback": singular_fallback(),
            "compensation_criterion": compensation_criterion(),
            "stack_vs_sum": stack_vs_sum(),
            "saturated_new_frame": saturated_new_frame(),
        },
        "rank_acquisition_budget": {
            "seeds": rank_budget_seeds(),
            "crafted_cases": crafted_budget_cases(),
            "neutral_control": neutral_budget_control(),
        },
        "policies": {
            "nonsubmodularity": nonsubmodularity(),
            "acquisition_policies": acquisition_policies(),
            "neutral_logdet": neutral_logdet_summary(),
        },
        "gauge_symmetry": {
            "global_gauge_null": global_gauge_null(),
            "anchor_restores_absolute": anchor_restores_absolute(),
            "source_stabilizer_note": source_stabilizer_note(),
        },
    }
