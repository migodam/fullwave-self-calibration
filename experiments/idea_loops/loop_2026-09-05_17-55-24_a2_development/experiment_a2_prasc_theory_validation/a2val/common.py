"""Shared linear-algebra helpers for the A2 validation package.

The helpers here are deliberately small; later experiments (E2/E3/E5) can add
their own helpers without modifying E1.
"""

from __future__ import annotations

import numpy as np


def _as_2d(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X)
    if X.ndim == 0:
        raise ValueError("X must have at least one dimension")
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    return X


def orth_proj(X, rcond=None):
    """Orthogonal projector onto the column range of ``X``.

    Projection is computed from the compact SVD ``X = U S V^H``:

        P_col = U[:, :r] U[:, :r]^H,

    with ``r`` the numerical rank (singular values above a tolerance).  The
    default tolerance follows ``numpy.linalg.pinv``: ``rcond * s[0]`` with
    ``rcond = max(shape) * eps``.  For complex ``X`` the projector is
    Hermitian (``P^H = P``).

    Parameters
    ----------
    X : (n, p) array_like
        Matrix whose column range is projected onto.
    rcond : float or None
        Relative singular-value cutoff; if None, the numpy pinv default is
        used.

    Returns
    -------
    P : (n, n) ndarray
        Orthogonal projector onto col-range(X).
    """
    X = _as_2d(X)
    n = X.shape[0]
    if n == 0:
        return np.zeros((0, 0), dtype=np.complex128 if np.iscomplexobj(X) else np.float64)
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    dtype = np.complex128 if np.iscomplexobj(X) else np.float64
    if s.size == 0 or s[0] <= 0.0:
        return np.zeros((n, n), dtype=dtype)
    if rcond is None:
        rcond = max(X.shape) * np.finfo(s.dtype).eps
    tol = rcond * s[0]
    rank = int(np.count_nonzero(s > tol))
    if rank == 0:
        return np.zeros((n, n), dtype=dtype)
    uu = u[:, :rank]
    return (uu @ uu.conj().T).astype(dtype, copy=False)


def proj_null(X, rcond=None):
    """Orthogonal projector onto the null space of ``X^H`` (left null range).

    Equivalently the complement of the column-range projector:

        P_null = I - P_col(X).

    Parameters are the same as :func:`orth_proj`.
    """
    n = _as_2d(X).shape[0]
    P = orth_proj(X, rcond=rcond)
    I = np.eye(n, dtype=P.dtype)
    return I - P


def numerical_rank(X, rcond=None):
    """Numerical rank of ``X`` under an SVD cutoff (helper for later use)."""
    X = _as_2d(X)
    s = np.linalg.svd(X, compute_uv=False)
    if s.size == 0 or s[0] <= 0.0:
        return 0
    if rcond is None:
        rcond = max(X.shape) * np.finfo(s.dtype).eps
    return int(np.count_nonzero(s > rcond * s[0]))

