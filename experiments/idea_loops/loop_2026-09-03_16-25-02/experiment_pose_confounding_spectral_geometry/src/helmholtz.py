"""2D scalar Helmholtz contrast-source forward model and analytic Jacobians.

Conventions implemented exactly (documented again because every numerical
claim depends on them):

Domain / grid
  D = [-0.5, 0.5]^2, world-fixed uniform N x N cell-center grid, cell side
  h = 1/N, cell area h^2.  Cell centers

      z_{i,j} = (-0.5 + (i + 0.5) h , -0.5 + (j + 0.5) h),

  flattened row-major with index n = i*N + j  (i = row = x index).

Green function (k_b real positive)
  g(r, r') = (i/4) H_0^{(1)}(k_b |r - r'|)
  grad_r g(r, r') = -(i k_b/4) H_1^{(1)}(k_b R) (r - r')/R   (R = |r - r'| > 0)
  and the gradient is the zero vector at R = 0 (never evaluated there).

Self-cell singular quadrature (domain propagator diagonal): integrate g over
the disk of area h^2, radius a = h/sqrt(pi):

      I_self = (i/4) 2 pi int_0^a r H_0^{(1)}(k_b r) dr
             = (i pi a / (2 k_b)) H_1^{(1)}(k_b a) - 1/k_b^2.

The lower-endpoint term -1/k_b^2 is required: with x = k_b r,
d[x H_1^{(1)}(x)]/dx = x H_0^{(1)}(x) and
lim_{x->0} x H_1^{(1)}(x) = -2i/pi, so the endpoint at r = 0 contributes
-1/k_b^2 after the 2 pi angular factor.  Without that term the purported
integral tends to 1/k_b^2 instead of zero as h -> 0 (h -> 0 means a -> 0).
The diagonal Green entry of G_D is k_b^2 * I_self (NOT multiplied again by
h^2; I_self already carries the cell-area factor).  Off-diagonal entries are
k_b^2 h^2 g(...).  This mirrors the stated rule
      G_D = k_b^2 * (Green_matrix_NxN * h^2),
      diagonal entries = k_b^2 * I_self.

Operators (dense complex128)
  G_D = k_b^2 h^2 * point-Green matrix, diagonal k_b^2 g_self   (N^2 x N^2)
  G_S = k_b^2 h^2 * Green matrix from grid cells to receivers   (n_rx x N^2)
  D_chi = diag(chi), chi in R^{N^2}
  M = I - D_chi G_D
  E_inc[n] = q g(z_n, s_tx),  q = 1.0
  J = solve(M, D_chi E_inc)
  E_tot = E_inc + G_D J
  F(chi) = G_S J

Matrix convention used everywhere in this module:
  green_matrix(src, dst) returns G[d, s] = g(dst_d, src_s), i.e. rows are the
  first argument ("field/observation" point) and columns the second argument
  ("source/integration" point).  G_S rows are therefore receivers r_a and
  columns grid centers z_n:  G_S[a, n] = k_b^2 h^2 g(r_a, z_n).

Full-wave map Jacobian (map chi -> data F, pose fixed)
  A_t = G_S solve(M, diag(E_tot))                (n_rx x N^2)
Born approximation
  A_born,t = G_S diag(E_inc);  F_born = G_S (chi * E_inc)

Pose parametrisation x_t = (p_x, p_y, theta)
  r_a = p + R(theta) s_rx[a],  s_tx = p + R(theta) s_tx_body,
  R(theta) = [[cos t, -sin t], [sin t, cos t]]
  d r_a / d px = (1,0), d r_a / d py = (0,1),
  d r_a / d theta = J (r_a - p),  J = [[0,-1],[1,0]]   (same for s_tx)

Pose Jacobian B_t (n_rx x 3), world-fixed grid so D_x G_D = 0:
  B_t[:, l] = (D_xG_S_l) J  +  G_S solve(M, D_chi (D_xE_inc_l))
  D_xG_S_l[a,n] = k_b^2 h^2 grad_1 g(r_a, z_n) . (d r_a / d x_l)
  D_xE_inc_l[n] = q grad_2 g(z_n, s_tx) . (d s_tx / d x_l)

SIGN NOTE (verified by the Family 1 finite-difference checks): the
mathematical derivative of the implemented incident field E_inc[n] =
q g(z_n, s_tx) with respect to the transmitter position s_tx is

      grad_s g(z, s) = +(i k_b/4) H_1^{(1)}(k_b R) (z - s)/R
                     = - grad_z g(z, s)
                     = - green_grad_first(src=s, dst=z).

Writing s = z - (z - s), d g(z, s)/ds = -d g(z, s)/dz, and the numerical
transmitter-position sweep confirms the plus sign (a relative discrepancy of
exactly 2 if the minus sign is used).  The code below evaluates grad_s with
the explicit helper green_grad_source, which applies that sign.

LU factorisation of M is reused across all right-hand sides that share M.
"""

from __future__ import annotations

import numpy as np
from functools import lru_cache
from scipy.linalg import lu_factor, lu_solve, svdvals
from scipy.special import hankel1

__all__ = [
    "make_grid",
    "self_cell_green",
    "green_matrix",
    "green_grad_first",
    "green_grad_source",
    "pose_geometry",
    "build_operators",
    "build_AB",
    "whiten_realify",
    "born_forward",
    "forward_measurements",
]

SELF_CELL_FORMULA = (
    "I_self = (i*pi*a/(2*k_b))*H_1^(1)(k_b*a) - 1/k_b^2, a = h/sqrt(pi)"
)
SELF_CELL_FORMULA_VERSION = (
    "self-cell v2 corrected 2026-09-03: complete equal-area disk integral "
    "(adds the -1/k_b^2 lower-endpoint term omitted by v1)"
)

_J90 = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=float)
_SQRT2 = np.sqrt(2.0)


# --------------------------------------------------------------------------
# Grid and Green functions
# --------------------------------------------------------------------------

def make_grid(N: int) -> tuple[np.ndarray, float]:
    """Return (points, h) for the D = [-0.5, 0.5]^2 N x N cell-center grid.

    points has shape (N*N, 2); index n = i*N + j holds the cell center
    (x_i, y_j) with x_i, y_j = -0.5 + (0..N-1 + 0.5) / N.
    """
    if not isinstance(N, (int, np.integer)) or int(N) < 1:
        raise ValueError(f"N must be a positive integer, got {N!r}")
    N = int(N)
    h = 1.0 / N
    centers = -0.5 + (np.arange(N, dtype=float) + 0.5) * h
    X, Y = np.meshgrid(centers, centers, indexing="ij")  # X = x index i, Y = y index j
    points = np.empty((N * N, 2), dtype=float)
    points[:, 0] = X.ravel(order="C")
    points[:, 1] = Y.ravel(order="C")
    return points, h


def self_cell_green(k_b: float, h: float) -> complex:
    """Disk-integrated self-cell Green value I_self (area units, no h^2 factor).

    Complete equal-area disk integral of g(r, r') = (i/4) H_0^{(1)}(k_b |r-r'|)
    over the disk of area h^2, radius a = h/sqrt(pi):

      I_self = (i/4) 2 pi int_0^a r H_0^{(1)}(k_b r) dr
             = (i pi a / (2 k_b)) H_1^{(1)}(k_b a) - 1/k_b^2.

    The final -1/k_b^2 comes from the lower endpoint x -> 0 of
    d[x H_1^{(1)}(x)]/dx = x H_0^{(1)}(x) with lim_{x->0} x H_1^{(1)}(x) =
    -2i/pi, multiplied by the 2 pi angular factor.  Omitting it would leave a
    cell integral that tends to 1/k_b^2 rather than zero as h -> 0.
    """
    k_b = float(k_b)
    a = float(h) / np.sqrt(np.pi)
    return (1j * np.pi * a / (2.0 * k_b)) * hankel1(1, k_b * a) - 1.0 / k_b**2


def _pair_diff_norm(src: np.ndarray, dst: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """diff[d, s, :] = dst[d] - src[s]; R[d, s] = |diff[d, s, :]|."""
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    diff = dst[:, None, :] - src[None, :, :]
    R = np.linalg.norm(diff, axis=2)
    return diff, R


def green_matrix(
    src: np.ndarray,
    dst: np.ndarray,
    k_b: float,
    include_self: bool = False,
    h: float | None = None,
) -> np.ndarray:
    """Point Green matrix G[d, s] = g(dst[d], src[s]).

    src: (N_s, 2) source/second-argument points (matrix columns).
    dst: (N_d, 2) field/first-argument points (matrix rows).

    Off-diagonal (R > 0) entries are the point values (i/4) H0(k_b R).
    Coincident pairs (R = 0) return 0, since g(0) is singular; with
    include_self=True they are replaced by the disk-integrated self-cell value
    g_self (requires h).  Note the self value already contains the cell-area
    factor, so callers that multiply off-diagonal point values by h^2 must not
    multiply the diagonal self entry by h^2 again.
    """
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    if src.shape[1] != 2 or dst.shape[1] != 2:
        raise ValueError("src and dst must have shape (..., 2)")
    _, R = _pair_diff_norm(src, dst)
    safe = R > 0.0
    G = np.zeros(R.shape, dtype=np.complex128)
    if safe.any():
        G[safe] = (1j / 4.0) * hankel1(0, float(k_b) * R[safe])
    if include_self:
        if h is None:
            raise ValueError("include_self=True requires the cell side h")
        G[~safe] = self_cell_green(k_b, h)
    return G


def green_grad_first(src: np.ndarray, dst: np.ndarray, k_b: float) -> np.ndarray:
    """Gradient of the Green function with respect to its first argument.

    Returns V[d, s, :] = grad_r g(r, r')|_(r = dst[d], r' = src[s])
                       = -(i k_b/4) H_1^{(1)}(k_b R) (dst[d] - src[s]) / R,
    and the zero vector wherever R = 0 (never evaluated there).

    In particular, with src = s and dst = z this function returns
    grad_z g(z, s) = -(i k_b/4) H_1^{(1)}(k_b R) (z - s)/R.  The derivative
    with respect to the source position is the opposite vector,
    grad_s g(z, s) = -grad_z g(z, s) = +(i k_b/4) H_1^{(1)}(k_b R) (z - s)/R;
    use green_grad_source for that sign (see SIGN NOTE in the module docstring).
    """
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    if src.shape[1] != 2 or dst.shape[1] != 2:
        raise ValueError("src and dst must have shape (..., 2)")
    diff, R = _pair_diff_norm(src, dst)
    V = np.zeros(diff.shape, dtype=np.complex128)
    safe = R > 0.0
    if safe.any():
        H = hankel1(1, float(k_b) * R[safe])
        V[safe] = (
            (-(1j * k_b) / 4.0)
            * H[:, None]
            * (diff[safe] / R[safe][:, None])
        )
    return V


def green_grad_source(z: np.ndarray, s: np.ndarray, k_b: float) -> np.ndarray:
    """Gradient of g(z, s) = (i/4) H_0^{(1)}(k_b |z - s|) wrt the source s.

    Returns V[n, :] = grad_s g(z_n, s)
                    = +(i k_b/4) H_1^{(1)}(k_b R) (z_n - s) / R
                    = -grad_z g(z_n, s)
                    = -green_grad_first(src=s, dst=z)[:, 0, :],
    with z an (N_z, 2) array of field points and s a single (2,) source point.
    The zero vector is returned wherever R = 0 (never evaluated there).
    """
    s = np.atleast_2d(np.asarray(s, dtype=float))
    z = np.atleast_2d(np.asarray(z, dtype=float))
    if z.shape[1] != 2 or s.shape[1] != 2:
        raise ValueError("z and s must have shape (..., 2)")
    return -green_grad_first(s, z, k_b)[:, 0, :]


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def pose_geometry(
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
) -> dict[str, np.ndarray]:
    """World positions and pose derivatives for transmitter/receiver offsets.

    poses: (T, 3), columns (p_x, p_y, theta).
    rx_offsets: (n_rx, 2) body coordinates of receivers.
    tx_offset: (2,) body coordinate of the transmitter.

    Returns dict with
      rx_world (T, n_rx, 2), tx_world (T, 2),
      drx_world (T, 3, n_rx, 2) derivative per pose component l = px, py, theta,
      dtx_world (T, 3, 2).
    """
    poses = np.atleast_2d(np.asarray(poses, dtype=float))
    rx_offsets = np.atleast_2d(np.asarray(rx_offsets, dtype=float))
    tx_offset = np.asarray(tx_offset, dtype=float).reshape(2)
    if poses.shape[1] != 3:
        raise ValueError("poses must have shape (T, 3)")
    T = poses.shape[0]
    n_rx = rx_offsets.shape[0]
    p = poses[:, :2]                       # (T, 2)
    theta = poses[:, 2]
    c, s = np.cos(theta), np.sin(theta)
    R = np.empty((T, 2, 2), dtype=float)
    R[:, 0, 0] = c
    R[:, 0, 1] = -s
    R[:, 1, 0] = s
    R[:, 1, 1] = c

    rx_world = p[:, None, :] + np.einsum("tij,aj->tai", R, rx_offsets)
    tx_world = p + np.einsum("tij,j->ti", R, tx_offset)

    drx = np.zeros((T, 3, n_rx, 2), dtype=float)
    dtx = np.zeros((T, 3, 2), dtype=float)
    drx[:, 0, :, :] = (1.0, 0.0)           # d/d p_x
    drx[:, 1, :, :] = (0.0, 1.0)           # d/d p_y
    dtx[:, 0, :] = (1.0, 0.0)
    dtx[:, 1, :] = (0.0, 1.0)
    drx[:, 2, :, :] = np.einsum("ij,taj->tai", _J90, rx_world - p[:, None, :])
    dtx[:, 2, :] = np.einsum("ij,tj->ti", _J90, tx_world - p)
    return {
        "rx_world": rx_world,
        "tx_world": tx_world,
        "drx_world": drx,
        "dtx_world": dtx,
    }


# --------------------------------------------------------------------------
# Domain operator
# --------------------------------------------------------------------------

@lru_cache(maxsize=16)
def _domain_core(N: int, k_b: float) -> tuple[np.ndarray, float, np.ndarray]:
    """(points, h, G_D) with the disk-regularised domain propagator."""
    points, h = make_grid(N)
    raw = green_matrix(points, points, k_b)          # point values, diag = 0
    G_D = (k_b**2) * (h**2) * raw                    # off-diagonal quadrature
    np.fill_diagonal(G_D, (k_b**2) * self_cell_green(k_b, h))
    return points, h, G_D


# --------------------------------------------------------------------------
# Core operators
# --------------------------------------------------------------------------

def _Dchi_times(M_diag: np.ndarray, A: np.ndarray) -> np.ndarray:
    """diag(chi) @ A using only the diagonal vector chi (row scaling)."""
    return M_diag[:, None] * A


def build_operators(
    chi: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    N: int,
    k_b: float,
) -> dict:
    """Evaluate the full contrast-source operators for a fixed (chi, poses).

    Returns a dict (see module docstring for the conventions):
      grid metadata, points, h, G_D, M,
      G_S_list, E_inc_list, J_list, E_tot_list, F_list (per pose),
      state residual per pose, sigma_min(M), ||M||_2,
      sigma_min(M)/||M||_2, max |E_inc|, max |E_tot|.
    """
    chi = np.asarray(chi, dtype=float).reshape(-1)
    poses = np.atleast_2d(np.asarray(poses, dtype=float))
    points, h, G_D = _domain_core(int(N), float(k_b))
    S = points.shape[0]
    if chi.shape[0] != S:
        raise ValueError(f"chi length {chi.shape[0]} != N^2 = {S}")
    T = poses.shape[0]
    geom = pose_geometry(poses, rx_offsets, tx_offset)

    M = np.eye(S, dtype=np.complex128) - _Dchi_times(chi, G_D)
    lu, piv = lu_factor(M, check_finite=False)

    G_S_list = []
    E_inc_list = []
    J_list = []
    E_tot_list = []
    F_list = []
    residual_list = []
    max_abs_E_inc = 0.0
    max_abs_E_tot = 0.0

    for t in range(T):
        rx_t = geom["rx_world"][t]
        tx_t = geom["tx_world"][t]
        G_S = (k_b**2) * (h**2) * green_matrix(points, rx_t, k_b)
        # E_inc[n] = q g(z_n, s_tx), q = 1.0
        E_inc = green_matrix(tx_t[None, :], points, k_b)[:, 0]
        rhs = chi * E_inc                          # D_chi E_inc
        J = lu_solve((lu, piv), rhs, check_finite=False)
        E_tot = E_inc + G_D @ J
        F = G_S @ J
        res = np.linalg.norm(M @ J - rhs, ord=2)
        denom = np.linalg.norm(rhs, ord=2)
        residual_list.append(float(res / denom) if denom > 0.0 else np.nan)
        G_S_list.append(G_S)
        E_inc_list.append(E_inc)
        J_list.append(J)
        E_tot_list.append(E_tot)
        F_list.append(F)
        max_abs_E_inc = max(max_abs_E_inc, float(np.max(np.abs(E_inc))))
        max_abs_E_tot = max(max_abs_E_tot, float(np.max(np.abs(E_tot))))

    sv = svdvals(M, check_finite=False)
    sigma_min = float(sv[-1])
    M_norm = float(sv[0])
    max_field = max_abs_E_tot
    return {
        "N": int(N),
        "k_b": float(k_b),
        "h": h,
        "points": points,
        "chi": chi,
        "poses": poses,
        "geom": geom,
        "G_D": G_D,
        "M": M,
        "G_S_list": G_S_list,
        "E_inc_list": E_inc_list,
        "J_list": J_list,
        "E_tot_list": E_tot_list,
        "F_list": F_list,
        "state_residual_list": residual_list,
        "max_state_residual": float(max(residual_list)),
        "sigma_min_M": sigma_min,
        "M_norm": M_norm,
        "M_min_sigma": sigma_min,
        "sigma_min_over_M_norm": float(sigma_min / M_norm),
        "max_abs_E_inc": max_abs_E_inc,
        "max_abs_E_tot": max_abs_E_tot,
        "max_field": max_field,
    }


def _pose_G_S_Einc(geom: dict, points: np.ndarray, h: float, k_b: float, t: int):
    """G_S_t and E_inc_t for pose t (q = 1.0)."""
    rx_t = geom["rx_world"][t]
    tx_t = geom["tx_world"][t]
    G_S = (k_b**2) * (h**2) * green_matrix(points, rx_t, k_b)
    E_inc = green_matrix(tx_t[None, :], points, k_b)[:, 0]
    return G_S, E_inc


def build_AB(
    chi: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    N: int,
    k_b: float,
):
    """Full-wave map and pose Jacobians at a fixed (chi, poses).

    Returns (A_complex, B_complex, forward_F, diagnostics):
      A_complex: (T*n_rx, N^2), rows grouped by pose: block t = A_t.
      B_complex: (T*n_rx, 3T), block t has columns 3t:3t+3 = B_t.
      forward_F: (T*n_rx,) complex, full-wave scattered-field data F(chi).
      diagnostics: list over poses with A_t, B_t, F_t, E_inc_t, E_tot_t,
                   norm(A_t), sigma_max/norm etc. for downstream families.

    One LU factorisation of M is used for the N^2 + 3T right-hand sides.
    """
    ops = build_operators(chi, poses, rx_offsets, tx_offset, N, k_b)
    points, h, G_D = ops["points"], ops["h"], ops["G_D"]
    geom = ops["geom"]
    M = ops["M"]
    chi_v = ops["chi"]
    lu, piv = lu_factor(M, check_finite=False)
    S = points.shape[0]
    poses = np.atleast_2d(np.asarray(poses, dtype=float))
    T = poses.shape[0]
    n_rx = np.atleast_2d(np.asarray(rx_offsets, dtype=float)).shape[0]
    rows = T * n_rx

    A = np.empty((rows, S), dtype=np.complex128)
    B = np.zeros((rows, 3 * T), dtype=np.complex128)
    F_full = np.empty((rows,), dtype=np.complex128)
    diagnostics = []

    for t in range(T):
        G_S_t = ops["G_S_list"][t]
        E_inc_t = ops["E_inc_list"][t]
        J_t = ops["J_list"][t]
        E_tot_t = ops["E_tot_list"][t]
        F_t = ops["F_list"][t]
        sl = slice(t * n_rx, (t + 1) * n_rx)

        X = lu_solve((lu, piv), np.diag(E_tot_t), check_finite=False)  # (S, S)
        A_t = G_S_t @ X

        rx_t = geom["rx_world"][t]
        tx_t = geom["tx_world"][t]
        drx_t = geom["drx_world"][t]      # (3, n_rx, 2)
        dtx_t = geom["dtx_world"][t]      # (3, 2)

        # grad_1 g(r_a, z_n): first argument = receiver r_a (row), second = z_n.
        grad1 = green_grad_first(points, rx_t, k_b)   # (n_rx, S, 2)
        # grad_2 g(z_n, s_tx) = grad_s g = -grad_z g   (see SIGN NOTE); the
        # explicit helper applies the corrected plus sign for source motion.
        grad2 = green_grad_source(points, tx_t, k_b)  # (S, 2)

        B_t = np.empty((n_rx, 3), dtype=np.complex128)
        for l in range(3):
            # D_xG_S_l[a, n] = k_b^2 h^2 grad_1 g(r_a, z_n) . (d r_a / d x_l)
            DXGS_l = (k_b**2) * (h**2) * np.einsum(
                "asd,ad->as", grad1, drx_t[l]
            )
            # D_xE_inc_l[n] = q grad_2 g(z_n, s_tx) . (d s_tx / d x_l), q = 1
            DXE_l = np.einsum("sd,d->s", grad2, dtx_t[l])
            Y = lu_solve((lu, piv), chi_v * DXE_l, check_finite=False)
            B_t[:, l] = DXGS_l @ J_t + G_S_t @ Y

        A[sl, :] = A_t
        B[sl, 3 * t : 3 * t + 3] = B_t
        F_full[sl] = F_t
        diagnostics.append(
            {
                "pose": poses[t],
                "A_t": A_t,
                "B_t": B_t,
                "F_t": F_t,
                "E_inc_t": E_inc_t,
                "E_tot_t": E_tot_t,
                "norm_A_t": float(np.linalg.norm(A_t, ord=2)),
                "sigma_min_A_t": float(
                    np.linalg.svd(A_t, compute_uv=False)[-1]
                    if min(A_t.shape) > 0
                    else np.nan
                ),
                "norm_B_t": float(np.linalg.norm(B_t, ord=2)),
            }
        )

    return A, B, F_full, diagnostics


def whiten_realify(A: np.ndarray, B: np.ndarray, W: np.ndarray | None = None):
    """Whiten with Sigma^{-1/2} = W (identity noise if W is None), then realify.

    A_hat = W A, B_hat = W B;
    A_R = sqrt(2) [Re(A_hat); Im(A_hat)], B_R = sqrt(2) [Re(B_hat); Im(B_hat)].
    Real parameters chi and X; real data dimension 2*(T*n_rx).
    """
    A = np.asarray(A)
    B = np.asarray(B)
    if W is not None:
        W = np.asarray(W)
        A_hat = W @ A
        B_hat = W @ B
    else:
        A_hat, B_hat = A, B
    A_R = _SQRT2 * np.vstack([np.real(A_hat), np.imag(A_hat)])
    B_R = _SQRT2 * np.vstack([np.real(B_hat), np.imag(B_hat)])
    return A_R, B_R


def born_forward(
    chi: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    N: int,
    k_b: float,
):
    """Born scattered field and Born map Jacobian.

    Returns (F_born, A_born), stacked over poses as (T*n_rx, N^2) for A_born
    and (T*n_rx,) for F_born, with
      A_born,t = G_S diag(E_inc),  F_born,t = G_S (chi * E_inc).
    """
    chi = np.asarray(chi, dtype=float).reshape(-1)
    poses = np.atleast_2d(np.asarray(poses, dtype=float))
    points, h, _ = _domain_core(int(N), float(k_b))
    geom = pose_geometry(poses, rx_offsets, tx_offset)
    T = poses.shape[0]
    n_rx = np.atleast_2d(np.asarray(rx_offsets, dtype=float)).shape[0]
    rows = T * n_rx
    S = points.shape[0]
    A_born = np.empty((rows, S), dtype=np.complex128)
    F_born = np.empty((rows,), dtype=np.complex128)
    for t in range(T):
        G_S_t, E_inc_t = _pose_G_S_Einc(geom, points, h, float(k_b), t)
        sl = slice(t * n_rx, (t + 1) * n_rx)
        A_born[sl, :] = G_S_t * E_inc_t[None, :]      # G_S diag(E_inc)
        F_born[sl] = G_S_t @ (chi * E_inc_t)
    return F_born, A_born


def forward_measurements(
    chi: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    N: int,
    k_b: float,
) -> np.ndarray:
    """Full-wave scattered-field data vector F(chi, X), shape (T*n_rx,)."""
    ops = build_operators(chi, poses, rx_offsets, tx_offset, N, k_b)
    return np.concatenate([F for F in ops["F_list"]])
