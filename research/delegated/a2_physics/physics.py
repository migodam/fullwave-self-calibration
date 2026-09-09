"""Reusable 2D scalar full-wave physical core for the A2 PRASC-SOM package.

This module implements the *physical forward model and analytic parameter
derivatives* only.  It makes no research claims and does not run any
experiment at import time.

Conventions (fixed by ``research/delegated/a2_physics/TASK.md`` and the
theorem package sections 2, 8 and 9):

Time / Green function
  Time convention ``e^{-i omega t}``.  The 2D free-space outgoing Green
  function is

      g(r, r') = (i/4) H_0^{(1)}(k |r - r'|),      e^{-i omega t}

  whose large-argument form is a growing phase ``exp(+i k r)``, i.e. an
  outgoing wave under ``e^{-i omega t}``.  ``point_green`` exposes the point
  value so the convention is directly testable.

Domain / grid
  World-fixed ``D = [-0.5, 0.5]^2`` metres, uniform ``N x N`` cell-centre
  grid, cell side ``h = 1/N``.  Points are flattened row-major with index
  ``n = i*N + j``, where ``(x_i, y_j)`` are the cell centres.  The inverse
  model uses ``N=16``, truth evaluation ``N=32``, and ``N=8`` is accepted
  only for smoke tests (enforced by ``Config`` validation).

Operators
  Both domain (D) and receiver (S) scattering operators carry the
  ``k^2`` frequency factor and the cell-area integration:

      D[a, b] = k^2 h^2 g(z_a, z_b)          (a != b)
      D[a, a] = k^2 I_self                    (equal-area disk integral)
      S[r, n] = k^2 h^2 g(rx_r, z_n)

  with the corrected equal-area disk self-cell value

      I_self = (i pi a/(2 k)) H_1^{(1)}(k a) - 1/k^2,   a = h/sqrt(pi).

  The ``-1/k^2`` lower-endpoint term is essential; omitting it makes the
  cell integral tend to ``1/k^2`` rather than zero as ``h -> 0``.

Material / dispersion
  Nine unnormalised Gaussian basis functions on the Cartesian square of
  ``{-0.25, 0, 0.25}`` metres:

      Phi[n, j] = exp(-|z_n - c_j|^2 / (2 width^2)),   width = 0.16 m.

  Real coefficients ``alpha`` (length 9).  ``epsilon_r - 1 = Phi @ alpha``
  and the declared Ohmic law ``sigma = 0.005*(epsilon_r - 1)`` S/m gives
  the frequency-dependent complex contrast

      chi_omega = (Phi @ alpha) * (1 + i*0.005/(omega*epsilon0)).

  ``c0 = 299792458``, ``epsilon0 = 8.8541878128e-12``, ``omega = k*c0``.

State equations
  ``M j = b`` with ``M = I - diag(chi) D`` and ``b = chi * E_inc``, where
  ``E_inc[n] = g(z_n, tx)`` is the unit point-source incident field.  The
  internal total field is ``E_tot = E_inc + D j`` and the scattered data are
  ``S j``.  Raw data are the TOTAL field at receivers:

      total = direct(rx, tx) + S j,    direct(rx, tx) = g(rx, tx).

  ``incident`` and ``scattered`` are retained separately in the return dict.

Geometry / pose
  Three nominal poses (px, py, theta):
      (-0.15, -0.12, 0), (0.15, 0, 0.65), (0, 0.15, 1.3).
  Pose 0 is anchored.  One shared unknown ``x = (dx, dy, dtheta)`` adds the
  same world translation and rotation angle to poses 1 and 2:

      pose_eff[t] = pose_nominal[t] + x      for t in {1, 2}, else nominal.

  This is NOT six independent pose unknowns.  Receivers and transmitters
  co-move rigidly with each pose, so all intra-pose receiver/transmitter
  distances (and therefore the direct term) are invariant under ``x``.  The
  pose metric is ``diag(1, 1, R_eff^2)`` with ``R_eff = 1.5`` m.

  Receiver body frame: radius 1.5 m, 12 angles
      full:    2*pi*m/12
      limited: linspace(-pi/4, pi/4, 12)
  Transmitter body frame: radius 1.9 m, two angles {-pi/3, +pi/3}.

Derivatives
  For the reduced physical model the exact tangents are

      A = S M^{-1} diag(E_tot) T,        T = d chi / d alpha,
      B = H_S + S M^{-1} H_D + C_x,

  with receiver term ``H_S = (dS/dx) j``, material term
  ``H_D = chi * (dE_inc/dx)`` and the direct-incident term ``C_x`` (exactly
  zero here because receivers and transmitters co-move rigidly, but kept in
  the code so a general geometry cannot silently drop it).  All
  derivatives are analytic; the domain matrix ``D`` is pose-independent so
  ``dD/dx = 0``.

Work accounting
  One LU factorisation is shared by all illuminations/poses at a given
  frequency.  Every individual right-hand side solved through an LU is
  charged, not each batch.  Factorisations are counted and timed
  separately.  ``adjoint_gradient`` charges every forward state solve and
  every adjoint right-hand side it performs.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1

__all__ = [
    "C0",
    "EPS0",
    "SELF_CELL_FORMULA",
    "Config",
    "Model",
    "make_grid",
    "point_green",
    "self_cell_green",
    "green_matrix",
    "green_grad_first",
    "green_grad_source",
    "pose_geometry",
    "gaussian_basis",
    "realify",
    "realify_jacobian",
]

C0 = 299792458.0
EPS0 = 8.8541878128e-12
SQRT2 = float(np.sqrt(2.0))

SELF_CELL_FORMULA = (
    "I_self = (i*pi*a/(2*k))*H_1^(1)(k*a) - 1/k^2, a = h/sqrt(pi)"
)

_J90 = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=float)
_ALLOWED_APERTURES = ("full", "limited")
_ALLOWED_N = (8, 16, 32)


# ---------------------------------------------------------------------------
# Grid, Green functions and basis
# ---------------------------------------------------------------------------

def make_grid(N: int) -> tuple[np.ndarray, float]:
    """Cell-centre grid for ``[-0.5, 0.5]^2`` with side ``h = 1/N``."""
    if not isinstance(N, (int, np.integer)) or int(N) < 1:
        raise ValueError(f"N must be a positive integer, got {N!r}")
    N = int(N)
    h = 1.0 / N
    centres = -0.5 + (np.arange(N, dtype=float) + 0.5) * h
    X, Y = np.meshgrid(centres, centres, indexing="ij")
    points = np.empty((N * N, 2), dtype=float)
    points[:, 0] = X.ravel(order="C")
    points[:, 1] = Y.ravel(order="C")
    return points, h


def point_green(R: float | np.ndarray, k: float) -> np.ndarray:
    """``g(R) = (i/4) H_0^{(1)}(k R)`` for positive distance(s) ``R``."""
    R = np.asarray(R, dtype=float)
    k = float(k)
    out = np.zeros(R.shape, dtype=np.complex128)
    safe = R > 0.0
    if safe.any():
        out[safe] = (1j / 4.0) * hankel1(0, k * R[safe])
    return out


def self_cell_green(k: float, h: float) -> complex:
    """Complete equal-area disk self-cell integral of ``g`` (area units).

    ``a = h/sqrt(pi)`` is the radius of a disk of area ``h^2``.  The
    ``-1/k^2`` endpoint term is included (see module docstring).
    """
    k = float(k)
    a = float(h) / np.sqrt(np.pi)
    return (1j * np.pi * a / (2.0 * k)) * hankel1(1, k * a) - 1.0 / k**2


def _pair_diff_norm(
    src: np.ndarray, dst: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    diff = dst[:, None, :] - src[None, :, :]
    R = np.linalg.norm(diff, axis=2)
    return diff, R


def green_matrix(
    src: np.ndarray, dst: np.ndarray, k: float
) -> np.ndarray:
    """Point Green matrix ``G[d, s] = g(dst[d], src[s])``.

    Coincident pairs return 0; the singular point value is not used.
    Callers apply the self-cell regularisation and ``k^2`` factors.
    """
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    if src.shape[1] != 2 or dst.shape[1] != 2:
        raise ValueError("src and dst must have shape (..., 2)")
    _, R = _pair_diff_norm(src, dst)
    G = np.zeros(R.shape, dtype=np.complex128)
    safe = R > 0.0
    if safe.any():
        G[safe] = (1j / 4.0) * hankel1(0, float(k) * R[safe])
    return G


def green_grad_first(src: np.ndarray, dst: np.ndarray, k: float) -> np.ndarray:
    """Gradient of ``g`` with respect to its first (field) argument.

    Returns ``V[d, s, :] = grad_r g(r, r')`` at ``r = dst[d]``,
    ``r' = src[s]``; zero vector for coincident pairs.
    """
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    if src.shape[1] != 2 or dst.shape[1] != 2:
        raise ValueError("src and dst must have shape (..., 2)")
    diff, R = _pair_diff_norm(src, dst)
    V = np.zeros(diff.shape, dtype=np.complex128)
    safe = R > 0.0
    if safe.any():
        H = hankel1(1, float(k) * R[safe])
        V[safe] = (
            (-(1j * float(k)) / 4.0)
            * H[:, None]
            * (diff[safe] / R[safe][:, None])
        )
    return V


def green_grad_source(z: np.ndarray, s: np.ndarray, k: float) -> np.ndarray:
    """``grad_s g(z_n, s) = -grad_z g(z_n, s)`` for one source point ``s``.

    Returns ``(n_z, 2)``; zero vector wherever ``z_n == s``.
    """
    z = np.atleast_2d(np.asarray(z, dtype=float))
    s = np.atleast_2d(np.asarray(s, dtype=float))
    if z.shape[1] != 2 or s.shape[1] != 2:
        raise ValueError("z and s must have shape (..., 2)")
    return -green_grad_first(s, z, k)[:, 0, :]


def gaussian_basis(
    points: np.ndarray, width: float, centres: tuple[float, float, float]
) -> np.ndarray:
    """Unnormalised Gaussian columns ``exp(-|x-c|^2/(2 width^2))``.

    The nine centres are the Cartesian product ``centres x centres``.
    """
    width = float(width)
    if width <= 0.0:
        raise ValueError("basis width must be positive")
    cols = []
    for cx in centres:
        for cy in centres:
            c = np.array([float(cx), float(cy)], dtype=float)
            d = points - c[None, :]
            cols.append(np.exp(-np.einsum("ij,ij->i", d, d) / (2.0 * width**2)))
    return np.stack(cols, axis=1)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def pose_geometry(
    poses: np.ndarray,
    rx_body: np.ndarray,
    tx_body: np.ndarray,
) -> dict[str, np.ndarray]:
    """World positions and per-pose-parameter derivatives.

    ``poses``: ``(T, 3)`` ``(px, py, theta)``.
    ``rx_body``: ``(R, 2)``, ``tx_body``: ``(U, 2)`` body offsets.
    Returns world receiver/transmitter positions and derivatives with
    respect to the three parameters of each pose:
      ``drx_pose`` ``(T, 3, R, 2)``, ``dtx_pose`` ``(T, 3, U, 2)``.
    """
    poses = np.atleast_2d(np.asarray(poses, dtype=float))
    rx_body = np.atleast_2d(np.asarray(rx_body, dtype=float))
    tx_body = np.atleast_2d(np.asarray(tx_body, dtype=float))
    if poses.shape[1] != 3:
        raise ValueError("poses must have shape (T, 3)")
    T = poses.shape[0]
    R = rx_body.shape[0]
    U = tx_body.shape[0]
    p = poses[:, :2]
    theta = poses[:, 2]
    c, s = np.cos(theta), np.sin(theta)
    rot = np.empty((T, 2, 2), dtype=float)
    rot[:, 0, 0] = c
    rot[:, 0, 1] = -s
    rot[:, 1, 0] = s
    rot[:, 1, 1] = c

    rx_world = p[:, None, :] + np.einsum("tij,aj->tai", rot, rx_body)
    tx_world = p[:, None, :] + np.einsum("tij,uj->tui", rot, tx_body)

    drx = np.zeros((T, 3, R, 2), dtype=float)
    dtx = np.zeros((T, 3, U, 2), dtype=float)
    drx[:, 0, :, 0] = 1.0
    drx[:, 1, :, 1] = 1.0
    dtx[:, 0, :, 0] = 1.0
    dtx[:, 1, :, 1] = 1.0
    drx[:, 2] = np.einsum("ij,taj->tai", _J90, rx_world - p[:, None, :])
    dtx[:, 2] = np.einsum("ij,tuj->tui", _J90, tx_world - p[:, None, :])
    return {
        "rx_world": rx_world,
        "tx_world": tx_world,
        "drx_pose": drx,
        "dtx_pose": dtx,
    }


# ---------------------------------------------------------------------------
# Config / model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Config:
    """Fixed physical/experimental configuration for the A2 core.

    ``material_basis`` optionally replaces the nine Gaussian columns with a
    caller-supplied ``(N*N, q)`` spatial basis (e.g. pixel identity for
    later imaging).  The same real-coefficient Ohmic dispersion law is
    applied, so ``q`` may differ from 9.
    """

    N: int = 16
    aperture: str = "full"
    kmax: float = 4 * np.pi
    material_basis: np.ndarray | None = None
    basis_centres: tuple[float, float, float] = (-0.25, 0.0, 0.25)
    basis_width: float = 0.16
    loss_coefficient: float = 0.005
    receiver_radius: float = 1.5
    transmitter_radius: float = 1.9
    receiver_angles_full: tuple[float, ...] | None = None
    transmitter_body_angles: tuple[float, float] = (-np.pi / 3.0, np.pi / 3.0)
    pose_metric_rotation_weight: float = 1.5

    def __post_init__(self) -> None:
        if int(self.N) not in _ALLOWED_N:
            raise ValueError(
                "N must be 8 (smoke only), 16 (inverse) or 32 (truth); "
                f"got {self.N!r}"
            )
        if self.aperture not in _ALLOWED_APERTURES:
            raise ValueError(
                f"aperture must be one of {_ALLOWED_APERTURES}; got "
                f"{self.aperture!r}"
            )
        if not np.isfinite(self.kmax) or self.kmax <= 0.0:
            raise ValueError(f"kmax must be a positive finite number")
        if self.material_basis is not None:
            mb = np.asarray(self.material_basis)
            if mb.ndim != 2 or mb.shape[0] != int(self.N) ** 2:
                raise ValueError(
                    "material_basis must have shape (N*N, q) with rows "
                    f"N*N={int(self.N) ** 2}"
                )


class Model:
    """Full-wave forward model for the fixed A2 physical specification."""

    def __init__(self, config: Config | None = None) -> None:
        self.config = config if config is not None else Config()
        cfg = self.config
        self.points, self.h = make_grid(int(cfg.N))
        self.n_cells = self.points.shape[0]

        if cfg.material_basis is None:
            self.basis = gaussian_basis(
                self.points, cfg.basis_width, cfg.basis_centres
            )
        else:
            mb = np.asarray(cfg.material_basis, dtype=float)
            if mb.shape[0] != self.n_cells:
                raise ValueError(
                    "material_basis rows do not match N*N for this model"
                )
            self.basis = mb.copy()
        self.n_alpha = self.basis.shape[1]

        self.ks = cfg.kmax * np.array([0.25, 0.5, 0.75, 1.0])
        self.n_frequencies = len(self.ks)

        if cfg.receiver_angles_full is None:
            rx_angles = (
                2.0 * np.pi * np.arange(12) / 12.0
                if cfg.aperture == "full"
                else np.linspace(-np.pi / 4.0, np.pi / 4.0, 12)
            )
        else:
            rx_angles = np.asarray(cfg.receiver_angles_full, dtype=float)
            if cfg.aperture == "limited":
                raise ValueError(
                    "receiver_angles_full is only meaningful for aperture='full'"
                )
        if rx_angles.ndim != 1 or rx_angles.shape[0] < 1:
            raise ValueError("receiver angles must be a 1-D array")
        self.rx_body = cfg.receiver_radius * np.stack(
            (np.cos(rx_angles), np.sin(rx_angles)), axis=1
        )
        tx_angles = np.asarray(cfg.transmitter_body_angles, dtype=float)
        self.tx_body = cfg.transmitter_radius * np.stack(
            (np.cos(tx_angles), np.sin(tx_angles)), axis=1
        )

        self.nominal_poses = np.array(
            [[-0.15, -0.12, 0.0], [0.15, 0.0, 0.65], [0.0, 0.15, 1.3]],
            dtype=float,
        )
        self.n_pose = 3
        self.n_tx = self.tx_body.shape[0]
        self.n_rx = self.rx_body.shape[0]
        w = cfg.pose_metric_rotation_weight
        self.pose_metric = np.diag([1.0, 1.0, float(w) ** 2])

        self._domain_cache: dict[int, np.ndarray] = {}
        self._work: dict[str, float | int] = {}

    # -- geometry -----------------------------------------------------------

    def _effective_geometry(self, x: np.ndarray) -> dict[str, np.ndarray]:
        """Effective poses (anchor plus shared ``x``) and x-derivatives."""
        eff = self.nominal_poses.copy()
        eff[1:, :] = eff[1:, :] + np.asarray(x, dtype=float)[None, :]
        base = pose_geometry(eff, self.rx_body, self.tx_body)
        # derivative of effective pose t wrt shared x is identity for t>=1,
        # zero for the anchored pose 0.
        drx_x = np.zeros_like(base["drx_pose"])
        dtx_x = np.zeros_like(base["dtx_pose"])
        drx_x[1:] = base["drx_pose"][1:]
        dtx_x[1:] = base["dtx_pose"][1:]
        return {
            "poses_eff": eff,
            "rx_world": base["rx_world"],
            "tx_world": base["tx_world"],
            "drx_x": drx_x,
            "dtx_x": dtx_x,
        }

    # -- cached domain operator ---------------------------------------------

    def _domain_operator(self, freq_idx: int) -> np.ndarray:
        if freq_idx in self._domain_cache:
            return self._domain_cache[freq_idx]
        k = float(self.ks[freq_idx])
        raw = green_matrix(self.points, self.points, k)
        D = (k**2) * (self.h**2) * raw
        np.fill_diagonal(D, (k**2) * self_cell_green(k, self.h))
        self._domain_cache[freq_idx] = D
        return D

    # -- work bookkeeping ---------------------------------------------------

    def _reset_work(self) -> None:
        self._work = {
            "factorizations": 0,
            "factorization_seconds": 0.0,
            "rhs_solves_total": 0,
            "rhs_solves": {
                "forward_state": 0,
                "material": 0,
                "pose": 0,
                "adjoint": 0,
            },
            "rhs_solve_seconds": 0.0,
            "operator_products": 0,
        }

    def _tick_product(self, count: int = 1) -> None:
        self._work["operator_products"] += int(count)

    def _factorize(self, M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        t0 = time.perf_counter()
        lu, piv = lu_factor(M, check_finite=False)
        self._work["factorizations"] += 1
        self._work["factorization_seconds"] += time.perf_counter() - t0
        return lu, piv

    def _solve(
        self,
        lu_piv: tuple[np.ndarray, np.ndarray],
        rhs: np.ndarray,
        kind: str,
        trans: int = 0,
    ) -> np.ndarray:
        n_col = 1 if rhs.ndim == 1 else rhs.shape[1]
        t0 = time.perf_counter()
        out = lu_solve(lu_piv, rhs, trans=trans, check_finite=False)
        self._work["rhs_solve_seconds"] += time.perf_counter() - t0
        self._work["rhs_solves_total"] += n_col
        self._work["rhs_solves"][kind] += n_col
        return out

    # -- validation helpers -------------------------------------------------

    def _validate_alpha(self, alpha: np.ndarray) -> np.ndarray:
        a = np.asarray(alpha, dtype=float).reshape(-1)
        if a.shape[0] != self.n_alpha:
            raise ValueError(
                f"alpha length {a.shape[0]} != n_alpha {self.n_alpha}"
            )
        if not np.all(np.isfinite(a)):
            raise ValueError("alpha must be finite")
        return a

    def _validate_x(self, x: np.ndarray) -> np.ndarray:
        v = np.asarray(x, dtype=float).reshape(3)
        if not np.all(np.isfinite(v)):
            raise ValueError("x must be a finite length-3 vector")
        return v

    def _normalize_freq_ids(
        self, freq_ids: int | list[int] | np.ndarray | None
    ) -> list[int]:
        if freq_ids is None:
            return list(range(self.n_frequencies))
        arr = np.asarray(freq_ids).reshape(-1)
        ids = [int(i) for i in arr]
        if any(i < 0 or i >= self.n_frequencies for i in ids):
            raise ValueError(
                f"freq_ids must be in 0..{self.n_frequencies - 1}"
            )
        return sorted(set(ids))

    # -- core computation ----------------------------------------------------

    def _compute(
        self,
        alpha: np.ndarray,
        x: np.ndarray,
        freq_ids: list[int],
        jacobian: bool,
    ) -> dict:
        cfg = self.config
        S = self.n_cells
        R = self.n_rx
        U = self.n_tx
        T = self.n_pose
        geom = self._effective_geometry(x)
        n_states = len(freq_ids) * T * U
        n_rows = n_states * R

        total = np.empty(n_rows, dtype=np.complex128)
        scattered = np.empty(n_rows, dtype=np.complex128)
        incident = np.empty(n_rows, dtype=np.complex128)
        states: list[dict] = []
        blocks: dict[str, list] = {
            "freq_idx": [],
            "pose_idx": [],
            "illum_idx": [],
            "row_start": [],
            "row_stop": [],
        }

        lu_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        chi_cache: dict[int, np.ndarray] = {}
        tmat_cache: dict[int, np.ndarray] = {}
        m_cache: dict[int, np.ndarray] = {}
        row = 0

        for fi in freq_ids:
            k = float(self.ks[fi])
            omega = k * C0
            loss_factor = 1.0 + 1j * float(cfg.loss_coefficient) / (omega * EPS0)
            D = self._domain_operator(fi)
            chi = (self.basis @ alpha) * loss_factor
            Tmat = self.basis * loss_factor
            chi_cache[fi] = chi
            tmat_cache[fi] = Tmat

            M = np.eye(S, dtype=np.complex128)
            M -= chi[:, None] * D  # diag(chi) @ D (row scaling)
            m_cache[fi] = M
            lu_cache[fi] = self._factorize(M)
            lu_piv = lu_cache[fi]

            for p in range(T):
                rx_world = geom["rx_world"][p]
                tx_world = geom["tx_world"][p]
                drx_x = geom["drx_x"][p]  # (3, R, 2)
                dtx_x = geom["dtx_x"][p]  # (3, U, 2)

                S_mat = (k**2) * (self.h**2) * green_matrix(
                    self.points, rx_world, k
                )
                # receiver-side Green gradient for S derivative
                grad_rx = green_grad_first(self.points, rx_world, k)  # (R,S,2)
                dS_dx = []
                for l in range(3):
                    dS_dx.append(
                        (k**2)
                        * (self.h**2)
                        * np.einsum("asd,ad->as", grad_rx, drx_x[l])
                    )

                for u in range(U):
                    tx_u = tx_world[u]
                    E_inc = green_matrix(tx_u[None, :], self.points, k)[:, 0]
                    b = chi * E_inc
                    j = self._solve(lu_piv, b, "forward_state")
                    Etot = E_inc + D @ j
                    sca = S_mat @ j
                    direct = green_matrix(tx_u[None, :], rx_world, k)[:, 0]
                    tot = direct + sca
                    self._tick_product(3)  # M@j residual, D@j, S_mat@j
                    res = M @ j - b
                    res_norm = float(np.linalg.norm(res))
                    b_norm = float(np.linalg.norm(b))
                    res_rel = (
                        res_norm / b_norm if b_norm > 0.0 else float(np.nan)
                    )

                    # db/dx = chi * dE_inc/dx  (transmitter motion only)
                    grad_src = green_grad_source(self.points, tx_u, k)  # (S,2)
                    db_dx = chi[:, None] * np.einsum(
                        "nd,ld->nl", grad_src, dtx_x[:, u, :]
                    )
                    # direct-incident derivative via relative rx-tx motion
                    grad_at_rx = green_grad_first(
                        tx_u[None, :], rx_world, k
                    )[:, 0, :]  # (R, 2)
                    rel_dx = drx_x - dtx_x[:, u, None, :]  # (3,R,2)
                    direct_dx = np.einsum("ad,lad->la", grad_at_rx, rel_dx)  # (3,R)

                    start = row
                    stop = row + R
                    row = stop
                    total[start:stop] = tot
                    scattered[start:stop] = sca
                    incident[start:stop] = direct

                    states.append(
                        {
                            "freq_idx": fi,
                            "pose_idx": p,
                            "illum_idx": u,
                            "k": k,
                            "omega": omega,
                            "wavelength": 2.0 * np.pi / k,
                            "pose": geom["poses_eff"][p].copy(),
                            "rx_body": self.rx_body.copy(),
                            "tx_body": self.tx_body.copy(),
                            "rx_world": rx_world,
                            "tx_world": tx_world,
                            "M": M,
                            "D": D,
                            "S": S_mat,
                            "chi": chi,
                            "Tmat": Tmat,
                            "j": j,
                            "b": b,
                            "E_inc": E_inc,
                            "E_total": Etot,
                            "direct": direct,
                            "scattered": sca,
                            "total": tot,
                            "dS_dx": dS_dx,
                            "db_dx": db_dx,
                            "direct_derivative_dx": direct_dx,  # (3, R)
                            "state_residual_abs": res_norm,
                            "state_residual_rel": res_rel,
                            "row_start": start,
                            "row_stop": stop,
                        }
                    )
                    blocks["freq_idx"].append(fi)
                    blocks["pose_idx"].append(p)
                    blocks["illum_idx"].append(u)
                    blocks["row_start"].append(start)
                    blocks["row_stop"].append(stop)

        A = None
        B = None
        if jacobian:
            A = np.empty((n_rows, self.n_alpha), dtype=np.complex128)
            B = np.zeros((n_rows, 3), dtype=np.complex128)
            for st, state in enumerate(states):
                fi = state["freq_idx"]
                lu_piv = lu_cache[fi]
                Tmat = tmat_cache[fi]
                j = state["j"]
                S_mat = state["S"]
                Etot = state["E_total"]
                db_dx = state["db_dx"]

                # material tangent: X = M^{-1} diag(E_tot) Tmat  (q RHS)
                X = self._solve(
                    lu_piv, Etot[:, None] * Tmat, "material"
                )
                A_block = S_mat @ X
                # pose tangent through b:  Y = M^{-1} db/dx  (3 RHS)
                Y = self._solve(lu_piv, db_dx, "pose")
                B_block = S_mat @ Y
                self._tick_product(5)  # S@X, S@Y, and 3 dS_dx @ j
                for l in range(3):
                    B_block[:, l] += state["dS_dx"][l] @ j + state[
                        "direct_derivative_dx"
                    ][l, :]
                st_state = states[st]
                st_state["A_block"] = A_block
                st_state["B_block"] = B_block
                sl = slice(st_state["row_start"], st_state["row_stop"])
                A[sl, :] = A_block
                B[sl, :] = B_block

        result = {
            "total": total,
            "scattered": scattered,
            "incident": incident,
            "states": states,
            "blocks": blocks,
            "lu_cache": lu_cache,
            "geometry": geom,
        }
        if jacobian:
            result["A"] = A
            result["B"] = B
        return result

    # -- public API ----------------------------------------------------------

    def forward(
        self,
        alpha: np.ndarray,
        x: np.ndarray,
        freq_ids: int | list[int] | np.ndarray | None = None,
        jacobian: bool = False,
    ) -> dict:
        """Run the full-wave forward model and optionally exact tangents.

        Stable order: frequency, pose, illumination, receiver.
        Returns dict with ``total``, ``scattered``, ``incident``, ``A`` and
        ``B`` (complex, only when ``jacobian=True``), ``blocks``, per-state
        ``states`` and ``work``.
        """
        alpha = self._validate_alpha(alpha)
        x = self._validate_x(x)
        ids = self._normalize_freq_ids(freq_ids)
        self._reset_work()
        t0 = time.perf_counter()
        core = self._compute(alpha, x, ids, bool(jacobian))
        wall = time.perf_counter() - t0
        work = self._work.copy()
        work["wall_seconds"] = wall
        work["frequencies"] = len(ids)
        work["states"] = len(core["states"])
        out = {
            "total": core["total"],
            "scattered": core["scattered"],
            "incident": core["incident"],
            "states": core["states"],
            "blocks": core["blocks"],
            "work": work,
        }
        if jacobian:
            out["A"] = core["A"]
            out["B"] = core["B"]
        return out

    def adjoint_gradient(
        self,
        alpha: np.ndarray,
        x: np.ndarray,
        y: np.ndarray,
        sigma: float = 1.0,
        freq_ids: int | list[int] | np.ndarray | None = None,
    ) -> dict:
        """Gradient of ``0.5 ||realify(total - y, sigma)||^2`` by adjoints.

        Charges every forward-state solve and every adjoint right-hand side
        in ``work``.  Uses the existing per-frequency LU factorisation
        (``trans=2``), so no additional factorisation is performed.
        """
        alpha = self._validate_alpha(alpha)
        x = self._validate_x(x)
        sigma = float(sigma)
        if not np.isfinite(sigma) or sigma <= 0.0:
            raise ValueError("sigma must be a positive finite number")
        ids = self._normalize_freq_ids(freq_ids)
        self._reset_work()
        t0 = time.perf_counter()
        core = self._compute(alpha, x, ids, jacobian=False)
        y_arr = np.asarray(y, dtype=np.complex128).reshape(-1)
        if y_arr.shape[0] != core["total"].shape[0]:
            raise ValueError(
                "y length must match the selected forward rows "
                f"({core['total'].shape[0]}), got {y_arr.shape[0]}"
            )
        res = core["total"] - y_arr
        grad_alpha = np.zeros(self.n_alpha, dtype=np.complex128)
        grad_x = np.zeros(3, dtype=np.complex128)
        scale = 2.0 / sigma**2

        for st, state in enumerate(core["states"]):
            fi = state["freq_idx"]
            lu_piv = core["lu_cache"][fi]
            w = res[state["row_start"] : state["row_stop"]]
            S_mat = state["S"]
            v = S_mat.conj().T @ w
            lam = self._solve(lu_piv, v, "adjoint", trans=2)
            self._tick_product(1)  # S^H w

            Tmat = state["Tmat"]
            Etot = state["E_total"]
            # A^H w = Tmat^H (conj(E_tot) .* lam)
            grad_alpha += Tmat.conj().T @ (np.conj(Etot) * lam)
            self._tick_product(1)

            # B^H w terms
            C = state["direct_derivative_dx"]  # (3, R)
            b_x = np.einsum("la,a->l", C.conj(), w)  # direct term (3,)
            self._tick_product(1)
            for l in range(3):
                u = state["dS_dx"][l].conj().T @ w  # receiver sampling
                b_x[l] += np.conj(state["j"]) @ u
            self._tick_product(3)
            b_x += np.conj(state["db_dx"]).T @ lam  # material RHS in b
            self._tick_product(1)
            grad_x += b_x

        loss = 0.5 * float(np.sum(np.abs(realify(res, sigma)) ** 2))
        wall = time.perf_counter() - t0
        work = self._work.copy()
        work["wall_seconds"] = wall
        work["frequencies"] = len(ids)
        work["states"] = len(core["states"])
        return {
            "grad_alpha": scale * np.real(grad_alpha),
            "grad_x": scale * np.real(grad_x),
            "loss": loss,
            "work": work,
        }


# ---------------------------------------------------------------------------
# Realification
# ---------------------------------------------------------------------------

def realify(z: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    """``sqrt(2)/sigma * [real(z); imag(z)]`` for complex data ``z``.

    With proper complex noise of per-real-component variance ``sigma^2/2``,
    the realified vector has identity covariance.
    """
    z = np.asarray(z)
    sigma = float(sigma)
    if sigma <= 0.0 or not np.isfinite(sigma):
        raise ValueError("sigma must be a positive finite number")
    return (SQRT2 / sigma) * np.concatenate([np.real(z), np.imag(z)])


def realify_jacobian(J: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    """Realified Jacobian ``sqrt(2)/sigma * [Re(J); Im(J)]``."""
    J = np.asarray(J)
    sigma = float(sigma)
    if sigma <= 0.0 or not np.isfinite(sigma):
        raise ValueError("sigma must be a positive finite number")
    return (SQRT2 / sigma) * np.vstack([np.real(J), np.imag(J)])
