"""Matrix-free FFT dyadic convolution and GMRES for the parent 3D voxel VIE.

This is known-algorithm engineering for the *same* occupied regular-grid dipole
VIE implemented in ``research/trispace_self_calibration/a3_research/maxwell3d.py``.
The parent file is read-only reference material; no parent source is modified.

Conventions are inherited exactly from the parent module (which is imported
read-only and used as the oracle in tests):

* time convention exp(-i omega t);
* point-major, component-minor vector ordering ``(p*3 + c)``;
* voxel polarizability
  ``alpha = a0/(1 - 1j*k**3*a0/(6*pi))`` with
  ``a0 = 3*h**3*fill*(eps-1)/(eps+2)``;
* subcell fill from the same ``fill_quadrature`` grid;
* physical currents solve ``(I - diag(alpha) K) p = alpha*E_inc``, i.e. the
  exact physical state with no additional free-current nuisance freedom
  (r_free = 0).

The solver embeds the occupied cells in the parent bounding rectangular grid
and performs the dyadic convolution with a zero-padded FFT.  The linear
convolution padding factor is 2 along every axis.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

# Read-only import of the parent 3D Maxwell reference/voxelisation code.
_PARENT_DIR = (
    Path(__file__).resolve().parents[3]
    / "research"
    / "trispace_self_calibration"
    / "a3_research"
)
sys.path.insert(0, str(_PARENT_DIR))

import maxwell3d as _parent  # noqa: E402

DipoleVIE = _parent.DipoleVIE
dipole_kernel = _parent.dipole_kernel
voxelize = _parent.voxelize
illuminations = _parent.illuminations
receivers = _parent.receivers
treams_field = _parent.treams_field


def bounding_grid(
    centers: Sequence[Sequence[float]],
    radii: Sequence[float],
    spacing: float,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int], np.ndarray]:
    """Return the parent rectangular embedding grid for the occupied scene.

    ``lo``/``hi`` are the voxel integer bounds used by ``voxelize``.  Returns
    ``(lo, hi, rect_shape, full_points)`` where ``full_points`` are all cell
    centres on the rectangular grid (occupied and empty), in the same
    row-major order as ``np.meshgrid(..., indexing='ij')``.
    """
    centers_arr = np.asarray(centers, dtype=float)
    radii_arr = np.asarray(radii, dtype=float)
    lo = np.floor(
        np.min(centers_arr - radii_arr[:, None], axis=0) / spacing
    ).astype(np.intp)
    hi = np.ceil(
        np.max(centers_arr + radii_arr[:, None], axis=0) / spacing
    ).astype(np.intp)
    shape = tuple(int(x) for x in (hi - lo))
    axes = [np.arange(lo[i], hi[i]) + 0.5 for i in range(3)]
    xyz = np.stack(
        np.meshgrid(*[a * spacing for a in axes], indexing="ij"), axis=-1
    ).reshape(-1, 3)
    return lo, hi, shape, xyz


def dyadic_kernel_values(
    rect_shape: tuple[int, int, int],
    pad_shape: tuple[int, int, int],
    spacing: float,
    k: float,
) -> np.ndarray:
    """Dyadic ``k**2*G`` on the padded lag grid, shape ``pad_shape + (3,3)``.

    For axis length ``L`` and padded length ``2*L``, padded index ``q`` maps
    to the physical integer lag ``q`` for ``q < L`` and ``q - 2*L`` otherwise.
    Values at the never-used ``q = L`` (integer lag ``-L``) are zero-filled.
    The zero-displacement self-cell entry is exactly zero, as in the parent
    ``dipole_kernel``.
    """
    Lx, Ly, Lz = rect_shape
    Nx, Ny, Nz = pad_shape
    L = np.asarray(rect_shape)
    N = np.asarray(pad_shape)
    q = np.stack(
        np.meshgrid(
            np.arange(Nx), np.arange(Ny), np.arange(Nz), indexing="ij"
        ),
        axis=-1,
    )
    lag = np.where(
        q < L[None, None, None, :],
        q,
        q - N[None, None, None, :],
    )
    dr = lag.astype(np.float64) * spacing
    d = np.linalg.norm(dr, axis=-1)
    safe = np.where(d > 0.0, d, 1.0)
    u = dr / safe[..., None]
    uu = u[..., None, :] * u[..., :, None]
    eye = np.eye(3)
    g = np.exp(1j * k * safe) / (4.0 * np.pi * safe)
    out = g[..., None, None] * (
        k * k * (eye - uu)
        + (1j * k / safe - 1.0 / safe**2)[..., None, None] * (eye - 3.0 * uu)
    )
    out[d == 0] = 0.0
    return out


class FFTDyadicKernel:
    """Zero-padded FFT implementation of the translation-invariant dyadic kernel."""

    def __init__(
        self,
        centers: Sequence[Sequence[float]],
        radii: Sequence[float],
        spacing: float,
        k: float,
        fill_quadrature: int = 1,
    ) -> None:
        self.centers = np.asarray(centers, dtype=float)
        self.radii = np.asarray(radii, dtype=float)
        self.spacing = float(spacing)
        self.k = float(k)
        self.fill_quadrature = int(fill_quadrature)

        lo, hi, rect_shape, full_points = bounding_grid(
            self.centers, self.radii, self.spacing
        )
        self.lo = lo
        self.hi = hi
        self.rect_shape = rect_shape
        self.pad_shape = tuple(2 * n for n in rect_shape)
        self.full_points = full_points
        self.n_rect = int(np.prod(rect_shape))
        self.n_pad = int(np.prod(self.pad_shape))

        points, labels, fill = voxelize(
            self.centers,
            self.radii,
            self.spacing,
            self.fill_quadrature,
        )
        self.points = points
        self.labels = labels
        self.fill = fill
        self.n_vox = int(points.shape[0])

        index = np.rint(
            (points - (lo + 0.5) * self.spacing) / self.spacing
        ).astype(np.intp)
        if index.size and (
            np.any(index < 0) or np.any(index >= np.asarray(rect_shape))
        ):
            raise ValueError("active voxel outside embedding rectangle")
        self.active_flat = np.ravel_multi_index(index.T, rect_shape)
        self.active_unravel = np.unravel_index(
            self.active_flat, rect_shape
        )
        if not np.array_equal(points, full_points[self.active_flat]):
            raise RuntimeError("active voxel index embedding failed")

        values = dyadic_kernel_values(
            self.rect_shape, self.pad_shape, self.spacing, self.k
        )
        self.kernel_values = values
        # spectra[field_component][source_component] on the padded grid
        self.spectra = [
            [
                np.fft.fftn(values[..., i, j], axes=(0, 1, 2))
                for j in range(3)
            ]
            for i in range(3)
        ]

    @property
    def kernel_spectra_bytes(self) -> int:
        return 9 * self.n_pad * np.dtype(complex).itemsize

    def convolve(self, source: np.ndarray) -> np.ndarray:
        """Linear dyadic convolution on the full rectangular grid.

        ``source`` has shape ``rect_shape + (3,)`` and gives the three vector
        components on every rectangular-grid cell.  The returned array has
        shape ``rect_shape + (3,)`` and equals the exact dyadic kernel action
        on all rectangular cells (occupied and empty).
        """
        if source.shape != self.rect_shape + (3,):
            raise ValueError(
                f"source has shape {source.shape}; expected "
                f"{self.rect_shape + (3,)}"
            )
        uhat = [
            np.fft.fftn(source[..., j], s=self.pad_shape, axes=(0, 1, 2))
            for j in range(3)
        ]
        result = []
        for i in range(3):
            acc = (
                self.spectra[i][0] * uhat[0]
                + self.spectra[i][1] * uhat[1]
                + self.spectra[i][2] * uhat[2]
            )
            full = np.fft.ifftn(acc, axes=(0, 1, 2))
            result.append(full[tuple(slice(0, n) for n in self.rect_shape)])
        return np.stack(result, axis=-1)

    def to_full_grid(self, vector: np.ndarray) -> np.ndarray:
        """Embed an occupied-grid vector ``(P,3)`` into ``rect_shape+(3,)``."""
        grid = np.zeros(self.rect_shape + (3,), dtype=complex)
        grid[self.active_unravel] = vector
        return grid

    def from_full_grid(self, grid: np.ndarray) -> np.ndarray:
        """Extract an occupied-grid vector ``(P,3)`` from a full grid."""
        return grid[self.active_unravel]


class FFTVIE:
    """Matrix-free FFT dipole VIE with GMRES for physical currents.

    This is the numerical replacement for ``DipoleVIE`` at occupancies where
    the parent dense LU factorization is not tractable.  The discretized
    operator and right-hand sides are byte-for-byte the parent definitions.
    """

    def __init__(
        self,
        centers: Sequence[Sequence[float]],
        radii: Sequence[float],
        spacing: float,
        k: float,
        fill_quadrature: int = 6,
        gmres_rtol: float = 1e-9,
        gmres_atol: float = 0.0,
        restart: int = 60,
        maxiter: int | None = None,
        wall_limit_seconds: float | None = None,
    ) -> None:
        self.kernel = FFTDyadicKernel(
            centers, radii, spacing, k, fill_quadrature=fill_quadrature
        )
        self.k = float(k)
        self.spacing = float(spacing)
        self.eps_shape = None
        self._alpha_values = None
        self._alpha_components = None
        self.incident = self._build_incident()
        self.ndof = 3 * self.kernel.n_vox
        self.gmres_rtol = float(gmres_rtol)
        self.gmres_atol = float(gmres_atol)
        self.restart = int(restart)
        self.maxiter = (
            int(maxiter) if maxiter is not None else min(1500, 5 * self.ndof)
        )
        self.wall_limit_seconds = wall_limit_seconds

        self.solve_started = 0.0
        self.matvecs = 0
        self.true_check_matvecs = 0
        self.iterations: list[int] = []
        self.residual_history: list[list[float]] = []
        self.solve_wall_seconds: list[float] = []
        self.gmres_flags: list[int] = []
        self.true_residuals: list[float] = []
        self.true_relative_residuals: list[float] = []

    @property
    def centers(self) -> np.ndarray:
        return self.kernel.centers

    @property
    def radii(self) -> np.ndarray:
        return self.kernel.radii

    @property
    def points(self) -> np.ndarray:
        return self.kernel.points

    @property
    def labels(self) -> np.ndarray:
        return self.kernel.labels

    @property
    def fill(self) -> np.ndarray:
        return self.kernel.fill

    @property
    def rect_shape(self) -> tuple[int, int, int]:
        return self.kernel.rect_shape

    @property
    def pad_shape(self) -> tuple[int, int, int]:
        return self.kernel.pad_shape

    def _build_incident(self) -> np.ndarray:
        points = self.kernel.points
        return np.stack(
            [
                np.exp(1j * self.k * (points @ direction))[:, None] * pol
                for direction, pol in illuminations()
            ],
            axis=-1,
        ).reshape(3 * len(points), -1)

    def set_eps(self, eps: Sequence[complex]) -> None:
        """Compute parent-identical polarizability for the material state."""
        e = np.asarray(eps, dtype=complex)[self.kernel.labels]
        a0 = 3.0 * self.spacing**3 * self.kernel.fill * (e - 1.0) / (e + 2.0)
        alpha = a0 / (1.0 - 1j * self.k**3 * a0 / (6.0 * np.pi))
        self.eps_shape = tuple(eps)
        self._alpha_values = alpha
        self._alpha_components = np.repeat(alpha, 3)

    def _check_wall_limit(self) -> None:
        if self.wall_limit_seconds is None or self.solve_started == 0.0:
            return
        if time.perf_counter() - self.solve_started > self.wall_limit_seconds:
            raise RuntimeError(
                f"single GMRES solve exceeded {self.wall_limit_seconds:g} s"
            )

    def apply_operator(self, x: np.ndarray, count: bool = True) -> np.ndarray:
        """Matrix-free action of ``A = I - diag(alpha) K`` on flat currents."""
        if count:
            self.matvecs += 1
            self._check_wall_limit()
        x = np.asarray(x, dtype=complex)
        xr = x.reshape(-1, 3)
        embedded = self.kernel.to_full_grid(xr)
        convolved = self.kernel.convolve(embedded)
        kv = self.kernel.from_full_grid(convolved).reshape(-1)
        return x - self._alpha_components * kv

    def solve(self, eps: Sequence[complex]) -> np.ndarray:
        """Solve all four physical-illumination current vectors with GMRES."""
        if (
            self._alpha_values is None
            or self.eps_shape != tuple(eps)
        ):
            self.set_eps(eps)
        b_all = self._alpha_components[:, None] * self.incident
        solution = np.empty((self.ndof, 4), dtype=complex)
        operator = LinearOperator(
            (self.ndof, self.ndof),
            matvec=lambda x: self.apply_operator(x, count=True),
            dtype=complex,
        )
        for column in range(4):
            b = b_all[:, column]
            hist: list[float] = []
            iters = {"n": 0}

            def callback(rnorm: float, *args: Any) -> None:
                iters["n"] += 1
                hist.append(float(rnorm))

            self.solve_started = time.perf_counter()
            t0 = self.solve_started
            x, flag = gmres(
                operator,
                b,
                rtol=self.gmres_rtol,
                atol=self.gmres_atol,
                restart=self.restart,
                maxiter=self.maxiter,
                callback=callback,
                callback_type="pr_norm",
            )
            self.solve_wall_seconds.append(time.perf_counter() - t0)
            self.gmres_flags.append(int(flag))
            self.iterations.append(int(iters["n"]))
            self.residual_history.append(list(hist))
            solution[:, column] = x

            # Independent true-residual audit using the same matrix-free
            # operator (one extra matvec per RHS, accounted separately).
            self.true_check_matvecs += 1
            true_res = b - self.apply_operator(x, count=False)
            self.true_residuals.append(float(np.linalg.norm(true_res)))
            self.true_relative_residuals.append(
                float(np.linalg.norm(true_res) / np.linalg.norm(b))
            )
        self.solve_started = 0.0
        return solution

    def field(self, eps: Sequence[complex], rx: np.ndarray) -> np.ndarray:
        """Scattered field at ``rx`` for all four illuminations."""
        p = self.solve(eps)
        return (
            dipole_kernel(rx, self.kernel.points, self.k) @ p
        ).reshape(len(rx), 3, -1)

    def memory_estimate(self) -> dict[str, float]:
        """Conservative peak-memory estimate in bytes, by component."""
        npad = self.kernel.n_pad
        nrect = self.kernel.n_rect
        ndof = self.ndof
        complex_bytes = np.dtype(complex).itemsize
        kernel_spectra = 9 * npad * complex_bytes
        # source copies, forward spectra, spectral accumulators, inverse
        # results, and FFT internal scratch are all padded-size arrays.
        fft_working = 6 * npad * complex_bytes
        rect_arrays = 3 * nrect * complex_bytes
        gmres_vectors = (self.restart + 4) * ndof * complex_bytes
        total = kernel_spectra + fft_working + rect_arrays + gmres_vectors
        return {
            "kernel_spectra_bytes": float(kernel_spectra),
            "fft_working_bytes": float(fft_working),
            "rect_arrays_bytes": float(rect_arrays),
            "gmres_vectors_bytes": float(gmres_vectors),
            "estimated_peak_bytes": float(total),
        }

    def solve_stats_row(self, eps: Sequence[complex]) -> dict[str, Any]:
        """Compact per-solve statistics after :meth:`solve`."""
        return {
            "scene_eps": [complex(e) for e in eps],
            "voxels": self.kernel.n_vox,
            "dof": self.ndof,
            "rect_shape": list(self.rect_shape),
            "pad_shape": list(self.pad_shape),
            "padded_cells": self.kernel.n_pad,
            "gmres_rtol": self.gmres_rtol,
            "restart": self.restart,
            "matvecs": self.matvecs,
            "true_check_matvecs": self.true_check_matvecs,
            "iterations": list(self.iterations),
            "residual_history": [
                [float(value) for value in history]
                for history in self.residual_history
            ],
            "gmres_flags": list(self.gmres_flags),
            "solve_wall_seconds": [float(x) for x in self.solve_wall_seconds],
            "total_solve_wall_seconds": float(sum(self.solve_wall_seconds)),
            "true_relative_residuals": [
                float(x) for x in self.true_relative_residuals
            ],
            "gmres_converged": all(
                flag == 0 for flag in self.gmres_flags
            ),
            "memory_estimate": {
                key: float(value)
                for key, value in self.memory_estimate().items()
            },
        }


def relative_field_error(
    predicted: np.ndarray, reference: np.ndarray
) -> float:
    """Frobenius relative field error over receivers/components/illuminations."""
    return float(
        np.linalg.norm(np.asarray(predicted) - np.asarray(reference))
        / np.linalg.norm(np.asarray(reference))
    )


def pointwise_errors(
    predicted: np.ndarray, reference: np.ndarray
) -> dict[str, Any]:
    """Preserve every pointwise complex field error with its shape."""
    err = np.asarray(predicted) - np.asarray(reference)
    return {
        "shape": list(err.shape),
        "errors_re_im": [
            [float(z.real), float(z.imag)] for z in err.reshape(-1)
        ],
        "max_abs": float(np.max(np.abs(err))),
        "rms": float(np.sqrt(np.mean(np.abs(err) ** 2))),
    }


__all__ = [
    "FFTDyadicKernel",
    "FFTVIE",
    "bounding_grid",
    "dyadic_kernel_values",
    "pointwise_errors",
    "relative_field_error",
]
