"""Tangent-enabled matrix-free FFT dipole VIE adapter (material derivatives).

This module extends the read-only ``a3_maxwell_fft/maxwell_fft.py`` GMRES
solver with the parent ``DipoleVIE.currents(eps, derivatives=True)`` interface
so that the parent 3D calibration estimator can be driven with an injected FFT
adapter.  All electromagnetic conventions (time sign, ordering, voxelisation,
Clausius-Mossotti + radiative correction, subcell fill, physical currents with
``r_free=0``) are the read-only parent/FFT definitions; nothing here modifies
parent or a3_maxwell_fft sources.

The tangent right-hand side is the analytic derivative of the discrete
physical-current equation ``(I - diag(alpha) K) p = alpha E_inc`` for each
known region ``j``:

    A dp_j = d(alpha_j)/d(eps_j) * mask_j * (E_inc + K p),

with the actual occupied-cell polarizability (fill and radiative correction)
differentiated.  Tangent columns are solved with the same matrix-free operator
and are verified against the dense parent and finite differences by the test
suite.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

_FFT_DIR = Path(__file__).resolve().parent.parent / "a3_maxwell_fft"
sys.path.insert(0, str(_FFT_DIR))

from maxwell_fft import FFTVIE  # noqa: E402  (read-only base)

_PARENT_DIR = (
    Path(__file__).resolve().parents[3]
    / "research"
    / "trispace_self_calibration"
    / "a3_research"
)
sys.path.insert(0, str(_PARENT_DIR))

import maxwell3d as _parent  # noqa: E402

dipole_kernel = _parent.dipole_kernel
voxelize = _parent.voxelize
illuminations = _parent.illuminations
receivers = _parent.receivers
treams_field = _parent.treams_field


class TangentFFTVIE(FFTVIE):
    """FFT GMRES VIE with parent-compatible ``currents(eps, derivatives=...)``.

    ``work`` is kept ledger-compatible with ``maxwell3d.DipoleVIE.work``
    (``factorizations``, ``rhs_columns``, ``wall_seconds``) so the parent
    estimator's work-ledger arithmetic is unchanged.  Additional matrix-free
    solve/matvec counters are exposed under ``fft_work`` and in
    :meth:`work_row` for separate accounting.
    """

    def __init__(
        self,
        centers: Sequence[Sequence[float]],
        radii: Sequence[float],
        spacing: float,
        k: float,
        fill_quadrature: int = 6,
        gmres_rtol: float = 1e-11,
        gmres_atol: float = 0.0,
        restart: int = 120,
        maxiter: int | None = None,
        wall_limit_seconds: float | None = None,
    ) -> None:
        super().__init__(
            centers,
            radii,
            spacing,
            k,
            fill_quadrature=fill_quadrature,
            gmres_rtol=gmres_rtol,
            gmres_atol=gmres_atol,
            restart=restart,
            maxiter=maxiter,
            wall_limit_seconds=wall_limit_seconds,
        )
        self.work: dict[str, float] = {
            "factorizations": 0.0,
            "rhs_columns": 0.0,
            "wall_seconds": 0.0,
        }
        self.fft_work: dict[str, float] = {
            "solves": 0.0,
            "matvecs": 0.0,
            "true_check_matvecs": 0.0,
            "tangent_field_convolutions": 0.0,
            "operator_updates": 0.0,
        }
        self.last_stats: dict[str, Any] = {}

    # -- column solves -----------------------------------------------------
    def _solve_columns(self, rhs: np.ndarray) -> np.ndarray:
        """Solve the current matrix-free operator for each column of ``rhs``.

        Returns a complex array with one solution per RHS column.  Every GMRES
        solve and every independent true-residual matvec is recorded in the
        same counter scheme as the base solver.
        """
        rhs = np.asarray(rhs, dtype=complex)
        if rhs.ndim == 1:
            rhs = rhs[:, None]
        n_cols = rhs.shape[1]
        solution = np.empty_like(rhs)
        operator = LinearOperator(
            (self.ndof, self.ndof),
            matvec=lambda x: self.apply_operator(x, count=True),
            dtype=complex,
        )
        for col in range(n_cols):
            hist: list[float] = []
            iters = {"n": 0}

            def callback(rnorm: float, *args: Any) -> None:
                iters["n"] += 1
                hist.append(float(rnorm))

            self.solve_started = time.perf_counter()
            t0 = self.solve_started
            x, flag = gmres(
                operator,
                rhs[:, col],
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
            solution[:, col] = x
            self.true_check_matvecs += 1
            true_res = rhs[:, col] - self.apply_operator(x, count=False)
            self.true_residuals.append(float(np.linalg.norm(true_res)))
            self.true_relative_residuals.append(
                float(np.linalg.norm(true_res) / max(np.linalg.norm(rhs[:, col]), 1e-300))
            )
            limit = max(10 * self.gmres_rtol, 1e-9)
            if flag != 0 or not np.isfinite(self.true_relative_residuals[-1]) or self.true_relative_residuals[-1] > limit:
                self.solve_started = 0.0
                raise RuntimeError(f"GMRES failed: flag={flag}, true relative residual={self.true_relative_residuals[-1]:.3g}")
        self.solve_started = 0.0
        return solution

    # -- parent-compatible adapter ----------------------------------------
    def currents(
        self, eps: Sequence[complex], derivatives: bool = False
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        """Parent-compatible currents and (optionally) material tangent.

        ``eps`` are complex region permittivities exactly as consumed by the
        parent adapter (e.g. real search value plus fixed imaginary loss).
        Returns ``p`` with shape ``(3*P, n_illuminations)`` and, when
        ``derivatives`` is true, ``(p, dp)`` where ``dp`` has shape
        ``(3*P, n_illuminations, len(eps))`` and orders the two trailing axes
        illumination-major then region, exactly as the dense parent does.
        """
        if self._alpha_values is None or self.eps_shape != tuple(eps):
            self.set_eps(eps)
            self.fft_work["operator_updates"] += 1
        start = time.perf_counter()
        before = (
            self.matvecs,
            self.true_check_matvecs,
            len(self.gmres_flags),
        )

        b_current = self._alpha_components[:, None] * self.incident
        p = self._solve_columns(b_current)
        n_tangent_columns = p.shape[1] * int(len(eps)) if derivatives else 0

        if derivatives:
            # Total field (incident + self-field) used by the analytic tangent.
            # p is (3*P, n_ill); each illumination is a (P,3) vector field in
            # point-major/component-minor order.
            per_ill = p.T.reshape(
                p.shape[1], self.kernel.n_vox, 3
            )
            kp_per_ill = np.stack(
                [
                    self.kernel.from_full_grid(
                        self.kernel.convolve(
                            self.kernel.to_full_grid(field)
                        )
                    )
                    for field in per_ill
                ],
                axis=0,
            )
            kp_flat = kp_per_ill.reshape(p.shape[1], -1).T
            self.fft_work["tangent_field_convolutions"] += p.shape[1]
            total = self.incident + kp_flat

            e = np.asarray(eps, dtype=complex)[self.kernel.labels]
            a0 = (
                3.0
                * self.spacing**3
                * self.kernel.fill
                * (e - 1.0)
                / (e + 2.0)
            )
            alpha = a0 / (1.0 - 1j * self.k**3 * a0 / (6.0 * np.pi))
            denom = (1.0 - 1j * self.k**3 * a0 / (6.0 * np.pi)) ** 2
            dalpha = (
                9.0 * self.spacing**3 * self.kernel.fill / (e + 2.0) ** 2
            ) / denom

            n_regions = int(len(eps))
            rhs = np.stack(
                [
                    np.repeat(dalpha * (self.kernel.labels == j), 3)[:, None]
                    * total
                    for j in range(n_regions)
                ],
                axis=-1,
            )
            dp = self._solve_columns(
                rhs.reshape(self.ndof, n_tangent_columns)
            )
            dp = dp.reshape(p.shape + (n_regions,))
        else:
            dp = None

        # Matrix-free GMRES performs no LU factorization.
        self.work["rhs_columns"] += float(
            p.shape[1] + n_tangent_columns
        )
        self.work["wall_seconds"] += time.perf_counter() - start
        self.fft_work["solves"] += float(
            len(self.gmres_flags) - before[2]
        )
        self.fft_work["matvecs"] += float(self.matvecs - before[0])
        self.fft_work["true_check_matvecs"] += float(
            self.true_check_matvecs - before[1]
        )
        self.solve_started = 0.0
        return (p, dp) if derivatives else p


    def fft_work_delta(self, before: dict[str, float]) -> dict[str, float]:
        """Deltas of the matrix-free counters for one model activity window."""
        return {
            key: float(self.fft_work[key] - before[key])
            for key in before
        }

    def work_row(self) -> dict[str, Any]:
        """Compact solve/memory accounting, excluding per-solve histories."""
        memory = self.memory_estimate()
        flags = list(self.gmres_flags)
        return {
            "work": dict(self.work),
            "fft_work": dict(self.fft_work),
            "solves": float(len(flags)),
            "matvecs": float(self.matvecs),
            "true_check_matvecs": float(self.true_check_matvecs),
            "iterations": list(self.iterations),
            "gmres_flags": flags,
            "gmres_converged": bool(flags and all(f == 0 for f in flags)),
            "true_relative_residuals": list(self.true_relative_residuals),
            "solve_wall_seconds": list(self.solve_wall_seconds),
            "memory_estimate_bytes": memory,
            "memory_estimate_mib": memory["estimated_peak_bytes"] / 2**20,
        }


def _model_fft_work(model: TangentFFTVIE) -> dict[str, float]:
    return dict(model.fft_work)


__all__ = [
    "TangentFFTVIE",
    "dipole_kernel",
    "receivers",
    "treams_field",
    "voxelize",
]
