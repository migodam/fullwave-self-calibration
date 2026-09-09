"""Unit tests for the corrected Family-1 Helmholtz conventions.

Run from the experiment root:
    .venv/bin/python -m unittest discover -s src -p "test_*.py" -v

Tests:
  1. self_cell_green equals a direct numerical radial (disk) quadrature of
     (i/4) H_0^{(1)}(k_b r) for many cell sizes;
  2. the complete I_self tends to zero as h -> 0 (the old rule, which omitted
     the -1/k_b^2 lower endpoint, saturates near 1/k_b^2 and fails this test);
  3. green_grad_source (used for the incident-field source derivative) matches
     a centred finite difference of g(z, s) in the source coordinate s -- an
     implementation with the old sign fails this at relative error ~2;
  4. grad_s g = -grad_z g for the implemented helpers;
  5. a pose finite difference of the full forward map against build_AB B,
     which would also fail if the source-gradient sign regressed.

Deterministic seeds: 12345 (source-gradient samples), 54321 (pose FD).
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.special import hankel1

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import helmholtz as hh  # noqa: E402

KB = 2.0 * np.pi
SELF_CELL_NS = [8, 12, 16, 24, 32, 40, 64, 128]
QUAD_NODES = 2048


def _disk_quadrature(k_b: float, h: float, nodes: int = QUAD_NODES) -> complex:
    """2 pi int_0^a r (i/4) H_0^{(1)}(k_b r) dr by Gauss-Legendre on [0, a]."""
    a = float(h) / np.sqrt(np.pi)
    x, w = np.polynomial.legendre.leggauss(nodes)
    r = 0.5 * a * (x + 1.0)          # map [-1, 1] -> [0, a]
    integrand = (1j / 4.0) * (2.0 * np.pi) * r * hankel1(0, k_b * r)
    return complex(np.dot(w, integrand) * (0.5 * a))


class TestSelfCellGreen(unittest.TestCase):
    """Corrected complete equal-area disk integral I_self."""

    def test_matches_radial_disk_quadrature(self) -> None:
        for N in SELF_CELL_NS:
            h = 1.0 / N
            analytic = hh.self_cell_green(KB, h)
            numeric = _disk_quadrature(KB, h)
            self.assertLess(
                abs(analytic - numeric),
                1e-11,
                msg=f"N={N}: analytic {analytic} vs quadrature {numeric}",
            )

    def test_complete_formula_tends_to_zero_under_refinement(self) -> None:
        mags = [abs(hh.self_cell_green(KB, 1.0 / N)) for N in SELF_CELL_NS]
        # Old (invalidated) rule gives |I| ~ 1/k_b^2 ~ 0.025 for every N here;
        # the complete integral starts near 5e-3 at N=8 and decays to ~4e-5.
        self.assertLess(max(mags), 1e-2)
        self.assertLess(mags[-1], 0.1 * mags[0])
        for coarse, fine in zip(mags, mags[1:]):
            self.assertLess(
                fine, coarse,
                msg="|I_self| must decrease over the refinement ladder",
            )
        # Domain-propagator diagonal k_b^2 I_self -> 0 (not -> 1).
        self.assertLess(abs(KB ** 2 * hh.self_cell_green(KB, 1.0 / 128.0)), 1e-2)


class TestSourceGradientSign(unittest.TestCase):
    """grad_s g = -grad_z g = +(i k_b/4) H_1(k_b R)(z - s)/R."""

    def test_green_grad_source_matches_finite_difference_in_s(self) -> None:
        rng = np.random.default_rng(12345)
        delta = 1e-6
        worst = 0.0
        samples = 0
        for _ in range(64):
            z = rng.uniform(-0.45, 0.45, 2)
            s = rng.uniform(-0.45, 0.45, 2)
            if np.linalg.norm(z - s) < 0.08:
                continue
            axis = int(rng.integers(0, 2))
            ep = np.zeros(2)
            ep[axis] = delta
            gp = hh.green_matrix((s + ep)[None, :], z[None, :], KB)[0, 0]
            gm = hh.green_matrix((s - ep)[None, :], z[None, :], KB)[0, 0]
            fd = (gp - gm) / (2.0 * delta)
            analytic = hh.green_grad_source(z[None, :], s, KB)[0, axis]
            denom = abs(analytic)
            if denom == 0.0:
                continue
            worst = max(worst, abs(fd - analytic) / denom)
            samples += 1
        self.assertGreater(samples, 20)
        # Correct sign errors ~1e-9 at delta=1e-6; old (wrong) sign ~2.0.
        self.assertLess(worst, 1e-5)

    def test_grad_source_equals_minus_grad_z(self) -> None:
        rng = np.random.default_rng(2026)
        worst = 0.0
        for _ in range(24):
            z = rng.uniform(-0.45, 0.45, (12, 2))
            s = rng.uniform(-0.45, 0.45, 2)
            if np.min(np.linalg.norm(z - s, axis=1)) < 0.05:
                continue
            grad_s = hh.green_grad_source(z, s, KB)
            grad_z = hh.green_grad_first(s[None, :], z, KB)[:, 0, :]
            scale = max(np.linalg.norm(grad_s), np.finfo(float).eps)
            worst = max(worst, float(np.linalg.norm(grad_s + grad_z) / scale))
        self.assertLess(worst, 1e-12)

    def test_pose_fd_would_fail_with_old_sign(self) -> None:
        """Full-map pose FD vs B: guards sign integration in build_AB."""
        from family1_pilot import CONFIG, build_poses, make_chi0

        cfg = dict(CONFIG)
        N = 16
        points, _ = hh.make_grid(N)
        chi0 = make_chi0(points, cfg)
        poses = build_poses(cfg)
        X0 = poses.reshape(-1)
        rx_offsets = np.asarray(cfg["rx_offsets"], dtype=float)
        tx_offset = np.asarray(cfg["tx_offset"], dtype=float)
        _, B, _, _ = hh.build_AB(
            chi0, poses, rx_offsets, tx_offset, N, cfg["k_b"]
        )
        rng = np.random.default_rng(54321)
        dX = rng.standard_normal(3 * cfg["T"])
        dX /= np.linalg.norm(dX)
        step = 1e-4
        Fp = hh.forward_measurements(
            chi0, (X0 + step * dX).reshape(cfg["T"], 3),
            rx_offsets, tx_offset, N, cfg["k_b"],
        )
        Fm = hh.forward_measurements(
            chi0, (X0 - step * dX).reshape(cfg["T"], 3),
            rx_offsets, tx_offset, N, cfg["k_b"],
        )
        fd = (Fp - Fm) / (2.0 * step)
        ref = B @ dX
        rel = float(
            np.linalg.norm(fd - ref) / max(np.linalg.norm(ref), np.finfo(float).eps)
        )
        # A sign regression in the source-gradient term gives relative errors
        # of order 1; the correct implementation is truncation-limited ~1e-10.
        self.assertLess(rel, 1e-4)


if __name__ == "__main__":
    unittest.main()
