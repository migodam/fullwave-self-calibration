#!/usr/bin/env python3
"""Unit/smoke checks for the A2 E4 solver package (no final seeds)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOLVER_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SOLVER_DIR))
sys.path.insert(0, str(ROOT / "research/delegated/a2_physics"))
sys.path.insert(0, str(ROOT / "research/trispace_self_calibration/a2_research"))

from physics import Config, Model  # noqa: E402


class A2SolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Model(Config(N=16, aperture="full"))
        rng = np.random.default_rng(11)
        cls.alpha = np.array([0.5, 0.45, 0.55, 0.5, 0.48, 0.52, 0.47, 0.5, 0.51])
        cls.x = np.array([0.01, -0.02, 0.03])
        fw = cls.model.forward(cls.alpha, cls.x, jacobian=False)
        cls.sigma = 0.005
        cls.y = fw["total"] + (cls.sigma / np.sqrt(2)) * (
            rng.standard_normal(fw["total"].shape)
            + 1j * rng.standard_normal(fw["total"].shape)
        )

    def test_final_seed_gate(self):
        from common import path_for_seed

        with self.assertRaises(RuntimeError):
            path_for_seed(1001)

    def test_reduced_full_basis_matches_exact(self):
        from common import CostLedger
        from reduced import reduced_evaluate

        m = self.model
        U = np.eye(m.n_cells, dtype=np.complex128)
        out = reduced_evaluate(
            m, self.alpha, self.x, U, (0,), CostLedger(), need_jacobian=True
        )
        exact = m.forward(self.alpha, self.x, [0], jacobian=True)
        self.assertLess(
            np.linalg.norm(out["total"] - exact["total"]) / np.linalg.norm(exact["total"]),
            1e-12,
        )
        self.assertLess(
            np.linalg.norm(out["A"] - exact["A"]) / np.linalg.norm(exact["A"]), 1e-11
        )
        self.assertLess(
            np.linalg.norm(out["B"] - exact["B"]) / np.linalg.norm(exact["B"]), 1e-11
        )

    def test_phaseless_gradient_finite_difference(self):
        from common import CostLedger
        from objective import intensity_nll_and_weights, phaseless_gradient

        m = self.model
        led = CostLedger()
        g = phaseless_gradient(
            m, self.alpha, self.x, self.y, self.sigma, (0,), led
        )
        ga = np.concatenate([g["grad_alpha"], g["grad_x"]])
        eps = 1e-6
        fd = np.zeros(12)
        for i in range(12):
            if i < 9:
                p = self.alpha.copy()
                p[i] += eps
                f1 = m.forward(p, self.x, [0], jacobian=False)["total"]
                l1 = intensity_nll_and_weights(f1, self.y[:72], self.sigma)[0]
                p[i] -= 2 * eps
                f2 = m.forward(p, self.x, [0], jacobian=False)["total"]
                l2 = intensity_nll_and_weights(f2, self.y[:72], self.sigma)[0]
            else:
                xp = self.x.copy()
                xp[i - 9] += eps
                f1 = m.forward(self.alpha, xp, [0], jacobian=False)["total"]
                l1 = intensity_nll_and_weights(f1, self.y[:72], self.sigma)[0]
                xp[i - 9] -= 2 * eps
                f2 = m.forward(self.alpha, xp, [0], jacobian=False)["total"]
                l2 = intensity_nll_and_weights(f2, self.y[:72], self.sigma)[0]
            fd[i] = (l1 - l2) / (2 * eps)
        self.assertLess(
            np.linalg.norm(ga - fd) / np.linalg.norm(fd), 2e-5
        )

    def test_passive_certificate_available(self):
        from cert import reduced_certification
        from common import CostLedger
        from reduced import build_sensing_basis, reduced_evaluate

        led = CostLedger()
        basis = build_sensing_basis(self.model, self.x, (0,), 8, led)
        out = reduced_evaluate(
            self.model,
            self.alpha,
            self.x,
            basis["U"],
            (0,),
            led,
            need_jacobian=True,
        )
        cert = reduced_certification(
            self.model,
            self.alpha,
            self.y[:72],
            self.sigma,
            (0,),
            out,
            {},
        )
        self.assertTrue(cert["available"])
        # conservative certificates are expected to refuse low-rank charts;
        # that refusal is reported, not silently replaced by true-error tests
        self.assertFalse(cert["pass"])

    def test_direct_solve_returns_finite_and_ledgered(self):
        from solver import solve

        res = solve(
            self.model,
            self.y,
            self.sigma,
            np.full(9, 0.4),
            np.zeros(3),
            "direct",
            200,
            {"maxls": 10, "ftol": 1e-8, "init_units": 0},
        )
        self.assertTrue(np.all(np.isfinite(res["alpha_est"])))
        self.assertTrue(np.all(np.isfinite(res["x_est"])))
        self.assertGreaterEqual(res["ledger"]["units"], 0)
        self.assertLessEqual(res["ledger"]["units"], 200.0 + 1e-6)
        self.assertEqual(res["ledger"]["full_rhs_solves"], res["ledger"]["units"])

    def test_intensity_nll_stable(self):
        from objective import intensity_nll_and_weights

        y = self.y[:72]
        mu = self.model.forward(self.alpha, self.x, [0], jacobian=False)["total"]
        nll, w = intensity_nll_and_weights(mu, y, self.sigma)
        self.assertTrue(np.isfinite(nll))
        self.assertTrue(np.all(np.isfinite(w)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
