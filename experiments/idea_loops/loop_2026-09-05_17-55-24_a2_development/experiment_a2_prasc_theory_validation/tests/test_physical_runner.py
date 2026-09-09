"""Light smoke test for the OPTIONAL physical second-stage runner.

Runs ONE N=8 full-aperture forward from the shared physical core and checks the
realified tangent shapes.  It deliberately does NOT run the 24 validation seeds
(that is the runner's job, not a unit test).

Run from the experiment root (pytest is not installed in the venv):

    python tests/test_physical_runner.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

PHYSICS_DIR = (
    "/Volumes/migodam's-external-brain/Research/Inv_SLAM/"
    "research/delegated/a2_physics"
)
if PHYSICS_DIR not in sys.path:
    sys.path.insert(0, PHYSICS_DIR)
import physics  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def test_physical_forward_realified_shapes():
    cfg = physics.Config(N=8, aperture="full")
    model = physics.Model(cfg)
    alpha = np.full(9, 0.2)
    x = np.zeros(3)
    out = model.forward(alpha, x, jacobian=True)
    assert out["A"].shape == (288, 9)
    assert out["B"].shape == (288, 3)
    A_r = physics.realify_jacobian(out["A"])
    B_r = physics.realify_jacobian(out["B"])
    assert A_r.shape == (576, 9)
    assert B_r.shape == (576, 3)
    assert out["work"]["rhs_solves_total"] > 0


if __name__ == "__main__":
    test_physical_forward_realified_shapes()
    print("test_physical_runner: OK")
