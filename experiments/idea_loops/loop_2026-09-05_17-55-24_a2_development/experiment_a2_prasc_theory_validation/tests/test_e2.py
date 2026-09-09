"""Assertion runner for experiment E2 (pytest-compatible, self-runnable).

Run from the experiment root directly (pytest is not installed in the venv):

    python tests/test_e2.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import e2


def test_exact_t5a_t5b():
    ec = e2.exact_controls()
    for name in ("equality_case", "strict_loss_case",
                 "complete_hiding_case", "loss_without_rank_drop"):
        rec = ec[name]
        assert rec["loss_residual_fro"] <= 1e-12, (name, rec["loss_residual_fro"])
        assert rec["rank_drop_predicted"] == rec["rank_drop_actual"], name


def test_exact_expected_values():
    ec = e2.exact_controls()
    assert np.allclose(ec["equality_case"]["j0"], np.eye(2), atol=1e-12)
    assert np.allclose(ec["equality_case"]["j1"], np.eye(2), atol=1e-12)
    assert np.allclose(ec["strict_loss_case"]["j0"], np.eye(2), atol=1e-12)
    assert np.allclose(ec["strict_loss_case"]["j1"], np.diag([0.0, 1.0]), atol=1e-12)
    assert np.allclose(ec["complete_hiding_case"]["j1"], np.zeros((2, 2)), atol=1e-12)
    assert abs(ec["loss_without_rank_drop"]["j0_scalar"] - 2.0) < 1e-12
    assert abs(ec["loss_without_rank_drop"]["j1_scalar"] - 1.0) < 1e-12


def test_nonmonotone_rho():
    rec = e2.exact_controls()["nonmonotone_rho"]
    expect = [0.75, 1.0, 0.0]
    for rr, val in zip(rec["rho_records"], expect):
        assert len(rr["rho"]) == 1
        assert abs(rr["rho"][0] - val) <= 1e-12, (rr["rho"], val)
    assert rec["rho_sequence"][0] < rec["rho_sequence"][1]
    assert rec["rho_sequence"][2] < rec["rho_sequence"][1]


def test_rank_uncertainty():
    rec = e2.exact_controls()["rank_uncertainty"]
    assert abs(rec["j_x_0"] - 1.0) < 1e-12
    assert all(abs(x) < 1e-12 for x in rec["j_x_nonzero"])


def test_saturation_zero_up_to_backward_error():
    rec = e2.exact_controls()["saturation_control"]
    assert rec["backward_scaled_residual"] < 100.0
    assert rec["j_x_fro"] <= 100.0 * rec["backward_scale"]
    # No relative-to-self rank claim: record the absolute singular values.
    assert all(sv < 1e-12 for sv in rec["bv_absolute_singular_values"])


def test_rank_threshold_control():
    rec = e2.exact_controls()["rank_threshold_control"]
    assert rec["rank_absolute"] == 2
    assert rec["rank_relative"] == 2
    assert abs(rec["ratios_to_sigma1"][2] - 1e-14) < 1e-30


def test_projector_gap_closure():
    rec = e2.exact_controls()["projector_gap_closure"]
    assert all(r["rank_top_one_projector"] == 1 for r in rec["records"])
    assert abs(rec["fro_P_plus_minus_0_05"] - np.sqrt(2.0)) < 1e-12
    assert abs(rec["left_right_limits"]["fro_distance"] - np.sqrt(2.0)) < 1e-12


def test_real_safe_space():
    rec = e2.real_safe_space()
    assert rec["dim_kerF"] == rec["bound"]["d_minus_rankF"]
    assert rec["bound"]["holds"]
    assert rec["safe_loss_fro"] <= 1e-10
    w = rec["diff_eigs"]
    assert min(w) >= -1e-10
    assert max(w) > 1e-10


def test_complex_safe_kernel():
    rec = e2.complex_safe_space()
    assert rec["jc_invariance_residual_fro"] <= 1e-10
    assert rec["max_abs_F_Z_S"] <= 1e-10
    assert rec["safe_loss_fro_after_S"] <= 1e-10
    assert rec["real_dim_S"] == 2
    assert rec["unsafe"]["decreases_loewner"]
    assert rec["unsafe"]["z_in_ker_F"]
    assert not rec["unsafe"]["in_ker_F_Jc"]


def test_utility_projection():
    rec = e2.utility_projection()
    assert rec["maximizer_dominates"]
    assert rec["trace_claimed_maximizer"] >= rec["max_random_trace"] - 1e-8


def test_random_tangents():
    rec = e2.random_tangent_checks()
    assert rec["n_checks"] == 36
    assert rec["all_t5b_matches"]
    assert rec["max_scaled_residual"] < 100.0


if __name__ == "__main__":
    failures = 0
    total = 0
    for name, fn in sorted(list(globals().items())):
        if name.startswith("test_") and callable(fn):
            total += 1
            try:
                fn()
                print(f"PASS {name}")
            except Exception as exc:  # noqa: BLE001 - assertion runner
                failures += 1
                print(f"FAIL {name}: {type(exc).__name__}: {exc}")
    print(f"\n{total - failures}/{total} tests passed")
    raise SystemExit(1 if failures else 0)
