"""Assertion runner for experiment E3 (pytest-compatible, self-runnable).

Run from the experiment root directly (pytest is not installed in the venv):

    python tests/test_e3.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import e3


def test_exact_small_dual_spectra():
    rec = e3.exact_small_dual_spectra()
    assert rec["all_pass"]
    for label in ("pi/4", "pi/2"):
        r = rec["records"][label]
        assert r["map_residual"] <= 1e-10, (label, r["map_residual"])
        assert r["pose_residual"] <= 1e-10, (label, r["pose_residual"])
        assert np.allclose(r["map_spectrum"], r["expected_map_spectrum"],
                           atol=1e-10)
        assert np.allclose(r["pose_spectrum"], r["expected_pose_spectrum"],
                           atol=1e-10)
    p4 = rec["records"]["pi/4"]
    assert np.allclose(p4["c"], [np.cos(np.pi / 4), 0.0], atol=1e-12)
    assert p4["zero_multiplicity_predicted"] == 0
    assert p4["unit_multiplicity_map_predicted"] == 1
    assert p4["unit_multiplicity_pose_predicted"] == 1
    assert np.allclose(p4["map_spectrum"], [0.5, 1.0], atol=1e-12)
    assert np.allclose(p4["pose_spectrum"], [0.5, 1.0], atol=1e-12)
    p2 = rec["records"]["pi/2"]
    assert np.allclose(p2["map_spectrum"], [1.0, 1.0], atol=1e-12)
    assert p2["unit_multiplicity_map_predicted"] == 2


def test_random_nuisance_theorem3():
    rec = e3.theorem3_fixture_check()
    assert rec["pass"]
    assert rec["max_spectral_residual"] <= 1e-10
    assert rec["zero_multiplicity_predicted"] == rec["zero_multiplicity_observed"]
    assert rec["unit_map_multiplicity_predicted"] == rec["unit_map_multiplicity_observed"]
    assert rec["unit_pose_multiplicity_predicted"] == rec["unit_pose_multiplicity_observed"]
    assert rec["raw_spectra_differ_from_normalized"]


def test_rho_undefined_control():
    rec = e3.rho_undefined_control()
    assert rec["pass"]
    assert rec["undefined_map"] == 1
    assert rec["undefined_pose"] == 1
    assert len(rec["map_spectrum"]) == rec["support_map"] == 2
    assert len(rec["pose_spectrum"]) == rec["support_pose"] == 1
    assert np.allclose(rec["map_spectrum"], 1.0, atol=1e-12)
    assert np.allclose(rec["pose_spectrum"], 1.0, atol=1e-12)


def test_scaling_control():
    rec = e3.scaling_control()
    assert rec["pass"]
    for r in rec["records"]:
        assert max(r["map_abs_delta_vs_unscaled"],
                   r["pose_abs_delta_vs_unscaled"]) <= 1e-10
        assert abs(r["J_x_fro_ratio_to_s2"] - 1.0) <= 1e-10


def test_b_zero_control():
    rec = e3.b_zero_control()
    assert rec["pass"]
    assert rec["map_retention_all_one"]
    assert rec["pose_support_dim"] == 0
    assert rec["J_x_fro"] == 0.0


def test_pose_prior_loewner_sandwich():
    rec = e3.pose_prior_checks()
    assert rec["sandwich"]["K_e<=K_eL"]
    assert rec["sandwich"]["K_eL<=K0"]
    assert rec["sandwich"]["min_eig_K_eL-K_e"] >= -1e-10
    assert rec["sandwich"]["min_eig_K0-K_eL"] >= -1e-10


def test_pose_prior_monotonicity():
    rec = e3.pose_prior_checks()["monotonicity"]
    assert rec["pass"]
    assert min(rec["min_eig_differences"].values()) >= -1e-10


def test_pose_prior_variational():
    rec = e3.pose_prior_checks()["variational"]
    assert rec["pass"]
    assert rec["max_residual"] <= 1e-8


def test_pose_prior_augmented_angles():
    rec = e3.pose_prior_checks()["augmented_principal_angles"]
    assert rec["augmented_residual"] <= 1e-10
    assert rec["closed_form_vs_augmented_residual"] <= 1e-10


def test_pose_prior_lambda_zero_recovers_ke():
    rec = e3.pose_prior_checks()["lambda_zero_recovery"]
    assert rec["fro_residual_vs_K_e"] <= 1e-10


def test_theorem7_analytic_mc_variance():
    rec = e3.theorem7_checks()
    assert rec["fixture"]["rank_B_v"] == rec["fixture"]["dims"]["p"]
    assert rec["pass"]
    assert rec["mc"]["all_variance_within_3se"]
    assert rec["mc"]["all_bias_within_3se"]
    assert max(r["sample_cov_vs_Jx_inv_fro_rel"]
               for r in rec["mc"]["records"]) < 0.2


def test_theorem7_deterministic_bound():
    rec = e3.theorem7_checks()["deterministic_bound"]
    assert rec["bound_satisfied"]
    assert rec["worst_analytic_Mx_error"] <= rec["bound"] * (1 + 1e-12) + 1e-12


def test_adaptive_rank_caveat_labeled():
    rec = e3.adaptive_rank_caveat()
    assert rec["B_v_ranks_full"]
    assert 0.05 < rec["fraction_selected_N1"] < 0.95
    assert rec["fro_diff_pooled_vs_majority_inverse_info"] > 0.0


def test_confounding_sweep_and_crossover():
    rec = e3.confounding_sweep()
    assert rec["crossover_eps_param"] == 1.0
    for row in rec["records"]:
        assert row["truncation_better"] == (row["eps_param"] < 1.0)
        if row["eps_param"] == 1.0:
            assert np.isclose(row["full_ls_risk"], 2.0)
            assert np.isclose(row["truncated_risk"], 2.0)
        assert row["full_ls_risk"] == 1.0 + row["eps_param"] ** (-2)
    assert rec["map_observation"]["map_information_retention"] == 1.0
    assert rec["map_observation"]["pose_truncation_bias"] == 1.0


def test_small_residual_large_bias_numbers():
    slb = e3.confounding_sweep()["small_residual_large_bias"]
    assert np.isclose(slb["full_ls_risk_exact"], 1.0 + 1e-6)
    assert np.isclose(slb["truncated_mse_exact"], 1.0 + 1e-6)
    assert slb["truncation_bias"] == 1.0
    assert np.isclose(slb["residual_y2_under_truncation"], 1e-3)


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
