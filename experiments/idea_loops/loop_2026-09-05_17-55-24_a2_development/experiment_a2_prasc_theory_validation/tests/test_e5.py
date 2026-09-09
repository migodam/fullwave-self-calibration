"""Assertion runner for experiment E5 (pytest-compatible, self-runnable).

Run from the experiment root directly (pytest is not installed in the venv):

    python tests/test_e5.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import e5


def test_minimal_pair_values():
    rec = e5.minimal_pair()
    assert rec["pass"]
    assert abs(float(np.linalg.norm(rec["J1"]))) <= 1e-12
    assert abs(float(np.linalg.norm(rec["J2_contrast"]))) <= 1e-12
    assert np.isclose(rec["J_stack_contrast"].flat[0], 2.0, atol=1e-12)
    assert np.isclose(rec["sum_J_contrast"].flat[0], 0.0, atol=1e-12)
    assert rec["J_stack_contrast"].flat[0] > rec["sum_J_contrast"].flat[0]
    assert abs(rec["J_stack_duplicate_fro"]) <= 1e-12


def test_innovation_matches_direct_stack():
    rec = e5.random_innovation()
    assert rec["pass"]
    assert rec["formula_direct_max_residual"] <= 1e-12
    assert rec["formula_direct_fro_residual"] <= 1e-12
    assert rec["equality_V_fro"] <= 1e-10
    assert rec["equality_J_change_fro"] <= 1e-12
    assert rec["ker_J_new_dim"] == rec["ker_J_cap_ker_V_dim"]
    assert rec["ker_subspace_distance"] <= 1e-10
    sweep = e5.random_innovation_sweep()
    assert sweep["worst_fro_residual"] <= 1e-12


def test_singular_fallback_psd_and_direct_match():
    rec = e5.singular_fallback()
    assert rec["pass"]
    assert rec["rank_A"] == 2 and rec["rank_G"] == 2
    assert rec["G_min_eig"] <= 1e-10
    assert rec["J_new_psd"] >= -1e-10
    # Direct == variational/augmented stack; the naive pinv formula is wrong.
    assert rec["naive_pinv_wrong_fro_difference"] > 1e-10
    assert rec["naive_pinv_formula_invalid"]


def test_compensation_criterion():
    rec = e5.compensation_criterion()
    assert rec["pass"]
    assert rec["necessary_bound_ok"]
    assert rec["visible"]["stacked_visible"]
    assert rec["visible"]["full_col_rank_D"]
    assert rec["visible"]["max_per_frame_fro"] <= 1e-10
    assert not rec["hidden"]["stacked_visible"]
    assert rec["hidden"]["rank_D"] == 0
    assert rec["hidden"]["max_per_frame_fro"] <= 1e-10


def test_stack_vs_sum_psd_and_equality():
    rec = e5.stack_vs_sum()
    assert rec["pass"]
    assert rec["random_three_frame"]["min_eig_stack_minus_sum"] >= -1e-10
    assert rec["equality_trivial_common_c"]["gap_fro"] <= 1e-12
    nt = rec["equality_nontrivial_orthogonal_residual"]
    assert nt["min_eig_stack_minus_sum"] >= -1e-10
    assert nt["fro_gap"] <= 1e-10
    assert nt["J_stack_min_eig"] > 1e-8  # nontrivial equality, nonzero J
    assert nt["max_map_direction_deviation"] <= 1e-10


def test_saturated_new_frame_no_innovation():
    rec = e5.saturated_new_frame()
    assert rec["pass"]
    assert rec["a_clean_fro"] <= 1e-12
    assert rec["b_clean_fro"] <= 1e-12
    assert rec["J_change_fro"] <= 1e-10
    assert rec["V_fro"] <= 1e-12
    assert rec["I_acq_fro"] <= 1e-12


def test_rank_acquisition_budget_identity():
    rec = e5.rank_budget_seeds()
    assert rec["pass"]
    assert rec["worst_identity_fro_residual"] <= 1e-12
    for r in rec["records"]:
        assert r["identity_fro_residual"] <= 1e-12
        assert abs(r["min_eig_J_final_minus_J_original"]
                   - r["min_eig_I_acq_minus_L_rank"]) <= 1e-10


def test_crafted_budget_compensation_and_degradation():
    rec = e5.crafted_budget_cases()
    assert rec["pass"]
    assert rec["compensate"]["I_minus_L_min_eig"] >= -1e-8
    assert abs(rec["compensate"]["J_final_minus_original_min_eig"]) <= 1e-8
    assert rec["compensate"]["identity_fro"] <= 1e-12
    assert rec["noncompensate"]["I_minus_L_min_eig"] < -1e-8
    assert rec["noncompensate"]["J_final_minus_original_min_eig"] < -1e-8
    assert rec["noncompensate"]["identity_fro"] <= 1e-12


def test_neutral_rank_budget_control():
    rec = e5.neutral_budget_control()
    assert rec["pass"]
    assert rec["neutral_enlargement"]
    assert rec["L_rank_max_eig"] <= 1e-10
    assert rec["I_acq_min_eig"] >= -1e-10
    assert rec["min_eig_J_final_minus_J_original"] >= -1e-10
    assert rec["identity_fro_residual"] <= 1e-12


def test_nonsubmodularity_inequality():
    rec = e5.nonsubmodularity()
    assert rec["pass"]
    for r in rec["records"]:
        assert r["submodularity_violated"]
        assert r["g1_plus_g2"] < r["g12_plus_g_empty"]
        assert r["submodularity_gap"] > 1e-12


def test_acquisition_policies_run():
    rec = e5.acquisition_policies()
    assert rec["pass"]
    assert len(rec["greedy"]["subset"]) == 3
    assert len(rec["pair_lookahead"]["subset"]) == 3
    assert len(rec["exhaustive"]["best_subset"]) == 3
    assert len(rec["random"]["rows"]) == 12
    assert rec["exhaustive"]["evaluations"] == 56
    assert rec["greedy"]["logdet"] <= rec["exhaustive"]["best_logdet"] + 1e-9
    assert rec["pair_lookahead"]["logdet"] <= \
        rec["exhaustive"]["best_logdet"] + 1e-9


def test_global_gauge_null_persists():
    rec = e5.global_gauge_null()
    assert rec["pass"]
    assert rec["records"][0]["min_eig_J"] <= 1e-10   # L=2
    assert rec["records"][1]["min_eig_J"] <= 1e-10   # L=4
    assert rec["records"][0]["tangent_null_dim"] == 1
    assert rec["records"][1]["tangent_null_dim"] == 1


def test_anchor_restores_absolute():
    rec = e5.anchor_restores_absolute()
    assert rec["pass"]
    for r in rec["records"]:
        assert r["min_eig_J"] > 1e-10
        assert r["tangent_null_dim"] == 0


def test_source_stabilizer_is_not_claimed_passed():
    rec = e5.source_stabilizer_note()
    assert rec["status"] == "not_yet_run"
    assert rec["pass"] is None
    assert len(rec["reason"]) > 20


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
