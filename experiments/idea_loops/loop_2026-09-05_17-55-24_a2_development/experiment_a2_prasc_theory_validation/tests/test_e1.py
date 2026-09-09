"""Unit tests for experiment E1 (pytest-compatible, self-runnable).

Run from the experiment root either with pytest or directly:

    python tests/test_e1.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import common
from a2val import e1


def test_orth_proj_and_null():
    rng = np.random.default_rng(7)
    X = rng.standard_normal((9, 4)) + 1j * rng.standard_normal((9, 4))
    P = common.orth_proj(X)
    assert np.allclose(P @ X, X, atol=1e-12)
    assert np.allclose(P @ P, P, atol=1e-12)
    assert np.allclose(P, P.conj().T, atol=1e-12)
    N = common.proj_null(X)
    assert np.allclose(N @ X, 0, atol=1e-12)
    assert np.allclose(P + N, np.eye(X.shape[0]), atol=1e-12)


def test_phase_only_jph_zero():
    j = e1.phaseless_fisher("phase_only", 1.0)
    assert abs(j - 0.0) < 1e-8, j


def test_scale_family_jph_one():
    j = e1.phaseless_fisher("scale_family", 0.0)
    assert abs(j - 1.0) < 1e-6, j


def test_phase_only_jcoh_two():
    j = e1.analytic_values("phase_only", 1.0)["j_coh_analytic"]
    assert abs(j - 2.0) < 1e-12, j


def test_total_field_contraction():
    a = e1.analytic_values("total_field_ref", np.pi / 2.0)
    assert 0.0 < a["j_ph_analytic"] < a["j_coh_analytic"]


def test_s_lambda_matches_logpdf_fd():
    # Compare phaseless_score_lambda with a central difference of the ncx2
    # logpdf with respect to the noncentrality at fixed scale.
    lam = 2.0
    scale = 0.5
    h = 1e-5
    for t in (0.05, 0.3, 1.0, 3.0, 7.5):
        w = t / scale
        analytic = float(e1.phaseless_score_lambda(w, lam))
        fd = (e1.stats.ncx2.logpdf(t, 2, nc=lam + h, scale=scale)
              - e1.stats.ncx2.logpdf(t, 2, nc=lam - h, scale=scale)) / (2 * h)
        assert abs(analytic - fd) < 1e-5, (t, analytic, fd)


def test_s_lambda_expectation_zero():
    rng = np.random.default_rng(42)
    lam = 2.0
    w = rng.noncentral_chisquare(2, lam, size=200_000)
    s = e1.phaseless_score_lambda(w, lam)
    se = np.std(s, ddof=1) / np.sqrt(s.size)
    assert abs(np.mean(s)) < 6 * se, (np.mean(s), se)


def test_scale_coherent_and_phaseless_empirical_fisher_one():
    rec = e1.mc_score_check("scale_family", 0.0, seed=101, n_samples=2000)
    assert abs(rec["j_ph_emp"] - 1.0) < 0.2
    assert abs(rec["j_coh_emp"] - 1.0) < 0.2
    assert rec["mc_contraction_ok"]


def test_mismatch_fisher_value():
    a = e1.analytic_values("mismatch_intensity_noise", np.pi / 2.0)
    expected = (2.0 * 0.7 * 1.0) ** 2 / 1.0e-3  # (d|mu|^2/dx)^2 / sigma_z2
    assert abs(a["j_mismatch_analytic"] - expected) < 1e-9
    assert a["j_mismatch_analytic"] > a["j_coh_analytic"]


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(list(globals().items())):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception as exc:  # noqa: BLE001 - assertion runner
                failures += 1
                print(f"FAIL {name}: {type(exc).__name__}: {exc}")
    print(f"\n{len([k for k in globals() if k.startswith('test_') and callable(globals()[k])]) - failures}/{len([k for k in globals() if k.startswith('test_') and callable(globals()[k])])} tests passed")
    raise SystemExit(1 if failures else 0)
