"""Deterministic algebra checks, not electromagnetic validation or novelty proof."""
import json
from pathlib import Path
import numpy as np


def estimator(j, w):
    return np.linalg.solve(j.T @ w @ j, j.T @ w)


def main():
    rng = np.random.default_rng(9330)
    results = []
    for seed in range(12):
        j = rng.normal(size=(16, 5))
        q = rng.normal(size=(16, 3))
        w = np.linalg.inv(np.eye(16) + q @ q.T)
        l = estimator(j, w)
        t = rng.normal(size=(2, 5))
        e = rng.normal(size=(16, 3))
        a = t @ l @ e
        _, singular, vh = np.linalg.svd(a, full_matrices=False)
        worst = np.linalg.norm(a @ vh[0]) ** 2
        assert np.allclose(worst, singular[0] ** 2)
        assert np.allclose(l @ j, np.eye(5))
        # Exact covariance from an equiprobable zero-mean finite noise law.
        noise = np.sqrt(16) * np.concatenate((np.eye(16), -np.eye(16)), axis=1)
        d = e[:, 0]
        errors = t @ l @ (d[:, None] + noise)
        empirical = np.mean(np.sum(errors * errors, axis=0))
        risk = np.linalg.norm(t @ l, 'fro') ** 2 + np.linalg.norm(t @ l @ d) ** 2
        assert np.allclose(empirical, risk, rtol=1e-12)
        shortcut = np.trace(t @ np.linalg.inv(j.T @ w @ j) @ t.T)
        actual = np.linalg.norm(t @ l, 'fro') ** 2
        results.append(dict(case=seed, risk_identity_error=float(abs(empirical-risk)),
                            variance_shortcut_error=float(abs(shortcut-actual))))

    # Exact safe suppression uses the profiled task tangent, not raw B.
    n = rng.normal(size=(12, 3))
    b = rng.normal(size=(12, 2))
    bv = b - n @ np.linalg.lstsq(n, b, rcond=None)[0]
    u = rng.normal(size=12)
    u -= bv @ np.linalg.lstsq(bv, u, rcond=None)[0]
    u /= np.linalg.norm(u)
    j = np.column_stack((b, n))
    original = estimator(j, np.eye(12))[:2]
    changed = estimator(j, np.eye(12) - .95 * np.outer(u, u))[:2]
    safe_error = float(np.linalg.norm(changed-original))
    assert safe_error < 1e-11

    b = np.array([[1.], [0.]])
    n = np.ones((2, 1))
    w = np.diag([1., 1./10.])
    information = float((b.T @ w @ b - b.T @ w @ n @
                         np.linalg.solve(n.T @ w @ n, n.T @ w @ b)).item())
    assert np.allclose(information, 1./11.)
    l = estimator(np.column_stack((b, n)), w)[0]
    assert np.allclose(l @ l, 2.)
    # Identical means and sharp midpoint two-world risk lower bound.
    h = np.array([1., -2.])
    b = rng.normal(size=(8, 2))
    assert np.array_equal(b @ h, np.zeros(8) + b @ h)
    bound = float(h @ h / 4.)
    midpoint = h / 2.
    assert np.allclose(max(midpoint @ midpoint, (midpoint-h) @ (midpoint-h)), bound)
    report = dict(passed=True, random_cases=results, safe_suppression_error=safe_error,
                  raw_orthogonality_counterexample_information=information,
                  actual_white_noise_variance=float(l@l), two_world_bound=bound,
                  scope='Fixed linear real-whitened identifiable models; no data-dependent certificate.')
    out = Path(__file__).resolve().parent / 'results/a3_3_task_risk_checks.json'
    out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
