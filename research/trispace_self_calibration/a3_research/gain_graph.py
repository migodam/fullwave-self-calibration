"""Raw-likelihood local gain quotient for sparse complex transfer measurements.

Graph incidence/closure is established mathematics; this module claims no novelty.
No noisy cross-ratios or fitted-residual model-error certificates are used.
"""
import numpy as np
from scipy.linalg import solve_triangular


def complex_real_matrix(a):
    a = np.asarray(a)
    return np.block([[a.real, -a.imag], [a.imag, a.real]])


def real_parameter_jacobian(a):
    a = np.asarray(a)
    return np.concatenate((a.real, a.imag), axis=0)


def incidence(edges):
    edges = list(edges)
    vertices = sorted({('r', int(r)) for r, t in edges} |
                      {('t', int(t)) for r, t in edges})
    index = {v: i for i, v in enumerate(vertices)}
    b = np.zeros((len(edges), len(vertices)))
    for i, (r, t) in enumerate(edges):
        b[i, index['r', int(r)]] = 1.
        b[i, index['t', int(t)]] = -1.
    return b


def span(a, relative_tolerance=1e-10):
    a = np.asarray(a, dtype=float)
    if a.shape[1] == 0:
        return np.zeros((a.shape[0], 0))
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    if len(s) == 0 or s[0] == 0:
        return u[:, :0]
    return u[:, s > relative_tolerance*s[0]]


def visible_geometry(h, edges, geometry, material=None, covariance=None,
                     electronic=None, target_scale=None):
    """Return all target singular values, including nonidentifiable directions.

    geometry/material/electronic derivatives have REAL parameter columns.
    Gain parameters are complex and unrestricted at this local reference.
    covariance is real noise covariance in [Re observations, Im observations] order.
    target_scale maps dimensionless target increments to physical increments.
    """
    h = np.asarray(h, complex).reshape(-1)
    b = incidence(edges)
    if len(h) != b.shape[0]:
        raise ValueError('one field per edge required')
    gg = complex_real_matrix(h[:, None]*b)
    extras = [real_parameter_jacobian(v) for v in (material, electronic) if v is not None]
    nuisance = np.concatenate([gg]+extras, axis=1)
    target = real_parameter_jacobian(geometry)
    if target_scale is not None:
        target = target @ np.diag(np.asarray(target_scale))
    if covariance is not None:
        chol = np.linalg.cholesky(covariance)
        nuisance = solve_triangular(chol, nuisance, lower=True)
        target = solve_triangular(chol, target, lower=True)
    q = span(nuisance)
    visible = target-q@(q.T@target)
    vals = np.linalg.svd(visible, compute_uv=False)
    vals = np.pad(vals, (0, max(0, target.shape[1]-len(vals))))
    scale = max(np.linalg.norm(target, 2), np.finfo(float).tiny)
    rank = int(np.count_nonzero(vals > 1e-9*scale))
    return {'visible': visible, 'singular_values': vals,
            'target_rank': rank, 'target_dimension': target.shape[1],
            'nuisance_rank': q.shape[1],
            'gain_cycle_dimension_complex': int(len(h)-np.linalg.matrix_rank(b)),
            'nonzero_graph_assumption': bool(np.all(np.abs(h)>0))}


def branch_bound(predictions, beta, delta=.05):
    """Conditional bank guarantee; caller must separately establish coverage.

    Predictions are frozen REAL WHITENED vectors for independent new data.
    beta bounds deterministic prediction mismatch, not a fitted residual.
    """
    mu = np.asarray(predictions, float)
    if mu.ndim != 2 or len(mu) < 2 or beta < 0 or not 0 < delta < 1:
        raise ValueError('need >=2 candidate vectors, beta>=0, 0<delta<1')
    distances = np.linalg.norm(mu[:, None]-mu[None, :], axis=-1)
    np.fill_diagonal(distances, np.inf)
    risks = np.minimum(1., np.exp(-np.maximum(distances-2*beta, 0)**2/8).sum(axis=1))
    return {'worst_conditional_bound': float(risks.max()),
            'minimum_distance': float(distances.min()),
            'unlock_conditional_on_coverage': bool(risks.max() <= delta),
            'global_coverage_established': False}
