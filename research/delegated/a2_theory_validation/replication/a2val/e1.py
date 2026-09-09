"""Experiment E1: matched coherent vs phaseless information contraction.

Theory
------

Parent observations are proper circular complex Gaussian

    Y ~ CN(mu(theta), sigma^2),      E |Y - mu|^2 = sigma^2,

so each of real and imaginary parts has variance sigma^2/2.  The phaseless
observation is the intensity I = |Y|^2, distributed as the Rice-squared /
noncentral chi-square law:

    I = (sigma^2/2) * chi2_2(lambda),
    lambda = 2 |mu|^2 / sigma^2   (noncentrality of the chi-square parent),

i.e. scipy.stats.ncx2(df=2, nc=lambda, scale=sigma^2/2).

For a real scalar parameter theta at fixed sigma^2 the coherent Fisher and
score are

    J_coh    = (2/sigma^2) |d mu/d theta|^2,
    s_coh(Y) = (2/sigma^2) Re[conj(Y - mu) (d mu/d theta)].

When sigma^2 = exp(theta) (scale family, mu == 0) the covariance-score term
must be added (never use only the mean-Jacobian Gram): with
u = d log(sigma^2)/d theta,

    s_coh(Y) = u * (|Y|^2/sigma^2 - 1),       J_coh = u^2.

For fixed scale and intensity t, write w = t/(sigma^2/2) = 2 t/sigma^2.  The
intensity score with respect to lambda is, for lambda > 0,

    s_lambda(t) = -1/2 + (1/2) sqrt(w/lambda)
                  * I1(sqrt(lambda w)) / I0(sqrt(lambda w)),

and the phaseless score for theta is s_ph(t) = s_lambda(t) * d lambda/d theta.
For lambda < 1e-12 the limit s_lambda -> -1/2 + w/4 is used (equivalently a
finite difference of the ncx2 logpdf would give the same derivative).  For the
scale family (lambda = 0 identically) s_ph is the derivative of the ncx2
logpdf with respect to theta, computed by a symmetric finite difference in the
scale/log-sigma^2 coordinate.

The exact phaseless Fisher is

    J_ph = integral_0^inf s_ph(t)^2 f_I(t) dt,

evaluated by scipy.integrate.quad on [ppf(1e-12), ppf(1 - 1e-12)] and checked
against a full-domain quadrature to confirm the tails are negligible.

Matched-model identity (Theorem 1 of the A2 package):

    s_ph(I) = E[s_coh(Y) | I],        J_ph <= J_coh,

with equality iff the coherent score is measurable with respect to I.
Experiment E1 checks this identity by equal-count binning of I.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, special, stats


# ---------------------------------------------------------------------------
# Model registry
# ---------------------------------------------------------------------------

#: Default model parameters (merged with per-call ``params`` overrides).
DEFAULTS = {
    "phase_only": {"a": 1.0, "sigma2": 1.0},
    "nuisance_phase": {"a": 1.0, "sigma2": 1.0, "eta": 0.7},
    "scale_family": {},  # theta = log sigma2, sigma2 = exp(theta)
    "total_field_ref": {"a": 1.0, "c": 0.7, "sigma2": 1.0},
    "mismatch_intensity_noise": {
        "a": 1.0,
        "c": 0.7,
        "sigma2": 1.0,
        "sigma_z2": 1.0e-3,
    },
}

FIXED_SCALE_MODELS = ("phase_only", "nuisance_phase", "total_field_ref", "mismatch_intensity_noise")
CONTROL_MODELS = ("phase_only", "nuisance_phase", "scale_family", "total_field_ref")


def _params(model: str, params=None) -> dict:
    d = dict(DEFAULTS[model])
    if params:
        d.update(params)
    return d


def model_mu(model: str, theta: float, params=None) -> complex:
    """Mean ``mu(theta)`` of the coherent parent for the named control model."""
    p = _params(model, params)
    if model == "scale_family":
        return 0.0 + 0.0j
    if model == "phase_only":
        return p["a"] * np.exp(1j * theta)
    if model == "nuisance_phase":
        return p["a"] * np.exp(1j * (theta + p["eta"]))
    if model in ("total_field_ref", "mismatch_intensity_noise"):
        return p["c"] + p["a"] * np.exp(1j * theta)
    raise ValueError(f"unknown model {model!r}")


def model_dmu(model: str, theta: float, params=None) -> complex:
    """Derivative of ``mu`` with respect to the target parameter theta."""
    p = _params(model, params)
    if model == "scale_family":
        return 0.0 + 0.0j
    if model == "phase_only":
        return 1j * p["a"] * np.exp(1j * theta)
    if model == "nuisance_phase":
        return 1j * p["a"] * np.exp(1j * (theta + p["eta"]))
    if model in ("total_field_ref", "mismatch_intensity_noise"):
        return 1j * p["a"] * np.exp(1j * theta)
    raise ValueError(f"unknown model {model!r}")


def model_dmu_nuisance(model: str, theta: float, params=None) -> complex:
    """Derivative of ``mu`` with respect to the nuisance phase eta (model 2)."""
    if model != "nuisance_phase":
        raise ValueError(f"nuisance derivative only defined for nuisance_phase, got {model!r}")
    # d/deta of a*exp(i*(x+eta)) equals d/dx of the same expression.
    return model_dmu(model, theta, params=params)


def model_sigma2(model: str, theta: float, params=None) -> float:
    p = _params(model, params)
    if model == "scale_family":
        return float(np.exp(theta))
    return float(p["sigma2"])


def model_logsigma2_deriv(model: str, theta: float, params=None) -> float:
    """u = d log(sigma2)/d theta (0 for fixed-scale mean models, 1 for scale)."""
    p = _params(model, params)
    if model == "scale_family":
        return 1.0
    _ = p
    return 0.0


def model_nc(model: str, theta: float, params=None) -> float:
    """Noncentrality lambda = 2|mu|^2/sigma^2 of the intensity law."""
    sigma2 = model_sigma2(model, theta, params=params)
    if sigma2 <= 0:
        raise ValueError("sigma2 must be positive")
    return float(2.0 * abs(model_mu(model, theta, params=params)) ** 2 / sigma2)


def model_dnc_dtheta(model: str, theta: float, params=None) -> float:
    """d lambda/d theta for fixed-scale mean models (d lambda/dx formula)."""
    sigma2 = model_sigma2(model, theta, params=params)
    mu = model_mu(model, theta, params=params)
    dmu = model_dmu(model, theta, params=params)
    # lambda = 2|mu|^2/sigma2 -> dlambda/dtheta = (4/sigma2) Re[conj(mu) dmu]
    value = float((4.0 / sigma2) * np.real(np.conj(mu) * dmu))
    # Snap derivatives that are zero in exact arithmetic (pure phase and
    # nuisance-phase models have Re[conj(mu)*i*mu] = 0) to exactly zero so the
    # phaseless score is identically zero up to floating-point representation.
    guard = 1e-12 * max(1.0, (4.0 / sigma2) * abs(mu) * abs(dmu))
    return 0.0 if abs(value) < guard else value


# ---------------------------------------------------------------------------
# Coherent experiment
# ---------------------------------------------------------------------------


def coherent_score(y, mu, dmu, sigma2=1.0, u=0.0):
    """Coherent score for the complex-Gaussian likelihood.

    Parameters
    ----------
    y : complex scalar or ndarray
        Observed coherent sample(s).
    mu, dmu : complex scalar or ndarray
        Mean and its derivative with respect to the target parameter.
    sigma2 : float
        Per-component variance scale (E|y-mu|^2 = sigma2).
    u : float
        d log(sigma2)/d theta when the covariance itself depends on theta.

    Returns
    -------
    ndarray / scalar with
        s_coh = u*(|y-mu|^2/sigma2 - 1)
                + (2/sigma2)*Re[conj(y-mu)*dmu].
    """
    y = np.asarray(y)
    mu = np.asarray(mu)
    dmu = np.asarray(dmu)
    resid = y - mu
    scale_term = u * (np.abs(resid) ** 2 / sigma2 - 1.0)
    mean_term = (2.0 / sigma2) * np.real(np.conj(resid) * dmu)
    return scale_term + mean_term


def coherent_fisher(mu, dmu, sigma2=1.0, u=0.0):
    """Analytic coherent Fisher for one scalar parameter.

    For fixed covariance this is J_coh = (2/sigma2)*|dmu|^2; for a pure scale
    model (dmu = 0) it is u^2.  The cross term between the circular-Gaussian
    mean and covariance scores vanishes, so the two contributions add.
    """
    mu = complex(mu)
    dmu = complex(dmu)
    return float(u * u + (2.0 / sigma2) * abs(dmu) ** 2)


# ---------------------------------------------------------------------------
# Phaseless (intensity) experiment
# ---------------------------------------------------------------------------


def _bessel_ratio_small(z):
    """I1(z)/I0(z) from the convergent series, safe for very small z."""
    # I1/I0 = z/2 * (1 - z^2/8 + z^4/192 - z^6/3072 + z^8/46080 - ...)
    z2 = z * z
    series = 1.0 - z2 / 8.0 + z2**2 / 192.0 - z2**3 / 3072.0 + z2**4 / 46080.0
    return 0.5 * z * series


def phaseless_score_lambda(w, lam: float):
    """Intensity score with respect to lambda at fixed scale.

    ``w = t/(sigma2/2)`` is the intensity normalised by the chi-square scale.
    For lambda > 0:

        s_lambda = -1/2 + (1/2) sqrt(w/lambda)
                   * I1(sqrt(lambda w)) / I0(sqrt(lambda w)).

    For lambda < 1e-12 the lambda -> 0 limit s_lambda = -1/2 + w/4 is used
    (this is also the finite-difference limit of the ncx2 logpdf).  The
    Bessel ratio uses a small-argument series when the product ``lambda*w`` is
    tiny to avoid cancellation.
    """
    scalar_input = np.ndim(w) == 0
    w = np.atleast_1d(np.asarray(w, dtype=float))
    lam = float(lam)
    if lam < 1e-12:
        result = -0.5 + 0.25 * w
        return float(result[0]) if scalar_input else result
    out = np.full_like(w, np.nan)
    z = np.sqrt(lam * w)
    if w.size == 0:
        return out
    if np.any(w <= 0.0):
        out = out.copy()
        zero_mask = w <= 0.0
        # t = 0 belongs to the support; derivative d/dlambda log f_I(0) = -1/2.
        out[zero_mask] = -0.5
    pos = w > 0.0
    if np.any(pos):
        zp = z[pos]
        wp = w[pos]
        small = zp < 1e-6
        idx_pos = np.flatnonzero(pos)
        if np.any(small):
            # sqrt(w/lambda) * (I1/I0) ~ w/2 * (1 + O(z^2)) in this regime;
            # compute the product stably from the small-z series.
            sub = small
            ratio = _bessel_ratio_small(zp[sub])
            out[idx_pos[sub]] = -0.5 + 0.5 * np.sqrt(wp[sub] / lam) * ratio
        mid = ~small
    if np.any(mid):
        ratio = special.i1(zp[mid]) / special.i0(zp[mid])
        out[idx_pos[mid]] = -0.5 + 0.5 * np.sqrt(wp[mid] / lam) * ratio
    return float(out[0]) if scalar_input else out


def _phaseless_scale_score(t, sigma2, h=1e-6):
    """d/d(theta) log f_I(t) for the scale model via symmetric logpdf
    finite differences in the chi-square scale coordinate
    (scale = sigma2/2, theta = log sigma2 = log 2 + log scale)."""
    s = sigma2 / 2.0
    s_hi = s * (1.0 + h)
    s_lo = s * (1.0 - h)
    denom = float(np.log(s_hi) - np.log(s_lo))  # ~ 2h
    f_hi = stats.ncx2.logpdf(t, 2, nc=0.0, scale=s_hi)
    f_lo = stats.ncx2.logpdf(t, 2, nc=0.0, scale=s_lo)
    return (f_hi - f_lo) / denom


def phaseless_score(model: str, t, theta: float, params=None, scale_fd_h=1e-6):
    """Phaseless (intensity) score d/d(theta) log f_I(t) for a control model.

    For fixed-scale mean models the score is
    ``s_lambda(t/(sigma2/2), lambda) * d lambda/d theta``.  For the scale
    family (lambda = 0 for all theta) the score is the logpdf derivative with
    respect to theta computed by a central finite difference of
    ``scipy.stats.ncx2.logpdf`` at ``scale*(1 +/- 1e-6)``.
    """
    if model not in CONTROL_MODELS:
        raise ValueError(f"exact intensity likelihood not defined for model {model!r}")
    t = np.asarray(t, dtype=float)
    sigma2 = model_sigma2(model, theta, params=params)
    if model == "scale_family":
        if np.ndim(t) == 0:
            return float(_phaseless_scale_score(float(t), sigma2, h=scale_fd_h))
        return np.array([_phaseless_scale_score(float(v), sigma2, h=scale_fd_h) for v in t])
    lam = model_nc(model, theta, params=params)
    dlam = model_dnc_dtheta(model, theta, params=params)
    w = 2.0 * t / sigma2  # t/(sigma2/2)
    return dlam * phaseless_score_lambda(w, lam)


def _phaseless_integrand(model, theta, params, scale_fd_h, t):
    sigma2 = model_sigma2(model, theta, params=params)
    lam = model_nc(model, theta, params=params)
    pdf = stats.ncx2.pdf(t, 2, nc=lam, scale=sigma2 / 2.0)
    s = phaseless_score(model, t, theta, params=params, scale_fd_h=scale_fd_h)
    return s * s * pdf


def phaseless_fisher(model: str, theta: float, params=None, alpha=1e-12,
                     scale_fd_h=1e-6, epsabs=1e-12, epsrel=1e-10, limit=300):
    """Exact phaseless Fisher J_ph by quadrature of s_ph(t)^2 f_I(t).

    The integration domain is [ppf(alpha), ppf(1-alpha)] from the induced
    ncx2 intensity law.  ``phaseless_fisher_detail`` additionally computes the
    full-domain integral and reports the tail check.
    """
    return phaseless_fisher_detail(model, theta, params=params, alpha=alpha,
                                   scale_fd_h=scale_fd_h, epsabs=epsabs,
                                   epsrel=epsrel, limit=limit)["j_ph_bounded"]


def phaseless_fisher_detail(model: str, theta: float, params=None, alpha=1e-12,
                            scale_fd_h=1e-6, epsabs=1e-12, epsrel=1e-10,
                            limit=300) -> dict:
    sigma2 = model_sigma2(model, theta, params=params)
    lam = model_nc(model, theta, params=params)
    dist = stats.ncx2(df=2, nc=lam, scale=sigma2 / 2.0)
    lo = max(float(dist.ppf(alpha)), 0.0)
    hi = float(dist.ppf(1.0 - alpha))

    def integrand(t):
        return _phaseless_integrand(model, theta, params, scale_fd_h, t)

    j_b, err_b = integrate.quad(integrand, lo, hi, epsabs=epsabs,
                                epsrel=epsrel, limit=limit)
    # Full-domain cross check (upper tail to +inf and lower [0, lo]).
    j_hi, err_hi = integrate.quad(integrand, hi, np.inf, epsabs=epsabs,
                                  epsrel=epsrel, limit=limit)
    j_lo, err_lo = integrate.quad(integrand, 0.0, lo, epsabs=epsabs,
                                  epsrel=epsrel, limit=limit)
    j_full = j_b + j_hi + j_lo
    rel_diff = abs(j_full - j_b) / j_full if j_full > 0 else 0.0
    return {
        "j_ph_bounded": float(j_b),
        "j_ph_bounded_err": float(err_b),
        "j_ph_full": float(j_full),
        "j_ph_full_err": float(err_hi + err_lo),
        "rel_diff_full_vs_bounded": float(rel_diff),
        "lower_tail_mass": float(dist.cdf(lo)),
        "upper_tail_mass": float(1.0 - dist.cdf(hi)),
        "lambda": float(lam),
        "quad_lo": float(lo),
        "quad_hi": float(hi),
    }


# ---------------------------------------------------------------------------
# Analytic summary and Monte Carlo score checks
# ---------------------------------------------------------------------------


def analytic_values(model: str, theta: float, params=None) -> dict:
    """Analytic E1 quantities for a control / mismatch model."""
    p = _params(model, params)
    sigma2 = model_sigma2(model, theta, params=params)
    mu = model_mu(model, theta, params=params)
    dmu = model_dmu(model, theta, params=params)
    u = model_logsigma2_deriv(model, theta, params=params)
    out = {
        "model": model,
        "theta": float(theta),
        "sigma2": float(sigma2),
        "mu": complex(mu),
        "j_coh_raw_analytic": float(coherent_fisher(mu, dmu, sigma2, u)),
    }
    if model == "nuisance_phase":
        # 2x2 Fisher for (x, eta) with dmu/dx == dmu/deta; efficient target
        # information is J_xx - J_xeta J_etaeta^{-1} J_etax = 0.
        j = out["j_coh_raw_analytic"]
        out["j_eff_coh_analytic"] = 0.0
        out["j_coh_analytic"] = 0.0
        out["j_ph_analytic"] = 0.0
        out["note"] = "j_coh_analytic is the nuisance-eliminated efficient coherent Fisher; raw J_xx in j_coh_raw_analytic."
    elif model == "scale_family":
        out["j_coh_analytic"] = float(u * u)
        out["j_ph_analytic"] = float(u * u)  # intensity is sufficient for the scale score
    elif model == "mismatch_intensity_noise":
        out["j_coh_analytic"] = out["j_coh_raw_analytic"]
        m_prime = float(-2.0 * p["c"] * p["a"] * np.sin(theta))
        out["intensity_map"] = float(abs(mu) ** 2)
        out["d_intensity_dtheta"] = m_prime
        out["sigma_z2"] = float(p["sigma_z2"])
        # Mismatched direct-intensity Gaussian sensor: Z ~ N(|mu|^2, sigma_z2)
        # ignores the coherent speckle fluctuation of |Y|^2.  NOT an exact
        # induced-intensity likelihood, hence not covered by the theorem.
        out["j_mismatch_analytic"] = float(m_prime * m_prime / p["sigma_z2"])
    else:
        out["j_coh_analytic"] = out["j_coh_raw_analytic"]
        out["j_ph_analytic"] = phaseless_fisher(model, theta, params=params)
    return out


def _wls_slope_intercept_r2(x, y, w):
    """Weighted least-squares y ~ slope*x + intercept and R^2."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    wt = w.sum()
    xm = np.average(x, weights=w)
    ym = np.average(y, weights=w)
    sxx = float(np.sum(w * (x - xm) ** 2))
    sxy = float(np.sum(w * (x - xm) * (y - ym)))
    # Degenerate regressor: all bin means of the phaseless score are zero (or
    # at numerical noise level), e.g. models whose intensity is independent of
    # theta.  A slope is then meaningless.
    if sxx <= 1e-300 or np.all(np.abs(x - xm) < 1e-10):
        return float("nan"), float("nan"), float("nan")
    slope = sxy / sxx
    intercept = float(ym - slope * xm)
    yhat = slope * x + intercept
    ss_res = float(np.sum(w * (y - yhat) ** 2))
    ss_tot = float(np.sum(w * (y - ym) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-300 else float("nan")
    return slope, intercept, r2


def _binned_check(I, s_coh, s_ph, n_bins=40):
    """Equal-count binning estimate of s_ph(I) vs E[s_coh(Y)|I]."""
    n = len(I)
    edges = np.quantile(I, np.linspace(0.0, 1.0, n_bins + 1))
    # Make the last edge strictly larger than max(I).
    edges[-1] += 1e-12 * max(1.0, float(np.max(I)))
    idx = np.clip(np.searchsorted(edges, I, side="right") - 1, 0, n_bins - 1)
    centers = []
    counts = []
    means_coh = []
    means_ph = []
    var_coh = []
    for k in range(n_bins):
        mask = idx == k
        counts.append(int(mask.sum()))
        if counts[-1] == 0:
            continue
        centers.append(float(np.mean(I[mask])))
        means_coh.append(float(np.mean(s_coh[mask])))
        means_ph.append(float(np.mean(s_ph[mask])))
        var_coh.append(float(np.var(s_coh[mask], ddof=1)) if counts[-1] > 1 else float("nan"))
    centers = np.asarray(centers)
    means_coh = np.asarray(means_coh)
    means_ph = np.asarray(means_ph)
    counts = np.asarray(counts)
    var_coh = np.asarray(var_coh)
    diff = means_coh - means_ph
    max_abs_diff = float(np.max(np.abs(diff))) if diff.size else float("nan")
    wrmse = float(np.sqrt(np.sum(counts * diff**2) / np.sum(counts))) if counts.size else float("nan")
    slope, intercept, r2 = _wls_slope_intercept_r2(means_ph, means_coh, counts)
    se_bin = np.sqrt(var_coh / counts)
    std_diff = np.abs(diff) / se_bin if counts.size else np.array([])
    max_std_bin_diff = float(np.nanmax(std_diff)) if std_diff.size else float("nan")
    return {
        "n_bins_used": int(counts.size),
        "max_bin_mean_diff": max_abs_diff,
        "bin_weighted_rmse": wrmse,
        "max_std_bin_diff": max_std_bin_diff,
        "identity_within_4_per_bin_se": bool(max_std_bin_diff <= 4.0) if np.isfinite(max_std_bin_diff) else False,
        "regression_slope": slope,
        "regression_intercept": intercept,
        "regression_r2": r2,
        "bin_centers": centers,
        "bin_counts": counts,
        "bin_mean_s_coh": means_coh,
        "bin_mean_s_ph": means_ph,
    }


def _empirical_moments(v):
    v = np.asarray(v, dtype=float)
    n = v.size
    m = float(np.mean(v))
    se_m = float(np.std(v, ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    m2 = float(np.mean(v**2))
    se_m2 = float(np.std(v**2, ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    return m, se_m, m2, se_m2


def draw_complex_gaussian(rng, mu, sigma2, n):
    """n draws Y ~ CN(mu, sigma2), real/imag variance sigma2/2 each."""
    sigma = np.sqrt(sigma2)
    eps = (sigma / np.sqrt(2.0)) * (
        rng.standard_normal(n) + 1j * rng.standard_normal(n)
    )
    return mu + eps


def mc_score_check(model: str, theta: float, seed: int = 101, n_samples: int = 2000,
                   n_bins: int = 40, params=None, scale_fd_h=1e-6):
    """Monte Carlo score checks for one model / seed.

    Checks
    ------
    * s_ph(I) = E[s_coh(Y) | I] via 40 equal-count bins of I;
    * empirical Fisher means mean(s_coh^2), mean(s_ph^2) with MC standard
      errors and J_ph_emp <= J_coh_emp within ~3 MC standard errors;
    * zero score expectations mean(s_coh), mean(s_ph).

    For ``nuisance_phase`` the raw target coherent score and an identically
    equal nuisance score are computed; the nuisance-eliminated efficient
    coherent score s_eff = s_x - s_eta is zero and its empirical Fisher is
    reported as ``j_coh_emp`` (raw in ``j_coh_raw_emp``).  The score identity
    is still checked with the raw coherent score.
    """
    if model not in CONTROL_MODELS:
        raise ValueError(f"mc_score_check defined only for CONTROL_MODELS, got {model!r}")
    rng = np.random.default_rng(seed)
    sigma2 = model_sigma2(model, theta, params=params)
    mu = model_mu(model, theta, params=params)
    dmu = model_dmu(model, theta, params=params)
    u = model_logsigma2_deriv(model, theta, params=params)
    y = draw_complex_gaussian(rng, complex(mu), sigma2, n_samples)
    I = np.abs(y) ** 2
    s_coh = np.asarray(coherent_score(y, mu, dmu, sigma2, u))

    if model == "scale_family":
        s_ph = np.array([phaseless_score(model, t, theta, params=params,
                                         scale_fd_h=scale_fd_h) for t in I])
    else:
        s_ph = np.asarray(phaseless_score(model, I, theta, params=params))

    bin_res = _binned_check(I, s_coh, s_ph, n_bins=n_bins)
    m_coh, se_coh, f_coh, se_f_coh = _empirical_moments(s_coh)
    m_ph, se_ph, f_ph, se_f_ph = _empirical_moments(s_ph)
    # 4 sigma per-seed window (across 40 seed-model runs ~3sigma excursions are
    # expected by chance; aggregate means are also reported in the JSON).
    score_expectation_ok = bool(
        abs(m_coh) <= 4.0 * se_coh and abs(m_ph) <= 4.0 * se_ph
    )
    rec = {
        "seed": int(seed),
        "n_samples": int(n_samples),
        "j_ph_emp": f_ph,
        "j_ph_emp_se": se_f_ph,
        "j_coh_emp": f_coh,
        "j_coh_emp_se": se_f_coh,
        "mean_coh": m_coh,
        "mean_coh_se": se_coh,
        "mean_ph": m_ph,
        "mean_ph_se": se_ph,
        "score_expectation_ok": score_expectation_ok,
        "max_bin_mean_diff": bin_res["max_bin_mean_diff"],
        "bin_weighted_rmse": bin_res["bin_weighted_rmse"],
        "max_std_bin_diff": bin_res["max_std_bin_diff"],
        "identity_within_4_per_bin_se": bin_res["identity_within_4_per_bin_se"],
        "regression_slope": bin_res["regression_slope"],
        "regression_intercept": bin_res["regression_intercept"],
        "regression_r2": bin_res["regression_r2"],
    }
    if model == "nuisance_phase":
        dmu_n = model_dmu_nuisance(model, theta, params=params)
        s_n = np.asarray(coherent_score(y, mu, dmu_n, sigma2, u))
        s_eff = s_coh - s_n  # analytically zero by construction
        _, _, f_eff, se_f_eff = _empirical_moments(s_eff)
        rec["j_coh_raw_emp"] = f_coh
        rec["j_coh_raw_emp_se"] = se_f_coh
        rec["j_eff_coh_emp"] = f_eff
        rec["j_eff_coh_emp_se"] = se_f_eff
        rec["j_coh_emp"] = f_eff
        rec["j_coh_emp_se"] = se_f_eff
    # Contraction check: J_ph_emp <= J_coh_emp + 3*sqrt(se_ph^2 + se_coh^2).
    diff = rec["j_ph_emp"] - rec["j_coh_emp"]
    diff_se = float(np.sqrt(rec["j_ph_emp_se"] ** 2 + rec["j_coh_emp_se"] ** 2))
    # A tiny absolute slack keeps exactly-zero phaseless scores (phase-only,
    # nuisance-phase models) from failing on floating-point roundoff.
    slack = 1e-12 * max(1.0, abs(rec["j_ph_emp"]) + abs(rec["j_coh_emp"]))
    rec["mc_contraction_ok"] = bool(diff <= 3.0 * diff_se + slack)
    return rec
