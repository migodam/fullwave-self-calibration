"""Model-conditional passive-medium state/derivative certificates (parent derivation)."""
from dataclasses import dataclass
import numpy as np
from scipy.special import j1

@dataclass(frozen=True)
class PassiveBound:
    available: bool
    alpha: float
    delta: float
    weighted_inverse_bound: float
    euclidean_inverse_bound: float
    sqrt_u: np.ndarray
    reason: str

def passive_bound(u, tau, k, h, margin=1e-12):
    u = np.asarray(u, dtype=float)
    if u.ndim != 1 or not np.all(np.isfinite(u)) or np.min(u) <= 0 or tau <= 0:
        return PassiveBound(False, float('nan'), float('nan'), np.inf, np.inf,
                            np.sqrt(np.maximum(u, 0)), 'positive_common_loss_required')
    a = h/np.sqrt(np.pi)
    delta = float(k*k*h*h/4 - np.pi*k*a*j1(k*a)/2)
    alpha = float(tau/(1+tau*tau)-max(delta,0)*np.max(u)-margin)
    if alpha <= 0:
        return PassiveBound(False, alpha, delta, np.inf, np.inf, np.sqrt(u),
                            'passivity_margin_unavailable')
    bound = 1/(np.sqrt(1+tau*tau)*alpha)
    return PassiveBound(True, alpha, delta, bound,
                        bound*np.sqrt(np.max(u)/np.min(u)), np.sqrt(u), 'analytic_bound')

def state_bound(cert, residual, S, sigma=1.0):
    """Bound whitened complex output error, per state RHS; not sqrt2-real error."""
    if not cert.available:
        return np.inf
    r = np.asarray(residual)
    rw = r/cert.sqrt_u if r.ndim == 1 else r/cert.sqrt_u[:,None]
    norm_r = np.linalg.norm(rw, axis=0)
    sensing = np.linalg.norm(np.asarray(S)*cert.sqrt_u[None,:], 'fro')/sigma
    return sensing*cert.weighted_inverse_bound*norm_r

def derivative_bound(cert, state_residual, derivative_residual, Mv, S, Sv=None, sigma=1):
    if not cert.available:
        return np.inf
    weighted_state_error = cert.weighted_inverse_bound*np.linalg.norm(state_residual/cert.sqrt_u)
    weighted_Mv = np.asarray(Mv)*cert.sqrt_u[None,:]/cert.sqrt_u[:,None]
    weighted_deriv_error = cert.weighted_inverse_bound*(
        np.linalg.norm(derivative_residual/cert.sqrt_u)
        +np.linalg.norm(weighted_Mv,'fro')*weighted_state_error)
    answer = np.linalg.norm(S*cert.sqrt_u[None,:],'fro')*weighted_deriv_error
    if Sv is not None:
        answer += np.linalg.norm(Sv*cert.sqrt_u[None,:],'fro')*weighted_state_error
    return float(answer/sigma)

def objective_interval(residual, error_bound, prior=0.0):
    r = np.linalg.norm(residual)
    return (0.5*max(0.0,r-error_bound)**2+prior,
            0.5*(r+error_bound)**2+prior)

def gradient_error_bound(residual, jacobian, residual_error, jacobian_error):
    return (jacobian_error*np.linalg.norm(residual)
            +(np.linalg.norm(jacobian,'fro')+jacobian_error)*residual_error)
