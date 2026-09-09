"""Standard fixed-design task risk, not a new statistical theorem.

Explicit nuisance identifiability checks; no inverse weighted information
shortcut under true white noise. Action scores are surrogates unless error
coverage, linearization, and independence are separately justified.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import svd,solve


def weighted_task_risk(J,W,T,error_mean=None,error_factor=None):
    J=np.asarray(J,float); T=np.atleast_2d(T)
    s=svd(J,compute_uv=False)
    if s[-1]<=1e-10*s[0]:
        return {'risk':float('inf'),'variance':float('inf'),'bias2':float('inf'),'identified':False}
    L=solve(J.T@W@J,J.T@W,assume_a='sym')
    TL=T@L; variance=float(np.linalg.norm(TL)**2)
    bias2=0.
    if error_mean is not None: bias2+=float(np.linalg.norm(TL@error_mean)**2)
    if error_factor is not None:
        # Stochastic sample-model MSE, covariance EE^T; not worst-case ellipsoid.
        bias2+=float(np.linalg.norm(TL@error_factor)**2)
    return {'risk':variance+bias2,'variance':variance,'bias2':bias2,'identified':True,'L':L}


def precision_choices(error_samples):
    E=np.asarray(error_samples,float); mean=E.mean(axis=0); centered=E-mean
    cov=centered.T@centered/max(1,len(E)-1)
    vals,vecs=np.linalg.eigh(cov); j=np.argmax(vals)
    m=E.shape[1]; I=np.eye(m)
    return {'raw':I,'isotropic':I/(1+np.trace(cov)/m),
            'rank_one':np.linalg.inv(I+vals[j]*np.outer(vecs[:,j],vecs[:,j])),
            'sampled_eem':np.linalg.inv(I+cov)},mean,cov

@dataclass
class Evidence:
    geometry_risk: float
    tolerance2: float
    coverage_ok: bool
    fidelity_supported: bool
    electronics_ambiguity: bool=False
    reference_available: bool=False
    branch_ambiguity: bool=False
    acquisition_available: bool=False
    refinement_available: bool=False
    downweight_risk: float=float('inf')
    refined_risk: float=float('inf')


def action(e):
    if e.electronics_ambiguity:
        return 'add_electronics_reference' if e.reference_available else 'reject'
    if e.branch_ambiguity or not e.coverage_ok:
        return 'acquire_new_data' if e.acquisition_available else 'reject'
    if not e.fidelity_supported:
        return 'refine_model' if e.refinement_available else 'reject'
    if e.geometry_risk<=e.tolerance2: return 'accept'
    if e.refinement_available and e.refined_risk<min(e.geometry_risk,e.downweight_risk): return 'refine_model'
    if e.downweight_risk<e.geometry_risk and e.downweight_risk<=e.tolerance2: return 'downweight'
    return 'acquire_new_data' if e.acquisition_available else 'reject'
