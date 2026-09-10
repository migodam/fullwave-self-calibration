"""Registered Q5 V2 diagnostics, NOT a final population test or certificate."""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k]='1'
from pathlib import Path
import sys, json, hashlib, time, platform, resource
import numpy as np
from scipy.optimize import least_squares
import scipy
from multipole import Cluster, directions
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
import maxwell3d as mx
CENTERS=np.array([[-.06,0,0],[.055,.02,0]])
RADII=[.035,.025]
LOSS=np.array([.03,.05])
STARTS=([2.,3.,0.],[1.6,2.2,-1.],[3.8,4.8,1.])


def encode(z):
    a=np.asarray(z)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def profile(y, f, sigma, reference=None, sigma_r=.01):
    """Exact scalar annulus projection, proper-complex real whitening."""
    if sigma<=0 or sigma_r<=0:
        raise ValueError('Standard deviations must be positive')
    v=np.asarray(f,complex)/sigma
    b=np.asarray(y,complex)/sigma
    if reference is not None:
        v=np.r_[v,1/sigma_r+0j];b=np.r_[b,reference/sigma_r]
    energy=float(np.vdot(v,v).real)
    if energy<=0:raise ValueError('Zero model and no gain reference')
    g0=np.vdot(v,b)/energy
    g=np.clip(abs(g0),.75,1.25)*np.exp(1j*np.angle(g0))
    r=b-g*v
    return np.sqrt(2)*np.r_[r.real,r.imag],g


def fit_model(model,rx,y,sigma,reference=None,extra=None):
    extra_points,extra_index=None,None
    if extra is not None:extra_points,extra_index,extra_y=extra
    def predicted(theta):
        f=model.forward(theta,rx)
        if extra is not None:
            # Same shift and same gain as the base acquisition.
            f=np.r_[f,model.forward(theta,extra_points)[extra_index]]
        return f
    yy=y if extra is None else np.r_[y,extra_y]
    fits=[];begin=time.perf_counter()
    for initial in STARTS:
        counter=[0]
        def residual(t):
            counter[0]+=1
            return profile(yy,predicted(t),sigma,reference)[0]
        f=least_squares(residual,initial,bounds=([1.5,2,-2],[4,5,2]),
                        max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8)
        fits.append({'theta':f.x.tolist(),'rss_real_whitened':float(f.fun@f.fun),
                     'nfev':int(f.nfev),'actual_residual_evaluations':counter[0],
                     'optimizer_success':bool(f.success),'status':int(f.status)})
    best=min(fits,key=lambda r:r['rss_real_whitened'])
    theta=np.array(best['theta']);res,g=profile(yy,predicted(theta),sigma,reference)
    return {'theta':theta.tolist(),'gain':encode(g),'starts':fits,
            'rss_real_whitened':float(res@res),'seconds':time.perf_counter()-begin}


