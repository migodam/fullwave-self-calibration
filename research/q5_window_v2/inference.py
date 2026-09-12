"""Correctly weighted inference utilities; no optimizer novelty is claimed."""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from scipy.special import spherical_jn, spherical_yn


def profile(y, f, sigma, reference=None, sigma_reference=.01, annulus=(.75, 1.25)):
    """Return sqrt(2)-realified residual and constrained complex GLS gain.

    sigma is the proper-complex standard deviation, E|noise|^2=sigma^2.
    A reference is an independent proper-complex observation z=g+noise.
    """
    y=np.asarray(y,complex).ravel(); f=np.asarray(f,complex).ravel()
    if y.shape != f.shape or not np.isfinite(sigma) or sigma<=0:
        raise ValueError('same-shaped data/field and a positive noise scale required')
    lo,hi=annulus
    if not (0<lo<=hi):raise ValueError('invalid gain annulus')
    v=f/sigma; b=y/sigma
    if reference is not None:
        if sigma_reference<=0:raise ValueError('invalid reference noise')
        v=np.r_[v,1/sigma_reference]; b=np.r_[b,reference/sigma_reference]
    den=np.vdot(v,v).real
    if den==0:raise ValueError('zero observation operator')
    c=np.vdot(v,b)/den
    gain=np.clip(abs(c),lo,hi)*np.exp(1j*np.angle(c))
    r=b-gain*v
    return np.sqrt(2)*np.r_[r.real,r.imag],complex(gain)


def modal_amplitude(epsilon, size):
    """SciPy point evaluation for estimation/tests, not an interval certificate."""
    if size<=0 or epsilon<=0:raise ValueError('positive material and size required')
    m=np.sqrt(epsilon); z=m*size
    j=spherical_jn(1,size); y=spherical_yn(1,size)
    dj=spherical_jn(1,size,True)+j/size
    dy=spherical_yn(1,size,True)+y/size
    ji=spherical_jn(1,z); di=spherical_jn(1,z,True)+ji/z
    q=(m*ji*dj-j*di)/(m*ji*dy-y*di)
    return abs(q/np.sqrt(1+q*q))


def amplitude_intervals(w,z,B,Ba,domains=((1.5,4.),(2.,5.)),sizes=(.2,.15),annulus=(.75,1.25)):
    """Exact-form set-membership reduction with numerical inverse evaluation.

    Outputs are numerical diagnostics, NOT outward-rounded machine certificates.
    Unknowns: common positive amplitude a and independent monotone modal epsilons.
    Errors satisfy |w_i-a h_i|<=B_i and |z-a|<=Ba, including deterministic bias.
    """
    w=np.asarray(w,float); B=np.asarray(B,float)
    if w.shape!=B.shape or len(w)!=len(domains) or np.any(B<0) or Ba<0:
        raise ValueError('invalid error bounds or dimensions')
    if not np.all(np.isfinite(np.r_[w,B,z,Ba])):raise ValueError('nonfinite input')
    hL=np.array([modal_amplitude(a,s) for (a,b),s in zip(domains,sizes)])
    hU=np.array([modal_amplitude(b,s) for (a,b),s in zip(domains,sizes)])
    alo=max(annulus[0],z-Ba,float(np.max((w-B)/hU)))
    ahi=min(annulus[1],z+Ba,float(np.min((w+B)/hL)))
    if not 0<alo<=ahi:return {'consistent':False,'status':'inconsistent_or_roundoff'}
    out=[]
    for i,((a,b),s) in enumerate(zip(domains,sizes)):
        def inv(v):
            if v<=hL[i]:return a
            if v>=hU[i]:return b
            return brentq(lambda e:modal_amplitude(e,s)-v,a,b,xtol=1e-13)
        left=inv((w[i]-B[i])/ahi);right=inv((w[i]+B[i])/alo)
        out.append([left,right])
    box=np.array(out)
    return {'consistent':True,'gain_interval':[alo,ahi],'material_intervals':out,
            'estimate':box.mean(axis=1).tolist(),'halfwidth':((box[:,1]-box[:,0])/2).tolist(),
            'status':'numerical_set_interval_diagnostic_not_machine_certificate'}
