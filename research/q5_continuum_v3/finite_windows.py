"""Finite-set task windows and exact algebra for bounded amplitude readouts.

The deterministic theorem is about compatible observation sets, not Fisher
information. Numerical inverse functions below must be replaced by a validated
monotone inverse when their returned endpoints are used as certificates.
"""
from __future__ import annotations
import numpy as np


def profile_gain(y, f, sigma, reference=None, reference_sigma=None,
                 annulus=(.75,1.25)):
    """Exact annulus GLS minimizer for a SINGLE common complex gain.

    Proper complex noise convention: E|noise|^2=sigma^2. This minimizes the
    complex GLS objective, not a marginal likelihood; there is no logdet.
    """
    y,f=np.asarray(y,complex).ravel(),np.asarray(f,complex).ravel()
    if y.shape!=f.shape or sigma<=0:raise ValueError('Invalid field data or noise')
    numerator=np.vdot(f,y)/sigma**2;denominator=np.vdot(f,f).real/sigma**2
    if reference is not None:
        if reference_sigma is None or reference_sigma<=0:raise ValueError('Reference variance required')
        numerator+=reference/reference_sigma**2;denominator+=1/reference_sigma**2
    if not (0<annulus[0]<=annulus[1]):raise ValueError('Invalid gain annulus')
    if denominator<=0:raise ValueError('Zero experiment vector')
    z=numerator/denominator
    phase=z/abs(z) if z!=0 else 1+0j
    return phase*np.clip(abs(z),*annulus)


def real_whiten(residual,sigma):
    z=np.asarray(residual,complex).ravel()
    if sigma<=0:raise ValueError('Positive complex noise std required')
    return np.sqrt(2)*np.r_[z.real,z.imag]/sigma


def finite_tube_task_window(means,tasks,error_radii,tolerance):
    """EXACT finite-list algebra evaluated in floating point, not continuum proof.

    Observation tubes are closed Euclidean balls, radius_i. Task loss is l_inf.
    Return an incompatible-task pair whose tubes intersect, or finite-list pass.
    A finite-list pass does not establish a continuum parameter window.
    """
    mu=np.asarray(means);q=np.asarray(tasks,float)
    if q.ndim==1:q=q[:,None]
    r=np.broadcast_to(np.asarray(error_radii,float),(len(mu),))
    if len(q)!=len(mu) or tolerance<0 or np.any(r<0):raise ValueError('Invalid finite experiment')
    for i in range(len(mu)):
        for j in range(i):
            if np.max(abs(q[i]-q[j]))>2*tolerance and np.linalg.norm(mu[i]-mu[j])<=r[i]+r[j]:
                return {'finite_list_pass':False,'witness_indices':[j,i],'continuum_claim':False}
    return {'finite_list_pass':True,'witness_indices':None,'continuum_claim':False}


def magnitude_feasible_intervals(readouts,reference,channel_bounds,reference_bound,
                                 h_lower,h_upper,inverses,annulus=(.75,1.25)):
    """Coordinate projections of the exact rectangular magnitude feasible set.

    Readout class: R_i=a*h_i(e_i)+u_i, Z=a+v, |u_i|<=B_i, |v|<=B_a.
    h_i are positive strictly increasing. Bounds returned by arbitrary floating
    inverses are diagnostic; the algebra is exact but does not certify an inverse.
    """
    R,B,lo,hi=[np.asarray(t,float) for t in (readouts,channel_bounds,h_lower,h_upper)]
    if not (R.shape==B.shape==lo.shape==hi.shape) or np.any(B<0) or np.any(lo<=0) or np.any(hi<lo):
        raise ValueError('Invalid amplitude model')
    if reference_bound<0:raise ValueError('Negative reference bound')
    amin,amax=annulus
    left=max(amin,reference-reference_bound,np.max((R-B)/hi))
    right=min(amax,reference+reference_bound,np.min((R+B)/lo))
    if left>right:return None
    lower=[inv(np.clip((R[i]-B[i])/right,lo[i],hi[i])) for i,inv in enumerate(inverses)]
    upper=[inv(np.clip((R[i]+B[i])/left,lo[i],hi[i])) for i,inv in enumerate(inverses)]
    return {'amplitude_interval':[float(left),float(right)],'material_lower':lower,
            'material_upper':upper,'material_midrange':((np.array(lower)+upper)/2).tolist(),
            'coordinate_halfwidths':((np.array(upper)-lower)/2).tolist(),
            'validated_inverse_required_for_certificate':True}
