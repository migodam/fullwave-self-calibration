"""Implement the class-M amplitude estimator; no class-C guarantee is implied."""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from scipy.special import spherical_jn, spherical_yn
from multipole import angular, quadrature, sphere_t

DOMAINS=((1.5,4.0,.2),(2.0,5.0,.15))


def readout_operator():
    """18 directions and a Cartesian row acting only on tangential components."""
    directions,weights=quadrature(3,6)
    _,p,_=angular(directions,1)
    x=2.0
    dh=spherical_jn(1,x,True)+spherical_jn(1,x)/x+1j*(
        spherical_yn(1,x,True)+spherical_yn(1,x)/x)
    c=-np.sqrt(6*np.pi)
    row=-weights[:,None]*p[:,:,1].conj()/(c*dh)
    return directions,row


def project_modal_field(field):
    """field shape (...,18,3), Cartesian representation of measured tangent field."""
    field=np.asarray(field,complex)
    if field.shape[-2:]!=(18,3) or not np.all(np.isfinite(field)):
        raise ValueError('finite (...,18,3) fields required')
    return np.sum(field*readout_operator()[1],axis=(-2,-1))


def material_from_amplitude(y,reference,object_index):
    """Invert a separately read modal amplitude with a positive real reference.

    This numerical root finder implements the certified estimator. Its floating
    point root tolerance is an implementation tolerance, not interval arithmetic.
    """
    if object_index not in (0,1):
        raise ValueError('object_index must be 0 or 1')
    if not np.isreal(reference) or not np.isfinite(reference) or reference<=0:
        raise ValueError('positive finite real amplitude reference required')
    if not np.isfinite(y):raise ValueError('nonfinite modal observation')
    lo,hi,x=DOMAINS[object_index]
    def h(e):return float(abs(sphere_t(e,x,1)[4]))
    v=float(np.clip(abs(y)/float(np.real(reference)),h(lo),h(hi)))
    if v<=h(lo):return lo
    if v>=h(hi):return hi
    return float(brentq(lambda e:h(e)-v,lo,hi,xtol=5e-13))


def require_continuum_certificate(metadata):
    """Fail closed before interpreting a class-C diagnostic as certification."""
    required=('continuous_forward_enclosure','roundoff_enclosure','complete_pair_cover')
    if not all(metadata.get(k) is True for k in required):
        raise ValueError('UNRESOLVED: class-C continuous enclosure or complete cover missing')
    # Metadata is a gate, never itself proof that an enclosure is valid.
    return True
