"""Outward set-membership estimator for class M, conditional on supplied bounds.

Exact decimal strings/Fractions only. It does NOT verify hardware/noise bounds.
The derivative floor 1/10000 is certified on both fixed material domains by
interval_certificate.py. No claim of a new set-inversion algorithm is made.
"""
from fractions import Fraction as F
from interval_certificate import rat,quantities,bounds,outward

DOMAINS=((F(3,2),F(4)),(F(2),F(5)))
SIZES=(F(1,5),F(3,20))


def exact(x):
    if isinstance(x,float):raise TypeError('use an exact decimal string or Fraction')
    return F(x)


def amplitude(e,x):
    return quantities(rat(e.numerator,e.denominator),rat(x.numerator,x.denominator))['h']


def inverse_enclosure(target,domain,size,tol=F(1,10**12)):
    lo,hi=domain;target=F(target)
    # Globally proved, deliberately conservative derivative floor.
    slope=F(1,10000)
    while hi-lo>tol:
        mid=(lo+hi)/2;hl,hu=bounds(amplitude(mid,size))
        if target<hl:hi=mid
        elif target>hu:lo=mid
        else:
            radius=max(target-hl,hu-target)/slope
            lo=max(lo,mid-radius);hi=min(hi,mid+radius)
            break
    return lo,hi


def estimate(w,z,B,Ba,delta=F(1,10)):
    if len(w)!=2 or len(B)!=2:raise ValueError('exactly two class-M channels required')
    w=list(map(exact,w));B=list(map(exact,B));z=exact(z);Ba=exact(Ba);delta=exact(delta)
    if any(x<0 for x in B) or Ba<0 or delta<0:raise ValueError('nonnegative bounds required')
    aL=max(F(3,4),z-Ba);aU=min(F(5,4),z+Ba)
    for wi,bi,dom,x in zip(w,B,DOMAINS,SIZES):
        hL=amplitude(dom[0],x);hU=amplitude(dom[1],x)
        aL=max(aL,bounds(rat((wi-bi).numerator,(wi-bi).denominator)/hU)[0])
        aU=min(aU,bounds(rat((wi+bi).numerator,(wi+bi).denominator)/hL)[1])
    if aL>aU:return {'status':'certified_inconsistent_conditional_on_error_bounds'}
    boxes=[]
    for wi,bi,dom,x in zip(w,B,DOMAINS,SIZES):
        left=inverse_enclosure((wi-bi)/aU,dom,x)[0]
        right=inverse_enclosure((wi+bi)/aL,dom,x)[1]
        if left>right:return {'status':'certified_inconsistent_conditional_on_error_bounds'}
        boxes.append((left,right))
    passed=all((hi-lo)/2<=delta for lo,hi in boxes)
    return {'status':'certified_task_error_conditional_on_bounds' if passed else 'unresolved_task_tolerance',
            'gain_outer_interval':[str(aL),str(aU)],
            'material_outer_intervals':[[str(lo),str(hi)] for lo,hi in boxes],
            'midpoint_estimates_exact':[str((lo+hi)/2) for lo,hi in boxes],
            'halfwidth_upper_decimal':[outward((hi-lo)/2,True) for lo,hi in boxes],
            'scope':'class M only; assumes certified complete input error bounds, known geometry and modal readout'}
