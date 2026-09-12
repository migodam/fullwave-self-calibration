"""Certified numerical inversion of the restricted monotone modal statistics.
This is an implementation of the proved scalar estimators, not a new GLS method.
Floating complex observations are treated as their exact binary rational values.
The returned interval covers the clipped statistic inverse, not automatically
true material: the measurement/model budget must also be applied.
"""
from fractions import Fraction as F
import json
from pathlib import Path
from certify_modal import I,q_enclosure,boundary
from certify_ratio import magnetic,compact
HERE=Path(__file__).resolve().parents[1]
DOMAINS=[(F(3,2),F(4)),(F(2),F(5))]
SIZES=[F(1,5),F(3,20)]

def norm2(z):return F(float(z.real))**2+F(float(z.imag))**2

def inverse_statistic(target,lower,upper,oracle,derivative_lower,tolerance=F(1,10**10)):
    """Exact interval bisection, including a derivative-safe ambiguous-sign step."""
    target=F(target);lower,upper=F(lower),F(upper);d=F(derivative_lower)
    if d<=0 or tolerance<=0:raise ValueError('Positive slope and tolerance required')
    for _ in range(256):
        if upper-lower<=tolerance:return I(lower,upper)
        mid=(lower+upper)/2;image=oracle(I(mid))
        if image.hi<target:lower=mid
        elif image.lo>target:upper=mid
        else:
            radius=max(target-image.lo,image.hi-target)/d
            new_lo=max(lower,mid-radius);new_hi=min(upper,mid+radius)
            if new_hi-new_lo>=upper-lower:
                raise ArithmeticError('Enclosure precision insufficient; do not return a false narrow bracket')
            lower,upper=new_lo,new_hi
    raise ArithmeticError('Bisection cap reached')

def amplitude_readout(Y,Z,mode):
    if mode not in (0,1):raise ValueError('Only the two certified class-M modes are supported')
    z2=norm2(Z)
    if z2<=0:raise ValueError('Reference magnitude is zero')
    c=json.loads((HERE/'results/modal_certificate.json').read_text())['results'][mode]
    m=F(c['slope_lower']);H=F(c['amplitude_upper']);hmin=F(c['q_lower'])/(1+H*H)
    x=SIZES[mode];bd=boundary(x)
    def square_amplitude(e):
        q,_,_=q_enclosure(e,x,bd);q=compact(q)
        return q*q/(1+q*q)
    return inverse_statistic(norm2(Y)/z2,*DOMAINS[mode],square_amplitude,2*hmin*m)

def ratio_readout(YE,YM,mode):
    if mode not in (0,1):raise ValueError('Only the two certified class-M modes are supported')
    y2=norm2(YE)
    if y2<=0:raise ValueError('Electric channel is zero')
    target=(F(float(YM.real))*F(float(YE.real))+F(float(YM.imag))*F(float(YE.imag)))/y2
    c=json.loads((HERE/'results/ratio_certificate.json').read_text())['results'][mode]
    x=SIZES[mode];bd=boundary(x)
    def ratio_real(e):
        q,_,_=q_enclosure(e,x,bd);b,_=magnetic(e,x,bd);q,b=compact(q),compact(b)
        return (b/q)*(1+q*b)/(1+b*b)
    return inverse_statistic(target,*DOMAINS[mode],ratio_real,F(c['real_ratio_slope_lower']))
