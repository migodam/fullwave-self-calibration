"""Validated special functions and spherical quadrature; static core imports."""
from fractions import Fraction
from functools import lru_cache
import math
import numpy as np
from certified_modal import I
from interval_core import R,C,PI,stack,sincos,pi_integer


def j_regular(z,order):
    z=C.of(z);zz=z*z;lead=C(1);out=[]
    for l in range(order+1):
        if l:lead=lead*z/(2*l+1)
        term=lead;total=lead
        for s in range(30):
            term=-term*zz/(2*(s+1)*(2*l+2*s+3));total=total+term
        nxt=-term*zz/(2*31*(2*l+63))
        ratio=z.abs2()/(2*32*(2*l+65))
        if np.any(ratio.hi>=1):raise ValueError('Bessel remainder ratio >= 1')
        rem=nxt.abs()/(1-ratio)
        out.append(total+C(R(-rem.hi,rem.hi),R(-rem.hi,rem.hi)))
    return out


def h_outgoing(z,order):
    z=R.of(z);s,c=sincos(z)
    out=[C(s/z,-c/z)]
    if order:out.append(C(s/z.square()-c/z,-c/z.square()-s/z))
    for l in range(1,order):out.append(out[-1]*(2*l+1)/z-out[-2])
    return out


@lru_cache(None)
def norm_lm(l,m):
    return R.of((I.of(Fraction((2*l+1)*math.factorial(l-m),math.factorial(l+m)))/(4*pi_integer())).sqrt())


@lru_cache(None)
def gauss_certified(n):
    """n disjoint sign-changing rational brackets prove all n Legendre roots."""
    from numpy.polynomial.legendre import leggauss
    def leg(x):
        p0=I.of(1);p1=I.of(x)
        for k in range(1,n):p0,p1=p1,((2*k+1)*x*p1-k*p0)/(k+1)
        return p0,p1
    roots=[];weights=[];brackets=[]
    for guess in leggauss(n)[0]:
        a=Fraction(float(guess))-Fraction(1,10**12);b=a+Fraction(2,10**12)
        pa=leg(I.of(a))[1];pb=leg(I.of(b))[1]
        if not (pa.hi<0<pb.lo or pb.hi<0<pa.lo):raise ArithmeticError('Root bracket not verified')
        signa=1 if pa.lo>0 else -1
        for _ in range(145):
            mid=(a+b)/2;pm=leg(I.of(mid))[1]
            if pm.lo>0:sign=1
            elif pm.hi<0:sign=-1
            else:break
            if sign==signa:a=mid
            else:b=mid
            if b-a<Fraction(1,10**48):break
        if brackets and brackets[-1][1]>=a:raise ArithmeticError('Overlapping brackets')
        brackets.append((a,b));x=I.hull(a,b);p0,p1=leg(x)
        derivative=n*(p0-x*p1)/(1-x**2)
        w=2/((1-x**2)*derivative**2)
        if w.lo<=0:raise ArithmeticError('Weight is not positive')
        roots.append(R.of(x));weights.append(R.of(w))
    return stack(roots),stack(weights)


def spherical_quadrature(n):
    u,w=gauss_certified(n);phi=PI*R.of(np.arange(2*n))/n
    ss,cc=sincos(phi);s=(1-u.square()).sqrt()
    xx=s[:,None]*cc[None,:];yy=s[:,None]*ss[None,:]
    zz=R(np.broadcast_to(u.lo[:,None],xx.lo.shape),np.broadcast_to(u.hi[:,None],xx.hi.shape))
    pts=stack([xx,yy,zz],axis=-1)
    ww=w[:,None]*(PI/n);ww=R(np.broadcast_to(ww.lo,xx.lo.shape),np.broadcast_to(ww.hi,xx.hi.shape))
    return R(pts.lo.reshape(-1,3),pts.hi.reshape(-1,3)),R(ww.lo.ravel(),ww.hi.ravel())
