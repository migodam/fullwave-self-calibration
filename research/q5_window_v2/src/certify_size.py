"""Uniform size derivative certificate; optional class-M manufacturing budget.
Uses the exact rational kernel. A size band is covered, never sampled.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
from certify_modal import I,series,decimal_bound
from certify_ratio import compact

def sincos(x):
    s=x*x
    def horner(offset):
        p=I(0)
        for n in reversed(range(20)):p=p*s+F((-1)**n,factorial(2*n+offset))
        if offset==1:p=p*x
        tail=x.hi**(40+offset)/factorial(40+offset)
        return compact(p+I(-tail,tail))
    return horner(1),horner(0)

def derivative(e,x):
    s=e*x*x; sx=x*x
    J,Jp=series(s,1,0),series(s,1,1)
    v=series(s,2,0);vx=2*e*x*series(s,2,1)
    u=e*x*J;ux=e*(J+2*s*Jp)
    j=x*series(sx,1,0);dj=series(sx,2,0)
    si,co=sincos(x)
    y=-co/(x*x)-si/x;dy=si/(x*x)+co/(x*x*x)-co/x
    jx=dj-j/x;djx=-jx/x-j+j/(x*x)
    yx=dy-y/x;dyx=-yx/x-y+y/(x*x)
    n=u*dj-j*v;d=u*dy-y*v
    nx=ux*dj+u*djx-jx*v-j*vx
    dx=ux*dy+u*dyx-yx*v-y*vx
    assert d.lo>0
    qp=(nx*d-n*dx)/(d*d)
    return max(abs(qp.lo),abs(qp.hi))

def main():
    out=[]
    # Cover eps by width 1/100; size interval itself is included in each box.
    for x0,lo,hi in [(F(1,5),F(3,2),F(4)),(F(3,20),F(2),F(5))]:
        dx=x0*F(1,20000);x=I(x0-dx,x0+dx);L=F(0)
        for i in range(int((hi-lo)*100)):
            L=max(L,derivative(I(lo+F(i,100),lo+F(i+1,100)),x))
        out.append({'x0':str(x0),'relative_size_uncertainty':'1/20000',
            'epsilon_boxes':int((hi-lo)*100),'size_interval':[str(x.lo),str(x.hi)],
            'amplitude_size_derivative_upper':decimal_bound(L,upper=True),
            'gain_weighted_amplitude_bias_upper':decimal_bound(F(5,4)*L*dx,upper=True),
            'fits_1e_6_allocation':F(5,4)*L*dx < F(1,10**6)})
    p=Path(__file__).resolve().parents[1]/'results/size_certificate.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps({'scope':'class M lossless modes only; size ka uncertain, all other shape/constitutive assumptions exact','results':out},indent=2)+'\n')
    print(p.read_text())
if __name__=='__main__':main()
