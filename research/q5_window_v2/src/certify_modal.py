"""Rational interval certificate: no floating-point arithmetic in proof kernel.
The two size parameters are 1/5 and 3/20. All arithmetic uses Fraction.
Polynomial tails are bounded analytically, not by sample differences.
"""
from fractions import Fraction as F
from math import factorial
from dataclasses import dataclass
import argparse,hashlib,json,time
from pathlib import Path

@dataclass(frozen=True)
class I:
    lo:F
    hi:F
    def __init__(self,lo,hi=None):
        a=F(lo);b=a if hi is None else F(hi)
        if a>b:raise ValueError('reversed interval')
        object.__setattr__(self,'lo',a);object.__setattr__(self,'hi',b)
    def __add__(self,b):
        b=as_i(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-as_i(b)
    def __rsub__(self,b):return as_i(b)+-self
    def __mul__(self,b):
        b=as_i(b);a=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi];return I(min(a),max(a))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=as_i(b)
        if b.lo<=0<=b.hi:raise ValueError('denominator contains zero')
        return self*I(1/b.hi,1/b.lo)
    def __rtruediv__(self,b):return as_i(b)/self

def as_i(a):return a if isinstance(a,I) else I(a)
def falling(n,d):return factorial(n)//factorial(n-d)

def series(s,r,d,N=16):
    """d-th derivative of sum (-1)^n (2n+2)^r s^n/(2n+3)!.
    r in {1,2}, d<=2, 0<=s<=1. From n>=16 successive absolute
    derivative terms have ratio <=1/2. The tail is <=2*first omitted.
    """
    assert r in (1,2) and 0<=d<=2 and 0<=s.lo<=s.hi<=1 and N>=4
    coeff=[F((-1)**n*(2*n+2)**r*falling(n,d),factorial(2*n+3)) for n in range(d,N)]
    p=I(0)
    for c in coeff[::-1]:p=p*s+c
    rem=2*F((2*N+2)**r*falling(N,d),factorial(2*N+3))*s.hi**(N-d)
    return p+I(-rem,rem)

def trig(x,N=20):
    assert 0<=x<=1
    sn=sum((F((-1)**n,factorial(2*n+1))*x**(2*n+1) for n in range(N)),F(0))
    cs=sum((F((-1)**n,factorial(2*n))*x**(2*n) for n in range(N)),F(0))
    rs=x**(2*N+1)/factorial(2*N+1);rc=x**(2*N)/factorial(2*N)
    return I(sn-rs,sn+rs),I(cs-rc,cs+rc)

def boundary(x):
    s,c=trig(x)
    j=I(x)*series(I(x*x),1,0); dj=series(I(x*x),2,0)
    y=-c/(x*x)-s/x;dy=s/(x*x)+c/(x**3)-c/x
    return j,dj,y,dy

def q_enclosure(e,x,bd):
    j,dj,y,dy=bd;s=e*(x*x)
    J=series(s,1,0);J1=series(s,1,1);D=series(s,2,0);D1=series(s,2,1)
    u=e*x*J;u1=x*(J+s*J1);v=D;v1=x*x*D1
    n=u*dj-j*v; den=u*dy-y*v
    n1=u1*dj-j*v1;den1=u1*dy-y*v1
    q=n/den; q1=(n1*den-n*den1)/(den*den)
    return q,q1,den

def decimal_bound(v,places=16,upper=False):
    """Directed decimal export by integer division, no float conversion."""
    t=v*10**places;n=t.numerator//t.denominator
    if upper and F(n)<t:n+=1
    sign='-' if n<0 else '';n=abs(n)
    return sign+str(n//10**places)+'.'+str(n%10**places).zfill(places)

def certify(stepden=1000):
    ans=[];start=time.perf_counter()
    for xmin,lo,hi in [(F(1,5),F(3,2),F(4)),(F(3,20),F(2),F(5))]:
        count=int((hi-lo)*stepden);assert F(count,stepden)==hi-lo
        bd=boundary(xmin);ml=None;H=F(0);qlo=None;dlo=None
        for n in range(count):
            left=lo+F(n,stepden);right=lo+F(n+1,stepden)
            q,qp,den=q_enclosure(I(left,right),xmin,bd)
            if q.lo<=0 or qp.lo<=0 or den.lo<=0:raise AssertionError(('positivity failure',n))
            # h'=q'/(1+q^2)^(3/2) >= q'_lo/(1+q_hi^2)^2.
            slope=qp.lo/(1+q.hi*q.hi)**2
            ml=slope if ml is None else min(ml,slope)
            H=max(H,q.hi);qlo=q.lo if qlo is None else min(qlo,q.lo)
            dlo=den.lo if dlo is None else min(dlo,den.lo)
        ans.append(dict(x=str(xmin),domain=[str(lo),str(hi)],boxes=count,step=str(F(1,stepden)),
                    slope_lower=decimal_bound(ml),amplitude_upper=decimal_bound(H,upper=True),
                    q_lower=decimal_bound(qlo),denominator_lower=decimal_bound(dlo),
                    first_left=str(lo),last_right=str(right),coverage='exact rational consecutive closed boxes'))
        print(ans[-1],flush=True)
    return dict(arithmetic='fractions.Fraction; exact rational interval kernel',series_terms=16,trig_terms=20,
                endpoint_minimum_claim=False,midpoint_lipschitz_used=False,results=ans,
                seconds=time.perf_counter()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stepden',type=int,default=1000);ap.add_argument('--output',type=Path)
    ar=ap.parse_args();out=certify(ar.stepden)
    if ar.output:
        if ar.output.exists():raise FileExistsError(ar.output)
        ar.output.write_text(json.dumps(out,indent=2)+'\n')
