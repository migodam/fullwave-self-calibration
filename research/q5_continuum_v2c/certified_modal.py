"""Integer fixed-point interval certificate for the restricted lossless class M.

All interval endpoints are integers / 10**60. No floating-point or external
interval library enters a proof calculation. Trigonometric polynomials use a
Taylor remainder, sqrt uses integer isqrt, and export is outward rounded.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import factorial, isqrt
import argparse, hashlib, json, time
from pathlib import Path

S=10**60

def ceildiv(a,b):
    return -((-a)//b)

@dataclass(frozen=True)
class I:
    lo: int
    hi: int
    def __post_init__(self):
        if self.lo>self.hi:raise ValueError('Reversed interval')
    @staticmethod
    def of(v):
        if isinstance(v,I):return v
        if isinstance(v,float):raise TypeError('Float inputs are forbidden')
        v=Fraction(v);return I(v.numerator*S//v.denominator,
                               ceildiv(v.numerator*S,v.denominator))
    @staticmethod
    def hull(a,b):
        return I(I.of(a).lo,I.of(b).hi)
    def __add__(self,b):
        b=I.of(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-I.of(b)
    def __rsub__(self,b):return I.of(b)+-self
    def __mul__(self,b):
        b=I.of(b);v=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return I(min(v)//S,ceildiv(max(v),S))
    __rmul__=__mul__
    def inv(self):
        if self.lo<=0<=self.hi:raise ZeroDivisionError('Interval contains zero')
        return I(S*S//self.hi,ceildiv(S*S,self.lo))
    def __truediv__(self,b):return self*I.of(b).inv()
    def __rtruediv__(self,b):return I.of(b)*self.inv()
    def __pow__(self,n):
        if not isinstance(n,int) or n<0:raise ValueError('Nonnegative integer power required')
        if n==0:return I.of(1)
        if n==1:return self
        if n==2:
            a,b=self.lo,self.hi
            low=0 if a<=0<=b else min(a*a,b*b)
            return I(low//S,ceildiv(max(a*a,b*b),S))
        v=I.of(1);base=self
        while n:
            if n&1:v=v*base
            n//=2
            if n:base=base**2
        return v
    def sqrt(self):
        if self.lo<0:raise ValueError('Negative sqrt')
        a=isqrt(self.lo*S);b=isqrt(self.hi*S)
        return I(a,b if b*b==self.hi*S else b+1)
    def absmax(self):return I(max(abs(self.lo),abs(self.hi)),max(abs(self.lo),abs(self.hi)))
    def export(self,digits=24):
        d=10**(60-digits)
        def text(v):
            sign='-' if v<0 else '';v=abs(v)
            return sign+str(v//10**digits)+'.'+str(v%10**digits).zfill(digits)
        return [text(self.lo//d),text(ceildiv(self.hi,d))]


def trig(z,cosine=False):
    z=I.of(z)
    if z.absmax().hi>S:raise ValueError('Taylor implementation domain is |z| <= 1')
    n=18;shift=0 if cosine else 1
    p=I.of(Fraction((-1)**n,factorial(2*n+shift)))
    zz=z**2
    for j in range(n-1,-1,-1):
        p=p*zz+Fraction((-1)**j,factorial(2*j+shift))
    if not cosine:p=p*z
    power=2*n+2 if cosine else 2*n+3
    rem=z.absmax()**power/factorial(power)
    return p+I(-rem.hi,rem.hi)


def radial(z):
    """j1, j1', D_j, D_j', D_j'', y1, D_y; identities audited in Q5."""
    z=I.of(z);s=trig(z);c=trig(z,True)
    j=s/z**2-c/z
    jp=s/z+2*c/z**2-2*s/z**3
    jpp=-2*jp/z-(1-2/z**2)*j
    dj=jp+j/z
    djp=-jp/z-j+j/z**2
    djpp=-jpp/z-jp+2*jp/z**2-2*j/z**3
    y=-c/z**2-s/z
    dy=-c/z+s/z**2+c/z**3
    return j,jp,dj,djp,djpp,y,dy


def modal(eps,x):
    eps=I.of(eps);x=I.of(x);m=eps.sqrt();z=m*x
    j,jp,dj,djp,_,y,dy=radial(x)
    ji,jip,di,dip,_,_,_=radial(z)
    n=m*ji*dj-j*di;d=m*ji*dy-y*di
    mp=1/(2*m);zp=x*mp
    ne=(mp*ji+m*jip*zp)*dj-j*dip*zp
    de=(mp*ji+m*jip*zp)*dy-y*dip*zp
    q=n/d;qe=(ne*d-n*de)/d**2
    h=q/(1+q**2).sqrt()
    he=qe/((1+q**2)*(1+q**2).sqrt())
    return {'q':q,'denominator':d,'q_prime':qe,'h':h,'h_prime':he}


def certify(step_den=1000, export_boxes=None):
    if not isinstance(step_den,int) or step_den<=0 or step_den%2:
        raise ValueError('step_den must be a positive even integer')
    output={'arithmetic':'integer endpoints / 10^60; sqrt=isqrt; sin/cos Taylor degree 37/36 with rigorous remainder',
            'input_policy':'Exact rational endpoints and ka; binary floats rejected',
            'slope_infimum_at_endpoint_claimed':False,'classes':[]}
    file=Path(export_boxes).open('w') if export_boxes else None
    try:
        for x,lo,hi in [(Fraction(1,5),Fraction(3,2),Fraction(4)),
                        (Fraction(3,20),Fraction(2),Fraction(5))]:
            first=int(lo*step_den);last=int(hi*step_den)
            lower={};upper={};hashes=hashlib.sha256();prev=lo
            for index in range(first,last):
                a,b=Fraction(index,step_den),Fraction(index+1,step_den)
                assert a==prev and b>a;prev=b
                v=modal(I.hull(a,b),I.of(x))
                for key in ['denominator','q','q_prime','h_prime']:
                    if v[key].lo<=0:raise ArithmeticError(f'Cannot certify {key} on {a},{b}')
                for key,val in v.items():
                    lower[key]=min(lower.get(key,val.lo),val.lo)
                    upper[key]=max(upper.get(key,val.hi),val.hi)
                row={'ka':str(x),'a':str(a),'b':str(b),'bounds':{k:t.export() for k,t in v.items()}}
                line=json.dumps(row,sort_keys=True,separators=(',',':'))+'\n'
                hashes.update(line.encode())
                if file:file.write(line)
            assert prev==hi
            # Since q>0 and h'>0 on the whole interval, h's maximum is at hi.
            # This does NOT assert that h' attains its minimum at hi.
            endpoint=modal(I.of(hi),I.of(x))['h']
            result={'ka_exact':str(x),'epsilon_exact':[str(lo),str(hi)],
                    'coverage':{'first_index':first,'last_index_exclusive':last,
                                'denominator':step_den,'boxes':last-first,'gaps':0},
                    'bounds':{k:I(lower[k],upper[k]).export() for k in lower},
                    'amplitude_max_endpoint_enclosure':endpoint.export(),
                    'per_box_sha256':hashes.hexdigest()}
            output['classes'].append(result)
    finally:
        if file:file.close()
    output['certificate_pass']=True
    return output


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default='results/modal_certificate.json')
    p.add_argument('--boxes');p.add_argument('--step-den',type=int,default=1000)
    a=p.parse_args();out=Path(a.output)
    if out.exists():raise FileExistsError(f'Refuse to overwrite {out}')
    start=time.perf_counter();result=certify(a.step_den,a.boxes)
    result['wall_seconds']=time.perf_counter()-start
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
