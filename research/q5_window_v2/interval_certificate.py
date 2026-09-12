"""Q5-V2: enclosing real Mie amplitudes with exact rational domains.

Run: python interval_certificate.py --output results/interval_certificate.json
No float input enters interval arithmetic. Every material box is evaluated,
not sampled. The trusted computing base is Python integers and mpmath.iv.
Exports use integer-directed decimal rounding of exact dyadic endpoints.
This certifies a restricted mathematical model, not laboratory error bounds.
"""
from __future__ import annotations
import argparse, hashlib, json, time, math
from fractions import Fraction
from pathlib import Path
import mpmath
from mpmath import iv

iv.dps = 45

def frac_raw(raw):
    sign, man, exp, _ = raw
    if man == 0:
        if exp != 0: raise ArithmeticError('nonfinite endpoint')
        return Fraction(0)
    return (-1 if sign else 1) * Fraction(man) * (Fraction(2) ** exp)

def bounds(x):
    return tuple(frac_raw(t) for t in x._mpi_)

def rat(p, q=1):
    return iv.mpf(p) / iv.mpf(q)

def interval(lo: Fraction, hi: Fraction):
    a=rat(lo.numerator,lo.denominator)
    b=rat(hi.numerator,hi.denominator)
    return iv.mpf([a.a,b.b])

def outward(f: Fraction, upper: bool, digits=24):
    scale=10**digits
    n=f.numerator*scale
    v=-((-n)//f.denominator) if upper else n//f.denominator
    s='-' if v<0 else ''
    v=abs(v)
    return f'{s}{v//scale}.{v%scale:0{digits}d}'

def export(x):
    lo,hi=bounds(x)
    ans=[outward(lo,False),outward(hi,True)]
    assert Fraction(ans[0])<=lo<=hi<=Fraction(ans[1])
    return ans

class Jet:
    """Value, first derivative, second derivative; no factorial scaling."""
    def __init__(self,v,d=0,dd=0):
        self.v=v if hasattr(v,'_mpi_') else rat(v)
        self.d=d if hasattr(d,'_mpi_') else rat(d)
        self.dd=dd if hasattr(dd,'_mpi_') else rat(dd)
    @staticmethod
    def asjet(x): return x if isinstance(x,Jet) else Jet(x)
    def __add__(self,b):
        b=self.asjet(b); return Jet(self.v+b.v,self.d+b.d,self.dd+b.dd)
    __radd__=__add__
    def __neg__(self): return Jet(-self.v,-self.d,-self.dd)
    def __sub__(self,b): return self+-self.asjet(b)
    def __rsub__(self,b): return self.asjet(b)+-self
    def __mul__(self,b):
        b=self.asjet(b)
        return Jet(self.v*b.v,self.d*b.v+self.v*b.d,self.dd*b.v+2*self.d*b.d+self.v*b.dd)
    __rmul__=__mul__
    def reciprocal(self):
        if bounds(self.v)[0]<=0<=bounds(self.v)[1]: raise ArithmeticError('zero denominator')
        return Jet(1/self.v,-self.d/self.v**2,2*self.d**2/self.v**3-self.dd/self.v**2)
    def __truediv__(self,b): return self*self.asjet(b).reciprocal()
    def __rtruediv__(self,b): return self.asjet(b)*self.reciprocal()
    def sqrt(self):
        v=iv.sqrt(self.v)
        return Jet(v,self.d/(2*v),self.dd/(2*v)-self.d**2/(4*v**3))
    def sin(self): return Jet(iv.sin(self.v),iv.cos(self.v)*self.d,iv.cos(self.v)*self.dd-iv.sin(self.v)*self.d**2)
    def cos(self): return Jet(iv.cos(self.v),-iv.sin(self.v)*self.d,-iv.sin(self.v)*self.dd-iv.cos(self.v)*self.d**2)
    def __pow__(self,n):
        if not isinstance(n,int): raise TypeError('integer powers only')
        if n<0: return self.reciprocal()**(-n)
        r=Jet(1)
        for _ in range(n): r=r*self
        return r

def j1(z): return z.sin()/z**2-z.cos()/z

def y1(z): return -z.cos()/z**2-z.sin()/z

def dj(z): return z.sin()/z+z.cos()/z**2-z.sin()/z**3

def dy(z): return -z.cos()/z+z.sin()/z**2+z.cos()/z**3

def material_series(eps,x,kind):
    """Enclose P=m*j1(mx) or Q=D_j(mx), including two eps derivatives.

    For n>=16, eps*x*x<=1/5 and r<=2, successive absolute
    differentiated terms have ratio <1/100. The omitted tail is
    bounded by (100/99) times its first term. See PROOFS.md.
    """
    assert bounds(eps*x*x)[1] <= Fraction(1,5)
    N=16; sums=[rat(0),rat(0),rat(0)]
    powers=[rat(1)]
    for _ in range(N+1): powers.append(powers[-1]*eps)
    for n in range(N):
        c=rat((-1)**n*(2*n+2),math.factorial(2*n+3))
        deg=n+1 if kind=='P' else n
        fac=x**(2*n+1) if kind=='P' else (2*n+2)*x**(2*n)
        for r in range(3):
            if deg>=r:
                sums[r]+=c*fac*math.prod(range(deg-r+1,deg+1))*powers[deg-r]
    total=Jet(*sums)
    n=N; coef=rat(2*n+2,math.factorial(2*n+3))
    deg=n+1 if kind=='P' else n
    fac=x**(2*n+1) if kind=='P' else (2*n+2)*x**(2*n)
    rs=[]
    for r in range(3):
        falling=math.prod(range(deg-r+1,deg+1))
        radius=coef*fac*falling*eps**(deg-r)*rat(100,99)
        upper=bounds(radius)[1]
        rs.append(interval(-upper,upper))
    return Jet(total.v+rs[0],total.d+rs[1],total.dd+rs[2])

def quantities(eps,x):
    xx=Jet(x)
    p=material_series(eps,x,'P'); s=material_series(eps,x,'Q')
    n=p*dj(xx)-j1(xx)*s
    d=p*dy(xx)-y1(xx)*s
    q=n/d; h=q/(1+q*q).sqrt()
    return {'den':d.v,'q':q.v,'qprime':q.d,'h':h.v,'hprime':h.d,'hsecond':h.dd}

def sweep(left,right,count,x):
    extrema={}; trace=hashlib.sha256(); old=None; start=time.perf_counter()
    for k in range(count):
        lo=left+(right-left)*Fraction(k,count)
        hi=left+(right-left)*Fraction(k+1,count)
        assert old is None or old==lo
        old=hi
        box=interval(lo,hi)
        assert bounds(box)[0]<=lo<=hi<=bounds(box)[1]
        r=quantities(box,x)
        for name,v in r.items():
            a,b=bounds(v)
            if name not in extrema: extrema[name]=[a,b]
            else: extrema[name]=[min(extrema[name][0],a),max(extrema[name][1],b)]
        trace.update(json.dumps([str(lo),str(hi),{n:export(v) for n,v in r.items()}],sort_keys=True).encode())
    assert old==right
    out={n:[outward(a,False),outward(b,True)] for n,(a,b) in extrema.items()}
    return {'domain':[str(left),str(right)],'boxes':count,'coverage':'consecutive exact rational boxes including endpoints',
            'x_enclosure':export(x),'enclosures':out,'all_box_sha256':trace.hexdigest(),
            'positive_q':extrema['q'][0]>0,'positive_hprime':extrema['hprime'][0]>0,
            'negative_hsecond':extrema['hsecond'][1]<0,'denominator_excludes_zero':not(extrema['den'][0]<=0<=extrema['den'][1]),
            'elapsed_seconds':time.perf_counter()-start}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('results/interval_certificate.json'))
    args=p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    rows=[sweep(Fraction(3,2),Fraction(4),2500,rat(1,5)),sweep(Fraction(2),Fraction(5),3000,rat(3,20))]
    # Exact finite secant boundary for delta=1/10 in amplitude-only experiment.
    # ba <= 0.003291 conservatively encloses the inherited reference noise event.
    # b_noise uses sqrt(log(1000)) directly in interval arithmetic.
    noise=[rat(2666408132,10**15),rat(1804798621,10**15)]
    delta=rat(1,10);amin=rat(3,4);ba=rat(3291,10**6)
    budgets=[]
    for x,upper,s in zip([rat(1,5),rat(3,20)],[4,5],noise):
        hu=quantities(rat(upper),x)['h'];hl=quantities(rat(upper)-2*delta,x)['h']
        critical=(amin*hu-(amin+2*ba)*hl)/2
        bn=s*iv.sqrt(iv.ln(rat(1000)))
        budgets.append({'h_upper':export(hu),'h_upper_minus_2delta':export(hl),'amplitude_noise_bound':export(bn),
                        'critical_total_modal_bound':export(critical),'sharp_remaining_model_budget':export(critical-bn)})
    result={'mpmath':mpmath.__version__,'iv_dps':iv.dps,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'separately observed lossless spheres; exact x; no geometry uncertainty or hardware certification',
            'domains':rows,'amplitude_only_sharp_budgets':budgets,
            'passed':all(r['positive_q'] and r['positive_hprime'] and r['negative_hsecond'] and r['denominator_excludes_zero'] for r in rows)}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
