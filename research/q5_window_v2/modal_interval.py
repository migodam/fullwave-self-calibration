"""Continuum enclosure for the restricted lossless electric-dipole experiment.

No float inputs are used in the certificate. mpmath.iv supplies directed
rounding; all exported bounds are exact dyadic rationals (not nearest floats).
A positive-power series removes cancellation and the size-parameter sqrt.
This is a reproducible interval certificate, not a formally verified library.
"""
from __future__ import annotations
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import argparse, hashlib, json, time
from functools import lru_cache
import mpmath as mp
iv = mp.iv


def rat(q):
    q = Q(q)
    return iv.mpf(q.numerator) / q.denominator


def endpoint(t):
    sign, man, exp, bc = t
    return Q((-1 if sign else 1)*man)*(Q(2)**exp)


def bounds(x):
    return tuple(endpoint(t) for t in x._mpi_)


def enc(a, b):
    # Preserve enclosure when rational endpoints are rounded to iv dyadics.
    aa, bb = rat(a), rat(b)
    return iv.mpf([aa.a, bb.b])


def export(x):
    lo, hi = bounds(x)
    return [str(lo), str(hi)]


def falling(n, d):
    return factorial(n)//factorial(n-d) if n >= d else 0


def series(e, x, kind, derivative=0, xderivative=0, terms=18):
    """m*j1(mx) (kind J) or Dj(mx) (kind D), with uniform tail.

    For e<=5 and x<=1/5, successive absolute terms after index 18,
    including the supported derivatives, have ratio <1/2. Thus twice
    the first omitted absolute term bounds the entire remainder.
    """
    assert kind in ('J', 'D') and 0 <= derivative <= 2 and 0 <= xderivative <= 1
    elo,ehi=bounds(e); xlo,xhi=bounds(x)
    assert 0 < elo <= ehi <= 5 and 0 < xlo <= xhi <= Q(201,1000)
    out=rat(0)
    def term(k, ee, xx, absolute=False):
        ep=k+1 if kind=='J' else k
        xp=2*k+1 if kind=='J' else 2*k
        c=Q(2*(k+1),factorial(2*k+3))
        if kind=='D': c*=2*k+2
        c*=falling(ep,derivative)*falling(xp,xderivative)
        if not c: return rat(0)
        if not absolute and k%2: c=-c
        return rat(c)*ee**(ep-derivative)*xx**(xp-xderivative)
    # Horner evaluation of the same polynomial reduces repeated powers.
    k0=next(k for k in range(terms) if falling(k+1 if kind=='J' else k,derivative) and falling(2*k+1 if kind=='J' else 2*k,xderivative))
    coeff=[]
    for k in range(k0,terms):
        ep=k+1 if kind=='J' else k; xp=2*k+1 if kind=='J' else 2*k
        c=Q(2*(k+1),factorial(2*k+3))
        if kind=='D': c*=2*k+2
        c*=falling(ep,derivative)*falling(xp,xderivative)
        coeff.append(rat(-c if k%2 else c))
    v=coeff[-1]; z=e*x*x
    for c in reversed(coeff[:-1]): v=c+z*v
    out=v*e**((k0+1 if kind=='J' else k0)-derivative)*x**((2*k0+1 if kind=='J' else 2*k0)-xderivative)
    first=term(terms,rat(ehi),rat(xhi),True)
    # An explicit rational upper bound for all k >= terms. Each polynomial
    # derivative-factor ratio is <= ((k+3)/(k-1))**2, and the coefficient
    # ratio is <= 2*(k+2)/(k+1)/((2*k+4)*(2*k+5)).
    k=terms
    ratio=Q(2*(k+2),k+1)/((2*k+4)*(2*k+5))*Q(k+3,k-1)**2*5*xhi*xhi
    assert ratio < Q(1,2)
    radius=bounds(2*first)[1]
    return out+enc(-radius,radius)


@lru_cache(maxsize=8)
def outer(xlo,xhi,dps):
    x=enc(xlo,xhi)
    return tuple(series(rat(1),x,kind,xderivative=xd) for kind,xd in [('J',0),('D',0),('J',1),('D',1)])


def enclose(e, x):
    # j1(x) = J(1,x); the outer y1,Dy have stable signs at these x.
    j,dj,jx,djx=outer(*bounds(x),iv.dps)
    y=-iv.cos(x)/x**2-iv.sin(x)/x
    dy=iv.sin(x)/x**2+iv.cos(x)/x**3-iv.cos(x)/x
    J=series(e,x,'J'); D=series(e,x,'D')
    Je=series(e,x,'J',1); De=series(e,x,'D',1)
    N=J*dj-j*D; den=J*dy-y*D
    assert bounds(den)[0]>0
    q=N/den
    qp=((Je*dj-j*De)*den-N*(Je*dy-y*De))/den**2
    assert bounds(q)[0]>0 and bounds(qp)[0]>0
    h=q/iv.sqrt(1+q*q)
    hp=qp/(1+q*q)**rat(Q(3,2))
    # x derivative, used only for an explicit radius/frequency error budget.
    yx=iv.sin(x)/x**2+2*iv.cos(x)/x**3-iv.cos(x)/x+iv.sin(x)/x**2
    dyx=iv.sin(x)/x+2*iv.cos(x)/x**2-3*iv.sin(x)/x**3-3*iv.cos(x)/x**4
    Jx=series(e,x,'J',xderivative=1); Dx=series(e,x,'D',xderivative=1)
    Nx=Jx*dj+J*djx-jx*D-j*Dx
    denx=Jx*dy+J*dyx-yx*D-y*Dx
    qx=(Nx*den-N*denx)/den**2
    # |dt/dx|=|qx|/(1+q^2) for real q, t=iq/(1-iq).
    tx=abs(qx)/(1+q*q)
    return {'den':den,'q':q,'qp':qp,'h':h,'hp':hp,'abs_tx':tx}


def certify(outdir: Path, cells_per_unit=1000):
    iv.dps=60
    outdir.mkdir(parents=True,exist_ok=True)
    trace=outdir/'modal_cover.jsonl'
    if trace.exists(): raise FileExistsError(trace)
    result={'backend':'mpmath.iv','backend_version':mp.__version__,'dps':iv.dps,
            'bound_export':'exact signed dyadic rational strings',
            'endpoint_minimum_claim':False,'midpoint_Lipschitz_used':False,
            'scope':'class M, fixed exact rational size, no class C transfer','regions':[]}
    start=time.perf_counter()
    with trace.open('w') as f:
        for a,b,x in [(Q(3,2),Q(4),Q(1,5)),(Q(2),Q(5),Q(3,20))]:
            n=int((b-a)*cells_per_unit)
            extrema={};last=a
            for i in range(n):
                lo=a+(b-a)*Q(i,n);hi=a+(b-a)*Q(i+1,n)
                assert lo==last and hi>lo
                last=hi
                z=enclose(enc(lo,hi),rat(x))
                rec={'lo':str(lo),'hi':str(hi),'x':str(x),'bounds':{k:export(v) for k,v in z.items()}}
                f.write(json.dumps(rec,separators=(',',':'))+'\n')
                for k,v in z.items():
                    l,u=bounds(v)
                    if k not in extrema: extrema[k]=[l,u]
                    else: extrema[k]=[min(l,extrema[k][0]),max(u,extrema[k][1])]
            assert last==b
            result['regions'].append({'material':[str(a),str(b)],'size':str(x),'size_enclosure':export(rat(x)),
                'cells':n,'step':str((b-a)/n),'coverage':'adjacent closed rational intervals, exact first/last check',
                'bounds':{k:[str(v[0]),str(v[1])] for k,v in extrema.items()},
                'display_only':{k:[float(v[0]),float(v[1])] for k,v in extrema.items()}})
    result['trace_sha256']=hashlib.sha256(trace.read_bytes()).hexdigest()
    result['seconds']=time.perf_counter()-start
    target=outdir/'modal_certificate.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({**{k:v for k,v in result.items() if k!='regions'},'regions':[r['display_only'] for r in result['regions']]},indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path(__file__).parent/'results');p.add_argument('--cells-per-unit',type=int,default=1000)
    args=p.parse_args()
    if args.cells_per_unit<1: p.error('cells-per-unit must be positive')
    certify(args.out,args.cells_per_unit)
