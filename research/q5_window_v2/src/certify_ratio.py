"""Extra restricted proof: finite material inversion from calibrated E/M ratio.
Known geometry and two individually resolved modes; no transfer to class C.
"""
from certify_modal import *

def magnetic(e,x,bd):
    j,dj,y,dy=bd;s=e*x*x
    J=series(s,1,0);J1=series(s,1,1);D=series(s,2,0);D1=series(s,2,1)
    u=x*J;u1=x**3*J1;v=D;v1=x*x*D1
    n=u*dj-j*v;d=u*dy-y*v;n1=u1*dj-j*v1;d1=u1*dy-y*v1
    return n/d,(n1*d-n*d1)/(d*d)

def compact(v):
    den=10**30
    low=(v.lo*den).numerator//(v.lo*den).denominator
    high=-((-v.hi*den).numerator//(-v.hi*den).denominator)
    return I(F(low,den),F(high,den))

def cert():
    out=[]
    for x,lo,hi in [(F(1,5),F(3,2),F(4)),(F(3,20),F(2),F(5))]:
        bd=boundary(x);ml=None;R=F(0);count=int((hi-lo)*1000)
        for k in range(count):
            e=I(lo+F(k,1000),lo+F(k+1,1000));a,ap,_=q_enclosure(e,x,bd);b,bp=magnetic(e,x,bd)
            a,ap,b,bp=map(compact,(a,ap,b,bp))
            if b.lo<=0:raise AssertionError('magnetic positivity failed')
            u=b/a;up=(bp*a-b*ap)/(a*a)
            v=1+a*b;vp=ap*b+a*bp;w=1+b*b;wp=2*b*bp
            rp=(up*v+u*vp)/w-u*v*wp/(w*w)
            # |tM/tE| <= b/a * sqrt(1+a^2) <= b_hi/a_lo*(1+a_hi^2)
            R=max(R,b.hi/a.lo*(1+a.hi*a.hi))
            if rp.lo<=0:raise AssertionError(('ratio monotonicity failed',x,k))
            ml=rp.lo if ml is None else min(ml,rp.lo)
        out.append({'x':str(x),'domain':[str(lo),str(hi)],'boxes':count,
            'real_ratio_slope_lower':decimal_bound(ml),'complex_ratio_modulus_upper':decimal_bound(R,upper=True)})
    return {'arithmetic':'Fraction intervals; explicit outward rational compression to denominator 10^30','scope':'two resolved E/M modes, known geometry, shared complex gain',
            'results':out,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
    d=cert();p=Path(__file__).resolve().parents[1]/'results/ratio_certificate.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d,indent=2))
