"""Audit a sharp finite-secant assertion found in a pre-existing branch report.
Amplitude experiment only; NOT a necessity theorem for coherent complex data.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json,time
from certify_modal import I,series,boundary,q_enclosure,decimal_bound
from certify_ratio import compact
HERE=Path(__file__).resolve().parents[1]

def sqrt_bounds(v):
    """Outward sqrt bounds on a nonnegative rational interval."""
    if v.lo<0:raise ValueError('negative radicand')
    Q=10**30
    def floorroot(x):return F(isqrt((x*Q*Q).numerator//(x*Q*Q).denominator),Q)
    lo=floorroot(v.lo);hi=floorroot(v.hi)
    if hi*hi<v.hi:hi+=F(1,Q)
    assert lo*lo<=v.lo and hi*hi>=v.hi
    return I(lo,hi)

def h(e,x):
    q,_,_=q_enclosure(I(e),x,boundary(x));q=compact(q)
    return compact(q/sqrt_bounds(1+q*q))

def q_second(e,x,bd):
    j,dj,y,dy=bd;s=e*x*x
    J,J1,J2=[compact(series(s,1,d)) for d in range(3)]
    D,D1,D2=[compact(series(s,2,d)) for d in range(3)]
    u=e*x*J;up=x*(J+s*J1);upp=x**3*(2*J1+s*J2)
    v=D;vp=x*x*D1;vpp=x**4*D2
    n=u*dj-j*v;den=u*dy-y*v
    np=up*dj-j*vp;dp=up*dy-y*vp
    npp=upp*dj-j*vpp;dpp=upp*dy-y*vpp
    n,den,np,dp,npp,dpp=map(compact,(n,den,np,dp,npp,dpp))
    return compact((npp*den-n*dpp)/(den*den)-2*dp*(np*den-n*dp)/(den*den*den))

def main():
    start=time.perf_counter();out=[];delta=F(1,10);ba=F('0.003291');da=2*ba
    sigma=[F('0.000002666408132'),F('0.000001804798621')];gamma=F('2.629')
    for i,(x,lo,hi) in enumerate([(F(1,5),F(3,2),F(4)),(F(3,20),F(2),F(5))]):
        bd=boundary(x);upper=None
        for n in range(int((hi-lo)*1000)):
            box=I(lo+F(n,1000),lo+F(n+1,1000));qpp=q_second(box,x,bd)
            if qpp.hi>=0:raise AssertionError(('concavity unresolved',i,n))
            upper=qpp.hi if upper is None else max(upper,qpp.hi)
        hu=h(hi,x);hl=h(lo,x);hd=h(hi-2*delta,x)
        margin=compact(F(3,4)*hu-(F(3,4)+da)*hd)
        allowance=margin.lo/2-sigma[i]*gamma
        assert F(3,4)*hu.lo>(F(3,4)+da)*hl.hi
        out.append({'x':str(x),'domain':[str(lo),str(hi)],'boxes':int((hi-lo)*1000),
          'q_second_derivative_upper':decimal_bound(upper,upper=True),
          'h_concave':'q>0, qprime>0, qsecond<0 imply hsecond<0',
          'finite_secant_separation_lower':decimal_bound(margin.lo),
          'finite_secant_separation_upper':decimal_bound(margin.hi,upper=True),
          'noise_bound':decimal_bound(sigma[i]*gamma,upper=True),
          'additional_amplitude_bias_lower':decimal_bound(allowance),
          'other_channel_overlap_condition':'satisfied for each channel',
          'h_upper_endpoint_enclosure':[decimal_bound(hu.lo),decimal_bound(hu.hi,upper=True)]})
    # Fixed raw amplitude error experiment. It is distinct from normalized noise.
    raw=F(1,10**6)
    def c2(z):return (z**4-z*z+1)/(z**6)
    windows=[]
    for row in out:
        margin=F(row['finite_secant_separation_lower']);threshold=(2*raw/margin)**2
        brackets=[]
        for bound in ('finite_secant_separation_lower','finite_secant_separation_upper'):
            threshold_here=(2*raw/F(row[bound]))**2
            lo,hi=F(1,5),F(100)
            for _ in range(60):
                mid=(lo+hi)/2
                if c2(mid)>=threshold_here:lo=mid
                else:hi=mid
            brackets.append((lo,hi))
        lo,hi=brackets[0][0],brackets[1][1]
        windows.append({'size_parameter':row['x'],'threshold_c_squared':str(threshold),
            'certified_nominal_kR_upper_lower':decimal_bound(lo),
            'transition_upper':decimal_bound(hi,upper=True),
            'note':'both separation endpoints propagated; true amplitude-experiment cutoff enclosed; no necessity claim for coherent fields'})
    result={'arithmetic':'Fraction intervals, outward sqrt, exact closed-box coverage',
            'class':'M-amplitude: bounded additive amplitude observations and amplitude reference',
            'reference_total_error':str(ba),'target_error':str(delta),'results':out,
            'fixed_raw_error_window':{'raw_absolute_amplitude_error':str(raw),'common_exterior_lower_kR':'1/5','results':windows},
            'seconds':time.perf_counter()-start}
    p=HERE/'results/secant_certificate_final.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(result,indent=2)+'\n');print(p.read_text())
if __name__=='__main__':main()
