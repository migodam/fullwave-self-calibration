"""Covered finite secants and radial windows, ideal separate amplitude readout.

This independently checks the amplitude-window formula appearing in another
branch report. Negative results here do NOT apply to full complex data.
The radial experiment fixes RAW absolute error; the 18-direction certificate
at kR=2 is a separate implementation, not silently extended to every radius.
"""
from fractions import Fraction as F
from pathlib import Path
import math,json,hashlib
from interval_certificate import I,SCALE,modal_intervals
HERE=Path(__file__).resolve().parent


def root(iv):
    if iv.lo<0:raise ValueError('negative square root')
    lo=math.isqrt(iv.lo*SCALE);hi=math.isqrt(iv.hi*SCALE)
    if hi*hi<iv.hi*SCALE:hi+=1
    return I(lo,hi)


def amplitude(e,x):
    q=modal_intervals(I.of(e),x)[0]
    return q/root(1+q*q)


def radial_squared(z):
    z=F(z)
    return z**-2-z**-4+z**-6


def window_bracket(budget,raw=F(1,10**6)):
    """Enclose unique positive c(z)=raw/budget root by exact rational bisection."""
    lo,hi=F(1,5),F(100)
    target=(I.of(raw)/budget)**2
    assert radial_squared(lo)>F(target.hi,SCALE)
    assert radial_squared(hi)<F(target.lo,SCALE)
    for _ in range(100):
        if hi-lo<F(1,10**12):break
        mid=(lo+hi)/2;v=radial_squared(mid)
        if v>F(target.hi,SCALE):lo=mid
        elif v<F(target.lo,SCALE):hi=mid
        else:raise ArithmeticError('rounding overlap; increase lattice precision')
    return [str(lo),str(hi)],I.of(lo,hi).export()


def main():
    out=HERE/'results/finite_windows.json'
    if out.exists():raise FileExistsError(out)
    amin=F(3,4);amax=F(5,4);Ba=F('0.003291');delta=F(1,10)
    da=min(2*Ba,amax-amin)
    rows=[]
    for i,(lower,upper,x) in enumerate([(F(3,2),F(4),F(1,5)),(F(2),F(5),F(3,20))]):
        hU=amplitude(upper,x);hV=amplitude(upper-2*delta,x);hL=amplitude(lower,x)
        # Verify the other-coordinate overlap construction even for the full gain annulus.
        assert hU.lo>0 and (hU/hL).lo>I.of(amax/amin).hi
        budget=(amin*hU-(amin+da)*hV)/2
        assert budget.lo>0
        fraction,decimal=window_bracket(budget)
        noise=I.of([F('2.666408132e-6'),F('1.804798621e-6')][i])*root(I.of(F('6.908')))
        bias=budget-noise
        rows.append({'object':i+1,'max_total_normalized_amplitude_error':budget.export(),
                     'radial_endpoint_exact_enclosure':fraction,'radial_endpoint_decimal_enclosure':decimal,
                     'single_read_model_budget_safe_noise_radius':bias.export()})
    # exp(6.908)>1000 certifies the complex 0.001 tail event without uncertified log.
    assert sum((F('6.908')**j/F(math.factorial(j)) for j in range(80)),F(0))>1000
    # Real reference: a rational central-normal-probability enclosure.
    z=F('3.291');acc=I.of(0)
    for j in range(100):
        acc+=F((-1)**j)*z**(2*j+1)/F(2**j*math.factorial(j)*(2*j+1))
    next_term=z**201/F(2**100*math.factorial(100)*201)
    acc+=I.of(-next_term,next_term)
    central=root(2/I.of(F(333,106),F(355,113)))*acc
    real_tail=1-central
    assert real_tail.hi<I.of(F(1,1000)).lo
    result={'reference_real_gaussian_tail_upper':real_tail.export()[1],
            'scope':'ideal separated amplitude experiment, not full-complex impossibility and not class C',
            'reference_bound_exact':str(Ba),'raw_projected_error_exact':'1/1000000',
            'common_standoff_domain':'kR > 1/5; ideal complete modal projection',
            'derivative_identity':'(c^2)\'= -2*((z^2-1)^2+2)/z^7 < 0',
            'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'provenance':'posthoc independent verification of finite-amplitude formula seen in branch report d5e4c2eb; no new algorithm claim'}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
