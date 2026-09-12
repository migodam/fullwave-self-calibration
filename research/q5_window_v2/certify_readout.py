"""Exact rational readout tail/budget and a finite calibration counterexample.

The certificate concerns ideal radial field samples, registered lossless class M,
and the explicit instrument-tolerance model in MANUSCRIPT_SECTIONS.md. It is NOT
proof that an actual probe array attains those tolerances.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json
from dyadic_interval import I,SCALE,modal
HERE=Path(__file__).resolve().parent

def odd_factorial(n):
    out=1
    for j in range(1,n+1,2):out*=j
    return out

def upper_tail(x,first):
    """Sum_{l>=first} 3*l*(l+1)*x**(2*l+1)/(2*l-1)!!, geometric majorant."""
    if not (0<x<=F(1,4)) or first<2:raise ValueError('Unsupported tail domain')
    first_term=3*first*(first+1)*x**(2*first+1)/odd_factorial(2*first-1)
    ratio=F(first+2,first)*x*x/(2*first+1)
    return first_term/(1-ratio)

def interval(f):
    f=F(f);return I.rational(f.numerator,f.denominator)

def lower_exp(z,terms=30):
    return sum((z**n/F(factorial(n)) for n in range(terms)),F(0))

def certify(out):
    if out.exists():raise FileExistsError(out)
    # All constants are chosen outward from the two stored full-cover certificates.
    nominal=json.loads((HERE/'results/modal_certificate.json').read_text())
    sizecert=json.loads((HERE/'results/modal_size_100ppm.json').read_text())
    constants=[dict(x=F(1,5),m=F(439,10**6),H=F(271,10**5),H_size=F(3034,10**6),Lx=F(49597,10**6),sigma=F(26664081,10**13)),
               dict(x=F(3,20),m=F(133,10**6),H=F(1305,10**6),H_size=F(1546,10**6),Lx=F(34613,10**6),sigma=F(18047986,10**13))]
    def dyadic(v,key):return F(int(v[key]),1<<v['denominator_power_of_two'])
    for c,n,s in zip(constants,nominal['rows'],sizecert['rows']):
        assert c['m']<=dyadic(n['bounds']['hprime'],'lo_dyadic')
        assert c['H']>=dyadic(n['bounds']['h'],'hi_dyadic')
        assert c['H_size']>=dyadic(s['bounds']['h'],'hi_dyadic')
        assert c['Lx']>=dyadic(s['bounds']['hx'],'hi_dyadic')
    normal_t=F(7,2);complex_t=F(263,100)
    # Exact exp series lower bounds and Mills' inequality prove all three
    # scalar tail probabilities < 0.001; sqrt(2*pi)>5/2 is conservative.
    pc=1/lower_exp(complex_t**2)
    pa=2/(normal_t*F(5,2)*lower_exp(normal_t**2/2))
    assert pc<F(1,1000) and pa<F(1,1000)
    rows=[]
    for reps,sqrt_reps in [(1,1),(16,4)]:
        ref_noise=normal_t*F(1,1000)/sqrt_reps
        ref_drift=F(1,1000);ba=ref_noise+ref_drift
        for i,c in enumerate(constants):
            x=c['x'];xmax=x*F(10001,10000)
            tail=upper_tail(xmax,15);all_high=upper_tail(xmax,2)
            calibration=F(5,4)*F(1,1000)*F(1001,1000)*(c['H_size']+all_high)
            illumination=F(5,4)*F(1,1000)*c['H_size']
            sizing=F(5,4)*c['Lx']*x/F(10000)
            other=F(1,2*10**6)
            bdet=calibration+illumination+sizing+other+F(5,4)*F(1001,1000)*tail
            bn=complex_t*c['sigma']/sqrt_reps
            bound=(bn+bdet+c['H']*ba)/((F(3,4)-ba)*c['m'])
            budget=F(1,10)*(F(3,4)-ba)*c['m']-bn-c['H']*ba
            quantities={'alias_tail':tail,'all_higher_modes_majorant':all_high,'node_calibration':calibration,
              'incident_modal_calibration':illumination,'size_bias':sizing,'other_residual_cap':other,
              'total_deterministic_amplitude_bias':bdet,'modal_noise_event':bn,'reference_noise_event':ref_noise,
              'reference_drift':ref_drift,'reference_total_event':ba,'allowed_deterministic_bias':budget,
              'material_error_bound':bound}
            rows.append({'sphere':i+1,'repetitions':reps,'pass_delta_0p1':bound<=F(1,10),
                         'quantities':{k:interval(v).record() for k,v in quantities.items()}})
    qa=modal(I.rational(47,10),I.rational(3,20))['q']
    qb=modal(I.rational(491,100),I.rational(3,20))['q']
    eta=(qb-qa)/((2*qa*qb)**2+(qa+qb)**2).sqrt()
    assert eta.lo>0 and eta.hi<=I.rational(125,10000).lo
    result={'kind':'rational_tail_budget_and_finite_modal_worlds','status':'proved_with_declared_readout_assumptions',
      'hardware_attainment':'unverified','readout':{'ntheta':8,'nphi':4,'kR':1,'cancelled_degrees':list(range(2,15)),
        'samples_two_spheres_with_16_repeats':1024,'background_subtracted_acquisitions_if_doubled':2048,
        'amplitude_reference_repetitions':16},
      'noise_tail_bounds':{'proper_complex_one_channel':interval(pc).record(),'real_reference':interval(pa).record(),
                          'union_upper':interval(2*pc+pa).record()},
      'budget_rows':rows,'competing_worlds':{'epsilon1_both':3,'epsilon2_A_rational':[47,10],
        'epsilon2_B_rational':[491,100],'common_gain_both':1,'reference_both':1,
        'relative_second_channel_uncertainty_cap_rational':[125,10000],
        'eta_required':eta.record(),'identity':'kappa_A=2*t_B/(t_A+t_B); kappa_B=2*t_A/(t_A+t_B)',
        'scope':'class M with one unreferenced modal transfer factor, not simultaneous lossy class C'},
      'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'dyadic_interval.py']}}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':certify(HERE/'results/readout_certificate.json')
