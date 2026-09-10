"""Outward rational final budget exports; floating budget file remains diagnostic."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
from certify_modal import decimal_bound
HERE=Path(__file__).resolve().parents[1]

def certificate():
    c=json.loads((HERE/'results/modal_certificate.json').read_text())['results']
    rc=json.loads((HERE/'results/ratio_certificate.json').read_text())['results']
    # exp(2.629^2)>1000, established by a positive rational Taylor partial sum.
    gamma=F(2629,1000);s=gamma*gamma
    assert sum((s**n/F(factorial(n)) for n in range(41)),F(0))>1000
    sigmas=[F('0.000002666409'),F('0.000001804800')]
    slopes=[F(v['slope_lower']) for v in c];caps=[F(v['amplitude_upper']) for v in c]
    ba=gamma/F(4000)+F(1,2000)
    out=[];counts=[];legacy=[]
    for i,(m,H,sigma) in enumerate(zip(slopes,caps,sigmas)):
        bias=F(5,4)*(H/F(2000)+F(3,10000)*caps[1-i])+F(21,10**7)
        B=sigma*gamma/4+bias
        error=(B+H*ba)/((F(3,4)-ba)*m)
        assert error<F(1,10)
        out.append({'material':i+1,'total_reference_error_upper':decimal_bound(ba,upper=True),
          'modal_bias_upper':decimal_bound(bias,upper=True),
          'absolute_material_error_upper':decimal_bound(error,upper=True)})
        b_a=F('0.003290527')
        allowed=F(1,10)*(F(3,4)-b_a)*m-sigma*gamma-H*b_a
        legacy.append(decimal_bound(allowed))
        r=F(rc[i]['complex_ratio_modulus_upper']);mr=F(rc[i]['real_ratio_slope_lower'])
        qmin=F(c[i]['q_lower']);emin=qmin/(1+H*H);noise=sigma*gamma
        n=((noise*(1+r)+F(1,10)*mr*noise)/(F(1,10)*mr*F(3,4)*emin))**2
        nceil=-((-n.numerator)//n.denominator)
        assert F(nceil)>=n
        counts.append(nceil)
    return {'arithmetic':'exact rational; directed decimal export',
      'noise':'proper CN; 16 independent repeats/channel; union bound for three complex radial events',
      'event_probability_lower':'997/1000','gamma_squared_greater_than_log_1000_proof':'41 positive rational Taylor terms of exp(gamma^2)',
      'sigma_upper_bounds':[str(v) for v in sigmas],
      'conditional_reference_readout':out,
      'inherited_real_reference_allowed_modal_bias_lower':legacy,
      'EM_ratio_noise_only_sufficient_repeat_counts':counts,
      'hardware_attainability':'NOT established by this arithmetic certificate'}

if __name__=='__main__':
    p=HERE/'results/certified_budgets.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(certificate(),indent=2)+'\n');print(p.read_text())
