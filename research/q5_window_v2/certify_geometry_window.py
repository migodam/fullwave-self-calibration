"""Finite radial-calibration window for the exact dipole tensor class G.

This is NOT a material certificate and NOT a two-lossy-sphere certificate.
Covers continuous nominal-radius bins and all finite separated radius pairs.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
from dyadic_interval import I,SCALE
from certify_readout import lower_exp,interval
HERE=Path(__file__).resolve().parent

def orientation_polynomial(z,w):
    # |<T(z,n),T(w,n)>|^2 - |<T(z,n),T(w,m)>|^2 for n perpendicular m.
    return (3*w**4*z**4-7*w**4*z**2+3*w**4+24*w**3*z**3-7*w**2*z**4
            +27*w**2*z**2+9*w**2+3*z**4+9*z**2+27)

def certify(out):
    if out.exists():raise FileExistsError(out)
    q=F(600);rho=F(5);gap=F(1,25);wpose=F(1,10);nsub=30;rows=[]
    for i in range(3,40):
        nominal=(F(i,10),F(i+1,10));lo=nominal[0]-wpose;hi=nominal[1]+wpose
        edges=[lo+(hi-lo)*F(j,nsub) for j in range(nsub+1)]
        minimum=None;count=0;bad_orientation=0
        for j in range(nsub):
            for k in range(j,nsub):
                a,b,c,d=edges[j],edges[j+1],edges[k],edges[k+1]
                if d-a<gap:continue
                z=I.box(a.numerator,a.denominator,b.numerator,b.denominator)
                w=I.box(c.numerator,c.denominator,d.numerator,d.denominator)
                p=orientation_polynomial(z,w)
                if p.lo<=0:bad_orientation+=1;lower=I.rational(0)
                else:
                    # Domain restriction |z-w|>=gap is used explicitly, not
                    # replaced by a zero lower bound on overlapping boxes.
                    sep=max(gap,c-b)
                    num=interval(2*sep*sep*((a+c)**2+(a*c)**2)).sqrt()
                    den=interval((b**4+b*b+3)*(d**4+d*d+3)).sqrt()
                    lower=num/den
                minimum=lower.lo if minimum is None else min(minimum,lower.lo);count+=1
        threshold=interval(2*rho/q)
        rows.append({'nominal_z_bin':[[v.numerator,v.denominator] for v in nominal],
            'covered_actual_z_domain':[[v.numerator,v.denominator] for v in [lo,hi]],
            'parameter_cells':nsub,'covered_separated_pairs':count,
            'orientation_unproved_boxes':bad_orientation,
            'projective_sine_lower':I(minimum,minimum).record(),
            'finite_radius_error_certified':minimum>threshold.hi})
    # ||n||^2/sigma^2 is Gamma(shape=9, scale=1). Exact positive exp-series lower
    # bound yields a rigorous upper bound on its integer-shape tail at 25.
    from math import factorial
    tail=sum((F(25)**j/F(factorial(j)) for j in range(9)),F(0))/lower_exp(F(25),80)
    assert tail<F(3,1000)
    result={'status':'finite_radial_geometry_certificate_only','class':'G_exact_regular_electric_dipole_tensor',
      'assumptions':{'proper_complex_entries':9,'sigma_each_entry':1,'received_signal_norm_lower':600,
        'noise_norm_event_upper':5,'radial_pose_halfwidth_in_kR':.1,'radial_error_tolerance_in_kR':.02,
        'material_and_common_gain_profiled':True,'orientation_may_vary':True},
      'finite_separation_formula':'sin^2=2*(z-w)^2*((z+w)^2+z^2*w^2)/((z^4+z^2+3)*(w^4+w^2+3))',
      'noise_event_failure_upper':interval(tail).record(),'rows':rows,
      'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'dyadic_interval.py']},
      'not_claimed':['No material recovery','No universal optimal standoff','No fixed-power optimum',
                     'Uncertified bins are not proved impossible','No transfer to class C']}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print('certified nominal bins',[r['nominal_z_bin'] for r in rows if r['finite_radius_error_certified']])
    print('all orientation boxes verified',all(r['orientation_unproved_boxes']==0 for r in rows))
if __name__=='__main__':certify(HERE/'results/geometry_window_certificate.json')
