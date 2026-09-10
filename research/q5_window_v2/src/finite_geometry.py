"""Finite projective chord of the exact electric-dipole tensor.
All bounds below use exact rational arithmetic. This is class D, not class C.
"""
from fractions import Fraction as F
import json
from pathlib import Path
from certify_modal import decimal_bound

def chord_squared(z,w):
    z,w=F(z),F(w)
    return 2*(z-w)**2*((z+w)**2+z*z*w*w)/((z**4+z*z+3)*(w**4+w*w+3))

def covered_chord_lower_squared(a,b,delta):
    """Covers all z,w in [a,b] with |z-w|>=delta. Conservative, no grid."""
    a,b,delta=map(F,(a,b,delta))
    if not 0<a<b or not 0<delta<=b-a:raise ValueError('invalid range domain')
    return 2*delta*delta*(4*a*a+a**4)/(b**4+b*b+3)**2

def main():
    h=F(1,20);tol=F(1,50);relative_error=F(1,200)
    out={'class':'D: isolated electric-dipole tensor; known direction; profiled scalar',
      'range_halfwidth_z':str(h),'target_error_z':str(tol),'total_error_over_min_signal_norm':str(relative_error),
      'formula':'2*(z-w)^2*((z+w)^2+z^2*w^2)/((z^4+z^2+3)*(w^4+w^2+3))',
      'noise_convention':'deterministic joint data-space ball; NOT fixed-power comparison',
      'not_certified_is_not_impossible':True,'designs':[]}
    for z in [F(1,5),F(2,5),F(3,5),F(4,5),F(1),F(3,2),F(2),F(3),F(5)]:
        lb=covered_chord_lower_squared(z-h,z+h,2*tol)
        out['designs'].append({'nominal_kR':str(z),'lower_squared_sine':decimal_bound(lb),
            'threshold_squared':str(4*relative_error**2),'finite_range_certificate':lb>4*relative_error**2})
    # Entire design interval, NOT only the listed design points:
    # z0 in [0.6,1.5] implies all true z in [0.55,1.55].
    lb=covered_chord_lower_squared(F(11,20),F(31,20),2*tol)
    out['whole_design_interval']={'nominal_interval':['3/5','3/2'],'lower_squared_sine':decimal_bound(lb),
             'passed':lb>4*relative_error**2}
    lbs=[covered_chord_lower_squared(F(60+i,100)-h,F(61+i,100)+h,2*tol) for i in range(90)]
    out['refined_whole_design_interval']={'nominal_interval':['3/5','3/2'],'design_boxes':90,'step':'1/100','lower_squared_sine':decimal_bound(min(lbs)),'passed':min(lbs)>4*relative_error**2,'coverage':'each interval bound covers every nominal design and every parameter pair, not node samples'}
    p=Path(__file__).resolve().parents[1]/'results/finite_geometry.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
