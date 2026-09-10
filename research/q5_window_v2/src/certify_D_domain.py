"""Exact compact class-D domain and material-counterexample check."""
from fractions import Fraction as F
from pathlib import Path
import json
from certify_modal import I,boundary,q_enclosure,decimal_bound

def certificate():
    x=F(1,10);bd=boundary(x)
    q,qp,_=q_enclosure(I(2,F(23,10)),x,bd)
    assert q.lo>0 and qp.lo>0
    q0,_,_=q_enclosure(I(2),x,bd);q1,_,_=q_enclosure(I(F(23,10)),x,bd)
    g2=(q0*q0/(1+q0*q0))/(q1*q1/(1+q1*q1))
    assert F(3,4)**2<g2.lo<=g2.hi<F(5,4)**2
    return {'class':'D: epsilon in [2,23/10], ka=1/10, gain modulus [3/4,5/4]',
      'q_lower':decimal_bound(q.lo),'qprime_lower':decimal_bound(qp.lo),
      'counterexample_materials':['2','23/10'],'counterexample_gains':['1','t_E(2)/t_E(23/10)'],
      'second_gain_abs_squared':[decimal_bound(g2.lo),decimal_bound(g2.hi,upper=True)],
      'equality':'exact coefficient compensation at the same range, not a small optimized residual',
      'scope':'projected electric dipole only; not a class-C ambiguity'}

if __name__=='__main__':
    p=Path(__file__).resolve().parents[1]/'results/D_certificate.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(certificate(),indent=2)+'\n');print(p.read_text())
