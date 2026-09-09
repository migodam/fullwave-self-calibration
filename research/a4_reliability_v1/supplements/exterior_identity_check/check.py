"""Post-frozen algebraic checks; no fitting, test-set reuse or threshold tuning."""
import json,sys
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from modal import dyad
z=s.symbols('z',positive=True,real=True);x=s.symbols('x',positive=True)
a=1+s.I/z-1/z**2; b=-1-3*s.I/z+3/z**2;l=a+b
W=s.factor(a*s.diff(l,z)-l*s.diff(a,z))
H=s.factor(s.expand_complex(2*a*s.conjugate(a)+l*s.conjugate(l)))
opt=s.factor(s.diff(x*(x+4)/(x*x+x+3)**2,x))
errs=[]
n=np.array([.3,-.5,.7]);n/=np.linalg.norm(n)
for zz in np.logspace(-3,3,301):
 y=dyad(n*zz/10,10);den=zz**4+zz**2+3
 for kind,expected in [('radiative',3*(zz**2+1)/den),('static',(zz**4+3*zz**2)/den)]:
  empirical=np.linalg.norm(y-dyad(n*zz/10,10,fidelity=kind))**2/np.linalg.norm(y)**2
  errs.append(abs(empirical-expected))
assert max(errs)<1e-12
out={'Wronskian':str(W),'squared_matrix_norm':str(H),'range_objective_derivative_x':str(opt),'exterior_identity_max_abs_error':max(errs),'z_values':301,'checks':602,'scope':'analytical model pairs, not a continuum DDA error bound'}
Path(__file__).with_name('result.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
