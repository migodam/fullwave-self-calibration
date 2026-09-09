"""Small independent algebra checks for the information-versus-bias distinction."""
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(8460);checks=[]
for case in range(10):
    nl=rng.normal(size=(20,5));nh=rng.normal(size=(12,5))
    bl=rng.normal(size=(20,3));bh=rng.normal(size=(12,3))
    profile=lambda b,n:b.T@(b-n@np.linalg.lstsq(n,b,rcond=None)[0])
    il=profile(bl,nl);ij=profile(np.r_[bl,bh],np.r_[nl,nh])
    h=nl.T@nl;q=bh-nh@np.linalg.solve(h,nl.T@bl)
    w=np.linalg.inv(np.eye(12)+nh@np.linalg.solve(h,nh.T))
    err=np.linalg.norm(ij-il-q.T@w@q)/np.linalg.norm(ij)
    assert err<1e-12 and np.linalg.eigvalsh(ij-il).min()>-1e-10
    j=rng.normal(size=(30,7));e=rng.normal(size=(30,4));t=np.diag(np.arange(1,8))[:3]
    m=t@np.linalg.pinv(j)@e;u,s,vh=np.linalg.svd(m,full_matrices=False)
    top=vh[0];bias=np.linalg.norm(m@top)**2
    assert abs(bias-s[0]**2)<1e-10
    checks.append(dict(case=case,information_identity_error=float(err)))
a=2.;good=1/(1+a*a);bad=good+a*a*9/(1+a*a)**2
assert good<1<bad
result=dict(passed=True,checks=checks,scalar_correct_model_mse=good,
            scalar_misspecified_mse=bad,old_mse=1.,
            scope='algebra only; no claim of known electromagnetic discrepancy envelope')
(OUT/'results/frequency_information_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
