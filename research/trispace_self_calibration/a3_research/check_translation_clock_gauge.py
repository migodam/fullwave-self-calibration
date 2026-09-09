"""Independent full multiple-scattering covariance check; no inversion claim."""
import json
from pathlib import Path
import numpy as np
from maxwell3d import treams_field,receivers,illuminations

OUT=Path(__file__).resolve().parent
centers=np.array([[-.2,0,0],[.2,.03,.04]])
radii=[.12,.1];eps=[2.5+.03j,3.2+.04j];rx=receivers()
directions=np.array([d for d,p in illuminations()])
difference=directions[1:]-directions[0]
assert np.linalg.matrix_rank(difference)==1
tetra=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)
assert np.linalg.matrix_rank(tetra[1:]-tetra[0])==3
rows=[]
for k in [3.,9.,18.]:
    base,_=treams_field(centers,radii,eps,k,rx,lmax=6)
    for name,t in [('transverse',np.array([0.,.08,0.])),
                   ('common_clock',np.array([.05,0.,.05])),
                   ('non_gauge_translation',np.array([.06,0.,0.]))]:
        shifted,_=treams_field(centers+t,radii,eps,k,rx+t,lmax=6)
        phase=np.exp(1j*k*(directions@t))
        covariance_error=float(np.linalg.norm(shifted-base*phase)/np.linalg.norm(base))
        assert covariance_error<1e-11
        c=float(directions[0]@t)
        compensated=shifted*np.exp(-1j*k*c)
        common_error=float(np.linalg.norm(compensated-base)/np.linalg.norm(base))
        in_family=bool(np.linalg.norm(difference@t)<1e-14)
        if in_family:assert common_error<1e-11
        individual=shifted*np.exp(-1j*k*(directions@t))
        individual_error=float(np.linalg.norm(individual-base)/np.linalg.norm(base))
        assert individual_error<1e-11
        rows.append(dict(k=k,translation_name=name,translation=t.tolist(),
            common_clock_gauge=in_family,field_covariance_error=covariance_error,
            fixed_compensating_common_clock_error=common_error,
            individual_delay_compensation_error=individual_error))
result=dict(passed=True,rows=rows,two_direction_difference_rank=1,
            tetrahedral_direction_difference_rank=3,
            scope='full-wave covariance check; non-gauge row is not a proof of absence of every other ambiguity')
(OUT/'results/translation_clock_gauge_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
