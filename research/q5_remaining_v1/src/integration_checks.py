"""Read-only artifact checks plus independent checks of integration-critical formulas."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.special import spherical_jn as jn, spherical_yn as yn

ROOT=Path(__file__).resolve().parents[1]
W=ROOT/'pipeline/boundary_20260910_171833/experiment_q5_boundary_reconstruction'
s=json.loads((W/'results/FINAL_SUMMARY.json').read_text())
bad=[r['path'] for r in s['artifacts'] if hashlib.sha256((W/r['path']).read_bytes()).hexdigest()!=r['sha256']]
assert not bad
loss=np.array([.03,.05]);q0,q1=30.,40.
g0=1.;g1=(q0+1j)/(q1+1j)
assert .75<=abs(g1)<=1.25
assert np.max(abs(g0*(q0+1j)*loss-g1*(q1+1j)*loss))<1e-14
assert np.all(1+q0*loss>=[1.5,2]) and np.all(1+q1*loss<=[4,5])
rows=[]
for z in (.2,.4,.7):
    j=jn(1,z);jp=jn(1,z,True);jpp=-2*jp/z-(1-2/z**2)*j
    dj=jp+j/z;dy=yn(1,z,True)+yn(1,z)/z
    first=-jp/z-j+j/z**2
    second=-jpp/z-jp+2*jp/z**2-2*j/z**3
    closed=-np.sin(z)/z-3*np.cos(z)/z**2+7*np.sin(z)/z**3+12*np.cos(z)/z**4-12*np.sin(z)/z**5
    assert abs(second-closed)<1e-10
    rows.append({'z':z,'wrong_first_identity_error':float(abs(first-(-dy+2*dj/z))),
                 'wrong_second_identity_error':float(abs(second-(-jpp/z-2*jp/z+2*jp/z**2-2*j/z**3))),
                 'correct_second_vs_worker_code_error':float(abs(second-closed))})
result={'artifact_hashes_pass':True,'artifacts_checked':len(s['artifacts']),
        'worker_constant_categories':s['counts'],'admissible_born_pair':{'epsilon0':(1+q0*loss).tolist(),
        'epsilon1':(1+q1*loss).tolist(),'gain1_modulus':abs(g1)},'derivative_checks':rows,
        'pipeline_exit_code':0,'native_development_rounds':1,'native_turn_cap':18,
        'native_evaluator_decision':'abandon','native_writeup_executed':False,
        'scientific_status':'incomplete; bounded verification integrated, not paper acceptance'}
p=ROOT/'results/integration_audit.json'
if p.exists():raise RuntimeError('Preserve prior integration audit')
p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='worker_constant_categories'},indent=2))
