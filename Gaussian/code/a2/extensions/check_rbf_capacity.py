"""Linear feasibility separates logit-family capacity from optimization failure."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'code'))
from a2.extensions.probability import all_states
from a2.extensions.probability_run import centres
c=centres();rc=np.array([[-.085,0],[.085,0],[0,.085]]);P=np.c_[np.ones(8),np.exp(-np.sum((c[:,None]-rc[None])**2,axis=2)/(2*.075**2))]
rows=[]
for z in all_states():
 A=-(2*z[:,None]-1)*P;o=linprog(np.zeros(4),A_ub=A,b_ub=-np.ones(8),bounds=[(None,None)]*4,method='highs');feasible=o.success
 if feasible:assert np.min((2*z-1)*(P@o.x))>1-1e-7
 rows.append(dict(state=z.tolist(),separable=bool(feasible),solver_status=o.status))
r=dict(scope='Post-hoc exact finite LP diagnostic: whether fixed Gaussian logits can concentrate on each binary pattern, not a Bayesian optimality test',feature_rank=int(np.linalg.matrix_rank(P)),separable_count=sum(x['separable'] for x in rows),total=256,states=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());p=ROOT/'runs/a2/extensions/rbf_capacity.json';p.write_text(json.dumps(r,indent=2));print({k:v for k,v in r.items() if k!='states'})
