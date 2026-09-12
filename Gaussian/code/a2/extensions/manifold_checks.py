import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'code'))
from a2.extensions.manifold import *
r=Path(__file__).resolve().parents[3]/'runs/a2/extensions/manifold';r.mkdir(parents=True,exist_ok=True)
x=np.random.default_rng(1).normal(size=(100,2));S=covariance(np.log(.03),np.log(.05),.4);H=SYM_BASIS[2]*.01
v=covariance_direction(x,.5,(.01,-.02),S,H);f=fd_covariance_direction(x,.5,(.01,-.02),S,H)
o={'peak_amplitude_covariance_fd_relative':float(np.linalg.norm(v-f)/np.linalg.norm(v)),'spd_min_eigenvalue':float(np.linalg.eigvalsh(retract_spd(S,H))[0]),'old_angle_chart_relative':chart_equivalence(x,.5,.01,-.02,np.log(.03),np.log(.05),.4),'trace_term_used':False}
(r/'checks.json').write_text(json.dumps(o,indent=2));print(json.dumps(o))
