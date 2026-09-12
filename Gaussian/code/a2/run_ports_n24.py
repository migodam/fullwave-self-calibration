import json,hashlib,sys,platform
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.ports import GaussianComponent,LocalTNetwork,render,ownership_blocks
cfg={'n':24,'frequency_hz':2.25e9,'n_tx':6,'n_rx':20,'rank':3,'enrich':1,'full_rtol':1e-9,'port_rtol':1e-8,'components':[[.32,.018,-.065,0,.045,.05,.15],[.28,.015,.065,.012,.05,.04,-.18],[.22,.012,0,.075,.038,.048,.25]]}
cs=[GaussianComponent(a+1j*b,(x,y),(sx,sy),ang) for a,b,x,y,sx,sy,ang in cfg['components']];v=VIE(Geometry(n=24,n_tx=6,n_rx=20),cfg['frequency_hz']);chi=render(cs,v.points);f=v.forward(chi,rtol=cfg['full_rtol'],blocks=ownership_blocks(cs,v.points));net=LocalTNetwork(v,cs);et=v.E+v.D.matmat(f['current']);cis=[x[:,None]*et for x in net.chis];p=net.solve_ports(3,rtol=cfg['port_rtol'],enrich=1)
o={**cfg,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'environment':{'python':sys.version,'platform':platform.platform()},'branch':'matrix_free_local_response','identity_errors':net.full_component_identity_error(cis,v.E),'network_dimension':p['network_dimension'],'scattered_relative_error':float(np.linalg.norm(p['scattered']-f['scattered'])/np.linalg.norm(f['scattered'])),'full_current_state_residual_relative':float(np.linalg.norm(f['current']-chi[:,None]*(v.E+v.D.matmat(f['current'])))/np.linalg.norm(f['current']))};out=ROOT/'runs/a2/ports';out.mkdir(parents=True,exist_ok=True);(out/'ports_n24_regression.json').write_text(json.dumps(o,indent=2));print(json.dumps(o))
