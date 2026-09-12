from pathlib import Path
import sys,json,time,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np,torch
torch.set_num_threads(2)
from a2.physics import Geometry,VIE
from a2.gpu.torch_vie import TorchVIE,TorchMultiVIE
out=Path(__file__).resolve().parents[3]/'runs/a2/gpu';out.mkdir(parents=True,exist_ok=True);results=[]
for device,dtype in [('cpu',torch.complex128)]+([('cuda',torch.complex64)] if torch.cuda.is_available() else []):
 v=[VIE(Geometry(n=16,n_tx=4,n_rx=16),f) for f in [1.25e9,2e9,2.75e9]];backs=[TorchVIE(q,device,dtype) for q in v];b=TorchMultiVIE(backs);chi=.4*np.exp(-np.sum(v[0].points**2,1)/.006)*(1+.15j);c=torch.tensor(chi,dtype=dtype,device=device,requires_grad=True)
 f=b.scattered(c);ref=torch.stack([q.scattered(c) for q in backs]);err=float(torch.linalg.vector_norm(f-ref)/torch.linalg.vector_norm(ref));g1=torch.autograd.grad(.5*f.abs().square().sum(),c)[0];g2=torch.autograd.grad(.5*ref.abs().square().sum(),c)[0];ge=float(torch.linalg.vector_norm(g1-g2)/torch.linalg.vector_norm(g2));assert err<1e-5 and ge<1e-4
 results.append(dict(device=device,field_error=err,gradient_error=ge,counts=b.counts))
r=dict(results=results,hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('torch_vie.py')]});(out/'batch_checks.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
