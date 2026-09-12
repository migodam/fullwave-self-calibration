from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np,torch
torch.set_num_threads(2)
from a2.physics import Geometry,VIE
from a2.gpu.torch_vie import TorchVIE
ROOT=Path(__file__).resolve().parents[3];out=ROOT/'runs/a2/gpu';out.mkdir(parents=True,exist_ok=True)
np.random.seed(42);torch.manual_seed(42);v=VIE(Geometry(n=16,n_tx=4,n_rx=16),2.25e9);chi=np.exp(-np.sum(v.points**2,1)/.008)*(.6+.09j);ref=v.forward(chi,rtol=1e-11)['scattered'];records=[]
for device,dtype in [('cpu',torch.complex128)]+([('cuda',torch.complex64)] if torch.cuda.is_available() else []):
 t=TorchVIE(v,device,dtype);c=torch.tensor(chi,dtype=dtype,device=device,requires_grad=True);f=t.scattered(c);err=np.linalg.norm(f.detach().cpu().numpy()-ref)/np.linalg.norm(ref);target=torch.tensor(ref*.8,dtype=dtype,device=device);loss=.5*abs(f-target).square().sum();loss.backward();d=torch.tensor(np.random.normal(size=chi.size)+1j*np.random.normal(size=chi.size),dtype=dtype,device=device);d*=.1
 eps=2e-3 if device=='cuda' else 1e-5
 with torch.no_grad():
  fp=t.forward(c+eps*d)[0];fm=t.forward(c-eps*d)[0];fd=(.5*abs(fp-target).square().sum()-.5*abs(fm-target).square().sum())/(2*eps);ad=(c.grad.conj()*d).sum().real;grad_error=float(abs(fd-ad)/max(abs(fd),abs(ad),1e-20))
 z=torch.randn(t.N,2,dtype=dtype,device=device);w=torch.randn_like(z);aa=(t.D(z).conj()*w).sum();bb=(z.conj()*t.D(w,True)).sum();adj=float(abs(aa-bb)/max(abs(aa),abs(bb),1e-20));assert err<(1e-4 if device=='cuda' else 1e-7);assert grad_error<(5e-3 if device=='cuda' else 1e-5);assert adj<(1e-5 if device=='cuda' else 1e-10)
 records.append(dict(device=device,scattered_relative=float(err),gradient_relative=grad_error,adjoint_relative=adj,counts=t.counts,last=t.last,peak_bytes=torch.cuda.max_memory_allocated() if device=='cuda' else None))
r=dict(records=records,torch=torch.__version__,gpu=torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('torch_vie.py')]},scope='N16 direct CPU reference, implicit gradient and adjoint checks; not large imaging evidence');(out/'backend_checks.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
