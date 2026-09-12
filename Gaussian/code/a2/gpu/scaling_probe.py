"""CUDA state-size scaling; full solve costs, not an inversion advantage claim."""
from pathlib import Path
import sys,os,time,json,hashlib,platform
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='2'
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np,torch
from a2.physics import VIE,Geometry
from a2.gpu.torch_vie import TorchVIE
torch.set_num_threads(2);ROOT=Path(__file__).resolve().parents[3];out=ROOT/'runs/a2/gpu';out.mkdir(parents=True,exist_ok=True);records=[]
for n in [64,128,256]:
 torch.cuda.empty_cache();torch.cuda.reset_peak_memory_stats();tic=time.perf_counter();v=VIE(Geometry(n=n,n_tx=12,n_rx=64,aperture='half'),2.75e9);b=TorchVIE(v,'cuda',torch.complex64);p=v.points;c0=(.25*np.exp(-np.sum((p-[.04,.02])**2,1)/.008)+.15*np.exp(-np.sum((p+[.08,.03])**2,1)/.003))*(1+.15j);c=torch.tensor(c0,dtype=torch.complex64,device='cuda',requires_grad=True);torch.cuda.synchronize();setup=time.perf_counter()-tic;times=[]
 for rep in range(4):
  tic=time.perf_counter();f=b.scattered(c);loss=.5*f.abs().square().sum();loss.backward();torch.cuda.synchronize();times.append(time.perf_counter()-tic);c.grad=None
 records.append(dict(n=n,state_dimension=n*n,measurements_complex=12*64,setup_seconds=setup,warmup_seconds=times[0],forward_adjoint_seconds=times[1:],peak_gpu_bytes=torch.cuda.max_memory_allocated(),counts=b.counts,last=b.last));print(records[-1],flush=True);del b,v,c,f,loss
r=dict(records=records,torch=torch.__version__,gpu=torch.cuda.get_device_name(0),python=platform.python_version(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Single frequency known material; includes physical forward and implicit adjoint, no image optimizer');(out/'scaling.json').write_text(json.dumps(r,indent=2))
