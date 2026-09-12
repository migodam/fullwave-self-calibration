"""CUDA full-wave representation campaign; a bounded comparison, not a SOM test.

Default production run: 30 fixed cases x {Gaussian K=16,64,144; voxel N128}.
The N192 cell-integrated fields are generated once per case; optimisation uses a
separate N128 point-kernel model, deliberately exposing discretisation mismatch.
"""
import argparse, hashlib, json, os, platform, time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='2'
import numpy as np
import torch
torch.set_num_threads(2)
import torch.nn.functional as F
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
import sys; sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.gpu.torch_vie import TorchVIE,TorchMultiVIE

DATA_BACKENDS={}
INVERSE_BACKENDS={}
FREQ=(1.25e9,2.0e9,2.75e9); LOSS=1+.15j; SIDE=.4
OUT=ROOT/'runs/a2/gpu/representation_campaign'; FIG=ROOT/'figures/a2/gpu'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic_json(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,default=str));os.replace(tmp,path)
def atomic_torch(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp');torch.save(obj,tmp);os.replace(tmp,path)
def geom(n): return Geometry(n=n,n_tx=12,n_rx=64,aperture='half')
def grid_t(n,dev):
    x=torch.linspace(-.2+.2/n,.2-.2/n,n,device=dev);xx,yy=torch.meshgrid(x,x,indexing='ij');return torch.stack((xx.reshape(-1),yy.reshape(-1)),1)
def inv_softplus(x): return torch.log(torch.expm1(torch.clamp(x,min=1e-7)))
def smooth(x,n):
    z=x.reshape(n,n)
    # Same discrete pixel-difference penalty for every method; no TV is used.
    return ((z[1:]-z[:-1]).square().mean()+(z[:,1:]-z[:,:-1]).square().mean())
def receiver_masks(device):
    rx=torch.arange(64,device=device);return (rx%2==0),(rx%2==1)

def truth_real(points,case):
    """Three fixed families, with no Gaussian-centre grid used for generation."""
    rng=np.random.default_rng(20260915+case); fam=case//10; x,y=points[:,0],points[:,1]
    if fam==0: # 16--32 off-grid Gaussian inclusions
        z=torch.zeros_like(x)
        for _ in range(int(rng.integers(16,33))):
            c=torch.tensor(rng.uniform(-.135,.135,2),device=points.device);s=torch.tensor(rng.uniform(.006,.022,2),device=points.device);a=float(rng.uniform(.08,.30))
            z=z+a*torch.exp(-.5*(((points-c)/s).square().sum(1)))
    elif fam==1: # sharp disks and rectangles
        z=torch.zeros_like(x)
        for j in range(8):
            c=torch.tensor(rng.uniform(-.13,.13,2),device=points.device);a=float(rng.uniform(.18,.55))
            if j%2: z=z+a*((x-c[0]).abs()<rng.uniform(.012,.035))*((y-c[1]).abs()<rng.uniform(.012,.035))
            else: z=z+a*(((x-c[0]).square()+(y-c[1]).square())<rng.uniform(.010,.030)**2)
    else: # thin curving material ribbons with explicitly empty holes
        z=.38*torch.exp(-((y-.055*torch.sin(21*x+.6*case))/.012).square())+.30*torch.exp(-((y+.055*torch.sin(17*x-.3*case))/.010).square())
        for _ in range(4):
            c=torch.tensor(rng.uniform(-.11,.11,2),device=points.device);z=z-.45*torch.exp(-((x-c[0]).square()+(y-c[1]).square())/.00028)
        z=torch.clamp(z,min=0)
    return z

class GaussianRep(torch.nn.Module):
    def __init__(self,k,points,case):
        super().__init__();self.k=k;self.register_buffer('points',points)
        m=int(np.ceil(np.sqrt(k)));a=torch.linspace(-.14,.14,m,device=points.device);xx,yy=torch.meshgrid(a,a,indexing='ij');c=torch.stack((xx.reshape(-1),yy.reshape(-1)),1)[:k]
        # Fixed uniform-grid initialization is independent of the truth.
        self.raw_amp=torch.nn.Parameter(torch.full((k,),float(inv_softplus(torch.tensor(.025))),device=points.device))
        self.raw_center=torch.nn.Parameter(torch.atanh(torch.clamp(c/.18,-.999,.999)))
        raw_diag=float(inv_softplus(torch.tensor((.018-.004)/.04)))
        self.raw_l=torch.nn.Parameter(torch.stack((torch.full((k,),raw_diag,device=points.device),torch.zeros(k,device=points.device),torch.full((k,),raw_diag,device=points.device)),1))
    def render(self):
        c=.18*torch.tanh(self.raw_center);d=self.points[:,None,:]-c[None]
        l00=.004+.04*F.softplus(self.raw_l[:,0]);l10=.02*self.raw_l[:,1];l11=.004+.04*F.softplus(self.raw_l[:,2])
        # q=||L^{-1}(x-centre)||², L lower triangular, hence SPD covariance LLᵀ.
        u0=d[:,:,0]/l00;u1=(d[:,:,1]-l10*u0)/l11;q=u0.square()+u1.square()
        return (F.softplus(self.raw_amp)[None]*torch.exp(-.5*q)).sum(1)

class VoxelRep(torch.nn.Module):
    def __init__(self,initial): super().__init__();self.raw=torch.nn.Parameter(inv_softplus(initial.detach().clamp_min(1e-7)))
    def render(self): return F.softplus(self.raw)

def make_backends(n,device):
    key=(n,str(device))
    if key not in INVERSE_BACKENDS:INVERSE_BACKENDS[key]=[TorchMultiVIE([TorchVIE(VIE(geom(n),hz),device=device,dtype=torch.complex64) for hz in FREQ])]
    return INVERSE_BACKENDS[key]
def predict(backends,real):
    chi=real.to(torch.complex64)*torch.as_tensor(LOSS,dtype=torch.complex64,device=real.device)
    return list(backends[0].scattered(chi).unbind(0))
def counts(backends):
    keys=('forward_rhs','adjoint_rhs','linear_iterations','solves');return {k:int(sum(b.counts[k] for b in backends)) for k in keys}
def subtract(a,b): return {k:a[k]-b[k] for k in a}

def make_data(case,data_n,device):
    t=time.perf_counter();pts=grid_t(data_n,device);real=truth_real(pts,case);key=(data_n,str(device));first_setup=key not in DATA_BACKENDS
    if first_setup:DATA_BACKENDS[key]=[TorchVIE(VIE(geom(data_n),hz,cell_integrated=True),device=device,dtype=torch.complex64) for hz in FREQ]
    backs=DATA_BACKENDS[key];before=counts(backs);fields=[]
    with torch.no_grad():
        for b in backs:
            fields.append(b.scattered(real.to(torch.complex64)*torch.as_tensor(LOSS,dtype=torch.complex64,device=device)).detach())
    return real.detach(),fields,{'seconds_including_first_setup':time.perf_counter()-t,'first_setup':first_setup,'counts':subtract(counts(backs),before)}

def evaluate(backends,real,data,held,n):
    with torch.no_grad():
        p=predict(backends,real);den=sum(torch.sum(torch.abs(y[held])**2) for y in data);num=sum(torch.sum(torch.abs(q[held]-y[held])**2) for q,y in zip(p,data))
        # N128 material difference is reported after data truth was bilinearly sampled.
        return float(torch.sqrt(num/den).cpu()),p
def iou(a,b):
    aa=a>.12;bb=b>.12;return float((aa&bb).sum().float()/torch.clamp((aa|bb).sum(),min=1))

def run_method(case,label,rep,data,truth_n,backends,args,config_hash,run_dir):
    ck=run_dir/f'case_{case:02d}_{label}.pt';js=run_dir/f'case_{case:02d}_{label}.json';train,held=receiver_masks(rep.render().device)
    opt=torch.optim.Adam(rep.parameters(),lr=args.lr);traj=[];start=0;prior_seconds=0.;prior_counts={k:0 for k in counts(backends)}
    if ck.exists():
        saved=torch.load(ck,map_location=rep.render().device,weights_only=False)
        if saved['config_hash']!=config_hash: raise RuntimeError(f'checkpoint config mismatch: {ck}')
        rep.load_state_dict(saved['model']);opt.load_state_dict(saved['optimizer']);traj=saved['trajectory'];start=saved['step']
    if js.exists():
        previous=json.loads(js.read_text())
        if previous.get('config_hash')==config_hash and previous['steps_completed']==args.steps:return previous
    if ck.exists():
        prior_seconds=saved.get('cumulative_seconds',0.);prior_counts=saved.get('cumulative_counts',prior_counts)
    before=counts(backends);t=time.perf_counter()
    den=sum(torch.sum(torch.abs(y[train])**2) for y in data).detach()
    for step in range(start,args.steps):
        opt.zero_grad(set_to_none=True);real=rep.render();p=predict(backends,real)
        fit=sum(torch.sum(torch.abs(q[train]-y[train])**2) for q,y in zip(p,data))/den
        reg=1e-4*smooth(real,args.n);loss=fit+reg;loss.backward();opt.step()
        traj.append({'step':step+1,'loss':float(loss.detach().cpu()),'fit':float(fit.detach().cpu()),'smooth':float(reg.detach().cpu())})
        atomic_torch(ck,{'config_hash':config_hash,'step':step+1,'model':rep.state_dict(),'optimizer':opt.state_dict(),'trajectory':traj,'cumulative_seconds':prior_seconds+time.perf_counter()-t,'cumulative_counts':{k:prior_counts[k]+counts(backends)[k]-before[k] for k in before}})
    torch.cuda.synchronize() if rep.render().device.type=='cuda' else None
    seconds=prior_seconds+time.perf_counter()-t;real=rep.render().detach();held_error,p=evaluate(backends,real,data,held,args.n)
    truth_small=F.interpolate(truth_n.reshape(1,1,int(np.sqrt(truth_n.numel())),int(np.sqrt(truth_n.numel()))),size=(args.n,args.n),mode='bilinear',align_corners=False)[0,0].reshape(-1)
    result={'config_hash':config_hash,'parameter_count':sum(p.numel() for p in rep.parameters()),'case':case,'method':label,'steps_completed':args.steps,'material_l2_relative':float(torch.linalg.vector_norm(real-truth_small)/torch.linalg.vector_norm(truth_small)),'held_clean_field_relative':held_error,'threshold_iou_0p12':iou(real,truth_small),'actual_optimization_seconds':seconds,'counts_delta':{k:prior_counts[k]+counts(backends)[k]-before[k] for k in before},'trajectory':traj,'regularizer':'1e-4 times mean squared adjacent-pixel differences; same discrete penalty for Gaussian and voxel'}
    np.savez_compressed(run_dir/f'case_{case:02d}_{label}_arrays.npz',truth=truth_small.cpu().numpy(),reconstruction=real.cpu().numpy())
    atomic_json(js,result);print('fit_done',case,label,result['material_l2_relative'],result['held_clean_field_relative'],seconds,flush=True);return result

def plot_case(case,truth,results,n):
    names=['truth']+[r['method'] for r in results];arr=[truth.cpu().numpy()]+[np.load(OUT/f'case_{case:02d}_{r["method"]}_arrays.npz')['reconstruction'] for r in results]
    fig,ax=plt.subplots(1,len(arr),figsize=(3*len(arr),3),sharex=True,sharey=True,constrained_layout=True);vmax=max(float(np.max(z)) for z in arr)
    for a,z,name in zip(np.atleast_1d(ax),arr,names):
        im=a.imshow(z.reshape(n,n).T,origin='lower',extent=(-20,20,-20,20),vmin=0,vmax=vmax,cmap='viridis');a.set(title=name,xlabel='x (cm)',ylabel='y (cm)')
    fig.colorbar(im,ax=np.atleast_1d(ax).tolist(),label='real contrast, common scale');fig.savefig(FIG/f'representation_case_{case:02d}.png',dpi=160);plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cases',type=int,nargs='+',default=list(range(30)));ap.add_argument('--steps',type=int,default=80);ap.add_argument('--n',type=int,default=128);ap.add_argument('--data-n',type=int);ap.add_argument('--lr',type=float,default=.03);ap.add_argument('--device',default='cuda');ap.add_argument('--run-name',default='main');args=ap.parse_args()
    if not torch.cuda.is_available() and args.device.startswith('cuda'): raise RuntimeError('CUDA is required for this campaign')
    if args.n<8 or any(c<0 or c>=30 for c in args.cases): raise ValueError('n>=8 and cases must be in 0..29')
    global OUT,FIG
    OUT=OUT/args.run_name;FIG=FIG/args.run_name
    device=torch.device(args.device);data_n=args.data_n or (192 if args.n==128 else args.n);OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    sources=[Path(__file__),Path(__file__).with_name('torch_vie.py'),ROOT/'code/a2/physics.py'];cfg={'case_universe':list(range(30)),'lr':args.lr,'noise':'noiseless development comparison','source_hashes':{p.name:sha(p) for p in sources},'steps':args.steps,'inverse_n':args.n,'truth_n':data_n,'frequencies_hz':FREQ,'geometry':'12 Tx, 64 half-aperture Rx; even Rx train and odd Rx held','loss_factor':[1,.15],'methods':['gaussian_K16','gaussian_K64','gaussian_K144','voxel'],'seed_rule':'20260915 + case','scope':'ordinary positive Gaussian representation versus voxel; no SOM contribution. Optimizer, regularizer and step count remain a declared development setting, not a fairness proof.'}
    config_hash=hashlib.sha256(json.dumps(cfg,sort_keys=True,default=str).encode()).hexdigest();existing=OUT/'frozen_config.json'
    if existing.exists() and json.loads(existing.read_text())['config_hash']!=config_hash:raise RuntimeError('Run name already contains different configuration; use another --run-name')
    atomic_json(OUT/'frozen_config.json',{'config':cfg,'config_hash':config_hash,'source_hashes':{str(p.relative_to(ROOT)):sha(p) for p in sources},'torch':torch.__version__,'cuda':torch.cuda.get_device_name(device) if device.type=='cuda' else None,'platform':platform.platform()})
    summary=[]
    for case in args.cases:
        print('case_start',case,flush=True)
        torch.manual_seed(20260915+case);np.random.seed(20260915+case);data_truth,data,data_meta=make_data(case,data_n,device);np.savez_compressed(OUT/f'case_{case:02d}_data.npz',truth=data_truth.cpu().numpy(),fields=np.stack([x.cpu().numpy() for x in data]));truth_inv=F.interpolate(data_truth.reshape(1,1,data_n,data_n),size=(args.n,args.n),mode='bilinear',align_corners=False)[0,0].reshape(-1)
        setup=time.perf_counter();backs=make_backends(args.n,device);pts=grid_t(args.n,device);base64=GaussianRep(64,pts,case);setup_seconds=time.perf_counter()-setup
        runs=[]
        for k in (16,64,144): runs.append(run_method(case,f'gaussian_K{k}',GaussianRep(k,pts,case),data,truth_inv,backs,args,config_hash,OUT))
        runs.append(run_method(case,'voxel',VoxelRep(base64.render()),data,truth_inv,backs,args,config_hash,OUT));plot_case(case,truth_inv,runs,args.n)
        rec={'case':case,'family':['off_grid_gaussian_mixture','sharp_multi_inclusion','curves_with_holes'][case//10],'data_generation':data_meta,'inverse_setup_seconds':setup_seconds,'methods':runs,'cuda_peak_bytes':int(torch.cuda.max_memory_allocated(device)) if device.type=='cuda' else None};atomic_json(OUT/f'case_{case:02d}_summary.json',rec);summary.append(rec)
    atomic_json(OUT/'summary_last_invocation.json',{'config_hash':config_hash,'cases':summary});print(json.dumps({'completed_cases':args.cases,'config_hash':config_hash},indent=2))
if __name__=='__main__': main()
