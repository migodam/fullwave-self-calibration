"""Native SPD full-wave tangent and finite-secant comparison, frozen common budget."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[key]='1'
from pathlib import Path
import sys,json,time,hashlib,platform
import numpy as np
from scipy.linalg import expm
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import VIE,Geometry
from a2.extensions.manifold import SYM_BASIS,peak_gaussian
OUT=ROOT/'runs/a2/extensions/manifold_v2';OUT.mkdir(parents=True,exist_ok=True)
FIG=ROOT/'figures/a2/extensions';FIG.mkdir(parents=True,exist_ok=True)
CFG=dict(seed=20260913,n_inverse=32,n_data=64,frequencies=[1.5e9,2.75e9],n_tx=6,n_rx=24,noise_relative=.01,families=['separated','close','sharp'],cases_per_family=10,max_outer=20,max_rhs=3600,trust_radius=.7,methods=['full_lm','spectral','weak_secant','random_secant'],scope='2D scalar; two fixed Gaussian components; known loss ratio; fixed common ROI initialization; half receive/full transmit aperture; independent grid and quadrature')
def sqrtm(S):
 w,u=np.linalg.eigh(S);return (u*np.sqrt(w))@u.T

def material(points,state,derivatives=False):
 chi=np.zeros(len(points));cols=[]
 for a,mu,S in state:
  root=sqrtm(S);inv=np.linalg.inv(S);q=points-mu;g=peak_gaussian(points,a,mu,S);chi+=g
  if derivatives:
   cols.extend([g,g*(q@inv[:,0])*.03,g*(q@inv[:,1])*.03])
   for E in SYM_BASIS:
    H=root@E@root;cols.append(.5*g*np.einsum('ni,ij,jk,kl,nl->n',q,inv,H,inv,q))
 return (chi*(1+.15j),np.asarray(cols).T*(1+.15j)) if derivatives else chi*(1+.15j)

def retract(state,step):
 out=[]
 for k,(a,mu,S) in enumerate(state):
  d=step[6*k:6*k+6];root=sqrtm(S);H=sum(d[3+i]*E for i,E in enumerate(SYM_BASIS));Sp=root@expm(H)@root
  ap=a*np.exp(d[0]);mp=mu+.03*d[1:3];w=np.linalg.eigvalsh(Sp)
  if not (.015<ap<3 and np.max(abs(mp))<.16 and w.min()>.009**2 and w.max()<.12**2):return None
  out.append((ap,mp,Sp))
 return out

def initial():return [(.35,np.array([-.055,0.]),np.eye(2)*.04**2),(.35,np.array([.055,0.]),np.eye(2)*.04**2)]

def truth(points,family,rng):
 centers=np.array([[-.065,.015],[.06,-.02]]) if family!='close' else np.array([[-.022,.008],[.022,-.008]])
 centers+=rng.normal(0,.007,(2,2));amps=rng.uniform(.4,.85,2);width=rng.uniform(.023,.039,(2,2));angles=rng.uniform(-1,1,2);out=np.zeros(len(points))
 for mu,a,w,t in zip(centers,amps,width,angles):
  R=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]]);q=(points-mu)@R
  out+=a*((np.sum((q/(w*1.5))**2,axis=1)<1) if family=='sharp' else np.exp(-.5*np.sum((q/w)**2,axis=1)))
 return out*(1+.15j)

def real(z):return np.r_[z.real,z.imag]
def make_solvers(n,integrated=False):return [VIE(Geometry(n=n,n_tx=6,n_rx=24,aperture='half'),f,cell_integrated=integrated) for f in CFG['frequencies']]
class Model:
 def __init__(self,solvers):self.v=solvers;self.rhs=0;self.forward_calls=0;self.tangent_calls=0
 def field(self,state):
  fw=[v.forward(material(v.points,state),rtol=1e-8) for v in self.v];self.rhs+=12;self.forward_calls+=1
  return np.concatenate([o['scattered'].ravel() for o in fw]),fw
 def jac(self,state,fw):
  blocks=[]
  for v,o in zip(self.v,fw):
   c,dc=material(v.points,state,True);blocks.append(np.column_stack([v.material_tangent(c,d,forward=o,rtol=1e-8)['scattered'].ravel() for d in dc.T]))
  self.rhs+=144;self.tangent_calls+=1;return np.concatenate(blocks)

def fit(solvers,y,mask,method,seed):
 model=Model(solvers);state=initial();tic=time.perf_counter();f,fw=model.field(state);scale=np.linalg.norm(y[mask]);r=real((y-f)[mask])/scale;hist=[];rng=np.random.default_rng(seed);damping=.003
 for iteration in range(CFG['max_outer']):
  # Every method pays actual forward/tangent RHS count; at most five trial fields.
  if model.rhs+144+60>CFG['max_rhs']:break
  Jc=model.jac(state,fw);J=np.r_[Jc[mask].real,Jc[mask].imag]/scale;U,s,Vh=np.linalg.svd(J,full_matrices=False)
  rank=max(1,int(np.sum(s>.08*s[0])));rank=len(s) if method=='full_lm' else rank;V=Vh[:rank].T;Js=J@V
  def visible(rr):return V@np.linalg.solve(Js.T@Js+damping*np.eye(rank),Js.T@rr)
  proposals=[visible(r)];probe_count=0
  if method in ('weak_secant','random_secant'):
   if method=='weak_secant':direction=Vh[min(rank,len(s)-1)].copy()
   else:direction=rng.normal(size=len(s));direction/=np.linalg.norm(direction)
   for sign in (-1,1):
    h=sign*.5*direction;sp=retract(state,h)
    if sp is None:continue
    fp,_=model.field(sp);probe_count+=1;delta=real((fp-f)[mask])/scale
    proposals.append(h+visible(r-delta))
  best=(float(r@r),state,f,fw,r);accepted=False
  for step in proposals:
   norm=np.linalg.norm(step)
   if norm>CFG['trust_radius']:step*=CFG['trust_radius']/norm
   sp=retract(state,step)
   if sp is None:continue
   fp,fwp=model.field(sp);rp=real((y-fp)[mask])/scale;val=float(rp@rp)
   if val<best[0]:best=(val,sp,fp,fwp,rp);accepted=True
  hist.append(dict(iteration=iteration,loss=best[0],rank=rank,rhs=model.rhs,accepted=accepted,probe_fields=probe_count))
  if accepted:_,state,f,fw,r=best;damping=max(1e-5,damping*.7)
  else:damping=min(10,damping*4)
 return state,f,dict(seconds=time.perf_counter()-tic,rhs=model.rhs,forward_calls=model.forward_calls,tangent_calls=model.tangent_calls,history=hist,final_train_relative=float(np.linalg.norm((f-y)[mask])/scale))

def main():
 files=[Path(__file__),ROOT/'code/a2/extensions/manifold.py',ROOT/'code/a2/physics.py']
 (OUT/'frozen_config.json').write_text(json.dumps(dict(config=CFG,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},environment=dict(python=platform.python_version(),numpy=np.__version__)),indent=2))
 inv=make_solvers(32);data=make_solvers(64,True);mask=np.tile(np.repeat(np.arange(24)%2==0,6),2);records=[];images=[]
 for family in CFG['families']:
  for rep in range(10):
   seed=CFG['seed']+100*CFG['families'].index(family)+rep;rng=np.random.default_rng(seed);state_rng=rng.bit_generator.state
   ctruth=truth(data[0].points,family,rng);rng.bit_generator.state=state_rng;itruth=truth(inv[0].points,family,rng)
   tic=time.perf_counter();clean=np.concatenate([v.forward(ctruth,rtol=1e-9)['scattered'].ravel() for v in data]);data_seconds=time.perf_counter()-tic
   sigma=.01*np.linalg.norm(clean)/np.sqrt(len(clean));y=clean+sigma/np.sqrt(2)*(rng.normal(size=clean.size)+1j*rng.normal(size=clean.size));case=dict(family=family,rep=rep,data_seconds=data_seconds,methods={});im=[itruth.real]
   for method in CFG['methods']:
    st,f,meta=fit(inv,y,mask,method,seed);c=material(inv[0].points,st);meta.update(material_relative_l2=float(np.linalg.norm(c-itruth)/np.linalg.norm(itruth)),heldout_clean_relative=float(np.linalg.norm((f-clean)[~mask])/np.linalg.norm(clean[~mask])),train_clean_relative=float(np.linalg.norm((f-clean)[mask])/np.linalg.norm(clean[mask])))
    case['methods'][method]=meta;im.append(c.real)
   records.append(case);images.append(im);(OUT/'results.json').write_text(json.dumps(dict(config=CFG,records=records),indent=2));np.savez_compressed(OUT/'images.npz',images=np.asarray(images));print(family,rep,[(k,round(v['material_relative_l2'],3)) for k,v in case['methods'].items()],flush=True)
 plot(records,np.asarray(images))
def plot(records,images):
 idx=[next(i for i,r in enumerate(records) if r['family']==f) for f in CFG['families']];fig,axes=plt.subplots(3,5,figsize=(14,8),sharex=True,sharey=True);vmax=float(images[idx,0].max())
 for row,i in enumerate(idx):
  for col,a in enumerate(axes[row]):
   obj=a.imshow(images[i,col].reshape(32,32).T,origin='lower',extent=[-20,20,-20,20],vmin=0,vmax=vmax,cmap='viridis');a.set_title(('truth' if col==0 else CFG['methods'][col-1]) if row==0 else '')
   if col==0:a.set_ylabel(records[i]['family']+'\ny (cm)')
   if row==2:a.set_xlabel('x (cm)')
 fig.colorbar(obj,ax=axes.ravel().tolist(),label='real contrast',shrink=.7);fig.savefig(FIG/'manifold_reconstructions.png',dpi=170,bbox_inches='tight');plt.close(fig)
if __name__=='__main__':main()
