"""Independent v2 finite-patch Bayes check with N64 generated data and held Rx."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='2'
from pathlib import Path
import sys,json,hashlib,platform,time
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.extensions.probability import all_states,log_prior,normalized_logweights,q_probs,fit_product,fit_gaussian_logit
from a2.extensions.probability_run import centres,chi_for,patch_map
OUT=ROOT/'runs/a2/extensions/probability/v2';FIG=ROOT/'figures/a2/extensions';OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
CFG={'n_model':32,'n_data':64,'cell_integrated_data':True,'frequencies_hz':[1.5e9,2.75e9],'n_tx':4,'n_rx':24,'train_receivers':'even','held_receivers':'odd','noise_fraction':{'high':.03,'medium':.1,'low':.3},'per_group':10,'seed':20260913,'finite_dictionary':'8 fixed known Gaussian-centre Voronoi patches, binary labels; not blind Gaussian imaging','prior':'full-support Ising ring prior; empty state retained','empty_state_noise_scale':'fixed N32 all-material scattered-field RMS reference'}
def g(n,cell=False):return Geometry(n=n,n_tx=4,n_rx=24,aperture='half')
def mask():
 m=[]
 for _ in range(2):m.extend(np.repeat(np.arange(24)%2==0,4))
 return np.asarray(m)
def field_n64(z):
 o=[];wall=0.
 for f in CFG['frequencies_hz']:
  v=VIE(g(64),f,cell_integrated=True);q=v.forward(chi_for(v.points,z),rtol=1e-9);o.append(q['scattered'].ravel());wall+=q['wall_seconds']
 return np.concatenate(o),wall
def dist(q,z,F,y):
 pred=q@F;mean=q@all_states();mapz=all_states()[q.argmax()];return {'brier':float(np.mean((mean-z)**2)),'log_score':float(-np.mean(z*np.log(np.clip(mean,1e-15,1))+(1-z)*np.log(np.clip(1-mean,1e-15,1)))),'entropy':float(np.mean(-mean*np.log(np.clip(mean,1e-15,1))-(1-mean)*np.log(np.clip(1-mean,1e-15,1)))),'fit_field_relative':float(np.linalg.norm(pred[TR]-y[TR])/np.linalg.norm(y[TR])),'held_field_relative':float(np.linalg.norm(pred[HE]-y[HE])/np.linalg.norm(y[HE])),'posterior_patch_mean':mean.tolist(),'map_state':mapz.tolist()}
def main():
 cache_path=ROOT/'runs/a2/extensions/probability/state_fields_n32.npz'
 src=[ROOT/'code/a2/extensions/probability_v2.py',ROOT/'code/a2/extensions/probability.py',ROOT/'code/a2/extensions/probability_run.py',ROOT/'code/a2/physics.py']
 frozen={'config':CFG,'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in src},'state_cache_sha256':hashlib.sha256(cache_path.read_bytes()).hexdigest(),'environment':{'numpy':np.__version__,'platform':platform.platform()},'scope':'N64 cell-integrated data; N32 cached finite dictionary posterior. Held odd receivers are not used in fitting.'}
 (OUT/'frozen_config.json').write_text(json.dumps(frozen,indent=2))
 tic_all=time.perf_counter();tic_cache=time.perf_counter()
 cache=np.load(cache_path);states,F=cache['states'],cache['fields'];cache_load_seconds=time.perf_counter()-tic_cache
 lp=log_prior(states);lpn=normalized_logweights(lp);prior=np.exp(lpn);reference_rms=float(np.sqrt(np.mean(abs(F[-1])**2)));global TR,HE;TR=mask();HE=~TR
 # Actual RBF Gaussian-logit features over frozen patch centres.
 c=centres();rbf_c=np.array([[-.085,0.],[.085,0.],[0.,.085]]);P=np.c_[np.ones(8),np.exp(-np.sum((c[:,None]-rbf_c[None])**2,axis=2)/(2*.075**2))]
 rng=np.random.default_rng(CFG['seed']);records=[];examples={};n64cost=0.
 for group,frac in CFG['noise_fraction'].items():
  for rep in range(10):
   idx=int(rng.choice(256,p=prior));z=states[idx];truth,wall=field_n64(z);n64cost+=wall;sigma=frac*max(float(np.sqrt(np.mean(abs(truth)**2))),reference_rms);y=truth+sigma/np.sqrt(2)*(rng.normal(size=len(truth))+1j*rng.normal(size=len(truth)))
   ll=-np.sum(abs(F[:,TR]-y[TR])**2,axis=1)/(sigma*sigma);postlog=normalized_logweights(lpn+ll);post=np.exp(postlog);logev=float(np.log(np.exp(lpn+ll-(lpn+ll).max()).sum())+(lpn+ll).max())
   tic=time.perf_counter();op,ep,qp=fit_product(ll,lpn,states);product_seconds=time.perf_counter()-tic
   tic=time.perf_counter();og,eg,qg=fit_gaussian_logit(ll,lpn,states,P);gaussian_seconds=time.perf_counter()-tic
   qs={'exact':post,'product_vi':q_probs(qp,states),'gaussian_rbf_logit_vi':q_probs(qg,states)}
   row={'group':group,'rep':rep,'truth_index':idx,'truth_state':z.tolist(),'sigma':float(sigma),'sigma_reference_rms':reference_rms,'methods':{}}
   for name,q in qs.items():
    elbo=float(q@(ll+lpn)-np.sum(q*np.log(np.clip(q,1e-300,1))))
    logq=np.log(np.clip(q,1e-300,1))
    kl=float(np.sum(q*(logq-postlog)))
    m=dist(q,z,F,y)
    opt=None if name=='exact' else (op if name=='product_vi' else og)
    m.update({'elbo':elbo,'log_evidence_minus_elbo':float(logev-elbo),'kl_to_exact':kl,'elbo_kl_identity_error':float(abs((logev-elbo)-kl)),'optimizer_success':None if opt is None else bool(opt.success),'optimizer_message':None if opt is None else str(opt.message),'optimizer_seconds':0. if name=='exact' else (product_seconds if name=='product_vi' else gaussian_seconds)})
    row['methods'][name]=m
   records.append(row);examples.setdefault(group,row)
 result={'cache_cost_n32':{'status':'reused immutable 256-state pilot cache; v2 recomputes full-support prior','cache_load_seconds':cache_load_seconds,'state_cache_sha256':frozen['state_cache_sha256']},'n64_data_forward_wall_seconds':n64cost,'total_wall_seconds':time.perf_counter()-tic_all,'records':records,'model_grid_vs_data_grid':'Posterior likelihood uses the frozen N32 point-kernel state cache while data are newly generated at N64 with cell-integrated sources, receivers, and FFT Green kernel. Odd Rx are independent of fit channels.','gaussian_logit_design':{'intercept':True,'rbf_centres_m':rbf_c.tolist(),'rbf_sigma_m':.075,'features':'[1, exp(-||centre-c_l||^2/(2 sigma^2)) for three fixed RBF centres]'},'noise_scale_reference_rms':reference_rms}
 identity=max(r['methods'][n]['elbo_kl_identity_error'] for r in records for n in r['methods'])
 if identity > 2e-8: raise RuntimeError(f'ELBO/KL identity check failed: {identity}')
 result['max_elbo_kl_identity_error']=identity
 (OUT/'results.json').write_text(json.dumps(result,indent=2))
 # Three rows (noise groups), four comparable material maps with common scale.
 fig,ax=plt.subplots(3,4,figsize=(10,7));v=VIE(g(32),CFG['frequencies_hz'][0]);vmax=.55
 for i,(group,row) in enumerate(examples.items()):
  maps=[chi_for(v.points,row['truth_state']).real]
  for name in ('exact','product_vi','gaussian_rbf_logit_vi'):
   mean=np.asarray(row['methods'][name]['posterior_patch_mean']);ids,inside=patch_map(v.points);maps.append(np.where(inside,mean[ids]*.55,0).reshape(32,32))
  for j,(im,title) in enumerate(zip(maps,['truth','exact','product','Gaussian RBF'])):
   h=ax[i,j].imshow(im.reshape(32,32).T if im.ndim==1 else im.T,origin='lower',vmin=0,vmax=vmax,cmap='viridis');ax[i,j].set(title=title if i==0 else '',ylabel=group if j==0 else '',xticks=[],yticks=[])
 fig.colorbar(h,ax=ax.ravel().tolist(),label='posterior/material contrast, common scale');fig.tight_layout();fig.savefig(FIG/'probability_v2_reconstructions.png',dpi=170)
 print(json.dumps({'records':len(records),'n64_wall':n64cost},indent=2))
if __name__=='__main__':main()
