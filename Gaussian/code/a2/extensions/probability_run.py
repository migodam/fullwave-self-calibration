"""Frozen finite-patch full-wave Bayesian imaging experiment."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='2'
from pathlib import Path
import sys,json,hashlib,platform,time
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.extensions.probability import *
OUT=ROOT/'runs/a2/extensions/probability';FIG=ROOT/'figures/a2/extensions';OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
CFG={'patches':8,'n_development':32,'n_validation':64,'frequencies_hz':[1.5e9,2.75e9],'n_tx':4,'n_rx':24,'aperture':'half','roi_radius_m':.14,'chi_material':[0.,.55+.11j],'ising_eta':.65,'ising_field':-.12,'condition_nonempty_prior':True,'snr_fraction':{'high':.003,'medium':.01,'low':.03},'noise_per_group':10,'seed':20260912}
def geom(n):return Geometry(n=n,n_tx=4,n_rx=24,aperture='half')
def centres():
 a=np.linspace(0,2*np.pi,8,endpoint=False);return .085*np.c_[np.cos(a),np.sin(a)]
def patch_map(points):
 c=centres();d=((points[:,None]-c[None])**2).sum(2);ids=d.argmin(1);inside=(points**2).sum(1)<=CFG['roi_radius_m']**2;return ids,inside
def chi_for(points,state):
 ids,inside=patch_map(points);return np.where(inside,np.asarray(state)[ids]*CFG['chi_material'][1],0j)
def cache_fields(n=32):
 states=all_states();fields=[];tic=time.perf_counter();solve_seconds=[]
 for z in states:
  blocks=[]
  for f in CFG['frequencies_hz']:
   v=VIE(geom(n),f);o=v.forward(chi_for(v.points,z),rtol=1e-9);blocks.append(o['scattered'].ravel());solve_seconds.append(o['wall_seconds'])
  fields.append(np.concatenate(blocks))
 return states,np.asarray(fields),{'wall_seconds':time.perf_counter()-tic,'state_solves':512,'individual_forward_wall_seconds':solve_seconds}
def metrics(q,z,F,y,sigma):
 pred=q@F;mean=(q@all_states());mapz=all_states()[q.argmax()];return {'brier':float(np.mean((mean-z)**2)),'log_score':float(-np.mean(z*np.log(np.clip(mean,1e-15,1))+(1-z)*np.log(np.clip(1-mean,1e-15,1)))),'entropy':float(np.mean(-mean*np.log(np.clip(mean,1e-15,1))-(1-mean)*np.log(np.clip(1-mean,1e-15,1)))),'posterior_predictive_relative':float(np.linalg.norm(pred-y)/np.linalg.norm(y)),'posterior_patch_mean':mean.tolist(),'map_state':mapz.tolist()}
def main():
 source={'config':CFG,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dependency':{'numpy':np.__version__,'platform':platform.platform()},'scope':'Finite 8-patch epistemic label posterior. Every label realization is solved full-wave; no F(Echi) substitution.'};(OUT/'frozen_config.json').write_text(json.dumps(source,indent=2,default=str))
 states,F,cost=cache_fields();lp=log_prior(states,CFG['ising_eta'],CFG['ising_field']);lp[0]=-1e6;lpn=normalized_logweights(lp);prior=np.exp(lpn);np.savez_compressed(OUT/'state_fields_n32.npz',states=states,fields=F,log_prior=lp)
 # fixed Gaussian-logit design: intercept + x/y/radial functions over the eight known centres.
 c=centres();P=np.c_[np.ones(8),c[:,0]/.085,c[:,1]/.085,(c[:,0]**2-c[:,1]**2)/.085**2]
 rng=np.random.default_rng(CFG['seed']);records=[];reliability={'exact':[],'product_vi':[],'gaussian_logit_vi':[]};examples={}
 for group,frac in CFG['snr_fraction'].items():
  for rep in range(10):
   idx=int(rng.choice(256,p=prior));z=states[idx];truth=F[idx];sigma=frac*np.sqrt(np.mean(abs(truth)**2));y=truth+sigma/np.sqrt(2)*(rng.normal(size=len(truth))+1j*rng.normal(size=len(truth)));ll=-np.sum(abs(F-y[None])**2,axis=1)/(sigma*sigma);post=np.exp(normalized_logweights(lp+ll));logev=float(np.log(np.exp(lp+ll-(lp+ll).max()).sum())+(lp+ll).max()-np.log(np.exp(lp).sum()))
   op,ep,qp=fit_product(ll,lpn,states);og,eg,qg=fit_gaussian_logit(ll,lpn,states,P);qp_state=q_probs(qp,states);qg_state=q_probs(qg,states);exact_elbo=float(post@(ll+lpn)-np.sum(post*np.log(np.clip(post,1e-300,1))))
   # ELBO identity: log evidence = ELBO + KL(q||posterior).
   row={'group':group,'rep':rep,'truth_state':z.tolist(),'truth_index':idx,'sigma':float(sigma),'exact_log_evidence':logev,'exact_elbo':exact_elbo,'methods':{}}
   for name,q,e,o in [('exact',post,exact_elbo,None),('product_vi',qp_state,ep,op),('gaussian_logit_vi',qg_state,eg,og)]:
    m=metrics(q,z,F,y,sigma);m.update({'elbo':float(e),'kl_to_exact':float(logev-e),'optimizer_success':None if o is None else bool(o.success),'optimizer_message':None if o is None else str(o.message)});row['methods'][name]=m;reliability[name].extend(zip(m['posterior_patch_mean'],z.tolist()))
   records.append(row);examples.setdefault(group,row)
 # N64 different-grid field validation for the three chosen physical states.
 val=[]
 for group,row in examples.items():
  z=np.asarray(row['truth_state']);blocks=[]
  for f in CFG['frequencies_hz']:v=VIE(geom(64),f);blocks.append(v.forward(chi_for(v.points,z),rtol=1e-9)['scattered'].ravel())
  f64=np.concatenate(blocks);f32=F[int(row['truth_index'])];val.append({'group':group,'truth_state':z.tolist(),'n64_field_norm':float(np.linalg.norm(f64)),'n32_to_n64_scattered_relative':float(np.linalg.norm(f64-f32)/np.linalg.norm(f64))})
 result={'cache_cost':cost,'records':records,'n64_different_grid_validation':val,'known_patch_geometry':'Eight fixed Gaussian-centre Voronoi regions inside declared ROI; labels are not arbitrary pixels or blind geometry.','scope':'Prior-predictive calibration only; not repeated-noise calibration conditional on a fixed truth, generalization, or a novelty claim.'};(OUT/'results.json').write_text(json.dumps(result,indent=2))
 # Three physical maps and reliability curves.
 fig,ax=plt.subplots(1,3,figsize=(9,3));
 for a,(g,r) in zip(ax,examples.items()):
  v=VIE(geom(32),CFG['frequencies_hz'][0]);a.imshow(chi_for(v.points,r['truth_state']).real.reshape(32,32).T,origin='lower',cmap='viridis');a.set(title=g,xticks=[],yticks=[])
 fig.tight_layout();fig.savefig(FIG/'probability_true_positions.png',dpi=170)
 fig,ax=plt.subplots(figsize=(5,4));bins=np.linspace(0,1,6)
 for name,pairs in reliability.items():
  p,z=np.asarray(pairs).T;xs=[];ys=[]
  for lo,hi in zip(bins[:-1],bins[1:]):
   m=(p>=lo)&(p<hi if hi<1 else p<=hi)
   if m.any():xs.append(p[m].mean());ys.append(z[m].mean())
  ax.plot(xs,ys,'o-',label=name)
 ax.plot([0,1],[0,1],'k--');ax.set(xlabel='posterior probability bin mean',ylabel='prior-predictive empirical frequency');ax.legend();fig.tight_layout();fig.savefig(FIG/'probability_reliability.png',dpi=170)
 print(json.dumps({'records':len(records),'cache_seconds':cost['wall_seconds'],'validation':val},indent=2))
if __name__=='__main__':main()
