"""Frozen finite-catalog operational certificate probe; not continuous imaging."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
from scipy.sparse.linalg import gmres
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.certificate import gamma_for_fft_kernel
OUT=ROOT/'runs/a2/finite_catalog';FIG=ROOT/'figures/a2';OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
AMP=(.3,.45,.6);SIG=(.035,.045,.055);STAGES=(2,4,8,16)
def g(p,a,s,c=(0,0)):
 return a*np.exp(-np.sum((p-np.asarray(c))**2,axis=1)/(2*s*s))
def pair(p,a,s,ratio):
 d=ratio*s;sx=np.sqrt(s*s-d*d);sy=s;aa=a*s*s/(2*sx*sy)
 x,y=p[:,0],p[:,1]
 return aa*np.exp(-.5*((x+d)/sx)**2-.5*(y/sy)**2)+aa*np.exp(-.5*((x-d)/sx)**2-.5*(y/sy)**2)
def catalog(p):
 old=[('old',a,s,0.,g(p,a,s)) for a in AMP for s in SIG]
 new=old+[('split',a,s,r,pair(p,a,s,r)) for a in AMP for s in SIG for r in (.3,.7)]
 return old,new
def field(v,base,exact=True,restart=None):
 chi=base*(1+.3j*2.25e9/(v.k*299792458/(2*np.pi)))
 if exact:return v.forward(chi,rtol=1e-10)['scattered'],{'solves':v.E.shape[1],'info':[0]*v.E.shape[1]},chi
 A=v.system(chi);j=np.empty_like(v.E);info=[]
 for q in range(v.E.shape[1]): j[:,q],ii=gmres(A,chi*v.E[:,q],restart=restart,maxiter=1,rtol=0,atol=0);info.append(int(ii))
 return v.S@j,{'solves':v.E.shape[1],'info':info},chi,j
def main():
 rng=np.random.default_rng(2026091106);vs=[VIE(Geometry(n=32,n_tx=4,n_rx=16),f) for f in (1.5e9,2.75e9)];old,new=catalog(vs[0].points);allc=old+new
 # Analytic pair invariants: total mass=a*2pi*s^2 and second moments
 # E[x^2]=E[y^2]=s^2 under the normalized mixture.
 for s in SIG:
  for r in (.3,.7):
   d=r*s;sx=np.sqrt(s*s-d*d);sy=s;aa=s*s/(2*sx*sy)
   assert abs(2*aa*sx*sy-s*s)<1e-14 and abs(sx*sx+d*d-s*s)<1e-14 and abs(sy*sy-s*s)<1e-14
 # Exact fields are precomputed only as fixed finite-catalog evaluation data.
 refs=[]
 for c in allc:refs.append(np.concatenate([field(v,c[4])[0].ravel() for v in vs]))
 refs=np.asarray(refs);np.savez_compressed(OUT/'reference_catalog_fields.npz',fields=refs,labels=np.asarray([str(c[:4]) for c in allc]));true=[('old_mid',old[4]),('split_03',next(c for c in new if c[1:4]==(.45,.045,.3))),('split_07',next(c for c in new if c[1:4]==(.45,.045,.7)))]
 result={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Frozen finite same-operator catalog check. Exact finite catalog references are evaluation only; no statement about continuous material classes, independent solver, or global theorem.','stages':list(STAGES),'cases':[]};M=refs.shape[1]
 for label,tcase in true:
  ti=allc.index(tcase);truth=refs[ti];sig=.01*np.sqrt(np.mean(abs(truth)**2));budget=float(sig*np.sqrt(M+2*np.sqrt(M*np.log(100))+2*np.log(100)))
  case={'truth':label,'noise_sigma_complex':float(sig),'budget':budget,'seeds':[]}
  stage_cache={}
  for seed in range(10):
   nr=np.random.default_rng(seed);y=truth+sig/np.sqrt(2)*(nr.normal(size=M)+1j*nr.normal(size=M));exact_err=np.linalg.norm(refs-y,axis=1);eo=exact_err[:len(old)].min();en=exact_err[len(old):].min();exact='supported_new_structure' if eo>budget and en<=budget else ('both_compatible' if eo<=budget and en<=budget else 'undetermined')
   row={'seed':seed,'exact_finite_catalog_decision':exact,'exact_old_min':float(eo),'exact_new_min':float(en),'stages':[]}
   for stage in STAGES:
    errs=[] if seed==0 else [float(np.linalg.norm(z-y)) for z in stage_cache[stage]['fields']];bounds=[] if seed==0 else list(stage_cache[stage]['bounds']);infos=[];cost=0
    approx_fields=[]
    for c in ([] if seed else allc):
     ff=[];bb=[]
     for v in vs:
      out=field(v,c[4],False,stage);f,meta,chi,j=out;ff.append(f.ravel());gamma,gmeta=gamma_for_fft_kernel(chi,v.D)
      if gamma is None or not gmeta['applicable']: bb.append(np.inf)
      else:
       rp=v.E-j/chi[:,None]+v.D.matmat(j); sg=np.linalg.norm(v.S/np.sqrt(gamma)[None,:],'fro'); bb.append(float(sg*np.linalg.norm(rp/np.sqrt(gamma)[:,None],'fro')))
      infos+=meta['info'];cost+=meta['solves']
     z=np.concatenate(ff);approx_fields.append(z);errs.append(float(np.linalg.norm(z-y)));bounds.append(float(np.sqrt(sum(x*x for x in bb))))
    if seed==0: stage_cache[stage]={'fields':approx_fields,'bounds':bounds,'solver_rhs':cost,'infos':infos}
    errs=np.asarray(errs);bounds=np.asarray(bounds);ol=max(0.,float(np.min(errs[:len(old)]+bounds[:len(old)]*0)-np.max(bounds[:len(old)]))) # conservative union lower using candidatewise e-b
    ol=float(np.min(np.maximum(0,errs[:len(old)]-bounds[:len(old)]))); ou=float(np.min(errs[:len(old)]+bounds[:len(old)]))
    nu=float(np.min(errs[len(old):]+bounds[len(old):]));naive='supported_new_structure' if errs[:len(old)].min()>budget and errs[len(old):].min()<=budget else ('both_compatible' if errs[:len(old)].min()<=budget and errs[len(old):].min()<=budget else 'undetermined')
    cert='supported_new_structure' if ol>budget and nu<=budget else ('both_compatible' if ou<=budget and nu<=budget else 'undetermined')
    row['stages'].append({'restart':stage,'old_lower':ol,'old_upper':ou,'new_upper':nu,'certificate_decision':cert,'naive_decision':naive,'solver_rhs_precompute':stage_cache[stage]['solver_rhs'],'nonconverged_info_count_precompute':sum(i!=0 for i in stage_cache[stage]['infos']),'mean_field_error':float(errs.mean()),'mean_raw_bound':float(bounds.mean()),'candidate_field_errors':errs.tolist(),'candidate_raw_bounds':bounds.tolist()})
   case['seeds'].append(row)
  np.savez_compressed(OUT/f'approximate_fields_{label}.npz',**{f'restart_{k}':np.asarray(v['fields']) for k,v in stage_cache.items()})
  result['cases'].append(case);(OUT/'results.json').write_text(json.dumps(result,indent=2))
 # lightweight plots after all fixed records.
 fig,ax=plt.subplots(1,2,figsize=(10,3.5));
 for c in result['cases']:
  ax[0].plot(STAGES,[np.mean([s['stages'][i]['mean_field_error'] for s in c['seeds']]) for i in range(4)],marker='o',label=c['truth']);ax[1].plot(STAGES,[np.mean([s['stages'][i]['mean_raw_bound'] for s in c['seeds']]) for i in range(4)],marker='o',label=c['truth'])
 ax[0].set(xlabel='GMRES restart, one cycle',ylabel='mean candidate field error');ax[1].set(xlabel='GMRES restart, one cycle',ylabel='mean raw bound');[a.legend() for a in ax];fig.tight_layout();fig.savefig(FIG/'catalog_error_bound.png',dpi=160)
 print(json.dumps({'cases':len(result['cases']),'seeds':30,'catalog_size':len(allc)}))
if __name__=='__main__':main()
