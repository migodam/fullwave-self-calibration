"""Post-hoc weak-data diagnostic; it does not replace the frozen v2 run."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[k]='2'
from pathlib import Path
import sys,json,hashlib,platform,time
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.extensions.probability import all_states,log_prior,normalized_logweights,q_probs,fit_product,fit_gaussian_logit
from a2.extensions.probability_run import centres,chi_for
OUT=ROOT/'runs/a2/extensions/probability/weak';OUT.mkdir(parents=True,exist_ok=True)
CFG={'kind':'post-hoc weak-data diagnostic; does not replace frozen v2','n_model':32,'n_data':64,'cell_integrated_data':True,'frequency_hz':1.5e9,'n_tx':4,'n_rx':24,'complex_flatten_indices':[0,48],'noise_fractions':{'high':.3,'medium':.8,'low':1.5},'per_group':10,'seed':20260914,'prior':'full-support eight-patch Ising','dictionary':'8 known Gaussian-centre Voronoi patches, binary labels; not blind Gaussian imaging'}
def geom(n): return Geometry(n=n,n_tx=4,n_rx=24,aperture='half')
def n64field(z):
    v=VIE(geom(64),CFG['frequency_hz'],cell_integrated=True)
    return v.forward(chi_for(v.points,z),rtol=1e-9)
def dist(q,z,F,y,ix):
    pred=q@F;mean=q@all_states()
    return {'brier':float(np.mean((mean-z)**2)),'log_score':float(-np.mean(z*np.log(np.clip(mean,1e-15,1))+(1-z)*np.log(np.clip(1-mean,1e-15,1)))),'entropy':float(np.mean(-mean*np.log(np.clip(mean,1e-15,1))-(1-mean)*np.log(np.clip(1-mean,1e-15,1)))),'field_relative_to_noisy_y':float(np.linalg.norm(pred-y)/np.linalg.norm(y)),'posterior_patch_mean':mean.tolist(),'map_state':all_states()[q.argmax()].tolist()}
def main():
    cache_path=ROOT/'runs/a2/extensions/probability/state_fields_n32.npz'
    src=[Path(__file__),ROOT/'code/a2/extensions/probability.py',ROOT/'code/a2/extensions/probability_run.py',ROOT/'code/a2/physics.py']
    frozen={'config':CFG,'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in src},'state_cache_sha256':hashlib.sha256(cache_path.read_bytes()).hexdigest(),'environment':{'numpy':np.__version__,'platform':platform.platform()}}
    (OUT/'frozen_config.json').write_text(json.dumps(frozen,indent=2))
    tic=time.perf_counter();cache=np.load(cache_path);states,Fall=cache['states'],cache['fields'];ix=np.asarray(CFG['complex_flatten_indices']);F=Fall[:,ix];lp=log_prior(states);lpn=normalized_logweights(lp);prior=np.exp(lpn)
    c=centres();rc=np.array([[-.085,0.],[.085,0.],[0.,.085]]);P=np.c_[np.ones(8),np.exp(-np.sum((c[:,None]-rc[None])**2,axis=2)/(2*.075**2))]
    ref=n64field(states[-1]);reference_rms=float(np.sqrt(np.mean(abs(ref['scattered'].ravel()[ix])**2)))
    rng=np.random.default_rng(CFG['seed']);records=[];truths=[];clean=[];forward_seconds=ref['wall_seconds']
    for group,frac in CFG['noise_fractions'].items():
        for rep in range(CFG['per_group']):
            idx=int(rng.choice(256,p=prior));z=states[idx];o=n64field(z);forward_seconds+=o['wall_seconds'];truth=o['scattered'].ravel()[ix];sigma=frac*reference_rms
            y=truth+sigma/np.sqrt(2)*(rng.normal(size=2)+1j*rng.normal(size=2));ll=-np.sum(abs(F-y)**2,axis=1)/(sigma*sigma);postlog=normalized_logweights(lpn+ll);post=np.exp(postlog);logev=float(np.log(np.exp(lpn+ll-(lpn+ll).max()).sum())+(lpn+ll).max())
            t=time.perf_counter();op,ep,qp=fit_product(ll,lpn,states);tp=time.perf_counter()-t;t=time.perf_counter();og,eg,qg=fit_gaussian_logit(ll,lpn,states,P);tg=time.perf_counter()-t
            row={'group':group,'rep':rep,'truth_index':idx,'truth_state':z.tolist(),'sigma':sigma,'methods':{}}
            for name,q,opt,seconds in [('exact',post,None,None),('product_vi',q_probs(qp,states),op,tp),('gaussian_rbf_logit_vi',q_probs(qg,states),og,tg)]:
                logq=np.log(np.clip(q,1e-300,1));elbo=float(q@(ll+lpn)-np.sum(q*logq));kl=float(np.sum(q*(logq-postlog)));m=dist(q,z,F,y,ix);m.update({'elbo':elbo,'kl_to_exact':kl,'log_evidence_minus_elbo':logev-elbo,'elbo_kl_identity_error':abs((logev-elbo)-kl),'optimizer_success':None if opt is None else bool(opt.success),'optimizer_seconds':seconds});row['methods'][name]=m
            records.append(row);truths.append(z);clean.append(truth)
    maxerr=max(r['methods'][m]['elbo_kl_identity_error'] for r in records for m in r['methods'])
    if maxerr>2e-8: raise RuntimeError(maxerr)
    np.savez_compressed(OUT/'n64_clean_fields_and_truth.npz',truth_indices=np.asarray([r['truth_index'] for r in records]),truth_states=np.asarray(truths),clean_held_channels=np.asarray(clean),channel_indices=ix)
    result={'scope':CFG['kind'],'records':records,'reference_rms':reference_rms,'n64_data_forward_wall_seconds':forward_seconds,'total_wall_seconds':time.perf_counter()-tic,'max_elbo_kl_identity_error':maxerr,'gaussian_rbf_design':{'intercept':True,'centres_m':rc.tolist(),'sigma_m':.075}}
    (OUT/'results.json').write_text(json.dumps(result,indent=2));print(json.dumps({'records':len(records),'n64_seconds':forward_seconds,'identity':maxerr},indent=2))
if __name__=='__main__': main()
