"""Post-hoc RBF-VI initialization diagnosis for frozen probability_v2 data."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[k]='2'
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.extensions.probability import all_states,log_prior,normalized_logweights,q_probs,elbo_product,fit_product
from a2.extensions.probability_run import centres,chi_for
OUT=ROOT/'runs/a2/extensions/probability/v2/optimization_posthoc';OUT.mkdir(parents=True,exist_ok=True)
CFG={'scope':'post-hoc initialization diagnosis; does not replace frozen v2','seed_v2_reconstruction':20260913,'random_start_seed':20260915,'n_data':64,'cell_integrated_data':True,'frequencies_hz':[1.5e9,2.75e9],'channels':'all original v2 even-Rx training channels','starts':['zero','product-marginal logit projected to RBF design','random_0','random_1'],'selection':'maximum ELBO; product-marginal start is data-derived and is not truth-oracle'}
def geom(n): return Geometry(n=n,n_tx=4,n_rx=24,aperture='half')
def receiver_mask(): return np.concatenate([np.repeat(np.arange(24)%2==0,4) for _ in range(2)])
def n64(z):
    f=[];wall=0.
    for hz in CFG['frequencies_hz']:
        v=VIE(geom(64),hz,cell_integrated=True);o=v.forward(chi_for(v.points,z),rtol=1e-9)
        f.append(o['scattered'].ravel());wall+=o['wall_seconds']
    return np.concatenate(f),wall
def fit(loglike,lp,states,P,x0):
    t=time.perf_counter();o=minimize(lambda x:-elbo_product(loglike,lp,states,P@x)[0],x0,method='L-BFGS-B',options={'maxiter':200,'ftol':1e-12});e,q=elbo_product(loglike,lp,states,P@o.x)
    return o,e,q,time.perf_counter()-t
def main():
    cache_path=ROOT/'runs/a2/extensions/probability/state_fields_n32.npz';src=[Path(__file__),ROOT/'code/a2/extensions/probability.py',ROOT/'code/a2/extensions/probability_run.py',ROOT/'code/a2/physics.py']
    (OUT/'frozen_config.json').write_text(json.dumps({'config':CFG,'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in src},'cache_sha256':hashlib.sha256(cache_path.read_bytes()).hexdigest()},indent=2))
    C=np.load(cache_path);states,F=C['states'],C['fields'];lpn=normalized_logweights(log_prior(states));prior=np.exp(lpn);ref=float(np.sqrt(np.mean(abs(F[-1])**2)));tr=receiver_mask()
    c=centres();rc=np.array([[-.085,0.],[.085,0.],[0.,.085]]);P=np.c_[np.ones(8),np.exp(-np.sum((c[:,None]-rc[None])**2,axis=2)/(2*.075**2))]
    rng=np.random.default_rng(CFG['seed_v2_reconstruction']);rs=np.random.default_rng(CFG['random_start_seed']);records=[];idxs=[];zs=[];cleans=[];ys=[];sigmas=[];wall=0.
    for group,frac in {'high':.03,'medium':.1,'low':.3}.items():
        for rep in range(10):
            idx=int(rng.choice(256,p=prior));z=states[idx];clean,w=n64(z);wall+=w;sigma=frac*max(float(np.sqrt(np.mean(abs(clean)**2))),ref);y=clean+sigma/np.sqrt(2)*(rng.normal(size=len(clean))+1j*rng.normal(size=len(clean)))
            ll=-np.sum(abs(F[:,tr]-y[tr])**2,axis=1)/(sigma*sigma);postlog=normalized_logweights(lpn+ll);op,ep,qp=fit_product(ll,lpn,states)
            logits=np.log(np.clip(qp,1e-15,1-1e-15)/np.clip(1-qp,1e-15,1));data_x0=np.linalg.lstsq(P,np.clip(logits,-6,6),rcond=None)[0]
            candidates=[]
            for label,x0 in [('zero',np.zeros(P.shape[1])),('product_projected',data_x0),('random_0',rs.normal(size=P.shape[1])),('random_1',rs.normal(size=P.shape[1]))]:
                o,e,q,sec=fit(ll,lpn,states,P,x0);qq=q_probs(q,states);logq=np.log(np.clip(qq,1e-300,1));kl=float(np.sum(qq*(logq-postlog)));candidates.append({'start':label,'elbo':e,'kl_to_exact':kl,'brier':float(np.mean((q-z)**2)),'success':bool(o.success),'seconds':sec,'theta':o.x.tolist()})
            best=max(candidates,key=lambda a:a['elbo']);records.append({'group':group,'rep':rep,'truth_index':idx,'sigma':sigma,'product_elbo':ep,'product_marginals':qp.tolist(),'candidates':candidates,'selected':best})
            idxs.append(idx);zs.append(z);cleans.append(clean);ys.append(y);sigmas.append(sigma)
    np.savez_compressed(OUT/'v2_reconstructed_data.npz',truth_indices=np.asarray(idxs),truth_states=np.asarray(zs),clean_fields=np.asarray(cleans),noisy_fields=np.asarray(ys),sigmas=np.asarray(sigmas))
    result={'scope':CFG['scope'],'n64_reconstruction_wall_seconds':wall,'records':records,'selection_note':CFG['selection']};(OUT/'results.json').write_text(json.dumps(result,indent=2));print(json.dumps({'n':len(records),'wall':wall},indent=2))
if __name__=='__main__':main()
