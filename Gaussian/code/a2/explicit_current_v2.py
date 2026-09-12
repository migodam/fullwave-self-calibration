"""Correctly composed explicit-current TSOM-style bounded diagnostic.

Not a claim about traditional TSOM: deterministic/weak splits and rank rules
are declared heuristics for this custom scalar implementation.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='2'
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,svds
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.som import FREQUENCIES_HZ,SCALES,LOWER,UPPER,chi_and_partials,solve_stack,components
from a2.ports import render

RANKS=(16,32,64,128)
def split(v,train,y):
 S=v.S[train];U,s,Vh=la.svd(S,full_matrices=False);keep=s>=.02*s[0];Ur=U[:,keep];sr=s[keep];Vr=Vh.conj().T[:,keep]
 j=Vr@((Ur.conj().T@y)/sr[:,None]);return j,lambda x:x-Vr@(Vr.conj().T@x),{'retained':int(keep.sum()),'threshold':.02}
def basis(v,P,rank):
 n=v.D.npix;op=LinearOperator((n,n),matvec=lambda x:v.D.matvec(P(x)),rmatvec=lambda x:P(v.D.rmatvec(x)),dtype=complex)
 _,_,vh=svds(op,k=min(rank,n-2),which='LM',tol=1e-6);Q,_=la.qr(P(vh.conj().T),mode='economic');return Q[:,:rank]
def profile(vies,train,obs,theta,cache,ridge=1e-8):
 vals=[];curr=[];meta=[]
 for i,(v,f,y) in enumerate(zip(vies,FREQUENCIES_HZ,obs)):
  chi,_=chi_and_partials(theta,v,f);jd,P,_=split(v,train,y[train]);Q,DQ,SQ=cache[i];r0=jd-chi[:,None]*(v.E+v.D.matmat(jd));G=Q-chi[:,None]*DQ
  # data and state blocks both normalized, equal declared weight ratio.
  nd=max(np.linalg.norm(y[train]),1e-30);ns=max(np.linalg.norm(chi[:,None]*v.E),1e-30);A=np.vstack([SQ/nd,G/ns]);beta=np.empty((Q.shape[1],jd.shape[1]),complex)
  for t in range(jd.shape[1]):
   b=np.r_[-(v.S[train]@jd-y[train])[:,t]/nd,(-r0[:,t]/ns)];beta[:,t]=la.solve(A.conj().T@A+ridge*np.eye(Q.shape[1]),A.conj().T@b)
  j=jd+Q@beta;dr=v.S[train]@j-y[train];sr=j-chi[:,None]*(v.E+v.D.matmat(j));vals.append(np.r_[dr.ravel()/nd,sr.ravel()/ns]);curr.append(j);meta.append(float(np.linalg.norm(sr)/ns))
 return np.concatenate([np.r_[z.real,z.imag] for z in vals]),curr,max(meta)
def main():
 out=ROOT/'runs/a2/explicit_current_v2';out.mkdir(parents=True,exist_ok=True);src=ROOT/'runs/a2/explicit_current/frozen_development_config.json';cases=json.loads((ROOT/'runs/a2/explicit_current/results.json').read_text()) if (ROOT/'runs/a2/explicit_current/results.json').exists() else []
 cfg={'ranks':RANKS,'deterministic_relative_singular_threshold':.02,'ridge_relative':1e-8,'state_expand':.08,'state_stop':.012,'max_outer':30,'scope':'custom explicit-current diagnostic; fixed Gaussian 12 parameters/known loss; not traditional TSOM'};(out/'frozen_config.json').write_text(json.dumps({'config':cfg,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'v1_config':str(src)},indent=2))
 geom=Geometry(n=32,n_tx=6,n_rx=24,aperture='half');vies=[VIE(geom,f) for f in FREQUENCIES_HZ];train=np.arange(24)%2==0;held=~train;records=[]
 selected=cases if len(cases)<=3 else [cases[i] for i in (0,5,7)]
 for c in selected:
  theta=np.array(c['initial']);_,clean,_=solve_stack(np.array(c['truth']),vies,rtol=2e-9);rng=np.random.default_rng(c['noise_seed']);allf=np.concatenate([x['scattered'].ravel() for x in clean]);sig=.01*np.linalg.norm(allf)/np.sqrt(allf.size);obs=[x['scattered']+sig/np.sqrt(2)*(rng.normal(size=x['scattered'].shape)+1j*rng.normal(size=x['scattered'].shape)) for x in clean];rank=16;hist=[];t0=time.perf_counter()
  for it in range(30):
   cache=[]
   for v,y in zip(vies,obs):
    jd,P,_=split(v,train,y[train]);Q=basis(v,P,rank);cache.append((Q,v.D.matmat(Q),v.S[train]@Q))
   r,cur,state=profile(vies,train,obs,theta,cache);hist.append({'iteration':it,'rank':rank,'state_relative':state,'loss':float(np.mean(r*r))})
   if state>.08 and rank<128:rank*=2;continue
   if state<=.012:break
   # Profiled finite-difference LM step in the frozen 12 coordinates.
   J=[];eps=2e-4
   for q in range(12):
    d=np.zeros(12);d[q]=SCALES[q]*eps;rp,*_=profile(vies,train,obs,np.clip(theta+d,LOWER,UPPER),cache);rm,*_=profile(vies,train,obs,np.clip(theta-d,LOWER,UPPER),cache);J.append((rp-rm)/(2*eps))
   J=np.column_stack(J);delta=la.solve(J.T@J+1e-5*np.eye(12),-J.T@r);cand=np.clip(theta+SCALES*delta,LOWER,UPPER);rc,*_=profile(vies,train,obs,cand,cache)
   if np.mean(rc*rc)<np.mean(r*r):theta=cand
   else:break
  r,cur,state=profile(vies,train,obs,theta,cache);truth=render(components(c['truth'],FREQUENCIES_HZ[0]),vies[0].points).real;got=render(components(theta,FREQUENCIES_HZ[0]),vies[0].points).real;helderr=np.sqrt(sum(np.linalg.norm(v.S[held]@j-clean[i]['scattered'][held])**2 for i,(v,j) in enumerate(zip(vies,cur))))/np.sqrt(sum(np.linalg.norm(x['scattered'][held])**2 for x in clean));records.append({'id':c['id'],'parameters':theta.tolist(),'material_relative_l2':float(np.linalg.norm(got-truth)/np.linalg.norm(truth)),'held_scattered_relative':float(helderr),'state_relative':state,'history':hist,'wall_seconds':time.perf_counter()-t0,'diagnostic_three_case':len(cases)>3});(out/'results.json').write_text(json.dumps(records,indent=2));print(json.dumps(records[-1]),flush=True)
if __name__=='__main__':main()
