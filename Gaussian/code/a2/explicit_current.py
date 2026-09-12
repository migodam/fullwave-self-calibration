"""Bounded explicit-current Gaussian TSOM/CSI development branch.

This is deliberately separate from the state-eliminated tangent selector.  It
uses train-receiver data to build a minimum-norm deterministic current and a
measurement-null, D-informed current basis; Gaussian material is then fitted
to a declared state-consistency residual by variable projection of weak-current
coefficients.  The current operator is not treated as an additional datum.
"""
from __future__ import annotations
import os
for _x in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_x] = "2"
import argparse, hashlib, json, platform, sys, time
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
import scipy, scipy.linalg as la
from scipy.sparse.linalg import LinearOperator, svds

ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry, VIE
from a2.ports import render
from a2.som import FREQUENCIES_HZ, SCALES, LOWER, UPPER, components, chi_and_partials, solve_stack, fit_parameter_tangent_som, FitConfig

@dataclass(frozen=True)
class Config:
    n: int=32; n_tx: int=6; n_rx: int=24; aperture: str='half'; train_stride: int=2
    noise_fraction: float=.01; seed: int=20260911; max_outer: int=30; weak_initial: int=8
    weak_cap: int=32; state_rel_expand: float=.08; state_rel_stop: float=.012
    ridge: float=1e-8; theta_trust: float=.75; gradient_stop: float=1e-5
    snapshots_every: int=5; wall_snapshot_seconds: float=30.; tangent_first: int=4; tangent_second: int=0

def realify(x): return np.r_[np.asarray(x).real.ravel(),np.asarray(x).imag.ravel()]

def deterministic_current(vie, ytrain, train):
    """jdet=S_train^*(S_train S_train^*)^-1 y, independently for each Tx."""
    st=vie.S[train]; gram=st@st.conj().T
    return st.conj().T@la.solve(gram,ytrain,assume_a='her')

def weak_basis(vie, train, rank):
    """Qweak=orth((I-Sdag S)V_D); V_D is a matrix-free D right-SVD range."""
    st=vie.S[train]; gram=st@st.conj().T
    def project(x):
        return x-st.conj().T@la.solve(gram,st@x,assume_a='her')
    n=vie.D.npix
    op=LinearOperator((n,n),matvec=lambda x:vie.D.matvec(x),rmatvec=lambda x:vie.D.rmatvec(x),dtype=complex)
    _,s,vh=svds(op,k=min(rank,n-2),which='LM',tol=2e-6,maxiter=300)
    q,_=la.qr(project(vh.conj().T),mode='economic')
    return q[:,:min(rank,q.shape[1])], np.sort(s)[::-1]

def variable_project(vies, train, observed, theta, bases, cfg):
    """For fixed theta solve each weak gamma against full state residual."""
    residuals=[]; currents=[]; gammas=[]; meta=[]
    for vie,fr,y,q in zip(vies,FREQUENCIES_HZ,observed,bases):
        chi,_=chi_and_partials(theta,vie,fr); jdet=deterministic_current(vie,y[train],train)
        r0=jdet-chi[:,None]*(vie.E+vie.D.matmat(jdet)); g=q-chi[:,None]*vie.D.matmat(q)
        normal=g.conj().T@g+cfg.ridge*np.eye(q.shape[1]); gamma=la.solve(normal,-g.conj().T@r0)
        r=r0+g@gamma; currents.append(jdet+q@gamma); gammas.append(gamma); residuals.append(r)
        meta.append({'weak_rank':q.shape[1],'condition_normal':float(np.linalg.cond(normal)),'state_relative':float(np.linalg.norm(r)/max(np.linalg.norm(chi[:,None]*vie.E),1e-30))})
    return np.concatenate([realify(x) for x in residuals]),currents,gammas,meta

def state_jacobian_fd(vies,train,observed,theta,bases,cfg,base):
    """Profiled finite differences in only 12 material coordinates; logged cost."""
    j=np.empty((base.size,theta.size)); eps=2e-4
    for p in range(theta.size):
        d=np.zeros_like(theta);d[p]=SCALES[p]*eps
        rp,*_=variable_project(vies,train,observed,np.clip(theta+d,LOWER,UPPER),bases,cfg)
        rm,*_=variable_project(vies,train,observed,np.clip(theta-d,LOWER,UPPER),bases,cfg)
        j[:,p]=(rp-rm)/(2*eps)
    return j

def explicit_fit(case,vies,train,held,observed,clean,cfg):
    theta=np.asarray(case['initial'],float).copy(); rank=cfg.weak_initial; started=time.perf_counter(); snapshots=[]; history=[]; last_snapshot=0.; cost={'basis_svd_calls':0,'D_svd_rank_sum':0,'profile_evaluations':0,'profile_fd_evaluations':0}
    termination='iteration_cap'
    for it in range(cfg.max_outer):
        bases=[]; spectra=[]
        for vie in vies:
            q,s=weak_basis(vie,train,rank);bases.append(q);spectra.append(s.tolist());cost['basis_svd_calls']+=1;cost['D_svd_rank_sum']+=rank
        r,curr,gamma,m=variable_project(vies,train,observed,theta,bases,cfg);cost['profile_evaluations']+=1
        loss=float(np.mean(r*r)); state=max(x['state_relative'] for x in m)
        if state>cfg.state_rel_expand and rank<cfg.weak_cap:
            rank=min(2*rank,cfg.weak_cap);history.append({'iteration':it,'event':'expand_weak_rank','state_relative':state,'new_weak_rank':rank});continue
        jac=state_jacobian_fd(vies,train,observed,theta,bases,cfg,r);cost['profile_fd_evaluations']+=24
        grad=jac.T@r;ginf=float(np.linalg.norm(grad,np.inf))
        if state<=cfg.state_rel_stop: termination='state_discrepancy'; break
        if ginf<=cfg.gradient_stop: termination='gradient'; break
        delta=la.solve(jac.T@jac+1e-5*np.eye(theta.size),-grad,assume_a='pos');delta*=min(1.,cfg.theta_trust/max(np.linalg.norm(delta),1e-30))
        candidate=np.clip(theta+SCALES*delta,LOWER,UPPER);rc,*_=variable_project(vies,train,observed,candidate,bases,cfg);cost['profile_evaluations']+=1
        accepted=float(np.mean(rc*rc))<loss
        if accepted:theta=candidate
        else: termination='rejected_step'
        helderr=np.sqrt(sum(np.linalg.norm(vie.S[held]@cur-clean[i]['scattered'][held])**2 for i,(vie,cur) in enumerate(zip(vies,curr))))/np.sqrt(sum(np.linalg.norm(x['scattered'][held])**2 for x in clean))
        rec={'iteration':it,'loss':loss,'state_relative':state,'gradient_inf':ginf,'weak_rank':rank,'accepted':accepted,'held_scattered_relative':float(helderr),'normal_conditions':[x['condition_normal'] for x in m]};history.append(rec)
        if it%cfg.snapshots_every==0 or time.perf_counter()-last_snapshot>=cfg.wall_snapshot_seconds:
            snapshots.append({**rec,'wall_seconds':time.perf_counter()-started});last_snapshot=time.perf_counter()
        if not accepted:break
    # final physical current and output metrics
    bases=[weak_basis(v,train,rank)[0] for v in vies];cost['basis_svd_calls']+=len(vies);cost['D_svd_rank_sum']+=len(vies)*rank
    r,curr,_,meta=variable_project(vies,train,observed,theta,bases,cfg);cost['profile_evaluations']+=1
    true=render(components(case['truth'],FREQUENCIES_HZ[0]),vies[0].points).real;got=render(components(theta,FREQUENCIES_HZ[0]),vies[0].points).real
    helderr=np.sqrt(sum(np.linalg.norm(vie.S[held]@cur-clean[i]['scattered'][held])**2 for i,(vie,cur) in enumerate(zip(vies,curr))))/np.sqrt(sum(np.linalg.norm(x['scattered'][held])**2 for x in clean))
    return {'method':'explicit_current_tsom_csi','parameters':theta.tolist(),'material_relative_l2':float(np.linalg.norm(got-true)/np.linalg.norm(true)),'held_scattered_relative':float(helderr),'state_relative':max(x['state_relative'] for x in meta),'termination':termination,'history':history,'snapshots':snapshots,'cost':cost,'final_weak_rank':rank,'wall_seconds':time.perf_counter()-started}

def cases(cfg):
    rng=np.random.default_rng(cfg.seed);out=[]
    for i in range(10):
        strong=i>=5; close=i%2==1; a=.65 if strong else .38; sep=.048 if close else .105
        truth=np.array([a,-sep/2+rng.uniform(-.012,.012),rng.uniform(-.025,.025),np.log(.022+rng.uniform(0,.012)),np.log(.024+rng.uniform(0,.014)),rng.uniform(-.4,.4),a*.82,sep/2+rng.uniform(-.012,.012),rng.uniform(-.025,.025),np.log(.021+rng.uniform(0,.014)),np.log(.023+rng.uniform(0,.016)),rng.uniform(-.4,.4)])
        initial=truth.copy();initial[::6]*=.65;initial[1::6]+=rng.normal(0,.012,2);initial[2::6]+=rng.normal(0,.012,2);initial[3::6]+=rng.normal(0,.18,2);initial[4::6]+=rng.normal(0,.18,2)
        out.append({'id':i,'family':'strong' if strong else 'moderate','geometry':'close' if close else 'separated','truth':truth.tolist(),'initial':initial.tolist(),'noise_seed':cfg.seed+100+i})
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');p.add_argument('--resume',action='store_true');a=p.parse_args();cfg=Config();out=ROOT/'runs/a2/explicit_current';out.mkdir(parents=True,exist_ok=True)
    source={'config':asdict(cfg),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'development only; explicit current state consistency is not extra observation; frozen test unopened','environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'thread_cap':2}}
    (out/'frozen_development_config.json').write_text(json.dumps(source,indent=2))
    geom=Geometry(n=cfg.n,n_tx=cfg.n_tx,n_rx=cfg.n_rx,aperture=cfg.aperture);vies=[VIE(geom,f) for f in FREQUENCIES_HZ];train=np.arange(cfg.n_rx)%cfg.train_stride==0;held=~train; records=json.loads((out/'results.json').read_text()) if a.resume and (out/'results.json').exists() else []
    for case in cases(cfg)[:1 if a.pilot else 10]:
        if any(x['id']==case['id'] for x in records):continue
        _,clean,_=solve_stack(np.array(case['truth']),vies,rtol=2e-9);rng=np.random.default_rng(case['noise_seed']);allsca=np.concatenate([x['scattered'].ravel() for x in clean]);sig=cfg.noise_fraction*np.linalg.norm(allsca)/np.sqrt(allsca.size);obs=[x['scattered']+sig/np.sqrt(2)*(rng.normal(size=x['scattered'].shape)+1j*rng.normal(size=x['scattered'].shape)) for x in clean]
        explicit=explicit_fit(case,vies,train,held,obs,clean,cfg)
        # Same full-wave data, initial point, and 30-step cap; this is baseline, not an oracle.
        lm=fit_parameter_tangent_som('ordinary_lm',np.array(case['initial']),np.array(case['truth']),obs,clean,vies,train,held,sig,FitConfig(max_iter=cfg.max_outer,first_rank=12,second_rank=0))
        adaptive=fit_parameter_tangent_som('adaptive_tsomg',np.array(case['initial']),np.array(case['truth']),obs,clean,vies,train,held,sig,FitConfig(max_iter=cfg.max_outer,first_rank=cfg.tangent_first,second_rank=cfg.tangent_second,adaptive_rank=True))
        rec={**case,'noise_sigma':float(sig),'explicit_current':explicit,'gaussian_lm':lm,'adaptive_tangent_tsomg':adaptive};records.append(rec);(out/'results.json').write_text(json.dumps(records,indent=2));print(json.dumps({'id':case['id'],'explicit':explicit['material_relative_l2'],'lm':lm['material_relative_l2'],'adaptive':adaptive['material_relative_l2'],'wall':explicit['wall_seconds']+lm['wall_seconds']+adaptive['wall_seconds']},indent=2),flush=True)
if __name__=='__main__':main()
