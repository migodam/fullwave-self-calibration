"""Conditional finite-bank nonlinear gain-graph mechanism, not a recovery algorithm."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares
from maxwell3d import treams_field, receivers
from gain_graph import incidence

OUT=Path(__file__).resolve().parent
GRAPHS={
    'tree':[(r,0) for r in range(6)]+[(0,t) for t in range(1,4)],
    'three_cycles':[(0,t) for t in range(4)]+[(1,t) for t in range(3)]+[(2,0),(2,3)]}


class GainFit:
    def __init__(self,field,data,edges,sigma):
        self.h=np.array([field[r,t] for r,t in edges])
        self.y=np.array([data[r,t] for r,t in edges]);self.sigma=sigma
        self.vertices=sorted({('r',r) for r,t in edges}|{('t',t) for r,t in edges})
        self.anchor=('t',0);self.free=[v for v in self.vertices if v!=self.anchor]
        self.lookup={v:i for i,v in enumerate(self.free)};self.edges=edges
    def pack(self,g):
        v=np.array([g[x] for x in self.free]);return np.r_[v.real,v.imag]
    def field_jac(self,z):
        n=len(self.free);vals=z[:n]+1j*z[n:]
        g={self.anchor:1.+0j};g.update(zip(self.free,vals))
        pred=[];jac=np.zeros((len(self.edges),n),complex)
        for i,(r,t) in enumerate(self.edges):
            rv,tv=('r',r),('t',t)
            pred.append(self.h[i]*g[rv]*g[tv])
            for vertex,other in [(rv,tv),(tv,rv)]:
                if vertex!=self.anchor:jac[i,self.lookup[vertex]]=self.h[i]*g[other]
        jc=np.c_[jac,1j*jac]
        return np.array(pred),jc
    def residual(self,z):
        pred,_=self.field_jac(z);d=pred-self.y
        return np.sqrt(2)/self.sigma*np.r_[d.real,d.imag]
    def jac(self,z):
        _,j=self.field_jac(z);return np.sqrt(2)/self.sigma*np.r_[j.real,j.imag]
    def tree_start(self):
        ratios=self.y/self.h
        if np.any(np.abs(self.h)<1e-14) or np.any(np.abs(ratios)<1e-14):
            raise ValueError('nonzero edge field/data required for constructive start')
        g={self.anchor:1.+0j}
        for _ in self.vertices:
            for i,(r,t) in enumerate(self.edges):
                rv,tv=('r',r),('t',t)
                if tv in g and rv not in g:g[rv]=ratios[i]/g[tv]
                if rv in g and tv not in g:g[tv]=ratios[i]/g[rv]
        if len(g)!=len(self.vertices):raise ValueError('disconnected graph needs more anchors')
        return self.pack(g)


def check():
    rng=np.random.default_rng(8399)
    h=rng.normal(size=(6,4))+1j*rng.normal(size=(6,4))
    y=rng.normal(size=(6,4))+1j*rng.normal(size=(6,4))
    checks=[]
    for name,edges in GRAPHS.items():
        fit=GainFit(h,y,edges,.1);z=fit.tree_start()
        d=rng.normal(size=z.size);step=1e-6
        fd=(fit.residual(z+step*d)-fit.residual(z-step*d))/(2*step)
        err=np.linalg.norm(fd-fit.jac(z)@d)/np.linalg.norm(fd)
        assert err<1e-7
        tree_error=np.linalg.norm(fit.field_jac(z)[0]-fit.y)/np.linalg.norm(fit.y)
        if name=='tree':assert tree_error<1e-13
        checks.append(dict(graph=name,jacobian_relative_error=float(err),constructive_residual=float(tree_error)))
    (OUT/'results/nonlinear_gain_checks.json').write_text(json.dumps(dict(passed=True,checks=checks),indent=2)+'\n')
    print(json.dumps(checks),flush=True)


def run():
    start=time.perf_counter()
    centers=[[-.2,0,0],[.2,.05,0]];radii=[.12,.10]
    true_pose=np.array([.03,-.025,.02]);true_eps=np.array([2.4,3.1])
    bank=[('oracle_truth',true_pose,true_eps),
          ('pose_plus',true_pose+np.array([.12,0,0]),true_eps),
          ('pose_minus',true_pose-np.array([.12,0,0]),true_eps),
          ('material_plus',true_pose,true_eps+np.array([.5,-.4])),
          ('mixed',true_pose+np.array([0,.08,-.06]),true_eps+np.array([-.3,.4]))]
    fields=[]
    for name,x,eps in bank:
        f,_=treams_field(centers,radii,eps+.02j,9.,receivers(6)+x,lmax=4)
        fields.append(np.einsum('rct,c->rt',f,np.array([1.,2.,3.])/np.sqrt(14)))
    field_seconds=time.perf_counter()-start
    rows=[]
    for seed in [8301,8302]:
        rng=np.random.default_rng(seed)
        gr=np.exp(rng.normal(0,.08,6)+1j*rng.normal(0,.2,6))
        gt=np.exp(rng.normal(0,.08,4)+1j*rng.normal(0,.2,4))
        mu=gr[:,None]*fields[0]*gt[None,:]
        sigma=np.linalg.norm(mu)/np.sqrt(mu.size)*10**(-30/20)
        y=mu+sigma/np.sqrt(2)*(rng.normal(size=mu.shape)+1j*rng.normal(size=mu.shape))
        for graph,edges in GRAPHS.items():
            cases=[];t0=time.perf_counter()
            for i,(name,x,eps) in enumerate(bank):
                fit=GainFit(fields[i],y,edges,sigma);n=len(fit.free)
                starts=[fit.tree_start(),np.r_[np.ones(n),np.zeros(n)]]
                srng=np.random.default_rng(8350+i)
                starts.append(starts[1]+srng.normal(0,.1,2*n))
                trials=[]
                for z in starts:
                    t=time.perf_counter()
                    opt=least_squares(fit.residual,z,jac=fit.jac,max_nfev=200,
                        ftol=1e-11,xtol=1e-11,gtol=1e-9)
                    trials.append(dict(loss=float(opt.fun@opt.fun/2),status=int(opt.status),
                        nfev=opt.nfev,seconds=time.perf_counter()-t))
                cases.append(dict(candidate=name,pose=x.tolist(),eps=eps.tolist(),
                    best_loss=min(v['loss'] for v in trials),starts=trials))
            scores=np.array([c['best_loss'] for c in cases])
            rows.append(dict(seed=seed,graph=graph,edges=edges,
                cycles=len(edges)-int(np.linalg.matrix_rank(incidence(edges))),
                candidate_scores=cases,lowest_loss_candidate=bank[int(scores.argmin())][0],
                wrong_bank_best_loss=float(scores[1:].min()),oracle_candidate_loss=float(scores[0]),
                seconds=time.perf_counter()-t0,
                scope='oracle covered finite bank, gain-profiled training residual; no independent validation or global coverage guarantee'))
            print(json.dumps({k:rows[-1][k] for k in ['seed','graph','cycles','wrong_bank_best_loss','oracle_candidate_loss','seconds']}),flush=True)
    (OUT/'results/nonlinear_gain_graph.json').write_text(json.dumps(dict(field_generation_seconds=field_seconds,rows=rows),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--test',action='store_true');args=p.parse_args()
    if args.test:check()
    else:run()
