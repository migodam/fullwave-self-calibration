"""Exploratory algorithmic translation: residual enrichment, no truth decisions.

Not part of the frozen E4 method. Fixed seed set and rank cap declared here.
All exact states/Jacobians are evaluated only after recording certificate tests.
"""
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import qr,svd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/'research/delegated/a2_solver_audited'),str(ROOT/'research/delegated/a2_physics'),str(HERE)]
from physics import Config,Model
from common import CostLedger
from reduced import build_sensing_basis,reduced_evaluate
from cert import reduced_certification

def orth_append(U,Z,k):
    Z=Z-U@(U.conj().T@Z)
    v,s,_=svd(Z,full_matrices=False,check_finite=False)
    q=min(k,np.count_nonzero(s>max(1e-12,1e-10*s[0]))) if len(s) else 0
    return qr(np.c_[U,v[:,:q]],mode='economic')[0]
def run():
    records=[];start=time.perf_counter()
    for ap in ['full','limited']:
        m=Model(Config(N=16,aperture=ap))
        for seed in [61,62]:
            rng=np.random.default_rng(seed);a=rng.uniform(.15,.5,9) if seed==61 else rng.uniform(.7,1.4,9)
            x=np.array([.04,-.02,.03]);sigma=.002
            # Observations used only in local gradient test. This is a
            # declared reference-state test, not blind inversion evidence.
            y=m.forward(a,np.zeros(3))['total']
            for fi in [0,3]:
                l=CostLedger();initial=build_sensing_basis(m,x,(fi,),16,l,tsom_append=True,tsom_target_fraction=.5)['U']
                for policy in ['sensing_domain','residual_enriched']:
                    U=initial.copy();history=[];t0=time.perf_counter()
                    for iteration in range(9):
                        out=reduced_evaluate(m,a,x,U,(fi,),l,True)
                        c=reduced_certification(m,a,y[fi*72:(fi+1)*72],sigma,(fi,),out,{})
                        history.append({'rank':U.shape[1],**c})
                        if c['pass'] or policy=='sensing_domain' or U.shape[1]>=128:break
                        # Both state and tangent residuals are computable;
                        # their singular directions do not require true j.
                        Z=np.column_stack([np.c_[s['state_residual'],s['derivative_residual']] for s in out['states']])
                        new=orth_append(U,Z,min(16,128-U.shape[1]))
                        if new.shape[1]==U.shape[1]:break
                        U=new
                    # Evaluation-only full physical truth, never gate input.
                    fw=m.forward(a,x,[fi],jacobian=True)
                    er=np.sqrt(2)/sigma*np.linalg.norm(out['total']-fw['total'])
                    AB=np.c_[fw['A'],fw['B']];J=np.sqrt(2)/sigma*np.vstack([AB.real,AB.imag]);J[:,-3:]/=[2,2,3]
                    res=np.sqrt(2)/sigma*np.r_[(fw['total']-y[fi*72:(fi+1)*72]).real,(fw['total']-y[fi*72:(fi+1)*72]).imag]
                    ABt=np.c_[out['A'],out['B']];Jt=np.sqrt(2)/sigma*np.vstack([ABt.real,ABt.imag]);Jt[:,-3:]/=[2,2,3]
                    rt=np.sqrt(2)/sigma*np.r_[(out['total']-y[fi*72:(fi+1)*72]).real,(out['total']-y[fi*72:(fi+1)*72]).imag]
                    ge=np.linalg.norm(J.T@res-Jt.T@rt)
                    records.append({'aperture':ap,'seed':seed,'freq':fi,'policy':policy,'history':history,'actual_data_error':float(er),'actual_gradient_error':float(ge),'state_bound_valid':bool(er<=c['data_error_bound']*(1+1e-9)),'gradient_bound_valid':bool(ge<=c['gradient_error_bound']*(1+1e-9)),'wall_including_exact_diagnostic':time.perf_counter()-t0,'r_free':0})
                    print(ap,seed,fi,policy,'rank',U.shape[1],'admitted',c['pass'],flush=True)
    result={'status':'exploratory_model_conditional','runs':len(records),'records':records,'wall_seconds':time.perf_counter()-start,'decisions':'only computable residual/derivative certificates; exact values evaluation-only'}
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results/enrichment_probe.json').write_text(json.dumps(result,indent=2)+'\n')
    assert all(r['state_bound_valid'] and r['gradient_bound_valid'] for r in records)
if __name__=='__main__':run()
