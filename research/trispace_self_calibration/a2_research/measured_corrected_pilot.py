"""Post-audit, bounded model-adequacy pilot, NOT antenna self-calibration.

Retains raw worker outputs. Fixes internal wave conventions and train/test
gain leakage. Point-source illumination is only a horn approximation. The
data convention follows the primary descriptor (p. 1570, exp(+i omega t)),
independently cross-checked with TRAINING INCIDENT fields, not scattering test
residuals. This is post-hoc exploratory model checking, not preregistered
performance validation. Real permittivity only; no lossy recovery claim.
"""
from pathlib import Path
import importlib.util
import json
import time
import numpy as np
from scipy.special import hankel1
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parent
path=ROOT.parents[1]/"delegated/a2_measured/model.py"
spec=importlib.util.spec_from_file_location("mie",path)
mie=importlib.util.module_from_spec(spec); spec.loader.exec_module(mie)

def correct_weights(k,u,nmax,a=mie.A_RAD):
    n=np.arange(-nmax,nmax+1)
    return (1j)**n*np.exp(-1j*n*np.arctan2(u[1],u[0]))
mie.jacobi_anger=correct_weights

def main():
    start=time.perf_counter()
    view,rec,freq,total,inc=mie.load_data()
    sc=total-inc
    train=(view%2==0); test=~train
    _,_,src,rx=mie.geometry()
    distance=np.linalg.norm(rx[rec-1]-src[view-1],axis=1)
    signs={}
    # One scalar gain per frequency; this cannot remove angle-dependent phase.
    for conjugate in (False,True):
        ratios=[]
        for fq in range(1,9):
            ix=train&(freq==fq)
            pred=hankel1(0,2*np.pi*fq*1e9/mie.C0*distance[ix])
            if conjugate: pred=pred.conj()
            gain=np.vdot(pred,inc[ix])/np.vdot(pred,pred)
            ratios.append(float(np.linalg.norm(gain*pred-inc[ix])**2/np.linalg.norm(inc[ix])**2))
        signs[str(conjugate)]=ratios
    conjugate=True  # Primary descriptor, p. 1570: exp(+i omega t); conjugate H1.
    assert min((False,True),key=lambda b:np.mean(signs[str(b)])) == conjugate
    keys,_=mie.build_index(view,rec,freq,sc)
    tr_views=np.arange(2,37,2); te_views=np.arange(1,37,2)
    # Coordinates scaled in centimetres for stable finite differences.
    def predictions(z,views,fq):
        batch=mie.forward_view_batch(fq,views,(z[0]*.01,z[1]*.01),complex(z[2]))
        pred=np.concatenate([batch[int(v)][rec[keys[(int(v),fq)]]-1] for v in views])
        return pred.conj() if conjugate else pred
    def evaluate(z,views,gains=None):
        num=den=0.; fitted={}; per={}
        for fq in range(1,9):
            pred=predictions(z,views,fq)
            obs=np.concatenate([sc[keys[(int(v),fq)]] for v in views])
            gain=np.vdot(pred,obs)/np.vdot(pred,pred) if gains is None else gains[fq]
            n=float(np.linalg.norm(gain*pred-obs)**2); d=float(np.linalg.norm(obs)**2)
            num+=n; den+=d; fitted[fq]=complex(gain); per[fq]=n/d
        return num/den,fitted,per
    runs=[]
    for z0 in ([3,0,3],[0,3,3],[-3,0,3],[0,-3,3]):
        fit=minimize(lambda z:evaluate(z,tr_views)[0],z0,method="L-BFGS-B",
                     bounds=[(-6,6),(-6,6),(1.2,6)],options={"maxiter":80,"ftol":1e-10})
        runs.append(dict(start=z0,z=fit.x.tolist(),train=float(fit.fun),nit=int(fit.nit),status=int(fit.status)))
        print(json.dumps(runs[-1]),flush=True)
    best=min(runs,key=lambda r:r["train"])
    tr,gains,pertr=evaluate(best["z"],tr_views)
    te,_,perte=evaluate(best["z"],te_views,gains)
    result=dict(status="exploratory_model_adequacy_only", convention_selection="primary descriptor p.1570 states exp(+i omega t); entire H1 model conjugated; training incident cross-check agrees",conjugated=conjugate,
                incident_train_normalized_squared_residual_by_sign=signs,
                data_split="even source views training; odd source views evaluation; gain estimated only on training",
                model="plane-wave dielectric cylinder, nominal source/receiver radii, fixed 15 mm object radius; object centre and real permittivity fit",
                fits=runs,best_center_m=(np.asarray(best["z"][:2])*.01).tolist(),best_eps_real=best["z"][2],
                training_normalized_squared_residual=tr,test_normalized_squared_residual_with_training_gains=te,
                train_per_frequency=pertr,test_per_frequency=perte,
                gains_train={k:[v.real,v.imag] for k,v in gains.items()},
                known_antenna_offset_truth=False,preregistered=False,wall_seconds=time.perf_counter()-start)
    (ROOT/"results/measured_corrected_pilot.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
