"""Registered shape/material stress experiment. Default: derivative checks only."""
import argparse
import fcntl
import gc
import hashlib
import json
import os
from pathlib import Path
import sys
import time

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parent
DELEGATED=ROOT.parents[1]/'delegated'
sys.path.insert(0,str(DELEGATED/'a3_discrepancy_weighting'))
import discrepancy_weighting as dw
mf=dw.mf
OUT=ROOT/'results/shape_material_stress'
CASES=((8401,'ellipsoid',1.8),(8402,'ellipsoid',3.5),(8403,'box',1.8),(8404,'box',3.5))

class ShapeModel(mf.ExplicitFrequencyModel):
    def __init__(self,shape,n,frequencies,rx=None):
        self.frequencies=tuple(float(k) for k in frequencies)
        self.rx=mf.receivers() if rx is None else np.asarray(rx)
        self.models={}
        self.setup_seconds={}
        for k in dict.fromkeys(self.frequencies):
            start=time.perf_counter()
            self.models[k]=mf.ShapeVIE(shape,n,k)
            self.setup_seconds[str(k)]=time.perf_counter()-start
        self.cache_material=None
        self.currents={}
        self.calls=0
        self.field_seconds={str(k):0. for k in self.models}

def provenance():
    paths=[Path(__file__),ROOT/'SHAPE_MATERIAL_STRESS_PROTOCOL.md',
           Path(dw.__file__),Path(mf.__file__),ROOT/'nonspherical3d.py',
           ROOT/'nonspherical_calibration.py',ROOT/'maxwell3d.py',
           ROOT/'external/adda/src/seq/adda',
           DELEGATED/'a3_maxwell_refine/tangent_fft.py',
           DELEGATED/'a3_maxwell_fft/maxwell_fft.py']
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def test():
    OUT.mkdir(parents=True,exist_ok=True)
    hashes=provenance()
    z=np.array([2.3,.03,-.02,.04,.03]+[.01]*4+[.02]*4)
    errors={}
    for shape in ('ellipsoid','box'):
        model=ShapeModel(shape,8,mf.CHOICES['low_high'])
        y,j=model.field_jac(z)
        errs=[]
        for c in range(13):
            dz=np.eye(13)[c]*1e-5
            fd=(model.field_jac(z+dz)[0]-model.field_jac(z-dz)[0])/2e-5
            errs.append(float(np.linalg.norm(fd-j[...,c])/np.linalg.norm(fd)))
        assert max(errs)<1e-7
        errors[shape]=errs
        if shape=='ellipsoid':
            old=mf.ExplicitFrequencyModel(8,mf.CHOICES['low_high'])
            oy,oj=old.field_jac(z)
            assert np.array_equal(y,oy) and np.array_equal(j,oj)
    tagged=np.arange(576).reshape(4,12,3,4)*(1+2j)
    assert np.array_equal(mf.real(tagged)[dw.high_rows(144)],mf.real(tagged[3]))
    d,di=dw.mode_vectors(tagged[3],1.)
    for method in ('warm_raw','isotropic','rank1'):
        w,c=dw.whitening(d,di,method)
        assert np.linalg.norm(w@c@w.T-np.eye(288))/np.sqrt(288)<1e-8
    assert provenance()==hashes
    result=dict(passed=True,errors=errors,ellipsoid_bit_identical=True,
                high_block_ordering_passed=True,provenance=hashes)
    mf.write_json(OUT/'checks.json',result)
    print(json.dumps(dict(checks_passed=True,max_derivative_error=max(max(v) for v in errors.values()))),flush=True)

def npz_write(path,arrays):
    tmp=path.with_name(path.name+'.tmp')
    with tmp.open('wb') as stream:
        np.savez_compressed(stream,**arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp,path)

def case_data(seed,shape,epsilon,hashes):
    path=OUT/f'data_{seed}.npz'
    meta_path=OUT/f'data_{seed}.json'
    if meta_path.exists():
        meta=json.loads(meta_path.read_text())
        assert meta['provenance']==hashes
        if meta['status']=='ok':
            assert meta['data_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(path) as archive:
                return {k:archive[k].copy() for k in archive.files},meta
        raise RuntimeError('Incomplete or failed generation preserved; inspect before resume')
    meta=dict(seed=seed,shape=shape,epsilon=epsilon,status='started',provenance=hashes)
    mf.write_json(meta_path,meta)
    start=time.perf_counter()
    try:
        z,rng=mf.scene(seed)
        z[0]=epsilon
        rx=mf.receivers()
        held=mf.receivers(17,1.6)
        raw=[];heldraw=[];sources=[]
        for k in mf.CHOICES['low_high']:
            p,q,record=mf.reference(shape,64,k,epsilon+1j*mf.LOSS)
            raw.append(mf.independent_radiation(p,q,rx+z[1:4],k))
            heldraw.append(mf.independent_radiation(p,q,held+z[1:4],k))
            sources.append(record)
        raw,heldraw=np.array(raw),np.array(heldraw)
        factor=mf.electronics(z,mf.CHOICES['low_high'])[:,None,None,:]
        mean=raw*factor
        sigma=float(np.linalg.norm(mean[:3])/np.sqrt(mean[:3].size)*10**(-30/20))
        noise=(rng.normal(size=mean.shape)+1j*rng.normal(size=mean.shape))/np.sqrt(2)
        y=mean+sigma*noise
        ref_rng=np.random.default_rng(seed+10000)
        reference=mf.electronics(z,mf.LOW)+mf.REF_SIGMA/np.sqrt(2)*(ref_rng.normal(size=(3,4))+1j*ref_rng.normal(size=(3,4)))
        # Independent reference refinement is diagnostic and never changes y.
        ref_start=time.perf_counter()
        p96,q96,record96=mf.reference(shape,96,18.,epsilon+1j*mf.LOSS)
        refined=mf.independent_radiation(p96,q96,rx+z[1:4],18.)
        heldrefined=mf.independent_radiation(p96,q96,held+z[1:4],18.)
        sensitivity=float(np.linalg.norm(raw[3]-refined)/np.linalg.norm(refined))
        held_sensitivity=float(np.linalg.norm(heldraw[3]-heldrefined)/np.linalg.norm(heldrefined))
        arrays=dict(true=z,sigma=np.array(sigma),low=y[:3],low_high=y,reference=reference,
                    training_mean=mean,held_sensor_mean=heldraw[:3]*factor[:3],
                    held_structural_mean=heldraw[:3],noise=noise)
        assert np.array_equal(arrays['low'],arrays['low_high'][:3])
        npz_write(path,arrays)
        meta.update(status='ok',generation_seconds=time.perf_counter()-start,
                    reference_refinement_seconds=time.perf_counter()-ref_start,
                    high_reference_relative_difference=sensitivity,
                    held_high_reference_relative_difference=held_sensitivity,
                    reference_sensitivity_below_one_percent=bool(max(sensitivity,held_sensitivity)<.01),
                    source_records=sources,refined_source_record=record96,
                    data_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    except Exception as exc:
        meta.update(status='failed',error_type=type(exc).__name__,error=str(exc))
        mf.write_json(meta_path,meta)
        raise
    mf.write_json(meta_path,meta)
    return arrays,meta

def make_mode(shape,pilot,sigma):
    start=time.perf_counter()
    predictions=[]
    for n in (16,32):
        model=ShapeModel(shape,n,(18.,))
        predictions.append(model.field_jac(np.array(pilot))[0][0])
        del model
        gc.collect()
    d,di=dw.mode_vectors(predictions[1]-predictions[0],sigma)
    return dict(delta=d.tolist(),delta_i=di.tolist(),seconds=time.perf_counter()-start,
                pilot=pilot,trace_inflation=float(d@d))

def fit(row,shape,arrays,startpoint,w,save):
    choice='low' if row['method']=='low' else 'low_high'
    start=time.perf_counter()
    model=None
    row['current_stage']='fit_setup'
    save()
    try:
        model=ShapeModel(shape,32,mf.CHOICES[choice])
        row['setup_seconds']=time.perf_counter()-start
        row['current_stage']='fit_solve'
        save()
        data=dict(sigma=float(arrays['sigma']),observations={choice:arrays[choice]},
                  references={choice:arrays['reference']})
        cache={}
        def ev(z):
            if cache.get('key')!=z.tobytes():
                r,j=mf.residual_jac(z,model,data,choice,row['use_reference'])
                if w is not None:
                    idx=dw.high_rows(144)
                    r[idx]=w@r[idx];j[idx]=w@j[idx]
                cache.update(key=z.tobytes(),r=r,j=j)
            return cache['r'],cache['j']
        solve=time.perf_counter()
        try:
            opt=least_squares(lambda z:ev(z)[0],np.array(startpoint),jac=lambda z:ev(z)[1],
                bounds=(mf.LO,mf.HI),x_scale=mf.SCALE,max_nfev=35,ftol=1e-9,xtol=1e-9,gtol=1e-7)
        finally:
            row['solve_seconds']=time.perf_counter()-solve
        row.update(estimated=opt.x.tolist(),status='converged' if opt.status>0 else 'iteration_limit',
                   optimizer_status=int(opt.status),nfev=int(opt.nfev),objective=float(opt.cost),
                   evaluation_status='pending',training_work=model.ledger())
        save()
    finally:
        row['fit_attempt_seconds']=time.perf_counter()-start
        if model is not None:row['training_work']=model.ledger()
        del model
        gc.collect()

def evaluate(row,shape,arrays,save):
    row['evaluation_status']='started';row['current_stage']='evaluation'
    save()
    start=time.perf_counter()
    try:
        z=np.array(row['estimated']);truth=arrays['true']
        model=ShapeModel(shape,32,mf.LOW,mf.receivers(17,1.6))
        pred,_=model.field_jac(z)
        structural=pred/mf.electronics(z,mf.LOW)[:,None,None,:]
        row.update(pose_error_m=float(np.linalg.norm(z[1:4]-truth[1:4])),
                   material_relative_error=float(abs(z[0]-truth[0])/truth[0]),
                   sensor_low_band=mf.metrics(pred,arrays['held_sensor_mean']),
                   structural_low_band=mf.metrics(structural,arrays['held_structural_mean']))
        del model
        model=ShapeModel(shape,32,(18.,))
        high,_=model.field_jac(z)
        row['high_raw_fit']=mf.metrics(high[0],arrays['low_high'][3])
        row['evaluation_status']='ok';row['current_stage']='complete'
    finally:
        row['evaluation_attempt_seconds']=time.perf_counter()-start
    save()

def run():
    hashes=provenance()
    manifest=OUT/'manifest.json'
    if manifest.exists():assert json.loads(manifest.read_text())==hashes
    mf.write_json(manifest,hashes)
    test()
    path=OUT/'fits.json'
    rows=json.loads(path.read_text()) if path.exists() else []
    assert all(r['provenance']==hashes for r in rows)
    for seed,shape,epsilon in CASES:
        arrays,meta=case_data(seed,shape,epsilon,hashes)
        for ref in (False,True):
            for method in ('low','warm_raw','isotropic','rank1'):
                prior=[r for r in rows if (r['seed'],r['use_reference'],r['method'])==(seed,ref,method)]
                assert len(prior)<=1
                previous=prior[0] if prior else None
                action=mf.checkpoint_action(previous)
                if action=='skip_completed':continue
                if action=='resume_evaluation':
                    row=previous
                    row['attempt_history']=mf.preserve_interruption(row,action)
                else:
                    row=dict(seed=seed,shape=shape,epsilon=epsilon,use_reference=ref,method=method,
                             status='started',provenance=hashes,data_sha256=meta['data_sha256'],
                             scope='factorial known-support one-material development stress test')
                    if previous:
                        row['attempt_history']=mf.preserve_interruption(previous,action)
                        rows[rows.index(previous)]=row
                    else:rows.append(row)
                def save():mf.write_json(path,rows)
                save()
                start=time.perf_counter()
                try:
                    if action!='resume_evaluation':
                        w=None;initial=mf.BASE.copy()
                        if method!='low':
                            row['current_stage']='mode_construction'
                            pilot=next(r for r in rows if r['seed']==seed and r['use_reference']==ref and r['method']=='low')
                            if 'estimated' not in pilot:raise RuntimeError('Pilot fit failed; no replacement start')
                            initial=pilot['estimated']
                            row['historical_pilot_fit_seconds']=pilot.get('fit_attempt_seconds')
                            mode_path=OUT/f'mode_{seed}_{int(ref)}.json'
                            if mode_path.exists():
                                mode=json.loads(mode_path.read_text())
                                assert mode['provenance']==hashes and mode['pilot']==initial
                            else:
                                mode=make_mode(shape,initial,float(arrays['sigma']))
                                mode['provenance']=hashes
                                mf.write_json(mode_path,mode)
                            row['mode_construction_seconds_shared']=mode['seconds']
                            w,_=dw.whitening(mode['delta'],mode['delta_i'],method)
                        row['start']=list(initial)
                        fit(row,shape,arrays,initial,w,save)
                    evaluate(row,shape,arrays,save)
                except Exception as exc:
                    row.update(failure_stage=row.get('current_stage'),error_type=type(exc).__name__,error=str(exc))
                    if 'estimated' in row:row['evaluation_status']='failed'
                    else:row.update(status='failed',evaluation_status='not_reached')
                row['attempt_seconds']=time.perf_counter()-start
                if row.get('attempt_history'):row['uncheckpointed_interrupted_cost']='unknown, not zero'
                save()
                print(json.dumps({k:row.get(k) for k in ('seed','method','use_reference','status','evaluation_status','failure_stage')}),flush=True)
                gc.collect()
    assert provenance()==hashes

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.run:
        with (OUT/'benchmark.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            run()
    else:test()
