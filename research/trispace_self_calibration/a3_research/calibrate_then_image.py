"""Frozen-geometry low-band refits; true geometry is a labelled oracle only."""
import argparse
import fcntl
import gc
import hashlib
import json
from pathlib import Path
import time
import shape_material_stress as ss
import numpy as np
from scipy.optimize import least_squares

mf=ss.mf
ROOT=ss.ROOT
SOURCE=ss.OUT
OUT=ROOT/'results/calibrate_then_image'
IDX=np.array([0]+list(range(4,13)))
SOURCES=('warm_raw','isotropic','rank1','oracle_true_geometry')

def provenance():
    p=ss.provenance()
    for path in [Path(__file__),ROOT/'CALIBRATE_THEN_IMAGE_PROTOCOL.md',SOURCE/'fits.json',
                 *[SOURCE/f'data_{seed}.npz' for seed,_,_ in ss.CASES]]:
        p[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    return p

def load():
    rows=json.loads((SOURCE/'fits.json').read_text())
    assert len(rows)==32 and all(r.get('evaluation_status')=='ok' for r in rows)
    manifest=json.loads((SOURCE/'manifest.json').read_text())
    for path,digest in manifest.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
    return {(r['seed'],r['use_reference'],r['method']):r for r in rows}

def ev_factory(base,shape,n,arrays,use_ref):
    model=ss.ShapeModel(shape,n,mf.LOW)
    data=dict(sigma=float(arrays['sigma']),observations={'low':arrays['low']},references={'low':arrays['reference']})
    cache={}
    def ev(v):
        if cache.get('key')!=v.tobytes():
            z=base.copy();z[IDX]=v
            assert np.array_equal(z[1:4],base[1:4])
            r,j=mf.residual_jac(z,model,data,'low',use_ref)
            cache.update(key=v.tobytes(),r=r,j=j[:,IDX])
        return cache['r'],cache['j']
    return ev,model

def test():
    hashes=provenance();load()
    with np.load(SOURCE/'data_8403.npz') as data:
        arrays={k:data[k].copy() for k in data.files}
    base=np.array([2.3,.03,-.02,.04,.03]+[.01]*4+[.02]*4)
    records=[]
    for ref in (False,True):
        ev,model=ev_factory(base,'box',8,arrays,ref)
        r,j=ev(base[IDX]);errors=[]
        for c in range(10):
            step=np.eye(10)[c]*1e-5
            fd=(ev(base[IDX]+step)[0]-ev(base[IDX]-step)[0])/2e-5
            errors.append(float(np.linalg.norm(fd-j[:,c])/np.linalg.norm(fd)))
        assert max(errors)<1e-7
        records.append(dict(use_reference=ref,column_errors=errors))
    assert provenance()==hashes
    mf.write_json(OUT/'checks.json',dict(passed=True,records=records,provenance=hashes,
                                       fixed_geometry_exact=True))
    print(json.dumps(dict(checks_passed=True,max_error=max(max(r['column_errors']) for r in records))),flush=True)

def fit(row,shape,arrays,save):
    row['current_stage']='fit_setup';save()
    start=time.perf_counter();model=None
    try:
        base=np.array(row['start'])
        ev,model=ev_factory(base,shape,32,arrays,row['use_reference'])
        row['setup_seconds']=time.perf_counter()-start
        row['current_stage']='fit_solve';save()
        solve=time.perf_counter()
        try:
            opt=least_squares(lambda v:ev(v)[0],base[IDX],jac=lambda v:ev(v)[1],
                bounds=(mf.LO[IDX],mf.HI[IDX]),x_scale=mf.SCALE[IDX],max_nfev=35,
                ftol=1e-9,xtol=1e-9,gtol=1e-7)
        finally:row['solve_seconds']=time.perf_counter()-solve
        estimated=base.copy();estimated[IDX]=opt.x
        assert np.array_equal(estimated[1:4],base[1:4])
        row.update(estimated=estimated.tolist(),optimizer_status=int(opt.status),nfev=int(opt.nfev),
                   objective=float(opt.cost),status='converged' if opt.status>0 else 'iteration_limit',
                   training_work=model.ledger(),evaluation_status='pending')
        save()
    finally:
        row['fit_attempt_seconds']=time.perf_counter()-start
        if model is not None:row['training_work']=model.ledger()
        del model;gc.collect()

def run():
    hashes=provenance();prior=load()
    path=OUT/'fits.json';manifest=OUT/'manifest.json'
    if manifest.exists():assert json.loads(manifest.read_text())==hashes
    mf.write_json(manifest,hashes);test()
    rows=json.loads(path.read_text()) if path.exists() else []
    assert all(r['provenance']==hashes for r in rows)
    for seed,shape,epsilon in ss.CASES:
        with np.load(SOURCE/f'data_{seed}.npz') as archive:
            arrays={k:archive[k].copy() for k in archive.files}
        for ref in (False,True):
            pilot=prior[seed,ref,'low']
            for source in SOURCES:
                previous=[r for r in rows if (r['seed'],r['use_reference'],r['geometry_source'])==(seed,ref,source)]
                assert len(previous)<=1
                previous=previous[0] if previous else None
                action=mf.checkpoint_action(previous)
                if action=='skip_completed':continue
                if action=='resume_evaluation':
                    row=previous;row['attempt_history']=mf.preserve_interruption(row,action)
                else:
                    base=np.array(pilot['estimated'])
                    oracle=source=='oracle_true_geometry'
                    calibration=None if oracle else prior[seed,ref,source]
                    base[1:4]=arrays['true'][1:4] if oracle else np.array(calibration['estimated'])[1:4]
                    row=dict(seed=seed,shape=shape,epsilon=epsilon,use_reference=ref,geometry_source=source,
                        oracle=oracle,start=base.tolist(),status='started',provenance=hashes,
                        historical_pilot_fit_seconds=pilot.get('fit_attempt_seconds'),
                        historical_calibration_fit_seconds=None if oracle else calibration.get('fit_attempt_seconds'),
                        historical_mode_seconds=None if oracle else calibration.get('mode_construction_seconds_shared'),
                        scope='Fixed estimated geometry, low-band nuisance/material refit; oracle separately labelled')
                    if previous:
                        row['attempt_history']=mf.preserve_interruption(previous,action)
                        rows[rows.index(previous)]=row
                    else:rows.append(row)
                def save():mf.write_json(path,rows)
                save();start=time.perf_counter()
                try:
                    if action!='resume_evaluation':fit(row,shape,arrays,save)
                    ss.evaluate(row,shape,arrays,save)
                    assert np.array_equal(np.array(row['estimated'])[1:4],np.array(row['start'])[1:4])
                    row['fixed_geometry_exact']=True
                except Exception as exc:
                    row.update(failure_stage=row.get('current_stage'),error_type=type(exc).__name__,error=str(exc))
                    if 'estimated' in row:row['evaluation_status']='failed'
                    else:row.update(status='failed',evaluation_status='not_reached')
                row['attempt_seconds']=time.perf_counter()-start
                if row.get('attempt_history'):row['uncheckpointed_interrupted_cost']='unknown, not zero'
                save()
                print(json.dumps({k:row.get(k) for k in ('seed','use_reference','geometry_source','status','evaluation_status')}),flush=True)
                gc.collect()
    assert provenance()==hashes

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.run:
        with (SOURCE/'benchmark.lock').open('r+b') as source_lock,(OUT/'benchmark.lock').open('a') as lock:
            fcntl.flock(source_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            run()
    else:test()
