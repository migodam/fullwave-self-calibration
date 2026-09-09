"""Audited tuning/freeze/final runner; separate immutable records per budget."""
import argparse,hashlib,json,time,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE),str(ROOT/'research/delegated/a2_physics'),str(ROOT/'research/trispace_self_calibration/a2_research')]
from common import jsonable,Calibration,default_settings
from run_e4 import make_models,run_one
from generate import nominal_sigma,generate_seed_scene,stratum_for_seed,init_offsets
from reference import reference_metrics

METHODS=['direct','prasc','fixed_rank','phaseless','direct_control']
def dump(p,v):
    p.write_text(json.dumps(jsonable(v),indent=2)+'\n')
def code_hash():
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))}
def candidates(method):
    opts=[{'maxls':m,'ftol':f,'gtol':1e-6,'init_max_units':0,'max_reduced_probes':1} for m in [10,20] for f in [1e-8,1e-10]]
    if method=='fixed_rank':
        for o,r in zip(opts,[16,36,64,36]):o['fixed_rank']=r
    return opts
def tune():
    models=make_models()
    cal=json.loads((HERE.parent/'a2_solver/calibration_tune.json').read_text())
    # Keep previously measured tuning-only ratios, not recalibrated after tests.
    dump(HERE/'calibration_tune.json',cal)
    sig={a:nominal_sigma(models[(a,32)]) for a in ['full','limited']}
    dump(HERE/'noise_sigma_tune.json',{'sigma_by_aperture':sig})
    ref={a:reference_metrics(models[(a,16)],models[(a,32)],sig[a],a) for a in sig}
    assert all(r['full_rank'] for r in ref.values())
    dump(HERE/'reference_design.json',ref)
    records=[];target=HERE/'tuning_corrected.jsonl'
    if target.exists():raise RuntimeError('Do not overwrite tuning records')
    with target.open('w') as f:
        for seed in range(1,11):
            ap=stratum_for_seed(seed)[1];scene=generate_seed_scene(seed,models[(ap,32)],sig[ap])
            ii=seed%12;info={'index':ii,**init_offsets(seed)[ii]}
            for m in METHODS:
                for ci,cand in enumerate(candidates(m)):
                    rec=run_one(SimpleNamespace(budget=800),models,seed,scene,info,m,{**default_settings(),**cand},ref)
                    rec['candidate']=ci;records.append(rec);f.write(json.dumps(jsonable(rec))+'\n');f.flush()
            print(f'tuning seed {seed}/10 done',flush=True)
    settings={}
    for m in METHODS:
        rr=[r for r in records if r['method']==m]
        def score(i):
            z=[r for r in rr if r['candidate']==i]
            return (-np.mean([r['success'] for r in z]),np.mean([r['metrics']['loss_exact_final'] for r in z]),i)
        i=min(range(4),key=score);settings[m]={**default_settings(),**candidates(m)[i]}
    frozen={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'code_hash':code_hash(),'methods':METHODS,'budgets':[200,800],'reference_design':ref,'sigma_by_aperture':sig,'calibration':cal,'settings':settings,'seeds_tune':list(range(1,11)),'seeds_final':list(range(1001,1021)),'initialization':'constant_nonoracle_alpha_0.5','primary':'200 only; secondary800 cannot rescue','bootstrap_seed':20260906,'bootstrap_draws':10000,'alpha_one_sided':.05/3,'noninferiority':'per-aperture difference divided by frozen margin, upper bound <1'}
    dump(HERE/'frozen_parent.json',frozen)
    print('FROZEN',hashlib.sha256((HERE/'frozen_parent.json').read_bytes()).hexdigest(),flush=True)
def final(budget):
    frozen=json.loads((HERE/'frozen_parent.json').read_text())
    assert frozen['code_hash']==code_hash(),'Code changed after freeze'
    assert budget in frozen['budgets']
    models=make_models();target=HERE/f'final_{budget}.jsonl'
    if target.exists():raise RuntimeError('Final output exists; never append duplicate runs')
    with target.open('w') as f:
        for seed in frozen['seeds_final']:
            ap=stratum_for_seed(seed)[1];scene=generate_seed_scene(seed,models[(ap,32)],frozen['sigma_by_aperture'][ap],allow_final=True)
            for ii,off in enumerate(init_offsets(seed)):
                for m in METHODS:
                    rec=run_one(SimpleNamespace(budget=budget),models,seed,scene,{'index':ii,**off},m,frozen['settings'][m],frozen['reference_design'],frozen)
                    # Only numerical sensing rank was computed, never L_det.
                    rec['sensing_numerical_rank']=rec.pop('L_det',None);rec['L_det']=None
                    f.write(json.dumps(jsonable(rec))+'\n');f.flush()
            print(f'final budget {budget}, seed {seed}, 60 paired records',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['tune','final']);p.add_argument('--budget',type=int,default=200);args=p.parse_args()
    tune() if args.mode=='tune' else final(args.budget)
