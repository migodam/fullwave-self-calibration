"""Run fixed forty-scene confirmatory benchmark without inspecting interim outcomes."""
import gc
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import timed_rom as timed

OUT=Path(__file__).resolve().parent
sm=timed.seedmod


def case(seed):
    if seed not in range(9201,9241):raise ValueError('unregistered final seed')
    streams=np.random.SeedSequence(seed).spawn(3)
    mat,noise,pose=[np.random.default_rng(s) for s in streams]
    alpha=np.full(9,.08);n=int(mat.integers(2,4))
    alpha[mat.choice(9,n,replace=False)]+=mat.uniform(.35,.90,n)
    model=sm.model_for(16,9)
    background=model.forward(np.full(9,.08),np.zeros(3),list(sm.ctl.FREQ_IDS),jacobian=False)['total']
    sigma=float(np.sqrt(np.mean(abs(background)**2)/1000))
    truth=model.forward(alpha,np.zeros(3),jacobian=False)
    y=truth['total']+sigma/np.sqrt(2)*(noise.normal(size=truth['total'].shape)+1j*noise.normal(size=truth['total'].shape))
    angle=pose.uniform(0,2*np.pi);radius=pose.uniform(.04,.08)
    theta=pose.uniform(.01,.05)*(1 if pose.random()<.5 else -1)
    return dict(alpha_true=alpha,x_true=np.zeros(3),y=y,sigma=sigma,
        x0=np.array([radius*np.cos(angle),radius*np.sin(angle),theta]))


def main():
    dest=OUT/'results/rom_final.json'
    digest=hashlib.sha256((OUT/'ROM_FINAL_PROTOCOL.md').read_bytes()).hexdigest()
    source_paths=[OUT/'rom_final.py',OUT/'timed_rom.py',OUT/'rom_seed_comparison.py',
                  sm.REUSE/'controller.py',sm.REUSE/'cost.py']
    source_paths+=list((sm.REUSE.parent/'a3_rom').glob('*.py'))
    source_paths+=list((sm.REUSE.parent/'a2_physics').glob('*.py'))
    hashes={str(p.relative_to(OUT.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    code_digest=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    rows=json.loads(dest.read_text()) if dest.exists() else []
    if rows and any(r['protocol_sha256']!=digest or r['source_digest']!=code_digest for r in rows):
        raise RuntimeError('protocol or source changed after results: do not resume')
    sm.ctl.write_json(OUT/'results/rom_final_source_manifest.json',dict(protocol_sha256=digest,
        source_digest=code_digest,source_hashes=hashes))
    for seed in range(9201,9241):
        data=case(seed);offset=(seed-9201)%len(sm.METHODS)
        order=sm.METHODS[offset:]+sm.METHODS[:offset]
        for method in order:
            if any(r['seed']==seed and r['method']==method for r in rows):continue
            started=time.perf_counter()
            try:
                model=sm.model_for(16,9);setup=time.perf_counter()-started
                row=timed.solve(model,data,method,8.,setup,maxiter=200)
                row['joint_success']=row['pose_error_m']<=.05 and row['material_rmse']<=.05
            except Exception as exc:
                row=dict(method=method,status='exception',joint_success=False,
                    error_type=type(exc).__name__,error=str(exc),
                    actual_consumed_seconds=time.perf_counter()-started)
            row.update(seed=seed,protocol_sha256=digest,source_digest=code_digest)
            rows.append(row);sm.ctl.write_json(dest,rows)
            # Progress only: do not inspect outcomes before all scenes finish.
            print(json.dumps(dict(completed=len(rows),planned=240,seed=seed,method=method)),flush=True)
            if 'model' in locals():del model
            gc.collect()


if __name__=='__main__':main()
