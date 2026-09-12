"""Resume interrupted V2 only from immutable saved raw observations.

No thresholds, starts, model settings, or data are changed. The hosted process
was killed by its 200 s execution timeout after completing scenes 0 through 10.
Existing scene files are loaded and never overwritten.
"""
import argparse, hashlib, json, platform, resource, time
from pathlib import Path
import numpy as np
import recovery_v2 as r

def resume(output):
    output=Path(output)
    if (output/'summary.json').exists():raise FileExistsError(output/'summary.json')
    clock=time.perf_counter();rawfile=output/'raw_observations.json'
    raw_sha=hashlib.sha256(rawfile.read_bytes()).hexdigest();pack=json.loads(rawfile.read_text())
    raw=pack['scenes'];rx=np.array(pack['rx']);candidate_rx=np.array(pack['candidate_rx'])
    model=r.SphereCluster(order=3,ntheta=12);results=[];resumed=[]
    for row in raw:
        i=row['scene'];target=output/f'scene_{i:02d}.json'
        if target.exists():results.append(json.loads(target.read_text()));continue
        theta=np.array(row['truth']);g=r.complex_read(row['gain']);y=r.complex_read(row['y'])
        field=r.complex_read(row['noiseless_field']);extra=r.complex_read(row['candidate_y']);z=r.complex_read(row['reference'])
        sigma=row['sigma'];res={'scene':i,'truth':theta.tolist(),'methods':{},'oracle_diagnostics':{}}
        for name,kwargs in [('no_reference',{}),('noisy_reference_GLS',{'z':z}),
            ('fixed_EM',{'candidate_points':candidate_rx,'extra_index':0,'extra_y':extra[0]}),
            ('random_EM',{'candidate_points':candidate_rx,'extra_index':row['random_index'],'extra_y':extra[row['random_index']]})]:
            res['methods'][name]=r.annotate(r.fit(model,rx,y,sigma,**kwargs),theta,g)
        for name,kwargs in [('noiseless_unknown_gain',{}),('noiseless_known_gain',{'gain_known':g}),
                            ('noiseless_known_gain_geometry',{'gain_known':g,'geometry_known':theta[2]})]:
            res['oracle_diagnostics'][name]=r.annotate(r.fit(model,rx,g*field,sigma,**kwargs),theta,g)
        res['structural_error']=r.orthogonal_diagnostic(model,theta,rx,field,g)
        target.write_text(json.dumps(res,indent=2)+'\n');results.append(res);resumed.append(i)
        print('resumed',i,flush=True)
    convergence=[]
    for spacing in [.009,.0075]:
        fine=r.DipoleVIE(r.CENTERS,r.RADII,spacing,r.K,4)
        for row in raw[:4]:
            theta=np.array(row['truth']);shift=np.array([theta[2]*.001,0,0])
            f=fine.field(theta[:2]+1j*r.LOSS,rx+shift).ravel();fm=model.forward(theta,rx)
            coarse=r.complex_read(row['noiseless_field'])
            convergence.append({'scene':row['scene'],'spacing':spacing,'cells':len(fine.points),
                'relative_difference_to_coarse_DDA':float(np.linalg.norm(f-coarse)/np.linalg.norm(f)),
                'relative_difference_to_lmax3':float(np.linalg.norm(f-fm)/np.linalg.norm(f))})
        del fine
    order_checks=[];m4=r.SphereCluster(order=4,ntheta=14);m5=r.SphereCluster(order=5,ntheta=16)
    for row in raw[:4]:
        theta=np.array(row['truth']);f3=model.forward(theta,rx);f4=m4.forward(theta,rx);f5=m5.forward(theta,rx)
        order_checks.append({'scene':row['scene'],'l3_l4':float(np.linalg.norm(f3-f4)/np.linalg.norm(f4)),
            'l4_l5':float(np.linalg.norm(f4-f5)/np.linalg.norm(f5))})
    summary={};oracle={}
    for group,table in [('methods',summary),('oracle_diagnostics',oracle)]:
        for name in results[0][group]:
            valid=[s[group][name] for s in results if 'failure' not in s[group][name]]
            table[name]={'successes':sum(s['material_success'] for s in valid),'attempts':12,
                'median_material_errors':np.median([s['material_errors'] for s in valid],axis=0).tolist(),
                'median_geometry_error_mm':float(np.median([s['geometry_error_mm'] for s in valid])),
                'total_forward_evaluations':sum(s['forward_evaluations'] for s in valid),
                'total_wall_seconds':sum(s['wall_seconds'] for s in valid)}
    sources=[r.HERE/'recovery_v2.py',r.HERE/'sphere_cluster.py',Path(__file__),r.ROOT/'research/trispace_self_calibration/a3_research/maxwell3d.py']
    final={'kind':'registered_V2_independent_diagnostic_not_final_test','status':'executed_with_logged_timeout_resume',
        'summary':summary,'oracle_diagnostics':oracle,'convergence_diagnostics_not_bounds':convergence,
        'order_checks_not_bounds':order_checks,'original_process_timeout_seconds':200,
        'resumed_scenes':resumed,'resume_wall_seconds':time.perf_counter()-clock,
        'end_to_end_wall_seconds':None,'end_to_end_wall_note':'original process interrupted; no aggregate runtime superiority claim',
        'maxrss_kib_linux_resume_process':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'python':platform.python_version(),'numpy':np.__version__,
        'source_sha256':{str(p.relative_to(r.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'raw_sha256':raw_sha,'limitations':['No certified DDA continuum error','No class C covering finite bound',
          'No hardware experiment','Equal extra scalar count is not equal electronics/EM hardware cost']}
    assert hashlib.sha256(rawfile.read_bytes()).hexdigest()==raw_sha
    (output/'summary.json').write_text(json.dumps(final,indent=2)+'\n');print(json.dumps(final,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=r.HERE/'results/recovery_v2')
    resume(p.parse_args().out)
