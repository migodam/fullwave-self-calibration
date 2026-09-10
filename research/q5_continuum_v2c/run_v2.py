"""Registered V2-C experiment runner; no algorithm changes in delivery split."""
from recovery_common import *


def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'results/recovery_v2.json')
    args=p.parse_args();target=args.output
    if target.exists():raise FileExistsError('Existing evidence is immutable: '+str(target))
    rng=np.random.default_rng(2026091107)
    action_rng=np.random.default_rng(2026091108)
    # Pin all truths before drawing observations: reproducible, independent streams.
    truths=[(np.r_[rng.uniform([1.5,2],[4,5]),rng.uniform(-2,2)],
             rng.uniform(.75,1.25)*np.exp(1j*rng.uniform(-np.pi,np.pi))) for _ in range(12)]
    rx=directions(12,.6);val=directions(16,.6);struct=directions(36,.15)
    inv=Cluster(CENTERS,RADII,18,3)
    check5=Cluster(CENTERS,RADII,18,5)
    check7=Cluster(CENTERS,RADII,18,7)
    checkq=Cluster(CENTERS,RADII,18,3,22)
    dda=mx.DipoleVIE(CENTERS,RADII,.011,18,fill_quadrature=4)
    rows=[];start=time.perf_counter()
    source_paths=[Path(__file__),HERE/'recovery_common.py',HERE/'multipole.py',Path(mx.__file__)]
    source_hashes={str(s.relative_to(ROOT)):hashlib.sha256(s.read_bytes()).hexdigest() for s in source_paths}
    def save(complete=False):
        summary={method:{'material_successes':sum(r['methods'][method]['material_success'] for r in rows),
                         'median_material_errors':np.median([r['methods'][method]['material_errors'] for r in rows],axis=0).tolist() if rows else None}
                 for method in ('no_reference','noisy_reference_gls','fixed_em','random_em')}
        obj={'protocol':'Q5 V2 2026-09-11; material-only threshold 0.1',
             'complete':complete,'seed':2026091107,'extra_action_seed':2026091108,
             'fixed_extra_index':0,'truths_drawn_before_noise':True,
             'scope':'Registered development diagnostics; no continuum/model-error or population certificate.',
             'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                            'blas_threads':1,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
             'source_hashes':source_hashes,'dda_cells':len(dda.points),
             'rx':rx.tolist(),'extra_rx':val.tolist(),'structural_points':struct.tolist(),
             'summary':summary,'rows':rows,'wall_seconds':time.perf_counter()-start}
        target.parent.mkdir(parents=True,exist_ok=True)
        tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(target)
    for scene,(truth,gain) in enumerate(truths):
        gen=time.perf_counter()
        currents=dda.currents(truth[:2]+1j*LOSS)
        def field(points,shift=True):
            p=points+np.array([truth[2]*.001,0,0]) if shift else points
            return (mx.dipole_kernel(p,dda.points,18)@currents).ravel()
        yy=field(rx);vv=field(val);ss=field(struct,False)
        sigma=.01*np.linalg.norm(gain*yy)/np.sqrt(yy.size)
        noise=lambda n:sigma/np.sqrt(2)*(rng.normal(size=n)+1j*rng.normal(size=n))
        y=gain*yy+noise(len(yy));extras=gain*vv+noise(len(vv))
        reference=gain+.01/np.sqrt(2)*(rng.normal()+1j*rng.normal())
        idx=int(action_rng.integers(len(vv)))
        row={'scene':scene,'truth':truth.tolist(),'true_gain':encode(gain),'sigma_complex':sigma,
             'reference_sigma_complex':.01,'y':encode(y),'reference':encode(reference),
             'extra_observations':encode(extras),'true_base_field':encode(yy),
             'true_extra_field':encode(vv),'true_structural_field_fixed_world':encode(ss),
             'random_extra_index':idx,'generation_seconds':time.perf_counter()-gen,'methods':{}}
        for name,ref,extra in [('no_reference',None,None),('noisy_reference_gls',reference,None),
                              ('fixed_em',None,(val,0,extras[0])),('random_em',None,(val,idx,extras[idx]))]:
            result=fit_model(inv,rx,y,sigma,ref,extra)
            theta=np.array(result['theta']);g=complex(result['gain']['real'],result['gain']['imag'])
            err=np.abs(theta-truth)
            result.update(material_errors=err[:2].tolist(),material_success=bool(np.all(err[:2]<=.1)),
                          geometry_error_mm=float(err[2]),gain_relative_error=float(abs(g-gain)/abs(gain)),
                          heldout_noiseless_sensor_error=float(np.linalg.norm(g*inv.forward(theta,val)-gain*vv)/np.linalg.norm(gain*vv)),
                          structural_error_fixed_world=float(np.linalg.norm(inv.field(theta[:2]+1j*LOSS,struct).ravel()-ss)/np.linalg.norm(ss)),
                          scientific_acceptance='unresolved')
            row['methods'][name]=result
        if scene<4:
            conv={'scope':'solver differences, NOT upper error bounds','orders':{},'dda_grids':{}}
            baseline=inv.forward(truth,rx)
            for name,model in [('L5',check5),('L7',check7),('L3_q22',checkq)]:
                v=model.forward(truth,rx)
                conv['orders'][name]={'relative_to_L3':float(np.linalg.norm(v-baseline)/np.linalg.norm(v)),
                                      'relative_to_DDA011':float(np.linalg.norm(v-yy)/np.linalg.norm(v))}
            for spacing in [.009,.008]:
                begin=time.perf_counter()
                refined=mx.DipoleVIE(CENTERS,RADII,spacing,18,fill_quadrature=4)
                cc=refined.currents(truth[:2]+1j*LOSS)
                ff=(mx.dipole_kernel(rx+[truth[2]*.001,0,0],refined.points,18)@cc).ravel()
                conv['dda_grids'][str(spacing)]={'cells':len(refined.points),
                    'relative_to_DDA011':float(np.linalg.norm(ff-yy)/np.linalg.norm(ff)),
                    'relative_to_L7':float(np.linalg.norm(ff-check7.forward(truth,rx))/np.linalg.norm(check7.forward(truth,rx))),
                    'seconds':time.perf_counter()-begin}
                del refined
            row['convergence_diagnostics']=conv
        rows.append(row);save()
        print(json.dumps({'scene':scene,'material_success':{k:v['material_success'] for k,v in row['methods'].items()}}),flush=True)
    save(True)

if __name__=='__main__':main()
