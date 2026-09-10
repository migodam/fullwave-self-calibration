"""Explicitly post-hoc controls. Oracles are not deployable competing methods."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.optimize import least_squares
from multipole import Cluster
from run_v2 import fit_model,encode,STARTS,CENTERS,RADII,LOSS


def decode(v):return np.array(v['real'])+1j*np.array(v['imag'])


def oracle_fit(model,rx,y,sigma,gain,truth,kind):
    fits=[];before=time.perf_counter()
    pose_fixed=kind=='gain_geometry'
    def residual(p):
        t=np.r_[p,truth[2]] if pose_fixed else p
        f=model.forward(t,rx)
        g=gain if kind!='amplitude' else abs(gain)*np.exp(1j*np.angle(np.vdot(f,y)))
        r=(y-g*f)/sigma
        return np.sqrt(2)*np.r_[r.real,r.imag]
    for initial in STARTS:
        init=initial[:2] if pose_fixed else initial
        lo=[1.5,2] if pose_fixed else [1.5,2,-2]
        hi=[4,5] if pose_fixed else [4,5,2]
        fit=least_squares(residual,init,bounds=(lo,hi),max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8)
        t=np.r_[fit.x,truth[2]] if pose_fixed else fit.x
        fits.append({'theta':t.tolist(),'rss_real_whitened':float(fit.fun@fit.fun),
                     'nfev':int(fit.nfev),'optimizer_success':bool(fit.success)})
    best=min(fits,key=lambda v:v['rss_real_whitened']);theta=np.array(best['theta'])
    g=gain if kind!='amplitude' else abs(gain)*np.exp(1j*np.angle(np.vdot(model.forward(theta,rx),y)))
    return dict(theta=theta.tolist(),gain=encode(g),starts=fits,
                rss_real_whitened=best['rss_real_whitened'],seconds=time.perf_counter()-before)


def main():
    here=Path(__file__).parent;target=here/'results/attribution_controls.json'
    if target.exists():raise FileExistsError('Immutable output exists')
    source=here/'results/recovery_v2.json';data=json.loads(source.read_text())
    model=Cluster(CENTERS,RADII,18,3);truth_model=Cluster(CENTERS,RADII,18,7)
    rx=np.array(data['rx']);rows=[];begin=time.perf_counter()
    for row in data['rows']:
        truth=np.array(row['truth']);g=complex(decode(row['true_gain']));y=decode(row['y'])
        sigma=row['sigma_complex'];z=complex(decode(row['reference']));f_dda=decode(row['true_base_field'])
        noise=y-g*f_dda;model_data=g*truth_model.forward(truth,rx)+noise
        controls={
            'same_model_no_reference':fit_model(model,rx,model_data,sigma),
            'same_model_noisy_reference':fit_model(model,rx,model_data,sigma,z),
            'DDA_exact_amplitude':oracle_fit(model,rx,y,sigma,g,truth,'amplitude'),
            'DDA_exact_complex_gain':oracle_fit(model,rx,y,sigma,g,truth,'gain'),
            'DDA_exact_gain_geometry':oracle_fit(model,rx,y,sigma,g,truth,'gain_geometry'),
            'DDA_noiseless_exact_gain_geometry':oracle_fit(model,rx,g*f_dda,sigma,g,truth,'gain_geometry'),
            'DDA_wrong_reference':fit_model(model,rx,y,sigma,z*1.02*np.exp(.03j))}
        for name,v in controls.items():
            t=np.array(v['theta']);gg=complex(decode(v['gain']));er=np.abs(t-truth)
            v.update(material_errors=er[:2].tolist(),material_success=bool(np.all(er[:2]<=.1)),
                     geometry_error_mm=float(er[2]),gain_relative_error=float(abs(gg-g)/abs(g)),
                     scientific_acceptance='unresolved; oracle/post-hoc diagnostic, not certified recovery')
        rows.append({'scene':row['scene'],'controls':controls})
        print(row['scene'],{k:v['material_success'] for k,v in controls.items()},flush=True)
    summary={k:{'material_successes':sum(r['controls'][k]['material_success'] for r in rows),
                'median_material_errors':np.median([r['controls'][k]['material_errors'] for r in rows],axis=0).tolist(),
                'median_geometry_error_mm':float(np.median([r['controls'][k]['geometry_error_mm'] for r in rows])),
                'median_gain_relative_error':float(np.median([r['controls'][k]['gain_relative_error'] for r in rows]))}
             for k in rows[0]['controls']}
    obj={'scope':'Post-hoc mechanism controls, separate protocol; no population or method-superiority claim',
         'source_data_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'summary':summary,'rows':rows,'wall_seconds':time.perf_counter()-begin}
    target.write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
