"""Registered post-hoc diagnostic. It deliberately reuses development data."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
import json,time,hashlib
from scipy.optimize import least_squares
from sphere_solver import SphereSolver,receivers
from recovery import profile
HERE=Path(__file__).resolve().parents[1]

def main():
    outpath=HERE/'results/attribution_diagnostic.json'
    if outpath.exists():raise FileExistsError(outpath)
    rawpath=HERE/'results/observations_v2.npz';raw=np.load(rawpath)
    base=json.loads((HERE/'results/recovery_v2.json').read_text())
    solver=SphereSolver();rx=receivers();start=time.perf_counter();rows=[]
    for scene in range(12):
        truth=raw[f'{scene}_truth'];g=complex(raw[f'{scene}_gain']);y=raw[f'{scene}_y']
        sig=base['rows'][scene]['sigma'];row={'scene':scene,'oracles':{}}
        for name in ['true_gain','true_shift','true_gain_and_shift']:
            gain_fixed=name in ('true_gain','true_gain_and_shift')
            shift_fixed=name in ('true_shift','true_gain_and_shift')
            def theta(x):return np.r_[x,truth[2]] if shift_fixed else x
            def fun(x):
                f=solver.forward(theta(x),rx)
                if not gain_fixed:return profile(y,f,sig)[0]
                r=(y-g*f)/sig;return np.sqrt(2)*np.r_[r.real,r.imag]
            fits=[]
            for initial in ([2,3,0],[1.6,2.2,-1],[3.8,4.8,1]):
                init=initial[:2] if shift_fixed else initial
                lo=[1.5,2] if shift_fixed else [1.5,2,-2]
                hi=[4,5] if shift_fixed else [4,5,2]
                res=least_squares(fun,init,bounds=(lo,hi),max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8)
                t=theta(res.x)
                fits.append({'theta':t.tolist(),'rss':float(res.fun@res.fun),'nfev':res.nfev,
                    'optimizer_success':bool(res.success),'initial':init})
            best=min(fits,key=lambda x:x['rss']);err=abs(np.array(best['theta'])-truth)
            row['oracles'][name]={'all_starts':fits,'theta':best['theta'],'material_errors':err[:2].tolist(),
               'task_success':bool(np.all(err[:2]<=.1)),'geometry_error_mm':float(err[2])}
        rows.append(row)
    out={'scope':'POST-HOC DEVELOPMENT DIAGNOSTIC, not final or deployable performance',
         'source_data_sha256':hashlib.sha256(rawpath.read_bytes()).hexdigest(),
         'registration_sha256':hashlib.sha256((HERE/'ATTRIBUTION_DIAGNOSTIC_PROTOCOL.md').read_bytes()).hexdigest(),
         'rows':rows,'seconds':time.perf_counter()-start,
         'counts':{key:sum(row['oracles'][key]['task_success'] for row in rows) for key in rows[0]['oracles']}}
    outpath.write_text(json.dumps(out,indent=2)+'\n');print(out['counts'])
if __name__=='__main__':main()
