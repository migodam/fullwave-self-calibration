"""Full-Maxwell mechanism tests, plus exact modal/linear negative controls."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
import maxwell3d as mx
sys.path.insert(0,str(ROOT/'public_release/research/a4_reliability_v1/src'))
import modal


def profiled(y, a, gain):
    y=np.asarray(y); a=np.asarray(a)
    if gain=='entry':
        return np.where(np.abs(y)>1e-290,0.,a)
    out=a.copy()
    groups=[(...,)] if gain=='shared' else [(slice(None),slice(None),p) for p in range(y.shape[-1])]
    for idx in groups:
        v=y[idx].ravel(); w=a[idx].ravel()
        out[idx]=(w-v*np.vdot(v,w)/np.vdot(v,v)).reshape(a[idx].shape)
    return out


def run():
    out=HERE/'results';out.mkdir(exist_ok=True)
    paths=[Path(__file__),HERE/'PLAN_AND_PROTOCOL.md',Path(mx.__file__),Path(modal.__file__)]
    manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    target=out/'scale_gauge.json'
    if target.exists():raise RuntimeError('Preserve existing run; use a separately registered version')
    start=time.perf_counter()
    checks=[]
    # A4's exact full-wave radial response, not a Born approximation.
    for k in (4.,12.,24.):
        for radius in (.12,.4,1.):
            r=radius*np.array([.3,.4,np.sqrt(.75)])
            y=modal.sphere_field(r,k,2+.03j,.045).ravel()
            z=modal.sphere_field(r,k,4+.03j,.045).ravel()
            g=np.vdot(z,y)/np.vdot(z,z)
            mismatch=float(np.linalg.norm(y-g*z)/np.linalg.norm(y))
            assert mismatch<1e-12
            checks.append(dict(k=k,radius=radius,gain_compensated_relative_error=mismatch,
                               compensating_gain=[float(g.real),float(g.imag)]))
    rows=[]
    def save(complete=False):
        record=dict(complete=complete,source_hashes=manifest,modal_counterexamples=checks,
                    rows=rows,wall_seconds=time.perf_counter()-start,
                    scope='Deterministic local/finite-pair diagnostics, not recovery, global stability or novelty proof.')
        temporary=target.with_suffix('.tmp')
        temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(target)
    save()
    centers=[[-.06,0.,0.],[.055,.02,0.]];radii=[.035,.025]
    template=np.array([1+.03j,2+.05j]);k=18.
    for scale in (.01,.03,.1,.3,1.,2.):
        for radius in (.2,.6,2.):
            fields={}
            for order in (3,4):
                fs=[]
                for multiplier in (1.,np.exp(1e-4),np.exp(-1e-4),1.25):
                    y,_=mx.treams_field(centers,radii,1+scale*multiplier*template,k,
                                        mx.receivers(12,radius),lmax=order)
                    fs.append(y)
                fields[order]=fs
            y,plus,minus,other=fields[4]; derivative=(plus-minus)/2e-4
            result=dict(scale=scale,radius=radius,field_norm=float(np.linalg.norm(y)),
                        lmax_relative_change=float(np.linalg.norm(fields[3][0]-y)/np.linalg.norm(y)),gains={})
            for gain in ('shared','per_illumination','entry'):
                visible=profiled(y,derivative,gain)
                # min_g ||y-g*other||, not division by noisy data.
                finite=profiled(other,y,gain)
                low_derivative=(fields[3][1]-fields[3][2])/2e-4
                low_visible=profiled(fields[3][0],low_derivative,gain)
                result['gains'][gain]=dict(
                    absolute_log_scale_sensitivity=float(np.linalg.norm(visible)),
                    relative_log_scale_sensitivity=float(np.linalg.norm(visible)/np.linalg.norm(y)),
                    finite_pair_relative_residual=float(np.linalg.norm(finite)/np.linalg.norm(y)),
                    lmax3_relative_sensitivity=float(np.linalg.norm(low_visible)/np.linalg.norm(fields[3][0])))
            linear_control=float(np.linalg.norm(profiled(1.25*y,y,'shared'))/np.linalg.norm(y))
            assert linear_control<1e-12
            result['linear_scaling_control']=linear_control
            assert result['gains']['entry']['relative_log_scale_sensitivity']==0.
            rows.append(result);save()
            print(json.dumps(result),flush=True)
    assert manifest=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    save(True)


if __name__=='__main__':run()
