"""Frozen V2 independent-model diagnostic stream, not a final evaluation.

Outputs preserve raw observations. DDA and sphere-cluster differences are
numerical diagnostics, not rigorous continuum-model error bounds.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys, json, time, hashlib, platform, resource
from functools import lru_cache
import numpy as np
import scipy
from scipy.optimize import least_squares
from multipole import Cluster, boundary_residual
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
import maxwell3d as mx
CENTERS=np.array([[-.06,0,0],[.055,.02,0]])
RADII=np.array([.035,.025]); LOSS=np.array([.03,.05]); K=18.
STARTS=([2,3,0],[1.6,2.2,-1],[3.8,4.8,1])
BOUNDS=([1.5,2,-2],[4,5,2])


def encode(z):
    z=np.asarray(z)
    return np.stack([z.real,z.imag],axis=-1).tolist()


def profile(y,f,sigma,reference=None):
    """Annulus-constrained deterministic profiling, without a logdet term."""
    v=np.asarray(f,complex)/sigma; b=np.asarray(y,complex)/sigma
    if reference is not None:
        z,sr=reference
        v=np.r_[v,1/sr]; b=np.r_[b,z/sr]
    energy=np.vdot(v,v).real
    if energy<=0 or not np.isfinite(energy):
        raise ValueError('zero or invalid model energy')
    g0=np.vdot(v,b)/energy
    gain=float(np.clip(abs(g0),.75,1.25))*np.exp(1j*np.angle(g0))
    residual=b-v*gain
    return np.sqrt(2)*np.r_[residual.real,residual.imag],gain


def main(output=None):
    target=Path(output) if output else HERE/'results/independent_v2.json'
    if target.exists(): raise FileExistsError('Immutable output exists: '+str(target))
    target.parent.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(2026091107)
    action_rng=np.random.default_rng(2026091108)
    rx=mx.receivers(12,.6); pool=mx.receivers(16,.6)
    points=np.concatenate([rx,pool]); base_count=144
    model=Cluster(CENTERS,RADII,K,mx.illuminations(),order=3)
    dda=mx.DipoleVIE(CENTERS,RADII,.011,K,fill_quadrature=4)
    fine=mx.DipoleVIE(CENTERS,RADII,.009,K,fill_quadrature=4)
    high=[Cluster(CENTERS,RADII,K,mx.illuminations(),order=o) for o in (4,5)]
    quad=Cluster(CENTERS,RADII,K,mx.illuminations(),order=3,nt=24,nphi=48)
    @lru_cache(maxsize=128)
    def obs(shift):
        return model.observation_matrix(points+np.array([shift*.001,0,0]))
    def forward(theta):
        eps=tuple(np.asarray(theta[:2])+1j*LOSS)
        return (obs(float(theta[2]))@model.coefficients(eps)).ravel()
    sources=[Path(__file__),HERE/'multipole.py',Path(mx.__file__)]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    rows=[]; start=time.perf_counter()
    header={'protocol':'PROTOCOL.md + EXECUTION_FREEZE.md','seed':2026091107,
        'action_seed':2026091108,'source_hashes':hashes,
        'environment':{'python':platform.python_version(),'numpy':np.__version__,
                       'scipy':scipy.__version__,'platform':platform.platform()},
        'scope':'12 new diagnostic scenes; no population claim or certified continuum bound',
        'primary_success':'both absolute material errors <=0.1',
        'base_complex_observations':144,'candidate_complex_observations':192,
        'reference_complex_sigma':.01,'gain_annulus':[.75,1.25],
        'dda_cells':len(dda.points),'fine_dda_cells':len(fine.points),
        'source_boundary_check':boundary_residual(3+.05j,.63,5)}
    def save(complete=False):
        data={**header,'complete':complete,'rows':rows,'wall_seconds':time.perf_counter()-start,
              'peak_rss_kib_linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        temp=target.with_suffix('.partial'); temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(target)
    for scene in range(12):
        truth=np.r_[rng.uniform([1.5,2],[4,5]),rng.uniform(-2,2)]
        gain=rng.uniform(.75,1.25)*np.exp(1j*rng.uniform(-np.pi,np.pi))
        shifted=points+np.array([truth[2]*.001,0,0])
        tick=time.perf_counter()
        currents=dda.currents(truth[:2]+1j*LOSS)
        exact=(mx.dipole_kernel(shifted,dda.points,K)@currents).ravel()
        sigma=.01*np.linalg.norm(gain*exact[:base_count])/np.sqrt(base_count)
        y=gain*exact+sigma/np.sqrt(2)*(rng.normal(size=exact.size)+1j*rng.normal(size=exact.size))
        reference=gain+.01/np.sqrt(2)*(rng.normal()+1j*rng.normal())
        random_index=int(action_rng.integers(192))
        row={'scene':scene,'truth':truth.tolist(),'gain':encode(gain),'sigma':sigma,
             'reference':encode(reference),'fixed_index':0,'random_index':random_index,
             'noise_free_all_fields':encode(exact),'observed_all_fields':encode(y),
             'generation_seconds':time.perf_counter()-tick,'methods':{}}
        if scene<4:
            lo=forward(truth)
            f4,f5=[m.field(truth[:2]+1j*LOSS,shifted).ravel() for m in high]
            fq=quad.field(truth[:2]+1j*LOSS,shifted).ravel()
            ff=fine.field(truth[:2]+1j*LOSS,shifted).ravel()
            rel=lambda a,b:float(np.linalg.norm(a-b)/np.linalg.norm(b))
            row['convergence_diagnostics']={
                'L3_vs_L4':rel(lo,f4),'L4_vs_L5':rel(f4,f5),
                'quadrature_16x32_vs_24x48':rel(lo,fq),
                'DDA_011_vs_L5':rel(exact,f5),'DDA_009_vs_L5':rel(ff,f5),
                'DDA_011_vs_009':rel(exact,ff),
                'interpretation':'diagnostic differences, NOT rigorous error bounds'}
        configs=[('no_reference',None,None),('noisy_reference_GLS',None,(reference,.01)),
                 ('fixed_EM',0,None),('random_EM',random_index,None)]
        for name,extra,ref in configs:
            indices=np.r_[np.arange(base_count),base_count+extra] if extra is not None else np.arange(base_count)
            yy=y[indices]
            def residual(theta):return profile(yy,forward(theta)[indices],sigma,ref)[0]
            tick=time.perf_counter(); fits=[]
            for initial in STARTS:
                fit=least_squares(residual,initial,bounds=BOUNDS,max_nfev=80,
                                  ftol=1e-8,xtol=1e-8,gtol=1e-8)
                fits.append({'theta':fit.x.tolist(),'rss_real_white':float(fit.fun@fit.fun),
                    'nfev':fit.nfev,'njev':fit.njev,'status':int(fit.status),
                    'optimizer_success':bool(fit.success),'message':fit.message})
            best=min(fits,key=lambda v:v['rss_real_white']); theta=np.array(best['theta'])
            rb,estimated_gain=profile(yy,forward(theta)[indices],sigma,ref)
            error=np.abs(theta-truth)
            row['methods'][name]={'starts':fits,'theta':theta.tolist(),'gain':encode(estimated_gain),
                'material_errors':error[:2].tolist(),'geometry_error_mm':float(error[2]),
                'gain_relative_error':float(abs(estimated_gain-gain)/abs(gain)),
                'task_success':bool(np.all(error[:2]<=.1)),
                'rss_real_white':float(rb@rb),
                'sensor_prediction_error':float(np.linalg.norm(estimated_gain*forward(theta)-gain*exact)/np.linalg.norm(gain*exact)),
                'elapsed_seconds':time.perf_counter()-tick,'scientific_acceptance':'unresolved'}
        rows.append(row);save()
        print(json.dumps({'scene':scene,'success':{k:v['task_success'] for k,v in row['methods'].items()},
                          'errors':{k:v['material_errors'] for k,v in row['methods'].items()}}),flush=True)
    save(True)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output')
    main(parser.parse_args().output)
