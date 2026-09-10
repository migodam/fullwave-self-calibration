"""Independent Q5 development reconstruction; original evidence untouched."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import sys
import json
import hashlib
import time
from functools import lru_cache
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
import treams

HERE=Path(__file__).resolve().parents[1]
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
import maxwell3d as mx
CENTERS=np.array([[-.06,0,0],[.055,.02,0]])
RADII=[.035,.025]
LOSS=np.array([.03,.05])
K=18.


@lru_cache(maxsize=24)
def scattered(e1,e2,order):
    spheres=[treams.TMatrix.sphere(order,K,r,[treams.Material(e),treams.Material()])
             for r,e in zip(RADII,np.array([e1,e2])+1j*LOSS)]
    tm=treams.TMatrix.cluster(spheres,CENTERS).interaction.solve()
    waves=[]
    for direction,pol in mx.illuminations():
        inc=treams.plane_wave(K*direction,pol.tolist(),k0=K,material=tm.material)
        waves.append(tm@inc.expand(tm.basis))
    return waves


def forward(theta,rx,order=3):
    points=rx+np.array([theta[2]*.001,0,0])
    return np.stack([np.asarray(w.efield(points)) for w in scattered(float(theta[0]),float(theta[1]),order)],axis=-1).ravel()


def profile(y,f,sigma,reference=None):
    v=f/sigma;b=y/sigma
    if reference is not None:
        v=np.r_[v,100.+0j];b=np.r_[b,reference*100.]
    g=np.vdot(v,b)/np.vdot(v,v)
    g=np.clip(abs(g),.75,1.25)*np.exp(1j*np.angle(g))
    r=b-g*v
    return np.r_[r.real,r.imag],g


def main():
    target=HERE/'results/development_cycle1.json'
    if target.exists(): raise RuntimeError('Immutable output already exists')
    paths=[Path(__file__),HERE/'docs/DEVELOPMENT_PROTOCOL.md',Path(mx.__file__)]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    rng=np.random.default_rng(2026091101)
    rx=mx.receivers(12,.6);valid=mx.receivers(16,.6);struct=mx.receivers(36,.15)
    start=time.perf_counter();rows=[]
    dda=mx.DipoleVIE(CENTERS,RADII,.011,K,fill_quadrature=4)
    def save(complete=False):
        data={'complete':complete,'source_hashes':hashes,'seed':2026091101,'rows':rows,
              'wall_seconds':time.perf_counter()-start,'dda_cells':len(dda.points),
              'scope':'Development; independent discretization not continuum truth; all certification unresolved.'}
        tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(target)
    for scene in range(24):
        truth=np.r_[rng.uniform([1.5,2],[4,5]),rng.uniform(-2,2)]
        gain=rng.uniform(.75,1.25)*np.exp(1j*rng.uniform(-np.pi,np.pi))
        before=time.perf_counter()
        if scene<12:
            yy=forward(truth,rx,4);vv=forward(truth,valid,4);ss=forward(truth,struct,4)
        else:
            p=dda.currents(truth[:2]+1j*LOSS)
            def df(points,shift=True):
                pts=points+np.array([truth[2]*.001,0,0]) if shift else points
                return (mx.dipole_kernel(pts,dda.points,K)@p).ravel()
            # Kernel flatten is (receiver,component,illumination), same as Treams.
            yy=df(rx);vv=df(valid);ss=df(struct)
        generation=time.perf_counter()-before
        sigma=.01*np.linalg.norm(gain*yy)/np.sqrt(yy.size)
        y=gain*yy+sigma/np.sqrt(2)*(rng.normal(size=yy.size)+1j*rng.normal(size=yy.size))
        reference=gain+.01/np.sqrt(2)*(rng.normal()+1j*rng.normal())
        row={'scene':scene,'data_class':'same_model' if scene<12 else 'independent_DDA',
             'truth':truth.tolist(),'gain':[gain.real,gain.imag],'sigma':sigma,
             'data_generation_seconds':generation,'methods':{}}
        for name,ref in [('no_reference',None),('joint_reference',reference)]:
            before=time.perf_counter();fits=[]
            for initial in ([2,3,0],[1.6,2.2,-1],[3.8,4.8,1]):
                fit=least_squares(lambda t:profile(y,forward(t,rx),sigma,ref)[0],initial,
                    bounds=([1.5,2,-2],[4,5,2]),max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8)
                fits.append({'theta':fit.x.tolist(),'rss':float(fit.fun@fit.fun),
                             'nfev':fit.nfev,'optimizer_success':bool(fit.success)})
            best=min(fits,key=lambda f:f['rss']);t=np.array(best['theta'])
            residual,g=profile(y,forward(t,rx),sigma,ref)
            errors=np.abs(t-truth); ge=abs(g-gain)/abs(gain)
            row['methods'][name]={'starts':fits,'theta':t.tolist(),'gain':[g.real,g.imag],
                'material_errors':errors[:2].tolist(),'geometry_error_mm':float(errors[2]),
                'gain_relative_error':float(ge),'task_success':bool(np.all(errors[:2]<=.1) and errors[2]<=.5 and ge<=.05),
                'sensor_prediction_error':float(np.linalg.norm(g*forward(t,valid)-gain*vv)/np.linalg.norm(gain*vv)),
                'structural_field_error':float(np.linalg.norm(forward(t,struct)-ss)/np.linalg.norm(ss)),
                'rss':float(residual@residual),'heuristic_residual_pass':bool(residual@residual<=chi2.ppf(.99,len(residual)-5)),
                'scientific_acceptance':'unresolved','total_seconds':time.perf_counter()-before}
        rows.append(row);save()
        print(json.dumps({'scene':scene,'class':row['data_class'],'success':{n:m['task_success'] for n,m in row['methods'].items()}}),flush=True)
    save(True)


if __name__=='__main__':main()
