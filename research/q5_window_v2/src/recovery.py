"""Pre-registered V2 diagnostic stream; preserve raw observations and failures."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,time,json,hashlib,resource,platform
import numpy as np
from scipy.optimize import least_squares
import scipy
from sphere_solver import SphereSolver,receivers,CENTERS,RADII,LOSS,K
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
import maxwell3d as mx
HERE=Path(__file__).resolve().parents[1]


def profile(y,f,sigma,ref=None):
    v=f/sigma;b=y/sigma
    if ref is not None:v=np.r_[v,100.+0j];b=np.r_[b,100*ref]
    q=np.vdot(v,v).real
    if not q>0:raise ValueError('zero template')
    c=np.vdot(v,b);a=np.clip(abs(c)/q,.75,1.25)
    g=a*np.exp(1j*np.angle(c));r=b-g*v
    return np.sqrt(2)*np.r_[r.real,r.imag],g


def main():
    target=HERE/'results/recovery_v2.json'
    if target.exists():raise FileExistsError(target)
    if not (HERE/'results/validation.json').exists():raise RuntimeError('Run validation first')
    start=time.perf_counter();rng=np.random.default_rng(2026091107)
    indices=np.random.default_rng(2026091108).integers(192,size=12)
    rx=receivers();valid=receivers(16);struct=receivers(36,.15)
    solver=SphereSolver();dda=mx.DipoleVIE(CENTERS,RADII,.011,K,fill_quadrature=4)
    fine=mx.DipoleVIE(CENTERS,RADII,.009,K,fill_quadrature=4)
    L5=SphereSolver(5,28)
    sources=[Path(__file__),HERE/'src/sphere_solver.py',HERE/'IMPLEMENTATION_FREEZE.md',Path(mx.__file__)]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    rows=[];raw={};inputs=[]
    # Generate every observation before any reconstruction; no RNG-dependent fit path.
    for scene in range(12):
        t=np.r_[rng.uniform([1.5,2],[4,5]),rng.uniform(-2,2)]
        g=rng.uniform(.75,1.25)*np.exp(1j*rng.uniform(-np.pi,np.pi))
        before=time.perf_counter();p=dda.currents(t[:2]+1j*LOSS)
        def df(points,p=p,grid=dda):
            return (mx.dipole_kernel(points+[t[2]*.001,0,0],grid.points,K)@p).ravel()
        f=df(rx);extra=df(valid);structure=(mx.dipole_kernel(struct,dda.points,K)@p).ravel()
        sig=.01*np.linalg.norm(g*f)/np.sqrt(f.size)
        noise=lambda n:sig/np.sqrt(2)*(rng.normal(size=n)+1j*rng.normal(size=n))
        obs=g*f+noise(f.size);exobs=g*extra+noise(extra.size)
        ref=g+.01/np.sqrt(2)*(rng.normal()+1j*rng.normal())
        row={'scene':scene,'truth':t.tolist(),'gain':[g.real,g.imag],'sigma':float(sig),
             'random_index':int(indices[scene]),'methods':{},'generation_seconds':time.perf_counter()-before,
             'fixed_world_relative_forward_disagreement':float(np.linalg.norm(solver.forward(t,rx)-f)/np.linalg.norm(f))}
        if scene<4:
            pf=fine.currents(t[:2]+1j*LOSS);ff=df(rx,pf,fine)
            row['discretization_diagnostics']={'coarse_fine_DDA_relative':float(np.linalg.norm(f-ff)/np.linalg.norm(ff)),
                'fine_DDA_L5_relative':float(np.linalg.norm(ff-L5.forward(t,rx))/np.linalg.norm(L5.forward(t,rx))),
                'coarse_DDA_L5_relative':float(np.linalg.norm(f-L5.forward(t,rx))/np.linalg.norm(L5.forward(t,rx))),
                'label':'diagnostic differences, NOT error upper bounds'}
        inputs.append((t,g,sig,obs,exobs,ref,structure));rows.append(row)
        for n,v in [('truth',t),('gain',g),('y',obs),('extra_y',exobs),('reference',ref),('noiseless',f),('structure',structure)]:
            raw[f'{scene}_{n}']=v
    rawfile=HERE/'results/observations_v2.npz'
    if rawfile.exists():raise FileExistsError(rawfile)
    np.savez_compressed(rawfile,**raw)
    def save(done=False):
        output={'complete':done,'seed':2026091107,'random_index_seed':2026091108,'source_hashes':hashes,
            'observations_sha256':hashlib.sha256(rawfile.read_bytes()).hexdigest(),
            'numpy':np.__version__,'scipy':scipy.__version__,'python':platform.python_version(),
            'dda_cells':len(dda.points),'fine_dda_cells':len(fine.points),'rows':rows,
            'seconds':time.perf_counter()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'independent numerical development diagnostics; NOT continuum truth or final test',
            'continuum_certificate':'unresolved: no certified continuum VIE residual bound'}
        tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(output,indent=2)+'\n');tmp.replace(target)
    save()
    for scene,(truth,gain,sigma,y,ey,ref,structure) in enumerate(inputs):
        for name,idx,reference in [('no_reference',None,None),('noisy_reference_GLS',None,ref),('fixed_EM',0,None),('random_EM',int(indices[scene]),None)]:
            before=time.perf_counter();yy=y if idx is None else np.r_[y,ey[idx]]
            def field(t):
                base=solver.forward(t,rx)
                return base if idx is None else np.r_[base,solver.forward(t,valid)[idx]]
            fits=[]
            for initial in ([2,3,0],[1.6,2.2,-1],[3.8,4.8,1]):
                fit=least_squares(lambda t:profile(yy,field(t),sigma,reference)[0],initial,
                    bounds=([1.5,2,-2],[4,5,2]),max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8)
                fits.append({'theta':fit.x.tolist(),'rss':float(fit.fun@fit.fun),'nfev':fit.nfev,
                    'optimizer_success':bool(fit.success),'initial':initial})
            best=min(fits,key=lambda f:f['rss']);t=np.array(best['theta'])
            r,g=profile(yy,field(t),sigma,reference);errors=np.abs(t-truth)
            # Evaluate structure at FIXED world positions, never fitted receiver positions.
            st=solver.field(t[:2]+1j*LOSS,struct).ravel()
            rows[scene]['methods'][name]={'starts':fits,'theta':t.tolist(),'gain':[g.real,g.imag],
                'material_errors':errors[:2].tolist(),'task_success':bool(np.all(errors[:2]<=.1)),
                'geometry_error_mm':float(errors[2]),'gain_relative_error':float(abs(g-gain)/abs(gain)),
                'sensor_relative_residual':float(np.linalg.norm(y-g*solver.forward(t,rx))/np.linalg.norm(y)),
                'structural_field_error':float(np.linalg.norm(st-structure)/np.linalg.norm(structure)),
                'rss':float(r@r),'total_seconds':time.perf_counter()-before,'scientific_acceptance':'unresolved'}
            save()
        print(scene,{n:m['task_success'] for n,m in rows[scene]['methods'].items()},flush=True)
    save(True)

if __name__=='__main__':main()
