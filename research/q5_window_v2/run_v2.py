"""Frozen V2 independent-development stream; no selector and no pass retuning.

Run from any working directory. Results are immutable; --out chooses a fresh
folder. The primary target is both material absolute errors <= 0.1.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
from pathlib import Path
import sys,json,time,hashlib,resource,argparse,platform
import numpy as np
import scipy
from scipy.optimize import least_squares
from multipole import SphereCluster,ShiftReadout,CENTERS,RADII,LOSS,K
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
import maxwell3d as mx
STARTS=[[2.,3.,0.],[1.6,2.2,-1.],[3.8,4.8,1.]]


def pack(z):
    a=np.asarray(z).ravel();return np.c_[a.real,a.imag].tolist()


def real(z): return np.r_[np.real(z),np.imag(z)]
def rel(a,b): return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))


def profile(y,f,sigma,reference=None,reference_sigma=.01,phase_free=False):
    """Exact scalar annulus profiling; phase_free profiles reference phase too."""
    v=f/sigma;b=y/sigma
    if phase_free:
        if reference is None: raise ValueError('phase-free reference missing')
        cross=np.vdot(v,b);aa=(abs(cross)+abs(reference)/reference_sigma**2)/(np.vdot(v,v).real+1/reference_sigma**2)
        g=np.clip(aa,.75,1.25)*np.exp(1j*np.angle(cross))
        r=np.r_[real(b-g*v),(abs(reference)-abs(g))/reference_sigma]
    else:
        if reference is not None:
            v=np.r_[v,1/reference_sigma];b=np.r_[b,reference/reference_sigma]
        gg=np.vdot(v,b)/np.vdot(v,v).real
        g=np.clip(abs(gg),.75,1.25)*np.exp(1j*np.angle(gg))
        r=real(b-g*v)
    return np.sqrt(2)*r,complex(g)


def fit(y,readout,sigma,ref=None,phase_free=False,known_shift=None):
    fits=[];before=time.perf_counter()
    calls=0
    def theta(t): return np.r_[t,known_shift] if known_shift is not None else t
    def residual(t):
        nonlocal calls
        calls+=1
        return profile(y,readout(theta(t)),sigma,ref,phase_free=phase_free)[0]
    for init in STARTS:
        x0=init[:2] if known_shift is not None else init
        lower=[1.5,2.] if known_shift is not None else [1.5,2.,-2.]
        upper=[4.,5.] if known_shift is not None else [4.,5.,2.]
        opt=least_squares(residual,x0,bounds=(lower,upper),max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8)
        fits.append({'theta':theta(opt.x).tolist(),'rss':float(opt.fun@opt.fun),'nfev':int(opt.nfev),'status':int(opt.status),'optimizer_success':bool(opt.success)})
    best=min(fits,key=lambda v:v['rss']);t=np.array(best['theta'])
    r,g=profile(y,readout(t),sigma,ref,phase_free=phase_free)
    return {'starts':fits,'theta':t.tolist(),'gain':[g.real,g.imag], 'rss':float(r@r),'residual_real_dimension':int(r.size),
            'residual_evaluations_including_fd':calls,'fit_seconds':time.perf_counter()-before}


def decomposition(model,read,truth,gain,clean,sigma):
    """Order-dependent tangent diagnostic, not physical attribution or a bound."""
    f=read(truth);d=gain*(clean-f)
    scale=np.sqrt(2)/sigma
    gaincols=np.stack([real(f),real(1j*f)],axis=1)*scale
    jac=[]
    for axis,h in enumerate([1e-4,1e-4,1e-3]):
        plus=truth.copy();minus=truth.copy();plus[axis]+=h;minus[axis]-=h
        # Direct readout is valid outside a tiny edge step and avoids interpolation domain clipping.
        df=(model.forward(plus,read.rx)-model.forward(minus,read.rx))/(2*h)
        jac.append(real(gain*df)*scale)
    material=np.stack(jac[:2],axis=1);geometry=jac[2][:,None]
    blocks=[gaincols,geometry,material];residual=real(d)*scale
    initial=float(residual@residual);fractions=[];basis0=np.empty((residual.size,0))
    for block in blocks:
        block=block-basis0@(basis0.T@block)
        u,s,v=np.linalg.svd(block,full_matrices=False)
        keep=s>max(s[0]*1e-10,1e-12);u=u[:,keep]
        comp=u@(u.T@residual);fractions.append(float(comp@comp)/max(initial,1e-300));residual-=comp
        basis0=np.c_[basis0,u]
    fractions.append(float(residual@residual)/max(initial,1e-300))
    J=np.c_[material,geometry,gaincols];b=real(d)*scale
    bias_none=np.linalg.lstsq(J,b,rcond=None)[0]
    extra=np.zeros((2,5));extra[0,3]=np.sqrt(2)/.01;extra[1,4]=np.sqrt(2)/.01
    bias_ref=np.linalg.lstsq(np.r_[J,extra],np.r_[b,[0.,0.]],rcond=None)[0]
    return {'label':'ordered gain -> geometry -> material -> remainder; not unique causal attribution',
            'energy_fractions':fractions,'linearized_material_bias_no_ref':bias_none[:2].tolist(),
            'linearized_material_bias_ref':bias_ref[:2].tolist(),
            'true_world_forward_difference_relative':rel(clean,f)}


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    report=out/'recovery.json';rawpath=out/'observations.jsonl'
    if report.exists() or rawpath.exists(): raise FileExistsError('immutable experiment output already exists')
    before=time.perf_counter();rng=np.random.default_rng(2026091107);extra_rng=np.random.default_rng(2026091108)
    rx=mx.receivers(12,.6);pool=mx.receivers(16,.6);struct=mx.receivers(36,.15)
    model=SphereCluster(3,18,36);read=ShiftReadout(model,rx);read_extra=ShiftReadout(model,pool)
    structural_readout=model.readout(struct)
    dda=mx.DipoleVIE(CENTERS,RADII,.011,K,fill_quadrature=4)
    refined={h:mx.DipoleVIE(CENTERS,RADII,h,K,fill_quadrature=4) for h in (.009,.007)}
    m5=SphereCluster(5,18,36);q3=SphereCluster(3,26,52)
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'multipole.py',Path(mx.__file__)]}
    result={'complete':False,'registered_seed':2026091107,'independent_extra_seed':2026091108,'source_sha256':hashes,
            'scope':'V2 independent development; not final test, not continuum truth, no new-policy claim',
            'primary_methods':['no_reference','noisy_reference_GLS','fixed_EM_scalar','random_EM_scalar'],
            'secondary_ablations':['phase_discard_reference','known_shift_reference'],
            'success_rule':'both material errors <= 0.1; geometry/gain reported separately',
            'noise':'proper CN(0,sigma^2 I), sqrt(2) real whitening, reference proper CN variance 0.01^2',
            'dda_cells':len(dda.points),'refinement_cells':{str(h):len(v.points) for h,v in refined.items()},
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'cpu_count':os.cpu_count(),'blas_threads':1},'rows':[]}
    def save():
        result['wall_seconds']=time.perf_counter()-before
        result['peak_rss_MiB_linux']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        tmp=report.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(report)
    with rawpath.open('w') as raw:
        for scene in range(12):
            truth=np.r_[rng.uniform([1.5,2.],[4.,5.]),rng.uniform(-2.,2.)]
            gain=rng.uniform(.75,1.25)*np.exp(1j*rng.uniform(-np.pi,np.pi))
            start=time.perf_counter();p=dda.currents(truth[:2]+1j*LOSS)
            shift=np.array([truth[2]*.001,0,0])
            clean=(mx.dipole_kernel(rx+shift,dda.points,K)@p).ravel()
            candidate=(mx.dipole_kernel(pool+shift,dda.points,K)@p).ravel()
            fixed_structure=(mx.dipole_kernel(struct,dda.points,K)@p).ravel()
            sigma=.01*np.linalg.norm(gain*clean)/np.sqrt(clean.size)
            y=gain*clean+sigma/np.sqrt(2)*(rng.normal(size=clean.size)+1j*rng.normal(size=clean.size))
            reference=gain+.01/np.sqrt(2)*(rng.normal()+1j*rng.normal())
            index=int(extra_rng.integers(candidate.size))
            added={name:gain*candidate[idx]+sigma/np.sqrt(2)*(rng.normal()+1j*rng.normal()) for name,idx in [('fixed',0),('random',index)]}
            raw.write(json.dumps({'scene':scene,'truth':truth.tolist(),'gain':pack(gain),'sigma':sigma,'y':pack(y),
                'reference':pack(reference),'fixed_extra':pack(added['fixed']),'random_extra':pack(added['random']),'random_index':index},separators=(',',':'))+'\n');raw.flush()
            row={'scene':scene,'truth':truth.tolist(),'gain':[float(gain.real),float(gain.imag)],'sigma':float(sigma),
                 'generation_seconds':time.perf_counter()-start,'random_extra_index':index,'methods':{},'scientific_acceptance':'unresolved'}
            for name in result['primary_methods']+result['secondary_ablations']:
                ref=None;phase_free=False;known_shift=None;yy=y;rr=read
                if name in ('noisy_reference_GLS','phase_discard_reference','known_shift_reference'): ref=reference
                if name=='phase_discard_reference': phase_free=True
                if name=='known_shift_reference': known_shift=float(truth[2])
                if name in ('fixed_EM_scalar','random_EM_scalar'):
                    mode='fixed' if name.startswith('fixed') else 'random';idx=0 if mode=='fixed' else index
                    yy=np.r_[y,added[mode]]
                    rr=lambda t,idx=idx:np.r_[read(t),read_extra(t)[idx]]
                fit0=fit(yy,rr,sigma,ref,phase_free,known_shift)
                t=np.array(fit0['theta']);g=complex(*fit0['gain']);errors=np.abs(t-truth)
                pp=model.coefficients(*(t[:2]+1j*LOSS));pred_structure=(structural_readout@pp).ravel()
                fit0.update({'material_errors':errors[:2].tolist(),'material_success':bool(np.all(errors[:2]<=.1)),
                    'geometry_error_mm':float(errors[2]),'gain_relative_error':float(abs(g-gain)/abs(gain)),
                    'fixed_world_structural_error':rel(pred_structure,fixed_structure),'scientific_acceptance':'unresolved'})
                row['methods'][name]=fit0
            row['discrepancy_diagnostic']=decomposition(model,read,truth,gain,clean,sigma)
            if scene<4:
                fields={'.011':clean}
                for h,dd in refined.items():
                    pp=dd.currents(truth[:2]+1j*LOSS)
                    fields[str(h)]=(mx.dipole_kernel(rx+shift,dd.points,K)@pp).ravel()
                f5=m5.forward(truth,rx);f3=model.forward(truth,rx)
                row['convergence_diagnostics']={'scope':'differences, NOT certified error bounds',
                    'DDA_vs_L5':{h:rel(f,f5) for h,f in fields.items()},
                    'DDA_009_vs_007':rel(fields['0.009'],fields['0.007']),
                    'L3_vs_L5':rel(f3,f5),'quadrature18_vs26':rel(q3.forward(truth,rx),f3)}
            result['rows'].append(row);save()
            print(json.dumps({'scene':scene,'success':{n:m['material_success'] for n,m in row['methods'].items()},'seconds':result['wall_seconds']}),flush=True)
    result['complete']=True;result['observations_sha256']=hashlib.sha256(rawpath.read_bytes()).hexdigest()
    result['summary']={name:{'successes':sum(r['methods'][name]['material_success'] for r in result['rows']),
           'n':12,'median_material_errors':np.median([r['methods'][name]['material_errors'] for r in result['rows']],axis=0).tolist(),
           'median_geometry_error_mm':float(np.median([r['methods'][name]['geometry_error_mm'] for r in result['rows']])),
           'median_structural_error':float(np.median([r['methods'][name]['fixed_world_structural_error'] for r in result['rows']])),
           'total_fit_seconds':sum(r['methods'][name]['fit_seconds'] for r in result['rows'])} for name in result['primary_methods']+result['secondary_ablations']}
    save();print(json.dumps(result['summary'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=HERE/'results/v2');a=p.parse_args();run(a.out)
