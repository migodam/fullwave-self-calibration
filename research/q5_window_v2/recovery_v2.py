"""Frozen V2 independent DDA diagnostics, not a new acquisition algorithm.

Run from the repository with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.
All observations and all starts are retained. Never overwrites an existing run.
"""
import argparse, hashlib, json, os, platform, resource, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from sphere_cluster import SphereCluster, CENTERS, RADII, LOSS, K
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
from maxwell3d import DipoleVIE, receivers, dipole_kernel

STARTS=[[2.,3.,0.],[1.6,2.2,-1.],[3.8,4.8,1.]]
LOW=np.array([1.5,2.,-2.]);HIGH=np.array([4.,5.,2.])

def complex_record(z):
    a=np.asarray(z)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}

def complex_read(d):
    return np.array(d['real'])+1j*np.array(d['imag'])

def realwhite(z):
    a=np.asarray(z).ravel()
    return np.sqrt(2)*np.r_[a.real,a.imag]

def profiled_gain(f,y,sigma,z=None,sigma_r=.01):
    """Exact minimizer on the closed gain annulus, including both boundaries."""
    denom=np.vdot(f,f).real/sigma**2
    num=np.vdot(f,y)/sigma**2
    if z is not None:
        denom+=1/sigma_r**2;num+=z/sigma_r**2
    if denom<=0:raise ValueError('No gain-sensitive data')
    unconstrained=num/denom
    gain=np.clip(abs(unconstrained),.75,1.25)*np.exp(1j*np.angle(unconstrained))
    return gain

def fit(model,points,y,sigma,z=None,candidate_points=None,extra_index=None,extra_y=None,
        gain_known=None,geometry_known=None):
    calls=0;clock=time.perf_counter();all_starts=[]
    def decode(t):
        return np.r_[t,geometry_known] if geometry_known is not None else np.asarray(t)
    def calc(t):
        nonlocal calls
        calls+=1;theta=decode(t)
        f=model.forward(theta,points);yy=y
        if extra_index is not None:
            extra=model.forward(theta,candidate_points)[extra_index]
            f=np.r_[f,extra];yy=np.r_[y,extra_y]
        g=gain_known if gain_known is not None else profiled_gain(f,yy,sigma,z)
        r=(yy-g*f)/sigma
        if z is not None:r=np.r_[r,(z-g)/.01]
        return realwhite(r),g
    for initial in STARTS:
        x0=initial[:2] if geometry_known is not None else initial
        bounds=(LOW[:2],HIGH[:2]) if geometry_known is not None else (LOW,HIGH)
        try:
            sol=least_squares(lambda t:calc(t)[0],x0,bounds=bounds,max_nfev=80,
                              ftol=1e-8,xtol=1e-8,gtol=1e-8)
            rr,gg=calc(sol.x)
            all_starts.append({'initial':x0,'theta':decode(sol.x).tolist(),'gain':complex_record(gg),
                'cost':float(rr@rr),'nfev':int(sol.nfev),'termination':int(sol.status),
                'optimality':float(sol.optimality),'solver_success':bool(sol.success)})
        except Exception as exc:
            all_starts.append({'initial':x0,'failure':repr(exc)})
    valid=[a for a in all_starts if 'cost' in a]
    if not valid:return {'failure':'all starts failed','starts':all_starts,'wall_seconds':time.perf_counter()-clock}
    best=min(valid,key=lambda s:s['cost'])
    return {**best,'starts':all_starts,'forward_evaluations':calls,'wall_seconds':time.perf_counter()-clock}

def annotate(fit_result,truth,g):
    if 'failure' in fit_result:return fit_result
    errors=np.abs(np.array(fit_result['theta'])-truth)
    fit_result.update(material_errors=errors[:2].tolist(),geometry_error_mm=float(errors[2]),
        gain_error=float(abs(complex_read(fit_result['gain'])-g)),
        material_success=bool(np.all(errors[:2]<=.1)))
    return fit_result

def numerical_derivative(model,theta,rx,j,step=1e-4):
    h=np.zeros(3);h[j]=step
    return (model.forward(theta+h,rx)-model.forward(theta-h,rx))/(2*step)

def orthogonal_diagnostic(model,theta,rx,f_truth,g):
    f=model.forward(theta,rx);e=realwhite(g*(f_truth-f))
    gg=np.c_[realwhite(f),realwhite(1j*f)]
    xx=np.c_[realwhite(g*numerical_derivative(model,theta,rx,2))]
    aa=np.stack([realwhite(g*numerical_derivative(model,theta,rx,j)) for j in [0,1]],axis=1)
    remaining=e.copy();basis=np.empty((len(e),0));parts={}
    for name,columns in [('gain',gg),('geometry_after_gain',xx),('material_after_gain_geometry',aa)]:
        residual=columns-basis@(basis.T@columns)
        u,s,v=np.linalg.svd(residual,full_matrices=False)
        q=u[:,s>max(s[0]*1e-10,1e-15)]
        projection=q@(q.T@remaining);parts[name]=float(np.linalg.norm(projection)**2)
        remaining-=projection;basis=np.c_[basis,q]
    parts['orthogonal_remainder']=float(remaining@remaining)
    total=float(e@e)
    return {'definition':'ordered tangent projection; not causal unique attribution',
            'relative_field_difference':float(np.linalg.norm(f_truth-f)/np.linalg.norm(f_truth)),
            'squared_norm_fractions':{key:val/total for key,val in parts.items()},
            'sum_error':float(abs(sum(parts.values())-total)/max(total,1e-300))}

def run(output):
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    start=time.perf_counter();model=SphereCluster(order=3,ntheta=12)
    checks=json.loads((HERE/'results/solver_checks.json').read_text())
    if checks['status']!='passed_numerical_checks':raise RuntimeError('Solver checks must pass before recovery')
    rx=receivers(12,.6);candidate_rx=receivers(16,.6)
    dda=DipoleVIE(CENTERS,RADII,.011,K,4)
    rng=np.random.default_rng(2026091107);choice_rng=np.random.default_rng(2026091108)
    def cn(n,s):return s/np.sqrt(2)*(rng.standard_normal(n)+1j*rng.standard_normal(n))
    raw=[]
    for index in range(12):
        theta=rng.uniform(LOW,HIGH);g=rng.uniform(.75,1.25)*np.exp(1j*rng.uniform(-np.pi,np.pi))
        shift=np.array([theta[2]*.001,0.,0.]);p=dda.currents(theta[:2]+1j*LOSS)
        field=(dipole_kernel(rx+shift,dda.points,K)@p).ravel()
        candidate=(dipole_kernel(candidate_rx+shift,dda.points,K)@p).ravel()
        sigma=.01*np.linalg.norm(g*field)/np.sqrt(len(field))
        y=g*field+cn(len(field),sigma);extra_y=g*candidate+cn(len(candidate),sigma)
        z=g+cn(1,.01)[0];random_index=int(choice_rng.integers(len(candidate)))
        raw.append({'scene':index,'truth':theta.tolist(),'gain':complex_record(g),'sigma':float(sigma),
            'reference':complex_record(z),'sigma_reference':.01,'y':complex_record(y),
            'noiseless_field':complex_record(field),'candidate_y':complex_record(extra_y),
            'candidate_noiseless_field':complex_record(candidate),'random_index':random_index,'fixed_index':0})
    (output/'raw_observations.json').write_text(json.dumps({'seed':2026091107,'choice_seed':2026091108,
        'rx':rx.tolist(),'candidate_rx':candidate_rx.tolist(),'scenes':raw},indent=2)+'\n')
    generation_seconds=time.perf_counter()-start
    results=[]
    for row in raw:
        i=row['scene'];theta=np.array(row['truth']);g=complex_read(row['gain']);y=complex_read(row['y'])
        field=complex_read(row['noiseless_field']);extra=complex_read(row['candidate_y']);z=complex_read(row['reference'])
        sigma=row['sigma'];res={'scene':i,'truth':theta.tolist(),'methods':{},'oracle_diagnostics':{}}
        for name,kwargs in [('no_reference',{}),('noisy_reference_GLS',{'z':z}),
            ('fixed_EM',{'candidate_points':candidate_rx,'extra_index':0,'extra_y':extra[0]}),
            ('random_EM',{'candidate_points':candidate_rx,'extra_index':row['random_index'],'extra_y':extra[row['random_index']]})]:
            fitted=fit(model,rx,y,sigma,**kwargs)
            res['methods'][name]=annotate(fitted,theta,g)
        for name,kwargs in [('noiseless_unknown_gain',{}),('noiseless_known_gain',{'gain_known':g}),
                            ('noiseless_known_gain_geometry',{'gain_known':g,'geometry_known':theta[2]})]:
            fitted=fit(model,rx,g*field,sigma,**kwargs)
            res['oracle_diagnostics'][name]=annotate(fitted,theta,g)
        res['structural_error']=orthogonal_diagnostic(model,theta,rx,field,g)
        results.append(res)
        (output/f'scene_{i:02d}.json').write_text(json.dumps(res,indent=2)+'\n')
        print('scene',i,{name:r.get('material_errors',r.get('failure')) for name,r in res['methods'].items()},flush=True)
    convergence=[]
    for spacing in [.009,.0075]:
        fine=DipoleVIE(CENTERS,RADII,spacing,K,4)
        for row in raw[:4]:
            theta=np.array(row['truth']);shift=np.array([theta[2]*.001,0,0])
            f_fine=fine.field(theta[:2]+1j*LOSS,rx+shift).ravel()
            fm=model.forward(theta,rx);coarse=complex_read(row['noiseless_field'])
            convergence.append({'scene':row['scene'],'spacing':spacing,'cells':len(fine.points),
                'relative_difference_to_coarse_DDA':float(np.linalg.norm(f_fine-coarse)/np.linalg.norm(f_fine)),
                'relative_difference_to_lmax3':float(np.linalg.norm(f_fine-fm)/np.linalg.norm(f_fine))})
        del fine
    order_checks=[]
    m4=SphereCluster(order=4,ntheta=14);m5=SphereCluster(order=5,ntheta=16)
    for row in raw[:4]:
        theta=np.array(row['truth']);f3=model.forward(theta,rx);f4=m4.forward(theta,rx);f5=m5.forward(theta,rx)
        order_checks.append({'scene':row['scene'],'l3_l4':float(np.linalg.norm(f3-f4)/np.linalg.norm(f4)),
            'l4_l5':float(np.linalg.norm(f4-f5)/np.linalg.norm(f5))})
    summary={}
    for name in results[0]['methods']:
        valid=[r['methods'][name] for r in results if 'failure' not in r['methods'][name]]
        summary[name]={'successes':sum(r['material_success'] for r in valid),'attempts':12,
            'median_material_errors':np.median([r['material_errors'] for r in valid],axis=0).tolist(),
            'median_geometry_error_mm':float(np.median([r['geometry_error_mm'] for r in valid])),
            'total_forward_evaluations':sum(r['forward_evaluations'] for r in valid),
            'total_wall_seconds':sum(r['wall_seconds'] for r in valid)}
    final={'kind':'registered_V2_independent_diagnostic_not_final_test','status':'executed',
        'summary':summary,'convergence_diagnostics_not_bounds':convergence,'order_checks_not_bounds':order_checks,
        'data_generation_seconds':generation_seconds,'wall_seconds':time.perf_counter()-start,
        'maxrss_kib_linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'python':platform.python_version(),'numpy':np.__version__,
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),HERE/'sphere_cluster.py',ROOT/'research/trispace_self_calibration/a3_research/maxwell3d.py']},
        'raw_sha256':hashlib.sha256((output/'raw_observations.json').read_bytes()).hexdigest(),
        'limitations':['No certified DDA continuum error','No class C covering finite bound','No hardware experiment',
                      'Equal extra scalar count is not equal electronics/EM hardware cost']}
    (output/'summary.json').write_text(json.dumps(final,indent=2)+'\n');print(json.dumps(final,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=HERE/'results/recovery_v2')
    args=p.parse_args();run(args.out)
