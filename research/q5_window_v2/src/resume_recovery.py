"""Resume immutable observations after interruption, with identical hyperparameters."""
from recovery import *
import shutil

def main():
    target=HERE/'results/recovery_v2.json';d=json.loads(target.read_text())
    if d['complete']:raise RuntimeError('Already complete')
    snapshot=HERE/'results/recovery_interrupted.json'
    if not snapshot.exists():shutil.copyfile(target,snapshot)
    raw=np.load(HERE/'results/observations_v2.npz');start=time.perf_counter();prior=d['seconds']
    solver=SphereSolver();rx=receivers();valid=receivers(16);struct=receivers(36,.15)
    d['resume_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    d['interruption']='Tool wall-clock limit; saved data reused; no fit hyperparameter change'
    def save():
        d['seconds']=prior+time.perf_counter()-start
        d['max_rss_kib']=max(d['max_rss_kib'],resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(target)
    processed=0
    for row in d['rows']:
        s=row['scene'];truth=raw[f'{s}_truth'];gain=complex(raw[f'{s}_gain']);sigma=row['sigma']
        y=raw[f'{s}_y'];ey=raw[f'{s}_extra_y'];ref=complex(raw[f'{s}_reference']);structure=raw[f'{s}_structure']
        for name,idx,reference in [('no_reference',None,None),('noisy_reference_GLS',None,ref),('fixed_EM',0,None),('random_EM',row['random_index'],None)]:
            if name in row['methods']:continue
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
            t=np.array(min(fits,key=lambda f:f['rss'])['theta']);r,g=profile(yy,field(t),sigma,reference);err=np.abs(t-truth)
            st=solver.field(t[:2]+1j*LOSS,struct).ravel()
            row['methods'][name]={'starts':fits,'theta':t.tolist(),'gain':[g.real,g.imag],
                'material_errors':err[:2].tolist(),'task_success':bool(np.all(err[:2]<=.1)),
                'geometry_error_mm':float(err[2]),'gain_relative_error':float(abs(g-gain)/abs(gain)),
                'sensor_relative_residual':float(np.linalg.norm(y-g*solver.forward(t,rx))/np.linalg.norm(y)),
                'structural_field_error':float(np.linalg.norm(st-structure)/np.linalg.norm(structure)),
                'rss':float(r@r),'total_seconds':time.perf_counter()-before,'scientific_acceptance':'unresolved'}
            save();print(s,name,row['methods'][name]['task_success'],flush=True)
            processed+=1
            if processed>=8:return
    d['complete']=True;save()
if __name__=='__main__':main()
