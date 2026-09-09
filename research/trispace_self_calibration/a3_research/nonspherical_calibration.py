"""Known-support ellipsoid calibration with an independent ADDA reference.

DEVELOPMENT, not general 3D imaging. Reference N64 LDR; inverse N16/32
fill-weighted CM+RR. Three translations, one real epsilon, one common delay
and four frequency-shared complex illumination gains. No receiver gains here.
"""
import argparse
import gc
import json
import time
import numpy as np
from scipy.optimize import least_squares
from nonspherical3d import OUT, ShapeVIE, reference, independent_radiation
from maxwell3d import receivers, dipole_kernel

KS=np.array([6.,9.,18.]);LOSS=.03


class Model:
    def __init__(self,n,rx=None):
        self.models=[ShapeVIE('ellipsoid',n,k) for k in KS]
        self.rx=receivers() if rx is None else rx
        self.cache=None

    def field_jac(self,z):
        # epsilon, translation xyz, delay, four log amplitudes, four phases
        if self.cache is None or self.cache[0]!=z[0]:
            self.cache=(z[0],[m.currents([z[0]+1j*LOSS],True) for m in self.models])
        ys=[];js=[]
        for k,m,(p,dp) in zip(KS,self.models,self.cache[1]):
            rx=self.rx+z[1:4];op=dipole_kernel(rx,m.points,k)
            raw=(op@p).reshape(len(rx),3,4)
            j=np.zeros(raw.shape+(13,),complex)
            j[...,0]=(op@dp[:,:,0]).reshape(raw.shape)
            for d in range(3):
                h=np.eye(3)[d]*1e-5
                dop=(dipole_kernel(rx+h,m.points,k)-dipole_kernel(rx-h,m.points,k))/(2e-5)
                j[...,1+d]=(dop@p).reshape(raw.shape)
            factor=np.exp(z[5:9]+1j*z[9:13]+1j*k*z[4])
            y=raw*factor;j*=factor[None,None,:,None]
            j[...,4]=1j*k*y
            for t in range(4):j[:,:,t,5+t]=y[:,:,t];j[:,:,t,9+t]=1j*y[:,:,t]
            ys.append(y);js.append(j)
        return np.array(ys),np.array(js)


def reference_field(sources,z,rx):
    raw=np.array([independent_radiation(p,q,rx+z[1:4],k) for k,(p,q) in zip(KS,sources)])
    factor=np.exp(z[5:9]+1j*z[9:13]+1j*KS[:,None,None,None]*z[4])
    return raw*factor,raw


def metrics(a,b):
    return dict(relative_field_error=float(np.linalg.norm(a-b)/np.linalg.norm(b)),
                weighted_phase_rmse_rad=float(np.sqrt(np.sum(abs(b)**2*np.angle(a*b.conj())**2)/np.sum(abs(b)**2))))


def test():
    z=np.array([2.3,.03,-.02,.04,.03]+[.01]*4+[.02]*4)
    m=Model(8);y,j=m.field_jac(z);h=1e-5
    fd=[]
    for d in range(13):
        e=np.eye(13)[d]*h
        fd.append((m.field_jac(z+e)[0]-m.field_jac(z-e)[0])/(2*h))
    fd=np.stack(fd,-1)
    errors=[float(np.linalg.norm(fd[...,d]-j[...,d])/np.linalg.norm(fd[...,d])) for d in range(13)]
    assert max(errors)<1e-7
    row=dict(jacobian_column_errors=errors,passed=True,scope='13-parameter ellipsoid adapter finite-difference check')
    (OUT/'results/nonspherical_calibration_checks.json').write_text(json.dumps(row,indent=2)+'\n')
    print(json.dumps(row),flush=True)


def run(seed,ns,methods):
    rng=np.random.default_rng(seed)
    true=np.zeros(13);true[0]=2.5
    true[1:4]=rng.normal(size=3);true[1:4]*=.09/np.linalg.norm(true[1:4])
    true[4]=.06;true[5:9]=rng.normal(0,.06,4);true[9:13]=rng.normal(0,.12,4)
    sources=[];ref_seconds=0.
    for k in KS:
        p,q,meta=reference('ellipsoid',64,k,true[0]+1j*LOSS)
        sources.append((p,q));ref_seconds+=meta['reference_seconds']
    mean,_=reference_field(sources,true,receivers())
    sigma=np.linalg.norm(mean)/np.sqrt(mean.size)*10**(-30/20)
    y=mean+sigma/np.sqrt(2)*(rng.normal(size=mean.shape)+1j*rng.normal(size=mean.shape))
    ref0=np.exp(true[5:9]+1j*true[9:13]+1j*KS[:,None]*true[4])
    refnoise=np.random.default_rng(seed+10000)
    ref=ref0+.02/np.sqrt(2)*(refnoise.normal(size=ref0.shape)+1j*refnoise.normal(size=ref0.shape))
    dest=OUT/'results/nonspherical_calibration.json'
    rows=json.loads(dest.read_text()) if dest.exists() else []
    real=lambda a:np.r_[a.real.ravel(),a.imag.ravel()]
    lo=np.array([1.2]+[-.25]*3+[-.2]+[-.5]*4+[-1.]*4)
    hi=np.array([5.]+[.25]*3+[.2]+[.5]*4+[1.]*4)
    scale=np.array([1.]+[.1]*12)
    held=receivers(17,1.6);heldtrue,structtrue=reference_field(sources,true,held)
    for n in ns:
        for method in methods:
            if any(r['seed']==seed and r['n']==n and r['method']==method for r in rows):continue
            start=time.perf_counter();model=Model(n);base=np.zeros(13);base[0]=2.
            idx={'fixed_wrong_pose':[0,4]+list(range(5,13)),
                 'geometry_only':[0,1,2,3],
                 'unified':list(range(13)),
                 'unified_reference':list(range(13))}[method]
            cache={}
            def ev(v):
                if cache.get('key')!=v.tobytes():
                    z=base.copy();z[idx]=v
                    pred,j=model.field_jac(z);r=np.sqrt(2)/sigma*real(pred-y)
                    jc=j.reshape(-1,13)[:,idx];jr=np.sqrt(2)/sigma*np.r_[jc.real,jc.imag]
                    if method=='unified_reference':
                        rp=np.exp(z[5:9]+1j*z[9:13]+1j*KS[:,None]*z[4])
                        rj=np.zeros((3,4,13),complex);rj[...,4]=1j*KS[:,None]*rp
                        for t in range(4):rj[:,t,5+t]=rp[:,t];rj[:,t,9+t]=1j*rp[:,t]
                        r=np.r_[r,np.sqrt(2)/.02*real(rp-ref)]
                        jj=rj.reshape(-1,13)[:,idx];jr=np.r_[jr,np.sqrt(2)/.02*np.r_[jj.real,jj.imag]]
                    cache.update(key=v.tobytes(),r=r,j=jr)
                return cache['r'],cache['j']
            opt=least_squares(lambda v:ev(v)[0],base[idx],jac=lambda v:ev(v)[1],
                              bounds=(lo[idx],hi[idx]),x_scale=scale[idx],max_nfev=35,
                              ftol=1e-9,xtol=1e-9,gtol=1e-7)
            estimated=base.copy();estimated[idx]=opt.x
            fitting_seconds=time.perf_counter()-start
            model.rx=held;pred,_=model.field_jac(estimated)
            factor=np.exp(estimated[5:9]+1j*estimated[9:13]+1j*KS[:,None,None,None]*estimated[4])
            row=dict(seed=seed,n=n,spacing=.32/n,method=method,estimated=estimated.tolist(),true=true.tolist(),
                frequencies=KS.tolist(),optimizer_status=int(opt.status),nfev=int(opt.nfev),
                pose_error_m=float(np.linalg.norm(estimated[1:4]-true[1:4])),
                material_relative_error=float(abs(estimated[0]-true[0])/true[0]),
                delay_error_m=float(abs(estimated[4]-true[4])),
                reduced_chisquare=float(opt.fun@opt.fun/(len(opt.fun)-len(idx))),
                deployable_prediction=metrics(pred,heldtrue),reconstructed_structural_field=metrics(pred/factor,structtrue),
                fitting_seconds=fitting_seconds,total_seconds=time.perf_counter()-start,
                reference_generation_seconds=ref_seconds,work=[m.work for m in model.models],
                status='converged' if opt.status>0 else 'iteration_limit',
                scope='development one known-support homogeneous ellipsoid, independent N64 ADDA LDR; frequency-shared illumination gains; no free Rx gains, no general image')
            rows.append(row);dest.write_text(json.dumps(rows,indent=2)+'\n')
            print(json.dumps(row),flush=True)
            del model;gc.collect()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--test',action='store_true')
    p.add_argument('--seeds',nargs='+',type=int,default=[8101]);p.add_argument('--ns',nargs='+',type=int,default=[16,32])
    p.add_argument('--methods',nargs='+',default=['fixed_wrong_pose','geometry_only','unified','unified_reference'])
    args=p.parse_args()
    if args.test:test()
    else:
        for seed in args.seeds:run(seed,args.ns,args.methods)
