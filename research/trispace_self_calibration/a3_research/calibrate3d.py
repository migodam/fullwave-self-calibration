"""Independent-reference 3D DEVELOPMENT calibration with exact material tangent.

Known supports and incident plane waves; only Rx translation is calibrated.
This is two-region permittivity estimation, NOT high-dimensional 3D imaging.
Per-illumination gains are shared across frequencies; delay is a common length.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares
from maxwell3d import DipoleVIE,dipole_kernel,treams_field,receivers

KS=np.array([3.,6.,9.])
CENTERS=[[-.2,0,0],[.2,.03,.04]]
RADII=[.12,.10]
LOSS=np.array([.03,.04])


class JointModel:
    def __init__(self,spacing=.04):
        self.models=[DipoleVIE(CENTERS,RADII,spacing,k,fill_quadrature=6) for k in KS]
        self.rx=receivers()
        self.cached=None
        self.n_predict=0
    def field_jac(self,z):
        # z: eps(2), receiver translation(3), delay-length(1), logamp(4), phase(4)
        material=tuple(z[:2])
        if self.cached is None or self.cached[0]!=material:
            self.cached=(material,[m.currents(np.array(material)+1j*LOSS,True) for m in self.models])
        yy=[];jj=[]
        for k,m,(p,dp) in zip(KS,self.models,self.cached[1]):
            pts=self.rx+z[2:5]
            op=dipole_kernel(pts,m.points,k)
            raw=(op@p).reshape(len(pts),3,4)
            jac=np.zeros(raw.shape+(14,),complex)
            jac[...,:2]=(op@dp.reshape(len(p),-1)).reshape(raw.shape+(2,))
            for d in range(3):
                shift=np.eye(3)[d]*1e-5
                dop=(dipole_kernel(pts+shift,m.points,k)-dipole_kernel(pts-shift,m.points,k))/(2e-5)
                jac[...,2+d]=(dop@p).reshape(raw.shape)
            factor=np.exp(z[6:10]+1j*z[10:14]+1j*k*z[5])
            y=raw*factor
            jac=jac*factor[None,None,:,None]
            jac[...,5]=1j*k*y
            for t in range(4):
                jac[:,:,t,6+t]=y[:,:,t]
                jac[:,:,t,10+t]=1j*y[:,:,t]
            yy.append(y);jj.append(jac)
        self.n_predict+=1
        return np.array(yy),np.array(jj)


def real(a):return np.r_[a.real.ravel(),a.imag.ravel()]


def truth(seed):
    rng=np.random.default_rng(seed)
    z=np.zeros(14)
    z[:2]=[2.2+.4*rng.random(),2.8+.6*rng.random()]
    z[2:5]=rng.normal(size=3);z[2:5]*=.09/np.linalg.norm(z[2:5])
    z[5]=.06
    z[6:10]=rng.normal(0,.06,4);z[10:14]=rng.normal(0,.12,4)
    raw=np.array([treams_field(CENTERS,RADII,z[:2]+1j*LOSS,k,receivers()+z[2:5],8 if k>9 else 5)[0] for k in KS])
    y0=raw*np.exp(z[6:10]+1j*z[10:14]+1j*KS[:,None,None,None]*z[5])
    sigma=np.linalg.norm(y0)/np.sqrt(y0.size)*10**(-30/20)
    y=y0+sigma/np.sqrt(2)*(rng.normal(size=y0.shape)+1j*rng.normal(size=y0.shape))
    return z,y,y0,sigma


def fit(model,y,sigma,ztrue,method,max_nfev=45,reference=None,reference_sigma=.02):
    base=np.zeros(14);base[:2]=[2.,2.5]
    if method=='oracle':base[2:5]=ztrue[2:5]
    indices={'unified':list(range(14)), 'geometry_only':list(range(5)),
             'unified_reference':list(range(14)),
             'fixed_wrong_pose':[0,1,5]+list(range(6,14)),
             'oracle':[0,1,5]+list(range(6,14))}[method]
    lo=np.array([1.1,1.1]+[-.3]*3+[-.3]+[-.5]*4+[-1.]*4)
    hi=np.array([5.,5.]+[.3]*3+[.3]+[.5]*4+[1.]*4)
    scale=np.array([1.,1.]+[.1]*3+[.1]+[.1]*8)
    cached={}
    def evaluate(v):
        key=tuple(v)
        if cached.get('key')!=key:
            z=base.copy();z[indices]=v
            pred,j=model.field_jac(z)
            r=np.sqrt(2)/sigma*real(pred-y)
            jc=j.reshape(-1,14)[:,indices]
            jr=np.sqrt(2)/sigma*np.r_[jc.real,jc.imag]
            if reference is not None:
                refpred=np.exp(z[6:10]+1j*z[10:14]+1j*KS[:,None]*z[5])
                refj=np.zeros((len(KS),4,14),complex)
                refj[...,5]=1j*KS[:,None]*refpred
                for t in range(4):
                    refj[:,t,6+t]=refpred[:,t]
                    refj[:,t,10+t]=1j*refpred[:,t]
                r=np.r_[r,np.sqrt(2)/reference_sigma*real(refpred-reference)]
                rj=refj.reshape(-1,14)[:,indices]
                jr=np.r_[jr,np.sqrt(2)/reference_sigma*np.r_[rj.real,rj.imag]]
            cached.update(key=key,r=r,j=jr,z=z)
        return cached['r'],cached['j']
    t0=time.perf_counter()
    before=[m.work.copy() for m in model.models]
    opt=least_squares(lambda v:evaluate(v)[0],base[indices],jac=lambda v:evaluate(v)[1],
        bounds=(lo[indices],hi[indices]),x_scale=scale[indices],max_nfev=max_nfev,
        ftol=1e-8,xtol=1e-8,gtol=1e-7)
    z=base.copy();z[indices]=opt.x
    pred,j=model.field_jac(z)
    # Held-out Rx angles: evaluate independent reference at estimated material/pose.
    held=receivers(17,1.6)
    ptrue=np.array([treams_field(CENTERS,RADII,ztrue[:2]+1j*LOSS,k,held+ztrue[2:5],8 if k>9 else 5)[0] for k in KS])
    pest=np.array([treams_field(CENTERS,RADII,z[:2]+1j*LOSS,k,held+z[2:5],8 if k>9 else 5)[0] for k in KS])
    phase=np.angle(pest*np.conj(ptrue))
    phase_rmse=np.sqrt(np.sum(abs(ptrue)**2*phase**2)/np.sum(abs(ptrue)**2))
    n=len(opt.fun)
    work={key:sum(m.work[key]-b[key] for m,b in zip(model.models,before)) for key in before[0]}
    # Local information diagnostics in declared parameter units, NOT global coverage.
    jfit=opt.jac
    uj,sj,vhj=np.linalg.svd(jfit,full_matrices=False)
    local_rank=int(np.count_nonzero(sj>1e-10*sj[0]))
    local_sd=np.sqrt(np.sum((vhj.T/sj)**2,axis=1)) if local_rank==len(indices) else np.full(len(indices),np.inf)
    return dict(method=method,estimated=z.tolist(),pose_error_m=float(np.linalg.norm(z[2:5]-ztrue[2:5])),
        material_relative_error=float(np.linalg.norm(z[:2]-ztrue[:2])/np.linalg.norm(ztrue[:2])),
        delay_error_m=float(abs(z[5]-ztrue[5])),reduced_chisquare=float(opt.fun@opt.fun/(n-len(indices))),
        heldout_phase_rmse_rad=float(phase_rmse),heldout_field_relative_error=float(np.linalg.norm(pest-ptrue)/np.linalg.norm(ptrue)),
        nfev=int(opt.nfev),status=int(opt.status),elapsed_seconds=time.perf_counter()-t0,work=work,
        fitted_indices=indices,local_rank=local_rank,local_linearized_sd=[float(s) if np.isfinite(s) else None for s in local_sd],reference_used=reference is not None,
        scope='development independent 3D reference; known two-region support and world-fixed illumination; no coverage certificate')


def main():
    global KS
    parser=argparse.ArgumentParser();parser.add_argument('--seeds',nargs='+',type=int,default=[6101,6102])
    parser.add_argument('--spacing',type=float,default=.04);parser.add_argument('--max-nfev',type=int,default=45)
    parser.add_argument('--reference-study',action='store_true')
    parser.add_argument('--frequencies',nargs='+',type=float,default=[3.,6.,9.])
    args=parser.parse_args();KS=np.array(args.frequencies);rows=[]
    for seed in args.seeds:
        z,y,y0,sigma=truth(seed)
        ref_rng=np.random.default_rng(seed+10000)
        reference=np.exp(z[6:10]+1j*z[10:14]+1j*KS[:,None]*z[5])
        reference+=.02/np.sqrt(2)*(ref_rng.normal(size=reference.shape)+1j*ref_rng.normal(size=reference.shape))
        methods=['unified','unified_reference'] if args.reference_study else ['fixed_wrong_pose','geometry_only','unified','oracle']
        for method in methods:
            model=JointModel(args.spacing)
            row=fit(model,y,sigma,z,method,args.max_nfev,reference=reference if method=='unified_reference' else None)
            row.update(seed=seed,spacing=args.spacing,frequencies=KS.tolist(),true=z.tolist(),noise_sigma=float(sigma))
            rows.append(row);print(json.dumps(row),flush=True)
            label='reference3d_development_' if args.reference_study else 'calibration3d_development_'
            suffix='' if args.frequencies==[3.,6.,9.] else '_k'+'-'.join(str(int(k)) for k in KS)
            path=Path(__file__).resolve().parent/'results'/(label+str(args.spacing)+suffix+'.json')
            path.write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':main()
