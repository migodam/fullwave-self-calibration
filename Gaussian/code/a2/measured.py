"""Fixed-geometry Fresnel dielTM_dec8f Gaussian material fit (bounded)."""
from __future__ import annotations
import hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator,gmres
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import C0,Geometry,VIE,point_green

FREQS=(2,4,6); VIEW_IDS=np.arange(0,36,3); TRAIN_RX_PARITY=0
OUT=ROOT/'runs/a2/measured';FIG=ROOT/'figures/a2';OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)

def gauss_model(z,p,count):
    z=np.asarray(z); cols=[]; chi=0.
    for q in z.reshape(count,6):
        a,cx,cy,lsx,lsy,ang=q; sx,sy=np.exp(lsx),np.exp(lsy);d=p-np.array([cx,cy]);c,s=np.cos(ang),np.sin(ang);u=c*d[:,0]+s*d[:,1];v=-s*d[:,0]+c*d[:,1];g=np.exp(-.5*((u/sx)**2+(v/sy)**2));t=a*g;chi=chi+t
        cols += [g,t*(u*c/sx**2-v*s/sy**2),t*(u*s/sx**2+v*c/sy**2),t*u*u/sx**2,t*v*v/sy**2,t*u*v*(1/sy**2-1/sx**2)]
    return np.asarray(chi),np.column_stack(cols)

def disk_model(z,p):
    a,cx,cy,lr=z;rad=np.exp(lr);d=p-np.array([cx,cy]);rr=np.hypot(d[:,0],d[:,1]);edge=.001
    w=1/(1+np.exp(np.clip((rr-rad)/edge,-60,60)));chi=a*w;dw=w*(1-w)/edge
    # softened rasterization differentiates a nominally sharp disk only.
    return chi,np.c_[w,a*dw*d[:,0]/np.maximum(rr,1e-12),a*dw*d[:,1]/np.maximum(rr,1e-12),a*dw*rad]

class MeasuredProblem:
    def __init__(self):
        raw=np.loadtxt(ROOT/'data/dielTM_dec8f.exp');self.sha=hashlib.sha256((ROOT/'data/dielTM_dec8f.exp').read_bytes()).hexdigest();self.blocks=[]
        src=.720*np.c_[np.cos(np.arange(36)*np.pi/18),np.sin(np.arange(36)*np.pi/18)];rx=.760*np.c_[np.cos(np.arange(72)*np.pi/36),np.sin(np.arange(72)*np.pi/36)]
        for fg in FREQS:
            r=raw[raw[:,2]==fg]; vi=r[:,0].astype(int)-1;ri=r[:,1].astype(int)-1;k=2*np.pi*fg*1e9/C0
            # Established convention: conjugate records to exp(-iwt), then
            # estimate one complex point-source factor per view from back arc.
            inc=np.conj(r[:,5]+1j*r[:,6]);y=np.conj(r[:,3]+1j*r[:,4])-inc;direct=point_green(rx,src,k);coef=np.zeros(36,complex)
            for v in range(36):
                back=(vi==v)&(np.cos((ri*5-v*10)*np.pi/180)<-np.cos(np.pi/12));coef[v]=np.vdot(direct[ri[back],v],inc[back])/max(np.vdot(direct[ri[back],v],direct[ri[back],v]).real,1e-30)
            keep=np.isin(vi,VIEW_IDS); Y=np.full((72,len(VIEW_IDS)),np.nan+0j);valid=np.zeros_like(Y,dtype=bool)
            pos={v:j for j,v in enumerate(VIEW_IDS)}
            for rr,vv,yy in zip(ri[keep],vi[keep],y[keep]):Y[rr,pos[vv]]=yy;valid[rr,pos[vv]]=True
            train=valid & ((np.arange(72)[:,None]%2)==TRAIN_RX_PARITY)
            v=VIE(Geometry(n=64,side=.12,n_tx=len(VIEW_IDS),n_rx=72),frequency_hz=fg*1e9)
            v.tx=src[VIEW_IDS];v.rx=rx;v.E=point_green(v.points,v.tx,k)*coef[VIEW_IDS][None,:];v.S=k*k*v.h*v.h*point_green(v.rx,v.points,k);v.direct=point_green(v.rx,v.tx,k)
            self.blocks.append(dict(frequency_GHz=fg,vie=v,Y=Y,valid=valid,train=train,source_factor=coef[VIEW_IDS]))
    def evaluate(self,z,kind,gradient=True,profile_gain=True):
        total=0.; grad=None; detail=[]; images=[]; its=[]
        for b in self.blocks:
            v=b['vie'];chi,P=(gauss_model(z,v.points,kind) if kind!='disk' else disk_model(z,v.points));j,meta=v.solve(chi,v.E,rtol=3e-7,maxiter=120);pred=v.S@j;mask=b['train'];yy=b['Y'][mask];pp=pred[mask];gain=np.vdot(pp,yy)/max(np.vdot(pp,pp).real,1e-30) if profile_gain else 1+0j;res=gain*pred-b['Y'];den=max(np.vdot(yy,yy).real,1e-30);total+=np.vdot(res[mask],res[mask]).real/den;its.extend(meta['iterations']);images.append(chi)
            hmask=b['valid']&~mask;detail.append(dict(frequency_GHz=b['frequency_GHz'],gain=[gain.real,gain.imag],train_scattered_relative=float(np.linalg.norm(res[mask])/np.linalg.norm(b['Y'][mask])),held_receiver_scattered_relative=float(np.linalg.norm(res[hmask])/np.linalg.norm(b['Y'][hmask]))))
            if gradient:
                e=np.zeros_like(pred);e[mask]=np.conj(gain)*res[mask]/den
                A=v.system(chi);adj=np.empty_like(j)
                AH=LinearOperator(A.shape,matvec=A.rmatvec,dtype=complex)
                for q in range(j.shape[1]):
                    adj[:,q],info=gmres(AH,v.S.conj().T@e[:,q],rtol=3e-7,atol=0,maxiter=120)
                    if info: raise RuntimeError('adjoint GMRES failed')
                gc=2*np.real(np.sum(np.conj(adj)*(v.E+v.D.matmat(j)),axis=1));gg=P.T@gc;grad=gg if grad is None else grad+gg
        return (float(total),grad,detail,images,its) if gradient else (float(total),detail,images,its)

def fit(problem,count,start,kind='gaussian',maxiter=45,profile_gain=True):
    if kind=='disk':bounds=[(.05,6),(-.035,.035),(-.035,.035),(np.log(.008),np.log(.03))]
    else: bounds=sum(([ (.02,6),(-.04,.04),(-.04,.04),(np.log(.006),np.log(.04)),(np.log(.006),np.log(.04)),(-np.pi,np.pi)] for _ in range(count)),[])
    calls=0
    def fn(z):
        nonlocal calls;calls+=1;v,g,*_=problem.evaluate(z,'disk' if kind=='disk' else count,True,profile_gain);return v,g
    t=time.perf_counter();o=minimize(fn,start,jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':maxiter,'ftol':1e-8,'gtol':2e-5,'maxls':12});val,detail,images,its=problem.evaluate(o.x,'disk' if kind=='disk' else count,False,profile_gain)
    return dict(parameters=o.x.tolist(),objective=val,success=bool(o.success),message=str(o.message),iterations=int(o.nit),evaluations=calls,wall_seconds=time.perf_counter()-t,gmres_iterations=its,metrics=detail,image=images[0].reshape(64,64))

def main():
    p=MeasuredProblem();rng=np.random.default_rng(20260911);starts=[]
    for _ in range(10): starts.append(np.array([.7+rng.uniform(-.25,.25),rng.uniform(-.018,.018),rng.uniform(-.018,.018),np.log(rng.uniform(.010,.025)),np.log(rng.uniform(.010,.025)),rng.uniform(-.5,.5)]))
    results={'acquisition':{'frequencies_GHz':list(FREQS),'source_views_zero_based':VIEW_IDS.tolist(),'train_receivers':'even global receiver indices; odd held out','geometry':'fixed source radius 0.720 m, receiver radius 0.760 m, 36/72 published angular locations; no geometry calibration'},'data_sha256':p.sha,'models':{},'scope':'One published 2-D TM cylinder acquisition. Per-frequency complex gain profiling handles source normalization only; no hardware calibration, uncertainty model, 3-D claim, or diverse-target claim.'}
    for count in (1,2,4):
        fits=[]
        for s in starts:
            z=np.tile(s,count); # stagger centers deterministically
            z=z.reshape(count,6);z[:,1]+=np.linspace(-.008,.008,count);z=z.ravel();fits.append(fit(p,count,z,maxiter=45))
        best=min(fits,key=lambda x:x['objective']);results['models'][f'{count}_gaussian']={k:v for k,v in best.items() if k!='image'};results['models'][f'{count}_gaussian']['paired_start_objectives']=[x['objective'] for x in fits];results['models'][f'{count}_gaussian']['successful_starts']=sum(x['success'] for x in fits);results['models'][f'{count}_gaussian']['image']=best['image']
    diskstarts=[np.array([2.0,0,0,np.log(.015)]),np.array([2.5,.01,-.01,np.log(.014)]),np.array([3,-.01,.01,np.log(.016)])]
    diskfits=[fit(p,0,s,'disk',45) for s in diskstarts];db=min(diskfits,key=lambda x:x['objective']);results['known_shape_disk_oracle']={k:v for k,v in db.items() if k!='image'};results['known_shape_disk_oracle']['image']=db['image'];results['known_shape_disk_oracle']['warning']='Known-shape softened disk diagnostic only; not a fair general inverse baseline.'
    # Boundary-floor diagnostic relative to fitted disk proxy: distinguish its
    # narrow boundary band from interior/exterior mismatch on the same grid.
    ref=db['image'];x=np.linspace(-.06+.12/128,.06-.12/128,64);xx,yy=np.meshgrid(x,x,indexing='ij');rad=np.hypot(xx-db['parameters'][1],yy-db['parameters'][2]);R=np.exp(db['parameters'][3]);band=np.abs(rad-R)<=.003
    floor={}
    for q in (1,2,4):
        im=results['models'][f'{q}_gaussian']['image'];floor[str(q)]={'relative_to_disk':float(np.linalg.norm(im-ref)/max(np.linalg.norm(ref),1e-30)),'boundary_band_rmse':float(np.sqrt(np.mean((im[band]-ref[band])**2))),'nonboundary_rmse':float(np.sqrt(np.mean((im[~band]-ref[~band])**2)))}
    results['gaussian_boundary_floor_against_known_shape_disk']=floor
    imgs=[results['models'][f'{q}_gaussian'].pop('image') for q in (1,2,4)]+[results['known_shape_disk_oracle'].pop('image')]
    vmax=max(float(np.max(z)) for z in imgs);fig,ax=plt.subplots(1,4,figsize=(13,3.3),constrained_layout=True)
    for a,im,title in zip(ax,imgs,['1 Gaussian','2 Gaussians','4 Gaussians','known-shape disk oracle']):
        h=a.imshow(im.T,origin='lower',extent=(-60,60,-60,60),vmin=0,vmax=vmax,cmap='viridis');a.set(title=title,xlabel='x (mm)',ylabel='y (mm)')
    fig.colorbar(h,ax=ax,label='fitted real scalar contrast; common color scale');fig.savefig(FIG/'measured_reconstructions.png',dpi=180)
    (OUT/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps({k:{'objective':v['objective'],'held':[m['held_receiver_scattered_relative'] for m in v['metrics']]} for k,v in results['models'].items()},indent=2))
if __name__=='__main__':main()
