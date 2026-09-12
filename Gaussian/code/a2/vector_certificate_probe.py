"""Bounded 3-D vector dipole certificate check; no continuum Maxwell claim."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.certificate import certify

def dyadic(points,sources,k):
    dr=points[:,None,:]-sources[None,:,:];r=np.linalg.norm(dr,axis=-1)
    safe=np.where(r>0,r,1.);u=dr/safe[:,:,None];uu=u[:,:,:,None]*u[:,:,None,:]
    I=np.eye(3);g=np.exp(1j*k*safe)/(4*np.pi*safe)
    tensor=g[:,:,None,None]*(k*k*(I-uu)+(1j*k/safe-1/safe**2)[:,:,None,None]*(I-3*uu))
    tensor[r==0]=0
    return tensor.transpose(0,2,1,3).reshape(3*len(points),3*len(sources))

def main():
    rows=[]
    for n in (3,4):
        h=.21/n;axis=(np.arange(n)+.5)*h-.105
        points=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1).reshape(-1,3)
        for k in (4.,8.):
            K=dyadic(points,points,k);beta=k**3/(6*np.pi)
            gram=(K-K.conj().T)/(2j)+beta*np.eye(len(K))
            gram_min=float(np.linalg.eigvalsh(gram)[0])
            for kind in ('Gaussian','two_ellipsoids'):
                if kind=='Gaussian':base=.6*np.exp(-np.sum((points-np.array([.02,0.,0.]))**2,axis=1)/(2*.065**2))+.02
                else:
                    left=np.sum(((points-np.array([-.045,0.,0.]))/np.array([.04,.07,.05]))**2,axis=1)<1
                    right=np.sum(((points-np.array([.045,.015,0.]))/np.array([.03,.05,.06]))**2,axis=1)<1
                    base=.02+.6*left+.4*right
                chi=base*(1+.3j);vol=h**3
                inverse_alpha=1/(3*vol)+1/(vol*chi)-1j*beta
                L=np.diag(np.repeat(inverse_alpha,3))-K
                gamma=np.repeat(chi.imag/(vol*np.abs(chi)**2),3)
                diss_min=float(np.linalg.eigvalsh(-(L-L.conj().T)/(2j)-np.diag(gamma))[0])
                illuminations=[]
                for direction,pol in [(np.array([0,0,1]),np.array([1,0,0])),(np.array([0,0,1]),np.array([0,1,0])),(np.array([1,0,0]),np.array([0,0,1])),(np.array([1,0,0]),np.array([0,1,0]))]:
                    illuminations.append((np.exp(1j*k*(points@direction))[:,None]*pol).ravel())
                E=np.column_stack(illuminations)
                phi=np.arange(8)*np.pi/4;rx=np.c_[.5*np.cos(phi),.5*np.sin(phi),np.full(8,.2)]
                C=dyadic(rx,points,k);exact=np.linalg.solve(L,E);F=C@exact
                # Polynomial state range: no exact current used to construct it.
                snapshots=[];v=E.copy()
                for _ in range(8):snapshots.append(v);v=K@v/max(np.linalg.norm(K,ord=np.inf),1.)
                basis=np.linalg.qr(np.column_stack(snapshots),mode='reduced')[0]
                for rank in (8,16):
                    Q=basis[:,:rank];J=Q@np.linalg.solve(Q.conj().T@L@Q,Q.conj().T@E)
                    V=np.linalg.qr(np.column_stack((Q,C.conj().T)),mode='reduced')[0]
                    Z=V@np.linalg.solve(V.conj().T@L.conj().T@V,V.conj().T@C.conj().T)
                    cert=certify(L,E,C,J,Z,gamma)
                    error=float(np.linalg.norm(F-cert.corrected_field))
                    assert gram_min>-1e-10 and diss_min>-1e-9 and error<=cert.corrected_bound+1e-10
                    rows.append(dict(n=n,vector_unknowns=len(K),k=k,material=kind,rank=rank,
                                     radiation_gram_min=gram_min,dissipation_margin=diss_min,
                                     gamma_min=float(gamma.min()),corrected_error_relative=error/np.linalg.norm(F),
                                     bound_relative=cert.corrected_bound/np.linalg.norm(F),covered=True))
    out=ROOT/'runs/a2/vector_certificate';out.mkdir(parents=True,exist_ok=True)
    (out/'results.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='3-D vector radiatively corrected point-dipole finite systems; positive-loss materials only. No continuous Maxwell accuracy or inverse algorithm acceptance.',rows=rows),indent=2))
    print(json.dumps(dict(checks=len(rows),all_covered=all(r['covered'] for r in rows))))
if __name__=='__main__':main()
