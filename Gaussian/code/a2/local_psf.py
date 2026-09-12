"""Algorithm-specific local impulse response on the Gaussian material tangent.

The image impulse is first projected by G-dagger. This is not a global PSF or
an arbitrary-shape resolution guarantee. Input/metric convention is serialized.
"""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.som import FREQUENCIES_HZ,chi_and_partials,tangent_operators,folds,solve_stack

def main():
    geometry=Geometry(n=32,n_tx=4,n_rx=24,aperture='half')
    vies=[VIE(geometry,f) for f in FREQUENCIES_HZ]
    theta=np.array([.58,-.042,-.015,np.log(.031),np.log(.024),.18,.46,.046,.018,np.log(.026),np.log(.036),-.31])
    _,forward,_=solve_stack(theta,vies,rtol=2e-9)
    field=np.concatenate([f['scattered'].ravel() for f in forward]);sigma=.01*np.linalg.norm(field)/np.sqrt(field.size)
    train=np.arange(24)%2==0;A,B,_,_,_=tangent_operators(theta,vies,train,sigma,rtol=2e-9)
    v1,v2,meta=folds(A,B,first_target=6,second_target=2)
    material,partials=chi_and_partials(theta,vies[0],FREQUENCIES_HZ[0]);G=partials.real.T
    lam=.01*np.trace(A.T@A)/theta.size
    center=np.array([-.042,-.015]);index=int(np.argmin(np.sum((vies[0].points-center)**2,axis=1)))
    impulse=np.zeros(32**2);impulse[index]=1.;q=np.linalg.pinv(G,rcond=1e-10)@impulse
    outputs={}
    for name,V in [('full_LM',np.eye(theta.size)),('first_fold_6',v1),('two_fold_8',np.column_stack((v1,v2)))]:
        R=V@np.linalg.solve(V.T@A.T@A@V+lam*np.eye(V.shape[1]),V.T@A.T@A)
        outputs[name]=G@R@q
    fig,axes=plt.subplots(1,4,figsize=(12,3.5),layout='constrained')
    axes[0].imshow(material.real.reshape(32,32).T,origin='lower',extent=(-20,20,-20,20),cmap='viridis')
    axes[0].plot(center[0]*100,center[1]*100,'w+',ms=9);axes[0].set_title('Linearization material\nwhite +: projected impulse')
    limit=max(np.max(np.abs(x)) for x in outputs.values())
    for ax,(name,y) in zip(axes[1:],outputs.items()):
        im=ax.imshow(y.reshape(32,32).T,origin='lower',extent=(-20,20,-20,20),cmap='RdBu_r',vmin=-limit,vmax=limit)
        ax.set_title(name.replace('_',' '))
    for ax in axes:ax.set(xlabel='x (cm)',ylabel='y (cm)')
    fig.colorbar(im,ax=list(axes[1:]),label='signed local response, same input and scale',shrink=.8)
    fig.suptitle('Finite Gaussian-tangent impulse response; not a global resolution claim',fontsize=11)
    out=ROOT/'runs/a2/psf';out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'arrays.npz',A=A,B=B,G=G,material=material,**outputs)
    (out/'results.json').write_text(json.dumps(dict(scope='Algorithm-specific local response GRGdagger; fixed known material and projected pixel impulse, no global recovery guarantee',noise_convention='complex RMS sigma; real Fisher is 2 A.T A',n=32,side_m=.4,frequencies=list(FREQUENCIES_HZ),sigma=float(sigma),damping=float(lam),input_pixel=index,parameter_rank=int(np.linalg.matrix_rank(G)),fold=meta),indent=2))
    fig.savefig(ROOT/'figures/a2/local_psf.png',dpi=170);print('Saved local Gaussian-tangent PSF diagnostic')
if __name__=='__main__':main()
