import json,sys,numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.explicit_current_v2 import split
v=VIE(Geometry(n=24,n_tx=2,n_rx=16),1.5e9);rng=np.random.default_rng(19);y=rng.normal(size=(8,2))+1j*rng.normal(size=(8,2));_,P,_=split(v,np.arange(16)%2==0,y);x=rng.normal(size=v.D.npix)+1j*rng.normal(size=v.D.npix);z=rng.normal(size=v.D.npix)+1j*rng.normal(size=v.D.npix);lhs=np.vdot(v.D.matvec(P(x)),z);rhs=np.vdot(x,P(v.D.rmatvec(z)));err=float(abs(lhs-rhs)/max(abs(lhs),1e-30));assert err<1e-10
out=ROOT/'runs/a2/explicit_current_v2';(out/'composition_check.json').write_text(json.dumps({'PDP_adjoint_relative':err,'scope':'discrete composition check only'}));print(err)
