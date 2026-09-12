"""Independent local reconstruction of Pro A2 checks (original code unavailable)."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import sys,json,time,hashlib
import numpy as np,scipy.linalg as la
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'code'))
from a12_graph.core import Geometry,Operators,FREQ_GHZ,C0
from a2.certificate import gamma_diagonal,certify,galerkin,structure_decision
out=root/'runs/a2/certificate';out.mkdir(parents=True,exist_ok=True);fig=root/'figures/a2';fig.mkdir(exist_ok=True);t=time.perf_counter()
ops=Operators(Geometry(n=16,n_tx=4,n_rx=12),3);pts=ops.points;xx,yy=pts.T
materials=[.04+.6*np.exp(-((xx+.04)**2+(yy-.015)**2)/(.045**2)),.04+.5*np.exp(-((xx+.035)**2+yy**2)/(.04**2))+.4*np.exp(-((xx-.025)**2+yy**2)/(.032**2)),.04+.7*((abs(xx)<.055)&(abs(yy)<.038))]
bases={p:Operators(Geometry(n=16,n_tx=4,n_rx=12),p).Q for p in [1,3,6,10,15]};rows=[]
for fi,f in enumerate(FREQ_GHZ):
 D=ops.D[fi];E=ops.E[fi];C=ops.S[fi]
 for mi,a in enumerate(materials):
  chi=a*(1+.30j*2.25/f);L=np.diag(1/chi)-D;gam,reason=gamma_diagonal(chi,2*np.pi*f*1e9/C0,ops.h);full=la.solve(L,E);F=C@full;born=C@(chi[:,None]*E)
  for p,pv in zip([1,3,6,10],[3,6,10,15]):
   Q,V=bases[p],bases[pv];J=galerkin(L,E,Q);Z=galerkin(L,C.conj().T,V,True);cc=certify(L,E,C,J,Z,gam);actual=la.norm(F-cc.corrected_field);raw=la.norm(F-C@J);larger=C@galerkin(L,E,V)
   A=np.eye(len(chi))-chi[:,None]*D;JA=Q@la.solve(Q.conj().T@A@Q,Q.conj().T@(chi[:,None]*E))
   rows.append({'frequency_GHz':float(f),'material':mi,'R':Q.shape[1],'Rd':V.shape[1],'gamma_min':cc.gamma_min,'corrected_rel':actual/la.norm(F),'bound_rel':cc.corrected_bound/la.norm(F),'raw_rel':raw/la.norm(F),'raw_bound_rel':cc.raw_bound/la.norm(F),'covered':bool(actual<=cc.corrected_bound+1e-11*max(1,la.norm(F))),'nested_equivalence_rel':float(la.norm(cc.corrected_field-larger)/la.norm(F)),'A_vs_L_reduced_rel':float(la.norm(C@(JA-J))/la.norm(F)),'born_rel':float(la.norm(born-F)/la.norm(F))})
# Whole finite-dimensional Schur transfer with nonzero hidden input/output.
rng=np.random.default_rng(25);Z=rng.normal(size=(9,9))+1j*rng.normal(size=(9,9));L=(Z+Z.conj().T)/2-1j*(np.eye(9)+Z.conj().T@Z);B=rng.normal(size=(9,3))+1j*rng.normal(size=(9,3));C=rng.normal(size=(4,9))+1j*rng.normal(size=(4,9));aa,ab,ba,bb=L[:4,:4],L[:4,4:],L[4:,:4],L[4:,4:];Li=aa-ab@la.solve(bb,ba);Be=B[:4]-ab@la.solve(bb,B[4:]);Ce=C[:,:4]-C[:,4:]@la.solve(bb,ba);direct=C[:,4:]@la.solve(bb,B[4:]);full=C@la.solve(L,B);correct=Ce@la.solve(Li,Be)+direct;wrong=C[:,:4]@la.solve(Li,B[:4]);schur={'complete_relative_error':float(la.norm(correct-full)/la.norm(full)),'self_energy_only_relative_error':float(la.norm(wrong-full)/la.norm(full))}
# Changed support / nonlossy materials explicitly return not applicable.
_,zreason=gamma_diagonal(np.zeros(3),10,.02);ng,nreason=gamma_diagonal(np.ones(3),10,.02)
summary={'scope':'Independent reconstructed A2 discrete checks; NOT original Pro code rerun, no material inversion acceptance','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'schur':schur,'zero_contrast_reason':zreason,'lossless_applicable':bool(np.min(ng)>0),'decisions':[structure_decision(2,.5,1),structure_decision(.2,.5,1,old_witness=.7),structure_decision(2,3,1,new_lower=2)],'elapsed_seconds':time.perf_counter()-t}
(out/'results.json').write_text(json.dumps(summary,indent=2));assert all(x['covered'] for x in rows);assert max(x['nested_equivalence_rel'] for x in rows)<1e-10;assert schur['complete_relative_error']<1e-10
f,ax=plt.subplots(1,2,figsize=(10,4),constrained_layout=True);x=np.arange(len(rows));ax[0].semilogy(x,[z['bound_rel'] for z in rows],'.-',label='Certified corrected-field bound');ax[0].semilogy(x,[z['corrected_rel'] for z in rows],'.-',label='Actual corrected-field error');ax[0].set(xlabel='Fixed case index (36 checks)',ylabel='Relative scattered-field error',title='Coverage and conservatism are different');ax[0].legend(fontsize=8);ax[1].semilogy(x,[z['A_vs_L_reduced_rel'] for z in rows],'.');ax[1].set(xlabel='Fixed case index',ylabel='Relative field difference',title='A-Galerkin and L-Galerkin differ');f.savefig(fig/'certificate_checks.png',dpi=160);print(json.dumps({'checks':len(rows),'covered':sum(x['covered'] for x in rows),'max_nested_error':max(x['nested_equivalence_rel'] for x in rows),'schur':schur,'seconds':summary['elapsed_seconds']}))
