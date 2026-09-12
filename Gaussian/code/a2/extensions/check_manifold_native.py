"""Portable native-chart checks: derivatives, acceptance and counted RHS."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from a2.extensions.manifold_revised import initial,material,retract,make_solvers,fit,CFG
s=initial();v=make_solvers(16);c,J=material(v[0].points,s,True);errors=[]
for k in range(12):
 h=np.zeros(12);h[k]=1e-5;fd=(material(v[0].points,retract(s,h))-material(v[0].points,retract(s,-h)))/2e-5;errors.append(float(np.linalg.norm(fd-J[:,k])/np.linalg.norm(J[:,k])))
assert max(errors)<1e-7
field_errors=[]
for solver in v:
 c,dc=material(solver.points,s,True);fw=solver.forward(c,rtol=1e-11)
 for k in range(12):
  h=np.zeros(12);h[k]=1e-4
  fd=(solver.forward(material(solver.points,retract(s,h)),rtol=1e-11)['scattered']-solver.forward(material(solver.points,retract(s,-h)),rtol=1e-11)['scattered'])/2e-4
  analytic=solver.material_tangent(c,dc[:,k],forward=fw,rtol=1e-11)['scattered'];field_errors.append(float(np.linalg.norm(fd-analytic)/np.linalg.norm(analytic)))
assert max(field_errors)<1e-5
m0=retract(s,np.r_[np.array([.15,.1,-.1,.2,-.1,.1]),np.zeros(6)])
f=np.concatenate([q.forward(material(q.points,m0))['scattered'].ravel() for q in v]);mask=np.tile(np.repeat(np.arange(24)%2==0,6),2);CFG['max_outer']=2;out={}
for method in CFG['methods']:
 st,y,meta=fit(v,f,mask,method,42);loss=[row['loss'] for row in meta['history']];assert all(b<=a+1e-14 for a,b in zip(loss,loss[1:]));assert meta['rhs']==12*meta['forward_calls']+144*meta['tangent_calls'];out[method]=dict(rhs=meta['rhs'],loss=loss)
result=dict(field_fd_max=max(field_errors),fd_max=max(errors),methods=out,scope='N16 two-iteration mechanical check, not imaging evidence',source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('manifold_revised.py'),Path(__file__).with_name('manifold.py'),Path(__file__).parents[1]/'physics.py']});p=Path(__file__).resolve().parents[3]/'runs/a2/extensions/manifold_native_checks.json';p.write_text(json.dumps(result,indent=2));print(json.dumps(result))
