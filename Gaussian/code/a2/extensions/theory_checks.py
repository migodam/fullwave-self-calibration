"""Finite algebra checks, not a substitute for full-wave imaging validation."""
from pathlib import Path
import json,hashlib
import numpy as np
rng=np.random.default_rng(2)
x=rng.normal(size=(40,2));Sigma=np.array([[.3,.05],[.05,.2]]);H=np.array([[.03,.01],[.01,-.01]]);iv=np.linalg.inv(Sigma)
def g(T):return np.exp(-.5*np.einsum('ni,ij,nj->n',x,np.linalg.inv(T),x))
an=.5*g(Sigma)*np.einsum('ni,ij,jk,kl,nl->n',x,iv,H,iv,x);fd=(g(Sigma+1e-5*H)-g(Sigma-1e-5*H))/2e-5
norm=lambda T:g(T)/np.sqrt(np.linalg.det(T));an_norm=(an-.5*g(Sigma)*np.trace(iv@H))/np.sqrt(np.linalg.det(Sigma));fd_norm=(norm(Sigma+1e-5*H)-norm(Sigma-1e-5*H))/2e-5
# Fourier response of unit-mass symmetric children: exp(+t^2/2) cos(t).
q=np.array([.7,.3]);split=[]
for delta in [.12,.06,.03]:
 d=np.array([delta,0.]);Sp=Sigma-np.outer(d,d);ac=.5*np.sqrt(np.linalg.det(Sigma)/np.linalg.det(Sp));mass=2*ac*np.sqrt(np.linalg.det(Sp))/np.sqrt(np.linalg.det(Sigma));split.append(dict(delta=delta,mass_ratio=mass,covariance_error=float(np.linalg.norm(Sp+np.outer(d,d)-Sigma)),fourier_relative=float(abs(np.exp((q@d)**2/2)*np.cos(q@d)-1))))
states=np.array([[-1,-1],[-1,1],[1,-1],[1,1.]])
p=np.ones(4)/4;pc=np.array([.5,0,0,.5]);A=np.array([[1,2j],[.2+.1j,1]])
def moments(fields,weights):
 mean=weights@fields;z=fields-mean;return mean,np.einsum('i,ij,ik->jk',weights,z,z.conj()),np.einsum('i,ij,ik->jk',weights,z,z)
b1=moments(states@A.T,p);b2=moments(states@A.T,pc)
D=np.array([[.1+.15j,.02+.03j],[.02+.03j,.08+.12j]]);M=np.diag([.2+.04j,.3+.06j]);S=np.array([[1,.2j],[.3,1-.1j]]);E=np.array([1,.8+.3j]);B=np.linalg.solve(np.eye(2)-D@M,D);U=np.linalg.inv(np.eye(2)-M@D);E0=np.linalg.solve(np.eye(2)-D@M,E)
checks=[]
for delta in [.2,.1,.05]:
 dc=delta*(1+.2j)*states;fields=[]
 for z in dc:
  X=M+np.diag(z);fields.append(S@X@np.linalg.solve(np.eye(2)-D@X,E))
 fields=np.asarray(fields);C=np.einsum('i,ij,ik->jk',p,dc,dc);pred=S@M@E0+S@U@(C*B)@E0
 qmax=max(np.linalg.norm(B@np.diag(z),2) for z in dc);dmax=max(np.linalg.norm(np.diag(z),2) for z in dc);bound=np.linalg.norm(S@U,2)*dmax*qmax*qmax/(1-qmax)*np.linalg.norm(E0)
 error=np.linalg.norm(p@fields-pred);checks.append(dict(delta=delta,error=float(error),bound=float(bound),q=float(qmax),different_joint_mean=float(np.linalg.norm(p@fields-pc@fields))));assert error<=bound
conditional=np.array([[.9,.1],[.1,.9]]);label_prior=np.array([.5,.5]);static=float(label_prior@conditional.prod(axis=1));resampled=float((label_prior@conditional).prod());assert abs(static-resampled)>.1
out=dict(scope='Finite complex matrices and Fourier moments only; no Maxwell validation claim',peak_fd_relative=float(np.linalg.norm(an-fd)/np.linalg.norm(an)),normalized_fd_relative=float(np.linalg.norm(an_norm-fd_norm)/np.linalg.norm(an_norm)),incorrect_peak_trace_relative=float(np.linalg.norm(an-.5*g(Sigma)*np.trace(iv@H)-fd)/np.linalg.norm(fd)),split=split,born_mean_difference=float(np.linalg.norm(b1[0]-b2[0])),born_covariance_difference=float(np.linalg.norm(b1[1]-b2[1])),born_pseudocovariance_norm=float(np.linalg.norm(b1[2])),second_order=checks,static_likelihood=static,resampled_likelihood=resampled,misaligned_product_norm=float(np.linalg.norm(np.diag([1,0])@np.diag([0,1]))))
assert out['peak_fd_relative']<1e-7 and out['normalized_fd_relative']<1e-7
assert out['incorrect_peak_trace_relative']>.1
assert out['born_mean_difference']<1e-14 and out['born_covariance_difference']>0
assert checks[-1]['error']<checks[0]['error']/100
assert 3.8<np.log2(split[0]['fourier_relative']/split[1]['fourier_relative'])<4.2
out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
target=Path(__file__).resolve().parents[3]/'runs/a2/extensions/theory_checks.json';target.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
