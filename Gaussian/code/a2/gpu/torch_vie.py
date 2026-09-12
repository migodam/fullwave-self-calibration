"""Implicit full-wave PyTorch VIE with batched nonsymmetric BiCGSTAB.

Kernel is fixed. Parameters enter only chi. Every solve checks its true residual.
"""
import torch
from a2.physics import VIE
class TorchVIE:
 def __init__(self,vie:VIE,device='cpu',dtype=torch.complex128,tol=None):
  self.vie,self.device,self.dtype=vie,torch.device(device),dtype;self.n=vie.D.n;self.N=self.n**2;self.shape=vie.D.fft_shape;self.tol=tol or (2e-5 if dtype==torch.complex64 else 1e-9)
  k=torch.as_tensor(vie.D.kernel,dtype=dtype,device=device);self.fk=torch.fft.fft2(k,s=self.shape);self.fka=torch.fft.fft2(k.conj(),s=self.shape)
  self.E=torch.as_tensor(vie.E,dtype=dtype,device=device);self.S=torch.as_tensor(vie.S,dtype=dtype,device=device);self.last={};self.counts=dict(forward_rhs=0,adjoint_rhs=0,linear_iterations=0,solves=0)
 def D(self,x,adj=False):
  if x.ndim==1:x=x[:,None]
  z=torch.fft.ifft2(torch.fft.fft2(x.T.reshape(-1,self.n,self.n),s=self.shape)*(self.fka if adj else self.fk))
  return z[:,self.n-1:2*self.n-1,self.n-1:2*self.n-1].reshape(-1,self.N).T
 def A(self,chi,x,adj=False):return x-self.D(chi.conj()[:,None]*x,True) if adj else x-chi[:,None]*self.D(x)
 def solve_linear(self,chi,b,adj=False,maxiter=300):
  x=torch.zeros_like(b);r=b.clone();rh=r.clone();p=torch.zeros_like(b);v=p.clone();one=torch.ones(b.shape[1],dtype=b.dtype,device=b.device);rho_old=one.clone();alpha=one.clone();omega=one.clone();bn=torch.linalg.vector_norm(b,dim=0);active=bn>0
  dot=lambda a,b:(a.conj()*b).sum(0)
  safe=lambda z:torch.where(z.abs()>1e-30,z,one)
  it=0
  for it in range(1,maxiter+1):
   if not bool(active.any()):break
   rho=dot(rh,r);beta=(rho/safe(rho_old))*(alpha/safe(omega));p=r+beta*(p-omega*v);p[:,~active]=0;v=self.A(chi,p,adj);alpha=rho/safe(dot(rh,v));s=r-alpha*v;t=self.A(chi,s,adj);omega=dot(t,s)/safe(dot(t,t));x=x+alpha*p+omega*s;r=s-omega*t
   active=torch.linalg.vector_norm(r,dim=0)>self.tol*bn
   r[:,~active]=0;rho_old=rho
  true_r=b-self.A(chi,x,adj);rel=torch.linalg.vector_norm(true_r,dim=0)/torch.clamp(bn,min=1e-30)
  if not bool(torch.isfinite(rel).all()) or bool((rel>self.tol*1.5).any()):raise RuntimeError('BiCGSTAB true residual failed: '+str(rel.detach().cpu().tolist()))
  self.last=dict(iterations=it,relative_residual=rel.detach().cpu().tolist());self.counts['linear_iterations']+=it;self.counts['solves']+=1;self.counts['adjoint_rhs' if adj else 'forward_rhs']+=b.shape[1]
  return x
 def receive(self,j):return self.S@j
 def receive_adjoint(self,r):return self.S.conj().T@r
 def forward(self,chi):
  j=self.solve_linear(chi,chi[:,None]*self.E);return self.receive(j),j
 def scattered(self,chi):return _Implicit.apply(chi,self)
class _Implicit(torch.autograd.Function):
 @staticmethod
 def forward(ctx,chi,backend):
  f,j=backend.forward(chi);ctx.backend=backend;ctx.save_for_backward(chi,j);return f
 @staticmethod
 def backward(ctx,grad):
  chi,j=ctx.saved_tensors;b=ctx.backend;lam=b.solve_linear(chi,b.receive_adjoint(grad),adj=True);total=b.E+b.D(j);return (total.conj()*lam).sum(1),None

class TorchMultiVIE(TorchVIE):
 """Batch independent frequency state equations sharing the same chi."""
 def __init__(self,backends):
  b=backends[0];self.n=b.n;self.N=b.N;self.shape=b.shape;self.device=b.device;self.dtype=b.dtype;self.tol=b.tol;self.nf=len(backends);self.nt=b.E.shape[1]
  assert all(q.n==self.n and q.E.shape==b.E.shape and q.S.shape==b.S.shape for q in backends)
  self.E=torch.cat([q.E for q in backends],dim=1);self.S=torch.stack([q.S for q in backends]);self.fk=torch.stack([q.fk for q in backends])[:,None];self.fka=torch.stack([q.fka for q in backends])[:,None];self.last={};self.counts=dict(forward_rhs=0,adjoint_rhs=0,linear_iterations=0,solves=0)
 def D(self,x,adj=False):
  z=torch.fft.ifft2(torch.fft.fft2(x.T.reshape(self.nf,self.nt,self.n,self.n),s=self.shape)*(self.fka if adj else self.fk))
  return z[:,:,self.n-1:2*self.n-1,self.n-1:2*self.n-1].reshape(-1,self.N).T
 def receive(self,j):return self.S@j.T.reshape(self.nf,self.nt,self.N).transpose(1,2)
 def receive_adjoint(self,r):return (self.S.conj().transpose(1,2)@r).transpose(1,2).reshape(-1,self.N).T
