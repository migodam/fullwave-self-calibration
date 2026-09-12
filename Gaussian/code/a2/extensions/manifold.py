"""SPD Gaussian peak-amplitude chart utilities for A2 manifold experiments."""
import numpy as np
import scipy.linalg as la

SYM_BASIS=(np.array([[1.,0.],[0.,0.]]),np.array([[0.,0.],[0.,1.]]),np.array([[0.,1.],[1.,0.]])/np.sqrt(2))

def covariance(logsx,logsy,angle):
 c,s=np.cos(angle),np.sin(angle);r=np.array([[c,-s],[s,c]]);return r@np.diag(np.exp([2*logsx,2*logsy]))@r.T
def retract_spd(S,H):
 w,u=la.eigh(S);root=(u*np.sqrt(w))@u.T;inv=(u/np.sqrt(w))@u.T;return root@la.expm(inv@H@inv)@root
def peak_gaussian(points,a,mu,S):
 q=points-np.asarray(mu);inv=np.linalg.inv(S);return a*np.exp(-.5*np.einsum('ni,ij,nj->n',q,inv,q))
def covariance_direction(points,a,mu,S,H):
 q=points-np.asarray(mu);inv=np.linalg.inv(S);g=peak_gaussian(points,a,mu,S);return .5*g*np.einsum('ni,ij,jk,kl,nl->n',q,inv,H,inv,q)
def fd_covariance_direction(points,a,mu,S,H,eps=1e-5):
 return (peak_gaussian(points,a,mu,retract_spd(S,eps*H))-peak_gaussian(points,a,mu,retract_spd(S,-eps*H)))/(2*eps)
def chart_equivalence(points,a,cx,cy,lsx,lsy,angle):
 S=covariance(lsx,lsy,angle);old=a*np.exp(-.5*(((np.cos(angle)*(points[:,0]-cx)+np.sin(angle)*(points[:,1]-cy))/np.exp(lsx))**2+((-np.sin(angle)*(points[:,0]-cx)+np.cos(angle)*(points[:,1]-cy))/np.exp(lsy))**2));new=peak_gaussian(points,a,(cx,cy),S);return float(np.linalg.norm(old-new)/max(np.linalg.norm(old),1e-30))
