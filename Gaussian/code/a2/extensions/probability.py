"""Finite latent-material full-wave Bayes utilities; no mean-material shortcut."""
import numpy as np
from scipy.optimize import minimize

def all_states(k=8): return ((np.arange(1<<k)[:,None]>>np.arange(k))&1).astype(int)
def log_prior(states,eta=.65,field=-.12):
    # Ring adjacency of fixed spatial patches; fully normalized by enumeration.
    return eta*np.sum(states*np.roll(states,-1,axis=1)+(1-states)*(1-np.roll(states,-1,axis=1)),axis=1)+field*states.sum(1)
def normalized_logweights(x):
    x=x-x.max();return x-np.log(np.exp(x).sum())
def q_probs(q,states): return np.prod(np.where(states,q,1-q),axis=1)
def elbo_product(loglike,lp,states,logits):
    q=1/(1+np.exp(-np.clip(logits,-60,60)));w=q_probs(q,states); entropy=-np.sum(q*np.log(np.clip(q,1e-15,1))+(1-q)*np.log(np.clip(1-q,1e-15,1)))
    return float(w@(loglike+lp)+entropy),q
def fit_product(loglike,lp,states,x0=None):
    x0=np.zeros(states.shape[1]) if x0 is None else x0
    o=minimize(lambda x:-elbo_product(loglike,lp,states,x)[0],x0,method='L-BFGS-B',options={'maxiter':200,'ftol':1e-12})
    e,q=elbo_product(loglike,lp,states,o.x);return o,e,q
def fit_gaussian_logit(loglike,lp,states,P):
    o=minimize(lambda t:-elbo_product(loglike,lp,states,P@t)[0],np.zeros(P.shape[1]),method='L-BFGS-B',options={'maxiter':200,'ftol':1e-12})
    e,q=elbo_product(loglike,lp,states,P@o.x);return o,e,q
