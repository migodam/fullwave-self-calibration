"""Budgeted set-membership exclusion with explicitly retained unresolved cells.

Analytic bounds are exact-arithmetic bounds. NumPy evaluation has a small
conservative padding, but is NOT directed-rounding interval certification.
The reported operational status is conditional/numerically audited coverage.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
from time import perf_counter
import heapq
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
from modal import unit_pattern, profile_complex, tensor_proposal, dyad_jacobian, realify, profile_columns

@dataclass
class Block:
    y: np.ndarray
    k: float
    offset: np.ndarray

@dataclass
class Cell:
    center: np.ndarray
    half: np.ndarray
    def contains(self,r): return bool(np.all(abs(np.asarray(r)-self.center)<=self.half+1e-12))
    def farthest(self,r): return float(np.linalg.norm(abs(self.center-np.asarray(r))+self.half))
    def split(self):
        i=int(np.argmax(self.half)); h=self.half.copy(); h[i]*=.5
        shift=np.zeros(3); shift[i]=h[i]
        return Cell(self.center-shift,h),Cell(self.center+shift,h)

class ModalProblem:
    def __init__(self,blocks):
        self.blocks=blocks; self.evaluations=0
    def residual(self,r):
        self.evaluations+=1; res=[]
        for b in self.blocks:
            if np.linalg.norm(np.asarray(r)+b.offset)<1e-10:
                return np.ones(18*len(self.blocks))*1e8
            f=unit_pattern(np.asarray(r)+b.offset,b.k)
            _,v=profile_complex(b.y,f); res.extend([v.real,v.imag])
        return np.concatenate(res)
    def norm(self,r): return float(np.linalg.norm(self.residual(r)))
    def threshold(self,alpha=.01,beta=0.):
        return float(np.sqrt(chi2.ppf(1-alpha,18*len(self.blocks))/2)+beta)
    def lower_bound(self,cell):
        if any(np.linalg.norm(cell.center+b.offset)<1e-10 for b in self.blocks): return 0.
        rho=np.linalg.norm(cell.half); changes=[]
        for b in self.blocks:
            rmin=np.linalg.norm(np.maximum(abs(cell.center+b.offset)-cell.half,0.))
            dp=1. if rmin<=0 else min(1.,np.sqrt(3)*rho/rmin)
            changes.append(np.linalg.norm(b.y)*dp)
        # Bound projector variation, then use the reverse triangle inequality.
        pad=1e-10*(1+sum(np.linalg.norm(b.y) for b in self.blocks))
        return max(0.,self.norm(cell.center)-np.linalg.norm(changes)-pad)
    def fit(self,roots,max_nfev=160,extra_starts=None):
        starts=[]
        b=self.blocks[0]
        for p in tensor_proposal(b.y,b.k): starts.append(p-b.offset)
        if extra_starts: starts.extend(extra_starts)
        answers=[]; before=self.evaluations; t=perf_counter()
        for root in roots:
            choices=[root.center]
            choices.extend(p for p in starts if root.contains(p))
            best=None
            for p in choices:
                lo=root.center-root.half; hi=root.center+root.half
                p=np.maximum(lo+1e-10,np.minimum(hi-1e-10,p))
                fit=least_squares(self.residual,p,bounds=(lo,hi),max_nfev=max_nfev,
                                  ftol=1e-10,xtol=1e-10,gtol=1e-10)
                row={'r':fit.x,'norm':float(np.linalg.norm(fit.fun)),
                     'success':bool(fit.success),'nfev':fit.nfev}
                if best is None or row['norm']<best['norm']: best=row
            answers.append(best)
        # Keep a deterministic tie rule, not one that consults the true branch.
        answers.sort(key=lambda a:a['norm'])
        return answers,{'seconds':perf_counter()-t,'profile_evaluations':self.evaluations-before}
    def local_information(self,r):
        I=np.zeros((3,3))
        for b in self.blocks:
            from modal import dyad
            m=dyad(np.asarray(r)+b.offset,b.k).ravel()
            c,_=profile_complex(b.y,m)
            B=realify(c*dyad_jacobian(np.asarray(r)+b.offset,b.k))
            N=realify(np.column_stack([m,1j*m])); V,_=profile_columns(B,N)
            I+=V.T@V
        return I


def cover(problem,roots,estimate,tolerance=.015,alpha=.005,beta=0.,max_cells=12000,min_half=1e-5):
    """Never discard budget-exhausted cells. A local optimizer is not a bound."""
    t=perf_counter(); before=problem.evaluations; eps=problem.threshold(alpha,beta)
    queue=[]; serial=0; excluded=0; inside=[]; stalled=[]; examined=0
    for cell in roots:
        heapq.heappush(queue,(-cell.farthest(estimate),serial,cell)); serial+=1
    while queue and examined<max_cells:
        _,_,cell=heapq.heappop(queue); examined+=1
        if problem.lower_bound(cell)>eps:
            excluded+=1; continue
        if cell.farthest(estimate)<=tolerance:
            inside.append(cell); continue
        if np.max(cell.half)<=min_half:
            stalled.append(cell); continue
        for child in cell.split():
            heapq.heappush(queue,(-child.farthest(estimate),serial,child)); serial+=1
    pending=stalled+[x[2] for x in queue]
    feasible_estimate=problem.norm(estimate)<=eps
    accepted=bool(feasible_estimate and not pending and inside)
    cells=inside+pending
    return {'accepted':accepted,'status':'accept_conditional' if accepted else ('reject_model_or_prior' if not cells else 'unresolved'),
            'threshold':eps,'examined':examined,'excluded':excluded,'inside_cells':len(inside),
            'unresolved_cells':len(pending),'seconds':perf_counter()-t,
            'profile_evaluations':problem.evaluations-before,
            'enclosure_kind':'analytic_bound_float_checked_not_interval',
            'cells':cells}


def choose_offset(candidates,k,pool,method,rng,q_design=100.):
    """All scores use frozen training candidates, never the true receiver."""
    if method=='random': return np.asarray(pool[int(rng.integers(len(pool)))]),[]
    scores=[]
    for d in pool:
        if method=='branch':
            angles=[]
            for i in range(len(candidates)):
                for j in range(i):
                    f=unit_pattern(candidates[i]+d,k); g=unit_pattern(candidates[j]+d,k)
                    angles.append(max(0.,1-abs(np.vdot(f,g))**2))
            score=q_design**2*min(angles) if angles else 0.
        elif method=='fisher':
            from modal import theoretical_singular_values
            score=min(theoretical_singular_values(np.linalg.norm(candidates[0]+d),k,q_design))**2
        else: raise ValueError(method)
        scores.append(float(score))
    return np.asarray(pool[int(np.argmax(scores))]),scores
