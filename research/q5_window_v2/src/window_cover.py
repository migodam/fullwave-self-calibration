"""Finite full-world pair covering, with exact rational interval arithmetic.

An oracle must enclose the real observation vector over its entire input box.
Its proof includes the forward model, common nuisance parameters, units and
whitening. This module does NOT certify an arbitrary oracle, or estimate a
model-error bound from mesh differences. Class C has no such oracle yet.
"""
from fractions import Fraction as F
from dataclasses import dataclass
from typing import Callable, Sequence
from certify_modal import I

@dataclass(frozen=True)
class CertificateAssumptions:
    oracle_proof: str
    model_bound_proof: str
    noise_radius: F
    model_radius: F
    def validate(self):
        if not self.oracle_proof or not self.model_bound_proof:
            return False
        if self.noise_radius<0 or self.model_radius<0:
            raise ValueError('Error radii must be nonnegative')
        return True

def distance_lower_squared(a: Sequence[I], b: Sequence[I]) -> F:
    if len(a)!=len(b):raise ValueError('Observation dimensions differ')
    total=F(0)
    for u,v in zip(a,b):
        gap=max(F(0),u.lo-v.hi,v.lo-u.hi)
        total+=gap*gap
    return total

def finite_pair_cover(domain: Sequence[I], oracle: Callable,
                      material_indices: Sequence[int], tolerance: Sequence[F],
                      assumptions: CertificateAssumptions, max_nodes: int=20000):
    """Sufficient uniform L-infinity material certificate, not a local SVD.

    Both full worlds retain their own nuisance variables; parameters shared
    across acquisitions must be shared INSIDE each oracle call. On resource
    exhaustion all unprocessed boxes remain unresolved. Strict factor two is
    required for the coordinatewise-midrange estimator.
    """
    if not assumptions.validate():
        return {'status':'unresolved','reason':'missing proved oracle or model bound'}
    if max_nodes<1:raise ValueError('max_nodes must be positive')
    domain=tuple(domain);indices=tuple(material_indices)
    tol=tuple(map(F,tolerance))
    if len(indices)!=len(tol) or not indices or any(v<=0 for v in tol):
        raise ValueError('Invalid material tolerance')
    if any(i<0 or i>=len(domain) for i in indices):raise ValueError('Invalid index')
    scales=[v.hi-v.lo for v in domain]
    if any(v<=0 for v in scales):raise ValueError('Eliminate fixed parameters first')
    threshold=4*(assumptions.noise_radius+assumptions.model_radius)**2
    stack=[(domain,domain)];visited=separated=irrelevant=0
    while stack and visited<max_nodes:
        left,right=stack.pop();visited+=1
        possible_bad=any(max(abs(left[i].lo-right[i].hi),abs(left[i].hi-right[i].lo))>2*d
                         for i,d in zip(indices,tol))
        if not possible_bad:irrelevant+=1;continue
        if distance_lower_squared(oracle(left),oracle(right))>threshold:
            separated+=1;continue
        flat=left+right
        axis=max(range(len(flat)),key=lambda j:(flat[j].hi-flat[j].lo)/scales[j%len(domain)])
        cell=flat[axis];mid=(cell.lo+cell.hi)/2
        for half in (I(cell.lo,mid),I(mid,cell.hi)):
            new=list(flat);new[axis]=half
            stack.append((tuple(new[:len(domain)]),tuple(new[len(domain):])))
    return {'status':'certified' if not stack else 'unresolved',
            'visited':visited,'separated_boxes':separated,
            'irrelevant_boxes':irrelevant,'unprocessed_boxes':len(stack),
            'noise_plus_model_radius':str(assumptions.noise_radius+assumptions.model_radius),
            'scope':'conditional on the declared continuum enclosure and error proofs'}

def class_C_status():
    return {'status':'unresolved','reason':'coercivity exists, but a continuum field enclosure and actual residual upper bound are not yet supplied'}
