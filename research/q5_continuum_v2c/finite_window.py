"""Fail-closed set-membership task enclosure; not a novel acquisition selector.

A caller must supply a justified interval range for every complete state box,
including all shared nuisances, and justified componentwise model-error bounds.
The full class-C forward-range provider is NOT established by the two-world
certificate. Leaving proof_id=None therefore returns unresolved.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from typing import Callable,Sequence
from certified_modal import I,S

Box=tuple[tuple[F,F],...]

@dataclass(frozen=True)
class BoxModel:
    ranges: Callable[[Box],Sequence[I]]
    model_error: tuple[F,...]
    proof_id: str | None = None


def certify_task(model:BoxModel,observation:Sequence[F],noise_bounds:Sequence[F],
                 domain:Box,task_coordinates:Sequence[int],tolerance:F,
                 max_evaluations:int=10000):
    """Coordinatewise material error guarantee, conditional on the noise event.

    Accepted estimates are centers of a hull covering ALL non-excluded boxes,
    not merely optimizer neighborhoods. Exact rational subdivision has no gaps.
    The function deliberately does not treat finite optimized distances as a
    global lower bound. The proof_id records provenance, not an automatic proof.
    """
    if model.proof_id is None:
        return {'status':'unresolved','reason':'No justified complete-box/model-error provider'}
    y=tuple(map(F,observation));nb=tuple(map(F,noise_bounds));tol=F(tolerance)
    if not len(y)==len(nb)==len(model.model_error) or tol<0 or any(v<0 for v in nb+model.model_error):
        raise ValueError('Invalid dimensions or uncertainty bounds')
    root=tuple((F(a),F(b)) for a,b in domain)
    if any(a>b for a,b in root) or not task_coordinates:raise ValueError('Invalid domain/task')
    if any(i<0 or i>=len(root) for i in task_coordinates):raise ValueError('Invalid task coordinate')
    active=[];pending=[root];calls=0
    def compatible(box):
        vals=model.ranges(box)
        if len(vals)!=len(y) or not all(isinstance(v,I) for v in vals):
            raise ValueError('Provider must return proved integer interval ranges')
        for v,yy,n,b in zip(vals,y,nb,model.model_error):
            error=n+b
            if F(v.lo,S)>yy+error or F(v.hi,S)<yy-error:return False
        return True
    def result(status,boxes,reason=None):
        r={'status':status,'evaluations':calls,'remaining_boxes':len(boxes),'proof_id':model.proof_id}
        if boxes:
            hull=[(min(b[i][0] for b in boxes),max(b[i][1] for b in boxes)) for i in task_coordinates]
            r['task_hull']=[[str(a),str(b)] for a,b in hull]
            r['estimate']=[str((a+b)/2) for a,b in hull]
        if reason:r['reason']=reason
        return r
    while True:
        while pending:
            if calls>=max_evaluations:
                return result('unresolved',active+pending,'Box budget exhausted; untested boxes remain covered')
            box=pending.pop();calls+=1
            if compatible(box):active.append(box)
        if not active:return result('incompatible',[],'No state remains under the declared model/noise event')
        hull=[(min(b[i][0] for b in active),max(b[i][1] for b in active)) for i in task_coordinates]
        if all(b-a<=2*tol for a,b in hull):return result('conditional_task_certificate',active)
        # Split a widest box exactly. Ordering has no role in correctness.
        choices=[(max(b-a for a,b in box),j) for j,box in enumerate(active)]
        width,j=max(choices);box=active.pop(j)
        if width==0:return result('unresolved',active+[box],'Distinct compatible point states remain')
        axis=max(range(len(box)),key=lambda i:box[i][1]-box[i][0]);a,b=box[axis];mid=(a+b)/2
        left=list(box);right=list(box);left[axis]=(a,mid);right[axis]=(mid,b)
        pending.extend([tuple(left),tuple(right)])
