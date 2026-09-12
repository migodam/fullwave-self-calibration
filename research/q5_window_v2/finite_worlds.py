"""Search upper bounds on finite class C distances. NOT continuum certificates."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.special import ndtr,ndtri
from sphere_cluster import SphereCluster
from recovery_v2 import receivers,profiled_gain,realwhite,complex_record
HERE=Path(__file__).resolve().parent

def run(out):
    if out.exists():raise FileExistsError(out)
    model=SphereCluster(order=3,ntheta=12);rx=receivers(12,.6);rows=[]
    gA=.8+0j
    for e1 in [2.,3.,4.]:
        for e2 in [3.,4.,5.]:
            a=np.array([e1,e2,0.]);fa=model.forward(a,rx);muA=gA*fa
            sigma=.01*np.linalg.norm(muA)/np.sqrt(len(muA));e2b=e2-.21
            def objective(v,ref=False):
                b=np.array([v[0],e2b,v[1]]);fb=model.forward(b,rx)
                gb=profiled_gain(fb,muA,sigma,z=gA if ref else None)
                rr=(muA-gb*fb)/sigma
                if ref:rr=np.r_[rr,(gA-gb)/.01]
                return realwhite(rr)
            for ref in [False,True]:
                starts=[[max(1.5,e1-.15),0.],[e1,0.]]
                fits=[least_squares(lambda v:objective(v,ref),s,bounds=([1.5,-2],[4,2]),max_nfev=100,
                       ftol=1e-11,xtol=1e-11,gtol=1e-11) for s in starts]
                sol=min(fits,key=lambda x:np.linalg.norm(x.fun))
                b=np.array([sol.x[0],e2b,sol.x[1]]);fb=model.forward(b,rx)
                gb=profiled_gain(fb,muA,sigma,z=gA if ref else None)
                distance=np.linalg.norm(sol.fun)/np.sqrt(2)
                rows.append({'A':a.tolist(),'B':b.tolist(),'gain_A':complex_record(gA),'gain_B':complex_record(gb),
                    'sigma_sensor_common':sigma,'sigma_reference_common':.01 if ref else None,
                    'numerical_complex_whitened_distance_upper_candidate':float(distance),
                    'model_based_equal_prior_binary_error':float(ndtr(-distance/np.sqrt(2))),
                    'reference':ref})
    selected=[]
    high=SphereCluster(order=5,ntheta=18)
    for ref in [False,True]:
        row=min([r for r in rows if r['reference']==ref],key=lambda r:r['numerical_complex_whitened_distance_upper_candidate'])
        g=row['gain_B'];gb=g['real']+1j*g['imag'];fa=high.forward(row['A'],rx);fb=high.forward(row['B'],rx)
        rr=(gA*fa-gb*fb)/row['sigma_sensor_common']
        if ref:rr=np.r_[rr,(gA-gb)/.01]
        dh=np.linalg.norm(rr);target=np.sqrt(2)*ndtri(.75)
        selected.append({**row,'same_worlds_lmax5_distance':float(dh),
          'high_order_difference_NOT_error_bound':float(abs(dh-row['numerical_complex_whitened_distance_upper_candidate'])),
          'per_world_whitened_mean_error_cap_needed_for_p_at_least_quarter':float((target-dh)/2)})
    result={'status':'numerical_upper_candidates_only_continuum_verification_unresolved','candidate_count':len(rows),
        'all_candidates':rows,'selected':selected,
        'not_a_claim':['No optimization minimum is a separation lower bound','No truncation difference is a model-error bound',
                       'No worst-case statistical risk claim for continuum Maxwell without mean error certificate']}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':run(HERE/'results/finite_worlds_diagnostic.json')
