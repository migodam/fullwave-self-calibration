"""Local parent tests: dimensionality, pre-call budgets and chart guard."""
import json
import numpy as np
from controller import A2, required_rhs, bounds, pack, unpack, solve, chart_guard, OUT

def run():
    rows=[]
    for q in [9,49]:
        pts,_=A2.make_grid(8)
        axis=np.linspace(-.6,.6,int(np.sqrt(q)))
        centers=np.array([(i,j) for i in axis for j in axis])
        mb=np.exp(-np.sum((pts[:,None,:]-centers[None,:,:])**2,axis=2)/.12)
        m=A2.Model(A2.Config(N=8,material_basis=mb))
        a=np.full(q,.1); x=np.array([.02,-.01,.03])
        a2,x2=unpack(pack(a,x))
        assert len(bounds(q))==q+3 and np.allclose(a,a2) and np.allclose(x,x2)
        for kind in ['forward','jacobian','adjoint']:
            if kind=='adjoint':
                yy=m.forward(a,x,[0,1,3])['total']
                fw=m.adjoint_gradient(a,x,yy,1.,[0,1,3])
            else:
                fw=m.forward(a,x,[0,1,3],jacobian=kind=='jacobian')
            assert fw['work']['rhs_solves_total']==required_rhs(m,kind)
        full=m.forward(a,np.zeros(3))
        data=dict(y=full['total'],sigma=.01,x0=x,alpha_true=a,x_true=np.zeros(3),
                  true_scattered=m.forward(a,np.zeros(3),[0,1,3])['scattered'])
        for method in ['direct_gn','direct_adjoint']:
            out=solve(m,data,method,maxiter=2,budget=17)
            assert out['online_work']['full_rhs_columns']==0
            assert out['status']=='budget_stopped_before_call'
        rows.append(dict(q=q,dimension_and_pre_call_budget_pass=True))
    fw=dict(total=np.array([100.+0j]),scattered=np.array([.1+0j]))
    red=dict(total=np.array([100.01+0j]),max_state_res_rel=1e-5,per_freq={3:dict(rank_ok=True)})
    assert not chart_guard(red,fw)['passed']
    rows.append(dict(dominant_direct_field_cannot_mask_scattering_error=True))
    (OUT/'parent_checks.json').write_text(json.dumps(dict(passed=True,checks=rows),indent=2)+'\n')
    print(json.dumps(rows))

if __name__=='__main__': run()
