"""Registered finite-candidate acquisition development; no global certificate."""
from recovery_cycle1 import *


def main():
    output=HERE/'results/development_cycle2.json'
    if output.exists(): raise RuntimeError('Immutable output exists')
    old=json.loads((HERE/'results/development_cycle1.json').read_text())
    assert old['complete']
    rng=np.random.default_rng(2026091101); extra_rng=np.random.default_rng(2026091102)
    rx=mx.receivers(12,.6);valid=mx.receivers(16,.6);struct=mx.receivers(36,.15)
    grid=np.vstack([mx.receivers(24,r) for r in (.14,.2)])
    dda=mx.DipoleVIE(CENTERS,RADII,.011,K,fill_quadrature=4)
    rows=[];start=time.perf_counter()
    paths=[Path(__file__),HERE/'src/recovery_cycle1.py',HERE/'docs/CYCLE2_PROTOCOL.md']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    def save(complete=False):
        data={'complete':complete,'source_hashes':hashes,'rows':rows,'wall_seconds':time.perf_counter()-start,
              'scope':'Reused development scenes; pair-based upper-bound proxy, no certified acceptance.'}
        temp=output.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(output)
    for scene in range(24):
        truth=np.r_[rng.uniform([1.5,2],[4,5]),rng.uniform(-2,2)]
        gain=rng.uniform(.75,1.25)*np.exp(1j*rng.uniform(-np.pi,np.pi))
        if scene<12:
            yy=forward(truth,rx,4)
        else:
            p=dda.currents(truth[:2]+1j*LOSS)
            def df(points,shift=True):
                pts=points+np.array([truth[2]*.001,0,0]) if shift else points
                return (mx.dipole_kernel(pts,dda.points,K)@p).ravel()
            yy=df(rx)
        sigma=.01*np.linalg.norm(gain*yy)/np.sqrt(yy.size)
        y=gain*yy+sigma/np.sqrt(2)*(rng.normal(size=yy.size)+1j*rng.normal(size=yy.size))
        reference=gain+.01/np.sqrt(2)*(rng.normal()+1j*rng.normal())
        assert np.allclose(truth,old['rows'][scene]['truth'],rtol=0,atol=1e-14)
        if scene<12: continue
        before=time.perf_counter()
        t0=np.array(old['rows'][scene]['methods']['no_reference']['theta'])
        _,g0=profile(y,forward(t0,rx),sigma)
        alternatives=[]
        for shift in (-.25,.25):
            e1=t0[0]+shift
            if not 1.5<=e1<=4:continue
            fit=least_squares(lambda u:profile(y,forward([e1,*u],rx),sigma)[0],t0[1:],
                bounds=([2,-2],[5,2]),max_nfev=50,ftol=1e-8,xtol=1e-8,gtol=1e-8)
            ta=np.r_[e1,fit.x];_,ga=profile(y,forward(ta,rx),sigma)
            alternatives.append((ta,ga))
        assert alternatives
        predictions=g0*forward(t0,grid)
        score=np.min([np.abs(predictions-ga*forward(ta,grid))/sigma for ta,ga in alternatives],axis=0)
        refscore=min(abs(g0-ga)/.01 for ta,ga in alternatives)
        selected=int(np.argmax(score));random=int(extra_rng.integers(len(score)))
        decision='reference' if refscore>=score[selected] else 'selected_EM'
        design_seconds=time.perf_counter()-before
        rawgrid=gain*df(grid)+sigma/np.sqrt(2)*(extra_rng.normal(size=len(score))+1j*extra_rng.normal(size=len(score)))
        vv=df(valid);ss=df(struct,False)
        row={'scene':scene,'truth':truth.tolist(),'selected_index':selected,'random_index':random,
             'selected_nominal_receiver':grid[selected//12].tolist(),
             'selected_component':(selected%12)//4,'selected_illumination':selected%4,
             'em_score':float(score[selected]),'reference_score':float(refscore),'policy_action':decision,
             'design_seconds':design_seconds,'methods':{},'alternatives':[
                 {'theta':t.tolist(),'gain':[g.real,g.imag]} for t,g in alternatives],
             'base_data':[[v.real,v.imag] for v in y], 'reference_data':[reference.real,reference.imag],
             'sigma':sigma,'gain_truth':[gain.real,gain.imag]}
        for name,index in [('selected_EM',selected),('random_EM',random),('fixed_EM',0)]:
            before=time.perf_counter();point=grid[index//12:index//12+1];component=index%12
            def fields(t):return np.r_[forward(t,rx),forward(t,point)[component]]
            data=np.r_[y,rawgrid[index]];fits=[]
            for initial in ([2,3,0],[1.6,2.2,-1],[3.8,4.8,1]):
                fit=least_squares(lambda t:profile(data,fields(t),sigma)[0],initial,
                    bounds=([1.5,2,-2],[4,5,2]),max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8)
                fits.append({'theta':fit.x.tolist(),'rss':float(fit.fun@fit.fun),'nfev':fit.nfev,'optimizer_success':bool(fit.success)})
            best=min(fits,key=lambda f:f['rss']);t=np.array(best['theta']);r,g=profile(data,fields(t),sigma)
            errors=np.abs(t-truth);ge=abs(g-gain)/abs(gain)
            row['methods'][name]={'starts':fits,'theta':t.tolist(),'gain':[g.real,g.imag],
                'material_errors':errors[:2].tolist(),'geometry_error_mm':float(errors[2]),'gain_relative_error':float(ge),
                'task_success':bool(np.all(errors[:2]<=.1) and errors[2]<=.5 and ge<=.05),
                'sensor_prediction_error':float(np.linalg.norm(g*forward(t,valid)-gain*vv)/np.linalg.norm(gain*vv)),
                'structural_field_error':float(np.linalg.norm(forward([*t[:2],0],struct)-ss)/np.linalg.norm(ss)),
                'whitened_rss':float(2*r@r),'heuristic_residual_pass':bool(2*r@r<=chi2.ppf(.99,len(r)-5)),
                'scientific_acceptance':'unresolved','fit_seconds':time.perf_counter()-before,
                'extra_datum':[rawgrid[index].real,rawgrid[index].imag]}
        rows.append(row);save()
        print(json.dumps({'scene':scene,'action':decision,'success':{n:m['task_success'] for n,m in row['methods'].items()}}),flush=True)
    save(True)


if __name__=='__main__':main()
