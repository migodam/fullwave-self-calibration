"""Post-registration mechanism diagnostics, never new fair baselines.

Oracle gain/geometry and same-model replacements expose possible mechanisms;
none of their success counts establishes a practical algorithm advantage.
"""
from run_experiment import *
from scipy.linalg import svd
from scipy.special import spherical_jn,spherical_yn
from multipole import modes,angular,quadrature


def decode(v):
    a=np.asarray(v);return a[...,0]+1j*a[...,1]


def fit_task(fun,truth,known_geometry=False):
    initial=[s[:2] if known_geometry else s for s in STARTS]
    bounds=([1.5,2],[4,5]) if known_geometry else BOUNDS
    fits=[least_squares(fun,s,bounds=bounds,max_nfev=80,ftol=1e-8,xtol=1e-8,gtol=1e-8) for s in initial]
    fit=min(fits,key=lambda r:r.fun@r.fun)
    return {'theta':fit.x.tolist(),'errors':np.abs(fit.x[:2]-truth[:2]).tolist(),
            'task_success':bool(np.all(np.abs(fit.x[:2]-truth[:2])<=.1)),
            'rss':float(fit.fun@fit.fun)}


def modal_checks():
    dirs,w=quadrature(3,6);Y,P,C=angular(dirs,1)
    dh=spherical_jn(1,2,True)+spherical_jn(1,2)/2+1j*(spherical_yn(1,2,True)+spherical_yn(1,2)/2)
    c=-np.sqrt(6*np.pi)
    row=-w[:,None]*P[:,:,1].conj()/(dh*c)
    ans=[]
    for x,interval in [(.2,[1.5,2.75,4]),(.15,[2,3.5,5])]:
        model=Cluster([[0,0,0]],[x],1.,[(np.array([1,0,0]),np.array([0,0,1]))],order=8)
        for eps in interval:
            field=model.field([eps],2*dirs)[:,:,0]
            measured=np.sum(row*field)
            target=model.coefficients((eps,))[len(modes(8))+1,0]/c
            ans.append({'size':x,'eps':eps,'absolute_modal_error':float(abs(measured-target))})
    return {'squared_row_norm':float(np.sum(abs(row)**2)),
            'expected_squared_row_norm':224/1053,'modal_tests':ans}


def main():
    path=HERE/'results/mechanism_diagnostics.json'
    if path.exists():raise FileExistsError(path)
    data=json.loads((HERE/'results/independent_v2.json').read_text())
    points=mx.receivers(12,.6); model=Cluster(CENTERS,RADII,K,mx.illuminations(),order=3)
    @lru_cache(maxsize=128)
    def obs(x):return model.observation_matrix(points+np.array([x*.001,0,0]))
    def forward(theta):return (obs(float(theta[2]))@model.coefficients(tuple(theta[:2]+1j*LOSS))).ravel()
    rows=[]
    for r in data['rows']:
        truth=np.array(r['truth']);g=complex(decode(r['gain']));sigma=r['sigma']
        y=decode(r['observed_all_fields'])[:144];ftrue=decode(r['noise_free_all_fields'])[:144]
        z=complex(decode(r['reference']));noise=y-g*ftrue
        fmatched=forward(truth);same=g*fmatched+noise
        out={}
        for tag,yy in [('same_model',same),('independent',y)]:
            for ref in (False,True):
                out[tag+('_reference' if ref else '_no_reference')]=fit_task(
                    lambda t:profile(yy,forward(np.array(t)),sigma,(z,.01) if ref else None)[0],truth)
        out['oracle_geometry_noisy_reference']=fit_task(
            lambda t:profile(y,forward(np.r_[t,truth[2]]),sigma,(z,.01))[0],truth,True)
        def fixed_gain(t,yy,geometry=False):
            theta=np.r_[t,truth[2]] if geometry else np.array(t)
            res=(yy-g*forward(theta))/sigma
            return np.sqrt(2)*np.r_[res.real,res.imag]
        out['oracle_gain']=fit_task(lambda t:fixed_gain(t,y),truth)
        out['oracle_gain_and_geometry']=fit_task(lambda t:fixed_gain(t,y,True),truth,True)
        out['noise_free_oracle_gain_and_geometry']=fit_task(lambda t:fixed_gain(t,g*ftrue,True),truth,True)
        def amplitude_only(t):
            f=forward(np.array(t));ga=abs(g)*np.exp(1j*np.angle(np.vdot(f,y)))
            res=(y-ga*f)/sigma
            return np.sqrt(2)*np.r_[res.real,res.imag]
        out['oracle_gain_amplitude_only']=fit_task(amplitude_only,truth)
        # Projection of deterministic model discrepancy onto the LOCAL real tangent.
        eps=1e-4;steps=np.array([eps,eps,eps])
        columns=[]
        for j in range(3):
            v=np.eye(3)[j]*steps[j]
            columns.append(g*(forward(truth+v)-forward(truth-v))/(2*steps[j]))
        columns += [fmatched,1j*fmatched]
        J=np.stack(columns,axis=1);JR=np.r_[J.real,J.imag]
        discrepancy=g*(ftrue-fmatched);br=np.r_[discrepancy.real,discrepancy.imag]
        projected=JR@np.linalg.lstsq(JR,br,rcond=None)[0]
        rows.append({'scene':r['scene'],'interventions':out,
            'local_tangent_model_error_fraction':float(np.linalg.norm(projected)/np.linalg.norm(br)),
            'model_relative_error_at_truth':float(np.linalg.norm(discrepancy)/np.linalg.norm(g*ftrue))})
    result={'scope':'posthoc mechanism/oracle diagnostics, not fair acquisition baselines',
            'rows':rows,'modal_readout_checks':modal_checks(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    result['counts']={key:sum(r['interventions'][key]['task_success'] for r in rows) for key in rows[0]['interventions']}
    path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['counts'],indent=2))
if __name__=='__main__':main()
