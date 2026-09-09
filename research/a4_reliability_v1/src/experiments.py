"""Self-contained A4 experiments. No A3 endpoints imported or overwritten."""
from __future__ import annotations
from pathlib import Path
from time import perf_counter
import json,sys,hashlib,platform,os
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import cholesky,solve_triangular
from scipy.stats import beta as beta_distribution
from modal import *
from dda import grid_target,solve as dda_solve
from coverage import Block,Cell,ModalProblem,cover,choose_offset

ROOT=Path(__file__).resolve().parents[1]

def jsonable(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,complex): return [x.real,x.imag]
    if isinstance(x,Cell): return {'center':x.center.tolist(),'half':x.half.tolist()}
    raise TypeError(type(x).__name__)

def save(name,data):
    path=ROOT/'results'/name;path.write_text(json.dumps(data,default=jsonable,indent=2,allow_nan=False));return path

def complex_noise(rng,shape):
    return (rng.normal(size=shape)+1j*rng.normal(size=shape))/np.sqrt(2)


def spectral_experiment(seed=715):
    rng=np.random.default_rng(seed);rows=[]
    for z in np.logspace(-2,2,41):
        n=rng.normal(size=3);n/=np.linalg.norm(n);k=10.;r=z/k*n
        num=visible_singular_values(r,k);theory=theoretical_singular_values(z/k,k)
        rows.append({'z':z,'numerical':num,'theory':theory,
                     'relative_error':np.linalg.norm(num-theory)/np.linalg.norm(theory)})
    out={'seed':seed,'rows':rows,'near_log_slope':float(np.polyfit(np.log([r['z'] for r in rows[:7]]),np.log([r['numerical'][2] for r in rows[:7]]),1)[0]),
         'far_log_slope':float(np.polyfit(np.log([r['z'] for r in rows[-7:]]),np.log([r['numerical'][2] for r in rows[-7:]]),1)[0]),
         'maximum_formula_relative_error':max(x['relative_error'] for x in rows)}
    return out


def dda_experiment():
    rows=[];checks=[]
    shapes={'sphere':(.05,.05,.05),'near_sphere':(.053,.05,.047),'ellipsoid':(.07,.045,.035),'box':(.05,.04,.035)}
    for shape,axes in shapes.items():
      for eps in [1.5+.02j,6.+.1j]:
       for k in [6.,18.,36.]:
        for h in [.025,.05/3]:
         points,v=grid_target(shape,h,axes);s=dda_solve(points,v,k,eps)
         for R in [.12,.4]:
          n=np.array([.5,.4,.7]);n/=np.linalg.norm(n);r=R*n
          y=s.field(r);q=np.linalg.norm(y);m=dyad(r,k).ravel()
          _,res=profile_complex(y,m)
          B=realify(s.geometry_jacobian(r)/q)
          mat=s.material_jacobian(r)[:,None]/q
          G=np.column_stack([y.ravel()/q,1j*y.ravel()/q])
          vis,rank=profile_columns(B,realify(np.column_stack([mat,G])))
          sv=np.linalg.svd(vis,compute_uv=False)
          row={'shape':shape,'epsilon':[eps.real,eps.imag],'k':k,'h':h,'R':R,'n_cells':len(points),
               'pattern_residual':float(np.linalg.norm(res)/q),'profiled_singular_values':sv,
               'sphere_formula_min':float(theoretical_singular_values(R,k)[2]),'nuisance_rank':rank,
               'linear_residual':s.relative_linear_residual,'setup_seconds':s.setup_seconds,'solve_seconds':s.solve_seconds,
               'cost_scope':'same state reused at two receivers; setup/solve repeated in rows, count once per shape/eps/k/h'}
          if shape=='sphere':
           ref=sphere_field(r,k,eps,.05)
           row['mie_relative_field_error']=float(np.linalg.norm(y-ref)/np.linalg.norm(ref))
          rows.append(row)
    return {'kind':'independent analytic Mie vs in-house DDA for spheres; shape extensions NOT continuum-certified','rows':rows}


def make_branch_scene(rng,kind):
    nominal=np.array([.13,.065,.12])+rng.uniform(-.015,.015,3)
    half=np.full(3,.04 if kind=='tangent_two_world' else .014); sign=1 if kind=='easy' else (-1 if kind=='bank_missing' else int(rng.choice([-1,1])))
    true=sign*nominal+rng.uniform(-.009,.009,3)
    roots=[Cell(nominal,half)] if kind=='easy' else [Cell(nominal,half),Cell(-nominal,half)]
    epsilon=rng.uniform(1.4,6.)+1j*rng.uniform(.01,.15);radius=rng.uniform(.025,.045)
    gain=(1+rng.uniform(-.15,.15))*np.exp(1j*rng.uniform(-.5,.5));clock=rng.uniform(-.02,.02)
    ks=[6.,12.,24.];q=18. if kind=='low_snr' else 400.
    state={}
    if kind=='model_discrepancy':
        for k in ks:
            pts,v=grid_target('ellipsoid',.02,(.055,.038,.03));state[k]=dda_solve(pts,v,k,epsilon)
    def mean_at(k,d):
        x=true+d
        # The deliberate two-world control places the data exactly on a wrong
        # admissible modal response, outside the beta=0 discrepancy assumption.
        if kind=='tangent_two_world': x=true+np.array([.025,0,0])+d
        field=state[k].field(x) if state else sphere_field(x,k,epsilon,radius)
        return gain*np.exp(1j*k*clock)*field.ravel()
    base=[mean_at(k,np.zeros(3)) for k in ks]
    sigmas=[np.linalg.norm(m)/q for m in base]
    train=[Block(m/s+complex_noise(rng,9),k,np.zeros(3)) for m,s,k in zip(base,sigmas,ks)]
    valid=[Block(m/s+complex_noise(rng,9),k,np.zeros(3)) for m,s,k in zip(base,sigmas,ks)]
    # Potential acquisitions and their noises are pre-generated and shared by policies.
    pool=np.array([[.035,0,0],[-.035,0,0],[0,.035,0],[0,-.035,0],[0,0,.035],[0,0,-.035]])
    potential=[];potential_v=[]
    for d in pool:
        m=mean_at(12.,d)/sigmas[1]
        potential.append(Block(m+complex_noise(rng,9),12.,d))
        potential_v.append(Block(m+complex_noise(rng,9),12.,d))
    return {'true':true,'nominal':nominal,'roots':roots,'train':train,'valid':valid,'pool':pool,
            'potential':potential,'potential_v':potential_v,'epsilon':epsilon,'radius':radius,'kind':kind,
            'noise_q':q,'gain':gain,'clock':clock}


def branch_experiment(seed,n_per_kind=4,max_cells=6000):
    rng=np.random.default_rng(seed);rows=[];raw={};scene_idx=0
    kinds=['easy','antipodal','bank_missing','low_snr','model_discrepancy','tangent_two_world']
    for kind in kinds:
      for repeat in range(n_per_kind):
        s=make_branch_scene(rng,kind);sid=f'{kind}_{repeat}';scene_idx+=1
        for j,b in enumerate(s['train']):raw[f'{sid}_train_{j}']=b.y
        for j,b in enumerate(s['valid']):raw[f'{sid}_validation_{j}']=b.y
        for j,b in enumerate(s['potential']):raw[f'{sid}_potential_{j}']=b.y
        for j,b in enumerate(s['potential_v']):raw[f'{sid}_potential_validation_{j}']=b.y
        raw[f'{sid}_true']=s['true']
        # Freeze common starting candidates before independent validation.
        base=ModalProblem(s['train']);bank,bcost=base.fit(s['roots']);positive=next(a for a in bank if s['roots'][0].contains(a['r']))
        for policy in ['naive_multistart','algebraic_multistart','fisher_acquisition','random_acquisition','branch_acquisition','coverage_aware']:
          t=perf_counter();P=ModalProblem(s['train']);V=ModalProblem(s['valid']);fitcost=bcost.copy();acq=0;design=0.;offset=np.zeros(3)
          pool_bank=bank if policy!='naive_multistart' else [positive]
          estimate=pool_bank[0]['r'];present=any(np.linalg.norm(a['r']-s['true'])<=.015 for a in pool_bank)
          if policy in ['fisher_acquisition','random_acquisition','branch_acquisition','coverage_aware']:
            td=perf_counter(); method='branch' if policy in ['branch_acquisition','coverage_aware'] else policy.split('_')[0]
            candidates=[a['r'] for a in bank]
            offset,scores=choose_offset(candidates,12.,s['pool'],method,rng)
            idx=int(np.argmin(np.linalg.norm(s['pool']-offset,axis=1)));design=perf_counter()-td;acq=9
            P=ModalProblem(s['train']+[s['potential'][idx]]);V=ModalProblem(s['valid']+[s['potential_v'][idx]])
            pool_bank,extra=P.fit(s['roots'],extra_starts=candidates)
            estimate=pool_bank[0]['r'];fitcost={'seconds':bcost['seconds']+extra['seconds'],'profile_evaluations':bcost['profile_evaluations']+extra['profile_evaluations']}
          # Selection has not looked at validation. This is a frozen-geometry
          # feasibility test with nuisance scalar reprofiled, not noise-free gains.
          valnorm=V.norm(estimate); valok=valnorm<=V.threshold(.005)
          err=float(np.linalg.norm(estimate-s['true']));correct=err<=.015
          extra_record={}
          if policy=='coverage_aware':
            cov=cover(P,s['roots'],estimate,tolerance=.015,alpha=.005,max_cells=max_cells)
            accepted=bool(cov['accepted'] and valok)
            extra_record={key:value for key,value in cov.items() if key!='cells'}
            truth_not_eliminated=any(c.contains(s['true']) for c in cov['cells'])
            extra_record['true_in_outer_regions']=truth_not_eliminated
            extra_record['outer_cells']=cov['cells']
          else: accepted=bool(valok and P.norm(estimate)<=P.threshold(.005))
          represented_final=any(np.linalg.norm(a['r']-s['true'])<=.015 for a in pool_bank)
          row={'scene_id':sid,'kind':kind,'policy':policy,'true':s['true'],'estimate':estimate,
               'geometry_error_m':err,'correct':correct,'accepted':accepted,'wrong_accepted':accepted and not correct,
               'correct_rejected':correct and not accepted,'initial_bank_missing':not present,
               'final_bank_covered':represented_final,'wrong_selection_when_covered':represented_final and not correct,
               'validation_residual':valnorm,'validation_passed':bool(valok),'offset':offset,
               'q':s['noise_q'],'new_complex_acquisitions':acq,'training_complex_data':27+acq,
               'validation_complex_data':27+acq,'electronics_references':0,'model':'full Maxwell electric-l=1',
               'fit_seconds':fitcost['seconds'],'fit_profile_evaluations':fitcost['profile_evaluations'],
               'design_seconds':design,'additional_wall_seconds':perf_counter()-t,'coverage':extra_record}
          rows.append(row)
    return {'seed':seed,'independent_scenes':scene_idx,'noise_replicates_per_scene':1,'rows':rows},raw


def rayleigh_coefficient(eps,k,a):return (2/3)*k*k*a**3*(eps-1)/(eps+2)


def calibration_forward(theta,ks,a,fidelity):
    r=theta[:3];eps=theta[3]+.05j;g=theta[4]+1j*theta[5];clock=theta[6]
    out=[]
    for k in ks:
        if fidelity=='fine':y=sphere_field(r,k,eps,a)
        else:y=rayleigh_coefficient(eps,k,a)*dyad(r,k)
        out.append(g*np.exp(1j*k*clock)*y.ravel())
    return np.concatenate(out)


def reference_forward(theta,ks):return (theta[4]+1j*theta[5])*np.exp(1j*np.asarray(ks)*theta[6])


def calibration_experiment(seed,n_scenes=8,training_samples=32):
    """Strong prior baselines for joint 3D pose, material, shared gain and clock.

    Every method gets identical independent pilot and noisy reference channels.
    EEM construction cost is charged per scene (no favorable amortization).
    Conditional EEM regresses discrepancy on theta, a standard Gaussian linear
    conditional-error surrogate. Neither learned covariance is a certificate.
    """
    rng=np.random.default_rng(seed);rows=[];raw={}
    ks=[4.,8.,22.];a=.045;methods=['low_only','raw_all','isotropic','rank_one','sampled_eem','conditional_eem','coarse_to_fine','full_fine','risk_controller']
    for scene in range(n_scenes):
        nominal=np.array([.13,.07,.11])+rng.uniform(-.008,.008,3)
        true=np.r_[nominal+rng.uniform(-.007,.007,3),rng.uniform(1.5,6.),rng.uniform(.85,1.15),rng.uniform(-.2,.2),rng.uniform(-.01,.01)]
        mean=calibration_forward(true,ks,a,'fine'); sigma=np.concatenate([np.full(9,np.linalg.norm(mean[9*i:9*i+9])/180.) for i in range(3)])
        noise=complex_noise(rng,27)*sigma; y=mean+noise; pilot=mean+complex_noise(rng,27)*sigma
        rsigma=.015;refs=reference_forward(true,ks)+rsigma*complex_noise(rng,3)
        pilot_refs=reference_forward(true,ks)+rsigma*complex_noise(rng,3)
        init=np.r_[nominal,2.7,1.,0.,0.];bounds=(np.r_[nominal-.022,1.15,.4,-.6,-.05],np.r_[nominal+.022,7.,1.6,.6,.05])
        counts={'coarse':0,'fine':0};local_start=perf_counter()
        def prediction(theta,fidelity): counts[fidelity]+=1;return calibration_forward(theta,ks,a,fidelity)
        def residual(theta,fidelity,selected,yobs,transform=None,mean_error=None,refobs=None):
            v=prediction(theta,fidelity)[selected]-yobs[selected]
            v=realify(v/sigma[selected])
            if mean_error is not None:v+=mean_error(theta)
            if transform is not None:v=transform@v
            return np.r_[v,realify((reference_forward(theta,ks)-(refs if refobs is None else refobs))/rsigma)]
        low=np.arange(18);allrows=np.arange(27)
        pilotfit=least_squares(lambda t:residual(t,'coarse',low,pilot,refobs=pilot_refs),init,bounds=bounds,max_nfev=180,x_scale='jac')
        common_init=pilotfit.x.copy();pilot_seconds=perf_counter()-local_start;pilot_counts=counts.copy()
        # Prior training uses declared nominal/bounds, never true theta or test noise.
        tc=perf_counter();thetas=[];errors=[]
        for _ in range(training_samples):
            t=np.r_[nominal+rng.uniform(-.014,.014,3),rng.uniform(1.2,7.),rng.uniform(.7,1.3),rng.uniform(-.3,.3),rng.uniform(-.03,.03)]
            thetas.append(t);errors.append(realify((prediction(t,'fine')-prediction(t,'coarse'))/sigma))
        X=np.asarray(thetas);E=np.asarray(errors);emean=E.mean(0);ec=E-emean
        cov=ec.T@ec/(training_samples-1);m=len(emean);I=np.eye(m)
        # Only downweight high-frequency directions for isotropic/rank-one controls.
        hi=np.r_[np.arange(18,27),np.arange(45,54)]
        hp=np.zeros(m);hp[hi]=1.;masked=ec*hp
        hc=masked.T@masked/(training_samples-1);ev,eu=np.linalg.eigh(hc)
        covs={'raw_all':np.zeros_like(cov),'isotropic':np.diag(hp*np.trace(hc)/18),
              'rank_one':ev[-1]*np.outer(eu[:,-1],eu[:,-1]),'sampled_eem':cov}
        xc=X-X.mean(0);xs=xc.std(0);xr=xc/xs
        reg=np.linalg.lstsq(np.r_[xr,.2*np.eye(7)],np.r_[ec,np.zeros((7,m))],rcond=None)[0]
        cres=ec-xr@reg;ccov=cres.T@cres/(training_samples-1)
        covs['conditional_eem']=ccov
        transforms={key:solve_triangular(cholesky(I+C,lower=True),I,lower=True) for key,C in covs.items()}
        construction_seconds=perf_counter()-tc
        construction_counts={key:counts[key]-pilot_counts[key] for key in counts}
        before_screen=counts.copy()
        # A frozen surrogate controller: evaluate sample-model geometry/material
        # task risk at the independent pilot. Mean and covariance are from the
        # training model family; all candidate fine calls are charged.
        tr=perf_counter();t0=common_init
        def jfinite(f):
            out=[]
            for j in range(7):
                step=1e-5*max(abs(t0[j]),.01);v=np.zeros(7);v[j]=step
                out.append((f(t0+v)-f(t0-v))/(2*step))
            return np.column_stack(out)
        from risk import weighted_task_risk
        Rj=jfinite(lambda t:realify(reference_forward(t,ks)/rsigma))
        jc=jfinite(lambda t:realify(prediction(t,'coarse')/sigma));jf=jfinite(lambda t:realify(prediction(t,'fine')/sigma))
        Caug=np.zeros((60,60));Caug[:54,:54]=cov;ef=np.linalg.cholesky(Caug+1e-12*np.eye(60))
        T=np.zeros((4,7));T[:3,:3]=np.eye(3)/.015;T[3,3]=1/.15
        scores={}
        for key in ['raw_all','isotropic','rank_one','sampled_eem']:
            W=np.eye(60);W[:54,:54]=transforms[key].T@transforms[key]
            a_score=weighted_task_risk(np.r_[jc,Rj],W,T,np.r_[emean if key!='sampled_eem' else np.zeros(54),np.zeros(6)],ef)
            scores[key]=a_score['risk']
        scores['full_fine']=weighted_task_risk(np.r_[jf,Rj],np.eye(60),T)['risk']
        choice=min(scores,key=scores.get);screen_seconds=perf_counter()-tr
        raw[f's{scene}_truth']=true;raw[f's{scene}_y']=y;raw[f's{scene}_pilot']=pilot;raw[f's{scene}_references']=refs;raw[f's{scene}_pilot_references']=pilot_refs;raw[f's{scene}_training_theta']=X;raw[f's{scene}_training_errors']=E
        screen_counts={key:counts[key]-before_screen[key] for key in counts}
        for method in methods:
            start_counts=counts.copy();t=perf_counter();selected=low if method=='low_only' else allrows
            target=choice if method=='risk_controller' else method
            fidelity='fine' if target in ['full_fine','coarse_to_fine'] else 'coarse'
            transform=transforms.get(target)
            errorfun=(lambda t:emean) if target=='sampled_eem' else ((lambda t:emean+((t-X.mean(0))/xs)@reg) if target=='conditional_eem' else None)
            x0=common_init.copy()
            if target=='coarse_to_fine':
                warm=least_squares(lambda t:residual(t,'coarse',allrows,y),x0,bounds=bounds,max_nfev=180,x_scale='jac');x0=warm.x
            result=least_squares(lambda t:residual(t,fidelity,selected,y,transform,errorfun),x0,bounds=bounds,max_nfev=250,x_scale='jac',ftol=1e-9,xtol=1e-9,gtol=1e-9)
            elapsed=perf_counter()-t;est=result.x
            modecost=construction_seconds if method in ['isotropic','rank_one','sampled_eem','conditional_eem','risk_controller'] else 0.
            extra=screen_seconds if method=='risk_controller' else 0.
            counts_cost={key:counts[key]-start_counts[key]+pilot_counts[key]+(construction_counts[key] if modecost else 0)+(screen_counts[key] if extra else 0) for key in counts}
            ge=float(np.linalg.norm(est[:3]-true[:3]));me=float(abs(est[3]-true[3])/true[3]);loss=(ge/.015)**2+((est[3]-true[3])/.15)**2
            row={'scene':scene,'method':method,'selected_action':('refine_model' if target=='full_fine' else ('downweight' if target in transforms and target!='raw_all' else 'retain')),
                 'selected_method':target,'true':true,'estimate':est,'geometry_error_m':ge,'material_relative_error':me,'scaled_task_loss':float(loss),
                 'fit_success':bool(result.success),'nfev':result.nfev,'fit_seconds':elapsed,'pilot_seconds':pilot_seconds,'mode_seconds':modecost,
                 'screen_seconds':extra,'total_seconds':elapsed+pilot_seconds+modecost+extra,
                 'forward_calls_including_construction':counts_cost,'training_data_complex':len(selected),'pilot_complex':18,
                 'reference_complex':3,'pilot_reference_complex':3,'offline_error_sample_pairs':training_samples if modecost else 0,
                 'risk_scores':{key:(float(value) if np.isfinite(value) else None) for key,value in scores.items()} if method=='risk_controller' else {},
                 'risk_status':'surrogate, not calibrated coverage','covariance_construction_scope':'per scene, no amortization'}
            rows.append(row)
    return {'seed':seed,'independent_scenes':n_scenes,'rows':rows,'comparison_scope':'restricted radial full-wave Mie target vs Rayleigh coarse model; no A3 reproduction'},raw


def summarize_branch(out):
    rows=out['rows'];table=[]
    for policy in sorted({r['policy'] for r in rows}):
        rr=[r for r in rows if r['policy']==policy];N=len(rr);wa=sum(r['wrong_accepted'] for r in rr);cov=[r for r in rr if r['final_bank_covered']]
        table.append({'policy':policy,'scenes':N,'accepted':sum(r['accepted'] for r in rr),'correct':sum(r['correct'] for r in rr),
                      'wrong_accepted':wa,'wrong_accept_CP_upper_one_sided_95':float(beta_distribution.ppf(.95,wa+1,N-wa)) if wa<N else 1.,
                      'correct_rejected':sum(r['correct_rejected'] for r in rr),'initial_bank_missing':sum(r['initial_bank_missing'] for r in rr),
                      'covered_scenes':len(cov),'wrong_selected_given_covered_count':sum(r['wrong_selection_when_covered'] for r in cov),
                      'median_fit_seconds':float(np.median([r['fit_seconds'] for r in rr])),
                      'median_coverage_seconds':float(np.median([r['coverage'].get('seconds',0) for r in rr])),
                      'median_error_m':float(np.median([r['geometry_error_m'] for r in rr]))})
    return table

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['dev','frozen']);args=p.parse_args()
    if args.mode=='dev':
        a,raw=branch_experiment(81031,1,max_cells=2000);save('development_branch.json',a);np.savez_compressed(ROOT/'results/development_branch_data.npz',**raw)
        b,raw=calibration_experiment(82031,2,16);save('development_calibration.json',b);np.savez_compressed(ROOT/'results/development_calibration_data.npz',**raw)
        print(json.dumps({'branch':summarize_branch(a),'calibration_rows':len(b['rows'])},indent=2))
    else:
        manifest=json.loads((ROOT/'provenance/frozen_manifest.json').read_text())
        for rel,digest in manifest['source_hashes'].items():
            if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=digest:raise RuntimeError('Frozen source mismatch: '+rel)
        if (ROOT/'results/frozen_complete.json').exists():raise RuntimeError('Refusing to overwrite completed frozen results')
        t=perf_counter();save('frozen_spectral.json',spectral_experiment(91031));save('frozen_dda.json',dda_experiment())
        a,raw=branch_experiment(92031,4,max_cells=6000);save('frozen_branch.json',a);np.savez_compressed(ROOT/'results/frozen_branch_data.npz',**raw)
        b,raw=calibration_experiment(93031,8,32);save('frozen_calibration.json',b);np.savez_compressed(ROOT/'results/frozen_calibration_data.npz',**raw)
        save('frozen_complete.json',{'seconds':perf_counter()-t,'source_manifest_sha256':hashlib.sha256((ROOT/'provenance/frozen_manifest.json').read_bytes()).hexdigest(),
                                    'branch_scenes':a['independent_scenes'],'branch_rows':len(a['rows']),'calibration_scenes':b['independent_scenes'],'calibration_rows':len(b['rows']),
                                    'environment':{'platform':platform.platform(),'python':sys.version,'numpy':np.__version__,'threads':os.getenv('OPENBLAS_NUM_THREADS')}})
        print('FROZEN COLLECTION COMPLETE; no algorithm changes permitted for these results')
