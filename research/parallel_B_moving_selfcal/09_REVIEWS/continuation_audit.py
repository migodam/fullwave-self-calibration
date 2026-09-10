"""Limited fresh audit, NOT a replay of missing B1 recovery arrays/experiment.py.
Run: python continuation_audit.py --output NEW_DIRECTORY
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '05_CODE'))
import model

EXPECTED = {
    'model.py':'738251fa5462ea79a6671a71f4dcfeeb40068e7f3d0c9c18ad466a330e7af17a',
    'vsw.py':'8e960c2052fa126427ec203eeefabc955e6ea554369dda7af3cf5fbbf6dfe509',
    'dda.py':'bb401b4719e5e74e4dd5614ca566e85294e75e32c9d893540eb75e289958030b',
}


def relative(x, y):
    return float(np.linalg.norm(x-y) / max(np.linalg.norm(y), 1e-300))


def run(output: Path):
    output.mkdir(parents=True, exist_ok=False)
    result = {'scope':'Fresh limited numerical audit; NOT original 36-world recovery replay',
              'base_commit':'685532825a15930b06656bee52649c8ab5a18edf',
              'seed':2026091199,'checks':{},'failures':[],
              'missing_originals':['experiment.py','Moving_Array_B1_Research.zip',
                                   'raw S/I/X observations and optimizer starts']}
    start = time.perf_counter()
    def record(name, value, passed=True):
        result['checks'][name]={'value':value,'passed':bool(passed)}
        (output/'AUDIT.json').write_text(json.dumps(result,indent=2))
        print(name,value,'PASS' if passed else 'FAIL',flush=True)
    try:
        hashes={n:hashlib.sha256((HERE.parent/'05_CODE'/n).read_bytes()).hexdigest() for n in EXPECTED}
        record('source_hashes', hashes,hashes==EXPECTED)
        if hashes!=EXPECTED:
            raise RuntimeError('Source byte mismatch; refusing to attribute results to recorded source')
        rng=np.random.default_rng(result['seed'])
        s=np.array([2.,2.7]); factor=1.12
        z=np.linspace(-1,1,model.T)
        pp=np.concatenate([model.points(t,z[t]) for t in range(model.T)])
        inc=np.stack([model.incident(pp,k) for k in model.K])
        record('direct_cross_max_abs',float(np.max(np.abs(inc[:,:,0,:]))),np.max(np.abs(inc[:,:,0,:]))==0)
        valid=np.all(factor*s>=model.LOWER) and np.all(factor*s<=model.UPPER) and .65<=1/factor<=1.35
        record('finite_pair_domain_and_disjoint_material_targets',{'s':s.tolist(),'partner':(factor*s).tolist(),
            'gain_partner':1/factor,'material_changes':((factor-1)*s).tolist()},valid and np.all((factor-1)*s>.2))
        tau=np.linspace(-10,10,model.T); g=np.exp(.03+1j*np.linspace(-.1,.1,model.T))
        alias=relative(model.gain(g,tau+2500),model.gain(g,tau))
        record('delay_alias_relative_error',alias,alias<1e-12)
        q=3.; interval=4.; u=np.linspace(0,interval,101)
        cond=q*u-(q*u)**2/(q*interval)
        record('brownian_bridge_midpoint_variance',float(cond[50]),abs(cond[50]-q*interval/4)<1e-14)
        alpha=np.array([.19+.02j,.25+.04j]); coupling=.4+.17j; cg=.9+.2j
        inv=np.diag(1/alpha)-np.array([[0,coupling],[coupling,0]])
        response=cg*np.linalg.inv(inv); v=np.linalg.inv(response)
        recovered=-coupling/v[0,1]; ar=1/(recovered*np.diag(v))
        err=max(abs(recovered-cg),np.max(abs(ar-alpha)))
        record('conditional_interaction_identity_error',float(err),err<1e-12)
        eta=1e-6; perturb=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        perturb*=eta/np.linalg.norm(perturb,2)
        vv=np.linalg.inv(response+perturb); nr=np.linalg.norm(v,2)
        beta=nr*nr*eta/(1-nr*eta)
        bound=abs(cg)**2*beta/(abs(coupling)-abs(cg)*beta)
        actual=abs(-coupling/vv[0,1]-cg)
        record('conditional_interaction_noise_bound',{'actual':float(actual),'upper':float(bound)},actual<=bound)
        vector=rng.normal(size=31)+1j*rng.normal(size=31)
        data=rng.normal(size=31)+1j*rng.normal(size=31)
        g0=np.vdot(vector,data)/np.vdot(vector,vector).real
        gp=np.clip(abs(g0),.65,1.35)*np.exp(1j*np.angle(g0))
        trials=np.linspace(.65,1.35,21)[:,None]*np.exp(1j*np.linspace(-np.pi,np.pi,1001))[None,:]
        a=np.vdot(vector,vector).real; b=np.vdot(vector,data)
        j=lambda gg:np.vdot(data,data).real+a*abs(gg)**2-2*np.real(np.conj(gg)*b)
        record('annulus_profile_vs_dense_grid',{'profile':float(j(gp)),'grid':float(np.min(j(trials)))},j(gp)<=np.min(j(trials))+1e-12)
        born=model.Born()
        bf=born.field(s,z); bf2=born.field(s*factor,z)/factor
        br=relative(bf2[:,:,:,0,:],bf[:,:,:,0,:])
        direct_change=relative(bf2[:,:,:,2,:],bf[:,:,:,2,:])
        record('born_cross_scale_residual',br,br<1e-12)
        record('nonzero_copolar_path_breaks_that_scale_pair',direct_change,direct_change>1e-3)
        del born
        fw=model.Forward(order=4)
        f=fw.field(s,z); exact=fw.field(s,z,exact=True)
        interpolation=relative(f[:,:,:,0,:],exact[:,:,:,0,:])
        record('cross_interpolation_relative_error',interpolation,interpolation<1e-9)
        f0=fw.cross(s,np.zeros(model.T)); f1=fw.cross(s*factor,np.zeros(model.T))
        gains=np.sum(f1.conj()*f0,axis=(1,2,3))/np.sum(abs(f1)**2,axis=(1,2,3))
        scale_residual=relative(gains[:,None,None,None]*f1,f0)
        # Diagnostic only: not a covered lower bound after unknown geometry profiling.
        record('fullwave_fixed_geometry_scale_pair_residual',scale_residual,scale_residual>1e-8)
        allf=fw.field(s,np.zeros(model.T))
        ratio=float(np.linalg.norm(allf[:,:,:,2,:])/np.linalg.norm(allf[:,:,:,0,:]))
        record('nominal_copolar_cross_rms_ratio',ratio)
        order6=model.Forward(order=6)
        f6=order6.cross(s,np.zeros(model.T))
        order_error=relative(f0,f6)
        record('cross_order4_vs6_relative_difference',order_error,order_error<1e-4)
        del order6
        independent=model.Independent(spacing=.0065)
        dda=independent.field(s,np.zeros(model.T))[:,:,:,0,:]
        disagreement=relative(dda,f0)
        record('DDA_VSW_cross_disagreement_DIAGNOSTIC',disagreement)
        np.savez_compressed(output/'fresh_fields.npz',s=s,dz=z,born_cross=bf[:,:,:,0,:],
                            born_partner=bf2[:,:,:,0,:],vsw_cross=f0,vsw_partner=f1,
                            dda_cross=dda,vsw_order6=f6,gain_fit=gains)
    except Exception:
        result['failures'].append(traceback.format_exc())
    result['elapsed_seconds']=time.perf_counter()-start
    result['all_checks_passed']=not result['failures'] and all(c['passed'] for c in result['checks'].values())
    result['n_checks']=len(result['checks'])
    result['fresh_recovery_runs']=0
    result['original_recovery_counts_recomputed_from_raw']=False
    result['audit_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (output/'AUDIT.json').write_text(json.dumps(result,indent=2))
    if not result['all_checks_passed']:
        print(result['failures'])
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.output)
