"""M2: finite bounded-amplitude windows, not an all-coherent impossibility result.

Necessary/sufficient statement and proof are in MANUSCRIPT_ZH.md. This module
exports directed interval budgets and a numerical implementation of the
coordinatewise minimax estimator. The estimator's root routine is diagnostic;
only the interval trace and analytic inequalities supply certified bounds.
"""
from __future__ import annotations
from pathlib import Path
from math import factorial
import json,hashlib,time,argparse
import numpy as np
from scipy.optimize import brentq
from modal_interval import Q,iv,rat,enc,bounds,export,series,outer,falling,enclose
from multipole import mie
AMIN=Q(3,4);AMAX=Q(5,4);DELTA=Q(1,10)
REGIONS=[(Q(3,2),Q(4),Q(1,5)),(Q(2),Q(5),Q(3,20))]
SIGMA=[Q('0.0000026664081'),Q('0.0000018047986')]
ZA=Q(3290526731492,10**12)


def cdf_enclose(x):
    """Gaussian CDF from power series with a uniform geometric tail."""
    total=rat(0);xx=rat(x)
    for k in range(100):
        c=Q((-1)**k,(2*k+1)*2**k*factorial(k))
        total+=rat(c)*xx**(2*k+1)
    first=xx**201/rat(Q(201*2**100*factorial(100)))
    assert x*x/Q(2*101)<Q(1,2)
    r=bounds(2*first)[1]
    return rat(Q(1,2))+(total+enc(-r,r))/iv.sqrt(2*iv.pi)


def complex_series(e,x,kind,derivative=0,xderivative=0,terms=18):
    """Same entire series, enclosing eps in a complex strip |eps| <= 5.001."""
    k0=next(k for k in range(terms) if falling(k+1 if kind=='J' else k,derivative) and falling(2*k+1 if kind=='J' else 2*k,xderivative))
    def coeff(k):
        ep=k+1 if kind=='J' else k;xp=2*k+1 if kind=='J' else 2*k
        c=Q(2*(k+1),factorial(2*k+3))
        if kind=='D':c*=2*k+2
        return c*falling(ep,derivative)*falling(xp,xderivative)
    v=rat((-1)**(terms-1)*coeff(terms-1));s=e*x*x
    for k in range(terms-2,k0-1,-1):v=rat((-1)**k*coeff(k))+s*v
    v*=e**((k0+1 if kind=='J' else k0)-derivative)*x**((2*k0+1 if kind=='J' else 2*k0)-xderivative)
    E=Q(5001,1000);X=bounds(x)[1];k=terms
    assert bounds(abs(e))[1]<E and X<Q(201,1000)
    ratio=Q(2*(k+2),k+1)/((2*k+4)*(2*k+5))*Q(k+3,k-1)**2*E*X*X
    assert ratio<Q(1,2)
    ep=k+1 if kind=='J' else k;xp=2*k+1 if kind=='J' else 2*k
    radius=2*coeff(k)*E**(ep-derivative)*X**(xp-xderivative)
    return v+iv.mpc(enc(-radius,radius),enc(-radius,radius))


def physical_derivatives(erange,x0):
    x=enc(x0*(1-Q(3,10000)),x0*(1+Q(3,10000)))
    e=iv.mpc(enc(*erange),enc(0,Q(1,1000)))
    j,dj,jx,djx=outer(*bounds(x),iv.dps)
    y=-iv.cos(x)/x**2-iv.sin(x)/x
    dy=iv.sin(x)/x**2+iv.cos(x)/x**3-iv.cos(x)/x
    yx=2*iv.sin(x)/x**2+2*iv.cos(x)/x**3-iv.cos(x)/x
    dyx=iv.sin(x)/x+2*iv.cos(x)/x**2-3*iv.sin(x)/x**3-3*iv.cos(x)/x**4
    J=complex_series(e,x,'J');D=complex_series(e,x,'D')
    Je=complex_series(e,x,'J',1);De=complex_series(e,x,'D',1)
    Jx=complex_series(e,x,'J',xderivative=1);Dx=complex_series(e,x,'D',xderivative=1)
    N=J*dj-j*D;den=J*dy-y*D
    Ne=Je*dj-j*De;dene=Je*dy-y*De
    Nx=Jx*dj+J*djx-jx*D-j*Dx;denx=Jx*dy+J*dyx-yx*D-y*Dx
    I=iv.mpc(0,1);W=N+I*den
    assert bounds(abs(W))[0]>0
    te=abs((Ne*den-N*dene)/W**2)
    tx=abs((Nx*den-N*denx)/W**2)
    return te,tx


def hfloat(e,x): return float(abs(mie(e,x,1)[3]))


def feasible_intervals(r,z,B,Ba):
    """Exact formula, numerical scalar inversion. Raises on empty compatibility."""
    r=np.asarray(r,float);B=np.asarray(B,float)
    if r.shape!=(2,) or B.shape!=(2,) or np.any(B<0) or Ba<0:raise ValueError('invalid amplitude bounds')
    lower=[];upper=[]
    for lo,hi,x in REGIONS:
        lower.append(hfloat(float(lo),float(x)));upper.append(hfloat(float(hi),float(x)))
    amin=max(float(AMIN),z-Ba,*((r-B)/upper));amax=min(float(AMAX),z+Ba,*((r+B)/lower))
    if amin>amax or amax<=0:raise ValueError('empty compatibility set')
    intervals=[]
    for i,(lo,hi,x) in enumerate(REGIONS):
        fun=lambda e:hfloat(e,float(x))
        def inv(value):
            if value<=lower[i]:return float(lo)
            if value>=upper[i]:return float(hi)
            return float(brentq(lambda e:fun(e)-value,float(lo),float(hi),xtol=1e-12))
        intervals.append([inv(max(lower[i],(r[i]-B[i])/amax)),inv(min(upper[i],(r[i]+B[i])/amin))])
    v=np.array(intervals)
    return {'gain_interval':[amin,amax],'material_intervals':intervals,'estimate':v.mean(axis=1).tolist(),'radius':((v[:,1]-v[:,0])/2).tolist()}


def require_class_c_certificate(covered_lower_bound,model_error_upper_bound,noise_radius):
    """Fail closed. None is not zero; an optimizer minimum is not a lower bound."""
    if covered_lower_bound is None or model_error_upper_bound is None:
        return {'status':'unresolved','reason':'covered separation and a justified model-error bound are both required'}
    if min(covered_lower_bound,model_error_upper_bound,noise_radius)<0:raise ValueError('negative norm bound')
    margin=covered_lower_bound-2*(model_error_upper_bound+noise_radius)
    return {'status':'sufficient_margin' if margin>0 else 'certificate_not_established','margin':margin,
            'note':'inputs must already be valid covered bounds; this routine cannot certify them'}


def main(out,resume=False):
    iv.dps=60;out.mkdir(parents=True,exist_ok=True)
    trace=out/'finite_cover.jsonl'
    cached={}
    if trace.exists():
        if not resume:raise FileExistsError(trace)
        lines=trace.read_text().splitlines();valid=[]
        for j,line in enumerate(lines):
            try: rec=json.loads(line)
            except json.JSONDecodeError:
                if j!=len(lines)-1:raise
                trace.with_name('incomplete_tail.txt').write_text(line);break
            key=(rec['kind'],rec['channel'],tuple(rec['box']))
            if key in cached:raise ValueError('duplicate trace box')
            cached[key]=rec;valid.append(line)
        trace.write_text('\n'.join(valid)+'\n')
    start=time.perf_counter();cdf=cdf_enclose(ZA)
    assert bounds(cdf)[0]>=Q(9995,10000)
    ba_noise=ZA*Q(1,1000)
    bi_noise=[bounds(rat(s)*iv.sqrt(iv.ln(1000)))[1] for s in SIGMA]
    result={'scope':'necessary/sufficient only for the stated bounded-amplitude experiment',
            'delta':str(DELTA),'normal_quantile_upper':str(ZA),'CDF_at_upper':export(cdf),'event_failure_union_bound':'0.003',
            'noise_bounds':{'modal':[str(x) for x in bi_noise],'real_reference':str(ba_noise)},'regions':[]}
    seen=set()
    with trace.open('a' if resume else 'w') as log:
        for i,(lo,hi,x) in enumerate(REGIONS):
            n=int((hi-lo-2*DELTA)*1000);maxdiff=None;prev=lo
            for k in range(n):
                a=lo+Q(k,1000);b=a+Q(1,1000);assert a==prev;prev=b
                key=('finite_increment_derivative',i,(str(a),str(b)));seen.add(key)
                if key in cached:upper=Q(cached[key]['upper'])
                else:
                    hp0=enclose(enc(a,b),rat(x))['hp'];hp1=enclose(enc(a+2*DELTA,b+2*DELTA),rat(x))['hp']
                    upper=bounds(hp1-hp0)[1]
                    log.write(json.dumps({'kind':key[0],'channel':i,'box':[str(a),str(b)],'upper':str(upper)},separators=(',',':'))+'\n')
                assert upper<0
                maxdiff=upper if maxdiff is None else max(maxdiff,upper)
            assert prev==hi-2*DELTA
            h0=enclose(rat(lo),rat(x))['h'];hm=enclose(rat(hi-2*DELTA),rat(x))['h'];h1=enclose(rat(hi),rat(x))['h']
            # This overlap condition makes the OTHER material uninformative in the
            # least-favourable amplitude pair used for the necessity direction.
            assert bounds(h1/h0)[0]>AMAX/AMIN
            le=Q(0);lx=Q(0)
            count=int((hi-lo)*100)
            for k in range(count):
                a=lo+Q(k,100);b=a+Q(1,100)
                key=('complex_strip_and_size',i,(str(a),str(b)));seen.add(key)
                if key in cached:
                    eu=Q(cached[key]['abs_te'][1]);xu=Q(cached[key]['abs_tx'][1])
                else:
                    te,tx=physical_derivatives((a,b),x);eu=bounds(te)[1];xu=bounds(tx)[1]
                    log.write(json.dumps({'kind':key[0],'channel':i,'box':[str(a),str(b)],'abs_te':export(te),'abs_tx':export(tx)},separators=(',',':'))+'\n')
                le=max(le,eu);lx=max(lx,xu)
            hcap=bounds(h1)[1]
            leakage=AMAX*hcap*Q(2,1000)
            channel_gain=AMAX*hcap*Q(1,1000)
            size=AMAX*lx*x*Q(3,10000)
            loss=AMAX*le*Q(1,1000)
            mixed=Q(1,1000)*(leakage+size+loss)
            physical=leakage+channel_gain+size+loss+mixed
            scenarios={}
            for repeats,bias in [(1,Q(0)),(16,Q(5,10000))]:
                root=1 if repeats==1 else 4
                Ba=ba_noise/root+bias;d=min(2*Ba,AMAX-AMIN)
                margin=rat(AMIN)*(h1-hm)-rat(d)*hm
                allowed=margin/2-rat(bi_noise[i]/root)
                old_m=[Q('0.000429'),Q('0.000129')][i];old_H=[Q('0.00272'),Q('0.00131')][i]
                inherited=DELTA*(AMIN-Ba)*old_m-bi_noise[i]/root-old_H*Ba
                totalB=bi_noise[i]/root+physical
                dstar=(rat(AMIN)*(h1-hm)-rat(2*totalB))/hm
                scenarios[str(repeats)]={'repeats':repeats,'reference_bias':str(bias),'total_reference_error':str(Ba),
                    'finite_allowed_model_bias':export(allowed),'finite_allowed_model_bias_display':[float(z) for z in bounds(allowed)],
                    'inherited_ratio_allowed_model_bias':str(inherited),
                    'fixed_physical_budget':str(physical),'fixed_physical_budget_display':float(physical),
                    'physical_budget_pass':bool(physical<=bounds(allowed)[0]),
                    'max_total_reference_error_for_this_B':export(dstar/2)}
            result['regions'].append({'channel':i,'size':str(x),'material':[str(lo),str(hi)],'finite_cover_cells':n,
                'finite_derivative_upper':str(maxdiff),'endpoint_minimum':'proved for F_d(e)=amin*(h(e+0.2)-h(e))-d*h(e), all d>=0',
                'other_channel_overlap_ratio':export(h1/h0),'physical_strip_cells':count,'abs_te_max':str(le),'abs_tx_max':str(lx),
                'bias_components_display':{'modal_leakage':float(leakage),'relative_channel_gain':float(channel_gain),'size':float(size),'loss':float(loss),'mixed_terms':float(mixed)},
                'scenarios':scenarios})
    assert set(cached).issubset(seen)
    # Fixed unit-test stream, distinct from EM recovery experiments.
    rng=np.random.default_rng(2026091117);maxratio=0.
    B=np.array([float(bi_noise[i]/4+Q(result['regions'][i]['scenarios']['16']['fixed_physical_budget'])) for i in range(2)])
    Ba=float(ba_noise/4+Q(5,10000));empty_rejected=False
    for _ in range(100):
        truth=rng.uniform([1.5,2.],[4.,5.]);gain=rng.uniform(.75,1.25)
        rr=gain*np.array([hfloat(truth[i],float(REGIONS[i][2])) for i in range(2)])+rng.uniform(-B,B)
        zz=gain+rng.uniform(-Ba,Ba);fit=feasible_intervals(rr,zz,B,Ba)
        assert np.all(np.abs(np.array(fit['estimate'])-truth)<=np.array(fit['radius'])+1e-10)
        maxratio=max(maxratio,max(fit['radius'])/.1)
        assert max(fit['radius'])<=.1+1e-10
    try:feasible_intervals([1.,1.],.75,B,Ba)
    except ValueError:empty_rejected=True
    assert empty_rejected
    result['unit_checks']={'seed':2026091117,'valid_error_box_cases':100,'maximum_radius_over_target':maxratio,
        'empty_set_rejected':empty_rejected,'missing_class_C_bound':require_class_c_certificate(None,None,1.)}
    result['trace_sha256']=hashlib.sha256(trace.read_bytes()).hexdigest();result['seconds']=time.perf_counter()-start
    (out/'finite_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({**result,'regions':[{k:v for k,v in r.items() if k in ['channel','bias_components_display','scenarios']} for r in result['regions']]},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path(__file__).parent/'results/m2');p.add_argument('--resume',action='store_true');args=p.parse_args();main(args.out,args.resume)
