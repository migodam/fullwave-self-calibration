"""Frozen-protocol analysis; refuses incomplete final samples."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta,binom

OUT=Path(__file__).resolve().parent


def cp_upper(x,n,alpha=.025):
    return 1. if x==n else float(beta.ppf(1-alpha,x+1,n-x))


def cp_lower(x,n,alpha=.025):
    return 0. if x==0 else float(beta.ppf(alpha,x,n-x+1))


def median_interval(values):
    a=np.sort(values);n=len(a)
    choices=[k for k in range(1,n//2+1) if 2*binom.cdf(k-1,n,.5)<=.05]
    if not choices:return [0.,None]
    k=max(choices)
    return [float(a[k-1]),float(a[n-k])]


def analyze():
    rows=json.loads((OUT/'results/rom_final.json').read_text())
    methods=['direct_gn','direct_adjoint','generic_task','sensing_task','twofold_task','krylov_task']
    required={(s,m) for s in range(9201,9241) for m in methods}
    assert len(rows)==240 and {(r['seed'],r['method']) for r in rows}==required, 'final sample incomplete'
    assert len({r['protocol_sha256'] for r in rows})==1 and len({r['source_digest'] for r in rows})==1
    indexed={(r['seed'],r['method']):r for r in rows}
    wins=losses=0;ratios=[];eligible=[]
    for seed in range(9201,9241):
        a=indexed[seed,'twofold_task'];b=indexed[seed,'krylov_task']
        wins+=int(a['joint_success'] and not b['joint_success'])
        losses+=int(b['joint_success'] and not a['joint_success'])
        if a['joint_success'] and b['joint_success']:
            gap=abs(a['final_loss']-b['final_loss'])/max(abs(b['final_loss']),1e-14)
            if gap<=1e-4:
                ratios.append(a['actual_consumed_seconds']/b['actual_consumed_seconds'])
                eligible.append(seed)
    upper=cp_upper(wins,40)-cp_lower(losses,40)
    lower=cp_lower(wins,40)-cp_upper(losses,40)
    pvalue=float(binom.cdf(sum(r<=.8 for r in ratios),len(ratios),.5)) if ratios else None
    descriptive=[]
    for method in methods:
        group=[r for r in rows if r['method']==method]
        descriptive.append(dict(method=method,joint_successes=sum(r['joint_success'] for r in group),n=40,
            exceptions=sum(r['status']=='exception' for r in group),
            mean_seconds=float(np.mean([r['actual_consumed_seconds'] for r in group])),
            max_atomic_overrun=max(r.get('atomic_overrun_seconds',0.) for r in group),
            deadline_exits=sum(r['status']=='decision_deadline' for r in group),
            accepted_reduced_updates=sum(r.get('reduced_updates',0) for r in group)))
    for row in rows:
        assert all(e['seconds']<8. for e in row.get('events',[]))
    result=dict(protocol_sha256=rows[0]['protocol_sha256'],source_digest=rows[0]['source_digest'],
        n=40,descriptive=descriptive,
        primary_success=dict(twofold_only_wins=wins,krylov_only_wins=losses,
            observed_difference=(wins-losses)/40,one_sided_95_upper=upper,
            one_sided_95_lower=lower,excludes_ten_point_advantage=upper<.10),
        primary_speed=dict(eligible_n=len(ratios),eligible_seeds=eligible,paired_ratios=ratios,
            conditional_median_ratio=float(np.median(ratios)) if ratios else None,
            exact_95_median_interval=median_interval(ratios),
            p_value_against_median_ratio_at_most_point8=pvalue,
            excludes_conditional_median_twenty_percent_saving=pvalue is not None and pvalue<.05),
        scope='bounded same-discretization q9 population; success at decision deadline; speed conditional on successful objective-matched pairs; atomic overrun disclosed',
        inference_assumptions='Independent sampled scenes; runtime interpretation also assumes sufficiently stationary execution conditions. No broad SOM-family or high-dimensional claim.')
    (OUT/'results/rom_final_analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':analyze()
