"""Read-only E4 records analysis; no test-dependent estimator modification."""
import json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=ROOT/'research/delegated/a2_solver_audited'
DEST=HERE/'results';FIG=HERE/'figures'
def loadrows(p):return [json.loads(x) for x in p.read_text().splitlines() if x.strip()]
def analyze(budget):
    rows=loadrows(SOURCE/f'final_{budget}.jsonl');frozen=json.loads((SOURCE/'frozen_parent.json').read_text())
    methods=frozen['methods'];seeds=frozen['seeds_final'];ref=frozen['reference_design']
    assert len(rows)==1200
    by={(r['seed'],r['init_index'],r['method']):r for r in rows};assert len(by)==len(rows)
    assert set(by)=={(s,i,m) for s in seeds for i in range(12) for m in methods}
    assert all(r['budget']==budget and r['ledger']['units']<=budget+1e-8 for r in rows)
    assert all(r['r_free']==0 and r['L_det'] is None for r in rows)
    assert all(np.isfinite(r['metrics']['pose_error_lever_m']) for r in rows)
    tables=[]
    for m in methods:
        rr=[r for r in rows if r['method']==m]
        tables.append({'method':m,'n':len(rr),'success_n':sum(r['success'] for r in rr),'success_rate':np.mean([r['success'] for r in rr]),'pose_m':np.mean([r['metrics']['pose_error_lever_m'] for r in rr]),'task_map_rmse':np.mean([r['metrics']['task_map_rmse'] for r in rr]),'full_coeff_rmse':np.mean([r['metrics']['full_map_rmse'] for r in rr]),'mean_work':np.mean([r['ledger']['units'] for r in rr]),'mean_wall_s':np.mean([r['wall_seconds'] for r in rr]),'completed_budget_policy':sum(r['finished_budget_policy'] for r in rr),'reduced_moves':sum(r['reduced_moves'] or 0 for r in rr),'certificate_probes':sum(len(r['cert_history']) for r in rr),'fallback_events':sum(len(r['fallback_events']) for r in rr),'final_data_evaluations':sum((r['stage_evaluations'] or [0,0,0,0])[-1] for r in rr),'local_fullrank_but_recovery_failed':sum(r['false_calibrated'] for r in rr)})
    dif=[]
    for s in seeds:
        ds=[]
        for i in range(12):
            p,d=by[(s,i,'prasc')],by[(s,i,'direct')];q=ref[p['aperture']]
            ds.append([int(p['success'])-int(d['success']),(p['metrics']['pose_error_lever_m']-d['metrics']['pose_error_lever_m'])/q['Delta_x_noninf_m'],(p['metrics']['task_map_rmse']-d['metrics']['task_map_rmse'])/q['Delta_chi_noninf']])
        dif.append(np.mean(ds,axis=0))
    dif=np.array(dif);rng=np.random.default_rng(20260906);idx=rng.integers(0,20,(10000,20));boot=dif[idx].mean(axis=1)
    lb=float(np.quantile(boot[:,0],.05/3));up=np.quantile(boot[:,1:],1-.05/3,axis=0)
    gates=[lb>0,up[0]<1,up[1]<1]
    sec=[]
    for m in ['fixed_rank','phaseless','direct_control']:
        d=np.array([np.mean([int(by[(s,i,m)]['success'])-int(by[(s,i,'direct')]['success']) for i in range(12)]) for s in seeds]);b=d[idx].mean(axis=1)
        # Centered bootstrap null, exploratory secondary test only.
        pv=(1+np.count_nonzero((b-d.mean())>=d.mean()))/(len(b)+1)
        sec.append({'method_minus_direct':m,'success_difference':float(d.mean()),'ci95':np.quantile(b,[.025,.975]).tolist(),'p_one_sided_centered_bootstrap':float(pv)})
    order=sorted(range(3),key=lambda i:sec[i]['p_one_sided_centered_bootstrap']);last=0
    for j,k in enumerate(order):
        last=max(last,min(1,(3-j)*sec[k]['p_one_sided_centered_bootstrap']));sec[k]['holm_p']=last
    out={'budget':budget,'raw_sha256':hashlib.sha256((SOURCE/f'final_{budget}.jsonl').read_bytes()).hexdigest(),'validated_records':len(rows),'table':tables,'primary_comparison':{'mean_seed_differences':dif.mean(axis=0).tolist(),'one_sided_level':1-.05/3,'success_lower':lb,'pose_normalized_upper':float(up[0]),'map_normalized_upper':float(up[1]),'gates':[bool(v) for v in gates],'compound_pass':bool(all(gates)),'ci95_descriptive':np.quantile(boot,[.025,.975],axis=0).tolist(),'bootstrap_interval_sensitivity':{str(n):np.quantile(boot[:n],[.05/3,1-.05/3],axis=0).tolist() for n in [1000,5000,10000]},'all_zero_success_warning':bool(np.all(dif[:,0]==0))},'secondary_holm':sec,'estimator_mutations_after_freeze':False,'notes':'Mean lever-metric pose error and task RMSE are unconditional; primary budget200 only. A fullrank local Gram is not a covariance coverage gate.'}
    (DEST/f'e4_{budget}_summary.json').write_text(json.dumps(out,indent=2)+'\n')
    return out,rows
def plots(outputs):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(12,3.7));colors={'direct':'#324d69','prasc':'#0c8e82','fixed_rank':'#bb8553','phaseless':'#ac4f79','direct_control':'#888888'}
    for m in colors:
        for ax,key in zip(axs,['success_rate','pose_m','task_map_rmse']):
            vals=[next(t for t in o[0]['table'] if t['method']==m)[key] for o in outputs]
            ax.plot([200,800],vals,'o-',label=m,color=colors[m])
    for ax,title in zip(axs,['Joint recovery (240 runs/method)','Mean pose error (m)','Mean resolved-map RMSE']):ax.set(title=title,xlabel='Charged work cap');ax.set_xticks([200,800]);ax.grid(alpha=.2)
    axs[0].set_ylim(0,1)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=5,frameon=False);fig.tight_layout(rect=[0,.13,1,1]);fig.savefig(FIG/'e4_budget_comparison.png',dpi=180);plt.close(fig)
    # Basin view: values are empirical success over five scenes per stratum,
    # never interpolated between the twelve prespecified initializations.
    fig,axs=plt.subplots(2,2,figsize=(9,6.5))
    for bi,(out,rs) in enumerate(outputs):
        for mi,m in enumerate(['direct','prasc']):
            arr=np.array([[np.mean([r['success'] for r in rs if r['method']==m and r['init_index']==i+4*j]) for i in range(4)] for j in range(3)])
            ax=axs[bi,mi];ax.imshow(arr,vmin=0,vmax=1,cmap='viridis',aspect='auto');ax.set(title=f'{m}; budget {out["budget"]}',xticks=range(4),xticklabels=['0','90','180','270'],yticks=range(3),yticklabels=['0.125','0.5','1']);ax.set_xlabel('Translation direction (deg)');ax.set_ylabel('Initial metric error / wavelength')
            for j in range(3):
                for i in range(4):ax.text(i,j,f'{arr[j,i]:.2f}',ha='center',va='center',color='white')
    fig.tight_layout();fig.savefig(FIG/'e4_basin_map.png',dpi=180);plt.close(fig)
def main():
    DEST.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True);oo=[analyze(b) for b in [200,800]];plots(oo)
    lines=['# Frozen E4 results','', 'All 2,400 records verified: unique complete pairing, budgets, finite errors, r_free=0, no L_det substitution.','']
    for o,rs in oo:
        lines += [f'## Budget {o["budget"]}', '', '| Method | successes/240 | pose m | task map RMSE | work | wall s | reduced moves |','|---|---:|---:|---:|---:|---:|---:|']
        for t in o['table']:lines.append(f'| {t["method"]} | {t["success_n"]} | {t["pose_m"]:.5f} | {t["task_map_rmse"]:.5f} | {t["mean_work"]:.2f} | {t["mean_wall_s"]:.3f} | {t["reduced_moves"]} |')
        p=o['primary_comparison'];lines += ['',f'Simultaneous limits: success lower {p["success_lower"]:.6f}; normalized pose upper {p["pose_normalized_upper"]:.6f}; normalized map upper {p["map_normalized_upper"]:.6f}. Success gate: FAIL. Both noninferiority checks: VACUOUS (identical outputs), not scientific passes. Compound pass: {p["compound_pass"]}.','']
    (DEST/'E4_RESULTS.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
if __name__=='__main__':main()
