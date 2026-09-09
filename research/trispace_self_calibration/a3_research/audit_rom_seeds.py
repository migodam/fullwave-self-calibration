"""Audit serial seed development comparisons and make a descriptive figure."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
METHODS=['direct_gn','direct_adjoint','generic_task','sensing_task','twofold_task','krylov_task']


def audit():
    allgroups=[]
    for q,count in [(9,4),(49,2)]:
        rows=json.loads((OUT/f'results/rom_seed_q{q}.json').read_text())
        assert len(rows)==6*count and len({(r['seed'],r['method']) for r in rows})==len(rows)
        for method in METHODS:
            group=[r for r in rows if r['method']==method]
            assert len(group)==count
            audits=[e['offline_high_frequency'] for r in group for e in r['events']
                    if e.get('accepted_reduced')]
            largest=max([max(a['map_jac_relative_error'],a['pose_jac_relative_error']) for a in audits],default=0.)
            allgroups.append(dict(q=q,method=method,scenes=count,
                mean_cold_seconds=float(np.mean([r['cold_seconds'] for r in group])),
                cold_seconds=[r['cold_seconds'] for r in group],
                mean_pose_error_m=float(np.mean([r['pose_error'] for r in group])),
                mean_material_rmse=float(np.mean([r['material_rmse'] for r in group])),
                mean_loss=float(np.mean([r['final_loss'] for r in group])),
                mean_full_rhs=float(np.mean([r['online_work']['full_rhs_columns'] for r in group])),
                accepted_reduced=sum(r['reduced_steps'] for r in group),
                fallback_steps=sum(r['fallback_steps'] for r in group),
                max_audited_jacobian_error=largest,
                all_accepted_jacobians_within_one_percent=largest<=.01,
                statuses=[r['status'] for r in group]))
        assert all(r['online_work']['full_rhs_columns']<=8000 for r in rows)
    result=dict(groups=allgroups,
        scope='serial DEVELOPMENT; four q9 and two q49 scenes; RHS/iteration caps are not equal wall budgets',
        completion_check=True)
    (OUT/'results/rom_seed_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    labels=['Direct\nGN','Direct\nadjoint','Generic\ntask','Sensing\ntask','Twofold\ntask','Krylov\ntask']
    for ax,q in zip(axes,[9,49]):
        groups=[g for g in allgroups if g['q']==q]
        ax.bar(range(6),[g['mean_cold_seconds'] for g in groups],color=['#596b7d']*2+['#2d8795']*4)
        for x,g in enumerate(groups):
            jitter=np.linspace(-.12,.12,len(g['cold_seconds']))
            ax.scatter(x+jitter,g['cold_seconds'],color='#131c24',s=18,zorder=3)
        ax.set_xticks(range(6),labels,fontsize=8)
        ax.set_ylabel('Total online + setup time (s)')
        ax.set_title(f'{q} material coefficients; {groups[0]["scenes"]} development scenes')
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Matched chart seeds: serial development costs, not final superiority')
    fig.text(.5,.02,'Bars: means; dots: individual scenes. Common RHS/iteration ceilings are not equal wall budgets.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.055,1,.95))
    fig.savefig(OUT/'figures/rom_seed_serial.png',dpi=180)
    print(json.dumps(result,indent=2))


if __name__=='__main__':audit()
