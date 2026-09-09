"""Generate evidence tables/figures from completed raw records only."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
DELEGATED=OUT.parents[1]/'delegated'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})

def main():
    (OUT/'figures').mkdir(exist_ok=True)
    aggregate={}
    for fname in ['reuse_n16_q9_r128.json','reuse_n32_q49_r128.json']:
        path=DELEGATED/'a3_rom_reuse'/fname
        if not path.exists():continue
        rows=json.loads(path.read_text());summ=[]
        for method in sorted({r['method'] for r in rows}):
            rr=[r for r in rows if r['method']==method]
            item=dict(method=method,scenes=len(rr))
            for k in ['pose_error','material_rmse','phase_rms','final_loss','online_seconds','cold_seconds','audit_seconds','reduced_steps','reuse_hits','fallback_steps','rebuilds']:
                item[k]=float(np.mean([r[k] for r in rr]))
            item['full_rhs_mean']=float(np.mean([r['online_work']['full_rhs_columns'] for r in rr]))
            errors=[e['offline_high_frequency'] for r in rr for e in r['events'] if 'offline_high_frequency' in e]
            item['high_frequency_audited_updates']=len(errors)
            item['high_frequency_tangent_passes']=sum(e['map_jac_relative_error']<=.01 and e['pose_jac_relative_error']<=.01 for e in errors)
            item['max_map_jac_error']=max([e['map_jac_relative_error'] for e in errors],default=None)
            item['max_pose_jac_error']=max([e['pose_jac_relative_error'] for e in errors],default=None)
            summ.append(item)
        aggregate[fname]=summ
    path=DELEGATED/'a3_maxwell_refine'/'calibration_fft_raw.json'
    if path.exists():
        rows=json.loads(path.read_text())
        aggregate['calibration3d']=[{k:r.get(k) for k in ['seed','spacing','frequencies','method','pose_error_m','material_relative_error','heldout_phase_rmse_rad','reduced_chisquare','nfev','optimizer_status','status','work_ledger_match','total_elapsed_seconds']} for r in rows]
    (OUT/'results'/'second_pass_summary.json').write_text(json.dumps(aggregate,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(6.5,4))
    old=json.loads((DELEGATED/'a3_maxwell_fft'/'treams_audit_raw.json').read_text())
    if isinstance(old,dict):old=old.get('rows',old.get('cases',[]))
    fine=json.loads((OUT/'results'/'maxwell3d_fine_forward.json').read_text())
    for k in [9.,18.]:
        vals={(r['spacing'],r['relative_field_error_high']) for r in old+fine if r.get('scene')=='pair' and r['k']==k and r.get('status')=='ok'}
        vals=sorted(vals,reverse=True)
        ax.plot([r[0]*1000 for r in vals],[100*r[1] for r in vals],'-o',label=f'k={k:g} rad/m')
    ax.set(xlabel='Inverse-model cell spacing (mm)',ylabel='Relative vector-field error (%)',title='Independent sphere-cluster reference: refinement matters')
    ax.invert_xaxis();ax.legend();ax.grid(alpha=.2);fig.tight_layout()
    fig.savefig(OUT/'figures'/'maxwell_refinement.png');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(9,4))
    colors=['#354c7c','#37828a','#cb8640','#8c5c93']
    for ax,(fname,title) in zip(axs,[('reuse_n16_q9_r128.json','9 material parameters / 4 scenes'),('reuse_n32_q49_r128.json','49 material parameters / up to 2 scenes')]):
        values=aggregate.get(fname,[])
        if not values:continue
        names=[v['method'].replace('direct_','').replace('_task','') for v in values]
        ax.bar(names,[v['online_seconds'] for v in values],color=colors)
        ax.set(title=title,ylabel='Development online wall time (s)')
        ax.tick_params(axis='x',rotation=20)
    fig.suptitle('Actual chart reuse; these are NOT isolated final timings',fontsize=11)
    fig.tight_layout();fig.savefig(OUT/'figures'/'rom_reuse_development.png');plt.close(fig)
    print(json.dumps(aggregate,indent=2))

if __name__=='__main__':main()
