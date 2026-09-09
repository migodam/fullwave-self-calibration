"""Figures and compact evidence from completed independent-reference audits."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})


def main():
    summary={}
    fig,axs=plt.subplots(1,2,figsize=(9.6,3.9))
    for ax,shape in zip(axs,['ellipsoid','box']):
        rows=json.loads((OUT/f'results/nonspherical_{shape}.json').read_text())
        summary[shape]=rows
        for k in [9.,18.]:
            vals=sorted([r for r in rows if r['k']==k],key=lambda r:r['n'])
            ax.plot([r['spacing']*1000 for r in vals],[r['relative_field_disagreement']*100 for r in vals],'-o',label=f'k={k:g}')
        ax.set(title=shape.capitalize(),xlabel='Cell spacing (mm)',ylabel='ADDA / inverse field difference (%)',yscale='log')
        ax.invert_xaxis();ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Independent implementations, shared DDA approximation',fontsize=11)
    fig.tight_layout();fig.savefig(OUT/'figures/nonspherical_refinement.png');plt.close(fig)

    rows=json.loads((OUT/'results/prediction_audit3d.json').read_text())
    fig,ax=plt.subplots(figsize=(7.2,4))
    groups=[]
    for method,label in [('unified','Without reference'),('unified_reference','Noisy electronics reference')]:
        rr=[r for r in rows if r['spacing']==.01 and max(r['frequencies'])==9 and r['method']==method]
        groups.append([100*np.mean([r[key]['field_relative_error'] for r in rr]) for key in ['deployable_sensor_prediction','reconstructed_structural_field']])
    xx=np.arange(2)
    ax.bar(xx-.18,groups[0],.36,label='Without reference',color='#345e89')
    ax.bar(xx+.18,groups[1],.36,label='Noisy electronics reference',color='#bc7738')
    ax.set(xticks=xx,xticklabels=['Actual sensor prediction','Reconstructed structural field'],ylabel='Relative field error (%)',
           title='Good prediction does not ensure good recovered physics')
    ax.legend();ax.grid(axis='y',alpha=.2);fig.tight_layout()
    fig.savefig(OUT/'figures/prediction_vs_recovery.png');plt.close(fig)
    summary['low_band_prediction_field_error_percent']=groups

    path=OUT/'results/nonspherical_calibration.json'
    if path.exists():
        rows=json.loads(path.read_text());aggregate=[]
        for n in sorted({r['n'] for r in rows}):
            for method in ['fixed_wrong_pose','geometry_only','unified','unified_reference']:
                rr=[r for r in rows if r['n']==n and r['method']==method]
                if not rr:continue
                aggregate.append(dict(n=n,method=method,scenes=len(rr),
                    pose_mm=float(1000*np.mean([r['pose_error_m'] for r in rr])),
                    material_percent=float(100*np.mean([r['material_relative_error'] for r in rr])),
                    sensor_field_percent=float(100*np.mean([r['deployable_prediction']['relative_field_error'] for r in rr])),
                    sensor_phase_rad=float(np.mean([r['deployable_prediction']['weighted_phase_rmse_rad'] for r in rr])),
                    structural_field_percent=float(100*np.mean([r['reconstructed_structural_field']['relative_field_error'] for r in rr])),
                    converged=sum(r['optimizer_status']>0 for r in rr)))
        summary['ellipsoid_calibration']=aggregate
        fig,axs=plt.subplots(1,2,figsize=(9,4.2));vals=[r for r in aggregate if r['n']==32]
        labels=['Wrong pose','Geometry only','Joint','Joint + reference']
        for ax,key,ylabel in zip(axs,['pose_mm','material_percent'],['Position error (mm)','Material error (%)']):
            ax.bar(np.arange(len(vals)),[r[key] for r in vals],color=['#969da3','#c68e60','#335f89','#318d82'])
            ax.set(xticks=np.arange(len(vals)),xticklabels=labels,ylabel=ylabel,yscale='log')
            ax.tick_params(axis='x',rotation=20);ax.grid(axis='y',alpha=.2)
        fig.suptitle('Known-support ellipsoid: 2 development scenes, not general imaging',fontsize=11)
        fig.tight_layout();fig.savefig(OUT/'figures/ellipsoid_calibration.png');plt.close(fig)
    (OUT/'results/independent_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary.get('ellipsoid_calibration',[]),indent=2))


if __name__=='__main__':main()
