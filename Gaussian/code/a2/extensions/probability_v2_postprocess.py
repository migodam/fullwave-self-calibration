"""Non-solver postprocessing for the frozen probability_v2 execution."""
from pathlib import Path
import json, hashlib, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'code'))
OUT=ROOT/'runs/a2/extensions/probability/v2'
FIG=ROOT/'figures/a2/extensions'

def mean(rows,key): return float(np.mean([r[key] for r in rows]))
def main():
    result=json.loads((OUT/'results.json').read_text())
    pilot=json.loads((ROOT/'runs/a2/extensions/probability/results.json').read_text())
    # Exact enumeration has no optimizer.  Zero seconds in the solver record did
    # not mean it was free, so make that field explicitly unavailable.
    for row in result['records']:
        row['methods']['exact']['optimizer_seconds']=None
        row['methods']['exact']['optimizer_timing_note']='No separate timing: finite enumeration posterior arithmetic shares the cached-field likelihood evaluation.'
    methods=['exact','product_vi','gaussian_rbf_logit_vi']; groups=['high','medium','low']
    aggregate={}
    for method in methods:
        aggregate[method]={}
        for group in groups:
            rows=[r['methods'][method] for r in result['records'] if r['group']==group]
            aggregate[method][group]={
                'n':len(rows),'fit_field_relative_to_noisy_y':mean(rows,'fit_field_relative'),
                'held_field_relative_to_noisy_y':mean(rows,'held_field_relative'),
                'kl_to_exact':mean(rows,'kl_to_exact'),'brier':mean(rows,'brier'),
                'optimizer_failure_count':sum(x['optimizer_success'] is False for x in rows),
                'optimizer_mean_seconds':None if method=='exact' else mean(rows,'optimizer_seconds')}
    result['n32_cache_construction_cost_from_pilot']=pilot['cache_cost']
    result['postprocess_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['field_error_note']='Fit and held errors are against independently noise-corrupted N64 data channels; they are not clean-manifold error and should not be compared as such.'
    (OUT/'results.json').write_text(json.dumps(result,indent=2))
    (OUT/'metrics.json').write_text(json.dumps({'aggregate':aggregate,'max_elbo_kl_identity_error':result['max_elbo_kl_identity_error'],'timing_note':'N32 construction is inherited from the pilot; N64 forward cost and end-to-end v2 wall time are recorded in results.json.'},indent=2))
    # Four columns use the same material scale. Coordinates are physical cm.
    from probability_run import centres,chi_for,patch_map,geom,CFG as PILOT_CFG
    from a2.physics import VIE
    examples={}
    for r in result['records']: examples.setdefault(r['group'],r)
    v=VIE(geom(32),1.5e9); extent=(-20,20,-20,20); vmax=.55
    fig,ax=plt.subplots(3,4,figsize=(11,8),sharex=True,sharey=True,constrained_layout=True)
    im=None
    for i,group in enumerate(groups):
        row=examples[group];maps=[chi_for(v.points,row['truth_state']).real]
        ids,inside=patch_map(v.points)
        for name in ('exact','product_vi','gaussian_rbf_logit_vi'):
            q=np.asarray(row['methods'][name]['posterior_patch_mean'])
            maps.append(np.where(inside,q[ids]*.55,0).reshape(32,32))
        for j,(data,title) in enumerate(zip(maps,['truth','exact','product VI','Gaussian-RBF VI'])):
            im=ax[i,j].imshow(data.reshape(32,32).T if data.ndim==1 else data.T,origin='lower',extent=extent,vmin=0,vmax=vmax,cmap='viridis')
            if i==0: ax[i,j].set_title(title)
            if j==0: ax[i,j].set_ylabel(f'{group}, state {row["truth_index"]}\ny (cm)')
            if i==2: ax[i,j].set_xlabel('x (cm)')
    fig.colorbar(im,ax=ax.ravel().tolist(),label='material contrast / posterior mean; common scale')
    fig.savefig(FIG/'probability_v2_reconstructions.png',dpi=170);plt.close(fig)
    # This is prior-predictive calibration: each point pools separately drawn
    # states and noise, not repeated noise for a fixed material realization.
    fig,ax=plt.subplots(figsize=(5,4)); bins=np.linspace(0,1,6)
    for method in methods:
        pairs=[]
        for r in result['records']:
            p=r['methods'][method]['posterior_patch_mean'];pairs.extend(zip(p,r['truth_state']))
        p,z=np.asarray(pairs).T;xs=[];ys=[]
        for lo,hi in zip(bins[:-1],bins[1:]):
            keep=(p>=lo)&(p<(hi if hi<1 else hi+1e-12))
            if keep.any(): xs.append(float(p[keep].mean()));ys.append(float(z[keep].mean()))
        ax.plot(xs,ys,'o-',label=method)
    ax.plot([0,1],[0,1],'k--',lw=1);ax.set(xlabel='posterior probability bin mean',ylabel='prior-predictive empirical frequency');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(FIG/'probability_v2_reliability.png',dpi=170)
    print(json.dumps(aggregate,indent=2))
if __name__=='__main__': main()
