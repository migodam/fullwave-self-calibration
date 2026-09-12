"""Read-only extension-result summary. Refuses partial manifold_v3 campaigns."""
from pathlib import Path
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
RUN=ROOT/'runs/a2/extensions'; FIG=ROOT/'figures/a2/extensions'; DOC=ROOT/'delegated/a2_extensions'
METHODS=['full_lm','spectral','weak_secant','random_secant']; BUDGETS=[1000,2000,3000,3400]
def load(p): return json.loads(Path(p).read_text())
def mean(xs): return None if not xs else float(np.mean(xs))
def at_budget(meta,b):
    h=[z for z in meta['history'] if z['rhs']<=b]
    return None if not h else h[-1]
def probability_summary():
    v2=load(RUN/'probability/v2/results.json'); weak=load(RUN/'probability/weak/results.json'); opt=load(RUN/'probability/v2/optimization_posthoc/results.json'); cap=load(RUN/'rbf_capacity.json')
    out={'frozen_v2_development':{'records':len(v2['records']),'max_elbo_kl_identity_error':v2['max_elbo_kl_identity_error'],'methods':{}},'posthoc_weak_data_diagnostic':{'records':len(weak['records']),'max_elbo_kl_identity_error':weak['max_elbo_kl_identity_error'],'methods':{}},'rbf_capacity_posthoc':{'feature_rank':cap['feature_rank'],'separable_count':cap['separable_count'],'total':cap['total'],'scope':cap['scope']},'posthoc_initialization_diagnostic':{}}
    for name in ['exact','product_vi','gaussian_rbf_logit_vi']:
        for src,key in [(v2,'frozen_v2_development'),(weak,'posthoc_weak_data_diagnostic')]:
            d={}
            for group in ['high','medium','low']:
                x=[r['methods'][name] for r in src['records'] if r['group']==group]
                d[group]={'n':len(x),'kl_to_exact_mean':mean([z['kl_to_exact'] for z in x]),'brier_mean':mean([z['brier'] for z in x]),'entropy_mean':mean([z['entropy'] for z in x]),'optimizer_failure_count':sum(z['optimizer_success'] is False for z in x),'held_error_to_noisy_y_mean':mean([z.get('held_field_relative',z.get('field_relative_to_noisy_y')) for z in x])}
            out[key]['methods'][name]=d
    status={tuple(x['state']):x['separable'] for x in cap['states']};base={(r['group'],r['rep']):r['methods']['gaussian_rbf_logit_vi']['kl_to_exact'] for r in v2['records']}
    rows=[]
    for r in opt['records']:
        b=base[(r['group'],r['rep'])];s=r['selected'];rows.append({'capacity_separable':bool(status[tuple(load(RUN/'probability/v2/results.json')['records'][['high','medium','low'].index(r['group'])*10+r['rep']]['truth_state'])]),'baseline_kl':b,'selected_kl':s['kl_to_exact'],'improved':s['kl_to_exact']<b,'selected_start':s['start']})
    for flag in [True,False]:
        x=[r for r in rows if r['capacity_separable']==flag];out['posthoc_initialization_diagnostic'][str(flag)]={'n':len(x),'improved_count':sum(r['improved'] for r in x),'baseline_kl_mean':mean([r['baseline_kl'] for r in x]),'selected_kl_mean':mean([r['selected_kl'] for r in x]),'selected_start_counts':{s:sum(r['selected_start']==s for r in x) for s in ['zero','product_projected','random_0','random_1']}}
    out['interpretation_guard']='An optimizer success flag only means the declared numerical optimizer terminated successfully; it does not establish posterior correctness. Weak and initialization results are post-hoc diagnostics and do not replace frozen v2.'
    return out
def manifold_summary(data):
    records=data['records'];out={'config':data['config'],'records':len(records),'budgets':{},'terminal':{},'paired_note':'Each family/budget cell retains one value per shared 10-case seed for each method. Missing means no recorded history point at or below that budget; no interpolation was used.'}
    for fam in data['config']['families']:
        rs=[r for r in records if r['family']==fam];out['budgets'][fam]={};out['terminal'][fam]={}
        for b in BUDGETS:
            cell={}
            for m in METHODS:
                vals=[]
                for r in rs:
                    h=at_budget(r['methods'][m],b)
                    vals.append(None if h is None else {'material_relative_l2':h['material_relative_l2'],'held_clean_relative':h['heldout_clean_relative'],'rhs':h['rhs']})
                valid=[v for v in vals if v is not None];cell[m]={'n_reached':len(valid),'material_l2_mean':mean([v['material_relative_l2'] for v in valid]),'held_clean_mean':mean([v['held_clean_relative'] for v in valid]),'paired_cases':vals}
            # Paired weak-vs-full differences only where both histories reached.
            pairs=[]
            for i in range(len(rs)):
                a,bw=cell['full_lm']['paired_cases'][i],cell['weak_secant']['paired_cases'][i]
                if a is not None and bw is not None:pairs.append({'weak_minus_full_material_l2':bw['material_relative_l2']-a['material_relative_l2'],'weak_minus_full_held_clean':bw['held_clean_relative']-a['held_clean_relative']})
            cell['paired_weak_minus_full']={'n':len(pairs),'material_l2_mean':mean([x['weak_minus_full_material_l2'] for x in pairs]),'held_clean_mean':mean([x['weak_minus_full_held_clean'] for x in pairs]),'cases':pairs}
            out['budgets'][fam][str(b)]=cell
        for m in METHODS:
            x=[r['methods'][m] for r in rs];out['terminal'][fam][m]={'n':len(x),'material_l2_mean':mean([z['material_relative_l2'] for z in x]),'held_clean_mean':mean([z['heldout_clean_relative'] for z in x]),'seconds_mean':mean([z['seconds'] for z in x]),'actual_rhs_mean':mean([z['rhs'] for z in x]),'actual_rhs_values':[z['rhs'] for z in x]}
    out['terminal_cost_note']='Terminal actual RHS and seconds are reported separately because methods often terminate at different RHS; this is not described as equal-RHS comparison.'
    return out
def plot_budget(summary):
    fig,ax=plt.subplots(2,3,figsize=(12,6.5))
    fig.subplots_adjust(top=.88,bottom=.12,wspace=.08,hspace=.22)
    colors=dict(full_lm='C0',spectral='C1',weak_secant='C2',random_secant='C3')
    for col,fam in enumerate(summary['config']['families']):
        for m in METHODS:
            x=[];y=[]
            for b in BUDGETS:
                v=summary['budgets'][fam][str(b)][m]
                if v['material_l2_mean'] is not None:x.append(b);y.append(v['material_l2_mean'])
            if x:ax[0,col].plot(x,np.asarray(y)*100,'o-',color=colors[m],label=m)
            x=[];y=[]
            for b in BUDGETS:
                v=summary['budgets'][fam][str(b)][m]
                if v['held_clean_mean'] is not None:x.append(b);y.append(v['held_clean_mean'])
            if x:ax[1,col].plot(x,np.asarray(y)*100,'o-',color=colors[m],label=m)
        ax[0,col].set_title(fam)
        for row in range(2):
            ax[row,col].set(xlabel='RHS budget',xticks=BUDGETS)
            if col==2: ax[row,col].yaxis.tick_right()
    ax[0,0].set_ylabel('mean material relative L2 (%)');ax[1,0].set_ylabel('mean held clean-field relative error (%)')
    fig.legend(*ax[0,0].get_legend_handles_labels(),loc='upper center',ncol=4,fontsize=8,bbox_to_anchor=(.5,.985))
    fig.text(.5,.025,'Each panel has its own y scale. Markers are recorded history points only; no extra solves or interpolation.',ha='center',fontsize=8)
    fig.savefig(FIG/'manifold_budget.png',dpi=170);plt.close(fig)
def plot_best_worst(data):
    imgs=np.load(RUN/'manifold_v3/images.npz')['images'];recs=data['records'];families=data['config']['families'];cols=['truth']+METHODS
    picks=[]
    for fam in families:
        ind=[i for i,r in enumerate(recs) if r['family']==fam];best=min(ind,key=lambda i:recs[i]['methods']['weak_secant']['material_relative_l2']);worst=max(ind,key=lambda i:recs[i]['methods']['weak_secant']['material_relative_l2']);picks.extend([best,worst])
    vmax=float(max(imgs[i,0].max() for i in picks));fig,axes=plt.subplots(len(picks),5,figsize=(12,2.3*len(picks)),sharex=True,sharey=True,constrained_layout=True)
    for row,i in enumerate(picks):
        tag='best' if row%2==0 else 'worst';r=recs[i]
        for j,a in enumerate(axes[row]):
            im=a.imshow(imgs[i,j].reshape(32,32).T,origin='lower',extent=(-20,20,-20,20),vmin=0,vmax=vmax,cmap='viridis');a.set_title(cols[j] if row==0 else '');
            if j==0:a.set_ylabel(f"{r['family']} {tag}\ny (cm)")
            if row==len(picks)-1:a.set_xlabel('x (cm)')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='real contrast; common scale');fig.savefig(FIG/'manifold_best_worst.png',dpi=170);plt.close(fig)
def plot_spectral_failures(data):
    imgs=np.load(RUN/'manifold_v3/images.npz')['images'];recs=data['records'];picks=[]
    for fam in data['config']['families']:
        ids=[i for i,r in enumerate(recs) if r['family']==fam];picks.append(max(ids,key=lambda i:recs[i]['methods']['spectral']['material_relative_l2']-recs[i]['methods']['full_lm']['material_relative_l2']))
    vmax=float(imgs[picks].max());fig,ax=plt.subplots(3,5,figsize=(13,8),constrained_layout=True,sharex=True,sharey=True);metadata=[]
    for row,i in enumerate(picks):
        r=recs[i];metadata.append({'record_index':i,'family':r['family'],'selection':'largest terminal spectral minus full_lm material error; post-hoc display','methods':{m:{k:r['methods'][m][k] for k in ['material_relative_l2','heldout_clean_relative','rhs']} for m in METHODS}})
        for col,a in enumerate(ax[row]):
            im=a.imshow(imgs[i,col].reshape(32,32).T,origin='lower',extent=(-20,20,-20,20),vmin=0,vmax=vmax,cmap='viridis')
            if col==0:a.set_title(f"{r['family']} truth");a.set_ylabel('y (cm)')
            else:
                m=METHODS[col-1];v=r['methods'][m];a.set_title(f"{m}\nL2 {100*v['material_relative_l2']:.2f}%, RHS {v['rhs']}",fontsize=9)
            if row==2:a.set_xlabel('x (cm)')
    fig.colorbar(im,ax=ax.ravel().tolist(),label='real contrast; common scale');fig.savefig(FIG/'manifold_spectral_failures.png',dpi=170);plt.close(fig)
    (RUN/'manifold_spectral_failure_selection.json').write_text(json.dumps(metadata,indent=2))
def plot_probability(p):
    fig,ax=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True);methods=['exact','product_vi','gaussian_rbf_logit_vi'];groups=['high','medium','low'];labels={'frozen_v2_development':['0.03','0.1','0.3'],'posthoc_weak_data_diagnostic':['0.3','0.8','1.5']}
    for panel,key,title in [(0,'frozen_v2_development','frozen v2 development'),(1,'posthoc_weak_data_diagnostic','weak-data post-hoc diagnostic')]:
        for m in methods:ax[panel].plot(labels[key],[p[key]['methods'][m][g]['kl_to_exact_mean'] for g in groups],'o-',label=m)
        ax[panel].set(yscale='symlog',ylabel='mean KL to exact',xlabel='noise fraction',title=title);ax[panel].legend(fontsize=7)
    cap=p['rbf_capacity_posthoc'];ax[2].bar(['LP separable','LP inseparable'],[cap['separable_count'],cap['total']-cap['separable_count']],color=['C2','C3']);ax[2].set(ylabel='of 256 finite states',title='fixed RBF feature capacity (post-hoc)')
    fig.savefig(FIG/'probability_summary.png',dpi=170);plt.close(fig)
def write_doc(p,m):
    lines=['# A2 extensions 汇总（自动读取结果）','',f"manifold_v3 已读取 {m['records']}/30 cases；预算曲线只取 history 中 `rhs <= budget` 的最后点，缺失不插值。",'',f"概率 frozen v2 为 {p['frozen_v2_development']['records']} 条；weak/posthoc 为 {p['posthoc_weak_data_diagnostic']['records']} 条。固定 RBF LP 可分 {p['rbf_capacity_posthoc']['separable_count']}/{p['rbf_capacity_posthoc']['total']}。",'',p['interpretation_guard'],'',m['terminal_cost_note'],'','详细机器可读数据：`runs/a2/extensions/SUMMARY.json`。图：`figures/a2/extensions/probability_summary.png`、`manifold_budget.png`、`manifold_best_worst.png`。']
    (DOC/'SUMMARY_ZH.md').write_text('\n'.join(lines)+'\n')
def main():
    data=load(RUN/'manifold_v3/results.json')
    if len(data.get('records',[]))!=30: raise SystemExit(f'manifold_v3 incomplete: {len(data.get("records",[]))}/30; summary not generated')
    p=probability_summary();m=manifold_summary(data);out={'probability':p,'manifold_v3':m};atomic=RUN/'SUMMARY.json';atomic.write_text(json.dumps(out,indent=2));plot_probability(p);plot_budget(m);plot_best_worst(data);plot_spectral_failures(data);DOC.mkdir(parents=True,exist_ok=True);write_doc(p,m);print(json.dumps({'written':str(atomic),'records':m['records']},indent=2))
if __name__=='__main__':main()
