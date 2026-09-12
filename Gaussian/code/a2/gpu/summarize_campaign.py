"""Strict read-only postprocessing for a completed RTX representation campaign.

It intentionally refuses partial synchronisations: output only follows a full,
hash-consistent 30-case, four-method, 120-step main128 result set.
"""
from pathlib import Path
import json, hashlib, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
RUN=ROOT/'runs/a2/gpu/representation_campaign/main128';FIG=ROOT/'figures/a2/gpu';DOC=ROOT/'delegated/a2_gpu'
METHODS=['gaussian_K16','gaussian_K64','gaussian_K144','voxel'];P={'gaussian_K16':96,'gaussian_K64':384,'gaussian_K144':864,'voxel':16384}
def read(p): return json.loads(Path(p).read_text())
def atomic_json(p,x):
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,indent=2));q.replace(p)
def validate():
    f=read(RUN/'frozen_config.json');cfg=f['config'];h=f['config_hash']
    if cfg.get('steps')!=120 or cfg.get('inverse_n')!=128: raise SystemExit('main128 frozen config is not the declared 120-step N128 campaign')
    records=[];problems=[]
    for name,expected in f['source_hashes'].items():
        path=ROOT/name.replace(chr(92),'/')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:problems.append('source hash mismatch '+name)
    for case in range(30):
        path=RUN/f'case_{case:02d}_summary.json'
        if not path.exists():problems.append(f'missing {path.name}');continue
        r=read(path); got={x.get('method'):x for x in r.get('methods',[])}
        if set(got)!=set(METHODS):problems.append(f'case {case}: methods {sorted(got)}');continue
        for m in METHODS:
            x=got[m]
            if x.get('config_hash')!=h:problems.append(f'case {case} {m}: config hash mismatch')
            if x.get('steps_completed')!=120:problems.append(f'case {case} {m}: steps={x.get("steps_completed")}')
            if x.get('parameter_count')!=P[m]:problems.append(f'case {case} {m}: p={x.get("parameter_count")}')
            if x['counts_delta']['forward_rhs']!=4356 or x['counts_delta']['adjoint_rhs']!=4320:problems.append(f'case {case} {m}: wrong RHS accounting')
            if len(x['trajectory'])!=120 or [q['step'] for q in x['trajectory']]!=list(range(1,121)):problems.append(f'case {case} {m}: incomplete trajectory')
            for key in ['material_l2_relative','held_clean_field_relative','threshold_iou_0p12','actual_optimization_seconds']:
                if not np.isfinite(x[key]):problems.append(f'case {case} {m}: nonfinite '+key)
            a=RUN/f'case_{case:02d}_{m}_arrays.npz'
            if not a.exists():problems.append(f'missing {a.name}');continue
            z=np.load(a)
            if not {'truth','reconstruction'}<=set(z.files):problems.append(f'{a.name}: missing truth/reconstruction')
            elif not (np.isfinite(z['truth']).all() and np.isfinite(z['reconstruction']).all()):problems.append(f'{a.name}: nonfinite image')
            elif z['truth'].size!=128**2 or z['reconstruction'].size!=128**2:problems.append(f'{a.name}: wrong N128 array size')
            else:
                computed=np.linalg.norm(z['truth']-z['reconstruction'])/np.linalg.norm(z['truth'])
                if abs(computed-x['material_l2_relative'])>1e-6:problems.append(f'{a.name}: metric does not match image')
        dp=RUN/f'case_{case:02d}_data.npz'
        if not dp.exists():problems.append(f'missing {dp.name}')
        records.append(r)
    if problems: raise SystemExit('campaign incomplete/invalid; main summary not generated:\n'+'\n'.join(problems[:20]))
    return cfg,h,records
def avg(x): return float(np.mean(x))
def data(records):
    rows=[]
    for r in records:
        for x in r['methods']:
            y=dict(x);y['family']=r['family'];y['data_seconds_including_first_setup']=r['data_generation'].get('seconds_including_first_setup');y['inverse_setup_seconds']=r.get('inverse_setup_seconds');y['method_seconds_including_checkpoint_io']=x['actual_optimization_seconds']
            rows.append(y)
    return rows
def summarize(cfg,h,records):
    rows=data(records);out={'config':cfg,'config_hash':h,'scope':'Bounded noiseless development representation comparison. Ordinary Gaussian components versus voxels; no SOM contribution. Training trajectories are fit loss only; no held intermediate frontier was stored.','methods':{},'families':{},'paired_k144_vs_voxel':{}}
    for m in METHODS:
        x=[r for r in rows if r['method']==m];out['methods'][m]={'parameter_count':P[m],'n':len(x),'material_l2_mean':avg([r['material_l2_relative'] for r in x]),'held_clean_mean':avg([r['held_clean_field_relative'] for r in x]),'iou_mean':avg([r['threshold_iou_0p12'] for r in x]),'seconds_including_checkpoint_io_mean':avg([r['method_seconds_including_checkpoint_io'] for r in x]),'counts_mean':{k:avg([r['counts_delta'][k] for r in x]) for k in x[0]['counts_delta']}}
    for fam in [records[i]['family'] for i in (0,10,20)]:
        rr=[r for r in records if r['family']==fam];out['families'][fam]={}
        for m in METHODS:
            x=[next(q for q in r['methods'] if q['method']==m) for r in rr];out['families'][fam][m]={'n':len(x),'material_l2_mean':avg([q['material_l2_relative'] for q in x]),'held_clean_mean':avg([q['held_clean_field_relative'] for q in x]),'iou_mean':avg([q['threshold_iou_0p12'] for q in x]),'seconds_including_checkpoint_io_mean':avg([q['actual_optimization_seconds'] for q in x])}
        pair=[]
        for r in rr:
            a=next(q for q in r['methods'] if q['method']=='gaussian_K144');b=next(q for q in r['methods'] if q['method']=='voxel');pair.append({'case':r['case'],'k144_minus_voxel_material_l2':a['material_l2_relative']-b['material_l2_relative'],'k144_minus_voxel_held_clean':a['held_clean_field_relative']-b['held_clean_field_relative'],'k144_win_material_l2':a['material_l2_relative']<b['material_l2_relative']})
        out['paired_k144_vs_voxel'][fam]={'n':len(pair),'k144_material_win_count':sum(q['k144_win_material_l2'] for q in pair),'delta_material_l2_mean':avg([q['k144_minus_voxel_material_l2'] for q in pair]),'delta_held_clean_mean':avg([q['k144_minus_voxel_held_clean'] for q in pair]),'cases':pair}
    out['paired_cost_excluding_interrupted_case18']={}
    rr=[r for r in records if r['case']!=18]
    for m in METHODS:
        ratios=[];times=[]
        for r in rr:
            a=next(q for q in r['methods'] if q['method']==m);v=next(q for q in r['methods'] if q['method']=='voxel');times.append(a['actual_optimization_seconds']);ratios.append(a['actual_optimization_seconds']/v['actual_optimization_seconds'])
        out['paired_cost_excluding_interrupted_case18'][m]={'n':len(rr),'seconds_including_checkpoint_io_mean':avg(times),'paired_seconds_ratio_to_voxel_median':float(np.median(ratios))}
    out['cost_accounting']={'completed_optimization_seconds_including_checkpoint_io_sum':sum(r['actual_optimization_seconds'] for r in rows),'recorded_data_generation_seconds_sum':sum(r['data_generation']['seconds_including_first_setup'] for r in records),'recorded_inverse_setup_seconds_sum':sum(r['inverse_setup_seconds'] for r in records),'peak_cuda_allocated_bytes_across_invocations':max(r['cuda_peak_bytes'] for r in records),'not_fully_metered':['interrupted incomplete K144 work','original case18 data generation and inverse setup before interruption','plotting, final field evaluation, inter-stage Python overhead, sync and downtime'],'interpretation':'Recorded costs are component totals, not complete uninterrupted end-to-end runtime. Cost ratios omit mixed pre/post-interruption case18; method order is fixed and checkpoint I/O remains included.'}
    out['cost_note']='actual_optimization_seconds measures the runner loop containing atomic checkpoints, so it includes checkpoint I/O. Data-generation and inverse setup costs are separate and are not silently added to each method.'
    return out,rows
def array(case,m): return np.load(RUN/f'case_{case:02d}_{m}_arrays.npz')
def representative(records):
    picks=[0,10,20];fig,ax=plt.subplots(3,5,figsize=(13,8),sharex=True,sharey=True,constrained_layout=True);vals=[]
    for c in picks:
        vals.append(array(c,METHODS[0])['truth'])
        vals += [array(c,m)['reconstruction'] for m in METHODS]
    vmax=max(float(x.max()) for x in vals)
    for row,c in enumerate(picks):
        fam=records[c]['family'];ims=[array(c,METHODS[0])['truth']]+[array(c,m)['reconstruction'] for m in METHODS]
        for col,(a,z,title) in enumerate(zip(ax[row],ims,['truth']+METHODS)):
            im=a.imshow(z.reshape(128,128).T,origin='lower',extent=(-20,20,-20,20),vmin=0,vmax=vmax,cmap='viridis')
            if row==0:a.set_title(title)
            if col==0:a.set_ylabel(f'{fam}\ny (cm)')
            if row==2:a.set_xlabel('x (cm)')
    fig.colorbar(im,ax=ax.ravel().tolist(),label='real contrast; common scale');fig.savefig(FIG/'representation_main128_representative.png',dpi=170);plt.close(fig)
def best_worst(records,summary):
    choices=[]
    for fam,cell in summary['paired_k144_vs_voxel'].items():
        xs=cell['cases'];choices += [(fam,min(xs,key=lambda q:q['k144_minus_voxel_material_l2'])['case'],'best K144-vs-voxel'),(fam,max(xs,key=lambda q:q['k144_minus_voxel_material_l2'])['case'],'worst K144-vs-voxel')]
    diffs=[];images=[]
    short={'off_grid_gaussian_mixture':'smooth mixture','sharp_multi_inclusion':'sharp inclusions','curves_with_holes':'curves + holes'}
    for _,c,_ in choices:
        t=array(c,'voxel')['truth'];g=array(c,'gaussian_K144')['reconstruction'];v=array(c,'voxel')['reconstruction'];images.append([t,g,v]);diffs.append(np.abs(g-t)-np.abs(v-t))
    lim=max(float(np.max(np.abs(q))) for q in diffs);vmax=max(float(z.max()) for row in images for z in row)
    fig,ax=plt.subplots(6,4,figsize=(12,16),sharex=True,sharey=True,constrained_layout=True)
    for row,((fam,c,label),ims,d) in enumerate(zip(choices,images,diffs)):
        for col,z in enumerate(ims):
            im=ax[row,col].imshow(z.reshape(128,128).T,origin='lower',extent=(-20,20,-20,20),vmin=0,vmax=vmax,cmap='viridis')
        err=ax[row,3].imshow(d.reshape(128,128).T,origin='lower',extent=(-20,20,-20,20),vmin=-lim,vmax=lim,cmap='coolwarm')
        ax[row,0].set_ylabel(f'{short[fam]}\ncase {c}, '+('best' if row%2==0 else 'worst')+' difference\ny (cm)')
        if row==0:
            for a,title in zip(ax[row],['truth','Gaussian K144','voxel','absolute-error difference']):a.set_title(title)
        if row==5:
            for a in ax[row]:a.set_xlabel('x (cm)')
    fig.colorbar(im,ax=ax[:,:3].ravel().tolist(),label='real contrast; common scale',shrink=.55)
    fig.colorbar(err,ax=ax[:,3].ravel().tolist(),label='|K144-truth| - |voxel-truth|; blue favors K144',shrink=.55)
    fig.savefig(FIG/'representation_main128_k144_vs_voxel_best_worst.png',dpi=160);plt.close(fig)

def cost_accuracy(rows):
    fig,ax=plt.subplots(1,2,figsize=(9,3.7),constrained_layout=True);colors=dict(zip(METHODS,['C0','C1','C2','C3']))
    for m in METHODS:
        x=[r for r in rows if r['method']==m];s=[20+100*r['threshold_iou_0p12'] for r in x]
        ax[0].scatter([r['method_seconds_including_checkpoint_io'] for r in x],[r['material_l2_relative'] for r in x],s=s,color=colors[m],alpha=.7,label=m)
        ax[1].scatter([r['method_seconds_including_checkpoint_io'] for r in x],[r['held_clean_field_relative'] for r in x],s=s,color=colors[m],alpha=.7,label=m)
    ax[0].set(xlabel='method seconds incl. checkpoint I/O',ylabel='material relative L2');ax[1].set(xlabel='method seconds incl. checkpoint I/O',ylabel='held clean-field relative error');ax[1].legend(fontsize=7);fig.savefig(FIG/'representation_main128_cost_accuracy.png',dpi=170);plt.close(fig)
def scaling():
    p=ROOT/'runs/a2/gpu/scaling_remote.json'
    if not p.exists(): return None
    d=read(p);r=d.get('records',[])
    if not r:return None
    fig,ax=plt.subplots(1,2,figsize=(8,3.2),constrained_layout=True);n=[x['n'] for x in r];sec=[avg(x['forward_adjoint_seconds']) for x in r];mem=[x['peak_gpu_bytes']/2**20 for x in r]
    ax[0].plot(n,sec,'o-');ax[0].set(xlabel='grid n',ylabel='forward+adjoint seconds');ax[1].plot(n,mem,'o-');ax[1].set(xlabel='grid n',ylabel='peak GPU MiB');fig.savefig(FIG/'representation_gpu_scaling.png',dpi=170);plt.close(fig)
    return {'scope':d.get('scope'),'records':[{'n':x['n'],'state_dimension':x['state_dimension'],'forward_adjoint_seconds_mean':avg(x['forward_adjoint_seconds']),'peak_gpu_bytes':x['peak_gpu_bytes']} for x in r]}
def doc(s):
    lines=['# RTX 4060 Gaussian/voxel campaign 候选结果段落','',f"完整性检查通过：30 cases × 4 methods，N128 inverse、N192 cell-integrated data、120 steps、统一 config hash `{s['config_hash'][:12]}`。",'', '该有界 noiseless development comparison 将普通正 Gaussian component 表示（p=96/384/864）与 voxel（p=16384）比较；它不是 SOM 结果，也不证明正则化或步数设置已达到公平充分。每方法秒数来自包含原子 checkpoint I/O 的优化循环，数据生成与 setup 单列。', '', '报告 material L2、held clean-field error、IoU、RHS/迭代计数及 Gaussian K144 相对 voxel 的 10-case family-paired 差值；训练轨迹只是 fit loss，未保存中途 held 轨迹，故不作 held 前沿结论。全局峰值显存不作为方法间比较。', '', '机器可读汇总：`runs/a2/gpu/representation_campaign/main128/SUMMARY.json`。']
    DOC.mkdir(parents=True,exist_ok=True);(DOC/'SUMMARY_ZH.md').write_text('\n'.join(lines)+'\n')
def main():
    cfg,h,recs=validate();FIG.mkdir(parents=True,exist_ok=True);s,rows=summarize(cfg,h,recs);s['scaling_remote']=scaling();atomic_json(RUN/'SUMMARY.json',s);representative(recs);best_worst(recs,s);cost_accuracy(rows);doc(s);print(json.dumps({'records':len(recs),'methods':len(rows),'written':str(RUN/'SUMMARY.json')},indent=2))
if __name__=='__main__': main()
