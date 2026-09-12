"""Root synthesis from the complete 30-case record; no additional fitting."""
from pathlib import Path
import collections,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.som import FREQUENCIES_HZ,components
from a2.ports import render
from a2.run_imaging_campaign import material_truth
out=ROOT/'runs/a2/imaging';r=json.loads((out/'results.json').read_text());assert len(r)==30
families=['gaussian_close','gaussian_strong','sharp_boundary'];labels=['Close Gaussians','Higher contrast','Disk with notch'];ms=['ordinary_lm','fixed_twofold_tsomg','adaptive_tsomg'];names=['LM','固定twofold','自适应'];flabels=['相邻Gaussian','较高对比度','尖边缺口']
v=VIE(Geometry(n=64,n_tx=12,n_rx=48,aperture='half'),FREQUENCIES_HZ[0]);fig,axes=plt.subplots(3,4,figsize=(11,8.5),layout='constrained')
lines=['**完整执行：30个案例×5方法＝150次反演。** 下表均为每组十例的中位数；拒步停止与达到迭代上限分开计数，不称数值崩溃或统计失败率。',['|目标|方法|材料L2|留出场|秒/例|拒步/达cap（各/10）|','|---|---|---:|---:|---:|---:|']]
lines=[lines[0],'']+lines[1];checks={}
for i,(fam,label,fl) in enumerate(zip(families,labels,flabels)):
    rows=[q for q in r if q['family']==fam];assert len(rows)==10
    rows.sort(key=lambda z:z['methods']['adaptive_tsomg']['material_relative_l2_independent_data']);z=rows[5]
    arr=[material_truth(z,v)]+[render(components(np.array(z['methods'][m]['parameters']),FREQUENCIES_HZ[0]),v.points).real for m in ms];vmax=max(a.max() for a in arr)
    for j,a in enumerate(arr):
        im=axes[i,j].imshow(a.reshape(64,64).T,origin='lower',extent=(-20,20,-20,20),vmin=0,vmax=vmax,cmap='magma')
        title=label+'\n'+z['id'] if j==0 else ['LM','Fixed twofold','Adaptive'][j-1]+f"\nL2={z['methods'][ms[j-1]]['material_relative_l2_independent_data']:.1%}"
        axes[i,j].set(xlabel='x (cm)',ylabel='y (cm)',title=title,xlim=(-12,12),ylim=(-12,12))
    fig.colorbar(im,ax=axes[i,:],shrink=.85,label='Real contrast')
    checks[fam]={}
    for m,n in zip(ms,names):
        a=[q['methods'][m] for q in rows];term=collections.Counter(q['termination'] for q in a)
        for q in a:assert np.isfinite(q['material_relative_l2_independent_data']) and np.isfinite(q['held_scattered_relative'])
        lines.append(f"|{fl}|{n}|{np.median([q['material_relative_l2_independent_data'] for q in a]):.2%}|{np.median([q['held_scattered_relative'] for q in a]):.2%}|{np.median([q['wall_seconds'] for q in a]):.2f}|{term['rejected_step']}/{term['iteration_cap']}|")
        checks[fam][m]={'terminations':dict(term),'state_residual_max':max(q['state_relative_residual'] for q in a)}
fig.savefig(ROOT/'figures/a2/imaging_typical_overview.png',dpi=165);plt.close(fig)
lines+=['','![三类目标的实际成像，左为真值](%s)'%(ROOT/'figures/a2/imaging_typical_overview.png'),'','图中每行选自该组自适应材料误差排序的中间案例，各方法使用同一个目标和初值；显示中央24cm区域，原计算区域为40cm，每行共用色标。这是按既定排序规则选择的展示，不用来单独判优。','',
'相邻目标上，自适应材料中位误差1.26%对LM的3.80%，留出场也较低；但只有6/10配对同时在材料和留出场上更好，耗时更长。较高对比度组LM材料更准且更快，自适应最终均扩至12个完整方向。因此值得研究的是早期限制方向是否改变非线性路径，而不是宣称内部传播增加信息或固定twofold已胜出。','',
'视觉上，相邻组的几张图很接近，不能凭外观宣布明显优势；较高对比度组固定twofold改变了目标形状/幅度；尖边真值有明确缺口，LM保留了部分凹陷却把它变成两个平滑团，自适应则几乎填平缺口。更快的尖边自适应伴随更高材料和留出场误差，不能称同精度加速。下一步需分离边界表示、优化路径与物理失配，不能只串联一个网络把图变锐。','',
'全部五方法、停止比例和配对胜数见 [完整统计](../runs/a2/imaging/SUMMARY_ZH.md)。各组最好/典型/最差完整图：[相邻](../figures/a2/imaging_gaussian_close.png)、[较高对比度](../figures/a2/imaging_gaussian_strong.png)、[尖边](../figures/a2/imaging_sharp_boundary.png)。这些是开发诊断；尖边原Gaussian代理参数误差已标为无效，没有进入上述材料分数。']
p=ROOT/'deliverables/RESEARCH_REPORT_ZH.md';s=p.read_text();a=s.index('<!-- IMAGING_SUMMARY_START -->');b=s.index('<!-- IMAGING_SUMMARY_END -->')+len('<!-- IMAGING_SUMMARY_END -->');s=s[:a]+'<!-- IMAGING_SUMMARY_START -->\n'+'\n'.join(lines)+'\n<!-- IMAGING_SUMMARY_END -->'+s[b:];p.write_text(s)
(out/'root_statistics_check.json').write_text(json.dumps({'cases':30,'fits':150,'scope':'Post-run root summary; no new success/failure thresholds.','checks':checks},indent=2))
manifest=json.loads((out/'POSTRUN_MANIFEST.json').read_text());manifest['frozen_runner_sha256']=json.loads((out/'frozen_config.json').read_text())['source_sha256'];manifest['postrun_current_python_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'code/a2').glob('*.py')};manifest['note']='Post-run only. Frozen runner digest is preserved separately; this does not backfill missing historical dependency snapshots.';(out/'POSTRUN_MANIFEST.json').write_text(json.dumps(manifest,indent=2));print({'cases':30,'fits':150})
