from pathlib import Path
import json,hashlib,sys,time,numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R/'code'))
from a2.physics import Geometry,VIE
from a2.ports import render
from a2.som import FREQUENCIES_HZ,components
from a2.run_imaging_campaign import material_truth
out=R/'runs/a2/imaging';fig=R/'figures/a2';
while not (out/'results.json').exists() or len(json.loads((out/'results.json').read_text()))<30: time.sleep(30)
r=json.loads((out/'results.json').read_text()); inv=VIE(Geometry(n=64,n_tx=12,n_rx=48,aperture='half'),FREQUENCIES_HZ[0]);methods=list(r[0]['methods']); summary={}
for row in r:
 for m,x in row['methods'].items():
  if row['family']=='sharp_boundary':x['gaussian_proxy_valid']=False;x['gaussian_proxy_invalid_reason']='non-Gaussian disk-notch truth; use material_relative_l2_independent_data'
  else:x['gaussian_proxy_valid']=True
for fam in sorted({x['family'] for x in r}):
 summary[fam]={}
 for m in methods:
  a=[x['methods'][m] for x in r if x['family']==fam];
  summary[fam][m]={'n':len(a),'material_median':float(np.median([z['material_relative_l2_independent_data'] for z in a])),'held_median':float(np.median([z['held_scattered_relative'] for z in a])),'wall_median':float(np.median([z['wall_seconds'] for z in a])),'failure_or_cap_fraction':float(np.mean([z['termination'] in ('iteration_cap','rejected_step') for z in a]))}
 base=[x['methods']['ordinary_lm']['material_relative_l2_independent_data'] for x in r if x['family']==fam]
 for m in methods: summary[fam][m]['paired_material_wins_vs_lm']=int(sum(z['material_relative_l2_independent_data']<b for z,b in zip([x['methods'][m] for x in r if x['family']==fam],base)))
(out/'results.json').write_text(json.dumps(r,indent=2));(out/'summary.json').write_text(json.dumps(summary,indent=2))
# compact typical + supplement
for fam in summary:
 rows=[x for x in r if x['family']==fam];rows.sort(key=lambda x:x['methods']['adaptive_tsomg']['material_relative_l2_independent_data']);picks=[rows[0],rows[len(rows)//2],rows[-1]]
 for compact,sel in [(True,[picks[1]]),(False,picks)]:
  ms=['ordinary_lm','fixed_twofold_tsomg','adaptive_tsomg'] if compact else methods;f,ax=plt.subplots(len(sel),1+len(ms),figsize=(10,3*len(sel)),constrained_layout=True);ax=np.atleast_2d(ax); vmax=max(material_truth(z,inv).max() for z in sel)
  for i,z in enumerate(sel):
   arr=[material_truth(z,inv)]+[render(components(np.array(z['methods'][m]['parameters']),FREQUENCIES_HZ[0]),inv.points).real for m in ms]
   for j,q in enumerate(arr):
    im=ax[i,j].imshow(q.reshape(64,64).T,origin='lower',extent=[-20,20,-20,20],vmin=0,vmax=vmax,cmap='magma');ax[i,j].set(xlabel='x (cm)',ylabel='y (cm)',title=(('typical' if compact else ('best','typical','worst')[i])+' '+z['id']) if j==0 else ms[j-1])
   f.colorbar(im,ax=ax[i,:],shrink=.7)
  f.savefig(fig/(f'imaging_{fam}_compact.png' if compact else f'imaging_{fam}.png'),dpi=150);plt.close(f)
lines=['# A2 独立网格开发成像汇总','','本次是 3 个固定形态模板 × 10 个随机平移/噪声/真值暖启动的局部恢复诊断；不是盲反演、30 个独立形状目标或冻结 900 例验收。N128 单元积分生成数据，N64 反演。','']
for fam,d in summary.items():
 lines+=['\n## '+fam,'','|方法|真值材料L2中位数|held中位数|时间中位数(s)|拒步停止或达cap|对LM材料胜数|','|---|---:|---:|---:|---:|---:|']+[f"|{m}|{v['material_median']:.2%}|{v['held_median']:.2%}|{v['wall_median']:.2f}|{v['failure_or_cap_fraction']:.0%}|{v['paired_material_wins_vs_lm']}/10|" for m,v in d.items()]
lines+=['','限制：LM 墙钟含 tangent/B/SVD 诊断，非优化纯 LM 基线；GSVD 控制未欧氏归一且 I 阻尼改变系数度量，非严格同参数度量消融。adaptive 的梯度缺口只是启发式，且最终常增加到更多方向，不能归因第二层或信息新增。sharp 的旧 Gaussian 参数代理诊断保留但无效，只使用独立网格材料 L2。']
(out/'SUMMARY_ZH.md').write_text('\n'.join(lines))
(out/'POSTRUN_MANIFEST.json').write_text(json.dumps({'postrun':True,'result_sha256':hashlib.sha256((out/'results.json').read_bytes()).hexdigest(),'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [out/'summary.json',out/'SUMMARY_ZH.md']},'note':'post-run manifest; frozen runner hash remains in frozen_config.json'},indent=2))
