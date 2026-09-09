"""Read-only analysis of frozen A4 endpoints. Never refits or edits raw results."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
def load(name): return json.loads((OUT/(name+'.json')).read_text())
def cp_upper(k,n): return float(beta.ppf(.95,k+1,n-k)) if k<n else 1.
def dumps(path,data): path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
manifest=json.loads((ROOT/'provenance/frozen_manifest.json').read_text())
checks={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in manifest['source_hashes'].items()}
assert all(checks.values()),checks
cal=load('frozen_calibration')['rows']; br=load('frozen_branch')['rows']; spec=load('frozen_spectral'); dda=load('frozen_dda')['rows']
methods=sorted({r['method'] for r in cal}); policies=sorted({r['policy'] for r in br}); kinds=sorted({r['kind'] for r in br})
cs=[]
for m in methods:
 rows=[r for r in cal if r['method']==m]
 cs.append(dict(method=m,n=len(rows),median_geometry_mm=float(np.median([r['geometry_error_m'] for r in rows])*1000),median_material_percent=float(np.median([r['material_relative_error'] for r in rows])*100),mean_scaled_task_loss=float(np.mean([r['scaled_task_loss'] for r in rows])),median_total_seconds=float(np.median([r['total_seconds'] for r in rows])),all_converged=all(r['fit_success'] for r in rows)))
bs=[]
for p in policies:
 for kind in ['all']+kinds:
  rows=[r for r in br if r['policy']==p and (kind=='all' or r['kind']==kind)]
  sums={s:sum(int(r[s]) for r in rows) for s in ['accepted','correct','wrong_accepted','correct_rejected','initial_bank_missing','final_bank_covered','wrong_selection_when_covered']}
  bs.append(dict(policy=p,kind=kind,n=len(rows),**sums,median_geometry_mm=float(np.median([r['geometry_error_m'] for r in rows])*1000),median_fit_seconds=float(np.median([r['fit_seconds'] for r in rows])),median_additional_wall_seconds=float(np.median([r['additional_wall_seconds'] for r in rows])),unresolved_cells=[r.get('coverage',{}).get('unresolved_cells',0) for r in rows]))
a=sorted([r for r in cal if r['method']=='risk_controller'],key=lambda r:r['scene']); b=sorted([r for r in cal if r['method']=='full_fine'],key=lambda r:r['scene'])
ratios=np.array([x['total_seconds']/y['total_seconds'] for x,y in zip(a,b)])
rng=np.random.default_rng(20260909); boot=np.median(rng.choice(ratios,size=(10000,len(ratios)),replace=True),axis=1)
paired={'controller_actions':[r['selected_method'] for r in a],'median_paired_cost_ratio':float(np.median(ratios)),'paired_ratios':ratios.tolist(),'descriptive_bootstrap_95_interval':np.quantile(boot,[.025,.975]).tolist(),'max_parameter_difference':float(max(np.max(abs(np.array(x['estimate'])-y['estimate'])) for x,y in zip(a,b))),'note':'Only eight scenes. Timing includes the independent pilot and charged construction/screening; hardware mode-synthesis cost unavailable.'}
valid=[r for r in br if r['policy']=='coverage_aware' and r['kind'] in ['easy','antipodal','bank_missing','low_snr']]
validstats={'n':len(valid),'accepted':sum(r['accepted'] for r in valid),'wrong_accepted':sum(r['wrong_accepted'] for r in valid),'zero_events_one_sided_95_upper_if_iid_population':cp_upper(sum(r['wrong_accepted'] for r in valid),len(valid)),'interpretation':'Stratified finite demonstration, not an iid population assurance or a calibrated sub-percent empirical failure rate.'}
dss=[]
for s in sorted({r['shape'] for r in dda}):
 rows=[r for r in dda if r['shape']==s]; ratio=np.array([min(r['profiled_singular_values'])/r['sphere_formula_min'] for r in rows])
 dss.append({'shape':s,'receiver_cases':len(rows),'min_singular_ratio_min_median_max':[float(f(ratio)) for f in [np.min,np.median,np.max]],'pattern_residual_min_max':[float(min(r['pattern_residual'] for r in rows)),float(max(r['pattern_residual'] for r in rows))],'nuisance_ranks':sorted({r['nuisance_rank'] for r in rows})})
ref=load('reference_v2_fresh'); rs={'fresh_scenes':len(ref),'median_material_percent':float(np.median([r['material_relative_error'] for r in ref])*100),'max_material_visible_norm_without_reference':max(r['material_visible_norm_without_reference'] for r in ref),'min_material_visible_norm_with_reference':min(r['material_visible_norm_with_reference'] for r in ref),'maximum_no_reference_profile_loss_range':max(r['no_reference_profile_loss_range'] for r in ref),'scope':'geometry known externally; three independent noisy complex gain references'}
summary={'calibration':cs,'branch':bs,'controller_vs_fine':paired,'exact_modal_coverage_strata':validstats,'dda':dss,'reference_v2':rs,'spectral':{k:v for k,v in spec.items() if k!='rows'},'completion':load('frozen_complete')}
dumps(OUT/'A4_ANALYSIS.json',summary)
# Source and raw-data integrity audit. No unrun test is reported as passed.
npz_info={}
for p in [OUT/'frozen_branch_data.npz',OUT/'frozen_calibration_data.npz']:
 with np.load(p) as z: npz_info[p.name]={'arrays':len(z.files),'all_finite':all(np.all(np.isfinite(z[k])) for k in z.files),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
assert len(cal)==72 and len(br)==144 and len(dda)==96
assert all(x['all_finite'] for x in npz_info.values())
dumps(OUT/'A4_AUDIT.json',{'frozen_sources_match':checks,'raw_npz':npz_info,'record_counts':{'calibration':len(cal),'branch':len(br),'dda_receiver_cases':len(dda),'spectral':len(spec['rows']),'diversity':len(load('frozen_diversity'))},'known_withdrawal':'frozen_reference.json material_visible_norm_without_reference has a row-index diagnostic bug; use separately preregistered reference_v2_fresh.json. Original fits are retained.','certification':'No interval arithmetic, continuum certificate, A3 rerun or hardware result is claimed.'})
# Figures use Matplotlib defaults; one independent figure per plot.
z=np.array([r['z'] for r in spec['rows']]); vals=np.array([r['theory'] for r in spec['rows']]); nums=np.array([r['numerical'] for r in spec['rows']])
fig,ax=plt.subplots(figsize=(7.6,4.8)); ax.loglog(z,vals[:,0],label='Tangential: analytic');ax.loglog(z,vals[:,2],label='Radial: analytic');ax.loglog(z[::4],nums[::4,2],'o',fillstyle='none',label='Profiled Jacobian');ax.axvline(1,linestyle=':',label='kR = 1');ax.set(xlabel='Dimensionless range kR (fixed k and relative SNR)',ylabel='Visible singular value',title='Positive radial information has two limiting degeneracies');ax.legend();ax.grid(True,which='both',alpha=.2);fig.tight_layout();fig.savefig(FIG/'01_modal_information.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.8))
for s in sorted({r['shape'] for r in dda}):
 rows=[r for r in dda if r['shape']==s]; ax.scatter([r['pattern_residual'] for r in rows],[min(r['profiled_singular_values'])/r['sphere_formula_min'] for r in rows],label=s,alpha=.75)
ax.axhline(1,linestyle=':');ax.set_xscale('log');ax.set(xlabel='Relative residual after best scalar modal fit',ylabel='Material-profiled minimum / ideal modal minimum',title='A small field mismatch is not a profiled-information guarantee');ax.legend();ax.grid(True,alpha=.2);fig.tight_layout();fig.savefig(FIG/'02_shape_stress.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,5));pos=np.arange(len(cs)); ax.barh(pos,[r['mean_scaled_task_loss'] for r in cs]);ax.set_yticks(pos,[r['method'] for r in cs]);ax.set_xscale('log');ax.set(xlabel='Mean scaled task loss (8 matched scenes; log scale)',title='Strong approximation-error and fine-model baselines');ax.invert_yaxis();ax.grid(True,axis='x',alpha=.2);fig.tight_layout();fig.savefig(FIG/'03_task_risk_baselines.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(7.8,4.7));ax.plot(range(1,9),ratios,'o-');ax.axhline(1,linestyle=':');ax.set(xlabel='Frozen scene',ylabel='Controller / full-fine charged wall time',title='Identical endpoints, extra controller cost');ax.grid(True,alpha=.2);fig.tight_layout();fig.savefig(FIG/'04_controller_cost.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(8.7,4.8));sel=[r for r in bs if r['policy']=='coverage_aware' and r['kind']!='all'];x=np.arange(len(sel));ax.bar(x-.18,[r['accepted'] for r in sel],width=.36,label='Accepted');ax.bar(x+.18,[r['wrong_accepted'] for r in sel],width=.36,label='Wrong and accepted');ax.set_xticks(x,[r['kind'].replace('_','\n') for r in sel]);ax.set(ylabel='Scene count (4 per stratum)',title='Rejection works only within its physical assumptions',ylim=(0,4.7));ax.legend();fig.tight_layout();fig.savefig(FIG/'05_coverage_strata.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(7.7,4.7));ax.semilogy([r['material_visible_norm_without_reference'] for r in ref],'o',label='Without reference (FD roundoff)');ax.semilogy([r['material_visible_norm_with_reference'] for r in ref],'s',label='With three noisy complex references');ax.set(xlabel='Fresh reference-audit scene',ylabel='Profiled material sensitivity',title='Reference removes a specific material-gain ambiguity');ax.legend();ax.grid(True,alpha=.2);fig.tight_layout();fig.savefig(FIG/'06_reference_rank.png',dpi=180);plt.close(fig)
lines=['# Executed A4 numerical summary','', '## Nine matched baselines','', '| Method | Median geometry (mm) | Median material (%) | Mean scaled task loss | Median charged wall time (s) |','|---|---:|---:|---:|---:|']
for r in cs: lines.append(f"| {r['method']} | {r['median_geometry_mm']:.4g} | {r['median_material_percent']:.4g} | {r['mean_scaled_task_loss']:.5g} | {r['median_total_seconds']:.5g} |")
lines+=['','## Coverage-aware policy by stratum','','| Stratum | Scenes | Accepted | Wrong accepted | Correct rejected |','|---|---:|---:|---:|---:|']
for r in sel: lines.append(f"| {r['kind']} | {r['n']} | {r['accepted']} | {r['wrong_accepted']} | {r['correct_rejected']} |")
lines+=['','## Paired controller/fine result','',json.dumps(paired,indent=2),'','## Exact-modal strata','',json.dumps(validstats,indent=2),'','## Reference correction, fresh sample','',json.dumps(rs,indent=2)]
(ROOT/'docs/NUMERICAL_TABLES.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(summary,indent=2)[:22000])
