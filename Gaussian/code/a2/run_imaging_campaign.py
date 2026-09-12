"""Independent-grid, bounded A2 imaging development campaign (not frozen test)."""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[k]='2'
import argparse, hashlib, json, sys, time
from dataclasses import asdict, dataclass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.ports import render
from a2.som import FREQUENCIES_HZ,SCALES,components,fit_parameter_tangent_som,FitConfig

@dataclass(frozen=True)
class CampaignConfig:
    inverse_n:int=64; data_n:int=128; n_tx:int=12; n_rx:int=48; aperture:str='half'; quad:int=4
    objects_per_family:int=10; noise_fraction:float=.01; seed:int=20260914; max_iter:int=30
    first_rank:int=6; second_rank:int=2; thread_cap:int=2

def sharp_field(points,center,radius,notch):
    d=points-np.asarray(center); disk=(d[:,0]**2+d[:,1]**2<=radius**2)
    # A disk with an intentionally non-Gaussian rectangular notch removed.
    cut=(np.abs(d[:,0]-radius*.35)<notch)&(d[:,1]>-notch)
    return .62*(disk & ~cut)

def make_cases(c):
    rng=np.random.default_rng(c.seed);out=[]
    for family in ('gaussian_close','gaussian_strong','sharp_boundary'):
      for i in range(c.objects_per_family):
        shift=rng.uniform(-.025,.025,2)
        if family=='gaussian_close':
          a,sep=.42,.046; truth=np.array([a,-sep/2+shift[0],shift[1],np.log(.023),np.log(.028),.15,a*.85,sep/2+shift[0],shift[1]+.012,np.log(.026),np.log(.022),-.23]); source='gaussian'
        elif family=='gaussian_strong':
          a,sep=.72,.09; truth=np.array([a,-sep/2+shift[0],shift[1],np.log(.020),np.log(.035),.31,a*.72,sep/2+shift[0],shift[1]-.015,np.log(.032),np.log(.021),-.18]); source='gaussian'
        else:
          truth=np.array([.38,shift[0],shift[1],np.log(.04),np.log(.04),0.,.22,shift[0]+.025,shift[1],np.log(.03),np.log(.03),0.]);source='sharp'
        init=truth.copy();init[::6]*=.62;init[1::6]+=rng.normal(0,.012,2);init[2::6]+=rng.normal(0,.012,2);init[3::6]+=rng.normal(0,.17,2);init[4::6]+=rng.normal(0,.17,2)
        out.append({'id':f'{family}_{i:02d}','family':family,'source':source,'truth_theta':truth.tolist(),'initial':init.tolist(),'shift':shift.tolist(),'noise_seed':int(c.seed+1000+len(out))})
    return out

def chi(case,vie,f):
    if case['source']=='gaussian': return render(components(np.asarray(case['truth_theta']),f),vie.points)
    return sharp_field(vie.points,case['shift'],.055,.014)*(1+.02j*2.25e9/f)

def data(case,data_vies,c):
    clean=[]
    for f,v in zip(FREQUENCIES_HZ,data_vies): clean.append(v.forward(chi(case,v,f),rtol=2e-8))
    rng=np.random.default_rng(case['noise_seed']);allv=np.concatenate([x['scattered'].ravel() for x in clean]);sig=c.noise_fraction*np.linalg.norm(allv)/np.sqrt(allv.size)
    obs=[x['scattered']+sig/np.sqrt(2)*(rng.normal(size=x['scattered'].shape)+1j*rng.normal(size=x['scattered'].shape)) for x in clean]
    return clean,obs,sig

def material_truth(case,vie): return (render(components(np.asarray(case['truth_theta']),FREQUENCIES_HZ[0]),vie.points).real if case['source']=='gaussian' else sharp_field(vie.points,case['shift'],.055,.014))
def evaluate(rec,case,inv_vie):
    got=render(components(np.asarray(rec['parameters']),FREQUENCIES_HZ[0]),inv_vie.points).real;true=material_truth(case,inv_vie)
    rec['material_relative_l2_independent_data'] = float(np.linalg.norm(got-true)/np.linalg.norm(true));return got,true

def figures(records,inv_vie,figdir):
    # Preregistered representative selector: min/median/max adaptive material error within family.
    for family in ('gaussian_close','gaussian_strong','sharp_boundary'):
      rows=[r for r in records if r['family']==family]; rows.sort(key=lambda r:r['methods']['adaptive_tsomg']['material_relative_l2_independent_data'])
      picks=[rows[0],rows[len(rows)//2],rows[-1]];methods=('ordinary_lm','fixed_twofold_tsomg','adaptive_tsomg','random_complement_tsomg','gsvd_ratio_tsomg')
      fig,ax=plt.subplots(3,6,figsize=(16,8),constrained_layout=True); vmax=max(np.max(material_truth(x,inv_vie)) for x in picks)
      for i,row in enumerate(picks):
        truth=material_truth(row,inv_vie).reshape(inv_vie.geometry.n,inv_vie.geometry.n);ax[i,0].imshow(truth,origin='lower',vmin=0,vmax=vmax,cmap='magma');ax[i,0].set_title(('best','typical','worst')[i]+' truth')
        for j,m in enumerate(methods,1):
          got=render(components(np.asarray(row['methods'][m]['parameters']),FREQUENCIES_HZ[0]),inv_vie.points).real.reshape(inv_vie.geometry.n,inv_vie.geometry.n);ax[i,j].imshow(got,origin='lower',vmin=0,vmax=vmax,cmap='magma');ax[i,j].set_title(m.replace('_',' ')+'\n%.1f%%'%(100*row['methods'][m]['material_relative_l2_independent_data']))
        for a in ax[i]:a.set_xticks([]);a.set_yticks([])
      fig.savefig(figdir/f'imaging_{family}.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4),constrained_layout=True); methods=list(records[0]['methods'])
    for m in methods:
      x=np.array([r['methods'][m]['material_relative_l2_independent_data'] for r in records]);t=np.array([r['methods'][m]['wall_seconds'] for r in records]);ax[0].scatter(t,x,label=m,s=15);ax[1].hist(x,bins=12,alpha=.45,label=m)
    ax[0].set(xlabel='wall seconds',ylabel='material relative L2');ax[1].set(xlabel='material relative L2',ylabel='cases');ax[0].legend(fontsize=6);ax[1].legend(fontsize=6);fig.savefig(figdir/'imaging_timing_failures.png',dpi=150);plt.close(fig)

def main():
 p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');p.add_argument('--resume',action='store_true');a=p.parse_args();c=CampaignConfig();out=ROOT/'runs/a2/imaging';figdir=ROOT/'figures/a2';out.mkdir(parents=True,exist_ok=True);figdir.mkdir(parents=True,exist_ok=True)
 source={'config':asdict(c),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'protocol':'independent N128 cell-integrated data / N64 point inverse; development only, not frozen900','representative_selector':'within-family adaptive material score min/median/max','methods':['ordinary_lm','fixed_twofold_tsomg','adaptive_tsomg','random_complement_tsomg','gsvd_ratio_tsomg']};(out/'frozen_config.json').write_text(json.dumps(source,indent=2))
 ig=Geometry(n=c.inverse_n,n_tx=c.n_tx,n_rx=c.n_rx,aperture=c.aperture);dg=Geometry(n=c.data_n,n_tx=c.n_tx,n_rx=c.n_rx,aperture=c.aperture);inv=[VIE(ig,f) for f in FREQUENCIES_HZ];dat=[VIE(dg,f,cell_integrated=True,quadrature_order=c.quad) for f in FREQUENCIES_HZ];train=np.arange(c.n_rx)%2==0;held=~train
 records=json.loads((out/'results.json').read_text()) if a.resume and (out/'results.json').exists() else []
 for case in make_cases(c)[:1 if a.pilot else 3*c.objects_per_family]:
  if any(r['id']==case['id'] for r in records):continue
  clean,obs,sig=data(case,dat,c); # inverse uses the independent data directly, not an N64-generated surrogate
  methods={};spec=[('ordinary_lm',FitConfig(max_iter=c.max_iter)),('fixed_twofold_tsomg',FitConfig(max_iter=c.max_iter,first_rank=c.first_rank,second_rank=c.second_rank)),('adaptive_tsomg',FitConfig(max_iter=c.max_iter,first_rank=4,second_rank=0,adaptive_rank=True)),('random_complement_tsomg',FitConfig(max_iter=c.max_iter,first_rank=c.first_rank,second_rank=c.second_rank)),('gsvd_ratio_tsomg',FitConfig(max_iter=c.max_iter,first_rank=c.first_rank,second_rank=c.second_rank))]
  method_dir=out/'checkpoints'/case['id'];method_dir.mkdir(parents=True,exist_ok=True)
  for name,cfg in spec:
    checkpoint=method_dir/(name+'.json')
    if a.resume and checkpoint.exists():
      r=json.loads(checkpoint.read_text())
    else:
      call='twofold_tsomg' if name=='fixed_twofold_tsomg' else name; r=fit_parameter_tangent_som(call,np.asarray(case['initial']),np.asarray(case['truth_theta']),obs,clean,inv,train,held,sig,cfg);r['method']=name;evaluate(r,case,inv[0]);checkpoint.write_text(json.dumps(r,indent=2))
      if r['wall_seconds']>300: print(json.dumps({'warning':'method exceeded 5 minutes; retained fixed cap','case':case['id'],'method':name,'wall_seconds':r['wall_seconds']}),flush=True)
    methods[name]=r
  rec={**case,'noise_sigma':float(sig),'methods':methods};records.append(rec);(out/'results.json').write_text(json.dumps(records,indent=2));print(json.dumps({'id':case['id'],'scores':{k:v['material_relative_l2_independent_data'] for k,v in methods.items()}},indent=2),flush=True)
 if len(records)>=3*c.objects_per_family: figures(records,inv[0],figdir)
if __name__=='__main__':main()
