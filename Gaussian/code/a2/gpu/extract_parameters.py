"""Post-hoc, CPU-only extraction of completed Gaussian parameter checkpoints."""
from pathlib import Path
import json,hashlib,numpy as np,torch
ROOT=Path(__file__).resolve().parents[3]
RUN=ROOT/'runs/a2/gpu/representation_campaign/main128'
def main():
    records=[]
    for case in range(30):
        for k in (16,64,144):
            p=RUN/f'case_{case:02d}_gaussian_K{k}.pt'
            # These are this campaign's own local checkpoints, not third-party uploads.
            d=torch.load(p,map_location='cpu',weights_only=False)
            assert d['step']==120
            m=d['model'];a=torch.nn.functional.softplus(m['raw_amp']).numpy();c=(.18*torch.tanh(m['raw_center'])).numpy();v=m['raw_l'];L=np.zeros((k,2,2));L[:,0,0]=(.004+.04*torch.nn.functional.softplus(v[:,0])).numpy();L[:,1,0]=(.02*v[:,1]).numpy();L[:,1,1]=(.004+.04*torch.nn.functional.softplus(v[:,2])).numpy();cov=L@L.transpose(0,2,1);scales=np.sqrt(np.linalg.eigvalsh(cov));active=a>.01
            records.append({'case':case,'method':f'gaussian_K{k}','config_hash':d['config_hash'],'checkpoint_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'amplitude':a.tolist(),'center_m':c.tolist(),'covariance_m2':cov.tolist(),'principal_sigma_m':scales.tolist(),'component_count':k,'active_amplitude_gt_0p01_count':int(active.sum()),'min_principal_sigma_m':float(scales[:,0].min()),'active_subgrid_sigma_count':int(((scales[:,0]<.4/128)&active).sum()),'active_anisotropy_gt_10_count':int(((scales[:,1]/scales[:,0]>10)&active).sum())})
    out={'scope':'Post-hoc representation diagnostic, not causal proof of aliasing or false structures. Active amplitude >0.01 and sigma<h thresholds selected after viewing images. No optimization is changed.','grid_spacing_m':.4/128,'records':records,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (RUN/'gaussian_parameters_posthoc.json').write_text(json.dumps(out,indent=2));print(json.dumps({'records':len(records),'subgrid_active_components':sum(r['active_subgrid_sigma_count'] for r in records)}))
if __name__=='__main__':main()
