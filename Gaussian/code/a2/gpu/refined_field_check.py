"""Post-hoc N256 point-kernel field reassessment; no optimization or selection."""
from pathlib import Path
import hashlib,json,time,numpy as np,torch
import torch.nn.functional as fn
ROOT=Path(__file__).resolve().parents[3]
import sys;sys.path.insert(0,str(ROOT/'code'))
from a2.gpu.representation_campaign import make_backends,grid_t,truth_real,LOSS,FREQ
RUN=ROOT/'runs/a2/gpu/representation_campaign/main128'
def main():
    start=time.perf_counter();n=256;dev=torch.device('cuda');pts=grid_t(n,dev);backs=make_backends(n,dev);setup=time.perf_counter()-start
    pp=json.loads((RUN/'gaussian_parameters_posthoc.json').read_text());params={(r['case'],r['method']):r for r in pp['records']};records=[]
    for case in range(30):
        y=torch.as_tensor(np.load(RUN/f'case_{case:02d}_data.npz')['fields'],device=dev);truth=truth_real(pts,case);held=torch.arange(1,64,2,device=dev)
        for name in ['gaussian_K16','gaussian_K64','gaussian_K144','voxel']:
            t=time.perf_counter()
            with torch.no_grad():
                if name=='voxel':
                    a=torch.as_tensor(np.load(RUN/f'case_{case:02d}_{name}_arrays.npz')['reconstruction'],device=dev);chi=fn.interpolate(a.reshape(1,1,128,128),size=(n,n),mode='nearest').reshape(-1)
                else:
                    p=params[(case,name)];a=torch.tensor(p['amplitude'],dtype=torch.float32,device=dev);c=torch.tensor(p['center_m'],dtype=torch.float32,device=dev);cov=torch.tensor(p['covariance_m2'],dtype=torch.float32,device=dev);L=torch.linalg.cholesky(cov);d=pts[:,None]-c[None];u0=d[:,:,0]/L[None,:,0,0];u1=(d[:,:,1]-L[None,:,1,0]*u0)/L[None,:,1,1];chi=(a[None]*torch.exp(-.5*(u0*u0+u1*u1))).sum(1)
                field=backs[0].scattered(chi.to(torch.complex64)*LOSS);err=float(torch.linalg.vector_norm(field[:,held]-y[:,held])/torch.linalg.vector_norm(y[:,held]));mat=float(torch.linalg.vector_norm(chi-truth)/torch.linalg.vector_norm(truth))
            old=json.loads((RUN/f'case_{case:02d}_{name}.json').read_text());records.append({'case':case,'method':name,'held_clean_N256':err,'held_clean_original_N128':old['held_clean_field_relative'],'material_relative_N256':mat,'seconds':time.perf_counter()-t})
        print('refined_case',case,flush=True)
    out={'kind':'POST-HOC independent-grid forward reassessment, no reoptimization','inverse_n':128,'reassessment_n':n,'data_n':192,'kernel':'N256 point/collocation with existing self term, compared to original N192 integrated data; same VIE software, not independent solver','gaussian_extension':'continuous Gaussian parameters evaluated on N256','voxel_extension':'piecewise constant N128 cells, nearest-neighbor 2x2 prolongation','scope':'Mesh sensitivity diagnostic; not a continuous-model error certificate or frozen confirmation test','parameter_json_sha256':hashlib.sha256((RUN/'gaussian_parameters_posthoc.json').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dependencies_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'code/a2/physics.py',ROOT/'code/a2/gpu/torch_vie.py',ROOT/'code/a2/gpu/representation_campaign.py']},'setup_seconds':setup,'total_seconds':time.perf_counter()-start,'counts':backs[0].counts,'records':records}
    (RUN/'refined_field_posthoc.json').write_text(json.dumps(out,indent=2));print(json.dumps({'records':len(records),'seconds':out['total_seconds']}))
if __name__=='__main__':main()
