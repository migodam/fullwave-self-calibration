"""DEVELOPMENT independent 3D discretizer validation, no inverse claims."""
import json
import time
from pathlib import Path
import numpy as np
from maxwell3d import DipoleVIE,treams_field,receivers,dipole_kernel


def main():
    rows=[]; rx=receivers()
    scenes=[('single',[[0.,0.,0.]],[.12],[2.4+.03j]),
            ('pair',[[-.2,0.,0.],[.2,.03,.04]],[.12,.10],[2.4+.03j,3.0+.04j])]
    for name,centers,radii,eps in scenes:
        for k in [3.,6.,9.]:
            refs=[]
            for l in [3,5]:
                ref,xs=treams_field(centers,radii,eps,k,rx,l)
                refs.append(ref)
            convergence=np.linalg.norm(refs[0]-refs[1])/np.linalg.norm(refs[1])
            for h in [.06,.04,.03]:
                start=time.perf_counter()
                model=DipoleVIE(centers,radii,h,k)
                pred=model.field(eps,rx)
                rel=np.linalg.norm(pred-refs[1])/np.linalg.norm(refs[1])
                kernel=model.kernel
                rec=np.linalg.norm(kernel-kernel.T)/max(np.linalg.norm(kernel),1e-30)
                row=dict(scene=name,k=k,spacing=h,voxels=len(model.points),reference_lmax=5,
                    reference_relative_convergence=float(convergence),relative_field_error=float(rel),
                    reciprocity_error=float(rec),xs_sca_ext=xs,
                    elapsed_seconds=time.perf_counter()-start,scope='development forward only')
                rows.append(row); print(json.dumps(row),flush=True)
                out=Path(__file__).resolve().parent/'results';out.mkdir(exist_ok=True)
                (out/'maxwell3d_forward_development.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':main()
