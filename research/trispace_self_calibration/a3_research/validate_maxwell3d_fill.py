"""Development correction: resolve boundary volume before calibration claims."""
import json
from pathlib import Path
import numpy as np
from maxwell3d import DipoleVIE,treams_field,receivers

if __name__=='__main__':
    rows=[]
    for name,centers,radii,eps in [('single',[[0,0,0]],[.12],[2.4+.03j]),('pair',[[-.2,0,0],[.2,.03,.04]],[.12,.10],[2.4+.03j,3+.04j])]:
        for k in [3.,6.,9.]:
            ref,_=treams_field(centers,radii,eps,k,receivers(),5)
            for h in [.06,.04,.03,.025]:
                m=DipoleVIE(centers,radii,h,k,fill_quadrature=6)
                y=m.field(eps,receivers())
                row=dict(scene=name,k=k,spacing=h,voxels=len(m.points),fill_quadrature=6,
                    relative_error=float(np.linalg.norm(y-ref)/np.linalg.norm(ref)),
                    volume_relative_error=float(abs(m.fill.sum()*h**3-sum(4*np.pi*np.asarray(radii)**3/3))/sum(4*np.pi*np.asarray(radii)**3/3)))
                rows.append(row);print(json.dumps(row),flush=True)
                (Path(__file__).resolve().parent/'results/maxwell3d_fill_development.json').write_text(json.dumps(rows,indent=2)+'\n')
