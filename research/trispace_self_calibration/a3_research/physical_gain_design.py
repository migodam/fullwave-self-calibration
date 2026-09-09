"""Development: sparse-link gain obstruction on an actual 3D Maxwell discretizer.

Known two-region material family; no global calibration or acquisition-win claim.
Topology choices below are fixed before seeing fields. No truth-based selection.
"""
import json
from pathlib import Path
import numpy as np
from calibrate3d import JointModel
from gain_graph import visible_geometry

if __name__=='__main__':
    model=JointModel(.04)
    z=np.zeros(14);z[:2]=[2.4,3.];z[2:5]=[.03,-.02,.01]
    y,j=model.field_jac(z)
    rows=[]
    for fi,k in enumerate([3.,6.,9.]):
        h=y[fi].reshape(-1,4)
        b=j[fi,...,2:5].reshape(-1,4,3)
        a=j[fi,...,:2].reshape(-1,4,2)
        nr,nt=h.shape
        tree=[(r,0) for r in range(nr)]+[(0,t) for t in range(1,nt)]
        designs={'tree':tree,'one_cycle':tree+[(1,1)],'two_cycles':tree+[(1,1),(2,1)],
                 'three_cycles':tree+[(1,1),(2,1),(3,1)],
                 'full':[(r,t) for r in range(nr) for t in range(nt)]}
        for name,edges in designs.items():
            rr,tt=np.array(edges).T
            yy=h[rr,tt];bb=b[rr,tt];aa=a[rr,tt]
            sigma=np.linalg.norm(h)/np.sqrt(h.size)*10**(-30/20)
            for unknown_material in [False,True]:
                out=visible_geometry(yy,edges,bb,material=aa if unknown_material else None,
                    covariance=np.eye(2*len(edges))*sigma*sigma/2,target_scale=[.01]*3)
                row={key:(value.tolist() if isinstance(value,np.ndarray) else value) for key,value in out.items() if key!='visible'}
                row.update(k=k,design=name,unknown_material=unknown_material,observed_complex_channels=len(edges),
                    pose_unit_m=.01,noise_sigma=float(sigma),scope='local development; physical-state q_material=2,r_free=0')
                rows.append(row)
    path=Path(__file__).resolve().parent/'results/physical_gain_design.json'
    path.write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))
