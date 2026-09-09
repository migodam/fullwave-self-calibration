"""Independent vector sphere-cluster check of far-field gain hiding at fixed SNR."""
import json
from pathlib import Path
import numpy as np
from maxwell3d import receivers,treams_field
from gain_graph import visible_geometry,span,real_parameter_jacobian

OUT=Path(__file__).resolve().parent
C=[[-.2,0,0],[.2,.03,.04]];RAD=[.12,.1];EPS=[2.4+.03j,3.+.04j]

def one(k,radius):
    rx=receivers(12,radius)
    order=8 if k>9 else 5
    def field(pts,eps=EPS):return treams_field(C,RAD,eps,k,pts,order)[0].reshape(36,4)
    h=field(rx); delta=1e-4
    b=np.stack([(field(rx+np.eye(3)[d]*delta)-field(rx-np.eye(3)[d]*delta))/(2*delta) for d in range(3)],axis=-1)
    # Independently half the displacement for a numerical derivative audit.
    bh=np.stack([(field(rx+np.eye(3)[d]*delta/2)-field(rx-np.eye(3)[d]*delta/2))/delta for d in range(3)],axis=-1)
    n=rx/radius
    leading=(1j*k-1/radius)*np.repeat(n,3,axis=0)[:,None,:]*h[:,:,None]
    sigma=float(np.linalg.norm(h)/np.sqrt(h.size)/np.sqrt(1000))
    scale=np.sqrt(2)/sigma
    edges=[(r,t) for r in range(36) for t in range(4)]
    vv=visible_geometry(h.ravel(),edges,b.reshape(-1,3))
    raw=real_parameter_jacobian(b.reshape(-1,3))*scale
    visible=vv['visible']*scale
    # Receiver-only gains (shared across illuminations, independent component).
    gg=np.zeros((144,36),complex)
    for r in range(36):gg[4*r:4*r+4,r]=h[r]
    from gain_graph import complex_real_matrix
    qg=span(complex_real_matrix(gg))
    rxvis=raw-qg@(qg.T@raw)
    gl=real_parameter_jacobian(leading.reshape(-1,3))*scale
    return dict(k=k,radius=radius,noise_sigma=sigma,
                derivative_halving_error=float(np.linalg.norm(b-bh)/np.linalg.norm(bh)),
                leading_gain_projector_residual=float(np.linalg.norm(gl-qg@(qg.T@gl))/np.linalg.norm(gl)),
                raw_singular_values=np.linalg.svd(raw,compute_uv=False).tolist(),
                rx_gain_profiled_singular_values=np.linalg.svd(rxvis,compute_uv=False).tolist(),
                rx_tx_gain_profiled_singular_values=np.linalg.svd(visible,compute_uv=False).tolist(),
                leading_subtraction_relative_norm=float(np.linalg.norm(b-leading)/np.linalg.norm(h)),
                relative_visible_norm=float(np.linalg.norm(vv['visible'])/np.linalg.norm(h)),
                target_rank=vv['target_rank'],noise_experiment='fixed 30dB global complex SNR; per-component variance sigma^2/2',
                scope='local receiver-translation derivative, material fixed, single-frequency free complex gains, independent Treams fields')

if __name__=='__main__':
    path=OUT/'results'/'far_field_gain.json'
    rows=json.loads(path.read_text()) if path.exists() else []
    for k in [9.,18.]:
        for radius in [1.3,2.6,5.2,10.4,20.8]:
            if any(r['k']==k and r['radius']==radius for r in rows):continue
            row=one(k,radius);rows.append(row)
            path.write_text(json.dumps(rows,indent=2)+'\n')
            print(json.dumps(row),flush=True)
