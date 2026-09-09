"""Post-fit audit separating deployable prediction and oracle-model diagnostics.

Frozen estimates only: no refitting on held-out receivers. Original metrics
are preserved, not overwritten. Noiseless independent held-out field is used
to evaluate mean-prediction error; it is not an acquired calibration reference.
"""
import gc
import json
import sys
import time
from pathlib import Path
import numpy as np
from maxwell3d import dipole_kernel, receivers, treams_field
import calibrate3d as cal

OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT.parents[1]/'delegated/a3_maxwell_refine'))
from tangent_fft import TangentFFTVIE


def metrics(pred,truth):
    phase=np.angle(pred*truth.conj())
    return dict(field_relative_error=float(np.linalg.norm(pred-truth)/np.linalg.norm(truth)),
                weighted_phase_rmse_rad=float(np.sqrt(np.sum(abs(truth)**2*phase**2)/np.sum(abs(truth)**2))))


def audit(row):
    z=np.array(row['estimated']);true=np.array(row['true']);ks=np.array(row['frequencies'])
    held=receivers(17,1.6);estimated=[];t0=time.perf_counter()
    for k in ks:
        model=TangentFFTVIE(cal.CENTERS,cal.RADII,row['spacing'],k,fill_quadrature=6)
        current=model.currents(z[:2]+1j*cal.LOSS,False)
        estimated.append((dipole_kernel(held+z[2:5],model.points,k)@current).reshape(17,3,4))
        del model;gc.collect()
    estimated=np.array(estimated)
    truth=np.array([treams_field(cal.CENTERS,cal.RADII,true[:2]+1j*cal.LOSS,k,held+true[2:5],8 if k>9 else 5)[0] for k in ks])
    electronics=lambda a:np.exp(a[6:10]+1j*a[10:14]+1j*ks[:,None,None,None]*a[5])
    return dict(seed=row['seed'],spacing=row['spacing'],frequencies=ks.tolist(),method=row['method'],
                deployable_sensor_prediction=metrics(estimated*electronics(z),truth*electronics(true)),
                reconstructed_structural_field=metrics(estimated,truth),
                preserved_oracle_parameter_transfer=dict(field_relative_error=row['heldout_field_relative_error'],
                    weighted_phase_rmse_rad=row['heldout_phase_rmse_rad']),
                audit_seconds=time.perf_counter()-t0,
                scope='frozen estimates; new receiver angles; actual inverse solver and fitted electronics; independent noiseless mean target; no refit')


if __name__=='__main__':
    source=OUT.parents[1]/'delegated/a3_maxwell_refine/calibration_fft_raw.json'
    destination=OUT/'results/prediction_audit3d.json'
    rows=json.loads(destination.read_text()) if destination.exists() else []
    key=lambda r:(r['seed'],r['spacing'],tuple(r['frequencies']),r['method'])
    for row in json.loads(source.read_text()):
        if key(row) in {key(r) for r in rows}:continue
        result=audit(row);rows.append(result)
        destination.write_text(json.dumps(rows,indent=2)+'\n')
        print(json.dumps(result),flush=True)
