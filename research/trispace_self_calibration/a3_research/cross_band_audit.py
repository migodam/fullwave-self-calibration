"""Frozen-estimate per-frequency audit; never refits or uses pooled frequency errors."""
import gc
import json
from pathlib import Path
import time
import numpy as np
from prediction_audit3d import TangentFFTVIE,metrics,cal,dipole_kernel,receivers,treams_field

OUT=Path(__file__).resolve().parent


def audit(row):
    z=np.array(row['estimated']);true=np.array(row['true']);held=receivers(17,1.6)
    per=[];started=time.perf_counter()
    for k in row['frequencies']:
        model=TangentFFTVIE(cal.CENTERS,cal.RADII,.01,k,fill_quadrature=6)
        current=model.currents(z[:2]+1j*cal.LOSS,False)
        pred=(dipole_kernel(held+z[2:5],model.points,k)@current).reshape(17,3,4)
        target=treams_field(cal.CENTERS,cal.RADII,true[:2]+1j*cal.LOSS,k,held+true[2:5],8 if k>9 else 5)[0]
        gain=lambda a:np.exp(a[6:10]+1j*a[10:14]+1j*k*a[5])
        per.append(dict(k=k,sensor=metrics(pred*gain(z),target*gain(true)),structural=metrics(pred,target)))
        del model;gc.collect()
    return dict(seed=row['seed'],method=row['method'],frequencies=row['frequencies'],spacing=.01,
        per_frequency=per,seconds=time.perf_counter()-started,
        scope='frozen DEVELOPMENT estimates; actual inverse predictions, no refit; not a matched-noise causal experiment')


def main():
    source=OUT.parents[1]/'delegated/a3_maxwell_refine/calibration_fft_raw.json'
    dest=OUT/'results/cross_band_audit.json'
    source_rows=[r for r in json.loads(source.read_text()) if r['spacing']==.01]
    assert len(source_rows)==8
    rows=json.loads(dest.read_text()) if dest.exists() else []
    key=lambda r:(r['seed'],r['method'],tuple(r['frequencies']))
    for row in source_rows:
        if key(row) in {key(r) for r in rows}:continue
        result=audit(row);rows.append(result)
        dest.write_text(json.dumps(rows,indent=2)+'\n')
        print(json.dumps(dict(completed=len(rows),planned=8,seed=row['seed'],method=row['method'])),flush=True)


if __name__=='__main__':main()
