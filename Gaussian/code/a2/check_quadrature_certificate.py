"""Independent dense Hermitian checks for the GL-kernel A2 extension."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[key]='1'
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,FFTGreen
from a2.certificate import gamma_for_fft_kernel

def main():
    rows=[]
    for n in (8,12):
        geom=Geometry(n=n);h=geom.side/n
        for kh in (.2,1.,2.9):
            for integrated in (False,True):
                D=FFTGreen(geom,kh/h,cell_integrated=integrated)
                chi=np.full(n*n,.5+.2j)
                gamma,meta=gamma_for_fft_kernel(chi,D)
                dense=D.block(np.arange(n*n),np.arange(n*n))
                imag=(dense-dense.conj().T)/(2j)
                spectral=float(np.linalg.eigvalsh(imag)[0])
                margin=spectral-meta['self_shift']
                rows.append(dict(n=n,integrated=integrated,
                                 hermitian_min=spectral,lower_margin=margin,**meta))
                assert margin>=-2e-12
    D=FFTGreen(Geometry(n=8),3.5/(.4/8),cell_integrated=True)
    _,failure=gamma_for_fft_kernel(np.full(64,.5+.2j),D)
    assert not failure['applicable']
    changed=FFTGreen(Geometry(n=8),1/(.4/8),cell_integrated=True)
    changed.kernel[7,8]-=5j
    _,tamper=gamma_for_fft_kernel(np.full(64,.5+.2j),changed)
    assert not tamper['applicable']
    changed=FFTGreen(Geometry(n=8),1/(.4/8),cell_integrated=True)
    changed._fk[0,0]+=1j
    _,buffer=gamma_for_fft_kernel(np.full(64,.5+.2j),changed)
    assert not buffer['applicable']
    out=ROOT/'runs/a2/quadrature_certificate';out.mkdir(parents=True,exist_ok=True)
    (out/'results.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows,unavailable=failure,changed_kernel_rejected=tamper,changed_fft_rejected=buffer),indent=2))
    print(json.dumps(dict(checks=len(rows),minimum_margin=min(x['lower_margin'] for x in rows),outside_condition=failure),indent=2))
if __name__=='__main__':main()
