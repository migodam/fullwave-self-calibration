"""Bounded verification and timing probe.  Run with BLAS threads <=2."""
from __future__ import annotations
import json, os, resource, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.ports import GaussianComponent,LocalTNetwork,render,ownership_blocks,tail_audit,pack,unpack,parameter_scales

OUT=ROOT/'runs/a2/ports'; OUT.mkdir(parents=True,exist_ok=True)
def rss_mb():
    # macOS ru_maxrss is bytes; Linux is KiB.
    raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw/(1024*1024) if sys.platform == 'darwin' else raw/1024
def main():
    comps=[GaussianComponent(.45+.02j,(-.07,0.),(.045,.055),.2),GaussianComponent(.34+.02j,(.07,.015),(.05,.04),-.3)]
    result={"scope":"Finite-grid scalar VIE numerical checks. Gaussian tails remain nonzero; response ports are not automatically power-S parameters.","blas_thread_limit":{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}}
    # Dense comparison is intentionally small and checks FFT forward plus adjoint.
    v=VIE(Geometry(n=16));chi=render(comps,v.points); blocks=ownership_blocks(comps,v.points)
    f=v.forward(chi,rtol=1e-10,blocks=blocks); d=v.dense_reference(chi)
    rng=np.random.default_rng(7);x=rng.normal(size=256)+1j*rng.normal(size=256);y=rng.normal(size=256)+1j*rng.normal(size=256)
    tangent=v.material_tangent(chi,.01*render([comps[0]],v.points),f,rtol=1e-10,blocks=blocks)
    h=2e-5; fp=v.forward(chi+h*.01*render([comps[0]],v.points),rtol=1e-10,blocks=blocks)['total'];fm=v.forward(chi-h*.01*render([comps[0]],v.points),rtol=1e-10,blocks=blocks)['total']
    result['n16_dense']={"forward_current_relative_error":float(np.linalg.norm(f['current']-d['current'])/np.linalg.norm(d['current'])),"forward_field_relative_error":float(np.linalg.norm(f['total']-d['total'])/np.linalg.norm(d['total'])),"adjoint_inner_product_relative_error":float(abs(np.vdot(v.D.matvec(x),y)-np.vdot(x,v.D.rmatvec(y)))/max(abs(np.vdot(v.D.matvec(x),y)),1e-30)),"tangent_fd_relative_error":float(np.linalg.norm((fp-fm)/(2*h)-tangent['total'])/np.linalg.norm(tangent['total'])),"gmres_iterations":f['solve']['iterations']}
    net=LocalTNetwork(v,comps); p=net.solve_ports(rank=3,rtol=1e-9,enrich=1)
    result['n16_response_ports']={"network_dimension":p['network_dimension'],"field_relative_error_vs_full":float(np.linalg.norm(p['total']-f['total'])/np.linalg.norm(f['total'])),"current_relative_error_vs_full":float(np.linalg.norm(p['current']-f['current'])/np.linalg.norm(f['current']))}
    result['parameterization']={"pack_roundtrip_max":float(np.max(abs(pack(unpack(pack(comps)))-pack(comps)))),"scales_positive":bool(np.all(parameter_scales(comps)>0))}
    result['tail_audit']=tail_audit(comps,v.points,v.geometry.side)
    timing=[]
    for n in (64,128):
        vv=VIE(Geometry(n=n,n_tx=4,n_rx=12)); cc=render(comps,vv.points); bb=ownership_blocks(comps,vv.points); tic=time.perf_counter(); ff=vv.forward(cc,rtol=2e-7,maxiter=180,blocks=bb)
        timing.append({"n":n,"unknowns":n*n,"wall_seconds":time.perf_counter()-tic,"rss_peak_mb":rss_mb(),"gmres_iterations":ff['solve']['iterations'],"operator_storage":"FFT kernel and work arrays; no dense D"})
    result['scalability_timing']=timing
    (OUT/'results.json').write_text(json.dumps(result,indent=2));(OUT/'README.md').write_text('# A2 Gaussian component/port probe\n\n`results.json` records finite-grid FFT-VIE checks, response-derived local-T port error and N64/N128 timings. It is a bounded numerical infrastructure check, not a low-rank theorem, continuum convergence result, or power-scattering-port claim.\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
