"""Bounded quadrature, grid-refinement, cylinder, and local-T port checks."""
from __future__ import annotations
import json, os, resource, sys, time
from pathlib import Path
import numpy as np
from scipy.special import jv, jvp, hankel1, h1vp
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.ports import GaussianComponent,LocalTNetwork,render,ownership_blocks
OUT=ROOT/'runs/a2/ports';OUT.mkdir(parents=True,exist_ok=True)
def rss_mb():
    x=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return x/(1024*1024) if sys.platform=='darwin' else x/1024

def components(separation):
    return [GaussianComponent(.32+.018j,(-separation,0.),(.045,.050),.15),
            GaussianComponent(.28+.015j,(separation,.012),(.050,.040),-.18),
            GaussianComponent(.22+.012j,(0.,.075),(.038,.048),.25)]

def scattered_rel(a,b): return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))

def cylinder_exact(k,rxy,radius,chi,mmax=50):
    """Scalar transmission-cylinder scattered field for e^(ik x) incidence."""
    ki=k*np.sqrt(1+chi); kr=k*radius; kir=ki*radius; out=np.zeros(len(rxy),complex)
    rr=np.linalg.norm(rxy,axis=1);phi=np.arctan2(rxy[:,1],rxy[:,0])
    for m in range(-mmax,mmax+1):
        inc=(1j)**m
        num=jv(m,kr)*ki*jvp(m,kir)-k*jvp(m,kr)*jv(m,kir)
        den=k*h1vp(m,kr)*jv(m,kir)-hankel1(m,kr)*ki*jvp(m,kir)
        out += inc*(num/den)*hankel1(m,k*rr)*np.exp(1j*m*phi)
    return out

def cylinder_control(n):
    g=Geometry(n=n,n_tx=1,n_rx=40);v=VIE(g,cell_integrated=True)
    radius=.075;chi0=.15+.01j; chi=chi0*(np.sum(v.points*v.points,axis=1)<=radius*radius)
    # Plane wave is used only for this analytic forward control.
    e=np.exp(1j*v.k*v.points[:,0])[:,None]; blocks=[np.flatnonzero(chi)]
    j,_=v.solve(chi,e,rtol=2e-8,maxiter=240,blocks=blocks);sca=v.S@j
    exact=cylinder_exact(v.k,v.rx,radius,chi0)
    return {"n":n,"scattered_relative_error_to_cylinder_series":scattered_rel(sca[:,0],exact),
            "scope":"pixelized disk plus finite-cell quadrature; this is a scalar 2-D convention/control, not a Maxwell validation."}

def main():
    result={"scope":"Bounded finite-grid checks. Cell-integrated is an independent quadrature option; results do not infer continuum convergence or a universal port rank.","threads":{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}}
    result['non_power_two_geometry_supported']=Geometry(n=96).n == 96
    # Same physical Gaussian scene, sampled on three grids, and point vs
    # cell-integrated quadrature.  N128 cell result is the numerical reference.
    conv=[]; cache={}
    cc=components(.065)
    for n in (32,64,128):
        for backend in ('point','cell_integrated'):
            v=VIE(Geometry(n=n,n_tx=4,n_rx=24),cell_integrated=(backend=='cell_integrated'))
            chi=render(cc,v.points); tic=time.perf_counter();f=v.forward(chi,rtol=2e-8,maxiter=240,blocks=ownership_blocks(cc,v.points))
            cache[n,backend]=f['scattered'];conv.append({"n":n,"backend":backend,"wall_seconds":time.perf_counter()-tic,"rss_peak_mb":rss_mb(),"gmres_iterations":f['solve']['iterations']})
    ref=cache[128,'cell_integrated']
    for row in conv: row['scattered_relative_to_n128_cell']=scattered_rel(cache[row['n'],row['backend']],ref)
    result['grid_convergence']=conv
    result['cylinder_analytic_control']=[cylinder_control(n) for n in (32,64,128)]
    # Three separation/overlap cases.  Exact component identity precedes port
    # truncation; every reported field error uses scattered field only.
    sweeps=[]
    for sep,label in ((.12,'separated'),(.065,'moderate_overlap'),(.025,'strong_overlap')):
        v=VIE(Geometry(n=16,n_tx=8,n_rx=32));cs=components(sep);chi=render(cs,v.points);full=v.dense_reference(chi);net=LocalTNetwork(v,cs)
        et=v.E+v.D.matmat(full['current']);cis=[xi[:,None]*et for xi in net.chis]
        item={"configuration":label,"center_separation":sep,"full_component_identity_relative_errors":net.full_component_identity_error(cis,v.E),"ports":[]}
        for rank in (3,6,12,24):
            p=net.solve_ports(rank,rtol=1e-9,enrich=2)
            item['ports'].append({"requested_ports_per_component":rank,"network_dimension":p['network_dimension'],"scattered_field_relative_error":scattered_rel(p['scattered'],full['scattered']),"current_relative_error":scattered_rel(p['current'],full['current'])})
        sweeps.append(item)
    result['local_t_response_port_sweep']=sweeps
    (OUT/'quadrature_ports_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
