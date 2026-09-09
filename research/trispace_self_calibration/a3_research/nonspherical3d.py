"""ADDA cross-implementation validation, with declared shared DDA approximation.

ADDA solves the reference currents (QMR, LDR), independently of our FFT/GMRES
CM+RR model. An RRC same-grid control tests units and the adapter, not physics
accuracy. Sphere/Treams and refinement controls test discretization separately.
All numerical lengths passed to ADDA share a common scale; its micrometre unit
is reinterpreted as metres with the same k*length and nondispersive epsilon.
"""
import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path
import numpy as np

from maxwell3d import dipole_kernel, illuminations, receivers, treams_field

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[1] / 'delegated' / 'a3_maxwell_refine'))
from tangent_fft import TangentFFTVIE
from maxwell_fft import FFTDyadicKernel, dyadic_kernel_values

ADDA = OUT / 'external/adda/src/seq/adda'
AXES = {'sphere': np.array([.16]*3),
        'ellipsoid': np.array([.16, .11, .075]),
        'box': np.array([.16, .11, .075])}


def shape_grid(shape, n, quadrature=1):
    axes = AXES[shape]
    h = 2*axes[0]/n
    dims = np.ceil(2*axes/h - 1e-10).astype(int)
    index = np.indices(dims).reshape(3, -1).T
    points = (index - (dims-1)/2)*h
    offsets = (np.arange(quadrature)+.5)/quadrature-.5
    fill = np.zeros(len(points))
    for off in np.stack(np.meshgrid(offsets, offsets, offsets, indexing='ij'), -1).reshape(-1, 3):
        pp = (points+off*h)/axes
        inside = (np.max(abs(pp), axis=1)<1) if shape=='box' else (np.sum(pp**2, axis=1)<1)
        fill += inside/quadrature**3
    take = fill>0
    return h, tuple(dims), points, index[take], fill[take]


class ShapeKernel(FFTDyadicKernel):
    def __init__(self, shape, n, k, quadrature=4):
        h, dims, full, idx, fill = shape_grid(shape, n, quadrature)
        self.spacing=h; self.k=k; self.rect_shape=dims
        self.pad_shape=tuple(2*v for v in dims)
        self.n_rect=int(np.prod(dims)); self.n_pad=int(np.prod(self.pad_shape))
        self.full_points=full; self.active_flat=np.ravel_multi_index(idx.T, dims)
        self.active_unravel=np.unravel_index(self.active_flat, dims)
        self.points=full[self.active_flat]; self.n_vox=len(fill)
        self.fill=fill; self.labels=np.zeros(len(fill), dtype=int)
        self.centers=np.array([[0.,0.,0.]])
        self.radii=np.array([max(AXES[shape])]) # bounding metadata, not shape geometry
        val=dyadic_kernel_values(dims, self.pad_shape, h, k)
        self.spectra=[[np.fft.fftn(val[...,i,j]) for j in range(3)] for i in range(3)]


class ShapeVIE(TangentFFTVIE):
    def __init__(self, shape, n, k, quadrature=4):
        h=2*AXES[shape][0]/n
        # Initialize solver bookkeeping on a tiny temporary sphere, then replace
        # the entire geometry/operator before any physical state is solved.
        super().__init__([[0,0,0]], [h], h, k, fill_quadrature=1,
                         gmres_rtol=1e-10, restart=80, maxiter=100,
                         wall_limit_seconds=120)
        self.kernel=ShapeKernel(shape, n, k, quadrature)
        self.ndof=3*self.kernel.n_vox
        self.incident=self._build_incident()


def read_vector(path):
    a=np.loadtxt(path, skiprows=1)
    return a[:,:3], a[:,4::2]+1j*a[:,5::2]


def independent_radiation(points, dipoles, rx, k):
    """CGS point-dipole expression, independently coded without parent kernel.

    p_parent=4*pi*P_ADDA. Returns scattered, not total, fields. This evaluator
    shares the analytic free-space Green function with the inverse model.
    """
    out=[]
    for r in rx:
        dr=r-points; d=np.linalg.norm(dr, axis=1); u=dr/d[:,None]
        dot=np.einsum('ni,nit->nt', u, dipoles)
        transverse=dipoles-u[:,:,None]*dot[:,None,:]
        longitudinal=3*u[:,:,None]*dot[:,None,:]-dipoles
        term=(k*k/d)[:,None,None]*transverse + (1/d**3-1j*k/d**2)[:,None,None]*longitudinal
        out.append(np.sum(np.exp(1j*k*d)[:,None,None]*term, axis=0))
    return np.array(out)


def reference(shape, n, k, epsilon, prescription='ldr'):
    """Cached independent executable runs; no silent use of incomplete output."""
    h,dims,full,idx,fill=shape_grid(shape,n,1)
    tag=f'{shape}_n{n}_k{k:g}_{epsilon.real:g}_{epsilon.imag:g}_{prescription}'
    folder=OUT/'results/adda'/tag
    folder.mkdir(parents=True,exist_ok=True)
    geom=folder/'shape.dat'
    if not geom.exists(): np.savetxt(geom,idx,fmt='%d',header='ADDA centered integer geometry')
    # ADDA recenters the occupied bounding box. Require this to be our intended
    # origin, rather than silently importing a displaced reference object.
    actual=(idx-(idx.min(0)+idx.max(0))/2)*h
    desired=full[np.ravel_multi_index(idx.T,dims)]
    assert np.max(abs(actual-desired))<1e-12
    m=np.sqrt(epsilon+0j)
    currents=[]; elapsed=0.; meta=[]
    for group,direction in enumerate([[0,0,1],[1,0,0]]):
        run=folder/f'prop{group}'
        command=[str(ADDA),'-shape','read',str(geom),'-lambda',str(2*np.pi/k),
                 '-dpl',str(2*np.pi/(k*h)),'-m',str(m.real),str(m.imag),
                 '-pol',prescription,'-int','poi','-no_vol_cor','-sym','no',
                 '-eps','10','-iter','qmr','-recalc_resid','-maxiter','500',
                 '-prop',*[str(v) for v in direction],'-store_dip_pol','-store_beam',
                 '-scat_matr','none','-dir',str(run)]
        stamp=folder/f'prop{group}_status.json'
        if not stamp.exists():
            start=time.perf_counter()
            try:
                done=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=180)
            except subprocess.TimeoutExpired as exc:
                stamp.write_text(json.dumps({'status':'timeout','command':command})+'\n')
                raise RuntimeError('ADDA reference timed out; failure preserved') from exc
            (folder/f'prop{group}_console.log').write_text(done.stdout)
            state=dict(status='ok' if done.returncode==0 else 'failed',returncode=done.returncode,
                       executable_sha256=hashlib.sha256(ADDA.read_bytes()).hexdigest(),
                       elapsed_seconds=time.perf_counter()-start,command=command)
            stamp.write_text(json.dumps(state,indent=2)+'\n')
        state=json.loads(stamp.read_text())
        if state['status']!='ok':raise RuntimeError(f'Unsuccessful ADDA reference {stamp}')
        residuals=[float(v) for v in re.findall(r'Final \(recalculated\) residual norm:\s*([\d.Ee+\-]+)', (run/'log').read_text())]
        if len(residuals)!=2 or not np.all(np.isfinite(residuals)) or max(residuals)>1e-9:
            raise RuntimeError(f'Independent ADDA true-residual audit failed: {run}')
        state['verified_true_residuals']=residuals
        elapsed+=state['elapsed_seconds']; meta.append(state)
        pp=[]; polarizations=[]
        for suffix in ['X','Y']:
            points,p=read_vector(run/f'DipPol-{suffix}')
            bpoints,beam=read_vector(run/f'IncBeam-{suffix}')
            assert np.max(abs(points-bpoints))<1e-12
            pol=beam*np.exp(-1j*k*(points@direction))[:,None]
            assert np.max(abs(pol-pol.mean(0)))<2e-8
            polarizations.append(pol.mean(0));pp.append(p)
        basis=np.array(polarizations).T
        target=np.array([illuminations()[2*group+i][1] for i in range(2)]).T
        coef=np.linalg.lstsq(basis,target,rcond=None)[0]
        assert np.linalg.norm(basis@coef-target)<1e-8
        currents.append(np.stack(pp,axis=-1)@coef)
    return points,np.concatenate(currents,axis=-1),dict(
        voxels=len(points),spacing=h,reference_seconds=elapsed,folder=str(folder),runs=meta)


def controls():
    k=9.;eps=2.5+.03j;rx=receivers(13,1.3)
    pts,p,meta=reference('sphere',8,k,eps,'rrc')
    ref=independent_radiation(pts,p,rx,k)
    via_parent=(dipole_kernel(rx,pts,k)@(4*np.pi*p).reshape(3*len(pts),4)).reshape(len(rx),3,4)
    radiation_error=float(np.linalg.norm(ref-via_parent)/np.linalg.norm(ref))
    model=ShapeVIE('sphere',8,k,quadrature=1)
    inv=model.field([eps],rx)
    same_error=float(np.linalg.norm(ref-inv)/np.linalg.norm(ref))
    exact=treams_field([[0,0,0]],[.16],[eps],k,rx,lmax=7)[0]
    sphere_error=float(np.linalg.norm(ref-exact)/np.linalg.norm(exact))
    # Exact state convention control is separate from coarse continuum error.
    assert radiation_error<1e-12 and same_error<1e-7
    p0,dp=model.currents([eps],derivatives=True)
    step=1e-5
    fd=(model.currents([eps+step])-model.currents([eps-step]))/(2*step)
    tangent_error=float(np.linalg.norm(fd-dp[:,:,0])/np.linalg.norm(fd))
    assert tangent_error<1e-7
    result=dict(radiation_adapter_error=radiation_error,same_grid_rrc_error=same_error,
                sphere_coarse_treams_error=sphere_error,material_tangent_fd_error=tangent_error,
                passed=True,scope='adapter/derivative controls; continuum discrepancy not a pass criterion')
    (OUT/'results/nonspherical_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


def study(shape,ks,ns):
    destination=OUT/'results'/f'nonspherical_{shape}.json'
    rows=json.loads(destination.read_text()) if destination.exists() else []
    eps=2.5+.03j;rx=receivers(17,1.3)
    for k in ks:
        previous=None
        for n in ns:
            pts,p,meta=reference(shape,n,k,eps)
            ref=independent_radiation(pts,p,rx,k)
            record=next((r for r in rows if r['n']==n and r['k']==k),None)
            if record is None:
                start=time.perf_counter();model=ShapeVIE(shape,n,k)
                inv=model.field([eps],rx)
                record=dict(shape=shape,n=n,k=k,epsilon=[eps.real,eps.imag],
                    relative_field_disagreement=float(np.linalg.norm(inv-ref)/np.linalg.norm(ref)),
                    reference_refinement_change=(None if previous is None else float(np.linalg.norm(ref-previous)/np.linalg.norm(ref))),
                    inverse_voxels=model.kernel.n_vox,reference_voxels=len(pts),spacing=meta['spacing'],
                    inverse_seconds=time.perf_counter()-start,reference_seconds=meta['reference_seconds'],
                    max_true_iterative_residual=max(model.true_relative_residuals),
                    memory_estimate_mib=model.memory_estimate()['estimated_peak_bytes']/2**20,
                    scope='development; independent ADDA LDR vs fill-weighted CM+RR; shared DDA family',status='ok')
                if shape=='sphere':
                    exact=treams_field([[0,0,0]],[.16],[eps],k,rx,lmax=9)[0]
                    record['reference_treams_error']=float(np.linalg.norm(ref-exact)/np.linalg.norm(exact))
                    record['inverse_treams_error']=float(np.linalg.norm(inv-exact)/np.linalg.norm(exact))
                rows.append(record)
                destination.write_text(json.dumps(rows,indent=2)+'\n')
                np.savez_compressed(OUT/'results/adda'/f'{shape}_n{n}_k{k:g}_fields.npz',rx=rx,reference=ref,inverse=inv)
            print(json.dumps(record),flush=True);previous=ref
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true')
    parser.add_argument('--shape',choices=AXES,default='ellipsoid')
    parser.add_argument('--ns',type=int,nargs='+',default=[12,20,32])
    parser.add_argument('--ks',type=float,nargs='+',default=[9.,18.])
    args=parser.parse_args()
    if args.controls: controls()
    else: study(args.shape,args.ks,args.ns)
