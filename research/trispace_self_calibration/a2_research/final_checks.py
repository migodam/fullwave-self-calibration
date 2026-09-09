"""Independent parent checks; never changes frozen optimization records."""
from pathlib import Path
import importlib.util
import json
import hashlib
import numpy as np
from scipy.special import hankel1, jv

ROOT = Path(__file__).resolve().parent
DELEGATED = ROOT.parents[1] / "delegated"

def efficient(J):
    # Square-root residualization avoids cancellation in an inverse Schur form.
    w, v = np.linalg.eigh(np.asarray(J))
    assert min(w) > -1e-12
    X = np.sqrt(np.maximum(w, 0))[:, None] * v.T
    A, B = X[:, :9], X[:, 9:]
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    U = U[:, s > 1e-12 * s[0]]
    V = B - U @ (U.T @ B)
    return V.T @ V

def fisher_check():
    source = DELEGATED / "a2_theory_validation/replication/results/e1_physical_results.json"
    records = json.loads(source.read_text())["per_seed"]
    out = []
    for r in records:
        jc, jp = efficient(r["J_coh"]), efficient(r["J_ph"])
        ev = np.linalg.eigvalsh(jc-jp)
        assert ev[0] >= -1e-10 * max(1, np.linalg.norm(jc))
        out.append(dict(seed=r["seed"], min_eigenvalue_difference=float(ev[0]),
                        trace_coherent=float(np.trace(jc)), trace_intensity=float(np.trace(jp)),
                        eigenvalues_coherent=np.linalg.eigvalsh(jc).tolist(),
                        eigenvalues_intensity=np.linalg.eigvalsh(jp).tolist()))
    return dict(status="pass", source=str(source), scope="cached independently quadrature-evaluated physical Fisher matrices; material eliminated by square-root projection; not new scenes", cases=out)

def measured_sign_check():
    raw=DELEGATED/"a2_literature/data/2001_iop_17_6_301/dielTM_dec8f.exp"
    assert hashlib.sha256(raw.read_bytes()).hexdigest()=="476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb"
    arr=np.loadtxt(raw)
    assert arr.shape==(14112,7) and np.isfinite(arr).all()
    assert len(np.unique(arr[:,:3],axis=0))==14112
    assert set(arr[:,0])==set(range(1,37)) and set(arr[:,2])==set(range(1,9))
    for v in range(1,37):
        for f in range(1,9):
            assert np.count_nonzero((arr[:,0]==v)&(arr[:,2]==f))==49
    path = DELEGATED / "a2_measured/model.py"
    spec = importlib.util.spec_from_file_location("audited_mie_worker", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    k = 2*np.pi*8e9/module.C0
    nm = int(np.ceil(k*module.A_RAD))+15
    n = np.arange(-nm, nm+1)
    theta = np.arange(8192)*2*np.pi/8192
    angular = np.exp(1j*n[:, None]*theta[None, :])
    z = k*.2
    j, h = jv(n, z), hankel1(n, z)
    dj = .5*(jv(n-1,z)-jv(n+1,z))*k
    dh = .5*(hankel1(n-1,z)-hankel1(n+1,z))*k
    # e^{-i omega t}: outward power is Im(E* d_r E)/(2 omega mu).
    # The common positive factor 1/(2 omega mu) is omitted for sign checks.
    outgoing = float(np.imag(np.conj(hankel1(0,z))*(-k*hankel1(1,z))))
    assert outgoing > 0
    flux = {}
    w = (1j)**n  # plane wave travelling in +x, consistent with outgoing H1
    for key, eps in [("lossless", 3+0j), ("passive_positive_imag", 3+.2j), ("active_negative_imag", 3-.2j)]:
        c = module.mie_coefficients(k, eps, module.A_RAD, nm)
        E = ((w*(j+c*h))[:, None]*angular).sum(axis=0)
        Er = ((w*(dj+c*dh))[:, None]*angular).sum(axis=0)
        flux[key] = float(2*np.pi*.2*np.mean(np.imag(np.conj(E)*Er)))
    assert abs(flux["lossless"]) < 1e-10
    assert flux["passive_positive_imag"] < 0 < flux["active_negative_imag"]
    # Old incident weights encode -k u, whereas H1 source phase has +k R.
    u = np.array([.6,.8]); phi=np.arctan2(u[1],u[0])
    coords=.007*np.array([np.cos(theta),np.sin(theta)])
    basis=jv(n,k*.007)[:,None]*angular
    old=module.jacobi_anger(k,u,nm)@basis
    correct=((1j)**n*np.exp(-1j*n*phi))@basis
    expected=np.exp(1j*k*(u@coords))
    olderr=float(np.max(np.abs(old-expected)))
    newerr=float(np.max(np.abs(correct-expected)))
    assert olderr>.1 and newerr<1e-12
    return dict(status="worker_fit_rejected", ingestion="parent raw checksum/shape/finiteness/unique full coverage passed", convention="H1 outgoing, e^-iwt, passive Im(epsilon)>0", outgoing_H1_flux=outgoing,
                integrated_flux_without_positive_factor=flux,
                old_incident_direction_error=olderr, corrected_incident_direction_error=newerr,
                holdout_audit="run_fits calls gains_and_residual separately on test views, profiling gains again; not held-out prediction",
                parameter_audit="c_x,c_y are object coordinates, not antenna offset ground truth")

def main():
    result = {"efficient_fisher":fisher_check(), "measured_model_audit":measured_sign_check()}
    rows=json.loads((ROOT/"results/passivity_checks.json").read_text())["cases"]
    table=[]
    for k in sorted(set(r["k"] for r in rows)):
        sub=[r for r in rows if r["k"]==k and r["available"]]
        table.append(dict(k=k, state_ratio_range=[min(r["state_ratio"] for r in sub),max(r["state_ratio"] for r in sub)], derivative_ratio_range=[min(r["derivative_ratio"] for r in sub),max(r["derivative_ratio"] for r in sub)]))
    result["certificate_tightness"]=table
    (ROOT/"results/final_parent_checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__": main()
