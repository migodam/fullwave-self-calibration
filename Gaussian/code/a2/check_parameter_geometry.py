"""Small exact linear-algebra checks for the A2 parameter SOM audit.

These do not test nonlinear inverse convergence or electromagnetic novelty.
"""
from pathlib import Path
import json, hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]

def main():
    rng = np.random.default_rng(2026091102)
    records = []
    for _ in range(30):
        A = rng.normal(size=(31, 12)); r = rng.normal(size=31)
        V = np.linalg.qr(rng.normal(size=(12, 7)))[0]
        lam = .1; H = A.T @ A + lam*np.eye(12); g = A.T @ r
        q = np.linalg.solve(H, g)
        qv = V @ np.linalg.solve(V.T @ H @ V, V.T @ g)
        model = lambda z: .5*np.linalg.norm(A@z-r)**2 + .5*lam*np.linalg.norm(z)**2
        gap = model(qv)-model(q); identity = .5*(qv-q)@H@(qv-q)
        residual = g-H@qv; bound = np.linalg.norm(residual)**2/(2*lam)
        # Parameter units transformed together with their physical metric.
        units = np.diag(np.exp(rng.uniform(-6,6,12)))
        P = np.diag(np.exp(rng.uniform(-2,2,12)))
        physical_A = A@P
        transformed_A = A@np.linalg.inv(units)@(units@P)
        records.append(dict(gap=float(gap),identity_error=float(abs(gap-identity)),
                            bound=float(bound),covered=bool(gap<=bound+1e-10),
                            coordinate_error=float(np.linalg.norm(physical_A-transformed_A))))
    # A decisive counterexample to using projected-gradient coverage as an
    # optimal-step certificate. H is SPD and g lies exactly in span(e1).
    H = np.array([[1.,.9],[.9,1.]])
    g = np.array([1.,0.]); V = np.array([[1.],[0.]])
    q = np.linalg.solve(H,g); qv = V@np.linalg.solve(V.T@H@V,V.T@g)
    deficit = 1.-np.linalg.norm(V.T@g)**2/np.linalg.norm(g)**2
    counter = dict(gradient_deficit=float(deficit),full_step=q.tolist(),
                   restricted_step=qv.tolist(),model_gap=float(.5*(q-qv)@H@(q-qv)))
    out = dict(scope='Finite dimensional algebra checks, no full-wave inverse acceptance',
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               rows=records,gradient_rule_counterexample=counter)
    target=ROOT/'runs/a2/parameter_geometry'; target.mkdir(parents=True,exist_ok=True)
    (target/'results.json').write_text(json.dumps(out,indent=2))
    assert all(x['covered'] and x['identity_error']<1e-10 and x['coordinate_error']<1e-10 for x in records)
    assert deficit == 0 and counter['model_gap'] > 2
    print(json.dumps({'checks':len(records),'counterexample':counter},indent=2))

if __name__ == '__main__': main()
