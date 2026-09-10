"""Registered physical Born / isolated-sphere / interacting-sphere controls."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scale_gauge_experiment import mx, profiled

HERE = Path(__file__).resolve().parent
CENTERS = np.array([[-.06, 0., 0.], [.055, .02, 0.]])
RADII = [.035, .025]
TEMPLATE = np.array([1+.03j, 2+.05j])
K = 18.
RX = mx.receivers(12, .6)


def field(scale, interaction):
    eps = 1+scale*TEMPLATE
    if interaction:
        return mx.treams_field(CENTERS, RADII, eps, K, RX, lmax=4)[0]
    return sum(mx.treams_field([c], [r], [e], K, RX, lmax=4)[0]
               for c, r, e in zip(CENTERS, RADII, eps))


def born(spacing):
    xyz, labels, fill = mx.voxelize(CENTERS, RADII, spacing, fill_quadrature=4)
    inc = np.stack([np.exp(1j*K*(xyz@d))[:, None]*p
                    for d, p in mx.illuminations()], axis=-1)
    currents = (spacing**3*fill*TEMPLATE[labels])[:, None, None]*inc
    return (mx.dipole_kernel(RX, xyz, K)@currents.reshape(3*len(xyz), -1)).reshape(12, 3, 4)


def main():
    target = HERE/'results/interaction_ablation.json'
    if target.exists():
        raise RuntimeError('Preserve completed evidence; register another version')
    start = time.perf_counter()
    paths = [Path(__file__), HERE/'ROUND2_PROTOCOL.md', Path(mx.__file__),
             HERE/'scale_gauge_experiment.py']
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    b0, b1 = born(.007), born(.005)
    quadrature = float(np.linalg.norm(b1-b0)/np.linalg.norm(b1))
    born_hidden = float(np.linalg.norm(profiled(b1, 1.25*b1, 'shared'))/np.linalg.norm(b1))
    rows = []
    for scale in (.01, .03, .1, .3, 1., 2.):
        row = {'scale': scale, 'methods': {}}
        central = {}
        for interaction, name in ((True, 'interacting'), (False, 'isolated_sum')):
            y = field(scale, interaction)
            central[name] = y
            ds = []
            for h in (1e-4, 5e-5):
                d = (field(scale*np.exp(h), interaction)-field(scale*np.exp(-h), interaction))/(2*h)
                ds.append(profiled(y, d, 'shared'))
            other = field(1.25*scale, interaction)
            row['methods'][name] = {
                'relative_sensitivity': float(np.linalg.norm(ds[0])/np.linalg.norm(y)),
                'absolute_sensitivity': float(np.linalg.norm(ds[0])),
                'finite_pair_residual': float(np.linalg.norm(profiled(other, y, 'shared'))/np.linalg.norm(y)),
                'derivative_step_change': float(np.linalg.norm(ds[1]-ds[0])/max(np.linalg.norm(ds[1]), 1e-30)),
                'relative_born_mismatch': float(np.linalg.norm(y-scale*b1)/np.linalg.norm(y))}
        row['relative_interaction_field_change'] = float(np.linalg.norm(central['interacting']-central['isolated_sum'])/np.linalg.norm(central['interacting']))
        rows.append(row)
    result = {'complete': True, 'source_hashes': hashes, 'rows': rows,
              'born_quadrature_change': quadrature, 'born_gain_hidden_residue': born_hidden,
              'wall_seconds': time.perf_counter()-start,
              'checks': {'born_scaling': born_hidden < 1e-12,
                         'derivative_steps': all(m['derivative_step_change'] < 1e-3 for r in rows for m in r['methods'].values())},
              'scope': 'Deterministic mechanism ablation, not recovery, solver certification or novelty.'}
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
