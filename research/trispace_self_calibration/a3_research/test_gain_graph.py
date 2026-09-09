"""Parent independent algebra audit, not a physical calibration benchmark."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from gain_graph import incidence, visible_geometry, complex_real_matrix, span, branch_bound


def run():
    rng = np.random.default_rng(5101)
    checks = []
    def check(name, value, threshold):
        checks.append(dict(name=name, value=float(value), threshold=float(threshold), passed=bool(value<=threshold)))
    edges = [(0, 0), (0, 1), (1, 0), (1, 1)]
    h = rng.normal(size=4)+1j*rng.normal(size=4)
    j = rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3))
    cov0 = rng.normal(size=(8, 8))
    cov = cov0@cov0.T+np.eye(8)
    result = visible_geometry(h, edges, j, covariance=cov)
    check('rectangle_has_one_complex_cycle', abs(result['gain_cycle_dimension_complex']-1), 0)
    check('one_cycle_at_most_two_real_geometry_directions', max(0,result['target_rank']-2), 0)
    for k in [1, 2, 3]:
        out = visible_geometry(h[:k], edges[:k], j[:k])
        check('forest_zero_information_'+str(k), np.linalg.norm(out['visible'])/np.linalg.norm(j[:k]), 1e-12)
    # Correctly transformed cycle covariance equals whitened gain projection.
    b = incidence(edges)
    z = null_space(b.T).T
    lraw = complex_real_matrix(z/h[None, :])
    chol = np.linalg.cholesky(cov)
    l = lraw@chol
    pcycle = l.T@np.linalg.solve(l@l.T, l)
    gain = np.linalg.solve(chol, complex_real_matrix(h[:, None]*b))
    q = span(gain)
    check('correlated_noise_cycle_projector', np.linalg.norm(pcycle-(np.eye(8)-q@q.T)), 1e-11)
    wrong = lraw.T@np.linalg.solve(lraw@lraw.T, lraw)
    check('negative_control_unwhitened_projection_disagrees', 1/np.linalg.norm(wrong-pcycle), 100)
    # Rank-one transfer geometry changes remain in the gain orbit.
    u = np.array([1+.2j, .7-.1j]); v=np.array([.8+.4j, 1.2-.1j])
    du=rng.normal(size=2)+1j*rng.normal(size=2)
    hv=(u[:,None]*v).ravel()
    dj=(du[:,None]*v).ravel()[:,None]
    out=visible_geometry(hv, edges, dj)
    check('rank_one_pose_hidden', np.linalg.norm(out['visible'])/np.linalg.norm(dj), 1e-12)
    # A new receiver needs a pair of edges to form a cycle.
    base=[(0,0),(0,1)]
    dims=[]
    for extra in [[],[(1,0)],[(1,1)],[(1,0),(1,1)]]:
        bb=incidence(base+extra)
        dims.append(len(base+extra)-np.linalg.matrix_rank(bb))
    check('pair_complementarity', np.linalg.norm(np.array(dims)-[0,0,0,1]), 0)
    # A3 mechanical translation / body phase-center gauge.
    from scipy.spatial.transform import Rotation
    rotations=[np.eye(3),Rotation.from_rotvec([0,0,.6]).as_matrix(),Rotation.from_rotvec([.5,0,0]).as_matrix()]
    two=np.vstack([np.c_[np.eye(3),r] for r in rotations[:2]])
    three=np.vstack([np.c_[np.eye(3),r] for r in rotations])
    check('two_3D_orientations_rank5', abs(np.linalg.matrix_rank(two)-5), 0)
    check('three_noncoaxial_orientations_rank6', abs(np.linalg.matrix_rank(three)-6), 0)
    r=1.3; delay=.07; g=.8+.2j; shift=.2; ks=np.array([2,5,12.])
    y=g*np.exp(1j*ks*(r+delay))/r
    yp=g*(r+shift)/r*np.exp(1j*ks*(r+shift+delay-shift))/(r+shift)
    check('all_frequency_range_clock_gain_gauge', np.linalg.norm(y-yp), 1e-12)
    bound=branch_bound([[0,0],[10,0]], beta=.1)
    check('conditional_not_global_certificate', int(bound['global_coverage_established']), 0)
    report={'scope':'parent algebra only; not final physical or SOM gate', 'seed':5101,'checks':checks,'passed':all(c['passed'] for c in checks)}
    path=Path(__file__).resolve().parent/'results'
    path.mkdir(exist_ok=True)
    (path/'parent_gain_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if not report['passed']:
        raise AssertionError('parent gain checks failed')


if __name__=='__main__':
    run()
