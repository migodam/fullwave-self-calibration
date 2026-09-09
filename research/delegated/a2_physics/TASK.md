# Bounded worker: full-wave physical core

You are the DeepSeek Flash implementation worker. Do not launch another worker
or AI Scientist. Do not edit provider/model configuration. Do not print secrets,
store private reasoning, or run destructive commands. Use apply_patch for files.
Work only under `research/delegated/a2_physics/`. Read the workspace AGENTS.md.

Implement and test a reusable Python physical core, not the final research claims.
Read `Theory/Questions/A2_PRASC_SOM_VALIDATION_PROTOCOL.md` and theorem package
sections 2, 8 and 9. Read the corrected older core at
`experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/src/helmholtz.py`.
Use the existing venv Python at
`experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python`.

## Fixed physical specification

- 2D scalar outgoing Helmholtz, e^(-i omega t); domain [-0.5,0.5]^2 metres.
- Inverse N=16, truth N=32 (N=8 is allowed ONLY for smoke tests).
- Highest k=4*pi rad/m, lambda_min=0.5 m. Frequencies kmax*[.25,.5,.75,1].
- Green g=i/4 H0^(1)(kr), D and S include k^2 times cell integration. Correct
  equal-area disk self-cell integral includes the negative 1/k^2 endpoint.
- Nine real coefficients on unnormalized Gaussian basis functions with centers
  {-0.25,0,0.25} squared and width 0.16 m. epsilon_r-1=Phi@alpha.
- Ohmic law: sigma=0.005*(epsilon_r-1) S/m, so contrast=(Phi@alpha)*
  (1+i*0.005/(omega*epsilon0)). Fixed conductivity law, not independently
  estimated conductivity. c0=299792458, epsilon0=8.8541878128e-12.
- Nominal poses (px,py,theta): (-.15,-.12,0), (.15,0,.65), (0,.15,1.3).
- First pose anchored. The same unknown x=(dx,dy,dtheta) adds world translation
  and rotation angle to poses 2 and 3. It is a shared extrinsic error, NOT six
  independent unknowns. Receivers and transmitters co-move rigidly.
- 12 receiver body positions radius1.5; full aperture 12 angles 2*pi*m/12;
  limited aperture 12 inclusive angles linspace(-pi/4,pi/4,12). Two transmitter
  body positions radius1.9 at angles {-pi/3,+pi/3}.
- Point-source incident field g(z,tx). Raw data is TOTAL field at receivers:
  g(rx,tx)+S j. Retain background and scattered components separately.
- For fixed material/pose, share one M=I-diag(contrast)D factorization among all
  illuminations/poses at that frequency. Multi-RHS solves are charged per RHS,
  not once per batch. Factorizations are separately counted and timed.

## API to implement in physics.py

Use NumPy/SciPy complex128. Importing must perform no simulation or file writes.
Provide:
- `Config(N=16, aperture='full', kmax=4*np.pi)` dataclass.
- `Model(config)` with `.points`, `.basis` (n,9), `.ks`, `.pose_metric`
  diag(1,1,1.5**2), and `.forward(alpha,x,freq_ids=None,jacobian=False)`.
- forward returns a dict with `total` complex flat vector, `scattered`,
  `incident`, `A` complex (m,9) and `B` complex(m,3) if requested, `blocks`
  giving frequency/pose/illumination and row slices, `work` counts, `states`.
  Stable order frequency, pose, illumination, receiver. freq_ids default all4.
- Each state record includes M, D, S, j, b, Etotal, source/receiver positions,
  derivatives dS/dx, db/dx and material derivative factor or enough to build
  reduced state solves. Cache pose-independent D and frequency factors.
- A is total physical material derivative S M^-1 diag(Etotal) T. B includes
  receiver and transmitter effects. Direct incident at receiver derivative must
  be accounted for (it is exactly zero for rigid co-motion in this model).
- `realify(z,sigma=1)` -> sqrt(2)/sigma*[real;imag], and analogous real Jacobian.
- Support optional material_basis matrix override to permit pixel imaging later,
  while standard A2 core has nine real unknowns. Config/model validation should
  prevent incompatible sizes.

Prefer analytic batched parameter derivatives. Also offer optional efficient
adjoint gradient of 0.5||realify(total-y,sigma)||^2 if straightforward, with cost
counts including every forward/adjoint RHS. Never count evaluation truth solves
as solver access and never hide derivative/certificate costs.

## Required checks and outputs

Create `test_physics.py` executable (unittest acceptable), checks analytic A/B
versus centered FD at nonzero pose, both apertures, all four frequencies; raw
total=scattered+incident; state residual; outgoing convention; S/D k^2 factors;
noise real variance; gauge-related rigid distance invariance; factor reuse.
Run checks and N16/32 timing. Tests may write `checks.json`, `summary.md`, and
`environment.json` under this folder. Do not run seeds1001+ or invent results.
Return only success/failure, 3–10 facts, paths, uncertainties and intervention needs.
