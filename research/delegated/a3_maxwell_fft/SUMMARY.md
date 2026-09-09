# A3 matrix-free FFT dyadic convolution + GMRES VIE worker

Status: bounded numerical development **success** (not a production or
scientific acceptance gate). No recursive delegation, no pipeline launch, no
installation, and no modification of parent code. All work lives in
`research/delegated/a3_maxwell_fft/`.

## What was implemented

`maxwell_fft.py` provides a reusable matrix-free route for the *same* occupied
regular-grid dipole VIE as the read-only parent
`research/trispace_self_calibration/a3_research/maxwell3d.py`:

- occupied cells are embedded in the parent bounding rectangular grid and the
  dyadic kernel is applied by zero-padded FFT (padding factor 2 per axis);
- component ordering, `exp(-i omega t)`, polarizability, and
  `fill_quadrature=6` subcell fill are inherited exactly from the parent;
- physical currents solve `(I - diag(alpha) K) p = alpha E_inc` for all four
  illuminations by restarted GMRES, i.e. `r_free = 0` (no added current
  nuisance freedom);
- every solve records matvecs, iterations, GMRES flags, wall time, grid and
  padding shape, true residual (independent post-solve check) and a
  documented conservative peak-memory estimate;
- a short independent Treams sphere/cluster convergence audit is provided at
  k = 9 and 18 with grid spacings .03, .02, .015.

## Verification evidence (small grids, against the dense parent)

All 24 checks in `checks.json` pass.

- FFT full-grid dyadic action vs dense parent kernel: relative error
  3e-16-9e-16 for origin, off-origin and parent pair grids.
- Occupied-grid FFT action vs dense occupied kernel: same level.
- Full-grid and occupied-grid reciprocity `u^T K v = v^T K u`: errors
  < 2e-15.
- Exact voxel embedding, labels, fill, polarizability and incident-field match
  with `DipoleVIE` (incident match exactly 0.0).
- GMRES currents vs dense LU currents: relative error 3e-11-8e-11;
  receiver-field relative error < 2e-11.
- Parent fill-study k=9, h=.03 field errors are reproduced to absolute
  difference < 2e-12 (single and pair), confirming that the FFT/GMRES route
  returns the parent VIE solution, not a modified model.

## Treams convergence audit (k=9: lmax5 vs 8; k=18: lmax8 vs 10)

Raw per-case data (all pointwise complex field errors, residuals, resolutions,
grid/padding and reference cross-sections) are preserved in
`treams_audit_raw.json`.

| scene | k | spacing | voxels | rect -> padded grid | dof | iterations (4 ill.) | solve wall s | est. peak MiB | true rel. residual | field rel. error |
|---|---|---|---|---|---|---|---|---|---|---|
| single | 9 | .030 | 408 | 8^3 -> 16^3 | 1224 | 16,16,16,16 | 0.024 | 2.5 | 6.3e-11 | 0.02076 |
| single | 9 | .020 | 1208 | 12^3 -> 24^3 | 3624 | 17,17,17,17 | 0.054 | 7.9 | 3.9e-11 | 0.01513 |
| single | 9 | .015 | 2632 | 16^3 -> 32^3 | 7896 | 17,17,17,17 | 0.154 | 17.8 | 4.5e-11 | 0.01268 |
| single | 18 | .030 | 408 | 8^3 -> 16^3 | 1224 | 21,21,21,21 | 0.031 | 2.5 | 4.9e-11 | 0.12565 |
| single | 18 | .020 | 1208 | 12^3 -> 24^3 | 3624 | 22,22,22,22 | 0.070 | 7.9 | 5.7e-11 | 0.06859 |
| single | 18 | .015 | 2632 | 16^3 -> 32^3 | 7896 | 22,22,22,22 | 0.197 | 17.8 | 9.5e-11 | 0.04487 |
| pair | 9 | .030 | 652 | [22,9,9] -> [44,18,18] | 1956 | 19,19,19,19 | 0.067 | 5.9 | 5.5e-11 | 0.02051 |
| pair | 9 | .020 | 1936 | [32,13,14] -> [64,26,28] | 5808 | 20,20,20,20 | 0.198 | 18.4 | 4.4e-11 | 0.01678 |
| pair | 9 | .015 | 4266 | [43,17,18] -> [86,34,36] | 12798 | 20,20,20,20 | 0.788 | 41.1 | 6.2e-11 | 0.01475 |
| pair | 18 | .030 | 652 | [22,9,9] -> [44,18,18] | 1956 | 34,34,33,33 | 0.117 | 5.9 | 9.0e-11 | 0.24839 |
| pair | 18 | .020 | 1936 | [32,13,14] -> [64,26,28] | 5808 | 37,37,36,36 | 0.357 | 18.4 | 8.3e-11 | 0.14325 |
| pair | 18 | .015 | 4266 | [43,17,18] -> [86,34,36] | 12798 | 38,38,37,37 | 1.428 | 41.1 | 8.6e-11 | 0.09452 |

## Findings

1. **FFT is an exact replacement for the dense dyadic operator**: verified to
   ~1e-15 on full and occupied grids, including off-origin embedding and
   reciprocity. Any remaining difference is roundoff, not padding aliasing.
2. **GMRES solves the physical current state**: true post-solve residuals are
   ~1e-10-1e-11, orders of magnitude below the 1e-2-0.25 discretization
   errors. GMRES tolerance does **not** remove discretization bias.
3. **Field error decreases monotonically with refinement in every audit
   scene/frequency** down to h=.015 at k=18, so this bounded audit does not
   reveal a non-converging finite-cell/polarizability failure. Errors are
   still substantial at high frequency on the finest audited grid (single
   4.5%, pair 9.5% at k=18), i.e. this is a numerical route to make such
   refinement affordable, not an accuracy certificate.
4. **Audit cost stayed tiny**: the hardest case (pair, k=18, h=.015) solved in
   about 1.4 s with estimated peak memory 41 MiB, so the 120 s / 4 GiB guards
   were not needed; all 12 cases completed (`status: ok`).
5. **Treams references are converged at the requested lmax pairs**: reference
   relative convergence is 7e-9-3e-8 at k=9 (lmax5 vs 8) and 1e-10-3e-10 at
   k=18 (lmax8 vs 10); low/high-lmax VIE errors agree to the same order.
6. No inverse sweep, formal runtime comparison, or physical-constant tuning
   was performed; parent hashes are recorded in `environment.json`.

## Artifacts and reproduction

- `maxwell_fft.py` - reusable matrix-free FFT dyadic kernel + GMRES solver.
- `test_maxwell_fft.py` - dense-oracle test suite (writes `checks.json`).
- `run_maxwell_fft.py` - validation + full Treams audit reproduction.
- `validation_raw.json`, `treams_audit_raw.json` - raw JSON evidence.
- `checks.json` - 24 check results (all passed).
- `environment.json` - interpreter/platform/package versions, command, and
  SHA-256 digests of the four parent read-only files.

Reproduce with:

```bash
research/trispace_self_calibration/a3_research/.venv3d/bin/python \
  research/delegated/a3_maxwell_fft/run_maxwell_fft.py
```

(the script reruns the checks and refuses to audit if any check fails; tests
can also be run directly with `test_maxwell_fft.py`).

## Uncertainties and scope

- Passes are bounded development checks of the implemented discrete model;
  they do not establish physical/continuum accuracy, journal acceptance, or
  that the remaining k=18 pair error has a finite-cell-only explanation.
- The memory estimate is a documented conservative formula, not an OS-measured
  peak-resident RSS.
- Reference quality is Treams multipole-cluster quality at the stated lmax
  pairs; it is an independent numerical reference, not measurement data.
