# Family-1 self-cell correction and revalidation report

Date: 2026-09-03 (UTC runs 09:26:02-09:26:17; SGT ~17:26).
Experiment: `experiment_pose_confounding_spectral_geometry`.
Audit gate: `context/PARENT_CORRECTIONS.md` (2026-09-03) and workshop rule 16
(`context/workshop.md`).

## 1. What the parent audit required

For the implemented 2D scalar Helmholtz Green function
`g(r,r') = (i/4) H_0^(1)(k_b |r-r'|)` with equal-area disk radius
`a = h/sqrt(pi)`, the domain-propagator diagonal must use the complete disk
integral

```
I_self = (i*pi*a/(2*k_b)) * H_1^(1)(k_b*a) - 1/k_b^2
[G_D]_nn = k_b^2 * I_self
```

The previous implementation omitted the `-1/k_b^2` lower-endpoint term
(`lim_{x->0} x H_1^(1)(x) = -2i/pi` under `d[x H_1(x)]/dx = x H_0(x)`), so its
alleged integral tended to `1/k_b^2` instead of zero as `h -> 0`. The
source-position gradient `grad_s g = -grad_z g` had already been corrected and
had to be preserved.

## 2. Invalidation archive (evidence preserved)

The pre-correction snapshot is preserved verbatim under
`archive/INVALIDATED_self_cell_2026-09-03T092339Z/` with a README and
`checksums.txt`:

- `src/helmholtz.py`, `src/family1_pilot.py`
- `results/family1_pilot_results.json`, `figures/family1_fd_convergence.png`
- `logs/environment.md` and the prior `notes/*` dumps

Archive time: 2026-09-03T09:23:39Z. Live file digests at archive time are
listed in the archive README. The old log remains in the live tree as
historical evidence but is superseded by
`logs/environment_self_cell_corrected.md`.

## 3. Source corrections (live tree)

`src/helmholtz.py`:

- `self_cell_green` now returns the complete equal-area disk integral
  `(i*pi*a/(2*k_b))*H_1^(1)(k_b*a) - 1/k_b^2`, with the derivation and the
  lower-endpoint limit documented.
- Module and function documentation updated (`SELF_CELL_FORMULA`,
  `SELF_CELL_FORMULA_VERSION = "self-cell v2 corrected 2026-09-03 ..."`).
- New exported helper `green_grad_source(z, s, k_b)` implements
  `grad_s g = +(i k_b/4) H_1^(1)(k_b R)(z - s)/R = -grad_z g`, and
  `build_AB` now calls it. A numerical check confirmed the refactor is
  bit-for-bit identical to the previously validated `-green_grad_first(...)`
  expression used for `grad2` (max abs diff 0.0), so the corrected sign was
  preserved.

`src/family1_pilot.py`: result JSON now stamps `self_cell_formula`,
`self_cell_formula_version`, the parent-correction gate, the archive path, the
exact command, runtime, platform/package versions, and source SHA-256 digests;
figure title identifies the corrected self-cell version.

## 4. Unit tests (new)

`src/test_family1_corrections.py`, run with
`.venv/bin/python -m unittest discover -s src -p "test_*.py" -v`
(5/5 pass, 5.43 s wall):

1. `self_cell_green` equals direct Gauss-Legendre radial disk quadrature
   (`abs err < 1e-11` enforced) for N in {8,12,16,24,32,40,64,128};
2. complete `I_self` tends to zero under refinement (monotone magnitude drop,
   `>=50x` reduction, `max |I_self| < 1e-2`; the old rule fails because it
   stays near `1/k_b^2 ~= 0.025`);
3. `green_grad_source` matches a centred finite difference of `g(z,s)` in the
   source coordinate `s` (rel err < 1e-5 enforced at `delta = 1e-6`);
4. `grad_s = -grad_z` algebraically for the helpers;
5. full-map pose FD against `B` (rel err < 1e-4 enforced at
   `delta = 1e-4`, seed 54321), guarding sign integration in `build_AB`.

Sign-test discrimination was checked by a scratch mutation of
`green_grad_source` to the old sign: worst relative error 2.000e+00 against
the finite difference versus 1.737e-09 for the corrected helper, so the old
sign fails the unit test.

## 5. Radial-quadrature validation of I_self

Runner: `src/validate_self_cell.py`; raw data:
`results/self_cell_quadrature_validation.json`; figure:
`figures/self_cell_quadrature_validation.png`. Method: 2048-node
Gauss-Legendre quadrature of `(i/4)*2*pi*r*H_0^(1)(k_b r)` over `[0, a]`,
k_b = 2 pi.

| N | h | abs I_self (analytic) | abs old invalid rule | abs quadrature error | abs [G_D]_nn |
|---:|---:|---:|---:|---:|---:|
| 8  | 0.125      | 5.123e-03 | 2.901e-02 | 3.5e-16 | 2.023e-01 |
| 16 | 0.0625     | 1.630e-03 | 2.666e-02 | 7.7e-17 | 6.433e-02 |
| 32 | 0.03125    | 5.003e-04 | 2.577e-02 | 1.7e-17 | 1.975e-02 |
| 40 | 0.025      | 3.399e-04 | 2.563e-02 | 2.3e-17 | 1.342e-02 |
| 64 | 0.015625   | 1.493e-04 | 2.547e-02 | 2.0e-18 | 5.895e-03 |
| 128| 0.0078125  | 4.358e-05 | 2.537e-02 | 1.0e-17 | 1.720e-03 |

Full ladder also includes N=12 and N=24 in the JSON. Max quadrature abs error
over the ladder: 3.13e-16. `tends_to_zero_under_refinement = true`:
`|I_self|` drops monotonically from 5.12e-3 (N=8) to 4.36e-5 (N=128) with a
first/last ratio of 117.6, while the invalidated rule remains at ~2.5e-2 (its
`1/k_b^2` limit). This is a finite-ladder diagnostic, not a continuum proof.

## 6. Family-1 pilot rerun from scratch

Runner: `src/family1_pilot.py` on the corrected harness; raw data:
`results/family1_pilot_results.json`; figure:
`figures/family1_fd_convergence.png`. Runtime: 2.13 s (recorded in JSON).

Selected results (identical seeds and config as the archived run):

- Max relative map-Jacobian FD error at FD step 1e-3: 2.25e-09 (seed 0),
  1.74e-08 (seed 1), 1.77e-09 (seed 2).
- Max relative pose-Jacobian FD error at FD step 1e-3: 8.38e-07 (seed 10),
  4.28e-06 (seed 11), 4.09e-06 (seed 12).
- Convergence slopes (5-point O(h^2) window): map seed 0 slope 1.976
  (R2=0.9998); pose seed 10 slope 2.0000 (R2=1.0000).
- State: `sigma_min(M)/||M||_2 = 0.7189`, `sigma_min(M) = 0.8329`, max state
  residual 9.34e-16, max `|E_tot| = 8.67e-02`.
- Born validity at s=1: `||D_chi G_D||_2 = 0.326`; Born/full-wave relative
  discrepancies at s = 0.01..0.2 range 1.77e-03..3.54e-02 (F) and
  1.89e-03..3.82e-02 (A Frobenius).

All raw per-seed curves, slopes, and the Born table are in the JSON.

## 7. Forward-model grid-refinement diagnostic

Runner: `src/family1_grid_refinement.py`; raw data:
`results/family1_grid_refinement.json`; figure:
`figures/family1_grid_refinement.png`. Same continuous two-Gaussian contrast
scene, same poses/receivers/k_b, grids N in {16,24,32,40} on Apple Silicon
CPU. Data dimension is fixed (24 complex measurements at shared receivers).

| N | h | rel diff vs previous | rel diff vs finest (N=40) | sigma_min(M)/||M|| |
|---:|---:|---:|---:|---:|
| 16 | 0.0625   | -          | 1.453e-03 | 0.71889 |
| 24 | 0.041667 | 9.505e-04  | 5.030e-04 | 0.71867 |
| 32 | 0.03125  | 3.430e-04  | 1.601e-04 | 0.71858 |
| 40 | 0.025    | 1.601e-04  | 0.0 (reference) | 0.71854 |

Status of the plan's `32x32` vs `40x40` sensitivity: **both completed, not
deferred** (N=32 1.18 s, N=40 3.71 s in this run). The successive differences
shrink monotonically over this ladder. Per the rules, this is labelled a
*forward-model refinement diagnostic*, not continuum convergence: the grids
are not a nested restriction of one continuum problem, and no continuum
transfer claim is made.

## 8. Environment, commands, seeds, tolerances

Apple Silicon CPU only (`macOS-26.6.2-arm64-arm-64bit`, `arm64`; no GPU, no
MPS, no CUDA). Virtualenv Python 3.12.13; numpy 2.5.2; scipy 1.18.1;
matplotlib 3.11.1.

Commands (all from the experiment root):

```bash
.venv/bin/python -m unittest discover -s src -p "test_*.py" -v
.venv/bin/python src/validate_self_cell.py
.venv/bin/python src/family1_grid_refinement.py
.venv/bin/python src/family1_pilot.py
```

Measured runtimes (recorded by the scripts where available): tests 5.43 s;
quadrature validation 5.12 s (JSON metadata); grid refinement 5.24 s total
(per-N 0.06/0.29/1.18/3.71 s); pilot 2.13 s.

Seeds: pilot map seeds {0,1,2}, pose seeds {10,11,12} (unchanged); source
gradient sample seed 12345; pose-FD unit-test seed 54321; identity sample
seed 2026. No random draws in the validators.

Tolerances used: unit-test quadrature agreement <1e-11; sign FD rel error
<1e-5 at delta 1e-6; pose FD rel error <1e-4 at delta 1e-4; pilot reports raw
FD errors with denominator guard `norm_guard_eps = eps_machine = 2.22e-16`
(no pass/fail threshold is applied inside the pilot); validator JSON records
raw numbers and a `tends_to_zero` flag with explicitly stated criteria.

## 9. SHA-256 digests (live corrected artifacts)

```text
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
e712bebe7f322773ef390ab042075afcedd8efe2369fb38cd2165324df0810c7  src/test_family1_corrections.py
42dd75c079d8b88a1a6d9ed2aa3f00535788fefc461837262d555d0743b611bc  src/validate_self_cell.py
b00cfe6fa6781ffa5c78fc6c04fc9ad59779a33f06876252eac8fe79064a645b  src/family1_grid_refinement.py
4e58d209a58aa8ba5062bfe8433030d1cfa6921efb59eaca5066d226e1f0d9e7  results/family1_pilot_results.json
62e1c128a6d7fa290fef0e526b5cbdf2b816ea7fe4a81dbaff8bc52585c3333e  results/self_cell_quadrature_validation.json
8677545c35437b139445dbc700d0bfbf25d46e5fc7a94295b93170b0da29607e  results/family1_grid_refinement.json
d538df048f07004afbd98e290180111362ec136ee1f8924dba5e37fa18b5e03d  figures/family1_fd_convergence.png
b577a2f38386df5cf7514829a96511a89e0cfbda1b3da85fbd0ac55670cef469  figures/self_cell_quadrature_validation.png
96e56db5694046c11f44fd6e5c0669206275ec9fa2715b449e102177ad04e044  figures/family1_grid_refinement.png
```

The source digests are also embedded in each result JSON's `source_sha256`.

## 10. Claims and limitations

- Corrected `self_cell_green` now matches the audited equal-area disk
  integral and its documented `1/k_b^2` lower endpoint; direct quadrature
  agreement is at roundoff (`~1e-16`).
- `grad_s = -grad_z` is preserved and unit-tested against finite differences
  in the source coordinate; the old sign fails that test at rel error 2.
- All Family-1 pilot checks rerun cleanly on the corrected formula; no
  mathematical ambiguity remained, so no blocking question was raised.
- Passes remain bounded checks of the implemented discrete harness: internal
  derivative consistency and self-cell quadrature fidelity. They are not a
  continuum-convergence proof, physical acceptance test, or novelty claim.
