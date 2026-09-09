# A3 reproducibility entry points and limits

Workspace root: `/Volumes/migodam's-external-brain/Research/Inv_SLAM`.
Run the commands below from that directory. They reuse the installed A3 CPU
environment; they do not install providers, CUDA, services or dependencies.
The complete project is still a working research package, not final TAP release.

## Environment and preserved inputs

- Numerical interpreter: `research/trispace_self_calibration/a3_research/.venv3d/bin/python`.
  Current scientific stack recorded during this stage: NumPy 2.5.3, SciPy
  1.16.3, Treams 0.4.7, Matplotlib 3.11.1, Apple Silicon CPU.
- Use `OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1` for
  serial timing. Do not run other heavy jobs during a timed comparison.
- Original Theory/Questions/A3.md and A3_2.md and A2 results remain preserved.
  Some A3 adapters import `research/delegated/a2_physics` and
  `research/delegated/a2_measured` read-only; this is not a self-contained wheel.
- ADDA revision, compiler options, static FFTW archive hash, unit/polarization
  conversions and independence limits are in `external/BUILD_PROVENANCE.md`.
  Use that record rather than substituting current upstream HEAD silently.
- Fresnel input: `research/delegated/a2_literature/data/2001_iop_17_6_301/dielTM_dec8f.exp`.
  SHA256 `476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb`.
  Preserve source/measurement conventions and training-only gains.
- Provider credentials are not experiment inputs or release artifacts. Do not
  bundle User/API.md, user configuration, virtual environments or credentials.

## Checks and short analyses

Prefix each script below with the numerical interpreter above. Tests write
their own check JSON; summaries never replace raw fits.

| Entry point | Purpose |
|---|---|
| `a3_research/nonspherical3d.py --controls` | Independent-code radiation/polarization controls; needs the local ADDA executable |
| `a3_research/nonspherical_calibration.py --test` | Known-ellipsoid physical Jacobian |
| `a3_research/measured_extension.py --test` | Cylindrical illumination/Graf identity |
| `a3_research/multifidelity_calibration.py --test` | Value/tangent anchor consistency and derivative |
| `a3_research/audit_multifidelity.py` | Recompute all ten paired endpoints and cost conclusions |
| `a3_research/rom_seed_comparison.py --test` | Four chart seeds, orthogonality and fixed-chart derivative |
| `a3_research/audit_rom_seeds.py` | Recompute 36 completed comparisons and descriptive figure |
| `a3_research/nonlinear_gain_graph.py --test` | Nonlinear complex-gain derivatives and tree construction |
| `a3_research/timed_rom_development.py` | Deadline implementation tests on inspected development data |
| `a3_research/analyze_rom_final.py` | Refuses incomplete final sample; analyzes only all 240 final rows |
| `a3_research/audit_matched_frequency.py` | Audits all 12 matched low/high/repeat fits and actual shared data |
| `a3_research/check_translation_clock_gauge.py` | Nine full-wave field-translation checks and illumination-direction ranks |
| `a3_research/check_frequency_information.py` | Standard information increment and biased-risk algebra controls |
| `a3_research/audit_discrepancy_weighting.py` | Requires all 16 new weighting outcomes; verifies frozen inputs and paired contrasts |
| `a3_research/audit_correct_model_frequency.py` | Verifies all 12 intentionally shared-model controls; separates post-hoc oracle mode alignment |
| `a3_research/shape_material_stress.py` | N8 two-shape derivative/equivalence and whitening checks, no full fits without --run |
| `a3_research/audit_shape_material_stress.py` | Requires all 32 registered stress outcomes, reports every case and reference sensitivity |

Here `a3_research/` in the table abbreviates
`research/trispace_self_calibration/a3_research/`, not a second directory.

## Longer experiments and checkpoint behavior

- `nonspherical_calibration.py --seeds 8101 8102 --ns 16 32` runs the declared
  sixteen inverse fits. `nonspherical3d.py --shape ... --ns ... --ks 9 18`
  supplies the separate forward-resolution studies; use the exact grids listed
  in the raw records (box boundary refinement is nonmonotone).
- `multifidelity_calibration.py` runs five controls on two already inspected
  scenes under MULTIFIDELITY_PROTOCOL.md. `measured_extension.py` runs the nine
  fold/method combinations with four training-selected starts.
- `rom_seed_comparison.py --q 9`, then `--q 49`, run the matched-seed development
  controls serially. Timings include setup; offline exact derivative audits do
  not enter chart admission or reported online time.
- `nonlinear_gain_graph.py` is a one-scene finite-bank mechanism control with
  60 gain fits, NOT independent-model inversion or a generated-bank success.
- `rom_final.py` runs the fixed forty-scene final protocol, 240 method fits.
  It saves each completed method and refuses resume if protocol/source digests
  change. During collection the log reports counts only. Do not inspect final
  outcomes to tune algorithms. Analyze only after all records exist.
- `research/delegated/a3_matched_frequency/matched_frequency.py --run` runs
  the matched acquisition control under MATCHED_FREQUENCY_PROTOCOL.md. Its 12
  existing endpoints and NPZ observations are immutable inputs to later work.
- `research/delegated/a3_discrepancy_weighting/discrepancy_weighting.py --run`
  runs 16 frozen-weight continuations under DISCREPANCY_WEIGHTING_PROTOCOL.md.
  Without `--run`, it executes bounded checks only. It imports matched sources
  read-only, hashes all inputs, stores each low-pilot-derived mode, and saves
  each endpoint before staged low/high evaluation. Completed failures are not
  selectively retried. Historical pilot cost plus new costs is not fresh timing.
- CORRECT_MODEL_FREQUENCY_PROTOCOL.md registers a separate shared-model
  attribution control. It must not be described as independent numerical
  validation; check its own result records before claiming execution.
- `research/delegated/a3_correct_model_frequency/correct_model_frequency.py --run`
  completed the 12 matched-noise shared-model controls. New means and all fits
  have immutable source/input manifests; original ADDA data are untouched.
- `shape_material_stress.py --run` executes four new support/material cases,
  32 pilot/continuation fits, and independent-reference N64/N96 sensitivity
  checks under SHAPE_MATERIAL_STRESS_PROTOCOL.md. All shape/material cases
  remain known-support scalar-permittivity inversions. Failed or incomplete
  reference generation is preserved and requires inspection before resume.

Several scripts resume by skipping existing rows; invoking them again is a
checkpoint check, not an independent rerun. A fresh reproduction needs a
separate workspace/output copy with preserved immutable inputs and no prior
result records. Do not delete original result files to manufacture a clean run.
Final source digests are in `results/rom_final_source_manifest.json`. Older
development artifacts do not retroactively acquire that version guarantee.

## Reading results correctly

Start with ACCEPTANCE.md, CONTRIBUTION_EVIDENCE.md and METRIC_SEMANTICS.md.
The current paper is PAPER_DRAFT_A3.md; the explanatory report is
`communication/A3_RESEARCH_REPORT_ZH.md`. The actual autonomous pipeline's
execution, rejected drift and failed writeup are in PIPELINE_AUDIT.md.
An isolated Codex review and recheck are in
`research/delegated/a3_current_critic/`; neither is external specialist review.
# A3_3 supplement checks (2026-09-08)

`check_a3_3_task_risk.py` produces `results/a3_3_task_risk_checks.json`:
12 fixed-linear risk identities, profiled orthogonality invariance, a raw-angle
counterexample and a two-world bound check. These are deterministic algebra,
not independent electromagnetic validation. Run with the existing `.venv3d`
Python; no additional dependencies or paid model calls.

`calibrate_then_image.py --run` follows CALIBRATE_THEN_IMAGE_PROTOCOL.md;
all32 fits are complete. `audit_calibrate_then_image.py` checks source/input
hashes, common nuisance starts, fixed geometry and oracle labels, producing
`results/calibrate_then_image_audit.json`. Original calibration sources and
datasets were not changed. The new Appendix H and source audit are separate
from frozen experiment provenance.
