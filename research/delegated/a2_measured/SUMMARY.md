# SUMMARY: bounded measured-data pilot (primary 2001 Fresnel, dielTM_dec8f)

> SUPERSEDED FIT INTERPRETATION: the parent rejected the initial wave/passivity
> implementation and test-gain leakage. Read the parent audit notice at the end
> and `research/trispace_self_calibration/a2_research/FINAL_REVIEW.md` before
> citing any fit/model-validation claim here. Raw ingestion remains valid.

Isolated worker pass, no child agents, no AI Scientist. Scope: ingestion QA +
one bounded full-wave 2D TM fit pilot. Main Codex owns all scientific claims.
No solver/physics/highdim code and no original downloaded data were modified;
all outputs live under `research/delegated/a2_measured/`.

## Gates, kept separate

| Gate | Status | Evidence |
|---|---|---|
| download (canonical 2001 files) | MET (A2 literature worker) | `a2_literature/public_data_audit.md`; SHA256 re-verified here |
| ingestion / data QA | MET | `results_ingestion.json`, `test_ingestion.py` (4/4 pass) |
| analytic model self-verification | MET | `results_model.json.self_checks` (conservation, boundary, PEC limit, truncation) |
| model fit vs measured data | BOUNDED, PARTIAL | residual 0.47 (3-param) / 0.25 (4-param); strong only below ~4 GHz |
| actual unknown-geometry recovery | NOT MET | no ground truth exists; only consistency with the published "about 30 mm" |

## Findings (concise)

1. Column semantics are resolved from the fulltext intro: (view, receiver,
   frequency in GHz, Re/Im total E, Re/Im incident E); scattered =
   total - incident; `exp(+i w t)`; TM = E parallel to the cylinder axis.
2. Published geometry: cylinder radius 15 mm, `eps_r = 3 +/- 0.3`,
   centre "about 30 mm" off-axis (direction unpublished), `de = 720 +/- 3 mm`,
   `dr = 760 +/- 3 mm`, 36 views x 49 receivers x 1-8 GHz.
3. Ingestion is clean: 14112 rows, all finite, zero duplicates, full
   frequency coverage, hash matches the prior worker's record.
4. The measured incident field is rotationally consistent across views to
   <1%, and its phase confirms the `exp(+i w t)` convention.
5. A point-source/line-source incident model is rejected: ~99% relative
   residual (horn beam pattern dominates the receiver arc); a full-ring
   equivalent-source fit is ill-conditioned in the near field. These
   unsuccessful checks are preserved in the code path (documented here).
6. The adopted forward model is an analytic Mie series for the dielectric
   cylinder under plane-wave illumination from the emitter direction, with a
   single complex gain per frequency profiled identically for every variant.
   The core is self-verified: Jacobi-Anger err 2e-15, boundary continuity
   7e-16, optical theorem 1 - 1e-15, PEC limit 1.2e-4 (at eps=1e8), series
   tail 2e-54.
7. Identifiability of the joint (cx, cy, eps_r) fit is demonstrated by the
   injected-perturbation pilot: injected centre (24, -8) mm + eps_r 2.7
   recovered to 8 um centre error and 3e-4 eps error under 5% noise
   (residual 0.5%). Injected errors are synthetic, not natural truth.
8. On the real data the recovered centre is stable across band choices at
   (-0.001, -0.027) m, i.e. radius ~26.5-27 mm from the axis - consistent
   with the published "about 30 mm". The direction is a gauge choice and is
   not claimed physical.
9. eps_r is NOT identifiable from this data under this model: it drifts
   5.2 -> 2.5 -> 1.2 as the fitted band grows (1-2, 1-3, 1-4 GHz) and the
   1-8 GHz compromise is 2.10 (3-param) or 1.81 + tan(delta)=0.33 (4-param),
   outside the published 3 +/- 0.3. The 4-param loss is an effective
   parameter absorbing model error, not a physical loss measurement.
10. Residuals are small at 1-2 GHz (2-8%) but large at 4-8 GHz (0.4-0.9),
    consistent with the plane-wave incident model breaking down where the
    horn beam narrows; this is model misspecification, not a data defect.

## Artifacts

- `METADATA.md` - sources, column semantics, geometry, unknowns, evidence levels
- `ingestion.py` / `test_ingestion.py` - parser + QA + executable tests
- `model.py` - Mie forward model, self-checks, nominal vs joint fits,
  low-frequency diagnostic, robustness pilot
- `results_ingestion.json`, `results_model.json`
- `plots/qa_amplitude_phase.png`, `plots/qa_incident_rotation_consistency.png`,
  `plots/fit_scattered_comparison.png`, `plots/fit_residual_and_stability.png`,
  `plots/robustness_basin.png`
- `sources/Ip01Introduction_Belkebir.pdf` (+ `.txt`)

## Limitations

- Plane-wave incident model with per-frequency scalar gain cannot represent
  the horn's angular pattern, near-field, receiver antenna factors, or any
  per-view phase-centre jitter; this dominates the 4-8 GHz residuals.
- Field units are unpublished; gains absorb absolute scale per frequency.
- The offset direction is a measurement-gauge ambiguity (target-rotation vs
  source-rotation equivalence up to reflection); no physical direction is
  claimed.
- Loss and eps_r values are effective/compensatory under the misspecified
  forward model; no uncertainty quantification beyond multi-start agreement.
- This is a bounded single-target pilot on one file; it is not a pixel-domain
  inverse-scattering reconstruction and proves nothing about the project's
  algorithm on measured data.

No algorithm-superiority claim is made from this exploratory dataset.
# Parent audit notice, 2026-09-06

This original worker report is retained for traceability, not adopted as the
final interpretation. Its initial Mie fit mixes propagation/passivity signs
and re-profiles gains on test views. Model-fit conclusions and the claimed
unknown-geometry validation are superseded by
`research/trispace_self_calibration/a2_research/FINAL_REVIEW.md` and the corrected
pilot recorded there. Data ingestion remains valid. The primary descriptor
does explicitly state exp(+i omega t); the issue is implementation consistency.
