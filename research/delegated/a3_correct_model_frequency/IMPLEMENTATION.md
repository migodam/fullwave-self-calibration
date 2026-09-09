# Correct-model frequency control implementation

Status: implemented for parent review on 2026-09-08. Lightweight checks pass.
**Zero physical solver calls, no N32 means generated, and no benchmark fits run.**
The discrepancy-weighting run and existing matched experiment remain untouched.

## Files and accepted bounded checks

- `correct_model_frequency.py`: registered 12-fit control, importing immutable
  matched-frequency helpers and their existing physics.
- `checks.json`: provenance, exact saved-input reuse, acquisition sharing,
  shape/electronic-convention and inherited checkpoint recovery checks.

The checks use deterministic, nonphysical field fixtures with the required
array shapes. They do not generate random noise or call an inverse model.
For both scenes they verify exact saved truth, absolute sigma, low noise,
extra noise, electronics references and receiver arrays; exact low observation
identity across all three choices; independent repeated k6 observation; and
the saved standardized extra-noise reuse for high and repeat blocks. They also
verify electronic-factor removal, prior nominal starts, and unchanged hashes
of imported sources, frozen matched results and input NPZs. The existing
validated physical derivative implementation is not changed or retested here.

## Sole mean substitution

At benchmark start only, N32 `ExplicitFrequencyModel` generates training means
for k3/k6/k9/k18 from the stored true parameters and training receivers. Another
N32 low-band model generates held-out sensor means at the stored held receivers.
The structural means divide out true electronics with the inherited convention.
These intentional shared-model means replace the independent ADDA means.

The assembly directly applies the saved arrays:

- low observations = new low means + saved sigma × saved low_noise;
- high observation = new k18 mean + saved sigma × saved extra_noise;
- independent repeated k6 = new k6 mean + saved sigma × saved extra_noise.

The saved noisy low-band electronic reference is copied without regeneration.
Sigma remains the original absolute scale even if relative SNR changes. There
are no additional gains, weights, noise draws, priors, tuning or warm starts.
The new mean/data archives are saved only in this delegated directory.

## Twelve fits and evaluation

The runner iterates two scenes × low/low+high/low+repeat × reference absent or
present. Every N32 fit starts from epsilon=2 and the other 12 parameters zero.
Bounds, scales, `max_nfev=35`, `ftol=xtol=1e-9`, `gtol=1e-7`, residual assembly,
and analytic/material and translation derivative implementation come from the
matched experiment. Truth is used only in generation and frozen evaluation.

Each optimizer endpoint is saved before any offline metrics. Evaluation uses
the same low-band held-out receiver channels and separate sensor/structural
field and phase errors, plus pose/material, delay and electronics errors.
The comparisons retain high-minus-low and high-minus-repeat contrasts within
this control, and each configuration's difference from the corresponding
independent-code matched configuration. Failures remain visible and are not
silently dropped from comparisons.

The experiment is explicitly a shared-N32-generator/inverse attribution
control. It supplies no independent Maxwell, continuum-accuracy, measured-data,
success-rate or population inference. Parent owns interpretation of whether
the resulting contrasts support a model-discrepancy explanation in these cases.

## Checkpoints, failures and timing

Generated archives use flushed atomic NPZ replacement. JSON records use the
inherited flushed atomic writer. A manifest binds the new source, protocol,
frozen matched inputs/results and imported physics sources. Every saved mean
archive has a content hash checked on reuse. Source/input changes, duplicate
fit states and changed cached data are refused.

Interrupted mean generation may restart with previous attempt metadata retained.
Completed generation failures remain recorded and propagate to the affected
fit configurations without automatic selective retry. Interrupted fits without
an estimate restart from the same nominal start, preserving their history.
Frozen estimates resume only evaluation. Completed fit or evaluation failures
are retained and skipped. Unknown interrupted work costs are not counted as zero.

Failures distinguish mean generation, fit setup, fit solve, parameter evaluation
and low-band evaluation. Available physical ledgers, stage/attempt wall times
and frozen endpoints remain preserved. A failed mean generation can report
elapsed time without unreturned internal work counters; no false zero-work
claim is made.

Data-generation setup/solve work and wall time are recorded once per scene and
marked shared across its six fits. Fit setup and solve, fit attempt wall time,
and frozen evaluation costs remain separate. The runner acquires its own
exclusive lock and the existing matched/discrepancy benchmark locks before
generating any mean. It therefore refuses a simultaneous experiment run; the
parent will still launch it only after discrepancy timing collection completes.
Existing lock files are opened without truncation or content changes.

## Commands

Executed successfully from the repository root:

```sh
research/trispace_self_calibration/a3_research/.venv3d/bin/python research/delegated/a3_correct_model_frequency/correct_model_frequency.py
```

Prepared full command, **not executed**:

```sh
research/trispace_self_calibration/a3_research/.venv3d/bin/python research/delegated/a3_correct_model_frequency/correct_model_frequency.py --run
```

The script uses the existing environment and sets single-thread CPU limits
before NumPy/SciPy imports. No API, service, installation or additional worker
is invoked. Full-run authorization and scientific review remain with the parent.
