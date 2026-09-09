# Matched-frequency development experiment implementation

Status: implemented and bounded checks passed on 2026-09-08. **No benchmark
fits have been launched.** Parent source review and explicit run authorization
remain the requested next step. This is an isolated implementation helper;
scientific interpretation remains with parent Codex.

## Files and verification

- `matched_frequency.py`: explicit-frequency model, deterministic acquisition
  assembly, reference checks, registered 12-fit runner and frozen-estimate audit.
- `checks.json`: N8 mixed-frequency derivative and assembly checks.
- `reference_data_checks.json`: checks using actual independent N64 ADDA means,
  plus complete per-frequency ADDA metadata and verified true residuals.

All 13 columns of the mixed `(3,6,9,18,6)` model were checked by centered finite
differences: maximum relative error **5.2783e-9**, below `1e-7`. The material
column error was `1.2446e-10`. Electronics-reference columns also passed.
The full real/imag residual, with reference absent and present, passed directional
checks at `5.70e-10` and `5.69e-10`. Adapter predictions and Jacobians exactly
equal the existing `nonspherical_calibration.Model(8)` on its original band.
The original global `KS` remains unchanged. Repeated-k physical outputs and
Jacobians are bit-identical; only four physical models serve five check blocks.

Actual N64 data sharing passed for both scenes. The absolute proper-complex
noise standard deviations are `0.0002617167877636827` (8101) and
`0.00025944607320512125` (8102). Hashes of low observations, shared electronics
reference, and standardized extra noise are in the JSON artifacts.

## Protocol implementation

The inverse remains the same 13-parameter known-support homogeneous ellipsoid:
one real permittivity with fixed imaginary part 0.03, three receiver translations,
one delay, four log amplitudes and four phases. `ShapeVIE`, dipole radiation,
reference-cache validation and field metrics are imported from existing sources.
No existing physics source or frozen estimate/result is edited.

The physical scenes reproduce the existing seeded scene construction before
making any new noise draws. Every scene makes one low-band noise draw and one
separate extra draw. Low observations are copied exactly into both extended
acquisitions. The extra standardized proper-complex noise is reused for high
k18 and repeated k6, with one common absolute sigma derived only from low-band
mean RMS at 30 dB. Repeated k6 is a new measurement and separate residual block.
Physics and material tangents are solved once per unique k and material state.

The same `(3,6,9)` electronics observations, from the established `seed+10000`
reference RNG and sigma 0.02, are supplied to every reference-present fit. All
choices hold the identical reference array; there is no k18 reference. Proper
complex residuals are whitened with `sqrt(2)/sigma` for real and imaginary parts.

Exactly 2 scenes × 3 choices × 2 reference conditions are prescribed, N32,
identical nominal start epsilon 2 and other parameters zero, original bounds
and scales, `max_nfev=35`, `ftol=xtol=1e-9`, `gtol=1e-7`. Truth is used only
to generate observations and evaluate endpoints, never as optimizer input.

Each estimate is persisted before evaluation. The evaluator constructs a fresh
N32 model using only k3/k6/k9 and the same held receivers `(count=17,radius=1.6)`.
Sensor and electronic-factor-removed structural fields are compared separately
against noiseless independently computed ADDA radiation means, with aggregate
and per-low-frequency field and phase metrics. High training-frequency error
does not enter low-band prediction statistics. Delay, material, pose, amplitude,
gain-phase and low-band electronic-response errors are also recorded.

Setup, optimizer solve and independent evaluation wall times are separated.
Per-frequency ledgers record physical RHS/FFT work, true residuals, projection
time, setup time and acquisition multiplicity. Cache-read cost is separated
from historical ADDA generation cost. Each fit is recorded before it begins;
failures retain error type, message, work where available, and frozen estimates
if optimization finished. Checkpoints are atomically replaced after flushing,
and an exclusive OS lock prevents simultaneous benchmark runners. On restart,
completed attempts (including fit or evaluation failures) are skipped. A
checkpoint containing a frozen estimate with unfinished evaluation resumes
only the independent evaluation, without an optimizer call. An incomplete
`started` checkpoint without an estimate restarts that interrupted fit from
the same nominal start. In either recovery case, the entire prior attempt
metadata is retained in `attempt_history`, with the recovery action and
timestamp. Uncheckpointed interrupted work has explicitly unknown cost; it
is not silently counted as zero. Unknown or duplicate checkpoint states are
refused. Runs from a different source hash are still refused. Completed failed
attempts are never selectively retried.

The lightweight checks additionally exercise pre-estimate interruption,
frozen-estimate evaluation recovery, preserved optimizer metadata and skipped
completed failures. Recovery checks use a synthetic evaluator and invoke no
optimizer or N32 physics. Fit and evaluation interruption histories remain
separate from the prescribed 12 distinct configuration records.

`comparisons.json` reports signed low+high minus low and low+high minus
low+repeat differences for both sensor/structural field and phase metrics,
retaining harmful and helping outcomes. An incomplete comparison is marked
explicitly. No significance or general high-frequency benefit claim is made.

## Reference preparation

The expected k3 N64 cache was initially absent. Parent explicitly authorized
this input-preparation step separately from the benchmark. It was generated
using the existing `reference('ellipsoid',64,3.,2.5+.03j)` function and existing
ADDA executable: 44,264 voxels, 0.005 spacing, two successful propagation runs,
about 2.08 seconds combined executable time, maximum verified recalculated
residual `9.7112e-11`. New cache files are confined to the previously missing
`a3_research/results/adda/ellipsoid_n64_k3_2.5_0.03_ldr/` directory.
Existing k6/k9/k18 caches were read and validated. The experiment runner
preflights all required caches and refuses missing files rather than launching
an unrequested ADDA generation job.

## Commands

From the repository root, executed successfully:

```sh
research/trispace_self_calibration/a3_research/.venv3d/bin/python research/delegated/a3_matched_frequency/matched_frequency.py
research/trispace_self_calibration/a3_research/.venv3d/bin/python research/delegated/a3_matched_frequency/matched_frequency.py --check-reference-data
```

The script sets OMP/OpenBLAS/MKL/Accelerate/NumExpr thread limits to one before
importing NumPy/SciPy. No installation, external API or new service is used.

After parent authorization only, the prepared full-run command is:

```sh
research/trispace_self_calibration/a3_research/.venv3d/bin/python research/delegated/a3_matched_frequency/matched_frequency.py --run
```

It re-runs bounded checks before fitting, then emits `fits.json`,
`comparisons.json`, reference-cost metadata and per-scene data NPZ/JSON files
within this delegated directory. Current checks establish implementation and
data-sharing correctness at N8; they do not establish N32 optimizer success,
the scientific outcome, or production/general-image acceptance.
