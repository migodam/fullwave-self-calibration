# A2 PRASC-SOM theory validation: independent replication summary

Bounded numerical sanity / negative-control supplement. All artifacts were
produced in this folder on 2026-09-06 by an isolated worker. No physics,
solver, global configuration, or original experiment package was modified;
the copied code was byte-identical (SHA-256 verified) to the source package
at `experiments/idea_loops/loop_2026-09-05_17-55-24_a2_development/
experiment_a2_prasc_theory_validation`, and the shared core
`research/delegated/a2_physics/physics.py` was imported read-only.

This file does **not** assert that a previous "round-2 replication" existed;
every number below is from a fresh execution in this directory.

## 1. Environment and exact executed commands

Environment: macOS Python 3.13.13 in the shared experiment venv; NumPy 2.5.2;
SciPy 1.18.1; Matplotlib 3.11.1; pytest not installed. All commands below ran
from this replication root with `OPENBLAS_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1`.

```sh
cd "/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_theory_validation/replication"
PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python"

# 1. 50-test suite (assertion runners, exit code 0 each)
"$PY" tests/test_e1.py
"$PY" tests/test_e2.py
"$PY" tests/test_e3.py
"$PY" tests/test_e5.py
"$PY" tests/test_physical_runner.py

# 2. Full experiment reruns to this folder's results/figures/
"$PY" scripts/run_e1.py
"$PY" scripts/run_e1_physical.py
"$PY" scripts/run_e2.py
"$PY" scripts/run_e3.py
"$PY" scripts/run_e5.py
"$PY" scripts/run_e3_physical_limited.py
"$PY" scripts/run_physical_optional.py
"$PY" scripts/compute_ci_stats.py

# 3. New threshold / negative-control ablations (this replication)
"$PY" scripts/run_ablations.py

# 4. E5 physical acquisition extension: 8 receiver configurations
"$PY" scripts/run_e5_physical_acquisition.py
```

Logs: `results/*_run.txt`, `results/run_*_log.txt`, and script-generated
`results/*_commands.txt`. Test suite: `results/replication_test_suite_report.txt`.

## 2. Run counts, sample sizes, and runtimes

| Component | Counts / sample sizes | Shell wall |
|---|---|---:|
| Test suite | 50/50 pass: E1 9, E2 11, E3 15, E5 14, physical smoke 1 | ~1-2 s total |
| E1 synthetic | 5 models (4 matched + mismatch), seeds 101-110 x 2000 | ~3 s |
| E1 physical Fisher | 10 seeds (101-110), N=8 full aperture, 288 complex / 576 realified rows per seed; every intensity row evaluated by 1-D ncx2 quadrature | ~15 s (13.37 s recorded) |
| E2 | exact controls + 36 seeded tangent checks (201-212) | <1 s |
| E3 synthetic | dual spectra; Theorem-7 MC seeds 301-310 x 1000 | ~1 s |
| E3 physical limited | 10 seeds x {full, limited} = 20 physical records | ~2 s |
| E5 synthetic | minimal pair, 24 budget records (401-412), policies | ~1 s |
| Physical optional | seeds 201-212 (phase1) and 401-412 (phase2) | ~1 s |
| Ablations | A1 144 rows + 36 T5 rows, A2 63 rows, A3 81 rows, A4 27 rows | ~2 s (0.75 s recorded) |
| E5 physical acquisition | 8 receiver configs, K=3, decision seed 501, truth seeds 1001/1002 | ~1 s (0.44 s recorded) |

## 3. Fresh rerun metrics (all in results/*.json and *.csv)

### E1 physical Fisher scene (optional script, rerun here)

- Loewner `J_coh - J_ph` within backward tolerance: pass (all 10 seeds and CI).
- Strict positive information `0 < tr(J_ph) < tr(J_coh)` on all 10 seeds: pass.
- Mean per-seed `tr(J_coh) = 1.3589`, `tr(J_ph) = 2.1925e-3`, mean
  `min_eig(J_coh-J_ph) = 1.756e-3` with 95% CI lower bound `1.701e-3`.
- Per-seed quadrature: all rows used quadrature (no lambda fast path);
  max quadrature error estimate recorded per seed.

### E3 physical full vs limited aperture (optional script, rerun here)

- Theorem-3 residuals: all 20 records <= 1e-8 (pass).
- Pose-prior Loewner sandwiches (Lambda=0, 0.1, 1): all 20 records pass.
- Paired full-minus-limited map-retention trace: mean `+0.02566`, 95% CI
  `[0.02354, 0.02778]` (limited aperture lowers retained map information).
- Paired pose-information difference: mean ~1e-27 with CI containing 0.
- Non-vacuity, explicit: after the declared C1 envelope (4 left singular
  directions of the realified tangent), visible pose `B_v` is numerically zero;
  Theorem-7 MC was evaluable in **0/20** records and was skipped (the
  script's `all_mc_within_3se=True` is therefore vacuous over 0 MC rows and is
  not counted as Theorem-7 evidence here).

### Physical optional tangent runner (E2/E3/E5), rerun here

- Phase 0: forward OK; realified `A_r` (576,9), `B_r` (576,3).
- Phase 1: T5b 24/24, max T5a backward-scaled residual 0.0101 (< 100),
  max Theorem-3 residual 3.17e-11 (<= 1e-6), prior sandwiches pass 12/12;
  Theorem 7 skipped 12/12 (rank-deficient visible `B_v`), as in E3 limited.
- Phase 2: innovation relative residual max 4.48e-16, budget identity max
  6.64e-17, `J_stack >= sum J_l` 12/12, all 12 rank-budget steps used the
  singular-fallback path.

## 4. New threshold / negative-control ablations (scripts/run_ablations.py)

Raw data: `results/ablations_results.json`, `results/ablations_summary.csv`;
figure `figures/ablation_2x2.png`. All threshold combinations are recorded;
none were pruned.

| Axis | Design | Result |
|---|---|---|
| A1 physical vs free-current envelope | 12 physical seeds; C0 = exact reduced state (no free-current block) vs declared envelope widths 2/4/8; rank/support thresholds 1e-12, 1e-10, 1e-8 | Physical `min eig(Jx)` positive on 12/12 seeds (range 6.73e-3 .. 1.09e-2). Envelope `min eig(Jx)` ~0 (|value| <= 8.5e-16) for C2/C4/C8 on all seeds (36 saturated records). T5 checks 36/36, max backward-scaled residual 0.0101. Zero envelope visibility is recorded as envelope failure only (G0-D), not physical nonidentifiability. |
| A2 relative-only vs absolute + bias gates | b_zero and small-residual/large-bias exact fixtures; 3 retention tolerances x 3 absolute gates x 3 residual gates | Relative-only gate false assurance in 9/9 b_zero combinations (rho=1 while `Jx=0`). Residual-only gate false assurance in 6/6 small-residual/large-bias combinations (residual 1e-3 vs pose bias 1.0). Confounding eps sweep crossover at eps=1 recorded in raw rows. |
| A3 shared-map vs independent-map | complementary pair + random 3-frame + hidden common compensator; 3 noise scales x 3 PSD gates x 3 visibility gates | Exact stack >= sum on all complementary records at every PSD gate. Independent-map (per-frame-sum) criterion reports no visibility in 27/27 complementary records while the shared stack is visible: recorded false assurance for the independent-map mismatch control. Exact examples preserved (minimal pair, compensation criterion, stack-vs-sum all pass). |
| A4 greedy vs pair-lookahead | fixed 8-candidate fixture; 3 regularisation gates x 3 pose-noise scales x 3 certification tolerances | Greedy gap range [0, 0.08552]. Exhaustive strictly better in 7/9 combos; greedy is optimal in 2/9; pair-lookahead recovers the exhaustive optimum in 3/9. Raw subsets/logdets for all 9 combos and 12 random seeds are in the JSON. |

## 5. E5 physical acquisition: 8 receiver configurations

Script: `scripts/run_e5_physical_acquisition.py`. Model support check:
`physics.Config(N=8, aperture='full', receiver_angles_full=<12 angles>)` is a
supported configuration, so the extension was run.

- Old model: full-ring configuration, frequencies 0 and 1 (shared map/pose).
- Candidate actions: the same new frequency (freq 2) measured with each of 8
  fixed 12-channel receiver configurations (full ring, rotated ring, 4 arcs,
  front half-arc, front/back twin arcs).
- Decision stage used only the fixed declared reference state (seed 501).
  Truth seeds 1001/1002 were used only to evaluate predeclared subsets;
  no truth tangent entered a decision and no nonlinear reconstruction was run.
- Reference exact shared-map logdet: greedy = pair = exhaustive at -4.40688
  (subset `{arc12_back, arc12_right, ring12_full}`), greedy gap 0.
- Incorrect independent-map criterion: exhaustive-independent selects
  `{ring12_full, ring12_rot15, arc12_back}` at -4.47085, i.e. an exact-criterion
  gap of 0.06397.
- Random best-of-12 (seed 409) reaches only -4.55130 at the reference.
- Truth transfer of the exact scores: Spearman rho 0.9943 (seed 1001) and
  0.9951 (seed 1002); max logdet drift 0.064 and 0.330.
- Charged work: 2340 RHS-equivalent solves, 1440 operator products for 3
  physical states (1 decision + 2 truth); all 168 triple logdet evaluations
  were offline matrix operations after precomputing tangent blocks.
- All forwards satisfy the physical state residual test (max ~5.1e-16).

## 6. Evidence tiers

- **T1 deterministic exact algebra** (no randomness): E2 exact controls, E3
  dual-spectra/confounding closed forms, E5 minimal pair/budget identities,
  A2 gate fixtures, A3 exact stack controls.
- **T2 seeded Monte Carlo**: E1 synthetic score checks (seeds 101-110 x 2000),
  E3 Theorem-7 MC (301-310 x 1000), E5 seeded budget/policy records.
- **T3 physical tangent smoke on the N=8 scalar core** at declared
  reference/truth states: E1 physical Fisher, E2/E3 physical optional, E3
  full-vs-limited, E5 physical acquisition; no nonlinear inversion.
- **T4 matched-budget nonlinear E4 comparison**: not run (outside this
  bounded supplement).

## 7. Protocol ledger for this replication

| Protocol element | Status | Evidence |
|---|---|---|
| Common conventions (matched phaseless law, shared map, no oracle rank, fixed anchors) | pass | synthetic modules + physical runner settings; raw settings in JSON |
| E1 scalar matched models + mismatch control | pass | fresh `e1_results.json` |
| E1 physical full-wave Fisher scene | pass (fresh) | `e1_physical_results.json`; strict positive info and Loewner PSD, N=10 |
| E2 nested-rank exact controls + 12 physical linearizations | pass | `e2_results.json`, `physical_optional_results.json` |
| E2 numerical-rank/saturation controls | pass | recorded absolute singular values; roundoff never counted as physical rank |
| E3 synthetic dual spectra, priors, Theorem-7 MC, confounding | pass | `e3_results.json`; MC within 3 SE |
| E3 physical full and limited aperture | executed with explicit non-vacuity | 20 records; Theorem-7 not evaluable (0/20, saturated `B_v`), priors/Theorem-3 pass |
| E3 nearly-parallel map/pose-range fixture | not separately implemented | partial only via b_zero/scaling controls |
| E4 primary matched-budget method comparison | not run | outside this package |
| E5 shared-map algebra, gauge, budget, policies (synthetic) | pass | `e5_results.json` |
| E5 physical 8-candidate acquisition comparison | executed tangent-level | `e5_physical_acquisition_results.json`; reference-decision + separate truth evaluation |
| E5 source-stabilizer (isotropic vs directional) | not run | preserved as `not_yet_run` |
| Publication evidence boundary (3D Maxwell / measured array) | not run | outside this package |

## 8. Limitations and non-claims

This replication is a bounded numerical sanity/negative-control supplement:
N=8 scalar full-wave tangents, small declared rank/current spaces, one
reference and two truth states for the physical acquisition diagnostic, no
nonlinear inversion, no E4 endpoint, no method-superiority, no universal
frame-count, basin, hardware, or publication-acceptance claim. Core theory and
final interpretation remain with the Codex orchestrator.
