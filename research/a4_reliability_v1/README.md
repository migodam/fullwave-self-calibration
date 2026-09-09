# A4: restricted full-wave phase-attribution research package

**Status: working research; not a submission-ready TAP method paper.**

Start with [the scientific verdict](docs/EXECUTIVE_SCIENTIFIC_VERDICT.md), then [PAPER_DRAFT_A4.md](PAPER_DRAFT_A4.md). The retained route is B: an exact geometry-information condition for anchored radial targets under calibrated electric-l=1 modal illumination. It is not universal tomography or a cleared novelty claim.

## What is included

| Need | File |
|---|---|
| Manuscript | `PAPER_DRAFT_A4.md` |
| Assumptions, proofs and counterexamples | `docs/THEORY_A4.md` |
| Algorithm / risk / six-action limits | `docs/ALGORITHM_A4.md` |
| Executed experiments and all comparison tables | `docs/EXPERIMENT_REPORT.md` |
| Prior-art and narrative audit | `docs/LITERATURE_AND_NARRATIVE_REVIEW.md`, `docs/REFERENCES.md` |
| Strict reviewer disposition | `docs/TAP_REVIEW_AND_CLOSURE.md` |
| No-hardware-yet protocol | `docs/HARDWARE_PROTOCOL.md` |
| Corrections / input boundary | `docs/PROVENANCE_AND_CORRECTIONS.md` |
| Frozen source/protocol hashes | `provenance/frozen_manifest.json` |
| Raw endpoints / noisy measurements | `results/frozen_*.json`, `results/frozen_*_data.npz` |
| Fresh diagnostic correction | `supplements/reference_audit_v2/` |
| Analysis / figures | `analyze_a4.py`, `results/A4_ANALYSIS.json`, `figures/` |

## Reproduce checks without touching frozen data

Tested with Python3.13.5; exact numerical-library versions are in `requirements.txt` and the environment record. Install into a virtual environment before running.

```sh
python -m pip install -r requirements.txt
sh run_tests.sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analyze_a4.py
python audit_package.py
```

The test launcher disables external pytest plugins because `src/coverage.py` shares its name with the third-party coverage tool. The default-import failure is preserved and explained. All 21 intended cases pass.

`analyze_a4.py` reads frozen endpoints and writes derived tables/plots; it does not rerun or overwrite frozen experiments. `audit_package.py` only checks integrity and stored-data invariants. Use `rerun_fresh.py` with a **new** output directory for a new experiment; it refuses to overwrite a destination. The historical frozen CLI refuses to overwrite an already completed frozen collection, but the safe fresh-run wrapper is preferred for new cases.

## Results that change the scientific decision

The risk selector matches full-fine estimates in 8/8 cases but has a median paired charged time ratio of 1.861. It is removed from main claims. Coverage retains unfinished cells and rejects several failure types but accepts all 4 explicitly indistinguishable wrong-world discrepancies when beta=0 is unjustified. Only the restricted EM positive condition remains a main-contribution candidate, with unresolved prior overlap.

## Repository and evidence boundary

Place this directory at `research/a4_reliability_v1/` in the supplied repository. This is a new standalone implementation, not the omitted A3 original source. No A3 frozen endpoints were changed; no remote commit/PR was made. The 21 checks,41 spectral points,48 DDA states,24 branch scenes,8 joint scenes and fresh16 reference scenes have different statistical units and must not be pooled into a single success rate.

All figures are generated from preserved data and use no fabricated measured results. This package includes no third-party article PDFs, font files, credentials or original user-file copies.

Example fresh run, with an explicit new output directory outside this package:

```sh
python rerun_fresh.py --output ../a4_fresh_20260909 --seed 20260909 --branch-per-kind 1 --calibration-scenes 2
```

The wrapper itself was syntax/CLI checked. The stored experiments were run using the frozen source entrypoints, not regenerated using this post-freeze packaging wrapper.
