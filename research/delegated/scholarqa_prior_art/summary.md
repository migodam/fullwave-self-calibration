# Prior-art synthesis summary (bounded to retrieved records)

## Findings

1. **No retrieved paper defines the proposed combination.** Across the 40
   candidate abstracts, none analyzes pose-nuisance-eliminated effective map
   information, principal-angle retention, a pose-map gauge, multifrequency
   pose compensation, or trajectory sensitivity for full-wave volumetric
   inverse-scattering SLAM. No abstract even mentions SLAM, pose, trajectory,
   Fisher/EFIM, or principal angles.
2. **Closest partial overlaps are joint-unknown inversions, not SLAM:**
   singular-value-optimized sensor geometry and frequency sampling
   (Capozzoli 2023); scatterer surface position/geometry parameter
   identification under observability conditions (Xu Zhang 2025);
   joint shape + unknown impedance via multifrequency continuation (Borges
   2021); statistics of random obstacle geometry (Zhiqi Sun 2026). Each keeps
   the sensors fixed and estimates a scene or model unknown.
3. **The only spectral-subspace substrate retrieved is generic SOM
   algebra** (Chen 2010; Zhong 2010/2011; Miao Wang 2024; Siampour 2024;
   Tingsen Zhang 2025): SVD of the induced-current to field operator
   partitions signal vs ambiguous subspaces. This is a static, known-geometry
   substrate — separated from, not supportive of, the pose-aware combination.
4. **Multifrequency in the retrieved record compensates contrast, boundary
   impedance, or noise-regularization — never pose.** (Miao Wang 2024,
   Borges 2021, Capozzoli 2023.)
5. **Feasibility constraints are explicit in retrieved work:** weak/strong
   contrast instability (Wei 2018), minimax-optimal logarithmic Bayesian
   contraction (Furuya 2024), Nyquist/polarization observability requirements
   (Xu Zhang 2025), and aspect-limited acquisition loss (Capozzoli 2023).
6. **Identifiability theory exists in adjacent settings** — Schiffer almost-
   sure uniqueness (Hongyu Liu 2024) and fixed-angle uniqueness under
   variable sound speed (Oksanen 2026) — but none addresses pose-map gauge
   structure or multi-pose trajectories.
7. **The retrieval is materially incomplete:** 8 of 10 search operations and
   all snippet searches failed with 429, so EFIM/Schur algebra, wave/RF SLAM,
   bilinear gauge, unknown-sensor inverse scattering, and robust-design
   literatures were not screened; no Tier A (full-text) evidence exists.

## Closest-overlap assessment (bounded to retrieved records)

Among the retrieved records the closest neighbors are Borges 2021 (joint
shape + unknown model parameter, multifrequency), Xu Zhang 2025 (position/
geometry parameter observability), and Capozzoli 2023 (information-driven
measurement geometry). None approaches pose-trajectory estimation or
nuisance elimination, so the closest-overlap distance is large: the proposed
formulation is unrepresented in this bundle. This is not a global novelty
claim.

## Artifacts

- `screened_papers.md` — 21 selected / 19 excluded, grouped and tiered.
- `claim_ledger.md` — 21 bounded claims (all Tier B; contextualizes/
  qualifies only) mapped to verified IDs/DOIs/arXiv.
- `verified_references.json` — canonical `scholarqa-cli verify` output;
  21/21 resolved, 0 unresolved.
- `SEARCH_PLAN.md` and `evidence.json` preserved unchanged.

## Uncertainties

- Abstract-only screening: several high-value SOM classics required OpenAlex
  abstract recovery because Semantic Scholar elides them; [30] is OpenAlex
  title-matched (S2 record has no year/abstract).
- arXiv-only items ([29], [32], [25]) are unrefereed; two are dated 2026.
- 19 exclusions may hide relevant variants, especially niche no-abstract
  SOM applications.
- The two successful queries bias the candidate set; the absent categories
  could contain the true closest prior work.

## Where Codex intervention is needed

- Decide whether to re-collect the four 429-failed queries (especially
  EFIM/Fisher, wave/RF SLAM, bilinear-gauge, unknown-sensor, robust-design
  formulations) with an API key or backoff before the novelty judgment.
- Interpret whether the SOM singular-spectrum substrate and the
  geometry/information tradeoffs above count as enabling antecedents or as
  distinct prior art.
- Make the final novelty judgment, which remains outside this worker's
  authority.
