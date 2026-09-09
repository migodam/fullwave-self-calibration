# Screened papers: full-wave inverse-scattering SLAM pose-confounding thesis

Screening date: 2026-09-03
Bundle: `research/delegated/scholarqa_prior_art/evidence.json` (40 candidates)
Queries that actually produced candidates: 2 of 6

- `subspace optimization inverse scattering singular spectrum pose` (20 papers)
- `multifrequency inverse scattering unknown measurement geometry` (20 papers)

All other query runs and every snippet search failed with Semantic Scholar
429. There are therefore **no retrieved snippets, no full-text passages, and
no SPECTER embeddings** in the bundle. Every substantive tier below is B
(abstract). Tier A is unavailable by construction. Metadata alone is Tier C.

## Selection result

21 of 40 candidates selected; 19 excluded as irrelevant or near-duplicate.
No candidate combines full-wave volumetric inverse scattering with SLAM,
trajectory/pose estimation, pose-nuisance elimination, EFIM, principal
angles, or a pose-map gauge structure. The categories below record the
nearest retrievable neighbors, not the proposed combination itself.

## 1. Closest overlap (partial only — none is the combination)

These are the retrieved papers whose unknowns come closest to a joint
"geometry + nuisance" estimation problem, but none estimates a moving
sensor trajectory or eliminates a pose nuisance from map information.

| Paper (evidence tier) | Why it is the nearest neighbor | Gap vs the thesis |
|---|---|---|
| Capozzoli et al. 2023, SVD optimization for multifrequency multimonostatic inverse scattering, URSI Radio Sci. Lett., DOI 10.46620/22-0065 (Tier B) | Optimizes transmitter/receiver locations and frequency sampling via singular values under the Born approximation; studies aspect-limited acquisition | Design-time sensor placement, not online pose estimation; no trajectory, no nuisance elimination |
| Xu Zhang et al. 2025, Parameter identification for electrically large curved surfaces, IEEE TAP, DOI 10.1109/TAP.2024.3524008 (Tier B) | Inverts surface geometry/position parameters from scattered waves under explicit observability conditions (Nyquist frequency sampling, polarization, mono/bistatic modes) | Inverts the *scatterer's* static shape/position, not the sensor's pose; no SLAM loop |
| Borges et al. 2021, Multifrequency inverse obstacle scattering with unknown impedance, Adv. Comput. Math., DOI 10.1007/s10444-021-09915-1 (Tier B) | Jointly recovers shape and an unknown impedance function via frequency continuation; shape recovery is robust even to wrong boundary conditions | The joint unknown is a boundary impedance, not a trajectory; fixed receivers |
| Zhiqi Sun et al. 2026, Inverse scattering for 3-D random obstacles with multifrequency data, arXiv:2607.13473 (Tier B) | Estimates reference geometry plus Karhunen-Loeve/covariance statistics of random shape fluctuations; proves far-field probability law determines the radial function in distribution | Geometry is random per-sample, but sensors are fixed; no pose/measurement-geometry uncertainty |

## 2. Enabling components

### 2a. SOM / current-space spectral machinery (generic algebra substrate)

These papers define the singular-spectrum partition of the induced-current to
field operator. This is the closest retrieved algebraic precedent for a
"retention spectrum," but it is generic: fixed arrays, known geometry, no pose.

| Paper (Tier B) | Contribution |
|---|---|
| Chen 2010, Subspace-based optimization method for solving inverse-scattering problems, IEEE TGRS, DOI 10.1109/TGRS.2009.2025122 | Foundational SOM: part of the contrast source is fixed by spectrum analysis, the complementary part by optimization; flexible signal/noise subspace partition gives noise robustness |
| Zhong et al. 2011, An FFT twofold subspace-based optimization method, IEEE TAP, DOI 10.1109/TAP.2010.2103027 | Replaces singular vectors with discrete Fourier bases for the current subspace, avoiding SVD; 2-D TM and 3-D validation |
| Zhong et al. 2010, An improved SOM and its implementation in 3-D inverse problems, IEEE TGRS, DOI 10.1109/TGRS.2010.2049744 | Cheaper construction of the ambiguous part and cheaper SVD of the current-to-field map, enabling 3-D |
| Miao Wang et al. 2024, Cross-correlated SOM, IEEE TAP, DOI 10.1109/TAP.2024.3450328 | Combines cross-correlated cost with subspace optimization; multifrequency data and per-frequency L-curve regularization improve high-contrast inversion; robustness gains cost extra compute |
| Siampour et al. 2024, Imaging through non-radiating subspace, J. Opt. Photonics Res., DOI 10.47852/bonviewjopr42022785 | SQP inside SOM; analyzes singular values of the induced-current to exterior-field operator to weigh non-radiating subspace reconstruction; FEM-BI/CDM/EFIE forward models |
| Tingsen Zhang et al. 2025, Kernel subspace optimization with shape priors, IEEE Antennas Wireless Propag. Lett., DOI 10.1109/LAWP.2025.3594841 | SVD-based low-dimensional shape parameterization maps pixel space to a prior space; kernel modeling of the state-residual map; robust under high noise |

### 2b. Multifrequency / geometry / 3-D reconstruction machinery

| Paper (Tier B) | Contribution |
|---|---|
| Wei 2018, Effects of multiple scattering on resolution of full-wave inverse-scattering solver, ICEAA, DOI 10.1109/ICEAA.2018.8520511 | Quantifies resolution behavior vs scatterer strength and noise in full-wave reconstructions |
| Haynes et al. 2024, Born inversion with object T-matrix and full bistatic spherical geometry, IEEE TAP, DOI 10.1109/TAP.2023.3329645 | Semi-analytic 3-D Born reconstruction; Toeplitz normal equations; FFT-accelerated CG |
| Zimmerling et al. 2025, Lippmann-Schwinger-Lanczos algorithm with unknown reflectivity and loss, arXiv:2511.15058 | Data-driven reduced-order model approximates internal fields from spectral or frequency-domain measurements, beyond the Born background-field replacement |
| Yutong Du et al. 2026, Coordinate-residual physics-driven neural network, arXiv:2608.09382 (S2 record lacks year/abstract; OpenAlex lookup) | Reconstructs 3-D contrast from normalized spatial coordinates without preliminary region selection; 2.10% mean error vs 7.97% CSI; ~5.5-12.1x speedups |
| Youzi He et al. 2026, Review of inverse-scattering imaging via transmission eigenfunctions, Mathematics, DOI 10.3390/math14101586 | Surveys spectral/eigenfunction imaging; open problems include super-resolution explanation and limited aperture |
| Sendrea et al. 2025, Non-iterative approximations via eigenfunction expansions, ACES, DOI 10.23919/ACES66556.2025.11052550 | Eigencurrents of a reference geometry approximate shape and material parameters; Austria profile to ~1e-2 error |

### 2c. Statistical / identifiability theory

| Paper (Tier B) | Contribution |
|---|---|
| Furuya et al. 2024, Consistency of the Bayes method for the inverse scattering problem, Inverse Problems, DOI 10.1088/1361-6420/ad3089 | Posterior consistency with explicit, minimax-optimal logarithmic contraction rate for refractive-index recovery from far-field data |
| Hongyu Liu et al. 2024, Solving Schiffer's problem almost surely, Commun. Anal. Comput., DOI 10.3934/cac.2024014 | Probabilistic resolution of single-far-field-measurement shape uniqueness |
| Oksanen et al. 2026, Fixed-angle inverse scattering with non-constant velocity, arXiv:2608.13670 | Uniqueness of top-order coefficients from finitely many plane waves plus complementary solutions under pseudoconvexity/no-caustics/spanning conditions |

### 2d. Adjacent but distinct modalities (explicitly separated)

These retrieve "information retention vs singular spectrum" or
"background/nuisance" themes in different physics. They are not usable as
substantive support for the wave-SLAM combination.

| Paper (Tier B) | What it shows | Why it stays separated |
|---|---|---|
| Rambo et al. 2024/2025, Information-theory optimization of small-angle scattering, DOI 10.1016/j.bpj.2025.06.031 | Shannon number equals the maximum number of significant singular values and marks the well/ill-conditioned transition of a band-limited inverse transform | Solution-state SAXS/SANS of biomolecules, not wave imaging; no geometry/pose unknowns |
| Hohage, Novikov & Sivkin 2024, Phase retrieval and phaseless inverse scattering with background information, Inverse Problems, DOI 10.1088/1361-6420/ad6fc6 | A known background scatterer supplies the missing phase for recovering an unknown potential | "Background" is a fixed a priori scatterer, not a trajectory/nuisance variable |

## 3. Counterevidence / limitations (against optimistic feasibility claims)

| Paper (Tier B) | Limitation recorded |
|---|---|
| Wei 2018 (DOI 10.1109/ICEAA.2018.8520511) | Weak scatterers are inversion-unstable, especially under noise; strong scatterers are hard to recover even at low noise; moderate contrast is the stable regime |
| Furuya et al. 2024 (DOI 10.1088/1361-6420/ad3089) | Information about the map contracts only logarithmically in sample size — a fundamental, minimax-optimal statistical hardness for IS map estimation |
| Xu Zhang et al. 2025 (DOI 10.1109/TAP.2024.3524008) | Unique parameter inversion requires Nyquist-satisfying frequency sampling, suitable polarization, and adequate mono/bistatic coverage; robustness is bounded by SNR and data volume |
| Capozzoli et al. 2023 (DOI 10.46620/22-0065) | Aspect-limited acquisition degrades reconstruction — restricted/perturbed viewing geometry is a real cost, not a free parameter |
| Miao Wang et al. 2024 (DOI 10.1109/TAP.2024.3450328) | Enhanced robustness comes at higher computational complexity — a constraint on any trajectory-sensitive online method |

## 4. Explicit separation from the proposed combination

- **Generic EFIM / Schur algebra:** not present in the retrieved bundle at
  all (expected queries 429'd). Nothing here may be cited for
  nuisance-eliminated Fisher information.
- **SOM / current-space work:** present (Section 2a), but it is static,
  known-geometry algebra. It is a spectral substrate only, not pose-aware
  SLAM. Do not cite it as if it analyzed a pose-map gauge.
- **Unknown-sensor inverse scattering:** absent. The nearest retrieved
  analogies ([20] unknown impedance, [29] random geometry, [38] background
  scatterer) treat unknowns other than sensor pose and must be labeled as
  analogies.
- **Wave/RF SLAM, bilinear pose-map gauge, robust trajectory design:** absent
  from the bundle; none of the 40 abstracts mentions SLAM, pose, trajectory,
  Fisher/EFIM, principal angles, or canonical correlations (verified by an
  abstract-level keyword scan; the only "gauge" hit is the AdS/CFT holography
  paper [23], excluded).

## 5. Irrelevant exclusions

| # | Paper | Reason for exclusion |
|---|---|---|
| 2 | Enhanced CC-SOM (IWS 2024) | Near-duplicate of Miao Wang 2024 |
| 3 | Choice of subspace for QMR (arXiv:2508.05793) | Generic Krylov linear algebra; no scattering or geometry |
| 5 | Multiresolution SOM (JOSA A 2011) | SOM variant; near-duplicate |
| 6 | Phaseless-data SOM (IEEE TGRS 2011) | SOM variant; phaseless niche |
| 8 | SOM with inhomogeneous background (Inv. Probl. 2010) | SOM variant; near-duplicate |
| 9 | Regularization-parameter role in SOM (APMC 2009) | Niche methodological note |
| 12 | USCTNet HSI reconstruction (arXiv:2509.10651) | RGB-to-hyperspectral; different modality |
| 15, 16 | SOM for anisotropic laminates | Niche application duplicates |
| 17 | Distorted-Born subspace optimization (PIERS 2016) | SOM variant; near-duplicate |
| 18 | Quantum inverse algorithm (JPC A 2023) | Quantum computing; different field |
| 19 | Nanoscale NDE inverse scattering (2010) | Niche metrology; no abstract in bundle |
| 23 | Inverse scattering in holographic bulk reconstruction (JHEP) | AdS/CFT; different domain |
| 26 | T-matrix plant stem model (CAMA 2024) | Niche agricultural radar application |
| 33 | Shape optimization for diffusion/DOT (2025) | Diffusion equation, not full-wave scattering; no pose |
| 34 | Inverse Compton x-ray geometry optimization (2025) | Accelerator physics |
| 35 | One- and two-photon light in inverse scattering (2023) | Quantum-optics resolution limits; modality mismatch |
| 36 | Biharmonic wave direct imaging (2025) | Different PDE model |
| 39 | Temporal domain derivative in acoustic obstacle scattering (2025) | Time-domain obstacle shape reconstruction; no pose |

## 6. Coverage gaps and omitted steps

- 429 gaps: all snippet searches (6) and paper searches for 4 of 6 queries
  failed (recorded in `evidence.json` `operation_errors`). The candidate set
  therefore over-represents the two queries that succeeded
  (subspace/SOM-oriented and multifrequency-oriented).
- No snippet or full-text retrieval was possible, so Tier A is absent and the
  screening is abstract-level.
- No SPECTER/embedding step was recorded in the bundle.
- EFIM/Schur, wave/RF SLAM, bilinear-gauge, unknown-sensor, and robust-design
  literatures could not be screened; their absence is a retrieval gap, not
  evidence of novelty.
- Four SOM classics had publisher-elided abstracts in Semantic Scholar; their
  abstracts were recovered from OpenAlex records ([7], [4], [13], and the
  excluded [8]). [30]'s abstract and arXiv ID come from an OpenAlex
  title lookup because the S2 record lacks both.
