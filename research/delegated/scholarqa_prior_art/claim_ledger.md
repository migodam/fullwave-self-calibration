# Claim ledger: full-wave inverse-scattering SLAM pose-confounding thesis

Provenance notes that apply to every entry:

- All substantive evidence is **Tier B** (abstract-level). No snippet or
  full-text passage was retrieved (snippet searches all 429'd), so **Tier A
  is absent**. `scholarqa-cli verify` confirmed bibliographic identity only
  (21/21 requested citations resolved, 0 unresolved); it does not confirm
  claim support.
- Abstracts for the SOM classics ([7], [4], [13]) were recovered from
  OpenAlex because Semantic Scholar elides them at the publisher's request.
  [30]'s abstract and arXiv ID come from an OpenAlex title lookup (S2 record
  has no year/abstract). Everything else is from `evidence.json`.
- No entry below is marked "supports" the proposed pose-nuisance-eliminated
  effective-map-information combination: nothing retrieved analyzes that
  combination. Entries are contextualizing or qualifying only.

## C1. SOM partitions the induced-current space by spectrum analysis

- **Claim:** Subspace-based optimization splits the contrast source into a
  part fixed by spectrum analysis of the induced-current to scattered-field
  mapping and a complementary part recovered by optimization, with a flexible
  signal/noise partition that yields noise robustness and fast convergence.
- **Paper:** Chen 2010, S2 `16b84bd4dfb430f5b1518d7a2f58367149c45497`,
  DOI 10.1109/TGRS.2009.2025122.
- **Evidence (Tier B, abstract):** "part of the contrast source is determined
  from the spectrum analysis without using any optimization, whereas the rest
  is determined by optimization"; "great flexibility in partitioning the
  space of induced current into two orthogonal complementary subspaces: the
  signal subspace and the noise subspace"; robust against noise, fast
  convergence; shares properties with contrast-source inversion.
- **Relationship:** Contextualizes — a generic singular-spectrum partition of
  the full-wave measurement operator, with fixed geometry and no pose.
- **Study context:** Electromagnetic inverse scattering, numerical
  simulations; dataset and frequency details not stated in the abstract.
- **Confidence:** High — direct abstract statement from the foundational SOM
  paper.

## C2. Fourier bases replace the SVD subspace to scale SOM

- **Claim:** FFT/improved variants replace the singular-vector current
  subspace with discrete Fourier bases or cheaper ambiguous-part
  constructions, avoiding the dominant SVD cost and extending SOM to 3-D.
- **Papers:** Zhong et al. 2011, S2 `e5eff7c7a6a23de2571cded519dd3f38585cd9a7`,
  DOI 10.1109/TAP.2010.2103027; Zhong et al. 2010, S2
  `863518871fbd222749c198b61a73065597b6f692`, DOI 10.1109/TGRS.2010.2049744.
- **Evidence (Tier B, abstracts):** FFT-TSOM uses "discrete Fourier bases to
  construct a current subspace that is a good approximation to the original
  current subspace spanned by singular vectors," avoiding SVD; the improved
  SOM "reduces ... the computational complexity of the singular-value
  decomposition of the mapping from the induced current to scattered
  fields," enabling 3-D problems.
- **Relationship:** Contextualizes — scalability of the spectral-subspace
  substrate in full-wave inversion; no pose content.
- **Study context:** 2-D TM and 3-D electromagnetic inverse scattering,
  numerical validation.
- **Confidence:** High — direct abstract statements.

## C3. Multifrequency + cross-correlated cost improves high-contrast robustness

- **Claim:** Combining the cross-correlated cost function with subspace
  optimization, and using multifrequency data with per-frequency L-curve
  regularization, improves inversion of high-contrast scatterers at higher
  computational cost.
- **Paper:** Miao Wang et al. 2024, S2
  `c0bc431a56b37f41d6f9fc27a9a7335a4aa6fd30`, DOI 10.1109/TAP.2024.3450328.
- **Evidence (Tier B, abstract):** CC-SOM "incorporat[es] the subspace
  optimization strategy" into CC-CSI; "multifrequency data are used to
  improve the inversion performance of high-contrast scatterers, where the
  L-curve method is introduced to select the regularization parameters of
  each frequency point"; SVD/FFT speedups; "advantages ... come at the cost
  of higher computational complexity."
- **Relationship:** Contextualizes (multifrequency compensation is for
  scatterer contrast, not pose); qualifies (cost).
- **Study context:** Synthetic and experimental electromagnetic inversion;
  high-contrast scatterers.
- **Confidence:** High — direct abstract statements.

## C4. The non-radiating subspace is governed by the operator's singular values

- **Claim:** In SOM, reconstruction quality through the non-radiating
  subspace is analyzed via the singular values of the induced-current to
  exterior-field mapping, and SQP minimization of radiating and non-radiating
  objectives improves imaging in complex background media.
- **Paper:** Siampour et al. 2024, S2
  `0ab5f694342580d436815deb0331b29bc517b955`,
  DOI 10.47852/bonviewjopr42022785.
- **Evidence (Tier B, abstract):** "We investigate the influence of
  non-radiating (NR) subspace reconstruction on imaging quality by analyzing
  the induced current-exterior field mapping operator's singular values";
  SQP minimizes radiating and non-radiating objectives; forward models
  include coupled dipole, FEM-boundary integral, and EFIE.
- **Relationship:** Contextualizes — singular-spectrum weighting of
  information-bearing vs invisible subspaces, still with fixed geometry.
- **Study context:** Numerical experiments for buried/complex dielectric
  objects; applications claimed in medical imaging, remote sensing, NDE.
- **Confidence:** High — direct abstract statements.

## C5. Shape-prior kernel SOM maps pixel space to a prior subspace

- **Claim:** A kernel SOM with SVD-based low-dimensional shape parameters
  maps the high-dimensional pixel space to a structured prior space,
  enforcing physical consistency and improving robust reconstruction under
  high noise and complex structures.
- **Paper:** Tingsen Zhang et al. 2025, S2
  `a1c65e723b23b9a33b8ddb158760f1dbabd2386e`, DOI 10.1109/LAWP.2025.3594841.
- **Evidence (Tier B, abstract):** "Leveraging singular value decomposition,
  this method employs low-dimensional shape parameters ... effectively
  mapping the high-dimensional pixel space to a structured prior space";
  kernel methods model the state-residual mapping; "two-timescale stochastic
  optimization"; "superior reconstruction performance under high-noise and
  complex structural scenarios."
- **Relationship:** Contextualizes — a subspace/prior compression of the map
  representation; no pose.
- **Study context:** Electromagnetic full-wave inverse scattering imaging
  (EFWISI), numerical experiments.
- **Confidence:** High — direct abstract statements.

## C6. Full-wave resolution depends strongly on scatterer strength

- **Claim:** In full-wave inverse-scattering reconstructions, weak scatterers
  are unstable especially under noise, strong scatterers are hard to recover
  even at low noise, and moderate contrast is the stable regime.
- **Paper:** Wei 2018, S2 `8f08f1a299b9fdab25e7cf3e62a3c9ab1c363b9c`,
  DOI 10.1109/ICEAA.2018.8520511.
- **Evidence (Tier B, abstract):** "for weak scatterers, inversion is
  unstable, especially when noise is high. For strong scatterers, one can
  hardly obtain a solution that is close to the exact one even when the
  noise level is low"; moderate scatterers are more stable; twofold SOM used
  in experiments.
- **Relationship:** Qualifies — bounds any feasibility claim for the thesis's
  full-wave map estimate across contrast regimes.
- **Study context:** Permittivity reconstruction from scattered fields for
  weak/moderate/strong scatterers under different noise levels.
- **Confidence:** High — direct abstract statements.

## C7. Sensor locations and frequency sampling can be singular-value optimized

- **Claim:** Under the Born approximation, transmitter/receiver locations and
  the frequency sampling step can be chosen by singular-value optimization,
  and aspect-limited acquisitions degrade the achievable performance.
- **Paper:** Capozzoli et al. 2023, S2
  `33f8e5c7ee0d77321c157bd1dfb222731e3b67d7`, DOI 10.46620/22-0065.
- **Evidence (Tier B, abstract):** "We optimize the transmitter and receiver
  locations, as well as the frequency sampling step, and we numerically study
  the performance of the method for aspect-limited acquisitions," over
  circular 2-D domains with multifrequency multimonostatic data.
- **Relationship:** Contextualizes (design-time geometry optimization, the
  closest retrieved geometry/information tradeoff); qualifies (limited
  aperture costs information).
- **Study context:** 2-D Born-approximation inverse scattering, numerical
  study.
- **Confidence:** High — direct abstract statement; the relevance to SLAM is
  analyst interpretation (no pose estimation in the paper).

## C8. Shape and an unknown impedance are jointly recoverable across frequencies

- **Claim:** Recursive linearization across frequencies recovers both the
  shape and the impedance function of an obstacle, and shape recovery remains
  accurate even under mismatched sound-hard/sound-soft boundary conditions.
- **Paper:** Borges et al. 2021, S2
  `09cc53da44b9181b3c8154b23cc8402aaacbe300`, DOI 10.1007/s10444-021-09915-1.
- **Evidence (Tier B, abstract):** RLA is "a continuation method in
  frequency"; recovers "the shape and impedance functions of the object";
  "one can recover the shape with high accuracy even when the measurements
  are generated by sound-hard or sound-soft objects"; convergence behavior
  questions remain open.
- **Relationship:** Contextualizes — joint recovery of a scene property and
  an unknown model parameter, and multi-frequency robustness; the nuisance is
  a boundary impedance, not a sensor pose.
- **Study context:** Acoustic obstacle scattering, plane waves from multiple
  directions and frequencies, numerical examples.
- **Confidence:** High — direct abstract statements.

## C9. Unique surface-parameter inversion requires specific observation conditions

- **Claim:** Unique inversion of geometry/position parameters of electrically
  large curved surfaces requires frequency sampling satisfying Nyquist,
  suitable polarization and mono/bistatic observation schemes, and robustness
  is limited by SNR and data volume.
- **Paper:** Xu Zhang et al. 2025, S2
  `b2d77880f9e2172666e6157f8322824b855365de`, DOI 10.1109/TAP.2024.3524008.
- **Evidence (Tier B, abstract):** "determines the necessary conditions of
  observation for obtaining a unique solution in parameter inversion";
  singly curved lengths invertible "when frequency points satisfy the
  Nyquist sampling theorem"; algebraic linearization of curvature/position
  terms; measurement schemes for monostatic/bistatic modes; "relationship
  between the robustness of the method and ... SNR and the amount of radar
  data."
- **Relationship:** Qualifies — observability constraints relevant to any
  wave-based geometry inference; the inverted position parameters belong to
  the scatterer, not the sensor.
- **Study context:** Radar parameter inversion of general curved surfaces,
  numerical simulations.
- **Confidence:** High — direct abstract statements.

## C10. Random obstacle geometry and its statistics are recoverable

- **Claim:** For 3-D random star-shaped obstacles, a Monte-Carlo
  multifrequency recursive-linearization scheme recovers the reference
  geometry and fluctuation statistics (Karhunen-Loeve eigenvalues,
  covariance hyperparameters, marginal distributions), and the far-field
  probability law uniquely determines the radial function in distribution.
- **Paper:** Zhiqi Sun et al. 2026, S2
  `4babffad2e2855db27b0bdde500c54cab5859eab`, arXiv:2607.13473.
- **Evidence (Tier B, abstract):** estimates "the reference geometry and key
  statistics of the shape fluctuation field including Karhunen-Loeve
  eigenvalues, covariance hyper-parameters ... representative marginal
  distributions"; "the probability law of the far-field data uniquely
  determines the radial function in distribution."
- **Relationship:** Contextualizes — statistical treatment of geometry
  uncertainty, but for random *object* shape with fixed sensors, not pose.
- **Study context:** 3-D acoustic inverse scattering with Gaussian and
  non-Gaussian shape fluctuations; Monte-Carlo numerical experiments.
- **Confidence:** High for the stated result; preprint (arXiv-only) not yet
  peer-reviewed, which tempers strength.

## C11. Single far-field measurement determines an obstacle almost surely

- **Claim:** Schiffer's problem (shape uniqueness from one far-field
  measurement) holds in more general settings in a probabilistic
  ("almost surely") sense.
- **Paper:** Hongyu Liu et al. 2024, S2
  `c08f4f1745f46afda7a7eb4ef1ec78f01ae7144f`, DOI 10.3934/cac.2024014.
- **Evidence (Tier B, abstract):** "We show that this conjecture holds true
  in more general settings in the probability sense," with implications for
  practice and broader inverse problems.
- **Relationship:** Contextualizes — a uniqueness/identifiability benchmark;
  says nothing about pose-map ambiguity or multi-pose trajectories.
- **Study context:** Inverse obstacle scattering theory; a short note with
  probabilistic perspective.
- **Confidence:** High — direct abstract statement.

## C12. Bayesian map estimation contracts only logarithmically

- **Claim:** A nonparametric Bayesian posterior for the refractive index from
  noisy far-field measurements is consistent with an explicit logarithmic
  contraction rate, and this rate is optimal in the minimax sense.
- **Paper:** Furuya et al. 2024, S2
  `8b9ad0c5d2a6ffd97aac30dea714fd9458efe65d`, DOI 10.1088/1361-6420/ad3089.
- **Evidence (Tier B, abstract):** "establish the consistency of the
  posterior distribution with an explicit contraction rate in terms of the
  sample size"; "the contraction rate is of a logarithmic type"; "such
  contraction rate is optimal in the statistical minimax sense"; two priors
  considered.
- **Relationship:** Qualifies — a fundamental, minimax-optimal statistical
  hardness for map information from far-field data (no pose nuisance in the
  model).
- **Study context:** Far-field inverse scattering; frequentist analysis of
  the Bayesian posterior over large noisy samples.
- **Confidence:** High — direct abstract statements.

## C13. Transmission eigenfunctions encode qualitative and quantitative scatterer information

- **Claim:** Transmission eigenfunctions carry intrinsic qualitative and
  quantitative information about scatterers and support imaging algorithms;
  open problems include explaining super-resolution and handling
  limited-aperture data.
- **Paper:** Youzi He et al. 2026, S2
  `0dac64ff7d2f2b0ed4d410ae6a514344da890ad8`, DOI 10.3390/math14101586.
- **Evidence (Tier B, abstract):** "the transmission eigenfunctions contain
  important qualitative and quantitative information about the unknown
  scatterers"; review covers "theoretical justification ... comparisons ...
  with traditional qualitative reconstruction methods"; open problems:
  "theoretical explanation of super-resolution effect, limited-aperture
  data, and the extension to more complex physical models."
- **Relationship:** Contextualizes — spectral/eigenfunction view of what the
  scattering operator retains; no pose.
- **Study context:** Review article; transmission-eigenfunction imaging
  literature.
- **Confidence:** Medium — review-level claims; no single study's data.

## C14. Reference-geometry eigencurrents give non-iterative reconstruction

- **Claim:** Eigencurrents of a reference geometry non-iteratively
  approximate an unknown scatterer's shape and electrical properties,
  reconstructing the Austria profile with errors on the order of 1e-2.
- **Paper:** Sendrea et al. 2025, S2
  `6d09c4e2fdc5f45e5650940c2db802b9db1185dd`,
  DOI 10.23919/ACES66556.2025.11052550.
- **Evidence (Tier B, abstract):** "utilizes eigencurrents of a reference
  geometry to approximate the shape and electrical properties of an unknown
  scatterer"; "validated by reconstructing the 'Austria' profile ... with
  errors on the order of 10^-2."
- **Relationship:** Contextualizes — spectral-basis reconstruction; no pose.
- **Study context:** Electromagnetic inverse scattering; Austria benchmark.
- **Confidence:** Medium — single-profile demonstration, conference-length
  record.

## C15. A data-driven ROM replaces the Born internal-field approximation

- **Claim:** The Lippmann-Schwinger-Lanczos reduced-order model, built from
  measured spectral or frequency-domain data, approximates internal fields
  and recovers unknown reflectivity and loss, improving accuracy and
  robustness versus the Born background-field replacement.
- **Paper:** Zimmerling et al. 2025, S2
  `5c9fcb23ce10e6faf767e39624751cb8a01868b9`, arXiv:2511.15058.
- **Evidence (Tier B, abstract):** "approximating the internal solutions
  through the lifting of states from a reduced-order model constructed
  directly from the measured data"; constructions "based on spectral data
  and another on frequency-domain measurements"; "Compared to the Born
  approximation ... our approach yields more accurate internal
  reconstructions and enables faster and more robust recovery of the
  contrast."
- **Relationship:** Contextualizes — data-driven internal-field machinery for
  unknown media; 1-D, no pose.
- **Study context:** One-dimensional attenuating media with unknown
  reflectivity and loss; numerical experiments.
- **Confidence:** Medium — preprint (arXiv-only), 1-D case.

## C16. Full-sphere quad-pol Born inversion has Toeplitz structure

- **Claim:** With a known object T-matrix and quad-pol bistatic measurements
  over a sphere, Born-approximation dielectric inversion reduces to a
  semi-analytic deconvolution whose normal matrix is analytic and Toeplitz,
  enabling FFT-accelerated conjugate-gradient solution.
- **Paper:** Haynes et al. 2024, S2
  `c4eafff1cf52243494cd66de98b4e6d1b75f2989`, DOI 10.1109/TAP.2023.3329645.
- **Evidence (Tier B, abstract):** "the dielectric contrast is computed by
  deconvolving the T-matrix backprojected image with the vector bistatic
  point target response"; "The normal matrix is given analytically and has
  Toeplitz structure that facilitates fast Fourier transform
  (FFT)-accelerated conjugate gradient solutions."
- **Relationship:** Contextualizes — structure of the 3-D Born measurement
  operator under full spherical sampling; fixed geometry, no pose.
- **Study context:** Electromagnetic inversion of electrically large
  dielectric collections; numerical validation.
- **Confidence:** High — direct abstract statements.

## C17. Coordinate-parameterized networks reconstruct 3-D contrast directly

- **Claim:** A physics-driven network that takes normalized spatial
  coordinates with a residual convolutional backbone reconstructs 3-D
  contrast without a preliminary reconstruction, achieving 2.10% mean error
  versus 7.97% (CSI) and 3.99% (L2/3-FBE-WCIE) and ~5.5-12.1x speedups.
- **Paper:** Yutong Du et al. 2026, S2
  `6269e09f30f4de6ed4ecb19f6ecf6dc4c2399a72`, arXiv:2608.09382 (OpenAlex).
- **Evidence (Tier B, abstract):** "directly reconstructs the unknown
  contrast distribution using normalized spatial coordinates and a residual
  convolutional network, without requiring a preliminary reconstruction";
  noise-free 3-D cases: "average relative error of 2.10%, compared with
  7.97% for CSI and 3.99% for L2/3-FBE-WCIE"; "approximately 5.5- and
  12.1-fold speedups"; 3-D Fresnel experiments.
- **Relationship:** Contextualizes — coordinate-parameterized 3-D full-wave
  reconstruction; "coordinates" are scene coordinates, not sensor pose.
- **Study context:** 3-D electromagnetic inverse scattering; synthetic and
  Fresnel experiments.
- **Confidence:** High for the stated benchmark; relevance to pose is
  analyst interpretation and is limited.

## C18. Shannon number equals significant singular values (distinct modality)

- **Claim:** In solution-state small-angle scattering, the Shannon number of
  a band-limited signal equals the maximum number of significant singular
  values and marks the well/ill-conditioned transition, allowing direct
  Shannon-limited inversion under detector oversampling.
- **Paper:** Rambo et al. 2024/2025, S2
  `8c64d4be9472bb988211915523f5146b58bb0917`, DOI 10.1016/j.bpj.2025.06.031.
- **Evidence (Tier B, abstract):** "the ill-conditioning of the inverse
  problem is directly related to the Shannon number"; "The limit corresponds
  to the maximum number of significant singular values that can be recovered
  in a SAXS experiment"; oversampling "enables a direct inverse Fourier
  transform"; hybrid AIC/Durbin-Watson scoring.
- **Relationship:** Contextualizes only — an information-retention vs
  singular-spectrum correspondence in a different modality (SAXS/SANS);
  explicitly separated, not usable as substantive support for wave SLAM.
- **Study context:** Solution-state X-ray/neutron small-angle scattering of
  biological macromolecules.
- **Confidence:** High for the modality; low relevance transfer to wave SLAM.

## C19. Fixed-angle uniqueness under variable sound speed

- **Claim:** With fixed-angle boundary data from finitely many plane waves
  and their complementary solutions, the highest-order coefficients of a
  second-order hyperbolic operator are uniquely recoverable under
  pseudoconvexity, no-caustics, and spanning conditions.
- **Paper:** Oksanen et al. 2026, S2
  `9fb62b8817283e8bb5bcbd55be81ac9025a437aa`, arXiv:2608.13670.
- **Evidence (Tier B, abstract):** "by measuring the boundary data of
  finitely many plane waves and their complementary solutions, one can
  uniquely recover the unknown coefficients of the highest order terms";
  requires pseudoconvexity, no-caustics, and spanning; Carleman estimates.
- **Relationship:** Contextualizes — formal identifiability for wave
  operators with non-constant velocity; no pose or SLAM.
- **Study context:** Mathematics of formally determined inverse problems for
  wave equations.
- **Confidence:** Medium — preprint (arXiv-only), pure-theory result.

## C20. Known background scatterers supply missing phase (distinct nuisance)

- **Claim:** Using an a priori known, sufficiently disjoint background
  scatterer, an iterative scheme recovers a compactly supported potential
  from a single phaseless differential cross section; relaxing disjointness
  requires additional monochromatic data.
- **Paper:** Hohage, Novikov & Sivkin 2024, S2
  `7dab2b76dfdb113d479e3acb1faf92dee85e3cb7`, DOI 10.1088/1361-6420/ad6fc6.
- **Evidence (Tier B, abstract):** "we propose an iterative scheme for
  finding the potential from measurements of a single differential
  scattering cross section corresponding to the sum of the unknown potential
  and a known background potential"; if disjointness is relaxed, "similar
  results ... from additional monochromatic measurements."
- **Relationship:** Contextualizes — a "known background" compensates an
  unknown; the background is a fixed a priori scatterer, not a trajectory.
- **Study context:** Multidimensional Schrodinger equation, fixed energy;
  numerical examples.
- **Confidence:** High — direct abstract statements.

## C21. Coverage claim: no retrieved record addresses the pose-SLAM combination

- **Claim:** Among the 40 candidate records in `evidence.json`, no abstract
  mentions SLAM, pose estimation, trajectory/sensor-position estimation,
  Fisher/EFIM, principal angles, canonical correlations, or a pose-map gauge
  structure; none analyzes pose-nuisance-eliminated map information.
- **Paper:** Negative finding across all 40 candidate IDs listed in
  `evidence.json` `candidate_paper_ids` (abstract-level keyword and content
  scan; the single "gauge" hit is the AdS/CFT paper [23], excluded).
- **Evidence (Tier B):** abstract-level scan of the full bundle.
- **Relationship:** Defines the gap — this is a bounded, non-exhaustive
  retrieval gap, not a global novelty claim.
- **Study context:** The two queries that succeeded were
  subspace/optimization-oriented and multifrequency-oriented, biasing the
  candidate set away from SLAM/EFIM literatures.
- **Confidence:** High for the retrieved 40 records; explicitly cannot be
  extended beyond the bundle because four search queries and all snippet
  searches failed with 429.
