# Targeted Prior-Art Addendum

Date checked: 2026-09-03. Scope: targeted discovery of the literatures that the
partial ScholarQA collection missed because of HTTP 429 responses. Evidence is
publisher/arXiv abstract or publisher metadata unless explicitly noted; it is
not an exhaustive search and cannot support a global novelty claim.

## Closest neighboring lines

1. **RF-SLAM map-information bounds.** Deutschmann, Li, Meyer, and Leitinger,
   *Posterior Cramér--Rao Bounds on Localization and Mapping Errors in
   Distributed MIMO SLAM*, Asilomar 2025,
   DOI 10.1109/IEEECONF67917.2025.11443670 / arXiv:2506.19957. The work derives a
   posterior mapping error bound for positions and orientations of specular
   surfaces jointly with mobile transceiver localization. This is the closest
   retrieved neighbor on map-vs-pose information. Its map is a finite
   geometric/specular-landmark model with parametric single/double-bounce
   paths, not a volumetric contrast field governed by the Lippmann--Schwinger
   equation; the abstract does not claim principal-angle retention or a
   low-rank pose defect.
   Source: https://arxiv.org/abs/2506.19957

2. **Radar autofocus with physical antenna-position ambiguity.** Mansour,
   Liu, Kamilov, and Boufounos, *Sparse Blind Deconvolution for Distributed
   Radar Autofocus Imaging*, IEEE TCI 4(4), 2018,
   DOI 10.1109/TCI.2018.2875375 / arXiv:1805.03269. It explicitly treats
   antenna position ambiguity on mobile radar platforms and maps each position
   error to an image-domain spatial-shift kernel, giving a multichannel blind
   deconvolution problem. This is a close physical and bilinear antecedent. It
   assumes a static radar image with sparse/piecewise-smooth priors rather than
   full-wave multiple-scattering volumetric contrast, and it does not provide
   the proposed Fisher/principal-angle spectrum.
   Sources: https://arxiv.org/abs/1805.03269 and
   https://www.merl.com/publications/TR2018-055

3. **Joint SAR image and phase-error estimation.** Önhon and Çetin,
   *A Sparsity-Driven Approach for Joint SAR Imaging and Phase Error
   Correction*, IEEE TIP 21(4), 2012, DOI 10.1109/TIP.2011.2179056, jointly
   alternates image formation and motion-induced phase-error correction.
   Scarnati and Gelb, *Joint image formation and two-dimensional autofocusing
   for synthetic aperture radar data*, JCP 374, 2018,
   DOI 10.1016/j.jcp.2018.07.059, jointly estimates a frequency- and
   azimuth-dependent phase correction and a piecewise-smooth image. These are
   direct antecedents to joint map/model-error estimation, but their nuisance
   variables are phase-error fields under SAR image models, not shared SE(2)
   poses in a nonlinear full-wave contrast-source model.
   Sources: https://research.sabanciuniv.edu/id/eprint/19013/ and
   https://doi.org/10.1016/j.jcp.2018.07.059

4. **Generic bilinear identifiability up to transformations.** Li, Lee, and
   Bresler, *A Unified Framework for Identifiability Analysis in Bilinear
   Inverse Problems with Applications to Subspace and Sparsity Models*,
   arXiv:1501.06120 (later IEEE TIT, DOI 10.1109/TIT.2016.2637933), supplies
   necessary/sufficient identifiability conditions up to transformation groups
   for structured bilinear inverse problems such as blind gain/phase
   calibration. It is an enabling theorem family, not a result specialized to
   Helmholtz pose tangents or Schur-eliminated map information.
   Source: https://arxiv.org/abs/1501.06120

5. **Radio-SLAM algorithms and geometric multipath maps.** The 2024--2025
   mmWave radio-SLAM literature jointly estimates user state and geometric
   landmarks and includes CRLB-based validation, but uses channel-parameter or
   point/surface representations rather than full-wave volumetric material
   inversion. Representative primary sources are
   DOI 10.1109/TCOMM.2024.3393977 (multipath identification, localization, and
   environment mapping) and DOI 10.1109/JSAC.2024.3413995 (end-to-end mmWave
   radio SLAM).
   Sources: https://ieeexplore.ieee.org/document/10509551/ and
   https://ieeexplore.ieee.org/document/10556695/

6. **FIM/Jacobian observability for SLAM-based sensor calibration.** Su, Kong,
   Sukkarieh, and Huang, *Necessary and Sufficient Conditions for Observability
   of SLAM-Based TDOA Sensor Array Calibration and Source Localization*, IEEE
   T-RO 37(5), 2021, DOI 10.1109/TRO.2021.3069140, derives necessary and
   sufficient full-column-rank conditions using the equivalence of FIM and
   Jacobian rank, including impossible motion/array configurations. This is a
   close methodological antecedent for rank-based joint map/nuisance
   observability. Its forward model is parametric TDOA with microphone geometry,
   offsets, and clock drift, not coherent Helmholtz scattering from a
   distributed material map; it does not study the retained-information
   spectrum after nuisance elimination.
   Source: https://opus.lib.uts.edu.au/handle/10453/157465

7. **Simultaneous inverse-scattering reconstruction and transmitter
   localization.** Karthik and Ghosh, *A Scalable Deep Learning Model for
   Simultaneous Reconstruction and Transmitter Localization in Inverse
   Scattering*, PIERS 2023, DOI 10.1109/PIERS59004.2023.10221374, is a mandatory
   direct-title neighbor. The publisher abstract and accessible introduction
   confirm that it jointly reconstructs contrast and localizes transmitters.
   Its stated multiple-incidence method processes each incidence independently
   with a scalable DNN and reduces a localization grid search from exponential
   to linear complexity. The accessible text does **not** establish that it
   derives a pose-eliminated Fisher operator, principal-angle retention law, or
   a finite-rank information defect. The remainder of the paper was behind a
   sign-in gate in this check, so its exact forward discretization, handling of
   multiple scattering, and any later theoretical analysis remain unverified.
   No manuscript may imply that jointly estimating an inverse-scattering map
   and transmitter location is new. The same authors' earlier titles on
   calibration-free single-incidence inversion
   (DOI 10.1109/PIERS53385.2021.9694867) and arbitrary transmitter
   configurations are part of the same mandatory follow-up cluster.
   Sources: https://ieeexplore.ieee.org/document/10221374,
   https://prague2023.piers.org/session.html?sid=S072, and
   https://doi.org/10.1109/PIERS59004.2023.10221374

## Bounded novelty consequence

The defensible contribution is not “joint imaging with unknown motion,” “Fisher
information for SLAM,” “bilinear identifiability,” or “subspace optimization”
individually. Each has clear prior art. A potentially differentiating synthesis,
still requiring full-text verification, is the explicit **whitened/realified
map-tangent versus physical-pose-tangent geometry for a nonlinear full-wave
volumetric contrast model**, with (i) pose-eliminated information expressed as
a finite-rank defect, (ii) no-prior generalized retention identified with
principal angles, (iii) Born-empty-background separated as a second-order
bilinear operator-uncertainty regime, and (iv) frequency/trajectory/rank-event
predictions tested in one reproducible solver.

Paper language must say “we formulate/test” and describe the search boundary;
it must not say “first,” “no prior work,” or imply exhaustive coverage.
