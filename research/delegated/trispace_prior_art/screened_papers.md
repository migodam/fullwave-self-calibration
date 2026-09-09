# Screened anchors: TriSpace SOM self-calibration

Date: 2026-09-04. Scope: bounded evidence cleanup for the 16 DOI/arXiv
anchors supplied by the requester; no broad re-search was run. Targeted
lookups were used only to recover abstract/landing text for anchors without
local descriptions.

Evidence tiers used here:

- `metadata` — bibliographic identity only (no abstract inspected).
- `abstract` — abstract/landing/intro/publisher text inspected (in this run
  or in the local addendum / ScholarQA bundle).
- `full-text` — article text actually inspected.
- `local-source` — the source documents or their digest characterize the
  paper in context, but this is a characterization, not full-text inspection
  of the paper.

No anchor was full-text inspected in this task. "Retrieval-bounded" applies
throughout: nothing below supports a global novelty statement.

## Focus object for the overlap columns

Unknown/unstable antenna geometry changes the sensing Green operator
`G_S(X)`; the claimed object is a canonical pose-equivalent induced-current
lift `delta j_eq` inside a **retained SOM/TSOM current space**, jointly using
(a) an irreducible data residual (data the restricted map cannot explain) and
(b) a state-equation inconsistency of the lifted current, for
self-calibration. A bare statement that geometry error is representable by a
current is treated as vacuous when the unrestricted sensing map has full row
rank, so overlap is assessed against the restricted lift, not the general
equivalence statement.

| # | Anchor | Mechanism (retrieved) | Overlap with focus | Decisive difference | Evidence |
|---|---|---|---|---|---|
| 1 | Chen, SOM, `10.1109/TGRS.2009.2025122` | Fixed-geometry contrast-source inversion: split the induced current `J` into a deterministic part fixed by the singular spectrum of the current-to-data operator `G_S` and a complementary part found by optimization. | Supplies the current-space `V_S^±` partition and retained current subspace that the claimed canonical lift must live in. SOM subspace reduction is antecedent machinery. | Geometry is known and fixed; `G_S` has no pose argument; no pose-equivalent lift, no residual pair used for self-calibration. | `local-source` + `metadata` (ScholarQA local verification); abstract-level characterization; no full text in this task. |
| 2 | Zhong and Chen, TSOM, `10.1088/0266-5611/25/8/085003` | Twofold SOM: applies a second current-space reduction through the internal/domain operator `G_D` after the `G_S` reduction. | Gives the `G_D`/`V_D` layer that defines state-operator structure and current-space state consistency for fixed geometry. | No unknown geometry; `V_S`, `V_D` are pose-independent; no `V_P`-type lift or calibration residual. | `local-source` (source-A ledger) + `metadata`; abstract-level; no full text in this task. |
| 3 | Hanabusa et al., CSI calibration, `10.1109/LGRS.2022.3169799` | Deep-learning data-calibration of measured frequency-domain responses before quantitative CSI for microwave subsurface imaging. | Establishes "calibrate data before contrast-source inversion" as prior art; calibration is explicitly in the CSI pipeline. | Calibration is a learned data-domain transformation, not a canonical current-space lift of a pose-dependent `G_S`; unknown physical geometry is not the calibrated variable. | `abstract` (publisher/repository text); `metadata`. |
| 4 | Cathers et al., microwave imaging calibration, `10.1109/OJAP.2023.3329356` | Calibration of an electromagnetic imaging system using 2-port error models of antennas and field propagation, replacing scalar calibration. | Systematic measurement-chain calibration for EM imaging is prior art; calibration accuracy affects subsequent inversion quality. | Hardware/system calibration performed on the measurement chain with known array geometry; not joint self-calibration during full-wave contrast inversion and not a pose-lift construction. | `abstract` (open-access text and catalog snippets); `metadata`. |
| 5 | Huang, Nammour, Symes, source-receiver extension FWI, `10.1190/geo2016-0301.1` | Full-waveform inversion with model extension: adds source and receiver extension degrees of freedom (extended modeling) and a matching objective so FWI tolerates kinematically inaccurate initial models. | Artificial source/receiver extension absorbs data-model mismatch — the same conceptual genus as adding geometry DOFs; receiver extension is specifically prior art. | Seismic wave-equation FWI, extension in source/receiver coordinates; no contrast-source current-space SOM/TSOM retention, no equivalent-current lift, no joint irreducible-data-plus-state residual pair. | `abstract` (landing/citation text); no full text in this task. |
| 6 | Metivier and Brossier, receiver-extension FWI, `10.1190/geo2020-0922.1` | Time-domain FWI receiver-extension ("relocalization"): receiver locations are introduced as artificial degrees of freedom, allowing FWI from crude initial models. | Direct antecedent for receiver-position unknowns absorbing data misfit and for treating geometry as an artificial extended DOF. | Receiver relocation is an artificial data-side unknown in acoustic/seismic FWI; no full-wave EM contrast source, no retained current subspace, no state-equation inconsistency of a lifted current used for calibration. | `abstract` (HAL/repository); no full text in this task. |
| 7 | da Silva et al., receiver-coordinate inaccuracies, `10.3997/2214-4609.2023101497` | Applies the receiver-extension FWI strategy to suppress 4D (time-lapse) noise caused by receiver-coordinate inaccuracies. | Receiver-coordinate error as a physical nuisance and receiver-extension correction is prior art in FWI. | Time-lapse seismic setting and empirical portability study; no volumetric EM/current-space formulation or residual-pair self-calibration. | `abstract` (EAGE proceedings text); `metadata`. |
| 8 | Karthik and Ghosh, simultaneous IS reconstruction + transmitter localization, `10.1109/PIERS59004.2023.10221374` | Scalable DNN for simultaneous contrast reconstruction and transmitter localization in inverse scattering. | Direct title-level and abstract-level neighbor: joint scene reconstruction plus unknown transmitter geometry is prior art. | Data-driven, per-incidence DNN with a localization search; the accessible text does not establish a retained SOM/TSOM current-space lift or a pose-defect residual analysis; remainder of paper behind a gate in the earlier check. | `abstract` + accessible intro (local addendum); full text not inspected. |
| 9 | Li, Lee, Bresler, bilinear blind calibration, `10.1109/TIT.2016.2637933` | Unified identifiability theory for bilinear inverse problems up to transformation groups, covering blind gain/phase calibration under subspace/sparsity models. | Bilinear structure and identifiability-up-to-gauge is prior art and applies to joint contrast/geometry problems; identifiability claims must respect this theory. | General abstract framework, not specialized to Helmholtz/Lippmann-Schwinger operators, pose-equivalent currents, or residual/state-defect coupling. | `local-source` (source-A ledger) + `metadata`; no full text in this task. |
| 10 | Mansour et al., radar autofocus, `10.1109/TCI.2018.2875375` | Sparse blind deconvolution for distributed radar autofocus: each antenna-position error is mapped to an image-domain spatial-shift kernel, yielding multichannel blind deconvolution. | Joint image formation with physical antenna-position ambiguity is prior art; mapping geometry error into an equivalent imaging-domain operator is a close conceptual antecedent to a lift. | Static radar image model with sparse/piecewise-smooth priors, linearized imaging; no nonlinear full-wave volumetric contrast, no current space, no state equation. | `abstract` (local addendum); no full text in this task. |
| 11 | Onhon and Cetin, joint SAR imaging and phase error, `10.1109/TIP.2011.2179056` | Sparsity-driven joint SAR image formation and phase-error correction (motion-induced phase errors). | Joint image + nuisance phase/geometry-model correction is prior art. | SAR image model and phase-error nuisance; not a full-wave contrast-source model, not current-space, no retained subspace/residual-pair geometry. | `abstract` (local addendum); no full text in this task. |
| 12 | Scarnati and Gelb, joint SAR imaging/autofocus, `10.1016/j.jcp.2018.07.059` | Joint image formation and two-dimensional autofocusing with frequency/azimuth-dependent phase correction. | Same genus as #10/#11: joint scene and model-error estimation is prior art. | SAR phase-error parameterization; no full-wave volumetric contrast or current-space SOM structure. | `abstract` (local addendum); no full text in this task. |
| 13 | Zhang et al., virtual antennas with SOM, arXiv `2312.17504` / `10.1109/TMTT.2024.3385996` | Frequency-domain zero-padding interpolation synthesizes extra (virtual) antenna samples with non-redundant scattered-field information in noise; MBA (linear) and SOM (nonlinear) then image the scatterers. | Uses the SOM current-space machinery with added synthetic measurement positions; shows that adding antenna samples can improve imaging. | Virtual antennas are interpolated data at known positions — aperture/data augmentation, not estimation of unknown physical geometry; no pose-dependent `G_S` or calibration residual. | `abstract` (arXiv and publisher text); `metadata`. |
| 14 | Sun, Kooij, Yarovoy, multifrequency CC-CSI, `10.1029/2017RS006505` | Cross-correlated CSI (CC-CSI) extended to multifrequency processing (MF-CC-CSI) with finite-difference frequency-domain forward modeling; cost function cross-correlates data-equation and state-equation errors. | Contrast-source framework jointly weighting data and state residuals is prior art; the idea that both residuals matter for current-space inversion is present in the CSI family. | Known geometry; calibration variables absent; the cross-correlated cost is an inversion criterion, not a residual pair discriminating a pose-lift from retained-current ambiguity. | `abstract` (publisher/repository text); no full text in this task. |
| 15 | Bevacqua et al., virtual experiments / equivalent contrast sources, `10.1109/TAP.2014.2382114` | Designs synthetic ("virtual") experiments by recombining real data and assumes focused contrast sources, recasting nonlinear inverse scattering algebraically (closed-form). | Virtual experiments and contrast-source manipulation/recombination in current space are prior art; equivalence-type constructions do not by themselves imply geometry calibration novelty. | Known array/scatterer geometry; virtual experiments are a solver/data-recombination device for the contrast, not a pose-equivalent lift for unknown `G_S`. | `abstract` (publisher/repository text); `metadata`. |
| 16 | Audibert and Haddar, limited-aperture GLSM, `10.1137/16M110112X` | Extends the generalized linear sampling method to limited-aperture far-field data where the operator factorization is non-symmetric, with a modified regularization functional. | Limited/partial aperture changes the measurement operator's structure; factorization/asymmetry of the measurement operator is analyzed. | Qualitative (sampling/indicator) method for target support; does not recover contrast or currents and does not estimate unknown aperture geometry. | `abstract` (HAL text); `metadata`. |

## Cross-cutting observations

1. The generic antecedents are confirmed as prior art: equivalent-source /
   contrast-source manipulation, self-calibration, joint image + calibration
   optimization, and source/receiver extension. None of these can be claimed
   as novel.
2. The only anchors that operate in the same current space as SOM/TSOM
   (`G_S`/`G_D` subspaces) assume fixed, known geometry (#1, #2, #13, #14,
   #15); the anchors that handle unknown geometry or calibration live either
   in data space, in a phase/shift-kernel image model, or on the
   measurement chain (#3-#12).
3. No retrieved anchor defines a pose-dependent SOM/TSOM spectral partition
   `V_S^±(X)`, `V_D^±(X)` or a `V_P` subspace **inside the retained current
   space** and then uses the joint (data residual + state inconsistency) of a
   canonical lift for self-calibration. That absence is a retrieval gap, not
   a novelty conclusion.
