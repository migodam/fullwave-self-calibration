# Closest-overlap analysis (retrieval-bounded)

Date: 2026-09-04. Object analyzed, in the exact restricted form:

1. unknown/unstable antenna geometry perturbs the sensing Green operator
   `G_S(X)`;
2. a canonical pose-equivalent induced-current lift is constructed **inside a
   retained SOM/TSOM current space** (not by bare equivalence in the full
   current space);
3. the construction retains an **irreducible data residual** and the
   lifted current's **state-equation inconsistency**;
4. those two objects are used jointly for self-calibration.

The bare equivalence statement (any geometry-induced data change can be
written as some current change) is treated as prior art and, worse, as
vacuous in the unrestricted case: when `G_S` is an `m x S`, `S >= m`,
full-row-rank matrix, every single-snapshot data vector lies in
`Ran G_S`, so every pose-induced perturbation is exactly representable. The
meaningful object must be the lift restricted to the retained
(TSOM-selected) current subspace, or to a shared/stacked pose-degree
structure that breaks the single-snapshot degeneracy. Overlap below is
assessed against that object.

## Closest conceptual overlaps

### 1. Receiver/source extension in FWI (#5 Huang et al.; #6 Metivier and
Brossier; #7 da Silva et al.)

Closest genus. Source-receiver extension and receiver relocalization already
(a) add artificial source/receiver degrees of freedom, (b) absorb data misfit
that would otherwise look like model error, and (c) estimate a corrected
geometry-dependent mapping as part of inversion. The receiver-extension work
in particular treats receiver-coordinate inaccuracy as the nuisance and uses
the misfit reduction as the driving signal.

Decisive difference from the claimed object:

- extension lives in source/receiver coordinates of a seismic wave-equation
  forward model, not in the induced-current (contrast-source) space of a
  full-wave EM/Lippmann-Schwinger model;
- there is no retained SOM/TSOM spectral subspace whose boundaries change
  with pose, and no minimum-norm lift `G_S(X)^+ dG_S[X] j` constrained to
  that subspace;
- the residual used is the overall data misfit (after extension), not a
  separately retained irreducible data residual paired with a
  state-equation inconsistency of the lifted current;
- "receiver position as an artificial DOF in FWI" is therefore prior art and
  must be cited, while the current-space residual-pair construction is the
  only part that remains retrieval-unresolved.

### 2. SOM/TSOM current-space substrate (#1 Chen; #2 Zhong-Chen; #14 CC-CSI;
#15 Bevacqua et al.)

Closest algebra. These papers establish that current/contrast-source
variables live in a current space, that `G_S` and `G_D` spectral structure
defines retained vs complementary current parts, and that data and state
residuals can both be incorporated (CC-CSI cross-correlates the two error
terms; virtual-experiment methods recombine measurements to shape contrast
sources).

Decisive difference:

- geometry is fixed and known in all of these papers; `V_S^±`, `V_D^±` do
  not depend on pose;
- no paper in this set defines a pose-induced perturbation of the
  current-to-data operator and lifts it into the *retained* current subspace
  with the resulting leakage/irreducible residual computed;
- CC-CSI's two error terms are an inversion objective under known geometry,
  not a calibration residual pair for an unknown operator;
- Bevacqua et al. use "virtual experiments" as a solver device to make the
  contrast sources focused, not to estimate an unknown sensing operator.

### 3. Joint scene estimation with unknown geometry/phase (autofocus and DNN
lines: #8 Karthik-Ghosh; #10 Mansour; #11 Onhon-Cetin; #12 Scarnati-Gelb)

Closest problem statement. Joint recovery of an image/contrast with unknown
transmitter positions, antenna-position shift kernels, or phase-error fields
is prior art, including in inverse scattering itself (Karthik-Ghosh).

Decisive difference:

- radar autofocus uses a linearized image/phase model or sparse
  blind-deconvolution representation, not a nonlinear full-wave volumetric
  contrast-current model;
- Karthik-Ghosh uses a data-driven DNN and a transmitter-localization search;
  the retrieved text does not establish a spectral current-space reduction or
  a pose-defect residual theory;
- none of them computes a SOM/TSOM-restricted pose-lift or a state-defect
  inconsistency, so the joint-inversion framing is antecedent but the
  restricted operator-level construction is not found in this set.

### 4. Calibration of the measurement chain or data (#3 Hanabusa; #4
Cathers)

Closest practical calibration antecedent inside microwave imaging: DNN-based
data conversion before CSI and 2-port error models for the whole imaging
system.

Decisive difference:

- these are calibration stages with known array geometry (systematic
  measurement errors, antenna/field propagation models), not pose-dependent
  `G_S` estimation interleaved with a retained-current reconstruction;
- they calibrate the data/measurement hardware, not a canonical lift of the
  geometry error into current space.

### 5. Identifiability/gauge theory (#9 Li-Lee-Bresler) and aperture
structure (#16 Audibert-Haddar)

Methodological context. Bilinear identifiability up to transformation groups
governs any contrast-times-operator bilinear structure, and limited-aperture
analysis shows that aperture-induced asymmetry of the measurement operator is
already studied in qualitative imaging.

Decisive difference:

- Li-Lee-Bresler is generic and not specialized to Helmholtz operators,
  contrast sources, or the pose-current residual pair;
- Audibert-Haddar is a fixed-geometry qualitative method; it supplies
  factorization/aperture context only.

## Bounded overlap verdict

Within this verified set:

- every *generic* ingredient (equivalent/contrast-source representation,
  self-calibration, joint image + geometry/phase estimation, source/receiver
  extension, data/state residual combination in CSI, virtual measurements)
  is prior art;
- the exact restricted object — a pose-equivalent current lift inside a
  retained SOM/TSOM current subspace whose irreducible data residual and
  state-equation inconsistency are jointly used for self-calibration — was
  not found in the 16 anchors or in the local digests/addendum;
- a **full row-rank vacuity warning** applies: without the retained-subspace
  restriction (or stacked/shared pose structure), the "pose subspace"
  construction is degenerate and any novelty claim built on bare equivalence
  would be overclaiming;
- this verdict is retrieval-bounded: it is not evidence of global absence and
  does not decide whether the restricted object is a supported recombination,
  a publishable narrow formulation, or a variant already present in literature
  outside this set.

## Residual risks for Codex judgment

1. Whether "irreducible data residual" is meant as
   `P_{ker G_S*} dG_S[δX]j` (sensing path leaves `Ran G_S`) or as
   retained-subspace leakage `P_{Ran U_T^⊥} δj_eq`; the two are different and
   have different vacuity properties.
2. Whether the state-equation inconsistency is measured on the lift only or
   on the full pose-perturbed induced current (which includes
   `M^{-1} D_χ dE^inc[δX]`, already inside `Ran G_S` under world-fixed grids).
3. Whether per-channel phase/gain calibration belongs to "pose" at all; the
   retrieved anchors #3/#4 show measurement-chain calibration is a separate,
   prior-art track.
4. Whether "inside a retained SOM/TSOM current space" means the retained
   `V_S^+`/`V_D^+` subspace of the *nominal* or the *perturbed* operator;
   the pose-dependence of the subspace itself (`V_S^±(X)`) is open and not
   resolved by any anchor.
