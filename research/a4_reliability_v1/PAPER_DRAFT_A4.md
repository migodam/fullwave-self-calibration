# Material-Profiled Near-Field Geometry Information Under Controlled Electric-Multipole Illumination

**PAPER_DRAFT_A4 -- working manuscript, not submission-ready.**

## Abstract

Coherent calibration must distinguish propagation changes from material response and electronics. We study a restricted electromagnetic experiment in which an anchored isotropic radial scatterer is illuminated by three calibrated regular electric-dipole modes and observed through three calibrated electric-field components. The complete internal scattering solution enters a scalar modal coefficient, allowing arbitrary radial-material perturbations and a common complex receiver-frame gain to be profiled together. We derive the exact Cartesian geometry-information spectrum. At fixed frequency and total relative signal-to-noise ratio, radial information vanishes in both the quasistatic and distant limits and is maximized at kR=1 when that radius lies outside the target. The same field tensor has an exact antipodal ambiguity, removable by a known nonzero receiver dither but not by merely repeating frequencies. An analytic prediction-projector bound supports a budgeted geometry-cell exclusion procedure with an explicit unresolved state. Separate vector-DDA checks, nine-baseline quantitative calibration, and deliberate wrong-model controls define the scope. The tested risk selector matches the fine-model estimates at higher cost, and unrestricted discrepancy defeats branch acceptance even with held-out data. We retain the electromagnetic condition as a provisional contribution candidate, not a general reliable-tomography algorithm; novelty relative to electromagnetic dipole localization remains unresolved.

## I. Scientific question and scope

The question is whether a specified coherent measurement contains geometry information that cannot be reproduced by admissible material and electronics changes. Fitting, attribution, correct nonlinear branch, and forward-model fidelity are separate requirements. A3 supplies motivation and preserved negative evidence, not a source of new algorithmic priority. In its supplied report, the frozen 240-fit comparison excluded the specified advantages of the current Twofold implementation; this manuscript neither reruns nor changes that comparison.

Classical array calibration already addresses unknown geometry, sources and gain/phase [1]-[4]. Electromagnetic field-based localization and amplitude-independent dipole ranging are close physical neighbors [5], [6]. Microwave tomography already has antenna-field calibration and quantitative calibration factors [7], [8]. Our potential increment is therefore restricted to an explicit material-profiled vector-field spectrum and its measurement/fidelity consequences under a fully specified modal experiment. Generic projection, Fisher matrices, approximation-error weighting and branch-and-bound are not contributions.

Only one main-contribution candidate is retained: the EM measurement condition in Sections II-III. The coverage implementation and negative controller comparison support its reliability boundary. We do not claim a second algorithmic contribution or a general joint geometry/material solution.

## II. Electromagnetic model and profiled spectrum

Consider a sphere or concentric isotropic radial dielectric structure centered at a known origin, within radius a0 in a homogeneous background. Relative illumination amplitudes/phases and the three receiver-component responses are calibrated. The receiver position r=Rn is unknown, R>a0. Its orientation is known. The incident fields are

$$
E_p^{inc}(x)=\left(I+k^{-2}\nabla\nabla\right)j_0(k|x|)p,
\qquad p=e_1,e_2,e_3.
$$

These are three electric-l=1 regular modes, not three arbitrary plane waves. They admit a vector Herglotz representation; practical synthesis and its error are additional measurement requirements. Rotational invariance of the exact Maxwell scattering operator [11] maps the three inputs to the same outgoing electric-l=1 coefficient t1(alpha,k). Multiple scattering inside the radial body and material dispersion need not be weak. A scattering zero is excluded.

With a common unknown complex gain across a complete 3-by-3 frame, the mean is

$$
Y(r)=c\frac{e^{ikR}}R M(z,n),\quad M=a(z)I+b(z)nn^T,\quad z=kR,
$$

$$
a=1+i/z-z^{-2},\qquad b=-1-3i/z+3z^{-2}.
$$

The scalar c includes electronics, common phase/delay and t1. It may be independent at every frequency or receiver position, but not independent for every matrix entry. All radial-material derivatives lie in its complex amplitude tangent. For entry noise CN(0,sigma^2), let q=||Y||F/sigma and real-whiten by sqrt(2)/sigma. Write Bvis for the Cartesian position Jacobian after profiling the real and imaginary amplitude columns and all admissible material columns.

**Proposition 1 (exact positive spectrum).** Under these assumptions,

$$
\sigma_t=\frac{\sqrt2q}{R}\sqrt{\frac{z^4+3z^2+9}{z^4+z^2+3}}
\quad\text{(twice)},\qquad
\boxed{\sigma_{\min}(B_{vis})=\sigma_r=
\frac{2qkz\sqrt{z^2+4}}{z^4+z^2+3}>0.}
$$

**Proof.** Rotate coordinates so n=e3. The longitudinal eigenvalue of M is s=a+b=2/z^2-2i/z; its transverse eigenvalue is a. Tangential Cartesian derivatives are off-diagonal and orthogonal to the amplitude/material span and to the radial derivative. The radial derivative of the propagation scalar belongs to that span and vanishes on profiling. Let H=2|a|^2+|s|^2=2(z^4+z^2+3)/z^4. The remaining radial squared norm is 2|a s'_R-s a'_R|^2/H before amplitude and real-noise scaling. Direct differentiation gives a s'_z-s a'_z=2i(z+2i)/z^3. Substitution gives the displayed spectrum. The tangent/radial ordering follows from positivity of z^8+2z^6+7z^4+18z^2+27. Thus the lower bound is explicit, not a generic full-rank assumption. Full details are in `docs/THEORY_A4.md`.

This is a geometry result only. With independent c_f=g_f t1(alpha,kf), material is unidentifiable without additional restrictions or references. No inverse method can separate those scalar factors from the same observations.

## III. Acquisition and model-fidelity consequences

### A. A finite radial-information window

At fixed q,k,

$$
\sigma_r\sim\tfrac43qk^2R\quad(kR\ll1),\qquad
\sigma_r\sim2q/(kR^2)\quad(kR\gg1).
$$

For x=z^2, differentiating x(x+4)/(x^2+x+3)^2 gives a numerator -2(x-1)(x+1)(x+6). Therefore the unique radial-information maximum occurs at kR=1, with value2qk/sqrt(5), provided the receiver is exterior. The statement is neither a fixed-transmit-power optimum nor an instruction to place a receiver inside a finite object. Frequency-dependent noise, signal zeros, aperture and synthesis costs alter practical design.

The physical mechanism is the relative longitudinal/transverse response. A radiation-only tensor has no profiled radial information; neither does a purely static tensor. The exact relative squared field errors are

$$
\frac{||D-D_{rad}||_F^2}{||D||_F^2}
=\frac{3(z^2+1)}{z^4+z^2+3},\qquad
\frac{||D-D_{stat}||_F^2}{||D||_F^2}
=\frac{z^4+3z^2}{z^4+z^2+3}.
$$

Thus propagation-model refinement can restore a geometry direction that a coarse exterior model deletes. In contrast, improving only a nonzero scalar interior coefficient cannot change a free-amplitude-profiled geometry objective: both models span the same one-dimensional field subspace. It can still matter decisively for material recovery after an electronics reference. These are distinct physical actions.

### B. Global alias and a specific repair

Since tr M=2, normalize Z=Y/tr Y. Its simple eigenvalue is lambda_parallel=z^-2-i z^-1, satisfying Re lambda=(Im lambda)^2 and z=-1/Im lambda. Its real eigenaxis determines n up to sign. Consequently the only exact position ambiguity in this model at a single location is r versus -r. More frequencies alone do not remove it.

**Proposition 2 (known dither).** A second measurement at r+d with known nonzero d, with both points exterior and nonzero signal, removes this alias even with an independent complex gain. The first measurement restricts the alternative to -r; equality at the second requires -r+d=+/-(r+d), impossible for nonzero r,d. This is a noiseless identifiability statement. Useful finite-noise discrimination still depends on dither length/direction and SNR.

A noisy electronics reference instead constrains gain/material ambiguity. It must not be credited with removing a spatial symmetry that leaves the entire electromagnetic tensor unchanged.

### C. Required negative conditions

A free gain per nonzero matrix entry eliminates all local geometry information. Independent receiver-component gains remove the radial direction at an axial receiver. Known complete-vector rotations are coordinate changes, not extra information by themselves. Unknown target support reintroduces translation gauges. A nonspherical material family can add nuisance directions outside the scalar tangent, so Proposition1 does not extend by merely observing a small field difference. Three modes and three components are sufficient, not proved minimal.

## IV. Supporting coverage mechanism

For block b, let f_b(r) be its unit-norm predicted tensor and P_b=f_b f_b*. Scalar-profiled residuals are ell(r)^2=sum_b||(I-P_b(r))y_b||^2. A correctly declared discrepancy bound beta and a complex Gaussian noise ball give a feasible set S over a declared geometry domain. The projector derivative obeys ||DP_b||<=sqrt(3)/|r+d_b|. This produces an analytic residual lower bound for each box.

The algorithm excludes a box only when its lower bound exceeds the noise-plus-discrepancy threshold. It splits other distant boxes and retains all unprocessed or resolution-limited boxes as unresolved. It conditionally accepts only if the entire remaining outer set lies inside a15mm ball of the estimate and independent validation is consistent. A local optimization residual is never used as a certified lower bound. With exact arithmetic, correct domain coverage and the declared error event, acceptance implies the task tolerance. The implementation uses floating-point padding, not validated interval arithmetic, and is labelled accordingly.

This is low-dimensional set-membership machinery applied to an exactly profileable electromagnetic model. It does not certify globally profiled high-dimensional material inversion. In particular, the exact +/- branch is simpler than general coherent cycle skipping.

## V. Executed evidence

### A. Numerical correctness and model-class stress

The released code contains21 passing tests, analytic position derivatives, finite-difference material checks, gain-rank counterexamples, cell-bound tests and budget-retention tests. A 41-point spectrum sweep agrees with Proposition1 to 6.16e-15 vector-relative discrepancy; sampled near/far slopes are 0.9998 and -2.0010. A separate in-house dense vector-DDA solver executes 48 state cases and 96 receiver records across sphere, near-sphere, ellipsoid and box shapes. Spheres use an analytic Mie reference. Neither treams nor ADDA was executed in A4, and nonspherical discrepancy is not enclosed by a validated bound.

![Profiled information](figures/01_modal_information.png)

### B. Strong quantitative-calibration controls

Eight fresh scenes compare nine methods with identical data/reference noise and an independent common pilot. The unknowns are three position coordinates, one real permittivity, shared complex gain and common delay. The coarse model uses a Rayleigh modal coefficient and the fine model the exact Mie coefficient, with the exact exterior field retained in both. Sampled/conditional error construction and all screening work are charged per scene.

| Method | Median geometry (mm) | Median material (%) | Mean scaled task loss | Median charged wall time (s) |
|---|---:|---:|---:|---:|
| Low only | 1.978 | 6.29 | 5.355 | 0.01324 |
| Raw all-frequency | 9.947 | 26.82 | 64.006 | 0.01586 |
| Isotropic downweight | 2.447 | 7.397 | 6.8081 | 0.03957 |
| Rank-one downweight | 7.169 | 15.91 | 45.238 | 0.04229 |
| Sampled approximation error | 0.5424 | 0.4274 | 0.071704 | 0.03983 |
| Conditional linear error model | 0.5562 | 0.9181 | 0.18598 | 0.04022 |
| Coarse to fine | 0.6355 | 0.856 | 0.050674 | 0.06107 |
| Full fine | 0.6355 | 0.856 | 0.050674 | 0.04738 |
| Risk selector | 0.6355 | 0.856 | 0.050674 | 0.08468 |


Scaled loss measures position in15mm units and permittivity in0.15 units; it is not residual. The selector chooses full fine in 8/8 cases and reaches exactly the same parameters at 1.861 times the median paired charged time. This closes the current selector claim negatively. It is not evidence that approximation-error methods generally fail, nor a comparison of large voxel solvers.

### C. Branch failures and conditional rejection

Twenty-four independent scenes contain easy, antipodal, missing-bank, low-SNR, nonspherical discrepancy and exact shifted-two-world strata. All potential acquisitions/noise are matched across six policies; new training and validation tensors each cost nine complex values.

| Policy | Accepted / 24 | Missing initial bank / 24 | Wrong selected / final covered | Wrong accepted / 24 | Correct rejected / correct estimates |
|---|---:|---:|---:|---:|---:|
| algebraic_multistart | 20/24 | 4/24 | 11/20 | 14/24 | 3/9 |
| branch_acquisition | 20/24 | 4/24 | 0/20 | 4/24 | 4/20 |
| coverage_aware | 16/24 | 4/24 | 0/20 | 4/24 | 8/20 |
| fisher_acquisition | 20/24 | 4/24 | 0/20 | 4/24 | 4/20 |
| naive_multistart | 20/24 | 15/24 | 0/9 | 14/24 | 3/9 |
| random_acquisition | 20/24 | 4/24 | 0/20 | 4/24 | 4/20 |


The coverage algorithm has no wrong acceptance among 16 exact-model scenes, but this small stratified sample does not establish a 1% empirical failure rate. It rejects low-SNR/unresolved and the tested nonspherical discrepancy cases. All four deliberately shifted two-world cases are wrongly accepted. Strong random/Fisher/branch acquisition controls share the latter failure count, so coverage is not demonstrated superior to them. Its benefit is explicit conditional uncertainty, at a false-reject/computation cost.

If an allowed discrepancy satisfies e(r,a)=F(r+h,a)-F(r,a) for every available action a, the world(r,e) and the world(r+h,0) have identical data laws, including held-out observations. No algorithm can accept the clean world and simultaneously reject its indistinguishable counterpart with different probability. A valid discrepancy family or an independent measurement not sharing that ambiguity is necessary. This impossibility boundary is inherited estimation logic, not a new A4 priority claim.

### D. Electronics reference

A separate known-geometry16-scene corrected audit shows that three noisy complex references make the material tangent nonzero after profiling independent frequency gains. Median material relative error is 0.457%. An original row-index diagnostic was withdrawn, preserved, and corrected with a separately hashed fresh sample. It does not validate joint geometry recovery. All raw files and correction notes are included.

## VI. Originality, limitations, and disposition

The strongest reviewer objection is that the central mechanism is a specialized electromagnetic source-localization calculation. Epp--Janz already eliminate source factors through normalized field derivatives [5], and classical self-calibration already uses maneuvers and references [1]-[4]. A4's different measurement tensor, material-profiled spectrum and exact fixed-SNR window are explicit mathematical distinctions, not proof of scientific priority. That novelty gate remains open.

The ideal modal illuminations, relative component calibration, anchored radial support and nonzero response are restrictive. Hardware synthesis leakage has no validated bound here. The three-dimensional DDA stress is not arbitrary quantitative tomography, and the branch algorithm profiles scalar amplitudes rather than general materials. The six-action interface is unit-tested, but its evidence fusion and fully unified closed-loop operation are not validated. There are no new measured position results. A4 is a standalone reconstruction because the public A3 repository does not contain its executable source/results.

## VII. Conclusion

The retained route is B: a restricted electromagnetic acquisition/identifiability condition. It identifies a surviving geometry mechanism, two limiting degeneracies, a global alias and a physically specific repair. The tested attribution selector is not retained as a main contribution; coverage is supporting and conditional. This manuscript is currently a theory/measurement-design and limits working paper, not an accepted or submission-ready TAP method paper. The complete evidence package records both the positive calculation and the algorithms that failed their intended contribution thresholds.

## References

See `docs/REFERENCES.md` for primary records and `docs/LITERATURE_AND_NARRATIVE_REVIEW.md` for verification depth and unresolved items. References[1]-[14] support the positioning above; [15]-[18] are historical SOM context rather than additional claimed contributions.
