# Exact full-wave object/receiver/clock translation gauge

Parent derivation, 2026-09-08. Translation covariance of scattering is an
established physical symmetry; the following reconstruction is an assumption
audit, not a claim of first discovery. It identifies what our known-support
inverse experiments are anchoring and why broader imaging needs a gauge check.

## Physical setup

The background is spatially homogeneous and the scatterer is described by a
compact material distribution chi. For illumination l, the world-fixed
incident field is p_l exp(i k d_l·r), with known unit direction d_l and fixed
transverse polarization p_l. The receiving component/direction is unchanged
by translating a receiver. Frequencies share electronic coefficients a_l,
while ell is an unknown common delay length. The scattered data are

$$
\mu_{krl}=a_l e^{ik\ell}E^s_{kl}[\chi](x_r).
$$

Assume the physical forward solution exists uniquely. The admissible material
family must include the translated distribution under consideration; fixed
known supports can deliberately break this symmetry. No Born approximation
or weak-scattering assumption is used.

## Exact symmetry and proof

Define chi_t(r)=chi(r-t). Translation invariance of the background Green
kernel, together with the incident plane-wave phase, gives

$$
E^s_{kl}[\chi_t](x_r+t)=e^{ikd_l\cdot t}E^s_{kl}[\chi](x_r).
$$

One direct proof is to translate the volume-integral equation: the transformed
current is j_t(r)=exp(i k d_l·t)j(r-t). Changing integration variable preserves
the interaction kernel and reproduces the incident field at r. Uniqueness then
identifies this current with the full multiple-scattering solution. Substitution
in the radiation integral yields the displayed field identity.

If a scalar c satisfies d_l·t=c for every illumination, then

$$
(\chi,x_r,\ell,a_l)\mapsto
(\chi_t,x_r+t,\ell-c,a_l)
$$

leaves every measurement unchanged at every frequency. Thus continuous
translation gauges include the family satisfying

$$
D t=0,\qquad
D=\begin{bmatrix}(d_2-d_1)^T\\\vdots\\(d_L-d_1)^T\end{bmatrix},
\qquad c=d_1^Tt.
$$

This constructed family is not an exhaustive classification of the full inverse
problem's ambiguities. Zero scattering channels, further electronic/material
symmetries or phase-branch ambiguities may exist. Unknown frequency-independent
illumination phases can additionally create discrete aliases at commensurate
frequencies; the continuous family above already establishes nonidentifiability.

For the two distinct directions used in our current 3D experiments, d1=e_z
and d2=e_x (two polarizations each), D has rank1. Hence t_x=t_z with arbitrary
t_y is a two-dimensional family. Extra polarizations with the same directions
do not increase this direction-difference rank. Four affinely independent
directions, for example tetrahedral directions, remove this particular common-
clock translation family, but do NOT prove complete material/array identifiability.

If each illumination has an independent unknown delay length ell_l, every
translation is compensated by ell_l -> ell_l-d_l·t. More frequencies do not
remove that exact gauge. If the common clock is known instead, d_l·t=0 for
all l is sufficient for an unchanged-gain continuous translation family.

## Meaning for our results and algorithm design

1. Receiver-translation recovery with known object support uses a physical
   spatial anchor. It is not evidence of absolute object and array localization
   when both supports and positions can translate freely.
2. A fixed material basis can break an underlying physical gauge by excluding
   translated materials. Apparent finite-dimensional rank must therefore be
   interpreted together with that prior, not as an intrinsic full-field guarantee.
3. A consistent larger inverse problem must fix a justified gauge or estimate
   relative geometry/identified combinations. Gauge fixing chooses coordinates;
   it does not add physical information or prove absolute location recovery.
4. A genuinely known spatial reference, an appropriate phase/clock reference,
   or additional calibrated directions may constrain this family. Which remedy
   is physically available must be declared, not smuggled in through truth.
5. Removing this symmetry is necessary for absolute attribution but not sufficient
   for nonlinear branch coverage or fine material recovery. Those remain separate.

Numerical Maxwell covariance checks are prepared in
`check_translation_clock_gauge.py`; the executed status is determined by its
result artifact, not by this prose. This proof does not certify global novelty.
