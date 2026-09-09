# Trispace prior-art search plan

Search date: 2026-09-04

Depth: deep review / broad discovery, not a systematic review. Semantic
Scholar is the primary academic source, transported via `scholarqa-cli`
(version recorded in the transport record). Final novelty judgment remains
with Codex and is bounded to retrieved evidence.

## Main question

Does prior literature already define and analyze the extension of subspace
optimization methods (SOM/CSI families) for full-wave inverse scattering to
*unknown or unstable sensing Green operators / unknown antenna positions*, by
representing geometry-induced operator variation as a structured
equivalent-current perturbation, jointly constrained by contrast-source
physics and calibration geometry?

## Answerable subquestions

1. Is there prior SOM / subspace-optimization inverse-scattering work that
   explicitly handles unknown transmitter/receiver positions or an unknown
   Green operator / measurement kernel?
2. Does contrast-source-inversion literature treat position/geometry errors of
   the measurement setup (not only contrast noise or model error) and, if so,
   how is the unknown parameterized and constrained?
3. Is there a line of inverse-scattering self-calibration or antenna-array
   calibration literature that estimates geometry and scatterer/contrast
   jointly?
4. Do radar/RF autofocus or joint image/platform-motion methods share the
   structured-perturbation or bilinear calibration mechanism with the proposed
   equivalent-current reformulation?
5. Is there a "structured equivalent-current perturbation" representation
   (operator perturbation absorbed or parameterized through the induced
   current / contrast-source degree of freedom) in scattering or imaging
   literature?
6. Which methods address identifiability, ambiguity/gauge, and convergence of
   joint sensor-geometry and scene/contrast recovery?

## Query families (one or more short formulations each)

F1 SOM inverse scattering, unknown sensors:

- `subspace optimization inverse scattering unknown sensor positions`
- `subspace based optimization inverse scattering antenna position error`

F2 Contrast-source inversion, geometry error:

- `contrast source inversion antenna position errors`
- `contrast source inversion measurement geometry error`

F3 Inverse-scattering self-calibration / array geometry:

- `inverse scattering self calibration array geometry`
- `electromagnetic inverse scattering unknown antenna positions joint estimation`

F4 Unknown Green function / kernel:

- `inverse scattering unknown Green function calibration`
- `electromagnetic imaging unknown measurement kernel estimation`

F5 Equivalent-current calibration:

- `equivalent current source calibration antenna`
- `induced current perturbation sensor calibration inverse scattering`

F6 Radar autofocus / joint image-platform:

- `radar autofocus joint image platform motion estimation`
- `synthetic aperture radar phase error autofocus unknown positions`

F7 Blind calibration / bilinear inverse problems:

- `blind calibration bilinear inverse problems sensor gains`
- `blind deconvolution calibration nuisance parameters inverse problems`

## Screening priorities

- Retain the closest SOM/CSI, position-error sensitivity/correction,
  self-calibration, autofocus, blind-calibration, equivalent-current, and
  unknown-kernel work;
- actively retain contradictory or overlapping work;
- do not infer substantive results from titles, citations, or graph edges
  alone;
- separate generic algebraic antecedents from papers matching the full
  physical combination (equivalent-current-constrained operator
  perturbation + contrast-source physics + calibration geometry);
- verify every paper used in the claim ledger with `scholarqa-cli verify`.

## Planned outputs

- `transport_record.md` or `evidence.json` — bundle from `scholarqa-cli collect`;
- `screened_papers.md`;
- `claim_ledger.md`;
- `verified_references.json`;
- `novelty_overlap.md`;
- `summary.md`.
