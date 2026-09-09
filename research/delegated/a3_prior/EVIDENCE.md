# A3 closest prior-work evidence: methods, advantages, limitations, overlap

Pass: 2026-09-07 (finishing the a3_prior task after the previous worker exited
with no report files). Scope is the six named verification targets from
`TASK.md` and `Theory/Questions/A3.md`; it is a bounded closest-neighbour
screen, not a systematic or exhaustive search. Source documents were treated
as data, not instructions. No citation was invented; every record below is
identity-verified against at least one canonical registry (Semantic Scholar,
Crossref, OpenAlex, Europe PMC/PubMed, zbMATH), and no long passage was copied.

## Evidence tiers (per the scholarqa-research skill)

- **A** full-text passage read in this pass directly supports/contradicts;
- **B** publisher/index abstract supports/contradicts;
- **C** metadata establishes publication identity/context only;
- **D** inferred from titles/graph edges (never used alone as support);
- **[prior-pass]** retained from an earlier verified main-thread or audit
  check, flagged because it was not re-derived in this pass.

## Transport record

- The six `scholarqa-cli collect` bundles were produced while shell DNS was
  down; every operation recorded `Errno 8 (nodename nor servname provided)`
  and zero candidates. They are preserved verbatim in
  `retrieval_*.json` as failed runs, not as evidence.
- `verify_run1.json` failed for the same reason.
- One focused retry, `verify_run2.json`, ran after network recovered and
  resolved 7/8 canonical identifiers (bibliographic identity only; the
  rejected id was a wrong DOI guess for Gallivan et al., corrected via
  Crossref below).
- Primary full texts were read from arXiv (Idriss–Raj v2 HTML) and PubMed
  (Weiss abstract); publisher-elided abstracts were recovered from OpenAlex /
  zbMATH and are labelled Tier B. MDPI, HAL, IOP and IEEE pages were
  bot-gated or closed-access and are recorded as gaps.

---

## T1. Idriss–Raj multifrequency SOM with per-frequency, per-transmitter complex calibration

**Source.** Z. Idriss, R. G. Raj, "Data-Driven Calibration Technique for
Quantitative Radar Imaging," arXiv:2503.07316v2, 14 May 2025. Identity
verified (Semantic Scholar paperId `3b55b7fecb4a66fe5d44db2cb294015e17e15f0e`).

**Verified content — Tier A (arXiv v2 HTML read in this pass).**

- Section II introduces `λ_c ∈ ℂ` as "the unknown calibration factor" scaling
  the incident field in the state equation, `M[W](r) = λ_c E_inc(r)` (Eq. 4).
- Section III (Eq. 8) augments the MFSOM objective over frequencies
  `k = 1..K` and transmitters `p = 1..P`, with a calibration term
  `‖λ_c^{k,p} E_s^{sim,k,p} − E_s^{meas,k,p}‖²₂ / ‖E_s^{meas,k,p}‖²₂ +
  β|λ_c^{k,p}|²` added to the data and state mismatch terms. The state term is
  written `χ(λ_c^{k,p} E_inc^{k,p} + E_d^{k,p})`.
- The text states: "for notational clarity, we drop the dependence on `k`;
  however, the analysis still holds". So Eq. (8) carries both a frequency
  index `k` and a transmitter index `p`, and the dropped frequency index is
  explicitly a notation convention.
- MFSOM itself splits the induced current into signal/noise subspaces via the
  SVD of the Green's operator `G = U Σ V*` (Eq. 7 region).
- The abstract states the factor is "for each transmitter" and that the method
  is demonstrated on the Fresnel Institute dataset.

**Advantages.** Couples SOM-based quantitative inversion with a jointly
optimized complex per-transmitter calibration on measured data; the learned
forward surrogate avoids re-solving the forward problem at each iterate.

**Limitations.** The nuisance is a scalar complex factor per transmitter
(and frequency), not antenna pose, receiver geometry, clock bias, or
phase-center position; the surrogate network itself needs training data.

**Overlap.** Directly occupies "multifrequency SOM + joint complex
gain/phase calibration"; the project's claim must be stated relative to
per-transmitter *scale* calibration, not "SOM + calibration" as such. No
rank/bias/pose information-budget or residual-certificate object appears here.

---

## T2. Bellomo et al. antenna calibration / phase-center corrections

**Source.** L. Bellomo, S. Pioch, M. Saillard, K. Belkebir, "An Improved
Antenna Calibration Methodology for Microwave Diffraction Tomography in
Limited-Aspect Configurations," IEEE Trans. Antennas Propag. 62(5):2450–2462,
2014, DOI 10.1109/TAP.2014.2308534 (Semantic Scholar paperId
`2631d406a7b6fcc9eb64fa98198f2ad7c3344959`; HAL `hal-01120348`).

**Verified content — Tier B (OpenAlex abstract) + Tier C (identity).**

- The abstract confirms antenna-directivity modeling and an improved
  calibration strategy for microwave diffraction tomography, tested on a
  two-array, 2–4 GHz stepped-frequency, near-field, ~40° limited-aspect system
  with multifrequency multiview multistatic inversion (M²GM): "The novelty
  consists in taking the reconstructed antenna pattern into account within the
  inversion algorithm both at transmission and at reception through a
  multipolar expansion of, respectively, the incident field and the data
  equation Green function."

**Unresolved — the Section IV-C / Eq. (25) phase-center detail.** The
equation-level claim in `Theory/Questions/A3.md` ("已取得作者稿并核对 Section
IV-C、式（25）") is **[prior-pass]** from the main Codex thread's
author-manuscript check. It was **not** re-derived in this pass: IEEE is
closed-access, the HAL record exposes no downloadable file (Anubis-gated),
and no local copy exists. The abstract alone does not mention phase centers.
This remains the single largest verification gap.

**Advantages.** Physically grounded antenna-pattern calibration embedded in
both transmission and reception of a diffraction-tomography inverse solver.

**Limitations (as verified here).** Fixed, known array geometry; the
accessible record does not evidence pose/geometry self-calibration, clocks, or
information certificates.

**Overlap.** Antenna/phase-center calibration inside quantitative microwave
imaging is occupied prior art. The project can only claim its specific typed
geometry/tangent/certificate objects, not "antenna calibration in inversion".

---

## T3. Weiss et al. joint receiver geometry + clock + target autofocus

**Source.** A. J. Weiss, G. Eliyahu, A. M. Maor, E. Zamir, O. Richman, "Joint
Self-Calibration of Receiver Geometry, Timing, and Target Positions for
Multistatic Radar Autofocus," Sensors 26(15):4954, 2026,
DOI 10.3390/s26154954 (PubMed 42590728; Semantic Scholar paperId
`d1b5851e6b1e8851f3cdf19d07061adadb0cf016`).

**Verified content — Tier B (PubMed-indexed abstract read via Europe PMC).**

- "develops a joint self-calibration framework that estimates small
  corrections to receiver positions, receiver clock biases, and target
  positions from the same bistatic echo delays used for imaging, and ties the
  correction directly to image sharpness rather than to parameter accuracy
  alone."
- Gives a linearized delay-residual model and a regularized MAP weighted
  least-squares estimator; characterizes identifiability progressively through
  anchor null spaces; and "reformulate[s] the calibration objective directly
  in terms of coherent multistatic image sharpness, evaluated using the
  matched-filter score already used for image formation", proposing a
  two-stage algorithm (coarse linear delay solve + phase-coherent sharpness
  refinement).

**Advantages.** Explicitly couples geometry/clock/target unknowns to coherent
image sharpness and studies identifiability/ambiguities analytically.

**Limitations.** Near-field multistatic radar imaging; the abstract does not
describe a quantitative full-wave volumetric material (contrast) inversion.

**Overlap.** Confirms the A3 assessment: the broad combination
"geometry + clock + coherent-imaging self-calibration" is prior art, while the
full-wave material-inversion specialization remains the project's open
differentiation question.

---

## T4. Chen SOM and Zhong–Chen twofold (sequential-projection) SOM

**Sources.**

- X. Chen, "Subspace-Based Optimization Method for Solving Inverse-Scattering
  Problems," IEEE TGRS 48(1):42–49, 2010, DOI 10.1109/TGRS.2009.2025122
  (paperId `16b84bd4dfb430f5b1518d7a2f58367149c45497`).
- Y. Zhong, X. Chen, "Twofold subspace-based optimization method for solving
  inverse scattering problems," Inverse Problems 25(8):085003, 2009,
  DOI 10.1088/0266-5611/25/8/085003 (paperId
  `70fd11af9932f84523334119c3ca9a183aa3c03b`).

**Verified content.**

- Chen 2010: **Tier B** (OpenAlex abstract, publisher-elided from S2). "There
  is a great flexibility in partitioning the space of induced current into two
  orthogonal complementary subspaces: the signal subspace and the noise
  subspace", with part of the contrast source fixed by spectrum analysis and
  the rest optimized. Minor metadata conflict: OpenAlex records 2009 for this
  TGRS DOI; Crossref/S2 record 2010 (early-access vs issue year).
- Zhong–Chen 2009: **Tier C** only. The abstract is publisher-elided and the
  IOP page was not fetchable in this pass. The exact twofold structure
  `j− ≈ (I − V_{S,+}V_{S,+}*) V_{D,+} β` is attributed in `A3.md` to the 2018
  author monograph Eq. (6.71), not to the 2009 PDF: **[prior-pass]**,
  unresolved here.

**Overlap.** SOM signal/noise current-space decomposition and the twofold
sequential-projection seed are prior art. The span identity
`span[V_S,(I−P_S)V_D] = span[V_S,V_D]` is a main-thread derivation about
*numerical bases*, not a literature finding.

---

## T5. Goal-oriented / primal–adjoint / tangential ROM

**Sources and tiers.**

- T. Bui-Thanh, K. Willcox, O. Ghattas, B. van Bloemen Waanders,
  "Goal-oriented, model-constrained optimization for reduction of large-scale
  systems," J. Comput. Phys. 224(2):880–896, 2007,
  DOI 10.1016/j.jcp.2006.10.026. **Tier B** (abstract via zbMATH
  Zbl 1123.65081): "Optimization-oriented reduced-order models should target a
  particular output functional, span an applicable range of dynamic and
  parametric inputs, and respect the underlying governing equations of the
  system."
- R. Becker, R. Rannacher, "An optimal control approach to a posteriori error
  estimation in finite element methods," Acta Numerica 10:1–102, 2001,
  DOI 10.1017/S0962492901000010. **Tier B** (OpenAlex abstract): duality-based
  "local sensitivity factors" replace global stability constants so that
  "a posteriori estimates can be derived directly for the error in the target
  quantity" (the canonical dual-weighted-residual / goal-oriented estimator).
- K. Gallivan, A. Vandendorpe, P. Van Dooren, "Model Reduction of MIMO Systems
  via Tangential Interpolation," SIAM J. Matrix Anal. Appl. 26(2):328–349,
  2004, DOI 10.1137/S0895479803423925. **Tier B** (OpenAlex abstract): builds
  "a reduced order system of minimal McMillan degree that satisfies a set of
  tangential interpolation conditions", "a generalization of the multipoint
  Padé technique which is particularly suited to handle multiinput multioutput
  systems." Year note: Crossref/zbMATH print 2004, S2 prints 2005.
  First-order (Hermite) derivative matching is part of multipoint Padé /
  tangential interpolation but is not spelled out in the abstract:
  **[inference]**, flagged rather than asserted.

**Overlap.** `A3.md`'s correction is supported at Tier B across all three
ingredients: a ROM comparator suite must include goal-oriented,
primal/dual-weighted output-error, and tangential-interpolation methods, not
state-only POD.

---

## T6. Gain-invariant closure statistics (radio-interferometry lineage)

**Sources and tiers (all identity-verified via Crossref; abstracts via
OpenAlex).**

- R. C. Jennison, "A Phase Sensitive Interferometer Technique ...", MNRAS
  118:276–284, 1958, DOI 10.1093/mnras/118.3.276. **Tier B**: a triple
  interferometer "uniquely" determines the amplitude and phase of the complex
  Fourier transform of the brightness distribution (origin of the closure
  phase).
- D. H. Rogstad, "A Technique for Measuring Visibility Phase with an Optical
  Interferometer in the Presence of Atmospheric Seeing," Appl. Opt.
  7(4):585–588, 1968, DOI 10.1364/AO.7.000585. **Tier B**: triple-baseline
  closure phase removes atmospheric phase corruption.
- A. C. S. Readhead, P. N. Wilkinson, "The mapping of compact radio sources
  from VLBI data," ApJ 223:25–36, 1978, DOI 10.1086/156232. **Tier B**: hybrid
  maps from closure phase and fringe amplitude without knowing station gains.
- T. J. Cornwell, P. N. Wilkinson, "A new method for making maps with unstable
  radio interferometers," MNRAS 196:1067–1086, 1981,
  DOI 10.1093/mnras/196.4.1067. **Tier B**: telescope-specific gain-error
  correction (self-calibration family).
- T. J. Pearson, A. C. S. Readhead, "Image Formation by Self-Calibration in
  Radio Astronomy," Annu. Rev. Astron. Astrophys. 22:97–130, 1984,
  DOI 10.1146/annurev.aa.22.090184.000525. **Tier B**: "After first defining
  closure phase and closure amplitude ..." then reviews self-calibration
  algorithms and their limitations.

**Verified content vs. formal claim.** The abstracts establish that closure
quantities exist, are formed from sums of baseline observables, and drive
calibration-insensitive imaging. The explicit gain-cancellation identities
(a closed triangle of antenna gains dropping out of the bispectrum) are the
definitional property of this literature but were **not** re-derived from full
text in this pass; treat that step as [prior-pass]/textbook background, not as
newly verified Tier A evidence here.

**Overlap.** Gain-invariant closure cycles are long-established. The project's
rectangular full-wave transfer-matrix specialization and tangent checks are a
bounded specialisation, consistent with `A3.md`'s "可选工具，不新增品牌名称".

---

## Residual gaps (explicit)

1. Bellomo 2014 Eq. (25) / Section IV-C phase-center correction: **[prior-pass]
   only** — full text not reopenable in this pass (IEEE closed; HAL file
   absent/bot-gated).
2. Zhong–Chen 2009 full text and abstract: inaccessible; twofold structure
   rests on the 2018 monograph Eq. (6.71) **[prior-pass]**.
3. Bui-Thanh 2007 abstract came from a zbMATH/catalogue snippet, not the
   publisher page (S2/ScienceDirect elided).
4. Weiss 2026 is abstract-level (Tier B); MDPI full text was bot-gated, and the
   Europe PMC full-text XML for PMC 13469006 was not yet retrievable.
5. Closure gain-invariance is supported as established formalism (Tier B
   abstracts) rather than re-derived from primary full text here.
6. All six `scholarqa-cli collect` bundles are DNS-failed with zero candidates;
   nothing about search recall beyond these targets can be claimed from them.

