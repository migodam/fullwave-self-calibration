# A2 literature audit: bounded evidence and measured-data audit

Worker: isolated DeepSeek Pro pass, 6 September 2026.
Scope source: `Theory/Questions/A2_PRASC_SOM_THEOREM_PACKAGE_EN.md` and
`A2_PRASC_SOM_VALIDATION_PROTOCOL.md` (intended claim scope only; they are not
authoritative literature). Final novelty/theory judgment stays with the main
Codex thread.

Evidence-level labels used throughout:

- `[fulltext]` — a primary document (PDF/HTML) was actually read in this pass
  or is inherited from the prior verified pass and re-checked here.
- `[abstract]` — the publisher/indexed abstract was read.
- `[metadata]` — bibliographic identity only (title/authors/venue/year/DOI),
  resolved from Crossref, arXiv, Semantic Scholar, or HAL.
- `[inference]` — a conclusion drawn by this worker from the above; not a
  statement contained verbatim in a source.
- `[prior-pass]` — retained from the theorem package's own access notes,
  flagged because it was not re-derived in this pass.

## 1. Intended claim scope (as declared in the theorem package)

The manuscript claims a geometry self-calibration framework that preserves
coherent phase, residual-certified physical states with `r_free = 0`, a
SOM-informed reduced basis that is not the classical independent-current SOM,
and joint map/pose diagnostics with limitations. Two stated negative targets:

1. broad phaseless claims implying "no current constraint" is impossible or
   that a current subspace necessarily absorbs all phase information;
2. claims that no prior self-calibration/calibration work exists.

Neither the theorem package nor this audit claims novelty ("first") anywhere.
The verdict the audit can support is a **conditional recommendation**: TGRS- and
TAP-scoped primary sources show dense prior art on SOM, phaseless SOM, and
antenna/calibration methodology, so any publication claim must be scoped to
the typed bias/rank/acquisition certificates, not to SOM-with-calibration as a
whole. Nothing in this audit establishes a universal basin, an exclusive SOM
advantage, or superiority over a matched direct joint solver.

## 2. Verification method

- ONE ScholarQA `verify` attempt (as instructed), 8/8 records resolved, no 429
  error this time. Output is bibliographic identity only.
- Canonical fallbacks: Crossref REST API, arXiv export API, HAL API, Europe
  PMC, IOPscience supplementary-data pages, and one author-institution PDF.
- No citation was invented; every DOI in `references.bib` was resolved through
  at least one canonical registry.

## 3. Primary publication metadata ledger

| Reference | Verified identity | Evidence level |
|---|---|---|
| Xudong Chen, "Subspace-Based Optimization Method for Solving Inverse-Scattering Problems," IEEE TGRS 48(1), 42-49, 2010. DOI 10.1109/TGRS.2009.2025122 | Crossref + ScholarQA agree | `[metadata]`; equation-level SOM claims rely on the 2018 monograph `[prior-pass]` |
| Y. Zhong, X. Chen, "Twofold subspace-based optimization method for solving inverse scattering problems," Inverse Problems 25, 085003, 2009. DOI 10.1088/0266-5611/25/8/085003 | Crossref + ScholarQA agree | `[metadata]`; original PDF retrieval previously failed `[prior-pass]` |
| Li Pan, Yu Zhong, Xudong Chen, Swee Ping Yeo, "Subspace-Based Optimization Method for Inverse Scattering Problems Utilizing Phaseless Data," IEEE TGRS 49(3), 981-987, 2011. DOI 10.1109/TGRS.2010.2070512 | Crossref + ScholarQA agree | `[metadata]`; abstract screened previously `[prior-pass]` |
| K. Xu, L. Wu, X. Ye, X. Chen, "Deep Learning-Based Inversion Methods for Solving Inverse Scattering Problems With Phaseless Data," IEEE TAP 68(11), 7457-7470, 2020. DOI 10.1109/TAP.2020.2998171 | Crossref + PDF from NUS author page (2.0 MB) agree | `[fulltext]` — PDF read in this pass |
| Z. Idriss, R. G. Raj, "Data-Driven Calibration Technique for Quantitative Radar Imaging," arXiv:2503.07316 (v2), 2025; IEEE CISA 2025 | arXiv API + full HTML read | `[fulltext]` — v2 HTML read in this pass |
| L. Bellomo, S. Pioch, M. Saillard, K. Belkebir, "An Improved Antenna Calibration Methodology for Microwave Diffraction Tomography in Limited-Aspect Configurations," IEEE TAP 62(5), 2450-2462, 2014. DOI 10.1109/TAP.2014.2308534 | Crossref + ScholarQA + HAL record agree; HAL has metadata only, no file | `[metadata]`; phase-center Eq. (25) claim is `[prior-pass]` and not re-verifiable here |
| G. Huang, R. Nammour, W. Symes, "Full-waveform inversion via source-receiver extension," Geophysics 82(3), R153-R171, 2017. DOI 10.1190/geo2016-0301.1 | Crossref + ScholarQA agree | `[metadata]` + publisher-blocked abstract; full PDF blocked `[prior-pass]` |
| Q. Wang, A. H. Paulus, T. F. Eibert, "Phase-Corrected Near-Field Microwave Imaging via Inverse Source Reconstruction with Modulated Signals," arXiv:2605.03875, 2026; EuCAP 2026, DOI 10.23919/EuCAP68105.2026.11612585 | arXiv API + Crossref + ScholarQA agree | `[abstract]` |
| T. Lu, K. Agarwal, Y. Zhong, X. Chen, "Through-Wall Imaging: Application of Subspace-Based Optimization Method," PIER 102, 351-366, 2010. DOI 10.2528/PIER10020903 | Crossref + ScholarQA agree | `[abstract]` |
| Q. Meng, K. Xu, F. Shen, B. Zhang, D. Ye, J. Huangfu, C. Li, L. Ran, "Microwave Imaging under Oblique Illumination," Sensors 16(7), 1046, 2016. DOI 10.3390/s16071046 | Europe PMC + Crossref agree | `[fulltext]` — PMC4970093 read in this pass |

Fresnel measured-data primary records (for the data audit):

- K. Belkebir, M. Saillard, "Special section: Testing inversion algorithms
  against experimental data," Inverse Problems 17, 1565-1571, 2001.
  DOI 10.1088/0266-5611/17/6/301 `[metadata]`
- K. Belkebir, M. Saillard, "Testing inversion algorithms against experimental
  data: inhomogeneous targets," Inverse Problems 21, S1-S3, 2005.
  DOI 10.1088/0266-5611/21/6/S01 `[metadata]`
- J.-M. Geffrin, P. Sabouroux, C. Eyraud, "Free space experimental scattering
  database continuation: experimental set-up and measurement precision,"
  Inverse Problems 21, S117-S130, 2005. DOI 10.1088/0266-5611/21/6/S09
  `[metadata]` (Crossref confirms DOI S09 -> pages S117-S130)
- J.-M. Geffrin, P. Sabouroux, "Continuing with the Fresnel database:
  experimental setup and improvements in 3D scattering measurements," Inverse
  Problems 25, 024001, 2009. DOI 10.1088/0266-5611/25/2/024001 `[metadata]`

## 4. Closest overlaps, with what each does and does not establish

### 4.1 SOM / TSOM lineage (Chen 2010; Zhong-Chen 2009)

`[metadata]` here, `[prior-pass]` equation-level. Chen 2010 is the original
deterministic/ambiguous current decomposition; Zhong-Chen 2009 is TSOM. Both
predate this manuscript. `[inference]` Any SOM-based reduced basis must be
distinguished from these; the theorem package's separation of `L_det`,
`r_num`, and `r_free` is compatible with, and motivated by, this lineage. No
prior-art gap exists at the level of "SOM decomposition exists."

### 4.2 Phaseless SOM retains a current constraint (Pan et al. 2011)

`[metadata]` here, `[prior-pass]` abstract/monograph. The monograph account
(Eqs. 8.28-8.30) shows phaseless SOM keeps the deterministic/ambiguous
partition and uses total-field intensity. `[inference]` This directly opposes a
broad "phaseless means no current constraint" reading: the published phaseless
SOM still has a data-determined major current. It does not speak to the
manuscript's matched-information inequality (Theorem 1), which is a separate
statistical statement needing its own validation (protocol E1).

### 4.3 Learning phaseless inversion with current-based schemes (Xu et al. 2020)

`[fulltext]` NUS PDF read. The paper solves phaseless ISPs with U-net CNN
training schemes whose inputs include phaseless total field (DIS), dominant
induced currents retrieved by Levenberg-Marquardt (PD-DICs), and a phaseless
contrast-source inversion scheme (PD-CSI). Its main result is that phaseless
deep reconstruction outperforms its baselines for that benchmark. `[inference]`
It is evidence that phaseless inversion with structured currents is an active,
occupied area, but it performs no geometry self-calibration and no
information-theoretic guarantee; it does not pre-empt the manuscript's typed
certificates, nor does it support any "phaseless has no current structure"
claim.

### 4.4 Joint SOM + complex transmitter calibration already exists (Idriss-Raj 2025)

`[fulltext]` arXiv:2503.07316v2 HTML read in this pass. Section III uses the
multiple-frequency subspace-based optimization method (MFSOM): the current is
split into signal/noise subspaces by SVD of the Green's operator (Eq. 6-7),
and the objective (Eq. 8) jointly optimizes a calibration factor per
transmitter together with the current and contrast. Fig. 6(a) states the
factor is allowed to be **complex** "to account for phase factors" and that a
purely real factor fails; factors differ per transmitter; it is validated on
the Fresnel Institute measured dataset. `[inference]` This is the strongest
closest overlap: SOM-based quantitative inversion with jointly estimated
complex gain/phase calibration on measured data already exists. What it lacks
(and therefore what remains distinguishable) is the manuscript's
rank/bias/pose information budget, residual-certified physical states, and
geometry (not merely per-transmitter scale) calibration. No "first
SOM + calibration" claim is tenable.

### 4.5 Antenna calibration with phase-center handling (Bellomo et al. 2014)

`[metadata]` this pass. The publisher-visible abstract (search index) states
the work underlines properly modeling antenna directivity and proposes an
improved calibration strategy for microwave diffraction tomography. The
theorem package `[prior-pass]` cites author-manuscript Section IV-C and
Eq. (25) for explicit phase-center correction and constrained spectral/physical
calibration modes. This pass could not re-open the full text (IEEE closed
access, HAL record has no attached file), so the equation-level claim must be
re-verified before any publication sentence rests on it. `[inference]` At
minimum it establishes prior art on *antenna* calibration inside diffraction
tomography; it does not establish the manuscript's pose/geometry certificates.

### 4.6 Source-receiver extension (Huang-Nammour-Symes 2017)

`[metadata]` this pass; publisher abstracts are elided and the full PDF is
blocked `[prior-pass]`. Source-receiver extension is the standard
conceptual reference for joint data/geometry extension in FWI.
`[inference]` The manuscript's "joint map/pose diagnostics" should be
positioned against extension/relocalization methods rather than claimed as
unoccupied. Equation-level equivalence to the manuscript remains unverified
and is a residual gap.

### 4.7 Inverse source vs inverse scattering (Wang-Paulus-Eibert 2026)

`[abstract]` arXiv + EuCAP record. This work uses inverse source
reconstruction for passive near-field microwave imaging, with phase correction
for coherent multi-frequency superposition. It is an explicit, current example
of the inverse-source/inverse-scattering distinction the task asked to keep
separate. It does not address map/pose self-calibration.

### 4.8 SOM application examples

- Through-wall: Lu et al. 2010, PIER 102, 351-366 `[abstract]`. SOM with a
  layered-medium model for through-wall imaging; numerical only.
- Oblique illumination: Meng et al. 2016, Sensors 16(7), 1046 `[fulltext]`
  PMC4970093 read. Experimental 2D SOM under oblique illumination, with the
  stated finding that cross-polarization neglect costs little for 2D-like
  objects.

Both confirm SOM is established in application contexts; neither calibrates
geometry or proves information guarantees.

## 5. Support for the two intended disproofs

- **"Phaseless has no current constraint"** — Pan et al. 2011 (phaseless SOM
  keeps the major-current partition, `[prior-pass]` monograph) and Xu et al.
  2020 (phaseless schemes built around dominant induced currents, `[fulltext]`)
  together falsify the broad reading. They do not by themselves prove the
  manuscript's sharper statement that phaseless observation of a matched
  experiment never has *larger efficient information*; that needs protocol E1.
- **"No prior self-calibration"** — Idriss-Raj 2025 (joint MFSOM + complex
  per-transmitter calibration on measured data, `[fulltext]`) and Bellomo et
  al. 2014 (antenna calibration in diffraction tomography, `[metadata]` +
  `[prior-pass]` phase-center detail) close that gap. Any claim must be about
  the specific geometry/pose self-calibration with certified residuals, not
  about calibration existing at all.

## 6. Firstness and recommendation

No source supports a novelty claim at the level of SOM, phaseless SOM,
self-calibration, antenna calibration, or source-receiver extension. The only
defensible conditional recommendation is: **scope the contribution to the
typed rank/bias/acquisition information budget and its certificates**, and
test it against a matched direct joint solver and against Idriss-Raj-style
SOM-with-calibration. TGRS (imaging method) and TAP (calibration method) both
have explicit evidence gates in the validation protocol that are not met by
metadata or by elementary matrix checks.

## 7. Measured-data gate

MET for the 2001 opus TM long dielectric cylinders (official IOP-served
files), with hash/format/provenance in `public_data_audit.md`. The 2005 opus
TM foam targets were also downloaded, with an explicit access-route caveat
(canonical IOP domain returned a bot challenge; the IOP beta mirror served
public pre-signed links). No generated data was substituted and no inversion
was implemented. Details and SHA256 sums: `data/SHA256SUMS.txt` and
`data/README.md`.

## 8. Uncertainties and residual gaps

1. Bellomo 2014 equation-level phase-center content could not be re-opened;
   it remains `[prior-pass]`.
2. Huang-Nammour-Symes 2017 full text remains blocked; abstract is
   publisher-elided. Only metadata identity is fresh.
3. Chen 2010 / Zhong-Chen 2009 / Pan et al. 2011 abstracts are
   publisher-elided via the public API; equation-level attribution comes from
   the 2018 author monograph `[prior-pass]`.
4. The 2005 S09 data route on the canonical domain is bot-gated; the beta
   mirror was used and is recorded as such.
5. No abstract of this audit constitutes novelty search exhaustiveness; it is
   a bounded closest-overlap screen.

## 9. Deliverables

- `EVIDENCE.md` (this file)
- `references.bib` — verified BibTeX records
- `public_data_audit.md` — measured-data access, provenance, format, hashes
- `data/2001_iop_17_6_301/` — 8 `.exp` files (2001 opus)
- `data/2005_iop_21_6_S09/` — 8 `.exp` files (2005 opus, via IOP beta mirror)
- `data/SHA256SUMS.txt`, `data/README.md`
