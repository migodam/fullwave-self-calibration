A2 literature audit complete. All three deliverables plus the raw-data mirror are in `research/delegated/a2_literature/`.

Findings (worker pass; novelty judgment stays with you):

1. ScholarQA verify ran once as instructed: 8/8 records resolved, no 429 this time. All DOIs independently cross-checked against Crossref/arXiv.
2. Xu/Wu/Ye/Chen TAP2020 DOI resolved: 10.1109/TAP.2020.2998171, TAP 68(11):7457-7470; fulltext read from the NUS author PDF. It solves phaseless ISPs with current-based schemes (PD-DICs, PD-CSI) — supports the "phaseless still has structured currents" disproof, does no geometry calibration.
3. Idriss/Raj 2503.07316 v2 fulltext read: MFSOM plus jointly estimated **complex per-transmitter** calibration factors on Fresnel measured data (Eqs. 6-8, Fig. 6(a)). This is the strongest closest overlap — an "SOM + calibration exists" claim is not tenable; only the typed rank/bias/pose certificates remain distinguishable.
4. Bellomo 2014 confirmed (TAP 62(5):2450-2462); HAL record has no file, so the phase-center Eq. (25) detail stays unverified prior-pass.
5. Huang/Nammour/Symes 2017 confirmed (Geophysics 82(3):R153-R171); fulltext still blocked, abstract publisher-elided.
6. Wang/Paulus/Eibert 2605.03875 confirmed as EuCAP 2026 (10.23919/EuCAP68105.2026.11612585); abstract-level inverse-source vs inverse-scattering distinction.
7. SOM application examples confirmed: Lu et al. 2010 (PIER 102:351-366, abstract) and Meng et al. 2016 (Sensors 16(7):1046, fulltext via PMC4970093).
8. Measured-data gate MET: official 2001-opus TM long-dielectric-cylinder files (`dielTM_dec8f.exp` etc.) downloaded from IOP-served public links on the canonical domain, no CAPTCHA. Hashes and format audit in the files below.
9. 2005-opus TM foam files also mirrored, but via the IOP beta domain because the canonical S09 data page is Radware bot-gated; flagged as secondary. Old `fresnel.fr` data paths are 404 and the old iop.org 2005 link is Cloudflare-challenged (not bypassed).
10. No source supports any "first" claim at the level of SOM, phaseless SOM, self-calibration, or antenna calibration; the evidence supports only the conditional recommendation scoped to the bias/rank/acquisition certificates.

Deliverables:

- [EVIDENCE.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_literature/EVIDENCE.md)
- [references.bib](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_literature/references.bib)
- [public_data_audit.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_literature/public_data_audit.md)
- [data/SHA256SUMS.txt](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_literature/data/SHA256SUMS.txt) and [data/README.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_literature/data/README.md)

Uncertainties needing your call: Bellomo phase-center content could not be re-opened (prior-pass only); Huang 2017 fulltext remains blocked; Chen/Zhong/Pan abstracts are publisher-elided so equation-level attribution rests on the 2018 monograph; the 2005 data route via the IOP beta mirror needs your sign-off before treating it as first-class evidence (the 2001 canonical route is clean and already satisfies the gate).