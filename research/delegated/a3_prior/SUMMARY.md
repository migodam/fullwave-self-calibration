# A3 prior-work evidence extraction — summary

Status: task finished with explicit gaps; no global novelty judgment made
(that stays with the main Codex thread).

## What is now verified

1. Idriss–Raj arXiv:2503.07316v2, Eq. (8): the complex calibration factor is
   indexed `λ_c^{k,p}` (frequency × transmitter) inside a K×P MFSOM objective,
   and the paper says the `k`-dependence is dropped only "for notational
   clarity". Tier A. The A3 claim is confirmed verbatim.
2. Idriss–Raj calibrates per transmitter on the Fresnel measured dataset
   (Tier A, abstract of the same full text).
3. Weiss et al., Sensors 26(15):4954 (2026): joint receiver-geometry +
   clock-bias + target-position self-calibration tied to coherent
   matched-filter image sharpness, with identifiability analysis. Tier B
   (PubMed abstract). Confirms the A3 claim that this broad combination is
   occupied, and that it is not quantitative full-wave material inversion.
4. Chen 2010 SOM signal/noise current-subspace partition: Tier B.
   Zhong–Chen 2009 TSOM identity: Tier C; its equation-level twofold
   structure remains [prior-pass] via the 2018 monograph Eq. (6.71).
5. Goal-oriented ROM (Bui-Thanh et al. 2007), primal–dual/dual-weighted
   residual estimation (Becker–Rannacher 2001), tangential-interpolation MIMO
   reduction (Gallivan et al. 2004): all Tier B. The A3 correction that a
   fair ROM comparator must go beyond state-only POD is supported.
6. Closure-phase/closure-amplitude lineage (Jennison 1958; Rogstad 1968;
   Readhead–Wilkinson 1978; Cornwell–Wilkinson 1981; Pearson–Readhead 1984):
   Tier B as established calibration-insensitive observables. The project's
   rectangular-MIMO cycles are a specialisation, not an invention.
7. Bibliographic identity: `verify_run2.json` resolved 7/8 canonical
   identifiers; the remaining 6 records were resolved via Crossref/OpenAlex/
   zbMATH/PubMed.

## Verification gaps (marked, not hidden)

- Bellomo 2014 Section IV-C / Eq. (25) phase-center detail: **[prior-pass]
  only**. Full text unreachable this pass (IEEE closed; HAL file absent /
  bot-gated; no local copy). The abstract supports antenna calibration in
  diffraction tomography but does not mention phase centers.
- Zhong–Chen 2009 abstract/full text inaccessible; twofold structure rests on
  the monograph, not the 2009 PDF.
- Weiss 2026 and Bui-Thanh 2007 are abstract-level; Bui-Thanh's abstract came
  from zbMATH rather than the publisher page.
- The six `scholarqa-cli collect` bundles are DNS-failed with zero candidates,
  so no recall claim beyond the named targets is possible.
- Closure gain-invariance is cited as established formalism (Tier B), not
  re-derived from full text here.

## Artifacts

- `EVIDENCE.md` — per-target methods/advantages/limitations/overlap and gaps.
- `CLAIM_LEDGER.md` — 16 bounded claims with tiers, excerpts, relationship,
  confidence.
- `REFERENCES.md` — 13 verified records (stable IDs only; no fabricated
  BibTeX).
- `retrieval_*.json` — six preserved ScholarQA collect bundles (all
  DNS-failed, zero candidates).
- `verify_run1.json` (failed), `verify_run2.json` (7/8 resolved).

## Where Codex intervention is needed

- Re-open the Bellomo author manuscript (or a library copy) to confirm
  Section IV-C / Eq. (25) before any publication sentence rests on it.
- Decide whether Tier B closure evidence suffices for the A3 "可选工具"
  wording, or whether a closure full text should be re-derived.
- Keep the final novelty judgment in the main thread; these reports only
  bound what the retrieved neighbours already occupy.

