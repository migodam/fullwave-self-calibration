# Full-Wave Self-Calibration · A1–A5 / Q5 research archive

## Latest public snapshot: Q5, 10 September 2026

The automatic verification has now ended and is integrated into the existing
Q5 manuscript, theory, Chinese report and GPT prompt. The same reading entry and
ZIP below replace the earlier pre-completion snapshot; no second packet is needed.

**[Start here: Q5 ChatGPT reading package](Q5_CHATGPT_START_HERE.md)** — theory,
English working manuscript, Chinese report, next research prompt, source code,
registered experiment parameters and immutable numerical results.

Q5 is scientifically incomplete. Two development cycles are complete; the
candidate-selection advantage is withdrawn. Independent-model recovery was
1/12 without reference and 6/12 with reference; the new policy reached 3/12.
No final 128-scene test was generated. Supplied A5, parent reconstructions and
provisional worker calculations are kept distinct. The older A1–A4 entries
below remain available as history, not the latest verdict.

Coherent electromagnetic inverse scattering, TriSpace SOM, map–pose coupling,
unknown arrays, geometry/electronic calibration, and model discrepancy.

**Working research, not a submission-ready or accepted TAP paper.** This archive
preserves theories, implementations, experiments, corrections and failed routes.
Historical positive claims must be read alongside later audits. The supplied
A4 package is included intact; packaging checks are not independent validation
of all scientific claims.

## 给 GPT / ChatGPT 的入口

Start with **[GPT_READING_GUIDE.md](GPT_READING_GUIDE.md)** for reading order and
direct raw-file links. Use [RESEARCH_NAVIGATION.md](RESEARCH_NAVIGATION.md) for
the A1→A4 timeline and SOM/SLAM topic map. Every exported file is enumerated in
[FILE_INDEX.md](FILE_INDEX.md), with source hashes in [PUBLIC_MANIFEST.json](PUBLIC_MANIFEST.json).

| Need | Entry |
| --- | --- |
| Latest supplied A4 English manuscript | [A4 paper](research/a4_reliability_v1/PAPER_DRAFT_A4.md) |
| A4 theory, proofs, code and experiments | [A4 package](research/a4_reliability_v1/README.md) |
| Last parent-audited A3 English manuscript | [A3 paper](research/trispace_self_calibration/a3_research/PAPER_DRAFT_A3.md) |
| GPT Pro total prompt through A3_3 | [Remaining theory/design questions](research/trispace_self_calibration/a3_research/GPT_PRO_REMAINING_QUESTIONS_ZH.md) — historical input to A4, not a post-A4 rewrite |
| A4 input question | [Q4](Theory/Questions/Q4.md) |
| 中文汇报 | [A3 report](communication/A3_RESEARCH_REPORT_ZH.md), [A2 report](communication/A2_RESEARCH_REPORT_ZH.md) |
| Original theory context | [SOM / SLAM context](Theory/SOM_SLAM_THEORY_CONTEXT.md) |
| Original Q/A documents | [Theory/Questions](Theory/Questions) |
| Gaussian research workspace | [Gaussian A1 + A1_2 + A2](Gaussian/deliverables/START_HERE.md) |
| Code entrypoints and dependencies | [Code and reproducibility](CODE_AND_REPRODUCIBILITY.md) |
| Boundaries and checks of this upload | [Publication scope](PUBLICATION_SCOPE.md) |

## Layout

```text
Theory/                         Original research context and Q/A inputs
Gaussian/                       Gaussian A1/A1_2/A2 theory, code and evidence
communication/                  Chinese research reports
research/a4_reliability_v1/      Supplied A4 package, source and results
research/trispace_self_calibration/
  a2_research/                  A2 paper, theory audit, analysis and results
  a3_research/                  A3 physics/solver code, protocols and results
research/delegated/              Supporting implementations, tests and reviews
experiments/idea_loops/          Historical autonomous idea-validation iterations
```

Code and numerical results are now included. This is not a single installable
application: historical paths may refer to local environments, excluded solver
dependencies or original measured datasets. Credentials, private user files,
raw model transcripts, virtual environments, third-party article PDFs and
bundled external solver source are not published. No open-source license is
assigned by this upload; third-party rights are not overridden.
