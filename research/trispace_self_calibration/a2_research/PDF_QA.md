# PDF rendering QA

Final review PDF: `output/pdf/trispace_a2_research_draft.pdf`, 12 A4 pages.
SHA-256: `f6bfe5441d0d9a49a0654906dd2d682ff378664cea007ba3c699a4ab9c4009e0`.
Markdown source hash is recorded in `output/pdf/build_manifest.json`.

Generated with Tectonic/XeTeX and inspected via Poppler PNG rendering, not just
text extraction. Initial build issues (Python delimiter and unavailable font)
were corrected before PDF creation succeeded. The final build reports no
overfull/underfull boxes, missing characters, or undefined references. Absolute
local font/image path warnings are retained as portability limitations.

Inspected manuscript text, all equation-heavy pages, four tables, both figures,
appendices and reference page. First full render had 13 pages; the final revision
integrates references after Appendix D and has 12. Latest affected pages and
equation/table pages were re-rendered and visually checked. No clipped formulas,
overlapping table text, missing panels or placeholder citations were observed.
The scientific figure has a shared contrast color scale, explicit exploratory
seed, and a known-pose reference label, not an "upper bound" claim. The all-zero
recovery graph uses the full [0,1] probability axis.

Review-format decisions: single-column equations, readable tables, long proofs
in Appendices A--D, a dedicated reconstruction figure page, sequential page
numbers and article-style references. This is not an IEEE submission template
or accessibility-tagged final publication PDF. Visual QA establishes rendering
quality only; scientific acceptance is separately recorded in FINAL_REVIEW.md.
