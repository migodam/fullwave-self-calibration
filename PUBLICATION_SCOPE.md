# Public archive scope and release checks

Date: 9 September 2026. User-authorized public expansion of the existing repository.

## Included

Authored theory notes and supplied Q/A documents, English working papers, Chinese
reports, Python/shell source, requirements, protocols, numerical JSON/JSONL,
synthetic arrays and generated figures. A1–A3 idea loops include outcome records,
parent reviews and invalidated implementations so failure history is not erased.

The entire 66-file A4 package is expanded from the user-supplied archive into
`research/a4_reliability_v1/`, without executing its code during extraction or
changing its frozen sources/results. Top-level A4 documents in Theory/Questions
are also retained as supplied; packaged documents use their original relative
figure paths and should be preferred for reading.

## Excluded

- Private `User/` content, API keys/configuration and environment files.
- Virtual environments, caches, binaries and heavy solver work products.
- Raw model transcripts / worker logs and raw literature retrieval bundles.
- Third-party article PDFs, external solver source and original measured datasets
  whose redistribution terms have not been reviewed.
- Files outside the explicit research roots/type allowlist; files above32MiB.

Source paths and provider names mentioned in authored audit reports can remain as
provenance. Their presence does not supply a credential or reproduce the original
machine environment. Historical references to excluded files are not necessarily
resolvable from this clone.

## What was checked in this publication pass

- Scoped export checked for common secret/token/private-key patterns and raw
  transcript markers; no candidate files were flagged by those checks.
- Source bytes are preserved and hashes are listed in PUBLIC_MANIFEST.json.
- A4's read-only package audit passed frozen source hashes, result counts,
  stored parameter-loss relationships, cost accounting and raw-array checks.
- A4 pytest could not start in the available A3 environment (pytest absent).
  The package's reported 21 passes are historical, not a new successful test run.
- New root navigation links are checked locally; historical documents can retain
  missing local/dependency references. The GitHub remote commit is checked after push.

These checks do not certify scientific novelty, continuum Maxwell accuracy,
global recovery, hardware calibration, or TAP readiness. The new A4 theorem
and experiments still require substantive independent review. Publishing source
does not upgrade a supplied result to an independently replicated result.

No open-source license is assigned by this upload. References identify prior work;
they do not transfer third-party rights. A later licensed software release can be
prepared after ownership and dependency review.
