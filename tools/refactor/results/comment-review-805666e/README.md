# Comment readability review

This follow-up rewrites internal comments without changing package behavior.
It starts from `805666e` (the illustrated-tour checkpoint), which already carries
all refactor qualification described in the handoff. Source bytes now differ
because comments changed; the earlier qualification remains tied to its named
commits, not to a claim that current sources are byte-identical.

The review covers handwritten package R/C++/headers, tests, permanent diagnostic
scripts, bundled data scripts and vignette code comments. Clear existing comments
were kept. Generated Rcpp files and frozen historical comparison evidence remain
verbatim. The main changes replace internal jargon with the job being done,
identify what cell numbers refer to, add small examples, and explain why ordering,
fit reuse and error-cleanup rules matter. Statistical formulas and necessary
C++ safety details remain explicit.

Validation against the preceding checkpoint:

- All 51 changed R files have identical parsed executable expressions, using
  `parse(..., keep.source = FALSE)` and `identical()`.
- All 21 changed C++/header files have identical Clang raw tokens after excluding
  comments and whitespace. The parser consumes the complete Clang output.
- Rcpp attribute comments and public roxygen documentation blocks are unchanged;
  only internal `@noRd` descriptions were rewritten. Generated bindings,
  NAMESPACE, DESCRIPTION and CI configuration are unchanged.
- Checker controls accept comment-only changes and detect changed numbers and
  strings in both languages.
- `air format --check R` and `git diff --check` pass.
- The walkthrough renders with Quarto. Its result-order explanation now says
  failed-QC rows are added back and results sorted by p-value/response, matching
  the code; the two affected diagrams were regenerated and visually inspected.

`validation.json` records per-file hashes, token counts, review scope and checks.
`verification-scripts.tar.gz` retains the exact local verification scripts and log.
The scripts use the recorded Mac checkout/temporary paths and baseline; adjust
those paths if reproducing elsewhere. No numerical, resource or full package test
rerun is claimed for this comment-only change.
