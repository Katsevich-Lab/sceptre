# Final RunPlan and unified FitCache qualification — 2026-10-19

Implementation: `e545b0891292ba82399271ec35c2e1f6818d0cbe`, following the
separate typed-RunPlan checkpoint `e86eb48`. The working saved schema remains 3,
with explicit adapters for old list plans and physical response-only cache
layouts. [Exact current field tables](../schema3/OBJECT_FIELDS.md) are generated
from the actual S4 class definitions.

| Check | Terminal result |
|---|---|
| `R CMD build .` | Passed, including vignette generation. |
| `R CMD check --no-manual --no-build-vignettes` | **Status: OK**, including examples, installed tests and vignette R code. |
| Final installed-package tests | **5,069 passed; 0 failures, warnings or skips.** |
| Source-suite checkpoint | 205 cases, 5,069 expectations; 0 failures/errors/warnings/skips. |
| Matched-budget reference comparison | Seven workflows, 231 tolerance-zero comparisons; all equal. |
| BiocCheck 1.47.26 | **0 errors, 1 warning, 6 notes.** |
| Independent integration review | No actionable findings; [scope and report](INDEPENDENT_REVIEW.md). |

Environment: R 4.6.0, macOS arm64, ondisc 1.2.0. The exact final archive checksum,
158 source-file comparisons, terminal check logs, test output and JSON totals
are in [validation/](validation/). The source-suite and parity checkpoint was
frozen before one nonnumerical namespace qualification of the RunPlan slot
setter. Its original hashes are preserved. The rebuilt archive and final
installed-package suite cover the exact final source. There are no native-source
changes.

The final check does not include manual PDF generation or rebuilding vignette
outputs during check; the build generated them and check executed their R code.
BiocCheck's R-version recommendation and interactive bioc-help check were
disabled. The warning concerns the existing lack of Bioconductor infrastructure
dependencies. Notes concern console output, condition messages, nonlocal
assignment, function length, line width and indentation; these are not claimed
resolved.

[The numerical evidence](validation/README.md) distinguishes equal-budget
reference comparisons from the old calibration-budget overwrite. New cache tests
also count actual fitting calls: warm runs perform no perturbation refits while
still generating full CRT banks. Real low-MOI NT fixtures exercise both shared
and per-entry incompatible-context misses. Other tests cover independent family
invalidation, compact entries, single-fit unsupported-model fallback, legitimate
warnings, transactional failure, saved cache/plan migration and serialization.

[Eugenio Mattei's PR #233](PR233.md) was inspected by title, author, pinned diff,
body and inline review. Its compact exact reconstruction, independent invalidation,
calibration bypass and fresh draws were adopted. Current ownership uses one
FitCache with shared context and per-group input hashes, without repeated NT
indices. Aliased/nonfinite or inexact fits use their original probabilities and
remain uncached. Actual-fit warnings are retained as diagnostics and are not
replayed on warm hits. [The RunPlan assessment](../runplan/DECISION.md) records why
S4 is a reasonable consistency choice rather than an automatic correctness or
immutability improvement.

Qualification is local to the recorded platform/runtime. It establishes no
large-screen speed or memory benchmark, cross-platform arithmetic portability,
all-historical-object guarantee or Nextflow compatibility. Active backing files
remain immutable by contract; no external-tamper monitoring was added. No push,
PR publication, merge or release was performed. The later book-grounded audit is
outside this completed implementation scope.
