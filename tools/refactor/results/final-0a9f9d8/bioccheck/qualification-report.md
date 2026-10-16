# Actual final candidate BiocCheck qualification

The actual `0a9f9d8` check completed with **0 errors, 1 warning and 6 notes**, compared with **0/1/5** in the unchanged original `21f9ea0` actual check. The original run was reused after its pinned checker/input provenance was reverified; it was not rerun. Historical `32b2076` remains separately preserved at **0/1/7**. This is completed diagnostic qualification, not Bioconductor submission approval.

The final candidate is `0a9f9d879db6a52e3daef237daa1711b9d542976`. Its tarball SHA256 is `e442622f397248ba40788ef088b55c393d616a07602044c8cf0c82cbc40e5413`. Its successful full R CMD check record binds the same tarball, returns zero, ends `Status: OK`, and includes vignette rebuilding. The original and candidate preclean native compile flags match. All current source Git blobs, packaged R/src/NAMESPACE bytes and DLL/tarball hashes were reverified.

Both actual checks used the same unmodified R driver and isolated **BiocCheck 1.48.1**, **R 4.6.0**, **Bioconductor 3.23**, BiocManager 1.30.27 and BiocVersion 3.23.1. Exact commands, runtime versions and flags are included. All planned checks were enabled, including new-package diagnostics. The documented `install="check:<preclean install.log>"` and explicit library reuse the verified fresh installed package. `quit-with-status=FALSE` enables result serialization; it does not erase findings. No R packages were installed or updated for this candidate requalification.

## Findings and dispositions

The warning is inherited: no Bioconductor dependencies were detected. Adding a nominal dependency to silence it is outside this refactor. The inherited notes concern the R minimum-version suggestion, unverifiable Bioc-devel subscription, existing cat/print output, function lengths and long lines. No membership or submission-readiness claim is made.

One added NOTE remains: the nested `paste0` progress-message construction at `R/execution.R:41`. Original code precomputed `log_file_pattern` at `R/medium_level_functs_v2.R:218/380`. This remains an explicitly deferred, nonblocking source-hygiene follow-up; no additional source fix was made here.

The historical indentation NOTE is absent. Reported long lines are **145 original → 472 historical32 → 146 final**. All 146 final locations are inventoried: 119 inherited/moved handwritten exact trimmed lines, 8 inherited generated lines, 11 changed canonical generated binding lines and 8 introduced/rewritten handwritten lines. Their kinds are 109 code/literal, 18 roxygen and 19 generated. Exact trimmed-text matching describes text reuse, not behavioral equivalence; no remaining location is hidden. Preserve diagnostic text and canonical Rcpp generation in any later cleanup.

Functions over 50 physical lines are **48 original → 41 historical32 → 43 final**; the longest are **283 → 274 → 276**. The formatter changes physical-line metrics, which are guidance rather than independent proof of maintainability or numerical correctness.

## Evidence and scope

The bundle includes raw RDS/JSON results, full checker logs, exact commands/environment, all finding dispositions and formatting locations, verified fresh build/check records, pinned checker/dependency setup, reused original raw results and the immutable historical manifest. No production or user/system-library edits were made. This lane is local macOS arm64 qualification; the coordinator owns other platform/resource/scientific acceptance.

Official references: [release checker](https://bioconductor.org/packages/release/bioc/html/BiocCheck.html), [supported release/R pair](https://bioconductor.org/install/), [checker interpretation/options](https://bioconductor.org/packages/release/bioc/vignettes/BiocCheck/inst/doc/BiocCheck.html).
