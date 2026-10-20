# Independent candidate diagnostic-instance review

Candidate: `01bb39ccf8373ae8965a2679ee4813a63ea2d85d`.
Actual submitted baseline: `21f9ea098b69c444a884a49bc254de711968e3b5`.
Reviewed archive: `/private/tmp/sceptre-bioc-cleaned-01bb39c/sceptre_0.99.0.tar.gz`.
Verified SHA-256: `606363fce32238354fc0d126d79f515973da7c8287ca43b14ce53ce2d00df4b7`.

**The candidate passes this independent source-instance comparison: no new flagged function sites, long-line instances, display calls, or coding-practice categories were found.** This is a bounded scanner, source-lineage, and cleanup-quality review; it is not a claim of Bioconductor acceptance, numerical parity, or an independently rerun full package check.

## Independent evidence and version handling

The supplied archive was extracted under `sceptre/`. A separate comparison verified 180tracked package source/test/documentation files byte-for-byte against frozen 01bb39c. Built DESCRIPTION retains every original DCF field value after ordinary DCF whitespace normalization and adds standard build metadata only. See `archive-git-source-manifest.json` and `built-description-metadata.json`.

I ran focused scanner functions in separate fresh R processes using the unmodified installed BiocCheck 1.47.26 and isolated 1.49.30 packages. Their DESCRIPTION hashes, library paths, function-length scanner code/hashes, and session information are saved. I did not substitute or patch a scanner, apply the owner's allowlist, change check options, or run another full BiocCheck/R CMD check.

| Scanner | Submitted long functions | Candidate long functions | Submitted/candidate long lines |
|---|---:|---:|---:|
| BiocCheck 1.49.30, the original submission checker |49|45|145/108|
| BiocCheck 1.47.26, the previous local checker |48|44|145/108|

The discrepancy is now resolved exactly. The 52-line `StatQQBand$compute_group` method in `R/qq_plot_helpers.R` is present at submitted 36–87 and candidate 57–108, with 46 coding lines in both. The older scanner stores results keyed by function name and overwrites this site when it reaches another `compute_group`; 1.49.30 retains both sites. This is an inherited site and checker behavior, not a new function or an unexplained waiver. See `checker-version-instance-difference.json`, both complete version-specific inventories, and both saved `getFunctionLengths-*.R` implementations.

The owner-produced completed 1.49.30 full log was inspected and copied with its invocation metadata. It uses the original official action's `no-check-R-ver` and `no-check-bioc-help` settings and the exact candidate archive. Its source findings agree with this independent inventory: the inherited no-Bioconductor-dependency warning and three note categories (display calls, function lengths, line lengths); man-page examples and coding practice pass. The coordinator's separate all-default 1.47.26 comparison must retain environment/network/version/subscription outcomes distinctly. This report does not add or recommend any disabling flags.

## All 108remaining long lines: exact text and multiplicity

`long-line-lineage-independent.tsv` contains every candidate flagged line, full text, width, exact candidate path/line and submitted path/line/width. **Matches are exact complete source-line text, including indentation and whitespace; no whitespace normalization was used.** An additional raw-byte comparison confirms the 108 matches, including all leading/trailing whitespace, excluding line terminators only; see `long-line-raw-byte-proof.json`. Each submitted occurrence is consumed once. No candidate line reuses a submitted occurrence, and removed baseline lines are not used as numerical credit for different new lines. All 108match unique actual submitted occurrences out of 145.

## All 45remaining long-function sites

`function-lineage-independent.tsv` lists all 45 candidate flagged functions with physical/coding lengths, exact start/end, submitted identity/start/end/length, and review basis. Every candidate maps to a unique already-flagged submitted site. Ordinary same-name sites were checked against the original source inventory, not merely total counts. Anonymous methods/callbacks were matched by their enclosing S4 method or plotting helper, never by line-number-suffixed labels alone.

The following renamed successors remain flagged. Their functional/site lineage was reviewed independently, using the submitted source and the earlier independent instance audit before the owner's final allowlist. It is not asserted from the owner's CSV:

| Candidate | Submitted predecessor | Reviewed responsibility |
|---|---|---|
| `analyze_association_pairs` | `run_sceptre_analysis_high_level` | Complete-family orchestration, full draws, permutation/CRT dispatch and assembly |
| `run_permutation_analysis` | `run_perm_test_in_memory` | Response-major permutation execution |
| `analyze_response_block` | `analyze_given_response_ids` | Nested response worker and per-pair native results |
| `run_crt_analysis` | `run_crt_in_memory_v2` | Guide-major CRT coordination and shared-response preparation |
| `analyze_guide_block` | `perform_association_analysis` | Nested guide worker and fitting-scope branch |
| `test_crt_shared_fit` | `crt_glm_factored_out` | Low-level shared-response CRT worker |
| `test_crt_group_fit` | `discovery_ntcells_crt` | Low-level treatment-plus-NT CRT worker |
| `construct_negative_control_pairs` | `construct_negative_control_pairs_v2` | Calibration group/pair sampler coordination |
| `pair_results_from_native` | `construct_data_frame_v2` | Native test-output conversion and optional diagnostic/null columns |
| `integrate_guides` | `apply_grouping_to_result` | Union/singleton/Bonferroni guide-target integration |

Additionally, same-name `set_matrix_accessibility` moves from `R/preprocessing_functions.R` to `R/storage.R`; the detailed S4print method and nested `grnas_to_plot` plotting callback retain their submitted responsibilities despite shifted anonymous labels. The table's exact baseline locations and one-to-one matching rule prevent one old monolith from excusing multiple newly flagged helpers. Functional/site lineage is not a claim of whole-implementation AST or numerical equivalence to the submitted 21revision; those are different qualification questions.

## All 9display operations and removed new categories

`display-lineage-independent.tsv` maps all nine display operations. Current `R/output_resources.R:18` retains the detailed summary print formerly at submitted `R/s4_analysis_functs_2.R:723`. Current S4print cat lines 90,92,112,243,294,305,321,334 map to submitted 91,93,107,234,266,272,285,295. There is no new display call, and the cleanup does not substitute the shorter `show()` output for detailed serialization.

Both independent scanners find zero `1:...` sequence, paste-in-signaler, and superassignment instances in the candidate, as in submitted 21. The output-level membership uses `seq_len(3L)`. Progress formatting is computed before `message()`. Warning collection uses an explicitly local environment while preserving warning propagation and diagnostics. Migration stage helpers return the updated object, with explicit ordered rebinding. No warnings were muffled to remove these findings.

## Cleanup quality, visibility and runnable contracts

The cleanup extracts coherent included helpers: manifest-envelope validation; work-unit execution/fit packaging; collector proposal ordering, cache merging and finalization; legacy component/request/result migration; per-entry fit validation; storage counting/order restoration; and preparation-stage checks. The reviewed implementations remain normally spaced and readable. A parser scan found no added semicolon packing, suppression calls, or source/sys.source loading. The cleanup changes no `src/`, tests, dependency declarations, NAMESPACE or .Rbuildignore relative to 169c16e. The parsed RcppExports wrappers are exactly AST-identical to 169c16e; only layout changes.

Relative to actual submission21, build ignores add development evidence/doc directories, not a replacement executable implementation. The extracted archive includes the production helpers; it does not source ignored tools code at runtime. New digest usage supplies input/execution identity certificates, and the BH/testthat metadata changes serve real compiler/test requirements; no artificial Bioconductor dependency was added to silence the inherited warning.

All three exported distributed API Rd files contain real, unguarded runnable examples with public import, setup, guide assignment, QC, prepare, every work unit, collect, and get_result. They do not use tests-only fixtures, private budget knobs, dontrun/donttest, or comment-only placeholders. Preparation happens once and collection follows proposals in reversed arrival order, making the lifecycle contract concrete. `runnable-example-structure.tsv` records this independent static review. The coordinator separately reports successful execution in a fresh exact-archive library (two units and a current 418-row result for each example); that execution was not rerun by this reviewer.

## Result and limits

No material diagnostic-instance or cleanup-evasion issue remains in this reviewed candidate. The saved independent instance summary records 45 inherited flagged function sites, 108 exact inherited long-line instances, 9 inherited display operations, and zero new source-instance/category findings. Full numerical, collector, installed-example, package-check, and pipeline qualification are separate evidence owned by the coordinator/other reviewers. Original-vs-local runtime/checker/options differences are preserved, and no overall acceptance or broad-platform claim follows from this source-instance pass.
