# F1/F3/F4 independent source review

Exact candidate `9afed2983522089a0d12450296881b348d7dbc41` (d7fee2d plus the data.table-safe follow-up), base `3d9fafc`. Reviewed all five owned files, actual plotting/result producers, transient-view callers and numerical factor formula. Read-only; no R process/build/test/benchmark.

**Outcome: no P0/P1/material P2 source finding.** Runtime checks on this exact follow-up remain a separate gate; earlier 45 source-evaluated checks covered the preceding pin and are not relabeled as expanded final-runtime evidence.

## F1 diagnostic invariance

The former whole-row na.omit dropped valid discovery rows for amount 2 missing skew-normal parameters and amount 3 padded null draws. The replacement requires completeness only of p_value, log_2_fold_change and significant: the primary values actually consumed by the QQ/volcano/common-family summary. It preserves existing significance flags/threshold construction, plot clipping/downsampling/RNG and panel family semantics; it does not re-test hypotheses or replace significance by a p threshold. Rows missing a primary value remain excluded. No new finite-value restriction silently excludes existing infinite effects.

Optional diagnostic fields do not select rows anymore. The explicit three `$` vectors avoid data.table j-symbol semantics and make the meaning independent of package-awareness/source-evaluation context. Actual union/singleton result producers can retain data.table class, so this matters beyond a fabricated dataframe fixture. Subsetting creates the display view without mutating the facade/retained table alias.

The new independent primary-result fixture checks known p/effect/significance values and discovery count. Amounts 1/2/3 are compared through actual built QQ bands/points, volcano data/threshold and text, rather than only row counts or unevaluated ggplots. The follow-up covers both data.frame and data.table input, actual aliases, serialized facade/result nonmutation and caller RNG; missing-p/effect/flag/failed-QC examples remain separate. The expanded runtime is honestly deferred under the resource quiet window.

## F3 transient field retirement

Removed control_fit_pos was not read by production or pinned consumers; it was reconstructed only for tests. The fit-to-analysis/source/treatment/n_trt view still carries every actually used value. There is no saved-slot/API change or replacement alias. Updated tests independently derive complement positions and assert actual source/design mapping for NT discovery, NT calibration and complement controls, instead of merely dropping the lost assertions. The return contract now accurately says the control complement is implicit rather than materialized.

## F4 score contract correction

The factor is `(Lambda^(-1/2) * t(U)) %*% t(wZ)`. Its rows are eigenbasis coordinates in eigen() order, with count equal to design coefficient count; they are not original coefficient coordinates/order. The amended prose correctly distinguishes those statements. Numeric expressions/eigen order are unchanged. This is a correction to earlier misleading documentation, not a new transform or a claim of exact row identity under arbitrary eigenvector conventions.

## Remaining scope

Run the tiny expanded plot checks and the normal integrated contracts on the coherent installed source. These fixes do not resolve native R-longjmp ownership, B6 precision/resources or final architecture coverage; those remain separately assigned. No scientific method/default or extra display feature is required by this review.

## Source SHA256

- `R/plotting_functions.R`: `530af0bafdd5a752cb0c6b3f410bb5034ac71ce15b7822cdc04e1b93851f0b35`
- `R/analysis_context.R`: `2faabf7efc68b7eea5e4c7c5661c42c93053061234833a30c8f451350b380105`
- `R/score_preparation.R`: `9c945ceca4fe851e9f5192e5da6c169466f031ea1a2cf706660a5d03128949b3`
- `tests/testthat/test-refactor-indexing.R`: `da2b278441b0121213df35bb58be2201d73987af87e94a9597b7c61b14b45dee`
- `tests/testthat/test-refactor-plot-diagnostic-na.R`: `3b6f59e3b1b06791c704e64e987204bb3881157ceacfb43d56cf87dc441c9b39`
