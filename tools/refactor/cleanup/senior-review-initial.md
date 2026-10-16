# Independent cleanup-plan review

Outcome: CHANGES REQUESTED before affected implementation resumes. The proposed responsibility split is supported; the required amendments below concern safety, precise dispositions and migration dependencies. No P0 found. This is a plan review, not implementation or qualification approval.

Reviewed source pin: c7d62502757f2a17741e2ef18b9e4d006bf063dc. All 53 entries in review-coverage.json match both that Git source and current production files by SHA256. This covers 34 handwritten R modules, 14 C++ implementations, three headers and two generated binding files. No production files were changed, tests run, or unsafe examples executed.

## Required amendments

1. P1: Make S2's calibration bound refer to actual post-QC nonempty NT assignment groups, not metadata guide-label counts. The existing validator counts metadata rows, while construct_negative_control_pairs samples names from the pruned assignment list. QC can remove an entire NT guide, making a previously valid/default calibration_group_size exceed the active count and reach unchecked iterate_over_combinations writes. Require a finite scalar integer group size in [1, length(current NT groups)] before combinations, and validate positive integral calibration-pair requests before narrowing. Include a public post-QC group-loss regression. Do not silently clamp the group size: that changes the calibration treatment definition.

2. P2: Replace S1/S2's generic zero rejection/native-precondition prose with an entrypoint contract table. Fisher-Yates allows 0 <= M <= n_tot and B >= 0; safe zero-size banks/samples are not automatically errors. The hybrid algorithm has stricter positive control/minimum-size prerequisites. Check sums/products in a wide type before narrowing; do not impose an invented arbitrary cardinality cap. For sample_combinations, explicitly choose an informative error for nonfinite/nonpositive p_hat or an unrepresentable estimate, rather than inventing a candidate-rate fallback. Check the skew-normal minimum when the fitter is actually invoked; do not reject unused short B2 metadata from a stage-1-only result. Validate draw-bank construction/bounds at the producer and inexpensive boundary dimensions; avoid a full bank rescan per pair. Preserve S1's M=min(old_upper_bound,n_cells) only for complement calibration, leaving every formerly safe M/draw stream identical.

3. P2: Resolve contradictory context/factory/intercept instructions before CL3. CLEANUP_PLAN I4 retains and renames the useful treatment lookup closure, while core C05 says remove its factory; choose one concrete target and mark the other proposal superseded. More materially, execution E7 deliberately preserves the NT one-column drop/error while B3 proposes correcting it and the target context promises a matrix design. Recommend adopting the bounded drop=FALSE correction, with separate valid intercept-only tests, before or atomically with CL3. Do not first add deliberate vector/error compatibility machinery and then remove it in CL6. Saved full-source design and cache invalidation remain unchanged.

4. P2: Complete the retained/deleted private-symbol crosswalk before its corresponding deletion/rename unit. The slow CRT sampler is still undecided; the actual native R entry name and renamed count helpers also need a clear destination. Recommend retiring the unused slow CRT sampler from production, retaining its immutable historical oracle/fixtures as archive evidence. Use an explicit reviewed retained-sampler mapping for candidate comparisons; never route both actual and expected through the frozen reference or silently filter by what functions happen to exist. Preserve three real result adapters, pipeline parts, determine_grnas_in_use's restricted/trial mode, and pod/file helpers with their pinned actual consumers. Do not add aliases to satisfy tests alone.

5. P2: Give B1-B6 an adopted/deferred/owner-choice ledger, not merely an open list. Proposed dispositions appear below. Each adopted observable correction needs its own rationale, affected domain and independent regression, with frozen normal-case fixtures unchanged. Deferred limitations need an honest completion statement; they must not be declared absent.

6. P2: Make the ownership-copy resource exception explicit. O1's roughly16 bytes/nnz transient copy is an engineering estimate, not measured peak RSS or qualification. Require a meaningful triplet conversion measurement including sorting workspace. If it crosses the existing general resource gate, record and review that bounded correctness cost as an exception or investigate it; do not relabel the old resource report as passing the new path. Existing general cases must still meet their agreed gates.

## Recommended behavior dispositions under current cleanup authority

- S1 and precisely bounded S2: adopt safety fixes. Error on invalid/unrepresentable prerequisites; preserve safe draw streams and do not choose a new sampler. A zero empirical candidate rate gets an informative failure, not an invented probability floor. If successful zero-rate calibration is wanted, that requires Eugene's scientific choice.
- B1: adopt stable cell/feature ordering and actual feature-ID-to-local-row matching, preserving stored entries/explicit zeros and the separately characterized duplicate policy. No automatic duplicate coalescing. Any proposed summation/coalescing policy needs an explicit additional decision.
- B2: adopt removal of the cis early exit and safe coordinate arithmetic. This corrects the stated distance predicate while preserving supplied row order; it selects no new statistical method.
- B3: adopt preserved matrix dimensions for already accepted intercept-only formulas. Test NT calibration/discovery and both traversal modes with an independent view/calculation. Do not infer expected numbers from the new context helper itself.
- B4: adopt documented QC option forwarding; correct bin expansion to run only beyond max_single_bin after a reproducer; make volcano colors follow the existing significant flags rather than define a new tie rule. Per-panel QQ selection may be a bounded correction after a focused reproducer establishes the intended existing one-band-per-panel contract. Defer Bonferroni guide/target pooling or labeling interpretation, and any new global downsampling cap; ask Eugene only if changing those semantics now is required.
- B5: adopt the wrong mixture-range variable correction as a named invalid-input bug fix. Defer broad new validation/error precedence and altered numerical fallback policy. Guarantee connection cleanup without redesigning malformed-file support. If narrowing estimate_theta's catch is proposed, first enumerate actual numerical exception classes and preserve their existing estimate/fallback arithmetic; propagation of nonnumerical failures can then be reviewed as a separate bounded exception-policy fix.
- B6: retain and document the current backend limitation for this refactor. Existing in-memory aggregation may be made sparse only after exact ordering/duplicates/numerical/resource evidence. Adding backed aggregation is a new backend capability and needs Eugene's choice if sought now.
- O1: adopt ownership corrections with independent baseline copies, caller serialization and measured sparse-copy cost; returned mathematics/schema/order remain unchanged.

## Design and qualification assessment

The plan simplifies underlying logic: it eliminates unused modes and fallbacks, makes count coordinates direct, and establishes one transient view/design authority. Keeping separate permutation/CRT traversal, guide-specific NT fits, compact response caches, active mixture fallback and pipeline facade adapters is justified by actual behavior and consumers. No replacement architecture is required.

Coverage is broad and source-pinned; it is not evidence of bug absence. Keep the existing real-consumer, fixed-bank RNG, saved-object, cache-operation-count, serial/fork and release-resource gates. Correctness regressions should compare independent intended calculations for newly repaired domains; ordinary supported domains retain the pinned original oracle. Finish with a final source/deletion review after the final immutable build. This review does not itself authorize new scientific choices or claim pristine source.

## Reviewed artifact hashes

- `tools/refactor/CLEANUP_PLAN.md`: `8a13a89021fba8df8f175fbfcd9205541257d7f5835be7e0251184fe30bd3646`
- `tools/refactor/cleanup/review-coverage.json`: `294d286459c2809b2782296a37a0190d14483d4218c8fd19c6e83f2b0e2f2a76`
- `tools/refactor/cleanup/core-review.md`: `aac2c5d77af64e574280ec0e96fa77c7eb51095ed3dea7a649294506613e38ad`
- `tools/refactor/cleanup/core-inventory.json`: `d336a65d7e98cb5148346f58c97bf6f90aa8d5de68e781c4c50d98bd63411ea0`
- `tools/refactor/cleanup/execution-proposal.md`: `a19037a1910dfc90bc3717fbf7f3261fe4a8c08771bfe5d8e30372a9d1844f32`
- `tools/refactor/cleanup/execution-inventory.json`: `eab7742dd09a803983cdd1a672e6fe0898c6d598f095d8cc0561582cd7695854`
- `tools/refactor/cleanup/legacy-findings.json`: `c8a64b5d8d58ea0cfd547fac085cd158485ac661fb5f6102c45db22db8595a8e`
- `tools/refactor/cleanup/legacy-review.md`: `6773caf877ee2f200d5dc7ee3d96fcbdce62a35f08afda27ea0c73f871e2c366`
- `tools/refactor/cleanup/native-review.md`: `41601307675492521bc50c33e1e488be86b940fd01c3fa416645d9ae8aa39d95`
