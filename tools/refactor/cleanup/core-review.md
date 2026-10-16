# Core-module audit for the plan-review checkpoint

Audit source: `c7d62502757f2a17741e2ef18b9e4d006bf063dc`. Read-only coverage remains **52 named functions** in nine modules, **28 existing source contracts** and the same source hashes/caller evidence. Later planning dispositions below do not change the audited production behavior. No production edits, tests or resource measurements are performed by this amendment.

The accompanying [core-inventory.json](core-inventory.json) separates pinned `executed_contract` evidence from adopted planning dispositions and explicitly labeled superseded proposals. Held implementation `8f1b8f4` and worker changes `4aac60c`/`13ff940` remain unintegrated at the audit pin.

## Senior-review amendments: final targets

- **C05:** Retain the useful lean treatment lookup closure, renamed `make_treatment_index_lookup`. Its result remains treatment positions/count, with no full per-pair view construction or all-guide cell-map bank. Remove only the identity/fallback routes and their obsolete caller reconstruction. The earlier proposal to remove the factory is superseded.
- **C07/B3:** Adopt the bounded `drop=FALSE` repair for already accepted intercept-only formulas **before or atomically with CL3**. The common design and guide-specific slices remain matrices. Preserve the full source-row saved design and strict cache invalidation comparison. The earlier deferral/vector-error compatibility proposal is superseded; do not implement an intermediate error-preservation route.
- **C02/O1:** An approximately 16-byte/nnz owned numeric payload is an engineering estimate, not measured peak RSS or sorting workspace. Follow [triplet-resource-plan.md](triplet-resource-plan.md); crossing general resource gates requires investigation or an explicit reviewed correctness-cost exception. Earlier resource reports do not qualify this changed path.

These resolve required amendments 3 and 6 of `/private/tmp/sceptre-cleanup-plan-senior-review.md`; the review remains changes-requested for its other plan amendments. This is adopted planning scope, not implementation/qualification approval.

## Material findings and current dispositions

### C01 — Canonical NT coordinates are undone and restored inside owned sparse path

**P1 maintainability/coordinate ownership**. Evidence: `R/analysis_context.R:51-60`; `R/storage.R:97-111`; `src/pairwise_qc.cpp:28-39`.

Use one canonical sparse kernel with NT groups in analysis positions, convert owned native positions to zero once outside loops, remove old private v3 entry point rather than a test-only adapter. Keep external ondisc NT-local adaptation only in its branch.

**Validation:** Noncontiguous/reordered source maps, NT pool not sorted, response-sorted scan with exact request restoration/repeats, stored zeros/fractions, complement/NT totals, input nonmutation; frozen public science/RNG unchanged.

**Resource cost:** O(total position counts) native conversion/copy per scan, replacing redundant R group copies/match and per-response redirection; no dense realization, no extra summary matrix.

### C02 — Triplet conversion reorders caller-owned sparse slots

**P1 ownership contract**. Evidence: `R/storage.R:313-324`.

ADOPTED OWNERSHIP TARGET, PENDING IMPLEMENTATION/QUALIFICATION: use owned sorting vectors while preserving returned values/stored entries/dimnames and stored-FALSE-to-one semantics. The held implementation is not integrated at the audit pin. Qualify the new path using triplet-resource-plan.md.

**Validation:** Unsorted dgT/lgT input, row/column orientation, serialized input unchanged, exact old sceptre output on independent copies; ondisc comparison must restore column names because dependency drops them.

**Resource cost:** Roughly 16 bytes per numeric stored entry describes only an engineering estimate for owned i/j/x payload. It is not measured peak RSS, sorting workspace, or qualification. Measure total/conversion-phase RSS and timing on matched fresh builds; investigate a crossed general gate or obtain an explicit reviewed bounded-correctness exception. Earlier reports do not qualify this new path.

### C03 — Full source-row cached design mislabeled analysis-cell order

**P1 contract accuracy**. Evidence: `R/fit_cache.R:6-10`; `R/s4_analysis_functs_1.R:138-144`; `R/s4_analysis_functs_1.R:175`; `R/association_execution.R:91-97`.

Correct source contract. Preserve full source-row S4 design and strict full-matrix comparison; distinguish it from transient analysis/common/guide-fit design.

**Validation:** Noncontiguous cells and parameter change outside current analysis view; unchanged invalidation decisions and saved schema.

**Resource cost:** Documentation correction zero runtime; context composition should remove duplicate subsetting without persisting additional full designs.

### C04 — Five private aliases have no production package or bounded executable downstream consumer

**P2 obsolete surface**. Evidence: `R/precomputation_functions.R:3-15`; `R/result_processing.R:163-165`; `R/analysis_context.R:154-161`; `tests/testthat/test-refactor-indexing.R:61`; `tests/testthat/test-refactor-indexing.R:85`.

Delete perform_response_precomputation, perform_grna_precomputation, compute_precomputation_pieces, construct_data_frame_v2, get_idx_vector_calibration_check; remove empty file and retarget tests/docs to actual core. Historical proposal prose is not compatibility authority.

**Validation:** Re-run bounded pinned consumer search, exported API/formals/method snapshots, component tests and actual native numerical fixtures; no golden regeneration.

**Resource cost:** Removes wrappers/call layers; no material allocations.

### C05 — Retain the lean treatment lookup closure by purpose; remove identity and fallback routes

**P2 unnecessary indirection**. Evidence: `R/analysis_context.R:113-121`; `R/analysis_context.R:151`; `R/analysis_context.R:165-167`; `R/association_workers.R:143-147`.

ADOPTED PLAN: retain the useful per-run treatment lookup closure, renamed make_treatment_index_lookup, with its lean trt_idxs/n_trt result and existing shared scope. Remove native_observed_positions identity coercion and legacy_fit_source_cells reconstruction fallback with migrated callers. Keep observed one-based R versus synthetic zero-based bank conversion at the real Rcpp boundary. The earlier factory-removal proposal is superseded.

**Validation:** All four worker modes retain treatment order/n_trt, aligned source/design arrays, response/guide fit and bank sharing. Lookup tests exercise its actual membership semantics, including ordered NT calibration, not only a closure identity. Keep native observed-versus-synthetic score fixture; remove passthrough-only assertions.

**Resource cost:** The closure resolves treatment positions only. No full source/view construction per pair, no guide-by-cell view bank, and no additional fitted-vector cache.

### C06 — Assignment-coefficient and reduced-score modes are test-only

**P2 unused branches**. Evidence: `R/association_workers.R:135-139`; `R/association_workers.R:216-220`; `R/fitting.R:54-72`; `R/score_preparation.R:34-58`; `tests/testthat/test-refactor-fitting.R:12-16`; `tests/testthat/test-refactor-fitting.R:31-35`.

After deleting generic aliases, prefer probabilities-only assignment fitter and full current score preparation; remove unused flags/branches unless review identifies a concrete current supported consumer. Do not keep complexity solely for hypothetical extensions or self-authored tests.

**Validation:** Current production fit counts/probabilities/score arrays exact, all scientific fixtures unchanged; direct inference coordinate fixture uses current full score record.

**Resource cost:** No added allocation; do not couple fitting to score/D or introduce an extension registry to replace unused flags.

### C07 — Context is partly duplicated and its options field is unused

**P1 single context authority**. Evidence: `R/analysis_context.R:123-146`; `R/analysis.R:14-56`; `R/association_execution.R:78-97`; `R/association_execution.R:257-275`.

ADOPTED PLAN: compose one transient common fit view and aligned matrix design once, with drop=FALSE. Adopt the bounded B3 intercept-only dimension correction before or atomically with CL3; guide-specific NT slices also remain matrices. Choose one authority for options/index/cache and remove duplicate driver fallback routes. Preserve the full source-row S4 design and strict full-matrix cache comparison. The earlier decision to defer the repair is superseded; no temporary vector/error compatibility route is retained.

**Validation:** Independently construct intended retained-source and treatment/NT rows for NT calibration and NT discovery, exercising accepted intercept-only formulas under permutations and CRT. Assert matrix dimensions and numerical results from the independent view/calculation, not expectations generated by the new context helper. Preserve frozen normal-case fixtures and cold/warm, serial/fork fit/bank counts.

**Resource cost:** Reuse one common matrix view/design; form guide-specific matrices only at their existing reuse scope. No all-guide view bank, persistent design duplication, or per-response fitted-vector cache. The B3 repair is a prerequisite to or part of CL3.

### C08 — Overlapping NT saved-schema test freezes malformed unsupported input

**P2 test authority**. Evidence: `tests/testthat/test-refactor-indexing.R:84-99`; `R/analysis_context.R:20-30`; `R/cellwise_qc_functs.R:98-116`; `R/s4_analysis_functs_1.R:93-99`.

Retain real saved low-MOI NT-local schema and ordered pool; test nonempty disjoint groups after QC. Separately test supported high-MOI union deduplication under complement. Do not require out-of-range cumulative runs for overlapping low-MOI groups.

**Validation:** Saved-object reopening/getter/assignment fixtures plus valid partition maps. If schema normalization is ever changed, treat that as a genuine serialized boundary decision.

**Resource cost:** No feature or matrix-layout change.

### C09 — Cache compactness, miss-only merge and invalidation are real semantics

**retain deliberate invariants**. Evidence: `R/fit_cache.R:27-85`; `R/association_execution.R:129-139`; `R/association_execution.R:213`; `R/association_execution.R:300-325`; `R/s4_analysis_functs_1.R:497-500`.

Keep setter/QC-owned invalidation and response-ID compact records; document merge as append, not upsert. Do not add assignment cache, opaque fingerprints, fitted-vector caches or PR233 behavior.

**Validation:** Cold/warm fit operation counts, parameter/control/low-MOI assignment/analysis-cell invalidations, no invalidation for mere report/QC-threshold changes with same view.

**Resource cost:** Existing compact saved cache only; no full response x cell state.

### C10 — Three result adapters have actual pinned pipeline ::: consumers

**retain genuine consumer boundary**. Evidence: `R/result_processing.R:145-161`; `sceptre-pipeline@38dc760:bin/process_association_analysis_results.R:34-36`; `sceptre-pipeline@38dc760:bin/run_association_analysis_nuclear.R:66`.

Retain process_calibration_result, process_discovery_result, apply_grouping_to_result as thin facade adapters until coordinated pipeline migration. Keep core finalization on plain tables. Add contracts stating complete-family standard correction versus existing nuclear raw shard behavior.

**Validation:** Pinned standard/nuclear artifact tests, Bonferroni integration, passing/failed rows, output amounts 1/2/3, diagnostic alpha collision and strict significance threshold.

**Resource cost:** No new cross-pod combine or global nuclear matrix/manifest.

### C11 — Similar count names do not imply interchangeable definitions

**retain observed numerical/input semantics**. Evidence: `src/sparse_matrix_functs.cpp:40-96`; `src/pairwise_qc.cpp:25-40`; `R/storage.R:265-287`; `tests/testthat/test-refactor-storage.R:16-28`; `tests/testthat/test-refactor-storage.R:55-77`.

Keep loaded-value >0.5 cell count, stored-position pair count, stored-FALSE conversion, sparse fractional values, and existing backed integer count/threshold coercion distinct. No scientific normalization disguised as source cleanup.

**Validation:** Existing stored-zero/fraction/threshold representation cases plus optional ondisc lane with real initialized handles; assertions should reflect supported representation differences.

**Resource cost:** No extra scans or dense realization.

### C12 — 28 source contracts exist but coverage does not establish accuracy; 24 assigned-module functions lack them

**P2 contract completion**. Evidence: `R/storage.R`; `R/fitting.R`; `R/score_preparation.R`; `R/fit_cache.R`; `R/pair_qc.R`; `R/workflow_state.R`; `R/analysis_context.R`; `R/result_processing.R`; `R/precomputation_functions.R`.

Correct C03, label threshold/index validation as caller preconditions, then add concise @noRd contracts only to retained named helpers in context/results. Remove obsolete functions instead of documenting them into permanence.

**Validation:** Per-function caller/shape/order/side-effect audit in inventory.json; parse/generated-doc/API checks supplement but cannot substitute for that review.

**Resource cost:** No runtime cost; avoid ceremonial anonymous-closure documentation.

## Complete named-function coverage

Caller names and source spans in this table are the historical audit evidence; the final column is the current planning target.

| Module/function at audit pin | Actual package callers | Current classification and disposition |
| --- | --- | --- |
| `R/analysis_context.R:4 prepare_cell_map` | prepare_index_context, update_grna_assignments_given_qc | necessary coordinate core: Retain central source-view validation; source map uniqueness differs from duplicate pair request policy. |
| `R/analysis_context.R:14 source_membership_to_analysis` | update_grna_assignments_given_qc | necessary coordinate core: Retain; do not replace ordered NT concatenation with this set-membership rule. |
| `R/analysis_context.R:20 legacy_nt_fit_assignments` | update_indiv_grna_assignments_for_nt_cells | genuine saved-state boundary: Retain saved low-MOI schema and order; rename boundary clearly if useful; replace malformed-overlap characterization with supported partition test. |
| `R/analysis_context.R:32 calibration_group_membership` | get_idx_vector_calibration_check, group_treatment_analysis | necessary coordinate core: Retain common actual union logic; test supported low-MOI partition and high-MOI overlap independently. |
| `R/analysis_context.R:39 resolve_pair_rows` | prepare_index_context | necessary request identity core: Retain identity map; invalid missing labels are precondition violations, not silently supported behavior. |
| `R/analysis_context.R:47 prepare_index_context` | compute_pairwise_qc_information, prepare_analysis_context, run_crt_in_memory_v2, run_perm_test_in_memory | necessary normalization boundary: Retain normalization at saved-facade boundary; avoid converting back for owned sparse kernel. |
| `R/analysis_context.R:79 group_treatment_analysis` | make_fit_view, prepared_index_factory | necessary coordinate core: Retain shared membership rule; keep distinct from fit-view positions. |
| `R/analysis_context.R:87 make_fit_view` | discovery_ntcells_crt, discovery_ntcells_perm_test, run_crt_in_memory_v2, run_perm_test_in_memory | necessary fit-view core: Retain explicit source/analysis/fit maps and order/duplicate semantics at existing shared scopes. The context selects a common matrix design and guide-specific workers use drop=FALSE, adopting B3 before or atomically with CL3; never construct full per-pair or all-guide view banks. |
| `R/analysis_context.R:113 prepared_index_factory` | run_crt_in_memory_v2, run_perm_test_in_memory | retain lean lookup; rename by purpose: ADOPTED PLAN: rename to make_treatment_index_lookup and retain its useful lean per-run closure returning trt_idxs/n_trt. Preserve membership/order semantics and shared scope; do not materialize full per-pair or all-guide cell views. Earlier factory-removal proposal is superseded. |
| `R/analysis_context.R:123 prepare_analysis_context` | run_sceptre_analysis_high_level | necessary boundary; incomplete composition: ADOPTED PLAN: prepare one common fit view and aligned matrix design with drop=FALSE; remove duplicate option/index/cache authorities. Adopt accepted intercept-only formula repair before or atomically with CL3, including guide-specific matrix slices. Preserve the saved full source-row matrix/cache comparison; no temporary vector/error compatibility machinery. |
| `R/analysis_context.R:151 native_observed_positions` | crt_glm_factored_out, discovery_ntcells_crt, discovery_ntcells_perm_test, perm_test_glm_factored_out | remove private passthrough: Remove helper and identity-wrapper tests; keep explicit coordinate contract and meaningful observed-vs-synthetic native score fixture. |
| `R/analysis_context.R:154 get_idx_vector_calibration_check` | none | remove private alias: Delete; test calibration_group_membership directly and remove obsolete private-wrapper claims. |
| `R/analysis_context.R:165 legacy_fit_source_cells` | crt_glm_factored_out | remove unsupported private fallback: Remove fallback and obsolete optional argument route with callers; if a real downstream worker consumer is discovered, localize named adapter there. |
| `R/fit_cache.R:13 response_cache_parameters_changed` | set_analysis_parameters | contract correction required: Correct text to full source-row schema; do not alter saved slot or replace comparison with view hashes. |
| `R/fit_cache.R:27 response_cache_assignments_changed` | assign_grnas | necessary cache policy: Retain supported-setter ownership and explain why control/view combinations differ. |
| `R/fit_cache.R:42 invalidate_response_cache` | assign_grnas, run_qc_pt_1, set_analysis_parameters | necessary cache policy: Retain central mutation point; event documents ownership, not generic extensible event framework. |
| `R/fit_cache.R:57 lookup_response_fit` | crt_glm_factored_out, run_crt_in_memory_v2, run_perm_test_in_memory | necessary compact cache interface: Retain alongside explicit invalidation guarantees, not implied self-validating cache. |
| `R/fit_cache.R:68 merge_response_fits` | run_crt_in_memory_v2, run_perm_test_in_memory | necessary compact cache interface: Retain simple append under actual driver invariant; test warm-cache/repartition behavior, not hypothetical key framework. |
| `R/fit_cache.R:82 commit_response_fits` | run_calibration_check_pt_2, run_discovery_analysis, run_power_check | necessary facade commit: Retain post-success facade boundary and compact saved schema. |
| `R/fitting.R:16 fit_response_current` | discovery_ntcells_crt, discovery_ntcells_perm_test, perform_response_precomputation, run_crt_in_memory_v2, run_perm_test_in_memory | necessary numerical core: Retain arithmetic and named failure contract; input design is already selected by caller, not full-source slot. |
| `R/fitting.R:54 fit_assignment_current` | crt_glm_factored_out, discovery_ntcells_crt, perform_grna_precomputation | test-only alternate branch decision: Prefer probabilities-only core after deleting unused aliases; otherwise name a real retained consumer for the flag. |
| `R/pair_qc.R:11 collect_pair_qc_metrics` | compute_pairwise_qc_information | necessary measurement boundary: Retain semantic measurement stage required by accepted QC architecture; do not call it an external compatibility shim. |
| `R/pair_qc.R:29 evaluate_pair_eligibility` | compute_qc_metrics | necessary policy boundary: Retain one rule owner; clarify parameter preconditions instead of claiming comprehensive scalar validation. |
| `R/precomputation_functions.R:3 perform_response_precomputation` | none | remove private alias: Delete file alias and update developer guide to actual fitter; do not preserve for historical baseline naming. |
| `R/precomputation_functions.R:7 perform_grna_precomputation` | none | remove private alias: Delete with obsolete return-mode decision; keep actual numerical core. |
| `R/precomputation_functions.R:11 compute_precomputation_pieces` | none | remove private alias: Delete; retarget native coordinate fixture to prepare_score_current and actual compact record. |
| `R/result_processing.R:3 pair_results_from_native` | construct_data_frame_v2, run_crt_in_memory_v2, run_perm_test_in_memory | necessary result core: Retain; document row/schema/type contracts and actual data.frame/data.table behavior, test no accidental scientific column changes. |
| `R/result_processing.R:64 integrate_guides_current` | apply_grouping_to_result, finalize_association_result | necessary result core: Retain supported strategy semantics; document mapping cardinality/preconditions and no extra correction for calibration. |
| `R/result_processing.R:123 assemble_pair_results` | run_sceptre_analysis_high_level | necessary family boundary: Retain family assembly separate from inference; test failed rows and diagnostic schema. |
| `R/result_processing.R:131 finalize_association_result` | process_calibration_result, process_discovery_result | necessary family boundary: Retain one complete-family standard finalizer; full outputs/Bonferroni/alpha-collision tests are material. |
| `R/result_processing.R:145 process_calibration_result` | run_calibration_check_pt_2; pinned pipeline calls verified | genuine executable consumer boundary: Retain thin adapter until coordinated pipeline migration; genuine caller, not baseline/test justification. |
| `R/result_processing.R:151 process_discovery_result` | run_discovery_analysis; pinned pipeline calls verified | genuine executable consumer boundary: Retain adapter and both script contracts; do not introduce global nuclear correction. |
| `R/result_processing.R:157 apply_grouping_to_result` | run_power_check; pinned pipeline calls verified | genuine executable consumer boundary: Retain until coordinated external migration; source-contract needed, no generic alias cleanup. |
| `R/result_processing.R:163 construct_data_frame_v2` | none | remove private alias: Delete; no need to keep renamed private test-only aliases or old version suffix. |
| `R/score_preparation.R:13 compute_D_matrix` | prepare_score_current | necessary numerical core: Retain exact algebra/eigen ordering; independently test shape and finite reference calculation. |
| `R/score_preparation.R:34 prepare_score_current` | compute_precomputation_pieces, crt_glm_factored_out, discovery_ntcells_crt, discovery_ntcells_perm_test, run_perm_test_in_memory | test-only alternate branch decision: Prefer pruning dead reduced branch/flag without coupling fit to D; defer new inference-specific preparation until an actual method. |
| `R/storage.R:12 storage_response_index` | assign_grnas_to_cells_thresholding, load_row | necessary adapter: Retain at the storage/ID boundary; document preconditions, not invented validation. |
| `R/storage.R:26 storage_read_response` | crt_glm_factored_out, discovery_ntcells_crt, load_row, run_crt_in_memory_v2, run_perm_test_in_memory | necessary adapter: Retain; no eager full assay or claim of arbitrary block-read efficiency. |
| `R/storage.R:46 storage_threshold_row` | assign_grnas_to_cells_thresholding | necessary adapter: Retain backend boundary and representation characterization. |
| `R/storage.R:65 storage_owned_groups` | storage_group_counts, storage_pair_counts | necessary defensive adapter: Keep only where the actually called kernel mutates R vectors; remove redundant sparse-pair copies in canonical proposal. |
| `R/storage.R:94 storage_pair_counts` | collect_pair_qc_metrics | simplify internal coordinates: Canonical sparse kernel consumes analysis groups directly, converts positions once without R mutation; keep NT-local conversion only for external ondisc. |
| `R/storage.R:162 storage_group_counts` | construct_negative_control_pairs_v2 | necessary adapter: Retain batched behavior and explicit allocation contract; no per-cell R loops. |
| `R/storage.R:192 storage_create_backed_pair` | write_matrices_to_disk | genuine compatibility boundary: Keep private adapter: public ondisc constructor fixes ID zero and cannot preserve existing import identity. |
| `R/storage.R:219 storage_reopen_backed_pair` | read_ondisc_backed_sceptre_object | genuine compatibility boundary: Retain for public saved-object reader and existing .odm layout, no serialized portable pointers. |
| `R/storage.R:244 storage_cell_covariates` | compute_cell_covariates | necessary adapter: Retain native boundary; do not unify its positive-value count with stored-position pair counts. |
| `R/storage.R:270 set_matrix_accessibility` | assign_grnas_to_cells_mixture, assign_grnas_to_cells_thresholding, compute_cell_covariates, plot_assign_grnas, plot_grna_count_distributions, process_covariates_and_matrices | ownership correction adopted; qualification pending: ADOPTED OWNERSHIP TARGET, PENDING QUALIFICATION: own sorting vectors and preserve exact returned semantics. Execute triplet-resource-plan.md after approval; 16 bytes/nnz is a payload estimate, not measured total RSS/workspace or an automatic resource exception. |
| `R/storage.R:354 convert_pointer_to_index_vector_v2` | set_matrix_accessibility | rename private helper: Rename sparse_pointer_to_indices; only two callers in storage.R, no bounded external executable consumer. |
| `R/workflow_state.R:9 workflow_steps` | workflow_initial_state, workflow_status, workflow_transition | genuine saved-state boundary: Retain fixed saved names/order and sibling semantics. |
| `R/workflow_state.R:20 workflow_initial_state` | init_sceptre_object | necessary state policy: Retain one initialization owner; public import behavior unchanged. |
| `R/workflow_state.R:36 workflow_transition` | perform_status_check_and_update | necessary state policy: Retain one owner; test reruns and siblings, do not add workflow framework. |
| `R/workflow_state.R:60 workflow_status` | print(sceptre_object) method in s4_helpers.R:87 | genuine saved-state display boundary: Retain unknown fields for display of supported historical objects; correct anonymous print-method caller attribution. |
| `R/workflow_state.R:76 perform_status_check_and_update` | assign_grnas, run_calibration_check_pt_1, run_discovery_analysis, run_power_check, run_qc_pt_1, set_analysis_parameters | necessary facade boundary: Retain supported wrapper boundary; no state transition injected in numerical core. |

## Implementation sequence after overall plan acceptance

1. Record senior plan-review amendments and precise genuine boundaries before implementation resumes.
2. Accept contract accuracy and final retained/deleted/renamed symbol crosswalk.
3. One atomic canonical-kernel+adapter+generated-wrapper+focused-test change; no obsolete test-only v3 shim.
4. Separate ownership correction with independent output parity/nonmutation and the explicit triplet resource protocol; crossed gates require investigation or reviewed exception.
5. Before or atomically with CL3 adopt bounded drop=FALSE B3 repair, then context/worker migration with retained make_treatment_index_lookup closure and original fit/bank sharing.
6. Complete source contracts for retained helpers and actual consumer/package/API/scientific checks.
7. Requalify affected resource paths on matched fresh release builds; do not reuse old reports as qualification of new ownership path.

Meaningful workloads use warm-up plus at least three measured repetitions, fixed inputs/dependencies/BLAS/threads/workers/budgets/cache. Existing general gates remain median elapsed <=+10% and peak worker RSS increase <=max(10%,64 MiB). The dedicated numeric/logical unsorted triplet experiment in triplet-resource-plan.md measures total and conversion-phase RSS including sorting workspace. Preserve initial crossings/outliers; investigate or record an explicitly reviewed ownership-correctness exception without relabeling the result a general-gate pass. Other general cases must still pass their agreed gates.

**Preserve:** Public exports/formals/methods/56-slot schema, supported statistical fits/QC/RNG/draw order/guide integration/family correction, ondisc optionality/.odm layout, backend support or pipeline standard/nuclear artifact contracts.

**Named bounded corrections:** B3 repairs the previously failing intercept-only accepted-formula domain; O1 removes accidental caller mutation with unchanged returned semantics. They require independent regressions, frozen normal-case fixtures and the stated resource review.

## Retrieval limits and compatibility authority

Does not search arbitrary unpublished user scripts, all prior chats, remote repositories or execute complete regression/resource suite. Test string mentions and historical documentation are not executable consumer authority. Call graph is lexical, not dynamic proof.

External executable search covered local sceptre-pipeline, sceptre-book, sceptre-team and sceptre3-project-v2 at the recorded revisions where available; sceptre-team has no committed HEAD. All four pipeline references matched pinned blobs. Historical prose and test-only calls are not authority for retaining obsolete internals. Existing supported saved-object/NT-local schemas, actual pipeline result/part/file consumers, native pointer/base boundary, and optional ondisc file interfaces remain genuine boundaries. Scientific semantics and fit/cache/draw reuse invariants remain requirements independently of old private symbol names.

No new tests or timing claims were produced during this planning amendment. Earlier ownership-regression results belong to the held proposal, not the audited production tree.
