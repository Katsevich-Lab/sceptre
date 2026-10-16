# Remaining-source audit checklist

Pinned source: `c63526329f386bb01ddb6fe7fc0a31d2cbd23077`. Read-only inventory; production implementation is paused. Paused documentation checkpoint: `a2e66f0`.

21 included R files plus three testthat helper files. Fourteen R modules excluded by explicit ownership. Static caller candidates cover every R module, every testthat R file, and pipeline bin R files at 38dc760a5b4105d1c3f339cd45f2c98d971087f4. Dynamic dispatch/other external private consumers require consumer-inventory/manual checks; absence here is not proof of global non-use.

## Coverage

| File | Named top-level functions | Public | Existing private @noRd | Checkpoint private drafts | Review |
|---|---:|---:|---:|---:|---|
| `R/assign_grna_functions.R` | 5 | 0 | 0 | 5 | source read; generated binding inventory only where marked |
| `R/aux_functions.R` | 4 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/cellwise_qc_functs.R` | 3 | 0 | 0 | 3 | source read; generated binding inventory only where marked |
| `R/check_functions.R` | 7 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/data.R` | 0 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/getters_and_setters.R` | 6 | 4 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/import_functs.R` | 17 | 3 | 0 | 14 | source read; generated binding inventory only where marked |
| `R/neg_control_functions.R` | 3 | 0 | 0 | 3 | source read; generated binding inventory only where marked |
| `R/ondisc_functs.R` | 4 | 2 | 0 | 2 | source read; generated binding inventory only where marked |
| `R/pair_constructor_functs.R` | 3 | 3 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/pairwise_qc_functs.R` | 2 | 0 | 0 | 2 | source read; generated binding inventory only where marked |
| `R/plotting_functions.R` | 14 | 8 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/preprocessing_functions.R` | 4 | 0 | 0 | 4 | source read; generated binding inventory only where marked |
| `R/qq_plot_helpers.R` | 3 | 0 | 3 | 0 | source read; generated binding inventory only where marked |
| `R/RcppExports.R` | 21 | 0 | 2 | 0 | source read; generated binding inventory only where marked |
| `R/s4_analysis_functs_1.R` | 5 | 3 | 0 | 2 | source read; generated binding inventory only where marked |
| `R/s4_analysis_functs_2.R` | 8 | 5 | 0 | 3 | source read; generated binding inventory only where marked |
| `R/s4_class.R` | 0 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/s4_helpers.R` | 1 | 0 | 0 | 1 | source read; generated binding inventory only where marked |
| `R/sceptre.R` | 0 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `R/test_functs.R` | 2 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `tests/testthat/helper-mock-data.R` | 8 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `tests/testthat/helper-refactor-native.R` | 5 | 0 | 0 | 0 | source read; generated binding inventory only where marked |
| `tests/testthat/helper-refactor.R` | 25 | 0 | 0 | 0 | source read; generated binding inventory only where marked |

{
  "included_R_files": 21,
  "included_test_helpers": 3,
  "excluded_R_files": 14,
  "handwritten_top_level_functions": 91,
  "exported_top_level_functions": 28,
  "internal_top_level_functions": 63,
  "baseline_internal_noRd": 3,
  "checkpoint_internal_contracts": 39,
  "remaining_missing_internal_contracts": 21
}

The three S4 show/print/plot methods, three QQ ggproto methods, named local helpers, anonymous callbacks and three test factory bindings are separately included in JSON. Anonymous/local callbacks inherit an enclosing contract. No boilerplate is proposed.

## Explicit ownership exclusions

- `R/storage.R`
- `R/fitting.R`
- `R/score_preparation.R`
- `R/fit_cache.R`
- `R/pair_qc.R`
- `R/workflow_state.R`
- `R/execution.R`
- `R/association_execution.R`
- `R/association_workers.R`
- `R/mixture_model_functs.R`
- `R/analysis.R`
- `R/analysis_context.R`
- `R/result_processing.R`
- `R/precomputation_functions.R`

## R/assign_grna_functions.R

Assignment/pooling stage; source-cell memberships until QC. Pipeline pods defer pooling. All 5 internal contracts drafted in paused checkpoint.

- `assign_grnas_to_cells` (line 4): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:327.
- `assign_grnas_to_cells_thresholding` (line 84): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/assign_grna_functions.R:39.
- `assign_grnas_to_cells_maximum` (line 104): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/assign_grna_functions.R:46.
- `process_initial_assignment_list` (line 136): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/assign_grna_functions.R:75; /Users/ekatsevi/code/research/sceptre-pipeline/bin/process_grna_assignments.R:40. Manifest consumers: sceptre-pipeline/bin/process_grna_assignments.R.
- `determine_grnas_in_use` (line 203): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/assign_grna_functions.R:21; tests/testthat/test-aux_functions.R:175; /Users/ekatsevi/code/research/sceptre-pipeline/bin/output_grna_info.R:16. Manifest consumers: sceptre-pipeline/bin/output_grna_info.R.

## R/aux_functions.R

Default formula selection, aligned import summaries, deterministic worker partitions, and log-directory ownership. Four undocumented private functions.

- `auto_construct_formula_object` (line 1): missing named internal contract. Call candidates: R/s4_analysis_functs_1.R:106; R/s4_analysis_functs_1.R:290; tests/testthat/test-aux_functions.R:202; tests/testthat/test-aux_functions.R:206.
  Contract: Data frame -> formula. Excludes mito, optionally guide totals; logs count-named columns, omits other columns with >=15 distinct values. Column order determines formula order. Empty eligible sets/non-syntactic names need an explicit supported-input policy.
- `auto_compute_cell_covariates` (line 50): missing named internal contract. Call candidates: R/import_functs.R:964; tests/testthat/test-aux_functions.R:35; tests/testthat/test-aux_functions.R:60; tests/testthat/test-aux_functions.R:78; +3 more in JSON.
  Contract: Aligned response/guide matrices and extra covariates -> source-cell-ordered data frame with prefixed response/guide summaries and response column IDs as row names. Column scans allocate converted sparse storage; no backed-store creation.
- `partition_response_ids` (line 87): missing named internal contract. Call candidates: R/association_execution.R:195; R/association_execution.R:317; R/association_execution.R:406; R/mixture_model_functs.R:109; +1 more in JSON.
  Contract: Character response IDs + parallel/core options -> named list of ID partitions. Uses fixed-seed permutation only above the 2*n_processors threshold, otherwise preserves input order. Reads/restores parallelly option; preserves caller RNG with withr.
- `get_log_dir` (line 136): missing named internal contract. Call candidates: R/execution.R:7; R/execution.R:17; R/mixture_model_functs.R:57; R/mixture_model_functs.R:103.
  Contract: Character base path -> base/sceptre_logs/ path; creates directories recursively. Actual filesystem side effect, not merely formatting.

## R/cellwise_qc_functs.R

Raw source indices -> retained analysis positions -> ordered NT-local compatibility groups. Repeated QC must rebuild from raw assignments. Three drafted contracts.

- `determine_cells_to_retain` (line 1): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:466.
- `determine_cells_to_retain::exclude_cells_by_clipping_counts` (line 13): missing named internal contract. Call candidates: R/cellwise_qc_functs.R:17; R/cellwise_qc_functs.R:22.
- `update_grna_assignments_given_qc` (line 78): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:503.
- `update_indiv_grna_assignments_for_nt_cells` (line 124): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/cellwise_qc_functs.R:111.

## R/check_functions.R

Seven validation entrypoints. Error text/order and accepted input sets are behavior; documented preconditions must not be confused with new validation. Repeated processor guards and unused formal are cleanup candidates.

- `check_import_data_inputs` (line 1): missing named internal contract. Call candidates: R/import_functs.R:114; R/import_functs.R:162; R/import_functs.R:624; tests/testthat/test-check_functions.R:38; +35 more in JSON.
  Contract: Aligned feature-by-cell stores, guide metadata, MOI and optional cell covariates -> NULL/error. Validates supported classes, dimensions, duplicate feature IDs and supplied barcode alignment, metadata and covariates. Does not reorder/coerce counts or validate all matrix values.
- `check_import_data_inputs::check_matrix_class` (line 54): missing named internal contract. Call candidates: R/check_functions.R:75; R/check_functions.R:80.
- `check_set_analysis_parameters` (line 227): missing named internal contract. Call candidates: R/s4_analysis_functs_1.R:121.
  Contract: S4 facade plus pairs/formula/view/test options -> NULL/error. Checks pair membership, duplicate complete rows, formula variables, control compatibility, strategy and backed-permutations restriction. Error ordering is observable; missing required pair columns are accessed before the intended column-name guard.
- `check_parallel_supported` (line 410): missing named internal contract. Call candidates: R/s4_analysis_functs_1.R:265; R/s4_analysis_functs_2.R:70; R/s4_analysis_functs_2.R:246; R/s4_analysis_functs_2.R:349; +6 more in JSON.
  Contract: Logical scalar parallel and OS string -> NULL/error; explicit Windows fork restriction. os_type parameter supports deterministic platform tests.
- `check_assign_grna_inputs` (line 425): missing named internal contract. Call candidates: R/s4_analysis_functs_1.R:308.
  Contract: Object, method, method-specific named hyperparameters and processor count -> NULL/error. Validates option membership/ranges; current perturbation-range upper bound mistakenly reads pi_guess_range. Do not silently tighten accepted values in a refactor.
- `check_run_qc_inputs` (line 565): missing named internal contract. Call candidates: R/s4_analysis_functs_1.R:442.
  Contract: QC cutoffs/percentile vectors and initial guide list -> NULL/error. Requires at least one assignment. Does not validate additional_cells_to_remove or p_mito_threshold; scalar/NA/integer checks are incomplete.
- `check_calibration_check_inputs` (line 611): missing named internal contract. Call candidates: R/s4_analysis_functs_2.R:143.
  Contract: Post-QC facade, pair count, processors -> NULL/error. Uses NT metadata-row count, not surviving assignment-group count. Rejects zero pair count but does not validate group size or all invalid count forms.
- `check_discovery_analysis_inputs` (line 652): missing named internal contract. Call candidates: R/s4_analysis_functs_2.R:266; R/s4_analysis_functs_2.R:369.
  Contract: Prepared pair table, controls, metadata, results/counts/processors -> NULL/error; missing calibration result emits a message. calibration_check_run is an unused lazy formal, omitted by current callers.

## R/data.R

Six bundled dataset declarations/documentation blocks; no functions. Preserve dataset object schemas/IDs and source attribution. No transformation implementation here.


## R/getters_and_setters.R

Four exported getters and two private S4-list-payload setters. Saved shape and raw-versus-QC assignment views are material; setter docs are missing.

- `get_response_matrix` (line 37): existing public roxygen; accuracy findings tracked separately. Call candidates: R/analysis_context.R:133; R/cellwise_qc_functs.R:57; R/cellwise_qc_functs.R:82; R/check_functions.R:238; +18 more in JSON. Manifest consumers: sceptre-pipeline/bin/prepare_association_analysis_nuclear.R.
  Contract: S4 -> first response_matrix list payload, retaining backed handle or matrix without eager read. Stripped serialized objects require reopening before this is meaningful.
- `get_grna_matrix` (line 45): existing public roxygen; accuracy findings tracked separately. Call candidates: R/assign_grna_functions.R:12; R/assign_grna_functions.R:141; R/check_functions.R:239; R/getters_and_setters.R:124; +8 more in JSON.
  Contract: S4 -> first grna_matrix list payload with source row/column order unchanged.
- `get_cell_covariates` (line 53): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-import_functs.R:59; tests/testthat/test-import_functs.R:89; tests/testthat/test-import_functs.R:161; tests/testthat/test-import_functs.R:185; +1 more in JSON.
  Contract: S4 -> stored source-cell-ordered covariate data frame; no QC subsetting.
- `get_grna_assignments` (line 106): existing public roxygen; accuracy findings tracked separately. Call candidates: R/s4_analysis_functs_2.R:619; /Users/ekatsevi/code/research/sceptre-pipeline/bin/process_grna_assignments.R:42. Manifest consumers: sceptre-pipeline/bin/process_grna_assignments.R.
  Contract: Assigned S4 -> sparse logical guide-by-source-cell membership matrix; optional QC selects cells_in_use. Native decrement acts on unlist scratch storage; row order follows initial guide list; must preserve raw assignments and repeated-call behavior.
- `set_response_matrix` (line 136): missing named internal contract. Call candidates: R/import_functs.R:831; R/ondisc_functs.R:76; R/ondisc_functs.R:89.
  Contract: S4 + store -> facade with first response_matrix list element replaced. One-element list is a saved-object compatibility boundary, not an arbitrary nested data model.
- `set_grna_matrix` (line 141): missing named internal contract. Call candidates: R/import_functs.R:832; R/ondisc_functs.R:77; R/ondisc_functs.R:88.
  Contract: S4 + store -> facade with first grna_matrix list element replaced; same saved-object convention.

## R/import_functs.R

Three public import facades, fourteen private helpers, and local parsing closures. Caller order, source identities, vector aggregation, backed RNG/files, scratch CSC ownership and canonical paths require distinct contracts. Fourteen drafted checkpoint contracts.

- `import_data` (line 70): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/helper-refactor.R:357; tests/testthat/test-assign_grnas.R:60; tests/testthat/test-assign_grnas.R:141; tests/testthat/test-assign_grnas.R:221; +33 more in JSON.
- `import_data_memory` (line 105): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:93; R/import_functs.R:578; R/import_functs.R:798.
- `import_data_use_ondisc` (line 152): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:82.
- `import_data_from_cellranger` (line 279): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-import_functs.R:44; tests/testthat/test-import_functs.R:74; tests/testthat/test-import_functs.R:101; tests/testthat/test-import_functs.R:119; +6 more in JSON.
- `import_data_from_cellranger_memory` (line 307): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:289.
- `import_data_from_cellranger_use_ondisc` (line 593): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:296.
- `import_data_from_parse` (line 720): existing public roxygen; accuracy findings tracked separately. Call candidates: none statically located.
- `import_data_from_parse::load_matrix` (line 750): missing named internal contract. Call candidates: R/import_functs.R:775; R/import_functs.R:776.
- `init_sceptre_object` (line 813): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:142; R/import_functs.R:194; R/import_functs.R:646.
- `create_csc_matrix` (line 851): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:511; R/import_functs.R:768.
- `get_mtx_metadata` (line 862): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:450; R/import_functs.R:751.
- `write_matrices_to_disk` (line 881): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:186.
- `prepare_directory_to_write` (line 894): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:171; R/import_functs.R:601; R/import_functs.R:887; R/ondisc_functs.R:91; +4 more in JSON.
- `prepare_ondisc_file_path` (line 919): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/ondisc_functs.R:60; R/ondisc_functs.R:64; R/ondisc_functs.R:68.
- `check_path_has_no_tilde` (line 940): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:914; R/import_functs.R:935.
- `process_covariates_and_matrices` (line 958): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:135; R/import_functs.R:174.
- `collapse_grna_target_data_frame` (line 990): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:129; R/import_functs.R:618.
- `collapse_grna_matrix` (line 1008): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/import_functs.R:125.

## R/neg_control_functions.R

Three private functions, prepared NT summaries and native deterministic sampling. Group IDs encode guide order with &; requested-count shortfalls are messages. Three drafted contracts.

- `construct_negative_control_pairs_v2` (line 1): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_2.R:151.
- `compute_calibration_group_size` (line 128): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_2.R:122.
- `get_undercover_group_names` (line 143): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/neg_control_functions.R:99.

## R/ondisc_functs.R

Two exported saved-object read/write functions plus two private pipeline text/pod helpers. Object payload is stripped/restored around unchanged backing identities. Two drafted private contracts.

- `read_ondisc_backed_sceptre_object` (line 54): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-refactor-contracts.R:345; tests/testthat/test-refactor-contracts.R:358; tests/testthat/test-refactor-contracts.R:360; /Users/ekatsevi/code/research/sceptre-pipeline/bin/assign_grnas.R:19; +10 more in JSON. Manifest consumers: sceptre-pipeline/bin/assign_grnas.R, sceptre-pipeline/bin/output_grna_info.R, sceptre-pipeline/bin/prepare_association_analyses.R, sceptre-pipeline/bin/prepare_association_analysis_nuclear.R, sceptre-pipeline/bin/process_association_analysis_results.R, sceptre-pipeline/bin/process_grna_assignments.R, sceptre-pipeline/bin/run_association_analysis.R, sceptre-pipeline/bin/run_association_analysis_nuclear.R, sceptre-pipeline/bin/run_qc.R, sceptre-pipeline/bin/run_qc_nuclear.R, sceptre-pipeline/bin/set_analysis_parameters.R.
- `write_ondisc_backed_sceptre_object` (line 84): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-refactor-contracts.R:342.
- `get_id_vect` (line 99): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: /Users/ekatsevi/code/research/sceptre-pipeline/bin/output_grna_info.R:18; /Users/ekatsevi/code/research/sceptre-pipeline/bin/prepare_association_analyses.R:39. Manifest consumers: sceptre-pipeline/bin/output_grna_info.R, sceptre-pipeline/bin/prepare_association_analyses.R.
- `write_vector` (line 109): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: /Users/ekatsevi/code/research/sceptre-pipeline/bin/output_grna_info.R:22; /Users/ekatsevi/code/research/sceptre-pipeline/bin/output_grna_info.R:24; /Users/ekatsevi/code/research/sceptre-pipeline/bin/prepare_association_analyses.R:57; /Users/ekatsevi/code/research/sceptre-pipeline/bin/prepare_association_analyses.R:58; +1 more in JSON. Manifest consumers: sceptre-pipeline/bin/output_grna_info.R, sceptre-pipeline/bin/prepare_association_analyses.R, sceptre-pipeline/bin/prepare_association_analysis_nuclear.R.

## R/pair_constructor_functs.R

Three exported pair constructors; no private top-level functions. Ordering, data.table class/ownership and whole-target PC exclusion must stay visible.

- `construct_cis_pairs` (line 58): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-pair_constructor_functs.R:16; tests/testthat/test-parallel_mac.R:21; tests/testthat/test-parallel_mac.R:72.
  Contract: Object + coordinate table + positive controls/distance -> data.frame(grna_target,response_id). Excludes all positive-control targets and NTs, walks chromosomes/targets in observed order, computes midpoint spans then native distance matches. Loads bundled GRCh38 data only into local environment by default.
- `construct_trans_pairs` (line 238): existing public roxygen; accuracy findings tracked separately. Call candidates: none statically located.
  Contract: Object + exclusions -> data.table target/response Cartesian product in expand.grid order. Retains observed target order, response store row order; data.table::setDT currently mutates caller positive_control_pairs class by reference when excluding PC pairs.
- `construct_positive_control_pairs` (line 315): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-pair_constructor_functs.R:13; tests/testthat/test-parallel_mac.R:20; tests/testthat/test-parallel_mac.R:71.
  Contract: Object -> data.frame target/response name matches, unique in guide metadata order; does not fit models or change the facade.

## R/pairwise_qc_functs.R

Two private S4 adapters for measurement/policy; stored sparse positions differ from value>0.5 covariate nonzeros. Response-sorted repeated rows and NT summary shape preserved. Two drafted contracts.

- `compute_pairwise_qc_information` (line 1): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:512.
- `compute_qc_metrics` (line 61): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:516.

## R/plotting_functions.R

Eight public plots and six undocumented private plot helpers, plus two named local bin/histogram closures. All produce presentation objects; forwarding, grouping, RNG, threshold/tie and ID-view behavior must not change accidentally.

- `get_my_theme` (line 1): missing named internal contract. Call candidates: R/plotting_functions.R:184; R/plotting_functions.R:365; R/plotting_functions.R:392; R/plotting_functions.R:441; +9 more in JSON.
  Contract: Numeric text size -> ggplot2 theme; trivial helper needs only a brief shared plotting contract.
- `plot_grna_count_distributions` (line 65): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-plotting_functions.R:97; tests/testthat/test-plotting_functions.R:104; /Users/ekatsevi/code/research/sceptre-pipeline/bin/process_grna_assignments.R:46. Manifest consumers: sceptre-pipeline/bin/process_grna_assignments.R.
  Contract: Object plus selected guides/bin threshold -> list of ggplots or combined cowplot. Reads source-cell guide rows, uses fixed-seed guide selection, clips values >=10000. Binning branch is incorrect when explicit threshold exceeds observed counts above 10.
- `plot_grna_count_distributions::grna_expressions_to_binned_factor` (line 113): missing named internal contract. Call candidates: R/plotting_functions.R:160.
- `plot_assign_grnas` (line 280): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-plotting_functions.R:119; /Users/ekatsevi/code/research/sceptre-pipeline/bin/process_grna_assignments.R:47. Manifest consumers: sceptre-pipeline/bin/process_grna_assignments.R.
  Contract: Assigned object -> panels/combined plot from raw source assignments and UMI counts. Low-MOI maximum excludes flagged cells for one panel; other summaries remain pre-QC. Fixed-seed sampling/jitter; no scientific result mutation.
- `plot_run_calibration_check` (line 531): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-plotting_functions.R:140.
  Contract: Calibrated object -> QQ/band/effect-summary panels or combined cowplot, reading processed calibration fields. Public return_indiv_plots prose reverses actual TRUE-list/FALSE-combined behavior.
- `compare_calibration_and_discovery_results` (line 645): missing named internal contract. Call candidates: R/plotting_functions.R:863; R/plotting_functions.R:874.
  Contract: Two result frames + p threshold/scale/display options -> QQ ggplot. Equal per-panel group sample sizes are required by the QQ band; caller drops unequal calibration data. No facade writes.
- `make_volcano_plot` (line 723): missing named internal contract. Call candidates: R/plotting_functions.R:888.
  Contract: Result frame + threshold/limits -> ggplot; clips displayed p/effect values only. Current rejection color uses strict p < threshold, so the maximum significant p-value can be colored non-rejected.
- `plot_run_discovery_analysis` (line 836): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-plotting_functions.R:163.
  Contract: Discovery object -> QQ/volcano/text panels or combined plot. Drops NA rows, hides calibration QQ overlay if row counts differ, derives threshold from largest significant p. Fixed-seed downsampling via helper.
- `plot_covariates` (line 992): existing public roxygen; accuracy findings tracked separately. Call candidates: /Users/ekatsevi/code/research/sceptre-pipeline/bin/run_qc_nuclear.R:36; /Users/ekatsevi/code/research/sceptre-pipeline/bin/run_qc.R:46. Manifest consumers: sceptre-pipeline/bin/run_qc.R, sceptre-pipeline/bin/run_qc_nuclear.R.
  Contract: Imported/QC facade -> histograms using stored thresholds once QC was run; source-cell covariate rows, not retained-only view. No scientific state changes.
- `plot_covariates::make_histogram` (line 1005): missing named internal contract. Call candidates: R/plotting_functions.R:1042; R/plotting_functions.R:1048; R/plotting_functions.R:1056.
- `plot_run_qc` (line 1134): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-plotting_functions.R:186; /Users/ekatsevi/code/research/sceptre-pipeline/bin/run_qc.R:47. Manifest consumers: sceptre-pipeline/bin/run_qc.R.
  Contract: QC facade -> cellwise plot and optional pairwise plot/list. Currently hardcodes pairwise downsample/size/transparency rather than forwarding user arguments.
- `plot_cellwise_qc` (line 1165): missing named internal contract. Call candidates: R/plotting_functions.R:1144; /Users/ekatsevi/code/research/sceptre-pipeline/bin/run_qc_nuclear.R:37. Manifest consumers: sceptre-pipeline/bin/run_qc_nuclear.R.
  Contract: QC facade -> ggplot of overlapping per-rule/source-cell removal percentages. Metrics order is positional; high MOI omits multiple-guide bar. Applies the same theme twice.
- `plot_pairwise_qc` (line 1207): missing named internal contract. Call candidates: R/plotting_functions.R:1146.
  Contract: QC facade -> single ggplot of discovery pair counts/eligibility, fixed-seed row sampling. return_indiv_plots formal is unused. Repeated pair rows remain distinct plot rows.
- `plot_run_power_check` (line 1322): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-plotting_functions.R:201.
  Contract: Power object -> clipped positive/negative-control p-value jitter ggplot; negative controls are sampled per p-value bin; fixed seed. Reads results, no statistical recomputation.
- `downsample_result_data_frame` (line 1388): missing named internal contract. Call candidates: R/plotting_functions.R:885; R/plotting_functions.R:1338.
  Contract: Result frame -> grouped tibble sampled separately within p-value bins under fixed seed. downsample_pairs is a per-bin cap, not a global cap; grouping/order are observable.
- `plot_response_grna_target_pair` (line 1446): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/test-plotting_functions.R:208.
  Contract: QC facade + response/group -> normalized source counts through analysis-cell memberships and NT/complement controls, ggplot with optional result p annotation. Bonferroni guide-level groups versus target-level result IDs are not explicitly handled by singleton-only branch.

## R/preprocessing_functions.R

Four private helpers covering full-row finite design, CSC cell summaries, transient native banks and integration joins. Four drafted contracts; preserve fixed native draw order and no pointer persistence.

- `convert_covariate_df_to_design_matrix` (line 1): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/mixture_model_functs.R:25; R/s4_analysis_functs_1.R:138; tests/testthat/test-preprocessing_functs.R:80; tests/testthat/test-preprocessing_functs.R:89; +4 more in JSON.
- `compute_cell_covariates` (line 36): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/aux_functions.R:57; R/aux_functions.R:66.
- `get_synthetic_permutation_idxs` (line 76): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/analysis.R:19.
- `update_dfs_based_on_grouping_strategy` (line 140): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:178.
- `update_dfs_based_on_grouping_strategy::update_targeting_pairs` (line 153): missing named internal contract. Call candidates: R/preprocessing_functions.R:177; R/preprocessing_functions.R:181.

## R/qq_plot_helpers.R

Three private @noRd helper blocks exist, but layer return/data/scale/rank/grouping/true sampling-cap semantics need refinement. Three named ggproto methods and two local transform closures inherit enclosing contracts.

- `StatQQBand::setup_data` (line 8): missing named internal contract. Call candidates: none statically located.
- `StatQQBand::compute_group` (line 36): missing named internal contract. Call candidates: none statically located.
- `stat_qq_band` (line 99): existing @noRd; QQ semantics need refinement. Call candidates: R/plotting_functions.R:569; R/plotting_functions.R:587; R/plotting_functions.R:670.
  Contract: ggplot layer arguments -> ribbon LayerInstance backed by StatQQBand. y arrives transformed; inverse/transform around beta order-statistic bands, rank ties first, panel grouping requirements, deterministic occupied-grid downsampling.
- `StatQQPoints::compute_group` (line 134): missing named internal contract. Call candidates: none statically located.
- `stat_qq_points` (line 191): existing @noRd; QQ semantics need refinement. Call candidates: R/plotting_functions.R:563; R/plotting_functions.R:581; R/plotting_functions.R:669.
  Contract: ggplot layer arguments -> point LayerInstance backed by StatQQPoints. Inverse/clip/re-transform y and expected quantiles, deterministic grid sampling. Existing max_pts_to_plot prose overstates a hard global cap.
- `revlog_trans` (line 230): existing @noRd; QQ semantics need refinement. Call candidates: R/plotting_functions.R:588; R/plotting_functions.R:589; R/plotting_functions.R:692; R/plotting_functions.R:693; +2 more in JSON.
  Contract: Positive logarithm base -> scales transformation object with inverse, log breaks and positive domain. Inner trans/inv closures inherit this contract; no standalone boilerplate needed.
- `revlog_trans::trans` (line 231): missing named internal contract. Call candidates: none statically located.
- `revlog_trans::inv` (line 234): missing named internal contract. Call candidates: none statically located.

## R/RcppExports.R

Twenty-one generated native-call wrappers. Inventory only; do not hand-document/edit generated code. Native owners must document C++ entrypoint argument/index/ownership contracts at their source.

- `estimate_theta` (line 4): generated; document native source, not this wrapper. Call candidates: R/fitting.R:10.
- `fisher_yates_sampler` (line 14): generated; document native source, not this wrapper. Call candidates: R/preprocessing_functions.R:94; R/preprocessing_functions.R:104; R/preprocessing_functions.R:129.
- `hybrid_fisher_iwor_sampler` (line 25): generated; document native source, not this wrapper. Call candidates: R/preprocessing_functions.R:115.
- `crt_index_sampler` (line 29): generated; document native source, not this wrapper. Call candidates: none statically located.
- `crt_index_sampler_fast` (line 33): generated; document native source, not this wrapper. Call candidates: R/association_workers.R:140; R/association_workers.R:221.
- `run_low_level_test_full_v4` (line 37): generated; document native source, not this wrapper. Call candidates: R/association_workers.R:19; R/association_workers.R:87; R/association_workers.R:161; R/association_workers.R:244.
- `run_reduced_em_algo_cpp` (line 41): generated; document native source, not this wrapper. Call candidates: R/mixture_model_functs.R:162.
- `sample_combinations_v2` (line 45): generated; document native source, not this wrapper. Call candidates: R/neg_control_functions.R:41.
- `iterate_over_combinations` (line 49): generated; document native source, not this wrapper. Call candidates: R/neg_control_functions.R:29.
- `increment_matrix` (line 53): generated; document native source, not this wrapper. Call candidates: R/neg_control_functions.R:98.
- `sample_undercover_pairs_v2` (line 57): generated; document native source, not this wrapper. Call candidates: R/neg_control_functions.R:56; R/neg_control_functions.R:82; tests/testthat/test-refactor-qc.R:20.
- `compute_n_trt_cells_matrix` (line 61): generated; document native source, not this wrapper. Call candidates: R/storage.R:111; tests/testthat/test-refactor-storage.R:37.
- `compute_genes_within_distance` (line 65): generated; document native source, not this wrapper. Call candidates: R/pair_constructor_functs.R:140.
- `evaluate_pair_eligibility_cpp` (line 69): generated; document native source, not this wrapper. Call candidates: R/pair_qc.R:10.
- `compute_nt_nonzero_matrix_and_n_ok_pairs_v3` (line 73): generated; document native source, not this wrapper. Call candidates: R/storage.R:79; tests/testthat/test-refactor-storage.R:22.
- `load_csr_row` (line 77): generated; document native source, not this wrapper. Call candidates: R/storage.R:13.
- `obtain_pointer_vector` (line 81): generated; document native source, not this wrapper. Call candidates: R/import_functs.R:852; R/storage.R:224; R/storage.R:229.
- `compute_cell_covariates_cpp` (line 85): generated; document native source, not this wrapper. Call candidates: R/storage.R:153.
- `compute_n_grnas_per_cell_vector` (line 89): generated; document native source, not this wrapper. Call candidates: R/assign_grna_functions.R:146.
- `increment_vector` (line 93): generated; document native source, not this wrapper. Call candidates: R/getters_and_setters.R:116; R/import_functs.R:508; R/import_functs.R:509; R/import_functs.R:764; +1 more in JSON.
- `threshold_count_matrix` (line 97): generated; document native source, not this wrapper. Call candidates: R/storage.R:25; tests/testthat/test-refactor-storage.R:19; tests/testthat/test-refactor-storage.R:20.

## R/s4_analysis_functs_1.R

Three public facades and two private staged QC boundaries consumed by pipeline. Facade flags/caches/slot writes remain compatible. Two drafted contracts.

- `set_analysis_parameters` (line 67): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/helper-refactor.R:368; tests/testthat/helper-refactor.R:439; tests/testthat/helper-refactor.R:447; tests/testthat/helper-refactor.R:515; +36 more in JSON. Manifest consumers: sceptre-pipeline/bin/set_analysis_parameters.R.
- `assign_grnas` (line 255): existing public roxygen; accuracy findings tracked separately. Call candidates: R/s4_analysis_functs_2.R:641; tests/testthat/helper-refactor.R:376; tests/testthat/helper-refactor.R:451; tests/testthat/test-assign_grnas.R:69; +30 more in JSON. Manifest consumers: sceptre-pipeline/bin/assign_grnas.R, sceptre-pipeline/bin/process_grna_assignments.R.
- `run_qc` (line 406): existing public roxygen; accuracy findings tracked separately. Call candidates: R/s4_analysis_functs_2.R:646; R/s4_analysis_functs_2.R:652; tests/testthat/helper-refactor.R:377; tests/testthat/helper-refactor.R:453; +33 more in JSON. Manifest consumers: sceptre-pipeline/bin/run_qc.R.
- `run_qc_pt_1` (line 428): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:415. Manifest consumers: sceptre-pipeline/bin/run_qc_nuclear.R.
- `run_qc_pt_2` (line 508): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_1.R:424; /Users/ekatsevi/code/research/sceptre-pipeline/bin/run_association_analysis_nuclear.R:52. Manifest consumers: sceptre-pipeline/bin/run_association_analysis_nuclear.R.

## R/s4_analysis_functs_2.R

Five public facades and three private automatic/prepared calibration boundaries. Pipeline retains raw pod output for merge. Three drafted contracts.

- `run_calibration_check` (line 60): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/helper-refactor.R:418; tests/testthat/test-parallel_mac.R:39; tests/testthat/test-plotting_functions.R:79; tests/testthat/test-refactor-cache.R:8; +7 more in JSON.
- `run_calibration_check_pt_1` (line 89): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_2.R:72. Manifest consumers: sceptre-pipeline/bin/prepare_association_analyses.R.
- `run_calibration_check_pt_2` (line 165): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_2.R:78. Manifest consumers: sceptre-pipeline/bin/run_association_analysis.R.
- `run_power_check` (line 237): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/helper-refactor.R:422; tests/testthat/test-parallel_mac.R:44; tests/testthat/test-plotting_functions.R:80; tests/testthat/test-s4_analysis_functs_2.R:408; +1 more in JSON. Manifest consumers: sceptre-pipeline/bin/run_association_analysis.R.
- `run_discovery_analysis` (line 340): existing public roxygen; accuracy findings tracked separately. Call candidates: tests/testthat/helper-refactor.R:424; tests/testthat/helper-refactor.R:434; tests/testthat/helper-refactor.R:435; tests/testthat/helper-refactor.R:438; +17 more in JSON. Manifest consumers: sceptre-pipeline/bin/run_association_analysis.R, sceptre-pipeline/bin/run_association_analysis_nuclear.R.
- `get_result` (line 457): existing public roxygen; accuracy findings tracked separately. Call candidates: R/s4_analysis_functs_2.R:607; tests/testthat/test-s4_analysis_functs_2.R:508.
- `write_outputs_to_directory` (line 534): existing public roxygen; accuracy findings tracked separately. Call candidates: none statically located.
- `skip_assign_grnas_and_run_qc` (line 628): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/s4_analysis_functs_2.R:97; R/s4_analysis_functs_2.R:247; R/s4_analysis_functs_2.R:350.

## R/s4_class.R

Persistent 56-slot class schema; commented-out class unions are stale, but apparently unused slots cannot be removed under saved-object compatibility. @import Matrix must remain.


## R/s4_helpers.R

show/print/plot S4 methods plus one private row loader. Public help contains old assignment names; print has two dead locals. One drafted load_row contract, public methods/text unchanged.

- `<anonymous@83>::get_mark` (line 88): missing named internal contract. Call candidates: R/s4_helpers.R:95.
- `load_row` (line 363): baseline missing; drafted in paused checkpoint a2e66f0. Call candidates: R/mixture_model_functs.R:78; R/plotting_functions.R:97; R/plotting_functions.R:325; R/plotting_functions.R:1484; +2 more in JSON.

## R/sceptre.R

Package-level registration/globalVariables declarations and public examples. Duplicate global variable names and narrative style are low-risk metadata cleanup; preserve registration/import semantics.


## R/test_functs.R

Two unexported dimname helpers shipped in namespace solely for testing/no known callers. Candidate for relocation/removal after private-consumer review.

- `set_colnames` (line 1): missing named internal contract. Call candidates: tests/testthat/helper-mock-data.R:203.
  Contract: Test-only dimname adapter: matrix-like x + names/NULL -> value with column names replaced. Called only from helper-mock-data within inspected scope; currently unnecessarily ships in package namespace.
- `set_rownames` (line 6): missing named internal contract. Call candidates: none statically located.
  Contract: Test-only dimname adapter; no located callers in R/, testthat tests/helpers or pinned pipeline bin scripts. Removal requires a declared private-consumer policy, not merely assuming global non-use.

## tests/testthat/helper-mock-data.R

Eight named test helpers; four have detailed roxygen prose and four detailed ordinary comments. Random residual matrix/covariate generation consumes R RNG; malformed-parameter checks have an NT-pattern typo. Keep fixture inputs/distributions separate from production cleanup.

- `.make_mock_target_and_chr_only` (line 13): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-mock-data.R:119; tests/testthat/test-mock-data-functions.R:13; tests/testthat/test-mock-data-functions.R:43.
- `.make_mock_cross_target_chr_with_distances` (line 44): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-mock-data.R:123; tests/testthat/test-mock-data-functions.R:47.
- `make_mock_grna_target_data` (line 101): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-assign_grnas.R:5; tests/testthat/test-aux_functions.R:9; tests/testthat/test-check_functions.R:7; tests/testthat/test-mock-data-functions.R:122; +6 more in JSON.
- `.make_mock_patterned_matrix` (line 175): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-mock-data.R:280; tests/testthat/helper-mock-data.R:293; tests/testthat/test-mock-data-functions.R:171; tests/testthat/test-mock-data-functions.R:177; +8 more in JSON.
- `.make_mock_matrix_pattern_list` (line 242): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-mock-data.R:398; tests/testthat/helper-mock-data.R:416; tests/testthat/helper-mock-data.R:555; tests/testthat/test-mock-data-functions.R:330; +2 more in JSON.
- `make_mock_grna_matrices` (line 341): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-assign_grnas.R:9; tests/testthat/test-aux_functions.R:15; tests/testthat/test-aux_functions.R:54; tests/testthat/test-check_functions.R:22; +10 more in JSON.
- `make_mock_response_matrices` (line 533): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-assign_grnas.R:16; tests/testthat/test-aux_functions.R:10; tests/testthat/test-aux_functions.R:30; tests/testthat/test-check_functions.R:17; +5 more in JSON.
- `make_mock_extra_covariates_data_frames` (line 616): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-aux_functions.R:22; tests/testthat/test-check_functions.R:29; tests/testthat/test-check_functions.R:420; tests/testthat/test-check_functions.R:436; +4 more in JSON.

## tests/testthat/helper-refactor-native.R

Five direct named helpers plus cached bridge factory. Production calls use getFromNamespace; test-only sourceCpp side effect/cache and immutable sampler provenance. Current bridge cache is session-singleton; path/cache policy is a harness issue, not runtime redesign.

- `refactor_native_function` (line 1): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor-native.R:46; tests/testthat/helper-refactor-native.R:50.
- `refactor_native_platform` (line 17): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-native.R:105.
- `refactor_native_run_case` (line 37): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-native.R:10; tests/testthat/test-refactor-native.R:29; tests/testthat/test-refactor-native.R:45; tests/testthat/test-refactor-native.R:47.
- `refactor_native_sample` (line 49): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-native.R:60; tests/testthat/test-refactor-native.R:61; tests/testthat/test-refactor-native.R:97; tests/testthat/test-refactor-native.R:109.
- `refactor_native_reference` (line 55): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-native.R:87.

## tests/testthat/helper-refactor.R

Twenty-five direct named helpers plus two cached factories, local diagnostic/capture helpers. Reads/hashes/validates oracle files and compiler configuration; native compilation is generator-owned, never automatic fixture regeneration in tests. Frozen inputs and source provenance are immutable.

- `refactor_fixture_dir` (line 2): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:136; tests/testthat/helper-refactor.R:137; tests/testthat/helper-refactor.R:211; tests/testthat/helper-refactor.R:228; +27 more in JSON.
- `refactor_private` (line 6): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:476; tests/testthat/helper-refactor.R:485; tests/testthat/helper-refactor.R:486; tests/testthat/test-refactor-contracts.R:348; +16 more in JSON.
- `refactor_oracle_file` (line 10): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:30; tests/testthat/helper-refactor.R:30; tests/testthat/helper-refactor.R:37; tests/testthat/helper-refactor.R:57.
- `refactor_oracle_makevars` (line 20): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:55.
- `refactor_oracle_toolchain::binary_state` (line 35): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:63; tests/testthat/helper-refactor.R:83.
- `refactor_oracle_environment` (line 107): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:139; tests/testthat/helper-refactor.R:203; tests/testthat/helper-refactor.R:256; tests/testthat/test-refactor-contracts.R:144; +6 more in JSON.
- `refactor_oracle_environment_differences` (line 121): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:129; tests/testthat/helper-refactor.R:222; tests/testthat/test-refactor-contracts.R:94; tests/testthat/test-refactor-contracts.R:96.
- `refactor_oracle_environment_differences::named` (line 123): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:124; tests/testthat/helper-refactor.R:124.
- `refactor_oracle_environment_differences::field` (line 127): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:129; tests/testthat/helper-refactor.R:129.
- `refactor_frozen_environment` (line 135): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:197; tests/testthat/helper-refactor.R:259; tests/testthat/test-refactor-contracts.R:18; tests/testthat/test-refactor-contracts.R:54; +3 more in JSON.
- `refactor_oracle_input_files` (line 153): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:226.
- `refactor_oracle_output_files` (line 157): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:212; tests/testthat/test-refactor-contracts.R:165; tests/testthat/test-refactor-contracts.R:166.
- `refactor_oracle_hashes` (line 163): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:170; tests/testthat/helper-refactor.R:228; tests/testthat/helper-refactor.R:233; tests/testthat/helper-refactor.R:238; +1 more in JSON.
- `refactor_oracle_verify_source` (line 167): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:236; tests/testthat/test-refactor-contracts.R:202; tests/testthat/test-refactor-contracts.R:205.
- `refactor_expected` (line 246): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-contracts.R:232; tests/testthat/test-refactor-contracts.R:271; tests/testthat/test-refactor-contracts.R:286; tests/testthat/test-refactor-contracts.R:293; +4 more in JSON.
- `refactor_api` (line 267): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-contracts.R:226.
- `refactor_api::describe_formals` (line 272): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:284.
- `refactor_make_inputs` (line 290): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: none statically located.
- `refactor_case_definitions` (line 331): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: none statically located.
- `refactor_expand_inputs` (line 343): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: none statically located.
- `refactor_import` (line 356): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:367; tests/testthat/test-refactor-contracts.R:357.
- `refactor_prepare` (line 366): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:416; tests/testthat/helper-refactor.R:433; tests/testthat/helper-refactor.R:458; tests/testthat/helper-refactor.R:508; +6 more in JSON.
- `refactor_plain` (line 385): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:402.
- `refactor_observe` (line 394): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:417; tests/testthat/helper-refactor.R:421; tests/testthat/helper-refactor.R:426; tests/testthat/helper-refactor.R:464; +9 more in JSON.
- `refactor_capture` (line 405): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:415; tests/testthat/helper-refactor.R:432; tests/testthat/helper-refactor.R:500; tests/testthat/test-refactor-contracts.R:284; +3 more in JSON.
- `refactor_run_case` (line 414): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-contracts.R:234; tests/testthat/test-refactor-contracts.R:302; tests/testthat/test-refactor-contracts.R:303; tests/testthat/test-refactor-contracts.R:325.
- `refactor_cache_observations` (line 431): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-contracts.R:270.
- `refactor_index_observations` (line 474): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-contracts.R:259.
- `refactor_edge_observations` (line 498): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/test-refactor-contracts.R:292.
- `refactor_edge_observations::outcome` (line 499): test helper; plain/roxygen coverage varies, module notes apply. Call candidates: tests/testthat/helper-refactor.R:508; tests/testthat/helper-refactor.R:509; tests/testthat/helper-refactor.R:513.

## Finding register

Priorities concern cleanup planning and existing source debt; no new production regression or unexecuted reproducer is claimed. A compatibility refactor must not absorb behavior fixes without explicit scope approval.

### L01 (P2) — Internal contracts are absent for 60/63 named private functions at c635263; three QQ helpers have partial @noRd blocks. Anonymous callbacks do not need boilerplate.

Evidence: `all included handwritten R modules`.

Proposed action: Review paused 39-function checkpoint, then document the remaining 21 functions and shared QQ/method/module contracts with types, views, order, duplicate policy, returns and effects.

Behavior boundary: Documentation/comments only; parsed expressions, exports, public formals and saved schema must remain identical.

Uncertainty: Draft accuracy still requires senior review; no claim of completed full-source cleanup.

Validation: AST identity, complete function checklist, no generated help/export churn for internal @noRd work.

### L02 (P2) — Public plot return prose reverses TRUE-list/FALSE-combined behavior in several functions; generic plot help uses obsolete run_assign_grnas/plot_run_assign_grnas names; ymax QQ prose says less instead of greater.

Evidence: `R/s4_helpers.R:311-312`, `R/plotting_functions.R:258-264,498-504,800-805,973-979,1105-1110`, `R/qq_plot_helpers.R:178-187`.

Proposed action: Prepare one explicitly reviewed public-doc synchronization patch, generating only matching Rd text with unchanged aliases/formals.

Behavior boundary: Documentation correction preserves current results; unlike @noRd work it intentionally changes help text and needs a narrow generated-file whitelist.

Uncertainty: Review every affected return variant; plot_run_qc returns a single cellwise plot if no discovery pairs even with list request.

Validation: API inventory exact; Rd usages/aliases unchanged; evaluate return-class branches without changing code.

### L03 (P3) — Commented-out class unions, duplicate globalVariables entries and broken numbered/narrative comments remain outside or before the checkpoint.

Evidence: `R/s4_class.R:1-7`, `R/sceptre.R:1-35`, `ten checkpoint files`.

Proposed action: Remove dead comments/duplicates and retain concise explanations of compatibility/index decisions.

Behavior boundary: Comment deletion can be AST-identical; globalVariables deduplication is metadata behavior and should be separately reviewed. Never remove S4 slots.

Uncertainty: Some apparently unused slots, including mitochondrial_gene, persist for saved-object compatibility.

Validation: AST diff for comments; preserve 56-slot/type inventory and NAMESPACE registration.

### L04 (P3) — Confirmed unused print locals, unused check_discovery_analysis_inputs calibration_check_run formal, unused private plot return flag, redundant theme addition and repeated processor guards.

Evidence: `R/s4_helpers.R:104,106`, `R/check_functions.R:239,656`, `R/plotting_functions.R:1193,1198,1212`.

Proposed action: Remove proven dead locals; review unused private formals/call arguments and factor shared processor validation only if messages/order remain identical.

Behavior boundary: Executable cleanup, not documentation. Removing formals changes private call shape; preserve known pipeline/private consumers or explicitly approve that change.

Uncertainty: No located caller uses calibration_check_run; dynamic/private external use is not globally excluded. Redundant getter call removal can change malformed-object errors.

Validation: Current callers/consumer inventory; supported workflow/result/print compatibility; focused invalid-input message/order checks if validators are consolidated.

### L05 (P3) — The NT assignment helper is a one-call forwarder and the matrix writer repeats directory preparation already done by its sole import caller; v2 naming remains in handwritten calibration adapter.

Evidence: `R/cellwise_qc_functs.R:124-126`, `R/import_functs.R:881-892`, `R/s4_analysis_functs_2.R:151-155`.

Proposed action: Consider direct component calls/canonical private names after caller review; avoid adding compatibility wrappers solely to hide avoidable indirection.

Behavior boundary: Keep algorithm, native entrypoints, saved slot layout and externally consumed pipeline functions stable.

Uncertainty: Private downstream consumers beyond the inspected repositories may exist.

Validation: Static/dynamic consumer search, unchanged input ownership/source maps/backed import output and RNG characterization.

### L06 (P2) — construct_trans_pairs calls setDT on caller-supplied positive_control_pairs, changing a data.frame class by reference.

Evidence: `R/pair_constructor_functs.R:274-278`.

Proposed action: Decide explicit ownership policy; use a local owned conversion if removing this side effect is approved.

Behavior boundary: Removing caller mutation is observable behavior and must not be mislabeled AST-only cleanup.

Uncertainty: Supported numeric pair results should match, but callers relying on incidental class mutation would differ.

Validation: Serialize/class-check caller before/after and compare returned class, row order and values under every exclusion mode.

### L07 (P2) — Vector aggregation materializes a dense guide matrix; the R-backed import path does not perform vector_id aggregation although the in-memory path does.

Evidence: `R/import_functs.R:1008-1027`, `R/import_functs.R:152-205`.

Proposed action: Document current memory/aggregation boundaries now; scope later storage-aware aggregation separately rather than burying it in docs or claiming parity.

Behavior boundary: Sparse aggregation or backed parity changes execution and potentially functionality/numeric summation order; outside pure cleanup.

Uncertainty: The intended cross-path vector contract needs product/statistical owner confirmation.

Validation: Existing vector import cases plus an approved paired in-memory/backed workload before any future behavior change.

### L08 (P2) — CSC construction directly adopts zero-based scratch vectors without sorting/coalescing/shape validation; Cell Ranger modality indexing subtracts a minimum and assumes contiguous feature blocks.

Evidence: `R/import_functs.R:503-513,851-859`.

Proposed action: Make parser preconditions explicit; propose an owned validated triplet/index boundary only after examining interleaved modalities/duplicate coordinate policy.

Behavior boundary: Sorting, coalescing or accepting interleaved features can change behavior; no silent normalization in this cleanup.

Uncertainty: Standard bundled inputs satisfy assumptions; frequency/validity of interleaved external Cell Ranger inputs is unverified.

Validation: Approved parser edge fixtures: dimensions/order/duplicate policy, explicit zeros, interleaved modality IDs, caller scratch ownership; existing import fixtures unchanged.

### L09 (P3) — Matrix Market metadata reader assumes a well-formed comment/header line and closes its connection only on the success path.

Evidence: `R/import_functs.R:862-878`.

Proposed action: Separate metadata parsing from file lifetime; use guaranteed connection cleanup and clearly state accepted spacing/header contract.

Behavior boundary: Success-path resource cleanup can preserve values; validation/error improvements alter malformed-input behavior and require explicit scope.

Uncertainty: Do not invent support for malformed Matrix Market variants; gzip behavior already exercised.

Validation: Existing gz/plain imports plus small malformed/truncated-header/connection-lifetime probes if approved.

### L10 (P2) — Calibration group size is not validated; zero-only pair-count guard accepts other invalid forms; QC omits validation of additional cell indices and mito cutoff, and alpha is not checked by set-parameters validator.

Evidence: `R/check_functions.R:565-649`, `R/s4_analysis_functs_1.R:406-425`, `R/s4_analysis_functs_2.R:89-149`.

Proposed action: Inventory exact accepted/precondition sets; decide a separate input-validation issue or explicitly approved cleanup stage.

Behavior boundary: New rejection/coercion changes public invalid-input behavior; not allowed to slip into a compatibility refactor.

Uncertainty: Severity for malformed sizes depends on native downstream constraints; no new crash claim without a reproducer.

Validation: Targeted scalar/NA/integer/range cases plus immutable existing edge behavior to distinguish intended fixes from accidental changes.

### L11 (P2) — g_pert_guess_range upper-bound validation reads pi_guess_range > 20 instead of the perturbation range.

Evidence: `R/check_functions.R:521-525`.

Proposed action: Track a separate precise validation bug fix; preserve the legacy check while only documenting refactor contracts unless approved.

Behavior boundary: Changing accepted input set is an observable bug fix.

Uncertainty: Source-level defect confirmed; downstream practical effect for values >20 was not executed in this audit.

Validation: Small validator test with valid pi range and perturbation upper bound >20, plus valid-boundary cases after approval.

### L12 (P2) — Pair validator subsets grna_target/response_id to check NAs before checking that those columns exist.

Evidence: `R/check_functions.R:247-274`.

Proposed action: Decide whether the intended named-column error should precede NA checking; separate error-contract fix.

Behavior boundary: Reordering changes error text/class on malformed public input.

Uncertainty: Valid scientific workflows unaffected; stored edge cases may not cover this malformed table.

Validation: Missing-column, NA-value and valid-pair cases with explicit expected error precedence.

### L13 (P2) — plot_run_qc ignores user downsample_pairs, point_size and transparency, passing constants to plot_pairwise_qc.

Evidence: `R/plotting_functions.R:1134-1152`.

Proposed action: Track/approve a distinct argument-forwarding fix; document actual current boundary while planning.

Behavior boundary: Nondefault exported plotting behavior changes; default output can remain the same.

Uncertainty: Confirmed by source and unused-formal analysis; no plot output was regenerated during inventory.

Validation: A small forwarding probe or ggplot layer-data comparison for nondefault values; default visualization unchanged.

### L14 (P2) — Volcano coloring uses p < p_thresh although p_thresh is the maximum significant p-value, excluding that boundary from the rejection color.

Evidence: `R/plotting_functions.R:734`, `R/plotting_functions.R:851-855`.

Proposed action: Separate visual significance-boundary bug fix with explicitly agreed tie convention.

Behavior boundary: Presentation behavior changes while scientific significant flags/results stay fixed.

Uncertainty: Depends on whether the intended plotted rejection set follows significant flags exactly; confirm with owner.

Validation: A two-point/tied-threshold plot-data test against agreed significant flags.

### L15 (P2) — Guide-count binning branches on max count >10 even when a larger supplied threshold makes max_single_bin exceed the data; log2 can receive a negative argument.

Evidence: `R/plotting_functions.R:115-131`.

Proposed action: Track/approve a focused binning boundary fix; retain current formulas during documentation-only work.

Behavior boundary: Fixes a valid nondefault plotting input, not refactor-only behavior preservation.

Uncertainty: Source arithmetic demonstrates risk; no reproduction claimed beyond inspection.

Validation: Guide max count between 10 and a high explicit threshold, and thresholds/defaults around the branch.

### L16 (P3) — Downsample helpers cap rows per occupied bin/grid rather than enforce a global max; result-frame helper returns a grouped tibble.

Evidence: `R/plotting_functions.R:1388-1396`, `R/qq_plot_helpers.R:74-85,91-96,179-185`.

Proposed action: Correct contracts and misleading count prose; preserve sampling/RNG/order/grouping in cleanup.

Behavior boundary: Ungrouping or imposing a true global cap changes values/classes/performance and is a separate change.

Uncertainty: Exact maximal QQ occupied-bin count depends on monotonicity/scales; do not promise a universal total cap.

Validation: Existing fixed-seed plots and grouped result class/order; typed internal docs sufficient for the first stage.

### L17 (P2) — Response/group plot treats only singleton integration as guide-level; Bonferroni uses guide-level assignment groups but target-level result IDs.

Evidence: `R/plotting_functions.R:1453,1467-1474,1531-1542`.

Proposed action: Clarify current supported plot input/annotation behavior and separately scope Bonferroni plotting semantics.

Behavior boundary: Changing guide-versus-target selection/annotation changes exported functionality.

Uncertainty: A source-level ID mismatch is visible; intended pooling/plot labeling requires owner decision.

Validation: Approved Bonferroni guide/target example, both control views, matching p-value annotation/result IDs.

### L18 (P3) — Assignment getter constructs CSR slots manually and decrements unlist scratch through a native mutator.

Evidence: `R/getters_and_setters.R:106-133`.

Proposed action: Audit and explicitly test scratch ownership before simplifying it; preserve source/QC views and repeated calls.

Behavior boundary: An ownership fix should not silently change sparse class, row/column IDs or raw assignments.

Uncertainty: No mutation bug is asserted; single-guide alias behavior needs a focused probe.

Validation: Caller raw serialization unchanged across one-guide/repeated get calls and optional QC columns.

### L19 (P3) — Production namespace contains test-only set_colnames and unreferenced set_rownames; NT mock-pattern validation checks length(non_nt_patterns) instead of nt_patterns.

Evidence: `R/test_functs.R:1-9`, `tests/testthat/helper-mock-data.R:203,381-384`.

Proposed action: Relocate the used dimname helper into test support and consider removing the unreferenced helper after private-consumer review; treat mock validation correction separately.

Behavior boundary: Package-private symbol removal and test-generation behavior are distinct from runtime numeric changes. Never regenerate frozen fixture inputs under this cleanup.

Uncertainty: Caller absence is scoped to R/tests/pinned pipeline, not every external user.

Validation: Full existing test helpers/fixtures; explicit old input checksums; no automatic fixture regeneration.

### L20 (P3) — QQ setup requires equal group sizes within each panel, then selects group == min(group) globally; panels with distinct group IDs may lose their confidence band.

Evidence: `R/qq_plot_helpers.R:8-37`.

Proposed action: Record the per-panel grouping contract; scope any panel-local selection correction separately.

Behavior boundary: Changes ggplot data returned for some facet inputs; not pure cleanup.

Uncertainty: Potential facet-specific issue inferred from source; needs a small ggplot-layer reproducer before being labeled confirmed runtime failure.

Validation: Two panels with distinct group IDs and equal per-panel sample sizes, plus existing single-panel QQ behavior.

## Suggested order for the parent addendum (proposal only)

1. Senior-review the function/module checklist and paused checkpoint accuracy before resuming cleanup.
2. Apply reviewed internal contracts/comment cleanup with AST identity and no generated-help/API changes. Keep named closures brief and shared contracts explicit.
3. In a separate approved executable-cleanup batch, remove only proven unused locals/indirection/test-only helpers and address explicit ownership, preserving public and known pipeline/private contracts.
4. Keep validation/plot/parser behavior bugs as individually approved fixes or follow-up issues; do not silently change invalid-input acceptance, selection order, sampling, output grouping or sparse zeros.
5. Reuse the preserved 3,522 installed checks, 112 native checks, 597 latest-harness checks and baseline-warning package qualification for documentation-only work. Require targeted meaningful tests for each approved executable change, then the relevant existing gates; no fixture regeneration from candidate code.

Artifacts: `inventory.json` contains every function/callback, caller candidate and checkpoint contract text; `findings.json` contains the material register. This audit does not implement the plan or claim external private-caller completeness.

## Source hash coverage and behavior split

SHA256 is recorded for all 24 included files in `source-sha256.json` and `source-sha256.tsv`; the Git source archive hash is in `inventory.json`. Every named function line was checked against its declaration.

- AST or documentation only: L01, L02, L03, L16.
- executable engineering requires reviewed behavior boundary: L04, L05, L06, L09, L18, L19.
- valid input behavior bugs confirmed by source: L13, L14, L15.
- invalid input or validation behavior changes: L10, L11, L12.
- potential or policy dependent behavior requires reproducer or owner decision: L07, L08, L17, L20.

The identified mutation, range-check, plot-forwarding, volcano-boundary and histogram-branch snippets were cross-checked against pinned original 21f9ea0; they already existed there. No new refactor regression is asserted. See `original-behavior-crosscheck.json`. No runtime reproducer or changed output is claimed for source-only findings.
