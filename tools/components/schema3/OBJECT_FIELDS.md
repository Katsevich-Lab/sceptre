# Final schema 3 object fields

Generated from `methods::getSlots()` by `generate-field-tables.R`. Declared
types are exact. These tables describe the live schema, not legacy aliases.

`CountStore` is the union of `dgRMatrix`, `BackedCountStore` and `NULL`.
`CountHandle` is the union of `externalptr` and `NULL`. Neither union names
an optional dependency class; loading the package does not load ondisc.

## sceptre_object

| Field | Exact declared type | Meaning |
|---|---|---|
| `data` | `ScreenData` | Imported screen data and backing-store identity. |
| `spec` | `AnalysisSpec` | Resolved requested scientific and configuration choices. |
| `preparation` | `PreparationRecord` | Durable design, source assignments, QC payloads and completed preparation records. |
| `fits` | `FitCache` | Compact response and CRT perturbation fits with independently validated actual fitting identities. |
| `results` | `ResultRecord` | Completed tables with their own provenance, validity identity and history. |
| `schema_version` | `integer` | Saved schema version; scalar integer 3. |

## ScreenData

| Field | Exact declared type | Meaning |
|---|---|---|
| `response_names` | `character` | Original response labels aligned with response rows, used for mitochondrial summaries; empty when legacy labels are unknown. |
| `response_matrix` | `CountStore` | Typed CountStore: row-accessible dgRMatrix, or BackedCountStore owning a descriptor and optional live handle; NULL only before import or for legacy unavailable counts. |
| `grna_matrix` | `CountStore` | Typed CountStore with the same ordered source cells as responses; backed stores retain dimensions and labels when closed. |
| `covariate_data_frame` | `data.frame` | One row per source cell: imported count summaries and supplied extra covariates. |
| `grna_target_data_frame` | `data.frame` | Canonical imported guide-to-target metadata, preserving supplied row order and multiplicity. |
| `low_moi` | `logical` | Whether the screen is analyzed as low MOI. |
| `import_grna_assignment_info` | `list` | Source-cell-aligned maximum-guide labels and maximum-guide UMI fractions computed at import or supported replacement. |
| `identity` | `character` | Content/metadata fingerprint assigned at import or supported replacement; active backing files are immutable by contract. |

## AnalysisSpec

| Field | Exact declared type | Meaning |
|---|---|---|
| `discovery_pairs` | `data.frame` | Canonical requested discovery target-response rows, without derived guide expansion. |
| `positive_control_pairs` | `data.frame` | Canonical requested positive-control target-response rows, without derived guide expansion. |
| `formula_object` | `formula` | Resolved requested model formula; the actual compiled design belongs to preparation. |
| `side_code` | `integer` | Requested test direction: -1 left, 0 both, 1 right. |
| `resampling_approximation` | `character` | Requested tail policy: skew_normal or no_approximation. |
| `control_group_complement` | `logical` | TRUE selects complement controls; FALSE selects nontargeting controls. |
| `run_permutations` | `logical` | TRUE selects permutation resampling; FALSE selects CRT. |
| `n_nonzero_trt_thresh` | `integer` | Minimum treatment observations for pair eligibility. |
| `n_nonzero_cntrl_thresh` | `integer` | Minimum control observations for pair eligibility. |
| `grna_integration_strategy` | `character` | Requested guide integration: union, singleton or bonferroni. |
| `grna_assignment_method` | `character` | Resolved guide-assignment choice: maximum, thresholding or mixture. |
| `grna_assignment_hyperparameters` | `list` | Method-specific requested assignment settings, including any assignment-model formula. |
| `multiple_testing_alpha` | `numeric` | Requested reporting significance threshold; also determines empirical resampling precision. |
| `multiple_testing_method` | `character` | Requested complete-family multiple-testing correction method. |
| `cellwise_qc_thresholds` | `list` | Requested UMI/nonzero quantile ranges and mitochondrial-fraction cutoff. |
| `additional_cells_to_remove` | `integer` | Additional source-cell positions explicitly excluded by the requested QC configuration. |
| `perturbation_fit_method` | `character` | Optional experimental CRT logistic fitter; empty means the reference glm.fit and preserves historical default certificates. |

## PreparationRecord

| Field | Exact declared type | Meaning |
|---|---|---|
| `design_matrix` | `matrix` | Actual compiled model matrix, one row per source cell, with design columns, contrasts and ordering retained. |
| `nt_nonzero_matrix` | `matrix` | Cached NT-guide-by-response stored-entry counts used by calibration. |
| `control_nonzero_totals` | `integer` | Cached per-response stored-entry totals: all selected cells for complement controls, or summed NT-group counts for NT controls (overlaps count again). |
| `discovery_pairs_with_info` | `data.frame` | Compiled discovery rows with treatment/control counts and eligibility decisions, including failed-QC rows. |
| `positive_control_pairs_with_info` | `data.frame` | Compiled positive-control rows with treatment/control counts and eligibility decisions. |
| `guide_assignments` | `list` | Named individual-guide memberships in canonical source-cell coordinates, durable before QC. |
| `source_groups` | `list` | Target/guide groups and individual NT memberships in source-cell coordinates. |
| `analysis_groups` | `list` | Groups remapped to analysis/NT-fit coordinates after cell QC, including the ordered NT-to-analysis mapping. |
| `guide_counts_per_cell` | `integer` | Cached per-cell assignment counts used for low-MOI exclusions and assignment diagnostics. |
| `assignment_exclusions` | `integer` | Source cells excluded because assignments are absent or ambiguous under the chosen assignment method. |
| `analysis_cells` | `integer` | Ordered source-cell positions retained by cell QC. |
| `cell_removal_metrics` | `integer` | Named counts of QC exclusion reasons; reasons may overlap. |
| `stages` | `list` | Completed design/assignment/cell_qc/pair_qc records. An entry exists only after validation; its identity and prerequisites determine readiness. |
| `guide_groups` | `character` | Derived group labels aligned one-for-one with canonical guide metadata rows; no duplicated metadata table. |
| `discovery_pairs` | `data.frame` | Compiled discovery target-to-guide-group expansion, preserving requested row multiplicity. |
| `positive_control_pairs` | `data.frame` | Compiled positive-control target-to-guide-group expansion. |

## FitCache

| Field | Exact declared type | Meaning |
|---|---|---|
| `response_fits` | `list` | Named response fits: each contains fitted_coefs and scalar theta. Identifiers, coefficient/design compatibility and dispersion are validated. |
| `response_identity` | `list` | Actual data identity, ordered fitting cells, ordered design columns, design-value fingerprint and fitter version for response fits. |
| `perturbation_fits` | `list` | Named CRT logistic fits: compact fitted_coefs, five fit diagnostics and three hashes of actual ordered fit cells, group design and binary treatment. No cell-length vectors or sampled draws. |
| `perturbation_identity` | `list` | Shared ordered analysis source cells, design columns/value hash, logistic fitter method and arithmetic runtime. NT group-specific fitting inputs are certified by each entry without repeated cell vectors. |

## ResultRecord

| Field | Exact declared type | Meaning |
|---|---|---|
| `records` | `list` | Completed run entries keyed by analysis name. Each owns its table, frozen spec/run plan, scope, identity and previous run history; invalidated_by explains historical status. |
| `import_metadata` | `list` | Boundary-only migration diagnostics, such as old schema version, discarded fit count or a temporary legacy reopen hint; not workflow state. |

## RunPlan

| Field | Exact declared type | Meaning |
|---|---|---|
| `analysis` | `character` | Canonical scalar analysis mode: run_calibration_check, run_power_check or run_discovery_analysis. |
| `pairs` | `data.frame` | Owned prepared hypothesis rows and complete logical QC status, preserving order, failed-QC rows, repeated requests and table representation. Public result commits require the complete prepared family. |
| `B1` | `integer` | Realized scalar nonnegative integer stage-one resampling budget. |
| `B2` | `integer` | Realized scalar nonnegative integer stage-two resampling budget; the total budget is checked before native execution. |
| `B3` | `integer` | Realized scalar nonnegative integer stage-three resampling budget, frozen for this run independently of later calibration or reporting choices. |
| `n_calibration_pairs` | `integer` | Requested positive scalar count for calibration, which may exceed the selected qualifying family; integer(0) for other analyses. |
| `calibration_group_size` | `integer` | Positive scalar number of NT guides per calibration group; integer(0) for other analyses. |
| `identity` | `list` | Frozen data, ordered cells, actual design, assignments, choices and canonical execution-payload identity; object-context matching is validated separately. |

## CellView

| Field | Exact declared type | Meaning |
|---|---|---|
| `fit_to_analysis` | `integer` | Ordered positions in the selected analysis cells used for this fit. |
| `source_cells` | `integer` | The corresponding ordered original count-matrix columns. |
| `trt_fit_pos` | `integer` | Treatment positions within this fit. |
| `n_trt` | `integer` | Validated treatment-position count; transient numerical convenience. |

## BackedCountStore

| Field | Exact declared type | Meaning |
|---|---|---|
| `descriptor` | `list` | Authoritative backing path, anchored content_md5, integer dimensions, dimnames and native integer_id; validated on creation and reopening. |
| `handle` | `CountHandle` | Process-local external pointer when open, NULL when closed; never used as content identity. |

