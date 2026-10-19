# Schema 1 and 2 field decisions

This is a migration/ownership inventory, not a compatibility facade. Authority:
current `R/s4_class.R`, `R/schema_migration.R`, lifecycle/cache/result code; the
physical fixtures `sceptre-0.99.0-21f9ea0.rds` and
`sceptre-schema2-4aad80b.rds`; schema 1 source at `97d95da`; and schema 2 source
at `4aad80b`. Physical attributes were inspected because the loaded S4 class
registration cannot describe an older serialized layout.

The schema 1 fixture has **56 payload attributes plus `class`**. The schema 2
fixture has **11 coordinator attributes and 59 component attributes**, excluding
their `class` attributes. The 56 rows below account for every schema 1 field and
its schema 2 location. The following 14 rows account for all additional schema 2
attributes, for a total of 70. Nested list metadata is inventoried separately.

Notation: `D` = `data`, `S` = `spec`, `P` = `preparation`, `F` = `fits`, `R` =
`results`; bare paths denote coordinator fields. The schema 3 destination refers
to current authority, not a retained alias. “Historical snapshot” means saved
result provenance only, never a working setting or reusable cache identity.

## All 56 schema 1 fields

| Schema 1 field | Schema 2 path | Schema 3 authority / decision at migration |
|---|---|---|
| `response_matrix` | D.response_matrix | D.response_matrix, direct typed CountStore. Unwrap at most one old list element; convert ondisc state/descriptors at the boundary. |
| `grna_matrix` | D.grna_matrix | D.grna_matrix, same typed-store policy, aligned source cells. |
| `covariate_data_frame` | D.covariate_data_frame | Retained canonical source-cell covariates. |
| `covariate_matrix` | P.covariate_matrix | P.design_matrix; durable compiled values/attributes, source-cell rows. |
| `grna_target_data_frame` | D.grna_target_data_frame | Retained canonical guide metadata; remove schema 1 derived grna_group column. Group labels have separate derived authority. |
| `low_moi` | D.low_moi | Retained screen property. |
| `covariate_names` | D.covariate_names | Removed redundant name vector; current covariate column names are authoritative. Not archived separately. |
| `discovery_pairs` | S.discovery_pairs | S.discovery_pairs becomes canonical requested target-response rows; old expanded rows go to P.discovery_pairs. Prefer saved schema 2 originals; otherwise reconstruct ordered expansion blocks with multiplicity. |
| `positive_control_pairs` | S.positive_control_pairs | Same split/recovery policy into canonical S requests and compiled P pairs. |
| `formula_object` | S.formula_object | Retained resolved response-model formula. Actual compiled design governs fit identity; provenance snapshots store formula language without mutable environments. |
| `side_code` | S.side_code | Retained resolved test direction. |
| `resampling_approximation` | S.resampling_approximation | Retained scientific choice; realized budgets move to run plans. |
| `control_group_complement` | S.control_group_complement | Retained resolved complement-versus-NT control choice. |
| `run_permutations` | S.run_permutations | Retained resolved permutation-versus-CRT choice. |
| `n_nonzero_trt_thresh` | S.n_nonzero_trt_thresh | Retained pair-QC setting. |
| `n_nonzero_cntrl_thresh` | S.n_nonzero_cntrl_thresh | Retained pair-QC setting. |
| `B1` | S.B1 | Removed working slot. New transient plans derive it; old value enters migrated result effective_snapshot where a result exists. |
| `B2` | S.B2 | Same transient-plan/historical-snapshot policy. |
| `B3` | S.B3 | Same policy; no calibration/discovery shared mutable budget. Old value is not proof of the budget actually used by every old result. |
| `grna_integration_strategy` | S.grna_integration_strategy | Retained resolved union/singleton/Bonferroni choice. |
| `grna_assignment_method` | S.grna_assignment_method | Retained resolved assignment choice. |
| `grna_assignment_hyperparameters` | S.grna_assignment_hyperparameters | Retained method-specific list; fresh stages fingerprint its values. |
| `multiple_testing_alpha` | S.multiple_testing_alpha | Retained working reporting/testing choice; each new result freezes its plan/spec. Empirical B3 derives from this value. |
| `multiple_testing_method` | S.multiple_testing_method | Retained working correction choice; complete-family correction is validated on commit. |
| `cell_removal_metrics` | P.cell_removal_metrics | Retained durable QC summary. |
| `cellwise_qc_thresholds` | S.cellwise_qc_thresholds | Retained resolved cell-QC list. |
| `mitochondrial_gene` | D.mitochondrial_gene | Removed unused field. Existing mitochondrial summaries remain in covariates; missing import display labels are not fabricated from IDs. |
| `M_matrix` | P.M_matrix | P.nt_nonzero_matrix; retained NT-group-by-response nonzero summary. |
| `n_nonzero_tot_vector` | P.n_nonzero_tot_vector | P.control_nonzero_totals; retained response totals for the existing control universe. |
| `discovery_pairs_with_info` | P.discovery_pairs_with_info | Retained compiled pair-QC table including failed rows; completeness validated against P.discovery_pairs. |
| `positive_control_pairs_with_info` | P.positive_control_pairs_with_info | Same retained pair-QC authority for positive controls. |
| `negative_control_pairs` | S.negative_control_pairs | Removed working slot. New selected calibration family belongs to its transient plan; old values are historical effective_snapshot only. |
| `initial_grna_assignment_list` | P.initial_grna_assignment_list | P.guide_assignments; named durable guide-to-source-cell memberships before QC. |
| `grna_assignments_raw` | P.grna_assignments_raw | P.source_groups; durable grouped source-cell memberships. |
| `grna_assignments` | P.grna_assignments | P.analysis_groups; retained QC-derived memberships with explicit NT coordinate interpretation. |
| `grnas_per_cell` | P.grnas_per_cell | P.guide_counts_per_cell; retain the per-cell summary rather than recomputing a potentially costly aggregation for every plot. |
| `cells_w_zero_or_twoplus_grnas` | P.cells_w_zero_or_twoplus_grnas | P.assignment_exclusions; source-coordinate ambiguous/low-total assignment exclusions. |
| `cells_in_use` | P.cells_in_use | P.analysis_cells; unique ordered retained source-cell positions. |
| `n_discovery_pairs` | P.n_discovery_pairs | Removed counter. Derive total from canonical S.discovery_pairs via pair_count(). |
| `n_ok_discovery_pairs` | P.n_ok_discovery_pairs | Removed counter. Derive from pair-QC rows; Bonferroni counts eligible response-target families, not expanded guides. |
| `n_positive_control_pairs` | P.n_positive_control_pairs | Removed counter. Derive total from canonical S.positive_control_pairs. |
| `n_ok_positive_control_pairs` | P.n_ok_positive_control_pairs | Removed counter. Same passing-family policy as discovery. |
| `calibration_group_size` | S.calibration_group_size | Removed working slot; selected size is transient and frozen with the corresponding result. Old value only a historical effective snapshot. |
| `n_calibration_pairs` | S.n_calibration_pairs | Removed working slot; requested count belongs to the calibration plan, distinct from the number of available sampled pairs. Old value historical only. |
| `import_grna_assignment_info` | D.import_grna_assignment_info | Retained max_grna/max_grna_frac_umis import summaries aligned to source cells. |
| `mean_cells_per_grna` | P.mean_cells_per_grna | Removed redundant scalar. mean_assignment_size() derives mean lengths of durable guide_assignments. |
| `response_precomputations` | F.response_precomputations | New F.entries cache has explicit actual fitting identity. Old entries are discarded as unverified; only discarded_unverified_fit_count is recorded. |
| `last_function_called` | last_function_called | Removed call journal; not retained as live state or used to infer completion. |
| `functs_called` | functs_called | Removed Boolean authority. Used only as evidence when reconstructing preparation stages whose durable payloads validate. Never certifies old results/fits; raw flag vector is not retained. |
| `calibration_result` | R.calibration_result | R.records$run_calibration_check$table. Import as historical, preserve recognized history/provenance, invalidate currentness. |
| `power_result` | R.power_result | R.records$run_power_check$table; same historical-only migration. |
| `discovery_result` | R.discovery_result | R.records$run_discovery_analysis$table; same historical-only migration. |
| `nuclear` | nuclear | Removed pipeline mode; no corresponding working field. Not retained separately. |
| `nf_pipeline` | nf_pipeline | Removed pipeline/raw-result mode; public analyses finalize complete families. Old saved result scope may be retained as historical metadata only. |
| `integer_id` | D.integer_id | No screen-level ID slot. Backed descriptors own import IDs; an unanchored old ID can temporarily enter R.import_metadata$legacy_backend_id for reopening. It never establishes immutable content. |
| `elements_to_analyze` | elements_to_analyze | Removed pipeline subset scope. Current public result commits require their complete prepared family; no corresponding live slot. |

## All 14 additional schema 2 attributes

| Schema 2 path | Schema 3 decision |
|---|---|
| `data` | Retained coordinator reference, now strictly ScreenData. |
| `spec` | Retained coordinator reference, now resolved scientific choices only. |
| `preparation` | Retained coordinator reference with renamed domain payloads and rebuilt stage identities. |
| `fits` | Reference now typed ResponseFitCache, replacing FitRecord; old fits not reused. |
| `results` | Retained ResultRecord reference; its entries now own tables rather than referring to separate sibling slots. |
| `schema_version` | Set to 3L by explicit migration. Unknown versions/layouts error; no implicit slot-based upgrade. |
| `D.response_names` | Retained display annotations when saved. Schema 1 lacked them; leave unknown/empty instead of substituting canonical count IDs as biological labels. |
| `D.identity` | Recomputed content identity from recognized imported payloads/descriptors. Old identity is not adopted as a new cache certificate. |
| `D.backend` | Removed duplicate descriptor container. Each typed BackedCountStore owns its descriptor and process-local handle. Its response/grna descriptors are consumed during migration. |
| `S.requested` | Removed second settings authority. Saved original discovery/positive-control pairs recover canonical requests when available. Other working choices come from resolved S slots; old result provenance retains its own saved snapshot. |
| `P.grna_target_data_frame` | Removed duplicate full metadata table. Extract its derived grna_group vector into P.guide_groups; canonical metadata remains in D. |
| `P.stages` | Old complete/coordinates/cells/identity records are not trusted or copied. Rebuild current identity-only stages after payload validation and recognizable request recovery. |
| `F.identity` | Discard old fitting identity, including its incomplete design/version assumptions. New cache identity includes ordered actual fitting cells/design columns/values and method. |
| `R.records` | Migrate known analysis records/history into entries owning historical tables, saved provenance, legacy scope and an invalidation reason. Old currentness assertions do not survive. |

## Nested schema 2 metadata and containers

| Old nested field(s) | Retention or derivation decision |
|---|---|
| Count-store lists, element 1 | Unwrapped once; multiple-element wrappers error. Sparse normalized stores become direct dgRMatrix; ondisc objects/descriptors become BackedCountStore. Empty recognizable stores remain unavailable, not fabricated count data. |
| D.backend$response / $grna | Descriptor path, content_md5, dimensions, dimnames and integer_id move together into each typed store. A live handle is process-local; serialized output closes it. |
| D.import_grna_assignment_info$max_grna / $max_grna_frac_umis | Kept with their owning import summary list. |
| S.requested$discovery_pairs / $positive_control_pairs | Prefer saved originals, strip derived grna_group, normalize ID columns to character. Preserve repeated requests. |
| S.requested$side / $grna_integration_strategy / $formula_object / $resampling_approximation / $control_group / $resampling_mechanism / $multiple_testing_method / $multiple_testing_alpha | No independent working copies. Resolved S fields remain authoritative; retained old result provenance may contain these exact old snapshots. |
| S.grna_assignment_hyperparameters and S.cellwise_qc_thresholds entries | Retained method/QC settings in their existing lists; current signatures include their resolved values. See ARCHITECTURE.md for actual method-specific keys. |
| P.initial_grna_assignment_list$<guide ID> | Retain under guide_assignments with source coordinates. |
| P.grna_assignments_raw$grna_group_idxs / $indiv_nt_grna_idxs | Retain under source_groups. |
| P.grna_assignments$grna_group_idxs / $indiv_nt_grna_idxs / optional $all_nt_idxs | Retain under analysis_groups, preserving source-to-analysis and NT-list ordering. |
| P.stages$design / $assignment / $cell_qc / $pair_qc; their complete, coordinates, cells, identity members | Drop old completion/identity metadata and duplicate cells, then reconstruct validated current stage records. Current analysis_cells is the sole retained-cell authority. |
| F.response_precomputations$<response ID>, fitted_coefs/theta | Not reused or archived as fits. Count discarded entries in import_metadata. New runs create independently validated entries. |
| F.identity$data / $cells / $design / $method | Drop old cache certificate. Current identity additionally records explicit columns and actual selected design values. |
| R.records$<analysis>$complete / $status / $identity / $invalidation_reason | Drop old completion/currentness assertions. New migrated identity is NULL and invalidated_by is legacy provenance unknown. |
| R.records$<analysis>$scope | Keep old scope when supplied, otherwise unknown; migrated entries remain historical regardless. |
| R.records$<analysis>$provenance | Preserve under provenance$saved_provenance, including its schema_version, package_version, requested/effective/pairs snapshots and any contained assignment/QC/method details. Do not reinterpret those snapshots as verified effective history. |
| R.records$<analysis>$history, old history entry result/table | Recursively migrate recognized histories; preserve each historical table/provenance. Flattened current-entry history does not duplicate a nested history chain inside each old entry. |
| S.B1/B2/B3/calibration_group_size/n_calibration_pairs/negative_control_pairs at import | Preserve recognizable values as each migrated result's provenance$effective_snapshot, explicitly unknown as evidence of that result's actual original run. |
| Top-level and component class attributes | Read physical attributes and construct current registered classes; do not transplant obsolete class layouts. |
| Unrecognized attributes or nonanalysis R.records entries (including an older records$legacy catch-all) | No generic lossless archive is implemented. These are not copied into a second live schema. Keep the original RDS if such implementation-specific metadata must be retained. |

## Migration guarantees and limits

Recognized data, preparation values, result tables and saved result histories
remain readable. Old result records are explicitly historical with unknown
provenance, even if old flags/status said current. Old response fits are
intentionally discarded rather than blessed by an incomplete identity.

For union, empty or already unexpanded requests, migration reconstructs canonical
requests directly. Singleton/Bonferroni expansion recovery walks ordered metadata
blocks and preserves duplicate requests; it does not simply unique rows. Recovery
certainty is saved per family as saved_original, reconstructed or unknown. Unknown
recovery leaves tables readable but establishes no current preparation stages.
When recognizable, saved flags and validated design/assignment/QC payloads can
establish fresh stage identities; they do not resurrect results.

Schema 1 has no saved response display-label authority. Mitochondrial summaries
already in covariates are retained, but missing original labels are not inferred
from canonical feature IDs. Removed logging/pipeline fields and unrecognized
attributes are deliberately outside the retention guarantee. No migration claims
Nextflow compatibility, historical RNG replay or verified old per-run budgets.

New schema 3 additions have one authority: S.additional_cells_to_remove;
P.guide_groups and canonical-to-compiled pair tables; typed backing descriptors;
ResponseFitCache entries/identity; result entries/provenance/history and
import_metadata. Transient plans/CellViews are not separate persisted working
fields. See OBJECT_FIELDS.md for the exhaustive current slot/type listing.

## Earlier schema 3 run-plan values

Earlier development snapshots stored each result's realized plan as a plain
list. The current owner is the S4 `RunPlan` value within that same result's
provenance. Its scalar fields, selected family and original frozen identity move
together; the canonical payload keeps the old field order and explicitly named
`NULL` calibration fields for hash compatibility. History entries are converted
at the same boundary. Unknown fields or malformed values are retained as
unrecognized historical provenance, with a diagnostic, rather than silently
dropped or re-signed as a valid plan. No legacy list-field forwarding layer or
`RunPlan` dollar-access generic is installed.
