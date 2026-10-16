# Private-symbol dispositions and caller migration

Planning amendment dated 2026-10-16, based on inspected production source
`c7d62502757f2a17741e2ef18b9e4d006bf063dc`. These are concrete cleanup targets,
not implemented changes or qualification claims. Public exports/formals/methods,
56-slot saved objects, ordinary-case science/RNG and real consumer contracts
remain protected. **No compatibility aliases are proposed.** A removed private
name is removed from source, generated bindings and current executable tests;
historical source/oracles/reports keep their original names and hashes.

Scope: bounded reread of the source/callers below, the prior core/execution/native
inventories, `tests/testthat/helper-refactor-native.R`, native fixture inputs and
bridge, and the 11 pinned pipeline-script rows. No test or workflow was run for
this amendment. No claim covers arbitrary unpublished private-namespace users.
Pipeline source is pinned at `38dc760a5b4105d1c3f339cd45f2c98d971087f4`; book
contracts at `63704b0cca4b334f5b957da928ce718d38184f20`. Source inspection and
historical executed qualification remain distinct.

## Resolved conflicting proposals

- Retain the useful per-run treatment lookup closure, renamed
  `prepared_index_factory` → `make_treatment_index_lookup`. This supersedes core
  review C05's factory-deletion alternative. Its implementation uses the one
  canonical index; workers do not reconstruct indices or cache all guide views.
- Adopt preserved matrix dimensions (`drop=FALSE`) for already accepted
  intercept-only designs **before or atomically with CL3**. This supersedes the
  original execution E7 proposal to reproduce the NT vector/error. Saved
  full-source design/cache rules remain unchanged; independent valid NT
  calibration/discovery regressions cover both traversal modes. Existing audit
  evidence honestly records the former error; ordinary fixtures stay frozen.
- Retire production `crt_index_sampler`, the uncalled slow CRT implementation.
  Keep its immutable original source/config/full-bank archive as historical
  evidence. It has no candidate substitute. Active comparisons are enumerated
  explicitly below; function existence never decides coverage.
- The held checkpoint `8f1b8f4` named its replacement pair kernel
  `compute_pair_counts_sparse`. That provisional name is superseded by
  `sparse_pair_counts`; do not integrate or retain a second alias under it.

## R aliases, context and numerical components

| Current private symbol / owner | Final disposition and owner | Actual caller/test migration |
| --- | --- | --- |
| `perform_response_precomputation`, `R/precomputation_functions.R:3` | Delete forwarder; retained numerical destination `fit_response` in `R/fitting.R`. Delete precomputation_functions.R after all three aliases are removed. | No runtime caller found; developer source-tracing prose names the actual fitter. |
| `perform_grna_precomputation`, same file :7 | Delete forwarder; destination `fit_assignment_probabilities`, `R/fitting.R`. | No runtime caller. Do not retain coefficient-return mode solely for an alias/test. |
| `compute_precomputation_pieces`, same file :11 | Delete forwarder; destination `prepare_score_components`, `R/score_preparation.R`, consuming the actual compact fit record. | Retarget `test-refactor-indexing.R:61` from historical coefficient/theta arguments to real component inputs; retain independent observed-versus-synthetic coordinate fixture. |
| `fit_response_current`, `R/fitting.R:16` | Rename `fit_response`, same owner; arithmetic unchanged. | The two group-fit workers, shared-view driver cache misses and component tests. Remove alias caller above. |
| `fit_assignment_current`, `R/fitting.R:54` | Rename `fit_assignment_probabilities`, same owner; remove unused `return_fitted_values` flag/coefficient-only mode. | The two CRT workers currently request TRUE (:135,:216). Update focused fitter tests to retained probability output; no replacement compatibility mode. |
| `prepare_score_current`, `R/score_preparation.R:34` | Rename `prepare_score_components`, same owner; remove unused reduced-mode flag and retain full current score record. | The group-fit permutation worker, both CRT workers and shared permutation preparation in its driver; update direct tests. No new score registry. |
| `construct_data_frame_v2`, `R/result_processing.R:163` | Delete alias; retained `pair_results_from_native`, same file :3. | No runtime caller found; drivers already use actual core. |
| `get_idx_vector_calibration_check`, `R/analysis_context.R:154` | Delete alias; retained `calibration_group_membership`, same file :32. | Retarget `test-refactor-indexing.R:85`; test supported low-MOI partitions and high-MOI complement unions separately. |
| `legacy_nt_fit_assignments`, `R/analysis_context.R:20` | Rename `encode_nt_fit_assignments`, same owner; preserve valid saved NT-local schema and ordered pool, using explicit membership-to-pool positions. | `update_grna_assignments_given_qc` calls it directly. Supported saved partition/getter tests replace malformed overlapping cumulative positions. |
| `update_indiv_grna_assignments_for_nt_cells`, `R/cellwise_qc_functs.R:124` | Delete one-call forwarding alias. | Its actual caller :111 calls `encode_nt_fit_assignments` directly; indexing test :91 migrates to actual boundary. |
| `prepared_index_factory`, `R/analysis_context.R:113` | Rename `make_treatment_index_lookup`, same owner; retain lean closure. | Both renamed drivers prepare it once. The shared-fit permutation worker receives the accessor; the CRT guide block resolves plain treatment inputs before calling its shared-fit worker. C05 deletion superseded. |
| `native_observed_positions`, `R/analysis_context.R:151` | Delete integer passthrough. | Four workers pass integer R one-based treatment positions directly to the actual Rcpp boundary, which alone subtracts one. Keep meaningful native base fixture; delete identity-wrapper-only assertion. |
| `legacy_fit_source_cells`, `R/analysis_context.R:165` | Delete optional direct-worker fallback. | Remove its sole call, `crt_glm_factored_out:146`, atomically with explicit source-cell/row worker inputs. |
| `prepare_cell_map`, `source_membership_to_analysis`, `calibration_group_membership`, `resolve_pair_rows`, `prepare_index_context`, `group_treatment_analysis`, `make_fit_view`, `prepare_analysis_context`, `R/analysis_context.R` | Retain names and coordinate ownership. Context owns one transient common view/design; saved design stays full-source. | Actual QC/context/worker callers remain. Normalization happens once; no sparse NT round trip or all-guide view bank. |
| `integrate_guides_current`, `R/result_processing.R:64` | Rename `integrate_guides`, same owner. | Actual `apply_grouping_to_result` and `finalize_association_result` calls migrate; three real facade result adapters below remain. |
| `set_colnames`, `R/test_functs.R:1` | Remove from package; use ordinary `colnames(x) <- NULL` directly in test support. | Sole inspected executable use `tests/testthat/helper-mock-data.R:203`. Do not regenerate frozen inputs. |
| `set_rownames`, `R/test_functs.R:6` | Delete unused function; remove empty test_functs.R. | No executable package/test/pinned-consumer caller found. |

## Drivers, workers and explicit effects

| Current private symbol / owner | Final symbol / owner | Caller migration and retained role |
| --- | --- | --- |
| `run_sceptre_analysis_high_level`, `R/analysis.R:3` | `analyze_association_pairs`, same file | Three calls `R/s4_analysis_functs_2.R:175,278,381`; remove redundant calibration_check argument. S4 public functions/formals stay unchanged. Boundary creates global permutation bank, performs initial GC, dispatches and assembles family. |
| `run_perm_test_in_memory`, `R/association_execution.R:43` | `run_permutation_analysis`, same file | Direct call from renamed boundary, replacing string do.call. `(context, execution, synthetic_idxs)`; delete nullable context reconstruction, raw duplicate fields and unread thresholds. Preserve response traversal and shared fits/scores. |
| `run_crt_in_memory_v2`, same file :223 | `run_crt_analysis`, same file | Direct call `(context, execution)`; preserve response prepass, then guide traversal, full bank per guide and cache merging. |
| Local `analyze_given_response_ids`, same file :98 | `analyze_response_block`, same owner | Permutation block executor callback only; retains request order, response reads and miss-only cache collection. |
| Local `run_precomp_on_given_responses`, same file :277 | `prepare_response_block`, same owner | CRT prepass callback; preserve cache-hit response read and log creation lifecycle. |
| Local `perform_association_analysis`, same file :329 | `analyze_guide_block`, same owner | CRT guide callback; preserves file append, guide order and response-row resolution. |
| `perm_test_glm_factored_out`, `R/association_workers.R:2` | `test_permutation_shared_fit`, same file | Permutation driver :146; explicit response/score, groups, treatment lookup, bank and options. |
| `discovery_ntcells_perm_test`, same file :43 | `test_permutation_group_fit`, same file | Driver :160; explicit prepared index and common-view arrays; active group view and response fit remain per group. Delete index=NULL reconstruction/raw NT inputs. |
| `crt_glm_factored_out`, same file :111 | `test_crt_shared_fit`, same file | CRT driver :356; explicit store row indices/source cells/treatment/cache. One assignment fit/full bank reused for every response. Delete optional rematching/source-cell fallback. |
| `discovery_ntcells_crt`, same file :184 | `test_crt_group_fit`, same file | CRT driver :376; explicit rows and active guide fit view; one assignment fit/bank followed by guide-specific response fits. Delete partial index reconstruction. |
| `get_id_from_idx`, `R/association_execution.R:2` | Delete combined helper; `report_worker_progress` and `collect_worker_garbage` owned by `R/execution.R` | Three association loops (:108,286,338) and mixture :69 perform plain ID lookup, report, then explicit maintenance. Preserve message/file bytes, first/every-five reports, every-200 GC even with progress off. |
| `analysis_log_file`, `announce_parallel_analysis`, `execute_analysis_blocks`, `bind_analysis_results`, `R/execution.R` | Retain focused helpers and names | Real log create/append, exact announcements, serial/fork block order/core count and rbindlist(fill=TRUE) semantics. Not generic extension wrappers. |
| `assign_grnas_to_cells_mixture`, `R/mixture_model_functs.R:1` | Retain name/owner | Actual assignment dispatch :27 and public mixture workflows. Local `analyze_given_grna_ids` → `assign_guide_block` for named responsibility; preserve names/order/fork behavior. |
| `obtain_em_assignments`, same file :129 | Retain name/owner; remove unused `use_glm` flag/false mode | Sole actual call :79 always uses TRUE (:84). Keep Poisson initialization, warning policy, EM, active low-count/nonconvergence backup threshold. |
| `get_random_starting_guesses`, same file :187 | Retain name/owner | Mixture :17 and aux test :160; fixed with_seed(4), pi draws before perturbation draws, caller RNG restored. |

Current executable selectors in `tools/refactor/measure.py`, current walkthrough
prose and direct component test lookups migrate with these names. Historical
measurement JSON and pinned source/oracle names do not change. No test-only alias
is retained to make a namespace lookup succeed.

## Native boundary, canonical counts and scratch names

For every renamed Rcpp export, regenerate **both** `R/RcppExports.R` and
`src/RcppExports.cpp` from source. Its generated R function and registered
`_sceptre_<name>` destination use the final name; old registration disappears.
Do not hand-maintain either file or add an old-name forwarding registration.

| Current native/private symbol / owner | Final symbol / owner | Actual migration / coordinate and ownership reason |
| --- | --- | --- |
| `run_low_level_test_full_v4`, `src/low_level_full_test.cpp:27` | `run_association_test`, same file | Four R workers; test-refactor-indexing.R:68; helper-refactor-native.R:46. Retain actual Rcpp record/pointer boundary: observed R one-based positions converted once; bank entries already zero-based; same result fields. |
| `compute_observed_full_statistic_v2`, same file :9 | Delete unused global wrapper | Actual core `sceptre::inference::observed_score`, inference_scores.cpp/inference.h, remains. No executable runtime consumer found. |
| `compute_null_full_statistics`, same file :15 | Delete unused global wrapper | Actual core `sceptre::inference::null_score_segment` remains; no pointer wrapper required outside actual R entry. |
| `estimate_log_fold_change_v2`, `src/inference_scores.cpp:98` | Delete unused global wrapper and declaration in shared header | Actual core `sceptre::inference::estimate_effect` remains. No direct runtime consumer found. |
| `sceptre::inference::legacy_result_record`, inference_stages.cpp:51 / inference.h:57 | Rename `sceptre::inference::association_result_record`, same owners | Actual Rcpp entry and `tools/refactor/native-components.cpp:49` migrate. Record preserves real result schema; name does not imply optional old algorithm. |
| `observed_positions_zero_based`, `estimate_effect`, `observed_score`, `null_score_segment`, `evaluate_staged_tail`, inference owners | Retain names and inference ownership | Actual entry/stages and native-component harness consume them; preserve const full bank, selected stage segments, sort/fallback arithmetic. |
| `compute_nt_nonzero_matrix_and_n_ok_pairs_v3`, `src/pairwise_qc.cpp:10` | Replace with `sparse_pair_counts`, new owner `src/sparse_counts.cpp` | `R/storage.R:134` migrates. Prepared treatment/NT groups use analysis positions; owned C++ vectors convert one-based positions once. Preserve response-sorted fused scan and adapter request restoration, stored zeros and NT-group-sum totals. Delete superseded pairwise_qc.cpp implementation. Held `compute_pair_counts_sparse` also retires, with no alias. |
| `compute_n_trt_cells_matrix`, `src/negative_control_functions.cpp:155` | `sparse_group_counts`, `src/sparse_counts.cpp` | `R/storage.R:173` migrates. Move count ownership out of calibration sampling; same all-response group matrix/selection contract; nonmutating owned positions. |
| `load_nonzero_posits`, shared_low_level_functions.cpp/.h | `load_stored_positions`, private to sparse_counts.cpp (anonymous namespace) | The two canonical count kernels are its only production callers. Not an Rcpp export. Stored positions include explicit zeros; direct source/analysis coordinate contract. Delete obsolete shared implementation/header after tail declarations relocate. |
| `compute_colwise_max`, `src/sparse_matrix_functs.cpp:92` | Delete unused unexported helper | No runtime/test/pinned consumer found; do not invent maximum-assignment use. |
| `synth_idx_list_to_matrix`, `synth_idx_list_to_r_list`, `print_synth_idx_list_row`, `test`, `src/check_functions.cpp` | Delete all four and empty file | No runtime callers; test-only native bridge materializes banks independently. No generated exports exist for these four. |
| `increment_vector`, `src/sparse_matrix_functs.cpp:130` | Delete after actual callers migrate to owned R arithmetic. | B1 removed the importer calls. The sole remaining assignment getter now uses `unlist(...) - 1L` to create its owned zero-based Matrix slot. Independent single/multiple/empty-guide, repeated-call, slot and serialized-caller tests preserve the getter's existing QC dimension drop. This supersedes the proposed private `offset_indices_in_place` rename: no production caller needs a native vector mutator. Native retirement removes its export/registration without adding an alias. |
| `increment_matrix`, `src/negative_control_functions.cpp:87` | `increment_group_indices_in_place`, same file | Sole R caller neg_control_functions.R:98 migrates. Mutates scratch combination matrix from zero-based to one-based after native sampling, before label construction; no new generic matrix operation. |
| `convert_pointer_to_index_vector_v2`, `R/storage.R:354` | `sparse_pointer_to_indices`, same file | Two real storage calls :306,:308 migrate, including held 13ff940 naming proposal. No old-name alias. |
| `load_csr_row`, `load_sparse_vector`, `obtain_pointer_vector`, sparse_matrix_functs.cpp | Retain names/owner | storage read and conversion; importer CSC construction; cell covariate loading. Real representation primitives, not duplicate counting semantics. |
| `compute_cell_covariates_cpp`, sparse_matrix_functs.cpp:40 | Retain name/owner, use owned mitochondrial coordinate conversion | Actual storage_cell_covariates caller. Preserve loaded-value >0.5 counting; not interchangeable with stored-position pair counts. Remove R defensive copies only when own kernel stops mutating input. |
| `compute_n_grnas_per_cell_vector`, `threshold_count_matrix`, sparse_matrix_functs.cpp | Retain names/owner | Actual getter/threshold storage callers and public workflows; preserve count/coercion semantics. |
| `sample_combinations_v2`, negative_control_functions.cpp:30 | Rename `sample_calibration_groups`, same owner | construct_negative_control_pairs_v2 :41 migrates; same fixed stream/sampling order with separately reviewed safety prerequisites. |
| `sample_undercover_pairs_v2`, same file :100 | Rename `sample_calibration_pairs`, same owner | Actual R calls :56,:82 and test-refactor-qc.R:20 migrate. Shared eligibility predicate retained; no new rate floor/heuristic. |
| `iterate_over_combinations`, same file :66 | Retain name/owner | Actual calibration group construction :29; post-QC active-group bounds checked before native write. |
| `construct_negative_control_pairs_v2`, `R/neg_control_functions.R:1` | Rename `construct_negative_control_pairs`, same owner | Actual run_calibration_check_pt_1 call R/s4_analysis_functs_2.R:151 migrates; pipeline calls the retained facade part, not this helper. |
| `compute_empirical_p_value`, `fit_and_evaluate_skew_normal` | Retain implementations in inference_tail.cpp; declarations owned by inference.h | Stages call them; remove their declarations from obsolete shared_low_level_functions.h. Keep actual cross-file boundary. |
| `fit_skew_normal_funct`, `check_sn_tail`, `check_for_outliers` | Retain as tail-local helpers in anonymous namespace of inference_tail.cpp | Only tail implementation calls them. Remove unnecessary shared declarations; do not alter tail tests/fallbacks under an ownership rename. |

`legacy_pair_eligible` / `legacy_pair_eligibility_value` in pair_eligibility.h are
retained: “legacy” identifies the deliberately fixed statistical policy, unlike
an obsolete version suffix. `estimate_theta` and native EM entrypoints retain
their active algorithms; no new numerical fallback or unreviewed exception
policy is implied. R storage count/read/backed adapters retain names and own
backend selection. `storage_owned_groups` stays **only at external ondisc calls**
that still mutate vectors; canonical owned sparse kernels do not need R copies.

## Explicit sampler mapping and independent evidence

The immutable input contains exactly four named historical configurations. The
mapping below is authoritative, not a `names()`/`exists()`-selected subset.

| Frozen configuration | Frozen original function and arguments | Candidate production function | Required comparison |
| --- | --- | --- | --- |
| `fisher_yates` | `fisher_yates_sampler(n_tot=13L,M=4L,B=17L)` | `fisher_yates_sampler`, generate_samples_functions.cpp:30 | Candidate **full 17-draw bank** versus independent original full bank; retain range/membership/determinism/RNG laws. |
| `hybrid` | `hybrid_fisher_iwor_sampler(N=9L,m=2L,M=5L,B=17L)` | `hybrid_fisher_iwor_sampler`, same file :66 | Candidate **full 17-draw bank** versus independent original full bank; retain every requested prefix law. |
| `crt_fast_production` | `crt_index_sampler_fast(fitted_probabilities=c(0,.15,.4,.85,1),B=5498L)` | Rename to `sample_crt_indices`, same file :159 | Candidate **full 5,498-draw bank**, not prefix, versus independent original fast sampler full bank. Both CRT workers migrate. |
| `crt_legacy` | `crt_index_sampler(fitted_probabilities=c(0,.15,.4,.85,1),B=17L)` | **None: delete production slow sampler, generated R wrapper and registration.** | Frozen original slow oracle full bank versus unchanged archived `native-expected.rds$sampler_draws$crt_legacy` on its recorded platform/toolchain; source/config hashes checked on every platform. This is archive integrity, not candidate parity. |

The active candidate map is a literal three-entry mapping:

```r
c(fisher_yates = "fisher_yates_sampler",
  hybrid = "hybrid_fisher_iwor_sampler",
  crt_fast_production = "sample_crt_indices")
```

Assert the frozen input has exactly the four named configurations above and
that the active map has exactly the three retained keys. Missing candidate
symbols are hard failures. No `exists` filtering, namespace-intersection
selection or fallback to a reference function is allowed. Frozen function names
remain unchanged in inputs; the map selects candidate names separately.

For each active key, actual is produced by the **installed candidate namespace**
using its mapped name and the unchanged arguments, then fully materialized by
the test-only bridge. Expected is independently produced by
`refactor_native_reference_sample(config$function_name, config$args)`, compiled
from the frozen original include. Never send candidate expected/actual both
through the frozen reference. Candidate-vs-stored-bank checks on the recorded
platform also remain for the three active keys. The retired slow key runs only
the explicitly named oracle-versus-archive check, not production laws/parity.
Do not silently substitute fast CRT for it or relabel archival coverage as an
active fourth sampler.

Immutable authority: original sceptre
`21f9ea098b69c444a884a49bc254de711968e3b5`, original source
`src/generate_samples_functions.cpp`, Git blob
`5a839ba737db9d6ef0e5e031adb7486d3692ff91`, SHA256
`deb2dee4338f5779588408f4d80caa9a42d64e1df5e2a5c519c8361dd1c06cfd`.
`tests/testthat/fixtures/refactor/native-sampler-reference.inc`, original fixture
inputs/expected/manifest and historical reports are not regenerated. The bridge
keeps all four historical oracle dispatch cases. Its source belongs to test
support and never includes candidate implementation source. This resolves the
slow-sampler decision without weakening the original independent full-bank gate.

## Real adapters and consumers retained by name

| Symbol / owner | Retained consumer and contract |
| --- | --- |
| `process_calibration_result`, result_processing.R:145 | Package calibration part; pipeline bin/process_association_analysis_results.R:34. Core standard complete-family finalization remains behind this facade adapter. |
| `apply_grouping_to_result`, same file :157 | Package power check; pipeline process_association_analysis_results.R:35. Real guide-integration adapter, not obsolete alias. |
| `process_discovery_result`, same file :151 | Package discovery; pipeline process_association_analysis_results.R:36 and run_association_analysis_nuclear.R:66. Preserve standard versus nuclear raw/pod behavior; no new global nuclear correction. |
| `run_calibration_check_pt_1`, s4_analysis_functs_2.R:89 | Pipeline prepare_association_analyses.R:27; creates prepared calibration pair family before pod work. |
| `run_calibration_check_pt_2`, same file :165 | Pipeline run_association_analysis.R:27. Preserve nf_pipeline output boundary, result.rds and precomputations.rds contract. |
| `run_qc_pt_1`, `run_qc_pt_2`, s4_analysis_functs_1.R | Pipeline run_qc_nuclear.R:33 and run_association_analysis_nuclear.R:52. Keep separated nuclear preparation/pod QC phases. |
| `process_initial_assignment_list`, assign_grna_functions.R:136 | Pipeline process_grna_assignments.R:40 and package assignment path; retain actual post-pod assignment processing. |
| `determine_grnas_in_use`, same file :203 | Package assignment dispatch and pipeline output_grna_info.R:16 passes **trial as restricted_grnas**. Retain FALSE/default and TRUE modes, standard target restriction, nuclear fixed-seed sample of at most 30 targeting metadata rows plus NT, and current order. The test default alone is not coverage of this real pipeline mode. |
| `get_id_vect`, ondisc_functs.R:99 | output_grna_info.R:18 and prepare_association_analyses.R:39. Keep actual pod assignment order/IDs and current rounding semantics; do not silently replace with a new chunking rule. |
| `write_vector`, same file :109 | output_grna_info.R:22,24; prepare_association_analyses.R:57,58; prepare_association_analysis_nuclear.R:28. Retain text filenames/content contract; connection cleanup is a separately bounded fix, not a helper deletion. |
| `partition_response_ids`, aux_functions.R:87; `get_log_dir`, :136 | Package association/mixture block scheduling and log paths. Retain ordering, minimum-block behavior and canonical path handling; no new scheduler. |
| `load_row`, s4_helpers.R:363 | Mixture/plot/current callers; real ID-to-storage boundary; retain and route through storage. |
| `storage_create_backed_pair`, `storage_reopen_backed_pair`, storage.R | Public backed import/reopen paths, original .odm IDs and process-local handle reattachment; needed ondisc private-constructor boundary. |
| `read_ondisc_backed_sceptre_object`, `write_ondisc_backed_sceptre_object`, ondisc_functs.R | Public API retained. Pinned pipeline scripts actually reopen saved objects; no serialized pointer redesign. |

Exported `run_power_check`, `run_discovery_analysis`, `assign_grnas`, setters,
getters, plots and public imports retain names/formals even where a pipeline
script uses `:::` spelling. The 11 script hashes/arguments/direct slots/fixed
filenames remain in consumer-inventory.json. Standard worker result.rds /
precomputations.rds, standard gather's complete-family outputs and nuclear
result_<pod>.parquet/.rds are separate genuine contracts. No rename here
requires a pipeline source edit. Full paired script/affected book qualification
follows implementation; this document does not claim those runs passed anew.

## Atomic migration/readback gates

- Land source definitions, actual callers, direct tests/tool lookups and generated
  bindings coherently. Search active R/src/test/tool execution for each removed
  name; distinguish intentional historical fixture/provenance/report mentions.
- All three active sampler full-bank comparisons remain candidate versus the
  independent original. The retired fourth configuration has the explicitly
  separate oracle/archive integrity check. Missing symbols or mappings fail.
- Canonical pair/group kernel, loader ownership and storage callers land as one
  unit; no historical v3 or held-checkpoint alias. Preserve ondisc adaptation
  only where its actual native contract requires it.
- Retain exact public API/S4/saved-object, supported indexing/ordering/cache/RNG,
  real consumer and resource gates. Newly adopted intercept-only behavior uses
  independent intended calculations, without refreshing ordinary goldens.
- Recheck final source/deletion crosswalk against the final immutable build;
  doc adoption or old 28df390/90664a7 evidence is not final qualification.
