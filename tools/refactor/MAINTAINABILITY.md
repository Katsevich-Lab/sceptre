# Maintainability scorecard

Freeze this baseline before RF1. Compatibility remains a hard gate. This
scorecard supplies several independent engineering observations; it does not
combine them into one score or demand a decrease in every number at every stage.
Production R and handwritten C++ are measured separately, excluding generated
Rcpp wrappers, tests and the frozen test oracle. Analyzer versions and rules are
pinned in the result files.

## Change locality and rule ownership

These are fixed source-tracing exercises, not timed human experiments or claims
that every listed function must be edited. Repeat the same questions after the
relevant boundaries stabilize and at RF9. Record functions/files actually read
and edited, participant familiarity/tooling, ambiguities and the resulting diff.
Do not infer comprehension time from static counts.

| Fixed exercise | Baseline source-tracing set | Intended stage observation |
| --- | --- | --- |
| Locate and change the legacy pairwise eligibility rule consistently with calibration | `compute_qc_metrics` in `R/pairwise_qc_functs.R`; `construct_negative_control_pairs_v2` in `R/neg_control_functions.R`; `sample_undercover_pairs_v2` in `src/negative_control_functions.cpp` | One owner for the rule; calibration and discovery use the same implementation |
| Trace all response fitting paths before adding an output-compatible batched fitter | `perform_response_precomputation` in `R/precomputation_functions.R`; `discovery_ntcells_perm_test` and `discovery_ntcells_crt` in `R/association_testing_functions_v2.R`; `run_perm_test_in_memory` and `run_crt_in_memory_v2` in `R/medium_level_functs_v2.R` | Fitter/score preparation boundary and explicit fit views; retain NT-specific fits and traversal reuse |
| Trace the response-row and grouped-nonzero operations before changing their storage implementation | `load_row` in `R/s4_helpers.R`; `compute_pairwise_qc_information` in `R/pairwise_qc_functs.R`; `construct_negative_control_pairs_v2` in `R/neg_control_functions.R`; direct CSR reads in both CRT workhorses and CRT driver | Statistical consumers use adapter operations; sparse and backed kernels remain legitimate distinct implementations |
| Determine whether a response fit remains valid after formula, assignment or QC-cell changes | `set_analysis_parameters`, `assign_grnas`, `run_qc_pt_1` in `R/s4_analysis_functs_1.R`; response-cache collection in both drivers; three facade writebacks in `R/s4_analysis_functs_2.R` | Central ownership of invalidation/merging and a short explanation of each fit's actual view |

The reviewed baseline has two independently maintained eligibility predicates:
`R/pairwise_qc_functs.R:133–134` and
`src/negative_control_functions.cpp:137`. Response-cache clearing has three
event-specific sites (`R/s4_analysis_functs_1.R:180,322,501`). Those are separate
responsibilities despite sharing a file. Count future ownership by responsibility,
not arbitrary filenames. Cache writeback/merging is distinct from invalidation.

Raw-storage and coordinate-candidate site inventories are diagnostics. They need
reviewed classification before being called architectural violations: imports,
adapter implementations and native kernels legitimately touch representations;
statistical orchestration should not. Distinct sparse/HDF5 kernels and response/
guide traversal plans are not duplicate scientific logic merely because they
look similar.

The reproducible initial R-only site scan records six private ondisc calls in
four files, 51 lines touching raw representation slots in nine files, 17 narrow
view-remapping candidate lines in four files, three reset events and three
facade cache writebacks. Each matching line is retained in
`results/architecture-sites.json`; these counts are not all boundary violations
or a complete inventory of coordinate conversions. The baseline and candidate
site reports are identical before structural work.

## Fixed independent-testability inventory

Keep this denominator fixed. Do not improve it by adding trivial helpers.

| Core operation and current entry | Needs an entire `sceptre_object`? | Other preparation burden |
| --- | --- | --- |
| Response fit: `perform_response_precomputation` | No | Ordinary expression vector and design matrix |
| Assignment fit: `perform_grna_precomputation` | No | Ordinary membership indices and design matrix |
| Score pieces: `compute_precomputation_pieces` | No | Ordinary vectors/matrices and compact fitted coefficients |
| Legacy eligibility: `compute_qc_metrics` | Yes | Pair tables, thresholds, grouping strategy and object mutation |
| Staged inference: `run_low_level_test_full_v4` | No | Ordinary score/effect inputs **plus a native resampling pointer**; the RF0 bridge does not change its production interface |
| Final discovery processing: `process_discovery_result` | Yes | Grouping, family and correction settings taken from the object |
| Response row: `load_row` | No | Explicit supported matrix representation and row ID |

Thus five of seven entries can be called without the entire facade, and four of
seven accept only ordinary explicit inputs in the sense above. The inference
pointer burden remains visible. This is a scoped entrypoint inventory, not a
claim that these operations have complete independent production test coverage.

## Complexity, coupling and duplication

Use per-function cyclomatic complexity, nesting/length and named hot spots;
report maxima and upper percentiles rather than relying on package averages.
R cyclocomp and C++ Lizard results have distinct analyzer definitions and are not
directly comparable or additive. Cognitive complexity is unmeasured until a
trustworthy pinned language analyzer is available. No percentage target is set
before inspecting the baseline.

The first measurements use cyclocomp **1.1.2** for R and Lizard **1.24.1** with
its `nd` extension for C++. There are 184 R function literals measured as
exclusive units (nested function bodies measured separately): max 49, p90 9.70,
p95 13.85. The 122 named top-level owners have inclusive p95 22.75. Handwritten
C++ has 43 functions: max 17, p90 8.80, p95 13.60. Maximum nesting is 7 under
the documented R AST rule and 5 under Lizard's different C++ rule.

For the named drivers, exclusive/inclusive R complexity is 7/15 for
`run_perm_test_in_memory` and 9/23 for `run_crt_in_memory_v2`; inspecting only
the outer functions would hide their worker bodies. Native staged inference
has Lizard complexity 9 and skew-normal fitting/evaluation 14. Detailed
per-function/source hashes and summaries are in `results/maintainability-baseline/`.
Baseline and candidate detailed metrics match exactly before structural work.

Reproduce with analyzer dependencies in separate development-only libraries:

```sh
R_LIBS=/path/to/cyclocomp-library PYTHONPATH=/path/to/lizard-library \
  python3 tools/refactor/measure.py /path/to/source /path/to/output \
  --source-commit recorded-commit
python3 tools/refactor/measure-sites.py /path/to/source /path/to/site-inventory.json
```

The analyzer's source-commit argument is only a label. Per-file and aggregate
source hashes establish actual content identity. The Python driver resolves
`measure-r.R` beside itself. Neither script evaluates the measured package.

Static direct-call coupling, when recorded, excludes dynamic calls and S4
dispatch and therefore is not a complete runtime dependency graph. Module edge
counts can grow with useful decomposition; inspect cycles, forbidden crossings
and unusually broad responsibilities. Reviewed duplicate-rule ownership above
is measured separately from token/block clone detection, which is not yet
measured. No test-coverage percentage or LOC reduction proves readability.

At RF1 inspect remapping sites and the effort to explain source/analysis/NT/fit
views. At RF2–3 inspect storage leaks and locality; RF4 eligibility ownership;
RF5 fitting/scoring testability and cache validity; RF6 duplicated orchestration
and driver nesting; RF7 native score/effect/inference boundaries; RF8 component
coupling and wrapper burden. At RF9 repeat the entire scorecard and the fixed
human exercises, including a responsibility/dependency map and deleted paths.
Cheap automated metrics can be refreshed at every reviewed commit; human
exercises belong at RF0, stable major boundaries and RF9. Human timing has not
been collected for this baseline.

Each completed stage should improve its targeted responsibility. An adapter can
add calls; temporary compatibility wrappers can add lines; splitting functions
can raise summed cyclomatic complexity because each function contributes one.
Record any worsening with its benefit and, if temporary, removal stage. Final
acceptance requires overall architectural improvement and retained behavior.

## RF9 source walkthrough

The coordinator repeated the four fixed source-tracing questions on the
integrated refactor, checked actual callers and removed routes, and had an
independent senior engineer review the source. This is a source/design exercise,
not a timed comprehension study or an implementation of any deferred method.
No human timing, cognitive-complexity or clone-coverage result is claimed.

| Fixed exercise | Refactored ownership and concrete change path | Removed coordination burden |
| --- | --- | --- |
| Change the legacy eligibility rule consistently with calibration | `src/pair_eligibility.h` owns both valid-count policy and missing-value compatibility; its R vector boundary and the unchanged candidate sampler call it. Measurements live separately in `R/pair_qc.R`/storage. | Discovery and calibration no longer maintain independent threshold predicates. A selected future policy needing new inputs still requires explicitly defining and passing them. |
| Trace response fitting before adding an output-compatible fitter | `R/fitting.R` owns response/assignment fits; `R/score_preparation.R` owns derived arrays. Workers consume explicit fit views and call these modules. | The numerical implementation is independent of storage, S4 state, scheduling and caches. NT fits retain their correct guide-specific view. No speculative batch/fitter registry was added. |
| Change response-row/grouped-nonzero storage | `R/storage.R` owns reads, thresholds, pair/group kernels, backing creation and reopening. Actual assignment/QC/calibration/CRT consumers call it. | Three direct CRT CSR reads and backend branches in statistical consumers are gone. Sparse/HDF5 kernels remain legitimate separate implementations. |
| Decide response-cache validity after supported mutations | `R/fit_cache.R` owns invalidation tests, resets, lookups/merges and facade writeback; setters provide the actual change event. | Three reset implementations and three direct writebacks become calls to one owner. Control-view and actual-design/contrast changes now conservatively clear incompatible fits. |

The seven-operation denominator remains fixed. Current fit, assignment, score,
eligibility, staged inference, finalization and row-read operations can all be
tested without an entire `sceptre_object`. The latter three interfaces that
previously imposed additional burdens are now `evaluate_pair_eligibility()`
(plain aligned measurements/options), native `evaluate_staged_tail()` (explicit
arrays and a const ordinary draw bank) and `finalize_association_result()`
(tables and settings). The legacy Rcpp entrypoint still accepts its pointer for
compatibility; the new private component test does not require one. Thus seven
of seven now have facade-independent explicit-input component boundaries,
compared with five of seven facade-independent and four of seven ordinary-input
entries at baseline. This is still an interface inventory, not complete coverage
of all numerical states.

The before/after dependency map is described in `ARCHITECTURE.md`. Public S4
wrappers now prepare one transient context, dispatch to separate permutation and
CRT plans, and commit returned results/cache values. Results no longer obtain
statistical settings from the S4 object except in thin compatibility adapters.
Storage never applies QC policy, fitting never reads a store or chooses a work
plan, and native inference never loads an S4 object or serialized handle.

Deletion and migration checks found no production callers of the removed
precomputation/table aliases or optional worker reconstruction routes. The lean
`make_treatment_index_lookup()` remains an actual canonical-index consumer;
real pipeline result adapters also remain. The old combined native inference
body and old storage/indexing branches
were replaced, rather than left alongside an unused new abstraction. New
function/module counts can increase as these responsibilities become explicit;
the final metric report must retain that cost along with any complexity gains.

## RF9 checkpoint measurements before the reviewed cleanup

The 90664a7 production checkpoint was measured with the same pinned analyzers and
rules as RF0. Detailed results and per-file hashes are in
`results/maintainability-refactor/`; the handwritten-source SHA256 is
`1c55411ee3e352056834972d78ec23d5bac557729f203a1bf2a22679ccd2d36d`.
The commit label is `90664a7`, whose production R/native files are identical to
`28df390`. Subsequent cleanup and native retirement require a fresh final-source
measurement; these numbers remain historical. Documentation and validation-tool
changes do not enter this metric.

| Within-language metric | Original | Refactored |
| --- | ---: | ---: |
| R exclusive function units | 184 | 226 |
| R exclusive complexity, maximum / p90 / p95 | 49 / 9.70 / 13.85 | 49 / 9.00 / 13.00 |
| R named top-level owners | 122 | 164 |
| R inclusive named-owner complexity, p95 | 22.75 | 16.85 |
| R control nesting, maximum | 7 | 7 |
| Handwritten C++ functions | 43 | 52 |
| C++ complexity, maximum / p90 / p95 | 17 / 8.80 / 13.60 | 17 / 8.00 / 11.80 |
| C++ nesting, maximum | 5 | 5 |

Permutation driver exclusive/inclusive complexity changes from 7/15 to 10/14;
CRT changes from 9/23 to 11/18. The outer drivers gain explicit compatibility
guards/context handling while nested orchestration becomes smaller. Maxima and
maximum nesting remain unchanged; the refactor does not claim every function
was simplified. More named functions and compatibility delegates are a cost of
the explicit boundaries. The ownership/caller migrations and seven-operation
testability walkthrough above establish the architectural improvement; static
complexity alone does not establish readability or comprehension time.

The identical site scanner now records six private ondisc calls in one file
(previously four), 29 raw-storage candidate lines in three files (previously
51/nine), one narrow view-remapping candidate in one file (previously 17/four),
one cache-reset owner and one facade-writeback owner (previously three each).
`results/architecture-sites-refactor.json` retains every matching line. The
remaining raw-storage candidates are the adapter and legitimate import/output
construction, not statistical consumers. These deliberately narrow patterns are
not a complete coordinate or dependency inventory. Cognitive complexity, clone
coverage, human timing and a complete dynamic call graph remain unmeasured.
