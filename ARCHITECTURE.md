# Internal architecture

This is the lasting guide to how SCEPTRE's internal code fits together and where
to make changes. Keep it current as the package develops. The public functions,
defaults and saved `sceptre_object` remain the user interface; its 56 fields (S4
slots) hold data, settings, workflow status, reusable fits and results.
For tests, source checks and the distinction between permanent package material
and temporary project records, see [Contributing](CONTRIBUTING.md).

## Ownership and data flow

| Responsibility | Files | Inputs and outputs |
| --- | --- | --- |
| Public workflow and saved state | `R/s4_analysis_functs_1.R`, `R/s4_analysis_functs_2.R`, `R/workflow_state.R` | Validate supported requests, enforce prerequisites, commit successful outputs |
| Analysis preparation | `R/analysis.R`, `R/analysis_context.R` | Temporary lists of settings, requested pairs and ordered cell selections |
| Sparse/backed operations | `R/storage.R`, `src/sparse_counts.cpp`, `src/sparse_matrix_functs.cpp` | R 1-based source columns; representation operations with declared count conventions |
| QC measurements and eligibility | `R/pair_qc.R`, `src/pair_eligibility.h` | Aligned counts and the existing shared threshold predicate |
| Calibration candidate construction | `R/neg_control_functions.R`, `src/negative_control_functions.cpp` | Existing candidate stream, RNG and acceptance rule |
| Model fitting | `R/fitting.R` | Response/design or treatment/design arrays to compact fits/probabilities |
| Score preparation | `R/score_preparation.R` | Ordered response/design/compact fit to `mu`, `a`, `w`, `D` |
| Cache validity and updates | `R/fit_cache.R` | Supported mutation events, keyed compact fits, saving fits in the returned object |
| Distinct association plans | `R/association_execution.R`, `R/association_workers.R` | Response traversal for permutation; guide traversal for CRT |
| Worker mechanics | `R/execution.R` | Partition dispatch, progress/logs and result binding |
| Native inference | `src/inference_scores.cpp`, `src/inference_tail.cpp`, `src/inference_stages.cpp` | Numeric arrays and shared random draws; checked Rcpp entry in `src/low_level_full_test.cpp` |
| Results and finalization | `R/result_processing.R` | Native records to raw tables, failed-QC rows, guide integration and family correction |

![Main SCEPTRE modules, with the requests and returned values between them.](tools/architecture/module-map.svg)

An arrow points from the part asking for work to the part doing it; labels show
what is passed in and returned. The map groups related jobs and shows the main
calls, rather than every helper call or their chronological order. Workflow
functions also ask fit reuse to clear invalid fits and save completed fits.
These are groups of files in one R package, not separate installed packages.


The plans share worker mechanics, not a universal statistical inner loop.
Permutation analysis loads each response and reuses its preparation where valid.
CRT prepares shared response fits, then fits one assignment model and generates
the full draw bank once per active guide/view, reusing it across that guide's
responses. NT-control discovery retains guide-specific response fits.

## Coordinates and lifetime

`analysis_to_source` maps R analysis positions to original matrix columns.
`nt_to_analysis` retains the ordered NT-only view, which can be nonmonotone.
For example, suppose original columns 2, 4 and 7 survive QC. Their analysis
positions are 1, 2 and 3, so `analysis_to_source = c(2, 4, 7)`. If a fit uses
analysis positions 3 and 1, its counts and design rows must both use original
columns 7 and 2, in that order. Within that fit, they are positions 1 and 2.

Temporary guide-group memberships use analysis positions. Some saved NT-only
assignments instead number cells within the NT selection; `prepare_index_context()`
converts those once before calculations begin.

`make_fit_view()` supplies ordered source columns and local treatment positions.
Response extraction and design subsetting use the same view. Raw
source memberships are normalized in analysis order, preserving the old mask
semantics rather than the order of the raw vector. Request-row maps preserve
identity independently of response-sorted count scans and execution traversal.

Storage takes R 1-based source positions. The Rcpp `run_association_test()` entry
also accepts R 1-based observed treatment positions, converts them once, and
calls private components using zero-based positions. Synthetic draw banks are
already zero-based and never pass through that conversion.

Contexts, score arrays, assignment probabilities and native draw pointers are
transient. No full store copy or per-pair cell map enters the saved cache. Only
compact response coefficients/dispersion remain in its existing list. Its view
is the common analysis universe for complement controls or the NT-only
calibration universe for NT controls; guide-specific NT-discovery fits are not
persisted. Setters clear incompatible fits when assignments, analysis cells,
actual designs or control views change. Assignment probabilities and their
random draws are currently prepared during a run rather than saved in this list.

The saved `covariate_matrix` remains in full source-cell order. The transient
context selects its common design once, and guide-specific NT fits select their
own aligned rows from that view. Cache identity checks the actual source design
and fit-cell scope; changing a fitted response's name alone cannot establish
that an old coefficient/dispersion record is reusable.
Score vectors share that fit-cell order. `D` has one column per fit cell and
one row per eigencomponent in `eigen()` order; its row count equals the number
of design coefficients, but its rows are not original coefficient coordinates.

Native count kernels and calibration-group construction allocate their R
outputs before native compute owners. Calibration-pair and association returns
instead protect an R continuation before computation, then transfer the exact
selected native result into a raw R conversion leaf with cleanup on success or
R error. Compute owners end before that leaf; it allocates the actual accepted
rows/selected null segment, with no unused-stage buffer or bank copy. The
controlled ownership probes qualify those declared payload/scratch boundaries.

The private ondisc route is deliberate: backing pointers, identity-preserving
construction, reattachment and private kernels are confined to `storage.R`.
In-memory workflows do not need ondisc. Import-format parsers and construction
of the exported assignment matrix retain their format/output-specific sparse
construction code. Existing HDF5 and sparse native scans remain separate.
Dense import forces general shape before sparsification so both triangles and
near-symmetric values survive; the supported coercion can add one dense input
copy. Sparse conversion owns columns before sorting, reusing only constructor-
owned columns. In-memory vector aggregation prepares sparse row access and
sums each selected guide group in dense column blocks using the original
platform's dense accumulator. Blocks target 2^20 elements; one complete selected
column is the minimum and can exceed that target. Existing dense vector-by-cell
output/transpose workspace is additional. The complete sparse guide input is
not materialized densely. File-coordinate duplicates are
rejected by the owned Cell Ranger/Parse constructors; the external ondisc file
parser remains a separate boundary.
Response access materializes one complete source row before selecting a cell
view. Sparse QC computes the all-response NT summary. Backed regular QC can
request a fused NT summary; the actual backed nuclear/pod path omits that
summary and counts its requested response block. These are concrete
I/O and allocation boundaries, not promises of arbitrary out-of-core scalability.

## Where to make a later change

- **Output-compatible fitter:** start in `fitting.R`, whose response fitter
  currently executes Poisson `glm.fit()` plus a clipped dispersion estimate,
  and whose assignment fitter executes logistic `glm.fit()`. Score preparation
  owns the current projection and array orientation. Use the actual fit view;
  an all-cell response fit cannot replace a treatment-plus-NT fit. A different
  fitter or effective settings will also need a concrete cache-identity change
  and scientific validation when introduced.
- **QC policy:** measurements come from `collect_pair_qc_metrics()`. The current
  threshold rule lives in `pair_eligibility.h`; the vector adapter and native
  calibration sampler call it. A replacement requiring side, family or new
  measurements must define those inputs and calibration behavior together.
  Do not change only the discovery filter or replace the candidate stream.
- **Inference method:** aligned response/treatment inputs, design and current
  score/effect arrays are available in the workers. Implement the selected runner there and reuse raw
  result conversion/finalization where justified. SPA or Wilcoxon need not
  generate this method's `D` or draw bank. This guide does not select their
  statistical definitions or add empty extension hooks.
- **Storage adapter:** implement the small operations in `storage.R` and test
  values, order, repeats, shapes, explicit-zero/fractional semantics and caller
  nonmutation. Counts describe requested pairs; they need not allocate a full
  Cartesian table. Keep representation branching out of statistical workers.
  Enabling backed CRT later also requires separate parity/resource validation
  and review of the existing validation guard and distributed work plan.

SCEPTRE's sparse count kernels convert coordinates once into native-owned
scratch and leave R inputs unchanged. External ondisc kernels that decrement
indices receive disposable copies at the storage boundary. Current counters
count stored sparse positions, including explicit zeros. Storage conversion of
logical sparse input gives every stored entry the legacy value one, including
stored FALSE; vector-count aggregation instead follows Matrix's logical values.
Ondisc construction coerces fractional values to integers. Those legacy rules
are characterized rather than silently normalized.

## Results and component compatibility

`pair_results_from_native()` preserves output amounts 1/2/3 and legacy fields,
including `se_fold_change`. `assemble_pair_results()` restores failed-QC rows
and the output ordering. `integrate_guides()` and
`finalize_association_result()` take plain tables and explicit settings; neither
needs an S4 object. Standard calibration/discovery integrate and correct the
complete family after collection; power analysis applies grouping alone.
Nuclear pipeline shards remain raw Parquet without a
global correction or significance column.
Stage-1 and empirical-tail diagnostics, and padded null-draw columns, can be
missing in valid results. Discovery plotting selects complete p-value, effect
and significance fields; optional diagnostic missingness does not remove rows.

Nextflow calls internal QC/calibration phase functions and result processors
directly. Preserve those entrypoints and their intermediate table/file formats
until changes are coordinated with the pipeline. Standard collection combines
all worker outputs before adjustment; nuclear tasks retain their raw outputs.
`make_treatment_index_lookup()` keeps one set of cell lookups, rather than
storing a cell selection for every guide at once. Callers supply the prepared
inputs explicitly.

## Qualified companion versions

| Companion | Qualified source | Supported scope |
| --- | --- | --- |
| ondisc | `589803edea2a66ad15e7f07186bee16c7739ab8c` | Backed permutation analyses; see the optional installation instructions in README |
| sceptre-pipeline | `38dc760a5b4105d1c3f339cd45f2c98d971087f4` | Standard and nuclear workflows using the retained private phase functions and artifact formats |

The recorded pipeline comparison ran locally with matched arrival/seed settings.
It does not qualify containers, cluster schedulers, array execution or uncontrolled
arrival ordering. Keep the tested version matrix current when changing either
repository. Private phase functions are compatibility entry points for this
consumer, not new public exports. Coordinate their retirement with a tested
pipeline migration. Nuclear artifacts remain raw per-pod Parquet; standard
collection applies complete-family correction. Neither gains backed CRT here.

## Trace one association

Start with `run_discovery_analysis()` in `s4_analysis_functs_2.R`. It owns workflow
checks and successful saving successful results and fits; `analyze_association_pairs()` builds
the transient context, filters the eligible family and dispatches its plan.
`run_permutation_analysis()` visits responses, whereas `run_crt_analysis()` first
prepares response fits and then visits guides. The four functions in
`association_workers.R` express shared versus guide-specific fits explicitly.
They use storage reads and aligned fit/score inputs, then call `run_association_test()`
in `low_level_full_test.cpp`; score/effect arithmetic and selected tail segments
belong to the inference components, not the public workflow. Result assembly restores
failed-QC requests before complete-family guide integration/correction, or
retains raw output for the pipeline collector.

For a small numerical trace, the observed/synthetic-coordinate test in
`tests/testthat/test-refactor-indexing.R` independently specifies cell membership,
builds score components and checks the same native entry point. Its supplied
draws separate the one-based observed boundary from zero-based synthetic
positions without relying on a sampler. Permanent tests cover declared stages, effect formulas and sampler laws.

## Verification for contributors

Run targeted tests for the boundary being changed, then the package suite at a
coherent integration checkpoint:

```r
testthat::test_local(filter = "refactor-indexing")
devtools::check()
```

Use `SCEPTRE_REQUIRE_ONDISC=true` for required backed coverage. Run private native
component checks with `Rscript tools/checks/test-native-components.R .` from a
source checkout; [tools/checks/README.md](tools/checks/README.md) describes the
other current-source safety checks used by CI. Routine package checks use
current code, explicit mathematical expectations and current-mode comparisons;
they require no old package or historical comparison outputs. Regression tests
remain useful after the change that introduced them is complete.
Use separate source trees for concurrent compilation/check/install jobs.

Preserve supported workflows and the public API. Numerical or method changes
need their own review and validation, including the cells used for each fit
and whether saved fits remain reusable.
