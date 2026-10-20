## Refactor review: empirical two-sided p-values

- Two-sided empirical p-values now cap twice the smaller inclusive add-one tail
  at 1. All-tied finite scores therefore return 1 instead of 2. One-sided
  p-values, values already at or below 1, and null-score ordering are unchanged.

## Refactor review: invalid observed scores

- Association comparisons with a nonfinite observed score (NaN or either
  infinity, including scores from zero or negative adjusted variance) now remain
  in the tested family with p-value 1. No null tail is evaluated for that pair.
  At output levels 2 and 3, stage 0 identifies this fallback and `z_orig` retains
  the original score. Tail parameters and level-3 null-score columns are missing
  for these pairs; an entirely stage-0 result has one missing `z_null_1` column.
  Finite scores and null-score semantics are unchanged. Input and runtime errors
  still raise errors rather than being converted to p-values.

# sceptre (development version)

- Legacy saved-table migration normalizes table containers while preserving
  columns and metadata. Backed imports choose smaller HDF5 chunks for stores
  with 2–1000 stored entries; stores with fewer than two entries fail before
  either backing file is created. The optional rhdf5 test dependency is declared.

- The distributed execution API now has runnable end-to-end examples. Internal
  validation and migration helpers are shorter, and warning capture uses local
  state without superassignment. These changes preserve numerical execution,
  result ownership and runtime certification.
- Public result extraction and completed-result admission now own nested tables
  and historical provenance. Editing returned data.tables cannot change saved
  results, earlier objects or later history; original RunPlan certificates are
  preserved. Failed or interrupted summary export restores the caller's sinks.
- Generated default formulas support nonsyntactic covariate names and use a
  minimal environment. User-supplied formula environments remain intact.
  Correction settings are validated before analysis: alpha must be one finite
  number in `(0, 1]`, and the method must be supported by `stats::p.adjust()`.
- Calibration now derives its resampling budget from its own requested family.
  In the empirical path this intentionally corrects the former leakage of
  discovery/power budgets into calibration. Scientific parity comparisons use
  the same explicit budgets when comparing the old and current implementations.
- Worker failures now signal `sceptre_worker_error` with block, analysis, phase
  and original cause metadata. Interrupts retain their condition; missing child
  replies remain failures and cannot produce partial completed results.
- Experimental `prepare_execution()`, `execute_work_unit()` and
  `collect_execution()` support complete-family distributed work without legacy
  pipeline slots. Global plans and draw scopes remain frozen; validated raw
  proposals are integrated and adjusted once. Runtime certificates bind actual
  R functions, S4 definitions/dispatch and loaded native library bytes.
- Backed operations diagnose incompatible optional ondisc interfaces before
  file/native work. Saved-object help explains closed readRDS inspection,
  process-local reopening, verified relocation and immutable active backing.
  The declared testthat minimum is now 3.1.8; additional R/Bioconductor support
  gates retain separate ondisc-free checks.

- Repeated CRT analyses can reuse compact treatment-model fits when ordered
  cells, design and treatment membership agree. Calibration bypasses this cache
  and CRT resamples stay fresh. Fitting warnings occur on actual fits; warm
  reuse retains their diagnostics without re-emitting them.

- Internal analysis code now has shared cell-coordinate and storage boundaries,
  separate fitting, score, QC policy, inference and result-processing modules,
  and common worker utilities. Public workflows retain their API; saved objects
  use schema 3 components with explicit migration from earlier layouts.
- Response fits are invalidated when the actual design or control-cell view
  changes. This prevents NT-only calibration fits from being reused for a later
  complement-control analysis and preserves agreement with a fresh fit.

- Count-file import now handles unsorted entries and interleaved Cell Ranger
  gene/guide rows correctly. In-memory Cell Ranger and Parse import reject
  repeated `(cell, feature)` coordinates with a clear error. This rule concerns
  file coordinates, not repeated guide metadata or existing R triplet inputs;
  the external ondisc file parser has its own behavior.
- Dense symmetric and nearly symmetric input matrices retain all supplied
  values. Sparse conversions and summary/getter operations no longer modify
  arrays shared with the caller's counts, indices or metadata.
- Cis-pair construction handles unsorted gene positions and avoids overflow
  when calculating genomic midpoints and distances.
- Intercept-only formulas now work with non-targeting controls. Calibration
  handles overlapping NT groups safely and checks the guides actually retained
  after QC before sampling groups. Invalid native requests fail clearly before
  unsafe allocation or indexing; file connections and native buffers are cleaned
  up on the checked error paths.
- Mixture-assignment input validation checks the intended perturbation bound.
  The known statistical inconsistency described below is retained pending a
  separately validated replacement. Dispersion fitting retains its numerical
  fallback while allowing unrelated system/allocation failures to propagate.
- QC plots honor their supplied options; large guide-count thresholds no longer
  break histogram construction. Volcano colors follow stored significance calls,
  QQ confidence bands appear in each facet, and missing optional diagnostics no
  longer remove valid discovery rows from plots.
- Guide-vector aggregation avoids making the whole sparse input dense while
  retaining its existing values, ordering and output shape.
- Parallel analysis reports a failed worker's original error and block number
  before combining results. A child exiting without replying cannot silently drop a block of work;
  successful empty results remain supported.

# sceptre 0.99.0 (2026-06-24)

This release prepares `sceptre` for submission to Bioconductor.

-   On Windows, setting `parallel = TRUE` now raises a clear error explaining that sceptre's parallel mode is unsupported on Windows (it is fork-based). `ondisc`/HDF5 file paths are also normalized to avoid Windows-specific issues.

-   Core detection for parallel processing now uses the `parallelly` package instead of `parallel::detectCores()`, giving more reliable core counts in containers and on computing clusters.

-   The plotting functions were updated for compatibility with ggplot2 4.0.

-   The bundled `import_data_from_parse()` example data is now a documented, reproducible subset of a public Parse Biosciences dataset (CC BY 4.0), generated by a script in `inst/scripts/`.

-   The documentation now describes a known fitting inconsistency in the
    historical `"mixture"` gRNA-assignment method. The implementation estimates
    the perturbation parameter under a multiplicative parameterization but
    applies it additively in the subsequent E-step. Consequently, the resulting
    parameter estimates are not maximum-likelihood estimates for the mixture
    model used to compute the fitted probabilities. The historical behavior is
    retained because it has performed well in benchmarks; the `"maximum"` and
    `"thresholding"` methods are unaffected. We plan to replace the historical
    mixture implementation after a statistically coherent and empirically
    validated alternative is available. Progress is tracked in
    [issue #220](https://github.com/Katsevich-Lab/sceptre/issues/220).

# sceptre 0.10.3 (2025-08-28)

In Version 0.10.3, we added standard errors to fold change estimates; see new `fold_change` and `se_fold_change` columns of association analysis outputs. In addition, we fixed bugs reported in issues [#162](https://github.com/Katsevich-Lab/sceptre/issues/162), [#164](https://github.com/Katsevich-Lab/sceptre/issues/164), and [#187](https://github.com/Katsevich-Lab/sceptre/issues/187) as well as a bug reported via email that arose when `grna_target_data_frame` contained `NA` values in the `grna_id` or `grna_target` columns.

# sceptre 0.10.2 (2025-01-17)

Version 0.10.2 is a minor update in which we replaced the larger, real example
data by smaller, simulated example data. Furthermore, the example data now live 
in the `sceptre` package, rather than in the external `sceptredata` package. 
Finally, the R-CMD-CHECK badge was restored on the website.

# sceptre 0.10.1 (2024-04-30)

Version 0.10.1 includes minor updates to v0.10.0. A new gene position data frame 
`gene_position_data_frame_grch37` was added.

# sceptre 0.10.0 (2024-04-02)

Version 0.10.0 is a major update to the `sceptre` package. This version provides support for the analysis of large-scale single-cell CRISPR screen data. It also includes several other, minor updates.

-   We have introduced the `ondisc`-backed `sceptre_object`, which is a special kind of `sceptre_object` in which the data are stored on-disk as opposed to in-memory.

-    We have made the `sceptre` R package compatible with the `sceptre` Nextflow pipeline.

-   We have added the function `plot_response_grna_target_pair()`, which creates a violin plot of the expression level of a specified gene across treatment and control groups of a specified target.

-   We have added the function `get_grna_assignments()`, which facilitates obtaining the gRNA-to-cell assignments.

-   We have updated the `sceptre` e-book, adding two new parts: a part on at-scale `sceptre` and a part on the methodology underlying `sceptre`.

-   We have made the discovery pairs data frame and the positive control pairs data frame optional arguments to the `set_analysis_parameters()` function.

-   We have added a comprehensive suite of unit tests to help verify correctness of the code.

-   We have fleshed out the man pages, for example by adding a runable example to each.

-   We have moved the example data within the `sceptre` package into the companion `sceptredata` package.

-   We have issued minor bug fixes.

Please note that v0.10.0 is a higher version number than v0.9.2. Also, note that you will need to recreate your `sceptre_object` and rerun your analysis to use version 0.10.0. However, you should be able to use the exact same code to do so (assuming you currently are using v0.9.0 or higher).

# sceptre 0.9.2 (2023-12-08)

Version 0.9.2 is a minor update to version 0.9.0.

1. Add an `n_processors` argument to the functions that enable parallelization to allow users to select the number of processors to use. (The default, `n_processors = "auto"`, selects the number of processors to use based on the number of processors available on the machine.)
2. Add a `log_dir` argument, enabling users to specify the directory in which to write the log files.
3. Accelerate the `plot_assign_grnas()` function.
4. Fix the functionality for identifying mitochondrial genes; now, genes prefixed by "MT-" or "mt-" are considered mitochondrial.

# sceptre 0.9.1 (2023-10-24)

Version 0.9.1 is a minor update to version 0.9.0.

-   We have added an experimental `import_data_from_parse()` function to import data from the output of Parse Biosciences CRISPR Detect.

-   We have added support for the `"bonferroni"` gRNA integration strategy.

# sceptre 0.9.0 (2023-10-20)

Version 0.9.0 is a total rework of the `sceptre` package. The new version of the package has a fresh user interface and is faster, more memory-efficient, and more fully featured than previous versions. We summarize key updates here.

-   We have added a `sceptre_object` class to represent the single-cell CRISPR screen data.
-   We have unified low-MOI and high-MOI analysis into a single interface.
-   We have written a [manual](https://timothy-barry.github.io/sceptre-book/) to guide users through the entire process of analyzing their single-cell CRISPR screen data.
-   We have added a new mixture model method for assigning gRNAs to cells in a principled way.
-   We have made the experimental high-MOI functionality (from version 0.3.0) the default functionality for high-MOI analysis.
-   We have added functionality to carry out cell-wise QC.
-   Mac and Linux users now can run `sceptre` in parallel.
-   We have added a suite of plotting functions to help users visualize the output of the different steps of the `sceptre` pipeline.
-   We have added helper functions to facilitate the construction of *cis* and *trans* discovery sets.
-   We have added a function to import data into `sceptre` from the output of one or more calls to CellRanger count.

# sceptre 0.3.0 (2023-07-13)

Version 0.3.0 introduces a new, experimental high MOI function. We expect the experimental high MOI function to be faster, more memory efficient, and more powerful than the current high MOI function on most datasets. The current high MOI function likely will be deprecated in the next version of the package in favor of the experimental function. Please let us know about your experience using the experimental high MOI function, in particular whether you run into any bugs.

We also have added a new plotting function, namely `plot_resampling_distribution`. Small changes to the API of the `run_sceptre_lowmoi` function are detailed in the function documentation.

# sceptre 0.2.0 (2023-04-03)

Version 0.2.0 is our biggest update yet. We have added functionality for low MOI analysis! The low MOI module is based new statistical methods and computational algorithms.

# sceptre 0.1.0 (2022-03-10)

Version 0.1.0 is a major update to `sceptre`. Usability and speed are improved considerably.

### Usability

-   The function `run_sceptre_gRNA_gene_pair`, which was redundant, is now deprecated.

-   `run_sceptre_high_moi`(previously called `run_sceptre_in_memory`) is simpler to use: the function now has only four required arguments: `gene_matrix` (previously called `expression_matrix`), `gRNA_matrix` (previously called `expression_matrix`), `covariate_matrix`, and `gene_gRNA_pairs`. The formerly required arguments `storage_dir` and `side` are now set to `tempdir()` and "both" by default. Additionally, the argument `pod_sizes` is removed entirely (and handled internally).

-   `run_sceptre_high_moi` has the additional optional arguments `full_output` and `parallel`. `full_output` controls the complexity of the data frame outputted by the method. When `full_output` is set to FALSE (the default), `run_sceptre_high_moi` outputs a data frame with four columns only, all of which are easy to interpret: `gene_id`, `gRNA_id`, `p_value`, and `z_value`. `parallel` controls whether the function is parallelized (TRUE; default) or not (FALSE).

-   `run_sceptre_high_moi` now accepts a raw (i.e., unthresholded) gRNA matrix or a user-thresholded gRNA matrix.

-   A new auxiliary function `combine_gRNAs` combines gRNAs that target the same chromosomal site.

-   Numerous checks have been added to `run_sceptre_high_moi` ensure that the input is valid. For example, `run_sceptre_high_moi` checks that gene IDs and gRNA IDs in `gene_gRNA_pairs` are in fact subsets of the row names of `gene_matrix` and `gRNA_matrix`, respectively.

### Speed

Two accelerations have been implemented to improve speed. These accelerations do not affect the API of the package.

-   First, the test statistic used in the conditional randomization test has changed. Previously, the test statistic was a z-score derived from a *Wald* test of a fitted negative binomial GLM. Now, the test statistic is a z-score derived from a *score* test of the same negative binomial GLM, which is asymptotically equivalent to the former but more robust in finite samples. Additionally, this score test-based *z*-score is computed via an explicit formula, sidestepping the need to fit a GLM, as was done previously. Overall, the new test statistic is faster to compute and more robust than the previous test statistic.

-   The synthetic perturbation indicators are now generated as part of the gRNA precomputation, factoring out this somewhat time-intensive step from the pairwise tests of association.

# sceptre 0.0.1 (Initial release)
