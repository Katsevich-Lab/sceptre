# sceptre 0.99.2

## New features

-   A `sceptre_object` now stores its data in a `SingleCellExperiment` (response counts in the `counts` assay, gRNA counts in the `"grna"` alternative experiment, cell covariates in `colData`), which replaces the slots `response_matrix`, `grna_matrix` and `covariate_data_frame`. `get_response_matrix()`, `get_grna_matrix()` and `get_cell_covariates()` return the same values as before. `updateObject()` converts objects saved by earlier versions; until it is called on such an object, the pipeline functions, `get_response_matrix()`, `get_grna_matrix()` and `get_cell_covariates()` stop with a message naming it. A `sceptre_object` records the version of sceptre that created it.

-   New `import_data_from_sce()` imports a `SingleCellExperiment` (gRNA counts in an alternative experiment, gRNA targets and the regions they target in its `rowData`, the latter in the columns `target_chr`, `target_start` and `target_end`).

-   Disk-backed analysis no longer uses `ondisc`: with `backed = TRUE` and `directory_to_write`, `import_data()`, `import_data_from_cellranger()` and `import_data_from_sce()` keep the response and gRNA counts on disk in `screen.h5mu`, a MuData file that Bioconductor's MuData and Python's mudata can read. New `write_sceptre_object()` and `read_sceptre_object()` save and load a `sceptre_object`, standard or disk-backed, as a directory holding `screen.h5mu` and `sceptre_object.rds`. An existing `screen.h5mu`, which a disk-backed `sceptre_object` may be reading, is never replaced: the directory must not already hold one, unless the `sceptre_object` being written reads its counts from it.

-   With `backed = TRUE`, `import_data_from_cellranger()` holds the matrix of one directory in memory at a time and writes the response counts to `screen.h5mu` through temporary files in `directory_to_write`, so the whole screen is never in memory at once.

-   `set_analysis_parameters()` gains `seed` (default `NULL`: drawn from R's random number generator, so `set.seed()` makes an analysis reproducible); the seed is recorded and printed by `print()` and in `analysis_summary.txt`. In sceptre 0.99.2, `seed = 4` reproduces the results of 0.99.1 and earlier, which used the fixed seed 4, apart from the bug fixes below.

-   `assign_grnas()`, `run_calibration_check()`, `run_power_check()` and `run_discovery_analysis()` run in parallel through BiocParallel: the new argument `BPPARAM` takes a BiocParallel back end (for example `BiocParallel::MulticoreParam(4)`, or `BiocParallel::SnowParam(4)` on Windows); results do not depend on the back end or number of workers.

-   The `sceptre_object` class is documented (`?"sceptre_object-class"`) and has a validity method.

-   The vignette uses BiocStyle.

## Changes

-   Plots that show a random subset of points or gRNAs now draw it from R's random number generator, so `set.seed()` makes them reproducible; previously they always used the seed 4.

-   The default formula gives a message naming each covariate it leaves out (non-count covariates with 15 or more distinct values).

-   sceptre now requires R 4.5.0 or later (its dependencies do).

## Bug fixes

-   With `resampling_approximation = "no_approximation"`, `run_calibration_check()` no longer changes the number of resamples used by a later `run_discovery_analysis()` or `run_power_check()`, so their p-values no longer depend on whether the calibration check ran first.

-   The maximum gRNA assignment method requires `min_grna_n_umis_threshold` to be at least 1; 0 assigned gRNAs to cells with no gRNA UMIs.

-   In high MOI, `set_analysis_parameters(control_group = "nt_cells")` gives a warning that the complement set is used as the control group; previously the complement was used silently.

-   `import_data_from_cellranger()` no longer discards a user-supplied `batch` column of `extra_covariates`: with one directory it is kept (previously it was dropped), and with several directories, where `batch` records the directory of each cell, it is an error (previously it was overwritten). The levels of `batch` (`b1`, `b2`, ...) follow the order of `directories`; previously, with ten or more directories, they were sorted as text.

-   Importing a response or gRNA matrix with missing, negative, infinite or non-integer entries is an error; previously such counts were accepted (non-integer counts then gave many warnings in the analysis).

-   `construct_cis_pairs()` stops with an error when the chromosomes of the gRNA targets and of the responses are named differently (for example, `"chr1"` and `"1"`); previously it returned no pairs.

-   `import_data_from_cellranger()` reads compressed (`.gz`) Cell Ranger files without the R.utils package, which sceptre does not depend on; previously it stopped with an error when R.utils was not installed.

-   When `grna_target_data_frame` gives `vector_id` together with `chr`, `start` and `end`, each vector is now one unit whose region spans those of its gRNAs. Previously a vector was listed once per gRNA: under every `grna_integration_strategy` the calibration check's default `calibration_group_size` was too large, and in addition, with `"singleton"` the vector's pairs were analyzed and reported repeatedly and with `"bonferroni"` its p-values were too large. The gRNAs of a vector must target one chromosome; other columns whose values differ between them are dropped, with a message. A factor `vector_id` no longer leaves the collapsed gRNA counts without row names, which made the analysis fail.

-   `run_qc(additional_cells_to_remove =)` checks its argument: it must be a vector of integer indices between 1 and the number of cells, and a logical vector or cell IDs are an error. Previously a logical vector was read as the indices 0 and 1, so it removed only the first cell while the QC metrics reported every cell marked `TRUE` as removed, and out-of-range indices were silently ignored.

-   `import_data()` and `import_data_from_sce()` require `response_names` to have one element per response. Previously a vector of another length was accepted: a shorter one computed `response_p_mito` from the wrong responses, and a longer one with mitochondrial names beyond the last response could crash R.

## Deprecated and removed

The deprecated arguments and functions below still work in this version, with a warning; they will be made defunct in the next Bioconductor release and removed in the one after. Reading `.odm` files and the slots of the Nextflow pipeline are removed in this version.

-   The arguments `parallel` and `n_processors` of `assign_grnas()`, `run_calibration_check()`, `run_power_check()` and `run_discovery_analysis()` are deprecated in favor of `BPPARAM`; `parallel = TRUE` still runs on `n_processors` forked workers, with a warning.

-   The argument `use_ondisc` is deprecated in favor of `backed`, and `read_ondisc_backed_sceptre_object()` and `write_ondisc_backed_sceptre_object()` in favor of `read_sceptre_object()` and `write_sceptre_object()`; `.odm` files written by earlier versions cannot be read, so re-import the data. The `ondisc` package is no longer used.

-   The `sceptre_object` slots used by the Nextflow pipeline (`nf_pipeline`, `nuclear`, `integer_id`, `elements_to_analyze`) are removed.

# sceptre 0.99.1

-   Two-sided empirical p-values are now capped at 1. Previously, when the observed test statistic tied with many of the resampled statistics (for example, a response with the same count in every cell), the reported p-value could exceed 1, up to 2 ([issue #237](https://github.com/Katsevich-Lab/sceptre/issues/237)). One-sided p-values, p-values at or below 1 and significance calls are unchanged. With `grna_integration_strategy = "bonferroni"`, if two or more of a target's gRNAs have p-values capped at 1, the combined result now takes its fold change and nonzero-cell counts from the first of these gRNAs rather than from the one with the smallest uncapped p-value; its p-value (1) and significance call are unchanged.

-   In high MOI with `resampling_mechanism = "permutations"`, the calibration check no longer requests more permuted cells than exist when non-targeting gRNA groups overlap. Previously the resampling step could read and write outside its memory in this case. The high-MOI default, the CRT, was not affected.

-   The calibration check now compares `calibration_group_size` with the number of non-targeting gRNAs that still have cells after QC, and stops with a clear error if the size is too large. Previously it counted the non-targeting gRNAs in `grna_target_data_frame`; a group size larger than the number of non-targeting gRNAs with cells could make the calibration check write outside allocated memory and then fail with a misleading error, and with NT-cell controls and one remaining non-targeting gRNA, the calibration check ran without control cells.

-   `assign_grnas(method = "mixture")` now checks that both entries of `g_pert_guess_range` lie in [0, 20], as documented; previously it checked `pi_guess_range` by mistake, so a larger `g_pert_guess_range` was accepted. Values such as `c(10, 5000)` now give an error (the default, `log(c(10, 5000))`, is within range).

-   With `parallel = TRUE`, an error in a parallel worker is now reported with its original message, and a worker that exits without returning results (for example, because it ran out of memory) now causes an error. Previously the original error was replaced by an unrelated one, and a worker that exited left its results missing with only a warning.

-   `plot_run_qc()` now applies its `downsample_pairs`, `point_size` and `transparency` arguments to the pairwise QC plot; previously they were ignored.

-   `plot_grna_count_distributions()` no longer fails when `threshold` is above 11 and the largest plotted count lies between 11 and `threshold - 2`.

-   In `plot_run_discovery_analysis()`, the volcano plot now colors points by the `significant` column of the discovery result. Previously the significant pair with the largest p-value was colored as not significant, and with no discoveries all points had a missing color (grey); they now have the color of non-significant points.

-   `plot_run_discovery_analysis()` no longer drops valid pairs from its plots when the discovery analysis was run with `output_amount = 2` or `3`. Previously rows with missing optional diagnostic columns were removed, so the plots could show far fewer pairs than the result table, and in some settings none.

-   `import_data()` no longer reorders the slots of a `dgTMatrix` or `lgTMatrix` supplied as the response or gRNA matrix; the caller's matrix is left unchanged.

-   `set_analysis_parameters()` now discards saved response-model fits when the control group or the design matrix changes. Previously, in low MOI, after running the calibration check with NT-cell controls and then switching to `control_group = "complement"` with the same formula, the discovery analysis reused fits estimated on the NT cells only (and vice versa), which changed its p-values. Changing factor contrasts with the same formula also left stale fits in place. Analyses that call `set_analysis_parameters()` once are unaffected by this change.

-   `construct_cis_pairs()` now finds every gene within the distance threshold when the gene positions are not sorted within a chromosome. Previously it could omit nearby genes. This affects the bundled `gene_position_data_frame_grch37`, which is not sorted, and unsorted custom position tables; `gene_position_data_frame_grch38` is sorted and gives the same pairs as before.

-   `construct_cis_pairs()` no longer overflows integer arithmetic for genomic coordinates above about one billion.

-   A formula with only an intercept (`formula_object = ~ 1`) now works with NT-cell controls in the calibration check and discovery analysis; previously these failed with an error.

-   `construct_trans_pairs()` with `pairs_to_exclude = "pc_pairs"` no longer converts the caller's `positive_control_pairs` data frame to a `data.table`.

-   `import_data_from_cellranger()` (with `use_ondisc = FALSE`) and `import_data_from_parse()` now read Matrix Market files with eminem, a C++ parser provided by the Bioconductor package assorthead, and accept general coordinate matrices with integer or real values, as Cell Ranger and Parse write them; other Matrix Market variants, such as symmetric or pattern matrices, are now an error. They now import files whose entries are not sorted by cell; previously, counts from such files could be lost or assigned to the wrong cell or feature. Files sorted by cell and then by feature, such as the bundled example files, give the same result as before.

-   `import_data_from_cellranger()` (with `use_ondisc = FALSE`) and `import_data_from_parse()` now stop with an error, naming the file and line, when a matrix file is malformed: for example, when an entry lies outside the dimensions declared in the file's header or a coordinate is not an integer (such as `2.0`). Previously such entries could be dropped or moved silently to another cell. Both functions also stop when a matrix file lists the same (cell, feature) entry more than once, and `import_data_from_cellranger()` stops when a matrix file and its `features.tsv` file have different numbers of rows.

-   `import_data_from_cellranger()` (with `use_ondisc = FALSE`) now assigns counts to the correct rows when gene and gRNA features are interleaved in `features.tsv`, rather than listed in two contiguous blocks.

-   `import_data()` no longer loses entries of a dense (base R) square response or gRNA matrix that is symmetric or nearly symmetric.

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
