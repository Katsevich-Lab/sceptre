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

Please note that v0.10.0 is a higher version number than v0.9.2.

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
