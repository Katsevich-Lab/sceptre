# Write outputs to directory

`write_outputs_to_directory()` writes the outputs of a `sceptre`
analysis to a directory on disk. `write_outputs_to_directory()` writes
several files to the specified directory: a text-based summary of the
analysis (`analysis_summary.txt`), the various plots (`*.png`), the
calibration check, power check, discovery analysis results
(`results_run_calibration_check.rds`, `results_run_power_check.rds`, and
`results_run_discovery_analysis.rds`, respectively), and the binary
gRNA-to-cell assignment matrix (`grna_assignment_matrix.rds`). See
[Section 8 of the introductory chapter in the
manual](https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory)
for more information about this function.

## Usage

``` r
write_outputs_to_directory(sceptre_object, directory)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- directory:

  a string giving the file path to a directory on disk in which to write
  the results

## Value

No return value; called for its side effect of writing files to
`directory`.

## Examples

``` r
data(highmoi_example_data)
data(grna_target_data_frame_highmoi)
# import data
sceptre_object <- import_data(
  response_matrix = highmoi_example_data$response_matrix,
  grna_matrix = highmoi_example_data$grna_matrix,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = highmoi_example_data$extra_covariates,
  response_names = highmoi_example_data$gene_names
)
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
discovery_pairs <- construct_cis_pairs(sceptre_object,
  positive_control_pairs = positive_control_pairs,
  distance_threshold = 5e6
)
sceptre_object |>
  set_analysis_parameters(
    side = "left",
    resampling_mechanism = "permutations",
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_calibration_check() |>
  run_power_check() |>
  run_discovery_analysis() |>
  write_outputs_to_directory(paste0(tempdir(), "/sceptre_outputs"))
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Constructing negative control pairs.
#>  ✓
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000253631 (1 of 98)
#> Analyzing pairs containing response ENSG00000100053 (5 of 98)
#> Analyzing pairs containing response ENSG00000100325 (10 of 98)
#> Analyzing pairs containing response ENSG00000253963 (15 of 98)
#> Analyzing pairs containing response ENSG00000100314 (20 of 98)
#> Analyzing pairs containing response ENSG00000177993 (25 of 98)
#> Analyzing pairs containing response ENSG00000099956 (30 of 98)
#> Analyzing pairs containing response ENSG00000253920 (35 of 98)
#> Analyzing pairs containing response ENSG00000203280 (40 of 98)
#> Analyzing pairs containing response ENSG00000187792 (45 of 98)
#> Analyzing pairs containing response ENSG00000211666 (50 of 98)
#> Analyzing pairs containing response ENSG00000236611 (55 of 98)
#> Analyzing pairs containing response ENSG00000253546 (60 of 98)
#> Analyzing pairs containing response ENSG00000099917 (65 of 98)
#> Analyzing pairs containing response ENSG00000253889 (70 of 98)
#> Analyzing pairs containing response ENSG00000099889 (75 of 98)
#> Analyzing pairs containing response ENSG00000241973 (80 of 98)
#> Analyzing pairs containing response ENSG00000225783 (85 of 98)
#> Analyzing pairs containing response ENSG00000279548 (90 of 98)
#> Analyzing pairs containing response ENSG00000100068 (95 of 98)
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000224277 (1 of 5)
#> Analyzing pairs containing response ENSG00000226772 (5 of 5)
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000100218 (1 of 97)
#> Analyzing pairs containing response ENSG00000099958 (5 of 97)
#> Analyzing pairs containing response ENSG00000211638 (10 of 97)
#> Analyzing pairs containing response ENSG00000211685 (15 of 97)
#> Analyzing pairs containing response ENSG00000220891 (20 of 97)
#> Analyzing pairs containing response ENSG00000187905 (25 of 97)
#> Analyzing pairs containing response ENSG00000235954 (30 of 97)
#> Analyzing pairs containing response ENSG00000244486 (35 of 97)
#> Analyzing pairs containing response ENSG00000225783 (40 of 97)
#> Analyzing pairs containing response ENSG00000224277 (45 of 97)
#> Analyzing pairs containing response ENSG00000286326 (50 of 97)
#> Analyzing pairs containing response ENSG00000211672 (55 of 97)
#> Analyzing pairs containing response ENSG00000203280 (60 of 97)
#> Analyzing pairs containing response ENSG00000233521 (65 of 97)
#> Analyzing pairs containing response ENSG00000100319 (70 of 97)
#> Analyzing pairs containing response ENSG00000211674 (75 of 97)
#> Analyzing pairs containing response ENSG00000099917 (80 of 97)
#> Analyzing pairs containing response ENSG00000100053 (85 of 97)
#> Analyzing pairs containing response ENSG00000229770 (90 of 97)
#> Analyzing pairs containing response ENSG00000253920 (95 of 97)
# files written to "sceptre_outputs" in tempdir()
```
