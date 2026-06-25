# Get result

`get_result()` returns a data frame containing the result of a
calibration check, power check, or discovery analysis. We pass as
arguments `sceptre_object` and `analysis`, where the latter is a string
indicating the function whose results we are querying. The output is a
data frame, the rows of which correspond to target-response pairs, and
the columns of which are as follows: `response_id`, `grna_target`,
`n_nonzero_trt`, `n_nonzero_cntrl`, `pass_qc` (a `TRUE`/`FALSE` value
indicating whether the pair passes pairwise QC), `p_value`,
`fold_change`, `se_fold_change` (standard error for the fold change
estimate), `log_2_fold_change`, and `significant` (a `TRUE`/`FALSE`
value indicating whether the pair is called as significant). The p-value
contained within the `p_value` column is a raw (i.e.,
non-multiplicity-adjusted) p-value. See [Section 8 of the introductory
chapter in the
manual](https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory)
for more information about this function.

## Usage

``` r
get_result(sceptre_object, analysis)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- analysis:

  a string indicating the name of the analysis whose results we are
  querying, one of `"run_calibration_check"`, `"run_power_check"`, or
  `"run_discovery_analysis"`.

## Value

a data frame containing the results of the analysis

## Note

If `output_amount` is set to `2` or `3` in
[`run_calibration_check()`](https://katsevich-lab.github.io/sceptre/reference/run_calibration_check.md),
[`run_power_check()`](https://katsevich-lab.github.io/sceptre/reference/run_power_check.md),
or
[`run_discovery_analysis()`](https://katsevich-lab.github.io/sceptre/reference/run_discovery_analysis.md),
then the result data frame contains additional columns; see [Chapter 6
in the
manual](https://timothy-barry.github.io/sceptre-book/run-calibration-check.html#sec-run_calibration_check_output_amount)
for more information.

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
pc_result <- sceptre_object |>
  set_analysis_parameters(
    side = "left",
    resampling_mechanism = "permutations",
    positive_control_pairs = positive_control_pairs
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_power_check() |>
  get_result("run_power_check")
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> The calibration check (`run_calibration_check()`) should be run before the power check.
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000224277 (1 of 5)
#> Analyzing pairs containing response ENSG00000226772 (5 of 5)
```
