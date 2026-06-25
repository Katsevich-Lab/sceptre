# Run power check

`run_power_check()` runs the power check. The power check entails
applying `sceptre` to analyze positive control pairs — pairs for which
we know there is an association between the target and the response — to
ensure that `sceptre` is capable of detecting true associations. See
[Chapter 6 of the
manual](https://timothy-barry.github.io/sceptre-book/run-power-check-and-discovery-analysis.html)
for more detailed information about this function.

## Usage

``` r
run_power_check(
  sceptre_object,
  output_amount = 1,
  print_progress = TRUE,
  parallel = FALSE,
  n_processors = "auto",
  log_dir = tempdir()
)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- output_amount:

  (optional; default `1`) an integer taking values 1, 2, or 3 specifying
  the amount of information to return. `1` returns the least amount of
  information and `3` the most.

- print_progress:

  (optional; default `TRUE`) a logical indicating whether to print
  progress updates

- parallel:

  (optional; default `FALSE`) a logical indicating whether to run the
  function in parallel. `parallel = TRUE` is recommended only on Mac; it
  is not supported on Windows and may behave unreliably on Linux
  clusters.

- n_processors:

  (optional; default `"auto"`) an integer specifying the number of
  processors to use if `parallel` is set to `TRUE`. The default,
  `"auto"`, uses half the physical cores. The fraction may be tuned via
  the `parallelly.availableCores.fraction` R option.

- log_dir:

  (optional; default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)) a string
  indicating the directory in which to write the log files (ignored if
  `parallel = FALSE`)

## Value

an updated `sceptre_object` in which the power check has been carried
out

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

# set analysis parameters, assign grnas, run qc
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
sceptre_object <- sceptre_object |>
  set_analysis_parameters(
    side = "left",
    resampling_mechanism = "permutations",
    positive_control_pairs = positive_control_pairs
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_power_check()
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> The calibration check (`run_calibration_check()`) should be run before the power check.
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000224277 (1 of 5)
#> Analyzing pairs containing response ENSG00000226772 (5 of 5)
```
