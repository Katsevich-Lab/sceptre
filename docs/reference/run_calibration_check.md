# Run calibration check

`run_calibration_check()` runs the calibration check. The calibration
check involves applying sceptre to analyze negative control
target-response pairs — pairs for which we know there is no association
between the target and response — to ensure control of the false
discovery rate. The calibration check enables us to verify that the
discovery set that sceptre ultimately produces is not contaminated by
excess false positives. See [Chapter 5 of the
manual](https://timothy-barry.github.io/sceptre-book/run-calibration-check.html)
for more detailed information about this function.

## Usage

``` r
run_calibration_check(
  sceptre_object,
  n_calibration_pairs = NULL,
  calibration_group_size = NULL,
  print_progress = TRUE,
  parallel = FALSE,
  n_processors = "auto",
  log_dir = tempdir(),
  output_amount = 1
)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- n_calibration_pairs:

  (optional) the number of negative control pairs to construct and test
  for association

- calibration_group_size:

  (optional) the number of negative control gRNAs to randomly assemble
  to form each negative control target

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

- output_amount:

  (optional; default `1`) an integer taking values 1, 2, or 3 specifying
  the amount of information to return. `1` returns the least amount of
  information and `3` the most.

## Value

an updated `sceptre_object` in which the calibration check has been
carried out

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

# set analysis parameters, assign grnas, run qc, run calibration check
sceptre_object <- sceptre_object |>
  set_analysis_parameters(
    side = "left",
    resampling_mechanism = "permutations"
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_calibration_check(
    n_calibration_pairs = 500,
    calibration_group_size = 2
  )
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Constructing negative control pairs.
#>  ✓
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000253631 (1 of 96)
#> Analyzing pairs containing response ENSG00000100053 (5 of 96)
#> Analyzing pairs containing response ENSG00000100325 (10 of 96)
#> Analyzing pairs containing response ENSG00000253963 (15 of 96)
#> Analyzing pairs containing response ENSG00000100314 (20 of 96)
#> Analyzing pairs containing response ENSG00000177993 (25 of 96)
#> Analyzing pairs containing response ENSG00000099956 (30 of 96)
#> Analyzing pairs containing response ENSG00000253920 (35 of 96)
#> Analyzing pairs containing response ENSG00000203280 (40 of 96)
#> Analyzing pairs containing response ENSG00000187792 (45 of 96)
#> Analyzing pairs containing response ENSG00000211666 (50 of 96)
#> Analyzing pairs containing response ENSG00000236611 (55 of 96)
#> Analyzing pairs containing response ENSG00000253546 (60 of 96)
#> Analyzing pairs containing response ENSG00000099917 (65 of 96)
#> Analyzing pairs containing response ENSG00000253889 (70 of 96)
#> Analyzing pairs containing response ENSG00000099889 (75 of 96)
#> Analyzing pairs containing response ENSG00000241973 (80 of 96)
#> Analyzing pairs containing response ENSG00000225783 (85 of 96)
#> Analyzing pairs containing response ENSG00000279548 (90 of 96)
#> Analyzing pairs containing response ENSG00000100068 (95 of 96)
```
